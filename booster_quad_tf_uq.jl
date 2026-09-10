module BoosterUQ

using ChainRulesCore
using GTPSA
using LinearAlgebra
using SciBmad
import ReverseDiff
import Turing

export BoosterInferenceData,
    BoosterSimulator,
    GTPSASensitivity,
    H_BPM_NAMES,
    H_CORRECTOR_NAMES,
    H_QUAD_LABELS,
    QUAD_FACTOR_NAMES,
    V_BPM_NAMES,
    V_CORRECTOR_NAMES,
    V_QUAD_LABELS,
    booster_model,
    lognormal_quad_prior,
    orbit_and_quad_jacobian,
    predict_orbits,
    predict_orbits_and_jacobian,
    set_corrector_currents!,
    set_quad_factors!

# The translated lattice defines hundreds of names used by its DefExprs. Keep
# those names private; the simulator replaces the inferred magnet strengths
# with direct values after copying the lattice.
module _BoosterLatticeTemplate
include(joinpath(@__DIR__, "..", "booster_lattice", "booster_run.jl"))
end

const H_QUAD_LABELS = (
    :A2, :A4, :A6, :A8, :B2, :B4, :B6, :B8, :C2, :C4, :C6, :C8,
    :D2, :D4, :D6, :D8, :E2, :E4, :E6, :E8, :F2, :F4, :F6, :F8,
)
const V_QUAD_LABELS = (
    :A1, :A3, :A5, :A7, :B1, :B3, :B5, :B7, :C1, :C3, :C5, :C7,
    :D1, :D3, :D5, :D7, :E1, :E3, :E5, :E7, :F1, :F3, :F5, :F7,
)
const H_CORRECTOR_NAMES = Tuple(Symbol("DHC", label) for label in H_QUAD_LABELS)
const V_CORRECTOR_NAMES = Tuple(Symbol("DVC", label) for label in V_QUAD_LABELS)

const H_BPM_LABELS = (
    :A6, :A8, :B4, :B6, :C2, :C6, :D2, :D8, :E2, :E4, :E6, :E8,
    :F2, :F4, :F8,
)
const V_BPM_LABELS = (
    :A1, :A3, :A5, :A7, :B1, :B5, :B7, :C1, :C3, :C5, :D3, :D5,
    :E3, :E5, :E7, :F1, :F3, :F5,
)
const H_BPM_NAMES = Tuple("PUEH" * string(label) for label in H_BPM_LABELS)
const V_BPM_NAMES = Tuple("PUEV" * string(label) for label in V_BPM_LABELS)
const QUAD_FACTOR_NAMES = Tuple(vcat(
    ["QH_" * string(label) for label in H_QUAD_LABELS],
    ["QV_" * string(label) for label in V_QUAD_LABELS],
))

const N_QH = length(H_QUAD_LABELS)
const N_QV = length(V_QUAD_LABELS)
const N_QUADS = N_QH + N_QV
const N_HC = length(H_CORRECTOR_NAMES)
const N_VC = length(V_CORRECTOR_NAMES)
const N_BPMS = length(H_BPM_NAMES) + length(V_BPM_NAMES)

Base.@kwdef struct GTPSASensitivity{AD<:Turing.ADTypes.AbstractADType}
    closed_orbit_adtype::AD = Turing.AutoForwardDiff()
    reltol::Float64 = 1e-13
    abstol::Float64 = 1e-13
    maxiter::Int = 100
    warm_start::Bool = false
end

"""
    BoosterSimulator(; sensitivity=GTPSASensitivity())

An independently mutable Booster lattice with one cached first-order GTPSA
workspace. Use one simulator per process or concurrently executing chain.
"""
mutable struct BoosterSimulator{L,S<:GTPSASensitivity,D,T}
    lattice::L
    sensitivity::S
    descriptor::D
    phase_variables::Vector{T}
    quad_parameters::Vector{T}
    quad_indices::Vector{Int}
    quad_base_strengths::Vector{Float64}
    hcorrector_indices::Vector{Int}
    vcorrector_indices::Vector{Int}
    bpm_by_index::Dict{Int,Tuple{Int,Int}}
    closed_orbit_guess::Matrix{Float64}
    coasting_beam::Bool
end

function BoosterSimulator(; sensitivity::GTPSASensitivity=GTPSASensitivity())
    template = _BoosterLatticeTemplate.booster
    lattice = Beamline(
        deepcopy(collect(template.line));
        species_ref=template.species_ref,
        p_over_q_ref=template.p_over_q_ref,
    )
    index = Dict(Symbol(element.name) => i for (i, element) in enumerate(lattice.line))

    quad_indices = [
        [index[Symbol("QH", label)] for label in H_QUAD_LABELS]
        [index[Symbol("QV", label)] for label in V_QUAD_LABELS]
    ]
    quad_base_strengths = Float64[
        [
            _BoosterLatticeTemplate.CBLH(
                getfield(_BoosterLatticeTemplate, Symbol("IQH", label))
            )() for label in H_QUAD_LABELS
        ]
        [
            (label in (:D5, :F5) ?
             _BoosterLatticeTemplate.CBLVear :
             _BoosterLatticeTemplate.CBLV)(
                getfield(_BoosterLatticeTemplate, Symbol("IQV", label))
            )() for label in V_QUAD_LABELS
        ]
    ]

    # A zero index represents the two requested horizontal correctors which do
    # not exist or are turned off in the physical lattice.
    hcorrector_indices = [get(index, name, 0) for name in H_CORRECTOR_NAMES]
    vcorrector_indices = [index[name] for name in V_CORRECTOR_NAMES]

    bpm_names = (H_BPM_NAMES..., V_BPM_NAMES...)
    h_bpm_names = Set(H_BPM_NAMES)
    bpm_by_index = Dict(
        index[Symbol(name)] => (row, name in h_bpm_names ? 1 : 3)
        for (row, name) in enumerate(bpm_names)
    )

    descriptor = Descriptor(4, 1, N_QUADS, 1)
    phase_variables = GTPSA.vars(descriptor)
    quad_parameters = GTPSA.params(descriptor)

    simulator = BoosterSimulator(
        lattice,
        sensitivity,
        descriptor,
        phase_variables,
        quad_parameters,
        quad_indices,
        quad_base_strengths,
        hcorrector_indices,
        vcorrector_indices,
        bpm_by_index,
        zeros(1, 6),
        false,
    )
    set_corrector_currents!(simulator, zeros(N_HC), zeros(N_VC))
    set_quad_factors!(simulator, ones(N_QH), ones(N_QV))
    simulator.coasting_beam =
        SciBmad.coast_check(lattice, sensitivity.closed_orbit_adtype)
    return simulator
end

function set_corrector_currents!(
    simulator::BoosterSimulator,
    currents_hc::AbstractVector,
    currents_vc::AbstractVector,
)
    length(currents_hc) == N_HC ||
        throw(ArgumentError("currents_hc must contain $N_HC values"))
    length(currents_vc) == N_VC ||
        throw(ArgumentError("currents_vc must contain $N_VC values"))

    calibration = _BoosterLatticeTemplate.CorCalib
    for (element_index, current) in zip(simulator.hcorrector_indices, currents_hc)
        iszero(element_index) && continue
        simulator.lattice.line[element_index].Bn0L = calibration * current
    end
    for (element_index, current) in zip(simulator.vcorrector_indices, currents_vc)
        simulator.lattice.line[element_index].Bs0L = calibration * current
    end
    return nothing
end

function set_quad_factors!(
    simulator::BoosterSimulator,
    vars_qh::AbstractVector,
    vars_qv::AbstractVector;
    tpsa::Bool=false,
)
    length(vars_qh) == N_QH ||
        throw(ArgumentError("vars_qh must contain $N_QH values"))
    length(vars_qv) == N_QV ||
        throw(ArgumentError("vars_qv must contain $N_QV values"))

    for (i, factor) in enumerate(Iterators.flatten((vars_qh, vars_qv)))
        value = tpsa ? factor + simulator.quad_parameters[i] : factor
        simulator.lattice.line[simulator.quad_indices[i]].Bn1L =
            simulator.quad_base_strengths[i] * value
    end
    return nothing
end

"""
    orbit_and_quad_jacobian(simulator, currents_hc, currents_vc, vars_qh, vars_qv)

Return the 33 selected BPM readings in millimetres and their `33 × 48`
Jacobian with respect to the quadrupole factors.
"""
function orbit_and_quad_jacobian(
    simulator::BoosterSimulator,
    currents_hc::AbstractVector,
    currents_vc::AbstractVector,
    vars_qh::AbstractVector{<:Real},
    vars_qv::AbstractVector{<:Real},
)
    # GTPSA has a process-global current descriptor. With one simulator per
    # process, selecting the cached descriptor here is sufficient; no lock or
    # save/restore wrapper is needed.
    GTPSA.desc_current = simulator.descriptor

    set_corrector_currents!(simulator, currents_hc, currents_vc)
    set_quad_factors!(simulator, vars_qh, vars_qv)

    sensitivity = simulator.sensitivity
    sensitivity.warm_start || fill!(simulator.closed_orbit_guess, 0.0)
    result = find_closed_orbit(
        simulator.lattice;
        v0=simulator.closed_orbit_guess,
        coasting_beam=simulator.coasting_beam,
        autodiff=sensitivity.closed_orbit_adtype,
        reltol=sensitivity.reltol,
        abstol=sensitivity.abstol,
        maxiter=sensitivity.maxiter,
        warn=false,
        batch=Val(false),
    )
    all(iszero, result.sol.retcode) ||
        throw(ErrorException("closed-orbit finder did not converge"))
    closed_orbit = result.v0
    sensitivity.warm_start && (simulator.closed_orbit_guess .= closed_orbit)

    set_quad_factors!(simulator, vars_qh, vars_qv; tpsa=true)
    zero_tps = zero(first(simulator.phase_variables))
    phase_seed = vcat(simulator.phase_variables, zero_tps, zero_tps)
    one_turn = Bunch(
        closed_orbit + phase_seed',
        species=simulator.lattice.species_ref,
        p_over_q_ref=simulator.lattice.p_over_q_ref,
    )
    track!(one_turn, simulator.lattice)
    one_turn_jacobian =
        GTPSA.jacobian(one_turn.coords.v[1, :]; include_params=true)[1:4, :]
    R = @view one_turn_jacobian[:, 1:4]
    B = @view one_turn_jacobian[:, 5:end]
    closed_orbit_response = (I - R) \ B

    transverse_response = closed_orbit_response * simulator.quad_parameters
    parameterized_closed_orbit = closed_orbit +
        vcat(transverse_response, zero_tps, zero_tps)'
    response_bunch = Bunch(
        parameterized_closed_orbit,
        species=simulator.lattice.species_ref,
        p_over_q_ref=simulator.lattice.p_over_q_ref,
    )

    orbit_mm = zeros(N_BPMS)
    jacobian_mm = zeros(N_BPMS, N_QUADS)
    for (element_index, element) in enumerate(simulator.lattice.line)
        track!(response_bunch, element)
        haskey(simulator.bpm_by_index, element_index) || continue
        row, plane = simulator.bpm_by_index[element_index]
        orbit_mm[row] = 1e3 * scalar(response_bunch.coords.v[1, plane])
        jacobian_mm[row, :] .= 1e3 .* GTPSA.jacobian(
            response_bunch.coords.v[1, :]; include_params=true
        )[plane, 5:end]
    end
    return orbit_mm, jacobian_mm
end

"""
    predict_orbits_and_jacobian(simulator, currents_hc, currents_vc, quad_factors)

Evaluate all corrector settings. The Jacobian rows follow `vec(predictions)`.
"""
function predict_orbits_and_jacobian(
    simulator::BoosterSimulator,
    currents_hc::AbstractMatrix,
    currents_vc::AbstractMatrix,
    quad_factors::AbstractVector{<:Real},
)
    size(currents_hc, 2) == N_HC ||
        throw(ArgumentError("currents_hc must have $N_HC columns"))
    size(currents_vc, 2) == N_VC ||
        throw(ArgumentError("currents_vc must have $N_VC columns"))
    size(currents_hc, 1) == size(currents_vc, 1) ||
        throw(ArgumentError("current matrices must have the same number of rows"))
    length(quad_factors) == N_QUADS ||
        throw(ArgumentError("quad_factors must contain $N_QUADS values"))

    predictions = Matrix{Float64}(undef, size(currents_hc, 1), N_BPMS)
    jacobian = Matrix{Float64}(undef, length(predictions), N_QUADS)
    linear_indices = LinearIndices(predictions)
    vars_qh = @view quad_factors[1:N_QH]
    vars_qv = @view quad_factors[N_QH+1:N_QUADS]

    for setting in axes(currents_hc, 1)
        orbit, orbit_jacobian = orbit_and_quad_jacobian(
            simulator,
            @view(currents_hc[setting, :]),
            @view(currents_vc[setting, :]),
            vars_qh,
            vars_qv,
        )
        predictions[setting, :] .= orbit
        for bpm in axes(predictions, 2)
            jacobian[linear_indices[setting, bpm], :] .=
                @view orbit_jacobian[bpm, :]
        end
    end
    return predictions, jacobian
end

function predict_orbits(
    simulator::BoosterSimulator,
    currents_hc::AbstractMatrix,
    currents_vc::AbstractMatrix,
    quad_factors::AbstractVector{<:Real},
)
    predictions, _ = predict_orbits_and_jacobian(
        simulator, currents_hc, currents_vc, quad_factors
    )
    return predictions
end

function ChainRulesCore.rrule(
    ::typeof(predict_orbits),
    simulator::BoosterSimulator,
    currents_hc::AbstractMatrix,
    currents_vc::AbstractMatrix,
    quad_factors::AbstractVector{<:Real},
)
    predictions, jacobian = predict_orbits_and_jacobian(
        simulator, currents_hc, currents_vc, quad_factors
    )
    function pullback(prediction_tangent)
        prediction_tangent = unthunk(prediction_tangent)
        factor_tangent = prediction_tangent isa AbstractZero ?
            ZeroTangent() :
            ProjectTo(quad_factors)(jacobian' * vec(prediction_tangent))
        return NoTangent(), NoTangent(), NoTangent(), NoTangent(), factor_tangent
    end
    return predictions, pullback
end

# ReverseDiff needs an explicit registration to consume the ChainRules rule.
ReverseDiff.@grad_from_chainrules predict_orbits(
    simulator::BoosterSimulator,
    currents_hc::AbstractMatrix,
    currents_vc::AbstractMatrix,
    quad_factors::ReverseDiff.TrackedArray,
)
# DynamicPPL's product transform tracks the vector elements individually.
ReverseDiff.@grad_from_chainrules predict_orbits(
    simulator::BoosterSimulator,
    currents_hc::AbstractMatrix,
    currents_vc::AbstractMatrix,
    quad_factors::AbstractVector{<:ReverseDiff.TrackedReal},
)

struct BoosterInferenceData
    currents_hc::Matrix{Float64}
    currents_vc::Matrix{Float64}
    observed_orbit_mm::Matrix{Float64}
    noise_std_mm::Matrix{Float64}
end

function BoosterInferenceData(
    currents_hc::AbstractVector,
    currents_vc::AbstractVector,
    observed_orbit_mm::AbstractVector;
    noise_std_mm,
)
    return BoosterInferenceData(
        reshape(currents_hc, 1, :),
        reshape(currents_vc, 1, :),
        reshape(observed_orbit_mm, 1, :);
        noise_std_mm,
    )
end

function BoosterInferenceData(
    currents_hc::AbstractMatrix,
    currents_vc::AbstractMatrix,
    observed_orbit_mm::AbstractMatrix;
    noise_std_mm,
)
    hc = Matrix{Float64}(currents_hc)
    vc = Matrix{Float64}(currents_vc)
    observed = Matrix{Float64}(observed_orbit_mm)
    size(hc, 2) == N_HC || throw(ArgumentError("currents_hc has wrong width"))
    size(vc) == (size(hc, 1), N_VC) ||
        throw(ArgumentError("currents_vc has wrong size"))
    size(observed) == (size(hc, 1), N_BPMS) ||
        throw(ArgumentError("observed_orbit_mm has wrong size"))

    noise = noise_std_mm isa Real ?
        fill(Float64(noise_std_mm), size(observed)) :
        Matrix{Float64}(noise_std_mm)
    size(noise) == size(observed) ||
        throw(ArgumentError("noise_std_mm must be scalar or match observations"))
    all(isfinite, noise) && all(>(0), noise) ||
        throw(ArgumentError("noise_std_mm must be positive and finite"))
    return BoosterInferenceData(hc, vc, observed, noise)
end

function lognormal_quad_prior(log_std; center::Symbol=:median)
    stds = log_std isa Real ? fill(Float64(log_std), N_QUADS) : Float64.(log_std)
    length(stds) == N_QUADS ||
        throw(ArgumentError("log_std must contain $N_QUADS values"))
    all(isfinite, stds) && all(>(0), stds) ||
        throw(ArgumentError("log_std values must be positive and finite"))
    locations = center === :median ? zeros(N_QUADS) :
        center === :mean ? -0.5 .* abs2.(stds) :
        throw(ArgumentError("center must be :median or :mean"))
    return Turing.product_distribution(Turing.LogNormal.(locations, stds))
end

Turing.@model function booster_model(
    data::BoosterInferenceData,
    simulator::BoosterSimulator,
    quad_prior=lognormal_quad_prior(0.05),
)
    quad_factors ~ quad_prior
    predicted_orbit_mm = predict_orbits(
        simulator, data.currents_hc, data.currents_vc, quad_factors
    )
    for index in eachindex(data.observed_orbit_mm)
        Turing.@addlogprob! Turing.logpdf(
            Turing.Normal(predicted_orbit_mm[index], data.noise_std_mm[index]),
            data.observed_orbit_mm[index],
        )
    end
    return predicted_orbit_mm
end

end # module BoosterUQ

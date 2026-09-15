using Beamlines
include("booster_conversions.jl")
include("booster_lattice.jl")

booster = Beamline(
    BOOSTER;
    species_ref=Species("proton"),
    p_over_q_ref=reference_rigidity_expression(),
)


ele_index = Dict(getproperty.(booster.line, :name) .=> getproperty.(booster.line, :beamline_index))


for ele in vcat(booster.line)
    if ele.kind == "SBend"
        ele.x1_limit = -0.08
        ele.x2_limit =  0.08
        ele.y1_limit = -0.033
        ele.y2_limit =  0.033
        ele.aperture_shape = ApertureShape.Rectangular
    elseif (ele.name != "SPTMD3") && (ele.kind in ["Drift", "Quadrupole", "Sextupole"] || ele.name == "SPTMD6")
        ele.x1_limit = -0.0742
        ele.x2_limit =  0.0742
        ele.y1_limit = -0.0742
        ele.y2_limit =  0.0742
    end
end

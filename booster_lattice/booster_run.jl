using Beamlines
include("booster_conversions.jl")
include("booster_lattice.jl")

# Injection bump kicker calibrations (Booster Injection Kicker Data by Kip)
IKHC1  = 0.; IKHC3 = 0.; IKHC7 = 0.; IKHD1 = 0.
IKHCAL = 0.016 / 1200   # Bdl/I [Tm/A]
BDOT   = 0.

booster = Beamline(BOOSTER, species_ref=Species("proton"), p_over_q_ref=BRHO)


ele_index = Dict(getproperty.(booster.line, :name) .=> getproperty.(booster.line, :beamline_index))


for ele in vcat(booster.line)
    if ele.kind == "SBend"
        ele.Bn0 = BDIPO
        ele.Bn1 = B10
        ele.Bn2 = B20
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

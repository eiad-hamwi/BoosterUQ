# BoosterUQ

Julia packages needed: `Beamlines`, `BeamTracking`, and `SciBmad`

Example that calculates horizontal closed orbit at BPM positions:
```
using SciBmad

include("./booster_lattice/booster_run.jl")

CORRECTOR_CURRENTS[:DHCA2] = 17.4

t = twiss(booster)

h_bpm_indices = findall(element -> (element.kind == "BPM" && contains(element.name, "H")), booster.line)


#=
using Plots
plot(t.s, t.x)
scatter!(t.s[h_bpm_indices], t.x[h_bpm_indices])
=#
```

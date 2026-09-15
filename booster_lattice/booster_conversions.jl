using Beamlines

# Booster operating-point inputs and magnet transfer functions.
#
#  MAD file originally developed by M.Tanaka and A.Luccio for
#  Slow Extracted Beam from the Booster.
#
#  MODIFICATIONS APRIL 1992 BY MMB FOR NONLINEAR EFFECTS
#  
#  MODIFICATIONS oct 1992 BY AUL FOR NONLINEAR EFFECTS 
#  
#  MODIFICATIONS oct 2001 BY KAB FOR HIGH PRECISION TUNE PREDICTIONS
#  
#  10/01 - KAB

inch = 0.0254; # in meters

# %%
# Main magnet current (dipoles + quads)
IDIPO    = 1255.8 ;

# Tune quadrupoles
IQHC     =  209.2
IQVC     = -117.4;

# Chrom sextupoles
ISH      = 0.0  ;
ISV      = 0.0  ;

# Resonance drive sextupoles
ISEBC8F8 =  0. ;
ISEBB4E4 =  0. ;

# D3 extracion bump
IC7  =  0.;   
ID1  =  0.;  
ID4  =  0.;  
ID7  =  0.;
IE1  =  0.; 

# Fast kickers
IKHC1 = 0.;
IKHC3 = 0.;
IKHC7 = 0.;
IKHD1 = 0.;

# Others
BDOT     = 0.   ;
IQHBD    = 3.4 ;   # value based on 1992/3 data
IQVBD    = 4.8 ;   # value based on 1992/3 data
CBRHO    = 0.0 ;

D03kick  = 0.   ;

# Local tune trims used by the polarized-proton optics. 
const QV_TRIM_CURRENTS = Dict{Symbol,Float64}(
    :B5 => 0.0, :B7 => 0.0,
    :C1 => 0.0, :C3 => 0.0, :C5 => 0.0, :C7 => 0.0,
    :D1 => 0.0, :D3 => 0.0,
)
const QH_TRIM_CURRENTS = Dict{Symbol,Float64}(
    :B6 => 0.0, :B8 => 0.0,
    :C2 => 0.0, :C4 => 0.0, :C6 => 0.0, :C8 => 0.0,
    :D2 => 0.0, :D4 => 0.0,
)

# %%
# 1/2 integer stop band correctors
QVSTR1 = 0.;
QHSTR1 = 0.;
QVSTR2 = 0.;
QHSTR2 = 0.;

#---------------------------------------------------

# %%
# Dipole parameters
LEND = 2.42;
# ANGD  CONSTANT= 0.174533;
# a little more precise (it is 10 degrees, after all) - kab
ANGD = 0.174532925;
# RHO = LEND/ANGD;
RHO = 13.8656;  #Bend radius of Booster main dipole


# Dipole field from current, including saturation.
bdip(_IDIP) = (
  0.0009122 
  + _IDIP * ( 2.371E-4
  + _IDIP * ( 1.717E-8
  + _IDIP * (-2.412E-11 
  + _IDIP * ( 1.836e-14 
  + _IDIP * (-7.88e-18
  + _IDIP * ( 1.891e-21 
  + _IDIP * (-2.351e-25
  + _IDIP *   1.163e-29
         )))))))
)

brho(idipo, bdot=0) = bdip(idipo) * RHO - CBRHO * bdot

# Merged Horner form of K10 * BRHO.
function b10(_IDIP, bdot=0)
    _α = CBRHO * bdot / RHO
    return (
    -3.2574662e-6 + _α * 0.003571
    + _IDIP * (-8.48982844e-7 + _α * 2.52e-6
    + _IDIP * (-6.532006010e-10 - _α * 6.145e-9
    + _IDIP * (1.4933551414e-12 + _α * 7.113e-12
    + _IDIP * (-1.5816306318e-15 - _α * 4.531e-15
    + _IDIP * (7.843050448e-19 + _α * 1.666e-18
    + _IDIP * (-1.939736294e-23 - _α * 3.523e-22
    + _IDIP * (-2.37341900484e-25 + _α * 3.972e-26
    + _IDIP * (1.882277001034e-28 - _α * 1.847e-30
    + _IDIP * (-8.99582648e-32
    + _IDIP * (3.089792104e-35
    + _IDIP * (-7.84830113e-39
    + _IDIP * (1.45747595e-42
    + _IDIP * (-1.9186619e-46
    + _IDIP * (1.6928098e-50
    + _IDIP * (-8.961733e-55
    + _IDIP * (2.148061e-59
    ))))))))))))))))
    )
end

SBDOT = 0.;

b20(idipo, bdot=0) = (-0.4438 + 0.31764 * SBDOT * brho(idipo, bdot)) / RHO


function dipole_strengths(idipo, bdot=0)
    return (
        Bn0=bdip(idipo),
        Bn1=b10(idipo, bdot),
        Bn2=b20(idipo, bdot),
        p_over_q_ref=brho(idipo, bdot),
    )
end

function dipole_field_expressions()
    return (
        Bn0=DefExpr(() -> dipole_strengths(IDIPO, BDOT).Bn0),
        Bn1=DefExpr(() -> dipole_strengths(IDIPO, BDOT).Bn1),
        Bn2=DefExpr(() -> dipole_strengths(IDIPO, BDOT).Bn2),
    )
end

reference_rigidity_expression() =
    DefExpr(() -> dipole_strengths(IDIPO, BDOT).p_over_q_ref)


# %%
# Quadrupoles
# CK1: CONSTANT = 8.616694e-4;
CKC = 0.20022;    #quad trim coil vs. main coil
# To match change in power supply response to back emf due to bdot we add
# a bdot term, with a calibration coefficient
effective_qh_current(idipo, iqhc, bdot=0) = idipo + CKC * (iqhc + bdot * IQHBD)
effective_qv_current(idipo, iqvc, bdot=0) = idipo + CKC * (-iqvc + bdot * IQVBD)
LENQH = 0.493;
LENQV = 0.504;
LENACQ = 0.17;
# Adjust adjacent drifts to keep the circumference fixed.
LVUS = (0.504-LENQV)/2
LVDS = (0.504-LENQV)/2
LHUS = (0.493-LENQH)/2
LHDS = (0.493-LENQH)/2

# Fifth-order saturation fits to the measured integrated gradients.
b1lh(_IQH) = (
    0.001818 
    + _IQH * ( 9.080e-4 
    + _IQH * ( 6.657e-9 
    + _IQH * (-7.225e-12 
    + _IQH * ( 3.239e-15 
    - _IQH *   5.07e-19))))
)

cblh(iqh, idipo, bdot=0) = (1 - 0.00004179 * bdot / bdip(idipo)) * b1lh(iqh)

b1lv(_IQV) = (
    0.002099 
    + _IQV * ( 9.257e-4 
    + _IQV * ( 1.164e-8 
    + _IQV * (-1.046e-11 
    + _IQV * ( 4.057e-15 
    - _IQV *   5.75e-19))))
)

cblv(iqv, idipo, bdot=0; ear=false) = -(
    1 - (ear ? 0.000062913 : 0.000041942) * bdot / bdip(idipo)
) * 1.0030 * b1lv(iqv)


function quad_strengths(idipo, iqhc, iqvc, bdot=0)
    iqh = effective_qh_current(idipo, iqhc, bdot)
    iqv = effective_qv_current(idipo, iqvc, bdot)
    return (
        horizontal=cblh(iqh, idipo, bdot),
        vertical=cblv(iqv, idipo, bdot),
        vertical_ear=cblv(iqv, idipo, bdot; ear=true),
    )
end


# Four stop-band supplies drive alternating quadrupoles according to Booster
# TN149 and AGS TN465. Small trim coils have two turns.
const STOPBAND_CALIBRATION = 0.4

function stopband_current(plane::Symbol, label::Symbol)
    text = String(label)
    length(text) == 2 || throw(ArgumentError("invalid quadrupole label $label"))
    sector = Int(text[1]) - Int('A')
    position = parse(Int, text[2:2])
    0 <= sector <= 5 || throw(ArgumentError("invalid quadrupole label $label"))
    sector_sign = iseven(sector) ? 1 : -1

    if plane === :H
        position in (2, 8) && return sector_sign * QHSTR1
        position == 4 && return sector_sign * QHSTR2
        position == 6 && return -sector_sign * QHSTR2
    elseif plane === :V
        position in (1, 7) && return sector_sign * QVSTR1
        position == 3 && return sector_sign * QVSTR2
        position == 5 && return -sector_sign * QVSTR2
    else
        throw(ArgumentError("plane must be :H or :V"))
    end
    throw(ArgumentError("label $label does not belong to plane $plane"))
end

function quad_current(
    plane::Symbol,
    label::Symbol;
    idipo=IDIPO,
    iqhc=IQHC,
    iqvc=IQVC,
    bdot=BDOT,
)
    if plane === :H
        return effective_qh_current(idipo, iqhc, bdot) +
            CKC * get(QH_TRIM_CURRENTS, label, 0.0) +
            STOPBAND_CALIBRATION * stopband_current(plane, label)
    elseif plane === :V
        return effective_qv_current(idipo, iqvc, bdot) +
            CKC * get(QV_TRIM_CURRENTS, label, 0.0) +
            STOPBAND_CALIBRATION * stopband_current(plane, label)
    end
    throw(ArgumentError("plane must be :H or :V"))
end

function quad_integrated_gradient_value(plane::Symbol, label::Symbol)
    current = quad_current(plane, label)
    plane === :H && return cblh(current, IDIPO, BDOT)
    return cblv(current, IDIPO, BDOT; ear=label in (:D5, :F5))
end

quad_integrated_gradient(plane::Symbol, label::Symbol) =
    DefExpr(() -> quad_integrated_gradient_value(plane, label))


# AC quadrupole based on AGS AC skew quads.
IACQA4 = 0.;
# Sextupoles
LENS = 0.1;

const SEXTUPOLE_CALIBRATION = 0.013144

function chromatic_sextupole_integrated_field(plane::Symbol, current::Real)
    plane === :H && return SEXTUPOLE_CALIBRATION * current
    plane === :V && return -SEXTUPOLE_CALIBRATION * current
    throw(ArgumentError("plane must be :H or :V"))
end

function sextupole_integrated_field(
    name::Symbol,
    ish::Real,
    isv::Real,
    isebc8f8::Real=0,
    isebb4e4::Real=0,
)
    name in (:SHB4, :SHE4) &&
        return chromatic_sextupole_integrated_field(:H, ish - isebb4e4)
    name === :SHC8 &&
        return chromatic_sextupole_integrated_field(:H, ish + isebc8f8)
    name === :SHF8 &&
        return chromatic_sextupole_integrated_field(:H, ish - isebc8f8)

    text = String(name)
    startswith(text, "SH") &&
        return chromatic_sextupole_integrated_field(:H, ish)
    startswith(text, "SV") &&
        return chromatic_sextupole_integrated_field(:V, isv)
    throw(ArgumentError("invalid sextupole name $name"))
end

sextupole_integrated_field(name::Symbol) = DefExpr(
    () -> sextupole_integrated_field(
        name, ISH, ISV, ISEBC8F8, ISEBB4E4
    )
)

# Slow correctors
#   kab, 8/31/07
# Calib. of correctors from BoosterTN 224, R.Thern
#  Bdl/A = 9.75e-5 Tm/A
#  so a kick (in radians) is (1/Brho)*(Bdl/A)*amps
# 
CorCalib = 0.0000975;

const CORRECTOR_CURRENTS = let currents = Dict{Symbol,Float64}()
    for sector in 'A':'F', position in (2, 4, 6, 8)
        currents[Symbol("DHC", sector, position)] = 0.0
    end
    for sector in 'A':'F', position in (1, 3, 5, 7)
        currents[Symbol("DVC", sector, position)] = 0.0
    end
    currents
end

corrector_integrated_field(current::Real) = CorCalib * current
corrector_integrated_field(name::Symbol) =
    DefExpr(() -> corrector_integrated_field(CORRECTOR_CURRENTS[name]))

using Beamlines

# %% [markdown]
#  booster.conversions
#    parameter conversions and polynomials used in Booster lattice
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
#  Expected Variables:
#    1.  IDIPO               = for main dipole current
#    2.  BRHO                = used in tune and chromaticity calculations
#    3.  IQHC & IQVC         = tune trim power supplies
#    4.  ISH & ISV           = chrom. sextupoles
#    5.  ISEBC8F8 & ISEBB4E4 = drive sextupoles for 13/3 resonance
#    6.  D03kick             = D3 thin septum
#    7.  IC7, ID1,ID2, ID4, ID6, ID7, IE1  = for D3/D6 orbit bump
#             call file = BOOSTER.BMP to use these
#    8.  BDOT
#
#  Lines which can be used:
#  
#   BOOSTER
#   
#   B_D3 
#   
#   D3toD6
#  
#  10/01 - KAB

inch = 0.0254; # in meters

# %%
# ------------Basic Relationships---------------------
#  EE    = EK + E0;
#  PC    = sqrt(EE*EE - E0*E0);
#  GAMMA = EE/E0;
#  BETA  = sqrt(1.-1./(GAMMA*GAMMA));
# ---------------------------------------------------

# %%
# Main magnet current (dipoles + quads)
IDIPO    = 1680.0 ;

# Tune quadrupoles
IQHC     = -212.0;
IQVC     = -282.4;

# Chrom sextupoles
ISH      = 0.1  ;
ISV      = 6.3  ;

# Resonance drive sextupoles
ISEBC8F8 = -93. ;
ISEBB4E4 =  90. ;

# D3 extracion bump
IC7  = 127.94 #-32.54;   
ID1  =  39.61 # 61.49;  
ID4  = -87.25 #-21.02;  
ID7  =  79.58 # 49.81;
IE1  = -89.1 # -6.1; 

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
ANGD6    = 0.   ;

IQVTB5   = 0.;
IQVTB7   = 0.;
IQVTC1   = 0.;
IQVTC3   = 0.;
IQVTC5   = 0.;
IQVTC7   = 0.;
IQVTD1   = 0.;
IQVTD3   = 0.;

IQHTB6   = 0.;
IQHTB8   = 0.;
IQHTC2   = 0.;
IQHTC4   = 0.;
IQHTC6   = 0.;
IQHTC8   = 0.;
IQHTD2   = 0.;
IQHTD4   = 0.;

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


# %%
# Transfer function, dipole field from current.
# use expansion to include saturation
_BDIP(IDIP) = DefExpr(
  () -> begin
  _IDIP = IDIP()
  (
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
  end
);

BDIPO = _BDIP(() -> IDIPO)

BRHO  = DefExpr(() -> BDIPO()*RHO-CBRHO*BDOT);

# %%
# Bend
# K10 = -2.E-3/RHO;
# K10 = -3.918E-3/RHO;   #6/20/91 from R.Thern Tech Note 190
# 8th order fit to data in tech note 190
K10 = DefExpr(
  () -> ( -0.003571 
  + IDIPO * (-2.52e-6
  + IDIPO * ( 6.145e-9
  + IDIPO * (-7.113e-12
  + IDIPO * ( 4.531e-15
  + IDIPO * (-1.666e-18
  + IDIPO * (3.523e-22
  + IDIPO * (-3.972e-26
  + IDIPO * 1.847e-30))))))))/RHO);

# B10 = K10 * (BDIPO*RHO - CBRHO*BDOT); merged Horner: coef[n](IDIPO^n) uses (P_K10*BDIPO)_n - _α*(P_K10)_n for n≤8, else (P_K10*BDIPO)_n; _α = CBRHO*BDOT/RHO
_B10(IDIP) = let _α = CBRHO * BDOT / RHO;
  DefExpr(
  () -> begin 
  _IDIP = IDIP();
    (
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
  );
end

B10 = _B10(() -> IDIPO)


SBDOT = 0.;
# K20 = -0.4438/RHO;     #6/20/91 from R.Thern Tech Note 190

K20 = DefExpr(() ->(-0.4438+0.31764*SBDOT)/RHO);   #MMB 6/20/91 from R.Thern Tech Note 190 and SYL TN 147

_B20(IDIP) = let __BDIP = _BDIP(IDIP)
  DefExpr(
    () -> begin
      _BRHO = __BDIP() * RHO - CBRHO * BDOT
      (-0.4438 + 0.31764 * SBDOT * _BRHO) / RHO
    end
  );
end

B20 = _B20(() -> IDIPO);


# %%
# Quadrupoles
# CK1: CONSTANT = 8.616694e-4;
CKC = 0.20022;    #quad trim coil vs. main coil
# To match change in power supply response to back emf due to bdot we add
# a bdot term, with a calibration coefficient
IQH = DefExpr(() ->IDIPO + CKC*(IQHC + BDOT*IQHBD))
IQV = DefExpr(() ->IDIPO + CKC*(IQVC + BDOT*IQVBD))
# old lengths
LENQH = 0.493;
LENQV = 0.504;
# ac quad based on AGS AC Skew quads
LENACQ = 0.17;
# lengths corrected for saturation fit to values in E.Bleser Report
# LENQH = DefExpr(() ->1.0*(0.4998-IQH()*1.569e-05+IQH()^2*7.900e-09-IQH()^3*1.215e-12));
# LENQV = DefExpr(() ->1.0*(0.5125-IQV()*1.946e-05+IQV()^2*9.609e-09-IQV()^3*1.404e-12));
# lengths corrected using fit to Wuzeng Meng quadrupole models
# LENQH = DefExpr(() ->1.0*(0.4974+IQH()*4.863e-08+IQH()^2*1.539e-11-IQH()^3*1.166e-13));
# LENQV = DefExpr(() ->1.0*(0.5091-IQV()*9.765e-08+IQV()^2*1.049e-10-IQV()^3*1.336e-13));
# of course we need to adjust the drifts on either side to keep Circ. Constant
LVUS = (0.504-LENQV)/2
LVDS = (0.504-LENQV)/2
LHUS = (0.493-LENQH)/2
LHDS = (0.493-LENQH)/2

# %% [markdown]
# -------------------------------------------------------- 
#  MODIFICATIONS oct 1992 BY AUL FOR NONLINEAR EFFECTS
#  
#  quad integrated gradient in tesla. This is a linear aproximation valid
#  
#  in the range 600-800 amps. See file luccio/mad_spline.dat --- old stuff; kab
#
#  New stuff: 10/18/01 Kevin Brown
#  
#   added effective length as a function of current (see LENQH & V above)
#   
#    effective length was never measured, so we had to use Opera3D model to calculate
#    
#   added 5th order polynomial to include saturation effects.
#   
#    for shorts this was fudged by 0.9996 to match real data
#    
#    for longs this was fudged by 1.0024 to match real data
#    
#   added BDot dependence based on Opera2D transient model, adjusted to match
#   
#   real data.
#
#   short quads ( current in amps )
#   
#  5th order, match to field measurements

# %%
B1LH(IQH)   = DefExpr(
  () -> begin
    _IQH = IQH()
    (
    0.001818 
    + _IQH * ( 9.080e-4 
    + _IQH * ( 6.657e-9 
    + _IQH * (-7.225e-12 
    + _IQH * ( 3.239e-15 
    - _IQH *   5.07e-19))))
    )
  end
  );

# for 2D transient model, coeff. is 5e-5
# CKH  = DefExpr(() ->(1-0.00005*BDOT/BDIPO())*1.00*B1LH()/(BRHO*LENQH);
# for 3D transient model
CBLH(IQH) = DefExpr(() -> (1-0.00004179*(BDOT/BDIPO()))*1.00*B1LH(IQH)());
# fit to experimental data, Note sign change 
# 0.00013;
# CKH = (1+0.0001313*BDOT/BDIPO())*1.00*B1LH()/(BRHO*LENQH);

#  long quads ( current in amps )
# 5th order, match to field measurements
B1LV(IQV)   = DefExpr(
  () -> begin
    _IQV = IQV()
    (
    0.002099 
    + _IQV * ( 9.257e-4 
    + _IQV * ( 1.164e-8 
    + _IQV * (-1.046e-11 
    + _IQV * ( 4.057e-15 
    - _IQV *   5.75e-19))))
    )
  end
  );

# for 2D transient model, coeff. is 5e-5
# CKV   = DefExpr(() -> -(1-0.00005*BDOT/BDIPO())*1.00*B1LV()/(BRHO*LENQV));
# for 3D transient model
CBLV(IQV) = DefExpr(() -> -(1-0.000041942*(BDOT/BDIPO()))*1.0030*B1LV(IQV)());
CBLVear(IQV) = DefExpr(() -> -(1-0.000062913*(BDOT/BDIPO()))*1.0030*B1LV(IQV)());


# %% [markdown]
#  Stop band correctors are configured as follows per Booster TN149 and AGS TN 465
#  
#   4 power supplies for 1/2 integer to 2 turns on main quadrupoles
#   
#   p.s. 1 QVSTR1 = +A1, +A7, -B1, -B7, ...
#   
#   p.s. 2 QHSTR1 = +A2, +A8, -B2, -B8, ...
#   
#   p.s. 3 QVSTR2 = +A3, -A5, -B3, +B5, ...
#   
#   p.s. 4 QHSTR2 = +A4, -A6, -B4, +B6, ...
#
#  IQVT and IQHT are trims for pol. proton test, small beta at C5 foil 7/6/07 KAB
#
#  small trim coils are 2 turns

# %%
ckc2 = 0.4;

IQVA1   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*( QVSTR1));
IQVA3   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*( QVSTR2));
IQVA5   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*(-QVSTR2));
IQVA7   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*( QVSTR1));
IQVB1   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*(-QVSTR1));
IQVB3   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*(-QVSTR2));
IQVB5   = DefExpr(() -> IDIPO + CKC*(-IQVC + IQVTB5 + BDOT*IQVBD) + ckc2*( QVSTR2));
IQVB7   = DefExpr(() -> IDIPO + CKC*(-IQVC + IQVTB7 + BDOT*IQVBD) + ckc2*(-QVSTR1));
IQVC1   = DefExpr(() -> IDIPO + CKC*(-IQVC + IQVTC1 + BDOT*IQVBD) + ckc2*( QVSTR1));
IQVC3   = DefExpr(() -> IDIPO + CKC*(-IQVC + IQVTC3 + BDOT*IQVBD) + ckc2*( QVSTR2));
IQVC5   = DefExpr(() -> IDIPO + CKC*(-IQVC + IQVTC5 + BDOT*IQVBD) + ckc2*(-QVSTR2));
IQVC7   = DefExpr(() -> IDIPO + CKC*(-IQVC + IQVTC7 + BDOT*IQVBD) + ckc2*( QVSTR1));
IQVD1   = DefExpr(() -> IDIPO + CKC*(-IQVC + IQVTD1 + BDOT*IQVBD) + ckc2*(-QVSTR1));
IQVD3   = DefExpr(() -> IDIPO + CKC*(-IQVC + IQVTD3 + BDOT*IQVBD) + ckc2*(-QVSTR2));
IQVD5   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*( QVSTR2));
IQVD7   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*(-QVSTR1));
IQVE1   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*( QVSTR1));
IQVE3   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*( QVSTR2));
IQVE5   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*(-QVSTR2));
IQVE7   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*( QVSTR1));
IQVF1   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*(-QVSTR1));
IQVF3   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*(-QVSTR2));
IQVF5   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*( QVSTR2));
IQVF7   = DefExpr(() -> IDIPO + CKC*(-IQVC + BDOT*IQVBD) + ckc2*(-QVSTR1));

IQHA2   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*( QHSTR1));
IQHA4   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*( QHSTR2));
IQHA6   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*(-QHSTR2));
IQHA8   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*( QHSTR1));
IQHB2   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*(-QHSTR1));
IQHB4   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*(-QHSTR2));
IQHB6   = DefExpr(() -> IDIPO + CKC*(IQHC + IQHTB6 + BDOT*IQHBD) + ckc2*( QHSTR2));
IQHB8   = DefExpr(() -> IDIPO + CKC*(IQHC + IQHTB8 + BDOT*IQHBD) + ckc2*(-QHSTR1));
IQHC2   = DefExpr(() -> IDIPO + CKC*(IQHC + IQHTC2 + BDOT*IQHBD) + ckc2*( QHSTR1));
IQHC4   = DefExpr(() -> IDIPO + CKC*(IQHC + IQHTC4 + BDOT*IQHBD) + ckc2*( QHSTR2));
IQHC6   = DefExpr(() -> IDIPO + CKC*(IQHC + IQHTC6 + BDOT*IQHBD) + ckc2*(-QHSTR2));
IQHC8   = DefExpr(() -> IDIPO + CKC*(IQHC + IQHTC8 + BDOT*IQHBD) + ckc2*( QHSTR1));
IQHD2   = DefExpr(() -> IDIPO + CKC*(IQHC + IQHTD2 + BDOT*IQHBD) + ckc2*(-QHSTR1));
IQHD4   = DefExpr(() -> IDIPO + CKC*(IQHC + IQHTD4 + BDOT*IQHBD) + ckc2*(-QHSTR2));
IQHD6   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*( QHSTR2));
IQHD8   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*(-QHSTR1));
IQHE2   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*( QHSTR1));
IQHE4   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*( QHSTR2));
IQHE6   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*(-QHSTR2));
IQHE8   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*( QHSTR1));
IQHF2   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*(-QHSTR1));
IQHF4   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*(-QHSTR2));
IQHF6   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*( QHSTR2));
IQHF8   = DefExpr(() -> IDIPO + CKC*(IQHC + BDOT*IQHBD) + ckc2*(-QHSTR1));


# %%
# for ac quadrupole
IACQA4 = 0.;
# ACQA4: QUADRUPOLE, L:=LENACQ, K1:=K1ACQA4;



# %%

# %%
# fit to experimental data, Note sign change 
# 0.00019
# CKV   := -(1+0.0001318*(BDOT/BDIPO))*1.0030*B1LV/(BRHO*LENQV);
# CKVear := -(1+0.00095*(BDOT/BDIPO))*1.0030*B1LV/(BRHO*LENQV);
# -------------------------------------------------------- 
#  
# Sextupoles
# 
LENS = 0.1;
#LENS: CONSTANT = 0.0631;         #E.Bleser Tech Note # 182
SCON = DefExpr(() -> 0.013144/(LENS*BRHO()));    #E.Bleser Tech Note # 182
K2SH = DefExpr(() -> SCON());
K2SV = DefExpr(() -> -SCON());



# %%
# Slow correctors
#   kab, 8/31/07
# Calib. of correctors from BoosterTN 224, R.Thern
#  Bdl/A = 9.75e-5 Tm/A
#  so a kick (in radians) is (1/Brho)*(Bdl/A)*amps
# 
CorCalib = 0.0000975;
#  currents are called as DHCA2i, etc, initialized to 0
DHCA2i=0.;
DHCA4i=0.;
DHCA6i=0.;
DHCA8i=0.;
DHCB2i=0.;
DHCB4i=0.;
DHCB6i=0.;
DHCB8i=0.;
DHCC2i=0.;
DHCC4i=0.;
DHCC6i=0.;
DHCC8i=0.;
DHCD2i=0.;
DHCD4i=0.;
DHCD6i=0.;
DHCD8i=0.;
DHCE2i=0.;
DHCE4i=0.;
DHCE6i=0.;
DHCE8i=0.;
DHCF2i=0.;
DHCF4i=0.;
DHCF6i=0.;
DHCF8i=0.;

DVCA1i=0.;
DVCA3i=0.;
DVCA5i=0.;
DVCA7i=0.;
DVCB1i=0.;
DVCB3i=0.;
DVCB5i=0.;
DVCB7i=0.;
DVCC1i=0.;
DVCC3i=0.;
DVCC5i=0.;
DVCC7i=0.;
DVCD1i=0.;
DVCD3i=0.;
DVCD5i=0.;
DVCD7i=0.;
DVCE1i=0.;
DVCE3i=0.;
DVCE5i=0.;
DVCE7i=0.;
DVCF1i=0.;
DVCF3i=0.;
DVCF5i=0.;
DVCF7i=0.;

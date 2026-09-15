using Beamlines


@elements begin
# %%
#Drift parameters
L007 =  Drift(L=0.069600);
# L011 = Drift(L=0.111825);  # these are now defined based on EffL of quads
# L012 = Drift(L=0.117325);  # see quadrupoles section below
L014 =   Drift(L=0.138050);
L019 =   Drift(L=0.186080);
L0218 =  Drift(L=0.217555);  #7/26/91 for new position of SPTMC3
L028 =   Drift(L=0.276675);
L02835 = Drift(L=0.2835);
L02869 = Drift(L=0.2869);
# L029 =   Drift(L=0.289875);  # these are now defined based on EffL of quads
# L031 =   Drift(L=0.295375);  # see quadrupoles section below
L030 =   Drift(L=0.293725);
L032 =   Drift(L=0.305784);
L03025 = Drift(L=0.3025);
L057 =   Drift(L=0.570400);
L0605 =  Drift(L=0.60500);
L222 =   Drift(L=2.218011);
L271 =   Drift(L=2.715375);
L2757 =  Drift(L=2.756936);   #7/26/91 for new position of SPTMC3
L328 =   Drift(L=3.280275);

# %%
# Dipole parameters
LDH =  Drift(L=LEND);
LDHH = Drift(L=LEND/2);
# ionization profile monitors
LDH1 = Drift(L=1.5 - (34.89inch-0.28));
LDH2 = Drift(L=0.32);
LDH3 = Drift(L=0.32);
LDH4 = Drift(); # 0.28); # LDH1 + LDH2 + LDH3 + LDH4 == LEND
IPMV =  Marker();  # VERTICAL IPM
IPMSK = Marker();  # SKEW FIELD CORRECTOR
IPMH =  Marker();  # HORIZONTAL IPM

# %%
# Bend

DIPOLE_FIELDS = dipole_field_expressions();

DHA1 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHA2 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHA4 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHA5 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHA7 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHA8 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);

DHB1 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHB2 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHB4 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHB5 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHB7 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHB8 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);

# a study of a broken vacuum chamber coil - Booster History 1992
C7BMPED=1.0;
DHC1 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHC2 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHC4 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHC5 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHC7 = SBend(L=LEND, angle=ANGD*C7BMPED; DIPOLE_FIELDS...);
DHC8 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);

DHD1 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHD2 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHD4 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHD5 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHD7 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHD8 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);

DHE1 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHE2 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHE4 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHE5 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHE7 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHE8 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);

DHF1 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHF2 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHF4 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHF5 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHF7 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
DHF8 = SBend(L=LEND, angle=ANGD; DIPOLE_FIELDS...);
end

# %%
# Bumps for D3 and D6 using thin kicks in middle of magnets
# bumpcal=4.215e-4;
bumpcal=1.6*4.215e-4;

@elements begin
# %%
# Quadrupoles

L011=Drift(L=0.111825+LVUS);
L012=Drift(L=0.117325+LHUS);
L029=Drift(L=0.289875+LVDS);
L031=Drift(L=0.295375+LHDS);
L031s=Drift(L=(0.295375+LHDS-LENACQ)/2.);
 
QVA1 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :A1));
QVA3 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :A3));
QVA5 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :A5));
QVA7 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :A7));

QHA2 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :A2));
QHA4 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :A4));
QHA6 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :A6));
QHA8 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :A8));

QVB1 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :B1));
QVB3 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :B3));
QVB5 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :B5));
QVB7 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :B7));

QHB2 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :B2));
QHB4 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :B4));
QHB6 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :B6));
QHB8 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :B8));

QVC1 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :C1));
QVC3 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :C3));
QVC5 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :C5));
QVC7 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :C7));

QHC2 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :C2));
QHC4 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :C4));
QHC6 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :C6));
QHC8 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :C8));

QVD1 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :D1));
QVD3 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :D3));
QVD5 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :D5)); # ear chamber
QVD7 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :D7));

QHD2 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :D2));
QHD4 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :D4));
QHD6 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :D6));
QHD8 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :D8));

QVE1 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :E1));
QVE3 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :E3));
QVE5 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :E5));
QVE7 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :E7));

QHE2 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :E2));
QHE4 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :E4));
QHE6 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :E6));
QHE8 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :E8));

QVF1 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :F1));
QVF3 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :F3));
QVF5 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :F5)); # ear chamber
QVF7 = Quadrupole(L=LENQV, Bn1L=quad_integrated_gradient(:V, :F7));

QHF2 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :F2));
QHF4 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :F4));
QHF6 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :F6));
QHF8 = Quadrupole(L=LENQH, Bn1L=quad_integrated_gradient(:H, :F8));

# AC Quad
ACQA4= Quadrupole(L=LENACQ, Bn1L=DefExpr(() -> 6.7e-3 * IACQA4 * LENQV));
end
# %%
# Sextupoles
SV_Bn2L=DefExpr(() -> chromatic_sextupole_integrated_field(:V, ISV)) # -0.013144 * ISV
SH_Bn2L=DefExpr(() -> chromatic_sextupole_integrated_field(:H, ISH)) #  0.013144 * ISH

@elements begin
SVA1= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHA2= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVA3= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHA4= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVA5= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHA6= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVA7= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHA8= Sextupole(L=LENS, Bn2L=SH_Bn2L);

SVB1= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHB2= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVB3= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHB4= Sextupole(L=LENS, Bn2L=sextupole_integrated_field(:SHB4));  # for 13/3 extraction
SVB5= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHB6= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVB7= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHB8= Sextupole(L=LENS, Bn2L=SH_Bn2L);

SVC1= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHC2= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVC3= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHC4= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVC5= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHC6= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVC7= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHC8= Sextupole(L=LENS, Bn2L=sextupole_integrated_field(:SHC8));

SVD1= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHD2= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVD3= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHD4= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVD5= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHD6= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVD7= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHD8= Sextupole(L=LENS, Bn2L=SH_Bn2L);

SVE1= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHE2= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVE3= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHE4= Sextupole(L=LENS, Bn2L=sextupole_integrated_field(:SHE4));
SVE5= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHE6= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVE7= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHE8= Sextupole(L=LENS, Bn2L=SH_Bn2L);

SVF1= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHF2= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVF3= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHF4= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVF5= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHF6= Sextupole(L=LENS, Bn2L=SH_Bn2L);
SVF7= Sextupole(L=LENS, Bn2L=SV_Bn2L);
SHF8= Sextupole(L=LENS, Bn2L=sextupole_integrated_field(:SHF8));   # for 13/3 extraction

# %%
# Slow correctors
#   kab, 8/31/07
# Calib. of correctors from BoosterTN 224, R.Thern
#  Bdl/A = 9.75e-5 Tm/A
#  so a kick (in radians) is (1/Brho)*(Bdl/A)*amps
 
DVCA1=VKicker(Bs0L=corrector_integrated_field(:DVCA1));
DVCA3=VKicker(Bs0L=corrector_integrated_field(:DVCA3));
DVCA5=VKicker(Bs0L=corrector_integrated_field(:DVCA5));
DVCA7=VKicker(Bs0L=corrector_integrated_field(:DVCA7));

DHCA2=HKicker(Bn0L=corrector_integrated_field(:DHCA2));
DHCA4=HKicker(Bn0L=corrector_integrated_field(:DHCA4));
DHCA6=HKicker(Bn0L=corrector_integrated_field(:DHCA6));
DHCA8=HKicker(Bn0L=corrector_integrated_field(:DHCA8));

DVCB1=VKicker(Bs0L=corrector_integrated_field(:DVCB1));
DVCB3=VKicker(Bs0L=corrector_integrated_field(:DVCB3));
DVCB5=VKicker(Bs0L=corrector_integrated_field(:DVCB5));
DVCB7=VKicker(Bs0L=corrector_integrated_field(:DVCB7));

DHCB2=HKicker(Bn0L=corrector_integrated_field(:DHCB2));
DHCB4=HKicker(Bn0L=corrector_integrated_field(:DHCB4));
DHCB6=HKicker(Bn0L=corrector_integrated_field(:DHCB6));
DHCB8=HKicker(Bn0L=corrector_integrated_field(:DHCB8));

DVCC1=VKicker(Bs0L=corrector_integrated_field(:DVCC1));
DVCC3=VKicker(Bs0L=corrector_integrated_field(:DVCC3));
DVCC5=VKicker(Bs0L=corrector_integrated_field(:DVCC5));
DVCC7=VKicker(Bs0L=corrector_integrated_field(:DVCC7));

DHCC2=HKicker(Bn0L=corrector_integrated_field(:DHCC2));
DHCC4=HKicker(Bn0L=corrector_integrated_field(:DHCC4));
DHCC6=HKicker(Bn0L=corrector_integrated_field(:DHCC6));
DHCC8=HKicker(Bn0L=corrector_integrated_field(:DHCC8));

DVCD1=VKicker(Bs0L=corrector_integrated_field(:DVCD1));
DVCD3=VKicker(Bs0L=corrector_integrated_field(:DVCD3));
DVCD5=VKicker(Bs0L=corrector_integrated_field(:DVCD5));
DVCD7=VKicker(Bs0L=corrector_integrated_field(:DVCD7));

DHCD2=HKicker(Bn0L=corrector_integrated_field(:DHCD2));
DHCD4=HKicker(Bn0L=corrector_integrated_field(:DHCD4));
DHCD6=HKicker();                                    #missing
DHCD8=HKicker(Bn0L=corrector_integrated_field(:DHCD8));

DVCE1=VKicker(Bs0L=corrector_integrated_field(:DVCE1));
DVCE3=VKicker(Bs0L=corrector_integrated_field(:DVCE3));
DVCE5=VKicker(Bs0L=corrector_integrated_field(:DVCE5));
DVCE7=VKicker(Bs0L=corrector_integrated_field(:DVCE7));

DHCE2=HKicker(Bn0L=corrector_integrated_field(:DHCE2));
DHCE4=HKicker(Bn0L=corrector_integrated_field(:DHCE4));
DHCE6=HKicker(Bn0L=corrector_integrated_field(:DHCE6));
DHCE8=HKicker(Bn0L=corrector_integrated_field(:DHCE8));

DVCF1=VKicker(Bs0L=corrector_integrated_field(:DVCF1));
DVCF3=VKicker(Bs0L=corrector_integrated_field(:DVCF3));
DVCF5=VKicker(Bs0L=corrector_integrated_field(:DVCF5));
DVCF7=VKicker(Bs0L=corrector_integrated_field(:DVCF7));

DHCF2=HKicker(Bn0L=corrector_integrated_field(:DHCF2));
DHCF4=HKicker(Bn0L=corrector_integrated_field(:DHCF4));
DHCF6=HKicker();                                    #missing
DHCF8=HKicker(Bn0L=corrector_integrated_field(:DHCF8));

# %%
# Beam position monitors

PUEVA1 = LineElement(kind="BPM");
PUEHA2 = LineElement(kind="BPM");
PUEVA3 = LineElement(kind="BPM");
PUEHA4 = LineElement(kind="BPM");
PUEVA5 = LineElement(kind="BPM");
PUEHA6 = LineElement(kind="BPM");
PUEVA7 = LineElement(kind="BPM");
PUEHA8 = LineElement(kind="BPM");

PUEVB1 = LineElement(kind="BPM");
PUEHB2 = LineElement(kind="BPM");
PUEVB3 = LineElement(kind="BPM");
PUEHB4 = LineElement(kind="BPM");
PUEVB5 = LineElement(kind="BPM");
PUEHB6 = LineElement(kind="BPM");
PUEVB7 = LineElement(kind="BPM");
PUEHB8 = LineElement(kind="BPM");

PUEVC1 = LineElement(kind="BPM");
PUEHC2 = LineElement(kind="BPM");
PUEVC3 = LineElement(kind="BPM");
PUEHC4 = LineElement(kind="BPM");
PUEVC5 = LineElement(kind="BPM");
PUEHC6 = LineElement(kind="BPM");
PUEVC7 = LineElement(kind="BPM");
PUEHC8 = LineElement(kind="BPM");

PUEVD1 = LineElement(kind="BPM");
PUEHD2 = LineElement(kind="BPM");
PUEVD3 = LineElement(kind="BPM");
PUEHD4 = LineElement(kind="BPM");
PUEVD5 = LineElement(kind="BPM");
PUEHD6 = LineElement(kind="BPM");      #missing
PUEVD7 = LineElement(kind="BPM");
PUEHD8 = LineElement(kind="BPM");

PUEVE1 = LineElement(kind="BPM");
PUEHE2 = LineElement(kind="BPM");
PUEVE3 = LineElement(kind="BPM");
PUEHE4 = LineElement(kind="BPM");
PUEVE5 = LineElement(kind="BPM");
PUEHE6 = LineElement(kind="BPM");
PUEVE7 = LineElement(kind="BPM");
PUEHE8 = LineElement(kind="BPM");

PUEVF1 = LineElement(kind="BPM");
PUEHF2 = LineElement(kind="BPM");
PUEVF3 = LineElement(kind="BPM");
PUEHF4 = LineElement(kind="BPM");
PUEVF5 = LineElement(kind="BPM");
PUEHF6 = LineElement(kind="BPM");       #missing
PUEVF7 = LineElement(kind="BPM");
PUEHF8 = LineElement(kind="BPM");

MONH= Marker();
MONV= Marker();

CAVITY= Marker(); #RFCavity(voltage=0.00, harmon=3);

IJFOIL= Marker();     # Proton Injection Foil
SPTMC3= Marker();     # Heavy Ion Injection Septum

START_OF_BEXT= Marker();  # This is just upstream of F2 main mag.
START_OF_BTA = Marker();  # This is just upstream of F6 SEPTUM (SEPTF6)
end


# %%
# Proton Injection Fast Kickers
# Heavy Ion Injection Fast Kickers 
IKHCAL = 0.016 / 1200;

@elements begin
IJKDHC1=HKicker(Bn0L= DefExpr(() -> IKHCAL * IKHC1)); #         T7_KHC1 +T16_KHC1;
IJKDHC3=HKicker(Bn0L= DefExpr(() -> IKHCAL * IKHC3)); #T2_KHC3 +T7_KHC3 +T16_KHC3;
IJKDHC7=HKicker(Bn0L= DefExpr(() -> IKHCAL * IKHC7)); #T2_KHC7 +T7_KHC7 +T16_KHC7;
IJKDHD1=HKicker(Bn0L= DefExpr(() -> IKHCAL * IKHD1)); #T2_KHD1 +T7_KHD1 +T16_KHD1;
end
#HKICKER, IJKDHC3, KICK= 0.716/161.29;    //1/17/92 to match experiment
#HKICKER, IJKDHD1, KICK= 0.380/161.29;    //1/17/92 to match experiment

# Extraction Kickers

F3Kick = 0.0;

@elements begin
X1DHF3=HKicker(Kn0L= DefExpr(() -> F3Kick));
X2DHF3=HKicker(Kn0L= DefExpr(() -> F3Kick));
X3DHF3=HKicker(Kn0L= DefExpr(() -> F3Kick));
X4DHF3=HKicker(Kn0L= DefExpr(() -> F3Kick));
SPTMD3=HKicker(  L = 34.89inch,
                 Kn0L= DefExpr(() -> D03kick),
                 y1_limit=-0.096, y2_limit=0.096, 
                 x1_limit=-0.096, x2_limit=0.05, 
                 aperture_shape = ApertureShape.Rectangular);
end
# %% [markdown]
#     probably should use an RBEND for the beamline def
#     beam enters normal to face of magnet, well not exactly
#     to properly model, the Booster equilibrium orbit is distorted and beam enters
#     D3 at some angle. That angle will come from the orbit distortion and hopefully
#     doesn't change much from setup  to setup = however getting that number could
#     be useful for a more realistic simulation
#     The D3 bend is nominally 3 mRad

# %%
L_D3S=0.67;


# the D6 septum
@elements begin

SPTMF6= Marker();     # Extraction Septum
SPTMB6= Marker();     # SEB Septum
SPTMD6= Marker(x1_limit = -3.0625inch, x2_limit = 3.0625inch);     # SEB Septum

# collimator for SEB (AUL & ST 12/23/92)
COLLB6= Marker();
end
# %%
# Special Sections

LJC1 = [L028, IJKDHC1, L030];
#LJC3 = [L032, IJKDHC3, L222, SPTMC3, L019];     #Heavy Ion Injection
LJC3 = [L032, IJKDHC3, L2757, SPTMC3, L0218];     #7/26/91 Heavy Ion Injection
LJC7 = [L028, IJKDHC7, L030];
LJD1 = [L028, IJKDHD1, L030];
LIJF = [L02835, IJFOIL, L02869];                #Proton Injection
RFA3 = [LDHH, CAVITY, LDHH];                    #Proton cavity


# %%
#Extraction fast kicker
 
S0= 0.130;   # distance along missing dipole to upstream edge of 1st kicker
@elements begin
LEK1=Drift(L=S0+0.2731);
LEK2=Drift(L=0.6096);
LEK3=Drift(L=0.5004);
LEK4=Drift(L=0.5391);
LEK5=Drift(L=0.4978-S0);
end
LXF3 = [LEK1,X1DHF3,LEK2,X2DHF3,LEK3,X3DHF3,LEK4,X4DHF3,LEK5];
LXF6 = [START_OF_BTA, SPTMF6, LDH];                           #Extraction septum

# %%
# Sublines


LC1 = vcat(LJC1, [DVCC1, L007, SVC1, L014, PUEVC1, L011, QVC1, L029, DHC1])
LC2 = [L057, DHCC2, L007, SHC2, L014, PUEHC2, L012, QHC2, L031, DHC2]
LC3 = vcat([L057, DVCC3, L007, SVC3, L014, PUEVC3, L011, QVC3], LJC3);
LC4 = [      DHCC4, L007, SHC4, L014, PUEHC4, L012, QHC4, L031, DHC4]   #7/26/91
LC5 = [L057, DVCC5, L007, SVC5, L014, PUEVC5, L011, QVC5, L029, DHC5]
LC6 = vcat(LIJF, [DHCC6, L007, SHC6, L014, PUEHC6, L012, QHC6, L031, LDH])
LC7 = vcat(LJC7, [DVCC7, L007, SVC7, L014, PUEVC7, L011, QVC7, L029, DHC7])
LC8 = [L057, DHCC8, L007, SHC8, L014, PUEHC8, L012, QHC8, L031, DHC8]


LD1 = vcat(LJD1, [DVCD1, L007, SVD1, L014, PUEVD1, L011, QVD1, L029, DHD1]);
LD2 =     [L057,  DHCD2, L007, SHD2, L014, PUEHD2, L012, QHD2, L031, DHD2]

LIPM  = [LDH1, IPMV, LDH2, IPMSK, LDH3, IPMH, SPTMD3, LDH4]; 
LD3  = vcat([L057, DVCD3, L007, SVD3, L014, PUEVD3, L011, QVD3, L029], LIPM);

LD4 = [L057, DHCD4, L007, SHD4, L014, PUEHD4, L012, QHD4, L031, DHD4]
LD5 = [L057, DVCD5, L007, SVD5, L014, PUEVD5, L011, QVD5, L029, DHD5]

LD6     =[L057,DHCD6,L007,SHD6,L014,PUEHD6,L012,QHD6,L031,SPTMD6,LDH]

LD7 = [L057, DVCD7, L007, SVD7, L014, PUEVD7, L011, QVD7, L029, DHD7]
LD8 = [L057, DHCD8, L007, SHD8, L014, PUEHD8, L012, QHD8, L031, DHD8]


LE1 = [L057, DVCE1, L007, SVE1, L014, PUEVE1, L011, QVE1, L029, DHE1]
LE2 = [L057, DHCE2, L007, SHE2, L014, PUEHE2, L012, QHE2, L031, DHE2]
LE3 = [L057, DVCE3, L007, SVE3, L014, PUEVE3, L011, QVE3, L029, LDH]
LE4 = [L057, DHCE4, L007, SHE4, L014, PUEHE4, L012, QHE4, L031, DHE4] 
LE5 = [L057, DVCE5, L007, SVE5, L014, PUEVE5, L011, QVE5, L029, DHE5]
LE6 = [L057, DHCE6, L007, SHE6, L014, PUEHE6, L012, QHE6, L031, LDH]
LE7 = [L057, DVCE7, L007, SVE7, L014, PUEVE7, L011, QVE7, L029, DHE7]
LE8 = [L057, DHCE8, L007, SHE8, L014, PUEHE8, L012, QHE8, L031, DHE8]


LF1 = [L057, DVCF1, L007, SVF1, L014, PUEVF1, L011, QVF1, L029, DHF1]
LF2 = [L057, DHCF2, L007, SHF2, L014, PUEHF2, L012, QHF2, START_OF_BEXT, L031, DHF2]
LF3 = vcat([L057, DVCF3, L007, SVF3, L014, PUEVF3, L011, QVF3, L029], LXF3);
LF4 = [L057, DHCF4, L007, SHF4, L014, PUEHF4, L012, QHF4, L031, DHF4]
LF5 = [L057, DVCF5, L007, SVF5, L014, PUEVF5, L011, QVF5, L029, DHF5]
LF6 = vcat([L057, DHCF6, L007, SHF6, L014, PUEHF6, L012, QHF6, L031], LXF6);
LF7 = [L057, DVCF7, L007, SVF7, L014, PUEVF7, L011, QVF7, L029, DHF7]
LF8 = [L057, DHCF8, L007, SHF8, L014, PUEHF8, L012, QHF8, L031, DHF8]


LA1 = [L057, DVCA1, L007, SVA1, L014, PUEVA1, L011, QVA1, L029, DHA1]
LA2 = [L057, DHCA2, L007, SHA2, L014, PUEHA2, L012, QHA2, L031, DHA2]
LA3 = vcat([L057, DVCA3, L007, SVA3, L014, PUEVA3, L011, QVA3, L029], RFA3)
LA4 = [L057, DHCA4, L007, SHA4, L014, PUEHA4, L012, QHA4, L031s, ACQA4, L031s, DHA4]
LA5 = [L057, DVCA5, L007, SVA5, L014, PUEVA5, L011, QVA5, L029, DHA5]
LA6 = [L057, DHCA6, L007, SHA6, L014, PUEHA6, L012, QHA6, L031, LDH]
LA7 = [L057, DVCA7, L007, SVA7, L014, PUEVA7, L011, QVA7, L029, DHA7]
LA8 = [L057, DHCA8, L007, SHA8, L014, PUEHA8, L012, QHA8, L031, DHA8]


LB1 = [L057, DVCB1, L007, SVB1, L014, PUEVB1, L011, QVB1, L029, DHB1]
LB2 = [L057, DHCB2, L007, SHB2, L014, PUEHB2, L012, QHB2, L031, DHB2]
LB3 = [L057, DVCB3, L007, SVB3, L014, PUEVB3, L011, QVB3, L029, LDH];
LB4 = [L057, DHCB4, L007, SHB4, L014, PUEHB4, L012, QHB4, L031, DHB4] 
LB5 = [L057, DVCB5, L007, SVB5, L014, PUEVB5, L011, QVB5, L029, DHB5]
LB6 = [L057, DHCB6, L007, SHB6, L014, PUEHB6, L012, QHB6, L031, LDHH, SPTMB6, LDHH];
LB7 = [L057, DVCB7, L007, SVB7, L014, PUEVB7, L011, QVB7, L029, DHB7]
LB8 = [L057, DHCB8, L007, SHB8, L014, PUEHB8, L012, QHB8, L031, DHB8]




LLF = vcat( [L057, DVCB1, L007, SVB1, L014, PUEVB1, L011, QVB1, L029, DHB1],
            [L057, DHCB2, L007, SHB2, L014, PUEHB2, L012, QHB2, L031, DHB2],
            [L057, DVCB3, L007, SVB3, L014, PUEVB3, L011, QVB3, L029, LDH,
             L057, DHCB4, L007, SHB4, L014, PUEHB4, L012, QHB4, L031, DHB4],
            [L057, DVCB5, L007, SVB5, L014, PUEVB5, L011, QVB5, L029, DHB5],
            [L057, DHCB6, L007, SHB6, L014, PUEHB6, L012, QHB6, L031, LDHH]);

# temp for SEB--------------------------------- 

# %%
# Superperiods and Machine

Asuper = vcat(LA1, LA2, LA3, LA4, LA5, LA6, LA7, LA8);
Bsuper = vcat(LB1, LB2, LB3, LB4, LB5, LB6, LB7, LB8);
Csuper = vcat(LC1, LC2, LC3, LC4, LC5, LC6, LC7, LC8);
Dsuper = vcat(LD1, LD2, LD3, LD4, LD5, LD6, LD7, LD8);
Esuper = vcat(LE1, LE2, LE3, LE4, LE5, LE6, LE7, LE8);
Fsuper = vcat(LF1, LF2, LF3, LF4, LF5, LF6, LF7, LF8);




BOOSTER = vcat(Asuper, Bsuper, Csuper, Dsuper, Esuper, Fsuper);

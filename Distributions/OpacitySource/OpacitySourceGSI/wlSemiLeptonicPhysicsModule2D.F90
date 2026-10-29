!!! modified date: June 17, 2020
!!! Authors: Codes wirtten by Gang Guo, based on the formalism by
!!! Andreas Lohs & Gang Guo
!!! neutrino CC opacity is given by 2D integrals
!!! processes to be considered in this routine
!!! ReactionIndex=1: v + n -> e- + p  (opacity in 1/cm)
!!! ReactionIndex=2: vb + p -> e+ + n (opacity in 1/cm)
!!! ReactionIndex=3: vb + p + e- -> n (opacity in 1/cm)
!!! ReactionIndex=4: n -> vb + e- + p [emissivity in v/(s cm^3 MeV)]

!=======================================================================
!  wlSemiLeptonicPhysicsModule2D
!
!  Physics shared by all 2D integration methods of the semi-leptonic CC
!  rates (Guo et al. 2020): constants and form-factor parameters, the
!  angle-integrated amplitude (Calc_AmpCoeffs, Calc_Ampsq), kinematic
!  windows and thermal cutoffs (Range_pn(_D), Threshold_Cut, Ebounds(_D),
!  Q0_Window_Envelope, Collinear_Roots), Fermi-Dirac factors (FD) and the
!  Gauss-Legendre rule (gauleg). No quadrature layout lives here: see
!  IntegrationMethods/ and the front end wlSemiLeptonicOpacityModule2D.
!=======================================================================
MODULE wlSemiLeptonicPhysicsModule2D

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mpi, Vud, mn, mp, &
   massA, massV, gamma_p, gamma_n

  IMPLICIT NONE
  PRIVATE

  !===========================================================================
  !  Semi-leptonic charged-current (CC) neutrino rates via 2D integration
  !
  !  ReactionIndex:
  !    1 : nu capture on neutron
  !    2 : anti-nu capture on proton
  !    3 : inverse neutron decay
  !    4 : neutron decay emissivity
  !
  !  WhichCorrection:
  !    0 : LO only
  !    1 : LO + weak magnetism (WM)
  !    2 : LO + WM + pseudoscalar (PS)
  !    3 : LO + WM + PS + form-factor dependence
  !
  !  Notes:
  !   - For ReactionIndex=4, final-state neutrino blocking is not included.
  !===========================================================================

  ! If you want to reproduce exactly the GSI numbers you need:
  ! REAL(DP),  PARAMETER :: Gw_MeV=1.166d-11,Vud=0.97427d0, F2wm0 = 3.706d0, &
  !      ga =1.2723d0,gv=1.d0, Mpi=139.57d0, Mnp=938.919d0, &
  !      Dnp=1.293d0, massA=1.0d3,massV=840.d0 , GfVud2 = (Gw_MeV*Vud)**2
  REAL(DP), PARAMETER :: F2wm0 = gamma_p - gamma_n - 1.0d0
  REAL(DP), PARAMETER :: Tfac  = 30.0d0           ! thermal cutoff, relative to threshold (Threshold_Cut)
  REAL(DP), PARAMETER :: TfacMin = 30.0d0         ! min. cutoff above rest mass (dilute matter)
  REAL(DP), PARAMETER :: FD_MaxExponent = 700.0d0  ! EXP overflow guard in FD
  REAL(DP), PARAMETER :: gasq  = ga**2
  REAL(DP), PARAMETER :: gvsq  = gv**2
  REAL(DP), PARAMETER :: gva   = gv*ga

  ! Internal tolerances/constants used in bounds logic
  REAL(DP), PARAMETER :: ROOT_REL_TOL      = 1.0d-8   ! root matching in Ebounds (capture)
  REAL(DP), PARAMETER :: Q0_FALLBACK_TFAC  = 50.0d0   ! thermal cap factor for open q0 windows in Ebounds_D

  ! Max Gauss-Legendre points per dimension. Fixed-size locals keep the
  ! integrators free of automatic arrays (which GPU offload compilers handle
  ! poorly); nE inputs are clamped to this value.
  INTEGER, PARAMETER :: nQuadMax = 256

  ! Per-run constants for Calc_Ampsq (hoisted out of the per-node hot path).
  ! MassNuc is the average bare nucleon mass entering the WM/PS/FF hadronic
  ! currents (Guo Eq. 26); the *massSq combinations are the mass parameters
  ! of the leading-order form-factor block.
  REAL(DP), PARAMETER :: MassNuc  = 0.5d0*( mn + mp )
  REAL(DP), PARAMETER :: AmassSq  = massA**2
  REAL(DP), PARAMETER :: VmassSq  = 1.0d0/(1.0d0/massV**2-F2wm0/MassNuc**2/8.0d0)
  REAL(DP), PARAMETER :: AVmassSq = 1.0d0/(0.5d0/massA**2+0.5d0/massV**2-F2wm0/MassNuc**2/16.0d0)
  REAL(DP), PARAMETER :: AFmassSq = 1.0d0/(0.5d0/massA**2+0.5d0/massV**2+1.0d0/MassNuc**2/16.0d0)
  REAL(DP), PARAMETER :: VFmassSq = 1.0d0/(1.0d0/massV**2-(F2wm0-1.0d0)/16.0d0/MassNuc**2)
  REAL(DP), PARAMETER :: FmassSq  = 1.0d0/(1.0d0/massV**2+1.0d0/8.0d0/MassNuc**2)
  REAL(DP), PARAMETER :: MassMy   = 0.25d0*AmassSq+mpi**2

  ! Per-call amplitude coefficients (depend on reaction/thermo state only,
  ! not on the quadrature node); filled once per call by Calc_AmpCoeffs.
  TYPE :: AmpCoeffsType
    REAL(DP) :: dU2, dmf, Qmass
    REAL(DP) :: gvf, gaf, f2sq
    REAL(DP) :: cplA, cplB, cplC, cplD
    REAL(DP) :: cplA2, cplB2, cplC2, cplD2
  END TYPE AmpCoeffsType

  ! Hidden output-unit control for opacity channels.
  ! Output is in cm^-1 (hbarc in MeV cm); output_in_km_inverse = .false.
  ! makes opacity_unit_factor exactly 1. NOTE: the km switch is unused and
  ! would be wrong if enabled (1/cm -> 1/km needs a factor 1e5, cm_to_km is
  ! a length conversion), and only the Energy scheme applies it.
  LOGICAL,  PARAMETER :: output_in_km_inverse = .false.
  REAL(DP), PARAMETER :: cm_to_km             = 1.0d-5
  REAL(DP), PARAMETER :: opacity_unit_factor  = MERGE(cm_to_km, 1.0d0, output_in_km_inverse)

  PUBLIC :: AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq
  PUBLIC :: Range_pn, Range_pn_D, Threshold_Cut
  PUBLIC :: Ebounds, Ebounds_D, Q0_Window_Envelope, Collinear_Roots
  PUBLIC :: Calc_PseudoscalarMoments
  PUBLIC :: FD, gauleg
  PUBLIC :: FDPair
  PUBLIC :: nQuadMax, Tfac, TfacMin, Q0_FALLBACK_TFAC
  PUBLIC :: opacity_unit_factor

CONTAINS

!=======================================================================
!  Calc_AmpCoeffs: per-call amplitude coefficients (node-independent)
!=======================================================================
SUBROUTINE Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, C )

  REAL(DP), INTENT(IN)  :: Mass2, Mass4, U2, U4, anti
  INTEGER , INTENT(IN)  :: IncludeCorrections(3)
  TYPE(AmpCoeffsType), INTENT(OUT) :: C
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  C%dU2   = U2 - U4
  C%dmf   = Mass2 - Mass4
  C%Qmass = 0.5d0*( Mass2**2 - Mass4**2 )

  IF(IncludeCorrections(1) .eq. 1) THEN
     C%gvf  = gv*F2wm0
     C%gaf  = ga*F2wm0
     C%f2sq = F2wm0**2
  ELSE
     C%gvf  = 0.0d0
     C%gaf  = 0.0d0
     C%f2sq = 0.0d0
  END IF

  C%cplA = (gvsq+2.0d0*anti*gva+gasq)+anti*2.0d0*C%gaf*Mass2/MassNuc*(1.0d0-C%dmf/2.0d0/Mass2)
  C%cplB = (gvsq-2.0d0*anti*gva+gasq)-anti*2.0d0*C%gaf*Mass2/MassNuc*(1.0d0-C%dmf/2.0d0/Mass2)
  C%cplC = C%f2sq/MassNuc**2
  C%cplD =-C%f2sq/MassNuc**2

  C%cplA2 = gvsq/VmassSq+gasq/AmassSq+2.0d0*anti*gva/AVmassSq+anti*2.0d0*&
      C%gaf/AFmassSq*Mass2/MassNuc*(1.0d0-C%dmf/2.0d0/Mass2)
  C%cplB2 = gvsq/VmassSq+gasq/AmassSq-2.0d0*anti*gva/AVmassSq-anti*2.0d0*&
      C%gaf/AFmassSq*Mass2/MassNuc*(1.0d0-C%dmf/2.0d0/Mass2)
  C%cplC2 = C%f2sq/MassNuc**2/FmassSq
  C%cplD2 =-C%f2sq/MassNuc**2/FmassSq

END SUBROUTINE Calc_AmpCoeffs

!!! calculate the amplitudes, keep the following unchanged
SUBROUTINE Calc_Ampsq( Enu, E2, E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, Ampsq )

  REAL(DP), INTENT(IN)  :: Enu, E2, E3, Mass2, Mass3, Mass4, U2, U4
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3)
  TYPE(AmpCoeffsType), INTENT(IN) :: Coeffs
  REAL(DP), INTENT(OUT) :: Ampsq
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: P1,P2,P3,P4,Pmax,Pmin,Pmax2,Pmin2,Pmax3,Pmin3
  REAL(DP) :: Ia,Ib,Ic,Id,Ie,Ig,Ih,Ij,Ik,Il,A,B,Ct,D,E,del0,del1,del2,eps0,eps1,alp0,alp1,alp2,bet0,bet1
  REAL(DP) :: cplE,cplG,cplH,cplJ,cplK,cplL,JFF,KFF,LFF
  REAL(DP) :: PaD1,PaD3,PaD5,PaD7,PaD9,PaDm1,PaDm3
  REAL(DP) :: PbD1,PbD3,PbD5,PbD7,PbDm1
  REAL(DP) :: PcD1,PcD3,PcD5,PcD7,PcD9,PcDm1
  REAL(DP) :: E2f,E4f,E1,Qmass,dmf,dU2,E4,P2sq,P3sq,P4sq
  REAL(DP) :: Iam,Iap,Ibm,Ibp,Icm,Icp,Idp,Iem,Iep,Igp
  REAL(DP) :: gvf,f2sq
  REAL(DP) :: tm3sq,tbet0,raw_ampsq
  REAL(DP) :: txi1,Ia2,tgam0,Ib2,txi3,Ic2,Id2,Ie2,Ig2,Ih2,Ij2,Ik2,Il2
  REAL(DP) :: cplE2,cplG2,cplH2,cplJ2,cplK2,cplL2
  REAL(DP) :: ZAP, Moments(7), IlAP, IlPP, IkAP, tIkAP, IkPP, tIkPP
  LOGICAL :: ValidPseudoscalar
  REAL(DP) :: IjAP,tIjAP,IGPP,tIgPP,cplLAP,cplKAP,cplJAP,cplLPP,cplKPP,cplGPP
  REAL(DP) :: ampsq_ff0,ampsq_pseudo
  REAL(DP) :: vb1,vc1,vb3,vc3,vg2,vg1,vg0,vgm1,vgm2,vh0,vh1,vh2
  REAL(DP) :: deltaP,S2,S3,S4,S5,S6,S7,S8,S9

  Ampsq = 0.0d0

  ! Per-call quantities hoisted to Calc_AmpCoeffs (MassNuc, the *massSq
  ! form-factor constants and MassMy are module parameters now)
  dU2   = Coeffs%dU2
  dmf   = Coeffs%dmf
  Qmass = Coeffs%Qmass
  gvf   = Coeffs%gvf
  f2sq  = Coeffs%f2sq

  E1 = Enu
  E4 = E1 + E2 - E3
  E4f = E4 - U4
  E2f = E2 - U2
  ! Off-shell particle (outside the kinematic window): no contribution.
  ! The radicands are tested before the SQRT so that no NaN is ever
  ! produced (FPE traps in debug builds, NaN tests removed by fast-math).
  ! The tested values are the ones passed to SQRT: with FMA contraction a
  ! recomputed radicand can round differently from the tested one.
  IF( E2f < Mass2 .OR. E4f < Mass4 ) RETURN
  P3sq = E3**2 - Mass3**2
  P2sq = E2f**2 - Mass2**2
  P4sq = E4f**2 - Mass4**2
  IF( P3sq < 0.0d0 .OR. P2sq < 0.0d0 .OR. P4sq < 0.0d0 ) RETURN

  P1 = E1
  P3 = SQRT( P3sq )
  P2 = SQRT( P2sq )
  P4 = SQRT( P4sq )

  IF ((P1+P2)>(P3+P4)) THEN
     Pmax = P3+P4
  ELSE
     Pmax = P1+P2
  endif
  IF ((ABS(P1-P2))>(ABS(P3-P4))) THEN
     Pmin = ABS(P1-P2)
  ELSE
     Pmin = ABS(P3-P4)
  endif

  ! No overlap of the angular ranges: no allowed configuration. Pmin = 0
  ! (P1 = P2 and P3 = P4 exactly) is a measure-zero point where the 1/Pmin
  ! terms below are singular; skip it.
  IF( Pmax <= Pmin .OR. Pmin <= 0.0d0 ) RETURN

  ! b^n - a^n = (b-a)*S_n with S_n = b*S_{n-1} + a^(n-1) (all terms
  ! positive, no cancellation for narrow windows)
  deltaP = Pmax - Pmin
  S2 = Pmax + Pmin
  S3 = Pmax*S2 + Pmin**2
  S4 = Pmax*S3 + Pmin**3
  S5 = Pmax*S4 + Pmin**4
  S6 = Pmax*S5 + Pmin**5
  S7 = Pmax*S6 + Pmin**6
  S8 = Pmax*S7 + Pmin**7
  S9 = Pmax*S8 + Pmin**8
  PaD1  = deltaP
  PaD3  = deltaP*S3
  PaD5  = deltaP*S5
  PaD7  = deltaP*S7
  PaD9  = deltaP*S9
  PaDm1 = -deltaP/(Pmax*Pmin)
  PaDm3 = -deltaP*S3/(Pmax*Pmin)**3

  A = E1*E2f+0.5d0*(P1**2+P2**2)
  B = E3*E4f+0.5d0*(P3**2+P4**2)
  Ia = pi**2/15.0d0*(3.0d0*PaD5-10.0d0*(A+B)*PaD3+60.0d0*A*B*PaD1)

  IF ((P1+P4)>(P3+P2)) THEN
    Pmax2 = P3+P2
  ELSE
    Pmax2 = P1+P4
  endif
  IF ((ABS(P1-P4))>(ABS(P3-P2))) THEN
    Pmin2 = ABS(P1-P4)
  ELSE
    Pmin2 = ABS(P3-P2)
  endif

  IF( Pmax2 <= Pmin2 .OR. Pmin2 <= 0.0d0 ) RETURN

  deltaP = Pmax2 - Pmin2
  S2 = Pmax2 + Pmin2
  S3 = Pmax2*S2 + Pmin2**2
  S4 = Pmax2*S3 + Pmin2**3
  S5 = Pmax2*S4 + Pmin2**4
  S6 = Pmax2*S5 + Pmin2**5
  S7 = Pmax2*S6 + Pmin2**6
  PbD1  = deltaP
  PbD3  = deltaP*S3
  PbD5  = deltaP*S5
  PbD7  = deltaP*S7
  PbDm1 = -deltaP/(Pmax2*Pmin2)

  Ct = E1*E4f-0.5d0*(P1**2+P4**2)
  D = E2f*E3-0.5d0*(P2**2+P3**2)
  Ib = pi**2/15.0d0*(3.0d0*PbD5+10.0d0*(Ct+D)*PbD3+60.0d0*Ct*D*PbD1)

  del0 = -(P1**2-P2**2)*(P3**2-P4**2)/4.0d0
  del1 = E1*E3 + (-P1**2+P2**2-P3**2+P4**2)/4.0d0
  del2 = -0.25d0
  eps0 = A   ! identical to A above
  eps1 = -0.5d0
  Ic = pi**2*4.0d0*(del2*eps1**2/7.0d0*PaD7+(2.0d0*del2*eps0*eps1+del1*eps1**2)/5.0d0*PaD5+(del2*eps0**2+2.0d0*del1*eps0&
       *eps1+del0*eps1**2)/3.0d0*PaD3+(del1*eps0**2+2.0d0*del0*eps0*eps1)*PaD1-del0*eps0**2*PaDm1)

  IF ((P1+P3)>(P2+P4)) THEN
    Pmax3 = P2+P4
  ELSE
    Pmax3 = P1+P3
  endif
  IF ((ABS(P1-P3))>(ABS(P2-P4))) THEN
    Pmin3 = ABS(P1-P3)
  ELSE
    Pmin3 = ABS(P2-P4)
  endif

  IF( Pmax3 <= Pmin3 .OR. Pmin3 <= 0.0d0 ) RETURN

  deltaP = Pmax3 - Pmin3
  S2 = Pmax3 + Pmin3
  S3 = Pmax3*S2 + Pmin3**2
  S4 = Pmax3*S3 + Pmin3**3
  S5 = Pmax3*S4 + Pmin3**4
  S6 = Pmax3*S5 + Pmin3**5
  S7 = Pmax3*S6 + Pmin3**6
  S8 = Pmax3*S7 + Pmin3**7
  S9 = Pmax3*S8 + Pmin3**8
  PcD1  = deltaP
  PcD3  = deltaP*S3
  PcD5  = deltaP*S5
  PcD7  = deltaP*S7
  PcD9  = deltaP*S9
  PcDm1 = -deltaP/(Pmax3*Pmin3)

  alp0 = (P1**2-P3**2)*(P2**2-P4**2)/4.0d0
  alp1 = E1*E2f + (P1**2+P2**2-P3**2-P4**2)/4.0d0
  alp2 = 0.25d0
  E = E1*E3 - 0.5d0*(P1**2+P3**2)
  bet0 = E   ! identical to E
  bet1 = 0.5d0
  Id = pi**2*4.0d0*(alp2*bet1**2/7.0d0*PcD7+(2.0d0*alp2*bet0*bet1+alp1*bet1**2)/5.0d0*PcD5+(alp2*bet0**2+2.0d0*alp1*bet0&
       *bet1+alp0*bet1**2)/3.0d0*PcD3+(alp1*bet0**2+2.0d0*alp0*bet0*bet1)*PcD1-alp0*bet0**2*PcDm1)

  Ie = pi**2/15.0d0*(3.0d0*PaD5-20.0d0*A*PaD3+60.0d0*A**2*PaD1)

  Ig = pi**2/15.0d0*(3.0d0*PcD5+20.0d0*E*PcD3+60.0d0*E**2*PcD1)

  Ih = pi**2*4.0d0*(alp2*bet1/5.0d0*PcD5+(alp2*bet0+alp1*bet1)/3.0d0*PcD3+(alp1*bet0+alp0*bet1)*PcD1-alp0*bet0*PcDm1)

  Ij = pi**2/15.0d0*(-10.0d0*PaD3+60.0d0*A*PaD1)

  Ik = pi**2/15.0d0*(10.0d0*PcD3+60.0d0*E*PcD1)

  Il = pi**2/15.0d0*(60.0d0*PaD1)

  cplE = F2sq/MassNuc**2*(-0.5d0*Mass3**2+dU2*(E3-E1)-0.5d0*dU2**2)
  cplG = gvf*Mass2/MassNuc*(2.0d0-dmf/Mass2)+0.5d0*F2sq/MassNuc**2*(Mass2*Mass4-Qmass+0.25d0*Mass3**2-dU2*(E1+E2f)&
       -dU2**2/4.0d0)
  cplH = 0.5d0*F2sq/MassNuc**2*(2.0d0*Qmass+Mass3**2+dU2*(3.0d0*E1-E3+2.0d0*E4f))

  JFF = -Mass3**2*(E1+0.5d0*E2f+0.5d0*E4f)+Qmass*(E3-3.0d0*E1) +0.5d0*dU2*(E4f*(3.0d0*E3-5.0d0*E1)+E2f*(E3+E1)+E3**2-E1**2&
       -2.0d0*Qmass) + dU2**2*(E1-E3-E4f)+0.5d0*dU2**3
  JFF = JFF*dU2
  cplJ = gvf*dmf/MassNuc*0.5d0*(Mass3**2-dU2*(E1+E3))+F2sq/MassNuc**2*0.5d0*JFF

  kFF = -(Mass2+3.0d0*Mass4)*Mass2*0.25d0*Mass3**2 + Qmass**2 + Qmass*0.25d0*Mass3**2-Mass3**4/8.0d0 +dU2*(0.5d0*Qmass&
       *(3.0d0*E1-E2f+E3+3.0d0*E4f)+0.25d0*Mass3**2*(2.0d0*E2f+E3+E1)+Mass2*Mass4*(E3-E1)) +dU2**2*( 0.25d0*(Mass2**2-Mass2&
       *Mass4-3.0d0*Qmass)+E4f*0.5d0*(2.0d0*E1-E2f+E3+E4f)+E2f*(0.5d0*E1-E3)+0.5d0*E1**2) +dU2**3*0.25d0*( -E1+2.0d0*E2f-E3&
       -2.0d0*E4f ) + 0.125d0*dU2**4
  cplK = (gasq-gvsq)*Mass2*Mass4+gvf*Mass2/MassNuc*0.5d0*(-3.0d0*Mass3**2+4.0d0*dU2*(E3-E1)-dU2**2 +dmf/Mass2*(2.0d0*Qmass&
       +Mass3**2+dU2*(E4f+2.0d0*E1-E3)))+0.5d0*F2sq/MassNuc**2*KFF

  LFF = Mass3**2*(Mass2+Mass4)**2-4.0d0*Qmass**2+dU2*(-Mass3**2*(E2f+E4f+E1)-2.0d0*E3*(Mass2**2+Mass2*Mass4) +2.0d0*Qmass&
       *(E2f-3.0d0*E4f-E1)) +2.0d0*dU2**2*( E2f*E3+E4f*(E2f-E4f-E1) + Qmass ) + dU2**3*(-E2f+E4f+E1)
  LFF = LFF*dU2*E1*0.25d0
  cplL = gvf*Mass2/MassNuc*dU2*E1*(Mass3**2-dU2*E3-0.5d0*dmf/Mass2*(Qmass+0.5d0*Mass3**2+dU2*E4f-0.5d0*dU2**2))   +0.5d0&
       *F2sq/MassNuc**2*LFF

  raw_ampsq = Coeffs%cplA*Ia + Coeffs%cplB*Ib + Coeffs%cplC*Ic + Coeffs%cplD*Id + cplE*Ie + cplG*Ig + cplH*Ih + cplJ*Ij &
       + cplK*Ik + cplL*Il

  ampsq_ff0    = 0.0d0
  ampsq_pseudo = 0.0d0
  ! Apply crossing symmetry minus sign for 3-to-1 and 1-to-3 processes
  IF (ReactionIndex == 3 .OR. ReactionIndex == 4) THEN
     Ampsq = max(0.0d0, -raw_ampsq)
  ELSE
     Ampsq = max(0.0d0, raw_ampsq)
  END IF

!!! adding form-factor dependences
  tm3sq = Mass3**2 + 2.0d0*dU2*(E1-E3) + dU2**2
  IF(IncludeCorrections(3) .eq. 1) THEN
    tbet0 = (P1**2-P2**2)*(P4**2-P3**2)/4.0d0
    txi1 = 0.5d0*tm3sq-E1*E3+(P1**2-P2**2+P3**2-P4**2)/4.0d0

    Iap = pi**2/15.0d0*(5.0d0/7.0d0*3.0d0*PaD7-3.0d0/5.0d0*10.0d0*(A+B)*PaD5+60.0d0/3.0d0*A*B*PaD3)
    Iam = pi**2/15.0d0*(5.0d0/3.0d0*3.0d0*PaD3-3.0d0*10.0d0*(A+B)*PaD1-60.0d0*A*B*PaDm1)
    Ia2 = 4.0d0*( 0.5d0*Iap + 2.0d0*txi1*Ia - 2.0d0*tbet0*Iam )

    tgam0 = (P1**2-P4**2)*(P2**2-P3**2)/4.0d0

    Ibp = pi**2/15.0d0*(5.0d0/7.0d0*3.0d0*PbD7+3.0d0/5.0d0*10.0d0*(Ct+D)*PbD5+60.0d0/3.0d0*Ct*D*PbD3)
    Ibm = pi**2/15.0d0*(5.0d0/3.0d0*3.0d0*PbD3+3.0d0*10.0d0*(Ct+D)*PbD1-60.0d0*Ct*D*PbDm1)
    Ib2 = 4.0d0*( 0.5d0*Ibp + 2.0d0*txi1*Ib - 2.0d0*tgam0*Ibm )

    Icp = pi**2*4.0d0*(del2*eps1**2/9.0d0*PaD9+(2.0d0*del2*eps0*eps1+del1*eps1**2)/7.0d0*PaD7+(del2*eps0**2+2.0d0*del1*eps0&
       *eps1+del0*eps1**2)/5.0d0*PaD5+(del1*eps0**2+2.0d0*del0*eps0*eps1)/3.0d0*PaD3+del0*eps0**2*PaD1)
    Icm =  pi**2*4.0d0*(del2*eps1**2/5.0d0*PaD5 + (2.0d0*del2*eps0*eps1+del1*eps1**2)/3.0d0*PaD3+(del2*eps0**2+2.0d0*del1&
       *eps0*eps1+del0*eps1**2)*PaD1-(del1*eps0**2+2.0d0*del0*eps0*eps1)*PaDm1-1.0d0/3.0d0*del0*eps0**2*PaDm3)

    Ic2 = 4.0d0*( 0.5d0*Icp + 2.0d0*txi1*Ic - 2.0d0*tbet0*Icm )
    ! <(p1.p3) q*^2>_phi = <p1.p3><q*^2> - 2 Var(p1vec.p3vec), Var = p1perp^2 p3perp^2/2:
    ! Ic2 -= 16 pi^2 Int dp (A-p^2/2)^2 p1perp^2 p3perp^2, with (u=p^2)
    ! p1perp^2 = -(u - b1 + c1/u)/4, b1 = 2(P1^2+P2^2), c1 = (P1^2-P2^2)^2 (same for 3,4)
    vb1 = 2.0d0*(P1**2+P2**2);  vc1 = (P1**2-P2**2)**2
    vb3 = 2.0d0*(P3**2+P4**2);  vc3 = (P3**2-P4**2)**2
    vg2 = 1.0d0; vg1 = -(vb1+vb3); vg0 = vb1*vb3+vc1+vc3
    vgm1 = -(vb1*vc3+vb3*vc1); vgm2 = vc1*vc3
    vh0 = A**2; vh1 = -A; vh2 = 0.25d0
    Ic2 = Ic2 - pi**2*( (vh2*vg2)                    *PaD9/9.0d0 &
                      + (vh2*vg1+vh1*vg2)            *PaD7/7.0d0 &
                      + (vh2*vg0+vh1*vg1+vh0*vg2)    *PaD5/5.0d0 &
                      + (vh2*vgm1+vh1*vg0+vh0*vg1)   *PaD3/3.0d0 &
                      + (vh2*vgm2+vh1*vgm1+vh0*vg0)  *PaD1       &
                      - (vh1*vgm2+vh0*vgm1)          *PaDm1      &
                      - (vh0*vgm2)                   *PaDm3/3.0d0 )


    txi3 = 0.5d0*tm3sq-E1*E3+(P1**2+P3**2)/2.0d0
    Idp = pi**2*4.0d0*(alp2*bet1**2/9.0d0*PcD9+(2.0d0*alp2*bet0*bet1+alp1*bet1**2)/7.0d0*PcD7+(alp2*bet0**2+2.0d0*alp1*bet0&
       *bet1+alp0*bet1**2)/5.0d0*PcD5+(alp1*bet0**2+2.0d0*alp0*bet0*bet1)/3.0d0*PcD3+alp0*bet0**2*PcD1)
    Id2 = 4.0d0*( -Idp + 2.0d0*txi3*Id )

    Iep = pi**2/15.0d0*(5.0d0/7.0d0*3.0d0*PaD7-3.0d0/5.0d0*20.0d0*A*PaD5+60.0d0/3.0d0*A**2*PaD3)
    Iem = pi**2/15.0d0*(5.0d0/3.0d0*3.0d0*PaD3-3.0d0*20.0d0*A*PaD1-60.0d0*A**2*PaDm1)
    Ie2 = 4.0d0*( 0.5d0*Iep + 2.0d0*txi1*Ie - 2.0d0*tbet0*Iem )

    Igp = pi**2/15.0d0*(5.0d0/7.0d0*3.0d0*PcD7+3.0d0/5.0d0*20.0d0*E*PcD5+60.0d0/3.0d0*E**2*PcD3)
    Ig2 = 4.0d0*( -Igp + 2.0d0*txi3*Ig )

    IH2 = 4.0d0*( tm3sq*Ih - 2.0d0*Id )
    Ij2 = 4.0d0*( tm3sq*Ij - 2.0d0*Ih )
    Ik2 = 4.0d0*( tm3sq*Ik - 2.0d0*IG )
    Il2 = 4.0d0*( tm3sq*Il - 2.0d0*Ik )

    cplE2 = F2sq/FmassSq/MassNuc**2*(-0.5d0*Mass3**2+dU2*(E3-E1)-0.5d0*dU2**2)
    cplG2 = gvf/VFmassSq*Mass2/MassNuc*(2.0d0-dmf/Mass2)+0.5d0*F2sq/FmassSq/MassNuc**2*(Mass2*Mass4-Qmass+0.25d0*Mass3**2&
       -dU2*(E1+E2f)-dU2**2/4.0d0)
    cplH2 = 0.5d0*F2sq/FmassSq/MassNuc**2*(2.0d0*Qmass+Mass3**2+dU2*(3.0d0*E1-E3+2.0d0*E4f))
    cplJ2 = gvf/VFmassSq*dmf/MassNuc*0.5d0*(Mass3**2-dU2*(E1+E3))+F2sq/FmassSq/MassNuc**2*0.5d0*JFF
    cplK2 = (gasq/AmassSq-gvsq/VmassSq)*Mass2*Mass4+gvf/VFmassSq*Mass2/MassNuc*0.5d0*(-3.0d0*Mass3**2+4.0d0*dU2*(E3-E1)&
       -dU2**2+dmf/Mass2*(2.0d0*Qmass+Mass3**2+dU2*(E4f+2.0d0*E1-E3)))+0.5d0*F2sq/MassNuc**2*KFF/FmassSq
    cplL2 = gvf/VFmassSq*Mass2/MassNuc*dU2*E1*(Mass3**2-dU2*E3-0.5d0*dmf/Mass2*(Qmass+0.5d0*Mass3**2+dU2*E4f-0.5d0*dU2**2))&
       +0.5d0*F2sq/MassNuc**2*LFF/FmassSq

    Ampsq_ff0 = Coeffs%cplA2*Ia2 + Coeffs%cplB2*Ib2 + Coeffs%cplC2*Ic2 + Coeffs%cplD2*Id2 + cplE2*Ie2 + cplG2*Ig2 + cplH2&
       *Ih2 + cplJ2*Ij2 + cplK2*Ik2 + cplL2*Il2
     IF(Ampsq+ampsq_ff0>=0.0d0) Ampsq=Ampsq+ampsq_ff0
  END IF

!!! adding pseudoscalar terms w/o form factor dependence
  IF(IncludeCorrections(2) .eq. 1) THEN

    ZAP = mpi**2-tm3sq+2.0d0*E1*E3-P1**2-P3**2

    CALL Calc_PseudoscalarMoments( Pmin3, Pmax3, ZAP, Moments, ValidPseudoscalar )
    IF( .NOT. ValidPseudoscalar )THEN
      ! A negative amplitude reports an invalid denominator to the caller.
      Ampsq = -1.0_DP
      RETURN
    END IF

    IF( IncludeCorrections(3) == 0 )THEN
      IlAP = 4.0_DP*pi**2*Moments(1)
      IlPP = 4.0_DP*pi**2*Moments(5)
      tIkAP = 2.0_DP*pi**2*Moments(2)
      tIkPP = 2.0_DP*pi**2*Moments(6)
      tIjAP = 4.0_DP*pi**2*alp0*Moments(4)
      tIgPP = pi**2*Moments(7)
    ELSE
      IlAP = 16.0_DP*pi**2/AmassSq*(MassMy*Moments(1)-PcD1)
      IlPP = 16.0_DP*pi**2/AmassSq*(MassMy*Moments(5)-Moments(1))
      tIkAP = 8.0_DP*pi**2/AmassSq*(MassMy*Moments(2)-PcD3/3.0_DP)
      tIkPP = 8.0_DP*pi**2/AmassSq*(MassMy*Moments(6)-Moments(2))
      tIjAP = 16.0_DP*pi**2/AmassSq*alp0*(MassMy*Moments(4)-1.0_DP/Pmin3+1.0_DP/Pmax3)
      tIgPP = 4.0_DP*pi**2/AmassSq*(MassMy*Moments(7)-Moments(3))
    END IF

      IkAP = E*IlAP+tIkAP
      IkPP = E*IlPP+tIkPP
      IjAP = alp1*IlAP+0.5d0*tIkAP+tIjAP
      IgPP = E**2*IlPP+2.0d0*E*tIkPP+tIgPP

      cplLAP = 2.0d0*MassNuc*gasq*dU2*E1*(2.0d0*Mass2*Mass3**2-dmf*(Qmass+0.5d0*Mass3**2)-dU2*(2.0d0*E3*Mass2+E4f*dmf)+0.5d0&
       *dU2**2*dmf)
      cplKAP = 2.0d0*MassNuc*gasq*(Mass2*(dU2**2-Mass3**2)+dmf*dU2*(E2f+E1))
      cplJAP = 2.0d0*MassNuc*gasq*dmf*(Mass3**2-dU2*(E1+E3))

      cplLPP = 2.0d0*MassNuc**2*gasq*dU2*E1*(Mass3**4-Mass3**2*dmf**2+dU2*(-Mass3**2*(3.0d0*E3-2.0d0*E1)+dmf**2*E3)+dU2**2&
       *(Mass3**2+2.0d0*E3*(E3-E1))-dU2**3*E3)
      cplKPP = 2.0d0*MassNuc**2*gasq*(-0.5d0*Mass3**4+0.5d0*dmf**2*Mass3**2+dU2*Mass3**2*(E3-3.0d0*E1)+dU2**2*(2.0d0*E1*E3&
       -0.5d0*dmf**2)+dU2**3*(E4f-E2f)-0.5d0*dU2**4)
      cplGPP = 2.0d0*MassNuc**2*gasq*(Mass3**2-dU2**2)

      ampsq_pseudo = cplLAP*IlAP+cplKAP*IkAP+cplJAP*IjAP + cplLPP*IlPP +cplKPP*IkPP+cplGPP*IgPP
      IF(Ampsq+ampsq_pseudo>=0.0d0) Ampsq=Ampsq+ampsq_pseudo
    END IF

  ! Crossing (Guo Eq. 39): for the 3-body channels EVERY term of |M|^2 (LO, WM, FF, PS)
  ! changes sign. Apply the sign to the full sum once, then clip once.
  IF (ReactionIndex == 3 .OR. ReactionIndex == 4) THEN
    Ampsq = MAX(0.0d0, -(raw_ampsq + ampsq_ff0 + ampsq_pseudo))
  END IF

END SUBROUTINE Calc_Ampsq

PURE SUBROUTINE Calc_PseudoscalarMoments( A, B, Z, Moments, Valid )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
  REAL(DP), INTENT(IN) :: A, B, Z
  REAL(DP), INTENT(OUT) :: Moments(7)
  LOGICAL, INTENT(OUT) :: Valid
  REAL(DP) :: Width, A2, B2, DA, DB, H2, HPos2, R, R2, Delta
  REAL(DP) :: C1, C3, C5, P3, P5, Power, Ratio, V1, V3, V5
  REAL(DP) :: S1, S3, S5, T1, T3, T5, Root, Argument, Y
  INTEGER :: K

  ! Moments are integrals of Q**n/(Q**2+Z)**r over [A,B]:
  ! (n,r) = (0,1),(2,1),(4,1),(-2,1),(0,2),(2,2),(4,2).
  Moments = 0.0_DP
  Valid = .FALSE.
  IF( A <= 0.0_DP .OR. B < A ) RETURN
  A2 = A*A
  B2 = B*B
  DA = A2+Z
  DB = B2+Z
  IF( DA <= 0.0_DP .OR. DB <= 0.0_DP ) RETURN
  Valid = .TRUE.
  Width = B-A
  IF( Width == 0.0_DP ) RETURN
  R = A/B
  R2 = R*R
  Delta = Width/B
  C1 = Delta
  C3 = Delta*(1.0_DP+R+R2)
  H2 = Delta/A
  HPos2 = Width*(A2+A*B+B2)/3.0_DP

  IF( Z == 0.0_DP )THEN
    Moments = [ H2, Width, HPos2, C3/(3.0_DP*A*A2), &
                C3/(3.0_DP*A*A2), H2, Width ]
    RETURN
  END IF

  IF( ABS(Z) <= 0.125_DP*A2 )THEN
    ! Geometric series in -Z/A**2; the omitted term is below roundoff.
    ! Cj = 1-R**j is accumulated without subtracting nearly equal powers.
    S1 = 0.0_DP
    S3 = 0.0_DP
    T1 = 0.0_DP
    T3 = 0.0_DP
    P3 = R*R2
    Power = 1.0_DP
    Ratio = -Z/A2
    DO K = 0, 20
      V1 = Power*C1/REAL(2*K+1,DP)
      V3 = Power*C3/REAL(2*K+3,DP)
      S1 = S1+V1
      S3 = S3+V3
      T1 = T1+REAL(K+1,DP)*V1
      T3 = T3+REAL(K+1,DP)*V3
      C1 = C3
      C3 = C3+P3*Delta*(1.0_DP+R)
      P3 = P3*R2
      Power = Power*Ratio
      IF( ABS(Power)*REAL(K+2,DP) <= EPSILON(1.0_DP)/64.0_DP ) EXIT
    END DO
    Moments(1) = S1/A
    Moments(4) = S3/(A*A2)
    Moments(5) = T3/(A*A2)
    Moments(6) = T1/A
    Moments(2) = Width-Z*Moments(1)
    Moments(3) = HPos2-Z*Moments(2)
    Moments(7) = Moments(2)-Z*Moments(6)
  ELSE IF( Z > 0.0_DP .AND. B2 <= 0.125_DP*Z )THEN
    ! The complementary series avoids cancellation for large positive Z.
    S1 = 0.0_DP
    S3 = 0.0_DP
    S5 = 0.0_DP
    T1 = 0.0_DP
    T3 = 0.0_DP
    T5 = 0.0_DP
    P5 = R*R2*R2
    C5 = C3+R*R2*Delta*(1.0_DP+R)
    Power = 1.0_DP
    Ratio = -B2/Z
    DO K = 0, 20
      V1 = Power*C1/REAL(2*K+1,DP)
      V3 = Power*C3/REAL(2*K+3,DP)
      V5 = Power*C5/REAL(2*K+5,DP)
      S1 = S1+V1
      S3 = S3+V3
      S5 = S5+V5
      T1 = T1+REAL(K+1,DP)*V1
      T3 = T3+REAL(K+1,DP)*V3
      T5 = T5+REAL(K+1,DP)*V5
      C1 = C3
      C3 = C5
      C5 = C5+P5*Delta*(1.0_DP+R)
      P5 = P5*R2
      Power = Power*Ratio
      IF( ABS(Power)*REAL(K+2,DP) <= EPSILON(1.0_DP)/64.0_DP ) EXIT
    END DO
    Moments(1) = B/Z*S1
    Moments(2) = B*B2/Z*S3
    Moments(3) = B*B2*B2/Z*S5
    Moments(4) = (H2-Moments(1))/Z
    Moments(5) = B/Z/Z*T1
    Moments(6) = B*B2/Z/Z*T3
    Moments(7) = B*B2*B2/Z/Z*T5
  ELSE
    IF( Z > 0.0_DP )THEN
      Root = SQRT(Z)
      Moments(1) = ATAN2(Root*Width,Z+A*B)/Root
    ELSE
      Root = SQRT(-Z)
      ! log(1+x) with the rounding of 1+x compensated for small x.
      Argument = 2.0_DP*Root*Width/((A-Root)*(B+Root))
      Y = 1.0_DP+Argument
      IF( Y == 1.0_DP )THEN
        Moments(1) = Argument/(2.0_DP*Root)
      ELSE
        Moments(1) = LOG(Y)*(Argument/(Y-1.0_DP))/(2.0_DP*Root)
      END IF
    END IF
    Moments(2) = Width-Z*Moments(1)
    Moments(3) = HPos2-Z*Moments(2)
    Moments(4) = (H2-Moments(1))/Z
    Moments(5) = (Moments(1)+Width*(Z-A*B)/DA/DB)/(2.0_DP*Z)
    Moments(6) = Moments(1)-Z*Moments(5)
    Moments(7) = Moments(2)-Z*Moments(6)
  END IF
END SUBROUTINE Calc_PseudoscalarMoments


SUBROUTINE Ebounds( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, xq0_min, xq0_max )

  ! Allowed q0 = Enu - E3 range for captures at fixed E2 (Guo Eq. C1),
  ! evaluated in closed form by Q0_Window_Envelope (s3 = -1).
  ! On an empty window a zero-width interval is returned (the E3
  ! integration in Integral_2D then skips this E2 node).

  REAL(DP), INTENT(IN)  :: E2   ! on-shell: E2 = SQRT(P2**2+Mass2**2)+U2 (set by caller)
  REAL(DP), INTENT(IN)  :: Enu, T, P2, Mass2, Mass3, Mass4, U2, U4
  REAL(DP), INTENT(OUT) :: xq0_min, xq0_max
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  LOGICAL  :: found

  ! Captures: q0 <= Enu - Mass3 always, so that is the (never-active) cap
  CALL Q0_Window_Envelope( -1.0d0, Enu, P2, Mass3, Mass4, U4-E2, Enu-Mass3, &
                           xq0_min, xq0_max, found )
  IF( .NOT. found ) THEN
    xq0_min = Enu - Mass3   ! zero-width window: E3 collapses to Mass3
    xq0_max = xq0_min
  END IF

END SUBROUTINE Ebounds

SUBROUTINE Range_pn( Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, p2min, p2max )

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2
  REAL(DP), INTENT(OUT) :: p2min, p2max
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: tmp,mpt,Epr,Esq,Equ,p2a,p2b,p20,F0,Finf,Fmax,Fmin

  ! Thermal cutoff mu-U+Tfac*T, but at least TfacMin*T above the rest
  ! mass: mu-U+Tfac*T alone drops below Mass2 in dilute matter (SQRT of a
  ! negative number). Unchanged wherever mu-U >= Mass2-(Tfac-TfacMin)*T.
  p2min = 0.0d0
  p2max = SQRT( MAX(Mu2 - U2 + Tfac*T, Mass2 + TfacMin*T)**2 - Mass2**2 )

  mpt = Mass4 + Mass3
  F0 = SQRT(Enu**2 + mpt**2) + U4 - Mass2 - U2
  Finf = -Enu + U4 - U2
  IF(mpt < Mass2) THEN
    tmp = (1.d0+mpt/Mass2)/(1.d0-mpt**2/Mass2**2)*Enu
    Fmin = SQRT( (tmp-Enu)**2+mpt**2 )+U4-SQRT(tmp**2+Mass2**2)-U2
    Fmax = max( F0, Finf )
  ELSE
    Fmin = Finf
    Fmax = F0
  END IF

  IF( Fmax <= Enu) THEN
    CALL Threshold_Cut( p2min, Mass2, Tfac*T, p2max )
    RETURN
  END IF
  IF( Fmin >= Enu) THEN
    p2max = -1000.0d0
    RETURN
  END IF

  Epr = Enu - U4 + U2
  Esq = Enu**2 + mpt**2 - Mass2**2 - Epr**2
  Equ = Esq**2 + 4.d0*Mass2**2*( Enu**2-Epr**2 )
  IF(Equ>0.0d0 .and. Enu**2 .ne. Epr**2) THEN
    p2a = ( Enu*Esq - ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    p2b = ( Enu*Esq + ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    IF(Epr<0.d0 .and. p2b>0.d0) THEN
      p2min=p2b
    ELSE IF(Epr>=0.d0 .and. p2a>0.d0) THEN
      p2min=p2a
    END IF

  ELSE IF(Enu**2 .eq. Epr**2) THEN
    p20 = -(4.0d0*Epr**2*Mass2**2-Esq**2)/(4.0d0*Esq*Enu)
    IF(p20>0.d0) p2min=p20
  END IF

  CALL Threshold_Cut( p2min, Mass2, Tfac*T, p2max )

END SUBROUTINE Range_pn

!=======================================================================
!  Threshold_Cut: make the thermal cutoff relative to the kinematic
!  threshold. pmax is raised (never lowered) so that the integration
!  covers at least dE above the threshold energy SQRT(pmin**2+Mass**2).
!  The mu-relative cutoff alone truncates rates whose threshold lies far
!  above mu (e.g. muonic inverse decay / decay), where the whole
!  integrand sits above mu - U + Tfac*T.
!=======================================================================
SUBROUTINE Threshold_Cut( pmin, Mass, dE, pmax )

  REAL(DP), INTENT(IN)    :: pmin, Mass, dE
  REAL(DP), INTENT(INOUT) :: pmax
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: Ecut

  Ecut = SQRT( pmin**2 + Mass**2 ) + dE
  pmax = MAX( pmax, SQRT( Ecut**2 - Mass**2 ) )

END SUBROUTINE Threshold_Cut

SUBROUTINE Ebounds_D( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max, ReactionIndex, xq0_min, xq0_max )

  REAL(DP), INTENT(IN)  :: E2   ! on-shell: E2 = SQRT(P2**2+Mass2**2)+U2 (set by caller)
  REAL(DP), INTENT(IN)  :: Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max
  INTEGER,  INTENT(IN)  :: ReactionIndex
  REAL(DP), INTENT(OUT) :: xq0_min, xq0_max
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  ! Allowed q0 = Enu + E3 range for inverse decay / decay at fixed E2
  ! (Guo Eq. C16), evaluated in closed form by Q0_Window_Envelope (s3 = +1).
  ! If the window never closes as |q| -> infinity (the paper's "no upper
  ! bound" case), q0 is capped thermally via q0_cap. On an empty window a
  ! zero-width interval is returned.

  REAL(DP) :: q0_cap
  LOGICAL  :: found

  ! Thermal cap on q0 for the open-window case:
  !  reaction 3: q0 = Enu + E3, absorbed lepton cut off by its occupation.
  !              Measured from MAX(Mu3, Mass3): with Mu3 < Mass3 (e.g. no
  !              muons) a cap from Mu3 or 0 can fall below E3 = Mass3 and
  !              truncate the window (breaks detailed balance with reaction 4).
  !  reaction 4: q0 = E4 - E2, emitting nucleon capped consistently with p4max
  IF( ReactionIndex == 3 ) THEN
    q0_cap = Enu + MAX(Mu3, Mass3) + Q0_FALLBACK_TFAC*T
  ELSE
    q0_cap = SQRT(p4max**2 + Mass4**2) + U4 - E2
  END IF

  CALL Q0_Window_Envelope( +1.0d0, Enu, P2, Mass3, Mass4, U4-E2, q0_cap, &
                           xq0_min, xq0_max, found )
  IF( .NOT. found ) THEN
    xq0_min = Enu + Mass3   ! zero-width window: E3 collapses to Mass3
    xq0_max = xq0_min
  END IF

END SUBROUTINE Ebounds_D

!=======================================================================
!  Q0_Window_Envelope: allowed q0 range over momentum transfer |q|,
!  computed in closed form.
!
!  With y = E4 = q0 - U4mE2 and W = Enu - U4mE2 (so E3 = s3*(y - W);
!  s3 = -1: capture, q0 = Enu - E3;  s3 = +1: decay, q0 = Enu + E3),
!  the lepton and nucleon constraints at fixed |q| read
!    p3 in [ |Enu - q|, Enu + q ]   and   p4 in [ |P2 - q|, P2 + q ].
!  At fixed y both p3 and p4 are fixed, and some |q| satisfies both iff
!  the lengths (Enu, P2, p3, p4) close a quadrilateral, i.e.
!    p3 + p4 >= |Enu - P2|   and   |p3 - p4| <= Enu + P2.
!  The allowed y set is therefore closed and its extrema can only sit at
!    E4 = Mass4,  E3 = Mass3,  or a collinear boundary |p3 +/- p4| = P,
!  P in { |Enu - P2|, Enu + P2 }. Each collinear boundary is a two-body
!  configuration with E4 - s3*E3 = W and invariant s = W**2 - P**2,
!  solved exactly. Any real collinear root with E3 > 0 is a genuine
!  allowed configuration (|p3 +/- p4| = D or S implies the quadrilateral
!  condition, and y**2 >= Mass4**2 follows from the quadratic), so the
!  roots only need the sign test on E3 (squaring admits E3 < 0); allowed
!  non-boundary roots lie inside the hull and do not move the extrema.
!  The full momentum test is needed only for the at-rest candidates: on a
!  collinear root next to p4 = 0 or p3 = 0 it would amplify the round-off
!  in y through SQRT(y**2 - Mass**2) and reject the true bound.
!  For decays the window stays open as y -> infinity
!  (p3 - p4 -> -W) iff |W| < Enu + P2; q0 is then capped by the
!  caller-supplied thermal cap q0_cap.
!=======================================================================
SUBROUTINE Q0_Window_Envelope( s3, Enu, P2, Mass3, Mass4, U4mE2, q0_cap, &
                               xq0_min, xq0_max, found )

  REAL(DP), INTENT(IN)  :: s3, Enu, P2, Mass3, Mass4, U4mE2, q0_cap
  REAL(DP), INTENT(OUT) :: xq0_min, xq0_max
  LOGICAL , INTENT(OUT) :: found
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: W, D, S, tol, y(2)
  INTEGER  :: iP, ir, nRoot

  W   = Enu - U4mE2
  D   = ABS( Enu - P2 )
  S   = Enu + P2
  tol = 1.0d-12*( ABS(W) + S + Mass3 + Mass4 + 1.0d0 )

  found   = .FALSE.
  xq0_min =  HUGE(1.0d0)
  xq0_max = -HUGE(1.0d0)

  ! Candidates for the extrema of y: a particle at rest (full test) ...
  IF( Allowed( Mass4 ) )          CALL Accept( Mass4 )
  IF( Allowed( W + s3*Mass3 ) )   CALL Accept( W + s3*Mass3 )

  ! ... or the collinear configurations (both roots, cancellation-free)
  DO iP = 1, 2
    CALL Collinear_Roots( W, MERGE( D, S, iP == 1 ), Mass3, Mass4, nRoot, y )
    DO ir = 1, nRoot
      CALL AcceptRoot( y(ir) )
    END DO
  END DO

  IF( .NOT. found ) RETURN

  xq0_min = xq0_min + U4mE2
  xq0_max = xq0_max + U4mE2

  ! Window open as q0 -> infinity (decay only) -> apply thermal cap
  IF( s3 > 0.0d0 .AND. ABS(W) < S ) xq0_max = q0_cap
  xq0_max = MIN( xq0_max, q0_cap )
  IF( xq0_max < xq0_min ) xq0_max = xq0_min   ! cap below window: zero width

CONTAINS

  SUBROUTINE Accept( y )
    REAL(DP), INTENT(IN) :: y
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    found   = .TRUE.
    xq0_min = MIN( xq0_min, y )
    xq0_max = MAX( xq0_max, y )
  END SUBROUTINE Accept

  SUBROUTINE AcceptRoot( y )
    ! collinear root: genuine iff E4 > 0 and E3 = s3*(y - W) > 0
    REAL(DP), INTENT(IN) :: y
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    IF( y > 0.0d0 .AND. s3*( y - W ) > 0.0d0 ) CALL Accept( y )
  END SUBROUTINE AcceptRoot

  LOGICAL FUNCTION Allowed( y )
    REAL(DP), INTENT(IN) :: y
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP) :: E3, p3, p4
    E3 = s3*( y - W )
    Allowed = ( E3 >= Mass3 - tol .AND. y >= Mass4 - tol )
    IF( .NOT. Allowed ) RETURN
    p3 = SQRT( MAX( E3**2 - Mass3**2, 0.0d0 ) )
    p4 = SQRT( MAX( y**2  - Mass4**2, 0.0d0 ) )
    Allowed = ( p3 + p4 >= D - tol .AND. ABS(p3 - p4) <= S + tol )
  END FUNCTION Allowed

END SUBROUTINE Q0_Window_Envelope

!=======================================================================
!  Collinear_Roots: two-body configurations with E4 - s3*E3 = W and
!  |p3 +/- p4| = P (invariant s = W**2 - P**2), returned as y = E4 - U4.
!  Both roots of the quadratic, cancellation-free. sInv and the Kallen
!  function lam are factored so that lam = 0 exactly at the two-body
!  threshold (no spurious lam < 0 from round-off). Squaring admits
!  E3 < 0: callers apply their own sign test.
!=======================================================================
SUBROUTINE Collinear_Roots( W, P, Mass3, Mass4, nRoot, y )

  REAL(DP), INTENT(IN)  :: W, P, Mass3, Mass4
  INTEGER , INTENT(OUT) :: nRoot
  REAL(DP), INTENT(OUT) :: y(2)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: sInv, a, lam, t, K

  nRoot = 0
  y     = 0.0d0
  sInv = ( W - P )*( W + P )
  a    = sInv + Mass4**2 - Mass3**2
  lam  = ( sInv - (Mass4 + Mass3)**2 )*( sInv - (Mass4 - Mass3)**2 )
  IF( lam < 0.0d0 ) RETURN
  t = W*a + SIGN( P*SQRT(lam), W*a )
  K = 0.25d0*( a**2 + 4.0d0*(P*Mass4)**2 )
  IF( sInv /= 0.0d0 ) THEN
    nRoot = nRoot + 1
    y(nRoot) = 0.5d0*t/sInv
  END IF
  IF( t /= 0.0d0 ) THEN
    nRoot = nRoot + 1
    y(nRoot) = 2.0d0*K/t
  END IF

END SUBROUTINE Collinear_Roots

SUBROUTINE Range_pn_D( Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4, ReactionIndex, p2min, p2max, p4max )

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4
  INTEGER , INTENT(IN)  :: ReactionIndex
  REAL(DP), INTENT(OUT) :: p2min, p2max, p4max
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: tmp,mpt,Epr,Esq,Equ,p2a,p2b,p20
  REAL(DP) :: F0, Fmax, E4cut

  p2min = 0.0d0
  ! Thermal cutoffs as in Range_pn
  p4max = SQRT( MAX(Mu4 - U4 + Tfac*T, Mass4 + TfacMin*T)**2 - Mass4**2 )
  SELECT CASE(ReactionIndex)
  CASE (3)
    p2max = SQRT( MAX(Mu2 - U2 + Tfac*T, Mass2 + TfacMin*T)**2 - Mass2**2 )
  CASE (4)
    ! MAX(...,0): with FMA contraction x**2 - Mass2**2 can round to -eps
    ! when x = Mass2, and the SQRT would return NaN.
    p2min = SQRT( MAX( (max(Mass2,Mu2-Tfac*T-U2))**2 - Mass2**2, 0.0d0 ) )
    p2max = SQRT( MAX( (SQRT(p4max**2+Mass4**2)+U4-U2)**2-Mass2**2, 0.0d0 ) )
  CASE DEFAULT
    ! Invalid input: no I/O or STOP in device-callable code; empty range
    p2min = 0.0d0
    p2max = -1000.0d0
    RETURN
  END SELECT

  mpt = Mass4 - Mass3
  IF( mpt > Mass2) THEN
    tmp = Mass2/(mpt-Mass2)*Enu
  ELSE
    ! F increases with p2 here: test the p2 -> infinity limit, since the
    ! thermal p2max is extended past the threshold below
    tmp = MAX( p2max, 1.0d6 )
  END IF
  Fmax = SQRT( (tmp+Enu)**2+mpt**2 )+U4-SQRT(tmp**2+Mass2**2)-U2

  IF(Enu>=Fmax) THEN
    p2max = -1000.0d0
    RETURN
  END IF

  F0 = SQRT(Enu**2 + mpt**2)+U4-Mass2-U2
  IF( .NOT. ( F0 >= Enu .and. U4 >= U2 ) ) THEN
    Epr = Enu - U4 + U2
    Esq = Enu**2 + mpt**2 - Mass2**2 - Epr**2
    Equ = Esq**2 + 4.d0*( (Mass2*Enu)**2 - (Epr*Mass2)**2 )
    IF(Equ>0.0d0.and. Enu**2 .ne. Epr**2) THEN
      p2a = -( Enu*Esq - ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
      p2b = -( Enu*Esq + ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
      IF(Epr<0.d0 .and. p2b>0.d0) THEN
        p2min=p2b
      ELSE IF(Epr>=0.0d0 .and. p2a>0.d0) THEN
        p2min=p2a
      END IF

    ELSE IF(Enu**2 .eq. Epr**2) THEN
      p20 = (4.0d0*Epr**2*Mass2**2-Esq**2)/(4.0d0*Esq*Enu)
      IF(p20>0.d0) p2min=p20
    END IF
  END IF

  ! Thermal cutoff relative to the kinematic threshold (see Threshold_Cut).
  ! Reaction 4: E4 >= E2(p2min) + Enu + Mass3; cut on E4, and on E2 <= E4.
  IF( ReactionIndex == 3 ) THEN
    CALL Threshold_Cut( p2min, Mass2, Tfac*T, p2max )
  ELSE
    E4cut = SQRT( p2min**2 + Mass2**2 ) + U2 + Enu + Mass3 + Tfac*T
    IF( E4cut - U4 > SQRT( p4max**2 + Mass4**2 ) ) THEN
      p4max = SQRT( (E4cut - U4)**2 - Mass4**2 )
      p2max = MAX( p2max, SQRT( (E4cut - U2)**2 - Mass2**2 ) )
    END IF
  END IF

END SUBROUTINE Range_pn_D

!=======================================================================
!
!     NUMERICAL RECIPES: Gauss-Legendre integration
!
!=======================================================================
SUBROUTINE gauleg(x1, x2, x, w, n)

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: x1, x2
  REAL(DP), INTENT(OUT), dimension(n) :: x, w
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
  !-----------------------------------------------------------------------
  !
  !	input:
  !	x1   lower integration boundary
  !	x2   upper integration boundary
  !	n    discretization
  !
  !	output:
  !	x   position of integrand evaluation
  !	w   weight
  !
  !-----------------------------------------------------------------------
  INTEGER :: i,j,m
  REAL(DP), PARAMETER :: eps=3.d-14
  REAL(DP) :: p1,p2,p3,pp,xl,xm,z,z1

  m = (n+1)/2
  xm = 0.5d0*(x2 + x1)
  xl = 0.5d0*(x2 - x1)

  DO i=1,m
    z = COS(pi*(i - 0.25d0)/(n + 0.5d0))
    DO
      p1 = 1.d0
      p2 = 0.d0
      DO j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.d0*j - 1.d0)*z*p2 - (j - 1.d0)*p3)/j
      END DO
      pp = n * (z * p1 - p2) / (z*z - 1.d0)
      z1 = z
      z = z1 - p1/pp
      IF (ABS(z - z1) <= eps) EXIT
    END DO
    x(i)        = xm - xl*z
    x(n + 1 - i)= xm + xl*z
    w(i)        = 2.d0 * xl / ((1.d0 - z*z) * pp * pp)
    w(n + 1 - i)= w(i)
  END DO

END SUBROUTINE gauleg

PURE ELEMENTAL FUNCTION FD(E, mu, T)
  REAL(DP), INTENT(IN) :: E, mu, T
  REAL(DP) :: FD, x, z
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  IF( T <= 0.0d0 )THEN
    FD = MERGE( 1.0d0, 0.0d0, E < mu )
    RETURN
  END IF

  x = ( E - mu ) / T
  z = EXP( -ABS(x) )
  FD = MERGE( z, 1.0d0, x >= 0.0d0 ) * ( 1.0d0 / ( 1.0d0 + z ) )
END FUNCTION

PURE ELEMENTAL SUBROUTINE FDPair(E, mu, T, f, vacancy)
  REAL(DP), INTENT(IN)  :: E, mu, T
  REAL(DP), INTENT(OUT) :: f, vacancy
  REAL(DP) :: x, z, inverse
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  IF( T <= 0.0d0 )THEN
    f = MERGE( 1.0d0, 0.0d0, E < mu )
    vacancy = 1.0d0-f
    RETURN
  END IF
  x = ( E-mu ) / T
  z = EXP( -ABS(x) )
  inverse = 1.0d0 / ( 1.0d0+z )
  f       = MERGE( z, 1.0d0, x >= 0.0d0 ) * inverse
  vacancy = MERGE( 1.0d0, z, x >= 0.0d0 ) * inverse
END SUBROUTINE FDPair


END MODULE wlSemiLeptonicPhysicsModule2D

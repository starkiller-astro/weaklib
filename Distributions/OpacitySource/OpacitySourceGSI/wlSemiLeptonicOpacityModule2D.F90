!!! modified date: June 17, 2020
!!! Authors: Codes wirtten by Gang Guo, based on the formalism by
!!! Andreas Lohs & Gang Guo
!!! neutrino CC opacity is given by 2D integrals
!!! processes to be considered in this routine
!!! ReactionIndex=1: v + n -> e- + p  (oapcity in 1/km)
!!! ReactionIndex=2: vb + p -> e+ + n (opacity in 1/km)
!!! ReactionIndex=3: vb + p + e- -> n (opacity in 1/km)
!!! ReactionIndex=4: n -> vb + e- + p [emissivity in v/(s cm^3 MeV)]

!!! nE defined in module: the No. of grids used in gaussian
!!! quadrature for each dimention
!!! seems nE needs be >~30 to reach an accuracy of 5%

MODULE wlSemiLeptonicOpacityModule2D

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mpi, Vud, mn, mp, &
   massA, massV, gamma_p, gamma_n, hbarc, hbar
  ! USE wlEosConstantsModule, ONLY: &
  !  pi, hbarc, hbar, mn, mp

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
  REAL(DP), PARAMETER :: Tfac  = 100.0d0
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

  ! Hidden output-unit control for opacity channels
  ! Default preserves current behavior: km^-1.
  ! Set .false. to return cm^-1 without changing public API.
  LOGICAL,  PARAMETER :: output_in_km_inverse = .false.
  REAL(DP), PARAMETER :: cm_to_km             = 1.0d-5
  REAL(DP), PARAMETER :: opacity_unit_factor  = MERGE(cm_to_km, 1.0d0, output_in_km_inverse)

  !========================================================================
! Routine to calculate the neutrino CC opacity via 2D integrals
! Inputs:
! (1) WhichCorrection=0,1,2,3 for LO, LO+weak magnetism(WM),
! LO+WM+pseduoscalar(PS), LO+WM+PS+form factor dependencies
! (2) ReactionIndex=ReactionIndex=1,2,3,4 for nu capture on n, nub capture on p, inverse
! neutron decay, neutron decay emissivity
! (3) Enu(NP): neutrino energy array in MeV with dimension NP
! (4) xT: Temperature in MeV
! (5) xMul, xml: relativistic chemical potentials & masses of
! e^-(mu^-) in MeV
! (6) xMun,xMup,xmn,xmp,xUn,xUp: relativistic chemical
! potentials, masses, & potentials of neutron and proton in MeV
! relativsitic disperstion relation for nucleons,
! E2=SQRT(p^2+Mass2^2)+U2
! Outputs:
! (7) OpaA(NP): opacities in km^-1; for neutron decay, emissivities
! in (s cm^3 MeV)^-1
!     **** blocking factor of final state neutrinos in neutron decay
!     is not considered ****
!=============================================================
  PUBLIC :: Opacity_CC_2D
  PUBLIC :: Opacity_CC_2D_PrecomputedNodes
  PUBLIC :: gauleg
  ! Building blocks exported for the 2-phase GPU integration layer
  ! (wlSemiLeptonicOpacityIntegrationModule2D)
  PUBLIC :: AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq
  PUBLIC :: Range_pn, Range_pn_D, Ebounds, Ebounds_D, FD
  PUBLIC :: opacity_unit_factor

CONTAINS

SUBROUTINE Opacity_CC_2D(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp, nE)

  ! Back-compatible entry point: computes the [0,1] Gauss-Legendre nodes
  ! locally. In performance-critical (GPU / in-flight) code precompute the
  ! nodes once with gauleg and call Opacity_CC_2D_PrecomputedNodes instead.
  INTEGER , INTENT(IN)  :: ReactionIndex, WhichCorrection, nE
  REAL(DP), INTENT(IN)  :: Enu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp
  REAL(DP), INTENT(OUT) :: OpaA
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: xa(nQuadMax), wxa(nQuadMax)
  INTEGER  :: nEuse

  nEuse = MIN(nE, nQuadMax)
  CALL gauleg(0.d0, 1.d0, xa, wxa, nEuse)
  CALL Opacity_CC_2D_PrecomputedNodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
       xml, xmn, xmp, xUn, xUp, nEuse, xa, wxa)

END SUBROUTINE Opacity_CC_2D

SUBROUTINE Opacity_CC_2D_PrecomputedNodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp, nE, xa, wxa)

  ! Same physics as Opacity_CC_2D, but takes precomputed Gauss-Legendre
  ! nodes/weights on [0,1] (from CALL gauleg(0.d0, 1.d0, xa, wxa, nE)).
  INTEGER , INTENT(IN)  :: ReactionIndex, WhichCorrection, nE
  REAL(DP), INTENT(IN)  :: Enu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)
  REAL(DP), INTENT(OUT) :: OpaA
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: anti
  INTEGER  :: IncludeCorrections(3)

  SELECT CASE(WhichCorrection)
  CASE (0)
     IncludeCorrections = (/ 0, 0, 0 /)
  CASE (1)
     IncludeCorrections = (/ 1, 0, 0 /)
  CASE (2)
     IncludeCorrections = (/ 1, 1, 0 /)
  CASE (3)
     IncludeCorrections = (/ 1, 1, 1 /)
  CASE DEFAULT
     ! Invalid input: no I/O or STOP in device-callable code; return zero
     OpaA = 0.0d0
     RETURN
  END SELECT

  SELECT CASE(ReactionIndex)
  CASE (1)
    anti = 1.0d0
    CALL Integral_2D(Enu, xT, xmn, xml, xmp, xUn, xUp, xMun, xMul, xMup, anti, 1, IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (2)
    anti = -1.0d0
    CALL Integral_2D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, -xMul, xMun, anti, 2, IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (3)
    anti = -1.0d0
    CALL Integral_2D_D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 3, IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (4)
    anti = -1.0d0
    CALL Integral_2D_D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 4, IncludeCorrections, nE, xa, wxa, OpaA)
  CASE DEFAULT
    OpaA = 0.0d0
  END SELECT

END SUBROUTINE Opacity_CC_2D_PrecomputedNodes

! For a more efficient parallelization this whole thing can be probably taken out, and we make 4 different integrations, one per reaction
SUBROUTINE Integral_2D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, IncludeCorrections, nE, xa, wxa, res)   ! captures

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)   ! Gauss-Legendre nodes/weights on [0,1]
  REAL(DP), INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: P2,xq0_min,xq0_max,E3_min,E3_max,E2_min,E2_max
  REAL(DP) :: wE2, wE3
  INTEGER  :: i, j
  REAL(DP) :: xamp, E2, E3, E4
  REAL(DP) :: xf2,xf3,xf4
  REAL(DP) :: p2min, p2max
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  CALL Range_pn(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, p2min, p2max)
  IF( p2min > p2max ) RETURN
  E2_min = SQRT(p2min**2 + Mass2**2) + U2
  E2_max = SQRT(p2max**2 + Mass2**2) + U2

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )

  ! Loop over incoming particle (i.e. particle 2)
  DO i=1, nE
    E2 = xa(i)*(E2_max - E2_min) + E2_min
    wE2 = wxa(i)*(E2_max - E2_min)
    P2 = SQRT( (E2 - U2)**2 - Mass2**2 )
    CALL Ebounds( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, xq0_min, xq0_max )
    E3_min = Enu - xq0_max
    E3_max = Enu - xq0_min
    IF( E3_max - E3_min <= 0.0d0 ) CYCLE   ! empty kinematic window at this E2
    xf2 = FD(E2, Mu2, T)
  ! Loop over outgoing particle (i.e. particle 3)
    DO j=1,nE
      E3 = xa(j)*(E3_max - E3_min) + E3_min
      wE3 = wxa(j)*(E3_max - E3_min)
      E4 = Enu + E2 - E3
      CALL Calc_ampsq( Enu, E2, E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, xamp)
      xf3 = FD(E3, Mu3, T)
      xf4 = FD(E4, Mu4, T)
      res = res + xf2*(1.0d0-xf3)*(1.0d0-xf4)*xamp*wE2*wE3/Enu**2
    END DO
  END DO

  ! This is in 1/km
  res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
  ! This is in 1/km or 1/cm, depending on the value of output_in_km_inverse
  res = res*opacity_unit_factor

END SUBROUTINE Integral_2D

! For a more efficient parallelization this whole thing can be probably taken out, and we make 4 different integrations, one per reaction
SUBROUTINE Integral_2D_D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, IncludeCorrections, nE, xa, wxa, res)   ! decay/inverse decay

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)   ! Gauss-Legendre nodes/weights on [0,1]
  REAL(DP), INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: P2,xq0_min,xq0_max,E3_min,E3_max,E2_min,E2_max
  REAL(DP) :: wE2, wE3
  INTEGER :: i,j
  REAL(DP) :: xamp
  REAL(DP) :: E2, E3, E4
  REAL(DP) :: xf2,xf3,xf4
  REAL(DP) :: p2min, p2max, p4max
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  CALL Range_pn_D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4, ReactionIndex, p2min, p2max, p4max)
  IF( p2min > p2max ) RETURN
  E2_min = SQRT(p2min**2 + Mass2**2) + U2
  E2_max = SQRT(p2max**2 + Mass2**2) + U2

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )

  ! Loop over incoming particle (i.e. particle 2)
  DO i=1, nE
    E2 = xa(i)*(E2_max - E2_min) + E2_min
    wE2 = wxa(i)*(E2_max - E2_min)
    P2 = SQRT( (E2-U2)**2 - Mass2**2 )
    CALL Ebounds_D( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max, ReactionIndex, xq0_min, xq0_max )
    E3_min = xq0_min - Enu
    E3_max = xq0_max - Enu
    IF( E3_max - E3_min <= 0.0d0 ) CYCLE   ! empty kinematic window at this E2
    xf2 = FD(E2, Mu2, T)
  ! Loop over outgoing particle (i.e. particle 3)
    DO j=1,nE
      E3 = xa(j)*(E3_max - E3_min) + E3_min
      wE3 = wxa(j)*(E3_max - E3_min)
      E4 = Enu + E3 + E2
      CALL Calc_ampsq(Enu, E2, -E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, xamp)
      xf3 = FD(E3, Mu3, T)
      xf4 = FD(E4, Mu4, T)
      IF(ReactionIndex.eq.3) THEN ! inverse decay
        res = res + xf2*xf3*(1.0d0-xf4)*xamp*wE2*wE3/Enu**2
      ELSE IF(ReactionIndex.eq.4) THEN ! n decay, *blocking of neutrinos is not added*
        res = res + (1.0d0-xf2)*(1.0d0-xf3)*xf4*xamp*wE2*wE3/Enu**2
      END IF
    END DO
  END DO

  IF(ReactionIndex.eq.3) THEN
  ! This is in 1/km
    res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
  ! This is in 1/km or 1/cm, depending on the value of output_in_km_inverse
    res = res*opacity_unit_factor
  ELSE IF(ReactionIndex.eq.4) THEN
    res = res*(Gw_MeV*Vud)**2/32.0d0/(pi**7)*Enu**2/hbarc**3/hbar
  END IF

END SUBROUTINE Integral_2D_D

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
  REAL(DP) :: E2f,E4f,E1,Qmass,dmf,dU2,E4
  REAL(DP) :: Iam,Iap,Ibm,Ibp,Icm,Icp,Idp,Iem,Iep,Igp
  REAL(DP) :: gvf,f2sq
  REAL(DP) :: tm3sq,tbet0,raw_ampsq
  REAL(DP) :: txi1,Ia2,tgam0,Ib2,txi3,Ic2,Id2,Ie2,Ig2,Ih2,Ij2,Ik2,Il2
  REAL(DP) :: cplE2,cplG2,cplH2,cplJ2,cplK2,cplL2
  REAL(DP) :: YAP,ZAP,Inte1,Inte2,IlAP,IlPP,IkAP,tIkAP,IkPP,tIkPP
  REAL(DP) :: IjAP,tIjAP,IGPP,tIgPP,TanArg,LogArg,cplLAP,cplKAP,cplJAP,cplLPP,cplKPP,cplGPP
  REAL(DP) :: ampsq_ff0,ampsq_pseudo       

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
  P1 = E1
  P3 = SQRT( E3**2 - Mass3**2 )
  P2 = SQRT( (E2 - U2)**2 - Mass2**2 )
  P4 = SQRT( (E4 - U4)**2 - Mass4**2 )

  IF( P1 .ne. P1) RETURN
  IF( P2 .ne. P2) RETURN
  IF( P3 .ne. P3) RETURN
  IF( P4 .ne. P4) RETURN

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

  PaD1  = Pmax - Pmin
  PaD3  = Pmax**3 - Pmin**3
  PaD5  = Pmax**5 - Pmin**5
  PaD7  = Pmax**7 - Pmin**7
  PaD9  = Pmax**9 - Pmin**9
  PaDm1 = 1.0d0/Pmax - 1.0d0/Pmin
  PaDm3 = 1.0d0/Pmax**3 - 1.0d0/Pmin**3

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

  PbD1  = Pmax2 - Pmin2
  PbD3  = Pmax2**3 - Pmin2**3
  PbD5  = Pmax2**5 - Pmin2**5
  PbD7  = Pmax2**7 - Pmin2**7
  PbDm1 = 1.0d0/Pmax2 - 1.0d0/Pmin2

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

  PcD1  = Pmax3 - Pmin3
  PcD3  = Pmax3**3 - Pmin3**3
  PcD5  = Pmax3**5 - Pmin3**5
  PcD7  = Pmax3**7 - Pmin3**7
  PcD9  = Pmax3**9 - Pmin3**9
  PcDm1 = 1.0d0/Pmax3 - 1.0d0/Pmin3

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

    YAP = 0.25d0*AmassSq+tm3sq-2.0d0*E1*E3+P1**2+P3**2
    ZAP = mpi**2-tm3sq+2.0d0*E1*E3-P1**2-P3**2

    ! The branches below cover ZAP>0 and ZAP<0 only; initialize so that the
    ! measure-zero case ZAP==0 gives a vanishing pseudoscalar contribution
    ! instead of reading undefined values.
    Inte1 = 0.0d0;  Inte2 = 0.0d0
    ILAp  = 0.0d0;  IlPP  = 0.0d0
    tIkAP = 0.0d0;  tIkPP = 0.0d0
    tIjAP = 0.0d0;  tIgPP = 0.0d0

    IF(IncludeCorrections(3) .eq. 0) THEN
      IF( ZAP >0.0d0 ) THEN
        Inte1 = ( atan(Pmax3/sqrt(ZAP))-atan(Pmin3/sqrt(ZAP)) )/sqrt(ZAP) 
        Inte2 = 0.5d0*( Pmax3/(Pmax3**2+ZAP)-Pmin3/(Pmin3**2+ZAp) + Inte1 )/ZAP 
        ILAp = 4.0d0*pi**2*Inte1
        Ilpp = 4.0d0*pi**2*inte2
        tIkAp = 2.0d0*pi**2*(PcD1-ZAP*inte1)
        tIkPP = 2.0d0*pi**2*(-ZAP*Inte2+Inte1)
        tIjAp = 4.0d0*pi**2*alp0*( -1.0d0/ZAP*PcDm1  -1.0d0/ZAP*inte1 )
        tIgPP = pi**2*( ZAP**2*Inte2 - 2.0d0*ZAP*Inte1 +PcD1 )
      ELSE IF(ZAP<0.0d0) THEN
        ZAP = -ZAP
        TanArg = sqrt(ZAP)*PcD1/(Pmax3*Pmin3-ZAP)
        LogArg = abs((1.0d0+TanArg)/(1.0d0-TanArg))
        Inte1 = 0.5d0*log(LogArg)/sqrt(ZAP)
        Inte2 = -0.5d0*( Pmax3/(Pmax3**2-ZAP)-Pmin3/(Pmin3**2-ZAp) + Inte1 )/ZAp 
        ILAp = 4.0d0*pi**2*Inte1
        Ilpp = 4.0d0*pi**2*inte2
        tIkAp = 2.0d0*pi**2*(PcD1+ZAP*inte1)
        tIkPP = 2.0d0*pi**2*(ZAP*Inte2+Inte1)
        tIjAp = 4.0d0*pi**2*alp0*( 1.0d0/ZAP*PcDm1  +1.0d0/ZAP*inte1 )
        tIgPP = pi**2*( ZAP**2*Inte2+2.0d0*ZAP*Inte1  +PcD1)
      END IF

    ELSE IF(IncludeCorrections(3) .eq. 1) THEN
      IF( ZAP>0.0d0) THEN
        Inte1 = ( atan(Pmax3/sqrt(ZAP))-atan(Pmin3/sqrt(ZAP)) )/sqrt(ZAP)
        Inte2 = 0.5d0*( Pmax3/(Pmax3**2+ZAP)-Pmin3/(Pmin3**2+ZAp) + Inte1 )/ZAP 
        ILAp = 16.0d0*pi**2/AmassSq*(MassMy*Inte1-PcD1)
        Ilpp = 16.0d0*pi**2/AmassSq*( MassMy*inte2 - inte1 )
        tIkAp = 8.0d0*pi**2/AmassSq*( MassMy*PcD1-PcD3/3.0d0-MassMy*ZAP*inte1 ) 
        tIkPP = 8.0d0*pi**2/AmassSq*( -PcD1-ZAP*MassMy*Inte2+(MassMy+ZAP)*Inte1 ) 
        tIjAp = 16.0d0*pi**2/AmassSq*alp0*( -(MassMy-ZAP)/ZAP*PcDm1-MassMy/ZAP*inte1 )
        tIgPP = 4.0d0*pi**2/AmassSq*( ZAP**2*MassMy*Inte2 - (2.0d0*ZAP*MassMy+ZAp**2)*Inte1+(ZAP+MassMy)*PcD1-PcD3/3.0d0 ) 
      ELSE IF (ZAP<0.0d0) THEN
        ZAP = -ZAP
        TanArg = sqrt(ZAP)*PcD1/(Pmax3*Pmin3-ZAP)
        LogArg = abs((1.0d0+TanArg)/(1.0d0-TanArg))
        Inte1 = 0.5d0*log(LogArg)/sqrt(ZAP)
        Inte2 = -0.5d0*( Pmax3/(Pmax3**2-ZAP)-Pmin3/(Pmin3**2-ZAp) + Inte1 )/ZAp 
        ILAp = 16.0d0*pi**2/AmassSq*(MassMy*Inte1-PcD1)
        Ilpp = 16.0d0*pi**2/AmassSq*( MassMy*inte2 - inte1 )
        tIkAp = 8.0d0*pi**2/AmassSq*( MassMy*PcD1-PcD3/3.0d0+MassMy*ZAP*inte1 ) 
        tIkPP = 8.0d0*pi**2/AmassSq*( -PcD1+ZAP*MassMy*Inte2+(MassMy-ZAP)*Inte1 ) 
        tIjAp = 16.0d0*pi**2/AmassSq*alp0*( (MassMy+ZAP)/ZAP*PcDm1+MassMy/ZAP*inte1 )
        tIgPP = 4.0d0*pi**2/AmassSq*( ZAP**2*MassMy*Inte2 - (-2.0d0*ZAP*MassMy+ZAp**2)*Inte1 +(-ZAP+MassMy)*PcD1&
       -PcD3/3.0d0 ) 
      END IF
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

END SUBROUTINE Calc_Ampsq

SUBROUTINE Ebounds( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, xq0_min, xq0_max )

  ! Allowed q0 = Enu - E3 range for captures at fixed E2 (Guo Eq. C1),
  ! evaluated by the direct window-envelope engine (s3 = -1).
  ! On an empty window a zero-width interval is returned (the E3
  ! integration in Integral_2D then skips this E2 node).

  REAL(DP), INTENT(INOUT) :: E2
  REAL(DP), INTENT(IN)  :: Enu, T, P2, Mass2, Mass3, Mass4, U2, U4
  REAL(DP), INTENT(OUT) :: xq0_min, xq0_max
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: qRoots(4), tmp1, tmp2, xA, xB, xC, disc
  INTEGER  :: nRoots
  LOGICAL  :: found

  E2 = SQRT( P2**2 + Mass2**2 ) + U2

  ! Window-closure roots, Eq. C11 (guarded; used as breakpoint candidates)
  nRoots = 0
  tmp1 = U4 - E2 - Enu
  tmp2 = Enu**2 + Mass3**2 - Mass4**2 - P2**2
  xC   = 4.d0*tmp1**2*(P2**2 + Mass4**2) - (tmp2 - tmp1**2)**2

  xA = 4.d0*( tmp1**2 - (Enu-P2)**2 )
  xB = 4.d0*( (tmp2 - tmp1**2)*Enu - (tmp2 + tmp1**2)*P2 )
  disc = xB**2 - 4.d0*xA*xC
  IF( disc >= 0.0d0 .AND. ABS(xA) > 0.0d0 ) THEN
    qRoots(nRoots+1) = ABS( ( xB + SQRT(disc))/(2.d0*xA) )
    qRoots(nRoots+2) = ABS( (-xB + SQRT(disc))/(2.d0*xA) )
    nRoots = nRoots + 2
  END IF

  xA = 4.d0*( tmp1**2 - (Enu+P2)**2 )
  xB = 4.d0*( -(tmp2 - tmp1**2)*Enu - (tmp2 + tmp1**2)*P2 )
  disc = xB**2 - 4.d0*xA*xC
  IF( disc >= 0.0d0 .AND. ABS(xA) > 0.0d0 ) THEN
    qRoots(nRoots+1) = ABS( ( xB + SQRT(disc))/(2.d0*xA) )
    qRoots(nRoots+2) = ABS( (-xB + SQRT(disc))/(2.d0*xA) )
    nRoots = nRoots + 2
  END IF

  ! Captures: q0 <= Enu - Mass3 always, so that is the (never-active) cap
  CALL Q0_Window_Envelope( -1.0d0, Enu, P2, Mass3, Mass4, U4-E2, qRoots, nRoots, &
                           Enu-Mass3, xq0_min, xq0_max, found )
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

  p2min = 0.0d0
  p2max = SQRT( (Tfac*T + Mu2 - U2)**2 - Mass2**2 )

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

  IF( Fmax <= Enu) RETURN
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

END SUBROUTINE Range_pn

SUBROUTINE Ebounds_D( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max, ReactionIndex, xq0_min, xq0_max )

  REAL(DP), INTENT(INOUT) :: E2
  REAL(DP), INTENT(IN)  :: Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max
  INTEGER,  INTENT(IN)  :: ReactionIndex
  REAL(DP), INTENT(OUT) :: xq0_min, xq0_max
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  ! Allowed q0 = Enu + E3 range for inverse decay / decay at fixed E2
  ! (Guo Eq. C16), evaluated by the direct window-envelope engine (s3 = +1).
  ! If the window never closes as |q| -> infinity (the paper's "no upper
  ! bound" case), q0 is capped thermally via q0_cap. On an empty window a
  ! zero-width interval is returned.

  REAL(DP) :: qRoots(4), tmp1, tmp2, xA, xB, xC, disc
  REAL(DP) :: q0_cap
  INTEGER  :: nRoots
  LOGICAL  :: found

  E2 = SQRT( P2**2 + Mass2**2 ) + U2

  ! Thermal cap on q0 for the open-window case:
  !  reaction 3: q0 = Enu + E3, absorbed lepton cut off by its occupation
  !  reaction 4: q0 = E4 - E2, emitting nucleon capped consistently with p4max
  IF( ReactionIndex == 3 ) THEN
    q0_cap = Enu + MAX(Mu3, 0.0d0) + Q0_FALLBACK_TFAC*T
  ELSE
    q0_cap = SQRT(p4max**2 + Mass4**2) + U4 - E2
  END IF

  ! Window-closure roots, Eq. C25 (guarded; used as breakpoint candidates)
  nRoots = 0
  tmp1 = U4 - E2 - Enu
  tmp2 = Enu**2 + Mass3**2 - Mass4**2 - P2**2
  xC   = 4.d0*tmp1**2*(P2**2 + Mass4**2) - (tmp2 - tmp1**2)**2

  xA = 4.d0*( tmp1**2 - (Enu+P2)**2 )
  xB = 4.d0*( (tmp2 - tmp1**2)*Enu + (tmp2 + tmp1**2)*P2 )
  disc = xB**2 - 4.d0*xA*xC
  IF( disc >= 0.0d0 .AND. ABS(xA) > 0.0d0 ) THEN
    qRoots(nRoots+1) = ABS( ( xB + SQRT(disc))/(2.d0*xA) )
    qRoots(nRoots+2) = ABS( (-xB + SQRT(disc))/(2.d0*xA) )
    nRoots = nRoots + 2
  END IF

  xA = 4.d0*( tmp1**2 - (Enu-P2)**2 )
  xB = 4.d0*( -(tmp2 - tmp1**2)*Enu + (tmp2 + tmp1**2)*P2 )
  disc = xB**2 - 4.d0*xA*xC
  IF( disc >= 0.0d0 .AND. ABS(xA) > 0.0d0 ) THEN
    qRoots(nRoots+1) = ABS( ( xB + SQRT(disc))/(2.d0*xA) )
    qRoots(nRoots+2) = ABS( (-xB + SQRT(disc))/(2.d0*xA) )
    nRoots = nRoots + 2
  END IF

  CALL Q0_Window_Envelope( +1.0d0, Enu, P2, Mass3, Mass4, U4-E2, qRoots, nRoots, &
                           q0_cap, xq0_min, xq0_max, found )
  IF( .NOT. found ) THEN
    xq0_min = Enu + Mass3   ! zero-width window: E3 collapses to Mass3
    xq0_max = xq0_min
  END IF

END SUBROUTINE Ebounds_D

!=======================================================================
!  Q0_Window_Envelope: allowed q0 range over momentum transfer |q|
!
!  At fixed |q| the kinematically allowed q0 interval is
!    [ WinLo(|q|), WinHi(|q|) ]              (Guo Eqs. C1 / C16)
!  built from the lepton-side branches
!    qL_{a,b}(|q|) = Enu + s3*SQRT( (Enu +/- |q|)^2 + Mass3^2 )
!    (s3 = -1: capture, q0 = Enu - E3;  s3 = +1: decay, q0 = Enu + E3)
!  and the nucleon-side branches
!    qH_{a,b}(|q|) = SQRT( (P2 +/- |q|)^2 + Mass4^2 ) + U4 - E2 .
!  The q0 bounds are the extrema of the envelope over all |q| with a
!  nonempty window. Extrema can only occur at branch stationary points
!  (|q| = Enu, |q| = P2), at branch crossings (kinks of the min/max
!  envelopes), or at window opening/closure points. Closure points and
!  crossings are refined by bisection, with the analytic closure roots
!  (Eqs. C11/C25) supplied by the caller as cheap breakpoint candidates.
!  If the window never closes (possible for decay kinematics), the
!  caller-supplied thermal cap q0_cap bounds q0 from above.
!=======================================================================
SUBROUTINE Q0_Window_Envelope( s3, Enu, P2, Mass3, Mass4, U4mE2, qRootsIn, nRootsIn, &
                               q0_cap, xq0_min, xq0_max, found )

  REAL(DP), INTENT(IN)  :: s3, Enu, P2, Mass3, Mass4, U4mE2, q0_cap
  REAL(DP), INTENT(IN)  :: qRootsIn(4)
  INTEGER , INTENT(IN)  :: nRootsIn
  REAL(DP), INTENT(OUT) :: xq0_min, xq0_max
  LOGICAL , INTENT(OUT) :: found
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  INTEGER , PARAMETER :: nCandMax = 8, nScanMax = 2*nCandMax + 2, nBisect = 60
  REAL(DP) :: qCand(nCandMax), qScan(nScanMax)
  LOGICAL  :: allowedScan(nScanMax), sgnScan(nScanMax), sgnA
  REAL(DP) :: q_big, qa, qb, qm, q_eval, tmpq
  INTEGER  :: i, j, k, nCand, nScan

  q_big = P2 + Enu + MAX(0.0d0, -U4mE2) + ABS(q0_cap) + 10.0d0

  ! Candidate |q| breakpoints: branch stationary points and closure roots
  nCand = 0
  CALL AddCand( 0.0d0 )
  CALL AddCand( Enu )
  CALL AddCand( P2 )
  DO i = 1, nRootsIn
    CALL AddCand( qRootsIn(i) )
  END DO

  ! Sort candidates ascending (insertion sort, nCand <= 7)
  DO i = 2, nCand
    tmpq = qCand(i)
    j = i - 1
    DO WHILE( j >= 1 )
      IF( qCand(j) <= tmpq ) EXIT
      qCand(j+1) = qCand(j)
      j = j - 1
    END DO
    qCand(j+1) = tmpq
  END DO

  ! Scan points: candidates, midpoints between neighbours, and the far probe
  nScan = 0
  DO i = 1, nCand
    nScan = nScan + 1;  qScan(nScan) = qCand(i)
    IF( i < nCand ) THEN
      nScan = nScan + 1;  qScan(nScan) = 0.5d0*(qCand(i) + qCand(i+1))
    END IF
  END DO
  nScan = nScan + 1;  qScan(nScan) = 0.5d0*(qCand(nCand) + q_big)
  nScan = nScan + 1;  qScan(nScan) = q_big

  ! Envelope extrema over allowed scan points
  found   = .FALSE.
  xq0_min =  HUGE(1.0d0)
  xq0_max = -HUGE(1.0d0)
  DO i = 1, nScan
    allowedScan(i) = ( SignTest(1, qScan(i)) )
    IF( allowedScan(i) ) THEN
      found = .TRUE.
      CALL UpdateEnvelope( qScan(i) )
    END IF
  END DO

  ! Refine sign-change points by bisection and include them:
  !   k = 1: window opening/closure (gap = WinHi - WinLo)
  !   k = 2: lower-branch crossing  (qL_min - qH_min), kink of WinLo
  !   k = 3: upper-branch crossing  (qL_max - qH_max), kink of WinHi
  DO k = 1, 3
    DO i = 1, nScan
      sgnScan(i) = SignTest(k, qScan(i))
    END DO
    DO i = 1, nScan-1
      IF( sgnScan(i) .NEQV. sgnScan(i+1) ) THEN
        qa   = qScan(i)
        qb   = qScan(i+1)
        sgnA = sgnScan(i)
        DO j = 1, nBisect
          qm = 0.5d0*(qa + qb)
          IF( SignTest(k, qm) .EQV. sgnA ) THEN
            qa = qm
          ELSE
            qb = qm
          END IF
        END DO
        IF( k == 1 ) THEN
          ! take the endpoint on the allowed side of the boundary
          IF( allowedScan(i) ) THEN
            q_eval = qa
          ELSE
            q_eval = qb
          END IF
          found = .TRUE.
          CALL UpdateEnvelope( q_eval )
        ELSE
          ! envelope kink: include only if inside the allowed set
          qm = 0.5d0*(qa + qb)
          IF( SignTest(1, qm) ) THEN
            found = .TRUE.
            CALL UpdateEnvelope( qm )
          END IF
        END IF
      END IF
    END DO
  END DO

  IF( .NOT. found ) RETURN

  ! Window still open at the far probe -> unbounded above, apply thermal cap
  IF( allowedScan(nScan) ) xq0_max = q0_cap
  xq0_max = MIN( xq0_max, q0_cap )
  IF( xq0_max < xq0_min ) xq0_max = xq0_min   ! cap below window: zero width

CONTAINS

  SUBROUTINE AddCand( q )
    REAL(DP), INTENT(IN) :: q
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    IF( q /= q ) RETURN                        ! NaN guard
    IF( q < 0.0d0 .OR. q > q_big ) RETURN
    IF( nCand >= nCandMax ) RETURN
    nCand = nCand + 1
    qCand(nCand) = q
  END SUBROUTINE AddCand

  PURE SUBROUTINE Branches( q, qLmin, qLmax, qHmin, qHmax )
    REAL(DP), INTENT(IN)  :: q
    REAL(DP), INTENT(OUT) :: qLmin, qLmax, qHmin, qHmax
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP) :: qLa, qLb
    qLa   = Enu + s3*SQRT( (Enu + q)**2 + Mass3**2 )
    qLb   = Enu + s3*SQRT( (Enu - q)**2 + Mass3**2 )
    qLmin = MIN( qLa, qLb )
    qLmax = MAX( qLa, qLb )
    qHmin = SQRT( (P2 - q)**2 + Mass4**2 ) + U4mE2
    qHmax = SQRT( (P2 + q)**2 + Mass4**2 ) + U4mE2
  END SUBROUTINE Branches

  FUNCTION SignTest( kind, q ) RESULT( res )
    INTEGER , INTENT(IN) :: kind
    REAL(DP), INTENT(IN) :: q
    LOGICAL :: res
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP) :: qLmin, qLmax, qHmin, qHmax
    CALL Branches( q, qLmin, qLmax, qHmin, qHmax )
    SELECT CASE( kind )
    CASE(1);  res = ( MIN(qLmax,qHmax) - MAX(qLmin,qHmin) >= 0.0d0 )
    CASE(2);  res = ( qLmin - qHmin >= 0.0d0 )
    CASE DEFAULT;  res = ( qLmax - qHmax >= 0.0d0 )
    END SELECT
  END FUNCTION SignTest

  SUBROUTINE UpdateEnvelope( q )
    REAL(DP), INTENT(IN) :: q
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP) :: qLmin, qLmax, qHmin, qHmax
    CALL Branches( q, qLmin, qLmax, qHmin, qHmax )
    xq0_min = MIN( xq0_min, MAX(qLmin, qHmin) )
    xq0_max = MAX( xq0_max, MIN(qLmax, qHmax) )
  END SUBROUTINE UpdateEnvelope

END SUBROUTINE Q0_Window_Envelope

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
  REAL(DP) :: F0, Fmax

  p2min = 0.0d0
  p4max = SQRT( (Tfac*T + Mu4 - U4)**2 - Mass4**2 )
  SELECT CASE(ReactionIndex)
  CASE (3)
    p2max = SQRT( (Tfac*T + Mu2 - U2)**2 - Mass2**2 )
  CASE (4)
    p2min = SQRT( (max(Mass2,Mu2-Tfac*T-U2))**2 - Mass2**2 )
    p2max = SQRT( (SQRT(p4max**2+Mass4**2)+U4-U2)**2-Mass2**2 )
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
    tmp = p2max
  END IF
  Fmax = SQRT( (tmp+Enu)**2+mpt**2 )+U4-SQRT(tmp**2+Mass2**2)-U2

  IF(Enu>=Fmax) THEN
    p2max = -1000.0d0
    RETURN
  END IF

  F0 = SQRT(Enu**2 + mpt**2)+U4-Mass2-U2
  IF( F0 >= Enu .and. U4 >= U2) RETURN

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
  REAL(DP) :: FD
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  FD = 1.0d0 / (EXP((E - mu)/T) + 1.0d0)
END FUNCTION

END MODULE wlSemiLeptonicOpacityModule2D

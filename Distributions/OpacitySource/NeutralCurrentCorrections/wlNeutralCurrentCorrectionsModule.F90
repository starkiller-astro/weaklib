MODULE wlNeutralCurrentCorrectionsModule

  USE wlKindModule, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeNCWeakMagnetismCorrection
  PUBLIC :: ComputeNCManyBodyCorrection

  INTERFACE ComputeNCManyBodyCorrection
    MODULE PROCEDURE ComputeNCManyBodyCorrection_Scalar
    MODULE PROCEDURE ComputeNCManyBodyCorrection_Vector
  END INTERFACE 

  CONTAINS


  SUBROUTINE ComputeNCWeakMagnetismCorrection &
    (iE_B, iE_E, ENu, Xi_Nu_N, Xi_Nu_P, Xi_NuBar_N, Xi_NuBar_P )

  IMPLICIT none


  INTEGER,  INTENT(in)  :: iE_B, iE_E
  REAL(DP), INTENT(in)  :: ENu       (iE_B:iE_E)
  REAL(DP), INTENT(out) :: Xi_Nu_N   (iE_B:iE_E)
  REAL(DP), INTENT(out) :: Xi_Nu_P   (iE_B:iE_E)
  REAL(DP), INTENT(out) :: Xi_NuBar_N(iE_B:iE_E)
  REAL(DP), INTENT(out) :: Xi_NuBar_P(iE_B:iE_E)

  REAL(DP), PARAMETER   :: twothd    = 2.0d0/3.0d0

  REAL(DP), PARAMETER   :: sin2W     = 0.23116d+00   ! weak mixing angle

  REAL(DP), PARAMETER   :: mu_p      = 1.79285d0     ! anomalous proton magnetic moment

  REAL(DP), PARAMETER   :: mu_n      = -1.91304d0    ! anomalous neutron magnetic moment

  REAL(DP), PARAMETER   :: ga        = 1.26d+00      ! Gamow-Teller beta decay constant. This is
!                                                       not unity as the oupling for Gamow-Teller
!                                                       decay is renormalized.

  REAL(DP), PARAMETER   :: mb        = 931.494d+00   ! baryon mass [MeV]

  REAL(DP), PARAMETER   :: cv_p      = 0.5d0 - 2.d0 * sin2W    ! vector current neutrino-proton coupling constant

  REAL(DP), PARAMETER   :: cv_n      = - 0.5d0       ! vector current neutrino-neutron coupling constant

  REAL(DP), PARAMETER   :: ca_p      = ga/2.d0       ! axial vector current neutrino-proton coupling constant

  REAL(DP), PARAMETER   :: ca_n      = -ga/2.d0      ! axial vector current neutrino-neutron coupling constant

  REAL(DP), PARAMETER   :: CV_p2     = CV_p * CV_p

  REAL(DP), PARAMETER   :: CV_n2     = CV_n * CV_n

  REAL(DP), PARAMETER   :: CA_p2     = CA_p * CA_p

  REAL(DP), PARAMETER   :: CA_n2     = CA_n * CA_n

  REAL(DP), PARAMETER   :: F2_p      = 0.5d0 * ( mu_p - mu_n ) - 2.d0 * sin2W * mu_p    ! proton form factor

  REAL(DP), PARAMETER   :: F2_n      = -0.5d0 * ( mu_p - mu_n ) - 2.d0 * sin2W * mu_n    ! neutron form factor

  REAL(DP), PARAMETER   :: CV_F2CA_p = (CV_p + F2_p) * CA_p

  REAL(DP), PARAMETER   :: CV_F2CA_n = (CV_n + F2_n) * CA_n

  REAL(DP), PARAMETER   :: CVF2_p    = CV_p * F2_p

  REAL(DP), PARAMETER   :: CVF2_n    = CV_n * F2_n

  REAL(DP), PARAMETER   :: F22_p     = F2_p * F2_p

  REAL(DP), PARAMETER   :: F22_n     = F2_n * F2_n

  REAL(DP), PARAMETER   :: CV5CA_p   = ( 2.d0/3.d0 ) * ( CV_p2 + 5.d0 * CA_p2 )

  REAL(DP), PARAMETER   :: CV5CA_n   = ( 2.d0/3.d0 ) * ( CV_n2 + 5.d0 * CA_n2 )

  INTEGER               :: iE            ! energy index counter

  REAL(DP)              :: e             ! ENu/mc2
  REAL(DP)              :: e2            ! e**2
  REAL(DP)              :: e3            ! e**3
  REAL(DP)              :: zeta          ! 1 + 2e
  REAL(DP)              :: zeta3         ! zeta**3
  REAL(DP)              :: ln_zeta       ! ln(zeta)

  REAL(DP)              :: chi_wm_rec    ! cross section correction due to weak magnetism and recoil
  REAL(DP)              :: chi_wm_rec1   ! cross section correction 1 due to weak magnetism and recoil
  REAL(DP)              :: chi_wm_rec2   ! cross section correction 2 due to weak magnetism and recoil
  REAL(DP)              :: chi_wm_rec3   ! cross section correction 3 due to weak magnetism and recoil
  REAL(DP)              :: chi_wm_rec4   ! cross section correction 4 due to weak magnetism and recoil
  REAL(DP)              :: chi_wm_rec5   ! cross section correction 5 due to weak magnetism and recoil
  REAL(DP)              :: chi_rec       ! cross section correction due to recoil
  REAL(DP)              :: chi_rec1      ! cross section correction 1 due to recoil
  REAL(DP)              :: chi_rec2      ! cross section correction 2 due to recoil

!-----------------------------------------------------------------------
!
!         \\\\\ GENERAL NEUTRINO-NUCLEON WEAK MAGNETISM /////
!         \\\\\       AND RECOIL CORRECTION TEMRS       /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Dimensionless energies
!-----------------------------------------------------------------------

  DO iE = iE_B, iE_E
    e                = ENu(iE)/mb
    e2               = e * e
    e3               = e2 * e
    zeta             = 1.d0 + 2.d0 * e
    zeta3            = zeta * zeta * zeta
    ln_zeta          = DLOG( zeta )

!-----------------------------------------------------------------------
!  chi_wm_rec terms
!-----------------------------------------------------------------------

    chi_wm_rec1      = ( ( e - 1.d0 )/( 2.d0 * e3 ) ) * ln_zeta           &
&                    + ( 3.d0 + 12.d0 * e + 9.d0 * e2 - 10.d0 * e3 )/( 3.d0 * e2 * zeta3 )
    chi_wm_rec2      = ( ( e + 1.d0 )/( 2.d0 * e3 ) ) * ln_zeta           &
&                    - ( 3.d0 + 18.d0 * e + 27.d0 * e2 + 10.d0 * e3 )/( 3.d0 * e2 * zeta3 )
    chi_wm_rec3      = ln_zeta/e2 - ( 2.d0 + 10.d0 * e + 28.d0 * e2/3.d0 )/( e * zeta3 )
    chi_wm_rec4      = ln_zeta/e2 - twothd * ( 3.d0 + 15.d0 * e + 22.d0 * e2 )/( e * zeta3 )
    chi_wm_rec5      = ln_zeta/( 4.d0 * e2 ) + ( - 3.d0 - 15.d0 * e - 22.d0 * e2 + 8.d0 * e3 )/( 6.d0 * e * zeta3 )

!-----------------------------------------------------------------------
!
!     \\\\\ GENERAL NEUTRINO-NUCLEON RECOIL CORRECTION TEMRS /////
!
!-----------------------------------------------------------------------

    chi_rec1         = ( ( e + 1.d0 )/e3 ) * ln_zeta - 2.d0/e2
    chi_rec2         = ( ( - 1.d0 - e + 2.d0 * e2 ) * ln_zeta + 2.d0 * e )/( e3 * zeta )

!-----------------------------------------------------------------------
!
!        \\\\\ NEUTRINO-PROTON WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

    chi_wm_rec       = ( CV_p2 * chi_wm_rec1 + CA_p2 * chi_wm_rec2        &
&                    + CV_F2CA_p * chi_wm_rec3 + CVF2_p * chi_wm_rec4     &
&                    + F22_p * chi_wm_rec5 )/CV5CA_p
    chi_rec          = ( CV_p2 * chi_rec1 + CA_p2 * chi_rec2 )/CV5CA_p
    Xi_Nu_P(iE)      = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!
!       \\\\\ NEUTRINO-NEUTRON WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

    chi_wm_rec       = ( CV_n2 * chi_wm_rec1 + CA_n2 * chi_wm_rec2        &
&                    + CV_F2CA_n * chi_wm_rec3 + CVF2_n * chi_wm_rec4     &
&                    + F22_n * chi_wm_rec5 )/CV5CA_n
    chi_rec          = ( CV_n2 * chi_rec1 + CA_n2 * chi_rec2 )/CV5CA_n
    Xi_Nu_N(iE)      = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!
!      \\\\\ ANTINEUTRINO-PROTON WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

    chi_wm_rec       = ( CV_p2 * chi_wm_rec1 + CA_p2 * chi_wm_rec2        &
&                    - CV_F2CA_p * chi_wm_rec3 + CVF2_p * chi_wm_rec4     &
&                    + F22_p * chi_wm_rec5 )/CV5CA_p
    chi_rec          = ( CV_p2 * chi_rec1 + CA_p2 * chi_rec2 )/CV5CA_p
    Xi_NuBar_P(iE)   = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!
!     \\\\\ ANTINEUTRINO-NEUTRON WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

    chi_wm_rec       = ( CV_n2 * chi_wm_rec1 + CA_n2 * chi_wm_rec2        &
&                    - CV_F2CA_n * chi_wm_rec3 + CVF2_n * chi_wm_rec4     &
&                    + F22_n * chi_wm_rec5 )/CV5CA_n
    chi_rec          = ( CV_n2 * chi_rec1 + CA_n2 * chi_rec2 )/CV5CA_n
    Xi_NuBar_N(iE)   = chi_wm_rec/chi_rec

  END DO

  END SUBROUTINE ComputeNCWeakMagnetismCorrection

  SUBROUTINE ComputeNCManyBodyCorrection_Scalar &
    ( brydns, tmev, Yp, S_tot )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: brydns ! baryon density [fm^{-3}]
    REAL(DP), INTENT(in)  :: tmev   ! temperature [MeV]
    REAL(DP), INTENT(in)  :: Yp     ! proton fraction, needs to become Ye+Ym later!
    REAL(DP), INTENT(out) :: S_tot


    REAL(DP), PARAMETER   :: A_0 = 920.d0
    REAL(DP), PARAMETER   :: B_0 = 3.05d0
    REAL(DP), PARAMETER   :: C_0 = 6140.d0
    REAL(DP), PARAMETER   :: D_0 = 1.5d+13
    REAL(DP), PARAMETER   :: ga = 1.26d0   ! axial charge of the nucleon
    REAL(DP), PARAMETER   :: ga52 = 5.d0 * ga**2
    REAL(DP)              :: A             ! fitting formula parameter
    REAL(DP)              :: B             ! fitting formula parameter
    REAL(DP)              :: C             ! fitting formula parameter
    REAL(DP)              :: S_A           ! fitting formula parameter
    REAL(DP)              :: S_V           ! fitting formula parameter

    A       = A_0 * ( brydns * ( 1.d0 - Yp + Yp**2 ) )/( tmev**1.22d0 )
    B       = b_0/( tmev**0.75d0 )
    C       = C_0 * ( brydns * Yp * ( 1.d0 - Yp ) )/ ( tmev**0.5d0 ) + D_0 * brydns**4/( tmev**6 )
    S_A     = 1.d0/( 1.d0 + A * ( 1.d0 + B * exp(-C) ) )
    S_V     = 1.d0

    S_tot   = ( ga52 * S_A + ( 1.d0 - Yp ) * S_V )/( ga52 + ( 1.d0 - Yp ) )

    RETURN

  END SUBROUTINE ComputeNCManyBodyCorrection_Scalar

  SUBROUTINE ComputeNCManyBodyCorrection_Vector &
    ( iP_B, iP_E, D, T, Yp, S_tot )

    INTEGER,  INTENT(in)  :: iP_B, iP_E
    REAL(DP), INTENT(in)  :: D(iP_B:iP_E), T(iP_B:iP_E), Yp(iP_B:iP_E)
    REAL(DP), INTENT(out) :: S_tot(iP_B:iP_E)

    INTEGER :: iP

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: D, T, Yp ) &
    !$OMP MAP( from: S_tot )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( D, T, Yp ) &
    !$ACC COPYOUT( S_tot )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = iP_B, iP_E

      CALL ComputeNCManyBodyCorrection_Scalar &
             ( D(iP), T(iP), Yp(iP), S_tot(iP) )

    END DO

  END SUBROUTINE ComputeNCManyBodyCorrection_Vector

END MODULE

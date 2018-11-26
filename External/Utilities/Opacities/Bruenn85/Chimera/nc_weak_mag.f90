SUBROUTINE nc_weak_mag( e_nu, xi_p_wm, xi_n_wm, xib_p_wm, xib_n_wm, nez )
!-----------------------------------------------------------------------
!
!    File:         nc_weak_mag
!    Module:       nc_weak_mag
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/05/07
!
!    Purpose:
!      To compute the weak magnetism correction for neutral current
!       neutrio scattering on neutrons and protons
!
!    Subprograms called:
!       none
!
!    Input arguments:
!  e_nu           : neutrino energy [MeV]
!
!    Output arguments:
!  xi_p_wm        : weak magnetism correction for neutrino-proton scattering
!  xi_n_wm        : weak magnetism correction for neutrino-neutron scattering
!  xib_p_wm       : weak magnetism correction for antineutrino-proton scattering
!  xib_n_wm       : weak magnetism correction for antineutrino-neutron scattering
!
!    Modules used:
!  kind_module
!  numerical_module
!  physcnst_module
!
!  prb_cntl_module
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: half, one, twothd
USE physcnst_module, ONLY: sin2W, mu_p, mu_n, mp, mn, ga, mb, &
&              CV_n,CV_p,CA_p,CA_n,F2_p,F2_n,     &
&              CV_p2,CV_n2,CA_p2,CA_n2,           &
&              CV_F2CA_n,CV_F2CA_p,CVF2_p,CVF2_n, &
&              F22_p,F22_n,CV5CA_p,CV5CA_n

USE prb_cntl_module, ONLY: inc_wm

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                       :: nez  ! energy grid size

REAL(double), INTENT(in), DIMENSION(nez)  :: e_nu ! energy grid

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(double), INTENT(out), DIMENSION(nez) :: xi_p_wm       ! weak magnetism correction for neutrino-proton scattering
REAL(double), INTENT(out), DIMENSION(nez) :: xi_n_wm       ! weak magnetism correction for neutrino-neutron scattering
REAL(double), INTENT(out), DIMENSION(nez) :: xib_p_wm      ! weak magnetism correction for antineutrino-proton scattering
REAL(double), INTENT(out), DIMENSION(nez) :: xib_n_wm      ! weak magnetism correction for antineutrino-neutron scattering

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                         :: k             ! energy index counter

REAL(double)                    :: e             ! e_nu/mc2
REAL(double)                    :: e2            ! e**2
REAL(double)                    :: e3            ! e**3
REAL(double)                    :: zeta          ! 1 + 2e
REAL(double)                    :: zeta3         ! zeta**3
REAL(double)                    :: ln_zeta       ! ln(zeta)

REAL(double)                    :: chi_wm_rec    ! cross section correction due to weak magnetism and recoil
REAL(double)                    :: chi_wm_rec1   ! cross section correction 1 due to weak magnetism and recoil
REAL(double)                    :: chi_wm_rec2   ! cross section correction 2 due to weak magnetism and recoil
REAL(double)                    :: chi_wm_rec3   ! cross section correction 3 due to weak magnetism and recoil
REAL(double)                    :: chi_wm_rec4   ! cross section correction 4 due to weak magnetism and recoil
REAL(double)                    :: chi_wm_rec5   ! cross section correction 5 due to weak magnetism and recoil
REAL(double)                    :: chi_rec       ! cross section correction due to recoil
REAL(double)                    :: chi_rec1      ! cross section correction 1 due to recoil
REAL(double)                    :: chi_rec2      ! cross section correction 2 due to recoil

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set weak magnetism corrections to zero if inc_wm = 0
!-----------------------------------------------------------------------

IF ( inc_wm == 0 ) THEN
  xi_p_wm          = one
  xi_n_wm          = one
  xib_p_wm         = one
  xib_n_wm         = one
  RETURN
END IF ! inc_wm == 0

!-----------------------------------------------------------------------
!
!         \\\\\ GENERAL NEUTRINO-NUCLEON WEAK MAGNETISM /////
!         \\\\\       AND RECOIL CORRECTION TEMRS       /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Dimensionless energies
!-----------------------------------------------------------------------

DO k = 1,nez
  e                = e_nu(k)/mb
  e2               = e * e
  e3               = e2 * e
  zeta             = 1.d0 + 2.d0 * e
  zeta3            = zeta * zeta * zeta
  ln_zeta          = DLOG( zeta )

!-----------------------------------------------------------------------
!  chi_wm_rec terms
!-----------------------------------------------------------------------

  chi_wm_rec1      = ( ( e - 1.d0 )/( 2.d0 * e3 ) ) * ln_zeta           &
&                  + ( 3.d0 + 12.d0 * e + 9.d0 * e2 - 10.d0 * e3 )/( 3.d0 * e2 * zeta3 )
  chi_wm_rec2      = ( ( e + 1.d0 )/( 2.d0 * e3 ) ) * ln_zeta           &
&                  - ( 3.d0 + 18.d0 * e + 27.d0 * e2 + 10.d0 * e3 )/( 3.d0 * e2 * zeta3 )
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
&                  + CV_F2CA_p * chi_wm_rec3 + CVF2_p * chi_wm_rec4     &
&                  + F22_p * chi_wm_rec5 )/CV5CA_p
  chi_rec          = ( CV_p2 * chi_rec1 + CA_p2 * chi_rec2 )/CV5CA_p
  xi_p_wm(k)       = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!
!       \\\\\ NEUTRINO-NEUTRON WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

  chi_wm_rec       = ( CV_n2 * chi_wm_rec1 + CA_n2 * chi_wm_rec2        &
&                  + CV_F2CA_n * chi_wm_rec3 + CVF2_n * chi_wm_rec4     &
&                  + F22_n * chi_wm_rec5 )/CV5CA_n
  chi_rec          = ( CV_n2 * chi_rec1 + CA_n2 * chi_rec2 )/CV5CA_n
  xi_n_wm(k)       = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!
!      \\\\\ ANTINEUTRINO-PROTON WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

  chi_wm_rec       = ( CV_p2 * chi_wm_rec1 + CA_p2 * chi_wm_rec2        &
&                  - CV_F2CA_p * chi_wm_rec3 + CVF2_p * chi_wm_rec4     &
&                  + F22_p * chi_wm_rec5 )/CV5CA_p
  chi_rec          = ( CV_p2 * chi_rec1 + CA_p2 * chi_rec2 )/CV5CA_p
  xib_p_wm(k)      = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!
!     \\\\\ ANTINEUTRINO-NEUTRON WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

  chi_wm_rec       = ( CV_n2 * chi_wm_rec1 + CA_n2 * chi_wm_rec2        &
&                  - CV_F2CA_n * chi_wm_rec3 + CVF2_n * chi_wm_rec4     &
&                  + F22_n * chi_wm_rec5 )/CV5CA_n
  chi_rec          = ( CV_n2 * chi_rec1 + CA_n2 * chi_rec2 )/CV5CA_n
  xib_n_wm(k)      = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

END DO ! k = 1,nez

RETURN
END SUBROUTINE nc_weak_mag

SUBROUTINE cc_weak_mag_weaklib( e_nu, xi_n_wm, xib_p_wm, nez )
!-----------------------------------------------------------------------
!
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/05/07
!
!    Purpose:
!      To compute the weak magnetism correction for neutrino absorption
!       on neutrons and antineutrino absorption on protons.
!
!    Subprograms called:
!       none
!
!    Input arguments:
!  e_nu           : neutrino energy [MeV]
!  nez            : number of energy groups
!
!    Output arguments:
!  xi_n_wm        : weak magnetism correction for neutrino absorption on neutrons
!  xib_p_wm       : weak magnetism correction for antineutrino absorption on protons
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: half, one, third
USE physcnst_module, ONLY: sin2W, mu_p, mu_n, mp, mn, ga

USE prb_cntl_module, ONLY: icc_wm

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nez          ! number of energy groups
REAL(dp), INTENT(in), DIMENSION(nez) :: e_nu ! neutrino energy [MeV]

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(out), DIMENSION(nez) :: xi_n_wm  ! weak magnetism correction for neutrino absorption on neutrons
REAL(dp), INTENT(out), DIMENSION(nez) :: xib_p_wm ! weak magnetism correction for antineutrino absorption on protons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                         :: k             ! energy zone counter

REAL(dp), PARAMETER         :: CV      = 1.d0    ! charged vector current coupling constant
REAL(dp), PARAMETER         :: CA      = 1.26    ! charged axial vector current coupling constant
REAL(dp), PARAMETER         :: F2      = 3.706d0 ! neutron-proton form factor
REAL(dp), PARAMETER         :: CV_2    = CV * CV
REAL(dp), PARAMETER         :: CA_2    = CA * CA
REAL(dp), PARAMETER         :: CV_F2CA = (CV + F2)*CA
REAL(dp), PARAMETER         :: CVF2    = CV * F2
REAL(dp), PARAMETER         :: F22     = F2 * F2
REAL(dp), PARAMETER         :: CV3CA   = ( CV_2 + 3.d0 * CA_2 )
REAL(dp), PARAMETER         :: m       = 0.5d0 * ( mp + mn )

REAL(dp)                    :: e             ! e_nu/mc2
REAL(dp)                    :: e2            ! e**2
REAL(dp)                    :: zeta          ! 1 + 2e
REAL(dp)                    :: zeta3         ! zeta**3
REAL(dp)                    :: ln_zeta       ! ln(zeta)

REAL(dp)                    :: chi_wm_rec    ! cross section correction due to weak magnetism and recoil
REAL(dp)                    :: chi_wm_rec1   ! cross section correction 1 due to weak magnetism and recoil
REAL(dp)                    :: chi_wm_rec2   ! cross section correction 2 due to weak magnetism and recoil
REAL(dp)                    :: chi_wm_rec3   ! cross section correction 3 due to weak magnetism and recoil
REAL(dp)                    :: chi_wm_rec4   ! cross section correction 4 due to weak magnetism and recoil
REAL(dp)                    :: chi_wm_rec5   ! cross section correction 5 due to weak magnetism and recoil
REAL(dp)                    :: chi_rec       ! cross section correction due to recoil
REAL(dp)                    :: chi_rec1      ! cross section correction 1 due to recoil
REAL(dp)                    :: chi_rec2      ! cross section correction 2 due to recoil

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set weak magnetism corrections to zero if inc_wm = 0
!-----------------------------------------------------------------------

!IF ( icc_wm == 0 ) THEN
!  xi_n_wm          = one
!  xib_p_wm         = one
!  RETURN
!END IF ! icc_wm == 0

!-----------------------------------------------------------------------
!
!         \\\\\ GENERAL NEUTRINO-NUCLEON WEAK MAGNETISM /////
!         \\\\\       AND RECOIL CORRECTION TEMRS       /////
!
!-----------------------------------------------------------------------

DO k = 1,nez

!-----------------------------------------------------------------------
!  Dimensionless energies
!-----------------------------------------------------------------------

  e                = e_nu(k)/m
  e2               = e * e
  zeta             = 1.d0 + 2.d0 * e
  zeta3            = zeta * zeta * zeta
  ln_zeta          = DLOG( zeta )
 
!-----------------------------------------------------------------------
!  chi_wm_rec terms
!-----------------------------------------------------------------------

  chi_wm_rec1      = 1.d0 + 4.d0 * e + 16.d0 * e2 * third
  chi_wm_rec2      = 3.d0 * ( 1.d0 + 4.d0 * e * third ) * ( 1.d0 + 4.d0 * e * third )
  chi_wm_rec3      = 4.d0 * e * ( 1.d0 + 4.d0 * e * third )
  chi_wm_rec4      = 8.d0 * third * e2
  chi_wm_rec5      = 5.d0 * third * e2 * ( 1.d0 + 2.d0 * e/5.d0 )

!-----------------------------------------------------------------------
!
!     \\\\\ GENERAL NEUTRINO-NUCLEON RECOIL CORRECTION TEMRS /////
!
!-----------------------------------------------------------------------

  chi_rec1         = 1.d0/e - ( 1.d0/( 2.d0 * e2 ) ) * ln_zeta
  chi_rec2         = ( zeta * ln_zeta - 2.d0 * e + 4.d0 * e2 )/( 2.d0 * e2 * zeta )

!-----------------------------------------------------------------------
!
!       \\\\\ NEUTRINO ABSORPTION WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

  chi_wm_rec       = ( CV_2 * chi_wm_rec1 + CA_2 * chi_wm_rec2        &
&                  + CV_F2CA * chi_wm_rec3 + CVF2 * chi_wm_rec4     &
&                  + F22 * chi_wm_rec5 )/( zeta3 * CV3CA )
  chi_rec          = ( CV_2 * chi_rec1 + CA_2 * chi_rec2 )/CV3CA
  xi_n_wm(k)       = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!
!     \\\\\ ANTINEUTRINO ABSORPTION WEAK MAGNETISM CORRECTION /////
!
!-----------------------------------------------------------------------

  chi_wm_rec       = ( CV_2 * chi_wm_rec1 + CA_2 * chi_wm_rec2        &
&                  - CV_F2CA * chi_wm_rec3 + CVF2 * chi_wm_rec4     &
&                  + F22 * chi_wm_rec5 )/( zeta3 * CV3CA )
  chi_rec          = ( CV_2 * chi_rec1 + CA_2 * chi_rec2 )/CV3CA
  xib_p_wm(k)      = chi_wm_rec/chi_rec

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

END DO ! k = 1,nez

RETURN
END SUBROUTINE cc_weak_mag_weaklib

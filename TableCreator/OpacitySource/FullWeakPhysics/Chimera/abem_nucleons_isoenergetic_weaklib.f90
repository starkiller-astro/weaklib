SUBROUTINE abem_nucleons_isoenergetic_weaklib &
           ( n, e_in, rho, t, ye, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
             cmpe, absornp, emitnp, nez)
!-----------------------------------------------------------------------
!
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      10/23/18
!
!    Purpose:
!      To call calculate the neutrino absorption and emission rates.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!  abemfrnpetanp : calculates rates of neutrino absorption and emission
!                   when chemical potential data are available
!  abemfrnp      : calculates rates of neutrino absorption and emission
!                   when chemical potential data are not available
!
!    Input arguments:
!  n             : neutrino type
!                  (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  e_in          : neutrino energy [MeV]
!  rho           : matter density [g cm^{-3}]
!  t             : matter temperature [K]
!  ye            : matter electron fraction
!  xneut         : free neutron mass fraction
!  xprot         : free proton mass fraction
!  xh            : heavy nucleus mass fraction
!  ah            : heavy nucleus mass number
!  zh            : heavy nucleus charge number
!  cmpn          : free neutron chemical potential
!                  (excluding rest mass) [MeV]
!  cmpp          : free proton chemical potential
!                  (excluding rest mass) [MeV]
!  cmpe          : electron chemical potential
!                  (including rest mass) [MeV]
!
!  EmAb_nucleons_recoil : Include corrections to nucleons EmAb due to Reddy98
!  EmAb_nucleons_weak_magnetism : Include weak magnetism corrections to EmAb on free nucleons
!
!    Output arguments:
!  absornp       : inverse mean free path for absorption on free
!                  nucleons (/cm)
!  emitnp        : inverse mean free path for emission from free
!                  nucleons (/cm)
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: &
      zero, one, epsilon, pi
USE physcnst_module, ONLY: &
      cvel, hbar, gv, ga, mn, mp, me, Gw, kmev, dmnp, rmu

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)          :: n             ! neutrino flavor index
INTEGER, INTENT(in)          :: nez           ! number of energy groups

REAL(dp), DIMENSION(nez), INTENT(in) :: e_in
                            ! zone centered incoming neutrino energy [MeV]
REAL(dp), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(dp), INTENT(in)    :: t             ! temperature [K]
REAL(dp), INTENT(in)    :: ye            ! electron fraction
REAL(dp), INTENT(in)    :: xneut         ! free neutron mass fraction
REAL(dp), INTENT(in)    :: xprot         ! free proton mass fraction
REAL(dp), INTENT(in)    :: xh            ! heavy nuclei mass fraction
REAL(dp), INTENT(in)    :: ah            ! heavy nuclei mass number
REAL(dp), INTENT(in)    :: zh            ! heavy nuclei charge number
REAL(dp), INTENT(in)    :: cmpn          ! neutron chemical porential
REAL(dp), INTENT(in)    :: cmpp          ! proton chemical porential
REAL(dp), INTENT(in)    :: cmpe          ! electron chemical porential

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), DIMENSION(nez), INTENT(out) :: absornp
                 ! inverse mean free path for absorption on free nucleons
REAL(dp), DIMENSION(nez), INTENT(out) :: emitnp
                 ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

!i_abemetanp        = .true.

INTEGER                     :: k             ! energy group index

LOGICAL                     :: i_abemetanp

REAL(dp)                :: xneutp        ! free neutron mass fraction
REAL(dp)                :: xprotp        ! free proton mass fraction
REAL(dp)                :: xhe           ! helium mass fraction
REAL(dp)                :: yep           ! electron fraction

REAL(dp)                :: tmev          ! temperature [MeV]
REAL(dp)                :: m_trgt_i      ! mass of the initial target particle [MeV]
REAL(dp)                :: m_trgt_f      ! mass of the final target particle [MeV]
REAL(dp)                :: m_lep         ! mass of the final lepton [MeV]
REAL(dp)                :: cmp_trgt_i    ! chemical potential of the initial target particle [MeV]
REAL(dp)                :: cmp_trgt_f    ! chemical potential of the transformed target particle [MeV]
REAL(dp)                :: cmp_lep       ! chemcal potential of the secondary lepton [MeV]
REAL(dp)                :: ab_r0_nu      ! zero moment of he inverse mean free path for neutrino absorption on free neutrons
REAL(dp)                :: ab_r1_nu      ! first moment of he inverse mean free path for neutrino absorption on free neutrons
REAL(dp)                :: e_out_e       ! mean energy of the emitted electron
REAL(dp)                :: ab_r0_nub     ! zero moment of he inverse mean free path for antineutrino absorption on free protons
REAL(dp)                :: ab_r1_nub     ! first moment of he inverse mean free path for antineutrino absorption on free protons
REAL(dp)                :: e_out_p       ! mean energy of the emitted positron
REAL(dp)                :: etam          ! - ( e_in + dmnp + cmpn - cmpp - cmpe )/tmev
REAL(dp)                :: etap          ! - ( e_in + cmpp + cmpe - dmnp - cmpn )/tmev

REAL(dp)                :: xi_n_wm       ! weak magnetism correction for antineutrino absorption on neutrons
REAL(dp)                :: xib_p_wm      ! weak magnetism correction for antineutrino-proton scattering

REAL(dp), PARAMETER     :: x_min = 1.d-30 ! minimum mass fraction fraction
REAL(dp), PARAMETER     :: rho_etanp = 1.d+10 ! minimum density to comopute approximate nucleon blocking factors

REAL(dp)                :: fexp

EXTERNAL fexp

!-----------------------------------------------------------------------
!
!             \\\\\ ADJUST XNEUT AND XPROT IF NECESSARY ////
!
!  Adjust xneut and xprot if ye > yep or if ye < 0.03 to cover regions
!   out of range of the EOS.
!
!-----------------------------------------------------------------------

xneutp             = xneut
xprotp             = xprot
i_abemetanp        = .true.

IF ( ye > 0.5d0 ) THEN

  yep              = DMIN1( ye, 0.99d0 )
  xhe              = DMAX1( one - xneut - xprot - xh, zero )
  i_abemetanp      = .false.

!-----------------------------------------------------------------------
!
!  Are there enough nucleons to shuffle between neutrons and protons
!   to reproduce yep?
!  If so, adjust xneut and xprot to reproduce yep.
!
!-----------------------------------------------------------------------

  IF ( xneut + xprot + 0.5d0 * xhe + xh * zh/( ah + epsilon ) > yep ) THEN

    xprotp         = DMAX1( yep - 0.5d0 * xhe - xh * zh/( ah + epsilon ), x_min )
    xneutp         = DMAX1( xneut + xprot - xprotp, x_min )

  ELSE ! if not, do the best possible

    xprotp         = xneut + xprot
    xneutp         = 1.d-30

  END IF ! xneut + xprot + 0.5d0 * xhe + xh * zh/( ah + epsilon ) > yep

END IF ! ye > 0.5d0

IF ( ye < 0.03d0 ) THEN
  xprotp           = ye
  xneutp           = one - ye
  i_abemetanp      = .false.
END IF ! ye < 0.03d0

!-----------------------------------------------------------------------
!
!           \\\\\ COMPUTE EMISSION AND ABSORPTION RATES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  IF rho < 1.d+10 use abemfrnp instead of abemfrnpetanp to avoid
!   potential numerical problems
!-----------------------------------------------------------------------

  IF ( rho < rho_etanp ) i_abemetanp = .false.

  IF ( i_abemetanp ) THEN

!-----------------------------------------------------------------------
!  Neutrino and antineutrino absorption on neutrons and protons with
!   approximate nucleon blocking but no recoil or thermal motions
!-----------------------------------------------------------------------

  DO k = 1,nez
    CALL abemfrnpetanp &
         ( n, e_in(k), rho, t, xneutp, xprotp, cmpn, cmpp, &
           cmpe, absornp(k), emitnp(k) )

  END DO

  ELSE

!-----------------------------------------------------------------------
!  Neutrino and antineutrino absorption on neutrons and protons with
!   no recoil, thermal motions, or nucleon blocking
!-----------------------------------------------------------------------

  DO k = 1,nez

    CALL abemfrnp &
         ( n, e_in(k), rho, t, xneutp, xprotp, cmpe, &
           absornp(k), emitnp(k) )

  END DO
  END IF ! i_abemetanp

RETURN

  END SUBROUTINE abem_nucleons_isoenergetic_weaklib

SUBROUTINE abem_cal( n, e_in, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
& cmpe, absornp, emitnp, ye_cube, nez, nse, eos )
!-----------------------------------------------------------------------
!
!    File:         abem_cal
!    Module:       abem_cal
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/06
!
!    Purpose:
!      To call calculate the neutrino absorption and emission rates.
!      (1) No recoil or thermal motions: Subroutine abemfrnpetanp is called if
!       ye is <= 0.5. If ye > 0.5, the free neutron and proton mass fractions
!       are adjusted to reflect ye > 0.5 (the LS EOS cannot handle ye > 0.5),
!       and subroutine abemfrnp is called to compute the absorption and
!       emission rates.
!      (2) Recoil and thermal motions Included: Subroutine nu_N_absr_momts is
!       called. If ye > 0.5, the free neutron and proton mass fractions
!       are adjusted to reflect ye > 0.5, and neutron and proton chemical
!       potentials are computed.
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
!  n             : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  e_in          : neutrino energy [MeV]
!  rho           : matter density [g cm^{-3}]
!  t             : matter temperature [K]
!  xneut         : free neutron mass fraction
!  xprot         : free proton mass fraction
!  xh            : heavy nucleus mass fraction
!  ah            : heavy nucleus mass number
!  zh            : heavy nucleus charge number
!  cmpn          : free neutron chemical potential (excluding rest mass) [MeV]
!  cmpp          : free proton chemical potential (excluding rest mass) [MeV]
!  cmpe          : electron chemical potential (including rest mass) [MeV]
!  ye_cube       : electron fraction at the cube corner
!
!    Output arguments:
!  absor         : absorption inverse mean free path (/cm)
!  emitnp        : emission inverse mean free path (/cm)!  
!
!    Modules used:
!  kind_module
!  numerical_module
!  physcnst_module,
! 
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, one, epsilon, pi
USE physcnst_module, ONLY: cvel, hbar, gv, ga, mn, mp, me, Gw, kmev,    &
& dmnp, rmu

USE prb_cntl_module, ONLY: i_aeps

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)          :: n             ! neutrino flavor index
INTEGER, INTENT(in)          :: nez           ! number of energy groups
INTEGER, INTENT(in)          :: nse           ! NSE flag
CHARACTER(LEN=1), INTENT(in) :: eos           ! EoS flag

REAL(double), DIMENSION(nez), INTENT(in) :: e_in ! zone centered incoming neutrino energy [MeV]
REAL(double), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(double), INTENT(in)    :: t             ! temperature [K]
REAL(double), INTENT(in)    :: ye_cube       ! electron fraction at the cube corner
REAL(double), INTENT(in)    :: xneut         ! free neutron mass fraction
REAL(double), INTENT(in)    :: xprot         ! free proton mass fraction
REAL(double), INTENT(in)    :: xh            ! heavy nuclei mass fraction
REAL(double), INTENT(in)    :: ah            ! heavy nuclei mass number
REAL(double), INTENT(in)    :: zh            ! heavy nuclei charge number
REAL(double), INTENT(in)    :: cmpn          ! neutron chemical porential
REAL(double), INTENT(in)    :: cmpp          ! proton chemical porential
REAL(double), INTENT(in)    :: cmpe          ! electron chemical porential

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(double), DIMENSION(nez), INTENT(out) :: absornp ! inverse mean free path for absorption on free nucleons
REAL(double), DIMENSION(nez), INTENT(out) :: emitnp  ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                    :: i_abemetanp

INTEGER                    :: k             ! energy group index

REAL(double)               :: xneutp        ! free neutron mass fraction
REAL(double)               :: xprotp        ! free proton mass fraction
REAL(double)               :: xhe           ! helium mass fraction
REAL(double)               :: yep           ! electron fraction

REAL(double)               :: tmev          ! temperature [MeV]
REAL(double)               :: m_trgt_i      ! mass of the initial target particle [MeV]
REAL(double)               :: m_trgt_f      ! mass of the final target particle [MeV]
REAL(double)               :: m_lep         ! mass of the final lepton [MeV]
REAL(double)               :: cmp_trgt_i    ! chemical potential of the initial target particle [MeV]
REAL(double)               :: cmp_trgt_f    ! chemical potential of the transformed target particle [MeV]
REAL(double)               :: cmp_lep       ! chemcal potential of the secondary lepton [MeV]
REAL(double)               :: ab_r0_nu      ! zero moment of he inverse mean free path for neutrino absorption on free neutrons
REAL(double)               :: ab_r1_nu      ! first moment of he inverse mean free path for neutrino absorption on free neutrons
REAL(double)               :: e_out_e       ! mean energy of the emitted electron
REAL(double)               :: ab_r0_nub     ! zero moment of he inverse mean free path for antineutrino absorption on free protons
REAL(double)               :: ab_r1_nub     ! first moment of he inverse mean free path for antineutrino absorption on free protons
REAL(double)               :: e_out_p       ! mean energy of the emitted positron
REAL(double)               :: etam          ! - ( e_in + dmnp + cmpn - cmpp - cmpe )/tmev 
REAL(double)               :: etap          ! - ( e_in + cmpp + cmpe - dmnp - cmpn )/tmev 

REAL(double)               :: xi_n_wm(nez)  ! weak magnetism correction for antineutrino absorption on neutrons
REAL(double)               :: xib_p_wm(nez) ! weak magnetism correction for antineutrino-proton scattering

REAL(double), PARAMETER    :: x_min = 1.d-30 ! minimum mass fraction fraction
REAL(double), PARAMETER    :: rho_etanp = 1.d+10 ! minimum density to comopute approximate nucleon blocking factors

REAL(double), EXTERNAL     :: fexp          ! floored/capped exponential function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( n > 2 ) THEN
  absornp          = zero
  emitnp           = zero
  RETURN
END IF ! n > 2

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

IF ( ye_cube > 0.5d0 .and. eos == 'L' ) THEN
  yep              = DMIN1( ye_cube, 0.99d0 )
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

END IF ! ye_cube > 0.5d0

IF ( ye_cube < 0.03d0 ) THEN
  xprotp           = ye_cube
  xneutp           = one - ye_cube
  i_abemetanp      = .false.
END IF ! ye_cube < 0.03d0

!-----------------------------------------------------------------------
!
!           \\\\\ COMPUTE EMISSION AND ABSORPTION RATES ////
!
!-----------------------------------------------------------------------

IF ( i_aeps == 0  .or.  rho < 1.d+09 .or. nse == 0 ) THEN

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
      CALL abemfrnpetanp( n, e_in(k), rho, t, xneutp, xprotp, cmpn, cmpp, &
&      cmpe, absornp(k), emitnp(k) )
    END DO ! k = 1,nez
  ELSE

!-----------------------------------------------------------------------
!  Neutrino and antineutrino absorption on neutrons and protons with
!   no recoil, thermal motions, or nucleon blocking
!-----------------------------------------------------------------------

    DO k = 1,nez 
      CALL abemfrnp( n, e_in(k), rho, t, xneutp, xprotp, cmpe, absornp(k),&
&      emitnp(k) )
    END DO ! k = 1, nez
  END IF ! i_abemetanp
ELSE ! i_aeps /= 0

!-----------------------------------------------------------------------
!  Neutrino absorption on neutrons with recoil, thermal motions, and
!   nucleon blocking
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( k, tmev, m_trgt_i, m_trgt_f, m_lep, cmp_trgt_i, etap,    &
!$OMP          cmp_trgt_f, cmp_lep, ab_r0_nu, ab_r1_nu, e_out_e, etam,  &
!$OMP          ab_r0_nub, ab_r1_nub, e_out_p ) &
!$OMP SHARED( nez, n, t, e_in, cmpn, cmpp, cmpe, emitnp, absornp )
  IF ( n == 1 ) THEN
    tmev           = kmev * t
    m_trgt_i       = mn
    m_trgt_f       = mp
    m_lep          = me
    cmp_trgt_i     = cmpn + dmnp + mn
    cmp_trgt_f     = cmpp + dmnp + mp
    cmp_lep        = cmpe
!$OMP DO SCHEDULE(static)
    DO k = 1,nez
      CALL nu_N_absr_momts( e_in(k), tmev, m_trgt_i, m_trgt_f, m_lep,   &
&      cmp_trgt_i, cmp_trgt_f, cmp_lep, ab_r0_nu, ab_r1_nu, e_out_e )
      absornp(k)   = ab_r0_nu
      etam         = - ( e_in(k) + dmnp + cmpn - cmpp - cmpe )/tmev
      emitnp(k)    = fexp(etam) * absornp(k)
    END DO ! k = 1,nez
!$OMP END DO
  
!-----------------------------------------------------------------------
!  Antieutrino absorption on protons with recoil, thermal motions, and
!   nucleon blocking
!-----------------------------------------------------------------------

  ELSE IF ( n == 2 ) THEN
    tmev           = kmev * t
    m_trgt_i       = mp
    m_trgt_f       = mn
    m_lep          = me
    cmp_trgt_i     = cmpp + dmnp + mp
    cmp_trgt_f     = cmpn + dmnp + mn
    cmp_lep        = - cmpe
!$OMP DO SCHEDULE(static)
    DO k = 1,nez
      CALL nu_N_absr_momts( e_in(k), tmev, m_trgt_i, m_trgt_f, m_lep,   &
&      cmp_trgt_i, cmp_trgt_f, cmp_lep, ab_r0_nub, ab_r1_nub, e_out_p )
      absornp(k)   = ab_r0_nub
      etap         = - ( e_in(k) + cmpp + cmpe - dmnp - cmpn )/tmev
      emitnp(k)    = fexp(etap) * absornp(k)
    END DO ! k = 1,nez
!$OMP END DO

!-----------------------------------------------------------------------
!  Absorption and emission on free nucleons is zero for mu and tau
!   neutrinos and antineutrinos
!-----------------------------------------------------------------------

  ELSE
    absornp        = zero
    emitnp         = zero
  END IF ! n == 1
!$OMP END PARALLEL

END IF ! i_aeps == 0

!-----------------------------------------------------------------------
!  Weak magnetism corrections for neutrino and antineutrino absorptions
!   on nucleons.
!-----------------------------------------------------------------------

CALL cc_weak_mag( e_in, xi_n_wm, xib_p_wm, nez )

IF ( n == 1 ) THEN
  absornp          = xi_n_wm * absornp
  emitnp           = xi_n_wm * emitnp
ELSE IF ( n == 2 ) THEN
  absornp          = xib_p_wm * absornp
  emitnp           = xib_p_wm * emitnp
END IF ! n == 1

RETURN
END SUBROUTINE abem_cal

SUBROUTINE abem_cal_weaklib &
           ( n, e_in, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
             cmpe, absornp, emitnp, nez )
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
!    Output arguments:
!  absornp       : inverse mean free path for absorption on free
!                  nucleons (/cm)
!  emitnp        : inverse mean free path for emission from free
!                  nucleons (/cm)
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
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

REAL(double), DIMENSION(nez), INTENT(in) :: e_in 
                            ! zone centered incoming neutrino energy [MeV]
REAL(double), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(double), INTENT(in)    :: t             ! temperature [K]
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

REAL(double), DIMENSION(nez), INTENT(out) :: absornp 
                 ! inverse mean free path for absorption on free nucleons
REAL(double), DIMENSION(nez), INTENT(out) :: emitnp 
                 ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                    :: i_abemetanp

INTEGER                    :: k             ! energy group index

REAL(double), PARAMETER    :: x_min = 1.d-30 
                           ! minimum mass fraction fraction
REAL(double), PARAMETER    :: rho_etanp = 1.d+10 
                           ! minimum density to comopute approximate 
                           ! nucleon blocking factors

i_abemetanp        = .true.

!-----------------------------------------------------------------------
!
!           \\\\\ COMPUTE EMISSION AND ABSORPTION RATES ////
!
!-----------------------------------------------------------------------

  IF ( rho < rho_etanp ) i_abemetanp = .false.

  IF ( i_abemetanp ) THEN

!-----------------------------------------------------------------------
!  Neutrino and antineutrino absorption on neutrons and protons with
!   approximate nucleon blocking but no recoil or thermal motions
!-----------------------------------------------------------------------

    DO k = 1,nez

      CALL abemfrnpetanp &
           ( n, e_in(k), rho, t, xneut, xprot, cmpn, cmpp, &
             cmpe, absornp(k), emitnp(k) )

    END DO ! k = 1,nez
  ELSE

!-----------------------------------------------------------------------
!  Neutrino and antineutrino absorption on neutrons and protons with
!   no recoil, thermal motions, or nucleon blocking
!-----------------------------------------------------------------------

    DO k = 1,nez 

      CALL abemfrnp &
      ( n, e_in(k), rho, t, xneut, xprot, cmpe, &
        absornp(k), emitnp(k) )

    END DO ! k = 1, nez
  END IF ! i_abemetanp

  RETURN

  END SUBROUTINE abem_cal_weaklib

SUBROUTINE abemrgn_weaklib &
       ( n, e_in, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
         cmpe, absor, emit, ye, nez )
!--------------------------------------------------------------------
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      10/29/18
!
!    Purpose:
!      To call the subroutines that compute the absorption and
!       emission inverse mean free paths and to store the total
!       inverse mean free paths in arrays, absor and emit,
!       respectively.
!
!    Subprograms called:
!  abem_cal_weaklib
!                : computes neutrino emission and absorption 
!                  on free neutrons and protons
!  abemnc_weaklib
!                : computes neutrino emission and absorption 
!                  on heavy nuclei using the FFN formalizm
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
!  ye            : electron fraction
!  nez           : size of neutrino energy array e_in
!
!    Output arguments:
!  absor         : absorption inverse mean free path (/cm)
!  emit          : emission inverse mean free path (/cm)! 
!--------------------------------------------------------------------

USE kind_module, ONLY: &
      double
USE numerical_module, ONLY: &
      zero, one, epsilon
USE physcnst_module, ONLY: &
      kmev, dmnp
USE prb_cntl_module, ONLY: &
      iaence

IMPLICIT NONE

!--------------------------------------------------------------------
!        Input variables.
!--------------------------------------------------------------------

INTEGER, INTENT(in)          :: n       ! neutrino flavor index
INTEGER, INTENT(in)          :: nez     ! number of energy groups

REAL(double), DIMENSION(nez), INTENT(in) :: e_in 
                      ! zone centered incoming neutrino energy [MeV]
REAL(double), INTENT(in)    :: rho      ! density (g/cm^3)
REAL(double), INTENT(in)    :: t        ! temperature [K]
REAL(double), INTENT(in)    :: ye       ! electron fraction
REAL(double), INTENT(in)    :: xneut    ! free neutron mass fraction
REAL(double), INTENT(in)    :: xprot    ! free proton mass fraction
REAL(double), INTENT(in)    :: xh       ! heavy nuclei mass fraction
REAL(double), INTENT(in)    :: ah       ! heavy nuclei mass number
REAL(double), INTENT(in)    :: zh       ! heavy nuclei charge number
REAL(double), INTENT(in)    :: cmpn     ! neutron chemical porential
REAL(double), INTENT(in)    :: cmpp     ! proton chemical porential
REAL(double), INTENT(in)    :: cmpe     ! electron chemical porential

!--------------------------------------------------------------------
!        Output variables.
!--------------------------------------------------------------------

REAL(double), DIMENSION(nez), INTENT(out) :: absor 
                      ! inverse mean free path for absorption
REAL(double), DIMENSION(nez), INTENT(out) :: emit 
                      ! inverse mean free path for emission 

!--------------------------------------------------------------------
!        Local variables
!--------------------------------------------------------------------
REAL(double), DIMENSION(nez)  :: absrnp, absrnc, emisnp, emisnc

!--------------------------------------------------------------------
!  n-neutrino - free nucleon absorption and emission inverse mean
!   free paths (/cm)
!--------------------------------------------------------------------

  CALL abem_cal_weaklib &
        ( n, e_in, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
          cmpe, absrnp, emisnp, nez )

!--------------------------------------------------------------------
!  n-neutrino - nuclei absorption and emission inverse mean free
!   paths (/cm).
!--------------------------------------------------------------------

  CALL abemnc_weaklib &
        ( n, nez, e_in, rho, t, xh, ah, zh, cmpn, cmpp,  &
          cmpe, absrnc, emisnc )

  absor = absrnp + absrnc
  emit  = emisnp + emisnc

RETURN

END SUBROUTINE abemrgn_weaklib

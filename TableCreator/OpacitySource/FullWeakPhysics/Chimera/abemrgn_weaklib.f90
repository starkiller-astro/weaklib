SUBROUTINE abemrgn_weaklib &
       ( n, e_in, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
         cmpe, absor, emit, ye, nez, & 
         EmAb_nucleons_isoenergetic, EmAb_nuclei_FFN, &
         EmAb_nucleons_recoil, EmAb_nucleons_weak_magnetism )
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
!  ye            : electron fraction
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
!  EmAb_nucleons_isoenergetic : Include isoenergetic EmAb on free nucleons using Bruenn85 formalism
!  EmAb_nuclei_FFN : Include EmAb on nuclei using FFN formalism
!  EmAb_nucleons_recoil : Include corrections to nucleons EmAb due to Reddy98
!  EmAb_nucleons_weak_magnetism : Include weak magnetism corrections to EmAb on free nucleons
!
!    Output arguments:
!  absor         : absorption inverse mean free path (/cm)
!  emit          : emission inverse mean free path (/cm)!
!--------------------------------------------------------------------

USE wlKindModule, ONLY: &
      dp
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

REAL(dp), DIMENSION(nez), INTENT(in) :: e_in 
                      ! zone centered incoming neutrino energy [MeV]
REAL(dp), INTENT(in)    :: rho      ! density (g/cm^3)
REAL(dp), INTENT(in)    :: t        ! temperature [K]
REAL(dp), INTENT(in)    :: ye       ! electron fraction
REAL(dp), INTENT(in)    :: xneut    ! free neutron mass fraction
REAL(dp), INTENT(in)    :: xprot    ! free proton mass fraction
REAL(dp), INTENT(in)    :: xh       ! heavy nuclei mass fraction
REAL(dp), INTENT(in)    :: ah       ! heavy nuclei mass number
REAL(dp), INTENT(in)    :: zh       ! heavy nuclei charge number
REAL(dp), INTENT(in)    :: cmpn     ! neutron chemical porential
REAL(dp), INTENT(in)    :: cmpp     ! proton chemical porential
REAL(dp), INTENT(in)    :: cmpe     ! electron chemical porential

INTEGER, INTENT(in)         :: EmAb_nucleons_isoenergetic   ! Flag to calculate isoenergetic EmAb on free nucleons 
                                                            ! using Bruenn85
INTEGER, INTENT(in)         :: EmAb_nuclei_FFN              ! Flag to calculate EmAb on nuclei using FFN formalism
INTEGER, INTENT(in)         :: EmAb_nucleons_recoil         ! Flag to recoil, nucleon final-state blocking, 
                                                            !and special relativity corrections
INTEGER, INTENT(in)         :: EmAb_nucleons_weak_magnetism ! Flag to include weak_magnetism corrections

!--------------------------------------------------------------------
!        Output variables.
!--------------------------------------------------------------------

REAL(dp), DIMENSION(nez), INTENT(out) :: absor 
                      ! inverse mean free path for absorption
REAL(dp), DIMENSION(nez), INTENT(out) :: emit 
                      ! inverse mean free path for emission 

!--------------------------------------------------------------------
!        Local variables
!--------------------------------------------------------------------
REAL(dp), DIMENSION(nez)  :: ab_nucleons, ab_nuclei, em_nucleons, em_nuclei

!--------------------------------------------------------------------
!  n-neutrino - free nucleon absorption and emission inverse mean
!   free paths (/cm)
!--------------------------------------------------------------------

  IF(EmAb_nucleons_isoenergetic .gt. 0) THEN

    CALL abem_cal_weaklib &
         ( n, e_in, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
           cmpe, ab_nucleons, em_nucleons, nez, EmAb_nucleons_recoil, &
           EmAb_nucleons_weak_magnetism )

  ENDIF

!--------------------------------------------------------------------
!  n-neutrino - nuclei absorption and emission inverse mean free
!   paths (/cm).
!--------------------------------------------------------------------

  IF(EmAb_nuclei_FFN .gt. 0) THEN

    CALL abemnc_weaklib &
         ( n, nez, e_in, rho, t, xh, ah, zh, cmpn, cmpp,  &
           cmpe, ab_nuclei, em_nuclei )

  ENDIF

  absor = ab_nucleons + ab_nuclei
  emit  = em_nucleons + em_nuclei

RETURN

END SUBROUTINE abemrgn_weaklib

PROGRAM wlTestElasticOpacity

  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  USE wlLeptonEOSModule
  USE wlHelmMuonIOModuleHDF
  USE wlEosConstantsModule, ONLY: kmev
  USE wlCalculateAbEmOpacityModule, ONLY: &
    ElasticAbsorptionOpacityNuMu
  IMPLICIT NONE

  INTEGER                        :: i
#ifdef EOSMODE_3D
    TYPE(EquationOfStateTableType) :: EOSTable
#elif defined(EOSMODE_4D)
    TYPE(EquationOfState4DTableType) :: EOSTable
#elif defined(EOSMODE_COMPOSE)
    TYPE(EquationOfStateCompOSETableType) :: EOSTable
#endif
  TYPE(MuonTableType) :: MuonTable
  REAL(DP) :: D, T, Ye, Ym, OpElNu, OpElNuBar
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: E
  REAL(DP) :: Emin, Emax, dE
  INTEGER  :: iE, nE

  CALL ReadEquationOfStateTableHDF( EOSTable, "BaryonsPlusHelmPlusMuonsEOS.h5" )
  CALL ReadMuonTableHDF( MuonTable, "BaryonsPlusHelmPlusMuonsEOS.h5" )

  ! See Fischer 2020 Table II for Thermo State and Figure 1 for opacity
  nE = 100
  Emin = 10.0d0
  Emax = 300.0d0
  dE = (Emax - Emin) / DBLE(nE)

  ALLOCATE( E(nE) )

  DO iE=1,nE
    E(iE) = Emin + (iE-1)*dE
  ENDDO

  ! Case (a) of Table II in FIscher 2020
  D  = 1.0d13
  T  = 10.0d0 / kmev
  Ye = 0.2d0 !Does not matter...
  Ym = 1.0d-4

  ! ! Case (b) of Table II in FIscher 2020
  ! D  = 1.0d14
  ! T  = 25.0d0 / kmev
  ! Ye = 0.15d0 !Does not matter...
  ! Ym = 0.05

  ! Print headers
  WRITE(*,'(A16, A18, A20, A10, A10, A22, A24)') 'Energy (MeV)', 'Density (g/cc)', &
    'Temperature (MeV)', 'Ye', 'Ym', 'Opacity Nu (1/cm)', 'Opacity NuBar (1/cm)'

  DO iE=1,nE
    CALL ElasticAbsorptionOpacityNuMu(D, T, Ye, Ym, E(iE), &
      EOSTable, MuonTable, OpElNu, OpElNuBar)

    ! Print values under each header
    WRITE(*,'(F16.6, ES18.6, F20.6, F10.4, F10.4, ES22.6, ES24.6)') &
      E(iE), D, T*kmev, Ye, Ym, OpElNu, OpElNuBar

  ENDDO

  DEALLOCATE(E)

END PROGRAM wlTestElasticOpacity

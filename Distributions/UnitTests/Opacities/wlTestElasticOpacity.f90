PROGRAM wlTestElasticOpacity

  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  USE wlLeptonEOSTableModule
  USE wlLeptonPhotonGasEOS
  USE wlHelmIOModuleHDF
  USE wlEosConstantsModule, ONLY: &
    rmu, mp, me, mn, mmu, cm3fm3, kmev
  USE wlCalculateAbEmOpacityModule, ONLY: &
    ElasticAbsorptionOpacityNum, &
    ElasticAbsorptionOpacityNue, &
    CalculateHorowitzWeakMagRecoil

  IMPLICIT NONE

  LOGICAL, PARAMETER :: IncludeElasticWeakMagRecoil = .TRUE.
  INTEGER                        :: i
#ifdef EOSMODE_3D
    TYPE(EquationOfStateTableType) :: EOSTable
#elif defined(EOSMODE_4D)
    TYPE(EquationOfState4DTableType) :: EOSTable
#elif defined(EOSMODE_COMPOSE)
    TYPE(EquationOfStateCompOSETableType) :: EOSTable
#endif
  TYPE(HelmTableType) :: HelmTableMuons
  TYPE(LeptonGasType) :: MuonGasState

  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(LeptonGasType) :: ElectronGasState

  REAL(DP) :: D, T, Ye, Ym
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: E
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: OpElNue , OpElNueBar
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: OpElNum, OpElNumBar
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: WeakMagCorrLep, WeakMagCorrLepBar
  REAL(DP) :: Emin, Emax, dE
  INTEGER  :: iE, nE

  ! LOCAL VARIABLES
  REAL(DP) :: Yp, Xp, Xn, dYp
  REAL(DP) :: Mup, Mun, Mumu, Mue
  REAL(DP) :: Mp_eff, Mn_eff, Up, Un, n_n, n_p
  INTEGER  :: iEOS_Rho, iEOS_T, iEOS_Yp, iDV
  
  CALL ReadEquationOfStateTableHDF( EOSTable, "BaryonsPlusPhotonsPlusLeptonsEOS.h5" )
  CALL ReadHelmholtzTableHDF( HelmTableElectrons, "BaryonsPlusPhotonsPlusLeptonsEOS.h5", "HelmTableElectrons" )
  CALL ReadHelmholtzTableHDF( HelmTableMuons, "BaryonsPlusPhotonsPlusLeptonsEOS.h5", "HelmTableMuons" )

  ! See Fischer 2020 Table II for Thermo State and Figure 1 for opacity
  nE = 100
  Emin = 10.0d0
  Emax = 300.0d0
  dE = (Emax - Emin) / DBLE(nE)

  ALLOCATE( E(nE) )
  ALLOCATE( OpElNue(nE) )
  ALLOCATE( OpElNum(nE) )
  ALLOCATE( OpElNueBar(nE) )
  ALLOCATE( OpElNumBar(nE) )
  ALLOCATE( WeakMagCorrLep(nE) )
  ALLOCATE( WeakMagCorrLepBar(nE) )

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

  ! Calculate all relevant quantities
  iEOS_Rho = EOSTable % TS % Indices % iRho
  iEOS_T   = EOSTable % TS % Indices % iT
  iEOS_Yp  = EOSTable % TS % Indices % iYe ! Ye is actually Yp, but for consistency reasons it was not changed

  Yp = Ye + Ym

  iDV = EOSTable % DV % Indices % iNeutronChemicalPotential
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Mun )

  iDV = EOSTable % DV % Indices % iProtonChemicalPotential
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Mup )

#ifdef EOSMODE_3D
  Mn_eff = 0.0d0
  Mp_eff = 0.0d0
  Up     = 0.0d0
  Un     = 0.0d0
#else
  iDV = EOSTable % DV % Indices % iNeutronEffMass
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Mn_eff )

  iDV = EOSTable % DV % Indices % iNeutronSelfEnergy
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Un )

  iDV = EOSTable % DV % Indices % iProtonEffMass
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Mp_eff )

  iDV = EOSTable % DV % Indices % iProtonSelfEnergy
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Up )
#endif

  iDV = EOSTable % DV % Indices % iProtonMassFraction
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Xp )

  iDV = EOSTable % DV % Indices % iNeutronMassFraction
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T)   % Values,   &
        EOSTable % TS % States(iEOS_Yp)  % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Xn )


  ElectronGasState % t   = T
  ElectronGasState % rho = D
  ElectronGasState % yL  = Ye
  CALL LeptonGasEOS(HelmTableElectrons, ElectronGasState)
  Mue = ElectronGasState % mu

  MuonGasState % t   = T
  MuonGasState % rho = D
  MuonGasState % yL  = Ym
  CALL LeptonGasEOS(HelmTableMuons, MuonGasState)
  Mumu = MuonGasState % mu

  WRITE(*,*) 't     =', T * kmev
  WRITE(*,*) 'ye    =', ye
  WRITE(*,*) 'ym    =', ym
  WRITE(*,*) 'd     =', d
  WRITE(*,*) 'yp    =', yp
  WRITE(*,*) 'chemmu=', Mumu
  WRITE(*,*) 'cheme =', Mue
  WRITE(*,*) 'chemn =', Mun - mn
  WRITE(*,*) 'chemp =', Mup - mp
  WRITE(*,*) 'xn    =', xn
  WRITE(*,*) 'xp    =', xp
  WRITE(*,*) 'un    =', un
  WRITE(*,*) 'up    =', up
  WRITE(*,*) 'massn =', Mn_eff
  WRITE(*,*) 'massp =', Mp_eff
  WRITE(*,*) 'nB    =', D * cm3fm3 / rmu
  WRITE(*,*) 'phin  =', Mun - Mn_eff - Un
  WRITE(*,*) 'phip  =', Mup - Mp_eff - Up

  DO iE=1,nE
    ! Initialize Recoil correction!
    IF (IncludeElasticWeakMagRecoil) THEN
      CALL CalculateHorowitzWeakMagRecoil(E(iE), WeakMagCorrLep(iE), WeakMagCorrLepBar(iE))
    ELSE
      WeakMagCorrLep(iE)    = 1.0d0
      WeakMagCorrLepBar(iE) = 1.0d0
    ENDIF

    CALL ElasticAbsorptionOpacityNue (D, T, E(iE), &
      Mun, Mn_eff, Un, Xn, Mup, Mp_eff, Up, Xp, Mue , &
      OpElNue(iE) , OpElNueBar(iE))

      OpElNue(iE)    = OpElNue(iE)    * WeakMagCorrLep(iE)
      OpElNueBar(iE) = OpElNueBar(iE) * WeakMagCorrLepBar(iE)

    CALL ElasticAbsorptionOpacityNum(D, T, E(iE), &
      Mun, Mn_eff, Un, Xn, Mup, Mp_eff, Up, Xp, Mumu, &
      OpElNum(iE), OpElNumBar(iE))

      OpElNum(iE)    = OpElNum(iE)    * WeakMagCorrLep(iE)
      OpElNumBar(iE) = OpElNumBar(iE) * WeakMagCorrLepBar(iE)

  ENDDO

  ! Print headers
  WRITE(*,*)
  WRITE(*,'(A13, A14, A11, A9, A9, A20, A23, A20, A23)') 'Energy (MeV)', 'Rho (g/cc)', &
    'T (MeV)', 'Ye', 'Ym', 'Opacity Nue (1/cm)', 'Opacity NueBar (1/cm)', &
    'Opacity Num (1/cm)', 'Opacity NumBar (1/cm)'
  DO iE=1,nE
    ! Print values under each header
    WRITE(*,'(F13.6, ES14.6, F11.6, F9.4, F9.4, ES20.8, ES23.8, ES20.8, ES23.8)') &
      E(iE), D, T*kmev, Ye, Ym, OpElNue(iE), OpElNueBar(iE), OpElNum(iE), OpElNumBar(iE)
  ENDDO

  DEALLOCATE(E)

END PROGRAM wlTestElasticOpacity

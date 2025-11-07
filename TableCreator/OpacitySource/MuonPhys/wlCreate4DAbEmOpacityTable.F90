PROGRAM wlCreate4DAbEmOpacityTable

  USE wlKindModule, ONLY: dp  
  USE HDF5
  USE wlEosConstantsModule, ONLY: kmev
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  USE wlLeptonEOSTableModule
  USE wlLeptonPhotonGasEOS
  USE wlHelmIOModuleHDF
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
  Opacity_CC_2D
  USE wlCalculateAbEmOpacityModule, ONLY: &
    ElasticAbsorptionOpacityNum, &
    ElasticAbsorptionOpacityNue, &
    CalculateHorowitzWeakMagRecoil
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DeAllocateOpacityTable
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlOpacityFieldsModule, ONLY: &
    iNu_e, iNu_e_bar

  implicit none

  INTEGER, PARAMETER  ::  nE_2D = 40

  INTEGER :: i, j, k, l, nThermoPoints, iE2, iE3

  INTEGER :: icount, nCount
  REAL(DP), allocatable :: OpaA_2D(:,:,:,:,:,:), WeakMagCorrLep(:), WeakMagCorrLepBar(:)
  REAL(DP), allocatable :: T(:), Rho(:), Ye(:), Ym(:), E(:)
  INTEGER :: nE, nRho, nT, nYe, nYm, iE, iRho, iT, iYe, iYm, ierr
  INTEGER :: iEOS_Rho, iEOS_T, iEOS_Yp, iRho_Op, iT_Op, iYe_Op
  REAL(DP) :: xYp, xMumu, xMue, xMun, xMup, xXn, xXp, xUn, xUp, xMn_eff, xMp_eff, xMl, xMul

#ifdef EOSMODE_3D
    TYPE(EquationOfStateTableType) :: EOSTable
#elif defined(EOSMODE_4D)
    TYPE(EquationOfState4DTableType) :: EOSTable
#elif defined(EOSMODE_COMPOSE)
    TYPE(EquationOfStateCompOSETableType) :: EOSTable
#endif
  TYPE(HelmTableType) :: HelmTableMuons
  TYPE(HelmTableType) :: HelmTableElectrons

  TYPE(OpacityTableType) :: OpacityTable
  CHARACTER(len=128) :: FileName_EmAb

  CALL MPI_INIT( ierr )

  CALL InitializeHDF( )
  CALL ReadEquationOfStateTableHDF( EOSTable, "BaryonsPlusPhotonsPlusLeptonsEOS.h5" )
  CALL ReadHelmholtzTableHDF( HelmTableElectrons, "BaryonsPlusPhotonsPlusLeptonsEOS.h5", "HelmTableElectrons" )
  CALL ReadHelmholtzTableHDF( HelmTableMuons, "BaryonsPlusPhotonsPlusLeptonsEOS.h5", "HelmTableMuons" )

  iEOS_Rho = EOSTable % TS % Indices % iRho
  iEOS_T   = EOSTable % TS % Indices % iT
  iEOS_Yp  = EOSTable % TS % Indices % iYe 

  ! Read Opacity Table
  FileName_EmAb = "wl-Op-SFHo-15-25-50-E40-EmAb.h5"
  CALL ReadOpacityTableHDF( OpacityTable,   &
       FileName_EmAb_Option                 &
       = FileName_EmAb, &
       EquationOfStateTableName_Option      &
       = "BaryonsPlusPhotonsPlusLeptonsEOS.h5",  &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )
  WRITE(*,*) 'Done Reading Tables'

  iRho_Op   = OpacityTable % TS % Indices % iRho
  iT_Op     = OpacityTable % TS % Indices % iT
  iYe_Op    = OpacityTable % TS % Indices % iYe
            
  nE = SIZE( OpacityTable % EnergyGrid % Values )
  nRho = SIZE( OpacityTable % TS % States(iRho_Op) % Values )
  nT = SIZE( OpacityTable % TS % States(iT_Op) % Values )
  nYe = SIZE( OpacityTable % TS % States(iYe_Op) % Values )
  nYm = 30

  ALLOCATE( T(nT) )
  ALLOCATE( Rho(nRho) )
  ALLOCATE( Ye(nYe) )
  ALLOCATE( Ym(nYm) )
  ALLOCATE( E(nE), WeakMagCorrLep(nE), WeakMagCorrLepBar(nE) )

  T = OpacityTable % TS % States(iT_Op) % Values
  Rho = OpacityTable % TS % States(iRho_Op) % Values
  Ye = OpacityTable % TS % States(iYe_Op) % Values
  Ym = logspace_array(1.0d-5, 1.0d-1, nYm)
  E = OpacityTable % EnergyGrid % Values

  WRITE(*,*) Ym
  ALLOCATE( OpaA_2D(nYm, nYe, nT, nRho, nE, 2) )

  DO iE = 1, nE
    CALL CalculateHorowitzWeakMagRecoil(E(iE), WeakMagCorrLep(iE), WeakMagCorrLepBar(iE))
  END DO

  iCount = 0
  nCount = nT * nRho * nYe * nYm * nE
  !$OMP PARALLEL DO COLLAPSE(5) &
  !$OMP PRIVATE(iT, iRho, iYe, iYm, iE, xMumu, xMue, xMun, &
  !$OMP xMup, xXn, xXp, xUn, xUp, xMn_eff, xMp_eff)
  DO iT=1, nT
  DO iRho=1, nRho
  DO iYe=1, nYe
  DO iYm=1, nYm
  DO iE=1, nE
    CALL ApplyEOS(T(iT), Rho(iRho), Ye(iYe), Ym(iYm), &
      xMumu, xMue, xMun, xMup, xXn, xXp, xUn, xUp, xMn_eff, xMp_eff)

    CALL Opacity_CC_2D(3, 1, E(iE), OpaA_2D(iYm, iYe, iT, iRho, iE, 1), &
          T(i)*kmev, xMul, xMun, xMup, xMl, xMn_eff, xMp_eff, xUn, xUp, nE_2D)

    CALL Opacity_CC_2D(3, 2, E(iE), OpaA_2D(iYm, iYe, iT, iRho, iE, 2), &
          T(i)*kmev, xMul, xMun, xMup, xMl, xMn_eff, xMp_eff, xUn, xUp, nE_2D)

    OpaA_2D(iYm, iYe, iT, iRho, iE, 1) = OpaA_2D(iYm, iYe, iT, iRho, iE, 1) * WeakMagCorrLep(iE)
    OpaA_2D(iYm, iYe, iT, iRho, iE, 2) = OpaA_2D(iYm, iYe, iT, iRho, iE, 2) * WeakMagCorrLepBar(iE)

    IF ( MOD(iCount, 10000) .EQ. 0 ) THEN
      WRITE(*,*) iCount, '/', nCount
    ENDIF
    iCount = iCount + 1

  END DO
  END DO
  END DO
  END DO
  END DO
  !$OMP END PARALLEL DO

  OPEN(200, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_2D
  CLOSE(200)

  DEALLOCATE(OpaA_2D)
  DEALLOCATE(T, Rho, Ye, Ym, E)

CONTAINS 

pure function logspace_array(x_start, y_end, n_points) result(res)
    ! Arguments
    real(dp), intent(in) :: x_start        ! The actual starting value (e.g., 1.0)
    real(dp), intent(in) :: y_end          ! The actual ending value (e.g., 1000.0)
    integer, intent(in)  :: n_points       ! The number of points to generate

    ! Function Result
    real(dp), allocatable :: res(:)         ! The array of log-spaced numbers

    ! Local variables
    real(dp) :: log_x, log_y, step
    integer  :: i

    ! Handle edge case: if N <= 0, return an empty array
    if (n_points <= 0) then
      allocate (res(0))
      return
    end if

    ! 1. Calculate the starting and ending exponents (log10(x) and log10(y))
    log_x = log10(x_start)
    log_y = log10(y_end)

    ! Handle edge case: if N = 1, return an array with just the end point
    if (n_points == 1) then
      allocate (res(1))
      res(1) = y_end
      return
    end if

    ! 2. Calculate the linear step size for the exponents
    step = (log_y - log_x) / real(n_points - 1, kind=dp)

    ! 3. Allocate the array and calculate the values using an array constructor
    allocate (res(n_points))
    
    ! Calculate 10**exponent, where exponent is linearly spaced
    ! The exponent is log_x + (i - 1) * step for i = 1 to n_points
    res = 10.0_dp ** ( &
           log_x + real((/ (i, i=1, n_points) /) - 1, kind=dp) * step &
          )

  end function logspace_array


SUBROUTINE ApplyEOS(T, D, Ye, Ym, Mumu, Mue, Mun, Mup, Xn, Xp, Un, Up, Mn_eff, Mp_eff)

  REAL(DP), INTENT(IN)  :: T, D, Ye, Ym
  REAL(DP), INTENT(OUT) :: Mumu, Mue, Mun, Mup, Xn, Xp, Un, Up, Mn_eff, Mp_eff

  TYPE(LeptonGasType) :: MuonGasState
  TYPE(LeptonGasType) :: ElectronGasState
  REAL(DP) :: Yp
  INTEGER  :: iDV

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

#ifdef EOSMODE_3D
  Mumu = 0.0d0
  iDV = EOSTable % DV % Indices % iElectronChemicalPotential
  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        EOSTable % DV % Offsets(iDV), &
        EOSTable % DV % Variables(iDV) % Values(:,:,:), Mue )
#else

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

#endif

END SUBROUTINE ApplyEOS


END PROGRAM wlCreate4DAbEmOpacityTable


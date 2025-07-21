PROGRAM wlTestSpeed

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, &
   Vud, massA, massV, gamma_p, gamma_n
  USE HDF5
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  USE wlLeptonEOSModule
  USE wlMuonEOS
  USE wlElectronPhotonEOS
  USE wlHelmMuonIOModuleHDF
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D
  USE wlSemiLeptonicOpacityModule4D, ONLY: &
    Opacity_CC_4D
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

  IMPLICIT NONE

  INTEGER , PARAMETER :: NP = 60, nOp = 2, iApprox = 0 ! from 0 to 3
  INTEGER , PARAMETER :: nOutTerminal = 100
  REAL(DP), PARAMETER :: masse = me , massm = mmu
  REAL(DP), PARAMETER :: massn = mn, massp = mp

  REAL(DP) :: EnuA(NP), Error(NP)
#ifdef EOSMODE_3D
    TYPE(EquationOfStateTableType) :: EOSTable
#elif defined(EOSMODE_4D)
    TYPE(EquationOfState4DTableType) :: EOSTable
#elif defined(EOSMODE_COMPOSE)
    TYPE(EquationOfStateCompOSETableType) :: EOSTable
#endif
  TYPE(MuonTableType) :: MuonTable
  TYPE(HelmTableType) :: HelmTable

  TYPE(OpacityTableType) :: OpacityTable
  REAL(dp), DIMENSION(2) :: Offset_Em

  INTEGER :: i, k, l, iCount, nThermoPoints, ierr
  INTEGER :: iEOS_Rho, iEOS_T, iEOS_Yp
  LOGICAL, PARAMETER :: DoMuons = .true.

  REAL(DP), ALLOCATABLE :: OpaA_El(:,:,:), OpaA_Table(:,:,:)
  REAL(DP), ALLOCATABLE :: OpaA_2D(:,:,:), OpaA_4D(:,:,:)
  REAL(DP), ALLOCATABLE :: OpaA_2D_OLD(:,:,:), OpaA_4D_OLD(:,:,:)
  REAL(DP), ALLOCATABLE :: xT(:), xD(:), xYe(:), xYm(:)

  REAL(DP) :: xMumu, xMue, xMun, xMup, xXn, xXp, xUn, xUp, xMn_eff, xMp_eff, xMl, xMul

  REAL(DP) :: buffer
  REAL(DP) :: t1, t2, t_start, t_end, t_2D, t_4D, t_El, t_Tab
  REAL(DP) :: t_2D_OLD, t_4D_OLD, err

  CALL MPI_INIT( ierr )

  CALL ReadEquationOfStateTableHDF( EOSTable, "BaryonsPlusHelmPlusMuonsEOS.h5" )
  CALL ReadHelmholtzTableHDF( HelmTable, "BaryonsPlusHelmPlusMuonsEOS.h5"  )
  CALL ReadMuonTableHDF( MuonTable, "BaryonsPlusHelmPlusMuonsEOS.h5" )

  iEOS_Rho = EOSTable % TS % Indices % iRho
  iEOS_T   = EOSTable % TS % Indices % iT
  iEOS_Yp  = EOSTable % TS % Indices % iYe   

  ! Initialize neutrino energies
  DO l = 1, NP
    EnuA(l) = l * 5.d0
  END DO

  CALL CalculateHorowitzWeakMagRecoil(EnuA, NP)

  ! READ thermodynamic data
  OPEN(UNIT=123, FILE=trim(adjustl('ThermoConditions.dat')), STATUS='OLD', ACTION='READ')
  READ(123,*) nThermoPoints
  READ(123,*)

  ! You can also set nThermoPoints to a smaller value for quick checks
  nThermoPoints = 50
  ALLOCATE(OpaA_El(NP, nThermoPoints, nOp))
  ALLOCATE(OpaA_2D(NP, nThermoPoints, nOp))
  ALLOCATE(OpaA_4D(NP, nThermoPoints, nOp))
  ALLOCATE(OpaA_2D_OLD(NP, nThermoPoints, nOp))
  ALLOCATE(OpaA_4D_OLD(NP, nThermoPoints, nOp))
  ALLOCATE(OpaA_Table(NP, nThermoPoints, nOp))
  ALLOCATE(xT(nThermoPoints) , xD(nThermoPoints), &
           xYe(nThermoPoints), xYm(nThermoPoints))

  DO i = 1, nThermoPoints
    READ(123,*) xT(i), xD(i), xYe(i), xYm(i), &
      buffer, buffer, buffer, buffer, buffer, buffer
  END DO
  CLOSE(123)

    
  ! Read Opacity Table
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable,   &
       FileName_EmAb_Option                 &
       = "wl-Op-SFHo-15-25-50-E40-EmAb.h5", &
       EquationOfStateTableName_Option      &
       = "BaryonsPlusHelmPlusMuonsEOS.h5",  &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )
  WRITE(*,*) 'Done Reading Tables'

  Offset_Em = OpacityTable % EmAb % Offsets

  ASSOCIATE( TableEm1  => OpacityTable % EmAb % Opacity(1) % Values, &
             TableEm2  => OpacityTable % EmAb % Opacity(2) % Values, &
             iRho_Op   => OpacityTable % TS % Indices % iRho, &
             iT_Op     => OpacityTable % TS % Indices % iT, &
             iYe_Op    => OpacityTable % TS % Indices % iYe )

  CALL CPU_TIME(t_start)

  OpaA_El = 0.d0
  OpaA_2D = 0.d0
  OpaA_2D_OLD = 0.d0
  OpaA_4D = 0.d0
  OpaA_4D_OLD = 0.d0

  t_El  = 0.0d0
  t_Tab = 0.0d0
  t_2D  = 0.0d0
  t_4D  = 0.0d0
  t_2D_OLD = 0.0d0
  t_4D_OLD = 0.0d0
  iCount = 0
  DO i = 1, nThermoPoints

    CALL ApplyEOS(xT(i), xD(i), xYe(i), xYm(i), &
          xMumu, xMue, xMun, xMup, xXn, xXp, xUn, xUp, xMn_eff, xMp_eff)

    IF (DoMuons) THEN
      xMul = xMumu
      xMl = massm
    ELSE
      xMul = xMue
      xMl = masse
    END IF

    DO k = 1, nOp

      CALL CPU_TIME(t1)
      CALL Opacity_CC_2D_GSI(iApprox, k, NP, EnuA, OpaA_2D_OLD(:, i, k), &
              xT(i), xMul, xMun, xMup, xMl, massn, massp, xUn, xUp)
      CALL CPU_TIME(t2)
      OpaA_2D_OLD(:, i, k) = OpaA_2D_OLD(:, i, k) * 1.0d5 ! 1/km to 1/cm
      t_2D_OLD = t_2D_OLD + t2 - t1

      CALL CPU_TIME(t1)
      CALL Opacity_CC_4D_GSI(iApprox, k, NP, EnuA, OpaA_4D_OLD(:, i, k), &
              xT(i), xMul, xMun, xMup, xMl, massn, massp, xUn, xUp)
      CALL CPU_TIME(t2)
      OpaA_4D_OLD(:, i, k) = OpaA_4D_OLD(:, i, k) * 1.0d5 ! 1/km to 1/cm
      t_4D_OLD = t_4D_OLD + t2 - t1
    
    END DO

    DO l=1, NP

      CALL CPU_TIME(t1)
      ! interpolate electron neutrino EmAb opacity
      CALL LogInterpolateSingleVariable & 
            ( LOG10( EnuA(l) ), LOG10( xD(i) ),  &
              LOG10( xT(i) ), xYe(i),            & 
              LOG10( OpacityTable % EnergyGrid % Values ),        &
              LOG10( OpacityTable % TS % States(iRho_Op) % Values ), &
              LOG10( OpacityTable % TS % States(iT_Op) % Values ),   &
              OpacityTable % TS % States(iYe_Op) % Values,           &
              Offset_Em(iNu_e), TableEm1, OpaA_Table(l, i, 1) )

      ! interpolate electron antineutrino EmAb opacity
      CALL LogInterpolateSingleVariable & 
            ( LOG10( EnuA(l) ), LOG10( xD(i) ),  &
              LOG10( xT(i) ), xYe(i),            & 
              LOG10( OpacityTable % EnergyGrid % Values ),        &
              LOG10( OpacityTable % TS % States(iRho_Op) % Values ), &
              LOG10( OpacityTable % TS % States(iT_Op) % Values ),   &
              OpacityTable % TS % States(iYe_Op) % Values,           &
              Offset_Em(iNu_e_bar), TableEm2, OpaA_Table(l, i, 2) )
      CALL CPU_TIME(t2)
      t_Tab = t_Tab + t2 - t1

      CALL CPU_TIME(t1)
      IF (DoMuons) THEN
        CALL ElasticAbsorptionOpacityNum (xD(i), xT(i), EnuA(l), l, &
          xMun, xMn_eff, xUn, xXn, xMup, xMp_eff, xUp, xXp, xMumu , &
          OpaA_El(l, i, 1), OpaA_El(l, i, 2) , IncludeWeakMagRecoil_Option=.true.)
      ELSE
        CALL ElasticAbsorptionOpacityNue (xD(i), xT(i), EnuA(l), l, &
          xMun, xMn_eff, xUn, xXn, xMup, xMp_eff, xUp, xXp, xMue , &
          OpaA_El(l, i, 1), OpaA_El(l, i, 2), IncludeWeakMagRecoil_Option=.true.)
      END IF
      CALL CPU_TIME(t2)
      t_El = t_El + t2 - t1

      DO k = 1, nOp
        CALL CPU_TIME(t1)
        CALL Opacity_CC_2D(iApprox, k, EnuA(l), OpaA_2D(l, i, k), &
              xT(i), xMul, xMun, xMup, xMl, massn, massp, xUn, xUp)
        CALL CPU_TIME(t2)
        t_2D = t_2D + t2 - t1

        CALL CPU_TIME(t1)
        CALL Opacity_CC_4D(iApprox, k, EnuA(l), OpaA_4D(l, i, k), &
              xT(i), xMul, xMun, xMup, xMl, massn, massp, xUn, xUp)
        CALL CPU_TIME(t2)
        t_4D = t_4D + t2 - t1
      END DO

      iCount = iCount + 1
      IF (MOD(iCount, nOutTerminal) == 0) WRITE(*,*) 'Done with', iCount
    END DO
  END DO
  
  END ASSOCIATE

  CALL CPU_TIME(t_end)

  WRITE(*,'(/,A,f13.6)') 'Total wall‑clock time :', t_end - t_start
  WRITE(*,'(A,f13.6)')   ' Elastic Opacity      :', t_El
  WRITE(*,'(A,f13.6)')   ' Table   Opacity      :', t_Tab
  WRITE(*,'(A,f13.6)')   ' OLD 2D kernels       :', t_2D_old
  WRITE(*,'(A,f13.6)')   ' NEW 2D kernels       :', t_2D
  WRITE(*,'(A,f13.6)')   ' OLD 4D kernels       :', t_4D_old
  WRITE(*,'(A,f13.6)')   ' NEW 4D kernels       :', t_4D
  
  OPEN (200, file='OpaA_El.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_2D
  CLOSE(200)

  OPEN (201, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  WRITE(201) OpaA_2D
  CLOSE(201)

  OPEN (202, file='OpaA_2D_OLD.bin', status='replace', form='unformatted', access='stream')
  WRITE(202) OpaA_2D_OLD
  CLOSE(202)

  OPEN (203, file='OpaA_4D.bin', status='replace', form='unformatted', access='stream')
  WRITE(203) OpaA_4D
  CLOSE(203)

  OPEN (204, file='OpaA_4D_OLD.bin', status='replace', form='unformatted', access='stream')
  WRITE(204) OpaA_4D_OLD
  CLOSE(204)

  DEALLOCATE(OpaA_2D, OpaA_2D_OLD, OpaA_4D, OpaA_4D_OLD, OpaA_El, OpaA_Table)
  DEALLOCATE(xT, xD, xYe, xYm)

CONTAINS

SUBROUTINE ApplyEOS(T, D, Ye, Ym, Mumu, Mue, Mun, Mup, Xn, Xp, Un, Up, Mn_eff, Mp_eff)

  REAL(DP), INTENT(IN)  :: T, D, Ye, Ym
  REAL(DP), INTENT(OUT) :: Mumu, Mue, Mun, Mup, Xn, Xp, Un, Up, Mn_eff, Mp_eff

  TYPE(MuonStateType) :: MuonState
  TYPE(ElectronPhotonStateType) :: ElectronPhotonState
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

  ! Muons
  MuonState % t     = T
  MuonState % rhoym = D * Ym
  CALL FullMuonEOS(MuonTable, MuonState)
  Mumu = MuonState % mu

  ! Electrons
  ElectronPhotonState % t   = T
  ElectronPhotonState % rho = D
  ElectronPhotonState % ye  = Ye
  CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)
  Mue = ElectronPhotonState % mue

END SUBROUTINE ApplyEOS

END PROGRAM wlTestSpeed



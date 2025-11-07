PROGRAM wlCompareOpacities

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, kmev, &
   Vud, massA, massV, gamma_p, gamma_n, kmev_inv
  USE HDF5
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

  INTEGER , PARAMETER :: NP = 40, nOp = 2, nApprox = 4, nE_2D = 50
  INTEGER , PARAMETER :: nOutTerminal = 10000
  LOGICAL , PARAMETER :: UseReddy = .TRUE.
  LOGICAL , PARAMETER :: IncldudeGSI_2D = .TRUE.
  LOGICAL , PARAMETER :: IncldudeGSI_4D = .FALSE.
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
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(HelmTableType) :: HelmTableMuons

  TYPE(OpacityTableType) :: OpacityTable
  REAL(dp), DIMENSION(2) :: Offset_Em

  INTEGER :: i, k, l, j, iCount, nThermoPoints, ierr
  INTEGER :: iEOS_Rho, iEOS_T, iEOS_Yp

  REAL(DP) :: WeakMagCorrLep(NP), WeakMagCorrLepBar(NP)
  REAL(DP) :: WeakMagCorrLep_WL(NP), WeakMagCorrLepBar_WL(NP)
  REAL(DP), ALLOCATABLE :: OpaA_El(:,:,:,:,:,:), OpaA_El_bare(:,:,:,:,:,:), OpaA_Table(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: OpaA_2D(:,:,:,:,:,:), OpaA_4D(:,:,:,:,:,:)

  REAL(DP) :: xT, xD, xYe, xMumu, xMue, xMun, xMup, xXn, xXp, xUn, xUp, xMn_eff, xMp_eff, xMl, xMul

  REAL(DP) :: buffer
  REAL(DP) :: t1, t2, t_start, t_end, t_2D_OLD, t_2D, t_4D, t_El, t_Tab
  REAL(DP) :: err
  CHARACTER(len=128) :: FileName_EmAb, FileName_EOS

  INTEGER :: nPoints(3), iRho, iT, iYe, iYe_table
  INTEGER :: iRho_min, iT_max, iYe_array(3)
  INTEGER :: nRho, nT, nYe

  iYe_array = (/ 1, 5, 10 /)
  iRho_min = 161
  iT_max = 60

  CALL MPI_INIT( ierr )

  CALL InitializeHDF( )
#ifdef EOSMODE_3D
  FileName_EOS = "wl-EOS-SFHo-15-25-50.h5"
  CALL ReadEquationOfStateTableHDF( EOSTable, FileName_EOS )
#else
  FileName_EOS = "BaryonsPlusPhotonsPlusLeptonsEOS.h5"
  CALL ReadEquationOfStateTableHDF( EOSTable, FileName_EOS )
  CALL ReadHelmholtzTableHDF( HelmTableElectrons, FileName_EOS, "HelmTableElectrons"  )
#endif 

  ! Read Opacity Table
  IF (UseReddy) THEN
    FileName_EmAb = "wl-Op-SFHo-15-25-50-E40-EmAb-Reddy.h5"
  ELSE
    FileName_EmAb = "wl-Op-SFHo-15-25-50-E40-EmAb.h5"
  ENDIF
  CALL ReadOpacityTableHDF( OpacityTable,   &
       FileName_EmAb_Option                 &
       = FileName_EmAb, &
       EquationOfStateTableName_Option      &
       = FileName_EOS,  &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )
  WRITE(*,*) 'Done Reading Tables'

  iEOS_Rho = EOSTable % TS % Indices % iRho
  iEOS_T   = EOSTable % TS % Indices % iT
  iEOS_Yp  = EOSTable % TS % Indices % iYe 
  nPoints(:)  = EOSTable % TS % nPoints

  ! Initialize neutrino energies NP must be equal to OpacityTable % EnergyGrid % nPointsE
  DO l = 1, NP
    EnuA(l) = OpacityTable % EnergyGrid % Values(l)
    CALL CalculateHorowitzWeakMagRecoil(EnuA(l), WeakMagCorrLep(l), WeakMagCorrLepBar(l))
  END DO

  CALL cc_weak_mag_weaklib( EnuA, WeakMagCorrLep_WL, WeakMagCorrLepBar_WL, NP )

  DO l = 1, NP
    WRITE(*,'(3F6.2)') EnuA(l), WeakMagCorrLep(l)/WeakMagCorrLep_WL(l), WeakMagCorrLepBar(l)/WeakMagCorrLepBar_WL(l)
  END DO
  
  WeakMagCorrLep = WeakMagCorrLep_WL
  WeakMagCorrLepBar = WeakMagCorrLepBar_WL

  nRho = nPoints(1) - iRho_min
  nT = iT_max
  nYe = 3

  ! You can also set nThermoPoints to a smaller value for quick checks
  ! nThermoPoints = 100
  ALLOCATE(OpaA_El(NP, nRho, nT, nYe, 2, nOp))
  ALLOCATE(OpaA_El_bare(NP, nRho, nT, nYe, 2, nOp))
  ALLOCATE(OpaA_2D(NP, nRho, nT, nYe, nApprox, nOp))
  ALLOCATE(OpaA_4D(NP, nRho, nT, nYe, nApprox, nOp))
  ALLOCATE(OpaA_Table(NP, nRho, nT, nYe, nOp))
  WRITE(*,*) NP, nRho, nT, nYe, 2, nOp

  Offset_Em = OpacityTable % EmAb % Offsets

  ASSOCIATE( TableEm1  => OpacityTable % EmAb % Opacity(1) % Values, &
             TableEm2  => OpacityTable % EmAb % Opacity(2) % Values, &
             iRho_Op   => OpacityTable % TS % Indices % iRho, &
             iT_Op     => OpacityTable % TS % Indices % iT, &
             iYe_Op    => OpacityTable % TS % Indices % iYe )

  CALL CPU_TIME(t_start)

  OpaA_El = 0.d0
  OpaA_El_bare = 0.d0
  OpaA_2D = 0.d0

  t_El  = 0.0d0
  t_Tab = 0.0d0
  t_2D  = 0.0d0
  t_4D  = 0.0d0
  iCount = 0
#ifdef USE_OMP
  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(iYe_table, xD, xT, xYe, xMumu, xMue, xMun, xMup, &
  !$OMP xXn, xXp, xUn, xUp, xMn_eff, xMp_eff, xMul, xMl, l, k, j, t1, t2), &
  !$OMP REDUCTION(+:t_Tab, t_El, t_2D, t_4D, iCount)
#endif
  DO iRho = 1, nRho
    DO iT = 1, nT
      DO iYe = 1, 3

        iYe_table = iYe_array(iYe)

        xD = EOSTable % TS % States(1) % Values(iRho_min + iRho)
        xT = EOSTable % TS % States(2) % Values(iT)
        xYe = EOSTable % TS % States(3) % Values(iYe_table)

        CALL ApplyEOS(xT, xD, xYe, 0.0d0, &
              xMumu, xMue, xMun, xMup, xXn, xXp, xUn, xUp, xMn_eff, xMp_eff)

        xMul = xMue
        xMl = masse

        DO l=1, NP

          CALL CPU_TIME(t1)
          ! interpolate electron neutrino EmAb opacity
          CALL LogInterpolateSingleVariable & 
                ( LOG10( EnuA(l) ), LOG10( xD ),  &
                  LOG10( xT ), xYe,            & 
                  LOG10( OpacityTable % EnergyGrid % Values ),        &
                  LOG10( OpacityTable % TS % States(iRho_Op) % Values ), &
                  LOG10( OpacityTable % TS % States(iT_Op) % Values ),   &
                  OpacityTable % TS % States(iYe_Op) % Values,           &
                  Offset_Em(iNu_e), TableEm1, OpaA_Table(l, iRho, iT, iYe, 1) )

          ! interpolate electron antineutrino EmAb opacity
          CALL LogInterpolateSingleVariable & 
                ( LOG10( EnuA(l) ), LOG10( xD ),  &
                  LOG10( xT ), xYe,            & 
                  LOG10( OpacityTable % EnergyGrid % Values ),        &
                  LOG10( OpacityTable % TS % States(iRho_Op) % Values ), &
                  LOG10( OpacityTable % TS % States(iT_Op) % Values ),   &
                  OpacityTable % TS % States(iYe_Op) % Values,           &
                  Offset_Em(iNu_e_bar), TableEm2, OpaA_Table(l, iRho, iT, iYe, 2) )
          CALL CPU_TIME(t2)
          t_Tab = t_Tab + t2 - t1

          CALL CPU_TIME(t1)
          CALL ElasticAbsorptionOpacityNue ( xD, xT, EnuA(l), &
            xMun, xMn_eff, xUn, xXn, xMup, xMp_eff, xUp, xXp, xMue, &
            OpaA_El(l, iRho, iT, iYe, 1, 1), OpaA_El(l, iRho, iT, iYe, 1, 2) )
          
          CALL ElasticAbsorptionOpacityNue ( xD, xT, EnuA(l), &
            xMun, mn, 0.0d0, 1.0d0 - xYe, xMup, mp, 0.0d0, xYe, xMue, &
            OpaA_El_bare(l, iRho, iT, iYe, 1, 1), OpaA_El_bare(l, iRho, iT, iYe, 1, 2) )

          OpaA_El(l, iRho, iT, iYe, 2, 1) = OpaA_El(l, iRho, iT, iYe, 1, 1) * WeakMagCorrLep(l)
          OpaA_El(l, iRho, iT, iYe, 2, 2) = OpaA_El(l, iRho, iT, iYe, 1, 2) * WeakMagCorrLepBar(l)
          
          OpaA_El_bare(l, iRho, iT, iYe, 2, 1) = OpaA_El_bare(l, iRho, iT, iYe, 1, 1) * WeakMagCorrLep(l)
          OpaA_El_bare(l, iRho, iT, iYe, 2, 2) = OpaA_El_bare(l, iRho, iT, iYe, 1, 2) * WeakMagCorrLepBar(l)
          
          CALL CPU_TIME(t2)
          t_El = t_El + t2 - t1

          IF (IncldudeGSI_2D) THEN
            DO k = 1, nOp
              DO j = 1, nApprox
                CALL CPU_TIME(t1)
                CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, iRho, iT, iYe, j, k), &
                      xT*kmev, xMul, xMun, xMup, xMl, xMn_eff, xMp_eff, xUn, xUp, nE_2D)
                CALL CPU_TIME(t2)
                t_2D = t_2D + t2 - t1
              END DO
            END DO
          ENDIF

          IF (IncldudeGSI_4D) THEN
            DO k = 1, nOp
              DO j = 1, nApprox
                CALL CPU_TIME(t1)
                CALL Opacity_CC_4D(j-1, k, EnuA(l), OpaA_4D(l, iRho, iT, iYe, j, k), &
                      xT*kmev, xMul, xMun, xMup, xMl, xMn_eff, xMp_eff, xUn, xUp)
                CALL CPU_TIME(t2)
                t_4D = t_4D + t2 - t1
              END DO
            END DO
          ENDIF

          iCount = iCount + 1
          IF (MOD(iCount, nOutTerminal) == 0) THEN
            WRITE(*,*) 'Done with', iCount
          ENDIF
        END DO

      END DO
    END DO
  END DO
#ifdef USE_OMP
  !$OMP END PARALLEL DO
#endif

  END ASSOCIATE

  CALL CPU_TIME(t_end)

  WRITE(*,'(/,A,f13.6)') 'Total wall‑clock time :', t_end - t_start
  WRITE(*,'(A,f13.6)')   ' Elastic Opacity      :', t_El
  WRITE(*,'(A,f13.6)')   ' Table   Opacity      :', t_Tab
  WRITE(*,'(A,f13.6)')   ' NEW 2D kernels       :', t_2D
  WRITE(*,'(A,f13.6)')   ' NEW 4D kernels       :', t_4D
  
  OPEN (200, file='OpaCompare_El.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_El
  CLOSE(200)

  OPEN (200, file='OpaCompare_El_bare.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_El_bare
  CLOSE(200)

  IF (IncldudeGSI_2D) THEN
    OPEN (200, file='OpaCompare_2D.bin', status='replace', form='unformatted', access='stream')
    WRITE(200) OpaA_2D
    CLOSE(200)
  ENDIF

  IF (IncldudeGSI_4D) THEN
    OPEN (200, file='OpaCompare_4D.bin', status='replace', form='unformatted', access='stream')
    WRITE(200) OpaA_4D
    CLOSE(200)
  ENDIF

  IF (UseReddy) THEN
    OPEN (200, file='OpaCompare_Table_Reddy.bin', status='replace', form='unformatted', access='stream')
  ELSE
    OPEN (200, file='OpaCompare_Table.bin', status='replace', form='unformatted', access='stream')
  ENDIF
  WRITE(200) OpaA_Table
  CLOSE(200)

  DEALLOCATE(OpaA_2D)
  DEALLOCATE(OpaA_El)
  DEALLOCATE(OpaA_El_bare)
  DEALLOCATE(OpaA_Table)

CONTAINS

SUBROUTINE ApplyEOS(T, D, Ye, Ym, Mumu, Mue, Mun, Mup, Xn, Xp, Un, Up, Mn_eff, Mp_eff)

  REAL(DP), INTENT(IN)  :: T, D, Ye, Ym
  REAL(DP), INTENT(OUT) :: Mumu, Mue, Mun, Mup, Xn, Xp, Un, Up, Mn_eff, Mp_eff
  REAL(DP), PARAMETER   :: dmnp = 1.29333922d0

  TYPE(LeptonGasType) :: ElectronGasState
  TYPE(LeptonGasType) :: MuonGasState
  REAL(DP) :: Yp, min_M, OS_M_new
  INTEGER  :: iDV

  Yp = Ye + Ym

#ifdef EOSMODE_3D
  ! Shift chemical potential
  ! Apply the shift for proton chemical potential
  iDV = EOSTable % DV % Indices % iProtonChemicalPotential
  ASSOCIATE( Mp_T  => EOSTable % DV % Variables(iDV) % Values(:,:,:), &
             OS_Mp => EOSTable % DV % Offsets(iDV) )

  IF ( OS_Mp > 0.0d0 ) THEN
    min_M = -0.5d0 * OS_Mp
  ELSE
    min_M = MINVAL( 10.0d0**Mp_T )
  ENDIF
  OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + mp + dmnp )
  Mp_T     = LOG10( 10.0d0**Mp_T - OS_Mp + mp + dmnp + OS_M_new)
  OS_Mp    = OS_M_new

  CALL LogInterpolateSingleVariable   &
    ( D, T, Yp, &
      EOSTable % TS % States(iEOS_Rho) % Values, &
      EOSTable % TS % States(iEOS_T) % Values,   &
      EOSTable % TS % States(iEOS_Yp) % Values,  &
      OS_Mp, &
      Mp_T, Mup )
  END ASSOCIATE

  ! Apply the shift for neutron chemical potential
  iDV = EOSTable % DV % Indices % iNeutronChemicalPotential
  ASSOCIATE( Mn_T  => EOSTable % DV % Variables(iDV) % Values(:,:,:), &
             OS_Mn => EOSTable % DV % Offsets(iDV) )

  IF ( OS_Mn > 0.0d0 ) THEN
    min_M = -0.5d0 * OS_Mn
  ELSE
    min_M = MINVAL( 10.0d0**Mn_T )
  ENDIF
  OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + mp + dmnp )
  Mn_T     = LOG10( 10.0d0**Mn_T - OS_Mn + mp + dmnp + OS_M_new)
  OS_Mn    = OS_M_new

  CALL LogInterpolateSingleVariable   &
      ( D, T, Yp, &
        EOSTable % TS % States(iEOS_Rho) % Values, &
        EOSTable % TS % States(iEOS_T) % Values,   &
        EOSTable % TS % States(iEOS_Yp) % Values,  &
        OS_Mn, &
        Mn_T, Mun )

  END ASSOCIATE
#else
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
#endif

#ifdef EOSMODE_3D
  Mn_eff = mn
  Mp_eff = mp
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

  ! Electrons
  ElectronGasState % t   = T
  ElectronGasState % rho = D
  ElectronGasState % yL  = Ye
  CALL LeptonGasEOS(HelmTableElectrons, ElectronGasState)
  Mue = ElectronGasState % mu

  ! Electrons
  MuonGasState % t   = T
  MuonGasState % rho = D
  MuonGasState % yL  = Ye
  CALL LeptonGasEOS(HelmTableMuons, MuonGasState)
  Mumu = MuonGasState % mu

#endif

END SUBROUTINE ApplyEOS

END PROGRAM wlCompareOpacities

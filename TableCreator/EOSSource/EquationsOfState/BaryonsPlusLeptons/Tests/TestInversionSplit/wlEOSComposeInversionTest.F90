PROGRAM wlEOSComposeInversionTest

  ! Compares temperature inversions (from E, P, and S) obtained with:
  !   1) Legacy3D : legacy 3D table (baryons + electrons + photons combined
  !                 at table creation), inverted with wlEOSInversionModule
  !   2) Combined : BaryonsPlusPhotonsPlusLeptonsEOS.h5, leptons/photons added
  !                 at the table nodes, then interpolated together
  !                 (wlEOSComponentsCombinedInversionModule)
  !   3) Separate : BaryonsPlusPhotonsPlusLeptonsEOS.h5, baryons interpolated,
  !                 leptons/photons evaluated at the exact point
  !                 (wlEOSComponentsSeparateInversionModule)
  !
  ! Test A (round trip): X is computed with each method's own forward EOS and
  !   inverted with the same method -> measures the inversion itself.
  ! Test B (common target): X is computed with the Separate forward EOS and
  !   inverted with every method -> measures the temperature each method
  !   would return for the same thermodynamic state.

  USE wlKindModule, ONLY: &
    dp
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType, &
    EquationOfStateCompOSETableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_3D_Custom_Point
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log
  USE wlEOSInversionModule, ONLY: &
    InitializeEOSInversion, &
    DescribeEOSInversionError, &
    ComputeTemperatureWith_DEY_Single_Guess_Error, &
    ComputeTemperatureWith_DEY_Single_NoGuess_Error, &
    ComputeTemperatureWith_DPY_Single_Guess_Error, &
    ComputeTemperatureWith_DPY_Single_NoGuess_Error, &
    ComputeTemperatureWith_DSY_Single_Guess_Error, &
    ComputeTemperatureWith_DSY_Single_NoGuess_Error
  USE wlEOSComponentsCombinedInversionModule, ONLY: &
    InitializeInversion_Comb => InitializeEOSComponentsInversion, &
    DescribeError_Comb       => DescribeEOSComponentsInversionError, &
    DEYpYl_Guess_Comb        => ComputeTemperatureWith_DEYpYl_Single_Guess_Error, &
    DEYpYl_NoGuess_Comb      => ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error, &
    DPYpYl_Guess_Comb        => ComputeTemperatureWith_DPYpYl_Single_Guess_Error, &
    DPYpYl_NoGuess_Comb      => ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error, &
    DSYpYl_Guess_Comb        => ComputeTemperatureWith_DSYpYl_Single_Guess_Error, &
    DSYpYl_NoGuess_Comb      => ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error
  USE wlEOSComponentsSeparateInversionModule, ONLY: &
    InitializeInversion_Sep => InitializeEOSComponentsInversion, &
    DescribeError_Sep       => DescribeEOSComponentsInversionError, &
    DEYpYl_Guess_Sep        => ComputeTemperatureWith_DEYpYl_Single_Guess_Error, &
    DEYpYl_NoGuess_Sep      => ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error, &
    DPYpYl_Guess_Sep        => ComputeTemperatureWith_DPYpYl_Single_Guess_Error, &
    DPYpYl_NoGuess_Sep      => ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error, &
    DSYpYl_Guess_Sep        => ComputeTemperatureWith_DSYpYl_Single_Guess_Error, &
    DSYpYl_NoGuess_Sep      => ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error
  USE wlHelmIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: &
    GetPhotonLeptonGasEOS

  IMPLICIT NONE

  INTEGER, PARAMETER :: nPoints = 2**12

  INTEGER, PARAMETER :: iLegacy = 1, iCombined = 2, iSeparate = 3
  INTEGER, PARAMETER :: iVar_E = 1, iVar_P = 2, iVar_S = 3
  INTEGER, PARAMETER :: iGoodGuess = 1, iBadGuess = 2, iNoGuess = 3

  CHARACTER(8),  PARAMETER :: MethodNames(3) = [ 'Legacy3D', 'Combined', 'Separate' ]
  CHARACTER(1),  PARAMETER :: VarNames(3)    = [ 'E', 'P', 'S' ]
  CHARACTER(10), PARAMETER :: GuessNames(3)  = [ 'Good Guess', 'Bad Guess ', 'No Guess  ' ]

  ! --- Muon fraction of the sampled points: Ym = Ym_over_Yp * Yp.
  !     The legacy table has no muons, so keep this zero (or tiny) for
  !     a like-for-like comparison with Legacy3D.
  REAL(dp), PARAMETER :: Ym_over_Yp = 0.0d0

  CHARACTER(len=128), PARAMETER :: &
    BaryonPlusHelmTableName = 'BaryonsPlusPhotonsPlusLeptonsEOS.h5', &
    LegacyTableName         = '3DEOSTable.h5'

  TYPE(EquationOfStateCompOSETableType) :: EOSBaryonTable
  TYPE(EquationOfStateTableType)        :: EOSLegacyTable
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(HelmTableType) :: HelmTableMuons

  ! --- New (split) table: baryons only, P stored linearly, E and S as log10
  REAL(dp), ALLOCATABLE :: Ds_bary(:), Ts_bary(:), Yps_bary(:)
  REAL(dp), ALLOCATABLE :: Xs_bary(:,:,:,:)
  ! --- Baryons + leptons + photons summed at the nodes (log10), Combined forward
  REAL(dp), ALLOCATABLE :: Xs_full(:,:,:,:)
  REAL(dp) :: OS_bary(3)
  ! --- Legacy table (log10)
  REAL(dp), ALLOCATABLE :: Ds_leg(:), Ts_leg(:), Yps_leg(:)
  REAL(dp), ALLOCATABLE :: Xs_leg(:,:,:,:)
  REAL(dp) :: OS_leg(3)

  INTEGER  :: iDV_bary(3), iDV_leg(3)
  INTEGER  :: iD_bary, iT_bary, iYp_bary
  INTEGER  :: SizeDs, SizeTs, SizeYps
  INTEGER  :: iD, iT, iYp, iP, iVar, iGuess, iMethod, iRef, n_rndm
  REAL(dp) :: Ymu_temp, P_LeptPhot, E_LeptPhot, S_LeptPhot
  REAL(dp) :: Amp, minT, maxT, tCPU

  REAL(dp), DIMENSION(nPoints) :: &
    D, T, Yp, Ye, Ym, T_Good, T_Bad, T_Inv, dT, RelDiff, &
    rndm_D, rndm_T, rndm_Yp
  INTEGER,  DIMENSION(nPoints) :: Error_Inv
  REAL(dp) :: X(nPoints,3,3) ! (point, variable, method)

  CALL InitializeHDF( )

  ! --- New table: Helmholtz electrons and muons, and baryons ---
  CALL ReadHelmholtzTableHDF( HelmTableElectrons, BaryonPlusHelmTableName, "HelmTableElectrons" )
  CALL ReadHelmholtzTableHDF( HelmTableMuons, BaryonPlusHelmTableName, &
      "HelmTableMuons", eos_MinD = 1.0d7 )
  CALL ReadEquationOfStateTableHDF( EOSBaryonTable, BaryonPlusHelmTableName )

  ! --- Legacy table ---
  CALL ReadEquationOfStateTableHDF( EOSLegacyTable, LegacyTableName )

  CALL FinalizeHDF( )

  ! -------------------------------------------------------------------
  ! --- Unpack new table ----------------------------------------------
  ! -------------------------------------------------------------------

  iD_bary  = EOSBaryonTable % TS % Indices % iRho
  iT_bary  = EOSBaryonTable % TS % Indices % iT
  iYp_bary = EOSBaryonTable % TS % Indices % iYe

  SizeDs  = EOSBaryonTable % TS % nPoints(iD_bary )
  SizeTs  = EOSBaryonTable % TS % nPoints(iT_bary )
  SizeYps = EOSBaryonTable % TS % nPoints(iYp_bary)

  ALLOCATE( Ds_bary(SizeDs), Ts_bary(SizeTs), Yps_bary(SizeYps) )
  Ds_bary  = EOSBaryonTable % TS % States(iD_bary ) % Values
  Ts_bary  = EOSBaryonTable % TS % States(iT_bary ) % Values
  Yps_bary = EOSBaryonTable % TS % States(iYp_bary) % Values

  iDV_bary(iVar_E) = EOSBaryonTable % DV % Indices % iInternalEnergyDensity
  iDV_bary(iVar_P) = EOSBaryonTable % DV % Indices % iPressure
  iDV_bary(iVar_S) = EOSBaryonTable % DV % Indices % iEntropyPerBaryon

  ALLOCATE( Xs_bary(SizeDs,SizeTs,SizeYps,3) )
  ALLOCATE( Xs_full(SizeDs,SizeTs,SizeYps,3) )
  DO iVar = 1, 3
    OS_bary(iVar) = EOSBaryonTable % DV % Offsets(iDV_bary(iVar))
    Xs_bary(:,:,:,iVar) = EOSBaryonTable % DV % Variables(iDV_bary(iVar)) % Values
  END DO

  ! --- Build full EOS (baryons + leptons + photons at the nodes) ---
  DO iYp = 1, SizeYps
    DO iT = 1, SizeTs
      DO iD = 1, SizeDs

        Ymu_temp = Yps_bary(iYp) * Ym_over_Yp

        CALL GetPhotonLeptonGasEOS( &
          Ds_bary(iD), Ts_bary(iT), Yps_bary(iYp) - Ymu_temp, Ymu_temp, &
          HelmTableElectrons, HelmTableMuons, P_LeptPhot, E_LeptPhot, S_LeptPhot)

        Xs_full(iD,iT,iYp,iVar_E) = LOG10( 10.0d0**( Xs_bary(iD,iT,iYp,iVar_E) ) + E_LeptPhot )
        Xs_full(iD,iT,iYp,iVar_P) = LOG10(           Xs_bary(iD,iT,iYp,iVar_P)   + P_LeptPhot )
        Xs_full(iD,iT,iYp,iVar_S) = LOG10( 10.0d0**( Xs_bary(iD,iT,iYp,iVar_S) ) + S_LeptPhot )

        IF ( ANY( Xs_full(iD,iT,iYp,:) .NE. Xs_full(iD,iT,iYp,:) ) ) THEN
          WRITE(*,*) 'NaN in full EOS', iD, iT, iYp, Xs_bary(iD,iT,iYp,:), &
            E_LeptPhot, P_LeptPhot, S_LeptPhot
          STOP
        END IF

      END DO
    END DO
  END DO

  ! -------------------------------------------------------------------
  ! --- Unpack legacy table -------------------------------------------
  ! -------------------------------------------------------------------

  ASSOCIATE( TS => EOSLegacyTable % TS )
  ALLOCATE( Ds_leg (TS % nPoints(TS % Indices % iRho)) )
  ALLOCATE( Ts_leg (TS % nPoints(TS % Indices % iT  )) )
  ALLOCATE( Yps_leg(TS % nPoints(TS % Indices % iYe )) )
  Ds_leg  = TS % States(TS % Indices % iRho) % Values
  Ts_leg  = TS % States(TS % Indices % iT  ) % Values
  Yps_leg = TS % States(TS % Indices % iYe ) % Values
  END ASSOCIATE

  iDV_leg(iVar_E) = EOSLegacyTable % DV % Indices % iInternalEnergyDensity
  iDV_leg(iVar_P) = EOSLegacyTable % DV % Indices % iPressure
  iDV_leg(iVar_S) = EOSLegacyTable % DV % Indices % iEntropyPerBaryon

  ALLOCATE( Xs_leg(SIZE(Ds_leg),SIZE(Ts_leg),SIZE(Yps_leg),3) )
  DO iVar = 1, 3
    OS_leg(iVar) = EOSLegacyTable % DV % Offsets(iDV_leg(iVar))
    Xs_leg(:,:,:,iVar) = EOSLegacyTable % DV % Variables(iDV_leg(iVar)) % Values
  END DO

  ! -------------------------------------------------------------------
  ! --- Initialize the three inversions -------------------------------
  ! -------------------------------------------------------------------

  CALL InitializeEOSInversion &
         ( Ds_leg, Ts_leg, Yps_leg, &
           10.0d0**( Xs_leg(:,:,:,iVar_E) ) - OS_leg(iVar_E), &
           10.0d0**( Xs_leg(:,:,:,iVar_P) ) - OS_leg(iVar_P), &
           10.0d0**( Xs_leg(:,:,:,iVar_S) ) - OS_leg(iVar_S), &
           Verbose_Option = .TRUE. )

  CALL InitializeInversion_Comb &
         ( Ds_bary, Ts_bary, Yps_bary, &
           10.0d0**( Xs_full(:,:,:,iVar_E) ) - OS_bary(iVar_E), &
           10.0d0**( Xs_full(:,:,:,iVar_P) ) - OS_bary(iVar_P), &
           10.0d0**( Xs_full(:,:,:,iVar_S) ) - OS_bary(iVar_S), &
           HelmTableElectrons, HelmTableMuons, &
           Verbose_Option = .TRUE. )

  CALL InitializeInversion_Sep &
         ( Ds_bary, Ts_bary, Yps_bary, &
           10.0d0**( Xs_full(:,:,:,iVar_E) ) - OS_bary(iVar_E), &
           10.0d0**( Xs_full(:,:,:,iVar_P) ) - OS_bary(iVar_P), &
           10.0d0**( Xs_full(:,:,:,iVar_S) ) - OS_bary(iVar_S), &
           HelmTableElectrons, HelmTableMuons, &
           Verbose_Option = .TRUE. )

  ! -------------------------------------------------------------------
  ! --- Random sampling of the table ----------------------------------
  ! -------------------------------------------------------------------

  WRITE(*,*)
  WRITE(*,'(A4,A10,I10.10)') '', 'nPoints = ', nPoints

  n_rndm = nPoints

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )
  ASSOCIATE &
    ( minD => EOSBaryonTable % TS % MinValues(iD_bary), &
      maxD => EOSBaryonTable % TS % MaxValues(iD_bary) )
  D(:) = 10.0d0**( LOG10(minD) + ( LOG10(1.0d12) - LOG10(minD) ) * rndm_D )
  PRINT*, "Min/Max D  = ", MINVAL( D ), MAXVAL( D )
  PRINT*, "             ", minD, maxD
  END ASSOCIATE

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )
  minT = EOSBaryonTable % TS % MinValues(iT_bary)
  maxT = EOSBaryonTable % TS % MaxValues(iT_bary)
  T(:) = 10.0d0**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )
  PRINT*, "Min/Max T  = ", MINVAL( T ), MAXVAL( T )
  PRINT*, "             ", minT, maxT

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Yp )
  ASSOCIATE &
    ( minYp => EOSBaryonTable % TS % MinValues(iYp_bary), &
      maxYp => EOSBaryonTable % TS % MaxValues(iYp_bary) )
  Yp(:) = minYp + ( maxYp - minYp ) * rndm_Yp
  PRINT*, "Min/Max Yp = ", MINVAL( Yp ), MAXVAL( Yp )
  PRINT*, "             ", minYp, maxYp
  END ASSOCIATE

  Ym = Yp * Ym_over_Yp
  Ye = Yp - Ym

  ! -------------------------------------------------------------------
  ! --- Forward EOS with each method ----------------------------------
  ! -------------------------------------------------------------------

  DO iMethod = 1, 3
    DO iVar = 1, 3
      DO iP = 1, nPoints
        X(iP,iVar,iMethod) = ForwardEOS( iMethod, iVar, D(iP), T(iP), Yp(iP), Ye(iP), Ym(iP) )
      END DO
    END DO
  END DO

  DO iRef = iLegacy, iSeparate, iSeparate - iLegacy
    WRITE(*,*)
    WRITE(*,'(A,A,A)') ' Forward EOS differences w.r.t. ', TRIM( MethodNames(iRef) ), &
      ':  max / mean |X - X_Ref| / |X_Ref|'
    DO iVar = 1, 3
      DO iMethod = 1, 3
        IF ( iMethod == iRef ) CYCLE
        RelDiff = ABS( X(:,iVar,iMethod) - X(:,iVar,iRef) ) / ABS( X(:,iVar,iRef) )
        WRITE(*,'(A4,A1,A3,A8,2ES12.3)') '', VarNames(iVar), '', MethodNames(iMethod), &
          MAXVAL( RelDiff ), SUM( RelDiff ) / DBLE( nPoints )
      END DO
    END DO
  END DO

  ! -------------------------------------------------------------------
  ! --- Initial guesses -----------------------------------------------
  ! -------------------------------------------------------------------

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  Amp = 0.001_dp
  dT  = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T
  T_Good = MIN( MAX( T + dT, 1.0001_dp * minT ), 0.9999_dp * maxT )

  Amp = 0.1_dp * ( 0.9999_dp * maxT - 1.0001_dp * minT )
  dT  = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T
  T_Bad = MIN( MAX( T + dT, 1.0001_dp * minT ), 0.9999_dp * maxT )

  ! -------------------------------------------------------------------
  ! --- Test A: round trip with each method's own forward EOS ---------
  ! -------------------------------------------------------------------

  WRITE(*,*)
  WRITE(*,'(A)') ' Test A: round trip (X from the method''s own forward EOS)'
  CALL WriteHeader( )
  DO iVar = 1, 3
    DO iGuess = 1, 3
      DO iMethod = 1, 3
        CALL InvertAllPoints( iMethod, iVar, iGuess, X(:,iVar,iMethod), T_Inv, Error_Inv, tCPU )
        CALL ReportInversion( iMethod, iVar, iGuess, X(:,iVar,iMethod), T_Inv, Error_Inv, tCPU )
      END DO
    END DO
    WRITE(*,*)
  END DO

  ! -------------------------------------------------------------------
  ! --- Test B: common target X from the Separate forward EOS ---------
  ! -------------------------------------------------------------------

  WRITE(*,'(A)') ' Test B: common target (X from the Separate forward EOS)'
  CALL WriteHeader( )
  DO iVar = 1, 3
    DO iGuess = 1, 3
      DO iMethod = 1, 3
        CALL InvertAllPoints( iMethod, iVar, iGuess, X(:,iVar,iSeparate), T_Inv, Error_Inv, tCPU )
        CALL ReportInversion( iMethod, iVar, iGuess, X(:,iVar,iSeparate), T_Inv, Error_Inv, tCPU )
      END DO
    END DO
    WRITE(*,*)
  END DO

  DEALLOCATE( Ds_bary, Ts_bary, Yps_bary, Xs_bary, Xs_full )
  DEALLOCATE( Ds_leg, Ts_leg, Yps_leg, Xs_leg )

CONTAINS


  ! Forward EOS: X(D,T,Yp) for the given method and variable
  REAL(dp) FUNCTION ForwardEOS( iMethod, iVar, D_P, T_P, Yp_P, Ye_P, Ym_P )

    INTEGER,  INTENT(in) :: iMethod, iVar
    REAL(dp), INTENT(in) :: D_P, T_P, Yp_P, Ye_P, Ym_P

    INTEGER  :: iD_P, iT_P, iYp_P
    REAL(dp) :: X_P, dD_P, dT_P, dYp_P
    REAL(dp) :: P_LP, E_LP, S_LP

    SELECT CASE ( iMethod )

    CASE ( iLegacy )

      ! Legacy table has no muons: all positive charge is balanced by electrons
      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, Ds_leg, Ts_leg, Yps_leg, OS_leg(iVar), &
               Xs_leg(:,:,:,iVar), X_P )

    CASE ( iCombined )

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, Ds_bary, Ts_bary, Yps_bary, OS_bary(iVar), &
               Xs_full(:,:,:,iVar), X_P )

    CASE DEFAULT

      IF ( iVar == iVar_P ) THEN

        iD_P  = Index1D_Log( D_P , Ds_bary  )
        iT_P  = Index1D_Log( T_P , Ts_bary  )
        iYp_P = Index1D_Lin( Yp_P, Yps_bary )

        ! Baryon P: linear in (log D, log T, Yp) if it is <= 0 at some corner
        ! of the cell, log-trilinear otherwise (as in wlSoundSpeedModule)
        ASSOCIATE( p => Xs_bary(iD_P:iD_P+1,iT_P:iT_P+1,iYp_P:iYp_P+1,iVar_P) )
        IF ( MINVAL( p ) <= 0.0_dp ) THEN
          dD_P  = LOG10( D_P / Ds_bary(iD_P) ) / LOG10( Ds_bary(iD_P+1) / Ds_bary(iD_P) )
          dT_P  = LOG10( T_P / Ts_bary(iT_P) ) / LOG10( Ts_bary(iT_P+1) / Ts_bary(iT_P) )
          dYp_P = ( Yp_P - Yps_bary(iYp_P) ) / ( Yps_bary(iYp_P+1) - Yps_bary(iYp_P) )
          X_P &
            = ( 1.0_dp - dYp_P ) * ( ( 1.0_dp - dT_P ) * ( ( 1.0_dp - dD_P ) * p(1,1,1) + dD_P * p(2,1,1) )   &
                                 +            dT_P   * ( ( 1.0_dp - dD_P ) * p(1,2,1) + dD_P * p(2,2,1) ) ) &
            +            dYp_P   * ( ( 1.0_dp - dT_P ) * ( ( 1.0_dp - dD_P ) * p(1,1,2) + dD_P * p(2,1,2) )   &
                                 +            dT_P   * ( ( 1.0_dp - dD_P ) * p(1,2,2) + dD_P * p(2,2,2) ) )
        ELSE
          CALL LogInterpolateSingleVariable_3D_Custom_Point &
                 ( D_P, T_P, Yp_P, &
                   Ds_bary(iD_P:iD_P+1), Ts_bary(iT_P:iT_P+1), Yps_bary(iYp_P:iYp_P+1), &
                   0.0_dp, LOG10( p ), X_P )
        END IF
        END ASSOCIATE

      ELSE

        CALL LogInterpolateSingleVariable_3D_Custom_Point &
               ( D_P, T_P, Yp_P, Ds_bary, Ts_bary, Yps_bary, OS_bary(iVar), &
                 Xs_bary(:,:,:,iVar), X_P )

      END IF

      CALL GetPhotonLeptonGasEOS( &
        D_P, T_P, Ye_P, Ym_P, &
        HelmTableElectrons, HelmTableMuons, P_LP, E_LP, S_LP )

      SELECT CASE ( iVar )
      CASE ( iVar_E )
        X_P = X_P + E_LP
      CASE ( iVar_P )
        X_P = X_P + P_LP
      CASE DEFAULT
        X_P = X_P + S_LP
      END SELECT

    END SELECT

    ForwardEOS = X_P

  END FUNCTION ForwardEOS


  ! Invert all sampled points with the given method, variable, and guess mode
  SUBROUTINE InvertAllPoints( iMethod, iVar, iGuess, X_Target, T_Out, Error_Out, tCPU )

    INTEGER,  INTENT(in)  :: iMethod, iVar, iGuess
    REAL(dp), INTENT(in)  :: X_Target(:)
    REAL(dp), INTENT(out) :: T_Out(:)
    INTEGER,  INTENT(out) :: Error_Out(:)
    REAL(dp), INTENT(out) :: tCPU

    INTEGER  :: i
    REAL(dp) :: T_Guess, tBegin, tEnd

    T_Out     = 0.0_dp
    Error_Out = 0

    CALL CPU_TIME( tBegin )

    DO i = 1, nPoints

      SELECT CASE ( iGuess )
      CASE ( iGoodGuess )
        T_Guess = T_Good(i)
      CASE ( iBadGuess )
        T_Guess = T_Bad(i)
      CASE DEFAULT
        T_Guess = 0.0_dp
      END SELECT

      CALL InvertOnePoint &
             ( iMethod, iVar, iGuess /= iNoGuess, D(i), X_Target(i), &
               Yp(i), Ye(i), Ym(i), T_Guess, T_Out(i), Error_Out(i) )

    END DO

    CALL CPU_TIME( tEnd )
    tCPU = tEnd - tBegin

  END SUBROUTINE InvertAllPoints


  SUBROUTINE InvertOnePoint &
    ( iMethod, iVar, UseGuess, D_P, X_P, Yp_P, Ye_P, Ym_P, T_Guess, T_P, Error_P )

    INTEGER,  INTENT(in)  :: iMethod, iVar
    LOGICAL,  INTENT(in)  :: UseGuess
    REAL(dp), INTENT(in)  :: D_P, X_P, Yp_P, Ye_P, Ym_P, T_Guess
    REAL(dp), INTENT(out) :: T_P
    INTEGER,  INTENT(out) :: Error_P

    SELECT CASE ( iMethod )

    CASE ( iLegacy )

      ASSOCIATE( Xs => Xs_leg(:,:,:,iVar), OS => OS_leg(iVar) )
      SELECT CASE ( iVar )
      CASE ( iVar_E )
        IF ( UseGuess ) THEN
          CALL ComputeTemperatureWith_DEY_Single_Guess_Error &
                 ( D_P, X_P, Yp_P, Ds_leg, Ts_leg, Yps_leg, Xs, OS, T_P, T_Guess, Error_P )
        ELSE
          CALL ComputeTemperatureWith_DEY_Single_NoGuess_Error &
                 ( D_P, X_P, Yp_P, Ds_leg, Ts_leg, Yps_leg, Xs, OS, T_P, Error_P )
        END IF
      CASE ( iVar_P )
        IF ( UseGuess ) THEN
          CALL ComputeTemperatureWith_DPY_Single_Guess_Error &
                 ( D_P, X_P, Yp_P, Ds_leg, Ts_leg, Yps_leg, Xs, OS, T_P, T_Guess, Error_P )
        ELSE
          CALL ComputeTemperatureWith_DPY_Single_NoGuess_Error &
                 ( D_P, X_P, Yp_P, Ds_leg, Ts_leg, Yps_leg, Xs, OS, T_P, Error_P )
        END IF
      CASE DEFAULT
        IF ( UseGuess ) THEN
          CALL ComputeTemperatureWith_DSY_Single_Guess_Error &
                 ( D_P, X_P, Yp_P, Ds_leg, Ts_leg, Yps_leg, Xs, OS, T_P, T_Guess, Error_P )
        ELSE
          CALL ComputeTemperatureWith_DSY_Single_NoGuess_Error &
                 ( D_P, X_P, Yp_P, Ds_leg, Ts_leg, Yps_leg, Xs, OS, T_P, Error_P )
        END IF
      END SELECT
      END ASSOCIATE

    CASE ( iCombined )

      ASSOCIATE( Xs => Xs_bary(:,:,:,iVar), OS => OS_bary(iVar) )
      SELECT CASE ( iVar )
      CASE ( iVar_E )
        IF ( UseGuess ) THEN
          CALL DEYpYl_Guess_Comb &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, T_Guess, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        ELSE
          CALL DEYpYl_NoGuess_Comb &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        END IF
      CASE ( iVar_P )
        IF ( UseGuess ) THEN
          CALL DPYpYl_Guess_Comb &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, T_Guess, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        ELSE
          CALL DPYpYl_NoGuess_Comb &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        END IF
      CASE DEFAULT
        IF ( UseGuess ) THEN
          CALL DSYpYl_Guess_Comb &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, T_Guess, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        ELSE
          CALL DSYpYl_NoGuess_Comb &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        END IF
      END SELECT
      END ASSOCIATE

    CASE DEFAULT

      ASSOCIATE( Xs => Xs_bary(:,:,:,iVar), OS => OS_bary(iVar) )
      SELECT CASE ( iVar )
      CASE ( iVar_E )
        IF ( UseGuess ) THEN
          CALL DEYpYl_Guess_Sep &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, T_Guess, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        ELSE
          CALL DEYpYl_NoGuess_Sep &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        END IF
      CASE ( iVar_P )
        IF ( UseGuess ) THEN
          CALL DPYpYl_Guess_Sep &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, T_Guess, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        ELSE
          CALL DPYpYl_NoGuess_Sep &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        END IF
      CASE DEFAULT
        IF ( UseGuess ) THEN
          CALL DSYpYl_Guess_Sep &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, T_Guess, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        ELSE
          CALL DSYpYl_NoGuess_Sep &
                 ( D_P, X_P, Ye_P, Ym_P, Ds_bary, Ts_bary, Yps_bary, Xs, OS, T_P, &
                   Error_P, HelmTableElectrons, HelmTableMuons )
        END IF
      END SELECT
      END ASSOCIATE

    END SELECT

  END SUBROUTINE InvertOnePoint


  SUBROUTINE WriteHeader( )

    WRITE(*,'(A4,A3,A12,A10,4A12,A8,A10)') &
      '', 'Var', 'Guess', 'Method', &
      'max dT/T', 'mean dT/T', 'max res.', 'mean res.', 'nFail', 'CPU [s]'

  END SUBROUTINE WriteHeader


  ! Relative temperature error w.r.t. the sampled T, and relative residual
  ! of the method's own forward EOS at the recovered T, over successful points
  SUBROUTINE ReportInversion( iMethod, iVar, iGuess, X_Target, T_Out, Error_Out, tCPU )

    INTEGER,  INTENT(in) :: iMethod, iVar, iGuess
    REAL(dp), INTENT(in) :: X_Target(:), T_Out(:)
    INTEGER,  INTENT(in) :: Error_Out(:)
    REAL(dp), INTENT(in) :: tCPU

    INTEGER  :: i, nFail, nOK, iFirstFail
    REAL(dp) :: ErrT(nPoints), Residual(nPoints)

    ErrT       = 0.0_dp
    Residual   = 0.0_dp
    iFirstFail = 0

    DO i = 1, nPoints
      IF ( Error_Out(i) == 0 ) THEN
        ErrT(i) = ABS( T(i) - T_Out(i) ) / T(i)
        Residual(i) &
          = ABS( ForwardEOS( iMethod, iVar, D(i), T_Out(i), Yp(i), Ye(i), Ym(i) ) &
                 - X_Target(i) ) / ABS( X_Target(i) )
      ELSE IF ( iFirstFail == 0 ) THEN
        iFirstFail = i
      END IF
    END DO

    nFail = COUNT( Error_Out /= 0 )
    nOK   = MAX( nPoints - nFail, 1 )

    WRITE(*,'(A4,A3,A12,A10,4ES12.3,I8,F10.4)') &
      '', VarNames(iVar), GuessNames(iGuess), MethodNames(iMethod), &
      MAXVAL( ErrT ), SUM( ErrT ) / DBLE( nOK ), &
      MAXVAL( Residual ), SUM( Residual ) / DBLE( nOK ), nFail, tCPU

    IF ( iFirstFail > 0 ) THEN
      WRITE(*,'(A8,A,I6,A,4ES13.5)') '', 'first failure: iP = ', iFirstFail, &
        '  D, X, T, Yp = ', D(iFirstFail), X_Target(iFirstFail), T(iFirstFail), Yp(iFirstFail)
      IF ( iMethod == iLegacy ) THEN
        CALL DescribeEOSInversionError( Error_Out(iFirstFail) )
      ELSE IF ( iMethod == iCombined ) THEN
        CALL DescribeError_Comb( Error_Out(iFirstFail) )
      ELSE
        CALL DescribeError_Sep( Error_Out(iFirstFail) )
      END IF
    END IF

  END SUBROUTINE ReportInversion


END PROGRAM wlEOSComposeInversionTest

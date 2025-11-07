PROGRAM wlEOSComposeInversionTestOnePointEnergy

  USE wlKindModule, ONLY: &
    dp
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType, &
    EquationOfStateCompOSETableType, &
    EquationOfState4DTableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_3D_Custom, &
    LogInterpolateSingleVariable_3D_Custom_Point
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log
#ifdef EOSMODE_3D
  USE wlEOSInversionModule
#else
#ifdef INTERPOLATION_SPLIT_TABLE_COMBINED
  USE wlEOSComponentsCombinedInversionModule
#else
  USE wlEOSComponentsSeparateInversionModule
#endif
#endif
  USE wlHelmIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: &
    GetPhotonLeptonGasEOS
  USE wlExtPhysicalConstantsModule, ONLY: &
    kmev, rmu, kmev_inv, ergmev, cvel
    
  IMPLICIT NONE

  INTEGER :: &
    n_rndm, &
    iP, &
    iMaxError, &
    ierr
  INTEGER :: &
    Error_E, &
    Error_P, &
    Error_S
  REAL(dp) :: &
    tBegin, &
    tEnd, &
    tCPU, &
    T_Guess, &
    Amp, &
    E_T, S_T, P_T
  REAL(dp) &
    D, T, Yp, Ye, Ym, E, P, S, T_E, T_P, T_S, dT, &
    E_int_error, P_int_error, S_int_error, &
    rndm_D, &
    rndm_T, &
    rndm_Yp, &
    Error
#ifdef EOSMODE_3D
    TYPE(EquationOfStateTableType) :: EOSTable
#elif defined(EOSMODE_4D)
    TYPE(EquationOfState4DTableType) :: EOSTable
#elif defined(EOSMODE_COMPOSE)
    TYPE(EquationOfStateCompOSETableType) :: EOSTable
#endif
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(HelmTableType) :: HelmTableMuons

  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Ds_bary, Ts_bary, Yps_bary, &
    Ds_full, Ts_full, Yps_full
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    Ps_bary, Es_bary, Ss_bary, &
    Ps_full, Es_full, Ss_full
  REAL(DP) :: &
    OS_P, OS_E, OS_S
  INTEGER :: SizeDs, SizeTs, SizeYps, &
    iD_bary, iT_bary, iYp_bary, iP_bary, iS_bary, iE_bary
  
  INTEGER  :: iD, iT, iYp
  REAL(DP) :: P_bary, E_bary, S_bary
  REAL(DP) :: P_LeptPhot, E_LeptPhot, S_LeptPhot
  REAL(DP) :: Ymu_temp, Yp_over_Ymu

  CHARACTER(len=128) :: BaryonPlusHelmTableName
  
  CALL InitializeHDF( )
  Yp_over_Ymu = 0.0d0

#ifdef EOSMODE_COMPOSE
  BaryonPlusHelmTableName = 'BaryonsPlusPhotonsPlusLeptonsEOS.h5'
  
  ! read in helmholtz table
  CALL ReadHelmholtzTableHDF( HelmTableElectrons, BaryonPlusHelmTableName, "HelmTableElectrons" )

  ! read in muon table
  CALL ReadHelmholtzTableHDF( HelmTableMuons, BaryonPlusHelmTableName, "HelmTableMuons" )

  ! read in baryon table -------------------------------
  CALL ReadEquationOfStateTableHDF( EOSTable, BaryonPlusHelmTableName )
  
#else
  BaryonPlusHelmTableName = '3DEOSTable.h5'
  ! read in baryon table -------------------------------
  CALL ReadEquationOfStateTableHDF( EOSTable, BaryonPlusHelmTableName )
#endif

  CALL FinalizeHDF( )

  iD_bary  = EOSTable % TS % Indices % iRho
  iT_bary  = EOSTable % TS % Indices % iT
  iYp_bary = EOSTable % TS % Indices % iYe

  SizeDs  = EOSTable % TS % nPoints(iD_bary )
  SizeTs  = EOSTable % TS % nPoints(iT_bary )
  SizeYps = EOSTable % TS % nPoints(iYp_bary)

  ALLOCATE( Ds_bary(SizeDs) )
  Ds_bary = EOSTable % TS % States(iD_bary) % Values

  ALLOCATE( Ts_bary(SizeTs) )
  Ts_bary = EOSTable % TS % States(iT_bary) % Values

  ALLOCATE( Yps_bary(SizeYps) )
  Yps_bary = EOSTable % TS % States(iYp_bary) % Values

  iP_bary = EOSTable % DV % Indices % iPressure
  iS_bary = EOSTable % DV % Indices % iEntropyPerBaryon
  iE_bary = EOSTable % DV % Indices % iInternalEnergyDensitY

  OS_P = EOSTable % DV % Offsets(iP_bary)
  OS_S = EOSTable % DV % Offsets(iS_bary)
  OS_E = EOSTable % DV % Offsets(iE_bary)

  WRITE(*,*) OS_P
  ALLOCATE( Ps_bary(SizeDs, SizeTs, SizeYps) )
  ALLOCATE( Ss_bary(SizeDs, SizeTs, SizeYps) )
  ALLOCATE( Es_bary(SizeDs, SizeTs, SizeYps) )

  Ps_bary = EOSTable % DV % Variables(iP_bary ) % Values
  Ss_bary = EOSTable % DV % Variables(iS_bary ) % Values
  Es_bary = EOSTable % DV % Variables(iE_bary ) % Values

  ALLOCATE( Ds_full(SizeDs) )
  ALLOCATE( Ts_full(SizeTs) )
  ALLOCATE( Yps_full(SizeYps) )
  
  Ds_full = Ds_bary
  Ts_full = Ts_bary
  Yps_full = Yps_bary
  
#ifdef EOSMODE_COMPOSE
  
  ALLOCATE( Ps_full(SizeDs, SizeTs, SizeYps) )
  ALLOCATE( Ss_full(SizeDs, SizeTs, SizeYps) )
  ALLOCATE( Es_full(SizeDs, SizeTs, SizeYps) )
           
  ! Build full EOS for comparison purposes
  DO iD=1,SizeDs
    DO iT=1,SizeTs
      DO iYp=1,SizeYps

        Ymu_temp = Yps_bary(iYp) / Yp_over_Ymu
        Ymu_temp = 0.0d0

        ! Now add electron photon and muon component
        CALL GetPhotonLeptonGasEOS( &
          Ds_bary(iD), Ts_bary(iT), Yps_bary(iYp) - Ymu_temp, Ymu_temp, &
          HelmTableElectrons, HelmTableMuons, P_LeptPhot, E_LeptPhot, S_LeptPhot)
                
        Es_full(iD,iT,iYp) = LOG10( 10.0d0**( Es_bary(iD,iT,iYp) ) + E_LeptPhot)
        Ps_full(iD,iT,iYp) = LOG10( Ps_bary(iD,iT,iYp) + P_LeptPhot)
        Ss_full(iD,iT,iYp) = LOG10( 10.0d0**( Ss_bary(iD,iT,iYp) ) + S_LeptPhot)

        IF (Es_full(iD,iT,iYp) .NE. Es_full(iD,iT,iYp)) THEN
            WRITE(*,*) 'E', iD,iT,iYp, Es_bary(iD,iT,iYp), E_LeptPhot
          STOP
        ENDIF
        IF (Ps_full(iD,iT,iYp) .NE. Ps_full(iD,iT,iYp)) THEN
          WRITE(*,*) 'P', iD,iT,iYp, Ps_bary(iD,iT,iYp), P_LeptPhot
          STOP
        ENDIF

        IF (Ss_full(iD,iT,iYp) .NE. Ss_full(iD,iT,iYp)) THEN
          WRITE(*,*) 'S', iD,iT,iYp, Ss_bary(iD,iT,iYp), S_LeptPhot
          STOP
        ENDIF
                
      ENDDO
    ENDDO
  ENDDO
  
  CALL InitializeEOSComponentsInversion &
         ( Ds_bary, Ts_bary, Yps_bary, &
           10.0d0**( Es_full ) - OS_E, &
           10.0d0**( Ps_full ) - OS_P, &
           10.0d0**( Ss_full ) - OS_S, &
           HelmTableElectrons, &
           HelmTableMuons, &
           Verbose_Option = .TRUE. )

#else
  
  CALL InitializeEOSInversion &
         ( Ds_bary, Ts_bary, Yps_bary, &
           10.0d0**( Es_bary ) - OS_E, &
           10.0d0**( Ps_bary ) - OS_P, &
           10.0d0**( Ss_bary ) - OS_S, &
           Verbose_Option = .TRUE. )

#endif

  ! --- Random Sampling of EOSTable Table ---

  ! --- Initialize Points ---

  ! This crashes
  D  = 1.4627E+14
  E  = 5.7620E+19
  Ye = 0.43603993
  Ym = 0.0d0

  ! This does not crash but is messed up
  D  = 2.937198130e+13
  E  = 3.542233480e+19
  Ye = 0.43603993
  Ym = 0.0d0

  ! Messed up in a different way
  D  = 2.447488340e+13
  E  = 3.363482380e+19
  Ye = 0.43603993
  Ym = 0.0d0

  ! This crashes
  D  = 1.6079E+14
  E  = 6.0164E+19
  Ye = 0.43639
  Ym = 0.0d0

  ! This does not 
  T_Guess = 3.894870830e+10

#ifdef EOSMODE_COMPOSE
  CALL ComputeTemperatureWith_DEYpYl_Single_Guess &
          ( D, E, Ye, Ym, &
          Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, &
          T_E, T_Guess, Error_E, HelmTableElectrons, HelmTableMuons )
#else
    CALL ComputeTemperatureWith_DEY_Single_Guess &
          ( D, E, Ye, &
          Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, &
          T_E, T_Guess, Error_E )
#endif

  IF( Error_E .NE. 0 ) THEN 
    WRITE(*,*) 'Inversion Not Succesfull with guess'
    WRITE(*,*) iP, D, E, Ye, Ym
#ifdef EOSMODE_COMPOSE
    CALL DescribeEOSComponentsInversionError( Error_E )
#else
    CALL DescribeEOSInversionError( Error_E )
#endif
  ELSE
    WRITE(*,*) 'Inversion Succesfull with guess'
    WRITE(*,*) D, E, Ye, Ym, T_E
  ENDIF

#ifdef EOSMODE_COMPOSE
  CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess &
          ( D, E, Ye, Ym, &
          Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, &
          T_E, Error_E, HelmTableElectrons, HelmTableMuons )
#else
    CALL ComputeTemperatureWith_DEY_Single_NoGuess &
          ( D, E, Ye, &
          Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, &
          T_E, Error_E )
#endif

  IF( Error_E .NE. 0 ) THEN 
    WRITE(*,*) 'Inversion Not Succesfull without guess'
    WRITE(*,*) iP, D, E, Ye, Ym
#ifdef EOSMODE_COMPOSE
    CALL DescribeEOSComponentsInversionError( Error_E )
#else
    CALL DescribeEOSInversionError( Error_E )
#endif
  ELSE
    WRITE(*,*) 'Inversion Succesfull without guess'
    WRITE(*,*) D, E, Ye, Ym, T_E
  ENDIF

END PROGRAM wlEOSComposeInversionTestOnePointEnergy

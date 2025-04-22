MODULE wlEosTemperatureSeparateInversionModule

  USE wlKindModule, ONLY: &
    dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log
  USE wlMuonEOS, ONLY: &
    MuonStateType, FullMuonEOS
  USE wlElectronPhotonEOS, ONLY: &
    ElectronPhotonStateType, ElectronPhotonEOS
  USE wlLeptonEOSModule, ONLY: &
    HelmholtzTableType, MuonEOSType
  USE wlHelmMuonIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF, ReadMuonTableHDF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeLeptonTables
  PUBLIC :: InvertTemperatureWith_DEYpYl_Guess
  PUBLIC :: InvertTemperatureWith_DEYpYl_NoGuess
  PUBLIC :: InvertTemperatureWith_DPYpYl_Guess
  PUBLIC :: InvertTemperatureWith_DPYpYl_NoGuess
  PUBLIC :: InvertTemperatureWith_DSYpYl_Guess
  PUBLIC :: InvertTemperatureWith_DSYpYl_NoGuess
  
  TYPE(HelmholtzTableType) :: HelmholtzTable 
  TYPE(MuonEOSType) :: MuonTable
  
CONTAINS

  SUBROUTINE InitializeLeptonTables( HelmholtzTableName, MuonTableName )

    CHARACTER(len=*), INTENT(IN) :: HelmholtzTableName, MuonTableName

    CALL ReadHelmholtzTableHDF( HelmholtzTable, HelmholtzTableName )
    CALL ReadMuonTableHDF( MuonTable, MuonTableName )

  END SUBROUTINE InitializeLeptonTables

  REAL(dp) FUNCTION InverseLogInterp( x_a, x_b, y_a, y_b, y, OS )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x_a, x_b, y_a, y_b, y, OS

    InverseLogInterp &
      = 10.0_dp**( LOG10( x_a ) + LOG10( x_b/x_a ) &
                 * LOG10( (y+OS)/(y_a+OS) ) / LOG10( (y_b+OS)/(y_a+OS) ) )

    RETURN
  END FUNCTION InverseLogInterp
  
  SUBROUTINE InvertTemperatureWith_DEYpYl_Guess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs,  OS, T, T_Guess, Error )

#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error

    REAL(dp) :: Yp
    
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y
    INTEGER  :: SizeDs, SizeTs, SizeYps
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: Xb_a, Xb_b, Xb_c, Xb_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Yps_i(2)
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)
    
    ! Electron and Muon quantities
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState

    Yp = Ye + Ym
    
    T = 0.0_dp
    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYps = SIZE( Yps )

    iD = Index1D_Log( D,       Ds )
    iT = Index1D_Log( T_Guess, Ts )
    iYp = Index1D_Lin( Yp,       Yps )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Yps_i = Yps(iYp:iYp+1)
    ! -------------------------------------------------------------------

    ! --- First Check if Initial Guess Gives a Solution ---
    i_a = iT
    T_a = Ts(i_a)
    Xbs_a = Xs(iD:iD+1,i_a,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)
     
    X_a = Xb_a + ElectronPhotonState % e + MuonState % e
    f_a = X - X_a
    
    i_b = i_a + 1
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)
    
    X_b = Xb_b + ElectronPhotonState % e + MuonState % e 

    f_b = X - X_b
    
    IF ( f_a * f_b <= 0.0_dp ) THEN
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xbs_a = Xs(iD:iD+1,i_a,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)
    
    X_a = Xb_a + ElectronPhotonState % e + MuonState % e

    f_a = X - X_a
    
    i_b = SizeTs
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    X_b = Xb_b + ElectronPhotonState % e + ElectronPhotonState % e

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      i_c = MIN( MAX( i_a + 1, iT ), i_b - 1 )
      DO WHILE ( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Xbs_c = Xs(iD:iD+1,i_c,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t   = T_c
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye  = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t     = T_c
        MuonState % rhoym = D * Ym
        CALL FullMuonEOS(MuonTable, MuonState)

        X_c = Xb_c + ElectronPhotonState % e + MuonState % e 

        f_c = X - X_c

        IF ( f_a * f_c <= 0.0_dp ) THEN
          i_b = i_c
          T_b = T_c
          X_b = X_c
          f_b = f_c
        ELSE
          i_a = i_c
          T_a = T_c
          X_a = X_c
          f_a = f_c
        END IF

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )

      END DO

    ELSE

      ! --- Pick Root Nearest Initial Guess ---
      d_i = SizeTs
      i_c = MIN( MAX( 1, iT ), SizeTs - 1 )
      T_c = T_a
      X_c = X_a
      DO i = 2, SizeTs - 1

        T_i = Ts(i)
        Xbs_i = Xs(iD:iD+1,i,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_i, Xb_i )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_i
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t     = T_i
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)
        
        X_i = Xb_i + ElectronPhotonState % e + MuonState % e 

        f_c = X - X_c
        f_b = X - X_i

        d_c = ABS( i - i_c )

        ! WRITE(*,*) i, X, X_i, X_c

        IF ( f_c * f_b <= 0.0_dp .AND. d_c < d_i ) THEN
          d_i = d_c
          i_a = i - 1
          T_a = T_c
          X_a = X_c
          i_b = i
          T_b = T_i
          X_b = X_i
        ELSE
          T_c = T_i
          X_c = X_i
        END IF

      END DO

      IF ( d_i >= SizeTs ) THEN 
          WRITE(*,*) 'Error in Inversion', X, X_a, X_b, T_a, T_b
          Error = 13
        ENDIF
    END IF

    IF ( Error .NE. 0 ) THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
    END IF

  END SUBROUTINE InvertTemperatureWith_DEYpYl_Guess
  
  SUBROUTINE InvertTemperatureWith_DEYpYl_NoGuess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, OS, T, Error )

#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    
    REAL(dp), INTENT(out) :: T
    INTEGER,  INTENT(out) :: Error

    REAL(dp) :: Yp
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y
    INTEGER  :: SizeDs, SizeTs, SizeYps
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: Xb_a, Xb_b, Xb_c, Xb_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Yps_i(2)
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)
    
    ! Electron and Muon quantities
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState

    Yp = Ye + Ym

    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYps = SIZE( Yps )

    iD = Index1D_Log( D, Ds )
    iYp = Index1D_Lin( Yp, Yps )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Yps_i = Yps(iYp:iYp+1)

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xbs_a = Xs(iD:iD+1,i_a,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)

    X_a = Xb_a + ElectronPhotonState % e + MuonState % e 

    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
    
    X_b = Xb_b + ElectronPhotonState % e + MuonState % e 

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      DO WHILE ( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Xbs_c = Xs(iD:iD+1,i_c,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_c
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_c
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)
 
        X_c = Xb_c + ElectronPhotonState % e + MuonState % e 

        f_c = X - X_c

        IF ( f_a * f_c <= 0.0_dp ) THEN
          i_b = i_c
          T_b = T_c
          X_b = X_c
          f_b = f_c
        ELSE
          i_a = i_c
          T_a = T_c
          X_a = X_c
          f_a = f_c
        END IF

      END DO

    ELSE

      ! --- Pick Highest Temperature Root ---
      DO i = SizeTs - 1, 2, -1

        T_i = Ts(i)
        Xbs_i = Xs(iD:iD+1,i,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_i, Xb_i )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_i
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_i
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)

        X_i = Xb_i + ElectronPhotonState % e + MuonState % e 

        f_a = X - X_i
        f_b = X - X_b

        IF ( f_a * f_b <= 0.0_dp ) THEN
          i_a = i
          T_a = T_i
          X_a = X_i
          EXIT
        ELSE
          i_b = i
          T_b = T_i
          X_b = X_i
        END IF

      END DO

      f_a = X - X_a
      f_b = X - X_b

      IF ( f_a * f_b > 0.0_dp ) THEN
        Error = 13
      ENDIF
    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
    END IF

  END SUBROUTINE InvertTemperatureWith_DEYpYl_NoGuess

  
  SUBROUTINE InvertTemperatureWith_DPYpYl_Guess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, OS, T, T_Guess, Error )

#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error

    REAL(dp) :: Yp
    REAL(dp) :: OS_local
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y
    INTEGER  :: SizeDs, SizeTs, SizeYps
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: Xb_a, Xb_b, Xb_c, Xb_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Yps_i(2)
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)
    
    ! Electron and Muon quantities
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState

    Yp = Ye + Ym
    
    T = 0.0_dp
    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYps = SIZE( Yps )

    iD = Index1D_Log( D,       Ds )
    iT = Index1D_Log( T_Guess, Ts )
    iYp = Index1D_Lin( Yp,       Yps )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Yps_i = Yps(iYp:iYp+1)
    ! -------------------------------------------------------------------

    ! --- First Check if Initial Guess Gives a Solution ---
    i_a = iT
    T_a = Ts(i_a)
    ! Find Local offset
    OS_local =  MINVAL( Xs(iD:iD+1,i_a,iYp:iYp+1) )
    IF (OS_local .lt. 0.0_dp) THEN
      OS_local = -1.1d0*OS_local
    ELSE
      OS_local = 0.0_dp
    ENDIF
    Xbs_a = LOG10(Xs(iD:iD+1,i_a,iYp:iYp+1) + OS_local)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)
    
    X_a = Xb_a + ElectronPhotonState % p + MuonState % p

    f_a = X - X_a
    
    i_b = i_a + 1
    T_b = Ts(i_b)
    OS_local =  MINVAL( Xs(iD:iD+1,i_b,iYp:iYp+1) )
    IF (OS_local .lt. 0.0_dp) THEN
      OS_local = -1.1d0*OS_local
    ELSE
      OS_local = 0.0_dp
    ENDIF
    Xbs_b = LOG10(Xs(iD:iD+1,i_b,iYp:iYp+1) + OS_local)

    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    X_b = Xb_b + ElectronPhotonState % p + MuonState % p

    f_b = X - X_b
    
    IF ( f_a * f_b <= 0.0_dp ) THEN
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    OS_local =  MINVAL( Xs(iD:iD+1,i_a,iYp:iYp+1) )
    IF (OS_local .lt. 0.0_dp) THEN
      OS_local = -1.1d0*OS_local
    ELSE
      OS_local = 0.0_dp
    ENDIF
    Xbs_a = LOG10(Xs(iD:iD+1,i_a,iYp:iYp+1) + OS_local)

    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    X_a = Xb_a + ElectronPhotonState % p + MuonState % p

    f_a = X - X_a
    
    i_b = SizeTs
    T_b = Ts(i_b)
    OS_local =  MINVAL( Xs(iD:iD+1,i_b,iYp:iYp+1) )
    IF (OS_local .lt. 0.0_dp) THEN
      OS_local = -1.1d0*OS_local
    ELSE
      OS_local = 0.0_dp
    ENDIF
    Xbs_b = LOG10(Xs(iD:iD+1,i_b,iYp:iYp+1) + OS_local)

    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    X_b = Xb_b + ElectronPhotonState % p + MuonState % p

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      i_c = MIN( MAX( i_a + 1, iT ), i_b - 1 )
      DO WHILE ( i_b > i_a + 1 )

        T_c = Ts(i_c)
        OS_local =  MINVAL( Xs(iD:iD+1,i_c,iYp:iYp+1) )
        IF (OS_local .lt. 0.0_dp) THEN
          OS_local = -1.1d0*OS_local
        ELSE
          OS_local = 0.0_dp
        ENDIF
        Xbs_c = LOG10(Xs(iD:iD+1,i_c,iYp:iYp+1) + OS_local)
            
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t   = T_c
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye  = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t     = T_c
        MuonState % rhoym = D * Ym
        CALL FullMuonEOS(MuonTable, MuonState)

        X_c = Xb_c + ElectronPhotonState % p + MuonState % p

        f_c = X - X_c

        IF ( f_a * f_c <= 0.0_dp ) THEN
          i_b = i_c
          T_b = T_c
          X_b = X_c
          f_b = f_c
        ELSE
          i_a = i_c
          T_a = T_c
          X_a = X_c
          f_a = f_c
        END IF

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )

      END DO

    ELSE

      ! --- Pick Root Nearest Initial Guess ---
      d_i = SizeTs
      i_c = MIN( MAX( 1, iT ), SizeTs - 1 )
      T_c = T_a
      X_c = X_a
      DO i = 2, SizeTs - 1

        T_i = Ts(i)
        OS_local =  MINVAL( Xs(iD:iD+1,i,iYp:iYp+1) )
        IF (OS_local .lt. 0.0_dp) THEN
          OS_local = -1.1d0*OS_local
        ELSE
          OS_local = 0.0_dp
        ENDIF
        Xbs_i = LOG10(Xs(iD:iD+1,i,iYp:iYp+1) + OS_local)

        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_i, Xb_i )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_i
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_i
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)

        X_i = Xb_i + ElectronPhotonState % p + MuonState % p

        f_c = X - X_c
        f_b = X - X_i

        d_c = ABS( i - i_c )

        ! WRITE(*,*) i, X, X_i, X_c

        IF ( f_c * f_b <= 0.0_dp .AND. d_c < d_i ) THEN
          d_i = d_c
          i_a = i - 1
          T_a = T_c
          X_a = X_c
          i_b = i
          T_b = T_i
          X_b = X_i
        ELSE
          T_c = T_i
          X_c = X_i
        END IF

      END DO

      IF ( d_i >= SizeTs ) THEN 
          Error = 13
      ENDIF
    END IF

    IF ( Error .NE. 0 ) THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
    END IF

  END SUBROUTINE InvertTemperatureWith_DPYpYl_Guess
  
  SUBROUTINE InvertTemperatureWith_DPYpYl_NoGuess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, OS, T, Error )

#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    
    REAL(dp), INTENT(out) :: T
    INTEGER,  INTENT(out) :: Error

    REAL(dp) :: Yp
    REAL(dp) :: OS_local
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y
    INTEGER  :: SizeDs, SizeTs, SizeYps
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: Xb_a, Xb_b, Xb_c, Xb_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Yps_i(2)
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)

    ! Electron and Muon quantities
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState

    Yp = Ye + Ym

    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYps = SIZE( Yps )

    iD = Index1D_Log( D, Ds )
    iYp = Index1D_Lin( Yp, Yps )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Yps_i = Yps(iYp:iYp+1)

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    OS_local =  MINVAL( Xs(iD:iD+1,i_a,iYp:iYp+1) )
    IF (OS_local .lt. 0.0_dp) THEN
      OS_local = -1.1d0*OS_local
    ELSE
      OS_local = 0.0_dp
    ENDIF
    Xbs_a = LOG10(Xs(iD:iD+1,i_a,iYp:iYp+1) + OS_local)

    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
   
    X_a = Xb_a + ElectronPhotonState % p + MuonState % p

    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    OS_local =  MINVAL( Xs(iD:iD+1,i_b,iYp:iYp+1) )
    IF (OS_local .lt. 0.0_dp) THEN
      OS_local = -1.1d0*OS_local
    ELSE
      OS_local = 0.0_dp
    ENDIF
    Xbs_b = LOG10(Xs(iD:iD+1,i_b,iYp:iYp+1) + OS_local)

    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
  
    X_b = Xb_b + ElectronPhotonState % p + MuonState % p

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      DO WHILE ( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        OS_local =  MINVAL( Xs(iD:iD+1,i_c,iYp:iYp+1) )
        IF (OS_local .lt. 0.0_dp) THEN
          OS_local = -1.1d0*OS_local
        ELSE
          OS_local = 0.0_dp
        ENDIF
        Xbs_c = LOG10(Xs(iD:iD+1,i_c,iYp:iYp+1) + OS_local)

        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_c
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_c
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)

        X_c = Xb_c + ElectronPhotonState % p + MuonState % p

        f_c = X - X_c

        IF ( f_a * f_c <= 0.0_dp ) THEN
          i_b = i_c
          T_b = T_c
          X_b = X_c
          f_b = f_c
        ELSE
          i_a = i_c
          T_a = T_c
          X_a = X_c
          f_a = f_c
        END IF

      END DO

    ELSE

      ! --- Pick Highest Temperature Root ---
      DO i = SizeTs - 1, 2, -1

        T_i = Ts(i)
        OS_local =  MINVAL( Xs(iD:iD+1,i,iYp:iYp+1) )
        IF (OS_local .lt. 0.0_dp) THEN
          OS_local = -1.1d0*OS_local
        ELSE
          OS_local = 0.0_dp
        ENDIF
        Xbs_i = LOG10(Xs(iD:iD+1,i,iYp:iYp+1) + OS_local)

        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS_local, Xbs_i, Xb_i )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_i
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_i
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)

        X_i = Xb_i + ElectronPhotonState % p + MuonState % p

        f_a = X - X_i
        f_b = X - X_b

        IF ( f_a * f_b <= 0.0_dp ) THEN
          i_a = i
          T_a = T_i
          X_a = X_i
          EXIT
        ELSE
          i_b = i
          T_b = T_i
          X_b = X_i
        END IF

      END DO

      f_a = X - X_a
      f_b = X - X_b

      IF ( f_a * f_b > 0.0_dp ) Error = 13

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS_local )
    END IF

  END SUBROUTINE InvertTemperatureWith_DPYpYl_NoGuess

  
  SUBROUTINE InvertTemperatureWith_DSYpYl_Guess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, OS, T, T_Guess, Error )

#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error

    REAL(dp) :: Yp
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y
    INTEGER  :: SizeDs, SizeTs, SizeYps
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: Xb_a, Xb_b, Xb_c, Xb_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Yps_i(2)
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)

    ! Electron and Muon quantities
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState

    Yp = Ye + Ym
    
    T = 0.0_dp
    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYps = SIZE( Yps )

    iD = Index1D_Log( D,       Ds )
    iT = Index1D_Log( T_Guess, Ts )
    iYp = Index1D_Lin( Yp,       Yps )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Yps_i = Yps(iYp:iYp+1)
    ! -------------------------------------------------------------------

    ! --- First Check if Initial Guess Gives a Solution ---
    i_a = iT
    T_a = Ts(i_a)
    Xbs_a = Xs(iD:iD+1,i_a,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    X_a = Xb_a + ElectronPhotonState % s + MuonState % s

    f_a = X - X_a
    
    i_b = i_a + 1
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    X_b = Xb_b + ElectronPhotonState % s + MuonState % s

    f_b = X - X_b
    
    IF ( f_a * f_b <= 0.0_dp ) THEN
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X,  OS )
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xbs_a = Xs(iD:iD+1,i_a,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    X_a = Xb_a + ElectronPhotonState % s + MuonState % s

    f_a = X - X_a
    
    i_b = SizeTs
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    X_b = Xb_b + ElectronPhotonState % s + MuonState % s

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      i_c = MIN( MAX( i_a + 1, iT ), i_b - 1 )
      DO WHILE ( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Xbs_c = Xs(iD:iD+1,i_c,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t   = T_c
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye  = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t     = T_c
        MuonState % rhoym = D * Ym
        CALL FullMuonEOS(MuonTable, MuonState)
    
        X_c = Xb_c + ElectronPhotonState % s + MuonState % s

        f_c = X - X_c

        IF ( f_a * f_c <= 0.0_dp ) THEN
          i_b = i_c
          T_b = T_c
          X_b = X_c
          f_b = f_c
        ELSE
          i_a = i_c
          T_a = T_c
          X_a = X_c
          f_a = f_c
        END IF

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )

      END DO

    ELSE

      ! --- Pick Root Nearest Initial Guess ---
      d_i = SizeTs
      i_c = MIN( MAX( 1, iT ), SizeTs - 1 )
      T_c = T_a
      X_c = X_a
      DO i = 2, SizeTs - 1

        T_i = Ts(i)
        Xbs_i = Xs(iD:iD+1,i,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_i, Xb_i )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_i
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_i
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)

        X_i = Xb_i + ElectronPhotonState % s + MuonState % s

        f_c = X - X_c
        f_b = X - X_i

        d_c = ABS( i - i_c )

        IF ( f_c * f_b <= 0.0_dp .AND. d_c < d_i ) THEN
          d_i = d_c
          i_a = i - 1
          T_a = T_c
          X_a = X_c
          i_b = i
          T_b = T_i
          X_b = X_i
        ELSE
          T_c = T_i
          X_c = X_i
        END IF

      END DO

      IF ( d_i >= SizeTs ) THEN 
          Error = 13
      ENDIF
    END IF

    IF ( Error .NE. 0 ) THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X,  OS )
    END IF

  END SUBROUTINE InvertTemperatureWith_DSYpYl_Guess
  
  SUBROUTINE InvertTemperatureWith_DSYpYl_NoGuess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs,  OS, T, Error )

#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  ::  OS
    
    REAL(dp), INTENT(out) :: T
    INTEGER,  INTENT(out) :: Error

    REAL(dp) :: Yp
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y
    INTEGER  :: SizeDs, SizeTs, SizeYps
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: Xb_a, Xb_b, Xb_c, Xb_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Yps_i(2)
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)
    
    ! Electron and Muon quantities
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState

    Yp = Ye + Ym

    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYps = SIZE( Yps )

    iD = Index1D_Log( D, Ds )
    iYp = Index1D_Lin( Yp, Yps )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Yps_i = Yps(iYp:iYp+1)

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xbs_a = Xs(iD:iD+1,i_a,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_a
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_a
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
   
    X_a = Xb_a + ElectronPhotonState % s + MuonState % s

    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    ElectronPhotonState % t = T_b
    ElectronPhotonState % rho = D
    ElectronPhotonState % ye = Ye
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    MuonState % t = T_b
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
  
    X_b = Xb_b + ElectronPhotonState % s + MuonState % s

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      DO WHILE ( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Xbs_c = Xs(iD:iD+1,i_c,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_c
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_c
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)

        X_c = Xb_c + ElectronPhotonState % s + MuonState % s

        f_c = X - X_c

        IF ( f_a * f_c <= 0.0_dp ) THEN
          i_b = i_c
          T_b = T_c
          X_b = X_c
          f_b = f_c
        ELSE
          i_a = i_c
          T_a = T_c
          X_a = X_c
          f_a = f_c
        END IF

      END DO

    ELSE

      ! --- Pick Highest Temperature Root ---
      DO i = SizeTs - 1, 2, -1

        T_i = Ts(i)
        Xbs_i = Xs(iD:iD+1,i,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_i, Xb_i )

        ! Calculate electron and muon contribution
        ElectronPhotonState % t = T_i
        ElectronPhotonState % rho = D
        ElectronPhotonState % ye = Ye
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_i
        MuonState % rhoym = D * Ym
        
        CALL FullMuonEOS(MuonTable, MuonState)

        X_i = Xb_i + ElectronPhotonState % s + MuonState % s

        f_a = X - X_i
        f_b = X - X_b

        IF ( f_a * f_b <= 0.0_dp ) THEN
          i_a = i
          T_a = T_i
          X_a = X_i
          EXIT
        ELSE
          i_b = i
          T_b = T_i
          X_b = X_i
        END IF

      END DO

      f_a = X - X_a
      f_b = X - X_b

      IF ( f_a * f_b > 0.0_dp ) Error = 13

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X,  OS )
    END IF

  END SUBROUTINE InvertTemperatureWith_DSYpYl_NoGuess

END MODULE wlEosTemperatureSeparateInversionModule
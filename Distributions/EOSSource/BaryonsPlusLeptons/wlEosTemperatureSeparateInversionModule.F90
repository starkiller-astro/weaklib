MODULE wlEosTemperatureSeparateInversionModule

  USE wlKindModule, ONLY: &
    dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: &
    GetPhotonLeptonGasEOS, &
    GetPhotonLeptonGasEOS_dT

  IMPLICIT NONE
  PRIVATE

  INTEGER,  PARAMETER :: iVar_E = 1
  INTEGER,  PARAMETER :: iVar_P = 2
  INTEGER,  PARAMETER :: iVar_S = 3

  ! --- Tolerances for the temperature refinement within a table cell ---
  INTEGER,  PARAMETER :: MaxIter_Refine = 30
  REAL(dp), PARAMETER :: Tol_LnT        = 1.0d-12
  REAL(dp), PARAMETER :: Tol_X          = 1.0d-13
  REAL(dp), PARAMETER :: ln10           = 2.302585092994045684_dp

  PUBLIC :: InvertTemperatureWith_DEYpYl_Guess
  PUBLIC :: InvertTemperatureWith_DEYpYl_NoGuess
  PUBLIC :: InvertTemperatureWith_DPYpYl_Guess
  PUBLIC :: InvertTemperatureWith_DPYpYl_NoGuess
  PUBLIC :: InvertTemperatureWith_DSYpYl_Guess
  PUBLIC :: InvertTemperatureWith_DSYpYl_NoGuess
  
CONTAINS

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

  ! Given a bracketing table interval [ Ts(i_a), Ts(i_a+1) ], solve
  !   X_bary(D,T,Yp) + X_LeptPhot(D,T,Ye,Ym) = X
  ! for T, where X_bary is evaluated exactly as in the forward call
  ! (trilinear interpolation of log10(X+OS) in (log D, log T, Yp)) and
  ! X_LeptPhot is evaluated at the exact point. Safeguarded Newton in ln T.
  ! On input T is the initial guess, on output the refined temperature.
  SUBROUTINE RefineTemperatureInCell &
    ( iVar, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
      HelmTableElectrons, HelmTableMuons )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)    :: iVar, iD, i_a, iYp
    REAL(dp), INTENT(in)    :: D, X, Ye, Ym
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)    :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(inout) :: T
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons

    INTEGER  :: iL_T, iter
    REAL(dp) :: Xbs(2,2,2), L(2)
    REAL(dp) :: OS_cell, Yp, dD, dY
    REAL(dp) :: L_a, dL, LnT_a, dLnT
    REAL(dp) :: LnT, LnT_lo, LnT_hi, LnT_new, LnT_best
    REAL(dp) :: f, f_lo, f_hi, f_best, dfdLnT
    LOGICAL  :: LinearCell

    Yp = Ye + Ym

    ! --- Baryonic table values on the cell ---
    LinearCell = .FALSE.
    OS_cell    = OS
    IF ( iVar == iVar_P ) THEN
      ! Pressure is stored linearly and can be negative. As in the forward
      ! interpolation (wlSoundSpeedModule), it is interpolated linearly in
      ! (log D, log T, Yp) if it is <= 0 at some corner of the cell, and
      ! log-trilinearly otherwise
      OS_cell    = 0.0_dp
      LinearCell = MINVAL( Xs(iD:iD+1,i_a:i_a+1,iYp:iYp+1) ) <= 0.0_dp
      IF ( LinearCell ) THEN
        Xbs = Xs(iD:iD+1,i_a:i_a+1,iYp:iYp+1)
      ELSE
        Xbs = LOG10( Xs(iD:iD+1,i_a:i_a+1,iYp:iYp+1) )
      END IF
    ELSE
      Xbs = Xs(iD:iD+1,i_a:i_a+1,iYp:iYp+1)
    END IF

    dD = LOG10( D / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
    dY = ( Yp - Yps(iYp) ) / ( Yps(iYp+1) - Yps(iYp) )

    ! --- Bilinear interpolation in (log D, Yp) on the two T slices ---
    DO iL_T = 1, 2
      L(iL_T) &
        =   ( 1.0_dp - dY ) * ( ( 1.0_dp - dD ) * Xbs(1,iL_T,1) + dD * Xbs(2,iL_T,1) ) &
          +            dY   * ( ( 1.0_dp - dD ) * Xbs(1,iL_T,2) + dD * Xbs(2,iL_T,2) )
    END DO

    L_a   = L(1)
    dL    = L(2) - L(1)
    LnT_a = LOG( Ts(i_a) )
    dLnT  = LOG( Ts(i_a+1) ) - LnT_a

    ! --- Check that the consistent interpolant brackets the root ---
    LnT_lo = LnT_a
    LnT_hi = LnT_a + dLnT
    CALL ComputeCellResidual &
           ( iVar, LnT_lo, D, X, Ye, Ym, L_a, dL, LnT_a, dLnT, OS_cell, LinearCell, &
             HelmTableElectrons, HelmTableMuons, f_lo, dfdLnT )
    CALL ComputeCellResidual &
           ( iVar, LnT_hi, D, X, Ye, Ym, L_a, dL, LnT_a, dLnT, OS_cell, LinearCell, &
             HelmTableElectrons, HelmTableMuons, f_hi, dfdLnT )

    IF ( f_lo == 0.0_dp ) THEN
      T = Ts(i_a)
      RETURN
    ELSE IF ( f_hi == 0.0_dp ) THEN
      T = Ts(i_a+1)
      RETURN
    ELSE IF ( f_lo * f_hi > 0.0_dp ) THEN
      ! Keep the initial guess
      RETURN
    END IF

    LnT = MIN( MAX( LOG( T ), LnT_lo ), LnT_hi )

    ! --- Best iterate so far (returned if the iteration does not converge) ---
    LnT_best = LnT_lo
    f_best   = f_lo
    IF ( ABS( f_hi ) < ABS( f_best ) ) THEN
      LnT_best = LnT_hi
      f_best   = f_hi
    END IF

    DO iter = 1, MaxIter_Refine

      CALL ComputeCellResidual &
             ( iVar, LnT, D, X, Ye, Ym, L_a, dL, LnT_a, dLnT, OS_cell, LinearCell, &
               HelmTableElectrons, HelmTableMuons, f, dfdLnT )

      IF ( ABS( f ) < ABS( f_best ) ) THEN
        LnT_best = LnT
        f_best   = f
      END IF

      IF ( ABS( f ) <= Tol_X * ABS( X ) ) EXIT

      ! --- Update bracket ---
      IF ( f * f_lo < 0.0_dp ) THEN
        LnT_hi = LnT
        f_hi   = f
      ELSE
        LnT_lo = LnT
        f_lo   = f
      END IF

      ! --- Newton step. A step below tolerance means convergence: accept it
      !     without the bracket test, since near the root (at roundoff level)
      !     it can fall just outside the bracket and would trigger a spurious
      !     bisection jump. Otherwise bisect if the step leaves the bracket.
      IF ( dfdLnT /= 0.0_dp ) THEN
        LnT_new = LnT - f / dfdLnT
        IF ( ABS( LnT_new - LnT ) < Tol_LnT ) THEN
          LnT_best = LnT_new
          EXIT
        END IF
      END IF
      IF ( dfdLnT == 0.0_dp .OR. LnT_new <= LnT_lo .OR. LnT_new >= LnT_hi ) THEN
        LnT_new = 0.5_dp * ( LnT_lo + LnT_hi )
      END IF

      ! --- Bracket collapsed: nothing left to gain ---
      IF ( LnT_hi - LnT_lo < Tol_LnT ) EXIT

      LnT = LnT_new

    END DO

    T = EXP( LnT_best )

  END SUBROUTINE RefineTemperatureInCell


  ! Residual f = X_bary + X_LeptPhot - X and its derivative w.r.t. ln T
  ! inside one table temperature cell. X_bary + OS_cell is 10**L, with L
  ! linear in ln T, or L itself if LinearCell
  SUBROUTINE ComputeCellResidual &
    ( iVar, LnT, D, X, Ye, Ym, L_a, dL, LnT_a, dLnT, OS_cell, LinearCell, &
      HelmTableElectrons, HelmTableMuons, f, dfdLnT )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iVar
    REAL(dp), INTENT(in)  :: LnT, D, X, Ye, Ym
    REAL(dp), INTENT(in)  :: L_a, dL, LnT_a, dLnT, OS_cell
    LOGICAL,  INTENT(in)  :: LinearCell
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons
    REAL(dp), INTENT(out) :: f, dfdLnT

    REAL(dp) :: T, Xb_plus_OS
    REAL(dp) :: P_LeptPhot, E_LeptPhot, S_LeptPhot
    REAL(dp) :: dPdT_LeptPhot, dEdT_LeptPhot, dSdT_LeptPhot

    T = EXP( LnT )

    IF ( LinearCell ) THEN
      Xb_plus_OS = L_a + dL * ( LnT - LnT_a ) / dLnT
      dfdLnT     = dL / dLnT
    ELSE
      Xb_plus_OS = 10.0_dp**( L_a + dL * ( LnT - LnT_a ) / dLnT )
      dfdLnT     = ln10 * Xb_plus_OS * dL / dLnT
    END IF

    CALL GetPhotonLeptonGasEOS_dT( &
      D, T, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot, &
      dPdT_LeptPhot, dEdT_LeptPhot, dSdT_LeptPhot)

    SELECT CASE ( iVar )
    CASE ( iVar_E )
      f      = Xb_plus_OS - OS_cell + E_LeptPhot - X
      dfdLnT = dfdLnT + T * dEdT_LeptPhot
    CASE ( iVar_P )
      f      = Xb_plus_OS - OS_cell + P_LeptPhot - X
      dfdLnT = dfdLnT + T * dPdT_LeptPhot
    CASE DEFAULT
      f      = Xb_plus_OS - OS_cell + S_LeptPhot - X
      dfdLnT = dfdLnT + T * dSdT_LeptPhot
    END SELECT

  END SUBROUTINE ComputeCellResidual


  ! Baryon pressure on one table temperature slice, interpolated in
  ! (log D, Yp): linearly if it is <= 0 at some corner, log-bilinearly
  ! otherwise (same choice as the forward interpolation)
  REAL(dp) FUNCTION BaryonPressureSlice( LogD, Yp, LogDs_i, Yps_i, Ps )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: LogD, Yp, LogDs_i(2), Yps_i(2), Ps(2,2)

    REAL(dp) :: dD, dY

    dD = ( LogD - LogDs_i(1) ) / ( LogDs_i(2) - LogDs_i(1) )
    dY = ( Yp   - Yps_i  (1) ) / ( Yps_i  (2) - Yps_i  (1) )

    IF ( MINVAL( Ps ) <= 0.0_dp ) THEN
      BaryonPressureSlice &
        =   ( 1.0_dp - dY ) * ( ( 1.0_dp - dD ) * Ps(1,1) + dD * Ps(2,1) ) &
          +            dY   * ( ( 1.0_dp - dD ) * Ps(1,2) + dD * Ps(2,2) )
    ELSE
      BaryonPressureSlice &
        = 10.0_dp**(   ( 1.0_dp - dY ) * (   ( 1.0_dp - dD ) * LOG10( Ps(1,1) )   &
                                         +            dD   * LOG10( Ps(2,1) ) ) &
                     +            dY   * (   ( 1.0_dp - dD ) * LOG10( Ps(1,2) )   &
                                         +            dD   * LOG10( Ps(2,2) ) ) )
    END IF

  END FUNCTION BaryonPressureSlice
  
  SUBROUTINE InvertTemperatureWith_DEYpYl_Guess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs,  OS, T, T_Guess, Error, HelmTableElectrons, HelmTableMuons )

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
    REAL(dp) :: P_LeptPhot, E_LeptPhot, S_LeptPhot
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)
    
    ! Electron and Muon quantities
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons

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
     
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_a = Xb_a + E_LeptPhot
    f_a = X - X_a
    
    i_b = i_a + 1
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)
    
    X_b = Xb_b + E_LeptPhot 
    f_b = X - X_b
    
    IF ( f_a * f_b <= 0.0_dp ) THEN
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
      CALL RefineTemperatureInCell &
             ( iVar_E, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xbs_a = Xs(iD:iD+1,i_a,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)
    
    X_a = Xb_a + E_LeptPhot
    f_a = X - X_a
    
    i_b = SizeTs
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_b = Xb_b + E_LeptPhot
    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      i_c = MIN( MAX( i_a + 1, iT ), i_b - 1 )
      DO WHILE ( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Xbs_c = Xs(iD:iD+1,i_c,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        CALL GetPhotonLeptonGasEOS( &
          D, T_c, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_c = Xb_c + E_LeptPhot 
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
        CALL GetPhotonLeptonGasEOS( &
          D, T_i, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)
        
        X_i = Xb_i + E_LeptPhot 

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
      CALL RefineTemperatureInCell &
             ( iVar_E, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
    END IF

  END SUBROUTINE InvertTemperatureWith_DEYpYl_Guess
  
  SUBROUTINE InvertTemperatureWith_DEYpYl_NoGuess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, OS, T, Error, HelmTableElectrons, HelmTableMuons )

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
    REAL(dp) :: P_LeptPhot, E_LeptPhot, S_LeptPhot
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)
    
    ! Electron and Muon quantities
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons

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
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_a = Xb_a + E_LeptPhot 
    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)
    
    X_b = Xb_b + E_LeptPhot 

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      DO WHILE ( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Xbs_c = Xs(iD:iD+1,i_c,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        CALL GetPhotonLeptonGasEOS( &
          D, T_c, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_c = Xb_c + E_LeptPhot 
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
        CALL GetPhotonLeptonGasEOS( &
          D, T_i, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_i = Xb_i + E_LeptPhot 

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
      CALL RefineTemperatureInCell &
             ( iVar_E, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
    END IF

  END SUBROUTINE InvertTemperatureWith_DEYpYl_NoGuess

  
  SUBROUTINE InvertTemperatureWith_DPYpYl_Guess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, OS, T, T_Guess, Error, HelmTableElectrons, HelmTableMuons )

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
    REAL(dp) :: P_LeptPhot, E_LeptPhot, S_LeptPhot
    
    ! Electron and Muon quantities
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons

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
    Xb_a = BaryonPressureSlice &
             ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i_a,iYp:iYp+1) )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)
    
    X_a = Xb_a + P_LeptPhot
    f_a = X - X_a
    
    i_b = i_a + 1
    T_b = Ts(i_b)
    Xb_b = BaryonPressureSlice &
             ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i_b,iYp:iYp+1) )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_b = Xb_b + P_LeptPhot
    f_b = X - X_b
    
    IF ( f_a * f_b <= 0.0_dp ) THEN
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
      CALL RefineTemperatureInCell &
             ( iVar_P, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xb_a = BaryonPressureSlice &
             ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i_a,iYp:iYp+1) )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_a = Xb_a + P_LeptPhot
    f_a = X - X_a
    
    i_b = SizeTs
    T_b = Ts(i_b)
    Xb_b = BaryonPressureSlice &
             ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i_b,iYp:iYp+1) )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_b = Xb_b + P_LeptPhot
    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      i_c = MIN( MAX( i_a + 1, iT ), i_b - 1 )
      DO WHILE ( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Xb_c = BaryonPressureSlice &
                 ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i_c,iYp:iYp+1) )

        ! Calculate electron and muon contribution
        CALL GetPhotonLeptonGasEOS( &
          D, T_c, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_c = Xb_c + P_LeptPhot
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
        Xb_i = BaryonPressureSlice &
                 ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i,iYp:iYp+1) )

        ! Calculate electron and muon contribution
        CALL GetPhotonLeptonGasEOS( &
          D, T_i, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_i = Xb_i + P_LeptPhot

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
      CALL RefineTemperatureInCell &
             ( iVar_P, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
    END IF

  END SUBROUTINE InvertTemperatureWith_DPYpYl_Guess
  
  SUBROUTINE InvertTemperatureWith_DPYpYl_NoGuess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, OS, T, Error, HelmTableElectrons, HelmTableMuons )

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
    REAL(dp) :: P_LeptPhot, E_LeptPhot, S_LeptPhot

    ! Electron and Muon quantities
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons
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
    Xb_a = BaryonPressureSlice &
             ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i_a,iYp:iYp+1) )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_a = Xb_a + P_LeptPhot
    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Xb_b = BaryonPressureSlice &
             ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i_b,iYp:iYp+1) )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_b = Xb_b + P_LeptPhot
    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      DO WHILE ( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Xb_c = BaryonPressureSlice &
                 ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i_c,iYp:iYp+1) )

        ! Calculate electron and muon contribution
        CALL GetPhotonLeptonGasEOS( &
          D, T_c, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_c = Xb_c + P_LeptPhot
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
        Xb_i = BaryonPressureSlice &
                 ( LogD, Yp, LogDs_i, Yps_i, Xs(iD:iD+1,i,iYp:iYp+1) )

        ! Calculate electron and muon contribution
        CALL GetPhotonLeptonGasEOS( &
          D, T_i, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_i = Xb_i + P_LeptPhot

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
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
      CALL RefineTemperatureInCell &
             ( iVar_P, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
    END IF

  END SUBROUTINE InvertTemperatureWith_DPYpYl_NoGuess

  
  SUBROUTINE InvertTemperatureWith_DSYpYl_Guess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, OS, T, T_Guess, Error, HelmTableElectrons, HelmTableMuons )

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
    REAL(dp) :: P_LeptPhot, E_LeptPhot, S_LeptPhot
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)

    ! Electron and Muon quantities
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons
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
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_a = Xb_a + S_LeptPhot
    f_a = X - X_a
    
    i_b = i_a + 1
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_b = Xb_b + S_LeptPhot
    f_b = X - X_b
    
    IF ( f_a * f_b <= 0.0_dp ) THEN
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X,  OS )
      CALL RefineTemperatureInCell &
             ( iVar_S, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xbs_a = Xs(iD:iD+1,i_a,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_a, Xb_a )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_a = Xb_a + S_LeptPhot
    f_a = X - X_a
    
    i_b = SizeTs
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_b = Xb_b + S_LeptPhot
    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      i_c = MIN( MAX( i_a + 1, iT ), i_b - 1 )
      DO WHILE ( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Xbs_c = Xs(iD:iD+1,i_c,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        CALL GetPhotonLeptonGasEOS( &
          D, T_c, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_c = Xb_c + S_LeptPhot
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
        CALL GetPhotonLeptonGasEOS( &
          D, T_i, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_i = Xb_i + S_LeptPhot

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
      CALL RefineTemperatureInCell &
             ( iVar_S, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
    END IF

  END SUBROUTINE InvertTemperatureWith_DSYpYl_Guess
  
  SUBROUTINE InvertTemperatureWith_DSYpYl_NoGuess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs,  OS, T, Error, HelmTableElectrons, HelmTableMuons )

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
    REAL(dp) :: P_LeptPhot, E_LeptPhot, S_LeptPhot
    REAL(dp) :: Xbs_a(2,2), Xbs_b(2,2), Xbs_c(2,2), Xbs_i(2,2)
    
    ! Electron and Muon quantities
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons
    
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
    CALL GetPhotonLeptonGasEOS( &
      D, T_a, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_a = Xb_a + S_LeptPhot
    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Xbs_b = Xs(iD:iD+1,i_b,iYp:iYp+1)
    
    CALL LogInterpolateSingleVariable &
           ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_b, Xb_b )

    ! Calculate electron and muon contribution
    CALL GetPhotonLeptonGasEOS( &
      D, T_b, Ye, Ym, &
      HelmTableElectrons, HelmTableMuons, &
      P_LeptPhot, E_LeptPhot, S_LeptPhot)

    X_b = Xb_b + S_LeptPhot
    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      DO WHILE ( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Xbs_c = Xs(iD:iD+1,i_c,iYp:iYp+1)
        
        CALL LogInterpolateSingleVariable &
               ( LogD, Yp, LogDs_i, Yps_i,  OS, Xbs_c, Xb_c )

        ! Calculate electron and muon contribution
        CALL GetPhotonLeptonGasEOS( &
          D, T_c, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_c = Xb_c + S_LeptPhot
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
        CALL GetPhotonLeptonGasEOS( &
          D, T_i, Ye, Ym, &
          HelmTableElectrons, HelmTableMuons, &
          P_LeptPhot, E_LeptPhot, S_LeptPhot)

        X_i = Xb_i + S_LeptPhot

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
      CALL RefineTemperatureInCell &
             ( iVar_S, D, X, Ye, Ym, iD, i_a, iYp, Ds, Ts, Yps, Xs, OS, T, &
               HelmTableElectrons, HelmTableMuons )
    END IF

  END SUBROUTINE InvertTemperatureWith_DSYpYl_NoGuess

END MODULE wlEosTemperatureSeparateInversionModule
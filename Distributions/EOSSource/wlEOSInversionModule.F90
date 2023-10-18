MODULE wlEOSInversionModule

  USE wlKindModule, ONLY: &
    dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_2D_Custom_Point
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeEOSInversion
  PUBLIC :: ComputeTemperatureWith_DEY
  PUBLIC :: ComputeTemperatureWith_DEY_Many
  PUBLIC :: ComputeTemperatureWith_DEY_Single
  PUBLIC :: ComputeTemperatureWith_DEY_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DEY_Single_Guess_Error
  PUBLIC :: ComputeTemperatureWith_DEY_Single_Guess_NoError
  PUBLIC :: ComputeTemperatureWith_DEY_Single_NoGuess
  PUBLIC :: ComputeTemperatureWith_DEY_Single_NoGuess_Error
  PUBLIC :: ComputeTemperatureWith_DEY_Single_NoGuess_NoError
  PUBLIC :: ComputeTemperatureWith_DPY
  PUBLIC :: ComputeTemperatureWith_DPY_Many
  PUBLIC :: ComputeTemperatureWith_DPY_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DPY_Single_Guess_Error
  PUBLIC :: ComputeTemperatureWith_DPY_Single_Guess_NoError
  PUBLIC :: ComputeTemperatureWith_DPY_Single_NoGuess
  PUBLIC :: ComputeTemperatureWith_DPY_Single_NoGuess_Error
  PUBLIC :: ComputeTemperatureWith_DPY_Single_NoGuess_NoError
  PUBLIC :: ComputeTemperatureWith_DSY
  PUBLIC :: ComputeTemperatureWith_DSY_Many
  PUBLIC :: ComputeTemperatureWith_DSY_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DSY_Single_Guess_Error
  PUBLIC :: ComputeTemperatureWith_DSY_Single_Guess_NoError
  PUBLIC :: ComputeTemperatureWith_DSY_Single_NoGuess
  PUBLIC :: ComputeTemperatureWith_DSY_Single_NoGuess_Error
  PUBLIC :: ComputeTemperatureWith_DSY_Single_NoGuess_NoError
  PUBLIC :: DescribeEOSInversionError

  LOGICAL  :: InversionInitialized
  REAL(dp), PUBLIC :: MinD, MaxD
  REAL(dp), PUBLIC :: MinT, MaxT
  REAL(dp), PUBLIC :: MinY, MaxY
  REAL(dp), PUBLIC :: MinE, MaxE
  REAL(dp), PUBLIC :: MinP, MaxP
  REAL(dp), PUBLIC :: MinS, MaxS

#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET( &
    !$OMP   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$OMP   MinY, MaxY, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#elif defined(WEAKLIB_OACC)
    !$ACC DECLARE CREATE( &
    !$ACC   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$ACC   MinY, MaxY, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#endif

  INTERFACE ComputeTemperatureWith_DEY
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Many
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DEY

  INTERFACE ComputeTemperatureWith_DEY_Single
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DEY_Single

  INTERFACE ComputeTemperatureWith_DEY_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_Guess_NoError
  END INTERFACE ComputeTemperatureWith_DEY_Single_Guess

  INTERFACE ComputeTemperatureWith_DEY_Single_NoGuess
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DEY_Single_NoGuess

  INTERFACE ComputeTemperatureWith_DPY
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Many
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DPY

  INTERFACE ComputeTemperatureWith_DPY_Single
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DPY_Single

  INTERFACE ComputeTemperatureWith_DPY_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_Guess_NoError
  END INTERFACE ComputeTemperatureWith_DPY_Single_Guess

  INTERFACE ComputeTemperatureWith_DPY_Single_NoGuess
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DPY_Single_NoGuess

  INTERFACE ComputeTemperatureWith_DSY
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Many
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DSY

  INTERFACE ComputeTemperatureWith_DSY_Single
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DSY_Single

  INTERFACE ComputeTemperatureWith_DSY_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_Guess_NoError
  END INTERFACE ComputeTemperatureWith_DSY_Single_Guess

  INTERFACE ComputeTemperatureWith_DSY_Single_NoGuess
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DSY_Single_NoGuess


CONTAINS


  SUBROUTINE InitializeEOSInversion( Ds, Ts, Ys, Es, Ps, Ss, Verbose_Option )

    REAL(dp), INTENT(in) :: Ds(1:)      , Ts(1:)      , Ys(1:)
    REAL(dp), INTENT(in) :: Es(1:,1:,1:), Ps(1:,1:,1:), Ss(1:,1:,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    MinD = MINVAL( Ds ); MaxD = MAXVAL( Ds )
    MinT = MINVAL( Ts ); MaxT = MAXVAL( Ts )
    MinY = MINVAL( Ys ); MaxY = MAXVAL( Ys )
    MinE = MINVAL( Es ); MaxE = MAXVAL( Es )
    MinP = MINVAL( Ps ); MaxP = MAXVAL( Ps )
    MinS = MINVAL( Ss ); MaxS = MAXVAL( Ss )

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'InitializeEOSInversion'
      WRITE(*,*)
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max D   [g cm^-3] = ', MinD, MaxD
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max T         [K] = ', MinT, MaxT
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max Y             = ', MinY, MaxY
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max E  [erg g^-1] = ', MinE, MaxE
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max P [dyn cm^-2] = ', MinP, MaxP
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max S       [k_B] = ', MinS, MaxS
      WRITE(*,*)
    END IF

    InversionInitialized = .TRUE.

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET UPDATE TO( &
    !$OMP   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$OMP   MinY, MaxY, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#elif defined(WEAKLIB_OACC)
    !$ACC UPDATE DEVICE( &
    !$ACC   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$ACC   MinY, MaxY, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#endif

  END SUBROUTINE InitializeEOSInversion


  INTEGER FUNCTION CheckInputError &
    ( D, X, Y, MinX, MaxX )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: D, X, Y, MinX, MaxX

    CheckInputError = 0

    IF ( .NOT. InversionInitialized ) THEN
      CheckInputError = 10
    ELSE IF ( D /= D .OR. X /= X .OR. Y /= Y ) THEN
      CheckInputError = 11
    ELSE IF ( D < MinD .OR. D > MaxD ) THEN
      CheckInputError = 01
#if EOS_DEBUG
      WRITE(*,*) 'ComputeTemperature ERROR 01'
      WRITE(*,*) 'input D, E/P/S, Y :', D, X, Y
      WRITE(*,*) 'Table Min/Max D:', MinD, MaxD
#endif
    ELSE IF ( X < MinX .OR. X > MaxX ) THEN
      CheckInputError = 02
#if EOS_DEBUG
      WRITE(*,*) 'ComputeTemperature ERROR 02'
      WRITE(*,*) 'input D, E/P/S, Y :', D, X, Y
      WRITE(*,*) 'Table Min/Max E/P/S:', MinX, MaxX
#endif
    ELSE IF ( Y < MinY .OR. Y > MaxY ) THEN
      CheckInputError = 03
#if EOS_DEBUG
      WRITE(*,*) 'ComputeTemperature ERROR 03'
      WRITE(*,*) 'input D, E/P/S, Y :', D, X, Y
      WRITE(*,*) 'Table Min/Max Y:', MinY, MaxY
#endif
    END IF

  END FUNCTION CheckInputError


  SUBROUTINE DescribeEOSInversionError( Error )

    INTEGER, INTENT(in) :: Error

    CHARACTER(64) :: ErrorString(00:13)

    IF( Error > 13 ) STOP 'ERROR in EOSInversionError flag'

    ErrorString(00) = 'Returned Successfully'
    ErrorString(01) = 'First Argument (D) Outside Table Bounds'
    ErrorString(02) = 'Second Argument (E, P, or S) Outside Table Bounds'
    ErrorString(03) = 'Third Argument (Y) Outside Table Bounds'
    ErrorString(10) = 'EOS Inversion Not Initialized'
    ErrorString(11) = 'NAN in Argument(s)'
    ErrorString(13) = 'Unable to Find Any Root'

    WRITE(*,*)
    WRITE(*,*) '  wlEOSInversionModule ERROR: ' // TRIM( ErrorString(Error) )
    WRITE(*,*)

  END SUBROUTINE DescribeEOSInversionError


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


  SUBROUTINE ComputeTemperatureWith_DXY_Guess &
    ( D, X, Y, Ds, Ts, Ys, Xs, OS, T, T_Guess, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error

    INTEGER  :: iD, iT, iY
    INTEGER  :: SizeDs, SizeTs, SizeYs
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Ys_i(2)
    REAL(dp) :: Xs_a(2,2), Xs_b(2,2), Xs_c(2,2), Xs_i(2,2)

    ! -------------------------------------------------------------------

    T = 0.0_dp
    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D_Log( D,       Ds )
    iT = Index1D_Log( T_Guess, Ts )
    iY = Index1D_Lin( Y,       Ys )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iY = MIN( MAX( 1, iY ), SizeYs - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Ys_i = Ys(iY:iY+1)

    ! -------------------------------------------------------------------

    ! --- First Check if Initial Guess Gives a Solution ---

    i_a = iT
    T_a = Ts(i_a)
    Xs_a = Xs(iD:iD+1,i_a,iY:iY+1)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Y, LogDs_i, Ys_i, OS, Xs_a, X_a )

    f_a = X - X_a

    i_b = i_a + 1
    T_b = Ts(i_b)
    Xs_b = Xs(iD:iD+1,i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Y, LogDs_i, Ys_i, OS, Xs_b, X_b )

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xs_a = Xs(iD:iD+1,i_a,iY:iY+1)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Y, LogDs_i, Ys_i, OS, Xs_a, X_a )

    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Xs_b = Xs(iD:iD+1,i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Y, LogDs_i, Ys_i, OS, Xs_b, X_b )

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      i_c = MIN( MAX( i_a + 1, iT ), i_b - 1 )
      DO WHILE ( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Xs_c = Xs(iD:iD+1,i_c,iY:iY+1)

        CALL LogInterpolateSingleVariable_2D_Custom_Point &
               ( LogD, Y, LogDs_i, Ys_i, OS, Xs_c, X_c )

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
        Xs_i = Xs(iD:iD+1,i,iY:iY+1)

        CALL LogInterpolateSingleVariable_2D_Custom_Point &
               ( LogD, Y, LogDs_i, Ys_i, OS, Xs_i, X_i )

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

      IF ( d_i >= SizeTs ) Error = 13

    END IF

    IF ( Error .NE. 0 ) THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DXY_Guess


  SUBROUTINE ComputeTemperatureWith_DXY_NoGuess &
    ( D, X, Y, Ds, Ts, Ys, Xs, OS, T, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    INTEGER,  INTENT(out) :: Error

    INTEGER  :: iD, iY
    INTEGER  :: SizeDs, SizeTs, SizeYs
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Ys_i(2)
    REAL(dp) :: Xs_a(2,2), Xs_b(2,2), Xs_c(2,2), Xs_i(2,2)

    ! -------------------------------------------------------------------
    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D_Log( D, Ds )
    iY = Index1D_Lin( Y, Ys )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iY = MIN( MAX( 1, iY ), SizeYs - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Ys_i = Ys(iY:iY+1)

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xs_a = Xs(iD:iD+1,i_a,iY:iY+1)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Y, LogDs_i, Ys_i, OS, Xs_a, X_a )

    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Xs_b = Xs(iD:iD+1,i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Y, LogDs_i, Ys_i, OS, Xs_b, X_b )

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      DO WHILE ( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Xs_c = Xs(iD:iD+1,i_c,iY:iY+1)

        CALL LogInterpolateSingleVariable_2D_Custom_Point &
               ( LogD, Y, LogDs_i, Ys_i, OS, Xs_c, X_c )

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
        Xs_i = Xs(iD:iD+1,i,iY:iY+1)

        CALL LogInterpolateSingleVariable_2D_Custom_Point &
               ( LogD, Y, LogDs_i, Ys_i, OS, Xs_i, X_i )

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
    END IF

  END SUBROUTINE ComputeTemperatureWith_DXY_NoGuess


  SUBROUTINE ComputeTemperatureWith_DEY_Many &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T, UseInitialGuess_Option, Error_Option )

    REAL(dp), INTENT(in)    :: D (1:), E (1:), Y (1:)
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(inout) :: T(1:)
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(1:)

    REAL(dp) :: T_Guess
    LOGICAL  :: UseInitialGuess
    INTEGER  :: Error(SIZE(D))
    INTEGER  :: i

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

    Error = 0

    DO i = 1, SIZE( D )
      IF ( UseInitialGuess ) THEN
        T_Guess = T(i)
        CALL ComputeTemperatureWith_DEY_Single_Guess_Error &
               ( D(i), E(i), Y(i), Ds, Ts, Ys, Es, OS, T(i), T_Guess, Error(i) )
      ELSE
        CALL ComputeTemperatureWith_DEY_Single_NoGuess_Error &
               ( D(i), E(i), Y(i), Ds, Ts, Ys, Es, OS, T(i), Error(i) )
      END IF
    END DO

    IF( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE ComputeTemperatureWith_DEY_Many


  SUBROUTINE ComputeTemperatureWith_DEY_Single_Guess_Error &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T, T_Guess, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , E     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error

    T = 0.0_dp
    Error = CheckInputError( D, E, Y, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_Guess &
             ( D, E, Y, Ds, Ts, Ys, Es, Os, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Single_Guess_Error


  SUBROUTINE ComputeTemperatureWith_DEY_Single_Guess_NoError &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T, T_Guess )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , E     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess

    INTEGER  :: Error

    T = 0.0_dp
    Error = CheckInputError( D, E, Y, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_Guess &
             ( D, E, Y, Ds, Ts, Ys, Es, Os, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Single_Guess_NoError


  SUBROUTINE ComputeTemperatureWith_DEY_Single_NoGuess_Error &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , E     , Y
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out)   :: Error

    T = 0.0_dp
    Error = CheckInputError( D, E, Y, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_NoGuess &
             ( D, E, Y, Ds, Ts, Ys, Es, Os, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Single_NoGuess_Error


  SUBROUTINE ComputeTemperatureWith_DEY_Single_NoGuess_NoError &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , E     , Y
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T

    INTEGER  :: Error

    T = 0.0_dp
    Error = CheckInputError( D, E, Y, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_NoGuess &
             ( D, E, Y, Ds, Ts, Ys, Es, Os, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Single_NoGuess_NoError


  SUBROUTINE ComputeTemperatureWith_DPY_Many &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T, UseInitialGuess_Option, Error_Option )

    REAL(dp), INTENT(in)    :: D (1:), P (1:), Y (1:)
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(inout) :: T(1:)
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(1:)

    REAL(dp) :: T_Guess
    LOGICAL  :: UseInitialGuess
    INTEGER  :: Error(SIZE(D))
    INTEGER  :: i

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

    Error = 0

    DO i = 1, SIZE( D )
      IF ( UseInitialGuess ) THEN
        T_Guess = T(i)
        CALL ComputeTemperatureWith_DPY_Single_Guess_Error &
               ( D(i), P(i), Y(i), Ds, Ts, Ys, Ps, OS, T(i), T_Guess, Error(i) )
      ELSE
        CALL ComputeTemperatureWith_DPY_Single_NoGuess_Error &
               ( D(i), P(i), Y(i), Ds, Ts, Ys, Ps, OS, T(i), Error(i) )
      END IF
    END DO

    IF( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE ComputeTemperatureWith_DPY_Many


  SUBROUTINE ComputeTemperatureWith_DPY_Single_Guess_Error &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T, T_Guess, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , P     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error

    T = 0.0_dp
    Error = CheckInputError( D, P, Y, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_Guess &
             ( D, P, Y, Ds, Ts, Ys, Ps, Os, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Single_Guess_Error


  SUBROUTINE ComputeTemperatureWith_DPY_Single_Guess_NoError &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T, T_Guess )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , P     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess

    INTEGER  :: Error

    T = 0.0_dp
    Error = CheckInputError( D, P, Y, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_Guess &
             ( D, P, Y, Ds, Ts, Ys, Ps, Os, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Single_Guess_NoError


  SUBROUTINE ComputeTemperatureWith_DPY_Single_NoGuess_Error &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , P     , Y
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out)   :: Error

    T = 0.0_dp
    Error = CheckInputError( D, P, Y, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_NoGuess &
             ( D, P, Y, Ds, Ts, Ys, Ps, Os, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Single_NoGuess_Error


  SUBROUTINE ComputeTemperatureWith_DPY_Single_NoGuess_NoError &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , P     , Y
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T

    INTEGER  :: Error

    T = 0.0_dp
    Error = CheckInputError( D, P, Y, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_NoGuess &
             ( D, P, Y, Ds, Ts, Ys, Ps, Os, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Single_NoGuess_NoError


  SUBROUTINE ComputeTemperatureWith_DSY_Many &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T, UseInitialGuess_Option, Error_Option )

    REAL(dp), INTENT(in)    :: D (1:), S (1:), Y (1:)
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(inout) :: T(1:)
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(1:)

    REAL(dp) :: T_Guess
    LOGICAL  :: UseInitialGuess
    INTEGER  :: Error(SIZE(D))
    INTEGER  :: i

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

    Error = 0

    DO i = 1, SIZE( D )
      IF ( UseInitialGuess ) THEN
        T_Guess = T(i)
        CALL ComputeTemperatureWith_DSY_Single_Guess_Error &
               ( D(i), S(i), Y(i), Ds, Ts, Ys, Ss, OS, T(i), T_Guess, Error(i) )
      ELSE
        CALL ComputeTemperatureWith_DSY_Single_NoGuess_Error &
               ( D(i), S(i), Y(i), Ds, Ts, Ys, Ss, OS, T(i), Error(i) )
      END IF
    END DO

    IF( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE ComputeTemperatureWith_DSY_Many


  SUBROUTINE ComputeTemperatureWith_DSY_Single_Guess_Error &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T, T_Guess, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , S     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error

    T = 0.0_dp
    Error = CheckInputError( D, S, Y, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_Guess &
             ( D, S, Y, Ds, Ts, Ys, Ss, Os, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Single_Guess_Error


  SUBROUTINE ComputeTemperatureWith_DSY_Single_Guess_NoError &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T, T_Guess )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , S     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess

    INTEGER  :: Error

    T = 0.0_dp
    Error = CheckInputError( D, S, Y, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_Guess &
             ( D, S, Y, Ds, Ts, Ys, Ss, Os, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Single_Guess_NoError


  SUBROUTINE ComputeTemperatureWith_DSY_Single_NoGuess_Error &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , S     , Y
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out)   :: Error

    T = 0.0_dp
    Error = CheckInputError( D, S, Y, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_NoGuess &
             ( D, S, Y, Ds, Ts, Ys, Ss, Os, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Single_NoGuess_Error


  SUBROUTINE ComputeTemperatureWith_DSY_Single_NoGuess_NoError &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , S     , Y
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)    :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T

    INTEGER  :: Error

    T = 0.0_dp
    Error = CheckInputError( D, S, Y, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXY_NoGuess &
             ( D, S, Y, Ds, Ts, Ys, Ss, Os, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Single_NoGuess_NoError


END MODULE wlEOSInversionModule

MODULE wlEOSInversionModule

  USE wlKindModule, ONLY: &
    dp
  USE wlInterpolationModule, ONLY: &
    Index1D, &
    LogInterpolateSingleVariable

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeEOSInversion
  PUBLIC :: ComputeTemperatureWith_DEY
  PUBLIC :: ComputeTemperatureWith_DEY_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DPY
  PUBLIC :: ComputeTemperatureWith_DPY_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DSY
  PUBLIC :: DescribeEOSInversionError

  LOGICAL  :: InversionInitialized
  REAL(dp) :: MinD, MaxD
  REAL(dp) :: MinT, MaxT
  REAL(dp) :: MinY, MaxY
  REAL(dp) :: MinE, MaxE
  REAL(dp) :: MinP, MaxP
  REAL(dp) :: MinS, MaxS

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
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single_NoGuess
  END INTERFACE ComputeTemperatureWith_DEY

  INTERFACE ComputeTemperatureWith_DPY
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Many
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single_NoGuess
  END INTERFACE ComputeTemperatureWith_DPY

  INTERFACE ComputeTemperatureWith_DSY
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Many
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single_NoGuess
  END INTERFACE ComputeTemperatureWith_DSY

CONTAINS


  SUBROUTINE InitializeEOSInversion( Ds, Ts, Ys, Es, Ps, Ss, Verbose_Option )

    REAL(dp), INTENT(in) :: Ds(:)
    REAL(dp), INTENT(in) :: Ts(:)
    REAL(dp), INTENT(in) :: Ys(:)
    REAL(dp), INTENT(in) :: Es(:,:,:)
    REAL(dp), INTENT(in) :: Ps(:,:,:)
    REAL(dp), INTENT(in) :: Ss(:,:,:)
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
      WRITE(*,'(A6,A24,2ES10.3E2)') '', 'Min/Max D   [g cm^-3] = ', MinD, MaxD
      WRITE(*,'(A6,A24,2ES10.3E2)') '', 'Min/Max T         [K] = ', MinT, MaxT
      WRITE(*,'(A6,A24,2ES10.3E2)') '', 'Min/Max Y             = ', MinY, MaxY
      WRITE(*,'(A6,A24,2ES10.3E2)') '', 'Min/Max E  [erg g^-1] = ', MinE, MaxE
      WRITE(*,'(A6,A24,2ES10.3E2)') '', 'Min/Max P [dyn cm^-2] = ', MinP, MaxP
      WRITE(*,'(A6,A24,2ES10.3E2)') '', 'Min/Max S       [k_B] = ', MinS, MaxS
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


  SUBROUTINE ComputeTemperatureWith_DEY_Many &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T, UseInitialGuess_Option, Error_Option )

    REAL(dp), INTENT(in)    :: D(:)
    REAL(dp), INTENT(in)    :: E(:)
    REAL(dp), INTENT(in)    :: Y(:)
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Es(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(inout) :: T(:)
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(:)

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
        CALL ComputeTemperatureWith_DEY_Single_Guess &
               ( D(i), E(i), Y(i), Ds, Ts, Ys, Es, OS, T(i), T_Guess, &
                 Error_Option = Error(i) )

      ELSE

        CALL ComputeTemperatureWith_DEY_Single_NoGuess &
               ( D(i), E(i), Y(i), Ds, Ts, Ys, Es, OS, T(i), &
                 Error_Option = Error(i) )

      END IF

    END DO

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Many


  SUBROUTINE ComputeTemperatureWith_DEY_Single_Guess &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T, T_Guess, Error_Option )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D
    REAL(dp), INTENT(in)    :: E
    REAL(dp), INTENT(in)    :: Y
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Es(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    REAL(dp), INTENT(in)    :: T_Guess
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER  :: Error
    INTEGER  :: iD, iT, iY, lo, hi
    INTEGER  :: i, d_c, d_i
    INTEGER  :: SizeDs, SizeTs, SizeYs
    INTEGER  :: i_a, i_b, i_c
    REAL(dp) :: T_a, T_b, T_c
    REAL(dp) :: E_a, E_b, E_c
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: T_i, E_i
    REAL(dp) :: Ds_i(1:2), Ts_i(1:3), Ys_i(1:2)
    REAL(dp) :: Ts_a(1:2), Ts_b(1:2), Ts_c(1:2)
    REAL(dp), DIMENSION(1:2,1:2,1:2) :: Es_a, Es_b, Es_c
    REAL(dp), DIMENSION(1:2,1:3,1:2) :: Es_i

    ! --- Initial Error Check -------------------------------------------

    Error = 0

    IF( .NOT. InversionInitialized )THEN
      Error = 10
    END IF

    IF( D < MinD .OR. D > MaxD )THEN
      Error = 01
    END IF

    IF( E < MinE .OR. E > MaxE )THEN
      Error = 02
    END IF

    IF( Y < MinY .OR. Y > MaxY )THEN
      Error = 03
    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D( D, Ds, SizeDs )
    iT = Index1D( T_Guess, Ts, SizeTs )
    iY = Index1D( Y, Ys, SizeYs )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iY = MIN( MAX( 1, iY ), SizeYs - 1 )

    Ds_i(1:2) = Ds(iD:iD+1)
    Ys_i(1:2) = Ys(iY:iY+1)

    ! -------------------------------------------------------------------

    ! --- First Check if Initial Guess Gives a Solution ---

    i_a = iT
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Es_a(1:2,1:2,1:2) = Es(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Es_a, E_a )

    f_a = E - E_a

    i_b = i_a + 1
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Es_b(1:2,1:2,1:2) = Es(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Es_b, E_b )

    f_b = E - E_b

    IF( f_a * f_b < 0.0_dp )THEN
      T = InverseLogInterp( Ts(i_a), Ts(i_b), E_a, E_b, E, OS )
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Es_a(1:2,1:2,1:2) = Es(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Es_a, E_a )

    f_a = E - E_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Es_b(1:2,1:2,1:2) = Es(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Es_b, E_b )

    f_b = E - E_b

    IF( f_a * f_b < 0.0_dp )THEN

      i_c = MIN( MAX( i_a+1, iT ), i_b-1 )

      DO WHILE( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Ts_c(1:2) = Ts(i_c:i_c+1)
        Es_c(1:2,1:2,1:2) = Es(iD:iD+1,i_c:i_c+1,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_c, Y, Ds_i, Ts_c, Ys_i, OS, Es_c, E_c )

        f_c = E - E_c

        IF( f_a * f_c < 0.0_dp )THEN

          i_b = i_c
          T_b = T_c
          E_b = E_c
          f_b = f_c

        ELSE

          i_a = i_c
          T_a = T_c
          E_a = E_c
          f_a = f_c

        END IF

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )

      END DO

    ELSE

      ! --- Pick Root Nearest Initial Guess ---
      d_i = SizeTs
      i_c = MIN( MAX( 1, iT ), SizeTs - 1 )
      T_c = T_a
      E_c = E_a
      DO i = 2, SizeTs - 1

        lo = i-1
        hi = i+1
        T_i = Ts(i)
        Ts_i(1:3) = Ts(lo:hi)
        Es_i(1:2,1:3,1:2) = Es(iD:iD+1,lo:hi,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_i, Y, Ds_i, Ts_i, Ys_i, OS, Es_i, E_i )

        f_c = E - E_c
        f_b = E - E_i

        d_c = ABS( i - i_c )

        IF( f_c * f_b < 0.0_dp .AND. d_c < d_i )THEN
          d_i = d_c
          i_a = i - 1
          T_a = T_c
          E_a = E_c
          i_b = i
          T_b = T_i
          E_b = E_i
        ELSE
          T_c = T_i
          E_c = E_i
        END IF

      END DO

      IF( d_i >= SizeTs )THEN
        Error = 13
      END IF

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, E_a, E_b, E, OS )
    END IF

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Single_Guess


  SUBROUTINE ComputeTemperatureWith_DEY_Single_NoGuess &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T, Error_Option )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D
    REAL(dp), INTENT(in)    :: E
    REAL(dp), INTENT(in)    :: Y
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Es(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER  :: Error
    INTEGER  :: iD, iY, lo, hi
    INTEGER  :: i
    INTEGER  :: SizeDs, SizeTs, SizeYs
    INTEGER  :: i_a, i_b, i_c
    REAL(dp) :: T_a, T_b, T_c
    REAL(dp) :: E_a, E_b, E_c
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: T_i, E_i
    REAL(dp) :: Ds_i(1:2), Ts_i(1:3), Ys_i(1:2)
    REAL(dp) :: Ts_a(1:2), Ts_b(1:2), Ts_c(1:2)
    REAL(dp), DIMENSION(1:2,1:2,1:2) :: Es_a, Es_b, Es_c
    REAL(dp), DIMENSION(1:2,1:3,1:2) :: Es_i

    ! --- Initial Error Check -------------------------------------------

    Error = 0

    IF( .NOT. InversionInitialized )THEN
      Error = 10
    END IF

    IF( D < MinD .OR. D > MaxD )THEN
      Error = 01
    END IF

    IF( E < MinE .OR. E > MaxE )THEN
      Error = 02
    END IF

    IF( Y < MinY .OR. Y > MaxY )THEN
      Error = 03
    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D( D, Ds, SizeDs )
    iY = Index1D( Y, Ys, SizeYs )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iY = MIN( MAX( 1, iY ), SizeYs - 1 )

    Ds_i(1:2) = Ds(iD:iD+1)
    Ys_i(1:2) = Ys(iY:iY+1)

    i_a = 1
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Es_a(1:2,1:2,1:2) = Es(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Es_a, E_a )

    f_a = E - E_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Es_b(1:2,1:2,1:2) = Es(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Es_b, E_b )

    f_b = E - E_b

    IF( f_a * f_b < 0.0_dp )THEN

      DO WHILE( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Ts_c(1:2) = Ts(i_c:i_c+1)
        Es_c(1:2,1:2,1:2) = Es(iD:iD+1,i_c:i_c+1,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_c, Y, Ds_i, Ts_c, Ys_i, OS, Es_c, E_c )

        f_c = E - E_c

        IF( f_a * f_c < 0.0_dp )THEN

          i_b = i_c
          T_b = T_c
          E_b = E_c
          f_b = f_c

        ELSE

          i_a = i_c
          T_a = T_c
          E_a = E_c
          f_a = f_c

        END IF

      END DO

    ELSE

      ! --- Pick Highest Temperature Root ---
      DO i = SizeTs - 1, 2, -1

        lo = i-1
        hi = i+1
        T_i = Ts(i)
        Ts_i(1:3) = Ts(lo:hi)
        Es_i(1:2,1:3,1:2) = Es(iD:iD+1,lo:hi,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_i, Y, Ds_i, Ts_i, Ys_i, OS, Es_i, E_i )

        f_a = E - E_i
        f_b = E - E_b

        IF( f_a * f_b < 0.0_dp )THEN
          i_a = i
          T_a = T_i
          E_a = E_i
          EXIT
        ELSE
          i_b = i
          T_b = T_i
          E_b = E_i
        END IF

      END DO

      f_a = E - E_a
      f_b = E - E_b

      IF( f_a * f_b >= 0.0_dp )THEN

        Error = 13

      END IF

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, E_a, E_b, E, OS )
    END IF

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Single_NoGuess


  SUBROUTINE ComputeTemperatureWith_DPY_Many &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T, UseInitialGuess_Option, Error_Option )

    REAL(dp), INTENT(in)    :: D(:)
    REAL(dp), INTENT(in)    :: P(:)
    REAL(dp), INTENT(in)    :: Y(:)
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Ps(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(inout) :: T(:)
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(:)

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
        CALL ComputeTemperatureWith_DPY_Single_Guess &
               ( D(i), P(i), Y(i), Ds, Ts, Ys, Ps, OS, T(i), T_Guess, &
                 Error_Option = Error(i) )

      ELSE

        CALL ComputeTemperatureWith_DPY_Single_NoGuess &
               ( D(i), P(i), Y(i), Ds, Ts, Ys, Ps, OS, T(i), &
                 Error_Option = Error(i) )

      END IF

    END DO

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Many


  SUBROUTINE ComputeTemperatureWith_DPY_Single_Guess &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T, T_Guess, Error_Option )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D
    REAL(dp), INTENT(in)    :: P
    REAL(dp), INTENT(in)    :: Y
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Ps(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    REAL(dp), INTENT(in)    :: T_Guess
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER  :: Error
    INTEGER  :: iD, iT, iY, lo, hi
    INTEGER  :: i, d_c, d_i
    INTEGER  :: SizeDs, SizeTs, SizeYs
    INTEGER  :: i_a, i_b, i_c
    REAL(dp) :: T_a, T_b, T_c
    REAL(dp) :: P_a, P_b, P_c
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: T_i, P_i
    REAL(dp) :: Ds_i(1:2), Ts_i(1:3), Ys_i(1:2)
    REAL(dp) :: Ts_a(1:2), Ts_b(1:2), Ts_c(1:2)
    REAL(dp), DIMENSION(1:2,1:2,1:2) :: Ps_a, Ps_b, Ps_c
    REAL(dp), DIMENSION(1:2,1:3,1:2) :: Ps_i

    ! --- Initial Error Check -------------------------------------------

    Error = 0

    IF( .NOT. InversionInitialized )THEN
      Error = 10
    END IF

    IF( D < MinD .OR. D > MaxD )THEN
      Error = 01
    END IF

    IF( P < MinP .OR. P > MaxP )THEN
      Error = 02
    END IF

    IF( Y < MinY .OR. Y > MaxY )THEN
      Error = 03
    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D( D, Ds, SizeDs )
    iT = Index1D( T_Guess, Ts, SizeTs )
    iY = Index1D( Y, Ys, SizeYs )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iY = MIN( MAX( 1, iY ), SizeYs - 1 )

    Ds_i(1:2) = Ds(iD:iD+1)
    Ys_i(1:2) = Ys(iY:iY+1)

    ! -------------------------------------------------------------------

    ! --- First Check if Initial Guess Gives a Solution ---

    i_a = iT
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Ps_a(1:2,1:2,1:2) = Ps(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Ps_a, P_a )

    f_a = P - P_a

    i_b = i_a + 1
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Ps_b(1:2,1:2,1:2) = Ps(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Ps_b, P_b )

    f_b = P - P_b

    IF( f_a * f_b < 0.0_dp )THEN
      T = InverseLogInterp( Ts(i_a), Ts(i_b), P_a, P_b, P, OS )
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Ps_a(1:2,1:2,1:2) = Ps(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Ps_a, P_a )

    f_a = P - P_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Ps_b(1:2,1:2,1:2) = Ps(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Ps_b, P_b )

    f_b = P - P_b

    IF( f_a * f_b < 0.0_dp )THEN

      i_c = MIN( MAX( i_a+1, iT ), i_b-1 )

      DO WHILE( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Ts_c(1:2) = Ts(i_c:i_c+1)
        Ps_c(1:2,1:2,1:2) = Ps(iD:iD+1,i_c:i_c+1,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_c, Y, Ds_i, Ts_c, Ys_i, OS, Ps_c, P_c )

        f_c = P - P_c

        IF( f_a * f_c < 0.0_dp )THEN

          i_b = i_c
          T_b = T_c
          P_b = P_c
          f_b = f_c

        ELSE

          i_a = i_c
          T_a = T_c
          P_a = P_c
          f_a = f_c

        END IF

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )

      END DO

    ELSE

      ! --- Pick Root Nearest Initial Guess ---
      d_i = SizeTs
      i_c = MIN( MAX( 1, iT ), SizeTs - 1 )
      T_c = T_a
      P_c = P_a
      DO i = 2, SizeTs - 1

        lo = i-1
        hi = i+1
        T_i = Ts(i)
        Ts_i(1:3) = Ts(lo:hi)
        Ps_i(1:2,1:3,1:2) = Ps(iD:iD+1,lo:hi,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_i, Y, Ds_i, Ts_i, Ys_i, OS, Ps_i, P_i )

        f_c = P - P_c
        f_b = P - P_i

        d_c = ABS( i - i_c )

        IF( f_c * f_b < 0.0_dp .AND. d_c < d_i )THEN
          d_i = d_c
          i_a = i - 1
          T_a = T_c
          P_a = P_c
          i_b = i
          T_b = T_i
          P_b = P_i
        ELSE
          T_c = T_i
          P_c = P_i
        END IF

      END DO

      IF( d_i >= SizeTs )THEN
        Error = 13
      END IF

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, P_a, P_b, P, OS )
    END IF

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Single_Guess


  SUBROUTINE ComputeTemperatureWith_DPY_Single_NoGuess &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T, Error_Option )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D
    REAL(dp), INTENT(in)    :: P
    REAL(dp), INTENT(in)    :: Y
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Ps(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER  :: Error
    INTEGER  :: iD, iY, lo, hi
    INTEGER  :: i
    INTEGER  :: SizeDs, SizeTs, SizeYs
    INTEGER  :: i_a, i_b, i_c
    REAL(dp) :: T_a, T_b, T_c
    REAL(dp) :: P_a, P_b, P_c
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: T_i, P_i
    REAL(dp) :: Ds_i(1:2), Ts_i(1:3), Ys_i(1:2)
    REAL(dp) :: Ts_a(1:2), Ts_b(1:2), Ts_c(1:2)
    REAL(dp), DIMENSION(1:2,1:2,1:2) :: Ps_a, Ps_b, Ps_c
    REAL(dp), DIMENSION(1:2,1:3,1:2) :: Ps_i

    ! --- Initial Error Check -------------------------------------------

    Error = 0

    IF( .NOT. InversionInitialized )THEN
      Error = 10
    END IF

    IF( D < MinD .OR. D > MaxD )THEN
      Error = 01
    END IF

    IF( P < MinP .OR. P > MaxP )THEN
      Error = 02
    END IF

    IF( Y < MinY .OR. Y > MaxY )THEN
      Error = 03
    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D( D, Ds, SizeDs )
    iY = Index1D( Y, Ys, SizeYs )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iY = MIN( MAX( 1, iY ), SizeYs - 1 )

    Ds_i(1:2) = Ds(iD:iD+1)
    Ys_i(1:2) = Ys(iY:iY+1)

    i_a = 1
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Ps_a(1:2,1:2,1:2) = Ps(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Ps_a, P_a )

    f_a = P - P_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Ps_b(1:2,1:2,1:2) = Ps(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Ps_b, P_b )

    f_b = P - P_b

    IF( f_a * f_b < 0.0_dp )THEN

      DO WHILE( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Ts_c(1:2) = Ts(i_c:i_c+1)
        Ps_c(1:2,1:2,1:2) = Ps(iD:iD+1,i_c:i_c+1,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_c, Y, Ds_i, Ts_c, Ys_i, OS, Ps_c, P_c )

        f_c = P - P_c

        IF( f_a * f_c < 0.0_dp )THEN

          i_b = i_c
          T_b = T_c
          P_b = P_c
          f_b = f_c

        ELSE

          i_a = i_c
          T_a = T_c
          P_a = P_c
          f_a = f_c

        END IF

      END DO

    ELSE

      ! --- Pick Highest Temperature Root ---
      DO i = SizeTs - 1, 2, -1

        lo = i-1
        hi = i+1
        T_i = Ts(i)
        Ts_i(1:3) = Ts(lo:hi)
        Ps_i(1:2,1:3,1:2) = Ps(iD:iD+1,lo:hi,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_i, Y, Ds_i, Ts_i, Ys_i, OS, Ps_i, P_i )

        f_a = P - P_i
        f_b = P - P_b

        IF( f_a * f_b < 0.0_dp )THEN
          i_a = i
          T_a = T_i
          P_a = P_i
          EXIT
        ELSE
          i_b = i
          T_b = T_i
          P_b = P_i
        END IF

      END DO

      f_a = P - P_a
      f_b = P - P_b

      IF( f_a * f_b >= 0.0_dp )THEN

        Error = 13

      END IF

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, P_a, P_b, P, OS )
    END IF

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Single_NoGuess


  SUBROUTINE ComputeTemperatureWith_DSY_Many &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T, UseInitialGuess_Option, Error_Option )

    REAL(dp), INTENT(in)    :: D(:)
    REAL(dp), INTENT(in)    :: S(:)
    REAL(dp), INTENT(in)    :: Y(:)
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Ss(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(inout) :: T(:)
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(:)

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
        CALL ComputeTemperatureWith_DSY_Single_Guess &
               ( D(i), S(i), Y(i), Ds, Ts, Ys, Ss, OS, T(i), T_Guess, &
                 Error_Option = Error(i) )

      ELSE

        CALL ComputeTemperatureWith_DSY_Single_NoGuess &
               ( D(i), S(i), Y(i), Ds, Ts, Ys, Ss, OS, T(i), &
                 Error_Option = Error(i) )

      END IF

    END DO

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Many


  SUBROUTINE ComputeTemperatureWith_DSY_Single_Guess &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T, T_Guess, Error_Option )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D
    REAL(dp), INTENT(in)    :: S
    REAL(dp), INTENT(in)    :: Y
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Ss(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    REAL(dp), INTENT(in)    :: T_Guess
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER  :: Error
    INTEGER  :: iD, iT, iY, lo, hi
    INTEGER  :: i, d_c, d_i
    INTEGER  :: SizeDs, SizeTs, SizeYs
    INTEGER  :: i_a, i_b, i_c
    REAL(dp) :: T_a, T_b, T_c
    REAL(dp) :: S_a, S_b, S_c
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: T_i, S_i
    REAL(dp) :: Ds_i(1:2), Ts_i(1:3), Ys_i(1:2)
    REAL(dp) :: Ts_a(1:2), Ts_b(1:2), Ts_c(1:2)
    REAL(dp), DIMENSION(1:2,1:2,1:2) :: Ss_a, Ss_b, Ss_c
    REAL(dp), DIMENSION(1:2,1:3,1:2) :: Ss_i

    ! --- Initial Error Check -------------------------------------------

    Error = 0

    IF( .NOT. InversionInitialized )THEN
      Error = 10
    END IF

    IF( D < MinD .OR. D > MaxD )THEN
      Error = 01
    END IF

    IF( S < MinS .OR. S > MaxS )THEN
      Error = 02
    END IF

    IF( Y < MinY .OR. Y > MaxY )THEN
      Error = 03
    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D( D, Ds, SizeDs )
    iT = Index1D( T_Guess, Ts, SizeTs )
    iY = Index1D( Y, Ys, SizeYs )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iY = MIN( MAX( 1, iY ), SizeYs - 1 )

    Ds_i(1:2) = Ds(iD:iD+1)
    Ys_i(1:2) = Ys(iY:iY+1)

    ! -------------------------------------------------------------------

    ! --- First Check if Initial Guess Gives a Solution ---

    i_a = iT
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Ss_a(1:2,1:2,1:2) = Ss(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Ss_a, S_a )

    f_a = S - S_a

    i_b = i_a + 1
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Ss_b(1:2,1:2,1:2) = Ss(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Ss_b, S_b )

    f_b = S - S_b

    IF( f_a * f_b < 0.0_dp )THEN
      T = InverseLogInterp( Ts(i_a), Ts(i_b), S_a, S_b, S, OS )
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Ss_a(1:2,1:2,1:2) = Ss(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Ss_a, S_a )

    f_a = S - S_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Ss_b(1:2,1:2,1:2) = Ss(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Ss_b, S_b )

    f_b = S - S_b

    IF( f_a * f_b < 0.0_dp )THEN

      i_c = MIN( MAX( i_a+1, iT ), i_b-1 )

      DO WHILE( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Ts_c(1:2) = Ts(i_c:i_c+1)
        Ss_c(1:2,1:2,1:2) = Ss(iD:iD+1,i_c:i_c+1,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_c, Y, Ds_i, Ts_c, Ys_i, OS, Ss_c, S_c )

        f_c = S - S_c

        IF( f_a * f_c < 0.0_dp )THEN

          i_b = i_c
          T_b = T_c
          S_b = S_c
          f_b = f_c

        ELSE

          i_a = i_c
          T_a = T_c
          S_a = S_c
          f_a = f_c

        END IF

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )

      END DO

    ELSE

      ! --- Pick Root Nearest Initial Guess ---
      d_i = SizeTs
      i_c = MIN( MAX( 1, iT ), SizeTs - 1 )
      T_c = T_a
      S_c = S_a
      DO i = 2, SizeTs - 1

        lo = i-1
        hi = i+1
        T_i = Ts(i)
        Ts_i(1:3) = Ts(lo:hi)
        Ss_i(1:2,1:3,1:2) = Ss(iD:iD+1,lo:hi,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_i, Y, Ds_i, Ts_i, Ys_i, OS, Ss_i, S_i )

        f_c = S - S_c
        f_b = S - S_i

        d_c = ABS( i - i_c )

        IF( f_c * f_b < 0.0_dp .AND. d_c < d_i )THEN
          d_i = d_c
          i_a = i - 1
          T_a = T_c
          S_a = S_c
          i_b = i
          T_b = T_i
          S_b = S_i
        ELSE
          T_c = T_i
          S_c = S_i
        END IF

      END DO

      IF( d_i >= SizeTs )THEN
        Error = 13
      END IF

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, S_a, S_b, S, OS )
    END IF

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Single_Guess


  SUBROUTINE ComputeTemperatureWith_DSY_Single_NoGuess &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T, Error_Option )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D
    REAL(dp), INTENT(in)    :: S
    REAL(dp), INTENT(in)    :: Y
    REAL(dp), INTENT(in)    :: Ds(:)
    REAL(dp), INTENT(in)    :: Ts(:)
    REAL(dp), INTENT(in)    :: Ys(:)
    REAL(dp), INTENT(in)    :: Ss(:,:,:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER  :: Error
    INTEGER  :: iD, iY, lo, hi
    INTEGER  :: i
    INTEGER  :: SizeDs, SizeTs, SizeYs
    INTEGER  :: i_a, i_b, i_c
    REAL(dp) :: T_a, T_b, T_c
    REAL(dp) :: S_a, S_b, S_c
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: T_i, S_i
    REAL(dp) :: Ds_i(1:2), Ts_i(1:3), Ys_i(1:2)
    REAL(dp) :: Ts_a(1:2), Ts_b(1:2), Ts_c(1:2)
    REAL(dp), DIMENSION(1:2,1:2,1:2) :: Ss_a, Ss_b, Ss_c
    REAL(dp), DIMENSION(1:2,1:3,1:2) :: Ss_i

    ! --- Initial Error Check -------------------------------------------

    Error = 0

    IF( .NOT. InversionInitialized )THEN
      Error = 10
    END IF

    IF( D < MinD .OR. D > MaxD )THEN
      Error = 01
    END IF

    IF( S < MinS .OR. S > MaxS )THEN
      Error = 02
    END IF

    IF( Y < MinY .OR. Y > MaxY )THEN
      Error = 03
    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D( D, Ds, SizeDs )
    iY = Index1D( Y, Ys, SizeYs )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iY = MIN( MAX( 1, iY ), SizeYs - 1 )

    Ds_i(1:2) = Ds(iD:iD+1)
    Ys_i(1:2) = Ys(iY:iY+1)

    i_a = 1
    T_a = Ts(i_a)
    Ts_a(1:2) = Ts(i_a:i_a+1)
    Ss_a(1:2,1:2,1:2) = Ss(iD:iD+1,i_a:i_a+1,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_a, Y, Ds_i, Ts_a, Ys_i, OS, Ss_a, S_a )

    f_a = S - S_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Ts_b(1:2) = Ts(i_b-1:i_b)
    Ss_b(1:2,1:2,1:2) = Ss(iD:iD+1,i_b-1:i_b,iY:iY+1)

    CALL LogInterpolateSingleVariable &
           ( D, T_b, Y, Ds_i, Ts_b, Ys_i, OS, Ss_b, S_b )

    f_b = S - S_b

    IF( f_a * f_b < 0.0_dp )THEN

      DO WHILE( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Ts_c(1:2) = Ts(i_c:i_c+1)
        Ss_c(1:2,1:2,1:2) = Ss(iD:iD+1,i_c:i_c+1,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_c, Y, Ds_i, Ts_c, Ys_i, OS, Ss_c, S_c )

        f_c = S - S_c

        IF( f_a * f_c < 0.0_dp )THEN

          i_b = i_c
          T_b = T_c
          S_b = S_c
          f_b = f_c

        ELSE

          i_a = i_c
          T_a = T_c
          S_a = S_c
          f_a = f_c

        END IF

      END DO

    ELSE

      ! --- Pick Highest Temperature Root ---
      DO i = SizeTs - 1, 2, -1

        lo = i-1
        hi = i+1
        T_i = Ts(i)
        Ts_i(1:3) = Ts(lo:hi)
        Ss_i(1:2,1:3,1:2) = Ss(iD:iD+1,lo:hi,iY:iY+1)

        CALL LogInterpolateSingleVariable &
               ( D, T_i, Y, Ds_i, Ts_i, Ys_i, OS, Ss_i, S_i )

        f_a = S - S_i
        f_b = S - S_b

        IF( f_a * f_b < 0.0_dp )THEN
          i_a = i
          T_a = T_i
          S_a = S_i
          EXIT
        ELSE
          i_b = i
          T_b = T_i
          S_b = S_i
        END IF

      END DO

      f_a = S - S_a
      f_b = S - S_b

      IF( f_a * f_b >= 0.0_dp )THEN

        Error = 13

      END IF

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, S_a, S_b, S, OS )
    END IF

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Single_NoGuess


  REAL(dp) FUNCTION InverseLogInterp( x_a, x_b, y_a, y_b, y, OS )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x_a, x_b, y_a, y_b, y, OS

    InverseLogInterp &
      = 10.d0**( LOG10( x_a ) + LOG10( x_b/x_a ) &
                 * LOG10( (y+OS)/(y_a+OS) ) / LOG10( (y_b+OS)/(y_a+OS) ) )

    RETURN
  END FUNCTION InverseLogInterp


  SUBROUTINE DescribeEOSInversionError( Error )

    INTEGER, INTENT(in) :: Error

    CHARACTER(64) :: ErrorString(00:13)

    ErrorString(00) = 'Returned Successfully'
    ErrorString(01) = 'First Argument (D) Outside Table Bounds'
    ErrorString(02) = 'Second Argument (E, P, or S) Outside Table Bounds'
    ErrorString(03) = 'Third Argument (Y) Outside Table Bounds'
    ErrorString(10) = 'EOS Inversion Not Initialized'
    ErrorString(13) = 'Unable to Find Any Root'

    WRITE(*,*)
    WRITE(*,*) '  wlEOSInversionModule ERROR: ' // TRIM( ErrorString(Error) )
    WRITE(*,*)

  END SUBROUTINE DescribeEOSInversionError


  SUBROUTINE WriteVector( N, Vec, FileName )

    INTEGER,          INTENT(in) :: N
    REAL(dp),         INTENT(in) :: Vec(N)
    CHARACTER(LEN=*), INTENT(in) :: FileName

    INTEGER :: FUNIT

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) N, Vec(1:N)

    CLOSE( FUNIT )

  END SUBROUTINE WriteVector


END MODULE wlEOSInversionModule

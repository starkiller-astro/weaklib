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
  PUBLIC :: ComputeTemperatureWith_DPY
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
    MODULE PROCEDURE ComputeTemperatureWith_DEY_Single
  END INTERFACE ComputeTemperatureWith_DEY

  INTERFACE ComputeTemperatureWith_DPY
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Many
    MODULE PROCEDURE ComputeTemperatureWith_DPY_Single
  END INTERFACE ComputeTemperatureWith_DPY

  INTERFACE ComputeTemperatureWith_DSY
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Many
    MODULE PROCEDURE ComputeTemperatureWith_DSY_Single
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

    LOGICAL :: UseInitialGuess
    INTEGER :: Error(SIZE(D))
    INTEGER :: i

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

    Error = 0

    DO i = 1, SIZE( D )

      CALL ComputeTemperatureWith_DEY_Single &
             ( D(i), E(i), Y(i), Ds, Ts, Ys, Es, OS, T(i), &
               UseInitialGuess_Option = UseInitialGuess, &
               Error_Option = Error(i) )

    END DO

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Many


  SUBROUTINE ComputeTemperatureWith_DEY_Single &
    ( D, E, Y, Ds, Ts, Ys, Es, OS, T, UseInitialGuess_Option, Error_Option )
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
    REAL(dp), INTENT(inout) :: T
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    LOGICAL  :: UseInitialGuess
    LOGICAL  :: LocalRoot(SIZE(Ts))
    INTEGER  :: Error
    INTEGER  :: iD, iY, i, d_i, lo, hi, SizeTs
    INTEGER  :: i_a, i_b, i_c
    REAL(dp) :: E_a, E_b, E_c
    REAL(dp) :: f_a, f_b, f_c, T_0
    REAL(dp) :: EvsT(SIZE(Ts))

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

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
      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF
      RETURN
    END IF

    iD = MIN( MAX( 1, Index1D( D, Ds, SIZE( Ds ) ) ), SIZE( Ds ) - 1 )
    iY = MIN( MAX( 1, Index1D( Y, Ys, SIZE( Ys ) ) ), SIZE( Ys ) - 1 )

    ! -------------------------------------------------------------------

    IF( UseInitialGuess )THEN

      ! --- First Check if Initial Guess Gives a Solution ---

      i_a = MIN( MAX( 1, Index1D( T, Ts, SIZE( Ts ) ) ), SIZE( Ts ) - 1 )

      CALL LogInterpolateSingleVariable &
             ( D, Ts(i_a), Y, Ds(iD:iD+1), Ts(i_a:i_a+1), Ys(iY:iY+1), &
               OS, Es(iD:iD+1,i_a:i_a+1,iY:iY+1), E_a )

      f_a = E - E_a

      i_b = i_a + 1

      CALL LogInterpolateSingleVariable &
             ( D, Ts(i_b), Y, Ds(iD:iD+1), Ts(i_b-1:i_b), Ys(iY:iY+1), &
               OS, Es(iD:iD+1,i_b-1:i_b,iY:iY+1), E_b )

      f_b = E - E_b

      IF( f_a * f_b < 0.0_dp )THEN
        T = InverseLogInterp( Ts(i_a), Ts(i_b), E_a, E_b, E, OS )
        IF( PRESENT( Error_Option ) )THEN
          Error_Option = Error
        END IF
        RETURN
      END IF

    END IF

    ! -------------------------------------------------------------------

    i_a = 1

    CALL LogInterpolateSingleVariable &
           ( D, Ts(i_a), Y, Ds(iD:iD+1), Ts(i_a:i_a+1), Ys(iY:iY+1), &
             OS, Es(iD:iD+1,i_a:i_a+1,iY:iY+1), E_a )

    f_a = E - E_a

    i_b = SIZE( Ts )

    CALL LogInterpolateSingleVariable &
           ( D, Ts(i_b), Y, Ds(iD:iD+1), Ts(i_b-1:i_b), Ys(iY:iY+1), &
             OS, Es(iD:iD+1,i_b-1:i_b,iY:iY+1), E_b )

    f_b = E - E_b

    IF( f_a * f_b < 0.0_dp )THEN

      DO WHILE( i_b > i_a + 1 )

        IF( UseInitialGuess )THEN
          i_c = MIN( MAX( i_a+1, Index1D( T, Ts, SIZE( Ts ) ) ), i_b-1 )
          UseInitialGuess = .FALSE.
        ELSE
          i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        END IF

        CALL LogInterpolateSingleVariable &
               ( D, Ts(i_c), Y, Ds(iD:iD+1), Ts(i_c:i_c+1), Ys(iY:iY+1), &
                 OS, Es(iD:iD+1,i_c:i_c+1,iY:iY+1), E_c )

        f_c = E - E_c

        IF( f_a * f_c < 0.0_dp )THEN

          i_b = i_c
          E_b = E_c
          f_b = f_c

        ELSE

          i_a = i_c
          E_a = E_c
          f_a = f_c

        END IF

      END DO

    ELSE

      SizeTs = SIZE( Ts )

      DO i = 1, SizeTs
        lo = MAX( i-1, 1 )
        hi = MIN( i+1, SizeTs )
        CALL LogInterpolateSingleVariable &
               ( D, Ts(i), Y, Ds(iD:iD+1), Ts(lo:hi), Ys(iY:iY+1), &
                 OS, Es(iD:iD+1,lo:hi,iY:iY+1), EvsT(i) )
      END DO

      LocalRoot = .FALSE.
      DO i = 1, SizeTs - 1
        IF( (E-EvsT(i)) * (E-EvsT(i+1)) < 0.0_dp )THEN
          LocalRoot(i) = .TRUE.
        END IF
      END DO

      IF( ANY( LocalRoot ) )THEN
        ! --- Pick Highest Temperature Root ---
        DO i = 1, SizeTs - 1
          IF( LocalRoot(i) )THEN
            i_a = i
            i_b = i + 1
          END IF
        END DO
        IF( UseInitialGuess )THEN
          ! --- Pick Root Nearest Initial Guess ---
          i_c = MIN( MAX( 1, Index1D( T, Ts, SizeTs ) ), SizeTs-1 )
          d_i = SizeTs
          DO i = 1, SizeTs - 1
            IF( LocalRoot(i) )THEN
              IF( ABS(i-i_c) < d_i )THEN
                i_a = i
                i_b = i + 1
                d_i = ABS(i-i_c)
              END IF
            END IF
          END DO
        END IF
        E_a = EvsT(i_a)
        E_b = EvsT(i_b)
      ELSE
        T = 0.0_dp
        Error = 13
        IF( PRESENT( Error_Option ) )THEN
          Error_Option = Error
        END IF
        RETURN
      END IF

    END IF

    T = InverseLogInterp( Ts(i_a), Ts(i_b), E_a, E_b, E, OS )

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEY_Single


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

    LOGICAL :: UseInitialGuess
    INTEGER :: Error(SIZE(D))
    INTEGER :: i

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

    Error = 0

    DO i = 1, SIZE( D )

      CALL ComputeTemperatureWith_DPY_Single &
             ( D(i), P(i), Y(i), Ds, Ts, Ys, Ps, OS, T(i), &
               UseInitialGuess_Option = UseInitialGuess, &
               Error_Option = Error(i) )

    END DO

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Many


  SUBROUTINE ComputeTemperatureWith_DPY_Single &
    ( D, P, Y, Ds, Ts, Ys, Ps, OS, T, UseInitialGuess_Option, Error_Option )
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
    REAL(dp), INTENT(inout) :: T
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    LOGICAL :: UseInitialGuess
    INTEGER :: Error

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

    Error = 13

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPY_Single


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

    LOGICAL :: UseInitialGuess
    INTEGER :: Error(SIZE(D))
    INTEGER :: i

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

    Error = 0

    DO i = 1, SIZE( D )

      CALL ComputeTemperatureWith_DSY_Single &
             ( D(i), S(i), Y(i), Ds, Ts, Ys, Ss, OS, T(i), &
               UseInitialGuess_Option = UseInitialGuess )

    END DO

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Many


  SUBROUTINE ComputeTemperatureWith_DSY_Single &
    ( D, S, Y, Ds, Ts, Ys, Ss, OS, T, UseInitialGuess_Option, Error_Option )
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
    REAL(dp), INTENT(inout) :: T
    LOGICAL,  INTENT(in),  OPTIONAL :: UseInitialGuess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    LOGICAL :: UseInitialGuess
    INTEGER :: Error

    IF( PRESENT( UseInitialGuess_Option ) )THEN
      UseInitialGuess = UseInitialGuess_Option
    ELSE
      UseInitialGuess = .FALSE.
    END IF

    Error = 13

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSY_Single


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

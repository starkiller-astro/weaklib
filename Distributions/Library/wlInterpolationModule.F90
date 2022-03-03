MODULE wlInterpolationModule

  USE wlKindModule, ONLY: dp
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    GetIndexAndDelta_Lin, &
    GetIndexAndDelta_Log, &
    LinearInterp1D_1DArray_Point, &
    LinearInterp2D_2DArray_Point, &
    LinearInterp3D_3DArray_Point, &
    LinearInterp4D_4DArray_Point, &
    LinearInterp5D_5DArray_Point, &
    LinearInterp2D_3DArray_1DAligned_Point, &
    LinearInterp3D_4DArray_1DAligned_Point, &
    LinearInterp4D_5DArray_1DAligned_Point, &
    LinearInterp2D_4DArray_2DAligned_Point, &
    LinearInterp3D_5DArray_2DAligned_Point, &
    LinearInterpDeriv3D_3DArray_Point, &
    LinearInterpDeriv4D_4DArray_Point, &
    LinearInterpDeriv2D_4DArray_2DAligned_Point

#if defined(WEAKLIB_OACC)
  USE openacc, ONLY: acc_async_sync
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: LogInterpolateSingleVariable
  PUBLIC :: LogInterpolateDifferentiateSingleVariable

  PUBLIC :: LogInterpolateSingleVariable_2D_Custom
  PUBLIC :: LogInterpolateSingleVariable_2D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_3D_Custom
  PUBLIC :: LogInterpolateSingleVariable_3D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_4D_Custom
  PUBLIC :: LogInterpolateSingleVariable_4D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_1D3D_Custom
  PUBLIC :: LogInterpolateSingleVariable_1D3D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom_Aligned
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_3D_Custom
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_3D_Custom_Point
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_2D2D_Custom
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Point
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P

  PUBLIC :: SumLogInterpolateSingleVariable_2D2D_Custom_Aligned

  REAL(dp), PARAMETER :: One = 1.0_dp
  REAL(dp), PARAMETER :: ln10 = LOG(10.d0)

  INTERFACE LogInterpolateSingleVariable
    MODULE PROCEDURE LogInterpolateSingleVariable_2D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_2D_Custom_Point
    MODULE PROCEDURE LogInterpolateSingleVariable_3D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_3D_Custom_Point
    MODULE PROCEDURE LogInterpolateSingleVariable_4D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_4D_Custom_Point
    MODULE PROCEDURE LogInterpolateSingleVariable_1D3D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_1D3D_Custom_Point
    MODULE PROCEDURE LogInterpolateSingleVariable_2D2D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_2D2D_Custom_Point
    MODULE PROCEDURE LogInterpolateSingleVariable_2D2D_Custom_Aligned
    MODULE PROCEDURE LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point
  END INTERFACE LogInterpolateSingleVariable

  INTERFACE LogInterpolateDifferentiateSingleVariable
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_3D_Custom
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_3D_Custom_Point
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_2D2D_Custom
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Point
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P
  END INTERFACE LogInterpolateDifferentiateSingleVariable

CONTAINS


  SUBROUTINE LogInterpolateSingleVariable_2D_Custom &
    ( X, Y, Xs, Ys, OS, Table, Interpolant, Error_Option, GPU_Option )

    REAL(dp), INTENT(in)  :: X (:), Y (:)
    REAL(dp), INTENT(in)  :: Xs(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option
    LOGICAL,  INTENT(in),  OPTIONAL :: GPU_Option

    INTEGER :: iP, Error
    LOGICAL :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

    Error = 0
    IF( .NOT. SIZE(X) == SIZE(Y) )THEN
      Error = 1
      IF( PRESENT( Error_Option ) ) Error_Option = Error
      RETURN
    END IF

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD IF( do_gpu )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR IF( do_gpu ) &
    !$ACC PRESENT( X, Xs, Y, Ys, OS, Table, Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, SIZE( X )

      CALL LogInterpolateSingleVariable_2D_Custom_Point &
             ( X(iP), Y(iP), Xs, Ys, OS, Table, Interpolant(iP) )

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_2D_Custom


  SUBROUTINE LogInterpolateSingleVariable_2D_Custom_Point &
    ( X, Y, Xs, Ys, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: X, Y
    REAL(dp), INTENT(in)  :: Xs(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:)
    REAL(dp), INTENT(out) :: Interpolant

    INTEGER  :: iX, iY
    REAL(dp) :: dX, dY

    CALL GetIndexAndDelta_Lin( X, Xs, iX, dX )
    CALL GetIndexAndDelta_Lin( Y, Ys, iY, dY )

    CALL LinearInterp2D_2DArray_Point &
           ( iX, iY, dX, dY, OS, Table, Interpolant )

  END SUBROUTINE LogInterpolateSingleVariable_2D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom &
    ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, OS, Table, Interpolant, GPU_Option )

    REAL(dp), INTENT(in)  :: LogE (:), LogD (:), LogT (:), Y (:)
    REAL(dp), INTENT(in)  :: LogEs(:), LogDs(:), LogTs(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option

    INTEGER  :: i, j
    INTEGER  :: iD, iT, iY, iE
    REAL(dp) :: dD, dT, dY, dE
    LOGICAL  :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) IF( do_gpu ) &
    !$OMP PRIVATE( iE, iD, dD, iT, dE, dT, iY, dY )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF( do_gpu ) &
    !$ACC PRIVATE( iE, iD, dD, iT, dE, dT, iY, dY ) &
    !$ACC PRESENT( LogE, LogEs, LogD, LogDs, LogT, LogTs, Y, Ys, OS, Table, Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iE, iD, dD, iT, dE, dT, iY, dY )
#endif
    DO j = 1, SIZE( LogD )
      DO i = 1, SIZE( LogE )

        CALL GetIndexAndDelta_Lin( LogE(i), LogEs, iE, dE )
        CALL GetIndexAndDelta_Lin( LogD(j), LogDs, iD, dD )
        CALL GetIndexAndDelta_Lin( LogT(j), LogTs, iT, dT )
        CALL GetIndexAndDelta_Lin(    Y(j),    Ys, iY, dY )

        CALL LinearInterp4D_4DArray_Point &
               ( iE, iD, iT, iY, dE, dD, dT, dY, OS, Table, Interpolant(i,j) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom


  SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom_Point &
    ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogE (:), LogD    , LogT    , Y
    REAL(dp), INTENT(in)  :: LogEs(:), LogDs(:), LogTs(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)

    INTEGER  :: i
    INTEGER  :: iD, iT, iY, iE
    REAL(dp) :: dD, dT, dY, dE

    CALL GetIndexAndDelta_Lin( LogD, LogDs, iD, dD )
    CALL GetIndexAndDelta_Lin( LogT, LogTs, iT, dT )
    CALL GetIndexAndDelta_Lin(    Y,    Ys, iY, dY )

    DO i = 1, SIZE( LogE )

      CALL GetIndexAndDelta_Lin( LogE(i), LogEs, iE, dE )

      CALL LinearInterp4D_4DArray_Point &
             ( iE, iD, iT, iY, dE, dD, dT, dY, OS, Table, Interpolant(i) )

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogE (:), LogT (:), LogX (:)
    REAL(dp), INTENT(in)  :: LogEs(:), LogTs(:), LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, l, ij, i0, j0, SizeE
    INTEGER  :: iT, iX, iE1, iE2
    REAL(dp) :: dT, dX, dE1, dE2
    LOGICAL  :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

    IF( PRESENT( ASYNC_Option ) )THEN
      async_flag = ASYNC_Option
    ELSE
#if defined(WEAKLIB_OMP_OL)
      async_flag = -1
#elif defined(WEAKLIB_OACC)
      async_flag = acc_async_sync
#endif
    END IF

    SizeE = SIZE( LogE )

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) IF( do_gpu ) &
    !$OMP PRIVATE( iE1, dE1, iE2, dE2, iT, dT, iX, dX, &
    !$OMP          i0, j0, i, j )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF( do_gpu ) ASYNC( async_flag ) &
    !$ACC PRIVATE( iE1, dE1, iE2, dE2, iT, dT, iX, dX, &
    !$ACC          i0, j0, i, j ), &
    !$ACC PRESENT( LogE, LogEs, LogT, LogTs, LogX, LogXs, OS, Table, Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iE1, dE1, iE2, dE2, iT, dT, iX, dX, &
    !$OMP          i0, j0, i, j )
#endif
    DO l = 1, SIZE( LogT )
      DO ij = 1, SizeE*(SizeE+1)/2 ! Collapsed triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        CALL GetIndexAndDelta_Lin( LogE(i), LogEs, iE1, dE1 )
        CALL GetIndexAndDelta_Lin( LogE(j), LogEs, iE2, dE2 )
        CALL GetIndexAndDelta_Lin( LogT(l), LogTs, iT, dT )
        CALL GetIndexAndDelta_Lin( LogX(l), LogXs, iX, dX )

        CALL LinearInterp4D_4DArray_Point &
               ( iE1, iE2, iT, iX, dE1, dE2, dT, dX, OS, Table, Interpolant(i,j,l) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Point &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogE (:), LogT    , LogX
    REAL(dp), INTENT(in)  :: LogEs(:), LogTs(:), LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)

    INTEGER  :: i, j, ij, i0, j0, SizeE
    INTEGER  :: iT, iX, iE1, iE2
    REAL(dp) :: dT, dX, dE1, dE2

    SizeE = SIZE( LogE )

    CALL GetIndexAndDelta_Lin( LogT, LogTs, iT, dT )
    CALL GetIndexAndDelta_Lin( LogX, LogXs, iX, dX )

    DO ij = 1, SizeE*(SizeE+1)/2 ! Collapsed triangular loop: DO j = 1, SizeE ; DO i = 1, j
      j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
      i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
      IF ( i0 > j0 ) THEN
        j = SizeE - j0 + 1
        i = SizeE - i0 + 2
      ELSE
        j = j0
        i = i0
      END IF

      CALL GetIndexAndDelta_Lin( LogE(i), LogEs, iE1, dE1 )
      CALL GetIndexAndDelta_Lin( LogE(j), LogEs, iE2, dE2 )

      CALL LinearInterp4D_4DArray_Point &
             ( iE1, iE2, iT, iX, dE1, dE2, dT, dX, OS, Table, Interpolant(i,j) )

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Aligned &
    ( LogT, LogX, LogTs, LogXs, OS, Table, Interpolant, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogT (:), LogX (:)
    REAL(dp), INTENT(in)  :: LogTs(:), LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, k, ij, i0, j0, SizeE
    INTEGER  :: iT, iX
    REAL(dp) :: dT, dX
    LOGICAL  :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

    IF( PRESENT( ASYNC_Option ) )THEN
      async_flag = ASYNC_Option
    ELSE
#if defined(WEAKLIB_OMP_OL)
      async_flag = -1
#elif defined(WEAKLIB_OACC)
      async_flag = acc_async_sync
#endif
    END IF

    SizeE = SIZE( Interpolant, DIM = 1 )

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE IF( do_gpu ) &
    !$OMP PRIVATE( iT, dT, iX, dX )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG IF( do_gpu ) ASYNC( async_flag ) &
    !$ACC PRIVATE( iT, dT, iX, dX ) &
    !$ACC PRESENT( LogT, LogTs, LogX, LogXs, OS, Table, Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iT, dT, iX, dX, &
    !$OMP          i0, j0, i, j )
#endif
    DO k = 1, SIZE( LogT )

      !CALL GetIndexAndDelta_Lin( LogT(k), LogTs, iT, dT )
      !CALL GetIndexAndDelta_Lin( LogX(k), LogXs, iX, dX )
      iT = Index1D_Lin( LogT(k), LogTs )
      iX = Index1D_Lin( LogX(k), LogXs )
      dT = ( LogT(k) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
      dX = ( LogX(k) - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( i0, j0, i, j )
#elif defined(WEAKLIB_OACC)
      !$ACC LOOP VECTOR &
      !$ACC PRIVATE( i0, j0, i, j )
#endif
      DO ij = 1, SizeE*(SizeE+1)/2 ! Collapsed triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        CALL LinearInterp2D_4DArray_2DAligned_Point &
               ( i, j, iT, iX, dT, dX, OS, Table, Interpolant(i,j,k) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Aligned


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
    ( LogT, LogX, LogTs, LogXs, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogT    , LogX
    REAL(dp), INTENT(in)  :: LogTs(:), LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)

    INTEGER  :: i, j
    INTEGER  :: iT, iX
    REAL(dp) :: dT, dX

    CALL GetIndexAndDelta_Lin( LogT, LogTs, iT, dT )
    CALL GetIndexAndDelta_Lin( LogX, LogXs, iX, dX )

    DO j = 1, SIZE( Interpolant, DIM = 1 )
      DO i = 1, j

        CALL LinearInterp2D_4DArray_2DAligned_Point &
               ( i, j, iT, iX, dT, dX, OS, Table, Interpolant(i,j) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point


  SUBROUTINE SumLogInterpolateSingleVariable_2D2D_Custom_Aligned &
    ( LogD, LogT, LogDs, LogTs, Alpha, OS, Table, Interpolant, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogD (:,:), LogT (:)
    REAL(dp), INTENT(in)  :: LogDs(:)  , LogTs(:)
    REAL(dp), INTENT(in)  :: Alpha(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, k, l, ij, i0, j0, SizeE
    INTEGER  :: iD(SIZE(Alpha)), iT
    REAL(dp) :: dD(SIZE(Alpha)), dT
    REAL(dp) :: Interp, SumInterp
    LOGICAL  :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

    IF( PRESENT( ASYNC_Option ) )THEN
      async_flag = ASYNC_Option
    ELSE
#if defined(WEAKLIB_OMP_OL)
      async_flag = -1
#elif defined(WEAKLIB_OACC)
      async_flag = acc_async_sync
#endif
    END IF

    SizeE = SIZE( Interpolant, DIM = 1 )

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE IF( do_gpu ) &
    !$OMP PRIVATE( iT, dT, iD, dD )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG IF( do_gpu ) ASYNC( async_flag ) &
    !$ACC PRIVATE( iT, dT, iD, dD ) &
    !$ACC PRESENT( LogT, LogTs, LogD, LogDs, OS, Table, Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iT, dT, iD, dD, Interp, SumInterp, &
    !$OMP          i0, j0, i, j )
#endif
    DO k = 1, SIZE( LogT )

      !CALL GetIndexAndDelta_Lin( LogT(k), LogTs, iT, dT )
      !CALL GetIndexAndDelta_Lin( LogD(k), LogDs, iD, dD )
      iT = Index1D_Lin( LogT(k), LogTs )
      dT = ( LogT(k) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
      DO l = 1, SIZE( Alpha )
        iD(l) = Index1D_Lin( LogD(l,k), LogDs )
        dD(l) = ( LogD(l,k) - LogDs(iD(l)) ) / ( LogDs(iD(l)+1) - LogDs(iD(l)) )
      END DO

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( i0, j0, i, j, Interp, SumInterp )
#elif defined(WEAKLIB_OACC)
      !$ACC LOOP VECTOR &
      !$ACC PRIVATE( i0, j0, i, j, Interp, SumInterp )
#endif
      DO ij = 1, SizeE*(SizeE+1)/2 ! Collapsed triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        SumInterp = Zero
        DO l = 1, SIZE( Alpha )
          CALL LinearInterp2D_4DArray_2DAligned_Point &
                 ( i, j, iD(l), iT, dD(l), dT, OS, Table, Interp )
          SumInterp = SumInterp + Alpha(l) * Interp
        END DO
        Interpolant(i,j,k) = SumInterp

      END DO
    END DO

  END SUBROUTINE SumLogInterpolateSingleVariable_2D2D_Custom_Aligned


  SUBROUTINE LogInterpolateSingleVariable_3D_Custom &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant, Error_Option, GPU_Option )

    REAL(dp), INTENT(in)  :: D (:), T (:), Y (:)
    REAL(dp), INTENT(in)  :: Ds(:), Ts(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option
    LOGICAL,  INTENT(in),  OPTIONAL :: GPU_Option

    INTEGER :: iP, Error
    LOGICAL :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

    Error = 0
    IF( .NOT. ALL( [ SIZE(T), SIZE(Y) ] == SIZE(D) ) )THEN
      Error = 1
      IF( PRESENT( Error_Option ) ) Error_Option = Error
      RETURN
    END IF

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD IF( do_gpu )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR IF( do_gpu ) &
    !$ACC PRESENT( D, Ds, T, Ts, Y, Ys, OS, Table, Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, SIZE( D )

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D(iP), T(iP), Y(iP), Ds, Ts, Ys, OS, Table, Interpolant(iP) )

    END DO

    IF( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE LogInterpolateSingleVariable_3D_Custom


  SUBROUTINE LogInterpolateSingleVariable_3D_Custom_Point &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D    , T    , Y
    REAL(dp), INTENT(in)  :: Ds(:), Ts(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:)
    REAL(dp), INTENT(out) :: Interpolant

    INTEGER  :: iD, iT, iY
    REAL(dp) :: dD, dT, dY

    CALL GetIndexAndDelta_Log( D, Ds, iD, dD )
    CALL GetIndexAndDelta_Log( T, Ts, iT, dT )
    CALL GetIndexAndDelta_Lin( Y, Ys, iY, dY )

    CALL LinearInterp3D_3DArray_Point &
           ( iD, iT, iY, dD, dT, dY, OS, Table, Interpolant )

  END SUBROUTINE LogInterpolateSingleVariable_3D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_4D_Custom &
      ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, &
        OS, Table, Interpolant, Error_Option, GPU_Option )

    REAL(dp), INTENT(in)  :: LogE (:), LogD (:), LogT (:),  Y(:)
    REAL(dp), INTENT(in)  :: LogEs(:), LogDs(:), LogTs(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option
    LOGICAL,  INTENT(in),  OPTIONAL :: GPU_Option

    INTEGER :: iP, Error
    LOGICAL :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

    Error = 0
    IF( .NOT. ALL( [ SIZE(LogD), SIZE(LogT), SIZE(Y) ] == SIZE(LogE) ) )THEN
      Error = 1
      IF( PRESENT( Error_Option ) ) Error_Option = Error
      RETURN
    END IF

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD IF( do_gpu )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR IF( do_gpu ) &
    !$ACC PRESENT( LogE, LogEs, LogD, LogDs, LogT, LogTs, Y, Ys, OS, Table, Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, SIZE( LogE )

      CALL LogInterpolateSingleVariable_4D_Custom_Point &
             ( LogE(iP), LogD(iP), LogT(iP), Y(iP), &
               LogEs, LogDs, LogTs, Ys, OS, Table, Interpolant(iP) )

    END DO

    IF( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE LogInterpolateSingleVariable_4D_Custom


  SUBROUTINE LogInterpolateSingleVariable_4D_Custom_Point &
    ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogE    , LogD    , LogT    , Y
    REAL(dp), INTENT(in)  :: LogEs(:), LogDs(:), LogTs(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant

    INTEGER  :: iD, iT, iY, iE
    REAL(dp) :: dD, dT, dY, dE

    CALL GetIndexAndDelta_Lin( LogE, LogEs, iE, dE )
    CALL GetIndexAndDelta_Lin( LogD, LogDs, iD, dD )
    CALL GetIndexAndDelta_Lin( LogT, LogTs, iT, dT )
    CALL GetIndexAndDelta_Lin(    Y,    Ys, iY, dY )

    CALL LinearInterp4D_4DArray_Point &
           ( iE, iD, iT, iY, dE, dD, dT, dY, OS, Table, Interpolant )

  END SUBROUTINE LogInterpolateSingleVariable_4D_Custom_Point


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D_Custom &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant, Derivative )

    REAL(dp), INTENT(in)  :: D (:), T (:), Y (:)
    REAL(dp), INTENT(in)  :: Ds(:), Ts(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)
    REAL(dp), INTENT(out) :: Derivative(:,:)

    INTEGER  :: iD, iT, iY, iP
    REAL(dp) :: dD, dT, dY
    REAL(dp) :: aD, aT, aY

    DO iP = 1, SIZE( D )

      CALL GetIndexAndDelta_Log( D(iP), Ds, iD, dD )
      CALL GetIndexAndDelta_Log( T(iP), Ts, iT, dT )
      CALL GetIndexAndDelta_Lin( Y(iP), Ys, iY, dY )
      aD = One / ( D(iP) * LOG10( Ds(iD+1) / Ds(iD) ) )
      aT = One / ( T(iP) * LOG10( Ts(iT+1) / Ts(iT) ) )
      aY = ln10 / ( Ys(iY+1) - Ys(iY) )

      CALL LinearInterpDeriv3D_3DArray_Point &
             ( iD, iT, iY, dD, dT, dY, aD, aT, aY, OS, Table, Interpolant(iP), &
               Derivative(iP,1), Derivative(iP,2), Derivative(iP,3) )

    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D_Custom


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant, Derivative )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D    , T    , Y
    REAL(dp), INTENT(in)  :: Ds(:), Ts(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:)
    REAL(dp), INTENT(out) :: Interpolant
    REAL(dp), INTENT(out) :: Derivative(:)

    INTEGER  :: iD, iT, iY
    REAL(dp) :: dD, dT, dY
    REAL(dp) :: aD, aT, aY

    CALL GetIndexAndDelta_Log( D, Ds, iD, dD )
    CALL GetIndexAndDelta_Log( T, Ts, iT, dT )
    CALL GetIndexAndDelta_Lin( Y, Ys, iY, dY )
    aD = One / ( D * LOG10( Ds(iD+1) / Ds(iD) ) )
    aT = One / ( T * LOG10( Ts(iT+1) / Ts(iT) ) )
    aY = ln10 / ( Ys(iY+1) - Ys(iY) )

    CALL LinearInterpDeriv3D_3DArray_Point &
           ( iD, iT, iY, dD, dT, dY, aD, aT, aY, OS, Table, Interpolant, &
             Derivative(1), Derivative(2), Derivative(3) )

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D_Custom_Point


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant, &
      DerivativeT, DerivativeX, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogE (:), LogT (:), LogX (:)
    REAL(dp), INTENT(in)  :: LogEs(:), LogTs(:), LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:)
    REAL(dp), INTENT(out) :: DerivativeT(:,:,:)
    REAL(dp), INTENT(out) :: DerivativeX(:,:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, l, ij, i0, j0, SizeE
    INTEGER  :: iT, iX, iE1, iE2
    REAL(dp) :: dT, dX, dE1, dE2
    REAL(dp) :: aT, aX, dI1, dI2
    LOGICAL  :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

    IF( PRESENT( ASYNC_Option ) )THEN
      async_flag = ASYNC_Option
    ELSE
#if defined(WEAKLIB_OMP_OL)
      async_flag = -1
#elif defined(WEAKLIB_OACC)
      async_flag = acc_async_sync
#endif
    END IF

    SizeE = SIZE( Interpolant, DIM = 1 )

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) IF( do_gpu ) &
    !$OMP PRIVATE( iE1, dE1, iE2, dE2, iT, dT, aT, iX, dX, aX, dI1, dI2, &
    !$OMP          i0, j0, i, j )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF( do_gpu ) ASYNC( async_flag ) &
    !$ACC PRIVATE( iE1, dE1, iE2, dE2, iT, dT, aT, iX, dX, aX, dI1, dI2, &
    !$ACC          i0, j0, i, j ), &
    !$ACC PRESENT( LogE, LogEs, LogT, LogTs, LogX, LogXs, OS, &
    !$ACC          Table, Interpolant, DerivativeT, DerivativeX )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iE1, dE1, iE2, dE2, iT, dT, aT, iX, dX, aX, dI1, dI2, &
    !$OMP          i0, j0, i, j )
#endif
    DO l = 1, SIZE( LogT )
      DO ij = 1, SizeE*(SizeE+1)/2 ! Collapsed triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        CALL GetIndexAndDelta_Lin( LogE(i), LogEs, iE1, dE1 )
        CALL GetIndexAndDelta_Lin( LogE(j), LogEs, iE2, dE2 )
        CALL GetIndexAndDelta_Lin( LogT(l), LogTs, iT , dT  )
        CALL GetIndexAndDelta_Lin( LogX(l), LogXs, iX , dX  )
        aT = One / ( LogTs(iT+1) - LogTs(iT) ) / 10.0d0**LogT(l)
        aX = One / ( LogXs(iX+1) - LogXs(iX) ) / 10.0d0**LogX(l)

        CALL LinearInterpDeriv4D_4DArray_Point &
               ( iE1, iE2, iT, iX, dE1, dE2, dT, dX, One, One, aT, aX, &
                 OS, Table, Interpolant(i,j,l), dI1, dI2, DerivativeT(i,j,l), DerivativeX(i,j,l) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Point &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant, &
      DerivativeT, DerivativeX )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogE (:), LogT    , LogX
    REAL(dp), INTENT(in)  :: LogEs(:), LogTs(:), LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)
    REAL(dp), INTENT(out) :: DerivativeT(:,:)
    REAL(dp), INTENT(out) :: DerivativeX(:,:)

    INTEGER  :: i, j
    INTEGER  :: iT, iX, iE1, iE2
    REAL(dp) :: dT, dX, dE1, dE2
    REAL(dp) :: aT, aX, dI1, dI2

    CALL GetIndexAndDelta_Lin( LogT, LogTs, iT, dT )
    CALL GetIndexAndDelta_Lin( LogX, LogXs, iX, dX )
    aT = One / ( LogTs(iT+1) - LogTs(iT) ) / 10.0d0**LogT
    aX = One / ( LogXs(iX+1) - LogXs(iX) ) / 10.0d0**LogX

    DO j = 1, SIZE( LogE )
      CALL GetIndexAndDelta_Lin( LogE(j), LogEs, iE2, dE2 )
      DO i = 1, j
        CALL GetIndexAndDelta_Lin( LogE(i), LogEs, iE1, dE1 )

        CALL LinearInterpDeriv4D_4DArray_Point &
               ( iE1, iE2, iT, iX, dE1, dE2, dT, dX, One, One, aT, aX, &
                 OS, Table, Interpolant(i,j), dI1, dI2, DerivativeT(i,j), DerivativeX(i,j) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Point


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned &
    ( LogT, LogX, LogTs, LogXs, OS, Table, Interpolant, &
      DerivativeT, DerivativeX, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogT (:), LogX (:)
    REAL(dp), INTENT(in)  :: LogTs(:), LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:)
    REAL(dp), INTENT(out) :: DerivativeT(:,:,:)
    REAL(dp), INTENT(out) :: DerivativeX(:,:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, k, ij, i0, j0, SizeE
    INTEGER  :: iT, iX
    REAL(dp) :: dT, aT, dX, aX
    LOGICAL  :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

    IF( PRESENT( ASYNC_Option ) )THEN
      async_flag = ASYNC_Option
    ELSE
#if defined(WEAKLIB_OMP_OL)
      async_flag = -1
#elif defined(WEAKLIB_OACC)
      async_flag = acc_async_sync
#endif
    END IF

    SizeE = SIZE( Interpolant, DIM = 1 )

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE IF( do_gpu ) &
    !$OMP PRIVATE( iT, dT, aT, iX, dX, aX )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG IF( do_gpu ) ASYNC( async_flag ) &
    !$ACC PRIVATE( iT, dT, aT, iX, dX, aX ) &
    !$ACC PRESENT( LogT, LogTs, LogX, LogXs, OS, Table, &
    !$ACC          DerivativeT, DerivativeX, Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iT, dT, aT, iX, dX, aX, &
    !$OMP          i0, j0, i, j )
#endif
    DO k = 1, SIZE( LogT )

      !CALL GetIndexAndDelta_Lin( LogT(k), LogTs, iT, dT )
      !CALL GetIndexAndDelta_Lin( LogX(k), LogXs, iX, dX )
      iT = Index1D_Lin( LogT(k), LogTs )
      iX = Index1D_Lin( LogX(k), LogXs )
      dT = ( LogT(k) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
      dX = ( LogX(k) - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )
      aT = One / ( LogTs(iT+1) - LogTs(iT) ) / 10.0d0**LogT(k)
      aX = One / ( LogXs(iX+1) - LogXs(iX) ) / 10.0d0**LogX(k)

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( i0, j0, i, j )
#elif defined(WEAKLIB_OACC)
      !$ACC LOOP VECTOR &
      !$ACC PRIVATE( i0, j0, i, j )
#endif
      DO ij = 1, SizeE*(SizeE+1)/2 ! Collapsed triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        CALL LinearInterpDeriv2D_4DArray_2DAligned_Point &
               ( i, j, iT, iX, dT, dX, aT, aX, &
                 OS, Table, Interpolant(i,j,k), DerivativeT(i,j,k), DerivativeX(i,j,k) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P &
    ( LogT, LogX, LogTs, LogXs, OS, Table, Interpolant, &
      DerivativeT, DerivativeX )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogT    , LogX
    REAL(dp), INTENT(in)  :: LogTs(:), LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)
    REAL(dp), INTENT(out) :: DerivativeT(:,:)
    REAL(dp), INTENT(out) :: DerivativeX(:,:)

    INTEGER  :: i, j
    INTEGER  :: iT, iX
    REAL(dp) :: dT, aT, dX, aX

    CALL GetIndexAndDelta_Lin( LogT, LogTs, iT, dT )
    CALL GetIndexAndDelta_Lin( LogX, LogXs, iX, dX )
    aT = One / ( LogTs(iT+1) - LogTs(iT) ) / 10.0d0**LogT
    aX = One / ( LogXs(iX+1) - LogXs(iX) ) / 10.0d0**LogX

    DO j = 1, SIZE( Interpolant, DIM = 1 )
      DO i = 1, j

        CALL LinearInterpDeriv2D_4DArray_2DAligned_Point &
               ( i, j, iT, iX, dT, dX, aT, aX, &
                 OS, Table, Interpolant(i,j), DerivativeT(i,j), DerivativeX(i,j) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P


END MODULE wlInterpolationModule

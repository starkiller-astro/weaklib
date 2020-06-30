MODULE wlInterpolationModule

  USE wlKindModule, ONLY: dp
#if defined(WEAKLIB_OACC)
  USE openacc, ONLY: acc_async_sync
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: locate
  PUBLIC :: Index1D
  PUBLIC :: Index1D_Lin
  PUBLIC :: TriLinear
  PUBLIC :: LogInterpolateSingleVariable
  PUBLIC :: LogInterpolateDifferentiateSingleVariable
  PUBLIC :: LogInterpolateSingleVariable_2D_Custom
  PUBLIC :: LogInterpolateSingleVariable_2D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_3D_Custom
  PUBLIC :: LogInterpolateSingleVariable_3D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_1D3D
  PUBLIC :: LogInterpolateSingleVariable_1D3D_Custom
  PUBLIC :: LogInterpolateSingleVariable_1D3D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_2D2D
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom_Point
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom_Aligned
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_2D2D_Custom
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Point
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P
  PUBLIC :: LogInterpolateOpacity_2D1D2D
  PUBLIC :: LogInterpolateOpacity_2D1D2D_Custom

  REAL(dp), PARAMETER :: One = 1.0_dp
  REAL(dp), PARAMETER :: ln10 = LOG(10.d0)

  INTERFACE LogInterpolateSingleVariable
    MODULE PROCEDURE LogInterpolateSingleVariable_2D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_2D_Custom_Point
    MODULE PROCEDURE LogInterpolateSingleVariable_3D
    MODULE PROCEDURE LogInterpolateSingleVariable_3D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_3D_Custom_Point
    MODULE PROCEDURE LogInterpolateSingleVariable_1D3D
    MODULE PROCEDURE LogInterpolateSingleVariable_1D3D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_1D3D_Custom_Point
    MODULE PROCEDURE LogInterpolateSingleVariable_4D_Custom_Point
  END INTERFACE LogInterpolateSingleVariable

  INTERFACE LogInterpolateDifferentiateSingleVariable
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_3D
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_3D_Custom
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_3D_Custom_Point
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_4D
  END INTERFACE LogInterpolateDifferentiateSingleVariable

CONTAINS


  SUBROUTINE locate( xx, n, x, j )

    INTEGER,  INTENT(in)  :: n
    INTEGER,  INTENT(out) :: j
    REAL(dp), INTENT(in)  :: x,xx(n)

    INTEGER :: jl,jm,ju

    jl = 0
    ju = n+1
    DO WHILE ( ju - jl > 1 )
      jm = (ju+jl)/2
      IF ((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) THEN
        jl = jm
      ELSE
        ju = jm
      END IF
    END DO

    IF (x.eq.xx(1)) THEN
      j = 1
    ELSEIF (x.eq.xx(n)) THEN
      j = n-1  
    ELSE
      j = jl
    END IF

  END SUBROUTINE locate


  INTEGER FUNCTION Index1D( x, xx, n )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    INTEGER :: il, im, iu

    il = 0
    iu = n+1
    DO WHILE ( iu - il > 1 )
      im = (iu+il)/2
      IF ((xx(n).ge.xx(1)).eqv.(x.ge.xx(im))) THEN
        il = im
      ELSE
        iu = im
      END IF
    END DO

    IF (x.eq.xx(1)) THEN
      Index1D = 1
    ELSEIF (x.eq.xx(n)) THEN
      Index1D = n-1
    ELSE
      Index1D = il
    END IF

    RETURN
  END FUNCTION Index1D


  FUNCTION Index1D_Lin( x, xx, n ) &
      RESULT( Index1D )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n
    INTEGER :: Index1D

    Index1D &
      = FLOOR( 1 + (n-1)*(x-xx(1))/(xx(n)-xx(1)) + 1.d-12 )

    RETURN
  END FUNCTION Index1D_Lin


  INTEGER FUNCTION Index1D_Log( x, xx, n )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    Index1D_Log &
      = FLOOR( 1 + (n-1)*LOG10(x/xx(1))/LOG10(xx(n)/xx(1)) + 1.d-12 )

    RETURN
  END FUNCTION Index1D_Log


  REAL(dp) FUNCTION Linear &
    ( p0, p1, dX1 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p0, p1, dX1

    REAL(dp) :: ddX1

    ddX1 = One - dX1

    Linear =  ddX1 * p0 &
             + dX1 * p1 

    RETURN
  END FUNCTION Linear


  REAL(dp) FUNCTION BiLinear &
    ( p00, p10, p01, p11, dX1, dX2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p00, p10, p01, p11, dX1, dX2

    REAL(dp) :: ddX1, ddX2

    ddX1 = One - dX1
    ddX2 = One - dX2

    BiLinear =  ddX1 * ( ddX2 * p00 + dX2 * p01 ) &
               + dX1 * ( ddX2 * p10 + dX2 * p11 )

    RETURN
  END FUNCTION BiLinear


  REAL(dp) FUNCTION TriLinear &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2, dX3 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX1, dX2, dX3

    REAL(dp) :: ddX1, ddX2, ddX3

    ddX1 = One - dX1
    ddX2 = One - dX2
    ddX3 = One - dX3

    TriLinear                                        &
      = ddX3                                         &
         * (   ddX2 * ( ddX1 * p000 + dX1 * p100 )   &
             +  dX2 * ( ddX1 * p010 + dX1 * p110 ) ) &
      +  dX3                                         &
         * (   ddX2 * ( ddX1 * p001 + dX1 * p101 )   &
             +  dX2 * ( ddX1 * p011 + dX1 * p111 ) )

    RETURN
  END FUNCTION TriLinear


  REAL(dp) FUNCTION dTriLineardX1 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX2, dX3 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX2, dX3

    REAL(dp) :: ddX2, ddX3

    ddX2 = One - dX2
    ddX3 = One - dX3

    dTriLineardX1 &
      = ddX3 &
          * ( - ddX2 * p000 + ddX2 * p100   &
              -  dX2 * p010 +  dX2 * p110 ) &
       + dX3 &
          * ( - ddX2 * p001 + ddX2 * p101   &
              -  dX2 * p011 +  dX2 * p111 )

    RETURN
  END FUNCTION dTriLineardX1


  REAL(dp) FUNCTION dTriLineardX2 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX3 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX1, dX3

    REAL(dp) :: ddX1, ddX3

    ddX1 = One - dX1
    ddX3 = One - dX3

    dTriLineardX2 &
      =  ddX3 &
           * ( - ddX1 * p000 - dX1 * p100   &
               + ddX1 * p010 + dX1 * p110 ) &
        + dX3 &
           * ( - ddX1 * p001 - dX1 * p101   &
               + ddX1 * p011 + dX1 * p111 )

    RETURN
  END FUNCTION dTriLineardX2


  REAL(dp) FUNCTION dTriLineardX3 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX1, dX2

    REAL(dp) :: ddX1, ddX2

    ddX1 = One - dX1
    ddX2 = One - dX2

    dTriLineardX3 &
      = - ddX1 * ddX2 * p000 - dX1 * ddX2 * p100 &
        - ddX1 *  dX2 * p010 - dX1 *  dX2 * p110 &
        + ddX1 * ddX2 * p001 + dX1 * ddX2 * p101 &
        + ddX1 *  dX2 * p011 + dX1 *  dX2 * p111

    RETURN
  END FUNCTION dTriLineardX3


  REAL(dp) FUNCTION TetraLinear &
    ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3, dX4 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3, dX4

    REAL(dp) :: ddX1, ddX2, ddX3, ddX4

    ddX1 = 1.0_dp - dX1
    ddX2 = 1.0_dp - dX2
    ddX3 = 1.0_dp - dX3
    ddX4 = 1.0_dp - dX4

    TetraLinear                                                    &
      = ddX4                                                       &
         * (   ddX3 * (  ddX2 * ( ddX1 * p0000 + dX1 * p1000 )     &
                        + dX2 * ( ddX1 * p0100 + dX1 * p1100 ) )   &
             +  dX3 * (  ddX2 * ( ddX1 * p0010 + dX1 * p1010 )     &
                        + dX2 * ( ddX1 * p0110 + dX1 * p1110 ) ) ) &
      +  dX4                                                       &
         * (   ddX3 * (  ddX2 * ( ddX1 * p0001 + dX1 * p1001 )     &
                        + dX2 * ( ddX1 * p0101 + dX1 * p1101 ) )   &
             +  dX3 * (  ddX2 * ( ddX1 * p0011 + dX1 * p1011 )     &
                        + dX2 * ( ddX1 * p0111 + dX1 * p1111 ) ) )

    RETURN
  END FUNCTION TetraLinear


  SUBROUTINE LogInterpolateSingleVariable_2D_Custom &
    ( X, Y, Xs, Ys, OS, Table, Interpolant )

    REAL(dp), INTENT(in)  :: X (:), Y (:)
    REAL(dp), INTENT(in)  :: Xs(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)

    INTEGER  :: iP

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

    INTEGER  :: &
      iX, iY
    REAL(dp) :: &
      dX, dY, &
      p00, p10, p01, p11

    iX = Index1D( X, Xs, SIZE( Xs ) )
    iY = Index1D( Y, Ys, SIZE( Ys ) )

    dX = LOG10( X / Xs(iX) ) / LOG10( Xs(iX+1) / Xs(iX) )
    dY = LOG10( Y / Ys(iY) ) / LOG10( Ys(iY+1) / Ys(iY) )

    p00 = Table( iX  , iY   )
    P10 = Table( iX+1, iY   )
    p01 = Table( iX  , iY+1 )
    p11 = Table( iX+1, iY+1 )

    Interpolant &
      = 10.0d0**( BiLinear( p00, p10, p01, p11, dX, dY ) ) - OS

  END SUBROUTINE LogInterpolateSingleVariable_2D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_3D &
               ( x1, x2, x3, Coordinate1, Coordinate2, Coordinate3, &
                 LogInterp, Offset, Table, Interpolant, MaskVar )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp 
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    REAL(dp), INTENT(in)                :: Offset
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: MaskVar
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: work_mask
    INTEGER                             :: Masksize
    REAL(dp), DIMENSION(:), INTENT(out) :: Interpolant 

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: delta
    INTEGER :: i, j, k, il1, il2, il3
  
    epsilon = 1.d-200

    ! Check the usage of x2
    ! x2 is used to size the run of the do loop because the input is expected 
    ! to be a series of 3-tuple rho, T, Ye points; 
    !     so SIZE(x1) = SIZE(x2) = SIZE(x3) 
    ! If this subroutine is being used to make an array, make sure the other 
    ! inputs are the same size (i.e. x1(1) = x(2) = x(3) etc.)
    !-------------------------------

    Masksize = SIZE( x2 )
    ALLOCATE( work_mask( Masksize ) )

    IF ( PRESENT(MaskVar) ) THEN
      work_mask = MaskVar
    ELSE
      work_mask = .true.
    END IF

    DO i = 1, Masksize
      
      IF ( .not.work_mask(i) ) CYCLE

      CALL locate( Coordinate1, SIZE( Coordinate1 ), x1(i), il1 )
      CALL locate( Coordinate2, SIZE( Coordinate2 ), x2(i), il2 )
      CALL locate( Coordinate3, SIZE( Coordinate3 ), x3(i), il3 )

      p000 = ( Table( il1  , il2  , il3   ) )
      p100 = ( Table( il1+1, il2  , il3   ) )
      p010 = ( Table( il1  , il2+1, il3   ) )
      p110 = ( Table( il1+1, il2+1, il3   ) )
      p001 = ( Table( il1  , il2  , il3+1 ) )
      p101 = ( Table( il1+1, il2  , il3+1 ) )
      p011 = ( Table( il1  , il2+1, il3+1 ) )
      p111 = ( Table( il1+1, il2+1, il3+1 ) )

      IF ( LogInterp(1) == 1 ) THEN 
        delta(1) = LOG10( x1(i) / Coordinate1(il1) ) &
                     / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        delta(1) = ( x1(i) - Coordinate1(il1) ) &
                     / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN 
        delta(2) = LOG10( x2(i) / Coordinate2(il2) ) &
                     / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        delta(2) = ( x2(i) - Coordinate2(il2) ) &
                     / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN 
        delta(3) = LOG10( x3(i) / Coordinate3(il3) ) &
                     / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        delta(3) = ( x3(i) - Coordinate3(il3) ) &
                     / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      Interpolant(i) &
        = 10.d0**( &
            TriLinear &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, delta(1), delta(2), delta(3) ) ) - Offset
   
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_3D

  
  SUBROUTINE LogInterpolateSingleVariable_1D3D &
    ( x1, x2, x3, x4, Coordinate1, Coordinate2, Coordinate3, Coordinate4, &
      LogInterp, Offset, Table, Interpolant )

    REAL(dp), INTENT(in)  :: x1(:)
    REAL(dp), INTENT(in)  :: x2(:)
    REAL(dp), INTENT(in)  :: x3(:)
    REAL(dp), INTENT(in)  :: x4(:)
    REAL(dp), INTENT(in)  :: Coordinate1(:)
    REAL(dp), INTENT(in)  :: Coordinate2(:)
    REAL(dp), INTENT(in)  :: Coordinate3(:)
    REAL(dp), INTENT(in)  :: Coordinate4(:)
    INTEGER,  INTENT(in)  :: LogInterp(4)
    REAL(dp), INTENT(in)  :: Offset
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)

    INTEGER  :: &
      i, j, il1, il2, il3, il4
    REAL(dp) :: &
      alpha(4), delta(4)
    REAL(dp) :: &
      p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111,&
      p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    DO j = 1, SIZE( x2 )

      il2 = Index1D( x2(j), Coordinate2, SIZE( Coordinate2 ) )
      il3 = Index1D( x3(j), Coordinate3, SIZE( Coordinate3 ) )
      il4 = Index1D( x4(j), Coordinate4, SIZE( Coordinate4 ) )

      IF ( LogInterp(2) == 1 ) THEN
        delta(2) = LOG10( x2(j) / Coordinate2(il2) ) &
                     / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        delta(2) = ( x2(j) - Coordinate2(il2) ) &
                     / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        delta(3) = LOG10( x3(j) / Coordinate3(il3) ) &
                     / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        delta(3) = ( x3(j) - Coordinate3(il3) ) &
                     / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      IF ( LogInterp(4) == 1 ) THEN
        delta(4) = LOG10( x4(j) / Coordinate4(il4) ) &
                     / LOG10( Coordinate4(il4+1) / Coordinate4(il4) )
      ELSE
        delta(4) = ( x4(j) - Coordinate4(il4) ) &
                     / ( Coordinate4(il4+1) - Coordinate4(il4) )
      END IF

      DO i = 1, SIZE( x1 )

        il1 = Index1D( x1(i), Coordinate1, SIZE( Coordinate1 ) )

        IF ( LogInterp(1) == 1 ) THEN
          delta(1) = LOG10( x1(i) / Coordinate1(il1) ) &
                      / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
        ELSE
          delta(1) = ( x1(i) - Coordinate1(il1) ) &
                       / ( Coordinate1(il1+1) - Coordinate1(il1) )
        END IF

        p0000 = ( Table( il1  , il2  , il3  , il4   ) )
        p0001 = ( Table( il1  , il2  , il3  , il4+1 ) )
        p0010 = ( Table( il1  , il2  , il3+1, il4   ) )
        p0011 = ( Table( il1  , il2  , il3+1, il4+1 ) )
        p0100 = ( Table( il1  , il2+1, il3  , il4   ) )
        p0101 = ( Table( il1  , il2+1, il3  , il4+1 ) )
        p0110 = ( Table( il1  , il2+1, il3+1, il4   ) )
        p0111 = ( Table( il1  , il2+1, il3+1, il4+1 ) )
        p1000 = ( Table( il1+1, il2  , il3  , il4   ) )
        p1001 = ( Table( il1+1, il2  , il3  , il4+1 ) )
        p1010 = ( Table( il1+1, il2  , il3+1, il4   ) )
        p1011 = ( Table( il1+1, il2  , il3+1, il4+1 ) )
        p1100 = ( Table( il1+1, il2+1, il3  , il4   ) )
        p1101 = ( Table( il1+1, il2+1, il3  , il4+1 ) )
        p1110 = ( Table( il1+1, il2+1, il3+1, il4   ) )
        p1111 = ( Table( il1+1, il2+1, il3+1, il4+1 ) )

        Interpolant(i,j) &
          = TetraLinear &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                delta(1), delta(2), delta(3), delta(4) )

      END DO ! i

    END DO ! j

    DO j = 1, SIZE( x2 )
      DO i = 1, SIZE( X1 )

        Interpolant(i,j) &
          = 10.0_DP**( Interpolant(i,j) ) - Offset

      END DO
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_1D3D


  SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom &
    ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, OS, Table, Interpolant, GPU_Option )

    REAL(dp), INTENT(in)  :: LogE(:)
    REAL(dp), INTENT(in)  :: LogD(:)
    REAL(dp), INTENT(in)  :: LogT(:)
    REAL(dp), INTENT(in)  :: Y(:)
    REAL(dp), INTENT(in)  :: LogEs(:)
    REAL(dp), INTENT(in)  :: LogDs(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option

    INTEGER  :: i, j
    INTEGER  :: iD, iT, iY, iE(SIZE(LogE))
    REAL(dp) :: dD, dT, dY, dE(SIZE(LogE))
    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111
    LOGICAL  :: do_gpu

    IF( PRESENT( GPU_Option ) )THEN
      do_gpu = GPU_Option
    ELSE
      do_gpu = .FALSE.
    END IF

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: iE, dE )
#elif defined(WEAKLIB_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( iE, dE )
#endif

    ! --- Precompute Indices and Differences ---

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( iE, dE, LogE, LogEs )
#endif
    DO i = 1, SIZE( LogE )
      iE(i) = Index1D_Lin( LogE(i), LogEs, SIZE( LogEs ) )
      dE(i) = ( LogE(i) - LogEs(iE(i)) ) / ( LogEs(iE(i)+1) - LogEs(iE(i)) )
    END DO

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( iD, dD, iT, dT, iY, dY )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( iD, dD, iT, dT, iY, dY ) &
    !$ACC PRESENT( iE, dE, LogE, LogEs, LogD, LogDs, LogT, LogTs, Y, Ys, OS, Table, Interpolant )
#endif
    DO j = 1, SIZE( LogD )

      iY = Index1D_Lin( Y(j), Ys, SIZE( Ys ) )
      dY = ( Y(j) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

      iT = Index1D_Lin( LogT(j), LogTs, SIZE( LogTs ) )
      dT = ( LogT(j) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )

      iD = Index1D_Lin( LogD(j), LogDs, SIZE( LogDs ) )
      dD = ( LogD(j) - LogDs(iD) ) / ( LogDs(iD+1) - LogDs(iD) )

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      !$OMP          p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111 )
#elif defined(WEAKLIB_OACC)
      !$ACC LOOP VECTOR &
      !$ACC PRIVATE( p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      !$ACC          p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111 )
#endif
      DO i = 1, SIZE( LogE )

        p0000 = TABLE( iE(i)  , iD  , iT  , iY   )
        p1000 = TABLE( iE(i)+1, iD  , iT  , iY   )
        p0100 = TABLE( iE(i)  , iD+1, iT  , iY   )
        p1100 = TABLE( iE(i)+1, iD+1, iT  , iY   )
        p0010 = TABLE( iE(i)  , iD  , iT+1, iY   )
        p1010 = TABLE( iE(i)+1, iD  , iT+1, iY   )
        p0110 = TABLE( iE(i)  , iD+1, iT+1, iY   )
        p1110 = TABLE( iE(i)+1, iD+1, iT+1, iY   )
        p0001 = TABLE( iE(i)  , iD  , iT  , iY+1 )
        p1001 = TABLE( iE(i)+1, iD  , iT  , iY+1 )
        p0101 = TABLE( iE(i)  , iD+1, iT  , iY+1 )
        p1101 = TABLE( iE(i)+1, iD+1, iT  , iY+1 )
        p0011 = TABLE( iE(i)  , iD  , iT+1, iY+1 )
        p1011 = TABLE( iE(i)+1, iD  , iT+1, iY+1 )
        p0111 = TABLE( iE(i)  , iD+1, iT+1, iY+1 )
        p1111 = TABLE( iE(i)+1, iD+1, iT+1, iY+1 )

        Interpolant(i,j) &
          = 10.0d0**( &
              TetraLinear &
                ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                  p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                  dE(i), dD, dT, dY ) ) - OS

      END DO ! i

    END DO ! j

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: iE, dE )
#elif defined(WEAKLIB_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( iE, dE )
#endif

  END SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom


  SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom_Point &
    ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogE(:)
    REAL(dp), INTENT(in)  :: LogD
    REAL(dp), INTENT(in)  :: LogT
    REAL(dp), INTENT(in)  :: Y
    REAL(dp), INTENT(in)  :: LogEs(:)
    REAL(dp), INTENT(in)  :: LogDs(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)

    INTEGER  :: i
    INTEGER  :: iD, iT, iY, iE
    REAL(dp) :: dD, dT, dY, dE
    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    iY = Index1D_Lin( Y, Ys, SIZE( Ys ) )
    dY = ( Y - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

    iT = Index1D_Lin( LogT, LogTs, SIZE( LogTs ) )
    dT = ( LogT - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )

    iD = Index1D_Lin( LogD, LogDs, SIZE( LogDs ) )
    dD = ( LogD - LogDs(iD) ) / ( LogDs(iD+1) - LogDs(iD) )

    DO i = 1, SIZE( LogE )

      iE = Index1D_Lin( LogE(i), LogEs, SIZE( LogEs ) )
      dE = ( LogE(i) - LogEs(iE) ) / ( LogEs(iE+1) - LogEs(iE) )

      p0000 = TABLE( iE  , iD  , iT  , iY   )
      p1000 = TABLE( iE+1, iD  , iT  , iY   )
      p0100 = TABLE( iE  , iD+1, iT  , iY   )
      p1100 = TABLE( iE+1, iD+1, iT  , iY   )
      p0010 = TABLE( iE  , iD  , iT+1, iY   )
      p1010 = TABLE( iE+1, iD  , iT+1, iY   )
      p0110 = TABLE( iE  , iD+1, iT+1, iY   )
      p1110 = TABLE( iE+1, iD+1, iT+1, iY   )
      p0001 = TABLE( iE  , iD  , iT  , iY+1 )
      p1001 = TABLE( iE+1, iD  , iT  , iY+1 )
      p0101 = TABLE( iE  , iD+1, iT  , iY+1 )
      p1101 = TABLE( iE+1, iD+1, iT  , iY+1 )
      p0011 = TABLE( iE  , iD  , iT+1, iY+1 )
      p1011 = TABLE( iE+1, iD  , iT+1, iY+1 )
      p0111 = TABLE( iE  , iD+1, iT+1, iY+1 )
      p1111 = TABLE( iE+1, iD+1, iT+1, iY+1 )

      Interpolant(i) &
        = 10.0d0**( &
            TetraLinear &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dE, dD, dT, dY ) ) - OS

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_2D2D &
               ( x1, x2, x3, x4, Coordinate1, Coordinate2, Coordinate3, &
                 Coordinate4, LogInterp, Offset, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: x1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x4
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate4
    INTEGER,  DIMENSION(4),       INTENT(in)  :: LogInterp
    REAL(dp),                     INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:,:,:),   INTENT(out) :: Interpolant

    INTEGER :: &
      i, ip, j, k, il1, il2, il3, il4, &
      SizeC1, SizeC2, SizeC3, SizeC4, &
      SizeX1, SizeX2, SizeX3
    REAL(dp), DIMENSION(4) :: &
      delta
    REAL(dp) :: &
      p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    REAL(dp), PARAMETER :: kmev = 8.61733d-11

    SizeX1 = SIZE( x1 )  ! ep
    SizeX2 = SIZE( x2 )  ! e
    SizeX3 = SIZE( x3 )  ! T

    SizeC1 = SIZE( Coordinate1 )  ! ep
    SizeC2 = SIZE( Coordinate2 )  ! e 
    SizeC3 = SIZE( Coordinate3 )  ! T
    SizeC4 = SIZE( Coordinate4 )  ! eta

    DO k = 1, SizeX3

      IF ( LogInterp(3) == 1 ) THEN
        il3 = Index1D_Log( x3(k), Coordinate3, SizeC3 )
        delta(3) &
          = LOG10( x3(k) / Coordinate3(il3) ) &
            / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        il3 = Index1D_Lin( x3(k), Coordinate3, SizeC3 )
        delta(3) &
          = ( x3(k) - Coordinate3(il3) ) &
            / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      IF ( LogInterp(4) == 1 ) THEN
        il4 = Index1D_Log( x4(k), Coordinate4, SizeC4 ) 
        delta(4) &
          = LOG10( x4(k) / Coordinate4(il4) ) &
            / LOG10( Coordinate4(il4+1) / Coordinate4(il4) )
      ELSE
        il4 = Index1D_Lin( x4(k), Coordinate4, SizeC4 ) 
        delta(4) &
          = ( x4(k) - Coordinate4(il4) ) &
            / ( Coordinate4(il4+1) - Coordinate4(il4) )
      END IF
             
      DO j = 1, SizeX2

        IF ( LogInterp(2) == 1 ) THEN
          il2 = Index1D_Log( x2(j), Coordinate2, SizeC2 )
          delta(2) &
            = LOG10( x2(j) / Coordinate2(il2) ) &
              / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
        ELSE
          il2 = Index1D_Lin( x2(j), Coordinate2, SizeC2 )
          delta(2) &
            = ( x2(j) - Coordinate2(il2) ) &
              / ( Coordinate2(il2+1) - Coordinate2(il2) )
        END IF

        DO i = 1, j

          IF ( LogInterp(1) == 1 ) THEN
            il1 = Index1D_Log( x1(i), Coordinate1, SizeC1 )
            delta(1) &
              = LOG10( x1(i) / Coordinate1(il1) ) &
                / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
          ELSE
            il1 = Index1D_Lin( x1(i), Coordinate1, SizeC1 )
            delta(1) &
              = ( x1(i) - Coordinate1(il1) ) &
                / ( Coordinate1(il1+1) - Coordinate1(il1) )
          END IF

          p0000 = ( Table( il1  , il2  , il3  , il4   ) )
          p0001 = ( Table( il1  , il2  , il3  , il4+1 ) )
          p0010 = ( Table( il1  , il2  , il3+1, il4   ) )
          p0011 = ( Table( il1  , il2  , il3+1, il4+1 ) )
          p0100 = ( Table( il1  , il2+1, il3  , il4   ) )
          p0101 = ( Table( il1  , il2+1, il3  , il4+1 ) )
          p0110 = ( Table( il1  , il2+1, il3+1, il4   ) )
          p0111 = ( Table( il1  , il2+1, il3+1, il4+1 ) )
          p1000 = ( Table( il1+1, il2  , il3  , il4   ) )
          p1001 = ( Table( il1+1, il2  , il3  , il4+1 ) )
          p1010 = ( Table( il1+1, il2  , il3+1, il4   ) )
          p1011 = ( Table( il1+1, il2  , il3+1, il4+1 ) )
          p1100 = ( Table( il1+1, il2+1, il3  , il4   ) )
          p1101 = ( Table( il1+1, il2+1, il3  , il4+1 ) )
          p1110 = ( Table( il1+1, il2+1, il3+1, il4   ) )
          p1111 = ( Table( il1+1, il2+1, il3+1, il4+1 ) ) 

          Interpolant(i,j,k) &
            = TetraLinear &
                ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                  p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                  delta(1), delta(2), delta(3), delta(4) )

        END DO ! i
      END DO ! j
    END DO ! k
  

    DO k = 1,SizeX3
      DO j = 1, SizeX2
        DO i = 1, j

          Interpolant(i,j,k) &
            = 10.d0**( Interpolant(i,j,k) ) - Offset

        END DO ! i

!        DO ip = j, SizeX1
!
!          Interpolant(ip,j,k) &
!            = Interpolant(j,ip,k) * &
!              EXP( ( x1(j) - x2(ip) )/( kmev * x3(k) ) )
!
!        END DO ! ip

      END DO ! j
    END DO ! k


  END SUBROUTINE LogInterpolateSingleVariable_2D2D


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogE(:)
    REAL(dp), INTENT(in)  :: LogT(:)
    REAL(dp), INTENT(in)  :: LogX(:)
    REAL(dp), INTENT(in)  :: LogEs(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, l, ij, i0, j0, SizeE
    INTEGER  :: iT, iX, iE(SIZE(LogE))
    REAL(dp) :: dT, dX, dE(SIZE(LogE))
    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111
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
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: iE, dE )
#elif defined(WEAKLIB_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( iE, dE )
#endif

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC( async_flag ) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( iE, dE, LogE, LogEs )
#endif
    DO i = 1, SizeE
      iE(i) = Index1D_Lin( LogE(i), LogEs, SIZE( LogEs ) )
      dE(i) = ( LogE(i) - LogEs(iE(i)) ) / ( LogEs(iE(i)+1) - LogEs(iE(i)) )
    END DO

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( iT, dT, iX, dX )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG ASYNC( async_flag ) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( iT, dT, iX, dX ) &
    !$ACC PRESENT( iE, dE, LogE, LogEs, LogT, LogTs, LogX, LogXs, OS, Table, Interpolant )
#endif
    DO l = 1, SIZE( LogT )

      iT = Index1D_Lin( LogT(l), LogTs, SIZE( LogTs ) )
      dT = ( LogT(l) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )

      iX = Index1D_Lin( LogX(l), LogXs, SIZE( LogXs ) )
      dX = ( LogX(l) - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( i0, j0, i, j, &
      !$OMP          p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      !$OMP          p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111 )
#elif defined(WEAKLIB_OACC)
      !$ACC LOOP VECTOR &
      !$ACC PRIVATE( i0, j0, i, j, &
      !$ACC          p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      !$ACC          p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111 )
#endif
      DO ij = 1, SizeE*(SizeE+1)/2 ! Fused triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        p0000 = TABLE( iE(i)  , iE(j)  , iT  , iX   )
        p1000 = TABLE( iE(i)+1, iE(j)  , iT  , iX   )
        p0100 = TABLE( iE(i)  , iE(j)+1, iT  , iX   )
        p1100 = TABLE( iE(i)+1, iE(j)+1, iT  , iX   )
        p0010 = TABLE( iE(i)  , iE(j)  , iT+1, iX   )
        p1010 = TABLE( iE(i)+1, iE(j)  , iT+1, iX   )
        p0110 = TABLE( iE(i)  , iE(j)+1, iT+1, iX   )
        p1110 = TABLE( iE(i)+1, iE(j)+1, iT+1, iX   )
        p0001 = TABLE( iE(i)  , iE(j)  , iT  , iX+1 )
        p1001 = TABLE( iE(i)+1, iE(j)  , iT  , iX+1 )
        p0101 = TABLE( iE(i)  , iE(j)+1, iT  , iX+1 )
        p1101 = TABLE( iE(i)+1, iE(j)+1, iT  , iX+1 )
        p0011 = TABLE( iE(i)  , iE(j)  , iT+1, iX+1 )
        p1011 = TABLE( iE(i)+1, iE(j)  , iT+1, iX+1 )
        p0111 = TABLE( iE(i)  , iE(j)+1, iT+1, iX+1 )
        p1111 = TABLE( iE(i)+1, iE(j)+1, iT+1, iX+1 )

        Interpolant(i,j,l) &
          = 10.0d0**( &
              TetraLinear &
                ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                  p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                  dE(i), dE(j), dT, dX ) ) - OS

      END DO

    END DO

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: iE, dE )
#elif defined(WEAKLIB_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( iE, dE )
#endif

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Point &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogE(:)
    REAL(dp), INTENT(in)  :: LogT
    REAL(dp), INTENT(in)  :: LogX
    REAL(dp), INTENT(in)  :: LogEs(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)

    INTEGER  :: i, j
    INTEGER  :: iT, iX, iE(SIZE(LogE))
    REAL(dp) :: dT, dX, dE(SIZE(LogE))
    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    iT = Index1D_Lin( LogT, LogTs, SIZE( LogTs ) )
    dT = ( LogT - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )

    iX = Index1D_Lin( LogX, LogXs, SIZE( LogXs ) )
    dX = ( LogX - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )

    DO i = 1, SIZE( LogE )
      iE(i) = Index1D_Lin( LogE(i), LogEs, SIZE( LogEs ) )
      dE(i) = ( LogE(i) - LogEs(iE(i)) ) / ( LogEs(iE(i)+1) - LogEs(iE(i)) )
    END DO

    DO j = 1, SIZE( LogE )
      DO i = 1, j

        p0000 = TABLE( iE(i)  , iE(j)  , iT  , iX   )
        p1000 = TABLE( iE(i)+1, iE(j)  , iT  , iX   )
        p0100 = TABLE( iE(i)  , iE(j)+1, iT  , iX   )
        p1100 = TABLE( iE(i)+1, iE(j)+1, iT  , iX   )
        p0010 = TABLE( iE(i)  , iE(j)  , iT+1, iX   )
        p1010 = TABLE( iE(i)+1, iE(j)  , iT+1, iX   )
        p0110 = TABLE( iE(i)  , iE(j)+1, iT+1, iX   )
        p1110 = TABLE( iE(i)+1, iE(j)+1, iT+1, iX   )
        p0001 = TABLE( iE(i)  , iE(j)  , iT  , iX+1 )
        p1001 = TABLE( iE(i)+1, iE(j)  , iT  , iX+1 )
        p0101 = TABLE( iE(i)  , iE(j)+1, iT  , iX+1 )
        p1101 = TABLE( iE(i)+1, iE(j)+1, iT  , iX+1 )
        p0011 = TABLE( iE(i)  , iE(j)  , iT+1, iX+1 )
        p1011 = TABLE( iE(i)+1, iE(j)  , iT+1, iX+1 )
        p0111 = TABLE( iE(i)  , iE(j)+1, iT+1, iX+1 )
        p1111 = TABLE( iE(i)+1, iE(j)+1, iT+1, iX+1 )

        Interpolant(i,j) &
          = 10.0d0**( &
              TetraLinear &
                ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                  p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                  dE(i), dE(j), dT, dX ) ) - OS

      END DO
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Aligned &
    ( LogT, LogX, LogTs, LogXs, OS, Table, Interpolant, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogT(:)
    REAL(dp), INTENT(in)  :: LogX(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, k, ij, i0, j0, SizeE
    INTEGER  :: iT, iX
    REAL(dp) :: dT, dX
    REAL(dp) :: p00, p10, p01, p11
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
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( iT, dT, iX, dX )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG ASYNC( async_flag ) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( iT, dT, iX, dX ) &
    !$ACC PRESENT( LogT, LogTs, LogX, LogXs, OS, Table, Interpolant )
#endif
    DO k = 1, SIZE( LogT )

      iT = Index1D_Lin( LogT(k), LogTs, SIZE( LogTs ) )
      dT = ( LogT(k) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )

      iX = Index1D_Lin( LogX(k), LogXs, SIZE( LogXs ) )
      dX = ( LogX(k) - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( i0, j0, i, j, p00, p10, p01, p11 )
#elif defined(WEAKLIB_OACC)
      !$ACC LOOP VECTOR &
      !$ACC PRIVATE( i0, j0, i, j, p00, p10, p01, p11 )
#endif
      DO ij = 1, SizeE*(SizeE+1)/2 ! Fused triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        p00 = Table(i,j,iT  ,iX  )
        p10 = Table(i,j,iT+1,iX  )
        p01 = Table(i,j,iT  ,iX+1)
        p11 = Table(i,j,iT+1,iX+1)

        Interpolant(i,j,k) &
          = 10.0d0**( BiLinear( p00, p10, p01, p11, dT, dX ) ) - OS

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

    REAL(dp), INTENT(in)  :: LogT
    REAL(dp), INTENT(in)  :: LogX
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)

    INTEGER  :: i, j
    INTEGER  :: iT, iX
    REAL(dp) :: dT, dX
    REAL(dp) :: p00, p10, p01, p11

    iT = Index1D_Lin( LogT, LogTs, SIZE( LogTs ) )
    dT = ( LogT - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )

    iX = Index1D_Lin( LogX, LogXs, SIZE( LogXs ) )
    dX = ( LogX - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )

    DO j = 1, SIZE( Interpolant, DIM = 1 )
      DO i = 1, j

        p00 = Table(i,j,iT  ,iX  )
        p10 = Table(i,j,iT+1,iX  )
        p01 = Table(i,j,iT  ,iX+1)
        p11 = Table(i,j,iT+1,iX+1)

        Interpolant(i,j) &
          = 10.0d0**( BiLinear( p00, p10, p01, p11, dT, dX ) ) - OS

      END DO
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point


  SUBROUTINE LogInterpolateOpacity_2D1D2D &
    ( E, X1, X2, TabE, TabX1, TabX2, LogInterpE, LogInterpX, OS, &
      TabOpacity, Opacity )

    REAL(dp), INTENT(in)  :: E(:), X1(:), X2(:)
    REAL(dp), INTENT(in)  :: TabE(:), TabX1(:), TabX2(:)
    INTEGER,  INTENT(in)  :: LogInterpE
    INTEGER,  INTENT(in)  :: LogInterpX(2)
    REAL(dp), INTENT(in)  :: OS(:)
    REAL(dp), INTENT(in)  :: TabOpacity(:,:,:,:,:)
    REAL(dp), INTENT(out) :: Opacity(:,:,:,:)

    INTEGER  :: iX, iX1, iX2, l, k, kp
    REAL(dp) :: dX1, dX2
    REAL(dp) :: &
      p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111
    INTEGER,  ALLOCATABLE :: iE(:)
    REAL(dp), ALLOCATABLE :: dE(:)

    ALLOCATE( iE(SIZE(E)), dE(SIZE(E)) )

    DO k = 1, SIZE( E )

      IF( LogInterpE == 1 )THEN
        iE(k) = Index1D_Log( E(k), TabE, SIZE( TabE ) )
        dE(k) = LOG10( E(k) / TabE(iE(k)) ) &
                  / LOG10( TabE(iE(k)+1) / TabE(iE(k)) )
      ELSE
        iE(k) = Index1D_Lin( E(k), TabE, SIZE( TabE ) )
        dE(k) = ( E(k) - TabE(iE(k)) ) &
                  / ( TabE(iE(k)+1) - TabE(iE(k)) )
      END IF

    END DO

    DO iX = 1, SIZE( X1 )

      IF( LogInterpX(1) == 1 )THEN
        iX1 = Index1D_Log( X1(iX), TabX1, SIZE( TabX1 ) )
        dX1 = LOG10( X1(iX) / TabX1(iX1) ) &
                / LOG10( TabX1(iX1+1) / TabX1(iX1) )
      ELSE
        iX1 = Index1D_Lin( X1(iX), TabX1, SIZE( TabX1 ) )
        dX1 = ( X1(iX) - TabX1(iX1) ) &
                / ( TabX1(iX1+1) - TabX1(iX1) )
      END IF

      IF( LogInterpX(2) == 1 )THEN
        iX2 = Index1D_Log( X2(iX), TabX2, SIZE( TabX2 ) )
        dX2 = LOG10( X2(iX) / TabX2(iX2) ) &
                / LOG10( TabX2(iX2+1) / TabX2(iX2) )
      ELSE
        iX2 = Index1D_Lin( X2(iX), TabX2, SIZE( TabX2 ) )
        dX2 = ( X2(iX) - TabX2(iX2) ) &
                / ( TabX2(iX2+1) - TabX2(iX2) )
      END IF

      DO l  = 1, SIZE( TabOpacity, DIM = 3 )
      DO k  = 1, SIZE( E )
      DO kp = 1, k

        p0000 = TabOpacity( iE(kp)  , iE(k)  , l, iX1  , iX2   )
        p0001 = TabOpacity( iE(kp)  , iE(k)  , l, iX1  , iX2+1 )
        p0010 = TabOpacity( iE(kp)  , iE(k)  , l, iX1+1, iX2   )
        p0011 = TabOpacity( iE(kp)  , iE(k)  , l, iX1+1, iX2+1 )
        p0100 = TabOpacity( iE(kp)  , iE(k)+1, l, iX1  , iX2   )
        p0101 = TabOpacity( iE(kp)  , iE(k)+1, l, iX1  , iX2+1 )
        p0110 = TabOpacity( iE(kp)  , iE(k)+1, l, iX1+1, iX2   )
        p0111 = TabOpacity( iE(kp)  , iE(k)+1, l, iX1+1, iX2+1 )
        p1000 = TabOpacity( iE(kp)+1, iE(k)  , l, iX1  , iX2   )
        p1001 = TabOpacity( iE(kp)+1, iE(k)  , l, iX1  , iX2+1 )
        p1010 = TabOpacity( iE(kp)+1, iE(k)  , l, iX1+1, iX2   )
        p1011 = TabOpacity( iE(kp)+1, iE(k)  , l, iX1+1, iX2+1 )
        p1100 = TabOpacity( iE(kp)+1, iE(k)+1, l, iX1  , iX2   )
        p1101 = TabOpacity( iE(kp)+1, iE(k)+1, l, iX1  , iX2+1 )
        p1110 = TabOpacity( iE(kp)+1, iE(k)+1, l, iX1+1, iX2   )
        p1111 = TabOpacity( iE(kp)+1, iE(k)+1, l, iX1+1, iX2+1 )

        Opacity(kp,k,l,iX) &
          = TetraLinear &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dE(kp), dE(k), dX1, dX2 )

        Opacity(kp,k,l,iX) = 10.0d0**( Opacity(kp,k,l,iX) ) - OS(l)

      END DO
      END DO
      END DO

    END DO

    DEALLOCATE( iE, dE )

  END SUBROUTINE LogInterpolateOpacity_2D1D2D


  SUBROUTINE LogInterpolateOpacity_2D1D2D_Custom &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant )

    REAL(dp), INTENT(in)  :: LogE(:)
    REAL(dp), INTENT(in)  :: LogT(:)
    REAL(dp), INTENT(in)  :: LogX(:)
    REAL(dp), INTENT(in)  :: LogEs(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
    REAL(dp), INTENT(in)  :: OS(:)
    REAL(dp), INTENT(in)  :: Table(:,:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:,:)

    INTEGER  :: i, j, k, kp, l, m, iT, iX
    INTEGER  :: p1, p2, p3, p4
    INTEGER  :: iE(SIZE(LogE))
    REAL(dp) :: dT, dX
    REAL(dp) :: dE(SIZE(LogE))
    REAL(dp) :: p(0:1,0:1,0:1,0:1)

    DO i = 1, SIZE( LogE )
      iE(i) = Index1D_Lin( LogE(i), LogEs, SIZE( LogEs ) )
      dE(i) = ( LogE(i) - LogEs(iE(i)) ) / ( LogEs(iE(i)+1) - LogEs(iE(i)) )
    END DO

    DO m = 1, SIZE( LogT )

      iT = Index1D_Lin( LogT(m), LogTs, SIZE( LogTs ) )
      dT = ( LogT(m) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )

      iX = Index1D_Lin( LogX(m), LogXs, SIZE( LogXs ) )
      dX = ( LogX(m) - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )

      DO l = 1, SIZE( TABLE, DIM = 3 )
      DO j = 1, SIZE( LogE )
      DO i = 1, j

        kp = iE(i)
        k  = iE(j)

        DO p4 = 0, 1
        DO p3 = 0, 1
        DO p2 = 0, 1
        DO p1 = 0, 1

          p(p1,p2,p3,p4) = TABLE(kp+p1,k+p2,l,iT+p3,iX+p4)

        END DO
        END DO
        END DO
        END DO

        Interpolant(i,j,l,m) &
          = TetraLinear &
              ( p(0,0,0,0), p(1,0,0,0), p(0,1,0,0), p(1,1,0,0), &
                p(0,0,1,0), p(1,0,1,0), p(0,1,1,0), p(1,1,1,0), &
                p(0,0,0,1), p(1,0,0,1), p(0,1,0,1), p(1,1,0,1), &
                p(0,0,1,1), p(1,0,1,1), p(0,1,1,1), p(1,1,1,1), &
                dE(i), dE(j), dT, dX )

        Interpolant(i,j,l,m) = 10.0d0**( Interpolant(i,j,l,m) ) - OS(l)

      END DO
      END DO
      END DO

    END DO

  END SUBROUTINE LogInterpolateOpacity_2D1D2D_Custom


  SUBROUTINE LogInterpolateSingleVariable_3D_Custom &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant, Error_Option )

    REAL(dp), INTENT(in)  :: D (:), T (:), Y (:)
    REAL(dp), INTENT(in)  :: Ds(:), Ts(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER  :: &
      iP, iD, iT, iY, Error
    REAL(dp) :: &
      dD, dT, dY, &
      p000, p100, p010, p110, &
      p001, p101, p011, p111

    Error = 0

    IF( .NOT. ALL( [ SIZE(T), SIZE(Y) ] == SIZE(D) ) )THEN

      Error = 1

      IF( PRESENT( Error_Option ) )THEN
        Error_Option = Error
      END IF

      RETURN

    END IF

    DO iP = 1, SIZE( D )

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D(iP), T(iP), Y(iP), Ds, Ts, Ys, OS, Table, Interpolant(iP) )

    END DO

    IF( PRESENT( Error_Option ) )THEN
      Error_Option = Error
    END IF

  END SUBROUTINE LogInterpolateSingleVariable_3D_Custom


  SUBROUTINE LogInterpolateSingleVariable_3D_Custom_Point &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D,  T,  Y
    REAL(dp), INTENT(in)  :: Ds(:), Ts(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:)
    REAL(dp), INTENT(out) :: Interpolant

    INTEGER  :: &
      iD, iT, iY
    REAL(dp) :: &
      dD, dT, dY, &
      p000, p100, p010, p110, &
      p001, p101, p011, p111
    INTEGER :: &
      SizeDs, SizeTs, SizeYs

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYs = SIZE( Ys )

    iD = Index1D( D, Ds, SizeDs )
    iT = Index1D( T, Ts, SizeTs )
    iY = Index1D( Y, Ys, SizeYs )

    dD = LOG10( D / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
    dT = LOG10( T / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
    dY = ( Y - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

    p000 = Table( iD  , iT  , iY   )
    p100 = Table( iD+1, iT  , iY   )
    p010 = Table( iD  , iT+1, iY   )
    p110 = Table( iD+1, iT+1, iY   )
    p001 = Table( iD  , iT  , iY+1 )
    p101 = Table( iD+1, iT  , iY+1 )
    p011 = Table( iD  , iT+1, iY+1 )
    p111 = Table( iD+1, iT+1, iY+1 )

    Interpolant &
      = 10.0d0**( &
          TriLinear &
            ( p000, p100, p010, p110, &
              p001, p101, p011, p111, &
              dD, dT, dY ) ) - OS

  END SUBROUTINE LogInterpolateSingleVariable_3D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_4D &
               ( x1, x2, x3, x4, Coordinate1, Coordinate2, Coordinate3, &
                 Coordinate4, LogInterp, Offset, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: x1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: x4
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Coordinate4
    INTEGER,  DIMENSION(4),       INTENT(in)  :: LogInterp
    REAL(dp),                     INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),       INTENT(out) :: Interpolant

    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111,&
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111
    REAL(dp), DIMENSION(4) :: alpha, delta
    INTEGER :: i, il1, il2, il3, il4

    IF ( ( SIZE(x1) .NE. SIZE(x2) ) .OR. ( SIZE(x2) .NE. SIZE(x3) ) &
         .OR. ( SIZE(x3) .NE. SIZE(x4) ) .OR. ( SIZE(x4) .NE. SIZE (x1) ) ) &
    THEN
      WRITE(*,*) &
        'ERROR: describe arrays (of interpolation point) have diff size.'
      RETURN
    END IF
    
    DO i = 1, SIZE( x1 )

      CALL locate( Coordinate1, SIZE( Coordinate1 ), x1(i), il1 )
      CALL locate( Coordinate2, SIZE( Coordinate2 ), x2(i), il2 )
      CALL locate( Coordinate3, SIZE( Coordinate3 ), x3(i), il3 )
      CALL locate( Coordinate4, SIZE( Coordinate4 ), x4(i), il4 )

      p0000 = ( Table( il1  , il2  , il3  , il4   ) )
      p0001 = ( Table( il1  , il2  , il3  , il4+1 ) )
      p0010 = ( Table( il1  , il2  , il3+1, il4   ) )
      p0011 = ( Table( il1  , il2  , il3+1, il4+1 ) )
      p0100 = ( Table( il1  , il2+1, il3  , il4   ) )
      p0101 = ( Table( il1  , il2+1, il3  , il4+1 ) )
      p0110 = ( Table( il1  , il2+1, il3+1, il4   ) )
      p0111 = ( Table( il1  , il2+1, il3+1, il4+1 ) )
      p1000 = ( Table( il1+1, il2  , il3  , il4   ) )
      p1001 = ( Table( il1+1, il2  , il3  , il4+1 ) )
      p1010 = ( Table( il1+1, il2  , il3+1, il4   ) )
      p1011 = ( Table( il1+1, il2  , il3+1, il4+1 ) )
      p1100 = ( Table( il1+1, il2+1, il3  , il4   ) )
      p1101 = ( Table( il1+1, il2+1, il3  , il4+1 ) )
      p1110 = ( Table( il1+1, il2+1, il3+1, il4   ) )
      p1111 = ( Table( il1+1, il2+1, il3+1, il4+1 ) )

      IF ( LogInterp(1) == 1 ) THEN
        delta(1) = LOG10( x1(i) / Coordinate1(il1) ) &
                    / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        delta(1) = ( x1(i) - Coordinate1(il1) ) &
                     / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        delta(2) = LOG10( x2(i) / Coordinate2(il2) ) &
                     / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        delta(2) = ( x2(i) - Coordinate2(il2) ) &
                     / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        delta(3) = LOG10( x3(i) / Coordinate3(il3) ) &
                     / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        delta(3) = ( x3(i) - Coordinate3(il3) ) &
                     / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      IF ( LogInterp(4) == 1 ) THEN
        delta(4) = LOG10( x4(i) / Coordinate4(il4) ) &
                     / LOG10( Coordinate4(il4+1) / Coordinate4(il4) )
      ELSE
        delta(4) = ( x4(i) - Coordinate4(il4) ) &
                     / ( Coordinate4(il4+1) - Coordinate4(il4) )
      END IF

      Interpolant(i) &
      = 10.d0**( &
        TetraLinear(p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      delta(1), delta(2), delta(3), delta(4) ) ) - Offset

    END DO


  END SUBROUTINE LogInterpolateSingleVariable_4D


  SUBROUTINE LogInterpolateSingleVariable_4D_Custom &
               ( E, D, T, Y, Es, Ds, Ts, Ys, OS, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: E,  D,  T,  Y
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Es, Ds, Ts, Ys
    REAL(dp),                     INTENT(in)  :: OS
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),       INTENT(out) :: Interpolant

    INTEGER  :: &
      iP, iE, iD, iT, iY
    REAL(dp) :: &
      p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111, &
      dE, dD, dT, dY

    IF( .NOT. ALL( [ SIZE(D), SIZE(T), SIZE(Y) ] == SIZE(E) ) )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'LogInterpolateSingleVariable_4D_Custom'
      WRITE(*,'(A4,A)') &
        '', 'ERROR: arrays of interpolation points have different sizes'
      WRITE(*,*)
      RETURN
    END IF

    DO iP = 1, SIZE( E )

      iE = Index1D( E(iP), Es, SIZE( Es ) )
      iD = Index1D( D(iP), Ds, SIZE( Ds ) )
      iT = Index1D( T(iP), Ts, SIZE( Ts ) )
      iY = Index1D( Y(iP), Ys, SIZE( Ys ) )

      p0000 = Table( iE  , iD  , iT  , iY   )
      p1000 = Table( iE+1, iD  , iT  , iY   )
      p0100 = Table( iE  , iD+1, iT  , iY   )
      p1100 = Table( iE+1, iD+1, iT  , iY   )
      p0010 = Table( iE  , iD  , iT+1, iY   )
      p1010 = Table( iE+1, iD  , iT+1, iY   )
      p0110 = Table( iE  , iD+1, iT+1, iY   )
      p1110 = Table( iE+1, iD+1, iT+1, iY   )
      p0001 = Table( iE  , iD  , iT  , iY+1 )
      p1001 = Table( iE+1, iD  , iT  , iY+1 )
      p0101 = Table( iE  , iD+1, iT  , iY+1 )
      p1101 = Table( iE+1, iD+1, iT  , iY+1 )
      p0011 = Table( iE  , iD  , iT+1, iY+1 )
      p1011 = Table( iE+1, iD  , iT+1, iY+1 )
      p0111 = Table( iE  , iD+1, iT+1, iY+1 )
      p1111 = Table( iE+1, iD+1, iT+1, iY+1 )

      dE = LOG10( E(iP) / Es(iE) ) / LOG10( Es(iE+1) / Es(iE) )
      dD = LOG10( D(iP) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
      dT = LOG10( T(iP) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
      dY = ( Y(iP) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

      Interpolant(iP) &
        = 10.0d0**( &
            TetraLinear &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dE, dD, dT, dY ) ) - OS

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_4D_Custom


  SUBROUTINE LogInterpolateSingleVariable_4D_Custom_Point &
    ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogE
    REAL(dp), INTENT(in)  :: LogD
    REAL(dp), INTENT(in)  :: LogT
    REAL(dp), INTENT(in)  :: Y
    REAL(dp), INTENT(in)  :: LogEs(:)
    REAL(dp), INTENT(in)  :: LogDs(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant

    INTEGER  :: iD, iT, iY, iE
    REAL(dp) :: dD, dT, dY, dE
    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    iE = Index1D_Lin( LogE, LogEs, SIZE( LogEs ) )
    iD = Index1D_Lin( LogD, LogDs, SIZE( LogDs ) )
    iT = Index1D_Lin( LogT, LogTs, SIZE( LogTs ) )
    iY = Index1D_Lin( Y, Ys, SIZE( Ys ) )

    p0000 = TABLE( iE  , iD  , iT  , iY   )
    p1000 = TABLE( iE+1, iD  , iT  , iY   )
    p0100 = TABLE( iE  , iD+1, iT  , iY   )
    p1100 = TABLE( iE+1, iD+1, iT  , iY   )
    p0010 = TABLE( iE  , iD  , iT+1, iY   )
    p1010 = TABLE( iE+1, iD  , iT+1, iY   )
    p0110 = TABLE( iE  , iD+1, iT+1, iY   )
    p1110 = TABLE( iE+1, iD+1, iT+1, iY   )
    p0001 = TABLE( iE  , iD  , iT  , iY+1 )
    p1001 = TABLE( iE+1, iD  , iT  , iY+1 )
    p0101 = TABLE( iE  , iD+1, iT  , iY+1 )
    p1101 = TABLE( iE+1, iD+1, iT  , iY+1 )
    p0011 = TABLE( iE  , iD  , iT+1, iY+1 )
    p1011 = TABLE( iE+1, iD  , iT+1, iY+1 )
    p0111 = TABLE( iE  , iD+1, iT+1, iY+1 )
    p1111 = TABLE( iE+1, iD+1, iT+1, iY+1 )

    dE = ( LogE - LogEs(iE) ) / ( LogEs(iE+1) - LogEs(iE) )
    dD = ( LogD - LogDs(iD) ) / ( LogDs(iD+1) - LogDs(iD) )
    dT = ( LogT - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
    dY = ( Y - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

    Interpolant &
      = 10.0d0**( &
          TetraLinear &
            ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
              p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
              dE, dD, dT, dY ) ) - OS

  END SUBROUTINE LogInterpolateSingleVariable_4D_Custom_Point


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D &
    ( x1, x2, x3, Coordinate1, Coordinate2, Coordinate3, LogInterp, Offset, &
      Table, Interpolant, Derivative, MaskVar)     

    REAL(dp), DIMENSION(:), INTENT(in)     :: x1
    REAL(dp), DIMENSION(:), INTENT(in)     :: x2
    REAL(dp), DIMENSION(:), INTENT(in)     :: x3
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate3
    INTEGER, DIMENSION(3), INTENT(in)      :: LogInterp 
    REAL(dp), INTENT(in)                   :: Offset
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    REAL(dp), DIMENSION(:), INTENT(out)    :: Interpolant
    REAL(dp), DIMENSION(:,:), INTENT(out)  :: Derivative 
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: MaskVar
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: work_mask
    INTEGER                             :: Masksize
    
    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: alpha, delta
    INTEGER :: i, j, k, il1, il2, il3

    epsilon = 1.d-200

    Masksize = SIZE( x2 )
    ALLOCATE( work_mask( Masksize ) )

    IF ( PRESENT(MaskVar) ) THEN
      work_mask = MaskVar
    ELSE
      work_mask = .true.
    END IF

    DO i = 1, Masksize

      IF ( .not.work_mask(i) ) CYCLE

      CALL locate( Coordinate1, SIZE( Coordinate1 ), x1(i), il1 )
      CALL locate( Coordinate2, SIZE( Coordinate2 ), x2(i), il2 )
      CALL locate( Coordinate3, SIZE( Coordinate3 ), x3(i), il3 )

      p000 = ( Table( il1  , il2  , il3   ) )
      p100 = ( Table( il1+1, il2  , il3   ) )
      p010 = ( Table( il1  , il2+1, il3   ) )
      p110 = ( Table( il1+1, il2+1, il3   ) )
      p001 = ( Table( il1  , il2  , il3+1 ) )
      p101 = ( Table( il1+1, il2  , il3+1 ) )
      p011 = ( Table( il1  , il2+1, il3+1 ) )
      p111 = ( Table( il1+1, il2+1, il3+1 ) )

      IF ( LogInterp(1) == 1 ) THEN
        alpha(1) &
          = 1.0d0 / ( x1(i) * LOG10( Coordinate1(il1+1) / Coordinate1(il1) ) )
        delta(1) &
          = LOG10( x1(i) / Coordinate1(il1) ) &
              / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        alpha(1) &
          = ln10 / ( Coordinate1(il1+1) - Coordinate1(il1) )
        delta(1) &
          = ( x1(i) - Coordinate1(il1) ) &
              / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        alpha(2) &
          = 1.0d0 / ( x2(i) * LOG10( Coordinate2(il2+1) / Coordinate2(il2) ) )
        delta(2) &
          = LOG10( x2(i) / Coordinate2(il2) ) &
              / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        alpha(2) &
          = ln10 / ( Coordinate2(il2+1) - Coordinate2(il2) )
        delta(2) &
          = ( x2(i) - Coordinate2(il2) ) &
              / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        alpha(3) &
          = 1.0d0 / ( x3(i) * LOG10( Coordinate3(il3+1) / Coordinate3(il3) ) )
        delta(3) &
          = LOG10( x3(i) / Coordinate3(il3) ) &
              / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        alpha(3) &
          = ln10 / ( Coordinate3(il3+1) - Coordinate3(il3) )
        delta(3) &
          = ( x3(i) - Coordinate3(il3) ) &
              / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      Interpolant(i) &
        = 10.d0**( &
            (1.0_dp - delta(3)) &
              * (   (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &
                  +           delta(1)  * (1.0_dp - delta(2)) * p100   &
                  + (1.0_dp - delta(1)) *           delta(2)  * p010   &
                  +           delta(1)  *           delta(2)  * p110 ) &
            +         delta(3) &
              * (   (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                  +           delta(1)  * (1.0_dp - delta(2)) * p101   &
                  + (1.0_dp - delta(1)) *           delta(2)  * p011   &
                  +           delta(1)  *           delta(2)  * p111 ) ) &
          - Offset

      Derivative(i,1) &
        = ( (Interpolant(i) ) * alpha(1) & 
            * ( (1.0_dp - delta(3)) * ( (delta(2) - 1.0_dp) * p000   &
                                    +  ( 1.0_dp - delta(2)) * p100   &
                                    -             delta(2)  * p010   &
                                    +             delta(2)  * p110 ) &
                         + delta(3) * ( (delta(2) - 1.0_dp) * p001   &
                                    +  ( 1.0_dp - delta(2)) * p101   &
                                    -             delta(2)  * p011   &
                                    +             delta(2)  * p111 ) ) )

      Derivative(i,2) &
        = ( ( Interpolant(i) ) * alpha(2) &
              * ( (1.0_dp - delta(3) ) * ( (delta(1) - 1.0_dp) * p000   &
                                       -             delta(1)  * p100   &
                                       +  ( 1.0_dp - delta(1)) * p010   &
                                       +             delta(1)  * p110 ) &
                            + delta(3) * ( (delta(1) - 1.0_dp) * p001   &
                                       -             delta(1)  * p101   &
                                       +   (1.0_dp - delta(1)) * p011   &
                                       +             delta(1)  * p111 ) ) )

      Derivative(i,3) &
        = ( ( Interpolant(i) ) * alpha(3) &
            * ( ( (delta(1) - 1.0_dp)) * (1.0_dp - delta(2)) * p000   &
                -            delta(1)  * (1.0_dp - delta(2)) * p100   &
                - ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                -            delta(1)  *           delta(2)  * p110   &
                +  (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                +            delta(1)  *           delta(2)  * p111 ) )

    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D_Custom &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant, Derivative )

    REAL(dp), INTENT(in)  :: D (:), T (:), Y (:)
    REAL(dp), INTENT(in)  :: Ds(:), Ts(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:)
    REAL(dp), INTENT(out) :: Derivative(:,:)

    INTEGER  :: &
      iP, iD, iT, iY
    REAL(dp) :: &
      dD, dT, dY, &
      aD, aT, aY, &
      p000, p100, p010, p110, &
      p001, p101, p011, p111

    DO iP = 1, SIZE( D )

      iD = Index1D( D(iP), Ds, SIZE( Ds ) )
      iT = Index1D( T(iP), Ts, SIZE( Ts ) )
      iY = Index1D( Y(iP), Ys, SIZE( Ys ) )

      dD = LOG10( D(iP) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
      dT = LOG10( T(iP) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
      dY = ( Y(iP) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

      aD = One / ( D(iP) * LOG10( Ds(iD+1) / Ds(iD) ) )
      aT = One / ( T(iP) * LOG10( Ts(iT+1) / Ts(iT) ) )
      aY = ln10 / ( Ys(iY+1) - Ys(iY) )

      p000 = Table( iD  , iT  , iY   )
      p100 = Table( iD+1, iT  , iY   )
      p010 = Table( iD  , iT+1, iY   )
      p110 = Table( iD+1, iT+1, iY   )
      p001 = Table( iD  , iT  , iY+1 )
      p101 = Table( iD+1, iT  , iY+1 )
      p011 = Table( iD  , iT+1, iY+1 )
      p111 = Table( iD+1, iT+1, iY+1 )

      Interpolant(iP) &
        = 10.0d0**( &
            TriLinear &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, dD, dT, dY ) ) - OS

      Derivative(iP,1) &
        = Interpolant(iP) * aD &
            * dTriLineardX1 &
                ( p000, p100, p010, p110, p001, p101, p011, p111, dT, dY )

      Derivative(iP,2) &
        = Interpolant(iP) * aT &
            * dTriLineardX2 &
                ( p000, p100, p010, p110, p001, p101, p011, p111, dD, dY )

      Derivative(iP,3) &
        = Interpolant(iP) * aY &
            * dTriLineardX3 &
                ( p000, p100, p010, p110, p001, p101, p011, p111, dD, dT )

    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D_Custom


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant, Derivative )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D, T, Y
    REAL(dp), INTENT(in)  :: Ds(:), Ts(:), Ys(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:)
    REAL(dp), INTENT(out) :: Interpolant
    REAL(dp), INTENT(out) :: Derivative(:)

    INTEGER  :: &
      iD, iT, iY
    REAL(dp) :: &
      dD, dT, dY, &
      aD, aT, aY, &
      p000, p100, p010, p110, &
      p001, p101, p011, p111

    iD = Index1D( D, Ds, SIZE( Ds ) )
    iT = Index1D( T, Ts, SIZE( Ts ) )
    iY = Index1D( Y, Ys, SIZE( Ys ) )

    dD = LOG10( D / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
    dT = LOG10( T / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
    dY = ( Y - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

    aD = One / ( D * LOG10( Ds(iD+1) / Ds(iD) ) )
    aT = One / ( T * LOG10( Ts(iT+1) / Ts(iT) ) )
    aY = ln10 / ( Ys(iY+1) - Ys(iY) )

    p000 = Table( iD  , iT  , iY   )
    p100 = Table( iD+1, iT  , iY   )
    p010 = Table( iD  , iT+1, iY   )
    p110 = Table( iD+1, iT+1, iY   )
    p001 = Table( iD  , iT  , iY+1 )
    p101 = Table( iD+1, iT  , iY+1 )
    p011 = Table( iD  , iT+1, iY+1 )
    p111 = Table( iD+1, iT+1, iY+1 )

    Interpolant &
      = 10.0d0**( &
          TriLinear &
            ( p000, p100, p010, p110, &
              p001, p101, p011, p111, dD, dT, dY ) ) - OS

    Derivative(1) &
      = Interpolant * aD &
          * dTriLineardX1 &
              ( p000, p100, p010, p110, p001, p101, p011, p111, dT, dY )

    Derivative(2) &
      = Interpolant * aT &
          * dTriLineardX2 &
              ( p000, p100, p010, p110, p001, p101, p011, p111, dD, dY )

    Derivative(3) &
      = Interpolant * aY &
          * dTriLineardX3 &
              ( p000, p100, p010, p110, p001, p101, p011, p111, dD, dT )

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D_Custom_Point


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant, &
      DerivativeT, DerivativeX, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogE(:)
    REAL(dp), INTENT(in)  :: LogT(:)
    REAL(dp), INTENT(in)  :: LogX(:)
    REAL(dp), INTENT(in)  :: LogEs(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:,:)
    REAL(dp), INTENT(out) :: DerivativeT(:,:,:)
    REAL(dp), INTENT(out) :: DerivativeX(:,:,:)
    LOGICAL,  INTENT(in), OPTIONAL :: GPU_Option
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, l, ij, i0, j0, SizeE
    INTEGER  :: iT, iX, iE(SIZE(LogE))
    REAL(dp) :: dT, dX, dE(SIZE(LogE))
    REAL(dp) :: aT, aX
    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111
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
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: iE, dE )
#elif defined(WEAKLIB_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( iE, dE )
#endif

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC( async_flag ) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( iE, dE, LogE, LogEs )
#endif
    DO i = 1, SizeE
      iE(i) = Index1D_Lin( LogE(i), LogEs, SIZE( LogEs ) )
      dE(i) = ( LogE(i) - LogEs(iE(i)) ) / ( LogEs(iE(i)+1) - LogEs(iE(i)) )
    END DO

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( iT, dT, aT, iX, dX, aX )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG ASYNC( async_flag ) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( iT, dT, aT, iX, dX, aX ) &
    !$ACC PRESENT( iE, dE, LogE, LogEs, LogT, LogTs, LogX, LogXs, OS, &
    !$ACC          Table, Interpolant, DerivativeT, DerivativeX )
#endif
    DO l = 1, SIZE( LogT )

      iT = Index1D_Lin( LogT(l), LogTs, SIZE( LogTs ) )
      dT = ( LogT(l) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
      aT = One / ( LogTs(iT+1) - LogTs(iT) ) / 10.0d0**LogT(l)

      iX = Index1D_Lin( LogX(l), LogXs, SIZE( LogXs ) )
      dX = ( LogX(l) - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )
      aX = One / ( LogXs(iX+1) - LogXs(iX) ) / 10.0d0**LogX(l)

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( i0, j0, i, j, &
      !$OMP          p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      !$OMP          p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111 )
#elif defined(WEAKLIB_OACC)
      !$ACC LOOP VECTOR &
      !$ACC PRIVATE( i0, j0, i, j, &
      !$ACC          p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      !$ACC          p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111 )
#endif
      DO ij = 1, SizeE*(SizeE+1)/2 ! Fused triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        p0000 = TABLE( iE(i)  , iE(j)  , iT  , iX   )
        p1000 = TABLE( iE(i)+1, iE(j)  , iT  , iX   )
        p0100 = TABLE( iE(i)  , iE(j)+1, iT  , iX   )
        p1100 = TABLE( iE(i)+1, iE(j)+1, iT  , iX   )
        p0010 = TABLE( iE(i)  , iE(j)  , iT+1, iX   )
        p1010 = TABLE( iE(i)+1, iE(j)  , iT+1, iX   )
        p0110 = TABLE( iE(i)  , iE(j)+1, iT+1, iX   )
        p1110 = TABLE( iE(i)+1, iE(j)+1, iT+1, iX   )
        p0001 = TABLE( iE(i)  , iE(j)  , iT  , iX+1 )
        p1001 = TABLE( iE(i)+1, iE(j)  , iT  , iX+1 )
        p0101 = TABLE( iE(i)  , iE(j)+1, iT  , iX+1 )
        p1101 = TABLE( iE(i)+1, iE(j)+1, iT  , iX+1 )
        p0011 = TABLE( iE(i)  , iE(j)  , iT+1, iX+1 )
        p1011 = TABLE( iE(i)+1, iE(j)  , iT+1, iX+1 )
        p0111 = TABLE( iE(i)  , iE(j)+1, iT+1, iX+1 )
        p1111 = TABLE( iE(i)+1, iE(j)+1, iT+1, iX+1 )

        Interpolant(i,j,l) &
          = 10.0d0**( &
              TetraLinear &
                ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                  p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                  dE(i), dE(j), dT, dX ) ) - OS

        DerivativeT(i,j,l) &
          = Interpolant(i,j,l) * aT &
              * ( Trilinear( p0000, p1000, p0100, p1100, p0001, p1001, p0101, p1101, &
                             dE(i), dE(j), dX ) &
                  - Trilinear( p0010, p1010, p0110, p1110, p0011, p1011, p0111, p1111, &
                               dE(i), dE(j), dX ) )

        DerivativeX(i,j,l) &
          = Interpolant(i,j,l) * aX &
              * ( Trilinear( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                             dE(i), dE(j), dT ) &
                  - Trilinear( p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                               dE(i), dE(j), dT ) )

      END DO
    END DO

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: iE, dE )
#elif defined(WEAKLIB_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( iE, dE )
#endif

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Point &
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant, &
      DerivativeT, DerivativeX )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogE(:)
    REAL(dp), INTENT(in)  :: LogT
    REAL(dp), INTENT(in)  :: LogX
    REAL(dp), INTENT(in)  :: LogEs(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)
    REAL(dp), INTENT(out) :: DerivativeT(:,:)
    REAL(dp), INTENT(out) :: DerivativeX(:,:)

    INTEGER  :: i, j
    INTEGER  :: iT, iX, iE(SIZE(LogE))
    REAL(dp) :: dT, dX, dE(SIZE(LogE))
    REAL(dp) :: aT, aX
    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    iT = Index1D_Lin( LogT, LogTs, SIZE( LogTs ) )
    dT = ( LogT - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
    aT = One / ( LogTs(iT+1) - LogTs(iT) ) / 10.0d0**LogT

    iX = Index1D_Lin( LogX, LogXs, SIZE( LogXs ) )
    dX = ( LogX - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )
    aX = One / ( LogXs(iX+1) - LogXs(iX) ) / 10.0d0**LogX

    DO i = 1, SIZE( LogE )
      iE(i) = Index1D_Lin( LogE(i), LogEs, SIZE( LogEs ) )
      dE(i) = ( LogE(i) - LogEs(iE(i)) ) / ( LogEs(iE(i)+1) - LogEs(iE(i)) )
    END DO

    DO j = 1, SIZE( LogE )
    DO i = 1, j

      p0000 = TABLE( iE(i)  , iE(j)  , iT  , iX   )
      p1000 = TABLE( iE(i)+1, iE(j)  , iT  , iX   )
      p0100 = TABLE( iE(i)  , iE(j)+1, iT  , iX   )
      p1100 = TABLE( iE(i)+1, iE(j)+1, iT  , iX   )
      p0010 = TABLE( iE(i)  , iE(j)  , iT+1, iX   )
      p1010 = TABLE( iE(i)+1, iE(j)  , iT+1, iX   )
      p0110 = TABLE( iE(i)  , iE(j)+1, iT+1, iX   )
      p1110 = TABLE( iE(i)+1, iE(j)+1, iT+1, iX   )
      p0001 = TABLE( iE(i)  , iE(j)  , iT  , iX+1 )
      p1001 = TABLE( iE(i)+1, iE(j)  , iT  , iX+1 )
      p0101 = TABLE( iE(i)  , iE(j)+1, iT  , iX+1 )
      p1101 = TABLE( iE(i)+1, iE(j)+1, iT  , iX+1 )
      p0011 = TABLE( iE(i)  , iE(j)  , iT+1, iX+1 )
      p1011 = TABLE( iE(i)+1, iE(j)  , iT+1, iX+1 )
      p0111 = TABLE( iE(i)  , iE(j)+1, iT+1, iX+1 )
      p1111 = TABLE( iE(i)+1, iE(j)+1, iT+1, iX+1 )

      Interpolant(i,j) &
        = 10.0d0**( &
            TetraLinear &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dE(i), dE(j), dT, dX ) ) - OS

      DerivativeT(i,j) &
        = Interpolant(i,j) * aT &
            * ( Trilinear( p0010, p1010, p0110, p1110, &
                           p0011, p1011, p0111, p1111, dE(i), dE(j), dX ) &
              - Trilinear( p0000, p1000, p0100, p1100, &
                           p0001, p1001, p0101, p1101, dE(i), dE(j), dX ) )

      DerivativeX(i,j) &
        = Interpolant(i,j) * aX &
            * ( Trilinear( p0001, p1001, p0101, p1101, &
                           p0011, p1011, p0111, p1111, dE(i), dE(j), dT ) &
              - Trilinear( p0000, p1000, p0100, p1100, &
                           p0010, p1010, p0110, p1110, dE(i), dE(j), dT ) )

    END DO
    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Point


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned &
    ( LogT, LogX, LogTs, LogXs, OS, Table, Interpolant, &
      DerivativeT, DerivativeX, GPU_Option, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogT(:)
    REAL(dp), INTENT(in)  :: LogX(:)
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
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
    REAL(dp) :: p00, p10, p01, p11
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
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( iT, dT, aT, iX, dX, aX )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG ASYNC( async_flag ) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( iT, dT, aT, iX, dX, aX ) &
    !$ACC PRESENT( LogT, LogTs, LogX, LogXs, OS, Table, Interpolant )
#endif
    DO k = 1, SIZE( LogT )

      iT = Index1D_Lin( LogT(k), LogTs, SIZE( LogTs ) )
      dT = ( LogT(k) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
      aT = One / ( LogTs(iT+1) - LogTs(iT) ) / 10.0d0**LogT(k)

      iX = Index1D_Lin( LogX(k), LogXs, SIZE( LogXs ) )
      dX = ( LogX(k) - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )
      aX = One / ( LogXs(iX+1) - LogXs(iX) ) / 10.0d0**LogX(k)

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD &
      !$OMP PRIVATE( i0, j0, i, j, p00, p10, p01, p11 )
#elif defined(WEAKLIB_OACC)
      !$ACC LOOP VECTOR &
      !$ACC PRIVATE( i0, j0, i, j, p00, p10, p01, p11 )
#endif
      DO ij = 1, SizeE*(SizeE+1)/2 ! Fused triangular loop: DO j = 1, SizeE ; DO i = 1, j
        j0 = MOD( (ij-1) / ( SizeE + 1 ), SizeE + 1 ) + 1
        i0 = MOD( (ij-1)                , SizeE + 1 ) + 1
        IF ( i0 > j0 ) THEN
          j = SizeE - j0 + 1
          i = SizeE - i0 + 2
        ELSE
          j = j0
          i = i0
        END IF

        p00 = Table(i,j,iT  ,iX  )
        p10 = Table(i,j,iT+1,iX  )
        p01 = Table(i,j,iT  ,iX+1)
        p11 = Table(i,j,iT+1,iX+1)

        Interpolant(i,j,k) &
          = 10.0d0**( BiLinear( p00, p10, p01, p11, dT, dX ) ) - OS

        DerivativeT(i,j,k) &
          = Interpolant(i,j,k) * aT &
              * ( Linear( p00, p01, dX ) - Linear( p10, p11, dX ) )

        DerivativeX(i,j,k) &
          = Interpolant(i,j,k) * aX &
              * ( Linear( p00, p10, dT ) - Linear( p01, p11, dX ) )

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

    REAL(dp), INTENT(in)  :: LogT
    REAL(dp), INTENT(in)  :: LogX
    REAL(dp), INTENT(in)  :: LogTs(:)
    REAL(dp), INTENT(in)  :: LogXs(:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(:,:,:,:)
    REAL(dp), INTENT(out) :: Interpolant(:,:)
    REAL(dp), INTENT(out) :: DerivativeT(:,:)
    REAL(dp), INTENT(out) :: DerivativeX(:,:)

    INTEGER  :: i, j
    INTEGER  :: iT, iX
    REAL(dp) :: dT, aT, dX, aX
    REAL(dp) :: p00, p10, p01, p11

    iT = Index1D_Lin( LogT, LogTs, SIZE( LogTs ) )
    dT = ( LogT - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
    aT = One / ( LogTs(iT+1) - LogTs(iT) ) / 10.0d0**LogT

    iX = Index1D_Lin( LogX, LogXs, SIZE( LogXs ) )
    dX = ( LogX - LogXs(iX) ) / ( LogXs(iX+1) - LogXs(iX) )
    aX = One / ( LogXs(iX+1) - LogXs(iX) ) / 10.0d0**LogX

    DO j = 1, SIZE( Interpolant, DIM = 1 )
      DO i = 1, j

        p00 = Table(i,j,iT  ,iX  )
        p10 = Table(i,j,iT+1,iX  )
        p01 = Table(i,j,iT  ,iX+1)
        p11 = Table(i,j,iT+1,iX+1)

        Interpolant(i,j) &
          = 10.0d0**( BiLinear( p00, p10, p01, p11, dT, dX ) ) - OS

        DerivativeT(i,j) &
          = Interpolant(i,j) * aT &
              * ( Linear( p00, p01, dX ) - Linear( p10, p11, dX ) )

        DerivativeX(i,j) &
          = Interpolant(i,j) * aX &
              * ( Linear( p00, p10, dT ) - Linear( p01, p11, dX ) )

      END DO
    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D &
   ( x1, x2, x3, x4, Coordinate1, Coordinate2, Coordinate3, Coordinate4, &
     LogInterp, Offset, Table, Interpolant, Derivative, debug )     

    REAL(dp), DIMENSION(:), INTENT(in)     :: x1
    REAL(dp), DIMENSION(:), INTENT(in)     :: x2
    REAL(dp), DIMENSION(:), INTENT(in)     :: x3
    REAL(dp), DIMENSION(:), INTENT(in)     :: x4
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate3
    REAL(dp), DIMENSION(:), INTENT(in)     :: Coordinate4
    INTEGER, DIMENSION(4), INTENT(in)      :: LogInterp
    REAL(dp), INTENT(in)                   :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in) :: Table
    REAL(dp), DIMENSION(:), INTENT(out)    :: Interpolant
    REAL(dp), DIMENSION(:,:), INTENT(out)  :: Derivative
    LOGICAL                                :: debug

    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111,&
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111
    REAL(dp), DIMENSION(4) :: alpha, delta
    INTEGER :: i, j, k, l, il1, il2, il3, il4
    
    IF ( ( SIZE(x1) .NE. SIZE(x2) ) .OR. ( SIZE(x2) .NE. SIZE(x3) ) &
         .OR. ( SIZE(x3) .NE. SIZE(x4) ) .OR. ( SIZE(x4) .NE. SIZE (x1) ) ) &
    THEN
      WRITE(*,*) &
        'ERROR: describe arrays (of interpolation point) have diff size.'
      RETURN
    END IF

    IF (debug) THEN
      WRITE(*,*) ' Now we are in 4D differentiate routine'
      WRITE(*,*) ' SIZE(x1) (number of loop) is ', SIZE( x1 )
    END IF

    DO i = 1, SIZE( x1 )

      CALL locate( Coordinate1, SIZE( Coordinate1 ), x1(i), il1 )
      CALL locate( Coordinate2, SIZE( Coordinate2 ), x2(i), il2 )
      CALL locate( Coordinate3, SIZE( Coordinate3 ), x3(i), il3 )
      CALL locate( Coordinate4, SIZE( Coordinate4 ), x4(i), il4 )

      IF (debug) THEN
        WRITE(*,*) ' The located postion is ', il1, il2, il3, il4
      END IF
 
      p0000 = ( Table( il1  , il2  , il3  , il4   ) )
      p0001 = ( Table( il1  , il2  , il3  , il4+1 ) )
      p0010 = ( Table( il1  , il2  , il3+1, il4   ) )
      p0011 = ( Table( il1  , il2  , il3+1, il4+1 ) )
      p0100 = ( Table( il1  , il2+1, il3  , il4   ) )
      p0101 = ( Table( il1  , il2+1, il3  , il4+1 ) )
      p0110 = ( Table( il1  , il2+1, il3+1, il4   ) )
      p0111 = ( Table( il1  , il2+1, il3+1, il4+1 ) )
      p1000 = ( Table( il1+1, il2  , il3  , il4   ) )
      p1001 = ( Table( il1+1, il2  , il3  , il4+1 ) )
      p1010 = ( Table( il1+1, il2  , il3+1, il4   ) )
      p1011 = ( Table( il1+1, il2  , il3+1, il4+1 ) )
      p1100 = ( Table( il1+1, il2+1, il3  , il4   ) )
      p1101 = ( Table( il1+1, il2+1, il3  , il4+1 ) )
      p1110 = ( Table( il1+1, il2+1, il3+1, il4   ) )
      p1111 = ( Table( il1+1, il2+1, il3+1, il4+1 ) )

      IF (debug) THEN
        WRITE(*,*) ' 16 corners are loaded. Then are:'
        WRITE(*,*) p0000, p0001, &
                   p0010, p0011, p0100, p0101, p0110, p0111, p1000, &
                   p1001, p1010, p1011, p1100, p1101, p1110, p1111
      END IF

      IF ( LogInterp(1) == 1 ) THEN
        alpha(1) &
          = 1.0_dp / ( x1(i) * LOG10( Coordinate1(il1+1) / Coordinate1(il1) ) )
        delta(1) &
          = LOG10( x1(i) / Coordinate1(il1) ) &
              / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
        alpha(1) &
          = ln10 / ( Coordinate1(il1+1) - Coordinate1(il1) )
        delta(1) &
          = ( x1(i) - Coordinate1(il1) ) &
              / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        alpha(2) &
          = 1.0_dp / ( x2(i) * LOG10( Coordinate2(il2+1) / Coordinate2(il2) ) )
        delta(2) &
          = LOG10( x2(i) / Coordinate2(il2) ) &
              / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
        alpha(2) &
          = ln10 / ( Coordinate2(il2+1) - Coordinate2(il2) )
        delta(2) &
          = ( x2(i) - Coordinate2(il2) ) &
              / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        alpha(3) &
          = 1.0_dp / ( x3(i) * LOG10( Coordinate3(il3+1) / Coordinate3(il3) ) )
        delta(3) &
          = LOG10( x3(i) / Coordinate3(il3) ) &
              / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
        alpha(3) &
          = ln10 / ( Coordinate3(il3+1) - Coordinate3(il3) )
        delta(3) &
          = ( x3(i) - Coordinate3(il3) ) &
              / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      IF ( LogInterp(4) == 1 ) THEN
        alpha(4) &
          = 1.0_dp / ( x4(i) * LOG10( Coordinate4(il4+1) / Coordinate4(il4) ) )
        delta(4) &
          = LOG10( x4(i) / Coordinate4(il4) ) &
              / LOG10( Coordinate4(il4+1) / Coordinate4(il4) )
      ELSE
        alpha(4) &
          = ln10 / ( Coordinate4(il4+1) - Coordinate4(il4) )
        delta(4) &
          = ( x4(i) - Coordinate4(il4) ) &
              / ( Coordinate4(il4+1) - Coordinate4(il4) )
      END IF
     
      IF (debug) THEN
        WRITE(*,*) '  Alpha and delta are calculated.'
        WRITE(*,*) ' delta = ', delta
        WRITE(*,*) ' alpha = ', alpha
      END IF

      Interpolant(i) &
        = 10.d0**( &
            (1.0_dp - delta(4)) &
              * (   (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0000   &
                  + (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1000   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0100   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  *           delta(1)  * p1100   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0010   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1010   &
                  +           delta(3)  * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0110   &
                  +           delta(3)  * &
                                    delta(2)  *           delta(1)  * p1110 ) &
            +         delta(4)  &
              * (   (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0001   &
                  + (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1001   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0101   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  *           delta(1)  * p1101   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0011   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1011   &
                  +           delta(3)  * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0111   &
                  +           delta(3)  * &
                                    delta(2)  *           delta(1)  * p1111 ) ) &
          - Offset

      IF (debug) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'Interpolant is calculated =', Interpolant(i)
        WRITE(*,*) 'LOG Interpolant is calculated =', LOG10( Interpolant(i) )
        WRITE(*,*) ''
      END IF

      Derivative(i,1) &  ! E
        = ( Interpolant(i) ) * alpha(1) &
          * ( (1.0_dp - delta(4)) * ( (1.0_dp - delta(3) ) * &
                                      ( ( delta(2) - 1.0_dp )  * p0000   &
                                      + ( 1.0_dp   - delta(2)) * p1000   &
                                                   - delta(2)  * p0100   &
                                                   + delta(2)  * p1100 ) &
                                    +           delta(3)   * &
                                      ( ( delta(2) - 1.0_dp )  * p0010   &
                                      - ( delta(2) - 1.0_dp )  * p1010   &
                                                   - delta(2)  * p0110   &
                                        +            delta(2)  * p1110 ))&
                      + delta(4) * ( (1.0_dp - delta(3) )  * &
                                     ( ( delta(2) - 1.0_dp )   * p0001   &
                                     + ( 1.0_dp   - delta(2))  * p1001   &
                                                  - delta(2)   * p0101   &
                                                  + delta(2)   * p1101 ) &
                                     +           delta(3)  * &
                                      ( ( delta(2) - 1.0_dp )  * p0011   &
                                      + ( 1.0_dp - delta(2) )  * p1011   &
                                                  - delta(2)   * p0111   &
                                                  + delta(2)   * p1111)) )
     
      Derivative(i,2) &  ! rho
        = ( Interpolant(i) ) * alpha(2) &
          * ( &
            (1.0_dp - delta(4)) &
              * (   (1.0_dp - delta(3)) * &
                                    ( -1.0_dp + delta(1)) * p0000   &
                  + (1.0_dp - delta(3)) * &
                                    (         - delta(1)) * p1000   &
                  + (1.0_dp - delta(3)) * &
                                    (  1.0_dp - delta(1)) * p0100   &
                  + (1.0_dp - delta(3)) * &
                                                delta(1)  * p1100   &
                  +           delta(3)  * &
                                    ( -1.0_dp + delta(1)) * p0010   &
                  +           delta(3)  * &
                                    (         - delta(1)) * p1010   &
                  +           delta(3)  * &
                                    (  1.0_dp - delta(1)) * p0110   &
                  +           delta(3)  * &
                                                delta(1)  * p1110 ) &
            +         delta(4)  &
              * (   (1.0_dp - delta(3)) * &
                                    ( -1.0_dp + delta(1)) * p0001   &
                  + (1.0_dp - delta(3)) * &
                                    (         - delta(1)) * p1001   &
                  + (1.0_dp - delta(3)) * &
                                    (  1.0_dp - delta(1)) * p0101   &
                  + (1.0_dp - delta(3)) * &
                                                delta(1)  * p1101   &
                  +           delta(3)  * &
                                    ( -1.0_dp - delta(1)) * p0011   &
                  +           delta(3)  * &
                                    (         - delta(1)) * p1011   &
                  +           delta(3)  * &
                                    (  1.0_dp - delta(1)) * p0111   &
                  +           delta(3)  * &
                                                delta(1)  * p1111 ) )

      Derivative(i,3) &  ! T
        = ( Interpolant(i) )  * alpha(3) &
          * ( (1.0_dp - delta(4)) &
              * (      (  (delta(2) - 1.0_dp) * (1.0_dp - delta(1)) * p0000   &
                  +       (delta(2) - 1.0_dp) *           delta(1)  * p1000   &
                  +                 delta(2)  * (delta(1) - 1.0_dp) * p0100   &
                  -                 delta(2)  *           delta(1)  * p1100 ) &
                  +    (  (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0010   &
                  +       (1.0_dp - delta(2)) *           delta(1)  * p1010   &
                  +                 delta(2)  * (1.0_dp - delta(1)) * p0110   &
                  +                 delta(2)  *           delta(1)  * p1110)) &
            +           delta(4)  &
              * (      (  (delta(2) - 1.0_dp) * (1.0_dp - delta(1)) * p0001   &
                  +       (delta(2) - 1.0_dp) *           delta(1)  * p1001   &
                  +                 delta(2)  * (delta(1) - 1.0_dp) * p0101   &
                  -                 delta(2)  *           delta(1)  * p1101 ) &
                  +    (  (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0011   &
                  +       (1.0_dp - delta(2)) *           delta(1)  * p1011   &
                  +                 delta(2)  * (1.0_dp - delta(1)) * p0111   &
                  +                 delta(2)  *           delta(1)  * p1111)) ) 

      Derivative(i,4) &  ! Ye
        = ( Interpolant(i) ) * alpha(4) &
          * ( - (   (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0000   &
                  + (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1000   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0100   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  *           delta(1)  * p1100   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0010   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1010   &
                  +           delta(3)  * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0110   &
                  +           delta(3)  * &
                                    delta(2)  *           delta(1)  * p1110 ) &
            +   (   (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0001   &
                  + (1.0_dp - delta(3)) * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1001   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  * (1.0_dp - delta(1)) * p0101   &
                  + (1.0_dp - delta(3)) * &
                                    delta(2)  *           delta(1)  * p1101   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) * (1.0_dp - delta(1)) * p0011   &
                  +           delta(3)  * &
                          (1.0_dp - delta(2)) *           delta(1)  * p1011   &
                  +           delta(3)  * &
                                   delta(2)  * (1.0_dp - delta(1)) * p0111   &
                  +           delta(3)  * &
                                    delta(2)  *           delta(1)  * p1111 ) )
    
    IF (debug) THEN
      WRITE(*,*) ' '
      WRITE(*,*) 'Derivative is calculated =', Derivative(i,:)
      WRITE(*,*) ''
    END IF

    END DO

    IF (debug) THEN
      WRITE(*,*) 'End of differentiate routine'
      WRITE(*,*) ''
    END IF

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D


END MODULE wlInterpolationModule

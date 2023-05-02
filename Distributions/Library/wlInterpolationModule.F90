MODULE wlInterpolationModule

  USE wlKindModule, ONLY: dp
  !USE wlInterpolationUtilitiesModule, ONLY: &
  !  Index1D_Lin, &
  !  GetIndexAndDelta_Lin, &
  !  GetIndexAndDelta_Log, &
  !  LinearInterp1D_1DArray_Point, &
  !  LinearInterp2D_2DArray_Point, &
  !  LinearInterp3D_3DArray_Point, &
  !  LinearInterp4D_4DArray_Point, &
  !  LinearInterp5D_5DArray_Point, &
  !  LinearInterp2D_3DArray_1DAligned_Point, &
  !  LinearInterp3D_4DArray_1DAligned_Point, &
  !  LinearInterp4D_5DArray_1DAligned_Point, &
  !  LinearInterp2D_4DArray_2DAligned_Point, &
  !  LinearInterp3D_5DArray_2DAligned_Point, &
  !  LinearInterpDeriv3D_3DArray_Point, &
  !  LinearInterpDeriv4D_4DArray_Point, &
  !  LinearInterpDeriv2D_4DArray_2DAligned_Point

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

  INTERFACE Index1D_Lin
    MODULE PROCEDURE Index1D_Lin1
    MODULE PROCEDURE Index1D_Lin2
  END INTERFACE Index1D_Lin

  INTERFACE Index1D_Log
    MODULE PROCEDURE Index1D_Log1
    MODULE PROCEDURE Index1D_Log2
  END INTERFACE Index1D_Log

  INTERFACE LinearInterp_Array_Point
    MODULE PROCEDURE LinearInterp1D_1DArray_Point
    MODULE PROCEDURE LinearInterp2D_2DArray_Point
    MODULE PROCEDURE LinearInterp3D_3DArray_Point
    MODULE PROCEDURE LinearInterp4D_4DArray_Point
    MODULE PROCEDURE LinearInterp5D_5DArray_Point
  END INTERFACE LinearInterp_Array_Point

  INTERFACE LinearInterpDeriv_Array_Point
    MODULE PROCEDURE LinearInterpDeriv3D_3DArray_Point
    MODULE PROCEDURE LinearInterpDeriv4D_4DArray_Point
  END INTERFACE LinearInterpDeriv_Array_Point

CONTAINS


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


  INTEGER FUNCTION Index1D_Lin1( x, xx, n )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    Index1D_Lin1 = 1 + FLOOR( (n-1)*(x-xx(1))/(xx(n)-xx(1)) )
    Index1D_Lin1 = MAX( 1, MIN( n-1, Index1D_Lin1 ) )

    RETURN
  END FUNCTION Index1D_Lin1


  INTEGER FUNCTION Index1D_Lin2( x, xx )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x, xx(:)
    INTEGER :: lo, hi

    lo = LBOUND(xx,1)
    hi = UBOUND(xx,1)

    Index1D_Lin2 = lo + FLOOR( (hi-lo)*(x-xx(lo))/(xx(hi)-xx(lo)) )
    Index1D_Lin2 = MAX( lo, MIN( hi-1, Index1D_Lin2 ) )

    RETURN
  END FUNCTION Index1D_Lin2


  INTEGER FUNCTION Index1D_Log1( x, xx, n )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    Index1D_Log1 = 1 + FLOOR( (n-1)*LOG10(x/xx(1))/LOG10(xx(n)/xx(1)) )
    Index1D_Log1 = MAX( 1, MIN( n-1, Index1D_Log1 ) )

    RETURN
  END FUNCTION Index1D_Log1


  INTEGER FUNCTION Index1D_Log2( x, xx )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x, xx(:)
    INTEGER :: lo, hi

    lo = LBOUND(xx,1)
    hi = UBOUND(xx,1)

    Index1D_Log2 = lo + FLOOR( (hi-lo)*LOG10(x/xx(lo))/LOG10(xx(hi)/xx(lo)) )
    Index1D_Log2 = MAX( lo, MIN( hi-1, Index1D_Log2 ) )

    RETURN
  END FUNCTION Index1D_Log2


  SUBROUTINE GetIndexAndDelta_Lin( Y, Ys, iY, dY )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    
    REAL(dp), INTENT(in)  :: Y, Ys(1:)
    INTEGER,  INTENT(out) :: iY
    REAL(dp), INTENT(out) :: dY

    iY = Index1D_Lin( Y, Ys )
    dY = ( Y - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

  END SUBROUTINE GetIndexAndDelta_Lin


  SUBROUTINE GetIndexAndDelta_Log( Y, Ys, iY, dY )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    
    REAL(dp), INTENT(in)  :: Y, Ys(1:)
    INTEGER,  INTENT(out) :: iY
    REAL(dp), INTENT(out) :: dY

    iY = Index1D_Log( Y, Ys )
    dY = LOG10( Y / Ys(iY) ) / LOG10( Ys(iY+1) / Ys(iY) )

  END SUBROUTINE GetIndexAndDelta_Log


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

    Linear = ddX1 * p0 + dX1 * p1 

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

    BiLinear &
      = Linear &
      ( Linear( p00, p10, dX1 ), &
        Linear( p01, p11, dX1 ), &
        dX2 )

    RETURN
  END FUNCTION BiLinear


  REAL(dp) FUNCTION dBiLineardX1 &
    ( p00, p10, p01, p11, dX2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p00, p10, p01, p11, dX2

    dBiLineardX1 &
      = Linear(p10, p11, dX2 ) &
      - Linear(p00, p01, dX2 )

    RETURN
  END FUNCTION dBiLineardX1


  REAL(dp) FUNCTION dBiLineardX2 &
    ( p00, p10, p01, p11, dX1 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p00, p10, p01, p11, dX1

    dBiLineardX2 &
      = Linear(p01, p11, dX1 ) &
      - Linear(p00, p10, dX1 )

    RETURN
  END FUNCTION dBiLineardX2


  REAL(dp) FUNCTION TriLinear &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2, dX3 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2, dX3

    TriLinear &
      = Linear &
      ( BiLinear( p000, p100, p010, p110, dX1, dX2 ), &
        BiLinear( p001, p101, p011, p111, dX1, dX2 ), &
        dX3 )

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
      p000, p100, p010, p110, p001, p101, p011, p111, dX2, dX3

    dTrilineardX1 &
      = Bilinear(p100, p110, p101, p111, dX2, dX3) &
      - Bilinear(p000, p010, p001, p011, dX2, dX3)

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
      p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX3

    dTrilineardX2 &
      = Bilinear(p010, p110, p011, p111, dX1, dX3) &
      - Bilinear(p000, p100, p001, p101, dX1, dX3)

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
      p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2

    dTrilineardX3 &
      = Bilinear(p001, p101, p011, p111, dX1, dX2) &
      - Bilinear(p000, p100, p010, p110, dX1, dX2)

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

    TetraLinear &
      = Linear &
      ( TriLinear( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, dX1, dX2, dX3 ), &
        TriLinear( p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, dX1, dX2, dX3 ), &
        dX4 )

    RETURN
  END FUNCTION TetraLinear


  REAL(dp) FUNCTION dTetraLineardX1 &
    ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX2, dX3, dX4 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX2, dX3, dX4

    dTetraLineardX1 &
      = TriLinear( p1000, p1100, p1010, p1110, p1001, p1101, p1011, p1111, dX2, dX3, dX4 ) &
      - TriLinear( p0000, p0100, p0010, p0110, p0001, p0101, p0011, p0111, dX2, dX3, dX4 )

    RETURN
  END FUNCTION dTetraLineardX1


  REAL(dp) FUNCTION dTetraLineardX2 &
    ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX3, dX4 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX3, dX4

    dTetraLineardX2 &
      = TriLinear( p0100, p1100, p0110, p1110, p0101, p1101, p0111, p1111, dX1, dX3, dX4 ) &
      - TriLinear( p0000, p1000, p0010, p1010, p0001, p1001, p0011, p1011, dX1, dX3, dX4 )

    RETURN
  END FUNCTION dTetraLineardX2


  REAL(dp) FUNCTION dTetraLineardX3 &
    ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX4 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX4

    dTetraLineardX3 &
      = TriLinear( p0010, p1010, p0110, p1110, p0011, p1011, p0111, p1111, dX1, dX2, dX4 ) &
      - TriLinear( p0000, p1000, p0100, p1100, p0001, p1001, p0101, p1101, dX1, dX2, dX4 )

    RETURN
  END FUNCTION dTetraLineardX3


  REAL(dp) FUNCTION dTetraLineardX4 &
    ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3

    dTetraLineardX4 &
      = TriLinear( p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, dX1, dX2, dX3 ) &
      - TriLinear( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, dX1, dX2, dX3 )

    RETURN
  END FUNCTION dTetraLineardX4


  REAL(dp) FUNCTION PentaLinear &
    ( p00000, p10000, p01000, p11000, p00100, p10100, p01100, p11100, &
      p00010, p10010, p01010, p11010, p00110, p10110, p01110, p11110, &
      p00001, p10001, p01001, p11001, p00101, p10101, p01101, p11101, &
      p00011, p10011, p01011, p11011, p00111, p10111, p01111, p11111, &
      dX1, dX2, dX3, dX4, dX5 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p00000, p10000, p01000, p11000, p00100, p10100, p01100, p11100, &
      p00010, p10010, p01010, p11010, p00110, p10110, p01110, p11110, &
      p00001, p10001, p01001, p11001, p00101, p10101, p01101, p11101, &
      p00011, p10011, p01011, p11011, p00111, p10111, p01111, p11111, &
      dX1, dX2, dX3, dX4, dX5

    PentaLinear &
      = Linear &
      ( TetraLinear( p00000, p10000, p01000, p11000, p00100, p10100, p01100, p11100, &
                     p00010, p10010, p01010, p11010, p00110, p10110, p01110, p11110, &
                     dX1, dX2, dX3, dX4 ), &
        TetraLinear( p00001, p10001, p01001, p11001, p00101, p10101, p01101, p11101, &
                     p00011, p10011, p01011, p11011, p00111, p10111, p01111, p11111, &
                     dX1, dX2, dX3, dX4 ), &
        dX5 )

    RETURN
  END FUNCTION PentaLinear


  SUBROUTINE LinearInterp1D_1DArray_Point &
    ( iY1, dY1, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iY1
    REAL(dp), INTENT(in)  :: dY1, OS, Table(:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p0, p1

    p0 = Table(iY1  )
    p1 = Table(iY1+1)

    Interpolant &
      = 10.0d0**( &
          Linear &
            ( p0, p1, &
              dY1 ) ) - OS

  END SUBROUTINE LinearInterp1D_1DArray_Point


  SUBROUTINE LinearInterp2D_2DArray_Point &
    ( iY1, iY2, dY1, dY2, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iY1, iY2
    REAL(dp), INTENT(in)  :: dY1, dY2, OS, Table(1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p00, p10, p01, p11

    p00 = Table(iY1  , iY2  )
    p10 = Table(iY1+1, iY2  )
    p01 = Table(iY1  , iY2+1)
    p11 = Table(iY1+1, iY2+1)

    Interpolant &
      = 10.0d0**( &
          BiLinear &
            ( p00, p10, p01, p11, &
              dY1, dY2 ) ) - OS

  END SUBROUTINE LinearInterp2D_2DArray_Point


  SUBROUTINE LinearInterp2D_3DArray_1DAligned_Point &
    ( iX1, iY1, iY2, dY1, dY2, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iX1, iY1, iY2
    REAL(dp), INTENT(in)  :: dY1, dY2, OS, Table(1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p00, p10, p01, p11

    p00 = Table(iX1, iY1  , iY2  )
    p10 = Table(iX1, iY1+1, iY2  )
    p01 = Table(iX1, iY1  , iY2+1)
    p11 = Table(iX1, iY1+1, iY2+1)

    Interpolant &
      = 10.0d0**( &
          BiLinear &
            ( p00, p10, p01, p11, &
              dY1, dY2 ) ) - OS

  END SUBROUTINE LinearInterp2D_3DArray_1DAligned_Point


  SUBROUTINE LinearInterp2D_4DArray_2DAligned_Point &
    ( iX1, iX2, iY1, iY2, dY1, dY2, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iX1, iX2, iY1, iY2
    REAL(dp), INTENT(in)  :: dY1, dY2, OS, Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p00, p10, p01, p11

    p00 = Table(iX1, iX2, iY1  , iY2  )
    p10 = Table(iX1, iX2, iY1+1, iY2  )
    p01 = Table(iX1, iX2, iY1  , iY2+1)
    p11 = Table(iX1, iX2, iY1+1, iY2+1)

    Interpolant &
      = 10.0d0**( &
          BiLinear &
            ( p00, p10, p01, p11, &
              dY1, dY2 ) ) - OS

  END SUBROUTINE LinearInterp2D_4DArray_2DAligned_Point


  SUBROUTINE LinearInterp3D_3DArray_Point &
    ( iY1, iY2, iY3, dY1, dY2, dY3, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iY1, iY2, iY3
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, OS, Table(1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p000, p100, p010, p110, p001, p101, p011, p111

    p000 = Table(iY1  , iY2  , iY3   )
    p100 = Table(iY1+1, iY2  , iY3   )
    p010 = Table(iY1  , iY2+1, iY3   )
    p110 = Table(iY1+1, iY2+1, iY3   )
    p001 = Table(iY1  , iY2  , iY3+1 )
    p101 = Table(iY1+1, iY2  , iY3+1 )
    p011 = Table(iY1  , iY2+1, iY3+1 )
    p111 = Table(iY1+1, iY2+1, iY3+1 )

    Interpolant &
      = 10.0d0**( &
          TriLinear &
            ( p000, p100, p010, p110, &
              p001, p101, p011, p111, &
              dY1, dY2, dY3 ) ) - OS

  END SUBROUTINE LinearInterp3D_3DArray_Point


  SUBROUTINE LinearInterp3D_4DArray_1DAligned_Point &
    ( iX1, iY1, iY2, iY3, dY1, dY2, dY3, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iX1, iY1, iY2, iY3
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, OS, Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p000, p100, p010, p110, p001, p101, p011, p111

    p000 = Table(iX1, iY1  , iY2  , iY3   )
    p100 = Table(iX1, iY1+1, iY2  , iY3   )
    p010 = Table(iX1, iY1  , iY2+1, iY3   )
    p110 = Table(iX1, iY1+1, iY2+1, iY3   )
    p001 = Table(iX1, iY1  , iY2  , iY3+1 )
    p101 = Table(iX1, iY1+1, iY2  , iY3+1 )
    p011 = Table(iX1, iY1  , iY2+1, iY3+1 )
    p111 = Table(iX1, iY1+1, iY2+1, iY3+1 )

    Interpolant &
      = 10.0d0**( &
          TriLinear &
            ( p000, p100, p010, p110, &
              p001, p101, p011, p111, &
              dY1, dY2, dY3 ) ) - OS

  END SUBROUTINE LinearInterp3D_4DArray_1DAligned_Point


  SUBROUTINE LinearInterp3D_5DArray_2DAligned_Point &
    ( iX1, iX2, iY1, iY2, iY3, dY1, dY2, dY3, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iX1, iX2, iY1, iY2, iY3
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, OS, Table(1:,1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p000, p100, p010, p110, p001, p101, p011, p111

    p000 = Table(iX1, iX2, iY1  , iY2  , iY3   )
    p100 = Table(iX1, iX2, iY1+1, iY2  , iY3   )
    p010 = Table(iX1, iX2, iY1  , iY2+1, iY3   )
    p110 = Table(iX1, iX2, iY1+1, iY2+1, iY3   )
    p001 = Table(iX1, iX2, iY1  , iY2  , iY3+1 )
    p101 = Table(iX1, iX2, iY1+1, iY2  , iY3+1 )
    p011 = Table(iX1, iX2, iY1  , iY2+1, iY3+1 )
    p111 = Table(iX1, iX2, iY1+1, iY2+1, iY3+1 )

    Interpolant &
      = 10.0d0**( &
          TriLinear &
            ( p000, p100, p010, p110, &
              p001, p101, p011, p111, &
              dY1, dY2, dY3 ) ) - OS

  END SUBROUTINE LinearInterp3D_5DArray_2DAligned_Point


  SUBROUTINE LinearInterp4D_4DArray_Point &
    ( iY1, iY2, iY3, iY4, dY1, dY2, dY3, dY4, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iY1, iY2, iY3, iY4
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, dY4, OS, Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    p0000 = Table(iY1  , iY2  , iY3  , iY4   )
    p1000 = Table(iY1+1, iY2  , iY3  , iY4   )
    p0100 = Table(iY1  , iY2+1, iY3  , iY4   )
    p1100 = Table(iY1+1, iY2+1, iY3  , iY4   )
    p0010 = Table(iY1  , iY2  , iY3+1, iY4   )
    p1010 = Table(iY1+1, iY2  , iY3+1, iY4   )
    p0110 = Table(iY1  , iY2+1, iY3+1, iY4   )
    p1110 = Table(iY1+1, iY2+1, iY3+1, iY4   )
    p0001 = Table(iY1  , iY2  , iY3  , iY4+1 )
    p1001 = Table(iY1+1, iY2  , iY3  , iY4+1 )
    p0101 = Table(iY1  , iY2+1, iY3  , iY4+1 )
    p1101 = Table(iY1+1, iY2+1, iY3  , iY4+1 )
    p0011 = Table(iY1  , iY2  , iY3+1, iY4+1 )
    p1011 = Table(iY1+1, iY2  , iY3+1, iY4+1 )
    p0111 = Table(iY1  , iY2+1, iY3+1, iY4+1 )
    p1111 = Table(iY1+1, iY2+1, iY3+1, iY4+1 )

    Interpolant &
      = 10.0d0**( &
          TetraLinear &
            ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
              p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
              dY1, dY2, dY3, dY4 ) ) - OS

  END SUBROUTINE LinearInterp4D_4DArray_Point


  SUBROUTINE LinearInterp4D_5DArray_1DAligned_Point &
    ( iX1, iY1, iY2, iY3, iY4, dY1, dY2, dY3, dY4, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iX1, iY1, iY2, iY3, iY4
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, dY4, OS, Table(1:,1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    p0000 = Table(iX1, iY1  , iY2  , iY3  , iY4   )
    p1000 = Table(iX1, iY1+1, iY2  , iY3  , iY4   )
    p0100 = Table(iX1, iY1  , iY2+1, iY3  , iY4   )
    p1100 = Table(iX1, iY1+1, iY2+1, iY3  , iY4   )
    p0010 = Table(iX1, iY1  , iY2  , iY3+1, iY4   )
    p1010 = Table(iX1, iY1+1, iY2  , iY3+1, iY4   )
    p0110 = Table(iX1, iY1  , iY2+1, iY3+1, iY4   )
    p1110 = Table(iX1, iY1+1, iY2+1, iY3+1, iY4   )
    p0001 = Table(iX1, iY1  , iY2  , iY3  , iY4+1 )
    p1001 = Table(iX1, iY1+1, iY2  , iY3  , iY4+1 )
    p0101 = Table(iX1, iY1  , iY2+1, iY3  , iY4+1 )
    p1101 = Table(iX1, iY1+1, iY2+1, iY3  , iY4+1 )
    p0011 = Table(iX1, iY1  , iY2  , iY3+1, iY4+1 )
    p1011 = Table(iX1, iY1+1, iY2  , iY3+1, iY4+1 )
    p0111 = Table(iX1, iY1  , iY2+1, iY3+1, iY4+1 )
    p1111 = Table(iX1, iY1+1, iY2+1, iY3+1, iY4+1 )

    Interpolant &
      = 10.0d0**( &
          TetraLinear &
            ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
              p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
              dY1, dY2, dY3, dY4 ) ) - OS

  END SUBROUTINE LinearInterp4D_5DArray_1DAligned_Point


  SUBROUTINE LinearInterp5D_5DArray_Point &
    ( iY1, iY2, iY3, iY4, iY5, dY1, dY2, dY3, dY4, dY5, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iY1, iY2, iY3, iY4, iY5
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, dY4, dY5, OS, Table(1:,1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    REAL(dp) :: p00000, p00001, p00010, p00011, p00100, p00101, p00110, p00111, &
                p01000, p01001, p01010, p01011, p01100, p01101, p01110, p01111, &
                p10000, p10001, p10010, p10011, p10100, p10101, p10110, p10111, &
                p11000, p11001, p11010, p11011, p11100, p11101, p11110, p11111

    p00000 = Table(iY1,   iY2  , iY3  , iY4  , iY5   )
    p01000 = Table(iY1,   iY2+1, iY3  , iY4  , iY5   )
    p00100 = Table(iY1,   iY2  , iY3+1, iY4  , iY5   )
    p01100 = Table(iY1,   iY2+1, iY3+1, iY4  , iY5   )
    p00010 = Table(iY1,   iY2  , iY3  , iY4+1, iY5   )
    p01010 = Table(iY1,   iY2+1, iY3  , iY4+1, iY5   )
    p00110 = Table(iY1,   iY2  , iY3+1, iY4+1, iY5   )
    p01110 = Table(iY1,   iY2+1, iY3+1, iY4+1, iY5   )
    p00001 = Table(iY1,   iY2  , iY3  , iY4  , iY5+1 )
    p01001 = Table(iY1,   iY2+1, iY3  , iY4  , iY5+1 )
    p00101 = Table(iY1,   iY2  , iY3+1, iY4  , iY5+1 )
    p01101 = Table(iY1,   iY2+1, iY3+1, iY4  , iY5+1 )
    p00011 = Table(iY1,   iY2  , iY3  , iY4+1, iY5+1 )
    p01011 = Table(iY1,   iY2+1, iY3  , iY4+1, iY5+1 )
    p00111 = Table(iY1,   iY2  , iY3+1, iY4+1, iY5+1 )
    p01111 = Table(iY1,   iY2+1, iY3+1, iY4+1, iY5+1 )
    p10000 = Table(iY1+1, iY2  , iY3  , iY4  , iY5   )
    p11000 = Table(iY1+1, iY2+1, iY3  , iY4  , iY5   )
    p10100 = Table(iY1+1, iY2  , iY3+1, iY4  , iY5   )
    p11100 = Table(iY1+1, iY2+1, iY3+1, iY4  , iY5   )
    p10010 = Table(iY1+1, iY2  , iY3  , iY4+1, iY5   )
    p11010 = Table(iY1+1, iY2+1, iY3  , iY4+1, iY5   )
    p10110 = Table(iY1+1, iY2  , iY3+1, iY4+1, iY5   )
    p11110 = Table(iY1+1, iY2+1, iY3+1, iY4+1, iY5   )
    p10001 = Table(iY1+1, iY2  , iY3  , iY4  , iY5+1 )
    p11001 = Table(iY1+1, iY2+1, iY3  , iY4  , iY5+1 )
    p10101 = Table(iY1+1, iY2  , iY3+1, iY4  , iY5+1 )
    p11101 = Table(iY1+1, iY2+1, iY3+1, iY4  , iY5+1 )
    p10011 = Table(iY1+1, iY2  , iY3  , iY4+1, iY5+1 )
    p11011 = Table(iY1+1, iY2+1, iY3  , iY4+1, iY5+1 )
    p10111 = Table(iY1+1, iY2  , iY3+1, iY4+1, iY5+1 )
    p11111 = Table(iY1+1, iY2+1, iY3+1, iY4+1, iY5+1 )

    Interpolant &
      = 10.0d0**( &
          PentaLinear &
            ( p00000, p01000, p00100, p01100, p00010, p01010, p00110, p01110, &
              p00001, p01001, p00101, p01101, p00011, p01011, p00111, p01111, &
              p10000, p11000, p10100, p11100, p10010, p11010, p10110, p11110, &
              p10001, p11001, p10101, p11101, p10011, p11011, p10111, p11111, &
              dY1, dY2, dY3, dY4, dY5 ) ) - OS

  END SUBROUTINE LinearInterp5D_5DArray_Point


  SUBROUTINE LinearInterpDeriv2D_4DArray_2DAligned_Point &
    ( iX1, iX2, iY1, iY2, dY1, dY2, aY1, aY2, OS, Table, &
      Interpolant, dIdY1, dIdY2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: iX1, iX2, iY1, iY2
    REAL(dp), INTENT(in) :: dY1, dY2, aY1, aY2, OS, Table(1:,1:,1:,1:)
    REAL(dp) :: Interpolant, dIdY1, dIdY2

    REAL(dp) :: p00, p10, p01, p11

    p00 = Table(iX1, iX2, iY1  , iY2  )
    p10 = Table(iX1, iX2, iY1+1, iY2  )
    p01 = Table(iX1, iX2, iY1  , iY2+1)
    p11 = Table(iX1, iX2, iY1+1, iY2+1)

    Interpolant &
      = 10.0d0**( &
          BiLinear &
            ( p00, p10, p01, p11, &
              dY1, dY2) ) - OS

    dIdY1 &
      = (Interpolant + OS) * aY1 &
          * dBiLineardX1( p00, p10, p01, p11, dY2 )

    dIdY2 &
      = (Interpolant + OS) * aY2 &
          * dBiLineardX2( p00, p10, p01, p11, dY1 )

  END SUBROUTINE LinearInterpDeriv2D_4DArray_2DAligned_Point


  SUBROUTINE LinearInterpDeriv3D_3DArray_Point &
    ( iY1, iY2, iY3, dY1, dY2, dY3, aY1, aY2, aY3, OS, Table, &
      Interpolant, dIdY1, dIdY2, dIdY3 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iY1, iY2, iY3
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, aY1, aY2, aY3, OS, Table(1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant, dIdY1, dIdY2, dIdY3

    REAL(dp) :: p000, p100, p010, p110, p001, p101, p011, p111

    p000 = Table(iY1  , iY2  , iY3   )
    p100 = Table(iY1+1, iY2  , iY3   )
    p010 = Table(iY1  , iY2+1, iY3   )
    p110 = Table(iY1+1, iY2+1, iY3   )
    p001 = Table(iY1  , iY2  , iY3+1 )
    p101 = Table(iY1+1, iY2  , iY3+1 )
    p011 = Table(iY1  , iY2+1, iY3+1 )
    p111 = Table(iY1+1, iY2+1, iY3+1 )

    Interpolant &
      = 10.0d0**( &
          TriLinear &
            ( p000, p100, p010, p110, &
              p001, p101, p011, p111, &
              dY1, dY2, dY3 ) ) - OS

    dIdY1 &
      = (Interpolant + OS) * aY1 &
          * dTriLineardX1 &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, &
                dY2, dY3 )

    dIdY2 &
      = (Interpolant + OS) * aY2 &
          * dTriLineardX2 &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, &
                dY1, dY3 )

    dIdY3 &
      = (Interpolant + OS) * aY3 &
          * dTriLineardX3 &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, &
                dY1, dY2 )

  END SUBROUTINE LinearInterpDeriv3D_3DArray_Point


  SUBROUTINE LinearInterpDeriv4D_4DArray_Point &
    ( iY1, iY2, iY3, iY4, dY1, dY2, dY3, dY4, aY1, aY2, aY3, aY4, OS, Table, &
      Interpolant, dIdY1, dIdY2, dIdY3, dIdY4 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iY1, iY2, iY3, iY4
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, dY4, aY1, aY2, aY3, aY4, OS, Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant, dIdY1, dIdY2, dIdY3, dIdY4

    REAL(dp) :: p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
                p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    p0000 = Table(iY1  , iY2  , iY3  , iY4   )
    p1000 = Table(iY1+1, iY2  , iY3  , iY4   )
    p0100 = Table(iY1  , iY2+1, iY3  , iY4   )
    p1100 = Table(iY1+1, iY2+1, iY3  , iY4   )
    p0010 = Table(iY1  , iY2  , iY3+1, iY4   )
    p1010 = Table(iY1+1, iY2  , iY3+1, iY4   )
    p0110 = Table(iY1  , iY2+1, iY3+1, iY4   )
    p1110 = Table(iY1+1, iY2+1, iY3+1, iY4   )
    p0001 = Table(iY1  , iY2  , iY3  , iY4+1 )
    p1001 = Table(iY1+1, iY2  , iY3  , iY4+1 )
    p0101 = Table(iY1  , iY2+1, iY3  , iY4+1 )
    p1101 = Table(iY1+1, iY2+1, iY3  , iY4+1 )
    p0011 = Table(iY1  , iY2  , iY3+1, iY4+1 )
    p1011 = Table(iY1+1, iY2  , iY3+1, iY4+1 )
    p0111 = Table(iY1  , iY2+1, iY3+1, iY4+1 )
    p1111 = Table(iY1+1, iY2+1, iY3+1, iY4+1 )

    Interpolant &
      = 10.0d0**( &
          TetraLinear &
            ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
              p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
              dY1, dY2, dY3, dY4 ) ) - OS

    dIdY1 &
      = (Interpolant + OS) * aY1 &
          * dTetraLineardX1 &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dY2, dY3, dY4 )

    dIdY2 &
      = (Interpolant + OS) * aY2 &
          * dTetraLineardX2 &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dY1, dY3, dY4 )

    dIdY3 &
      = (Interpolant + OS) * aY3 &
          * dTetraLineardX3 &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dY1, dY2, dY4 )

    dIdY4 &
      = (Interpolant + OS) * aY4 &
          * dTetraLineardX4 &
              ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
                p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
                dY1, dY2, dY3 )

  END SUBROUTINE LinearInterpDeriv4D_4DArray_Point


  SUBROUTINE LogInterpolateSingleVariable_2D_Custom &
    ( X, Y, Xs, Ys, OS, Table, Interpolant, Error_Option )

    REAL(dp), INTENT(in)  :: X (1:), Y (1:)
    REAL(dp), INTENT(in)  :: Xs(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER :: iP, Error

    Error = 0
    IF( .NOT. SIZE(X) == SIZE(Y) )THEN
      Error = 1
      IF( PRESENT( Error_Option ) ) Error_Option = Error
      RETURN
    END IF

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: X, Xs, Y, Ys, OS, Table ) &
    !$OMP MAP( from: Interpolant )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( X, Xs, Y, Ys, OS, Table ) &
    !$ACC COPYOUT( Interpolant )
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

    REAL(dp), INTENT(in)  :: X     , Y
    REAL(dp), INTENT(in)  :: Xs(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    INTEGER  :: iX, iY
    REAL(dp) :: dX, dY
    REAL(dp) :: p00, p10, p01, p11
    INTEGER :: loX, hiX
    INTEGER :: loY, hiY

    loX = LBOUND(Xs,1)
    hiX = UBOUND(Xs,1)
    loY = LBOUND(Ys,1)
    hiY = UBOUND(Ys,1)

    !CALL GetIndexAndDelta_Lin( X, Xs, iX, dX )
    !iX = Index1D_Lin( X, Xs )
    iX = MAX( loX, MIN( hiX-1, loX + FLOOR( (hiX-loX)*(X-Xs(loX))/(Xs(hiX)-Xs(loX)) ) ) )
    dX = ( X - Xs(iX) ) / ( Xs(iX+1) - Xs(iX) )
    !CALL GetIndexAndDelta_Lin( Y, Ys, iY, dY )
    !iY = Index1D_Lin( Y, Ys )
    iY = MAX( loY, MIN( hiY-1, loY + FLOOR( (hiY-loY)*(Y-Ys(loY))/(Ys(hiY)-Ys(loY)) ) ) )
    dY = ( Y - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

    !CALL LinearInterp2D_2DArray_Point &
    !       ( iX, iY, dX, dY, OS, Table, Interpolant )

    p00 = Table(iX  , iY  )
    p10 = Table(iX+1, iY  )
    p01 = Table(iX  , iY+1)
    p11 = Table(iX+1, iY+1)
    Interpolant &
      = 10.0d0 ** (   ( One - dY ) * ( ( One - dX ) * p00 + dX * p10 ) &
                    +         dY   * ( ( One - dX ) * p01 + dX * p11 ) ) - OS
    !Interpolant &
    !  = 10.0d0**( &
    !      BiLinear &
    !        ( p00, p10, p01, p11, &
    !          dY1, dY2 ) ) - OS

  END SUBROUTINE LogInterpolateSingleVariable_2D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom &
    ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, OS, Table, Interpolant )

    REAL(dp), INTENT(in)  :: LogE (1:), LogD (1:), LogT (1:), Y (1:)
    REAL(dp), INTENT(in)  :: LogEs(1:), LogDs(1:), LogTs(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:)

    INTEGER  :: i, j
    INTEGER  :: iD, iT, iY, iE
    REAL(dp) :: dD, dT, dY, dE

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iE, iD, dD, iT, dE, dT, iY, dY ) &
    !$OMP MAP( to: LogE, LogEs, LogD, LogDs, LogT, LogTs, Y, Ys, OS, Table ) &
    !$OMP MAP( from: Interpolant )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC PRIVATE( iE, iD, dD, iT, dE, dT, iY, dY ) &
    !$ACC COPYIN( LogE, LogEs, LogD, LogDs, LogT, LogTs, Y, Ys, OS, Table ) &
    !$ACC COPYOUT( Interpolant )
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

    REAL(dp), INTENT(in)  :: LogE (1:), LogD     , LogT     , Y
    REAL(dp), INTENT(in)  :: LogEs(1:), LogDs(1:), LogTs(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:)

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
    ( LogE, LogT, LogX, LogEs, LogTs, LogXs, OS, Table, Interpolant, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogE (1:), LogT (1:), LogX (1:)
    REAL(dp), INTENT(in)  :: LogEs(1:), LogTs(1:), LogXs(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:,1:)
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, l, ij, i0, j0, SizeE
    INTEGER  :: iT, iX, iE1, iE2
    REAL(dp) :: dT, dX, dE1, dE2

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
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iE1, dE1, iE2, dE2, iT, dT, iX, dX, &
    !$OMP          i0, j0, i, j ) &
    !$OMP MAP( to: LogE, LogEs, LogT, LogTs, LogX, LogXs, OS, Table ) &
    !$OMP MAP( from: Interpolant )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) ASYNC( async_flag ) &
    !$ACC PRIVATE( iE1, dE1, iE2, dE2, iT, dT, iX, dX, &
    !$ACC          i0, j0, i, j ), &
    !$ACC COPYIN( LogE, LogEs, LogT, LogTs, LogX, LogXs, OS, Table ) &
    !$ACC COPYOUT( Interpolant )
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

    REAL(dp), INTENT(in)  :: LogE (1:), LogT     , LogX
    REAL(dp), INTENT(in)  :: LogEs(1:), LogTs(1:), LogXs(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:)

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
    ( LogT, LogX, LogTs, LogXs, OS, Table, Interpolant, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogT (1:), LogX (1:)
    REAL(dp), INTENT(in)  :: LogTs(1:), LogXs(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:,1:)
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, k, ij, i0, j0, SizeE
    INTEGER  :: iT, iX
    REAL(dp) :: dT, dX
    REAL(dp) :: p00, p10, p01, p11

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
    !$OMP PRIVATE( iT, dT, iX, dX ) &
    !$OMP MAP( to: LogT, LogTs, LogX, LogXs, OS, Table ) &
    !$OMP MAP( from: Interpolant )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG ASYNC( async_flag ) &
    !$ACC PRIVATE( iT, dT, iX, dX ) &
    !$ACC COPYIN( LogT, LogTs, LogX, LogXs, OS, Table ) &
    !$ACC COPYOUT( Interpolant )
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
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( i, j, p00, p10, p01, p11 )
      DO j = 1, SizeE
        DO i = 1, SizeE
          IF ( i <= j ) THEN
            p00 = Table(i,j,iT  ,iX  )
            p10 = Table(i,j,iT+1,iX  )
            p01 = Table(i,j,iT  ,iX+1)
            p11 = Table(i,j,iT+1,iX+1)
            Interpolant(i,j,k) &
              = 10.0d0 ** (   ( One - dX ) * ( ( One - dT ) * p00 + dT * p10 ) &
                            +         dX   * ( ( One - dT ) * p01 + dT * p11 ) ) - OS
          END IF
        END DO
      END DO
#else
#if defined(WEAKLIB_OACC)
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
        !j = FLOOR( 0.5*( -1.0 + SQRT( 1.0 + 8.0*(ij-1) ) ) + EPSILON(1.0) ) + 1
        !i = (ij-1) - j*(j-1)/2 + 1
        CALL LinearInterp2D_4DArray_2DAligned_Point &
               ( i, j, iT, iX, dT, dX, OS, Table, Interpolant(i,j,k) )
      END DO
#endif
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Aligned


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
    ( LogT, LogX, LogTs, LogXs, OS, Table, Interpolant )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: LogT     , LogX
    REAL(dp), INTENT(in)  :: LogTs(1:), LogXs(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:)

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
    ( LogD, LogT, LogDs, LogTs, Alpha, OS, Table, Interpolant, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogD (1:,1:), LogT (1:)
    REAL(dp), INTENT(in)  :: LogDs(1:)   , LogTs(1:)
    REAL(dp), INTENT(in)  :: Alpha(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:,1:)
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, k, l, ij, i0, j0, SizeE
    INTEGER  :: iD(SIZE(Alpha)), iT
    REAL(dp) :: dD(SIZE(Alpha)), dT
    REAL(dp) :: Interp, SumInterp
    REAL(dp) :: p00, p10, p01, p11

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
    !$OMP PRIVATE( iT, dT, iD, dD ) &
    !$OMP MAP( to: LogT, LogTs, LogD, LogDs, OS, Table, Alpha ) &
    !$OMP MAP( from: Interpolant )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG ASYNC( async_flag ) &
    !$ACC PRIVATE( iT, dT, iD, dD ) &
    !$ACC COPYIN( LogT, LogTs, LogD, LogDs, OS, Table, Alpha ) &
    !$ACC COPYOUT( Interpolant )
#elif defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iT, dT, iD, dD, Interp, SumInterp, &
    !$OMP          i0, j0, i, j )
#endif
    DO k = 1, SIZE( LogT )

      !CALL GetIndexAndDelta_Lin( LogT(k), LogTs, iT, dT )
      !DO l = 1, SIZE( Alpha )
      !  CALL GetIndexAndDelta_Lin( LogD(l,k), LogDs, iD(l), dD(l) )
      !END DO
      iT = Index1D_Lin( LogT(k), LogTs )
      dT = ( LogT(k) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
      DO l = 1, SIZE( Alpha )
        iD(l) = Index1D_Lin( LogD(l,k), LogDs )
        dD(l) = ( LogD(l,k) - LogDs(iD(l)) ) / ( LogDs(iD(l)+1) - LogDs(iD(l)) )
      END DO

#if defined(WEAKLIB_OMP_OL)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( i, j, Interp, SumInterp, &
      !$OMP          p00, p10, p01, p11 )
      DO j = 1, SizeE
        DO i = 1, SizeE
          IF ( i <= j ) THEN
            SumInterp = 0.0d0
            DO l = 1, SIZE( Alpha )
              p00 = Table(i,j,iD(l)  ,iT  )
              p10 = Table(i,j,iD(l)+1,iT  )
              p01 = Table(i,j,iD(l)  ,iT+1)
              p11 = Table(i,j,iD(l)+1,iT+1)
              Interp &
                = 10.0d0 ** (   ( One - dT ) * ( ( One - dD(l) ) * p00 + dD(l) * p10 ) &
                              +         dT   * ( ( One - dD(l) ) * p01 + dD(l) * p11 ) ) - OS
              SumInterp = SumInterp + Alpha(l) * Interp
            END DO
            Interpolant(i,j,k) = SumInterp
          END IF
        END DO
      END DO
#else
#if defined(WEAKLIB_OACC)
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

        SumInterp = 0.0d0
        DO l = 1, SIZE( Alpha )
          CALL LinearInterp2D_4DArray_2DAligned_Point &
                 ( i, j, iD(l), iT, dD(l), dT, OS, Table, Interp )
          SumInterp = SumInterp + Alpha(l) * Interp
        END DO
        Interpolant(i,j,k) = SumInterp
      END DO
#endif
    END DO

  END SUBROUTINE SumLogInterpolateSingleVariable_2D2D_Custom_Aligned


  SUBROUTINE LogInterpolateSingleVariable_3D_Custom &
    ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant, Error_Option )

    REAL(dp), INTENT(in)  :: D (1:), T (1:), Y (1:)
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER :: iP, Error

    Error = 0
    IF( .NOT. ALL( [ SIZE(T), SIZE(Y) ] == SIZE(D) ) )THEN
      Error = 1
      IF( PRESENT( Error_Option ) ) Error_Option = Error
      RETURN
    END IF

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: D, Ds, T, Ts, Y, Ys, OS, Table ) &
    !$OMP MAP( from: Interpolant )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( D, Ds, T, Ts, Y, Ys, OS, Table ) &
    !$ACC COPYOUT( Interpolant )
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

    REAL(dp), INTENT(in)  :: D     , T     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant

    INTEGER  :: iD, iT, iY
    REAL(dp) :: dD, dT, dY
    REAL(dp) :: p000, p100, p010, p110, p001, p101, p011, p111
    INTEGER  :: loD, hiD
    INTEGER  :: loT, hiT
    INTEGER  :: loY, hiY

    !CALL GetIndexAndDelta_Log( D, Ds, iD, dD )
    !CALL GetIndexAndDelta_Log( T, Ts, iT, dT )
    !CALL GetIndexAndDelta_Lin( Y, Ys, iY, dY )

    !CALL LinearInterp3D_3DArray_Point &
    !       ( iD, iT, iY, dD, dT, dY, OS, Table, Interpolant )

    loD = LBOUND(Ds,1)
    hiD = UBOUND(Ds,1)
    loT = LBOUND(Ts,1)
    hiT = UBOUND(Ts,1)
    loY = LBOUND(Ys,1)
    hiY = UBOUND(Ys,1)

    iD = MAX( loD, MIN( hiD-1, loD + FLOOR( (hiD-loD)*LOG10(D/Ds(loD))/LOG10(Ds(hiD)/Ds(loD)) ) ) )
    dD = LOG10( D / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )

    iT = MAX( loT, MIN( hiT-1, loT + FLOOR( (hiT-loT)*LOG10(T/Ts(loT))/LOG10(Ts(hiT)/Ts(loT)) ) ) )
    dT = LOG10( T / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )

    iY = MAX( loY, MIN( hiY-1, loY + FLOOR( (hiY-loY)*(Y-Ys(loY))/(Ys(hiY)-Ys(loY)) ) ) )
    dY = ( Y - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

    p000 = Table(iD  , iT  , iY  )
    p100 = Table(iD+1, iT  , iY  )
    p010 = Table(iD  , iT+1, iY  )
    p110 = Table(iD+1, iT+1, iY  )
    p001 = Table(iD  , iT  , iY+1)
    p101 = Table(iD+1, iT  , iY+1)
    p011 = Table(iD  , iT+1, iY+1)
    p111 = Table(iD+1, iT+1, iY+1)

    Interpolant &
      = 10.0d0 ** (   ( One - dY ) * (   ( One - dT ) * ( ( One - dD ) * p000 + dD * p100 ) &
                                       +         dT   * ( ( One - dD ) * p010 + dD * p110 ) ) &
                    +         dY   * (   ( One - dT ) * ( ( One - dD ) * p001 + dD * p101 ) &
                                       +         dT   * ( ( One - dD ) * p011 + dD * p111 ) ) ) &
        - OS

    !Interpolant &
    !  = 10.0d0**( &
    !      TriLinear &
    !        ( p000, p100, p010, p110, &
    !          p001, p101, p011, p111, &
    !          dY1, dY2, dY3 ) ) - OS

  END SUBROUTINE LogInterpolateSingleVariable_3D_Custom_Point


  SUBROUTINE LogInterpolateSingleVariable_4D_Custom &
      ( LogE, LogD, LogT, Y, LogEs, LogDs, LogTs, Ys, &
        OS, Table, Interpolant, Error_Option )

    REAL(dp), INTENT(in)  :: LogE (1:), LogD (1:), LogT (1:),  Y(1:)
    REAL(dp), INTENT(in)  :: LogEs(1:), LogDs(1:), LogTs(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER :: iP, Error

    Error = 0
    IF( .NOT. ALL( [ SIZE(LogD), SIZE(LogT), SIZE(Y) ] == SIZE(LogE) ) )THEN
      Error = 1
      IF( PRESENT( Error_Option ) ) Error_Option = Error
      RETURN
    END IF

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: LogE, LogEs, LogD, LogDs, LogT, LogTs, Y, Ys, OS, Table ) &
    !$OMP MAP( from: Interpolant )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( LogE, LogEs, LogD, LogDs, LogT, LogTs, Y, Ys, OS, Table ) &
    !$ACC COPYOUT( Interpolant )
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

    REAL(dp), INTENT(in)  :: LogE     , LogD     , LogT     , Y
    REAL(dp), INTENT(in)  :: LogEs(1:), LogDs(1:), LogTs(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
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

    REAL(dp), INTENT(in)  :: D (1:), T (1:), Y (1:)
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:)
    REAL(dp), INTENT(out) :: Derivative(1:,1:)

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

    REAL(dp), INTENT(in)  :: D     , T     , Y
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Ys(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant
    REAL(dp), INTENT(out) :: Derivative(1:)

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
      DerivativeT, DerivativeX, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogE (1:), LogT (1:), LogX (1:)
    REAL(dp), INTENT(in)  :: LogEs(1:), LogTs(1:), LogXs(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:,1:)
    REAL(dp), INTENT(out) :: DerivativeT(1:,1:,1:)
    REAL(dp), INTENT(out) :: DerivativeX(1:,1:,1:)
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, l, ij, i0, j0, SizeE
    INTEGER  :: iT, iX, iE1, iE2
    REAL(dp) :: dT, dX, dE1, dE2
    REAL(dp) :: aT, aX, dI1, dI2

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
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( iE1, dE1, iE2, dE2, iT, dT, aT, iX, dX, aX, dI1, dI2, &
    !$OMP          i0, j0, i, j ) &
    !$OMP MAP( to: LogE, LogEs, LogT, LogTs, LogX, LogXs, OS, Table ) &
    !$OMP MAP( from: Interpolant, DerivativeT, DerivativeX )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) ASYNC( async_flag ) &
    !$ACC PRIVATE( iE1, dE1, iE2, dE2, iT, dT, aT, iX, dX, aX, dI1, dI2, &
    !$ACC          i0, j0, i, j ), &
    !$ACC COPYIN( LogE, LogEs, LogT, LogTs, LogX, LogXs, OS, Table ) &
    !$ACC COPYOUT( Interpolant, DerivativeT, DerivativeX )
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

    REAL(dp), INTENT(in)  :: LogE (1:), LogT     , LogX
    REAL(dp), INTENT(in)  :: LogEs(1:), LogTs(1:), LogXs(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:)
    REAL(dp), INTENT(out) :: DerivativeT(1:,1:)
    REAL(dp), INTENT(out) :: DerivativeX(1:,1:)

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
      DerivativeT, DerivativeX, ASYNC_Option )

    REAL(dp), INTENT(in)  :: LogT (1:), LogX (1:)
    REAL(dp), INTENT(in)  :: LogTs(1:), LogXs(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:,1:)
    REAL(dp), INTENT(out) :: DerivativeT(1:,1:,1:)
    REAL(dp), INTENT(out) :: DerivativeX(1:,1:,1:)
    INTEGER,  INTENT(in), OPTIONAL :: ASYNC_Option

    INTEGER  :: async_flag
    INTEGER  :: i, j, k, ij, i0, j0, SizeE
    INTEGER  :: iT, iX
    REAL(dp) :: dT, aT, dX, aX

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
    !$OMP PRIVATE( iT, dT, aT, iX, dX, aX ) &
    !$OMP MAP( to: LogT, LogTs, LogX, LogXs, OS, Table ) &
    !$OMP MAP( from: DerivativeT, DerivativeX, Interpolant )
#elif defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG ASYNC( async_flag ) &
    !$ACC PRIVATE( iT, dT, aT, iX, dX, aX ) &
    !$ACC COPYIN( LogT, LogTs, LogX, LogXs, OS, Table ) &
    !$ACC COPYOUT( DerivativeT, DerivativeX, Interpolant )
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

    REAL(dp), INTENT(in)  :: LogT     , LogX
    REAL(dp), INTENT(in)  :: LogTs(1:), LogXs(1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: Table(1:,1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant(1:,1:)
    REAL(dp), INTENT(out) :: DerivativeT(1:,1:)
    REAL(dp), INTENT(out) :: DerivativeX(1:,1:)

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

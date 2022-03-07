MODULE wlInterpolationUtilitiesModule

  USE wlKindModule, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Index1D
  PUBLIC :: Index1D_Lin
  PUBLIC :: Index1D_Log

  PUBLIC :: GetIndexAndDelta_Lin
  PUBLIC :: GetIndexAndDelta_Log

  PUBLIC :: Linear
  PUBLIC :: BiLinear
  PUBLIC :: TriLinear
  PUBLIC :: TetraLinear
  PUBLIC :: PentaLinear

  PUBLIC :: LinearInterp_Array_Point
  PUBLIC :: LinearInterp1D_1DArray_Point
  PUBLIC :: LinearInterp2D_2DArray_Point
  PUBLIC :: LinearInterp3D_3DArray_Point
  PUBLIC :: LinearInterp4D_4DArray_Point
  PUBLIC :: LinearInterp5D_5DArray_Point

  PUBLIC :: LinearInterp2D_3DArray_1DAligned_Point
  PUBLIC :: LinearInterp3D_4DArray_1DAligned_Point
  PUBLIC :: LinearInterp4D_5DArray_1DAligned_Point

  PUBLIC :: LinearInterp2D_4DArray_2DAligned_Point
  PUBLIC :: LinearInterp3D_5DArray_2DAligned_Point

  PUBLIC :: LinearInterpDeriv_Array_Point
  PUBLIC :: LinearInterpDeriv3D_3DArray_Point
  PUBLIC :: LinearInterpDeriv4D_4DArray_Point

  PUBLIC :: LinearInterpDeriv2D_4DArray_2DAligned_Point

  REAL(dp), PARAMETER :: One = 1.0_dp

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


END MODULE wlInterpolationUtilitiesModule

MODULE wlInterpolationModule

  USE wlKindModule, ONLY: dp
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  implicit none
  private

  PUBLIC :: LogInterpolateSingleVariable
  PUBLIC :: LogInterpolateAllVariables
  PUBLIC :: LogInterpolateDifferentiateSingleVariable
  PUBLIC :: LogInterpolateDifferentiateAllVariables
  PUBLIC :: MonotonicityCheck
  PUBLIC :: GetGamma1
  PUBLIC :: ComputeTempFromIntEnergy
  PUBLIC :: ComputeTempFromIntEnergy_Lookup
  PUBLIC :: ComputeTempFromIntEnergy_Bisection
  PUBLIC :: ComputeTempFromIntEnergy_Secant
  PUBLIC :: ComputeTempFromEntropy
  PUBLIC :: ComputeTempFromPressure
  PUBLIC :: ComputeTempFromPressure_Bisection
  PUBLIC :: EOSTableQuery
  PUBLIC :: LogInterpolateSingleVariable_1D3D
  PUBLIC :: LogInterpolateSingleVariable_1D3D_Custom
  PUBLIC :: LogInterpolateSingleVariable_2D2D
  PUBLIC :: LogInterpolateSingleVariable_2D2D_Custom

  REAL(dp), PARAMETER :: ln10 = LOG(10.d0)

  INTERFACE LogInterpolateSingleVariable
    MODULE PROCEDURE LogInterpolateSingleVariable_3D
    MODULE PROCEDURE LogInterpolateSingleVariable_3D_Custom
    MODULE PROCEDURE LogInterpolateSingleVariable_4D
    MODULE PROCEDURE LogInterpolateSingleVariable_4D_Custom
  END INTERFACE LogInterpolateSingleVariable

  INTERFACE LogInterpolateAllVariables
    MODULE PROCEDURE LogInterpolateAllVariables_3D
    MODULE PROCEDURE LogInterpolateAllVariables_3D_Custom
  END INTERFACE LogInterpolateAllVariables

  INTERFACE LogInterpolateDifferentiateSingleVariable
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_3D
    MODULE PROCEDURE LogInterpolateDifferentiateSingleVariable_4D
  END INTERFACE LogInterpolateDifferentiateSingleVariable

CONTAINS


  SUBROUTINE locate( xx, n, x, j )

    INTEGER, INTENT(in)      :: n
    INTEGER, INTENT(out)     :: j
    REAL(dp), INTENT(in)     :: x,xx(n)
    INTEGER                  :: jl,jm,ju

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


  PURE INTEGER FUNCTION Index1D( x, xx, n )

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


  PURE INTEGER FUNCTION Index1D_Lin( x, xx, n )

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    Index1D_Lin &
      = FLOOR( 1 + (n-1)*(x-xx(1))/(xx(n)-xx(1)) + 1.d-12 )

    RETURN
  END FUNCTION Index1D_Lin


  PURE INTEGER FUNCTION Index1D_Log( x, xx, n )

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    Index1D_Log &
      = FLOOR( 1 + (n-1)*LOG10(x/xx(1))/LOG10(xx(n)/xx(1)) + 1.d-12 )

    RETURN
  END FUNCTION Index1D_Log


  PURE REAL(dp) FUNCTION TriLinear &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2, dX3 )

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX1, dX2, dX3

    REAL(dp) :: ddX1, ddX2, ddX3

    ddX1 = 1.0_dp - dX1
    ddX2 = 1.0_dp - dX2
    ddX3 = 1.0_dp - dX3

    TriLinear                                        &
      = ddX3                                         &
         * (   ddX2 * ( ddX1 * p000 + dX1 * p100 )   &
             +  dX2 * ( ddX1 * p010 + dX1 * p110 ) ) &
      +  dX3                                         &
         * (   ddX2 * ( ddX1 * p001 + dX1 * p101 )   &
             +  dX2 * ( ddX1 * p011 + dX1 * p111 ) )

    RETURN
  END FUNCTION TriLinear


  PURE REAL(dp) FUNCTION TetraLinear &
    ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3, dX4 )

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
   
    END DO

  END SUBROUTINE LogInterpolateSingleVariable_3D

  
  SUBROUTINE LogInterpolateSingleVariable_1D3D &
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
    REAL(dp), DIMENSION(:,:),     INTENT(out) :: Interpolant

    INTEGER :: &
      i, j, il1, il2, il3, il4
    REAL(dp), DIMENSION(4) :: &
      alpha, delta
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
               ( LogX1, LogX2, LogX3, LinX4, LogCoordsX1, LogCoordsX2, &
                 LogCoordsX3, LinCoordsX4, Offset, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LinX4
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LinCoordsX4
    REAL(dp),                     INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:,:),     INTENT(out) :: Interpolant

    INTEGER :: i, j
    INTEGER :: iX1, iX2, iX3, iX4
    INTEGER :: p1, p2, p3, p4
    REAL(dp), DIMENSION(4) :: dX
    REAL(dp), DIMENSION(0:1,0:1,0:1,0:1) :: p

    DO j = 1, SIZE( LogX2 )

      iX4   = Index1D_Lin( LinX4(j), LinCoordsX4, SIZE( LinCoordsX4 ) )
      dX(4) = ( LinX4(j) - LinCoordsX4(iX4) ) &
              / ( LinCoordsX4(iX4+1) - LinCoordsX4(iX4) )

      iX3   = Index1D_Lin( LogX3(j), LogCoordsX3, SIZE( LogCoordsX3 ) )
      dX(3) = ( LogX3(j) - LogCoordsX3(iX3) ) &
              / ( LogCoordsX3(iX3+1) - LogCoordsX3(iX3) )

      iX2   = Index1D_Lin( LogX2(j), LogCoordsX2, SIZE( LogCoordsX2 ) )
      dX(2) = ( LogX2(j) - LogCoordsX2(iX2) ) &
              / ( LogCoordsX2(iX2+1) - LogCoordsX2(iX2) )

      DO i = 1, SIZE( LogX1 )

        iX1   = Index1D_Lin( LogX1(i), LogCoordsX1, SIZE( LogCoordsX1 ) )
        dX(1) = ( LogX1(i) - LogCoordsX1(iX1) ) &
                / ( LogCoordsX1(iX1+1) - LogCoordsX1(iX1) )

        DO p4 = 0, 1
          DO p3 = 0, 1
            DO p2 = 0, 1
              DO p1 = 0, 1

                p(p1,p2,p3,p4) &
                  = TABLE(iX1+p1,iX2+p2,iX3+p3,iX4+p4)

              END DO
            END DO
          END DO
        END DO

        Interpolant(i,j) &
          = TetraLinear &
              ( p(0,0,0,0), p(1,0,0,0), p(0,1,0,0), p(1,1,0,0), &
                p(0,0,1,0), p(1,0,1,0), p(0,1,1,0), p(1,1,1,0), &
                p(0,0,0,1), p(1,0,0,1), p(0,1,0,1), p(1,1,0,1), &
                p(0,0,1,1), p(1,0,1,1), p(0,1,1,1), p(1,1,1,1), &
                dX(1), dX(2), dX(3), dX(4) )

      END DO ! i
    END DO ! j

    Interpolant(:,:) &
      = 10**( Interpolant(:,:) ) - Offset

  END SUBROUTINE LogInterpolateSingleVariable_1D3D_Custom


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
      i, j, k, il1, il2, il3, il4, &
      SizeC1, SizeC2, SizeC3, SizeC4, &
      SizeX1, SizeX2, SizeX3
    REAL(dp), DIMENSION(4) :: &
      delta
    REAL(dp) :: &
      p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, &
      p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111

    SizeX1 = SIZE( x1 )
    SizeX2 = SIZE( x2 )
    SizeX3 = SIZE( x3 )

    SizeC1 = SIZE( Coordinate1 )
    SizeC2 = SIZE( Coordinate2 )
    SizeC3 = SIZE( Coordinate3 )
    SizeC4 = SIZE( Coordinate4 )

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

        DO i = 1, SizeX1

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
        DO i = 1, SizeX1

          Interpolant(i,j,k) &
            = 10.d0**( Interpolant(i,j,k) ) - Offset

        END DO ! i
      END DO ! j
    END DO ! k

  END SUBROUTINE LogInterpolateSingleVariable_2D2D


  SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom &
               ( LogX1, LogX2, LogX3, LogX4, LogCoordsX1, LogCoordsX2, &
                 LogCoordsX3, LogCoordsX4, Offset, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX4
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX4
    REAL(dp),                     INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:,:,:),   INTENT(out) :: Interpolant

    INTEGER :: i, j, k
    INTEGER :: iX1, iX2, iX3, iX4
    INTEGER :: p1, p2, p3, p4
    REAL(dp), DIMENSION(4) :: dX
    REAL(dp), DIMENSION(0:1,0:1,0:1,0:1) :: p

    DO k = 1, SIZE( LogX3 )

      iX4 &
        = Index1D_Lin( LogX4(k), LogCoordsX4, SIZE( LogCoordsX4 ) )
      dX(4) &
        = ( LogX4(k) - LogCoordsX4(iX4) ) &
            / ( LogCoordsX4(iX4+1) - LogCoordsX4(iX4) )

      iX3 &
        = Index1D_Lin( LogX3(k), LogCoordsX3, SIZE( LogCoordsX3 ) )
      dX(3) &
        = ( LogX3(k) - LogCoordsX3(iX3) ) &
            / ( LogCoordsX3(iX3+1) - LogCoordsX3(iX3) )

      DO j = 1, SIZE( LogX2 )

        iX2 &
          = Index1D_Lin( LogX2(j), LogCoordsX2, SIZE( LogCoordsX2 ) )
        dX(2) &
          = ( LogX2(j) - LogCoordsX2(iX2) ) &
              / ( LogCoordsX2(iX2+1) - LogCoordsX2(iX2) )

        DO i = 1, SIZE( LogX1 )

          iX1 &
            = Index1D_Lin( LogX1(i), LogCoordsX1, SIZE( LogCoordsX1 ) )
          dX(1) &
            = ( LogX1(i) - LogCoordsX1(iX1) ) &
                / ( LogCoordsX1(iX1+1) - LogCoordsX1(iX1) )

          DO p4 = 0, 1
            DO p3 = 0, 1
              DO p2 = 0, 1
                DO p1 = 0, 1

                  p(p1,p2,p3,p4) &
                    = TABLE(iX1+p1,iX2+p2,iX3+p3,iX4+p4)

                END DO
              END DO
            END DO
          END DO

          Interpolant(i,j,k) &
            = TetraLinear &
                ( p(0,0,0,0), p(1,0,0,0), p(0,1,0,0), p(1,1,0,0), &
                  p(0,0,1,0), p(1,0,1,0), p(0,1,1,0), p(1,1,1,0), &
                  p(0,0,0,1), p(1,0,0,1), p(0,1,0,1), p(1,1,0,1), &
                  p(0,0,1,1), p(1,0,1,1), p(0,1,1,1), p(1,1,1,1), &
                  dX(1), dX(2), dX(3), dX(4) )

        END DO ! i
      END DO ! j
    END DO ! k

    Interpolant(:,:,:) &
      = 10**( Interpolant(:,:,:) ) - Offset

  END SUBROUTINE LogInterpolateSingleVariable_2D2D_Custom


  SUBROUTINE LogInterpolateSingleVariable_3D_Custom &
               ( D, T, Y, Ds, Ts, Ys, OS, Table, Interpolant )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: D,  T,  Y
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Ds, Ts, Ys
    REAL(dp),                   INTENT(in)  :: OS
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),     INTENT(out) :: Interpolant

    INTEGER  :: &
      iP, iD, iT, iY
    REAL(dp) :: &
      dD, dT, dY, &
      p000, p100, p010, p110, &
      p001, p101, p011, p111

    IF( .NOT. ALL( [ SIZE(T), SIZE(Y) ] == SIZE(D) ) )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'LogInterpolateSingleVariable_3D_Custom'
      WRITE(*,'(A4,A)') &
        '', 'ERROR: arrays of interpolation points have different sizes'
      WRITE(*,*)
      RETURN
    END IF

    DO iP = 1, SIZE( D )

      iD = Index1D( D(iP), Ds, SIZE( Ds ) )
      iT = Index1D( T(iP), Ts, SIZE( Ts ) )
      iY = Index1D( Y(iP), Ys, SIZE( Ys ) )

      dD = LOG10( D(iP) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
      dT = LOG10( T(iP) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
      dY = ( Y(iP) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

      p000 = ( Table( iD  , iT  , iY   ) )
      p100 = ( Table( iD+1, iT  , iY   ) )
      p010 = ( Table( iD  , iT+1, iY   ) )
      p110 = ( Table( iD+1, iT+1, iY   ) )
      p001 = ( Table( iD  , iT  , iY+1 ) )
      p101 = ( Table( iD+1, iT  , iY+1 ) )
      p011 = ( Table( iD  , iT+1, iY+1 ) )
      p111 = ( Table( iD+1, iT+1, iY+1 ) )

      Interpolant(iP) &
        = 10.0d0**( &
            TriLinear &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, dD, dT, dY ) ) - OS

    END DO

  END SUBROUTINE LogInterpolateSingleVariable_3D_Custom


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

    WRITE(*,*) 'TetraLinear   '
  END SUBROUTINE LogInterpolateSingleVariable_4D_Custom


  SUBROUTINE LogInterpolateAllVariables_3D &
               ( x1, x2, x3, LogInterp, TS, DV, Interpolants, MaskVar )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp 
    TYPE(ThermoStateType), INTENT(in) :: TS
    TYPE(DependentVariablesType), INTENT(in) :: DV
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: MaskVar

    REAL(dp), DIMENSION(:,:), INTENT(out) :: Interpolants 

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: delta
    INTEGER :: i, j, dim1, dim2, dim3
    INTEGER, DIMENSION(:), ALLOCATABLE :: il1, il2, il3
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: work_mask
    INTEGER                             :: Masksize

    epsilon = 1.d-200

    Masksize = SIZE( x2 )
    dim1     = SIZE( TS % States(1) % Values )
    dim2     = SIZE( TS % States(2) % Values )
    dim3     = SIZE( TS % States(3) % Values )

    ALLOCATE( work_mask( Masksize ), delta( 3, Masksize ), &
              il1( Masksize), il2( Masksize ), il3( Masksize ) )

    IF ( PRESENT(MaskVar) ) THEN
      work_mask = MaskVar
    ELSE
      work_mask = .true.
    END IF

    DO i = 1, Masksize

      IF ( .not.work_mask(i) ) CYCLE

      ASSOCIATE( Coordinate1 => TS % States(1) % Values, &    
                 Coordinate2 => TS % States(2) % Values, &    
                 Coordinate3 => TS % States(3) % Values )   

      CALL locate( Coordinate1, dim1, x1(i), il1(i) )
      CALL locate( Coordinate2, dim2, x2(i), il2(i) )
      CALL locate( Coordinate3, dim3, x3(i), il3(i) )

      IF ( LogInterp(1) == 1 ) THEN
        delta(1,i) = LOG10( x1(i) / Coordinate1(il1(i)) ) / LOG10( Coordinate1(il1(i)+1) / Coordinate1(il1(i)) )
      ELSE
        delta(1,i) = ( x1(i) - Coordinate1(il1(i)) ) / ( Coordinate1(il1(i)+1) - Coordinate1(il1(i)) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        delta(2,i) = LOG10( x2(i) / Coordinate2(il2(i)) ) / LOG10( Coordinate2(il2(i)+1) / Coordinate2(il2(i)) )
      ELSE
        delta(2,i) = ( x2(i) - Coordinate2(il2(i)) ) / ( Coordinate2(il2(i)+1) - Coordinate2(il2(i)) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        delta(3,i) = LOG10( x3(i) / Coordinate3(il3(i)) ) / LOG10( Coordinate3(il3(i)+1) / Coordinate3(il3(i)) )
      ELSE
        delta(3,i) = ( x3(i) - Coordinate3(il3(i)) ) / ( Coordinate3(il3(i)+1) - Coordinate3(il3(i)) )
      END IF

      END ASSOCIATE

    END DO ! i = 1, Masksize

    DO j = 1, DV % nVariables

      ASSOCIATE( Table => DV % Variables(j) % Values(:,:,:), &
                 Offset => DV % Offsets(j) )
      DO i = 1, Masksize

        IF ( .not.work_mask(i) ) CYCLE

        p000 = ( Table( il1(i)  , il2(i)  , il3(i)   ) )
        p100 = ( Table( il1(i)+1, il2(i)  , il3(i)   ) )
        p010 = ( Table( il1(i)  , il2(i)+1, il3(i)   ) )
        p110 = ( Table( il1(i)+1, il2(i)+1, il3(i)   ) )
        p001 = ( Table( il1(i)  , il2(i)  , il3(i)+1 ) )
        p101 = ( Table( il1(i)+1, il2(i)  , il3(i)+1 ) )
        p011 = ( Table( il1(i)  , il2(i)+1, il3(i)+1 ) )
        p111 = ( Table( il1(i)+1, il2(i)+1, il3(i)+1 ) )

        Interpolants(i,j) &
          = 10.d0**( &
                (1.0_dp - delta(3,i)) * ( (1.0_dp - delta(1,i)) * (1.0_dp - delta(2,i)) * p000   &
                                       +            delta(1,i)  * (1.0_dp - delta(2,i)) * p100   &
                                       + ( 1.0_dp - delta(1,i)) *           delta(2,i)  * p010   &
                                       +            delta(1,i)  *           delta(2,i)  * p110 ) &
                        + delta(3,i)  * ( (1.0_dp - delta(1,i)) * (1.0_dp - delta(2,i)) * p001   &
                                       +            delta(1,i)  * (1.0_dp - delta(2,i)) * p101   &
                                       +  (1.0_dp - delta(1,i)) *           delta(2,i)  * p011   &
                                       +            delta(1,i)  *           delta(2,i)  * p111 ) &

                   ) - Offset

        END DO ! i = 1, Masksize

      END ASSOCIATE

    END DO ! j = 1, nVariables
        

  END SUBROUTINE LogInterpolateAllVariables_3D


  SUBROUTINE LogInterpolateAllVariables_3D_Custom &
               ( D, T, Y, Ds, Ts, Ys, DV, Interpolants )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: D,  T,  Y
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Ds, Ts, Ys
    TYPE(DependentVariablesType), INTENT(in)  :: DV
    REAL(dp), DIMENSION(:,:),     INTENT(out) :: Interpolants

    INTEGER :: &
      iP, iV, iD, iT, iY
    REAL(dp) :: &
      dD, dT, dY, &
      p000, p100, p010, p110, &
      p001, p101, p011, p111

    IF( .NOT. ALL( [ SIZE(T), SIZE(Y) ] == SIZE(D) ) )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'LogInterpolateAllVariables_3D_Custom'
      WRITE(*,'(A4,A)') &
        '', 'ERROR: arrays of interpolation points have different sizes'
      WRITE(*,*)
      RETURN
    END IF

    DO iP = 1, SIZE( D )

      iD = Index1D( D(iP), Ds, SIZE( Ds ) )
      iT = Index1D( T(iP), Ts, SIZE( Ts ) )
      iY = Index1D( Y(iP), Ys, SIZE( Ys ) )

      dD = LOG10( D(iP) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )
      dT = LOG10( T(iP) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )
      dY = ( Y(iP) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )

      DO iV = 1, DV % nVariables

        ASSOCIATE &
          ( Table => DV % Variables(iV) % Values(:,:,:), &
            OS    => DV % Offsets  (iV) )

        p000 = Table( iD  , iT  , iY   )
        p100 = Table( iD+1, iT  , iY   )
        p010 = Table( iD  , iT+1, iY   )
        p110 = Table( iD+1, iT+1, iY   )
        p001 = Table( iD  , iT  , iY+1 )
        p101 = Table( iD+1, iT  , iY+1 )
        p011 = Table( iD  , iT+1, iY+1 )
        p111 = Table( iD+1, iT+1, iY+1 )

        Interpolants(iV, iP) &
          = 10.0d0**( &
              TriLinear &
                ( p000, p100, p010, p110, &
                  p001, p101, p011, p111, dD, dT, dY ) ) - OS

        END ASSOCIATE ! Table, etc.

      END DO

    END DO

  END SUBROUTINE LogInterpolateAllVariables_3D_Custom


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D                    &
               ( x1, x2, x3, Coordinate1, Coordinate2, Coordinate3,          &
                 LogInterp, Offset, Table, Interpolant, Derivative, MaskVar)     

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

  SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D &
               ( x1, x2, x3, x4, Coordinate1, Coordinate2, Coordinate3, &
                 Coordinate4, LogInterp, Offset, Table, Interpolant, &
                 Derivative, debug )     

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


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D_Custom &
       ( E, D, T, Y, Es, Ds, Ts, Ys, OS, Table, Interpolant, Derivatives )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: E,  D,  T,  Y
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Es, Ds, Ts, Ys
    REAL(dp),                     INTENT(in)  :: OS
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),       INTENT(out) :: Interpolant
    REAL(dp), DIMENSION(:,:),     INTENT(out) :: Derivatives

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_4D_Custom


  SUBROUTINE LogInterpolateDifferentiateAllVariables &
               ( x1, x2, x3, LogInterp, TS, DV, Interpolants, Derivatives, MaskVar )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp 
    TYPE(ThermoStateType), INTENT(in) :: TS
    TYPE(DependentVariablesType), INTENT(in) :: DV
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: MaskVar

    REAL(dp), DIMENSION(:,:), INTENT(out) :: Interpolants 

    REAL(dp), DIMENSION(:,:,:), INTENT(out) :: Derivatives 
    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(:,:), ALLOCATABLe :: alpha, delta
    INTEGER :: i, j, dim1, dim2, dim3
    INTEGER, DIMENSION(:), ALLOCATABLE :: il1, il2, il3
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: work_mask
    INTEGER                             :: Masksize

    epsilon = 1.d-200

    Masksize = SIZE( x2 )
    dim1     = SIZE( TS % States(1) % Values )
    dim2     = SIZE( TS % States(2) % Values )
    dim3     = SIZE( TS % States(3) % Values )

    ALLOCATE( work_mask( Masksize ), alpha( 3, Masksize ), delta( 3, Masksize ), &
              il1( Masksize), il2( Masksize ), il3( Masksize ) )

    IF ( PRESENT(MaskVar) ) THEN
      work_mask = MaskVar
    ELSE
      work_mask = .true.
    END IF

    DO i = 1, Masksize

      IF ( .not.work_mask(i) ) CYCLE

      ASSOCIATE( Coordinate1 => TS % States(1) % Values, &
                 Coordinate2 => TS % States(2) % Values, &
                 Coordinate3 => TS % States(3) % Values )

      CALL locate( Coordinate1, dim1, x1(i), il1(i) )
      CALL locate( Coordinate2, dim2, x2(i), il2(i) )
      CALL locate( Coordinate3, dim3, x3(i), il3(i) )

      IF ( LogInterp(1) == 1 ) THEN
        alpha(1,i) = ( 1.0d0 ) / ( x1(i) * LOG10( Coordinate1(il1(i)+1) / Coordinate1(il1(i)) ) )
        delta(1,i) = LOG10( x1(i) / Coordinate1(il1(i)) ) / LOG10( Coordinate1(il1(i)+1) / Coordinate1(il1(i)) )
      ELSE
        alpha(1,i) = ( ln10 ) / ( Coordinate1(il1(i)+1) - Coordinate1(il1(i)) )
        delta(1,i) = ( x1(i) - Coordinate1(il1(i)) ) / ( Coordinate1(il1(i)+1) - Coordinate1(il1(i)) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        alpha(2,i) = ( 1.0d0 ) / ( x2(i) * LOG10( Coordinate2(il2(i)+1) / Coordinate2(il2(i)) ) )
        delta(2,i) = LOG10( x2(i) / Coordinate2(il2(i)) ) / LOG10( Coordinate2(il2(i)+1) / Coordinate2(il2(i)) )
      ELSE
        alpha(2,i) = ( ln10 ) / ( Coordinate2(il2(i)+1) - Coordinate2(il2(i)) )
        delta(2,i) = ( x2(i) - Coordinate2(il2(i)) ) / ( Coordinate2(il2(i)+1) - Coordinate2(il2(i)) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        alpha(3,i) = ( 1.0d0 ) / ( x3(i) * LOG10( Coordinate3(il3(i)+1) / Coordinate3(il3(i)) ) )
        delta(3,i) = LOG10( x3(i) / Coordinate3(il3(i)) ) / LOG10( Coordinate3(il3(i)+1) / Coordinate3(il3(i)) )
      ELSE
        alpha(3,i) = ( ln10 ) / ( Coordinate3(il3(i)+1) - Coordinate3(il3(i)) )
        delta(3,i) = ( x3(i) - Coordinate3(il3(i)) ) / ( Coordinate3(il3(i)+1) - Coordinate3(il3(i)) )
      END IF

      END ASSOCIATE

    END DO ! i = 1, Masksize

      DO j = 1, DV % nVariables

        ASSOCIATE( Table => DV % Variables(j) % Values(:,:,:), &
                   Offset => DV % Offsets(j) )
        DO i = 1, Masksize

          IF ( .not.work_mask(i) ) CYCLE

          p000 = ( Table( il1(i)  , il2(i)  , il3(i)   ) )
          p100 = ( Table( il1(i)+1, il2(i)  , il3(i)   ) )
          p010 = ( Table( il1(i)  , il2(i)+1, il3(i)   ) )
          p110 = ( Table( il1(i)+1, il2(i)+1, il3(i)   ) )
          p001 = ( Table( il1(i)  , il2(i)  , il3(i)+1 ) )
          p101 = ( Table( il1(i)+1, il2(i)  , il3(i)+1 ) )
          p011 = ( Table( il1(i)  , il2(i)+1, il3(i)+1 ) )
          p111 = ( Table( il1(i)+1, il2(i)+1, il3(i)+1 ) )

          Interpolants(i,j) &
            = 10.d0**( &
                  (1.0_dp - delta(3,i)) * ( (1.0_dp - delta(1,i)) * (1.0_dp - delta(2,i)) * p000   &
                                       +            delta(1,i)  * (1.0_dp - delta(2,i)) * p100   &
                                       + ( 1.0_dp - delta(1,i)) *           delta(2,i)  * p010   &
                                       +            delta(1,i)  *           delta(2,i)  * p110 ) &
                          + delta(3,i)  * ( (1.0_dp - delta(1,i)) * (1.0_dp - delta(2,i)) * p001   &
                                       +            delta(1,i)  * (1.0_dp - delta(2,i)) * p101   &
                                       +  (1.0_dp - delta(1,i)) *           delta(2,i)  * p011   &
                                       +            delta(1,i)  *           delta(2,i)  * p111 ) &
  
                     ) - Offset
  
          Derivatives(i,1,j) &
            = ( (Interpolants(i,j) ) * alpha(1,i) &
                * ( (1.0_dp - delta(3,i)) * ( (delta(2,i) - 1.0_dp) * p000   &
                                        +  ( 1.0_dp - delta(2,i)) * p100   &
                                        -             delta(2,i)  * p010   &
                                        +             delta(2,i)  * p110 ) &
                             + delta(3,i) * ( (delta(2,i) - 1.0_dp) * p001   &
                                        +  ( 1.0_dp - delta(2,i)) * p101   &
                                        -             delta(2,i)  * p011   &
                                        +             delta(2,i)  * p111 ) ) )
    
          Derivatives(i,2,j) &
            = ( ( Interpolants(i,j) ) * alpha(2,i) &
                * ( (1.0_dp - delta(3,i) ) * ( (delta(1,i) - 1.0_dp) * p000   &
                                         -             delta(1,i)  * p100   & 
                                         +  ( 1.0_dp - delta(1,i)) * p010   & 
                                         +             delta(1,i)  * p110 ) & 
                              + delta(3,i) * ( (delta(1,i) - 1.0_dp) * p001   & 
                                         -             delta(1,i)  * p101   & 
                                         +   (1.0_dp - delta(1,i)) * p011   & 
                                         +             delta(1,i)  * p111 ) ) )
  
          Derivatives(i,3,j) &
            = ( ( Interpolants(i,j) ) * alpha(3,i) &
                                       * ( ( (delta(1,i) - 1.0_dp)) * (1.0_dp - delta(2,i)) * p000   &
                                           -            delta(1,i)  * (1.0_dp - delta(2,i)) * p100   &
                                           - ( 1.0_dp - delta(1,i)) *           delta(2,i)  * p010   &
                                           -            delta(1,i)  *           delta(2,i)  * p110   &
                                           +  (1.0_dp - delta(1,i)) * (1.0_dp - delta(2,i)) * p001   &
                                           +            delta(1,i)  * (1.0_dp - delta(2,i)) * p101   &
                                           +  (1.0_dp - delta(1,i)) *           delta(2,i)  * p011   &
                                           +            delta(1,i)  *           delta(2,i)  * p111 ) )
          END DO ! i = 1, Masksize
  
        END ASSOCIATE

      END DO ! j = 1, nVariables

  END SUBROUTINE LogInterpolateDifferentiateAllVariables 

  SUBROUTINE GetGamma1( x1, x2, x3, Coordinate1, Coordinate2, &
                Coordinate3, LogInterp, TS, DV, Gamma1 ) 

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    INTEGER, DIMENSION(3), INTENT(in)  :: LogInterp
    TYPE(ThermoStateType), INTENT(in) :: TS
    TYPE(DependentVariablesType), INTENT(in) :: DV

    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: Gamma1

    INTEGER :: i
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Interpolant 
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative 

    ALLOCATE( Interpolant( SIZE(x1) ) )
    ALLOCATE( Gamma1( SIZE(x1) ) )
    ALLOCATE( Derivative( SIZE(x1), 3 ) )

      CALL LogInterpolateDifferentiateSingleVariable( x1, x2, x3,                 &
                                    TS % States(1) % Values(:),        &
                                    TS % States(2) % Values(:),        &
                                    TS % States(3) % Values(:),        &
                                    LogInterp, DV % Offsets(1),        &
                                    DV % Variables(1) % Values(:,:,:), &
                                    Interpolant(:), Derivative(:,:) )
       
      DO i = 1, SIZE(x1) 
        Gamma1(i) =  ( x1(i)/Interpolant(i) ) * Derivative(i, 1 ) 
      END DO

  END SUBROUTINE GetGamma1 
  
  SUBROUTINE MonotonicityCheck ( Table, Nrho, NT, NYe, Axis, Repaired )

    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    INTEGER, DIMENSION(:,:,:), INTENT(in) :: Repaired
    INTEGER, INTENT(in) :: Nrho 
    INTEGER, INTENT(in) :: NT
    INTEGER, INTENT(in) :: NYe
    INTEGER, INTENT(in) :: Axis

    INTEGER :: i, j, k, count

    97 FORMAT ("Table not monotonic in rho at (Nrho, NT, NYe) = ", 3(1x,i4) )
    98 FORMAT ("Table not monotonic in T at (Nrho, NT, NYe) = ", 3(1x,i4) )
    99 FORMAT ("Table not monotonic in Ye at (Nrho, NT, NYe) = ", 3(1x,i4) )
 
    count = 0
    
    SELECT CASE ( Axis ) 
     
    CASE( 1 )
      DO k = 1, NYe
        DO j = 1, NT  
          DO i = 2, Nrho - 1

            IF ( ( ( Table(i+1, j, k) - Table(i, j, k) ) * &
                 ( Table(i, j, k) - Table(i-1, j, k) ) ) < 0. ) THEN
              WRITE (*,97) i, j, k
              WRITE (*,*) "Repaired =", Repaired(i,j,k), Repaired(i+1,j,k), Repaired(i-1,j,k), &
                Repaired(i,j+1,k), Repaired(i,j-1,k), Repaired(i,j,k+1), Repaired(i,j,k-1)
              count = count + 1
            END IF
          END DO
        END DO
      END DO

    CASE( 2 )
      DO k = 1, NYe
        DO j = 2, NT - 1 
          DO i = 1, Nrho

            IF ( ( ( Table(i, j+1, k) - Table(i, j, k) ) * &
                 ( Table(i, j, k) - Table(i, j-1, k) ) ) < 0.) THEN 
              WRITE (*,98) i, j, k
              WRITE (*,*) "Repaired =", Repaired(i,j,k), Repaired(i+1,j,k), Repaired(i-1,j,k), &
                Repaired(i,j+1,k), Repaired(i,j-1,k), Repaired(i,j,k+1), Repaired(i,j,k-1)
              count = count + 1
            END IF
          END DO
        END DO
      END DO

   CASE( 3 )
      DO k = 2, NYe - 1
        DO j = 1, NT
          DO i = 1, Nrho

            IF ( ( ( Table(i, j, k+1) - Table(i, j, k) ) * &
                 ( Table(i, j, k) - Table(i, j, k-1) ) ) < 0. ) &
              WRITE (*, 99) i, j, k
              WRITE (*,*) "Repaired =", Repaired(i,j,k), Repaired(i+1,j,k), Repaired(i-1,j,k), &
                Repaired(i,j+1,k), Repaired(i,j-1,k), Repaired(i,j,k+1), Repaired(i,j,k-1)
              count = count + 1
          END DO
        END DO
      END DO

    CASE DEFAULT
      WRITE (*,*) "Invalid Axis", Axis
      STOP

    END SELECT
    WRITE (*,*) count, " Non-monotonic out of " , NYe*NT*Nrho

  END SUBROUTINE MonotonicityCheck 


  SUBROUTINE ComputeTempFromIntEnergy &
               ( rho, e_int, ye, density_table, temp_table, ye_table, &
                 LogInterp, energy_table, Offset, Temperature )

    REAL(dp), INTENT(in)                    :: rho
    REAL(dp), INTENT(in)                    :: e_int
    REAL(dp), INTENT(in)                    :: ye
    REAL(dp), DIMENSION(:), INTENT(in)      :: density_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: ye_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: temp_table
    INTEGER, DIMENSION(3), INTENT(in)       :: LogInterp 
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: energy_table
    REAL(dp), INTENT(in) :: Offset

    REAL(dp), DIMENSION(1)                  :: eibuff ! internal energy buffer
    INTEGER                                 :: nPoints
    INTEGER                                 :: i, j
    REAL(dp), DIMENSION(:), ALLOCATABLE :: rhobuff 
    REAL(dp), DIMENSION(:), ALLOCATABLE :: energy_array
    REAL(dp), DIMENSION(:), ALLOCATABLE :: yebuff

    REAL(dp), DIMENSION(1), INTENT(out)     :: Temperature

  nPoints = SIZE(temp_table)

  ALLOCATE( energy_array( nPoints ), rhobuff( nPoints ), yebuff( nPoints) )

  rhobuff(1:nPoints) = rho
  eibuff(1) = e_int
  yebuff(1:nPoints) = ye

  CALL LogInterpolateSingleVariable                                     &
         ( rhobuff, temp_table, yebuff,                                 &
           density_table,                                               &
           temp_table,                                                  &
           ye_table,                                                    &
           LogInterp,                                                   &
           Offset,                                                      &
           energy_table(:,:,:), energy_array )

    DO j = 1, SIZE( eibuff )

      CALL locate( energy_array, SIZE( energy_array ), eibuff(j), i )
      IF ( i == SIZE(energy_array) ) THEN
        STOP 'energy too high'
      END IF

      IF ( i == 0 ) THEN
        Temperature(j) = 0.d0
        CYCLE
      END IF
      Temperature(j) = 10.d0**( &
             LOG10( temp_table(i) ) + LOG10( temp_table(i+1) / temp_table(i) ) &
                        * LOG10( ( ( eibuff(j) + Offset ) / ( energy_array(i) + Offset ) ) )             &
                        / LOG10( ( energy_array(i+1) + Offset ) / ( energy_array(i) + Offset ) ) )
    END DO

    DEALLOCATE( energy_array, rhobuff, yebuff )

  END SUBROUTINE ComputeTempFromIntEnergy


  SUBROUTINE ComputeTempFromIntEnergy_Lookup &
    ( D, E, Y, D_table, T_table, Y_table, LogInterp, E_table, Offset, T )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: D
    REAL(dp), DIMENSION(:),     INTENT(in)  :: E
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Y
    REAL(dp), DIMENSION(:),     INTENT(in)  :: D_table
    REAL(dp), DIMENSION(:),     INTENT(in)  :: T_table
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Y_table
    INTEGER,  DIMENSION(3),     INTENT(in)  :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: E_table
    REAL(dp),                   INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:),     INTENT(out) :: T

    LOGICAL, PARAMETER :: Debug = .TRUE.
    INTEGER  :: i, Error
    INTEGER  :: ilD, ilY, ilE, ilE1, ilE2, ilE3, ilE4, ilEa, ilEb
    INTEGER  :: nPtsE
    REAL(dp) :: logE, tmpE
    REAL(dp), DIMENSION(:), ALLOCATABLE :: ptsD, ptsT, ptsY, ptsE

    DO i = 1, SIZE( D )

      Error = 0

      ilD = Index1D( D(i), D_table, SIZE( D_table ) )
      ilY = Index1D( Y(i), Y_table, SIZE( Y_table ) )

      logE = LOG10( E(i) + Offset )
      ilE1 = Index1D( logE, E_table(ilD,  :,ilY  ), SIZE( T_Table ) )
      ilE2 = Index1D( logE, E_table(ilD+1,:,ilY  ), SIZE( T_Table ) )
      ilE3 = Index1D( logE, E_table(ilD,  :,ilY+1), SIZE( T_Table ) )
      ilE4 = Index1D( logE, E_table(ilD+1,:,ilY+1), SIZE( T_Table ) )

      ilEa = MAX( MINVAL( [ilE1,ilE2,ilE3,ilE4] ) - 1, 1 )
      ilEb = MIN( MAXVAL( [ilE1,ilE2,ilE3,ilE4] ) + 2, SIZE( T_Table ) )

      nPtsE = SIZE( T_Table(ilEa:ilEb) )
      ALLOCATE( ptsD(nPtsE), ptsT(nPtsE), ptsY(nPtsE), ptsE(nPtsE) )
      ptsD = D(i)
      ptsT = T_Table(ilEa:ilEb)
      ptsY = Y(i)

      CALL LogInterpolateSingleVariable &
             ( ptsD, ptsT, ptsY, D_Table(ilD:ilD+1), T_Table(ilEa:ilEb), &
               Y_Table(ilY:ilY+1), LogInterp, Offset, &
               E_Table(ilD:ilD+1,ilEa:ilEb,ilY:ilY+1), ptsE )

      tmpE = E(i)

      IF( (ptsE(1)-E(i))*(ptsE(nPtsE)-E(i)) > 0.0_DP )THEN

        IF( Debug )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A)') &
            '', 'Warning: ComputeTempFromIntEnergy_Lookup'
          WRITE(*,'(A6,A20,ES10.4E2,A9,ES10.4E2)') &
            '', 'No Root Between T = ', ptsT(1), ' and T = ', ptsT(nPtsE)
          WRITE(*,*)
          WRITE(*,*) '  nPtsE = ', nPtsE
          WRITE(*,*) '  ptsE  = ', ptsE
          WRITE(*,*) '    ia  = ', ilEa
          WRITE(*,*) '    ib  = ', ilEb
          WRITE(*,*) '    Ta  = ', T_Table(ilEa)
          WRITE(*,*) '    Tb  = ', T_Table(ilEb)
          WRITE(*,*) '    Ea  = ', ptsE(1)
          WRITE(*,*) '    Eb  = ', ptsE(nPtsE)
          WRITE(*,*) '     i  = ', i
          WRITE(*,*) '     E  = ', E(i)
          WRITE(*,*) '     D  = ', D(i)
          WRITE(*,*) '     Y  = ', Y(i)
          WRITE(*,*)
          STOP
        END IF

        ! --- Reset Energy Density ---

        IF( E(i) < ptsE(1) .AND. E(i) < ptsE(nPtsE) )THEN

          tmpE = 0.5 * ( ptsE(1) + ptsE(2) )

        ELSE

          tmpE = 0.5 * ( ptsE(nPtsE-1) + ptsE(nPtsE) )

        END IF

        Error = 1

      END IF

      ilE = Index1D( tmpE, ptsE, nPtsE )

      T(i) &
        = 10.d0**( LOG10( ptsT(ilE) ) &
                   + LOG10( ptsT(ilE+1)/ptsT(ilE) ) &
                     * LOG10( (tmpE+Offset)/(ptsE(ilE)+Offset) ) &
                       / LOG10( (ptsE(ilE+1)+Offset)/(ptsE(ilE)+Offset) ) )

      IF( Error == 1 .AND. Debug )THEN
        WRITE(*,*)
        WRITE(*,*) '  T = ', T(i)
        WRITE(*,*)
      END IF

      DEALLOCATE( ptsD, ptsT, ptsY, ptsE )

    END DO

  END SUBROUTINE ComputeTempFromIntEnergy_Lookup


  SUBROUTINE ComputeTempFromIntEnergy_Bisection &
    ( D, E, Y, D_table, T_table, Y_table, LogInterp, E_table, Offset, T )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: D
    REAL(dp), DIMENSION(:),     INTENT(in)  :: E
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Y
    REAL(dp), DIMENSION(:),     INTENT(in)  :: D_table
    REAL(dp), DIMENSION(:),     INTENT(in)  :: T_table
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Y_table
    INTEGER,  DIMENSION(3),     INTENT(in)  :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: E_table
    REAL(dp),                   INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:),     INTENT(out) :: T

    LOGICAL  :: Converged
    INTEGER  :: i, Iter
    INTEGER,  PARAMETER :: MaxIter = 128
    REAL(DP), PARAMETER :: Tol = 1.0d-10
    REAL(dp), DIMENSION(1) :: a, b, c, ab, E_a, E_b, E_c, f_a, f_b, f_c

    DO i = 1, SIZE( D )

      a = T_table(1)
      CALL LogInterpolateSingleVariable &
             ( [ D(i) ] , a, [ Y(i) ], D_table, T_table, Y_table, &
               LogInterp, Offset, E_table, E_a )
      f_a = E(i) - E_a

      b = T_table(SIZE(T_table))
      CALL LogInterpolateSingleVariable &
             ( [ D(i) ] , b, [ Y(i) ], D_table, T_table, Y_table, &
               LogInterp, Offset, E_table, E_b )
      f_b = E(i) - E_b

      IF( ALL( f_a*f_b > 0.0_dp ) )THEN
        WRITE(*,*)
        WRITE(*,'(A4,A)') &
          '', 'Error: ComputeTempFromIntEnergy_Bisection'
        WRITE(*,'(A6,A20,ES10.4E2,A9,ES10.4E2)') &
          '', 'No Root Between T = ', a, ' and T = ', b
        WRITE(*,*)
        STOP
      END IF

      ab = b - a

      Converged = .FALSE.
      Iter      = 0
      DO WHILE ( .NOT. Converged )

        Iter = Iter + 1

        ab = 0.5_dp * ab
        c  = a + ab

        CALL LogInterpolateSingleVariable &
               ( [ D(i) ] , c, [ Y(i) ], D_table, T_table, Y_table, &
                 LogInterp, Offset, E_table, E_c )

        f_c = E(i) - E_c

        IF( ALL( f_a * f_c < 0.0_dp ) )THEN

          b   = c
          f_b = f_c

        ELSE

          a   = c
          f_a = f_c

        END IF

        IF( ALL( ABS( f_c ) / E(i) < Tol ) ) &
          Converged = .TRUE.

        IF( Iter > MaxIter )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A)') &
            '', 'ComputeTempFromIntEnergy_Bisection'
          WRITE(*,'(A6,A21,I4.4,A11)') &
            '', 'No Convergence After ', Iter, ' Iterations'
          WRITE(*,*)
        END IF

      END DO

      T(i) = c(1)

    END DO

  END SUBROUTINE ComputeTempFromIntEnergy_Bisection


  SUBROUTINE ComputeTempFromIntEnergy_Secant &
    ( D, E, Y, D_table, T_table, Y_table, LogInterp, E_table, Offset, T )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: D
    REAL(dp), DIMENSION(:),     INTENT(in)  :: E
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Y
    REAL(dp), DIMENSION(:),     INTENT(in)  :: D_table
    REAL(dp), DIMENSION(:),     INTENT(in)  :: T_table
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Y_table
    INTEGER,  DIMENSION(3),     INTENT(in)  :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: E_table
    REAL(dp),                   INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:),     INTENT(out) :: T

    LOGICAL  :: Converged
    INTEGER :: i, Iter
    INTEGER,  PARAMETER :: MaxIter = 128
    REAL(DP), PARAMETER :: Tol = 1.0d-10
    REAL(dp), DIMENSION(1) :: T_L, T_H, a, b, c, E_a, E_b, f_a, f_b

    T_L = T_table(1)
    T_H = T_table(SIZE(T_table))

    DO i = 1, SIZE( D )

      a = T_L
      CALL LogInterpolateSingleVariable &
             ( [ D(i) ] , a, [ Y(i) ], D_table, T_table, Y_table, &
               LogInterp, Offset, E_table, E_a )
      f_a = E(i) - E_a

      b = T_H
      CALL LogInterpolateSingleVariable &
             ( [ D(i) ] , b, [ Y(i) ], D_table, T_table, Y_table, &
               LogInterp, Offset, E_table, E_b )
      f_b = E(i) - E_b

      IF( ALL( f_a*f_b > 0.0_dp ) )THEN
        WRITE(*,*)
        WRITE(*,'(A4,A)') &
          '', 'Error: ComputeTempFromIntEnergy_Secant'
        WRITE(*,'(A6,A20,ES10.4E2,A9,ES10.4E2)') &
          '', 'No Root Between T = ', a, ' and T = ', b
        WRITE(*,*)
        STOP
      END IF

      Converged = .FALSE.
      Iter      = 0
      DO WHILE ( .NOT. Converged )

        Iter = Iter + 1

        IF( ALL( ABS( f_a ) > ABS( f_b ) ) )THEN
          ! -- Swap (a,b)
          a = a + b
          b = a - b
          a = a - b
          ! -- Swap (f_a,f_b)
          f_a = f_a + f_b
          f_b = f_a - f_b
          f_a = f_a - f_b
        END IF

        c = ( b * f_a - a * f_b ) / ( f_a - f_b )
        b = a; f_b = f_a

        a = c
        CALL LogInterpolateSingleVariable &
               ( [ D(i) ] , a, [ Y(i) ], D_table, T_table, Y_table, &
                 LogInterp, Offset, E_table, E_a )
        f_a = E(i) - E_a

        IF( ALL( ABS( f_a ) / E(i) < Tol ) ) &
          Converged = .TRUE.

        IF( Iter > MaxIter )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A)') &
            '', 'ComputeTempFromIntEnergy_Secant'
          WRITE(*,'(A6,A21,I4.4,A11)') &
            '', 'No Convergence After ', Iter, ' Iterations'
          WRITE(*,*)
        END IF

      END DO

      T(i) = c(1)

    END DO

  END SUBROUTINE ComputeTempFromIntEnergy_Secant


  SUBROUTINE ComputeTempFromEntropy &
               ( rho, s, ye, density_table, temp_table, ye_table, &
                 LogInterp, entropy_table, Offset, Temperature )

    REAL(dp), INTENT(in)                    :: rho
    REAL(dp), INTENT(in)                    :: s
    REAL(dp), INTENT(in)                    :: ye
    REAL(dp), DIMENSION(:), INTENT(in)      :: density_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: ye_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: temp_table
    INTEGER, DIMENSION(3), INTENT(in)       :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: entropy_table
    REAL(dp), INTENT(in) :: Offset

    REAL(dp), DIMENSION(1)                  :: sbuff ! entropy buffer
    INTEGER                                 :: nPoints
    INTEGER                                 :: i, j
    REAL(dp), DIMENSION(:), ALLOCATABLE :: rhobuff
    REAL(dp), DIMENSION(:), ALLOCATABLE :: entropy_array
    REAL(dp), DIMENSION(:), ALLOCATABLE :: yebuff

    REAL(dp), DIMENSION(1), INTENT(out)     :: Temperature

  nPoints = SIZE(temp_table)

  ALLOCATE( entropy_array( nPoints ), rhobuff( nPoints ), yebuff( nPoints) )

  rhobuff(1:nPoints) = rho
  sbuff(1) = s
  yebuff(1:nPoints) = ye

  CALL LogInterpolateSingleVariable                                     &
         ( rhobuff, temp_table, yebuff,                                 &
           density_table,                                               &
           temp_table,                                                  &
           ye_table,                                                    &
           LogInterp,                                                   &
           Offset,                                                      &
           entropy_table(:,:,:), entropy_array )

    DO j = 1, SIZE( sbuff )

      CALL locate( entropy_array, SIZE( entropy_array ), sbuff(j), i )
      IF ( i == SIZE(entropy_array) ) THEN
        STOP
      END IF

      IF ( i == 0 ) THEN
        Temperature(j) = 0.d0
        CYCLE
      END IF
      Temperature(j) = 10.d0**( &
             LOG10( temp_table(i) ) + LOG10( temp_table(i+1) / temp_table(i) ) &
                        * LOG10( ( sbuff(j) / entropy_array(i) ) )             &
                        / LOG10( entropy_array(i+1) / entropy_array(i) ) )
    END DO

  DEALLOCATE( entropy_array, rhobuff, yebuff )

  END SUBROUTINE ComputeTempFromEntropy

  SUBROUTINE ComputeTempFromPressure &
               ( rho, p, ye, density_table, temp_table, ye_table, &
                 LogInterp, pressure_table, Offset, Temperature )

    REAL(dp), INTENT(in)                    :: rho
    REAL(dp), INTENT(in)                    :: p
    REAL(dp), INTENT(in)                    :: ye
    REAL(dp), DIMENSION(:), INTENT(in)      :: density_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: ye_table
    REAL(dp), DIMENSION(:), INTENT(in)      :: temp_table
    INTEGER, DIMENSION(3), INTENT(in)       :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: pressure_table
    REAL(dp), INTENT(in) :: Offset

    REAL(dp), DIMENSION(1)                  :: pbuff ! pressure buffer
    INTEGER                                 :: nPoints
    INTEGER                                 :: i, j
    REAL(dp), DIMENSION(:), ALLOCATABLE :: rhobuff
    REAL(dp), DIMENSION(:), ALLOCATABLE :: pressure_array
    REAL(dp), DIMENSION(:), ALLOCATABLE :: yebuff

    REAL(dp), DIMENSION(1), INTENT(out)     :: Temperature

  nPoints = SIZE(temp_table)

  ALLOCATE( pressure_array( nPoints ), rhobuff( nPoints ), yebuff( nPoints) )

  rhobuff(1:nPoints) = rho
  pbuff(1) = p
  yebuff(1:nPoints) = ye

  CALL LogInterpolateSingleVariable                                     &
         ( rhobuff, temp_table, yebuff,                                 &
           density_table,                                               &
           temp_table,                                                  &
           ye_table,                                                    &
           LogInterp,                                                   &
           Offset,                                                      &
           pressure_table(:,:,:), pressure_array )

    DO j = 1, SIZE( pbuff )

      CALL locate( pressure_array, SIZE( pressure_array ), pbuff(j), i )
      IF ( i == SIZE(pressure_array) ) THEN
        STOP
      END IF

      IF ( i == 0 ) THEN
        Temperature(j) = 0.d0
        CYCLE
      END IF
      Temperature(j) = 10.d0**( &
             LOG10( temp_table(i) ) + LOG10( temp_table(i+1) / temp_table(i) ) &
                        * LOG10( ( pbuff(j) / pressure_array(i) ) )             &
                        / LOG10( pressure_array(i+1) / pressure_array(i) ) )
    END DO

  DEALLOCATE( pressure_array, rhobuff, yebuff )

  END SUBROUTINE ComputeTempFromPressure


  SUBROUTINE ComputeTempFromPressure_Bisection &
    ( D, P, Y, D_table, T_table, Y_table, LogInterp, P_table, Offset, T )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: D
    REAL(dp), DIMENSION(:),     INTENT(in)  :: P
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Y
    REAL(dp), DIMENSION(:),     INTENT(in)  :: D_table
    REAL(dp), DIMENSION(:),     INTENT(in)  :: T_table
    REAL(dp), DIMENSION(:),     INTENT(in)  :: Y_table
    INTEGER,  DIMENSION(3),     INTENT(in)  :: LogInterp
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: P_table
    REAL(dp),                   INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:),     INTENT(out) :: T

    LOGICAL  :: Converged
    INTEGER  :: i, Iter
    INTEGER,  PARAMETER :: MaxIter = 128
    REAL(DP), PARAMETER :: Tol = 1.0d-10
    REAL(dp), DIMENSION(1) :: a, b, c, ab, P_a, P_b, P_c, f_a, f_b, f_c

    DO i = 1, SIZE( D )

      a = T_table(1)
      CALL LogInterpolateSingleVariable &
             ( [ D(i) ] , a, [ Y(i) ], D_table, T_table, Y_table, &
               LogInterp, Offset, P_table, P_a )
      f_a = P(i) - P_a

      b = T_table(SIZE(T_table))
      CALL LogInterpolateSingleVariable &
             ( [ D(i) ] , b, [ Y(i) ], D_table, T_table, Y_table, &
               LogInterp, Offset, P_table, P_b )
      f_b = P(i) - P_b

      IF( ALL( f_a*f_b > 0.0_dp ) )THEN
        WRITE(*,*)
        WRITE(*,'(A4,A)') &
          '', 'Warning: ComputeTempFromPressure_Bisection'
        WRITE(*,'(A6,A27,ES10.4E2,A9,ES10.4E2)') &
          '', 'No Unique Root Between T = ', a, ' and T = ', b
        WRITE(*,'(A6,A6,ES10.4E2)') '', 'P_a = ', P_a
        WRITE(*,'(A6,A6,ES10.4E2)') '', 'P_i = ', P(i)
        WRITE(*,'(A6,A6,ES10.4E2)') '', 'P_b = ', P_b
        WRITE(*,*)
      END IF

      ab = b - a

      Converged = .FALSE.
      Iter      = 0
      DO WHILE ( .NOT. Converged )

        Iter = Iter + 1

        ab = 0.5_dp * ab
        c  = a + ab

        CALL LogInterpolateSingleVariable &
               ( [ D(i) ] , c, [ Y(i) ], D_table, T_table, Y_table, &
                 LogInterp, Offset, P_table, P_c )

        f_c = P(i) - P_c

        IF( ALL( f_a * f_c < 0.0_dp ) )THEN

          b   = c
          f_b = f_c

        ELSE

          a   = c
          f_a = f_c

        END IF

        IF( ALL( ABS( f_c ) / P(i) < Tol ) .OR. ALL( ab / a < Tol ) ) &
          Converged = .TRUE.

        IF( Iter > MaxIter )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A)') &
            '', 'ComputeTempFromPressure_Bisection'
          WRITE(*,'(A6,A21,I4.4,A11)') &
            '', 'No Convergence After ', Iter, ' Iterations'
          WRITE(*,*)
        END IF

      END DO

      T(i) = c(1)

    END DO

  END SUBROUTINE ComputeTempFromPressure_Bisection


  SUBROUTINE EOSTableQuery &
               ( rho, T, Ye, LogInterp, TS, DV, Interpolants )

    INTEGER                                    :: i, j
    REAL(dp), DIMENSION(:), INTENT(in)         :: rho
    REAL(dp), DIMENSION(:), INTENT(in)         :: T
    REAL(dp), DIMENSION(:), INTENT(in)         :: Ye
    INTEGER, DIMENSION(3), INTENT(in)          :: LogInterp
    TYPE(ThermoStateType), INTENT(in)          :: TS
    TYPE(DependentVariablesType), INTENT(in)   :: DV
    REAL(dp), DIMENSION(:,:), INTENT(out)      :: Interpolants

!    CALL InitializeHDF( )

!    CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

    CALL LogInterpolateAllVariables( rho, T, Ye, LogInterp, &
                                     TS, DV, Interpolants )
    DO i = 1, SIZE(rho)
      WRITE(*,*) 'Rho=', rho(i), 'T=', T(i), 'Ye=', Ye(i)
      DO j = 1, DV % nVariables
        WRITE(*,*) DV % Names(j), Interpolants(i,j)
      END DO
    END DO

  END SUBROUTINE EOSTableQuery

END MODULE wlInterpolationModule

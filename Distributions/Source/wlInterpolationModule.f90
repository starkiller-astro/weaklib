MODULE wlInterpolationModule

  USE wlKindModule, ONLY: dp

  implicit none

  PUBLIC LogInterpolateSingleVariable
  PUBLIC locate 

CONTAINS

  SUBROUTINE locate( xx, n, x, j )

    INTEGER, INTENT(in)      :: n
    INTEGER, INTENT(out)     :: j
    REAL(dp), INTENT(in)     ::  x,xx(n)
    INTEGER                  ::  jl,jm,ju

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

  SUBROUTINE LogInterpolateSingleVariable( x1, x2, x3, Coordinate1, Coordinate2, &
                                           Coordinate3, LogInterp, Table, Interpolant )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    LOGICAL, DIMENSION(3), INTENT(in)  :: LogInterp 
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table

    REAL(dp), DIMENSION(:), INTENT(out) :: Interpolant 

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111
    REAL(dp), DIMENSION(3) :: delta
    INTEGER :: il1, il2, il3
  
    DO i = 1, SIZE(x1)  
  
      CALL Locate( Coordinate1, x1(i), SIZE(Coordinate1), il1 ) 
      CALL Locate( Coordinate2, x2(i), SIZE(Coordinate2), il2 )
      CALL Locate( Coordinate3, x3(i), SIZE(Coordinate3), il3 )
    
      p000 = Table( il1  , il2  , il3   )
      p100 = Table( il1+1, il2  , il3   )
      p010 = Table( il1  , il2+1, il3   )
      p110 = Table( il1+1, il2+1, il3   )
      p001 = Table( il1  , il2  , il3+1 )
      p101 = Table( il1+1, il2  , il3+1 )
      p011 = Table( il1  , il2+1, il3+1 )
      p111 = Table( il1+1, il2+1, il3+1 )
  
      IF LogInterp = (/.true.,.true.,.false./) THEN

      delta(1) = LOG10( x1(i) - Coordinate1(il1) ) / LOG10( Coordinate1(il1+1) - Coordinate1(il1) )
      delta(2) = LOG10( x2(i) - Coordinate2(il2) ) / LOG10( Coordinate2(il2+1) - Coordinate2(il2) )
      delta(3) = ( x3(i) - Coordinate3(il3) ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
    
      ELSEIF LogInterp = (/.true.,.false.,.true./) THEN

      delta(1) = LOG10( x1(i) - Coordinate1(il1) ) / LOG10( Coordinate1(il1+1) - Coordinate1(il1) )
      delta(2) = ( x2(i) - Coordinate2(il2) ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      delta(3) = LOG10( x3(i) - Coordinate3(il3) ) / LOG10( Coordinate3(il3+1) - Coordinate3(il3) )
  
      ELSEIF LogInterp = (/.false.,.true.,.true./) THEN

      delta(1) = ( x1(i) - Coordinate1(il1) ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      delta(2) = LOG10( x2(i) - Coordinate2(il2) ) / LOG10( Coordinate2(il2+1) - Coordinate2(il2) )
      delta(3) = LOG10( x3(i) - Coordinate3(il3) ) / LOG10( Coordinate3(il3+1) - Coordinate3(il3) )
  
      END IF

      Interpolant(i) &
        = 10.d0**( &
            + (1.0_dp - delta(3)) * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &                
                                   +            delta(1)  * (1.0_dp - delta(2)) * p100   &
                                   + ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                                   +            delta(1)  *           delta(2)  * p110 ) &
                        delta(3)  * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                                   +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                                   +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                                   +            delta(1)  *           delta(2)  * p111 ) &
 
                 ) 
    END DO 

  END SUBROUTINE LogInterpolateSingleVariable

END MODULE wlInterpolationModule


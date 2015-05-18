MODULE wlInterpolationModule

  USE wlKindModule, ONLY: dp

  implicit none

  PUBLIC LogInterpolateSingleVariable
  PUBLIC locate 
  PUBLIC MonotonicityCheck

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

  SUBROUTINE LogInterpolateFine1D( i, j, k, iMinGradient, delta, Table, Interpolant )

    INTEGER, INTENT(in) :: i, j, k
    REAL(dp), INTENT(in) :: delta
    INTEGER, INTENT(in) :: iMinGradient
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    
    REAL(dp), INTENT(out) :: Interpolant 

    REAL(dp) :: p0, p1, epsilon, offset 
   
    epsilon = 1.d-200
 
    SELECT CASE( iMinGradient ) 

      CASE(1) 
        offset = -2.d0*MIN( 0.d0, Table( i-1, j, k ), Table( i+1, j, k ) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i-1, j, k ) + offset + epsilon  )
        p1 = LOG10( Table( i+1, j, k ) + offset + epsilon )

      CASE(2)
        offset = -2.d0*MIN( 0.d0, Table( i, j-1, k ), Table( i, j+1, k ) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i, j-1, k ) + offset + epsilon  )
        p1 = LOG10( Table( i, j+1, k ) + offset + epsilon  )

      CASE(3)
        offset = -2.d0*MIN( 0.d0, Table( i, j, k-1 ), Table( i, j, k+1 ) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i, j, k-1 ) + offset + epsilon  )
        p1 = LOG10( Table( i, j, k+1 ) + offset + epsilon  )
        
    END SELECT

    Interpolant = 10.d0**( delta * p1 + ( 1.d0 - delta ) * p0 ) - offset 
    WRITE (*,*) "Interpolant=", Interpolant 
    WRITE (*,*) "Delta=", delta 

  END SUBROUTINE LogInterpolateFine1D

  SUBROUTINE LogInterpolateCoarse1D( i, j, k, iMinGradient, iLimits, delta, Table, Interpolant )

    INTEGER, INTENT(in) :: i, j, k
    REAL(dp), INTENT(in) :: delta
    INTEGER, INTENT(in) :: iMinGradient
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    INTEGER, DIMENSION(:,:,:,:), INTENT(in) :: iLimits
    
    REAL(dp), INTENT(out) :: Interpolant 

    REAL(dp) :: p0, p1, offset, epsilon
    
    epsilon = 1.d-200

    SELECT CASE( iMinGradient ) 

      CASE(1) 
        offset = -2.d0*MIN( 0.d0, Table(iLimits(1,i,j,k ),j,k), Table(iLimits(2,i,j,k),j,k) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( iLimits(1,i,j,k), j, k ) + offset + epsilon  )
        p1 = LOG10( Table( iLimits(2,i,j,k), j, k ) + offset + epsilon  )

      CASE(2)
        offset = -2.d0*MIN( 0.d0, Table(i,iLimits(1,i,j,k ),k), Table(i,iLimits(2,i,j,k),k) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i, iLimits(1,i,j,k), k ) + offset + epsilon  )
        p1 = LOG10( Table( i, iLimits(2,i,j,k), k ) + offset + epsilon  )

      CASE(3)
        offset = -2.d0*MIN( 0.d0, Table(i,j,iLimits(1,i,j,k )), Table(i,j,iLimits(2,i,j,k)) )
        WRITE (*,*) "Offset=", offset
        p0 = LOG10( Table( i, j, iLimits(1,i,j,k) ) + offset + epsilon  )
        p1 = LOG10( Table( i, j, iLimits(2,i,j,k) ) + offset + epsilon  )
        
    END SELECT
    WRITE (*,*) "p0, p1 =", p0, p1
    Interpolant = 10.d0**( delta * p1 + ( 1.d0 - delta ) * p0 ) - offset 
    WRITE (*,*) "Interpolant=", Interpolant 
    WRITE (*,*) "Delta=", delta 

  END SUBROUTINE LogInterpolateCoarse1D

  SUBROUTINE LogInterpolateSingleVariable( x1, x2, x3, Coordinate1, Coordinate2, &
                                           Coordinate3, LogInterp, Offset, Table, Interpolant )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    LOGICAL, DIMENSION(3), INTENT(in)  :: LogInterp 
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    REAL(dp), INTENT(in) :: Offset
    
    REAL(dp), DIMENSION(:), INTENT(out) :: Interpolant 

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: delta
    INTEGER :: i, il1, il2, il3
  
    epsilon = 1.d-200

    DO i = 1, SIZE(x1)  
  
      CALL locate( Coordinate1, SIZE(Coordinate1), x1(i), il1 ) 
      CALL locate( Coordinate2, SIZE(Coordinate2), x2(i), il2 )
      CALL locate( Coordinate3, SIZE(Coordinate3), x3(i), il3 )

    !    WRITE (*,*) "Offset=", Offset

      p000 = ( Table( il1  , il2  , il3   ) )
      p100 = ( Table( il1+1, il2  , il3   ) )
      p010 = ( Table( il1  , il2+1, il3   ) )
      p110 = ( Table( il1+1, il2+1, il3   ) )
      p001 = ( Table( il1  , il2  , il3+1 ) )
      p101 = ( Table( il1+1, il2  , il3+1 ) )
      p011 = ( Table( il1  , il2+1, il3+1 ) )
      p111 = ( Table( il1+1, il2+1, il3+1 ) )

     ! WRITE (*,*) "p000 =", p000

      IF ( LogInterp(1) ) THEN 
      delta(1) = LOG10( x1(i) / Coordinate1(il1) ) / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
      delta(1) = ( x1(i) - Coordinate1(il1) ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) ) THEN 
      delta(2) = LOG10( x2(i) / Coordinate2(il2) ) / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
      delta(2) = ( x2(i) - Coordinate2(il2) ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) ) THEN 
      delta(3) = LOG10( x3(i) / Coordinate3(il3) ) / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
      delta(3) = ( x3(i) - Coordinate3(il3) ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF
     ! WRITE (*,*) "Deltas = ", delta
      Interpolant(i) &
        = 10.d0**( &
              (1.0_dp - delta(3)) * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &                
                                   +            delta(1)  * (1.0_dp - delta(2)) * p100   &
                                   + ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                                   +            delta(1)  *           delta(2)  * p110 ) &
                      + delta(3)  * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                                   +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                                   +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                                   +            delta(1)  *           delta(2)  * p111 ) &
 
                 ) - Offset 
    END DO 

  END SUBROUTINE LogInterpolateSingleVariable

  SUBROUTINE LoneCellLogInterpolateSingleVariable( x1, x2, x3, Coordinate1, Coordinate2, &
                                           Coordinate3, LogInterp, Table, Interpolant )

    INTEGER, INTENT(in) :: x1
    INTEGER, INTENT(in) :: x2
    INTEGER, INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    LOGICAL, DIMENSION(3), INTENT(in)  :: LogInterp 
    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    

    REAL(dp), INTENT(out) :: Interpolant 

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111
    REAL(dp), DIMENSION(3) :: delta
    INTEGER :: il1, il2, il3
    REAL(dp) :: PreInterpolant 
 

      il1 = x1 - 1
      il2 = x2 - 1
      il3 = x3 - 1

      p000 = LOG10( Table( il1  , il2  , il3   ) )
      p100 = LOG10( Table( il1+2, il2  , il3   ) )
      p010 = LOG10( Table( il1  , il2+2, il3   ) )
      p110 = LOG10( Table( il1+2, il2+2, il3   ) )
      p001 = LOG10( Table( il1  , il2  , il3+2 ) )
      p101 = LOG10( Table( il1+2, il2  , il3+2 ) )
      p011 = LOG10( Table( il1  , il2+2, il3+2 ) )
      p111 = LOG10( Table( il1+2, il2+2, il3+2 ) )

      WRITE (*,*) p000

      IF ( LogInterp(1) ) THEN 
      delta(1) = LOG10( Coordinate1(x1) / Coordinate1(il1) ) / LOG10( Coordinate1(il1+2) / Coordinate1(il1) )
      ELSE
      delta(1) = ( Coordinate1(x1) - Coordinate1(il1) ) / ( Coordinate1(il1+2) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) ) THEN 
      delta(2) = LOG10( Coordinate2(x2) / Coordinate2(il2) ) / LOG10( Coordinate2(il2+2) / Coordinate2(il2) )
      ELSE
      delta(2) = ( Coordinate2(x2) - Coordinate2(il2) ) / ( Coordinate2(il2+2) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) ) THEN 
      delta(3) = LOG10( Coordinate3(x3) / Coordinate3(il3) ) / LOG10( Coordinate3(il3+2) / Coordinate3(il3) )
      ELSE
      delta(3) = ( Coordinate3(x3) - Coordinate3(il3) ) / ( Coordinate3(il3+2) - Coordinate3(il3) )
      END IF

WRITE (*,*) delta
      Interpolant &
        = 10.d0**( &
           (  (1.0_dp - delta(3)) * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p000   &          
                                   +            delta(1)  * (1.0_dp - delta(2)) * p100   &
                                   + ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                                   +            delta(1)  *           delta(2)  * p110 ) &
                      + delta(3)  * ( (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                                   +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                                   +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                                   +            delta(1)  *           delta(2)  * p111 ) &
                 ) ) 
WRITE (*,*) Interpolant

  END SUBROUTINE LoneCellLogInterpolateSingleVariable
  
  SUBROUTINE MonotonicityCheck ( Table, Nrho, NT, NYe, Axis, Repaired )

    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    LOGICAL, DIMENSION(:,:,:), INTENT(in) :: Repaired
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

          END DO
        END DO
      END DO

    CASE DEFAULT
      WRITE (*,*) "Invalid Axis", Axis
      STOP

    END SELECT
    WRITE (*,*) count, " Non-monotonic out of " , NYe*NT*Nrho

  END SUBROUTINE MonotonicityCheck 

END MODULE wlInterpolationModule


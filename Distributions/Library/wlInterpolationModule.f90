MODULE wlInterpolationModule

  USE wlKindModule, ONLY: dp
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  implicit none

  PUBLIC LogInterpolateSingleVariable
  PUBLIC LogInterpolateAllVariables
  PUBLIC LogInterpolateDifferentiateSingleVariable
  PUBLIC LogInterpolateDifferentiateAllVariables
  PUBLIC locate 
  PUBLIC MonotonicityCheck
  PUBLIC GetGamma1

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
    INTEGER :: i, j, k, il1, il2, il3
  
    epsilon = 1.d-200
  
    ! Check the usage of x2
    ! x2 is used to size the run of the do loop because the input is expected 
    ! to be a series of 3-tuple rho, T, Ye points; so SIZE(x1) = SIZE(x2) = SIZE(x3) 
    !-------------------------------

    DO i = 1, SIZE(x2) 
 
      CALL locate( Coordinate1, SIZE(Coordinate1), x1(i), il1 )
      CALL locate( Coordinate2, SIZE(Coordinate2), x2(i), il2 )
      CALL locate( Coordinate3, SIZE(Coordinate3), x3(i), il3 )

      p000 = ( Table( il1  , il2  , il3   ) )
      p100 = ( Table( il1+1, il2  , il3   ) )
      p010 = ( Table( il1  , il2+1, il3   ) )
      p110 = ( Table( il1+1, il2+1, il3   ) )
      p001 = ( Table( il1  , il2  , il3+1 ) )
      p101 = ( Table( il1+1, il2  , il3+1 ) )
      p011 = ( Table( il1  , il2+1, il3+1 ) )
      p111 = ( Table( il1+1, il2+1, il3+1 ) )

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

  SUBROUTINE LogInterpolateAllVariables( x1, x2, x3, LogInterp, TS, DV, Interpolants )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    LOGICAL, DIMENSION(3), INTENT(in)  :: LogInterp 
    TYPE(ThermoStateType), INTENT(in) :: TS
    TYPE(DependentVariablesType), INTENT(in) :: DV

    REAL(dp), DIMENSION(:,:), INTENT(out) :: Interpolants 

    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: delta
    INTEGER :: i, j, il1, il2, il3
  
    epsilon = 1.d-200

    DO i = 1, SIZE(x2)

      ASSOCIATE( Coordinate1 => TS % States(1) % Values, &    
                 Coordinate2 => TS % States(2) % Values, &    
                 Coordinate3 => TS % States(3) % Values )   

      CALL locate( Coordinate1, SIZE(Coordinate1), x1(i), il1 )
      CALL locate( Coordinate2, SIZE(Coordinate2), x2(i), il2 )
      CALL locate( Coordinate3, SIZE(Coordinate3), x3(i), il3 )


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

      END ASSOCIATE

      DO j = 1, DV % nVariables

      
        ASSOCIATE( Table => DV % Variables(j) % Values(:,:,:), &
                   Offset => DV % Offsets(j) )
      
        p000 = ( Table( il1  , il2  , il3   ) )
        p100 = ( Table( il1+1, il2  , il3   ) )
        p010 = ( Table( il1  , il2+1, il3   ) )
        p110 = ( Table( il1+1, il2+1, il3   ) )
        p001 = ( Table( il1  , il2  , il3+1 ) )
        p101 = ( Table( il1+1, il2  , il3+1 ) )
        p011 = ( Table( il1  , il2+1, il3+1 ) )
        p111 = ( Table( il1+1, il2+1, il3+1 ) )
  
        Interpolants(i,j) &
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
        
        END ASSOCIATE
    
      END DO
    END DO

  END SUBROUTINE LogInterpolateAllVariables 

  SUBROUTINE LogInterpolateDifferentiateSingleVariable( x1, x2, x3, Coordinate1, Coordinate2, &
                                           Coordinate3, LogInterp, Offset, Table,        &
                                           Interpolant, Derivative )     

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
    REAL(dp), DIMENSION(:,:), INTENT(out) :: Derivative 
    
    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: alpha, delta
    INTEGER :: i, j, k, il1, il2, il3

    epsilon = 1.d-200


    DO i = 1, SIZE(x2)

      CALL locate( Coordinate1, SIZE(Coordinate1), x1(i), il1 )
      CALL locate( Coordinate2, SIZE(Coordinate2), x2(i), il2 )
      CALL locate( Coordinate3, SIZE(Coordinate3), x3(i), il3 )

      p000 = ( Table( il1  , il2  , il3   ) )
      p100 = ( Table( il1+1, il2  , il3   ) )
      p010 = ( Table( il1  , il2+1, il3   ) )
      p110 = ( Table( il1+1, il2+1, il3   ) )
      p001 = ( Table( il1  , il2  , il3+1 ) )
      p101 = ( Table( il1+1, il2  , il3+1 ) )
      p011 = ( Table( il1  , il2+1, il3+1 ) )
      p111 = ( Table( il1+1, il2+1, il3+1 ) )

      IF ( LogInterp(1) ) THEN
      alpha(1) = ( 1.0d0 ) / ( x1(i) * LOG10( Coordinate1(il1+1) / Coordinate1(il1) ) )
      delta(1) = LOG10( x1(i) / Coordinate1(il1) ) / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
      alpha(1) = ( 1.0d0 ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      delta(1) = ( x1(i) - Coordinate1(il1) ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) ) THEN
      alpha(2) = ( 1.0d0 ) / ( x2(i) * LOG10( Coordinate2(il2+1) / Coordinate2(il2) ) )
      delta(2) = LOG10( x2(i) / Coordinate2(il2) ) / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
      alpha(2) = ( 1.0d0 ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      delta(2) = ( x2(i) - Coordinate2(il2) ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) ) THEN
      alpha(3) = ( 1.0d0 ) / ( x3(i) * LOG10( Coordinate3(il3+1) / Coordinate3(il3) ) )
      delta(3) = LOG10( x3(i) / Coordinate3(il3) ) / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
      alpha(3) = ( 1.0d0 ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
      delta(3) = ( x3(i) - Coordinate3(il3) ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

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
              * ( (1.0_dp - delta(3) ) * ( (1.0_dp - delta(1)) * p000   &
                                       +             delta(1)  * p100   &
                                       -  ( 1.0_dp - delta(1)) * p010   &
                                       -             delta(1)  * p110 ) &
                            + delta(3) * ( (1.0_dp - delta(1)) * p001   &
                                       +             delta(1)  * p101   &
                                       -   (1.0_dp - delta(1)) * p011   &
                                       -             delta(1)  * p111 ) ) )

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

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable

  SUBROUTINE LogInterpolateDifferentiateAllVariables( x1, x2, x3, LogInterp, TS, DV, Interpolants, Derivatives )

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    LOGICAL, DIMENSION(3), INTENT(in)  :: LogInterp 
    TYPE(ThermoStateType), INTENT(in) :: TS
    TYPE(DependentVariablesType), INTENT(in) :: DV

    REAL(dp), DIMENSION(:,:), INTENT(out) :: Interpolants 

    REAL(dp), DIMENSION(:,:,:), INTENT(out) :: Derivatives 
    REAL(dp) :: p000, p100, p010, p001, p011, p101, p110, p111, epsilon
    REAL(dp), DIMENSION(3) :: alpha, delta
    INTEGER :: i, j, il1, il2, il3

    epsilon = 1.d-200

    DO i = 1, SIZE(x2)

      ASSOCIATE( Coordinate1 => TS % States(1) % Values, &
                 Coordinate2 => TS % States(2) % Values, &
                 Coordinate3 => TS % States(3) % Values )

      CALL locate( Coordinate1, SIZE(Coordinate1), x1(i), il1 )
      CALL locate( Coordinate2, SIZE(Coordinate2), x2(i), il2 )
      CALL locate( Coordinate3, SIZE(Coordinate3), x3(i), il3 )

      IF ( LogInterp(1) ) THEN
      alpha(1) = ( 1.0d0 ) / ( x1(i) * LOG10( Coordinate1(il1+1) / Coordinate1(il1) ) )
      delta(1) = LOG10( x1(i) / Coordinate1(il1) ) / LOG10( Coordinate1(il1+1) / Coordinate1(il1) )
      ELSE
      alpha(1) = ( 1.0d0 ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      delta(1) = ( x1(i) - Coordinate1(il1) ) / ( Coordinate1(il1+1) - Coordinate1(il1) )
      END IF

      IF ( LogInterp(2) ) THEN
      alpha(2) = ( 1.0d0 ) / ( x2(i) * LOG10( Coordinate2(il2+1) / Coordinate2(il2) ) )
      delta(2) = LOG10( x2(i) / Coordinate2(il2) ) / LOG10( Coordinate2(il2+1) / Coordinate2(il2) )
      ELSE
      alpha(2) = ( 1.0d0 ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      delta(2) = ( x2(i) - Coordinate2(il2) ) / ( Coordinate2(il2+1) - Coordinate2(il2) )
      END IF

      IF ( LogInterp(3) ) THEN
      alpha(3) = ( 1.0d0 ) / ( x3(i) * LOG10( Coordinate3(il3+1) / Coordinate3(il3) ) )
      delta(3) = LOG10( x3(i) / Coordinate3(il3) ) / LOG10( Coordinate3(il3+1) / Coordinate3(il3) )
      ELSE
      alpha(3) = ( 1.0d0 ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
      delta(3) = ( x3(i) - Coordinate3(il3) ) / ( Coordinate3(il3+1) - Coordinate3(il3) )
      END IF

      END ASSOCIATE

      DO j = 1, DV % nVariables

        ASSOCIATE( Table => DV % Variables(j) % Values(:,:,:), &
                   Offset => DV % Offsets(j) )

        p000 = ( Table( il1  , il2  , il3   ) )
        p100 = ( Table( il1+1, il2  , il3   ) )
        p010 = ( Table( il1  , il2+1, il3   ) )
        p110 = ( Table( il1+1, il2+1, il3   ) )
        p001 = ( Table( il1  , il2  , il3+1 ) )
        p101 = ( Table( il1+1, il2  , il3+1 ) )
        p011 = ( Table( il1  , il2+1, il3+1 ) )
        p111 = ( Table( il1+1, il2+1, il3+1 ) )

        Interpolants(i,j) &
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

        Derivatives(i,1,j) &
          = ( (Interpolants(i,j) ) * alpha(1) &
              * ( (1.0_dp - delta(3)) * ( (delta(2) - 1.0_dp) * p000   &
                                      +  ( 1.0_dp - delta(2)) * p100   &
                                      -             delta(2)  * p010   &
                                      +             delta(2)  * p110 ) &
                           + delta(3) * ( (delta(2) - 1.0_dp) * p001   &
                                      +  ( 1.0_dp - delta(2)) * p101   &
                                      -             delta(2)  * p011   &
                                      +             delta(2)  * p111 ) ) )
  
        Derivatives(i,2,j) &
          = ( ( Interpolants(i,j) ) * alpha(2) &
                * ( (1.0_dp - delta(3) ) * ( (1.0_dp - delta(1)) * p000   &
                                         +             delta(1)  * p100   &
                                         -  ( 1.0_dp - delta(1)) * p010   &
                                         -             delta(1)  * p110 ) &
                              + delta(3) * ( (1.0_dp - delta(1)) * p001   &
                                         +             delta(1)  * p101   &
                                         -   (1.0_dp - delta(1)) * p011   &
                                         -             delta(1)  * p111 ) ) )

        Derivatives(i,3,j) &
          = ( ( Interpolants(i,j) ) * alpha(3) &
                                     * ( ( (delta(1) - 1.0_dp)) * (1.0_dp - delta(2)) * p000   &
                                         -            delta(1)  * (1.0_dp - delta(2)) * p100   &
                                         - ( 1.0_dp - delta(1)) *           delta(2)  * p010   &
                                         -            delta(1)  *           delta(2)  * p110   &
                                         +  (1.0_dp - delta(1)) * (1.0_dp - delta(2)) * p001   &
                                         +            delta(1)  * (1.0_dp - delta(2)) * p101   &
                                         +  (1.0_dp - delta(1)) *           delta(2)  * p011   &
                                         +            delta(1)  *           delta(2)  * p111 ) )
  
  
        END ASSOCIATE

      END DO
    END DO

  END SUBROUTINE LogInterpolateDifferentiateAllVariables 

  SUBROUTINE GetGamma1( x1, x2, x3, Coordinate1, Coordinate2, &
                Coordinate3, LogInterp, TS, DV, Gamma1 ) 

    REAL(dp), DIMENSION(:), INTENT(in) :: x1
    REAL(dp), DIMENSION(:), INTENT(in) :: x2
    REAL(dp), DIMENSION(:), INTENT(in) :: x3
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate1
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate2
    REAL(dp), DIMENSION(:), INTENT(in) :: Coordinate3
    LOGICAL, DIMENSION(3), INTENT(in)  :: LogInterp
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
               ( E_Internal, Energy_Table, Temp_Table, Offset, Temperature ) 

    REAL(dp), DIMENSION(:), INTENT(in) :: E_Internal
    REAL(dp), DIMENSION(:), INTENT(in) :: Temp_Table 
    REAL(dp), DIMENSION(:), INTENT(in) :: Energy_Table
    REAL(dp), INTENT(in) :: Offset
    REAL(dp), DIMENSION(:), INTENT(out) :: Temperature
    REAL(dp) :: delta, epsilon
    INTEGER :: i, j
 
    epsilon = 1.d-200
    DO j = 1, SIZE( E_Internal )

      CALL locate( Energy_Table, SIZE( Energy_Table ), E_Internal(j), i )
      IF ( i == SIZE(Energy_Table) ) THEN
        STOP
      END IF

      IF ( i == 0 ) THEN
        Temperature(j) = 0.d0
        CYCLE
      END IF
      Temperature(j) =  Temp_Table(i) + ( Temp_Table(i+1) - Temp_Table(i) ) &
                        * ( ( E_internal(j) - Energy_Table(i) )             &
                        / ( Energy_Table(i+1) - Energy_Table(i) ) )
    END DO

  END SUBROUTINE ComputeTempFromIntEnergy

END MODULE wlInterpolationModule

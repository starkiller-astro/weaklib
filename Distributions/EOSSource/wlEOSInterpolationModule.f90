MODULE wlEOSInterpolationModule

  USE wlKindModule,                 ONLY: &
    dp
  USE wlInterpolationModule,        ONLY: &
    LogInterpolateDifferentiateSingleVariable_3D_Custom
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D, &
    TriLinear
  USE wlThermoStateModule,          ONLY: &
    ThermoStateType
  USE wlDependentVariablesModule,   ONLY: &
    DependentVariablesType
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: LogInterpolateAllVariables
  PUBLIC :: LogInterpolateDifferentiateAllVariables
  PUBLIC :: GetGamma1
  PUBLIC :: EOSTableQuery
  PUBLIC :: MonotonicityCheck

  INTERFACE LogInterpolateAllVariables
    MODULE PROCEDURE LogInterpolateAllVariables_3D
    MODULE PROCEDURE LogInterpolateAllVariables_3D_Custom
  END INTERFACE LogInterpolateAllVariables

  REAL(dp), PARAMETER :: ln10 = LOG(10.d0)

CONTAINS


  SUBROUTINE LogInterpolateAllVariables_3D &
    ( x1, x2, x3, LogInterp, TS, DV, Interpolants, MaskVar )

    REAL(dp), DIMENSION(:),           INTENT(in)  :: x1
    REAL(dp), DIMENSION(:),           INTENT(in)  :: x2
    REAL(dp), DIMENSION(:),           INTENT(in)  :: x3
    INTEGER,  DIMENSION(3),           INTENT(in)  :: LogInterp 
    TYPE(ThermoStateType),            INTENT(in)  :: TS
    TYPE(DependentVariablesType),     INTENT(in)  :: DV
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(in)  :: MaskVar
    REAL(dp), DIMENSION(:,:),         INTENT(out) :: Interpolants 

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

      il1(i)  = Index1D( x1(i), Coordinate1, dim1 )
      il2(i)  = Index1D( x2(i), Coordinate2, dim2 )
      il3(i)  = Index1D( x3(i), Coordinate3, dim3 )

      IF ( LogInterp(1) == 1 ) THEN
        delta(1,i) = LOG10( x1(i) / Coordinate1(il1(i)) ) &
                       / LOG10( Coordinate1(il1(i)+1) / Coordinate1(il1(i)) )
      ELSE
        delta(1,i) = ( x1(i) - Coordinate1(il1(i)) ) &
                       / ( Coordinate1(il1(i)+1) - Coordinate1(il1(i)) )
      END IF

      IF ( LogInterp(2) == 1 ) THEN
        delta(2,i) = LOG10( x2(i) / Coordinate2(il2(i)) ) &
                       / LOG10( Coordinate2(il2(i)+1) / Coordinate2(il2(i)) )
      ELSE
        delta(2,i) = ( x2(i) - Coordinate2(il2(i)) ) &
                       / ( Coordinate2(il2(i)+1) - Coordinate2(il2(i)) )
      END IF

      IF ( LogInterp(3) == 1 ) THEN
        delta(3,i) = LOG10( x3(i) / Coordinate3(il3(i)) ) &
                       / LOG10( Coordinate3(il3(i)+1) / Coordinate3(il3(i)) )
      ELSE
        delta(3,i) = ( x3(i) - Coordinate3(il3(i)) ) &
                       / ( Coordinate3(il3(i)+1) - Coordinate3(il3(i)) )
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


  SUBROUTINE LogInterpolateDifferentiateAllVariables &
    ( x1, x2, x3, LogInterp, TS, DV, Interpolants, Derivatives, MaskVar )

    REAL(dp), DIMENSION(:),           INTENT(in)  :: x1
    REAL(dp), DIMENSION(:),           INTENT(in)  :: x2
    REAL(dp), DIMENSION(:),           INTENT(in)  :: x3
    INTEGER,  DIMENSION(3),           INTENT(in)  :: LogInterp 
    TYPE(ThermoStateType),            INTENT(in)  :: TS
    TYPE(DependentVariablesType),     INTENT(in)  :: DV
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(in)  :: MaskVar
    REAL(dp), DIMENSION(:,:),         INTENT(out) :: Interpolants 

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

    ALLOCATE &
      ( work_mask( Masksize ), alpha( 3, Masksize ), delta( 3, Masksize ), &
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

      il1(i) = Index1D( x1(i), Coordinate1, dim1 )
      il2(i) = Index1D( x2(i), Coordinate2, dim2 )
      il3(i) = Index1D( x3(i), Coordinate3, dim3 )

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


  SUBROUTINE GetGamma1 &
    ( x1, x2, x3, Coordinate1, Coordinate2, Coordinate3, LogInterp, &
      TS, DV, Gamma1 ) 

    REAL(dp), DIMENSION(:),              INTENT(in)  :: x1
    REAL(dp), DIMENSION(:),              INTENT(in)  :: x2
    REAL(dp), DIMENSION(:),              INTENT(in)  :: x3
    REAL(dp), DIMENSION(:),              INTENT(in)  :: Coordinate1
    REAL(dp), DIMENSION(:),              INTENT(in)  :: Coordinate2
    REAL(dp), DIMENSION(:),              INTENT(in)  :: Coordinate3
    INTEGER,  DIMENSION(3),              INTENT(in)  :: LogInterp
    TYPE(ThermoStateType),               INTENT(in)  :: TS
    TYPE(DependentVariablesType),        INTENT(in)  :: DV
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: Gamma1

    INTEGER :: i
    REAL(dp), DIMENSION(:),   ALLOCATABLE :: Interpolant 
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative 

    ALLOCATE( Interpolant( SIZE(x1)    ) )
    ALLOCATE( Gamma1     ( SIZE(x1)    ) )
    ALLOCATE( Derivative ( SIZE(x1), 3 ) )

    CALL LogInterpolateDifferentiateSingleVariable_3D_Custom &
           ( x1, x2, x3, &
             TS % States(1) % Values(:), &
             TS % States(2) % Values(:), &
             TS % States(3) % Values(:), &
             DV % Offsets(1), &
             DV % Variables(1) % Values(:,:,:), &
             Interpolant(:), Derivative(:,:) )

    DO i = 1, SIZE(x1) 
      Gamma1(i) =  ( x1(i) / Interpolant(i) ) * Derivative(i,1) 
    END DO

  END SUBROUTINE GetGamma1 


  SUBROUTINE EOSTableQuery &
    ( rho, T, Ye, LogInterp, TS, DV, Interpolants, Verbose_Option )

    INTEGER                                   :: i, j
    REAL(dp), DIMENSION(:),       INTENT(in)  :: rho
    REAL(dp), DIMENSION(:),       INTENT(in)  :: T
    REAL(dp), DIMENSION(:),       INTENT(in)  :: Ye
    INTEGER,  DIMENSION(3),       INTENT(in)  :: LogInterp
    TYPE(ThermoStateType),        INTENT(in)  :: TS
    TYPE(DependentVariablesType), INTENT(in)  :: DV
    REAL(dp), DIMENSION(:,:),     INTENT(out) :: Interpolants
    LOGICAL,                INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL        :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    CALL LogInterpolateAllVariables &
           ( rho, T, Ye, LogInterp, TS, DV, Interpolants )

    IF( Verbose )THEN
      DO i = 1, SIZE(rho)
        WRITE(*,*) 'Rho=', rho(i), 'T=', T(i), 'Ye=', Ye(i)
        DO j = 1, DV % nVariables
          WRITE(*,*) DV % Names(j), Interpolants(i,j)
        END DO
      END DO
    END IF

  END SUBROUTINE EOSTableQuery


  SUBROUTINE MonotonicityCheck( Table, Nrho, NT, NYe, Axis, Repaired )

    REAL(dp), DIMENSION(:,:,:), INTENT(in) :: Table
    INTEGER,                    INTENT(in) :: Nrho 
    INTEGER,                    INTENT(in) :: NT
    INTEGER,                    INTENT(in) :: NYe
    INTEGER,                    INTENT(in) :: Axis
    INTEGER,  DIMENSION(:,:,:), INTENT(in) :: Repaired

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


END MODULE wlEOSInterpolationModule

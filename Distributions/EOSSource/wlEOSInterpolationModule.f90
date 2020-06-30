MODULE wlEOSInterpolationModule

  USE wlKindModule,                 ONLY: &
    dp
  USE wlInterpolationModule,        ONLY: &
    locate, &
    Index1D, &
    TriLinear, &
    LogInterpolateSingleVariable, &
    LogInterpolateDifferentiateSingleVariable
  USE wlThermoStateModule,          ONLY: &
    ThermoStateType
  USE wlDependentVariablesModule,   ONLY: &
    DependentVariablesType
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeTempFromIntEnergy
  PUBLIC :: ComputeTempFromIntEnergy_Lookup
  PUBLIC :: ComputeTempFromIntEnergy_Bisection
  PUBLIC :: ComputeTempFromIntEnergy_Secant
  PUBLIC :: ComputeTempFromEntropy
  PUBLIC :: ComputeTempFromPressure
  PUBLIC :: ComputeTempFromPressure_Bisection
  PUBLIC :: ComputeTempForVector
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
    INTEGER  :: i, j, Error
    INTEGER  :: ilD, ilY, ilE, ilE1, ilE2, ilE3, ilE4, ilEa, ilEb
    INTEGER  :: nPtsE
    LOGICAL :: LocalRoot
    REAL(dp) :: logE, tmpE, f_a, f_b, T_a, T_b, E_a, E_b
    REAL(dp) :: Ds(1:2), Ys(1:2)
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Ts, EvsT
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Es

    DO i = 1, SIZE( D )

      Error = 0

      ilD = Index1D( D(i), D_table, SIZE( D_table ) )
      ilY = Index1D( Y(i), Y_table, SIZE( Y_table ) )

      Ds(1:2) = D_Table(ilD:ilD+1)
      Ys(1:2) = Y_Table(ilY:ilY+1)

      logE = LOG10( E(i) + Offset )
      ilE1 = Index1D( logE, E_table(ilD,  :,ilY  ), SIZE( T_Table ) )
      ilE2 = Index1D( logE, E_table(ilD+1,:,ilY  ), SIZE( T_Table ) )
      ilE3 = Index1D( logE, E_table(ilD,  :,ilY+1), SIZE( T_Table ) )
      ilE4 = Index1D( logE, E_table(ilD+1,:,ilY+1), SIZE( T_Table ) )

      ilEa = MAX( MINVAL( [ilE1,ilE2,ilE3,ilE4] ) - 1, 1 )
      ilEb = MIN( MAXVAL( [ilE1,ilE2,ilE3,ilE4] ) + 2, SIZE( T_Table ) )

      nPtsE = ilEb - ilEa + 1
      ALLOCATE( Ts(1:nPtsE), Es(1:2,1:nPtsE,1:2), EvsT(1:nPtsE) )
      Ts(1:nPtsE) = T_Table(ilEa:ilEb)
      Es(1:2,1:nPtsE,1:2) = E_Table(ilD:ilD+1,ilEa:ilEb,ilY:ilY+1)

      DO j = 1, nPtsE
        CALL LogInterpolateSingleVariable &
               ( D(i), Ts(j), Y(i), Ds, Ts, Ys, Offset, Es, EvsT(j) )
      END DO

      ! --- Pick Highest Temperature Root ---
      DO j = nPtsE - 1, 1, -1

        f_a = E(i) - EvsT(j)
        f_b = E(i) - EvsT(j+1)

        LocalRoot = ( f_a * f_b < 0.0_DP )

        IF ( LocalRoot ) THEN
          ilE = j
          EXIT
        END IF

      END DO

      tmpE = E(i)

      T_a = Ts(1)
      T_b = Ts(nPtsE)

      E_a = EvsT(1)
      E_b = EvsT(nPtsE)

      f_a = E(i) - E_a
      f_b = E(i) - E_b

      IF( .not. LocalRoot )THEN

        IF( Debug )THEN
          WRITE(*,*)
          WRITE(*,'(A4,A)') &
            '', 'Warning: ComputeTempFromIntEnergy_Lookup'
          WRITE(*,'(A6,A20,ES10.4E2,A9,ES10.4E2)') &
            '', 'No Root Between T = ', Ts(1), ' and T = ', Ts(nPtsE)
          WRITE(*,*)
          WRITE(*,*) '  nPtsE = ', nPtsE
          WRITE(*,*) '  EvsT  = ', EvsT
          WRITE(*,*) '    ia  = ', ilEa
          WRITE(*,*) '    ib  = ', ilEb
          WRITE(*,*) '    Ta  = ', T_a
          WRITE(*,*) '    Tb  = ', T_b
          WRITE(*,*) '    Ea  = ', E_a
          WRITE(*,*) '    Eb  = ', E_b
          WRITE(*,*) '     i  = ', i
          WRITE(*,*) '     E  = ', E(i)
          WRITE(*,*) '     D  = ', D(i)
          WRITE(*,*) '     Y  = ', Y(i)
          WRITE(*,*)
          STOP
        END IF

        ! --- Reset Energy Density ---

        IF( E(i) < EvsT(1) .AND. E(i) < EvsT(nPtsE) )THEN

          ilE = 1
          tmpE = 0.5 * ( EvsT(1) + EvsT(2) )

        ELSE

          ilE = nPtsE - 1
          tmpE = 0.5 * ( EvsT(nPtsE-1) + EvsT(nPtsE) )

        END IF

        Error = 1

      END IF

      !ilE = Index1D( tmpE, EvsT, nPtsE )

      T(i) &
        = 10.d0**( LOG10( Ts(ilE) ) &
                   + LOG10( Ts(ilE+1)/Ts(ilE) ) &
                     * LOG10( (tmpE+Offset)/(EvsT(ilE)+Offset) ) &
                       / LOG10( (EvsT(ilE+1)+Offset)/(EvsT(ilE)+Offset) ) )

      IF( Error == 1 .AND. Debug )THEN
        WRITE(*,*)
        WRITE(*,*) '  T = ', T(i)
        WRITE(*,*)
      END IF

      DEALLOCATE( Ts, Es, EvsT )

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

      CALL locate( Coordinate1, dim1, x1(i), il1(i) )
      CALL locate( Coordinate2, dim2, x2(i), il2(i) )
      CALL locate( Coordinate3, dim3, x3(i), il3(i) )

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


  SUBROUTINE ComputeTempForVector &
    ( rho, alt, ye, EOSTable, InputFlag, InputMask, Temperature )

    REAL(dp), DIMENSION(:), INTENT(in)         :: rho
    REAL(dp), DIMENSION(:), INTENT(in)         :: alt
    REAL(dp), DIMENSION(:), INTENT(in)         :: ye
    TYPE(EquationOfStateTableType), INTENT(in) :: EOSTable
    INTEGER, DIMENSION(:), INTENT(in)          :: InputFlag
    LOGICAL, DIMENSION(:), INTENT(in)          :: InputMask
    REAL(dp), DIMENSION(:), INTENT(out)        :: Temperature

    INTEGER                                    :: i, ni, nf
    REAL(dp), DIMENSION(1)                     :: tbuff

    ni = LBOUND( rho, DIM = 1 )
    nf = UBOUND( rho, DIM = 1 )

    ASSOCIATE &
      ( WL_EINT    => EOSTable % DV % Indices % iInternalEnergyDensity, &
        WL_ENTROPY => EOSTable % DV % Indices % iEntropyPerBaryon,      &
        WL_PRESS   => EOSTable % DV % Indices % iPressure )

    DO i = ni, nf

      IF ( .not. InputMask(i) ) CYCLE

      IF ( InputFlag(i) == WL_EINT ) THEN

        CALL ComputeTempFromIntEnergy( rho(i), alt(i), ye(i),      &
           EOSTable % TS % States(1) % Values,                     &
           EOSTable % TS % States(2) % Values,                     &
           EOSTable % TS % States(3) % Values,                     &
           EOSTable % TS % LogInterp,                              &
           EOSTable % DV % Variables(WL_EINT) % Values(:,:,:),     &
           EOSTable % DV % Offsets(WL_EINT), tbuff )

      ELSEIF ( InputFlag(i) == WL_ENTROPY ) THEN

        CALL ComputeTempFromEntropy( rho(i), alt(i), ye(i),        &
           EOSTable % TS % States(1) % Values,                     &
           EOSTable % TS % States(2) % Values,                     &
           EOSTable % TS % States(3) % Values,                     &
           EOSTable % TS % LogInterp,                              &
           EOSTable % DV % Variables(WL_ENTROPY) % Values(:,:,:),  &
           EOSTable % DV % Offsets(WL_ENTROPY), tbuff )

      ELSEIF ( InputFlag(i) == WL_PRESS ) THEN

        CALL ComputeTempFromIntEnergy( rho(i), alt(i), ye(i),      &
           EOSTable % TS % States(1) % Values,                     &
           EOSTable % TS % States(2) % Values,                     &
           EOSTable % TS % States(3) % Values,                     &
           EOSTable % TS % LogInterp,                              &
           EOSTable % DV % Variables(WL_PRESS) % Values(:,:,:),    &
           EOSTable % DV % Offsets(WL_PRESS), tbuff )

      ELSE

        WRITE(*,*) &
          "Invalid thermodynamic variable flag in ComputeTempForVector: ", &
          InputMask(i)
        STOP

      END IF
      Temperature(i) = tbuff(1) 

    END DO

    END ASSOCIATE

  END SUBROUTINE ComputeTempForVector


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

    CALL LogInterpolateDifferentiateSingleVariable &
           ( x1, x2, x3, &
             TS % States(1) % Values(:), &
             TS % States(2) % Values(:), &
             TS % States(3) % Values(:), &
             LogInterp, DV % Offsets(1), &
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

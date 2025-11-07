MODULE wlEOSInterpolationModule_Old

  USE wlKindModule,                 ONLY: &
    dp
  USE wlInterpolationModule,        ONLY: &
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

      i = Index1D( eibuff(j), energy_array, SIZE( energy_array ) )
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

      i  = Index1D( sbuff(j), entropy_array, SIZE( entropy_array ) )
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

      i  = Index1D( pbuff(j), pressure_array, SIZE( pressure_array ) )
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


END MODULE wlEOSInterpolationModule_Old

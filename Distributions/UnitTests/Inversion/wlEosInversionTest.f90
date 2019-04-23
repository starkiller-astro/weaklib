PROGRAM wlEosInversionTest

  USE wlKindModule, ONLY: &
    dp
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    ComputeTempFromIntEnergy, &
    ComputeTempFromEntropy
  USE wlEOSInversionModule, ONLY: &
    InitializeEOSInversion, &
    ComputeTemperatureWith_DEY, &
    ComputeTemperatureWith_DSY, &
    DescribeEOSInversionError

  IMPLICIT NONE

  INTEGER :: &
    n_rndm, &
    iP
  INTEGER, PARAMETER :: &
    nPoints = 2**16, &
    iD = 1, iT = 2, iY = 3
  INTEGER, DIMENSION(nPoints) :: &
    Error_E, &
    Error_P, &
    Error_S
  REAL(dp) :: &
    tBegin, &
    tEnd, &
    Amp
  REAL(dp), DIMENSION(nPoints) :: &
    D, T, Y, E, P, S, T_E, T_P, T_S, dT, &
    rndm_D, &
    rndm_T, &
    rndm_Y, &
    Error
  TYPE(EquationOfStateTableType) :: &
    EosTab

  CALL InitializeHDF( )
  CALL ReadEquationOfStateTableHDF( EosTab, "EquationOfStateTable.h5" )
  CALL FinalizeHDF( )

  associate &
    ( iDtab => EosTab % TS % Indices % iRho, &
      iTtab => EosTab % TS % Indices % iT, &
      iYtab => EosTab % TS % Indices % iYe, &
      iEtab => EosTab % DV % Indices % iInternalEnergyDensity, &
      iPtab => EosTab % DV % Indices % iPressure, &
      iStab => EosTab % DV % Indices % iEntropyPerBaryon )

  associate &
    ( Dtab => EosTab % TS % States(iDtab) % Values, &
      Ttab => EosTab % TS % States(iTtab) % Values, &
      Ytab => EosTab % TS % States(iYtab) % Values, &
      Etab => EosTab % DV % Variables(iEtab) % Values, &
      Ptab => EosTab % DV % Variables(iPtab) % Values, &
      Stab => EosTab % DV % Variables(iStab) % Values, &
      OS_E => EosTab % DV % Offsets(iEtab), &
      OS_P => EosTab % DV % Offsets(iPtab), &
      OS_S => EosTab % DV % Offsets(iStab) )

  CALL InitializeEOSInversion &
         ( Dtab, Ttab, Ytab, &
           10.0d0**( Etab ) - OS_E, &
           10.0d0**( Ptab ) - OS_P, &
           10.0d0**( Stab ) - OS_S, &
           Verbose_Option = .TRUE. )

  end associate

  end associate

  WRITE(*,*)
  WRITE(*,'(A4,A10,I10.10)') '', 'nPoints = ', nPoints

  ! --- Random Sampling of EOS Table ---

  n_rndm = nPoints

  ! --- Initialize Density Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )

  associate &
    ( minD => EosTab % TS % MinValues(iD), &
      maxD => EosTab % TS % MaxValues(iD) )

  D(:) = 10**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  PRINT*, "Min/Max D = ", MINVAL( D ), MAXVAL( D )
  PRINT*, "            ", minD, maxD

  end associate

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  associate &
    ( minT => EosTab % TS % MinValues(iT), &
      maxT => EosTab % TS % MaxValues(iT) )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  PRINT*, "Min/Max T = ", MINVAL( T ), MAXVAL( T )
  PRINT*, "            ", minT, maxT

  end associate

  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Y )

  associate &
    ( minY => EosTab % TS % MinValues(iY), &
      maxY => EosTab % TS % MaxValues(iY) )

  Y(:) = minY + ( maxY - minY ) * rndm_Y

  PRINT*, "Min/Max Y = ", MINVAL( Y ), MAXVAL( Y )
  PRINT*, "            ", minY, maxY

  end associate

  associate &
    ( iDtab => EosTab % TS % Indices % iRho, &
      iTtab => EosTab % TS % Indices % iT, &
      iYtab => EosTab % TS % Indices % iYe, &
      iEtab => EosTab % DV % Indices % iInternalEnergyDensity, &
      iPtab => EosTab % DV % Indices % iPressure, &
      iStab => EosTab % DV % Indices % iEntropyPerBaryon )

  associate &
    ( Dtab => EosTab % TS % States(iDtab) % Values, &
      Ttab => EosTab % TS % States(iTtab) % Values, &
      Ytab => EosTab % TS % States(iYtab) % Values, &
      Etab => EosTab % DV % Variables(iEtab) % Values, &
      Ptab => EosTab % DV % Variables(iPtab) % Values, &
      Stab => EosTab % DV % Variables(iStab) % Values, &
      OS_E => EosTab % DV % Offsets(iEtab), &
      OS_P => EosTab % DV % Offsets(iPtab), &
      OS_S => EosTab % DV % Offsets(iStab) )

  ! --- Compute Internal Energy, Pressure, and Entropy ---

  CALL LogInterpolateSingleVariable &
         ( D, T, Y, Dtab, Ttab, Ytab, OS_E, Etab, E )

  WRITE(*,*)
  WRITE(*,*) "Min/Max E = ", MINVAL( E ), MAXVAL( E )

  CALL LogInterpolateSingleVariable &
         ( D, T, Y, Dtab, Ttab, Ytab, OS_P, Ptab, P )

  WRITE(*,*) "Min/Max P = ", MINVAL( P ), MAXVAL( P )

  CALL LogInterpolateSingleVariable &
         ( D, T, Y, Dtab, Ttab, Ytab, OS_S, Stab, S )

  WRITE(*,*) "Min/Max S = ", MINVAL( S ), MAXVAL( S )

  ! -------------------------------------------------------------------
  ! --- Recover Temperature from Internal Energy ----------------------
  ! -------------------------------------------------------------------

  Amp = 0.001_dp ! --- Perturbation Amplitude

  PRINT*
  PRINT*, "Perturbation Amplitude: ", Amp

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  associate &
    ( minT => 1.0001 * EosTab % TS % MinValues(iT), &
      maxT => 0.9999 * EosTab % TS % MaxValues(iT) )

  T_E = MIN( MAX( T + dT, minT ), maxT )

  end associate

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DEY &
         ( D, E, Y, Dtab, Ttab, Ytab, Etab, OS_E, T_E, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_E )

  CALL CPU_TIME( tEnd )

  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 )THEN
      CALL DescribeEOSInversionError( Error_E(iP) )
    END IF
  END DO

  Error = ABS( T - T_E ) / T

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEY:"
  PRINT*, "CPU_TIME = ", tEnd - tBegin
  PRINT*, MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "MINVAL(T) = ", MINVAL( T_E )

  ! -------------------------------------------------------------------

  CALL CPU_TIME( tBegin )

  DO iP = 1, nPoints

    CALL ComputeTempFromIntEnergy &
           ( D(iP), E(iP), Y(iP), Dtab, Ttab, Ytab, [1,1,0], &
             Etab, OS_E, T_E(iP:iP) )

  END DO

  CALL CPU_TIME( tEnd )

  Error = ABS( T - T_E ) / T

  PRINT*
  PRINT*, "ComputeTempFromIntEnergy:"
  PRINT*, "CPU_TIME = ", tEnd - tBegin
  PRINT*, MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "MINVAL(T) = ", MINVAL( T_E )

  STOP

  ! -------------------------------------------------------------------
  ! --- Recover Temperature from Entropy Per Baryon -------------------
  ! -------------------------------------------------------------------

  associate &
    ( minT => 1.0001 * EosTab % TS % MinValues(iT), &
      maxT => 0.9999 * EosTab % TS % MaxValues(iT) )

  T_S = MIN( MAX( T + dT, minT ), maxT )

  end associate

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DSY &
         ( D, S, Y, Dtab, Ttab, Ytab, Stab, OS_S, T_S, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_S )

  CALL CPU_TIME( tEnd )

  Error = ABS( T - T_S ) / T

  PRINT*
  PRINT*, "ComputeTemperatureWith_DSY:"
  PRINT*, "CPU_TIME = ", tEnd - tBegin
  PRINT*, MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "MINVAL(T) = ", MINVAL( T_S )

  ! -------------------------------------------------------------------

  CALL CPU_TIME( tBegin )

  DO iP = 1, nPoints

    CALL ComputeTempFromEntropy &
           ( D(iP), S(iP), Y(iP), Dtab, Ttab, Ytab, [1,1,0], &
             Stab, OS_S, T_S(iP:iP) )

  END DO

  CALL CPU_TIME( tEnd )

  Error = ABS( T - T_S ) / T

  PRINT*
  PRINT*, "ComputeTempFromEntropy:"
  PRINT*, "CPU_TIME = ", tEnd - tBegin
  PRINT*, MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "MINVAL(T) = ", MINVAL( T_S )

  end associate ! --- Dtab, etc.

  end associate ! --- iDtab, etc.

END PROGRAM wlEosInversionTest

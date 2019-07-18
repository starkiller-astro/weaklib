PROGRAM wlEosPerformanceTest

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
    LogInterpolateDifferentiateSingleVariable
  USE wlEOSInterpolationModule, ONLY: &
    LogInterpolateAllVariables

  IMPLICIT NONE

  INTEGER :: &
    n_rndm, &
    iPoint, iV, iDim
  INTEGER, DIMENSION(3) :: &
    LogInterp = [ 1, 1, 0 ]
  INTEGER, PARAMETER :: &
    iD = 1, iT = 2, iY = 3, &
    nPoints = 2**18, &
    nVariables = 15
  REAL(dp) :: &
    tBegin, &
    tEnd
  REAL(dp), DIMENSION(nPoints) :: &
    D, T, Y, &
    rndm_D, &
    rndm_T, &
    rndm_Y, &
    Interpolant1, &
    Interpolant2, &
    Interpolant3, &
    Error
  REAL(dp), DIMENSION(nPoints,3) :: &
    Derivative1, &
    Derivative2, &
    Derivative3, &
    ErrorD
  REAL(dp), DIMENSION(nPoints,nVariables) :: &
    Interpolants1
  REAL(dp), DIMENSION(nVariables,nPoints) :: &
    Interpolants2
  TYPE(EquationOfStateTableType) :: &
    EosTab

  CALL InitializeHDF( )
  CALL ReadEquationOfStateTableHDF( EosTab, "EquationOfStateTable.h5" )
  CALL FinalizeHDF( )

  WRITE(*,*)
  WRITE(*,'(A4,A10,I10.10)') '', 'nPoints = ', nPoints
  WRITE(*,*)

  n_rndm = nPoints

  ! --- Initialize Density Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )

  associate &
    ( minD => EosTab % TS % MinValues(iD), &
      maxD => EosTab % TS % MaxValues(iD) )

  D(:) = 10**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  end associate

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  associate &
    ( minT => EosTab % TS % MinValues(iT), &
      maxT => EosTab % TS % MaxValues(iT) )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  end associate

  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Y )

  associate &
    ( minY => EosTab % TS % MinValues(iY), &
      maxY => EosTab % TS % MaxValues(iY) )

  Y(:) = minY + ( maxY - minY ) * rndm_Y

  end associate

  ! --- Single Variable ---

  ASSOCIATE &
    ( iDtab => EosTab % TS % Indices % iRho, &
      iTtab => EosTab % TS % Indices % iT,   &
      iYtab => EosTab % TS % Indices % iYe,  &
      iPtab => EosTab % DV % Indices % iPressure )

  ASSOCIATE &
    ( Dtab => EosTab % TS % States   (iDtab) % Values, &
      Ttab => EosTab % TS % States   (iTtab) % Values, &
      Ytab => EosTab % TS % States   (iYtab) % Values, &
      Ptab => EosTab % DV % Variables(iPtab) % Values, &
      OS   => EosTab % DV % Offsets(iPtab) )

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateSingleVariable &
           ( D, T, Y, Dtab, Ttab, Ytab, LogInterp, OS, Ptab, Interpolant1 )

  CALL CPU_TIME( tEnd )

  WRITE(*,*)
  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_3D: ', tEnd - tBegin
  WRITE(*,*)

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateSingleVariable &
           ( D, T, Y, Dtab, Ttab, Ytab, OS, Ptab, Interpolant2 )

  CALL CPU_TIME( tEnd )

  WRITE(*,*)
  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_3D_Custom: ', tEnd - tBegin
  WRITE(*,*)

  CALL CPU_TIME( tBegin )

  DO iPoint = 1, nPoints

    CALL LogInterpolateSingleVariable &
           ( D(iPoint), T(iPoint), Y(iPoint), Dtab, Ttab, Ytab, OS, Ptab, &
             Interpolant3(iPoint) )

  END DO

  CALL CPU_TIME( tEnd )

  WRITE(*,*)
  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_3D_Custom_Point: ', tEnd - tBegin
  WRITE(*,*)

  Error = ABS( ( Interpolant1 - Interpolant2 ) / Interpolant1 )

  WRITE(*,*) MAXVAL( Error ), MINVAL( Error )

  Error = ABS( ( Interpolant1 - Interpolant3 ) / Interpolant1 )

  WRITE(*,*) MAXVAL( Error ), MINVAL( Error )

  ! --- Single Variable (Interpolate and Differentiate) ---

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateDifferentiateSingleVariable &
         ( D, T, Y, Dtab, Ttab, Ytab, LogInterp, OS, Ptab, Interpolant1, &
           Derivative1 )

  CALL CPU_TIME( tEnd )

  WRITE(*,*)
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateDifferentiateSingleVariable_3D: ', &
    tEnd - tBegin
  WRITE(*,*)

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateDifferentiateSingleVariable &
         ( D, T, Y, Dtab, Ttab, Ytab, OS, Ptab, Interpolant2, &
           Derivative2 )

  CALL CPU_TIME( tEnd )

  WRITE(*,*)
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateDifferentiateSingleVariable_3D_Custom: ', &
    tEnd - tBegin
  WRITE(*,*)

  CALL CPU_TIME( tBegin )

  DO iPoint = 1, nPoints

    CALL LogInterpolateDifferentiateSingleVariable &
           ( D(iPoint), T(iPoint), Y(iPoint), Dtab, Ttab, Ytab, OS, Ptab, &
             Interpolant3(iPoint), Derivative3(iPoint,:) )

  END DO

  CALL CPU_TIME( tEnd )

  WRITE(*,*)
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateDifferentiateSingleVariable_3D_Custom_Point: ', &
    tEnd - tBegin
  WRITE(*,*)

  Error = ABS( Interpolant1 - Interpolant2 ) / ABS( Interpolant1 )
  WRITE(*,*) MAXVAL( Error ), MINVAL( Error )
  DO iDim = 1, 3
    ErrorD(:,iDim) &
      = ABS( ( Derivative1(:,iDim) - Derivative2(:,iDim) ) ) &
          / MAX( ABS( Derivative1(:,iDim) ), SQRT( TINY( 1.0_dp ) ) )
    WRITE(*,*) '  ', MAXVAL( ErrorD(:,iDim) ), MINVAL( ErrorD(:,iDim) )
  END DO

  Error = ABS( Interpolant1 - Interpolant3 ) / ABS( Interpolant1 )
  WRITE(*,*) MAXVAL( Error ), MINVAL( Error )
  DO iDim = 1, 3
    ErrorD(:,iDim) &
      = ABS( ( Derivative1(:,iDim) - Derivative3(:,iDim) ) ) &
          / MAX( ABS( Derivative1(:,iDim) ), SQRT( TINY( 1.0_dp ) ) )
    WRITE(*,*) '  ', MAXVAL( ErrorD(:,iDim) ), MINVAL( ErrorD(:,iDim) )    
  END DO

  END ASSOCIATE ! Dtab, etc.

  END ASSOCIATE ! iDtab, etc.

  ! --- All Variables ---

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateAllVariables &
         ( D, T, Y, LogInterp, EosTab % TS, EosTab % DV, Interpolants1 )

  CALL CPU_TIME( tEnd )

  WRITE(*,*)
  WRITE(*,'(A4,A40,ES10.4E2)') &
    '', 'LogInterpolateAllVariables_3D: ', tEnd - tBegin
  WRITE(*,*)

  ASSOCIATE &
    ( iDtab => EosTab % TS % Indices % iRho, &
      iTtab => EosTab % TS % Indices % iT,   &
      iYtab => EosTab % TS % Indices % iYe )

  ASSOCIATE &
    ( Dtab => EosTab % TS % States(iDtab) % Values, &
      Ttab => EosTab % TS % States(iTtab) % Values, &
      Ytab => EosTab % TS % States(iYtab) % Values )

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateAllVariables &
         ( D, T, Y, Dtab, Ttab, Ytab, EosTab % DV, Interpolants2 )

  CALL CPU_TIME( tEnd )

  END ASSOCIATE ! Dtab, etc.

  END ASSOCIATE ! iDtab, etc.

  WRITE(*,*)
  WRITE(*,'(A4,A40,ES10.4E2)') &
    '', 'LogInterpolateAllVariables_3D_Custom: ', tEnd - tBegin
  WRITE(*,*)

  DO iV = 1, nVariables
    Error = ABS( Interpolants1(:,iV) - Interpolants2(iV,:) ) &
              / ABS( Interpolants1(:,iV) )
    WRITE(*,*) iV, EosTab % DV % Names(iV), MAXVAL( Error ), MINVAL( Error )
  END DO
  WRITE(*,*)

END PROGRAM wlEosPerformanceTest

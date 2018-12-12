PROGRAM wlOpacityPerformanceTest

  USE wlKindModule, ONLY: &
    dp
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable

  IMPLICIT NONE

  INTEGER :: &
    n_rndm, &
    iPoint
  INTEGER, DIMENSION(4) :: &
    LogInterp = [ 1, 1, 1, 0 ]
  INTEGER, PARAMETER :: &
    iD = 1, iT = 2, iY = 3, &
    nPoints = 2**20
  REAL(dp) :: &
    tBegin, &
    tEND
  REAL(dp), DIMENSION(nPoints) :: &
    E, D, T, Y, &
    rndm_E, &
    rndm_D, &
    rndm_T, &
    rndm_Y, &
    Interpolant1, &
    Interpolant2
  TYPE(OpacityTableType) :: &
    OpTab

  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpTab, "OpacityTable.h5" )
  CALL FinalizeHDF( )

  WRITE(*,*)
  WRITE(*,'(A4,A12,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min E = ', OpTab % EnergyGrid % MinValue, &
    ' / ', OpTab % EnergyGrid % MaxValue
  WRITE(*,'(A4,A12,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min D = ', OpTab % EOSTable % TS % MinValues(iD), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iD)
  WRITE(*,'(A4,A12,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min T = ', OpTab % EOSTable % TS % MinValues(iT), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iT)
  WRITE(*,'(A4,A12,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min Y = ', OpTab % EOSTable % TS % MinValues(iY), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iY)
  WRITE(*,*)
  WRITE(*,'(A4,A14,I10.10,A7)') &
    '', 'Interpolating ', nPoints, ' points'
  WRITE(*,*)

  n_rndm = nPoints

  ! --- Initialize Energy Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_E )

  ASSOCIATE &
    ( minE => OpTab % EnergyGrid % MinValue, &
      maxE => OpTab % EnergyGrid % MaxValue )

  E(:) = 10**( LOG10(minE) + ( LOG10(maxE) - LOG10(minE) ) * rndm_E )

  END ASSOCIATE

  ! --- Initialize Density Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )

  ASSOCIATE &
    ( minD => OpTab % EOSTable % TS % MinValues(iD), &
      maxD => OpTab % EOSTable % TS % MaxValues(iD) )

  D(:) = 10**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  END ASSOCIATE

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  ASSOCIATE &
    ( minT => OpTab % EOSTable % TS % MinValues(iT), &
      maxT => OpTab % EOSTable % TS % MaxValues(iT) )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  END ASSOCIATE

  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Y )

  ASSOCIATE &
    ( minY => OpTab % EOSTable % TS % MinValues(iY), &
      maxY => OpTab % EOSTable % TS % MaxValues(iY) )

  Y(:) = minY + ( maxY - minY ) * rndm_Y

  END ASSOCIATE

  ASSOCIATE &
    ( Etab => OpTab % EnergyGrid % Values,                  &
      Dtab => OpTab % EOSTable % TS % States(iD)  % Values, &
      Ttab => OpTab % EOSTable % TS % States(iT)  % Values, &
      Ytab => OpTab % EOSTable % TS % States(iY)  % Values, &
      Ctab => OpTab % EmAb % Absorptivity(1) % Values, &
      OS   => OpTab % EmAb % Offsets(1) )

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateSingleVariable &
         ( E, D, T, Y, Etab, Dtab, Ttab, Ytab, LogInterp, OS, Ctab, Interpolant1 )

  CALL CPU_TIME( tEND )

  WRITE(*,*)
  WRITE(*,'(A4,A40,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_4D: ', tEND - tBegin
  WRITE(*,*)

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateSingleVariable &
         ( E, D, T, Y, Etab, Dtab, Ttab, Ytab, OS, Ctab, Interpolant2 )

  CALL CPU_TIME( tEND )

  WRITE(*,*)
  WRITE(*,'(A4,A40,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_4D_Custom: ', tEND - tBegin
  WRITE(*,*)

  END ASSOCIATE

  WRITE(*,*) MAXVAL( ABS( Interpolant1 - Interpolant2 ) )

END PROGRAM wlOpacityPerformanceTest

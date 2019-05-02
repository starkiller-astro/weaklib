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
    iP
  INTEGER, DIMENSION(4) :: &
    LogInterp = [ 1, 1, 1, 0 ]
  INTEGER, PARAMETER :: &
    iD = 1, iT = 2, iY = 3, &
    nPoints = 2**20, &
    nPointsX = 2**16, &
    nPointsE = 2**05
  REAL(dp) :: &
    tBegin, &
    tEND
  REAL(dp), DIMENSION(nPointsX) :: &
    D, T, Y, &
    rndm_D, &
    rndm_T, &
    rndm_Y
  REAL(dp), DIMENSION(nPointsE) :: &
    E, LogE, &
    rndm_E
  REAL(dp), DIMENSION(nPointsE,nPointsX) :: &
    Interpolant1, &
    Interpolant2, &
    Interpolant3
  TYPE(OpacityTableType) :: &
    OpTab

  CALL InitializeHDF( )

  CALL ReadOpacityTableHDF &
         ( OpTab, "wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5" )

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

  n_rndm = nPointsE

  ! --- Initialize Energy Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_E )

  ASSOCIATE &
    ( minE => OpTab % EnergyGrid % MinValue, &
      maxE => OpTab % EnergyGrid % MaxValue )

  E(:) = 10**( LOG10(minE) + ( LOG10(maxE) - LOG10(minE) ) * rndm_E )

  END ASSOCIATE

  n_rndm = nPointsX

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
      Ctab => OpTab % EmAb % Opacity(1) % Values, &
      OS   => OpTab % EmAb % Offsets(1) )

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateSingleVariable &
         ( E, D, T, Y, Etab, Dtab, Ttab, Ytab, LogInterp, OS, Ctab, &
           Interpolant1 )

  CALL CPU_TIME( tEND )

  WRITE(*,*)
  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D: ', tEND - tBegin
  WRITE(*,*)

  CALL CPU_TIME( tBegin )

  CALL LogInterpolateSingleVariable &
         ( LOG10( E    ), LOG10( D    ), LOG10( T    ), Y, &
           LOG10( Etab ), LOG10( Dtab ), LOG10( Ttab ), Ytab, &
           OS, Ctab, Interpolant2 )

  CALL CPU_TIME( tEND )

  WRITE(*,*)
  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D_Custom: ', tEND - tBegin
  WRITE(*,*)

  CALL CPU_TIME( tBegin )

  LogE = LOG10( E )

  ASSOCIATE &
    ( LogEtab => LOG10( OpTab % EnergyGrid % Values ), &
      LogDtab => LOG10( OpTab % EOSTable % TS % States(iD)  % Values ), &
      LogTtab => LOG10( OpTab % EOSTable % TS % States(iT)  % Values ) )

  DO iP = 1, nPointsX

    CALL LogInterpolateSingleVariable &
           ( LogE,    LOG10( D(iP) ), LOG10( T(iP) ), Y(iP), &
             LogEtab, LogDtab,        LogTtab,        Ytab,  &
             OS, Ctab, Interpolant3(:,iP) )

  END DO

  END ASSOCIATE ! LogEtab, etc.

  CALL CPU_TIME( tEND )

  WRITE(*,*)
  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D_Custom_Point: ', tEND - tBegin
  WRITE(*,*)

  END ASSOCIATE

  WRITE(*,*) MINVAL( Interpolant1 ), MAXVAL( Interpolant1 )
  WRITE(*,*) MINVAL( Interpolant2 ), MAXVAL( Interpolant2 )
  WRITE(*,*) MINVAL( Interpolant3 ), MAXVAL( Interpolant3 )
  WRITE(*,*) MAXVAL( ABS( Interpolant1 - Interpolant2 ) )
  WRITE(*,*) MAXVAL( ABS( Interpolant1 - Interpolant3 ) )

END PROGRAM wlOpacityPerformanceTest

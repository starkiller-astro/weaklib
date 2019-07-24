PROGRAM wlOpacityAlignTest

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
    LogInterpolateSingleVariable_2D_Custom_Point, &
    LogInterpolateSingleVariable_2D2D_Custom_Point, &
    LogInterpolateSingleVariable_2D2D_Custom_Aligned

  IMPLICIT NONE

  INTEGER :: &
    i, j, k, l, iE, iX, n_rndm, Shape_H1(4)
  INTEGER, PARAMETER :: &
    iD = 1, iT = 2, iY = 3, &
    nPointsX = 2**16, &
    nPointsE = 2**05
  REAL(dp) :: &
    tBegin, tEnd
  REAL(dp), DIMENSION(nPointsE) :: &
    E, LogE
  REAL(dp), DIMENSION(nPointsX) :: &
    T, LogT, &
    Eta, LogEta, &
    rndm_T, &
    rndm_Eta
  REAL(dp), DIMENSION(nPointsE,nPointsE,nPointsX) :: &
    Interpolated_H1_0, &
    Interpolated_H1_1, &
    Interpolated_H1_2
  REAL(dp), ALLOCATABLE :: &
    Tabulated_H1(:,:,:,:), &
    Aligned_H1(:,:,:,:)
  TYPE(OpacityTableType) :: &
    OpTab

  CALL InitializeHDF( )

  CALL ReadOpacityTableHDF &
         ( OpTab, FileName_NES_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5', &
           Verbose_Option = .TRUE. )

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
  WRITE(*,'(A4,A11,I10.10)') &
    '', 'nPointsX = ', nPointsX
  WRITE(*,'(A4,A11,I10.10)') &
    '', 'nPointsE = ', nPointsE
  WRITE(*,*)

  ! --- Energy Points ---

  ASSOCIATE &
    ( minE => OpTab % EnergyGrid % MinValue, &
      maxE => OpTab % EnergyGrid % MaxValue )

  DO iE = 1, nPointsE

    E(iE) = 10**( LOG10(minE) + ( LOG10(maxE) - LOG10(minE) ) &
                                * DBLE(iE) / DBLE(nPointsE+1) )

  END DO

  LogE = LOG10( E )

  END ASSOCIATE

  n_rndm = nPointsX

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  ASSOCIATE &
    ( minT => OpTab % EOSTable % TS % MinValues(iT), &
      maxT => OpTab % EOSTable % TS % MaxValues(iT) )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  LogT = LOG10( T )

  END ASSOCIATE

  ! --- Initialize Degeneracy Parameter Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Eta )

  ASSOCIATE &
    ( minEta => MINVAL( OpTab % EtaGrid % Values ), &
      maxEta => MAXVAL( OpTab % EtaGrid % Values ) )

  Eta(:) = 10**( LOG10(minEta) + ( LOG10(maxEta) - LOG10(minEta) ) * rndm_Eta )

  LogEta = LOG10( Eta )

  END ASSOCIATE

  ASSOCIATE &
    ( LogEs   => LOG10( OpTab % EnergyGrid      % Values ), &
         Es   =>        OpTab % EnergyGrid      % Values  , &
      LogTs   => LOG10( OpTab % TS % States(iT) % Values ), &
         Ts   =>        OpTab % TS % States(iT) % Values  , &
      LogEtas => LOG10( OpTab % EtaGrid         % Values ), &
         Etas =>        OpTab % EtaGrid         % Values  , &
      OS      => OpTab % Scat_NES % Offsets(1,1) )

  Shape_H1 = SHAPE( OpTab % Scat_NES % Kernel(1) % Values(:,:,1,:,:) )

  ALLOCATE( Tabulated_H1(Shape_H1(1),Shape_H1(2),Shape_H1(3),Shape_H1(4)) )

  Tabulated_H1 = OpTab % Scat_NES % Kernel(1) % Values(:,:,1,:,:)

  CALL CPU_TIME( tBegin )

  Interpolated_H1_0 = 0.0d0

  DO iX = 1, nPointsX

    CALL LogInterpolateSingleVariable_2D2D_Custom_Point &
           ( LogE, LogT(iX), LogEta(iX), LogEs, LogTs, LogEtas, &
             OS, Tabulated_H1, Interpolated_H1_0(:,:,iX) )

  END DO

  CALL CPU_TIME( tEnd )

  WRITE(*,*)
  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'Interpolate H1 (Original): ', tEnd - tBegin

  ALLOCATE( Aligned_H1(SIZE(LogTs),SIZE(LogEtas),nPointsE,nPointsE) )

  CALL CPU_TIME( tBegin )

  DO l = 1, SIZE( LogEtas )
  DO k = 1, SIZE( LogTs   )
  DO j = 1, SIZE( E )
  DO i = 1, SIZE( E )

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( E(i), E(j), Es, Es, OS, Tabulated_H1(:,:,k,l), Aligned_H1(k,l,i,j) )

    Aligned_H1(k,l,i,j) = LOG10( Aligned_H1(k,l,i,j) + OS )

  END DO
  END DO
  END DO
  END DO

  CALL CPU_TIME( tEnd )

  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'Align H1: ', tEnd - tBegin

  ! --- Initial Test of Aligned Interpolation ---

  CALL CPU_TIME( tBegin )

  Interpolated_H1_1 = 0.0d0

  DO j = 1, nPointsE
  DO i = 1, j

    DO iX = 1, nPointsX

      CALL LogInterpolateSingleVariable_2D_Custom_Point &
             ( T(iX), Eta(iX), Ts, Etas, OS, Aligned_H1(:,:,i,j), &
               Interpolated_H1_1(i,j,iX) )

    END DO

  END DO
  END DO

  CALL CPU_TIME( tEnd )

  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'Interpolate H1 (Aligned 1): ', tEnd - tBegin
  WRITE(*,*)

  PRINT*, "  diff 1 = ", MAXVAL( ABS( Interpolated_H1_0 - Interpolated_H1_1 ) )
  PRINT*, "       0 = ", Interpolated_H1_0(1,1,1)
  PRINT*, "       1 = ", Interpolated_H1_1(1,1,1)

  ! --- Aligned Interpolation ---

  CALL CPU_TIME( tBegin )

  Interpolated_H1_2 = 0.0d0

  CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
         ( LogT, LogEta, LogTs, LogEtas, OS, Aligned_H1, Interpolated_H1_2 )

  CALL CPU_TIME( tEnd )

  WRITE(*,'(A4,A48,ES10.4E2)') &
    '', 'Interpolate H1 (Aligned 2): ', tEnd - tBegin
  WRITE(*,*)

  print*, "  diff 2 = ", MAXVAL( ABS( Interpolated_H1_0 - Interpolated_H1_2 ) )
  PRINT*, "       0 = ", Interpolated_H1_0(1,1,1)
  PRINT*, "       1 = ", Interpolated_H1_2(1,1,1)

  END ASSOCIATE ! LogEs, etc.

  DEALLOCATE( Tabulated_H1 )
  DEALLOCATE( Aligned_H1 )

END PROGRAM wlOpacityAlignTest

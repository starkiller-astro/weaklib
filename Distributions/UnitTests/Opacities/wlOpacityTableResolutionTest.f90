PROGRAM wlOpacityTableResolutionTest

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateSingleVariable_1D3D_Custom_Point, &
    LogInterpolateSingleVariable_2D2D_Custom_Point
  USE wlOpacityFieldsModule, ONLY: &
    iHi0, iHii0, iHi1, iHii1
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    AllocateOpacityTable, &
    DeAllocateOpacityTable, &
    DescribeOpacityTable
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF, &
    OpenFileHDF, &
    CloseFileHDF, &
    WriteHDF, &
    OpenGroupHDF, &
    CloseGroupHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF, &
    WriteOpacityTableHDF
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DescribeGrid, &
    MakeLogGrid
  USE wlExtPhysicalConstantsModule, ONLY: kMeV, ca, cv
  USE wlExtNumericalModule, ONLY: pi, half, twpi, zero
  USE HDF5

  IMPLICIT NONE
 
  !-------- Table filename --------------------------------------------------
  CHARACTER(256) :: HighResOpTableBase   = "wl-Op-SFHo-25-40-100-E40-B85"
  CHARACTER(256) :: LowResOpTableBase    = "wl-Op-SFHo-25-20-100-E40-B85"
  CHARACTER(256) :: HightResEOSTableName = "wl-EOS-SFHo-25-40-100.h5"
  INTEGER, DIMENSION(4) :: TableFlags = [1,1,1,1] ! [EmAb, Iso, NES, Pair]

  !-------- variables for reading opacity table -----------------------------
  TYPE(OpacityTableType) :: OpacityTableHi, OpacityTableLo, ErrorTable
  REAL(dp), DIMENSION(2) :: Offset_NES

  !-------- variables for reading parameters data ---------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: Inte_rho, Inte_T, &
                                             Inte_Ye, Inte_cmpe, database
  CHARACTER(LEN=100)                      :: Format1, Format2, Format3
  CHARACTER(LEN=30)                       :: a
  INTEGER                                 :: i, datasize, icmpe

  !-------- variables for output -------------------------------------------
  INTEGER(HID_T)                          :: file_id, group_id
  INTEGER(HSIZE_T)                        :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d
  INTEGER(HSIZE_T), DIMENSION(3)          :: datasize3d

  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: SumNES_nue, SumNES_nuebar, &
                                             SumNES_mutau, SumNES_mutaubar

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: InterpolantNES_nue, &
                                             InterpolantNES_nuebar, &
                                             InterpolantNES_mutau, &
                                             InterpolantNES_mutaubar
  
  CHARACTER(256)                          :: WriteTableName

  !-------- local variables -------------------------------------------------
  CHARACTER(128)                          :: FileNameHi(4), FileNameLo(4)

  INTEGER                                 :: ii, jj
  INTEGER                                 :: j_rho, k_t, l_ye, i_eta, i_r, t_m
  INTEGER, DIMENSION(4)                   :: LogInterp
  REAL(dp)                                :: cparpe  = (cv+ca)**2
  REAL(dp)                                :: cparne  = (cv-ca)**2
  REAL(dp)                                :: cparpmt = (cv+ca-2.d0)**2
  REAL(dp)                                :: cparnmt = (cv-ca)**2
  REAL(dp)                                :: sum_NES_nue, sum_NES_nuebar, &
                                             sum_NES_mutau, sum_NES_mutaubar, &
                                             root2p, root2n
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: NES0_nue, NES0_mutau 
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: NES0_nuebar, NES0_mutaubar
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: roots, widths

  INTEGER                                 :: nOpac_EmAb
  INTEGER                                 :: nOpac_Iso 
  INTEGER                                 :: nMom_Iso  
  INTEGER                                 :: nOpac_NES 
  INTEGER                                 :: nMom_NES  
  INTEGER                                 :: nOpac_Pair
  INTEGER                                 :: nMom_Pair 
  INTEGER                                 :: nPointsE, nPointsEta

  REAL(dp)                                :: rho, T, ye, eta
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: TableValue1D
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: InterEmAb, InterIso
  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: TableValue2D
  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: InterH0i, InterH0ii
  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: InterJ0i, InterJ0ii

  REAL(dp) :: OffRatio

  OffRatio = 1.0_dp - 1.0d-1

  IF( TableFlags(1) == 1 )THEN
    FileNameHi(1) = TRIM(HighResOpTableBase)//'-EmAb.h5'
    FileNameLo(1) = TRIM(LowResOpTableBase)//'-EmAb.h5'
    nOpac_EmAb = 2
  ELSE
    FileNameHi(1) = ''
    FileNameLo(1) = ''
    nOpac_EmAb = 0
  END IF

  IF( TableFlags(2) == 1 )THEN
    FileNameHi(2) = TRIM(HighResOpTableBase)//'-Iso.h5'
    FileNameLo(2) = TRIM(LowResOpTableBase)//'-Iso.h5'
    nOpac_Iso = 2
    nMom_Iso  = 2
  ELSE
    FileNameHi(2) = ''
    FileNameLo(2) = ''
    nOpac_Iso = 0
    nMom_Iso  = 0
  END IF

  IF( TableFlags(3) == 1 )THEN
    FileNameHi(3) = TRIM(HighResOpTableBase)//'-NES.h5'
    FileNameLo(3) = TRIM(LowResOpTableBase)//'-NES.h5'
    nOpac_NES = 1
    nMom_NES  = 4
  ELSE
    FileNameHi(3) = ''
    FileNameLo(3) = ''
    nOpac_NES = 0
    nMom_NES  = 0
  END IF

  IF( TableFlags(4) == 1 )THEN
    FileNameHi(4) = TRIM(HighResOpTableBase)//'-Pair.h5'
    FileNameLo(4) = TRIM(LowResOpTableBase)//'-Pair.h5'
    nOpac_Pair = 1
    nMom_Pair  = 4
  ELSE
    FileNameHi(4) = ''
    FileNameLo(4) = ''
    nOpac_Pair = 0
    nMom_Pair  = 0
  END IF

  !------------------------------------------------------
  ! read in the reference high resolution opacity table
  !------------------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTableHi,             &
       FileName_EmAB_Option                             &
       = TRIM(FileNameHi(1)),                           &
       FileName_Iso_Option                              &
       = TRIM(FileNameHi(2)),                           &
       FileName_NES_Option                              &
       = TRIM(FileNameHi(3)),                           &
       FileName_Pair_Option                             &
       = TRIM(FileNameHi(4)),                           &
       EquationOfStateTableName_Option                  &
       = TRIM(HightResEOSTableName),                    &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )
  
  !------------------------------------------------------
  !    read in the low resoltion opacity table
  !------------------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTableLo,             &
       FileName_EmAB_Option                             &
       = TRIM(FileNameLo(1)),                           &
       FileName_Iso_Option                              &
       = TRIM(FileNameLo(2)),                           &
       FileName_NES_Option                              &
       = TRIM(FileNameLo(3)),                           &
       FileName_Pair_Option                             &
       = TRIM(FileNameLo(4)),                           &
       EquationOfStateTableName_Option                  &
       = TRIM(HightResEOSTableName),                    &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )

  !------------------------------------------------------
  !   create an empty table for interpolate error
  !   it has same resolution as the high-res. table
  !------------------------------------------------------
  nPointsE = OpacityTableHi % nPointsE
  nPointsEta = OpacityTableHi % nPointsEta

  CALL InitializeHDF( )
  CALL AllocateOpacityTable( ErrorTable,                &
       nOpac_EmAb, nOpac_Iso, nMom_Iso, &
       nOpac_NES, nMom_NES, nOpac_Pair, nMom_Pair, &
       nPointsE, nPointsEta, &
       EquationOfStateTableName_Option &
       = TRIM(HightResEOSTableName) )
  CALL FinalizeHDF( )

   ! set up ErrorTable % EnergyGrid
   ErrorTable % EnergyGrid = OpacityTableHi % EnergyGrid
   ErrorTable % EnergyGrid % maxValue  &
                           = OpacityTableHi % EnergyGrid % Values(nPointsE) * OffRatio
   ErrorTable % EnergyGrid % Values(nPointsE)    &
                           = OpacityTableHi % EnergyGrid % Values(nPointsE) * OffRatio

   ! set up ErrorTable % EtaGrid
   IF( nPointsEta > 0 )THEN
     ErrorTable % EtaGrid = OpacityTableHi % EtaGrid
     ErrorTable % EtaGrid % maxValue  &
                             = OpacityTableHi % EtaGrid % Values(nPointsEta) * OffRatio
     ErrorTable % EtaGrid % Values(nPointsEta)    &
                             = OpacityTableHi % EtaGrid % Values(nPointsEta) * OffRatio
   END IF
   ! set up OpacityTableTypeEmAb
   IF( nOpac_EmAb .gt. 0 ) THEN
     ErrorTable % EmAb % nOpacities = OpacityTableHi % EmAb % nOpacities
     ErrorTable % EmAb % nPoints    = OpacityTableHi % EmAb % nPoints
     ErrorTable % EmAb % Names      = OpacityTableHi % EmAb % Names
     ErrorTable % EmAb % Units      = (/'DIMENSIONLESS','DIMENSIONLESS'/)
   END IF

   ! set up OpacityTableTypeScat Iso
   IF( nOpac_Iso .gt. 0 ) THEN
     ErrorTable % Scat_Iso % nOpacities= OpacityTableHi % Scat_Iso % nOpacities
     ErrorTable % Scat_Iso % nMoments  = OpacityTableHi % Scat_Iso % nMoments
     ErrorTable % Scat_Iso % nPoints   = OpacityTableHi % Scat_Iso % nPoints
     ErrorTable % Scat_Iso % Names     = OpacityTableHi % Scat_Iso % Names
     ErrorTable % Scat_Iso % Units     = (/'DIMENSIONLESS','DIMENSIONLESS'/)
   END IF

   ! set up OpacityTableTypeScat NES
   IF( nOpac_NES .gt. 0 ) THEN
     ErrorTable % Scat_NES % nOpacities = OpacityTableHi % Scat_NES % nOpacities
     ErrorTable % Scat_NES % nMoments   = OpacityTableHi % Scat_NES % nMoments
     ErrorTable % Scat_NES % nPoints    = OpacityTableHi % Scat_NES % nPoints
     ErrorTable % Scat_NES % Names      = OpacityTableHi % Scat_NES % Names
     ErrorTable % Scat_NES % Units      = (/'DIMENSIONLESS'/)
   END IF

   ! set up OpacityTableTypeScat Pair
   IF( nOpac_Pair .gt. 0 ) THEN
     ErrorTable % Scat_Pair % nOpacities = OpacityTableHi % Scat_Pair % nOpacities
     ErrorTable % Scat_Pair % nMoments   = OpacityTableHi % Scat_Pair % nMoments
     ErrorTable % Scat_Pair % nPoints    = OpacityTableHi % Scat_Pair % nPoints
     ErrorTable % Scat_Pair % Names      = OpacityTableHi % Scat_Pair % Names
     ErrorTable % Scat_Pair % Units      = (/'DIMENSIONLESS'/)
   END IF

   ! allocate local variables
   ALLOCATE( TableValue1D(nPointsE) )
   ALLOCATE(    InterEmAb(nPointsE) )
   ALLOCATE(    InterIso (nPointsE) )
   ALLOCATE( TableValue2D(nPointsE,nPointsE) )
   ALLOCATE(     InterH0i(nPointsE,nPointsE) )
   ALLOCATE(    InterH0ii(nPointsE,nPointsE) )
   ALLOCATE(     InterJ0i(nPointsE,nPointsE) )
   ALLOCATE(    InterJ0ii(nPointsE,nPointsE) )

  !----------------------------------------------------------------------------
  !   do interpolation for high-res. at low-res. table
  !----------------------------------------------------------------------------
  WRITE(*,*) 'Interpolating on low-resolution table ...'

  ASSOCIATE &
  ( iEOS_Rho      => OpacityTableLo % EOSTable % TS % Indices % iRho, &
    iEOS_T        => OpacityTableLo % EOSTable % TS % Indices % iT,   &
    iEOS_Ye       => OpacityTableLo % EOSTable % TS % Indices % iYe,  &
    iRho          => OpacityTableLo % TS % Indices % iRho, &
    iT            => OpacityTableLo % TS % Indices % iT,   &
    iYe           => OpacityTableLo % TS % Indices % iYe,  &
    LogInterp     => OpacityTableLo % EOSTable % TS % LogInterp )

  IF( TableFlags(1) + TableFlags(2) .ge. 1 )THEN

    WRITE(*,*) 'EmAb and Iso interpolating ... '

    DO l_ye = 1, ErrorTable % nPointsTS(iYe) - 1
      ye = ErrorTable % TS % States (iYe) % Values (l_ye)
      DO k_t = 1, ErrorTable % nPointsTS(iT) - 1
        T = ErrorTable % TS % States (iT) % Values (k_t)
        DO j_rho = 1, ErrorTable % nPointsTS(iRho) - 1
          rho = ErrorTable % TS % States (iRho) % Values (j_rho)
          IF( TableFlags(1) == 1 )THEN
            !--- EmAb Interpolation ---
            DO i_r = 1, nOpac_EmAb
              ! table value in high-res. table
              TableValue1D = OpacityTableHi % EmAb % Opacity(i_r) % &
                             Values(:,j_rho,k_t,l_ye)
              TableValue1D = 10**TableValue1D &
                             - OpacityTableHi % EmAb % Offsets(i_r)
              ! interpolate @ low-res. table
              CALL LogInterpolateSingleVariable_1D3D_Custom_Point          &
                   ( LOG10( ErrorTable % EnergyGrid % Values ),            &
                     LOG10( rho ), LOG10( T ), ye,                         &
                     LOG10( OpacityTableLo % EnergyGrid % Values ),        &
                     LOG10( OpacityTableLo % TS % States(iRho) % Values ), &
                     LOG10( OpacityTableLo % TS % States(iT) % Values ),   &
                     OpacityTableLo % TS % States(iYe) % Values,           &
                     OpacityTableLo % EmAb % Offsets(i_r),                 &
                     OpacityTableLo % EmAb % Opacity(i_r) % Values,        &
                     InterEmAb )
              ErrorTable % EmAb % Opacity(i_r) % Values(:,j_rho,k_t,l_ye)  &
                = LOG10( MAX( ABS( TableValue1D - InterEmAb ) &
                              / MAX(TableValue1D, 1.0d-30), 1.0d-30) )
            END DO
          END IF
          IF( TableFlags(2) == 1 )THEN
            !--- Iso Interpolation ---
            DO i_r = 1, nOpac_Iso
              DO t_m = 1, nMom_Iso
                ! table value in high-res. table
                TableValue1D = OpacityTableHi % Scat_Iso % Kernel(i_r) %      &
                               Values(:,t_m,j_rho,k_t,l_ye)
                TableValue1D = 10**TableValue1D &
                               - OpacityTableHi % Scat_Iso % Offsets(i_r,t_m)
                ! interpolate @ low-res. table
                CALL LogInterpolateSingleVariable_1D3D_Custom_Point           &
                     ( LOG10( ErrorTable % EnergyGrid % Values ),             &
                       LOG10( rho ), LOG10( T ), ye,                          &
                       LOG10( OpacityTableLo % EnergyGrid % Values ),         &
                       LOG10( OpacityTableLo % TS % States(iRho) % Values ),  &
                       LOG10( OpacityTableLo % TS % States(iT) % Values ),    &
                       OpacityTableLo % TS % States(iYe) % Values,            &
                       OpacityTableLo % Scat_Iso % Offsets(i_r,t_m),          &
                       OpacityTableLo % Scat_Iso % Kernel(i_r)                &
                       % Values(:,t_m,:,:,:),                                 &
                       InterIso )
                ErrorTable % Scat_Iso % Kernel(i_r) % &
                  Values(:,t_m,j_rho,k_t,l_ye)        &
                  = LOG10( MAX( ABS( TableValue1D - InterIso ) &
                                / MAX(TableValue1D, 1.0d-30), 1.0d-30) )
              END DO ! nMom
            END DO ! nOpac
          END IF
        END DO ! rho
      END DO ! T
    END DO ! ye
  
  END IF  
  
  IF( TableFlags(3) .eq. 1 )THEN

    WRITE(*,*) ' NES interpolating ... '

    DO i_eta = 1, ErrorTable % nPointsEta
      eta = ErrorTable % EtaGrid % Values(i_eta)
      DO k_t = 1, ErrorTable % nPointsTS(iT) - 1
        T    = ErrorTable % TS % States (iT) % Values (k_t)
        DO i_r = 1, nOpac_NES
          DO t_m = 1, nMom_NES
            ! table value in high-res. table
            TableValue2D = OpacityTableHi % Scat_NES % Kernel(i_r) % Values &
                         (:,:,t_m,k_t,i_eta)
            TableValue2D = 10**TableValue2D &
                        - OpacityTableHi % Scat_NES % Offsets(i_r,t_m)
            ! interpolate @ low-res. table
              CALL LogInterpolateSingleVariable_2D2D_Custom_Point        &
                     ( LOG10(ErrorTable % EnergyGrid % Values),          &
                       LOG10(T), LOG10(eta),                             &
                       LOG10(OpacityTableLo % EnergyGrid % Values),      &
                       LOG10(OpacityTableLo % TS % States(iT) % Values), &
                       LOG10(OpacityTableLo % EtaGrid % Values),         &
                       OpacityTableLo % Scat_NES % Offsets(i_r,t_m),     &
                       OpacityTableLo % Scat_NES % Kernel(i_r) %         &
                       Values(:,:,t_m,:,:), InterH0i )
              ErrorTable % Scat_NES % Kernel(i_r) % Values     &
                       ( :, :, t_m, k_t, i_eta )               &
                  = LOG10( MAX( ABS( TableValue2D - InterH0i ) &
                                / MAX(TableValue2D, 1.0d-30), 1.0d-30) )
          END DO
        END DO
      END DO
    END DO

  END IF

  IF( TableFlags(4) .eq. 1 )THEN

    WRITE(*,*) ' Pair interpolating ... '

    DO i_eta = 1, ErrorTable % nPointsEta
      eta = ErrorTable % EtaGrid % Values(i_eta)
      DO k_t = 1, ErrorTable % nPointsTS(iT) - 1
        T    = ErrorTable % TS % States (iT) % Values (k_t)
        DO i_r = 1, nOpac_Pair
          DO t_m = 1, nMom_Pair
            ! table value in high-res. table
            TableValue2D = OpacityTableHi % Scat_Pair % Kernel(i_r) % Values &
                         (:,:,t_m,k_t,i_eta)
            TableValue2D = 10**TableValue2D &
                        - OpacityTableHi % Scat_Pair % Offsets(i_r,t_m)
            ! interpolate @ low-res. table
              CALL LogInterpolateSingleVariable_2D2D_Custom_Point        &
                     ( LOG10(ErrorTable % EnergyGrid % Values),          &
                       LOG10(T), LOG10(eta),                             &
                       LOG10(OpacityTableLo % EnergyGrid % Values),      &
                       LOG10(OpacityTableLo % TS % States(iT) % Values), &
                       LOG10(OpacityTableLo % EtaGrid % Values),         &
                       OpacityTableLo % Scat_Pair % Offsets(i_r,t_m),    &
                       OpacityTableLo % Scat_Pair % Kernel(i_r) %        &
                       Values(:,:,t_m,:,:), InterJ0i )
              ErrorTable % Scat_Pair % Kernel(i_r) % Values    &
                       ( :, :, t_m, k_t, i_eta )               &
                  = LOG10( MAX( ABS( TableValue2D - InterJ0i ) &
                                /MAX(TableValue2D, 1.0d-30), 1.0d-30) )
          END DO
        END DO
      END DO
    END DO

  END IF

  END ASSOCIATE

  CALL DescribeOpacityTable( ErrorTable )

  IF( nOpac_EmAb > 0 ) THEN
    WriteTableName = 'InterpolatedError_EmAb.h5'
    CALL InitializeHDF( )
    WRITE(*,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_EmAb_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_Iso > 0 ) THEN
    WriteTableName = 'InterpolatedError_Iso.h5'
    CALL InitializeHDF( )
    WRITE(*,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_Iso_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_NES > 0 ) THEN
    WriteTableName = 'InterpolatedError_NES.h5'
    CALL InitializeHDF( )
    WRITE(*,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_NES_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_Pair > 0 ) THEN
    WriteTableName = 'InterpolatedError_Pair.h5'
    CALL InitializeHDF( )
    WRITE(*,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_Pair_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  CALL DeAllocateOpacityTable( OpacityTableHi )
  CALL DeAllocateOpacityTable( OpacityTableLo )
  CALL DeAllocateOpacityTable( ErrorTable )
  DEALLOCATE( InterEmAb, InterIso, TableValue1D )
  DEALLOCATE( InterH0i, InterH0ii, InterJ0i, InterJ0ii, TableValue2D )

END PROGRAM wlOpacityTableResolutionTest

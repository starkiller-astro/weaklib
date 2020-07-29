PROGRAM wlOpacityTableResolutionTest

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateSingleVariable_1D3D_Custom_Point, &
    LogInterpolateSingleVariable_2D2D_Custom_Point, &
    LinearInterp_Array_Point, &
    GetIndexAndDelta    
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
  USE, INTRINSIC :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  IMPLICIT NONE

  !-------- Table filename --------------------------------------------------
  !CHARACTER(256) :: HighResOpTableBase   = "wl-Op-SFHo-25-40-100-E40-B85"
  !CHARACTER(256) :: LowResOpTableBase    = "wl-Op-SFHo-25-20-100-E40-B85"

  CHARACTER(256) :: HighResOpTableBase     = "wl-Op-SFHo-25-40-100-E40-B85"
  CHARACTER(256) :: LowResOpTableBase      = "wl-Op-SFHo-15-25-50-E40-B85"

  CHARACTER(256) :: HighResOpTableBaseBrem = "wl-Op-SFHo-25-40-100-E40-HR98"
  CHARACTER(256) :: LowResOpTableBaseBrem  = "wl-Op-SFHo-15-25-50-E40-HR98"

  CHARACTER(256) :: HighResEOSTableName    = "wl-EOS-SFHo-25-40-100.h5"

  INTEGER, DIMENSION(5) :: TableFlags      = [1,1,1,1,1] ! [EmAb, Iso, NES, Pair, Brem]

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
  CHARACTER(128)                          :: FileNameHi(5), FileNameLo(5)

  INTEGER                                 :: ii, jj
  INTEGER                                 :: k, kp
  INTEGER                                 :: idxRho, idxT, idxYe, idxEta 
  REAL(dp)                                :: dRho, dT, dYe, dEta
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
  INTEGER                                 :: nOpac_Brem
  INTEGER                                 :: nMom_Brem 
  INTEGER                                 :: nPointsE, nPointsEta

  REAL(dp)                                :: rho, T, ye, eta
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: TableValue1D
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: InterEmAb, InterIso
  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: TableValue2D
  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: InterH0i, InterH0ii
  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: InterJ0i, InterJ0ii

  REAL(dp)   :: InterpValue
  REAL(dp)   :: ReferenceValue

  REAL(dp)   :: SMAPE_pt_max !maximum pointwise symmetric mean absolute percentage error
  REAL(dp)   :: max_rel_e_Ep, max_rel_e_E, max_rel_e_rho, max_rel_e_T, max_rel_e_Ye, max_rel_e_eta
  REAL(dp)   :: SMAPE    !symmetric mean absolute percentage error

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

  IF( TableFlags(5) == 1 )THEN
    FileNameHi(5) = TRIM(HighResOpTableBaseBrem)//'-Brem.h5'
    FileNameLo(5) = TRIM(LowResOpTableBaseBrem)//'-Brem.h5'
    nOpac_Brem = 1
    nMom_Brem  = 1
  ELSE
    FileNameHi(5) = ''
    FileNameLo(5) = ''
    nOpac_Brem = 0
    nMom_Brem  = 0
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
       FileName_Brem_Option                             &
       = TRIM(FileNameHi(5)),                           &
       EquationOfStateTableName_Option                  &
       = TRIM(HighResEOSTableName),                     &
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
       FileName_Brem_Option                             &
       = TRIM(FileNameLo(5)),                           &
       EquationOfStateTableName_Option                  &
       = TRIM(HighResEOSTableName),                     &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )

  !------------------------------------------------------
  !   create an empty table for interpolate error
  !   it has same resolution as the high-res. table
  !------------------------------------------------------
  nPointsE   = OpacityTableHi % nPointsE
  nPointsEta = OpacityTableHi % nPointsEta

  CALL InitializeHDF( )
  CALL AllocateOpacityTable( ErrorTable,                &
       nOpac_EmAb, nOpac_Iso, nMom_Iso, &
       nOpac_NES, nMom_NES, nOpac_Pair, nMom_Pair, nOpac_Brem, nMom_Brem, &
       nPointsE, nPointsEta, &
       EquationOfStateTableName_Option &
       = TRIM(HighResEOSTableName) )
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

   ! set up OpacityTableTypeBrem 
   IF( nOpac_Brem .gt. 0 ) THEN
     ErrorTable % Scat_Brem % nOpacities = OpacityTableHi % Scat_Brem % nOpacities
     ErrorTable % Scat_Brem % nMoments   = OpacityTableHi % Scat_Brem % nMoments
     ErrorTable % Scat_Brem % nPoints    = OpacityTableHi % Scat_Brem % nPoints
     ErrorTable % Scat_Brem % Names      = OpacityTableHi % Scat_Brem % Names
     ErrorTable % Scat_Brem % Units      = (/'DIMENSIONLESS'/)
   END IF

   ! allocate local variables
!   ALLOCATE( TableValue1D(nPointsE) )
!   ALLOCATE(    InterEmAb(nPointsE) )
!   ALLOCATE(    InterIso (nPointsE) )
!   ALLOCATE( TableValue2D(nPointsE,nPointsE) )
!   ALLOCATE(     InterH0i(nPointsE,nPointsE) )
!   ALLOCATE(    InterH0ii(nPointsE,nPointsE) )
!   ALLOCATE(     InterJ0i(nPointsE,nPointsE) )
!   ALLOCATE(    InterJ0ii(nPointsE,nPointsE) )

  !----------------------------------------------------------------------------
  !   do interpolation for high-res. at low-res. table
  !----------------------------------------------------------------------------
  WRITE(stdout,*) 'Checking table resolution: Interpolate from LoRes table at coordinates of HiRes table and compare'
  FLUSH(stdout)

  ASSOCIATE &
  ( iEOS_Rho      => OpacityTableLo % EOSTable % TS % Indices % iRho, &
    iEOS_T        => OpacityTableLo % EOSTable % TS % Indices % iT,   &
    iEOS_Ye       => OpacityTableLo % EOSTable % TS % Indices % iYe,  &
    iRho          => OpacityTableLo % TS % Indices % iRho, &
    iT            => OpacityTableLo % TS % Indices % iT,   &
    iYe           => OpacityTableLo % TS % Indices % iYe,  &
    LogInterp     => OpacityTableLo % EOSTable % TS % LogInterp )

  IF( TableFlags(1) .eq. 1 )THEN
  
    WRITE(stdout,*) ' Checking EmAb, Opacity(1)... '

    WRITE(stdout,'(A,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % EmAb % Opacity(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % EmAb % Opacity(1) % Values)

    SMAPE_pt_max = 0.0d0
    SMAPE        = 0.0d0
    DO i_r = 1, 1 !nOpac_Iso
      DO k = 1, OpacityTableHi % nPointsE
        DO l_ye = 1, OpacityTableHi % nPointsTS(iYe) - 1
          DO k_t = 1, OpacityTableHi % nPointsTS(iT) - 1
            DO j_rho = 1, OpacityTableHi % nPointsTS(iRho) - 1

              rho = OpacityTableHi % TS % States (iRho) % Values (j_rho)
              T   = OpacityTableHi % TS % States (iT) % Values (k_t)  
              Ye  = OpacityTableHi % TS % States (iYe) % Values (l_ye)  

              ReferenceValue = 10.d0**(OpacityTableHi % EmAb % Opacity(i_r) % Values &
                               (k,j_rho,k_t,l_ye)) - OpacityTableHi % EmAb % Offsets(i_r)

                CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableLo % TS % States (iRho) % Values), idxRho, dRho )
                CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT, dT )
                CALL GetIndexAndDelta( Ye,         OpacityTableLo % TS % States (iYe) % Values, idxYe, dYe )

                InterpValue = LinearInterp_Array_Point( k, idxRho, idxT, idxYe, dRho, dT, dYe, &
                              OpacityTableLo % EmAb % Offsets(i_r),       &
                              OpacityTableLo % EmAb % Opacity(i_r) % Values(:,:,:,:))

                IF((ABS(ReferenceValue)+ABS(InterpValue)) == 0.0d0) THEN
                  ErrorTable % EmAb % Opacity(i_r) % Values( k, j_rho, k_t, l_ye ) = 0.0d0
                ELSE
                  ErrorTable % EmAb % Opacity(i_r) % Values( k, j_rho, k_t, l_ye )               &
                    = 2.0d0 *ABS(ReferenceValue - InterpValue)/(ABS(ReferenceValue)+ABS(InterpValue))
                END IF

                IF(ErrorTable % EmAb % Opacity(i_r) % Values( k, j_rho, k_t, l_ye ) .gt. SMAPE_pt_max) THEN
                  SMAPE_pt_max  = ErrorTable % EmAb % Opacity(i_r) % Values( k, j_rho, k_t, l_ye )
                  max_rel_e_E   = OpacityTableHi % EnergyGrid % Values(k)
                  max_rel_e_rho = rho
                  max_rel_e_T   = T
                  max_rel_e_Ye  = Ye
                END IF
      
                SMAPE = SMAPE + ErrorTable % EmAb % Opacity(i_r) % Values( k, j_rho, k_t, l_ye )

              END DO
            END DO
          END DO
        END DO

      SMAPE = SMAPE / SIZE(ErrorTable % EmAb % Opacity(1) % Values)

      WRITE(stdout,'(A,5ES17.4)') 'EmAb table LoRes maximum SMAPE_pt, E, rho, T, Ye', & 
                  SMAPE_pt_max, max_rel_e_E, max_rel_e_rho, max_rel_e_T, max_rel_e_Ye
      WRITE(stdout,'(A,ES17.4,I4)') 'EmAb table LoRes SMAPE, Moment', SMAPE, t_m

    END DO

    ErrorTable % EmAb % Opacity(2) % Values( :, :, :, : ) = 0.0d0

  END IF

  FLUSH(stdout)

  IF( TableFlags(2) .eq. 1 )THEN
  
    WRITE(stdout,*) ' Checking Scat_Iso, Opacity(1) ... '

    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % Scat_Iso % Kernel(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % Scat_Iso % Kernel(1) % Values)

    DO t_m = 1, nMom_Iso 
      SMAPE_pt_max = 0.0d0
      SMAPE        = 0.0d0
      DO i_r = 1, 1 !nOpac_Iso
        DO k = 1, OpacityTableHi % nPointsE
          DO l_ye = 1, OpacityTableHi % nPointsTS(iYe) - 1
            DO k_t = 1, OpacityTableHi % nPointsTS(iT) - 1
              DO j_rho = 1, OpacityTableHi % nPointsTS(iRho) - 1

                rho = OpacityTableHi % TS % States (iRho) % Values (j_rho)
                T   = OpacityTableHi % TS % States (iT) % Values (k_t)  
                Ye  = OpacityTableHi % TS % States (iYe) % Values (l_ye)  

                ReferenceValue = 10.d0**(OpacityTableHi % Scat_Iso % Kernel(i_r) % Values &
                                 (k,t_m,j_rho,k_t,l_ye)) - OpacityTableHi % Scat_Iso % Offsets(i_r,t_m)

                CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableLo % TS % States (iRho) % Values), idxRho, dRho )
                CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT, dT )
                CALL GetIndexAndDelta( Ye,         OpacityTableLo % TS % States (iYe) % Values, idxYe, dYe )

                InterpValue = LinearInterp_Array_Point( k, idxRho, idxT, idxYe, dRho, dT, dYe, &
                              OpacityTableLo % Scat_Iso % Offsets(i_r,t_m),       &
                              OpacityTableLo % Scat_Iso % Kernel(i_r) % Values(:,t_m,:,:,:))

                IF((ABS(ReferenceValue)+ABS(InterpValue)) == 0.0d0) THEN
                  ErrorTable % Scat_Iso % Kernel(i_r) % Values( k, t_m, j_rho, k_t, l_ye ) = 0.0d0
                ELSE
                  ErrorTable % Scat_Iso % Kernel(i_r) % Values( k, t_m, j_rho, k_t, l_ye )               &
                    = 2.0d0 *ABS(ReferenceValue - InterpValue)/(ABS(ReferenceValue)+ABS(InterpValue))
                END IF

                IF(ErrorTable % Scat_Iso % Kernel(i_r) % Values( k, t_m, j_rho, k_t, l_ye ) .gt. SMAPE_pt_max) THEN
                  SMAPE_pt_max  = ErrorTable % Scat_Iso % Kernel(i_r) % Values( k, t_m, j_rho, k_t, l_ye )
                  max_rel_e_E   = OpacityTableHi % EnergyGrid % Values(k)
                  max_rel_e_rho = rho
                  max_rel_e_T   = T
                  max_rel_e_Ye  = Ye
                END IF
      
                SMAPE = SMAPE + ErrorTable % Scat_Iso % Kernel(i_r) % Values( k, t_m, j_rho, k_t, l_ye )

              END DO
            END DO
          END DO
        END DO
      END DO

      SMAPE = SMAPE / SIZE(ErrorTable % Scat_Iso % Kernel(1) % Values) * nMom_Iso

      WRITE(stdout,'(A,ES17.4,I4,4ES17.4)') 'Iso table LoRes maximum SMAPE_pt, Moment, E, rho, T, Ye', & 
                  SMAPE_pt_max, t_m, max_rel_e_E, max_rel_e_rho, max_rel_e_T, max_rel_e_Ye
      WRITE(stdout,'(A,ES17.4,I4)') 'Iso table LoRes SMAPE, Moment', SMAPE, t_m

    END DO

    ErrorTable % Scat_Iso % Kernel(2) % Values( :, :, :, :, : ) = 0.0d0

  END IF
  
  FLUSH(stdout)

  IF( TableFlags(3) .eq. 1 )THEN
  
    WRITE(stdout,*) ' Checking Scat_NES ... '

    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % Scat_NES % Kernel(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % Scat_NES % Kernel(1) % Values)

    DO t_m = 1, nMom_NES 
      SMAPE_pt_max = 0.0d0
      SMAPE        = 0.0d0
      DO i_r = 1, nOpac_NES
        DO kp = 1, OpacityTableHi % nPointsE
          DO k = 1, OpacityTableHi % nPointsE
            DO i_eta = 1, OpacityTableHi % nPointsEta - 1
              DO k_t = 1, OpacityTableHi % nPointsTS(iT) - 1

                eta = OpacityTableHi % EtaGrid % Values(i_eta)
                T   = OpacityTableHi % TS % States (iT) % Values (k_t)  

                ReferenceValue = 10.d0**(OpacityTableHi % Scat_NES % Kernel(i_r) % Values &
                                 (kp,k,t_m,k_t,i_eta)) - OpacityTableHi % Scat_NES % Offsets(i_r,t_m)

                CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT, dT )
                CALL GetIndexAndDelta( LOG10(eta), LOG10(OpacityTableLo % EtaGrid % Values), idxEta, dEta )

                InterpValue = LinearInterp_Array_Point( kp, k, idxT, idxEta, dT, dEta, &
                              OpacityTableLo % Scat_NES % Offsets(i_r,t_m),       &
                              OpacityTableLo % Scat_NES % Kernel(i_r) % Values(:,:,t_m,:,:))

                IF((ABS(ReferenceValue)+ABS(InterpValue)) == 0.0d0) THEN
                  ErrorTable % Scat_NES % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta ) = 0.0d0
                ELSE
                  ErrorTable % Scat_NES % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta )               &
                    = 2.0d0 *ABS(ReferenceValue - InterpValue)/(ABS(ReferenceValue)+ABS(InterpValue))
                END IF

                IF(ErrorTable % Scat_NES % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta ) .gt. SMAPE_pt_max) THEN
                  SMAPE_pt_max  = ErrorTable % Scat_NES % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta )
                  max_rel_e_Ep  = OpacityTableHi % EnergyGrid % Values(kp)
                  max_rel_e_E   = OpacityTableHi % EnergyGrid % Values(k)
                  max_rel_e_eta = eta
                  max_rel_e_T   = T
                END IF
      
                SMAPE = SMAPE + ErrorTable % Scat_NES % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta )

              END DO
            END DO
          END DO
        END DO
      END DO

      SMAPE = SMAPE / SIZE(ErrorTable % Scat_NES % Kernel(1) % Values) * nMom_NES

      WRITE(stdout,'(A,ES17.4,I4,4ES17.4)') 'NES table LoRes maximum SMAPE_pt, Moment, Ep, E, eta, T', & 
                  SMAPE_pt_max, t_m, max_rel_e_Ep, max_rel_e_E, max_rel_e_eta, max_rel_e_T
      WRITE(stdout,'(A,ES17.4,I4)') 'NES table LoRes SMAPE, Moment', SMAPE, t_m

    END DO

  END IF

  FLUSH(stdout)

  IF( TableFlags(4) .eq. 1 )THEN

    WRITE(stdout,*) ' Checking Scat_Pair, zeroth moments ... '

    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % Scat_Pair % Kernel(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % Scat_Pair % Kernel(1) % Values)

    DO t_m = 1, 2 !nMom_Pair first order moments not stored in tables yet
      SMAPE_pt_max = 0.0d0
      SMAPE        = 0.0d0
      DO i_r = 1, nOpac_Pair
        DO kp = 1, OpacityTableHi % nPointsE
          DO k = 1, OpacityTableHi % nPointsE
            DO i_eta = 1, OpacityTableHi % nPointsEta - 1
              DO k_t = 1, OpacityTableHi % nPointsTS(iT) - 1

                eta = OpacityTableHi % EtaGrid % Values(i_eta)
                T   = OpacityTableHi % TS % States (iT) % Values (k_t)  

                ReferenceValue = 10.d0**(OpacityTableHi % Scat_Pair % Kernel(i_r) % Values &
                                 (kp,k,t_m,k_t,i_eta)) - OpacityTableHi % Scat_Pair % Offsets(i_r,t_m)

                CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT, dT )
                CALL GetIndexAndDelta( LOG10(eta), LOG10(OpacityTableLo % EtaGrid % Values), idxEta, dEta )

                InterpValue = LinearInterp_Array_Point( kp, k, idxT, idxEta, dT, dEta, &
                              OpacityTableLo % Scat_Pair % Offsets(i_r,t_m),       &
                              OpacityTableLo % Scat_Pair % Kernel(i_r) % Values(:,:,t_m,:,:))

                IF((ABS(ReferenceValue)+ABS(InterpValue)) == 0.0d0) THEN
                  ErrorTable % Scat_Pair % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta ) = 0.0d0
                ELSE
                  ErrorTable % Scat_Pair % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta )               &
                    = 2.0d0 *ABS(ReferenceValue - InterpValue)/(ABS(ReferenceValue)+ABS(InterpValue))
                END IF

                IF(ErrorTable % Scat_Pair % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta ) .gt. SMAPE_pt_max) THEN
                  SMAPE_pt_max  = ErrorTable % Scat_Pair % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta )
                  max_rel_e_Ep  = OpacityTableHi % EnergyGrid % Values(kp)
                  max_rel_e_E   = OpacityTableHi % EnergyGrid % Values(k)
                  max_rel_e_eta = eta
                  max_rel_e_T   = T
                END IF
      
                SMAPE = SMAPE + ErrorTable % Scat_Pair % Kernel(i_r) % Values( kp, k, t_m, k_t, i_eta )

              END DO
            END DO
          END DO
        END DO
      END DO

      SMAPE = SMAPE / SIZE(ErrorTable % Scat_Pair % Kernel(1) % Values) * nMom_Pair

      WRITE(stdout,'(A,ES17.4,I4,4ES17.4)') 'Pair table LoRes maximum SMAPE_pt, Moment, Ep, E, eta, T', & 
                  SMAPE_pt_max, t_m, max_rel_e_Ep, max_rel_e_E, max_rel_e_eta, max_rel_e_T
      WRITE(stdout,'(A,ES17.4,I4)') 'Pair table LoRes SMAPE, Moment', SMAPE, t_m

    END DO

    !
    ErrorTable % Scat_Pair % Kernel(1) % Values( :, :, 3, :, : ) = 0.0d0
    ErrorTable % Scat_Pair % Kernel(1) % Values( :, :, 4, :, : ) = 0.0d0

  END IF

  FLUSH(stdout)

  IF( TableFlags(5) .eq. 1 )THEN

    WRITE(stdout,*) ' Checking Scat_Brem ... '

    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % Scat_Brem % Kernel(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % Scat_Brem % Kernel(1) % Values)
  
    SMAPE_pt_max = 0.0d0
    SMAPE        = 0.0d0

    DO kp = 1, OpacityTableHi % nPointsE 
      DO k = 1, OpacityTableHi % nPointsE 
        DO k_t = 1, OpacityTableHi % nPointsTS(iT) - 1
          DO j_rho = 1, OpacityTableHi % nPointsTS(iRho) - 1 
            DO i_r = 1, nOpac_Brem
              DO t_m = 1, nMom_Brem
                rho = OpacityTableHi % TS % States (iRho) % Values (j_rho)
                T   = OpacityTableHi % TS % States (iT) % Values (k_t)

                ! table value in high-res. table
                ReferenceValue = 10.d0**(OpacityTableHi % Scat_Brem % Kernel(i_r) % Values &
                                 (kp,k,t_m,j_rho,k_t)) - OpacityTableHi % Scat_Brem % Offsets(i_r,t_m)

                CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT, dT )
                CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableLo % TS % States (iRho) % Values), idxRho, dRho )

                InterpValue = LinearInterp_Array_Point( kp, k, idxRho, idxT, dRho, dT, &
                              OpacityTableLo % Scat_Brem % Offsets(i_r,t_m),       &
                              OpacityTableLo % Scat_Brem % Kernel(i_r) % Values(:,:,t_m,:,:))

                IF((ABS(ReferenceValue)+ABS(InterpValue)) == 0.0d0) THEN
                  ErrorTable % Scat_Brem % Kernel(i_r) % Values( kp, k, t_m, j_rho, k_t ) = 0.0d0
                ELSE
                  ErrorTable % Scat_Brem % Kernel(i_r) % Values( kp, k, t_m, j_rho, k_t )               &
                    = 2.0d0 *ABS(ReferenceValue - InterpValue)/(ABS(ReferenceValue)+ABS(InterpValue))
                END IF


                IF(ErrorTable % Scat_Brem % Kernel(i_r) % Values( kp, k, t_m, j_rho, k_t ) .gt. SMAPE_pt_max) THEN
                  SMAPE_pt_max  = ErrorTable % Scat_Brem % Kernel(i_r) % Values( kp, k, t_m, j_rho, k_t )
                  max_rel_e_Ep  = OpacityTableHi % EnergyGrid % Values(kp)
                  max_rel_e_E   = OpacityTableHi % EnergyGrid % Values(k)
                  max_rel_e_rho = rho
                  max_rel_e_T   = T
                END IF
      
                SMAPE = SMAPE + ErrorTable % Scat_Brem % Kernel(i_r) % Values( kp, k, t_m, j_rho, k_t )

            END DO
          END DO
        END DO
      END DO
    END DO
  END DO

  SMAPE = SMAPE / SIZE(ErrorTable % Scat_Brem % Kernel(1) % Values) * nMom_Brem

  WRITE(stdout,'(A,5ES17.4)') 'Brem table LoRes maximum SMAPE_pt, Ep, E, rho, T', & 
                  SMAPE_pt_max, max_rel_e_Ep, max_rel_e_E, max_rel_e_rho, max_rel_e_T
  WRITE(stdout,'(A,ES17.4,I4)') 'Brem table LoRes SMAPE', SMAPE
 

  END IF

  END ASSOCIATE

  !CALL DescribeOpacityTable( ErrorTable )

  IF( nOpac_EmAb > 0 ) THEN
    WriteTableName = 'InterpolatedError_EmAb.h5'
    CALL InitializeHDF( )
    WRITE(stdout,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_EmAb_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_Iso > 0 ) THEN
    WriteTableName = 'InterpolatedError_Iso.h5'
    CALL InitializeHDF( )
    WRITE(stdout,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_Iso_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_NES > 0 ) THEN
    WriteTableName = 'InterpolatedError_NES.h5'
    CALL InitializeHDF( )
    WRITE(stdout,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_NES_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_Pair > 0 ) THEN
    WriteTableName = 'InterpolatedError_Pair.h5'
    CALL InitializeHDF( )
    WRITE(stdout,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_Pair_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_Brem > 0 ) THEN
    WriteTableName = 'InterpolatedError_Brem.h5'
    CALL InitializeHDF( )
    WRITE(stdout,*) 'Write Error into file = ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( ErrorTable, TRIM(WriteTableName), WriteOpacity_Brem_Option = .true. )
    CALL FinalizeHDF( )
  END IF

!  CALL DeAllocateOpacityTable( OpacityTableHi )
!  CALL DeAllocateOpacityTable( OpacityTableLo )
!  CALL DeAllocateOpacityTable( ErrorTable )
!   DEALLOCATE( InterEmAb, InterIso, TableValue1D )
!   DEALLOCATE( InterH0i, InterH0ii, InterJ0i, InterJ0ii, TableValue2D )

  FLUSH(stdout)

END PROGRAM wlOpacityTableResolutionTest

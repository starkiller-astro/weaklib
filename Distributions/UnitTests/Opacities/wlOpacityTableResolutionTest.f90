PROGRAM wlOpacityTableResolutionTest

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LinearInterp_Array_Point, &
    GetIndexAndDelta    
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
  USE B85_scattIso
  USE B85_scattNES
  USE B85_pair
  USE HR98_Bremsstrahlung
  USE prb_cntl_module, ONLY: &
      i_aeps, iaefnp, rhoaefnp, iaence, iaenct, roaenct, &
      edmpa, edmpe, iaenca
  USE, INTRINSIC :: iso_fortran_env, ONLY : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_nan

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
  TYPE(OpacityTableType) :: OpacityTableHi, OpacityTableLo
  REAL(dp), DIMENSION(2) :: Offset_NES

  !-------- variables for reading parameters data ---------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: Inte_rho, Inte_T, &
                                             Inte_Ye, Inte_cmpe, database
  CHARACTER(LEN=100)                      :: Format1, Format2, Format3
  !CHARACTER(LEN=30)                       :: a
  INTEGER                                 :: i, datasize, icmpe

  !-------- variables for output -------------------------------------------
  INTEGER(HID_T)                          :: file_id, group_id
  INTEGER(HSIZE_T)                        :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d
  INTEGER(HSIZE_T), DIMENSION(3)          :: datasize3d
  INTEGER(HSIZE_T), DIMENSION(4)          :: datasize4d
  
  CHARACTER(256)                          :: WriteTableName

  !-------- local variables -------------------------------------------------
  CHARACTER(128)                          :: FileNameHi(5), FileNameLo(5)

  INTEGER                                 :: ii, jj
  INTEGER                                 :: k, kp
  INTEGER                                 :: idxRho, idxT, idxYe, idxEta 
  REAL(dp)                                :: dRho, dT, dYe, dEta
  INTEGER                                 :: idxRho_Lo, idxT_Lo, idxYe_Lo, idxEta_Lo 
  REAL(dp)                                :: dRho_Lo, dT_Lo, dYe_Lo, dEta_Lo
  INTEGER                                 :: i_r, t_m

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

  REAL(dp)   :: InterpValue, InterpLoValue, InterpHiValue
  REAL(dp)   :: ReferenceValue

  INTEGER, PARAMETER :: n_rows = 297
  INTEGER, PARAMETER :: n_cols = 4

  REAL(dp), dimension(n_rows,n_cols) :: TS_profile_D15
  INTEGER :: n, n_r

  REAL(dp), DIMENSION(:,:,:,:),   ALLOCATABLE :: EmAb_Interp
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Scat_Iso_Interp
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Scat_NES_Interp
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Scat_Pair_Interp
  REAL(dp), DIMENSION(:,:,:,:),   ALLOCATABLE :: Scat_Brem_Interp

  REAL(dp), DIMENSION(:,:),       ALLOCATABLE :: EOS_quantities
  REAL(dp), DIMENSION(:),         ALLOCATABLE :: Scat_Eta

  !local variables for table building

  REAL(dp)                :: energy, TMeV, Z, A, &
                             chem_e, chem_n, chem_p, xheavy, xn, &
                             xp, xhe, bb, minvar 

  REAL(dp), DIMENSION(:),   ALLOCATABLE :: absor, emit
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: cok
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: H0i, H0ii, H1i, H1ii
  REAL(dp)                              :: j0i, j0ii, j1i, j1ii
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: s_a
  REAL(dp), PARAMETER                   :: brem_rho_min = 1.0d+07 !switch Bremsstrahlung off below rho_min
  REAL(dp), PARAMETER                   :: brem_rho_max = 1.0d+15 !switch Bremsstrahlung off above rho_max

  INTEGER :: n_errors

  n_errors = 0

  WRITE(stdout,*) 'Reading in radial rho, T, Ye profile at bounce'

  OPEN (UNIT=99, FILE='f_profile.d', STATUS='old', ACTION='read')

  DO n=1, n_rows
      READ(99,*) TS_profile_D15(n,1), TS_profile_D15(n,2), TS_profile_D15(n,3), TS_profile_D15(n,4)
  END DO 

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

  nPointsE   = OpacityTableHi % nPointsE
  nPointsEta = OpacityTableHi % nPointsEta

   ! set up OpacityTableTypeEmAb
   IF( nOpac_EmAb .gt. 0 ) THEN
     ALLOCATE(absor(nPointsE))
     ALLOCATE(emit(nPointsE))
     ALLOCATE(EmAb_Interp(nPointsE,n_rows,nOpac_EmAb,3))
     ALLOCATE(EOS_quantities(n_rows,9))
   END IF

   ! set up OpacityTableTypeScat Iso
   IF( nOpac_Iso .gt. 0 ) THEN
     ALLOCATE(cok(nPointsE,2))
     ALLOCATE(Scat_Iso_Interp(nPointsE,n_rows,nOpac_Iso,nMom_Iso,3))
     IF(.not. ALLOCATED(EOS_quantities)) ALLOCATE(EOS_quantities(n_rows,9))
   END IF

   ! set up OpacityTableTypeScat NES
   IF( nOpac_NES .gt. 0 ) THEN
     ALLOCATE(H0i(nPointsE,nPointsE))
     ALLOCATE(H0ii(nPointsE,nPointsE))
     ALLOCATE(H1i(nPointsE,nPointsE))
     ALLOCATE(H1ii(nPointsE,nPointsE))
     ALLOCATE(Scat_NES_Interp(nPointsE,nPointsE,n_rows,nMom_NES,3))
     ALLOCATE(Scat_Eta(n_rows))
   END IF

   ! set up OpacityTableTypeScat Pair
   IF( nOpac_Pair .gt. 0 ) THEN
     ALLOCATE(Scat_Pair_Interp(nPointsE,nPointsE,n_rows,nMom_Pair,3))
     IF(.not. ALLOCATED(Scat_Eta)) ALLOCATE(Scat_Eta(n_rows))
   END IF

   ! set up OpacityTableTypeBrem 
   IF( nOpac_Brem .gt. 0 ) THEN
     ALLOCATE(s_a(nPointsE,nPointsE))
     ALLOCATE(Scat_Brem_Interp(nPointsE,nPointsE,n_rows,3))
   END IF


  WRITE(stdout,*) 'Checking table resolution: Calculate reference value opacities from direct call to'
  WRITE(stdout,*) 'table building routines for a radial (rho,T,Ye) profile and then interpolate the'
  WRITE(stdout,*) 'opacities from high and low resolution tables'

  ASSOCIATE &
  ( iEOS_Rho  => OpacityTableHi % EOSTable % TS % Indices % iRho, &
    iEOS_T    => OpacityTableHi % EOSTable % TS % Indices % iT,   &
    iEOS_Ye   => OpacityTableHi % EOSTable % TS % Indices % iYe,  &
    Indices   => OpacityTableHi % EOSTable % DV % Indices,        &
    DVOffs    => OpacityTableHi % EOSTable % DV % Offsets,        &
    DVar      => OpacityTableHi % EOSTable % DV % Variables,      &
    iRho      => OpacityTableHi % TS % Indices % iRho,            &
    iT        => OpacityTableHi % TS % Indices % iT,              &
    iYe       => OpacityTableHi % TS % Indices % iYe,             &
    LogInterp => OpacityTableHi % EOSTable % TS % LogInterp )



  WRITE(stdout,'(A,2ES17.4)') 'HiRes min/max Ye', minval(OpacityTableHi % TS % States (iYe) % Values), &
                                                  maxval(OpacityTableHi % TS % States (iYe) % Values) 

  IF( TableFlags(1) .eq. 1 )THEN
  
    WRITE(stdout,*) ' Checking EmAb... '

    WRITE(stdout,'(A,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % EmAb % Opacity(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % EmAb % Opacity(1) % Values)

    DO n_r = 1, n_rows

      rho = TS_profile_D15(n_r,2)
      T   = TS_profile_D15(n_r,3)
      Ye  = TS_profile_D15(n_r,4)

      CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableHi % TS % States (iRho) % Values), idxRho, dRho )
      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableHi % TS % States (iT) % Values), idxT, dT )
      CALL GetIndexAndDelta( Ye,         OpacityTableHi % TS % States (iYe) % Values, idxYe, dYe )

      CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableLo % TS % States (iRho) % Values), idxRho_Lo, dRho_Lo )
      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT_Lo, dT_Lo )
      CALL GetIndexAndDelta( Ye,         OpacityTableLo % TS % States (iYe) % Values, idxYe_Lo, dYe_Lo )
               

      !interpolate derived quantities from EOS
      chem_e = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,           &
                                                DVOffs(Indices % iElectronChemicalPotential), &
                                                DVar(Indices % iElectronChemicalPotential) % Values)
      EOS_quantities(n_r,1) = chem_e

      chem_p = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,         &
                                                DVOffs(Indices % iProtonChemicalPotential), &
                                                DVar(Indices % iProtonChemicalPotential) % Values)
      EOS_quantities(n_r,2) = chem_p


      chem_n = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,          &
                                                DVOffs(Indices % iNeutronChemicalPotential), &
                                                DVar(Indices % iNeutronChemicalPotential) % Values)
      EOS_quantities(n_r,3) = chem_n


      xp = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,    &
                                            DVOffs(Indices % iProtonMassFraction), &
                                            DVar(Indices % iProtonMassFraction) % Values)
      EOS_quantities(n_r,4) = xp


      xn = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,     &
                                            DVOffs(Indices % iNeutronMassFraction), &
                                            DVar(Indices % iNeutronMassFraction) % Values)
      EOS_quantities(n_r,5) = xn


      xhe = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,   &
                                             DVOffs(Indices % iAlphaMassFraction), &
                                             DVar(Indices % iAlphaMassFraction) % Values)
      EOS_quantities(n_r,6) = xhe


      xheavy = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,   &
                                                DVOffs(Indices % iHeavyMassFraction), &
                                                DVar(Indices % iHeavyMassFraction) % Values)
      EOS_quantities(n_r,7) = xheavy


      Z = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,   &
                                           DVOffs(Indices % iHeavyChargeNumber), &
                                           DVar(Indices % iHeavyChargeNumber) % Values)
      EOS_quantities(n_r,8) = Z


      A = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe, &
                                           DVOffs(Indices % iHeavyMassNumber), &
                                           DVar(Indices % iHeavyMassNumber) % Values)
      EOS_quantities(n_r,9) = A


      bb  = (chem_e + chem_p - chem_n)/(T*kMev)

      iaefnp = 1
      i_aeps = 0
      rhoaefnp = HUGE(1.d0) ! (?) 
      iaence = 1
      edmpe = 3.d0
      iaenca = 1 
      edmpa = 3.d0
      iaenct = 0
      roaenct = TINY(1.d0)

      DO i_r = 1, nOpac_EmAb

        CALL abemrgn_weaklib &
             ( i_r, OpacityTableHi % EnergyGrid % Values, &
               rho, T, xn, xp, xheavy, &
               A, Z, chem_n, chem_p, chem_e, & 
               absor, emit, ye, nPointsE )
                

        DO k = 1, OpacityTableHi % nPointsE

          ReferenceValue = absor(k) + emit(k)
 
          InterpHiValue = LinearInterp_Array_Point( k, idxRho, idxT, idxYe, dRho, dT, dYe, &
                                                    OpacityTableHi % EmAb % Offsets(i_r),  &
                                                    OpacityTableHi % EmAb % Opacity(i_r) % Values(:,:,:,:))

          InterpLoValue = LinearInterp_Array_Point( k, idxRho_Lo, idxT_Lo, idxYe_Lo, dRho_Lo, dT_Lo, dYe_Lo, &
                                                    OpacityTableLo % EmAb % Offsets(i_r),  &
                                                    OpacityTableLo % EmAb % Opacity(i_r) % Values(:,:,:,:))
                

          EmAb_Interp(k,n_r,i_r,1) = ReferenceValue
          EmAb_Interp(k,n_r,i_r,2) = InterpLoValue
          EmAb_Interp(k,n_r,i_r,3) = InterpHiValue

          IF(ieee_is_nan(ReferenceValue) .or. ieee_is_nan(InterpLoValue) .or. ieee_is_nan(InterpHiValue)) &
            n_errors = n_errors + 1 

        END DO
      END DO
    END DO

  CALL InitializeHDF( )
  CALL OpenFileHDF( 'EmAb_table_resolution.h5', .true., file_id )

  CALL OpenGroupHDF( 'EmAb_Interp', .true., file_id, group_id )

  datasize1d(1) = nPointsE
  CALL WriteHDF( "Energy", OpacityTableHi % EnergyGrid % Values(:), group_id, datasize1d )
  datasize1d(1) = n_rows
  CALL WriteHDF( "Radius", TS_profile_D15(:,1), group_id, datasize1d )
  CALL WriteHDF( "rho",    TS_profile_D15(:,2), group_id, datasize1d )
  CALL WriteHDF( "T",      TS_profile_D15(:,3), group_id, datasize1d )
  CALL WriteHDF( "Ye",     TS_profile_D15(:,4), group_id, datasize1d )

  CALL WriteHDF( "ElectronChemicalPotential", EOS_quantities(:,1), group_id, datasize1d )
  CALL WriteHDF( "ProtonChemicalPotential",   EOS_quantities(:,2), group_id, datasize1d )
  CALL WriteHDF( "NeutronChemicalPotential",  EOS_quantities(:,3), group_id, datasize1d )
  CALL WriteHDF( "ProtonMassFraction",        EOS_quantities(:,4), group_id, datasize1d )
  CALL WriteHDF( "NeutronMassFraction",       EOS_quantities(:,5), group_id, datasize1d )
  CALL WriteHDF( "AlphaMassFraction",         EOS_quantities(:,6), group_id, datasize1d )
  CALL WriteHDF( "HeavyMassFraction",         EOS_quantities(:,7), group_id, datasize1d )
  CALL WriteHDF( "HeavyChargeNumber",         EOS_quantities(:,8), group_id, datasize1d )
  CALL WriteHDF( "HeavyMassNumber",           EOS_quantities(:,9), group_id, datasize1d )

  datasize2d = [nPointsE,n_rows]

  CALL WriteHDF("EmAb_Ref_nue",         EmAb_Interp(:,:,1,1), group_id, datasize2d)
  CALL WriteHDF("EmAb_LoInt_nue",       EmAb_Interp(:,:,1,2), group_id, datasize2d)
  CALL WriteHDF("EmAb_HiInt_nue",       EmAb_Interp(:,:,1,3), group_id, datasize2d)

  CALL WriteHDF("EmAb_Ref_nuebar",         EmAb_Interp(:,:,2,1), group_id, datasize2d)
  CALL WriteHDF("EmAb_LoInt_nuebar",       EmAb_Interp(:,:,2,2), group_id, datasize2d)
  CALL WriteHDF("EmAb_HiInt_nuebar",       EmAb_Interp(:,:,2,3), group_id, datasize2d)

  CALL CloseGroupHDF( group_id )   

  END IF

  IF( TableFlags(2) .eq. 1 )THEN
  
    WRITE(stdout,*) ' Checking Scat_Iso... '

    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % Scat_Iso % Kernel(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % Scat_Iso % Kernel(1) % Values)

    DO n_r = 1, n_rows

      rho = TS_profile_D15(n_r,2)
      T   = TS_profile_D15(n_r,3)
      Ye  = TS_profile_D15(n_r,4)

      CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableHi % TS % States (iRho) % Values), idxRho, dRho )
      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableHi % TS % States (iT) % Values), idxT, dT )
      CALL GetIndexAndDelta( Ye,         OpacityTableHi % TS % States (iYe) % Values, idxYe, dYe )

      CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableLo % TS % States (iRho) % Values), idxRho_Lo, dRho_Lo )
      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT_Lo, dT_Lo )
      CALL GetIndexAndDelta( Ye,         OpacityTableLo % TS % States (iYe) % Values, idxYe_Lo, dYe_Lo )
               

      !interpolate derived quantities from EOS
      chem_e = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,           &
                                                DVOffs(Indices % iElectronChemicalPotential), &
                                                DVar(Indices % iElectronChemicalPotential) % Values)
      EOS_quantities(n_r,1) = chem_e

      chem_p = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,         &
                                                DVOffs(Indices % iProtonChemicalPotential), &
                                                DVar(Indices % iProtonChemicalPotential) % Values)
      EOS_quantities(n_r,2) = chem_p


      chem_n = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,          &
                                                DVOffs(Indices % iNeutronChemicalPotential), &
                                                DVar(Indices % iNeutronChemicalPotential) % Values)
      EOS_quantities(n_r,3) = chem_n


      xp = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,    &
                                            DVOffs(Indices % iProtonMassFraction), &
                                            DVar(Indices % iProtonMassFraction) % Values)
      EOS_quantities(n_r,4) = xp


      xn = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,     &
                                            DVOffs(Indices % iNeutronMassFraction), &
                                            DVar(Indices % iNeutronMassFraction) % Values)
      EOS_quantities(n_r,5) = xn


      xhe = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,   &
                                             DVOffs(Indices % iAlphaMassFraction), &
                                             DVar(Indices % iAlphaMassFraction) % Values)
      EOS_quantities(n_r,6) = xhe


      xheavy = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,   &
                                                DVOffs(Indices % iHeavyMassFraction), &
                                                DVar(Indices % iHeavyMassFraction) % Values)
      EOS_quantities(n_r,7) = xheavy


      Z = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,   &
                                           DVOffs(Indices % iHeavyChargeNumber), &
                                           DVar(Indices % iHeavyChargeNumber) % Values)
      EOS_quantities(n_r,8) = Z


      A = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe, &
                                           DVOffs(Indices % iHeavyMassNumber), &
                                           DVar(Indices % iHeavyMassNumber) % Values)
      EOS_quantities(n_r,9) = A


      bb  = (chem_e + chem_p - chem_n)/(T*kMev)

      iaefnp = 1
      i_aeps = 0
      rhoaefnp = HUGE(1.d0) ! (?) 
      iaence = 1
      edmpe = 3.d0
      iaenca = 1 
      edmpa = 3.d0
      iaenct = 0
      roaenct = TINY(1.d0)

      DO i_r = 1, nOpac_Iso

        CALL scatical_weaklib &
           ( i_r, OpacityTableHi % EnergyGrid % Values, &
             nPointsE, rho, T, xn, xp, xhe, xheavy, A, Z, cok )
                
        DO t_m = 1, nMom_Iso
          DO k = 1, OpacityTableHi % nPointsE

            ReferenceValue = cok(k,t_m)
 
            InterpHiValue = LinearInterp_Array_Point( k, idxRho, idxT, idxYe, dRho, dT, dYe, &
                                                      OpacityTableHi % Scat_Iso % Offsets(i_r,t_m),  &
                                                      OpacityTableHi % Scat_Iso % Kernel(i_r) % Values(:,t_m,:,:,:))

            InterpLoValue = LinearInterp_Array_Point( k, idxRho_Lo, idxT_Lo, idxYe_Lo, dRho_Lo, dT_Lo, dYe_Lo, &
                                                      OpacityTableLo % Scat_Iso % Offsets(i_r,t_m),  &
                                                      OpacityTableLo % Scat_Iso % Kernel(i_r) % Values(:,t_m,:,:,:))
                

            Scat_Iso_Interp(k,n_r,i_r,t_m,1) = ReferenceValue
            Scat_Iso_Interp(k,n_r,i_r,t_m,2) = InterpLoValue
            Scat_Iso_Interp(k,n_r,i_r,t_m,3) = InterpHiValue

            IF(ieee_is_nan(ReferenceValue) .or. ieee_is_nan(InterpLoValue) .or. ieee_is_nan(InterpHiValue)) &
              n_errors = n_errors + 1 
          END DO
        END DO
      END DO
    END DO

  CALL InitializeHDF( )
  CALL OpenFileHDF( 'Iso_table_resolution.h5', .true., file_id )

  CALL OpenGroupHDF( 'Iso_Interp', .true., file_id, group_id )

  datasize1d(1) = nPointsE
  CALL WriteHDF( "Energy", OpacityTableHi % EnergyGrid % Values(:), group_id, datasize1d )
  datasize1d(1) = n_rows
  CALL WriteHDF( "Radius", TS_profile_D15(:,1), group_id, datasize1d )
  CALL WriteHDF( "rho",    TS_profile_D15(:,2), group_id, datasize1d )
  CALL WriteHDF( "T",      TS_profile_D15(:,3), group_id, datasize1d )
  CALL WriteHDF( "Ye",     TS_profile_D15(:,4), group_id, datasize1d )

  CALL WriteHDF( "ElectronChemicalPotential", EOS_quantities(:,1), group_id, datasize1d )
  CALL WriteHDF( "ProtonChemicalPotential",   EOS_quantities(:,2), group_id, datasize1d )
  CALL WriteHDF( "NeutronChemicalPotential",  EOS_quantities(:,3), group_id, datasize1d )
  CALL WriteHDF( "ProtonMassFraction",        EOS_quantities(:,4), group_id, datasize1d )
  CALL WriteHDF( "NeutronMassFraction",       EOS_quantities(:,5), group_id, datasize1d )
  CALL WriteHDF( "AlphaMassFraction",         EOS_quantities(:,6), group_id, datasize1d )
  CALL WriteHDF( "HeavyMassFraction",         EOS_quantities(:,7), group_id, datasize1d )
  CALL WriteHDF( "HeavyChargeNumber",         EOS_quantities(:,8), group_id, datasize1d )
  CALL WriteHDF( "HeavyMassNumber",           EOS_quantities(:,9), group_id, datasize1d )

  datasize3d = [nPointsE,n_rows,nMom_Iso]

  CALL WriteHDF("Iso_Ref_nue_Mom0",      Scat_Iso_Interp(:,:,1,1,1), group_id, datasize3d)
  CALL WriteHDF("Iso_LoInt_nue_Mom0",    Scat_Iso_Interp(:,:,1,1,2), group_id, datasize3d)
  CALL WriteHDF("Iso_HiInt_nue_Mom0",    Scat_Iso_Interp(:,:,1,1,3), group_id, datasize3d)

  CALL WriteHDF("Iso_Ref_nuebar_Mom0",   Scat_Iso_Interp(:,:,2,1,1), group_id, datasize3d)
  CALL WriteHDF("Iso_LoInt_nuebar_Mom0", Scat_Iso_Interp(:,:,2,1,2), group_id, datasize3d)
  CALL WriteHDF("Iso_HiInt_nuebar_Mom0", Scat_Iso_Interp(:,:,2,1,3), group_id, datasize3d)

  CALL WriteHDF("Iso_Ref_nue_Mom1",      Scat_Iso_Interp(:,:,1,2,1), group_id, datasize3d)
  CALL WriteHDF("Iso_LoInt_nue_Mom1",    Scat_Iso_Interp(:,:,1,2,2), group_id, datasize3d)
  CALL WriteHDF("Iso_HiInt_nue_Mom1",    Scat_Iso_Interp(:,:,1,2,3), group_id, datasize3d)

  CALL WriteHDF("Iso_Ref_nuebar_Mom1",   Scat_Iso_Interp(:,:,2,2,1), group_id, datasize3d)
  CALL WriteHDF("Iso_LoInt_nuebar_Mom1", Scat_Iso_Interp(:,:,2,2,2), group_id, datasize3d)
  CALL WriteHDF("Iso_HiInt_nuebar_Mom1", Scat_Iso_Interp(:,:,2,2,3), group_id, datasize3d)

  CALL CloseGroupHDF( group_id ) 

  END IF

  IF( TableFlags(3) .eq. 1 )THEN
  
    WRITE(stdout,*) ' Checking Scat_NES ... '

    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % Scat_NES % Kernel(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % Scat_NES % Kernel(1) % Values)

    DO n_r = 1, n_rows

      rho = TS_profile_D15(n_r,2)
      T   = TS_profile_D15(n_r,3)
      Ye  = TS_profile_D15(n_r,4)

      CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableHi % TS % States (iRho) % Values), idxRho, dRho )
      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableHi % TS % States (iT) % Values), idxT, dT )
      CALL GetIndexAndDelta( Ye,         OpacityTableHi % TS % States (iYe) % Values, idxYe, dYe )

      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT_Lo, dT_Lo )
               

      !interpolate chem_e from EOS to calculate eta
      chem_e = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,           &
                                                DVOffs(Indices % iElectronChemicalPotential), &
                                                DVar(Indices % iElectronChemicalPotential) % Values)

      TMeV = T * kMeV
      eta = chem_e / TMeV
      Scat_Eta(n_r) = eta

      CALL GetIndexAndDelta( LOG10(eta), LOG10(OpacityTableHi % EtaGrid % Values), idxEta, dEta )
      CALL GetIndexAndDelta( LOG10(eta), LOG10(OpacityTableLo % EtaGrid % Values), idxEta_Lo, dEta_Lo )

      CALL scatergn_weaklib &
               ( nPointsE, OpacityTableHi % EnergyGrid % Values, &
                 TMeV, eta, H0i, H0ii, H1i, H1ii )

      DO kp = 1, nPointsE
        DO k = 1, nPointsE

          Scat_NES_Interp(kp,k,n_r,1,1) = H0i (k,kp)  
          Scat_NES_Interp(kp,k,n_r,2,1) = H0ii(k,kp)  
          Scat_NES_Interp(kp,k,n_r,3,1) = H1i (k,kp)  
          Scat_NES_Interp(kp,k,n_r,4,1) = H1ii(k,kp)  

          DO t_m = 1, nMom_NES
            InterpLoValue = LinearInterp_Array_Point( kp, k, idxT_Lo, idxEta_Lo, dT_Lo, dEta_Lo, &
                                                      OpacityTableLo % Scat_NES % Offsets(1,t_m),       &
                                                      OpacityTableLo % Scat_NES % Kernel(1) % Values(:,:,t_m,:,:))
            InterpHiValue = LinearInterp_Array_Point( kp, k, idxT, idxEta, dT, dEta, &
                                                      OpacityTableHi % Scat_NES % Offsets(1,t_m),       &
                                                      OpacityTableHi % Scat_NES % Kernel(1) % Values(:,:,t_m,:,:))

            Scat_NES_Interp(kp,k,n_r,t_m,2) = InterpLoValue
            Scat_NES_Interp(kp,k,n_r,t_m,3) = InterpHiValue

            IF(ieee_is_nan(Scat_NES_Interp(kp,k,n_r,t_m,1)) .or. ieee_is_nan(InterpLoValue) .or. ieee_is_nan(InterpHiValue)) &
              n_errors = n_errors + 1 
          END DO
        END DO
      END DO

    END DO

  CALL InitializeHDF( )
  CALL OpenFileHDF( 'NES_table_resolution.h5', .true., file_id )

  CALL OpenGroupHDF( 'NES_Interp', .true., file_id, group_id )

  datasize1d(1) = nPointsE
  CALL WriteHDF( "Energy", OpacityTableHi % EnergyGrid % Values(:), group_id, datasize1d )
  datasize1d(1) = n_rows
  CALL WriteHDF( "Radius", TS_profile_D15(:,1), group_id, datasize1d )
  CALL WriteHDF( "rho",    TS_profile_D15(:,2), group_id, datasize1d )
  CALL WriteHDF( "T",      TS_profile_D15(:,3), group_id, datasize1d )
  CALL WriteHDF( "Ye",     TS_profile_D15(:,4), group_id, datasize1d )
  CALL WriteHDF( "Eta",    Scat_Eta(:),         group_id, datasize1d )

  datasize3d = [nPointsE,nPointsE,n_rows]

  CALL WriteHDF("NES_H0i_Ref",    Scat_NES_Interp(:,:,:,1,1), group_id, datasize3d)
  CALL WriteHDF("NES_H0i_LoInt",  Scat_NES_Interp(:,:,:,1,2), group_id, datasize3d)
  CALL WriteHDF("NES_H0i_HiInt",  Scat_NES_Interp(:,:,:,1,3), group_id, datasize3d)

  CALL WriteHDF("NES_H0ii_Ref",   Scat_NES_Interp(:,:,:,2,1), group_id, datasize3d)
  CALL WriteHDF("NES_H0ii_LoInt", Scat_NES_Interp(:,:,:,2,2), group_id, datasize3d)
  CALL WriteHDF("NES_H0ii_HiInt", Scat_NES_Interp(:,:,:,2,3), group_id, datasize3d)

  CALL WriteHDF("NES_H1i_Ref",    Scat_NES_Interp(:,:,:,3,1), group_id, datasize3d)
  CALL WriteHDF("NES_H1i_LoInt",  Scat_NES_Interp(:,:,:,3,2), group_id, datasize3d)
  CALL WriteHDF("NES_H1i_HiInt",  Scat_NES_Interp(:,:,:,3,3), group_id, datasize3d)

  CALL WriteHDF("NES_H1ii_Ref",   Scat_NES_Interp(:,:,:,4,1), group_id, datasize3d)
  CALL WriteHDF("NES_H1ii_LoInt", Scat_NES_Interp(:,:,:,4,2), group_id, datasize3d)
  CALL WriteHDF("NES_H1ii_HiInt", Scat_NES_Interp(:,:,:,4,3), group_id, datasize3d)

  CALL CloseGroupHDF( group_id ) 

  END IF

  IF( TableFlags(4) .eq. 1 )THEN

    WRITE(stdout,*) ' Checking Scat_Pair ... '

    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % Scat_Pair % Kernel(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % Scat_Pair % Kernel(1) % Values)

    DO n_r = 1, n_rows

      rho = TS_profile_D15(n_r,2)
      T   = TS_profile_D15(n_r,3)
      Ye  = TS_profile_D15(n_r,4)

      CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableHi % TS % States (iRho) % Values), idxRho, dRho )
      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableHi % TS % States (iT) % Values), idxT, dT )
      CALL GetIndexAndDelta( Ye,         OpacityTableHi % TS % States (iYe) % Values, idxYe, dYe )

      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT_Lo, dT_Lo )
               

      !interpolate chem_e from EOS to calculate eta
      chem_e = LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,           &
                                                DVOffs(Indices % iElectronChemicalPotential), &
                                                DVar(Indices % iElectronChemicalPotential) % Values)

      TMeV = T * kMeV
      eta = chem_e / TMeV
      Scat_Eta(n_r) = eta

      CALL GetIndexAndDelta( LOG10(eta), LOG10(OpacityTableHi % EtaGrid % Values), idxEta, dEta )
      CALL GetIndexAndDelta( LOG10(eta), LOG10(OpacityTableLo % EtaGrid % Values), idxEta_Lo, dEta_Lo )

      DO kp = 1, nPointsE
        DO k = 1, nPointsE

          CALL paircal_weaklib( OpacityTableHi % EnergyGrid % Values(k),  &
                                OpacityTableHi % EnergyGrid % Values(kp), &
                                chem_e, T, j0i, j0ii, j1i, j1ii )

          Scat_Pair_Interp(kp,k,n_r,1,1) = j0i 
          Scat_Pair_Interp(kp,k,n_r,2,1) = j0ii  
          Scat_Pair_Interp(kp,k,n_r,3,1) = j1i   
          Scat_Pair_Interp(kp,k,n_r,4,1) = j1ii  

          DO t_m = 1, nMom_Pair
            InterpLoValue = LinearInterp_Array_Point( kp, k, idxT_Lo, idxEta_Lo, dT_Lo, dEta_Lo, &
                                                      OpacityTableLo % Scat_Pair % Offsets(1,t_m),       &
                                                      OpacityTableLo % Scat_Pair % Kernel(1) % Values(:,:,t_m,:,:))
            InterpHiValue = LinearInterp_Array_Point( kp, k, idxT, idxEta, dT, dEta, &
                                                      OpacityTableHi % Scat_Pair % Offsets(1,t_m),       &
                                                      OpacityTableHi % Scat_Pair % Kernel(1) % Values(:,:,t_m,:,:))

            Scat_Pair_Interp(kp,k,n_r,t_m,2) = InterpLoValue
            Scat_Pair_Interp(kp,k,n_r,t_m,3) = InterpHiValue

            IF(ieee_is_nan(Scat_Pair_Interp(kp,k,n_r,t_m,1)) .or. ieee_is_nan(InterpLoValue) .or. ieee_is_nan(InterpHiValue)) &
              n_errors = n_errors + 1 
          END DO
        END DO
      END DO

    END DO

  CALL InitializeHDF( )
  CALL OpenFileHDF( 'Pair_table_resolution.h5', .true., file_id )

  CALL OpenGroupHDF( 'Pair_Interp', .true., file_id, group_id )

  datasize1d(1) = nPointsE
  CALL WriteHDF( "Energy", OpacityTableHi % EnergyGrid % Values(:), group_id, datasize1d )
  datasize1d(1) = n_rows
  CALL WriteHDF( "Radius", TS_profile_D15(:,1), group_id, datasize1d )
  CALL WriteHDF( "rho",    TS_profile_D15(:,2), group_id, datasize1d )
  CALL WriteHDF( "T",      TS_profile_D15(:,3), group_id, datasize1d )
  CALL WriteHDF( "Ye",     TS_profile_D15(:,4), group_id, datasize1d )
  CALL WriteHDF( "Eta",    Scat_Eta(:),         group_id, datasize1d )

  datasize3d = [nPointsE,nPointsE,n_rows]

  CALL WriteHDF("Pair_j0i_Ref",    Scat_Pair_Interp(:,:,:,1,1), group_id, datasize3d)
  CALL WriteHDF("Pair_j0i_LoInt",  Scat_Pair_Interp(:,:,:,1,2), group_id, datasize3d)
  CALL WriteHDF("Pair_j0i_HiInt",  Scat_Pair_Interp(:,:,:,1,3), group_id, datasize3d)

  CALL WriteHDF("Pair_j0ii_Ref",   Scat_Pair_Interp(:,:,:,2,1), group_id, datasize3d)
  CALL WriteHDF("Pair_j0ii_LoInt", Scat_Pair_Interp(:,:,:,2,2), group_id, datasize3d)
  CALL WriteHDF("Pair_j0ii_HiInt", Scat_Pair_Interp(:,:,:,2,3), group_id, datasize3d)

  CALL WriteHDF("Pair_j1i_Ref",    Scat_Pair_Interp(:,:,:,3,1), group_id, datasize3d)
  CALL WriteHDF("Pair_j1i_LoInt",  Scat_Pair_Interp(:,:,:,3,2), group_id, datasize3d)
  CALL WriteHDF("Pair_j1i_HiInt",  Scat_Pair_Interp(:,:,:,3,3), group_id, datasize3d)

  CALL WriteHDF("Pair_j1ii_Ref",   Scat_Pair_Interp(:,:,:,4,1), group_id, datasize3d)
  CALL WriteHDF("Pair_j1ii_LoInt", Scat_Pair_Interp(:,:,:,4,2), group_id, datasize3d)
  CALL WriteHDF("Pair_j1ii_HiInt", Scat_Pair_Interp(:,:,:,4,3), group_id, datasize3d)

  CALL CloseGroupHDF( group_id ) 

  END IF

  IF( TableFlags(5) .eq. 1 )THEN

    WRITE(stdout,*) ' Checking Scat_Brem ... '

    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of LoResTable', shape(OpacityTableLo % Scat_Brem % Kernel(1) % Values)
    WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of HiResTable', shape(OpacityTableHi % Scat_Brem % Kernel(1) % Values)

    DO n_r = 1, n_rows

      rho = TS_profile_D15(n_r,2)
      T   = TS_profile_D15(n_r,3)
      Ye  = TS_profile_D15(n_r,4)

      CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableHi % TS % States (iRho) % Values), idxRho, dRho )
      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableHi % TS % States (iT) % Values), idxT, dT )

      CALL GetIndexAndDelta( LOG10(rho), LOG10(OpacityTableLo % TS % States (iRho) % Values), idxRho_Lo, dRho_Lo )
      CALL GetIndexAndDelta( LOG10(T),   LOG10(OpacityTableLo % TS % States (iT) % Values), idxT_Lo, dT_Lo )

      CALL bremcal_weaklib(nPointsE, OpacityTableHi % EnergyGrid % Values, & 
                           rho, T, s_a)

      DO kp = 1, nPointsE
        DO k = 1, nPointsE

          IF (rho < brem_rho_min .or. rho > brem_rho_max) THEN

            Scat_Brem_Interp(kp,k,n_r,1) = 0.d0
             
          ELSE
            Scat_Brem_Interp(kp,k,n_r,1) = s_a(k,kp) 
          END IF

          IF(ieee_is_nan(s_a(k,kp))) n_errors = n_errors + 1

          DO t_m = 1, nMom_Brem
            InterpLoValue = LinearInterp_Array_Point( kp, k, idxRho_Lo, idxT_Lo, dRho_Lo, dT_Lo, &
                                                      OpacityTableLo % Scat_Brem % Offsets(1,t_m),       &
                                                      OpacityTableLo % Scat_Brem % Kernel(1) % Values(:,:,t_m,:,:))
            InterpHiValue = LinearInterp_Array_Point( kp, k, idxRho, idxT, dRho, dT, &
                                                      OpacityTableHi % Scat_Brem % Offsets(1,t_m),       &
                                                      OpacityTableHi % Scat_Brem % Kernel(1) % Values(:,:,t_m,:,:))

            Scat_Brem_Interp(kp,k,n_r,2) = InterpLoValue
            Scat_Brem_Interp(kp,k,n_r,3) = InterpHiValue
            IF(ieee_is_nan(InterpLoValue) .or. ieee_is_nan(InterpHiValue)) &
              n_errors = n_errors + 1 

          END DO
        END DO
      END DO

    END DO

  CALL InitializeHDF( )
  CALL OpenFileHDF( 'Brem_table_resolution.h5', .true., file_id )

  CALL OpenGroupHDF( 'Brem_Interp', .true., file_id, group_id )

  datasize1d(1) = nPointsE
  CALL WriteHDF( "Energy", OpacityTableHi % EnergyGrid % Values(:), group_id, datasize1d )
  datasize1d(1) = n_rows
  CALL WriteHDF( "Radius", TS_profile_D15(:,1), group_id, datasize1d )
  CALL WriteHDF( "rho",    TS_profile_D15(:,2), group_id, datasize1d )
  CALL WriteHDF( "T",      TS_profile_D15(:,3), group_id, datasize1d )
  CALL WriteHDF( "Ye",     TS_profile_D15(:,4), group_id, datasize1d )

  datasize3d = [nPointsE,nPointsE,n_rows]

  CALL WriteHDF("Brem_s_a_Ref",    Scat_Brem_Interp(:,:,:,1), group_id, datasize3d)
  CALL WriteHDF("Brem_s_a_LoInt",  Scat_Brem_Interp(:,:,:,2), group_id, datasize3d)
  CALL WriteHDF("Brem_s_a_HiInt",  Scat_Brem_Interp(:,:,:,3), group_id, datasize3d)


  CALL CloseGroupHDF( group_id ) 

  END IF

  END ASSOCIATE

  IF(MAXVAL(TableFlags) > 0) THEN

    CALL DeAllocateOpacityTable( OpacityTableHi )
    CALL DeAllocateOpacityTable( OpacityTableLo )
  END IF

  IF(ALLOCATED(EmAb_Interp))      DEALLOCATE(EmAb_Interp)
  IF(ALLOCATED(Scat_Iso_Interp))  DEALLOCATE(Scat_Iso_Interp)
  IF(ALLOCATED(Scat_NES_Interp))  DEALLOCATE(Scat_NES_Interp)
  IF(ALLOCATED(Scat_Pair_Interp)) DEALLOCATE(Scat_Pair_Interp)
  IF(ALLOCATED(Scat_Brem_Interp)) DEALLOCATE(Scat_Brem_Interp)

  IF(ALLOCATED(EOS_quantities))   DEALLOCATE(EOS_quantities)
  IF(ALLOCATED(Scat_Eta))         DEALLOCATE(Scat_Eta)

  IF(ALLOCATED(absor))            DEALLOCATE(absor)
  IF(ALLOCATED(emit))             DEALLOCATE(emit)
  IF(ALLOCATED(cok))              DEALLOCATE(cok)
  IF(ALLOCATED(H0i))              DEALLOCATE(H0i)
  IF(ALLOCATED(H0ii))             DEALLOCATE(H0ii)
  IF(ALLOCATED(H1i))              DEALLOCATE(H1i)
  IF(ALLOCATED(H1ii))             DEALLOCATE(H1ii)
  IF(ALLOCATED(s_a))              DEALLOCATE(s_a)

  IF(n_errors > 0) THEN
    WRITE(stdout,*) 'FAILED table resolution test, nans were found in table building routines or interpolated values!'
  ELSE
    WRITE(stdout,*) 'PASSED table resolution test.'
  ENDIF

END PROGRAM wlOpacityTableResolutionTest

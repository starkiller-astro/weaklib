PROGRAM wlInterpolateBrem

  USE wlKindModule, ONLY: dp
  USE wlInterpolationUtilitiesModule, ONLY: &
    LinearInterp_Array_Point, &
    GetIndexAndDelta_Lin, &
    LinearInterp2D_4DArray_2DAligned_Point
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
  USE HDF5
  USE, INTRINSIC :: iso_fortran_env, ONLY : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_nan

  IMPLICIT NONE

  !-------- Table filename --------------------------------------------------

  CHARACTER(256) :: HighResOpTableBaseBrem = "wl-Op-SFHo-15-25-50-E40-HR98"

  CHARACTER(256) :: HighResEOSTableName    = "wl-EOS-SFHo-15-25-50.h5"

  !-------- variables for reading opacity table -----------------------------
  TYPE(OpacityTableType) :: OpacityTable

  !-------- variables for reading parameters data ---------------------------
  INTEGER                                 :: i, datasize, icmpe

  !-------- variables for output -------------------------------------------
  INTEGER(HID_T)                          :: file_id, group_id
  INTEGER(HSIZE_T)                        :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d
  INTEGER(HSIZE_T), DIMENSION(3)          :: datasize3d
  INTEGER(HSIZE_T), DIMENSION(4)          :: datasize4d
  
  !-------- local variables -------------------------------------------------
  CHARACTER(128)                          :: FileName

  INTEGER                                 :: ii, jj
  INTEGER                                 :: k, kp
  INTEGER                                 :: idxRho, idxT, idxYe, idxEta 
  REAL(dp)                                :: dRho, dT, dYe, dEta
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

  REAL(dp)                                :: xp, xn
  REAL(dp), PARAMETER                     :: coef_np  = 28.d0/3.d0

  REAL(dp)   :: InterpValue, InterpValue_xp, InterpValue_xn, InterpValue_xpxn

  INTEGER    :: idxRho_xp, idxRho_xn, idxRho_xpxn
  REAL(dp)   :: dRho_xp, dRho_xn, dRho_xpxn

  INTEGER, PARAMETER :: n_rows = 213
  INTEGER, PARAMETER :: n_cols = 4

  REAL(dp), dimension(n_rows,n_cols) :: TS_profile
  INTEGER :: n, n_r

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Scat_Brem_Interp
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Scat_Brem_Interp_Decomp

  REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: EOS_quantities


  REAL(dp), PARAMETER                   :: brem_rho_min = 1.0d+07 !switch Bremsstrahlung off below rho_min
  REAL(dp), PARAMETER                   :: brem_rho_max = 1.0d+15 !switch Bremsstrahlung off above rho_max

  INTEGER :: n_errors

  n_errors = 0

  WRITE(stdout,*) 'Reading in radial rho, T, Ye profile'

  OPEN (UNIT=99, FILE='profile.d', STATUS='old', ACTION='read')

  DO n=1, n_rows
      READ(99,*) TS_profile(n,1), TS_profile(n,2), TS_profile(n,3), TS_profile(n,4)
  END DO 


  FileName = TRIM(HighResOpTableBaseBrem)//'-Brem.h5'

  !------------------------------------------------------
  ! read in the reference high resolution opacity table
  !------------------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable,             &
       FileName_EmAB_Option                             &
       = '',                           &
       FileName_Iso_Option                              &
       = '',                           &
       FileName_NES_Option                              &
       = '',                           &
       FileName_Pair_Option                             &
       = '',                           &
       FileName_Brem_Option                             &
       = TRIM(FileName),                           &
       EquationOfStateTableName_Option                  &
       = TRIM(HighResEOSTableName),                     &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )
  

  nOpac_EmAb = 0

  nOpac_Iso = 0
  nMom_Iso  = 0

  nOpac_NES = 0
  nMom_NES  = 0

  nOpac_Pair = 0
  nMom_Pair  = 0

  nOpac_Brem = OpacityTable % nOpacities_Brem
  nMom_Brem  = OpacityTable % nMoments_Brem

  nPointsE   = OpacityTable % nPointsE
  nPointsEta = OpacityTable % nPointsEta


   ! set up OpacityTableTypeBrem 
  ALLOCATE(Scat_Brem_Interp(nPointsE,nPointsE,n_rows))
  ALLOCATE(Scat_Brem_Interp_Decomp(nPointsE,nPointsE,n_rows))


  ASSOCIATE &
  ( iEOS_Rho  => OpacityTable % EOSTable % TS % Indices % iRho, &
    iEOS_T    => OpacityTable % EOSTable % TS % Indices % iT,   &
    iEOS_Ye   => OpacityTable % EOSTable % TS % Indices % iYe,  &
    Indices   => OpacityTable % EOSTable % DV % Indices,        &
    DVOffs    => OpacityTable % EOSTable % DV % Offsets,        &
    DVar      => OpacityTable % EOSTable % DV % Variables,      &
    iRho      => OpacityTable % TS % Indices % iRho,            &
    iT        => OpacityTable % TS % Indices % iT,              &
    iYe       => OpacityTable % TS % Indices % iYe,             &
    LogInterp => OpacityTable % EOSTable % TS % LogInterp )



  WRITE(stdout,'(A,2ES17.4)') 'OpTab min/max Ye', minval(OpacityTable % TS % States (iYe) % Values), &
                                                  maxval(OpacityTable % TS % States (iYe) % Values) 




  WRITE(stdout,*) ' Interpolating Scat_Brem kernel... '

  WRITE(stdout,'(A,I4,I4,I4,I4,I4)') 'shape of Opacity Table', shape(OpacityTable % Scat_Brem % Kernel(1) % Values)

    DO n_r = 1, n_rows

      rho = TS_profile(n_r,2)
      T   = TS_profile(n_r,3)
      Ye  = TS_profile(n_r,4)


      CALL GetIndexAndDelta_Lin( LOG10(rho), LOG10(OpacityTable % TS % States (iRho) % Values), idxRho, dRho )
      CALL GetIndexAndDelta_Lin( LOG10(T),   LOG10(OpacityTable % TS % States (iT) % Values), idxT, dT )
      CALL GetIndexAndDelta_Lin( Ye,         OpacityTable % TS % States (iYe) % Values, idxYe, dYe )

      CALL LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,    &
                                            DVOffs(Indices % iProtonMassFraction), &
                                            DVar(Indices % iProtonMassFraction) % Values, &
                                            xp)

      CALL LinearInterp_Array_Point(idxRho, idxT, idxYe, dRho, dT, dYe,     &
                                            DVOffs(Indices % iNeutronMassFraction), &
                                            DVar(Indices % iNeutronMassFraction) % Values, &
                                            xn)

      IF(rho*xp > minval(OpacityTable % TS % States (iRho) % Values)) THEN
        CALL GetIndexAndDelta_Lin( LOG10(rho*xp), LOG10(OpacityTable % TS % States (iRho) % Values), idxRho_xp, dRho_xp )
      ELSE
        idxRho_xp = -1
        dRho_xp = 0.d0
      END IF

      IF(rho*xn > minval(OpacityTable % TS % States (iRho) % Values)) THEN
        CALL GetIndexAndDelta_Lin( LOG10(rho*xn), LOG10(OpacityTable % TS % States (iRho) % Values), idxRho_xp, dRho_xp )
      ELSE
        idxRho_xn = -1
        dRho_xn = 0.d0
      END IF

      IF(rho*SQRT( ABS( xn * xp ) + 1d-100 ) > minval(OpacityTable % TS % States (iRho) % Values)) THEN
        CALL GetIndexAndDelta_Lin( LOG10(rho*SQRT( ABS( xn * xp ) + 1d-100 )), & 
                               LOG10(OpacityTable % TS % States (iRho) % Values), idxRho_xpxn, dRho_xpxn )
      ELSE
        idxRho_xpxn = -1
        dRho_xpxn = 0.d0
      END IF

      DO kp = 1, nPointsE
        DO k = 1, nPointsE


          DO t_m = 1, nMom_Brem
            CALL LinearInterp2D_4DArray_2DAligned_Point( kp, k, idxRho, idxT, dRho, dT, &
                                                         OpacityTable % Scat_Brem % Offsets(1,t_m),       &
                                                         OpacityTable % Scat_Brem % Kernel(1) % Values(:,:,t_m,:,:), &
                                                         InterpValue )

            IF(idxRho_xp > 0) THEN
              CALL LinearInterp2D_4DArray_2DAligned_Point( kp, k, idxRho_xp, idxT, dRho_xp, dT, &
                                                           OpacityTable % Scat_Brem % Offsets(1,t_m),       &
                                                           OpacityTable % Scat_Brem % Kernel(1) % Values(:,:,t_m,:,:), &
                                                           InterpValue_xp)
            ELSE
              InterpValue_xp = 0.d0
            END IF

            IF(idxRho_xn > 0) THEN
              CALL LinearInterp2D_4DArray_2DAligned_Point( kp, k, idxRho_xn, idxT, dRho_xn, dT, &
                                                           OpacityTable % Scat_Brem % Offsets(1,t_m),       &
                                                           OpacityTable % Scat_Brem % Kernel(1) % Values(:,:,t_m,:,:), &
                                                           InterpValue_xn)
            ELSE
              InterpValue_xpxn = 0.d0
            END IF

            IF(idxRho_xpxn > 0) THEN
              CALL LinearInterp2D_4DArray_2DAligned_Point( kp, k, idxRho_xpxn, idxT, dRho_xpxn, dT, &
                                                           OpacityTable % Scat_Brem % Offsets(1,t_m),       &
                                                           OpacityTable % Scat_Brem % Kernel(1) % Values(:,:,t_m,:,:), &
                                                           InterpValue_xpxn)
            ELSE
              InterpValue_xpxn = 0.d0
            END IF


            Scat_Brem_Interp(kp,k,n_r) = InterpValue

            Scat_Brem_Interp_Decomp(kp,k,n_r) = InterpValue_xp + InterpValue_xn + coef_np * InterpValue_xpxn

            IF(ieee_is_nan(InterpValue))      n_errors = n_errors + 1 
            IF(ieee_is_nan(InterpValue_xp))   n_errors = n_errors + 1 
            IF(ieee_is_nan(InterpValue_xn))   n_errors = n_errors + 1 
            IF(ieee_is_nan(InterpValue_xpxn)) n_errors = n_errors + 1 

          END DO
        END DO
      END DO

    END DO

  CALL InitializeHDF( )
  CALL OpenFileHDF( 'InterpolatedBremOutput.h5', .true., file_id )

  CALL OpenGroupHDF( 'Brem_Interp', .true., file_id, group_id )

  datasize1d(1) = nPointsE
  CALL WriteHDF( "Energy", OpacityTable % EnergyGrid % Values(:), group_id, datasize1d )
  datasize1d(1) = n_rows
  CALL WriteHDF( "Radius", TS_profile(:,1), group_id, datasize1d )
  CALL WriteHDF( "rho",    TS_profile(:,2), group_id, datasize1d )
  CALL WriteHDF( "T",      TS_profile(:,3), group_id, datasize1d )
  CALL WriteHDF( "Ye",     TS_profile(:,4), group_id, datasize1d )

  datasize3d = [nPointsE,nPointsE,n_rows]

  CALL WriteHDF("Brem_s_a_interp",        Scat_Brem_Interp(:,:,:), group_id, datasize3d)
  CALL WriteHDF("Brem_s_a_interp_decomp", Scat_Brem_Interp_Decomp(:,:,:), group_id, datasize3d)


  CALL CloseGroupHDF( group_id ) 


  END ASSOCIATE

  CALL DeAllocateOpacityTable( OpacityTable )

  IF(ALLOCATED(Scat_Brem_Interp))        DEALLOCATE(Scat_Brem_Interp)
  IF(ALLOCATED(Scat_Brem_Interp_Decomp)) DEALLOCATE(Scat_Brem_Interp_Decomp)

  IF(ALLOCATED(EOS_quantities))   DEALLOCATE(EOS_quantities)

  IF(n_errors > 0) THEN
    WRITE(stdout,*) 'FAILED table resolution test, nans were found in table building routines or interpolated values!'
  ELSE
    WRITE(stdout,*) 'PASSED table resolution test.'
  ENDIF

END PROGRAM wlInterpolateBrem

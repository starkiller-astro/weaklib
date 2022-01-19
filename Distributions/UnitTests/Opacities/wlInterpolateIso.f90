PROGRAM wlInterpolateIso

  USE wlKindModule, ONLY: dp
  USE wlExtNumericalModule, ONLY: frpi
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_1D3D_Custom
  USE wlOpacityFieldsModule, ONLY: &
    iNu_e, iNu_e_bar
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DeAllocateOpacityTable
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF, &
    OpenFileHDF, &
    CloseFileHDF, &
    WriteHDF, &
    OpenGroupHDF, &
    CloseGroupHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DescribeGrid, &
    MakeLogGrid
  USE HDF5

  IMPLICIT NONE

  !--------- parameters for creating energy grid ----------------------------
  INTEGER, PARAMETER :: Inte_nPointE = 80
  REAL(dp)           :: Inte_Emin = 1.0d-1
  REAL(dp)           :: Inte_Emax = 3.0d02
  TYPE(GridType)     :: Inte_E

  !-------- variables for reading opacity table -----------------------------
  TYPE(OpacityTableType)   :: OpacityTable
  REAL(dp), DIMENSION(2,2) :: Offset_Iso

  !-------- variables for reading parameters data ---------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Inte_r, Inte_rho, Inte_T, &
                                         Inte_Ye, database
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3
  CHARACTER(LEN=30)                   :: a
  INTEGER                             :: i, datasize

  !-------- variables for interpolation -------------------------------------
  INTEGER                               :: iE
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: bufferIso10, bufferIso11, &
                                           bufferIso20, bufferIso21

  !-------- variables for output --------------------------------------------
  INTEGER(HID_T)                        :: file_id, group_id
  INTEGER(HSIZE_T)                      :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)        :: datasize2d

  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: InterpolantIso1
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: InterpolantIso2

  CHARACTER(LEN=128)                      :: EosTableName, OpTableName, &
                                             ProfileName

  CHARACTER(LEN=30)                     :: outfilename = &
                                           'InterpolatedIsoOutput.h5'

  !--------------------------------------------------
  !   take in profile and table names
  !--------------------------------------------------
  OPEN(12, FILE="dataList.txt", FORM = "formatted", &
           ACTION = 'read')
  READ(12,*) ProfileName
  READ(12,*) EosTableName
  READ(12,*) OpTableName ! AbEm's name, skipt
  READ(12,*) OpTableName
  CLOSE(12, STATUS = 'keep')

  !------------------------------------------------
  !   interpolated energy 
  !------------------------------------------------
  Format1 = "(A2,I5)"
  Format2 = "(4A12)"
  Format3 = "(4ES12.3)"

  CALL AllocateGrid( Inte_E, Inte_nPointE )

  Inte_E % Unit = 'MeV                  '
  Inte_E % Name = 'Intepolated Energy   '
  Inte_E % MinValue  = Inte_Emin
  Inte_E % MaxValue  = Inte_Emax
  Inte_E % LogInterp = 1
  Inte_E % nPoints   = Inte_nPointE

  CALL MakeLogGrid &
          ( Inte_E % MinValue, Inte_E % MaxValue, &
            Inte_E % nPoints, Inte_E % Values )

  !-------------------------------------------------------
  !    read in profile ( rho, T, Ye )
  !-------------------------------------------------------
  OPEN(1, FILE = TRIM(ProfileName), FORM = "formatted", &
          ACTION = 'read')
  READ( 1, Format1 ) a, datasize
  READ( 1, Format2 )

  ALLOCATE( database( datasize * 4) )
  ALLOCATE( Inte_r  ( datasize ) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_T  ( datasize ) )
  ALLOCATE( Inte_Ye ( datasize ) )
  ALLOCATE( bufferIso10    ( Inte_nPointE, datasize ) )
  ALLOCATE( bufferIso11    ( Inte_nPointE, datasize ) )
  ALLOCATE( bufferIso20    ( Inte_nPointE, datasize ) )
  ALLOCATE( bufferIso21    ( Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantIso1( Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantIso2( Inte_nPointE, datasize ) )

  READ( 1, Format3 ) database
  CLOSE( 1, STATUS = 'keep')  

  DO i = 1, datasize  
    Inte_r(i)   = database(i*4-3)
    Inte_rho(i) = database(i*4-2)
    Inte_T(i)   = database(i*4-1)
    Inte_Ye(i)  = database(i*4)
  END DO 

  !---------------------------------------------------------------
  !    read in the opacity table
  !---------------------------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, &
       FileName_Iso_Option                &
       = TRIM(OpTableName),               &
       EquationOfStateTableName_Option    &
       = TRIM(EosTableName),              &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )

  Offset_Iso(1:2,1:2) = OpacityTable % Scat_Iso % Offsets(1:2,1:2)

  !-------------------------------------------------------------------------
  !   do interpolation
  !-------------------------------------------------------------------------
  ASSOCIATE&
  ( TableES10  => OpacityTable % Scat_Iso % Kernel(iNu_e) % Values(:,1,:,:,:),    &
    TableES11  => OpacityTable % Scat_Iso % Kernel(iNu_e) % Values(:,2,:,:,:),    &
    TableES20  => OpacityTable % Scat_Iso % Kernel(iNu_e_bar) % Values(:,1,:,:,:),&
    TableES21  => OpacityTable % Scat_Iso % Kernel(iNu_e_bar) % Values(:,2,:,:,:),&
    Energy     => Inte_E % Values,                    &
    iRho       => OpacityTable % TS % Indices % iRho, &
    iT         => OpacityTable % TS % Indices % iT,   &
    iYe        => OpacityTable % TS % Indices % iYe )
  
  CALL LogInterpolateSingleVariable_1D3D_Custom         &
         ( LOG10( Energy ), LOG10( Inte_rho ),          &
           LOG10( Inte_T ), Inte_Ye,                    &
           LOG10( OpacityTable % EnergyGrid % Values ), &
           LOG10( OpacityTable % TS % States(iRho) % Values ), &
           LOG10( OpacityTable % TS % States(iT) % Values ),   &
           OpacityTable % TS % States(iYe) % Values,           &
           Offset_Iso(iNu_e,1), TableES10, bufferIso10 )

  CALL LogInterpolateSingleVariable_1D3D_Custom         &
         ( LOG10( Energy ), LOG10( Inte_rho ),          &
           LOG10( Inte_T ), Inte_Ye,                    &
           LOG10( OpacityTable % EnergyGrid % Values ), &
           LOG10( OpacityTable % TS % States(iRho) % Values ), &
           LOG10( OpacityTable % TS % States(iT) % Values ),   &
           OpacityTable % TS % States(iYe) % Values,           &
           Offset_Iso(iNu_e,2), TableES11, bufferIso11 )

  DO iE = 1, Inte_nPointE
    InterpolantIso1(iE,:) = frpi * Energy(iE)**2 &
                            * ( bufferIso10(iE,:) - bufferIso11(iE,:) / 3.0d0 )
                     ! (A41) in Bruenn 85
  END DO

  CALL LogInterpolateSingleVariable_1D3D_Custom         &
         ( LOG10( Energy ), LOG10( Inte_rho ),          &
           LOG10( Inte_T ), Inte_Ye,                    &
           LOG10( OpacityTable % EnergyGrid % Values ), &
           LOG10( OpacityTable % TS % States(iRho) % Values ), &
           LOG10( OpacityTable % TS % States(iT) % Values ),   &
           OpacityTable % TS % States(iYe) % Values,           &
           Offset_Iso(iNu_e_bar,1), TableES20, bufferIso20 )

  CALL LogInterpolateSingleVariable_1D3D_Custom         &
         ( LOG10( Energy ), LOG10( Inte_rho ),          &
           LOG10( Inte_T ), Inte_Ye,                    &
           LOG10( OpacityTable % EnergyGrid % Values ), &
           LOG10( OpacityTable % TS % States(iRho) % Values ), &
           LOG10( OpacityTable % TS % States(iT) % Values ),   &
           OpacityTable % TS % States(iYe) % Values,           &
           Offset_Iso(iNu_e_bar,2), TableES21, bufferIso21 )

  DO iE = 1, Inte_nPointE
    InterpolantIso2(iE,:) = frpi * Energy(iE)**2 &
                            * ( bufferIso20(iE,:) - bufferIso21(iE,:) / 3.0d0 )
                     ! (A41) in Bruenn 85
  END DO

  END ASSOCIATE ! Table

  CALL DeAllocateOpacityTable( OpacityTable )
  
  !-------------------------------------------------------------------------
  !   write out
  !-------------------------------------------------------------------------
  CALL InitializeHDF( )
  CALL OpenFileHDF( outfilename, .true., file_id )

  CALL OpenGroupHDF( 'ProfileInfo', .true., file_id, group_id )

  datasize1d(1) = Inte_E % nPoints
  CALL WriteHDF( "Energy", Inte_E % Values(:), group_id, datasize1d ) 
  datasize1d(1) = datasize
  CALL WriteHDF( "Radius", Inte_r, group_id, datasize1d )
  CALL WriteHDF( "Density", Inte_rho, group_id, datasize1d )
  CALL WriteHDF( "Temperature", Inte_T, group_id, datasize1d )
  CALL WriteHDF( "Electron Fraction", Inte_ye, group_id, datasize1d )

  CALL CloseGroupHDF( group_id )

  CALL OpenGroupHDF( 'Opacities', .true., file_id, group_id )

  datasize2d(2) = datasize 
  datasize2d(1) = Inte_E % nPoints
  CALL WriteHDF( "Iso Neutrino", InterpolantIso1, group_id, datasize2d )
  CALL WriteHDF( "Iso Antineutrino", InterpolantIso2, group_id, datasize2d )

  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )
  CALL FinalizeHDF( )

  PRINT*, 'Result was written into ', outfilename
  DEALLOCATE( Inte_r, Inte_rho, Inte_T, Inte_Ye, database )
  DEALLOCATE( bufferIso10, bufferIso11, bufferIso20, bufferIso21 )
  DEALLOCATE( InterpolantIso1, InterpolantIso2 )
 
END PROGRAM wlInterpolateIso

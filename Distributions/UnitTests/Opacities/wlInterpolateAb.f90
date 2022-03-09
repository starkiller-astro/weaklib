PROGRAM wlInterpolateAb

  USE wlKindModule, ONLY: dp
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

  !--------- parameters for creating energy grid --------------------------
  INTEGER, PARAMETER     :: Inte_nPointE = 80
  REAL(dp)               :: Inte_Emin = 1.0d-1
  REAL(dp)               :: Inte_Emax = 3.0d02
  TYPE(GridType)         :: Inte_E

  !-------- variables for reading opacity table ---------------------------
  TYPE(OpacityTableType) :: OpacityTable
  REAL(dp), DIMENSION(2) :: Offset_Em

  !-------- variables for reading parameters data -------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Inte_r, Inte_rho, Inte_T, &
                                         Inte_Ye, database
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3
  CHARACTER(LEN=30)                   :: a
  INTEGER                             :: i, datasize

  !-------- variables for output ------------------------------------------
  INTEGER(HID_T)                          :: file_id, group_id
  INTEGER(HSIZE_T)                        :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d

  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: InterpolantEm1, &
                                             InterpolantEm2

  CHARACTER(LEN=128)                      :: EosTableName, OpTableName, &
                                             ProfileName

  CHARACTER(LEN=30)                       :: Outfilename = &
                                            'InterpolatedAbOutput.h5'

  !--------------------------------------------------
  !   take in profile and table names
  !--------------------------------------------------
  OPEN(12, FILE="dataList.txt", FORM = "formatted", &
           ACTION = 'read')
  READ(12,*) ProfileName
  READ(12,*) EosTableName
  READ(12,*) OpTableName
  CLOSE(12, STATUS = 'keep')

  !--------------------------------------------------
  !   interpolated energy 
  !--------------------------------------------------
  Format1 = "(A2,I5)"
  Format2 = "(4A12)"
  Format3 = "(4ES12.3)"

  CALL AllocateGrid( Inte_E, Inte_nPointE )

  Inte_E % Unit = 'MeV                  '
  Inte_E % Name = 'Intepolated Energy   '
  Inte_E % MinValue = Inte_Emin
  Inte_E % MaxValue = Inte_Emax
  Inte_E % LogInterp = 1
  Inte_E % nPoints = Inte_nPointE

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
  ALLOCATE( InterpolantEm1( Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantEm2( Inte_nPointE, datasize ) )

  READ( 1, Format3 ) database
  CLOSE( 1, STATUS = 'keep')  

  DO i = 1, datasize  
    Inte_r(i)   = database(i*4-3)
    Inte_rho(i) = database(i*4-2)
    Inte_T(i)   = database(i*4-1)
    Inte_Ye(i)  = database(i*4)
  END DO 

  !--------------------------------------------
  !    read in the opacity table
  !--------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, &
       FileName_EmAb_Option               &
       = TRIM(OpTableName),               &
       EquationOfStateTableName_Option    &
       = TRIM(EosTableName),              &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )

  Offset_Em = OpacityTable % EmAb % Offsets

  !-----------------------------------------------------------------------
  !   do interpolation
  !-----------------------------------------------------------------------
  ASSOCIATE( TableEm1  => OpacityTable % EmAb % Opacity(1) % Values, &
             TableEm2  => OpacityTable % EmAb % Opacity(2) % Values, &
             Energy    => Inte_E % Values, &
             iRho      => OpacityTable % TS % Indices % iRho, &
             iT        => OpacityTable % TS % Indices % iT,   &
             iYe       => OpacityTable % TS % Indices % iYe   )

  ! interpolate electron neutrino EmAb opacity
  CALL LogInterpolateSingleVariable_1D3D_Custom & 
         ( LOG10( Energy ), LOG10( Inte_rho ),  &
           LOG10( Inte_T ), Inte_Ye,            & 
           LOG10( OpacityTable % EnergyGrid % Values ),        &
           LOG10( OpacityTable % TS % States(iRho) % Values ), &
           LOG10( OpacityTable % TS % States(iT) % Values ),   &
           OpacityTable % TS % States(iYe) % Values,           &
           Offset_Em(iNu_e), TableEm1, InterpolantEm1 )

  ! interpolate electron anti-neutrino EmAb opacity
  CALL LogInterpolateSingleVariable_1D3D_Custom &
         ( LOG10(Energy), LOG10(Inte_rho),      &
           LOG10(Inte_T), Inte_Ye,              &
           LOG10(OpacityTable % EnergyGrid % Values),         &
           LOG10(OpacityTable % TS % States(iRho) % Values),  &
           LOG10(OpacityTable % TS % States(iT) % Values),    &
           OpacityTable % TS % States(iYe) % Values,          &
           Offset_Em(iNu_e_bar), TableEm2, InterpolantEm2 )  

  END ASSOCIATE ! Table

  CALL DeAllocateOpacityTable( OpacityTable )
  
  !----------------------------------------------------------------------
  !   write out
  !----------------------------------------------------------------------
  CALL InitializeHDF( )
  CALL OpenFileHDF( Outfilename, .true., file_id )

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
  CALL WriteHDF( "EmAb_Electron", InterpolantEm1, group_id, datasize2d )
  CALL WriteHDF( "EmAb_ElecAnti", InterpolantEm2, group_id, datasize2d )

  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )
  CALL FinalizeHDF( )
 
  PRINT*, ''
  PRINT*, 'Result was written into ',Outfilename
  DEALLOCATE( Inte_r, Inte_rho, Inte_T, Inte_Ye, database )
  DEALLOCATE( InterpolantEm1, InterpolantEm2 )
 
END PROGRAM wlInterpolateAb

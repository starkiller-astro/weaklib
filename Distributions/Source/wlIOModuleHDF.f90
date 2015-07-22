MODULE wlIOModuleHDF

  USE wlKindModule, ONLY: dp
  USE HDF5 
  USE wlThermoStateModule
  USE wlDependentVariablesModule
  USE wlEquationOfStateTableModule

  implicit none
  PRIVATE
  INTEGER                                     :: hdferr

  PUBLIC InitializeHDF
  PUBLIC FinalizeHDF
  PUBLIC OpenFileHDF
  PUBLIC CloseFileHDF
  PUBLIC OpenGroupHDF
  PUBLIC CloseGroupHDF
  PUBLIC WriteHeaderHDF
  PUBLIC WriteEOSTableHDF
  PUBLIC WriteThermoStateHDF
  PUBLIC WriteDependentVariablesHDF
  PUBLIC ReadThermoStateHDF
  PUBLIC ReadDependentVariablesHDF
  PUBLIC ReadDimensionsHDF
  PUBLIC ReadCHIMERAHDF 
  PUBLIC ReadNumberVariablesHDF
  PUBLIC WriteEquationOfStateTableHDF
  PUBLIC ReadEquationOfStateTableHDF

CONTAINS

  SUBROUTINE InitializeHDF( )

    CALL h5open_f( hdferr )

  END SUBROUTINE InitializeHDF 

  SUBROUTINE FinalizeHDF( )

    CALL h5close_f( hdferr )

  END SUBROUTINE FinalizeHDF 

  SUBROUTINE OpenFileHDF( FileName, NewFile, file_id )  

    CHARACTER(len=*), INTENT(in)               :: FileName
    LOGICAL, INTENT(in)                        :: NewFile
    INTEGER(HID_T), INTENT(out)                :: file_id

    IF ( NewFile ) THEN

      CALL h5fcreate_f( TRIM( FileName ), H5F_ACC_TRUNC_F, file_id, hdferr)

    ELSE

      CALL h5fopen_f( TRIM( FileName ), H5F_ACC_RDONLY_F, file_id, hdferr)
    
    END IF

  END SUBROUTINE OpenFileHDF

  SUBROUTINE CloseFileHDF( file_id )  

    INTEGER(HID_T), INTENT(in)                  :: file_id

    CALL h5fclose_f( file_id, hdferr )

  END SUBROUTINE CloseFileHDF

  SUBROUTINE OpenGroupHDF( GroupName, NewGroup, file_id, group_id )  

    CHARACTER(len=*), INTENT(in)                :: GroupName
    LOGICAL, INTENT(in)                         :: NewGroup
    INTEGER(HID_T), INTENT(in)                  :: file_id
    INTEGER(HID_T), INTENT(out)                 :: group_id

    IF ( NewGroup ) THEN

      CALL h5gcreate_f( file_id, TRIM( GroupName ), group_id, hdferr )

    ELSE

      CALL h5gopen_f( file_id, TRIM( GroupName ), group_id, hdferr )
    
    END IF

  END SUBROUTINE OpenGroupHDF

  SUBROUTINE CloseGroupHDF( group_id )  

    INTEGER(HID_T), INTENT(in)                  :: group_id

    CALL h5gclose_f( group_id, hdferr )

  END SUBROUTINE CloseGroupHDF

  SUBROUTINE WriteHeaderHDF( file_id )

    INTEGER(HID_T), INTENT(in)                  :: file_id
 
  END SUBROUTINE WriteHeaderHDF 

  SUBROUTINE WriteEOSTableHDF( ThermoState, file_id )

    INTEGER(HID_T), INTENT(in)                  :: file_id
    TYPE(ThermoStateType), INTENT(in)           :: ThermoState

    INTEGER(HID_T) :: group_id
   
    CALL h5gcreate_f( file_id, "EOSTable", group_id, hdferr ) 

    CALL WriteThermoStateHDF( ThermoState, group_id ) 
 
    CALL h5gclose_f( group_id, hdferr ) 

  END SUBROUTINE WriteEOSTableHDF

  SUBROUTINE Write1dHDF_double( name, values, group_id, datasize, &
               desc_option, unit_option)

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:), INTENT(in)          :: values
   
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
  
    
    CALL h5screate_simple_f( 1, datasize, dataspace_id, hdferr )

    CALL h5dcreate_f( group_id, name, H5T_NATIVE_DOUBLE, &
           dataspace_id, dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_DOUBLE, &
           values, datasize, hdferr )

    CALL h5sclose_f( dataspace_id, hdferr ) 

    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write1dHDF_double

  SUBROUTINE Read1dHDF_double( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                    :: name
    INTEGER(HID_T)                               :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)   :: datasize
    REAL(dp), DIMENSION(:), INTENT(out)          :: values
    
    INTEGER(HID_T)                               :: dataset_id
  
    CALL h5dopen_f( group_id, name, dataset_id, hdferr )
    CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, &
                   values, datasize, hdferr )
    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Read1dHDF_double
  
  SUBROUTINE Write3dHDF_double( name, values, group_id, datasize, &
               desc_option, unit_option)

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(3), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:,:,:), INTENT(in)      :: values
   
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
  
    
    CALL h5screate_simple_f( 3, datasize, dataspace_id, hdferr )

    CALL h5dcreate_f( group_id, name, H5T_NATIVE_DOUBLE, &
           dataspace_id, dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_DOUBLE, &
           values, datasize, hdferr )

    CALL h5sclose_f( dataspace_id, hdferr ) 

    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write3dHDF_double

  SUBROUTINE Read3dHDF_double( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                     :: name
    INTEGER(HID_T)                               :: group_id
    INTEGER(HSIZE_T), DIMENSION(3), INTENT(in)   :: datasize
    REAL(dp), DIMENSION(:,:,:), INTENT(out)      :: values
    
    INTEGER(HID_T)                               :: dataset_id
  
    CALL h5dopen_f( group_id, name, dataset_id, hdferr )
    CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, &
                   values, datasize, hdferr )
    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Read3dHDF_double

  SUBROUTINE Write3dHDF_integer( name, values, group_id, datasize, &
               desc_option, unit_option)

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(3), INTENT(in)  :: datasize
    INTEGER, DIMENSION(:,:,:), INTENT(in)      :: values
   
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
  
    
    CALL h5screate_simple_f( 3, datasize, dataspace_id, hdferr )

    CALL h5dcreate_f( group_id, name, H5T_NATIVE_INTEGER, &
           dataspace_id, dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_INTEGER, &
           values, datasize, hdferr )

    CALL h5sclose_f( dataspace_id, hdferr ) 

    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write3dHDF_integer

  SUBROUTINE Read3dHDF_integer( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                     :: name
    INTEGER(HID_T)                               :: group_id
    INTEGER(HSIZE_T), DIMENSION(3), INTENT(in)   :: datasize
    INTEGER, DIMENSION(:,:,:), INTENT(out)      :: values
    
    INTEGER(HID_T)                               :: dataset_id
  
    CALL h5dopen_f( group_id, name, dataset_id, hdferr )
    CALL h5dread_f( dataset_id, H5T_NATIVE_INTEGER, &
                   values, datasize, hdferr )
    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Read3dHDF_integer
  
  SUBROUTINE Write1dHDF_integer( name, values, group_id, datasize, &
               desc_option, unit_option)

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)  :: datasize
    INTEGER, DIMENSION(:), INTENT(in)           :: values
   
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
  
    
    CALL h5screate_simple_f( 1, datasize, dataspace_id, hdferr )

    CALL h5dcreate_f( group_id, name, H5T_NATIVE_INTEGER, &
           dataspace_id, dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_INTEGER, &
           values, datasize, hdferr )

    CALL h5sclose_f( dataspace_id, hdferr ) 

    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write1dHDF_integer

  SUBROUTINE Read1dHDF_integer( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                    :: name
    INTEGER(HID_T)                               :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)   :: datasize
    INTEGER, DIMENSION(:), INTENT(out)           :: values
    
    INTEGER(HID_T)                               :: dataset_id
  
    CALL h5dopen_f( group_id, name, dataset_id, hdferr )
    CALL h5dread_f( dataset_id, H5T_NATIVE_INTEGER, &
                   values, datasize, hdferr )
    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Read1dHDF_integer

  SUBROUTINE Write1dHDF_string( name, values, group_id, datasize, &
               desc_option, unit_option)

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)  :: datasize
    CHARACTER(len=*), DIMENSION(:), INTENT(in)  :: values
   
    INTEGER(HSIZE_T)                            :: sizechar
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
  
    
    CALL h5screate_simple_f( 1, datasize, dataspace_id, hdferr )
    sizechar = LEN( values(1) )

    CALL h5tset_size_f( H5T_NATIVE_CHARACTER, sizechar, hdferr )
    CALL h5dcreate_f( group_id, name, H5T_NATIVE_CHARACTER, &
           dataspace_id, dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_CHARACTER, &
           values, datasize, hdferr )

    CALL h5sclose_f( dataspace_id, hdferr ) 

    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write1dHDF_string

  SUBROUTINE Read1dHDF_string( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                      :: name
    INTEGER(HID_T)                                :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)    :: datasize
    CHARACTER(len=*), DIMENSION(:), INTENT(inout) :: values

    INTEGER(HSIZE_T)                              :: sizechar
    INTEGER(HID_T)                                :: dataset_id
  
    sizechar = LEN( values(1) )
    CALL h5tset_size_f( H5T_NATIVE_CHARACTER, sizechar, hdferr )
    CALL h5dopen_f( group_id, name, dataset_id, hdferr )
    CALL h5dread_f( dataset_id, H5T_NATIVE_CHARACTER, &
                   values, datasize, hdferr )
    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Read1dHDF_string
  
  SUBROUTINE WriteThermoStateHDF( TS, group_id )

    TYPE(ThermoStateType), INTENT(in)           :: TS
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER                                     :: i

    datasize1d(1) = 3
    CALL Write1dHDF_integer( "Dimensions", TS % nPoints(:), &
                             group_id, datasize1d )
    
    CALL Write1dHDF_string( "Names", TS % Names(:), &
                             group_id, datasize1d )

    CALL Write1dHDF_string( "Units", TS % Units(:), &
                             group_id, datasize1d )
    DO i = 1, 3
      datasize1d(1) = TS % nPoints(i)
      CALL Write1dHDF_double( TS % Names(i), TS % States(i) % Values(:), &
                              group_id, datasize1d )
    END DO

  END SUBROUTINE WriteThermoStateHDF

  SUBROUTINE WriteDependentVariablesHDF( DV, group_id )

    TYPE(DependentVariablesType), INTENT(in)    :: DV
    INTEGER(HID_T), INTENT(in)                  :: group_id
    
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER                                     :: i

    datasize1d = SIZE( DV % Names )
    CALL Write1dHDF_string( "Names", DV % Names(:), &
                             group_id, datasize1d )

    CALL Write1dHDF_string( "Units", DV % Units(:), &
                             group_id, datasize1d )

    datasize1d = 3
    CALL Write1dHDF_integer( "Dimensions", DV % nPoints(:), &
                             group_id, datasize1d )
    datasize1d = 1 
    CALL Write1dHDF_integer( "nVariables", (/DV % nVariables/), &
                             group_id, datasize1d )
    datasize1d = DV % nVariables 
    CALL Write1dHDF_double( "Offsets", (/DV % Offsets/), &
                             group_id, datasize1d )
    DO i = 1, SIZE( DV % Names ) 
      datasize3d = SHAPE( DV % Variables(i) % Values ) 
      CALL Write3dHDF_double( DV % Names(i), DV % Variables(i) % Values(:,:,:), &
                              group_id, datasize3d )
    END DO

    datasize3d = SHAPE( DV % Repaired )
    CALL Write3dHDF_integer( "Repaired", DV % Repaired(:,:,:), &
                              group_id, datasize3d )

  END SUBROUTINE WriteDependentVariablesHDF

  SUBROUTINE ReadThermoStateHDF( TS, file_id )

    TYPE(ThermoStateType), INTENT(inout)        :: TS
    INTEGER(HID_T), INTENT(in)                  :: file_id

    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER                                     :: i
    
    CALL OpenGroupHDF( "ThermoState", .false., file_id, group_id )

    CALL Read1dHDF_string( "Names", TS % Names(:), &
                              group_id, datasize1d )

    CALL Read1dHDF_string( "Units", TS % Units(:), &
                              group_id, datasize1d )

    DO i = 1, 3
      datasize1d(1) = TS % nPoints(i)
      CALL Read1dHDF_double( TS % Names(i), TS % States(i) % Values(:), &
                              group_id, datasize1d )
      TS % minValues(i) = MINVAL( TS % States(i) % Values(:) )                     
      TS % maxValues(i) = MAXVAL( TS % States(i) % Values(:) )                     
    END DO

    CALL CloseGroupHDF( group_id )

  END SUBROUTINE ReadThermoStateHDF

  SUBROUTINE ReadDependentVariablesHDF( DV, file_id )

    TYPE(DependentVariablesType), INTENT(inout) :: DV
    INTEGER(HID_T), INTENT(in)                  :: file_id

    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER                                     :: i

    CALL OpenGroupHDF( "DependentVariables", .false., file_id, group_id )

    datasize1d = SIZE( DV % Names )
    CALL Read1dHDF_string( "Names", DV % Names(:), &
                              group_id, datasize1d )

    CALL Read1dHDF_string( "Units", DV % Units(:), &
                              group_id, datasize1d )

    DO i = 1, SIZE( DV % Names ) 
      datasize3d = SHAPE( DV % Variables(i) % Values ) 
      CALL Read3dHDF_double( DV % Names(i), DV % Variables(i) % Values(:,:,:), &
                              group_id, datasize3d )
    END DO

    DO i = 1, SIZE( DV % Names ) 
      datasize1d = SIZE( DV % Names )
      CALL Read1dHDF_double( "Offsets", DV % Offsets(:), &
                              group_id, datasize1d )
    END DO

    datasize3d = SHAPE( DV % Repaired )
    CALL Read3dHDF_integer( "Repaired", DV % Repaired(:,:,:), &
                              group_id, datasize3d )

    CALL CloseGroupHDF( group_id )

  END SUBROUTINE ReadDependentVariablesHDF
  
  SUBROUTINE ReadDimensionsHDF ( Dimensions, group_id ) 

    INTEGER(HID_T), INTENT(in)                  :: group_id
    INTEGER, DIMENSION(:), INTENT(inout)        :: Dimensions
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d

    datasize1d(1) = SIZE( Dimensions )
    CALL Read1dHDF_integer( "Dimensions", Dimensions(:), group_id, datasize1d ) 

  END SUBROUTINE ReadDimensionsHDF

  SUBROUTINE ReadNumberVariablesHDF ( nVariables, group_id )

    INTEGER(HID_T), INTENT(in)                  :: group_id
    INTEGER, INTENT(inout)                      :: nVariables
    INTEGER, DIMENSION(1)                       :: nVarTemp  
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d

    datasize1d(1) = 1
    CALL Read1dHDF_integer( "nVariables", nVarTemp(:), group_id, datasize1d )
    nVariables = nVarTemp(1)

  END SUBROUTINE ReadNumberVariablesHDF

  SUBROUTINE WriteEquationOfStateTableHDF( EOSTable )
!  SUBROUTINE WriteEquationOfStateTableHDF( EOSTable, Description  )

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
!    CHARACTER(len=*), INTENT(in) :: Description
    
!    CHARACTER(len=23), PARAMETER :: BaseFileName = "EquationOfStateTable.h5"  
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id

    
    CALL OpenFileHDF( "EquationOfStateTable.h5", .true., file_id )
!    CALL OpenFileHDF( TRIM(ADJUSTL(Description//BaseFileName)), .true., file_id )

    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( EOSTable % TS, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
    CALL WriteDependentVariablesHDF( EOSTable % DV, group_id )
    CALL CloseGroupHDF( group_id )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteEquationOfStateTableHDF
  
  SUBROUTINE ReadEquationOfStateTableHDF( EOSTable, FileName )

    CHARACTER(len=*), INTENT(in)                  :: FileName
    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER                                       :: nVariables
    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "DependentVariables", .false., file_id, group_id )

    CALL ReadDimensionsHDF( nPoints, group_id )
    CALL ReadNumberVariablesHDF( nVariables, group_id )
    CALL CloseGroupHDF( group_id )

    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

    CALL ReadThermoStateHDF( EOSTable % TS, file_id )

    CALL ReadDependentVariablesHDF( EOSTable % DV, file_id )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE ReadEquationOfStateTableHDF

  SUBROUTINE ReadCHIMERAHDF( Rho, T, Ye, E_Int, Entropy, NSE, FileName)

    CHARACTER(len=*), INTENT(in)                :: FileName
    REAL(dp), DIMENSION(1,240,722), INTENT(out) :: Rho
    REAL(dp), DIMENSION(1,240,722), INTENT(out) :: T 
    REAL(dp), DIMENSION(1,240,722), INTENT(out) :: Ye
    REAL(dp), DIMENSION(1,240,722), INTENT(out) :: E_Int
    REAL(dp), DIMENSION(1,240,722), INTENT(out) :: Entropy
    INTEGER, DIMENSION(1,240,722), INTENT(out) :: NSE

    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER(HID_T)                              :: file_id
    INTEGER(HID_T)                              :: group_id

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "fluid", .false., file_id, group_id )

    datasize3d = SHAPE(Rho)
    CALL Read3dHDF_double( "rho_c", Rho(:,:,:), &
                              group_id, datasize3d ) 

    CALL Read3dHDF_double( "t_c", T(:,:,:), &
                              group_id, datasize3d ) 

    CALL Read3dHDF_double( "ye_c", Ye(:,:,:), &
                              group_id, datasize3d ) 

    CALL Read3dHDF_double( "e_int", E_Int(:,:,:), &
                              group_id, datasize3d ) 

    CALL Read3dHDF_double( "entropy", Entropy(:,:,:), &
                              group_id, datasize3d ) 

    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "abundance", .false., file_id, group_id )

    datasize3d = SHAPE(NSE)
    CALL Read3dHDF_integer( "nse_c", NSE(:,:,:), &
                              group_id, datasize3d )

    CALL CloseGroupHDF( group_id )

  END SUBROUTINE ReadCHIMERAHDF


END MODULE wlIOModuleHDF


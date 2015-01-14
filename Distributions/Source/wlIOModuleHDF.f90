MODULE wlIOModuleHDF

  USE wlKindModule, ONLY: dp
  USE HDF5 
  USE wlThermoStateModule
  USE wlDependentVariablesModule

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
  PUBLIC LoadThermoStateHDF

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
    ! Write EOS name, EOS parameters w/ names and comments
 
  END SUBROUTINE WriteHeaderHDF 

  SUBROUTINE WriteEOSTableHDF( ThermoState, file_id )

    INTEGER(HID_T), INTENT(in)                  :: file_id
    TYPE(ThermoStateType), INTENT(in)           :: ThermoState

    INTEGER(HID_T) :: group_id
   
    CALL h5gcreate_f( file_id, "EOSTable", group_id, hdferr ) 

    CALL WriteThermoStateHDF( ThermoState, group_id ) 
 
    ! CALL WriteDependentVariablesHDF 

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
    CALL Write1dHDF_integer( "Dimensions", TS % nValues(:), &
                             group_id, datasize1d )
    
    CALL Write1dHDF_string( "Names", TS % Names(:), &
                             group_id, datasize1d )
    DO i = 1, 3
      datasize1d(1) = TS % nValues(i)
      CALL Write1dHDF_double( TS % Names(i), TS % States(i) % Values(:), &
                              group_id, datasize1d )
    END DO

  END SUBROUTINE WriteThermoStateHDF

  SUBROUTINE WriteDependentVariablesHDF( DV, group_id )

    TYPE(DependentVariablesType), INTENT(in)    :: DV
    INTEGER(HID_T), INTENT(in)                  :: group_id
    
    INTEGER, DIMENSION(12)                      :: nValues
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER                                     :: i

    datasize1d = SIZE( DV % Names )
    CALL Write1dHDF_string( "Names", DV % Names(:), &
                             group_id, datasize1d )
    nValues = (/10,10,10,10,10,10,10,10,10,10,10,10/)
     
    CALL Write1dHDF_integer( "Dimensions", nValues(:), &
                             group_id, datasize1d )
    DO i = 1, SIZE( DV % Names ) 
      datasize3d = SHAPE( DV % Variables(i) % Values ) 
      CALL Write3dHDF_double( DV % Names(i), DV % Variables(i) % Values(:,:,:), &
                              group_id, datasize3d )
    END DO

  END SUBROUTINE WriteDependentVariablesHDF

  SUBROUTINE ReadThermoStateHDF( TS, group_id )

    TYPE(ThermoStateType), INTENT(inout)        :: TS
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER                                     :: i
    
    CALL Read1dHDF_string( "Names", TS % Names(:), &
                              group_id, datasize1d )
    DO i = 1, 3
      datasize1d(1) = TS % nValues(i)
      CALL Read1dHDF_double( TS % Names(i), TS % States(i) % Values(:), &
                              group_id, datasize1d )
      TS % minValues(i) = MINVAL( TS % States(i) % Values(:) )                     
      TS % maxValues(i) = MAXVAL( TS % States(i) % Values(:) )                     
    END DO

  END SUBROUTINE ReadThermoStateHDF

  SUBROUTINE ReadDependentVariablesHDF( DV, group_id )

    TYPE(DependentVariablesType), INTENT(inout) :: DV
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER                                     :: i


    datasize1d = SIZE( DV % Names )
    CALL Read1dHDF_string( "Names", DV % Names(:), &
                              group_id, datasize1d )
    !DO i = 1, 3
    !  datasize1d(1) = DV % nValues(i)
    !  CALL Read1dHDF_double( DV % Names(i), DV % States(i) % Values(:), &
    !                          group_id, datasize1d )
      !DV % minValues(i) = MINVAL( DV % States(i) % Values(:) )
      !DV % maxValues(i) = MAXVAL( DV % States(i) % Values(:) )
    !END DO

    DO i = 1, 3
      datasize3d = SHAPE( DV % Variables(i) % Values ) 
      CALL Read3dHDF_double( DV % Names(i), DV % Variables(i) % Values(:,:,:), &
                              group_id, datasize3d )
    END DO

  END SUBROUTINE ReadDependentVariablesHDF
  
  SUBROUTINE ReadDimensionsHDF ( Dimensions, group_id ) 

    INTEGER(HID_T), INTENT(in)                  :: group_id
    INTEGER, DIMENSION(:), INTENT(inout)        :: Dimensions
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d

    datasize1d(1) = SIZE( Dimensions )
    CALL Read1dHDF_integer( "Dimensions", Dimensions(:), group_id, datasize1d ) 

  END SUBROUTINE ReadDimensionsHDF
  
  SUBROUTINE LoadThermoStateHDF( TS, file_id )

    TYPE(ThermoStateType), INTENT(inout)        :: TS
    INTEGER, DIMENSION(3)                       :: npts
    INTEGER(HID_T)                              :: group_id
    INTEGER(HID_T)                              :: file_id

  CALL OpenGroupHDF( "ThermoState", .false., file_id, group_id )

  CALL ReadDimensionsHDF( npts, group_id )
  CALL AllocateThermoState( TS, npts )

  TS % nValues(1:3) = npts(1:3)

  CALL ReadThermoStateHDF( TS, group_id )
  CALL CloseGroupHDF( group_id )

  END SUBROUTINE LoadThermoStateHDF



END MODULE wlIOModuleHDF


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
  PUBLIC DescribeEquationOfStateTable
  PUBLIC ReadEquationOfStateTableParallelHDF

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
  
  SUBROUTINE read_1d_slab_int(name, value, group_id, datasize)
    CHARACTER(*), INTENT(IN)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(1), INTENT(IN)  :: datasize
    INTEGER, DIMENSION(:), INTENT(OUT)          :: value

    INTEGER(HID_T)                              :: dataset_id
    INTEGER                                     :: error

    CALL h5dopen_f(group_id, name, dataset_id, error)
    CALL h5dread_f(dataset_id, H5T_NATIVE_INTEGER, &
                    value, datasize, error)
    CALL h5dclose_f(dataset_id, error)

  END SUBROUTINE read_1d_slab_int

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
    INTEGER, DIMENSION(1)                       :: buffer

  WRITE (*,*) "Starting HDF TS write "

    datasize1d(1) = 3
    CALL Write1dHDF_integer( "Dimensions", TS % nPoints(:), &
                             group_id, datasize1d )
    
    CALL Write1dHDF_integer( "LogInterp", TS % LogInterp(:), &
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
 
  WRITE (*,*) "Starting HDF indicies write "
    datasize1d = 1
    buffer(1) = TS % Indices % iRho
    CALL Write1dHDF_integer( "iRho", buffer, &
                             group_id, datasize1d )
    
    buffer(1) = TS % Indices % iT
    CALL Write1dHDF_integer( "iT", buffer, &
                             group_id, datasize1d )

    buffer(1) = TS % Indices % iYe
    CALL Write1dHDF_integer( "iYe", buffer, &
                             group_id, datasize1d )

  WRITE (*,*) "HDF indicies write successful "

  END SUBROUTINE WriteThermoStateHDF

  SUBROUTINE WriteDependentVariablesHDF( DV, group_id )

    TYPE(DependentVariablesType), INTENT(in)    :: DV
    INTEGER(HID_T), INTENT(in)                  :: group_id
    
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer

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

    datasize1d = 1
    buffer(1) = DV % Indices % iPressure 
    CALL Write1dHDF_integer( "iPressure", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iEntropyPerBaryon 
    CALL Write1dHDF_integer( "iEntropyPerBaryon", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iInternalEnergyDensity 
    CALL Write1dHDF_integer( "iInternalEnergyDensity", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iElectronChemicalPotential 
    CALL Write1dHDF_integer( "iElectronChemicalPotential", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iProtonChemicalPotential 
    CALL Write1dHDF_integer( "iProtonChemicalPotential", buffer, & 
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iNeutronChemicalPotential 
    CALL Write1dHDF_integer( "iNeutronChemicalPotential", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iProtonMassFraction 
    CALL Write1dHDF_integer( "iProtonMassFraction", buffer, & 
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iNeutronMassFraction 
    CALL Write1dHDF_integer( "iNeutronMassFraction", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iAlphaMassFraction 
    CALL Write1dHDF_integer( "iAlphaMassFraction", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iHeavyMassFraction 
    CALL Write1dHDF_integer( "iHeavyMassFraction", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iHeavyChargeNumber 
    CALL Write1dHDF_integer( "iHeavyChargeNumber", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iHeavyMassNumber 
    CALL Write1dHDF_integer( "iHeavyMassNumber", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iHeavyBindingEnergy 
    CALL Write1dHDF_integer( "iHeavyBindingEnergy", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iThermalEnergy 
    CALL Write1dHDF_integer( "iThermalEnergy", buffer, &
                             group_id, datasize1d )

    buffer(1) = DV % Indices % iGamma1 
    CALL Write1dHDF_integer( "iGamma1", buffer, &
                             group_id, datasize1d )

  END SUBROUTINE WriteDependentVariablesHDF

  SUBROUTINE ReadThermoStateHDF( TS, file_id )

    TYPE(ThermoStateType), INTENT(inout)        :: TS
    INTEGER(HID_T), INTENT(in)                  :: file_id

    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer 
    
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

    datasize1d(1) = 1
    CALL Read1dHDF_integer( "iRho", buffer, group_id, datasize1d )
    TS % Indices % iRho = buffer(1)

    CALL Read1dHDF_integer( "iT", buffer, group_id, datasize1d )
    TS % Indices % iT = buffer(1)

    CALL Read1dHDF_integer( "iYe", buffer, group_id, datasize1d )
    TS % Indices % iYe = buffer(1)

    datasize1d(1) = SIZE( TS % LogInterp )
    CALL Read1dHDF_integer( "LogInterp", TS % LogInterp(:), group_id, datasize1d )

    CALL CloseGroupHDF( group_id )

  END SUBROUTINE ReadThermoStateHDF

  SUBROUTINE ReadDependentVariablesHDF( DV, file_id )

    TYPE(DependentVariablesType), INTENT(inout) :: DV
    INTEGER(HID_T), INTENT(in)                  :: file_id

    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer

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

    datasize1d(1) = 1
    CALL Read1dHDF_integer( "iPressure", buffer, &
                             group_id, datasize1d )
    DV % Indices % iPressure = buffer(1)

    CALL Read1dHDF_integer( "iEntropyPerBaryon", buffer, &
                             group_id, datasize1d )
    DV % Indices % iEntropyPerBaryon = buffer(1)

    CALL Read1dHDF_integer( "iInternalEnergyDensity", buffer, &
                             group_id, datasize1d )
    DV % Indices % iInternalEnergyDensity = buffer(1)

    CALL Read1dHDF_integer( "iElectronChemicalPotential", buffer, &
                             group_id, datasize1d )
    DV % Indices % iElectronChemicalPotential = buffer(1)

    CALL Read1dHDF_integer( "iProtonChemicalPotential", buffer, &
                             group_id, datasize1d )
    DV % Indices % iProtonChemicalPotential = buffer(1)

    CALL Read1dHDF_integer( "iNeutronChemicalPotential", buffer, &
                             group_id, datasize1d )
    DV % Indices % iNeutronChemicalPotential = buffer(1)

    CALL Read1dHDF_integer( "iProtonMassFraction", buffer, &
                             group_id, datasize1d )
    DV % Indices % iProtonMassFraction = buffer(1)

    CALL Read1dHDF_integer( "iNeutronMassFraction", buffer, &
                             group_id, datasize1d )
    DV % Indices % iNeutronMassFraction = buffer(1)

    CALL Read1dHDF_integer( "iAlphaMassFraction", buffer, &
                             group_id, datasize1d )
    DV % Indices % iAlphaMassFraction = buffer(1)

    CALL Read1dHDF_integer( "iHeavyMassFraction", buffer, &
                             group_id, datasize1d )
    DV % Indices % iHeavyMassFraction = buffer(1)

    CALL Read1dHDF_integer( "iHeavyChargeNumber", buffer, &
                             group_id, datasize1d )
    DV % Indices % iHeavyChargeNumber = buffer(1)

    CALL Read1dHDF_integer( "iHeavyMassNumber", buffer, &
                             group_id, datasize1d )
    DV % Indices % iHeavyMassNumber = buffer(1)

    CALL Read1dHDF_integer( "iHeavyBindingEnergy", buffer, &
                             group_id, datasize1d )
    DV % Indices % iHeavyBindingEnergy = buffer(1)

    CALL Read1dHDF_integer( "iThermalEnergy", buffer, &
                             group_id, datasize1d )
    DV % Indices % iThermalEnergy = buffer(1)

    CALL Read1dHDF_integer( "iGamma1", buffer, &
                             group_id, datasize1d )
    DV % Indices % iGamma1 = buffer(1)

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

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id
    
    CALL OpenFileHDF( "EquationOfStateTable.h5", .true., file_id )

    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( EOSTable % TS, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
    CALL WriteDependentVariablesHDF( EOSTable % DV, group_id )
    CALL CloseGroupHDF( group_id )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteEquationOfStateTableHDF

  SUBROUTINE DescribeEquationOfStateTable( EOSTable )

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    INTEGER :: i

    DO i = 1, 3
      WRITE (*,*) "Independent Variable", i, "=", &
                  EOSTable % TS % Names(i), "Units:", EOSTable % TS % Units(i)
      WRITE (*,*) "Independent Variable", i, " Minimum = " , EOSTable % TS % minValues(i)
      WRITE (*,*) "Independent Variable", i, "Maximum =" , EOSTable % TS % maxValues(i)
      WRITE (*,*) "Number of Independent Variable", i, "Points =" , EOSTable % TS % nPoints(i) 
      IF ( EOSTable % TS % LogInterp(i) == 1 ) THEN
         WRITE (*,*) "Independent Variable ", i, "Grid Logarithmically Spaced" 
         ELSE
         WRITE (*,*) "Independent Variable", i, "Grid Linearly Spaced" 
      END IF  
    END DO

    DO i = 1, EOSTable % nVariables
      WRITE (*,*) "Dependent Variable", i, "=", &
                  EOSTable % DV % Names(i), "Units:", EOSTable % DV % Units(i)
    END DO

  END SUBROUTINE DescribeEquationOfStateTable
  
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

    CALL DescribeEquationOfStateTable( EOSTable )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE ReadEquationOfStateTableHDF

  SUBROUTINE ReadEquationOfStateTableParallelHDF( EOSTable, FileName, rootproc, COMMUNICATOR )

    USE MPI
    
    implicit none

    CHARACTER(len=*), INTENT(in)                  :: FileName
    INTEGER, INTENT(in)                           :: rootproc
    INTEGER, INTENT(in)                           :: COMMUNICATOR
    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER, DIMENSION(4)                         :: buffer
    INTEGER                                       :: nStates, nVariables, i
    INTEGER                                       :: i_count
    INTEGER                                       :: myid, ierr 
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id

    CALL MPI_COMM_RANK( COMMUNICATOR, myid , ierr )

    IF ( myid == rootproc ) THEN

      CALL OpenFileHDF( FileName, .false., file_id )

      CALL OpenGroupHDF( "DependentVariables", .false., file_id, group_id )

      CALL ReadDimensionsHDF( nPoints, group_id )
      CALL ReadNumberVariablesHDF( nVariables, group_id )
      CALL CloseGroupHDF( group_id )
    
      buffer(1) = nPoints(1)
      buffer(2) = nPoints(2)
      buffer(3) = nPoints(3)
      buffer(4) = nVariables

  WRITE (*,*) "in: process", myid, "buffer(4)", buffer(4) 

    END IF

    i_count = SIZE(buffer) 

    CALL MPI_BCAST( buffer, i_count, MPI_INTEGER, rootproc, COMMUNICATOR, ierr )

    IF ( myid /= rootproc ) THEN
      
      nPoints(1) = buffer(1)
      nPoints(2) = buffer(2)
      nPoints(3) = buffer(3)
      nVariables = buffer(4) 

  WRITE (*,*) "out process", myid, "buffer(4)", nVariables 

    END IF

    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

    IF ( myid == rootproc ) THEN

      CALL ReadThermoStateHDF( EOSTable % TS, file_id )

      CALL ReadDependentVariablesHDF( EOSTable % DV, file_id )

      CALL CloseFileHDF( file_id )

    END IF

    i_count = PRODUCT(nPoints) 
    
    nStates = 3

    DO i= 1, nStates
      CALL MPI_BCAST(EOSTable % TS % States(i) % Values(:), nPoints(i), &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    END DO

    CALL MPI_BCAST(EOSTable % TS % Names(:), nStates,     &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % Units(:), nStates,     &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % minValues(:), nStates, &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % maxValues(:), nStates, &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )

  WRITE (*,*) "process", myid, "test TS data", EOSTable % TS % States(1) % Values(10)
  WRITE (*,*) "process", myid, "test TS data", EOSTable % TS % Names(1) 
  WRITE (*,*) "process", myid, "test TS data", EOSTable % TS % Units(1) 

    DO i= 1, nVariables  
      CALL MPI_BCAST(EOSTable % DV % Variables(i) % Values(:,:,:), i_count,  &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    END DO

    CALL MPI_BCAST(EOSTable % DV % Names(:), nVariables,                   &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Units(:), nVariables,                   &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Offsets(:), nVariables,                 &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Repaired(:,:,:), i_count,               &
                   MPI_INTEGER, rootproc, COMMUNICATOR, ierr )

  END SUBROUTINE ReadEquationOfStateTableParallelHDF

  SUBROUTINE ReadCHIMERAHDF( Rho, T, Ye, E_Int, Entropy, NSE, imax, nx, ny, &
                             nz, FileName)

    CHARACTER(len=*), INTENT(in)                :: FileName
    INTEGER, INTENT(out) :: imax, nx, ny, nz
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Rho
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: T 
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Ye
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: E_Int
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Entropy
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: NSE
    INTEGER, DIMENSION(2) :: indices

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER(HID_T)                              :: file_id
    INTEGER(HID_T)                              :: group_id

    CALL OpenFileHDF( FileName, .false., file_id )
 
    CALL OpenGroupHDF( "mesh", .false., file_id, group_id )

    datasize1d(1) = 2
    CALL read_1d_slab_int('radial_index_bound', indices, group_id, &
           datasize1d)
    imax = indices(2)
    nx = imax + 2

    CALL read_1d_slab_int('theta_index_bound', indices, group_id, &
           datasize1d)
    ny = indices(2)

    CALL read_1d_slab_int('phi_index_bound', indices, group_id, &
           datasize1d)
    nz = indices(2)

    CALL CloseGroupHDF( group_id )

    ALLOCATE( Rho( nx, ny, nz ), T( nx, ny, nz ), Ye( nx, ny, nz ),           &
              E_Int( nx, ny, nz ), Entropy( nx, ny, nz ), NSE( nx + 1, ny, nz ) ) 

    CALL OpenGroupHDF( "fluid", .false., file_id, group_id )

    datasize3d = (/nx,ny,nz/)
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


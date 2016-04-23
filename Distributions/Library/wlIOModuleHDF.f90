MODULE wlIOModuleHDF

  USE wlKindModule, ONLY: dp
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  USE HDF5 

  IMPLICIT NONE
  PRIVATE

  INTEGER :: hdferr

  PUBLIC InitializeHDF
  PUBLIC FinalizeHDF
  PUBLIC OpenFileHDF
  PUBLIC CloseFileHDF
  PUBLIC OpenGroupHDF
  PUBLIC CloseGroupHDF
  PUBLIC WriteHeaderHDF
  PUBLIC Write1dHDF_double
  PUBLIC Read1dHDF_double
  PUBLIC Write3dHDF_double
  PUBLIC Read3dHDF_double
  PUBLIC Write3dHDF_integer
  PUBLIC Read3dHDF_integer
  PUBLIC read_1d_slab_int
  PUBLIC Write1dHDF_integer
  PUBLIC Read1dHDF_integer
  PUBLIC Write1dHDF_string
  PUBLIC Read1dHDF_string
  PUBLIC WriteThermoStateHDF
  PUBLIC WriteDependentVariablesHDF
  PUBLIC ReadThermoStateHDF
  PUBLIC ReadDependentVariablesHDF
  PUBLIC ReadDimensionsHDF
  PUBLIC ReadNumberVariablesHDF

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
  
  SUBROUTINE Write3dHDF_double &
               ( name, values, group_id, datasize, desc_option, unit_option )   

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

    CALL h5dcreate_f( group_id, name, H5T_NATIVE_DOUBLE, dataspace_id, &
                      dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_DOUBLE, values, datasize, hdferr )

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
    CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, values, datasize, hdferr )
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

    datasize1d(1) = 3
    CALL Write1dHDF_integer &
           ( "Dimensions", TS % nPoints(:), group_id, datasize1d )
    
    CALL Write1dHDF_integer &
           ( "LogInterp", TS % LogInterp(:), group_id, datasize1d )
    
    CALL Write1dHDF_string &
           ( "Names", TS % Names(:), group_id, datasize1d )

    CALL Write1dHDF_string &
           ( "Units", TS % Units(:), group_id, datasize1d )

    DO i = 1, 3
      datasize1d(1) = TS % nPoints(i)
      CALL Write1dHDF_double &
             ( TS % Names(i), TS % States(i) % Values(:), group_id, datasize1d )
    END DO
 
    datasize1d = 1
    buffer(1) = TS % Indices % iRho
    CALL Write1dHDF_integer( "iRho", buffer, group_id, datasize1d )
    
    buffer(1) = TS % Indices % iT
    CALL Write1dHDF_integer( "iT",   buffer, group_id, datasize1d )

    buffer(1) = TS % Indices % iYe
    CALL Write1dHDF_integer( "iYe",  buffer, group_id, datasize1d )

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

END MODULE wlIOModuleHDF


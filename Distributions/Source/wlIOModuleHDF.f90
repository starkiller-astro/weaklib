MODULE wlIOModuleHDF

  USE wlKindModule, ONLY: dp
  USE HDF5 
  USE wlThermoStateModule

  implicit none
  PRIVATE
  INTEGER        :: hdferr
  INTEGER        :: rank




  PUBLIC ReadGrid 
  PUBLIC WriteGrid

CONTAINS

  SUBROUTINE InitializeHDF( )

    CALL h5open_f(hdferr)

  END SUBROUTINE InitializeHDF 

  SUBROUTINE OpenFileHDF( FileName, NewFile, fileID )  

    CHARACTER(len=*), INTENT(in) :: FileName
    LOGICAL, INTENT(in) :: NewFile
    INTEGER(HID_T), INTENT(out) :: fileID

    IF ( NewFile ) THEN

      CALL h5fcreate_f( TRIM( FileName ), H5F_ACC_TRUNC_F, fileID, hdferr)

    ELSE

      CALL h5fopen_f( TRIM( FileName ), H5F_ACC_RDONLY_F, fileID, hdferr)
    
    END IF

  END SUBROUTINE OpenFileHDF

  SUBROUTINE WriteHeaderHDF( fileID )

    INTEGER(HID_T), INTENT(in) :: fileID
    ! Write EOS name, EOS parameters w/ names and comments
    
 
  END SUBROUTINE WriteHeaderHDF 

  SUBROUTINE WriteEOSTableHDF( ThermoState, fileID )

    INTEGER(HID_T), INTENT(in) :: fileID
    TYPE(ThermoStateType), INTENT(in)         :: ThermoState

    INTEGER(HID_T) :: groupID
   
    CALL h5gcreate_f( fileID, "EOSTable", groupID, hdferr) 

    CALL WriteThermoStateHDF( ThermoState, groupID ) 
 
    ! CALL WriteDependentVariablesHDF 

    CALL h5gclose_f( groupID, hdferr) 

  END SUBROUTINE WriteEOSTableHDF

  SUBROUTINE WriteThermoStateHDF( ThermoState, groupID )
  
    INTEGER(HID_T), INTENT(in) :: groupID
    TYPE(ThermoStateType), INTENT(in)         :: ThermoState

    INTEGER                                   :: i
    INTEGER, INTENT(in)                       :: rank
    INTEGER, DIMENSION(rank), INTENT(in)      :: griddims
    INTEGER(HID_T) :: dset_id
    INTEGER(HID_T) :: dspace_id

    DO i = 1, rank
       CALL h5dcreate_f(file_id, TRIM( ThermoState1 % Names(i) ),&
         H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)

       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,&
         Thermostate1 % States(i) % Values(:), griddims(i), hdferr)

       CALL h5dclose_f(dset_id, hdferr)

     END DO


  END SUBROUTINE WriteThermoStateHDF

  SUBROUTINE WriteGridHDF( )  ! Gets passed griddims, buffer, returns h5grid
   DO i = 1, griddims( )
     CALL h5screate_simple_f(3, griddims, space, hdferr) 
     CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, buffer, griddims, hdferr)
     CALL h5dclose_f(dset0, hdferr)
   END DO
  END SUBROUTINE WriteGrid

END MODULE wlIOModuleHDF


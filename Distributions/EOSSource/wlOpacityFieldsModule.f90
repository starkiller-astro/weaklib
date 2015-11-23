MODULE wlOpacityFieldsModule
!-------------------------------------------------------------
!   Module: wlOpacityFieldsModule
!   Author: R. Landfield
!   Date:   11/6/14
!-------------------------------------------------------------

  USE HDF5
  
!-------------------------------------------------------------
! Allocate Opacity Fields    
!-------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

!Char List
  INTEGER, PUBLIC, PARAMETER :: N_OPACITY_FIELDS = 2
  CHARACTER(32), PUBLIC, DIMENSION(N_OPACITY_FIELDS) ::&
                 OpacityFieldNames = (/'Absorption                      ', & 
                                       'Emission                        '/)

  INTEGER, PUBLIC, PARAMETER :: Absorption                   = 1
  INTEGER, PUBLIC, PARAMETER :: Emission                     = 2
  REAL(8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: OpacityData
  REAL(8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: Metadata
  CHARACTER(LEN=15), PARAMETER  :: filename   = "opacitytable.h5"
  INTEGER(HID_T), PUBLIC     :: file_id               ! HDF5 File identifier 
  INTEGER(HID_T), PUBLIC     :: dataset_id            ! HDF5 dataset identifier
  INTEGER(HID_T), PUBLIC     :: dataspace_id          ! HDF5 dataspace identifier
  INTEGER                    :: hdferr                ! Error flag
  INTEGER(HSIZE_T), ALLOCATABLE, DIMENSION(1:4) ::&
              bdims = (/dim3, dim0, dim1, dim2/)      ! size Opacity read buffer
  INTEGER(HSIZE_T), DIMENSION(1:2)   ::&
              mddims = (/dim4, dim5/)                 ! size metadata write buffer

  PUBLIC :: &
    AllocateOpacityFields, &
    DeAllocateOpacityFields, &
    ReadOpacityFields, &
    ReadOpacityMetadata 
  
CONTAINS

  Subroutine ReadOpacityMetadata
  !-------------------------------------------------------------
  !   Open the file and the metadata dataset.
  !-------------------------------------------------------------
    CALL h5open_f(hdferr)
    CALL h5fopen_f(filename, H5F_ACC_RDONLY, file_id, hdferr);
    CALL h5dopen_f(file_id, "metadata", dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, MetaData(:,1), mddims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)
   metadata(1,1) = Nrho
   metadata(1,2) = Nt
   metadata(1,3) = Nye
   metadata(1,4) = rhomin ! Min rho
   metadata(1,5) = rhomax !Max rho
   metadata(1,6) = Tmin   !Min T  
   metadata(1,7) = Tmax   !Max T  
   metadata(1,8) = Yemin  !Min Ye  
   metadata(1,9) = Yemax  !Max Ye
   metadata(1,10) = 12.d0    ! Number of EOS fields
   metadata(1,11) = 2.d0     !Number of Opacity fields
   metadata(1,12) = 1.d0     !Number of Neutrino Species
   metadata(1,13) = nez   !Number of Energy groups
   metadata(1,14) = Emin
   metadata(1,15) = Emax

  End Subroutine ReadOpacityMetadata
  
  Subroutine AllocateOpacityFields
  !-------------------------------------------------------------
  ! Initialize HDF5 FORTRAN interface.
  !-------------------------------------------------------------
    CALL h5open_f(filename, H5F_ACC_RDONLY, file_id, hdferr)
      ALLOCATE(OpacityData(Nrho,Nt,Nye,N_OPACITY_FIELDS))


    End Subroutine AllocateEosFields


  Subroutine ReadOpacityFields
  !-------------------------------------------------------------
  !   Open the file and the dataset.
  !-------------------------------------------------------------

    CALL h5open_f(hdferr)
    CALL h5fopen_f(filename, H5F_ACC_RDONLY, file_id, hdferr)

      DO i=1,N_OPACITY_FIELDS;
        CALL h5dopen_f(file_id, OpacityFieldNames(i), dataset_id, hdferr)
        CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE,&
                OpacityData(:,:,:,:,i), bdims, hdferr)
        CALL h5dclose_f(dataset_id, hdferr) 
      END DO

    CALL h5fclose_f(filename, H5F_ACC_RDONLY, file_id, hdferr);

    CALL h5close_f(hdferr)

  End Subroutine ReadOpacityFields

  Subroutine DeAllocateEosFields
  !-------------------------------------------------------------
  !   Close/release resources.
  !-------------------------------------------------------------
    h5dclose(dataset);
    h5sclose(dataspace);
    h5sclose(memspace);
    h5fclose(file);
  End Subroutine DeAllocateEosFields

End Module EosFieldsModule


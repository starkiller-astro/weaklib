MODULE EosFieldsModule
!-------------------------------------------------------------
!   Module: EosFieldsModule
!   Author: R. Landfield
!   Date:   11/6/14
!-------------------------------------------------------------

USE opacitytable.h5
USE HDF5


!-------------------------------------------------------------
! Allocate EOS Fields    
!-------------------------------------------------------------

IMPLICIT NONE
PRIVATE

!Char List
INTEGER, PUBLIC, PARAMETER :: N_EOS_FIELDS = 12
CHARACTER(32), PUBLIC, DIMENSION(N_EOS_FIELDS) ::&
                  EosFieldNames = (/'Pressure                        ', & 
                                    'InternalEnergyDensity           ', &
                                    'EntropyPerBaryon                ', &
                                    'NeutronChemicalPotential        ', &
                                    'ProtonChemicalPotential         ', &
                                    'ElectronChemicalPotential       ', &
                                    'NeutronMassFraction             ', &
                                    'ProtonMassFraction              ', &
                                    'HeliumMassFraction              ', &
                                    'HeavyMassFraction               ', &
                                    'HeavyMassNumber                 ', &
                                    'HeavyChargeNumber               '/)

INTEGER, PUBLIC, PARAMETER :: Pressure                   = 1
INTEGER, PUBLIC, PARAMETER :: InternalEnergyDensity      = 2
INTEGER, PUBLIC, PARAMETER :: EntropyPerBaryon           = 3
INTEGER, PUBLIC, PARAMETER :: NeutronChemicalPotential   = 4
INTEGER, PUBLIC, PARAMETER :: ProtonChemicalPotential    = 5
INTEGER, PUBLIC, PARAMETER :: ElectronChemicalPotential  = 6
INTEGER, PUBLIC, PARAMETER :: NeutronMassFraction        = 7
INTEGER, PUBLIC, PARAMETER :: ProtonMassFraction         = 8
INTEGER, PUBLIC, PARAMETER :: HeliumMassFraction         = 9
INTEGER, PUBLIC, PARAMETER :: HeavyMassFraction          = 10
INTEGER, PUBLIC, PARAMETER :: HeavyMassNumber            = 11
INTEGER, PUBLIC, PARAMETER :: HeavyChargeNumber          = 12
REAL, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) :: EosData
REAL, PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: Metadata
CHARACTER(LEN=15), PARAMETER  :: filename   = "opacitytable.h5"
INTEGER(HID_T), PUBLIC     :: file_id               ! HDF5 File identifier 
INTEGER(HID_T), PUBLIC     :: dataset_id            ! HDF5 dataset identifier
INTEGER(HID_T), PUBLIC     :: dataspace_id          ! HDF5 dataspace identifier
INTEGER                    :: hdferr                ! Error flag
INTEGER(HSIZE_T), ALLOCATABLE, DIMENSION(1:3) ::&
            dims = (/dim0, dim1, dim2/)             ! Size EOS read buffers
INTEGER(HSIZE_T), ALLOCATABLE, DIMENSION(1:4) ::&
            bdims = (/dim3, dim0, dim1, dim2/)      ! size Opacity read buffer
INTEGER(HSIZE_T), DIMENSION(1:2)   ::&
            mddims = (/dim4, dim5/)                 ! size metadata write buffer

PUBLIC :: &
  AllocateEosFields, &
  DeAllocateEosFields, &
  ReadEosFields, &
  ReadEosMetadata 
  
CONTAINS

  Subroutine ReadEosMetadata
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

  End Subroutine ReadEosMetadata
  
  Subroutine AllocateEosFields
  !-------------------------------------------------------------
  ! Initialize HDF5 FORTRAN interface.
  !-------------------------------------------------------------
    CALL h5open_f(filename, H5F_ACC_RDONLY, file_id, hdferr)
      ALLOCATE(GridData( 
      ALLOCATE(EosData(Nrho,Nt,Nye,N_EOS_FIELDS))


  End Subroutine AllocateEosFields


  Subroutine ReadEosFields
  !-------------------------------------------------------------
  !   Open the file and the dataset.
  !-------------------------------------------------------------

    CALL h5dopen_f(file_id, EosFieldNames(1), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,1), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr) 

    CALL h5dopen_f(file_id, EosFieldNames(2), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,2), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(3), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,3), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(4), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,4), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(5), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,5), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(6), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,6), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(7), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,7), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(8), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,8), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(9), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,9), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(10), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,10), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(11), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,11), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)

    CALL h5dopen_f(file_id, EosFieldNames(12), dataset_id, hdferr)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, EosData(:,:,:,12), dims, hdferr)
    CALL h5dclose_f(dataset_id, hdferr)



  End Subroutine ReadEosFields

  Subroutine DeAllocateEosFields
  !-------------------------------------------------------------
  !   Close/release resources.
  !-------------------------------------------------------------
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Fclose(file);
  End Subroutine DeAllocateEosFields

End Module EosFieldsModule


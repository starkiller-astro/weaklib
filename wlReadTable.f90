PROGRAM wlreadtable_get
       
!--------------------------------------------------
!    File:         wlreadtable_get
!    Module:       wlreadtable_get
!    AuthorType:   R. Landfield
!    Type:         Program
!
!    Date:         10/23/14
!
!    Purpose:
!        To query HDF5 precalculated table of EOS state variables
!        and Weak interaction values produced by wltable_get.f90     
!
!--------------------------------------------------

USE HDF5
USE wlEosFieldsModule
USE wlOpacityFieldsModule
USE wlStateModule
REAL(KIND=double), DIMENSION(Nrho*Nt*Nye), INTENT(inout) :: rho       ! density [g/cm**3]
REAL(KIND=double), DIMENSION(Nrho*Nt*Nye), INTENT(inout) :: ye        ! electron fraction
REAL(KIND=double), DIMENSION(Nrho*Nt*Nye), INTENT(inout) :: t         ! temperature [K]
REAL(KIND=double), DIMENSION(12, Nrho*Nt*Nye), INTENT(out) :: EosOutput   ! (rho, T, Ye, E, P, S,..) 
REAL(KIND=double), DIMENSION(2, Nrho*Nt*Nye), INTENT(out) :: OpacityOutput   ! (rho, T, Ye, E, P, S,..) 
INTEGER(HID_T), PUBLIC     :: file_id               ! HDF5 File identifier 
INTEGER(HID_T), PUBLIC     :: dataset_id            ! HDF5 dataset identifier
INTEGER(HID_T), PUBLIC     :: dataspace_id          ! HDF5 dataspace identifier
INTEGER                    :: hdferr                ! Error flag



!--------------------------------------------------
! Initialize HDF5
!--------------------------------------------------

   CALL h5open_f(hdferr)

!   CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)
!   At some point, we'll need to have a rho, t, y query call
!, i.e. a subroutiine above this that calls this program
   CALL h5fopen_f(opacitytable.h5, H5F_ACC_RDONLY_F, file, hdferr)

   CALL AllocateState(:,:)
     ! Inputs = ?
   CALL ReadEosMetadata
     ! Inputs = opacity.h5 
     ! Outputs = Eos Metadata
   CALL AllocateEosFields
     ! Inputs = Eos Metadata
     ! Outputs = none, it just allocates the table 
   CALL ReadEosFields
     ! Inputs = opacityfields.h5, Eos Metadata
     ! Outputs = local buffers
! Write Field values
   CALL ReadOpacityMetadata

   CALL AllocateOpacityFields

   CALL ReadOpacityFields
! Write Field values
   CALL DeAllocateState(:,:)
   
   CALL DeAllocateEosFields

   CALL DeAllocateOpacityFields
!allocate Nrho, Nye, Nt

!allocate each EOS variable from the AllEOSFields(nvar), nvar 1-12

   REAL, INTENT(in), ALLOCATABLE, DIMENSION(:,:) 
   CALL h5fclose_f(file , hdferr)



END PROGRAM wlreadtable_get

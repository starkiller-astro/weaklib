PROGRAM wlDependentVariablesTest
 
  USE wlDependentVariablesModule
  USE HDF5
  USE wlThermoStateModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, OpenFileHDF, OpenGroupHDF,         &
                           WriteDependentVariablesHDF, CloseGroupHDF,        & 
                           ReadDependentVariablesHDF, CloseFileHDF,          &
                           FinalizeHDF, ReadDimensionsHDF
  implicit none

  INTEGER :: i
  INTEGER, DIMENSION(12) :: npts
  INTEGER :: j
  TYPE(DependentVariablesType) :: DV
  TYPE(DependentVariablesType) :: DV2
  INTEGER(HID_T) :: file_id
  INTEGER(HID_T) :: group_id

print*,"1"

  CALL AllocateDependentVariables( DV, nValues = (/10,10,10,10,10,10,10,10,10,10,10,10/), nVariables = 12 )
  !CALL AllocateDependentVariables( DV, nValues = (/10,10,10/), nVariables = 3 )
  !CALL AllocateDependentVariables( DV, nVariables = 3 )

print*,"2"

  DV % Names(1:12) = (/'Pressure                        ', &
                       'Entropy Per Baryon              ', &
                       'Internal Energy Density         ', &
                       'Neutron Chemical Potential      ', &
                       'Electron Chemical Potential     ', &
                       'Proton Chemical Potential       ', &
                       'Neutron Mass Fraction           ', &
                       'Proton Mass Fraction            ', &
                       'Helium Mass Fraction            ', &
                       'Heavy Mass Fraction             ', &
                       'Heavy Mass Number               ', &
                       'Heavy Charge Number             '/)

  !DV % nValues(1:3) = (/10,10,10/)


print*,"3"

  DO i = 1, SIZE( DV % Variables )
    WRITE (*,*) SHAPE( DV % Variables(i) % Values )
  END DO

  DO i = 1, SIZE( DV % Variables )
    DV % Variables(i) % Values = i 
  END DO



print*,"4"

  CALL InitializeHDF( )
  CALL OpenFileHDF( "DependentVariablesFile.h5", .true., file_id )
  CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
  CALL WriteDependentVariablesHDF( DV, group_id )
  CALL CloseGroupHDF( group_id )
  CALL CloseFileHDF( file_id )
  CALL DeAllocateDependentVariables( DV )

print*,"5"

  CALL OpenFileHDF( "DependentVariablesFile.h5", .false., file_id )
  !CALL LoadDependentVariablesHDF( DV2, file_id )
  CALL OpenGroupHDF( "DependentVariables", .false., file_id, group_id )
  CALL ReadDimensionsHDF( npts, group_id )
  CALL AllocateDependentVariables( DV2, npts, nVariables = 12 )

  !DV2 % nValues(1:3) = npts(1:3)

  CALL ReadDependentVariablesHDF( DV2, group_id )
  CALL CloseGroupHDF( group_id )
  CALL CloseFileHDF( file_id )

  DO j = 1,12
    WRITE(*,*) TRIM( DV2 % Names(j) )
    !WRITE(*,*) DV2 % nValues(j)
    WRITE(*,*) DV2 % Variables(j) % Values(:,:,:)
  END DO

  CALL DeAllocateDependentVariables( DV2 )

  CALL FinalizeHDF( )
 
END PROGRAM wlDependentVariablesTest

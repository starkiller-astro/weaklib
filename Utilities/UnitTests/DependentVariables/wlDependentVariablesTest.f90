PROGRAM wlDependentVariablesTest
 
  USE wlDependentVariablesModule
  USE HDF5
  USE wlThermoStateModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, OpenFileHDF, OpenGroupHDF,         &
                           WriteDependentVariablesHDF,   &
                           CloseGroupHDF, CloseFileHDF, FinalizeHDF,         &
                           ReadDimensionsHDF
  implicit none

  INTEGER :: i
  TYPE(DependentVariablesType) :: DV
  INTEGER(HID_T) :: file_id
  INTEGER(HID_T) :: group_id

print*,"1"

  CALL AllocateDependentVariables( DV, nValues = (/10,10,10/), nVariables = 3 )
  !CALL AllocateDependentVariables( DV, nVariables = 3 )

print*,"2"

  DV % Names(1:3) = (/'Pressure                        ', &
                      'Entropy Per Baryon              ', &
                      'Internal Energy Density         '/)

  !DV % nValues(1:3) = (/10,10,10/)


print*,"3"

  DO i = 1, SIZE( DV % Variables )
    WRITE (*,*) SHAPE( DV % Variables(i) % Values )
  END DO

print*,"4"

  CALL InitializeHDF( )
  CALL OpenFileHDF( "DependentVariablesFile.h5", .true., file_id )
  CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
  CALL WriteDependentVariablesHDF( DV, group_id )
  CALL CloseGroupHDF( group_id )
  CALL CloseFileHDF( file_id )



  CALL DeAllocateDependentVariables( DV )
 
END PROGRAM wlDependentVariablesTest

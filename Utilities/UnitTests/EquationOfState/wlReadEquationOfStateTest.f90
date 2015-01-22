PROGRAM wlReadEquationOfStateTest
 
  USE wlDependentVariablesModule
  USE HDF5
  USE wlThermoStateModule
  USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, OpenFileHDF, OpenGroupHDF,          &
                           WriteDependentVariablesHDF, CloseGroupHDF,         & 
                           ReadDependentVariablesHDF, CloseFileHDF,           &
                           FinalizeHDF, ReadDimensionsHDF,                    &
                           LoadDependentVariablesHDF, ReadNumberVariablesHDF, &
                           WriteThermoStateHDF, LoadThermoStateHDF,           &
                           ReadThermoStateHDF                                    

  implicit none

  INTEGER                        :: i
  INTEGER, DIMENSION(3)          :: nPoints
  INTEGER                        :: nVariables
  INTEGER                        :: j
  TYPE(EquationOfStateTableType) :: EOSTable
  !TYPE(DependentVariablesType)   :: DV, DV2 
  INTEGER(HID_T)                 :: file_id
  INTEGER(HID_T)                 :: group_id


  !nPoints = (/10,10,10/)
  !nVariables = 12
  CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

!STOP

  CALL InitializeHDF( )

  CALL OpenFileHDF( "EquationOfStateTable.h5", .false., file_id )

  CALL LoadThermoStateHDF( EOSTable % TS, file_id ) 

  CALL LoadDependentVariablesHDF( EOSTable % DV, file_id )

  CALL CloseFileHDF( file_id )
!  CALL DeAllocateThermoState( EOSTable % TS )
!  CALL DeAllocateDependentVariables( EOSTable % DV )

  DO i = 1, SIZE( EOSTable % DV % Variables )
    WRITE (*,*) SHAPE( EOSTable % DV % Variables(i) % Values )
  END DO

  DO i = 1, SIZE( EOSTable % DV % Variables )
    EOSTable % DV % Variables(i) % Values = i
  END DO

  DO j = 1,3
    WRITE(*,*) TRIM( EOSTable % TS % Names(j) )
    WRITE(*,*) EOSTable % TS % nPoints(j)
    WRITE(*,*) EOSTable % TS % minValues(j), EOSTable % TS % maxValues(j)
    WRITE(*,*) EOSTable % TS % States(j) % Values(:)
  END DO

  DO j = 1,12
    WRITE(*,*)
    WRITE(*,*) TRIM( EOSTable % DV % Names(j) ) , j
    WRITE(*,*)
    WRITE(*,*) EOSTable % DV % Variables(j) % Values(:,:,:)
  END DO

  CALL DeAllocateEquationOfStateTable( EOSTable )
  CALL FinalizeHDF( )

END PROGRAM wlReadEquationOfStateTest

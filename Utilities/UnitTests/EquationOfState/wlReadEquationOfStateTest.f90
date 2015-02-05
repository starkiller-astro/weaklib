PROGRAM wlReadEquationOfStateTest
 
  USE HDF5
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF, & 
                           ReadEquationOfStateTableHDF 

  implicit none

  INTEGER                        :: i
  TYPE(EquationOfStateTableType) :: EOSTable

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  Write (*,*) EOSTable % nPoints
  Write (*,*) EOSTable % nVariables

  DO i = 1, SIZE( EOSTable % DV % Variables )
    WRITE (*,*) SHAPE( EOSTable % DV % Variables(i) % Values )
  END DO

  DO i = 1, SIZE( EOSTable % DV % Variables )
    EOSTable % DV % Variables(i) % Values = i
  END DO

  DO i = 1,3
    WRITE(*,*) TRIM( EOSTable % TS % Names(i) )
    WRITE(*,*) EOSTable % TS % nPoints(i)
    WRITE(*,*) EOSTable % TS % minValues(i), EOSTable % TS % maxValues(i)
    WRITE(*,*) EOSTable % TS % States(i) % Values(:)
  END DO

  DO i = 1,13
    WRITE(*,*)
    WRITE(*,*) TRIM( EOSTable % DV % Names(i) ) , i
    WRITE(*,*)
    WRITE(*,*) EOSTable % DV % Variables(i) % Values(:,:,:)
  END DO

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlReadEquationOfStateTest

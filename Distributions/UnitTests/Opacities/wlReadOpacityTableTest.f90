PROGRAM wlReadOpacityTableTest

  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DescribeOpacityTable
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF_New

  IMPLICIT NONE

  TYPE(OpacityTableType) :: OpacityTable

  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF_New &
          ( OpacityTable, "temp_EmAb.h5", &
            ReadOpacity_EmAb_Option = .TRUE. )
  CALL FinalizeHDF( ) 

  PRINT*
  PRINT*, "Describe OpacityTable"

  CALL DescribeOpacityTable( OpacityTable )

END PROGRAM wlReadOpacityTableTest

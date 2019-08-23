PROGRAM wlReadOpacityTableTest

  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DescribeOpacityTable
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF

  IMPLICIT NONE

  TYPE(OpacityTableType) :: OpacityTable

  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF &
          ( OpacityTable, FileName_EmAb_Option = "temp_EmAb.h5" )
  CALL FinalizeHDF( ) 

  PRINT*
  PRINT*, "Describe OpacityTable"

  CALL DescribeOpacityTable( OpacityTable )

END PROGRAM wlReadOpacityTableTest

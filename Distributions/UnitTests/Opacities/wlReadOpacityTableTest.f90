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
  CALL ReadOpacityTableHDF( OpacityTable, "OpTa_Grey_5quad_lowEOS.h5" )
  CALL FinalizeHDF( ) 

  PRINT*, "Describe OpacityTable"

  CALL DescribeOpacityTable( OpacityTable )

END PROGRAM wlReadOpacityTableTest

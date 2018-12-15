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
          ( OpacityTable, "temp_Pair.h5", &
            ReadOpacity_Pair_Option = .TRUE. )
  CALL FinalizeHDF( ) 

  PRINT*
  PRINT*, "Describe OpacityTable"

  CALL DescribeOpacityTable( OpacityTable )

END PROGRAM wlReadOpacityTableTest

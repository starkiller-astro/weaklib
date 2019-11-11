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
          ( OpacityTable, &
           FileName_EmAb_Option &
           = "wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5", &
           EquationOfStateTableName_Option = "wl-EOS-SFHo-15-25-50.h5" )
  CALL FinalizeHDF( ) 

  PRINT*
  PRINT*, "Describe OpacityTable"

  CALL DescribeOpacityTable( OpacityTable )

END PROGRAM wlReadOpacityTableTest

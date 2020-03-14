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
           = "", &
           FileName_Iso_Option &
           = "", &
           FileName_NES_Option &
           = "wl-Op-SFHo-25-20-100-E40-B85-NES.h5", &
           FileName_Pair_Option &
           = "", &
           EquationOfStateTableName_Option &
           = "wl-EOS-SFHo-25-40-100.h5" )

  CALL FinalizeHDF( ) 

  PRINT*
  PRINT*, "Describe OpacityTable"

  CALL DescribeOpacityTable( OpacityTable )

END PROGRAM wlReadOpacityTableTest

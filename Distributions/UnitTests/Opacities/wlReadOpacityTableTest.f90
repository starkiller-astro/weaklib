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
           = "wl-Op-LS220-15-25-50-Lower-T-E40-B85-AbEm.h5", &
           FileName_Iso_Option &
           = "", &
           FileName_NES_Option &
           = "", &
           FileName_Pair_Option &
           = "", &
           EquationOfStateTableName_Option &
           = "wl-EOS-LS220-15-25-50-Lower-T-rewrite.h5" )

  CALL FinalizeHDF( ) 

  PRINT*
  PRINT*, "Describe OpacityTable"

  CALL DescribeOpacityTable( OpacityTable )

END PROGRAM wlReadOpacityTableTest

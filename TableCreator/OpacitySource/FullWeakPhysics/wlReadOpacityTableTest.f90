PROGRAM wlReadOpacityTableTest

  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DescribeOpacityTable
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF

  IMPLICIT NONE

  TYPE(OpacityTableType) :: OpacityTable

  CALL InitializeHDF( )

  CALL ReadOpacityTableHDF &
          ( OpacityTable, &
           FileName_EmAb_Option &
           = "wl-Op-SFHo-15-25-50-E40-EmAb.h5", &
           FileName_Iso_Option &
           = "wl-Op-SFHo-15-25-50-E40-Iso.h5", &
           FileName_NNS_Option &
           = "wl-Op-SFHo-15-25-50-E40-NNS.h5", &
           FileName_NES_Option &
           = "wl-Op-SFHo-15-25-50-E40-NES.h5", &
           FileName_Pair_Option &
           = "wl-Op-SFHo-15-25-50-E40-Pair.h5", &
           FileName_Brem_Option &
           = "wl-Op-SFHo-15-25-50-E40-Brem.h5", &
           EquationOfStateTableName_Option &
           = "wl-EOS-SFHo-15-25-50.h5" )

  CALL DescribeOpacityTable( OpacityTable )

  CALL FinalizeHDF( ) 

END PROGRAM wlReadOpacityTableTest

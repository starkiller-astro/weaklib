PROGRAM wlReadOpacityTableTest
   
  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlGridModule, ONLY: MakeLogGrid
  USE wlThermoStateModule
  USE wlDependentVariablesModule
  USE wlOpacityTableModule

  USE wlIOModuleHDF, ONLY:&
          InitializeHDF, FinalizeHDF
  USE wlOpacityTableModule
  USE wlOpacityTableIOModuleHDF, ONLY: &
          ReadOpacityTableHDF
  USE wlExtPhysicalConstantsModule

  implicit none

  TYPE(OpacityTableType)  :: OpacityTable

  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( ) 

  PRINT*, "Describe OpacityTable"

  CALL DescribeOpacityTable( OpacityTable )

END PROGRAM wlReadOpacityTableTest

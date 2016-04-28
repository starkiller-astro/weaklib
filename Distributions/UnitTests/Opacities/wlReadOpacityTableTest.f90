PROGRAM wlReadOpacityTableTest
   
  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlGridModule, ONLY: MakeLogGrid
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  USE wlIOModuleHDF!, ONLY:&
       !   InitializeHDF, FinalizeHDF
  USE wlOpacityTableModule
  USE wlOpacityTableIOModuleHDF!, ONLY: &
        !  ReadOpacityTableHDF
  USE wlExtPhysicalConstantsModule

  implicit none

  TYPE(OpacityTableType)  :: OpacityTable

  CALL InitializeHDF( )

  CALL ReadOpacityTableHDF( OpacityTable, "TestTable.h5" )

!  Write (*,*) OpacityTable % EOSTable % nPoints
!  Write (*,*) OpacityTable % EOSTable % nVariables
!  Write (*,*) OpacityTable % EnergyGrid
!  Write (*,*) OpacityTable % EcapEm % Names

  CALL FinalizeHDF( )  

END PROGRAM wlReadOpacityTableTest

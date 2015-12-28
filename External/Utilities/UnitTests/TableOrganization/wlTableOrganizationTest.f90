PROGRAM wlTableOrganizationTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF

  implicit none

  TYPE(EquationOfStateTableType) :: OldEOSTable
  TYPE(EquationOfStateTableType) :: NewEOSTable
  TYPE(DVIDType)                 :: NewDVID
  INTEGER                        :: NewnVariables


    NewnVariables = 14

    NewDVID % iPressure = 1
    NewDVID % iEntropyPerBaryon = 3
    NewDVID % iInternalEnergyDensity = 2
    NewDVID % iElectronChemicalPotential = 6
    NewDVID % iProtonChemicalPotential = 5
    NewDVID % iNeutronChemicalPotential = 4
    NewDVID % iProtonMassFraction = 8
    NewDVID % iNeutronMassFraction = 7 
    NewDVID % iAlphaMassFraction = 0
    NewDVID % iHeavyMassFraction = 9
    NewDVID % iHeavyChargeNumber = 11
    NewDVID % iHeavyMassNumber = 10
    NewDVID % iHeavyBindingEnergy = 13
    NewDVID % iThermalEnergy = 14
    NewDVID % iGamma1 = 12


  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( OldEOSTable, "EquationOfStateTable.h5" )

  CALL MatchTableStructure( OldEOSTable, NewEOSTable, NewDVID, NewnVariables )

  CALL DescribeEquationOfStateTable( NewEOSTable )
  WRITE (*,*) "P=", NewEOSTable % DV % Indices % iPressure
  WRITE (*,*) "S=", NewEOSTable % DV % Indices % iEntropyPerBaryon
  WRITE (*,*) "E=", NewEOSTable % DV % Indices % iInternalEnergyDensity
  WRITE (*,*) "EChem=", NewEOSTable % DV % Indices % iElectronChemicalPotential
  WRITE (*,*) "PChem=", NewEOSTable % DV % Indices % iProtonChemicalPotential
  WRITE (*,*) "NChemE=", NewEOSTable % DV % Indices % iNeutronChemicalPotential
  WRITE (*,*) "PFrac=", NewEOSTable % DV % Indices % iProtonMassFraction
  WRITE (*,*) "NFrac=", NewEOSTable % DV % Indices % iNeutronMassFraction
  !WRITE (*,*) "AFrac=", NewEOSTable % DV % Indices % iAlphaMassFraction
  WRITE (*,*) "HFrac=", NewEOSTable % DV % Indices % iHeavyMassFraction
  WRITE (*,*) "HChar=", NewEOSTable % DV % Indices % iHeavyChargeNumber
  WRITE (*,*) "HMass#=", NewEOSTable % DV % Indices % iHeavyMassNumber
  WRITE (*,*) "HBE=", NewEOSTable % DV % Indices % iHeavyBindingEnergy
  WRITE (*,*) "TE=", NewEOSTable % DV % Indices % iThermalEnergy
  WRITE (*,*) "G=", NewEOSTable % DV % Indices % iGamma1 

  CALL DeAllocateEquationOfStateTable( NewEOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlTableOrganizationTest

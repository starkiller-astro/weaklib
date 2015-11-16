PROGRAM wlTableOrganizationTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF

  implicit none

  TYPE(EquationOfStateTableType) :: EOSTable
  TYPE(DVIDType) :: LocalDVID


    LocalDVID % iPressure = 1
    LocalDVID % iEntropyPerBaryon = 3
    LocalDVID % iInternalEnergyDensity = 2
    LocalDVID % iElectronChemicalPotential = 6
    LocalDVID % iProtonChemicalPotential = 5
    LocalDVID % iNeutronChemicalPotential = 4
    LocalDVID % iProtonMassFraction = 8
    LocalDVID % iNeutronMassFraction = 7 
    LocalDVID % iAlphaMassFraction = 15
    LocalDVID % iHeavyMassFraction = 9
    LocalDVID % iHeavyChargeNumber = 11
    LocalDVID % iHeavyMassNumber = 10
    LocalDVID % iHeavyBindingEnergy = 13
    LocalDVID % iThermalEnergy = 14
    LocalDVID % iGamma1 = 12


  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  CALL MatchTableStructure( EOSTable, LocalDVID )

  CALL DescribeEquationOfStateTable( EOSTable )
  WRITE (*,*) "P=", EOSTable % DV % Indices % iPressure
  WRITE (*,*) "S=", EOSTable % DV % Indices % iEntropyPerBaryon
  WRITE (*,*) "E=", EOSTable % DV % Indices % iInternalEnergyDensity
  WRITE (*,*) "EChem=", EOSTable % DV % Indices % iElectronChemicalPotential
  WRITE (*,*) "PChem=", EOSTable % DV % Indices % iProtonChemicalPotential
  WRITE (*,*) "NChemE=", EOSTable % DV % Indices % iNeutronChemicalPotential
  WRITE (*,*) "PFrac=", EOSTable % DV % Indices % iProtonMassFraction
  WRITE (*,*) "NFrac=", EOSTable % DV % Indices % iNeutronMassFraction
  WRITE (*,*) "AFrac=", EOSTable % DV % Indices % iAlphaMassFraction
  WRITE (*,*) "HFrac=", EOSTable % DV % Indices % iHeavyMassFraction
  WRITE (*,*) "HChar=", EOSTable % DV % Indices % iHeavyChargeNumber
  WRITE (*,*) "HMass#=", EOSTable % DV % Indices % iHeavyMassNumber
  WRITE (*,*) "HBE=", EOSTable % DV % Indices % iHeavyBindingEnergy
  WRITE (*,*) "TE=", EOSTable % DV % Indices % iThermalEnergy
  WRITE (*,*) "G=", EOSTable % DV % Indices % iGamma1 

 
! The routine will go through a do loop, and will check if the indicies match.
! If they don't, it will write the DV data to a buffer, put what index the local DVID 
! wants in that position to that positon, and write what's in the buffer to the 
! position that the DV index position the local DVID has it in.
! If the local DVID doesn't use that particular variable, it's not re-written.  
! In this test, the local indicies will be populated 

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlTableOrganizationTest

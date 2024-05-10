PROGRAM wlCreateEquationOfStateTable
    
    USE wlKindModule, ONLY: dp
    USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
    USE HDF5
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
    USE wlEOSIOModuleHDF
    USE wlExtNumericalModule, ONLY: zero
    USE wlCompOSEInterface, ONLY : ReadnPointsFromCompOSE, ReadCompOSETable, &
	RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE
    USE wlHelmholtzIOModuleHDF, ONLY : WriteHelmholtzTableHDF
    USE wlElectronEOSModule, ONLY: HelmholtzEOSType, iTempMax, iDenMax, &
        AllocateHelmEOS, DeAllocateHelmEOS, ReadHelmEOSdat
	
    IMPLICIT NONE
    
    INTEGER                        :: iVars
    INTEGER, DIMENSION(3)          :: nPointsBaryon
    INTEGER, DIMENSION(2)          :: nPointsHelm
    INTEGER                        :: nVariables
    TYPE(EquationOfStateTableType) :: EOSBaryonTable
    TYPE(HelmholtzEOSType) :: HelmEOSTable
    
    CHARACTER(len=128) :: CompOSEFilePath, HelmDatFilePath, BaryonEOSTableName, HelmEOSTableName
    
    REAL(dp) :: Minimum_Value
    
    BaryonEOSTableName = 'BaryonsEOS.h5'
    HelmEOSTableName = 'HelmEOS.h5'
    CompOSEFilePath = '../../../../../../../Compare_different_tables/SFHo_no_ele/'
    CALL ReadnPointsFromCompOSE(CompOSEFilePath, nPointsBaryon(1), nPointsBaryon(2), nPointsBaryon(3))
    nVariables = 15
    
    ! Fill the rho, temperature, ye, and compose EOS table arrays
    CALL ReadCompOSETable( CompOSEFilePath, nPointsBaryon(1), nPointsBaryon(2), nPointsBaryon(3) )
    
    ! -------------------- NOW DO BARYONIC EOS ----------------------------------------------
    PRINT*, "Allocate Baryonic EOS"
    CALL AllocateEquationOfStateTable( EOSBaryonTable, nPointsBaryon , nVariables )
    
    EOSBaryonTable % TS % Names(1:3) = (/'Density                         ',&
    'Temperature                     ',&
    'Proton Fraction                 '/)
    
    EOSBaryonTable % TS % Indices % iRho = 1
    EOSBaryonTable % TS % Indices % iT   = 2
    EOSBaryonTable % TS % Indices % iYe  = 3
    
    PRINT*, "Allocate Independent Variable Units " 
    
    EOSBaryonTable % TS % Units(1:3) = (/'Grams per cm^3                  ', &
    'K                               ', &
    '                                '/) 
    
    ! DO WE NEED THIS?
    EOSBaryonTable % TS % minValues(1:3) =  (/MINVAL(RhoCompOSE), MINVAL(TempCompOSE), MINVAL(YpCompOSE)/)
    EOSBaryonTable % TS % maxValues(1:3) =  (/MAXVAL(RhoCompOSE), MAXVAL(TempCompOSE), MAXVAL(YpCompOSE)/)
    
    PRINT*, "Label Grid Type"
    EOSBaryonTable % TS % LogInterp(1:3) =  (/1, 1, 0/)
    
    PRINT*, "Allocate Names " 
    EOSBaryonTable % DV % Names(1:15) = (/'Pressure                        ', &
    'Entropy Per Baryon              ', &
    'Internal Energy Density         ', &
    'Electron Chemical Potential     ', &
    'Proton Chemical Potential       ', &
    'Neutron Chemical Potential      ', &
    'Proton Mass Fraction            ', &
    'Neutron Mass Fraction           ', &
    'Alpha Mass Fraction             ', &
    'Heavy Mass Fraction             ', &
    'Heavy Charge Number             ', &
    'Heavy Mass Number               ', &
    'Heavy Binding Energy            ', &
    'Thermal Energy                  ', &
    'Gamma1                          '/)
    
    PRINT*, "Set Dependent Variable Identifier Indicies " 
    
    EOSBaryonTable % DV % Indices % iPressure = 1
    EOSBaryonTable % DV % Indices % iEntropyPerBaryon = 2
    EOSBaryonTable % DV % Indices % iInternalEnergyDensity = 3
    EOSBaryonTable % DV % Indices % iElectronChemicalPotential = 4
    EOSBaryonTable % DV % Indices % iProtonChemicalPotential = 5
    EOSBaryonTable % DV % Indices % iNeutronChemicalPotential = 6
    EOSBaryonTable % DV % Indices % iProtonMassFraction = 7
    EOSBaryonTable % DV % Indices % iNeutronMassFraction = 8
    EOSBaryonTable % DV % Indices % iAlphaMassFraction = 9
    EOSBaryonTable % DV % Indices % iHeavyMassFraction = 10
    EOSBaryonTable % DV % Indices % iHeavyChargeNumber = 11
    EOSBaryonTable % DV % Indices % iHeavyMassNumber = 12
    EOSBaryonTable % DV % Indices % iHeavyBindingEnergy = 13
    EOSBaryonTable % DV % Indices % iThermalEnergy = 14
    EOSBaryonTable % DV % Indices % iGamma1 = 15
    
    WRITE (*,*) "iRho, iGamma1", EOSBaryonTable % TS % Indices % iRho, EOSBaryonTable % DV % Indices % iGamma1
    
    PRINT*, "Allocate Dependent Variable Units " 
    EOSBaryonTable % DV % Units(1:15) = (/'Dynes per cm^2                  ', &
    'k_b per baryon                  ', &
    'erg per gram                    ', &
    'MeV                             ', &
    'MeV                             ', &
    'MeV                             ', &
    '                                ', &
    '                                ', &
    '                                ', &
    '                                ', &
    '                                ', &
    '                                ', &
    'MeV                             ', &
    'MeV                             ', &
    '                                '/)
    
    PRINT*, "Begin Populating EOSBaryonTable" 
    ! FROM HERE ON YOU HAVE TO CONVERT THE COMPOSE EOS TO WEAKLIB
    EOSBaryonTable % TS % States(1) % Values = RhoCompOSE
    EOSBaryonTable % TS % States(2) % Values = TempCompOSE
    EOSBaryonTable % TS % States(3) % Values = YpCompOSE
    
    ! thermo and compo state
    EOSBaryonTable % DV % Variables(1) % Values(:,:,:) = EOSCompOSE(:,:,:,1)
    EOSBaryonTable % DV % Variables(2) % Values(:,:,:) = EOSCompOSE(:,:,:,2)
    EOSBaryonTable % DV % Variables(3) % Values(:,:,:) = EOSCompOSE(:,:,:,3)
    EOSBaryonTable % DV % Variables(4) % Values(:,:,:) = EOSCompOSE(:,:,:,4)
    EOSBaryonTable % DV % Variables(5) % Values(:,:,:) = EOSCompOSE(:,:,:,5)
    EOSBaryonTable % DV % Variables(6) % Values(:,:,:) = EOSCompOSE(:,:,:,6)
    EOSBaryonTable % DV % Variables(7) % Values(:,:,:) = EOSCompOSE(:,:,:,7)
    EOSBaryonTable % DV % Variables(8) % Values(:,:,:) = EOSCompOSE(:,:,:,8)
    EOSBaryonTable % DV % Variables(9) % Values(:,:,:) = EOSCompOSE(:,:,:,9)
    EOSBaryonTable % DV % Variables(10) % Values(:,:,:) = EOSCompOSE(:,:,:,10)
    EOSBaryonTable % DV % Variables(11) % Values(:,:,:) = EOSCompOSE(:,:,:,11)
    EOSBaryonTable % DV % Variables(12) % Values(:,:,:) = EOSCompOSE(:,:,:,12)
    EOSBaryonTable % DV % Variables(13) % Values(:,:,:) = EOSCompOSE(:,:,:,13)
    EOSBaryonTable % DV % Variables(14) % Values(:,:,:) = EOSCompOSE(:,:,:,14)
    EOSBaryonTable % DV % Variables(15) % Values(:,:,:) = EOSCompOSE(:,:,:,15)
    
    PRINT*, "Offset negative quantities"
    
    DO iVars = 1, EOSBaryonTable % DV % nVariables
        Minimum_Value = MINVAL(EOSBaryonTable % DV % Variables(iVars) % Values)
        IF ( Minimum_Value .LT. zero) THEN
            EOSBaryonTable % DV % Offsets(iVars) = -1.1d0*Minimum_Value
            EOSBaryonTable % DV % Variables(iVars) % Values =  &
        EOSBaryonTable % DV % Variables(iVars) % Values -1.1d0*Minimum_Value
        ELSE
            EOSBaryonTable % DV % Offsets(iVars) = zero
        END IF
    END DO
    
    ! ------------- NOW DO ELECTRON EOS ------------------ !
    nPointsHelm = (/ iTempMax, iDenMax /)
    PRINT*, "Allocate Helmholtz EOS"
    CALL AllocateHelmEOS( HelmEOSTable, nPointsHelm )
    
    HelmDatFilePath = '../helm_table.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmEOSTable )
    
    ! NOW CREATE BARYONIC FILE
    
    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Baryonic EOS"
    CALL WriteEquationOfStateTableHDF( EOSBaryonTable, BaryonEOSTableName )
        
    ! NOW ADD Helmholtz EOS TO PREVIOUSLY CREATED H5 FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Appending Helmholtz EOS to HDF file"
    CALL WriteHelmholtzTableHDF( HelmEOSTable, HelmEOSTableName, .true. )
    CALL FinalizeHDF( )
        
    WRITE (*,*) "HDF write successful"
    
    CALL DeAllocateEquationOfStateTable( EOSBaryonTable )
    CALL DeAllocateHelmEOS( HelmEOSTable )
    
END PROGRAM wlCreateEquationOfStateTable

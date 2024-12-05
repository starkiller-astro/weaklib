PROGRAM wlCreateEquationOfStateTable
    
    USE wlKindModule, ONLY: dp
    USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
    USE HDF5
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
    USE wlEOSIOModuleHDF
    USE wlCompOSEInterface, ONLY : ReadnPointsFromCompOSE, ReadCompOSETable, &
        ReadCompOSEHDFTable, RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE
    USE wlHelmMuonIOModuleHDF, ONLY : WriteHelmholtzTableHDF, WriteMuonTableHDF
    USE wlLeptonEOSModule, ONLY: HelmholtzTableType, MuonEOSType, &
        iTempMax, iDenMax, &
        nTempMuon, nDenMuon, &
        ReadHelmEOSdat, ReadMuonEOSdat, &
        AllocateHelmholtzTable, DeallocateHelmholtzTable, &
        AllocateMuonEOS, DeAllocateMuonEOS 

    IMPLICIT NONE
    
    INTEGER                        :: iVars
    INTEGER, DIMENSION(3)          :: nPointsBaryon
    INTEGER, DIMENSION(2)          :: nPointsHelm, nPointsDenon
    INTEGER                        :: nVariables
    TYPE(EquationOfStateTableType) :: EOSTable
    TYPE(HelmholtzTableType) :: HelmEOSTable
    TYPE(MuonEOSType) :: MuonEOSTable
    
    CHARACTER(len=128) :: CompOSEFilePath, CompOSEFHDF5Path, &
        HelmDatFilePath, MuonDatFilePath, &
        BaryonEOSTableName, HelmEOSTableName, MuonEOSTableName
    
    REAL(dp) :: Minimum_Value
    LOGICAL  :: ReadFullTable, RedHDF5Table, ResetNegativePressure
    INTEGER  :: iLepton, iRho, iTemp, iYe

    ReadFullTable = .false.
    RedHDF5Table = .false.
    ResetNegativePressure = .false.
    
    IF (ReadFullTable) THEN
        ! for two separate tables
        ! path of the original compose table (.thermo and .compo for now, will add hdf5 option later)
        CompOSEFilePath = 'SFHo/'
        CompOSEFHDF5Path = 'SFHo/eoscompose.h5'
        IF (RedHDF5Table) THEN
            BaryonEOSTableName = 'FullEOS_interpolated.h5'
        ELSE
            BaryonEOSTableName = 'FullEOS.h5'
        ENDIF
        HelmEOSTableName = 'HelmEOS.h5'
        MuonEOSTableName = 'MuonEOS.h5'
        iLepton = 1
    ELSE 
        ! for one big table
        ! path of the original compose table (.thermo and .compo for now, will add hdf5 option later)
        CompOSEFilePath = 'SFHo_no_ele/'
        CompOSEFHDF5Path = 'SFHo_no_ele/eoscompose.h5'
        IF (RedHDF5Table) THEN
            BaryonEOSTableName = 'BaryonsPlusHelmPlusMuonsEOS_interpolated.h5'
        ELSE
            BaryonEOSTableName = 'BaryonsPlusHelmPlusMuonsEOS.h5'
        ENDIF
        HelmEOSTableName = BaryonEOSTableName
        MuonEOSTableName = BaryonEOSTableName
        iLepton = 0
    ENDIF 
    
    nVariables = 15

    IF (RedHDF5Table) THEN
        ! This is to read the HDF5 file
        CALL ReadCompOSEHDFTable( CompOSEFHDF5Path, nPointsBaryon(1), nPointsBaryon(2), nPointsBaryon(3), iLepton )
    ELSE 
        ! This is to read the text .thermo and .compo files
        CALL ReadnPointsFromCompOSE(CompOSEFilePath, nPointsBaryon(1), nPointsBaryon(2), nPointsBaryon(3))    
        ! Fill the rho, temperature, ye, and compose EOS table arrays
        CALL ReadCompOSETable( CompOSEFilePath, nPointsBaryon(1), nPointsBaryon(2), nPointsBaryon(3) )
    ENDIF 

    ! -------------------- NOW DO BARYONIC EOS ----------------------------------------------
    PRINT*, "Allocate Baryonic EOS"
    WRITE(*,*) nPointsBaryon
    CALL AllocateEquationOfStateTable( EOSTable, nPointsBaryon , nVariables )
    
    EOSTable % TS % Names(1:3) = (/'Density                         ',&
    'Temperature                     ',&
    'Proton Fraction                 '/)
    
    EOSTable % TS % Indices % iRho = 1
    EOSTable % TS % Indices % iT   = 2
    EOSTable % TS % Indices % iYe  = 3
    
    PRINT*, "Allocate Independent Variable Units " 
    
    EOSTable % TS % Units(1:3) = (/'Grams per cm^3                  ', &
    'K                               ', &
    '                                '/) 
    
    ! DO WE NEED THIS?
    EOSTable % TS % minValues(1:3) =  (/MINVAL(RhoCompOSE), MINVAL(TempCompOSE), MINVAL(YpCompOSE)/)
    EOSTable % TS % maxValues(1:3) =  (/MAXVAL(RhoCompOSE), MAXVAL(TempCompOSE), MAXVAL(YpCompOSE)/)
    
    PRINT*, "Label Grid Type"
    EOSTable % TS % LogInterp(1:3) =  (/1, 1, 0/)
    
    PRINT*, "Allocate Names " 
    EOSTable % DV % Names(1:15) = (/'Pressure                        ', &
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
    'Sound Speed                     ', &
    'Gamma1                          '/)
    
    PRINT*, "Set Dependent Variable Identifier Indicies "
    
    EOSTable % DV % Indices % iPressure = 1
    EOSTable % DV % Indices % iEntropyPerBaryon = 2
    EOSTable % DV % Indices % iInternalEnergyDensity = 3
    EOSTable % DV % Indices % iElectronChemicalPotential = 4
    EOSTable % DV % Indices % iProtonChemicalPotential = 5
    EOSTable % DV % Indices % iNeutronChemicalPotential = 6
    EOSTable % DV % Indices % iProtonMassFraction = 7
    EOSTable % DV % Indices % iNeutronMassFraction = 8
    EOSTable % DV % Indices % iAlphaMassFraction = 9
    EOSTable % DV % Indices % iHeavyMassFraction = 10
    EOSTable % DV % Indices % iHeavyChargeNumber = 11
    EOSTable % DV % Indices % iHeavyMassNumber = 12
    EOSTable % DV % Indices % iHeavyBindingEnergy = 13
    EOSTable % DV % Indices % iThermalEnergy = 14
    EOSTable % DV % Indices % iGamma1 = 15
    
    PRINT*, "Allocate Dependent Variable Units " 
    EOSTable % DV % Units(1:15) = (/'Dynes per cm^2                  ', &
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
    'cm per s                        ', &
    '                                '/)
    
    PRINT*, "Begin Populating EOSTable" 
    ! FROM HERE ON YOU HAVE TO CONVERT THE COMPOSE EOS TO WEAKLIB
    EOSTable % TS % States(1) % Values = RhoCompOSE
    EOSTable % TS % States(2) % Values = TempCompOSE
    EOSTable % TS % States(3) % Values = YpCompOSE
    
    ! thermo and compo state
    EOSTable % DV % Variables(1) % Values(:,:,:) = EOSCompOSE(:,:,:,1)
    EOSTable % DV % Variables(2) % Values(:,:,:) = EOSCompOSE(:,:,:,2)
    EOSTable % DV % Variables(3) % Values(:,:,:) = EOSCompOSE(:,:,:,3)
    EOSTable % DV % Variables(4) % Values(:,:,:) = EOSCompOSE(:,:,:,4)
    EOSTable % DV % Variables(5) % Values(:,:,:) = EOSCompOSE(:,:,:,5)
    EOSTable % DV % Variables(6) % Values(:,:,:) = EOSCompOSE(:,:,:,6)
    EOSTable % DV % Variables(7) % Values(:,:,:) = EOSCompOSE(:,:,:,7)
    EOSTable % DV % Variables(8) % Values(:,:,:) = EOSCompOSE(:,:,:,8)
    EOSTable % DV % Variables(9) % Values(:,:,:) = EOSCompOSE(:,:,:,9)
    EOSTable % DV % Variables(10) % Values(:,:,:) = EOSCompOSE(:,:,:,10)
    EOSTable % DV % Variables(11) % Values(:,:,:) = EOSCompOSE(:,:,:,11)
    EOSTable % DV % Variables(12) % Values(:,:,:) = EOSCompOSE(:,:,:,12)
    EOSTable % DV % Variables(13) % Values(:,:,:) = EOSCompOSE(:,:,:,13)
    EOSTable % DV % Variables(14) % Values(:,:,:) = EOSCompOSE(:,:,:,14)
    EOSTable % DV % Variables(15) % Values(:,:,:) = EOSCompOSE(:,:,:,15)
    
    DEALLOCATE(RhoCompOSE)
    DEALLOCATE(TempCompOSE)
    DEALLOCATE(YpCompOSE)
    DEALLOCATE(EOSCompOSE)
    
    PRINT*, "Offset negative quantities"
    
    IF (ResetNegativePressure) THEN
        DO iRho=1,nPointsBaryon(1)
         DO iTemp=1,nPointsBaryon(2)
          DO iYe=1,nPointsBaryon(3)
            IF (EOSTable % DV % Variables(1) % Values(iRho,iTemp,iYe) .lt. 0.0d0) THEN
                EOSTable % DV % Variables(1) % Values(iRho,iTemp,iYe) = 1.0d0
            ENDIF
          ENDDO
         ENDDO
        ENDDO
    ENDIF
    
    EOSTable % DV % Offsets(:) = 0.0_dp
    
    DO iVars = 2, EOSTable % DV % Indices % iNeutronChemicalPotential
       Minimum_Value = MINVAL(EOSTable % DV % Variables(iVars) % Values)
       IF ( Minimum_Value .eq. 0.0_dp) THEN
           EOSTable % DV % Offsets(iVars) = 1.0d-100
           EOSTable % DV % Variables(iVars) % Values =  &
               EOSTable % DV % Variables(iVars) % Values + & 
               EOSTable % DV % Offsets(iVars)
       ELSE IF ( Minimum_Value .lt. 0.0_dp) THEN
           EOSTable % DV % Offsets(iVars) = -1.1_dp * Minimum_Value
           EOSTable % DV % Variables(iVars) % Values =  &
               EOSTable % DV % Variables(iVars) % Values + & 
               EOSTable % DV % Offsets(iVars)
       ELSE
           EOSTable % DV % Offsets(iVars) = 0.0_dp
       END IF
       
       WRITE(*,*) iVars, EOSTable % DV % Offsets(iVars), Minimum_Value
   END DO
   
   !EOSTable % DV % Variables(EOSTable % DV % Indices % iPressure) % Values = &
   !    LOG10(EOSTable % DV % Variables(EOSTable % DV % Indices % iPressure) % Values)
   EOSTable % DV % Variables(EOSTable % DV % Indices % iEntropyPerBaryon) % Values = &
       LOG10(EOSTable % DV % Variables(EOSTable % DV % Indices % iEntropyPerBaryon) % Values)
   EOSTable % DV % Variables(EOSTable % DV % Indices % iInternalEnergyDensity) % Values = &
       LOG10(EOSTable % DV % Variables(EOSTable % DV % Indices % iInternalEnergyDensity) % Values)
       
   EOSTable % DV % Variables(EOSTable % DV % Indices % iProtonChemicalPotential) % Values = &
       LOG10(EOSTable % DV % Variables(EOSTable % DV % Indices % iProtonChemicalPotential) % Values)
   EOSTable % DV % Variables(EOSTable % DV % Indices % iNeutronChemicalPotential) % Values = &
       LOG10(EOSTable % DV % Variables(EOSTable % DV % Indices % iNeutronChemicalPotential) % Values)
                
    ! ------------- NOW DO ELECTRON EOS ------------------ !
    nPointsHelm = (/ iTempMax, iDenMax /)
    PRINT*, "Allocate Helmholtz EOS"
    CALL AllocateHelmholtzTable( HelmEOSTable, nPointsHelm )
    
    HelmDatFilePath = '../helm_table.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmEOSTable )
    
    ! ------------- NOW DO MUON EOS ------------------ !
    nPointsDenon = (/ nTempMuon, nDenMuon /)
    PRINT*, "Allocate Muon EOS"
    CALL AllocateMuonEOS( MuonEOSTable, nPointsDenon )
    
    MuonDatFilePath = '../muons_fixedrho.dat'
    CALL ReadMuonEOSdat( MuonDatFilePath, MuonEOSTable )
    
    ! NOW CREATE BARYONIC FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Baryonic EOS"
    CALL WriteEquationOfStateTableHDF( EOSTable, BaryonEOSTableName )
    CALL FinalizeHDF( )
    
    ! NOW ADD Helmholtz EOS TO PREVIOUSLY CREATED H5 FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Appending Helmholtz EOS to HDF file"
    CALL WriteHelmholtzTableHDF( HelmEOSTable, HelmEOSTableName, ReadFullTable )
    CALL FinalizeHDF( )
    
    ! NOW ADD Muon EOS TO PREVIOUSLY CREATED H5 FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Appending Muon EOS to HDF file"
    CALL WriteMuonTableHDF( MuonEOSTable, MuonEOSTableName, ReadFullTable )
    CALL FinalizeHDF( )
    
    WRITE (*,*) "HDF write successful"
    
    CALL DeAllocateEquationOfStateTable( EOSTable )
    CALL DeallocateHelmholtzTable( HelmEOSTable )
    CALL DeAllocateMuonEOS( MuonEOSTable )

    ! Now Create Electron EOS

END PROGRAM wlCreateEquationOfStateTable

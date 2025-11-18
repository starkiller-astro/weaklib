PROGRAM wlCreateEquationOfStateTable
    
    USE wlKindModule, ONLY: dp
    USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
    USE HDF5
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
    USE wlEOSIOModuleHDF
    USE wlCompOSEInterface, ONLY : ReadnPointsFromCompOSE, ReadCompOSETable, &
        ReadCompOSEHDFTable, RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE
    USE wlHelmIOModuleHDF, ONLY : WriteHelmholtzTableHDF
    USE wlLeptonEOSTableModule, ONLY: &
        HelmTableType, &
        iTempMax, iDenMax, &
        ReadHelmEOSdat, &
        AllocateHelmholtzTable, DeallocateHelmholtzTable
    USE wlEosConstantsModule, ONLY: &
      cvel, ergmev, cm3fm3, kmev_inv, rmu, mn, me, mp, mmu

    IMPLICIT NONE
    
    ! ! PARAMETERS OF MUON TABLE, NEED TO KNOW THIS IN ADVANCE
    ! INTEGER,  PARAMETER :: nTempMuons = 101
    ! REAL(DP), PARAMETER :: tlo_Muons  = 8.0e+00_dp
    ! REAL(DP), PARAMETER :: thi_Muons  = 1.3e+01_dp

    ! INTEGER,  PARAMETER :: nDenMuons  = 421
    ! REAL(DP), PARAMETER :: dlo_Muons  = -5.05e+00_dp
    ! REAL(DP), PARAMETER :: dhi_Muons  = 15.95e+00_dp

    ! PARAMETERS OF ELECTRON TABLE, NEED TO KNOW THIS IN ADVANCE
    ! THIS IS FOR THE FIRST MUON TABLE WE MADE
    INTEGER,  PARAMETER :: nTempMuons = 270
    REAL(DP), PARAMETER :: tlo_Muons  =  9.763597279797683_dp
    REAL(DP), PARAMETER :: thi_Muons  = 12.45359726262338_dp

    INTEGER, PARAMETER  :: nDenMuons  = 1501
    REAL(DP), PARAMETER :: dlo_Muons  =  3.220259248823259_dp
    REAL(DP), PARAMETER :: dhi_Muons  = 15.85594401137048_dp

    ! PARAMETERS OF ELECTRON TABLE, NEED TO KNOW THIS IN ADVANCE
    INTEGER,  PARAMETER :: nTempElectrons = 201
    REAL(DP), PARAMETER :: tlo_Electrons  = 3.0_dp
    REAL(DP), PARAMETER :: thi_Electrons  = 13.0_dp 

    INTEGER,  PARAMETER :: nDenElectrons  = 541
    REAL(DP), PARAMETER :: dlo_Electrons  = -12.0_dp
    REAL(DP), PARAMETER :: dhi_Electrons  = 15.0_dp

    INTEGER                        :: iVars
    INTEGER, DIMENSION(3)          :: nPointsBaryon
    INTEGER, DIMENSION(2)          :: nPointsHelm, nPointsDenon
    INTEGER                        :: nVariables
    TYPE(EquationOfStateCompOSETableType) :: EOSTable
    TYPE(HelmTableType) :: HelmTableElectrons
    TYPE(HelmTableType) :: HelmTableMuons
    
    CHARACTER(len=128) :: CompOSEFilePath, CompOSEFHDF5Path, BaryonEOSTableName, &
        HelmDatFilePath, HelmEOSTableName, MuonEOSTableName
    
    REAL(dp) :: Minimum_Value, Add_to_energy
    LOGICAL  :: ReadFullTable, RedHDF5Table, ResetNegativePressure
    INTEGER  :: iLepton, iRho, iTemp, iYe

    ReadFullTable = .false.
    RedHDF5Table  = .false.
    ResetNegativePressure = .false.
    nVariables = 19
    
    Add_to_energy = 8.9d0*ergmev/rmu
    Add_to_energy = 2.0d0*ergmev/rmu
    Add_to_energy = + cvel**2 - ergmev * mn / rmu + 8.9d0 * ergmev / rmu
    Add_to_energy = 0.0d0

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
            BaryonEOSTableName = 'BaryonsPlusPhotonsPlusLeptonsEOS_interpolated.h5'
        ELSE
            BaryonEOSTableName = 'BaryonsPlusPhotonsPlusLeptonsEOS.h5'
        ENDIF
        HelmEOSTableName = BaryonEOSTableName
        MuonEOSTableName = BaryonEOSTableName
        iLepton = 0
    ENDIF 

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
    EOSTable % DV % Names(1:19) = (/'Pressure                        ', &
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
    'Thermal Energy                  ', &
    'Heavy Binding Energy            ', &
    'Gamma1                          ', &
    'Neutron Effective Mass          ', &
    'Proton Effective Mass           ', &
    'Neutron Vector Self Energy      ', &
    'Proton Vector Self Energy       '/)
    
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
    EOSTable % DV % Indices % iNeutronEffMass = 16
    EOSTable % DV % Indices % iProtonEffMass = 17
    EOSTable % DV % Indices % iNeutronSelfEnergy = 18
    EOSTable % DV % Indices % iProtonSelfEnergy = 19
    
    PRINT*, "Allocate Dependent Variable Units " 
    EOSTable % DV % Units(1:19) = (/'Dynes per cm^2                  ', &
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
    '                                ', &
    'MeV                             ', &
    'MeV                             ', &
    'MeV                             ', &
    'MeV                             '/)
    
    PRINT*, "Begin Populating EOSTable" 
    ! FROM HERE ON YOU HAVE TO CONVERT THE COMPOSE EOS TO WEAKLIB
    EOSTable % TS % States(1) % Values = RhoCompOSE
    EOSTable % TS % States(2) % Values = TempCompOSE
    EOSTable % TS % States(3) % Values = YpCompOSE
    
    ! thermo and compo state
    EOSTable % DV % Variables(1)  % Values(:,:,:) = EOSCompOSE(:,:,:,1)
    EOSTable % DV % Variables(2)  % Values(:,:,:) = EOSCompOSE(:,:,:,2)
    EOSTable % DV % Variables(3)  % Values(:,:,:) = EOSCompOSE(:,:,:,3) + Add_to_energy
    EOSTable % DV % Variables(4)  % Values(:,:,:) = EOSCompOSE(:,:,:,4)
    EOSTable % DV % Variables(5)  % Values(:,:,:) = EOSCompOSE(:,:,:,5)
    EOSTable % DV % Variables(6)  % Values(:,:,:) = EOSCompOSE(:,:,:,6)
    EOSTable % DV % Variables(7)  % Values(:,:,:) = EOSCompOSE(:,:,:,7)
    EOSTable % DV % Variables(8)  % Values(:,:,:) = EOSCompOSE(:,:,:,8)
    EOSTable % DV % Variables(9)  % Values(:,:,:) = EOSCompOSE(:,:,:,9)
    EOSTable % DV % Variables(10) % Values(:,:,:) = EOSCompOSE(:,:,:,10)
    EOSTable % DV % Variables(11) % Values(:,:,:) = EOSCompOSE(:,:,:,11)
    EOSTable % DV % Variables(12) % Values(:,:,:) = EOSCompOSE(:,:,:,12)
    EOSTable % DV % Variables(13) % Values(:,:,:) = EOSCompOSE(:,:,:,13)
    EOSTable % DV % Variables(14) % Values(:,:,:) = EOSCompOSE(:,:,:,14)
    EOSTable % DV % Variables(15) % Values(:,:,:) = EOSCompOSE(:,:,:,15)
    EOSTable % DV % Variables(16) % Values(:,:,:) = EOSCompOSE(:,:,:,16)
    EOSTable % DV % Variables(17) % Values(:,:,:) = EOSCompOSE(:,:,:,17)
    EOSTable % DV % Variables(18) % Values(:,:,:) = EOSCompOSE(:,:,:,18)
    EOSTable % DV % Variables(19) % Values(:,:,:) = EOSCompOSE(:,:,:,19)
    
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
    
    WRITE(*,*) 'These are the offsets'
    DO iVars = 2, nVariables
        Minimum_Value = MINVAL(EOSTable % DV % Variables(iVars) % Values)
        IF ( Minimum_Value .lt. 0.0_dp) THEN
            EOSTable % DV % Offsets(iVars) = -1.1_dp * Minimum_Value
            EOSTable % DV % Variables(iVars) % Values =  &
                EOSTable % DV % Variables(iVars) % Values + & 
                EOSTable % DV % Offsets(iVars)
        ELSE IF ( Minimum_Value .gt. 0.0_dp) THEN
            EOSTable % DV % Offsets(iVars) = 0.0_dp
        ELSE
          WRITE(*,*) 'Minimum is zero', iVars, Minimum_Value
        ENDIF

        WRITE(*,*) iVars, EOSTable % DV % Offsets(iVars), Minimum_Value

    END DO
   
    ! Now log everything before putting it in the table (except for pressure)
    DO iVars = 2, nVariables
      EOSTable % DV % Variables(iVars) % Values = &
          LOG10( MAX(EOSTable % DV % Variables(iVars) % Values, 1.0d-99) )
    END DO
                  
    ! ------------- NOW DO ELECTRON EOS ------------------ !
    nPointsHelm = (/ nDenElectrons, nTempElectrons /)
    PRINT*, "Allocate Helmholtz EOS for Electrons"
    CALL AllocateHelmholtzTable( HelmTableElectrons, nPointsHelm )

    WRITE(*,*) 'HelmTableElectrons Dimensions = ', &
               HelmTableElectrons % nPointsDen, HelmTableElectrons % nPointsTemp
    
    HelmDatFilePath = '../electron_table_p256_q800.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmTableElectrons, me, &
       tlo_Electrons, thi_Electrons, dlo_Electrons, dhi_Electrons )
    
    ! ------------- NOW DO MUON EOS ------------------ !
    nPointsHelm = (/ nDenMuons, nTempMuons /)
    PRINT*, "Allocate Helmholtz EOS for Muons"
    CALL AllocateHelmholtzTable( HelmTableMuons, nPointsHelm )
    
    WRITE(*,*) 'HelmTableMuons Dimensions = ', &
               HelmTableMuons % nPointsDen, HelmTableMuons % nPointsTemp
    
    HelmDatFilePath = '/home/lbocciol/MuonProject/muon_table_p256_q800.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmTableMuons, mmu, &
       tlo_Muons, thi_Muons, dlo_Muons, dhi_Muons)
    
    ! NOW CREATE BARYONIC FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Baryonic EOS"
    CALL WriteEquationOfStateTableHDF( EOSTable, BaryonEOSTableName )
    CALL FinalizeHDF( )
    
    ! NOW ADD Helmholtz EOS TO PREVIOUSLY CREATED H5 FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Appending Helmholtz EOS for Electrons to HDF file"
    CALL WriteHelmholtzTableHDF( HelmTableElectrons, HelmEOSTableName, "HelmTableElectrons", ReadFullTable )
    CALL FinalizeHDF( )
    
    ! NOW ADD Muon EOS TO PREVIOUSLY CREATED H5 FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Appending Helmholtz EOS for Muons to HDF file"
    CALL WriteHelmholtzTableHDF( HelmTableMuons, MuonEOSTableName, "HelmTableMuons", ReadFullTable )
    CALL FinalizeHDF( )
    
    WRITE (*,*) "HDF write successful"
    
    CALL DeAllocateEquationOfStateTable( EOSTable )
    CALL DeallocateHelmholtzTable( HelmTableElectrons )
    CALL DeallocateHelmholtzTable( HelmTableMuons )

    ! Now Create Electron EOS

END PROGRAM wlCreateEquationOfStateTable

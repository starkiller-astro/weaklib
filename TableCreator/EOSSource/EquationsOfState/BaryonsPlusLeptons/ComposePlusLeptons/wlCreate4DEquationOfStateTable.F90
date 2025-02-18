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
    USE wlEosConstantsModule, ONLY: cvel, ergmev, cm3fm3, kmev_inv, rmu, mn, me, mp

    IMPLICIT NONE
    
    INTEGER                          :: iVars
    INTEGER, DIMENSION(3)            :: nPointsCompose
    INTEGER, DIMENSION(4)            :: nPoints
    INTEGER                          :: nVariables
    TYPE(EquationOfState4DTableType) :: EOSTable
    TYPE(HelmholtzTableType)         :: HelmEOSTable
    TYPE(MuonEOSType)                :: MuonEOSTable
    
    CHARACTER(len=128) :: CompOSEFilePath, CompOSEFHDF5Path, &
        HelmDatFilePath, MuonDatFilePath, &
        BaryonEOSTableName, HelmEOSTableName, MuonEOSTableName
    
    REAL(dp) :: Minimum_Value, Add_to_energy
    LOGICAL  :: RedHDF5Table, ResetNegativePressure
    INTEGER  :: iLepton, iRho, iTemp, iYe

    REAL(DP), allocatable, DIMENSION(1) :: YmGrid

    RedHDF5Table = .true.
    ResetNegativePressure = .false.
    
    Add_to_energy = 8.9d0*ergmev/rmu
    Add_to_energy = 2.0d0*ergmev/rmu
    Add_to_energy = + cvel**2 - ergmev * mn / rmu + 8.9d0 * ergmev / rmu

    ! path of the original compose table
    CompOSEFilePath = 'SFHo_no_ele/'
    CompOSEFHDF5Path = 'SFHo_no_ele/eoscompose.h5'
    IF (RedHDF5Table) THEN
        BaryonEOSTableName = 'BaryonsPlusHelmPlusMuonsEOS_interpolated.h5'
    ELSE
        BaryonEOSTableName = 'BaryonsPlusHelmPlusMuonsEOS.h5'
    ENDIF
    iLepton = 0
    
    nVariables = 17

    IF (RedHDF5Table) THEN
        ! This is to read the HDF5 file
        CALL ReadCompOSEHDFTable( CompOSEFHDF5Path, nPointsCompose(1), &
            nPointsCompose(2), nPointsCompose(3), iLepton )
    ELSE 
        ! This is to read the text .thermo and .compo files
        CALL ReadnPointsFromCompOSE(CompOSEFilePath, nPointsCompose(1), &
            nPointsCompose(2), nPointsCompose(3))    
        ! Fill the rho, temperature, ye, and compose EOS table arrays
        CALL ReadCompOSETable( CompOSEFilePath, nPointsCompose(1), &
            nPointsCompose(2), nPointsCompose(3) )
    ENDIF

    ! At this point you have filled the Compose EOS, now read in the Helmholtz and Muon EOS
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

    ! Set up desired grid
    nPoints(1) = nPointsCompose(1)
    nPoints(2) = nPointsCompose(2)
    nPoints(3) = nPointsCompose(3)
    nPoints(4) = 30

    ALLOCATE( YmGrid(nPoints(4)) )

    ! Set up logarithmic grid, we will see exactly how
    YmGrid = 0.0

    ! -------------------- NOW DO BARYONIC EOS ----------------------------------------------
    PRINT*, "Allocate Baryonic EOS"
    WRITE(*,*) nPoints
    CALL AllocateEquationOfState4DTable( EOSTable, nPoints , nVariables )
    
    EOSTable % TS % Names(1:4) = (/'Density            ',&
                                   'Temperature        ',&
                                   'Electron Fraction  ',&
                                   'Muon Fraction      '/)
    
    EOSTable % TS % Indices % iRho = 1
    EOSTable % TS % Indices % iT   = 2
    EOSTable % TS % Indices % iYe  = 3
    EOSTable % TS % Indices % iYm  = 4
    
    PRINT*, "Allocate Independent Variable Units " 
    
    EOSTable % TS % Units(1:4) = (/'Grams per cm^3                  ', &
    'K                               ', &
    '                                ',&
    '                                '/) 
    
    ! DO WE NEED THIS?
    EOSTable % TS % minValues(1:4) =  &
        (/MINVAL(RhoCompOSE), MINVAL(TempCompOSE), MINVAL(YpCompOSE), MINVAL(YmGrid)/)
    EOSTable % TS % maxValues(1:4) =  &
        (/MAXVAL(RhoCompOSE), MAXVAL(TempCompOSE), MAXVAL(YpCompOSE), MAXVAL(YmGrid)/)
    
    PRINT*, "Label Grid Type"
    EOSTable % TS % LogInterp(1:4) =  (/1, 1, 0, 1/)
    
    PRINT*, "Allocate Names " 
    EOSTable % DV % Names(1:17) = (/'Pressure                        ', &
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
    EOSTable % DV % Indices % iNeutronSelfEnergy = 16
    EOSTable % DV % Indices % iProtonSelfEnergy = 17
    
    PRINT*, "Allocate Dependent Variable Units " 
    EOSTable % DV % Units(1:17) = (/'Dynes per cm^2                  ', &
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
    'MeV                             '/)
    
    PRINT*, "Begin Populating EOSTable" 
    
    ! Now the real stuff
    ! This will look like 4 loops where you have to interpolate the Compose table onto the new grid
    ! By now having Ye + Ym = Yp. You only need to do 1D interpolations though (along Yp). 
    ! A more sophisticated interpolation scheme might be desirable. You will encounter issues 
    ! Where MAX(Ye) + MAX(Ym) will be above the Compose table. We should decide how to handle this.
    ! A straightforward way would be to do somthing like 
    ! IF (Yp > MAXVAL(YpCompOSE)) THEN
    !   Yp = Ye
    ! ELSE
    !   Yp = Ye + Ym
    ! ENDIF
    ! This is reasonable since MAX(Ym) ~ 0.1 and you only have muons inside the PNS where Ye < 0.4
    ! and typically MAXVAL(YpCompOSE) ~ 0.55 so you should be quite safe.



    ! NOW CREATE BARYONIC FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Baryonic EOS"
    CALL WriteEquationOfStateTableHDF( EOSTable, BaryonEOSTableName )
    CALL FinalizeHDF( )

    WRITE (*,*) "HDF write successful"
    
    CALL DeAllocateEquationOfStateTable( EOSTable )
    CALL DeallocateHelmholtzTable( HelmEOSTable )
    CALL DeAllocateMuonEOS( MuonEOSTable )

    ! Now Create Electron EOS

END PROGRAM wlCreateEquationOfStateTable

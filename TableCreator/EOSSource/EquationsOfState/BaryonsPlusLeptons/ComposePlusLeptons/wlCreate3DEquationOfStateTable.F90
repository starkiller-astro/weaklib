PROGRAM wlCreateEquationOfStateTable
    
    USE wlKindModule, ONLY: dp
    USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
    USE HDF5
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
    USE wlEOSIOModuleHDF
    USE wlInterpolationUtilitiesModule, ONLY: &
      LinearInterp_Array_Point, &
      GetIndexAndDelta_Lin, GetIndexAndDelta_Log
    USE wlCompOSEInterface, ONLY : &
      ReadnPointsFromCompOSE, ReadCompOSETable, &
      ReadCompOSEHDFTable, RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE
    USE wlHelmIOModuleHDF, ONLY : &
      WriteHelmholtzTableHDF
    USE wlLeptonEOSTableModule, ONLY: &
      HelmTableType, &
      iTempMax, iDenMax, &
      ReadHelmEOSdat, &
      AllocateHelmholtzTable, DeallocateHelmholtzTable
    USE wlLeptonPhotonGasEOS, ONLY: &
      PhotonGasType, LeptonGasType, &
      PhotonGasEOS , LeptonGasEOS , &
      GetPhotonLeptonGasEOS
    USE wlEosConstantsModule, ONLY: &
      cvel, ergmev, cm3fm3, kmev_inv, rmu, mn, me, mp, mmu
    USE wlSoundSpeedModule

    IMPLICIT NONE
    
    ! PARAMETERS OF MUON TABLE, NEED TO KNOW THIS IN ADVANCE
    INTEGER, PARAMETER  :: nTempMuons = 270
    REAL(DP), PARAMETER :: tlo_Muons  =  9.763597279797683_dp
    REAL(DP), PARAMETER :: thi_Muons  = 12.45359726262338_dp

    INTEGER, PARAMETER  :: nDenMuons  = 1501
    REAL(DP), PARAMETER :: dlo_Muons  =  3.220259248823259_dp
    REAL(DP), PARAMETER :: dhi_Muons  = 15.85594401137048_dp

    ! PARAMETERS OF ELECTRON TABLE, NEED TO KNOW THIS IN ADVANCE
    INTEGER, PARAMETER  :: nTempElectrons = 541
    REAL(DP), PARAMETER :: tlo_Electrons  = 3.0_dp
    REAL(DP), PARAMETER :: thi_Electrons  = 13.0_dp 

    INTEGER, PARAMETER  :: nDenElectrons  = 201
    REAL(DP), PARAMETER :: dlo_Electrons  = -12.0_dp
    REAL(DP), PARAMETER :: dhi_Electrons  = 15.0_dp
        
    INTEGER                          :: iVars
    INTEGER, DIMENSION(2)            :: nPoints2D
    INTEGER, DIMENSION(3)            :: nPointsCompose
    INTEGER, DIMENSION(3)            :: nPoints
    INTEGER                          :: nVariables
    TYPE(EquationOfStateCompOSETableType) :: EOSTable
    TYPE(HelmTableType)              :: HelmTableElectrons
    TYPE(HelmTableType)              :: HelmTableMuons
    TYPE(PhotonGasType)              :: PhotonGasState
    TYPE(LeptonGasType)              :: ElectronGasState
    TYPE(LeptonGasType)              :: MuonGasState
    
    CHARACTER(len=128) :: CompOSEFilePath, CompOSEFHDF5Path, &
        HelmDatFilePath, EOSTableName
    
    REAL(dp) :: Minimum_Value, Add_to_energy
    LOGICAL  :: RedHDF5Table
    INTEGER  :: iLepton, i, iRho, iTemp, iYe, iCount

    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: LogEnergy, LogEntropy

    REAL(DP) :: D, T, Yp, Ye
    REAL(DP) :: Xbary, dYp, OS_P, OS_E, OS_S, LocalOffset, Gamma, Cs

    RedHDF5Table = .false.
    
    Add_to_energy = 8.9d0*ergmev/rmu
    Add_to_energy = 2.0d0*ergmev/rmu
    Add_to_energy = + cvel**2 - ergmev * mn / rmu + 8.9d0 * ergmev / rmu
    Add_to_energy = 0.0d0

    ! path of the original compose table
    CompOSEFilePath = 'SFHo_no_ele/'
    CompOSEFHDF5Path = 'SFHo_no_ele/eoscompose.h5'
    IF (RedHDF5Table) THEN
        EOSTableName = '3DEOSTable_Highres.h5'
    ELSE
        EOSTableName = '3DEOSTable.h5'
    ENDIF
    iLepton = 0
    
    nVariables = 19

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
    nPoints2D = (/ nTempElectrons, nDenElectrons /)
    PRINT*, "Allocate Helmholtz EOS for Electrons"
    CALL AllocateHelmholtzTable( HelmTableElectrons, nPoints2D )
    
    HelmDatFilePath = '../electron_table_p256_q800.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmTableElectrons, me, &
       tlo_Electrons, thi_Electrons, dlo_Electrons, dhi_Electrons )
    
    ! ------------- NOW DO MUON EOS ------------------ !
    nPoints2D = (/ nTempMuons, nDenMuons /)
    PRINT*, "Allocate Helmholtz EOS for Muons"
    CALL AllocateHelmholtzTable( HelmTableMuons, nPoints2D )
    
    HelmDatFilePath = '../muon_table_p256_q800.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmTableMuons, mmu, &
       tlo_Muons, thi_Muons, dlo_Muons, dhi_Muons )

    ! Set up desired grid
    nPoints(1) = nPointsCompose(1)
    nPoints(2) = nPointsCompose(2)
    nPoints(3) = nPointsCompose(3)

    ! -------------------- NOW DO BARYONIC EOS ----------------------------------------------
    PRINT*, "Allocate Baryonic EOS"
    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )
             
    EOSTable % TS % States(1) % Values = RhoCompOSE
    EOSTable % TS % States(2) % Values = TempCompOSE
    EOSTable % TS % States(3) % Values = YpCompOSE ! This will actually become Ye!!!!

    EOSTable % TS % Names(1:3) = (/'Density            ',&
                                   'Temperature        ',&
                                   'Electron Fraction  '/)
    
    EOSTable % TS % Indices % iRho = 1
    EOSTable % TS % Indices % iT   = 2
    EOSTable % TS % Indices % iYe  = 3
    
    PRINT*, "Allocate Independent Variable Units " 
    
    EOSTable % TS % Units(1:3) = (/'Grams per cm^3                  ', &
    'K                               ', &
    '                                '/) 
    
    ! DO WE NEED THIS?
    EOSTable % TS % minValues(1:3) =  &
        (/MINVAL(RhoCompOSE), MINVAL(TempCompOSE), MINVAL(YpCompOSE)/)
    EOSTable % TS % maxValues(1:3) =  &
        (/MAXVAL(RhoCompOSE), MAXVAL(TempCompOSE), MAXVAL(YpCompOSE)/)
    
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
    

    iCount = 0
    !$OMP PARALLEL DO PRIVATE(PhotonGasState, ElectronGasState, Ye) &
    !$OMP SHARED(iCount) &
    !$OMP COLLAPSE(3)
    DO iRho=1,nPoints(1)
      DO iTemp=1,nPoints(2)
        DO iYe=1,nPoints(3)

          ! Safely increment iCount without locking the whole loop
          !$OMP ATOMIC
          iCount = iCount + 1

          ! Print progress only every 1% of completion
          IF (MOD(iCount, INT(0.01d0 * nPoints(1) * nPoints(2) * nPoints(3))) .EQ. 0) THEN
            !$OMP CRITICAL
            WRITE(*,*) "Progress:", DBLE(iCount) / DBLE(nPoints(1) * nPoints(2) * nPoints(3)) * 100.0d0, "%"
            !$OMP END CRITICAL
          ENDIF

          Ye = EOSTable % TS % States(3) % Values(iYe)

          PhotonGasState % T   = EOSTable % TS % States(2) % Values(iTemp)
          PhotonGasState % rho = EOSTable % TS % States(1) % Values(iRho)
          CALL PhotonGasEOS(PhotonGasState)
              
          ElectronGasState % T    = EOSTable % TS % States(2) % Values(iTemp)
          ElectronGasState % rho  = EOSTable % TS % States(1) % Values(iRho)
          ElectronGasState % yL   = Ye

          CALL LeptonGasEOS(HelmTableElectrons, ElectronGasState)

          ! PRESSURE
          EOSTable % DV % Variables(1) % Values(iRho,iTemp,iYe) = &
          EOSCompOSE(iRho,iTemp,iYe,1) + ElectronGasState % p + PhotonGasState % p

          ! ENTROPY
          EOSTable % DV % Variables(2) % Values(iRho,iTemp,iYe) = &
          EOSCompOSE(iRho,iTemp,iYe,2) + ElectronGasState % s + PhotonGasState % s

          ! ENERGY
          EOSTable % DV % Variables(3) % Values(iRho,iTemp,iYe) = &
          EOSCompOSE(iRho,iTemp,iYe,3) + ElectronGasState % e + PhotonGasState % e + Add_to_energy

          ! ELECTRON CHEMICAL POTENTIAL
          EOSTable % DV % Variables(4) % Values(iRho,iTemp,iYe) = &
            ElectronGasState % mu

        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! thermo and compo state
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
    ! Gamma calculated later
    EOSTable % DV % Variables(16) % Values(:,:,:) = EOSCompOSE(:,:,:,16)
    EOSTable % DV % Variables(17) % Values(:,:,:) = EOSCompOSE(:,:,:,17)
    EOSTable % DV % Variables(18) % Values(:,:,:) = EOSCompOSE(:,:,:,18)
    EOSTable % DV % Variables(19) % Values(:,:,:) = EOSCompOSE(:,:,:,19)

    ! GAMMA
    OS_P = 0.0_dp

    OS_S = MINVAL( EOSCompOSE(:,:,:,2) )
    IF (OS_S .lt. 0.0_dp) THEN
      OS_S = -1.1d0*OS_S
    ELSE IF (OS_S .eq. 0.0_dp) THEN
      OS_S = 1.0d-10
    ELSE
      OS_S = 0.0_dp
    ENDIF

    OS_E = MINVAL( EOSCompOSE(:,:,:,3) )
    IF (OS_E .lt. 0.0_dp) THEN
      OS_E = -1.1d0*OS_E
    ELSE IF (OS_E .eq. 0.0_dp) THEN
      OS_E = 1.0d-10
    ELSE
      OS_E = 0.0_dp
    ENDIF

    ALLOCATE( LogEnergy(nPointsCompose(1), nPointsCompose(2), nPointsCompose(3)) )
    ALLOCATE( LogEntropy(nPointsCompose(1), nPointsCompose(2), nPointsCompose(3)) )

    LogEnergy  = LOG10(EOSCompOSE(:,:,:,3) + OS_E)
    LogEntropy = LOG10(EOSCompOSE(:,:,:,2) + OS_S)

    iCount = 0
    !$OMP PARALLEL DO PRIVATE(D, T, Ye, Gamma, Cs) &
    !$OMP SHARED(EOSTable, RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE, &
    !$OMP LogEnergy, LogEntropy, OS_P, OS_E, OS_S, HelmTableElectrons, HelmTableMuons,iCount) &
    !$OMP COLLAPSE(3)
    DO iRho=1,nPoints(1)
      DO iTemp=1,nPoints(2)
        DO iYe=1,nPoints(3)

            ! Safely increment iCount without locking the whole loop
            !$OMP ATOMIC
            iCount = iCount + 1

            ! Print progress only every 1% of completion
            IF (MOD(iCount, INT(0.01d0 * nPoints(1) * nPoints(2) * nPoints(3)  )) .EQ. 0) THEN
              !$OMP CRITICAL
              WRITE(*,*) "Progress:", DBLE(iCount) / DBLE( nPoints(1) * nPoints(2) * nPoints(3) ) * 100.0d0, "%"
              !$OMP END CRITICAL
            ENDIF

            D  = EOSTable % TS % States(1) % Values(iRho)
            T  = EOSTable % TS % States(2) % Values(iTemp)
            Ye = EOSTable % TS % States(3) % Values(iYe)
          
            CALL CalculateSoundSpeed( D, T, Ye, 0.0d0, RhoCompOSE(:), &
                  TempCompOSE(:), YpCompOSE(:), &
                  EOSCompOSE(:,:,:,1), OS_P, &
                  LogEnergy (:,:,:), OS_E, &
                  LogEntropy(:,:,:), OS_S, &
                  HelmTableElectrons, HelmTableMuons, Gamma, Cs, .FALSE.)

            EOSTable % DV % Variables(15) % Values(iRho,iTemp,iYe) = Gamma
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    WRITE(*,*) 'DONE'

    ! FINAL HOUSEKEEPING
    DEALLOCATE(LogEnergy)
    DEALLOCATE(LogEntropy)
    DEALLOCATE(RhoCompOSE)
    DEALLOCATE(TempCompOSE)
    DEALLOCATE(YpCompOSE)
    DEALLOCATE(EOSCompOSE)
    
    ! Calculate maxima and minima
    DO i = 1, EOSTable % DV % nVariables
      EOSTable % DV % maxValues(i) = MAXVAL( EOSTable % DV % Variables(i) % Values)
      EOSTable % DV % minValues(i) = MINVAL( EOSTable % DV % Variables(i) % Values)
    ENDDO

    PRINT*, "Offset negative quantities"

    EOSTable % DV % Offsets(:) = 0.0_dp
    
    WRITE(*,*) 'These are the offsets'
    DO iVars = 1, nVariables
        Minimum_Value = MINVAL(EOSTable % DV % Variables(iVars) % Values)
        IF ( Minimum_Value .lt. 0.0_dp) THEN
            EOSTable % DV % Offsets(iVars) = -1.1_dp * Minimum_Value
            EOSTable % DV % Variables(iVars) % Values =  &
                EOSTable % DV % Variables(iVars) % Values + & 
                EOSTable % DV % Offsets(iVars)
        ELSE IF ( Minimum_Value .gt. 0.0_dp) THEN
            EOSTable % DV % Offsets(iVars) = 0.0_dp
        ELSE
          WRITE(*,*) 'Minimum is zero', iVars
          EOSTable % DV % Variables(iVars) % Values = & 
              MAX(EOSTable % DV % Variables(iVars) % Values, 1.0e-99_dp)
        ENDIF

        WRITE(*,*) iVars, EOSTable % DV % Offsets(iVars), Minimum_Value

        EOSTable % DV % Variables(iVars) % Values = &
          LOG10(EOSTable % DV % Variables(iVars) % Values)

    END DO
   
    ! NOW CREATE BARYONIC FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: 3D EOS"
    CALL WriteEquationOfStateTableHDF( EOSTable, EOSTableName )
    CALL FinalizeHDF( )

    WRITE (*,*) "HDF write successful"
    
    CALL DeAllocateEquationOfStateTable( EOSTable )
    CALL DeallocateHelmholtzTable( HelmTableElectrons )
    CALL DeallocateHelmholtzTable( HelmTableMuons )

    ! Now Create Electron EOS

END PROGRAM wlCreateEquationOfStateTable

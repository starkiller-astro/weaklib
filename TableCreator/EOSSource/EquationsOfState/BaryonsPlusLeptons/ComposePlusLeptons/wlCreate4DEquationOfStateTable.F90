PROGRAM wlCreateEquationOfStateTable
    
    USE wlKindModule, ONLY: dp
    USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
    USE HDF5
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
    USE wlEOSIOModuleHDF
    USE wlInterpolationUtilitiesModule, ONLY: &
      LinearInterp_Array_Point, &
      GetIndexAndDelta_Lin
    USE wlCompOSEInterface, ONLY : ReadnPointsFromCompOSE, ReadCompOSETable, &
      ReadCompOSEHDFTable, RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE
    USE wlLeptonEOSTableModule, ONLY: &
      HelmTableType, &
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
    INTEGER, PARAMETER  :: nTempMuons = 101
    REAL(DP), PARAMETER :: tlo_Muons  = 8.0e+00_dp
    REAL(DP), PARAMETER :: thi_Muons  = 1.3e+01_dp

    INTEGER, PARAMETER  :: nDenMuons  = 421
    REAL(DP), PARAMETER :: dlo_Muons  = -5.05e+00_dp
    REAL(DP), PARAMETER :: dhi_Muons  = 15.95e+00_dp

    ! PARAMETERS OF ELECTRON TABLE, NEED TO KNOW THIS IN ADVANCE
    INTEGER, PARAMETER  :: nTempElectrons = 201
    REAL(DP), PARAMETER :: tlo_Electrons  = 3.0_dp
    REAL(DP), PARAMETER :: thi_Electrons  = 13.0_dp 

    INTEGER, PARAMETER  :: nDenElectrons  = 541
    REAL(DP), PARAMETER :: dlo_Electrons  = -12.0_dp
    REAL(DP), PARAMETER :: dhi_Electrons  = 15.0_dp

    INTEGER                          :: iVars
    INTEGER, DIMENSION(2)            :: nPointsHelm
    INTEGER, DIMENSION(3)            :: nPointsCompose
    INTEGER, DIMENSION(4)            :: nPoints
    INTEGER                          :: nVariables
    TYPE(EquationOfState4DTableType) :: EOSTable
    TYPE(HelmTableType)              :: HelmTableElectrons
    TYPE(HelmTableType)              :: HelmTableMuons
    TYPE(PhotonGasType)              :: PhotonGasState
    TYPE(LeptonGasType)              :: ElectronGasState
    TYPE(LeptonGasType)              :: MuonGasState
    
    CHARACTER(len=128) :: CompOSEFilePath, CompOSEFHDF5Path, &
        HelmDatFilePath, MuonDatFilePath, &
        BaryonEOSTableName
    
    REAL(dp) :: Minimum_Value, Add_to_energy
    LOGICAL  :: RedHDF5Table
    INTEGER  :: iLepton, i, iRho, iTemp, iYe, iYm, iYp, iCount

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: YmGrid
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: LogEnergy, LogEntropy

    REAL(DP) :: D, T, Yp, Ye, Ym
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
        BaryonEOSTableName = '4DEOSTableHighres.h5'
    ELSE
        BaryonEOSTableName = '4DEOSTable.h5'
    ENDIF
    iLepton = 0
    
    nVariables = 20

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
    nPointsHelm = (/ nDenElectrons, nTempElectrons /)
    PRINT*, "Allocate Helmholtz EOS for Electrons"
    CALL AllocateHelmholtzTable( HelmTableElectrons, nPointsHelm )
    
    HelmDatFilePath = '../electron_table_p256_q800.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmTableElectrons, me, &
       tlo_Electrons, thi_Electrons, dlo_Electrons, dhi_Electrons )
    
    ! ------------- NOW DO MUON EOS ------------------ !
    nPointsHelm = (/ nDenMuons, nTempMuons /)
    PRINT*, "Allocate Helmholtz EOS for Muons"
    CALL AllocateHelmholtzTable( HelmTableMuons, nPointsHelm )
    
    HelmDatFilePath = '../muon_table_p256_q800.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmTableMuons, mmu, &
       tlo_Muons, thi_Muons, dlo_Muons, dhi_Muons )

    ! Set up desired grid
    nPoints(1) = nPointsCompose(1)
    nPoints(2) = nPointsCompose(2)
    nPoints(3) = nPointsCompose(3)
    nPoints(4) = 35

    ALLOCATE( YmGrid(nPoints(4)) )

    ! Set up logarithmic grid, we will see exactly how
     CALL MakeLogGrid( 1.0d-6, 0.1_dp, nPoints(4), YmGrid )

    ! -------------------- NOW DO BARYONIC EOS ----------------------------------------------
    PRINT*, "Allocate Baryonic EOS"
    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )
             
    EOSTable % TS % States(1) % Values = RhoCompOSE
    EOSTable % TS % States(2) % Values = TempCompOSE
    EOSTable % TS % States(3) % Values = YpCompOSE ! This will actually become Ye!!!!
    EOSTable % TS % States(4) % Values = YmGrid

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
    EOSTable % DV % Names(1:20) = (/'Pressure                        ', &
    'Entropy Per Baryon              ', &
    'Internal Energy Density         ', &
    'Electron Chemical Potential     ', &
    'Muon Chemical Potential         ', &
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
    EOSTable % DV % Indices % iMuonChemicalPotential = 5 ! careful because there is a mismatch with the EOSCompOSE array because of this
    EOSTable % DV % Indices % iProtonChemicalPotential = 6
    EOSTable % DV % Indices % iNeutronChemicalPotential = 7
    EOSTable % DV % Indices % iProtonMassFraction = 8
    EOSTable % DV % Indices % iNeutronMassFraction = 9
    EOSTable % DV % Indices % iAlphaMassFraction = 10
    EOSTable % DV % Indices % iHeavyMassFraction = 11
    EOSTable % DV % Indices % iHeavyChargeNumber = 12
    EOSTable % DV % Indices % iHeavyMassNumber = 13
    EOSTable % DV % Indices % iHeavyBindingEnergy = 14
    EOSTable % DV % Indices % iThermalEnergy = 15
    EOSTable % DV % Indices % iGamma1 = 16
    EOSTable % DV % Indices % iNeutronEffMass = 17
    EOSTable % DV % Indices % iProtonEffMass = 18
    EOSTable % DV % Indices % iNeutronSelfEnergy = 19
    EOSTable % DV % Indices % iProtonSelfEnergy = 20
    
    PRINT*, "Allocate Dependent Variable Units " 
    EOSTable % DV % Units(1:20) = (/'Dynes per cm^2                  ', &
    'k_b per baryon                  ', &
    'erg per gram                    ', &
    'MeV                             ', &
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

    iCount = 0
    !$OMP PARALLEL DO PRIVATE(ElectronGasState, PhotonGasState, MuonGasState, Ye, Ym, Yp, iYp, dYp, LocalOffset, Xbary) &
    !$OMP SHARED(iCount) &
    !$OMP COLLAPSE(3)
    DO iRho=1,nPoints(1)
      DO iTemp=1,nPoints(2)
        DO iYe=1,nPoints(3)-1
          DO iYm=1,nPoints(4)

            ! Safely increment iCount without locking the whole loop
            !$OMP ATOMIC
            iCount = iCount + 1

            ! Print progress only every 1% of completion
            IF (MOD(iCount, INT(0.01d0 * nPoints(1) * nPoints(2) * (nPoints(3)-1) * nPoints(4))) .EQ. 0) THEN
              !$OMP CRITICAL
              WRITE(*,*) "Progress:", DBLE(iCount) / DBLE(nPoints(1) * nPoints(2) * (nPoints(3)-1) * nPoints(4)) * 100.0d0, "%"
              !$OMP END CRITICAL
            ENDIF

            Ye = EOSTable % TS % States(3) % Values(iYe)
            Ym = EOSTable % TS % States(4) % Values(iYm)
            IF (Ye + Ym > YpCompOSE(nPoints(3))) THEN
              Ym = 0.0d0
            ENDIF
            Yp = Ye + Ym

            ! Now get index for baryonic table
            CALL GetIndexAndDelta_Lin( Yp, YpCompOSE, iYp, dYp )

            PhotonGasState % T   = EOSTable % TS % States(2) % Values(iTemp)
            PhotonGasState % rho = EOSTable % TS % States(1) % Values(iRho)
            CALL PhotonGasEOS(PhotonGasState)
                
            ElectronGasState % T    = EOSTable % TS % States(2) % Values(iTemp)
            ElectronGasState % rho  = EOSTable % TS % States(1) % Values(iRho)
            ElectronGasState % yL   = Ye

            CALL LeptonGasEOS(HelmTableElectrons, ElectronGasState)

            MuonGasState % T    = EOSTable % TS % States(2) % Values(iTemp)
            MuonGasState % rho  = EOSTable % TS % States(1) % Values(iRho)
            MuonGasState % yL   = Ym

            CALL LeptonGasEOS(HelmTableMuons, MuonGasState)

            ! PRESSURE
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,1) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF
            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,1) + LocalOffset), Xbary )
          
            EOSTable % DV % Variables(1) % Values(iRho,iTemp,iYe,iYm) = &
              Xbary + MuonGasState % p + PhotonGasState % p + ElectronGasState % p

            IF (EOSTable % DV % Variables(1) % Values(iRho,iTemp,iYe,iYm) < 0.0d0) THEN
              WRITE(*,*) iYe, iYp, Ye, Ym, Yp, dYp, Xbary, MuonGasState % p, PhotonGasState % p + ElectronGasState % p
              STOP
            ENDIF

            ! ENTROPY
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,2) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF
            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,2) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(2) % Values(iRho,iTemp,iYe,iYm) = &
              Xbary + MuonGasState % s + PhotonGasState % s + ElectronGasState % s

            ! ENERGY
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,3) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF
            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,3) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(3) % Values(iRho,iTemp,iYe,iYm) = &
              Xbary + MuonGasState % e + PhotonGasState % e + ElectronGasState % e + Add_to_energy

            ! ELECTRON CHEMICAL POTENTIAL
            EOSTable % DV % Variables(4) % Values(iRho,iTemp,iYe,iYm) = &
              ElectronGasState % mu

            ! MUON CHEMICAL POTENTIAL
            EOSTable % DV % Variables(5) % Values(iRho,iTemp,iYe,iYm) = &
              MuonGasState % mu

            ! PROTON CHEMICAL POTENTIAL
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,5) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,5) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(6) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! NEUTRON CHEMICAL POTENTIAL
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,6) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,6) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(7) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! PROTON MASS FRACTION
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,7) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            ! This is a workaround for now
            CALL LinearInterp_Array_Point( 1, dYp, 0.0d0, &
              EOSCompOSE(iRho,iTemp,iYp:iYp+1,7), Xbary )
            Xbary = LOG10( MAX(Xbary, 1.0d-99) )

            ! CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
            !   LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,7) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(8) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! NEUTRON MASS FRACTION
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,8) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            ! This is a workaround for now
            CALL LinearInterp_Array_Point( 1, dYp, 0.0d0, &
              EOSCompOSE(iRho,iTemp,iYp:iYp+1,8), Xbary )
            Xbary = LOG10( MAX(Xbary, 1.0d-99) )

            ! CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
            !   LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,8) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(9) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! ALPHA MASS FRACTION
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,9) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF
            
            ! This is a workaround for now
            CALL LinearInterp_Array_Point( 1, dYp, 0.0d0, &
              EOSCompOSE(iRho,iTemp,iYp:iYp+1,9), Xbary )
            Xbary = LOG10( MAX(Xbary, 1.0d-99) )

            ! CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
            !   LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,9) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(10) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! HEAVY MASS FRACTION
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,10) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            ! This is a workaround for now
            CALL LinearInterp_Array_Point( 1, dYp, 0.0d0, &
              EOSCompOSE(iRho,iTemp,iYp:iYp+1,10), Xbary )
            Xbary = LOG10( MAX(Xbary, 1.0d-99) )

            ! CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
            !   LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,10) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(11) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! HEAVY CHARGE NUMBER
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,11) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,11) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(12) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! HEAVY MASS NUMBER
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,12) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,12) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(13) % Values(iRho,iTemp,iYe,iYm) = Xbary
            ! HEAVY BINDING ENERGY
            EOSTable % DV % Variables(14) % Values(iRho,iTemp,iYe,iYm) = 1.0d-99
            ! THERMAL ENERGY
            EOSTable % DV % Variables(15) % Values(iRho,iTemp,iYe,iYm) = 1.0d-99

            ! NEUTRON EFF MASS
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,16) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,16) + LocalOffset), Xbary )

            EOSTable % DV % Variables(17) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! PROTON EFF MASS
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,17) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,17) + LocalOffset), Xbary )

            EOSTable % DV % Variables(18) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! PROTON SELF ENERGY
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,18) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,18) + LocalOffset), Xbary )

            EOSTable % DV % Variables(19) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! NEUTRON SELF ENERGY
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYp:iYp+1,19) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYp:iYp+1,19) + LocalOffset), Xbary )

            EOSTable % DV % Variables(20) % Values(iRho,iTemp,iYe,iYm) = Xbary

          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! HANDLE LAST POINT AD HOC
    DO iRho=1,nPoints(1)
      DO iTemp=1,nPoints(2)
        DO iYm=1,nPoints(4)

          PhotonGasState % T   = EOSTable % TS % States(2) % Values(iTemp)
          PhotonGasState % rho = EOSTable % TS % States(1) % Values(iRho)
          CALL PhotonGasEOS(PhotonGasState)
              
          ElectronGasState % T    = EOSTable % TS % States(2) % Values(iTemp)
          ElectronGasState % rho  = EOSTable % TS % States(1) % Values(iRho)
          ElectronGasState % yL   = EOSTable % TS % States(3) % Values(nPoints(3))

          CALL LeptonGasEOS(HelmTableElectrons, ElectronGasState)

          MuonGasState % T    = EOSTable % TS % States(2) % Values(iTemp)
          MuonGasState % rho  = EOSTable % TS % States(1) % Values(iRho)
          MuonGasState % yL   = EOSTable % TS % States(4) % Values(iYm)

          CALL LeptonGasEOS(HelmTableMuons, MuonGasState)
          
          EOSTable % DV % Variables(1) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),1) + &
            ElectronGasState % p + MuonGasState % p + PhotonGasState % p
          EOSTable % DV % Variables(2) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),2) + &
            ElectronGasState % s + MuonGasState % s + PhotonGasState % s
          EOSTable % DV % Variables(3) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),3) + &
            ElectronGasState % e + MuonGasState % e + PhotonGasState % e + Add_to_energy
          EOSTable % DV % Variables(4) % Values(iRho,iTemp,nPoints(3),iYm) = &
            ElectronGasState % mu
          EOSTable % DV % Variables(5) % Values(iRho,iTemp,nPoints(3),iYm) = &
            MuonGasState % mu
          EOSTable % DV % Variables(6) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),5)
          EOSTable % DV % Variables(7) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),6)
          EOSTable % DV % Variables(8) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),7)
          EOSTable % DV % Variables(9) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),8)
          EOSTable % DV % Variables(10) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),9)
          EOSTable % DV % Variables(11) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),10)
          EOSTable % DV % Variables(12) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),11)
          EOSTable % DV % Variables(13) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),12)
          EOSTable % DV % Variables(14) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),13)
          EOSTable % DV % Variables(15) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),14)
          EOSTable % DV % Variables(16) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),15)
          EOSTable % DV % Variables(17) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),16)
          EOSTable % DV % Variables(18) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),17)
          EOSTable % DV % Variables(19) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),18)
          EOSTable % DV % Variables(20) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),19)
        ENDDO
      ENDDO
    ENDDO

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
    !$OMP PARALLEL DO PRIVATE(D, T, Ye, Ym, Gamma, Cs) &
    !$OMP SHARED(EOSTable, RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE, &
    !$OMP LogEnergy, LogEntropy, OS_P, OS_E, OS_S, HelmTableElectrons, HelmTableMuons, iCount) &
    !$OMP COLLAPSE(3)
    DO iRho=1,nPoints(1)
      DO iTemp=1,nPoints(2)
        DO iYe=1,nPoints(3)
          DO iYm=1,nPoints(4)

            ! Safely increment iCount without locking the whole loop
            !$OMP ATOMIC
            iCount = iCount + 1

            ! Print progress only every 1% of completion
            IF (MOD(iCount, INT(0.01d0 * nPoints(1) * nPoints(2) * nPoints(3) * nPoints(4)   )) .EQ. 0) THEN
              !$OMP CRITICAL
              WRITE(*,*) "Progress:", DBLE(iCount) / DBLE( nPoints(1) * nPoints(2) * nPoints(3) * nPoints(4)) * 100.0d0, "%"
              !$OMP END CRITICAL
            ENDIF

            D  = EOSTable % TS % States(1) % Values(iRho)
            T  = EOSTable % TS % States(2) % Values(iTemp)
            Ye = EOSTable % TS % States(3) % Values(iYe)
            Ym = EOSTable % TS % States(4) % Values(iYm)
            IF (Ye + Ym > YpCompOSE(nPoints(3))) THEN
              Ym = 0
            ENDIF
          
            CALL CalculateSoundSpeed( D, T, Ye, Ym, RhoCompOSE(:), &
                  TempCompOSE(:), YpCompOSE(:), &
                  EOSCompOSE(:,:,:,1), OS_P, &
                  LogEnergy (:,:,:), OS_E, &
                  LogEntropy(:,:,:), OS_S, &
                  HelmTableElectrons, HelmTableMuons, Gamma, Cs, .FALSE.)

            EOSTable % DV % Variables(15) % Values(iRho,iTemp,iYe,iYm) = Gamma
          ENDDO
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
            ! Sometimes mass fractions are very small and negative, just set them to very small number
            IF (Minimum_Value .gt. -1.0d-15) THEN
              EOSTable % DV % Variables(iVars) % Values = MAX(1.0d-99, EOSTable % DV % Variables(iVars) % Values)
              EOSTable % DV % Offsets(iVars) = 0.0d0
            ELSE
              EOSTable % DV % Offsets(iVars) = -1.1_dp * Minimum_Value
              EOSTable % DV % Variables(iVars) % Values =  &
                  EOSTable % DV % Variables(iVars) % Values + & 
                  EOSTable % DV % Offsets(iVars)
            ENDIF
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
    WRITE (*,*) "Starting HDF write: Baryonic EOS"
    CALL WriteEquationOfStateTableHDF( EOSTable, BaryonEOSTableName )
    CALL FinalizeHDF( )

    WRITE (*,*) "HDF write successful"
    
    CALL DeAllocateEquationOfStateTable( EOSTable )
    CALL DeallocateHelmholtzTable( HelmTableElectrons )
    CALL DeallocateHelmholtzTable( HelmTableMuons )

    ! Now Create Electron EOS

END PROGRAM wlCreateEquationOfStateTable

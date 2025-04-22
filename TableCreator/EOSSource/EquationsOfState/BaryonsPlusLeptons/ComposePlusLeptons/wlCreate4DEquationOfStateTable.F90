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
    USE wlCompOSEInterface, ONLY : ReadnPointsFromCompOSE, ReadCompOSETable, &
        ReadCompOSEHDFTable, RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE
    USE wlHelmMuonIOModuleHDF, ONLY : WriteHelmholtzTableHDF, WriteMuonTableHDF
    USE wlLeptonEOSModule, ONLY: HelmholtzTableType, MuonEOSType, &
        iTempMax, iDenMax, &
        nTempMuon, nDenMuon, &
        ReadHelmEOSdat, ReadMuonEOSdat, &
        AllocateHelmholtzTable, DeallocateHelmholtzTable, &
        AllocateMuonEOS, DeAllocateMuonEOS 
    USE wlElectronPhotonEOS, ONLY: &
      ElectronPhotonEOS, ElectronPhotonStateType
    USE wlMuonEOS, ONLY: &
      FullMuonEOS, MuonStateType
    USE wlEosConstantsModule, ONLY: cvel, ergmev, cm3fm3, kmev_inv, rmu, mn, me, mp
    USE wlSoundSpeedModule

    IMPLICIT NONE
    
    INTEGER                          :: iVars
    INTEGER, DIMENSION(2)            :: nPoints2D
    INTEGER, DIMENSION(3)            :: nPointsCompose
    INTEGER, DIMENSION(4)            :: nPoints
    INTEGER                          :: nVariables
    TYPE(EquationOfState4DTableType) :: EOSTable
    TYPE(HelmholtzTableType)         :: HelmholtzTable
    TYPE(MuonEOSType)                :: MuonTable
    TYPE(ElectronPhotonStateType)    :: ElectronPhotonState
    TYPE(MuonStateType)              :: MuonState
    
    CHARACTER(len=128) :: CompOSEFilePath, CompOSEFHDF5Path, &
        HelmDatFilePath, MuonDatFilePath, &
        BaryonEOSTableName
    
    REAL(dp) :: Minimum_Value, Add_to_energy
    LOGICAL  :: RedHDF5Table
    INTEGER  :: iLepton, i, iRho, iTemp, iYe, iYm, iCount

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: YmGrid
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: LogEnergy, LogEntropy

    REAL(DP) :: D, T, Yp, Ye, Ym
    REAL(DP) :: Xbary, dYp, OS_P, OS_E, OS_S, LocalOffset, Gamma, Cs

    RedHDF5Table = .false.
    
    Add_to_energy = 8.9d0*ergmev/rmu
    Add_to_energy = 2.0d0*ergmev/rmu
    Add_to_energy = + cvel**2 - ergmev * mn / rmu + 8.9d0 * ergmev / rmu

    ! path of the original compose table
    CompOSEFilePath = 'SFHo_no_ele/'
    CompOSEFHDF5Path = 'SFHo_no_ele/eoscompose.h5'
    IF (RedHDF5Table) THEN
        BaryonEOSTableName = '4DEOSTableHighres.h5'
    ELSE
        BaryonEOSTableName = '4DEOSTable.h5'
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
    nPoints2D = (/ iTempMax, iDenMax /)
    PRINT*, "Allocate Helmholtz EOS"
    CALL AllocateHelmholtzTable( HelmholtzTable, nPoints2D )
    
    HelmDatFilePath = '../helm_table.dat'
    CALL ReadHelmEOSdat( HelmDatFilePath, HelmholtzTable )
    
    ! ------------- NOW DO MUON EOS ------------------ !
    nPoints2D = (/ nTempMuon, nDenMuon /)
    PRINT*, "Allocate Muon EOS"
    CALL AllocateMuonEOS( MuonTable, nPoints2D )
    
    MuonDatFilePath = '../muons_fixedrho.dat'
    CALL ReadMuonEOSdat( MuonDatFilePath, MuonTable )

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
    CALL AllocateEquationOfState4DTable( EOSTable, nPoints , nVariables )
             
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

    iCount = 0
    !$OMP PARALLEL DO PRIVATE(ElectronPhotonState, MuonState, Ye, Ym, Yp, dYp, LocalOffset, Xbary) &
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
              Ym = 0
            ENDIF
            Yp = Ye + Ym
            dYp = Yp - Ye

            ElectronPhotonState % rho = EOSTable % TS % States(1) % Values(iRho)
            ElectronPhotonState % t   = EOSTable % TS % States(2) % Values(iTemp)
            ElectronPhotonState % ye  = Ye
            CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

            MuonState % t     = EOSTable % TS % States(2) % Values(iTemp)
            MuonState % rhoym = EOSTable % TS % States(1) % Values(iRho) * Ym
            CALL FullMuonEOS(MuonTable, MuonState)

            ! PRESSURE
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,1) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF
            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,1) + LocalOffset), Xbary )
          
            EOSTable % DV % Variables(1) % Values(iRho,iTemp,iYe,iYm) = &
            Xbary + MuonState % p + ElectronPhotonState % p

            IF (iRho == 200 .and. iYe == 48 .and. iYm == 1) THEN
              WRITE(*,*) Xbary , MuonState % p , ElectronPhotonState % p, Ym
              WRITE(*,*) iRho,iTemp,iYe,iYm
            ENDIF

            ! ENTROPY
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,2) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF
            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,2) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(2) % Values(iRho,iTemp,iYe,iYm) = &
              Xbary + MuonState % s + ElectronPhotonState % s

            IF (iYm == 1 .and. MuonState % s / (Xbary + ElectronPhotonState % s) > 0.1d0) THEN
              WRITE(*,*) Xbary , MuonState % s , ElectronPhotonState % s, Ym
              WRITE(*,*) iRho,iTemp,iYe,iYm
            ENDIF

            ! ENERGY
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,3) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF
            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,3) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(3) % Values(iRho,iTemp,iYe,iYm) = &
              Xbary + MuonState % e + ElectronPhotonState % e + Add_to_energy

            ! ELECTRON CHEMICAL POTENTIAL
            EOSTable % DV % Variables(4) % Values(iRho,iTemp,iYe,iYm) = &
              ElectronPhotonState % mue

            ! PROTON CHEMICAL POTENTIAL
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,5) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,5) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(5) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! NEUTRON CHEMICAL POTENTIAL
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,6) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,6) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(6) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! PROTON MASS FRACTION
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,7) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,7) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(7) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! NEUTRON MASS FRACTION
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,8) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,8) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(8) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! ALPHA MASS FRACTION
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,9) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,9) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(9) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! HEAVY MASS FRACTION
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,10) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,10) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(10) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! HEAVY CHARGE NUMBER
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,11) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,11) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(11) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! HEAVY MASS NUMBER
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,12) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,12) + LocalOffset), Xbary )
            
            EOSTable % DV % Variables(12) % Values(iRho,iTemp,iYe,iYm) = Xbary
            ! HEAVY BINDING ENERGY
            EOSTable % DV % Variables(13) % Values(iRho,iTemp,iYe,iYm) = Xbary
            ! THERMAL ENERGY
            EOSTable % DV % Variables(14) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! PROTON SELF ENERGY
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,16) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,16) + LocalOffset), Xbary )

            EOSTable % DV % Variables(16) % Values(iRho,iTemp,iYe,iYm) = Xbary

            ! NEUTRON SELF ENERGY
            LocalOffset = MINVAL( EOSCompOSE(iRho,iTemp,iYe:iYe+1,17) )
            IF (LocalOffset .lt. 0.0_dp) THEN
                LocalOffset = -1.1d0*LocalOffset
            ELSE IF (LocalOffset .eq. 0.0_dp) THEN
                LocalOffset = 1.0d-10
            ELSE
                LocalOffset = 0.0_dp
            ENDIF

            CALL LinearInterp_Array_Point( 1, dYp, LocalOffset, &
              LOG10(EOSCompOSE(iRho,iTemp,iYe:iYe+1,17) + LocalOffset), Xbary )

            EOSTable % DV % Variables(17) % Values(iRho,iTemp,iYe,iYm) = Xbary

          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! HANDLE LAST POINT AD HOC
    DO iRho=1,nPoints(1)
      DO iTemp=1,nPoints(2)
        DO iYm=1,nPoints(4)

          ElectronPhotonState % rho = EOSTable % TS % States(1) % Values(iRho)
          ElectronPhotonState % t   = EOSTable % TS % States(2) % Values(iTemp)
          ElectronPhotonState % ye  = EOSTable % TS % States(3) % Values(nPoints(3))
          CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

          EOSTable % DV % Variables(1) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),1) + ElectronPhotonState % p
          EOSTable % DV % Variables(2) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),2) + ElectronPhotonState % s
          EOSTable % DV % Variables(3) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),3) + ElectronPhotonState % e
          EOSTable % DV % Variables(4) % Values(iRho,iTemp,nPoints(3),iYm) = &
            ElectronPhotonState % mue
          EOSTable % DV % Variables(5) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),5)
          EOSTable % DV % Variables(6) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),6)
          EOSTable % DV % Variables(7) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),7)
          EOSTable % DV % Variables(8) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),8)
          EOSTable % DV % Variables(9) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),9)
          EOSTable % DV % Variables(10) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),10)
          EOSTable % DV % Variables(11) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),11)
          EOSTable % DV % Variables(12) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),12)
          EOSTable % DV % Variables(13) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),13)
          EOSTable % DV % Variables(14) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),14)
          EOSTable % DV % Variables(15) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),15)
          EOSTable % DV % Variables(16) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),16)
          EOSTable % DV % Variables(17) % Values(iRho,iTemp,nPoints(3),iYm) = &
            EOSCompOSE(iRho,iTemp,nPoints(3),17)
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
    !$OMP LogEnergy, LogEntropy, OS_P, OS_E, OS_S, HelmholtzTable, MuonTable, iCount) &
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
                  HelmholtzTable, MuonTable, Gamma, Cs, .FALSE.)

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
            EOSTable % DV % Offsets(iVars) = -1.1_dp * Minimum_Value
            EOSTable % DV % Variables(iVars) % Values =  &
                EOSTable % DV % Variables(iVars) % Values + & 
                EOSTable % DV % Offsets(iVars)
        ELSE IF ( Minimum_Value .gt. 0.0_dp) THEN
            EOSTable % DV % Offsets(iVars) = 0.0_dp
        ELSE
          WRITE(*,*) 'Minimum is zero', iVars
          EOSTable % DV % Variables(iVars) % Values = & 
              MAX(EOSTable % DV % Variables(iVars) % Values, 1.0e-30_dp)
        ENDIF

        WRITE(*,*) iVars, EOSTable % DV % Offsets(iVars), Minimum_Value

        EOSTable % DV % Variables(iVars) % Values = &
          LOG10(EOSTable % DV % Variables(iVars) % Values)

    END DO
   
    ! NOW CREATE BARYONIC FILE
    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Baryonic EOS"
    CALL WriteEquationOfState4DTableHDF( EOSTable, BaryonEOSTableName )
    CALL FinalizeHDF( )

    WRITE (*,*) "HDF write successful"
    
    CALL DeAllocateEquationOfState4DTable( EOSTable )
    CALL DeallocateHelmholtzTable( HelmholtzTable )
    CALL DeAllocateMuonEOS( MuonTable )

    ! Now Create Electron EOS

END PROGRAM wlCreateEquationOfStateTable

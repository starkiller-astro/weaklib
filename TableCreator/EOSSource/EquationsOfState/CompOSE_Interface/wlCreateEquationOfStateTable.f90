PROGRAM wlCreateEquationOfStateTable
    
    USE wlKindModule, ONLY: dp
    USE HDF5
    USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
    USE wlEOSIOModuleHDF
    USE wlExtNumericalModule, ONLY: zero
    USE wlCompOSEInterface, ONLY : ReadnPointsFromCompOSE, ReadCompOSETable, &
	RhoCompOSE, TempCompOSE, YpCompOSE, EOSCompOSE
    
    implicit none
    
    INTEGER                        :: iVars
    INTEGER, DIMENSION(3)          :: nPointsBaryon, nPointsElectron
    INTEGER                        :: nVariables
    TYPE(EquationOfStateTableType) :: EOSBaryonTable
    TYPE(EquationOfStateTableType) :: EOSElectronTable
    
    CHARACTER(len=128) :: CompOSEFilePath
    
    REAL(dp) :: Minimum_Value
    
    ! These are the parameters for the helmoltz EOS 
    INTEGER, PARAMETER :: EOSiTMax=541, EOSiDenMax=201
    INTEGER :: istat=0
    INTEGER :: iDen, iT
    REAL(dp) :: Helm_f(EOSiDenMax,EOSiTMax), Helm_df_drho(EOSiDenMax,EOSiTMax), &
          Helm_df_dt(EOSiDenMax,EOSiTMax), Helm_d2f_d2rho(EOSiDenMax,EOSiTMax), &
          Helm_d2f_d2t(EOSiDenMax,EOSiTMax), Helm_d2f_drhodt(EOSiDenMax,EOSiTMax), & 
          Helm_d3f_d2rhodt(EOSiDenMax,EOSiTMax), Helm_d3f_drhod2t(EOSiDenMax,EOSiTMax), &
          Helm_d4f_d2rhod2t(EOSiDenMax,EOSiTMax)
    
    !!
    !!     f --  Helmholtz free energy
    !!     fd --  derivative of f wrt density
    !!     ft --  derivative of f wrt temperature
    !!     fdd --  second derivative of f wrt density
    !!     ftt --  second derivative of f wrt temperature
    !!     fdt --  second derivative of f wrt density and temperature
    !!     fddt --  third derivative of f wrt density^2 and temperature
    !!     fdtt --  third derivative of f wrt density and temperature^2 e.g. dF/(dd)(dt^2)
    !!     fddtt --  fourth derivative of f wrt density^2 and temperature^2
    !!     dpdf --  pressure derivative 
    !!     dpdfd -- 
    !!     dpdft --  
    !!     dpdfdt --  
    !!     ef --  electron chemical potential
    !!     efd --  
    !!     eft --  
    !!     efdt --  
    !!     xf --  number density
    !!     xfd --  
    !!     xft --  
    !!     xfdt --
    
    ! START THE CODE
    
    CompOSEFilePath = '../SFHo'
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
	
    nVariables = 15	
    nPointsElectron = (\ nPointsBaryon(1), nPointsBaryon(2), 1/)
    PRINT*, "Allocate Electron EOS"
    CALL AllocateEquationOfStateTable( EOSElectronTable, nPointsElectron , nVariables )
    
    EOSBaryonTable % TS % Names(1:3) = (/'Density                         ',&
    'Temperature                     ',&
    'Proton Fraction                 '/)
    
    EOSBaryonTable % TS % Indices % iRho = 1
    EOSBaryonTable % TS % Indices % iT   = 2
    EOSBaryonTable % TS % Indices % iYe  = 3
    
    OPEN(UNIT=1234,FILE='helm_table.dat',STATUS='old',IOSTAT=istat)
    
    IF (istat .ne. 0) THEN
        WRITE(*,*) 'Cannot open helm_table.dat!'
        STOP
    ENDIF
    
    !..read the helmholtz free energy table
    do iT=1,EOSiTMax
        do iDen=1,EOSiDenMax
            read(uniDentEos,*) Helm_f(iDen,iT), Helm_df_drho(iDen,iT), Helm_df_dt(iDen,iT),&
            Helm_d2f_d2rho(iDen,iT), Helm_d2f_d2t(iDen,iT), Helm_d2f_drhodt(iDen,iT), & 
            Helm_d3f_d2rhodt(iDen,iT), Helm_d3f_drhod2t(iDen,iT), Helm_d4f_d2rhod2t(iDen,iT)
        enddo
    enddo
    
    !..read the pressure derivative with density table
    DO iT=1,EOSiTMax
        DO iDen=1,EOSiDenMax
            read(uniDentEos,*) Helm_dp_df(iDen,iT), Helm_d2p_dfdrho(iDen,iT),&
            Helm_d2p_dfdt(iDen,iT), Helm_d3p_dfdrhodt(iDen,iT) ! not sure about the nomenclature, whiDench ones are 2nd or 3rd deriDenvs
        ENDDO
    ENDDO
    
    !..read the electron chemical potential table
    DO iT=1,EOSiTMax
        DO iDen=1,EOSiDenMax
            READ(uniDentEos,*) Helm_mue(iDen,iT), Helm_dmue_drho(iDen,iT),&
            Helm_dmue_dt(iDen,iT), Helm_d2mue_drhodt(iDen,iT)
        ENDDO
    ENDDO
    
    !..read the number density table
    DO iT=1,EOSiTMax
        DO iDen=1,EOSiDenMax
            READ(uniDentEos,*) Helm_ne(iDen,iT), Helm_dne_drho(iDen,iT),&
            Helm_dne_dt(iDen,iT), Helm_d2ne_drhodt(iDen,iT)
        ENDDO
    ENDDO
    
    ! NOW WRITE H5 FILES
    
    CALL InitializeHDF( )
    
    WRITE (*,*) "Starting HDF write "
    CALL WriteEquationOfStateTableHDF( EOSBaryonTable )
    
    CALL FinalizeHDF( )
    
    WRITE (*,*) "HDF write successful"
    
    CALL DeAllocateEquationOfStateTable( EOSBaryonTable )
    
END PROGRAM wlCreateEquationOfStateTable

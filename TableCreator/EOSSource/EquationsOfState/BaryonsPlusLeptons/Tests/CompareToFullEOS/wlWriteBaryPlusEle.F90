PROGRAM wlWriteBaryPlusEle
	
	USE wlKindModule, ONLY: dp
	USE wlExtNumericalModule, ONLY: zero, one, pi
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
	USE wlEOSIOModuleHDF
	USE wlLeptonEOSModule
	USE wlHelmMuonIOModuleHDF
	USE wlElectronEOS
	USE wlExtPhysicalConstantsModule, ONLY: &
		kmev, rmu, kmev_inv, ergmev, me, cvel
	
	IMPLICIT NONE
	
	TYPE(ElectronStateType) :: ElectronState
	TYPE(HelmholtzEOSType) :: HelmholtzTable
    TYPE(EquationOfStateTableType) :: EOSBaryonTable, EOSBaryonPlusEleTable, EOSFullTable
	
	INTEGER :: nRho, nTemp, nYp, iAbar, iZbar, iPressure, iRho, iTemp, iYp, iXa, iXn, iXh, &
		iEnergy, iEntropy, iMue, iMun, iMup, iRho_transition, ics2, iGamma
	REAL(dp) :: abar_total, zbar_total, abar, zbar, heavy_frac, p_frac, n_frac, alpha_frac
	
	REAL(dp) :: press_baryons, press_electrons, press_total
	REAL(dp) :: s_baryons, s_electrons, s_total
	REAL(dp) :: eps_baryons, eps_electrons, eps_total
	REAL(dp) :: cs2_baryons, cs2_electrons, cs2_total
	REAL(dp) :: gamma_baryons, gamma_electrons, gamma_total
	REAL(dp) :: mue_baryons, mue_electrons, mue_total
	REAL(dp) :: mup_baryons, mup_total, mun_baryons, mun_total

	REAL(dp) :: ematch
	REAL(dp), ALLOCATABLE :: eos_table(:,:,:,:)
	
	LOGICAL :: include_ion_contribution = .false.
	LOGICAL :: do_coulomb = .false.
	
    CHARACTER(len=128) :: BaryonEOSTableName, BaryonPlusEleName, FullEOSTableName
	
    BaryonEOSTableName = 'BaryonsPlusHelmEOS.h5'
	FullEOSTableName = 'FullEOS.h5' 
	IF (do_coulomb) THEN
		BaryonPlusEleName = 'BaryonPlusEleEOS_Coul.h5' 
	ELSE
		BaryonPlusEleName = 'BaryonPlusEleEOS_NoCoul.h5' 
	ENDIF
	
	! read in helmholtz table
	CALL ReadHelmholtzTableHDF( HelmholtzTable, BaryonEOSTableName )
	
	! read in baryon table -------------------------------
    CALL ReadEquationOfStateTableHDF( EOSBaryonTable, BaryonEOSTableName )

	! read in full table -------------------------------
    CALL ReadEquationOfStateTableHDF( EOSFullTable, FullEOSTableName )
	
	nRho = EOSBaryonTable % TS % nPoints(1)
	nTemp = EOSBaryonTable % TS % nPoints(2)
	nYp = EOSBaryonTable % TS % nPoints(3)
	
    iXn = EOSBaryonTable % DV % Indices % iNeutronMassFraction
    iXa = EOSBaryonTable % DV % Indices % iAlphaMassFraction
    iXh = EOSBaryonTable % DV % Indices % iHeavyMassFraction
	iAbar = EOSBaryonTable % DV % Indices % iHeavyMassNumber
	iZbar = EOSBaryonTable % DV % Indices % iHeavyChargeNumber
	iPressure = EOSBaryonTable % DV % Indices % iPressure
	iEntropy = EOSBaryonTable % DV % Indices % iEntropyPerBaryon
	iEnergy = EOSBaryonTable % DV % Indices % iInternalEnergyDensity
	iMue = EOSBaryonTable % DV % Indices % iElectronChemicalPotential
	iMup = EOSBaryonTable % DV % Indices % iProtonChemicalPotential
	iMun = EOSBaryonTable % DV % Indices % iNeutronChemicalPotential
	ics2 = EOSBaryonTable % DV % Indices % iThermalEnergy
	iGamma = EOSBaryonTable % DV % Indices % iGamma1

    PRINT*, "Allocate Baryon + Ele EOS"
    CALL AllocateEquationOfStateTable( EOSBaryonPlusEleTable, (/nRho, nTemp, nYp/) , 15 )
    
    EOSBaryonPlusEleTable % TS % Names(1:3) = (/'Density                         ',&
    'Temperature                     ',&
    'Proton Fraction                 '/)
    
    EOSBaryonPlusEleTable % TS % Indices % iRho = 1
    EOSBaryonPlusEleTable % TS % Indices % iT   = 2
    EOSBaryonPlusEleTable % TS % Indices % iYe  = 3
    
    PRINT*, "Allocate Independent Variable Units " 
    
    EOSBaryonPlusEleTable % TS % Units(1:3) = (/'Grams per cm^3                  ', &
    'K                               ', &
    '                                '/) 
    
    PRINT*, "Label Grid Type"
    EOSBaryonPlusEleTable % TS % LogInterp(1:3) =  (/1, 1, 0/)
    
    PRINT*, "Allocate Names " 
    EOSBaryonPlusEleTable % DV % Names(1:15) = (/'Pressure                        ', &
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
    
    EOSBaryonPlusEleTable % DV % Indices % iPressure = 1
    EOSBaryonPlusEleTable % DV % Indices % iEntropyPerBaryon = 2
    EOSBaryonPlusEleTable % DV % Indices % iInternalEnergyDensity = 3
    EOSBaryonPlusEleTable % DV % Indices % iElectronChemicalPotential = 4
    EOSBaryonPlusEleTable % DV % Indices % iProtonChemicalPotential = 5
    EOSBaryonPlusEleTable % DV % Indices % iNeutronChemicalPotential = 6
    EOSBaryonPlusEleTable % DV % Indices % iProtonMassFraction = 7
    EOSBaryonPlusEleTable % DV % Indices % iNeutronMassFraction = 8
    EOSBaryonPlusEleTable % DV % Indices % iAlphaMassFraction = 9
    EOSBaryonPlusEleTable % DV % Indices % iHeavyMassFraction = 10
    EOSBaryonPlusEleTable % DV % Indices % iHeavyChargeNumber = 11
    EOSBaryonPlusEleTable % DV % Indices % iHeavyMassNumber = 12
    EOSBaryonPlusEleTable % DV % Indices % iHeavyBindingEnergy = 13
    EOSBaryonPlusEleTable % DV % Indices % iThermalEnergy = 14
    EOSBaryonPlusEleTable % DV % Indices % iGamma1 = 15
    
    PRINT*, "Allocate Dependent Variable Units " 
    EOSBaryonPlusEleTable % DV % Units(1:15) = (/'Dynes per cm^2                  ', &
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

    PRINT*, "Begin Populating EOSBaryonPlusEleTable with do_coulomb = ", do_coulomb 
    EOSBaryonPlusEleTable % TS % States(1) % Values = EOSBaryonTable % TS % States(1) % Values
    EOSBaryonPlusEleTable % TS % States(2) % Values = EOSBaryonTable % TS % States(2) % Values
    EOSBaryonPlusEleTable % TS % States(3) % Values = EOSBaryonTable % TS % States(3) % Values
    
	DO iRho=1,nRho
		DO iTemp=1,nTemp
			DO iYp=1,nYp
				
				abar = EOSBaryonTable % DV % Variables(iAbar) % Values(iRho,iTemp,iYp)
				zbar = EOSBaryonTable % DV % Variables(iZbar) % Values(iRho,iTemp,iYp)
				heavy_frac = EOSBaryonTable % DV % Variables(iXh) % Values(iRho,iTemp,iYp)
				alpha_frac = EOSBaryonTable % DV % Variables(iXa) % Values(iRho,iTemp,iYp)
				p_frac = EOSBaryonTable % TS % States(3) % Values(iYp)
				n_frac = EOSBaryonTable % DV % Variables(iXn) % Values(iRho,iTemp,iYp)
				
				abar_total = abar * heavy_frac + p_frac + n_frac + 4.0d0*alpha_frac
				zbar_total = zbar * heavy_frac + p_frac + 2.0d0*alpha_frac
				
				! Initialize temperature, density, yp, Zbar and Abar
				ElectronState % t = EOSBaryonTable % TS % States(2) % Values(iTemp)
				ElectronState % rho = EOSBaryonTable % TS % States(1) % Values(iRho)
				ElectronState % abar = abar_total
				ElectronState % zbar = zbar_total
				ElectronState % y_e = EOSBaryonTable % TS % States(3) % Values(iYp)
				
				! calculate electron quantities
				CALL FullHelmEOS(1, HelmholtzTable, ElectronState, include_ion_contribution, do_coulomb)
			
				press_baryons = EOSBaryonTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iPressure)
				press_electrons = ElectronState % p
				
				s_baryons = EOSBaryonTable % DV % Variables(iEntropy) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iEntropy)
				s_electrons = ElectronState % s * rmu * kmev_inv / ergmev
				
				eps_baryons = EOSBaryonTable % DV % Variables(iEnergy) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iEnergy)
				eps_electrons = ElectronState % e - 0.511 / rmu * ergmev * ElectronState % y_e

				mue_electrons = ElectronState % eta * ElectronState % t*kmev + 0.511d0
				mup_baryons = EOSBaryonTable % DV % Variables(iMup) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iMup)
				mun_baryons = EOSBaryonTable % DV % Variables(imun) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(imun)
				
				cs2_baryons = EOSBaryonTable % DV % Variables(ics2) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(ics2)
				cs2_electrons = ( ElectronState % cs )**2.0d0
			
				Gamma_baryons = EOSBaryonTable % DV % Variables(iGamma) % Values(iRho,iTemp,iYp) - &
					EOSBaryonTable % DV % Offsets(iGamma)
				Gamma_electrons = ElectronState % gam1
				
				! thermo and compo state
				EOSBaryonPlusEleTable % DV % Variables(1) % Values(iRho,iTemp,iYp) = press_baryons + press_electrons
				EOSBaryonPlusEleTable % DV % Variables(2) % Values(iRho,iTemp,iYp) = s_baryons + s_electrons
				EOSBaryonPlusEleTable % DV % Variables(3) % Values(iRho,iTemp,iYp) = eps_baryons + eps_electrons
				EOSBaryonPlusEleTable % DV % Variables(4) % Values(iRho,iTemp,iYp) = mue_electrons
				EOSBaryonPlusEleTable % DV % Variables(5) % Values(iRho,iTemp,iYp) = mup_baryons
				EOSBaryonPlusEleTable % DV % Variables(6) % Values(iRho,iTemp,iYp) = mun_baryons
				EOSBaryonPlusEleTable % DV % Variables(7) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(7) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(8) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(8) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(9) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(9) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(10) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(10) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(11) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(11) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(12) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(12) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(13) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(13) % Values(iRho,iTemp,iYp)
					
				!EOSBaryonPlusEleTable % DV % Variables(14) % Values(iRho,iTemp,iYp) = cs2_baryons + cs2_electrons
				!EOSBaryonPlusEleTable % DV % Variables(15) % Values(iRho,iTemp,iYp) = Gamma_baryons + Gamma_electrons
				
			END DO
		END DO
	END DO
		
	ALLOCATE( eos_table(nRho, nTemp, nYp, 3) )
	eos_table(:,:,:,1) = LOG10( EOSBaryonPlusEleTable % DV % Variables(1) % Values(:,:,:) )
	eos_table(:,:,:,2) = EOSBaryonPlusEleTable % DV % Variables(2) % Values(:,:,:)
	eos_table(:,:,:,3) = LOG10( EOSBaryonPlusEleTable % DV % Variables(3) % Values(:,:,:) - &
			1.1d0*MINVAL(EOSBaryonPlusEleTable % DV % Variables(3) % Values(:,:,:)) )
	
	CALL derivatives_production(nRho, nTemp, nYp, &
		LOG10( EOSBaryonPlusEleTable % TS % States(1) % Values(:) ), &
		LOG10( EOSBaryonPlusEleTable % TS % States(2) % Values(:) ), &
		EOSBaryonPlusEleTable % TS % States(2) % Values(:), &
		eos_table, &
		EOSBaryonPlusEleTable % DV % Variables(14) % Values, &
		EOSBaryonPlusEleTable % DV % Variables(15) % Values )
	
	DO iRho=1,nRho
	  DO iTemp=1,nTemp
		DO iYp=1,nYp
		  EOSBaryonPlusEleTable % DV % Variables(14) % Values(iRho,iTemp,iYp) = &
			EOSBaryonPlusEleTable % DV % Variables(15) % Values(iRho,iTemp,iYp) * &
			EOSBaryonPlusEleTable % DV % Variables(1) % Values(iRho,iTemp,iYp) / &
			(EOSBaryonPlusEleTable % DV % Variables(1) % Values(iRho,iTemp,iYp) + &
			 (EOSBaryonPlusEleTable % DV % Variables(3) % Values(iRho,iTemp,iYp) + &
			0.511d0 / rmu * ergmev * EOSBaryonPlusEleTable % TS % States(3) % Values(iYp) + cvel**2) * &
			EOSBaryonPlusEleTable % TS % States(1) % Values(iRho) )
		END DO
	  END DO
	END DO

	EOSBaryonPlusEleTable % DV % Variables(14) % Values = &
		EOSBaryonPlusEleTable % DV % Variables(14) % Values * cvel**2

	WRITE(*,*) EOSBaryonPlusEleTable % DV % Variables(15) % Values(50,50,50), EOSFullTable % DV % Variables(15) % Values(50,50,50)
	WRITE(*,*) EOSBaryonPlusEleTable % DV % Variables(14) % Values(50,50,50), EOSFullTable % DV % Variables(14) % Values(50,50,50)
	
	EOSBaryonPlusEleTable % DV % Offsets(:) = zero
    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Baryonic + Electrons EOS"
    CALL WriteEquationOfStateTableHDF( EOSBaryonPlusEleTable, BaryonPlusEleName )
    CALL FinalizeHDF( )
	
END PROGRAM wlWriteBaryPlusEle
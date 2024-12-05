PROGRAM wlCompareToFullEOS
	
	USE wlKindModule, ONLY: dp
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
	USE wlEOSIOModuleHDF	
	USE wlLeptonEOSModule
	USE wlHelmMuonIOModuleHDF
	USE wlElectronPhotonEOS
	USE wlExtPhysicalConstantsModule, ONLY: kmev, rmu, kmev_inv, ergmev, me
	
	IMPLICIT NONE
	
	TYPE(ElectronPhotonStateType) :: ElectronState
	TYPE(HelmholtzTableType) :: HelmholtzTable
    TYPE(EquationOfStateTableType) :: EOSBaryonTable, EOSFullTable, EOSBaryonPlusEleTable, EOSEleTable
	
	INTEGER :: nRho, nTemp, nYp, iAbar, iZbar, iBE, iPressure, iRho, iTemp, iYp, iXp, iXa, iXn, iXh, &
		iEnergy, iEntropy, iMue, iMun, iMup, iRho_transition, ics2, iGamma1
	REAL(dp) :: abar_total, zbar_total, abar, zbar, heavy_frac, p_frac, n_frac, alpha_frac
	
	REAL(dp) :: press_baryons, press_electrons, press_total
	REAL(dp) :: s_baryons, s_electrons, s_total
	REAL(dp) :: eps_baryons, eps_electrons, eps_total
	REAL(dp) :: mue_baryons, mue_electrons, mue_total
	REAL(dp) :: mup_baryons, mup_total, mun_baryons, mun_total
	
	REAL(DP) :: energy_shift
	REAL(dp), ALLOCATABLE :: eos_table(:,:,:,:)
	
	LOGICAL :: RedHDF5Table
	LOGICAL :: include_ion_contribution = .false.
	LOGICAL :: do_coulomb
	INTEGER :: igamma = 1
	
    CHARACTER(len=128) :: BaryonEOSTableName, ElectronEOSTableName, BaryonPlusEleName, FullEOSTableName

    RedHDF5Table = .true.
	do_coulomb = .true.
	
	IF (RedHDF5Table) THEN
		FullEOSTableName = 'FullEOS_interpolated.h5'
		BaryonEOSTableName = 'BaryonsPlusHelmEOS_interpolated.h5'
	ELSE
		FullEOSTableName = 'FullEOS.h5'
		BaryonEOSTableName = 'BaryonsPlusHelmEOS.h5'
	ENDIF

	IF (do_coulomb) THEN
		IF (RedHDF5Table) THEN
			BaryonPlusEleName = 'BaryonPlusEleEOS_Coul_interpolated.h5'
		ELSE
			BaryonPlusEleName = 'BaryonPlusEleEOS_Coul.h5'
		ENDIF
		ElectronEOSTableName = 'EOSEle_Coul.h5' 
	ELSE
		IF (RedHDF5Table) THEN
			BaryonPlusEleName = 'BaryonPlusEleEOS_NoCoul_interpolated.h5'
		ELSE
			BaryonPlusEleName = 'BaryonPlusEleEOS_NoCoul.h5'
		ENDIF
		ElectronEOSTableName = 'EOSEle_NoCoul.h5' 
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
	
    iXp = EOSBaryonTable % DV % Indices % iProtonMassFraction
    iXn = EOSBaryonTable % DV % Indices % iNeutronMassFraction
    iXa = EOSBaryonTable % DV % Indices % iAlphaMassFraction
    iXh = EOSBaryonTable % DV % Indices % iHeavyMassFraction
	iAbar = EOSBaryonTable % DV % Indices % iHeavyMassNumber
	iZbar = EOSBaryonTable % DV % Indices % iHeavyChargeNumber
	iBE = EOSBaryonTable % DV % Indices % iHeavyBindingEnergy
	iPressure = EOSBaryonTable % DV % Indices % iPressure
	iEntropy = EOSBaryonTable % DV % Indices % iEntropyPerBaryon
	iEnergy = EOSBaryonTable % DV % Indices % iInternalEnergyDensity
	iMue = EOSBaryonTable % DV % Indices % iElectronChemicalPotential
	iMup = EOSBaryonTable % DV % Indices % iProtonChemicalPotential
	iMun = EOSBaryonTable % DV % Indices % iNeutronChemicalPotential
	ics2 = EOSBaryonTable % DV % Indices % iThermalEnergy
	iGamma1 = EOSBaryonTable % DV % Indices % iGamma1
	
	! Initialize electron and combined tables 
	EOSBaryonPlusEleTable = EOSBaryonTable
	EOSEleTable = EOSBaryonTable
	
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
				ElectronState % abar = abar
				ElectronState % zbar = zbar
				ElectronState % ye = EOSBaryonTable % TS % States(3) % Values(iYp)
				
				! calculate electron quantities
				CALL FullHelmEOS(1, HelmholtzTable, ElectronState, include_ion_contribution, do_coulomb)
			
				press_baryons = EOSBaryonTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iPressure)
				press_electrons = ElectronState % p
				press_total = EOSFullTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) - &
				EOSFullTable % DV % Offsets(iPressure)
				
				s_baryons = EOSBaryonTable % DV % Variables(iEntropy) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iEntropy)
				s_electrons = ElectronState % s * rmu * kmev_inv / ergmev
				s_total = EOSFullTable % DV % Variables(iEntropy) % Values(iRho,iTemp,iYp) - &
				EOSFullTable % DV % Offsets(iEntropy)
				
				eps_baryons = EOSBaryonTable % DV % Variables(iEnergy) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iEnergy)
				eps_electrons = ElectronState % e + 0.511 / rmu * ergmev * ElectronState % ye
				eps_total = EOSFullTable % DV % Variables(iEnergy) % Values(iRho,iTemp,iYp) - &
				EOSFullTable % DV % Offsets(iEnergy)
				
				mue_baryons = EOSBaryonTable % DV % Variables(iMue) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iMue)
				mue_electrons = ElectronState % eta * ElectronState % t*kmev + 0.511d0
				mue_total = EOSFullTable % DV % Variables(iMue) % Values(iRho,iTemp,iYp) - &
				EOSFullTable % DV % Offsets(iMue)
				
				mup_baryons = EOSBaryonTable % DV % Variables(iMup) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iMup)
				mup_total = EOSFullTable % DV % Variables(iMup) % Values(iRho,iTemp,iYp) - &
				EOSFullTable % DV % Offsets(iMup)
				
				mun_baryons = EOSBaryonTable % DV % Variables(imun) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(imun)
				mun_total = EOSFullTable % DV % Variables(imun) % Values(iRho,iTemp,iYp) - &
				EOSFullTable % DV % Offsets(imun)
				
				EOSBaryonPlusEleTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) = press_baryons + press_electrons
				EOSBaryonPlusEleTable % DV % Variables(iEntropy) % Values(iRho,iTemp,iYp) = s_baryons + s_electrons
				EOSBaryonPlusEleTable % DV % Variables(iEnergy) % Values(iRho,iTemp,iYp) = eps_baryons + eps_electrons
				EOSBaryonPlusEleTable % DV % Variables(iMue) % Values(iRho,iTemp,iYp) = mue_electrons
				EOSBaryonPlusEleTable % DV % Variables(iMup) % Values(iRho,iTemp,iYp) = mup_baryons
				EOSBaryonPlusEleTable % DV % Variables(iMun) % Values(iRho,iTemp,iYp) = mun_baryons
				EOSBaryonPlusEleTable % DV % Variables(iXp) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(iXp) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(iXn) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(iXn) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(iXa) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(iXa) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(iXh) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(iXh) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(iAbar) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(iAbar) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(iZbar) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(iZbar) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(iBE) % Values(iRho,iTemp,iYp) = & 
					EOSBaryonTable % DV % Variables(iBE) % Values(iRho,iTemp,iYp)
				EOSBaryonPlusEleTable % DV % Variables(ics2) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSBaryonPlusEleTable % DV % Variables(iGamma1) % Values(iRho,iTemp,iYp) = 0.0d0

				EOSEleTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) = press_electrons
				EOSEleTable % DV % Variables(iEntropy) % Values(iRho,iTemp,iYp) = s_electrons
				EOSEleTable % DV % Variables(iEnergy) % Values(iRho,iTemp,iYp) = eps_electrons
				EOSEleTable % DV % Variables(iMue) % Values(iRho,iTemp,iYp) = mue_electrons
				EOSEleTable % DV % Variables(iMup) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iMun) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iXp) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iXn) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iXa) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iXh) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iAbar) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iZbar) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iBE) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(ics2) % Values(iRho,iTemp,iYp) = 0.0d0
				EOSEleTable % DV % Variables(iGamma1) % Values(iRho,iTemp,iYp) = 0.0d0

			END DO
		END DO
	END DO
	
	ALLOCATE( eos_table(nRho, nTemp, nYp, 3) )
	eos_table(:,:,:,1) = LOG10(EOSBaryonPlusEleTable % DV % Variables(iPressure) % Values(:,:,:))
	eos_table(:,:,:,2) = EOSBaryonPlusEleTable % DV % Variables(iEntropy) % Values(:,:,:)
	energy_shift = -1.1d0 * MINVAL(EOSBaryonPlusEleTable % DV % Variables(iEnergy) % Values)
	eos_table(:,:,:,3) = LOG10(EOSBaryonPlusEleTable % DV % Variables(iEnergy) % Values(:,:,:) + energy_shift)
	
	CALL derivatives_production(igamma, nRho, nTemp, nYp, &
		LOG10(EOSBaryonPlusEleTable % TS % States(1) % Values(:)), &
		LOG10(EOSBaryonPlusEleTable % TS % States(2) % Values(:)), &
		EOSBaryonPlusEleTable % TS % States(3) % Values(:), &
		eos_table, energy_shift, &
		EOSBaryonPlusEleTable % DV % Variables(ics2) % Values, &
		EOSBaryonPlusEleTable % DV % Variables(iGamma1) % Values )

	eos_table(:,:,:,1) = LOG10(EOSEleTable % DV % Variables(iPressure) % Values(:,:,:))
	eos_table(:,:,:,2) = EOSEleTable % DV % Variables(iEntropy) % Values(:,:,:)
	energy_shift = -1.1d0 * MINVAL(EOSEleTable % DV % Variables(iEnergy) % Values)
	eos_table(:,:,:,3) = LOG10(EOSEleTable % DV % Variables(iEnergy) % Values(:,:,:) + energy_shift)
	
	CALL derivatives_production(igamma, nRho, nTemp, nYp, &
		LOG10(EOSEleTable % TS % States(1) % Values(:)), &
		LOG10(EOSEleTable % TS % States(2) % Values(:)), &
		EOSEleTable % TS % States(3) % Values(:), &
		eos_table, energy_shift, &
		EOSEleTable % DV % Variables(ics2) % Values, &
		EOSEleTable % DV % Variables(iGamma1) % Values )

	eos_table(:,:,:,1) = LOG10(EOSBaryonTable % DV % Variables(iPressure) % Values(:,:,:))
	eos_table(:,:,:,2) = EOSBaryonTable % DV % Variables(iEntropy) % Values(:,:,:)
	energy_shift = -1.1d0 * MINVAL(EOSBaryonTable % DV % Variables(iEnergy) % Values)
	eos_table(:,:,:,3) = LOG10(EOSBaryonTable % DV % Variables(iEnergy) % Values(:,:,:) + energy_shift)
	
	CALL derivatives_production(igamma, nRho, nTemp, nYp, &
		LOG10(EOSBaryonTable % TS % States(1) % Values(:)), &
		LOG10(EOSBaryonTable % TS % States(2) % Values(:)), &
		EOSBaryonTable % TS % States(3) % Values(:), &
		eos_table, energy_shift, &
		EOSBaryonTable % DV % Variables(ics2) % Values, &
		EOSBaryonTable % DV % Variables(iGamma1) % Values )

	! Do some housekeeping
	DO iRho=1,nRho
		DO iTemp=1,nTemp
			DO iYp=1,nYp	
				
				IF (EOSBaryonTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) .le. 0.0d0) THEN
					EOSBaryonTable % DV % Variables(ics2) % Values(iRho,iTemp,iYp) = 0.0d0
					EOSBaryonTable % DV % Variables(iGamma1) % Values(iRho,iTemp,iYp) = 0.0d0
				END IF

				IF (EOSEleTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) .le. 0.0d0) THEN
					EOSEleTable % DV % Variables(ics2) % Values(iRho,iTemp,iYp) = 0.0d0
					EOSEleTable % DV % Variables(iGamma1) % Values(iRho,iTemp,iYp) = 0.0d0
				END IF
				
				IF (EOSBaryonPlusEleTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) .le. 0.0d0) THEN
					EOSBaryonPlusEleTable % DV % Variables(ics2) % Values(iRho,iTemp,iYp) = 0.0d0
					EOSBaryonPlusEleTable % DV % Variables(iGamma1) % Values(iRho,iTemp,iYp) = 0.0d0
				END IF
				
			END DO
		END DO
	END DO

    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Baryonic + Electrons EOS"
    CALL WriteEquationOfStateTableHDF( EOSBaryonPlusEleTable, BaryonPlusEleName )
    CALL FinalizeHDF( )
	
	CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Electrons EOS"
    CALL WriteEquationOfStateTableHDF( EOSEleTable, ElectronEOSTableName )
    CALL FinalizeHDF( )

	DEALLOCATE(eos_table)
	
END PROGRAM wlCompareToFullEOS				

PROGRAM wlCompareToFullEOS
	
	USE wlKindModule, ONLY: dp
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
	USE wlEOSIOModuleHDF	
  USE wlLeptonEOSModule, ONLY: &
    HelmholtzTableType, MuonEOSType
  USE wlElectronPhotonEOS, ONLY: &
    ElectronPhotonEOS, ElectronPhotonStateType
  USE wlHelmholtzEOS, ONLY: &
    FullHelmEOS, HelmholtzStateType
  USE wlHelmMuonIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF
	USE wlExtPhysicalConstantsModule, ONLY: kmev, rmu, kmev_inv, ergmev, me
	USE wlExtEOSWrapperModule, ONLY : wlGetElectronEOS

	IMPLICIT NONE
	
	TYPE(ElectronPhotonStateType) :: ElectronPhotonState
	TYPE(HelmholtzStateType) :: HelmholtzState
	TYPE(HelmholtzTableType) :: HelmholtzTable
  TYPE(EquationOfStateTableType) :: EOSBaryonTable, EOSFullTable, EOSBaryonPlusEleTable, EOSEleTable
	
	INTEGER :: nRho, nTemp, nYp, iAbar, iZbar, iBE, iPressure, iRho, iTemp, iYp, iXp, iXa, iXn, iXh, &
		iEnergy, iEntropy, iMue, iMun, iMup, iRho_transition, ics2, iGamma1
	REAL(dp) :: abar_total, zbar_total, abar, zbar, heavy_frac, p_frac, n_frac, alpha_frac
	
	REAL(dp) :: press_e_BCK, entrop_e_BCK, energ_e_BCK, chem_e_BCK

  REAL(DP) :: press_electrons_helm
	REAL(dp) :: press_baryons, press_electrons, press_total
	REAL(dp) :: s_baryons, s_electrons, s_total
	REAL(dp) :: eps_baryons, eps_electrons, eps_total
	REAL(dp) :: mue_baryons, mue_electrons, mue_total
	REAL(dp) :: mup_baryons, mup_total, mun_baryons, mun_total
		
	LOGICAL :: RedHDF5Table
	
  CHARACTER(len=128) :: BaryonEOSTableName, FullEOSTableName

  RedHDF5Table = .false.
	
	IF (RedHDF5Table) THEN
		FullEOSTableName = 'FullEOS_interpolated.h5'
		BaryonEOSTableName = 'BaryonsPlusHelmPlusMuonsEOS_interpolated.h5'
	ELSE
		FullEOSTableName = 'FullEOS.h5'
		BaryonEOSTableName = 'BaryonsPlusHelmPlusMuonsEOS.h5'
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
				ElectronPhotonState % t = EOSBaryonTable % TS % States(2) % Values(iTemp)
				ElectronPhotonState % rho = EOSBaryonTable % TS % States(1) % Values(iRho)
				ElectronPhotonState % ye = EOSBaryonTable % TS % States(3) % Values(iYp)
				
        HelmholtzState % t = EOSBaryonTable % TS % States(2) % Values(iTemp)
				HelmholtzState % rho = EOSBaryonTable % TS % States(1) % Values(iRho)
				HelmholtzState % ye = EOSBaryonTable % TS % States(3) % Values(iYp)
				HelmholtzState % abar = abar_total
				HelmholtzState % zbar = zbar_total

				! calculate electron quantities
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)
        CALL FullHelmEOS(1, HelmholtzTable, HelmholtzState)
        CALL wlGetElectronEOS(ElectronPhotonState % rho, &
            ElectronPhotonState % t, ElectronPhotonState % ye, &
            press_e_BCK, entrop_e_BCK, energ_e_BCK, chem_e_BCK)

        entrop_e_BCK = entrop_e_BCK
				press_baryons = EOSBaryonTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iPressure)
				press_electrons = ElectronPhotonState % pele
        press_electrons_helm = HelmholtzState % pele

				s_baryons = EOSBaryonTable % DV % Variables(iEntropy) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iEntropy)
				s_electrons = ElectronPhotonState % sele
		
				eps_baryons = EOSBaryonTable % DV % Variables(iEnergy) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iEnergy)
				eps_electrons = ElectronPhotonState % eele

				mue_baryons = EOSBaryonTable % DV % Variables(iMue) % Values(iRho,iTemp,iYp) - &
				EOSBaryonTable % DV % Offsets(iMue)
				mue_electrons = ElectronPhotonState % mue

        IF ( ( ABS(entrop_e_BCK - s_electrons)/entrop_e_BCK > 3.0d-2 ) &
             .AND. (HelmholtzState % rho < 1.0d15) ) THEN
          WRITE(*,*) iRho, iTemp, iYp, &
          entrop_e_BCK, s_electrons, ABS(entrop_e_BCK - s_electrons)/entrop_e_BCK
        ENDIF

        ! IF ( (iRho == 1) .AND. (iTemp == 76) .AND. (iYp == 49) ) THEN
        !   WRITE(*,*) iRho, iTemp, iYp, &
        !   entrop_e_BCK, s_electrons, ABS(entrop_e_BCK - s_electrons)/entrop_e_BCK
        ! ENDIF

			END DO
		END DO
	END DO
		
END PROGRAM wlCompareToFullEOS				

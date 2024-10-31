PROGRAM wlCalculateSoundSpeed
	
	USE wlKindModule, ONLY: dp
	USE wlExtNumericalModule, ONLY: zero, one, pi
    USE wlEquationOfStateTableModule
    USE wlIOModuleHDF
	USE wlEOSIOModuleHDF	
	USE wlLeptonEOSModule, ONLY: &
		HelmholtzEOSType, MuonEOSType
	USE wlElectronEOS, ONLY: &
		FullHelmEOS, MinimalHelmEOS_rt, ElectronStateType
	USE wlMuonEOS, ONLY: &
		FullMuonEOS, MuonStateType
	USE wlHelmMuonIOModuleHDF, ONLY: &
		ReadHelmholtzTableHDF, ReadMuonTableHDF
	USE wlExtPhysicalConstantsModule, ONLY: kmev, rmu, kmev_inv, ergmev, me, cvel
	USE wlGammaSoundSpeed
	USE wlBolligSoundSpeed
	
	IMPLICIT NONE
	
	TYPE(HelmholtzEOSType) :: HelmholtzTable
	TYPE(ElectronStateType) :: ElectronState
	TYPE(MuonEOSType) :: MuonTable
	TYPE(MuonStateType) :: MuonState
    TYPE(EquationOfStateTableType) :: EOSBaryonTable, EOSFullTable, EOSBaryonPlusEleTable
	
	INTEGER :: nRho, nTemp, nYp, iAbar, iZbar, iBE, iPressure, iRho, iTemp, iYp, iXp, iXa, iXn, iXh, &
		iEnergy, iEntropy, iMue, iMun, iMup, iRho_transition, ics2, iGamma1
	REAL(dp) :: abar_total, zbar_total, abar, zbar, heavy_frac, p_frac, n_frac, alpha_frac
	
	REAL(dp) :: press_baryons, press_electrons, press_total
	REAL(dp) :: s_baryons, s_electrons, s_total
	REAL(dp) :: eps_baryons, eps_electrons, eps_total
	REAL(dp) :: mue_baryons, mue_electrons, mue_total
	REAL(dp) :: mup_baryons, mup_total, mun_baryons, mun_total
	
	REAL(DP) :: OS_P, OS_E, OS_S, OS_V
	REAL(DP) :: D, T, Yp, Ye, Ymu, cs2, Gamma
	REAL(dp), ALLOCATABLE :: eos_table(:,:,:,:), Pe(:,:,:), Ee(:,:,:), &
		P_T(:,:,:), V_T(:,:,:), D_T(:), T_T(:), Yp_T(:)
	
	LOGICAL :: RedHDF5Table
	INTEGER :: igamma = 3
	INTEGER :: iP_bary, iS_bary, iE_bary
	
	REAL(DP) :: tBegin, tEnd
	
    CHARACTER(len=128) :: BaryonEOSTableName, FullEOSTableName, BaryonPlusEleName

    RedHDF5Table = .true.
	
	IF (RedHDF5Table) THEN
		FullEOSTableName = '/mnt/c/Users/Lucab/GDrive/Research/Muons_project/weaklib_tables/FullEOS_interpolated.h5'
		BaryonEOSTableName = '/mnt/c/Users/Lucab/GDrive/Research/Muons_project/weaklib_tables/BaryonsPlusHelmPlusMuonsEOS_interpolated.h5'
	ELSE
		FullEOSTableName = '/mnt/c/Users/Lucab/GDrive/Research/Muons_project/weaklib_tables/FullEOS.h5'
		BaryonEOSTableName = '/mnt/c/Users/Lucab/GDrive/Research/Muons_project/weaklib_tables/BaryonsPlusHelmPlusMuonsEOS.h5'
	ENDIF
	
	BaryonPlusEleName = 'csBaryonPlusEle.h5'
	
	! read in helmholtz table
	CALL ReadHelmholtzTableHDF( HelmholtzTable, BaryonEOSTableName )
	
	! read in baryon table -------------------------------
    CALL ReadEquationOfStateTableHDF( EOSBaryonTable, BaryonEOSTableName )

	! read in muon table
	CALL ReadMuonTableHDF( MuonTable, BaryonEOSTableName )

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
	
	iP_bary = EOSBaryonTable % DV % Indices % iPressure
	iS_bary = EOSBaryonTable % DV % Indices % iEntropyPerBaryon
	iE_bary = EOSBaryonTable % DV % Indices % iInternalEnergyDensitY

	OS_P = EOSBaryonTable % DV % Offsets(iP_bary)
	OS_S = EOSBaryonTable % DV % Offsets(iS_bary)
	OS_E = EOSBaryonTable % DV % Offsets(iE_bary)
	
	ALLOCATE( Pe(nRho, nTemp, nYp) )
	ALLOCATE( Ee(nRho, nTemp, nYp) )
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
				ElectronState % y_e = EOSBaryonTable % TS % States(3) % Values(iYp)
				
				! calculate electron quantities
				CALL FullHelmEOS(1, HelmholtzTable, ElectronState, .false., .false.)

				Ee(iRho,iTemp,iYp) = ElectronState % e + 0.511d0 / rmu * ergmev * ElectronState % y_e
				Pe(iRho,iTemp,iYp) = ElectronState % p
				
				press_baryons = 10.0d0**(EOSBaryonTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp)) - &
					EOSBaryonTable % DV % Offsets(iPressure)
				press_electrons = ElectronState % p
				press_total = 10.0d0**(EOSFullTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp)) - &
					EOSFullTable % DV % Offsets(iPressure)
				
				s_baryons = 10.0d0**(EOSBaryonTable % DV % Variables(iEntropy) % Values(iRho,iTemp,iYp)) - &
					EOSBaryonTable % DV % Offsets(iEntropy)
				s_electrons = ElectronState % s * rmu * kmev_inv / ergmev
				s_total = 10.0d0**(EOSFullTable % DV % Variables(iEntropy) % Values(iRho,iTemp,iYp)) - &
					EOSFullTable % DV % Offsets(iEntropy)
				
				eps_baryons = 10.0d0**(EOSBaryonTable % DV % Variables(iEnergy) % Values(iRho,iTemp,iYp)) - &
					EOSBaryonTable % DV % Offsets(iEnergy)
				eps_electrons = ElectronState % e + 0.511 / rmu * ergmev * ElectronState % y_e
				eps_total = 10.0d0**(EOSFullTable % DV % Variables(iEnergy) % Values(iRho,iTemp,iYp)) - &
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

			END DO
		END DO
	END DO
	
	ALLOCATE( eos_table(nRho, nTemp, nYp, 3) )
	ALLOCATE( P_T(nRho, nTemp, nYp) )
	ALLOCATE( V_T(nRho, nTemp, nYp) )
	ALLOCATE( D_T(nRho) )
	ALLOCATE( T_T(nTemp) )
	ALLOCATE( Yp_T(nYp) )

	WRITE(*,*) 'start stcol'
	CALL CPU_TIME( tBegin )
	! SOME DEFINITION
	CALL CalculateSeparateSoundSpeed(nRho, nTemp, nYp, &
		10.0d0**(EOSBaryonTable % DV % Variables(iPressure) % Values(:,:,:)) - OS_P, Pe, &
		10.0d0**(EOSBaryonTable % DV % Variables(iEnergy) % Values(:,:,:)) - OS_E, Ee, &
		OS_E, &
		EOSBaryonTable % TS % States(1) % Values(:),  &
		EOSBaryonTable % TS % States(2) % Values(:),  &
		EOSBaryonTable % TS % States(3) % Values(:),  &
		EOSBaryonPlusEleTable % DV % Variables(iGamma1) % Values, &
		EOSBaryonPlusEleTable % DV % Variables(ics2) % Values )
	CALL CPU_TIME( tEnd )
	WRITE(*,*) 'end stcol:', tEnd - tBegin 

	! ! ANOTHER DEFINITION BUT THIS ONE IS FOR THE FULL TABLE SO PROB NOT NEEDED
	! eos_table(:,:,:,1) = EOSBaryonPlusEleTable % DV % Variables(iPressure) % Values(:,:,:)
	! eos_table(:,:,:,2) = EOSBaryonPlusEleTable % DV % Variables(iEntropy) % Values(:,:,:)
	! eos_table(:,:,:,3) = EOSBaryonPlusEleTable % DV % Variables(iEnergy) % Values(:,:,:)
	
	! CALL derivatives_production(igamma, nRho, nTemp, nYp, &
		! LOG10(EOSBaryonPlusEleTable % TS % States(1) % Values(:)), &
		! LOG10(EOSBaryonPlusEleTable % TS % States(2) % Values(:)), &
		! EOSBaryonPlusEleTable % TS % States(3) % Values(:), &
		! eos_table, OS_E, &
		! EOSBaryonPlusEleTable % DV % Variables(ics2) % Values, &
		! EOSBaryonPlusEleTable % DV % Variables(iGamma1) % Values )

	D_T = EOSBaryonTable % TS % States(1) % Values(:)
	T_T = EOSBaryonTable % TS % States(2) % Values(:)
	Yp_T = EOSBaryonTable % TS % States(3) % Values(:)
	P_T = EOSBaryonTable % DV % Variables(iPressure) % Values(:,:,:)
	V_T = EOSBaryonTable % DV % Variables(iEnergy) % Values(:,:,:)
	OS_V = OS_E
	
	! Do some housekeeping
	WRITE(*,*) 'start bollig'
	CALL CPU_TIME( tBegin )
	DO iRho=1,nRho
		DO iTemp=1,nTemp
			DO iYp=1,nYp	
				
				! BOLLIG DEFINITION
				D = EOSBaryonTable % TS % States(1) % Values(iRho)
				T = EOSBaryonTable % TS % States(2) % Values(iTemp)
				Yp = EOSBaryonTable % TS % States(3) % Values(iYp)
				
				Ye = Yp
				Ymu = 0.0d0
				
				CALL CalculateBolligSoundSpeed( D, T, Yp, Ye, Ymu, D_T, T_T, Yp_T, P_T, OS_P, V_T, OS_V, &
					'Energy', HelmholtzTable, MuonTable, Gamma, cs2)
								
				EOSBaryonPlusEleTable % DV % Variables(ics2) % Values(iRho,iTemp,iYp) = cs2
				EOSBaryonPlusEleTable % DV % Variables(iGamma1) % Values(iRho,iTemp,iYp) = Gamma			
				
				IF (EOSBaryonPlusEleTable % DV % Variables(iPressure) % Values(iRho,iTemp,iYp) .le. 0.0d0) THEN
					EOSBaryonPlusEleTable % DV % Variables(ics2) % Values(iRho,iTemp,iYp) = 0.0d0
					EOSBaryonPlusEleTable % DV % Variables(iGamma1) % Values(iRho,iTemp,iYp) = 0.0d0
				END IF
				
			END DO
		END DO
	END DO
	CALL CPU_TIME( tEnd )
	WRITE(*,*) 'end bollig:', tEnd - tBegin 

    CALL InitializeHDF( )
    WRITE (*,*) "Starting HDF write: Baryonic + Electrons EOS"
    CALL WriteEquationOfStateTableHDF( EOSBaryonPlusEleTable, BaryonPlusEleName )
    CALL FinalizeHDF( )

	DEALLOCATE(eos_table)
	
CONTAINS

	SUBROUTINE CalculateSeparateSoundSpeed(nrho, ntemp, nye, Pb_in, Pe_in, Eb_in, Ee_in, Eb_shift, rho, temp, ye, Gamma, cs2)
	
		INTEGER, INTENT(IN) :: nrho, ntemp, nye
		REAL(DP), INTENT(IN) :: Pb_in(nrho,ntemp,nye), Pe_in(nrho,ntemp,nye), Eb_in(nrho,ntemp,nye), Ee_in(nrho,ntemp,nye)
		
		REAL(DP) :: Eb_shift
		REAL(DP), INTENT(OUT) :: Gamma(nrho,ntemp,nye), cs2(nrho,ntemp,nye)
		REAL(DP) :: Pb(nrho,ntemp,nye), Pe(nrho,ntemp,nye)
		REAL(DP) :: Eb(nrho,ntemp,nye), Ee(nrho,ntemp,nye)
		REAL(DP) :: dlnEbdlnT(nrho,ntemp,nye), dlnEedlnT(nrho,ntemp,nye)
		REAL(DP) :: dlnPbdlnT(nrho,ntemp,nye), dlnPedlnT(nrho,ntemp,nye)
		REAL(DP) :: dlnPbdlnrho(nrho,ntemp,nye), dlnPedlnrho(nrho,ntemp,nye)
		REAL(DP) :: rho(nrho), temp(ntemp), ye(nye)
	
		integer :: i,j,k
		REAL(DP) :: f1, f2, x1, x2
		
		do k=1,nye
			do i=1,nrho
				do j=1,ntemp
				
					if (Pb_in(i,j,k) .lt. 0.0d0) then
						Pb(i,j,k) = 1.0d-30
					else
						Pb(i,j,k) = Pb_in(i,j,k)
					endif
					if (Pe_in(i,j,k) .lt. 0.0d0) then
						Pe(i,j,k) = 1.0d-30
					else
						Pe(i,j,k) = Pe_in(i,j,k)
					endif					
					
				enddo
			enddo
		enddo
		
		Eb = Eb_in
		Ee = Ee_in
		
		do k=1,nye
			do i=1,nrho
				! Energy
				do j=2,ntemp-1
					x1 = LOG10(temp(j-1))
					f1 = LOG10(Ee(i,j-1,k))
					x2 = LOG10(temp(j+1))
					f2 = LOG10(Ee(i,j+1,k))
					dlnEedlnT(i,j,k) = (f2-f1)/(x2-x1)
				enddo
				
				! boundaries: one-sided derivative
				x2 = LOG10(temp(2))
				x1 = LOG10(temp(1))
				f2 = LOG10(Ee(i,2,k))
				f1 = LOG10(Ee(i,1,k))
				dlnEedlnT(i,1,k) = (f2-f1)/(x2-x1)

				x2 = LOG10(temp(ntemp))
				x1 = LOG10(temp(ntemp-1))
				f2 = LOG10(Ee(i,ntemp,k))
				f1 = LOG10(Ee(i,ntemp-1,k))
				dlnEedlnT(i,1,k) = (f2-f1)/(x2-x1)
				
				do j=2,ntemp-1
					x1 = LOG10(temp(j-1))
					f1 = LOG10(Eb(i,j-1,k))
					x2 = LOG10(temp(j+1))
					f2 = LOG10(Eb(i,j+1,k))
					dlnEbdlnT(i,j,k) = (f2-f1)/(x2-x1)
				enddo
				
				! boundaries: one-sided derivative
				x2 = LOG10(temp(2))
				x1 = LOG10(temp(1))
				f2 = LOG10(Eb(i,2,k))
				f1 = LOG10(Eb(i,1,k))
				dlnEbdlnT(i,1,k) = (f2-f1)/(x2-x1) 

				x2 = LOG10(temp(ntemp))
				x1 = LOG10(temp(ntemp-1))
				f2 = LOG10(Eb(i,ntemp,k))
				f1 = LOG10(Eb(i,ntemp-1,k))
				dlnEbdlnT(i,1,k) = (f2-f1)/(x2-x1)
				
				! Pressure
				do j=2,ntemp-1
					x1 = LOG10(temp(j-1))
					f1 = LOG10(Pe(i,j-1,k))
					x2 = LOG10(temp(j+1))
					f2 = LOG10(Pe(i,j+1,k))
					dlnPedlnT(i,j,k) = (f2-f1)/(x2-x1)
				enddo
				
				! boundaries: one-sided derivative
				x2 = LOG10(temp(2))
				x1 = LOG10(temp(1))
				f2 = LOG10(Pe(i,2,k))
				f1 = LOG10(Pe(i,1,k))
				dlnPedlnT(i,1,k) = (f2-f1)/(x2-x1)

				x2 = LOG10(temp(ntemp))
				x1 = LOG10(temp(ntemp-1))
				f2 = LOG10(Pe(i,ntemp,k))
				f1 = LOG10(Pe(i,ntemp-1,k))
				dlnPedlnT(i,1,k) = (f2-f1)/(x2-x1)
				
				do j=2,ntemp-1
					x1 = LOG10(temp(j-1))
					f1 = LOG10(Pb(i,j-1,k))
					x2 = LOG10(temp(j+1))
					f2 = LOG10(Pb(i,j+1,k))
					dlnPbdlnT(i,j,k) = (f2-f1)/(x2-x1)
				enddo
				
				! boundaries: one-sided derivative
				x2 = LOG10(temp(2))
				x1 = LOG10(temp(1))
				f2 = LOG10(Pb(i,2,k))
				f1 = LOG10(Pb(i,1,k))
				dlnPbdlnT(i,1,k) = (f2-f1)/(x2-x1)

				x2 = LOG10(temp(ntemp))
				x1 = LOG10(temp(ntemp-1))
				f2 = LOG10(Pb(i,ntemp,k))
				f1 = LOG10(Pb(i,ntemp-1,k))
				dlnPbdlnT(i,1,k) = (f2-f1)/(x2-x1)	
			enddo
		enddo		
	
		do k=1,nye
			do j=1,ntemp
				do i=2,nrho-1
					x1 = LOG10(rho(i-1))
					x2 = LOG10(rho(i+1))
					f1 = LOG10(Pe(i-1,j,k))
					f2 = LOG10(Pe(i+1,j,k))
					dlnPedlnrho(i,j,k) = (f2-f1)/(x2-x1)
				enddo
				
				! boundaries: one-sided derivative
				x1 = LOG10(rho(1))
				x2 = LOG10(rho(2))
				f1 = LOG10(Pe(1,j,k))
				f2 = LOG10(Pe(2,j,k))
				dlnPedlnrho(i,j,k) = (f2-f1)/(x2-x1)

				x1 = LOG10(rho(nrho-1))
				x2 = LOG10(rho(nrho))
				f1 = LOG10(Pe(nrho-1,j,k))
				f2 = LOG10(Pe(nrho,j,k))
				dlnPedlnrho(i,j,k) = (f2-f1)/(x2-x1)

				do i=2,nrho-1
					x1 = LOG10(rho(i-1))
					x2 = LOG10(rho(i+1))
					f1 = LOG10(Pb(i-1,j,k))
					f2 = LOG10(Pb(i+1,j,k))
					dlnPbdlnrho(i,j,k) = (f2-f1)/(x2-x1)
				enddo
				
				! boundaries: one-sided derivative
				x1 = LOG10(rho(1))
				x2 = LOG10(rho(2))
				f1 = LOG10(Pb(1,j,k))
				f2 = LOG10(Pb(2,j,k))
				dlnPbdlnrho(i,j,k) = (f2-f1)/(x2-x1)

				x1 = LOG10(rho(nrho-1))
				x2 = LOG10(rho(nrho))
				f1 = LOG10(Pb(nrho-1,j,k))
				f2 = LOG10(Pb(nrho,j,k))
				dlnPbdlnrho(i,j,k) = (f2-f1)/(x2-x1)
			enddo
		enddo
	
		! Handle special cases where the baryonic pressures and/or energies are negative 
		do k=1,nye
			do j=1,ntemp
				do i=1,nrho		
					IF (Pb(i,j,k) == 0.0d0) THEN
						dlnPbdlnrho(i,j,k) = 0.0d0
						dlnPbdlnT(i,j,k) = 0.0d0
					ENDIF

					IF (Eb(i,j,k) == 0.0d0) THEN
						dlnEbdlnT(i,j,k) = 0.0d0
					ENDIF

				enddo
			enddo
		enddo	
		
		do k=1,nye
			do j=1,ntemp
				do i=1,nrho		
					Gamma(i,j,k) = (Pb(i,j,k) * dlnPbdlnrho(i,j,k) + Pe(i,j,k) * dlnPedlnrho(i,j,k) + &
						(Pb(i,j,k) * dlnPbdlnT(i,j,k) + Pe(i,j,k) * dlnPedlnT(i,j,k))**2.0d0 / &
						((Eb(i,j,k) * dlnEbdlnT(i,j,k) + Ee(i,j,k) * dlnEedlnT(i,j,k)) * rho(i) ) ) / &
						(Pb(i,j,k) + Pe(i,j,k))
	
					cs2(i,j,k) = Gamma(i,j,k) * (Pb(i,j,k) + Pe(i,j,k)) / rho(i) 
					cs2(i,j,k) = Gamma(i,j,k) * (Pb(i,j,k) + Pe(i,j,k)) / &
					(Pb(i,j,k) + Pe(i,j,k) + (Eb(i,j,k) - Eb_shift + Ee(i,j,k) - &
					0.511d0 / rmu * ergmev * ye(k) + cvel**2) * rho(i) ) * cvel**2
										
				enddo
				! i = 160
				! WRITE(*,*) Pb(i,j,k) * dlnPbdlnrho(i,j,k) + Pe(i,j,k) * dlnPedlnrho(i,j,k), &
				! (Pb(i,j,k) * dlnPbdlnT(i,j,k) + Pe(i,j,k) * dlnPedlnT(i,j,k))**2.0d0, &
				! (Eb(i,j,k) * dlnEbdlnT(i,j,k) + Ee(i,j,k) * dlnEedlnT(i,j,k)) * rho(i) , &
				! Gamma(i,j,k), cs2(i,j,k)

				! WRITE(*,*) cs2(i,j,k), Gamma(i,j,k), Eb(i,j,k) , dlnEbdlnT(i,j,k), Ee(i,j,k) , dlnEedlnT(i,j,k), Eb_shift

			enddo
		enddo
		! WRITE(*,*) MAXVAL(cs2), MINVAL(cs2)
		! WRITE(*,*) MAXVAL(Gamma), MINVAL(Gamma)
		! WRITE(*,*)
		! WRITE(*,*) MAXVAL(Pb), MAXVAL(Pe)
		! WRITE(*,*) MINVAL(Pb), MINVAL(Pe)
		! WRITE(*,*)
		! WRITE(*,*) MAXVAL(Eb), MAXVAL(Ee)
		! WRITE(*,*) MINVAL(Eb), MINVAL(Ee)
		! WRITE(*,*)
		
	END SUBROUTINE CalculateSeparateSoundSpeed
	
END PROGRAM wlCalculateSoundSpeed				

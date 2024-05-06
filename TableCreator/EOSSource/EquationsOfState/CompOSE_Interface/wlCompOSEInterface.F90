MODULE wlCompOSEInterface

  USE wlKindModule, ONLY: dp
  USE wlExtPhysicalConstantsModule, ONLY: cvel, ergmev, cm3fm3, kmev_inv, rmu
  USE wlExtNumericalModule, ONLY: zero, half, one, pi

  IMPLICIT NONE
  PRIVATE
	
  PUBLIC :: ReadnPointsFromCompOSE, ReadCompOSETable
	
  
  REAL(dp), PUBLIC, DIMENSION(:,:,:,:), ALLOCATABLE :: EOSCompOSE

  ! output baryon density (fm^-3), temperature (kelvin), and proton (or charge) fraction
  REAL(dp), DIMENSION(:), ALLOCATABLE :: nbCompOSE
  REAL(dp), PUBLIC, DIMENSION(:), ALLOCATABLE :: RhoCompOSE
  REAL(dp), PUBLIC, DIMENSION(:), ALLOCATABLE :: TempCompOSE
  REAL(dp), PUBLIC, DIMENSION(:), ALLOCATABLE :: YpCompOSE

  INTEGER, PARAMETER :: nVariablesCompOSE = 15
  INTEGER, PARAMETER :: iPressCompOSE = 1
  INTEGER, PARAMETER :: iEntropyCompOSE = 2
  INTEGER, PARAMETER :: iInternalEnergyDensityCompOSE = 3
  INTEGER, PARAMETER :: iElectronChemPotCompOSE = 4
  INTEGER, PARAMETER :: iProtonChemPotCompOSE = 5
  INTEGER, PARAMETER :: iNeutronChemPotCompOSE = 6
  INTEGER, PARAMETER :: iProtonMassFractionCompOSE = 7
  INTEGER, PARAMETER :: iNeutronMassFractionCompOSE = 8
  INTEGER, PARAMETER :: iAlphaMassFractionCompOSE = 9
  INTEGER, PARAMETER :: iHeavyMassFractionCompOSE = 10
  INTEGER, PARAMETER :: iHeavyChargeNumberCompOSE = 11
  INTEGER, PARAMETER :: iHeavyMassNumberCompOSE = 12
  INTEGER, PARAMETER :: iHeavyBindingEnergyCompOSE = 13
  INTEGER, PARAMETER :: iThermalEnergyCompOSE = 14
  INTEGER, PARAMETER :: iGamma1CompOSE = 15

CONTAINS

  SUBROUTINE ReadnPointsFromCompOSE( CompOSEFilePath, nRho, nTemp, nYp )
	
	! input path of the directory with the CompOSE output
	CHARACTER(len=128), INTENT(IN) :: CompOSEFilePath
	
	! Output variables
    INTEGER, INTENT(INOUT) :: nRho, nTemp, nYp
	
	
	! get density grid
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.nb', & 
		FORM = "formatted", ACTION = 'read')
	READ(123,*) 
	READ(123,*) nRho

	CLOSE(123)
	
	! get temperature grid
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.t', & 
		FORM = "formatted", ACTION = 'read')
	READ(123,*)
	READ(123,*) nTemp

	CLOSE(123)
	
	! get yp grid
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.yq', & 
		FORM = "formatted", ACTION = 'read')
	READ(123,*)
	READ(123,*) nYp

	CLOSE(123)
		
  END SUBROUTINE ReadnPointsFromCompOSE

  SUBROUTINE ReadCompOSETable( CompOSEFilePath, nRho, nTemp, nYp )
	
	! input path of the directory with the CompOSE output
	CHARACTER(len=128), INTENT(IN) :: CompOSEFilePath
	
	! output number of points in the table
    INTEGER, INTENT(IN) :: nRho, nTemp, nYp
	
	! Local Variables
    INTEGER :: iLepton, iRho, iT, iYp, i_tot
	REAL(dp) :: NeutronMass, ProtonMass
	REAL(dp) :: Q1, Q2, Q3, Q4, Q5, Q6, Q7, TotalInternalEnergy
	REAL(dp) :: Xp, Xn, Xa, Xh, Abar, Zbar
		
	! Get proton and neutron mass
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.thermo', & 
	FORM = "formatted", ACTION = 'read')
	READ(123,*) NeutronMass, ProtonMass, iLepton
	CLOSE(123)
	
	! In the future you will want to catch whether you are reading in leptons or NOT!!!!!
	!IF (iLepton .ne. zero) THEN
	!  WRITE(*,*) 'The table you are trying to read already contains leptons!'
	!  STOP
	!END IF
				
	! Get density. DECIDE HOW TO CALCULATE THE MASS DENSITY SINCE IT IS NOT OBVIOUS
	ALLOCATE( nbCompOSE( nRho ) )
	ALLOCATE( RhoCompOSE( nRho ) )
	
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.nb', & 
		FORM = "formatted", ACTION = 'read')
	READ(123,*)
	READ(123,*)
	DO iRho=1,nRho
	   READ(123,*) nbCompOSE(iRho)
	END DO
	CLOSE(123)	
	! Convert to grams per cm^3
	RhoCompOSE(:) = nbCompOSE(:)*rmu/cm3fm3
	
	! get temperature
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.t', & 
		FORM = "formatted", ACTION = 'read')
	READ(123,*)
	READ(123,*)
	
	ALLOCATE( TempCompOSE( nTemp ) )
	
	DO iT=1,nTemp
	   READ(123,*) TempCompOSE(iT)
	END DO
		
	CLOSE(123)
	! Convert to kelvin
	TempCompOSE(:) = TempCompOSE(:)*kmev_inv

	! get yp grid
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.yq', & 
		FORM = "formatted", ACTION = 'read')
	READ(123,*)
	READ(123,*)
	
	ALLOCATE( YpCompOSE( nYp ) )
	
	DO iYp=1,nYp
	   READ(123,*) YpCompOSE(iYp)
	END DO

	CLOSE(123)
	
	! Allocate table
	ALLOCATE( EOSCompOSE(nRho,nTemp,nYp,nVariablesCompOSE) )
	
	! get thermal state
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.thermo', & 
		FORM = "formatted", ACTION = 'read')
	READ(123,*) NeutronMass, ProtonMass, iLepton

	! Read CompOSE EOS and convert to weaklib units
	! for the internal energy you have to be careful. YOu have to convert baryon density to mass density. It seems to me that
	! the best way of doing this is to set a fixed conversion factor, i.e. the nucleon mass. Notice that Q7 technically uses the 
	! neutron mass, so we have to understand if we like that or not
    DO i_tot=1,nTemp*nRho*nYp
	
	  READ(123,*) iT, iRho, iYp, Q1, Q2, Q3, Q4, Q5, Q6, Q7
	  
	  EOSCompOSE(iRho,iT,iYp,iPressCompOSE) = Q1 * nbCompOSE(iRho) * ergmev / cm3fm3
	  EOSCompOSE(iRho,iT,iYp,iEntropyCompOSE) = Q2
	  TotalInternalEnergy = (one + Q7) * nbCompOSE(iRho) * NeutronMass * ergmev / cm3fm3
	  
	  EOSCompOSE(iRho,iT,iYp,iInternalEnergyDensityCompOSE) = TotalInternalEnergy / RhoCompOSE(iRho) - cvel**2.0d0
	  ! gamma is defined as P/(rho eps), but maybe eps should be the thermal energy? I don't think so
	  
	  ! now handle the chemical potentials
	  EOSCompOSE(iRho,iT,iYp,iProtonChemPotCompOSE) = (Q3 + Q4 + one)*NeutronMass
	  EOSCompOSE(iRho,iT,iYp,iNeutronChemPotCompOSE) = (Q3 + one)*NeutronMass
	  EOSCompOSE(iRho,iT,iYp,iElectronChemPotCompOSE) = Q5*NeutronMass - Q4*NeutronMass
	  
	END DO
	
	CLOSE(123)
	
	! get composition from custom made table (see python script to convert eos.compo into eos.compo.WL
	OPEN(123, FILE=trim(adjustl(CompOSEFilePath))//'/eos.compo.wl', & 
		FORM = "formatted", ACTION = 'read')

	! Read CompOSE EOS and convert to weaklib units
	! for the internal energy you have to be careful. YOu have to convert baryon density to mass density. It seems to me that
	! the best way of doing this is to set a fixed conversion factor, i.e. the nucleon mass. Notice that Q7 technically uses the 
	! neutron mass, so we have to understand if we like that or not
    DO i_tot=1,nTemp*nRho*nYp
	
	  READ(123,*) iT, iRho, iYp, Xp, Xn, Xa, Xh, Abar, Zbar
	  
	  EOSCompOSE(iRho,iT,iYp,iProtonMassFractionCompOSE) = Xp
	  EOSCompOSE(iRho,iT,iYp,iNeutronMassFractionCompOSE) = Xn
	  EOSCompOSE(iRho,iT,iYp,iAlphaMassFractionCompOSE) = Xa
	  EOSCompOSE(iRho,iT,iYp,iHeavyMassFractionCompOSE) = Xh
	  EOSCompOSE(iRho,iT,iYp,iHeavyMassNumberCompOSE) = Abar
	  EOSCompOSE(iRho,iT,iYp,iHeavyChargeNumberCompOSE) = Zbar
	  
	  ! Finally set other things to zero for now, CompOSE does not tell you what they are
	  EOSCompOSE(iRho,iT,iYp,iHeavyBindingEnergyCompOSE) = zero
	  EOSCompOSE(iRho,iT,iYp,iThermalEnergyCompOSE) = zero
	  EOSCompOSE(iRho,iT,iYp,iGamma1CompOSE) = zero
	  
	END DO

  END SUBROUTINE ReadCompOSETable

END MODULE wlCompOSEInterface

MODULE wlOpacityTableModule
!-----------------------------------------------------------------------
!
!    File:         wlOpacityTableModule.f90
!    Module:       wlOpacityTableModule
!    Type:         Module w/ Subroutines
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      2/22/16
!    WeakLib ver:  
!
!    Purpose:
!      Allocate/DeAllocate table for opacity considered EoS table.
!
!   
!
!    CONTAINS:
!    
!
!    Modules used:
!       wlKindModule
!       HDF5
!       wlThermoStateModule
!       wlDependentVariablesModule
!       wlGridModule, ONLY: makeEnergyGrid EnergyGridType 
!  
!-----------------------------------------------------------------------
!  NOTE: Only Type A interaction applied. Type B and Type C interaction 
!        needs to be added for future use.
!-----------------------------------------------------------------------
!                         Three Opacity Type
!
! OpacityType A for  ABEM( rho, T, Ye, E)
!
!                e- + p/A <--> v_e + n/A*
!
! OpacityType B for  ISO( e, rho, T, Ye, l)
!                        
!                v_i/anti(v_i) + A --> v_i/anti(v_i) + A 
!                v_i/anti(v_i) + e+/e-/n/p  <-->  v_i/anti(v_i) + e+/e-/n/p
!
! OpacityType C for  NISO( e_in, e_out, rho, T, Ye, l)
!
!                e+ + e-  <--> v_i + anti(v_i);   i=e, muon, tau
!                N + N   <--> N + N + v_i + anti(v_i)
!                 
!-----------------------------------------------------------------------
 
  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlGridModule, ONLY: EnergyGridType 
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF, ONLY: ReadEquationOfStateTableHDF

  implicit none
  PRIVATE
   
!---------------------------------------------
!
! OpacityType A for  ABEM( rho, T, Ye, E)
!
! OpacityType B for  ISO( e, rho, T, Ye, l)
!                        
! OpacityType C for  NISO( e_in, e_out, rho, T, Ye, l)
!         
!---------------------------------------------
  TYPE, PUBLIC :: OpacityTypeA  
    CHARACTER(LEN=32) :: Names
    CHARACTER(LEN=32) :: Units
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: Values 
  END TYPE 

  TYPE, PUBLIC :: OpacityTypeB  
    CHARACTER(LEN=32) :: Names
    CHARACTER(LEN=32) :: Units
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: Values 
  END TYPE 

  TYPE, PUBLIC :: OpacityTypeC  
    CHARACTER(LEN=32) :: Names
    CHARACTER(LEN=32) :: Units
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: Values 
  END TYPE 

!=======================================
! OpacityTableType
!=======================================

  TYPE, PUBLIC :: OpacityTableType
  
    TYPE(OpacityTypeA), ALLOCATABLE, DIMENSION(:) :: ECAPEM 
    TYPE(EnergyGridType)                          :: EnergyGrid
    TYPE(EquationOfStateTableType)       :: EOSTable           
  END TYPE 

 
   PUBLIC AllocateOpacityTable 
   PUBLIC DeAllocateOpacityTable
   PUBLIC totalECapEm

CONTAINS


!==========================================================================
! Public Subroutine for OpacityTable
!==========================================================================
 SUBROUTINE AllocateOpacityTable( OpacityTable, nSpeciesA, nPointsE ) ! 'EOSTable' refers to the EOSTable we want to connect

     TYPE(OpacityTableType), INTENT(inout)  :: OpacityTable
     INTEGER, INTENT(in)                    :: nPointsE
     INTEGER, INTENT(in)                    :: nSpeciesA

  WRITE(*,*) "Readin EOSTable" ! make the readin EOSTable in readOpacityTable surbroutine later on

   CALL InitializeHDF( )
   CALL ReadEquationOfStateTableHDF( OpacityTable % EOSTable, &
                                  & "EquationOfStateTable.h5" )

  WRITE(*,*) 'Allocating OpacityTypeA' 
     CALL AllocateOpacityTypeA &
          ( OpacityTable % ECAPEM, nSpeciesA, &
           &OpacityTable % EOSTable % nPoints, nPointsE )  ! just for TypeA
  WRITE(*,*) 'Allocating EnergyGrid' 
     CALL AllocateEnergyGrid &
          ( OpacityTable % EnergyGrid, nPointsE)

  WRITE(*,*) 'Allocation Complete'

   
  END SUBROUTINE AllocateOpacityTable


  SUBROUTINE DeAllocateOpacityTable( OpacityTable )

    TYPE(OpacityTableType) :: OpacityTable

    CALL DeAllocateOpacityTypeA( OpacityTable % ECAPEM)
    CALL DeAllocateEquationOfStateTable( OpacityTable % EOSTable) 
    CALL DeAllocateEnergyGrid( OpacityTable % EnergyGrid ) 
 
  END SUBROUTINE DeAllocateOpacityTable


!==========================================================================
! subroutine for OpacityTable
!==========================================================================
  SUBROUTINE AllocateOpacityTypeA( ECAPEM, nSpeciesA, nPointsTS, nPointsE )  ! For OpacityTypeA (D ='4')

    TYPE(OpacityTypeA), ALLOCATABLE, DIMENSION(:), INTENT(inout) :: ECAPEM
    INTEGER, DIMENSION(3), INTENT(in) :: nPointsTS
    INTEGER, INTENT(in)               :: nPointsE
    INTEGER, INTENT(in)               :: nSpeciesA !species of interactions
 
    INTEGER :: i

    ALLOCATE( ECAPEM ( nSpeciesA ) )

    Do i = 1, nSpeciesA
       ALLOCATE( ECAPEM(i) % Values( nPointsE, nPointsTS(1), &
                  nPointsTS(2),nPointsTS(3))  ) ! E, rho, T, Ye
    END DO
 
  END SUBROUTINE AllocateOpacityTypeA

!---------------------------------------------

  SUBROUTINE DeAllocateOpacityTypeA( ECAPEM ) 

     TYPE(OpacityTypeA), DIMENSION(:), INTENT(inout) :: ECAPEM 
 
    INTEGER :: i
  
    DO i = 1, SIZE ( ECAPEM )
       DEALLOCATE( ECAPEM(i) % Values )
    END DO

  END SUBROUTINE DeAllocateOpacityTypeA

!==========================================================================
! subroutine for EnergyGrid
!==========================================================================
  SUBROUTINE AllocateEnergyGrid( EnergyGrid, nPointsE )

    TYPE(EnergyGridType), INTENT(inout)      :: EnergyGrid
    INTEGER, INTENT(in)                         :: nPointsE

    EnergyGrid % nPointsE = nPointsE 

    ALLOCATE( EnergyGrid % Values(1:nPointsE) ) 

  END SUBROUTINE AllocateEnergyGrid

!---------------------------------------------

  SUBROUTINE DeAllocateEnergyGrid( EnergyGrid )
    
    TYPE(EnergyGridType)            :: EnergyGrid
 
    DEALLOCATE( EnergyGrid % Values )  

  END SUBROUTINE DeAllocateEnergyGrid
 


!========================Function=============================
REAL FUNCTION totalECapEm(energy, rho, T, Z, A, chem_e, chem_n, chem_p, xheavy, xn, xp )

    USE wlKindModule, ONLY: dp
    USE wlExtPhysicalConstantsModule ! kmev, therm1, therm2, dmnp, me, rmu (mean baryon mass [g])
  
	REAL(dp) :: energy, rho, T, Z, A, chem_e, chem_n, chem_p, xheavy, xn, xp   
        REAL(dp) :: TMeV, n, qpri, nhn, npz, etapn, jnucleon, jnuclear
        TMeV = T * kmev                      ! kmev = 8.61733d-11 [MeV K^{-1}]
        n = a - z
        qpri = chem_n - chem_p + 3.0         ! [MeV] 3MeV: energy of the 1f5/2 level
        
 
        if(n.le.34.0) nhn= 6.0
        if(n.gt.34.0.and.n.le.40.0) nhn= 40.0 - n
        if(n.gt.40.0) nhn= 0.0

        if(z.le.20.0) npz= 0.0
        if(z.gt.20.0.and.z.le.28.0) npz= z - 20.0
        if(z.gt.28.0) npz= 8.0

        etapn = rho * ( xn - xp ) / ( rmu * ( EXP( (chem_n-chem_p)/TMeV ) - 1 ) ) 

        jnucleon = therm1 * etapn * (energy+dmnp)**2 &
                   * SQRT( 1.0 - ( me / (energy+dmnp) )**2 ) &
                   / ( EXP( (energy+dmnp-chem_e) / TMeV ) + 1 )
   
        jnuclear = therm2 * rho * xheavy * npz * nhn * (energy+qpri)**2&
                   * SQRT( 1.0 - ( me / (energy+qpri) )**2 ) &
                   / ( rmu * a * ( EXP( (energy+qpri-chem_e) / TMeV ) + 1 ))

        totalECapEm = jnucleon + jnuclear

    RETURN 
END FUNCTION

END MODULE wlOpacityTableModule

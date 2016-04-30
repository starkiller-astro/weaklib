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
!       wlGridModule, ONLY: EnergyGridType 
!       wlEquationOfStateTableModule
!       wlIOModuleHDF
!       wlEOSIOModuleHDF, ONLY: ReadEquationOfStateTableHDF
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
  USE wlEnergyGridModule, ONLY: &
    EnergyGridType, &
    AllocateEnergyGrid, &
    DeallocateEnergyGrid, &
    DescribeEnergyGrid
  USE wlOpacityFieldsModule, ONLY: &
    OpacityTypeA, &
    OpacityTypeB, &
    OpacityTypeC, &
    AllocateOpacity, &
    DeallocateOpacity
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType, &
    DeAllocateEquationOfStateTable

  IMPLICIT NONE
  PRIVATE

!=======================================
! OpacityTableType
!=======================================

  TYPE, PUBLIC :: OpacityTableType
    INTEGER                        :: nOpacitiesA
    INTEGER                        :: nOpacitiesB
    INTEGER                        :: nMomentsB
    INTEGER                        :: nOpacitiesC
    INTEGER                        :: nMomentsC
    INTEGER                        :: nPointsE
    INTEGER, DIMENSION(3)          :: nPointsTS
    TYPE(EnergyGridType)           :: EnergyGrid
    TYPE(EquationOfStateTableType) :: EOSTable
    TYPE(OpacityTypeA)             :: &
      ecap       ! -- Electron Capture
    TYPE(OpacityTypeB)             :: &
      scatt_Iso  ! -- Isoenergenic Scattering
    TYPE(OpacityTypeC)             :: &
      scatt_nIso ! -- Non-Isoenergenic Scattering
  END TYPE OpacityTableType

  PUBLIC :: AllocateOpacityTable
  PUBLIC :: DeAllocateOpacityTable
  PUBLIC :: DescribeOpacityTable

CONTAINS


!==========================================================================
! Public Subroutine for OpacityTable
!==========================================================================


  SUBROUTINE AllocateOpacityTable &
               ( OpTab, nOpacA, nOpacB, nMomB, nOpacC, nMomC, nPointsE )

    TYPE(OpacityTableType), INTENT(inout) :: OpTab
    INTEGER, INTENT(in)                   :: nOpacA
    INTEGER, INTENT(in)                   :: nOpacB, nMomB
    INTEGER, INTENT(in)                   :: nOpacC, nMomC
    INTEGER, INTENT(in)                   :: nPointsE

    CALL AllocateEnergyGrid( OpTab % EnergyGrid, nPointsE )

    CALL InitializeHDF( )
    CALL ReadEquationOfStateTableHDF &
           ( OpTab % EOSTable, "EquationOfStateTable.h5" )
    CALL FinalizeHDF( )

    OpTab % nPointsE  = nPointsE
    OpTab % nPointsTS = OpTab % EOSTable % TS % nPoints

    ASSOCIATE( nPoints => OpTab % EOSTable % nPoints )

    CALL AllocateOpacity &
           ( OpTab % ecap,       [ nPointsE, nPoints ], &
             nOpacities = nOpacA )

    CALL AllocateOpacity &
           ( OpTab % scatt_Iso,  [ nPointsE, nPoints ], &
             nMoments = nMomB, nOpacities = nOpacB )

    CALL AllocateOpacity &
           ( OpTab % scatt_nIso, [ nPointsE, nPoints ], &
             nMoments = nMomC, nOpacities = nOpacC )

    END ASSOCIATE ! nPoints

  END SUBROUTINE AllocateOpacityTable


  SUBROUTINE DeAllocateOpacityTable( OpTab )

    TYPE(OpacityTableType) :: OpTab

    CALL DeAllocateOpacity( OpTab % ecap )
    CALL DeAllocateEquationOfStateTable( OpTab % EOSTable )
    CALL DeAllocateEnergyGrid( OpTab % EnergyGrid ) 
 
  END SUBROUTINE DeAllocateOpacityTable


  SUBROUTINE DescribeOpacityTable( OpTab )

    TYPE(OpacityTableType) :: OpTab

    WRITE(*,*)
    WRITE(*,'(A2,A)') ' ', 'DescribeOpacityTable'

    CALL DescribeEnergyGrid( OpTab % EnergyGrid )

  END SUBROUTINE DescribeOpacityTable


!!$  SUBROUTINE AllocateEmptyOpacityTable &
!!$               ( OpacityTable, nSpeciesA, nPointsE, nPoints, nVariables )
!!$
!!$    TYPE(OpacityTableType), INTENT(inout)  :: OpacityTable
!!$    INTEGER, DIMENSION(3), INTENT(in)      :: nPoints
!!$    INTEGER, INTENT(in)                    :: nPointsE
!!$    INTEGER, INTENT(in)                    :: nSpeciesA
!!$    INTEGER, INTENT(in)                    :: nVariables
!!$
!!$    WRITE(*,*)'Allocating EOSTable'
!!$    CALL AllocateEquationOfStateTable( OpacityTable % EOSTable, nPoints, nVariables )    
!!$   
!!$    WRITE(*,*) 'Allocating OpacityTypeA'
!!$    CALL AllocateOpacityTypeA &
!!$           ( OpacityTable % ECAPEM, nSpeciesA, &
!!$             OpacityTable % EOSTable % nPoints, nPointsE )  ! just for TypeA
!!$
!!$    WRITE(*,*) 'Allocating EnergyGrid'
!!$    CALL AllocateEnergyGrid( OpacityTable % EnergyGrid, nPointsE )
!!$
!!$    WRITE(*,*) 'Allocation Complete'
!!$
!!$  END SUBROUTINE AllocateEmptyOpacityTable



!==========================================================================
! subroutine for OpacityTable
!==========================================================================

!!$  SUBROUTINE AllocateOpacityTypeA( ECAPEM, nSpeciesA, nPointsTS, nPointsE )  ! For OpacityTypeA (D ='4')
!!$
!!$    TYPE(OpacityTypeA), ALLOCATABLE, DIMENSION(:), INTENT(inout) :: ECAPEM
!!$    INTEGER, DIMENSION(3), INTENT(in) :: nPointsTS
!!$    INTEGER, INTENT(in)               :: nPointsE
!!$    INTEGER, INTENT(in)               :: nSpeciesA !species of interactions
!!$ 
!!$    INTEGER :: i
!!$
!!$    ALLOCATE( ECAPEM ( nSpeciesA ) )
!!$
!!$    Do i = 1, nSpeciesA
!!$       ALLOCATE( ECAPEM(i) % Values( nPointsE, nPointsTS(1), &
!!$                  nPointsTS(2),nPointsTS(3))  ) ! E, rho, T, Ye
!!$    END DO
!!$ 
!!$  END SUBROUTINE AllocateOpacityTypeA

!---------------------------------------------

!!$  SUBROUTINE DeAllocateOpacityTypeA( ECAPEM ) 
!!$
!!$    TYPE(OpacityTypeA), DIMENSION(:), INTENT(inout) :: ECAPEM 
!!$ 
!!$    INTEGER :: i
!!$  
!!$    DO i = 1, SIZE ( ECAPEM )
!!$      DEALLOCATE( ECAPEM(i) % Values )
!!$    END DO
!!$
!!$  END SUBROUTINE DeAllocateOpacityTypeA

END MODULE wlOpacityTableModule

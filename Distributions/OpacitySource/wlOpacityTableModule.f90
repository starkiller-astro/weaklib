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
    DeallocateOpacity, &
    DescribeOpacity
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType, &
    DeAllocateEquationOfStateTable, &
    AllocateEquationOfStateTable

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
      thermEmAb  ! -- Thermal Emission and Absorption
    TYPE(OpacityTypeB)             :: &
      scatt_Iso  ! -- Isoenergenic Scattering
    TYPE(OpacityTypeC)             :: &
      scatt_nIso ! -- Non-Isoenergenic Scattering
  END TYPE OpacityTableType

  PUBLIC :: AllocateOpacityTable
  PUBLIC :: DeAllocateOpacityTable
  PUBLIC :: DescribeOpacityTable

CONTAINS


  SUBROUTINE AllocateOpacityTable &
               ( OpTab, nOpacA, nOpacB, nMomB, nOpacC, nMomC, nPointsE )

    TYPE(OpacityTableType), INTENT(inout) :: OpTab
    INTEGER, INTENT(in)                   :: nOpacA
    INTEGER, INTENT(in)                   :: nOpacB, nMomB
    INTEGER, INTENT(in)                   :: nOpacC, nMomC
    INTEGER, INTENT(in)                   :: nPointsE

    CALL AllocateEnergyGrid( OpTab % EnergyGrid, nPointsE )

!    CALL InitializeHDF( )
    WRITE(*,*)
    WRITE(*,'(A2,A)') ' ', 'Reading EquationOfStateTable.h5 '

    CALL ReadEquationOfStateTableHDF &
           ( OpTab % EOSTable, "EquationOfStateTable.h5" )
!    CALL FinalizeHDF( )

    OpTab % nOpacitiesA = nOpacA
    OpTab % nOpacitiesB = nOpacB
    OpTab % nMomentsB = nMomB
    OpTab % nOpacitiesC = nOpacC
    OpTab % nMomentsC = nMomC
    OpTab % nPointsE  = nPointsE
    OpTab % nPointsTS = OpTab % EOSTable % TS % nPoints

    ASSOCIATE( nPoints => OpTab % EOSTable % nPoints )

    CALL AllocateOpacity &
           ( OpTab % thermEmAb,  [ nPointsE, nPoints ], &
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

    CALL DeAllocateOpacity( OpTab % thermEmAb )
    CALL DeAllocateEquationOfStateTable( OpTab % EOSTable )
    CALL DeAllocateEnergyGrid( OpTab % EnergyGrid ) 
 
  END SUBROUTINE DeAllocateOpacityTable


  SUBROUTINE DescribeOpacityTable( OpTab )

    TYPE(OpacityTableType) :: OpTab

    WRITE(*,*)
    WRITE(*,'(A2,A)') ' ', 'DescribeOpacityTable'

    CALL DescribeEnergyGrid( OpTab % EnergyGrid )
    CALL DescribeOpacity( OpTab % thermEmAb )
    CALL DescribeOpacity( OpTab % scatt_Iso )
    CALL DescribeOpacity( OpTab % scatt_nIso )

  END SUBROUTINE DescribeOpacityTable

END MODULE wlOpacityTableModule

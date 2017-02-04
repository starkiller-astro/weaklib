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
 
  USE wlKindModule, ONLY: &
    dp
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DeallocateGrid, &
    DescribeGrid
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
    INTEGER                        :: nOpacitiesB, nMomentsB
    INTEGER                        :: nOpacitiesB_NES, nMomentsB_NES
    INTEGER                        :: nOpacitiesC, nMomentsC
    INTEGER                        :: nPointsE, nPointsEta
    INTEGER, DIMENSION(3)          :: nPointsTS
    TYPE(GridType)                 :: EnergyGrid
    TYPE(GridType)                 :: EtaGrid     ! -- eletron chemical potential / kT
    TYPE(EquationOfStateTableType) :: EOSTable
    TYPE(OpacityTypeA)             :: &
      thermEmAb  ! -- Thermal Emission and Absorption
    TYPE(OpacityTypeB)             :: &
      scatt_Iso  ! -- Isoenergenic Scattering
    TYPE(OpacityTypeB)             :: &
      scatt_NES  ! -- Inelastic Neutrino-Electron Scattering
    TYPE(OpacityTypeC)             :: &
      scatt_nIso ! -- Non-Isoenergenic Scattering
  END TYPE OpacityTableType

  PUBLIC :: AllocateOpacityTable
  PUBLIC :: DeAllocateOpacityTable
  PUBLIC :: DescribeOpacityTable

CONTAINS

  SUBROUTINE AllocateOpacityTable &
               ( OpTab, nOpacA, nOpacB, nMomB, nOpacB_NES, nMomB_NES, nOpacC, nMomC, nPointsE, nPointsEta )

    TYPE(OpacityTableType), INTENT(inout) :: OpTab
    INTEGER, INTENT(in)                   :: nOpacA
    INTEGER, INTENT(in)                   :: nOpacB, nMomB
    INTEGER, INTENT(in)                   :: nOpacB_NES, nMomB_NES
    INTEGER, INTENT(in)                   :: nOpacC, nMomC
    INTEGER, INTENT(in)                   :: nPointsE
    INTEGER, INTENT(in)                   :: nPointsEta
    INTEGER, DIMENSION(4)                 :: nPointsTemp

    WRITE(*,*)
    WRITE(*,'(A2,A)') ' ', 'Reading wl-EOS-LS220-20-40-100-Lower-T.h5 ... '

    CALL ReadEquationOfStateTableHDF &
           ( OpTab % EOSTable, "wl-EOS-LS220-20-40-100-Lower-T.h5" )

    WRITE(*,*) 'Read EOS sucessfully.'
    WRITE(*,*) 'Pass the parameter and allocate OpacityTable ... '

    OpTab % nOpacitiesA     = nOpacA
    OpTab % nOpacitiesB     = nOpacB
    OpTab % nMomentsB       = nMomB
    OpTab % nOpacitiesB_NES = nOpacB_NES
    OpTab % nMomentsB_NES   = nMomB_NES
    OpTab % nOpacitiesC     = nOpacC
    OpTab % nMomentsC       = nMomC
    OpTab % nPointsE        = nPointsE
    OpTab % nPointsEta      = nPointsEta
    OpTab % nPointsTS       = OpTab % EOSTable % TS % nPoints

    CALL AllocateGrid( OpTab % EnergyGrid, nPointsE   )
    CALL AllocateGrid( OpTab % EtaGrid,    nPointsEta )

    ASSOCIATE( nPoints => OpTab % EOSTable % nPoints )

    nPointsTemp(1:4) = [ nPointsE, nPoints ]

    CALL AllocateOpacity &
           ( OpTab % thermEmAb, nPointsTemp, &
             nOpacities = nOpacA )

    CALL AllocateOpacity &
           ( OpTab % scatt_Iso, nPointsTemp, &
             nMoments = nMomB, nOpacities = nOpacB )

    CALL AllocateOpacity &
           ( OpTab % scatt_nIso, nPointsTemp, &
             nMoments = nMomC, nOpacities = nOpacC )
  
    nPointsTemp(1:4) = [ nPointsE, nPointsE, nPoints(2), nPointsEta]

    CALL AllocateOpacity &
           ( OpTab % scatt_NES, nPointsTemp, &
             nMoments = nMomB_NES, nOpacities = nOpacB_NES )

    END ASSOCIATE ! nPoints

  END SUBROUTINE AllocateOpacityTable


  SUBROUTINE DeAllocateOpacityTable( OpTab )

    TYPE(OpacityTableType) :: OpTab
    
    WRITE(*,*)
    WRITE(*,*) 'DeAllocate OpacityTable'

    CALL DeAllocateOpacity( OpTab % thermEmAb )
    CALL DeAllocateOpacity( OpTab % scatt_Iso )
    CALL DeAllocateOpacity( OpTab % scatt_NES )
    CALL DeAllocateOpacity( OpTab % scatt_nIso )

    CALL DeAllocateGrid( OpTab % EnergyGrid ) 
    CALL DeAllocateGrid( OpTab % EtaGrid ) 

    CALL DeAllocateEquationOfStateTable( OpTab % EOSTable )
 
  END SUBROUTINE DeAllocateOpacityTable


  SUBROUTINE DescribeOpacityTable( OpTab )

    TYPE(OpacityTableType) :: OpTab

    WRITE(*,*)
    WRITE(*,'(A2,A)') ' ', 'DescribeOpacityTable'

    CALL DescribeGrid( OpTab % EnergyGrid )
    CALL DescribeGrid( OpTab % EtaGrid )

    CALL DescribeOpacity( OpTab % thermEmAb )
    CALL DescribeOpacity( OpTab % scatt_Iso )
    CALL DescribeOpacity( OpTab % scatt_NES )
    CALL DescribeOpacity( OpTab % scatt_nIso )

  END SUBROUTINE DescribeOpacityTable

END MODULE wlOpacityTableModule

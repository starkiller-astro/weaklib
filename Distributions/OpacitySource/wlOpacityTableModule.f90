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
  USE wlThermoStateModule, ONLY: &
    ThermoStateType, &
    AllocateThermoState, &
    DeAllocateThermoState, &
    CopyThermoState  

  IMPLICIT NONE
  PRIVATE

!=======================================
! OpacityTableType
!=======================================

  TYPE, PUBLIC :: OpacityTableType
    INTEGER                        :: nOpacitiesA
    INTEGER                        :: nOpacitiesB, nMomentsB
    INTEGER                        :: nOpacitiesB_NES, nMomentsB_NES
    INTEGER                        :: nOpacitiesB_TP, nMomentsB_TP
    INTEGER                        :: nOpacitiesC, nMomentsC
    INTEGER                        :: nPointsE, nPointsEta
    INTEGER, DIMENSION(3)          :: nPointsTS
    TYPE(GridType)                 :: EnergyGrid
    TYPE(GridType)                 :: EtaGrid     ! -- eletron chemical potential / kT
    TYPE(EquationOfStateTableType) :: EOSTable
    TYPE(ThermoStateType)          :: TS
    TYPE(OpacityTypeA)             :: &
      thermEmAb  ! -- Thermal Emission and Absorption
    TYPE(OpacityTypeB)             :: &
      scatt_Iso  ! -- Isoenergenic Scattering
    TYPE(OpacityTypeB)             :: &
      scatt_NES  ! -- Inelastic Neutrino-Electron Scattering
    TYPE(OpacityTypeB)             :: &
      scatt_TP   ! -- Thermal Production
    TYPE(OpacityTypeC)             :: &
      scatt_nIso ! -- Non-Isoenergenic Scattering
  END TYPE OpacityTableType

  PUBLIC :: AllocateOpacityTable
  PUBLIC :: DeAllocateOpacityTable
  PUBLIC :: DescribeOpacityTable

CONTAINS

  SUBROUTINE AllocateOpacityTable &
               ( OpTab, nOpacA, nOpacB, nMomB, nOpacB_NES, nMomB_NES, &
                 nOpacB_TP, nMomB_TP, nOpacC, nMomC, nPointsE, nPointsEta )

    TYPE(OpacityTableType), INTENT(inout) :: OpTab
    INTEGER, INTENT(in)                   :: nOpacA
    INTEGER, INTENT(in)                   :: nOpacB, nMomB
    INTEGER, INTENT(in)                   :: nOpacB_NES, nMomB_NES
    INTEGER, INTENT(in)                   :: nOpacB_TP, nMomB_TP
    INTEGER, INTENT(in)                   :: nOpacC, nMomC
    INTEGER, INTENT(in)                   :: nPointsE
    INTEGER, INTENT(in)                   :: nPointsEta
    INTEGER, DIMENSION(4)                 :: nPointsTemp

    WRITE(*,*)
    WRITE(*,*) '  1st Allocate OpacityTable ...          ' 
    WRITE(*,*) '      Reading EquationOfStateTable  ...  '
    WRITE(*,*) '      File: EquationOfStateTable.h5  ... '

    CALL ReadEquationOfStateTableHDF &
           ( OpTab % EOSTable, "EquationOfStateTable.h5" )

    WRITE(*,*) 'Read EOS sucessfully.'
    WRITE(*,*)
    WRITE(*,*) 'Pass the parameter and allocate OpacityTable ... '

    OpTab % nOpacitiesA     = nOpacA
    OpTab % nOpacitiesB     = nOpacB
    OpTab % nMomentsB       = nMomB
    OpTab % nOpacitiesB_NES = nOpacB_NES
    OpTab % nMomentsB_NES   = nMomB_NES
    OpTab % nOpacitiesB_TP  = nOpacB_TP
    OpTab % nMomentsB_TP    = nMomB_TP
    OpTab % nOpacitiesC     = nOpacC
    OpTab % nMomentsC       = nMomC
    OpTab % nPointsE        = nPointsE
    OpTab % nPointsEta      = nPointsEta
    OpTab % nPointsTS       = OpTab % EOSTable % TS % nPoints

    CALL AllocateGrid( OpTab % EnergyGrid, nPointsE   )
    CALL AllocateGrid( OpTab % EtaGrid,    nPointsEta )
    CALL AllocateThermoState( OpTab % TS, OpTab % EOSTable % TS % nPoints )

    CALL CopyThermoState( OpTab % TS, OpTab % EOSTable % TS )

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

    CALL AllocateOpacity &
           ( OpTab % scatt_TP, nPointsTemp, &
             nMoments = nMomB_TP, nOpacities = nOpacB_TP )

    END ASSOCIATE ! nPoints

    PRINT*, 'End allocte opacity table.'
    PRINT*

  END SUBROUTINE AllocateOpacityTable


  SUBROUTINE DeAllocateOpacityTable( OpTab )

    TYPE(OpacityTableType) :: OpTab
    
    WRITE(*,*)
    WRITE(*,*) 'DeAllocate OpacityTable'

    CALL DeAllocateOpacity( OpTab % thermEmAb )
    CALL DeAllocateOpacity( OpTab % scatt_Iso )
    CALL DeAllocateOpacity( OpTab % scatt_NES )
    CALL DeAllocateOpacity( OpTab % scatt_TP )
    CALL DeAllocateOpacity( OpTab % scatt_nIso )

    CALL DeAllocateThermoState( OpTab % TS )

    CALL DeAllocateGrid( OpTab % EnergyGrid ) 
    CALL DeAllocateGrid( OpTab % EtaGrid ) 

    CALL DeAllocateEquationOfStateTable( OpTab % EOSTable )
 
  END SUBROUTINE DeAllocateOpacityTable


  SUBROUTINE DescribeOpacityTable( OpTab )

    TYPE(OpacityTableType) :: OpTab
    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A2,A)') ' ', 'DescribeOpacityTable'

    ASSOCIATE( TS => OpTab % TS )

    WRITE(*,*)
    WRITE(*,'(A2,A)') ' ', 'DescribeThermoState'
    DO i = 1, 3

      WRITE(*,*)
      WRITE(*,'(A5,A21,I1.1,A3,A32)') &
        '', 'Independent Variable ', i, ' = ', TRIM( TS % Names(i) )
      WRITE(*,'(A7,A7,A20)') &
        '', 'Units: ', TRIM( TS % Units(i) )
      WRITE(*,'(A7,A12,ES12.6E2)') &
        '', 'Min Value = ', TS % minValues(i)
      WRITE(*,'(A7,A12,ES12.6E2)') &
        '', 'Max Value = ', TS % maxValues(i)
      WRITE(*,'(A7,A12,I4.4)') &
        '', 'nPoints   = ', TS % nPoints(i)

      IF ( TS % LogInterp(i) == 1 ) THEN
        WRITE (*,'(A7,A27)') &
          '', 'Grid Logarithmically Spaced'
      ELSE
        WRITE (*,'(A7,A20)') &
          '', 'Grid Linearly Spaced'
      END IF

    END DO

    END ASSOCIATE ! TS
    

    CALL DescribeGrid( OpTab % EnergyGrid )
    CALL DescribeGrid( OpTab % EtaGrid )

    CALL DescribeOpacity( OpTab % thermEmAb )
    CALL DescribeOpacity( OpTab % scatt_Iso )
    CALL DescribeOpacity( OpTab % scatt_NES )
    CALL DescribeOpacity( OpTab % scatt_TP )
    CALL DescribeOpacity( OpTab % scatt_nIso )

  END SUBROUTINE DescribeOpacityTable

END MODULE wlOpacityTableModule

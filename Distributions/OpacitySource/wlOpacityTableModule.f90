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
!-----------------------------------------------------------------------
!                         Three Opacity Type
!
! OpacityType EmAb for  ABEM( rho, T, Ye, E )
!
!                e- + p/A <--> v_e + n/A*
!
! OpacityType Scat for  ISO( e, l, rho, T, Ye )
!                        
!                v_i/anti(v_i) + A --> v_i/anti(v_i) + A 
!                v_i/anti(v_i) + e+/e-/n/p  <-->  v_i/anti(v_i) + e+/e-/n/p
!
!                  and  NES ( ep, e, l, T, eta )
!                       Pair( ep, e, l, T, eta )  
!-----------------------------------------------------------------------
 
  USE wlKindModule, ONLY: &
    dp
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DeallocateGrid, &
    DescribeGrid
  USE wlOpacityFieldsModule, ONLY: &
    OpacityTypeEmAb, &
    OpacityTypeScat, &
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
    INTEGER        :: nOpacities_EmAb
    INTEGER        :: nOpacities_Iso,  nMoments_Iso
    INTEGER        :: nOpacities_NES,  nMoments_NES
    INTEGER        :: nOpacities_Pair, nMoments_Pair
    INTEGER        :: nPointsE
    INTEGER        :: nPointsEta
    INTEGER        :: nPointsTS(3)
    TYPE(GridType) :: EnergyGrid
    TYPE(GridType) :: EtaGrid ! -- eletron chemical potential / kT
    TYPE(EquationOfStateTableType) :: EOSTable
    TYPE(ThermoStateType)          :: TS
    TYPE(OpacityTypeEmAb)          :: &
      EmAb       ! -- Corrected Absorption Opacity
    TYPE(OpacityTypeScat)          :: &
      Scat_Iso   ! -- Isoenergenic Scattering
    TYPE(OpacityTypeScat)          :: &
      Scat_NES   ! -- Inelastic Neutrino-Electron Scattering
    TYPE(OpacityTypeScat)          :: &
      Scat_Pair  ! -- Pair Production
  END TYPE OpacityTableType

  PUBLIC :: AllocateOpacityTable
  PUBLIC :: DeAllocateOpacityTable
  PUBLIC :: DescribeOpacityTable

CONTAINS

  SUBROUTINE AllocateOpacityTable &
    ( OpTab, nOpac_EmAb, nOpac_Iso, nMom_Iso, nOpac_NES, nMom_NES, &
      nOpac_Pair, nMom_Pair, nPointsE, nPointsEta, &
      EquationOfStateTableName_Option, Verbose_Option )

    TYPE(OpacityTableType), INTENT(inout)        :: OpTab
    INTEGER,                INTENT(in)           :: nOpac_EmAb
    INTEGER,                INTENT(in)           :: nOpac_Iso, nMom_Iso
    INTEGER,                INTENT(in)           :: nOpac_NES, nMom_NES
    INTEGER,                INTENT(in)           :: nOpac_Pair, nMom_Pair
    INTEGER,                INTENT(in)           :: nPointsE
    INTEGER,                INTENT(in)           :: nPointsEta
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    LOGICAL,                INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL        :: Verbose
    CHARACTER(256) :: EquationOfStateTableName
    INTEGER        :: nPointsTemp(5)

    IF( PRESENT( EquationOfStateTableName_Option ) )THEN
       EquationOfStateTableName = TRIM( EquationOfStateTableName_Option )
    ELSE
       EquationOfStateTableName = 'EquationOfStateTable.h5'
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'AllocateOpacityTable'
      WRITE(*,*)
      WRITE(*,'(A6,A9,A)') &
        '', 'Reading: ', TRIM( EquationOfStateTableName )
    END IF

    CALL ReadEquationOfStateTableHDF &
           ( OpTab % EOSTable, TRIM( EquationOfStateTableName ) )

    OpTab % nOpacities_EmAb = nOpac_EmAb
    OpTab % nOpacities_Iso  = nOpac_Iso
    OpTab % nMoments_Iso    = nMom_Iso
    OpTab % nOpacities_NES  = nOpac_NES
    OpTab % nMoments_NES    = nMom_NES
    OpTab % nOpacities_Pair = nOpac_Pair
    OpTab % nMoments_Pair   = nMom_Pair
    OpTab % nPointsE        = nPointsE
    OpTab % nPointsEta      = nPointsEta
    OpTab % nPointsTS       = OpTab % EOSTable % TS % nPoints

    CALL AllocateGrid( OpTab % EnergyGrid, nPointsE   )
    CALL AllocateGrid( OpTab % EtaGrid,    nPointsEta )
    CALL AllocateThermoState( OpTab % TS, OpTab % EOSTable % TS % nPoints )

    CALL CopyThermoState( OpTab % TS, OpTab % EOSTable % TS )

    ASSOCIATE( nPoints => OpTab % EOSTable % nPoints, &
               iT      => OpTab % EOSTable % TS % Indices % iT )

    nPointsTemp(1:4) = [ nPointsE, nPoints ]

    CALL AllocateOpacity &
           ( OpTab % EmAb, nPointsTemp(1:4), nOpacities = nOpac_EmAb )

    nPointsTemp(1:5) = [ nPointsE, nMom_Iso, nPoints ]
    CALL AllocateOpacity &
           ( OpTab % Scat_Iso, nPointsTemp(1:5), &
             nMoments = nMom_Iso, nOpacities = nOpac_Iso )
  
    nPointsTemp(1:5) = &
           [ nPointsE, nPointsE, nMom_NES, nPoints(iT), nPointsEta]

    CALL AllocateOpacity &
           ( OpTab % Scat_NES, nPointsTemp(1:5), &
             nMoments = nMom_NES, nOpacities = nOpac_NES )

    nPointsTemp(1:5) = &
           [ nPointsE, nPointsE, nMom_Pair, nPoints(iT), nPointsEta]

    CALL AllocateOpacity &
           ( OpTab % Scat_Pair, nPointsTemp(1:5), &
             nMoments = nMom_Pair, nOpacities = nOpac_Pair )

    END ASSOCIATE ! nPoints

  END SUBROUTINE AllocateOpacityTable


  SUBROUTINE DeAllocateOpacityTable( OpTab, Verbose_Option )

    TYPE(OpacityTableType) :: OpTab
    LOGICAL, OPTIONAL      :: Verbose_Option

    LOGICAL :: Verbose
    
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'Deallocating Opacity Table'
    END IF

    CALL DeAllocateOpacity( OpTab % EmAb ) 
    CALL DeAllocateOpacity( OpTab % Scat_Iso )
    CALL DeAllocateOpacity( OpTab % Scat_NES )
    CALL DeAllocateOpacity( OpTab % Scat_Pair )

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

    CALL DescribeOpacity( OpTab % EmAb ) 
    CALL DescribeOpacity( OpTab % Scat_Iso )
    CALL DescribeOpacity( OpTab % Scat_NES )
    CALL DescribeOpacity( OpTab % Scat_Pair )

  END SUBROUTINE DescribeOpacityTable

END MODULE wlOpacityTableModule

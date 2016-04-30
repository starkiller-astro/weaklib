MODULE wlOpacityFieldsModule

  USE wlKindModule, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PARAMETER :: iNu_e     = 1
  INTEGER, PUBLIC, PARAMETER :: iNu_e_bar = 2
  INTEGER, PUBLIC, PARAMETER :: iNu_x     = 3
  INTEGER, PUBLIC, PARAMETER :: iNu_x_bar = 4
  INTEGER, PUBLIC, PARAMETER :: nSpecies  = 4

!---------------------------------------------
!
! OpacityTypeA (Emission / Absorption)
!   Dependency ( E, rho, T, Ye )
!     E:   Neutrino Energy
!     rho: Mass Density
!     T:   Temperature
!     Ye:  Electron Fraction
!
!---------------------------------------------

  TYPE :: ValueTypeA
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Values
  END type ValueTypeA

  TYPE, PUBLIC :: OpacityTypeA
    INTEGER :: nOpacities
    INTEGER, DIMENSION(4) :: nPoints
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Names
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Species
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Units
    TYPE(ValueTypeA),  DIMENSION(:), ALLOCATABLE :: Absorptivity
  END TYPE

!---------------------------------------------
!
! OpacityType B (Elastic Scattering)
!   Dependency ( E, rho, T, Ye, l )
!     E:   Neutrino Energy
!     rho: Mass Density
!     T:   Temperature
!     Ye:  Electron Fraction
!     l:   Legendre Moment
!
!---------------------------------------------  

  TYPE :: ValueTypeB
    REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Values
  END type ValueTypeB

  TYPE, PUBLIC :: OpacityTypeB
    INTEGER :: nOpacities
    INTEGER :: nMoments
    INTEGER, DIMENSION(4) :: nPoints
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Names
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Species
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Units
    TYPE(ValueTypeB),  DIMENSION(:), ALLOCATABLE :: Kernel
  END TYPE

!---------------------------------------------
!
! OpacityType C (Inelastic Scattering)
!   Dependency ( E_in, E_out, rho, T, Ye, l )
!
!---------------------------------------------

  TYPE :: ValueTypeC
    REAL(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Values
  END type ValueTypeC

  TYPE, PUBLIC :: OpacityTypeC
    INTEGER :: nOpacities
    INTEGER :: nMoments
    INTEGER, DIMENSION(4) :: nPoints
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Names
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Species
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Units
    TYPE(ValueTypeC),  DIMENSION(:), ALLOCATABLE :: Kernel
  END TYPE

  PUBLIC :: AllocateOpacity
  PUBLIC :: DeallocateOpacity

  INTERFACE AllocateOpacity
    MODULE PROCEDURE AllocateOpacityTypeA
    MODULE PROCEDURE AllocateOpacityTypeB
    MODULE PROCEDURE AllocateOpacityTypeC
  END INTERFACE AllocateOpacity

  INTERFACE DeallocateOpacity
    MODULE PROCEDURE DeallocateOpacityTypeA
    MODULE PROCEDURE DeallocateOpacityTypeB
    MODULE PROCEDURE DeallocateOpacityTypeC
  END INTERFACE DeallocateOpacity

CONTAINS


  SUBROUTINE AllocateOpacityTypeA( Opacity, nPoints, nOpacities )

    TYPE(OpacityTypeA), INTENT(inout) :: &
      Opacity
    INTEGER, DIMENSION(4), INTENT(in) :: &
      nPoints
    INTEGER, INTENT(in) :: &
      nOpacities

    INTEGER :: i

    Opacity % nOpacities = nOpacities
    Opacity % nPoints    = nPoints

    ALLOCATE( Opacity % Names(nOpacities) )
    ALLOCATE( Opacity % Species(nOpacities) )
    ALLOCATE( Opacity % Units(nOpacities) )
    ALLOCATE( Opacity % Absorptivity(nOpacities) )

    DO i = 1, nOpacities
      ALLOCATE( Opacity % Absorptivity(i) % Values &
                  (nPoints(1), nPoints(2), nPoints(3), nPoints(4)) )
    END DO

  END SUBROUTINE AllocateOpacityTypeA


  SUBROUTINE DeallocateOpacityTypeA( Opacity )

    TYPE(OpacityTypeA), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % Absorptivity(i) % Values )
    END DO

    DEALLOCATE( Opacity % Absorptivity )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Species )
    DEALLOCATE( Opacity % Names )

  END SUBROUTINE DeallocateOpacityTypeA


  SUBROUTINE AllocateOpacityTypeB( Opacity, nPoints, nMoments, nOpacities )

    TYPE(OpacityTypeB), INTENT(inout) :: &
      Opacity
    INTEGER, DIMENSION(4), INTENT(in) :: &
      nPoints
    INTEGER, INTENT(in) :: &
      nMoments, &
      nOpacities

    INTEGER :: i

    Opacity % nOpacities = nOpacities
    Opacity % nMoments   = nMoments
    Opacity % nPoints    = nPoints

    ALLOCATE( Opacity % Names(nOpacities) )
    ALLOCATE( Opacity % Species(nOpacities) )
    ALLOCATE( Opacity % Units(nOpacities) )
    ALLOCATE( Opacity % Kernel(nOpacities) )

    DO i = 1, nOpacities
      ALLOCATE( Opacity % Kernel(i) % Values &
                  (nPoints(1), nPoints(2), nPoints(3), nPoints(4), &
                   nMoments) )
    END DO

  END SUBROUTINE AllocateOpacityTypeB


  SUBROUTINE DeallocateOpacityTypeB( Opacity )

    TYPE(OpacityTypeB), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % Kernel(i) % Values )
    END DO

    DEALLOCATE( Opacity % Kernel )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Species )
    DEALLOCATE( Opacity % Names )

  END SUBROUTINE DeallocateOpacityTypeB


  SUBROUTINE AllocateOpacityTypeC( Opacity, nPoints, nMoments, nOpacities )

    TYPE(OpacityTypeC), INTENT(inout) :: &
      Opacity
    INTEGER, DIMENSION(4), INTENT(in) :: &
      nPoints
    INTEGER, INTENT(in) :: &
      nMoments, &
      nOpacities

    INTEGER :: i

    Opacity % nOpacities = nOpacities
    Opacity % nMoments   = nMoments
    Opacity % nPoints    = nPoints

    ALLOCATE( Opacity % Names(nOpacities) )
    ALLOCATE( Opacity % Species(nOpacities) )
    ALLOCATE( Opacity % Units(nOpacities) )
    ALLOCATE( Opacity % Kernel(nOpacities) )

    DO i = 1, nOpacities
      ALLOCATE( Opacity % Kernel(i) % Values &
                  (nPoints(1), nPoints(1), nPoints(2), nPoints(3), &
                   nPoints(4), nMoments) )
    END DO

  END SUBROUTINE AllocateOpacityTypeC


  SUBROUTINE DeallocateOpacityTypeC( Opacity )

    TYPE(OpacityTypeC), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % Kernel(i) % Values )
    END DO

    DEALLOCATE( Opacity % Kernel )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Species )
    DEALLOCATE( Opacity % Names )

  END SUBROUTINE DeallocateOpacityTypeC


END MODULE wlOpacityFieldsModule

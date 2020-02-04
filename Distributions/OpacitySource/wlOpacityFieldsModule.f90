MODULE wlOpacityFieldsModule

  USE wlKindModule, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PARAMETER :: iNu_e     = 1
  INTEGER, PUBLIC, PARAMETER :: iNu_e_bar = 2
  INTEGER, PUBLIC, PARAMETER :: iNu_x     = 3
  INTEGER, PUBLIC, PARAMETER :: iNu_x_bar = 4
  INTEGER, PUBLIC, PARAMETER :: nSpecies  = 4

  ! --- NES Scattering Kernels ---

  INTEGER, PUBLIC, PARAMETER :: iHi0  = 1
  INTEGER, PUBLIC, PARAMETER :: iHii0 = 2
  INTEGER, PUBLIC, PARAMETER :: iHi1  = 3
  INTEGER, PUBLIC, PARAMETER :: iHii1 = 4

  ! --- Pair Kernels ---

  INTEGER, PUBLIC, PARAMETER :: iJi0  = 1
  INTEGER, PUBLIC, PARAMETER :: iJii0 = 2
  INTEGER, PUBLIC, PARAMETER :: iJi1  = 3
  INTEGER, PUBLIC, PARAMETER :: iJii1 = 4

  TYPE :: ValueType_3D
    REAL(dp), ALLOCATABLE :: Values(:,:,:)
  END type ValueType_3D

  TYPE :: ValueType_4D
    REAL(dp), ALLOCATABLE :: Values(:,:,:,:)
  END type ValueType_4D

  TYPE :: ValueType_5D
    REAL(dp), ALLOCATABLE :: Values(:,:,:,:,:)
  END type ValueType_5D

!---------------------------------------------
!
! OpacityTypeEmAb (Emission / Absorption)
!   Dependency ( E, rho, T, Ye )
!     E:   Neutrino Energy
!     rho: Mass Density
!     T:   Temperature
!     Ye:  Electron Fraction
!
!---------------------------------------------

  TYPE, PUBLIC :: OpacityTypeEmAb
    INTEGER                         :: nOpacities
    INTEGER                         :: nPoints(4)
    CHARACTER(LEN=32),  ALLOCATABLE :: Names(:)
    CHARACTER(LEN=32),  ALLOCATABLE :: Units(:)
    REAL(dp),           ALLOCATABLE :: Offsets(:)
    TYPE(ValueType_4D), ALLOCATABLE :: Opacity(:)
  END TYPE OpacityTypeEmAb

!---------------------------------------------
!
! OpacityTypeScat (Elastic Scattering)
!   Dependency ( E, l, rho, T, Ye )
!     E:   Neutrino Energy
!     l:   Legendre Moment
!     rho: Mass Density
!     T:   Temperature
!     Ye:  Electron Fraction
!
! OpacityTypeScat (Inelastic Neutrino-Electron Scattering)
!   Dependency ( E', E, l, T, Eta )
!     E':  Neutrino Energy
!     E:   Neutrino Energy
!     l:   Legendre Moment
!     T:   Temperature
!     Eta: Electron Chemical Pot. / Temperature
!
!---------------------------------------------

  TYPE, PUBLIC :: OpacityTypeScat
    INTEGER                         :: nOpacities
    INTEGER                         :: nMoments
    INTEGER                         :: nPoints(5)
    CHARACTER(LEN=32),  ALLOCATABLE :: Names(:)
    CHARACTER(LEN=32),  ALLOCATABLE :: Units(:)
    REAL(dp),           ALLOCATABLE :: Offsets(:,:)
    TYPE(ValueType_5D), ALLOCATABLE :: Kernel(:)
  END TYPE OpacityTypeScat

  PUBLIC :: AllocateOpacity
  PUBLIC :: DeallocateOpacity
  PUBLIC :: DescribeOpacity

  INTERFACE AllocateOpacity
    MODULE PROCEDURE AllocateOpacityTypeEmAb
    MODULE PROCEDURE AllocateOpacityTypeScat
  END INTERFACE AllocateOpacity

  INTERFACE DeallocateOpacity
    MODULE PROCEDURE DeallocateOpacityTypeEmAb
    MODULE PROCEDURE DeallocateOpacityTypeScat
  END INTERFACE DeallocateOpacity

  INTERFACE DescribeOpacity
    MODULE PROCEDURE DescribeOpacityTypeEmAb
    MODULE PROCEDURE DescribeOpacityTypeScat
  END INTERFACE DescribeOpacity

CONTAINS


  SUBROUTINE AllocateOpacityTypeEmAb( Opacity, nPoints, nOpacities )

    TYPE(OpacityTypeEmAb), INTENT(inout) :: Opacity
    INTEGER,               INTENT(in)    :: nPoints(4)
    INTEGER,               INTENT(in)    :: nOpacities

    INTEGER :: i

    Opacity % nOpacities = nOpacities
    Opacity % nPoints    = nPoints
  
    ALLOCATE( Opacity % Names(nOpacities) )
    ALLOCATE( Opacity % Units(nOpacities) )
    ALLOCATE( Opacity % Offsets(nOpacities) )
    ALLOCATE( Opacity % Opacity(nOpacities) )

    DO i = 1, nOpacities

      ALLOCATE( Opacity % Opacity(i) % Values &
                  (nPoints(1), nPoints(2), nPoints(3), nPoints(4)) )
    END DO

  END SUBROUTINE AllocateOpacityTypeEmAb


  SUBROUTINE DeallocateOpacityTypeEmAb( Opacity )

    TYPE(OpacityTypeEmAb), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % Opacity(i) % Values )
    END DO

    DEALLOCATE( Opacity % Opacity )
    DEALLOCATE( Opacity % Offsets )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Names )

  END SUBROUTINE DeallocateOpacityTypeEmAb
  

  SUBROUTINE DescribeOpacityTypeEmAb( Opacity )

    TYPE(OpacityTypeEmAb), INTENT(in) :: Opacity

    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A4,A)') ' ', 'Opacity Type EmAb'
    WRITE(*,'(A4,A)') ' ', '--------------'
    WRITE(*,'(A6,A13,I3.3)') &
      ' ', 'nOpacities = ', Opacity % nOpacities
    WRITE(*,'(A6,A13,4I5.4)') &
      ' ', 'nPoints    = ', Opacity % nPoints
    WRITE(*,'(A6,A13,I10.10)') &
      ' ', 'DOFs       = ', &
      Opacity % nOpacities * PRODUCT( Opacity % nPoints )

    DO i = 1, Opacity % nOpacities
      WRITE(*,*)
      WRITE(*,'(A6,A8,I3.3,A3,A)') &
        ' ', 'Opacity(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Units     = ', TRIM( Opacity % Units(i) )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % Opacity(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % Opacity(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Offset    = ', Opacity % Offsets(i)
    END DO
    WRITE(*,*)

  END SUBROUTINE DescribeOpacityTypeEmAb


  SUBROUTINE AllocateOpacityTypeScat( Opacity, nPoints, nMoments, nOpacities )

    INTEGER,               INTENT(in)    :: nMoments
    INTEGER,               INTENT(in)    :: nOpacities
    INTEGER,               INTENT(in)    :: nPoints(5)
    TYPE(OpacityTypeScat), INTENT(inout) :: Opacity

    INTEGER :: i

    Opacity % nOpacities = nOpacities
    Opacity % nMoments   = nMoments
    Opacity % nPoints    = nPoints

    ALLOCATE( Opacity % Names(nOpacities) )
    ALLOCATE( Opacity % Units(nOpacities) )
    ALLOCATE( Opacity % Offsets(nOpacities, nMoments) )
    ALLOCATE( Opacity % Kernel(nOpacities) )

    DO i = 1, nOpacities

      ALLOCATE &
        ( Opacity % Kernel(i) % Values &
            ( nPoints(1), nPoints(2), nPoints(3), nPoints(4), nPoints(5) ) )

    END DO

  END SUBROUTINE AllocateOpacityTypeScat


  SUBROUTINE DeallocateOpacityTypeScat( Opacity )

    TYPE(OpacityTypeScat), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % Kernel(i) % Values )
    END DO

    DEALLOCATE( Opacity % Kernel )
    DEALLOCATE( Opacity % Offsets )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Names )

  END SUBROUTINE DeallocateOpacityTypeScat


  SUBROUTINE DescribeOpacityTypeScat( Opacity )

    TYPE(OpacityTypeScat), INTENT(in) :: Opacity

    INTEGER :: i, l

    WRITE(*,*)
    WRITE(*,'(A4,A)') ' ', 'Opacity Type Scat'
    WRITE(*,'(A4,A)') ' ', '--------------'
    WRITE(*,'(A6,A13,I3.3)') &
      ' ', 'nOpacities = ', Opacity % nOpacities
    WRITE(*,'(A6,A13,I3.3)') &
      ' ', 'nMoments   = ', Opacity % nMoments
    WRITE(*,'(A6,A13,5I5.4)') &
      ' ', 'nPoints    = ', Opacity % nPoints
    WRITE(*,'(A6,A13,I10.10)') &
      ' ', 'DOFs       = ', &
      Opacity % nOpacities * PRODUCT( Opacity % nPoints )

    DO i = 1, Opacity % nOpacities
      WRITE(*,*)
      WRITE(*,'(A6,A8,I3.3,A3,A)') &
        ' ', 'Opacity(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Units     = ', TRIM( Opacity % Units(i) )

      DO l = 1, Opacity % nMoments
      WRITE(*,*)
         IF( Opacity % nMoments .eq. 4 )THEN
           WRITE(*,'(A8,A16,I3.3)') &
             ' ', 'For Moments l = ', l
           WRITE(*,'(A8,A12,ES12.4E3)') &
             ' ', 'Min Value = ', MINVAL( Opacity % Kernel(i) % Values(:,:,l,:,:) )
           WRITE(*,'(A8,A12,ES12.4E3)') &
             ' ', 'Max Value = ', MAXVAL( Opacity % Kernel(i) % Values(:,:,l,:,:) )
         ELSE
           WRITE(*,'(A8,A16,I3.3)') &
             ' ', 'For Moments l = ', l
           WRITE(*,'(A8,A12,ES12.4E3)') &
             ' ', 'Min Value = ', MINVAL( Opacity % Kernel(i) % Values(:,l,:,:,:) )
           WRITE(*,'(A8,A12,ES12.4E3)') &
             ' ', 'Max Value = ', MAXVAL( Opacity % Kernel(i) % Values(:,l,:,:,:) )
         END IF
         WRITE(*,'(A8,A12,ES12.4E3)') &
           ' ', 'Offset    = ', Opacity % Offsets(i,l)
      END DO ! l = nMoment
      WRITE(*,*)
    END DO
    WRITE(*,*)

  END SUBROUTINE DescribeOpacityTypeScat


END MODULE wlOpacityFieldsModule

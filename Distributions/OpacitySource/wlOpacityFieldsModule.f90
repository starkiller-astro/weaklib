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

  TYPE :: ValueType_3D
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Values
  END type ValueType_3D

  TYPE :: ValueType_4D
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Values
  END type ValueType_4D

  TYPE, PUBLIC :: OpacityTypeA
    INTEGER :: nOpacities
    INTEGER, DIMENSION(4) :: nPoints
    REAL(dp),             DIMENSION(:), ALLOCATABLE :: Offsets
    CHARACTER(LEN=32),    DIMENSION(:), ALLOCATABLE :: Names
    CHARACTER(LEN=32),    DIMENSION(:), ALLOCATABLE :: Species
    CHARACTER(LEN=32),    DIMENSION(:), ALLOCATABLE :: Units
    TYPE(ValueType_3D),   DIMENSION(:), ALLOCATABLE :: GreyOpacity_Number_FD
    TYPE(ValueType_3D),   DIMENSION(:), ALLOCATABLE :: GreyOpacity_Energy_FD
    TYPE(ValueType_3D),   DIMENSION(:), ALLOCATABLE :: GreyMoment_Number_FD
    TYPE(ValueType_3D),   DIMENSION(:), ALLOCATABLE :: GreyMoment_Energy_FD
    TYPE(ValueType_4D),   DIMENSION(:), ALLOCATABLE :: Absorptivity
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
! OpacityType B (Inelastic Neutrino-Electron Scattering)
!   Dependency ( E', E, T, Eta, l )
!     E':  Neutrino Energy
!     E:   Neutrino Energy
!     T:   Temperature
!     Eta: Electron Chemical Pot. / Temperature
!     l:   Legendre Moment
!
!---------------------------------------------  

  TYPE :: ValueType_5D
    REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Values
  END type ValueType_5D

  TYPE, PUBLIC :: OpacityTypeB
    INTEGER :: nOpacities
    INTEGER :: nMoments
    INTEGER, DIMENSION(4) :: nPoints
    REAL(dp),            DIMENSION(:,:), ALLOCATABLE :: Offsets
    CHARACTER(LEN=32),   DIMENSION(:), ALLOCATABLE :: Names
    CHARACTER(LEN=32),   DIMENSION(:), ALLOCATABLE :: Species
    CHARACTER(LEN=32),   DIMENSION(:), ALLOCATABLE :: Units
    TYPE(ValueType_4D),  DIMENSION(:), ALLOCATABLE :: GreyOpacity_Number_FD
    TYPE(ValueType_4D),  DIMENSION(:), ALLOCATABLE :: GreyOpacity_Energy_FD
!    TYPE(ValueType_4D),  DIMENSION(:), ALLOCATABLE :: GreyMoment_Number_FD
!    TYPE(ValueType_4D),  DIMENSION(:), ALLOCATABLE :: GreyMoment_Energy_FD
    TYPE(ValueType_5D),  DIMENSION(:), ALLOCATABLE :: Kernel
  END TYPE

!---------------------------------------------
!
! OpacityType C (Inelastic Scattering)
!   Dependency ( E_out, E_in, rho, T, Ye, l )
!
!---------------------------------------------

  TYPE :: ValueType_6D
    REAL(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Values
  END type ValueType_6D

  TYPE, PUBLIC :: OpacityTypeC
    INTEGER :: nOpacities
    INTEGER :: nMoments
    INTEGER, DIMENSION(4) :: nPoints
    REAL(dp),          DIMENSION(:), ALLOCATABLE :: Offsets
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Names
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Species
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Units
    TYPE(ValueType_6D),  DIMENSION(:), ALLOCATABLE :: Kernel
  END TYPE

  PUBLIC :: AllocateOpacity
  PUBLIC :: DeallocateOpacity
  PUBLIC :: DescribeOpacity

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

  INTERFACE DescribeOpacity
    MODULE PROCEDURE DescribeOpacityTypeA
    MODULE PROCEDURE DescribeOpacityTypeB
    MODULE PROCEDURE DescribeOpacityTypeC
  END INTERFACE DescribeOpacity

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
    ALLOCATE( Opacity % Offsets(nOpacities) )
    ALLOCATE( Opacity % GreyOpacity_Number_FD(nOpacities) )
    ALLOCATE( Opacity % GreyOpacity_Energy_FD(nOpacities) )
    ALLOCATE( Opacity % GreyMoment_Number_FD(nOpacities) )
    ALLOCATE( Opacity % GreyMoment_Energy_FD(nOpacities) )
    ALLOCATE( Opacity % Absorptivity(nOpacities) )
    
    DO i = 1, nOpacities

      ALLOCATE( Opacity % GreyOpacity_Number_FD(i) % Values &
                  ( nPoints(2), nPoints(3), nPoints(4)) )

      ALLOCATE( Opacity % GreyOpacity_Energy_FD(i) % Values &
                  ( nPoints(2), nPoints(3), nPoints(4)) )

      ALLOCATE( Opacity % GreyMoment_Energy_FD(i) % Values &
                  ( nPoints(2), nPoints(3), nPoints(4)) )

      ALLOCATE( Opacity % GreyMoment_Number_FD(i) % Values &
                  ( nPoints(2), nPoints(3), nPoints(4)) )

      ALLOCATE( Opacity % Absorptivity(i) % Values &
                  (nPoints(1), nPoints(2), nPoints(3), nPoints(4)) )
    END DO

  END SUBROUTINE AllocateOpacityTypeA


  SUBROUTINE DeallocateOpacityTypeA( Opacity )

    TYPE(OpacityTypeA), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % GreyOpacity_Number_FD(i) % Values )
      DEALLOCATE( Opacity % GreyOpacity_Energy_FD(i) % Values )
      DEALLOCATE( Opacity % GreyMoment_Number_FD(i) % Values )
      DEALLOCATE( Opacity % GreyMoment_Energy_FD(i) % Values )
      DEALLOCATE( Opacity % Absorptivity(i) % Values )
    END DO

    DEALLOCATE( Opacity % GreyOpacity_Number_FD )
    DEALLOCATE( Opacity % GreyOpacity_Energy_FD )
    DEALLOCATE( Opacity % GreyMoment_Number_FD )
    DEALLOCATE( Opacity % GreyMoment_Energy_FD )
    DEALLOCATE( Opacity % Absorptivity )
    DEALLOCATE( Opacity % Offsets )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Species )
    DEALLOCATE( Opacity % Names )

  END SUBROUTINE DeallocateOpacityTypeA


  SUBROUTINE DescribeOpacityTypeA( Opacity )

    TYPE(OpacityTypeA), INTENT(in) :: Opacity

    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A4,A)') ' ', 'Opacity Type A'
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
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Units     = ', TRIM( Opacity % Units(i) )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % Absorptivity(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % Absorptivity(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Offset     = ', Opacity % Offsets(i)
      WRITE(*,*)
      WRITE(*,'(A8,A22,I3.3,A3,A)') &
        ' ', 'GreyOpacity_Number_FD(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,3I5.4)') &
        ' ', 'Shape     = ', SHAPE( Opacity % GreyOpacity_Number_FD(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % GreyOpacity_Number_FD(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % GreyOpacity_Number_FD(i) % Values ) 
      WRITE(*,*)
      WRITE(*,'(A8,A22,I3.3,A3,A)') &
        ' ', 'GreyOpacity_Energy_FD(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,3I5.4)') &
        ' ', 'Shape     = ', SHAPE( Opacity % GreyOpacity_Energy_FD(i)% Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % GreyOpacity_Energy_FD(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % GreyOpacity_Energy_FD(i) % Values )
      WRITE(*,*)
      WRITE(*,'(A8,A21,I3.3,A3,A)') &
        ' ', 'GreyMoment_Number_FD(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,3I5.4)') &
        ' ', 'Shape     = ', SHAPE( Opacity % GreyMoment_Number_FD(i)% Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % GreyMoment_Number_FD(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % GreyMoment_Number_FD(i) % Values )
      WRITE(*,*)
      WRITE(*,'(A8,A21,I3.3,A3,A)') &
        ' ', 'GreyMoment_Energy_FD(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,3I5.4)') &
        ' ', 'Shape     = ', SHAPE( Opacity % GreyMoment_Energy_FD(i)% Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % GreyMoment_Energy_FD(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % GreyMoment_Energy_FD(i) % Values )

    END DO
    WRITE(*,*)

  END SUBROUTINE DescribeOpacityTypeA


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
    ALLOCATE( Opacity % Offsets(nOpacities, nMoments) )
    ALLOCATE( Opacity % GreyOpacity_Number_FD(nOpacities) )
    ALLOCATE( Opacity % GreyOpacity_Energy_FD(nOpacities) )
!    ALLOCATE( Opacity % GreyMoment_Number_FD(nOpacities) )
!    ALLOCATE( Opacity % GreyMoment_Energy_FD(nOpacities) )
    ALLOCATE( Opacity % Kernel(nOpacities) )

    DO i = 1, nOpacities

      ALLOCATE( Opacity % GreyOpacity_Number_FD(i) % Values &
                  ( nPoints(2), nPoints(3), nPoints(4), nMoments ) )

      ALLOCATE( Opacity % GreyOpacity_Energy_FD(i) % Values &
                  ( nPoints(2), nPoints(3), nPoints(4), nMoments ) )

!      ALLOCATE( Opacity % GreyMoment_Energy_FD(i) % Values &
!                  ( nPoints(2), nPoints(3), nPoints(4), nMoments ) )

!      ALLOCATE( Opacity % GreyMoment_Number_FD(i) % Values &
!                  ( nPoints(2), nPoints(3), nPoints(4), nMoments ) )

      ALLOCATE( Opacity % Kernel(i) % Values &
                  ( nPoints(1), nPoints(2), nPoints(3), nPoints(4), &
                    nMoments) )

    END DO

  END SUBROUTINE AllocateOpacityTypeB


  SUBROUTINE DeallocateOpacityTypeB( Opacity )

    TYPE(OpacityTypeB), INTENT(inout) :: Opacity

    INTEGER :: i

    DO i = 1, Opacity % nOpacities
      DEALLOCATE( Opacity % GreyOpacity_Number_FD(i) % Values )
      DEALLOCATE( Opacity % GreyOpacity_Energy_FD(i) % Values )
!      DEALLOCATE( Opacity % GreyMoment_Number_FD(i) % Values )
!      DEALLOCATE( Opacity % GreyMoment_Energy_FD(i) % Values )
      DEALLOCATE( Opacity % Kernel(i) % Values )
    END DO

    DEALLOCATE( Opacity % GreyOpacity_Number_FD )
    DEALLOCATE( Opacity % GreyOpacity_Energy_FD )
!    DEALLOCATE( Opacity % GreyMoment_Number_FD )
!    DEALLOCATE( Opacity % GreyMoment_Energy_FD )
    DEALLOCATE( Opacity % Kernel )
    DEALLOCATE( Opacity % Offsets )
    DEALLOCATE( Opacity % Units )
    DEALLOCATE( Opacity % Species )
    DEALLOCATE( Opacity % Names )

  END SUBROUTINE DeallocateOpacityTypeB


  SUBROUTINE DescribeOpacityTypeB( Opacity )

    TYPE(OpacityTypeB), INTENT(in) :: Opacity

    INTEGER :: i, l

    WRITE(*,*)
    WRITE(*,'(A4,A)') ' ', 'Opacity Type B'
    WRITE(*,'(A4,A)') ' ', '--------------'
    WRITE(*,'(A6,A13,I3.3)') &
      ' ', 'nOpacities = ', Opacity % nOpacities
    WRITE(*,'(A6,A13,I3.3)') &
      ' ', 'nMoments   = ', Opacity % nMoments
    WRITE(*,'(A6,A13,5I5.4)') &
      ' ', 'nPoints    = ', Opacity % nPoints, Opacity % nMoments
    WRITE(*,'(A6,A13,I10.10)') &
      ' ', 'DOFs       = ', &
      Opacity % nOpacities * Opacity % nMoments &
        * PRODUCT( Opacity % nPoints )

    DO i = 1, Opacity % nOpacities
      WRITE(*,*)
      WRITE(*,'(A6,A8,I3.3,A3,A)') &
        ' ', 'Opacity(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Units     = ', TRIM( Opacity % Units(i) )

      DO l = 1, Opacity % nMoments
      WRITE(*,*)
         WRITE(*,'(A8,A16,I3.3)') &
           ' ', 'For Moments l = ', l
         WRITE(*,'(A8,A12,ES12.4E3)') &
           ' ', 'Min Value = ', MINVAL( Opacity % Kernel(i) % Values(:,:,:,:,l) )
         WRITE(*,'(A8,A12,ES12.4E3)') &
           ' ', 'Max Value = ', MAXVAL( Opacity % Kernel(i) % Values(:,:,:,:,l) )
         WRITE(*,'(A8,A12,ES12.4E3)') &
           ' ', 'Offset     = ', Opacity % Offsets(i,l)
      END DO ! l = nMoment

      WRITE(*,*)
      WRITE(*,'(A8,A22,I3.3,A3,A)') &
        ' ', 'GreyOpacity_Number_FD(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,4I5.4)') &
        ' ', 'Shape     = ', SHAPE( Opacity % GreyOpacity_Number_FD(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % GreyOpacity_Number_FD(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % GreyOpacity_Number_FD(i) % Values ) 
      WRITE(*,*)
      WRITE(*,'(A8,A22,I3.3,A3,A)') &
        ' ', 'GreyOpacity_Energy_FD(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,4I5.4)') &
        ' ', 'Shape     = ', SHAPE( Opacity % GreyOpacity_Energy_FD(i)% Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % GreyOpacity_Energy_FD(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % GreyOpacity_Energy_FD(i) % Values )
      WRITE(*,*)
!      WRITE(*,'(A8,A21,I3.3,A3,A)') &
!        ' ', 'GreyMoment_Number_FD(',i,'): ', TRIM( Opacity % Names(i) )
!      WRITE(*,'(A8,A12,A)') &
!        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
!      WRITE(*,'(A8,A12,4I5.4)') &
!        ' ', 'Shape     = ', SHAPE( Opacity % GreyMoment_Number_FD(i)% Values )
!      WRITE(*,'(A8,A12,ES12.4E3)') &
!        ' ', 'Min Value = ', MINVAL( Opacity % GreyMoment_Number_FD(i) % Values )
!      WRITE(*,'(A8,A12,ES12.4E3)') &
!        ' ', 'Max Value = ', MAXVAL( Opacity % GreyMoment_Number_FD(i) % Values )
!      WRITE(*,*)
!      WRITE(*,'(A8,A21,I3.3,A3,A)') &
!        ' ', 'GreyMoment_Energy_FD(',i,'): ', TRIM( Opacity % Names(i) )
!      WRITE(*,'(A8,A12,A)') &
!        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
!      WRITE(*,'(A8,A12,4I5.4)') &
!        ' ', 'Shape     = ', SHAPE( Opacity % GreyMoment_Energy_FD(i)% Values )
!      WRITE(*,'(A8,A12,ES12.4E3)') &
!        ' ', 'Min Value = ', MINVAL( Opacity % GreyMoment_Energy_FD(i) % Values )
!      WRITE(*,'(A8,A12,ES12.4E3)') &
!        ' ', 'Max Value = ', MAXVAL( Opacity % GreyMoment_Energy_FD(i) % Values )
    END DO
    WRITE(*,*)

  END SUBROUTINE DescribeOpacityTypeB


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


  SUBROUTINE DescribeOpacityTypeC( Opacity )

    TYPE(OpacityTypeC), INTENT(in) :: Opacity

    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A4,A)') ' ', 'Opacity Type C'
    WRITE(*,'(A4,A)') ' ', '--------------'
    WRITE(*,'(A6,A13,I3.3)') &
      ' ', 'nOpacities = ', Opacity % nOpacities
    WRITE(*,'(A6,A13,I3.3)') &
      ' ', 'nMoments   = ', Opacity % nMoments
    WRITE(*,'(A6,A13,4I5.4)') &
      ' ', 'nPoints    = ', Opacity % nPoints
    WRITE(*,'(A6,A13,I10.10)') &
      ' ', 'DOFs       = ', &
      Opacity % nOpacities * Opacity % nMoments &
        * Opacity % nPoints(1) * PRODUCT( Opacity % nPoints )

    DO i = 1, Opacity % nOpacities
      WRITE(*,*)
      WRITE(*,'(A6,A8,I3.3,A3,A)') &
        ' ', 'Opacity(',i,'): ', TRIM( Opacity % Names(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Species   = ', TRIM( Opacity % Species(i) )
      WRITE(*,'(A8,A12,A)') &
        ' ', 'Units     = ', TRIM( Opacity % Units(i) )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Min Value = ', MINVAL( Opacity % Kernel(i) % Values )
      WRITE(*,'(A8,A12,ES12.4E3)') &
        ' ', 'Max Value = ', MAXVAL( Opacity % Kernel(i) % Values )
    END DO
    WRITE(*,*)

  END SUBROUTINE DescribeOpacityTypeC


END MODULE wlOpacityFieldsModule

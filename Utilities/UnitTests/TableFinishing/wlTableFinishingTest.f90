PROGRAM wlTableFinishingTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtTableFinishingModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF,  & 
                           ReadEquationOfStateTableHDF, & 
                           WriteEquationOfStateTableHDF
  implicit none

  INTEGER  :: i, j, k, l, MinGradient
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp) :: Interpolant

  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE            :: fails     
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE            :: LoneCells    
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE            :: LinOkX 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE            :: LinOkY 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE            :: LinOkZ 

  LogInterp = (/.true.,.true.,.false./)

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "StandardResEquationOfStateTableFirstPass.h5" )

    WRITE (*,*) "Table Read", EOSTable % nPoints 

    ALLOCATE( fails( EOSTable % nPoints(1), EOSTable % nPoints(2),              &
                     EOSTable % nPoints(3) ), LoneCells( EOSTable % nPoints(1), &
                     EOSTable % nPoints(2), EOSTable % nPoints(3) ),            &
                     LinOkX(EOSTable % nPoints(1), EOSTable % nPoints(2),       &   
                     EOSTable % nPoints(3)), LinOkY(EOSTable % nPoints(1),      &
                     EOSTable % nPoints(2), EOSTable % nPoints(3)),             &
                     LinOkZ(EOSTable % nPoints(1), EOSTable % nPoints(2),       &   
                     EOSTable % nPoints(3))  )
    WRITE (*,*) "Logicals Allocated" 

    fails(:,:,:) = EOSTable % DV % Variables(1) % Values(:,:,:) <= 0.0d0 


  CALL HoleCharacterize( fails, LinOkX, LinOkY, LinOkZ )

    DO k = 1, SIZE(fails, DIM=3)
      DO j = 1, SIZE(fails, DIM=2)
        DO i = 1, SIZE(fails, DIM=1)

          IF ( .not.fails(i,j,k) ) CYCLE
          IF ( fails(i,j,k) ) THEN
            CALL GradientCheck( i, j, k, LinOkX, LinOkY, LinOkZ, &
                                EOSTable % DV % Variables(1) % Values(:,:,:) , MinGradient)            
          WRITE (*,*) i, j, k, MinGradient
            IF ( MinGradient == 0 ) THEN
              CYCLE
            ELSE
            CALL LogLineInterpolateSingleVariable( i, j, k, &
                   EOSTable % DV % Variables(1) % Values(:,:,:), MinGradient, Interpolant)
            EOSTable % DV % Variables(1) % Values(i,j,k) = Interpolant
            END IF
          WRITE (*,*) Interpolant 
          END IF
        END DO
      END DO
    END DO

!  STOP
!  CALL LoneCellLocate( fails, LoneCells )
!
!    DO k = 2, SIZE(LoneCells, DIM=3) - 1
!      DO j = 2, SIZE(LoneCells, DIM=2) - 1
!        DO i = 2, SIZE(LoneCells, DIM=1) - 1
!
!          IF ( LoneCells(i,j,k) ) THEN
!            DO l = 1, EOSTable % nVariables 
!              CALL LoneCellLogInterpolateSingleVariable( i, j, k,                 &
!                                    EOSTable % TS % States(1) % Values,           &
!                                    EOSTable % TS % States(2) % Values,           &
!                                    EOSTable % TS % States(3) % Values,           &
!                                    LogInterp,                                    &
!                                    EOSTable % DV % Variables(l) % Values(:,:,:), &
!                                      Interpolant )
!              EOSTable % DV % Variables(l) % Values(i,j,k) = Interpolant
!            END DO
!          ELSE
!            CYCLE
!          END IF

!        END DO
!      END DO
!    END DO

  CALL WriteEquationOfStateTableHDF( EOSTable )

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlTableFinishingTest

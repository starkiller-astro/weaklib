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

  INTEGER  :: i, j, k, l, idim
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iMinGradient
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp) :: Interpolant
  REAL(dp) :: delta

  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: fail     
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: Repaired 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: LoneCells    
  LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: LinearOK 

  LogInterp = (/.true.,.true.,.false./)
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "StandardResEquationOfStateTable.h5" )

  ASSOCIATE( nPoints => EOSTable % nPoints )

  WRITE (*,*) "Table Read", nPoints 

  ALLOCATE( fail( nPoints(1), nPoints(2), nPoints(3) ),        &
            Repaired( nPoints(1), nPoints(2), nPoints(3) ),    & 
            LoneCells( nPoints(1), nPoints(2), nPoints(3) ),   & 
            LinearOK(3, nPoints(1), nPoints(2), nPoints(3) ),  &
            iMinGradient( nPoints(1), nPoints(2), nPoints(3)) ) 
            
  
  END ASSOCIATE

  WRITE (*,*) "Logicals Allocated" 

  Repaired = .false.

  fail(:,:,:) = EOSTable % DV % Variables(1) % Values(:,:,:) <= 0.0d0 

  
  ASSOCIATE( Pressure => EOSTable % DV % Variables(1) % Values(:,:,:) )

  CALL HoleCharacterize( fail, LinearOK, Pressure, iMinGradient )

  DO k = 1, SIZE(fail, DIM=3)
    DO j = 1, SIZE(fail, DIM=2)
      DO i = 1, SIZE(fail, DIM=1)

        IF ( .not.fail(i,j,k) ) CYCLE

        WRITE (*,*) i, j, k, iMinGradient(i,j,k)

        IF ( iMinGradient(i,j,k) == 0 ) CYCLE 

        delta = 0.5d0  !REPLACE 

        idim = iMinGradient(i,j,k)

        CALL LogLineInterpolateSingleVariable &
               ( i, j, k, idim, delta, Pressure, Interpolant )

        Pressure(i,j,k) = Interpolant

        Repaired(i,j,k) = .true.

        WRITE (*,*) Interpolant 

      END DO
    END DO
  END DO

  END ASSOCIATE ! Pressure

!  STOP
!  CALL LoneCellLocate( fail, LoneCells )
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

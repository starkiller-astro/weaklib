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

  INTEGER  :: i, j, k, l, m, nHoles
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iMinGradient
  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: iLimits
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp) :: InterpolantFine
  REAL(dp), DIMENSION(:), ALLOCATABLE :: InterpolantCoarse
  REAL(dp) :: DeltaFine
  REAL(dp), DIMENSION(:), ALLOCATABLE :: DeltaCoarse

  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: Fail 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: Repaired 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: LoneCells    

  LogInterp = (/.true.,.true.,.false./)
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "HighResEquationOfStateTable.h5" )

  ASSOCIATE( nPoints => EOSTable % nPoints )

  WRITE (*,*) "Table Read", nPoints 

  ALLOCATE( Fail( nPoints(1), nPoints(2), nPoints(3) ),         &
            Repaired( nPoints(1), nPoints(2), nPoints(3) ),     & 
            LoneCells( nPoints(1), nPoints(2), nPoints(3) ),    & 
            iMinGradient( nPoints(1), nPoints(2), nPoints(3) ), &  
            iLimits( 2, nPoints(1), nPoints(2), nPoints(3) ) ) 
  
  END ASSOCIATE

  WRITE (*,*) "Logicals Allocated" 

  Repaired = .false.

  Fail(:,:,:) = EOSTable % DV % Variables(1) % Values(:,:,:) <= 0.0d0 
  
  ASSOCIATE( Pressure => EOSTable % DV % Variables(1) % Values(:,:,:) )

  ! Find dimension (iMinGradient) with smallest gradient 
  !   for interpolation across single cell hole 
  
  CALL HoleCharacterizeFine( Fail, Pressure, iMinGradient )

  DO k = 1, SIZE(Fail, DIM=3)
    DO j = 1, SIZE(Fail, DIM=2)
      DO i = 1, SIZE(Fail, DIM=1)

        IF ( .not.Fail(i,j,k) ) CYCLE

        WRITE (*,*) i, j, k, iMinGradient(i,j,k)

        IF ( iMinGradient(i,j,k) == 0 ) CYCLE 

        DeltaFine = 0.5d0  !REPLACE 

        CALL LogInterpolateFine1D&
               ( i, j, k, iMinGradient(i,j,k), DeltaFine, Pressure, InterpolantFine )

        Pressure(i,j,k) = InterpolantFine

        Repaired(i,j,k) = .true.

        WRITE (*,*) InterpolantFine 

      END DO
    END DO
  END DO

  CALL HoleCharacterizeCoarse( Fail, Repaired, Pressure, iMinGradient, iLimits ) 

  DO k = 1, SIZE(Fail, DIM=3)
    DO j = 1, SIZE(Fail, DIM=2)
      DO i = 1, SIZE(Fail, DIM=1)

        IF ( .not.Fail(i,j,k) .or. Repaired(i,j,k) ) CYCLE

        WRITE (*,*) i, j, k, iMinGradient(i,j,k)

        IF ( iMinGradient(i,j,k) == 0 ) CYCLE

        nHoles = ( iLimits(2,i,j,k) - iLimits(1,i,j,k) - 1 )
          WRITE (*,*) "nHoles=", nHoles 
        WRITE (*,*) "iLimits =", iLimits(1,i,j,k), iLimits(2,i,j,k), "nHoles=", nHoles
        ALLOCATE( DeltaCoarse( nHoles ), InterpolantCoarse( nHoles ) )

        DO m = 1, nHoles 
          DeltaCoarse(m) = (m)/(nHoles+1.d0) 
          WRITE (*,*) "DC=", DeltaCoarse(m)
        END DO

        CALL LogInterpolateCoarse1D&
               ( i, j, k, iMinGradient(i,j,k), iLimits, DeltaCoarse, Pressure, InterpolantCoarse )
        
        DO l = 1, nHoles

          SELECT CASE( iMinGradient(i,j,k) )

            CASE(1)
              Pressure(i+l-1,j,k) = InterpolantCoarse(l)

              Repaired(i+l-1,j,k) = .true.

            CASE(2)
              Pressure(i,j+l-1,k) = InterpolantCoarse(l)

              Repaired(i,j+l-1,k) = .true.

            CASE(3)
              Pressure(i,j,k+l-1) = InterpolantCoarse(l)

              Repaired(i,j,k+l-1) = .true.

          END SELECT
        END DO

        WRITE (*,*) InterpolantCoarse

        DEALLOCATE( DeltaCoarse, InterpolantCoarse )
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

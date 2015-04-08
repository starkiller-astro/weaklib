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

  INTEGER  :: i, j, k, l
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iMinGradient
  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: iLimits
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp) :: InterpolantFine
  REAL(dp) :: InterpolantCoarse
  REAL(dp) :: DeltaFine
  REAL(dp) :: DeltaCoarse

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

        WRITE (*,*) "iLimits =", iLimits(1,i,j,k), iLimits(2,i,j,k)

        SELECT CASE( iMinGradient(i,j,k) )
          CASE(1)
            DeltaCoarse = DBLE( i - iLimits(1,i,j,k) ) & 
                            / DBLE( iLimits(2,i,j,k) - iLimits(1,i,j,k) )
          CASE(2)
            DeltaCoarse = DBLE( j - iLimits(1,i,j,k) ) & 
                            / DBLE( iLimits(2,i,j,k) - iLimits(1,i,j,k) )
          CASE(3)
            DeltaCoarse = DBLE( k - iLimits(1,i,j,k) ) & 
                            / DBLE( iLimits(2,i,j,k) - iLimits(1,i,j,k) )
        END SELECT

        CALL LogInterpolateCoarse1D&
               ( i, j, k, iMinGradient(i,j,k), iLimits, DeltaCoarse, Pressure, InterpolantCoarse )

        Pressure(i,j,k) = InterpolantCoarse
        ! Add other DV's
        Repaired(i,j,k) = .true.

        WRITE (*,*) InterpolantCoarse

      END DO
    END DO
  END DO

  END ASSOCIATE ! Pressure

  CALL WriteEquationOfStateTableHDF( EOSTable )

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlTableFinishingTest

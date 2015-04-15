PROGRAM wlTableFinishingTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtTableFinishingModule
  USE wlExtNumericalModule, ONLY: epsilon, zero
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF,  & 
                           ReadEquationOfStateTableHDF, & 
                           WriteEquationOfStateTableHDF
  implicit none

  INTEGER  :: i, j, k, l
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iMinGradient
  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: iLimits
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp) :: InterpolantFine, InterpolantCoarse
  REAL(dp) :: DeltaFine, DeltaCoarse
  REAL(dp) :: minvar
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: Fail 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: Repaired 
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE   :: LoneCells    

  LogInterp = (/.true.,.true.,.false./)
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "StandardResEquationOfStateTable.h5" )

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
  
  !ASSOCIATE( Pressure => EOSTable % DV % Variables(1) % Values(:,:,:) )

  ! Find dimension (iMinGradient) with smallest gradient 
  !   for interpolation across single cell hole 
  
  CALL HoleCharacterizeFine( Fail, EOSTable % DV % Variables(1) % Values, iMinGradient )

  DO k = 1, SIZE(Fail, DIM=3)
    DO j = 1, SIZE(Fail, DIM=2)
      DO i = 1, SIZE(Fail, DIM=1)

        IF ( .not.Fail(i,j,k) ) CYCLE

        WRITE (*,*) i, j, k, iMinGradient(i,j,k)

        IF ( iMinGradient(i,j,k) == 0 ) CYCLE 

        DeltaFine = 0.5d0 

        DO l = 1, EOSTable % nVariables

          CALL LogInterpolateFine1D&
                 ( i, j, k, iMinGradient(i,j,k), DeltaFine, &
                 EOSTable % DV % Variables(l) % Values, InterpolantFine )

          EOSTable % DV % Variables(l) % Values(i,j,k) = InterpolantFine

          Repaired(i,j,k) = .true.

        END DO

        WRITE (*,*) InterpolantFine 

      END DO
    END DO
  END DO

  CALL HoleCharacterizeCoarse( Fail, Repaired, &
         EOSTable % DV % Variables(1) % Values, iMinGradient, iLimits ) 

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

        DO l = 1, EOSTable % nVariables
          CALL LogInterpolateCoarse1D&
                 ( i, j, k, iMinGradient(i,j,k), iLimits, DeltaCoarse, &
                  EOSTable % DV % Variables(l) % Values, InterpolantCoarse )

          EOSTable % DV % Variables(l) % Values(i,j,k) = InterpolantCoarse
          Repaired(i,j,k) = .true.

        WRITE (*,*) InterpolantCoarse
        END DO

      END DO
    END DO
  END DO

  DO l = 1, EOSTable % nVariables
    WRITE (*,*) EOSTable % DV % Names(l)
    minvar = MINVAL( EOSTable % DV % Variables(l) % Values )
    WRITE (*,*) "minvar=", minvar
    EOSTable % DV % Offsets(l) = -2.d0 * MIN( 0.d0, minvar )
    WRITE (*,*) "Offset=", EOSTable % DV % Offsets(l)
    EOSTable % DV % Variables(l) % Values &
      = LOG10( EOSTable % DV % Variables(l) % Values &
               + EOSTable % DV % Offsets(l) + epsilon )        

  END DO

  CALL WriteEquationOfStateTableHDF( EOSTable )

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlTableFinishingTest

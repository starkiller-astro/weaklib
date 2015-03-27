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

  CALL ReadEquationOfStateTableHDF( EOSTable, "LowResEquationOfStateTable.h5" )

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

    !WRITE (*,*) "fails =", fails 

  CALL HoleCharacterize( fails, LinOkX, LinOkY, LinOkZ )

    !WRITE (*,*) "fails =", fails 
    !WRITE (*,*) LinOkY 
    !WRITE (*,*) "that was LinOkY" 

  STOP
  CALL LoneCellLocate( fails, LoneCells )

      DO k = 2, SIZE(LoneCells, DIM=3) - 1
        DO j = 2, SIZE(LoneCells, DIM=2) - 1
          DO i = 2, SIZE(LoneCells, DIM=1) - 1

            IF ( LoneCells(i,j,k) ) THEN
              DO l = 1, EOSTable % nVariables 
                CALL LoneCellLogInterpolateSingleVariable( i, j, k,                 &
                                      EOSTable % TS % States(1) % Values,           &
                                      EOSTable % TS % States(2) % Values,           &
                                      EOSTable % TS % States(3) % Values,           &
                                      LogInterp,                                    &
                                      EOSTable % DV % Variables(l) % Values(:,:,:), &
                                      Interpolant )
                EOSTable % DV % Variables(l) % Values(i,j,k) = Interpolant
              END DO
              WRITE (*,*) i, j, k
              WRITE (*,*) Interpolant 
              WRITE (*,*) EOSTable % DV % Variables(1) % Values(i,j,k)
            ELSE
              CYCLE
            END IF

          END DO
        END DO
      END DO
    !WRITE (*,*) EOSTable % DV % Variables(1) % Values  

  CALL WriteEquationOfStateTableHDF( EOSTable )

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlTableFinishingTest

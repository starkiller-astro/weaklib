PROGRAM wlTableFinishingTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtTableFinishingModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF, & 
                           ReadEquationOfStateTableHDF 

  implicit none

  INTEGER  :: i
  TYPE(EquationOfStateTableType) :: EOSTable


  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE            :: fails     
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE            :: LoneCells    



  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "LowResEquationOfStateTable.h5" )
    WRITE (*,*) "Table Read", EOSTable % nPoints 
    ALLOCATE( fails( EOSTable % nPoints(1), EOSTable % nPoints(2), &
                     EOSTable % nPoints(3) ), LoneCells( EOSTable % nPoints(1), &
                     EOSTable % nPoints(2), EOSTable % nPoints(3)  ) )
    WRITE (*,*) "Logicals Allocated" 
  fails(:,:,:) = EOSTable % DV % Variables(1) % Values(:,:,:) <= 0.0d0 

    WRITE (*,*) "fails =", fails 
  CALL LoneCellLocate( fails, LoneCells)

  WRITE (*,*) LoneCells

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlTableFinishingTest

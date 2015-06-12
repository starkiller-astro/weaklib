PROGRAM wlInversionTest

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
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp) :: E_Internal, Temperature
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Energy_Table, Temp_Table 

  !LogInterp = (/.true.,.true.,.false./)
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  ASSOCIATE( nPoints => EOSTable % nPoints )

  WRITE (*,*) "Table Read", nPoints 

  ALLOCATE( Energy_Table( nPoints(2) ), Temp_Table( nPoints(2) ) )
  
  WRITE (*,*) "Allocation Complete" 

  END ASSOCIATE
  
  ! Associate call variables with appropriate EOSTable-type structures
  ! Make any necessary dimensional handling changes to energy table
  ! Should I start with a "Dimension(:)" table, so that "locate" Subroutine 
  ! Can handle the table? DONE 

  !DO i=1, SIZE(nPoints(2))
  !  Energy_Table(i) = 10**( EOSTable % DV % Variables(3) % Values(10,i,10) ) 
  !  Temp_Table(i) = 10**( EOSTable % TS % States(2) % Values(i) ) 
  !END DO
  ! Internal energy for (10,:,10) in the standard res table runs from
  ! 8.22e32 to 1.1e42 
  E_Internal = 1.0d19  

  WRITE (*,*) "Internal Energy =", E_Internal

  Energy_Table = 10**( EOSTable % DV % Variables(3) % Values(30,:,5) ) 

  Temp_Table = ( EOSTable % TS % States(2) % Values ) 

  CALL TemperatureFinder( E_Internal, Energy_Table, Temp_table, &
                          EOSTable % DV % Offsets(3), Temperature )  

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInversionTest

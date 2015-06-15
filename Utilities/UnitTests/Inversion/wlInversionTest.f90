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
  !LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Energy_Table, Temp_Table, E_Internal, &
                                         Temperature
  !LogInterp = (/.true.,.true.,.false./)
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  ASSOCIATE( nPoints => EOSTable % nPoints )

  WRITE (*,*) "Table Read", nPoints 

  ALLOCATE( Energy_Table( nPoints(2) ), Temp_Table( nPoints(2) ),  &
            E_Internal( nPoints(2) - 1 ), Temperature( nPoints(2) - 1 ) ) 
  
  WRITE (*,*) "Allocation Complete" 

  END ASSOCIATE

  ! Internal energy for (10,:,10) in the standard res table runs from
  ! 8.22e32 to 1.1e42 


  Energy_Table = 10**( EOSTable % DV % Variables(3) % Values(25,:,25) ) 

  DO i = 1, SIZE(Energy_Table) - 1
  E_Internal(i) = Energy_Table(i) * 1.05d0  
  END DO

  Temp_Table = ( EOSTable % TS % States(2) % Values ) 

  CALL TemperatureFinder( E_Internal, Energy_Table, Temp_Table, &
                          EOSTable % DV % Offsets(3), Temperature )  
  DO i = 1, SIZE(E_Internal)
    WRITE (*,*) i, E_Internal(i), Temperature(i), Temp_Table(i), ( Temperature(i) - Temp_Table(i) ) / Temp_Table(i)
  END DO

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInversionTest

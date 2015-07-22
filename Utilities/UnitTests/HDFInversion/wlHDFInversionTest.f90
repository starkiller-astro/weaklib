PROGRAM wlHDFInversionTest

!  The goal of this test is to read a CHIMERA-slice HDF file (only the NSE zones)
!  and then to read the internal energy density of each zone, and then to convert 
!  (via inversion) that value into a temperature, and then to compare with the 
!  local temperature, both to determine table accuracy capability, but also to 
!  examine the effective phase space used by CHIMERA
!--------------------------------------------------------------------------------

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtTableFinishingModule
  USE wlExtNumericalModule, ONLY: epsilon, zero
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF,  & 
                           ReadEquationOfStateTableHDF, & 
                           WriteEquationOfStateTableHDF,&
                           ReadCHIMERAHDF
  implicit none

  INTEGER  :: i, j, k, l, TestUnit, ErrorUnit
  REAL(dp) :: maxnorm, L1Norm !, Temperature
  REAL(dp), DIMENSION(1) :: Temperature, Energy
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Norm
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Energy_Table

  REAL(dp), DIMENSION(1,240,722) :: Rho
  REAL(dp), DIMENSION(1,240,722) :: T
  REAL(dp), DIMENSION(1,240,722) :: Ye
  REAL(dp), DIMENSION(1,240,722) :: E_Int
  REAL(dp), DIMENSION(1,240,722) :: Entropy
  INTEGER, DIMENSION(1,240,722) :: NSE

!  nPoints = (/81,24,24/) ! Low Res
!  nPoints = (/81,500,24/) ! High Res in T only
!  nPoints = (/161,47,47/) ! Standard Res
!  nPoints = (/321,93,93/) ! High Res

  LogInterp = (/.true.,.true.,.false./)
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  WRITE (*,*) "Table read"

  OPEN( newunit = TestUnit, FILE="HDFInversionTableMap.d")
  OPEN( newunit = ErrorUnit, FILE="HDFInversionErrors.d")

  ASSOCIATE( nPoints   => EOSTable % nPoints,                 &
             TableRho  => EOSTable % TS % States(1) % Values, &
             TableTemp => EOSTable % TS % States(2) % Values, &
             TableYe   => EOSTable % TS % States(3) % Values  )


  ALLOCATE( Energy_Table( nPoints(2) ),      &
            Norm(1, 240, 722) )  

  WRITE (*,*) "Allocation Complete"
  

  CALL ReadCHIMERAHDF( Rho, T, Ye, E_Int, Entropy, NSE, &
                       "chimera_000000000_grid_1_01.h5") 

  WRITE (*,*) "Chimera read Complete"

  Norm = 0.0d0 

  maxnorm = 0.0d0

  DO i = 1, 240
    DO j = 1, 722
      IF ( NSE(1,i,j) == 0 ) THEN
        CYCLE
      END IF 
     
      CALL locate( TableRho, nPoints(1), Rho(1,i,j), k )

      CALL locate( TableYe, nPoints(3), Ye(1,i,j), l )


  WRITE (*,*) "Locate Complete", i, j, k, l, Rho(1,i,j)
!    CALL LogInterpolateSingleVariable &
!           ( Rho(1,i,j),        &
!             T(1,i,j),                                   &
!             Ye(1,i,j),     &
!             EOSTable % TS % States(1) % Values,           &
!             EOSTable % TS % States(2) % Values,           &
!             EOSTable % TS % States(3) % Values,           &
!             LogInterp,                                    &
!             EOSTable % DV % Offsets(3),                   &
!             EOSTable % DV % Variables(3) % Values(:,:,:), &
!             E_Internal )

      Energy_Table = 10.d0**( EOSTable % DV % Variables(3) % Values(k,:,l) ) &
                     - EOSTable % DV % Offsets(3) 

      Energy(1) = E_Int(1,i,j)

  WRITE (*,*) "Energy = ", Energy 

      CALL ComputeTempFromIntEnergy( Energy, Energy_Table, TableTemp, &
                          EOSTable % DV % Offsets(3), Temperature )  
  WRITE (*,*) "Temp = ", Temperature 

      Norm(1,i,j) = ABS( Temperature(1) - T(1,i,j) ) / T(1,i,j) 
      WRITE (ErrorUnit,'(2i4, 4es12.5)') i, j, E_Int(1,i,j), &
               Temperature, T(1,i,j), Norm(1,i,j)
      IF ( Norm(1,i,j) > maxnorm) THEN 
        maxnorm = Norm(1,i,j)
      END IF
    END DO
  END DO

 ! L1Norm = SUM( ABS( Norm( 1:nPoints(1), 1:( nPoints(2) - 1), 1:nPoints(3) - 1 ) ) ) &
  !           /( nPoints(1) * (nPoints(2) - 1 ) * ( nPoints(3) - 1 ) ) 

  END ASSOCIATE

  !WRITE (ErrorUnit,*) "L1Norm/N=", L1Norm 
  !WRITE (ErrorUnit,*) "maxnorm =", maxnorm 

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlHDFInversionTest

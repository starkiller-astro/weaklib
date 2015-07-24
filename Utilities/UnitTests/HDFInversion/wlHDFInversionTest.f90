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

  INTEGER  :: i, j, k, kk, l, TestUnit, ErrorUnit, ZoneLimit
  REAL(dp) :: maxnorm, L1Norm
  REAL(dp), DIMENSION(1) :: Temperature, Energy
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Norm
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Energy_Table, TestRho, TestYe

  REAL(dp), DIMENSION(722,240,1) :: Rho
  REAL(dp), DIMENSION(722,240,1) :: T
  REAL(dp), DIMENSION(722,240,1) :: Ye
  REAL(dp), DIMENSION(722,240,1) :: E_Int
  REAL(dp), DIMENSION(722,240,1) :: Entropy
  INTEGER, DIMENSION(722,240,1) :: NSE

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


  ALLOCATE( Energy_Table( nPoints(2) ), TestRho( nPoints(2) ), &
            Norm(722,240,1), TestYe( nPoints(2) ) )  

  WRITE (*,*) "Allocation Complete"
  

  CALL ReadCHIMERAHDF( Rho, T, Ye, E_Int, Entropy, NSE, &
                       "chimera_000000000_grid_1_01.h5") 

  WRITE (*,*) "Chimera read Complete"

  Norm = 0.0d0 

  maxnorm = 0.0d0
  
  k = 1
    DO j = 1, 240
      DO i = 1, 720 

      IF ( NSE(i,j,k) == 0 ) THEN
        CYCLE
      END IF 
     
      CALL locate( TableRho, nPoints(1), Rho(i,j,k), kk )

      IF ( kk  == 0 ) THEN
        CYCLE
      END IF 

      CALL locate( TableYe, nPoints(3), Ye(i,j,k), l )

      IF ( l == 0 ) THEN
        CYCLE
      END IF 


      !WRITE (*,*) "Locate Complete", i, j, kk, l, Rho(i,j,k)

      TestRho = Rho(i,j,k)
      TestYe = Ye(i,j,k)

      CALL LogInterpolateSingleVariable &
             ( TestRho,        &
               TableTemp,                                   &
               TestYe,     &
               EOSTable % TS % States(1) % Values,           &
               EOSTable % TS % States(2) % Values,           &
               EOSTable % TS % States(3) % Values,           &
               LogInterp,                                    &
               EOSTable % DV % Offsets(3),                   &
               EOSTable % DV % Variables(3) % Values(:,:,:), &
               Energy_Table )

!      Energy_Table = 10.d0**( EOSTable % DV % Variables(3) % Values(kk,:,l) ) &
!                     - EOSTable % DV % Offsets(3) 

      Energy(1) = E_Int(i,j,k)

      WRITE (*,*) "Energy = ", Energy 

      CALL ComputeTempFromIntEnergy( Energy, Energy_Table, TableTemp, &
                          EOSTable % DV % Offsets(3), Temperature )  
      WRITE (*,*) "Temp = ", Temperature 

      IF ( Temperature(1) < 1.d0  ) THEN
        CYCLE
      END IF

      Norm(i,j,k) = ABS( Temperature(1) - T(i,j,k) ) / T(i,j,k) 

      WRITE (*,*) "Norm = ", Norm(i,j,k) 

      WRITE (ErrorUnit,'(2i4, 6es12.5, 3i4)') i, j, E_Int(i,j,k), &
               Temperature, T(i,j,k), Norm(i,j,k), &
               MAXVAL(Energy_Table), MINVAL(Energy_Table), NSE(i,j,k), kk, l

      IF ( Norm(i,j,k) > maxnorm) THEN 

        maxnorm = Norm(i,j,k)

      END IF
    END DO
  END DO

  !L1Norm = SUM( ABS( Norm( 1:720 , 1:240, 1) ) ) &
  !           /( 720 * 240 ) 

  !WRITE (ErrorUnit,*) L1Norm, maxnorm 
  END ASSOCIATE

  !WRITE (ErrorUnit,*) "L1Norm/N=", L1Norm 
  !WRITE (ErrorUnit,*) "maxnorm =", maxnorm 

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlHDFInversionTest

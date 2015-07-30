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

  INTEGER  :: i, j, k, kk, l, TestUnit, ErrorUnit, ZoneLimit, nx, ny, nz, imax
  REAL(dp) :: maxnorm, L1Norm
  REAL(dp), DIMENSION(1) :: Temperature, Energy
  TYPE(EquationOfStateTableType) :: EOSTable
  LOGICAL, DIMENSION(3) :: LogInterp
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Norm
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Energy_Table, TestRho, TestT, TestYe

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Rho
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ye
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: E_Int
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Entropy
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NSE

!  nPoints = (/81,24,24/) ! Low Res
!  nPoints = (/81,500,24/) ! High Res in T only
!  nPoints = (/161,47,47/) ! Standard Res
!  nPoints = (/321,93,93/) ! High Res

  LogInterp = (/.true.,.true.,.false./)

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  OPEN( newunit = TestUnit, FILE="HDFInversionTableMap.d")
  OPEN( newunit = ErrorUnit, FILE="HDFInversionErrors.d")

  ASSOCIATE( nPoints   => EOSTable % nPoints,                 &
             TableRho  => EOSTable % TS % States(1) % Values, &
             TableTemp => EOSTable % TS % States(2) % Values, &
             TableYe   => EOSTable % TS % States(3) % Values  )
  
  CALL ReadCHIMERAHDF( Rho, T, Ye, E_Int, Entropy, NSE, imax, nx, ny, nz, &
                       "chimera_000018000_grid_1_01.h5") 

  ALLOCATE( Energy_Table( nPoints(2) ), TestRho( nPoints(2) ), &
            Norm(nx,ny,nz), TestYe( nPoints(2) ) )  

  WRITE (*,*) "imax, nx, ny, nz =", imax, nx, ny, nz
  WRITE (*,*) "Shape of Rho=", SHAPE(RHO)

  Norm = 0.0d0 

  maxnorm = 0.0d0
  
  DO k = 1, nz
    DO j = 1, ny 
      DO i = 1, imax 

      IF ( NSE(i,j,k) == 0 ) THEN
        CYCLE
      END IF 
     
      CALL locate( TableRho, nPoints(1), Rho(i,j,k), kk )

      IF ( kk  == 0 ) THEN
        WRITE (*,*) "Zone below density limit", i, j, k
        CYCLE
      END IF 

      CALL locate( TableYe, nPoints(3), Ye(i,j,k), l )

      IF ( l == 0 ) THEN
        WRITE (*,*) "Zone below ye limit", i, j, k
        CYCLE
      END IF 


      TestRho = Rho(i,j,k)
      !TestT = T(i,j,k)
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

      Energy(1) = E_Int(i,j,k)


      CALL ComputeTempFromIntEnergy( Energy, Energy_Table, TableTemp, &
                          EOSTable % DV % Offsets(3), Temperature )  

      IF ( Temperature(1) < 1.d0  ) THEN
        CYCLE
      END IF

      Norm(i,j,k) = ABS( Temperature(1) - T(i,j,k) ) / T(i,j,k) 

      WRITE (ErrorUnit,'(2i4, 4es12.5, i4)') i, j, E_Int(i,j,k), &
               Temperature, T(i,j,k), Norm(i,j,k), &
               NSE(i,j,k)

      IF ( Norm(i,j,k) > maxnorm) THEN 

        maxnorm = Norm(i,j,k)

      END IF
      END DO
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

PROGRAM wlInterpolationTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF, & 
                           ReadEquationOfStateTableHDF 

  implicit none

  INTEGER  :: i
  REAL(dp) :: rho
  REAL(dp) :: T  
  REAL(dp) :: Ye  
  INTEGER :: RhoLocation
  INTEGER :: TLocation
  INTEGER :: YeLocation
  LOGICAL, DIMENSION(3) :: LogInterp
  TYPE(EquationOfStateTableType) :: EOSTable


  REAL(dp) :: Interpolant

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "SmallEquationOfStateTable.h5" )

  LogInterp = (/.true.,.true.,.false./) 
  
  rho = ( EOSTable % TS % States(1) % Values(6) -         &
          EOSTable % TS % States(1) % Values(5) )/2.0d0 + &
          EOSTable % TS % States(1) % Values(5) 

  T = ( EOSTable % TS % States(2) % Values(6) -         &
        EOSTable % TS % States(2) % Values(5) )/2.0d0 + &
        EOSTable % TS % States(2) % Values(5) 

  Ye = ( EOSTable % TS % States(3) % Values(6) -         &
         EOSTable % TS % States(3) % Values(5) )/2.0d0 + &
         EOSTable % TS % States(3) % Values(5) 

  CALL LogInterpolateSingleVariable( rho, T, Ye,                                   &
                                     EOSTable % TS % States(1) % Values,           &
                                     EOSTable % TS % States(2) % Values,           &
                                     EOSTable % TS % States(3) % Values,           &
                                     LogInterp,                                    &
                                     EOSTable % DV % Variables(1) % Values(:,:,:), & 
                                     Interpolant )
   WRITE (*,*) Interpolant


  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInterpolationTest

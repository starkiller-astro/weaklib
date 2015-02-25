PROGRAM wlInterpolationTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF, & 
                           ReadEquationOfStateTableHDF 

  implicit none

  INTEGER  :: i
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  INTEGER :: RhoLocation
  INTEGER :: TLocation
  INTEGER :: YeLocation
  LOGICAL, DIMENSION(3) :: LogInterp
  TYPE(EquationOfStateTableType) :: EOSTable


  REAL(dp), DIMENSION(:), ALLOCATABLE :: Interpolant
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rand
  REAL(dp) :: Yemin, Yemax, logTmin, logTmax, logrhomin, logrhomax
  INTEGER :: i
  REAL(dp) :: LogT, logrho

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "SmallEquationOfStateTable.h5" )

  LogInterp = (/.true.,.true.,.false./)

  ALLOCATE( rho(10), T(10), Ye(10), rand(10,3), Interpolant(10) ) 

  CALL RANDOM_SEED( )
  CALL RANDOM_NUMBER( rand )
 
  Yemin = EOSTable % TS % minValues(3)
  Yemax = EOSTable % TS % maxValues(3)
  logTmin = LOG10( EOSTable % TS % minValues(2) )
  logTmax = LOG10( EOSTable % TS % maxValues(2) )
  logrhomin = LOG10( EOSTable % TS % minValues(1) )
  logrhomax = LOG10( EOSTable % TS % maxValues(1) )
  DO i=1,10
    Ye(i) = (Yemax - Yemin ) * rand(i,3) + Yemin   
    logT = (logTmax - logTmin ) * rand(i,2) + logTmin   
    T(i) = 10.d0**logT 
    logrho = (logrhomax - logrhomin ) * rand(i,1) + logrhomin   
    rho(i) = 10.d0**logrho 
  END DO 

  WRITE (*,*) "rho=", rho
  WRITE (*,*) "T=", T 
  WRITE (*,*) "Ye=", Ye 
  
!  rho = ( EOSTable % TS % States(1) % Values(6) -         &
!          EOSTable % TS % States(1) % Values(5) )/2.0d0 + &
!          EOSTable % TS % States(1) % Values(5) 

!  T = ( EOSTable % TS % States(2) % Values(6) -         &
!        EOSTable % TS % States(2) % Values(5) )/2.0d0 + &
!        EOSTable % TS % States(2) % Values(5) 

!  Ye = ( EOSTable % TS % States(3) % Values(6) -         &
!         EOSTable % TS % States(3) % Values(5) )/2.0d0 + &
!         EOSTable % TS % States(3) % Values(5) 

   EOSTable % DV % Variables(1) % Values(:,:,:) &
           = LOG10( EOSTable % DV % Variables(1) % Values(:,:,:) ) 

  CALL LogInterpolateSingleVariable( rho, T, Ye,                                   &
                                     EOSTable % TS % States(1) % Values,           &
                                     EOSTable % TS % States(2) % Values,           &
                                     EOSTable % TS % States(3) % Values,           &
                                     LogInterp,                                    &
                                     EOSTable % DV % Variables(1) % Values(:,:,:), & 
                                     Interpolant )
   WRITE (*,*) "Pressure =", Interpolant


  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInterpolationTest

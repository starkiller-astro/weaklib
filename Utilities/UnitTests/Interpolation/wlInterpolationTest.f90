PROGRAM wlInterpolationTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
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
  REAL(dp), DIMENSION(:), ALLOCATABLE :: DirectCall 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: press 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: entrop
  REAL(dp), DIMENSION(:), ALLOCATABLE :: energ
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_e
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_p
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_n 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn_prot 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn_neut 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn_alpha 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: z_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: a_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: be_heavy 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rand
  REAL(dp) :: Yemin, Yemax, logTmin, logTmax, logrhomin, logrhomax
  INTEGER :: i
  REAL(dp) :: LogT, logrho
  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  LOGICAL            :: fail        ! did EoS fail to converge

  LScompress = '220'
  LSFilePath = '../../../External/LS/Data'
  EOSFlag = "L"

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "SmallEquationOfStateTable.h5" )

  LogInterp = (/.true.,.true.,.false./)

  ALLOCATE( rho(10), T(10), Ye(10), rand(10,3), Interpolant(10), DirectCall(10),  &
            press(10), entrop(10), energ(10), chem_e(10), chem_p(10), chem_n(10), &
            xn_prot(10), xn_neut(10), xn_alpha(10), xn_heavy(10), z_heavy(10),    &
            a_heavy(10), be_heavy(10) ) 

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


   EOSTable % DV % Variables(1) % Values(:,:,:) &
           = LOG10( EOSTable % DV % Variables(1) % Values(:,:,:) ) 

  CALL LogInterpolateSingleVariable( rho, T, Ye,                                   &
                                     EOSTable % TS % States(1) % Values,           &
                                     EOSTable % TS % States(2) % Values,           &
                                     EOSTable % TS % States(3) % Values,           &
                                     LogInterp,                                    &
                                     EOSTable % DV % Variables(1) % Values(:,:,:), & 
                                     Interpolant )

  CALL wlExtInitializeEOS( LSFilePath, LScompress )

  DO i = 1, SIZE( rho )
    CALL wlGetFullEOS( rho(i), T(i), Ye(i), EOSFlag, fail, press(i), energ(i), &
                       entrop(i), chem_n(i), chem_p(i), chem_e(i), xn_neut(i), &
                       xn_prot(i), xn_alpha(i), xn_heavy(i), a_heavy(i),       &
                       z_heavy(i), be_heavy(i) ) 

    DirectCall(i) = press(i)
  END DO


  WRITE (*, '(4A22)' ) "rho=", "T=", "Ye=", "Pressure=" 
  DO i=1,10 
    WRITE (*, '(4E)' ) rho(i), T(i), Ye(i), Interpolant(i) 
  END DO

  WRITE (*, '(3A25)' ) "Interpolated Pressure= ", "Direct Call Pressure=", "Residual=" 
  DO i=1,10 
    WRITE (*, '(3E)' ) Interpolant(i), DirectCall(i), ABS( Interpolant(i) - DirectCall(i) ) 
  END DO

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInterpolationTest

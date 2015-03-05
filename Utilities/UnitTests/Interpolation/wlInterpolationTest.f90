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
  INTEGER, PARAMETER :: NumPoints = 100
  REAL(dp) :: LogT, logrho, L1norm, maxnorm 
  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  LOGICAL            :: fail        ! did EoS fail to converge

  99 FORMAT ("rho=", es12.5, " T=", es12.5, " Ye=" , es12.5, 3(1x,es12.5) )  


  LScompress = '220'
  LSFilePath = '../../../External/LS/Data'
  EOSFlag = "L"

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "LargeEquationOfStateTable.h5" )

  LogInterp = (/.true.,.true.,.false./)

  ALLOCATE( rho( NumPoints ), T( NumPoints ), Ye( NumPoints ), rand( NumPoints, 3 ), &
            Interpolant( NumPoints ), DirectCall( NumPoints ), press( NumPoints ),   &
            entrop( NumPoints ), energ( NumPoints ), chem_e( NumPoints ),            & 
            chem_p( NumPoints ), chem_n( NumPoints ), xn_prot( NumPoints ),          &
            xn_neut( NumPoints ), xn_alpha( NumPoints ), xn_heavy( NumPoints ),      &
            z_heavy( NumPoints ), a_heavy( NumPoints ), be_heavy( NumPoints ) ) 

  CALL RANDOM_SEED( )
  CALL RANDOM_NUMBER( rand )
 
  Yemin = EOSTable % TS % minValues(3)
  Yemax = EOSTable % TS % maxValues(3)
  logTmin = LOG10( EOSTable % TS % minValues(2) )
  logTmax = LOG10( EOSTable % TS % maxValues(2) )
  logrhomin = LOG10( EOSTable % TS % minValues(1) )
  logrhomax = LOG10( EOSTable % TS % maxValues(1) )

  DO i=1, NumPoints
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


  WRITE (*, '(4A22)' ) "rho=", " T=", " Ye=", "Pressure=" 

  DO i = 1, NumPoints 
    WRITE (*, '(4ES12.5)' ) rho(i), T(i), Ye(i), Interpolant(i) 
  END DO

!  WRITE (*, '(3A25)' ) "  Interp P= ", "Dir Call P=", "Residual=" 
  DO i = 1, NumPoints 
    WRITE (*, 99 ) rho(i), T(i), Ye(i), Interpolant(i), DirectCall(i), ABS( Interpolant(i) - DirectCall(i) ) 
  END DO

  L1norm = 0
  DO i =1, NumPoints 
    L1norm = L1norm + ABS( Interpolant(i) - DirectCall(i) ) !/ DirectCall(i)  
  END DO

  WRITE (*, '(A,ES12.5)' ) "L1norm =" , L1norm

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInterpolationTest

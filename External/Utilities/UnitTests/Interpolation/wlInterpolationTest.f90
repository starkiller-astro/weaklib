PROGRAM wlInterpolationTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
                  

  implicit none

  INTEGER  :: i, j, k, TestUnit1, TestUnit2, TestUnit3, TestUnit4
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  INTEGER :: NumGoodPoints
  TYPE(EquationOfStateTableType) :: EOSTable

  REAL(dp), DIMENSION(:), ALLOCATABLE :: Interpolant
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: DirectCall 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Interpolants 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Derivatives 
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
  REAL(dp), DIMENSION(:), ALLOCATABLE :: thermenergy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: gammaone
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rand
  REAL(dp) :: Yemin, Yemax, logrho, logTmin, logTmax, logrhomin, logrhomax, epsilon
  INTEGER, PARAMETER :: NumPoints = 1000
  REAL(dp), DIMENSION(15) :: L1norm, Maxnorm, L1norm2, Maxnorm2
  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  LOGICAL            :: fail        ! did EoS fail to converge

  99 FORMAT ("rho=", es12.5, " T=", es12.5, " Ye=" , es12.5, " Int=", es12.5, &
             " DC=", es12.5, " Diff=", es12.5 ) 
 4018 FORMAT (i4,3(es12.5))
  epsilon = 1.d-100

  OPEN( newunit = TestUnit1, FILE="SingleInterpolateTest.d")
  OPEN( newunit = TestUnit2, FILE="InterpolateAllTest.d")
  OPEN( newunit = TestUnit3, FILE="SingleDerivativeTest.d")
  OPEN( newunit = TestUnit4, FILE="AllDerivativesTest.d")

  LScompress = '220'
  LSFilePath = '../../../LS/Data'
  EOSFlag = "L"

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  ALLOCATE( rho( NumPoints ), T( NumPoints ), Ye( NumPoints ), rand( NumPoints, 3 ),  &
            Interpolant( NumPoints ), DirectCall( NumPoints, 15), press( NumPoints ), &
            entrop( NumPoints ), energ( NumPoints ), chem_e( NumPoints ),             & 
            chem_p( NumPoints ), chem_n( NumPoints ), xn_prot( NumPoints ),           &
            xn_neut( NumPoints ), xn_alpha( NumPoints ), xn_heavy( NumPoints ),       &
            z_heavy( NumPoints ), a_heavy( NumPoints ), be_heavy( NumPoints ),        &
            thermenergy( NumPoints ), gammaone( NumPoints ),                          &
            Interpolants( NumPoints, EOSTable % DV % nVariables  ),                   & 
            Derivative( NumPoints, 3 ),                                               & 
            Derivatives( NumPoints, 3, EOSTable % DV % nVariables  ) )

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

  CALL wlExtInitializeEOS( LSFilePath, LScompress )

  DO i = 1, SIZE( rho )

    CALL wlGetFullEOS( rho(i), T(i), Ye(i), EOSFlag, fail, press(i), energ(i), &
                       entrop(i), chem_n(i), chem_p(i), chem_e(i), xn_neut(i), &
                       xn_prot(i), xn_alpha(i), xn_heavy(i), a_heavy(i),       &
                       z_heavy(i), be_heavy(i), thermenergy(i), gammaone(i), 1, 1, 1 ) 

    DirectCall(i,1) = press(i)
    DirectCall(i,2) = entrop(i)
    DirectCall(i,3) = energ(i)
    DirectCall(i,4) = chem_n(i)
    DirectCall(i,5) = chem_p(i)
    DirectCall(i,6) = chem_e(i)
    DirectCall(i,7) = xn_neut(i)
    DirectCall(i,8) = xn_prot(i)
    DirectCall(i,9) = xn_alpha(i)
    DirectCall(i,10) = xn_heavy(i)
    DirectCall(i,11) = a_heavy(i)
    DirectCall(i,12) = z_heavy(i)
    DirectCall(i,13) = be_heavy(i)
    DirectCall(i,14) = thermenergy(i)
    DirectCall(i,15) = gammaone(i)

  END DO
  
  NumGoodPoints = 0

  DO i = 1, SIZE(rho)
    DO j = 1, EOSTable % DV % nVariables
      CALL LogInterpolateSingleVariable &
             ( rho, T, Ye,                                   &
               EOSTable % TS % States(1) % Values,           &
               EOSTable % TS % States(2) % Values,           &
               EOSTable % TS % States(3) % Values,           &
               EOSTable % TS % LogInterp,                    &
               EOSTable % DV % Offsets(j),                   &
               EOSTable % DV % Variables(j) % Values(:,:,:), & 
               Interpolant )

      WRITE (TestUnit1,*) "Interpolant=", Interpolanti(i),   &
                          "Direct Call=", DirectCall(i,j),   & 
                            "L1Norm=", ( DirectCall(i,j) -   & 
                          Interpolants(i,j) )/ DirectCall(i,j)
  
    END DO
  END DO
  
  CALL LogInterpolateAllVariables( rho, T, Ye, EOSTable % TS % LogInterp, &
                               EOSTable % TS, EOSTable % DV, Interpolants ) 
  DO i = 1, SIZE(rho)
    DO j = 1, EOSTable % DV % nVariables 
      WRITE (TestUnit2,*) "Interpolant =", Interpolants(i,j), &
                          "Direct Call =", DirectCall(i,j),   &
                          "L1Norm=", ( DirectCall(i,j) -      & 
                          Interpolants(i,j) )/ DirectCall(i,j)
    END DO
  END DO

  DO i = 1, EOSTable % DV % nVariables
    CALL LogInterpolateDifferentiateSingleVariable( rho, T, Ye,        &
                         EOSTable % TS % States(1) % Values(:),        &
                         EOSTable % TS % States(2) % Values(:),        &
                         EOSTable % TS % States(3) % Values(:),        &
                         EOSTable % TS % LogInterp,                    &
                         EOSTable % DV % Offsets(i),                   &
                         EOSTable % DV % Variables(i) % Values(:,:,:), &
                         Interpolant(:), Derivative(:,:) )
      Interpolants(:,i) = Interpolant(:)
      Derivatives(:,:,i) = Derivative(:,:)
    END DO

    DO j = 1, EOSTable % DV % nVariables 
      WRITE (TestUnit3,*) "Derivatives =", Derivatives(:,:,j)
    END DO

    DO j = 1, EOSTable % DV % nVariables 
      WRITE (TestUnit3,*) "Derivatives =", Derivatives(:,:,j)
    END DO

    DO i = 1, SIZE(rho)
      DO j = 1, EOSTable % DV % nVariables 
        WRITE (TestUnit1,*) "Interpolant =", Interpolants(i,j), "Direct Call =", DirectCall(i,j)
      END DO
      WRITE (TestUnit3,4018) i, Interpolants(i,1), DirectCall(i,1), ( DirectCall(i,1) - Interpolants(i,1) )/ DirectCall(i,1)
    END DO

  CALL LogInterpolateDifferentiateAllVariables( rho, T, Ye, EOSTable % TS % LogInterp, EOSTable % TS, EOSTable % DV, Interpolants, Derivatives ) 

    DO j = 1, EOSTable % DV % nVariables 
      WRITE (TestUnit4,*) "Derivatives =", Derivatives(:,:,j)
    END DO

    DO j = 1, EOSTable % DV % nVariables 
      WRITE (TestUnit4,*) "Derivatives =", Derivatives(:,:,j)
    END DO

    DO i = 1, SIZE(rho)
      DO j = 1, EOSTable % DV % nVariables 
        WRITE (TestUnit1,*) "Interpolant =", Interpolants(i,j), "Direct Call =", DirectCall(i,j)
      END DO
      WRITE (TestUnit4,4018) i, Interpolants(i,1), DirectCall(i,1), ( DirectCall(i,1) - Interpolants(i,1) )/ DirectCall(i,1)
    END DO

  CLOSE(TestUnit1)
  CLOSE(TestUnit2)
  CLOSE(TestUnit3)
  CLOSE(TestUnit4)

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInterpolationTest

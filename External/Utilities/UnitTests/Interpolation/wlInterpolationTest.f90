PROGRAM wlInterpolationTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
                  

  implicit none

  INTEGER  :: i, j, nout, TestUnit1, TestUnit2, TestUnit3
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  INTEGER :: RhoLocation, TLocation, YeLocation
  INTEGER :: NumGoodPoints
  INTEGER, DIMENSION(1) :: LocMax
  LOGICAL, DIMENSION(3) :: LogInterp
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
  REAL(dp) :: Yemin, Yemax, logTmin, logTmax, logrhomin, logrhomax, epsilon
  INTEGER, PARAMETER :: NumPoints = 1000
  REAL(dp) :: LogT, logrho
  REAL(dp), DIMENSION(15) :: L1norm, Maxnorm, L1norm2, Maxnorm2
  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  LOGICAL            :: fail        ! did EoS fail to converge

  99 FORMAT ("rho=", es12.5, " T=", es12.5, " Ye=" , es12.5, " Int=", es12.5, &
             " DC=", es12.5, " Diff=", es12.5 ) 
  epsilon = 1.d-100
  nout  = 3216

  OPEN( newunit = TestUnit1, FILE="InterpolateAllTest.d")
  OPEN( newunit = TestUnit2, FILE="SingleDerivativeTest.d")
  OPEN( newunit = TestUnit3, FILE="AllDerivativesTest.d")

  LScompress = '220'
  LSFilePath = '../../../LS/Data'
  EOSFlag = "L"

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )

  OPEN(nout, FILE="OutputFile")

  WRITE (nout,*) "Table Minimums"
  DO i = 1, EOSTable % DV % nVariables
    WRITE (nout,*) EOSTable % DV % Names(i) , MINVAL( EOSTable % DV % Variables(i) % Values) 
  END DO

  LogInterp = (/.true.,.true.,.false./)

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

  DO j = 1, EOSTable % DV % nVariables

    CALL LogInterpolateSingleVariable &
           ( rho, T, Ye,                                   &
             EOSTable % TS % States(1) % Values,           &
             EOSTable % TS % States(2) % Values,           &
             EOSTable % TS % States(3) % Values,           &
             LogInterp,                                    &
             EOSTable % DV % Offsets(j),                   &
             EOSTable % DV % Variables(j) % Values(:,:,:), & 
             Interpolant )

    WRITE (nout,*) EOSTable % DV % Names(j), " Interpolation Comparison"

    WRITE (nout,*) "Interpolant =", Interpolant, "Direct Call =", DirectCall(:,j) 

    L1norm(j) = SUM( ABS( Interpolant - DirectCall(:,j) )  & 
                  / ( ABS( DirectCall(:,j) ) + epsilon ),  &
                  MASK = DirectCall(:,1).gt.0.0d0 )

    NumGoodPoints = COUNT( DirectCall(:,1).gt.0.0d0 )

    Maxnorm(j) = MAXVAL( ABS( Interpolant - DirectCall(:,j) ) &
                   / ( ABS( DirectCall(:,j) ) + epsilon ) ,   &
                   MASK = DirectCall(:,1).gt.0.0d0 )

    LocMax = MAXLOC( ABS( Interpolant - DirectCall(:,j) )     &
               / ( ABS( DirectCall(:,j) ) + epsilon ) ,       &
               MASK = DirectCall(:,1).gt.0.0d0 )

    L1norm2(j) = SUM( ABS( Interpolant - DirectCall(:,j) ),  & 
                   MASK = DirectCall(:,1).gt.0.0d0 )

    Maxnorm2(j) = MAXVAL( ABS( Interpolant - DirectCall(:,j) ),  &
                    MASK = DirectCall(:,1).gt.0.0d0 )
    
    L1norm(j) = MIN( L1norm(j), L1norm2(j) )/NumPoints
    Maxnorm(j) = MIN( Maxnorm(j), Maxnorm2(j) )


  END DO
  !WRITE ( *, * ) "press        entropy      energy       chem_n       chem_p       chem_e       xn_neut      xn_prot      xn_alpha      xn_heavy      a_heavy     z_heavy      be_heavy"
  !WRITE ( *, '( 13(es12.5,x) )' ) L1norm(1)/Zonelimit, L1norm(2)/Zonelimit,   &
                                  !L1norm(3)/Zonelimit, L1norm(4)/Zonelimit,   &
                                  !L1norm(5)/Zonelimit, L1norm(6)/Zonelimit,   &
                                  !L1norm(7)/Zonelimit, L1norm(8)/Zonelimit,   &
                                  !L1norm(9)/Zonelimit, L1norm(10)/Zonelimit,  &
                                  !L1norm(11)/Zonelimit, L1norm(12)/Zonelimit, &
                                  !L1norm(13)/Zonelimit

  !WRITE (*,*) NumGoodPoints

  CLOSE(nout)

  
  CALL LogInterpolateAllVariables( rho, T, Ye, LogInterp, EOSTable % TS, EOSTable % DV, Interpolants ) 
  DO i = 1, SIZE(rho)
    DO j = 1, EOSTable % DV % nVariables 
      WRITE (TestUnit1,*) "Interpolant =", Interpolants(i,j), "Direct Call =", DirectCall(i,j)
    END DO
  END DO

  DO i = 1, EOSTable % DV % nVariables
    CALL LogInterpolateDifferentiateSingleVariable( rho, T, Ye, &
                         EOSTable % TS % States(1) % Values(:), &
                         EOSTable % TS % States(2) % Values(:), &
                         EOSTable % TS % States(3) % Values(:), &
                         LogInterp,                             &
                         EOSTable % DV % Offsets(i), &
                         EOSTable % DV % Variables(i) % Values(:,:,:), &
                         Interpolant(:), Derivative(:,:) )
      Interpolants(:,i) = Interpolant(:)
      Derivatives(:,:,i) = Derivative(:,:)
    END DO

  DO j = 1, EOSTable % DV % nVariables 
    WRITE (TestUnit2,*) "Derivatives =", Derivatives(:,:,j)
  END DO


  CALL LogInterpolateDifferentiateAllVariables( rho, T, Ye, LogInterp, EOSTable % TS, EOSTable % DV, Interpolants, Derivatives ) 

  DO j = 1, EOSTable % DV % nVariables 
    WRITE (TestUnit3,*) "Derivatives =", Derivatives(:,:,j)
  END DO

  CLOSE(TestUnit1)
  CLOSE(TestUnit2)
  CLOSE(TestUnit3)

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlInterpolationTest

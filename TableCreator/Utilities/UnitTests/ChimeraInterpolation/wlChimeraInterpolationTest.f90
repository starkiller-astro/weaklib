PROGRAM wlChimeraInterpolationTest

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlExtPhysicalConstantsModule, ONLY: kfm
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF

  implicit none

  INTEGER  :: i, j, unitout, ErrorUnit, TestUnit, ProfileUnit, ZoneLimit
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  INTEGER :: RhoLocation, TLocation, YeLocation
  INTEGER :: NumGoodPoints
  INTEGER, DIMENSION(1) :: LocMax
  INTEGER, DIMENSION(3) :: LogInterp
  TYPE(EquationOfStateTableType) :: EOSTable

  REAL(dp), DIMENSION(:), ALLOCATABLE :: Interpolant
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: TableProfile 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ChimeraProfile
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Radius 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: p 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: s
  REAL(dp), DIMENSION(:), ALLOCATABLE :: e_internal
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_e
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_p
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_n 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xp 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xa 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xhe 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: z_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: a_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: be_heavy 
  REAL(dp) :: minrho, epsilon, factor
  INTEGER, PARAMETER :: NumPoints = 720
  REAL(dp) :: LogT, logrho
  REAL(dp), DIMENSION(13) :: L1norm, Maxnorm, L1norm2, Maxnorm2
  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath
  CHARACTER(len=12)  :: FileName

  LOGICAL            :: fail        ! did EoS fail to converge

  99 FORMAT ("rho=", es12.5, " T=", es12.5, " Ye=" , es12.5, " Int=", es12.5, &
             " DC=", es12.5, " Diff=", es12.5 ) 
  epsilon = 1.d-100

  LScompress = '220'
  LSFilePath = '../../../External/LS/Data'
  EOSFlag = "L"

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )
  
!  minrho = EOSTable % TS % minValues(1) 
  minrho = 1.0e11 

  OPEN( newunit = unitout, FILE="StandardResInterpolationData10ms.d")
  OPEN( newunit = ErrorUnit, FILE="NewStandardTableErrors7mspb.d")
  OPEN( newunit = TestUnit, FILE="StandardResTableMap.d")

!  WRITE (unitout,*) "Table Minimums"
!  DO i = 1, EOSTable % DV % nVariables
!    WRITE (unitout,*) EOSTable % DV % Names(i) , MINVAL( EOSTable % DV % Variables(i) % Values) 
!  END DO

  LogInterp = (/1,1,0/)

  ALLOCATE( Radius( NumPoints ), rho( NumPoints ), T( NumPoints ), Ye( NumPoints ), &
            Interpolant( NumPoints ), ChimeraProfile( NumPoints, 13), p( NumPoints ), &
            TableProfile( NumPoints, 13 ), &
            s( NumPoints ), e_internal( NumPoints ), chem_e( NumPoints ),          & 
            chem_p( NumPoints ), chem_n( NumPoints ), xp( NumPoints ),             &
            xn( NumPoints ), xa( NumPoints ), xhe( NumPoints ),                    &
            z_heavy( NumPoints ), a_heavy( NumPoints ), be_heavy( NumPoints ) ) 

  FileName = 'Output10ms.d'
  OPEN( newunit = ProfileUnit, file = FileName )

  DO i = 1, NumPoints
    READ(ProfileUnit,'(17(es11.3,x))') Radius(i), rho(i), T(i), Ye(i), p(i), s(i), e_internal(i), &
                                    chem_n(i), chem_p(i), chem_e(i), xn(i), xp(i),  &
                                    xa(i), xhe(i), z_heavy(i), a_heavy(i), be_heavy(i)

    IF ( rho(i) > minrho ) ZoneLimit = i
    ChimeraProfile(i,1) = p(i)
    ChimeraProfile(i,2) = s(i)
    ChimeraProfile(i,3) = e_internal(i)
    ChimeraProfile(i,4) = chem_e(i)
    ChimeraProfile(i,5) = chem_p(i)
    ChimeraProfile(i,6) = chem_n(i)
    ChimeraProfile(i,7) = xp(i)
    ChimeraProfile(i,8) = xn(i)
    ChimeraProfile(i,9) = xa(i)
    ChimeraProfile(i,10) = xhe(i)
    ChimeraProfile(i,11) = z_heavy(i)
    ChimeraProfile(i,12) = a_heavy(i)
    ChimeraProfile(i,13) = be_heavy(i)

  END DO
 
  WRITE (*,*) "P =", ChimeraProfile(1,1)

  DO i = 1, 13 !EOSTable % DV % nVariables

    CALL LogInterpolateSingleVariable( rho(1:ZoneLimit), T(1:ZoneLimit), &
                                       Ye(1:ZoneLimit),                              &
                                       EOSTable % TS % States(1) % Values,           &
                                       EOSTable % TS % States(2) % Values,           &
                                       EOSTable % TS % States(3) % Values,           &
                                       LogInterp,                                    &
                                       EOSTable % DV % Offsets(i),                   &
                                       EOSTable % DV % Variables(i) % Values(:,:,:), & 
                                       Interpolant(1:ZoneLimit) )

    TableProfile(1:ZoneLimit,i) = Interpolant(1:ZoneLimit)

    L1norm(i) = SUM(ABS( Interpolant(1:ZoneLimit) &
                       -  ChimeraProfile(1:ZoneLimit,i) ) / &
                      ( ABS( ChimeraProfile(1:ZoneLimit,i) ) + epsilon ) )

!    Maxnorm(i) = MAXVAL( ABS( Interpolant(1:ZoneLimit) - ChimeraProfile(1:ZoneLimit,i) )  &
!                   / ( ABS( ChimeraProfile(1:ZoneLimit,i) ) + epsilon ) ) 

!    LocMax = MAXLOC( ABS( Interpolant(1:ZoneLimit) - ChimeraProfile(1:ZoneLimit,i) ) &
!               / ( ABS( ChimeraProfile(1:ZoneLimit,i) ) + epsilon ) ,  &
!               MASK = ChimeraProfile(1:ZoneLimit,1).gt.0.0d0 )

  END DO

  !WRITE ( *, * ) "press        entrop       energy       chem_e       chem_p       chem_n       xn_prot      xn_neut      xn_alpha      xn_heavy      z_heavy     a_heavy      be_heavy"
  !WRITE ( *, '( 13(es12.5,x) )' ) L1norm(1)/Zonelimit, L1norm(2)/Zonelimit,   &
  !                                L1norm(3)/Zonelimit, L1norm(4)/Zonelimit,   &
  !                                L1norm(5)/Zonelimit, L1norm(6)/Zonelimit,   &
  !                                L1norm(7)/Zonelimit, L1norm(8)/Zonelimit,   &
  !                                L1norm(9)/Zonelimit, L1norm(10)/Zonelimit,  &
  !                                L1norm(11)/Zonelimit, L1norm(12)/Zonelimit, &
  !                                L1norm(13)/Zonelimit

  !WRITE ( *,'( 13(es12.5,x) )' ) Maxnorm(1), Maxnorm(2), Maxnorm(3), Maxnorm(4), &
  !                               Maxnorm(5), Maxnorm(6), Maxnorm(7), Maxnorm(8), &
  !                               Maxnorm(9), Maxnorm(10), Maxnorm(11),           &
  !                               Maxnorm(12), Maxnorm(13)
  
  !WRITE ( ErrorUnit, * ) "press        entrop       energy       chem_e       chem_p       chem_n       xn_prot      xn_neut      xn_alpha      xn_heavy      z_heavy     a_heavy      be_heavy"
  WRITE ( ErrorUnit, '( 13(es12.5,x) )' ) L1norm(1)/Zonelimit, L1norm(2)/Zonelimit,   &
                                  L1norm(3)/Zonelimit, L1norm(4)/Zonelimit,   &
                                  L1norm(5)/Zonelimit, L1norm(6)/Zonelimit,   &
                                  L1norm(7)/Zonelimit, L1norm(8)/Zonelimit,   &
                                  L1norm(9)/Zonelimit, L1norm(10)/Zonelimit,  &
                                  L1norm(11)/Zonelimit, L1norm(12)/Zonelimit, &
                                  L1norm(13)/Zonelimit

  !WRITE ( ErrorUnit,'( 13(es12.5,x) )' ) Maxnorm(1), Maxnorm(2), Maxnorm(3), Maxnorm(4), &
  !                               Maxnorm(5), Maxnorm(6), Maxnorm(7), Maxnorm(8), &
  !                               Maxnorm(9), Maxnorm(10), Maxnorm(11),           &
  !                               Maxnorm(12), Maxnorm(13)
  


    CALL LogInterpolateSingleVariable( rho(1:ZoneLimit), T(1:ZoneLimit), &
                                       Ye(1:ZoneLimit),                              &
                                       EOSTable % TS % States(1) % Values,           &
                                       EOSTable % TS % States(2) % Values,           &
                                       EOSTable % TS % States(3) % Values,           &
                                       LogInterp,                                    &
                                       EOSTable % DV % Offsets(i),                   &
                                       EOSTable % DV % Variables(13) % Values(:,:,:), &
                                       Interpolant(1:ZoneLimit) )
! DO j = 1, 13
  DO i = 1, ZoneLimit
    WRITE (unitout,'(28es11.4)') Radius(i), rho(i), TableProfile(i,1), ChimeraProfile(i,1), &
            TableProfile(i,2), ChimeraProfile(i,2), TableProfile(i,3), ChimeraProfile(i,3), &
            TableProfile(i,4), ChimeraProfile(i,4), TableProfile(i,5), ChimeraProfile(i,5), &
            TableProfile(i,6), ChimeraProfile(i,6), TableProfile(i,7), ChimeraProfile(i,7), &
            TableProfile(i,8), ChimeraProfile(i,8), TableProfile(i,9), ChimeraProfile(i,9), &
            TableProfile(i,10), ChimeraProfile(i,10), TableProfile(i,11), ChimeraProfile(i,11), &
            TableProfile(i,12), ChimeraProfile(i,12), TableProfile(i,13), ChimeraProfile(i,13)  
  END DO

  CLOSE(unitout)

 ! DO i = 1 , SIZE( EOSTable % TS % States(2) % Values ) 
 !   WRITE(TestUnit,'(i4,7es22.15)' ) i, EOSTable % TS % States(2) % Values(i),  &
 !                       10**(EOSTable % DV % Variables(3) % Values(139,i,1))
 ! END DO

 ! WRITE (*,*) "Internal Energy Monotonicity Check"
 !
 ! CALL MonotonicityCheck( EOSTable % DV % Variables(3) % Values(:,:,:), &
 !                         EOSTable % nPoints(1), EOSTable % nPoints(2), &
 !                         EOSTable % nPoints(3), 2, EOSTable % DV % Repaired(:,:,:) )

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlChimeraInterpolationTest

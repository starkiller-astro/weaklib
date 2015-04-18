PROGRAM wlChimeraProfileReader1Dtest

  USE wlKindModule, ONLY: dp 
  USE wlIOModuleCHIMERA
  implicit none

  INTEGER  :: i, Maxzone, nout
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: s 

  nout = 3216
  Maxzone = 720

  ALLOCATE( r( Maxzone ), rho( Maxzone ), T( Maxzone ), Ye( Maxzone ), s( Maxzone ) ) 

  CALL ReadChimeraProfile1D( "10ms.d", Maxzone, r, rho, T, Ye, s, 46)

  OPEN(nout, FILE="Output10ms10ms10ms10ms10ms10ms10ms10ms10ms10ms.d")
  DO i = Maxzone, 1, -1
    WRITE (nout,'(5(es12.5,x))') r(i), rho(i), T(i), Ye(i), s(i) 
  END DO
  CLOSE(nout)

END PROGRAM wlChimeraProfileReader1Dtest

PROGRAM wlChimeraProfileReader1Dtest

  USE wlKindModule, ONLY: dp 
  USE wlIOModuleCHIMERA
  implicit none

  INTEGER  :: i, Maxzone
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  !REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  !REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  !REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 

  Maxzone = 540

  ALLOCATE( rho( Maxzone ), T( Maxzone ), Ye( Maxzone ) ) 

  CALL ReadChimeraProfile1D( "chimerafile", Maxzone, rho, T, Ye, 1)

  DO i = Maxzone, 1, -1
    WRITE (*,*) rho(i), T(i), Ye(i) 
  END DO

END PROGRAM wlChimeraProfileReader1Dtest

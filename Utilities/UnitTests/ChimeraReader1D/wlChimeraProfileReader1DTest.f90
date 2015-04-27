PROGRAM wlChimeraProfileReader1Dtest

  USE wlKindModule, ONLY: dp 
  USE wlIOModuleCHIMERA
  implicit none

  INTEGER  :: i, Maxzone, nout
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: p 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: s 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: e_internal
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_n
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_p
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_e
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xp 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xhe 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xa 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: a_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: z_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: be_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: u 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rstmss 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: vsound 

  nout = 3216
  Maxzone = 720

  ALLOCATE( r( Maxzone ), rho( Maxzone ), T( Maxzone ), Ye( Maxzone ),     &
              p( Maxzone ), s( Maxzone ), e_internal( Maxzone ),           &
              xn( Maxzone ), xp( Maxzone ), xhe( Maxzone ), xa( Maxzone ), &
              chem_n( Maxzone ), chem_p( Maxzone ), chem_e( Maxzone ),     &
              a_heavy( Maxzone ), z_heavy( Maxzone ), be_heavy( Maxzone ), &
              u( Maxzone ), rstmss( Maxzone ), vsound( Maxzone ) ) 

  CALL ReadChimeraProfile1D( "10ms.d", Maxzone, r, rho, T, Ye, p, s, e_internal, &
                               xn, xp, xhe, xa, chem_n, chem_p, chem_e, a_heavy,  &
                               z_heavy, be_heavy, u, rstmss, vsound, 6932)

  OPEN(nout, FILE="Output10ms.d")

  DO i = 1, Maxzone
!    WRITE (nout,'(20(es12.5,x))') r(i), rho(i), T(i), Ye(i), p(i), s(i), e_internal(i), &
!                                    xn(i), xp(i), xhe(i), xa(i), chem_n(i), chem_p(i),  &
!                                    chem_e(i), a_heavy(i), z_heavy(i), be_heavy(i),     &
!                                    u(i), rstmss(i), vsound(i) 

    WRITE (nout,'(17(es11.3,x))') r(i), rho(i), T(i), Ye(i), p(i), s(i), e_internal(i),   &
                                    chem_n(i), chem_p(i), chem_e(i), xn(i), xp(i),  &
                                    xa(i), xhe(i), z_heavy(i), a_heavy(i), be_heavy(i)
                                     


  END DO
  CLOSE(nout)

END PROGRAM wlChimeraProfileReader1Dtest

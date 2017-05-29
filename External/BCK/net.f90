SUBROUTINE net
!-----------------------------------------------------------------------
!
!    File:         net
!    Module:       net
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/17/11
!
!    Purpose:
!      To compute the nuclear statistical equilibrium of neutrons,
!       protons, alpha nuclei and 55Fe.
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Subprograms called:
!        none
!
!    Include files:
!  wlKindModule
!  wlExtNumericalModule
!  wlExtPhysicalConstantsModule
!
!  eos_bck_module
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE wlExtNumericalModule, ONLY: zero, third, half, one, epsilon, pi
USE wlExtPhysicalConstantsModule, ONLY: mb, hbarc, rmu, cm3fm3, kmev

USE eos_bck_module, ONLY: d=>dbck, t=>tbck, ye=>yebck, un, uhat, xpbck,  &
& xnbck, xabck, xhbck, ahbck, zabck, ph, sh, eh, b_energy=>b, pd, ed, sd, &
& upack, theta, b_hvy, nnc_nse_bck, ye_nnc_nse_min

IMPLICIT none

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, SAVE                     :: first = .true.

INTEGER, PARAMETER                :: ncnvge = 20        ! number of iterations to converge

INTEGER                           :: n                  ! do index
INTEGER                           :: n_flunk            ! number of attempts at convergece
INTEGER                           :: itrat              ! iteration index

INTEGER                           :: it_uhat            ! iteration index for uhat iteration
INTEGER                           :: it_un              ! iteration index for un iteration
INTEGER, PARAMETER                :: it_uhat_max = 50   ! maximum iteration index for uhat iteration
INTEGER, PARAMETER                :: it_un_max = 50     ! maximum iteration index for un iteration

REAL(dp), PARAMETER           :: tol = 1.d-07       ! tolerance

REAL(dp), PARAMETER, DIMENSION(nnc_nse_bck)                    &
&                                 :: a = (/ 1.d0,                  &
&                                           1.d0,                  &
&                                           4.d0,                  &
&                                          12.d0,                  &
&                                          16.d0,                  &
&                                          20.d0,                  &
&                                          24.d0,                  &
&                                          28.d0,                  &
&                                          32.d0,                  &
&                                          36.d0,                  &
&                                          40.d0,                  &
&                                          44.d0,                  &
&                                          48.d0,                  &
&                                          52.d0,                  &
&                                          56.d0,                  &
&                                          60.d0,                  &
&                                          56.d0 /)     ! nuclear mass number

REAL(dp), PARAMETER, DIMENSION(nnc_nse_bck)                    &
&                                 :: z = (/ 0.d0,                  &
&                                           1.d0,                  &
&                                           2.d0,                  &
&                                           6.d0,                  &
&                                           8.d0,                  &
&                                          10.d0,                  &
&                                          12.d0,                  &
&                                          14.d0,                  &
&                                          16.d0,                  &
&                                          18.d0,                  &
&                                          20.d0,                  &
&                                          22.d0,                  &
&                                          24.d0,                  &
&                                          26.d0,                  &
&                                          28.d0,                  &
&                                          30.d0,                  &
&                                          26.d0 /)     ! nuclear charge number

REAL(dp), PARAMETER, DIMENSION(nnc_nse_bck)                    &
&                                 :: z_over_a(1:nnc_nse_bck) = z(1:nnc_nse_bck)/a(1:nnc_nse_bck)
                                                        ! charge to mass ratios

REAL(dp), PARAMETER, DIMENSION(nnc_nse_bck)                    &
&                                 :: g = (/ 2.d0,                  &
&                                           2.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0,                  &
&                                           1.d0 /)     ! nuclear charge number
                                                       
REAL(dp), PARAMETER, DIMENSION(nnc_nse_bck)                    &
&                                 :: z2_over_a(1:nnc_nse_bck)  = z(1:nnc_nse_bck)*z(1:nnc_nse_bck)/a(1:nnc_nse_bck)
                                                        ! charge squared to mass ratios
                                                        
REAL(dp), DIMENSION(nnc_nse_bck), SAVE :: a52(1:nnc_nse_bck)    ! a**5/2
REAL(dp), DIMENSION(nnc_nse_bck), SAVE :: a52_1(1:nnc_nse_bck)  ! 1/a**5/2

REAL(dp), PARAMETER, DIMENSION(nnc_nse_bck)                    &
&                                 :: b = (/ 0.d0      ,            &
&                                           0.d0      ,            &
&                                          -7.073915d0,            &
&                                          -7.680144d0,            &
&                                          -7.976206d0,            &
&                                          -8.032240d0,            &
&                                          -8.260709d0,            &
&                                          -8.447744d0,            &
&                                          -8.493134d0,            &
&                                          -8.519909d0,            &
&                                          -8.551301d0,            &
&                                          -8.533518d0,            &
&                                          -8.572210d0,            &
&                                          -8.609598d0,            &
&                                          -8.642709d0,            &
&                                          -8.583273d0,            &
&                                          -8.790323d0 /)

REAL(dp)                      :: tolone             ! difference or x_tot from unity
REAL(dp)                      :: tolza              ! difference or z_over_a_sum from ye
REAL(dp)                      :: eps                ! measure of convergence

REAL(dp), SAVE                :: c0                 ! Saha constant
REAL(dp)                      :: therm              ! Saha constant
REAL(dp)                      :: therm_1            ! 1/therm
REAL(dp), PARAMETER           :: three9 = 0.999d0
REAL(dp), PARAMETER           :: five9 = 0.99999d0

REAL(dp)                      :: x_n                ! neutron mass fraction
REAL(dp)                      :: x_p                ! proton mass fraction
REAL(dp)                      :: x_a                ! alpha mass fraction
REAL(dp)                      :: x_Si               ! 28Si mass fraction
REAL(dp)                      :: x_Fe               ! 56Fe mass fraction
REAL(dp)                      :: x_Ni               ! 56Ni mass fraction
REAL(dp)                      :: up                 ! proton chemical potential
REAL(dp)                      :: u_a                ! alpha chemical potential
REAL(dp)                      :: u_Fe               ! 56Fe chemical potential
REAL(dp)                      :: u_Ni               ! 56Ni chemical potential

REAL(dp)                      :: x_tot              ! sum of the xn_net's
REAL(dp)                      :: z_over_a_sum       ! sum of the xn_net * z/a
REAL(dp)                      :: z_sum              ! sum of the xn_net * z
REAL(dp)                      :: a_sum              ! sum of the xn_net * a
REAL(dp)                      :: zza_sum            ! sum of the xn_net**2/a
REAL(dp)                      :: b_sum              ! sum of the xn_net * b
REAL(dp)                      :: a_1_sum            ! sum of the xn_net/a
REAL(dp)                      :: xhsum              ! mass fraction of heavy nuclei
REAL(dp)                      :: excited            ! energy of excited states
REAL(dp)                      :: rdenom             ! denominaor
REAL(dp)                      :: dun                ! correction to un
REAL(dp)                      :: duhat              ! correction to uhat

REAL(dp)                      :: uhat_min           ! minimum value of uhat to begin bisection iteration
REAL(dp)                      :: uhat_max           ! maximum value of uhat to begin bisection iteration
REAL(dp)                      :: un_min             ! minimum value of un to begin bisection iteration
REAL(dp)                      :: un_max             ! maximum value of un to begin bisection iteration

REAL(dp), DIMENSION(nnc_nse_bck) :: xn_net          ! composition mass fraction array

REAL(dp), EXTERNAL            :: fexp

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' un will not converge in net; un_min=',es11.3,' un_max=',es11.3, &
& ' un=',es11.3,' x_tot=',es11.3,' d=',es11.3,' t=',es11.3,' ye=',es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize nuclear parameters
!-----------------------------------------------------------------------

IF ( first ) THEN

  first                = .false.
  c0                   = ( mb/( 2.d0 * pi * hbarc**2 ) )**( 3.d0/2.d0 )
  a52(1:nnc_nse_bck)   = a(1:nnc_nse_bck)**2.5
  a52_1(1:nnc_nse_bck) = 1.d0/a52(1:nnc_nse_bck)

END IF ! first

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

n_flunk               = 0
therm                 = c0 * t * sqrt(t) / d
therm_1               = 1.d0/therm

!-----------------------------------------------------------------------
!  Begin cycle
!-----------------------------------------------------------------------

  100 CONTINUE

!-----------------------------------------------------------------------
!  Estimate initial values for un and uhat
!-----------------------------------------------------------------------

IF ( n_flunk == 0 ) THEN

!-----------------------------------------------------------------------
!  If n_flunk <= 1 and ye < ye_nnc_nse_min, compute initial un and uhat
!   assuming only Fe56 and neutrons if ye < ye_nnc_nse_min
!-----------------------------------------------------------------------

  IF ( ye < ye_nnc_nse_min ) THEN

    x_Fe              = ye /ye_nnc_nse_min
    x_n               = 1.d0 - x_n
    un                = t * DLOG( half * therm_1 * x_n )
    uhat              = ( un - b(17) - t/a(17) * DLOG( therm_1 * a52_1(17) * x_Fe ) ) &
&                     / ye_nnc_nse_min
    GO TO 50

!-----------------------------------------------------------------------
!  If n_flunk <= 1 and ye_nnc_nse_min < ye < 0.5d0, compute initial un
!   and uhat assuming only Fe56 and Ni56
!-----------------------------------------------------------------------

  ELSE IF ( ye_nnc_nse_min < ye  .and.  ye <= 0.5d0 ) THEN

    x_Fe              = DMAX1( ( 1.d0 - 2 * ye )/( 1.d0 - 2.d0 * ye_nnc_nse_min ), epsilon )
    x_Ni              = DMAX1( 1.d0 - x_Fe, epsilon )
    u_Fe              = t * DLOG( x_Fe * a52_1(17) * therm_1 )
    u_Ni              = t * DLOG( x_Ni * a52_1(15) * therm_1 )
    up                = 0.25d0 * ( 60.d0 * b(15) - 56.d0 * b(17) + 60.d0 * u_Ni/56.d0 - u_Fe )
    un                = - up + 2.d0 * b(15) + u_Ni/28.d0
    uhat              = un - up
    GO TO 50

!-----------------------------------------------------------------------
!  If n_flunk <= 1 and ye > 0.5, compute initial un and uhat assuming
!   only Ni56 and protons
!-----------------------------------------------------------------------

  ELSE ! ye > 0.5d0

    x_Ni             = DMIN1( 2.d0 * ( 1.d0 - ye), three9 )
    x_p              = 1.d0 - x_Ni
    up               = t * DLOG( 0.5d0 * therm_1 * x_p )
    uhat             = 2.d0 * ( b(15) + t/a(15) * DLOG( therm_1 * a52_1(15) * x_Ni ) - up )
    un               = uhat + up
    GO TO 50

  END IF ! ye < ye_nnc_nse_min

!-----------------------------------------------------------------------
!  If n_flunk = 1 and ye_nnc_nse_min < ye < 0.5, compute initial un and
!   uhat assuming only Fe56 and alphas
!-----------------------------------------------------------------------

ELSE IF ( n_flunk == 1 ) THEN

  IF ( ye <= 0.5d0 ) THEN

    x_Fe              = DMAX1( ( 1.d0 - 2 * ye )/( 1.d0 - 2.d0 * ye_nnc_nse_min ), epsilon )
    x_a               = DMAX1( 1.d0 - x_Fe, epsilon )
    u_Fe              = t * DLOG( x_Fe * a52_1(17) * therm_1 )
    u_a               = t * DLOG( x_a  * a52_1( 3) * therm_1 )
    un                = 0.25d0 * ( 56.d0 * b(17) - 52.d0 * b(3) + u_Fe - 13.d0 * u_a )
    up                = - un + 2.d0 * b(3) + 0.5d0 * u_a
    uhat              = un - up
    GO TO 50

!-----------------------------------------------------------------------
!  If n_flunk = 1 and ye > 0.5, compute initial un and uhat assuming
!   only protons and alphas
!-----------------------------------------------------------------------

  ELSE ! ye > 0.5d0

    x_a              = DMAX1( 0.5d0 * a(3) * ( 1.d0 - ye ), epsilon )
    x_p              = DMAX1( 1.d0 - x_a, epsilon )
    up               = t * DLOG( 0.5d0 * therm_1 * x_p )
    un               = 2.d0 * b(3) + 0.5d0 * t * DLOG( therm_1 * a52_1(3) * x_a ) - up
    uhat             = un - up
    GO TO 50

  END IF ! ye <= 0.5d0

!-----------------------------------------------------------------------
!  If n_flunk = 2 and ye < 0.5, compute initial un and uhat assuming
!   only neutrons and alphas
!-----------------------------------------------------------------------

ELSE IF ( n_flunk == 2   ) THEN

  IF ( ye < 0.5d0 ) THEN

    x_a              = 0.5d0 * a(3) * ye
    x_n              = DMAX1( 1.d0 - x_a, epsilon )
    un               = t * DLOG( 0.5d0 * therm_1 * x_n )
    uhat             = 2.d0 * un - 2.d0 * b(3) - 0.5d0 * t * DLOG( therm_1 * a52_1(3) * x_a )
    GO TO 50

!-----------------------------------------------------------------------
!  If n_flunk = 2 and ye > 0.5, compute initial un and uhat assuming
!   only protons and alphas
!-----------------------------------------------------------------------

  ELSE ! ye > 0.5d0

    x_a              = DMAX1( 0.5d0 * a(3) * ( 1.d0 - ye ), epsilon )
    x_p              = DMAX1( 1.d0 - x_a, epsilon )
    up               = t * DLOG( 0.5d0 * therm_1 * x_p )
    un               = 2.d0 * b(3) + 0.5d0 * t * DLOG( therm_1 * a52_1(3) * x_a ) - up
    uhat             = un - up
    GO TO 50

  END IF ! ye < 0.5d0

!-----------------------------------------------------------------------
!  If n_flunk = 3, compute initial un and uhat assuming only neutrons
!   and protons
!-----------------------------------------------------------------------

ELSE IF ( n_flunk == 3 ) THEN

  x_n                = 1.d0 - ye
  x_p                = ye
  uhat               = t * DLOG( x_n/x_p )
  un                 = t * DLOG( 0.5d0 * therm_1 * x_n )
  GO TO 50

!-----------------------------------------------------------------------
!  If n_flunk = 4, compute initial un and uhat assuming only 28Si
!   and protons
!-----------------------------------------------------------------------

ELSE IF ( n_flunk == 4 ) THEN

  x_Si               = DMIN1( 2.d0 * ( 1.d0 - ye ),three9 )
  x_p                = DMAX1( 1.d0 - x_Si, epsilon )
  up                 = t * DLOG( 0.5d0 * therm_1 * x_p )
  uhat               = 2.d0 * ( b(8) + t/a(8) * DLOG( therm_1 * a52_1(8) * x_Si ) - up )
  un                 = uhat + up
  GO TO 50

!-----------------------------------------------------------------------
!  If n_flunk = 5, compute initial un and uhat assuming only 56Fe
!   and protons
!-----------------------------------------------------------------------

ELSE IF ( n_flunk == 5 ) THEN

  x_Fe               = DMIN1( 2.d0 * ye, five9 )
  x_p                = DMAX1( 1.d0 - x_Fe, epsilon )
  up                 = t * DLOG( 0.5d0 * therm_1 * x_p )
  uhat               = ye_nnc_nse_min * ( b(17) + t/a(17) * DLOG( therm_1 * a52_1(17) * x_Fe ) - up )
  un                 = uhat + up
  GO TO 50

!-----------------------------------------------------------------------
!  If n_flunk = 6, in desparation, do a bisection iteration
!-----------------------------------------------------------------------

ELSE IF ( n_flunk == 6 ) THEN

  uhat_min           = -100.d0
  uhat_max           = 100.d0
  
  DO it_uhat = 1, it_uhat_max
    uhat             = 0.5d0 * ( uhat_min + uhat_max )
    un_min           = - 100.d0
    un_max           = 100.d0
    DO it_un = 1, it_un_max
      un             = 0.5d0 * ( un_min + un_max )
      DO n = 1, nnc_nse_bck
        xn_net(n)    = fexp( DLOG( g(n) * a52(n) * therm ) + ( a(n) * un - z(n) * uhat - a(n) * b(n) )/t )
      END DO ! n = 1, nnc_nse_bck
      x_tot          = SUM( xn_net(1:nnc_nse_bck) )
      tolone         = DLOG(x_tot)
      eps            = DABS(tolone)
      IF ( eps < tol ) EXIT
      IF ( it_un == it_un_max ) THEN
        WRITE (*,1001) un_min, un_max, un, x_tot, rmu * d/cm3fm3, t/kmev, ye
        STOP
      END IF ! it_un == it_un_max
      IF ( x_tot < 1.d0 ) THEN
        un_min       = un
      ELSE
        un_max       = un
      END IF ! x_tot < 1.d0
    END DO ! it_un = 1, it_un_max
    z_over_a_sum     = SUM( xn_net(1:nnc_nse_bck) * z_over_a(1:nnc_nse_bck) )
    tolza            = DLOG(z_over_a_sum/ye)
    eps              = DABS(tolza)
    IF ( eps < tol ) EXIT
    IF ( it_uhat == it_uhat_max ) THEN
      WRITE (*,1001) un_min, un_max, un, x_tot, rmu * d/cm3fm3, t/kmev, ye
      STOP
    END IF ! it_uhat == it_uhat_max
    IF ( z_over_a_sum < ye ) THEN
      uhat_max       = uhat
    ELSE
      uhat_min       = uhat
    END IF ! z_over_a_sum < ye
  END DO ! it_uhat = 1, it_uhat_max

  excited            = zero

  DO n = 1, nnc_nse_bck
    xn_net(n)        = fexp( DLOG( g(n) * a52(n) * therm ) + ( a(n) * un - z(n) * uhat - a(n) * b(n) )/t )
  END DO ! n = 1, nnc_nse_bck
  x_tot              = SUM( xn_net(1:nnc_nse_bck) )
  zza_sum            = SUM( xn_net(1:nnc_nse_bck) * z2_over_a(1:nnc_nse_bck) )
  b_sum              = SUM( xn_net(1:nnc_nse_bck) * b(1:nnc_nse_bck) )
  a_1_sum            = SUM( xn_net(1:nnc_nse_bck)/a(1:nnc_nse_bck) )
  GO TO 200

END IF ! n_flunk <= 1

!-----------------------------------------------------------------------
!  Iterate on un and uhat
!-----------------------------------------------------------------------

   50 CONTINUE

DO itrat = 1, ncnvge

  excited            = zero
  
  DO n = 1, nnc_nse_bck
    xn_net(n)        = fexp( DLOG( g(n) * a52(n) * therm ) + ( a(n) * un - z(n) * uhat - a(n) * b(n) )/t )
  END DO ! n = 1, nnc_nse_bck
  x_tot              = SUM( xn_net(1:nnc_nse_bck) )
  zza_sum            = SUM( xn_net(1:nnc_nse_bck) * z2_over_a(1:nnc_nse_bck) )
  z_sum              = SUM( xn_net(1:nnc_nse_bck) * z(1:nnc_nse_bck) )
  a_sum              = SUM( xn_net(1:nnc_nse_bck) * a(1:nnc_nse_bck) )
  b_sum              = SUM( xn_net(1:nnc_nse_bck) * b(1:nnc_nse_bck) )
  a_1_sum            = SUM( xn_net(1:nnc_nse_bck)/a(1:nnc_nse_bck) )
  z_over_a_sum       = SUM( xn_net(1:nnc_nse_bck) * z_over_a(1:nnc_nse_bck) )
  tolone             = DLOG(x_tot)
  tolza              = DLOG(z_over_a_sum/ye)
  eps                = DABS(tolone) + DABS(tolza)

!-----------------------------------------------------------------------
!  If eps < tol, done
!-----------------------------------------------------------------------

  IF ( eps < tol ) go to 200

!-----------------------------------------------------------------------
!  tolone >= 75.  .and.  tolza >= 75., try again with a different guess
!-----------------------------------------------------------------------

  IF ( tolone >= 75.  .and.  tolza >= 75. ) THEN
    n_flunk          = n_flunk + 1
    GO TO 100
  END IF ! tolone >= 75.  .and.  tolza >= 75.

!-----------------------------------------------------------------------
!  Compute the denominator of the solution matrix
!-----------------------------------------------------------------------

  rdenom             = ( ( z_sum/x_tot * z_sum/z_over_a_sum - a_sum/x_tot * zza_sum/z_over_a_sum ) )

!-----------------------------------------------------------------------
!  If rdenom = 0, try again with a different guess
!-----------------------------------------------------------------------

  IF ( rdenom == zero) THEN
    n_flunk          = n_flunk + 1
    GO TO 100
  END IF ! rdenom == zero

!-----------------------------------------------------------------------
!  Compute corrections
!-----------------------------------------------------------------------

  rdenom             = t/rdenom
  dun                = ( tolone * zza_sum/z_over_a_sum - tolza/x_tot * z_sum) * rdenom
  duhat              = ( z_sum/z_over_a_sum * tolone - a_sum/x_tot * tolza) * rdenom
  un                 = un + dun
  uhat               = uhat + duhat
  
  IF ( itrat == ncnvge ) THEN
    n_flunk          = n_flunk + 1
    GO TO 100
  END IF ! itrat == ncnvge

END DO ! itrat = 1, ncnvge

  200 CONTINUE

!-----------------------------------------------------------------------
!  Store mass fractions
!-----------------------------------------------------------------------

xnbck                = xn_net(1)
xpbck                = xn_net(2)
xabck                = xn_net(3)
xhbck                = SUM( xn_net(4:nnc_nse_bck) )

!-----------------------------------------------------------------------
!  Compute thermodynamic quantities
!-----------------------------------------------------------------------

upack                = zero
theta                = zero

xhsum                = xhbck
ahbck                = SUM( xn_net(4:nnc_nse_bck) * a(4:nnc_nse_bck) )
ahbck                = ahbck/xhsum
zabck                = SUM( xn_net(4:nnc_nse_bck) * z_over_a(4:nnc_nse_bck) )

eh                   = zero
ph                   = zero
sh                   = zero

b_energy             = b_sum
b_hvy                = SUM( xn_net(4:nnc_nse_bck) * b(4:nnc_nse_bck) ) / ( SUM(xn_net(4:nnc_nse_bck) ) + epsilon )
eh                   = b_sum  + eh
pd                   = d * t * a_1_sum
ed                   = 1.5 * t * a_1_sum
sd                   = 2.5 * a_1_sum - ( un - ye * uhat - b_sum + xhsum * excited )/t

RETURN

END SUBROUTINE net

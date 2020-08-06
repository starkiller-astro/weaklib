!calculate the zeroth order annihilation kernel from Hannestad & Raffelt 1998
!using:
!neutrino energies, rho, T
!The decomposition can later be taken into account by
!passing rho*xn, rho*xp, and rho*SQRT( ABS( xn * xp ))
!in separate calls and adding the contributions to the
!overall kernel
!calculates s_a, which is 1/2 \Phi^a_0_Brem, in subroutine bremcal_weaklib
!calls two subroutines: s_brem and g_brem
!https://iopscience.iop.org/article/10.1086/306303
!https://wwwth.mpp.mpg.de/members/raffelt/mypapers/199804.pdf

MODULE HR98_Bremsstrahlung

USE wlKindModule, ONLY: dp
USE wlExtPhysicalConstantsModule, ONLY: kMeV, cvel
USE wlExtNumericalModule, ONLY: pi, half, epsilon

IMPLICIT NONE


CONTAINS

subroutine bremcal_weaklib( nez, egrid, rho_x, t, s_a )

INTEGER,                  INTENT(in) :: nez       ! number of neutrino energies
REAL(dp), DIMENSION(nez), INTENT(in) :: egrid     ! neutrino energy grid

REAL(dp), INTENT(in)       :: rho_x                 ! density [g cm^{-3}]
REAL(dp), INTENT(in)       :: t                   ! temperature [K]
!REAL(dp), INTENT(in)       :: xn                  ! free neutron mass fraction
!REAL(dp), INTENT(in)       :: xp                  ! free proton mass fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), DIMENSION(nez,nez), INTENT(out) :: s_a  ! differential annihilation kernel, 1/2 \Phi^a_0

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER               :: i                    ! index to loop over nn, pp, and np contributions

REAL(dp), PARAMETER   :: tthird   = 2.d0/3.d0

REAL(dp), PARAMETER   :: C_AG_FnB = 3.79d+4/cvel             ! Raffelt constant
REAL(dp), PARAMETER   :: conv1    = 1.d0 / ( 2.d0 * pi )**3  ! 1/(hc)**3 in natural units
REAL(dp), PARAMETER   :: conv2    = 3.d0 * 4.d0 * pi               ! integrates ( 3 - costh ) over solid angle
REAL(dp), PARAMETER   :: coef     = conv1 * conv2 * C_AG_FnB
REAL(dp), PARAMETER   :: coef_np  = 28.d0/3.d0
REAL(dp)              :: rho_14          ! rho/1e+14
REAL(dp)              :: t_10            ! tmev/10
REAL(dp), DIMENSION(nez,nez)   :: x      ! (enu + enubar [MeV])/tmev
REAL(dp)              :: tmev            ! kt [MeV]
REAL(dp)              :: eta_star        ! degeneracy parameter
REAL(dp)              :: gamma           ! spin fluctuation rate
REAL(dp)              :: y               ! pion mass parameter
REAL(dp), DIMENSION(nez,nez)   :: sb     ! dimensionless fitting parameter
REAL(dp)              :: gb              ! dimensionless fitting parameter
!REAL(dp)              :: rho_x           ! density adjusted for the composition [g cm^{-3}]

INTEGER               :: k, kb           ! neutrino(anti) loop variables

!REAL(dp), DIMENSION(nez,nez) :: s_aa   ! differential absorption kernel

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Loop over nn, pp, and np contributions
!-----------------------------------------------------------------------

!DO i = 1, 3


!-----------------------------------------------------------------------
!  Adjust density to reflect composition
!-----------------------------------------------------------------------

!  IF ( i == 1 ) THEN
!    rho_x          = rho * xn
!  ELSE IF ( i == 2 ) THEN
!    rho_x          = rho * xp
!  ELSE IF ( i == 3 ) THEN
!    rho_x          = rho * SQRT( ABS( xn * xp ) + epsilon )
!  END IF ! i == 1

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

  rho_14           = rho_x/1.d+14
  tmev             = kmev * t
  t_10             = tmev/10.d0

!-----------------------------------------------------------------------
!  Compute eta_star, the neutron effective degeneracy parameter
!   (Equation 36 Hannestad & Raffelt)
!-----------------------------------------------------------------------

  eta_star         = 3.04d0 * rho_14**tthird / t_10

!-----------------------------------------------------------------------
!  Compute gamma, the spin fluctuation rate
!   (Equation 37 Hannestad & Raffelt)
!-----------------------------------------------------------------------

  gamma            = 8.6d0 * rho_14 / sqrt(t_10)

!-----------------------------------------------------------------------
!  Compute y, the pion mass parameter
!   (Equation 38 Hannestad & Raffelt)
!-----------------------------------------------------------------------

  y                = 1.94d0 / t_10

!-----------------------------------------------------------------------
!  Compute x, the dimensionless neutrino energy sum
!   (Equation 35 Hannestad & Raffelt)
!-----------------------------------------------------------------------

  DO kb = 1,nez
    DO k = 1,nez
      x(k,kb)      = ( egrid(k) + egrid(kb) ) / tmev
    END DO ! k = 1,nez
  END DO ! kb = 1,nez

!-----------------------------------------------------------------------
!  Compute the dimensionless fitting parameter, sb
!-----------------------------------------------------------------------

  CALL s_brem( nez, x, y, eta_star, sb)

!-----------------------------------------------------------------------
!  Compute the dimensionless fitting parameter, gb
!-----------------------------------------------------------------------

  CALL g_brem( y, eta_star, gb )

!-----------------------------------------------------------------------
!  Compute the differential absorption kernel, s_a
!-----------------------------------------------------------------------

  DO kb = 1,nez
    DO k = 1,nez
      s_a(k,kb)     = gamma/( x(k,kb)**2 + ( half * gamma * gb )**2 )         &
&                    * sb(k,kb)/tmev * coef * rho_14
    END DO ! k = 1,nez
  END DO ! kb = 1,nez

!END DO ! i = 1, 3

!-----------------------------------------------------------------------
!  Sum the contributions
!-----------------------------------------------------------------------

!DO kb = 1,nez
!  DO k = 1,nez
!    s_a(k,kb)   = s_aa(k,kb,1) + s_aa(k,kb,2) + coef_np * s_aa(k,kb,3)
!  END DO ! k = 1,nez
!END DO ! kb = 1,nez

RETURN

end subroutine bremcal_weaklib

SUBROUTINE s_brem( nez, x_p, y_p, eta_star_p, s_fit )
!-----------------------------------------------------------------------
!
!    File:         s_brem
!    Module:       s_brem
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/23/02
!
!    Purpose:
!      To compute s-component of the kernel for neutrino bremss-
!       trahlung and inelastic neutrino scattering in a nucleon
!       field.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  y           : m_{pi}^{2}/M_{N}c^{2}kT
!  eta_star    : neutron degeneracy parameter
!
!    Output arguments:
!  s_fit       : dimensionless quantity related to S_{sigma}
!
!    Modules:
!  kind_module
!  numerical_module
!  array_module
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE wlExtNumericalModule, ONLY: zero, half, one, third, pi, pi2, epsilon

!USE math_functions_module, ONLY: fexp

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------
INTEGER,                      INTENT(in)  :: nez         ! number of neutrino energy bins
REAL(dp), DIMENSION(nez,nez), INTENT(in)  :: x_p         ! w/tmev
REAL(dp),                     INTENT(in)  :: y_p         ! pion mass parameter
REAL(dp),                     INTENT(in)  :: eta_star_p  ! degeneracy parameter

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), DIMENSION(nez,nez), INTENT(out) :: s_fit       ! dimensionless fitting parameter

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp)              :: x                               ! use in place of x_p
REAL(dp)              :: y                               ! use in place of y_p
REAL(dp)              :: eta_star                        ! use in place of eta_star_p

REAL(dp), PARAMETER   :: x_min = 1.d-10
REAL(dp), PARAMETER   :: y_min = 1.d-10
REAL(dp), PARAMETER   :: eta_min = 1.d-10
REAL(dp), PARAMETER   :: fi_third = 5.d0/3.d0
REAL(dp), PARAMETER   :: fi_sixth = 5.d0/6.d0

REAL(dp)              :: pi1_2                           ! pi**(1/2)
REAL(dp)              :: pi1_8                           ! pi**(1/8)
REAL(dp)              :: s_ND_num                        ! numerator of s_ND
REAL(dp)              :: s_ND_denom                      ! denominator of s_ND
REAL(dp)              :: s_ND                            ! nondegenerate limit of s_fit
REAL(dp)              :: u                               ! y/(2*eta_star)
REAL(dp)              :: u2                              ! u**2
REAL(dp)              :: u_arg                           ! a function of u
REAL(dp)              :: f_u                             ! a function of u
REAL(dp)              :: s_D                             ! degenerate limit of s_fit
REAL(dp)              :: C_fit                           ! fitting function
REAL(dp)              :: F_denom                         ! denominator of F_fit
REAL(dp)              :: F_fit                           ! fitting function
REAL(dp)              :: G_fit                           ! fitting function
REAL(dp)              :: h_fit                           ! fitting function
REAL(dp)              :: p_fit                           ! fitting function

INTEGER               :: k, kb                           ! neutrino(anti) loop variables

REAL(KIND=dp)         :: fexp                            ! exponential

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

pi1_2              = pi**(1.d0/2.d0)
pi1_8              = pi**(1.d0/8.d0)

!-----------------------------------------------------------------------
!  Replace input variables for modification
!-----------------------------------------------------------------------

y                  = y_p
eta_star           = eta_star_p

!-----------------------------------------------------------------------
!  Prevent singular behavior of terms
!-----------------------------------------------------------------------

y                  = MAX( y       , y_min   )
eta_star           = MAX( eta_star, eta_min )

DO kb = 1,nez
  DO k = 1,nez
    x              = x_p(k,kb)
    IF ( x >= zero ) THEN
      x            = MAX( x       , x_min   )
     ELSE
      x            = MIN( x       ,-x_min   )
    END IF

!-----------------------------------------------------------------------
!  Compute S_ND, nondegenerate approximation
!   (Equation 45 Hannestad & Raffelt)
!-----------------------------------------------------------------------

    s_ND_num       = 2.d0 * pi1_2 * ( x + 2.d0 - exp(- y/12.d0 ) )**1.5d0 &
&                  * ( x**2 + 2.d0 * x * y + fi_third * y**2 + 1.d0 )
    s_ND_denom     = pi1_2 + ( pi1_8 + x + y )**4
    s_ND           = s_ND_num/s_ND_denom
    s_ND           = MAX( s_ND, epsilon )

IF ( x < zero ) S_ND = S_ND * fexp(-x)

!-----------------------------------------------------------------------
!  Compute S_D, degenerate approximation
!   (Equation 46 & 47 Hannestad & Raffelt)
!-----------------------------------------------------------------------

    u              = SQRT( y/( 2.d0 * eta_star ) ) + 1.d-10
    u2             = u * u
    u_arg          = u2/( 2.d0 * SQRT( 2.d0 * u2 + 4.d0 ) )
    f_u            = one - fi_sixth * u * ATAN( 2.d0/u )                 &
&                  + u2/( 3.d0 * ( u2 + 4.d0 ) )                          &
&                  + ATAN( one/u_arg ) * third * u_arg
    s_D            = 3.d0 * ( half * pi )**2.5d0 * eta_star**( -2.5d0 )   &
&                  * ( x**2 + 4.d0 * pi2 ) * x * f_u                      &
&                  / ( 4.d0 * pi2 * ( one - fexp(-x) ) )
    s_D            = MAX( s_D, epsilon )

!-----------------------------------------------------------------------
!  Compute F_fit
!   (Equation 50 Hannestad & Raffelt)
!-----------------------------------------------------------------------

    F_denom        = ( 3.d0 + ( x - 1.2d0 )**2 + x**(-4) )             &
&                  * ( one + eta_star**2 ) * ( one + y**4 )
    F_fit          = one + one/F_denom

!-----------------------------------------------------------------------
!  Compute G_fit
!   (Equation 50 Hannestad & Raffelt)
!-----------------------------------------------------------------------

    G_fit          = one - 0.0044d0 * x**1.1d0                            &
&                  * y/( 0.8d0 + 0.06d0 * y**1.05d0 )                     &
&                  * sqrt( eta_star )/( eta_star + 0.2d0 )

!-----------------------------------------------------------------------
!  Compute C_fit
!   (Equation 50 Hannestad & Raffelt)
!-----------------------------------------------------------------------

    h_fit          = 0.1d0 * eta_star                                     &
&                  / ( 2.39d0 + 0.1d0 * eta_star**1.1d0 )
    C_fit          = 1.1d0 * x**1.1d0 * h_fit                             &
&                  / ( 2.3d0 + h_fit * x**0.93d0 + 0.0001d0 * x**1.2d0 )  &
&                  * 30.d0/( 30.d0 + 0.005d0 * x**2.8d0 )

!-----------------------------------------------------------------------
!  Compute s_fit
!   (Equation 49 & 50 Hannestad & Raffelt)
!-----------------------------------------------------------------------

    p_fit          = 0.67d0 + 0.18d0 * y**0.4d0
    s_fit(k,kb)    = ( s_ND**(-p_fit) + s_D**(-p_fit) )**( - one/p_fit )  &
&                  * F_fit * ( one + C_fit * G_fit )

  END DO ! k = 1,nez
END DO ! kb = 1,nez

RETURN
END SUBROUTINE s_brem

SUBROUTINE g_brem( y_p, eta_star_p, g_fit )
!-----------------------------------------------------------------------
!
!    File:         g_brem
!    Module:       g_brem
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/23/02
!
!    Purpose:
!      To compute g_fit-component of the kernel for neutrino bremss-
!       trahlung and inelastic neutrino scattering in a nucleon
!       field.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  y           : m_{pi}^{2}/M_{N}c^{2}kT
!  eta_star    : neutron degeneracy parameter
!
!    Output arguments:
!  g           : dimensionless quantity related to S_{sigma}
!
!    Modules:
!  kind_module
!  numerical_module
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE wlExtNumericalModule, ONLY: one

!USE math_functions_module, ONLY: fexp
      
IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(IN)  :: y_p                  ! pion mass parameter
REAL(dp), INTENT(IN)  :: eta_star_p           ! degeneracy parameter

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(OUT) :: g_fit                ! dimensionless fitting parameter

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp)              :: y                    ! use in place of y_p
REAL(dp)              :: eta_star             ! use in place of eta_star_p

REAL(dp), PARAMETER   :: y_min = 1.d-10
REAL(dp), PARAMETER   :: eta_min = 1.d-10
REAL(dp)              :: y2                   ! y*y
REAL(dp)              :: eta_star_1           ! 1/eta_star
REAL(dp)              :: alpha1_denom
REAL(dp)              :: alpha1,alpha2,alpha3
REAL(dp)              :: p1,p2

REAL(KIND=dp)         :: fexp                 ! exponential

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Replace input variables for modification
!-----------------------------------------------------------------------

y                   = y_p
eta_star            = eta_star_p

!-----------------------------------------------------------------------
!  Prevent singular behavior of terms
!-----------------------------------------------------------------------

y                  = MAX( y       , y_min   )
eta_star           = MAX( eta_star, eta_min )

!-----------------------------------------------------------------------
!  Compute alpha1
!   (Equation 53 Hannestad & Raffelt)
!-----------------------------------------------------------------------

y2                 = y * y
eta_star_1         = one/eta_star
alpha1_denom       = 25.d0 * y2 + one
alpha1             = ( 0.5d0 + eta_star_1 )/( one + eta_star_1 )   &
&                  * one/alpha1_denom                              &
&                  + ( 0.5d0 + eta_star/15.6d0 ) * 25.d0 * y2      &
&                 / alpha1_denom

!-----------------------------------------------------------------------
!  Compute alpha2
!   (Equation 53 Hannestad & Raffelt)
!-----------------------------------------------------------------------

alpha2             = ( 0.63d0 + 0.04d0 * eta_star**1.45d0 )        &
&                  / ( one + 0.02d0 * eta_star**2.5d0 )

!-----------------------------------------------------------------------
!  Compute alpha3
!   (Equation 53 Hannestad & Raffelt)
!-----------------------------------------------------------------------

alpha3             = 1.2d0 * fexp( 0.6d0 * eta_star                &
&                  - 0.4d0 * eta_star**1.5d0 )

!-----------------------------------------------------------------------
!  Compute p1
!   (Equation 53 Hannestad & Raffelt)
!-----------------------------------------------------------------------

p1                 = ( 1.8d0 + 0.45d0 * eta_star )                 &
&                  / ( one + 0.15d0 * eta_star**1.5d0 )

!-----------------------------------------------------------------------
!  Compute p2
!   (Equation 53 Hannestad & Raffelt)
!-----------------------------------------------------------------------

p2                 = 2.3d0 - 0.05d0 * eta_star                       &
&                  / ( one + 0.025d0 * eta_star )

!-----------------------------------------------------------------------
!  Compute g_fit
!   (Equation 52 Hannestad & Raffelt)
!-----------------------------------------------------------------------

g_fit              = ( alpha1 + alpha2 * y**p1 )                   &
&                  / ( one + alpha3 * y**p2                        &
&                  + alpha2 * y**( p1 + 2.d0 )/13.75d0 )

RETURN
END SUBROUTINE g_brem


END MODULE HR98_Bremsstrahlung

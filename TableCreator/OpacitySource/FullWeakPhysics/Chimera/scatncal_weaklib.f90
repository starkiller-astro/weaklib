SUBROUTINE scatncal_weaklib &
           ( e_in, e_out, e_in_l, e_in_u, e_out_l, e_out_u, tmev, cmpt, sm, &
             phi0, phi1, cv_N, ca_N )
!-----------------------------------------------------------------------
!
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/07/02
!
!    Modified:     C. Y. Cardall for weaklib
!
!    Purpose:
!      To compute the first two Legendre moments of the
!       neutrino-nucleon scattering kernals.
!
!    Input variables passed through the calling statement
!   e_in      : zone-centered incident neutrino energy [MeV]
!   e_in_l    : inner boundary of the e_in energy zone [MeV]
!   e_in_u    : outer boundary of the e_in energy zone [MeV]
!   e_out     : zone-centered scattered neutrino energy [MeV]
!   e_out_l   : inner boundary of the e_out energy zone [MeV]
!   e_out_u   : outer boundary of the e_out energy zone [MeV]
!   tmev      : temperature [MeV]
!   cmpN      : nucleon chemical potential (including rest mass) [MeV]
!   sm        : rest mass of the nucleon [MeV]
!   cv_N      : neutrino-nucleon weak interaction coefficient
!   ca_N      : neutrino-nucleon weak interaction coefficient
!
!    Output variables passed through the calling statement
!   phi0      : zero moment of the neutrino-nucleon scatering kermel
!   phi1      : first moment of the neutrino-nucleon scatering kermel
!
!    Subprograms called:
!  n_difcs    : computes the neutrino-nucleon inelastic scattering rates
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: zero, half, one, pi
USE physcnst_module, ONLY: hbar, cvel
USE scat_n_module_weaklib

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(in)    :: e_in          ! zone centered incoming neutrino energy [MeV]
REAL(dp), INTENT(in)    :: e_out         ! zone centered incoming neutrino energy [MeV]
REAL(dp), INTENT(in)    :: e_in_l        ! inner zone edge of incoming neutrino energy [MeV]
REAL(dp), INTENT(in)    :: e_in_u        ! outer zone edge of incoming neutrino energy [MeV]
REAL(dp), INTENT(in)    :: e_out_l       ! inner zone edge of outgoing neutrino energy [MeV]
REAL(dp), INTENT(in)    :: e_out_u       ! outer zone edge of outgoing neutrino energy [MeV]
REAL(dp), INTENT(in)    :: tmev          ! temperature [MeV]
REAL(dp), INTENT(in)    :: cmpt          ! nucleon chamical potential
REAL(dp), INTENT(in)    :: sm            ! nucleon mass
REAL(dp), INTENT(in)    :: cv_N          ! neutrino-nucleon vector coupling constant
REAL(dp), INTENT(in)    :: ca_N          ! neutrino-nucleon axial vector coupling constant

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------
REAL(dp), INTENT(out)   :: phi0          ! zero angular moment of the neutrino-nucleon scattering kernal
REAL(dp), INTENT(out)   :: phi1          ! first angular moment of the neutrino-nucleon scattering kernal

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp)                :: d2_sigma        ! double differential cross section

REAL(dp), PARAMETER     :: mult = 5.d0     ! boundary multiplier

REAL(dp)                :: xu_a            ! upper limit of angular integration
REAL(dp)                :: xl_a            ! lower limit of angular integration
REAL(dp)                :: mid_a           ! midpoint of angular integration
REAL(dp)                :: width_a         ! half-width of angular integration
REAL(dp)                :: c_a             ! scaled points of angular quadrature
REAL(dp)                :: costh           ! points of angular quadrature

REAL(dp)                :: t_m             ! paremeter of the energy limits
REAL(dp)                :: prin            ! paremeter of the energy limits
REAL(dp)                :: radical         ! paremeter of the energy limits
REAL(dp)                :: xu_e            ! upper limit of energy integration
REAL(dp)                :: xl_e            ! lower limit of energy integration
REAL(dp)                :: mid_e           ! midpoint of energy integration
REAL(dp)                :: width_e         ! half-width of energy integration
REAL(dp)                :: c_e             ! scaled points of energy quadrature
REAL(dp)                :: e_out_q         ! points of energy quadrature
REAL(dp)                :: phi_e           ! energy integral for a given angle

REAL(dp)                :: xu_ein          ! upper limit of energy integration
REAL(dp)                :: xl_ein          ! lower limit of energy integration
REAL(dp)                :: mid_ein         ! midpoint of energy integration
REAL(dp)                :: width_ein       ! half-width of energy integration
REAL(dp)                :: c_ein           ! scaled points of energy quadrature
REAL(dp)                :: e_in_q          ! points of energy quadrature
REAL(dp)                :: phi0_ein        ! energy integral for a given angle
REAL(dp)                :: phi1_ein        ! energy integral for a given angle

REAL(dp)                :: xu_eot          ! upper limit of energy integration
REAL(dp)                :: xl_eot          ! lower limit of energy integration
REAL(dp)                :: xl_eot_s        ! lower limit of energy integration (temporary)
REAL(dp)                :: mid_eot         ! midpoint of energy integration
REAL(dp)                :: width_eot       ! half-width of energy integration
REAL(dp)                :: c_eot           ! scaled points of energy quadrature
REAL(dp)                :: e_ot            ! points of energy quadrature
REAL(dp)                :: phi0_eot        ! energy integral for a given angle
REAL(dp)                :: phi1_eot        ! energy integral for a given angle

INTEGER                     :: i_a             ! summation index of angular Gauss-Lagendre quadrature
INTEGER                     :: i_e             ! summation index of energy Gauss-Lagendre quadrature
INTEGER                     :: i_ein           ! summation index of energy Gauss-Lagendre quadrature
INTEGER                     :: i_eot           ! summation index of energy Gauss-Lagendre quadrature

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!              \\\\\ COMPUTE SCATTERING FUNCTIONS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                  ||||| Case 1: e_in = e_out |||||
!-----------------------------------------------------------------------

IF ( e_in == e_out ) THEN
print*, '            >>> e_in == e_out'

!-----------------------------------------------------------------------
!  Initialize scattering kernal and integration boundary parameter
!-----------------------------------------------------------------------

  t_m             = DMIN1( 2.d0 * mult * tmev/sm, 0.99d0 )
  phi0            = zero
  phi1            = zero

!-----------------------------------------------------------------------
!  Integrate over final neutrino angle
!-----------------------------------------------------------------------

  xu_a            = one
  xl_a            = - one
  mid_a           = half * ( xu_a + xl_a )
  width_a         = half * ( xu_a - xl_a)

  angle: DO i_a = 1,nquad_a

    c_a           = x_a(i_a) * width_a
    costh         = mid_a + c_a

!-----------------------------------------------------------------------
!  Upper and lower limits of energy integration
!-----------------------------------------------------------------------

    IF ( cmpt - sm < zero ) THEN

      prin        = one - t_m * costh
      radical     = dsqrt( 2.d0 * t_m * ( one - costh ) - t_m * t_m * ( one - costh * costh ) )
      xu_e        = e_in * ( prin + radical )/( one - t_m )
      xl_e        = e_in * ( prin - radical )/( one - t_m )
      xu_e        = min( xu_e, e_in_u )
      xl_e        = max( xl_e, e_in_l   )

    ELSE
    
      xu_e        = e_in_u
      xl_e        = e_in_l

    END IF

!-----------------------------------------------------------------------
!  Integrate over final neutrino energy
!-----------------------------------------------------------------------

    phi_e         = zero
    mid_e         = half * ( xu_e + xl_e )
    width_e       = half * ( xu_e - xl_e )

    energy: DO i_e = 1,nquad_e
      c_e         = x_e(i_e) * width_e
      e_out_q     = mid_e + c_e
         
      IF ( e_out_q <= zero ) CYCLE

!-----------------------------------------------------------------------
!  Get the differential cross section
!-----------------------------------------------------------------------

      CALL n_difcs( e_in, e_out_q, costh, tmev, sm, cmpt, cv_N, ca_N, d2_sigma )

!-----------------------------------------------------------------------
!  Sum outgoing neutrino energies.
!-----------------------------------------------------------------------

      phi_e       = phi_e + d2_sigma * wt_e(i_e)

    END DO energy

!-----------------------------------------------------------------------
!  Sum neutrino scattering angles
!-----------------------------------------------------------------------

    phi0          = phi0 + phi_e * wt_a(i_a) * width_e
    phi1          = phi1 + phi_e * wt_a(i_a) * width_e * costh

  END DO angle

!-----------------------------------------------------------------------
!  Assemble final result
!-----------------------------------------------------------------------

  phi0            = R_nes * C_nes * phi0 * width_a/( e_out**2 * ( e_in_u - e_in_l ) )
  phi1            = R_nes * C_nes * phi1 * width_a/( e_out**2 * ( e_in_u - e_in_l ) )

  RETURN

END IF ! e_in == e_out


!-----------------------------------------------------------------------
!    ||||| Case 2: e_in_l = e_out_u (contiguous energy zones) |||||
!-----------------------------------------------------------------------

IF ( e_in_l == e_out_u ) THEN
print*, '            >>> e_in_l == e_out_u'

!-----------------------------------------------------------------------
!  Initialize scattering kernal and integration boundary parameter
!-----------------------------------------------------------------------

  t_m             = DMIN1( 2.d0 * mult * tmev/sm, 0.99d0 )

  phi0            = zero
  phi1            = zero

!-----------------------------------------------------------------------
!  Integrate over final neutrino energy
!-----------------------------------------------------------------------

  xu_a            = one
  xl_a            = - one
  mid_a           = half * ( xu_a + xl_a )
  width_a         = half * ( xu_a - xl_a)

  angle2: DO i_a = 1,nquad_a

    c_a           = x_a(i_a) * width_a
    costh         = mid_a + c_a

!-----------------------------------------------------------------------
!  Upper and lower limits of energy integration
!-----------------------------------------------------------------------

    IF ( cmpt - sm < zero ) THEN

      prin        = one - t_m * costh
      radical     = DSQRT( 2.d0 * t_m * ( one - costh ) - t_m * t_m * ( one - costh * costh ) )
      xu_e        = e_in_l * ( prin + radical )/( one - t_m )
      xl_eot_s    = e_in_l * ( prin - radical )/( one - t_m )
      xu_e        = DMIN1( xu_e, e_in_u )
      xl_eot_s    = DMAX1( xl_eot_s, e_out_l   )

    ELSE
    
      xu_e        = e_in_u
      xl_eot_s    = e_out_l

    END IF

!-----------------------------------------------------------------------
!  Outer integral - integrate incoming neutrino energy from
!   e_in_l to e_in_u
!-----------------------------------------------------------------------

    phi0_ein      = zero
    phi1_ein      = zero
    xu_ein        = xu_e
    xl_ein        = e_in_l
    mid_ein       = half * ( xu_ein + xl_ein )
    width_ein     = half * ( xu_ein - xl_ein )
    outer_e: DO i_ein = 1,nquad_ein1
      c_ein       = x_ein1(i_ein) * width_ein 
      e_in_q      = mid_ein + c_ein

!-----------------------------------------------------------------------
!  Inner integral - integrate outgoing neutrino energy from
!   e_out_l to e_out_u
!-----------------------------------------------------------------------

      phi0_eot    = zero
      phi1_eot    = zero
      xu_eot      = e_out_u
      xl_eot      = xl_eot_s
      mid_eot     = half * ( xu_eot + xl_eot )
      width_eot   = half * ( xu_eot - xl_eot )
      inner_e: DO i_eot = 1,nquad_eot1
        c_eot     = x_eot1(i_eot) * width_eot
        e_ot      = mid_eot + c_eot

!-----------------------------------------------------------------------
!  Get the differential cross section
!-----------------------------------------------------------------------

        CALL n_difcs( e_in_q, e_ot, costh, tmev, sm, cmpt, cv_N, ca_N, d2_sigma )

!-----------------------------------------------------------------------
!  Sum outgoing neutrino energies.
!-----------------------------------------------------------------------

        phi0_eot  = phi0_eot + d2_sigma * wt_eot1(i_eot) * e_ot**2
        phi1_eot  = phi1_eot + d2_sigma * wt_eot1(i_eot) * e_ot**2

      END DO inner_e

!-----------------------------------------------------------------------
!  Sum incoming neutrino energies.
!-----------------------------------------------------------------------

      phi0_ein   = phi0_ein + phi0_eot * wt_ein1(i_ein) * width_eot * e_in_q**3
      phi1_ein   = phi1_ein + phi0_eot * wt_ein1(i_ein) * width_eot * e_in_q**3

    END DO outer_e

!-----------------------------------------------------------------------
!  Sum neutrino scattering angles
!-----------------------------------------------------------------------

    phi0         = phi0 + phi0_ein * wt_a(i_a) * width_ein
    phi1         = phi1 + phi1_ein * wt_a(i_a) * width_ein * costh 

  END DO angle2

  phi0            = R_nes * C_nes * phi0 * width_a/e_out**2                                      &
&                 / ( e_in**3  * ( e_in_u  - e_in_l  ) * e_out**2 * ( e_out_u - e_out_l ) )
  phi1            = R_nes * C_nes * phi1 * width_a/e_out**2                                      &
&                 / ( e_in**3  * ( e_in_u  - e_in_l  ) * e_out**2 * ( e_out_u - e_out_l ) )

  RETURN

END IF ! e_in_l == e_out_u


!-----------------------------------------------------------------------
!   ||||| Case 3: e_in_l > e_out_u (noncontiguous energy zones) |||||
!-----------------------------------------------------------------------
print*, '            >>> e_in_l > e_out'

phi0              = zero
phi1              = zero

!-----------------------------------------------------------------------
!  Integrate over final neutrino energy
!-----------------------------------------------------------------------

xu_a              = one
xl_a              = - one
mid_a             = half * ( xu_a + xl_a )
width_a           = half * ( xu_a - xl_a)

angle3: DO i_a = 1,nquad_a2

  c_a             = x_a2(i_a) * width_a
  costh           = mid_a + c_a

!-----------------------------------------------------------------------
!  Outer integral - integrate incoming neutrino energy from
!   e_in_l to e_in_u
!-----------------------------------------------------------------------

  phi0_ein        = zero
  phi1_ein        = zero
  xu_ein          = e_in_u
  xl_ein          = e_in_l
  mid_ein         = half * ( xu_ein + xl_ein )
  width_ein       = half * ( xu_ein - xl_ein )
  outer_e2: DO i_ein = 1,nquad_ein2
    c_ein         = x_ein2(i_ein) * width_ein 
    e_in_q        = mid_ein + c_ein

!-----------------------------------------------------------------------
!  Inner integral - integrate outgoing neutrino energy from
!   e_out_l to e_out_u
!-----------------------------------------------------------------------

    phi0_eot      = zero
    phi1_eot      = zero
    xu_eot        = e_out_u
    xl_eot        = e_out_l
    mid_eot       = half * ( xu_eot + xl_eot )
    width_eot     = half * ( xu_eot - xl_eot )
    inner_e2: DO i_eot = 1,nquad_eot2
      c_eot       = x_eot2(i_eot) * width_eot
      e_ot        = mid_eot + c_eot
      
!-----------------------------------------------------------------------
!  Calculate the differential cross section
!-----------------------------------------------------------------------

      CALL n_difcs( e_in_q, e_ot, costh, tmev, sm, cmpt, cv_N, ca_N, d2_sigma )

!-----------------------------------------------------------------------
!  Sum outgoing neutrino energies
!-----------------------------------------------------------------------

      phi0_eot  = phi0_eot + d2_sigma * wt_eot2(i_eot) * e_ot**2
      phi1_eot  = phi1_eot + d2_sigma * wt_eot2(i_eot) * e_ot**2

    END DO inner_e2

!-----------------------------------------------------------------------
!  Sum outgoing neutrino energies.
!-----------------------------------------------------------------------

    phi0_ein     = phi0_ein + phi0_eot * wt_ein2(i_ein) * width_eot * e_in_q**3
    phi1_ein     = phi1_ein + phi0_eot * wt_ein2(i_ein) * width_eot * e_in_q**3

   END DO outer_e2

!-----------------------------------------------------------------------
!  Sum neutrino scattering angles
!-----------------------------------------------------------------------

  phi0           = phi0 + phi0_ein * wt_a2(i_a) * width_ein
  phi1           = phi1 + phi1_ein * wt_a2(i_a) * width_ein * costh 

END DO angle3

phi0              = R_nes * C_nes * phi0 * width_a/e_out**2                                     &
&                 / ( e_in**3  * ( e_in_u  - e_in_l  ) * e_out**2 * ( e_out_u - e_out_l ) )
phi1              = R_nes * C_nes * phi1 * width_a/e_out**2                                     &
&                 / ( e_in**3  * ( e_in_u  - e_in_l  ) * e_out**2 * ( e_out_u - e_out_l ) )

RETURN
END SUBROUTINE scatncal_weaklib

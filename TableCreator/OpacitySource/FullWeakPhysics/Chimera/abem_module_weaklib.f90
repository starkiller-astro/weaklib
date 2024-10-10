MODULE abem_module_weaklib

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: zero, half, one, epsilon, pi
USE physcnst_module, ONLY: cvel, gv, ga
USE polylog_module_weaklib, ONLY: polylog

!-----------------------------------------------------------------------
!  Energy and angle quadrature points and weights
!-----------------------------------------------------------------------

INTEGER, PARAMETER           :: nleg_a = 64   ! number of points of ouitgoing lepton angular Gauss-Lagendre quadrature
REAL(dp), DIMENSION (nleg_a) :: x_a           ! scaled points of angular quadrature
REAL(dp), DIMENSION (nleg_a) :: wt_a          ! scaled weights of angular quadrature

INTEGER, PARAMETER           :: nleg_e = 64   ! number of points of ouitgoing lepton energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION (nleg_e) :: x_e           ! scaled points of energy quadrature
REAL(dp), DIMENSION (nleg_e) :: wt_e          ! scaled weights of energy quadrature

CONTAINS


SUBROUTINE init_quad_abem

  !-----------------------------------------------------------------------
  ! Get integration points and weights
  !-----------------------------------------------------------------------

  CALL gquad(nleg_a,x_a,wt_a,nleg_a)
  CALL gquad(nleg_e,x_e,wt_e,nleg_e)

END SUBROUTINE init_quad_abem


SUBROUTINE nu_N_absr_momts( enu_in, tmev, m_trgt_i, m_trgt_f, m_lep, &
& cmp_trgt_i, cmp_trgt_f, cmp_lep, ab_r0, ab_r1, e_out )
!-----------------------------------------------------------------------
!
!    File:         nu_N_absr_momts
!    Module:       nu_N_absr_momts
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/25/02
!
!    Purpose:
!      To compute moments of the neutrino emission and absorption
!       inverse mean free paths.
!
!    Input arguments:
!   enu_in     : absorbed or emitted neutrino [MeV]
!   tmev       : temperature [MeV]
!   m_trgt_i   : mass of the target particle [MeV]
!   m_trgt_f   : mass of the transformed target particle [MeV]
!   m_lep      : mass of the created lepton [MeV]
!   cmp_trgt_i : neutron chemical potential (including rest mass) [MeV]
!   cmp_trgt_f : proton chemical potential (including rest mass) [MeV]
!   cmp_lep    : electron chemical potential (including rest mass) [MeV]
!
!    Output arguments:
!   ab_r0      : zero angular moment of the inverse absorption mean free
!                 path (cm^{-1})
!   ab_r1      : first angular moment of the inverse absorption mean
!                 free path (cm^{-1})
!   e_out      : mean energy of the emitted lepton [MeV]
!
!    Subprograms called:
!      cc_difcs
!
!    Include files:
!  wlKindModule
!  numerical_module
!  physcnst_module
!
!  abem_module
!
!-----------------------------------------------------------------------

!USE wlKindModule, ONLY: dp
!USE numerical_module, ONLY: zero, half, one, epsilon
!USE physcnst_module, ONLY: cvel, gv, ga

!USE abem_module_weaklib, ONLY:  nleg_a, x_a, nleg_e, x_e, wt_e, wt_a

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(in)    :: enu_in        ! incoming neutrino energy
REAL(dp), INTENT(in)    :: tmev          ! temperature [MeV]
REAL(dp), INTENT(in)    :: m_trgt_i      ! mass of the initial target particle [MeV]
REAL(dp), INTENT(in)    :: m_trgt_f      ! mass of the final target particle [MeV]
REAL(dp), INTENT(in)    :: m_lep         ! mass of the final lepton [MeV]
REAL(dp), INTENT(in)    :: cmp_trgt_i    ! chemical potential of the initial target particle [MeV]
REAL(dp), INTENT(in)    :: cmp_trgt_f    ! chemical potential of the transformed target particle [MeV]
REAL(dp), INTENT(in)    :: cmp_lep       ! chemcal potential of the secondary lepton [MeV]

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(out)   :: ab_r0         ! zero angular moment of the inverse absorption mean free path (cm^{-1})
REAL(dp), INTENT(out)   :: ab_r1         ! first angular moment of the inverse absorption mean free path (cm^{-1})
REAL(dp), INTENT(out)   :: e_out         ! mean energy of the emitted lepton

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp), PARAMETER     :: m_cm = 1.d-2  ! meters/cm
REAL(dp)                :: m_trgt_i2     ! m_trgt_i^{2}
REAL(dp)                :: m_trgt_f2     ! m_trgt_f^{2}
REAL(dp)                :: m_lep2        ! m_lep/cm^{2}
REAL(dp)                :: dm_trgt       ! m_trgt_i - m_trgt_f
REAL(dp)                :: dm_trgtp      ! - dm_trgt

INTEGER                     :: i_a           ! summation index of angular Gauss-Lagendre quadrature
REAL(dp)                :: xu_a          ! upper limit of lepton angular quadrature
REAL(dp)                :: xl_a          ! lower limit of lepton angular quadrature
REAL(dp)                :: mid_a         ! midpoint of lepton angular quadrature
REAL(dp)                :: width_a       ! half-width of lepton angular quadrature
REAL(dp)                :: c_a           ! scaled points of angular quadrature
REAL(dp)                :: costh         ! points of angular quadrature

REAL(dp)                :: ab_r0_e       ! zero moment of the neutrino absorption inverse mean free path per angle
REAL(dp)                :: ab_r1_e       ! first moment of the neutrino absorption inverse mean free path per angle
REAL(dp)                :: e_out_e       ! parameter to estimate energy quadrature widths per angle

REAL(dp), PARAMETER     :: mult = 10.d0  ! boundary multiplier
REAL(dp)                :: t_m           ! paremeter of the energy limits
REAL(dp)                :: prin          ! paremeter of the energy limits
REAL(dp)                :: radical       ! paremeter of the energy limits

INTEGER                     :: i_e           ! summation index of ouitgoing lepton energy Gauss-Lagendre quadrature
REAL(dp)                :: xu_e          ! upper limit of lepton energy quadrature
REAL(dp)                :: xl_e          ! lower limit of lepton energy quadrature
REAL(dp)                :: mid_e         ! midpoint of lepton energy quadrature
REAL(dp)                :: width_e       ! half-width of lepton energy quadrature
REAL(dp)                :: c_e           ! scaled points of energy quadrature
REAL(dp)                :: e_lep_f       ! energy quadrature points

REAL(dp)                :: absr_eomega   ! absorption cross section per unit energy and solid angle

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

m_trgt_i2          = m_trgt_i * m_trgt_i
m_trgt_f2          = m_trgt_f * m_trgt_f
m_lep2             = m_lep * m_lep
dm_trgt            = m_trgt_i - m_trgt_f
dm_trgtp           = - dm_trgt

!-----------------------------------------------------------------------
!  Initialize for dp (angle and energy) integration
!-----------------------------------------------------------------------

ab_r0              = zero
ab_r1              = zero
e_out              = zero

!-----------------------------------------------------------------------
!
!           \\\\\ INTEGRATION OVER FINAL LEPTON ANGLE /////
!
!-----------------------------------------------------------------------

xu_a               = one
xl_a               = - one
mid_a              = half * ( xu_a + xl_a )
width_a            = half * ( xu_a - xl_a)
t_m                = 2.d0 * mult * tmev/m_trgt_i
t_m                = DMIN1( t_m, 0.999d0 )


DO i_a = 1,nleg_a
  c_a              = x_a(i_a) * width_a
  costh            = mid_a + c_a

!-----------------------------------------------------------------------
!  Initialize for final lepton energy integration
!-----------------------------------------------------------------------

  ab_r0_e          = zero
  ab_r1_e          = zero
  e_out_e          = zero

!-----------------------------------------------------------------------
!
!          \\\\\ INTEGRATION OVER GFINAL LEPTON ENERGY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine the energy quadrature limits
!-----------------------------------------------------------------------

  prin            = one - t_m * costh
  radical         = DSQRT( DMAX1( 2.d0 * t_m * ( one - costh ) - t_m * t_m * ( one - costh * costh ), epsilon ) )
  xu_e            = enu_in * ( prin + radical )/( one - t_m ) + dm_trgt
  xl_e            = enu_in * ( prin - radical )/( one - t_m ) + dm_trgt
  xl_e            = DMAX1( xl_e, zero   )

!      IF ( i_a == 2   .or.  i_a == 32  .or.  i_a == 64  .or.  i_a == 16  .or.  i_a == 48 ) &
!& WRITE (6,3005) t_m,prin,radical,xu_e,xl_e
! 3005 format (' t_m=',1pe10.3,' prin=',1pe10.3,' radical=',1pe10.3,' xu_e=',1pe10.3,' xl_e=',1pe10.3)

!-----------------------------------------------------------------------
!  Integrate over final lepton energy
!-----------------------------------------------------------------------

  mid_e            = half * ( xu_e + xl_e )
  width_e          = half * ( xu_e - xl_e )

  DO i_e = 1,nleg_e
    c_e            = x_e(i_e) * width_e
    e_lep_f        = mid_e + c_e
    IF ( e_lep_f <= zero ) CYCLE
    CALL cc_difcs( enu_in, e_lep_f, costh, tmev, m_trgt_i, dm_trgtp, &
&    cmp_trgt_i, cmp_trgt_f, cmp_lep, gv, ga, absr_eomega )

    absr_eomega    = DMAX1( absr_eomega, zero )
    ab_r0_e        = ab_r0_e + absr_eomega * wt_e(i_e)
    ab_r1_e        = ab_r1_e + absr_eomega * wt_e(i_e) * costh
    e_out_e        = e_out_e + absr_eomega * wt_e(i_e) * e_lep_f
  END DO ! i_e = 1,nleg_e

  ab_r0_e          = ab_r0_e * width_e
  ab_r1_e          = ab_r1_e * width_e
  e_out_e          = e_out_e * width_e

!-----------------------------------------------------------------------
!  Sum final result over angle and average the final electron energy
!-----------------------------------------------------------------------

  ab_r0            = ab_r0 + ab_r0_e * wt_a(i_a)
  ab_r1            = ab_r1 + ab_r1_e * wt_a(i_a)
  e_out            = e_out + e_out_e * wt_a(i_a)

END DO ! i_a = 1,nleg_a

e_out              = e_out/( ab_r0 + epsilon )
ab_r0              = m_cm * ab_r0 * width_a
ab_r1              = m_cm * ab_r1 * width_a

RETURN
END SUBROUTINE nu_N_absr_momts

!-------------------------------------------------------------------
      SUBROUTINE cc_difcs(e,ep,costh,T,s,dm,mun,mup,mue,gv,ga,fwqt)
!-------------------------------------------------------------------
!-    THIS ROUTINE CALCULATES DIFF CROSS SECTION /VOLUME
!-        FOR  CHARGED CURRENT REACTIONS
!     input: e = incoming neutrino energy
!            ep = outgoing neutrino energy
!            costh = cosine of the angle between in and out neutrino
!                    momentum vectors
!            T = temperature
!            s = mass of the target particle
!            mun = neutron chemical potential (including rest mass)
!            mup = proton chemical potential (  "   )
!            mue = electron chemical potential
!            dm = neutron - proton mass difference
!            gv = vector coupling
!            ga = axial coupling
!            fwqt = (1/V) dsigma/(d(costh) d(ep))
!    UNITS: if all energies are input in MeV the output fwqt has
!           dimensions of MeV^-1 m^-1
!-------------------------------------------------------------------

!USE wlKindModule, ONLY: dp
!USE numerical_module, ONLY: &
!      zero, one, epsilon, pi
!USE physcnst_module, ONLY: &
!      cvel, hbar, gv, ga, mn, mp, me, Gw, kmev, dmnp, rmu

      Implicit none

      real(dp) :: e,ep,costh,fwqt,qo,qu2,q,T,s,dm,gv,ga
      real(dp) :: forfac,mun,mup,mue,muh,efac,beta
      real(dp), parameter :: pi2=9.8696044d0
      real(dp), parameter :: G2= 6.79d-10
!
!                 G2 in MeV^-5 * m^-1
!

      real(dp) eminus,l,u,z,arg1,arg2,pd,fd
      real(dp) I0, I1, I2, uli1, uli2, uli3, bli1, bli2, bli3
      real(dp) zi1, zi2, zi3, I2U, I2D, impl, impt, impa, impva
      real(dp) r1, r2, r3, A, B, C

      real(dp)                :: fexp

      external fexp
!------------------------------
!--      KINEMATICAL FACTORS
!------------------------------
      qo = e - ep
      qu2 = 2.d0*E*EP*(costh-1.d0)
      q  = dsqrt(E*E + EP*EP - 2.d0*E*EP*costh)
      efac = (mue - ep)/T
      muh = mun - mup

      if(efac.gt.20.d0.or.qu2.ge.0.d0)then
         fwqt = 0.d0
         return
      else
         fd = 1.d0/(1.d0+fexp(efac))
      end if
!------------------------------
!     -      RESPONSE FUNCTIONS
!------------------------------

      beta = 1.d0 - 2.d0*s*dm/qu2

      eminus = 0.5d0*(-qo*beta + q*dsqrt(beta**2 - 4.d0*(s*s/qu2)))

      if(eminus.lt.s) eminus = s

      l = eminus/T
      u = mun/T
      z = (qo+muh)/T

      arg1 = l-u
      arg2 = l-u+z

      if(((arg1.gt.20.).and.(arg2.gt.20.)).or.(z.lt.-25.))then

         I0 = 0.d0
         I1 = 0.d0
         I2 = 0.d0
         pd = 0.d0

      else

         if (dabs(z).lt.1.d-8) then

            call polylog(arg1,uli1,uli2,uli3)

            zi1 = -1.d0/(1.d0+fexp(-arg1))
            zi2 = dlog(1.d0+fexp(arg1))
            zi3 = -uli2

            pd = 1.d0/(1.d0-0.5d0*z)

         else

            call polylog(arg1,uli1,uli2,uli3)
            call polylog(arg2,bli1,bli2,bli3)

            zi1 = (uli1 - bli1)/z
            zi2 = (uli2 - bli2)/z
            zi3 = (uli3 - bli3)/z
            pd = z/(1.d0-fexp(-z))

         endif

         I0 = T*(1.d0 + zi1)
         I1 = T*T*(u - 0.5d0*z + (zi2 + l*zi1))
         I2 = T*T*T*( u*(u-z) + (pi2+z*z)/3.d0 - 2.d0*zi3 &
            + l*(2.d0*zi2 + l*zi1) )

      endif

      if(q.gt.1.d-3)then

         forfac = 1.d0
!------------------------
!---  long vector:
!------------------------
         impl = qu2*(I2 + qo*I1 + 0.25d0*qu2*I0 )/(2.d0*pi*q**3)
!------------------------
!---  tran vector:
!------------------------
         impt = 0.5d0*impl + (s*s+0.5d0*qu2)*I0/(4.d0*pi*q)
!------------------------
!---  axial correction
!------------------------
         impa = s*s*I0/(2.d0*pi*q)
         impva = qu2*(qo*I0+2.d0*I1)/(8.d0*pi*q**3)

      else

         forfac = dsqrt(0.5*(costh-1.d0)/(E*EP))

         impl =  forfac*(I2 + qo*I1 + 0.25d0*qu2*I0 )/(2.d0*pi)

         impt = 0.5d0*impl + forfac*(s*s+0.5d0*qu2)*I0/(4.d0*pi)

         impa = forfac * s*s*I0/(2.d0*pi)

         impva = forfac*(qo*I0+2.d0*I1)/(8.d0*pi)

      end if

      R1  = (gv*gv+ga*ga)*(impt+impl)
      R2  = gv*gv*impt+ga*ga*(impt-impa)
      R3  = 2.d0*gv*ga*impva

      A = (2.d0*E*EP+0.5d0*qu2)/(q*q)
      B = E + EP

      if(q.gt.1.d-3)then
         A = (2.d0*E*EP+0.5d0*qu2)/(q*q)
         B = E + EP
         fwqt = (G2/pi2)*(costh-1.d0)*pd*(A*R1 + R2 + B*R3)
      else
         A = (2.d0*E*EP+0.5d0*qu2)/(q*q)
         B = E + EP
         fwqt = (G2/pi2)* pd *(A*R1 + R2 + B*R3)
      end if

      fwqt = fwqt*fd*ep**2

      return
      end

!-- moved to polylog_module_weaklib.f90
! !--------------------------------------------------
! SUBROUTINE POLYLOG(x,pl1,pl2,pl3)
! !--------------------------------------------------
!   use polylog_module_weaklib
!   !USE wlKindModule, ONLY: dp
!   !USE numerical_module, ONLY: &
!   !      zero, one, epsilon, pi
!   !USE physcnst_module, ONLY: &
!   !      cvel, hbar, gv, ga, mn, mp, me, Gw, kmev, dmnp, rmu

!   implicit none

!   real(dp) :: x,pl1,pl2,pl3
!   integer j
!   real(dp) :: pi2,arg,a,b
!   parameter (pi2=9.8696044d0)

!   real(dp)                :: fexp

!   external fexp

!   pl1 = x+dlog(1.d0+fexp(-x))

!   if (x.ge.4.6) then
!     pl2 = -(0.5d0*x*x + pi2/6.d0)/(1.d0+fexp(-1.5*x))
!     pl3 = -x*(x*x + pi2)/6.d0/(1.d0+fexp(-1.69*x))
!   else
!     arg=fexp(x)
!     j = int(10.d0*arg)
!     a = 10.d0*arg - dble(j)
!     b = 1.d0 - a
!     pl2 = pl2a(j)*b+pl2a(j+1)*a
!     pl3 = pl3a(j)*b+pl3a(j+1)*a
!   endif

!   return

! end

!       !------------------------------------------------------------------


END module abem_module_weaklib

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

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: &
      zero, one, epsilon, pi
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
         fd = 1.d0/(1.d0+dexp(efac))
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

            zi1 = -1.d0/(1.d0+dexp(-arg1))
            zi2 = dlog(1.d0+dexp(arg1))
            zi3 = -uli2

            pd = 1.d0/(1.d0-0.5d0*z)

         else

            call polylog(arg1,uli1,uli2,uli3)
            call polylog(arg2,bli1,bli2,bli3)

            zi1 = (uli1 - bli1)/z
            zi2 = (uli2 - bli2)/z
            zi3 = (uli3 - bli3)/z
            pd = z/(1.d0-dexp(-z))

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

!--------------------------------------------------
SUBROUTINE POLYLOG(x,pl1,pl2,pl3)
!--------------------------------------------------
  use polylog_module_weaklib
  USE wlKindModule, ONLY: dp
  USE numerical_module, ONLY: &
        zero, one, epsilon, pi
  !USE physcnst_module, ONLY: &
  !      cvel, hbar, gv, ga, mn, mp, me, Gw, kmev, dmnp, rmu

  implicit none

  real(dp) :: x,pl1,pl2,pl3
  integer j
  real(dp) :: pi2,arg,a,b
  parameter (pi2=9.8696044d0)

  pl1 = x+dlog(1.d0+dexp(-x))

  if (x.ge.4.6) then
    pl2 = -(0.5d0*x*x + pi2/6.d0)/(1.d0+dexp(-1.5*x))
    pl3 = -x*(x*x + pi2)/6.d0/(1.d0+dexp(-1.69*x))
  else
    arg=dexp(x)
    j = int(10.d0*arg)
    a = 10.d0*arg - dble(j)
    b = 1.d0 - a
    pl2 = pl2a(j)*b+pl2a(j+1)*a
    pl3 = pl3a(j)*b+pl3a(j+1)*a
  endif

  return

end

      !------------------------------------------------------------------

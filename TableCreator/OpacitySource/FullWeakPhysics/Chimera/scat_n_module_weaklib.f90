MODULE scat_n_module_weaklib

USE wlKindModule, ONLY: dp
USE polylog_module_weaklib, ONLY: polylog
USE numerical_module, ONLY: one, pi

!-----------------------------------------------------------------------
!  Constants
!-----------------------------------------------------------------------

REAL(dp), PARAMETER     :: m_cm = 1.d-2         ! meters/cm
REAL(dp)                :: C_nes                ! converts d sigma/ V d Omega dE3 from m^-1 to cm^-1
REAL(dp)                :: R_nes                ! converts d sigma/ V d Omega dE3 to R(E1,E3,cosQ)

!-----------------------------------------------------------------------
!  Gauss-Legendre quantities for scatncal
!-----------------------------------------------------------------------

INTEGER, PARAMETER               :: nquad_a = 32    ! number of points of angular Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_a) :: x_a             ! normalized points of angular Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_a) :: wt_a            ! weights of angular Gauss-Lagendre quadrature

INTEGER, PARAMETER               :: nquad_e = 32    ! number of points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_e) :: x_e             ! points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_e) :: wt_e            ! weights of energy Gauss-Lagendre quadrature

INTEGER, PARAMETER               :: nquad_ein1 = 8  ! number of points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_ein1) :: x_ein1       ! points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_ein1) :: wt_ein1      ! weights of energy Gauss-Lagendre quadrature

INTEGER, PARAMETER               :: nquad_eot1 = 8  ! number of points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_eot1) :: x_eot1       ! points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_eot1) :: wt_eot1      ! weights of energy Gauss-Lagendre quadrature

INTEGER, PARAMETER               :: nquad_ein2 = 4  ! number of points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_ein2) :: x_ein2       ! points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_ein2) :: wt_ein2      ! weights of energy Gauss-Lagendre quadrature

INTEGER, PARAMETER               :: nquad_eot2 = 4  ! number of points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_eot2) :: x_eot2       ! points of energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_eot2) :: wt_eot2      ! weights of energy Gauss-Lagendre quadrature

INTEGER, PARAMETER               :: nquad_a2 = 4    ! number of points of angular Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_a2) :: x_a2           ! normalized points of angular Gauss-Lagendre quadrature
REAL(dp), DIMENSION(nquad_a2) :: wt_a2          ! weights of angular Gauss-Lagendre quadrature

CONTAINS


SUBROUTINE init_quad_scat_n

  !-----------------------------------------------------------------------
  !  Initialize
  !
  !  C_nes converts d sigma/ V d Omega dE3 from /m to /cm
  !  R_nes converts d sigma/ V d cos theta to dE3 to d sigma/ V d Omega dE3
  !  Get quadrature points and weights
  !-----------------------------------------------------------------------

  C_nes                     = m_cm
!  R_nes                     = one/( 2.d0 * pi )
  R_nes                     = one / ( 2.d0 * pi )**2
                              !-- Additional 1/2pi for Chimera vs. thornado

  CALL gquad( nquad_a,    x_a,    wt_a,    nquad_a    )
  CALL gquad( nquad_e,    x_e,    wt_e,    nquad_e    )
  CALL gquad( nquad_ein1, x_ein1, wt_ein1, nquad_ein1 )
  CALL gquad( nquad_eot1, x_eot1, wt_eot1, nquad_eot1 )
  CALL gquad( nquad_ein2, x_ein2, wt_ein2, nquad_ein2 )
  CALL gquad( nquad_eot2, x_eot2, wt_eot2, nquad_eot2 )
  CALL gquad( nquad_a2,   x_a2,   wt_a2,   nquad_a2   )

END SUBROUTINE init_quad_scat_n

 
!-------------------------------------------------------------------   
      SUBROUTINE N_DIFCS(E,EP,costh,T,s,mun,cv,ca,fwqt)
!-------------------------------------------------------------------       
!-    THIS ROUTINE CALCULATES DIFF CROSS SECTION /VOLUME
!-        FOR  NEUTRAL CURRENT REACTIONS
!     input: E = incoming neutrino energy
!            EP = outgoing neutrino energy
!            costh = cosine of the angle between in and out neutrino 
!                    momentum vectors
!            T = temperature
!            s = mass of the target particle
!            mun = chemical potential (including rest mass of the target)
!            cv = vector coupling
!            ca = axial coupling
!            fwqt = (1/V) dsigma/(d(costh) d(EP))
!    UNITS: if all energies are input in MeV the output fwqt has 
!           dimensions of MeV^-1 m^-1
!-------------------------------------------------------------------       
      Implicit none 
      
      real(dp) :: E,EP,costh,fwqt,qo,qu2,q,T,s,mun,cv,ca
      real(dp), parameter :: pi=3.141592654d0, pi2=9.8696044d0, G2= 6.8945d-10
! 
!                 G2 in MeV^-5 * m^-1
!

      real(dp) :: eminus,l,u,z,arg1,arg2,pd
      real(dp) :: I0, I1, I2, uli1, uli2, uli3, bli1, bli2, bli3
      real(dp) :: zi1, zi2, zi3, impl, impt, impa, impva
      real(dp) :: r1, r2, r3, A, B

print*, '               >>> n_difcs'
print*, '               >>> E     ', E
print*, '               >>> EP    ', EP
print*, '               >>> costh ', costh
print*, '               >>> T     ', T
print*, '               >>> s     ', s
print*, '               >>> mun   ', mun
!------------------------------
!--      KINEMATICAL FACTORS          
!------------------------------
      qo = E - EP
      qu2 = 2.d0*E*EP*(costh-1.d0)
      q  = dsqrt(E*E + EP*EP - 2.d0*E*EP*costh)

!------------------------------
!-      RESPONSE FUNCTIONS
!------------------------------

        eminus = 0.5d0*(-qo + q*dsqrt(1.d0 - 4.d0*(s*s/qu2)))
print*, '               >>> eminus', eminus
        if(eminus.lt.s) eminus = s
print*, '               >>> eminus', eminus
   
        l = eminus/T 
        u = mun/T
        z = qo/T
print*, '               >>> l     ', l
print*, '               >>> u     ', u
print*, '               >>> z     ', z

        arg1 = l-u
        arg2 = l-u+z
print*, '               >>> arg1  ', arg1
print*, '               >>> arg2  ', arg2

        if(((arg1.gt.25.).and.(arg2.gt.25.)).or.(z.lt.-25.))then
print*, '               >>> 1'
           I0 = 0.d0
           I1 = 0.d0
           I2 = 0.d0
           pd = 0.d0
        else

           if (dabs(z).lt.1.d-3) then
print*, '               >>> 2a'

             call polylog(arg1,uli1,uli2,uli3)

             zi1 = -1.d0/(1.d0+dexp(-arg1))
             zi2 = dlog(1.d0+dexp(arg1))
             zi3 = -uli2
             pd = 1.d0/(1.d0-0.5d0*z)

           else
print*, '               >>> 2b'
           
             call polylog(arg1,uli1,uli2,uli3)
             call polylog(arg2,bli1,bli2,bli3)
print*, '               >>> uli1', uli1
print*, '               >>> uli2', uli2
print*, '               >>> uli3', uli3
print*, '               >>> bli1', bli1
print*, '               >>> bli2', bli2
print*, '               >>> bli3', bli3

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

!------------------------   
!---  long vector:
!------------------------  
print*, '               >>> I0', I0
print*, '               >>> I1', I1
print*, '               >>> I2', I2
        impl = qu2*(I2 + qo*I1 + 0.25d0*qu2*I0 )/(2.d0*pi*q**3)
!------------------------  
!---   tran vector:
!------------------------  
        impt = 0.5d0*impl + (s*s+0.5d0*qu2)*I0/(4.d0*pi*q)
!------------------------  
!---   axial correction
!------------------------  
        impa = s*s*I0/(2.d0*pi*q)
        impva = qu2*(qo*I0+2.d0*I1)/(8.d0*pi*q**3)
        
        R1  = (cv*cv+ca*ca)*(impt+impl)
        R2  = cv*cv*impt+ca*ca*(impt-impa)
        R3  = 2.d0*cv*ca*impva

        A = (2.d0*E*EP+0.5d0*qu2)/(q*q)
        B = E + EP

        fwqt = (G2/pi2)*(costh-1.d0)*pd*(A*R1 + R2 + B*R3) &
               *EP*EP
         
      return
      end

!-- moved to polylog_module_weaklib.f90
! !--------------------------------------------------   
!       SUBROUTINE POLYLOG(x,pl1,pl2,pl3)
! !--------------------------------------------------
!       use polylog_module_weaklib
!       implicit none
!       real(dp) :: x,pl1,pl2,pl3
!       integer j
!       real(dp) :: arg,a,b
!       real(dp), parameter :: pi2=9.8696044d0

!       pl1 = x+dlog(1.d0+dexp(-x))
    
!       if (x.ge.4.6) then
!         pl2 = -(0.5d0*x*x + pi2/6.d0)/(1.d0+dexp(-1.5*x))
!         pl3 = -x*(x*x + pi2)/6.d0/(1.d0+dexp(-1.69*x))
!       else
!         arg=dexp(x) 
!         j = int(10.d0*arg)
!         a = 10.d0*arg - dble(j)
!         b = 1.d0 - a
!         pl2 = pl2a(j)*b+pl2a(j+1)*a
!         pl3 = pl3a(j)*b+pl3a(j+1)*a
!       endif

!       return
!       end

! !------------------------------------------------------------------

END module scat_n_module_weaklib

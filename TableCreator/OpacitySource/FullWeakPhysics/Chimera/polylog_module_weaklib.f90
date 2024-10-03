MODULE polylog_module_weaklib

  USE wlKindModule, ONLY: dp

  implicit none

  integer, parameter :: nt=1000
  real(dp)   :: xa(0:nt), pl2a(0:nt), pl3a(0:nt) !Polylog table

CONTAINS

!-- Version previously in scat_n_module_weaklib.f90
!--------------------------------------------------   
      SUBROUTINE POLYLOG(x,pl1,pl2,pl3)
!--------------------------------------------------
!      use polylog_module_weaklib
      implicit none
      real(dp) :: x,pl1,pl2,pl3
      integer j
      real(dp) :: arg,a,b
      real(dp), parameter :: pi2=9.8696044d0

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

!-- Version previously in abem_module_weaklib.f90
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

end MODULE polylog_module_weaklib

!!! modified date: June 17, 2020
!!! Authors: Codes wirtten by Gang Guo
!!! neutrino CC opacity is given by 4D integrals using CUBA library
!!! processes to be considered in this routine
!!! opt=1: v + n -> e- + p  (oapcity in 1/km)
!!! opt=2: vb + p -> e+ + n (opacity in 1/km)
!!! opt=3: vb + p + e- -> n (opacity in 1/km)
!!! opt=4: n -> vb + e- + p [emissivity in v/(s cm^3 MeV)]

!!! globle modules for cuba high dimension integrals
!!! please adjust 'nstart' & 'nincrease' in the module
!!! to trade off between effiency & accuracy

module CUBA
  implicit none
  integer ndim, ncomp, nvec, last, seed, mineval, maxeval
  double precision epsrel, epsabs, userdata
  parameter (ndim = 4)
  parameter (ncomp = 1)
  parameter (userdata = 0)
  parameter (nvec = 1)
  parameter (epsrel = 1D-3)
  parameter (epsabs = 1D-5)
  parameter (last = 4)
  parameter (seed = 0)
  parameter (mineval = 0)
  parameter (maxeval = 1000000)

  integer nstart, nincrease, nbatch, gridno
  integer*8 spin
  character*(*) statefile
  parameter (nstart = 50000)	 !!! to be tuned by user
  parameter (nincrease = 2000)	 !!! to be tuned by user
  parameter (nbatch = 1000)
  parameter (gridno = 0)
  parameter (statefile = "")
  parameter (spin = -1)

  integer nnew
  double precision flatness
  parameter (nnew = 300000)
  parameter (flatness = 25D0)

  integer key1, key2, key3, maxpass
  double precision border, maxchisq, mindeviation
  integer ngiven, ldxgiven, nextra
  parameter (key1 = 47)
  parameter (key2 = 1)
  parameter (key3 = 1)
  parameter (maxpass = 5)
  parameter (border = 1D-2)
  parameter (maxchisq = 10D0)
  parameter (mindeviation = .25D0)
  parameter (ngiven = 0)
  parameter (ldxgiven = ndim)
  parameter (nextra = 0)

  integer key
  parameter (key = 9)
  external integrand,integrand_D
  double precision integral(ncomp), error(ncomp), prob(ncomp)
  integer:: verbose=0, nregions, neval, fail
end module CUBA

module global_4D
  implicit none
  real*8,parameter :: pi=3.1415927d0,Gf=1.166d-11,Vud=&
       0.97427d0, gA0 =1.2723d0,gV0=1.d0, Mpi=139.57,Mnp=938.919d0,&
       Dnp=1.293, MA=1.0d3,MV=840.d0 
  real*8:: me=0.511d0,mn=939.565d0,mp=938.272d0,&
           Un=-20.546d0,Up=-36.853d0
  real*8:: pe(0:3),pv(0:3),pn(0:3),pp(0:3),q(0:3),Enu,pn3,En,cthe,&
       q3,Jacob
  real*8:: ye,rho,Tem,nden,mun,mup,mue
  real*8:: jVA=1.d0,JAF=1.d0
  real*8:: pnmax, pnmin, Tfac=100.
  real*8:: ppmax,ppmin,pemax,pemin
  integer,parameter:: opt_pn=1
  integer:: opt_form,opt_pseudo,opt_wm,opt0
end module global_4D

!========================================================================
! Routine to calculate the neutrino CC opacity via 4D integrals
! (same inputs & output as the 2D subroutine)
! Inputs:
! (1) opt_had=0,1,2,3 for LO, LO+weak magnetism(WM),
! LO+WM+pseduoscalar(PS), LO+WM+PS+form factor dependencies
! (2) opt=opt0=1,2,3,4 for nu capture on n, nub capture on p, inverse
! neutron decay, neutron decay emissivity
! (3) EnuA(NP): neutrino energy array in MeV with dimension NP
! (4) xTem: Temperature in MeV
! (5) cheme, masse: relativistic chemical potentials & masses of
! e^-(mu^-) in MeV
! (6) chemn,chemp,massn,massp,xUn,xUp: relativistic chemical
! potentials, masses, & potentials of neutron and proton in MeV
! relativsitic disperstion relation for nucleons,
! En=sqrt(p^2+mn^2)+Un
! Outputs:
! (7) OpaA(NP): opacities in km^-1; for neutron decay, emissivities
! in (s cm^3 MeV)^-1
!     **** blocking factor of final state neutrinos in neutron decay
!     is not considered ****
!=============================================================
subroutine Opacity_CC_4D_GSI(opt_had,opt,NP,EnuA,OpaA,xTem,cheme,chemn,&
     chemp,masse,massn,massp,xUn,xUp)
  use global_4D
  implicit none
  integer:: opt_had,opt,NP,ii
  real*8:: EnuA(NP),OpaA(NP),xTem,cheme,chemn,chemp,masse,massn,massp,&
      xUn,xUp
  real*8,external:: Integration,integration_D

  ! CAREFUL. IN THIS ROUTINE A LOT OF NUMBERS ARE DEFINED AS REAL(4) INSTEAD OF REAL(8)
  ! E.G. 1.0 INSTEAD OF 1.0d0

  Tem = xTem
  mn = massn
  mp = massp
  me = masse
  mun = chemn
  mup = chemp
  mue = cheme
  Un = xUn
  Up = xUp
  jVA = 1.d0
  jAF = 1.d0

  select case(opt_had)
  case(0)
     opt_form = 0
     opt_pseudo = 0
     opt_wm = 0
  case (1)
     opt_wm = 1
     opt_form = 0
     opt_pseudo = 0
  case (2)
     opt_wm = 1
     opt_pseudo = 1
     opt_form = 0
  case (3)
     opt_wm = 1
     opt_pseudo = 1
     opt_form = 1
  end select

  Tfac = 100
  OpaA = 0.0
  opt0 = opt
  select case(opt0)
  case (1)
     do ii=1,NP
	Enu = EnuA(ii)
	OpaA(ii) = Integration()
     end do
  case (2)
     call nubp_4D()
     do ii=1,NP
	Enu = EnuA(ii)
	OpaA(ii) = Integration()
     end do
  case (3,4)
     call nubp_4D()
     mue = cheme
     do ii=1,NP
	Enu = EnuA(ii)
	OpaA(ii) = Integration_D()
     end do
  end select
end subroutine Opacity_CC_4D_GSI


function Integration()
  use global_4D
  use CUBA
  implicit none
  real*8:: Integration,tmp,tmp1,At,tmp2,tmp3,del,mpt

  Integration = 0.

  pnmin = 0.
  pnmax = sqrt( (Tfac*Tem + mun - Un)**2 - mn**2 )

  mpt = mp + me
  if(mpt<mn) then
     tmp = (1.d0+mpt/mn)/(1.d0-mpt**2/mn**2)*Enu
  else if(mpt>=mn) then
     tmp = pnmax
  end if
  tmp1 = sqrt( (tmp-Enu)**2+mpt**2 )+Up-sqrt(tmp**2+mn**2)-Un
  if(Enu<=tmp1) goto 12

  call Range_pn_4D()
  if( pnmin > pnmax ) goto 12
!==========    vegas  ==============================================
! by default vegas is always used
   call vegas(ndim, ncomp, integrand, userdata, nvec, &
      epsrel, epsabs, verbose, seed,  &
      mineval, maxeval, nstart, nincrease, nbatch,    &
      gridno, statefile, spin, &
      neval, fail, integral, error, prob)

  Integration = integral(1)/197.327d-18*2.d0
12  continue
  return
end function Integration

function Integration_D()
  use global_4D
  use CUBA
  implicit none
  real*8:: Integration_D,tmp,tmp1,At,tmp2,tmp3,del,mpt

  Integration_D = 0.0

  call Range_pn_dec_4D()
  if( pnmin > pnmax ) goto 12
!==========    vegas  ==============================================
  call vegas(ndim, ncomp, integrand_D, userdata, nvec, &
     epsrel, epsabs, verbose, seed,	 &
     mineval, maxeval, nstart, nincrease, nbatch,&
     gridno, statefile, spin,	 &
     neval, fail, integral, error, prob)

  if(opt0.eq.3) then
     Integration_D = integral(1)/197.327d-18*2.d0
  else if(opt0.eq.4) then
     Integration_D = integral(1)*Enu**2/(2.d0*pi**2)*2.d0&
	  *(1.d13/197.327)**3*(3.d23/197.327)
  end if
12  continue
  return
end function Integration_D


integer function integrand(ndim, xx, ncomp, ff)
  use global_4D
  implicit none
  integer:: ndim,ncomp
  real*8:: xx(*),ff(*)
  real*8,external:: amp_4D
  real*8:: xf2,xf3,xf4

  call initmomentum(xx)
  ff(1)=amp_4D()*q3*pn3/(16.d0*(2.d0*pi)**4*Enu**2*pp(0))*Jacob
  xf2 = 1.0/( exp((pn(0)+Un-mun)/Tem) + 1.0 )
  xf3 = 1.0/( exp((pe(0)-mue)/Tem) + 1.0 )
  xf4 = 1.0/( exp((pp(0)+Up-mup)/Tem) + 1.0 )
  ff(1) = ff(1)*xf2*(1.d0-xf4)*(1.d0-xf3)
  integrand = 0
end function integrand


subroutine initmomentum(x)
  use global_4D
  implicit none
  real*8:: x(4),En_min,En_max,tmp,tmp1,q3_max,q3_min,&
   cthe_max,cthe_min,cthe_n,Ep,pp3,q0,Ee,Pe3,tmp2,q3_lim(2),&
   q0_max,q0_min,tmp3,tmp4,cthe_v,phi_n

  En_min = max( sqrt(pnmin**2+mn**2)+Un, mp+Up-Enu+me)
  En_max = sqrt(pnmax**2+mn**2)+Un
  En = (En_max-En_min)*x(1) + En_min
  pn3 = sqrt( (En-Un)**2 - mn**2 )
  tmp = sqrt( (En-Up+Enu-me)**2 - mp**2 )
  tmp1 = pn3 + tmp
  tmp2 = Enu + sqrt( (Enu-mp-Up+En)**2 - me**2 )
  q3_max = min(tmp1, tmp2)
  call q3_sol(pn3,q3_lim)
  q3_min = q3_lim(1)
  q3_max = q3_lim(2)
  q3 = (q3_max-q3_min)*x(2) + q3_min

  tmp1 = Enu - sqrt( (Enu-q3)**2 + me**2 )! Lp_max
  tmp2 = Enu - sqrt( (Enu+q3)**2 + me**2 )! Lp_min
  tmp3 = sqrt( (pn3+q3)**2+mp**2) + Up - En	 ! hd_max
  tmp4 = sqrt( (pn3-q3)**2+mp**2) + Up - En	 ! hd_min
  q0_max = min(tmp1,tmp3)
  q0_min = max(tmp2,tmp4)

  tmp = mup-50.0*Tem-En
  if( tmp < q0_max ) then
     q0_min = max(q0_min, tmp)
  end if
  tmp = Enu-max(me,mue-50.0*Tem)
  if( q0_min<tmp ) then
     q0_max = min(q0_max, tmp)
  end if

  cthe_max = min(1.d0, ((q0_max+En-Up)**2-mp**2-pn3**2-q3**2)&
	    /(2.d0*pn3*q3) )
  cthe_min = max(-1.d0, ((q0_min+En-Up)**2-mp**2-pn3**2-q3**2)&
	   /(2.d0*pn3*q3) )
  cthe_n = (cthe_max - cthe_min)*x(3) + cthe_min
  phi_n = 2.d0*pi*x(4)

  pp3 = sqrt(pn3**2 + q3**2 + 2.d0*cthe_n*pn3*q3)
  Ep = sqrt(pp3**2 + mp**2) + Up
  q0 = Ep - En
  Ee = Enu - q0
  pe3 = sqrt(Ee**2 - me**2)

  cthe_v = (Enu**2+q3**2-pe3**2)/(2.d0*Enu*q3)
  cthe = (Enu**2+pe3**2-q3**2)/(2.d0*Enu*pe3)

  q(0:3) = (/q0-Up+Un,0.d0,0.d0,q3/)
  pv(0:3) = (/Enu,Enu*sqrt(1.d0-cthe_v**2),0.d0,Enu*cthe_v/)
  pe(0:3) = pv(0:3) - q(0:3)
  pe(0) = Enu - q0

  pn(0:3) = (/En-Un,pn3*sqrt(1.d0-cthe_n**2)*cos(phi_n),&
	   pn3*sqrt(1.d0-cthe_n**2)*sin(phi_n),pn3*cthe_n/)
  pp(0:3) = pn(0:3) + q(0:3)

  Jacob = (cthe_max-cthe_min)*(En_max-En_min)&
	 *(q3_max-q3_min)*(2.d0*pi)

12  continue

end subroutine initmomentum


subroutine nubp_4D()
  use global_4D
  implicit none
  real*8:: Mtmp, Utmp, mutmp
!  1) exchange mn <-> mp; and Un <-> Up, mun <-> mup
  Mtmp = Mn
  Mn = Mp
  Mp = Mtmp
  Utmp = Un
  Un = Up
  Up = Utmp
  mutmp = mun
  mun = mup
  mup = mutmp
!!  2) exchange mue -> -mue
  mue = -mue
!!  times -1 for VA & VF terms
  jVA = -1.d0
  jAF = -1.d0
end subroutine nubp_4D


subroutine Range_pn_4D()
  use global_4D
  implicit none
  real*8:: tmp,xpn3,mpt,Epr,Esq,Equ,p2a,p2b,p20,F0,Finf,Fmax,Fmin,pnmax0
  integer:: i

  pnmin = 0.
  pnmax = sqrt( (Tfac*Tem + mun - Un)**2 - mn**2 )
  pnmax0 = pnmax

  mpt = mp + me
  F0 = sqrt(Enu**2 + mpt**2) + Up - mn - Un
  Finf = -Enu + Up - Un
  if(mpt < mn) then
     tmp = (1.d0+mpt/mn)/(1.d0-mpt**2/mn**2)*Enu
     Fmin = sqrt( (tmp-Enu)**2+mpt**2 )+Up-sqrt(tmp**2+mn**2)-Un
     Fmax = max( F0, Finf )
  else
     Fmin = Finf
     Fmax = F0
  end if

  if( Fmax .le. Enu) goto 12
  if( Fmin .ge. Enu) then
     pnmax = -1000.0
     goto 12
  end if

  Epr = Enu - Up + Un
  Esq = Enu**2 + mpt**2 - mn**2 - Epr**2
  Equ = Esq**2 + 4.d0*mn**2*( Enu**2-Epr**2 )
  if(Equ.gt.0.d0.and. Enu**2.ne.Epr**2) then
     p2a = ( Enu*Esq - abs(Epr)*sqrt(Equ) )*0.5d0/(Enu**2-Epr**2)
     p2b = ( Enu*Esq + abs(Epr)*sqrt(Equ) )*0.5d0/(Enu**2-Epr**2)
     if(Epr<0.d0 .and. p2b>0.d0) then
	pnmin=p2b
     else if(Epr.ge.0.d0 .and. p2a>0.d0) then
	pnmin=p2a
     end if

  else if(Enu**2 .eq. Epr**2) then
     p20 = -(4.0*Epr**2*mn**2-Esq**2)/(4.0*Esq*Enu)
     if(p20>0.d0) pnmin=p20
  end if

12  continue

end subroutine Range_pn_4D


subroutine q3_sol(xPn3, q3_lim)
  use global_4D
  implicit none
  real*8:: xPn3,q3_lim(2),tmp1,tmp2,xA,xB,xC
  q3_lim = 0.
  tmp1 = Up-En-Enu
  tmp2 = Enu**2 + me**2 - mp**2 - xpn3**2
  xA = 4.d0*( tmp1**2-(Enu-xpn3)**2 )
  xB = 4.d0*( (tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*xpn3 )
  xC = 4.d0*tmp1**2*(xpn3**2+mp**2) - (tmp2-tmp1**2)**2

  q3_lim(1) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  q3_lim(2) = (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA)

  return
end subroutine q3_sol



integer function integrand_D(ndim, xx, ncomp, ff)
  use global_4D
  implicit none
  integer:: ndim,ncomp
  real*8:: xx(*),ff(*)
  real*8, external :: amp_4D
  real*8:: xf2,xf3,xf4
  call initmomentum_D(xx)
  ff(1)=amp_4D()*q3*pn3/(16.d0*(2.d0*pi)**4*Enu**2*pp(0))*Jacob
  xf2 = 1.0/( exp((pn(0)+Un-mun)/Tem) + 1.0 )
  xf3 = 1.0/( exp((pe(0)-mue)/Tem) + 1.0 )
  xf4 = 1.0/( exp((pp(0)+Up-mup)/Tem) + 1.0 )
  if(opt0.eq.3) then
     ff(1) = ff(1)*xf2*(1.d0-xf4)*xf3
  else if(opt0.eq.4) then
     ff(1) = ff(1)*(1.0-xf2)*xf4*(1.0-xf3)
  end if
  integrand_D = 0
end function integrand_D


subroutine initmomentum_D(x)
  use global_4D
  implicit none
  real*8:: x(4),En_min,En_max,tmp,tmp1,q3_max,q3_min,&
   cthe_max,cthe_min,cthe_n,Ep,pp3,q0,Ee,Pe3,tmp2,q3_lim(2),&
   q0_max,q0_min,tmp3,tmp4,cthe_v,phi_n
  real*8:: xq, lp_min,hd_max,hd_min
  integer:: i
  real*8:: xk,kmin,xAA,xQQ,xE2_min,xE2_max

  En_min = sqrt(pnmin**2+mn**2)+Un
  En_max = sqrt(pnmax**2+mn**2)+Un
  En = (En_max-En_min)*x(1) + En_min
  pn3 = sqrt( (En-Un)**2 - mn**2 )

  call q3_sol_D(pn3,q3_lim)
  q3_min = q3_lim(1)
  q3_max = q3_lim(2)
  q3 = (q3_max-q3_min)*x(2) + q3_min

  tmp1 = Enu + sqrt( (Enu+q3)**2 + me**2 )
  tmp2 = Enu + sqrt( (Enu-q3)**2 + me**2 )
  tmp3 = sqrt( (pn3+q3)**2+mp**2) + Up - En
  tmp4 = sqrt( (pn3-q3)**2+mp**2) + Up - En
  q0_max = min(tmp1,tmp3)
  q0_min = max(tmp2,tmp4)

  cthe_max = min(1.d0, ((q0_max+En-Up)**2-mp**2-pn3**2-q3**2)&
	    /(2.d0*pn3*q3))
  cthe_min = max(-1.d0, ((q0_min+En-Up)**2-mp**2-pn3**2-q3**2)&
	    /(2.d0*pn3*q3))


  cthe_n = (cthe_max - cthe_min)*x(3) + cthe_min
  phi_n = 2.d0*pi*x(4)

  pp3 = sqrt(pn3**2 + q3**2 + 2.d0*cthe_n*pn3*q3)
  Ep = sqrt(pp3**2 + mp**2) + Up
  q0 = Ep - En
  Ee = q0 - Enu
  pe3 = sqrt(Ee**2 - me**2)

  cthe_v = (Enu**2+q3**2-pe3**2)/(2.d0*Enu*q3)
  cthe = -(Enu**2+pe3**2-q3**2)/(2.d0*Enu*pe3)

  q(0:3) = (/q0-Up+Un,0.d0,0.d0,q3/)
  pv(0:3) = (/Enu,Enu*sqrt(1.d0-cthe_v**2),0.d0,Enu*cthe_v/)
  pe(0:3) = q(0:3) - pv(0:3)
  pe(0) = q0 - Enu

  pn(0:3) = (/En-Un,pn3*sqrt(1.d0-cthe_n**2)*cos(phi_n),&
      pn3*sqrt(1.d0-cthe_n**2)*sin(phi_n),pn3*cthe_n/)
  pp(0:3) = pn(0:3) + q(0:3)

  Jacob = (cthe_max-cthe_min)*(En_max-En_min)&
      *(q3_max-q3_min)*(2.d0*pi)
12  continue
end subroutine initmomentum_D

subroutine Range_pn_dec_4D()
  use global_4D
  implicit none
  real*8:: tmp,xpn3,mpt,Epr,Esq,Equ,p2a,p2b,p20
  real*8:: F0, xpn, xpn1, xF, Fmax
  real*8:: Lmax,Lmin,Hmin,Hmax,xq

  pnmin = 0.
  ppmax = sqrt( (Tfac*Tem + mup - Up)**2 - mp**2 )
  select case(opt0)
  case (3)
  pnmax = sqrt( (Tfac*Tem + mun - Un)**2 - mn**2 )
  case (4)
  pnmin = sqrt( (max(mn,mun-Tfac*Tem-Un))**2 - mn**2 )
  pnmax = sqrt( (sqrt(ppmax**2+mp**2)+Up-Un)**2-mn**2 )
  end select

  mpt = mp - me
  if( mpt .gt. mn) then
     tmp = mn/(mpt-mn)*Enu
  else
     tmp = pnmax
  end if
  Fmax = sqrt( (tmp+Enu)**2+mpt**2 )+Up-sqrt(tmp**2+mn**2)-Un

  if(Enu>=Fmax) then
     pnmax = -1000.
     goto 12
  end if

  F0 = sqrt(Enu**2 + mpt**2)+Up-mn-Un
  if( F0 .ge. Enu .and. Up .ge. Un) goto 12

  Epr = Enu - Up + Un
  Esq = Enu**2 + mpt**2 - mn**2 - Epr**2
  Equ = Esq**2 + 4.d0*( (mn*Enu)**2 - (Epr*mn)**2 )
  if(Equ.gt.0.d0.and. Enu**2.ne.Epr**2) then
     p2a = -( Enu*Esq - abs(Epr)*sqrt(Equ) )*0.5d0/(Enu**2-Epr**2)
     p2b = -( Enu*Esq + abs(Epr)*sqrt(Equ) )*0.5d0/(Enu**2-Epr**2)
     if(Epr<0.d0 .and. p2b>0.d0) then
	pnmin=p2b
     else if(Epr.ge.0.d0 .and. p2a>0.d0) then
	pnmin=p2a
     end if

  else if(Enu**2 .eq. Epr**2) then
     p20 = (4.0*Epr**2*mn**2-Esq**2)/(4.0*Esq*Enu)
     if(p20>0.d0) pnmin=p20
  end if
12  continue
end subroutine Range_pn_dec_4D


subroutine q3_sol_D(xPn3, q3_lim)
  use global_4D
  implicit none
  real*8:: xPn3,q3_lim(1:2),tmp1,tmp2,xA,xB,xC,xqmax,q3_lim2(1:2)
  real*8:: tmp3,tmp4,tmp5,tmp6

  q3_lim = 0.0
  q3_lim2 = 0.0

  tmp1 = Up-En-Enu
  tmp2 = Enu**2 + me**2 - mp**2 - xpn3**2
  xA = 4.d0*( tmp1**2-(Enu+xpn3)**2 )
  xB = 4.d0*( (tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*xpn3 )
  xC = 4.d0*tmp1**2*(xpn3**2+mp**2) - (tmp2-tmp1**2)**2
  q3_lim(1) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  q3_lim(2) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )

  tmp3 = Enu + sqrt( (Enu-q3_lim(1))**2 + me**2 ) &
    -(sqrt( (xpn3+q3_lim(1))**2+mp**2) + Up - En)
  tmp4 = Enu + sqrt( (Enu-q3_lim(2))**2 + me**2 ) &
    -(sqrt( (xpn3+q3_lim(2))**2+mp**2) + Up - En)
  tmp5 = Enu + sqrt( (Enu+q3_lim(1))**2 + me**2 ) &
    -(sqrt( (xpn3-q3_lim(1))**2+mp**2) + Up - En)
  tmp6 = Enu + sqrt( (Enu+q3_lim(2))**2 + me**2 ) &
    -(sqrt( (xpn3-q3_lim(2))**2+mp**2) + Up - En)

  if(abs(tmp3)>1.d-10 .and. abs(tmp5)>1.d-10) then
   q3_lim(1)=0.0
  end if
  if(abs(tmp4)>1.d-10 .and. abs(tmp6)>1.d-10) then
   q3_lim(2)=max( q3_lim(1) + 20.*Tem, Enu + (20.0*Tem + mue) )
  end if

  return
end subroutine q3_sol_D



function amp_4D()
  use global_4D
  implicit none
  real*8:: Mvv,Maa,Mff,Mva,Mvf,Maf,Mpp,Map,tmp1,tmp2,tmp3,&
     d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d55,amp_4D,gv,ga,f2,gp
  real*8,external:: cdot4
  real*8:: gasq,gvsq,gaf,gvf,gva,f2sq,gpsq,gap,pipole

  d1 = cdot4(pn,pv)
  d2 = cdot4(pp,pe)
  d3 = cdot4(pn,pe)
  d4 = cdot4(pp,pv)
  d5 = cdot4(pv,pe)
  d55 = cdot4(pn,pp)
  tmp1 = d1*d2
  tmp2 = d3*d4
  tmp3 = mn*mp*d5
  d6 = cdot4(pv,q)
  d7 = cdot4(pe,q)
  d8 = cdot4(pn,q)
  d9 = cdot4(pp,q)
  d10 = cdot4(q,q)

  gv = 1.0
  f2 = 3.706
  ga = ga0
  gp = 0.0 
  if(opt_form .eq. 1) then
     gv=( 1.d0-4.706d0*d10/(4.d0*Mnp**2) )/( 1.d0 - d10/(4.d0*Mnp**2))&
	/( 1.d0 - d10/MV**2 )**2
     f2 = 3.706d0/( 1.d0 - d10/(4.d0*Mnp**2) )/( 1.d0 - d10/MV**2 )**2
     ga = gA0/( 1.d0 - d10/MA**2 )**2
  end if
  if(opt_pseudo.eq.1) gp=2.d0*Mnp**2*ga/(mpi**2 - d10)  

  Mvv = 16.d0*(Gf*Vud*gv)**2*( tmp1 + tmp2 - tmp3 )
  Maa = 16.d0*(Gf*Vud*ga)**2*( tmp1 + tmp2 + tmp3 )
  Mva = jVA*32.d0*(Gf*Vud)**2*gv*ga*( tmp1 - tmp2 )
  Mvf = 8.d0*(Gf*Vud)**2*gv*f2/Mnp*(   &
    (d1*mp - d4*mn)*d7 + (d3*mp-d2*mn)*d6 + (d8*mp - d9*mn)*d5 )
  Maf = jAF*16.d0*(Gf*Vud)**2*ga*f2/Mnp*((d1*mp + d4*mn)*d7 &
       -(d3*mp+d2*mn)*d6)
  Mff = 2.d0*(Gf*Vud*f2/Mnp)**2*(      &
    d10*(d55*d5-2.d0*d1*d2-2.d0*d3*d4) &
    + 2.d0*d7*(d1*d9 + d4*d8 - d55*d6) + 2.d0*d6*(d3*d9 + d2*d8) &
    - mn*mp*(d5*d10 + 2.d0*d6*d7) )

  Mpp = 8.d0*(Gf*Gp*Vud/Mnp)**2*(2.d0*d6*d7-d10*d5)*(d55-mn*mp)
  Map = 16.d0*(Gf*Vud)**2*gA*gP/mnp*(  &
    d6*(d3*mp-d2*mn) + d7*(d1*mp-d4*mn) + d5*(d9*mn-d8*mp) )

  amp_4D = Mvv + Maa + Mva + opt_wm*(Mvf + Maf + Mff) &
          + opt_pseudo*(Mpp + Map) 

  return
end function amp_4D



function cdot4(xp1,xp2)
  implicit none
  real*8:: xp1(0:3),xp2(0:3),cdot4
  cdot4 = xp1(0)*xp2(0)-xp1(1)*xp2(1)-xp1(2)*xp2(2)-xp1(3)*xp2(3)
  return
end function cdot4

function cdot3(xp1,xp2)
  implicit none
  real*8:: xp1(0:3),xp2(0:3),cdot3
  cdot3 = xp1(1)*xp2(1)+xp1(2)*xp2(2)+xp1(3)*xp2(3)
  return
end function cdot3

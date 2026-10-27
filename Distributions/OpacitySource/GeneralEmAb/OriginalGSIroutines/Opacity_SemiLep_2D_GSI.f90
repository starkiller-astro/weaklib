!!! modified date: June 17, 2020
!!! Authors: Codes wirtten by Gang Guo, based on the formalism by
!!! Andreas Lohs & Gang Guo
!!! neutrino CC opacity is given by 2D integrals
!!! processes to be considered in this routine
!!! opt=1: v + n -> e- + p  (oapcity in 1/km)
!!! opt=2: vb + p -> e+ + n (opacity in 1/km)
!!! opt=3: vb + p + e- -> n (opacity in 1/km)
!!! opt=4: n -> vb + e- + p [emissivity in v/(s cm^3 MeV)]

!!! Ngrids defined in module: the No. of grids used in gaussian
!!! quadrature for each dimention
!!! seems Ngrids needs be >~30 to reach an accuracy of 5%

module CC_global
  implicit none
  real*8:: me,mn,mp,Un,Up
  real*8,parameter :: pi=3.1415927d0,Gf=1.166d-11,Vud=0.97427d0,gA0=&
       1.2723d0,gV0=1.d0,f2wm0=3.706d0
  real*8::Tem,mun,mup,mue
  real*8:: Enu,En,pnmax,pnmin,ppmax,pemax,pemin,Tfac=150.0
  integer:: anti,opt0,opt_form=0,opt_pseudo=0,opt_WM=0
  real*8:: massA=1.d3,massV=840.0,Mpi=139.57
  integer,parameter:: Ngrids=50      !!! to be tuned by user     
end module CC_global

!========================================================================
! Routine to calculate the neutrino CC opacity via 2D integrals
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
subroutine Opacity_CC_2D_GSI(opt_had,opt,NP,EnuA,OpaA,xTem,cheme,chemn,chemp,&
     masse,massn,massp,xUn,xUp) 
  use CC_global
  implicit none
  integer:: NP,i,opt,opt_had
  real*8:: EnuA(NP),OpaA(NP),xTem,cheme,chemn,chemp,phin,phip,masse,&
       massn,massp,xUn,xUp 

  ! CAREFUL. IN THIS ROUTINE A LOT OF NUMBERS ARE DEFINED AS REAL(4) INSTEAD OF REAL(8)
  ! E.G. 1.0 INSTEAD OF 1.0d0
       
  opt0 = opt

  select case(opt_had)
  case (0)
     opt_wm=0
     opt_pseudo=0
     opt_form=0
  case (1)
     opt_wm=1
     opt_pseudo=0
     opt_form=0
  case (2)
     opt_wm=1
     opt_pseudo=1
     opt_form=0
  case (3)
     opt_wm=1
     opt_pseudo=1
     opt_form=1
  end select

  me = masse

  Tem = xTem
  mue = cheme
  if(opt0.eq.2)then
     mue = -cheme
  end if
  select case(opt0)
  case (1)
     mn=massn
     mp=massp
     Un=xUn
     Up=xUp
     mun = chemn
     mup = chemp
     anti = 1
  case (2,3,4)
     mp=massn
     mn=massp
     Up=xUn
     Un=xUp
     mup = chemn
     mun = chemp
     anti = -1
  end select

  select case(opt0)
  case (1,2)
     do i=1,NP
	call Integral_2D(EnuA(i),OpaA(i))
     end do
  case (3,4)
     do i=1,NP
	call Integral_2D_D(EnuA(i),OpaA(i))
     end do
  end select

end subroutine Opacity_CC_2D_GSI

subroutine Integral_2D(xEnu,res)   ! captures
  use CC_global
  implicit none
  real*8:: xpn3,xq0_min,xq0_max,Ee_min,Ee_max,xEe,En_min,En_max
  real*8:: xa(Ngrids),wxa(Ngrids),wEn(Ngrids),wEe(Ngrids)
  integer:: i, j
  real*8:: xEnu,res,xamp,Ep
  real*8:: xf2,xf3,xf4

  Enu = xEnu
  res = 0.0
  call Range_pn()
  if( pnmin > pnmax ) goto 12
  En_min = sqrt(pnmin**2+mn**2) + Un
  En_max = sqrt(pnmax**2+mn**2) + Un
  call gauleg(0.d0,1.d0,xa,wxa,Ngrids)

  do i=1,Ngrids
     En = xa(i)*(En_max-En_min)+En_min
     wEn(i) = wxa(i)*(En_max-En_min)
     xpn3 = sqrt( (En-Un)**2-mn**2 )
     call Ebounds( xpn3, xq0_min, xq0_max)
     Ee_min = Enu - xq0_max
     Ee_max = Enu - xq0_min
     do j=1,Ngrids
	xEe = xa(j)*(Ee_max-Ee_min)+Ee_min
	wEe(j) = wxa(j)*(Ee_max-Ee_min)
	Ep = Enu+En-xEe
	call Calc_ampsq(En,xEe,xamp)
	xf2 = 1.0/( exp((En-mun)/Tem) + 1.0 )
	xf3 = 1.0/( exp((xEe-mue)/Tem) + 1.0 )
	xf4 = 1.0/( exp((Ep-mup)/Tem) + 1.0 )
	res = res + xf2*(1.0-xf3)*(1.0-xf4)*xamp*wEn(i)*wEe(j)/Enu**2
     end do
  end do
  res = res*(gf*Vud)**2/16.0/(pi**5)/197.327d-18
12 continue
end subroutine Integral_2D


subroutine Integral_2D_D(xEnu,res)   ! decay/inverse decay
  use CC_global
  implicit none
  real*8:: xpn3,xq0_min,xq0_max,Ee_min,Ee_max,xEe,En_min,En_max
  real*8:: xa(Ngrids),wxa(Ngrids),wEn(Ngrids),wEe(Ngrids)
  integer:: i,j
  real*8:: xEnu,res,xamp
  real*8:: E1,E2,E3,E4,E2f,E4f,mu2,mu4,mu3,U2,U4,Ep
  real*8:: xf2,xf3,xf4

  U2 = Un
  U4 = Up
  mu2 = mun
  mu4 = mup
  mu3 = mue

  Enu = xEnu
  E1 = xEnu
  res = 0.0
  call Range_pn_D()
  if( pnmin > pnmax ) goto 12
  En_min = sqrt(pnmin**2+mn**2) + Un
  En_max = sqrt(pnmax**2+mn**2) + Un

  call gauleg(0.d0,1.d0,xa,wxa,Ngrids)
  do i=1,Ngrids
     En = xa(i)*(En_max-En_min)+En_min
     wEn(i) = wxa(i)*(En_max-En_min)
     xpn3 = sqrt( (En-Un)**2-mn**2 )
     call Ebounds_D( xpn3, xq0_min, xq0_max )
     Ee_min = xq0_min - Enu
     Ee_max = xq0_max - Enu
     do j=1,Ngrids
	xEe = xa(j)*(Ee_max-Ee_min)+Ee_min
	wEe(j) = wxa(j)*(Ee_max-Ee_min)
	Ep = Enu+xEe+En
	call Calc_ampsq(En,-xEe,xamp)
	xf2 = 1.0/( exp((En-mun)/Tem) + 1.0 )
	xf3 = 1.0/( exp((xEe-mue)/Tem) + 1.0 )
	xf4 = 1.0/( exp((Ep-mup)/Tem) + 1.0 )
	if(opt0.eq.3) then	    ! inverse decay
	   res = res + xf2*xf3*(1.0-xf4)*xamp*wEn(i)*wEe(j)/Enu**2
	else if(opt0.eq.4) then	    ! n decay, *blocking of neutrinos is not added*
	   res = res + (1.0-xf2)*(1.0-xf3)*xf4*xamp*wEn(i)*wEe(j)/Enu**2
	end if
     end do
  end do

  if(opt0.eq.3) then
     res = res*(gf*Vud)**2/16.0/(pi**5)/197.327d-18
  else if(opt0.eq.4) then
     res = res*(gf*Vud)**2/32.0/(pi**7)*Enu**2*(1.d13/197.327)**3*&
          (3.d23/197.327) 
  end if

  res = -res

12 continue

end subroutine Integral_2D_D


!!! calculate the amplitudes, keep the following unchanged
subroutine Calc_Ampsq( E2, E3, ampsq )
  use CC_global
  implicit none
  real*8:: Ampsq, P1,P2,P3,P4,Pmax,Pmin,Pmax2,Pmin2,Pmax3,Pmin3
  real*8:: Ia,Ib,Ic,Id,Ie,Ig,Ih,Ij,Ik,Il,A,B,Ct,D,E,G,H,J,K,L,del0,&
       del1,del2,eps0,eps1,alp0,alp1,alp2,bet0,bet1,cplA,cplB,cplC,&
       cplD,cplE,cplG,cplH,cplJ,cplK,cplL,JFF,KFF,LFF 
  real*8:: m2f,E2f,m4f,E4f,E1,E3,Qmass,dmf,U2,U4,dU2,m3,E2,E4,m2,xE3
  real*8:: teps0, teps1
  real*8:: Iam,Iap,Ibm,Ibp,Icm,Icp,Idm,Idp,Iem,Iep,Igm,Igp,Ikm,Ikp
  real*8:: gasq=0.0,gvsq=0.0,gva=0.0,gaf=0.0,gvf=0.0,f2sq=0.0
  real*8:: AmassSq,VmassSq,AVmassSq,APmassSq,AFmassSq,VFmassSq,&
       FmassSq,tm3sq,tbet0 
  real*8:: xbet0,txi1,Ia2,tgam0,Ib2,txi3,Ic2,Id2,Ie2,Ig2,Ih2,Ij2,Ik2,Il2,&
       cplA2,cplB2,cplC2,cplD2,cplE2,cplG2,cplH2,cplJ2,cplK2,cplL2
  real*8::YAP,ZAP,MassMy,Inte1,Inte2,IlAP,IlPP,IkAP,tIkAP,IkPP,tIkPP,&
       IjAP,tIjAP,IGPP,tIgPP,TanArg,LogArg,cplLAP,cplKAP,cplJAP,&
       cplLPP,cplKPP,cplGPP 
  real*8:: ampsq_ff0,ampsq_pseudo       

  Ampsq = 0.0

  m2 = 0.5*(939.565+938.272)
  m2f = mn
  m4f = mp
  m3 = me
  U2 = Un
  U4 = Up
  dU2 = U2 - U4
  dmf = m2f - m4f
  Qmass = 0.5*( m2f**2 - m4f**2 )

  E1 = Enu
  E4 = E1 + E2 - E3
  E4f = E4 - U4
  E2f = E2 - U2
  P1 = E1
  P3 = sqrt( E3**2 - m3**2 )
  P2 = sqrt( (E2 - U2)**2 - m2f**2 )
  P4 = sqrt( (E4 - U4)**2 - m4f**2 )

  if( P1 .ne. P1) return
  if( P2 .ne. P2) return
  if( P3 .ne. P3) return
  if( P4 .ne. P4) return

  if ((P1+P2)>(P3+P4)) then
     Pmax = P3+P4
  else
     Pmax = P1+P2
  endif
  if ((abs(P1-P2))>(abs(P3-P4))) then
     Pmin = abs(P1-P2)
  else
     Pmin = abs(P3-P4)
  endif

  A = E1*E2f+0.5*(P1**2+P2**2)
  B = E3*E4f+0.5*(P3**2+P4**2)
  Ia = pi**2/15.0*(3.0*(Pmax**5-Pmin**5)-10.0*(A+B)*(Pmax**3-Pmin**3)+&
       60.0*A*B*(Pmax-Pmin)) 

  if ((P1+P4)>(P3+P2)) then
     Pmax2 = P3+P2
  else
     Pmax2 = P1+P4
  endif
  if ((abs(P1-P4))>(abs(P3-P2))) then
     Pmin2 = abs(P1-P4)
  else
     Pmin2 = abs(P3-P2)
  endif

  Ct = E1*E4f-0.5*(P1**2+P4**2)
  D = E2f*E3-0.5*(P2**2+P3**2)
  Ib = pi**2/15.0*(3.0*(Pmax2**5-Pmin2**5)+10.0*(Ct+D)*(Pmax2**3-&
       Pmin2**3)+60.0*Ct*D*(Pmax2-Pmin2)) 

  del0 = -(P1**2-P2**2)*(P3**2-P4**2)/4.0
  del1 = E1*E3 + (-P1**2+P2**2-P3**2+P4**2)/4.0
  del2 = -0.25
  eps0 = E1*E2f + (P1**2+P2**2)/2.0
  eps1 = -0.5
  Ic =	pi**2*4.0*(del2*eps1**2/7.0*(Pmax**7-Pmin**7)+(2.0*del2*eps0*&
       eps1+del1*eps1**2)/5.0*(Pmax**5-Pmin**5)+(del2*eps0**2+2.0*&
       del1*eps0*eps1+del0*eps1**2)/3.0*(Pmax**3-Pmin**3)+(del1*eps0**&
       2+2.0*del0*eps0*eps1)*(Pmax-Pmin)-del0*eps0**2*(1.0/Pmax-1.0/&
       Pmin)) 

  if ((P1+P3)>(P2+P4)) then
     Pmax3 = P2+P4
  else
     Pmax3 = P1+P3
  endif
  if ((abs(P1-P3))>(abs(P2-P4))) then
     Pmin3 = abs(P1-P3)
  else
     Pmin3 = abs(P2-P4)
  endif

  alp0 = (P1**2-P3**2)*(P2**2-P4**2)/4.0
  alp1 = E1*E2f + (P1**2+P2**2-P3**2-P4**2)/4.0
  alp2 = 0.25
  bet0 = E1*E3 - (P1**2+P3**2)/2.0
  bet1 = 0.5
  Id = pi**2*4.0*(alp2*bet1**2/7.0*(Pmax3**7-Pmin3**7)+(2.0*alp2*bet0*&
       bet1+alp1*bet1**2)/5.0*(Pmax3**5-Pmin3**5)+& 
       (alp2*bet0**2+2.0*alp1*bet0*bet1+alp0*bet1**2)/3.0*(Pmax3**3-&
       Pmin3**3)+(alp1*bet0**2+2.0*alp0*bet0*bet1)*& 
       (Pmax3-Pmin3)-alp0*bet0**2*(1.0/Pmax3-1.0/Pmin3))

  Ie = pi**2/15.0*(3.0*(Pmax**5-Pmin**5)-20.0*A*(Pmax**3-Pmin**3)+&
       60.0*A**2*(Pmax-Pmin)) 

  E = E1*E3 - 0.5*(P1**2+P3**2)
  Ig = pi**2/15.0*(3.0*(Pmax3**5-Pmin3**5)+20.0*E*(Pmax3**3-Pmin3**3)+&
       60.0*E**2*(Pmax3-Pmin3)) 

  Ih = pi**2*4.0*(alp2*bet1/5.0*(Pmax3**5-Pmin3**5)+(alp2*bet0+alp1*&
       bet1)/3.0*(Pmax3**3-Pmin3**3)+(alp1*bet0+alp0*bet1)& 
       *(Pmax3-Pmin3)-alp0*bet0*(1.0/Pmax3-1.0/Pmin3))

  Ij = pi**2/15.0*(-10.0*(Pmax**3-Pmin**3)+60.0*A*(Pmax-Pmin))

  Ik = pi**2/15.0*(10.0*(Pmax3**3-Pmin3**3)+60.0*E*(Pmax3-Pmin3))

  Il = pi**2/15.0*(60.0*(Pmax-Pmin))

  gasq = ga0**2
  gvsq = gv0**2
  gva = gv0*ga0
  if(opt_WM .eq. 1) then
     gvf = gv0*F2wm0
     gaf = ga0*F2wm0
     f2sq = f2Wm0**2
  else
     gvf = 0.0
     gaf = 0.0
     f2sq = 0.0
  end if

  xE3 = E3
  cplA = (gvsq+2.*anti*gva+gasq)+anti*2.0*gaf*m2f/m2*(1.0-dmf/2.0/m2f)
  cplB = (gvsq-2.*anti*gva+gasq)-anti*2.0*gaf*m2f/m2*(1.0-dmf/2.0/m2f)
  cplC = F2sq/m2**2
  cplD =-F2sq/m2**2
  cplE = F2sq/m2**2*(-0.5*m3**2+dU2*(xE3-E1)-0.5*dU2**2)
  cplG = gvf*m2f/m2*(2.0-dmf/m2f)+0.5*F2sq/m2**2*(m2f*m4f-Qmass+0.25*&
       m3**2-dU2*(E1+E2f)-dU2**2/4.0) 
  cplH = 0.5*F2sq/m2**2*(2.0*Qmass+m3**2+dU2*(3.0*E1-xE3+2.0*E4f))

  JFF = -m3**2*(E1+0.5*E2f+0.5*E4f)+Qmass*(xE3-3.0*E1)		  &
       +0.5*dU2*(E4f*(3.*xE3-5.*E1)+E2f*(xE3+E1)+xE3**2-E1**2-2.*Qmass)	 &
       + dU2**2*(E1-xE3-E4f)+0.5*dU2**3
  JFF = JFF*dU2
  cplJ = gvf*dmf/m2*0.5*(m3**2-dU2*(E1+xE3))+F2sq/m2**2*0.5*JFF

  kFF = -(m2f+3.*m4f)*m2f*0.25*m3**2 + Qmass**2 + Qmass*0.25*m3**2-m3**4/8.0 &
       +dU2*(0.5*Qmass*(3.*E1-E2f+xE3+3.*E4f)+0.25*m3**2*(2.*E2f+xE3+&
       E1)+m2f*m4f*(xE3-E1)) & 
       +dU2**2*( 0.25*(m2f**2-m2f*m4f-3.*Qmass)+E4f*0.5*(2.*E1-E2f+&
       xE3+E4f)+E2f*(0.5*E1-xE3)+0.5*E1**2)	 & 
       +dU2**3*0.25*( -E1+2.*E2f-xE3-2.*E4f ) + 0.125*dU2**4
  cplK = (gasq-gvsq)*m2f*m4f+gvf*m2f/m2*0.5*(-3.0*m3**2+4.0*dU2*(xE3-&
       E1)-dU2**2 & 
       +dmf/m2f*(2.0*Qmass+m3**2+dU2*(E4f+2.0*E1-xE3)))&
       +0.5*F2sq/m2**2*KFF

  LFF = m3**2*(m2f+m4f)**2-4.*Qmass**2+dU2*(-m3**2*(E2f+E4f+E1)-2.*&
       xE3*(m2f**2+m2f*m4f) & 
       +2.*Qmass*(E2f-3.*E4f-E1)) &
       +2.*dU2**2*( E2f*xE3+E4f*(E2f-E4f-E1) + Qmass ) + dU2**3*(-E2f+E4f+E1)
  LFF = LFF*dU2*E1*0.25
  cplL = gvf*m2f/m2*dU2*E1*(m3**2-dU2*xE3-0.5*dmf/m2f*(Qmass+0.5*m3**&
       2+dU2*E4f-0.5*dU2**2))   & 
       +0.5*F2sq/m2**2*LFF

  Ampsq = max(0.0,cplA*Ia+cplB*Ib+cplC*Ic+cplD*Id+cplE*Ie+cplG*Ig+cplH*Ih+&
       cplJ*Ij+cplK*Ik+cplL*Il) 

!!! adding form-factor dependences
  AmassSq = massA**2
  APmassSq = 2.0*AmassSq
  tm3sq = m3**2 + 2.0*dU2*(E1-E3) + dU2**2
  if(opt_form .eq. 1) then
     VmassSq = 1.0/(1.0/massV**2-3.706/m2**2/8.0)
     AVmassSq = 1.0/(0.5/massA**2+0.5/massV**2-3.706/m2**2/16.0)
     AFmassSq = 1.0/(0.5/massA**2+0.5/massV**2+1.0/m2**2/16.0)
     VFmassSq = 1.0/(1.0/massV**2-(3.706-1.0)/16.0/m2**2)
     FmassSq = 1.0/(1.0/massV**2+1.0/8.0/m2**2)

     tbet0 = (P1**2-P2**2)*(P4**2-P3**2)/4.0
     txi1 = 0.5*tm3sq-E1*E3+(P1**2-P2**2+P3**2-P4**2)/4.0

     Iap = pi**2/15.0*(5./7.*3.0*(Pmax**7-Pmin**7)-3.0/5.0*10.0*(A+B)*&
          (Pmax**5-Pmin**5)+60.0/3.0*A*B*(Pmax**3-Pmin**3)) 
     Iam = pi**2/15.0*(5./3.*3.0*(Pmax**3-Pmin**3)-3.0*10.0*(A+B)*&
          (Pmax-Pmin)-60.0*A*B*(1.0/Pmax-1.0/Pmin)) 
     Ia2 = 4.0*( 0.5*Iap + 2.*txi1*Ia - 2.0*tbet0*Iam )

     tgam0 = (P1**2-P4**2)*(P2**2-P3**2)/4.0

     Ibp = pi**2/15.0*(5.0/7.0*3.0*(Pmax2**7-Pmin2**7)+3.0/5.0*10.0*&
          (Ct+D)*(Pmax2**5-Pmin2**5)+60.0/3.0*Ct*D*(Pmax2**3-Pmin2**&
          3)) 
     Ibm = pi**2/15.0*(5.0/3.0*3.0*(Pmax2**3-Pmin2**3)+3.0*10.0*(Ct+&
          D)*(Pmax2-Pmin2)-60.0*Ct*D*(1.0/Pmax2-1.0/Pmin2)) 
     Ib2 = 4.0*( 0.5*Ibp + 2.*txi1*Ib - 2.0*tgam0*Ibm )

     Icp = pi**2*4.0*(del2*eps1**2/9.0*(Pmax**9-Pmin**9)+(2.0*del2*&
          eps0*eps1+del1*eps1**2)/7.0*(Pmax**7-Pmin**7)+(del2*eps0**2+&
	  2.0*del1*eps0*eps1+del0*eps1**2)/5.0*(Pmax**5-Pmin**5)+&
          (del1*eps0**2+2.0*del0*eps0*eps1)/3.0*(Pmax**3-Pmin**3)+&
          del0*eps0**2*(Pmax-Pmin))
     Icm =  pi**2*4.0*(del2*eps1**2/5.0*(Pmax**5-Pmin**5) + (2.0*del2*&
          eps0*eps1+del1*eps1**2)/3.0*(Pmax**3-Pmin**3)+(del2*eps0**2+& 
	  2.0*del1*eps0*eps1+del0*eps1**2)*(Pmax-Pmin)-(del1*eps0**2+&
          2.0*del0*eps0*eps1)*(1.0/Pmax-1.0/Pmin)-1.0/3.0*del0*eps0**&
          2*(1.0/Pmax**3-1.0/Pmin**3))
     Ic2 = 4.0*( 0.5*Icp + 2.*txi1*Ic - 2.0*tbet0*Icm )


     txi3 = 0.5*tm3sq-E1*E3+(P1**2+P3**2)/2.0
     Idp = pi**2*4.0*(alp2*bet1**2/9.0*(Pmax3**9-Pmin3**9)+(2.0*alp2*&
          bet0*bet1+alp1*bet1**2)/7.0*(Pmax3**7-Pmin3**7)+(alp2*bet0**&
          2+2.0*alp1*bet0*bet1+alp0*bet1**2)/5.0*(Pmax3**5-Pmin3**5)+&
          (alp1*bet0**2+2.0*alp0*bet0*bet1)/3.0*(Pmax3**3-Pmin3**3)+&
          alp0*bet0**2*(Pmax3-Pmin3)) 
     Id2 = 4.0*( -Idp + 2.*txi3*Id )

     Iep = pi**2/15.0*(5.0/7.0*3.0*(Pmax**7-Pmin**7)-3.0/5.0*20.0*A*&
          (Pmax**5-Pmin**5)+60.0/3.0*A**2*(Pmax**3-Pmin**3)) 
     Iem = pi**2/15.0*(5.0/3.0*3.0*(Pmax**3-Pmin**3)-3.0*20.0*A*(Pmax-&
          Pmin)-60.0*A**2*(1.0/Pmax-1.0/Pmin)) 
     Ie2 = 4.0*( 0.5*Iep + 2.*txi1*Ie - 2.0*tbet0*Iem )

     Igp = pi**2/15.0*(5.0/7.0*3.0*(Pmax3**7-Pmin3**7)+3.0/5.0*20.0*E*&
          (Pmax3**5-Pmin3**5)+60.0/3.0*E**2*(Pmax3**3-Pmin3**3)) 
     Ig2 = 4.0*( -Igp + 2.*txi3*Ig )

     IH2 = 4.0*( tm3sq*Ih - 2.0*Id )
     Ij2 = 4.0*( tm3sq*Ij - 2.0*Ih )
     Ik2 = 4.0*( tm3sq*Ik - 2.0*IG )
     Il2 = 4.0*( tm3sq*Il - 2.0*Ik )

     cplA2 = gvsq/VmassSq+gasq/AmassSq+2.0*anti*gva/AVmassSq+anti*2.0*&
          gaf/AFmassSq*m2f/m2*(1.0-dmf/2.0/m2f) 
     cplB2 = gvsq/VmassSq+gasq/AmassSq-2.0*anti*gva/AVmassSq-anti*2.0*&
          gaf/AFmassSq*m2f/m2*(1.0-dmf/2.0/m2f) 
     cplC2 = F2sq/m2**2/FmassSq
     cplD2 =-F2sq/m2**2/FmassSq
     cplE2 = F2sq/FmassSq/m2**2*(-0.5*m3**2+dU2*(xE3-E1)-0.5*dU2**2)
     cplG2 = gvf/VFmassSq*m2f/m2*(2.0-dmf/m2f)+0.5*F2sq/FmassSq/m2**2*&
          (m2f*m4f-Qmass+0.25*m3**2-dU2*(E1+E2f)-dU2**2/4.0) 
     cplH2 = 0.5*F2sq/FmassSq/m2**2*(2.0*Qmass+m3**2+dU2*(3.0*E1-xE3+&
          2.0*E4f)) 
     cplJ2 = gvf/VFmassSq*dmf/m2*0.5*(m3**2-dU2*(E1+xE3))+F2sq/&
          FmassSq/m2**2*0.5*JFF 
     cplK2 = (gasq/AmassSq-gvsq/VmassSq)*m2f*m4f+gvf/VFmassSq*m2f/m2*&
          0.5*(-3.0*m3**2+4.0*dU2*(xE3-E1)-dU2**2+dmf/m2f*(2.0*Qmass+&
          m3**2+dU2*(E4f+2.0*E1-xE3)))+0.5*F2sq/m2**2*KFF/FmassSq 
     cplL2 = gvf/VFmassSq*m2f/m2*dU2*E1*(m3**2-dU2*xE3-0.5*dmf/m2f*&
          (Qmass+0.5*m3**2+dU2*E4f-0.5*dU2**2))+0.5*F2sq/m2**2*LFF/&
          FmassSq 

     Ampsq_ff0 = cplA2*Ia2 + cplB2*Ib2 + cplC2*Ic2 + cplD2*Id2 + &
          cplE2*Ie2 + cplG2*Ig2 + cplH2*Ih2 + cplJ2*Ij2 + cplK2*Ik2 + &
          cplL2*Il2 
     if(Ampsq+ampsq_ff0.ge.0.0) Ampsq=Ampsq+ampsq_ff0           
  end if

!!! adding pseudoscalar terms w/o form factor dependence
  if(opt_pseudo .eq. 1) then

     YAP = 0.25*AmassSq+tm3sq-2.0*E1*E3+P1**2+P3**2
     ZAP = Mpi**2-tm3sq+2.0*E1*E3-P1**2-P3**2
     MassMy = 0.25*AmassSq+Mpi**2

     if(opt_form .eq. 0) then
	if( ZAP >0.0 ) then
	   Inte1 = ( atan(Pmax3/sqrt(ZAP))-atan(Pmin3/sqrt(ZAP)) )/&
         sqrt(ZAP) 
	   Inte2 = 0.5*( Pmax3/(Pmax3**2+ZAP)-Pmin3/(Pmin3**2+ZAp) + &
         Inte1 )/ZAP 
	   ILAp = 4.0*pi**2*Inte1
	   Ilpp = 4.0*pi**2*inte2
	   tIkAp = 2.0*pi**2*((Pmax3-Pmin3)-ZAP*inte1)
	   tIkPP = 2.0*pi**2*(-ZAP*Inte2+Inte1)
	   tIjAp = 4.0*pi**2*alp0*( -1.0/ZAP*(1.0/Pmax3-1.0/Pmin3)  &
		-1.0/ZAP*inte1 )
	   tIgPP = pi**2*( ZAP**2*Inte2 - 2.0*ZAP*Inte1 &
		+(Pmax3-Pmin3) )
	else if(ZAP.lt.0.0) then
	   ZAP = -ZAP
	   TanArg = sqrt(ZAP)*(Pmax3-Pmin3)/(Pmax3*Pmin3-ZAP)
	   LogArg = abs((1.0+TanArg)/(1.0-TanArg))
	   Inte1 = 0.5*log(LogArg)/sqrt(ZAP)
	   Inte2 = -0.5*( Pmax3/(Pmax3**2-ZAP)-Pmin3/(Pmin3**2-ZAp) + &
         Inte1 )/ZAp 
	   ILAp = 4.0*pi**2*Inte1
	   Ilpp = 4.0*pi**2*inte2
	   tIkAp = 2.0*pi**2*((Pmax3-Pmin3)+ZAP*inte1)
	   tIkPP = 2.0*pi**2*(ZAP*Inte2+Inte1)
	   tIjAp = 4.0*pi**2*alp0*( 1.0/ZAP*(1.0/Pmax3-1.0/Pmin3)  &
		+1.0/ZAP*inte1 )
	   tIgPP = pi**2*( ZAP**2*Inte2+2.0*ZAP*Inte1  &
		+(Pmax3-Pmin3))
	end if

     else if(opt_form .eq. 1) then
	if( ZAP>0.0) then
	   Inte1 = ( atan(Pmax3/sqrt(ZAP))-atan(Pmin3/sqrt(ZAP)) )/sqrt(ZAP)
	   Inte2 = 0.5*( Pmax3/(Pmax3**2+ZAP)-Pmin3/(Pmin3**2+ZAp) + &
         Inte1 )/ZAP 
	   ILAp = 16.0*pi**2/AmassSq*(MassMy*Inte1-(Pmax3-Pmin3))
	   Ilpp = 16.0*pi**2/AmassSq*( MassMy*inte2 - inte1 )
	   tIkAp = 8.0*pi**2/AmassSq*( MassMy*(Pmax3-Pmin3)-(Pmax3**3-&
         Pmin3**3)/3.0-MassMy*ZAP*inte1 ) 
	   tIkPP = 8.0*pi**2/AmassSq*( -(Pmax3-Pmin3)-ZAP*MassMy*&
         Inte2+(MassMy+ZAP)*Inte1 ) 
	   tIjAp = 16.0*pi**2/AmassSq*alp0*( -(MassMy-ZAP)/ZAP*(1.0/&
         Pmax3-1.0/Pmin3)-MassMy/ZAP*inte1 )
	   tIgPP = 4.0*pi**2/AmassSq*( ZAP**2*MassMy*Inte2 - (2.0*ZAP*&
         MassMy+ZAp**2)*Inte1+(ZAP+MassMy)*(Pmax3-Pmin3)-(Pmax3**3-&
         Pmin3**3)/3.0 ) 
	else if (ZAP.lt.0.0) then
	   ZAP = -ZAP
	   TanArg = sqrt(ZAP)*(Pmax3-Pmin3)/(Pmax3*Pmin3-ZAP)
	   LogArg = abs((1.0+TanArg)/(1.0-TanArg))
	   Inte1 = 0.5*log(LogArg)/sqrt(ZAP)
	   Inte2 = -0.5*( Pmax3/(Pmax3**2-ZAP)-Pmin3/(Pmin3**2-ZAp) + &
         Inte1 )/ZAp 
	   ILAp = 16.0*pi**2/AmassSq*(MassMy*Inte1-(Pmax3-Pmin3))
	   Ilpp = 16.0*pi**2/AmassSq*( MassMy*inte2 - inte1 )
	   tIkAp = 8.0*pi**2/AmassSq*( MassMy*(Pmax3-Pmin3)-(Pmax3**3-&
         Pmin3**3)/3.0+MassMy*ZAP*inte1 ) 
	   tIkPP = 8.0*pi**2/AmassSq*( -(Pmax3-Pmin3)+ZAP*MassMy*&
         Inte2+(MassMy-ZAP)*Inte1 ) 
	   tIjAp = 16.0*pi**2/AmassSq*alp0*( (MassMy+ZAP)/ZAP*(1.0/&
         Pmax3-1.0/Pmin3)+MassMy/ZAP*inte1 )
	   tIgPP = 4.0*pi**2/AmassSq*( ZAP**2*MassMy*Inte2 - (-2.0*&
         ZAP*MassMy+ZAp**2)*Inte1 +(-ZAP+MassMy)*(Pmax3-Pmin3)-(Pmax3**&
         3-Pmin3**3)/3.0 ) 
	end if
     end if

     IkAP = E*IlAP+tIkAP
     IkPP = E*IlPP+tIkPP
     IjAP = alp1*IlAP+0.5*tIkAP+tIjAP
     IgPP = E**2*IlPP+2.0*E*tIkPP+tIgPP

     cplLAP = 2.0*m2*gasq*dU2*E1*(2.0*m2f*m3**2-dmf*(Qmass+0.5*m3**2)-&
          dU2*(2.0*E3*m2f+E4f*dmf)+0.5*dU2**2*dmf) 
     cplKAP = 2.0*m2*gasq*(m2f*(dU2**2-m3**2)+dmf*dU2*(E2f+E1))
     cplJAP = 2.0*m2*gasq*dmf*(m3**2-dU2*(E1+E3))

     cplLPP = 2.0*m2**2*gasq*dU2*E1*(m3**4-m3**2*dmf**2+dU2*(-m3**2*&
          (3.0*E3-2.0*E1)+dmf**2*E3)+dU2**2*(m3**2+2.0*E3*(E3-E1))-&
          dU2**3*E3) 
     cplKPP = 2.0*m2**2*gasq*(-0.5*m3**4+0.5*dmf**2*m3**2+dU2*m3**2*&
          (E3-3.0*E1)+dU2**2*(2.0*E1*E3-0.5*dmf**2)+dU2**3*(E4f-E2f)-&
          0.5*dU2**4)  
     cplGPP = 2.0*m2**2*gasq*(m3**2-dU2**2)

     ampsq_pseudo = cplLAP*IlAP+cplKAP*IkAP+cplJAP*IjAP + cplLPP*IlPP +&
          cplKPP*IkPP+cplGPP*IgPP 
     if(Ampsq+ampsq_pseudo.ge.0.0) Ampsq=Ampsq+ampsq_pseudo
  end if


  return
end subroutine Calc_Ampsq


subroutine Ebounds( xpn3, xq0_min, xq0_max )
  use CC_global
  implicit none
  real*8:: xpn3, q3_lim(4), tmp1,tmp2,xA,xB,xC, xq3
  integer:: i
  real*8:: tmp3,tmp4,xq0_min,xq0_max,xq0, xq0_l,xq0_h

  En = sqrt( xpn3**2 + mn**2 ) + Un

  q3_lim = 0.
  tmp1 = Up-En-Enu
  tmp2 = Enu**2 + me**2 - mp**2 - xpn3**2

  xA = 4.d0*( tmp1**2-(Enu-xpn3)**2 )
  xB = 4.d0*( (tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*xpn3 )
  xC = 4.d0*tmp1**2*(xpn3**2+mp**2) - (tmp2-tmp1**2)**2
  if( xB**2 - 4.0*xA*xC .ge. 0.0) then
     q3_lim(1) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
     q3_lim(2) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  end if

  xA = 4.d0*( tmp1**2-(Enu+xpn3)**2 )
  xB = 4.d0*( -(tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*xpn3 )
  if( xB**2 - 4.0*xA*xC .ge. 0.0) then
     q3_lim(3) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
     q3_lim(4) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  end if

  xq0_min = 1.d30
  xq0_max = -1.d30
  xq0 = 0.0
  do i=1,4
     xq3 = q3_lim(i)
     tmp1 = Enu - sqrt( (Enu-xq3)**2 + me**2 )
     tmp2 = Enu - sqrt( (Enu+xq3)**2 + me**2 )
     tmp3 = sqrt( (xpn3+xq3)**2+mp**2) + Up - En
     tmp4 = sqrt( (xpn3-xq3)**2+mp**2) + Up - En

     if( abs((tmp1-tmp3)/tmp3)<1.d-6 .or. abs((tmp1-tmp4)/tmp4)<1.d-6 ) then
	xq0 = tmp1
     end if
     if( abs((tmp2-tmp3)/tmp3)<1.d-6 .or. abs((tmp2-tmp4)/tmp4)<1.d-6 ) then
	xq0 = tmp2
     end if
     xq0_min = min( xq0,  xq0_min)
     xq0_max = max( xq0,  xq0_max)
  end do

  xq0_l = Enu - me
  tmp1 = sqrt( (xpn3+Enu)**2+mp**2) + Up - En
  tmp2 = sqrt( (xpn3-Enu)**2+mp**2) + Up - En
  if( xq0_l .le. tmp1 .and. xq0_l .ge. tmp2) then
     xq0_max = max( xq0_l, xq0_max)
  end if

  xq0_h = mp + Up - En
  tmp1 = Enu - sqrt( (Enu-xpn3)**2 + me**2 )
  tmp2 = Enu - sqrt( (Enu+xpn3)**2 + me**2 )
  if( xq0_h .le. tmp1 .and. xq0_h .ge. tmp2) then
     xq0_min = min( xq0_h, xq0_min)
  end if

end subroutine Ebounds


subroutine Range_pn()
  use CC_global
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

12 continue

end subroutine Range_pn


subroutine Ebounds_D( xpn3, xq0_min, xq0_max )
  use CC_global
  implicit none
  real*8:: xpn3, q3_lim(4), tmp1,tmp2,xA,xB,xC, xq3
  integer:: i
  real*8:: tmp3,tmp4,xq0_min,xq0_max,xq0, xq0_l,xq0_h

  En = sqrt( xpn3**2 + mn**2 ) + Un

  q3_lim = 0.
  tmp1 = Up-En-Enu
  tmp2 = Enu**2 + me**2 - mp**2 - xpn3**2

  xA = 4.d0*( tmp1**2-(Enu+xpn3)**2 )
  xB = 4.d0*( (tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*xpn3 )
  xC = 4.d0*tmp1**2*(xpn3**2+mp**2) - (tmp2-tmp1**2)**2
  q3_lim(1) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  q3_lim(2) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )

  xA = 4.d0*( tmp1**2-(Enu-xpn3)**2 )
  xB = 4.d0*( -(tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*xpn3 )
  if( xB**2 - 4.0*xA*xC .ge. 0 ) then
     q3_lim(3) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
     q3_lim(4) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  end if

  xq0_min = 1.d10
  xq0_max = -1.d10
  do i=1,4
     xq3 = q3_lim(i)
     tmp1 = Enu + sqrt( (Enu+xq3)**2 + me**2 )
     tmp2 = Enu + sqrt( (Enu-xq3)**2 + me**2 )
     tmp3 = sqrt( (xpn3+xq3)**2+mp**2) + Up - En
     tmp4 = sqrt( (xpn3-xq3)**2+mp**2) + Up - En
     if( abs((tmp1-tmp3)/tmp3)<1.d-6 .or. abs((tmp1-tmp4)/tmp4)<1.d-6 ) then
	xq0 = tmp1
     end if
     if( abs((tmp2-tmp3)/tmp3)<1.d-6 .or. abs((tmp2-tmp4)/tmp4)<1.d-6 ) then
	xq0 = tmp2
     end if
     xq0_min = min( xq0,  xq0_min)
     xq0_max = max( xq0,  xq0_max)
  end do
  if( abs(xq0_max-xq0_min)<1.d-15 ) then
     if(opt0.eq.3) then
	xq0_max = max(Enu + mue + 50.0*Tem, xq0_min + 50.0*Tem)
     else if(opt0.eq.4) then
	xq0_max = max(xq0_min+50.0*Tem, sqrt(ppmax**2+mp**2)+Up-En)
     end if
  end if

  xq0_l = Enu + me
  tmp1 = sqrt( (xpn3+Enu)**2+mp**2) + Up - En
  tmp2 = sqrt( (xpn3-Enu)**2+mp**2) + Up - En
  if( xq0_l .le. tmp1 .and. xq0_l .ge. tmp2) then
     xq0_min = min( xq0_l, xq0_min)
  end if

  xq0_h = mp + Up - En
  tmp1 = Enu + sqrt( (Enu+xpn3)**2 + me**2 )
  tmp2 = Enu + sqrt( (Enu-xpn3)**2 + me**2 )
  if( xq0_h .le. tmp1 .and. xq0_h .ge. tmp2) then
     xq0_min = min( xq0_h, xq0_min)
  end if

end subroutine Ebounds_D

subroutine Range_pn_D()
  use CC_global
  implicit none
  real*8:: tmp,xpn3,mpt,Epr,Esq,Equ,p2a,p2b,p20
  integer:: i
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

12 continue
end subroutine Range_pn_D


!=======================================================================
!
!     NUMERICAL RECIPES: Gauss-Legendre integration
!
!=======================================================================
subroutine gauleg(x1,x2,x,w,n)
  implicit none
  integer :: n
  real*8 :: x1,x2
  real*8, dimension(n) :: x,w
  !-----------------------------------------------------------------------
  !
  !	input:
  !	x1   lower integration boundary
  !	x2   upper integration boundary
  !	n    discretization
  !
  !	output:
  !	x   position of integrand evaluation
  !	w   weight
  !
  !-----------------------------------------------------------------------
  integer :: i,j,m
  real*8, parameter :: eps=3.d-14,pi=3.141592654d0
  real*8 :: p1,p2,p3,pp,xl,xm,z,z1

  m = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)

  do i=1,m
     z = cos(pi*(i-0.25d0)/(n+0.5d0))

1    continue
     p1 = 1.d0
     p2 = 0.d0
     do j=1,n
	p3 = p2
	p2 = p1
	p1 = ((2.d0*j-1.d0)*z*p2 - (j-1.d0)*p3)/j
     enddo
     pp = n*(z*p1-p2)/(z*z-1.d0)
     z1 = z
     z = z1 - p1/pp
     if (abs(z-z1).gt.eps) goto 1
     x(i) = xm - xl*z
     x(n+1-i) = xm + xl*z
     w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
     w(n+1-i) = w(i)

  enddo
  return
end subroutine gauleg


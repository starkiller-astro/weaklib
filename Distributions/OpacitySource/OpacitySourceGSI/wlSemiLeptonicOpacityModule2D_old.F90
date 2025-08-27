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

MODULE wlSemiLeptonicOpacityModule2D_old

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, mpi, Vud, &
   massA, massV, gamma_p, gamma_n, hbarc, hbar

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PARAMETER :: F2wm0 = gamma_p - gamma_n - 1.0d0
  ! ! If you want to reproduce exactly the GSI numbers you need:
  ! REAL(DP),  PARAMETER :: Gw_MeV=1.166d-11,Vud=0.97427d0, F2wm0 = 3.706d0, &
  !      ga =1.2723d0,gv=1.d0, Mpi=139.57d0, Mnp=938.919d0, &
  !      massA = 1.0d3, massV = 840.d0, &
  !      Dnp=1.293d0, massA=1.0d3,massV=840.d0 , GfVud2 = (Gw_MeV*Vud)**2

  REAL(DP), PARAMETER :: Tfac = 150.0d0
  REAL(DP), PARAMETER :: gasq = ga**2, gvsq = gv**2, gva = gv*ga
  INTEGER,  PARAMETER :: Ngrids = 50      !!! to be tuned by user     

  REAL(DP) :: lepton_mass, neutron_mass, proton_mass, Un, Up
  REAL(DP) :: T, Mun, Mup, Mul, ml2, mn2, mp2
  REAL(DP) :: Enu, En, pnmax, pnmin, ppmax
  INTEGER  :: anti, opt0, opt_form=0, opt_pseudo=0, opt_WM=0

! #if defined(USE_OMP)
!   !$omp threadprivate(lepton_mass, neutron_mass, proton_mass, Un, Up, &
!   !$omp               T, Mun, Mup, Mul, ml2, mn2, mp2,                &
!   !$omp               Enu, En, pnmax, pnmin, ppmax,                   &
!   !$omp               anti, opt0, opt_form, opt_pseudo, opt_WM)
! #endif

!========================================================================
! Routine to calculate the neutrino CC opacity via 2D integrals
! Inputs:
! (1) opt_had=0,1,2,3 for LO, LO+weak magnetism(WM),
! LO+WM+pseduoscalar(PS), LO+WM+PS+form factor dependencies
! (2) opt=opt0=1,2,3,4 for nu capture on n, nub capture on p, inverse
! neutron decay, neutron decay emissivity
! (3) xEnu(NP): neutrino energy array in MeV with dimension NP
! (4) xT: Temperature in MeV
! (5) xMul, xml: relativistic chemical potentials & masses of
! e^-(mu^-) in MeV
! (6) xMun,xMup,xmn,xmp,xUn,xUp: relativistic chemical
! potentials, masses, & potentials of neutron and proton in MeV
! relativsitic disperstion relation for nucleons,
! En=SQRT(p^2+neutron_mass^2)+Un
! Outputs:
! (7) OpaA(NP): opacities in km^-1; for neutron decay, emissivities
! in (s cm^3 MeV)^-1
!     **** blocking factor of final state neutrinos in neutron decay
!     is not considered ****
!=============================================================
  PUBLIC :: Opacity_CC_2D

CONTAINS

SUBROUTINE Opacity_CC_2D(opt_had, opt, xEnu, OpaA, xT, xMul, xMun, xMup,&
     xml, xmn, xmp, xUn, xUp) 

  INTEGER :: opt,opt_had
  REAL(DP), INTENT(IN) :: xEnu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp 
  REAL(8), INTENT(OUT) :: OpaA

  opt0 = opt

  SELECT CASE(opt_had)
  CASE (0)
     opt_wm=0
     opt_pseudo=0
     opt_form=0
  CASE (1)
     opt_wm=1
     opt_pseudo=0
     opt_form=0
  CASE (2)
     opt_wm=1
     opt_pseudo=1
     opt_form=0
  CASE (3)
     opt_wm=1
     opt_pseudo=1
     opt_form=1
  END SELECT

  lepton_mass = xml

  T = xT
  Mul = xMul
  IF(opt0.eq.2) THEN
     Mul = -xMul
  END IF
  SELECT CASE(opt0)
  CASE (1)
     neutron_mass = xmn
     proton_mass = xmp
     Un = xUn
     Up = xUp
     Mun = xMun
     Mup = xMup
     anti = 1
  CASE (2,3,4)
     proton_mass = xmn
     neutron_mass = xmp
     Up = xUn
     Un = xUp
     Mup = xMun
     Mun = xMup
     anti = -1
  END SELECT

  mn2 = neutron_mass**2
  mp2 = proton_mass**2
  ml2 = lepton_mass**2

  SELECT CASE(opt0)
  CASE (1,2)
    CALL Integral_2D(xEnu, OpaA)
  CASE (3,4)
    CALL Integral_2D_D(xEnu, OpaA)
  END SELECT

END SUBROUTINE Opacity_CC_2D

SUBROUTINE Integral_2D(xEnu,res)   ! captures

  REAL(DP) :: xpn3,xq0_min,xq0_max,Ee_min,Ee_max,xEe,En_min,En_max
  REAL(DP) :: xa(Ngrids),wxa(Ngrids), wEn(Ngrids),wEe(Ngrids)
  INTEGER :: i, j
  REAL(DP) :: xEnu,res,xamp,Ep
  REAL(DP) :: xf2,xf3,xf4

  Enu = xEnu
  res = 0.0d0
  CALL Range_pn()
  IF( pnmin > pnmax ) RETURN
  En_min = SQRT(pnmin**2 + mn2) + Un
  En_max = SQRT(pnmax**2 + mn2) + Un

  CALL gauleg(0.d0, 1.d0, xa, wxa, Ngrids)

  DO i=1, Ngrids
    En = xa(i)*(En_max-En_min)+En_min
    wEn(i) = wxa(i)*(En_max-En_min)
    xpn3 = SQRT( (En - Un)**2 - mn2 )
    CALL Ebounds( xpn3, xq0_min, xq0_max)
    Ee_min = Enu - xq0_max
    Ee_max = Enu - xq0_min
    xf2 = 1.0d0/( EXP((En-Mun)/T) + 1.0d0 )
    DO j=1,Ngrids
      xEe = xa(j)*(Ee_max-Ee_min)+Ee_min
      wEe(j) = wxa(j)*(Ee_max-Ee_min)
      Ep = Enu+En-xEe
      CALL Calc_ampsq(En,xEe,xamp)
      xf3 = 1.0d0/( EXP((xEe-Mul)/T) + 1.0d0 )
      xf4 = 1.0d0/( EXP((Ep-Mup)/T) + 1.0d0 )
      res = res + xf2*(1.0d0-xf3)*(1.0d0-xf4)*xamp*wEn(i)*wEe(j)/Enu**2
    END DO
  END DO
  res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc*1.0d10

END SUBROUTINE Integral_2D

SUBROUTINE Integral_2D_D(xEnu,res)   ! decay/inverse decay

  REAL(DP) :: xpn3,xq0_min,xq0_max,Ee_min,Ee_max,xEe,En_min,En_max
  REAL(DP) :: xa(Ngrids),wxa(Ngrids), wEn(Ngrids), wEe(Ngrids)
  INTEGER :: i,j
  REAL(DP) :: xEnu,res,xamp
  REAL(DP) :: E1,mu2,mu4,mu3,U2,U4,Ep
  REAL(DP) :: xf2,xf3,xf4

  U2 = Un
  U4 = Up
  mu2 = Mun
  mu4 = Mup
  mu3 = Mul

  Enu = xEnu
  E1 = xEnu
  res = 0.0d0
  CALL Range_pn_D()
  IF( pnmin > pnmax ) RETURN
  En_min = SQRT(pnmin**2 + mn2) + Un
  En_max = SQRT(pnmax**2 + mn2) + Un

  CALL gauleg(0.d0, 1.d0, xa, wxa, Ngrids)

  DO i=1,Ngrids
    En = xa(i)*(En_max-En_min) + En_min
    wEn(i) = wxa(i)*(En_max-En_min)
    xpn3 = SQRT( (En-Un)**2 - mn2 )
    CALL Ebounds_D( xpn3, xq0_min, xq0_max )
    Ee_min = xq0_min - Enu
    Ee_max = xq0_max - Enu
    xf2 = 1.0d0/( EXP((En-Mun)/T) + 1.0d0 )
    DO j=1,Ngrids
      xEe = xa(j)*(Ee_max-Ee_min)+Ee_min
      wEe(j) = wxa(j)*(Ee_max-Ee_min)
      Ep = Enu+xEe+En
      CALL Calc_ampsq(En,-xEe,xamp)
      xf3 = 1.0d0/( EXP((xEe-Mul)/T) + 1.0d0 )
      xf4 = 1.0d0/( EXP((Ep-Mup)/T) + 1.0d0 )
      IF(opt0.eq.3) THEN ! inverse decay
        res = res + xf2*xf3*(1.0d0-xf4)*xamp*wEn(i)*wEe(j)/Enu**2
      ELSE IF(opt0.eq.4) THEN ! n decay, *blocking of neutrinos is not added*
        res = res + (1.0d0-xf2)*(1.0d0-xf3)*xf4*xamp*wEn(i)*wEe(j)/Enu**2
      END IF
    END DO
  END DO

  IF(opt0.eq.3) THEN
    res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc*1.0d10
  ELSE IF(opt0.eq.4) THEN
    res = res*(Gw_MeV*Vud)**2/32.0d0/(pi**7)*Enu**2/hbarc**3/hbar*1.0d10
  END IF

  res = -res

END SUBROUTINE Integral_2D_D


!!! calculate the amplitudes, keep the following unchanged
PURE ELEMENTAL SUBROUTINE Calc_Ampsq( E2, E3, Ampsq )

  REAL(DP), INTENT(IN)  :: E2, E3
  REAL(DP), INTENT(OUT) :: Ampsq

  REAL(DP) :: P1,P2,P3,P4,Pmax,Pmin,Pmax2,Pmin2,Pmax3,Pmin3
  REAL(DP) :: Ia,Ib,Ic,Id,Ie,Ig,Ih,Ij,Ik,Il,A,B,Ct,D,E,del0,&
       del1,del2,eps0,eps1,alp0,alp1,alp2,bet0,bet1,cplA,cplB,cplC,&
       cplD,cplE,cplG,cplH,cplJ,cplK,cplL,JFF,KFF,LFF 
  REAL(DP) :: m2f,E2f,m4f,E4f,E1,Qmass,dmf,U2,U4,dU2,m3,E4,m2,xE3
  REAL(DP) :: Iam,Iap,Ibm,Ibp,Icm,Icp,Idp,Iem,Iep,Igp
  REAL(DP) :: gaf,gvf,f2sq
  REAL(DP) :: AmassSq,VmassSq,AVmassSq,APmassSq,AFmassSq,VFmassSq,&
       FmassSq,tm3sq,tbet0 
  REAL(DP) :: txi1,Ia2,tgam0,Ib2,txi3,Ic2,Id2,Ie2,Ig2,Ih2,Ij2,Ik2,Il2,&
       cplA2,cplB2,cplC2,cplD2,cplE2,cplG2,cplH2,cplJ2,cplK2,cplL2
  REAL(DP) ::YAP,ZAP,MassMy,Inte1,Inte2,IlAP,IlPP,IkAP,tIkAP,IkPP,tIkPP,&
       IjAP,tIjAP,IGPP,tIgPP,TanArg,LogArg,cplLAP,cplKAP,cplJAP,&
       cplLPP,cplKPP,cplGPP 
  REAL(DP) :: ampsq_ff0,ampsq_pseudo       

  Ampsq = 0.0d0

  m2 = 0.5d0*( proton_mass + neutron_mass )
  m2f = neutron_mass
  m4f = proton_mass
  m3 = lepton_mass
  U2 = Un
  U4 = Up
  dU2 = U2 - U4
  dmf = m2f - m4f
  Qmass = 0.5d0*( m2f**2 - m4f**2 )

  E1 = Enu
  E4 = E1 + E2 - E3
  E4f = E4 - U4
  E2f = E2 - U2
  P1 = E1
  P3 = SQRT( E3**2 - m3**2 )
  P2 = SQRT( (E2 - U2)**2 - m2f**2 )
  P4 = SQRT( (E4 - U4)**2 - m4f**2 )

  IF( P1 .ne. P1) RETURN
  IF( P2 .ne. P2) RETURN
  IF( P3 .ne. P3) RETURN
  IF( P4 .ne. P4) RETURN

  IF ((P1+P2)>(P3+P4)) THEN
     Pmax = P3+P4
  ELSE
     Pmax = P1+P2
  endif
  IF ((ABS(P1-P2))>(ABS(P3-P4))) THEN
     Pmin = ABS(P1-P2)
  ELSE
     Pmin = ABS(P3-P4)
  endif

  A = E1*E2f+0.5d0*(P1**2+P2**2)
  B = E3*E4f+0.5d0*(P3**2+P4**2)
  Ia = pi**2/15.0d0*(3.0d0*(Pmax**5-Pmin**5)-10.0d0*(A+B)*(Pmax**3-Pmin**3)+&
       60.0d0*A*B*(Pmax-Pmin)) 

  IF ((P1+P4)>(P3+P2)) THEN
    Pmax2 = P3+P2
  ELSE
    Pmax2 = P1+P4
  endif
  IF ((ABS(P1-P4))>(ABS(P3-P2))) THEN
    Pmin2 = ABS(P1-P4)
  ELSE
    Pmin2 = ABS(P3-P2)
  endif

  Ct = E1*E4f-0.5d0*(P1**2+P4**2)
  D = E2f*E3-0.5d0*(P2**2+P3**2)
  Ib = pi**2/15.0d0*(3.0d0*(Pmax2**5-Pmin2**5)+10.0d0*(Ct+D)*(Pmax2**3-&
       Pmin2**3)+60.0d0*Ct*D*(Pmax2-Pmin2)) 

  del0 = -(P1**2-P2**2)*(P3**2-P4**2)/4.0d0
  del1 = E1*E3 + (-P1**2+P2**2-P3**2+P4**2)/4.0d0
  del2 = -0.25d0
  eps0 = E1*E2f + (P1**2+P2**2)/2.0d0
  eps1 = -0.5d0
  Ic = pi**2*4.0d0*(del2*eps1**2/7.0d0*(Pmax**7-Pmin**7)+(2.0d0*del2*eps0*&
       eps1+del1*eps1**2)/5.0d0*(Pmax**5-Pmin**5)+(del2*eps0**2+2.0d0*&
       del1*eps0*eps1+del0*eps1**2)/3.0d0*(Pmax**3-Pmin**3)+(del1*eps0**&
       2+2.0d0*del0*eps0*eps1)*(Pmax-Pmin)-del0*eps0**2*(1.0d0/Pmax-1.0d0/&
       Pmin)) 

  IF ((P1+P3)>(P2+P4)) THEN
    Pmax3 = P2+P4
  ELSE
    Pmax3 = P1+P3
  endif
  IF ((ABS(P1-P3))>(ABS(P2-P4))) THEN
    Pmin3 = ABS(P1-P3)
  ELSE
    Pmin3 = ABS(P2-P4)
  endif

  alp0 = (P1**2-P3**2)*(P2**2-P4**2)/4.0d0
  alp1 = E1*E2f + (P1**2+P2**2-P3**2-P4**2)/4.0d0
  alp2 = 0.25d0
  bet0 = E1*E3 - (P1**2+P3**2)/2.0d0
  bet1 = 0.5d0
  Id = pi**2*4.0d0*(alp2*bet1**2/7.0d0*(Pmax3**7-Pmin3**7)+(2.0d0*alp2*bet0*&
       bet1+alp1*bet1**2)/5.0d0*(Pmax3**5-Pmin3**5)+& 
       (alp2*bet0**2+2.0d0*alp1*bet0*bet1+alp0*bet1**2)/3.0d0*(Pmax3**3-&
       Pmin3**3)+(alp1*bet0**2+2.0d0*alp0*bet0*bet1)*& 
       (Pmax3-Pmin3)-alp0*bet0**2*(1.0d0/Pmax3-1.0d0/Pmin3))

  Ie = pi**2/15.0d0*(3.0d0*(Pmax**5-Pmin**5)-20.0d0*A*(Pmax**3-Pmin**3)+&
       60.0d0*A**2*(Pmax-Pmin)) 

  E = E1*E3 - 0.5d0*(P1**2+P3**2)
  Ig = pi**2/15.0d0*(3.0d0*(Pmax3**5-Pmin3**5)+20.0d0*E*(Pmax3**3-Pmin3**3)+&
       60.0d0*E**2*(Pmax3-Pmin3)) 

  Ih = pi**2*4.0d0*(alp2*bet1/5.0d0*(Pmax3**5-Pmin3**5)+(alp2*bet0+alp1*&
       bet1)/3.0d0*(Pmax3**3-Pmin3**3)+(alp1*bet0+alp0*bet1)& 
       *(Pmax3-Pmin3)-alp0*bet0*(1.0d0/Pmax3-1.0d0/Pmin3))

  Ij = pi**2/15.0d0*(-10.0d0*(Pmax**3-Pmin**3)+60.0d0*A*(Pmax-Pmin))

  Ik = pi**2/15.0d0*(10.0d0*(Pmax3**3-Pmin3**3)+60.0d0*E*(Pmax3-Pmin3))

  Il = pi**2/15.0d0*(60.0d0*(Pmax-Pmin))

  IF(opt_WM .eq. 1) THEN
     gvf = gv*F2wm0
     gaf = ga*F2wm0
     f2sq = F2wm0**2
  ELSE
     gvf = 0.0d0
     gaf = 0.0d0
     f2sq = 0.0d0
  END IF

  xE3 = E3
  cplA = (gvsq+2.0d0*anti*gva+gasq)+anti*2.0d0*gaf*m2f/m2*(1.0d0-dmf/2.0d0/m2f)
  cplB = (gvsq-2.0d0*anti*gva+gasq)-anti*2.0d0*gaf*m2f/m2*(1.0d0-dmf/2.0d0/m2f)
  cplC = F2sq/m2**2
  cplD =-F2sq/m2**2
  cplE = F2sq/m2**2*(-0.5d0*m3**2+dU2*(xE3-E1)-0.5d0*dU2**2)
  cplG = gvf*m2f/m2*(2.0d0-dmf/m2f)+0.5d0*F2sq/m2**2*(m2f*m4f-Qmass+0.25d0*&
       m3**2-dU2*(E1+E2f)-dU2**2/4.0d0) 
  cplH = 0.5d0*F2sq/m2**2*(2.0d0*Qmass+m3**2+dU2*(3.0d0*E1-xE3+2.0d0*E4f))

  JFF = -m3**2*(E1+0.5d0*E2f+0.5d0*E4f)+Qmass*(xE3-3.0d0*E1) &
       +0.5d0*dU2*(E4f*(3.0d0*xE3-5.0d0*E1)+E2f*(xE3+E1)+xE3**2-E1**2-2.0d0*Qmass) &
       + dU2**2*(E1-xE3-E4f)+0.5d0*dU2**3
  JFF = JFF*dU2
  cplJ = gvf*dmf/m2*0.5d0*(m3**2-dU2*(E1+xE3))+F2sq/m2**2*0.5d0*JFF

  kFF = -(m2f+3.0d0*m4f)*m2f*0.25d0*m3**2 + Qmass**2 + Qmass*0.25d0*m3**2-m3**4/8.0d0 &
       +dU2*(0.5d0*Qmass*(3.0d0*E1-E2f+xE3+3.0d0*E4f)+0.25d0*m3**2*(2.0d0*E2f+xE3+&
       E1)+m2f*m4f*(xE3-E1)) & 
       +dU2**2*( 0.25d0*(m2f**2-m2f*m4f-3.0d0*Qmass)+E4f*0.5d0*(2.0d0*E1-E2f+&
       xE3+E4f)+E2f*(0.5d0*E1-xE3)+0.5d0*E1**2) & 
       +dU2**3*0.25d0*( -E1+2.0d0*E2f-xE3-2.0d0*E4f ) + 0.125d0*dU2**4
  cplK = (gasq-gvsq)*m2f*m4f+gvf*m2f/m2*0.5d0*(-3.0d0*m3**2+4.0d0*dU2*(xE3-&
       E1)-dU2**2 & 
       +dmf/m2f*(2.0d0*Qmass+m3**2+dU2*(E4f+2.0d0*E1-xE3)))&
       +0.5d0*F2sq/m2**2*KFF

  LFF = m3**2*(m2f+m4f)**2-4.0d0*Qmass**2+dU2*(-m3**2*(E2f+E4f+E1)-2.0d0*&
       xE3*(m2f**2+m2f*m4f) & 
       +2.0d0*Qmass*(E2f-3.0d0*E4f-E1)) &
       +2.0d0*dU2**2*( E2f*xE3+E4f*(E2f-E4f-E1) + Qmass ) + dU2**3*(-E2f+E4f+E1)
  LFF = LFF*dU2*E1*0.25d0
  cplL = gvf*m2f/m2*dU2*E1*(m3**2-dU2*xE3-0.5d0*dmf/m2f*(Qmass+0.5d0*m3**&
       2+dU2*E4f-0.5d0*dU2**2))   & 
       +0.5d0*F2sq/m2**2*LFF

  Ampsq = max(0.0d0,cplA*Ia+cplB*Ib+cplC*Ic+cplD*Id+cplE*Ie+cplG*Ig+cplH*Ih+&
       cplJ*Ij+cplK*Ik+cplL*Il) 

!!! adding form-factor dependences
  AmassSq = massA**2
  APmassSq = 2.0d0*AmassSq
  tm3sq = m3**2 + 2.0d0*dU2*(E1-E3) + dU2**2
  IF(opt_form .eq. 1) THEN
    VmassSq = 1.0d0/(1.0d0/massV**2-F2wm0/m2**2/8.0d0)
    AVmassSq = 1.0d0/(0.5d0/massA**2+0.5d0/massV**2-F2wm0/m2**2/16.0d0)
    AFmassSq = 1.0d0/(0.5d0/massA**2+0.5d0/massV**2+1.0d0/m2**2/16.0d0)
    VFmassSq = 1.0d0/(1.0d0/massV**2-(F2wm0-1.0d0)/16.0d0/m2**2)
    FmassSq = 1.0d0/(1.0d0/massV**2+1.0d0/8.0d0/m2**2)

    tbet0 = (P1**2-P2**2)*(P4**2-P3**2)/4.0d0
    txi1 = 0.5d0*tm3sq-E1*E3+(P1**2-P2**2+P3**2-P4**2)/4.0d0

    Iap = pi**2/15.0d0*(5.0d0/7.0d0*3.0d0*(Pmax**7-Pmin**7)-3.0d0/5.0d0*10.0d0*(A+B)*&
        (Pmax**5-Pmin**5)+60.0d0/3.0d0*A*B*(Pmax**3-Pmin**3)) 
    Iam = pi**2/15.0d0*(5.0d0/3.0d0*3.0d0*(Pmax**3-Pmin**3)-3.0d0*10.0d0*(A+B)*&
        (Pmax-Pmin)-60.0d0*A*B*(1.0d0/Pmax-1.0d0/Pmin)) 
    Ia2 = 4.0d0*( 0.5d0*Iap + 2.0d0*txi1*Ia - 2.0d0*tbet0*Iam )

    tgam0 = (P1**2-P4**2)*(P2**2-P3**2)/4.0d0

    Ibp = pi**2/15.0d0*(5.0d0/7.0d0*3.0d0*(Pmax2**7-Pmin2**7)+3.0d0/5.0d0*10.0d0*&
        (Ct+D)*(Pmax2**5-Pmin2**5)+60.0d0/3.0d0*Ct*D*(Pmax2**3-Pmin2**&
        3)) 
    Ibm = pi**2/15.0d0*(5.0d0/3.0d0*3.0d0*(Pmax2**3-Pmin2**3)+3.0d0*10.0d0*(Ct+&
        D)*(Pmax2-Pmin2)-60.0d0*Ct*D*(1.0d0/Pmax2-1.0d0/Pmin2)) 
    Ib2 = 4.0d0*( 0.5d0*Ibp + 2.0d0*txi1*Ib - 2.0d0*tgam0*Ibm )

    Icp = pi**2*4.0d0*(del2*eps1**2/9.0d0*(Pmax**9-Pmin**9)+(2.0d0*del2*&
        eps0*eps1+del1*eps1**2)/7.0d0*(Pmax**7-Pmin**7)+(del2*eps0**2+&
        2.0d0*del1*eps0*eps1+del0*eps1**2)/5.0d0*(Pmax**5-Pmin**5)+&
        (del1*eps0**2+2.0d0*del0*eps0*eps1)/3.0d0*(Pmax**3-Pmin**3)+&
        del0*eps0**2*(Pmax-Pmin))
    Icm =  pi**2*4.0d0*(del2*eps1**2/5.0d0*(Pmax**5-Pmin**5) + (2.0d0*del2*&
        eps0*eps1+del1*eps1**2)/3.0d0*(Pmax**3-Pmin**3)+(del2*eps0**2+& 
        2.0d0*del1*eps0*eps1+del0*eps1**2)*(Pmax-Pmin)-(del1*eps0**2+&
        2.0d0*del0*eps0*eps1)*(1.0d0/Pmax-1.0d0/Pmin)-1.0d0/3.0d0*del0*eps0**&
        2*(1.0d0/Pmax**3-1.0d0/Pmin**3))

    Ic2 = 4.0d0*( 0.5d0*Icp + 2.0d0*txi1*Ic - 2.0d0*tbet0*Icm )


    txi3 = 0.5d0*tm3sq-E1*E3+(P1**2+P3**2)/2.0d0
    Idp = pi**2*4.0d0*(alp2*bet1**2/9.0d0*(Pmax3**9-Pmin3**9)+(2.0d0*alp2*&
        bet0*bet1+alp1*bet1**2)/7.0d0*(Pmax3**7-Pmin3**7)+(alp2*bet0**&
        2+2.0d0*alp1*bet0*bet1+alp0*bet1**2)/5.0d0*(Pmax3**5-Pmin3**5)+&
        (alp1*bet0**2+2.0d0*alp0*bet0*bet1)/3.0d0*(Pmax3**3-Pmin3**3)+&
        alp0*bet0**2*(Pmax3-Pmin3)) 
    Id2 = 4.0d0*( -Idp + 2.0d0*txi3*Id )

    Iep = pi**2/15.0d0*(5.0d0/7.0d0*3.0d0*(Pmax**7-Pmin**7)-3.0d0/5.0d0*20.0d0*A*&
        (Pmax**5-Pmin**5)+60.0d0/3.0d0*A**2*(Pmax**3-Pmin**3)) 
    Iem = pi**2/15.0d0*(5.0d0/3.0d0*3.0d0*(Pmax**3-Pmin**3)-3.0d0*20.0d0*A*(Pmax-&
        Pmin)-60.0d0*A**2*(1.0d0/Pmax-1.0d0/Pmin)) 
    Ie2 = 4.0d0*( 0.5d0*Iep + 2.0d0*txi1*Ie - 2.0d0*tbet0*Iem )

    Igp = pi**2/15.0d0*(5.0d0/7.0d0*3.0d0*(Pmax3**7-Pmin3**7)+3.0d0/5.0d0*20.0d0*E*&
        (Pmax3**5-Pmin3**5)+60.0d0/3.0d0*E**2*(Pmax3**3-Pmin3**3)) 
    Ig2 = 4.0d0*( -Igp + 2.0d0*txi3*Ig )

    IH2 = 4.0d0*( tm3sq*Ih - 2.0d0*Id )
    Ij2 = 4.0d0*( tm3sq*Ij - 2.0d0*Ih )
    Ik2 = 4.0d0*( tm3sq*Ik - 2.0d0*IG )
    Il2 = 4.0d0*( tm3sq*Il - 2.0d0*Ik )

    cplA2 = gvsq/VmassSq+gasq/AmassSq+2.0d0*anti*gva/AVmassSq+anti*2.0d0*&
        gaf/AFmassSq*m2f/m2*(1.0d0-dmf/2.0d0/m2f) 
    cplB2 = gvsq/VmassSq+gasq/AmassSq-2.0d0*anti*gva/AVmassSq-anti*2.0d0*&
        gaf/AFmassSq*m2f/m2*(1.0d0-dmf/2.0d0/m2f) 
    cplC2 = F2sq/m2**2/FmassSq
    cplD2 =-F2sq/m2**2/FmassSq
    cplE2 = F2sq/FmassSq/m2**2*(-0.5d0*m3**2+dU2*(xE3-E1)-0.5d0*dU2**2)
    cplG2 = gvf/VFmassSq*m2f/m2*(2.0d0-dmf/m2f)+0.5d0*F2sq/FmassSq/m2**2*&
        (m2f*m4f-Qmass+0.25d0*m3**2-dU2*(E1+E2f)-dU2**2/4.0d0) 
    cplH2 = 0.5d0*F2sq/FmassSq/m2**2*(2.0d0*Qmass+m3**2+dU2*(3.0d0*E1-xE3+&
        2.0d0*E4f)) 
    cplJ2 = gvf/VFmassSq*dmf/m2*0.5d0*(m3**2-dU2*(E1+xE3))+F2sq/&
        FmassSq/m2**2*0.5d0*JFF 
    cplK2 = (gasq/AmassSq-gvsq/VmassSq)*m2f*m4f+gvf/VFmassSq*m2f/m2*&
        0.5d0*(-3.0d0*m3**2+4.0d0*dU2*(xE3-E1)-dU2**2+dmf/m2f*(2.0d0*Qmass+&
        m3**2+dU2*(E4f+2.0d0*E1-xE3)))+0.5d0*F2sq/m2**2*KFF/FmassSq 
    cplL2 = gvf/VFmassSq*m2f/m2*dU2*E1*(m3**2-dU2*xE3-0.5d0*dmf/m2f*&
        (Qmass+0.5d0*m3**2+dU2*E4f-0.5d0*dU2**2))+0.5d0*F2sq/m2**2*LFF/&
        FmassSq 

    Ampsq_ff0 = cplA2*Ia2 + cplB2*Ib2 + cplC2*Ic2 + cplD2*Id2 + &
        cplE2*Ie2 + cplG2*Ig2 + cplH2*Ih2 + cplJ2*Ij2 + cplK2*Ik2 + &
        cplL2*Il2 
     IF(Ampsq+ampsq_ff0>=0.0d0) Ampsq=Ampsq+ampsq_ff0           
  END IF

!!! adding pseudoscalar terms w/o form factor dependence
  IF(opt_pseudo .eq. 1) THEN

    YAP = 0.25d0*AmassSq+tm3sq-2.0d0*E1*E3+P1**2+P3**2
    ZAP = mpi**2-tm3sq+2.0d0*E1*E3-P1**2-P3**2
    MassMy = 0.25d0*AmassSq+mpi**2

    IF(opt_form .eq. 0) THEN
      IF( ZAP >0.0d0 ) THEN
        Inte1 = ( atan(Pmax3/sqrt(ZAP))-atan(Pmin3/sqrt(ZAP)) )/&
            sqrt(ZAP) 
        Inte2 = 0.5d0*( Pmax3/(Pmax3**2+ZAP)-Pmin3/(Pmin3**2+ZAp) + &
            Inte1 )/ZAP 
        ILAp = 4.0d0*pi**2*Inte1
        Ilpp = 4.0d0*pi**2*inte2
        tIkAp = 2.0d0*pi**2*((Pmax3-Pmin3)-ZAP*inte1)
        tIkPP = 2.0d0*pi**2*(-ZAP*Inte2+Inte1)
        tIjAp = 4.0d0*pi**2*alp0*( -1.0d0/ZAP*(1.0d0/Pmax3-1.0d0/Pmin3)  &
        -1.0d0/ZAP*inte1 )
        tIgPP = pi**2*( ZAP**2*Inte2 - 2.0d0*ZAP*Inte1 &
        +(Pmax3-Pmin3) )
      ELSE IF(ZAP<0.0d0) THEN
        ZAP = -ZAP
        TanArg = sqrt(ZAP)*(Pmax3-Pmin3)/(Pmax3*Pmin3-ZAP)
        LogArg = abs((1.0d0+TanArg)/(1.0d0-TanArg))
        Inte1 = 0.5d0*log(LogArg)/sqrt(ZAP)
        Inte2 = -0.5d0*( Pmax3/(Pmax3**2-ZAP)-Pmin3/(Pmin3**2-ZAp) + &
            Inte1 )/ZAp 
        ILAp = 4.0d0*pi**2*Inte1
        Ilpp = 4.0d0*pi**2*inte2
        tIkAp = 2.0d0*pi**2*((Pmax3-Pmin3)+ZAP*inte1)
        tIkPP = 2.0d0*pi**2*(ZAP*Inte2+Inte1)
        tIjAp = 4.0d0*pi**2*alp0*( 1.0d0/ZAP*(1.0d0/Pmax3-1.0d0/Pmin3)  &
        +1.0d0/ZAP*inte1 )
        tIgPP = pi**2*( ZAP**2*Inte2+2.0d0*ZAP*Inte1  &
        +(Pmax3-Pmin3))
      END IF

    ELSE IF(opt_form .eq. 1) THEN
      IF( ZAP>0.0d0) THEN
        Inte1 = ( atan(Pmax3/sqrt(ZAP))-atan(Pmin3/sqrt(ZAP)) )/sqrt(ZAP)
        Inte2 = 0.5d0*( Pmax3/(Pmax3**2+ZAP)-Pmin3/(Pmin3**2+ZAp) + &
            Inte1 )/ZAP 
        ILAp = 16.0d0*pi**2/AmassSq*(MassMy*Inte1-(Pmax3-Pmin3))
        Ilpp = 16.0d0*pi**2/AmassSq*( MassMy*inte2 - inte1 )
        tIkAp = 8.0d0*pi**2/AmassSq*( MassMy*(Pmax3-Pmin3)-(Pmax3**3-&
            Pmin3**3)/3.0d0-MassMy*ZAP*inte1 ) 
        tIkPP = 8.0d0*pi**2/AmassSq*( -(Pmax3-Pmin3)-ZAP*MassMy*&
            Inte2+(MassMy+ZAP)*Inte1 ) 
        tIjAp = 16.0d0*pi**2/AmassSq*alp0*( -(MassMy-ZAP)/ZAP*(1.0d0/&
            Pmax3-1.0d0/Pmin3)-MassMy/ZAP*inte1 )
        tIgPP = 4.0d0*pi**2/AmassSq*( ZAP**2*MassMy*Inte2 - (2.0d0*ZAP*&
            MassMy+ZAp**2)*Inte1+(ZAP+MassMy)*(Pmax3-Pmin3)-(Pmax3**3-&
            Pmin3**3)/3.0d0 ) 
      ELSE IF (ZAP<0.0d0) THEN
        ZAP = -ZAP
        TanArg = sqrt(ZAP)*(Pmax3-Pmin3)/(Pmax3*Pmin3-ZAP)
        LogArg = abs((1.0d0+TanArg)/(1.0d0-TanArg))
        Inte1 = 0.5d0*log(LogArg)/sqrt(ZAP)
        Inte2 = -0.5d0*( Pmax3/(Pmax3**2-ZAP)-Pmin3/(Pmin3**2-ZAp) + &
            Inte1 )/ZAp 
        ILAp = 16.0d0*pi**2/AmassSq*(MassMy*Inte1-(Pmax3-Pmin3))
        Ilpp = 16.0d0*pi**2/AmassSq*( MassMy*inte2 - inte1 )
        tIkAp = 8.0d0*pi**2/AmassSq*( MassMy*(Pmax3-Pmin3)-(Pmax3**3-&
            Pmin3**3)/3.0d0+MassMy*ZAP*inte1 ) 
        tIkPP = 8.0d0*pi**2/AmassSq*( -(Pmax3-Pmin3)+ZAP*MassMy*&
            Inte2+(MassMy-ZAP)*Inte1 ) 
        tIjAp = 16.0d0*pi**2/AmassSq*alp0*( (MassMy+ZAP)/ZAP*(1.0d0/&
            Pmax3-1.0d0/Pmin3)+MassMy/ZAP*inte1 )
        tIgPP = 4.0d0*pi**2/AmassSq*( ZAP**2*MassMy*Inte2 - (-2.0d0*&
            ZAP*MassMy+ZAp**2)*Inte1 +(-ZAP+MassMy)*(Pmax3-Pmin3)-(Pmax3**&
            3-Pmin3**3)/3.0d0 ) 
      END IF
    END IF

      IkAP = E*IlAP+tIkAP
      IkPP = E*IlPP+tIkPP
      IjAP = alp1*IlAP+0.5d0*tIkAP+tIjAP
      IgPP = E**2*IlPP+2.0d0*E*tIkPP+tIgPP

      cplLAP = 2.0d0*m2*gasq*dU2*E1*(2.0d0*m2f*m3**2-dmf*(Qmass+0.5d0*m3**2)-&
            dU2*(2.0d0*E3*m2f+E4f*dmf)+0.5d0*dU2**2*dmf) 
      cplKAP = 2.0d0*m2*gasq*(m2f*(dU2**2-m3**2)+dmf*dU2*(E2f+E1))
      cplJAP = 2.0d0*m2*gasq*dmf*(m3**2-dU2*(E1+E3))

      cplLPP = 2.0d0*m2**2*gasq*dU2*E1*(m3**4-m3**2*dmf**2+dU2*(-m3**2*&
            (3.0d0*E3-2.0d0*E1)+dmf**2*E3)+dU2**2*(m3**2+2.0d0*E3*(E3-E1))-&
            dU2**3*E3) 
      cplKPP = 2.0d0*m2**2*gasq*(-0.5d0*m3**4+0.5d0*dmf**2*m3**2+dU2*m3**2*&
            (E3-3.0d0*E1)+dU2**2*(2.0d0*E1*E3-0.5d0*dmf**2)+dU2**3*(E4f-E2f)-&
            0.5d0*dU2**4)  
      cplGPP = 2.0d0*m2**2*gasq*(m3**2-dU2**2)

      ampsq_pseudo = cplLAP*IlAP+cplKAP*IkAP+cplJAP*IjAP + cplLPP*IlPP +&
            cplKPP*IkPP+cplGPP*IgPP 
      IF(Ampsq+ampsq_pseudo>=0.0d0) Ampsq=Ampsq+ampsq_pseudo
    END IF

END subroutine Calc_Ampsq

SUBROUTINE Ebounds( xpn3, xq0_min, xq0_max )

  REAL(DP) :: xpn3, q3_lim(4), tmp1,tmp2,xA,xB,xC, xq3
  INTEGER :: i
  REAL(DP) :: tmp3,tmp4,xq0_min,xq0_max,xq0, xq0_l,xq0_h

  En = SQRT( xpn3**2 + mn2 ) + Un

  q3_lim = 0.0d0
  tmp1 = Up-En-Enu
  tmp2 = Enu**2 + ml2 - mp2 - xpn3**2

  xA = 4.d0*( tmp1**2-(Enu-xpn3)**2 )
  xB = 4.d0*( (tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*xpn3 )
  xC = 4.d0*tmp1**2*(xpn3**2+mp2) - (tmp2-tmp1**2)**2
  IF( xB**2 - 4.0d0*xA*xC >= 0.0d0) THEN
    q3_lim(1) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
    q3_lim(2) = ABS( (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  END IF

  xA = 4.d0*( tmp1**2-(Enu+xpn3)**2 )
  xB = 4.d0*( -(tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*xpn3 )
  IF( xB**2 - 4.0d0*xA*xC >= 0.0d0) THEN
    q3_lim(3) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
    q3_lim(4) = ABS( (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  END IF

  xq0_min = 1.d30
  xq0_max = -1.d30
  xq0 = 0.0d0
  DO i=1,4
    xq3 = q3_lim(i)
    tmp1 = Enu - SQRT( (Enu-xq3)**2 + ml2 )
    tmp2 = Enu - SQRT( (Enu+xq3)**2 + ml2 )
    tmp3 = SQRT( (xpn3+xq3)**2+mp2) + Up - En
    tmp4 = SQRT( (xpn3-xq3)**2+mp2) + Up - En

    IF( ABS((tmp1-tmp3)/tmp3)<1.d-6 .or. ABS((tmp1-tmp4)/tmp4)<1.d-6 ) THEN
      xq0 = tmp1
    END IF
    IF( ABS((tmp2-tmp3)/tmp3)<1.d-6 .or. ABS((tmp2-tmp4)/tmp4)<1.d-6 ) THEN
      xq0 = tmp2
    END IF
    xq0_min = min( xq0,  xq0_min)
    xq0_max = max( xq0,  xq0_max)
  END DO

  xq0_l = Enu - lepton_mass
  tmp1 = SQRT( (xpn3+Enu)**2+mp2) + Up - En
  tmp2 = SQRT( (xpn3-Enu)**2+mp2) + Up - En
  IF( xq0_l <= tmp1 .and. xq0_l >= tmp2) THEN
     xq0_max = max( xq0_l, xq0_max)
  END IF

  xq0_h = proton_mass + Up - En
  tmp1 = Enu - SQRT( (Enu-xpn3)**2 + ml2 )
  tmp2 = Enu - SQRT( (Enu+xpn3)**2 + ml2 )
  IF( xq0_h <= tmp1 .and. xq0_h >= tmp2) THEN
     xq0_min = min( xq0_h, xq0_min)
  END IF

END SUBROUTINE Ebounds

SUBROUTINE Range_pn()

  REAL(DP) :: tmp,mpt,Epr,Esq,Equ,p2a,p2b,p20,F0,Finf,Fmax,Fmin,pnmax0

  pnmin = 0.0d0
  pnmax = SQRT( (Tfac*T + Mun - Un)**2 - mn2 )
  pnmax0 = pnmax

  mpt = proton_mass + lepton_mass
  F0 = SQRT(Enu**2 + mpt**2) + Up - neutron_mass - Un
  Finf = -Enu + Up - Un
  IF(mpt < neutron_mass) THEN
    tmp = (1.d0+mpt/neutron_mass)/(1.d0-mpt**2/mn2)*Enu
    Fmin = SQRT( (tmp-Enu)**2+mpt**2 )+Up-SQRT(tmp**2+mn2)-Un
    Fmax = max( F0, Finf )
  ELSE
    Fmin = Finf
    Fmax = F0
  END IF

  IF( Fmax <= Enu) RETURN
  IF( Fmin >= Enu) THEN
    pnmax = -1000.0d0
    RETURN
  END IF

  Epr = Enu - Up + Un
  Esq = Enu**2 + mpt**2 - mn2 - Epr**2
  Equ = Esq**2 + 4.d0*mn2*( Enu**2-Epr**2 )
  IF(Equ>0.0d0 .and. Enu**2 .ne. Epr**2) THEN
    p2a = ( Enu*Esq - ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    p2b = ( Enu*Esq + ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    IF(Epr<0.d0 .and. p2b>0.d0) THEN
      pnmin=p2b
    ELSE IF(Epr>=0.d0 .and. p2a>0.d0) THEN
      pnmin=p2a
    END IF

  ELSE IF(Enu**2 .eq. Epr**2) THEN
    p20 = -(4.0d0*Epr**2*mn2-Esq**2)/(4.0d0*Esq*Enu)
    IF(p20>0.d0) pnmin=p20
  END IF

END SUBROUTINE Range_pn

SUBROUTINE Ebounds_D( xpn3, xq0_min, xq0_max )

  REAL(DP) :: xpn3, q3_lim(4), tmp1,tmp2,xA,xB,xC, xq3
  INTEGER :: i
  REAL(DP) :: tmp3,tmp4,xq0_min,xq0_max,xq0, xq0_l,xq0_h

  En = SQRT( xpn3**2 + mn2 ) + Un

  q3_lim = 0.0d0
  tmp1 = Up-En-Enu
  tmp2 = Enu**2 + ml2 - mp2 - xpn3**2

  xA = 4.d0*( tmp1**2-(Enu+xpn3)**2 )
  xB = 4.d0*( (tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*xpn3 )
  xC = 4.d0*tmp1**2*(xpn3**2+mp2) - (tmp2-tmp1**2)**2
  q3_lim(1) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  q3_lim(2) = ABS( (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )

  xA = 4.d0*( tmp1**2-(Enu-xpn3)**2 )
  xB = 4.d0*( -(tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*xpn3 )
  IF( xB**2 - 4.0d0*xA*xC >= 0 ) THEN
     q3_lim(3) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
     q3_lim(4) = ABS( (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  END IF

  xq0_min = 1.d10
  xq0_max = -1.d10
  DO i=1,4
    xq3 = q3_lim(i)
    tmp1 = Enu + SQRT( (Enu+xq3)**2 + ml2 )
    tmp2 = Enu + SQRT( (Enu-xq3)**2 + ml2 )
    tmp3 = SQRT( (xpn3+xq3)**2+mp2) + Up - En
    tmp4 = SQRT( (xpn3-xq3)**2+mp2) + Up - En
    IF( ABS((tmp1-tmp3)/tmp3)<1.d-6 .or. ABS((tmp1-tmp4)/tmp4)<1.d-6 ) THEN
      xq0 = tmp1
    END IF
    IF( ABS((tmp2-tmp3)/tmp3)<1.d-6 .or. ABS((tmp2-tmp4)/tmp4)<1.d-6 ) THEN
      xq0 = tmp2
    END IF
    xq0_min = min( xq0,  xq0_min)
    xq0_max = max( xq0,  xq0_max)
  END DO
  IF( ABS(xq0_max-xq0_min)<1.d-15 ) THEN
    IF(opt0.eq.3) THEN
      xq0_max = max(Enu + Mul + 50.0d0*T, xq0_min + 50.0d0*T)
    ELSE IF(opt0.eq.4) THEN
      xq0_max = max(xq0_min+50.0d0*T, SQRT(ppmax**2+mp2)+Up-En)
    END IF
  END IF

  xq0_l = Enu + lepton_mass
  tmp1 = SQRT( (xpn3+Enu)**2+mp2) + Up - En
  tmp2 = SQRT( (xpn3-Enu)**2+mp2) + Up - En
  IF( xq0_l <= tmp1 .and. xq0_l >= tmp2) THEN
     xq0_min = min( xq0_l, xq0_min)
  END IF

  xq0_h = proton_mass + Up - En
  tmp1 = Enu + SQRT( (Enu+xpn3)**2 + ml2 )
  tmp2 = Enu + SQRT( (Enu-xpn3)**2 + ml2 )
  IF( xq0_h <= tmp1 .and. xq0_h >= tmp2) THEN
     xq0_min = min( xq0_h, xq0_min)
  END IF

END SUBROUTINE Ebounds_D

SUBROUTINE Range_pn_D()

  REAL(DP) :: tmp,mpt,Epr,Esq,Equ,p2a,p2b,p20
  REAL(DP) :: F0, Fmax

  pnmin = 0.0d0
  ppmax = SQRT( (Tfac*T + Mup - Up)**2 - mp2 )
  SELECT CASE(opt0)
  CASE (3)
    pnmax = SQRT( (Tfac*T + Mun - Un)**2 - mn2 )
  CASE (4)
    pnmin = SQRT( (max(neutron_mass,Mun-Tfac*T-Un))**2 - mn2 )
    pnmax = SQRT( (SQRT(ppmax**2+mp2)+Up-Un)**2-mn2 )
  END SELECT

  mpt = proton_mass - lepton_mass
  IF( mpt > neutron_mass) THEN
    tmp = neutron_mass/(mpt-neutron_mass)*Enu
  ELSE
    tmp = pnmax
  END IF
  Fmax = SQRT( (tmp+Enu)**2+mpt**2 )+Up-SQRT(tmp**2+mn2)-Un

  IF(Enu>=Fmax) THEN
    pnmax = -1000.0d0
    RETURN
  END IF

  F0 = SQRT(Enu**2 + mpt**2)+Up-neutron_mass-Un
  IF( F0 >= Enu .and. Up >= Un) RETURN

  Epr = Enu - Up + Un
  Esq = Enu**2 + mpt**2 - mn2 - Epr**2
  Equ = Esq**2 + 4.d0*( (neutron_mass*Enu)**2 - (Epr*neutron_mass)**2 )
  IF(Equ>0.0d0.and. Enu**2 .ne. Epr**2) THEN
    p2a = -( Enu*Esq - ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    p2b = -( Enu*Esq + ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    IF(Epr<0.d0 .and. p2b>0.d0) THEN
      pnmin=p2b
    ELSE IF(Epr>=0.0d0 .and. p2a>0.d0) THEN
      pnmin=p2a
    END IF

  ELSE IF(Enu**2 .eq. Epr**2) THEN
    p20 = (4.0d0*Epr**2*mn2-Esq**2)/(4.0d0*Esq*Enu)
    IF(p20>0.d0) pnmin=p20
  END IF

END SUBROUTINE Range_pn_D

!=======================================================================
!
!     NUMERICAL RECIPES: Gauss-Legendre integration
!
!=======================================================================
SUBROUTINE gauleg(x1, x2, x, w, n)

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: x1, x2
  REAL(DP), INTENT(OUT), dimension(n) :: x, w
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
  INTEGER :: i,j,m
  REAL(DP), PARAMETER :: eps=3.d-14
  REAL(DP) :: p1,p2,p3,pp,xl,xm,z,z1

  m = (n+1)/2
  xm = 0.5d0*(x2 + x1)
  xl = 0.5d0*(x2 - x1)

  DO i=1,m
    z = COS(pi*(i - 0.25d0)/(n + 0.5d0))
    DO
      p1 = 1.d0
      p2 = 0.d0
      DO j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.d0*j - 1.d0)*z*p2 - (j - 1.d0)*p3)/j
      END DO
      pp = n * (z * p1 - p2) / (z*z - 1.d0)
      z1 = z
      z = z1 - p1/pp
      IF (ABS(z - z1) <= eps) EXIT
    END DO
    x(i)        = xm - xl*z
    x(n + 1 - i)= xm + xl*z
    w(i)        = 2.d0 * xl / ((1.d0 - z*z) * pp * pp)
    w(n + 1 - i)= w(i)
  END DO
  
END SUBROUTINE gauleg

END MODULE wlSemiLeptonicOpacityModule2D_old

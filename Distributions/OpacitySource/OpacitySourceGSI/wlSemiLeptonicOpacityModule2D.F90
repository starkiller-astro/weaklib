!!! modified date: June 17, 2020
!!! Authors: Codes wirtten by Gang Guo, based on the formalism by
!!! Andreas Lohs & Gang Guo
!!! neutrino CC opacity is given by 2D integrals
!!! processes to be considered in this routine
!!! ReactionIndex=1: v + n -> e- + p  (oapcity in 1/km)
!!! ReactionIndex=2: vb + p -> e+ + n (opacity in 1/km)
!!! ReactionIndex=3: vb + p + e- -> n (opacity in 1/km)
!!! ReactionIndex=4: n -> vb + e- + p [emissivity in v/(s cm^3 MeV)]

!!! nE defined in module: the No. of grids used in gaussian
!!! quadrature for each dimention
!!! seems nE needs be >~30 to reach an accuracy of 5%

MODULE wlSemiLeptonicOpacityModule2D

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mpi, Vud, &
   massA, massV, gamma_p, gamma_n, hbarc, hbar
  ! USE wlEosConstantsModule, ONLY: &
  !  pi, hbarc, hbar

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PARAMETER :: F2wm0 = gamma_p - gamma_n - 1.0d0
  ! ! If you want to reproduce exactly the GSI numbers you need:
  ! REAL(DP),  PARAMETER :: Gw_MeV=1.166d-11,Vud=0.97427d0, F2wm0 = 3.706d0, &
  !      ga =1.2723d0,gv=1.d0, Mpi=139.57d0, Mnp=938.919d0, &
  !      Dnp=1.293d0, massA=1.0d3,massV=840.d0 , GfVud2 = (Gw_MeV*Vud)**2

  REAL(DP), PARAMETER :: Tfac = 150.0d0
  REAL(DP), PARAMETER :: gasq = ga**2, gvsq = gv**2, gva = gv*ga

  !========================================================================
! Routine to calculate the neutrino CC opacity via 2D integrals
! Inputs:
! (1) WhichCorrection=0,1,2,3 for LO, LO+weak magnetism(WM),
! LO+WM+pseduoscalar(PS), LO+WM+PS+form factor dependencies
! (2) ReactionIndex=ReactionIndex=1,2,3,4 for nu capture on n, nub capture on p, inverse
! neutron decay, neutron decay emissivity
! (3) Enu(NP): neutrino energy array in MeV with dimension NP
! (4) xT: Temperature in MeV
! (5) xMul, xml: relativistic chemical potentials & masses of
! e^-(mu^-) in MeV
! (6) xMun,xMup,xmn,xmp,xUn,xUp: relativistic chemical
! potentials, masses, & potentials of neutron and proton in MeV
! relativsitic disperstion relation for nucleons,
! E2=SQRT(p^2+Mass2^2)+U2
! Outputs:
! (7) OpaA(NP): opacities in km^-1; for neutron decay, emissivities
! in (s cm^3 MeV)^-1
!     **** blocking factor of final state neutrinos in neutron decay
!     is not considered ****
!=============================================================
  PUBLIC :: Opacity_CC_2D

CONTAINS

SUBROUTINE Opacity_CC_2D(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp, nE) 

  INTEGER , INTENT(IN)  :: ReactionIndex, WhichCorrection, nE
  REAL(DP), INTENT(IN)  :: Enu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp 
  REAL(DP), INTENT(OUT) :: OpaA

  REAL(DP) :: anti
  INTEGER  :: IncludeCorrections(3)

  SELECT CASE(WhichCorrection)
  CASE (0)
     IncludeCorrections(1) = 0
     IncludeCorrections(2) = 0
     IncludeCorrections(3) = 0
  CASE (1)
     IncludeCorrections(1) = 1
     IncludeCorrections(2) = 0
     IncludeCorrections(3) = 0
  CASE (2)
     IncludeCorrections(1) = 1
     IncludeCorrections(2) = 1
     IncludeCorrections(3) = 0
  CASE (3)
     IncludeCorrections(1) = 1
     IncludeCorrections(2) = 1
     IncludeCorrections(3) = 1
  END SELECT

  SELECT CASE(ReactionIndex)
  CASE (1)
    anti = 1.0d0
    CALL Integral_2D(Enu, xT, xmn, xml, xmp, xUn, xUp, xMun, xMul, xMup, anti, 1, IncludeCorrections, nE, OpaA)
  CASE (2)
    anti = -1.0d0
    CALL Integral_2D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, -xMul, xMun, anti, 2, IncludeCorrections, nE, OpaA)
  CASE (3)
    anti = -1.0d0
    CALL Integral_2D_D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 3, IncludeCorrections, nE, OpaA)
  CASE (4)
    anti = -1.0d0
    CALL Integral_2D_D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 4, IncludeCorrections, nE, OpaA)
  END SELECT

END SUBROUTINE Opacity_CC_2D

! For a more efficient parallelization this whole thing can be probably taken out, and we make 4 different integrations, one per reaction
SUBROUTINE Integral_2D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, IncludeCorrections, nE, res)   ! captures

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE
  REAL(DP), INTENT(OUT) :: res

  REAL(DP) :: P2,xq0_min,xq0_max,E3_min,E3_max,E2_min,E2_max
  REAL(DP) :: xa(nE), wxa(nE), wE2(nE), wE3(nE)
  INTEGER  :: i, j
  REAL(DP) :: xamp, E2, E3, E4
  REAL(DP) :: xf2,xf3,xf4
  REAL(DP) :: p2min, p2max

  res = 0.0d0
  CALL Range_pn(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, p2min, p2max)
  IF( p2min > p2max ) RETURN
  E2_min = SQRT(p2min**2 + Mass2**2) + U2
  E2_max = SQRT(p2max**2 + Mass2**2) + U2

  CALL gauleg(0.d0, 1.d0, xa, wxa, nE)

  ! Loop over incoming particle (i.e. particle 2)
  !!$omp parallel do private(i,j,E2,P2,xf2,E3,E3,xf3,xf4,xamp,E3_min,E3_max,xq0_min,xq0_max) reduction(+:res)
  DO i=1, nE
    E2 = xa(i)*(E2_max - E2_min) + E2_min
    wE2(i) = wxa(i)*(E2_max - E2_min)
    P2 = SQRT( (E2 - U2)**2 - Mass2**2 )
    CALL Ebounds( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, xq0_min, xq0_max )
    E3_min = Enu - xq0_max
    E3_max = Enu - xq0_min
    xf2 = FD(E2, Mu2, T)
  ! Loop over outgoing particle (i.e. particle 3)
    DO j=1,nE
      E3 = xa(j)*(E3_max - E3_min) + E3_min
      wE3(j) = wxa(j)*(E3_max - E3_min)
      E4 = Enu + E2 - E3
      CALL Calc_ampsq( Enu, T, E2, E3, Mass2, Mass3, Mass4, U2, U4, anti, ReactionIndex, IncludeCorrections, xamp)
      xf3 = FD(E3, Mu3, T)
      xf4 = FD(E4, Mu4, T)
      res = res + xf2*(1.0d0-xf3)*(1.0d0-xf4)*xamp*wE2(i)*wE3(j)/Enu**2
    END DO
  END DO
  !!$omp end parallel do

  res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc

END SUBROUTINE Integral_2D

! For a more efficient parallelization this whole thing can be probably taken out, and we make 4 different integrations, one per reaction
SUBROUTINE Integral_2D_D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, IncludeCorrections, nE, res)   ! decay/inverse decay

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE
  REAL(DP), INTENT(OUT) :: res

  REAL(DP) :: P2,xq0_min,xq0_max,E3_min,E3_max,E2_min,E2_max
  REAL(DP) :: xa(nE), wxa(nE), wE2(nE), wE3(nE)
  INTEGER :: i,j
  REAL(DP) :: xamp
  REAL(DP) :: E2, E3, E4
  REAL(DP) :: xf2,xf3,xf4
  REAL(DP) :: p2min, p2max, p4max

  res = 0.0d0
  CALL Range_pn_D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4, ReactionIndex, p2min, p2max, p4max)
  IF( p2min > p2max ) RETURN
  E2_min = SQRT(p2min**2 + Mass2**2) + U2
  E2_max = SQRT(p2max**2 + Mass2**2) + U2

  CALL gauleg(0.d0, 1.d0, xa, wxa, nE)

  ! Loop over incoming particle (i.e. particle 2)
  !!$omp parallel do private(i,j,E2,P2,xf2,E3,E3,xf3,xf4,xamp,E3_min,E3_max,xq0_min,xq0_max) reduction(+:res)
  DO i=1, nE
    E2 = xa(i)*(E2_max - E2_min) + E2_min
    wE2(i) = wxa(i)*(E2_max - E2_min)
    P2 = SQRT( (E2-U2)**2 - Mass2**2 )
    CALL Ebounds_D( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max, ReactionIndex, xq0_min, xq0_max )
    E3_min = xq0_min - Enu
    E3_max = xq0_max - Enu
    xf2 = FD(E2, Mu2, T)
  ! Loop over outgoing particle (i.e. particle 3)
    DO j=1,nE
      E3 = xa(j)*(E3_max - E3_min) + E3_min
      wE3(j) = wxa(j)*(E3_max - E3_min)
      E4 = Enu + E3 + E2
      CALL Calc_ampsq(Enu, T, E2, -E3, Mass2, Mass3, Mass4, U2, U4, anti, ReactionIndex, IncludeCorrections, xamp)
      xf3 = FD(E3, Mu3, T)
      xf4 = FD(E4, Mu4, T)
      IF(ReactionIndex.eq.3) THEN ! inverse decay
        res = res + xf2*xf3*(1.0d0-xf4)*xamp*wE2(i)*wE3(j)/Enu**2
      ELSE IF(ReactionIndex.eq.4) THEN ! n decay, *blocking of neutrinos is not added*
        res = res + (1.0d0-xf2)*(1.0d0-xf3)*xf4*xamp*wE2(i)*wE3(j)/Enu**2
      END IF
    END DO
  END DO
  !!$omp end parallel do

  IF(ReactionIndex.eq.3) THEN
    res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
  ELSE IF(ReactionIndex.eq.4) THEN
    res = res*(Gw_MeV*Vud)**2/32.0d0/(pi**7)*Enu**2/hbarc**3/hbar
  END IF

  res = -res

END SUBROUTINE Integral_2D_D


!!! calculate the amplitudes, keep the following unchanged
SUBROUTINE Calc_Ampsq( Enu, T, E2, E3, Mass2, Mass3, Mass4, U2, U4, anti, ReactionIndex, IncludeCorrections, Ampsq )

  REAL(DP), INTENT(IN)  :: Enu, T, E2, E3, Mass2, Mass3, Mass4, U2, U4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3)
  REAL(DP), INTENT(OUT) :: Ampsq

  REAL(DP) :: P1,P2,P3,P4,Pmax,Pmin,Pmax2,Pmin2,Pmax3,Pmin3
  REAL(DP) :: Ia,Ib,Ic,Id,Ie,Ig,Ih,Ij,Ik,Il,A,B,Ct,D,E,del0,&
       del1,del2,eps0,eps1,alp0,alp1,alp2,bet0,bet1,cplA,cplB,cplC,&
       cplD,cplE,cplG,cplH,cplJ,cplK,cplL,JFF,KFF,LFF 
  REAL(DP) :: E2f,E4f,E1,Qmass,dmf,dU2,E4,MassNuc
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

  MassNuc = 0.5d0*( Mass2 + Mass4 )
  dU2 = U2 - U4
  dmf = Mass2 - Mass4
  Qmass = 0.5d0*( Mass2**2 - Mass4**2 )

  E1 = Enu
  E4 = E1 + E2 - E3
  E4f = E4 - U4
  E2f = E2 - U2
  P1 = E1
  P3 = SQRT( E3**2 - Mass3**2 )
  P2 = SQRT( (E2 - U2)**2 - Mass2**2 )
  P4 = SQRT( (E4 - U4)**2 - Mass4**2 )

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

  IF(IncludeCorrections(1) .eq. 1) THEN
     gvf = gv*F2wm0
     gaf = ga*F2wm0
     f2sq = F2wm0**2
  ELSE
     gvf = 0.0d0
     gaf = 0.0d0
     f2sq = 0.0d0
  END IF

  cplA = (gvsq+2.0d0*anti*gva+gasq)+anti*2.0d0*gaf*Mass2/MassNuc*(1.0d0-dmf/2.0d0/Mass2)
  cplB = (gvsq-2.0d0*anti*gva+gasq)-anti*2.0d0*gaf*Mass2/MassNuc*(1.0d0-dmf/2.0d0/Mass2)
  cplC = F2sq/MassNuc**2
  cplD =-F2sq/MassNuc**2
  cplE = F2sq/MassNuc**2*(-0.5d0*Mass3**2+dU2*(E3-E1)-0.5d0*dU2**2)
  cplG = gvf*Mass2/MassNuc*(2.0d0-dmf/Mass2)+0.5d0*F2sq/MassNuc**2*(Mass2*Mass4-Qmass+0.25d0*&
       Mass3**2-dU2*(E1+E2f)-dU2**2/4.0d0) 
  cplH = 0.5d0*F2sq/MassNuc**2*(2.0d0*Qmass+Mass3**2+dU2*(3.0d0*E1-E3+2.0d0*E4f))

  JFF = -Mass3**2*(E1+0.5d0*E2f+0.5d0*E4f)+Qmass*(E3-3.0d0*E1) &
       +0.5d0*dU2*(E4f*(3.0d0*E3-5.0d0*E1)+E2f*(E3+E1)+E3**2-E1**2-2.0d0*Qmass) &
       + dU2**2*(E1-E3-E4f)+0.5d0*dU2**3
  JFF = JFF*dU2
  cplJ = gvf*dmf/MassNuc*0.5d0*(Mass3**2-dU2*(E1+E3))+F2sq/MassNuc**2*0.5d0*JFF

  kFF = -(Mass2+3.0d0*Mass4)*Mass2*0.25d0*Mass3**2 + Qmass**2 + Qmass*0.25d0*Mass3**2-Mass3**4/8.0d0 &
       +dU2*(0.5d0*Qmass*(3.0d0*E1-E2f+E3+3.0d0*E4f)+0.25d0*Mass3**2*(2.0d0*E2f+E3+&
       E1)+Mass2*Mass4*(E3-E1)) & 
       +dU2**2*( 0.25d0*(Mass2**2-Mass2*Mass4-3.0d0*Qmass)+E4f*0.5d0*(2.0d0*E1-E2f+&
       E3+E4f)+E2f*(0.5d0*E1-E3)+0.5d0*E1**2) & 
       +dU2**3*0.25d0*( -E1+2.0d0*E2f-E3-2.0d0*E4f ) + 0.125d0*dU2**4
  cplK = (gasq-gvsq)*Mass2*Mass4+gvf*Mass2/MassNuc*0.5d0*(-3.0d0*Mass3**2+4.0d0*dU2*(E3-&
       E1)-dU2**2 & 
       +dmf/Mass2*(2.0d0*Qmass+Mass3**2+dU2*(E4f+2.0d0*E1-E3)))&
       +0.5d0*F2sq/MassNuc**2*KFF

  LFF = Mass3**2*(Mass2+Mass4)**2-4.0d0*Qmass**2+dU2*(-Mass3**2*(E2f+E4f+E1)-2.0d0*&
       E3*(Mass2**2+Mass2*Mass4) & 
       +2.0d0*Qmass*(E2f-3.0d0*E4f-E1)) &
       +2.0d0*dU2**2*( E2f*E3+E4f*(E2f-E4f-E1) + Qmass ) + dU2**3*(-E2f+E4f+E1)
  LFF = LFF*dU2*E1*0.25d0
  cplL = gvf*Mass2/MassNuc*dU2*E1*(Mass3**2-dU2*E3-0.5d0*dmf/Mass2*(Qmass+0.5d0*Mass3**&
       2+dU2*E4f-0.5d0*dU2**2))   & 
       +0.5d0*F2sq/MassNuc**2*LFF

  Ampsq = max(0.0d0,cplA*Ia+cplB*Ib+cplC*Ic+cplD*Id+cplE*Ie+cplG*Ig+cplH*Ih+&
       cplJ*Ij+cplK*Ik+cplL*Il) 

!!! adding form-factor dependences
  AmassSq = massA**2
  APmassSq = 2.0d0*AmassSq
  tm3sq = Mass3**2 + 2.0d0*dU2*(E1-E3) + dU2**2
  IF(IncludeCorrections(3) .eq. 1) THEN
    VmassSq = 1.0d0/(1.0d0/massV**2-F2wm0/MassNuc**2/8.0d0)
    AVmassSq = 1.0d0/(0.5d0/massA**2+0.5d0/massV**2-F2wm0/MassNuc**2/16.0d0)
    AFmassSq = 1.0d0/(0.5d0/massA**2+0.5d0/massV**2+1.0d0/MassNuc**2/16.0d0)
    VFmassSq = 1.0d0/(1.0d0/massV**2-(F2wm0-1.0d0)/16.0d0/MassNuc**2)
    FmassSq = 1.0d0/(1.0d0/massV**2+1.0d0/8.0d0/MassNuc**2)

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
        gaf/AFmassSq*Mass2/MassNuc*(1.0d0-dmf/2.0d0/Mass2) 
    cplB2 = gvsq/VmassSq+gasq/AmassSq-2.0d0*anti*gva/AVmassSq-anti*2.0d0*&
        gaf/AFmassSq*Mass2/MassNuc*(1.0d0-dmf/2.0d0/Mass2) 
    cplC2 = F2sq/MassNuc**2/FmassSq
    cplD2 =-F2sq/MassNuc**2/FmassSq
    cplE2 = F2sq/FmassSq/MassNuc**2*(-0.5d0*Mass3**2+dU2*(E3-E1)-0.5d0*dU2**2)
    cplG2 = gvf/VFmassSq*Mass2/MassNuc*(2.0d0-dmf/Mass2)+0.5d0*F2sq/FmassSq/MassNuc**2*&
        (Mass2*Mass4-Qmass+0.25d0*Mass3**2-dU2*(E1+E2f)-dU2**2/4.0d0) 
    cplH2 = 0.5d0*F2sq/FmassSq/MassNuc**2*(2.0d0*Qmass+Mass3**2+dU2*(3.0d0*E1-E3+&
        2.0d0*E4f)) 
    cplJ2 = gvf/VFmassSq*dmf/MassNuc*0.5d0*(Mass3**2-dU2*(E1+E3))+F2sq/&
        FmassSq/MassNuc**2*0.5d0*JFF 
    cplK2 = (gasq/AmassSq-gvsq/VmassSq)*Mass2*Mass4+gvf/VFmassSq*Mass2/MassNuc*&
        0.5d0*(-3.0d0*Mass3**2+4.0d0*dU2*(E3-E1)-dU2**2+dmf/Mass2*(2.0d0*Qmass+&
        Mass3**2+dU2*(E4f+2.0d0*E1-E3)))+0.5d0*F2sq/MassNuc**2*KFF/FmassSq 
    cplL2 = gvf/VFmassSq*Mass2/MassNuc*dU2*E1*(Mass3**2-dU2*E3-0.5d0*dmf/Mass2*&
        (Qmass+0.5d0*Mass3**2+dU2*E4f-0.5d0*dU2**2))+0.5d0*F2sq/MassNuc**2*LFF/&
        FmassSq 

    Ampsq_ff0 = cplA2*Ia2 + cplB2*Ib2 + cplC2*Ic2 + cplD2*Id2 + &
        cplE2*Ie2 + cplG2*Ig2 + cplH2*Ih2 + cplJ2*Ij2 + cplK2*Ik2 + &
        cplL2*Il2 
     IF(Ampsq+ampsq_ff0>=0.0d0) Ampsq=Ampsq+ampsq_ff0           
  END IF

!!! adding pseudoscalar terms w/o form factor dependence
  IF(IncludeCorrections(2) .eq. 1) THEN

    YAP = 0.25d0*AmassSq+tm3sq-2.0d0*E1*E3+P1**2+P3**2
    ZAP = mpi**2-tm3sq+2.0d0*E1*E3-P1**2-P3**2
    MassMy = 0.25d0*AmassSq+mpi**2

    IF(IncludeCorrections(3) .eq. 0) THEN
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

    ELSE IF(IncludeCorrections(3) .eq. 1) THEN
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

      cplLAP = 2.0d0*MassNuc*gasq*dU2*E1*(2.0d0*Mass2*Mass3**2-dmf*(Qmass+0.5d0*Mass3**2)-&
            dU2*(2.0d0*E3*Mass2+E4f*dmf)+0.5d0*dU2**2*dmf) 
      cplKAP = 2.0d0*MassNuc*gasq*(Mass2*(dU2**2-Mass3**2)+dmf*dU2*(E2f+E1))
      cplJAP = 2.0d0*MassNuc*gasq*dmf*(Mass3**2-dU2*(E1+E3))

      cplLPP = 2.0d0*MassNuc**2*gasq*dU2*E1*(Mass3**4-Mass3**2*dmf**2+dU2*(-Mass3**2*&
            (3.0d0*E3-2.0d0*E1)+dmf**2*E3)+dU2**2*(Mass3**2+2.0d0*E3*(E3-E1))-&
            dU2**3*E3) 
      cplKPP = 2.0d0*MassNuc**2*gasq*(-0.5d0*Mass3**4+0.5d0*dmf**2*Mass3**2+dU2*Mass3**2*&
            (E3-3.0d0*E1)+dU2**2*(2.0d0*E1*E3-0.5d0*dmf**2)+dU2**3*(E4f-E2f)-&
            0.5d0*dU2**4)  
      cplGPP = 2.0d0*MassNuc**2*gasq*(Mass3**2-dU2**2)

      ampsq_pseudo = cplLAP*IlAP+cplKAP*IkAP+cplJAP*IjAP + cplLPP*IlPP +&
            cplKPP*IkPP+cplGPP*IgPP 
      IF(Ampsq+ampsq_pseudo>=0.0d0) Ampsq=Ampsq+ampsq_pseudo
    END IF

END subroutine Calc_Ampsq

SUBROUTINE Ebounds( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, xq0_min, xq0_max )

  REAL(DP), INTENT(INOUT) :: E2
  REAL(DP), INTENT(IN)  :: Enu, T, P2, Mass2, Mass3, Mass4, U2, U4
  REAL(DP), INTENT(OUT) :: xq0_min, xq0_max

  REAL(DP) :: q3_lim(4), tmp1, tmp2, xA, xB, xC, xq3
  INTEGER :: i
  REAL(DP) :: tmp3, tmp4, xq0, xq0_l, xq0_h

  E2 = SQRT( P2**2 + Mass2**2 ) + U2

  q3_lim = 0.0d0
  tmp1 = U4-E2-Enu
  tmp2 = Enu**2 + Mass3**2 - Mass4**2 - P2**2

  xA = 4.d0*( tmp1**2-(Enu-P2)**2 )
  xB = 4.d0*( (tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*P2 )
  xC = 4.d0*tmp1**2*(P2**2+Mass4**2) - (tmp2-tmp1**2)**2
  IF( xB**2 - 4.0d0*xA*xC >= 0.0d0) THEN
    q3_lim(1) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
    q3_lim(2) = ABS( (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  END IF

  xA = 4.d0*( tmp1**2-(Enu+P2)**2 )
  xB = 4.d0*( -(tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*P2 )
  IF( xB**2 - 4.0d0*xA*xC >= 0.0d0) THEN
    q3_lim(3) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
    q3_lim(4) = ABS( (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  END IF

  xq0_min = 1.d30
  xq0_max = -1.d30
  xq0 = 0.0d0
  DO i=1,4
    xq3 = q3_lim(i)
    tmp1 = Enu - SQRT( (Enu-xq3)**2 + Mass3**2 )
    tmp2 = Enu - SQRT( (Enu+xq3)**2 + Mass3**2 )
    tmp3 = SQRT( (P2+xq3)**2+Mass4**2) + U4 - E2
    tmp4 = SQRT( (P2-xq3)**2+Mass4**2) + U4 - E2

    IF( ABS((tmp1-tmp3)/tmp3)<1.d-6 .or. ABS((tmp1-tmp4)/tmp4)<1.d-6 ) THEN
      xq0 = tmp1
    END IF
    IF( ABS((tmp2-tmp3)/tmp3)<1.d-6 .or. ABS((tmp2-tmp4)/tmp4)<1.d-6 ) THEN
      xq0 = tmp2
    END IF
    xq0_min = min( xq0,  xq0_min)
    xq0_max = max( xq0,  xq0_max)
  END DO

  xq0_l = Enu - Mass3
  tmp1 = SQRT( (P2+Enu)**2+Mass4**2) + U4 - E2
  tmp2 = SQRT( (P2-Enu)**2+Mass4**2) + U4 - E2
  IF( xq0_l <= tmp1 .and. xq0_l >= tmp2) THEN
     xq0_max = max( xq0_l, xq0_max)
  END IF

  xq0_h = Mass4 + U4 - E2
  tmp1 = Enu - SQRT( (Enu-P2)**2 + Mass3**2 )
  tmp2 = Enu - SQRT( (Enu+P2)**2 + Mass3**2 )
  IF( xq0_h <= tmp1 .and. xq0_h >= tmp2) THEN
     xq0_min = min( xq0_h, xq0_min)
  END IF

END SUBROUTINE Ebounds

SUBROUTINE Range_pn( Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, p2min, p2max )

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2
  REAL(DP), INTENT(OUT) :: p2min, p2max

  REAL(DP) :: tmp,mpt,Epr,Esq,Equ,p2a,p2b,p20,F0,Finf,Fmax,Fmin

  p2min = 0.0d0
  p2max = SQRT( (Tfac*T + Mu2 - U2)**2 - Mass2**2 )

  mpt = Mass4 + Mass3
  F0 = SQRT(Enu**2 + mpt**2) + U4 - Mass2 - U2
  Finf = -Enu + U4 - U2
  IF(mpt < Mass2) THEN
    tmp = (1.d0+mpt/Mass2)/(1.d0-mpt**2/Mass2**2)*Enu
    Fmin = SQRT( (tmp-Enu)**2+mpt**2 )+U4-SQRT(tmp**2+Mass2**2)-U2
    Fmax = max( F0, Finf )
  ELSE
    Fmin = Finf
    Fmax = F0
  END IF

  IF( Fmax <= Enu) RETURN
  IF( Fmin >= Enu) THEN
    p2max = -1000.0d0
    RETURN
  END IF

  Epr = Enu - U4 + U2
  Esq = Enu**2 + mpt**2 - Mass2**2 - Epr**2
  Equ = Esq**2 + 4.d0*Mass2**2*( Enu**2-Epr**2 )
  IF(Equ>0.0d0 .and. Enu**2 .ne. Epr**2) THEN
    p2a = ( Enu*Esq - ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    p2b = ( Enu*Esq + ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    IF(Epr<0.d0 .and. p2b>0.d0) THEN
      p2min=p2b
    ELSE IF(Epr>=0.d0 .and. p2a>0.d0) THEN
      p2min=p2a
    END IF

  ELSE IF(Enu**2 .eq. Epr**2) THEN
    p20 = -(4.0d0*Epr**2*Mass2**2-Esq**2)/(4.0d0*Esq*Enu)
    IF(p20>0.d0) p2min=p20
  END IF

END SUBROUTINE Range_pn

SUBROUTINE Ebounds_D( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max, ReactionIndex, xq0_min, xq0_max )
  
  REAL(DP), INTENT(INOUT) :: E2
  REAL(DP), INTENT(IN)  :: Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max
  INTEGER,  INTENT(IN)  :: ReactionIndex
  REAL(DP), INTENT(OUT) :: xq0_min, xq0_max
  
  REAL(DP) :: q3_lim(4), tmp1, tmp2, xA, xB, xC, xq3
  INTEGER  :: i
  REAL(DP) :: tmp3, tmp4, xq0, xq0_l,xq0_h

  E2 = SQRT( P2**2 + Mass2**2 ) + U2

  q3_lim = 0.0d0
  tmp1 = U4-E2-Enu
  tmp2 = Enu**2 + Mass3**2 - Mass4**2 - P2**2

  xA = 4.d0*( tmp1**2-(Enu+P2)**2 )
  xB = 4.d0*( (tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*P2 )
  xC = 4.d0*tmp1**2*(P2**2+Mass4**2) - (tmp2-tmp1**2)**2
  q3_lim(1) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  q3_lim(2) = ABS( (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )

  xA = 4.d0*( tmp1**2-(Enu-P2)**2 )
  xB = 4.d0*( -(tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*P2 )
  IF( xB**2 - 4.0d0*xA*xC >= 0 ) THEN
     q3_lim(3) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
     q3_lim(4) = ABS( (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  END IF

  xq0_min = 1.d10
  xq0_max = -1.d10
  DO i=1,4
    xq3 = q3_lim(i)
    tmp1 = Enu + SQRT( (Enu+xq3)**2 + Mass3**2 )
    tmp2 = Enu + SQRT( (Enu-xq3)**2 + Mass3**2 )
    tmp3 = SQRT( (P2+xq3)**2+Mass4**2) + U4 - E2
    tmp4 = SQRT( (P2-xq3)**2+Mass4**2) + U4 - E2
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
    IF(ReactionIndex.eq.3) THEN
      xq0_max = max(Enu + Mu3 + 50.0d0*T, xq0_min + 50.0d0*T)
    ELSE IF(ReactionIndex.eq.4) THEN
      xq0_max = max(xq0_min+50.0d0*T, SQRT(p4max**2+Mass4**2)+U4-E2)
    END IF
  END IF

  xq0_l = Enu + Mass3
  tmp1 = SQRT( (P2+Enu)**2+Mass4**2) + U4 - E2
  tmp2 = SQRT( (P2-Enu)**2+Mass4**2) + U4 - E2
  IF( xq0_l <= tmp1 .and. xq0_l >= tmp2) THEN
     xq0_min = min( xq0_l, xq0_min)
  END IF

  xq0_h = Mass4 + U4 - E2
  tmp1 = Enu + SQRT( (Enu+P2)**2 + Mass3**2 )
  tmp2 = Enu + SQRT( (Enu-P2)**2 + Mass3**2 )
  IF( xq0_h <= tmp1 .and. xq0_h >= tmp2) THEN
     xq0_min = min( xq0_h, xq0_min)
  END IF

END SUBROUTINE Ebounds_D

SUBROUTINE Range_pn_D( Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4, ReactionIndex, p2min, p2max, p4max )

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4
  INTEGER , INTENT(IN)  :: ReactionIndex
  REAL(DP), INTENT(OUT) :: p2min, p2max, p4max

  REAL(DP) :: tmp,mpt,Epr,Esq,Equ,p2a,p2b,p20
  REAL(DP) :: F0, Fmax

  p2min = 0.0d0
  p4max = SQRT( (Tfac*T + Mu4 - U4)**2 - Mass4**2 )
  SELECT CASE(ReactionIndex)
  CASE (3)
    p2max = SQRT( (Tfac*T + Mu2 - U2)**2 - Mass2**2 )
  CASE (4)
    p2min = SQRT( (max(Mass2,Mu2-Tfac*T-U2))**2 - Mass2**2 )
    p2max = SQRT( (SQRT(p4max**2+Mass4**2)+U4-U2)**2-Mass2**2 )
  END SELECT

  mpt = Mass4 - Mass3
  IF( mpt > Mass2) THEN
    tmp = Mass2/(mpt-Mass2)*Enu
  ELSE
    tmp = p2max
  END IF
  Fmax = SQRT( (tmp+Enu)**2+mpt**2 )+U4-SQRT(tmp**2+Mass2**2)-U2

  IF(Enu>=Fmax) THEN
    p2max = -1000.0d0
    RETURN
  END IF

  F0 = SQRT(Enu**2 + mpt**2)+U4-Mass2-U2
  IF( F0 >= Enu .and. U4 >= U2) RETURN

  Epr = Enu - U4 + U2
  Esq = Enu**2 + mpt**2 - Mass2**2 - Epr**2
  Equ = Esq**2 + 4.d0*( (Mass2*Enu)**2 - (Epr*Mass2)**2 )
  IF(Equ>0.0d0.and. Enu**2 .ne. Epr**2) THEN
    p2a = -( Enu*Esq - ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    p2b = -( Enu*Esq + ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2-Epr**2)
    IF(Epr<0.d0 .and. p2b>0.d0) THEN
      p2min=p2b
    ELSE IF(Epr>=0.0d0 .and. p2a>0.d0) THEN
      p2min=p2a
    END IF

  ELSE IF(Enu**2 .eq. Epr**2) THEN
    p20 = (4.0d0*Epr**2*Mass2**2-Esq**2)/(4.0d0*Esq*Enu)
    IF(p20>0.d0) p2min=p20
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

PURE ELEMENTAL FUNCTION FD(E, mu, T)
  REAL(DP), INTENT(IN) :: E, mu, T
  REAL(DP) :: FD

  FD = 1.0d0 / (EXP((E - mu)/T) + 1.0d0)
END FUNCTION

END MODULE wlSemiLeptonicOpacityModule2D

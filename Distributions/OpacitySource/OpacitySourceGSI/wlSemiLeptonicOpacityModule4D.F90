!!! modified date: June 17, 2020
!!! Authors: Codes wirtten by Gang Guo
!!! neutrino CC opacity is given by 4D integrals using CUBA library
!!! processes to be considered in this routine
!!! opt=1: v + n -> e- + p  (oapcity in 1/km)
!!! opt=2: vb + p -> e+ + n (opacity in 1/km)
!!! opt=3: vb + p + e- -> n (opacity in 1/km)
!!! opt=4: n -> vb + e- + p [emissivity in v/(s cm^3 MeV)]

!!! globle modules for cuba high dimension integrals
!!! please adjust 'nstart_C' & 'nincrease_C' in the module
!!! to trade off between effiency & accuracy

MODULE wlSemiLeptonicOpacityModule4D

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, mpi, Vud, &
   massA, massV, gamma_p, gamma_n
   
  IMPLICIT NONE
  PRIVATE

  ! This part is for CUBA -----------------------------------------
  ! Integration settings
  INTEGER,  PARAMETER :: ndim_C     = 4
  INTEGER,  PARAMETER :: ncomp_C    = 1
  INTEGER,  PARAMETER :: nvec_C     = 1
  INTEGER,  PARAMETER :: last_C     = 4
  INTEGER,  PARAMETER :: seed_C     = 0
  INTEGER,  PARAMETER :: mineval_C  = 0
  INTEGER,  PARAMETER :: maxeval_C  = 1000000
  REAL(DP), PARAMETER :: epsrel_C   = 1.0d-3
  REAL(DP), PARAMETER :: epsabs_C   = 1.0d-5
  REAL(DP), PARAMETER :: userdata_C = 0.0d0

  ! Sampling configuration
  INTEGER,  PARAMETER :: nstart_C    = 10000     ! Tunable
  INTEGER,  PARAMETER :: nincrease_C = 1000      ! Tunable
  INTEGER,  PARAMETER :: nbatch_C    = 1000
  INTEGER,  PARAMETER :: gridno_C    = 0

  character(len=*), PARAMETER :: statefile_C = ""
  INTEGER(KIND=8),  PARAMETER :: spin_C = -1

  ! Suave-specific settings (IF used)
  INTEGER,  PARAMETER :: nnew_C      = 300000   ! More new points per cycle not used by vegas
  INTEGER,  PARAMETER :: nmin_C      = nnew_C / 2
  REAL(DP), PARAMETER :: flatness_C  = 3.0d0  ! Lower is better for sharp features not used by vegas

  ! Unused, only for Suave
  ! Cuhre-specific configuration 
  INTEGER,  PARAMETER :: key1_C         = 47 
  INTEGER,  PARAMETER :: key2_C         = 1
  INTEGER,  PARAMETER :: key3_C         = 1
  INTEGER,  PARAMETER :: maxpass_C      = 5
  REAL(DP), PARAMETER :: border_C       = 1.0d-2
  REAL(DP), PARAMETER :: maxchisq_C     = 10.0d0
  REAL(DP), PARAMETER :: mindeviation_C = 0.25d0
  INTEGER,  PARAMETER :: ngiven_C       = 0
  INTEGER,  PARAMETER :: ldxgiven_C     = ndim_C
  INTEGER,  PARAMETER :: nextra_C       = 0

  ! Vegas key_C
  INTEGER,  PARAMETER :: key_C = 9
  INTEGER,  PARAMETER :: verbose_C = 0
  ! This part is for CUBA -----------------------------------------

  ! These are constants  -----------------------------------------
  REAL(DP),  PARAMETER :: Mnp = (mp + mn)/2.0d0, F2wm0 = gamma_p - gamma_n - 1.0d0, &
       Dnp = mn - mp, GfVud2 = (Gw_MeV*Vud)**2
  ! ! If you want to reproduce exactly the GSI numbers you need:
  ! REAL(DP),  PARAMETER :: Gw_MeV=1.166d-11,Vud=0.97427d0, F2wm0 = 3.706d0, &
  !      ga =1.2723d0,gv=1.d0, Mpi=139.57d0, Mnp=938.919d0, &
  !      massA = 1.0d3, massV = 840.d0, &
  !      Dnp=1.293d0, massA=1.0d3,massV=840.d0 , GfVud2 = (Gw_MeV*Vud)**2
  REAL(DP),  PARAMETER :: Tfac = 100.0d0

  REAL(DP) :: jVA = 1.d0, JAF = 1.d0 ! They get modified later but only for cases 3 and 4, we'll see in the future
  REAL(DP):: pe(0:3), pv(0:3), pn(0:3), pp(0:3), q(0:3), &
    Enu, pn3, En, cthe, q3, Jacob

  INTEGER  :: opt_form, opt_pseudo, opt_wm, reaction_index

  ! These are the global variables needed by CUBA
  REAL(DP) :: T, Mup, Mul, Mun, Un, Up, proton_mass, neutron_mass, lepton_mass
  REAL(DP) :: pnmin, pnmax, ppmax, invT, mn2, mp2

#if defined(USE_OMP)
  !$omp threadprivate(pe, pv, pn, pp, q, Enu, pn3, En, cthe, q3, Jacob,  &
  !$omp               T, Mup, Mul, Mun, Un, Up,                          &
  !$omp               proton_mass, neutron_mass, lepton_mass,            &
  !$omp               pnmin, pnmax, ppmax, invT, mn2, mp2,               &
  !$omp               opt_form, opt_pseudo, opt_wm, reaction_index, jVA, jAF)
#endif

  PUBLIC :: Opacity_CC_4D

!========================================================================
! Routine to calculate the neutrino CC opacity via 4D integrals
! (same inputs & output as the 2D SUBROUTINE)
! Inputs:
! (1) opt_had=0,1,2,3 for LO, LO+weak magnetism(WM),
! LO+WM+pseduoscalar(PS), LO+WM+PS+form factor dependencies
! (2) opt=reaction_index=1,2,3,4 for nu capture on n, nub capture on p, inverse
! neutron decay, neutron decay emissivity
! (3) xEnu(NP): neutrino energy array in MeV with dimension NP
! (4) T: Temperature in MeV
! (5) Mul, lepton_mass: relativistic chemical potentials & masses of
! e^-(mu^-) in MeV
! (6) Mun,Mup,neutron_mass,proton_mass,Un,Up: relativistic chemical
! potentials, masses, & potentials of neutron and proton in MeV
! relativsitic disperstion relation for nucleons,
! En=SQRT(p^2+neutron_mass^2)+Un
! Outputs:
! (7) OpaA(NP): opacities in km^-1; for neutron decay, emissivities
! in (s cm^3 MeV)^-1
!     **** blocking factor of final state neutrinos in neutron decay
!     is not considered ****
!=============================================================

CONTAINS

SUBROUTINE Opacity_CC_4D(opt_had, opt, xEnu, OpaA, xT, xMul, xMun, &
     xMup, xml, xmn, xmp, xUn, xUp)

  INTEGER,  INTENT(IN)  :: opt_had, opt
  REAL(DP), INTENT(IN)  :: xEnu, xT, xMul, xMun, xMup, xml, xmn, xmp, xUn, xUp
  REAL(DP), INTENT(OUT) :: OpaA

  INTEGER :: ii

  T = xT
  Mul = xMul
  Mun = xMun
  Mup = xMup
  lepton_mass = xml
  neutron_mass = xmn
  proton_mass = xmp
  Un = xUn
  Up = xUp
  jVA = 1.d0
  jAF = 1.d0

  ! Precompute some stuff you are going to use thousands of times
  mn2 = neutron_mass**2
  mp2 = proton_mass**2
  invT = 1/T
  
  SELECT CASE(opt_had)
  CASE(0)
    opt_form = 0
    opt_pseudo = 0
    opt_wm = 0
  CASE (1)
    opt_wm = 1
    opt_form = 0
    opt_pseudo = 0
  CASE (2)
    opt_wm = 1
    opt_pseudo = 1
    opt_form = 0
  CASE (3)
    opt_wm = 1
    opt_pseudo = 1
    opt_form = 1
  END SELECT

  OpaA = 0.0d0
  reaction_index = opt
  SELECT CASE(reaction_index)
  CASE (1)
    Enu = xEnu
    OpaA = Integration()
  CASE (2)
    CALL nubp_4D()
    Enu = xEnu
    OpaA = Integration()
  CASE (3,4)
    CALL nubp_4D()
    Mul = xMul
    Enu = xEnu
    OpaA = Integration_D()
  END SELECT

END SUBROUTINE Opacity_CC_4D

REAL(dp) FUNCTION Integration()

  ! Local
  REAL(DP) :: tmp, tmp1, mpt
  REAL(DP) :: integral(ncomp_C), error(ncomp_C), prob(ncomp_C) 
  INTEGER  :: neval, fail

  Integration = 0.

  pnmin = 0.
  pnmax = SQRT( (Tfac*T + Mun - Un)**2 - mn2 )

  mpt = proton_mass + lepton_mass
  IF(mpt < neutron_mass) THEN
     tmp = (1.d0 + mpt/neutron_mass)/(1.d0 - mpt**2/mn2)*Enu
  ELSE IF(mpt >= neutron_mass) THEN
     tmp = pnmax
  END IF
  tmp1 = SQRT( (tmp-Enu)**2+mpt**2 )+Up-SQRT(tmp**2+mn2)-Un
  IF(Enu<=tmp1) RETURN

  CALL Range_pn_4D()
  IF( pnmin > pnmax ) RETURN

!==========    vegas  ==============================================
! by default vegas is always used
  CALL Vegas(ndim_C, ncomp_C, integrand, userdata_C, nvec_C, &
    epsrel_C, epsabs_C, verbose_C, seed_C,  &
    mineval_C, maxeval_C, nstart_C, nincrease_C, nbatch_C,    &
    gridno_C, statefile_C, spin_C, &
    neval, fail, integral, error, prob)

  Integration = integral(1)/197.327d-18*2.d0
  RETURN

END FUNCTION Integration

REAL(dp) FUNCTION Integration_D()

  ! Local
  REAL(DP) :: integral(ncomp_C), error(ncomp_C), prob(ncomp_C) 
  INTEGER  :: neval, fail

  Integration_D = 0.0d0

  CALL Range_pn_dec_4D()
  IF( pnmin > pnmax ) RETURN
!==========      ==============================================
  CALL Vegas(ndim_C, ncomp_C, integrand_D, userdata_C, nvec_C, &
    epsrel_C, epsabs_C, verbose_C, seed_C, &
    mineval_C, maxeval_C, nstart_C, nincrease_C, nbatch_C,&
    gridno_C, statefile_C, spin_C, &
    neval, fail, integral, error, prob)

  IF(reaction_index .eq. 3) THEN
     Integration_D = integral(1) / 197.327d-18*2.d0
  ELSE IF(reaction_index .eq. 4) THEN
     Integration_D = integral(1)*Enu**2 / (2.d0*pi**2)*2.d0 &
        *(1.d13/197.327d0)**3*(3.d23/197.327d0)
  END IF

  RETURN

END FUNCTION Integration_D

INTEGER FUNCTION integrand(ndim, xx, ncomp, ff)

  INTEGER  :: ndim, ncomp
  REAL(DP) :: xx(*),ff(*)
  REAL(DP) :: xf2,xf3,xf4

  CALL initmomentum(xx)
  ff(1) = amp_4D()*q3*pn3/(16.d0*(2.d0*pi)**4*Enu**2*pp(0))*Jacob
  xf2 = 1.0d0/( EXP((pn(0) + Un - Mun)*invT) + 1.0d0 )
  xf3 = 1.0d0/( EXP((pe(0) - Mul)*invT) + 1.0d0 )
  xf4 = 1.0d0/( EXP((pp(0) + Up - Mup)*invT) + 1.0d0 )
  ff(1) = ff(1)*xf2*(1.d0 - xf4)*(1.d0 - xf3)
  integrand = 0

END FUNCTION integrand

SUBROUTINE initmomentum(x)
  
  REAL(DP):: x(4),En_min,En_max,tmp,tmp1,q3_max,q3_min,&
   cthe_max,cthe_min,cthe_n,Ep,pp3,q0,Ee,Pe3,tmp2,q3_lim(2),&
   q0_max,q0_min,tmp3,tmp4,cthe_v,phi_n

  En_min = MAX( SQRT(pnmin**2 + mn2) + Un, &
          proton_mass + Up - Enu + lepton_mass)
  En_max = SQRT(pnmax**2 + mn2) + Un
  En = (En_max - En_min)*x(1) + En_min

  pn3 = SQRT( (En - Un)**2 - mn2 )
  tmp = SQRT( (En - Up + Enu - lepton_mass)**2 - mp2 )

  tmp1 = pn3 + tmp
  tmp2 = Enu + SQRT( (Enu - proton_mass - Up + En)**2 - lepton_mass**2 )
  q3_max = MIN(tmp1, tmp2)

  CALL q3_sol(pn3,q3_lim)
  
  q3_min = q3_lim(1)
  q3_max = q3_lim(2)
  q3 = (q3_max-q3_min)*x(2) + q3_min

  tmp1 = Enu - SQRT( (Enu - q3)**2 + lepton_mass**2 )! Lp_max
  tmp2 = Enu - SQRT( (Enu + q3)**2 + lepton_mass**2 )! Lp_min
  tmp3 = SQRT( (pn3 + q3)**2 + mp2) + Up - En ! hd_max
  tmp4 = SQRT( (pn3 - q3)**2 + mp2) + Up - En ! hd_min
  q0_max = MIN(tmp1, tmp3)
  q0_min = MAX(tmp2, tmp4)

  tmp = mup - 50.0d0*T - En
  IF( tmp < q0_max ) THEN
     q0_min = MAX(q0_min, tmp)
  END IF
  tmp = Enu - MAX(lepton_mass, Mul - 50.0d0*T)
  IF( q0_min<tmp ) THEN
     q0_max = MIN(q0_max, tmp)
  END IF

  cthe_max = MIN(1.d0, ((q0_max + En - Up)**2 - mp2 - pn3**2 - q3**2) &
    /(2.d0*pn3*q3) )
  cthe_min = MAX(-1.d0, ((q0_min + En - Up)**2 - mp2 - pn3**2 - q3**2) &
    /(2.d0*pn3*q3) )
  cthe_n = (cthe_max - cthe_min)*x(3) + cthe_min
  phi_n = 2.d0*pi*x(4)

  pp3 = SQRT(pn3**2 + q3**2 + 2.d0*cthe_n*pn3*q3)
  Ep = SQRT(pp3**2 + mp2) + Up
  q0 = Ep - En
  Ee = Enu - q0
  pe3 = SQRT(Ee**2 - lepton_mass**2)

  cthe_v = (Enu**2 + q3**2 - pe3**2)/(2.d0*Enu*q3)
  cthe = (Enu**2 + pe3**2 - q3**2)/(2.d0*Enu*pe3)

  q(0:3) = (/q0-Up+Un,0.d0,0.d0,q3/)
  pv(0:3) = (/Enu,Enu*SQRT(1.d0-cthe_v**2),0.d0,Enu*cthe_v/)
  pe(0:3) = pv(0:3) - q(0:3)
  pe(0) = Enu - q0

  pn(0:3) = (/En-Un,pn3*SQRT(1.d0 - cthe_n**2)*cos(phi_n), &
    pn3*SQRT(1.d0-cthe_n**2)*sin(phi_n),pn3*cthe_n/)
  pp(0:3) = pn(0:3) + q(0:3)

  Jacob = (cthe_max - cthe_min)*(En_max - En_min) &
    *(q3_max - q3_min)*(2.d0*pi)

END SUBROUTINE initmomentum

SUBROUTINE nubp_4D()

  REAL(DP):: Mtmp, Utmp, mutmp
  
!  1) exchange neutron_mass <-> proton_mass; and Un <-> Up, mun <-> mup
  Mtmp = neutron_mass
  neutron_mass = proton_mass
  proton_mass = Mtmp
  mn2 = neutron_mass**2
  mp2 = proton_mass**2
  Utmp = Un
  Un = Up
  Up = Utmp
  mutmp = mun
  mun = mup
  mup = mutmp
!!  2) exchange Mul -> -Mul
  Mul = -Mul
!!  times -1 for VA & VF terms
  jVA = -1.d0
  jAF = -1.d0

END SUBROUTINE nubp_4D

SUBROUTINE Range_pn_4D()
  
  REAL(DP):: tmp,mpt,Epr,Esq,Equ,p2a,p2b,p20,F0,Finf,Fmax,Fmin,pnmax0

  pnmin = 0.
  pnmax = SQRT( (Tfac*T + mun - Un)**2 - mn2 )
  pnmax0 = pnmax

  mpt = proton_mass + lepton_mass
  F0 = SQRT(Enu**2 + mpt**2) + Up - neutron_mass - Un
  Finf = -Enu + Up - Un
  IF(mpt < neutron_mass) THEN
     tmp = (1.d0 + mpt/neutron_mass)/(1.d0 - mpt**2/mn2)*Enu
     Fmin = SQRT( (tmp - Enu)**2 + mpt**2 ) + Up - SQRT(tmp**2 + mn2) - Un
     Fmax = MAX( F0, Finf )
  ELSE
     Fmin = Finf
     Fmax = F0
  END IF

  IF( Fmax .le. Enu) RETURN
  IF( Fmin .ge. Enu) THEN
     pnmax = -1000.0d0
     RETURN
  END IF

  Epr = Enu - Up + Un
  Esq = Enu**2 + mpt**2 - mn2 - Epr**2
  Equ = Esq**2 + 4.d0*mn2*( Enu**2 - Epr**2 )
  IF(Equ .gt. 0.d0 .and. Enu**2 .ne. Epr**2) THEN
     p2a = ( Enu*Esq - ABS(Epr)*SQRT(Equ) )*0.5d0 / (Enu**2 - Epr**2)
     p2b = ( Enu*Esq + ABS(Epr)*SQRT(Equ) )*0.5d0 / (Enu**2 - Epr**2)
     IF(Epr < 0.d0 .and. p2b > 0.d0) THEN
      pnmin=p2b
     ELSE IF(Epr .ge. 0.d0 .and. p2a > 0.d0) THEN
        pnmin=p2a
     END IF

  ELSE IF(Enu**2 .eq. Epr**2) THEN
     p20 = -(4.0d0*Epr**2*mn2 - Esq**2)/(4.0d0*Esq*Enu)
     IF(p20 > 0.d0) pnmin = p20
  END IF

END SUBROUTINE Range_pn_4D

SUBROUTINE q3_sol(xPn3, q3_lim)
  
  REAL(DP), INTENT(IN)  :: xPn3
  REAL(DP), INTENT(OUT) :: q3_lim(2)
  ! Local
  REAL(DP) :: tmp1, tmp2, xA, xB, xC

  q3_lim = 0.
  tmp1 = Up - En - Enu
  tmp2 = Enu**2 + lepton_mass**2 - mp2 - xpn3**2
  xA = 4.d0*( tmp1**2 - (Enu - xpn3)**2 )
  xB = 4.d0*( (tmp2 - tmp1**2)*Enu - (tmp2 + tmp1**2)*xpn3 )
  xC = 4.d0*tmp1**2*(xpn3**2 + mp2) - (tmp2 - tmp1**2)**2

  q3_lim(1) = ABS( (xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA) )
  q3_lim(2) = (-xB + SQRT(xB**2-4.d0*xA*xC))/(2.d0*xA)

END SUBROUTINE q3_sol

INTEGER FUNCTION integrand_D(ndim, xx, ncomp, ff)

  INTEGER:: ndim,ncomp
  REAL(DP):: xx(*),ff(*)
  REAL(DP):: xf2,xf3,xf4

  CALL initmomentum_D(xx)
  ff(1) = amp_4D()*q3*pn3/(16.d0*(2.d0*pi)**4*Enu**2*pp(0))*Jacob
  xf2 = 1.0d0/( EXP((pn(0)+Un-mun)*invT) + 1.0d0 )
  xf3 = 1.0d0/( EXP((pe(0)-Mul)*invT) + 1.0d0 )
  xf4 = 1.0d0/( EXP((pp(0)+Up-mup)*invT) + 1.0d0 )

  IF(reaction_index.eq.3) THEN
     ff(1) = ff(1)*xf2*(1.d0-xf4)*xf3
  ELSE IF(reaction_index.eq.4) THEN
     ff(1) = ff(1)*(1.0d0-xf2)*xf4*(1.0d0-xf3)
  END IF
  integrand_D = 0

END FUNCTION integrand_D

SUBROUTINE initmomentum_D(x)

  REAL(DP):: x(4), En_min, En_max, tmp1, q3_max, q3_min, &
   cthe_max, cthe_min, cthe_n, Ep, pp3, q0, Ee, Pe3, tmp2, q3_lim(2), &
   q0_max, q0_min, tmp3, tmp4, cthe_v, phi_n

  En_min = SQRT(pnmin**2 + mn2) + Un
  En_max = SQRT(pnmax**2 + mn2) + Un
  En = (En_max - En_min)*x(1) + En_min
  pn3 = SQRT( (En - Un)**2 - mn2 )

  CALL q3_sol_D(pn3, q3_lim)
  q3_min = q3_lim(1)
  q3_max = q3_lim(2)
  q3 = (q3_max-q3_min)*x(2) + q3_min

  tmp1 = Enu + SQRT( (Enu + q3)**2 + lepton_mass**2 )
  tmp2 = Enu + SQRT( (Enu - q3)**2 + lepton_mass**2 )
  tmp3 = SQRT( (pn3 + q3)**2 + mp2) + Up - En
  tmp4 = SQRT( (pn3 - q3)**2 + mp2) + Up - En
  q0_max = MIN(tmp1, tmp3)
  q0_min = MAX(tmp2, tmp4)

  cthe_max = MIN(1.d0, ((q0_max + En - Up)**2 &
      - mp2 - pn3**2 - q3**2) / (2.d0*pn3*q3))
  cthe_min = MAX(-1.d0, ((q0_min + En - Up)**2 &
      - mp2 - pn3**2 - q3**2) / (2.d0*pn3*q3))

  cthe_n = (cthe_max - cthe_min)*x(3) + cthe_min
  phi_n = 2.d0*pi*x(4)

  pp3 = SQRT(pn3**2 + q3**2 + 2.d0*cthe_n*pn3*q3)
  Ep = SQRT(pp3**2 + mp2) + Up
  q0 = Ep - En
  Ee = q0 - Enu
  pe3 = SQRT(Ee**2 - lepton_mass**2)

  cthe_v = (Enu**2 + q3**2 - pe3**2) / (2.d0*Enu*q3)
  cthe = -(Enu**2 + pe3**2 - q3**2) / (2.d0*Enu*pe3)

  q(0:3) = (/q0 - Up + Un, 0.d0, 0.d0, q3/)
  pv(0:3) = (/Enu, Enu*SQRT(1.d0 - cthe_v**2), 0.d0, Enu*cthe_v/)
  pe(0:3) = q(0:3) - pv(0:3)
  pe(0) = q0 - Enu

  pn(0:3) = (/En - Un, pn3*SQRT(1.d0 - cthe_n**2)*cos(phi_n),&
      pn3*SQRT(1.d0 - cthe_n**2)*sin(phi_n), pn3*cthe_n/)
  pp(0:3) = pn(0:3) + q(0:3)

  Jacob = (cthe_max - cthe_min)*(En_max - En_min)&
      *(q3_max - q3_min)*(2.d0*pi)

END SUBROUTINE initmomentum_D

SUBROUTINE Range_pn_dec_4D()
  
  REAL(DP):: tmp, mpt, Epr, Esq, Equ, p2a, p2b, p20
  REAL(DP):: F0, Fmax

  pnmin = 0.
  ppmax = SQRT( (Tfac*T + mup - Up)**2 - mp2 )
  SELECT CASE(reaction_index)
  CASE (3)
  pnmax = SQRT( (Tfac*T + mun - Un)**2 - mn2 )
  CASE (4)
  pnmin = SQRT( (MAX(neutron_mass, mun - Tfac*T - Un))**2 - mn2 )
  pnmax = SQRT( (SQRT(ppmax**2 + mp2) + Up - Un)**2 - mn2 )
  END SELECT

  mpt = proton_mass - lepton_mass
  IF( mpt .gt. neutron_mass) THEN
     tmp = neutron_mass/(mpt - neutron_mass)*Enu
  ELSE
     tmp = pnmax
  END IF
  Fmax = SQRT( (tmp + Enu)**2 + mpt**2 ) + Up - SQRT(tmp**2 + mn2) - Un

  IF(Enu>=Fmax) THEN
     pnmax = -1000.
     RETURN
  END IF

  F0 = SQRT(Enu**2 + mpt**2) + Up - neutron_mass - Un
  IF( F0 .ge. Enu .and. Up .ge. Un) RETURN

  Epr = Enu - Up + Un
  Esq = Enu**2 + mpt**2 - mn2 - Epr**2
  Equ = Esq**2 + 4.d0*( (neutron_mass*Enu)**2 - (Epr*neutron_mass)**2 )
  IF(Equ .gt. 0.d0 .and. Enu**2 .ne. Epr**2) THEN
     p2a = -( Enu*Esq - ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2 - Epr**2)
     p2b = -( Enu*Esq + ABS(Epr)*SQRT(Equ) )*0.5d0/(Enu**2 - Epr**2)
     IF(Epr < 0.d0 .and. p2b > 0.d0) THEN
      pnmin = p2b
     ELSE IF(Epr .ge. 0.d0 .and. p2a > 0.d0) THEN
      pnmin = p2a
     END IF

  ELSE IF(Enu**2 .eq. Epr**2) THEN
     p20 = (4.0d0*Epr**2*mn2 - Esq**2)/(4.0d0*Esq*Enu)
     IF(p20 > 0.d0) pnmin = p20
  END IF

END SUBROUTINE Range_pn_dec_4D

SUBROUTINE q3_sol_D(xPn3, q3_lim)

  REAL(DP), INTENT(IN)  :: xPn3
  REAL(DP), INTENT(OUT) :: q3_lim(1:2)
  REAL(DP) :: tmp1, tmp2, xA, xB, xC
  REAL(DP) :: tmp3, tmp4, tmp5, tmp6

  q3_lim = 0.0d0

  tmp1 = Up - En - Enu
  tmp2 = Enu**2 + lepton_mass**2 - mp2 - xpn3**2
  xA = 4.d0*( tmp1**2 - (Enu + xpn3)**2 )
  xB = 4.d0*( (tmp2 - tmp1**2)*Enu + (tmp2 + tmp1**2)*xpn3 )
  xC = 4.d0*tmp1**2*(xpn3**2 + mp2) - (tmp2 - tmp1**2)**2
  q3_lim(1) = ABS( (xB + SQRT(xB**2 - 4.d0*xA*xC))/(2.d0*xA) )
  q3_lim(2) = ABS( (-xB + SQRT(xB**2 - 4.d0*xA*xC))/(2.d0*xA) )

  tmp3 = Enu + SQRT( (Enu-q3_lim(1))**2 + lepton_mass**2 ) &
    -(SQRT( (xpn3+q3_lim(1))**2+mp2) + Up - En)
  tmp4 = Enu + SQRT( (Enu-q3_lim(2))**2 + lepton_mass**2 ) &
    -(SQRT( (xpn3+q3_lim(2))**2+mp2) + Up - En)
  tmp5 = Enu + SQRT( (Enu+q3_lim(1))**2 + lepton_mass**2 ) &
    -(SQRT( (xpn3-q3_lim(1))**2+mp2) + Up - En)
  tmp6 = Enu + SQRT( (Enu+q3_lim(2))**2 + lepton_mass**2 ) &
    -(SQRT( (xpn3-q3_lim(2))**2+mp2) + Up - En)

  IF(ABS(tmp3) > 1.d-10 .and. ABS(tmp5) > 1.d-10) THEN
   q3_lim(1)=0.0d0
  END IF

  IF(ABS(tmp4) > 1.d-10 .and. ABS(tmp6) > 1.d-10) THEN
   q3_lim(2)=MAX( q3_lim(1) + 20.0d0*T, Enu + (20.0d0*T + Mul) )
  END IF

END SUBROUTINE q3_sol_D

FUNCTION amp_4D()
  
  REAL(DP):: Mvv,Maa,Mff,Mva,Mvf,Maf,Mpp,Map,tmp1,tmp2,tmp3,&
     d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d55,amp_4D,gv_eff,ga_eff,f2,gp

  d1 = cdot4(pn,pv)
  d2 = cdot4(pp,pe)
  d3 = cdot4(pn,pe)
  d4 = cdot4(pp,pv)
  d5 = cdot4(pv,pe)
  d55 = cdot4(pn,pp)
  tmp1 = d1*d2
  tmp2 = d3*d4
  tmp3 = neutron_mass*proton_mass*d5
  d6 = cdot4(pv,q)
  d7 = cdot4(pe,q)
  d8 = cdot4(pn,q)
  d9 = cdot4(pp,q)
  d10 = cdot4(q,q)

  gv_eff = gv
  f2 = F2wm0
  ga_eff = ga
  gp = 0.0d0

  IF(opt_form .eq. 1) THEN
     gv_eff=( 1.d0-4.706d0*d10/(4.d0*Mnp**2) )/( 1.d0 - d10/(4.d0*Mnp**2)) &
      /( 1.d0 - d10/massV**2 )**2
     f2 = F2wm0/( 1.d0 - d10/(4.d0*Mnp**2) )/( 1.d0 - d10/massV**2 )**2
     ga_eff = ga/( 1.d0 - d10/massA**2 )**2
  END IF
  IF(opt_pseudo.eq.1) gp=2.d0*Mnp**2*ga_eff/(mpi**2 - d10)  

  Mvv = 16.d0*(Gw_MeV*Vud*gv_eff)**2*( tmp1 + tmp2 - tmp3 )
  Maa = 16.d0*(Gw_MeV*Vud*ga_eff)**2*( tmp1 + tmp2 + tmp3 )
  Mva = jVA*32.d0*GfVud2*gv_eff*ga_eff*( tmp1 - tmp2 )
  Mvf = 8.d0*GfVud2*gv_eff*f2/Mnp*(   &
    (d1*proton_mass - d4*neutron_mass)*d7 + (d3*proton_mass-d2*neutron_mass)*d6 + (d8*proton_mass - d9*neutron_mass)*d5 )
  Maf = jAF*16.d0*GfVud2*ga_eff*f2/Mnp*((d1*proton_mass + d4*neutron_mass)*d7 &
       -(d3*proton_mass+d2*neutron_mass)*d6)
  Mff = 2.d0*(Gw_MeV*Vud*f2/Mnp)**2*(      &
    d10*(d55*d5-2.d0*d1*d2-2.d0*d3*d4) &
    + 2.d0*d7*(d1*d9 + d4*d8 - d55*d6) + 2.d0*d6*(d3*d9 + d2*d8) &
    - neutron_mass*proton_mass*(d5*d10 + 2.d0*d6*d7) )

  Mpp = 8.d0*(Gw_MeV*Gp*Vud/Mnp)**2*(2.d0*d6*d7-d10*d5)*(d55-neutron_mass*proton_mass)
  Map = 16.d0*GfVud2*ga_eff*gP/mnp*(  &
    d6*(d3*proton_mass-d2*neutron_mass) + d7*(d1*proton_mass-d4*neutron_mass) + d5*(d9*neutron_mass-d8*proton_mass) )

  amp_4D = Mvv + Maa + Mva + opt_wm*(Mvf + Maf + Mff) &
          + opt_pseudo*(Mpp + Map) 

END FUNCTION amp_4D

REAL(DP) PURE FUNCTION cdot4(xp1,xp2)

  REAL(DP), INTENT(IN) :: xp1(0:3), xp2(0:3)

  cdot4 = xp1(0)*xp2(0) - xp1(1)*xp2(1) - xp1(2)*xp2(2) - xp1(3)*xp2(3)
  
END FUNCTION cdot4

END MODULE wlSemiLeptonicOpacityModule4D

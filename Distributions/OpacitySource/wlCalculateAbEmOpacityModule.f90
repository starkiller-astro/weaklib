MODULE wlCalculateAbEmOpacityModule

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    Gw, hbar, cvel, pi, ga, gv, Gw_MeV, &
    Vud, rmu, mp, me, mn, mmu, kmev
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlInterpolationUtilitiesModule, ONLY: &
    LinearInterp_Array_Point, &
    GetIndexAndDelta_Lin, GetIndexAndDelta_Log

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ElasticAbsorptionOpacityNum
  PUBLIC :: ElasticAbsorptionOpacityNue
  PUBLIC :: CalculateHorowitzWeakMagRecoil

  ! REAL(dp), PARAMETER  :: g2      = ( Gw/mp**2 )**2 * ( hbar * cvel )**2   ! square of the Fermi constant in units of cm^{2} MeV^{-2}
  ! REAL(dp), PARAMETER  :: C_Op_el = g2/pi*(gv**2 + 3*ga**2)
  REAL(dp), PARAMETER  :: C_Op_el = Gw_MeV**2 * Vud**2 * ( hbar * cvel )**2/pi*(gv**2 + 3*ga**2)

  ! Constants for Weak Magnetism Corrections
  REAL(dp), PARAMETER  :: Lambda_p = 1.793d0
  REAL(dp), PARAMETER  :: Lambda_n = -1.913d0
  REAL(dp), PARAMETER  :: sin2tw = 0.2325d0

CONTAINS

  SUBROUTINE ElasticAbsorptionOpacityNum(D, T, Enu, Mun, Mn_eff, Un, Xn, &
    Mup, Mp_eff, Up, Xp, Mumu, OpElNum, OpElNumBar)
    
    REAL(DP), INTENT(IN) :: D, T, Enu
    REAL(DP), INTENT(IN) :: Mun, Mn_eff, Un, Xn
    REAL(DP), INTENT(IN) :: Mup, Mp_eff, Up, Xp
    REAL(DP), INTENT(IN) :: Mumu
    REAL(DP), INTENT(OUT) :: OpElNum, OpElNumBar

    ! LOCAL VARIABLES
    REAL(DP) :: n_n, n_p

    n_n = Xn * D / rmu
    n_p = Xp * D / rmu

    ! num   + n -> mu- + p
    CALL ElasticAbsorptionOpacity(T, Enu, Mun, mn, Mn_eff, Un, n_n, &
      Mumu, mmu, Mup, mp, Mp_eff, Up, n_p, OpElNum)

    ! numbar + p -> mu+ + n
    CALL ElasticAbsorptionOpacity(T, Enu, Mup, mp, Mp_eff, Up, n_p, &
      -Mumu, mmu, Mun, mn, Mn_eff, Un, n_n, OpElNumBar)

  END SUBROUTINE ElasticAbsorptionOpacityNum


  SUBROUTINE ElasticAbsorptionOpacityNue(D, T, Enu, Mun, Mn_eff, Un, Xn, &
    Mup, Mp_eff, Up, Xp, Mue, OpElNue, OpElNueBar)
  
    REAL(DP), INTENT(IN)  :: D, T, Enu
    REAL(DP), INTENT(IN)  :: Mun, Mn_eff, Un, Xn
    REAL(DP), INTENT(IN)  :: Mup, Mp_eff, Up, Xp
    REAL(DP), INTENT(IN)  :: Mue
    REAL(DP), INTENT(OUT) :: OpElNue, OpElNueBar

    ! LOCAL VARIABLES
    REAL(DP) :: n_n, n_p

    n_n = Xn * D / rmu
    n_p = Xp * D / rmu

    ! nue   + n -> e- + p
    CALL ElasticAbsorptionOpacity(T, Enu, Mun, mn, Mn_eff, Un, n_n, &
      Mue, me, Mup, mp, Mp_eff, Up, n_p, OpElNue)

    ! nuebar + p -> e+ + n
    CALL ElasticAbsorptionOpacity(T, Enu, Mup, mp, Mp_eff, Up, n_p, &
      -Mue, me, Mun, mn, Mn_eff, Un, n_n, OpElNueBar)

  END SUBROUTINE ElasticAbsorptionOpacityNue
  
  SUBROUTINE ElasticAbsorptionOpacity(T, E1, Mu2, M2, M2_eff, U2, n2, &
    Mu3, M3, Mu4, M4, M4_eff, U4, n4, OpEl)

    ! For absorptivity of neutrinos we have: 
    ! 1 is neutrino, 2 is neutron, 3 is lepton, 4 is proton
    ! Notice that 
    REAL(DP), INTENT(IN)  :: T, E1
    REAL(DP), INTENT(IN)  :: Mu2, M2, M2_eff, U2, n2
    REAL(DP), INTENT(IN)  :: Mu3, M3
    REAL(DP), INTENT(IN)  :: Mu4, M4, M4_eff, U4, n4
    REAL(DP), INTENT(OUT) :: OpEl

    ! LOCAL
    REAL(DP) :: f3, baryons_term, phase_space_term, E3, phi2, phi4

    E3 = E1 + (M2-M4) + (U2-U4)

    IF (E3 < M3) THEN
        OpEl = 0.0d0
        RETURN
    ENDIF
      
    phi2 = Mu2 - M2_eff - U2
    phi4 = Mu4 - M4_eff - U4

    baryons_term = (n2 - n4)/( 1.0d0 - EXP((phi4 - phi2)/T/kmev))
    phase_space_term = E3**2 * SQRT(1.0d0 - (M3/E3)**2)
    f3 = 1.0d0/(EXP((E3 - Mu3)/T/kmev) + 1.0d0)

    OpEl = C_Op_el * baryons_term * phase_space_term * (1.0d0 - f3)

  END SUBROUTINE ElasticAbsorptionOpacity

  SUBROUTINE CalculateHorowitzWeakMagRecoil(E, WeakMagRecNuLep, WeakMagRecNuLepBar)
  ! This is Taken from Tobias' Fischer routines. In principle it can be
  ! Precalculated for each energy at the beginning...

    REAL(DP), INTENT(in) :: E
    REAL(DP), INTENT(out) :: WeakMagRecNuLep, WeakMagRecNuLepBar

    ! LOCAL
    REAL(DP) :: small_E, cv, ca, F2
    REAL(DP) :: temp_wm1, temp_wm2, temp_wm3
    INTEGER  :: iE

    ! Now for absorption of neutrinos on neutrons
    small_E = E/mn

    CALL CalculateNucleonFormFactors(small_E, 3, ca, cv, F2)

    temp_wm1 = cv*cv*(1.0d0 + 4.0d0*small_E + 16.0d0/3.0d0*small_E*small_E) &
      + 3.0d0*ca*ca*(1.0d0 + 4.0d0/3.0d0*small_E)**2 &
      + 8.0d0/3.0d0*cv*F2*small_E*small_E &
      + 5.0d0/3.0d0*small_E*small_E*(1.0d0 + 2.0d0/5.0d0*small_E)*F2*F2

    temp_wm2 = 4.0d0*(cv + F2)*ca*small_E*(1. + 4./3.0d0*small_E)
    temp_wm3 = (cv*cv + 3.00d0*ca*ca)*(1.0d0 + 2.0d0*small_E)**3

    WeakMagRecNuLep    = (temp_wm1 + temp_wm2)/temp_wm3

    ! Now for absorption of antineutrinos on protons
    small_E = E/mp

    CALL CalculateNucleonFormFactors(small_E, 4, ca, cv, F2)

    temp_wm1 = cv*cv*(1.0d0 + 4.0d0*small_E + 16.0d0/3.0d0*small_E*small_E) &
      + 3.0d0*ca*ca*(1.0d0 + 4.0d0/3.0d0*small_E)**2 &
      + 8.0d0/3.0d0*cv*F2*small_E*small_E &
      + 5.0d0/3.0d0*small_E*small_E*(1.0d0 + 2.0d0/5.0d0*small_E)*F2*F2

    temp_wm2 = 4.0d0*(cv + F2)*ca*small_E*(1. + 4./3.0d0*small_E)
    temp_wm3 = (cv*cv + 3.00d0*ca*ca)*(1.0d0 + 2.0d0*small_E)**3

    WeakMagRecNuLepBar = (temp_wm1 - temp_wm2)/temp_wm3

  END SUBROUTINE CalculateHorowitzWeakMagRecoil

  SUBROUTINE CalculateNucleonFormFactors(small_E, WhichReaction, ca, cv, F2)

    REAL(DP), INTENT(IN)  :: small_E
    INTEGER , INTENT(IN)  :: WhichReaction
    REAL(DP), INTENT(OUT) :: ca, cv, F2

    REAL(DP) :: tau, eta, G, Fp1, Fp2, Fn1, Fn2

    !-----------------------------------------------------------------------
    !
    !     Input:
    !     WhichReaction
    !     1 ... nu/nubar + p --> nu/nubar + p
    !     2 ... nu/nubar + n --> nu/nubar + n
    !     3 ... nu       + n --> l-       + p
    !     4 ... nubar    + p --> l+       + n
    !
    !     Output: cv, ca, and F2 as defined in Appendix B of
    !     Horowitz, Phys. Rev. D65, 043001 (2002)
    !
    !     We assume an average angle x=0, such that Q^2 = 2k^2/(1+k/M)
    !
    !     Luca : I think 1 and 2 are valid both for nu and nubar scattering, but 
    !     I might be wrong...
    !-----------------------------------------------------------------------

    tau = 0.5d0*small_E**2/(1.0d0 + small_E)
    eta = 1.0d0/(1.0d0 + 5.6d0*tau)
    G = 1.0d0/(1.0d0 + 4.97d0*tau)**2
    Fp1 = (1.0d0 + tau*(1.0d0 + Lambda_p))*G/(1.0d0 + tau)
    Fp2 = Lambda_p*G/(1.0d0 + tau)
    Fn1 = tau*Lambda_n*(1.0d0 - eta)*G/(1.0d0 + tau)
    Fn2 = Lambda_n*(1.0d0 + tau*eta)*G/(1.0d0 + tau)
      
    !.....form factors......................................................
    IF (WhichReaction .eq. 1) THEN
      cv = (0.5d0 - 2.0d0*sin2tw)*Fp1 - 0.5d0*Fn1
      ca = 0.5d0*ga/(1.0d0 + 3.53d0*tau)**2
      F2 = (0.5d0 - 2.0d0*sin2tw)*Fp2 - 0.5d0*Fn2
    ELSE IF (WhichReaction .eq. 2) THEN
      cv = (0.5d0 - 2.0d0*sin2tw)*Fn1 - 0.5d0*Fp1
      ca = -0.5d0*ga/(1. + 3.53d0*tau)**2
      F2 = (0.5d0 - 2.0d0*sin2tw)*Fn2 - 0.5d0*Fp2
    ELSE
      cv = Fp1 - Fn1
      ca = ga/(1.0d0 + 3.53d0*tau)**2
      F2 = Fp2 - Fn2
    ENDIF

  END SUBROUTINE CalculateNucleonFormFactors

END MODULE wlCalculateAbEmOpacityModule

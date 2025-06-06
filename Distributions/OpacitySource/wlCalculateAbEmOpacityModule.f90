MODULE wlCalculateAbEmOpacityModule

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    Gw, hbar, cvel, pi, ga, gv, &
    rmu, mp, me, mn, mmu, cm3fm3, kmev
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType, &
    EquationOfState4DTableType, &
    EquationOfStateCompOSETableType
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlInterpolationUtilitiesModule, ONLY: &
    LinearInterp_Array_Point, &
    GetIndexAndDelta_Lin, GetIndexAndDelta_Log
  USE wlMuonEOS, ONLY: &
    MuonStateType, FullMuonEOS
  USE wlLeptonEOSModule, ONLY: &
    MuonTableType

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ElasticAbsorptionOpacityNuMu

  REAL(dp), PARAMETER  :: g2      = ( Gw/mp**2 )**2 * ( hbar * cvel )**2   ! square of the Fermi constant in units of cm^{2} MeV^{-2}
  REAL(dp), PARAMETER  :: C_Op_el = g2/pi*(gv**2 + 3*ga**2)/cm3fm3

  ! Constants for Weak Magnetism Corrections
  REAL(dp), PARAMETER  :: Lambda_p = 1.793d0
  REAL(dp), PARAMETER  :: Lambda_n = -1.913d0
  REAL(dp), PARAMETER  :: sin2tw = 0.2325d0

  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: WeakMagRecNuLep, WeakMagRecNuLepBar
  LOGICAL :: WeakMagRecoilInitialized = .FALSE.

CONTAINS

  SUBROUTINE ElasticAbsorptionOpacityNuMu(D, T, Ye, Ym, Enu, nE, &
    EOSTable, MuonTable, OpElNuMu, OpElNuMuBar, IncludeWeakMagRecoil)
  
    REAL(DP), INTENT(IN)  :: D, T, Ye, Ym
    INTEGER , INTENT(IN)  :: nE
    REAL(DP), INTENT(IN), DIMENSION(:) :: Enu(nE)
    LOGICAL,  OPTIONAL  :: IncludeWeakMagRecoil

#ifdef EOSMODE_3D
    TYPE(EquationOfStateTableType), INTENT(in) :: EOSTable
#elif defined(EOSMODE_4D)
    TYPE(EquationOfState4DTableType), INTENT(in) :: EOSTable
#elif defined(EOSMODE_COMPOSE)
    TYPE(EquationOfStateCompOSETableType), INTENT(in) :: EOSTable
#endif
    TYPE(MuonTableType), INTENT(in) :: MuonTable

    REAL(DP), INTENT(OUT), DIMENSION(:) :: OpElNuMu(nE), OpElNuMuBar(nE)

    ! LOCAL VARIABLES
    REAL(DP) :: Yp, Xp, Xn, dYp
    REAL(DP) :: Mup, Mun, Mumu
    REAL(DP) :: Mp_eff, Mn_eff, Up, Un, n_n, n_p
    INTEGER  :: iEOS_Rho, iEOS_T, iEOS_Yp, iDV, iE
    LOGICAL  :: DoWeakMagRec

    TYPE(MuonStateType) :: MuonState
    
    IF ( PRESENT(IncludeWeakMagRecoil) ) THEN
      DoWeakMagRec = IncludeWeakMagRecoil
    ELSE
      DoWeakMagRec = .TRUE.
    END IF

    IF ((DoWeakMagRec) .AND. .NOT. WeakMagRecoilInitialized) THEN
      CALL CalculateHorowitzWeakMagRecoil(Enu, nE)
      WeakMagRecoilInitialized = .TRUE.
    ENDIF

    iEOS_Rho = EOSTable % TS % Indices % iRho
    iEOS_T   = EOSTable % TS % Indices % iT
    iEOS_Yp  = EOSTable % TS % Indices % iYe ! Ye is actually Yp, but for consistency reasons it was not changed

    Yp = Ye + Ym

    iDV = EOSTable % DV % Indices % iNeutronChemicalPotential
    CALL LogInterpolateSingleVariable   &
        ( D, T, Yp, &
          EOSTable % TS % States(iEOS_Rho) % Values, &
          EOSTable % TS % States(iEOS_T) % Values,   &
          EOSTable % TS % States(iEOS_Yp) % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Mun )

    iDV = EOSTable % DV % Indices % iProtonChemicalPotential
    CALL LogInterpolateSingleVariable   &
        ( D, T, Yp, &
          EOSTable % TS % States(iEOS_Rho) % Values, &
          EOSTable % TS % States(iEOS_T) % Values,   &
          EOSTable % TS % States(iEOS_Yp) % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Mup )

#ifdef EOSMODE_3D
    Mn_eff = 0.0d0
    Mp_eff = 0.0d0
    Up     = 0.0d0
    Un     = 0.0d0
#else
    iDV = EOSTable % DV % Indices % iNeutronEffMass
    CALL LogInterpolateSingleVariable   &
        ( D, T, Yp, &
          EOSTable % TS % States(iEOS_Rho) % Values, &
          EOSTable % TS % States(iEOS_T) % Values,   &
          EOSTable % TS % States(iEOS_Yp) % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Mn_eff )

    iDV = EOSTable % DV % Indices % iNeutronSelfEnergy
    CALL LogInterpolateSingleVariable   &
        ( D, T, Yp, &
          EOSTable % TS % States(iEOS_Rho) % Values, &
          EOSTable % TS % States(iEOS_T) % Values,   &
          EOSTable % TS % States(iEOS_Yp) % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Un )

    iDV = EOSTable % DV % Indices % iProtonEffMass
    CALL LogInterpolateSingleVariable   &
        ( D, T, Yp, &
          EOSTable % TS % States(iEOS_Rho) % Values, &
          EOSTable % TS % States(iEOS_T) % Values,   &
          EOSTable % TS % States(iEOS_Yp) % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Mp_eff )

    iDV = EOSTable % DV % Indices % iProtonSelfEnergy
    CALL LogInterpolateSingleVariable   &
        ( D, T, Yp, &
          EOSTable % TS % States(iEOS_Rho) % Values, &
          EOSTable % TS % States(iEOS_T) % Values,   &
          EOSTable % TS % States(iEOS_Yp) % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Up )
#endif

    iDV = EOSTable % DV % Indices % iProtonMassFraction
    CALL LogInterpolateSingleVariable   &
        ( D, T, Yp, &
          EOSTable % TS % States(iEOS_Rho) % Values, &
          EOSTable % TS % States(iEOS_T) % Values,   &
          EOSTable % TS % States(iEOS_Yp) % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Xp )

    iDV = EOSTable % DV % Indices % iNeutronMassFraction
    CALL LogInterpolateSingleVariable   &
        ( D, T, Yp, &
          EOSTable % TS % States(iEOS_Rho) % Values, &
          EOSTable % TS % States(iEOS_T)   % Values,   &
          EOSTable % TS % States(iEOS_Yp)  % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Xn )

    n_n = Xn * D * cm3fm3 / rmu
    n_p = Xp * D * cm3fm3 / rmu

    MuonState % t     = T
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    Mumu = MuonState % mu

    WRITE(*,*) 't     =', T * kmev
    WRITE(*,*) 'ye    =', ye
    WRITE(*,*) 'ym    =', ym
    WRITE(*,*) 'd     =', d
    WRITE(*,*) 'yp    =', yp
    WRITE(*,*) 'chemmu=', Mumu
    WRITE(*,*) 'chemn =', Mun - mn
    WRITE(*,*) 'chemp =', Mup - mp
    WRITE(*,*) 'xn    =', xn
    WRITE(*,*) 'xp    =', xp
    WRITE(*,*) 'un    =', un
    WRITE(*,*) 'up    =', up
    WRITE(*,*) 'massn =', Mn_eff
    WRITE(*,*) 'massp =', Mp_eff
    WRITE(*,*) 'nB    =', D * cm3fm3 / rmu
    WRITE(*,*) 'phin  =', Mun - Mn_eff - Un
    WRITE(*,*) 'phip  =', Mup - Mp_eff - Up

#if defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iE, OpElNuMu, OpElNuMuBar ) & 
    !$OMP SHARED( T, Enu, Mup, Mp_eff, Up, n_p, &
    !$OMP         Mumu, Mun, Mn_eff, Un, n_n )
#endif
    DO iE=1, nE

      ! numu   + n -> mu- + p
      CALL ElasticAbsorptionOpacity(T, Enu(iE), Mun, mn, Mn_eff, Un, n_n, &
        Mumu, mmu, Mup, mp, Mp_eff, Up, n_p, OpElNuMu(iE))

      ! numubar + p -> mu+ + n
      CALL ElasticAbsorptionOpacity(T, Enu(iE), Mup, mp, Mp_eff, Up, n_p, &
        -Mumu, mmu, Mun, mn, Mn_eff, Un, n_n, OpElNuMuBar(iE))

      OpElNuMu(iE)    = OpElNuMu(iE)   
      OpElNuMuBar(iE) = OpElNuMuBar(iE)

    ENDDO
#if defined(WEAKLIB_OMP)
    !$OMP END PARALLEL DO
#endif

  IF (DoWeakMagRec) THEN
    OpElNuMu(iE)    = OpElNuMu(iE)    * WeakMagRecNuLep   (iE)
    OpElNuMuBar(iE) = OpElNuMuBar(iE) * WeakMagRecNuLepBar(iE)
  ENDIF

  END SUBROUTINE ElasticAbsorptionOpacityNuMu

  SUBROUTINE ElasticAbsorptionOpacity(T, E1, Mu2, M2, M2_eff, U2, n2, &
    Mu3, M3, Mu4, M4, M4_eff, U4, n4, OpEl)

  ! For absorptivity of muon neutrinos we have: 
  ! 1 is neutrino, 2 is neutron, 3 is muon, 4 is proton
    REAL(DP), INTENT(IN)  :: T, E1
    REAL(DP), INTENT(IN)  :: Mu2, M2, M2_eff, U2, n2
    REAL(DP), INTENT(IN)  :: Mu3, M3
    REAL(DP), INTENT(IN)  :: Mu4, M4, M4_eff, U4, n4
    REAL(DP), INTENT(OUT) :: OpEl

    ! LOCAL
    REAL(DP) :: fmu, baryons_term, phase_space_term, E3, phi2, phi4

    E3 = E1 + (M2-M4) + (U2-U4)

    IF (E3 < M3) THEN
        OpEl = 0.0d0
        RETURN
    ENDIF
      
    phi2 = Mu2 - M2_eff - U2
    phi4 = Mu4 - M4_eff - U4

    baryons_term = (n2 - n4)/( 1.0d0 - EXP((phi4 - phi2)/T/kmev))
    phase_space_term = E3**2 * SQRT(1.0d0 - (M3/E3)**2)
    fmu = 1.0d0/(EXP((E3 - Mu3)/T) + 1.0d0)

    OpEl = C_Op_el * baryons_term * phase_space_term * (1.0d0 - fmu)

  END SUBROUTINE ElasticAbsorptionOpacity

  SUBROUTINE CalculateHorowitzWeakMagRecoil(E, nE)
  ! This is Taken from Tobias' Fischer routines. In principle it can be
  ! Precalculated for each energy at the beginning...

    INTEGER , INTENT(in)               :: nE
    REAL(DP), INTENT(in), DIMENSION(:) :: E(nE)

    ! LOCAL
    REAL(DP) :: small_E, cv, ca, F2
    REAL(DP) :: temp_wm1, temp_wm2, temp_wm3
    INTEGER  :: iE

    ALLOCATE( WeakMagRecNuLep(nE) )
    ALLOCATE( WeakMagRecNuLepBar(nE) )

#if defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( iE, small_E, cv, ca, F2, &
    !$OMP          temp_wm1, temp_wm2, temp_wm3 )
#endif
    DO iE=1, nE

      ! Now for absorption of neutrinos on neutrons
      small_E = E(iE)/mn

      CALL CalculateNucleonFormFactors(small_E, 3, ca, cv, F2)

      temp_wm1 = cv*cv*(1.0d0 + 4.0d0*small_E + 16.0d0/3.0d0*small_E*small_E) &
        + 3.0d0*ca*ca*(1.0d0 + 4.0d0/3.0d0*small_E)**2 &
        + 8.0d0/3.0d0*cv*F2*small_E*small_E &
        + 5.0d0/3.0d0*small_E*small_E*(1.0d0 + 2.0d0/5.0d0*small_E)*F2*F2

      temp_wm2 = 4.0d0*(cv + F2)*ca*small_E*(1. + 4./3.0d0*small_E)
      temp_wm3 = (cv*cv + 3.00d0*ca*ca)*(1.0d0 + 2.0d0*small_E)**3

      WeakMagRecNuLep(iE)    = (temp_wm1 + temp_wm2)/temp_wm3

      ! Now for absorption of antineutrinos on protons
      small_E = E(iE)/mp

      CALL CalculateNucleonFormFactors(small_E, 4, ca, cv, F2)

      temp_wm1 = cv*cv*(1.0d0 + 4.0d0*small_E + 16.0d0/3.0d0*small_E*small_E) &
        + 3.0d0*ca*ca*(1.0d0 + 4.0d0/3.0d0*small_E)**2 &
        + 8.0d0/3.0d0*cv*F2*small_E*small_E &
        + 5.0d0/3.0d0*small_E*small_E*(1.0d0 + 2.0d0/5.0d0*small_E)*F2*F2

      temp_wm2 = 4.0d0*(cv + F2)*ca*small_E*(1. + 4./3.0d0*small_E)
      temp_wm3 = (cv*cv + 3.00d0*ca*ca)*(1.0d0 + 2.0d0*small_E)**3

      WeakMagRecNuLepBar(iE) = (temp_wm1 - temp_wm2)/temp_wm3
    ENDDO
#if defined(WEAKLIB_OMP)
    !$OMP END PARALLEL DO
#endif

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

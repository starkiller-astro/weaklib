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
  USE wlMuonEOS, ONLY: &
    MuonStateType, FullMuonEOS
  USE wlLeptonEOSModule, ONLY: &
    MuonTableType

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ElasticAbsorptionOpacityNuMu
  PUBLIC :: ElasticAbsorptionOpacity

  REAL(dp), PARAMETER  :: g2      = ( Gw/mp**2 )**2 * ( hbar * cvel )**2   ! square of the Fermi constant in units of cm^{2} MeV^{-2}
  REAL(dp), PARAMETER  :: C_Op_el = g2/pi*(gv**2 + 3*ga**2)/cm3fm3

CONTAINS

  SUBROUTINE ElasticAbsorptionOpacityNuMu(D, T, Ye, Ym, Enu, &
    EOSTable, MuonTable, OpElNuMu, OpElNuMuBar)
  
  ! 1 is neutrino, 2 is neutron, 3 is muon, 4 is proton
    REAL(DP), INTENT(IN)  :: D, T, Ye, Ym, Enu
#ifdef EOSMODE_3D
    TYPE(EquationOfStateTableType), INTENT(in) :: EOSTable
#elif defined(EOSMODE_4D)
    TYPE(EquationOfState4DTableType), INTENT(in) :: EOSTable
#elif defined(EOSMODE_COMPOSE)
    TYPE(EquationOfStateCompOSETableType), INTENT(in) :: EOSTable
#endif
    TYPE(MuonTableType), INTENT(in) :: MuonTable
    REAL(DP), INTENT(OUT) :: OpElNuMu, OpElNuMuBar

    ! LOCAL VARIABLES
    REAL(DP) :: Yp, Xp, Xn
    REAL(DP) :: Mup, Mun, Mumu
    REAL(DP) :: Mp_eff, Mn_eff, Up, Un, n_n, n_p
    INTEGER  :: iEOS_Rho, iEOS_T, iEOS_Yp, iDV

    TYPE(MuonStateType) :: MuonState
    
    iEOS_Rho = EOSTable % TS % Indices % iRho
    iEOS_T = EOSTable % TS % Indices % iRho
    iEOS_Yp = EOSTable % TS % Indices % iYe ! Ye is actually Yp, but for consistency reasons it was not changed

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
          EOSTable % TS % States(iEOS_T) % Values,   &
          EOSTable % TS % States(iEOS_Yp) % Values,  &
          EOSTable % DV % Offsets(iDV), &
          EOSTable % DV % Variables(iDV) % Values(:,:,:), Xn )

    n_n = Xn * D * cm3fm3 / rmu
    n_p = Xp * D * cm3fm3 / rmu

    MuonState % t     = T
    MuonState % rhoym = D * Ym
    CALL FullMuonEOS(MuonTable, MuonState)

    Mumu = MuonState % mu

    CALL ElasticAbsorptionOpacity(T, Enu, Mup - mp, mp, Mp_eff, Up, n_p, &
      Mumu, mmu, Mun - mn, mn, Mn_eff, Un, n_n, OpElNuMu)

    CALL ElasticAbsorptionOpacity(T, Enu, Mun - mn, mn, Mn_eff, Un, n_n, &
      -Mumu, mmu, Mup - mp, mp, Mp_eff, Up, n_p, OpElNuMuBar)

  END SUBROUTINE ElasticAbsorptionOpacityNuMu

  SUBROUTINE ElasticAbsorptionOpacity(T, E1, Mu2, M2, M2_eff, U2, n2, &
    Mu3, M3, Mu4, M4, M4_eff, U4, n4, OpEl)

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

    ! WRITE(*,*) E3, E1, Mu2, Mu3, Mu4

  END SUBROUTINE ElasticAbsorptionOpacity


END MODULE wlCalculateAbEmOpacityModule

MODULE wlSoundSpeedModule
  
  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    kmev, rmu, kmev_inv, ergmev, cvel, cm3fm3, me, mmu
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: &
    PhotonGasType, LeptonGasType, &
    PhotonGasEOS , LeptonGasEOS , &
    GetPhotonLeptonGasEOS
  USE wlInterpolationUtilitiesModule, ONLY: &
    LinearInterpDeriv_Array_Point, &
    GetIndexAndDelta_Lin, GetIndexAndDelta_Log

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: CalculateSoundSpeed
  
  REAL(dp), PARAMETER :: ln10 = LOG(10.d0)
  
CONTAINS
  
  ! This one also calculates derivatives, but maybe you can provide derivatives ?
  SUBROUTINE CalculateSoundSpeed( D, T, Ye, Ym, D_T, T_T, Yp_T, &
    P_T, OS_P, E_T, OS_E, S_T, OS_S, &
    HelmTableElectrons, HelmTableMuons, Gamma, Cs, SeparateContributions)

    REAL(DP), INTENT(IN) :: D, T, Ye, Ym
    REAL(DP), INTENT(IN) :: D_T(1:), T_T(1:), Yp_T(1:)
    REAL(DP), INTENT(IN) :: P_T(1:,1:,1:), E_T(1:,1:,1:), S_T(1:,1:,1:)
  
    REAL(DP), INTENT(IN) :: OS_P, OS_E, OS_S
    LOGICAL, INTENT(IN)  :: SeparateContributions

    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons
    
    REAL(DP), INTENT(OUT)    :: Gamma, Cs

    REAL(DP) :: Pbary, Ebary, Sbary 
    REAL(DP) :: Pele , Eele , Sele
    REAL(DP) :: Pmu  , Emu  , Smu
    REAL(DP) :: Prad , Erad , Srad
    REAL(DP) :: Ptot , Etot , Stot

    REAL(DP) :: dPbarydD, dPbarydT
    REAL(DP) :: dEbarydD, dEbarydT
    REAL(DP) :: dSbarydD, dSbarydT

    REAL(DP) :: dPeledD , dPeledT
    REAL(DP) :: dEeledD , dEeledT
    REAL(DP) :: dSeledD , dSeledT
    
    REAL(DP) :: dPmudD  , dPmudT
    REAL(DP) :: dEmudD  , dEmudT
    REAL(DP) :: dSmudD  , dSmudT

    REAL(DP) :: dPraddD , dPraddT
    REAL(DP) :: dEraddD , dEraddT
    REAL(DP) :: dSraddD , dSraddT

    REAL(DP) :: dD, dT, dYp
    REAL(DP) :: aD, aT, aYp
    REAL(DP) :: Yp, Ye_over_Yp, Ym_over_Yp, h
    REAL(DP) :: D_cube, T_cube, Ye_cube, Ym_cube, P_cube, E_cube, S_cube
    REAL(DP) :: P_LeptPhot(2,2,2), E_LeptPhot(2,2,2), S_LeptPhot(2,2,2)
    REAL(DP) :: Ptot_T   (2,2,2), Etot_T(2,2,2)   , Stot_T(2,2,2)

    REAL(DP) :: dPdD, dPdT, dEdD, dEdT, dSdD, dSdT
    INTEGER  :: iD, iT, iYp, iL_D, iL_T, iL_Yp
    REAL(DP) :: eos_table(3,3,3,3), cs2_temp(3,3,3), gamma_temp(3,3,3)

    REAL(DP) :: LocalOffset

    TYPE(LeptonGasType) :: ElectronGasState
    TYPE(LeptonGasType) :: MuonGasState
    TYPE(PhotonGasType) :: PhotonGasState
      
    Yp = Ye + Ym

    IF (SeparateContributions) THEN
      ! PHOTONS
      PhotonGasState % T   = T
      PhotonGasState % rho = D
      
      CALL PhotonGasEOS(PhotonGasState)

      Prad = PhotonGasState % p
      Erad = PhotonGasState % e
      Srad = PhotonGasState % s 

      dPraddD = PhotonGasState % dpdr
      dPraddT = PhotonGasState % dpdt

      dEraddD = PhotonGasState % dedr
      dEraddT = PhotonGasState % dedt

      dSraddD = PhotonGasState % dsdr
      dSraddT = PhotonGasState % dsdt

      ! ELECTRONS
      ElectronGasState % t    = T
      ElectronGasState % rho  = D
      ElectronGasState % yL   = Ye

      ! calculate electron quantities
      CALL LeptonGasEOS(HelmTableElectrons, ElectronGasState)

      Pele = ElectronGasState % p
      Eele = ElectronGasState % e
      Sele = ElectronGasState % s 

      dPeledD = ElectronGasState % dpdr
      dPeledT = ElectronGasState % dpdt

      dEeledD = ElectronGasState % dedr
      dEeledT = ElectronGasState % dedt

      dSeledD = ElectronGasState % dsdr
      dSeledT = ElectronGasState % dsdt

      ! MUONS
      MuonGasState % t    = T
      MuonGasState % rho  = D
      MuonGasState % yL   = Ym

      ! calculate muon quantities
      CALL LeptonGasEOS(HelmTableMuons, MuonGasState)

      Pmu = MuonGasState % p
      Emu = MuonGasState % e
      Smu = MuonGasState % s 

      dPmudD = MuonGasState % dpdr
      dPmudT = MuonGasState % dpdt

      dEmudD = MuonGasState % dedr
      dEmudT = MuonGasState % dedt

      dSmudD = MuonGasState % dsdr
      dSmudT = MuonGasState % dsdt

      ! BARYONIC PART -----------------------------------------------------
      CALL GetIndexAndDelta_Log( D,  D_T,  iD,  dD  )
      CALL GetIndexAndDelta_Log( T,  T_T,  iT,  dT  )
      CALL GetIndexAndDelta_Lin( Yp, Yp_T, iYp, dYp )
      
      aD  = 1.0_dp / ( D * LOG10( D_T(iD+1) / D_T(iD) ) )
      aT  = 1.0_dp / ( T * LOG10( T_T(iT+1) / T_T(iT) ) )
      aYp = ln10   / ( Yp_T(iYp+1) - Yp_T(iYp) )
      
      LocalOffset = MINVAL( P_T(iD:iD+1,iT:iT+1,iYp:iYp+1) )
      IF (LocalOffset .lt. 0.0_dp) THEN
          LocalOffset = -1.1d0*LocalOffset
      ELSE
          LocalOffset = 0.0_dp
      ENDIF
      
      CALL LinearInterpDeriv_Array_Point &
          ( 1, 1, 1, dD, dT, dYp, aD, aT, aYp, LocalOffset, &
          LOG10( P_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + LocalOffset ), Pbary, &
          dPbarydD, dPbarydT )
      
      CALL LinearInterpDeriv_Array_Point &
          ( iD, iT, iYp, dD, dT, dYp, aD, aT, aYp, OS_E, E_T, Ebary, &
          dEbarydD, dEbarydT )

      CALL LinearInterpDeriv_Array_Point &
          ( iD, iT, iYp, dD, dT, dYp, aD, aT, aYp, OS_S, S_T, Sbary, &
          dSbarydD, dSbarydT )

      ! D and T in front take care of the derivative 
      ! wr2 logrho and logT. The rho multiplying the denominator in the second 
      ! one makes sure that you have erg/cm^3 instead of erg/g
      Gamma = (D*(dPbarydD + dPeledD + dPmudD + dPraddD) + &
            T*(dPbarydT + dPeledT + dPmudT + dPraddT)**2.0_DP / &
            (D*(dEbarydT + dEeledT + dEmudT + dEraddT)) ) / &
            (Pbary + Pele + Pmu + Prad)

      ! ! Another way of doing it
      ! Gamma = D*( (dPbarydD + dPeledD + dPmudD + dPraddD)   &
      !       -     (dSbarydD + dSeledD + dSmudD + dSraddD)   &
      !       *     (dPbarydT + dPeledT + dPmudT + dPraddT)   &
      !       /     (dSbarydT + dSeledT + dSmudT + dSraddT) ) &
      !       /     (Pbary + Pele + Pmu + Prad)

      ! relativistic definition with enthalpy
      h = (1.0_dp + (Ebary + Eele + Emu + Erad)/cvel**2 + (Pbary + Pele + Pmu + Prad)/D/cvel**2)
      
      Cs = SQRT(Gamma * (Pbary + Pele + Pmu + Prad) / (D*h))
      Cs = SQRT(Gamma * (Pbary + Pele + Pmu + Erad) / (D))

    ELSE

      Ye_over_Yp = Ye/Yp
      Ym_over_Yp = Ym/Yp

      CALL GetIndexAndDelta_Log( D,  D_T,  iD,  dD  )
      CALL GetIndexAndDelta_Log( T,  T_T,  iT,  dT  )
      CALL GetIndexAndDelta_Lin( Yp, Yp_T, iYp, dYp )
      
      aD = 1.0_dp / ( D * LOG10( D_T(iD+1) / D_T(iD) ) )
      aT = 1.0_dp / ( T * LOG10( T_T(iT+1) / T_T(iT) ) )
      aYp = ln10 / ( Yp_T(iYp+1) - Yp_T(iYp) )
    
      DO iL_T=1,2
        DO iL_D=1,2
          DO iL_Yp=1,2

            D_cube  = D_T(iD+iL_D-1)
            T_cube  = T_T(iT+iL_T-1)
            Ye_cube = Yp_T(iYp+iL_Yp-1) * Ye_over_Yp
            Ym_cube = Yp_T(iYp+iL_Yp-1) * Ym_over_Yp
            CALL GetPhotonLeptonGasEOS( &
              D_cube, T_cube, Ye_cube, Ym_cube, &
              HelmTableElectrons, HelmTableMuons, &
              P_cube, E_cube, S_cube)

            P_LeptPhot(iL_D,iL_T,iL_Yp) = P_cube
            E_LeptPhot(iL_D,iL_T,iL_Yp) = E_cube
            S_LeptPhot(iL_D,iL_T,iL_Yp) = S_cube

          END DO
        END DO
      END DO
    
      Ptot_T = LOG10(P_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + P_LeptPhot)
      CALL LinearInterpDeriv_Array_Point &
          ( 1, 1, 1, dD, dT, dYp, aD, aT, aYp, OS_P, Ptot_T, Ptot, &
          dPdD, dPdT )
      
      Etot_T = LOG10(10.00**E_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_LeptPhot )
      CALL LinearInterpDeriv_Array_Point &
          ( 1, 1, 1, dD, dT, dYp, aD, aT, aYp, OS_E, Etot_T, Etot, &
          dEdD, dEdT )

      Stot_T = LOG10(10.00**S_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_LeptPhot )
      CALL LinearInterpDeriv_Array_Point &
          ( 1, 1, 1, dD, dT, dYp, aD, aT, aYp, OS_S, Stot_T, Stot, &
          dSdD, dSdT )

      Gamma = (D*dPdD + T*dPdT**2.0_DP / &
          (D*dEdT) ) / Ptot

      ! ! Another way of doing it
      ! Gamma = D/Ptot * ( dPdD - (dSdD*dPdT/dSdT) )

      ! relativistic definition with enthalpy
      h = (1.0_dp + Etot/cvel**2 + Ptot/D/cvel**2)
      Cs = SQRT(Gamma * Ptot / (D*h))
      Cs = SQRT(Gamma * Ptot / (D))
      
      RETURN
      ! Below is how stellarcollapse does it...

      ! -----
      DO iL_T=-1,1
        DO iL_D=-1,1
          DO iL_Yp=-1,1

            D_cube  = D_T(iD+iL_D-1)
            T_cube  = T_T(iT+iL_T-1)
            Ye_cube = Yp_T(iYp+iL_Yp-1) * Ye_over_Yp
            Ym_cube = Yp_T(iYp+iL_Yp-1) * Ym_over_Yp
            CALL GetPhotonLeptonGasEOS( &
              D_cube, T_cube, Ye_cube, Ym_cube, &
              HelmTableElectrons, HelmTableMuons, &
              P_cube, E_cube, S_cube)
            
            eos_table(iL_D+2,iL_T+2,iL_Yp+2,1) = &
              LOG10(P_T(iD+iL_D,iT+iL_T,iYp+iL_Yp) + P_cube)
            eos_table(iL_D+2,iL_T+2,iL_Yp+2,2) = &
              LOG10(10.00**S_T(iD+iL_D,iT+iL_T,iYp+iL_Yp) + S_cube)
            eos_table(iL_D+2,iL_T+2,iL_Yp+2,3) = &
              LOG10(10.00**E_T(iD+iL_D,iT+iL_T,iYp+iL_Yp) + E_cube)
          END DO
        END DO
      END DO

      call derivatives_production(2, 3, 3, 3, &
          LOG10(D_T(iD-1:iD+1)), LOG10(T_T(iT-1:iT+1)), Yp_T(iYp-1:iYp+1), &
          eos_table, OS_E, cs2_temp, gamma_temp)
      Cs = SQRT(cs2_temp(2,2,2))
      ! WRITE(*,*) gamma_temp(2,2,2), Gamma
      Gamma = gamma_temp(2,2,2)

    ENDIF

  END SUBROUTINE CalculateSoundSpeed

! Bonus subroutines, currently unused
  ! This is taken from the stellarcollapse.org routines (Christian Ott's code)
  subroutine derivatives_production(igamma, nrho, ntemp, nye, logrho, logtemp, ye, &
    eos_table, energy_shift, cs2, gamma)
    
    integer, intent(in) :: igamma
    integer, intent(in) :: nrho, ntemp, nye
    real(dp), intent(in) :: eos_table(nrho,ntemp,nye,3)
    real(dp), intent(in) :: energy_shift
    real(dp), intent(in) :: logrho(nrho), logtemp(ntemp), ye(nye)
    real(dp), intent(out) :: cs2(nrho, ntemp, nye)
    real(dp), intent(out) :: gamma(nrho, ntemp, nye)
    
    integer i,j,k
    real(dp) :: dx,x1,x2,f1,f2,z,zz,h,cp,cv,beta_v,gamma_ad,kappa_t
    
    real(dp) :: dedT(nrho, ntemp, nye)
    real(dp) :: dsdlnT(nrho, ntemp, nye)
    real(dp) :: dsdlnrho(nrho, ntemp, nye)
    real(dp) :: dpdrhoT(nrho, ntemp, nye)
    real(dp) :: dedrhoT(nrho, ntemp, nye)
    real(dp) :: dpdT(nrho, ntemp, nye)
    real(dp) :: dpdrho(nrho, ntemp, nye)
    real(dp) :: dpde(nrho, ntemp, nye)
    integer  :: ipress = 1
    integer  :: ientropy = 2
    integer  :: ienergy = 3
    
    ! igamma: There are multiple ways to compute gamma1, via the entropy
    !         turns out to be pretty good.
      ! 1 -> Gamma1 = &
      !      dlnP/dlnrho|T,Y_e - ds/dlnrho|T,Y_e * (dlnP/dlnT / ds/dlnT)_rho,Y_e
      ! 2 -> Gamma1 = P/rho * ( dpdrho|e,Y_e + P/rho**2 dpde|rho,Y_e
    ! 3 -> Gamma1 = dlnP/dlnrho|T,Y_e + T * (dpdT|rho,Y_e)**2 / &
    !      (P*rho * dedT|rho,Y_e)
    
    !###########################################################################  
    ! dedT, e is erg in log10, T is in MeV in log 10
    dedT = 0.0d0
    
    do k=1,nye
      do i=1,nrho
        do j=2,ntemp-1
          x1 = logtemp(j-1)
          f1 = eos_table(i,j-1,k,ienergy)
          x2 = logtemp(j+1)
          f2 = eos_table(i,j+1,k,ienergy)
          dedT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)  &
          * 10.0d0**eos_table(i,j,k,ienergy)
        enddo
        
        ! boundaries: one-sided derivative
        x2 = logtemp(2)
        x1 = logtemp(1)
        f2 = eos_table(i,2,k,ienergy)
        f1 = eos_table(i,1,k,ienergy)
        dedT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1)  &
                * 10.0d0**eos_table(i,1,k,ienergy)
        x2 = logtemp(ntemp)
        x1 = logtemp(ntemp-1)
        f2 = eos_table(i,ntemp,k,ienergy)
        f1 = eos_table(i,ntemp-1,k,ienergy)
        dedT(i,ntemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(ntemp) &
        * 10.0d0**eos_table(i,ntemp,k,ienergy) 
      enddo
    enddo
    
    !########################################################################### 
    
    !dpdT, p is in dyn/cm^2 in log10, T is in MeV in log 10
    !dsdlnT, s is in k_B/baryon not in log 10, T is in meV in log 10
    dpdT = 0.0d0
    dsdlnT = 0.0d0
    
    do k=1,nye
      do i=1,nrho
        do j=2,ntemp-1
          x1 = logtemp(j-1)
          f1 = eos_table(i,j-1,k,ipress)
          x2 = logtemp(j+1)
          f2 = eos_table(i,j+1,k,ipress)
          dpdT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)  &
          * 10.0d0**eos_table(i,j,k,ipress)
          
          x1 = logtemp(j-1)
          f1 = log10(eos_table(i,j-1,k,ientropy))
          x2 = logtemp(j+1)
          f2 = log10(eos_table(i,j+1,k,ientropy))
          dsdlnT(i,j,k) = (f2-f1)/(x2-x1) * eos_table(i,j,k,ientropy)
          
        enddo
        
        ! boundaries: one-sided derivative
        x1 = logtemp(1)
        f1 = eos_table(i,1,k,ipress)
        x2 = logtemp(2)
        f2 = eos_table(i,2,k,ipress)
        dpdT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
        10.0d0**eos_table(i,1,k,ipress)
        
        x1 = logtemp(1)
        f1 = log10(eos_table(i,1,k,ientropy))
        x2 = logtemp(2)
        f2 = log10(eos_table(i,2,k,ientropy))
        dsdlnT(i,1,k) = (f2-f1)/(x2-x1) * eos_table(i,1,k,ientropy)
        
        
        x1 = logtemp(ntemp-1)
        f1 = eos_table(i,ntemp-1,k,ipress)
        x2 = logtemp(ntemp)
        f2 = eos_table(i,ntemp,k,ipress)
        dpdT(i,ntemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
        10.0d0**eos_table(i,1,k,ipress)
        
        x1 = logtemp(ntemp-1)
        f1 = log10(eos_table(i,ntemp-1,k,ientropy))
        x2 = logtemp(ntemp)
        f2 = log10(eos_table(i,ntemp,k,ientropy))
        dsdlnT(i,ntemp,k) = (f2-f1)/(x2-x1) & 
        * eos_table(i,ntemp,k,ientropy)
        
      enddo
    enddo
    
    !########################################################################### 
    ! dp/drho|T, p is in dyn/cm^2 in log 10, T is in MeV in log 10
      ! ds/dlnrho|T, s is in k_b/baryon not in log 10, T is in MeV in log 10
      ! de/drho|T, e is in erg in log 10, T is in MeV in log 10
      ! dp/drho|e = dp/drho|T + dp/dT * (-de/drho|T) / de/dT
    ! dp/de|rho = dp/dT / de/dT
    
    dpdrho   = 0.0d0
    dpde     = 0.0d0
    dsdlnrho = 0.0d0
    dpdrhoT  = 0.0d0
    dedrhoT  = 0.0d0
    
    do k=1,nye
      do j=1,ntemp
        do i=2,nrho-1
          x1 = logrho(i-1)
          x2 = logrho(i+1)
          f1 = eos_table(i-1,j,k,ipress)
          f2 = eos_table(i+1,j,k,ipress)
          dpdrhoT(i,j,k) = (f2-f1)/(x2-x1) /  10.0d0**logrho(i) &
          * 10.0d0**eos_table(i,j,k,ipress)
          
          x1 = logrho(i-1)
          x2 = logrho(i+1)
          f1 = eos_table(i-1,j,k,ienergy)
          f2 = eos_table(i+1,j,k,ienergy)
          dedrhoT(i,j,k) = (f2-f1)/(x2-x1) /  10.0d0**logrho(i) &
          * 10.0d0**eos_table(i,j,k,ienergy)
          
          x1 = logrho(i-1)
          x2 = logrho(i+1)
          f1 = log10(eos_table(i-1,j,k,ientropy))
          f2 = log10(eos_table(i+1,j,k,ientropy))
          dsdlnrho(i,j,k) = (f2-f1)/(x2-x1) * eos_table(i,j,k,ientropy)
        enddo
        
        ! boundaries: one-sided derivative
        x1 = logrho(1)
        x2 = logrho(2)
        f1 = eos_table(1,j,k,ipress)
        f2 = eos_table(2,j,k,ipress)
        dpdrhoT(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
        * 10.0d0**eos_table(1,j,k,ipress)
        
        f1 = eos_table(1,j,k,ienergy)
        f2 = eos_table(2,j,k,ienergy)
        dedrhoT(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
        * 10.0d0**eos_table(1,j,k,ienergy)
        
        x1 = logrho(nrho-1)
        x2 = logrho(nrho)
        f1 = eos_table(nrho-1,j,k,ipress)
        f2 = eos_table(nrho,j,k,ipress)
        dpdrhoT(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho) &
        * 10.0d0**eos_table(nrho,j,k,ipress)
        
        f1 = eos_table(nrho-1,j,k,ienergy)
        f2 = eos_table(nrho,j,k,ienergy)
        dedrhoT(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho) &
        * 10.0d0**eos_table(nrho,j,k,ienergy)
        
        x1 = logrho(1)
        x2 = logrho(2)
        f1 = log10(eos_table(1,j,k,ientropy))
        f2 = log10(eos_table(2,j,k,ientropy))
        dsdlnrho(1,j,k) = (f2-f1)/(x2-x1) * eos_table(1,j,k,ientropy)
        
        x1 = logrho(nrho-1)
        x2 = logrho(nrho)
        f1 = log10(eos_table(nrho-1,j,k,ientropy))
        f2 = log10(eos_table(nrho,j,k,ientropy))
        dsdlnrho(nrho,j,k) = (f2-f1)/(x2-x1) * eos_table(nrho,j,k,ientropy)
        
      enddo
    enddo
    
    do k=1,nye
      do j=1,ntemp
        do i=1,nrho
          dpdrho(i,j,k) = dpdrhoT(i,j,k) + dpdT(i,j,k) * &
          (-dedrhoT(i,j,k))/dedT(i,j,k)
          
          dpde(i,j,k) = dpdT(i,j,k) / dedT(i,j,k)
        enddo
      enddo
    enddo
    
    !########################################################################### 
    cs2 = 0.0d0
    gamma = 0.0d0
    
    do k=1,nye
      do j=1,ntemp
        do i=1,nrho
          z = 10.0d0**logrho(i)/10.0d0**eos_table(i,j,k,ipress)
          zz = 10.0d0**logtemp(j)/10.0d0**eos_table(i,j,k,ipress)
          
          if(igamma .eq. 1) then
            gamma(i,j,k) = z * dpdrhoT(i,j,k) & 
            - dsdlnrho(i,j,k) * zz * dpdT(i,j,k) / dsdlnT(i,j,k)
          endif
          
          if(igamma .eq. 2) then
            gamma(i,j,k) =  dpdrho(i,j,k) + 10.0d0**eos_table(i,j,k,ipress) / &
            (10.0d0**logrho(i))**2 * dpde(i,j,k)
            gamma(i,j,k) = z * gamma(i,j,k)
          endif
          
          if(igamma .eq. 3) then
            gamma(i,j,k) = 10.0d0**logrho(i)/10.0d0**eos_table(i,j,k,ipress) &
            * dpdrhoT(i,j,k) &
            + 10.0d0**logtemp(j)*dpdT(i,j,k)**2 / &
            (10.0d0**eos_table(i,j,k,ipress)*10.0d0**logrho(i) * dedT(i,j,k))
          endif
          
          cs2(i,j,k) = gamma(i,j,k) / z
          cs2(i,j,k) = gamma(i,j,k) * 10.0d0**eos_table(i,j,k,ipress) / &
          (10.0d0**eos_table(i,j,k,ipress) + (10.0d0**eos_table(i,j,k,ienergy) - energy_shift - &
          0.511d0 / rmu * ergmev * ye(k) + cvel**2) * 10.0d0**logrho(i) ) * cvel**2
          
          if(igamma .eq. 5) then
            h = 10.0d0**eos_table(i,j,k,ipress) + (10.0d0**eos_table(i,j,k,ienergy) - energy_shift - &
              0.511d0 / rmu * ergmev * ye(k) + cvel**2) * 10.0d0**logrho(i)
            
            kappa_t = 10.0d0**logrho(i) * dpdrho(i,j,k) * (cm3fm3/ergmev)
            kappa_t = 1.0d0 / kappa_t
            beta_v = dpdT(i,j,k) * (cm3fm3/ergmev) / kmev
            cv = 1.0d0 / (10.0d0**logrho(i)/rmu*cm3fm3) * dsdlnT(i,j,k) / LOG(10.0d0) 
            cp = cv + 10.0d0**logtemp(j)*kmev/(10.0d0**logrho(i)/rmu*cm3fm3) * beta_v * kappa_t * beta_v
            gamma_ad = cp/cv
                          
            cs2(i,j,k) = gamma_ad / (h * (cm3fm3/ergmev) * kappa_t )
            gamma(i,j,k) = gamma_ad / 10.0d0**eos_table(i,j,k,ipress)*ergmev/cm3fm3 / kappa_t
            cs2 = cs2 * cvel**2
            
          endif
        enddo
      enddo
    enddo

  end subroutine derivatives_production

END MODULE wlSoundSpeedModule

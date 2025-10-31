MODULE wlMuonEOS
    
    USE wlKindModule, ONLY: dp
    USE wlInterpolationUtilitiesModule, ONLY: &
        Index1D_Lin, &
        Index1D_Log
    USE wlLeptonEOSModule, ONLY: &
        MuonTableType
        
    IMPLICIT NONE
    PRIVATE
    
    TYPE, PUBLIC :: MuonStateType
        
        REAL(dp) :: t
        REAL(dp) :: rhoym
        REAL(dp) :: rho
        REAL(dp) :: mu
        REAL(dp) :: p
        REAL(dp) :: e
        REAL(dp) :: s
        REAL(dp) :: dlnPdlnrho
        REAL(dp) :: dlnPdlnT
        REAL(dp) :: dlnsdlnrho
        REAL(dp) :: dlnsdlnT
        REAL(dp) :: dlnedlnrho
        REAL(dp) :: dlnedlnT
        
    END TYPE MuonStateType
    
    PUBLIC :: FullMuonEOS
    
CONTAINS
    
    ! INTERPOLATION ROUTINES FROM TOBIAS
    SUBROUTINE FullMuonEOS(MuonTable, MuonState, CalculateDerivatives_Option)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
        TYPE(MuonTableType), INTENT(IN) :: MuonTable
        TYPE (MuonStateType), INTENT(INOUT) :: MuonState
        LOGICAL, DIMENSION(3), OPTIONAL, INTENT(IN)  :: CalculateDerivatives_Option
        LOGICAL :: CalculatePressureDerivatives, CalculateEntropyDerivatives, &
                   CalculateEnergyDerivatives
        
        !----locals ------------------------------------------------------------
        integer :: i
        
        real(dp) :: temp, rhoym, logt, logrhoym
        real(dp) :: dt,dd
        real(dp) :: rd,rt
        real(dp) :: sign_ym
        integer :: iT, iDen, iT_Max, iDen_Max
        
        IF ( PRESENT(CalculateDerivatives_Option) ) THEN
            CalculatePressureDerivatives = CalculateDerivatives_Option(1)
            CalculateEnergyDerivatives = CalculateDerivatives_Option(2)
            CalculateEntropyDerivatives = CalculateDerivatives_Option(3)
        ELSE
            CalculatePressureDerivatives = .FALSE.
            CalculateEnergyDerivatives = .FALSE.
            CalculateEntropyDerivatives = .FALSE.
        END IF
        
        !....exit early if the density is below a threshold that is consistent 
        !....with the EmAb on nucleons opacities .............................
        IF ( MuonState % rho < MuonTable % eos_minD ) THEN 
            MuonState % p  = 0.0d0
            MuonState % e  = 0.0d0
            MuonState % s  = 0.0d0
            MuonState % mu = 0.0d0
            RETURN
        END If


        !....convert input variables to log scale .............................
        temp  = MuonState % t
        rhoym = MuonState % rhoym
        sign_ym = SIGN(1.0d0, rhoym)
        rhoym = ABS(rhoym) 

        iT_Max = MuonTable % nPointsTemp
        iDen_Max = MuonTable % nPointsDen

        IF(  temp  < MuonTable % t(1) &
        .OR. temp  > MuonTable % t(iT_Max) ) THEN
            MuonState % p  = 0.0d0
            MuonState % e  = 0.0d0
            MuonState % s  = 0.0d0
            MuonState % mu = 0.0d0
!write(*,*) 'muon eos out of T table bounds'
!write(*,*) temp, MuonTable % t(1), MuonTable % t(iT_Max)
            RETURN
        ENDIF

        iT = Index1D_Log( temp, MuonTable % t(:) )

        IF(  rhoym  < MuonTable % rhoym(iT,1) &
        .OR. rhoym  > MuonTable % rhoym(iT,iDen_Max) ) THEN
            MuonState % p  = 0.0d0
            MuonState % e  = 0.0d0
            MuonState % s  = 0.0d0
            MuonState % mu = 0.0d0
!write(*,*) 'muon eos out of rhoym table bounds'
!write(*,*) rhoym, MuonTable % rhoym(iT,1), MuonTable % rhoym(iT,iDen_Max)
            RETURN
        ENDIF

        logt  = log10(temp)
        logrhoym = log10(rhoym)
        
        DO i = 1, iDen_Max
          IF (MuonTable % rhoym(iT,i) >= rhoym) THEN
            iDen = i
            EXIT
          ENDIF 
        END DO
        
        !.....determine interpolation weights ..................................
        dt = LOG10( MuonTable % t(iT+1) ) - LOG10( MuonTable % t(iT) )
        dd = LOG10( MuonTable % rhoym(iT,iDen+1) ) - LOG10( MuonTable % rhoym(iT,iDen) )
        
        rt = (logt - LOG10( MuonTable % t(iT) ))/dt
        rd = (logrhoym - LOG10( MuonTable % rhoym(iT,iDen) ))/dd
                
        !.....interpolation scheme.............................................
        MuonState % mu = & ! muon chemical potential (including rest mass)  [MeV]
            (1.-rt)*(1.-rd) * MuonTable % mu(iT  ,iDen  ) &
             +  rt *(1.-rd) * MuonTable % mu(iT+1,iDen  ) &
             +     rt * rd  * MuonTable % mu(iT+1,iDen+1) &
             + (1.-rt)* rd  * MuonTable % mu(iT  ,iDen+1)
        MuonState % mu = MuonState % mu * sign_ym
             
        MuonState % p = & ! pressure              [MeV/fm^3]
            (1.-rt)*(1.-rd) * MuonTable % p(iT  ,iDen  ) &
             +  rt *(1.-rd) * MuonTable % p(iT+1,iDen  ) &
             +     rt * rd  * MuonTable % p(iT+1,iDen+1) &
             + (1.-rt)* rd  * MuonTable % p(iT  ,iDen+1)
             
        MuonState % e = & ! internal energy       [MeV/fm^3]
            (1.-rt)*(1.-rd) * MuonTable % e(iT  ,iDen  ) &
             +  rt *(1.-rd) * MuonTable % e(iT+1,iDen  ) &
             +     rt * rd  * MuonTable % e(iT+1,iDen+1) &
             + (1.-rt)* rd  * MuonTable % e(iT  ,iDen+1)
             
        MuonState % s = & ! entropy
            (1.-rt)*(1.-rd) * MuonTable % s(iT  ,iDen  ) &
             +  rt *(1.-rd) * MuonTable % s(iT+1,iDen  ) &
             +     rt * rd  * MuonTable % s(iT+1,iDen+1) &
             + (1.-rt)* rd  * MuonTable % s(iT  ,iDen+1)
             
        IF (CalculatePressureDerivatives) THEN
            MuonState % dlnPdlnT = & ! Pressure derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnPdlnT(iT  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnPdlnT(iT+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnPdlnT(iT+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnPdlnT(iT  ,iDen+1)         

            MuonState % dlnPdlnrho = & ! Pressure derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnPdlnrho(iT  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnPdlnrho(iT+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnPdlnrho(iT+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnPdlnrho(iT  ,iDen+1)       
        
        ENDIF
             
        IF (CalculateEnergyDerivatives) THEN
            MuonState % dlnedlnT = & ! Energy derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnedlnT(iT  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnedlnT(iT+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnedlnT(iT+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnedlnT(iT  ,iDen+1)         

            MuonState % dlnedlnrho = & ! Energy derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnedlnrho(iT  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnedlnrho(iT+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnedlnrho(iT+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnedlnrho(iT  ,iDen+1)       
        
        ENDIF

        IF (CalculateEntropyDerivatives) THEN
            MuonState % dlnsdlnT = & ! Entropy derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnsdlnT(iT  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnsdlnT(iT+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnsdlnT(iT+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnsdlnT(iT  ,iDen+1)         

            MuonState % dlnsdlnrho = & ! Entropy derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnsdlnrho(iT  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnsdlnrho(iT+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnsdlnrho(iT+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnsdlnrho(iT  ,iDen+1)       
        
        ENDIF

        END SUBROUTINE FullMuonEOS

END MODULE wlMuonEOS    

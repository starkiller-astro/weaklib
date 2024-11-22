MODULE wlMuonEOS
    
    USE wlKindModule, ONLY: dp
    USE wlInterpolationUtilitiesModule, ONLY: &
        Index1D_Lin, &
        Index1D_Log
    USE wlLeptonEOSModule, ONLY: &
        MuonEOSType
        
    IMPLICIT NONE
    PRIVATE
    
    TYPE, PUBLIC :: MuonStateType
        
        REAL(dp) :: t
        REAL(dp) :: rhoym
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
        
        TYPE(MuonEOSType), INTENT(IN) :: MuonTable
        TYPE (MuonStateType), INTENT(INOUT) :: MuonState
        LOGICAL, DIMENSION(3), OPTIONAL, INTENT(IN)  :: CalculateDerivatives_Option
        LOGICAL :: CalculatePressureDerivatives, CalculateEntropyDerivatives, &
                   CalculateEnergyDerivatives
        
        !----locals ------------------------------------------------------------
        integer :: i
        
        real(dp) :: temp, rhoym, logt, logrhoymu
        real(dp) :: dt,dd
        real(dp) :: rd,rt
        integer :: it, it_1, iDen, it_max, iDen_max
        
        IF ( PRESENT(CalculateDerivatives_Option) ) THEN
            CalculatePressureDerivatives = CalculateDerivatives_Option(1)
            CalculateEnergyDerivatives = CalculateDerivatives_Option(2)
            CalculateEntropyDerivatives = CalculateDerivatives_Option(3)
        ELSE
            CalculatePressureDerivatives = .FALSE.
            CalculateEnergyDerivatives = .FALSE.
            CalculateEntropyDerivatives = .FALSE.
        END IF
        
        !....convert input variables to log scale .............................
        temp  = MuonState % t
        rhoym = MuonState % rhoym

        IF(  rhoym < MuonTable % rhoym(1) &
        .OR. temp  < MuonTable % t(1) &
        .OR. rhoym > MuonTable % rhoym(ubound(MuonTable % rhoym,1)) &
        .OR. temp  > MuonTable % t    (ubound(MuonTable % t,    1)) ) THEN
            MuonState % p  = 0.0d0
            MuonState % e  = 0.0d0
            MuonState % s  = 0.0d0
            MuonState % mu = 0.0d0
            RETURN
        ENDIF

        logt  = log10(temp)
        logrhoymu = log10(rhoym)
        
        it_max = MuonTable % nPointsTemp
        iDen_max = MuonTable % nPointsDen

        it = Index1D_Log( temp, MuonTable % t(:) )

        ! At these low temperatures muons don't really matter
        IF (it .eq. 1) THEN
            MuonState % p = 0.0d0
            MuonState % e = 0.0d0
            MuonState % s = 0.0d0
            MuonState % mu = 0.0d0
            RETURN
        ENDIF
        
        ! This is only if the grid is not log
        ! ! .....find t-index......................................................
        ! IF (temp .le. MuonTable % t(1)) THEN
            ! ! At these low temperatures muons don't really matter
            ! MuonState % p = 0.0d0
            ! MuonState % e = 0.0d0
            ! MuonState % s = 0.0d0
            ! MuonState % mu = 0.0d0
            ! RETURN
            
        ! ELSE IF (temp .ge. MuonTable % t(it_max)) THEN
            ! it = it_max
        ! ELSE
            ! DO i=1,it_max-1
                ! IF ( (temp .ge. MuonTable % t(i)) .and. (temp .lt. MuonTable % t(i+1)) ) THEN
                    ! it = i
                    ! EXIT
                ! END IF
            ! END DO
        ! END IF
        
        ! IF (it .ne. it_1) THEN
            ! WRITE(*,*) it, it_1
            ! STOP
        ! ENDIF
        
        iDen = Index1D_Log( rhoym, MuonTable % rhoym(:) )

        ! At these low densities muons don't really matter
        IF (iDen .eq. 1) THEN
            MuonState % p = 0.0d0
            MuonState % e = 0.0d0
            MuonState % s = 0.0d0
            MuonState % mu = 0.0d0
            RETURN
        ENDIF
        
        ! This is only if the grid is not log
        ! ! .....find d-index......................................................
        ! IF (rhoym .le. MuonTable % rhoym(it,1)) THEN
            ! ! At these low densities muons don't really matter
            ! MuonState % p = 0.0d0
            ! MuonState % e = 0.0d0
            ! MuonState % s = 0.0d0
            ! MuonState % mu = 0.0d0
            ! RETURN
        ! ELSE IF (rhoym .ge. MuonTable % rhoym(it,iDen_max)) THEN
            ! iDen = iDen_max
        ! ELSE
            ! DO i=iDen_max-1,1,-1
                ! IF((rhoym .ge. MuonTable % rhoym(it,i)) .and. &
                ! (rhoym .lt. MuonTable % rhoym(it,i+1)))THEN
                    ! iDen = i
                    ! exit
                ! END IF
            ! END DO
        ! END IF
        
        !.....determine interpolation weights ..................................
        dt = LOG10( MuonTable % t(it+1) ) - LOG10( MuonTable % t(it) )
        dd = LOG10( MuonTable % rhoym(iDen+1) ) - LOG10( MuonTable % rhoym(iDen) )
        
        rt = (logt - LOG10( MuonTable % t(it) ))/dt
        rd = (logrhoymu - LOG10( MuonTable % rhoym(iDen) ))/dd
                
        !.....interpolation scheme.............................................
        MuonState % mu = & ! muon chemical potential  [MeV]
            (1.-rt)*(1.-rd) * MuonTable % mu(it  ,iDen  ) &
             +  rt *(1.-rd) * MuonTable % mu(it+1,iDen  ) &
             +     rt * rd  * MuonTable % mu(it+1,iDen+1) &
             + (1.-rt)* rd  * MuonTable % mu(it  ,iDen+1)
             
        MuonState % p = & ! pressure              [MeV/fm^3]
            (1.-rt)*(1.-rd) * MuonTable % p(it  ,iDen  ) &
             +  rt *(1.-rd) * MuonTable % p(it+1,iDen  ) &
             +     rt * rd  * MuonTable % p(it+1,iDen+1) &
             + (1.-rt)* rd  * MuonTable % p(it  ,iDen+1)
             
        MuonState % e = & ! internal energy       [MeV/fm^3]
            (1.-rt)*(1.-rd) * MuonTable % e(it  ,iDen  ) &
             +  rt *(1.-rd) * MuonTable % e(it+1,iDen  ) &
             +     rt * rd  * MuonTable % e(it+1,iDen+1) &
             + (1.-rt)* rd  * MuonTable % e(it  ,iDen+1)
             
        MuonState % s = & ! entropy
            (1.-rt)*(1.-rd) * MuonTable % s(it  ,iDen  ) &
             +  rt *(1.-rd) * MuonTable % s(it+1,iDen  ) &
             +     rt * rd  * MuonTable % s(it+1,iDen+1) &
             + (1.-rt)* rd  * MuonTable % s(it  ,iDen+1)
             
        IF (CalculatePressureDerivatives) THEN
            MuonState % dlnPdlnT = & ! Pressure derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnPdlnT(it  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnPdlnT(it+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnPdlnT(it+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnPdlnT(it  ,iDen+1)         

            MuonState % dlnPdlnrho = & ! Pressure derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnPdlnrho(it  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnPdlnrho(it+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnPdlnrho(it+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnPdlnrho(it  ,iDen+1)       
        
        ENDIF
             
        IF (CalculateEnergyDerivatives) THEN
            MuonState % dlnedlnT = & ! Energy derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnedlnT(it  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnedlnT(it+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnedlnT(it+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnedlnT(it  ,iDen+1)         

            MuonState % dlnedlnrho = & ! Energy derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnedlnrho(it  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnedlnrho(it+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnedlnrho(it+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnedlnrho(it  ,iDen+1)       
        
        ENDIF

        IF (CalculateEntropyDerivatives) THEN
            MuonState % dlnsdlnT = & ! Entropy derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnsdlnT(it  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnsdlnT(it+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnsdlnT(it+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnsdlnT(it  ,iDen+1)         

            MuonState % dlnsdlnrho = & ! Entropy derivatives
                (1.-rt)*(1.-rd) * MuonTable % dlnsdlnrho(it  ,iDen  ) &
                 +  rt *(1.-rd) * MuonTable % dlnsdlnrho(it+1,iDen  ) &
                 +     rt * rd  * MuonTable % dlnsdlnrho(it+1,iDen+1) &
                 + (1.-rt)* rd  * MuonTable % dlnsdlnrho(it  ,iDen+1)       
        
        ENDIF

        END SUBROUTINE FullMuonEOS

END MODULE wlMuonEOS    

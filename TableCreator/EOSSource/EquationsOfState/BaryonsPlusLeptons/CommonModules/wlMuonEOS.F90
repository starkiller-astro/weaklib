MODULE wlMuonEOS
	
	USE wlKindModule, ONLY: dp
	USE wlExtNumericalModule, ONLY: zero, one, pi, half
	USE wlInterpolationUtilitiesModule, ONLY: &
		Index1D_Lin, &
		Index1D_Log
	USE wlLeptonEOSModule, ONLY: &
		MuonEOSType
		
	IMPLICIT NONE
	PRIVATE
	
	TYPE, PUBLIC :: MuonStateType
		
		REAL(dp) :: t
		REAL(dp) :: rhoymu
		REAL(dp) :: mu
		REAL(dp) :: p
		REAL(dp) :: e
		REAL(dp) :: s
		
	END TYPE MuonStateType
	
	PUBLIC :: FullMuonEOS
	
CONTAINS
	
	! INTERPOLATION ROUTINES FROM TOBIAS
	SUBROUTINE FullMuonEOS(MuonTable, MuonState)
		
		TYPE(MuonEOSType), INTENT(IN) :: MuonTable
        TYPE (MuonStateType), INTENT(INOUT) :: MuonState
		
		!----locals ------------------------------------------------------------
		integer :: i
		
		real(dp) :: temp, rhoymu, logt, logrhoymu
		real(dp) :: dt,dd
		real(dp) :: rd,rt
		integer :: it, it_1, iDen, it_max, iDen_max
		
		!....convert input variables to log scale .............................
		temp = MuonState % t
		rhoymu = MuonState % rhoymu
		logt = log10(temp)
		logrhoymu = log10(rhoymu)
		
		it_max = MuonTable % nPointsTemp
		iDen_max = MuonTable % nPointsMu

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
		
		iDen = Index1D_Log( rhoymu, MuonTable % rhoymu(:) )

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
		! IF (rhoymu .le. MuonTable % rhoymu(it,1)) THEN
			! ! At these low densities muons don't really matter
			! MuonState % p = 0.0d0
			! MuonState % e = 0.0d0
			! MuonState % s = 0.0d0
			! MuonState % mu = 0.0d0
			! RETURN
		! ELSE IF (rhoymu .ge. MuonTable % rhoymu(it,iDen_max)) THEN
			! iDen = iDen_max
		! ELSE
			! DO i=iDen_max-1,1,-1
				! IF((rhoymu .ge. MuonTable % rhoymu(it,i)) .and. &
				! (rhoymu .lt. MuonTable % rhoymu(it,i+1)))THEN
					! iDen = i
					! exit
				! END IF
			! END DO
		! END IF
		
		!.....determine interpolation weights ..................................
		dt = LOG10( MuonTable % t(it+1) ) - LOG10( MuonTable % t(it) )
		dd = LOG10( MuonTable % rhoymu(iDen+1) ) - LOG10( MuonTable % rhoymu(iDen) )
		
		rt = (logt - LOG10( MuonTable % t(it) ))/dt
		rd = (logrhoymu - LOG10( MuonTable % rhoymu(iDen) ))/dd
				
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
			 
	END SUBROUTINE FullMuonEOS

END MODULE wlMuonEOS	
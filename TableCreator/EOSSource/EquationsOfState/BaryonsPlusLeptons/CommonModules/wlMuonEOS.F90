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
		
		real(dp) :: temp, rhoymu, logt
		real(dp) :: dt,dd
		real(dp) :: rd,rt
		integer :: it, it_1, imu, it_max, imu_max
		
		!....convert input variables to log scale .............................
		temp = MuonState % t
		rhoymu = MuonState % rhoymu
		logt = log10(temp)
		
		it_max = MuonTable % nPointsTemp
		imu_max = MuonTable % nPointsMu

		it = Index1D_Log( temp, MuonTable % t(:) )

		! At these low temperatures muons don't really matter
		IF (it .eq. 1) THEN
			MuonState % p = 0.0d0
			MuonState % e = 0.0d0
			MuonState % s = 0.0d0
			MuonState % mu = 0.0d0
			RETURN
		ENDIF
		
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
		
		!.....find d-index......................................................
		IF (rhoymu .le. MuonTable % rhoymu(it,1)) THEN
			! At these low densities muons don't really matter
			MuonState % p = 0.0d0
			MuonState % e = 0.0d0
			MuonState % s = 0.0d0
			MuonState % mu = 0.0d0
			RETURN
		ELSE IF (rhoymu .ge. MuonTable % rhoymu(it,imu_max)) THEN
			imu = imu_max
		ELSE
			DO i=imu_max-1,1,-1
				IF((rhoymu .ge. MuonTable % rhoymu(it,i)) .and. &
				(rhoymu .lt. MuonTable % rhoymu(it,i+1)))THEN
					imu = i
					exit
				END IF
			END DO
		END IF
		
		!.....determine interpolation weights ..................................
		dt = LOG10( MuonTable % t(it+1) ) - LOG10( MuonTable % t(it) )
		dd = MuonTable % rhoymu(it,imu+1) - MuonTable % rhoymu(it,imu)
		
		rt = (logt - LOG10( MuonTable % t(it) ))/dt
		rd = (rhoymu - MuonTable % rhoymu(it,imu))/dd
				
		!.....interpolation scheme.............................................
		MuonState % mu = & ! muon chemical potential  [MeV]
					MuonTable % mu(imu  ) &
			 + rt * MuonTable % mu(imu+1) &
			 - rt * MuonTable % mu(imu  )

		MuonState % p = & ! pressure              [MeV/fm^3]
			(1.-rt)*(1.-rd) * MuonTable % p(it  ,imu  ) &
			 +  rt *(1.-rd) * MuonTable % p(it+1,imu  ) &
			 +     rt * rd  * MuonTable % p(it+1,imu+1) &
			 + (1.-rt)* rd  * MuonTable % p(it  ,imu+1)
			 
		MuonState % e = & ! internal energy       [MeV/fm^3]
			(1.-rt)*(1.-rd) * MuonTable % e(it  ,imu  ) &
			 +  rt *(1.-rd) * MuonTable % e(it+1,imu  ) &
			 +     rt * rd  * MuonTable % e(it+1,imu+1) &
			 + (1.-rt)* rd  * MuonTable % e(it  ,imu+1)
			 
		MuonState % s = & ! entropy
			(1.-rt)*(1.-rd) * MuonTable % s(it  ,imu  ) &
			 +  rt *(1.-rd) * MuonTable % s(it+1,imu  ) &
			 +     rt * rd  * MuonTable % s(it+1,imu+1) &
			 + (1.-rt)* rd  * MuonTable % s(it  ,imu+1)
			 
	END SUBROUTINE FullMuonEOS

END MODULE wlMuonEOS	
MODULE wlMuonEOS
	
	USE wlKindModule, ONLY: dp
	USE wlExtNumericalModule, ONLY: zero, one, pi, half
	USE wlLeptonEOSModule
	
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
		
		real :: temp, rhoymu, logt
		real :: dt,dd
		real :: rd,rt
		integer :: it,imu
		
		!....convert input variables to log scale .............................
		temp = MuonState % t
		rhoymu = MuonState % rhoymu
		logt = log10(temp)
		
		it_max = MuonTable % nPointsTemp
		imu_max = MuonTable % nPointsMu

		!.....find t-index......................................................
		IF (temp .le. MuonTable % t(1)) THEN
			it = 1
		ELSE IF (temp .ge. MuonTable % t(it_max)) THEN
			it = it_max
		ELSE
			DO i=1,it_max-1
				IF ( (temp .ge. MuonTable % t(i)) .and. (temp .lt. MuonTable % t(i+1)) ) THEN
					it = i
					EXIT
				END IF
			END DO
		END IF
		
		!.....find d-index......................................................
		IF (rhoymu .le. MuonTable % rhoymu(it,1)) THEN
			imu = 1
		ELSE IF (rhoymu .ge. MuonTable % rhoymu(it,imu_max)) THEN
			imu = imu_max
		ELSE
			DO i=1,imu_max-1
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
		
				REAL(dp) :: mu
		REAL(dp) :: p
		REAL(dp) :: e
		REAL(dp) :: s
		
		!.....interpolation scheme.............................................
		MuonState % mu = & ! muon chemical potential  [MeV]
			(1.-rt)*(1.-rd) * MuonTable % mu(it  ,imu  ) &
			 +  rt *(1.-rd) * MuonTable % mu(it+1,imu  ) &
			 +     rt * rd  * MuonTable % mu(it+1,imu+1) &
			 + (1.-rt)* rd  * MuonTable % mu(it  ,imu+1)

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
			 
	END SUBROUTINE InterpolateMuon

END MODULE wlMuonEOS	
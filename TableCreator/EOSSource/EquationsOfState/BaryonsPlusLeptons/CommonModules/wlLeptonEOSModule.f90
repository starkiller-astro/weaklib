MODULE wlLeptonEOSModule
	
	USE wlKindModule, ONLY: dp
	USE wlExtNumericalModule, ONLY: zero, one, pi
	
	IMPLICIT NONE
	PRIVATE
	
    INTEGER, PARAMETER, PUBLIC :: iTempMax=541, iDenMax=201 
    INTEGER, PARAMETER, PUBLIC :: nTempMuon=270, nDenMuon=1301 
	
 	TYPE, PUBLIC :: HelmholtzEOSType
		
		INTEGER :: nPointsDen !imax
		INTEGER :: nPointsTemp !jmax
		
		REAL(dp), DIMENSION(:), ALLOCATABLE :: t
		REAL(dp), DIMENSION(:), ALLOCATABLE :: d
		
		! Free Energy Table
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: f
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: fd
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ft
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: fdd
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ftt
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: fdt
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: fddt
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: fdtt
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: fddtt
		
		!  pressure derivative with density table
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dpdf
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dpdfd
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dpdft
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dpdfdt
		
		!  electron chemical potential table
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ef
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: efd
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: eft
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: efdt
		
		!  number density table
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xf
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xfd
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xft
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xfdt
		
		!  deltas needed for interpolation
		REAL(dp), DIMENSION(:), ALLOCATABLE :: dt
		REAL(dp), DIMENSION(:), ALLOCATABLE :: dt2
		REAL(dp), DIMENSION(:), ALLOCATABLE :: dti
		REAL(dp), DIMENSION(:), ALLOCATABLE :: dt2i
		
		REAL(dp), DIMENSION(:), ALLOCATABLE :: dd
		REAL(dp), DIMENSION(:), ALLOCATABLE :: dd2
		REAL(dp), DIMENSION(:), ALLOCATABLE :: ddi
		REAL(dp), DIMENSION(:), ALLOCATABLE :: dd2i
		
		! minimum and maximum possible densities
		REAL(dp) :: mintemp
		REAL(dp) :: maxtemp
		REAL(dp) :: mindens
		REAL(dp) :: maxdens
		
	END TYPE HelmholtzEOSType
	
	TYPE, PUBLIC :: MuonEOSType
			
		INTEGER :: nPointsMu
		INTEGER :: nPointsTemp
		
		REAL(dp), DIMENSION(:), ALLOCATABLE :: t
		REAL(dp), DIMENSION(:), ALLOCATABLE :: mu
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rhoymu
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: e
		REAL(dp), DIMENSION(:,:), ALLOCATABLE :: s

	END TYPE MuonEOSType

	PUBLIC AllocateHelmEOS
	PUBLIC DeAllocateHelmEOS
	PUBLIC ReadHelmEOSdat

	PUBLIC AllocateMuonEOS
	PUBLIC DeAllocateMuonEOS
	PUBLIC ReadMuonEOSdat
	
CONTAINS

	SUBROUTINE AllocateHelmEOS( HelmholtzEOS, nPoints )
		
		TYPE(HelmholtzEOSType)      :: HelmholtzEOS 
		INTEGER, DIMENSION(2), INTENT(IN) :: nPoints
				
		HelmholtzEOS % nPointsDen   = nPoints(1)
		HelmholtzEOS % nPointsTemp  = nPoints(2)
		
		ALLOCATE( HelmholtzEOS % d( nPoints(1) ) )
		ALLOCATE( HelmholtzEOS % t( nPoints(2) ) )
		
		ALLOCATE( HelmholtzEOS % f( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % fd( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % ft( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % fdd( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % ftt( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % fdt( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % fddt( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % fdtt( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % fddtt( nPoints(1), nPoints(2) ) )
		
		ALLOCATE( HelmholtzEOS % dpdf( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % dpdfd( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % dpdft( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % dpdfdt( nPoints(1), nPoints(2) ) )
		
		ALLOCATE( HelmholtzEOS % ef( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % efd( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % eft( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % efdt( nPoints(1), nPoints(2) ) )
		
		ALLOCATE( HelmholtzEOS % xf( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % xfd( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % xft( nPoints(1), nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % xfdt( nPoints(1), nPoints(2) ) )
		
		ALLOCATE( HelmholtzEOS % dt( nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % dt2( nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % dti( nPoints(2) ) )
		ALLOCATE( HelmholtzEOS % dt2i( nPoints(2) ) )
		
		ALLOCATE( HelmholtzEOS % dd( nPoints(1) ) )
		ALLOCATE( HelmholtzEOS % dd2( nPoints(1) ) )
		ALLOCATE( HelmholtzEOS % ddi( nPoints(1) ) )
		ALLOCATE( HelmholtzEOS % dd2i( nPoints(1) ) )
		
	END SUBROUTINE AllocateHelmEOS
	
	SUBROUTINE DeAllocateHelmEOS( HelmholtzEOS )
	
		TYPE(HelmholtzEOSType)      :: HelmholtzEOS 

		DeAllocate( HelmholtzEOS % t )
		DeAllocate( HelmholtzEOS % d )
		
		DeAllocate( HelmholtzEOS % f )
		DeAllocate( HelmholtzEOS % fd )
		DeAllocate( HelmholtzEOS % ft )
		DeAllocate( HelmholtzEOS % fdd )
		DeAllocate( HelmholtzEOS % ftt )
		DeAllocate( HelmholtzEOS % fdt )
		DeAllocate( HelmholtzEOS % fddt )
		DeAllocate( HelmholtzEOS % fdtt )
		DeAllocate( HelmholtzEOS % fddtt )
		
		DeAllocate( HelmholtzEOS % dpdf )
		DeAllocate( HelmholtzEOS % dpdfd )
		DeAllocate( HelmholtzEOS % dpdft )
		DeAllocate( HelmholtzEOS % dpdfdt )
		
		DeAllocate( HelmholtzEOS % ef )
		DeAllocate( HelmholtzEOS % efd )
		DeAllocate( HelmholtzEOS % eft )
		DeAllocate( HelmholtzEOS % efdt )
		
		DeAllocate( HelmholtzEOS % xf )
		DeAllocate( HelmholtzEOS % xfd )
		DeAllocate( HelmholtzEOS % xft )
		DeAllocate( HelmholtzEOS % xfdt )
		
		DeAllocate( HelmholtzEOS % dt )
		DeAllocate( HelmholtzEOS % dt2 )
		DeAllocate( HelmholtzEOS % dti )
		DeAllocate( HelmholtzEOS % dt2i )
		
		DeAllocate( HelmholtzEOS % dd )
		DeAllocate( HelmholtzEOS % dd2 )
		DeAllocate( HelmholtzEOS % ddi )
		DeAllocate( HelmholtzEOS % dd2i )	
	
	END SUBROUTINE DeAllocateHelmEOS

	SUBROUTINE AllocateMuonEOS( MuonEOS, nPoints )
		
		TYPE(MuonEOSType)      :: MuonEOS
		INTEGER, DIMENSION(2), INTENT(IN) :: nPoints

		MuonEOS % nPointsTemp  = nPoints(1)
		MuonEOS % nPointsMu   = nPoints(2)

		ALLOCATE( MuonEOS % t( nPoints(1) ) )
		ALLOCATE( MuonEOS % mu( nPoints(2) ) )
		ALLOCATE( MuonEOS % rhoymu( nPoints(1), nPoints(2) ) )
		ALLOCATE( MuonEOS % p( nPoints(1), nPoints(2) ) )
		ALLOCATE( MuonEOS % e( nPoints(1), nPoints(2) ) )
		ALLOCATE( MuonEOS % s( nPoints(1), nPoints(2) ) )

	END SUBROUTINE AllocateMuonEOS

	SUBROUTINE DeAllocateMuonEOS( MuonEOS )
		
		TYPE(MuonEOSType)      :: MuonEOS

		DEALLOCATE( MuonEOS % t )
		DEALLOCATE( MuonEOS % rhoymu )
		DEALLOCATE( MuonEOS % mu )
		DEALLOCATE( MuonEOS % p )
		DEALLOCATE( MuonEOS % e )
		DEALLOCATE( MuonEOS % s )

	END SUBROUTINE DeAllocateMuonEOS

	SUBROUTINE ReadHelmEOSdat(HelmDatFilePath, HelmholtzEOS)
		
		TYPE(HelmholtzEOSType), INTENT(INOUT) :: HelmholtzEOS 
		CHARACTER(len=128), INTENT(IN) :: HelmDatFilePath
		
		! Local variables
		REAL(dp) :: tlo, thi, tstp, tstpi, dlo, dhi, dstp, dstpi, tsav, dsav
		INTEGER :: istat=0
		INTEGER :: iDen, iT
		
        !..   read the helmholtz free energy table
        tlo   = 3.0d0
        thi   = 13.0d0
        tstp  = (thi - tlo)/float(HelmholtzEOS % nPointsTemp-1)
        tstpi = 1.0d0/tstp
        dlo   = -12.0d0
        dhi   = 15.0d0
        dstp  = (dhi - dlo)/float(HelmholtzEOS % nPointsDen-1)
        dstpi = 1.0d0/dstp
		
        do iT=1,HelmholtzEOS % nPointsTemp
			tsav = tlo + (iT-1)*tstp
			HelmholtzEOS % t(iT) = 10.0d0**(tsav)
		enddo
		do iDen=1,HelmholtzEOS % nPointsDen
			dsav = dlo + (iDen-1)*dstp
			HelmholtzEOS % d(iDen) = 10.0d0**(dsav)
		end do
		
		OPEN(UNIT=1234,FILE=TRIM(ADJUSTL(HelmDatFilePath)), STATUS='old', IOSTAT=istat)
		
		IF (istat .ne. 0) THEN
			WRITE(*,*) 'Cannot open HelmholtzEOS % table.dat!'
			STOP
		ENDIF
		
		!..read the helmholtz free energy table
		do iT=1, HelmholtzEOS % nPointsTemp
			do iDen=1, HelmholtzEOS % nPointsDen
				read(1234,*) HelmholtzEOS % f(iDen,iT), HelmholtzEOS % fd(iDen,iT), HelmholtzEOS % ft(iDen,iT),&
				HelmholtzEOS % fdd(iDen,iT), HelmholtzEOS % ftt(iDen,iT), HelmholtzEOS % fdt(iDen,iT), & 
				HelmholtzEOS % fddt(iDen,iT), HelmholtzEOS % fdtt(iDen,iT), HelmholtzEOS % fddtt(iDen,iT)
			enddo
		enddo
		
		!..read the pressure derivative with density table
		DO iT=1, HelmholtzEOS % nPointsTemp
			DO iDen=1, HelmholtzEOS % nPointsDen
				read(1234,*) HelmholtzEOS % dpdf(iDen,iT), HelmholtzEOS % dpdfd(iDen,iT),&
				HelmholtzEOS % dpdft(iDen,iT), HelmholtzEOS % dpdfdt(iDen,iT)
			ENDDO
		ENDDO
		
		!..read the electron chemical potential table
		DO iT=1, HelmholtzEOS % nPointsTemp
			DO iDen=1, HelmholtzEOS % nPointsDen
				READ(1234,*) HelmholtzEOS % ef(iDen,iT), HelmholtzEOS % efd(iDen,iT),&
				HelmholtzEOS % eft(iDen,iT), HelmholtzEOS % efdt(iDen,iT)
			ENDDO
		ENDDO
		
		!..read the number density table
		DO iT=1, HelmholtzEOS % nPointsTemp
			DO iDen=1, HelmholtzEOS % nPointsDen
				READ(1234,*) HelmholtzEOS % xf(iDen,iT), HelmholtzEOS % xfd(iDen,iT),&
				HelmholtzEOS % xft(iDen,iT), HelmholtzEOS % xfdt(iDen,iT)
			ENDDO
		ENDDO
		
		CLOSE(1234)
		
		!..   construct the temperature and density deltas and their inverses
		do iT = 1, HelmholtzEOS % nPointsTemp-1
			HelmholtzEOS % dt(iT)   = HelmholtzEOS % t(iT+1) - HelmholtzEOS % t(iT)
			HelmholtzEOS % dt2(iT)  = HelmholtzEOS % dt(iT) * HelmholtzEOS % dt(iT)
			HelmholtzEOS % dti(iT)  = 1.0d0/HelmholtzEOS % dt(iT)
			HelmholtzEOS % dt2i(iT) = 1.0d0/HelmholtzEOS % dt2(iT)
		end do
		
		do iDen = 1, HelmholtzEOS % nPointsDen-1
			HelmholtzEOS % dd(iDen)   = HelmholtzEOS % d(iDen+1) - HelmholtzEOS % d(iDen)
			HelmholtzEOS % dd2(iDen)  = HelmholtzEOS % dd(iDen) * HelmholtzEOS % dd(iDen)
			HelmholtzEOS % ddi(iDen)  = 1.0d0/HelmholtzEOS % dd(iDen)
			HelmholtzEOS % dd2i(iDen) = 1.0d0/HelmholtzEOS % dd2(iDen)
		end do
		
		! Set up the minimum and maximum possible densities.
		HelmholtzEOS % mintemp = 10.d0**tlo
		HelmholtzEOS % maxtemp = 10.d0**thi
		HelmholtzEOS % mindens = 10.d0**dlo
		HelmholtzEOS % maxdens = 10.d0**dhi
		
	END	SUBROUTINE ReadHelmEOSdat

	SUBROUTINE ReadMuonEOSdat(MuonDatFilePath, MuonEOS)
		
		TYPE(MuonEOSType), INTENT(INOUT) :: MuonEOS 
		CHARACTER(len=128), INTENT(IN) :: MuonDatFilePath
		
		! Local variables
		INTEGER :: istat=0
		INTEGER :: iMu, iT

        !..   read the muon table with Tobias' format
		OPEN(UNIT=1234,FILE=TRIM(ADJUSTL(MuonDatFilePath)), STATUS='old', IOSTAT=istat)
		IF (istat .ne. 0) THEN
			WRITE(*,*) 'Cannot open MuonEOS table.dat!'
			STOP
		ENDIF
		
		DO iT=1,MuonEOS % nPointsTemp
			IF (iT .eq. (nDenMuon+1)*iT) THEN
				READ(1234,*)
			ENDIF
			DO iMu=1,MuonEOS % nPointsMu
				READ(1234,*) MuonEOS % t(iT), MuonEOS % rhoymu(iT,iMu), &
					MuonEOS % p(iT,iMu), MuonEOS % e(iT,iMu), &
					MuonEOS % s(iT,iMu), MuonEOS % mu(iMu)

			ENDDO
		ENDDO
		
		! DO UNIT CONVERSION
		MuonEOS % t =
		
	END	SUBROUTINE ReadMuonEOSdat

END MODULE wlLeptonEOSModule
				
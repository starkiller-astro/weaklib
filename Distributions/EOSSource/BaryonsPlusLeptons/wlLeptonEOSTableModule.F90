MODULE wlLeptonEOSTableModule
  
  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    ergmev, cm3fm3, kmev, kmev_inv, rmu

  IMPLICIT NONE
  PRIVATE
  
  ! nPoints for classic HelmTable
  INTEGER, PARAMETER, PUBLIC :: iTempMax=541, iDenMax=201 

  TYPE, PUBLIC :: HelmTableType
    
    INTEGER :: nPointsDen  !imax
    INTEGER :: nPointsTemp !jmax
    
    REAL(DP), DIMENSION(:), ALLOCATABLE :: t
    REAL(DP), DIMENSION(:), ALLOCATABLE :: d
    
    ! Free Energy Table
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: f
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: fd
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ft
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: fdd
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ftt
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: fdt
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: fddt
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: fdtt
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: fddtt
    
    !  pressure derivative with density table
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dpdf
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dpdfd
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dpdft
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dpdfdt
    
    !  electron chemical potential table
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ef
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: efd
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: eft
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: efdt
    
    !  number density table
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xf
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xfd
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xft
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xfdt
    
    !  deltas needed for interpolation
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dt
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dt2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dti
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dt2i
    
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dd
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dd2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: ddi
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dd2i
    
    ! miniDenm and maxiDenm possible densities
    REAL(DP) :: mintemp
    REAL(DP) :: maxtemp
    REAL(DP) :: mindens
    REAL(DP) :: maxdens

    REAL(DP) :: eos_MinD
    REAL(DP) :: lepton_mass ! To distinguish electrons and muons
    
  END TYPE HelmTableType
  
  PUBLIC AllocateHelmholtzTable
  PUBLIC DeallocateHelmholtzTable
  PUBLIC ReadHelmEOSdat

CONTAINS

  SUBROUTINE AllocateHelmholtzTable( HelmTable, nPoints, eos_MinD )
    
    TYPE(HelmTableType), INTENT(INOUT) :: HelmTable 
    INTEGER, DIMENSION(2), INTENT(IN)  :: nPoints
    REAL(DP), INTENT(in), OPTIONAL     :: eos_MinD

    HelmTable % nPointsDen   = nPoints(1)
    HelmTable % nPointsTemp  = nPoints(2)
    
    IF ( PRESENT (eos_MinD) ) THEN
      HelmTable % eos_MinD = eos_MinD
    ELSE
      HelmTable % eos_MinD = 0.0_dp
    END IF

    ALLOCATE( HelmTable % d( nPoints(1) ) )
    ALLOCATE( HelmTable % t( nPoints(2) ) )
    
    ALLOCATE( HelmTable % f    ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % fd   ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % ft   ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % fdd  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % ftt  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % fdt  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % fddt ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % fdtt ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % fddtt( nPoints(1), nPoints(2) ) )
    
    ALLOCATE( HelmTable % dpdf  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % dpdfd ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % dpdft ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % dpdfdt( nPoints(1), nPoints(2) ) )
    
    ALLOCATE( HelmTable % ef  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % efd ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % eft ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % efdt( nPoints(1), nPoints(2) ) )
    
    ALLOCATE( HelmTable % xf  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % xfd ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % xft ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmTable % xfdt( nPoints(1), nPoints(2) ) )
    
    ALLOCATE( HelmTable % dt  ( nPoints(2) ) )
    ALLOCATE( HelmTable % dt2 ( nPoints(2) ) )
    ALLOCATE( HelmTable % dti ( nPoints(2) ) )
    ALLOCATE( HelmTable % dt2i( nPoints(2) ) )
    
    ALLOCATE( HelmTable % dd  ( nPoints(1) ) )
    ALLOCATE( HelmTable % dd2 ( nPoints(1) ) )
    ALLOCATE( HelmTable % ddi ( nPoints(1) ) )
    ALLOCATE( HelmTable % dd2i( nPoints(1) ) )
    
  END SUBROUTINE AllocateHelmholtzTable
  
  SUBROUTINE DeallocateHelmholtzTable( HelmTable )
  
    TYPE(HelmTableType)      :: HelmTable 

    DEALLOCATE( HelmTable % t )
    DEALLOCATE( HelmTable % d )
    
    DEALLOCATE( HelmTable % f )
    DEALLOCATE( HelmTable % fd )
    DEALLOCATE( HelmTable % ft )
    DEALLOCATE( HelmTable % fdd )
    DEALLOCATE( HelmTable % ftt )
    DEALLOCATE( HelmTable % fdt )
    DEALLOCATE( HelmTable % fddt )
    DEALLOCATE( HelmTable % fdtt )
    DEALLOCATE( HelmTable % fddtt )
    
    DEALLOCATE( HelmTable % dpdf )
    DEALLOCATE( HelmTable % dpdfd )
    DEALLOCATE( HelmTable % dpdft )
    DEALLOCATE( HelmTable % dpdfdt )
    
    DEALLOCATE( HelmTable % ef )
    DEALLOCATE( HelmTable % efd )
    DEALLOCATE( HelmTable % eft )
    DEALLOCATE( HelmTable % efdt )
    
    DEALLOCATE( HelmTable % xf )
    DEALLOCATE( HelmTable % xfd )
    DEALLOCATE( HelmTable % xft )
    DEALLOCATE( HelmTable % xfdt )
    
    DEALLOCATE( HelmTable % dt )
    DEALLOCATE( HelmTable % dt2 )
    DEALLOCATE( HelmTable % dti )
    DEALLOCATE( HelmTable % dt2i )
    
    DEALLOCATE( HelmTable % dd )
    DEALLOCATE( HelmTable % dd2 )
    DEALLOCATE( HelmTable % ddi )
    DEALLOCATE( HelmTable % dd2i )  

  END SUBROUTINE DeallocateHelmholtzTable

  SUBROUTINE ReadHelmEOSdat(HelmDatFilePath, HelmTable, lepton_mass, &
      t_low, t_high, d_low, d_high)
    
    CHARACTER(len=128) , INTENT(IN)    :: HelmDatFilePath
    TYPE(HelmTableType), INTENT(INOUT) :: HelmTable 
    REAL(DP), INTENT(IN) :: lepton_mass
    
    REAL(DP), INTENT(IN), OPTIONAL :: t_low, t_high, d_low, d_high
    ! Local variables
    REAL(DP) :: tlo, thi, dlo, dhi
    REAL(DP) :: tstp, tstpi, dstp, dstpi, tsav, dsav
    INTEGER  :: istat = 0
    INTEGER  :: iDen, iT
    
    HelmTable % lepton_mass = lepton_mass

    IF (PRESENT(t_low)) THEN
        tlo = t_low
    ELSE
        tlo = 3.0_dp ! Default value for tlo of standard Helmholtz table for electrons
    END IF
    
    IF (PRESENT(t_high)) THEN
        thi = t_high
    ELSE
        thi = 13.0_dp ! Default value for thi of standard Helmholtz table for electrons
    END IF
    
    IF (PRESENT(d_low)) THEN
        dlo = d_low
    ELSE
        dlo = -12.0_dp ! Default value for dlo of standard Helmholtz table for electrons
    END IF
    
    IF (PRESENT(d_high)) THEN
        dhi = d_high
    ELSE
        dhi = 15.0_dp ! Default value for dhi of standard Helmholtz table for electrons
    END IF

    !..   read the helmholtz free energy table
    tstp  = (thi - tlo) / REAL(HelmTable % nPointsTemp-1, DP)
    tstpi = 1.0_dp/tstp

    dstp  = (dhi - dlo) / REAL(HelmTable % nPointsDen -1, DP)
    dstpi = 1.0_dp/dstp

    DO iT=1,HelmTable % nPointsTemp
      tsav = tlo + REAL(iT-1, DP)*tstp
      HelmTable % t(iT) = 10.0_dp**(tsav)
    END DO
    DO iDen=1,HelmTable % nPointsDen
      dsav = dlo + REAL(iDen-1, DP)*dstp
      HelmTable % d(iDen) = 10.0_dp**(dsav)
    END DO
    
    OPEN(UNIT=1234,FILE=TRIM(ADJUSTL(HelmDatFilePath)), STATUS='old', IOSTAT=istat)
    
    IF (istat .ne. 0) THEN
      WRITE(*,*) 'Cannot open HelmTable % table.dat!'
      STOP
    ENDIF
    
    !..read the helmholtz free energy table
    DO iT=1, HelmTable % nPointsTemp
      DO iDen=1, HelmTable % nPointsDen
        read(1234,*) HelmTable % f(iDen,iT), HelmTable % fd(iDen,iT), HelmTable % ft(iDen,iT),&
        HelmTable % fdd(iDen,iT) , HelmTable % ftt(iDen,iT) , HelmTable % fdt(iDen,iT), & 
        HelmTable % fddt(iDen,iT), HelmTable % fdtt(iDen,iT), HelmTable % fddtt(iDen,iT)
      END DO
    END DO
    
    !..read the pressure derivative with density table
    DO iT=1, HelmTable % nPointsTemp
      DO iDen=1, HelmTable % nPointsDen
        read(1234,*) HelmTable % dpdf(iDen,iT), HelmTable % dpdfd(iDen,iT),&
        HelmTable % dpdft(iDen,iT), HelmTable % dpdfdt(iDen,iT)
      ENDDO
    ENDDO
    
    !..read the electron chemical potential table
    DO iT=1, HelmTable % nPointsTemp
      DO iDen=1, HelmTable % nPointsDen
        READ(1234,*) HelmTable % ef(iDen,iT), HelmTable % efd(iDen,iT),&
        HelmTable % eft(iDen,iT), HelmTable % efdt(iDen,iT)
      ENDDO
    ENDDO
    
    !..read the number density table
    DO iT=1, HelmTable % nPointsTemp
      DO iDen=1, HelmTable % nPointsDen
        READ(1234,*) HelmTable % xf(iDen,iT), HelmTable % xfd(iDen,iT),&
        HelmTable % xft(iDen,iT), HelmTable % xfdt(iDen,iT)
      ENDDO
    ENDDO
    
    CLOSE(1234)
    
    !..   construct the temperature and density deltas and their inverses
    DO iT = 1, HelmTable % nPointsTemp-1
      HelmTable % dt(iT)   = HelmTable % t(iT+1) - HelmTable % t(iT)
      HelmTable % dt2(iT)  = HelmTable % dt(iT) * HelmTable % dt(iT)
      HelmTable % dti(iT)  = 1.0_dp/HelmTable % dt(iT)
      HelmTable % dt2i(iT) = 1.0_dp/HelmTable % dt2(iT)
    END DO
    
    DO iDen = 1, HelmTable % nPointsDen-1
      HelmTable % dd(iDen)   = HelmTable % d(iDen+1) - HelmTable % d(iDen)
      HelmTable % dd2(iDen)  = HelmTable % dd(iDen) * HelmTable % dd(iDen)
      HelmTable % ddi(iDen)  = 1.0_dp/HelmTable % dd(iDen)
      HelmTable % dd2i(iDen) = 1.0_dp/HelmTable % dd2(iDen)
    END DO
    
    ! Set up the miniDenm and maxiDenm possible densities.
    HelmTable % mintemp = 10.0_dp**tlo
    HelmTable % maxtemp = 10.0_dp**thi
    HelmTable % mindens = 10.0_dp**dlo
    HelmTable % maxdens = 10.0_dp**dhi
    
  END  SUBROUTINE ReadHelmEOSdat

END MODULE wlLeptonEOSTableModule
        

MODULE wlLeptonEOSModule
  
  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    ergmev, cm3fm3, kmev, kmev_inv, rmu

  IMPLICIT NONE
  PRIVATE
  
    INTEGER, PARAMETER, PUBLIC :: iTempMax=541, iDenMax=201 
    INTEGER, PARAMETER, PUBLIC :: nTempMuon=270, nDenMuon=1501 
  
   TYPE, PUBLIC :: HelmholtzTableType
    
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
    
    ! miniDenm and maxiDenm possible densities
    REAL(dp) :: mintemp
    REAL(dp) :: maxtemp
    REAL(dp) :: mindens
    REAL(dp) :: maxdens
    
  END TYPE HelmholtzTableType
  
  TYPE, PUBLIC :: MuonEOSType
      
    INTEGER :: nPointsDen
    INTEGER :: nPointsTemp
    
    REAL(dp), DIMENSION(:),   ALLOCATABLE :: t
    REAL(dp), DIMENSION(:),   ALLOCATABLE :: rhoym
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: mu
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: e
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: s
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnPdlnrho
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnPdlnT
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnsdlnrho
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnsdlnT
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnedlnrho
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnedlnT

  END TYPE MuonEOSType

  PUBLIC AllocateHelmholtzTable
  PUBLIC DeallocateHelmholtzTable
  PUBLIC ReadHelmEOSdat

  PUBLIC AllocateMuonEOS
  PUBLIC DeAllocateMuonEOS
  PUBLIC ReadMuonEOSdat
  
CONTAINS

  SUBROUTINE AllocateHelmholtzTable( HelmholtzTable, nPoints )
    
    TYPE(HelmholtzTableType), INTENT(INOUT) :: HelmholtzTable 
    INTEGER, DIMENSION(2), INTENT(IN)     :: nPoints
        
    HelmholtzTable % nPointsDen   = nPoints(1)
    HelmholtzTable % nPointsTemp  = nPoints(2)
    
    ALLOCATE( HelmholtzTable % d( nPoints(1) ) )
    ALLOCATE( HelmholtzTable % t( nPoints(2) ) )
    
    ALLOCATE( HelmholtzTable % f    ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % fd   ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % ft   ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % fdd  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % ftt  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % fdt  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % fddt ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % fdtt ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % fddtt( nPoints(1), nPoints(2) ) )
    
    ALLOCATE( HelmholtzTable % dpdf  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % dpdfd ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % dpdft ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % dpdfdt( nPoints(1), nPoints(2) ) )
    
    ALLOCATE( HelmholtzTable % ef  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % efd ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % eft ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % efdt( nPoints(1), nPoints(2) ) )
    
    ALLOCATE( HelmholtzTable % xf  ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % xfd ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % xft ( nPoints(1), nPoints(2) ) )
    ALLOCATE( HelmholtzTable % xfdt( nPoints(1), nPoints(2) ) )
    
    ALLOCATE( HelmholtzTable % dt  ( nPoints(2) ) )
    ALLOCATE( HelmholtzTable % dt2 ( nPoints(2) ) )
    ALLOCATE( HelmholtzTable % dti ( nPoints(2) ) )
    ALLOCATE( HelmholtzTable % dt2i( nPoints(2) ) )
    
    ALLOCATE( HelmholtzTable % dd  ( nPoints(1) ) )
    ALLOCATE( HelmholtzTable % dd2 ( nPoints(1) ) )
    ALLOCATE( HelmholtzTable % ddi ( nPoints(1) ) )
    ALLOCATE( HelmholtzTable % dd2i( nPoints(1) ) )
    
  END SUBROUTINE AllocateHelmholtzTable
  
  SUBROUTINE DeallocateHelmholtzTable( HelmholtzTable )
  
    TYPE(HelmholtzTableType)      :: HelmholtzTable 

    DeAllocate( HelmholtzTable % t )
    DeAllocate( HelmholtzTable % d )
    
    DeAllocate( HelmholtzTable % f )
    DeAllocate( HelmholtzTable % fd )
    DeAllocate( HelmholtzTable % ft )
    DeAllocate( HelmholtzTable % fdd )
    DeAllocate( HelmholtzTable % ftt )
    DeAllocate( HelmholtzTable % fdt )
    DeAllocate( HelmholtzTable % fddt )
    DeAllocate( HelmholtzTable % fdtt )
    DeAllocate( HelmholtzTable % fddtt )
    
    DeAllocate( HelmholtzTable % dpdf )
    DeAllocate( HelmholtzTable % dpdfd )
    DeAllocate( HelmholtzTable % dpdft )
    DeAllocate( HelmholtzTable % dpdfdt )
    
    DeAllocate( HelmholtzTable % ef )
    DeAllocate( HelmholtzTable % efd )
    DeAllocate( HelmholtzTable % eft )
    DeAllocate( HelmholtzTable % efdt )
    
    DeAllocate( HelmholtzTable % xf )
    DeAllocate( HelmholtzTable % xfd )
    DeAllocate( HelmholtzTable % xft )
    DeAllocate( HelmholtzTable % xfdt )
    
    DeAllocate( HelmholtzTable % dt )
    DeAllocate( HelmholtzTable % dt2 )
    DeAllocate( HelmholtzTable % dti )
    DeAllocate( HelmholtzTable % dt2i )
    
    DeAllocate( HelmholtzTable % dd )
    DeAllocate( HelmholtzTable % dd2 )
    DeAllocate( HelmholtzTable % ddi )
    DeAllocate( HelmholtzTable % dd2i )  
  
  END SUBROUTINE DeallocateHelmholtzTable

  SUBROUTINE AllocateMuonEOS( MuonEOS, nPoints )
    
    TYPE(MuonEOSType)      :: MuonEOS
    INTEGER, DIMENSION(2), INTENT(IN) :: nPoints

    MuonEOS % nPointsTemp  = nPoints(1)
    MuonEOS % nPointsDen   = nPoints(2)

    ALLOCATE( MuonEOS % t( nPoints(1) ) )
    ALLOCATE( MuonEOS % rhoym( nPoints(2) ) )
    ALLOCATE( MuonEOS % mu( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % p( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % e( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % s( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % dlnPdlnrho( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % dlnPdlnT( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % dlnsdlnrho( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % dlnsdlnT( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % dlnedlnrho( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonEOS % dlnedlnT( nPoints(1), nPoints(2) ) )

  END SUBROUTINE AllocateMuonEOS

  SUBROUTINE DeAllocateMuonEOS( MuonEOS )
    
    TYPE(MuonEOSType)      :: MuonEOS

    DEALLOCATE( MuonEOS % t )
    DEALLOCATE( MuonEOS % rhoym )
    DEALLOCATE( MuonEOS % mu )
    DEALLOCATE( MuonEOS % p )
    DEALLOCATE( MuonEOS % e )
    DEALLOCATE( MuonEOS % s )
    DEALLOCATE( MuonEOS % dlnPdlnrho )
    DEALLOCATE( MuonEOS % dlnPdlnT )
    DEALLOCATE( MuonEOS % dlnsdlnrho )
    DEALLOCATE( MuonEOS % dlnsdlnT )
    DEALLOCATE( MuonEOS % dlnedlnrho )
    DEALLOCATE( MuonEOS % dlnedlnT )

  END SUBROUTINE DeAllocateMuonEOS

  SUBROUTINE ReadHelmEOSdat(HelmDatFilePath, HelmholtzTable)
    
    TYPE(HelmholtzTableType), INTENT(INOUT) :: HelmholtzTable 
    CHARACTER(len=128), INTENT(IN) :: HelmDatFilePath
    
    ! Local variables
    REAL(dp) :: tlo, thi, tstp, tstpi, dlo, dhi, dstp, dstpi, tsav, dsav
    INTEGER :: istat=0
    INTEGER :: iDen, iT
    
        !..   read the helmholtz free energy table
        tlo   = 3.0d0
        thi   = 13.0d0
        tstp  = (thi - tlo)/float(HelmholtzTable % nPointsTemp-1)
        tstpi = 1.0d0/tstp
        dlo   = -12.0d0
        dhi   = 15.0d0
        dstp  = (dhi - dlo)/float(HelmholtzTable % nPointsDen-1)
        dstpi = 1.0d0/dstp
    
        do iT=1,HelmholtzTable % nPointsTemp
      tsav = tlo + (iT-1)*tstp
      HelmholtzTable % t(iT) = 10.0d0**(tsav)
    enddo
    do iDen=1,HelmholtzTable % nPointsDen
      dsav = dlo + (iDen-1)*dstp
      HelmholtzTable % d(iDen) = 10.0d0**(dsav)
    end do
    
    OPEN(UNIT=1234,FILE=TRIM(ADJUSTL(HelmDatFilePath)), STATUS='old', IOSTAT=istat)
    
    IF (istat .ne. 0) THEN
      WRITE(*,*) 'Cannot open HelmholtzTable % table.dat!'
      STOP
    ENDIF
    
    !..read the helmholtz free energy table
    do iT=1, HelmholtzTable % nPointsTemp
      do iDen=1, HelmholtzTable % nPointsDen
        read(1234,*) HelmholtzTable % f(iDen,iT), HelmholtzTable % fd(iDen,iT), HelmholtzTable % ft(iDen,iT),&
        HelmholtzTable % fdd(iDen,iT), HelmholtzTable % ftt(iDen,iT), HelmholtzTable % fdt(iDen,iT), & 
        HelmholtzTable % fddt(iDen,iT), HelmholtzTable % fdtt(iDen,iT), HelmholtzTable % fddtt(iDen,iT)
      enddo
    enddo
    
    !..read the pressure derivative with density table
    DO iT=1, HelmholtzTable % nPointsTemp
      DO iDen=1, HelmholtzTable % nPointsDen
        read(1234,*) HelmholtzTable % dpdf(iDen,iT), HelmholtzTable % dpdfd(iDen,iT),&
        HelmholtzTable % dpdft(iDen,iT), HelmholtzTable % dpdfdt(iDen,iT)
      ENDDO
    ENDDO
    
    !..read the electron chemical potential table
    DO iT=1, HelmholtzTable % nPointsTemp
      DO iDen=1, HelmholtzTable % nPointsDen
        READ(1234,*) HelmholtzTable % ef(iDen,iT), HelmholtzTable % efd(iDen,iT),&
        HelmholtzTable % eft(iDen,iT), HelmholtzTable % efdt(iDen,iT)
      ENDDO
    ENDDO
    
    !..read the number density table
    DO iT=1, HelmholtzTable % nPointsTemp
      DO iDen=1, HelmholtzTable % nPointsDen
        READ(1234,*) HelmholtzTable % xf(iDen,iT), HelmholtzTable % xfd(iDen,iT),&
        HelmholtzTable % xft(iDen,iT), HelmholtzTable % xfdt(iDen,iT)
      ENDDO
    ENDDO
    
    CLOSE(1234)
    
    !..   construct the temperature and density deltas and their inverses
    do iT = 1, HelmholtzTable % nPointsTemp-1
      HelmholtzTable % dt(iT)   = HelmholtzTable % t(iT+1) - HelmholtzTable % t(iT)
      HelmholtzTable % dt2(iT)  = HelmholtzTable % dt(iT) * HelmholtzTable % dt(iT)
      HelmholtzTable % dti(iT)  = 1.0d0/HelmholtzTable % dt(iT)
      HelmholtzTable % dt2i(iT) = 1.0d0/HelmholtzTable % dt2(iT)
    end do
    
    do iDen = 1, HelmholtzTable % nPointsDen-1
      HelmholtzTable % dd(iDen)   = HelmholtzTable % d(iDen+1) - HelmholtzTable % d(iDen)
      HelmholtzTable % dd2(iDen)  = HelmholtzTable % dd(iDen) * HelmholtzTable % dd(iDen)
      HelmholtzTable % ddi(iDen)  = 1.0d0/HelmholtzTable % dd(iDen)
      HelmholtzTable % dd2i(iDen) = 1.0d0/HelmholtzTable % dd2(iDen)
    end do
    
    ! Set up the miniDenm and maxiDenm possible densities.
    HelmholtzTable % mintemp = 10.d0**tlo
    HelmholtzTable % maxtemp = 10.d0**thi
    HelmholtzTable % mindens = 10.d0**dlo
    HelmholtzTable % maxdens = 10.d0**dhi
    
  END  SUBROUTINE ReadHelmEOSdat

  SUBROUTINE ReadMuonEOSdat(MuonDatFilePath, MuonEOS)
    
    TYPE(MuonEOSType), INTENT(INOUT) :: MuonEOS 
    CHARACTER(len=128), INTENT(IN) :: MuonDatFilePath
    
    ! Local variables
    INTEGER :: istat=0
    INTEGER :: iDen, iT

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
      DO iDen=1,MuonEOS % nPointsDen
        READ(1234,*) MuonEOS % t(iT), MuonEOS % rhoym(iDen), &
          MuonEOS % p(iT,iDen), MuonEOS % e(iT,iDen), &
          MuonEOS % s(iT,iDen), MuonEOS % mu(iT,iDen)
      ENDDO
    ENDDO
    
    ! DO UNIT CONVERSION
    ! Based on Tobias's code, the chemical potential includes the muon rest mass
    MuonEOS % t = MuonEOS % t * kmev_inv !MeV to kelvin
    MuonEOS % rhoym = MuonEOS % rhoym * rmu/cm3fm3 !baryon number to g/cm^3
    MuonEOS % p = MuonEOS % p * ergmev / cm3fm3 !MeV/fm^3 to dyn
    MuonEOS % e = MuonEOS % e * ergmev / rmu 
    MuonEOS % s = MuonEOS % s / (kmev * ergmev / rmu)
    
    ! NOW CALCULATE DERIVATIVES AT CONSTANT RHO
    DO iDen=1,MuonEOS % nPointsDen
      MuonEOS % dlnsdlnT(:,iDen) = Gradient3pts( LOG10(MuonEOS % s(:,iDen)), LOG10(MuonEOS % t), MuonEOS % nPointsTemp)
      MuonEOS % dlnPdlnT(:,iDen) = Gradient3pts( LOG10(MuonEOS % p(:,iDen)), LOG10(MuonEOS % t), MuonEOS % nPointsTemp)
      MuonEOS % dlnedlnT(:,iDen) = Gradient3pts( LOG10(MuonEOS % e(:,iDen)), LOG10(MuonEOS % t), MuonEOS % nPointsTemp)
    ENDDO

    ! NOW CALCULATE DERIVATIVES AT CONSTANT T
    DO iT=1,MuonEOS % nPointsTemp
      MuonEOS % dlnsdlnrho(iT,:) = Gradient3pts( LOG10(MuonEOS % s(iT,:)), LOG10(MuonEOS % rhoym), MuonEOS % nPointsDen)
      MuonEOS % dlnPdlnrho(iT,:) = Gradient3pts( LOG10(MuonEOS % p(iT,:)), LOG10(MuonEOS % rhoym), MuonEOS % nPointsDen)
      MuonEOS % dlnedlnrho(iT,:) = Gradient3pts( LOG10(MuonEOS % e(iT,:)), LOG10(MuonEOS % rhoym), MuonEOS % nPointsDen)
    ENDDO
    
    WRITE(*,*) 'Bounds of the Muon Table'
    WRITE(*,*) MAXVAL(MuonEOS % t), MINVAL(MuonEOS % t)
    WRITE(*,*) MAXVAL(MuonEOS % rhoym), MINVAL(MuonEOS % rhoym)
    WRITE(*,*) MAXVAL(MuonEOS % p), MINVAL(MuonEOS % p)
    WRITE(*,*) MAXVAL(MuonEOS % e), MINVAL(MuonEOS % e)
    WRITE(*,*) MAXVAL(MuonEOS % s), MINVAL(MuonEOS % s)
    WRITE(*,*) MAXVAL(MuonEOS % mu), MINVAL(MuonEOS % mu)
    
  END SUBROUTINE ReadMuonEOSdat

  ! Tools to calculate derivatives, maybe worth moving into another MODULE
  FUNCTION Gradient(y, x, nX) RESULT(grad)
   
    INTEGER,  INTENT(IN) :: nX
    REAL(DP), INTENT(IN) :: y(nX), x(nX)
    
    INTEGER  :: i
    REAL(DP) :: grad(nX)

    grad(:) = 0.0d0
    !calculate derivativees with standard gradient formula
    DO i=2,nX-1
      grad(i) = (y(i+1) - y(i-1))/(x(i+1) - x(i-1))
    ENDDO
    
    ! calculate one-sided derivatives
    grad(1) = (y(2) - y(1))/ (x(2) - x(1))
    grad(nX) = (y(nX) - y(nX-1))/ (x(nX) - x(nX-1))

  END FUNCTION Gradient 

  FUNCTION Gradient3pts(y, x, nX) RESULT(grad)

    INTEGER,  INTENT(IN) :: nX
    REAL(DP), INTENT(IN) :: y(nX), x(nX)
    
    INTEGER  :: i
    REAL(DP) :: grad(nX)
    REAL(DP) :: h1, h2

    grad(:) = 0.0d0
    
    ! calculate derivatives with the 3 point formula
    DO i=2,nX-1
     
     h1 = x(i)   - x(i-1)
     h2 = x(i+1) - x(i)
    
     grad(i) = - y(i-1)*(h2/(h1*(h1+h2))) &
           + y(i  )*((h2-h1)/(h1*h2)) &
           + y(i+1)*((h1/(h2*(h1+h2))))
    
    ENDDO
            
    ! calculate one-sided derivatives
    grad(1) = (y(2) - y(1))/ (x(2) - x(1))
    grad(nX) = (y(nX) - y(nX-1))/ (x(nX) - x(nX-1))

  END FUNCTION Gradient3pts

  FUNCTION Gradient5pts(y, x, nX) RESULT(grad)

    INTEGER,  INTENT(IN) :: nX
    REAL(DP), INTENT(IN) :: y(nX), x(nX)
    
    INTEGER  :: i
    REAL(DP) :: grad(nX)
    REAL(DP) :: h1,h2,h3,h4,Ht1,Ht2

    grad(:) = 0.0d0
    
    ! calculate derivatives with the 5 point formula
    DO i=3,nX-2
     
     h1 = x(i-1) - x(i-2)
     h2 = x(i)   - x(i-1)
     h3 = x(i+1) - x(i)
     h4 = x(i+2) - x(i+1)
     Ht1 = h1+h2+h3
     Ht2 = Ht1+h4
    
     grad(i) = y(i-2)*(h2*h3*(h3+h4))/(h1*(h1+h2)*Ht1*Ht2) &
         - y(i-1)*((h1+h2)*h3*(h3+h4))/(h1*h2*(h2+h3)*(Ht2-h1)) &
         + y(i)*((h1+2.0d0*h2)*h3*(h3+h4)-(h1+h2)*h2*(2.0d0*h3+h4)) &
             /((h1+h2)*h2*h3*(h3+h4)) &
         + y(i+1)*((h2+h1)*h2*(h3+h4))/(Ht1*(h2+h3)*h3*h4) &
         - y(i+2)*(h2*(h1+h2)*h3)/(Ht2*(Ht2-h1)*(h3+h4)*h4)

    ENDDO
    
    ! simple gradient should be okay
    grad(2) = (y(3) - y(1)) / (x(3) - x(1))
    grad(nX-1) = (y(nX) - y(nX-2)) / (x(nX) - x(nX-2))
             
    ! calculate one-sided derivatives
    grad(1) = (y(2) - y(1)) / (x(2) - x(1))
    grad(nX) = (y(nX) - y(nX-1)) / (x(nX) - x(nX-1))

  END FUNCTION Gradient5pts

END MODULE wlLeptonEOSModule
        
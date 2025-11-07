MODULE wlLeptonEOSTableModule
  
  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    ergmev, cm3fm3, kmev, kmev_inv, rmu

  IMPLICIT NONE
  PRIVATE
  
    INTEGER, PARAMETER, PUBLIC :: iTempMax=541, iDenMax=201 
    INTEGER, PARAMETER, PUBLIC :: nTempMuon=270, nDenMuon=1301 
  
   TYPE, PUBLIC :: HelmTableType
    
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
    
  END TYPE HelmTableType
  
  TYPE, PUBLIC :: MuonTableType
      
    INTEGER :: nPointsDen
    INTEGER :: nPointsTemp
    
    REAL(dp), DIMENSION(:),   ALLOCATABLE :: t
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: mu
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rhoym
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: e
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: s
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnPdlnrho
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnPdlnT
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnsdlnrho
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnsdlnT
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnedlnrho
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dlnedlnT

    REAL(DP) :: eos_MinD

  END TYPE MuonTableType

  PUBLIC AllocateHelmholtzTable
  PUBLIC DeallocateHelmholtzTable
  PUBLIC ReadHelmEOSdat

  PUBLIC AllocateMuonTable
  PUBLIC DeAllocateMuonTable
  PUBLIC ReadMuonEOSdat
  
CONTAINS

  SUBROUTINE AllocateHelmholtzTable( HelmTable, nPoints )
    
    TYPE(HelmTableType), INTENT(INOUT) :: HelmTable 
    INTEGER, DIMENSION(2), INTENT(IN)     :: nPoints
        
    HelmTable % nPointsDen   = nPoints(1)
    HelmTable % nPointsTemp  = nPoints(2)
    
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

    DeAllocate( HelmTable % t )
    DeAllocate( HelmTable % d )
    
    DeAllocate( HelmTable % f )
    DeAllocate( HelmTable % fd )
    DeAllocate( HelmTable % ft )
    DeAllocate( HelmTable % fdd )
    DeAllocate( HelmTable % ftt )
    DeAllocate( HelmTable % fdt )
    DeAllocate( HelmTable % fddt )
    DeAllocate( HelmTable % fdtt )
    DeAllocate( HelmTable % fddtt )
    
    DeAllocate( HelmTable % dpdf )
    DeAllocate( HelmTable % dpdfd )
    DeAllocate( HelmTable % dpdft )
    DeAllocate( HelmTable % dpdfdt )
    
    DeAllocate( HelmTable % ef )
    DeAllocate( HelmTable % efd )
    DeAllocate( HelmTable % eft )
    DeAllocate( HelmTable % efdt )
    
    DeAllocate( HelmTable % xf )
    DeAllocate( HelmTable % xfd )
    DeAllocate( HelmTable % xft )
    DeAllocate( HelmTable % xfdt )
    
    DeAllocate( HelmTable % dt )
    DeAllocate( HelmTable % dt2 )
    DeAllocate( HelmTable % dti )
    DeAllocate( HelmTable % dt2i )
    
    DeAllocate( HelmTable % dd )
    DeAllocate( HelmTable % dd2 )
    DeAllocate( HelmTable % ddi )
    DeAllocate( HelmTable % dd2i )  
  
  END SUBROUTINE DeallocateHelmholtzTable

  SUBROUTINE AllocateMuonTable( MuonTable, nPoints, eos_MinD )
    
    TYPE(MuonTableType)      :: MuonTable
    INTEGER, DIMENSION(2), INTENT(IN) :: nPoints
    REAL(DP), INTENT(in), OPTIONAL :: eos_MinD


    IF ( PRESENT (eos_MinD) ) THEN
      MuonTable % eos_MinD = eos_MinD
    ELSE
      MuonTable % eos_MinD = 0.0d0
    END IF

    MuonTable % nPointsTemp  = nPoints(1)
    MuonTable % nPointsDen    = nPoints(2)

    ALLOCATE( MuonTable % t( nPoints(1) ) )
    ALLOCATE( MuonTable % mu( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % rhoym( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % p( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % e( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % s( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % dlnPdlnrho( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % dlnPdlnT( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % dlnsdlnrho( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % dlnsdlnT( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % dlnedlnrho( nPoints(1), nPoints(2) ) )
    ALLOCATE( MuonTable % dlnedlnT( nPoints(1), nPoints(2) ) )

  END SUBROUTINE AllocateMuonTable

  SUBROUTINE DeAllocateMuonTable( MuonTable )
    
    TYPE(MuonTableType)      :: MuonTable

    DEALLOCATE( MuonTable % t )
    DEALLOCATE( MuonTable % rhoym )
    DEALLOCATE( MuonTable % mu )
    DEALLOCATE( MuonTable % p )
    DEALLOCATE( MuonTable % e )
    DEALLOCATE( MuonTable % s )
    DEALLOCATE( MuonTable % dlnPdlnrho )
    DEALLOCATE( MuonTable % dlnPdlnT )
    DEALLOCATE( MuonTable % dlnsdlnrho )
    DEALLOCATE( MuonTable % dlnsdlnT )
    DEALLOCATE( MuonTable % dlnedlnrho )
    DEALLOCATE( MuonTable % dlnedlnT )

  END SUBROUTINE DeAllocateMuonTable

  SUBROUTINE ReadHelmEOSdat(HelmDatFilePath, HelmTable)
    
    TYPE(HelmTableType), INTENT(INOUT) :: HelmTable 
    CHARACTER(len=128), INTENT(IN) :: HelmDatFilePath
    
    ! Local variables
    REAL(dp) :: tlo, thi, tstp, tstpi, dlo, dhi, dstp, dstpi, tsav, dsav
    INTEGER :: istat=0
    INTEGER :: iDen, iT
    
        !..   read the helmholtz free energy table
        tlo   = 3.0d0
        thi   = 13.0d0
        tstp  = (thi - tlo)/real(HelmTable % nPointsTemp-1, dp)
        tstpi = 1.0d0/tstp
        dlo   = -12.0d0
        dhi   = 15.0d0
        dstp  = (dhi - dlo)/real(HelmTable % nPointsDen-1,  dp)
        dstpi = 1.0d0/dstp
    
        do iT=1,HelmTable % nPointsTemp
      tsav = tlo + (iT-1)*tstp
      HelmTable % t(iT) = 10.0d0**(tsav)
    enddo
    do iDen=1,HelmTable % nPointsDen
      dsav = dlo + (iDen-1)*dstp
      HelmTable % d(iDen) = 10.0d0**(dsav)
    end do
    
    OPEN(UNIT=1234,FILE=TRIM(ADJUSTL(HelmDatFilePath)), STATUS='old', IOSTAT=istat)
    
    IF (istat .ne. 0) THEN
      WRITE(*,*) 'Cannot open HelmTable % table.dat!'
      STOP
    ENDIF
    
    !..read the helmholtz free energy table
    do iT=1, HelmTable % nPointsTemp
      do iDen=1, HelmTable % nPointsDen
        read(1234,*) HelmTable % f(iDen,iT), HelmTable % fd(iDen,iT), HelmTable % ft(iDen,iT),&
        HelmTable % fdd(iDen,iT), HelmTable % ftt(iDen,iT), HelmTable % fdt(iDen,iT), & 
        HelmTable % fddt(iDen,iT), HelmTable % fdtt(iDen,iT), HelmTable % fddtt(iDen,iT)
      enddo
    enddo
    
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
    do iT = 1, HelmTable % nPointsTemp-1
      HelmTable % dt(iT)   = HelmTable % t(iT+1) - HelmTable % t(iT)
      HelmTable % dt2(iT)  = HelmTable % dt(iT) * HelmTable % dt(iT)
      HelmTable % dti(iT)  = 1.0d0/HelmTable % dt(iT)
      HelmTable % dt2i(iT) = 1.0d0/HelmTable % dt2(iT)
    end do
    
    do iDen = 1, HelmTable % nPointsDen-1
      HelmTable % dd(iDen)   = HelmTable % d(iDen+1) - HelmTable % d(iDen)
      HelmTable % dd2(iDen)  = HelmTable % dd(iDen) * HelmTable % dd(iDen)
      HelmTable % ddi(iDen)  = 1.0d0/HelmTable % dd(iDen)
      HelmTable % dd2i(iDen) = 1.0d0/HelmTable % dd2(iDen)
    end do
    
    ! Set up the miniDenm and maxiDenm possible densities.
    HelmTable % mintemp = 10.d0**tlo
    HelmTable % maxtemp = 10.d0**thi
    HelmTable % mindens = 10.d0**dlo
    HelmTable % maxdens = 10.d0**dhi
    
  END  SUBROUTINE ReadHelmEOSdat

  SUBROUTINE ReadMuonEOSdat(MuonDatFilePath, MuonTable)
    
    TYPE(MuonTableType), INTENT(INOUT) :: MuonTable 
    CHARACTER(len=128), INTENT(IN) :: MuonDatFilePath
    
    ! Local variables
    INTEGER :: istat=0
    INTEGER :: iDen, iT

        !..   read the muon table with Tobias' format
    OPEN(UNIT=1234,FILE=TRIM(ADJUSTL(MuonDatFilePath)), STATUS='old', IOSTAT=istat)
    IF (istat .ne. 0) THEN
      WRITE(*,*) 'Cannot open Muon table.dat!'
      STOP
    ENDIF
    
    DO iT=1,MuonTable % nPointsTemp
      IF (iT .eq. (nDenMuon+1)*iT) THEN
        READ(1234,*)
      ENDIF
      DO iDen=1,MuonTable % nPointsDen
        READ(1234,*) MuonTable % t(iT), MuonTable % rhoym(iT,iDen), &
          MuonTable % p(iT,iDen), MuonTable % e(iT,iDen), &
          MuonTable % s(iT,iDen), MuonTable % mu(iT,iDen)
      ENDDO
    ENDDO
    
    ! DO UNIT CONVERSION
    ! Based on Tobias's code, the chemical potential includes the muon rest mass
    MuonTable % t = MuonTable % t * kmev_inv !MeV to kelvin
    MuonTable % rhoym = MuonTable % rhoym * rmu/cm3fm3 !baryon number to g/cm^3
    MuonTable % p = MuonTable % p * ergmev / cm3fm3 !MeV/fm^3 to dyn
    MuonTable % e = MuonTable % e * ergmev / cm3fm3 / MuonTable % rhoym
    MuonTable % s = MuonTable % s / (kmev * ergmev / rmu)
    
    ! NOW CALCULATE DERIVATIVES AT CONSTANT RHO
    DO iDen=1, MuonTable % nPointsDen
      MuonTable % dlnsdlnT(:,iDen) = Gradient3pts( LOG10(MuonTable % s(:,iDen)), LOG10(MuonTable % t), MuonTable % nPointsTemp)
      MuonTable % dlnPdlnT(:,iDen) = Gradient3pts( LOG10(MuonTable % p(:,iDen)), LOG10(MuonTable % t), MuonTable % nPointsTemp)
      MuonTable % dlnedlnT(:,iDen) = Gradient3pts( LOG10(MuonTable % e(:,iDen)), LOG10(MuonTable % t), MuonTable % nPointsTemp)
    ENDDO

    ! NOW CALCULATE DERIVATIVES AT CONSTANT T
    DO iT=1,MuonTable % nPointsTemp
      MuonTable % dlnsdlnrho(iT,:) = Gradient3pts( LOG10(MuonTable % s(iT,:)), LOG10(MuonTable % rhoym), MuonTable % nPointsDen)
      MuonTable % dlnPdlnrho(iT,:) = Gradient3pts( LOG10(MuonTable % p(iT,:)), LOG10(MuonTable % rhoym), MuonTable % nPointsDen)
      MuonTable % dlnedlnrho(iT,:) = Gradient3pts( LOG10(MuonTable % e(iT,:)), LOG10(MuonTable % rhoym), MuonTable % nPointsDen)
    ENDDO
    
    WRITE(*,*) 'Bounds of the Muon Table'
    WRITE(*,*) MAXVAL(MuonTable % t), MINVAL(MuonTable % t)
    WRITE(*,*) MAXVAL(MuonTable % rhoym), MINVAL(MuonTable % rhoym)
    WRITE(*,*) MAXVAL(MuonTable % p), MINVAL(MuonTable % p)
    WRITE(*,*) MAXVAL(MuonTable % e), MINVAL(MuonTable % e)
    WRITE(*,*) MAXVAL(MuonTable % s), MINVAL(MuonTable % s)
    WRITE(*,*) MAXVAL(MuonTable % mu), MINVAL(MuonTable % mu)
    
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

END MODULE wlLeptonEOSTableModule
        

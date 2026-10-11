MODULE wlLeptonEOSTableModule

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    ergmev, cm3fm3, kmev, kmev_inv, rmu, pi, hbarc

  IMPLICIT NONE
  PRIVATE

  ! nPoints for classic HelmTable
  INTEGER, PARAMETER, PUBLIC :: iTempMax=201, iDenMax=541
  INTEGER, PARAMETER, PUBLIC :: nPairExtensionPoints=1025
  INTEGER, PARAMETER :: nPairQuadraturePoints=128

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

    ! Pair-symmetric thermodynamics.  These one-dimensional arrays regularize
    ! the logarithmic net-density table throughout the low-degeneracy regime,
    ! including rho*|Y_L|=0.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairTemperature
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairLogP0
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairDLogP0dT
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairLogA0
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairDLogA0dT

    ! miniDenm and maxiDenm possible densities
    REAL(DP) :: mintemp
    REAL(DP) :: maxtemp
    REAL(DP) :: mindens
    REAL(DP) :: maxdens

    REAL(DP) :: lepton_mass ! To distinguish electrons and muons

  END TYPE HelmTableType

  PUBLIC AllocateHelmholtzTable
  PUBLIC DeallocateHelmholtzTable
  PUBLIC ReadHelmEOSdat
  PUBLIC InitializePairSymmetricExtension

CONTAINS

  SUBROUTINE AllocateHelmholtzTable( HelmTable, nPoints )

    TYPE(HelmTableType), INTENT(INOUT) :: HelmTable
    INTEGER, DIMENSION(2), INTENT(IN)  :: nPoints

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

    ALLOCATE( HelmTable % PairTemperature( nPairExtensionPoints ) )
    ALLOCATE( HelmTable % PairLogP0      ( nPairExtensionPoints ) )
    ALLOCATE( HelmTable % PairDLogP0dT   ( nPairExtensionPoints ) )
    ALLOCATE( HelmTable % PairLogA0      ( nPairExtensionPoints ) )
    ALLOCATE( HelmTable % PairDLogA0dT   ( nPairExtensionPoints ) )

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

    DEALLOCATE( HelmTable % PairTemperature )
    DEALLOCATE( HelmTable % PairLogP0 )
    DEALLOCATE( HelmTable % PairDLogP0dT )
    DEALLOCATE( HelmTable % PairLogA0 )
    DEALLOCATE( HelmTable % PairDLogA0dT )

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

    CALL InitializePairSymmetricExtension( HelmTable )

  END  SUBROUTINE ReadHelmEOSdat

  SUBROUTINE InitializePairSymmetricExtension( HelmTable )

    TYPE(HelmTableType), INTENT(INOUT) :: HelmTable

    INTEGER :: iT
    REAL(DP) :: LogTMin, dLogT
    REAL(DP) :: Nodes(nPairQuadraturePoints)
    REAL(DP) :: Weights(nPairQuadraturePoints)

    LogTMin = LOG10( HelmTable % mintemp )
    dLogT = ( LOG10(HelmTable % maxtemp) - LogTMin ) &
            / REAL( nPairExtensionPoints - 1, DP )

    CALL ComputeGaussLegendreRule &
           ( nPairQuadraturePoints, Nodes, Weights )

    DO iT = 1, nPairExtensionPoints
      HelmTable % PairTemperature(iT) &
        = 10.0_DP**( LogTMin + REAL(iT-1,DP) * dLogT )
      CALL ComputePairSymmetricScales &
             ( HelmTable % PairTemperature(iT), HelmTable % lepton_mass, &
               Nodes, Weights, &
               HelmTable % PairLogP0(iT), HelmTable % PairLogA0(iT) )
    END DO

    CALL ComputeNodalDerivative &
           ( HelmTable % PairTemperature, HelmTable % PairLogP0, &
             HelmTable % PairDLogP0dT )
    CALL ComputeNodalDerivative &
           ( HelmTable % PairTemperature, HelmTable % PairLogA0, &
             HelmTable % PairDLogA0dT )

  END SUBROUTINE InitializePairSymmetricExtension

  SUBROUTINE ComputePairSymmetricScales &
    ( Temperature, LeptonMass, Nodes, Weights, LogP0, LogA0 )

    REAL(DP), INTENT(IN)  :: Temperature, LeptonMass
    REAL(DP), INTENT(IN)  :: Nodes(nPairQuadraturePoints)
    REAL(DP), INTENT(IN)  :: Weights(nPairQuadraturePoints)
    REAL(DP), INTENT(OUT) :: LogP0, LogA0

    REAL(DP), PARAMETER :: LogTwo = LOG( 2.0_DP )
    REAL(DP), PARAMETER :: LogThree = LOG( 3.0_DP )
    REAL(DP), PARAMETER :: LogPi = LOG( pi )

    INTEGER :: iQ
    REAL(DP) :: a, E, f, IntP, IntA, theta, x, z
    REAL(DP) :: LogK2, LogCommon

    theta = kmev * Temperature
    a     = LeptonMass / theta

    IF( a > 50.0_DP )THEN

      ! Maxwell--Boltzmann limit.  Work in logarithms because the symmetric
      ! muon-pair background is exponentially small at low temperature.
      LogK2 = 0.5_DP * LOG( pi / ( 2.0_DP * a ) ) - a &
              + LOG( 1.0_DP + 15.0_DP / ( 8.0_DP * a ) &
                     + 105.0_DP / ( 128.0_DP * a * a ) )
      LogCommon = LOG( 2.0_DP / ( pi * pi ) ) &
                  + 2.0_DP * LOG( LeptonMass ) + LOG( theta ) &
                  - 3.0_DP * LOG( hbarc ) + LogK2
      LogP0 = LogCommon + LOG( theta ) + LOG( ergmev )
      LogA0 = LogCommon

    ELSE

      IntP = 0.0_DP
      IntA = 0.0_DP
      DO iQ = 1, nPairQuadraturePoints
        z = 0.5_DP * ( Nodes(iQ) + 1.0_DP )
        x = z / ( 1.0_DP - z )
        E = SQRT( x * x + a * a )
        IF( E < 40.0_DP )THEN
          f = 1.0_DP / ( EXP(E) + 1.0_DP )
        ELSE
          f = EXP( -E )
        END IF
        IntP = IntP + 0.5_DP * Weights(iQ) &
                      * x**4 / E * f / ( 1.0_DP - z )**2
        IntA = IntA + 0.5_DP * Weights(iQ) &
                      * x*x * f * ( 1.0_DP - f ) / ( 1.0_DP - z )**2
      END DO

      LogP0 = LogTwo - LogThree - 2.0_DP * LogPi &
              + 4.0_DP * LOG(theta) - 3.0_DP * LOG(hbarc) &
              + LOG(IntP) + LOG(ergmev)
      LogA0 = LogTwo - 2.0_DP * LogPi &
              + 3.0_DP * LOG(theta) - 3.0_DP * LOG(hbarc) &
              + LOG(IntA)

    END IF

  END SUBROUTINE ComputePairSymmetricScales

  SUBROUTINE ComputeGaussLegendreRule( n, Nodes, Weights )

    INTEGER, INTENT(IN) :: n
    REAL(DP), DIMENSION(n), INTENT(OUT) :: Nodes, Weights

    INTEGER :: i, j, m
    REAL(DP) :: p1, p2, p3, pp, z, zOld

    m = ( n + 1 ) / 2
    DO i = 1, m
      z = COS( pi * ( REAL(i,DP) - 0.25_DP ) &
               / ( REAL(n,DP) + 0.5_DP ) )
      DO
        p1 = 1.0_DP
        p2 = 0.0_DP
        DO j = 1, n
          p3 = p2
          p2 = p1
          p1 = ( ( 2.0_DP * REAL(j,DP) - 1.0_DP ) * z * p2 &
                 - ( REAL(j,DP) - 1.0_DP ) * p3 ) / REAL(j,DP)
        END DO
        pp   = REAL(n,DP) * ( z * p1 - p2 ) / ( z*z - 1.0_DP )
        zOld = z
        z    = zOld - p1 / pp
        IF( ABS( z - zOld ) <= 8.0_DP * EPSILON(z) ) EXIT
      END DO

      ! Re-evaluate the derivative at the converged root.
      p1 = 1.0_DP
      p2 = 0.0_DP
      DO j = 1, n
        p3 = p2
        p2 = p1
        p1 = ( ( 2.0_DP * REAL(j,DP) - 1.0_DP ) * z * p2 &
               - ( REAL(j,DP) - 1.0_DP ) * p3 ) / REAL(j,DP)
      END DO
      pp = REAL(n,DP) * ( z * p1 - p2 ) / ( z*z - 1.0_DP )

      Nodes(i)       = -z
      Nodes(n+1-i)   =  z
      Weights(i)     = 2.0_DP / ( ( 1.0_DP - z*z ) * pp*pp )
      Weights(n+1-i) = Weights(i)
    END DO

  END SUBROUTINE ComputeGaussLegendreRule

  SUBROUTINE ComputeNodalDerivative( X, F, dFdX )

    REAL(DP), DIMENSION(:), INTENT(IN)  :: X, F
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dFdX

    INTEGER :: i, n
    REAL(DP) :: h0, h1

    n = SIZE(X)
    h0 = X(2) - X(1)
    h1 = X(3) - X(2)
    dFdX(1) = -( 2.0_DP*h0 + h1 ) / ( h0 * ( h0 + h1 ) ) * F(1) &
               + ( h0 + h1 ) / ( h0 * h1 ) * F(2) &
               - h0 / ( h1 * ( h0 + h1 ) ) * F(3)
    DO i = 2, n - 1
      h0 = X(i)   - X(i-1)
      h1 = X(i+1) - X(i)
      dFdX(i) = -h1 / ( h0 * ( h0 + h1 ) ) * F(i-1) &
                 + ( h1 - h0 ) / ( h0 * h1 ) * F(i) &
                 + h0 / ( h1 * ( h0 + h1 ) ) * F(i+1)
    END DO
    h0 = X(n-1) - X(n-2)
    h1 = X(n)   - X(n-1)
    dFdX(n) = h1 / ( h0 * ( h0 + h1 ) ) * F(n-2) &
               - ( h0 + h1 ) / ( h0 * h1 ) * F(n-1) &
               + ( h0 + 2.0_DP*h1 ) / ( h1 * ( h0 + h1 ) ) * F(n)

  END SUBROUTINE ComputeNodalDerivative

END MODULE wlLeptonEOSTableModule

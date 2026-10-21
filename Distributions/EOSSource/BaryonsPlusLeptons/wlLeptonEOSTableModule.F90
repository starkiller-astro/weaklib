MODULE wlLeptonEOSTableModule

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    ergmev, cm3fm3, kmev, kmev_inv, rmu, pi, hbarc_Default => hbarc

  IMPLICIT NONE
  PRIVATE

  ! nPoints for classic HelmTable
  INTEGER, PARAMETER, PUBLIC :: iTempMax=201, iDenMax=541
  INTEGER, PARAMETER, PUBLIC :: nPairExtensionPoints=1025
  INTEGER, PARAMETER :: nPairQuadraturePoints=128
  INTEGER, PARAMETER, PUBLIC :: nPairDensityTerms=7
  ! Explicit preset for the supplied muon Helm tables. These CGS inputs
  ! belong to this table model, not to the global weak-interaction constants.
  ! Order: rest energy, Boltzmann constant, hbar*c, number per mass.
  REAL(DP), PARAMETER, PUBLIC :: MuonHelmGeneratorConstants(4) = &
    [ 1.8835371982e-25_DP * 29979245800.0_DP**2, &
      1.3806490000000002e-16_DP, &
      1.0545718176461563e-27_DP * 29979245800.0_DP, &
      6.02214076e23_DP ]

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
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairD2LogP0dT2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairLogA0
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairDLogA0dT
    REAL(DP), DIMENSION(:), ALLOCATABLE :: PairD2LogA0dT2
    ! Normalized odd-density Taylor coefficients Q_j/A0 and their first
    ! two log-temperature derivatives, j=2..7. Q_1=A0 is stored above.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: PairDensitySeries

    ! miniDenm and maxiDenm possible densities
    REAL(DP) :: mintemp
    REAL(DP) :: maxtemp
    REAL(DP) :: mindens
    REAL(DP) :: maxdens

    REAL(DP) :: lepton_mass ! To distinguish electrons and muons
    ! Internal energy quantities use the existing runtime MeV-to-erg
    ! conversion, keeping API units unchanged. The gas itself uses the
    ! generator's physical constants. Defaults preserve legacy tables.
    REAL(DP) :: ThermalConstant = kmev
    REAL(DP) :: HbarC = hbarc_Default
    REAL(DP) :: NumberPerMass = 1.0_DP/rmu
    REAL(DP) :: DensityScale = 1.0_DP

  END TYPE HelmTableType

  PUBLIC AllocateHelmholtzTable
  PUBLIC DeallocateHelmholtzTable
  PUBLIC ReadHelmEOSdat
  PUBLIC InitializePairSymmetricExtension
  PUBLIC ComputeClassicalPairScales

CONTAINS

  SUBROUTINE ComputeClassicalPairScales &
    ( Temperature, LeptonMass, ThermalConstant, HbarC, LogA0, LA1, LA2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! Relativistic Maxwell--Boltzmann gas, including both charges. Only
    ! used for a=m/(kT)>=12. The 16-term K_2 asymptotic expansion is checked
    ! against the Bessel integral in the regression test. Logarithms retain
    ! the pair scale even when the physical pair density underflows.
    REAL(DP), INTENT(IN) :: Temperature, LeptonMass, ThermalConstant, HbarC
    REAL(DP), INTENT(OUT) :: LogA0, LA1, LA2
    INTEGER :: j
    REAL(DP) :: a, theta, term, series, series1, series2, q1, q2

    theta = ThermalConstant * Temperature
    a = LeptonMass / theta
    series = 1.0_DP
    series1 = 0.0_DP
    series2 = 0.0_DP
    term = 1.0_DP
    DO j = 1, 16
      term = term * ( 16.0_DP - REAL((2*j-1)**2,DP) ) &
                  / ( 8.0_DP * REAL(j,DP) * a )
      series = series + term
      series1 = series1 + REAL(j,DP) * term
      series2 = series2 + REAL(j*j,DP) * term
    END DO
    q1 = series1 / series
    q2 = series2 / series - q1*q1
    LogA0 = LOG(2.0_DP/(pi*pi)) + 2.0_DP*LOG(LeptonMass) &
              + LOG(theta) - 3.0_DP*LOG(HbarC) &
              + 0.5_DP*LOG(pi/(2.0_DP*a)) - a + LOG(series)
    LA1 = ( 1.5_DP + a + q1 ) / Temperature
    LA2 = ( -1.5_DP - 2.0_DP*a + q2 - q1 ) &
            / ( Temperature * Temperature )

  END SUBROUTINE ComputeClassicalPairScales

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
    ALLOCATE( HelmTable % PairD2LogP0dT2 ( nPairExtensionPoints ) )
    ALLOCATE( HelmTable % PairLogA0      ( nPairExtensionPoints ) )
    ALLOCATE( HelmTable % PairDLogA0dT   ( nPairExtensionPoints ) )
    ALLOCATE( HelmTable % PairD2LogA0dT2 ( nPairExtensionPoints ) )
    ALLOCATE( HelmTable % PairDensitySeries(3,nPairDensityTerms-1,nPairExtensionPoints) )

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
    DEALLOCATE( HelmTable % PairD2LogP0dT2 )
    DEALLOCATE( HelmTable % PairLogA0 )
    DEALLOCATE( HelmTable % PairDLogA0dT )
    DEALLOCATE( HelmTable % PairD2LogA0dT2 )
    DEALLOCATE( HelmTable % PairDensitySeries )

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
               HelmTable%ThermalConstant, HelmTable%HbarC, Nodes, Weights, &
               HelmTable % PairLogP0(iT), HelmTable % PairLogA0(iT), &
               HelmTable % PairDLogP0dT(iT), HelmTable % PairD2LogP0dT2(iT), &
               HelmTable % PairDLogA0dT(iT), HelmTable % PairD2LogA0dT2(iT) )
      CALL ComputePairDensitySeries &
             ( HelmTable%PairTemperature(iT), HelmTable%lepton_mass, HelmTable%ThermalConstant, Nodes, Weights, &
               HelmTable%PairDensitySeries(:,:,iT) )
    END DO


  END SUBROUTINE InitializePairSymmetricExtension

  SUBROUTINE ComputePairDensitySeries &
    ( Temperature, LeptonMass, ThermalConstant, Nodes, Weights, Series )

    REAL(DP), INTENT(IN) :: Temperature, LeptonMass, ThermalConstant
    REAL(DP), INTENT(IN) :: Nodes(nPairQuadraturePoints), Weights(nPairQuadraturePoints)
    REAL(DP), INTENT(OUT) :: Series(3,nPairDensityTerms-1)
    REAL(DP) :: c(0:2*nPairDensityTerms+1), sums(3,nPairDensityTerms)
    REAL(DP) :: a, z, x, E, f, weight, convolution, fact, r, rx
    INTEGER :: iQ, i, j, k, order

    a = LeptonMass/(ThermalConstant*Temperature)
    IF( a > 50.0_DP )THEN
      fact = 1.0_DP
      Series = 0.0_DP
      DO j = 2, nPairDensityTerms
        fact = fact*REAL((2*j-2)*(2*j-1),DP)
        Series(1,j-1) = 1.0_DP/fact
      END DO
      RETURN
    END IF
    sums = 0.0_DP
    DO iQ = 1, nPairQuadraturePoints
      z = 0.5_DP*(Nodes(iQ)+1.0_DP)
      x = z/(1.0_DP-z)
      E = SQRT(x*x+a*a)
      f = EXP(-E)/(1.0_DP+EXP(-E))
      c(0) = f
      ! Normalized Taylor coefficients, not large factorial derivatives.
      ! f'(eta)=f-f^2 supplies a stable power-series recurrence.
      DO i = 1, 2*nPairDensityTerms+1
        convolution = 0.0_DP
        DO k = 0, i-1
          convolution = convolution+c(k)*c(i-1-k)
        END DO
        c(i) = (c(i-1)-convolution)/REAL(i,DP)
      END DO
      weight = 0.5_DP*Weights(iQ)*x*x/(1.0_DP-z)**2
      DO j = 1, nPairDensityTerms
        order = 2*j-1
        sums(1,j) = sums(1,j)+weight*c(order)
        sums(2,j) = sums(2,j)+weight*E*REAL(order+1,DP)*c(order+1)
        sums(3,j) = sums(3,j)+weight*(E*E*REAL((order+1)*(order+2),DP)*c(order+2) &
                                      -E*REAL(order+1,DP)*c(order+1))
      END DO
    END DO
    DO j = 2, nPairDensityTerms
      r = sums(1,j)/sums(1,1)
      rx = (sums(2,j)-r*sums(2,1))/sums(1,1)
      Series(1,j-1) = r
      Series(2,j-1) = rx
      Series(3,j-1) = (sums(3,j)-r*sums(3,1)-2.0_DP*rx*sums(2,1))/sums(1,1)
    END DO

  END SUBROUTINE ComputePairDensitySeries

  SUBROUTINE ComputePairSymmetricScales &
    ( Temperature, LeptonMass, ThermalConstant, HbarC, Nodes, Weights, LogP0, LogA0, LP1, LP2, LA1, LA2 )

    REAL(DP), INTENT(IN)  :: Temperature, LeptonMass, ThermalConstant, HbarC
    REAL(DP), INTENT(IN)  :: Nodes(nPairQuadraturePoints)
    REAL(DP), INTENT(IN)  :: Weights(nPairQuadraturePoints)
    REAL(DP), INTENT(OUT) :: LogP0, LogA0
    REAL(DP), INTENT(OUT) :: LP1,LP2,LA1,LA2

    REAL(DP), PARAMETER :: LogTwo = LOG( 2.0_DP )
    REAL(DP), PARAMETER :: LogThree = LOG( 3.0_DP )
    REAL(DP), PARAMETER :: LogPi = LOG( pi )

    INTEGER :: iQ,j
    REAL(DP) :: a, E, f, IntP, IntA, theta, x, z
    REAL(DP) :: LogK2, LogCommon
    REAL(DP) :: IntU,IntUT,IntA1,IntA2,Weight,ff,Term,Series,Series1,Series2,q1,q2

    theta = ThermalConstant * Temperature
    a     = LeptonMass / theta

    IF( a > 50.0_DP )THEN

      ! Maxwell--Boltzmann limit.  Work in logarithms because the symmetric
      ! muon-pair background is exponentially small at low temperature.
      Series=1.0_DP
      Series1=0.0_DP
      Series2=0.0_DP
      Term=1.0_DP
      DO j=1,8
        Term=Term*(16.0_DP-REAL((2*j-1)**2,DP))/(8.0_DP*REAL(j,DP)*a)
        Series=Series+Term
        Series1=Series1+REAL(j,DP)*Term
        Series2=Series2+REAL(j*j,DP)*Term
      END DO
      q1=Series1/Series
      q2=Series2/Series-q1*q1
      LogK2 = 0.5_DP * LOG( pi / ( 2.0_DP * a ) ) - a + LOG(Series)
      LogCommon = LOG( 2.0_DP / ( pi * pi ) ) &
                  + 2.0_DP * LOG( LeptonMass ) + LOG( theta ) &
                  - 3.0_DP * LOG( hbarc ) + LogK2
      LogP0 = LogCommon + LOG( theta ) + LOG( ergmev )
      LogA0 = LogCommon
      LA1=(1.5_DP+a+q1)/Temperature
      LA2=(-1.5_DP-2.0_DP*a+q2-q1)/(Temperature*Temperature)
      LP1=LA1+1.0_DP/Temperature
      LP2=LA2-1.0_DP/(Temperature*Temperature)

    ELSE

      IntP = 0.0_DP
      IntA = 0.0_DP
      IntU = 0.0_DP
      IntUT = 0.0_DP
      IntA1 = 0.0_DP
      IntA2 = 0.0_DP
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
        Weight=0.5_DP*Weights(iQ)*x*x/(1.0_DP-z)**2
        ff=f*(1.0_DP-f)
        IntU=IntU+Weight*E*f
        IntUT=IntUT+Weight*E*E*ff
        IntA1=IntA1+Weight*E*ff*(1.0_DP-2.0_DP*f)
        IntA2=IntA2+Weight*ff*(E*E*(1.0_DP-6.0_DP*f+6.0_DP*f*f) &
                    -2.0_DP*E*(1.0_DP-2.0_DP*f))
      END DO

      LogP0 = LogTwo - LogThree - 2.0_DP * LogPi &
              + 4.0_DP * LOG(theta) - 3.0_DP * LOG(HbarC) &
              + LOG(IntP) + LOG(ergmev)
      LogA0 = LogTwo - 2.0_DP * LogPi &
              + 3.0_DP * LOG(theta) - 3.0_DP * LOG(HbarC) &
              + LOG(IntA)
      LP1=(1.0_DP+3.0_DP*IntU/IntP)/Temperature
      LP2=(3.0_DP*IntUT/IntP)/(Temperature*Temperature)-LP1*LP1
      LA1=(IntA1/IntA)/Temperature
      LA2=(IntA2/IntA)/(Temperature*Temperature)-LA1*LA1

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


END MODULE wlLeptonEOSTableModule

MODULE wlLeptonPhotonGasEOS

  USE wlKindModule, ONLY: dp
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType, ComputeClassicalPairScales, nPairDensityTerms
  USE wlHelmholtzEOS
  USE wlEosConstantsModule, ONLY: &
    pi, rmu, kerg, cvel, &
    kmev, ergmev, sigma_sb, asol

  IMPLICIT NONE
  PRIVATE

  ! A grand-pressure expansion through eta^14 supplies the regular Fermi
  ! pair limit. Join its Helmholtz free energy to the regularized table;
  ! all outputs, including mu and compressibility, use that same potential.
  ! These thresholds refer to asinh(n/A0), not to a floor on net density.
  REAL(DP), PARAMETER, PUBLIC :: PairSymmetricEtaMin = 0.5_DP
  REAL(DP), PARAMETER, PUBLIC :: PairSymmetricEtaMax = 1.0_DP

  ! Derivatives of one free-energy density F(T,n), n the signed net
  ! lepton number density. Store density derivatives in scaled form to
  ! avoid forming 1/n**2 near pair symmetry. Fn remains unscaled so that
  ! a representable chemical potential is not lost when n*Fn underflows.
  TYPE, PRIVATE :: LeptonFreeEnergyType
    REAL(DP) :: F, FT, FTT, Fn, nFTn, nnFnn
  END TYPE LeptonFreeEnergyType

  TYPE, PUBLIC :: LeptonGasType

    REAL(DP) :: T
    REAL(DP) :: rho
    REAL(DP) :: yL

    REAL(DP) :: p
    REAL(DP) :: dpdT
    REAL(DP) :: dpdr

    REAL(DP) :: e
    REAL(DP) :: dedT
    REAL(DP) :: dedr

    REAL(DP) :: s
    REAL(DP) :: dsdT
    REAL(DP) :: dsdr

    REAL(DP) :: mu

  END TYPE LeptonGasType

  TYPE, PUBLIC :: PhotonGasType

    REAL(DP) :: T
    REAL(DP) :: rho

    REAL(DP) :: p
    REAL(DP) :: dpdT
    REAL(DP) :: dpdr

    REAL(DP) :: e
    REAL(DP) :: dedT
    REAL(DP) :: dedr

    REAL(DP) :: s
    REAL(DP) :: dsdT
    REAL(DP) :: dsdr

  END TYPE PhotonGasType

  PUBLIC :: LeptonGasEOS, PhotonGasEOS, GetPhotonLeptonGasEOS
  PUBLIC :: ComputePairDegeneracy, ComputePairLeptonFraction

CONTAINS

  ! Wrapper to make this call easy
  SUBROUTINE GetPhotonLeptonGasEOS(D, T, Ye, Ym, &
    HelmTableElectrons, HelmTableMuons, P, E, S)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(IN)  :: D, T, Ye, Ym
    TYPE(HelmTableType), INTENT(IN) :: HelmTableElectrons
    TYPE(HelmTableType), INTENT(IN) :: HelmTableMuons
    REAL(DP), INTENT(OUT) :: P, E, S

    ! LOCAL
    TYPE(PhotonGasType) :: PhotonGasState
    TYPE(LeptonGasType) :: LeptonGasState

    P = 0.0_dp
    E = 0.0_dp
    S = 0.0_dp

    ! Initialize your states
    PhotonGasState % rho = D
    PhotonGasState % T   = T

    ! Do photons first
    CALL PhotonGasEOS(PhotonGasState)

    P = P + PhotonGasState % p
    E = E + PhotonGasState % e
    S = S + PhotonGasState % s

    ! Then electrons
    LeptonGasState % rho = D
    LeptonGasState % T   = T
    LeptonGasState % yL  = Ye

    CALL LeptonGasEOS(HelmTableElectrons, LeptonGasState)

    P = P + LeptonGasState % p
    E = E + LeptonGasState % e
    S = S + LeptonGasState % s

    ! Then muons
    LeptonGasState % rho = D
    LeptonGasState % T   = T
    LeptonGasState % yL  = Ym

    CALL LeptonGasEOS(HelmTableMuons, LeptonGasState)

    P = P + LeptonGasState % p
    E = E + LeptonGasState % e
    S = S + LeptonGasState % s

  END SUBROUTINE GetPhotonLeptonGasEOS

  SUBROUTINE PhotonGasEOS(PhotonGasState)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    TYPE(PhotonGasType), INTENT(INOUT) :: PhotonGasState

    ! Local variables
    REAL(DP) :: temp, den, tempi, deni
    REAL(DP) :: prad, erad, srad, &
                dpraddd, dpraddt, &
                deraddd, deraddt, &
                dsraddd, dsraddt

    temp = PhotonGasState % T
    den  = PhotonGasState % rho

    deni    = 1.0_dp/den
    tempi   = 1.0_dp/temp

    prad    = asol/3.0_dp * temp * temp * temp * temp
    dpraddd = 0.0_dp
    dpraddt = 4.0_dp * prad*tempi

    erad    = 3.0_dp * prad*deni
    deraddd = -erad*deni
    deraddt = 3.0_dp * dpraddt*deni

    srad    = (prad*deni + erad)*tempi
    dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
    dsraddt = (dpraddt*deni + deraddt - srad)*tempi

    PhotonGasState % p    = prad
    PhotonGasState % dpdT = dpraddt
    PhotonGasState % dpdr = dpraddd

    PhotonGasState % e    = erad
    PhotonGasState % dedT = deraddt
    PhotonGasState % dedr = deraddd

    PhotonGasState % s    = srad    / (kmev * ergmev / rmu)
    PhotonGasState % dsdT = dsraddt / (kmev * ergmev / rmu)
    PhotonGasState % dsdr = dsraddd / (kmev * ergmev / rmu)

  END SUBROUTINE PhotonGasEOS

  SUBROUTINE LeptonGasEOS(HelmTable, LeptonGasState)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    !..input arguments
    TYPE(HelmTableType), INTENT(IN)    :: HelmTable
    TYPE(LeptonGasType), INTENT(INOUT) :: LeptonGasState

    !..declare local variables

    REAL(DP) :: temp, den, yL, din, x, deni, tempi, &
    dplep_dt,  dplep_dd,  delep_dt, delep_dd, dslep_dd, dslep_dt,  &
    plep, elep, slep, etalep

    !..for the interpolations
    INTEGER  :: iat, jat, iCornerD, iCornerT, iCorner, iDerivative
    REAL(DP) :: free, df_d, df_t, df_tt, df_dt, df_dd
    REAL(DP) :: xt, xd, mxt, mxd, fi(36),   &
    si0t, si1t, si2t, si0mt, si1mt, si2mt,  &
    si0d, si1d, si2d, si0md, si1md, si2md,  &
    dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt,  &
    dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md,  &
    ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
    ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md

    ! extra variables I added that were previously allocated
    REAL(DP) :: smallt, smalld, hight, highd
    REAL(DP) :: tstp, dstp, dstpi, tstpi
    TYPE(HelmholtzStateType) :: HelmState
    TYPE(LeptonFreeEnergyType) :: TableFreeEnergy
    TYPE(LeptonFreeEnergyType) :: ReferenceFreeEnergy
    TYPE(LeptonGasType) :: ReferenceState
    REAL(DP) :: ReferenceDerivatives(9)
    REAL(DP) :: etaPair, sign_yL
    REAL(DP) :: DiluteWeight(6)
    TYPE(LeptonGasType) :: DiluteState
    TYPE(LeptonFreeEnergyType) :: DilutePotential

    temp = LeptonGasState % T
    den  = LeptonGasState % rho
    yL   = LeptonGasState % yL

    ! Handle negative Ym case
    sign_yL = SIGN(1.0d0, yL)
    yL = ABS(yL)

    ! Below the temperature table, retain the cold dilute limit.
    IF ( temp < HelmTable % mintemp ) THEN

      LeptonGasState % p    = 0.0_dp
      LeptonGasState % dpdT = 0.0_dp
      LeptonGasState % dpdr = 0.0_dp

      LeptonGasState % e    = HelmTable % lepton_mass / rmu * ergmev * yL
      LeptonGasState % dedT = 0.0_dp
      LeptonGasState % dedr = 0.0_dp

      LeptonGasState % s    = 0.0_dp
      LeptonGasState % dsdT = 0.0_dp
      LeptonGasState % dsdr = 0.0_dp

      IF( yL > 0.0_dp )THEN
        LeptonGasState % mu = sign_yL * HelmTable % lepton_mass
      ELSE
        LeptonGasState % mu = 0.0_dp
      END IF
      RETURN

    ENDIF

    IF ( (yL*den*HelmTable%DensityScale) > HelmTable % maxdens .OR. &
         temp   > HelmTable % maxtemp ) THEN

      STOP 'LeptonGasEOS: density or temperature is too high for the table'
      RETURN

    ENDIF

    ! In the nondegenerate massive-gas regime the analytic relativistic
    ! gas is more accurate than interpolating a cancelling rest-mass term.
    ! The common-potential transition is controlled by occupation, not by
    ! an arbitrary net-density floor or the lower edge of a particular file.
    CALL ComputeDilutePairWeight &
           ( HelmTable, temp, den, sign_yL*yL, DiluteWeight )
    IF( DiluteWeight(1) > 0.0_DP )THEN
      CALL EvaluatePairSymmetricExtension &
             ( HelmTable, temp, den, sign_yL*yL, DiluteState, &
               DilutePotential, Classical_Option=.TRUE. )
      IF( DiluteWeight(1) >= 1.0_DP )THEN
        LeptonGasState = DiluteState
        RETURN
      END IF
    END IF

    ! The logarithmic net-density coordinate is ill-conditioned near pair
    ! symmetry: the chemical potential is recovered by cancelling its rest
    ! mass against the tabulated degeneracy parameter.  Use the regular
    ! pair-symmetric representation throughout the nondegenerate interval,
    ! not merely below the first net-density table point.
    CALL ComputePairDegeneracy &
           ( HelmTable, temp, den, sign_yL * yL, etaPair )
    IF ( (yL*den*HelmTable%DensityScale) <= HelmTable % mindens .OR. &
         ABS(etaPair) <= PairSymmetricEtaMin ) THEN
      IF( ABS(etaPair) <= PairSymmetricEtaMax )THEN
        CALL EvaluateFermiPairExtension &
               ( HelmTable, temp, den, sign_yL*yL, LeptonGasState, TableFreeEnergy )
      ELSE
        CALL EvaluatePairSymmetricExtension &
               ( HelmTable, temp, den, sign_yL*yL, LeptonGasState, TableFreeEnergy )
      END IF
      IF( DiluteWeight(1) > 0.0_DP )THEN
        CALL BlendLeptonFreeEnergy &
               ( den*sign_yL*yL/rmu, DiluteWeight, DilutePotential, TableFreeEnergy )
        CALL SetLeptonStateFromFreeEnergy &
               ( temp, den, sign_yL*yL, TableFreeEnergy, LeptonGasState )
      END IF
      RETURN
    END IF

    smallt = HelmTable % mintemp
    smalld = HelmTable % mindens
    hight  = HelmTable % maxtemp
    highd  = HelmTable % maxdens

    tstp  = (LOG10(hight) - LOG10(smallt))/real(HelmTable % nPointsTemp -1, dp)
    dstp  = (LOG10(highd) - LOG10(smalld))/real(HelmTable % nPointsDen -1,  dp)

    tstpi = 1.0_dp / tstp
    dstpi = 1.0_dp / dstp

    din   = (yL*den*HelmTable%DensityScale)

    !..initialize
    deni    = 1.0_dp/den
    tempi   = 1.0_dp/temp

    !..ion section: removed in this

    !..enter the table with (yL*den*HelmTable%DensityScale)
    din = (yL*den*HelmTable%DensityScale)

    !..hash locate this temperature and density
    jat = INT( (LOG10(temp) - LOG10(smallt))*tstpi ) + 1
    jat = MAX( 1, MIN(jat, HelmTable % nPointsTemp-1) )
    iat = INT( (LOG10(din ) - LOG10(smalld))*dstpi ) + 1
    iat = MAX( 1, MIN(iat, HelmTable % nPointsDen-1) )

    fi(:) = 0.0_dp
    !..access the table locations only once
    fi(1)  = HelmTable % f(iat,jat)
    fi(2)  = HelmTable % f(iat+1,jat)
    fi(3)  = HelmTable % f(iat,jat+1)
    fi(4)  = HelmTable % f(iat+1,jat+1)
    fi(5)  = HelmTable % ft(iat,jat)
    fi(6)  = HelmTable % ft(iat+1,jat)
    fi(7)  = HelmTable % ft(iat,jat+1)
    fi(8)  = HelmTable % ft(iat+1,jat+1)
    fi(9)  = HelmTable % ftt(iat,jat)
    fi(10) = HelmTable % ftt(iat+1,jat)
    fi(11) = HelmTable % ftt(iat,jat+1)
    fi(12) = HelmTable % ftt(iat+1,jat+1)
    fi(13) = HelmTable % fd(iat,jat)
    fi(14) = HelmTable % fd(iat+1,jat)
    fi(15) = HelmTable % fd(iat,jat+1)
    fi(16) = HelmTable % fd(iat+1,jat+1)
    fi(17) = HelmTable % fdd(iat,jat)
    fi(18) = HelmTable % fdd(iat+1,jat)
    fi(19) = HelmTable % fdd(iat,jat+1)
    fi(20) = HelmTable % fdd(iat+1,jat+1)
    fi(21) = HelmTable % fdt(iat,jat)
    fi(22) = HelmTable % fdt(iat+1,jat)
    fi(23) = HelmTable % fdt(iat,jat+1)
    fi(24) = HelmTable % fdt(iat+1,jat+1)
    fi(25) = HelmTable % fddt(iat,jat)
    fi(26) = HelmTable % fddt(iat+1,jat)
    fi(27) = HelmTable % fddt(iat,jat+1)
    fi(28) = HelmTable % fddt(iat+1,jat+1)
    fi(29) = HelmTable % fdtt(iat,jat)
    fi(30) = HelmTable % fdtt(iat+1,jat)
    fi(31) = HelmTable % fdtt(iat,jat+1)
    fi(32) = HelmTable % fdtt(iat+1,jat+1)
    fi(33) = HelmTable % fddtt(iat,jat)
    fi(34) = HelmTable % fddtt(iat+1,jat)
    fi(35) = HelmTable % fddtt(iat,jat+1)
    fi(36) = HelmTable % fddtt(iat+1,jat+1)

    ! Interpolate the departure from the analytic pair gas. The original
    ! specific potential contains -P0/d and a classical log(d) term;
    ! interpolating those large terms directly need not retain convexity
    ! on the coarse first density interval. Subtract the SAME reference
    ! and all its derivatives at the corners, and add it back exactly.
    DO iCornerT=0,1
    DO iCornerD=0,1
      iCorner=1+iCornerD+2*iCornerT
      CALL EvaluatePairSymmetricExtension &
        (HelmTable,HelmTable%t(jat+iCornerT),HelmTable%d(iat+iCornerD)/HelmTable%DensityScale, &
         1.0_DP,ReferenceState,SpecificDerivatives_Option=ReferenceDerivatives)
      DO iDerivative=1,9
        fi(iCorner+4*(iDerivative-1))=fi(iCorner+4*(iDerivative-1))-ReferenceDerivatives(iDerivative)
      END DO
    END DO
    END DO

    !..various differences
    xt  = MAX( (temp - HelmTable % t(jat)) * HelmTable % dti(jat), 0.0_dp)
    xd  = MAX( (din  - HelmTable % d(iat)) * HelmTable % ddi(iat), 0.0_dp)
    mxt = 1.0_dp - xt
    mxd = 1.0_dp - xd

    !..the six density and six temperature basis functions
    si0t =   psi0(xt)
    si1t =   psi1(xt)*HelmTable % dt(jat)
    si2t =   psi2(xt)*HelmTable % dt2(jat)

    si0mt =  psi0(mxt)
    si1mt = -psi1(mxt)*HelmTable % dt(jat)
    si2mt =  psi2(mxt)*HelmTable % dt2(jat)

    si0d =   psi0(xd)
    si1d =   psi1(xd)*HelmTable % dd(iat)
    si2d =   psi2(xd)*HelmTable % dd2(iat)

    si0md =  psi0(mxd)
    si1md = -psi1(mxd)*HelmTable % dd(iat)
    si2md =  psi2(mxd)*HelmTable % dd2(iat)

    !..derivatives of the weight functions
    dsi0t =   dpsi0(xt )*HelmTable % dti(jat)
    dsi1t =   dpsi1(xt )
    dsi2t =   dpsi2(xt )*HelmTable % dt (jat)

    dsi0mt = -dpsi0(mxt)*HelmTable % dti(jat)
    dsi1mt =  dpsi1(mxt)
    dsi2mt = -dpsi2(mxt)*HelmTable % dt (jat)

    dsi0d =   dpsi0(xd )*HelmTable % ddi(iat)
    dsi1d =   dpsi1(xd )
    dsi2d =   dpsi2(xd )*HelmTable % dd (iat)

    dsi0md = -dpsi0(mxd)*HelmTable % ddi(iat)
    dsi1md =  dpsi1(mxd)
    dsi2md = -dpsi2(mxd)*HelmTable % dd (iat)

    !..second derivatives of the weight functions
    ddsi0t =   ddpsi0(xt)*HelmTable % dt2i(jat)
    ddsi1t =   ddpsi1(xt)*HelmTable % dti (jat)
    ddsi2t =   ddpsi2(xt)

    ddsi0mt =  ddpsi0(mxt)*HelmTable % dt2i(jat)
    ddsi1mt = -ddpsi1(mxt)*HelmTable % dti (jat)
    ddsi2mt =  ddpsi2(mxt)

    ddsi0d  =  ddpsi0(xd)*HelmTable % dd2i(iat)
    ddsi1d  =  ddpsi1(xd)*HelmTable % ddi (iat)
    ddsi2d  =  ddpsi2(xd)
    ddsi0md =  ddpsi0(mxd)*HelmTable % dd2i(iat)
    ddsi1md = -ddpsi1(mxd)*HelmTable % ddi (iat)
    ddsi2md =  ddpsi2(mxd)

    !..the free energy
    free  = h5( fi, &
      si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
      si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    !..derivative with respect to density
    df_d  = h5( fi, &
      si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
      dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

    !..derivative with respect to temperature
    df_t = h5( fi, &
      dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
      si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    !..derivative with respect to temperature**2
    df_tt = h5( fi, &
      ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
      si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    !..derivative with respect to temperature and density
    df_dt = h5( fi, &
      dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
      dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

    df_dd = h5( fi, &
      si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
      ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

    ! The tabulated Helmholtz potential is per unit net lepton mass
    ! density d=rho*|Y|. Include rest energy before differentiating.
    ! Chemical potential and compressibility must come from this same
    ! interpolant, not independent bicubic ef/dpdf tables at the join.
    CALL EvaluatePairSymmetricExtension &
      (HelmTable,temp,den,sign_yL*yL,ReferenceState,ReferenceFreeEnergy)
    TableFreeEnergy % F = din*free + ReferenceFreeEnergy%F
    TableFreeEnergy % FT = din*df_t + ReferenceFreeEnergy%FT
    TableFreeEnergy % FTT = din*df_tt + ReferenceFreeEnergy%FTT
    TableFreeEnergy % Fn = sign_yL &
      / HelmTable%NumberPerMass*(free+din*df_d) + ReferenceFreeEnergy%Fn
    TableFreeEnergy % nFTn = din*(df_t+din*df_dt)+ReferenceFreeEnergy%nFTn
    TableFreeEnergy % nnFnn = din*din*(2.0_DP*df_d+din*df_dd)+ReferenceFreeEnergy%nnFnn

    IF( ABS(etaPair) < PairSymmetricEtaMax .OR. &
        (yL*den*HelmTable%DensityScale) < HelmTable%d(2) )THEN
      CALL BlendPairSymmetricFreeEnergy &
        (HelmTable,temp,den,sign_yL*yL,etaPair,TableFreeEnergy)
    END IF
    IF( DiluteWeight(1) > 0.0_DP )THEN
      CALL BlendLeptonFreeEnergy &
             ( den*sign_yL*yL/rmu, DiluteWeight, DilutePotential, TableFreeEnergy )
    END IF
    CALL SetLeptonStateFromFreeEnergy &
      (temp,den,sign_yL*yL,TableFreeEnergy,LeptonGasState)

    ! There are some negatives in the table
    IF (LeptonGasState % p < 0.0_dp .OR. &
        LeptonGasState % e < 0.0_dp .OR. &
        LeptonGasState % s < 0.0_dp ) THEN

      WRITE(*,*) 'Something is negative, and it should not be'
      WRITE(*,*) LeptonGasState % rho, LeptonGasState % T, LeptonGasState % yL
      WRITE(*,*) den, temp, yL
      WRITE(*,*) LeptonGasState % p  , LeptonGasState % e, LeptonGasState % s

      HelmState % rho = LeptonGasState % rho
      HelmState % T   = LeptonGasState % T
      HelmState % ye  = LeptonGasState % yL
      CALL FullHelmEOS(1, HelmTable, HelmState)
      WRITE(*,*) HelmState % pele  , HelmState % eele, HelmState % sele
      WRITE(*,*) iat, jat, df_d, -df_t, yL
      WRITE(*,*) HelmTable % nPointsDen, HelmTable % nPointsTemp
      STOP

      ! Or if you do not want to STOP perhaps set everything to zero?
      ! But negative pressures is a bad sign
      LeptonGasState % p    = 0.0_dp
      LeptonGasState % dpdT = 0.0_dp
      LeptonGasState % dpdr = 0.0_dp

      LeptonGasState % e    = 0.0_dp
      LeptonGasState % dedT = 0.0_dp
      LeptonGasState % dedr = 0.0_dp

      LeptonGasState % s    = 0.0_dp
      LeptonGasState % dsdT = 0.0_dp
      LeptonGasState % dsdr = 0.0_dp

      LeptonGasState % mu   = 0.0_dp

    ENDIF

  END SUBROUTINE LeptonGasEOS

  SUBROUTINE ComputePairDegeneracy &
    ( HelmTable, Temperature, Density, SignedYL, Eta )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    TYPE(HelmTableType), INTENT(IN) :: HelmTable
    REAL(DP), INTENT(IN)  :: Temperature, Density, SignedYL
    REAL(DP), INTENT(OUT) :: Eta

    INTEGER :: jT
    REAL(DP) :: A0, dLogA0dT, d2LogA0dT2, LogA0, nNet

    CALL InterpolatePairScale &
           ( HelmTable % PairTemperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable % PairLogA0, &
             HelmTable % PairDLogA0dT, HelmTable % PairD2LogA0dT2, Temperature, jT, &
             A0, LogA0, dLogA0dT, d2LogA0dT2 )

    nNet = Density * SignedYL / rmu
    Eta = PairDegeneracyFromScale(nNet,A0,LogA0)

  END SUBROUTINE ComputePairDegeneracy


  REAL(DP) FUNCTION PairDegeneracyFromScale(nNet,A0,LogA0) RESULT(Eta)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: nNet,A0,LogA0
    IF( nNet == 0.0_DP )THEN
      Eta=0.0_DP
    ELSE IF( ABS(nNet) <= A0 )THEN
      Eta=ASINH(nNet/A0)
    ELSE
      ! asinh(n/A) without overflowing n/A, even if A is subnormal.
      Eta=SIGN(LOG(ABS(nNet))-LogA0+LOG(1.0_DP+HYPOT(1.0_DP,A0/nNet)),nNet)
    END IF
  END FUNCTION PairDegeneracyFromScale


  SUBROUTINE ComputePairLeptonFraction &
    ( HelmTable, Temperature, Density, Eta, SignedYL, dYLdEta )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    TYPE(HelmTableType), INTENT(IN) :: HelmTable
    REAL(DP), INTENT(IN)  :: Temperature, Density, Eta
    REAL(DP), INTENT(OUT) :: SignedYL, dYLdEta

    INTEGER :: jT
    REAL(DP) :: A0, dLogA0dT, d2LogA0dT2, LogA0, Scale
    REAL(DP) :: Series(3,nPairDensityTerms), value, slope
    INTEGER :: j

    CALL InterpolatePairScale &
           ( HelmTable % PairTemperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable % PairLogA0, &
             HelmTable % PairDLogA0dT, HelmTable % PairD2LogA0dT2, Temperature, jT, &
             A0, LogA0, dLogA0dT, d2LogA0dT2 )

    IF( Density <= 0.0_DP .OR. A0 <= 0.0_DP )THEN
      SignedYL = 0.0_DP
      dYLdEta  = 0.0_DP
      RETURN
    END IF

    Scale = rmu * A0 / Density
    IF( ABS(Eta) <= 1.0_DP )THEN
      CALL InterpolatePairDensitySeries(HelmTable,Temperature,Series)
      value = Series(1,nPairDensityTerms)
      slope = REAL(2*nPairDensityTerms-1,DP)*value
      DO j = nPairDensityTerms-1, 1, -1
        value = value*Eta*Eta+Series(1,j)
        slope = slope*Eta*Eta+REAL(2*j-1,DP)*Series(1,j)
      END DO
      SignedYL = Scale*value*Eta
      dYLdEta = Scale*slope
    ELSE
      SignedYL = Scale*SINH(Eta)
      dYLdEta = Scale*COSH(Eta)
    END IF

  END SUBROUTINE ComputePairLeptonFraction

  SUBROUTINE EvaluatePairSymmetricExtension &
    ( HelmTable, Temperature, Density, SignedYL, LeptonGasState, &
      FreeEnergy_Option, SpecificDerivatives_Option, Classical_Option )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    TYPE(HelmTableType), INTENT(IN)    :: HelmTable
    REAL(DP), INTENT(IN)               :: Temperature, Density, SignedYL
    TYPE(LeptonGasType), INTENT(INOUT) :: LeptonGasState
    TYPE(LeptonFreeEnergyType), INTENT(OUT), OPTIONAL :: FreeEnergy_Option
    REAL(DP), INTENT(OUT), OPTIONAL :: SpecificDerivatives_Option(9)
    LOGICAL, INTENT(IN), OPTIONAL :: Classical_Option

    INTEGER :: jT
    REAL(DP) :: A0, AMinusR, dLogA0dT, d2LogA0dT2, LogA0
    REAL(DP) :: dLogP0dT, d2LogP0dT2
    REAL(DP) :: F, F_T, F_TT, F_n, F_Tn
    REAL(DP) :: G, G_T, G_TT, eta, nNet, P0, P0_T, P0_TT
    REAL(DP) :: R, entropyUnit, LogP0, theta
    REAL(DP) :: d,ratioA,ratioN,delta,deltaT,deltaTT,pTT,K,K0,KT,KTT
    REAL(DP) :: LogZ2, LZ21, LZ22, LogD, LD1, LD2, LE1, LE2
    REAL(DP) :: CorrectionD, CorrectionE, CorrectionFn

    LeptonGasState%T = Temperature
    LeptonGasState%rho = Density
    LeptonGasState%yL = SignedYL

    CALL InterpolatePairScale &
           ( HelmTable % PairTemperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable % PairLogP0, &
             HelmTable % PairDLogP0dT, HelmTable % PairD2LogP0dT2, Temperature, jT, &
             P0, LogP0, dLogP0dT, d2LogP0dT2 )
    CALL InterpolatePairScale &
           ( HelmTable % PairTemperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable % PairLogA0, &
             HelmTable % PairDLogA0dT, HelmTable % PairD2LogA0dT2, Temperature, jT, &
             A0, LogA0, dLogA0dT, d2LogA0dT2 )

    IF( PRESENT(Classical_Option) )THEN
      IF( Classical_Option )THEN
        CALL ComputeClassicalPairScales &
               ( Temperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable%HbarC, LogA0, dLogA0dT, d2LogA0dT2 )
        A0 = EXP(LogA0)
        LogP0 = LogA0 + LOG(HelmTable%ThermalConstant*Temperature*ergmev)
        P0 = EXP(LogP0)
        dLogP0dT = dLogA0dT + 1.0_DP/Temperature
        d2LogP0dT2 = d2LogA0dT2 - 1.0_DP/(Temperature*Temperature)
      END IF
    END IF

    theta = HelmTable%ThermalConstant * Temperature
    nNet  = Density * SignedYL / rmu

    eta = PairDegeneracyFromScale(nNet,A0,LogA0)
    R = HYPOT(A0,nNet)

    IF( R > 0.0_DP )THEN
      AMinusR = -nNet * ( nNet / ( A0 + R ) )
    ELSE
      AMinusR = 0.0_DP
    END IF
    G       = nNet * eta + AMinusR
    G_T     = dLogA0dT * AMinusR
    IF( R > 0.0_DP )THEN
      G_TT = d2LogA0dT2 * AMinusR &
             - dLogA0dT*dLogA0dT * A0 * ( AMinusR / R )
    ELSE
      G_TT = 0.0_DP
    END IF

    P0_T  = P0 * dLogP0dT
    P0_TT = P0 * ( d2LogP0dT2 + dLogP0dT*dLogP0dT )

    F    = -P0 + ergmev * theta * G
    F_T  = -P0_T + ergmev * ( HelmTable%ThermalConstant * G + theta * G_T )
    F_TT = -P0_TT + ergmev * ( 2.0_DP * HelmTable%ThermalConstant * G_T &
                               + theta * G_TT )
    F_n  = ergmev * theta * eta
    IF( R > 0.0_DP )THEN
      F_Tn = ergmev * ( HelmTable%ThermalConstant * eta &
                        - theta * ( nNet / R ) * dLogA0dT )
    ELSE
      F_Tn = 0.0_DP
    END IF

    CorrectionD = 0.0_DP
    CorrectionE = 0.0_DP
    IF( PRESENT(Classical_Option) )THEN
      IF( Classical_Option )THEN
        ! First Fermi correction to the Legendre-transformed dilute gas:
        ! F_FD=F_MB+D+2D(n/A0)^2+O(f_max^2), D=theta*Z_2/4.
        ! Z_2=2 Z_1(T/2). Keep the density-dependent term in logarithms:
        ! D can underflow while D(n/A0)^2 remains representable.
        CALL ComputeClassicalPairScales &
               ( 0.5_DP*Temperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable%HbarC, LogZ2, LZ21, LZ22 )
        LogD = LOG(ergmev*theta/2.0_DP)+LogZ2
        LD1 = 1.0_DP/Temperature+0.5_DP*LZ21
        LD2 = -1.0_DP/(Temperature*Temperature)+0.25_DP*LZ22
        LE1 = LD1-2.0_DP*dLogA0dT
        LE2 = LD2-2.0_DP*d2LogA0dT2
        CorrectionD = EXP(LogD)
        CorrectionFn = 0.0_DP
        IF( nNet /= 0.0_DP )THEN
          CorrectionE = EXP(LogD+LOG(2.0_DP)+2.0_DP*(LOG(ABS(nNet))-LogA0))
          CorrectionFn = SIGN(EXP(LogD+LOG(4.0_DP)+LOG(ABS(nNet))-2.0_DP*LogA0),nNet)
        END IF
        F = F+CorrectionD+CorrectionE
        F_T = F_T+CorrectionD*LD1+CorrectionE*LE1
        F_TT = F_TT+CorrectionD*(LD1*LD1+LD2)+CorrectionE*(LE1*LE1+LE2)
        F_n = F_n+CorrectionFn
        F_Tn = F_Tn+CorrectionFn*LE1
      END IF
    END IF

    LeptonGasState % p    = P0 - ergmev * theta * AMinusR+CorrectionE-CorrectionD
    LeptonGasState % dpdT = nNet * F_Tn - F_T
    IF( R > 0.0_DP )THEN
      LeptonGasState % dpdr = ( ABS(SignedYL) / rmu ) * ergmev * theta &
                              * ABS( nNet / R )
    ELSE
      LeptonGasState % dpdr = 0.0_DP
    END IF
    LeptonGasState % dpdr = LeptonGasState%dpdr+2.0_DP*CorrectionE/Density

    LeptonGasState % e    = ( F - Temperature * F_T ) / Density
    LeptonGasState % dedT = -Temperature * F_TT / Density
    LeptonGasState % dedr = &
      ( nNet * ( F_n - Temperature * F_Tn ) &
        - ( F - Temperature * F_T ) ) / ( Density * Density )

    entropyUnit = kmev * ergmev / rmu
    LeptonGasState % s    = -F_T / Density / entropyUnit
    LeptonGasState % dsdT = -F_TT / Density / entropyUnit
    LeptonGasState % dsdr = &
      ( -nNet * F_Tn + F_T ) &
      / ( Density * Density ) / entropyUnit

    LeptonGasState % mu = F_n / ergmev

    IF( PRESENT(FreeEnergy_Option) )THEN
      FreeEnergy_Option%F = F
      FreeEnergy_Option%FT = F_T
      FreeEnergy_Option%FTT = F_TT
      FreeEnergy_Option%Fn = F_n
      FreeEnergy_Option%nFTn = nNet*F_Tn
      FreeEnergy_Option%nnFnn = 0.0_DP
      IF( R > 0.0_DP ) FreeEnergy_Option%nnFnn = ergmev*theta*nNet*(nNet/R)
      FreeEnergy_Option%nnFnn = FreeEnergy_Option%nnFnn+2.0_DP*CorrectionE
    END IF
    IF( PRESENT(SpecificDerivatives_Option) )THEN
      d=Density*ABS(SignedYL)*HelmTable%DensityScale
      ratioA=A0/R
      ratioN=nNet/R
      delta=-AMinusR
      deltaT=-dLogA0dT*ratioA*delta
      deltaTT=delta*(-d2LogA0dT2*ratioA &
        +dLogA0dT**2*ratioA*(ratioA-ratioN**2))
      pTT=P0_TT+ergmev*(2.0_DP*HelmTable%ThermalConstant*deltaT+theta*deltaTT)
      K0=nNet*ratioN
      K=ergmev*theta*K0
      KT=ergmev*K0*(HelmTable%ThermalConstant-theta*dLogA0dT*ratioA**2)
      KTT=ergmev*K0*(-2.0_DP*HelmTable%ThermalConstant*dLogA0dT*ratioA**2 &
        +theta*(-d2LogA0dT2*ratioA**2 &
                +dLogA0dT**2*(ratioA**4-2.0_DP*ratioA**2*ratioN**2)))
      SpecificDerivatives_Option(1)=F/d-HelmTable%lepton_mass*HelmTable%NumberPerMass*ergmev
      SpecificDerivatives_Option(2)=F_T/d
      SpecificDerivatives_Option(3)=F_TT/d
      SpecificDerivatives_Option(4)=LeptonGasState%p/(d*d)
      ! Match h5 packing: f,ft,ftt,fd,fdd,fdt,fddt,fdtt,fddtt.
      SpecificDerivatives_Option(5)=(K-2.0_DP*LeptonGasState%p)/(d*d*d)
      SpecificDerivatives_Option(6)=LeptonGasState%dpdT/(d*d)
      SpecificDerivatives_Option(7)=(KT-2.0_DP*LeptonGasState%dpdT)/(d*d*d)
      SpecificDerivatives_Option(8)=pTT/(d*d)
      SpecificDerivatives_Option(9)=(KTT-2.0_DP*pTT)/(d*d*d)
    END IF

  END SUBROUTINE EvaluatePairSymmetricExtension

  SUBROUTINE InterpolatePairDensitySeries &
    ( HelmTable, Temperature, Series )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(HelmTableType), INTENT(IN) :: HelmTable
    REAL(DP), INTENT(IN) :: Temperature
    REAL(DP), INTENT(OUT) :: Series(3,nPairDensityTerms)
    REAL(DP) :: h, x, y, f0, f1, m0, m1, q0, q1
    INTEGER :: jT, j, nT

    nT = SIZE(HelmTable%PairTemperature)
    jT = INT(LOG(Temperature/HelmTable%PairTemperature(1)) &
             /LOG(HelmTable%PairTemperature(nT)/HelmTable%PairTemperature(1))*REAL(nT-1,DP))+1
    jT = MAX(1,MIN(jT,nT-1))
    h = LOG(HelmTable%PairTemperature(jT+1)/HelmTable%PairTemperature(jT))
    x = MAX(0.0_DP,MIN(LOG(Temperature/HelmTable%PairTemperature(jT))/h,1.0_DP))
    y = 1.0_DP-x
    Series(:,1) = [1.0_DP,0.0_DP,0.0_DP]
    DO j = 2, nPairDensityTerms
      f0 = HelmTable%PairDensitySeries(1,j-1,jT)
      f1 = HelmTable%PairDensitySeries(1,j-1,jT+1)
      m0 = HelmTable%PairDensitySeries(2,j-1,jT)
      m1 = HelmTable%PairDensitySeries(2,j-1,jT+1)
      q0 = HelmTable%PairDensitySeries(3,j-1,jT)
      q1 = HelmTable%PairDensitySeries(3,j-1,jT+1)
      Series(1,j) = psi0(x)*f0+psi0(y)*f1+h*(psi1(x)*m0-psi1(y)*m1) &
                      +h*h*(psi2(x)*q0+psi2(y)*q1)
      Series(2,j) = (dpsi0(x)*f0-dpsi0(y)*f1)/h+dpsi1(x)*m0+dpsi1(y)*m1 &
                      +h*(dpsi2(x)*q0-dpsi2(y)*q1)
      Series(3,j) = (ddpsi0(x)*f0+ddpsi0(y)*f1)/(h*h) &
                      +(ddpsi1(x)*m0-ddpsi1(y)*m1)/h+ddpsi2(x)*q0+ddpsi2(y)*q1
    END DO

  END SUBROUTINE InterpolatePairDensitySeries

  SUBROUTINE EvaluateFermiPairExtension &
    ( HelmTable, Temperature, Density, SignedYL, State, Potential )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(HelmTableType), INTENT(IN) :: HelmTable
    REAL(DP), INTENT(IN) :: Temperature, Density, SignedYL
    TYPE(LeptonGasType), INTENT(OUT) :: State
    TYPE(LeptonFreeEnergyType), INTENT(OUT) :: Potential
    REAL(DP) :: Series(3,nPairDensityTerms), A0, LogA, LA1, LA2, P0, LogP, LP1, LP2
    REAL(DP) :: eta, eta2, ratio, value, slope, step, power, theta, nNet, etaT
    REAL(DP) :: pressure, pressureT, pressureTT, nT, LX, LXX, c, ct, ctt
    INTEGER :: jT, j, iteration

    CALL InterpolatePairScale &
           ( HelmTable%PairTemperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable%PairLogA0, &
             HelmTable%PairDLogA0dT, HelmTable%PairD2LogA0dT2, Temperature, jT, A0, LogA, LA1, LA2 )
    CALL InterpolatePairScale &
           ( HelmTable%PairTemperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable%PairLogP0, &
             HelmTable%PairDLogP0dT, HelmTable%PairD2LogP0dT2, Temperature, jT, P0, LogP, LP1, LP2 )
    IF( A0 <= TINY(1.0_DP) )THEN
      CALL EvaluatePairSymmetricExtension(HelmTable,Temperature,Density,SignedYL,State,Potential)
      RETURN
    END IF
    CALL InterpolatePairDensitySeries(HelmTable,Temperature,Series)
    nNet = Density*SignedYL/rmu
    ratio = nNet/A0
    eta = ASINH(ratio)
    ! Only the regular, bounded near-pair interval enters this inversion.
    ! This is not a nested thermodynamic solve: seven polynomial terms,
    ! with a positive slope and a few scalar Newton corrections.
    DO iteration = 1, 12
      eta2 = eta*eta
      value = Series(1,nPairDensityTerms)
      slope = REAL(2*nPairDensityTerms-1,DP)*value
      DO j = nPairDensityTerms-1, 1, -1
        value = value*eta2+Series(1,j)
        slope = slope*eta2+REAL(2*j-1,DP)*Series(1,j)
      END DO
      value = value*eta
      step = (value-ratio)/slope
      eta = eta-step
      IF( ABS(step) <= 8.0_DP*EPSILON(eta)*MAX(ABS(eta),TINY(eta)) ) EXIT
    END DO
    theta = HelmTable%ThermalConstant*Temperature
    LX = Temperature*LA1
    LXX = Temperature*Temperature*LA2+LX
    pressure = P0
    pressureT = P0*LP1
    pressureTT = P0*(LP1*LP1+LP2)
    nT = 0.0_DP
    power = eta
    DO j = 1, nPairDensityTerms
      c = ergmev*theta*A0/REAL(2*j,DP)
      ct = c*((1.0_DP+LX)*Series(1,j)+Series(2,j))/Temperature
      ctt = c*((LX*LX+LX+LXX)*Series(1,j) &
                 +(1.0_DP+2.0_DP*LX)*Series(2,j)+Series(3,j))/(Temperature*Temperature)
      pressure = pressure+c*Series(1,j)*power*eta
      pressureT = pressureT+ct*power*eta
      pressureTT = pressureTT+ctt*power*eta
      nT = nT+(LX*Series(1,j)+Series(2,j))*power
      power = power*eta*eta
    END DO
    etaT = -nT/(Temperature*slope)
    Potential%F = nNet*ergmev*theta*eta-pressure
    Potential%FT = nNet*ergmev*HelmTable%ThermalConstant*eta-pressureT
    Potential%FTT = -pressureTT+ergmev*theta*A0*slope*etaT*etaT
    Potential%Fn = ergmev*theta*eta
    Potential%nFTn = nNet*ergmev*HelmTable%ThermalConstant*(eta+Temperature*etaT)
    Potential%nnFnn = ergmev*theta*nNet*(ratio/slope)
    CALL SetLeptonStateFromFreeEnergy(Temperature,Density,SignedYL,Potential,State)

  END SUBROUTINE EvaluateFermiPairExtension

  SUBROUTINE SetLeptonStateFromFreeEnergy &
    (Temperature,Density,SignedYL,Potential,State)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: Temperature,Density,SignedYL
    TYPE(LeptonFreeEnergyType), INTENT(IN) :: Potential
    TYPE(LeptonGasType), INTENT(INOUT) :: State
    REAL(DP) :: nNet,EntropyUnit,InternalEnergy
    nNet = Density*SignedYL/rmu
    EntropyUnit = kmev*ergmev/rmu
    InternalEnergy = Potential%F-Temperature*Potential%FT
    State%T = Temperature
    State%rho = Density
    State%yL = SignedYL
    State%mu = Potential%Fn/ergmev
    State%p = nNet*Potential%Fn-Potential%F
    State%dpdT = Potential%nFTn-Potential%FT
    State%dpdr = Potential%nnFnn/Density
    State%e = InternalEnergy/Density
    State%dedT = -Temperature*Potential%FTT/Density
    State%dedr = (nNet*Potential%Fn-Temperature*Potential%nFTn &
      -InternalEnergy)/(Density*Density)
    State%s = -Potential%FT/(Density*EntropyUnit)
    State%dsdT = -Potential%FTT/(Density*EntropyUnit)
    State%dsdr = (Potential%FT-Potential%nFTn)/(Density*Density*EntropyUnit)
  END SUBROUTINE SetLeptonStateFromFreeEnergy


  SUBROUTINE ComputeDilutePairWeight &
    ( HelmTable, Temperature, Density, SignedYL, Weight )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(HelmTableType), INTENT(IN) :: HelmTable
    REAL(DP), INTENT(IN) :: Temperature, Density, SignedYL
    REAL(DP), INTENT(OUT) :: Weight(6)
    REAL(DP) :: a, LogA, LA1, LA2, A0, nNet, R, eta, rN, rZ
    REAL(DP) :: wa, wa1, wa2, wb, wb1, wb2, aT, aTT
    REAL(DP) :: bT, bTT, bZ, bZZ, bTZ, h, u1, u2

    Weight = 0.0_DP
    a = HelmTable%lepton_mass/(HelmTable%ThermalConstant*Temperature)
    IF( a <= 12.0_DP ) RETURN
    CALL ComputeClassicalPairScales &
           ( Temperature, HelmTable%lepton_mass, HelmTable%ThermalConstant, HelmTable%HbarC, LogA, LA1, LA2 )
    A0 = EXP(LogA)
    nNet = Density*SignedYL/rmu
    eta = ABS(PairDegeneracyFromScale(nNet,A0,LogA))
    h = HYPOT(1.0_DP,eta)
    IF( a-h <= 12.0_DP ) RETURN

    ! The maximum occupation is < exp(-12). Including the first Fermi
    ! correction leaves an O(exp(-24)) dilute error at the inner edge.
    ! z=log|n|. Store w,w_z,w_zz,w_T,w_TT,w_Tz.
    CALL SmoothTransition(a-12.0_DP,4.0_DP,wa,wa1,wa2)
    ! A smooth even weight is essential at n=0. A weight depending on
    ! |eta| would introduce a cusp into F and a spurious finite mu(0+).
    CALL SmoothTransition(a-h-12.0_DP,4.0_DP,wb,wb1,wb2)
    R = HYPOT(A0,nNet)
    rN = 0.0_DP
    rZ = 0.0_DP
    IF( R > 0.0_DP )THEN
      rN = ABS(nNet)/R
      rZ = rN*(A0/R)**2
    END IF
    aT = -a/Temperature
    aTT = 2.0_DP*a/(Temperature*Temperature)
    u1 = eta/h
    u2 = 1.0_DP/(h*h*h)
    bZ = -u1*rN
    bZZ = -u2*rN*rN-u1*rZ
    bT = u1*rN*LA1-a/Temperature
    bTT = aTT-(u2*rN*rN+u1*rZ)*LA1*LA1+u1*rN*LA2
    bTZ = (u2*rN*rN+u1*rZ)*LA1
    Weight(1) = wa*wb
    Weight(2) = wa*wb1*bZ
    Weight(3) = wa*(wb2*bZ*bZ+wb1*bZZ)
    Weight(4) = wa1*aT*wb+wa*wb1*bT
    Weight(5) = (wa2*aT*aT+wa1*aTT)*wb &
                  +2.0_DP*wa1*aT*wb1*bT+wa*(wb2*bT*bT+wb1*bTT)
    Weight(6) = wa1*aT*wb1*bZ+wa*(wb2*bT*bZ+wb1*bTZ)

  END SUBROUTINE ComputeDilutePairWeight

  SUBROUTINE BlendLeptonFreeEnergy &
    ( nNet, Weight, Reference, Potential )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: nNet, Weight(6)
    TYPE(LeptonFreeEnergyType), INTENT(IN) :: Reference
    TYPE(LeptonFreeEnergyType), INTENT(INOUT) :: Potential
    TYPE(LeptonFreeEnergyType) :: Difference
    REAL(DP) :: w, wz, wzz, wt, wtt, wtz

    w = Weight(1)
    wz = Weight(2)
    wzz = Weight(3)
    wt = Weight(4)
    wtt = Weight(5)
    wtz = Weight(6)
    Difference%F = Reference%F-Potential%F
    Difference%FT = Reference%FT-Potential%FT
    Difference%FTT = Reference%FTT-Potential%FTT
    Difference%Fn = Reference%Fn-Potential%Fn
    Difference%nFTn = Reference%nFTn-Potential%nFTn
    Difference%nnFnn = Reference%nnFnn-Potential%nnFnn
    Potential%F = Potential%F+w*Difference%F
    Potential%FT = Potential%FT+w*Difference%FT+wt*Difference%F
    Potential%FTT = Potential%FTT+w*Difference%FTT &
                     +2.0_DP*wt*Difference%FT+wtt*Difference%F
    Potential%Fn = Potential%Fn+w*Difference%Fn
    IF( nNet /= 0.0_DP ) Potential%Fn = Potential%Fn+(wz/nNet)*Difference%F
    Potential%nFTn = Potential%nFTn+w*Difference%nFTn &
                      +wz*Difference%FT+wt*(nNet*Difference%Fn)+wtz*Difference%F
    Potential%nnFnn = Potential%nnFnn+w*Difference%nnFnn &
                      +2.0_DP*wz*(nNet*Difference%Fn)+(wzz-wz)*Difference%F

  END SUBROUTINE BlendLeptonFreeEnergy

  SUBROUTINE SmoothTransition(Value,Width,W,W1,W2)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: Value,Width
    REAL(DP), INTENT(OUT) :: W,W1,W2
    REAL(DP) :: u
    u = MAX(0.0_DP,MIN(Value/Width,1.0_DP))
    W = u**3*(10.0_DP+u*(-15.0_DP+6.0_DP*u))
    W1 = 30.0_DP*u*u*(1.0_DP-u)**2/Width
    W2 = 60.0_DP*u*(1.0_DP-u)*(1.0_DP-2.0_DP*u)/(Width*Width)
  END SUBROUTINE SmoothTransition


  SUBROUTINE BlendPairSymmetricFreeEnergy &
    (HelmTable,Temperature,Density,SignedYL,Eta,TablePotential)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    TYPE(HelmTableType), INTENT(IN) :: HelmTable
    REAL(DP), INTENT(IN) :: Temperature,Density,SignedYL,Eta
    TYPE(LeptonFreeEnergyType), INTENT(INOUT) :: TablePotential
    TYPE(LeptonFreeEnergyType) :: Analytic,Difference
    TYPE(LeptonGasType) :: AnalyticState
    INTEGER :: jT
    REAL(DP) :: A0,LogA0,LA1,LA2,nNet,R,rN,rZ
    REAL(DP) :: wa,wa1,wa2,wd,wd1,wd2,waZ,waZZ,waT,waTT,waTZ
    REAL(DP) :: w,wZ,wZZ,wT,wTT,wTZ

    CALL EvaluateFermiPairExtension &
      (HelmTable,Temperature,Density,SignedYL,AnalyticState,Analytic)
    CALL InterpolatePairScale &
      (HelmTable%PairTemperature,HelmTable%lepton_mass,HelmTable%ThermalConstant,HelmTable%PairLogA0, &
       HelmTable%PairDLogA0dT, HelmTable % PairD2LogA0dT2,Temperature,jT,A0,LogA0,LA1,LA2)
    nNet = Density*SignedYL/rmu
    R = HYPOT(A0,nNet)
    IF( R <= 0.0_DP .OR. nNet == 0.0_DP )THEN
      TablePotential = Analytic
      RETURN
    END IF

    ! z=log|n|. Work with derivatives in z rather than 1/n**2.
    ! rN=d|eta|/dz; rZ=d rN/dz. Both are bounded.
    rN = ABS(nNet)/R
    rZ = rN*(A0/R)**2
    CALL SmoothTransition(ABS(Eta)-PairSymmetricEtaMin, &
      PairSymmetricEtaMax-PairSymmetricEtaMin,wa,wa1,wa2)
    waZ = wa1*rN
    waZZ = wa2*rN*rN+wa1*rZ
    waT = -wa1*rN*LA1
    waTT = wa2*(rN*LA1)**2+wa1*(rZ*LA1*LA1-rN*LA2)
    waTZ = -(wa2*rN*rN+wa1*rZ)*LA1
    CALL SmoothTransition(LOG(Density*ABS(SignedYL)*HelmTable%DensityScale/HelmTable%d(1)), &
      LOG(HelmTable%d(2)/HelmTable%d(1)),wd,wd1,wd2)
    w = wa*wd
    wZ = waZ*wd+wa*wd1
    wZZ = waZZ*wd+2.0_DP*waZ*wd1+wa*wd2
    wT = waT*wd
    wTT = waTT*wd
    wTZ = waTZ*wd+waT*wd1

    Difference%F = TablePotential%F-Analytic%F
    Difference%FT = TablePotential%FT-Analytic%FT
    Difference%FTT = TablePotential%FTT-Analytic%FTT
    Difference%Fn = TablePotential%Fn-Analytic%Fn
    Difference%nFTn = TablePotential%nFTn-Analytic%nFTn
    Difference%nnFnn = TablePotential%nnFnn-Analytic%nnFnn
    ! Differentiate F=Fa+w*(Ft-Fa), including every weight derivative.
    TablePotential%F = Analytic%F+w*Difference%F
    TablePotential%FT = Analytic%FT+w*Difference%FT+wT*Difference%F
    TablePotential%FTT = Analytic%FTT+w*Difference%FTT &
      +2.0_DP*wT*Difference%FT+wTT*Difference%F
    TablePotential%Fn = Analytic%Fn+w*Difference%Fn+(wZ/nNet)*Difference%F
    TablePotential%nFTn = Analytic%nFTn+w*Difference%nFTn &
      +wZ*Difference%FT+wT*(nNet*Difference%Fn)+wTZ*Difference%F
    TablePotential%nnFnn = Analytic%nnFnn+w*Difference%nnFnn &
      +2.0_DP*wZ*(nNet*Difference%Fn)+(wZZ-wZ)*Difference%F
  END SUBROUTINE BlendPairSymmetricFreeEnergy

  SUBROUTINE InterpolatePairScale &
    (TemperatureGrid,LeptonMass,ThermalConstant,LogScale,dLogScaledT,d2LogScaledT2,Temperature,jT_Option, &
     Scale,LogScaleValue,dLogdT,d2LogdT2)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: TemperatureGrid(:),LogScale(:),dLogScaledT(:),d2LogScaledT2(:)
    REAL(DP), INTENT(IN) :: Temperature, LeptonMass, ThermalConstant
    INTEGER, INTENT(OUT) :: jT_Option
    REAL(DP), INTENT(OUT) :: Scale,LogScaleValue,dLogdT,d2LogdT2
    INTEGER :: jT,nT
    REAL(DP) :: dT,x,y,f0,f1,m0,m1,q0,q1,value,a,a0,a1,hx,hxx
    nT=SIZE(TemperatureGrid)
    jT=INT((LOG10(Temperature)-LOG10(TemperatureGrid(1))) &
      /(LOG10(TemperatureGrid(nT))-LOG10(TemperatureGrid(1)))*REAL(nT-1,DP))+1
    jT=MAX(1,MIN(jT,nT-1))
    ! Interpolate in log T after removing the known Boltzmann exponent.
    ! In particular, do not interpolate -m/(kT) or differentiate it by
    ! differences: at the cold edge that destroyed entropy and heat capacity.
    dT=LOG(TemperatureGrid(jT+1)/TemperatureGrid(jT))
    x=MAX(0.0_DP,MIN(LOG(Temperature/TemperatureGrid(jT))/dT,1.0_DP))
    y=1.0_DP-x
    a=LeptonMass/(ThermalConstant*Temperature)
    a0=LeptonMass/(ThermalConstant*TemperatureGrid(jT))
    a1=LeptonMass/(ThermalConstant*TemperatureGrid(jT+1))
    f0=LogScale(jT)+a0
    f1=LogScale(jT+1)+a1
    m0=TemperatureGrid(jT)*dLogScaledT(jT)-a0
    m1=TemperatureGrid(jT+1)*dLogScaledT(jT+1)-a1
    q0=TemperatureGrid(jT)**2*d2LogScaledT2(jT) &
         +TemperatureGrid(jT)*dLogScaledT(jT)+a0
    q1=TemperatureGrid(jT+1)**2*d2LogScaledT2(jT+1) &
         +TemperatureGrid(jT+1)*dLogScaledT(jT+1)+a1
    value=psi0(x)*f0+psi0(y)*f1+dT*(psi1(x)*m0-psi1(y)*m1) &
      +dT*dT*(psi2(x)*q0+psi2(y)*q1)
    hx=(dpsi0(x)*f0-dpsi0(y)*f1)/dT+dpsi1(x)*m0+dpsi1(y)*m1 &
      +dT*(dpsi2(x)*q0-dpsi2(y)*q1)
    hxx=(ddpsi0(x)*f0+ddpsi0(y)*f1)/(dT*dT) &
      +(ddpsi1(x)*m0-ddpsi1(y)*m1)/dT+ddpsi2(x)*q0+ddpsi2(y)*q1
    dLogdT=(hx+a)/Temperature
    d2LogdT2=(hxx-hx-2.0_DP*a)/(Temperature*Temperature)
    LogScaleValue=value-a
    Scale=EXP(LogScaleValue)
    jT_Option=jT
  END SUBROUTINE InterpolatePairScale

  ! These are all the functions used fro the interpolation
  ! provided by Timmes
  ! quintic hermite polynomial functions
  ! psi0 and its derivatives
  REAL(DP) FUNCTION psi0(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    psi0 = z**3 * ( z * (-6.0_dp*z + 15.0_dp) -10.0_dp) + 1.0_dp
  END FUNCTION psi0

  REAL(DP) FUNCTION dpsi0(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    dpsi0 = z**2 * ( z * (-30.0_dp*z + 60.0_dp) - 30.0_dp)
  END FUNCTION dpsi0

  REAL(DP) FUNCTION ddpsi0(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    ddpsi0 = z* ( z*( -120.0_dp*z + 180.0_dp) -60.0_dp)
  END FUNCTION ddpsi0

  ! psi1 and its derivatives
  REAL(DP) FUNCTION psi1(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    psi1 = z* ( z**2 * ( z * (-3.0_dp*z + 8.0_dp) - 6.0_dp) + 1.0_dp)
  END FUNCTION psi1

  REAL(DP) FUNCTION dpsi1(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    dpsi1 = z*z * ( z * (-15.0_dp*z + 32.0_dp) - 18.0_dp) +1.0_dp
  END FUNCTION dpsi1

  REAL(DP) FUNCTION ddpsi1(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    ddpsi1 = z * (z * (-60.0_dp*z + 96.0_dp) -36.0_dp)
  END FUNCTION ddpsi1

  ! psi2  and its derivatives
  REAL(DP) FUNCTION psi2(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    psi2 = 0.5d0*z*z*( z* ( z * (-z + 3.0_dp) - 3.0_dp) + 1.0_dp)
  END FUNCTION psi2

  REAL(DP) FUNCTION dpsi2(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    dpsi2 = 0.5d0*z*( z*(z*(-5.0_dp*z + 12.0_dp) - 9.0_dp) + 2.0_dp)
  END FUNCTION dpsi2

  REAL(DP) FUNCTION ddpsi2(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    ddpsi2 = 0.5d0*(z*( z * (-20.0_dp*z + 36.0_dp) - 18.0_dp) + 2.0_dp)
  END FUNCTION ddpsi2


  ! biquintic hermite polynomial FUNCTION
  REAL(DP) FUNCTION h5(fi,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: fi(36)
    REAL(DP), INTENT(IN) :: w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md

    h5 =  fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
    + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
    + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
    + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
    + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
    + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
    + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
    + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
    + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
    + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
    + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
    + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
    + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
    + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
    + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
    + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
    + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
    + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt

  END FUNCTION h5

  ! cubic hermite polynomial FUNCTIONs
  ! psi0 & derivatives
  REAL(DP) FUNCTION xpsi0(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    xpsi0 = z * z * (2.0_dp*z - 3.0_dp) + 1.0
  END FUNCTION xpsi0

  REAL(DP) FUNCTION xdpsi0(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    xdpsi0 = z * (6.0_dp*z - 6.0_dp)
  END FUNCTION xdpsi0


  ! psi1 & derivatives
  REAL(DP) FUNCTION xpsi1(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    xpsi1 = z * ( z * (z - 2.0_dp) + 1.0_dp)
  END FUNCTION xpsi1

  REAL(DP) FUNCTION xdpsi1(z)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: z

    xdpsi1 = z * (3.0_dp*z - 4.0_dp) + 1.0_dp
  END FUNCTION xdpsi1

  ! bicubic hermite polynomial FUNCTION
  REAL(DP) FUNCTION h3(dfi,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: dfi(16)
    REAL(DP), INTENT(IN) :: w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md

    h3 =  dfi(1)  *w0d*w0t   +  dfi(2)  *w0md*w0t &
    + dfi(3)  *w0d*w0mt  +  dfi(4)  *w0md*w0mt &
    + dfi(5)  *w0d*w1t   +  dfi(6)  *w0md*w1t &
    + dfi(7)  *w0d*w1mt  +  dfi(8)  *w0md*w1mt &
    + dfi(9)  *w1d*w0t   +  dfi(10) *w1md*w0t &
    + dfi(11) *w1d*w0mt  +  dfi(12) *w1md*w0mt &
    + dfi(13) *w1d*w1t   +  dfi(14) *w1md*w1t &
    + dfi(15) *w1d*w1mt  +  dfi(16) *w1md*w1mt
  END FUNCTION h3

END MODULE wlLeptonPhotonGasEOS

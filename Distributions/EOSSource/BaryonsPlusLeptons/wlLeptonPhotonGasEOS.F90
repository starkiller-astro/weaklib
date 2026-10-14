MODULE wlLeptonPhotonGasEOS

  USE wlKindModule, ONLY: dp
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlHelmholtzEOS
  USE wlEosConstantsModule, ONLY: &
    pi, rmu, kerg, cvel, &
    kmev, ergmev, sigma_sb, asol

  IMPLICIT NONE
  PRIVATE

  ! The pair-symmetric representation agrees with the full Fermi gas through
  ! linear order in the net lepton density.  Use it directly for small eta
  ! and blend its Helmholtz free energy into the table over a finite interval.
  ! At the outer edge the first omitted relative correction is O(eta^2).
  REAL(DP), PARAMETER, PUBLIC :: PairSymmetricEtaMin = 5.0e-3_DP
  REAL(DP), PARAMETER, PUBLIC :: PairSymmetricEtaMax = 2.0e-2_DP

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
    INTEGER  :: iat, jat
    REAL(DP) :: free, df_d, df_t, df_tt, df_dt
    REAL(DP) :: xt, xd, mxt, mxd, fi(36),   &
    si0t, si1t, si2t, si0mt, si1mt, si2mt,  &
    si0d, si1d, si2d, si0md, si1md, si2md,  &
    dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt,  &
    dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md,  &
    ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt

    ! extra variables I added that were previously allocated
    REAL(DP) :: smallt, smalld, hight, highd
    REAL(DP) :: tstp, dstp, dstpi, tstpi
    TYPE(HelmholtzStateType) :: HelmState
    REAL(DP) :: etaPair, sign_yL

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

    IF ( yL*den > HelmTable % maxdens .OR. &
         temp   > HelmTable % maxtemp ) THEN

      STOP 'LeptonGasEOS: density or temperature is too high for the table'
      RETURN

    ENDIF

    ! The logarithmic net-density coordinate is ill-conditioned near pair
    ! symmetry: the chemical potential is recovered by cancelling its rest
    ! mass against the tabulated degeneracy parameter.  Use the regular
    ! pair-symmetric representation throughout the nondegenerate interval,
    ! not merely below the first net-density table point.
    CALL ComputePairDegeneracy &
           ( HelmTable, temp, den, sign_yL * yL, etaPair )
    IF ( yL*den <= HelmTable % mindens .OR. &
         ABS(etaPair) <= PairSymmetricEtaMin ) THEN
      CALL EvaluatePairSymmetricExtension &
             ( HelmTable, temp, den, sign_yL * yL, LeptonGasState )
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

    din   = yL * den

    !..initialize
    deni    = 1.0_dp/den
    tempi   = 1.0_dp/temp

    !..ion section: removed in this

    !..enter the table with yL*den
    din = yL*den

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

    !..now get the pressure derivative with density, chemical potential, and
    !..lepton positron number densities
    !..get the interpolation weight functions
    si0t   =  xpsi0(xt )
    si1t   =  xpsi1(xt ) * HelmTable % dt(jat)

    si0mt  =  xpsi0(mxt)
    si1mt  = -xpsi1(mxt) * HelmTable % dt(jat)

    si0d   =  xpsi0(xd )
    si1d   =  xpsi1(xd ) * HelmTable % dd(iat)

    si0md  =  xpsi0(mxd)
    si1md  = -xpsi1(mxd) * HelmTable % dd(iat)

    !..derivatives of weight functions
    dsi0t  = xdpsi0(xt ) * HelmTable % dti(jat)
    dsi1t  = xdpsi1(xt )

    dsi0mt = -xdpsi0(mxt) * HelmTable % dti(jat)
    dsi1mt =  xdpsi1(mxt)

    dsi0d  =  xdpsi0(xd ) * HelmTable % ddi(iat)
    dsi1d  =  xdpsi1(xd )

    dsi0md = -xdpsi0(mxd) * HelmTable % ddi(iat)
    dsi1md =  xdpsi1(mxd)

    !..look in the pressure derivative only once
    fi(1)  = HelmTable % dpdf  (iat,jat)
    fi(2)  = HelmTable % dpdf  (iat+1,jat)
    fi(3)  = HelmTable % dpdf  (iat,jat+1)
    fi(4)  = HelmTable % dpdf  (iat+1,jat+1)
    fi(5)  = HelmTable % dpdft (iat,jat)
    fi(6)  = HelmTable % dpdft (iat+1,jat)
    fi(7)  = HelmTable % dpdft (iat,jat+1)
    fi(8)  = HelmTable % dpdft (iat+1,jat+1)
    fi(9)  = HelmTable % dpdfd (iat,jat)
    fi(10) = HelmTable % dpdfd (iat+1,jat)
    fi(11) = HelmTable % dpdfd (iat,jat+1)
    fi(12) = HelmTable % dpdfd (iat+1,jat+1)
    fi(13) = HelmTable % dpdfdt(iat,jat)
    fi(14) = HelmTable % dpdfdt(iat+1,jat)
    fi(15) = HelmTable % dpdfdt(iat,jat+1)
    fi(16) = HelmTable % dpdfdt(iat+1,jat+1)

    !..pressure derivative with density
    dplep_dd  = h3(   fi, &
        si0t,   si1t,   si0mt,   si1mt, &
        si0d,   si1d,   si0md,   si1md)
    dplep_dd  = MAX( yL * dplep_dd, 0.0_dp )

    !..look in the lepton chemical potential table only once
    fi(1)  = HelmTable % ef  (iat,jat)
    fi(2)  = HelmTable % ef  (iat+1,jat)
    fi(3)  = HelmTable % ef  (iat,jat+1)
    fi(4)  = HelmTable % ef  (iat+1,jat+1)
    fi(5)  = HelmTable % eft (iat,jat)
    fi(6)  = HelmTable % eft (iat+1,jat)
    fi(7)  = HelmTable % eft (iat,jat+1)
    fi(8)  = HelmTable % eft (iat+1,jat+1)
    fi(9)  = HelmTable % efd (iat,jat)
    fi(10) = HelmTable % efd (iat+1,jat)
    fi(11) = HelmTable % efd (iat,jat+1)
    fi(12) = HelmTable % efd (iat+1,jat+1)
    fi(13) = HelmTable % efdt(iat,jat)
    fi(14) = HelmTable % efdt(iat+1,jat)
    fi(15) = HelmTable % efdt(iat,jat+1)
    fi(16) = HelmTable % efdt(iat+1,jat+1)

    !..lepton chemical potential etalep
    etalep  = h3( fi, &
        si0t,   si1t,   si0mt,   si1mt, &
        si0d,   si1d,   si0md,   si1md)

    !..look in the number density table only once
    fi(1)  = HelmTable % xf  (iat,jat)
    fi(2)  = HelmTable % xf  (iat+1,jat)
    fi(3)  = HelmTable % xf  (iat,jat+1)
    fi(4)  = HelmTable % xf  (iat+1,jat+1)
    fi(5)  = HelmTable % xft (iat,jat)
    fi(6)  = HelmTable % xft (iat+1,jat)
    fi(7)  = HelmTable % xft (iat,jat+1)
    fi(8)  = HelmTable % xft (iat+1,jat+1)
    fi(9)  = HelmTable % xfd (iat,jat)
    fi(10) = HelmTable % xfd (iat+1,jat)
    fi(11) = HelmTable % xfd (iat,jat+1)
    fi(12) = HelmTable % xfd (iat+1,jat+1)
    fi(13) = HelmTable % xfdt(iat,jat)
    fi(14) = HelmTable % xfdt(iat+1,jat)
    fi(15) = HelmTable % xfdt(iat,jat+1)
    fi(16) = HelmTable % xfdt(iat+1,jat+1)

    !..dplep_dd at high temperatures and low densities is below the
    !..floating point limit of the subtraction of two large terms.
    !..since dplep_dd doesn't enter the maxwell relations at all, use the
    !..bicubic interpolation done above instead of this one
    x        = din * din
    plep     = x   * df_d
    dplep_dt = x   * df_dt

    x        =  yL    * yL
    slep     = -df_t  * yL
    dslep_dt = -df_tt * yL
    dslep_dd = -df_dt * x

    elep     = yL*free + temp * slep
    delep_dt = temp * dslep_dt
    delep_dd = x * df_d + temp * dslep_dd

    !..maxwell relations; each is 0.0_dp IF the consistency is perfect
    ! x   = den * den
    ! dse = temp*dslep_dt/delep_dt - 1.0_dp
    ! dpe = (delep_dd*x + temp*dplep_dt)/pres - 1.0_dp
    ! dsp = -dslep_dd*x/dplep_dt - 1.0_dp

    LeptonGasState % p    = plep
    LeptonGasState % dpdT = dplep_dt
    LeptonGasState % dpdr = dplep_dd

    LeptonGasState % e    = elep + HelmTable % lepton_mass / rmu * ergmev * yL
    LeptonGasState % dedT = delep_dt
    LeptonGasState % dedr = delep_dd

    LeptonGasState % s    = slep     / (kmev * ergmev / rmu)
    LeptonGasState % dsdT = dslep_dt / (kmev * ergmev / rmu)
    LeptonGasState % dsdr = dslep_dd / (kmev * ergmev / rmu)

    ! do not forget to add back mass of the lepton!
    LeptonGasState % mu = etalep*temp*kmev + HelmTable % lepton_mass
    LeptonGasState % mu = sign_yL * LeptonGasState % mu

    IF( ABS(etaPair) < PairSymmetricEtaMax .OR. &
        yL*den < HelmTable % d(2) )THEN
      CALL BlendPairSymmetricThermodynamics &
             ( HelmTable, temp, den, sign_yL*yL, etaPair, &
               LeptonGasState )
    END IF

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
           ( HelmTable % PairTemperature, HelmTable % PairLogA0, &
             HelmTable % PairDLogA0dT, Temperature, jT, &
             A0, LogA0, dLogA0dT, d2LogA0dT2 )

    nNet = Density * SignedYL / rmu
    IF( A0 > 0.0_DP )THEN
      Eta = ASINH( nNet / A0 )
    ELSE IF( nNet == 0.0_DP )THEN
      Eta = 0.0_DP
    ELSE
      Eta = SIGN( LOG(2.0_DP*ABS(nNet)) - LogA0, nNet )
    END IF

  END SUBROUTINE ComputePairDegeneracy


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

    CALL InterpolatePairScale &
           ( HelmTable % PairTemperature, HelmTable % PairLogA0, &
             HelmTable % PairDLogA0dT, Temperature, jT, &
             A0, LogA0, dLogA0dT, d2LogA0dT2 )

    IF( Density <= 0.0_DP .OR. A0 <= 0.0_DP )THEN
      SignedYL = 0.0_DP
      dYLdEta  = 0.0_DP
      RETURN
    END IF

    Scale = rmu * A0 / Density
    SignedYL = Scale * SINH( Eta )
    dYLdEta  = Scale * COSH( Eta )

  END SUBROUTINE ComputePairLeptonFraction

  SUBROUTINE EvaluatePairSymmetricExtension &
    ( HelmTable, Temperature, Density, SignedYL, LeptonGasState )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    TYPE(HelmTableType), INTENT(IN)    :: HelmTable
    REAL(DP), INTENT(IN)               :: Temperature, Density, SignedYL
    TYPE(LeptonGasType), INTENT(INOUT) :: LeptonGasState

    INTEGER :: jT
    REAL(DP) :: A0, AMinusR, dLogA0dT, d2LogA0dT2, LogA0
    REAL(DP) :: dLogP0dT, d2LogP0dT2
    REAL(DP) :: F, F_T, F_TT, F_n, F_Tn
    REAL(DP) :: G, G_T, G_TT, eta, nNet, P0, P0_T, P0_TT
    REAL(DP) :: R, entropyUnit, LogP0, theta

    CALL InterpolatePairScale &
           ( HelmTable % PairTemperature, HelmTable % PairLogP0, &
             HelmTable % PairDLogP0dT, Temperature, jT, &
             P0, LogP0, dLogP0dT, d2LogP0dT2 )
    CALL InterpolatePairScale &
           ( HelmTable % PairTemperature, HelmTable % PairLogA0, &
             HelmTable % PairDLogA0dT, Temperature, jT, &
             A0, LogA0, dLogA0dT, d2LogA0dT2 )

    theta = kmev * Temperature
    nNet  = Density * SignedYL / rmu

    IF( A0 > 0.0_DP )THEN
      eta = ASINH( nNet / A0 )
      R   = HYPOT( A0, nNet )
    ELSE
      ! A0 can underflow only where the pair background is physically
      ! negligible.  Use the classical logarithmic limit without dividing by
      ! the underflowed value.
      IF( nNet == 0.0_DP )THEN
        eta = 0.0_DP
        R   = 0.0_DP
      ELSE
        eta = SIGN( LOG( 2.0_DP * ABS(nNet) ) - LogA0, nNet )
        R   = ABS(nNet)
      END IF
    END IF

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
    F_T  = -P0_T + ergmev * ( kmev * G + theta * G_T )
    F_TT = -P0_TT + ergmev * ( 2.0_DP * kmev * G_T &
                               + theta * G_TT )
    F_n  = ergmev * theta * eta
    IF( R > 0.0_DP )THEN
      F_Tn = ergmev * ( kmev * eta &
                        - theta * ( nNet / R ) * dLogA0dT )
    ELSE
      F_Tn = 0.0_DP
    END IF

    LeptonGasState % p    = P0 + ergmev * theta * ( R - A0 )
    LeptonGasState % dpdT = nNet * F_Tn - F_T
    IF( R > 0.0_DP )THEN
      LeptonGasState % dpdr = ABS(SignedYL) * ergmev * theta &
                              / rmu * ABS( nNet / R )
    ELSE
      LeptonGasState % dpdr = 0.0_DP
    END IF

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

    LeptonGasState % mu = theta * eta

  END SUBROUTINE EvaluatePairSymmetricExtension

  SUBROUTINE BlendPairSymmetricThermodynamics &
    ( HelmTable, Temperature, Density, SignedYL, Eta, TableState )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    TYPE(HelmTableType), INTENT(IN)    :: HelmTable
    REAL(DP), INTENT(IN)               :: Temperature, Density, SignedYL, Eta
    TYPE(LeptonGasType), INTENT(INOUT) :: TableState

    INTEGER :: jT
    REAL(DP) :: A0, dLogA0dT, d2LogA0dT2, LogA0
    REAL(DP) :: DeltaEta, EtaAbs, Eta_n, Eta_T, nNet, RatioN, R
    REAL(DP) :: u, w, w_n, w_rho, w_T
    REAL(DP) :: wa, waEta, wa_n, wa_T
    REAL(DP) :: wd, wdLog, wd_n
    REAL(DP) :: LogDensityWidth, NetMassDensity
    REAL(DP) :: Difference, OldValue
    TYPE(LeptonGasType) :: AnalyticState

    CALL EvaluatePairSymmetricExtension &
           ( HelmTable, Temperature, Density, SignedYL, AnalyticState )
    CALL InterpolatePairScale &
           ( HelmTable % PairTemperature, HelmTable % PairLogA0, &
             HelmTable % PairDLogA0dT, Temperature, jT, &
             A0, LogA0, dLogA0dT, d2LogA0dT2 )

    nNet = Density * SignedYL / rmu
    R = HYPOT( A0, nNet )
    IF( R <= 0.0_DP .OR. nNet == 0.0_DP )THEN
      TableState = AnalyticState
      RETURN
    END IF

    EtaAbs  = ABS(Eta)
    DeltaEta = PairSymmetricEtaMax - PairSymmetricEtaMin
    u = ( EtaAbs - PairSymmetricEtaMin ) / DeltaEta
    u = MAX( 0.0_DP, MIN( u, 1.0_DP ) )
    wa = u*u*u * ( 10.0_DP + u * ( -15.0_DP + 6.0_DP*u ) )
    waEta = 30.0_DP*u*u*(1.0_DP-u)*(1.0_DP-u) / DeltaEta

    RatioN = ABS(nNet) / R
    Eta_n  = SIGN( 1.0_DP, nNet ) / R
    Eta_T  = -RatioN * dLogA0dT

    wa_n  = waEta * Eta_n
    wa_T  = waEta * Eta_T

    NetMassDensity = Density * ABS(SignedYL)
    LogDensityWidth = LOG( HelmTable % d(2) / HelmTable % d(1) )
    u = LOG( NetMassDensity / HelmTable % d(1) ) / LogDensityWidth
    u = MAX( 0.0_DP, MIN( u, 1.0_DP ) )
    wd = u*u*u * ( 10.0_DP + u * ( -15.0_DP + 6.0_DP*u ) )
    wdLog = 30.0_DP*u*u*(1.0_DP-u)*(1.0_DP-u) / LogDensityWidth
    wd_n  = wdLog / nNet

    w    = wa * wd
    w_n  = wa_n*wd + wa*wd_n
    w_T  = wa_T * wd
    w_rho = w_n * nNet / Density

    OldValue = TableState % p
    Difference = OldValue - AnalyticState % p
    TableState % p = AnalyticState % p + w*Difference
    TableState % dpdT = AnalyticState % dpdT &
                         + w*(TableState % dpdT-AnalyticState % dpdT) &
                         + w_T*Difference
    TableState % dpdr = AnalyticState % dpdr &
                         + w*(TableState % dpdr-AnalyticState % dpdr) &
                         + w_rho*Difference

    OldValue = TableState % e
    Difference = OldValue - AnalyticState % e
    TableState % e = AnalyticState % e + w*Difference
    TableState % dedT = AnalyticState % dedT &
                         + w*(TableState % dedT-AnalyticState % dedT) &
                         + w_T*Difference
    TableState % dedr = AnalyticState % dedr &
                         + w*(TableState % dedr-AnalyticState % dedr) &
                         + w_rho*Difference

    OldValue = TableState % s
    Difference = OldValue - AnalyticState % s
    TableState % s = AnalyticState % s + w*Difference
    TableState % dsdT = AnalyticState % dsdT &
                         + w*(TableState % dsdT-AnalyticState % dsdT) &
                         + w_T*Difference
    TableState % dsdr = AnalyticState % dsdr &
                         + w*(TableState % dsdr-AnalyticState % dsdr) &
                         + w_rho*Difference

    TableState % mu = AnalyticState % mu &
                       + w*(TableState % mu-AnalyticState % mu)

  END SUBROUTINE BlendPairSymmetricThermodynamics

  SUBROUTINE InterpolatePairScale &
    ( TemperatureGrid, LogScale, dLogScaledT, Temperature, jT_Option, &
      Scale, LogScaleValue, dLogdT, d2LogdT2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), DIMENSION(:), INTENT(IN) :: TemperatureGrid
    REAL(DP), DIMENSION(:), INTENT(IN) :: LogScale, dLogScaledT
    REAL(DP), INTENT(IN) :: Temperature
    INTEGER, INTENT(OUT) :: jT_Option
    REAL(DP), INTENT(OUT) :: Scale, LogScaleValue, dLogdT, d2LogdT2

    INTEGER :: jT, nT
    REAL(DP) :: dT, f0, f1, h00, h01, h10, h11
    REAL(DP) :: m0, m1, value, x

    nT = SIZE(TemperatureGrid)
    jT = INT( ( LOG10(Temperature) - LOG10(TemperatureGrid(1)) ) &
              / ( LOG10(TemperatureGrid(nT)) &
                  - LOG10(TemperatureGrid(1)) ) * REAL(nT-1,DP) ) + 1
    jT = MAX( 1, MIN( jT, nT-1 ) )

    dT = TemperatureGrid(jT+1) - TemperatureGrid(jT)
    x  = ( Temperature - TemperatureGrid(jT) ) / dT
    x  = MAX( 0.0_DP, MIN( x, 1.0_DP ) )
    f0 = LogScale(jT)
    f1 = LogScale(jT+1)
    m0 = dLogScaledT(jT)
    m1 = dLogScaledT(jT+1)

    h00 =  2.0_DP*x**3 - 3.0_DP*x*x + 1.0_DP
    h10 =         x**3 - 2.0_DP*x*x + x
    h01 = -2.0_DP*x**3 + 3.0_DP*x*x
    h11 =         x**3 -        x*x
    value = h00*f0 + h10*dT*m0 + h01*f1 + h11*dT*m1

    dLogdT = ( ( 6.0_DP*x*x - 6.0_DP*x ) * f0 &
                + ( 3.0_DP*x*x - 4.0_DP*x + 1.0_DP ) * dT*m0 &
                + ( -6.0_DP*x*x + 6.0_DP*x ) * f1 &
                + ( 3.0_DP*x*x - 2.0_DP*x ) * dT*m1 ) / dT
    d2LogdT2 = ( ( 12.0_DP*x - 6.0_DP ) * f0 &
                   + ( 6.0_DP*x - 4.0_DP ) * dT*m0 &
                   + ( -12.0_DP*x + 6.0_DP ) * f1 &
                   + ( 6.0_DP*x - 2.0_DP ) * dT*m1 ) / ( dT*dT )

    Scale         = EXP(value)
    LogScaleValue = value
    jT_Option     = jT

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

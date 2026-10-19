PROGRAM wlTestAbEm2DVsGSI

  !=========================================================================
  !  Accuracy test for the semi-leptonic charged-current 2D kernels.
  !
  !  For every thermodynamic point in ThermoConditions.dat (up to
  !  nPointsMax), every neutrino energy, reaction (1-4) and correction
  !  level (0-3) it computes:
  !
  !    W : wlSemiLeptonicOpacityModule2D :: Opacity_CC_2D
  !    P : wlSemiLeptonicOpacityModule2D :: Opacity_CC_2D_PrecomputedNodes
  !    I : wlSemiLeptonicOpacityIntegrationModule2D wrappers
  !    G : original GSI reference code (Opacity_SemiLep_2D_GSI.f90)
  !
  !  and reports per (reaction, correction):
  !    - W vs P and W vs I (expected to agree bitwise)
  !    - W vs G (differences from GSI constants, Tfac, real(4) literals)
  !    - for reactions 3,4 also W vs |G| and the number of sign
  !      mismatches. The GSI code clips |M|^2 >= 0 BEFORE applying the
  !      crossing minus sign for 3->1 / 1->3 processes and then flips the
  !      sign of the result, so G <= 0 (and mostly 0) for reactions 3,4.
  !    - detailed balance between reaction 3 (inverse neutron decay
  !      absorption opacity, cm^-1) and reaction 4 (neutron decay
  !      emissivity, (s cm^3 MeV)^-1):
  !        j(E) = c * kappa_a(E) * E^2 / (2 pi^2 (hbar c)^3)
  !               * exp( -(E - mu_nubar)/T ),  mu_nubar = Mun - Mup - Mul
  !      This holds exactly for the integrands (same |M|^2, occupation
  !      factors related by mu_n = mu_p + mu_l + mu_nubar, and prefactor
  !      ratio c E^2/(2 pi^2 (hbar c)^3)), so the ratio should be one up
  !      to quadrature differences in the reaction-dependent bounds.
  !
  !  Units: W, P, I are in cm^-1 (reactions 1-3) and (s cm^3 MeV)^-1
  !  (reaction 4). G is in km^-1 for reactions 1-3 (converted here to
  !  cm^-1) and (s cm^3 MeV)^-1 for reaction 4 (no conversion).
  !
  !  The GSI code uses a hardcoded Ngrids = 50 Gauss-Legendre points per
  !  dimension, so nE = 50 is used for W, P and I as well.
  !=========================================================================

  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D, Opacity_CC_2D_PrecomputedNodes, gauleg
  USE wlSemiLeptonicOpacityIntegrationModule2D, ONLY: &
    NuAbsorptionOnNeutrons, NuBarAbsorptionOnProtons, &
    InverseNeutronDecay, NeutronEmissivity
  USE wlEosConstantsModule, ONLY: &
    pi, Gw_MeV, ga, mn, mp, me, mmu, &
    massA, massV, gamma_p, gamma_n, hbarc, cvel

  IMPLICIT NONE

  INTEGER , PARAMETER :: nPointsMax = 100   ! number of thermo points used
  INTEGER , PARAMETER :: nEnu = 10          ! number of neutrino energies
  INTEGER , PARAMETER :: nReac = 4, nCorr = 4
  INTEGER , PARAMETER :: nE = 50            ! = Ngrids hardcoded in GSI code
  REAL(DP), PARAMETER :: EnuMin = 2.0d0, EnuMax = 300.0d0
  REAL(DP), PARAMETER :: km_to_cm_inv = 1.0d-5
  ! W-vs-G "significant" samples: max(|W|,|G|) >= SigFrac * max_E |W| at the
  ! same (point, reaction, correction). Excludes the deep Boltzmann tails
  ! (1e-40 and below) where the two codes cut the phase space differently.
  REAL(DP), PARAMETER :: SigFrac = 1.0d-6
  INTEGER , PARAMETER :: nWorstDB = 5
  LOGICAL , PARAMETER :: DoMuons = .TRUE.

  ! GSI Constants, they will change the results!
  REAL(DP), PARAMETER :: pi_GSI = 3.1415927d0, Gw_MeV_GSI = 1.166d-11, &
      gA_GSI = 1.2723d0, Mnp_GSI = 938.919d0, &
      Dnp_GSI = 1.293d0, massA_GSI = 1.0d3, massV_GSI = 840.0d0, &
      F2wm0_GSI = 3.706d0

  REAL(DP) :: EnuA(nEnu), xa(nE), wxa(nE)
  REAL(DP) :: massl, mu_nubar, DBpref, DBratio
  INTEGER  :: IncludeCorrections(3)
  INTEGER  :: i, j, k, l, nThermoPoints, nPoints, nSamples, iS

  REAL(DP), ALLOCATABLE :: T(:), Rho(:), Ye(:), Ym(:), Mue(:), Mum(:), Mul(:), &
      Mun(:), Mup(:), Un(:), Up(:), EffMassn(:), EffMassp(:)
  REAL(DP), ALLOCATABLE :: Opa_W(:,:,:,:), Opa_P(:,:,:,:), Opa_I(:,:,:,:), Opa_G(:,:,:,:)
  REAL(DP), ALLOCATABLE :: DB(:,:,:)

  ! Statistics work arrays
  REAL(DP), ALLOCATABLE :: rWP(:), rWI(:), rWG(:), rWGs(:), rWAG(:), rDB(:)
  INTEGER  :: nWP, nWI, nWG, nWGs, nWAG, nDB
  REAL(DP) :: WmaxE, worst(nWorstDB)
  INTEGER  :: iworst(3,nWorstDB), m
  INTEGER  :: nNaN_W, nNaN_WI0, nExactWP, nExactWI, nSignMis, nOneZeroG, nOneZeroDB

  INTEGER  :: count_start, count_end, count_rate
  REAL(DP) :: t_W, t_P, t_I, t_G

  WRITE(*,*) 'CONSTANTS REL ERR (weaklib vs GSI)'
  WRITE(*,*) 'pi    ', ABS(pi - pi_GSI)/pi_GSI
  WRITE(*,*) 'Gw_MeV', ABS(Gw_MeV - Gw_MeV_GSI)/Gw_MeV_GSI
  WRITE(*,*) 'ga    ', ABS(ga - gA_GSI)/gA_GSI
  WRITE(*,*) 'Mnp   ', ABS(Mnp_GSI - (mn + mp)/2.0d0)/Mnp_GSI
  WRITE(*,*) 'Dnp   ', ABS(Dnp_GSI - (mn - mp))/Dnp_GSI
  WRITE(*,*) 'massA ', ABS(massA_GSI - massA)/massA_GSI
  WRITE(*,*) 'massV ', ABS(massV_GSI - massV)/massV_GSI
  WRITE(*,*) 'F2wm0 ', ABS(F2wm0_GSI - (gamma_p - gamma_n - 1.0d0))/F2wm0_GSI

  ! Log-spaced neutrino energies
  DO l = 1, nEnu
    EnuA(l) = EnuMin * (EnuMax/EnuMin)**( DBLE(l-1) / DBLE(nEnu-1) )
  END DO

  ! Read thermodynamic data
  ! Columns: T Rho Ye Ym Mue Mum Mun Mup Un Up M*n M*p
  ! (relativistic chemical potentials, MeV)
  OPEN(UNIT=123, FILE=TRIM(ADJUSTL('ThermoConditions.dat')), STATUS='OLD', ACTION='READ')
  READ(123,*) nThermoPoints
  READ(123,*)

  nPoints = MIN(nThermoPoints, nPointsMax)

  ALLOCATE(T(nPoints), Rho(nPoints), Ye(nPoints), Ym(nPoints))
  ALLOCATE(Mue(nPoints), Mum(nPoints), Mul(nPoints), Mun(nPoints), Mup(nPoints))
  ALLOCATE(Un(nPoints), Up(nPoints), EffMassn(nPoints), EffMassp(nPoints))

  DO i = 1, nPoints
    READ(123,*) T(i), Rho(i), Ye(i), Ym(i), Mue(i), Mum(i), &
        Mun(i), Mup(i), Un(i), Up(i), EffMassn(i), EffMassp(i)
  END DO
  CLOSE(123)

  IF (DoMuons) THEN
    Mul(:) = Mum(:)
    massl  = mmu
  ELSE
    Mul(:) = Mue(:)
    massl  = me
  END IF

  WRITE(*,'(/,A,I6,A,I4,A,I4)') ' Thermo points:', nPoints, '   energies:', nEnu, '   nE:', nE
  WRITE(*,'(A,L2)') ' DoMuons:', DoMuons

  ALLOCATE(Opa_W(nEnu, nPoints, nCorr, nReac))
  ALLOCATE(Opa_P(nEnu, nPoints, nCorr, nReac))
  ALLOCATE(Opa_I(nEnu, nPoints, nCorr, nReac))
  ALLOCATE(Opa_G(nEnu, nPoints, nCorr, nReac))
  ALLOCATE(DB(nEnu, nPoints, nCorr))

  Opa_W = 0.0d0
  Opa_P = 0.0d0
  Opa_I = 0.0d0
  Opa_G = 0.0d0
  DB    = 0.0d0

  CALL gauleg(0.d0, 1.d0, xa, wxa, nE)

  ! W: Opacity_CC_2D
  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
  DO k = 1, nReac
  DO j = 1, nCorr
  DO i = 1, nPoints
  DO l = 1, nEnu
    CALL Opacity_CC_2D(j-1, k, EnuA(l), Opa_W(l,i,j,k), T(i), Mul(i), Mun(i), Mup(i), &
          massl, EffMassn(i), EffMassp(i), Un(i), Up(i), nE)
  END DO
  END DO
  END DO
  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_W = REAL(count_end - count_start, DP) / REAL(count_rate, DP)

  ! P: Opacity_CC_2D_PrecomputedNodes
  CALL SYSTEM_CLOCK(count_start)
  DO k = 1, nReac
  DO j = 1, nCorr
  DO i = 1, nPoints
  DO l = 1, nEnu
    CALL Opacity_CC_2D_PrecomputedNodes(j-1, k, EnuA(l), Opa_P(l,i,j,k), T(i), Mul(i), &
          Mun(i), Mup(i), massl, EffMassn(i), EffMassp(i), Un(i), Up(i), nE, xa, wxa)
  END DO
  END DO
  END DO
  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_P = REAL(count_end - count_start, DP) / REAL(count_rate, DP)

  ! I: integration module wrappers. IncludeCorrections mapping is the
  ! same as in Opacity_CC_2D_PrecomputedNodes: (WM, PS, FF)
  CALL SYSTEM_CLOCK(count_start)
  DO j = 1, nCorr
    SELECT CASE(j-1)
    CASE (0)
      IncludeCorrections = (/ 0, 0, 0 /)
    CASE (1)
      IncludeCorrections = (/ 1, 0, 0 /)
    CASE (2)
      IncludeCorrections = (/ 1, 1, 0 /)
    CASE (3)
      IncludeCorrections = (/ 1, 1, 1 /)
    END SELECT
    DO i = 1, nPoints
    DO l = 1, nEnu
      CALL NuAbsorptionOnNeutrons(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
          EffMassn(i), EffMassp(i), Un(i), Up(i), &
          IncludeCorrections, nE, xa, wxa, xa, wxa, Opa_I(l,i,j,1))
      CALL NuBarAbsorptionOnProtons(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
          EffMassn(i), EffMassp(i), Un(i), Up(i), &
          IncludeCorrections, nE, xa, wxa, xa, wxa, Opa_I(l,i,j,2))
      CALL InverseNeutronDecay(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
          EffMassn(i), EffMassp(i), Un(i), Up(i), &
          IncludeCorrections, nE, xa, wxa, xa, wxa, Opa_I(l,i,j,3))
      CALL NeutronEmissivity(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
          EffMassn(i), EffMassp(i), Un(i), Up(i), &
          IncludeCorrections, nE, xa, wxa, xa, wxa, Opa_I(l,i,j,4))
    END DO
    END DO
  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_I = REAL(count_end - count_start, DP) / REAL(count_rate, DP)

  ! G: GSI reference (serial only, uses module globals)
  CALL SYSTEM_CLOCK(count_start)
  DO k = 1, nReac
  DO j = 1, nCorr
  DO i = 1, nPoints
    CALL Opacity_CC_2D_GSI(j-1, k, nEnu, EnuA, Opa_G(:,i,j,k), T(i), Mul(i), &
          Mun(i), Mup(i), massl, EffMassn(i), EffMassp(i), Un(i), Up(i))
    IF (k <= 3) Opa_G(:,i,j,k) = Opa_G(:,i,j,k) * km_to_cm_inv ! 1/km to 1/cm
  END DO
  END DO
  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_G = REAL(count_end - count_start, DP) / REAL(count_rate, DP)

  ! Detailed balance ratio j / ( c kappa_a E^2/(2 pi^2 (hbar c)^3) exp(-(E-mu)/T) )
  DO j = 1, nCorr
  DO i = 1, nPoints
    mu_nubar = Mun(i) - Mup(i) - Mul(i)
    DO l = 1, nEnu
      IF (Opa_W(l,i,j,3) > 0.0d0 .AND. Opa_W(l,i,j,4) > 0.0d0) THEN
        DBpref = cvel * EnuA(l)**2 / ( 2.0d0 * pi**2 * hbarc**3 )
        DB(l,i,j) = EXP( LOG(Opa_W(l,i,j,4)) - LOG(Opa_W(l,i,j,3)) - LOG(DBpref) &
                         + (EnuA(l) - mu_nubar) / T(i) )
      ELSE
        DB(l,i,j) = 0.0d0
      END IF
    END DO
  END DO
  END DO

  ! Per-sample output
  OPEN(UNIT=200, FILE='AbEm2DVsGSI.dat', STATUS='REPLACE', ACTION='WRITE')
  WRITE(200,'(A)') '# iPoint Reaction Correction Enu[MeV] T[MeV] W P I G DBratio(reac 4 only)'
  WRITE(200,'(A)') '# units: cm^-1 for reactions 1-3, (s cm^3 MeV)^-1 for reaction 4'
  DO k = 1, nReac
  DO j = 1, nCorr
  DO i = 1, nPoints
  DO l = 1, nEnu
    IF (k == 4) THEN
      DBratio = DB(l,i,j)
    ELSE
      DBratio = 0.0d0
    END IF
    WRITE(200,'(3I6,7ES17.7E3)') i, k, j-1, EnuA(l), T(i), &
        Opa_W(l,i,j,k), Opa_P(l,i,j,k), Opa_I(l,i,j,k), Opa_G(l,i,j,k), DBratio
  END DO
  END DO
  END DO
  END DO
  CLOSE(200)

  ! Summary
  nSamples = nEnu * nPoints
  ALLOCATE(rWP(nSamples), rWI(nSamples), rWG(nSamples), rWGs(nSamples), &
           rWAG(nSamples), rDB(nSamples))

  WRITE(*,'(/,A)') ' Relative differences |a-b|/max(|a|,|b|) over all (point, energy) samples'
  WRITE(*,'(A)')   ' (samples where both values are zero count as 0; NaNs are excluded and counted)'
  WRITE(*,'(A,I8)') ' Samples per (reaction, correction):', nSamples
  WRITE(*,'(A,ES8.1,A)') ' sig: only samples with max(|W|,|G|) >=', SigFrac, &
                         ' * max_E|W| at that point'
  WRITE(*,'(/,A)') ' Reac Corr |  W-P max   W-P med  exact | &
                   &W-I max   W-I med  exact NaN(W) NaN(W)&I=0 | &
                   &W-G max   W-G med  1zero  W-G max(sig) med(sig)  nSig | &
                   &W-|G| max  W-|G| med SignMis'

  DO k = 1, nReac
  DO j = 1, nCorr

    nWP = 0; nWI = 0; nWG = 0; nWGs = 0; nWAG = 0
    nNaN_W = 0; nNaN_WI0 = 0; nExactWP = 0; nExactWI = 0
    nSignMis = 0; nOneZeroG = 0

    DO i = 1, nPoints

    WmaxE = MAXVAL(ABS(Opa_W(:,i,j,k)), MASK = .NOT. (Opa_W(:,i,j,k) /= Opa_W(:,i,j,k)))

    DO l = 1, nEnu

      IF (IsNaN_DP(Opa_W(l,i,j,k))) THEN
        nNaN_W = nNaN_W + 1
        IF (Opa_I(l,i,j,k) == 0.0d0) nNaN_WI0 = nNaN_WI0 + 1
        CYCLE
      END IF

      IF (Opa_W(l,i,j,k) == Opa_P(l,i,j,k)) nExactWP = nExactWP + 1
      IF (Opa_W(l,i,j,k) == Opa_I(l,i,j,k)) nExactWI = nExactWI + 1

      IF (.NOT. IsNaN_DP(Opa_P(l,i,j,k))) THEN
        nWP = nWP + 1
        rWP(nWP) = RelDiff(Opa_W(l,i,j,k), Opa_P(l,i,j,k))
      END IF

      IF (.NOT. IsNaN_DP(Opa_I(l,i,j,k))) THEN
        nWI = nWI + 1
        rWI(nWI) = RelDiff(Opa_W(l,i,j,k), Opa_I(l,i,j,k))
      END IF

      IF (.NOT. IsNaN_DP(Opa_G(l,i,j,k))) THEN
        nWG = nWG + 1
        rWG(nWG) = RelDiff(Opa_W(l,i,j,k), Opa_G(l,i,j,k))
        IF ( (Opa_W(l,i,j,k) == 0.0d0) .NEQV. (Opa_G(l,i,j,k) == 0.0d0) ) &
          nOneZeroG = nOneZeroG + 1
        IF ( MAX(ABS(Opa_W(l,i,j,k)), ABS(Opa_G(l,i,j,k))) >= SigFrac * WmaxE &
             .AND. WmaxE > 0.0d0 ) THEN
          nWGs = nWGs + 1
          rWGs(nWGs) = rWG(nWG)
        END IF
        IF (k >= 3) THEN
          nWAG = nWAG + 1
          rWAG(nWAG) = RelDiff(Opa_W(l,i,j,k), ABS(Opa_G(l,i,j,k)))
          IF ( Opa_W(l,i,j,k) * Opa_G(l,i,j,k) < 0.0d0 .OR. &
               ( Opa_W(l,i,j,k) > 0.0d0 .AND. Opa_G(l,i,j,k) == 0.0d0 ) .OR. &
               ( Opa_W(l,i,j,k) == 0.0d0 .AND. Opa_G(l,i,j,k) /= 0.0d0 ) ) &
            nSignMis = nSignMis + 1
        END IF
      END IF

    END DO
    END DO

    IF (k >= 3) THEN
      WRITE(*,'(2I5,A,2ES10.2,I7,A,2ES10.2,I7,I7,I11,A,2ES10.2,I7,2ES13.2,I6,A,2ES11.2,I8)') &
        k, j-1, ' |', MaxOf(rWP, nWP), MedianOf(rWP, nWP), nExactWP, ' |', &
        MaxOf(rWI, nWI), MedianOf(rWI, nWI), nExactWI, nNaN_W, nNaN_WI0, ' |', &
        MaxOf(rWG, nWG), MedianOf(rWG, nWG), nOneZeroG, &
        MaxOf(rWGs, nWGs), MedianOf(rWGs, nWGs), nWGs, ' |', &
        MaxOf(rWAG, nWAG), MedianOf(rWAG, nWAG), nSignMis
    ELSE
      WRITE(*,'(2I5,A,2ES10.2,I7,A,2ES10.2,I7,I7,I11,A,2ES10.2,I7,2ES13.2,I6,A)') &
        k, j-1, ' |', MaxOf(rWP, nWP), MedianOf(rWP, nWP), nExactWP, ' |', &
        MaxOf(rWI, nWI), MedianOf(rWI, nWI), nExactWI, nNaN_W, nNaN_WI0, ' |', &
        MaxOf(rWG, nWG), MedianOf(rWG, nWG), nOneZeroG, &
        MaxOf(rWGs, nWGs), MedianOf(rWGs, nWGs), nWGs, ' |        -          -       -'
    END IF

  END DO
  END DO

  ! Detailed balance summary
  WRITE(*,'(/,A)') ' Detailed balance (reaction 4 vs reaction 3, W):'
  WRITE(*,'(A)')   '   ratio = j / ( c kappa_a E^2/(2 pi^2 (hbar c)^3) exp(-(E-mu_nubar)/T) ),'
  WRITE(*,'(A)')   '   mu_nubar = Mun - Mup - Mul; expected 1 (only samples with both > 0)'
  WRITE(*,'(A)')   ' Corr   nUsed  1zero   min ratio   max ratio   median |ratio-1|   max |ratio-1|'
  DO j = 1, nCorr
    nDB = 0; nOneZeroDB = 0
    DO i = 1, nPoints
    DO l = 1, nEnu
      IF (DB(l,i,j) > 0.0d0) THEN
        nDB = nDB + 1
        rDB(nDB) = DB(l,i,j)
      ELSE IF ( (Opa_W(l,i,j,3) > 0.0d0) .NEQV. (Opa_W(l,i,j,4) > 0.0d0) ) THEN
        nOneZeroDB = nOneZeroDB + 1
      END IF
    END DO
    END DO
    IF (nDB > 0) THEN
      WRITE(*,'(I5,2I7,2ES12.4,2ES18.4)') j-1, nDB, nOneZeroDB, &
        MINVAL(rDB(1:nDB)), MAXVAL(rDB(1:nDB)), &
        MedianOf(ABS(rDB(1:nDB) - 1.0d0), nDB), MAXVAL(ABS(rDB(1:nDB) - 1.0d0))
    ELSE
      WRITE(*,'(I5,2I7,A)') j-1, nDB, nOneZeroDB, '   (no samples)'
    END IF
  END DO

  ! Worst detailed-balance samples (all corrections)
  worst  = -1.0d0
  iworst = 0
  DO j = 1, nCorr
  DO i = 1, nPoints
  DO l = 1, nEnu
    IF (DB(l,i,j) > 0.0d0) THEN
      IF (ABS(DB(l,i,j) - 1.0d0) > worst(nWorstDB)) THEN
        m = nWorstDB
        DO WHILE (m > 1)
          IF (worst(m-1) >= ABS(DB(l,i,j) - 1.0d0)) EXIT
          worst(m) = worst(m-1); iworst(:,m) = iworst(:,m-1)
          m = m - 1
        END DO
        worst(m) = ABS(DB(l,i,j) - 1.0d0); iworst(:,m) = (/ l, i, j /)
      END IF
    END IF
  END DO
  END DO
  END DO
  WRITE(*,'(/,A)') ' Worst detailed-balance samples:'
  WRITE(*,'(A)')   ' iPoint Corr    Enu[MeV]     T[MeV]  mu_nubar      kappa_a(3)      j(4)       ratio'
  DO m = 1, nWorstDB
    IF (iworst(1,m) == 0) CYCLE
    l = iworst(1,m); i = iworst(2,m); j = iworst(3,m)
    WRITE(*,'(I7,I5,3F11.3,3ES12.4)') i, j-1, EnuA(l), T(i), Mun(i) - Mup(i) - Mul(i), &
      Opa_W(l,i,j,3), Opa_W(l,i,j,4), DB(l,i,j)
  END DO

  WRITE(*,'(/,A)')        ' Timings [s]'
  WRITE(*,'(A,F10.3)')    '   W (Opacity_CC_2D)           :', t_W
  WRITE(*,'(A,F10.3)')    '   P (PrecomputedNodes)        :', t_P
  WRITE(*,'(A,F10.3)')    '   I (integration wrappers)    :', t_I
  WRITE(*,'(A,F10.3)')    '   G (GSI reference)           :', t_G
  WRITE(*,'(/,A)')        ' Per-sample values written to AbEm2DVsGSI.dat'

  DEALLOCATE(rWP, rWI, rWG, rWGs, rWAG, rDB)
  DEALLOCATE(Opa_W, Opa_P, Opa_I, Opa_G, DB)
  DEALLOCATE(T, Rho, Ye, Ym, Mue, Mum, Mul, Mun, Mup, Un, Up, EffMassn, EffMassp)

CONTAINS

  LOGICAL FUNCTION IsNaN_DP(x)
    REAL(DP), INTENT(IN) :: x
    IsNaN_DP = (x /= x)
  END FUNCTION IsNaN_DP

  REAL(DP) FUNCTION RelDiff(a, b)
    REAL(DP), INTENT(IN) :: a, b
    REAL(DP) :: denom
    denom = MAX(ABS(a), ABS(b))
    IF (denom > 0.0d0) THEN
      RelDiff = ABS(a - b) / denom
    ELSE
      RelDiff = 0.0d0
    END IF
  END FUNCTION RelDiff

  REAL(DP) FUNCTION MaxOf(x, n)
    REAL(DP), INTENT(IN) :: x(:)
    INTEGER , INTENT(IN) :: n
    IF (n > 0) THEN
      MaxOf = MAXVAL(x(1:n))
    ELSE
      MaxOf = 0.0d0
    END IF
  END FUNCTION MaxOf

  REAL(DP) FUNCTION MedianOf(x, n)
    REAL(DP), INTENT(IN) :: x(:)
    INTEGER , INTENT(IN) :: n
    REAL(DP) :: w(n), tmp
    INTEGER  :: ii, jj
    IF (n <= 0) THEN
      MedianOf = 0.0d0
      RETURN
    END IF
    w = x(1:n)
    ! insertion sort (n is small)
    DO ii = 2, n
      tmp = w(ii)
      jj  = ii - 1
      DO WHILE (jj >= 1)
        IF (w(jj) <= tmp) EXIT
        w(jj+1) = w(jj)
        jj = jj - 1
      END DO
      w(jj+1) = tmp
    END DO
    IF (MOD(n,2) == 1) THEN
      MedianOf = w((n+1)/2)
    ELSE
      MedianOf = 0.5d0 * ( w(n/2) + w(n/2+1) )
    END IF
  END FUNCTION MedianOf

END PROGRAM wlTestAbEm2DVsGSI

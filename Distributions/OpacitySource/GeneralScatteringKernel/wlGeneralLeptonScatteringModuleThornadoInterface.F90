MODULE wlGeneralLeptonScatteringModuleThornadoInterface

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    pi, Gw_MeV, hbarc, sin2W, me, mmu, cvel
  USE wlGeneralLeptonScatteringModule, ONLY: &
    SelectLambdaFromProcessIndex, &
    Compute_Is_Integrals,         &
    A1_f, B1_f, C1_f, C3_f,       &
    gauleg, LegendrePolynomial

  IMPLICIT NONE
  PRIVATE

  INTEGER , PARAMETER, PUBLIC :: iProcessMax = 33
  REAL(DP), ALLOCATABLE    :: LambdaScatteringArray(:)
  REAL(DP), ALLOCATABLE    :: A_Scat(:, :, :, :)
  REAL(DP), ALLOCATABLE    :: B_Scat(:, :, :, :, :, :)
  REAL(DP), ALLOCATABLE    :: C_Scat(:, :, :, :, :, :)

  INTEGER                  :: nThetaScattering = 24
  REAL(DP), ALLOCATABLE    :: xa_cos_theta(:)
  REAL(DP), ALLOCATABLE    :: wa_cos_theta(:)

  INTEGER , PARAMETER      :: iABC_el = 1
  INTEGER , PARAMETER      :: iABC_mu = 2

  INTEGER , PARAMETER      :: nDistinctCases = 8
  REAL(DP), DIMENSION(8,2) :: MassPairDistinct
  INTEGER , PARAMETER      :: i_em_em = 1 ! particle 2 is e-  and particle 4 is e-
  INTEGER , PARAMETER      :: i_ep_ep = 2 ! particle 2 is e+  and particle 4 is e+
  INTEGER , PARAMETER      :: i_mm_mm = 3 ! particle 2 is mu- and particle 4 is mu-
  INTEGER , PARAMETER      :: i_mp_mp = 4 ! particle 2 is mu+ and particle 4 is mu+
  INTEGER , PARAMETER      :: i_em_mm = 5 ! particle 2 is e-  and particle 4 is mu-
  INTEGER , PARAMETER      :: i_mm_em = 6 ! particle 2 is mu- and particle 4 is e-
  INTEGER , PARAMETER      :: i_ep_mp = 7 ! particle 2 is e+  and particle 4 is mu+
  INTEGER , PARAMETER      :: i_mp_ep = 8 ! particle 2 is mu+ and particle 4 is e+

  INTEGER , DIMENSION(8), PARAMETER :: &
    iABC_m2 = (/ iABC_el, iABC_el, iABC_mu, iABC_mu, &
                 iABC_el, iABC_mu, iABC_el, iABC_mu /)
  INTEGER , DIMENSION(8), PARAMETER :: &
    iABC_m4 = (/ iABC_el, iABC_el, iABC_mu, iABC_mu, &
                 iABC_mu, iABC_el, iABC_mu, iABC_el /)

  INTEGER  :: iDistinctMap(iProcessMax)
  REAL(DP) :: lam1(iProcessMax), lam2(iProcessMax), lam3(iProcessMax)

  REAL(DP) :: conv_fac = 2.0d0*pi / ( (2.0d0 * pi)**3 * hbarc )

  PUBLIC :: InitGeneralScatteringKernels
  PUBLIC :: CalculateAllRout
  PUBLIC :: CalculateAllPhout

CONTAINS

  SUBROUTINE InitGeneralScatteringKernels(E1_array, E3_array, nE1, nE3, nTheta)

    REAL(DP), INTENT(IN) :: E1_array(nE1), E3_array(nE3)
    INTEGER , INTENT(IN) :: nE1, nE3, nTheta
    ! Here is where A, B, C are calculated
    ! Lambda are calculated
    ! GauLeg are calculated
    ! Potentially this is where the
    ! selection of which processes to include is done

    REAL(DP) :: E1, E3, costh, Delta, Delta5
    INTEGER  :: iE1, iE3, iTh, iProcess

    nThetaScattering = nTheta

    ALLOCATE( LambdaScatteringArray(iProcessMax) )
    ALLOCATE( A_Scat(nE1, nE3, nThetaScattering, 2) )
    ALLOCATE( B_Scat(nE1, nE3, nThetaScattering, 2, 2, 2) ) ! This could also be (nE1, nE2, nTheta, 3, 3) but more readable this way
    ALLOCATE( C_Scat(nE1, nE3, nThetaScattering, 2, 2, 3) )

    ! Take care of angular quadrature
    ALLOCATE( xa_cos_theta(nThetaScattering), wa_cos_theta(nThetaScattering) )
    CALL gauleg( -1.0d0, 1.0d0, xa_cos_theta, wa_cos_theta, nThetaScattering )

    ! Here we might want to iterate from iE_B to iE_S or whatever? Not sure
    DO iE1=1,nE1
      DO iE3=1,nE3
        DO iTh=1,nThetaScattering

          E1    = E1_array(iE1)
          E3    = E3_array(iE3)
          costh = xa_cos_theta(iTh)
          Delta  = SQRT(MAX(E1**2 - 2.0d0*E1*E3*costh + E3**2, 0.0d0))
          Delta5 = Delta**5

          A_Scat(iE1, iE3, iTh, 1) = A1_f(  E1,  E3, costh ) 
          A_Scat(iE1, iE3, iTh, 2) = A1_f( -E3, -E1, costh )
          
          B_Scat(iE1, iE3, iTh, iABC_el, iABC_el, 1) = B1_f(  E1,  E3, costh, 0.0d0 )
          B_Scat(iE1, iE3, iTh, iABC_el, iABC_el, 2) = B1_f( -E3, -E1, costh, 0.0d0 )

          B_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 1) = B1_f(  E1,  E3, costh, me**2 - mmu**2 )
          B_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 2) = B1_f( -E3, -E1, costh, me**2 - mmu**2 )

          B_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 1) = B1_f(  E1,  E3, costh, mmu**2 - me**2 )
          B_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 2) = B1_f( -E3, -E1, costh, mmu**2 - me**2 )
 
          B_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 1) = B1_f(  E1,  E3, costh, 0.0d0 )
          B_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 2) = B1_f( -E3, -E1, costh, 0.0d0 )
  
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_el, 1) = C1_f(  E1,  E3, costh, 0.0d0, me**2 )
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_el, 2) = C1_f( -E3, -E1, costh, 0.0d0, me**2 )
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_el, 3) = C3_f(  E1,  E3, costh, me, me )

          C_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 1) = C1_f(  E1,  E3, costh, me**2 - mmu**2, me**2 )
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 2) = C1_f( -E3, -E1, costh, me**2 - mmu**2, me**2 )
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 3) = C3_f(  E1,  E3, costh, me, mmu )

          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 1) = C1_f(  E1,  E3, costh, mmu**2 - me**2, mmu**2 )
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 2) = C1_f( -E3, -E1, costh, mmu**2 - me**2, mmu**2 )
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 3) = C3_f(  E1,  E3, costh, mmu, me ) ! Same as the one above actually
 
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 1) = C1_f(  E1,  E3, costh, 0.0d0, mmu**2 )
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 2) = C1_f( -E3, -E1, costh, 0.0d0, mmu**2 )
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 3) = C3_f(  E1,  E3, costh, mmu, mmu )

          ! Might as well build Delta into A, B, C while we are here
          Delta = SQRT(MAX(E1**2 - 2.0d0*E1*E3*costh + E3**2, 0.0d0))
          A_Scat(iE1, iE3, iTh, :      ) = A_Scat(iE1, iE3, iTh, :      ) / (16.0d0 * pi * Delta5)
          B_Scat(iE1, iE3, iTh, :, :, :) = B_Scat(iE1, iE3, iTh, :, :, :) / (16.0d0 * pi * Delta5)
          C_Scat(iE1, iE3, iTh, :, :, :) = C_Scat(iE1, iE3, iTh, :, :, :) / (16.0d0 * pi * Delta5)

        ENDDO
      ENDDO
    ENDDO

    ! This is to compute Distinct Is
    MassPairDistinct(i_em_em,1) = me
    MassPairDistinct(i_em_em,2) = me

    MassPairDistinct(i_ep_ep,1) = me
    MassPairDistinct(i_ep_ep,2) = me

    MassPairDistinct(i_mm_mm,1) = mmu
    MassPairDistinct(i_mm_mm,2) = mmu

    MassPairDistinct(i_mp_mp,1) = mmu
    MassPairDistinct(i_mp_mp,2) = mmu

    MassPairDistinct(i_em_mm,1) = me
    MassPairDistinct(i_em_mm,2) = mmu

    MassPairDistinct(i_mm_em,1) = mmu
    MassPairDistinct(i_mm_em,2) = me

    MassPairDistinct(i_ep_mp,1) = me
    MassPairDistinct(i_ep_mp,2) = mmu

    MassPairDistinct(i_mp_ep,1) = mmu
    MassPairDistinct(i_mp_ep,2) = me

    ! Now initialize Lambdas and map each process to the correct I
    DO iProcess = 1, iProcessMax
      CALL SelectLambdaFromProcessIndex( &
           iProcess, lam1(iProcess), lam2(iProcess), lam3(iProcess))
      CALL SelectiDistinctFromProcessIndex(&
           iProcess, iDistinctMap(iProcess))
    ENDDO

  END SUBROUTINE InitGeneralScatteringKernels

  SUBROUTINE SelectiDistinctFromProcessIndex(ProcessIndex, iDistinct)

    ! Same structure as SelectLambdaFromProcessIndex

    INTEGER, INTENT(IN)  :: ProcessIndex
    INTEGER, INTENT(OUT) :: iDistinct

    SELECT CASE(ProcessIndex)
    ! ============================================================
    ! ELASTIC NC SCATTERING ONLY (same in/out flavor): 1..24
    ! ============================================================

    ! 1) nu_e + e-
    CASE(1); iDistinct = i_em_em
    ! 2) nu_e + e+
    CASE(2); iDistinct = i_ep_ep
    ! 3) nu_bar_e + e-
    CASE(3); iDistinct = i_em_em
    ! 4) nu_bar_e + e+
    CASE(4); iDistinct = i_ep_ep

    ! 5) nu_mu + e-
    CASE(5); iDistinct = i_em_em
    ! 6) nu_mu + e+
    CASE(6); iDistinct = i_ep_ep
    ! 7) nu_bar_mu + e-
    CASE(7); iDistinct = i_em_em
    ! 8) nu_bar_mu + e+
    CASE(8); iDistinct = i_ep_ep

    ! 9) nu_tau + e-
    CASE(9); iDistinct = i_em_em
    ! 10) nu_tau + e+
    CASE(10); iDistinct = i_ep_ep
    ! 11) nu_bar_tau + e-
    CASE(11); iDistinct = i_em_em
    ! 12) nu_bar_tau + e+
    CASE(12); iDistinct = i_ep_ep

    ! 13) nu_e + mu-
    CASE(13); iDistinct = i_mm_mm
    ! 14) nu_e + mu+
    CASE(14); iDistinct = i_mp_mp
    ! 15) nu_bar_e + mu-
    CASE(15); iDistinct = i_mm_mm
    ! 16) nu_bar_e + mu+
    CASE(16); iDistinct = i_mp_mp

    ! 17) nu_mu + mu-
    CASE(17); iDistinct = i_mm_mm
    ! 18) nu_mu + mu+
    CASE(18); iDistinct = i_mp_mp
    ! 19) nu_bar_mu + mu-
    CASE(19); iDistinct = i_mm_mm
    ! 20) nu_bar_mu + mu+
    CASE(20); iDistinct = i_mp_mp

    ! 21) nu_tau + mu-
    CASE(21); iDistinct = i_mm_mm
    ! 22) nu_tau + mu+
    CASE(22); iDistinct = i_mp_mp
    ! 23) nu_bar_tau + mu-
    CASE(23); iDistinct = i_mm_mm
    ! 24) nu_bar_tau + mu+
    CASE(24); iDistinct = i_mp_mp

    ! ============================================================
    ! INELASTIC FLAVOR EXCHANGE (LFE) & CONVERSION (LFC): 25..32
    ! ============================================================
    
    ! LFE (t-channel W-exchange) 
    ! 25) nu_mu + e- -> nu_e + mu-
    CASE(25); iDistinct = i_em_mm
    ! 26) nu_e + mu- -> nu_mu + e-
    CASE(26); iDistinct = i_mm_em
    ! 27) nu_bar_mu + e+ -> nu_bar_e + mu+
    CASE(27); iDistinct = i_ep_mp
    ! 28) nu_bar_e + mu+ -> nu_bar_mu + e+
    CASE(28); iDistinct = i_mp_ep

    ! LFC (s-channel W-exchange)
    ! 29) nu_mu + mu+ -> nu_e + e+
    CASE(29); iDistinct = i_mp_ep
    ! 30) nu_e + e+ -> nu_mu + mu+
    CASE(30); iDistinct = i_ep_mp
    ! 31) nu_bar_mu + mu- -> nu_bar_e + e-
    CASE(31); iDistinct = i_mm_em
    ! 32) nu_bar_e + e- -> nu_bar_mu + mu-
    CASE(32); iDistinct = i_em_mm

    CASE(33); iDistinct = 0 ! Not sure about this one deal with it later

    CASE DEFAULT
      WRITE(*,*) 'Error: Unrecognized ProcessIndex: ', ProcessIndex
      STOP
    END SELECT

  END SUBROUTINE SelectiDistinctFromProcessIndex

  SUBROUTINE CalculateAllRout( E1, E3, costh, T, Mu_e, Mu_mu, iE1, iE3, iTh, Rout )

    REAL(DP), INTENT(IN)  :: E1, E3, costh, T, Mu_e, Mu_mu
    INTEGER , INTENT(IN)  :: iE1, iE3, iTh
    REAL(DP), INTENT(OUT) :: Rout(iProcessMax)
    REAL(DP) :: I0(nDistinctCases), R1(nDistinctCases)
    REAL(DP) :: I1(nDistinctCases), R2(nDistinctCases)
    REAL(DP) :: I2(nDistinctCases), R3(nDistinctCases)

    INTEGER  :: iDistinct, iProcess, Index2, Index4

    CALL CalculateDistinctI0I1I2( E1, E3, costh, T, Mu_e, Mu_mu, I0, I1, I2 )
    DO iDistinct = 1, nDistinctCases
      Index2 = iABC_m2( iDistinct )
      Index4 = iABC_m4( iDistinct )
      CALL CalculateR1R2R3( E1, E3, costh, &
          A_Scat(iE1, iE3, iTh, :), &
          B_Scat(iE1, iE3, iTh, Index2, Index4, :), &
          C_Scat(iE1, iE3, iTh, Index2, Index4, :), &
          I0(iDistinct), I1(iDistinct), I2(iDistinct), &
          R1(iDistinct), R2(iDistinct), R3(iDistinct) )
    ENDDO

    DO iProcess = 1, iProcessMax
      Rout(iProcess) = lam1(iProcess)*R1(iDistinctMap(iProcess)) + &
                       lam2(iProcess)*R2(iDistinctMap(iProcess)) + &
                       lam3(iProcess)*R3(iDistinctMap(iProcess))
    ENDDO

  END SUBROUTINE CalculateAllRout

  SUBROUTINE CalculateAllPhout( iE1, iE3, E1, E3, T, Mu_e, Mu_mu, Phout, nL )

    INTEGER , INTENT(IN)  :: iE1, iE3
    REAL(DP), INTENT(IN)  :: E1, E3, T, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: Phout(iProcessMax,nL)
    INTEGER , INTENT(IN)  :: nL

    REAL(DP) :: Rout(iProcessMax)
    REAL(DP) :: costh, Delta_mu, exponent, Pl_mu
    INTEGER  :: iTh, nTheta, iL, iProcess

    Phout = 0.0d0
    DO iTh=1, nThetaScattering
      costh = xa_cos_theta(iTh)
      CALL CalculateAllRout( E1, E3, costh, T, Mu_e, Mu_mu, iE1, iE3, iTh, Rout )

      ! --- Accumulate the integral for each moment ---
      DO iL = 1, nL
         Pl_mu = LegendrePolynomial(iL-1, costh)
        DO iProcess=1,iProcessMax
          Phout(iProcess,iL) = Phout(iProcess,iL) + Rout(iProcess) * Pl_mu * wa_cos_theta(iTh)
        END DO
      ENDDO
    ENDDO
    ! factors of 1/2 and 3/2 to match Bruenn and then factor of two to match thornado convention?
    Phout(:,1) = Phout(:,1) * conv_fac * 0.5d0 * 2.0d0
    Phout(:,2) = Phout(:,2) * conv_fac * 1.5d0 * 2.0d0

  END SUBROUTINE CalculateAllPhout

  SUBROUTINE CalculateR1R2R3( E1, E3, costh, A, B, C, I0_val, I1_val, I2_val, R1, R2, R3 )

    REAL(DP), INTENT(IN)  :: E1, E3, costh, I0_val, I1_val, I2_val
    REAL(DP), INTENT(IN)  :: A(2), B(2), C(3) ! Notice that A3 and B3 are identically zero
    REAL(DP), INTENT(OUT) :: R1, R2, R3

    REAL(DP) :: Delta

    ! --- Delta (Guo Eq. 4) ---
    Delta = SQRT(MAX(E1**2 - 2.0d0*E1*E3*costh + E3**2, 0.0d0))   ! Guo Eq. 4

    ! --- Notice that 1.0d0 / (16.0d0 * pi * Delta5) is built into A, B, C ---

    R1 = A(1)*I2_val + B(1)*I1_val + C(1)*I0_val   ! Guo Eq. 3
    R2 = A(2)*I2_val + B(2)*I1_val + C(2)*I0_val   ! Guo Eq. 3
    R3 =                             C(3)*I0_val   ! Guo Eq. 3 

  END SUBROUTINE CalculateR1R2R3

  SUBROUTINE SetMuDistinctCases(Mu_e, Mu_mu, MuPair)

    REAL(DP), INTENT(IN)  :: Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: MuPair(nDistinctCases,2)

    MuPair(i_em_em,1) = Mu_e
    MuPair(i_em_em,2) = Mu_e

    MuPair(i_ep_ep,1) = -Mu_e
    MuPair(i_ep_ep,2) = -Mu_e

    MuPair(i_mm_mm,1) = Mu_mu
    MuPair(i_mm_mm,2) = Mu_mu

    MuPair(i_mp_mp,1) = -Mu_mu
    MuPair(i_mp_mp,2) = -Mu_mu

    MuPair(i_em_mm,1) = Mu_e
    MuPair(i_em_mm,2) = Mu_mu

    MuPair(i_mm_em,1) = Mu_mu
    MuPair(i_mm_em,2) = Mu_e

    MuPair(i_ep_mp,1) = -Mu_e
    MuPair(i_ep_mp,2) = -Mu_mu

    MuPair(i_mp_ep,1) = -Mu_mu
    MuPair(i_mp_ep,2) = -Mu_e

  END SUBROUTINE SetMuDistinctCases

  SUBROUTINE CalculateDistinctI0I1I2( E1, E3, costh, T, Mu_e, Mu_mu, I0, I1, I2 )

    REAL(DP), INTENT(IN)  :: E1, E3, costh, T, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: I0(nDistinctCases), I1(nDistinctCases), I2(nDistinctCases)
    
    REAL(DP) :: Delta, xMu_l2, xMu_l4, m2, m4, disc, k_val, E_minus, Q
    REAL(DP) :: MuPairDistinct(nDistinctCases,2)
    INTEGER  :: iCase

    I0 = 0.0d0
    I1 = 0.0d0
    I2 = 0.0d0

    Delta = SQRT(MAX(E1**2 - 2.0d0*E1*E3*costh + E3**2, 0.0d0))   ! Guo Eq. 4
    IF (Delta < 1.0d-10) RETURN   ! pathological collinear case

    ! --- k (Guo Eq. 5): diverges at forward scattering mu=1 ---
    IF (1.0d0 - costh < 1.0d-10) RETURN   ! mu -> 1: E_minus -> inf, no phase space
    
    ! There are nDistinctCases different cases:

    CALL SetMuDistinctCases(Mu_e, Mu_mu, MuPairDistinct)
    DO iCase=1,nDistinctCases
      
      k_val = Q / (E1 * E3 * (1.0d0 - costh))   ! Guo Eq. 5
      disc = (1.0d0 + k_val)**2 + 2.0d0 * m2**2 / (E1 * E3 * (1.0d0 - costh)) ! Guo Eq. 4
      IF (disc < 0.0d0) RETURN

      m2     = MassPairDistinct(iCase,1)
      m4     = MassPairDistinct(iCase,1)
 
      xMu_l2 = MuPairDistinct(iCase,1)
      xMu_l4 = MuPairDistinct(iCase,2)

      Q       = 0.5d0 * (m4**2 - m2**2)
      E_minus = 0.5d0 * ( (E3 - E1)*(1.0d0 + k_val) + Delta*SQRT(disc) )   ! Guo Eq. 4
      E_minus = MAX(E_minus, m2)

      CALL Compute_Is_Integrals( E1, E3, E_minus, T, xMu_l2, xMu_l4, I0(iCase), I1(iCase), I2(iCase) )

    ENDDO

  END SUBROUTINE CalculateDistinctI0I1I2


END MODULE wlGeneralLeptonScatteringModuleThornadoInterface

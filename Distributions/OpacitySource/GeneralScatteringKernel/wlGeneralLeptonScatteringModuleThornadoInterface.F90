MODULE wlGeneralLeptonScatteringModuleThornadoInterface

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    pi, Gw_MeV, hbarc, sin2W, me, mmu, cvel
  USE wlGeneralLeptonScatteringModule, ONLY: &
    SelectLambdaFromProcessIndex, &
    Compute_Is_Integrals, FD,     &
    A1_f, B1_f, C1_f, C3_f,       &
    gauleg, LegendrePolynomial,   &
    gaulag, CalculateCollinearH0

  IMPLICIT NONE
  PRIVATE

  INTEGER , PARAMETER, PUBLIC :: iProcessMin_Default = 1
  INTEGER , PARAMETER, PUBLIC :: iProcessMax_Default = 34
  INTEGER                     :: iProcessMin
  INTEGER                     :: iProcessMax
  REAL(DP), PARAMETER         :: ElectronMass        = me
  REAL(DP), PARAMETER         :: MuonMass            = mmu
  REAL(DP), PARAMETER         :: E2_max_Default      = 100.0d0
  INTEGER , PARAMETER         :: nE2_Default         = 32
  INTEGER                     :: nE2_Lepton

  REAL(DP), ALLOCATABLE    :: LambdaScatteringArray(:)
  REAL(DP), ALLOCATABLE    :: A_Scat(:, :, :, :)
  REAL(DP), ALLOCATABLE    :: B_Scat(:, :, :, :, :, :)
  REAL(DP), ALLOCATABLE    :: C_Scat(:, :, :, :, :, :)

  REAL(DP), ALLOCATABLE    :: A_IMD(:, :, :, :)
  REAL(DP), ALLOCATABLE    :: B_IMD(:, :, :, :)
  REAL(DP), ALLOCATABLE    :: C_IMD(:, :, :, :)

  ! IMD threshold (Guo Eq. 16): kernel vanishes for k < k0
  REAL(DP), PARAMETER      :: k0_IMD = (MuonMass + ElectronMass) / (MuonMass - ElectronMass)

  INTEGER                  :: nThetaScattering
  REAL(DP), ALLOCATABLE    :: xa_cos_theta(:), wa_cos_theta(:)
  REAL(DP), ALLOCATABLE    :: xa_GLeg_ref(:), wa_GLeg_ref(:)
  REAL(DP), ALLOCATABLE    :: xa_GLag(:), wa_GLag(:)

  INTEGER , PARAMETER      :: iABC_el = 1
  INTEGER , PARAMETER      :: iABC_mu = 2

  INTEGER , PARAMETER      :: nDistinctCasesMax = 9
  INTEGER                  :: nDistinctCases
  REAL(DP), DIMENSION(nDistinctCasesMax,2) :: MassPairDistinct
  INTEGER , PARAMETER      :: i_em_em = 1 ! particle 2 is e-  and particle 4 is e-
  INTEGER , PARAMETER      :: i_ep_ep = 2 ! particle 2 is e+  and particle 4 is e+
  INTEGER , PARAMETER      :: i_mm_mm = 3 ! particle 2 is mu- and particle 4 is mu-
  INTEGER , PARAMETER      :: i_mp_mp = 4 ! particle 2 is mu+ and particle 4 is mu+
  INTEGER , PARAMETER      :: i_em_mm = 5 ! particle 2 is e-  and particle 4 is mu-
  INTEGER , PARAMETER      :: i_mm_em = 6 ! particle 2 is mu- and particle 4 is e-
  INTEGER , PARAMETER      :: i_ep_mp = 7 ! particle 2 is e+  and particle 4 is mu+
  INTEGER , PARAMETER      :: i_mp_ep = 8 ! particle 2 is mu+ and particle 4 is e+
  INTEGER , PARAMETER      :: i_IMD = 9 ! IMD: particle 2 is e-, particle 4 is mu- (applies to 33 and 34)

  INTEGER , DIMENSION(8), PARAMETER :: &
    iABC_m2 = (/ iABC_el, iABC_el, iABC_mu, iABC_mu, &
                 iABC_el, iABC_mu, iABC_el, iABC_mu /)
  INTEGER , DIMENSION(8), PARAMETER :: &
    iABC_m4 = (/ iABC_el, iABC_el, iABC_mu, iABC_mu, &
                 iABC_mu, iABC_el, iABC_mu, iABC_el /)

  INTEGER                             :: iDistinctMap(iProcessMax_Default)
  INTEGER , ALLOCATABLE, DIMENSION(:) :: ChosenToTrueCaseMap(:)
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: lam1(:), lam2(:), lam3(:)

  REAL(DP) :: conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc )

  PUBLIC :: InitGeneralScatteringKernels
  PUBLIC :: FinalizeGeneralScatteringKernels
  PUBLIC :: CalculateAllRout
  PUBLIC :: CalculateAllPhout
  PUBLIC :: CalculateAllRoutIntegrated

CONTAINS

  SUBROUTINE InitGeneralScatteringKernels(E1_array, E3_array, nE1, nE3, nTheta, &
                                          iProcessMin_Option, iProcessMax_Option, &
                                          E2_max_Option, nE2_Option)

    REAL(DP), INTENT(IN)           :: E1_array(nE1), E3_array(nE3)
    INTEGER , INTENT(IN)           :: nE1, nE3, nTheta
    INTEGER , INTENT(IN), OPTIONAL :: iProcessMin_Option, iProcessMax_Option
    REAL(DP), INTENT(IN), OPTIONAL :: E2_max_Option
    INTEGER , INTENT(IN), OPTIONAL :: nE2_Option

    REAL(DP) :: E1, E2, E3, costh, Delta, Delta5, E2_max
    INTEGER  :: iE1, iE2, iE3, iTh, iProcess, iCase
    
    REAL(DP) :: DeltaA, DeltaA5, Q_mu_e, Q_e_mu

    INTEGER, DIMENSION(nDistinctCasesMax) :: unique_cases
    INTEGER :: idx, current_val
    LOGICAL :: IncludeIMD

    nThetaScattering = nTheta

    IF ( PRESENT(nE2_Option) ) THEN
      nE2_Lepton = nE2_Option
    ELSE
      nE2_Lepton = nE2_Default
    ENDIF

    IF ( PRESENT(E2_max_Option) ) THEN
      E2_max = E2_max_Option
    ELSE
      E2_max = E2_max_Default
    ENDIF

    IF ( PRESENT(iProcessMin_Option) ) THEN
      iProcessMin = iProcessMin_Option
    ELSE
      iProcessMin = iProcessMin_Default
    ENDIF

    IF ( PRESENT(iProcessMax_Option) ) THEN
      iProcessMax = iProcessMax_Option
    ELSE
      iProcessMax = iProcessMax_Default
    ENDIF

    IncludeIMD = ( iProcessMax >= 33 )

    ALLOCATE( lam1(iProcessMax - iProcessMin + 1) )
    ALLOCATE( lam2(iProcessMax - iProcessMin + 1) )
    ALLOCATE( lam3(iProcessMax - iProcessMin + 1) )
    ALLOCATE( A_Scat(nE1, nE3, nThetaScattering, 2) )
    ALLOCATE( B_Scat(nE1, nE3, nThetaScattering, 2, 2, 2) ) ! This could also be (nE1, nE3, nTheta, 3, 3) but more readable this way
    ALLOCATE( C_Scat(nE1, nE3, nThetaScattering, 2, 2, 3) )
      
    IF (IncludeIMD) THEN
      ! Only one flavor combination (e- in, mu- out) and no C3 slot (lam3 = 0
      ! for both IMD processes), hence the small last dimension.
      ALLOCATE( A_IMD(nE1, nE3, nThetaScattering, 2) )
      ALLOCATE( B_IMD(nE1, nE3, nThetaScattering, 2) )
      ALLOCATE( C_IMD(nE1, nE3, nThetaScattering, 2) )
    ENDIF

    ! Take care of angular quadrature
    ALLOCATE( xa_cos_theta(nThetaScattering), wa_cos_theta(nThetaScattering) )
    CALL gauleg( -1.0d0, 1.0d0, xa_cos_theta, wa_cos_theta, nThetaScattering )

    Q_mu_e = 0.5d0 * (MuonMass**2 - ElectronMass**2)
    Q_e_mu = - Q_mu_e

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

          B_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 1) = B1_f(  E1,  E3, costh, Q_mu_e )
          B_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 2) = B1_f( -E3, -E1, costh, Q_mu_e )

          B_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 1) = B1_f(  E1,  E3, costh, Q_e_mu )
          B_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 2) = B1_f( -E3, -E1, costh, Q_e_mu )
 
          B_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 1) = B1_f(  E1,  E3, costh, 0.0d0 )
          B_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 2) = B1_f( -E3, -E1, costh, 0.0d0 )
  
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_el, 1) = C1_f(  E1,  E3, costh, 0.0d0, ElectronMass )
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_el, 2) = C1_f( -E3, -E1, costh, 0.0d0, ElectronMass )
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_el, 3) = C3_f(  E1,  E3, costh, ElectronMass, ElectronMass )

          C_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 1) = C1_f(  E1,  E3, costh, Q_mu_e, ElectronMass )
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 2) = C1_f( -E3, -E1, costh, Q_mu_e, ElectronMass )
          C_Scat(iE1, iE3, iTh, iABC_el, iABC_mu, 3) = C3_f(  E1,  E3, costh, ElectronMass, MuonMass )

          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 1) = C1_f(  E1,  E3, costh, Q_e_mu, MuonMass )
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 2) = C1_f( -E3, -E1, costh, Q_e_mu, MuonMass )
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_el, 3) = C3_f(  E1,  E3, costh, MuonMass, ElectronMass ) ! Same as the one above actually
 
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 1) = C1_f(  E1,  E3, costh, 0.0d0, MuonMass )
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 2) = C1_f( -E3, -E1, costh, 0.0d0, MuonMass )
          C_Scat(iE1, iE3, iTh, iABC_mu, iABC_mu, 3) = C3_f(  E1,  E3, costh, MuonMass, MuonMass )

          ! Might as well build Delta into A, B, C while we are here
          Delta = SQRT(MAX(E1**2 - 2.0d0*E1*E3*costh + E3**2, 0.0d0))
          A_Scat(iE1, iE3, iTh, :      ) = A_Scat(iE1, iE3, iTh, :      ) / (16.0d0 * pi * Delta5)
          B_Scat(iE1, iE3, iTh, :, :, :) = B_Scat(iE1, iE3, iTh, :, :, :) / (16.0d0 * pi * Delta5)
          C_Scat(iE1, iE3, iTh, :, :, :) = C_Scat(iE1, iE3, iTh, :, :, :) / (16.0d0 * pi * Delta5)

          IF (IncludeIMD) THEN
            ! IMD coefficients, the "1" set is process 33
            ! (nu_bar_e has energy E1), the "2" set is its E1 <-> E3 swap
            ! (process 34). Particle 2 is the electron in BOTH orientations.
            ! Notice that C is not needed since lam3 is always zero

            DeltaA  = SQRT(MAX(E1**2 + 2.0d0*E1*E3*costh + E3**2, 0.0d0))
            DeltaA5 = DeltaA**5

            A_IMD(iE1, iE3, iTh, 1) = A1_f(  E1, -E3, costh )        / (16.0d0 * pi * DeltaA5)
            A_IMD(iE1, iE3, iTh, 2) = A1_f(  E3, -E1, costh )        / (16.0d0 * pi * DeltaA5)

            B_IMD(iE1, iE3, iTh, 1) = B1_f(  E1, -E3, costh, Q_mu_e ) / (16.0d0 * pi * DeltaA5)
            B_IMD(iE1, iE3, iTh, 2) = B1_f(  E3, -E1, costh, Q_mu_e ) / (16.0d0 * pi * DeltaA5)

            C_IMD(iE1, iE3, iTh, 1) = C1_f(  E1, -E3, costh, Q_mu_e, ElectronMass ) / (16.0d0 * pi * DeltaA5)
            C_IMD(iE1, iE3, iTh, 2) = C1_f(  E3, -E1, costh, Q_mu_e, ElectronMass ) / (16.0d0 * pi * DeltaA5)
          ENDIF

        ENDDO
      ENDDO
    ENDDO

    ! --- Preallocate collinear grids (generated ONCE, reused every call) ---
    ALLOCATE( xa_GLeg_ref(nE2_Lepton), wa_GLeg_ref(nE2_Lepton) )
    ALLOCATE( xa_GLag(nE2_Lepton),     wa_GLag(nE2_Lepton) )

    ! Canonical GL panel on [-1,1]; affine-mapped to [0,E1] at use time.
    CALL gauleg( -1.0d0, 1.0d0, xa_GLeg_ref, wa_GLeg_ref, nE2_Lepton )

    ! Raw Gauss-Laguerre (alpha=0) nodes/weights for int_0^inf e^{-x} g(x) dx.
    CALL gaulag( xa_GLag, wa_GLag, nE2_Lepton, 0.0d0 )

    ! This is to compute Distinct Is
    MassPairDistinct(i_em_em,1) = ElectronMass
    MassPairDistinct(i_em_em,2) = ElectronMass

    MassPairDistinct(i_ep_ep,1) = ElectronMass
    MassPairDistinct(i_ep_ep,2) = ElectronMass

    MassPairDistinct(i_mm_mm,1) = MuonMass
    MassPairDistinct(i_mm_mm,2) = MuonMass

    MassPairDistinct(i_mp_mp,1) = MuonMass
    MassPairDistinct(i_mp_mp,2) = MuonMass

    MassPairDistinct(i_em_mm,1) = ElectronMass
    MassPairDistinct(i_em_mm,2) = MuonMass

    MassPairDistinct(i_mm_em,1) = MuonMass
    MassPairDistinct(i_mm_em,2) = ElectronMass

    MassPairDistinct(i_ep_mp,1) = ElectronMass
    MassPairDistinct(i_ep_mp,2) = MuonMass

    MassPairDistinct(i_mp_ep,1) = MuonMass
    MassPairDistinct(i_mp_ep,2) = ElectronMass

    MassPairDistinct(i_IMD,1) = ElectronMass
    MassPairDistinct(i_IMD,2) = MuonMass

    ! Now initialize Lambdas and map each process to the correct I
    nDistinctCases  = 0
    unique_cases(:) = -1
    DO iProcess = iProcessMin, iProcessMax
       ! Storing your index expression in a temporary variable to keep things scannable
       idx = iProcess - iProcessMin + 1

       CALL SelectLambdaFromProcessIndex( iProcess, &
                  lam1(idx), &
                  lam2(idx), &
                  lam3(idx))

       ! Notice that iDistinctMap is always called with the "True" iProcess
       CALL SelectiDistinctFromProcessIndex( iProcess, &
                  iDistinctMap(iProcess))

       current_val = iDistinctMap(iProcess)

       ! 2. Check if current_val is NOT present in the populated portion of unique_cases
       ! This basically looks at whether we have encountered current_val already, and if
       ! we have not, then increases nDistinctCases by 1.
       IF (.NOT. ANY(unique_cases(1:nDistinctCases) == current_val)) THEN
          nDistinctCases = nDistinctCases + 1
          unique_cases(nDistinctCases) = current_val
       END IF
    ENDDO

    ALLOCATE( ChosenToTrueCaseMap(nDistinctCases) )
    ! This is needed because in principle you might be needing Cases
    ! 3,4,5 for example (not sure if this specific combination can actually happend since
    ! user does not have complete freedom in choosing reactions, but can only
    ! choose iProcessMin and iProcessMas, cannot choose chunks by skipping reactions
    ! in between for example). But this part will work regardless
    DO iCase=1, nDistinctCases
      ChosenToTrueCaseMap(iCase) = unique_cases(iCase)
    ENDDO

    WRITE(*,*) 'Working with', nDistinctCases, 'Distinct Scattering Cases'

  END SUBROUTINE InitGeneralScatteringKernels

  SUBROUTINE FinalizeGeneralScatteringKernels()

    DEALLOCATE( lam1, lam2, lam3, A_Scat, B_Scat, C_Scat )
    DEALLOCATE( xa_cos_theta, wa_cos_theta )
    DEALLOCATE( xa_GLeg_ref, wa_GLeg_ref, xa_GLag, wa_GLag )
    DEALLOCATE( ChosenToTrueCaseMap )
    IF ( ALLOCATED(A_IMD) ) DEALLOCATE( A_IMD, B_IMD, C_IMD )

  END SUBROUTINE FinalizeGeneralScatteringKernels


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

    CASE(33); iDistinct = i_IMD
    CASE(34); iDistinct = i_IMD ! Same as 33 but need to switch nue_bar and nu_mu

    CASE DEFAULT
      WRITE(*,*) 'Error: Unrecognized ProcessIndex: ', ProcessIndex
      STOP
    END SELECT

  END SUBROUTINE SelectiDistinctFromProcessIndex

  SUBROUTINE CalculateAllPhout( iE1, iE3, E1, E3, T, Mu_e, Mu_mu, Phout, nL )

    INTEGER , INTENT(IN)  :: iE1, iE3
    REAL(DP), INTENT(IN)  :: E1, E3, T, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: Phout(iProcessMax - iProcessMin + 1,nL)
    INTEGER , INTENT(IN)  :: nL

    REAL(DP) :: costh, Pl_mu
    INTEGER  :: iTh, iL, iProcess, idx, iTrueCase, iCase

    ! Arrays to hold the integrated components for the distinct cases ONLY
    REAL(DP) :: Int_R1(nDistinctCasesMax, nL)
    REAL(DP) :: Int_R2(nDistinctCasesMax, nL)
    REAL(DP) :: Int_R3(nDistinctCasesMax, nL)

    ! Arrays to hold the angular integrand components for a specific theta
    REAL(DP) :: R1_th(nDistinctCasesMax)
    REAL(DP) :: R2_th(nDistinctCasesMax)
    REAL(DP) :: R3_th(nDistinctCasesMax)

    Phout  = 0.0d0
    Int_R1 = 0.0d0
    Int_R2 = 0.0d0
    Int_R3 = 0.0d0

    ! IF ( ABS(E3 - E1) < 1.0d-10) THEN
    !   IF (nL == 1)  THEN
    !     CALL CalculateAllPhi0_collinear( iE1, E1, T, Mu_e, Mu_mu, Phout(:,1) )
    !   ELSE
    !     CALL CalculateAllPhoutGeneral_collinear( iE1, E1, T, Mu_e, Mu_mu, Phout, nL )
    !   ENDIF
    !   RETURN
    ! ENDIF

    ! -------------------------------------------------------------------
    ! 1. Perform the angular integration over the DISTINCT cases ONLY
    ! -------------------------------------------------------------------
    DO iTh=1, nThetaScattering
      costh = xa_cos_theta(iTh)
      
      CALL CalculateAllR1R2R3( E1, E3, costh, T, Mu_e, Mu_mu, &
                               iE1, iE3, iTh, R1_th, R2_th, R3_th )

      DO iL = 1, nL
        Pl_mu = LegendrePolynomial(iL-1, costh) * wa_cos_theta(iTh)
        
        ! Only integrate the unique cases we actually need
        DO iCase = 1, nDistinctCases
          iTrueCase = ChosenToTrueCaseMap(iCase)
          Int_R1(iTrueCase, iL) = Int_R1(iTrueCase, iL) + R1_th(iTrueCase) * Pl_mu
          Int_R2(iTrueCase, iL) = Int_R2(iTrueCase, iL) + R2_th(iTrueCase) * Pl_mu
          Int_R3(iTrueCase, iL) = Int_R3(iTrueCase, iL) + R3_th(iTrueCase) * Pl_mu
        END DO
      ENDDO
    ENDDO

    ! -------------------------------------------------------------------
    ! 2. Map the integrated distinct cases back to the full iProcess array
    ! -------------------------------------------------------------------
    DO iL = 1, nL
      DO iProcess = iProcessMin, iProcessMax
        idx = iProcess - iProcessMin + 1
        iTrueCase = iDistinctMap(iProcess)

        ! Apply the Lambdas to the integrated values once
        Phout(idx, iL) = lam1(idx) * Int_R1(iTrueCase, iL) + &
                         lam2(idx) * Int_R2(iTrueCase, iL) + &
                         lam3(idx) * Int_R3(iTrueCase, iL)
      END DO
      
    ! factors of 1/2 and 3/2 to match Bruenn (i.e. ((iL-1) + 0.5d0))
      Phout(:, iL) = Phout(:, iL) * conv_fac * ((iL-1) + 0.5d0)
    END DO

  END SUBROUTINE CalculateAllPhout

  SUBROUTINE CalculateAllRout( E1, E3, costh, T, Mu_e, Mu_mu, iE1, iE3, iTh, Rout )

    REAL(DP), INTENT(IN)  :: E1, E3, costh, T, Mu_e, Mu_mu
    INTEGER , INTENT(IN)  :: iE1, iE3, iTh
    REAL(DP), INTENT(OUT) :: Rout(iProcessMax)

    ! Temporary arrays to hold the distinct R components for this theta
    REAL(DP) :: R1_th(nDistinctCasesMax)
    REAL(DP) :: R2_th(nDistinctCasesMax)
    REAL(DP) :: R3_th(nDistinctCasesMax)
    INTEGER  :: iProcess, idx, iTrueCase

    ! 1. Get the raw distinct R components for this specific angle/energy
    CALL CalculateAllR1R2R3( E1, E3, costh, T, Mu_e, Mu_mu, &
                             iE1, iE3, iTh, R1_th, R2_th, R3_th )

    ! 2. Safely initialize Rout
    Rout = 0.0d0

    ! 3. Map the distinct cases to the full iProcess array using Lambdas
    DO iProcess = iProcessMin, iProcessMax
      idx = iProcess - iProcessMin + 1
      
      ! Fetch the true underlying physics case ID for this process
      iTrueCase = iDistinctMap(iProcess)

      ! Combine them
      Rout(iProcess) = lam1(idx) * R1_th(iTrueCase) + &
                       lam2(idx) * R2_th(iTrueCase) + &
                       lam3(idx) * R3_th(iTrueCase)
    END DO

  END SUBROUTINE CalculateAllRout

  SUBROUTINE CalculateAllRoutIntegrated( iE1, iE3, E1, E3, T, Mu_e, Mu_mu, Rout )

    REAL(DP), INTENT(IN)  :: E1, E3, T, Mu_e, Mu_mu
    INTEGER , INTENT(IN)  :: iE1, iE3
    REAL(DP), INTENT(OUT) :: Rout(iProcessMax)

    ! Temporary arrays to hold the distinct R components for this theta
    REAL(DP) :: R1_th(nDistinctCasesMax)
    REAL(DP) :: R2_th(nDistinctCasesMax)
    REAL(DP) :: R3_th(nDistinctCasesMax)
    INTEGER  :: iProcess, idx, iTrueCase, iTh, iCase
    REAL(DP) :: costh

    ! Arrays to hold the integrated components for the distinct cases ONLY
    REAL(DP) :: Int_R1(nDistinctCasesMax)
    REAL(DP) :: Int_R2(nDistinctCasesMax)
    REAL(DP) :: Int_R3(nDistinctCasesMax)

    Int_R1(:) = 0.0d0
    Int_R2(:) = 0.0d0
    Int_R3(:) = 0.0d0
    DO iTh=1, nThetaScattering
      costh = xa_cos_theta(iTh)
      
      CALL CalculateAllR1R2R3( E1, E3, costh, T, Mu_e, Mu_mu, &
                               iE1, iE3, iTh, R1_th, R2_th, R3_th )

        ! Only integrate the unique cases we actually need
        DO iCase = 1, nDistinctCases
          iTrueCase = ChosenToTrueCaseMap(iCase)
          Int_R1(iTrueCase) = Int_R1(iTrueCase) + R1_th(iTrueCase) * wa_cos_theta(iTh)
          Int_R2(iTrueCase) = Int_R2(iTrueCase) + R2_th(iTrueCase) * wa_cos_theta(iTh)
          Int_R3(iTrueCase) = Int_R3(iTrueCase) + R3_th(iTrueCase) * wa_cos_theta(iTh)
        END DO

    ENDDO
    ! 2. Safely initialize Rout
    Rout = 0.0d0

    ! -------------------------------------------------------------------
    ! 2. Map the integrated distinct cases back to the full iProcess array
    ! -------------------------------------------------------------------
    DO iProcess = iProcessMin, iProcessMax
      idx = iProcess - iProcessMin + 1
      iTrueCase = iDistinctMap(iProcess)

      ! Apply the Lambdas to the integrated values once
      Rout(idx) = lam1(idx) * Int_R1(iTrueCase) + &
                  lam2(idx) * Int_R2(iTrueCase) + &
                  lam3(idx) * Int_R3(iTrueCase)
      
    END DO
    ! Clamp it to zero...
    Rout = MAX(Rout, 0.0d0)
    
  END SUBROUTINE CalculateAllRoutIntegrated

  SUBROUTINE CalculateAllR1R2R3( E1, E3, costh, T, Mu_e, Mu_mu, iE1, iE3, iTh, R1, R2, R3 )
    REAL(DP), INTENT(IN)  :: E1, E3, costh, T, Mu_e, Mu_mu
    INTEGER , INTENT(IN)  :: iE1, iE3, iTh
    REAL(DP), INTENT(OUT) :: R1(nDistinctCasesMax)
    REAL(DP), INTENT(OUT) :: R2(nDistinctCasesMax)
    REAL(DP), INTENT(OUT) :: R3(nDistinctCasesMax)

    REAL(DP) :: I0(nDistinctCasesMax), I1(nDistinctCasesMax), I2(nDistinctCasesMax)
    INTEGER  :: iCase, iTrueCase, Index2, Index4

    ! Ensure unused cases are safely zeroed out
    R1 = 0.0d0; R2 = 0.0d0; R3 = 0.0d0

    CALL CalculateDistinctI0I1I2( E1, E3, costh, T, Mu_e, Mu_mu, I0, I1, I2 )
    
    DO iCase = 1, nDistinctCases
      iTrueCase = ChosenToTrueCaseMap(iCase)

      IF (iTrueCase == i_IMD) THEN

        ! IMD uses its own coefficient tables (built at (E1,-E3) / (E3,-E1)
        ! with the DeltaA^5 normalization). R1 is process 33 (nu_bar_e has
        ! energy E1); R2 is the E1 <-> E3 swap, i.e. process 34. R3 never
        ! contributes because lam3 = 0 for both IMD processes.
        R1(iTrueCase) = A_IMD(iE1, iE3, iTh, 1)*I2(iCase) &
                      + B_IMD(iE1, iE3, iTh, 1)*I1(iCase) &
                      + C_IMD(iE1, iE3, iTh, 1)*I0(iCase)
        R2(iTrueCase) = A_IMD(iE1, iE3, iTh, 2)*I2(iCase) &
                      + B_IMD(iE1, iE3, iTh, 2)*I1(iCase) &
                      + C_IMD(iE1, iE3, iTh, 2)*I0(iCase)
        R3(iTrueCase) = 0.0d0

      ELSE

        Index2 = iABC_m2( iTrueCase )
        Index4 = iABC_m4( iTrueCase )
        
        ! We populate the R arrays by the TRUE case index so they match what Phout expects later
        CALL CalculateR1R2R3( E1, E3, costh, &
            A_Scat(iE1, iE3, iTh, :), &
            B_Scat(iE1, iE3, iTh, Index2, Index4, :), &
            C_Scat(iE1, iE3, iTh, Index2, Index4, :), &
            I0(iCase), I1(iCase), I2(iCase), &
            R1(iTrueCase), R2(iTrueCase), R3(iTrueCase) )
      ENDIF

    ENDDO

  END SUBROUTINE CalculateAllR1R2R3

  SUBROUTINE CalculateR1R2R3( E1, E3, costh, A, B, C, I0_val, I1_val, I2_val, R1, R2, R3 )

    REAL(DP), INTENT(IN)  :: E1, E3, costh, I0_val, I1_val, I2_val
    REAL(DP), INTENT(IN)  :: A(2), B(2), C(3) ! Notice that A3 and B3 are identically zero
    REAL(DP), INTENT(OUT) :: R1, R2, R3

    ! --- Delta (Guo Eq. 4) ---
    ! Delta = SQRT(MAX(E1**2 - 2.0d0*E1*E3*costh + E3**2, 0.0d0))   ! Guo Eq. 4
    ! --- Notice that 1.0d0 / (16.0d0 * pi * Delta5) is built into A, B, C ---
    ! In principle this function does not need costh as input but for clarity
    ! we leave it in

    R1 = A(1)*I2_val + B(1)*I1_val + C(1)*I0_val   ! Guo Eq. 3
    R2 = A(2)*I2_val + B(2)*I1_val + C(2)*I0_val   ! Guo Eq. 3
    R3 =                             C(3)*I0_val   ! Guo Eq. 3 

  END SUBROUTINE CalculateR1R2R3

  SUBROUTINE SetChemPotDistinctCases(Mu_e, Mu_mu, ChemPotPair)

    REAL(DP), INTENT(IN)  :: Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: ChemPotPair(nDistinctCasesMax,2)

    ChemPotPair(i_em_em,1) = Mu_e
    ChemPotPair(i_em_em,2) = Mu_e

    ChemPotPair(i_ep_ep,1) = -Mu_e
    ChemPotPair(i_ep_ep,2) = -Mu_e

    ChemPotPair(i_mm_mm,1) = Mu_mu
    ChemPotPair(i_mm_mm,2) = Mu_mu

    ChemPotPair(i_mp_mp,1) = -Mu_mu
    ChemPotPair(i_mp_mp,2) = -Mu_mu

    ChemPotPair(i_em_mm,1) = Mu_e
    ChemPotPair(i_em_mm,2) = Mu_mu

    ChemPotPair(i_mm_em,1) = Mu_mu
    ChemPotPair(i_mm_em,2) = Mu_e

    ChemPotPair(i_ep_mp,1) = -Mu_e
    ChemPotPair(i_ep_mp,2) = -Mu_mu

    ChemPotPair(i_mp_ep,1) = -Mu_mu
    ChemPotPair(i_mp_ep,2) = -Mu_e

    ChemPotPair(i_IMD,1) = Mu_e
    ChemPotPair(i_IMD,2) = Mu_mu

  END SUBROUTINE SetChemPotDistinctCases

  SUBROUTINE CalculateDistinctI0I1I2( E1, E3, costh, T, Mu_e, Mu_mu, I0, I1, I2 )

    REAL(DP), INTENT(IN)  :: E1, E3, costh, T, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: I0(nDistinctCases), I1(nDistinctCases), I2(nDistinctCases)
    
    REAL(DP) :: Delta, DeltaA, xMu_l2, xMu_l4, m2, m4, disc, k_val, E_minus, E_plus, Q
    REAL(DP) :: I0m, I1m, I2m, I0p, I1p, I2p
    REAL(DP) :: ChemPotPairDistinct(nDistinctCasesMax,2)
    INTEGER  :: iCase, iTrueCase

    I0 = 0.0d0
    I1 = 0.0d0
    I2 = 0.0d0

    Delta = SQRT(MAX(E1**2 - 2.0d0*E1*E3*costh + E3**2, 0.0d0))   ! Guo Eq. 4

    ! There are nDistinctCases different cases:
    CALL SetChemPotDistinctCases(Mu_e, Mu_mu, ChemPotPairDistinct)

    DO iCase=1,nDistinctCases
      
      iTrueCase = ChosenToTrueCaseMap(iCase)     
      m2        = MassPairDistinct(iTrueCase,1)
      m4        = MassPairDistinct(iTrueCase,2)

      xMu_l2 = ChemPotPairDistinct(iTrueCase,1)
      xMu_l4 = ChemPotPairDistinct(iTrueCase,2)

      Q = 0.5d0 * (m4**2 - m2**2)

      IF (iTrueCase == i_IMD) THEN

        ! ==============================================================
        ! IMD kinematics (Guo Eqs. 16-19)
        ! ==============================================================
        k_val = Q / (E1 * E3 * (1.0d0 - costh))   ! Guo Eq. 5
        IF (k_val < k0_IMD) CYCLE   ! Theta function, Guo Eq. 16

        DeltaA = SQRT(MAX(E1**2 + 2.0d0*E1*E3*costh + E3**2, 0.0d0))   ! Guo Eq. 17

        disc = (1.0d0 - k_val)**2 - 2.0d0 * m2**2 / (E1 * E3 * (1.0d0 - costh)) ! Guo Eq. 17
        IF (disc < 0.0d0) CYCLE

        E_minus = 0.5d0 * ( (E3 + E1)*(k_val - 1.0d0) - DeltaA*SQRT(disc) )   ! Guo Eq. 17
        E_plus  = 0.5d0 * ( (E3 + E1)*(k_val - 1.0d0) + DeltaA*SQRT(disc) )   ! Guo Eq. 17
        E_minus = MAX(E_minus, m2)
        E_plus  = MAX(E_plus , m2)
        IF (E_minus >= E_plus) CYCLE

        ! I_s = I_s(E_minus) - I_s(E_plus), both evaluated at (E1, -E3) (Guo Eq. 18).
        ! Note these are E1 <-> E3 symmetric, which is why one distinct case
        ! serves both process 33 and process 34.
        CALL Compute_Is_Integrals( E1, -E3, E_minus, T, xMu_l2, xMu_l4, I0m, I1m, I2m )
        CALL Compute_Is_Integrals( E1, -E3, E_plus , T, xMu_l2, xMu_l4, I0p, I1p, I2p )

        I0(iCase) = I0m - I0p
        I1(iCase) = I1m - I1p
        I2(iCase) = I2m - I2p

      ELSE

        IF (Delta < 1.0d-10) RETURN   ! pathological collinear case
        ! --- k (Guo Eq. 5): diverges at forward scattering mu=1 ---
        IF (1.0d0 - costh < 1.0d-10) RETURN

        k_val = Q / (E1 * E3 * (1.0d0 - costh))   ! Guo Eq. 5
        disc = (1.0d0 + k_val)**2 + 2.0d0 * m2**2 / (E1 * E3 * (1.0d0 - costh)) ! Guo Eq. 4
        IF (disc < 0.0d0) CYCLE

        E_minus = 0.5d0 * ( (E3 - E1)*(1.0d0 + k_val) + Delta*SQRT(disc) )   ! Guo Eq. 4
        E_minus = MAX(E_minus, m2)

        CALL Compute_Is_Integrals( E1, E3, E_minus, T, xMu_l2, xMu_l4, I0(iCase), I1(iCase), I2(iCase) )
      
      ENDIF

    ENDDO

  END SUBROUTINE CalculateDistinctI0I1I2

  SUBROUTINE CalculateAllPhi0_collinear( iE1, E1, xT, Mu_e, Mu_mu, Phi0 )

    INTEGER,  INTENT(IN)  :: iE1
    REAL(DP), INTENT(IN)  :: E1, xT, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: Phi0(iProcessMax - iProcessMin + 1)

    INTEGER  :: iE2, iProcess, iCase, iTrueCase, idx
    REAL(DP) :: E2, wE2, xMu_l2, xMu_l4, f2, f4, occ, xnode
    REAL(DP) :: val_H0I, val_H0II
    REAL(DP) :: ChemPotPairDistinct(nDistinctCasesMax,2)
    REAL(DP) :: integrand0I(nDistinctCasesMax), integrand0II(nDistinctCasesMax)

    integrand0I  = 0.0d0
    integrand0II = 0.0d0
    CALL SetChemPotDistinctCases(Mu_e, Mu_mu, ChemPotPairDistinct)

    DO iCase = 1, nDistinctCases
      iTrueCase = ChosenToTrueCaseMap(iCase)

      xMu_l2 = ChemPotPairDistinct(iTrueCase,1)
      xMu_l4 = ChemPotPairDistinct(iTrueCase,2)

      ! ================================================================
      ! PANEL 1: [0, E1] -- finite, use Gauss-Legendre.
      ! Affine-map the preallocated canonical [-1,1] grid to [0,E1].
      ! E2 = (E1/2)*(x_ref + 1),  weight scaled by (E1/2).
      ! ================================================================
      DO iE2 = 1, nE2_Lepton
        E2  = 0.5d0 * E1 * ( xa_GLeg_ref(iE2) + 1.0d0 )
        wE2 = 0.5d0 * E1 * wa_GLeg_ref(iE2)

        CALL CalculateCollinearH0(E1, E2, val_H0I, val_H0II)

        f2 = FD(E2, xMu_l2, xT)
        IF (ABS(xMu_l2 - xMu_l4) < 1.0d-10) THEN
          f4 = f2
        ELSE
          f4 = FD(E2, xMu_l4, xT)
        ENDIF
        occ = f2 * (1.0d0 - f4)

        integrand0I(iTrueCase)  = integrand0I(iTrueCase)  + wE2 * occ * val_H0I
        integrand0II(iTrueCase) = integrand0II(iTrueCase) + wE2 * occ * val_H0II
      ENDDO

      ! ================================================================
      ! PANEL 2: [E1, infinity) -- semi-infinite, use Gauss-Laguerre.
      ! Map raw GLag nodes (weight e^{-x}) via E2 = E1 + x*T.
      ! The built-in e^{-x} must be undone (multiply by e^{+x}) since the
      ! true integrand is f2*(1-f4)*H0, not a pure exponential; xT is the
      ! dE2 = xT dx Jacobian.
      ! ================================================================
      DO iE2 = 1, nE2_Lepton
        xnode = xa_GLag(iE2)
        E2    = E1 + xnode * xT
        wE2   = wa_GLag(iE2) * EXP(xnode) * xT

        CALL CalculateCollinearH0(E1, E2, val_H0I, val_H0II)

        f2 = FD(E2, xMu_l2, xT)
        IF (ABS(xMu_l2 - xMu_l4) < 1.0d-10) THEN
          f4 = f2
        ELSE
          f4 = FD(E2, xMu_l4, xT)
        ENDIF
        occ = f2 * (1.0d0 - f4)

        integrand0I(iTrueCase)  = integrand0I(iTrueCase)  + wE2 * occ * val_H0I
        integrand0II(iTrueCase) = integrand0II(iTrueCase) + wE2 * occ * val_H0II
      ENDDO

    ENDDO

    DO iProcess = iProcessMin, iProcessMax
      idx = iProcess - iProcessMin + 1
      iTrueCase = iDistinctMap(iProcess)

      Phi0(idx) = lam1(idx) * integrand0I(iTrueCase) + &
                  lam2(idx) * integrand0II(iTrueCase)
    END DO

    Phi0(:) = Phi0(:) / E1**4 * 4066026.8989943061

  END SUBROUTINE CalculateAllPhi0_collinear

  SUBROUTINE CalculateAllPhoutGeneral_collinear( iE1, E1, xT, Mu_e, Mu_mu, Phout, nL )
    
    INTEGER , INTENT(IN)  :: iE1
    REAL(DP), INTENT(IN)  :: E1, xT, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: Phout(iProcessMax - iProcessMin + 1, nL)
    INTEGER , INTENT(IN)  :: nL

    REAL(DP) :: ChemPotPairDistinct(nDistinctCasesMax,2)
    REAL(DP) :: integrand0I(nDistinctCasesMax), integrand0II(nDistinctCasesMax)
    REAL(DP) :: Phi0(iProcessMax - iProcessMin + 1)
    INTEGER  :: iL
    
    CALL CalculateAllPhi0_collinear( iE1, E1, xT, Mu_e, Mu_mu, Phi0 )

    Phout(:,1) = Phi0(:) * 0.5d0
    DO iL = 2, nL
      ! Here you have to reconstruct the kernel using eq. 42 from Mezzacappa & Bruenn 1993, 
      ! and then integrate over angle to get all other Phis. My guess is that from Phi0 you can
      ! directly calculate all Phis perhaps without having to do the Theta loop but maybe not.
      ! Don't have time now to implement this. Phi0 should be enough.
      Phout(:,iL) = 0.0d0
    ENDDO

  END SUBROUTINE CalculateAllPhoutGeneral_collinear

END MODULE wlGeneralLeptonScatteringModuleThornadoInterface

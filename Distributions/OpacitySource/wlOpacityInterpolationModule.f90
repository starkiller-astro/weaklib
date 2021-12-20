MODULE wlOpacityInterpolationModule

  USE wlKindModule, ONLY: &
    dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateOpacity_2D1D2D, &
    LinearInterp_Array_Point, &
    GetIndexAndDelta
  USE wlOpacityFieldsModule, ONLY: &
    iNu_e, iNu_e_bar, &
    iHi0, iHii0, iHi1, iHii1, &
    iJi0, iJii0, iJi1, iJii1

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: wlInterpolateOpacity_NES
  PUBLIC :: wlInterpolateOpacity_Pair
  PUBLIC :: wlInterpolateOpacity_Brem

  REAL(dp), PARAMETER :: cv = 0.96d+00 ! weak interaction constant
  REAL(dp), PARAMETER :: ca = 0.50d+00 ! weak interaction constant

CONTAINS


  SUBROUTINE wlInterpolateOpacity_NES &
    ( E, T, Eta, TabE, TabT, TabEta, OS, TabNES, PhiNES )

    ! --- T and TabT in MeV ---

    REAL(dp), INTENT(in)  :: E(:)
    REAL(dp), INTENT(in)  :: T(:)
    REAL(dp), INTENT(in)  :: Eta(:)
    REAL(DP), INTENT(in)  :: TabE(:)
    REAL(DP), INTENT(in)  :: TabT(:)
    REAL(DP), INTENT(in)  :: TabEta(:)
    REAL(DP), INTENT(in)  :: OS(:)
    REAL(DP), INTENT(in)  :: TabNES(:,:,:,:,:)
    REAL(DP), INTENT(out) :: PhiNES(:,:,:,:,:)

    INTEGER  :: nE, nX, iE, iEp, iX, iS, l
    REAL(DP) :: C_i(2), C_ii(2), minBF, maxBF
    REAL(DP), ALLOCATABLE :: NES(:,:,:,:)

    nE = SIZE( E )
    nX = SIZE( T )

    ALLOCATE( NES(nE,nE,4,nX) )

    CALL LogInterpolateOpacity_2D1D2D &
           ( E, T, Eta, TabE, TabT, TabEta, 1, [ 1, 1 ], OS, TabNES, NES )

    C_i (iNu_e)     = ( cv + ca )**2
    C_ii(iNu_e)     = ( cv - ca )**2

    C_i (iNu_e_bar) = ( cv - ca )**2
    C_ii(iNu_e_bar) = ( cv + ca )**2

    DO iX = 1, nX

      DO iE  = 1, nE
      DO iEp = 1, iE

        ! --- Electron Netrinos ---

        ! --- Phi_0 ---

        PhiNES(iEp,iE,1,iNu_e,iX) &
          = C_i(iNu_e) * NES(iEp,iE,iHi0,iX) &
              + C_ii(iNu_e) * NES(iEp,iE,iHii0,iX)

        ! --- Phi_1 ---

        PhiNES(iEp,iE,2,iNu_e,iX) &
          = C_i(iNu_e) * NES(iEp,iE,iHi1,iX) &
              + C_ii(iNu_e) * NES(iEp,iE,iHii1,iX)

        ! --- Electron Antinetrinos ---

        ! --- Phi_0 ---

        PhiNES(iEp,iE,1,iNu_e_bar,iX) &
          = C_i(iNu_e_bar) * NES(iEp,iE,iHi0,iX) &
              + C_ii(iNu_e_bar) * NES(iEp,iE,iHii0,iX)

        ! --- Phi_1 ---

        PhiNES(iEp,iE,2,iNu_e_bar,iX) &
          = C_i(iNu_e_bar) * NES(iEp,iE,iHi1,iX) &
              + C_ii(iNu_e_bar) * NES(iEp,iE,iHii1,iX)

      END DO
      END DO

      ! --- Enforce Detailed Balance ---

      minBF = + HUGE( 1.0_dp )
      maxBF = - HUGE( 1.0_dp )

      DO iS  = iNu_e, iNu_e_bar
      DO l   = 1, 2
      DO iE  = 1,    nE
      DO iEp = iE+1, nE

        PhiNES(iEp,iE,l,iS,iX) &
          = PhiNES(iE,iEp,l,iS,iX) * EXP( ( E(iE) - E(iEp) ) / T(iX) )

        minBF = MIN( minBF, EXP( ( E(iE) - E(iEp) ) / T(iX) ) )
        maxBF = MAX( maxBF, EXP( ( E(iE) - E(iEp) ) / T(iX) ) )

      END DO
      END DO
      END DO
      END DO

    END DO

    PRINT*
    PRINT*, "  NES: min/max BF = ", minBF, maxBF

    DEALLOCATE( NES )

  END SUBROUTINE wlInterpolateOpacity_NES


  SUBROUTINE wlInterpolateOpacity_Pair &
    ( E, T, Eta, TabE, TabT, TabEta, OS, TabPair, PhiPair )

    REAL(dp), INTENT(in)  :: E(:)
    REAL(dp), INTENT(in)  :: T(:)
    REAL(dp), INTENT(in)  :: Eta(:)
    REAL(DP), INTENT(in)  :: TabE(:)
    REAL(DP), INTENT(in)  :: TabT(:)
    REAL(DP), INTENT(in)  :: TabEta(:)
    REAL(DP), INTENT(in)  :: OS(:)
    REAL(DP), INTENT(in)  :: TabPair(:,:,:,:,:)
    REAL(DP), INTENT(out) :: PhiPair(:,:,:,:,:)

    INTEGER  :: nE, nX, iE, iEp, iX, iS, l
    REAL(DP) :: C_i(2), C_ii(2)
    REAL(DP), ALLOCATABLE :: Pair(:,:,:,:)

    nE = SIZE( E )
    nX = SIZE( T )

    ALLOCATE( Pair(nE,nE,4,nX) )

    CALL LogInterpolateOpacity_2D1D2D &
           ( E, T, Eta, TabE, TabT, TabEta, 1, [ 1, 1 ], OS, TabPair, Pair )

    C_i (iNu_e)     = ( cv + ca )**2
    C_ii(iNu_e)     = ( cv - ca )**2

    C_i (iNu_e_bar) = ( cv - ca )**2
    C_ii(iNu_e_bar) = ( cv + ca )**2

    DO iX = 1, nX

      DO iE  = 1, nE
      DO iEp = 1, iE

        ! --- Electron Netrinos ---

        ! --- Phi_0 ---

        PhiPair(iEp,iE,1,iNu_e,iX) &
          = C_i(iNu_e) * Pair(iEp,iE,iJi0,iX) &
              + C_ii(iNu_e) * Pair(iEp,iE,iJii0,iX)

        ! --- Phi_1 ---

        PhiPair(iEp,iE,2,iNu_e,iX) &
          = C_i(iNu_e) * Pair(iEp,iE,iJi1,iX) &
              + C_ii(iNu_e) * Pair(iEp,iE,iJii1,iX)

        ! --- Electron Antinetrinos ---

        ! --- Phi_0 ---

        PhiPair(iEp,iE,1,iNu_e_bar,iX) &
          = C_i(iNu_e_bar) * Pair(iEp,iE,iJi0,iX) &
              + C_ii(iNu_e_bar) * Pair(iEp,iE,iJii0,iX)

        ! --- Phi_1 ---

        PhiPair(iEp,iE,2,iNu_e_bar,iX) &
          = C_i(iNu_e_bar) * Pair(iEp,iE,iJi1,iX) &
              + C_ii(iNu_e_bar) * Pair(iEp,iE,iJii1,iX)

      END DO
      END DO

      ! --- Enforce Symmetry Relations (B85, Eq. (C64)) ---

      DO iE  = 1,    nE
      DO iEp = iE+1, nE

        ! --- Electron Neutrinos ---

        ! --- Phi_0 ---

        PhiPair(iEp,iE,1,iNu_e,iX) &
          = C_i(iNu_e) * Pair(iE,iEp,iJii0,iX) &
              + C_ii(iNu_e) * Pair(iE,iEp,iJi0,iX)

        ! --- Phi_1 ---

        PhiPair(iEp,iE,2,iNu_e,iX) &
          = C_i(iNu_e) * Pair(iE,iEp,iJii1,iX) &
              + C_ii(iNu_e) * Pair(iE,iEp,iJi1,iX)

        ! --- Electron Antineutrinos ---

        ! --- Phi_0 ---

        PhiPair(iEp,iE,1,iNu_e_bar,iX) &
          = C_i(iNu_e_bar) * Pair(iE,iEp,iJii0,iX) &
              + C_ii(iNu_e_bar) * Pair(iE,iEp,iJi0,iX)

        ! --- Phi_1 ---

        PhiPair(iEp,iE,2,iNu_e_bar,iX) &
          = C_i(iNu_e_bar) * Pair(iE,iEp,iJii1,iX) &
              + C_ii(iNu_e_bar) * Pair(iE,iEp,iJi1,iX)

      END DO
      END DO

    END DO

    DEALLOCATE( Pair )

  END SUBROUTINE wlInterpolateOpacity_Pair

  SUBROUTINE wlInterpolateOpacity_Brem &
    ( LogEp, LogE, LogD, LogT, Y, LogEps, LogEs, LogDs, LogTs, Ys, OS, Table, Phi0a_Brem)

    REAL(dp), INTENT(in) :: LogEp(:)
    REAL(dp), INTENT(in) :: LogE(:)
    REAL(dp), INTENT(in) :: LogD(:)
    REAL(dp), INTENT(in) :: LogT(:)
    REAL(dp), INTENT(in) :: Y(:)
    REAL(dp), INTENT(in) :: LogEps(:)
    REAL(dp), INTENT(in) :: LogEs(:)
    REAL(dp), INTENT(in) :: LogDs(:)
    REAL(dp), INTENT(in) :: LogTs(:)
    REAL(dp), INTENT(in) :: Ys(:)

    REAL(dp), INTENT(in) :: OS

    REAL(dp), INTENT(in) :: Table(:,:,:,:,:)

    REAL(dp), INTENT(inout) :: Phi0a_Brem(:,:,:,:,:)

    INTEGER  :: iD, iT, iY
    INTEGER  :: ii, jj, kk, ll, mm
    REAL(dp) :: dD, dT, dY

    DO mm = 1, SIZE(Y)
      CALL GetIndexAndDelta( Y(mm), Ys, iY, dY )
      DO ll = 1, SIZE(LogT)
        CALL GetIndexAndDelta( LogT(ll), LogTs, iT, dT )
        DO kk = 1, SIZE(LogD)
          CALL GetIndexAndDelta( LogD(kk), LogDs, iD, dD )
          DO jj = 1, SIZE(LogE)
            DO ii = 1, SIZE(LogEp)

              Phi0a_Brem(ii,jj,kk,ll,mm) = LinearInterp_Array_Point( ii, jj, iD, iT, iY, dD, dT, dY, OS, Table)

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE wlInterpolateOpacity_Brem

END MODULE wlOpacityInterpolationModule

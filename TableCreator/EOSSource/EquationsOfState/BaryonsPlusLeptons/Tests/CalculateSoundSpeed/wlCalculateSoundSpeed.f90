PROGRAM wlCalculateSoundSpeed

  ! Computes the adiabatic index Gamma1 and the sound speed from
  ! BaryonsPlusPhotonsPlusLeptonsEOS.h5 at every node of the reference table
  ! wl-EOS-SFHo-15-25-50.h5 (no muons) with
  !
  !   Separate : baryon table interpolated and differentiated, leptons and
  !              photons (value and derivatives) from Helmholtz at the exact
  !              point. The baryon pressure is interpolated log-trilinearly,
  !              except in cells where it is <= 0 at some corner: there it is
  !              interpolated linearly in (log D, log T, Yp)
  !   Combined : baryons + leptons + photons summed at the cell corners, then
  !              the log-trilinear interpolant of the total is differentiated
  !
  ! Both are computed by CalculateSoundSpeed in wlSoundSpeedModule.
  !
  ! The results and the reference values are written to SoundSpeed.bin
  ! (stream access, no record markers), in this order:
  !
  !   int32   nD, nT, nYe
  !   float64 D(nD), T(nT), Ye(nYe)
  !   float64 Gamma1(nD,nT,nYe,3)     [ Reference, Separate, Combined ]
  !   float64 Cs    (nD,nT,nYe,3)     [ Reference, Separate, Combined ]
  !   int32   NegativeP(nD,nT,nYe)    1 if the baryon P is <= 0 in the cell
  !
  ! Points outside the new table, with invalid reference values, or where
  ! a method fails (NaN or <= 0) are NaN. Plot with plot_sound_speed.py.

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: &
    IEEE_VALUE, &
    IEEE_QUIET_NAN, &
    IEEE_IS_NAN
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    INT32
  USE wlKindModule, ONLY: &
    dp
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType, &
    EquationOfStateCompOSETableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlHelmIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF
  USE wlInterpolationUtilitiesModule, ONLY: &
    GetIndexAndDelta_Lin, &
    GetIndexAndDelta_Log
  USE wlSoundSpeedModule, ONLY: &
    CalculateSoundSpeed

  IMPLICIT NONE

  INTEGER,  PARAMETER :: iReference = 1, iSeparate = 2, iCombined = 3
  CHARACTER(8), PARAMETER :: MethodNames(3) = [ 'Ref.    ', 'Separate', 'Combined' ]

  CHARACTER(len=128), PARAMETER :: &
    BaryonEOSTableName    = 'BaryonsPlusPhotonsPlusLeptonsEOS.h5', &
    ReferenceEOSTableName = 'wl-EOS-SFHo-15-25-50.h5', &
    OutputFileName        = 'SoundSpeed.bin'

  TYPE(EquationOfStateCompOSETableType) :: EOSBaryonTable
  TYPE(EquationOfStateTableType)        :: EOSReferenceTable
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(HelmTableType) :: HelmTableMuons

  ! --- New table (baryons only): P stored linearly, E and S as log10(X+OS)
  REAL(dp), ALLOCATABLE :: D_T(:), T_T(:), Yp_T(:)
  REAL(dp), ALLOCATABLE :: P_T(:,:,:), E_T(:,:,:), S_T(:,:,:)
  REAL(dp) :: OS_P, OS_E, OS_S

  ! --- Reference table: all variables stored as log10(X+OS)
  REAL(dp), ALLOCATABLE :: D_R(:), T_R(:), Ye_R(:)
  INTEGER  :: iP_R, iG_R
  REAL(dp) :: OS_P_R, OS_G_R

  REAL(dp),       ALLOCATABLE :: Gamma(:,:,:,:), Cs(:,:,:,:)
  INTEGER(INT32), ALLOCATABLE :: NegativeP(:,:,:)

  INTEGER  :: iD, iT, iYe, iMethod, iUnit
  INTEGER  :: nPoints, nInvalidRef, nInvalid(3), nSeparateCloser
  INTEGER  :: iD_T, iT_T, iYp_T
  REAL(dp) :: dD_T, dT_T, dYp_T
  REAL(dp) :: D, T, Ye, Ym, P_Ref, NaN
  REAL(dp) :: ErrG(3), ErrCs(3)
  REAL(dp) :: MaxErrG(3), SumErrG(3), MaxErrCs(3), SumErrCs(3)

  NaN = IEEE_VALUE( 0.0_dp, IEEE_QUIET_NAN )

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSBaryonTable, BaryonEOSTableName )
  CALL ReadHelmholtzTableHDF( HelmTableElectrons, BaryonEOSTableName, "HelmTableElectrons" )
  CALL ReadHelmholtzTableHDF( HelmTableMuons    , BaryonEOSTableName, "HelmTableMuons"     )

  CALL ReadEquationOfStateTableHDF( EOSReferenceTable, ReferenceEOSTableName )

  CALL FinalizeHDF( )

  ! --- Unpack new table ---

  ASSOCIATE( TState => EOSBaryonTable % TS, DV => EOSBaryonTable % DV )

  ALLOCATE( D_T (TState % nPoints(TState % Indices % iRho)) )
  ALLOCATE( T_T (TState % nPoints(TState % Indices % iT  )) )
  ALLOCATE( Yp_T(TState % nPoints(TState % Indices % iYe )) )
  D_T  = TState % States(TState % Indices % iRho) % Values
  T_T  = TState % States(TState % Indices % iT  ) % Values
  Yp_T = TState % States(TState % Indices % iYe ) % Values

  OS_P = DV % Offsets(DV % Indices % iPressure)
  OS_E = DV % Offsets(DV % Indices % iInternalEnergyDensity)
  OS_S = DV % Offsets(DV % Indices % iEntropyPerBaryon)

  ALLOCATE( P_T(SIZE(D_T),SIZE(T_T),SIZE(Yp_T)) )
  ALLOCATE( E_T(SIZE(D_T),SIZE(T_T),SIZE(Yp_T)) )
  ALLOCATE( S_T(SIZE(D_T),SIZE(T_T),SIZE(Yp_T)) )
  P_T = DV % Variables(DV % Indices % iPressure             ) % Values
  E_T = DV % Variables(DV % Indices % iInternalEnergyDensity) % Values
  S_T = DV % Variables(DV % Indices % iEntropyPerBaryon     ) % Values

  END ASSOCIATE

  ! --- Unpack reference table ---

  ASSOCIATE( TState => EOSReferenceTable % TS, DV => EOSReferenceTable % DV )

  ALLOCATE( D_R (TState % nPoints(TState % Indices % iRho)) )
  ALLOCATE( T_R (TState % nPoints(TState % Indices % iT  )) )
  ALLOCATE( Ye_R(TState % nPoints(TState % Indices % iYe )) )
  D_R  = TState % States(TState % Indices % iRho) % Values
  T_R  = TState % States(TState % Indices % iT  ) % Values
  Ye_R = TState % States(TState % Indices % iYe ) % Values

  iP_R   = DV % Indices % iPressure
  iG_R   = DV % Indices % iGamma1
  OS_P_R = DV % Offsets(iP_R)
  OS_G_R = DV % Offsets(iG_R)

  END ASSOCIATE

  ALLOCATE( Gamma    (SIZE(D_R),SIZE(T_R),SIZE(Ye_R),3) )
  ALLOCATE( Cs       (SIZE(D_R),SIZE(T_R),SIZE(Ye_R),3) )
  ALLOCATE( NegativeP(SIZE(D_R),SIZE(T_R),SIZE(Ye_R)  ) )
  Gamma     = NaN
  Cs        = NaN
  NegativeP = 0

  ! --- Evaluate at every reference node ---

  nPoints         = 0
  nInvalidRef     = 0
  nInvalid        = 0
  nSeparateCloser = 0
  MaxErrG  = 0.0_dp; SumErrG  = 0.0_dp
  MaxErrCs = 0.0_dp; SumErrCs = 0.0_dp

  Ym = 0.0_dp ! The reference table has no muons

  DO iYe = 1, SIZE(Ye_R)
    DO iT = 1, SIZE(T_R)
      DO iD = 1, SIZE(D_R)

        D  = D_R (iD)
        T  = T_R (iT)
        Ye = Ye_R(iYe)

        IF ( D  < D_T (1) .OR. D  > D_T (SIZE(D_T )) .OR. &
             T  < T_T (1) .OR. T  > T_T (SIZE(T_T )) .OR. &
             Ye < Yp_T(1) .OR. Ye > Yp_T(SIZE(Yp_T)) ) CYCLE

        Gamma(iD,iT,iYe,iReference) &
          = 10.0_dp**EOSReferenceTable % DV % Variables(iG_R) % Values(iD,iT,iYe) - OS_G_R
        P_Ref &
          = 10.0_dp**EOSReferenceTable % DV % Variables(iP_R) % Values(iD,iT,iYe) - OS_P_R

        ! --- Placeholder / unphysical reference values ---
        IF ( Gamma(iD,iT,iYe,iReference) <= 1.0d-20 .OR. P_Ref <= 0.0_dp ) THEN
          Gamma(iD,iT,iYe,iReference) = NaN
          nInvalidRef = nInvalidRef + 1
          CYCLE
        END IF

        Cs(iD,iT,iYe,iReference) = SQRT( Gamma(iD,iT,iYe,iReference) * P_Ref / D )

        DO iMethod = iSeparate, iCombined
          CALL CalculateSoundSpeed &
                 ( D, T, Ye, Ym, D_T, T_T, Yp_T, P_T, OS_P, E_T, OS_E, S_T, OS_S, &
                   HelmTableElectrons, HelmTableMuons, &
                   Gamma(iD,iT,iYe,iMethod), Cs(iD,iT,iYe,iMethod), iMethod == iSeparate )
        END DO

        ! --- Flag cells where the baryon pressure is <= 0 at some corner ---
        CALL GetIndexAndDelta_Log( D , D_T , iD_T , dD_T  )
        CALL GetIndexAndDelta_Log( T , T_T , iT_T , dT_T  )
        CALL GetIndexAndDelta_Lin( Ye, Yp_T, iYp_T, dYp_T )
        IF ( MINVAL( P_T(iD_T:iD_T+1,iT_T:iT_T+1,iYp_T:iYp_T+1) ) <= 0.0_dp ) &
          NegativeP(iD,iT,iYe) = 1

        ! --- Failed results (NaN or non-positive) are flagged as NaN ---
        DO iMethod = iSeparate, iCombined
          IF ( .NOT. ( Gamma(iD,iT,iYe,iMethod) > 0.0_dp ) .OR. &
               .NOT. ( Cs   (iD,iT,iYe,iMethod) > 0.0_dp ) ) THEN
            Gamma(iD,iT,iYe,iMethod) = NaN
            Cs   (iD,iT,iYe,iMethod) = NaN
            nInvalid(iMethod) = nInvalid(iMethod) + 1
          END IF
        END DO
        IF ( ANY( IEEE_IS_NAN( Gamma(iD,iT,iYe,iSeparate:iCombined) ) ) ) CYCLE

        nPoints = nPoints + 1

        ErrG  = ABS( Gamma(iD,iT,iYe,:) - Gamma(iD,iT,iYe,iReference) ) &
                / Gamma(iD,iT,iYe,iReference)
        ErrCs = ABS( Cs   (iD,iT,iYe,:) - Cs   (iD,iT,iYe,iReference) ) &
                / Cs   (iD,iT,iYe,iReference)

        MaxErrG  = MAX( MaxErrG , ErrG  )
        MaxErrCs = MAX( MaxErrCs, ErrCs )
        SumErrG  = SumErrG  + ErrG
        SumErrCs = SumErrCs + ErrCs

        IF ( ErrG(iSeparate) < ErrG(iCombined) ) nSeparateCloser = nSeparateCloser + 1

      END DO
    END DO
  END DO

  ! --- Output ---

  OPEN( NEWUNIT = iUnit, FILE = TRIM( OutputFileName ), ACCESS = 'STREAM', &
        FORM = 'UNFORMATTED', STATUS = 'REPLACE' )
  WRITE( iUnit ) INT( SIZE(D_R), INT32 ), INT( SIZE(T_R), INT32 ), INT( SIZE(Ye_R), INT32 )
  WRITE( iUnit ) D_R, T_R, Ye_R
  WRITE( iUnit ) Gamma
  WRITE( iUnit ) Cs
  WRITE( iUnit ) NegativeP
  CLOSE( iUnit )

  ! --- Summary ---

  WRITE(*,*)
  WRITE(*,'(A,I8,A,I6,A)') ' Compared at ', nPoints, ' reference nodes (', &
    nInvalidRef, ' skipped: invalid reference Gamma1 or P)'
  WRITE(*,'(A,I8,A)') ' Nodes in cells with baryon P <= 0: ', &
    COUNT( NegativeP == 1 ), ' (Separate interpolates P linearly there)'
  WRITE(*,'(A,2I8)') ' Failed results (NaN or <= 0), Separate / Combined: ', &
    nInvalid(iSeparate:iCombined)
  WRITE(*,*)
  WRITE(*,'(A)') ' Relative error w.r.t. reference'
  WRITE(*,'(A4,A10,4A14)') '', 'Method', 'max Gamma1', 'mean Gamma1', 'max cs', 'mean cs'
  DO iMethod = iSeparate, iCombined
    WRITE(*,'(A4,A10,4ES14.3)') '', MethodNames(iMethod), &
      MaxErrG (iMethod), SumErrG (iMethod) / DBLE( MAX( nPoints, 1 ) ), &
      MaxErrCs(iMethod), SumErrCs(iMethod) / DBLE( MAX( nPoints, 1 ) )
  END DO
  WRITE(*,'(A,F6.1,A)') ' Separate closer to the reference Gamma1 at ', &
    100.0_dp * DBLE( nSeparateCloser ) / DBLE( MAX( nPoints, 1 ) ), ' % of the points'
  WRITE(*,*)
  WRITE(*,'(A)') ' Results written to ' // TRIM( OutputFileName )

  DEALLOCATE( D_T, T_T, Yp_T, P_T, E_T, S_T )
  DEALLOCATE( D_R, T_R, Ye_R, Gamma, Cs, NegativeP )

END PROGRAM wlCalculateSoundSpeed

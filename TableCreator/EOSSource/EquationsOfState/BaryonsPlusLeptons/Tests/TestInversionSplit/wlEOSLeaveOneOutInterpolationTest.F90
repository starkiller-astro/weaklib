PROGRAM wlEOSLeaveOneOutInterpolationTest

  ! Leave-one-out accuracy test of the two ways of interpolating the
  ! baryon + lepton + photon EOS between temperature nodes:
  !
  !   Combined : log10( X_bary + X_LeptPhot ) built at the nodes and
  !              interpolated linearly in log T (same as the legacy 3D table)
  !   Separate : log10( X_bary ) interpolated linearly in log T, plus
  !              X_LeptPhot evaluated exactly at the point (Helmholtz)
  !
  ! Every other temperature node is dropped. At each dropped node T(i) both
  ! methods interpolate from T(i-1) and T(i+1), and are compared with the
  ! exact value there: baryon table node value + Helmholtz at T(i).
  ! Density and Yp stay on the nodes, so only the T interpolation is tested.
  ! The T spacing is doubled, so absolute errors are ~4x those of the full
  ! table; the comparison between the two methods is what matters.

  USE wlKindModule, ONLY: &
    dp
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateCompOSETableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlHelmIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: &
    GetPhotonLeptonGasEOS

  IMPLICIT NONE

  INTEGER, PARAMETER :: iVar_E = 1, iVar_P = 2, iVar_S = 3
  INTEGER, PARAMETER :: iCombined = 1, iSeparate = 2
  INTEGER, PARAMETER :: nBins = 6
  CHARACTER(1), PARAMETER :: VarNames(3) = [ 'E', 'P', 'S' ]

  ! --- Muon fraction: Ym = Ym_over_Yp * Yp
  REAL(dp), PARAMETER :: Ym_over_Yp = 0.0d0

  CHARACTER(len=128), PARAMETER :: &
    BaryonPlusHelmTableName = 'BaryonsPlusPhotonsPlusLeptonsEOS.h5'

  TYPE(EquationOfStateCompOSETableType) :: EOSBaryonTable
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(HelmTableType) :: HelmTableMuons

  REAL(dp), ALLOCATABLE :: Ds(:), Ts(:), Yps(:)
  REAL(dp), ALLOCATABLE :: Xs_bary(:,:,:,:)
  REAL(dp) :: OS(3)
  INTEGER  :: iDV(3)
  INTEGER  :: nD, nT, nYp

  INTEGER  :: iD, iT, iYp, iVar, iL, iBin, nPoints
  REAL(dp) :: Ye, Ym, w, LocalOffset
  REAL(dp) :: P_LeptPhot, E_LeptPhot, S_LeptPhot
  REAL(dp) :: X_LeptPhot(-1:1,3), X_bary(-1:1)
  REAL(dp) :: X_Exact, X_Interp(2), RelErr(2)
  REAL(dp) :: MaxErr(3,2), SumErr(3,2), MaxErrBin(nBins,3,2)
  INTEGER  :: nSeparateCloser(3)

  CALL InitializeHDF( )

  CALL ReadHelmholtzTableHDF( HelmTableElectrons, BaryonPlusHelmTableName, "HelmTableElectrons" )
  CALL ReadHelmholtzTableHDF( HelmTableMuons, BaryonPlusHelmTableName, &
      "HelmTableMuons", eos_MinD = 1.0d7 )
  CALL ReadEquationOfStateTableHDF( EOSBaryonTable, BaryonPlusHelmTableName )

  CALL FinalizeHDF( )

  ASSOCIATE( TState => EOSBaryonTable % TS )
  nD  = TState % nPoints(TState % Indices % iRho)
  nT  = TState % nPoints(TState % Indices % iT  )
  nYp = TState % nPoints(TState % Indices % iYe )
  ALLOCATE( Ds(nD), Ts(nT), Yps(nYp) )
  Ds  = TState % States(TState % Indices % iRho) % Values
  Ts  = TState % States(TState % Indices % iT  ) % Values
  Yps = TState % States(TState % Indices % iYe ) % Values
  END ASSOCIATE

  ! --- E and S are stored as log10(X + OS), P linearly ---
  iDV(iVar_E) = EOSBaryonTable % DV % Indices % iInternalEnergyDensity
  iDV(iVar_P) = EOSBaryonTable % DV % Indices % iPressure
  iDV(iVar_S) = EOSBaryonTable % DV % Indices % iEntropyPerBaryon

  ALLOCATE( Xs_bary(nD,nT,nYp,3) )
  DO iVar = 1, 3
    OS(iVar) = EOSBaryonTable % DV % Offsets(iDV(iVar))
    Xs_bary(:,:,:,iVar) = EOSBaryonTable % DV % Variables(iDV(iVar)) % Values
  END DO

  MaxErr          = 0.0_dp
  SumErr          = 0.0_dp
  MaxErrBin       = 0.0_dp
  nSeparateCloser = 0
  nPoints         = 0

  DO iYp = 1, nYp
    DO iD = 1, nD
      DO iT = 2, nT - 1, 2

        Ym = Ym_over_Yp * Yps(iYp)
        Ye = Yps(iYp) - Ym

        ! --- Position of T(i) in [ T(i-1), T(i+1) ] in log T ---
        w = LOG( Ts(iT) / Ts(iT-1) ) / LOG( Ts(iT+1) / Ts(iT-1) )

        ! --- Leptons and photons at T(i-1), T(i), T(i+1) ---
        DO iL = -1, 1
          CALL GetPhotonLeptonGasEOS &
                 ( Ds(iD), Ts(iT+iL), Ye, Ym, HelmTableElectrons, HelmTableMuons, &
                   P_LeptPhot, E_LeptPhot, S_LeptPhot )
          X_LeptPhot(iL,iVar_E) = E_LeptPhot
          X_LeptPhot(iL,iVar_P) = P_LeptPhot
          X_LeptPhot(iL,iVar_S) = S_LeptPhot
        END DO

        iBin = MIN( nBins, MAX( 1, INT( ( LOG10( Ts(iT) ) - 9.0_dp ) / 0.5_dp ) + 1 ) )

        DO iVar = 1, 3

          ! --- Linear baryon values at T(i-1), T(i), T(i+1) ---
          IF ( iVar == iVar_P ) THEN
            X_bary = Xs_bary(iD,iT-1:iT+1,iYp,iVar)
          ELSE
            X_bary = 10.0_dp**Xs_bary(iD,iT-1:iT+1,iYp,iVar) - OS(iVar)
          END IF

          X_Exact = X_bary(0) + X_LeptPhot(0,iVar)

          ! --- Combined: log of the total at the neighbours, interpolated ---
          IF ( iVar == iVar_P ) THEN
            X_Interp(iCombined) &
              = 10.0_dp**( ( 1.0_dp - w ) * LOG10( X_bary(-1) + X_LeptPhot(-1,iVar) ) &
                         +            w   * LOG10( X_bary(+1) + X_LeptPhot(+1,iVar) ) )
          ELSE
            X_Interp(iCombined) &
              = 10.0_dp**( ( 1.0_dp - w ) * LOG10( X_bary(-1) + OS(iVar) + X_LeptPhot(-1,iVar) ) &
                         +            w   * LOG10( X_bary(+1) + OS(iVar) + X_LeptPhot(+1,iVar) ) ) &
                - OS(iVar)
          END IF

          ! --- Separate: baryons interpolated, leptons/photons exact ---
          IF ( iVar == iVar_P ) THEN
            ! Local offset as in the forward pressure interpolation
            LocalOffset = MIN( X_bary(-1), X_bary(+1) )
            IF ( LocalOffset .lt. 0.0_dp ) THEN
              LocalOffset = -1.1d0*LocalOffset
            ELSE
              LocalOffset = 0.0_dp
            ENDIF
            X_Interp(iSeparate) &
              = 10.0_dp**( ( 1.0_dp - w ) * LOG10( X_bary(-1) + LocalOffset ) &
                         +            w   * LOG10( X_bary(+1) + LocalOffset ) ) &
                - LocalOffset + X_LeptPhot(0,iVar)
          ELSE
            X_Interp(iSeparate) &
              = 10.0_dp**( ( 1.0_dp - w ) * Xs_bary(iD,iT-1,iYp,iVar) &
                         +            w   * Xs_bary(iD,iT+1,iYp,iVar) ) &
                - OS(iVar) + X_LeptPhot(0,iVar)
          END IF

          RelErr = ABS( X_Interp - X_Exact ) / ABS( X_Exact )

          MaxErr(iVar,:) = MAX( MaxErr(iVar,:), RelErr )
          SumErr(iVar,:) = SumErr(iVar,:) + RelErr
          MaxErrBin(iBin,iVar,:) = MAX( MaxErrBin(iBin,iVar,:), RelErr )
          IF ( RelErr(iSeparate) < RelErr(iCombined) ) &
            nSeparateCloser(iVar) = nSeparateCloser(iVar) + 1

        END DO

        nPoints = nPoints + 1

      END DO
    END DO
  END DO

  WRITE(*,*)
  WRITE(*,'(A,I8,A)') ' Leave-one-out in T at ', nPoints, &
    ' dropped nodes: relative error w.r.t. exact value'
  WRITE(*,*)
  WRITE(*,'(A4,A3,4A15,A18)') '', 'Var', 'Combined max', 'Combined mean', &
    'Separate max', 'Separate mean', 'Separate closer'
  DO iVar = 1, 3
    WRITE(*,'(A4,A3,4ES15.3,F16.1,A2)') '', VarNames(iVar), &
      MaxErr(iVar,iCombined), SumErr(iVar,iCombined) / DBLE( nPoints ), &
      MaxErr(iVar,iSeparate), SumErr(iVar,iSeparate) / DBLE( nPoints ), &
      100.0_dp * DBLE( nSeparateCloser(iVar) ) / DBLE( nPoints ), ' %'
  END DO

  WRITE(*,*)
  WRITE(*,'(A)') ' Max relative error per log10 T bin'
  WRITE(*,'(A4,A3,A10,6A10)') '', 'Var', 'Method', &
    '[9,9.5)', '[9.5,10)', '[10,10.5)', '[10.5,11)', '[11,11.5)', '[11.5,..)'
  DO iVar = 1, 3
    WRITE(*,'(A4,A3,A10,6ES10.2)') '', VarNames(iVar), 'Combined', MaxErrBin(:,iVar,iCombined)
    WRITE(*,'(A4,A3,A10,6ES10.2)') '', VarNames(iVar), 'Separate', MaxErrBin(:,iVar,iSeparate)
  END DO

  DEALLOCATE( Ds, Ts, Yps, Xs_bary )

END PROGRAM wlEOSLeaveOneOutInterpolationTest

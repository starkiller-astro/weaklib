PROGRAM wlTestMuonChargedCurrent2D

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: mn, mp, mmu
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D, Opacity_CC_2D_PrecomputedNodes, gauleg
  USE wlSemiLeptonicOpacityIntegrationModule2D, ONLY: &
    CC2D_CallData, CC2D_RowData, &
    CC2D_PrepareCall, CC2D_RowSetup, CC2D_NodeKernel, &
    CC2D_ApplyPrefactor

  IMPLICIT NONE

  INTEGER, PARAMETER :: nStates = 2, nEnergies = 5
  INTEGER, PARAMETER :: nQ = 32, nQReference = 96
  INTEGER, PARAMETER :: nOrders = 5
  INTEGER, PARAMETER :: QuadratureOrder(nOrders) = [ 16, 24, 32, 48, 64 ]
  REAL(dp), PARAMETER :: RelativeTolerance = 1.0d-2
  REAL(dp), PARAMETER :: TinyRate = 1.0d-90
  REAL(dp), PARAMETER :: Energy(nEnergies) = &
    [ 1.0_dp, 10.0_dp, 30.0_dp, 80.0_dp, 180.0_dp ]
  REAL(dp), PARAMETER :: Temperature(nStates) = [ 38.3_dp, 15.2_dp ]
  REAL(dp), PARAMETER :: MuMuon(nStates) = [ 64.1_dp, 37.3_dp ]
  REAL(dp), PARAMETER :: MuNeutron(nStates) = [ 886.0_dp, 912.7_dp ]
  REAL(dp), PARAMETER :: MuProton(nStates) = [ 800.7_dp, 875.4_dp ]
  REAL(dp), PARAMETER :: UNeutron(nStates) = [ -24.9_dp, -6.1_dp ]
  REAL(dp), PARAMETER :: UProton(nStates) = [ -42.9_dp, -8.9_dp ]

  REAL(dp) :: x(nQ), w(nQ), xReference(nQReference), wReference(nQReference)
  REAL(dp) :: xOrder(nQReference,nOrders), wOrder(nQReference,nOrders)
  REAL(dp) :: Direct, Precomputed, TwoPhase, Reference
  REAL(dp) :: NodeContribution, RelativeError, MaxRelativeError
  REAL(dp) :: MaxRelativeErrorByOrder(nOrders)
  REAL(dp) :: MaxOpacityByReaction(3)
  INTEGER :: iState, iEnergy, iReaction, iOrder, iQ1, iQ2
  INTEGER :: IncludeCorrections(3)
  TYPE(CC2D_CallData) :: CallData
  TYPE(CC2D_RowData)  :: RowData

  CALL gauleg( 0.0_dp, 1.0_dp, x, w, nQ )
  CALL gauleg( 0.0_dp, 1.0_dp, xReference, wReference, nQReference )
  DO iOrder = 1, nOrders
    CALL gauleg &
           ( 0.0_dp, 1.0_dp, &
             xOrder(1:QuadratureOrder(iOrder),iOrder), &
             wOrder(1:QuadratureOrder(iOrder),iOrder), &
             QuadratureOrder(iOrder) )
  END DO

  MaxRelativeError = 0.0_dp
  MaxRelativeErrorByOrder = 0.0_dp
  MaxOpacityByReaction = 0.0_dp
  IncludeCorrections = [ 1, 1, 0 ]
  DO iState = 1, nStates
    DO iEnergy = 1, nEnergies
      DO iReaction = 1, 3
        CALL Opacity_CC_2D &
               ( 2, iReaction, Energy(iEnergy), Direct, &
                 Temperature(iState), MuMuon(iState), &
                 MuNeutron(iState), MuProton(iState), mmu, mn, mp, &
                 UNeutron(iState), UProton(iState), nQ )
        CALL Opacity_CC_2D_PrecomputedNodes &
               ( 2, iReaction, Energy(iEnergy), Precomputed, &
                 Temperature(iState), MuMuon(iState), &
                 MuNeutron(iState), MuProton(iState), mmu, mn, mp, &
                 UNeutron(iState), UProton(iState), nQ, x, w )
        CALL Opacity_CC_2D_PrecomputedNodes &
               ( 2, iReaction, Energy(iEnergy), Reference, &
                 Temperature(iState), MuMuon(iState), &
                 MuNeutron(iState), MuProton(iState), mmu, mn, mp, &
                 UNeutron(iState), UProton(iState), nQReference, &
                 xReference, wReference )

        CALL CC2D_PrepareCall &
               ( iReaction, Energy(iEnergy), Temperature(iState), &
                 MuMuon(iState), MuNeutron(iState), MuProton(iState), &
                 mmu, mn, mp, UNeutron(iState), UProton(iState), &
                 IncludeCorrections, CallData )
        TwoPhase = 0.0_dp
        DO iQ1 = 1, nQ
          CALL CC2D_RowSetup( CallData, x(iQ1), w(iQ1), RowData )
          DO iQ2 = 1, nQ
            CALL CC2D_NodeKernel &
                   ( CallData, RowData, x(iQ2), w(iQ2), NodeContribution )
            TwoPhase = TwoPhase + NodeContribution
          END DO
        END DO
        CALL CC2D_ApplyPrefactor( CallData, TwoPhase )

        IF( Direct /= Precomputed ) &
          ERROR STOP 'Direct and precomputed-node results differ.'
        IF( Precomputed /= TwoPhase ) &
          ERROR STOP 'Monolithic and two-phase results differ.'
        IF( Direct < 0.0_dp .OR. Direct /= Direct ) &
          ERROR STOP 'Muon charged-current opacity is invalid.'
        MaxOpacityByReaction(iReaction) &
          = MAX( MaxOpacityByReaction(iReaction), Direct )

        IF( Reference > TinyRate )THEN
          RelativeError = ABS( Precomputed - Reference ) / Reference
          MaxRelativeError = MAX( MaxRelativeError, RelativeError )
          DO iOrder = 1, nOrders
            CALL Opacity_CC_2D_PrecomputedNodes &
                   ( 2, iReaction, Energy(iEnergy), Direct, &
                     Temperature(iState), MuMuon(iState), &
                     MuNeutron(iState), MuProton(iState), mmu, mn, mp, &
                     UNeutron(iState), UProton(iState), &
                     QuadratureOrder(iOrder), &
                     xOrder(1:QuadratureOrder(iOrder),iOrder), &
                     wOrder(1:QuadratureOrder(iOrder),iOrder) )
            MaxRelativeErrorByOrder(iOrder) &
              = MAX( MaxRelativeErrorByOrder(iOrder), &
                     ABS( Direct - Reference ) / Reference )
          END DO
        END IF
      END DO
    END DO
  END DO

  WRITE(*,'(A,ES12.4E3)') &
    'maximum nQ=32 relative error against nQ=96: ', MaxRelativeError
  WRITE(*,'(A)') 'quadrature order and maximum relative error:'
  DO iOrder = 1, nOrders
    WRITE(*,'(I5,2X,ES12.4E3)') &
      QuadratureOrder(iOrder), MaxRelativeErrorByOrder(iOrder)
  END DO
  WRITE(*,'(A,3(1X,ES12.4E3))') &
    'maximum opacity by reaction:', MaxOpacityByReaction
  IF( MaxRelativeError > RelativeTolerance ) &
    ERROR STOP 'Muon charged-current quadrature did not converge.'
  IF( ANY( MaxOpacityByReaction(1:2) <= TinyRate ) ) &
    ERROR STOP 'A muon charged-current capture reaction was not exercised.'
  WRITE(*,'(A)') 'muon charged-current 2D test: PASS'

END PROGRAM wlTestMuonChargedCurrent2D

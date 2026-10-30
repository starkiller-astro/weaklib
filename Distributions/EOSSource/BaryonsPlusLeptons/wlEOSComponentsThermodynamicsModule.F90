MODULE wlEOSComponentsThermodynamicsModule

  USE wlKindModule, ONLY: dp
  USE wlLeptonEOSTableModule, ONLY: HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: GetPhotonLeptonGasEOS
  USE wlInterpolationUtilitiesModule, ONLY: &
    GetIndexAndDelta_Log, GetIndexAndDelta_Lin, &
    LinearInterp3D_3DArray_Point, LinearInterpDeriv_Array_Point

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeTotalEPS_CS
  PUBLIC :: EvaluatePhotonLeptonEOSPlane
  PUBLIC :: StoreCombinedEOSPlanes
  PUBLIC :: CombinedEOSPlaneCacheType

  ! Final inversion planes for one subsequent EOS query using the same
  ! unchanged Helm tables. No pointers, allocation or persistent device state.
  TYPE :: CombinedEOSPlaneCacheType
    LOGICAL :: Valid
    INTEGER :: iD, iT, iYp
    REAL(dp) :: D, Ye, Ym
    REAL(dp) :: DensityNodes(2), TemperatureNodes(2), YpNodes(2)
    REAL(dp) :: LeptonEPS(2,2,2,3)
  END TYPE CombinedEOSPlaneCacheType

CONTAINS

  SUBROUTINE EvaluatePhotonLeptonEOSPlane &
    ( DensityNodes, Temperature, YeNodes, YmNodes, &
      HelmTableElectrons, HelmTableMuons, LeptonEPS )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: DensityNodes(2), Temperature, YeNodes(2), YmNodes(2)
    TYPE(HelmTableType), INTENT(in) :: HelmTableElectrons, HelmTableMuons
    REAL(dp), INTENT(out) :: LeptonEPS(2,2,3)

    INTEGER :: iD, iY

    DO iD = 1, 2
      DO iY = 1, 2
        CALL GetPhotonLeptonGasEOS &
               ( DensityNodes(iD), Temperature, YeNodes(iY), YmNodes(iY), &
                 HelmTableElectrons, HelmTableMuons, &
                 LeptonEPS(iD,iY,2), LeptonEPS(iD,iY,1), LeptonEPS(iD,iY,3) )
      END DO
    END DO

  END SUBROUTINE EvaluatePhotonLeptonEOSPlane


  SUBROUTINE StoreCombinedEOSPlanes &
    ( D, Ye, Ym, Ds, Ts, Yps, iD, iTa, iTb, iYp, PlaneA, PlaneB, Cache )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: D, Ye, Ym, Ds(:), Ts(:), Yps(:)
    INTEGER, INTENT(in) :: iD, iTa, iTb, iYp
    REAL(dp), INTENT(in) :: PlaneA(2,2,3), PlaneB(2,2,3)
    TYPE(CombinedEOSPlaneCacheType), INTENT(out) :: Cache

    Cache%Valid = iTb == iTa+1
    IF( .NOT. Cache%Valid ) RETURN
    Cache%iD = iD
    Cache%iT = iTa
    Cache%iYp = iYp
    Cache%D = D
    Cache%Ye = Ye
    Cache%Ym = Ym
    Cache%DensityNodes = Ds(iD:iD+1)
    Cache%TemperatureNodes = Ts(iTa:iTb)
    Cache%YpNodes = Yps(iYp:iYp+1)
    Cache%LeptonEPS(:,1,:,:) = PlaneA
    Cache%LeptonEPS(:,2,:,:) = PlaneB

  END SUBROUTINE StoreCombinedEOSPlanes


  ! Combined EOS on the baryon grid: E [erg/g], P [dyn/cm^2], S [k_B/baryon],
  ! Gamma [dimensionless], Cs [cm/s]. D/T derivatives keep Ye and Ym fixed.
  SUBROUTINE ComputeTotalEPS_CS &
    ( D, T, Ye, Ym, Ds, Ts, Yps, Es, Ps, Ss, OS_E, OS_P, OS_S, &
      HelmTableElectrons, HelmTableMuons, E, P, S, &
      Gamma_Option, Cs_Option, dPdD_Option, dPdT_Option, &
      dEdD_Option, dEdT_Option, PlaneCache_Option )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: D, T, Ye, Ym, Ds(:), Ts(:), Yps(:)
    REAL(dp), INTENT(in) :: Es(:,:,:), Ps(:,:,:), Ss(:,:,:), OS_E, OS_P, OS_S
    TYPE(HelmTableType), INTENT(in) :: HelmTableElectrons, HelmTableMuons
    REAL(dp), INTENT(out) :: E, P, S
    REAL(dp), INTENT(out), OPTIONAL :: Gamma_Option, Cs_Option
    REAL(dp), INTENT(out), OPTIONAL :: dPdD_Option, dPdT_Option
    REAL(dp), INTENT(out), OPTIONAL :: dEdD_Option, dEdT_Option
    TYPE(CombinedEOSPlaneCacheType), INTENT(in), OPTIONAL :: PlaneCache_Option

    REAL(dp) :: Yp, YeNodes(2), YmNodes(2), dD, dT, dYp, aD, aT, aYp
    REAL(dp) :: LeptonEPS(2,2,2,3), LogEPS(2,2,2,3)
    REAL(dp) :: dPdD, dPdT, dEdD, dEdT, Gamma
    INTEGER :: iD, iT, iYp, iPlane
    LOGICAL :: ReusePlanes, NeedDerivatives

    Yp = Ye+Ym
    CALL GetIndexAndDelta_Log( D, Ds, iD, dD )
    CALL GetIndexAndDelta_Log( T, Ts, iT, dT )
    CALL GetIndexAndDelta_Lin( Yp, Yps, iYp, dYp )

    ReusePlanes = .FALSE.
    IF( PRESENT(PlaneCache_Option) )THEN
      IF( PlaneCache_Option%Valid )THEN
        ReusePlanes = PlaneCache_Option%iD == iD &
          .AND. PlaneCache_Option%iT == iT .AND. PlaneCache_Option%iYp == iYp &
          .AND. PlaneCache_Option%D == D .AND. PlaneCache_Option%Ye == Ye &
          .AND. PlaneCache_Option%Ym == Ym
        IF( ReusePlanes ) ReusePlanes = &
          ALL(PlaneCache_Option%DensityNodes == Ds(iD:iD+1)) &
          .AND. ALL(PlaneCache_Option%TemperatureNodes == Ts(iT:iT+1)) &
          .AND. ALL(PlaneCache_Option%YpNodes == Yps(iYp:iYp+1))
        IF( ReusePlanes ) LeptonEPS = PlaneCache_Option%LeptonEPS
      END IF
    END IF

    IF( .NOT. ReusePlanes )THEN
      YeNodes = Yps(iYp:iYp+1)*(Ye/Yp)
      YmNodes = Yps(iYp:iYp+1)*(Ym/Yp)
      DO iPlane = 1, 2
        CALL EvaluatePhotonLeptonEOSPlane &
               ( Ds(iD:iD+1), Ts(iT+iPlane-1), YeNodes, YmNodes, &
                 HelmTableElectrons, HelmTableMuons, LeptonEPS(:,iPlane,:,:) )
      END DO
    END IF

    ! E and S are stored as logarithms of shifted baryon quantities;
    ! pressure is stored linearly. Preserve those representations.
    LogEPS(:,:,:,1) = LOG10(10.0_dp**Es(iD:iD+1,iT:iT+1,iYp:iYp+1)+LeptonEPS(:,:,:,1))
    LogEPS(:,:,:,2) = LOG10(Ps(iD:iD+1,iT:iT+1,iYp:iYp+1)+LeptonEPS(:,:,:,2))
    LogEPS(:,:,:,3) = LOG10(10.0_dp**Ss(iD:iD+1,iT:iT+1,iYp:iYp+1)+LeptonEPS(:,:,:,3))

    NeedDerivatives = PRESENT(Gamma_Option) .OR. PRESENT(Cs_Option) &
      .OR. PRESENT(dPdD_Option) .OR. PRESENT(dPdT_Option) &
      .OR. PRESENT(dEdD_Option) .OR. PRESENT(dEdT_Option)

    IF( NeedDerivatives )THEN
      aD = 1.0_dp/(D*LOG10(Ds(iD+1)/Ds(iD)))
      aT = 1.0_dp/(T*LOG10(Ts(iT+1)/Ts(iT)))
      aYp = LOG(10.0_dp)/(Yps(iYp+1)-Yps(iYp))
      CALL LinearInterpDeriv_Array_Point &
             ( 1, 1, 1, dD, dT, dYp, aD, aT, aYp, OS_P, LogEPS(:,:,:,2), P, dPdD, dPdT )
      CALL LinearInterpDeriv_Array_Point &
             ( 1, 1, 1, dD, dT, dYp, aD, aT, aYp, OS_E, LogEPS(:,:,:,1), E, dEdD, dEdT )
      IF( PRESENT(dPdD_Option) ) dPdD_Option = dPdD
      IF( PRESENT(dPdT_Option) ) dPdT_Option = dPdT
      IF( PRESENT(dEdD_Option) ) dEdD_Option = dEdD
      IF( PRESENT(dEdT_Option) ) dEdT_Option = dEdT
      IF( PRESENT(Gamma_Option) .OR. PRESENT(Cs_Option) )THEN
        Gamma = (D*dPdD+T*dPdT**2.0_dp/(D*dEdT))/P
        IF( PRESENT(Gamma_Option) ) Gamma_Option = Gamma
        IF( PRESENT(Cs_Option) ) Cs_Option = SQRT(Gamma*P/D)
      END IF
    ELSE
      CALL LinearInterp3D_3DArray_Point &
             ( 1, 1, 1, dD, dT, dYp, OS_E, LogEPS(:,:,:,1), E )
      CALL LinearInterp3D_3DArray_Point &
             ( 1, 1, 1, dD, dT, dYp, OS_P, LogEPS(:,:,:,2), P )
    END IF
    CALL LinearInterp3D_3DArray_Point &
           ( 1, 1, 1, dD, dT, dYp, OS_S, LogEPS(:,:,:,3), S )

  END SUBROUTINE ComputeTotalEPS_CS

END MODULE wlEOSComponentsThermodynamicsModule

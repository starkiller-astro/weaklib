MODULE wlBolligSoundSpeed
	
	USE wlKindModule, ONLY: dp
	USE wlEosConstantsModule, ONLY: kmev, rmu, kmev_inv, ergmev, me, cvel, cm3fm3
	USE wlLeptonEOSModule, ONLY: &
		HelmholtzEOSType, MuonEOSType
	USE wlElectronEOS, ONLY: &
		FullHelmEOS, MinimalHelmEOS_rt, ElectronStateType
	USE wlMuonEOS, ONLY: &
		FullMuonEOS, MuonStateType
	USE wlInterpolationUtilitiesModule, ONLY: &
	GetIndexAndDelta_Lin, GetIndexAndDelta_Log, Trilinear, BiLinear, Linear

	IMPLICIT NONE
	PRIVATE
	
	PUBLIC :: CalculateBolligSoundSpeed
	
	REAL(dp), PARAMETER :: ln10 = LOG(10.d0)
	
CONTAINS
	
	! This one also calculates derivatives, but maybe you can provide derivatives ?
	SUBROUTINE CalculateBolligSoundSpeed( D, T, Yp, Ye, Ymu, D_T, T_T, Yp_T, P_T, OS_P, V_T, OS_V, &
		DependentVariable, HelmholtzTable, MuonTable, Gamma, cs2)

    REAL(DP), INTENT(IN)     :: D     , T     , Yp, Ye, Ymu
    REAL(DP), INTENT(IN)     :: D_T(1:), T_T(1:), Yp_T(1:)
    REAL(DP), INTENT(IN)     :: P_T(1:,1:,1:), V_T(1:,1:,1:)
	
    REAL(DP), INTENT(IN)     :: OS_P, OS_V
    CHARACTER(*), INTENT(IN) :: DependentVariable

	TYPE(HelmholtzEOSType), INTENT(IN) :: HelmholtzTable
	TYPE(MuonEOSType), INTENT(IN) :: MuonTable
	
	REAL(DP), INTENT(OUT)    :: Gamma, cs2

	REAL(DP) :: Pbary, Vbary
	REAL(DP) :: dPbarydD, dPbarydT
	REAL(DP) :: dVbarydD, dVbarydT
	REAL(DP) :: Pele, Vele, Eele, Sele
	REAL(DP) :: dPeledD, dPeledT
	REAL(DP) :: dVeledD, dVeledT
	REAL(DP) :: Pmu, Vmu
	REAL(DP) :: dPmudD, dPmudT
	REAL(DP) :: dVmudD, dVmudT
	REAL(DP) :: dD, dT, dYp
	REAL(DP) :: aD, aT, aYp

	INTEGER :: iD, iT, iYp

	TYPE(ElectronStateType) :: ElectronState
	TYPE(MuonStateType) :: MuonState
	
	LOGICAL :: UseEnergy
	
	IF (TRIM(DependentVariable) == 'Energy') THEN
		UseEnergy = .TRUE.
	ELSE IF (TRIM(DependentVariable) == 'Entropy') THEN
		UseEnergy = .FALSE.
	ELSE
		STOP 'Only Energy or Entropy expression are supported for Sound Speed calculation'
	ENDIF 
	
	! ELECTRON PART IS EASY -----------------------------------------------------
	! Initialize temperature, density, yp, Zbar and Abar
	ElectronState % t = T
	ElectronState % rho = D
	ElectronState % abar = 1.0d0 ! these are only used for ion contribution
	ElectronState % zbar = 1.0d0 ! these are only used for ion contribution
	ElectronState % y_e = Ye

	! calculate electron quantities
	CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

	Eele = ElectronState % e + me / rmu * ergmev * ElectronState % y_e ! add back mass to internal energy!
	Sele = ElectronState % s 
	Pele = ElectronState % p

	dPeledD = ElectronState % dpdr
	dPeledT = ElectronState % dpdt
	
	IF (UseEnergy) THEN

		Vele = ElectronState % e
		dVeledD = ElectronState % dedr
		dVeledT = ElectronState % dedt
	
	ELSE
	
		Vele = ElectronState % s
		dVeledD = ElectronState % dsdr
		dVeledT = ElectronState % dsdt
	
	ENDIF

	! BARYONIC PART -----------------------------------------------------
	CALL GetIndexAndDelta_Log( D, D_T, iD, dD )
	CALL GetIndexAndDelta_Log( T, T_T, iT, dT )
	CALL GetIndexAndDelta_Lin( Yp, Yp_T, iYp, dYp )
	
	aD = 1.0_dp / ( D * LOG10( D_T(iD+1) / D_T(iD) ) )
	aT = 1.0_dp / ( T * LOG10( T_T(iT+1) / T_T(iT) ) )
	aYp = ln10 / ( Yp_T(iYp+1) - Yp_T(iYp) )
  
	! is this linear or log derivative? I think it's log ??????
	CALL LinearInterpDeriv3D_2DArray_Point &
		   ( iD, iT, iYp, dD, dT, dYp, aD, aT, aYp, OS_P, P_T, Pbary, &
			 dPbarydD, dPbarydT )
	
	! is this linear or log derivative?
	CALL LinearInterpDeriv3D_2DArray_Point &
		   ( iD, iT, iYp, dD, dT, dYp, aD, aT, aYp, OS_V, V_T, Vbary, &
			 dVbarydD, dVbarydT )

	! MUON PART -----------------------------------------------------
	IF (( D * Ymu .lt. MuonTable % rhoymu(1) ) .or. (T .lt. MuonTable % t(1))) THEN
	
	  Pmu = 0.0d0
	  dPmudD = 0.0d0
	  dPmudT = 0.0d0

	  Vmu = 0.0d0
	  dVmudD = 0.0d0
	  dVmudT = 0.0d0
	
	ELSE
	
	  CALL GetIndexAndDelta_Log( D * Ymu, MuonTable % rhoymu(:), iD, dD )
	  CALL GetIndexAndDelta_Log( T, MuonTable % t(:), iT, dT )
	  
	  aD = 1.0_dp / ( D * LOG10( MuonTable % rhoymu(iD+1) / MuonTable % rhoymu(iD) ) )
	  aT = 1.0_dp / ( T * LOG10( MuonTable % t(iT+1) / MuonTable % t(iT) ) )

	  ! is this linear or log derivative?
	  CALL LinearInterpDeriv2D_2DArray_Point &
			 ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % p), Pmu, &
			   dPmudD, dPmudT )
	  
	  dPmudD = dPmudD * Ymu ! make sure the derivative is wr2 rho, not rhoymu

	  IF (UseEnergy) THEN
	  
		  ! make sure you handle rest mass correctly
		  CALL LinearInterpDeriv2D_2DArray_Point &
				 ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % e), Vmu, &
				   dVmudD, dVmudT )
				   
		  dVmudD = dVmudD * Ymu ! make sure the derivative is wr2 rho, not rhoymu

	  ELSE
	  
		  ! make sure you handle rest mass correctly
		  CALL LinearInterpDeriv2D_2DArray_Point &
				 ( iD, iT, dD, dT, aD, aT, 0.0_dp, MuonTable % s, Vmu, &
				   dVmudD, dVmudT )
				   
		  dVmudD = dVmudD * Ymu ! make sure the derivative is wr2 rho, not rhoymu

	  ENDIF
	ENDIF
	
	! Check that yu are doing this correctly, really check this like
	! do not trust me at all. D and T in front take care of the derivative 
	! wr2 logrho and logT. The rho multiplying the denominator in the second 
	! one makes sure that you have erg/cm^3 instead of erg/g
	Gamma = (D*(dPbarydD + dPeledD + dPmudD) + &
				T*(dPbarydT + dPeledT + dPmudT)**2.0_DP / &
				(D*(dVbarydT + dVeledT + dVmudT)) ) / &
				(Pbary + Pele + Pmu)
				
	! relativistic definition
	cs2 = Gamma * (Pbary + Pele + Pmu) / &
		(Pbary + Pele + Pmu + &
		(Vbary + Vele + Vmu + cvel**2.0_DP )*D ) * cvel**2

	END SUBROUTINE CalculateBolligSoundSpeed

  ! This one interpolates and calculates derivatives along first and second dimension
  SUBROUTINE LinearInterpDeriv3D_2DArray_Point &
	( iY1, iY2, iY3, dY1, dY2, dY3, aY1, aY2, aY3, OS, Table, &
      Interpolant, dIdY1, dIdY2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: iY1, iY2, iY3
    REAL(dp), INTENT(in)  :: dY1, dY2, dY3, aY1, aY2, aY3, OS, Table(1:,1:,1:)
    REAL(dp), INTENT(out) :: Interpolant, dIdY1, dIdY2

    REAL(dp) :: p000, p100, p010, p110, p001, p101, p011, p111

    p000 = Table(iY1  , iY2  , iY3   )
    p100 = Table(iY1+1, iY2  , iY3   )
    p010 = Table(iY1  , iY2+1, iY3   )
    p110 = Table(iY1+1, iY2+1, iY3   )
    p001 = Table(iY1  , iY2  , iY3+1 )
    p101 = Table(iY1+1, iY2  , iY3+1 )
    p011 = Table(iY1  , iY2+1, iY3+1 )
    p111 = Table(iY1+1, iY2+1, iY3+1 )

    Interpolant &
      = 10.0d0**( &
          TriLinear &
            ( p000, p100, p010, p110, &
              p001, p101, p011, p111, &
              dY1, dY2, dY3 ) ) - OS

    dIdY1 &
      = (Interpolant + OS) * aY1 &
          * dTriLineardX1 &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, &
                dY2, dY3 )

    dIdY2 &
      = (Interpolant + OS) * aY2 &
          * dTriLineardX2 &
              ( p000, p100, p010, p110, &
                p001, p101, p011, p111, &
                dY1, dY3 )

  END SUBROUTINE LinearInterpDeriv3D_2DArray_Point

  SUBROUTINE LinearInterpDeriv2D_2DArray_Point &
	  ( iY1, iY2, dY1, dY2, aY1, aY2, OS, Table, &
      Interpolant, dIdY1, dIdY2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in) :: iY1, iY2
    REAL(dp), INTENT(in) :: dY1, dY2, aY1, aY2, OS, Table(1:,1:)
    REAL(dp) :: Interpolant, dIdY1, dIdY2

    REAL(dp) :: p00, p10, p01, p11

    p00 = Table(iY1  , iY2  )
    p10 = Table(iY1+1, iY2  )
    p01 = Table(iY1  , iY2+1)
    p11 = Table(iY1+1, iY2+1)

    Interpolant &
      = 10.0d0**( &
          BiLinear &
            ( p00, p10, p01, p11, &
              dY1, dY2) ) - OS

    dIdY1 &
      = (Interpolant + OS) * aY1 &
          * dBiLineardX1( p00, p10, p01, p11, dY2 )

    dIdY2 &
      = (Interpolant + OS) * aY2 &
          * dBiLineardX2( p00, p10, p01, p11, dY1 )

  END SUBROUTINE LinearInterpDeriv2D_2DArray_Point

! Copied over from InterpolationModule
  REAL(dp) FUNCTION dBiLineardX1 &
    ( p00, p10, p01, p11, dX2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p00, p10, p01, p11, dX2

    dBiLineardX1 &
      = Linear(p10, p11, dX2 ) &
      - Linear(p00, p01, dX2 )

    RETURN
  END FUNCTION dBiLineardX1


  REAL(dp) FUNCTION dBiLineardX2 &
    ( p00, p10, p01, p11, dX1 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p00, p10, p01, p11, dX1

    dBiLineardX2 &
      = Linear(p01, p11, dX1 ) &
      - Linear(p00, p10, dX1 )

    RETURN
  END FUNCTION dBiLineardX2


  REAL(dp) FUNCTION dTriLineardX1 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX2, dX3 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, p001, p101, p011, p111, dX2, dX3

    dTrilineardX1 &
      = Bilinear(p100, p110, p101, p111, dX2, dX3) &
      - Bilinear(p000, p010, p001, p011, dX2, dX3)

    RETURN
  END FUNCTION dTriLineardX1


  REAL(dp) FUNCTION dTriLineardX2 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX3 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX3

    dTrilineardX2 &
      = Bilinear(p010, p110, p011, p111, dX1, dX3) &
      - Bilinear(p000, p100, p001, p101, dX1, dX3)

    RETURN
  END FUNCTION dTriLineardX2


  REAL(dp) FUNCTION dTriLineardX3 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2

    dTrilineardX3 &
      = Bilinear(p001, p101, p011, p111, dX1, dX2) &
      - Bilinear(p000, p100, p010, p110, dX1, dX2)

    RETURN
  END FUNCTION dTriLineardX3

END MODULE wlBolligSoundSpeed
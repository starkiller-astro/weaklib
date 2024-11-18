MODULE wlSoundSpeedModule
	
	USE wlKindModule, ONLY: dp
	USE wlEosConstantsModule, ONLY: kmev, rmu, kmev_inv, ergmev, me, mmu, cvel, cm3fm3
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
	
	PUBLIC :: CalculatewlSoundSpeed
	
	REAL(dp), PARAMETER :: ln10 = LOG(10.d0)
	
CONTAINS
	
	! This one also calculates derivatives, but maybe you can provide derivatives ?
	SUBROUTINE CalculatewlSoundSpeed( D, T, Ye, Ym, D_T, T_T, Yp_T, P_T, OS_P, E_T, OS_E, &
		HelmholtzTable, MuonTable, Gamma, cs2)

    REAL(DP), INTENT(IN)     :: D     , T     , Ye, Ym
    REAL(DP), INTENT(IN)     :: D_T(1:), T_T(1:), Yp_T(1:)
    REAL(DP), INTENT(IN)     :: P_T(1:,1:,1:), E_T(1:,1:,1:)
	
    REAL(DP), INTENT(IN)     :: OS_P, OS_E

	TYPE(HelmholtzEOSType), INTENT(IN) :: HelmholtzTable
	TYPE(MuonEOSType), INTENT(IN) :: MuonTable
	
	REAL(DP), INTENT(OUT)    :: Gamma, cs2

	REAL(DP) :: Pbary, Ebary
	REAL(DP) :: dPbarydD, dPbarydT
	REAL(DP) :: dEbarydD, dEbarydT
	REAL(DP) :: Pele, Eele, Sele
	REAL(DP) :: dPeledD, dPeledT
	REAL(DP) :: dEeledD, dEeledT
	REAL(DP) :: Pmu, Emu
	REAL(DP) :: dPmudD, dPmudT
	REAL(DP) :: dEmudD, dEmudT
	REAL(DP) :: dD, dT, dYp
	REAL(DP) :: aD, aT, aYp
	REAL(DP) :: h
	
	INTEGER :: iD, iT, iYp

	TYPE(ElectronStateType) :: ElectronState
	TYPE(MuonStateType) :: MuonState
		
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


	! BARYONIC PART -----------------------------------------------------
	CALL GetIndexAndDelta_Log( D, D_T, iD, dD )
	CALL GetIndexAndDelta_Log( T, T_T, iT, dT )
	CALL GetIndexAndDelta_Lin( Ye + Ym, Yp_T, iYp, dYp )
	
	aD = 1.0_dp / ( D * LOG10( D_T(iD+1) / D_T(iD) ) )
	aT = 1.0_dp / ( T * LOG10( T_T(iT+1) / T_T(iT) ) )
	aYp = ln10 / ( Yp_T(iYp+1) - Yp_T(iYp) )
  
	! is this linear or log derivative? I think it's log ??????
	CALL LinearInterpDeriv3D_2DArray_Point &
		   ( iD, iT, iYp, dD, dT, dYp, aD, aT, aYp, OS_P, P_T, Pbary, &
			 dPbarydD, dPbarydT )
	
	! is this linear or log derivative?
	CALL LinearInterpDeriv3D_2DArray_Point &
		   ( iD, iT, iYp, dD, dT, dYp, aD, aT, aYp, OS_E, E_T, Ebary, &
			 dEbarydD, dEbarydT )

	! MUON PART -----------------------------------------------------
	IF (( D * Ym .lt. MuonTable % rhoym(1) ) .or. (T .lt. MuonTable % t(1))) THEN
	
	  Pmu = 0.0d0
	  dPmudD = 0.0d0
	  dPmudT = 0.0d0

	  Emu = 0.0d0
	  dEmudD = 0.0d0
	  dEmudT = 0.0d0
	
	ELSE
	
	  CALL GetIndexAndDelta_Log( D * Ym, MuonTable % rhoym(:), iD, dD )
	  CALL GetIndexAndDelta_Log( T, MuonTable % t(:), iT, dT )
	  
	  aD = 1.0_dp / ( D * LOG10( MuonTable % rhoym(iD+1) / MuonTable % rhoym(iD) ) )
	  aT = 1.0_dp / ( T * LOG10( MuonTable % t(iT+1) / MuonTable % t(iT) ) )

	  ! is this linear or log derivative?
	  CALL LinearInterpDeriv2D_2DArray_Point &
			 ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % p), Pmu, &
			   dPmudD, dPmudT )
	  
	  dPmudD = dPmudD * Ym ! make sure the derivative is wr2 rho, not rhoym
	  
	  CALL LinearInterpDeriv2D_2DArray_Point &
			 ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % e), Emu, &
			   dEmudD, dEmudT )
			   
	  dEmudD = dEmudD * Ym ! make sure the derivative is wr2 rho, not rhoym
	  Emu = Emu + mmu / rmu * ergmev * Ym ! make sure you handle rest mass correctly

	ENDIF
	
	! Check that yu are doing this correctly, really check this like
	! do not trust me at all. D and T in front take care of the derivative 
	! wr2 logrho and logT. The rho multiplying the denominator in the second 
	! one makes sure that you have erg/cm^3 instead of erg/g
	Gamma = (D*(dPbarydD + dPeledD + dPmudD) + &
				T*(dPbarydT + dPeledT + dPmudT)**2.0_DP / &
				(D*(dEbarydT + dEeledT + dEmudT)) ) / &
				(Pbary + Pele + Pmu)
				
	! relativistic definition with enthalpy
	h = (1.0_dp + (Ebary + Eele + Emu)/cvel**2 + (Pbary + Pele + Pmu)/D/cvel**2)
	
	cs2 = Gamma * (Pbary + Pele + Pmu) / (D*h)

	END SUBROUTINE CalculatewlSoundSpeed

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

! Bonus subroutines, currently unused
	! This is taken from the stellarcollapse.org routines (Christian Ott's code)
	subroutine derivatives_production(igamma, nrho, ntemp, nye, logrho, logtemp, ye, &
		eos_table, energy_shift, cs2, gamma)
		
		integer, intent(in) :: igamma
		integer, intent(in) :: nrho, ntemp, nye
		real(dp), intent(in) :: eos_table(nrho,ntemp,nye,3)
		real(dp), intent(in) :: energy_shift
		real(dp), intent(in) :: logrho(nrho), logtemp(ntemp), ye(nye)
		real(dp), allocatable, intent(out) :: cs2(:,:,:)
		real(dp), allocatable, intent(out) :: gamma(:,:,:)
		
		integer i,j,k
		real(dp) :: dx,x1,x2,f1,f2,z,zz,h,cp,cv,beta_v,gamma_ad,kappa_t
		
		real(dp), allocatable :: dedT(:,:,:)
		real(dp), allocatable :: dsdlnT(:,:,:)
		real(dp), allocatable :: dsdlnrho(:,:,:)
		real(dp), allocatable :: dpdrhoT(:,:,:)
		real(dp), allocatable :: dedrhoT(:,:,:)
		real(dp), allocatable :: dpdT(:,:,:)
		real(dp), allocatable :: dpdrho(:,:,:)
		real(dp), allocatable :: dpde(:,:,:)
		integer :: ipress = 1
		integer :: ientropy = 2
		integer :: ienergy = 3
		
		! igamma: There are multiple ways to compute gamma1, via the entropy
		!         turns out to be pretty good.
			! 1 -> Gamma1 = &
			!      dlnP/dlnrho|T,Y_e - ds/dlnrho|T,Y_e * (dlnP/dlnT / ds/dlnT)_rho,Y_e
			! 2 -> Gamma1 = P/rho * ( dpdrho|e,Y_e + P/rho**2 dpde|rho,Y_e
		! 3 -> Gamma1 = dlnP/dlnrho|T,Y_e + T * (dpdT|rho,Y_e)**2 / &
		!      (P*rho * dedT|rho,Y_e)
		
		!###########################################################################  
		! dedT, e is erg in log10, T is in MeV in log 10
		allocate(dedT(nrho,ntemp,nye))
		dedT = 0.0d0
		
		do k=1,nye
			do i=1,nrho
				do j=2,ntemp-1
					x1 = logtemp(j-1)
					f1 = eos_table(i,j-1,k,ienergy)
					x2 = logtemp(j+1)
					f2 = eos_table(i,j+1,k,ienergy)
					dedT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)  &
					* 10.0d0**eos_table(i,j,k,ienergy)
				enddo
				
				! boundaries: one-sided derivative
				x2 = logtemp(2)
				x1 = logtemp(1)
				f2 = eos_table(i,2,k,ienergy)
				f1 = eos_table(i,1,k,ienergy)
				dedT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1)  &
                * 10.0d0**eos_table(i,1,k,ienergy)
				x2 = logtemp(ntemp)
				x1 = logtemp(ntemp-1)
				f2 = eos_table(i,ntemp,k,ienergy)
				f1 = eos_table(i,ntemp-1,k,ienergy)
				dedT(i,ntemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(ntemp) &
				* 10.0d0**eos_table(i,ntemp,k,ienergy) 
			enddo
		enddo
		
		!########################################################################### 
		
		!dpdT, p is in dyn/cm^2 in log10, T is in MeV in log 10
		!dsdlnT, s is in k_B/baryon not in log 10, T is in meV in log 10
		allocate(dpdT(nrho,ntemp,nye))
		allocate(dsdlnT(nrho,ntemp,nye))
		dpdT = 0.0d0
		dsdlnT = 0.0d0
		
		do k=1,nye
			do i=1,nrho
				do j=2,ntemp-1
					x1 = logtemp(j-1)
					f1 = eos_table(i,j-1,k,ipress)
					x2 = logtemp(j+1)
					f2 = eos_table(i,j+1,k,ipress)
					dpdT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)  &
					* 10.0d0**eos_table(i,j,k,ipress)
					
					x1 = logtemp(j-1)
					f1 = log10(eos_table(i,j-1,k,ientropy))
					x2 = logtemp(j+1)
					f2 = log10(eos_table(i,j+1,k,ientropy))
					dsdlnT(i,j,k) = (f2-f1)/(x2-x1) * eos_table(i,j,k,ientropy)
					
				enddo
				
				! boundaries: one-sided derivative
				x1 = logtemp(1)
				f1 = eos_table(i,1,k,ipress)
				x2 = logtemp(2)
				f2 = eos_table(i,2,k,ipress)
				dpdT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
				10.0d0**eos_table(i,1,k,ipress)
				
				x1 = logtemp(1)
				f1 = log10(eos_table(i,1,k,ientropy))
				x2 = logtemp(2)
				f2 = log10(eos_table(i,2,k,ientropy))
				dsdlnT(i,1,k) = (f2-f1)/(x2-x1) * eos_table(i,1,k,ientropy)
				
				
				x1 = logtemp(ntemp-1)
				f1 = eos_table(i,ntemp-1,k,ipress)
				x2 = logtemp(ntemp)
				f2 = eos_table(i,ntemp,k,ipress)
				dpdT(i,ntemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
				10.0d0**eos_table(i,1,k,ipress)
				
				x1 = logtemp(ntemp-1)
				f1 = log10(eos_table(i,ntemp-1,k,ientropy))
				x2 = logtemp(ntemp)
				f2 = log10(eos_table(i,ntemp,k,ientropy))
				dsdlnT(i,ntemp,k) = (f2-f1)/(x2-x1) & 
				* eos_table(i,ntemp,k,ientropy)
				
			enddo
		enddo
		
		!########################################################################### 
		! dp/drho|T, p is in dyn/cm^2 in log 10, T is in MeV in log 10
			! ds/dlnrho|T, s is in k_b/baryon not in log 10, T is in MeV in log 10
			! de/drho|T, e is in erg in log 10, T is in MeV in log 10
			! dp/drho|e = dp/drho|T + dp/dT * (-de/drho|T) / de/dT
		! dp/de|rho = dp/dT / de/dT
		
		allocate(dpdrho(nrho,ntemp,nye))
		allocate(dpde(nrho,ntemp,nye))
		allocate(dsdlnrho(nrho,ntemp,nye))
		allocate(dpdrhoT(nrho,ntemp,nye))
		allocate(dedrhoT(nrho,ntemp,nye))
		dpdrho   = 0.0d0
		dpde     = 0.0d0
		dsdlnrho = 0.0d0
		dpdrhoT  = 0.0d0
		dedrhoT  = 0.0d0
		
		do k=1,nye
			do j=1,ntemp
				do i=2,nrho-1
					x1 = logrho(i-1)
					x2 = logrho(i+1)
					f1 = eos_table(i-1,j,k,ipress)
					f2 = eos_table(i+1,j,k,ipress)
					dpdrhoT(i,j,k) = (f2-f1)/(x2-x1) /  10.0d0**logrho(i) &
					* 10.0d0**eos_table(i,j,k,ipress)
					
					x1 = logrho(i-1)
					x2 = logrho(i+1)
					f1 = eos_table(i-1,j,k,ienergy)
					f2 = eos_table(i+1,j,k,ienergy)
					dedrhoT(i,j,k) = (f2-f1)/(x2-x1) /  10.0d0**logrho(i) &
					* 10.0d0**eos_table(i,j,k,ienergy)
					
					x1 = logrho(i-1)
					x2 = logrho(i+1)
					f1 = log10(eos_table(i-1,j,k,ientropy))
					f2 = log10(eos_table(i+1,j,k,ientropy))
					dsdlnrho(i,j,k) = (f2-f1)/(x2-x1) * eos_table(i,j,k,ientropy)
				enddo
				
				! boundaries: one-sided derivative
				x1 = logrho(1)
				x2 = logrho(2)
				f1 = eos_table(1,j,k,ipress)
				f2 = eos_table(2,j,k,ipress)
				dpdrhoT(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
				* 10.0d0**eos_table(1,j,k,ipress)
				
				f1 = eos_table(1,j,k,ienergy)
				f2 = eos_table(2,j,k,ienergy)
				dedrhoT(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
				* 10.0d0**eos_table(1,j,k,ienergy)
				
				x1 = logrho(nrho-1)
				x2 = logrho(nrho)
				f1 = eos_table(nrho-1,j,k,ipress)
				f2 = eos_table(nrho,j,k,ipress)
				dpdrhoT(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho) &
				* 10.0d0**eos_table(nrho,j,k,ipress)
				
				f1 = eos_table(nrho-1,j,k,ienergy)
				f2 = eos_table(nrho,j,k,ienergy)
				dedrhoT(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho) &
				* 10.0d0**eos_table(nrho,j,k,ienergy)
				
				x1 = logrho(1)
				x2 = logrho(2)
				f1 = log10(eos_table(1,j,k,ientropy))
				f2 = log10(eos_table(2,j,k,ientropy))
				dsdlnrho(1,j,k) = (f2-f1)/(x2-x1) * eos_table(1,j,k,ientropy)
				
				x1 = logrho(nrho-1)
				x2 = logrho(nrho)
				f1 = log10(eos_table(nrho-1,j,k,ientropy))
				f2 = log10(eos_table(nrho,j,k,ientropy))
				dsdlnrho(nrho,j,k) = (f2-f1)/(x2-x1) * eos_table(nrho,j,k,ientropy)
				
			enddo
		enddo
		
		do k=1,nye
			do j=1,ntemp
				do i=1,nrho
					dpdrho(i,j,k) = dpdrhoT(i,j,k) + dpdT(i,j,k) * &
					(-dedrhoT(i,j,k))/dedT(i,j,k)
					
					dpde(i,j,k) = dpdT(i,j,k) / dedT(i,j,k)
				enddo
			enddo
		enddo
		
		!########################################################################### 
		allocate(gamma(nrho,ntemp,nye))
		allocate(cs2(nrho,ntemp,nye))
		cs2 = 0.0d0
		gamma = 0.0d0
		
		do k=1,nye
			do j=1,ntemp
				do i=1,nrho
					z = 10.0d0**logrho(i)/10.0d0**eos_table(i,j,k,ipress)
					zz = 10.0d0**logtemp(j)/10.0d0**eos_table(i,j,k,ipress)
					
					if(igamma .eq. 1) then
						gamma(i,j,k) = z * dpdrhoT(i,j,k) & 
						- dsdlnrho(i,j,k) * zz * dpdT(i,j,k) / dsdlnT(i,j,k)
					endif
					
					if(igamma .eq. 2) then
						gamma(i,j,k) =  dpdrho(i,j,k) + 10.0d0**eos_table(i,j,k,ipress) / &
						(10.0d0**logrho(i))**2 * dpde(i,j,k)
						gamma(i,j,k) = z * gamma(i,j,k)
					endif
					
					if(igamma .eq. 3) then
						gamma(i,j,k) = 10.0d0**logrho(i)/10.0d0**eos_table(i,j,k,ipress) &
						* dpdrhoT(i,j,k) &
						+ 10.0d0**logtemp(j)*dpdT(i,j,k)**2 / &
						(10.0d0**eos_table(i,j,k,ipress)*10.0d0**logrho(i) * dedT(i,j,k))
					endif
					
					cs2(i,j,k) = gamma(i,j,k) / z
					cs2(i,j,k) = gamma(i,j,k) * 10.0d0**eos_table(i,j,k,ipress) / &
					(10.0d0**eos_table(i,j,k,ipress) + (10.0d0**eos_table(i,j,k,ienergy) - energy_shift - &
					0.511d0 / rmu * ergmev * ye(k) + cvel**2) * 10.0d0**logrho(i) ) * cvel**2
					
					if(igamma .eq. 5) then
						h = 10.0d0**eos_table(i,j,k,ipress) + (10.0d0**eos_table(i,j,k,ienergy) - energy_shift - &
							0.511d0 / rmu * ergmev * ye(k) + cvel**2) * 10.0d0**logrho(i)
						
						kappa_t = 10.0d0**logrho(i) * dpdrho(i,j,k) * (cm3fm3/ergmev)
						kappa_t = 1.0d0 / kappa_t
						beta_v = dpdT(i,j,k) * (cm3fm3/ergmev) / kmev
						cv = 1.0d0 / (10.0d0**logrho(i)/rmu*cm3fm3) * dsdlnT(i,j,k) / LOG(10.0d0) 
						cp = cv + 10.0d0**logtemp(j)*kmev/(10.0d0**logrho(i)/rmu*cm3fm3) * beta_v * kappa_t * beta_v
						gamma_ad = cp/cv
													
						cs2(i,j,k) = gamma_ad / (h * (cm3fm3/ergmev) * kappa_t )
						gamma(i,j,k) = gamma_ad / 10.0d0**eos_table(i,j,k,ipress)*ergmev/cm3fm3 / kappa_t
						cs2 = cs2 * cvel**2
						
					endif
				enddo
			enddo
		enddo
		
		DEALLOCATE(dedT)
		DEALLOCATE(dsdlnT)
		DEALLOCATE(dsdlnrho)
		DEALLOCATE(dpdrhoT)
		DEALLOCATE(dedrhoT)
		DEALLOCATE(dpdT)
		DEALLOCATE(dpdrho)
		DEALLOCATE(dpde)
		
	end subroutine derivatives_production

END MODULE wlSoundSpeedModule
MODULE wlElectronEOS
	
	USE wlKindModule, ONLY: dp
	USE wlExtNumericalModule, ONLY: zero, one, pi, half
	USE wlLeptonEOSModule, ONLY: &
		HelmholtzEOSType
	USE wlExtPhysicalConstantsModule, ONLY: &
		ergmev, kmev, kmev_inv, rmu
		
	IMPLICIT NONE
	PRIVATE
	
	! set some parameters (potentially this could be changed? But why would you want to cahnge them?)
	LOGICAL,  PARAMETER :: input_is_constant = .true.
    REAL(dp), PARAMETER :: ttol = 1.0d-8
    REAL(dp), PARAMETER :: dtol = 1.0d-8

	INTEGER, PARAMETER :: max_newton = 100
    INTEGER, PARAMETER :: eos_input_rt = 1  ! rho, T are inputs
	INTEGER, PARAMETER :: eos_input_rh = 2  ! rho, h are inputs
	INTEGER, PARAMETER :: eos_input_tp = 3  ! T, p are inputs
	INTEGER, PARAMETER :: eos_input_rp = 4  ! rho, p are inputs
	INTEGER, PARAMETER :: eos_input_re = 5  ! rho, e are inputs
	INTEGER, PARAMETER :: eos_input_ps = 6  ! p, s are inputs
	INTEGER, PARAMETER :: eos_input_ph = 7  ! p, h are inputs
	INTEGER, PARAMETER :: eos_input_th = 8  ! T, h are inputs
	
	! these are used to allow for a generic interface to the
	! root finding
	INTEGER, PARAMETER :: itemp = 1
	INTEGER, PARAMETER :: idens = 2
	INTEGER, PARAMETER :: iener = 3
	INTEGER, PARAMETER :: ienth = 4
	INTEGER, PARAMETER :: ientr = 5
	INTEGER, PARAMETER :: ipres = 6

    ! Physical constants
    REAL(dp), PARAMETER :: h       = 6.62607535359d-27
    REAL(dp), PARAMETER :: hbar    = 0.5d0 * h/pi
    REAL(dp), PARAMETER :: qe      = 4.8032042712d-10
    REAL(dp), PARAMETER :: avo_eos = 6.022141d+23
    REAL(dp), PARAMETER :: clight  = 2.99792d+10
    REAL(dp), PARAMETER :: kerg    = 1.380648792741d-16
    REAL(dp), PARAMETER :: ev2erg_eos  = 1.602177d-12
    REAL(dp), PARAMETER :: kev     = kerg/ev2erg_eos
    REAL(dp), PARAMETER :: amu     = 1.0d0/avo_eos

    REAL(dp), PARAMETER :: ssol    = 5.67051d-5
    REAL(dp), PARAMETER :: asol    = 4.0d0 * ssol / clight
	
    ! Some other useful combinations of the constants
    REAL(dp), PARAMETER :: sioncon = (2.0d0 * pi * amu * kerg)/(h*h)
    REAL(dp), PARAMETER :: forth   = 4.0d0/3.0d0
    REAL(dp), PARAMETER :: forpi   = 4.0d0 * pi
    REAL(dp), PARAMETER :: forthpi = forth * pi
    REAL(dp), PARAMETER :: kergavo = kerg * avo_eos
    REAL(dp), PARAMETER :: ikavo   = 1.0d0/kergavo
    REAL(dp), PARAMETER :: asoli3  = asol/3.0d0
    REAL(dp), PARAMETER :: light2  = clight * clight

    ! Constants used for the Coulomb corrections
    REAL(dp), PARAMETER :: a1    = -0.898004d0
    REAL(dp), PARAMETER :: b1    =  0.96786d0
    REAL(dp), PARAMETER :: c1    =  0.220703d0
    REAL(dp), PARAMETER :: d1    = -0.86097d0
    REAL(dp), PARAMETER :: e1    =  2.5269d0
    REAL(dp), PARAMETER :: a2    =  0.29561d0
    REAL(dp), PARAMETER :: b2    =  1.9885d0
    REAL(dp), PARAMETER :: c2    =  0.288675d0
    REAL(dp), PARAMETER :: onethird = 1.0d0/3.0d0
    REAL(dp), PARAMETER :: esqu = qe * qe
	
	TYPE, PUBLIC :: ElectronStateType
		
		REAL(dp) :: rho
		REAL(dp) :: T
		REAL(dp) :: p
		REAL(dp) :: e
		REAL(dp) :: h
		REAL(dp) :: s
		
		REAL(dp) :: dpdT
		REAL(dp) :: dpdr
		REAL(dp) :: dedT
		REAL(dp) :: dedr
		REAL(dp) :: dhdT
		REAL(dp) :: dhdr
		REAL(dp) :: dsdT
		REAL(dp) :: dsdr
		REAL(dp) :: dpde
		REAL(dp) :: dpdr_e
		
		REAL(dp) :: cv
		REAL(dp) :: cp
		REAL(dp) :: xne
		REAL(dp) :: xnp
		REAL(dp) :: eta
		REAL(dp) :: detadt
		REAL(dp) :: pele
		REAL(dp) :: ppos
		REAL(dp) :: mu
		REAL(dp) :: mu_e
		REAL(dp) :: y_e
		REAL(dp) :: gam1
		REAL(dp) :: cs
		
		REAL(dp) :: abar
		REAL(dp) :: zbar
		
		REAL(dp) :: conductivity
		
	END TYPE ElectronStateType
	
	TYPE, PUBLIC :: HelmholtzStateType
			
		REAL(dp) :: rho
		REAL(dp) :: T
		REAL(dp) :: ye
		REAL(dp) :: abar
		REAL(dp) :: zbar
		
		REAL(dp) :: prad
		REAL(dp) :: pion
		REAL(dp) :: pele
		REAL(dp) :: pcoul
	
		REAL(dp) :: erad
		REAL(dp) :: eion
		REAL(dp) :: eele
		REAL(dp) :: ecoul
		
		REAL(dp) :: srad
		REAL(dp) :: sion
		REAL(dp) :: sele
		REAL(dp) :: scoul

		REAL(dp) :: dpraddd
		REAL(dp) :: dpiondd
		REAL(dp) :: dpepdd
		REAL(dp) :: dpcouldd

		REAL(dp) :: deraddd
		REAL(dp) :: deiondd
		REAL(dp) :: deepdd
		REAL(dp) :: decouldd
		
		REAL(dp) :: dsraddd
		REAL(dp) :: dsiondd
		REAL(dp) :: dsepdd
		REAL(dp) :: dscouldd
		
	END TYPE HelmholtzStateType
	
	PUBLIC :: FullHelmEOS, HelmEOSComponents
	
CONTAINS

    SUBROUTINE FullHelmEOS(input, HelmTable, state, include_ion_contribution, do_coulomb)
		
        !..input arguments
        INTEGER, INTENT(IN) :: input
		LOGICAL, INTENT(IN) :: include_ion_contribution, do_coulomb
        TYPE(HelmholtzEOSType), INTENT(IN) :: HelmTable
        TYPE (ElectronStateType), INTENT(INOUT) :: state
		
        !..rows to store EOS data
        REAL(dp) :: temp_row, &
		den_row, &
		abar_row, &
		zbar_row, &
		ye_row, &
		etot_row, &
		ptot_row, &
		cv_row, &
		cp_row,  &
		xne_row, &
		xnp_row, &
		etaele_row, &
		detadt_row, &
		pele_row, &
		ppos_row, &
		dpd_row,  &
		dpt_row, &
		dpa_row, &
		dpz_row,  &
		ded_row, &
		det_row, &
		dea_row,  &
		dez_row,  &
		stot_row, &
		dsd_row, &
		dst_row, &
		htot_row, &
		dhd_row, &
		dht_row, &
		dpe_row, &
		dpdr_e_row, &
		gam1_row, &
		cs_row
		
		!..declare local variables
		
		logical :: single_iter, double_iter, converged
		INTEGER :: var, dvar, var1, var2, iter
		REAL(dp) :: v_want
		REAL(dp) :: v1_want, v2_want
		REAL(dp) :: xnew, xtol, dvdx, smallx, error, v
		REAL(dp) :: v1, v2, dv1dt, dv1dr, dv2dt,dv2dr, delr, error1, error2, told, rold, tnew, rnew, v1i, v2i
		
		REAL(dp) :: x,y,z,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
		dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
		dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
		deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
		dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
		sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
		dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
		gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
		detadt,detadd,xnefer,dxnedt,dxnedd,s, &
		temp,den,abar,zbar,ytot1,ye,din
		
		
		!..for the interpolations
		INTEGER  :: iat,jat
		REAL(dp) :: free,df_d,df_t,df_tt,df_dt
		REAL(dp) :: xt,xd,mxt,mxd,fi(36), &
		si0t,si1t,si2t,si0mt,si1mt,si2mt, &
		si0d,si1d,si2d,si0md,si1md,si2md, &
		dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
		dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
		ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt
		
		!..for the coulomb corrections
		REAL(dp) :: dsdd,dsda,lami,inv_lami,lamida,lamidd,     &
		plasg,plasgdd,plasgdt,plasgda,plasgdz,     &
		ecoul,decouldd,decouldt,decoulda,decouldz, &
		pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
		scoul,dscouldd,dscouldt,dscoulda,dscouldz
				
		REAL(dp) :: p_temp, e_temp
		
		! extra variables I added that were previously allocated
		REAL(dp) :: smallt, smalld, hight, highd
		REAL(dp) :: tstp, dstp, dstpi, tstpi
		
		smallt = HelmTable % mintemp
		smalld = HelmTable % mindens
		hight = HelmTable % maxtemp
		highd = HelmTable % maxdens
				
		tstp  = (LOG10(hight) - LOG10(smallt))/float(HelmTable % nPointsTemp -1)
		dstp  = (LOG10(highd) - LOG10(smalld))/float(HelmTable % nPointsDen -1)
		
		tstpi = one / tstp
		dstpi = one / dstp
		
		temp_row = state % T
		den_row  = state % rho
		abar_row = state % abar
		zbar_row = state % zbar
		ye_row   = state % y_e
		
		! Initial setup for iterations
		
		single_iter = .false.
		double_iter = .false.
		
		IF (input .eq. eos_input_rt) THEN
			
			! Nothing to DO here.
		
		ELSEIF (input .eq. eos_input_rh) THEN
			
			single_iter = .true.
			v_want = state % h
			var  = ienth
			dvar = itemp
		
		ELSEIF (input .eq. eos_input_tp) THEN
			
			single_iter = .true.
			v_want = state % p
			var  = ipres
			dvar = idens
		
		ELSEIF (input .eq. eos_input_rp) THEN
			
			single_iter = .true.
			v_want = state % p
			var  = ipres
			dvar = itemp
		
		ELSEIF (input .eq. eos_input_re) THEN
			
			single_iter = .true.
			v_want = state % e
			var  = iener
			dvar = itemp
		
		ELSEIF (input .eq. eos_input_ps) THEN
			
			double_iter = .true.
			v1_want = state % p
			v2_want = state % s
			var1 = ipres
			var2 = ientr
		
		ELSEIF (input .eq. eos_input_ph) THEN
			
			double_iter = .true.
			v1_want = state % p
			v2_want = state % h
			var1 = ipres
			var2 = ienth
		
		ELSEIF (input .eq. eos_input_th) THEN
			
			single_iter = .true.
			v_want = state % h
			var  = ienth
			dvar = idens
			
		ENDIF
		
		ptot_row = 0.0d0
		dpt_row = 0.0d0
		dpd_row = 0.0d0
		dpa_row = 0.0d0
		dpz_row = 0.0d0
		dpe_row = 0.0d0
		dpdr_e_row = 0.0d0
		
		etot_row = 0.0d0
		det_row = 0.0d0
		ded_row = 0.0d0
		dea_row = 0.0d0
		dez_row = 0.0d0
		
		stot_row = 0.0d0
		dst_row = 0.0d0
		dsd_row = 0.0d0
		
		htot_row = 0.0d0
		dhd_row = 0.0d0
		dht_row = 0.0d0
		
		pele_row = 0.0d0
		ppos_row = 0.0d0
		
		xne_row = 0.0d0
		xnp_row = 0.0d0
		
		etaele_row = 0.0d0
		detadt_row = 0.0d0
		
		cv_row = 0.0d0
		cp_row = 0.0d0
		cs_row = 0.0d0
		gam1_row = 0.0d0
		
		converged = .false.
		
		IF (input .eq. eos_input_rt) converged = .true.
		
		DO iter = 1, max_newton
			
			temp  = temp_row
			den   =  den_row
			abar  = abar_row
			zbar  = zbar_row
			
			ytot1 = 1.0d0 / abar
			ye    = ye_row
			din   = ye * den
			
			!..initialize
			deni    = 1.0d0/den
			tempi   = 1.0d0/temp
			kt      = kerg * temp
			ktinv   = 1.0d0/kt
			
			!..radiation section:
			prad    = asoli3 * temp * temp * temp * temp
			dpraddd = 0.0d0
			dpraddt = 4.0d0 * prad*tempi
			
			erad    = 3.0d0 * prad*deni
			deraddd = -erad*deni
			deraddt = 3.0d0 * dpraddt*deni
			
			srad    = (prad*deni + erad)*tempi
			dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
			dsraddt = (dpraddt*deni + deraddt - srad)*tempi
			
			!..ion section:
			xni     = avo_eos * ytot1 * den
			dxnidd  = avo_eos * ytot1
			dxnida  = -xni * ytot1
			
			pion    = xni * kt
			dpiondd = dxnidd * kt
			dpiondt = xni * kerg
			
			eion    = 1.5d0 * pion*deni
			deiondd = (1.5d0 * dpiondd - eion)*deni
			deiondt = 1.5d0 * dpiondt*deni
			
			x       = abar*abar*SQRT(abar) * deni/avo_eos
			s       = sioncon * temp
			z       = x * s * SQRT(s)
			y       = LOG(z)
			sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
			dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
			- kergavo * deni * ytot1
			dsiondt = (dpiondt*deni + deiondt)*tempi -  &
			(pion*deni + eion) * tempi*tempi  &
			+ 1.5d0 * kergavo * tempi*ytot1
			x       = avo_eos*kerg/abar
			
			! IF (.NOT. include_ion_contribution) THEN
				! pion = zero
				! eion = zero
				! sion = zero

				! dpiondd = zero
				! dpiondt = zero
				! deiondd = zero
				! deiondt = zero

				! dsiondd = zero
				! dsiondt = zero
			! END IF
				
			!..electron-positron section:
			!..assume complete ionization
			xnem    = xni * zbar
			
			!..enter the table with ye*den
			din = ye*den
						
			!..hash locate this temperature and density
			jat = int((log10(temp) - LOG10(smallt))*tstpi) + 1
			jat = MAX(1,MIN(jat,HelmTable % nPointsTemp-1))
			iat = int((log10(din) - LOG10(smalld))*dstpi) + 1
			iat = MAX(1,MIN(iat,HelmTable % nPointsDen-1))
						
			fi(:) = 0.0d0
			!..access the table locations only once
			fi(1)  = HelmTable % f(iat,jat)
			fi(2)  = HelmTable % f(iat+1,jat)
			fi(3)  = HelmTable % f(iat,jat+1)
			fi(4)  = HelmTable % f(iat+1,jat+1)
			fi(5)  = HelmTable % ft(iat,jat)
			fi(6)  = HelmTable % ft(iat+1,jat)
			fi(7)  = HelmTable % ft(iat,jat+1)
			fi(8)  = HelmTable % ft(iat+1,jat+1)
			fi(9)  = HelmTable % ftt(iat,jat)
			fi(10) = HelmTable % ftt(iat+1,jat)
			fi(11) = HelmTable % ftt(iat,jat+1)
			fi(12) = HelmTable % ftt(iat+1,jat+1)
			fi(13) = HelmTable % fd(iat,jat)
			fi(14) = HelmTable % fd(iat+1,jat)
			fi(15) = HelmTable % fd(iat,jat+1)
			fi(16) = HelmTable % fd(iat+1,jat+1)
			fi(17) = HelmTable % fdd(iat,jat)
			fi(18) = HelmTable % fdd(iat+1,jat)
			fi(19) = HelmTable % fdd(iat,jat+1)
			fi(20) = HelmTable % fdd(iat+1,jat+1)
			fi(21) = HelmTable % fdt(iat,jat)
			fi(22) = HelmTable % fdt(iat+1,jat)
			fi(23) = HelmTable % fdt(iat,jat+1)
			fi(24) = HelmTable % fdt(iat+1,jat+1)
			fi(25) = HelmTable % fddt(iat,jat)
			fi(26) = HelmTable % fddt(iat+1,jat)
			fi(27) = HelmTable % fddt(iat,jat+1)
			fi(28) = HelmTable % fddt(iat+1,jat+1)
			fi(29) = HelmTable % fdtt(iat,jat)
			fi(30) = HelmTable % fdtt(iat+1,jat)
			fi(31) = HelmTable % fdtt(iat,jat+1)
			fi(32) = HelmTable % fdtt(iat+1,jat+1)
			fi(33) = HelmTable % fddtt(iat,jat)
			fi(34) = HelmTable % fddtt(iat+1,jat)
			fi(35) = HelmTable % fddtt(iat,jat+1)
			fi(36) = HelmTable % fddtt(iat+1,jat+1)
			
			!..various differences
			xt  = MAX( (temp - HelmTable % t(jat))*HelmTable % dti(jat), zero)
			xd  = MAX( (din - HelmTable % d(iat))*HelmTable % ddi(iat), zero)
			mxt = 1.0d0 - xt
			mxd = 1.0d0 - xd
						
			!..the six density and six temperature basis functions
			si0t =   psi0(xt)
			si1t =   psi1(xt)*HelmTable % dt(jat)
			si2t =   psi2(xt)*HelmTable % dt2(jat)
			
			si0mt =  psi0(mxt)
			si1mt = -psi1(mxt)*HelmTable % dt(jat)
			si2mt =  psi2(mxt)*HelmTable % dt2(jat)
			
			si0d =   psi0(xd)
			si1d =   psi1(xd)*HelmTable % dd(iat)
			si2d =   psi2(xd)*HelmTable % dd2(iat)
			
			si0md =  psi0(mxd)
			si1md = -psi1(mxd)*HelmTable % dd(iat)
			si2md =  psi2(mxd)*HelmTable % dd2(iat)
			
			!..derivatives of the weight functions
			dsi0t =   dpsi0(xt)*HelmTable % dti(jat)
			dsi1t =   dpsi1(xt)
			dsi2t =   dpsi2(xt)*HelmTable % dt(jat)
			
			dsi0mt = -dpsi0(mxt)*HelmTable % dti(jat)
			dsi1mt =  dpsi1(mxt)
			dsi2mt = -dpsi2(mxt)*HelmTable % dt(jat)
			
			dsi0d =   dpsi0(xd)*HelmTable % ddi(iat)
			dsi1d =   dpsi1(xd)
			dsi2d =   dpsi2(xd)*HelmTable % dd(iat)
			
			dsi0md = -dpsi0(mxd)*HelmTable % ddi(iat)
			dsi1md =  dpsi1(mxd)
			dsi2md = -dpsi2(mxd)*HelmTable % dd(iat)
			
			!..second derivatives of the weight functions
			ddsi0t =   ddpsi0(xt)*HelmTable % dt2i(jat)
			ddsi1t =   ddpsi1(xt)*HelmTable % dti(jat)
			ddsi2t =   ddpsi2(xt)
			
			ddsi0mt =  ddpsi0(mxt)*HelmTable % dt2i(jat)
			ddsi1mt = -ddpsi1(mxt)*HelmTable % dti(jat)
			ddsi2mt =  ddpsi2(mxt)
			
			!     ddsi0d =   ddpsi0(xd)*dd2i(iat)
			!     ddsi1d =   ddpsi1(xd)*HelmTable % ddi(iat)
			!     ddsi2d =   ddpsi2(xd)
			
			!     ddsi0md =  ddpsi0(mxd)*dd2i(iat)
			!     ddsi1md = -ddpsi1(mxd)*HelmTable % ddi(iat)
			!     ddsi2md =  ddpsi2(mxd)
			
			!..the free energy
			free  = h5( fi, &
				si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
				si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

			!..derivative with respect to density
			df_d  = h5( fi, &
				si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
				dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

			!..derivative with respect to temperature
			df_t = h5( fi, &
				dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
				si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

			!..derivative with respect to density**2
			!     df_dd = h5( &
			!               si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
			!               ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

			!..derivative with respect to temperature**2
			df_tt = h5( fi, &
				ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
				si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

			!..derivative with respect to temperature and density
			df_dt = h5( fi, &
				dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
				dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)
			
			!..now get the pressure derivative with density, chemical potential, and
			!..electron positron number densities
			!..get the interpolation weight functions
			si0t   =  xpsi0(xt)
			si1t   =  xpsi1(xt)*HelmTable % dt(jat)
			
			si0mt  =  xpsi0(mxt)
			si1mt  =  -xpsi1(mxt)*HelmTable % dt(jat)
			
			si0d   =  xpsi0(xd)
			si1d   =  xpsi1(xd)*HelmTable % dd(iat)
			
			si0md  =  xpsi0(mxd)
			si1md  =  -xpsi1(mxd)*HelmTable % dd(iat)
			
			!..derivatives of weight functions
			dsi0t  = xdpsi0(xt)*HelmTable % dti(jat)
			dsi1t  = xdpsi1(xt)
			
			dsi0mt = -xdpsi0(mxt)*HelmTable % dti(jat)
			dsi1mt = xdpsi1(mxt)
			
			dsi0d  = xdpsi0(xd)*HelmTable % ddi(iat)
			dsi1d  = xdpsi1(xd)
			
			dsi0md = -xdpsi0(mxd)*HelmTable % ddi(iat)
			dsi1md = xdpsi1(mxd)
			
			!..look in the pressure derivative only once
			fi(1)  = HelmTable % dpdf(iat,jat)
			fi(2)  = HelmTable % dpdf(iat+1,jat)
			fi(3)  = HelmTable % dpdf(iat,jat+1)
			fi(4)  = HelmTable % dpdf(iat+1,jat+1)
			fi(5)  = HelmTable % dpdft(iat,jat)
			fi(6)  = HelmTable % dpdft(iat+1,jat)
			fi(7)  = HelmTable % dpdft(iat,jat+1)
			fi(8)  = HelmTable % dpdft(iat+1,jat+1)
			fi(9)  = HelmTable % dpdfd(iat,jat)
			fi(10) = HelmTable % dpdfd(iat+1,jat)
			fi(11) = HelmTable % dpdfd(iat,jat+1)
			fi(12) = HelmTable % dpdfd(iat+1,jat+1)
			fi(13) = HelmTable % dpdfdt(iat,jat)
			fi(14) = HelmTable % dpdfdt(iat+1,jat)
			fi(15) = HelmTable % dpdfdt(iat,jat+1)
			fi(16) = HelmTable % dpdfdt(iat+1,jat+1)
			
			!..pressure derivative with density
			dpepdd  = h3(   fi, &
			si0t,   si1t,   si0mt,   si1mt, &
			si0d,   si1d,   si0md,   si1md)
			dpepdd  = MAX(ye * dpepdd,zero)
			
			!..look in the electron chemical potential table only once
			fi(1)  = HelmTable % ef(iat,jat)
			fi(2)  = HelmTable % ef(iat+1,jat)
			fi(3)  = HelmTable % ef(iat,jat+1)
			fi(4)  = HelmTable % ef(iat+1,jat+1)
			fi(5)  = HelmTable % eft(iat,jat)
			fi(6)  = HelmTable % eft(iat+1,jat)
			fi(7)  = HelmTable % eft(iat,jat+1)
			fi(8)  = HelmTable % eft(iat+1,jat+1)
			fi(9)  = HelmTable % efd(iat,jat)
			fi(10) = HelmTable % efd(iat+1,jat)
			fi(11) = HelmTable % efd(iat,jat+1)
			fi(12) = HelmTable % efd(iat+1,jat+1)
			fi(13) = HelmTable % efdt(iat,jat)
			fi(14) = HelmTable % efdt(iat+1,jat)
			fi(15) = HelmTable % efdt(iat,jat+1)
			fi(16) = HelmTable % efdt(iat+1,jat+1)
			
			!..electron chemical potential etaele
			etaele  = h3( fi, &
			si0t,   si1t,   si0mt,   si1mt, &
			si0d,   si1d,   si0md,   si1md)
			
			!..derivative with respect to density
			x       = h3( fi, &
			si0t,   si1t,   si0mt,   si1mt, &
			dsi0d,  dsi1d,  dsi0md,  dsi1md)
			detadd  = ye * x
			
			!..derivative with respect to temperature
			detadt  = h3( fi, &
			dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
			si0d,   si1d,   si0md,   si1md)
			
			
			!..look in the number density table only once
			fi(1)  = HelmTable % xf(iat,jat)
			fi(2)  = HelmTable % xf(iat+1,jat)
			fi(3)  = HelmTable % xf(iat,jat+1)
			fi(4)  = HelmTable % xf(iat+1,jat+1)
			fi(5)  = HelmTable % xft(iat,jat)
			fi(6)  = HelmTable % xft(iat+1,jat)
			fi(7)  = HelmTable % xft(iat,jat+1)
			fi(8)  = HelmTable % xft(iat+1,jat+1)
			fi(9)  = HelmTable % xfd(iat,jat)
			fi(10) = HelmTable % xfd(iat+1,jat)
			fi(11) = HelmTable % xfd(iat,jat+1)
			fi(12) = HelmTable % xfd(iat+1,jat+1)
			fi(13) = HelmTable % xfdt(iat,jat)
			fi(14) = HelmTable % xfdt(iat+1,jat)
			fi(15) = HelmTable % xfdt(iat,jat+1)
			fi(16) = HelmTable % xfdt(iat+1,jat+1)
			
			!..electron + positron number densities
			xnefer   = h3( fi, &
			si0t,   si1t,   si0mt,   si1mt, &
			si0d,   si1d,   si0md,   si1md)
			
			!..derivative with respect to density
			x        = h3( fi, &
			si0t,   si1t,   si0mt,   si1mt, &
			dsi0d,  dsi1d,  dsi0md,  dsi1md)
			x = MAX(x,zero)
			dxnedd   = ye * x
			
			!..derivative with respect to temperature
			dxnedt   = h3( fi, &
			dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
			si0d,   si1d,   si0md,   si1md)
			
			!..the desired electron-positron thermodynamic quantities
			
			!..dpepdd at high temperatures and low densities is below the
			!..floating point limit of the subtraction of two large terms.
			!..since dpresdd doesn't enter the maxwell relations at all, use the
			!..bicubic interpolation done above instead of this one
			x       = din * din
			pele    = x * df_d
			dpepdt  = x * df_dt
			!     dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
			s       = dpepdd/ye - 2.0d0 * din * df_d
			
			x       = ye * ye
			sele    = -df_t * ye
			dsepdt  = -df_tt * ye
			dsepdd  = -df_dt * x
			
			eele    = ye*free + temp * sele
			deepdt  = temp * dsepdt
			deepdd  = x * df_d + temp * dsepdd
			
			!..coulomb section:
			!..initialize
			pcoul    = 0.0d0
			dpcouldd = 0.0d0
			dpcouldt = 0.0d0
			dpcoulda = 0.0d0
			dpcouldz = 0.0d0
			ecoul    = 0.0d0
			decouldd = 0.0d0
			decouldt = 0.0d0
			decoulda = 0.0d0
			decouldz = 0.0d0
			scoul    = 0.0d0
			dscouldd = 0.0d0
			dscouldt = 0.0d0
			dscoulda = 0.0d0
			dscouldz = 0.0d0
			
			!..uniform background corrections only
			!..from yakovlev & shalybkov 1989
			!..lami is the average ion seperation
			!..plasg is the plasma coupling PARAMETER
			z        = forth * pi
			s        = z * xni
			dsdd     = z * dxnidd
			dsda     = z * dxnida
			
			lami     = 1.0d0/s**onethird
			inv_lami = 1.0d0/lami
			z        = -onethird * lami
			lamidd   = z * dsdd/s
			lamida   = z * dsda/s
			
			plasg    = zbar*zbar*esqu*ktinv*inv_lami
			z        = -plasg * inv_lami
			plasgdd  = z * lamidd
			plasgda  = z * lamida
			plasgdt  = -plasg*ktinv * kerg
			plasgdz  = 2.0d0 * plasg/zbar
			
			! TURN ON/OFF COULOMB
			IF ( do_coulomb ) THEN
				!...yakovlev & shalybkov 1989 equations 82, 85, 86, 87
				IF (plasg .ge. 1.0D0) THEN
					x        = plasg**(0.25d0)
					y        = avo_eos * ytot1 * kerg
					ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
					pcoul    = onethird * den * ecoul
					scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
					+ d1 * (LOG(plasg) - 1.0d0) - e1)
					
					y        = avo_eos*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
					decouldd = y * plasgdd
					decouldt = y * plasgdt + ecoul/temp
					decoulda = y * plasgda - ecoul/abar
					decouldz = y * plasgdz
					
					y        = onethird * den
					dpcouldd = onethird * ecoul + y*decouldd
					dpcouldt = y * decouldt
					dpcoulda = y * decoulda
					dpcouldz = y * decouldz
					
					y        = -avo_eos*kerg/(abar*plasg)* &
					(0.75d0*b1*x+1.25d0*c1/x+d1)
					dscouldd = y * plasgdd
					dscouldt = y * plasgdt
					dscoulda = y * plasgda - scoul/abar
					dscouldz = y * plasgdz
					
				!...yakovlev & shalybkov 1989 equations 102, 103, 104
				ELSE IF (plasg .lt. 1.0D0) THEN
					x        = plasg*SQRT(plasg)
					y        = plasg**b2
					z        = c2 * x - onethird * a2 * y
					pcoul    = -pion * z
					ecoul    = 3.0d0 * pcoul/den
					scoul    = -avo_eos/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)
					
					s        = 1.5d0*c2*x/plasg - onethird*a2*b2*y/plasg
					dpcouldd = -dpiondd*z - pion*s*plasgdd
					dpcouldt = -dpiondt*z - pion*s*plasgdt
					
					s        = 3.0d0/den
					decouldd = s * dpcouldd - ecoul/den
					decouldt = s * dpcouldt
					decoulda = s * dpcoulda
					decouldz = s * dpcouldz
					
					s        = -avo_eos*kerg/(abar*plasg)* &
					(1.5d0*c2*x-a2*(b2-1.0d0)*y)
					dscouldd = s * plasgdd
					dscouldt = s * plasgdt
					dscoulda = s * plasgda - scoul/abar
					dscouldz = s * plasgdz
				END IF
				
				IF (include_ion_contribution) THEN
					p_temp = prad + pion + pele + pcoul
					e_temp = erad + eion + eele + ecoul
				ELSE
					p_temp = prad + pele + pcoul
					e_temp = erad + eele + ecoul
				END IF
				
				! Disable Coulomb corrections IF they cause
				! the energy or pressure to go negative.				
				IF (p_temp .le. ZERO .or. e_temp .le. ZERO) THEN
					
					pcoul    = 0.0d0
					dpcouldd = 0.0d0
					dpcouldt = 0.0d0
					dpcoulda = 0.0d0
					dpcouldz = 0.0d0
					ecoul    = 0.0d0
					decouldd = 0.0d0
					decouldt = 0.0d0
					decoulda = 0.0d0
					decouldz = 0.0d0
					scoul    = 0.0d0
					dscouldd = 0.0d0
					dscouldt = 0.0d0
					dscoulda = 0.0d0
					dscouldz = 0.0d0
					
				END IF
			END IF
			
			!..sum all the components
			IF (include_ion_contribution) THEN
				pres    = prad + pion + pele + pcoul
				ener    = erad + eion + eele + ecoul
				entr    = srad + sion + sele + scoul
				
				dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd
				dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
				denerdd = deraddd + deiondd + deepdd + decouldd
				denerdt = deraddt + deiondt + deepdt + decouldt
				
				dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
				dentrdt = dsraddt + dsiondt + dsepdt + dscouldt
			ELSE
				pres    = prad + pele + pcoul
				ener    = erad + eele + ecoul
				entr    = srad + sele + scoul
				
				dpresdd = dpraddd + dpepdd + dpcouldd
				dpresdt = dpraddt + dpepdt + dpcouldt
				denerdd = deraddd + deepdd + decouldd
				denerdt = deraddt + deepdt + decouldt
				
				dentrdd = dsraddd + dsepdd + dscouldd
				dentrdt = dsraddt + dsepdt + dscouldt
			END IF
			
			!..the temperature and density exponents (c&g 9.81 9.82)
				!..the specific heat at constant volume (c&g 9.92)
				!..the third adiabatic exponent (c&g 9.93)
				!..the first adiabatic exponent (c&g 9.97)
				!..the second adiabatic exponent (c&g 9.105)
				!..the specific heat at constant pressure (c&g 9.98)
			!..and relativistic formula for the sound speed (c&g 14.29)
			zz    = pres*deni
			zzi   = den/pres
			chit  = temp/pres * dpresdt
			chid  = dpresdd*zzi
			cv    = denerdt
			x     = zz * chit/(temp * cv)
			gam3  = x + 1.0d0
			gam1  = chit*x + chid
			nabad = x/gam1
			gam2  = 1.0d0/(1.0d0 - nabad)
			cp    = cv * gam1/chid
			z     = 1.0d0 + (ener + light2)*zzi
			sound = clight * SQRT(gam1/z)
			
			!..maxwell relations; each is zero IF the consistency is perfect
			x   = den * den
			dse = temp*dentrdt/denerdt - 1.0d0
			dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
			dsp = -dentrdd*x/dpresdt - 1.0d0
			
			ptot_row = pres
			dpt_row = dpresdt
			dpd_row = dpresdd
			dpe_row = dpresdt / denerdt
			dpdr_e_row = dpresdd - dpresdt * denerdd / denerdt
			
			etot_row = ener
			det_row = denerdt
			ded_row = denerdd
			
			stot_row = entr
			dst_row = dentrdt
			dsd_row = dentrdd
			
			htot_row = ener + pres / den
			dhd_row = denerdd + dpresdd / den - pres / den**2
			dht_row = denerdt + dpresdt / den
			
			pele_row = pele
			ppos_row = 0.0d0
			
			xne_row = xnefer
			xnp_row = 0.0d0
			
			etaele_row = etaele
			detadt_row = detadt
			
			cv_row = cv
			cp_row = cp
			cs_row = sound
			gam1_row = gam1
			
			IF (converged) THEN
				
				EXIT
			
			ELSEIF (single_iter) THEN
				
				IF (dvar .eq. itemp) THEN
					
					x = temp_row
					smallx = smallt
					xtol = ttol
					
					IF (var .eq. ipres) THEN
						v    = ptot_row
					dvdx = dpt_row
					ELSEIF (var .eq. iener) THEN
						v    = etot_row
					dvdx = det_row
					ELSEIF (var .eq. ientr) THEN
						v    = stot_row
					dvdx = dst_row
					ELSEIF (var .eq. ienth) THEN
						v    = htot_row
					dvdx = dht_row
					ELSE
						EXIT
					ENDIF
				
				ELSE ! dvar == density
					
					x = den_row
					smallx = smalld
					xtol = dtol
					
					IF (var .eq. ipres) THEN
						v    = ptot_row
					dvdx = dpd_row
					ELSEIF (var .eq. iener) THEN
						v    = etot_row
					dvdx = ded_row
					ELSEIF (var .eq. ientr) THEN
						v    = stot_row
					dvdx = dsd_row
					ELSEIF (var .eq. ienth) THEN
						v    = htot_row
					dvdx = dhd_row
					ELSE
						EXIT
					ENDIF
					
				ENDIF
				
				! Now DO the calculation for the next guess for T/rho
				
				xnew = x - (v - v_want) / dvdx
				
				! Don't let the temperature/density change by more than a factor of two
				xnew = MAX(0.5 * x, MIN(xnew, 2.0 * x))
				
				! Don't let us freeze/evacuate
				xnew = MAX(smallx, xnew)
				
				! Store the new temperature/density
				
				IF (dvar .eq. itemp) THEN
				temp_row = xnew
				ELSE
					den_row  = xnew
				ENDIF
				
				! Compute the error from the last iteration
				
				error = ABS( (xnew - x) / x )
				
				IF (error .lt. xtol) converged = .true.
			
			ELSEIF (double_iter) THEN
				
				! figure out which variables we're using
				
				told = temp_row
				rold = den_row
				
				IF (var1 .eq. ipres) THEN
					v1    = ptot_row
					dv1dt = dpt_row
				dv1dr = dpd_row
				ELSEIF (var1 .eq. iener) THEN
					v1    = etot_row
					dv1dt = det_row
				dv1dr = ded_row
				ELSEIF (var1 .eq. ientr) THEN
					v1    = stot_row
					dv1dt = dst_row
				dv1dr = dsd_row
				ELSEIF (var1 .eq. ienth) THEN
					v1    = htot_row
					dv1dt = dht_row
				dv1dr = dhd_row
				ELSE
					EXIT
				ENDIF
				
				IF (var2 .eq. ipres) THEN
					v2    = ptot_row
					dv2dt = dpt_row
				dv2dr = dpd_row
				ELSEIF (var2 .eq. iener) THEN
					v2    = etot_row
					dv2dt = det_row
				dv2dr = ded_row
				ELSEIF (var2 .eq. ientr) THEN
					v2    = stot_row
					dv2dt = dst_row
				dv2dr = dsd_row
				ELSEIF (var2 .eq. ienth) THEN
					v2    = htot_row
					dv2dt = dht_row
				dv2dr = dhd_row
				ELSE
					EXIT
				ENDIF
				
				! Two functions, f and g, to iterate over
				v1i = v1_want - v1
				v2i = v2_want - v2
				
				!
					! 0 = f + dfdr * delr + dfdt * delt
					! 0 = g + dgdr * delr + dgdt * delt
				!
				
				! note that dfi/dT = - df/dT
				delr = (-v1i*dv2dt + v2i*dv1dt) / (dv2dr*dv1dt - dv2dt*dv1dr)
				
				rnew = rold + delr
				
				tnew = told + (v1i - dv1dr*delr) / dv1dt
				
				! Don't let the temperature or density change by more
				! than a factor of two
				tnew = MAX(HALF * told, MIN(tnew, 2.0d0 * told))
				rnew = MAX(HALF * rold, MIN(rnew, 2.0d0 * rold))
				
				! Don't let us freeze or evacuate
				tnew = MAX(smallt, tnew)
				rnew = MAX(smalld, rnew)
				
				! Store the new temperature and density
				den_row  = rnew
				temp_row = tnew
				
				! Compute the errors
				error1 = ABS( (rnew - rold) / rold )
				error2 = ABS( (tnew - told) / told )
				
				IF (error1 .LT. dtol .and. error2 .LT. ttol) converged = .true.
				
			ENDIF
			
		ENDDO
		
		state % T    = temp_row
		state % rho  = den_row
		
		state % p    = ptot_row
		state % dpdT = dpt_row
		state % dpdr = dpd_row
		
		
		state % dpde = dpe_row
		state % dpdr_e = dpdr_e_row
		
		state % e    = etot_row
		state % dedT = det_row
		state % dedr = ded_row
		
		state % s    = stot_row / (kmev * ergmev / rmu)
		state % dsdT = dst_row / (kmev * ergmev / rmu)
		state % dsdr = dsd_row / (kmev * ergmev / rmu)
		
		state % h    = htot_row
		state % dhdR = dhd_row
		state % dhdT = dht_row
		
		state % pele = pele_row
		state % ppos = ppos_row
		
		state % xne = xne_row
		state % xnp = xnp_row
		
		state % eta = etaele_row
		state % detadt = detadt_row
		
		state % cv   = cv_row
		state % cp   = cp_row
		state % gam1 = gam1_row
		state % cs   = cs_row
		
		! Take care of final housekeeping.
		
		! Count the positron contribution in the electron quantities.
		
		state % xne  = state % xne  + state % xnp
		state % pele = state % pele + state % ppos
		
		! Use the non-relativistic version of the sound speed, cs = SQRT(gam_1 * P / rho).
		! This replaces the relativistic version that comes out of helmeos.
		
		! state % cs = SQRT(state % gam1 * state % p / state % rho)
		
		IF (input_is_constant) THEN
			
			IF (input .eq. eos_input_rh) THEN
				
				state % h = v_want
			
			ELSEIF (input .eq. eos_input_tp) THEN
				
				state % p = v_want
			
			ELSEIF (input .eq. eos_input_rp) THEN
				
				state % p = v_want
			
			ELSEIF (input .eq. eos_input_re) THEN
				
				state % e = v_want
			
			ELSEIF (input .eq. eos_input_ps) THEN
				
				state % p = v1_want
				state % s = v2_want / (kmev * ergmev / rmu)
			
			ELSEIF (input .eq. eos_input_ph) THEN
				
				state % p = v1_want
				state % h = v2_want
			
			ELSEIF (input .eq. eos_input_th) THEN
				
				state % h = v_want
				
			ENDIF
			
		ENDIF
		
	END SUBROUTINE FullHelmEOS
	
   SUBROUTINE HelmEOSComponents(input, HelmTable, HelmState)
	   
		! In principle this should make the calculation very efficient, and only calculate electron, positrons, and photons,
		! and ions only if needed. But it's not complete yet, please disregard
        !..input arguments
        INTEGER, INTENT(IN) :: input
        TYPE(HelmholtzEOSType), INTENT(IN) :: HelmTable
        TYPE (HelmholtzStateType), INTENT(INOUT) :: HelmState
				
		!..declare local variables
		REAL(dp) :: x,y,z,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
		dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
		dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
		deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
		dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
		sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
		dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
		gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
		detadt,detadd,xnefer,dxnedt,dxnedd,s, &
		temp,den,abar,zbar,ytot1,ye,din
		
		
		!..for the interpolations
		INTEGER  :: iat,jat
		REAL(dp) :: free,df_d,df_t,df_tt,df_dt
		REAL(dp) :: xt,xd,mxt,mxd,fi(36), &
		si0t,si1t,si2t,si0mt,si1mt,si2mt, &
		si0d,si1d,si2d,si0md,si1md,si2md, &
		dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
		dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
		ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt
		
		!..for the coulomb corrections
		REAL(dp) :: dsdd,dsda,lami,inv_lami,lamida,lamidd,     &
		plasg,plasgdd,plasgdt,plasgda,plasgdz,     &
		ecoul,decouldd,decouldt,decoulda,decouldz, &
		pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
		scoul,dscouldd,dscouldt,dscoulda,dscouldz
				
		REAL(dp) :: p_temp, e_temp
		
		! extra variables I added that were previously allocated
		REAL(dp) :: smallt, smalld, hight, highd
		REAL(dp) :: tstp, dstp, dstpi, tstpi
		
		smallt = HelmTable % mintemp
		smalld = HelmTable % mindens
		hight = HelmTable % maxtemp
		highd = HelmTable % maxdens
				
		tstp  = (LOG10(hight) - LOG10(smallt))/float(HelmTable % nPointsTemp -1)
		dstp  = (LOG10(highd) - LOG10(smalld))/float(HelmTable % nPointsDen -1)
		
		tstpi = one / tstp
		dstpi = one / dstp
		
		IF (.not. (input .eq. eos_input_rt)) THEN
			WRITE(*,*) 'This routine is not optimized for anything other than temperature and density as inputs'
			STOP
		ENDIF 

		temp = HelmState % T
		den  = HelmState % rho
		abar = HelmState % abar
		zbar = HelmState % zbar
		ye   = HelmState % ye

		ytot1 = 1.0d0 / abar
		din   = ye * den
		
		!..initialize
		deni    = 1.0d0/den
		tempi   = 1.0d0/temp
		kt      = kerg * temp
		ktinv   = 1.0d0/kt
		
		!..radiation section:
		prad    = asoli3 * temp * temp * temp * temp
		dpraddd = 0.0d0
		dpraddt = 4.0d0 * prad*tempi
		
		erad    = 3.0d0 * prad*deni
		deraddd = -erad*deni
		deraddt = 3.0d0 * dpraddt*deni
		
		srad    = (prad*deni + erad)*tempi
		dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
		dsraddt = (dpraddt*deni + deraddt - srad)*tempi
		
		!..ion section:
		xni     = avo_eos * ytot1 * den
		dxnidd  = avo_eos * ytot1
		dxnida  = -xni * ytot1
		
		pion    = xni * kt
		dpiondd = dxnidd * kt
		dpiondt = xni * kerg
		
		eion    = 1.5d0 * pion*deni
		deiondd = (1.5d0 * dpiondd - eion)*deni
		deiondt = 1.5d0 * dpiondt*deni
		
		x       = abar*abar*SQRT(abar) * deni/avo_eos
		s       = sioncon * temp
		z       = x * s * SQRT(s)
		y       = LOG(z)
		sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
		dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
		- kergavo * deni * ytot1
		dsiondt = (dpiondt*deni + deiondt)*tempi -  &
		(pion*deni + eion) * tempi*tempi  &
		+ 1.5d0 * kergavo * tempi*ytot1
		x       = avo_eos*kerg/abar
			
		!..electron-positron section:
		!..assume complete ionization
		xnem    = xni * zbar
		
		!..enter the table with ye*den
		din = ye*den
					
		!..hash locate this temperature and density
		jat = int((log10(temp) - LOG10(smallt))*tstpi) + 1
		jat = MAX(1,MIN(jat,HelmTable % nPointsTemp-1))
		iat = int((log10(din) - LOG10(smalld))*dstpi) + 1
		iat = MAX(1,MIN(iat,HelmTable % nPointsDen-1))
					
		fi(:) = 0.0d0
		!..access the table locations only once
		fi(1)  = HelmTable % f(iat,jat)
		fi(2)  = HelmTable % f(iat+1,jat)
		fi(3)  = HelmTable % f(iat,jat+1)
		fi(4)  = HelmTable % f(iat+1,jat+1)
		fi(5)  = HelmTable % ft(iat,jat)
		fi(6)  = HelmTable % ft(iat+1,jat)
		fi(7)  = HelmTable % ft(iat,jat+1)
		fi(8)  = HelmTable % ft(iat+1,jat+1)
		fi(9)  = HelmTable % ftt(iat,jat)
		fi(10) = HelmTable % ftt(iat+1,jat)
		fi(11) = HelmTable % ftt(iat,jat+1)
		fi(12) = HelmTable % ftt(iat+1,jat+1)
		fi(13) = HelmTable % fd(iat,jat)
		fi(14) = HelmTable % fd(iat+1,jat)
		fi(15) = HelmTable % fd(iat,jat+1)
		fi(16) = HelmTable % fd(iat+1,jat+1)
		fi(17) = HelmTable % fdd(iat,jat)
		fi(18) = HelmTable % fdd(iat+1,jat)
		fi(19) = HelmTable % fdd(iat,jat+1)
		fi(20) = HelmTable % fdd(iat+1,jat+1)
		fi(21) = HelmTable % fdt(iat,jat)
		fi(22) = HelmTable % fdt(iat+1,jat)
		fi(23) = HelmTable % fdt(iat,jat+1)
		fi(24) = HelmTable % fdt(iat+1,jat+1)
		fi(25) = HelmTable % fddt(iat,jat)
		fi(26) = HelmTable % fddt(iat+1,jat)
		fi(27) = HelmTable % fddt(iat,jat+1)
		fi(28) = HelmTable % fddt(iat+1,jat+1)
		fi(29) = HelmTable % fdtt(iat,jat)
		fi(30) = HelmTable % fdtt(iat+1,jat)
		fi(31) = HelmTable % fdtt(iat,jat+1)
		fi(32) = HelmTable % fdtt(iat+1,jat+1)
		fi(33) = HelmTable % fddtt(iat,jat)
		fi(34) = HelmTable % fddtt(iat+1,jat)
		fi(35) = HelmTable % fddtt(iat,jat+1)
		fi(36) = HelmTable % fddtt(iat+1,jat+1)
		
		!..various differences
		xt  = MAX( (temp - HelmTable % t(jat))*HelmTable % dti(jat), zero)
		xd  = MAX( (din - HelmTable % d(iat))*HelmTable % ddi(iat), zero)
		mxt = 1.0d0 - xt
		mxd = 1.0d0 - xd
					
		!..the six density and six temperature basis functions
		si0t =   psi0(xt)
		si1t =   psi1(xt)*HelmTable % dt(jat)
		si2t =   psi2(xt)*HelmTable % dt2(jat)
		
		si0mt =  psi0(mxt)
		si1mt = -psi1(mxt)*HelmTable % dt(jat)
		si2mt =  psi2(mxt)*HelmTable % dt2(jat)
		
		si0d =   psi0(xd)
		si1d =   psi1(xd)*HelmTable % dd(iat)
		si2d =   psi2(xd)*HelmTable % dd2(iat)
		
		si0md =  psi0(mxd)
		si1md = -psi1(mxd)*HelmTable % dd(iat)
		si2md =  psi2(mxd)*HelmTable % dd2(iat)
		
		!..derivatives of the weight functions
		dsi0t =   dpsi0(xt)*HelmTable % dti(jat)
		dsi1t =   dpsi1(xt)
		dsi2t =   dpsi2(xt)*HelmTable % dt(jat)
		
		dsi0mt = -dpsi0(mxt)*HelmTable % dti(jat)
		dsi1mt =  dpsi1(mxt)
		dsi2mt = -dpsi2(mxt)*HelmTable % dt(jat)
		
		dsi0d =   dpsi0(xd)*HelmTable % ddi(iat)
		dsi1d =   dpsi1(xd)
		dsi2d =   dpsi2(xd)*HelmTable % dd(iat)
		
		dsi0md = -dpsi0(mxd)*HelmTable % ddi(iat)
		dsi1md =  dpsi1(mxd)
		dsi2md = -dpsi2(mxd)*HelmTable % dd(iat)
		
		!..second derivatives of the weight functions
		ddsi0t =   ddpsi0(xt)*HelmTable % dt2i(jat)
		ddsi1t =   ddpsi1(xt)*HelmTable % dti(jat)
		ddsi2t =   ddpsi2(xt)
		
		ddsi0mt =  ddpsi0(mxt)*HelmTable % dt2i(jat)
		ddsi1mt = -ddpsi1(mxt)*HelmTable % dti(jat)
		ddsi2mt =  ddpsi2(mxt)
		
		! ddsi0d =   ddpsi0(xd)*dd2i(iat)
		! ddsi1d =   ddpsi1(xd)*HelmTable % ddi(iat)
		! ddsi2d =   ddpsi2(xd)
	
		! ddsi0md =  ddpsi0(mxd)*dd2i(iat)
		! ddsi1md = -ddpsi1(mxd)*HelmTable % ddi(iat)
		! ddsi2md =  ddpsi2(mxd)
		
		!..the free energy
		free  = h5( fi, &
			si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
			si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

		!..derivative with respect to density
		df_d  = h5( fi, &
			si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
			dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

		!..derivative with respect to temperature
		df_t = h5( fi, &
			dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
			si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

		!..derivative with respect to density**2
		!     df_dd = h5( &
		!               si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
		!               ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

		!..derivative with respect to temperature**2
		df_tt = h5( fi, &
			ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
			si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

		!..derivative with respect to temperature and density
		df_dt = h5( fi, &
			dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
			dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)
		
		!..now get the pressure derivative with density, chemical potential, and
		!..electron positron number densities
		!..get the interpolation weight functions
		si0t   =  xpsi0(xt)
		si1t   =  xpsi1(xt)*HelmTable % dt(jat)
		
		si0mt  =  xpsi0(mxt)
		si1mt  =  -xpsi1(mxt)*HelmTable % dt(jat)
		
		si0d   =  xpsi0(xd)
		si1d   =  xpsi1(xd)*HelmTable % dd(iat)
		
		si0md  =  xpsi0(mxd)
		si1md  =  -xpsi1(mxd)*HelmTable % dd(iat)
		
		!..derivatives of weight functions
		dsi0t  = xdpsi0(xt)*HelmTable % dti(jat)
		dsi1t  = xdpsi1(xt)
		
		dsi0mt = -xdpsi0(mxt)*HelmTable % dti(jat)
		dsi1mt = xdpsi1(mxt)
		
		dsi0d  = xdpsi0(xd)*HelmTable % ddi(iat)
		dsi1d  = xdpsi1(xd)
		
		dsi0md = -xdpsi0(mxd)*HelmTable % ddi(iat)
		dsi1md = xdpsi1(mxd)
		
		!..look in the pressure derivative only once
		fi(1)  = HelmTable % dpdf(iat,jat)
		fi(2)  = HelmTable % dpdf(iat+1,jat)
		fi(3)  = HelmTable % dpdf(iat,jat+1)
		fi(4)  = HelmTable % dpdf(iat+1,jat+1)
		fi(5)  = HelmTable % dpdft(iat,jat)
		fi(6)  = HelmTable % dpdft(iat+1,jat)
		fi(7)  = HelmTable % dpdft(iat,jat+1)
		fi(8)  = HelmTable % dpdft(iat+1,jat+1)
		fi(9)  = HelmTable % dpdfd(iat,jat)
		fi(10) = HelmTable % dpdfd(iat+1,jat)
		fi(11) = HelmTable % dpdfd(iat,jat+1)
		fi(12) = HelmTable % dpdfd(iat+1,jat+1)
		fi(13) = HelmTable % dpdfdt(iat,jat)
		fi(14) = HelmTable % dpdfdt(iat+1,jat)
		fi(15) = HelmTable % dpdfdt(iat,jat+1)
		fi(16) = HelmTable % dpdfdt(iat+1,jat+1)
		
		!..pressure derivative with density
		dpepdd  = h3(   fi, &
		si0t,   si1t,   si0mt,   si1mt, &
		si0d,   si1d,   si0md,   si1md)
		dpepdd  = MAX(ye * dpepdd,zero)
		
		!..look in the electron chemical potential table only once
		fi(1)  = HelmTable % ef(iat,jat)
		fi(2)  = HelmTable % ef(iat+1,jat)
		fi(3)  = HelmTable % ef(iat,jat+1)
		fi(4)  = HelmTable % ef(iat+1,jat+1)
		fi(5)  = HelmTable % eft(iat,jat)
		fi(6)  = HelmTable % eft(iat+1,jat)
		fi(7)  = HelmTable % eft(iat,jat+1)
		fi(8)  = HelmTable % eft(iat+1,jat+1)
		fi(9)  = HelmTable % efd(iat,jat)
		fi(10) = HelmTable % efd(iat+1,jat)
		fi(11) = HelmTable % efd(iat,jat+1)
		fi(12) = HelmTable % efd(iat+1,jat+1)
		fi(13) = HelmTable % efdt(iat,jat)
		fi(14) = HelmTable % efdt(iat+1,jat)
		fi(15) = HelmTable % efdt(iat,jat+1)
		fi(16) = HelmTable % efdt(iat+1,jat+1)
		
		!..electron chemical potential etaele
		etaele  = h3( fi, &
		si0t,   si1t,   si0mt,   si1mt, &
		si0d,   si1d,   si0md,   si1md)
		
		!..derivative with respect to density
		x       = h3( fi, &
		si0t,   si1t,   si0mt,   si1mt, &
		dsi0d,  dsi1d,  dsi0md,  dsi1md)
		detadd  = ye * x
		
		!..derivative with respect to temperature
		detadt  = h3( fi, &
		dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
		si0d,   si1d,   si0md,   si1md)
		
		!..look in the number density table only once
		fi(1)  = HelmTable % xf(iat,jat)
		fi(2)  = HelmTable % xf(iat+1,jat)
		fi(3)  = HelmTable % xf(iat,jat+1)
		fi(4)  = HelmTable % xf(iat+1,jat+1)
		fi(5)  = HelmTable % xft(iat,jat)
		fi(6)  = HelmTable % xft(iat+1,jat)
		fi(7)  = HelmTable % xft(iat,jat+1)
		fi(8)  = HelmTable % xft(iat+1,jat+1)
		fi(9)  = HelmTable % xfd(iat,jat)
		fi(10) = HelmTable % xfd(iat+1,jat)
		fi(11) = HelmTable % xfd(iat,jat+1)
		fi(12) = HelmTable % xfd(iat+1,jat+1)
		fi(13) = HelmTable % xfdt(iat,jat)
		fi(14) = HelmTable % xfdt(iat+1,jat)
		fi(15) = HelmTable % xfdt(iat,jat+1)
		fi(16) = HelmTable % xfdt(iat+1,jat+1)
		
		!..electron + positron number densities
		xnefer   = h3( fi, &
		si0t,   si1t,   si0mt,   si1mt, &
		si0d,   si1d,   si0md,   si1md)
		
		!..derivative with respect to density
		x        = h3( fi, &
		si0t,   si1t,   si0mt,   si1mt, &
		dsi0d,  dsi1d,  dsi0md,  dsi1md)
		x = MAX(x,zero)
		dxnedd   = ye * x
		
		!..derivative with respect to temperature
		dxnedt   = h3( fi, &
		dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
		si0d,   si1d,   si0md,   si1md)
		
		!..the desired electron-positron thermodynamic quantities
		
		!..dpepdd at high temperatures and low densities is below the
		!..floating point limit of the subtraction of two large terms.
		!..since dpresdd doesn't enter the maxwell relations at all, use the
		!..bicubic interpolation done above instead of this one
		x       = din * din
		pele    = x * df_d
		dpepdt  = x * df_dt
		!     dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
		s       = dpepdd/ye - 2.0d0 * din * df_d
		
		x       = ye * ye
		sele    = -df_t * ye
		dsepdt  = -df_tt * ye
		dsepdd  = -df_dt * x
		
		eele    = ye*free + temp * sele
		deepdt  = temp * dsepdt
		deepdd  = x * df_d + temp * dsepdd
		
		!..coulomb section:
		!..initialize
		pcoul    = 0.0d0
		dpcouldd = 0.0d0
		dpcouldt = 0.0d0
		dpcoulda = 0.0d0
		dpcouldz = 0.0d0
		ecoul    = 0.0d0
		decouldd = 0.0d0
		decouldt = 0.0d0
		decoulda = 0.0d0
		decouldz = 0.0d0
		scoul    = 0.0d0
		dscouldd = 0.0d0
		dscouldt = 0.0d0
		dscoulda = 0.0d0
		dscouldz = 0.0d0
		
		!..uniform background corrections only
		!..from yakovlev & shalybkov 1989
		!..lami is the average ion seperation
		!..plasg is the plasma coupling PARAMETER
		z        = forth * pi
		s        = z * xni
		dsdd     = z * dxnidd
		dsda     = z * dxnida
		
		lami     = 1.0d0/s**onethird
		inv_lami = 1.0d0/lami
		z        = -onethird * lami
		lamidd   = z * dsdd/s
		lamida   = z * dsda/s
		
		plasg    = zbar*zbar*esqu*ktinv*inv_lami
		z        = -plasg * inv_lami
		plasgdd  = z * lamidd
		plasgda  = z * lamida
		plasgdt  = -plasg*ktinv * kerg
		plasgdz  = 2.0d0 * plasg/zbar
		
		! Calcuate Coulomb correction, it's up to the user to decide if you need to include them
		
		!...yakovlev & shalybkov 1989 equations 82, 85, 86, 87
		IF (plasg .ge. 1.0D0) THEN
			x        = plasg**(0.25d0)
			y        = avo_eos * ytot1 * kerg
			ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
			pcoul    = onethird * den * ecoul
			scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
			+ d1 * (LOG(plasg) - 1.0d0) - e1)
			
			y        = avo_eos*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
			decouldd = y * plasgdd
			decouldt = y * plasgdt + ecoul/temp
			decoulda = y * plasgda - ecoul/abar
			decouldz = y * plasgdz
			
			y        = onethird * den
			dpcouldd = onethird * ecoul + y*decouldd
			dpcouldt = y * decouldt
			dpcoulda = y * decoulda
			dpcouldz = y * decouldz
			
			y        = -avo_eos*kerg/(abar*plasg)* &
			(0.75d0*b1*x+1.25d0*c1/x+d1)
			dscouldd = y * plasgdd
			dscouldt = y * plasgdt
			dscoulda = y * plasgda - scoul/abar
			dscouldz = y * plasgdz
			
		!...yakovlev & shalybkov 1989 equations 102, 103, 104
		ELSE IF (plasg .lt. 1.0D0) THEN
			x        = plasg*SQRT(plasg)
			y        = plasg**b2
			z        = c2 * x - onethird * a2 * y
			pcoul    = -pion * z
			ecoul    = 3.0d0 * pcoul/den
			scoul    = -avo_eos/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)
			
			s        = 1.5d0*c2*x/plasg - onethird*a2*b2*y/plasg
			dpcouldd = -dpiondd*z - pion*s*plasgdd
			dpcouldt = -dpiondt*z - pion*s*plasgdt
			
			s        = 3.0d0/den
			decouldd = s * dpcouldd - ecoul/den
			decouldt = s * dpcouldt
			decoulda = s * dpcoulda
			decouldz = s * dpcouldz
			
			s        = -avo_eos*kerg/(abar*plasg)* &
			(1.5d0*c2*x-a2*(b2-1.0d0)*y)
			dscouldd = s * plasgdd
			dscouldt = s * plasgdt
			dscoulda = s * plasgda - scoul/abar
			dscouldz = s * plasgdz
		END IF
	
		p_temp = prad + pion + pele + pcoul
		e_temp = erad + eion + eele + ecoul
		
		! Disable Coulomb corrections IF they cause
		! the energy or pressure to go negative.				
		IF (p_temp .le. ZERO .or. e_temp .le. ZERO) THEN
			
			pcoul    = 0.0d0
			dpcouldd = 0.0d0
			dpcouldt = 0.0d0
			dpcoulda = 0.0d0
			dpcouldz = 0.0d0
			ecoul    = 0.0d0
			decouldd = 0.0d0
			decouldt = 0.0d0
			decoulda = 0.0d0
			decouldz = 0.0d0
			scoul    = 0.0d0
			dscouldd = 0.0d0
			dscouldt = 0.0d0
			dscoulda = 0.0d0
			dscouldz = 0.0d0
			
		END IF

		! pres    = prad + pion + pele + pcoul
		! ener    = erad + eion + eele + ecoul
		! entr    = srad + sion + sele + scoul
		
		! dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd
		! denerdd = deraddd + deiondd + deepdd + decouldd
		! dentrdd = dsraddd + dsiondd + dsepdd + dscouldd

		HelmState % prad = prad
		HelmState % pion = pion
		HelmState % pele = pele
		HelmState % pcoul = pcoul

		HelmState % erad = erad
		HelmState % eion = eion
		HelmState % eele = eele
		HelmState % ecoul = ecoul

		HelmState % srad = srad
		HelmState % sion = sion
		HelmState % sele = sele
		HelmState % scoul = scoul
		
		HelmState % dpraddd = dpraddd
		HelmState % dpiondd = dpiondd
		HelmState % dpepdd = dpepdd
		HelmState % dpcouldd = dpcouldd

		HelmState % deraddd = deraddd
		HelmState % deiondd = deiondd
		HelmState % deepdd = deepdd
		HelmState % decouldd = decouldd
		
		HelmState % dsraddd = dsraddd
		HelmState % dsiondd = dsiondd
		HelmState % dsepdd = dsepdd
		HelmState % dscouldd = dscouldd
		
	END SUBROUTINE HelmEOSComponents
	
	
	! These are all the functions used fro the interpolation
	! provided by Timmes
    ! quintic hermite polynomial functions
    ! psi0 and its derivatives
    REAL(dp) FUNCTION psi0(z)
		REAL(dp), INTENT(IN) :: z

		psi0 = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
	END FUNCTION psi0
	
    REAL(dp) FUNCTION dpsi0(z)
		REAL(dp), INTENT(IN) :: z

		dpsi0 = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
	END FUNCTION dpsi0
	
    REAL(dp) FUNCTION ddpsi0(z)
		REAL(dp), INTENT(IN) :: z

		ddpsi0 = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)
	END FUNCTION ddpsi0
	
    ! psi1 and its derivatives
    REAL(dp) FUNCTION psi1(z)
		REAL(dp), INTENT(IN) :: z

		psi1 = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
	END FUNCTION psi1
	
    REAL(dp) FUNCTION dpsi1(z)
		REAL(dp), INTENT(IN) :: z

		dpsi1 = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
	END FUNCTION dpsi1
	
    REAL(dp) FUNCTION ddpsi1(z)
		REAL(dp), INTENT(IN) :: z

		ddpsi1 = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
	END FUNCTION ddpsi1
	
    ! psi2  and its derivatives
    REAL(dp) FUNCTION psi2(z)
		REAL(dp), INTENT(IN) :: z

		psi2 = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
	END FUNCTION psi2
	
    REAL(dp) FUNCTION dpsi2(z)
		REAL(dp), INTENT(IN) :: z

		dpsi2 = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
	END FUNCTION dpsi2
	
    REAL(dp) FUNCTION ddpsi2(z)
		REAL(dp), INTENT(IN) :: z

		ddpsi2 = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)
	END FUNCTION ddpsi2
	
	
    ! biquintic hermite polynomial FUNCTION
    REAL(dp) FUNCTION h5(fi,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)
		REAL(dp), INTENT(IN) :: fi(36)
		REAL(dp), INTENT(IN) :: w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md
				
		h5 =  fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
		+ fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
		+ fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
		+ fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
		+ fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
		+ fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
		+ fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
		+ fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
		+ fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
		+ fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
		+ fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
		+ fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
		+ fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
		+ fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
		+ fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
		+ fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
		+ fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
		+ fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt

	END FUNCTION h5
	
    ! cubic hermite polynomial FUNCTIONs
    ! psi0 & derivatives
    REAL(dp) FUNCTION xpsi0(z)
		REAL(dp), INTENT(IN) :: z

		xpsi0 = z * z * (2.0d0*z - 3.0d0) + 1.0
	END FUNCTION xpsi0
	
    REAL(dp) FUNCTION xdpsi0(z)
		REAL(dp), INTENT(IN) :: z

		xdpsi0 = z * (6.0d0*z - 6.0d0)
	END FUNCTION xdpsi0
	
	
    ! psi1 & derivatives
    REAL(dp) FUNCTION xpsi1(z)
		REAL(dp), INTENT(IN) :: z

		xpsi1 = z * ( z * (z - 2.0d0) + 1.0d0)
	END FUNCTION xpsi1
	
    REAL(dp) FUNCTION xdpsi1(z)
		REAL(dp), INTENT(IN) :: z

		xdpsi1 = z * (3.0d0*z - 4.0d0) + 1.0d0
	END FUNCTION xdpsi1
	
    ! bicubic hermite polynomial FUNCTION
    REAL(dp) FUNCTION h3(dfi,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md)
		REAL(dp), INTENT(IN) :: dfi(16)
		REAL(dp), INTENT(IN) :: w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md

		h3 =  dfi(1)  *w0d*w0t   +  dfi(2)  *w0md*w0t &
		+ dfi(3)  *w0d*w0mt  +  dfi(4)  *w0md*w0mt &
		+ dfi(5)  *w0d*w1t   +  dfi(6)  *w0md*w1t &
		+ dfi(7)  *w0d*w1mt  +  dfi(8)  *w0md*w1mt &
		+ dfi(9)  *w1d*w0t   +  dfi(10) *w1md*w0t &
		+ dfi(11) *w1d*w0mt  +  dfi(12) *w1md*w0mt &
		+ dfi(13) *w1d*w1t   +  dfi(14) *w1md*w1t &
		+ dfi(15) *w1d*w1mt  +  dfi(16) *w1md*w1mt
	END FUNCTION h3
	
END MODULE wlElectronEOS
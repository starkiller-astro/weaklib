MODULE wlElectronPhotonEOS
  
  USE wlKindModule, ONLY: dp
  USE wlLeptonEOSModule, ONLY: &
    HelmholtzTableType
  USE wlEosConstantsModule, ONLY: &
    pi, rmu, kerg, cvel, me, &
    kmev, ergmev, sigma_sb, asol
    
  IMPLICIT NONE
  PRIVATE
  TYPE, PUBLIC :: ElectronPhotonStateType
    
    REAL(dp) :: T
    REAL(dp) :: rho
    REAL(dp) :: ye

    REAL(dp) :: p
    REAL(dp) :: pele
    REAL(dp) :: prad
    REAL(dp) :: dpdT
    REAL(dp) :: dpdr

    REAL(dp) :: e
    REAL(dp) :: eele
    REAL(dp) :: erad
    REAL(dp) :: dedT
    REAL(dp) :: dedr
    
    REAL(dp) :: s
    REAL(dp) :: sele
    REAL(dp) :: srad
    REAL(dp) :: dsdT
    REAL(dp) :: dsdr
  
    REAL(dp) :: mue
    
  END TYPE ElectronPhotonStateType

  PUBLIC :: ElectronPhotonEOS
  
CONTAINS

  SUBROUTINE ElectronPhotonEOS(HelmTable, ElectronPhotonState)
     
    !..input arguments
    TYPE(HelmholtzTableType), INTENT(IN) :: HelmTable
    TYPE (ElectronPhotonStateType), INTENT(INOUT) :: ElectronPhotonState

  !..declare local variables

  REAL(dp) :: x,deni,tempi,xni, &
  dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
  dpraddd,dpraddt,deraddd,deraddt, &
  deiondd,deiondt,dsraddd,dsraddt, &
  dse,dpe,dsp,kt,ktinv,prad,erad,srad, &
  pele,eele,sele,pres,ener,entr,dpresdd, &
  dpresdt,denerdd,denerdt,dentrdd,dentrdt,etaele, &
  detadt,detadd,xnefer,dxnedt,dxnedd,s, &
  temp,den,ye,din

  !..for the interpolations
  INTEGER  :: iat,jat
  REAL(dp) :: free,df_d,df_t,df_tt,df_dt
  REAL(dp) :: xt,xd,mxt,mxd,fi(36), &
  si0t,si1t,si2t,si0mt,si1mt,si2mt, &
  si0d,si1d,si2d,si0md,si1md,si2md, &
  dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
  dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
  ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt

  ! extra variables I added that were previously allocated
  REAL(dp) :: smallt, smalld, hight, highd
  REAL(dp) :: tstp, dstp, dstpi, tstpi

  smallt = HelmTable % mintemp
  smalld = HelmTable % mindens
  hight = HelmTable % maxtemp
  highd = HelmTable % maxdens
      
  tstp  = (LOG10(hight) - LOG10(smallt))/float(HelmTable % nPointsTemp -1)
  dstp  = (LOG10(highd) - LOG10(smalld))/float(HelmTable % nPointsDen -1)

  tstpi = 1.0_dp / tstp
  dstpi = 1.0_dp / dstp

  temp = ElectronPhotonState % T
  den  = ElectronPhotonState % rho
  ye   = ElectronPhotonState % ye
              
  din   = ye * den

  !..initialize
  deni    = 1.0d0/den
  tempi   = 1.0d0/temp
  kt      = kerg * temp
  ktinv   = 1.0d0/kt

  !..radiation section:
  prad    = asol/3.0d0 * temp * temp * temp * temp
  dpraddd = 0.0d0
  dpraddt = 4.0d0 * prad*tempi

  erad    = 3.0d0 * prad*deni
  deraddd = -erad*deni
  deraddt = 3.0d0 * dpraddt*deni

  srad    = (prad*deni + erad)*tempi
  dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
  dsraddt = (dpraddt*deni + deraddt - srad)*tempi

  !..ion section: removed in this
  
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
  xt  = MAX( (temp - HelmTable % t(jat))*HelmTable % dti(jat), 0.0_dp)
  xd  = MAX( (din - HelmTable % d(iat))*HelmTable % ddi(iat), 0.0_dp)
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
  dpepdd  = MAX(ye * dpepdd,0.0_dp)

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
  x = MAX(x,0.0_dp)
  dxnedd   = ye * x
  
  !..derivative with respect to temperature
  dxnedt   = h3( fi, &
  dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
  si0d,   si1d,   si0md,   si1md)

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

  pres    = prad + pele
  ener    = erad + eele
  entr    = srad + sele

  dpresdd = dpraddd + dpepdd
  dpresdt = dpraddt + dpepdt

  denerdd = deraddd + deepdd
  denerdt = deraddt + deepdt

  dentrdd = dsraddd + dsepdd
  dentrdt = dsraddt + dsepdt

  !..maxwell relations; each is 0.0_dp IF the consistency is perfect
  ! x   = den * den
  ! dse = temp*dentrdt/denerdt - 1.0d0
  ! dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
  ! dsp = -dentrdd*x/dpresdt - 1.0d0

  ElectronPhotonState % T    = temp
  ElectronPhotonState % rho  = den
  ElectronPhotonState % ye   = ye

  ElectronPhotonState % p    = pres
  ElectronPhotonState % pele = pele
  ElectronPhotonState % prad = prad
  ElectronPhotonState % dpdT = dpresdt
  ElectronPhotonState % dpdr = dpresdd

  ElectronPhotonState % e    = ener + me / rmu * ergmev * ye
  ElectronPhotonState % eele = eele + me / rmu * ergmev * ye
  ElectronPhotonState % erad = erad
  ElectronPhotonState % dedT = denerdt
  ElectronPhotonState % dedr = denerdd

  ElectronPhotonState % s    = entr    / (kmev * ergmev / rmu)
  ElectronPhotonState % sele = sele    / (kmev * ergmev / rmu)
  ElectronPhotonState % srad = srad    / (kmev * ergmev / rmu)
  ElectronPhotonState % dsdT = dentrdt / (kmev * ergmev / rmu)
  ElectronPhotonState % dsdr = dentrdd / (kmev * ergmev / rmu)

  ! do not forget to add back mass of the electron!
  ElectronPhotonState % mue = etaele*temp*kmev + me

END SUBROUTINE ElectronPhotonEOS

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
  
END MODULE wlElectronPhotonEOS

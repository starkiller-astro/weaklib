MODULE CC_module_weaklib

  USE wlKindModule, ONLY: dp
  
  USE HDF5
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  USE wlInterpolationUtilitiesModule, ONLY: &
    GetIndexAndDelta_Lin, &
    GetIndexAndDelta_Log, &
    LinearInterp_Array_Point
  
 
  USE, INTRINSIC :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  implicit none

  real(dp), parameter :: pi=3.1415927d0, &
                         Gf=1.166d-11, &
                         Vud=0.97427d0, &
                         gA0=1.2723d0, &
                         gV0=1.d0, &
                         f2wm0=3.706d0

   real(dp) :: Tem,mun,mup,mue,mn,mp,Un,Up,Enu,En, &
               pnmax,pnmin,ppmax,pemax,pemin

   real(dp) :: me

   real(dp), parameter :: Tfac=100.0d0,Tfac2=100.0d0

   real(dp), parameter :: Q = 1.2935d0

   integer:: anti,opt0

  CONTAINS



  SUBROUTINE Get_LS220_UnUp(D, T, Y, Un, Up)

    USE el_eos_module, ONLY: EPRESS, EU, ES
    USE eos_m4c_module, ONLY: PTOT, UTOT, STOT, XNUT, XPROT, XH, MUN, &
        MUPROT, A, X, GAM_S, BUNUC, VNOUT, VPOUT
    USE wlExtPhysicalConstantsModule, ONLY: dmnp, kmev, rmu, cm3fm3, ergmev, &
        asig
    USE e_p_eos_module

    implicit none

    real(dp), dimension(:), intent(in)        :: D, T, Y
    real(dp), dimension(:,:,:), intent(inout) :: Un, Up

    real(dp), parameter     :: kfm = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
    real(dp), parameter     :: kp  = ergmev/cm3fm3 ! ( erg/cm3 ) / ( mev/fm3 )
    real(dp), parameter     :: ku  = ergmev/rmu    ! ( # nucleons/gram )( erg/mev )
    real(dp), parameter     :: UTOT0 = 8.9d0       ! change in the zero of energy [MeV]
    real(dp), parameter     :: me = 0.510998d+00   ! electron mass [MeV

    !loop indices
    integer :: ii, jj, kk

    !original table size
    integer :: nD, nT, nY

    integer                 :: iflag       ! LS input type flag
    integer                 :: eosflg      ! LS output EoS type flag
    integer                 :: forflg      ! LS EoS forcing flag
    integer                 :: sf          ! eos success flag
    integer                 :: ieos        ! eos counter

    real(dp), dimension(4)  :: inpvar      ! input variables for LS EoS

    real(dp)                :: pe          ! raw electron pressure [dynes cm^{-2}]
    real(dp)                :: ee          ! raw electron energy [ergs cm^{-3}]
    real(dp)                :: entrop_e    ! entropy of electron/photon gas
    real(dp)                :: yeplus      ! positron fraction [unused]
    real(dp)                :: rel         ! relativity parameter [unused]
    real(dp)                :: tmev        ! temperature [MeV]
    real(dp)                :: brydns      ! LS input baryon density [fm^-3]

    real(dp)                :: xprev       ! LS input protron fraction
    real(dp)                :: pprev       ! previous value of proton density
    real(dp)                :: t_old       ! old value of temperature

    real(dp) :: ye_loc

    character(len=3)        :: LScompress
    character(len=128)      :: LSFilePath

    iflag     = 1
    forflg    = 0
    sf        = 1

    nD = size(D)
    nT = size(T)
    nY = size(Y)

    LScompress = '220'
    LSFilePath = '../../../Old/LS/Data/'

write(*,*) 5.2d-08 / kfm
write(*,*) 1d-12 / kfm


    CALL wlExtInitializeEOS( LSFilePath, LScompress )

    do ii = nD, nD
      do jj = nT, nT
        do kk = 45, 45

          brydns    = 11.937766417144367d0 !D(ii) * kfm

          tmev      = 181.97008586099830d0 !T(jj) * kmev
 
          t_old     = tmev

          ye_loc    = 3.0d-2
 
          pprev     = ye_loc * brydns !Y(kk) * brydns
  
          inpvar(1) = tmev
          inpvar(2) = 0.155d0
          inpvar(3) = -15.d0
          inpvar(4) = -10.d0

!write(*,*) ii, jj, kk
!write(*,*) D(ii), T(jj), Y(kk)
!write(*,*) brydns, tmev, pprev, xprev
!stop
          CALL inveos( inpvar, t_old, ye_loc, brydns, iflag, eosflg, forflg, sf,    &
               xprev, pprev )


          if ( sf == 0 ) then
write(*,*) 'Direct LS220 call failed at' ,ii ,jj ,kk   
Write(*,*) D(ii), T(jj), Y(kk)
          else
write(*,*) 'Direct LS220 call succeedes at' ,ii ,jj ,kk   
Write(*,*) D(ii), T(jj), Y(kk)
write(*,*) brydns, tmev
Write(*,*) VNOUT, VPOUT
            Un(ii,jj,kk) = VNOUT
            Up(ii,jj,kk) = VPOUT
          end if

         

        end do
      end do
    end do

write (*,*) 'Un', minval(Un), maxval(Un)
write (*,*) 'Up', minval(Up), maxval(Up)

  END SUBROUTINE Get_LS220_UnUp

  SUBROUTINE Interp_UnUp(filename, D, T, Y, Un, Up, massn, massp, which_EOS)

    character(len=*), intent(in)              :: filename 
    real(dp), dimension(:), intent(in)        :: D, T, Y
    real(dp), dimension(:,:,:), intent(inout) :: Un, Up
    real(dp), dimension(:,:,:), intent(inout) :: massn, massp
    integer, intent(in)                       :: which_EOS !1 for LS220, 2 for SFHo

    !interpolation indices and deltas
    integer :: iD, iT, iY
    real(dp) :: dD, dT, dY

    !loop indices
    integer :: ii, jj, kk

    !original table size
    integer :: nD, nT, nY

    !UnUp table size and data
    integer :: nDs, nTs, nYs
    real(dp), dimension(:), allocatable :: Ds, Ts, Ys
    real(dp), dimension(:,:,:), allocatable :: Uns, Ups
    real(dp), dimension(:,:,:), allocatable :: massns, massps

    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HSIZE_T), dimension(1) :: datasize1d    
    integer(HSIZE_T), dimension(3) :: datasize3d    

    integer, dimension(3) :: table_dimensions
   
    real(dp), dimension(4) :: offsets

    integer :: d_under, d_over, t_under, t_over

    real(dp), parameter :: mp =  938.272046d0  !proton restmass
    real(dp), parameter :: mn =  939.565379d0  !neutron restmass

    massn = mn
    massp = mp

    nD = size(D)
    nT = size(T)
    nY = size(Y)

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "ThermoState", .false., file_id, group_id )

    datasize1d = size(table_dimensions)
    if(which_EOS .eq. 1) then
      CALL ReadHDF("Dimensions_ls220", table_dimensions(:), group_id, datasize1d )
    else if (which_EOS .eq. 2) then
      CALL ReadHDF("Dimensions_sfho", table_dimensions(:), group_id, datasize1d )
    endif

    nDs = table_dimensions(1)
    nTs = table_dimensions(2)
    nYs = table_dimensions(3) 

    ALLOCATE(Ds (nDs))
    ALLOCATE(Ts (nTs))
    ALLOCATE(Ys (nYs))
    ALLOCATE(Uns(nDs,nTs,nYs))
    ALLOCATE(Ups(nDs,nTs,nYs))

    if(which_EOS .eq. 2) then

      ALLOCATE(massns(nDs,nTs,nYs))
      ALLOCATE(massps(nDs,nTs,nYs))

    endif

    datasize1d = nDs
    if(which_EOS .eq. 1) then
      CALL ReadHDF("Density_ls220", Ds(:), group_id, datasize1d)
    else if (which_EOS .eq. 2) then
      CALL ReadHDF("Density_sfho", Ds(:), group_id, datasize1d)
    endif

    datasize1d = nTs
    if(which_EOS .eq. 1) then
      CALL ReadHDF("Temperature_ls220", Ts(:), group_id, datasize1d)
    else if (which_EOS .eq. 2) then
      CALL ReadHDF("Temperature_sfho", Ts(:), group_id, datasize1d)
    endif
    
    datasize1d = nYs
    if(which_EOS .eq. 1) then
      CALL ReadHDF("Electron Fraction_ls220", Ys(:), group_id, datasize1d)
    else if (which_EOS .eq. 2) then
      CALL ReadHDF("Electron Fraction_sfho", Ys(:), group_id, datasize1d)
    endif

    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "DependentVariables", .false., file_id, group_id )

    datasize1D = size(offsets)
    if(which_EOS .eq. 1) then
      CALL ReadHDF("Offsets_ls220", offsets(:), group_id, datasize1d )
    else if (which_EOS .eq. 2) then
      CALL ReadHDF("Offsets_sfho", offsets(:), group_id, datasize1d )
    endif 

    datasize3d = [nDs,nTs,nYs]
    if(which_EOS .eq. 1) then
      CALL READHDF("Un_ls220", Uns(:,:,:), group_id, datasize3d)
      CALL READHDF("Up_ls220", Ups(:,:,:), group_id, datasize3d)
    else if (which_EOS .eq. 2) then
      CALL READHDF("Un_sfho", Uns(:,:,:), group_id, datasize3d)
      CALL READHDF("Up_sfho", Ups(:,:,:), group_id, datasize3d)
      CALL READHDF("massn_sfho", massns(:,:,:), group_id, datasize3d)
      CALL READHDF("massp_sfho", massps(:,:,:), group_id, datasize3d)
    endif

    CALL CloseGroupHDF( group_id )

    CALL CloseFileHDF( file_id )

write(*,*) 'min/max D', minval(D), maxval(D)
write(*,*) 'min/max Ds', minval(Ds), maxval(Ds)

write(*,*) 'min/max T', minval(T), maxval(T)
write(*,*) 'min/max Ts', minval(Ts), maxval(Ts)

write(*,*) 'min/max Y', minval(Y), maxval(Y)
write(*,*) 'min/max Ys', minval(Ys), maxval(Ys)

write (*,*) 'min/max Uns', minval(10**Uns)-offsets(1), maxval(10**Uns)-offsets(1)
write (*,*) 'min/max Ups', minval(10**Ups)-offsets(2), maxval(10**Ups)-offsets(2)

write (*,*) 'min/max Un', minval(Un), maxval(Un)
write (*,*) 'min/max Up', minval(Up), maxval(Up)
write (*,*) 'min/max massn', minval(massn), maxval(massn)
write (*,*) 'min/max massp', minval(massp), maxval(massp)

write(*,*) 'offsets', offsets(1), offsets(2)
if(which_EOS .eq. 2) then
  write(*,*) 'offsets', offsets(3), offsets(4)
endif

!write(*,*) nDs, nTs, nYs
!write(*,*) nD, nT, nY

d_under = 0
d_over  = 0
t_under = 0
t_over  = 0

    do ii = 1, nD
      do jj = 1, nT
        do kk = 1, nY
 
          if(D(ii) .le. minval(Ds)) then
            iD = 1
            dD = LOG10( D(ii) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )     
            d_under = d_under + 1
          else if (D(ii) .ge. maxval(Ds)) then
            iD = nDs-1
            dD = LOG10( D(ii) / Ds(iD) ) / LOG10( Ds(iD+1) / Ds(iD) )     
            d_over = d_over + 1
          else
            CALL GetIndexAndDelta_Log( D(ii), Ds, iD, dD )
          endif

          if(T(jj) .le. minval(Ts)) then
            iT = 1
            dT = LOG10( T(jj) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )     
            t_under = t_under + 1
          else if (T(jj) .ge. maxval(Ts)) then
            iT = nTs-1
            dT = LOG10( T(jj) / Ts(iT) ) / LOG10( Ts(iT+1) / Ts(iT) )     
            t_over = t_over + 1
          else
            CALL GetIndexAndDelta_Log( T(jj), Ts, iT, dT )
          endif

          if(Y(kk) .le. minval(Ys)) then
            iY = 1
            dY = ( Y(kk) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )
          else if (Y(kk) .ge. maxval(Ys)) then
            iY = nYs-1
            dY = ( Y(kk) - Ys(iY) ) / ( Ys(iY+1) - Ys(iY) )
          else
            CALL GetIndexAndDelta_Lin( Y(kk), Ys, iY, dY )
          endif

          CALL LinearInterp_Array_Point &
               (iD, iT, iY, dD, dT, dY, offsets(1), Uns, Un(ii,jj,kk))
  
          CALL LinearInterp_Array_Point &
               (iD, iT, iY, dD, dT, dY, offsets(2), Ups, Up(ii,jj,kk))

          if(which_EOS .eq. 2) then
            CALL LinearInterp_Array_Point &
                 (iD, iT, iY, dD, dT, dY, offsets(3), massns, massn(ii,jj,kk))
  
            CALL LinearInterp_Array_Point &
                 (iD, iT, iY, dD, dT, dY, offsets(4), massps, massp(ii,jj,kk))

          endif

        end do
      end do
    end do

write (*,*) 'Un interp min/max', minval(Un), maxval(Un)
write (*,*) 'Up interp min/max', minval(Up), maxval(Up)
write (*,*) 'Un interp min/max loc', minloc(Un), maxloc(Un)
write (*,*) 'Up interp min/max loc', minloc(Up), maxloc(Up)

write (*,*) 'massn interp min/max', minval(massn), maxval(massn)
write (*,*) 'massp interp min/max', minval(massp), maxval(massp)

!stop
  

    DEALLOCATE(Ds)
    DEALLOCATE(Ts)
    DEALLOCATE(Ys)

    DEALLOCATE(Uns)
    DEALLOCATE(Ups)

    if(which_EOS .eq. 2) then

      DEALLOCATE(massns)
      DEALLOCATE(massps)

    endif

  END SUBROUTINE Interp_UnUp

  SUBROUTINE CC_EmAb(e_nu, xTem, m, cheme, chemn, chemp, &
                     massn, massp, xUn, xUp,             &
                     species, inv_n_decay, absor, emit,  &
                     ab_inv_n_decay, em_inv_n_decay, nPointsE)

  real(dp), intent(in) :: xTem, m, cheme, chemn, chemp, &
                          massn, massp, xUn, xUp

  integer, intent(in)  :: species, inv_n_decay, nPointsE

  !process=1: v + n -> e- + p
  !process=2: vb + p -> e+ + n
  !process=3: vb + p + e- -> n

  real(dp), dimension(nPointsE), intent(in) :: e_nu ! neutrino energy grid [MeV]

  real(dp), dimension(nPointsE), intent(inout) :: absor, emit
  real(dp), dimension(nPointsE), intent(inout) :: ab_inv_n_decay
  real(dp), dimension(nPointsE), intent(inout) :: em_inv_n_decay

  integer :: k
  real(dp) :: op3, chemhat

  real(dp) :: fexp                 ! exponential

  external fexp

  chemhat = chemn-chemp

  me   = m
  Tem  = xTem

  if(species .eq. 1) then

    mue  = cheme
    mn   = massn
    mp   = massp
    Un   = xUn
    Up   = xUp
    mun  = chemn
    mup  = chemp
    anti = 1
    opt0 = 1
write(*,*) 'hello there', Tem, species, inv_n_decay
    do k = 1, nPointsE
      call Integral_2D(e_nu(k),absor(k))
      absor(k) = absor(k) / 1.0d2
      write(*,*) e_nu(k), absor(k)
      !detailed balance
      emit(k) = absor(k) / (fexp((e_nu(k)-cheme+chemhat+Q)/Tem))
    end do

  endif

  if(species .eq. 2) then

    if(inv_n_decay .eq. 1) then
      mp   = massn
      mn   = massp
      Up   = xUn
      Un   = xUp
      mup  = chemn
      mun  = chemp
      anti = -1
      opt0 = 3
      do k = 1,nPointsE
        call Integral_2D_D(e_nu(k),ab_inv_n_decay(k))
        ab_inv_n_decay(k) = ab_inv_n_decay(k) / 1.0d2
        !detailed balance
        em_inv_n_decay(k) = ab_inv_n_decay(k) / &
                            (fexp((e_nu(k)+cheme-chemhat-Q)/Tem))
      enddo

    else
      mue = -cheme
      opt0 = 2
      do k = 1, nPointsE
        call Integral_2D(e_nu(k),absor(k))
        absor(k) = absor(k) / 1.0d2
        !detailed balance
        emit(k)  = absor(k) / &
                            (fexp((e_nu(k)+cheme-chemhat-Q)/Tem))
      enddo

    endif
  endif

  END SUBROUTINE CC_EmAb


  SUBROUTINE Integral_2D(xEnu,res)   ! captures

     !use CC_global

    implicit none

!.....input.............................................................
    real(dp), intent(in) :: xEnu

!.....output............................................................
    real(dp), intent(out) :: res

!.....locals............................................................
    integer :: i, j
    real(dp) :: xamp,Ep
    real(dp) :: xf2,xf3,xf4
    real(dp) :: xpn3, xq0_min, xq0_max, Ee_min,Ee_max,xEe, En_min, En_max

    integer, parameter :: N1=64, N2=64

    real(dp), dimension(N1) :: Ena,wEn
    real(dp), dimension(N2) :: Eea,wEe

    real(dp) :: fexp                 ! exponential

    external fexp

    Enu = xEnu
    res = 0.0d0
    call Range_pn()
    !write(*,*) '2D int pnmin pnmax', pnmin, pnmax
    if ( pnmin > pnmax ) goto 12

    En_min = sqrt(pnmin**2+mn**2) + Un
    En_max = sqrt(pnmax**2+mn**2) + Un

    call gaulegg( En_min,En_max,Ena,wEn,N1 )

    do i=1,N1
      En = Ena(i)
      xpn3 = sqrt( (En-Un)**2 - mn**2 )
      call Ebounds( xpn3, xq0_min, xq0_max )
      Ee_min = Enu - xq0_max
      Ee_max = Enu - xq0_min
      call gaulegg( Ee_min,Ee_max,Eea,wEe,N2 )
      do j=1,N2
        xEe = Eea(j)
!write(*,*) j, En, xEe
        xf3 = 1.0d0/( fexp((xEe-mue)/Tem) + 1.0d0 )
        Ep = Enu+En-xEe
        call Calc_ampsq(En,xEe,xamp)
        if (xamp.lt.0.0d0) xamp = 0.0d0
        xf2 = 1.0d0/( fexp((En-mun)/Tem) + 1.0d0 )
        xf4 = 1.0d0/( fexp((Ep-mup)/Tem) + 1.0d0 )
        res = res + xf2*(1.0d0-xf3)*(1.0d0-xf4)*xamp*wEn(i)*wEe(j)/Enu**2
      end do
    end do

!      write(*,*) 'res=',res

    res = res*(gf*Vud)**2/16.0d0/(pi**5)/197.327d-18*1.d-3   ! m^-1

12    continue

  END SUBROUTINE Integral_2D

  SUBROUTINE Integral_2D_D(xEnu,res)   ! decay/inverse decay
     !use CC_global
    implicit none

    real(dp) :: xpn3, xq0_min, xq0_max, Ee_min,Ee_max,xEe, En_min, En_max

    integer,parameter :: N1=64, N2=64

    real(dp) :: Ena(N1),wEn(N1),Eea(N2),wEe(N2)
    integer  :: i, j
    real(dp) :: xEnu, res, xamp
    real(dp) :: E1,E2,E3,E4,E2f,E4f,mu2,mu4,mu3,U2,U4,Ep
    real(dp) :: xf2,xf3,xf4

    real(dp) :: fexp                 ! exponential

    external fexp

    U2 = Un
    U4 = Up
    mu2 = mun
    mu4 = mup
    mu3 = mue

    Enu = xEnu
    E1 = xEnu
    res = 0.0d0
    call Range_pn_D()
    if( pnmin > pnmax ) goto 12
    En_min = sqrt(pnmin**2+mn**2) + Un
    En_max = sqrt(pnmax**2+mn**2) + Un

    call gaulegg( En_min,En_max,Ena,wEn,N1 )

    do i=1,N1
      En = Ena(i)
      xpn3 = sqrt( (En-Un)**2-mn**2 )
      call Ebounds_D( xpn3, xq0_min, xq0_max )
      Ee_min = xq0_min - Enu
      Ee_max = xq0_max - Enu
      call gaulegg( Ee_min,Ee_max,Eea,wEe,N2 )
      do j=1,N2
        xEe = Eea(j)
        Ep = Enu+xEe+En
        call Calc_Ampsq_D(En,xEe,xamp)
        xf2 = 1.0d0/( fexp((En-mun)/Tem) + 1.0d0 )
        xf3 = 1.0d0/( fexp((xEe-mue)/Tem) + 1.0d0 )
        xf4 = 1.0d0/( fexp((Ep-mup)/Tem) + 1.0d0 )
        if(opt0.eq.3) then          ! inverse decay
          res = res + xf2*xf3*(1.0-xf4)*xamp*wEn(i)*wEe(j)/Enu**2
        else if(opt0.eq.4) then     ! decay, *blocking of neutrinos is not added*
          res = res + (1.0d0-xf2)*(1.0d0-xf3)*xf4*xamp*wEn(i)*wEe(j)/Enu**2
        end if
      end do
     end do

    if(opt0.eq.3) then ! inverse decay, opacity, m^-1
      res = res*(gf*Vud)**2/16.0d0/(pi**5)/197.327d-15
    else if(opt0.eq.4) then  ! neutron decay, emmisivity, 1/(s cm^3 MeV)
      res = res*(gf*Vud)**2/32.0d0/(pi**7)*Enu**2*(1.d13/197.327d0)**3 &
          * (3.d23/197.327d0)
    end if

12    continue

  END SUBROUTINE Integral_2D_D

! subroutines for v/vb captures

  SUBROUTINE Calc_ampsq( E2, E3, ampsq )
         !use CC_global
    implicit none

    real(dp) :: Ampsq, P1,P2,P3,P4,Pmax,Pmin,Pmax2,Pmin2,Pmax3,Pmin3
    real(dp) :: Ia,Ib,Ic,Id,Ie,Ig,Ih,Ij,Ik,Il,A,B,Ct,D,E,G,H,J,K,L, &
                del0,del1,del2,eps0,eps1,alp0,alp1,alp2,bet0,bet1, &
                cplA,cplB,cplC,cplD,cplE,cplG,cplH,cplJ,cplK,cplL,JFF,KFF,LFF
    real(dp) :: m2f,E2f,m4f,E4f,E1,E3,Qmass,dmf,U2,U4,dU2,m3,E2,E4,m2,xE3
    real(dp) :: teps0, teps1
    real(dp) :: Iam,Iap,Ibm,Ibp,Icm,Icp,Idm,Idp,Iem,Iep,Igm,Igp,Ikm,Ikp
    real(dp) :: gasq=0.0d0,gvsq=0.0d0,gva=0.0d0,gaf=0.0d0,gvf=0.0d0,f2sq=0.0d0

    m2 = 0.5d0*(939.565d0+938.272d0)   ! bare mass for the weak magnetism term in hardonic current
    m2f = mn                     ! actually particle-2
    m4f = mp                     ! actually particle-4
    m3 = me
    U2 = Un
    U4 = Up
    dU2 = U2 - U4
    dmf = m2f - m4f
    Qmass = 0.5d0*( m2f**2 - m4f**2 )

    E1 = Enu
    E4 = E1 + E2 - E3
    E4f = E4 - U4
    E2f = E2 - U2
    P1 = E1
    if(E3**2 .ge. m3**2) then
      P3 = sqrt( E3**2 - m3**2 )
    else
      P3 = 0.0d0
    endif
    if((E2 - U2)**2 .ge. m2f**2) then
      P2 = sqrt( (E2 - U2)**2 - m2f**2 )
    else
      P2 = 0.0d0
    endif
    if((E4 - U4)**2 .ge. m4f**2) then 
      P4 = sqrt( (E4 - U4)**2 - m4f**2 )
    else
      P4 = 0.0d0
    endif

    if( P1 .ne. P1) write(*,*) 'nan in hh1_a'
    if( P2 .ne. P2) write(*,*) 'nan in hh2_a'
    if( P3 .ne. P3) then
   ! debug. TF
   !        write(*,*) 'hh3_a'
      P3 = 0.
   !        write(*,*) 'E=',E1,E2,E3,E4
   !        write(*,*) 'U=',U2,U4
   !        write(*,*) 'mf=',m2f,m4f
   !        write(*,*) 'P=',P1,P2,P3,P4
   !        write(*,*)
   !        pause 'hh3_a'
    endif
    if( P4 .ne. P4) then
   ! debug. TF
   !        write(*,*) 'P4=',P4
      P4 = 0.  !P1+P2-P3
   !        write(*,*) 'E=',E1,E2,E3,E4
   !        write(*,*) 'U=',U2,U4
   !        write(*,*) 'mf=',m2f,m4f
   !        write(*,*) 'P=',P1,P2,P3,P4
   !        write(*,*)
   !        pause 'hh4_a'
    endif

    if ((P1+P2)>(P3+P4)) then
      Pmax = P3+P4
    else
      Pmax = P1+P2
    endif
    if ((abs(P1-P2))>(abs(P3-P4))) then
      Pmin = abs(P1-P2)
    else
      Pmin = abs(P3-P4)
    endif

    A = E1*E2f+0.5d0*(P1**2+P2**2)
    B = E3*E4f+0.5d0*(P3**2+P4**2)
    Ia = pi**2/15.0d0*(3.0d0*(Pmax**5-Pmin**5) &
       - 10.0d0*(A+B)*(Pmax**3-Pmin**3)+60.0d0*A*B*(Pmax-Pmin))

    if ((P1+P4)>(P3+P2)) then
      Pmax2 = P3+P2
    else
      Pmax2 = P1+P4
    endif
    if ((abs(P1-P4))>(abs(P3-P2))) then
      Pmin2 = abs(P1-P4)
    else
      Pmin2 = abs(P3-P2)
    endif

    Ct = E1*E4f-0.5d0*(P1**2+P4**2)
    D = E2f*E3-0.5d0*(P2**2+P3**2)
    Ib = pi**2/15.0d0*(3.0d0*(Pmax2**5-Pmin2**5) &
       + 10.0d0*(Ct+D)*(Pmax2**3-Pmin2**3)+60.0d0*Ct*D*(Pmax2-Pmin2))

    del0 = -(P1**2-P2**2)*(P3**2-P4**2)/4.0d0
    del1 = E1*E3 + (-P1**2+P2**2-P3**2+P4**2)/4.0d0
    del2 = -0.25d0
    eps0 = E1*E2f + (P1**2+P2**2)/2.0d0
    eps1 = -0.5d0
    Ic =  pi**2*4.0d0*(del2*eps1**2/7.0d0*(Pmax**7-Pmin**7) &
       + (2.0d0*del2*eps0*eps1+del1*eps1**2)/5.0d0*(Pmax**5-Pmin**5) &
       + (del2*eps0**2+2.0d0*del1*eps0*eps1+del0*eps1**2)/3.0d0 &
       * (Pmax**3-Pmin**3) &
       + (del1*eps0**2+2.0d0*del0*eps0*eps1)*(Pmax-Pmin) &
       - del0*eps0**2*(1.0d0/Pmax-1.0d0/Pmin))

    if ((P1+P3)>(P2+P4)) then
      Pmax3 = P2+P4
    else
      Pmax3 = P1+P3
    endif
    if ((abs(P1-P3))>(abs(P2-P4))) then
      Pmin3 = abs(P1-P3)
    else
      Pmin3 = abs(P2-P4)
    endif

    alp0 = (P1**2-P3**2)*(P2**2-P4**2)/4.0d0
    alp1 = E1*E2f + (P1**2+P2**2-P3**2-P4**2)/4.0d0
    alp2 = 0.25d0
    bet0 = E1*E3 - (P1**2+P3**2)/2.0d0
    bet1 = 0.5d0
    Id = pi**2*4.0d0*(alp2*bet1**2/7.0d0*(Pmax3**7-Pmin3**7) &
       + (2.0d0*alp2*bet0*bet1+alp1*bet1**2)/5.0d0*(Pmax3**5-Pmin3**5) &
       + (alp2*bet0**2+2.0*alp1*bet0*bet1+alp0*bet1**2)/3.0d0 &
       * (Pmax3**3-Pmin3**3) &
       + (alp1*bet0**2+2.0d0*alp0*bet0*bet1)*(Pmax3-Pmin3) &
       - alp0*bet0**2*(1.0d0/Pmax3-1.0d0/Pmin3))

    Ie = pi**2/15.0d0*(3.0d0*(Pmax**5-Pmin**5) &
       - 20.0d0*A*(Pmax**3-Pmin**3)+60.0d0*A**2*(Pmax-Pmin))

    E = E1*E3 - 0.5d0*(P1**2+P3**2)
    Ig = pi**2/15.0d0*(3.0d0*(Pmax3**5-Pmin3**5) &
       + 20.0d0*E*(Pmax3**3-Pmin3**3)+60.0d0*E**2*(Pmax3-Pmin3))

    Ih = pi**2*4.0d0*(alp2*bet1/5.0d0*(Pmax3**5-Pmin3**5) &
       + (alp2*bet0+alp1*bet1)/3.0d0*(Pmax3**3-Pmin3**3) &
       + (alp1*bet0+alp0*bet1)*(Pmax3-Pmin3) &
       - alp0*bet0*(1.0d0/Pmax3-1.0d0/Pmin3))

    Ij = pi**2/15.0d0*(-10.0d0*(Pmax**3-Pmin**3)+60.0d0*A*(Pmax-Pmin))

    Ik = pi**2/15.0d0*(10.0d0*(Pmax3**3-Pmin3**3)+60.0d0*E*(Pmax3-Pmin3))

    Il = pi**2/15.0d0*(60.0d0*(Pmax-Pmin))

    gasq = ga0**2
    gvsq = gv0**2
    gva = gv0*ga0
    gvf = gv0*F2wm0
    gaf = ga0*F2wm0
    f2sq = f2Wm0**2

    xE3 = E3
    cplA = (gvsq+2.*anti*gva+gasq) &
         + anti*2.0d0*gaf*m2f/m2*(1.0d0-dmf/2.0d0/m2f)
    cplB = (gvsq-2.*anti*gva+gasq) &
         - anti*2.0d0*gaf*m2f/m2*(1.0d0-dmf/2.0d0/m2f)
    cplC = F2sq/m2**2
    cplD =-F2sq/m2**2
    cplE = F2sq/m2**2*(-0.5d0*m3**2+dU2*(xE3-E1)-0.5d0*dU2**2)
    cplG = gvf*m2f/m2*(2.0d0-dmf/m2f) &
         + 0.5d0*F2sq/m2**2*(m2f*m4f-Qmass+0.25d0*m3**2 &
         - dU2*(E1+E2f)-dU2**2/4.0d0)
    cplH = 0.5d0*F2sq/m2**2*(2.0d0*Qmass+m3**2+dU2*(3.0d0*E1-xE3+2.0d0*E4f))

    JFF = -m3**2*(E1+0.5d0*E2f+0.5d0*E4f)+Qmass*(xE3-3.0d0*E1) &
        + 0.5d0*dU2*(E4f*(3.0d0*xE3-5.0d0*E1) &
        + E2f*(xE3+E1)+xE3**2-E1**2-2.0d0*Qmass) &
        + dU2**2*(E1-xE3-E4f)+0.5d0*dU2**3
    JFF = JFF*dU2
    cplJ = gvf*dmf/m2*0.5d0*(m3**2-dU2*(E1+xE3))+F2sq/m2**2*0.5d0*JFF

    kFF = -(m2f+3.*m4f)*m2f*0.25d0*m3**2 + Qmass**2 &
        + Qmass*0.25d0*m3**2-m3**4/8.0d0 &
        + dU2*(0.5d0*Qmass*(3.*E1-E2f+xE3+3.0d0*E4f) &
        + 0.25d0*m3**2*(2.0d0*E2f+xE3+E1)+m2f*m4f*(xE3-E1)) &
        + dU2**2*( 0.25d0*(m2f**2-m2f*m4f-3.0d0*Qmass) &
        + E4f*0.5d0*(2.0d0*E1-E2f+xE3+E4f) &
        + E2f*(0.5d0*E1-xE3)+0.5d0*E1**2 ) &
        + dU2**3*0.25d0*( -E1+2.0d0*E2f-xE3-2.0d0*E4f ) + 0.125d0*dU2**4
    cplK = (gasq-gvsq)*m2f*m4f &
         + gvf*m2f/m2*0.5d0*(-3.0d0*m3**2+4.0d0*dU2*(xE3-E1)-dU2**2 &
         + dmf/m2f*(2.0d0*Qmass+m3**2+dU2*(E4f+2.0d0*E1-xE3))) &
         + 0.5d0*F2sq/m2**2*KFF

    LFF = m3**2*(m2f+m4f)**2 &
        - 4.0d0*Qmass**2 + dU2*(-m3**2*(E2f+E4f+E1) &
        - 2.0d0*xE3*(m2f**2+m2f*m4f) &
        + 2.0d0*Qmass*(E2f-3.*E4f-E1)) &
        + 2.0d0*dU2**2*( E2f*xE3+E4f*(E2f-E4f-E1) + Qmass ) &
        + dU2**3*(-E2f+E4f+E1)
    LFF = LFF*dU2*E1*0.25d0
    cplL = gvf*m2f/m2*dU2*E1*(m3**2-dU2*xE3-0.5d0*dmf/m2f &
         * (Qmass+0.5d0*m3**2+dU2*E4f-0.5d0*dU2**2)) &
         + 0.5d0*F2sq/m2**2*LFF

    Ampsq = cplA*Ia &
          + cplB*Ib &
          + cplC*Ic &
          + cplD*Id &
          + cplE*Ie &
          + cplG*Ig &
          + cplH*Ih &
          + cplJ*Ij &
          + cplK*Ik &
          + cplL*Il
    return
  END SUBROUTINE Calc_ampsq


  SUBROUTINE Ebounds( xpn3, xq0_min, xq0_max )

    !use CC_global
    implicit none

    real(dp) :: xpn3, q3_lim(4), tmp1,tmp2,xA,xB,xC, xq3
    integer:: i
    real(dp) :: tmp3,tmp4,xq0_min,xq0_max,xq0, xq0_l,xq0_h

    En = sqrt( xpn3**2 + mn**2 ) + Un

    q3_lim = 0.0d0
    tmp1 = Up-En-Enu
    tmp2 = Enu**2 + me**2 - mp**2 - xpn3**2

    xA = 4.d0*( tmp1**2-(Enu-xpn3)**2 )
    xB = 4.d0*( (tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*xpn3 )
    xC = 4.d0*tmp1**2*(xpn3**2+mp**2) - (tmp2-tmp1**2)**2
    if( xB**2 - 4.0d0*xA*xC .ge. 0.0) then
      q3_lim(1) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
      q3_lim(2) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
    end if

    xA = 4.d0*( tmp1**2-(Enu+xpn3)**2 )
    xB = 4.d0*( -(tmp2-tmp1**2)*Enu - (tmp2+tmp1**2)*xpn3 )
    if( xB**2 - 4.0*xA*xC .ge. 0.0d0) then
      q3_lim(3) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
      q3_lim(4) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
    end if

    xq0_min = 1.d30
    xq0_max = -1.d30
    xq0 = 0.0d0
    do i=1,4
      xq3 = q3_lim(i)
      tmp1 = Enu - sqrt( (Enu-xq3)**2 + me**2 )
      tmp2 = Enu - sqrt( (Enu+xq3)**2 + me**2 )
      tmp3 = sqrt( (xpn3+xq3)**2+mp**2) + Up - En
      tmp4 = sqrt( (xpn3-xq3)**2+mp**2) + Up - En

      if( (abs((tmp1-tmp3)/tmp3)<1.d-6).or. &
          (abs((tmp1-tmp4)/tmp4)<1.d-6) ) then
        xq0 = tmp1
      end if
      if( (abs((tmp2-tmp3)/tmp3)<1.d-6).or. &
          (abs((tmp2-tmp4)/tmp4)<1.d-6) ) then
        xq0 = tmp2
      end if
      xq0_min = min( xq0,  xq0_min)
      xq0_max = max( xq0,  xq0_max)
    end do

    xq0_l = Enu - me
    tmp1 = sqrt( (xpn3+Enu)**2+mp**2) + Up - En
    tmp2 = sqrt( (xpn3-Enu)**2+mp**2) + Up - En
    if( xq0_l .le. tmp1 .and. xq0_l .ge. tmp2) then
      xq0_max = max( xq0_l, xq0_max)
    end if

    xq0_h = mp + Up - En
    tmp1 = Enu - sqrt( (Enu-xpn3)**2 + me**2 )
    tmp2 = Enu - sqrt( (Enu+xpn3)**2 + me**2 )
    if( xq0_h .le. tmp1 .and. xq0_h .ge. tmp2) then
      xq0_min = min( xq0_h, xq0_min)
    end if

   END SUBROUTINE Ebounds

   SUBROUTINE Range_pn()
     !use CC_global
     implicit none

     real(dp) :: tmp,tmp1,tmp2,tmp3,del,At,xpn3,mpt,tmp4
     real(dp) :: F0,Finf,Fmax,Fmin,pnmax0
     !real(dp),external :: test,test1,test2
     integer  :: i

     pnmin = 0.d0
     if((Tfac*Tem + mun - Un)**2 .ge. mn**2) then
       pnmax = sqrt( (Tfac*Tem + mun - Un)**2 - mn**2 )
     else 
       pnmax = 0.0d0
     endif
     pnmax0 = pnmax

     mpt = mp + me
     F0 = sqrt(Enu**2 + mpt**2) + Up - mn - Un
     Finf = -Enu + Up - Un
     if(mpt < mn) then
       tmp = (1.d0+mpt/mn)/(1.d0-mpt**2/mn**2)*Enu
       Fmin = sqrt( (tmp-Enu)**2+mpt**2 )+Up-sqrt(tmp**2+mn**2)-Un
       Fmax = max( F0, Finf )
     else
       Fmin = Finf
       Fmax = F0
     end if

     if( Fmax .le. Enu) goto 12
     if( Fmin .ge. Enu) then
       pnmax = -1000.0d0
       goto 12
     end if

     tmp = Enu - Up + Un
     At = Enu**2 + mpt**2 - mn**2 - tmp**2
     tmp1 = At**2 + 4.d0*( (mn*Enu)**2 - (tmp*mn)**2 )
     del = 1.d-5
     if(tmp1 > 0.d0 .and. Enu**2 .ne. tmp**2) then
       tmp2 = ( Enu*At - abs(tmp)*sqrt(tmp1) )*0.5d0/(Enu**2-tmp**2)
       tmp3 = ( Enu*At + abs(tmp)*sqrt(tmp1) )*0.5d0/(Enu**2-tmp**2)

       if(tmp2>0.d0) then
         if( (test1(tmp2-del)>0.d0 .and. test1(tmp2+del)<0.d0) &
            .or. (test1(tmp2-1.d0)>0.d0 .and. test1(tmp2+1.d0)<0.d0) )then
           pnmin = tmp2
         else if( (test1(tmp2-del)<0.d0 .and. test1(tmp2+del)>0.d0) &
                 .or. (test1(tmp2-1.d0)<0.d0 .and. test1(tmp2+1.d0)>0.d0) )then
           pnmax = tmp2
         end if
       end if

       if(tmp3 > 0.d0) then
         if( (test1(tmp3-del)>0.d0 .and. test1(tmp3+del)<0.d0) &
            .or. (test1(tmp3-1.d0)>0.d0 .and. test1(tmp3+1.d0)<0.d0) )then
           pnmin = tmp3
         else if( (test1(tmp3-del)<0.d0 .and. test1(tmp3+del)>0.d0) &
                 .or. (test1(tmp3-1.d0)<0.d0 .and. test1(tmp3+1.d0)>0.d0))then
           pnmax = tmp3
         end if
       end if

     else if(Enu**2 .eq. tmp**2) then
       tmp4 = -(4.0*tmp**2*mn**2-At**2)/(4.0*At*Enu)
       if(tmp4 > 0.0) then
         if( (test1(tmp4-del)>0.d0 .and. test1(tmp4+del)<0.d0) &
            .or. (test1(tmp4-1.d0)>0.d0 .and. test1(tmp4+1.d0)<0.d0) ) then
           pnmin = tmp4
         else if( (test1(tmp4-del)<0.d0 .and. test1(tmp4+del)>0.d0) &
                 .or. (test1(tmp4-1.d0)<0.d0 .and. test1(tmp4+1.d0)>0.d0))then
           pnmax = tmp4
         end if
       end if
     end if

12    continue

  END SUBROUTINE Range_pn

  FUNCTION test1(xPn3)
  !use CC_global
    implicit none
    real(dp) :: xpn3,test1,mpt

    mpt = mp + me
    test1 = sqrt((xpn3-Enu)**2 + mpt**2)-sqrt(xpn3**2+mn**2)+Up-Un-Enu
    return
  END FUNCTION test1

! subroutines for neutron decay and the inverse
  SUBROUTINE Calc_Ampsq_D( E2, E3, ampsq )
        !use CC_global
    implicit none
    real(dp) :: Ampsq, P1,P2,P3,P4,Pmax,Pmin,Pmax2,Pmin2,Pmax3,Pmin3
    real(dp) :: Ia,Ib,Ic,Id,Ie,Ig,Ih,Ij,Ik,Il,A,B,Ct,D,E,G,H,J,K,L, &
                del0,del1,del2,eps0,eps1,alp0,alp1,alp2,bet0,bet1, &
                cplA,cplB,cplC,cplD,cplE,cplG,cplH,cplJ,cplK,cplL,JFF,KFF,LFF
    real(dp) :: m2f,E2f,m4f,E4f,E1,E3,Qmass,dmf,U2,U4,dU2,m3,E2,E4,m2,xE3

    real(dp) :: gasq=0.d0,gvsq=0.d0,gva=0.d0,gaf=0.d0,gvf=0.d0,f2sq=0.d0

    m2 = 0.5d0*(939.565d0+938.272d0)   ! bare mass for the weak magnetism term in hardonic current

    m2f = mn     ! actually proton, particle-2
    m4f = mp     ! neutron, particle-4
    m3 = me
    U2 = Un
    U4 = Up
    dU2 = U2 - U4
    dmf = m2f - m4f
    Qmass = 0.5d0*( m2f**2 - m4f**2 )

    E1 = Enu
    E4 = E1 + E2 + E3
    E4f = E4 - U4
    E2f = E2 - U2
    P1 = E1
    if(E3**2 .ge. m3**2) then
      P3 = sqrt( E3**2 - m3**2 )
    else
      P3 = 0.0d0
    endif
    if((E2 - U2)**2 .ge. m2f**2) then
      P2 = sqrt( (E2 - U2)**2 - m2f**2 )
    else
      P2 = 0.0d0
    endif
    if((E4 - U4)**2 .ge. m4f**2) then
      P4 = sqrt( (E4 - U4)**2 - m4f**2 )
    else
      P4 = 0.0d0
    endif

    if( P1 .ne. P1) write(*,*) 'nan in hh1_b'
    if( P2 .ne. P2) write(*,*) 'nan in hh2_b'
    if( P3 .ne. P3) write(*,*) 'nan in hh3_b'
    if( P4 .ne. P4) write(*,*) 'nan in hh4_b'

    if ((P1+P2)>(P3+P4)) then
      Pmax = P3+P4
    else
      Pmax = P1+P2
    endif
    if ((abs(P1-P2))>(abs(P3-P4))) then
      Pmin = abs(P1-P2)
    else
      Pmin = abs(P3-P4)
    endif

    A = E1*E2f+0.5d0*(P1**2+P2**2)
    B = E3*E4f-0.5d0*(P3**2+P4**2)
    Ia = pi**2/15.0d0*(-3.0*(Pmax**5-Pmin**5) &
       - 10.0d0*(-A+B)*(Pmax**3-Pmin**3)+60.0d0*A*B*(Pmax-Pmin))

    if ((P1+P4)>(P3+P2)) then
      Pmax2 = P3+P2
    else
      Pmax2 = P1+P4
    endif
    if ((abs(P1-P4))>(abs(P3-P2))) then
      Pmin2 = abs(P1-P4)
    else
      Pmin2 = abs(P3-P2)
    endif

    D = E2f*E3+0.5d0*(P2**2+P3**2)
    Ct = E1*E4f-0.5d0*(P1**2+P4**2)
    Ib = pi**2/15.0d0*(-3.0*(Pmax2**5-Pmin2**5) &
       + 10.0d0*(-Ct+D)*(Pmax2**3-Pmin2**3)+60.0d0*Ct*D*(Pmax2-Pmin2))

    del0 = -(P1**2-P2**2)*(P3**2-P4**2)/4.0d0
    del1 = -E1*E3 + (-P1**2+P2**2-P3**2+P4**2)/4.0d0
    del2 = -0.25d0
    eps0 = E1*E2f + (P1**2+P2**2)/2.0d0
    eps1 = -0.5d0
    Ic = pi**2*4.0d0*(del2*eps1**2/7.0d0*(Pmax**7-Pmin**7) &
       + (2.0d0*del2*eps0*eps1+del1*eps1**2)/5.0d0*(Pmax**5-Pmin**5) &
       + (del2*eps0**2+2.0d0*del1*eps0*eps1+del0*eps1**2)/3.0d0 &
       * (Pmax**3-Pmin**3) &
       + (del1*eps0**2+2.0d0*del0*eps0*eps1) &
       * (Pmax-Pmin)-del0*eps0**2*(1.0d0/Pmax-1.0d0/Pmin))
    Ic = -Ic

    if ((P1+P3)>(P2+P4)) then
      Pmax3 = P2+P4
    else
      Pmax3 = P1+P3
    endif
    if ((abs(P1-P3))>(abs(P2-P4))) then
      Pmin3 = abs(P1-P3)
    else
      Pmin3 = abs(P2-P4)
    endif

    alp0 = (P1**2-P3**2)*(P2**2-P4**2)/4.0d0
    alp1 = E1*E2f + (P1**2+P2**2-P3**2-P4**2)/4.0d0
    alp2 = 0.25d0
    bet0 = -E1*E3 - (P1**2+P3**2)/2.0d0
    bet1 = 0.5d0
    Id = pi**2*4.0d0*(alp2*bet1**2/7.0d0*(Pmax3**7-Pmin3**7) &
       + (2.0d0*alp2*bet0*bet1+alp1*bet1**2)/5.0d0*(Pmax3**5-Pmin3**5) &
       + (alp2*bet0**2+2.0d0*alp1*bet0*bet1+alp0*bet1**2)/3.0d0 &
       * (Pmax3**3-Pmin3**3) &
       + (alp1*bet0**2+2.0d0*alp0*bet0*bet1)*(Pmax3-Pmin3) &
       - alp0*bet0**2*(1.0d0/Pmax3-1.0d0/Pmin3))

    Ie = pi**2/15.0d0*(3.0d0*(Pmax**5-Pmin**5) &
       - 20.0d0*A*(Pmax**3-Pmin**3)+60.0d0*A**2*(Pmax-Pmin))

    E = E1*E3 + 0.5d0*(P1**2+P3**2)
    Ig = pi**2/15.0d0*(3.0d0*(Pmax3**5-Pmin3**5) &
       - 20.0d0*E*(Pmax3**3-Pmin3**3)+60.0d0*E**2*(Pmax3-Pmin3))

    Ih = pi**2*4.0d0*(alp2*bet1/5.0d0*(Pmax3**5-Pmin3**5) &
       + (alp2*bet0+alp1*bet1)/3.0d0*(Pmax3**3-Pmin3**3) &
       + (alp1*bet0+alp0*bet1) &
       * (Pmax3-Pmin3)-alp0*bet0*(1.0d0/Pmax3-1.0d0/Pmin3))
    Ih = -Ih

    Ij = pi**2/15.0d0*(-10.0d0*(Pmax**3-Pmin**3)+60.0d0*A*(Pmax-Pmin))

    Ik = pi**2/15.0d0*(-10.0d0*(Pmax3**3-Pmin3**3)+60.0d0*E*(Pmax3-Pmin3))

    Il = pi**2/15.0d0*(60.0d0*(Pmax-Pmin))


    xE3 = -E3
    cplA = (gv0+anti*ga0)**2 &
         + anti*2.0d0*F2WM0*ga0*m2f/m2*(1.0d0-dmf/2.0d0/m2f)
    cplB = (gv0-anti*ga0)**2 &
         - anti*2.0d0*F2WM0*ga0*m2f/m2*(1.0d0-dmf/2.0d0/m2f)
    cplC = F2WM0**2/m2**2
    cplD =-F2WM0**2/m2**2
    cplE = F2WM0**2/m2**2*(-0.5d0*m3**2+dU2*(xE3-E1)-0.5d0*dU2**2)
    cplG = gv0*F2WM0*m2f/m2*(2.0d0-dmf/m2f) &
         + 0.5d0*F2WM0**2/m2**2*(m2f*m4f-Qmass+0.25d0*m3**2 &
         - dU2*(E1+E2f)-dU2**2/4.0d0)
    cplH = 0.5d0*F2WM0**2/m2**2 &
         * (2.0d0*Qmass+m3**2+dU2*(3.0d0*E1-xE3+2.0d0*E4f))

    JFF = -m3**2*(E1+0.5d0*E2f+0.5d0*E4f)+Qmass*(xE3-3.0d0*E1) &
        + 0.5d0*dU2*(E4f*(3.d0*xE3-5.d0*E1) &
        + E2f*(xE3+E1)+xE3**2-E1**2-2.d0*Qmass) &
        + dU2**2*(E1-xE3-E4f)+0.5d0*dU2**3
    JFF = JFF*dU2
    cplJ = gv0*F2WM0*dmf/m2*0.5d0*(m3**2-dU2*(E1+xE3)) &
         + F2WM0**2/m2**2*0.5d0*JFF

    kFF = -(m2f+3.d0*m4f)*m2f*0.25d0*m3**2 &
        + Qmass**2 &
        + Qmass*0.25d0*m3**2-m3**4/8.0d0 &
        + dU2*(0.5d0*Qmass*(3.*E1-E2f+xE3+3.d0*E4f) &
        + 0.25d0*m3**2*(2.d0*E2f+xE3+E1)+m2f*m4f*(xE3-E1)) &
        + dU2**2*( 0.25d0*(m2f**2-m2f*m4f-3.d0*Qmass) &
        + E4f*0.5d0*(2.d0*E1-E2f+xE3+E4f)+E2f*(0.5d0*E1-xE3)+0.5d0*E1**2 ) &
        + dU2**3*0.25d0*( -E1+2.d0*E2f-xE3-2.d0*E4f ) &
        + 0.125d0*dU2**4
    cplK = (ga0**2-gv0**2)*m2f*m4f &
         + gv0*F2WM0*m2f/m2*0.5d0*(-3.0d0*m3**2+4.0d0*dU2*(xE3-E1) &
         - dU2**2 &
         + dmf/m2f*(2.0d0*Qmass+m3**2+dU2*(E4f+2.0d0*E1-xE3))) &
         + 0.5d0*F2WM0**2/m2**2*KFF

    LFF = m3**2*(m2f+m4f)**2 &
        - 4.d0*Qmass**2+dU2*(-m3**2*(E2f+E4f+E1) &
        - 2.d0*xE3*(m2f**2+m2f*m4f) &
        + 2.d0*Qmass*(E2f-3.d0*E4f-E1)) &
        + 2.d0*dU2**2*( E2f*xE3+E4f*(E2f-E4f-E1) + Qmass ) &
        + dU2**3*(-E2f+E4f+E1)
    LFF = LFF*dU2*E1*0.25d0
    cplL = gv0*F2WM0*m2f/m2*dU2*E1 &
         * (m3**2-dU2*xE3-0.5d0*dmf/m2f &
         * (Qmass+0.5d0*m3**2+dU2*E4f-0.5d0*dU2**2)) &
         + 0.5d0*F2WM0**2/m2**2*LFF

    cplD = -cplD
    cplE = -cplE
    cplG = -cplG
    cplJ = -cplJ
    cplL = -cplL

    Ampsq = cplA*Ia &
          + cplB*Ib &
          + cplC*Ic &
          + cplD*Id &
          + cplE*Ie &
          + cplG*Ig &
          + cplH*Ih &
          + cplJ*Ij &
          + cplK*Ik &
          + cplL*Il

    return

  END SUBROUTINE Calc_Ampsq_D

  SUBROUTINE Ebounds_D( xpn3, xq0_min, xq0_max )
  !use CC_global
    implicit none

    real(dp) :: xpn3, q3_lim(4), tmp1,tmp2,xA,xB,xC, xq3
    integer  :: i
    real(dp) :: tmp3,tmp4,xq0_min,xq0_max,xq0, xq0_l,xq0_h

    En = sqrt( xpn3**2 + mn**2 ) + Un

    q3_lim = 0.d0
    tmp1 = Up-En-Enu
    tmp2 = Enu**2 + me**2 - mp**2 - xpn3**2

    xA = 4.d0*( tmp1**2-(Enu+xpn3)**2 )
    xB = 4.d0*( (tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*xpn3 )
    xC = 4.d0*tmp1**2*(xpn3**2+mp**2) - (tmp2-tmp1**2)**2
    q3_lim(1) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
    q3_lim(2) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )

    xA = 4.d0*( tmp1**2-(Enu-xpn3)**2 )
    xB = 4.d0*( -(tmp2-tmp1**2)*Enu + (tmp2+tmp1**2)*xpn3 )
    if( xB**2 - 4.0*xA*xC .ge. 0 ) then
      q3_lim(3) = abs( (xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
      q3_lim(4) = abs( (-xB + sqrt(xB**2-4.d0*xA*xC))/(2.d0*xA) )
    end if

    xq0_min = 1.d10
    xq0_max = -1.d10
    do i=1,4
      xq3 = q3_lim(i)
      tmp1 = Enu + sqrt( (Enu+xq3)**2 + me**2 )
      tmp2 = Enu + sqrt( (Enu-xq3)**2 + me**2 )
      tmp3 = sqrt( (xpn3+xq3)**2+mp**2) + Up - En
      tmp4 = sqrt( (xpn3-xq3)**2+mp**2) + Up - En
      if( (abs((tmp1-tmp3)/tmp3)<1.d-6).or. &
          (abs((tmp1-tmp4)/tmp4)<1.d-6) ) then
        xq0 = tmp1
      end if
      if( (abs((tmp2-tmp3)/tmp3)<1.d-6).or. &
          (abs((tmp2-tmp4)/tmp4)<1.d-6) ) then
        xq0 = tmp2
      end if
      xq0_min = min( xq0,  xq0_min)
      xq0_max = max( xq0,  xq0_max)
    end do
    if( abs(xq0_max-xq0_min)<1.d-15 ) then
      if(opt0.eq.3) then ! inverse decay
        xq0_max = max(Enu + mue + 50.0*Tem, xq0_min + 50.0*Tem)
      else if(opt0.eq.4) then ! decay
        xq0_max = max(xq0_min+50.0*Tem, sqrt(ppmax**2+mp**2)+Up-En)
      end if
    end if

    xq0_l = Enu + me
    tmp1 = sqrt( (xpn3+Enu)**2+mp**2) + Up - En
    tmp2 = sqrt( (xpn3-Enu)**2+mp**2) + Up - En
    if( xq0_l .le. tmp1 .and. xq0_l .ge. tmp2) then
      xq0_min = min( xq0_l, xq0_min)
    end if

    xq0_h = mp + Up - En
    tmp1 = Enu + sqrt( (Enu+xpn3)**2 + me**2 )
    tmp2 = Enu + sqrt( (Enu-xpn3)**2 + me**2 )
    if( xq0_h .le. tmp1 .and. xq0_h .ge. tmp2) then
      xq0_min = min( xq0_h, xq0_min)
    end if

  END SUBROUTINE Ebounds_D

  SUBROUTINE Range_pn_D()
  !use CC_global
    implicit none
    real(dp) :: tmp,tmp1,tmp2,tmp3,del,At,xpn3,mpt
    !real(dp),external:: test1_D
    integer  :: i
    real(dp) :: F0, xpn, xpn1, xF, Fmax
    real(dp) :: Lmax,Lmin,Hmin,Hmax,xq

    pnmin = 0.d0
    if((Tfac*Tem + mup - Up)**2 .ge. mp**2) then
      ppmax = sqrt( (Tfac*Tem + mup - Up)**2 - mp**2 )
    else
      ppmax = 0.0d0
    endif
    select case(opt0)
      case (3)  ! inverse decay
        if((Tfac*Tem + mun - Un)**2 .ge. mn**2) then
          pnmax = sqrt( (Tfac*Tem + mun - Un)**2 - mn**2 )
        else
          pnmax = 0.0d0
        endif
      case (4)  ! decay
        pnmin = sqrt( (max(mn,mun-Tfac*Tem-Un))**2 - mn**2 )
        pnmax = sqrt( (sqrt(ppmax**2+mp**2)+Up-Un)**2-mn**2 )
    end select

    pemax = sqrt( max( (me+20.d0*Tem)**2,(Tfac*Tem + mue)**2 )-me**2)
    pemin = sqrt( max( me**2,(-Tfac*Tem + mue)**2 ) - me**2 )

    mpt = mp - me
    if( mpt .gt. mn) then
      tmp = mn/(mpt-mn)*Enu
    else
      tmp = pnmax
    end if
    Fmax = sqrt( (tmp+Enu)**2+mpt**2 )+Up-sqrt(tmp**2+mn**2)-Un

    if(Enu>=Fmax) then
      pnmax = -1000.d0
      goto 12
    end if

    F0 = sqrt(Enu**2 + mpt**2)+Up-mn-Un
    if( F0 .ge. Enu .and. Up .ge. Un) goto 12


    del = 1.d-5
    tmp = Enu - Up + Un
    At = Enu**2 + mpt**2 - mn**2 - tmp**2
    if(tmp**2 .eq. Enu**2) then
      xpn =(4.d0*(Enu*mn)**2-(mpt**2-mn**2)**2)/( 4.d0*(mpt**2-mn**2)*Enu)
      if( test1_D(xpn-del)*test1_D(xpn+del) <= 0.d0 .and. xpn > 0.0d0) then
        pnmin = xpn
      end if

    else
      tmp1 = At**2 + 4.d0*( (mn*Enu)**2 - (tmp*mn)**2 )
      tmp2 =(Enu*At - abs(tmp)*sqrt(abs(tmp1)))*0.5d0/(tmp**2-Enu**2)
      tmp3 =(Enu*At + abs(tmp)*sqrt(abs(tmp1)))*0.5d0/(tmp**2-Enu**2)

      if(tmp2 > 0.d0) then
        if( (test1_D(tmp2-del)>0.d0 .and. test1_D(tmp2+del)<0.d0) &
        .or.(test1_D(tmp2-1.d0)>0.d0 .and. test1_D(tmp2+1.d0)<0.d0))then
          pnmax = tmp2
        else if( (test1_D(tmp2-del)<0.d0 .and. test1_D(tmp2+del)>0.d0) &
             .or.(test1_D(tmp2-1.d0)<0.d0 .and. test1_D(tmp2+1.d0)>0.d0))then
          pnmin = tmp2
        end if
      end if

      if(tmp3 > 0.d0) then
        if( (test1_D(tmp3-del)>0.d0 .and. test1_D(tmp3+del)<0.d0) &
        .or.(test1_D(tmp3-1.d0)>0.d0 .and. test1_D(tmp3+1.d0)<0.d0))then
          pnmax = tmp3
        else if( (test1_D(tmp3-del)<0.d0 .and. test1_D(tmp3+del)>0.d0) &
             .or.(test1_D(tmp3-1.d0)<0.d0 .and. test1_D(tmp3+1.d0)>0.d0))then
          pnmin = tmp3
        end if
      end if

    end if

12    continue

  END SUBROUTINE Range_pn_D




  FUNCTION test1_D(xPn3)
  !use CC_global
    implicit none
    real(dp) :: xpn3,test1_D,mpt

    mpt = mp - me
    test1_D = sqrt((xpn3+Enu)**2 + mpt**2) &
            - sqrt(xpn3**2+mn**2)+Up-Un-Enu
    return
  END FUNCTION test1_D

  SUBROUTINE gaulegg(x1,x2,x,w,n)

    implicit none

    integer  :: n
    real(dp) :: x1,x2,x(n),w(n)
    real(dp), parameter :: EPS = 3.0d-14
    integer  :: i,j,m
    real(dp) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)

    do 12 i=1,m
      z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
      p1=1.d0
      p2=0.d0
      do 11 j=1,n
        p3=p2
        p2=p1
        p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue

    return

  END SUBROUTINE gaulegg

END MODULE CC_module_weaklib

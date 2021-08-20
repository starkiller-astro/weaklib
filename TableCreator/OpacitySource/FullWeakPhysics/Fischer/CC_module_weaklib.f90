MODULE CC_global_weaklib

  USE wlKindModule, ONLY: dp

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

   integer:: anti,opt0                           

END MODULE         

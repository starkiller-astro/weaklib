PROGRAM wlGreyVariablesTEST

  USE wlKindModule, ONLY: dp
  USE GreyVariables, ONLY: &
        GreyMomentWithGaussianQuadrature, &
        GreyOpacityWithGaussianQuadrature

implicit none
    

   REAL(dp)                :: energy, rho, T, Z, A, chem_e, chem_n, &
                              chem_p, xheavy, xn, xp, &
                              bufferquad1, bufferquad2, bufferquad3,&
                              bufferquad4
   INTEGER                 :: nquad = 6 

PRINT*, "Gaussian Quadrature test"
   CALL GreyMomentWithGaussianQuadrature&
             ( nquad, 10.0_dp, bufferquad2,&
               "GreyMoment_Energy ", .TRUE. )

!   CALL GreyOpacityWithGaussianQuadrature&
!             ( nquad, 10.0_dp, &
!                    rho, T, Z, A, chem_e, chem_n,&
!                    chem_p, xheavy, xn, xp,&
!                    bufferquad3,"GreyOpacity_Energy ",.TRUE. )

  WRITE (*,*) "Done test"

!=============================================================

END PROGRAM wlGreyVariablesTEST

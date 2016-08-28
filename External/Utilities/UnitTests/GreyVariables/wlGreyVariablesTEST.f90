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
             ( nquad, 1.0_dp, bufferquad2,&
               "GreyMoment_Number ", .TRUE. )
   WRITE(*,*) "Expects 4.32833"

   CALL GreyMomentWithGaussianQuadrature&
             ( nquad, 1.0_dp, bufferquad2,&
               "GreyMoment_Energy ", .TRUE. )
   WRITE(*,*) "Expects 14.3894"

   CALL GreyOpacityWithGaussianQuadrature&
             ( nquad, 1.0_dp, &
                    rho, T, Z, A, chem_e, chem_n,&
                    chem_p, xheavy, xn, xp,&
                    bufferquad3,"GreyOpacity_Number ",.TRUE. )

   CALL GreyOpacityWithGaussianQuadrature&
             ( nquad, 1.0_dp, &
                    rho, T, Z, A, chem_e, chem_n,&
                    chem_p, xheavy, xn, xp,&
                    bufferquad3,"GreyOpacity_Energy ",.TRUE. )

  WRITE (*,*) "Done test"

!=============================================================

END PROGRAM wlGreyVariablesTEST

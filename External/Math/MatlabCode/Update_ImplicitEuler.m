function [ Nnew, nIter ]...
  = Update_ImplicitEuler( Nold, dt, R_In, R_Out, theta, N_g, Tol_N )

  [ Nnew, nIter ]...
    = Newton( Nold, dt, R_In, R_Out, theta, N_g, Tol_N );

end


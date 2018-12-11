function [ Nnew, nIter ]...
  = Update_SSPRK3( Nold, dt, R_In, R_Out, theta, N_g )

  nIter = 0;
  
  c21 = 3.0/4.0;
  c22 = 1.0/4.0;
  
  c31 = 1.0/3.0;
  c32 = 2.0/3.0;

  % Stage 1:
  N1...
    = ( eye( N_g ) + dt .* L_FUN( Nold, R_In, R_Out, theta, N_g ) ) * Nold;
  % Stage 2:
  N2...
    = c21 .* Nold + c22 .* ( ( eye( N_g ) + dt .* L_FUN( N1, R_In, R_Out, theta, N_g ) ) * N1 );
  % Stage 3:
  Nnew...
    = c31 .* Nold + c32 .* ( ( eye( N_g ) + dt .* L_FUN( N2, R_In, R_Out, theta, N_g ) ) * N2 );

end


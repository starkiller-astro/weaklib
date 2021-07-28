function [ Nnew, nIter ]...
  = Update_ForwardEular( Nold, dt, R_In, R_Out, theta, N_g )

  nIter = 0;

  Nnew...
    = Nold + dt .* L_FUN( Nold, R_In, R_Out, theta, N_g ) * Nold;

end


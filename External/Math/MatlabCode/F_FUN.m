function [ F ] = F_FUN( Nnew, Nold, dt, R_In, R_Out, theta, N_g )

  F = ( eye( N_g ) - dt * L_FUN( Nnew, R_In, R_Out, theta, N_g ) ) * Nnew - Nold;

end


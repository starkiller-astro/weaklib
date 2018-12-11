function [ L ] = L_FUN( N, R_In, R_Out, theta, N_g )

  k = diag( R_Out * ones( N_g, 1 ) + theta .* ( R_In - R_Out ) * N );
  L = R_In - k;
  
end


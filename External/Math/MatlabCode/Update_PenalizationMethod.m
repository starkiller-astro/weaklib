function [ Nnew, nIter ]...
  = Update_PenalizationMethod...
      ( Nold, N_Eq, dt, R_In, R_Out, theta, N_g, dV2, dV3 )

  dC = L_FUN( N_Eq, R_In, R_Out, theta, N_g )...
       - theta .* diag( N_Eq ) * ( R_In - R_Out );
  
  % Explicit Step:
  Nnew...
    = Nold...
        + dt .* ( L_FUN( Nold, R_In, R_Out, theta, N_g) * Nold...
                  + dC * ( N_Eq - Nold ) );
              
  % Implicit Step:
  Nnew = ( eye( N_g ) - dt .* dC ) \ ( Nnew - dt.* dC * N_Eq );
  
  nIter = 1;

end


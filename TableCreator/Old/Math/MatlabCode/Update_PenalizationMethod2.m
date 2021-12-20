function [ Nnew, nIter ]...
  = Update_PenalizationMethod2...
      ( Nold, N_Eq, dt, R_In, R_Out, theta, N_g, Correct )

  dC = L_FUN( N_Eq, R_In, R_Out, theta, N_g )...
       - theta .* diag( N_Eq ) * ( R_In - R_Out );

  N_0 = Nold;

  % Stage 1:

  % Explicit Step:
  N_1...
    = N_0...
        + dt .* ( L_FUN( N_0, R_In, R_Out, theta, N_g) * N_0...
                  + dC * ( N_Eq - N_0 ) );

  % Implicit Step:
  N_1 = ( eye( N_g ) - dt .* dC ) \ ( N_1 - dt.* dC * N_Eq );

  % Stage 2:

  % Explicit Step:
  N_2...
    = N_1...
        + dt .* ( L_FUN( N_1, R_In, R_Out, theta, N_g) * N_1...
                  + dC * ( N_Eq - N_1 ) );

  % Implicit Step:
  N_2 = ( eye( N_g ) - dt .* dC ) \ ( N_2 - dt.* dC * N_Eq );

  % Combine Steps:
  N_2 = 0.5 * ( N_0 + N_2 );
  
  if( Correct )

    % Correction Step:
    F = L_FUN( N_2, R_In, R_Out, theta, N_g) * N_2...
        + dC * ( N_Eq - N_2 );
    
    Nnew = ( eye( N_g ) + dt^2 .* dC * dC )...
           \ ( N_2 - dt^2 .* dC * ( F - dC * N_Eq ) );
       
  else
    
    Nnew = N_2;
    
  end

  nIter = 1;

end


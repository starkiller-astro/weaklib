function [ J ] = F_JAC( N, dt, R_In, R_Out, theta, N_g )

  J = eye( N_g )...
      - dt .* ( L_FUN( N, R_In, R_Out, theta, N_g )...
                - theta .* diag( N ) * ( R_In - R_Out ) );

end


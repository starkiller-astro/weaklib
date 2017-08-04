function [ Nnew, nIter ]...
  = Update_PenalizationMethod_dCconstant...
      ( Nold, N_Eq, dt, R_In, R_Out, theta, N_g )
% this function update N and nIter
% for fixed tempereature case
% with penalization method and P(Neq) = - lamp, lamp > 0

  dC = L_FUN( N_Eq, R_In, R_Out, theta, N_g )...
       - theta .* diag( N_Eq ) * ( R_In - R_Out );
 
  lamp = max(abs(eig(dC)));
  
  Nnew = Nold + ( dt / ( 1.0 + lamp * dt ) ) ...
      .* L_FUN( Nold, R_In, R_Out, theta, N_g) * Nold;
  
  nIter = 1;

end


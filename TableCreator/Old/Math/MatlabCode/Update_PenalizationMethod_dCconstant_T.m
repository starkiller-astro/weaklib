function [ Nnew, intEnew, Tnew, N_Eqnew, R_In_H, R_Out_H, nIter, dnormNES ]...
  = Update_PenalizationMethod_dCconstant_T...
      ( Nold, N_Eqold, intEold, eC, dV, h3c3, ...
    E, rho, Ye, eosD, eosT, eosY, intE_table, eosOSintE, ...
    Eta, NesR_0, NesOS, Me, eosOS,...
    chemTab, chemOS, dt, R_In, R_Out, theta, N_g )
% this function update N, intE, T, N_Eq, R_In_H, R_Out_H, nIter, dnormNES
% for changeable tempereature case
% with penalization method and P(Neq) = - lamp, lamp > 0
% Modified Update_PenalizationMethod_dCconstant.m
% and Update_PenalizationMethod_T.m

    W      = eC .* dV;
   C_0     = 1.6021773E-6; % MeV to erg
 
   dC      = L_FUN( N_Eqold, R_In, R_Out, theta, N_g )...
             - theta .* diag( N_Eqold ) * ( R_In - R_Out );   
   lamp    = max(abs(eig(dC)));
   
   %
   % Update N 
   Nnew    = Nold + ( dt / ( 1.0 + lamp * dt ) ) ...
             .* L_FUN( Nold, R_In, R_Out, theta, N_g) * Nold;
   
   %
   % Update specific internal energy
   RHS     = Nnew - Nold;
   intEnew = intEold - ( C_0 / rho ) * ( W' * RHS ) / h3c3;
   
   % Update Temperature
   Tnew    = ComputeTempFromIntEnergy_Bisection...
            (rho, intEnew, Ye, eosD, eosT, eosY, intE_table, eosOSintE );
   
   % Update NES Kernel
   [R_In_H, R_Out_H] = Update_NESKernel...
       ( eC, dV, rho, Tnew, Ye, E, eosD, eosT, eosY, ...
       Eta, NesR_0, NesOS, Me, eosOS);
   
   % Update N_Eq
   N_Eqnew = Update_Neq_FD...
       ( rho, Tnew, Ye, eosD, eosT, eosY, chemTab, chemOS, eC);
   
   nIter   = 1;
   
   dnormNES = max(max( abs( R_In_H - R_In ) ./ R_In ));
   
end


function [ Nnew, intEnew, Tnew, N_Eqnew, R_In_H, R_Out_H, nIter ]...
  = Update_ForwardEular_T( Nold, intEold, eC, dV, h3c3, ...
    E, rho, Ye, eosD, eosT, eosY, intE_table, eosOSintE, ...
    Eta, NesR_0, NesOS, Me, eosOS,...
    chemTab, chemOS, dt, R_In, R_Out, theta, N_g )
% this function update N, intE, T, N_Eq, R_In_H, R_Out_H, and nIter
% for changeable tempereature
% with forward Euler method

    W      = eC .* dV;
   C_0     = 1.6021773E-6; % MeV to erg
   
   % Update N 
   Nnew...
           = Nold + dt .* L_FUN( Nold, R_In, R_Out, theta, N_g ) * Nold;
   
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
   
  nIter = 1;


end


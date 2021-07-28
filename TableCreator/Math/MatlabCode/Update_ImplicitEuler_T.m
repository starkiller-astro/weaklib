 function [ Nnew, intEnew, Tnew, R_In_H, R_Out_H, nIter]...
  = Update_ImplicitEuler_T( Nold, intEold, eC, dV, h3c3, ...
    E, rho, Ye, eosD, eosT, eosY, intE_table, eosOSintE, ...
    Eta, NesR_0, NesOS, Me, eosOS,...
    dt, R_In, R_Out, theta, N_g, Tol_N)

   W = eC .* dV;
   C_0 = 1.6021773E-6; % MeV to erg 

   [ Nnew, nIter ]...
     = Newton( Nold, dt, R_In, R_Out, theta, N_g, Tol_N );
 
   % Update Specific Internal Energy
   RHS     = Nnew - Nold;
   intEnew = intEold - ( C_0 / rho ) * ( W' * RHS ) / h3c3;

   % Update Temperature
   Tnew = ComputeTempFromIntEnergy_Bisection...
    (rho, intEnew, Ye, eosD, eosT, eosY, intE_table, eosOSintE );

   % Update NES Kernel
   [R_In_H, R_Out_H] = Update_NESKernel...
    ( eC, dV, rho, Tnew, Ye, E, eosD, eosT, eosY, ...
    Eta, NesR_0, NesOS, Me, eosOS);

end


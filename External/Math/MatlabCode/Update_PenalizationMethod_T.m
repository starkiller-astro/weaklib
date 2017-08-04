function [ Nnew, intEnew, Tnew, N_Eqnew, R_In_H, R_Out_H, nIter ]...
    = Update_PenalizationMethod_T...
    ( Nold, N_Eq, intEold, eC, dV, h3c3, ...
    E, rho, Ye, eosD, eosT, eosY, intE_table, eosOSintE, ...
    Eta, NesR_0, NesOS, Me, eosOS,...
    chemTab, chemOS, dt, R_In, R_Out, theta, N_g)

  W = eC .* dV;
  C_0 = 1.6021773E-6; % MeV to erg

  dC = L_FUN( N_Eq, R_In, R_Out, theta, N_g )...
    - theta .* diag( N_Eq ) * ( R_In - R_Out );

% Explicit Step:

  Nnew...
    = Nold...
    + dt .* ( L_FUN( Nold, R_In, R_Out, theta, N_g) * Nold...
    + dC * ( N_Eq - Nold ) );

% Implicit Step:
  Nnew = ( eye( N_g ) - dt .* dC ) \ ( Nnew - dt.* dC * N_Eq );

%
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

% Update N_Eq
  N_Eqnew = Update_Neq_FD...
    ( rho, Tnew, Ye, eosD, eosT, eosY, chemTab, chemOS, eC);

  nIter = 1;
end


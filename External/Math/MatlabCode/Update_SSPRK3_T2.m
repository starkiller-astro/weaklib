function [ Nnew, intEnew, Tnew, R_In_H, R_Out_H, nIter ]...
    = Update_SSPRK3_T2( Nold, intEold, eC, dV, h3c3, ...
    E, rho, Ye, eosD, eosT, eosY, intE_table, eosOSintE, ...
    Eta, NesR_0, NesOS, Me, eosOS,...
    dt, R_In, R_Out, theta, N_g )

nIter = 0;

c21 = 3.0/4.0;
c22 = 1.0/4.0;

c31 = 1.0/3.0;
c32 = 2.0/3.0;

W = eC .* dV;
C_0 = 1.6021773E-6;
%
% Stage 1:
RHS   = L_FUN( Nold, R_In, R_Out, theta, N_g ) * Nold;
N1    = Nold + dt .* RHS;
intE1 = intEold - ( C_0 / rho ) * dt * ( W' * RHS ) / h3c3;

% Update Temperature
T1 = ComputeTempFromIntEnergy_Bisection...
    (rho, intE1, Ye, eosD, eosT, eosY, intE_table, eosOSintE );

% Update NES Kernel
[R_In, R_Out] = Update_NESKernel...
    ( eC, dV, rho, T1, Ye, E, eosD, eosT, eosY, ...
    Eta, NesR_0, NesOS, Me, eosOS);

%
% Stage 2:
RHS = L_FUN( N1, R_In, R_Out, theta, N_g ) * N1;
N2  = c21 .* Nold + c22 .* ( N1 + dt .* RHS );
intE2 = c21 * intEold + c22 * ( intE1 - ( C_0 / rho ) * dt * ( W' * RHS ) / h3c3 );

% Update Temperature
T2 = ComputeTempFromIntEnergy_Bisection...
    (rho, intE2, Ye, eosD, eosT, eosY, intE_table, eosOSintE );

% Update NES Kernel
[ R_In, R_Out ]...
  = Update_NESKernel...
      ( eC, dV, rho, T2, Ye, E, eosD, eosT, eosY, Eta, NesR_0, NesOS, Me, eosOS );

%
% Stage 3:
RHS = L_FUN( N2, R_In, R_Out, theta, N_g ) * N2;
Nnew = c31 .* Nold + c32 .* ( N2 + dt .* RHS );
intEnew = c31 * intEold + c32 * ( intE2 - ( C_0 / rho ) * dt * ( W' * RHS ) / h3c3 );

% Update Temperature
Tnew = ComputeTempFromIntEnergy_Bisection...
    (rho, intEnew, Ye, eosD, eosT, eosY, intE_table, eosOSintE );

% Update NES Kernel
[R_In_H, R_Out_H] = Update_NESKernel...
    ( eC, dV, rho, Tnew, Ye,  E, eosD, eosT, eosY, ...
    Eta, NesR_0, NesOS, Me, eosOS);

end


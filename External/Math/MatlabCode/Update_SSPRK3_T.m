function [ Nnew, intEnew, Tnew, R_In_H, R_Out_H, nIter ]...
    = Update_SSPRK3_T( Nold, intEold, eC, dV, h3c3, ...
    E, rho, Ye, eosD, eosT, eosY, intE_table, eosOSintE, ...
    Eta, NesR_0, NesOS, Me, eosOS,...
    dt, R_In, R_Out, theta, N_g )

nIter = 0;

c21 = 3.0/4.0;
c22 = 1.0/4.0;

c31 = 1.0/3.0;
c32 = 2.0/3.0;

%
% Stage 1:
N1...
    = Nold + dt .* L_FUN( Nold, R_In, R_Out, theta, N_g ) * Nold;

% Compute C(R_n, N_n)*dt
C_Fun = N1 - Nold;

% Update specific internal energy density
dEMeV = - sum( C_Fun .* eC .* dV )/h3c3; % [MeV cm-3]
intE1 = intEold + dEMeV * 1.6021773E-6/ rho; % [erg per gram]

% Update Temperature
T1 = ComputeTempFromIntEnergy_Bisection...
    (rho, intE1, Ye, eosD, eosT, eosY, intE_table, eosOSintE );

% Update NES Kernel
[R_In, R_Out] = Update_NESKernel...
    ( eC, dV, rho, T1, Ye, E, eosD, eosT, eosY, ...
    Eta, NesR_0, NesOS, Me, eosOS);

%
% Stage 2:
N2...
    = c21 .* Nold + c22 .* ( N1 + dt .* L_FUN( N1, R_In, R_Out, theta, N_g ) * N1 );

% Compute C(R_n, N_n)*dt
C_Fun = N2 - Nold;

% Update specific internal energy density
dEMeV = - sum( C_Fun .* eC .* dV )/h3c3; % [MeV cm-3]
intE2 = intEold + dEMeV * 1.6021773E-6/ rho; % [erg per gram]

% Update Temperature
T2 = ComputeTempFromIntEnergy_Bisection...
    (rho, intE2, Ye, eosD, eosT, eosY, intE_table, eosOSintE );

% Update NES Kernel
[R_In, R_Out] = Update_NESKernel...
    ( eC, dV, rho, T2, Ye, E, eosD, eosT, eosY, ...
    Eta, NesR_0, NesOS, Me, eosOS);

%
% Stage 3:
Nnew...
    = c31 .* Nold + c32 .* ( N2 + dt .* L_FUN( N2, R_In, R_Out, theta, N_g ) * N2 );

% Compute C(R_n, N_n)*dt
C_Fun = Nnew - Nold;

% Update specific internal energy density
dEMeV = - sum( C_Fun .* eC .* dV )/h3c3; % [MeV cm-3]
intEnew = intEold + dEMeV * 1.6021773E-6/ rho; % [erg per gram]

% Update Temperature
Tnew = ComputeTempFromIntEnergy_Bisection...
    (rho, intEnew, Ye, eosD, eosT, eosY, intE_table, eosOSintE );


% Update NES Kernel
[R_In_H, R_Out_H] = Update_NESKernel...
    ( eC, dV, rho, Tnew, Ye,  E, eosD, eosT, eosY, ...
    Eta, NesR_0, NesOS, Me, eosOS);

end


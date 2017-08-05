function [ Tnew ] = Update_Temperature...
    ( eC, dV, C_Fun, rho, T, Ye, eosD, eosT, eosY, intE)

c  = 2.99792d+10;  % [cm s-1]
h  = 4.13567d-21;  % Planck's constant [MeV s]
h3c3 = h^3 * c^3;  % [MeV3 cm3]

Cv = ComputeHeatCapacity(rho, T, Ye, eosD, eosT, eosY, intE);% [MeV T-1]


Tnew = T - sum( C_Fun .* eC .* dV )/Cv/h3c3; % [K]

end


function [ Neq ] = Update_Neq_FD...
    ( rho, T_K, Ye, eosD, eosT, eosY, chemTab, chemOS, eC)
% This function provides the equilibrium number density
% distribution function as a Fermi-Dirac function 
% for a given density, temperature, electron fraction,
% and chemical potential table 

  % find the neutrino chemical potential 
  chem = interpolateEos( rho, T_K, Ye, eosD, eosT, eosY, chemTab, chemOS );
   
  % exponent of the FD function
  kmev = 8.61733d-11;
  beta = ( eC - chem ) ./ ( T_K * kmev );
  
  %
  Neq  = 1.0 ./ ( exp( beta ) + 1.0 );


end


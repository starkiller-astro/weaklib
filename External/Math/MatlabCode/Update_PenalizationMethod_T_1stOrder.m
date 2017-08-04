function [ Nnew, intEnew, Tnew, N_Eqnew, R_In_H, R_Out_H, nIter, dnormNES ]...
    = Update_PenalizationMethod_T_1stOrder...
    ( Nold, N_Eqold, intEold, eC, dV, h3c3, ...
    E, rho, Told, Ye, eosD, eosT, eosY, intE_table, eosOSintE, ...
    Eta, NesR_0, NesOS, Me, eosOS,...
    chemmuTable, chemmuOS, dt, R_In, R_Out, theta, N_g )
% this function update N, intE, T, N_Eq, R_In_H, R_Out_H, nIter, dnormNES
% for changeable tempereature
% with first order penalization method and newton's method

  W = eC .* dV;       % [ Energy ^ 4 ] = [ MeV ^ 4 ]
  C_0 = 1.6021773E-6; % convert MeV to erg

  dC = L_FUN( N_Eqold, R_In, R_Out, theta, N_g )...
    - theta .* diag( N_Eqold ) * ( R_In - R_Out ); % Dimensionless

  A  =  eye( N_g) + dt .* dC; % Dimensionless
  tildeW = W.' * ( A \ dC ) .* dt; % [ s * MeV ^ 4 ]
  
  GUN = L_FUN( Nold, R_In, R_Out, theta, N_g) * Nold ...
        + dC * ( Nold - N_Eqold ); % [ N ]
    
  fun = rho * intEold / C_0  - W.' * ( A \ Nold ) / h3c3 ...
        + W.' * Nold / h3c3  - W.' * ( A \ GUN ) .* dt / h3c3 ; % [ MeV ]

  Tnew = Told;
  N_Eqnew = N_Eqold;
  intEnew = intEold;
  
  converged = false;
  nIter     = 0;
  Tol = 1E-12;
  
  while ( not( converged ) )

    Fun = rho * intEnew / C_0 + tildeW * N_Eqnew / h3c3 - fun;
    
    [ ~, ~, dintEdT, ~ ] = ...
        interpolateDifferentiateEos(rho, Tnew, Ye, eosD, eosT, eosY, intE_table, eosOSintE);
    
    kmev = 8.61733d-11;
    chemmu = interpolateEos( rho, Tnew, Ye, eosD, eosT, eosY, chemmuTable, chemmuOS );
    [ ~, ~, dchemdT, ~ ] = ...
        interpolateDifferentiateEos(rho, Tnew, Ye, eosD, eosT, eosY, chemmuTable, chemmuOS);
    
    Mid1 = power( ( exp( ( eC - chemmu ) / ( kmev * Tnew ) ) + theta ), -2.0 );
    Mid2 = kmev * Tnew* ( eC - chemmu ) .* exp( ( eC - chemmu )/( kmev * Tnew ) );
    Mid3 = ( eC - chemmu ) .* kmev - dchemdT/(kmev * Tnew ); 
    dNEqdT  = - Mid1 .* Mid2 .* Mid3;
    
    dFun =  rho * dintEdT / C_0 + tildeW * dNEqdT / h3c3;
    
    dT =  - dFun \ Fun;
       
    Tnew = Tnew + dT;
    
    intEnew = interpolateEos( rho, Tnew, Ye, eosD, eosT, eosY, intE_table, eosOSintE);
    
    N_Eqnew = Update_Neq_FD( rho, Tnew, Ye, eosD, eosT, eosY, chemmuTable, chemmuOS, eC);
    
    nIter = nIter + 1;
    
    if( abs( dT / Tnew ) < Tol )
      converged = true;
    end

  end

% Update Nnew  
  Nnew = A \ ( Nold + dt .* ( GUN  + dC * N_Eqnew ) );
  
% Update NES Kernel
  [R_In_H, R_Out_H] = Update_NESKernel...
    ( eC, dV, rho, Tnew, Ye, E, eosD, eosT, eosY, ...
    Eta, NesR_0, NesOS, Me, eosOS);

dnormNES = max(max( abs( R_In_H - R_In ) ./ R_In ));
end


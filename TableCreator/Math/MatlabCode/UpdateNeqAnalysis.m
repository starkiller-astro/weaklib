clear all

Model = '002';

theta = 1.0;

N_g = 40; % Number of Energy Groups

TableName = 'SFHo';

[ eosD, eosT, eosY, nD, nT, nY, P, S, intE, Me, Mp, Mn, Xp, Xn, Xa, Xh,...
    Zh, Ah, Eh, Eth, Gm, eosOS ] = ...
    readEosTable(['wl-EOS-' TableName '-25-40-100.h5']);
[ E, T, Eta, NesR_0, NesR_1, NesOS] = ...
    readNesScatteringOpacityTable(['wl-Op-' TableName '-25-40-100.h5' ]);

dmnp = 1.29333d+00;
Mmu = Mp - Mn + Me -dmnp;
MmuOS = abs( min( min( min( Mmu ) ) ) );
Mmu_log10 = log10( Mmu + MmuOS + 1e-100 );

[ rho, T0, Ye, eC, dV, R_In, R_Out, N_Eq, intE0] = InitializeNES_Interp_rhoTYe...
    (Model, N_g, E, eosD, eosT, eosY, Eta, NesR_0, NesOS(1), Me, eosOS(4), intE, eosOS(3),...
    Mmu, MmuOS);


% Multiply Rates with Momentum Space Volume Element:
R_In_H  = R_In  * diag( dV ) ;
R_Out_H = R_Out * diag( dV ) ;

isNeq( N_Eq, R_In_H, R_Out_H, theta, N_g );


%%% enforce the R_In and R_Out symmetry
[ R_In, R_Out ] = EnforceRinoutSym( R_Out, eC, T0 );

% Multiply Rates with Momentum Space Volume Element:
R_In_H  = R_In  * diag( dV ) ;
R_Out_H = R_Out * diag( dV ) ;


isNeq( N_Eq, R_In_H, R_Out_H, theta, N_g );



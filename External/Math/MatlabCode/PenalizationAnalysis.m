clear all
close all

Model = '003';

theta = 1.0;
N_g = 40; % Number of Energy Groups
[ eC, dV, R_In, R_Out, N_Eq ] = InitializeNES( Model, N_g );

% Multiply Rates with Momentum Space Volume Element:
R_In_H  = R_In  * diag( dV );
R_Out_H = R_Out * diag( dV );
dV3 = eC .* dV;

%
% Initial Distribution:
beta = 0.0;
N_0 = 7.5d-1 .* exp( - ( 1.0d2 - eC ).^2 / 1.0d2 );
N_0 = beta .* N_Eq + (1.0-beta) .* N_0;

Lmat = L_FUN( N_0,  R_In_H, R_Out_H, theta, N_g );
Pmat = L_FUN( N_Eq, R_In_H, R_Out_H, theta, N_g )...
       - diag( N_Eq ) * ( R_In_H - R_Out_H );
Kvec = ( R_In_H - R_Out_H ) * ( N_0 - N_Eq );

dt_K = min( 1.0 ./ abs( Kvec + 1.0d-16 ) );
dt   = min( [ 1.0 dt_K ] );

% Update Matrix from Nold to Nnew:
Wmat = ( eye(N_g) - dt .* Pmat )...
       \ ( eye(N_g) + dt .* ( Lmat - Pmat ) );

figure(1)
semilogx...
  ( eC, Lmat*N_0, '-o',...
    eC, Pmat*(N_0-N_Eq), '-x',...
    eC, Kvec, '-*' )

figure(2)
semilogx...
  ( eC, N_0,      '-o',...
    eC, Wmat*N_0, '-+' )

figure(3)
plot( eC, Wmat )

figure(4)
imagesc( Wmat )
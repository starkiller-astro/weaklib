clear all
% close all

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
beta = 0.9;
N_0 = 7.5d-1 .* exp( - ( 1.0d2 - eC ).^2 / 1.0d2 );
N_0 = beta .* N_Eq + (1.0-beta) .* rand(N_g,1) .* N_Eq;

Lmat   = L_FUN( N_0,  R_In_H, R_Out_H, theta, N_g );
Pmat   = L_FUN( N_Eq, R_In_H, R_Out_H, theta, N_g )...
         - diag( N_Eq ) * ( R_In_H - R_Out_H );
F_plus = R_In_H * N_0;
Kappa  = R_Out_H * ones(N_g,1)...
         - ( R_In_H - R_Out_H ) * N_0;
M_0    = F_plus ./ Kappa;
Kvec   = ( R_In_H - R_Out_H ) * ( N_0 - N_Eq );
Kvec   = F_plus ./ Kappa;

dt = [ 1.d-7 3.d-7 1.d-6 3.d-6 1.d-5 3.d-5 1.d-4 3.d-4 1.d-3 ]';

loglog( eC, N_0, '-k', 'linewidth', 2 )
hold on
for i = 1 : size(dt,1)
  loglog( eC, ( eye(N_g)-dt(i).*Pmat ) \ (eye(N_g)+dt(i).*Lmat)*N_0, '--' )
end
axis([1 300 1.d-10 10 ])
hold off
return

dt_K = min( 1.0 ./ abs( Kvec + 1.0d-16 ) );
dt   = min( [ 1.0 dt_K ] );

% Update Matrix from Nold to Nnew:
Wmat = ( eye(N_g) - dt .* Pmat )...
       \ ( eye(N_g) + dt .* ( Lmat - Pmat ) );

figure(1)
loglog(eC,N_0./abs(Lmat*N_0),'-x')
return
semilogx...
  ( eC, Lmat*N_0, '-o',...
    eC, Pmat*(N_0-N_Eq), '-x',...
    eC, Kvec, '-*' )

return

figure(2)
semilogx...
  ( eC, N_0,      '-o',...
    eC, Wmat*N_0, '-+' )

figure(3)
plot( eC, Wmat )

figure(4)
imagesc( Wmat )
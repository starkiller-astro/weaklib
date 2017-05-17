clear all
close all

Model = '004';

FileName = [ 'Relaxation_Analytic_dt_1m7_' Model '.mat' ];

N_g = 40; % Number of Energy Groups
[ eC, dV, R_In, R_Out, N_Eq ] = InitializeNES( Model, N_g );

% Multiply Rates with Momentum Space Volume Element:
R_In_H  = R_In  * diag( dV );
R_Out_H = R_Out * diag( dV );

N_t   = 10001;
Time = linspace(0.0, 1.0E-3, N_t)';
% Time  = logspace( -10, -3, N_t );

%
% Gaussian Initial Condition:
N_0 = 7.5d-1 .* exp( - ( 1.0d2 - eC ).^2 / 1.0d2 );

%
% Relaxation Matrix (Solving N' = A N):
A = L_FUN( N_0, R_In_H, R_Out_H, 0.0, N_g );

%
% Characteristic Decomposition:
[ R, D ] = eig( A );
lambda = diag( D );

%
% Initial Condition in Characteristic Variables:
W_0 = R \ N_0;

N_A = zeros( N_t, N_g );
W   = zeros( N_t, N_g );
tic
for i = 1 : N_t
  W(i,1:N_g)...
    = W_0(1:N_g) .* exp( lambda(1:N_g) .* Time(i) );
  N_A(i,1:N_g) = R * W(i,1:N_g)';
end
toc

%
% Save Solution:
Time_A    = Time;
Density_A = N_A;
save( FileName, 'Time_A', 'Density_A' );
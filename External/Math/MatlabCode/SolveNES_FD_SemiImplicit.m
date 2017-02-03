clear all
close all

Model = '001';

FileName = [ 'SemiImplicit_5m9_' Model '.mat' ];

N_g = 40; % Number of Energy Groups
[ eC, dV, R_In, R_Out ] = InitializeNES( Model, N_g );

% Multiply Rates with Momentum Space Volume Element:
R_In_H  = R_In  * diag( dV );
R_Out_H = R_Out * diag( dV );

%
% Computational Parameters:
t_end = 1.0d-4; % [ s ]
t     = 0.0d-0; % [ s ]
dt    = 5.0d-9; % [ s ] Time step
theta = 0.0;

%
% Gaussian Initial Condition:
N_0 = 7.5d-1 .* exp( - ( 1.0d2 - eC ).^2 / 1.0d2 );

%
Ntot_0 = sum( N_0 .* dV );

fig = figure;
figure( fig );
loglog( eC, N_0, '-ob', 'LineWidth', 2 )
axis( [ 1 3.0e2 1e-5 1.0e1 ] )
xlabel( '\epsilon [MeV]', 'FontSize', 20 )
ylabel( 'Particle Density [MeV^{-3}]', 'FontSize', 20 )
hold on

%
% Solve to End Time:
done  = false;
cycle = 0;
Nold  = N_0;
Nnew  = Nold;

tic
while ( not( done ) )

  cycle = cycle + 1;

  %
  % Update Matrix (Use Old State):
  F_In = R_In_H * Nold;
  k = R_Out_H * ones( N_g, 1 ) + theta .* ( R_In_H - R_Out_H ) * Nold;
  L = R_In_H - diag( k );
  
  if( t + dt > t_end )
    dt = t_end - t;
  end

  %
  % Record Detailed Solution:
  Time      (cycle,1) = t;
  Density   (cycle,:) = Nold;
  DensitySum(cycle,1) = sum( Nold .* dV );

  %
  % Backward Euler Solve:
  A = ( eye( N_g ) - dt * L );
  Nnew = A \ Nold;

  dN = max( abs( Nnew - Nold ) ./ Nold );

  Nold  = Nnew; 
  t     = t + dt;

  if ( t >= t_end )
    done = true;
  end

  if( mod(cycle, 100) == 1 )
    disp( fprintf( '  Cycle = %d, t = %d, dt = %d, ||dN|| = %d ', cycle, t, dt, dN ) );
  end

end
toc

figure( fig );
loglog( eC, Nnew, '-ok', 'LineWidth', 2 )
legend( 'Initial Condition', 'Final Solution',...
        'Location', 'SouthWest' )
hold off

Ntot = sum( Nnew .* dV );

disp( fprintf( 'Initial Particle Number = %d', Ntot_0 ) );
disp( fprintf( 'Final   Particle Number = %d', Ntot   ) );
disp( fprintf( '    Relative Difference = %d', abs( Ntot - Ntot_0 ) / Ntot_0 ) );

%
% Save Solution:
save( FileName, 'Time', 'Density', 'DensitySum' );
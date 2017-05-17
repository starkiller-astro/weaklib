clear all
close all

Model = '003';

Method = 'PenalizationEuler';

FileName = [ 'SemiImplicit_' Method '_dt_1m7_theta_0_' Model '.mat' ];

N_g = 40; % Number of Energy Groups

[ eC, dV, R_In, R_Out, N_Eq ] = InitializeNES( Model, N_g );

% Multiply Rates with Momentum Space Volume Element:
R_In_H  = R_In  * diag( dV );
R_Out_H = R_Out * diag( dV ); 

%
% Computational Parameters:
t_end  = 1.0d-4;  % [ s ]
t      = 0.0d-0;  % [ s ]
dt_max = 1.0d-7;  % [ s ] Maximum time step
dt     = 1.0d-7;  % [ s ] Initial time step
alpha  = 1.00;    % Factor to increase time step
theta  = 0.0;     % Include Fermi blocking with theta = 1.0
Tol_N  = 1.0d-14; % Tolerance for Newton iterations
ApplyCorrection = true; % For SIRK2 method


%
% Gaussian Initial Condition:
N_0 = 7.5d-1 .* exp( - ( 1.0d2 - eC ).^2 / 1.0d2 );
% N_0 = N_Eq;

%
% Initial Particle Number:
Ntot_0 = sum( N_0 .* dV );

%
% Solve to End Time:
done  = false;
cycle = 0;
Nold  = N_0;
Nnew  = Nold;

tic
while ( not( done ) )

  cycle = cycle + 1;
  
  dt = min( [ alpha * dt dt_max ] );
  
  if( t + dt > t_end )
    dt = t_end - t;
  end

  %
  % Record Detailed Solution:
  Time      (cycle,1) = t;
  TimeStep  (cycle,1) = dt;
  Density   (cycle,:) = Nold;
  DensitySum(cycle,1) = sum( Nold .* dV );

  switch Method
    case 'BackwardEuler'
      
      %
      % Backward Euler Solve:
      [ Nnew, Iterations(cycle,1) ]...
        = Update_ImplicitEuler...
            ( Nold, dt, R_In_H, R_Out_H, theta, N_g, Tol_N );
      
    case 'PenalizationEuler'
      
      %
      % Penalization Method Solve (Euler):
      [ Nnew, Iterations(cycle,1) ]...
        = Update_PenalizationMethod...
            ( Nold, N_Eq, dt, R_In_H, R_Out_H, theta, N_g );

    case 'PenalizationSIRK2'

      %
      % Penalization Method Solve (SIRK2):
      [ Nnew, Iterations(cycle,1) ]...
        = Update_PenalizationMethod2...
            ( Nold, N_Eq, dt, R_In_H, R_Out_H, theta, N_g,...
              ApplyCorrection );
          
    case 'SSPRK3'
        
      %
      % Third Order Explicit Runge-Kutta Method:
      
      [ Nnew, Iterations(cycle,1) ]...
        = Update_SSPRK3...
            ( Nold, dt, R_In_H, R_Out_H, theta, N_g );

    otherwise

      Nnew = Nold;
      Iterations(cycle,1) = 1;

  end

  dN = max( abs( Nnew - Nold ) ./ Nold );

  Nold = Nnew; 
  t    = t + dt;

  if ( t >= t_end )
    done = true;
  end

  if( mod(cycle, 100) == 1 )
    disp( fprintf( '  Cycle = %d, t = %d, dt = %d, ||dN|| = %d ', cycle, t, dt, dN ) );
  end

end
toc

Ntot = sum( Nnew .* dV );

disp( fprintf( 'Initial Particle Number = %d', Ntot_0 ) );
disp( fprintf( 'Final   Particle Number = %d', Ntot   ) );
disp( fprintf( '    Relative Difference = %d', abs( Ntot - Ntot_0 ) / Ntot_0 ) );

%
% Save Solution:
save( FileName, 'Time', 'TimeStep', 'Density', 'DensitySum', 'Iterations' );
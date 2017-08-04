function [ output_args ] = isNeq( N, R_In_H, R_Out_H, theta, N_g )
%Test whether the distribution function N is
%   the equilibrium number distribution function

output_args = false;

LN = L_FUN( N, R_In_H, R_Out_H, theta, N_g ) * N;

critical = ones( N_g, 1) * 1E-10;

logi =  ( abs(LN) < critical );

if ( isequal( logi, ones( N_g, 1) ) )
    output_args = true;
    disp( fprintf( 'Is N the equilibrium number density: YES' ) );
else
    disp( fprintf( 'Is N the equilibrium number density: NO!' ) );
    disp( fprintf( '      The maximum absolute value of L*N is   %d', max( abs( LN ) ) ) );    
end

end


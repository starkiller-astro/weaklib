function [ T ] = ComputeTempFromIntEnergy_Bisection...
    ( rho, E, Y, D1D, T1D, Y1D,  E_table, OS )
%ComputeTempFromIntEnergy_Bisection
%   E_table is internal energy density table

T = zeros( size( rho,1) );

for ii = 1 : size( rho, 1 )
    
    a = T1D(1)*1.001;
    E_a = interpolateEos( rho(ii), a, Y(ii), D1D, T1D, Y1D, E_table, OS );
    f_a = E(ii) - E_a;
    
    b = T1D(end-1);
    E_b = interpolateEos( rho(ii), b, Y(ii), D1D, T1D, Y1D, E_table, OS );
    f_b = E(ii) - E_b;
    
    if ( (f_a * f_b) > 0.0 )
        disp('Error: ComputeTempFromIntEnergy_Bisection');
        disp('Error: No Temperature found');
        break
    end
    
    ab = b - a;
    
    Converaged = false;
    Iter = 0;
    
    while ( not(Converaged) )
        Iter = Iter + 1;
        ab = 0.5 * ab;
        c = a + ab;
        E_c =  interpolateEos( rho(ii), c, Y(ii), D1D, T1D, Y1D, E_table, OS );
        f_c = E(ii) - E_c;
        
        if ( f_a*f_c < 0.0 )
            b = c;
            f_b = f_c;
        else
            a = c;
            f_a = f_c;
        end
        
        if( abs( f_c )/ E(ii) < 1.0E-10 ) 
            Converaged = true ;
        end
        
        if ( Iter > 128 )
            disp('ComputeTempFromIntEnergy_Bisection');
            disp('No Convergence After 128 iterations.');
        end
        
    end
    
    T(ii) = c;
    
end

end


function [ Nnew, nIter ] = Newton( Nold, dt, R_In, R_Out, theta, N_g, Tol )

  Nnew = Nold;
  
  converged = false;
  nIter     = 0;
  
  while ( not( converged ) )

    dN = - F_JAC( Nnew, dt, R_In, R_Out, theta, N_g )...
           \ F_FUN( Nnew, Nold, dt, R_In, R_Out, theta, N_g );
       
    Nnew = Nnew + dN;
    
    nIter = nIter + 1;
    
    if( norm( dN ./ Nnew ) < Tol )
      converged = true;
    end
    
  end

end


function [ dt ] = determine_dt_penalization( dtold, LL, N, dC, dt_max, diff_ratio )
% LL for L_FUN
% dC for P(Neq)
  
  C = LL*N;
  
  lamp = max(abs(eig(dC)));
  
  %% Boundary limit %%%%%%%%%%%%%%%%%%%%%%%%%%
  lim1sign = ( ones(size(C)) - sign( C ) )/ 2.0; % for L < 0
  lim1  = - N ./ C;
  lim1p = lim1 .* lim1sign; 
   
  lim2sign = ( ones(size(C)) + sign( C ) )/ 2.0; % for L > 0
  lim2  = ( ones( size(N) ) - N ) ./ C;
  lim2p = lim2 .* lim2sign; 
  
  lim = min( lim1p + lim2p );
  if ( lim < 1.0 / lamp)
      dt_0 = lim / ( 1.0 - lim * lamp );
  else
      dt_0 = dt_max;
  end
  
  dt = min( [dt_0 dtold*1.05 dt_max]);
  
    %% Consider the difference ratio between penalization method and forward Eular
    %  and normalized by N
  LLN = max( abs ( LL*N ./ N ) );
  
  dt_K5 = diff_ratio / (2.0 * LLN) + sqrt( diff_ratio / ( lamp * LLN ) + ...
      diff_ratio * diff_ratio / (4.0 * LLN * LLN ) );
  
  dt = min( [ dt dt_K5] );

end


function [ dt ] = determine_dt_forwardEular( dtold, LL, N, dt_max )
% LL for L_FUN
% dC for P(Neq)

  
  dt_lim1 =  ( ones( size(N) ) - N ) ./ ( LL * N ) / size(N,1) ; %  LL * N > 0
  dt_lim1sign = ( ones( size(N) ) + sign(  LL * N ) )/ 2.0;
  dt_lim1p = dt_lim1 .* dt_lim1sign;
  
  dt_lim2 = ( ( - N ) ./ ( LL * N ) ) / size(N,1); % LL * N < 0
  dt_lim2sign = ( ones( size(N) ) - sign(  LL * N) )/ 2.0;
  dt_lim2p = dt_lim2 .* dt_lim2sign;
  
  dt_0 = min( dt_lim1p + dt_lim2p ) ;

  dt = min( [ 1.05*dtold dt_0  dt_max]);
 
end


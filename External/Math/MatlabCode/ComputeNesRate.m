function [ R_0 ] = ComputeNesRate...
    ( Ep, E, D, T, Y, E1D, D1D, T1D, Y1D,...
    Eta1D, Op, OpOS, ecmpTable, ecmpOS )
%   Ep and E are array
%   D, T, T are number
   
%
%   Interpolating electron chemical potential

ecmp = interpolateEos( D, T, Y, D1D, T1D, Y1D, ecmpTable, ecmpOS );

%   Compute eta value
kmev = 8.61733d-11;
eta  = ecmp/T/kmev;

%   Interpolating
N_g = size(Ep,1);
R_0 = zeros(N_g, N_g);

for ii = 1:N_g
    R_0(:,ii) = interpolate4D( ...
        Ep, E(ii)*ones(N_g,1), T*ones(N_g,1), eta*ones(N_g,1),...
        E1D, E1D, T1D, Eta1D, [1,1,1,1], Op, OpOS );
end

% Unit Convertion cm-1 to s-1 and integral over mu'
c  = 2.99792d+10;  % cm s-1
R_0 = R_0 * c * 4 * pi; % Should be this one

end


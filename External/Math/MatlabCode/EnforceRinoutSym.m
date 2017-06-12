function [ R_In_low, R_Out_low ] = EnforceRinoutSym( R_out, eC, T_K )
% This function test the symmetry of R_In and R_Out

R_Out_low = tril(R_out);

kmev = 8.61733d-11;

for ii = 1:size(R_out,1)
    for jj = ii:size(R_out,1)
        beta = (eC(ii)-eC(jj))/(kmev*T_K);
        R_Out_low(ii,jj) = R_Out_low(jj,ii)*exp(beta);
    end
end

R_In_low = R_Out_low';

end


function [ R_In_H, R_Out_H ] = Update_NESKernel...
    ( eC, dV, rho, T, Y, E1D, D1D, T1D, Y1D, ...
    Eta1D, Op, OpOS, ecmpTable, ecmpOS)


R = ComputeNesRate...
    (eC, eC, rho, T, Y, E1D, D1D, T1D, Y1D,...
    Eta1D, Op, OpOS, ecmpTable, ecmpOS );

R_In  = R;
R_Out = R';
R_In_H  = R_In  * diag( dV ) ;
R_Out_H = R_Out * diag( dV ) ;

end


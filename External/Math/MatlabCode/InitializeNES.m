function [ eC, dV, R_In, R_Out, N_Eq ] = InitializeNES( Model, N_g )

  %
  %  Energy Bin Centers:
  [ eC ] = ReadData1D( '../Data/NES_RATES_EnergyBinCenter.dat', N_g );

  %
  % Energy Bin Widths:
  [ de ] = ReadData1D( '../Data/NES_RATES_EnergyBinWidth.dat', N_g );

  %
  % Energy Bin Volume:
  dV = ( (eC+0.5d0*de).^3 - (eC-0.5d0*de).^3 )/3.0d0;
  
  %
  % NES in Rates:
  [ R ] = ReadData2D( ['../Data/NES_RATES_R_In___' Model '.dat'], N_g, N_g );
  
  R_In  = R;
  R_Out = R';
  
  switch Model
    case '001'
      N_Eq = 1.0 ./ ( exp( (eC-145.254)./20.5399 ) + 1.0 );
    case '002'
      N_Eq = 1.0 ./ ( exp( (eC-045.835)./15.9751 ) + 1.0 );
    case '003'
      N_Eq = 1.0 ./ ( exp( (eC-020.183)./07.7141 ) + 1.0 );
    case '004'
      N_Eq = 1.0 ./ ( exp( (eC-009.118)./07.5830 ) + 1.0 );
    case '005'
      N_Eq = 1.0 ./ ( exp( (eC-003.886)./03.1448 ) + 1.0 );
    otherwise
      N_Eq = 1.0;
  end

end


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
  
  switch Model
    case '001'
      kT = 20.5399;
      N_Eq = 1.0 ./ ( exp( (eC-145.254) ./ kT ) + 1.0 );
    case '002'
      kT = 15.9751;
      N_Eq = 1.0 ./ ( exp( (eC-045.835) ./ kT ) + 1.0 );
    case '003'
      kT = 07.7141;
      N_Eq = 1.0 ./ ( exp( (eC-020.183) ./ kT ) + 1.0 );
    case '004'
      kT = 07.5830;
      N_Eq = 1.0 ./ ( exp( (eC-009.118) ./ kT ) + 1.0 );
    case '005'
      kT = 03.1448;
      N_Eq = 1.0 ./ ( exp( (eC-003.886) ./ kT ) + 1.0 );
    otherwise
      N_Eq = 1.0;
  end
  
  R_In  = R;
  for j = 1 : N_g
    for i = 1 : N_g
      if( i < j )
        R_In(i,j) = R_In(j,i)*exp( (eC(j)-eC(i))/kT );
      end
    end
  end
  R_Out = R_In';
end


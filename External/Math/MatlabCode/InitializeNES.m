function [ eC, dV, R_In, R_Out ] = InitializeNES( Model, N_g )

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

end


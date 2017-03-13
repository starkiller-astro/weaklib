function [ E, D, T, Y, R_0, R_1, OS ] = readIsoScatteringOpacityTable( opacityTableName )

  % Reads HDF5 Opacity Table

  disp( fprintf( 'INFO: Reading Isoenergetic Opacity from file: %s', opacityTableName ) );

  % Independent Variables:
  E = h5read( opacityTableName, '/EnergyGrid/Values' );
  D = h5read( opacityTableName, '/ThermoState/Density' );
  T = h5read( opacityTableName, '/ThermoState/Temperature' );
  Y = h5read( opacityTableName, '/ThermoState/Electron Fraction' );
  
  OS = h5read( opacityTableName, '/scatt_Iso/Offsets' );
  
  Tmp = h5read( opacityTableName, '/scatt_Iso/Kernel/Electron Iso-sca Kernel Momen' );
  R_0 = 10.^( Tmp(:,:,:,:,1) ) -OS(1);
  R_1 = 10.^( Tmp(:,:,:,:,2) ) -OS(2);

end


function [ E, T, Eta, R_0, R_1, OS] = readNesScatteringOpacityTable( opacityTableName )

  % Reads HDF5 Opacity Table

  disp( fprintf( 'INFO: Reading NES Scattering Opacity from file: %s', opacityTableName ) );

  % Independent Variables:
  E   = h5read( opacityTableName, '/EnergyGrid/Values' );
  T   = h5read( opacityTableName, '/ThermoState/Temperature' );
  Eta = h5read( opacityTableName, '/EtaGrid/Values' );
  
  OS = h5read( opacityTableName, '/scatt_NES/Offsets' );
  
  Tmp = h5read( opacityTableName, '/scatt_NES/Kernel/Electron non-Iso Kernel Momen' );
  R_0 = 10.^( Tmp(:,:,:,:,1) ) - OS(1);
  R_1 = 10.^( Tmp(:,:,:,:,2) ) - OS(2);
  
end


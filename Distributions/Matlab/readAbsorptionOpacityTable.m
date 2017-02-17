function [ E, D, T, Y, Chi ] = readAbsorptionOpacityTable( opacityTableName )

  % Reads HDF5 Opacity Table

  disp( fprintf( 'INFO: Reading Absorption Opacity from file: %s', opacityTableName ) );

  % Independent Variables:
  E = h5read( opacityTableName, '/EnergyGrid/Values' );
  D = h5read( opacityTableName, '/ThermoState/Density' );
  T = h5read( opacityTableName, '/ThermoState/Temperature' );
  Y = h5read( opacityTableName, '/ThermoState/Electron Fraction' );
  
  Chi = h5read( opacityTableName, '/thermEmAb/Absorptivity/Electron Capture Absorptivity' );

end


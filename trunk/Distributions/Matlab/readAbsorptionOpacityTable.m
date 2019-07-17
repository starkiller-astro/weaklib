function [ E, D, T, Y, Chi, OS ] = readAbsorptionOpacityTable( opacityTableName )

  % Reads HDF5 Opacity Table

  disp( fprintf( 'INFO: Reading Absorption Opacity from file: %s', opacityTableName ) );

  % Independent Variables:
  E = h5read( opacityTableName, '/EnergyGrid/Values' );
  D = h5read( opacityTableName, '/ThermoState/Density' );
  T = h5read( opacityTableName, '/ThermoState/Temperature' );
  Y = h5read( opacityTableName, '/ThermoState/Electron Fraction' );
  
  OS = h5read( opacityTableName, '/thermEmAb/Offsets' );
  
  Tmp = h5read( opacityTableName, '/thermEmAb/Absorptivity/Electron Capture Absorptivity' );
  Chi = 10.^(Tmp) - OS;

end


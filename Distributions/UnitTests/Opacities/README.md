# WeakLib Opacity UnitTests

## Implement Examples

- WeakLib Opacity table IO test:
    
    - `wlReadOpacityTableTest.f90`

- Interpolating on WeakLib Opacity tables:

    - `wlInterpolateAb.f90` for absorptivity and emissivity table
    - `wlInterpolateIso.f90` for isoenergetic scattering table
    - `wlInterpolateNES.f90` for neutrino-electron scattering table
    - `wlInterpolatePair.f90` for pair process table


## Step-by-step instruction for obtaining the comprehensive opacities with a given profile:

1.  Write the profile into the stardard format, as `ProfileExample.profile`.

    Keys:
    - Filename ends with `.profile`
    - First line of the profile gives the number of rho-T-Ye states
    - Star the data from the third line
    - Given the states in radius-density-temperature-electron fraction order
    - in unit of [Arbitrary], [g/cm^3], [K], [Dimensionless]

2.  Download, link or place the WeakLib EoS table and opacity tables with the profile under a same directory.
    
    Keys:
    - WeakLib tables can be found at https://code.ornl.gov/astro/weaklib-tables.
    - You can create a new directory under `./Executables` or any place you want to run the test.
    
3.  Copy `WriteTxtForInterpolation.sh` to the same directory and run it by
    ```$./WriteTxtForInterpolation.sh```
    
    A `dataList.txt` file should be generated.
    
    Keys:
    - Check `dataList.txt` to see if the files' name are all listed as profile, Eos table, opacity table
    
4.  Next, generate the interpolating executables.

    Keys:
    - Set your compile environment by 
    
    ```$source weaklib/Distributions/Workflow/SetEnvironment.sh machinename```
    
      Replace `weaklib` directory and `machinename` as needed
      
    - Then go to `Distributions/UnitTests/Opacities/Executables` make using
    
    ```$make Interpolations```
    
    to compile.
    - If everything works correctly, you should see four executables `wlInterpolate*_machinename` now.


5.  Link the executables to the same directory as the tables, and run the executables:
    
    - Run the executables by
    
    ```$./wlInterpolateAb_machinename```
    
    ```$./wlInterpolateIso_machinename```
    
    ```$./wlInterpolateNES_machinename```
    
    ```$./wlInterpolatePair_machinename```
    
    - Each of the run gives a `Interpolated*Output.h5` file which is the opacity rates.
    
6.  To plot the `Interpolated*Output.h5` result, use the `plotInterpolatedCompactHDF.m` matlab script.

    Keys:
    - Change the file directory and axis limit as needed.


## Ask For Help
- R. Chu : rchu@vols.utk.edu

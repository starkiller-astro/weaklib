# WeakLib Opacity UnitTests

## Implement Examples

- WeakLib Opacity table IO:

     wlReadOpacityTableTest.f90

- Interpolating on WeakLib Opacity tables:

     wlInterpolateAb.f90
     wlInterpolateIso.f90
     wlInterpolateNES.f90
     wlInterpolatePair.f90

- Step-by-step instruction for obtaining comprehensive opacities with a given profile:

1.  Write the profile into the stardard format, as `ProfileExample.profile`.

    Keys: 
         - Filename ends with `.profile`
         - First line of the profile gives the number of rho-T-Ye states
         - Star the data from the third line
         - Given the states in radius-density-temperature-electron fraction order
         - in unit of [Arbitrary], [g/cm^3], [K], [Dimensionless]

2.  d    



## Ask For Help
- R. Chu : rchu@vols.utk.edu

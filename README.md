# Welcome to WeakLib

WeakLib is a code library for astrophysical applications. It provides preprocessed equation of state (EoS) and neutrino opacity tables for use in neutrino transport calculations. WeakLib tables are intended to be usable in a straightforward manner by GENASIS, CHIMERA, FLASH, thornado, and other radiation hydrodynamics codes.



### Available Tables Download
[https://code.ornl.gov/astro/weaklib-tables](https://code.ornl.gov/astro/weaklib-tables)



### Build Your Own Tables
- #### Build EoS Table

- #### Build Opacity Table

  * Step 1: Have the EOS table ready.

     * Make sure it covers the domain and the physics with the resolution you desired

     * If the Eos table has a out-of-date structure or you want a different domain or 
       resultion, you can try to use
       weaklib/Distributions/UnitTests/EquationOfState/wlRewriteEquationOfStateTest.f90
       to rewrite the Eos table.
       Otherwise, build a new Eos table.

  * Step 2: Edit the opacity table creating driver:
    weaklib/External/Utilities/Opacities/Bruenn85/wlCreateOpacityTable.f90
    and compile the excutable.

  * Step 3: Move the excutable and Eos table to same directory and run the excutable.

  * Step 4: Test whether the opacity table is created correctly by runing interpolation test.
    Examples can be found under weaklib/Distributions/UnitTests/Opacities.


### Interpolate in WeakLib Opacity Tables

[Distributions/UnitTests/Opacities](Distributions/UnitTests/Opacities)


## Ask For Help
- R. Chu : rchu@vols.utk.edu

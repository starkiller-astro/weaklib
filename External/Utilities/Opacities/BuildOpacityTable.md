## Build Opacity Table

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
    Examples can be found under [/Distributions/UnitTests/Opacities](/Distributions/UnitTests/Opacities).

### Copy-To-Go Bash Script
  You can copy the following code and save it as '/weaklib/build.sh' and do '$source build.sh'.
  **Note**: change machine setting, 'summit', to the one you're running on. (Check [Distributions/Workflow/SetEnvironment.sh](Distributions/Workflow/SetEnvironment.sh) for the options)

```
cd External/Utilities/Opacities/Bruenn85/Executables
export WEAKLIB_DIR='../../../../../../weaklib'
source ../../../../../Distributions/Workflow/SetEnvironment.sh summit
[ -e TableCreatingLog ] && rm TableCreatingLog
make clobber
make wlCreateOpacityTable
[ -e wl-EOS-SFHo-15-25-50.h5 ] && rm wl-*.h5
ln -s ~/weaklib-tables/SFHo/LowRes/wl-EOS-SFHo-15-25-50.h5 .
./wlCreateOpacityTable_summit > TableCreatingLog
```

### Ask For Help
- R. Chu : rchu@vols.utk.edu

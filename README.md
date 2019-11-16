## Build Opacity Table

There are two ways to build a WeakLib opacity table.

(1) Build the table from the source code:

  * Have the EOS table ready.

     * Make sure it covers the domain and the physics with the resolution you desired

     * If the Eos table has a out-of-date structure or you want a different domain or 
       resultion, you can try to use
       weaklib/Distributions/UnitTests/EquationOfState/wlRewriteEquationOfStateTest.f90
       to rewrite the Eos table.
       Otherwise, build a new Eos table.

  * Edit the opacity table creating driver:
    weaklib/External/Utilities/Opacities/Bruenn85/wlCreateOpacityTable.f90
    and compile the excutable.

  * Move the excutable and Eos table to same directory and run the excutable.

  * Test whether the opacity table is created correctly by runing interpolation test.
    Examples can be found under weaklib/Distributions/UnitTests/Opacities.

(2) Send an email to rchu@vols.utk.edu with your requires.

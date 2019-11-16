## Build Opacity Table

There are two ways to build a WeakLib opacity table.

(1) Build the table from the source code:

  * Have the EOS table ready.

     * Make sure it covers the domain and the physics with the resolution you desired

     * If the Eos table has a out-of-date structure or you want a different domain or 
       resultion, you can try
       `weaklib/Distributions/UnitTests/EquationOfState/wlRewriteEquationOfStateTest.f90`
       to rewrite the Eos table.
       Otherwise, build a Eos table.

  * Edit the opacity table creating driver. 

  * Link the excuted driver and Eos table under same directory and run the excutable.

  * Test whether the opacity table is created correctly by runing interpolation test.

(2) Send an email to rchu@vols.utk.edu with your requires.

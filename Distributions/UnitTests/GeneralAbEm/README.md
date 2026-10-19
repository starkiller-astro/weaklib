# GeneralAbEm unit tests

Tests for the semi-leptonic charged-current absorption/emission kernels in
`Distributions/OpacitySource/OpacitySourceGSI`:

| Short name | Routine | Notes |
|---|---|---|
| W | `wlSemiLeptonicOpacityModule2D :: Opacity_CC_2D` (and `_PrecomputedNodes`) | 2D Gauss-Legendre, no CUBA |
| I | `wlSemiLeptonicOpacityIntegrationModule2D` wrappers | two-phase (GPU-oriented) driver, should match W bit for bit |
| 4D | `wlSemiLeptonicOpacityModule4D :: Opacity_CC_4D` | Vegas (CUBA) |
| G | `Opacity_SemiLep_2D_GSI.f90`, `Opacity_SemiLep_4D_GSI.f90` | original GSI reference code (unchanged physics, serial only) |

Reactions: 1 `nu + n -> l- + p`, 2 `nubar + p -> l+ + n`, 3 `nubar + l- + p -> n`
(inverse neutron decay), 4 `n -> p + l- + nubar` (emissivity).
Corrections: 0 LO, 1 +weak magnetism, 2 +pseudoscalar, 3 +form factors.

Units: the weaklib kernels return cm^-1 (reactions 1-3) and (s cm^3 MeV)^-1
(reaction 4). The GSI routines return km^-1 for reactions 1-3 (the tests
convert them) and (s cm^3 MeV)^-1 for reaction 4.

## Tests

| Test | Target | What it does |
|---|---|---|
| `wlTestAbEm2DVsGSI` | `all` | Accuracy test. For up to `nPointsMax` (default 100) points of `ThermoConditions.dat`, 10 log-spaced energies (2-300 MeV), reactions 1-4 and corrections 0-3 it computes W, W with precomputed nodes, I and G with nE = 50 (the GSI code has `Ngrids = 50` hardcoded). It prints per (reaction, correction) the max/median relative difference W-P, W-I (with NaN counts), W-G (all samples and "significant" samples, i.e. >= 1e-6 of the maximum over energies at that point), and for reactions 3,4 W-\|G\| and the number of sign mismatches. It also checks detailed balance between reactions 3 and 4, `j = c kappa_a E^2/(2 pi^2 (hbar c)^3) exp(-(E-mu_nubar)/T)` with `mu_nubar = Mun - Mup - Mul`, and lists the worst samples. Per-sample values go to `AbEm2DVsGSI.dat`. |
| `wlTestAbEm2DParallel` | `all` | Timing and consistency of W in serial, OpenMP, through the I wrappers and with the collapsed two-phase loop. Prints any mismatch above 1e-8 (relative). `nPointsMax` defaults to 200 (all 1000 points take about 5 min on 4 threads). |
| `wlTestAbEm2DVsTable` | `all` | Compares W against the tabulated EmAb opacities (needs the 3D SFHo EOS table and EmAb opacity table, so build with `EOSMODE=3D`). Also times both approaches with OpenMP: table interpolation (EOS chemical potentials + EmAb) vs the 2D calculation (W and the integration module), reporting wall time, time per (point, energy, species) and the 2D/table ratio. |
| `wlTestAbEm4DVsGSI` | `cuba` | W and 4D compared with the GSI 2D/4D codes. Reports sign mismatches, does not hide them. |
| `wlTestAbEm4DParallel` | `cuba` | Timing of the 4D kernels, serial vs OpenMP. |
| `wlCompareAbEmApproximations` | `compare` | Along `ThermoConditions.dat` (T, Rho, Ye only; the rest is recomputed from the EOS), compares W (and 4D when built with CUBA) with the elastic approximation (`wlCalculateAbEmOpacityModule`) and the EmAb table. Writes `.bin` files. |
| `wlCompareAbEmApproximationsTable` | `compare` | The same comparison on a sub-grid of the EOS/opacity table. Also compares the Horowitz weak-magnetism factor with Chimera's `cc_weak_mag_weaklib` (from `TableCreator`). OpenMP parallel. |

## Building

```
export WEAKLIB_DIR=/path/to/weaklib
export WEAKLIB_MACHINE=ubuntu_gnu        # or source Workflow/SetEnvironment.sh <machine>
cd Executables

make                                      # all: 2D tests, no CUBA, EOSMODE=3D
make cuba USE_CUBA=TRUE                   # 4D tests, need CUBA_DIR (dir with libcuba.a)
make compare EOSMODE=COMPOSE              # or EOSMODE=3D; add USE_CUBA=TRUE for 4D
make clobber
```

- **`make clobber` when you switch `EOSMODE` or `USE_CUBA`.** The library
  objects are built with `-DEOSMODE_x` / `-DWEAKLIB_CUBA`, and make will not
  rebuild them by itself.
- `$(CUBA_LIB)` is linked only when `USE_CUBA=TRUE`. `make cuba` without it
  stops with an error.
- The `compare` targets use `TableCreator/OpacitySource/FullWeakPhysics/Chimera`
  (through VPATH) for `cc_weak_mag_weaklib`.
- EOSMODE for `compare`: `3D` reads `wl-EOS-SFHo-15-25-50.h5`, using bare
  nucleon masses, U = 0 and relativistic chemical potentials rebuilt as in the
  table creator (`Mun = cmpn + dmnp + mn`, `Mup = cmpp + dmnp + mp`).
  `COMPOSE` reads `BaryonsPlusPhotonsPlusLeptonsEOS.h5` together with its
  Helmholtz electron and muon tables.

## Inputs (symlinks in `Executables/`)

| File | Needed by |
|---|---|
| `ThermoConditions.dat` (-> `../ThermoConditions.dat`) | all tests except `wlCompareAbEmApproximationsTable` |
| `wl-EOS-SFHo-15-25-50.h5` | `wlTestAbEm2DVsTable`, `compare` with `EOSMODE=3D` |
| `wl-Op-SFHo-15-25-50-E40-EmAb.h5` | `wlTestAbEm2DVsTable`, `compare` with `UseReddy = .FALSE.` |
| `wl-Op-SFHo-15-25-50-E40-EmAb-Reddy.h5` | `compare` (default `UseReddy = .TRUE.`) |
| `BaryonsPlusPhotonsPlusLeptonsEOS.h5` | `compare` with `EOSMODE=COMPOSE` |

`ThermoConditions.dat` holds 1000 rows. The first line is the row count and the
second the column names: `T Rho Ye Ym Mue Mum Mun Mup Un Up M*n M*p` (MeV,
g/cm^3, relativistic chemical potentials including rest mass and mean field,
mean-field potentials and effective masses).

## Running the CUBA tests

Run CUBA serially inside the OpenMP loops:

```
export CUBACORES=0
export CUBACORESMAX=0
OMP_NUM_THREADS=4 ./wlTestAbEm4DParallel_$WEAKLIB_MACHINE
```

Known issues in the 4D library module (see `OpacitySource/OpacitySourceGSI/README`):
- Its threadprivate block is guarded by `USE_OMP` rather than `WEAKLIB_OMP`.
- Its reaction-4 emissivity carries an extra factor 1e-5.

## Known differences with the GSI code

- G uses different constants (G_F, massA, Delta_np, ...), `Tfac = 150` and
  real(4) literals. Expect W-G at the 1e-3 level.
- For reactions 3,4, G clips |M|^2 >= 0 *before* applying the crossing minus
  sign and then flips the sign of the result. G is therefore <= 0 and mostly
  0 (it is nonzero only through the pseudoscalar and form-factor terms added
  after the clip). The sign-mismatch counts in `wlTestAbEm2DVsGSI` and
  `wlTestAbEm4DVsGSI` show this.

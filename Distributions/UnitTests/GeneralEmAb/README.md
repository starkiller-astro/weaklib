# GeneralEmAb unit tests

Tests for the semi-leptonic charged-current absorption/emission kernels in
`Distributions/OpacitySource/GeneralEmAb` (see the README there for the
layout and the three 2D integration methods ENERGY = "standard", MOMENTUM =
default, TAIL):

| Short name | Routine | Notes |
|---|---|---|
| W | `wlSemiLeptonicOpacityModule2D :: Opacity_CC_2D` after `InitOpacity_CC_2D(nE)` | 2D, default integration method (MOMENTUM unless built with `EMAB2D_METHOD=ENERGY\|TAIL`), no CUBA |
| I | `wlSemiLeptonicOpacityIntegrationModule2D` wrappers | two-phase (GPU-oriented) driver of the MOMENTUM method, should match W bit for bit (with the default method) |
| 4D | `wlSemiLeptonicOpacityModule4D :: Opacity_CC_4D` | Vegas (CUBA) |
| G | `GeneralEmAb/OriginalGSIroutines/Opacity_SemiLep_2D_GSI.f90`, `Opacity_SemiLep_4D_GSI.f90` | original GSI reference code (unchanged, serial only) |

Reactions: 1 `nu + n -> l- + p`, 2 `nubar + p -> l+ + n`, 3 `nubar + l- + p -> n`
(inverse neutron decay), 4 `n -> p + l- + nubar` (emissivity).
Corrections: 0 LO, 1 +weak magnetism, 2 +pseudoscalar, 3 +form factors.

Units: the weaklib kernels return cm^-1 (reactions 1-3) and (s cm^3 MeV)^-1
(reaction 4). The GSI routines return km^-1 for reactions 1-3 (the tests
convert them) and (s cm^3 MeV)^-1 for reaction 4.

## Tests

| Test | Target | What it does |
|---|---|---|
| `wlTestEmAb2DVsGSI` | `all` | Accuracy test. For up to `nPointsMax` (default 100) points of `ThermoConditions.dat`, 10 log-spaced energies (2-300 MeV), reactions 1-4 and corrections 0-3 it computes W (`Opacity_CC_2D`), P (the stateless `_Nodes` kernel), I and G with nE = 50 (the GSI code has `Ngrids = 50` hardcoded). It prints per (reaction, correction) the max/median relative difference W-P, W-I (with NaN counts), W-G (all samples and "significant" samples, i.e. >= 1e-6 of the maximum over energies at that point), and for reactions 3,4 W-\|G\| and the number of sign mismatches. It also checks detailed balance between reactions 3 and 4, `j = c kappa_a E^2/(2 pi^2 (hbar c)^3) exp(-(E-mu_nubar)/T)` with `mu_nubar = Mun - Mup - Mul`, and lists the worst samples. Per-sample values go to `EmAb2DVsGSI.dat`. |
| `wlTestEmAb2DParallel` | `all` | Timing and consistency of the MOMENTUM method (`InitOpacity_CC_2D(nE, iEmAb2D_Momentum)` + `Opacity_CC_2D`, so the check holds for any default) in serial, OpenMP, through the I wrappers and with the collapsed two-phase loop, for electrons and muons. Compares each against serial `Opacity_CC_2D`, prints any mismatch above 1e-8 (relative, `|a-b|/max(|a|,|b|)`) and the number of bit-identical samples (OpenMP and wrappers are expected to be bit-identical, the collapsed SUM only to round-off). `nPointsMax` defaults to 200 (about 2.5 min on 8 threads for both leptons). |
| `wlTestEmAb2DConvergence` | `all` | Accuracy vs the resolution parameter `nE` of one integration method, `iMethod` (default: the build's default method `iEmAb2D_Default`; set it to `iEmAb2D_Energy`/`_Momentum`/`_Tail` to test another). nE lists and benchmark are per method (MOMENTUM 6-96 vs 256, ENERGY 8-128 vs 256, TAIL 4-16 vs 24). For MOMENTUM nE is the number of outer nodes per segment in nucleon momentum, with `MAX(8,MIN(12,nE/2))` inner lepton nodes per segment (segments split at kinks and Fermi edges). For every `PointStride`-th row of `ThermoConditions.dat` (default all 1000), 12 log-spaced energies (1-300 MeV), electrons and muons, reactions 1-4 and corrections 0-3 it computes `Opacity_CC_2D` for `nE_List` (6-96) and compares with the benchmark `nE = 256`. Reports max/p99/median relative error (corrections 0-2 per reaction and lepton, correction 3 separately), zeros/NaNs, time per evaluation and speedup vs `nE = 50`, the smallest `nE` meeting 1e-2...1e-4 (by max and by p99) and the worst samples. Only significant samples count (ref >= 1e-6 of the max over energies and >= 1e-15 of the channel max). This is self-convergence: the floor of a few 1e-5 from the fixed inner rule and Tfac is not visible. Per-sample errors go to `EmAb2DConvergence.dat` (~35 MB). About 2.5 min on 8 threads. |
| `wlTestEmAb2DCompareSchemes` | `all` | Error vs cost of the three 2D integration methods, each selected with `InitOpacity_CC_2D(nE, method)` and evaluated with `Opacity_CC_2D`: A ENERGY (standard, nE x nE linear Gauss-Legendre in E2/E3), B MOMENTUM (`iEmAb2D_Momentum`: nE outer p2 nodes, kink-split inner rule) and C TAIL (`iEmAb2D_Tail`: no thermal cutoff, semi-infinite tail rule). Same samples and significance mask as `wlTestEmAb2DConvergence`. Each scheme runs its resolution list (`nE_A`, `nE_B`, `nE_C`) plus a high-resolution run `nE_Hi` (A 128, B 256, C 24; the cost of A and C grows like nE^2). The reference is C at 24 (`iRefScheme`): C has no Tfac and agrees with an independent 1024x1024 reference to ~1e-6 for corrections 0-2; B has a floor of a few 1e-5 from its inner rule (at most 12 nodes per segment). A "reference agreement" table gives the pairwise differences between the three high-resolution runs, i.e. how well the reference itself is known. Correction 3 does not converge in any scheme (~1e-2). For every run it prints max/p99/median error (corrections 0-2), max per reaction and per lepton, correction 3 max/median, zeros/NaNs, wall time and us/eval. It then lists all runs sorted by cost, and for 1e-2...1e-4 the cheapest (scheme, nE) meeting the tolerance by max and by p99, with the speedup vs B at nE = 50. `EmAb2DCompareSchemes.dat` (~60 MB) is written only if `WriteSamples`. About 12-13 min on 8 threads (A128, B256 and C24 take about half of it). |
| `wlTestEmAb2DFastSettings` | `all` | Cheapest quadrature setting meeting a physically weighted accuracy metric. Same samples as `wlTestEmAb2DConvergence` (`PointStride` default 1). Settings: B `Opacity_CC_2D_Momentum_Nodes` (nE2 in 6-24 outer x nE3 in 4-12 inner nodes per segment), A `Opacity_CC_2D_Energy_Nodes` (all nE2 x nE3 pairs of 8-24), C `Opacity_CC_2D_Tail_Nodes` (nE 4, 6, 8); reference C at nE = 16, plus a check run C24 (`DoCheck`). Metric: per (point, lepton, energy, correction) the reference species totals kappa_tot(nu_l) = kappa_1, kappa_tot(nubar_l) = kappa_2 + kappa_3; samples with kappa_tot,ref < `KappaFloor` (1e-10 cm^-1) are ignored; weighted error |kappa - ref|/kappa_tot,ref for reactions 1-3 and |j - j_ref|/j_ref * kappa_3,ref/kappa_tot,ref for reaction 4; tolerance `TolLow` (1e-2) for E <= E_split = max(mu,0) + `nT_Split` T (8, mu_nu = Mul + Mup - Mun, -mu_nu for nubar), `TolHigh` (0.2) above; NaN or a significant zero fails. A setting passes if the max normalised error (weighted error/tolerance) over all reactions, leptons and corrections 0-3 is <= 1. Prints per setting us/eval, max normalised error at low/high E, % failing, NaN/zero, the same on corrections 0-2 only, the strict max relative error of `wlTestEmAb2DConvergence` (corrections 0-2), the worst sample and PASS/FAIL (0-3 and 0-2); all runs sorted by cost with the cheapest passing setting; a sensitivity table (`TolLow` 3e-3/1e-2/3e-2 x `nT_Split` 5/8/12) from the stored values; and the worst samples of the cheapest passing setting. Correction 3 (clipped form factors) is not converged even in the reference (C24 vs C16 ~0.5 normalised at 300 MeV, vs ~5e-5 on corrections 0-2), which is why the 0-2 verdict is also given. All values are stored (~180 MB). `EmAb2DFastSettings.dat` only if `WriteSamples`. About 8-16 min on 8 threads (C16 and C24 take about half). |
| `wlTestEmAb2DVsTable` | `all` | Compares W against the tabulated EmAb opacities (needs the 3D SFHo EOS table and EmAb opacity table, so build with `EOSMODE=3D`). Also times both approaches with OpenMP: table interpolation (EOS chemical potentials + EmAb) vs the 2D calculation (W and the integration module), reporting wall time, time per (point, energy, species) and the 2D/table ratio. |
| `wlTestEmAb4DVsGSI` | `cuba` | W and 4D compared with the GSI 2D/4D codes. Reports sign mismatches, does not hide them. |
| `wlTestEmAb4DParallel` | `cuba` | Timing of the 4D kernels, serial vs OpenMP. |
| `wlCompareEmAbApproximations` | `compare` | Along `ThermoConditions.dat` (T, Rho, Ye only; the rest is recomputed from the EOS), compares W (and 4D when built with CUBA) with the elastic approximation (`wlCalculateEmAbOpacityModule`) and the EmAb table. Writes `.bin` files. |
| `wlCompareEmAbApproximationsTable` | `compare` | The same comparison on a sub-grid of the EOS/opacity table. Also compares the Horowitz weak-magnetism factor with Chimera's `cc_weak_mag_weaklib` (from `TableCreator`). OpenMP parallel. |

## Building

```
export WEAKLIB_DIR=/path/to/weaklib
export WEAKLIB_MACHINE=ubuntu_gnu        # or source Workflow/SetEnvironment.sh <machine>
cd Executables

make                                      # all: 2D tests, no CUBA, EOSMODE=3D
make cuba USE_CUBA=TRUE                   # 4D tests, need CUBA_DIR (dir with libcuba.a)
make EMAB2D_METHOD=ENERGY                 # default method of Opacity_CC_2D: MOMENTUM (default) | ENERGY | TAIL (make clobber first)
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
| `ThermoConditions.dat` (-> `../ThermoConditions.dat`) | all tests except `wlCompareEmAbApproximationsTable` |
| `wl-EOS-SFHo-15-25-50.h5` | `wlTestEmAb2DVsTable`, `compare` with `EOSMODE=3D` |
| `wl-Op-SFHo-15-25-50-E40-EmAb.h5` | `wlTestEmAb2DVsTable`, `compare` with `UseReddy = .FALSE.` |
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
OMP_NUM_THREADS=4 ./wlTestEmAb4DParallel_$WEAKLIB_MACHINE
```

Known issues in the 4D library module (see `OpacitySource/GeneralEmAb/README`):
- Its threadprivate block is guarded by `USE_OMP` rather than `WEAKLIB_OMP`.
- Its reaction-4 emissivity carries an extra factor 1e-5.

## Known differences with the GSI code

- G uses different constants (G_F, massA, Delta_np, ...), `Tfac = 150` and
  real(4) literals. Expect W-G at the 1e-3 level.
- For reactions 3,4, G clips |M|^2 >= 0 *before* applying the crossing minus
  sign and then flips the sign of the result. G is therefore <= 0 and mostly
  0 (it is nonzero only through the pseudoscalar and form-factor terms added
  after the clip). The sign-mismatch counts in `wlTestEmAb2DVsGSI` and
  `wlTestEmAb4DVsGSI` show this.

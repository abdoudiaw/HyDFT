# Input file reference (`hydft.in`)

All namelists are optional; absent values keep their defaults. Strings are
case-insensitive.

```
&grid     nx=64 ny=1 nz=1  lx=20 ly=1 lz=1  dealias=.true. /
&species  kind='ion'|'electron'
          gamma=2 kappa=0.1          ! ions: Γ, κ = a/λ_D
          rs=1.86 theta=-1           ! electrons: r_s and Γ (or θ > 0, which overrides Γ)
          mass=1
          terms='ideal,hartree,ry'   ! default per kind: ions as shown, electrons 'tf,gradient,hartree,ry'
          tfk_gamma=0.1111           ! gradient prefactor: 1/9 Kirzhnits, 1 von Weizsäcker
          closure='newtonian'|'maxwell'
          transport='constant'|'yukawa_fit'|'stanton_murillo'|'electron_fit'
          eta=0 xi=0 tau=0           ! η, ξ in m n₀ ω_p a², τ in 1/ω_p (constant); see docs/equations.md
          tau_model='ichimaru'|'constant'
          adiabatic_index=1.6667
          structure='hnc'|'file'|'none'  structure_file='...'  structure_kind='s'|'c' /
&hnc      nr=8192 rmax=80 tol=1e-9 mix=0.2 maxiter=20000 /
&time     dt=0 (auto) cfl=0.4 tend=50 nsteps_max=1e7 iout=0 iprobe=1 idiag=10 /
&initial  type='uniform'|'mode'|'random'|'gaussian' amplitude=1e-4 mode=1,0,0 width=1 seed=1234 /
&external type='none'|'mode'|'driven'|'gaussian' amplitude=0 mode=1,0,0 omega=0 width=1 ramp=0 /
&probe    mode=1,0,0 /
&linear   qmin=0.05 qmax=2.5 nq=50 omega_max=3 nomega=600 q_dsf=0.5,1.0,1.5 /
&output   prefix='out' /
```

Programs:

* `hydft input.in` — time-dependent run. Writes `<prefix>_diag.dat`
  (budgets), `<prefix>_probe.dat` (complex `δn_k/n₀` and `mom_k` of the probe
  mode every `iprobe` steps), and `<prefix>_fld_NNNNNN.bin` fields every `iout` steps.
* `hydft_linear input.in` — `<prefix>_dispersion.dat` (q, Re ω, Im ω, undamped ω, S(q))
  and `<prefix>_dsf.dat` (S(q,ω)/n₀ for each `q_dsf`).
* `hydft_hnc input.in` — standalone structure (`&hnc` namelist with
  `potential`, `gamma`, `kappa`, `rs`, `prefix`): `<prefix>_gr.dat`, `<prefix>_sk.dat`.

Python: `python/hydft_io.py` reads all outputs; `plot_dispersion.py`, `plot_dsf.py`;
`python/stanton_murillo.py` mirrors the Fortran transport fits (reduced and cgs units).
Notebook: `examples/electron_dsf/electron_dsf.ipynb` runs the electron example and
compares with Figs. 1–2 of the Sci. Rep. paper (reconstructs the supplement's η_l in Python).

## Adding a physics option

Free-energy term: extend `fe_term_t` (setup / add_mu / kernel / energy) in a
new file under `src/functional/` and add a `case` in `species_build`.
Likewise for a stress closure (`src/closure/stress.f90`), transport model
(`src/transport/transport.f90`), structure source (`src/structure/`), pair
potential (`pair_potential.f90`), external potential or initial condition.

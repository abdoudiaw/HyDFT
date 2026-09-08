# HyDFT — hydrodynamic density-functional theory for strongly coupled plasmas

`HyDFT` is a modern-Fortran (2018) code that solves the time-dependent,
nonlinear hydrodynamic equations of strongly coupled plasmas closed by dynamic
density functional theory (DDFT). One set of equations, one set of
abstractions, covers both

* the classical **viscoelastic density-functional (VEDF)** model for ions
  [Diaw & Murillo, *Phys. Rev. E* **92**, 013107 (2015)], and
* the **viscous quantum hydrodynamics (DDFT-QHD)** model for electrons
  [Diaw & Murillo, *Sci. Rep.* **7**, 15352 (2017)].

The same free-energy and stress-closure objects that drive the nonlinear solver
also generate the linear response (dispersion relations, susceptibility,
dynamic structure factor), so every simulation can be checked against the
analytic results of the two papers.

The code is organised in the spirit of
[CaNS](https://github.com/CaNS-World/CaNS): small modules, abstract derived
types at every physics extension point, periodic pseudo-spectral discretisation,
a single namelist input file, plain CMake build. MIT licensed.

## Status

Implemented and tested (`ctest`, 12 tests):

* HNC structure (Yukawa, Coulomb, Hansen–McDonald QSP), file input of MD `S(k)`;
* free-energy terms (classical ideal, Thomas–Fermi, Kirzhnits/vW gradient,
  Hartree, Ramakrishnan–Yussouff correlation) with linear kernels verified
  against finite differences of the nonlinear functional derivatives;
* linear response: PRE Eq. 54, the electron plasmon dispersion and DSF;
* the nonlinear solver (Newtonian and Maxwell closures, ions and electrons):
  mass/momentum conservation to round-off, single-mode frequency and damping
  within 1–2 % of the linear theory, driven-mode `χ(k,ω)` to 1e-5, relaxation
  to the DDFT ground state `μ[n] + v_ext = const`.

See `docs/equations.md` for the equations exactly as implemented (including
the choices made for the relaxation time) and `docs/input.md` for the input
reference. Examples for the papers' figures are in `examples/`.

---

## 1. The model

### 1.1 Moment equations from the BBGKY hierarchy

For a species of mass `m` interacting through a pair potential `v(r)` in an
external potential `v_ext(r,t)`, the first BBGKY equation yields the exact
(unclosed) moment equations for the density `n`, velocity `u` and pressure
tensor `P_ij`:

```
D_t n            = − n ∇·u
m n D_t u        = − ∇·P + n F_ext − C
D_t P_ij         = − P_ij ∇_k u_k − P_ik ∇_k u_j − P_kj ∇_k u_i − ∇_k Q_ijk − B_ij
```

with `D_t = ∂_t + u·∇`. The terms `C` and `B_ij` contain the two-body
distribution and therefore all correlation, field and collision effects.

### 1.2 Closure 1 — DDFT for the diagonal part

The DDFT closure replaces the two-body density by its equilibrium value for the
*instantaneous* one-body density, which forces the equilibrium of the fluid
equations to coincide with the density-functional ground state
`δF[n]/δn = μ`. Writing `P_ij = P δ_ij − Π_ij`, the momentum equation becomes

```
m n D_t u = − n ∇ μ[n] + n F_ext + ∇·Π,        μ[n] ≡ δF[n]/δn
```

where the total free-energy functional is a sum of terms:

| term | classical ions (VEDF) | quantum electrons (DDFT-QHD) |
|---|---|---|
| non-interacting | `F_id = β⁻¹ ∫ n (ln Λ³n − 1)` → `μ_id = β⁻¹ ln n` | Thomas–Fermi + Kirzhnits gradient correction `T_TFK[n]` (Eq. 17), `μ = α(n)/β + γ (3√2π²/8) β^{3/2} [ξ'(α)|∇n|² + 2ξ(α)∇²n]` (finite-T Bohm potential) |
| Hartree / mean field | `U = v ⋆ n`, `v(k) = 4πe²/(k²+k_D²)` (Yukawa, Boltzmann electrons) | `U = v ⋆ n`, `v(k) = 4πe²/k²` with neutralising background |
| correlation / exchange–correlation | Ramakrishnan–Yussouff: `F_cor = F_cor[n₀] − (2β)⁻¹ ∬ Δn Δn' c(|r−r'|)` → `μ_cor = − β⁻¹ c_cor ⋆ Δn` | same form with `c_ee`, `Ω_ex = Ω_H + Ω_xc` (Eq. 13) |

Here `Δn = n − n₀`, and `c_cor(k) = c(k) + β v(k)` is the part of the direct
correlation function beyond mean field, so that the Hartree term is not double
counted. `c(k)` comes from an Ornstein–Zernike/HNC calculation (Yukawa, Coulomb,
or a Hansen–McDonald quantum statistical potential), or from a file of `S(k)`
produced by molecular dynamics via `c(k) = [1 − 1/S(k)]/n₀`.

### 1.3 Closure 2 — relaxation of the dissipative stress

The off-diagonal stress `Π_ij` relaxes towards its Newtonian value on the
Maxwell time `τ`:

```
Newtonian (τ = 0):   Π = Π⁰ ≡ η (∇u + ∇uᵀ − ⅔ ∇·u I) + ξ ∇·u I
Maxwell   (τ > 0):   D_t Π = (Π⁰ − Π)/τ
```

Eliminating `Π` from the momentum equation gives the VEDF form of the paper,

```
(1 + τ D_t) [ m n D_t u + ∇P − n F_tot + n ∇ δF_cor/δn ] = ∇·Π⁰,
```

but `HyDFT` integrates `Π` as a state variable instead: it is equivalent, avoids
commuting `∇` with `D_t`, recovers Navier–Stokes exactly as `τ → 0`, and keeps
the high-frequency (elastic) response for `ωτ ≫ 1`. The relaxation time is
`ω_p τ = 3 η̄ / (1 − γμ + (4/15) E_c)` with `η̄ = (4η/3 + ξ)/(m n₀ ω_p a²)`
(Eq. 37–38).

### 1.4 The system that is integrated

Per species, conservative form:

```
∂_t n            + ∇·(n u)        = 0
∂_t (m n u)      + ∇·(m n u u)    = − n ∇ μ[n] + n F_ext + ∇·Π
∂_t Π            + u·∇Π           = (Π⁰ − Π)/τ            (only for the Maxwell closure)
```

### 1.5 Reduced units

Lengths in the Wigner–Seitz radius `a = (3/4πn₀)^{1/3}`, times in `1/ω_p`,
masses in `m`; hence `n₀ = 3/4π`, `k_B T = 1/(3Γ)` with `Γ = e²/(a k_B T)`,
`n₀ v(k) = 1/(k² + κ²)` with `κ = a/λ_D`, wave numbers `q = k a`. For
electrons the density parameter `r_s`, the coupling `Γ`, and the degeneracy
`θ = T/T_F` fix `β` and `α₀`.

### 1.6 Linear response (built from the same objects)

Linearising about `n₀` with each free-energy term contributing its second
functional derivative kernel `K_i(k) = δ²F_i/δn δn |_{n₀}`, and writing
`η_l = 4η/3 + ξ`:

```
χ(k,ω) = n₀ k² / [ − ω² + k² n₀ K(k) − i k² η_l ω / (m n₀ (1 − iωτ)) ]
S(k,ω) = − (1/π) Im χ(k,ω) / (1 − e^{−βħω})          (classical limit: k_B T/(πω) Im χ)
```

Special cases reproduced to round-off:

* classical Yukawa ions, `K = 1/(β n₀ S(q))`:
  `ω²/ω_p² = q²/(3Γ S(q)) − i η̄ q² ω / (1 − iωτ)`   (PRE Eq. 54)
* electrons with TFK + Hartree + `c_ee`:
  `ω_q² ≃ ω_p²[ −(q²/3Γ) n₀ C_ee(q) + λ̄_TF⁻² q² + (ν/4) λ̄_TF⁻⁴ q⁴ − η_l² q⁴/2 ]`  (Sci. Rep. Eq. 24)
  and the full DSF of Sci. Rep. Eq. 22.

Limits: `τ → 0` and `F_cor → 0` give Navier–Stokes; `τ → 0` gives Ying's
hydrodynamic Bloch equations; RPA `c = −βv`, `γ = 0`, `η = 0` gives Bohm–Gross.

---

## 2. Numerics

* Periodic box, Fourier pseudo-spectral in 1-D/2-D/3-D (arrays are always
  3-D, `ny = nz = 1` for 1-D), FFTW3, 2/3 dealiasing. Hartree and correlation
  terms are convolutions and are applied in `k`-space; the Bohm term uses
  spectral gradients.
* Low-storage third-order Runge–Kutta (Wray) for advection and forces;
  the stiff relaxation `(Π⁰−Π)/τ` is applied with an exact integrating factor
  per stage, so `τ → 0` costs nothing in stability.
* Time step from advective, acoustic (`c_s² = n₀ K(k)`), viscous and Bohm
  (`∝ k²`) constraints.
* Spatial operators live behind a single module so that a finite-volume,
  shock-capturing backend can be added later without touching the physics.

## 3. Code layout

```
src/core        kinds, utils, params (namelist), units, grid, fft (FFTW wrapper), spectral ops
src/special     Gauss–Legendre, Fermi–Dirac integrals I_p(α) and inverse of I_{1/2}
src/structure   pair potentials, OZ/HNC solver, structure sources (hnc | file)
src/functional  abstract free-energy term + ideal, TF, gradient, Hartree, RY correlation
src/transport   η, ξ, τ: constant | yukawa_fit (+ Ichimaru τ)
src/closure     stress closures: newtonian | maxwell
src/model       species factory, external potentials, model (RHS + RK3), IC, diagnostics, I/O, analysis
src/linear      χ(k,ω), ω(q), S(k,ω) from the same functional/closure objects
app/            hydft (simulation), hydft_linear (dispersion/DSF), hydft_hnc (structure)
test/           ctest programs
examples/       yukawa_iaw (PRE Figs. 2–4), electron_dsf (Sci. Rep. Figs. 1–2)
python/         readers and plotting
docs/           equations and input reference
```

Every physics choice (species kind, free-energy term, stress closure,
transport model, structure source, external potential) is an abstract derived
type selected by name in the input file; adding one means adding one file and
one `select case` in the factory.

## 4. Building

Requirements: gfortran ≥ 12 (or ifx), CMake ≥ 3.20, FFTW3; optional OpenMP,
HDF5.

```
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
build/hydft_linear examples/yukawa_iaw/g2_k0.1.in     # dispersion + DSF tables
build/hydft examples/yukawa_iaw/g2_k0.1_mode.in       # single-mode simulation
```

Debug builds (`-DCMAKE_BUILD_TYPE=Debug`) enable bounds checking and
floating-point traps.

## 5. References

1. A. Diaw and M. S. Murillo, *Generalized hydrodynamics model for strongly
   coupled plasmas*, Phys. Rev. E **92**, 013107 (2015).
2. A. Diaw and M. S. Murillo, *A viscous quantum hydrodynamics model based on
   dynamic density functional theory*, Sci. Rep. **7**, 15352 (2017).
3. U. M. B. Marconi and P. Tarazona, J. Chem. Phys. **110**, 8032 (1999);
   A. J. Archer and R. Evans, J. Chem. Phys. **121**, 4246 (2004) — DDFT closure.
4. T. V. Ramakrishnan and M. Yussouff, Phys. Rev. B **19**, 2775 (1979) —
   correlation functional.
5. J. Frenkel, *Kinetic Theory of Liquids* (Clarendon, 1946) — Maxwell
   relaxation / generalized hydrodynamics.

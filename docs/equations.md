# HyDFT — equations as implemented

All quantities are in reduced units: lengths in the Wigner–Seitz radius
`a = (3/4πn₀)^{1/3}`, times in `1/ω_p`, masses in `m`. Hence

    n₀ = 3/(4π),   e² = 1/3,   k_B T = 1/(3Γ),   β = 3Γ,   n₀ v(k) = 1/(k² + κ²)

and for electrons `ħ = 1/√(3 r_s)`, `θ = T/T_F = 2 r_s / (Γ (9π/4)^{2/3})`, with
`α₀` from `n₀ = c_n I_{1/2}(α₀)`, `c_n = √2 m^{3/2}/(π² ħ³ β^{3/2})`.

## Evolved system (per species, conservative form)

    ∂_t n        + ∇·(n u)        = 0
    ∂_t (m n u)  + ∇·(m n u u)    = − n ∇μ[n] + n F_ext + ∇·Π
    ∂_t Π        + u·∇Π           = (Π⁰ − Π)/τ            (maxwell closure only)

    Π⁰ = η (∇u + ∇uᵀ − ⅔ ∇·u I) + ξ ∇·u I,   F_ext = −∇v_ext

`newtonian` sets `Π = Π⁰`. The relaxation term is integrated per Runge–Kutta
stage with the exact solution for a source varying linearly across the stage,
so `τ → 0` recovers the Newtonian closure with no stability penalty.

## Free-energy terms and their linear kernels

`μ[n] = Σ_i δF_i/δn`;  `K_i(k) = δ²F_i/δnδn|_{n₀}` is what the linear tool uses.
`test_kernels` checks `K_i` against a finite difference of `μ_i` for every term.

| term | `μ_i[n]` | `K_i(k)` |
|---|---|---|
| `ideal` | `β⁻¹ ln n` | `1/(β n₀)` |
| `tf` | `α(n)/β` | `2/(β c_n I_{−1/2}(α₀))` |
| `gradient` | `−C [ (dξ/dn) |∇n|² + 2 ξ ∇²n ]`, `ξ = I'_{−1/2}/I²_{−1/2}`, `C = γ_g (3√2π²/8) ħ⁵β^{3/2}/m^{5/2}` | `2 C ξ(α₀) k²` |
| `hartree` | `v ⋆ n`, `v(k) = 4πe²/(k²+κ²)` (k = 0 dropped when κ = 0) | `v(k)` |
| `ry` | `−β⁻¹ c_cor ⋆ (n − n₀)` | `−c_cor(k)/β` |

`c_cor = c + βv` when a Hartree term is present (no double counting), else the
full `c(k)`. `c(k)` comes from the HNC solver or from a file of `S(k)` or `c(k)`.

Notes on the gradient term: the sign follows from `T_2 = C ∫ ξ |∇n|²` and gives
exactly the Bohm potential `−(ħ²/2m) ∇²√n/√n` for `γ_g = 1` at `T = 0`
(`γ_g = 1/9` is Kirzhnits). The Thomas–Fermi length is
`λ_TF² = 1/(4π e² ∂n/∂μ)`, which reduces to the Debye length classically; the
constants printed in Sci. Rep. Eq. 20 differ and were re-derived.

## Linear response (same objects)

With perturbations `∝ e^{i(k·r − ωt)}` and `η̄ = (4η/3 + ξ)/(m n₀ ω_p a²)`:

    χ(k,ω) = − k² n₀ / D,    D = −ω² + k² n₀ K(k)/m − i k² η̄ ω /(1 − iωτ)
    ω(k):  D = 0  (cubic in ω for τ > 0; the Re ω > 0 root is reported)
    S(k,ω)/n₀ = −(1/π n₀) Im χ / (1 − e^{−βħω})     (quantum)
              = −(1/π n₀ β ω) Im χ                   (classical)

`χ(k,0) = −1/K(k) < 0`. For classical Yukawa ions `n₀K = 1/(βS(k))` and
`D = 0` is PRE Eq. 54; for electrons it is the full form of Sci. Rep. Eq. 21/24.

## Transport

* `constant`: `η, ξ, τ` from the input.
* `yukawa_fit`: `η* = 0.482/G² + 0.629/G^{0.878} + 1.88×10⁻³ G` (Bastea 2005 OCP
  fit) at `G = Γ(1 + κ + κ²/2)e^{−κ}`; `ξ` from input.
* `tau_model = ichimaru`:
  `ω_p τ = 3Γ η̄ / (1 − γ_ad μ_ex + 4u/15)`, `u = E_c/(NkT)` from HNC,
  `μ_ex = u/3 + (Γ/9) du/dΓ` by finite difference of two extra HNC solves.
  Two choices to confirm against the original scripts: the factor `Γ` (PRE
  Eq. 37 prints `3η̄`; `Γ` arises from `n k_B T` in these units), and the use of
  the excess compressibility (the full one makes the modulus vanish near Γ ≈ 2).

Not built in yet: the electron viscosity interpolation (Conti–Vignale +
Stanton–Murillo) of the Sci. Rep. supplement, and the Stanton–Murillo
`η(Γ,κ), τ(Γ,κ)` collision-integral fits (need the `K11, K22` coefficients).

## Numerics

Periodic pseudo-spectral (FFTW r2c), 2/3 dealiasing of every product and of
the state after each stage; the k = 0 mode of the internal force `−n∇μ` is
removed so momentum is conserved to round-off. Low-storage RK3 (Wray). Time
step: `dt = cfl · min(1.7/ω_max, 2.5/(η̄ k_max²), dx/u_max)` with
`ω_max² = max_k k²[n₀K(k)/m + η̄/(mτ)]`.

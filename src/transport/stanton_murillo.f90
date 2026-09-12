! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Stanton & Murillo transport coefficients for screened one-component plasmas
!> and binary ionic mixtures [L. G. Stanton and M. S. Murillo, Phys. Rev. E 93,
!> 043203 (2016)], from the fits (C22-C24, Table IV) to the effective-Boltzmann
!> collision integrals
!>
!>   Omega_ij^(n,m) = sqrt(2 pi/(mu_ij T^3)) (Z_i Z_j e^2)^2 K_nm(g),   g = Z_i Z_j e^2/(lambda_eff T),
!>
!> with the effective screening length of Eq. (28)/(35), which gives (Eq. 54, 68)
!>
!>   g = Gamma sqrt(kappa^2 + 3 Gamma/(1 + 3 Gamma))                          (one species)
!>   g = Gamma_ij sqrt(kappa^2 + sum_k 3 x_k Gamma_kk/(1 + 3 (x_k/z_k)^{1/3} Gamma_kk))   (mixture)
!>
!> where kappa = a/lambda_e (electron screening), x_k = n_k/n_tot, z_k = Z_k n_k/sum Z_j n_j.
!> (Eq. 68 is printed with x_k^{-1}; expanding 1/(lambda_k^2 (1 + 3 Gamma_k^IS)) with
!> a_tot^2/lambda_k^2 = 3 x_k Gamma_kk gives x_k, which also keeps a trace species
!> from screening infinitely strongly.)
!>
!> Reduced single-species coefficients (Eqs. 56, 75, 82), a = ion-sphere radius:
!>   D*   = D/(wp a^2)      = sqrt(3 pi)/(12 Gamma^{5/2} K_11)
!>   eta* = eta/(m n wp a^2) = 5 sqrt(3 pi)/(36 Gamma^{5/2} K_22)
!>   K*   = K/(n wp a^2)    = 25 sqrt(3 pi)/(48 Gamma^{5/2} K_22)      (k_B = 1)
!>
!> Electron longitudinal viscosity of the Sci. Rep. 7, 15352 (2017) supplement, Eqs. (6-8):
!>   eta_l* = (4/3) (eta_CV* + theta eta_SM*)/(1 + theta),  eta_SM* at g = Gamma a/lambda_TF,
!>   eta_CV* = 1/sqrt(3 r_s) / (60 r_s^{-3/2} + 80 r_s^{-1} - 40 r_s^{-2/3} + 62 r_s^{-1/3}).
!>
!> Validity: the fits reproduce the collision integrals to ~1e-2; compared with MD the
!> effective-Boltzmann D*, eta* agree to ~10-30 % for Gamma <~ 10 and fail at strong coupling
!> (no viscosity minimum). Use the MD fit (transport = 'yukawa_fit') beyond Gamma ~ 10.
module hydft_stanton_murillo
  use hydft_kinds
  implicit none
  private
  public :: sm_knm, sm_k11, sm_k12, sm_k13, sm_k22
  public :: sm_g, sm_g_mixture, sm_diffusion, sm_viscosity, sm_conductivity
  public :: sm_electron_viscosity, sm_mixture
  public :: IDX_K11, IDX_K12, IDX_K13, IDX_K22

  integer, parameter :: IDX_K11 = 1, IDX_K12 = 2, IDX_K13 = 3, IDX_K22 = 4
  integer, parameter :: n_idx(4) = [1, 1, 1, 2], m_idx(4) = [1, 2, 3, 2]
  ! Table IV: weakly coupled (g < 1) a_1..a_5 and strongly coupled (g > 1) b_0..b_4
  real(rp), parameter :: a_tab(5, 4) = reshape([ &
    1.4660_rp, -1.7836_rp, 1.4313_rp, -0.55833_rp, 0.061162_rp, &
    0.52094_rp, 0.25153_rp, -1.1337_rp, 1.2155_rp, -0.43784_rp, &
    0.30346_rp, 0.23739_rp, -0.62167_rp, 0.56110_rp, -0.18046_rp, &
    0.85401_rp, -0.22898_rp, -0.60059_rp, 0.80591_rp, -0.30555_rp], [5, 4])
  real(rp), parameter :: b_tab(0:4, 4) = reshape([ &
    0.081033_rp, -0.091336_rp, 0.051760_rp, -0.50026_rp, 0.17044_rp, &
    0.20572_rp, -0.16536_rp, 0.061572_rp, -0.12770_rp, 0.066993_rp, &
    0.68375_rp, -0.38459_rp, 0.10711_rp, 0.10649_rp, 0.028760_rp, &
    0.43475_rp, -0.21147_rp, 0.11116_rp, 0.19665_rp, 0.15195_rp], [5, 4])
  real(rp), parameter :: sqrt3pi = 3.0699801238394655_rp   ! sqrt(3 pi)

contains

  !> Reduced collision integral K_nm(g), idx = IDX_K11 | IDX_K12 | IDX_K13 | IDX_K22 (Eqs. C22-C24).
  pure real(rp) function sm_knm(idx, g) result(k)
    integer, intent(in) :: idx
    real(rp), intent(in) :: g
    real(rp) :: s, lg, fact
    integer :: i
    fact = 1.0_rp
    do i = 2, m_idx(idx) - 1
      fact = fact*i
    end do
    if (g < 1.0_rp) then
      s = 0.0_rp
      do i = 5, 1, -1
        s = (s + a_tab(i, idx))*max(g, tiny(1.0_rp))
      end do
      k = -0.25_rp*n_idx(idx)*fact*log(s)
    else
      lg = log(g)
      k = (b_tab(0, idx) + b_tab(1, idx)*lg + b_tab(2, idx)*lg*lg)/(1.0_rp + b_tab(3, idx)*g + b_tab(4, idx)*g*g)
    end if
  end function sm_knm

  pure real(rp) function sm_k11(g)
    real(rp), intent(in) :: g
    sm_k11 = sm_knm(IDX_K11, g)
  end function sm_k11

  pure real(rp) function sm_k12(g)
    real(rp), intent(in) :: g
    sm_k12 = sm_knm(IDX_K12, g)
  end function sm_k12

  pure real(rp) function sm_k13(g)
    real(rp), intent(in) :: g
    sm_k13 = sm_knm(IDX_K13, g)
  end function sm_k13

  pure real(rp) function sm_k22(g)
    real(rp), intent(in) :: g
    sm_k22 = sm_knm(IDX_K22, g)
  end function sm_k22

  !> Plasma parameter of one ionic species with electron screening kappa = a/lambda_e (Eq. 54).
  pure real(rp) function sm_g(gamma, kappa) result(g)
    real(rp), intent(in) :: gamma, kappa
    g = gamma*sqrt(kappa*kappa + 3.0_rp*gamma/(1.0_rp + 3.0_rp*gamma))
  end function sm_g

  !> a_tot/lambda_eff for a mixture (the square root of Eq. 68): x = number fractions,
  !> z = charge fractions Z_k n_k/rho_tot, gkk = Z_k^2 e^2/(a_tot T). Multiply by Gamma_ij for g_ij.
  pure real(rp) function sm_g_mixture(kappa, x, z, gkk) result(s)
    real(rp), intent(in) :: kappa, x(:), z(:), gkk(:)
    integer :: k
    s = kappa*kappa
    do k = 1, size(x)
      if (x(k) <= 0.0_rp) cycle
      s = s + 3.0_rp*x(k)*gkk(k)/(1.0_rp + 3.0_rp*(x(k)/z(k))**(1.0_rp/3.0_rp)*gkk(k))
    end do
    s = sqrt(s)
  end function sm_g_mixture

  !> D* = D/(wp a^2), Eq. (56).
  pure real(rp) function sm_diffusion(gamma, kappa) result(d)
    real(rp), intent(in) :: gamma, kappa
    d = sqrt3pi/(12.0_rp*gamma**2.5_rp*sm_k11(sm_g(gamma, kappa)))
  end function sm_diffusion

  !> eta* = eta/(m n wp a^2), Eq. (75).
  pure real(rp) function sm_viscosity(gamma, kappa) result(eta)
    real(rp), intent(in) :: gamma, kappa
    eta = 5.0_rp*sqrt3pi/(36.0_rp*gamma**2.5_rp*sm_k22(sm_g(gamma, kappa)))
  end function sm_viscosity

  !> K* = K/(n k_B wp a^2), Eq. (82).
  pure real(rp) function sm_conductivity(gamma, kappa) result(kth)
    real(rp), intent(in) :: gamma, kappa
    kth = 25.0_rp*sqrt3pi/(48.0_rp*gamma**2.5_rp*sm_k22(sm_g(gamma, kappa)))
  end function sm_conductivity

  !> Electron longitudinal viscosity eta_l* (units n m wp a^2) of the Sci. Rep. supplement:
  !> Conti-Vignale at T = 0 interpolated with Stanton-Murillo at g = Gamma a/lambda_TF.
  pure subroutine sm_electron_viscosity(rs, gamma, theta, lambda_tf, eta_l, eta_cv, eta_sm)
    real(rp), intent(in) :: rs, gamma, theta, lambda_tf   !< lambda_tf in units of a
    real(rp), intent(out) :: eta_l
    real(rp), intent(out), optional :: eta_cv, eta_sm
    real(rp) :: cv, sm
    cv = 1.0_rp/sqrt(3.0_rp*rs)/(60.0_rp*rs**(-1.5_rp) + 80.0_rp/rs - 40.0_rp*rs**(-2.0_rp/3.0_rp) &
         + 62.0_rp*rs**(-1.0_rp/3.0_rp))
    sm = 5.0_rp*sqrt3pi/(36.0_rp*gamma**2.5_rp*sm_k22(gamma/lambda_tf))
    eta_l = 4.0_rp/3.0_rp*(cv + theta*sm)/(1.0_rp + theta)
    if (present(eta_cv)) eta_cv = cv
    if (present(eta_sm)) eta_sm = sm
  end subroutine sm_electron_viscosity

  !> Binary mixture, first-order Chapman-Enskog (Eqs. B8-B16). Inputs: Gamma_12 = Z_1 Z_2 e^2/(a_tot T),
  !> kappa = a_tot/lambda_e, number fraction x1 of species 1, charges z1, z2 and masses m1, m2 (any
  !> common unit). Outputs in reduced mixture units with a_tot and the hydrodynamic plasma frequency
  !> w_hp = sqrt(4 pi <Z>^2 e^2 n_tot/<m>), <Z> = sum x_k Z_k, <m> = sum x_k m_k:
  !> d12* = D_12/(w_hp a_tot^2), eta* = eta_tot/(rho w_hp a_tot^2), k* = K_tot/(n_tot k_B w_hp a_tot^2),
  !> kt = thermal diffusion ratio D_T/D_12 (dimensionless).
  pure subroutine sm_mixture(gamma12, kappa, x1, z1, z2, m1, m2, d12, eta, kth, kt)
    real(rp), intent(in) :: gamma12, kappa, x1, z1, z2, m1, m2
    real(rp), intent(out) :: d12, eta, kth, kt
    real(rp) :: x2, t, n, rhoc, zc(2), s, g11, g22, g12, mu12, m1f, m2f
    real(rp) :: o11_22, o22_22, o12_11, o12_12, o12_13, o12_22
    real(rp) :: eta1, eta2, k1, k2, aa, bb, cc, ee, p1, p2, q1, q2, q12, q12p, r1, r2, r12, r12p, s1, s2
    real(rp) :: zbar, mbar, whp, den
    ! units: a_tot = 1, e^2 = 1
    x2 = 1.0_rp - x1
    t = z1*z2/gamma12
    n = 3.0_rp/(4.0_rp*pi)
    rhoc = x1*z1 + x2*z2
    zc = [x1*z1, x2*z2]/rhoc
    s = sm_g_mixture(kappa, [x1, x2], zc, [z1*z1, z2*z2]/t)
    g11 = z1*z1/t*s; g22 = z2*z2/t*s; g12 = z1*z2/t*s
    mu12 = m1*m2/(m1 + m2)
    m1f = m1/(m1 + m2); m2f = m2/(m1 + m2)
    o11_22 = omega(0.5_rp*m1, z1*z1, sm_k22(g11))
    o22_22 = omega(0.5_rp*m2, z2*z2, sm_k22(g22))
    o12_11 = omega(mu12, z1*z2, sm_k11(g12))
    o12_12 = omega(mu12, z1*z2, sm_k12(g12))
    o12_13 = omega(mu12, z1*z2, sm_k13(g12))
    o12_22 = omega(mu12, z1*z2, sm_k22(g12))
    eta1 = 5.0_rp*t/(8.0_rp*o11_22);  eta2 = 5.0_rp*t/(8.0_rp*o22_22)
    k1 = 75.0_rp*t/(32.0_rp*m1*o11_22);  k2 = 75.0_rp*t/(32.0_rp*m2*o22_22)
    aa = o12_22/(5.0_rp*o12_11)
    bb = (5.0_rp*o12_12 - o12_13)/(5.0_rp*o12_11)
    cc = 2.0_rp*o12_12/(5.0_rp*o12_11) - 1.0_rp
    ee = t/(8.0_rp*m1f*m2f*o12_11)
    p1 = o11_22/(5.0_rp*m2f*o12_11);  p2 = o22_22/(5.0_rp*m1f*o12_11)
    q1 = p1*(6.0_rp*m2f**2 + 5.0_rp*m1f**2 - 4.0_rp*m1f**2*bb + 8.0_rp*m1f*m2f*aa)
    q2 = p2*(6.0_rp*m1f**2 + 5.0_rp*m2f**2 - 4.0_rp*m2f**2*bb + 8.0_rp*m1f*m2f*aa)
    q12p = 15.0_rp*ee*(p1 + p2 + (11.0_rp - 4.0_rp*bb - 8.0_rp*aa)*m1f*m2f)/(2.0_rp*(m1 + m2))
    q12 = 2.0_rp*p1*p2 + 3.0_rp*(m1f - m2f)**2*(5.0_rp - 4.0_rp*bb) + 4.0_rp*m1f*m2f*aa*(11.0_rp - 4.0_rp*bb)
    r1 = 2.0_rp/3.0_rp + m1f/m2f*aa;  r2 = 2.0_rp/3.0_rp + m2f/m1f*aa
    r12 = 4.0_rp*aa/(3.0_rp*m1f*m2f*ee) + ee/(2.0_rp*eta1*eta2)
    r12p = 4.0_rp/3.0_rp + ee/(2.0_rp*eta1) + ee/(2.0_rp*eta2) - 2.0_rp*aa
    s1 = m1f*p1 - m2f*(3.0_rp*(m2f - m1f) + 4.0_rp*m1f*aa)
    s2 = m2f*p2 - m1f*(3.0_rp*(m1f - m2f) + 4.0_rp*m2f*aa)
    den = x1*x1*q1 + x2*x2*q2 + x1*x2*q12
    d12 = 3.0_rp*t/(16.0_rp*n*mu12*o12_11)
    eta = (x1*x1*r1 + x2*x2*r2 + x1*x2*r12p)/(x1*x1*r1/eta1 + x2*x2*r2/eta2 + x1*x2*r12)
    kth = (x1*x1*q1*k1 + x2*x2*q2*k2 + x1*x2*q12p)/den
    kt = 5.0_rp*x1*x2*cc*(x1*s1 - x2*s2)/den
    ! reduced units of the mixture
    zbar = x1*z1 + x2*z2
    mbar = x1*m1 + x2*m2
    whp = sqrt(4.0_rp*pi*zbar*zbar*n/mbar)
    d12 = d12/whp
    eta = eta/(n*mbar*whp)
    kth = kth/(n*whp)
  contains
    pure real(rp) function omega(mu, zz, knm)
      real(rp), intent(in) :: mu, zz, knm
      omega = sqrt(2.0_rp*pi/(mu*t**3))*zz*zz*knm
    end function omega
  end subroutine sm_mixture

end module hydft_stanton_murillo

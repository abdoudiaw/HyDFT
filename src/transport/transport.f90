! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Transport coefficients in reduced units: shear viscosity eta and bulk
!> viscosity xi in units of m n0 wp a^2, Maxwell relaxation time tau in 1/wp.
!>
!> Models (input `transport`)
!>   constant        : eta, xi, tau from the input file.
!>   yukawa_fit      : eta from the OCP MD fit of Bastea [PRE 71, 056405 (2005)],
!>                     eta* = 0.482/G^2 + 0.629/G^0.878 + 1.88e-3 G, at the effective coupling
!>                     G = Gamma (1 + kappa + kappa^2/2) e^{-kappa}; xi from input. Covers the
!>                     viscosity minimum and strong coupling.
!>   stanton_murillo : eta from the effective-Boltzmann fit of Stanton & Murillo [PRE 93, 043203
!>                     (2016)], eta* = 5 sqrt(3 pi)/(36 Gamma^{5/2} K_22(g)), g = Gamma sqrt(kappa^2 +
!>                     3 Gamma/(1 + 3 Gamma)) (module hydft_stanton_murillo); xi from input.
!>                     Accurate for Gamma <~ 10 (weak to moderate coupling); no viscosity minimum.
!>   electron_fit    : electrons only; eta_l of the Sci. Rep. 7, 15352 supplement (Conti-Vignale at
!>                     T = 0 interpolated in theta with Stanton-Murillo at g = Gamma a/lambda_TF),
!>                     xi = 0 and eta = 3 eta_l/4.
!>
!> Relaxation time (input `tau_model`, Maxwell closure only)
!>   constant : tau from input.
!>   ichimaru : the Ichimaru/Kaw-Sen form (PRE Eq. 37)
!>              wp tau = 3 Gamma eta_l/(1 - gamma_ad mu_ex + 4 u/15), with u = E_c/(N kT)
!>              and mu_ex = u/3 + (Gamma/9) du/dGamma the excess compressibility, both
!>              from the HNC structure calculation. Notes: (i) Eq. 37 prints 3 eta_bar;
!>              the factor Gamma comes from n kT = n0 m wp^2 a^2/(3 Gamma) in these
!>              units. (ii) With the full compressibility 1 + mu_ex (Kaw & Sen 1998) the
!>              modulus vanishes near Gamma ~ 2; the excess form is positive for all
!>              Gamma and identical at strong coupling. Confirm against the original
!>              scripts before quantitative comparisons with the figures.
module hydft_transport
  use hydft_kinds
  use hydft_utils, only: lower
  use hydft_structure_source
  use hydft_units
  use hydft_stanton_murillo
  implicit none
  private
  public :: transport_t

  type :: transport_t
    character(len=16) :: model = 'constant'
    character(len=16) :: tau_model = 'constant'
    real(rp) :: eta = 0.0_rp
    real(rp) :: xi = 0.0_rp
    real(rp) :: tau = 0.0_rp
    real(rp) :: eta_l = 0.0_rp        !< 4 eta/3 + xi
    real(rp) :: eta_cv = 0.0_rp, eta_sm = 0.0_rp   !< electron_fit parts
    real(rp) :: elastic_modulus = 0.0_rp   !< 1 - gamma_ad mu_c + 4u/15 (dimensionless), when computed
  contains
    procedure :: setup => transport_setup
    procedure :: print => transport_print
  end type transport_t

contains

  subroutine transport_setup(tr, model, tau_model, eta_in, xi_in, tau_in, un, gamma_ad, st, maxwell)
    class(transport_t), intent(inout) :: tr
    character(len=*), intent(in) :: model, tau_model
    real(rp), intent(in) :: eta_in, xi_in, tau_in, gamma_ad
    type(units_t), intent(in) :: un
    type(structure_t), intent(in) :: st
    logical, intent(in) :: maxwell
    real(rp) :: geff, u, du, muc, gamma, kappa
    gamma = un%gamma
    kappa = un%kappa
    tr%model = lower(model)
    tr%tau_model = lower(tau_model)
    select case (trim(tr%model))
    case ('constant')
      tr%eta = eta_in
      tr%xi = xi_in
    case ('yukawa_fit')
      geff = gamma*(1.0_rp + kappa + 0.5_rp*kappa*kappa)*exp(-kappa)
      tr%eta = 0.482_rp/geff**2 + 0.629_rp/geff**0.878_rp + 1.88e-3_rp*geff
      tr%xi = xi_in
    case ('stanton_murillo')
      if (un%quantum) call fatal('transport = stanton_murillo is the ion (Yukawa) fit; use electron_fit for electrons')
      tr%eta = sm_viscosity(gamma, kappa)
      tr%xi = xi_in
    case ('electron_fit')
      if (.not. un%quantum) call fatal('transport = electron_fit needs kind = electron')
      call sm_electron_viscosity(un%rs, gamma, un%theta, un%lambda_tf, tr%eta_l, tr%eta_cv, tr%eta_sm)
      tr%eta = 0.75_rp*tr%eta_l
      tr%xi = 0.0_rp
    case default
      call fatal('unknown transport model: '//trim(model))
    end select
    tr%eta_l = 4.0_rp*tr%eta/3.0_rp + tr%xi
    tr%tau = 0.0_rp
    if (maxwell) then
      select case (trim(tr%tau_model))
      case ('constant')
        tr%tau = tau_in
      case ('ichimaru')
        if (.not. st%has_energy) call fatal('tau_model = ichimaru needs the HNC structure source (for E_c)')
        u = st%u_ex
        du = st%du_dgamma
        muc = u/3.0_rp + gamma*du/9.0_rp
        tr%elastic_modulus = 1.0_rp - gamma_ad*muc + 4.0_rp*u/15.0_rp
        if (tr%elastic_modulus <= 0.0_rp) call fatal('ichimaru tau: 1 - gamma mu + 4u/15 <= 0; use tau_model = constant')
        tr%tau = 3.0_rp*gamma*tr%eta_l/tr%elastic_modulus
      case default
        call fatal('unknown tau model: '//trim(tau_model))
      end select
    end if
  end subroutine transport_setup

  subroutine transport_print(tr)
    class(transport_t), intent(in) :: tr
    write(*,'(a,a,a,f9.5,a,f9.5,a,f9.5,a,f10.4)') '  transport ', trim(tr%model), '  eta =', tr%eta, '  xi =', tr%xi, &
      '  eta_l =', tr%eta_l, '  tau =', tr%tau
    if (tr%elastic_modulus /= 0.0_rp) write(*,'(a,f10.4)') '            1 - gamma mu + 4u/15 =', tr%elastic_modulus
    if (tr%model == 'electron_fit') write(*,'(a,f9.5,a,f9.5)') '            eta_CV =', tr%eta_cv, '  eta_SM =', tr%eta_sm
  end subroutine transport_print

end module hydft_transport

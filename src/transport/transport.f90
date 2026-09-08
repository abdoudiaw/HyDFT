! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Transport coefficients in reduced units: shear viscosity eta and bulk
!> viscosity xi in units of m n0 wp a^2, Maxwell relaxation time tau in 1/wp.
!>
!> Models
!>   constant   : eta, xi, tau from the input file.
!>   yukawa_fit : eta from the OCP fit of Bastea [PRE 71, 056405 (2005)],
!>                eta* = 0.482/G^2 + 0.629/G^0.878 + 1.88e-3 G, evaluated at the
!>                effective coupling G = Gamma (1 + kappa + kappa^2/2) e^{-kappa};
!>                xi from input; tau from the Ichimaru/Kaw-Sen form (PRE Eq. 37)
!>                wp tau = 3 Gamma eta_l/(1 - gamma_ad mu_ex + 4 u/15), with u = E_c/(N kT)
!>                and mu_ex = u/3 + (Gamma/9) du/dGamma the excess compressibility, both
!>                from the HNC structure calculation. Notes: (i) Eq. 37 prints 3 eta_bar;
!>                the factor Gamma comes from n kT = n0 m wp^2 a^2/(3 Gamma) in these
!>                units. (ii) With the full compressibility 1 + mu_ex (Kaw & Sen 1998) the
!>                modulus vanishes near Gamma ~ 2; the excess form is positive for all
!>                Gamma and identical at strong coupling. Confirm against the original
!>                scripts before quantitative comparisons with the figures.
module hydft_transport
  use hydft_kinds
  use hydft_utils, only: lower
  use hydft_structure_source
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
    real(rp) :: elastic_modulus = 0.0_rp   !< 1 - gamma_ad mu_c + 4u/15 (dimensionless), when computed
  contains
    procedure :: setup => transport_setup
    procedure :: print => transport_print
  end type transport_t

contains

  subroutine transport_setup(tr, model, tau_model, eta_in, xi_in, tau_in, gamma, kappa, gamma_ad, st, maxwell)
    class(transport_t), intent(inout) :: tr
    character(len=*), intent(in) :: model, tau_model
    real(rp), intent(in) :: eta_in, xi_in, tau_in, gamma, kappa, gamma_ad
    type(structure_t), intent(in) :: st
    logical, intent(in) :: maxwell
    real(rp) :: geff, u, du, muc
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
  end subroutine transport_print

end module hydft_transport

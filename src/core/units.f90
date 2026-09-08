! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Reduced units. Lengths in a = (3/4 pi n0)^{1/3}, time in 1/wp, mass m = 1,
!> so n0 = 3/(4 pi), e^2 = 1/3, k_B T = 1/(3 Gamma). For electrons
!> hbar = 1/sqrt(3 rs) and the reduced chemical potential alpha0 follows from
!> n0 = cn I_{1/2}(alpha0), cn = sqrt(2)/(pi^2 hbar^3 beta^{3/2}).
module hydft_units
  use hydft_kinds
  use hydft_fermi_dirac
  implicit none
  private
  public :: units_t

  type :: units_t
    logical  :: quantum = .false.
    real(rp) :: n0 = 3.0_rp/(4.0_rp*pi)
    real(rp) :: e2 = 1.0_rp/3.0_rp
    real(rp) :: gamma = 1.0_rp
    real(rp) :: kappa = 0.0_rp
    real(rp) :: beta = 3.0_rp
    real(rp) :: kT = 1.0_rp/3.0_rp
    real(rp) :: mass = 1.0_rp
    real(rp) :: rs = 1.0_rp
    real(rp) :: hbar = 0.0_rp
    real(rp) :: theta = huge(1.0_rp)
    real(rp) :: alpha0 = 0.0_rp
    real(rp) :: cn = 0.0_rp          !< n = cn I_{1/2}(alpha)
    real(rp) :: lambda_qsp = 0.0_rp  !< electron-pair thermal de Broglie length / a
  contains
    procedure :: init => units_init
    procedure :: print => units_print
  end type units_t

contains

  subroutine units_init(u, kind, gamma, kappa, rs, mass)
    class(units_t), intent(inout) :: u
    character(len=*), intent(in) :: kind
    real(rp), intent(in) :: gamma, kappa, rs, mass
    real(rp) :: tf
    u%gamma = gamma
    u%kappa = kappa
    u%mass = mass
    u%beta = 3.0_rp*gamma
    u%kT = 1.0_rp/u%beta
    u%rs = rs
    if (kind == 'electron') then
      u%quantum = .true.
      u%kappa = 0.0_rp
      u%hbar = 1.0_rp/sqrt(3.0_rp*rs)
      u%cn = sqrt2/(pi*pi*u%hbar**3*u%beta**1.5_rp)*u%mass**1.5_rp
      u%alpha0 = fd_alpha_from_ip12(u%n0/u%cn)
      ! T_F = hbar^2 (3 pi^2 n0)^{2/3}/(2 m)
      tf = u%hbar**2*(3.0_rp*pi*pi*u%n0)**(2.0_rp/3.0_rp)/(2.0_rp*u%mass)
      u%theta = u%kT/tf
      u%lambda_qsp = 2.0_rp*sqrt(pi*gamma/rs)
    else
      u%quantum = .false.
      u%hbar = 0.0_rp
      u%theta = huge(1.0_rp)
    end if
  end subroutine units_init

  subroutine units_print(u)
    class(units_t), intent(in) :: u
    write(*,'(a,f10.4,a,f10.4,a,f8.4)') '  units     Gamma =', u%gamma, '  beta =', u%beta, '  kappa =', u%kappa
    if (u%quantum) write(*,'(a,f8.4,a,es10.3,a,f10.4,a,f10.4)') '            rs =', u%rs, '  hbar =', u%hbar, '  theta =', u%theta, '  alpha0 =', u%alpha0
  end subroutine units_print

end module hydft_units

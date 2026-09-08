! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Pair potentials in reduced units: beta*v(r) with r in Wigner-Seitz radii.
!> Each potential provides the full beta*v(r) and a long-range/short-range
!> split, beta*v = beta*v_S + beta*v_L, with beta*v_L(r) = G [e^{-kappa r} - e^{-kappa' r}]/r
!> (kappa' = kappa + ksplit) whose Fourier transform is analytic. This makes
!> the Ornstein-Zernike solver work for Coulomb and Yukawa systems alike.
module hydft_pair_potential
  use hydft_kinds
  use hydft_utils, only: lower
  implicit none
  private
  public :: pair_potential_t, yukawa_t, coulomb_t, hansen_mcdonald_t, make_pair_potential

  type, abstract :: pair_potential_t
    real(rp) :: gamma = 1.0_rp     !< coupling prefactor of the 1/r tail (Gamma)
    real(rp) :: kappa = 0.0_rp     !< inverse screening length (0 for Coulomb)
    real(rp) :: ksplit = 1.0_rp    !< kappa' - kappa for the long/short split
    character(len=32) :: name = ''
  contains
    procedure(bv_iface), deferred :: beta_v       !< full beta v(r)
    procedure :: beta_v_long_r
    procedure :: beta_v_long_k
    procedure :: beta_v_short
    procedure :: beta_v_k                          !< beta v(k) of the 1/r tail (mean field)
  end type pair_potential_t

  abstract interface
    pure real(rp) function bv_iface(p, r)
      import :: pair_potential_t, rp
      class(pair_potential_t), intent(in) :: p
      real(rp), intent(in) :: r
    end function bv_iface
  end interface

  type, extends(pair_potential_t) :: yukawa_t
  contains
    procedure :: beta_v => yukawa_beta_v
  end type yukawa_t

  type, extends(pair_potential_t) :: coulomb_t
  contains
    procedure :: beta_v => coulomb_beta_v
  end type coulomb_t

  !> Hansen-McDonald quantum statistical potential for electrons:
  !> beta v = (G/r)(1 - e^{-2 pi r/L}) + ln2 exp(-4 pi r^2/(L^2 ln2)), L = thermal
  !> de Broglie wavelength of the electron pair in units of a.
  type, extends(pair_potential_t) :: hansen_mcdonald_t
    real(rp) :: lambda = 1.0_rp
    logical :: pauli = .true.
  contains
    procedure :: beta_v => hm_beta_v
  end type hansen_mcdonald_t

contains

  pure real(rp) function yukawa_beta_v(p, r) result(v)
    class(yukawa_t), intent(in) :: p
    real(rp), intent(in) :: r
    v = p%gamma*exp(-p%kappa*r)/r
  end function yukawa_beta_v

  pure real(rp) function coulomb_beta_v(p, r) result(v)
    class(coulomb_t), intent(in) :: p
    real(rp), intent(in) :: r
    v = p%gamma/r
  end function coulomb_beta_v

  pure real(rp) function hm_beta_v(p, r) result(v)
    class(hansen_mcdonald_t), intent(in) :: p
    real(rp), intent(in) :: r
    real(rp), parameter :: ln2 = 0.693147180559945309417232121458_rp
    if (r > 0.0_rp) then
      v = p%gamma*(1.0_rp - exp(-twopi*r/p%lambda))/r
    else
      v = p%gamma*twopi/p%lambda
    end if
    if (p%pauli) v = v + ln2*exp(-4.0_rp*pi*r*r/(p%lambda**2*ln2))
  end function hm_beta_v

  !> Long-range part in r space (finite at r = 0).
  pure real(rp) function beta_v_long_r(p, r) result(v)
    class(pair_potential_t), intent(in) :: p
    real(rp), intent(in) :: r
    real(rp) :: kp
    kp = p%kappa + p%ksplit
    if (r > 1.0e-12_rp) then
      v = p%gamma*(exp(-p%kappa*r) - exp(-kp*r))/r
    else
      v = p%gamma*(kp - p%kappa)
    end if
  end function beta_v_long_r

  !> Long-range part in k space: 4 pi G [1/(k^2+kappa^2) - 1/(k^2+kappa'^2)].
  pure real(rp) function beta_v_long_k(p, k) result(v)
    class(pair_potential_t), intent(in) :: p
    real(rp), intent(in) :: k
    real(rp) :: kp
    kp = p%kappa + p%ksplit
    v = 4.0_rp*pi*p%gamma*(1.0_rp/(k*k + p%kappa**2) - 1.0_rp/(k*k + kp*kp))
  end function beta_v_long_k

  pure real(rp) function beta_v_short(p, r) result(v)
    class(pair_potential_t), intent(in) :: p
    real(rp), intent(in) :: r
    v = p%beta_v(r) - p%beta_v_long_r(r)
  end function beta_v_short

  !> Mean-field kernel beta v(k) = 4 pi G/(k^2 + kappa^2) of the screened Coulomb tail.
  pure real(rp) function beta_v_k(p, k) result(v)
    class(pair_potential_t), intent(in) :: p
    real(rp), intent(in) :: k
    v = 4.0_rp*pi*p%gamma/(k*k + p%kappa**2)
  end function beta_v_k

  !> Factory.
  subroutine make_pair_potential(name, gamma, kappa, lambda, pot)
    character(len=*), intent(in) :: name
    real(rp), intent(in) :: gamma, kappa, lambda
    class(pair_potential_t), allocatable, intent(out) :: pot
    select case (trim(lower(name)))
    case ('yukawa')
      allocate(yukawa_t :: pot)
    case ('coulomb')
      allocate(coulomb_t :: pot)
    case ('hansen_mcdonald', 'hm', 'qsp')
      allocate(hansen_mcdonald_t :: pot)
      select type (pot)
      type is (hansen_mcdonald_t)
        pot%lambda = lambda
      end select
    case default
      call fatal('unknown pair potential: '//trim(name))
    end select
    pot%gamma = gamma
    pot%kappa = kappa
    pot%name = name
  end subroutine make_pair_potential

end module hydft_pair_potential

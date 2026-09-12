! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Linear response built from the same functional / transport / closure
!> objects as the nonlinear solver. With perturbations ~ exp(i(k.r - w t)):
!>
!>   chi(k,w) = -k^2 n0 / D(k,w),   D = -w^2 + k^2 n0 K(k) - i k^2 eta_bar w/(1 - i w tau)
!>   dispersion: D = 0
!>   S(k,w)/n0 = -(1/(pi n0)) Im chi/(1 - exp(-beta hbar w))    (quantum)
!>            = -(1/(pi n0 beta w)) Im chi                       (classical)
!>
!> chi is the density response to a potential energy delta v_ext; chi(k,0) = -1/K < 0.
!> For classical Yukawa ions K = 1/(beta n0 S(k)) and D = 0 is PRE Eq. 54;
!> for electrons D = 0 is the full form of Sci. Rep. Eq. 21/24.
module hydft_linear_response
  use hydft_kinds
  use hydft_species
  implicit none
  private
  public :: linear_t

  type :: linear_t
    real(rp) :: n0 = 3.0_rp/(4.0_rp*pi)
    real(rp) :: beta = 3.0_rp
    real(rp) :: hbar = 0.0_rp
    logical  :: quantum = .false.
    real(rp) :: eta_bar = 0.0_rp     !< eta_l/(m n0 wp a^2)
    real(rp) :: tau = 0.0_rp
    real(rp) :: mass = 1.0_rp
    type(species_t), pointer :: sp => null()
  contains
    procedure :: init => linear_init
    procedure :: c2 => linear_c2
    procedure :: denominator => linear_denominator
    procedure :: chi => linear_chi
    procedure :: dispersion => linear_dispersion
    procedure :: dsf => linear_dsf
  end type linear_t

contains

  subroutine linear_init(l, sp)
    class(linear_t), intent(inout) :: l
    type(species_t), intent(in), target :: sp
    l%sp => sp
    l%n0 = sp%un%n0
    l%beta = sp%un%beta
    l%hbar = sp%un%hbar
    l%quantum = sp%un%quantum
    l%eta_bar = sp%tr%eta_l
    l%tau = sp%stress%tau
    l%mass = sp%mass
  end subroutine linear_init

  !> Squared sound speed-like quantity n0 K(k)/m, so that w0^2 = k^2 c2.
  real(rp) function linear_c2(l, k) result(c2)
    class(linear_t), intent(in) :: l
    real(rp), intent(in) :: k
    c2 = l%n0*l%sp%fun%kernel(k)/l%mass
  end function linear_c2

  complex(rp) function linear_denominator(l, k, w) result(d)
    class(linear_t), intent(in) :: l
    real(rp), intent(in) :: k
    complex(rp), intent(in) :: w
    d = -w*w + k*k*l%c2(k) - i_unit*k*k*(l%eta_bar/l%mass)*w/(1.0_rp - i_unit*w*l%tau)
  end function linear_denominator

  complex(rp) function linear_chi(l, k, w) result(chi)
    class(linear_t), intent(in) :: l
    real(rp), intent(in) :: k, w
    chi = -k*k*l%n0/l%denominator(k, cmplx(w, 0.0_rp, rp))
  end function linear_chi

  !> Complex mode frequency: the root of D(k,w) = 0 with the largest real part
  !> (Re w > 0 branch of the propagating pair; Im w < 0 is damping).
  !> All roots of the cubic  -i tau w^3 + w^2 + i (k^2 eta_bar + c^2 k^2 tau) w - c^2 k^2 = 0
  !> (or the quadratic for tau = 0) are found by Durand-Kerner iteration.
  complex(rp) function linear_dispersion(l, k, all_roots) result(w)
    class(linear_t), intent(in) :: l
    real(rp), intent(in) :: k
    complex(rp), intent(out), optional :: all_roots(3)
    complex(rp) :: a(0:3), z(3), zn(3), num, den
    integer :: deg, i, j, it
    real(rp) :: c2k2, nu
    c2k2 = k*k*l%c2(k)
    nu = k*k*l%eta_bar/l%mass
    if (l%tau > 0.0_rp) then
      deg = 3
      a(3) = -i_unit*l%tau; a(2) = (1.0_rp, 0.0_rp); a(1) = i_unit*(nu + c2k2*l%tau); a(0) = -c2k2
    else
      deg = 2
      a(2) = (1.0_rp, 0.0_rp); a(1) = i_unit*nu; a(0) = -c2k2
    end if
    ! Durand-Kerner
    z(1) = cmplx(0.4_rp, 0.9_rp, rp); z(2) = z(1)**2; z(3) = z(1)**3
    do it = 1, 500
      do i = 1, deg
        num = poly(a, deg, z(i))
        den = a(deg)
        do j = 1, deg
          if (j /= i) den = den*(z(i) - z(j))
        end do
        zn(i) = z(i) - num/den
      end do
      if (maxval(abs(zn(1:deg) - z(1:deg))) < 1.0e-15_rp*max(1.0_rp, maxval(abs(zn(1:deg))))) then
        z = zn; exit
      end if
      z = zn
    end do
    ! polish with Newton
    do i = 1, deg
      do it = 1, 5
        z(i) = z(i) - poly(a, deg, z(i))/dpoly(a, deg, z(i))
      end do
    end do
    w = z(1)
    do i = 2, deg
      if (real(z(i)) > real(w)) w = z(i)
    end do
    if (present(all_roots)) then
      all_roots = (0.0_rp, 0.0_rp)
      all_roots(1:deg) = z(1:deg)
    end if
  contains
    pure complex(rp) function poly(a, deg, x)
      complex(rp), intent(in) :: a(0:3), x
      integer, intent(in) :: deg
      integer :: m
      poly = a(deg)
      do m = deg - 1, 0, -1
        poly = poly*x + a(m)
      end do
    end function poly
    pure complex(rp) function dpoly(a, deg, x)
      complex(rp), intent(in) :: a(0:3), x
      integer, intent(in) :: deg
      integer :: m
      dpoly = deg*a(deg)
      do m = deg - 1, 1, -1
        dpoly = dpoly*x + m*a(m)
      end do
    end function dpoly
  end function linear_dispersion

  !> Dynamic structure factor per particle, S(k,w)/n0 (units 1/wp).
  !> Im chi is odd in w, so -Im chi/w is used: it stays finite at w = 0, where
  !> Im chi -> -n0 eta_bar w/(m c2^2) and the Bose factor x/(1-exp(-x)) -> 1 + x/2.
  real(rp) function linear_dsf(l, k, w) result(s)
    class(linear_t), intent(in) :: l
    real(rp), intent(in) :: k, w
    real(rp) :: ratio, x
    if (abs(w) < 1.0e-6_rp) then
      ratio = l%n0*l%eta_bar/(l%mass*l%c2(k)**2)
    else
      ratio = -aimag(l%chi(k, w))/w
    end if
    if (l%quantum) then
      x = l%beta*l%hbar*w
      if (abs(x) < 1.0e-6_rp) then
        s = ratio/(pi*l%n0*l%beta*l%hbar)*(1.0_rp + 0.5_rp*x)
      else
        s = ratio/(pi*l%n0*l%beta*l%hbar)*x/(1.0_rp - exp(-x))
      end if
    else
      s = ratio/(pi*l%n0*l%beta)
    end if
  end function linear_dsf

end module hydft_linear_response

! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> The total free-energy functional: a list of terms. Provides the total
!> mu[n] on the grid, the total kernel K(k), and the total energy.
module hydft_functional
  use hydft_kinds
  use hydft_spectral
  use hydft_free_energy
  implicit none
  private
  public :: functional_t

  type :: functional_t
    integer :: nterms = 0
    type(fe_box_t), allocatable :: terms(:)
    real(rp) :: n0 = 3.0_rp/(4.0_rp*pi)
  contains
    procedure :: add => functional_add
    procedure :: setup => functional_setup
    procedure :: mu => functional_mu
    procedure :: kernel => functional_kernel
    procedure :: energy => functional_energy
    procedure :: has => functional_has
  end type functional_t

contains

  subroutine functional_add(f, term)
    class(functional_t), intent(inout) :: f
    class(fe_term_t), intent(in) :: term
    type(fe_box_t), allocatable :: tmp(:)
    integer :: i
    allocate(tmp(f%nterms + 1))
    do i = 1, f%nterms
      call move_alloc(f%terms(i)%t, tmp(i)%t)
    end do
    allocate(tmp(f%nterms + 1)%t, source=term)
    call move_alloc(tmp, f%terms)
    f%nterms = f%nterms + 1
  end subroutine functional_add

  subroutine functional_setup(f, spc)
    class(functional_t), intent(inout) :: f
    type(spec_t), intent(inout) :: spc
    integer :: i
    do i = 1, f%nterms
      call f%terms(i)%t%setup(spc)
    end do
  end subroutine functional_setup

  !> mu[n] = sum_i dF_i/dn on the grid.
  subroutine functional_mu(f, spc, n, nhat, mu)
    class(functional_t), intent(inout) :: f
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    real(rp), intent(out) :: mu(:,:,:)
    integer :: i
    mu = 0.0_rp
    do i = 1, f%nterms
      call f%terms(i)%t%add_mu(spc, n, nhat, mu)
    end do
  end subroutine functional_mu

  !> K(k) = sum_i K_i(k); n0 K(k) k^2 is the squared frequency of the undamped mode.
  real(rp) function functional_kernel(f, k) result(kk)
    class(functional_t), intent(in) :: f
    real(rp), intent(in) :: k
    integer :: i
    kk = 0.0_rp
    do i = 1, f%nterms
      kk = kk + f%terms(i)%t%kernel(k)
    end do
  end function functional_kernel

  real(rp) function functional_energy(f, spc, n, nhat) result(e)
    class(functional_t), intent(inout) :: f
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    integer :: i
    e = 0.0_rp
    do i = 1, f%nterms
      e = e + f%terms(i)%t%energy(spc, n, nhat)
    end do
  end function functional_energy

  logical function functional_has(f, name)
    class(functional_t), intent(in) :: f
    character(len=*), intent(in) :: name
    integer :: i
    functional_has = .false.
    do i = 1, f%nterms
      if (trim(f%terms(i)%t%name) == trim(name)) functional_has = .true.
    end do
  end function functional_has

end module hydft_functional

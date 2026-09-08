! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> External potential energy v_ext(r,t) per particle and its force -grad v_ext.
!>   none     : nothing
!>   mode     : A cos(k.r)                       (static single mode)
!>   driven   : A(t) cos(k.r - w t), A(t) = A min(1, t/ramp)   (for chi(k,w) checks)
!>   gaussian : -A exp(-|r - L/2|^2/(2 width^2))                (attractive well)
module hydft_external_potential
  use hydft_kinds
  use hydft_utils, only: lower
  use hydft_spectral
  implicit none
  private
  public :: vext_t

  type :: vext_t
    character(len=16) :: kind = 'none'
    real(rp) :: amp = 0.0_rp
    integer  :: mode(3) = [1, 0, 0]
    real(rp) :: omega = 0.0_rp
    real(rp) :: width = 1.0_rp
    real(rp) :: ramp = 0.0_rp
    real(rp) :: kvec(3) = 0.0_rp
    logical  :: active = .false.
    logical  :: time_dependent = .false.
    real(rp), allocatable :: v(:,:,:)
    real(rp), allocatable :: f(:,:,:,:)
  contains
    procedure :: setup => vext_setup
    procedure :: eval => vext_eval
  end type vext_t

contains

  subroutine vext_setup(e, spc, kind, amp, mode, omega, width, ramp)
    class(vext_t), intent(inout) :: e
    type(spec_t), intent(inout) :: spc
    character(len=*), intent(in) :: kind
    real(rp), intent(in) :: amp, omega, width, ramp
    integer, intent(in) :: mode(3)
    e%kind = lower(kind); e%amp = amp; e%mode = mode; e%omega = omega; e%width = width; e%ramp = ramp
    e%kvec = [twopi*mode(1)/spc%g%lx, twopi*mode(2)/spc%g%ly, twopi*mode(3)/spc%g%lz]
    if (spc%g%ny == 1) e%kvec(2) = 0.0_rp
    if (spc%g%nz == 1) e%kvec(3) = 0.0_rp
    select case (trim(e%kind))
    case ('none')
      e%active = .false.
    case ('mode', 'gaussian')
      e%active = .true.; e%time_dependent = .false.
    case ('driven')
      e%active = .true.; e%time_dependent = .true.
    case default
      call fatal('unknown external potential: '//trim(kind))
    end select
    if (allocated(e%v)) deallocate(e%v, e%f)
    allocate(e%v(spc%g%nx, spc%g%ny, spc%g%nz), e%f(spc%g%nx, spc%g%ny, spc%g%nz, 3))
    e%v = 0.0_rp; e%f = 0.0_rp
    if (e%active) call e%eval(spc, 0.0_rp)
  end subroutine vext_setup

  subroutine vext_eval(e, spc, t)
    class(vext_t), intent(inout) :: e
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: t
    integer :: i, j, k
    real(rp) :: a, ph, r2, x, y, z
    if (.not. e%active) return
    a = e%amp
    if (e%time_dependent .and. e%ramp > 0.0_rp) a = a*min(1.0_rp, t/e%ramp)
    select case (trim(e%kind))
    case ('mode', 'driven')
      do k = 1, spc%g%nz; do j = 1, spc%g%ny; do i = 1, spc%g%nx
        ph = e%kvec(1)*spc%g%x(i) + e%kvec(2)*spc%g%y(j) + e%kvec(3)*spc%g%z(k) - e%omega*t
        e%v(i,j,k) = a*cos(ph)
        e%f(i,j,k,:) = a*sin(ph)*e%kvec
      end do; end do; end do
    case ('gaussian')
      do k = 1, spc%g%nz; do j = 1, spc%g%ny; do i = 1, spc%g%nx
        x = spc%g%x(i) - 0.5_rp*spc%g%lx
        y = 0.0_rp; z = 0.0_rp
        if (spc%g%ny > 1) y = spc%g%y(j) - 0.5_rp*spc%g%ly
        if (spc%g%nz > 1) z = spc%g%z(k) - 0.5_rp*spc%g%lz
        r2 = x*x + y*y + z*z
        e%v(i,j,k) = -a*exp(-0.5_rp*r2/e%width**2)
      end do; end do; end do
      call spc%grad(e%v, e%f)
      e%f = -e%f
    end select
  end subroutine vext_eval

end module hydft_external_potential

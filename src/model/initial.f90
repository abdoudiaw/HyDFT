! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Initial conditions.
!>   uniform  : n = n0, u = 0
!>   mode     : n = n0 (1 + A cos(k.r)), u = 0
!>   random   : n = n0 (1 + A d), d a zero-mean, unit-rms dealiased random field
!>   gaussian : n = n0 (1 + A (g - <g>)), g = exp(-|r - L/2|^2/(2 width^2))
module hydft_initial
  use hydft_kinds
  use hydft_utils, only: lower
  use hydft_spectral
  use hydft_species
  implicit none
  private
  public :: set_initial_condition

contains

  subroutine set_initial_condition(spc, sp, kind, amp, mode, width, seed)
    type(spec_t), intent(inout) :: spc
    type(species_t), intent(inout) :: sp
    character(len=*), intent(in) :: kind
    real(rp), intent(in) :: amp, width
    integer, intent(in) :: mode(3), seed
    integer :: i, j, k, nseed
    integer, allocatable :: sd(:)
    real(rp) :: kv(3), ph, x, y, z, r2, mean, rms
    real(rp), allocatable :: d(:,:,:)
    allocate(d(spc%g%nx, spc%g%ny, spc%g%nz))
    sp%mom = 0.0_rp
    if (allocated(sp%pi)) sp%pi = 0.0_rp
    select case (trim(lower(kind)))
    case ('uniform')
      d = 0.0_rp
    case ('mode')
      kv = [twopi*mode(1)/spc%g%lx, twopi*mode(2)/spc%g%ly, twopi*mode(3)/spc%g%lz]
      if (spc%g%ny == 1) kv(2) = 0.0_rp
      if (spc%g%nz == 1) kv(3) = 0.0_rp
      do k = 1, spc%g%nz; do j = 1, spc%g%ny; do i = 1, spc%g%nx
        ph = kv(1)*spc%g%x(i) + kv(2)*spc%g%y(j) + kv(3)*spc%g%z(k)
        d(i,j,k) = cos(ph)
      end do; end do; end do
    case ('random')
      call random_seed(size=nseed)
      allocate(sd(nseed))
      sd = seed + 37*[(i, i = 1, nseed)]
      call random_seed(put=sd)
      call random_number(d)
      d = d - 0.5_rp
      call spc%dealias(d)
      mean = spc%mean(d)
      d = d - mean
      rms = sqrt(spc%mean(d*d))
      if (rms > 0.0_rp) d = d/rms
    case ('gaussian')
      do k = 1, spc%g%nz; do j = 1, spc%g%ny; do i = 1, spc%g%nx
        x = spc%g%x(i) - 0.5_rp*spc%g%lx
        y = 0.0_rp; z = 0.0_rp
        if (spc%g%ny > 1) y = spc%g%y(j) - 0.5_rp*spc%g%ly
        if (spc%g%nz > 1) z = spc%g%z(k) - 0.5_rp*spc%g%lz
        r2 = x*x + y*y + z*z
        d(i,j,k) = exp(-0.5_rp*r2/width**2)
      end do; end do; end do
      call spc%dealias(d)
      d = d - spc%mean(d)
    case default
      call fatal('unknown initial condition: '//trim(kind))
    end select
    sp%n = sp%un%n0*(1.0_rp + amp*d)
    if (minval(sp%n) <= 0.0_rp) call fatal('initial density is not positive; reduce the amplitude')
  end subroutine set_initial_condition

end module hydft_initial

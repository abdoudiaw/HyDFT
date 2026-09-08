! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Precision kinds and mathematical constants used throughout HyDFT.
module hydft_kinds
  use, intrinsic :: iso_fortran_env, only: real64, int64
  implicit none
  private
  public :: rp, ip, pi, twopi, sqrt2, sqrtpi, i_unit, stderr, stdout, fatal, warn

  integer, parameter :: rp = real64
  integer, parameter :: ip = int64
  real(rp), parameter :: pi     = 3.14159265358979323846264338327950288_rp
  real(rp), parameter :: twopi  = 2.0_rp*pi
  real(rp), parameter :: sqrt2  = 1.41421356237309504880168872420969808_rp
  real(rp), parameter :: sqrtpi = 1.77245385090551602729816748334114518_rp
  complex(rp), parameter :: i_unit = (0.0_rp, 1.0_rp)
  integer, parameter :: stderr = 0, stdout = 6

contains

  !> Print a message and stop with a non-zero code.
  subroutine fatal(msg)
    character(len=*), intent(in) :: msg
    write(stderr,'(a)') 'hydft: fatal: '//trim(msg)
    error stop 1
  end subroutine fatal

  subroutine warn(msg)
    character(len=*), intent(in) :: msg
    write(stderr,'(a)') 'hydft: warning: '//trim(msg)
  end subroutine warn

end module hydft_kinds

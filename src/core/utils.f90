! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Small helpers: string handling, interpolation, file units.
module hydft_utils
  use hydft_kinds
  implicit none
  private
  public :: lower, interp_linear, interp_linear_arr, open_file, trapz

contains

  !> Lower-case copy of a string.
  pure function lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, c
    t = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) t(i:i) = achar(c + 32)
    end do
  end function lower

  !> Piecewise-linear interpolation of y(x) at xq; x must be increasing.
  !> Beyond the table the end values are extrapolated linearly from the last two points.
  pure real(rp) function interp_linear(x, y, xq) result(yq)
    real(rp), intent(in) :: x(:), y(:), xq
    integer :: n, lo, hi, mid
    real(rp) :: w
    n = size(x)
    if (n == 1) then
      yq = y(1); return
    end if
    if (xq <= x(1)) then
      lo = 1
    else if (xq >= x(n)) then
      lo = n - 1
    else
      lo = 1; hi = n
      do while (hi - lo > 1)
        mid = (lo + hi)/2
        if (x(mid) <= xq) then
          lo = mid
        else
          hi = mid
        end if
      end do
    end if
    w = (xq - x(lo))/(x(lo+1) - x(lo))
    yq = (1.0_rp - w)*y(lo) + w*y(lo+1)
  end function interp_linear

  pure subroutine interp_linear_arr(x, y, xq, yq)
    real(rp), intent(in) :: x(:), y(:), xq(:)
    real(rp), intent(out) :: yq(:)
    integer :: i
    do i = 1, size(xq)
      yq(i) = interp_linear(x, y, xq(i))
    end do
  end subroutine interp_linear_arr

  !> Trapezoidal rule on a uniform grid with spacing h.
  pure real(rp) function trapz(f, h)
    real(rp), intent(in) :: f(:), h
    integer :: n
    n = size(f)
    trapz = h*(sum(f) - 0.5_rp*(f(1) + f(n)))
  end function trapz

  !> Open a file for reading or writing and return the unit; stops on error.
  integer function open_file(path, action) result(u)
    character(len=*), intent(in) :: path, action
    integer :: ios
    if (action == 'read') then
      open(newunit=u, file=path, status='old', action='read', iostat=ios)
    else
      open(newunit=u, file=path, status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) call fatal('cannot open file '//trim(path))
  end function open_file

end module hydft_utils

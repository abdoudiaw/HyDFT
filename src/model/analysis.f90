! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Post-processing helpers: fit frequency and damping of a damped
!> oscillation from a time series (zero crossings and peak envelope), and
!> the complex response amplitude of a driven mode.
module hydft_analysis
  use hydft_kinds
  implicit none
  private
  public :: fit_damped_mode, fit_driven_response

contains

  !> y(t) ~ A exp(gamma t) cos(omega t + phi) for t >= t_start.
  !> omega from the mean spacing of zero crossings (linearly interpolated),
  !> gamma from a least-squares line through ln|peak| vs t_peak.
  subroutine fit_damped_mode(t, y, t_start, omega, gamma, nzeros)
    real(rp), intent(in) :: t(:), y(:), t_start
    real(rp), intent(out) :: omega, gamma
    integer, intent(out) :: nzeros
    real(rp), allocatable :: tz(:), tp(:), lp(:)
    integer :: i, n, nz, np, i0
    real(rp) :: w, sx, sy, sxx, sxy, ymax, tmax
    n = size(t)
    allocate(tz(n), tp(n), lp(n))
    i0 = 1
    do while (i0 < n .and. t(i0) < t_start)
      i0 = i0 + 1
    end do
    nz = 0
    do i = i0, n - 1
      if (y(i) == 0.0_rp) cycle
      if (y(i)*y(i+1) < 0.0_rp) then
        nz = nz + 1
        w = y(i)/(y(i) - y(i+1))
        tz(nz) = t(i) + w*(t(i+1) - t(i))
      end if
    end do
    nzeros = nz
    omega = 0.0_rp; gamma = 0.0_rp
    if (nz < 3) return
    omega = pi*(nz - 1)/(tz(nz) - tz(1))
    ! peaks between consecutive zero crossings
    np = 0
    do i = 1, nz - 1
      ymax = 0.0_rp; tmax = 0.0_rp
      call peak_between(tz(i), tz(i+1), ymax, tmax)
      if (ymax > 0.0_rp) then
        np = np + 1
        tp(np) = tmax; lp(np) = log(ymax)
      end if
    end do
    if (np < 2) return
    sx = sum(tp(1:np)); sy = sum(lp(1:np)); sxx = sum(tp(1:np)**2); sxy = sum(tp(1:np)*lp(1:np))
    gamma = (np*sxy - sx*sy)/(np*sxx - sx*sx)
  contains
    subroutine peak_between(ta, tb, ymax, tmax)
      real(rp), intent(in) :: ta, tb
      real(rp), intent(out) :: ymax, tmax
      integer :: j
      ymax = 0.0_rp; tmax = 0.5_rp*(ta + tb)
      do j = 1, n
        if (t(j) > ta .and. t(j) < tb) then
          if (abs(y(j)) > ymax) then
            ymax = abs(y(j)); tmax = t(j)
          end if
        end if
      end do
    end subroutine peak_between
  end subroutine fit_damped_mode

  !> For a drive ~ exp(-i w t): the response coefficient a(t) = z(t) exp(i w t)
  !> averaged over t >= t_start, and its relative scatter (transient check).
  subroutine fit_driven_response(t, z, w, t_start, a, scatter)
    real(rp), intent(in) :: t(:), w, t_start
    complex(rp), intent(in) :: z(:)
    complex(rp), intent(out) :: a
    real(rp), intent(out) :: scatter
    integer :: i, cnt
    complex(rp) :: s
    real(rp) :: dev
    s = (0.0_rp, 0.0_rp); cnt = 0
    do i = 1, size(t)
      if (t(i) >= t_start) then
        s = s + z(i)*exp(i_unit*w*t(i)); cnt = cnt + 1
      end if
    end do
    a = s/max(cnt, 1)
    dev = 0.0_rp
    do i = 1, size(t)
      if (t(i) >= t_start) dev = max(dev, abs(z(i)*exp(i_unit*w*t(i)) - a))
    end do
    scatter = dev/max(abs(a), tiny(1.0_rp))
  end subroutine fit_driven_response

end module hydft_analysis

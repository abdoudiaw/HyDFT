! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Fermi-Dirac integrals: limits, derivative relations, inverse.
program test_fermi_dirac
  use hydft_kinds
  use hydft_fermi_dirac
  implicit none
  integer :: nfail, ia
  real(rp) :: a, v, ex, err, y, d, hh
  nfail = 0
  ! classical limit alpha = -20: I_p ~ Gamma(p+1) e^alpha
  a = -20.0_rp
  err = abs(fd_ip12(a)/(0.5_rp*sqrtpi*exp(a)) - 1.0_rp); call check(err < 1e-8_rp, 'I_1/2 classical', err)
  err = abs(fd_im12(a)/(sqrtpi*exp(a)) - 1.0_rp);        call check(err < 1e-8_rp, 'I_-1/2 classical', err)
  err = abs(fd_ip32(a)/(0.75_rp*sqrtpi*exp(a)) - 1.0_rp); call check(err < 1e-8_rp, 'I_3/2 classical', err)
  ! degenerate limit alpha = 250 (inside the table): Sommerfeld
  a = 250.0_rp
  ex = (2.0_rp/3.0_rp)*a**1.5_rp*(1.0_rp + pi**2/(8.0_rp*a*a))
  err = abs(fd_ip12(a)/ex - 1.0_rp); call check(err < 1e-8_rp, 'I_1/2 degenerate', err)
  ex = 2.0_rp*sqrt(a)*(1.0_rp - pi**2/(24.0_rp*a*a))
  err = abs(fd_im12(a)/ex - 1.0_rp); call check(err < 1e-8_rp, 'I_-1/2 degenerate', err)
  ! continuity across the table edge
  err = abs(fd_ip12(300.0_rp - 1e-9_rp)/fd_ip12(300.0_rp + 1e-9_rp) - 1.0_rp); call check(err < 1e-8_rp, 'continuity at amax', err)
  err = abs(fd_ip12(-30.0_rp - 1e-9_rp)/fd_ip12(-30.0_rp + 1e-9_rp) - 1.0_rp); call check(err < 1e-8_rp, 'continuity at amin', err)
  ! known value: I_{1/2}(0) = 0.6780938...
  err = abs(fd_ip12(0.0_rp) - 0.678093895_rp); call check(err < 1e-8_rp, 'I_1/2(0)', err)
  ! derivative relations by central differences
  hh = 1e-4_rp
  do ia = 0, 9
    a = -5.0_rp + 4.7_rp*ia
    d = (fd_ip12(a + hh) - fd_ip12(a - hh))/(2*hh)
    err = abs(d/(0.5_rp*fd_im12(a)) - 1.0_rp); call check(err < 1e-7_rp, 'dI_1/2 = I_-1/2 / 2', err)
    d = (fd_ip32(a + hh) - fd_ip32(a - hh))/(2*hh)
    err = abs(d/(1.5_rp*fd_ip12(a)) - 1.0_rp); call check(err < 1e-7_rp, 'dI_3/2 = 3 I_1/2 / 2', err)
    d = (fd_im12(a + hh) - fd_im12(a - hh))/(2*hh)
    err = abs(d/fd_dim12(a) - 1.0_rp); call check(err < 1e-6_rp, 'dI_-1/2 tabulated', err)
    d = (fd_dim12(a + hh) - fd_dim12(a - hh))/(2*hh)
    err = abs(d/fd_d2im12(a) - 1.0_rp); call check(err < 1e-5_rp, 'd2I_-1/2 tabulated', err)
    d = (fd_xi(a + hh) - fd_xi(a - hh))/(2*hh)
    err = abs(d/fd_dxi(a) - 1.0_rp); call check(err < 1e-5_rp, 'dxi/dalpha', err)
  end do
  ! inverse
  do ia = 0, 22
    a = -25.0_rp + 13.3_rp*ia
    y = fd_ip12(a)
    v = fd_alpha_from_ip12(y)
    err = abs(v - a); call check(err < 1e-9_rp*max(1.0_rp, abs(a)), 'inverse I_1/2', err)
  end do
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_fermi_dirac: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_fermi_dirac: all passed'
contains
  subroutine check(ok, name, val)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: name
    real(rp), intent(in) :: val
    if (ok) then
      write(*,'(a,a,es12.3)') '  pass  ', name, val
    else
      write(*,'(a,a,es12.3)') '  FAIL  ', name, val
      nfail = nfail + 1
    end if
  end subroutine check
end program test_fermi_dirac

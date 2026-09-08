! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> FFT round trip and mode amplitudes in 1-D and 3-D.
program test_fft
  use hydft_kinds
  use hydft_spectral
  implicit none
  type(spec_t) :: s
  real(rp), allocatable :: u(:,:,:), v(:,:,:)
  complex(rp), allocatable :: uhat(:,:,:)
  integer :: i, j, k, nfail
  real(rp) :: err
  complex(rp) :: a

  nfail = 0
  ! ---- 1-D
  call s%init(64, 1, 1, twopi, twopi, twopi)
  allocate(u(64,1,1), v(64,1,1), uhat(s%g%nkx,1,1))
  do i = 1, 64
    u(i,1,1) = 1.5_rp + 0.3_rp*cos(3.0_rp*s%g%x(i)) + 0.2_rp*sin(5.0_rp*s%g%x(i))
  end do
  call s%fwd(u, uhat)
  call s%bwd(uhat, v)
  err = maxval(abs(u - v))
  call check(err < 1.0e-13_rp, '1-D round trip', err)
  call check(abs(uhat(1,1,1) - 1.5_rp) < 1.0e-13_rp, '1-D mean', abs(uhat(1,1,1) - 1.5_rp))
  ! cos(3x) -> coefficient 0.15 at index 4; sin(5x) -> -0.1 i at index 6
  call check(abs(uhat(4,1,1) - 0.15_rp) < 1.0e-13_rp, '1-D cos mode', abs(uhat(4,1,1) - 0.15_rp))
  call check(abs(uhat(6,1,1) + 0.1_rp*i_unit) < 1.0e-13_rp, '1-D sin mode', abs(uhat(6,1,1) + 0.1_rp*i_unit))
  a = s%mode(u, 3, 0, 0)
  call check(abs(a - 0.15_rp) < 1.0e-13_rp, '1-D mode()', abs(a - 0.15_rp))
  deallocate(u, v, uhat)

  ! ---- 3-D
  call s%init(16, 12, 8, twopi, 4.0_rp, 1.0_rp)
  allocate(u(16,12,8), v(16,12,8), uhat(s%g%nkx,12,8))
  do k = 1, 8; do j = 1, 12; do i = 1, 16
    u(i,j,k) = 0.7_rp + cos(2.0_rp*s%g%x(i))*sin(twopi*s%g%y(j)/4.0_rp) + 0.4_rp*cos(twopi*2.0_rp*s%g%z(k))
  end do; end do; end do
  call s%fwd(u, uhat)
  call s%bwd(uhat, v)
  err = maxval(abs(u - v))
  call check(err < 1.0e-13_rp, '3-D round trip', err)
  a = s%mode(u, 0, 0, 2)
  call check(abs(a - 0.2_rp) < 1.0e-13_rp, '3-D z mode', abs(a - 0.2_rp))
  ! cos(2x) sin(ky y): amplitude 1 -> coefficient at (kx=2, ky=+1) is (1/2)(-i/2) = -i/4
  a = s%mode(u, 2, 1, 0)
  call check(abs(a + 0.25_rp*i_unit) < 1.0e-13_rp, '3-D xy mode', abs(a + 0.25_rp*i_unit))

  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_fft: ', nfail, ' failure(s)'
    error stop 1
  end if
  write(*,'(a)') 'test_fft: all passed'

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
end program test_fft

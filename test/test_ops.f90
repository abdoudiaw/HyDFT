! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Spectral derivative, Laplacian, divergence, convolution and dealiasing.
program test_ops
  use hydft_kinds
  use hydft_spectral
  implicit none
  type(spec_t) :: s
  integer, parameter :: nx = 32, ny = 16, nz = 1
  real(rp), allocatable :: u(:,:,:), du(:,:,:,:), lu(:,:,:), ex(:,:,:), v(:,:,:,:), dv(:,:,:), cu(:,:,:), kern(:,:,:)
  complex(rp), allocatable :: uhat(:,:,:)
  integer :: i, j, nfail
  real(rp) :: err, x, y, lx, ly

  nfail = 0
  lx = twopi; ly = 4.0_rp
  call s%init(nx, ny, nz, lx, ly, 1.0_rp)
  allocate(u(nx,ny,nz), du(nx,ny,nz,3), lu(nx,ny,nz), ex(nx,ny,nz), v(nx,ny,nz,3), dv(nx,ny,nz), cu(nx,ny,nz))
  allocate(kern(s%g%nkx,ny,nz), uhat(s%g%nkx,ny,nz))

  do j = 1, ny; do i = 1, nx
    x = s%g%x(i); y = s%g%y(j)
    u(i,j,1) = sin(2.0_rp*x)*cos(twopi*y/ly)
  end do; end do

  call s%grad(u, du)
  do j = 1, ny; do i = 1, nx
    x = s%g%x(i); y = s%g%y(j)
    ex(i,j,1) = 2.0_rp*cos(2.0_rp*x)*cos(twopi*y/ly)
  end do; end do
  err = maxval(abs(du(:,:,:,1) - ex)); call check(err < 1e-12_rp, 'd/dx', err)
  do j = 1, ny; do i = 1, nx
    x = s%g%x(i); y = s%g%y(j)
    ex(i,j,1) = -(twopi/ly)*sin(2.0_rp*x)*sin(twopi*y/ly)
  end do; end do
  err = maxval(abs(du(:,:,:,2) - ex)); call check(err < 1e-12_rp, 'd/dy', err)
  err = maxval(abs(du(:,:,:,3))); call check(err < 1e-14_rp, 'd/dz = 0 in 2-D', err)

  call s%lap(u, lu)
  ex = -(4.0_rp + (twopi/ly)**2)*u
  err = maxval(abs(lu - ex)); call check(err < 1e-11_rp, 'laplacian', err)

  v(:,:,:,1) = du(:,:,:,1); v(:,:,:,2) = du(:,:,:,2); v(:,:,:,3) = 0.0_rp
  call s%div(v, dv)
  err = maxval(abs(dv - ex)); call check(err < 1e-11_rp, 'div grad = lap', err)

  ! convolution with kernel 1/(1+k^2): eigenfunction test
  kern = 1.0_rp/(1.0_rp + s%g%k2)
  call s%fwd(u, uhat)
  call s%convolve(kern, uhat, cu)
  ex = u/(1.0_rp + 4.0_rp + (twopi/ly)**2)
  err = maxval(abs(cu - ex)); call check(err < 1e-13_rp, 'convolution', err)

  ! dealiasing: cos(12x) on nx=32 lies beyond 2/3*16 = 10.67 and must vanish
  do j = 1, ny; do i = 1, nx
    u(i,j,1) = cos(12.0_rp*s%g%x(i)) + cos(10.0_rp*s%g%x(i))
  end do; end do
  call s%dealias(u)
  do j = 1, ny; do i = 1, nx
    ex(i,j,1) = cos(10.0_rp*s%g%x(i))
  end do; end do
  err = maxval(abs(u - ex)); call check(err < 1e-13_rp, 'dealias removes k > 2/3 kmax', err)
  err = abs(s%g%kmax_dealiased - sqrt(100.0_rp + (twopi*5.0_rp/ly)**2)); call check(err < 1e-12_rp, 'kmax_dealiased', err)

  ! mean
  u = 3.0_rp + du(:,:,:,1)
  err = abs(s%mean(u) - 3.0_rp); call check(err < 1e-13_rp, 'mean', err)

  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_ops: ', nfail, ' failure(s)'
    error stop 1
  end if
  write(*,'(a)') 'test_ops: all passed'
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
end program test_ops

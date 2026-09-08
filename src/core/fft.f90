! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Thin FFTW3 wrapper: real-to-complex 3-D transforms on the grid.
!> Forward transforms are normalised so that spectral coefficients are
!> amplitudes (a plane wave of unit amplitude has coefficient 1).
!> This is the only file that needs to change to move to a distributed FFT.
module hydft_fft
  use, intrinsic :: iso_c_binding
  use hydft_kinds
  implicit none
  private
  public :: fft_t

  include 'fftw3.f03'

  type :: fft_t
    integer :: nx = 0, ny = 0, nz = 0, nkx = 0
    type(c_ptr) :: plan_fwd = c_null_ptr, plan_bwd = c_null_ptr
    real(c_double), pointer :: rwork(:,:,:) => null()
    complex(c_double_complex), pointer :: cwork(:,:,:) => null()
    type(c_ptr) :: rbuf = c_null_ptr, cbuf = c_null_ptr
    real(rp) :: norm = 1.0_rp
  contains
    procedure :: init => fft_init
    procedure :: fwd => fft_fwd
    procedure :: bwd => fft_bwd
    procedure :: destroy => fft_destroy
  end type fft_t

contains

  subroutine fft_init(f, nx, ny, nz)
    class(fft_t), intent(inout) :: f
    integer, intent(in) :: nx, ny, nz
    f%nx = nx; f%ny = ny; f%nz = nz; f%nkx = nx/2 + 1
    f%norm = 1.0_rp/real(nx*ny*nz, rp)
    f%rbuf = fftw_alloc_real(int(nx*ny*nz, c_size_t))
    f%cbuf = fftw_alloc_complex(int(f%nkx*ny*nz, c_size_t))
    call c_f_pointer(f%rbuf, f%rwork, [nx, ny, nz])
    call c_f_pointer(f%cbuf, f%cwork, [f%nkx, ny, nz])
    ! FFTW uses row-major ordering: the Fortran (nx,ny,nz) array is a C (nz,ny,nx) array.
    f%plan_fwd = fftw_plan_dft_r2c_3d(nz, ny, nx, f%rwork, f%cwork, FFTW_MEASURE)
    f%plan_bwd = fftw_plan_dft_c2r_3d(nz, ny, nx, f%cwork, f%rwork, FFTW_MEASURE)
    if (.not. c_associated(f%plan_fwd) .or. .not. c_associated(f%plan_bwd)) call fatal('fft: plan creation failed')
  end subroutine fft_init

  !> Forward transform: real u(nx,ny,nz) -> complex uhat(nkx,ny,nz), normalised.
  subroutine fft_fwd(f, u, uhat)
    class(fft_t), intent(inout) :: f
    real(rp), intent(in) :: u(:,:,:)
    complex(rp), intent(out) :: uhat(:,:,:)
    f%rwork = u
    call fftw_execute_dft_r2c(f%plan_fwd, f%rwork, f%cwork)
    uhat = f%cwork*f%norm
  end subroutine fft_fwd

  !> Backward transform: complex uhat -> real u. The input is preserved.
  subroutine fft_bwd(f, uhat, u)
    class(fft_t), intent(inout) :: f
    complex(rp), intent(in) :: uhat(:,:,:)
    real(rp), intent(out) :: u(:,:,:)
    f%cwork = uhat
    call fftw_execute_dft_c2r(f%plan_bwd, f%cwork, f%rwork)
    u = f%rwork
  end subroutine fft_bwd

  subroutine fft_destroy(f)
    class(fft_t), intent(inout) :: f
    if (c_associated(f%plan_fwd)) call fftw_destroy_plan(f%plan_fwd)
    if (c_associated(f%plan_bwd)) call fftw_destroy_plan(f%plan_bwd)
    if (c_associated(f%rbuf)) call fftw_free(f%rbuf)
    if (c_associated(f%cbuf)) call fftw_free(f%cbuf)
    f%plan_fwd = c_null_ptr; f%plan_bwd = c_null_ptr
    f%rbuf = c_null_ptr; f%cbuf = c_null_ptr
    nullify(f%rwork, f%cwork)
  end subroutine fft_destroy

end module hydft_fft

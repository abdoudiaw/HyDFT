! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Spectral operators on the periodic grid: gradient, divergence, Laplacian,
!> convolution with an isotropic kernel, and dealiasing. Everything the
!> physics modules need from the discretisation goes through this type, so a
!> different spatial backend can be swapped in behind the same interface.
module hydft_spectral
  use hydft_kinds
  use hydft_grid
  use hydft_fft
  implicit none
  private
  public :: spec_t

  type :: spec_t
    type(grid_t) :: g
    type(fft_t)  :: fft
    complex(rp), allocatable, private :: w1(:,:,:), w2(:,:,:)
  contains
    procedure :: init => spec_init
    procedure :: fwd => spec_fwd
    procedure :: bwd => spec_bwd
    procedure :: grad_hat => spec_grad_hat
    procedure :: grad => spec_grad
    procedure :: lap => spec_lap
    procedure :: div => spec_div
    procedure :: convolve => spec_convolve
    procedure :: dealias => spec_dealias
    procedure :: dealias_hat => spec_dealias_hat
    procedure :: mean => spec_mean
    procedure :: mode => spec_mode
  end type spec_t

contains

  subroutine spec_init(s, nx, ny, nz, lx, ly, lz, dealias)
    class(spec_t), intent(inout) :: s
    integer,  intent(in) :: nx, ny, nz
    real(rp), intent(in) :: lx, ly, lz
    logical,  intent(in), optional :: dealias
    call s%g%init(nx, ny, nz, lx, ly, lz, dealias)
    call s%fft%destroy()
    call s%fft%init(nx, ny, nz)
    if (allocated(s%w1)) deallocate(s%w1, s%w2)
    allocate(s%w1(s%g%nkx, ny, nz), s%w2(s%g%nkx, ny, nz))
  end subroutine spec_init

  subroutine spec_fwd(s, u, uhat)
    class(spec_t), intent(inout) :: s
    real(rp), intent(in) :: u(:,:,:)
    complex(rp), intent(out) :: uhat(:,:,:)
    call s%fft%fwd(u, uhat)
  end subroutine spec_fwd

  subroutine spec_bwd(s, uhat, u)
    class(spec_t), intent(inout) :: s
    complex(rp), intent(in) :: uhat(:,:,:)
    real(rp), intent(out) :: u(:,:,:)
    call s%fft%bwd(uhat, u)
  end subroutine spec_bwd

  !> Spectral gradient component idir (1,2,3) of uhat, returned in real space.
  subroutine spec_grad_hat(s, uhat, idir, du)
    class(spec_t), intent(inout) :: s
    complex(rp), intent(in) :: uhat(:,:,:)
    integer, intent(in) :: idir
    real(rp), intent(out) :: du(:,:,:)
    integer :: i, j, k
    select case (idir)
    case (1)
      do k = 1, s%g%nz; do j = 1, s%g%ny; do i = 1, s%g%nkx
        s%w1(i,j,k) = i_unit*s%g%kx(i)*uhat(i,j,k)*s%g%mask(i,j,k)
      end do; end do; end do
    case (2)
      do k = 1, s%g%nz; do j = 1, s%g%ny; do i = 1, s%g%nkx
        s%w1(i,j,k) = i_unit*s%g%ky(j)*uhat(i,j,k)*s%g%mask(i,j,k)
      end do; end do; end do
    case (3)
      do k = 1, s%g%nz; do j = 1, s%g%ny; do i = 1, s%g%nkx
        s%w1(i,j,k) = i_unit*s%g%kz(k)*uhat(i,j,k)*s%g%mask(i,j,k)
      end do; end do; end do
    case default
      call fatal('spectral: bad direction')
    end select
    call s%fft%bwd(s%w1, du)
  end subroutine spec_grad_hat

  !> Gradient of a real field: du(:,:,:,1:3).
  subroutine spec_grad(s, u, du)
    class(spec_t), intent(inout) :: s
    real(rp), intent(in) :: u(:,:,:)
    real(rp), intent(out) :: du(:,:,:,:)
    integer :: d
    call s%fft%fwd(u, s%w2)
    do d = 1, 3
      if (d > s%g%ndim) then
        du(:,:,:,d) = 0.0_rp
      else
        call s%grad_hat(s%w2, d, du(:,:,:,d))
      end if
    end do
  end subroutine spec_grad

  !> Laplacian of a real field.
  subroutine spec_lap(s, u, lu)
    class(spec_t), intent(inout) :: s
    real(rp), intent(in) :: u(:,:,:)
    real(rp), intent(out) :: lu(:,:,:)
    call s%fft%fwd(u, s%w2)
    s%w2 = -s%g%k2*s%w2*s%g%mask
    call s%fft%bwd(s%w2, lu)
  end subroutine spec_lap

  !> Divergence of a real vector field v(:,:,:,1:3).
  subroutine spec_div(s, v, dv)
    class(spec_t), intent(inout) :: s
    real(rp), intent(in) :: v(:,:,:,:)
    real(rp), intent(out) :: dv(:,:,:)
    integer :: d, i, j, k
    s%w1 = (0.0_rp, 0.0_rp)
    do d = 1, s%g%ndim
      call s%fft%fwd(v(:,:,:,d), s%w2)
      select case (d)
      case (1)
        do k = 1, s%g%nz; do j = 1, s%g%ny; do i = 1, s%g%nkx
          s%w1(i,j,k) = s%w1(i,j,k) + i_unit*s%g%kx(i)*s%w2(i,j,k)
        end do; end do; end do
      case (2)
        do k = 1, s%g%nz; do j = 1, s%g%ny; do i = 1, s%g%nkx
          s%w1(i,j,k) = s%w1(i,j,k) + i_unit*s%g%ky(j)*s%w2(i,j,k)
        end do; end do; end do
      case (3)
        do k = 1, s%g%nz; do j = 1, s%g%ny; do i = 1, s%g%nkx
          s%w1(i,j,k) = s%w1(i,j,k) + i_unit*s%g%kz(k)*s%w2(i,j,k)
        end do; end do; end do
      end select
    end do
    s%w1 = s%w1*s%g%mask
    call s%fft%bwd(s%w1, dv)
  end subroutine spec_div

  !> Convolution (kernel * u)(r) with kernel given on the half spectrum.
  subroutine spec_convolve(s, kernel, uhat, cu)
    class(spec_t), intent(inout) :: s
    real(rp), intent(in) :: kernel(:,:,:)
    complex(rp), intent(in) :: uhat(:,:,:)
    real(rp), intent(out) :: cu(:,:,:)
    s%w1 = kernel*uhat*s%g%mask
    call s%fft%bwd(s%w1, cu)
  end subroutine spec_convolve

  !> Apply the dealiasing mask to a real field in place.
  subroutine spec_dealias(s, u)
    class(spec_t), intent(inout) :: s
    real(rp), intent(inout) :: u(:,:,:)
    call s%fft%fwd(u, s%w2)
    s%w2 = s%w2*s%g%mask
    call s%fft%bwd(s%w2, u)
  end subroutine spec_dealias

  subroutine spec_dealias_hat(s, uhat)
    class(spec_t), intent(inout) :: s
    complex(rp), intent(inout) :: uhat(:,:,:)
    uhat = uhat*s%g%mask
  end subroutine spec_dealias_hat

  !> Volume mean of a real field.
  real(rp) function spec_mean(s, u)
    class(spec_t), intent(in) :: s
    real(rp), intent(in) :: u(:,:,:)
    spec_mean = sum(u)/real(s%g%nx*s%g%ny*s%g%nz, rp)
  end function spec_mean

  !> Complex amplitude of mode (mx,my,mz) (integer indices, mx >= 0) of a real field.
  complex(rp) function spec_mode(s, u, mx, my, mz)
    class(spec_t), intent(inout) :: s
    real(rp), intent(in) :: u(:,:,:)
    integer, intent(in) :: mx, my, mz
    integer :: jy, kz
    call s%fft%fwd(u, s%w2)
    jy = modulo(my, s%g%ny) + 1
    kz = modulo(mz, s%g%nz) + 1
    spec_mode = s%w2(mx+1, jy, kz)
  end function spec_mode

end module hydft_spectral

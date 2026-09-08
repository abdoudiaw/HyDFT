! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Periodic Cartesian grid and its wavenumbers. Arrays are always 3-D;
!> 1-D and 2-D runs set ny = nz = 1 (or nz = 1).
module hydft_grid
  use hydft_kinds
  implicit none
  private
  public :: grid_t

  type :: grid_t
    integer  :: nx = 1, ny = 1, nz = 1      !< real-space points
    integer  :: nkx = 1                      !< nx/2 + 1 (r2c half spectrum)
    integer  :: ndim = 1
    real(rp) :: lx = twopi, ly = twopi, lz = twopi
    real(rp) :: dx = 1.0_rp, dy = 1.0_rp, dz = 1.0_rp
    real(rp) :: volume = 1.0_rp
    real(rp), allocatable :: x(:), y(:), z(:)       !< node coordinates (cell centres)
    real(rp), allocatable :: kx(:), ky(:), kz(:)    !< wavenumbers, kx(1:nkx) >= 0
    real(rp), allocatable :: k2(:,:,:)              !< |k|^2 on the half spectrum
    real(rp), allocatable :: kabs(:,:,:)            !< |k|
    real(rp), allocatable :: mask(:,:,:)            !< 2/3-rule dealiasing mask (1 or 0)
    real(rp) :: kmax_dealiased = 0.0_rp             !< largest |k| kept after dealiasing
  contains
    procedure :: init => grid_init
  end type grid_t

contains

  subroutine grid_init(g, nx, ny, nz, lx, ly, lz, dealias)
    class(grid_t), intent(inout) :: g
    integer,  intent(in) :: nx, ny, nz
    real(rp), intent(in) :: lx, ly, lz
    logical,  intent(in), optional :: dealias
    integer :: i, j, k
    logical :: do_dealias
    real(rp) :: kxc, kyc, kzc

    do_dealias = .true.
    if (present(dealias)) do_dealias = dealias
    if (nx < 2) call fatal('grid: nx must be >= 2')
    if (mod(nx,2) /= 0) call fatal('grid: nx must be even')
    if (ny > 1 .and. mod(ny,2) /= 0) call fatal('grid: ny must be even (or 1)')
    if (nz > 1 .and. mod(nz,2) /= 0) call fatal('grid: nz must be even (or 1)')

    g%nx = nx; g%ny = ny; g%nz = nz
    g%nkx = nx/2 + 1
    g%ndim = 1
    if (ny > 1) g%ndim = 2
    if (nz > 1) g%ndim = 3
    g%lx = lx; g%ly = ly; g%lz = lz
    g%dx = lx/nx; g%dy = ly/ny; g%dz = lz/nz
    g%volume = lx
    if (ny > 1) g%volume = g%volume*ly
    if (nz > 1) g%volume = g%volume*lz

    if (allocated(g%x)) deallocate(g%x, g%y, g%z, g%kx, g%ky, g%kz, g%k2, g%kabs, g%mask)
    allocate(g%x(nx), g%y(ny), g%z(nz))
    do i = 1, nx
      g%x(i) = (i-1)*g%dx
    end do
    do j = 1, ny
      g%y(j) = (j-1)*g%dy
    end do
    do k = 1, nz
      g%z(k) = (k-1)*g%dz
    end do

    allocate(g%kx(g%nkx), g%ky(ny), g%kz(nz))
    do i = 1, g%nkx
      g%kx(i) = twopi*(i-1)/lx
    end do
    g%ky = 0.0_rp; g%kz = 0.0_rp
    if (ny > 1) then
      do j = 1, ny
        g%ky(j) = twopi*wrap(j-1, ny)/ly
      end do
    end if
    if (nz > 1) then
      do k = 1, nz
        g%kz(k) = twopi*wrap(k-1, nz)/lz
      end do
    end if

    allocate(g%k2(g%nkx,ny,nz), g%kabs(g%nkx,ny,nz), g%mask(g%nkx,ny,nz))
    kxc = twopi*(nx/2)/lx*(2.0_rp/3.0_rp)
    kyc = huge(1.0_rp); kzc = huge(1.0_rp)
    if (ny > 1) kyc = twopi*(ny/2)/ly*(2.0_rp/3.0_rp)
    if (nz > 1) kzc = twopi*(nz/2)/lz*(2.0_rp/3.0_rp)
    g%kmax_dealiased = 0.0_rp
    do k = 1, nz
      do j = 1, ny
        do i = 1, g%nkx
          g%k2(i,j,k) = g%kx(i)**2 + g%ky(j)**2 + g%kz(k)**2
          g%kabs(i,j,k) = sqrt(g%k2(i,j,k))
          g%mask(i,j,k) = 1.0_rp
          if (do_dealias) then
            if (abs(g%kx(i)) >= kxc .or. abs(g%ky(j)) >= kyc .or. abs(g%kz(k)) >= kzc) g%mask(i,j,k) = 0.0_rp
          else
            ! always drop the Nyquist modes, whose sign is undefined for real data
            if (i == g%nkx) g%mask(i,j,k) = 0.0_rp
            if (ny > 1 .and. j == ny/2+1) g%mask(i,j,k) = 0.0_rp
            if (nz > 1 .and. k == nz/2+1) g%mask(i,j,k) = 0.0_rp
          end if
          if (g%mask(i,j,k) > 0.0_rp) g%kmax_dealiased = max(g%kmax_dealiased, g%kabs(i,j,k))
        end do
      end do
    end do
  end subroutine grid_init

  !> Map FFT index 0..n-1 to signed integer wavenumber index.
  pure integer function wrap(m, n)
    integer, intent(in) :: m, n
    if (m <= n/2) then
      wrap = m
    else
      wrap = m - n
    end if
  end function wrap

end module hydft_grid

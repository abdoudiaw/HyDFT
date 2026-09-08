! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Hartree (mean-field) energy F_H = (1/2) int int n(r) v(|r-r'|) n(r') with
!> v(k) = 4 pi e^2/(k^2 + kappa^2): Yukawa for ions screened by Boltzmann
!> electrons, Coulomb (kappa = 0) for electrons in a neutralising background
!> (k = 0 mode removed). mu = v * n, K = v(k).
module hydft_fe_hartree
  use hydft_kinds
  use hydft_spectral
  use hydft_free_energy
  implicit none
  private
  public :: fe_hartree_t

  type, extends(fe_term_t) :: fe_hartree_t
    real(rp) :: e2 = 1.0_rp/3.0_rp
    real(rp) :: kappa = 0.0_rp
    real(rp), allocatable :: vk(:,:,:)     !< v(k) on the half spectrum
    real(rp), allocatable :: w(:,:,:)
  contains
    procedure :: setup => hartree_setup
    procedure :: add_mu => hartree_add_mu
    procedure :: kernel => hartree_kernel
    procedure :: energy => hartree_energy
    procedure :: beta_v_k => hartree_beta_v_k
  end type fe_hartree_t

contains

  subroutine hartree_setup(t, spc)
    class(fe_hartree_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    integer :: i, j, k
    t%name = 'hartree'
    if (allocated(t%vk)) deallocate(t%vk, t%w)
    allocate(t%vk(spc%g%nkx, spc%g%ny, spc%g%nz), t%w(spc%g%nx, spc%g%ny, spc%g%nz))
    do k = 1, spc%g%nz; do j = 1, spc%g%ny; do i = 1, spc%g%nkx
      t%vk(i,j,k) = t%kernel(spc%g%kabs(i,j,k))
    end do; end do; end do
  end subroutine hartree_setup

  subroutine hartree_add_mu(t, spc, n, nhat, mu)
    class(fe_hartree_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    real(rp), intent(inout) :: mu(:,:,:)
    call spc%convolve(t%vk, nhat, t%w)
    mu = mu + t%w
  end subroutine hartree_add_mu

  real(rp) function hartree_kernel(t, k) result(kk)
    class(fe_hartree_t), intent(in) :: t
    real(rp), intent(in) :: k
    if (k*k + t%kappa**2 > 0.0_rp) then
      kk = 4.0_rp*pi*t%e2/(k*k + t%kappa**2)
    else
      kk = 0.0_rp     ! neutralising background
    end if
  end function hartree_kernel

  !> beta v(k), the mean-field part of the direct correlation function (with sign reversed).
  real(rp) function hartree_beta_v_k(t, k) result(bv)
    class(fe_hartree_t), intent(in) :: t
    real(rp), intent(in) :: k
    bv = t%beta*t%kernel(k)
  end function hartree_beta_v_k

  real(rp) function hartree_energy(t, spc, n, nhat) result(f)
    class(fe_hartree_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    call spc%convolve(t%vk, nhat, t%w)
    f = 0.5_rp*spc%mean(n*t%w)*spc%g%volume
  end function hartree_energy

end module hydft_fe_hartree

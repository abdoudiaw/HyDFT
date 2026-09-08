! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Ramakrishnan-Yussouff correlation free energy (PRE Eq. 39, Sci. Rep. Eq. 13):
!>   F_cor = F_cor[n0] - (2 beta)^-1 int int dn(r) c_cor(|r-r'|) dn(r'),  dn = n - n0,
!>   mu = -beta^-1 c_cor * dn,   K = -c_cor(k)/beta.
!> c_cor(k) = c(k) + beta v(k) when an explicit Hartree term is present (so the
!> mean field is not double counted), otherwise the full c(k).
module hydft_fe_correlation_ry
  use hydft_kinds
  use hydft_spectral
  use hydft_free_energy
  use hydft_structure_source
  implicit none
  private
  public :: fe_ry_t

  type, extends(fe_term_t) :: fe_ry_t
    type(structure_t) :: st
    logical  :: subtract_mean_field = .true.
    real(rp) :: e2 = 1.0_rp/3.0_rp
    real(rp) :: kappa = 0.0_rp
    real(rp), allocatable :: ck(:,:,:)    !< c_cor(k) on the half spectrum
    real(rp), allocatable :: w(:,:,:)
  contains
    procedure :: setup => ry_setup
    procedure :: add_mu => ry_add_mu
    procedure :: kernel => ry_kernel
    procedure :: energy => ry_energy
    procedure :: c_cor => ry_c_cor
  end type fe_ry_t

contains

  !> c_cor(k), interpolated from the structure table.
  real(rp) function ry_c_cor(t, k) result(c)
    class(fe_ry_t), intent(in) :: t
    real(rp), intent(in) :: k
    real(rp) :: k1
    ! c(k) + beta v(k) is smooth at k -> 0 even for Coulomb; evaluate the table
    ! at k and add beta v at the same k. Below the first table point use the
    ! first table point (the sum is regular there).
    k1 = max(k, t%st%k(1))
    c = t%st%c_of_k(k1)
    if (t%subtract_mean_field) c = c + t%beta*4.0_rp*pi*t%e2/(k1*k1 + t%kappa**2)
  end function ry_c_cor

  subroutine ry_setup(t, spc)
    class(fe_ry_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    integer :: i, j, k
    t%name = 'ry'
    if (.not. allocated(t%st%k)) call fatal('ry term: no structure table')
    if (allocated(t%ck)) deallocate(t%ck, t%w)
    allocate(t%ck(spc%g%nkx, spc%g%ny, spc%g%nz), t%w(spc%g%nx, spc%g%ny, spc%g%nz))
    do k = 1, spc%g%nz; do j = 1, spc%g%ny; do i = 1, spc%g%nkx
      t%ck(i,j,k) = t%c_cor(spc%g%kabs(i,j,k))
    end do; end do; end do
  end subroutine ry_setup

  subroutine ry_add_mu(t, spc, n, nhat, mu)
    class(fe_ry_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    real(rp), intent(inout) :: mu(:,:,:)
    ! the k = 0 mode of dn = n - n0 vanishes for a mass-conserving run; drop it explicitly
    call spc%convolve(t%ck, nhat_no_mean(nhat), t%w)
    mu = mu - t%w/t%beta
  contains
    function nhat_no_mean(a) result(b)
      complex(rp), intent(in) :: a(:,:,:)
      complex(rp) :: b(size(a,1), size(a,2), size(a,3))
      b = a
      b(1,1,1) = a(1,1,1) - t%n0
    end function nhat_no_mean
  end subroutine ry_add_mu

  real(rp) function ry_kernel(t, k) result(kk)
    class(fe_ry_t), intent(in) :: t
    real(rp), intent(in) :: k
    kk = -t%c_cor(k)/t%beta
  end function ry_kernel

  real(rp) function ry_energy(t, spc, n, nhat) result(f)
    class(fe_ry_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    complex(rp) :: b(size(nhat,1), size(nhat,2), size(nhat,3))
    b = nhat
    b(1,1,1) = nhat(1,1,1) - t%n0
    call spc%convolve(t%ck, b, t%w)
    f = -0.5_rp*spc%mean((n - t%n0)*t%w)*spc%g%volume/t%beta
  end function ry_energy

end module hydft_fe_correlation_ry

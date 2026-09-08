! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Non-interacting free energy of the finite-temperature electron gas.
!>
!> fe_tf_t: Thomas-Fermi, T_TF = beta^-1 int [n alpha - (2/3) cn I_{3/2}(alpha)],
!>   n = cn I_{1/2}(alpha), cn = sqrt2 m^{3/2}/(pi^2 hbar^3 beta^{3/2});
!>   mu = alpha(n)/beta, K = 1/(beta dn/dalpha) = 2/(beta cn I_{-1/2}(alpha0)).
!>
!> fe_gradient_t: Kirzhnits/von Weizsaecker gradient correction (Sci. Rep. Eq. 17),
!>   T_2 = C int xi(alpha) |grad n|^2,  C = gamma_g (3 sqrt2 pi^2/8) hbar^5 beta^{3/2}/m^{5/2},
!>   xi = I'_{-1/2}/I_{-1/2}^2. Functional derivative (re-derived; gamma_g = 1 gives
!>   exactly the Bohm potential -(hbar^2/2m) lap(sqrt n)/sqrt n at T = 0):
!>   mu_2 = -C [ (dxi/dn) |grad n|^2 + 2 xi lap n ],   K_2 = 2 C xi(alpha0) k^2.
module hydft_fe_kinetic_tfk
  use hydft_kinds
  use hydft_spectral
  use hydft_free_energy
  use hydft_fermi_dirac
  implicit none
  private
  public :: fe_tf_t, fe_gradient_t

  type, extends(fe_term_t) :: fe_tf_t
    real(rp) :: cn = 1.0_rp
    real(rp) :: alpha0 = 0.0_rp
    real(rp), allocatable :: alpha(:,:,:)   !< work: alpha(n) on the grid
  contains
    procedure :: setup => tf_setup
    procedure :: add_mu => tf_add_mu
    procedure :: kernel => tf_kernel
    procedure :: energy => tf_energy
    procedure :: alpha_of_n => tf_alpha_of_n
  end type fe_tf_t

  type, extends(fe_term_t) :: fe_gradient_t
    real(rp) :: cn = 1.0_rp
    real(rp) :: alpha0 = 0.0_rp
    real(rp) :: cg = 0.0_rp        !< C above
    real(rp) :: gamma_g = 1.0_rp/9.0_rp
    real(rp), allocatable :: lapn(:,:,:), gn(:,:,:,:)
  contains
    procedure :: setup => grad_setup
    procedure :: add_mu => grad_add_mu
    procedure :: kernel => grad_kernel
    procedure :: energy => grad_energy
  end type fe_gradient_t

contains

  ! ------------------------------------------------------------------ TF
  subroutine tf_setup(t, spc)
    class(fe_tf_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    t%name = 'tf'
    if (allocated(t%alpha)) deallocate(t%alpha)
    allocate(t%alpha(spc%g%nx, spc%g%ny, spc%g%nz))
  end subroutine tf_setup

  subroutine tf_alpha_of_n(t, n, alpha)
    class(fe_tf_t), intent(in) :: t
    real(rp), intent(in) :: n(:,:,:)
    real(rp), intent(out) :: alpha(:,:,:)
    integer :: i, j, k
    if (minval(n) <= 0.0_rp) call fatal('tf term: density is not positive')
    do k = 1, size(n,3); do j = 1, size(n,2); do i = 1, size(n,1)
      alpha(i,j,k) = fd_alpha_from_ip12(n(i,j,k)/t%cn)
    end do; end do; end do
  end subroutine tf_alpha_of_n

  subroutine tf_add_mu(t, spc, n, nhat, mu)
    class(fe_tf_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    real(rp), intent(inout) :: mu(:,:,:)
    call t%alpha_of_n(n, t%alpha)
    mu = mu + t%alpha/t%beta
  end subroutine tf_add_mu

  real(rp) function tf_kernel(t, k) result(kk)
    class(fe_tf_t), intent(in) :: t
    real(rp), intent(in) :: k
    kk = 2.0_rp/(t%beta*t%cn*fd_im12(t%alpha0))
  end function tf_kernel

  real(rp) function tf_energy(t, spc, n, nhat) result(f)
    class(fe_tf_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    integer :: i, j, k
    real(rp) :: s
    call t%alpha_of_n(n, t%alpha)
    s = 0.0_rp
    do k = 1, size(n,3); do j = 1, size(n,2); do i = 1, size(n,1)
      s = s + n(i,j,k)*t%alpha(i,j,k) - (2.0_rp/3.0_rp)*t%cn*fd_ip32(t%alpha(i,j,k))
    end do; end do; end do
    f = s/real(size(n), rp)*spc%g%volume/t%beta
  end function tf_energy

  ! ------------------------------------------------------------ gradient
  subroutine grad_setup(t, spc)
    class(fe_gradient_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    t%name = 'gradient'
    if (allocated(t%lapn)) deallocate(t%lapn, t%gn)
    allocate(t%lapn(spc%g%nx, spc%g%ny, spc%g%nz), t%gn(spc%g%nx, spc%g%ny, spc%g%nz, 3))
  end subroutine grad_setup

  subroutine grad_add_mu(t, spc, n, nhat, mu)
    class(fe_gradient_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    real(rp), intent(inout) :: mu(:,:,:)
    integer :: i, j, k, d
    real(rp) :: a, xi, dxidn, g2
    if (minval(n) <= 0.0_rp) call fatal('gradient term: density is not positive')
    do d = 1, 3
      if (d <= spc%g%ndim) then
        call spc%grad_hat(nhat, d, t%gn(:,:,:,d))
      else
        t%gn(:,:,:,d) = 0.0_rp
      end if
    end do
    call spc%convolve(-spc%g%k2, nhat, t%lapn)
    do k = 1, size(n,3); do j = 1, size(n,2); do i = 1, size(n,1)
      a = fd_alpha_from_ip12(n(i,j,k)/t%cn)
      xi = fd_xi(a)
      dxidn = fd_dxi(a)/(0.5_rp*t%cn*fd_im12(a))      ! dxi/dalpha * dalpha/dn
      g2 = t%gn(i,j,k,1)**2 + t%gn(i,j,k,2)**2 + t%gn(i,j,k,3)**2
      mu(i,j,k) = mu(i,j,k) - t%cg*(dxidn*g2 + 2.0_rp*xi*t%lapn(i,j,k))
    end do; end do; end do
  end subroutine grad_add_mu

  real(rp) function grad_kernel(t, k) result(kk)
    class(fe_gradient_t), intent(in) :: t
    real(rp), intent(in) :: k
    kk = 2.0_rp*t%cg*fd_xi(t%alpha0)*k*k
  end function grad_kernel

  real(rp) function grad_energy(t, spc, n, nhat) result(f)
    class(fe_gradient_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    integer :: i, j, k, d
    real(rp) :: s, a
    do d = 1, 3
      if (d <= spc%g%ndim) then
        call spc%grad_hat(nhat, d, t%gn(:,:,:,d))
      else
        t%gn(:,:,:,d) = 0.0_rp
      end if
    end do
    s = 0.0_rp
    do k = 1, size(n,3); do j = 1, size(n,2); do i = 1, size(n,1)
      a = fd_alpha_from_ip12(n(i,j,k)/t%cn)
      s = s + fd_xi(a)*(t%gn(i,j,k,1)**2 + t%gn(i,j,k,2)**2 + t%gn(i,j,k,3)**2)
    end do; end do; end do
    f = t%cg*s/real(size(n), rp)*spc%g%volume
  end function grad_energy

end module hydft_fe_kinetic_tfk

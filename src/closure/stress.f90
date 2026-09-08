! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Stress closures for the dissipative (off-diagonal) stress Pi_ij.
!>   newtonian : Pi = Pi0(grad u) = eta (grad u + grad u^T - 2/3 div u I) + xi div u I
!>   maxwell   : Pi is a state variable, D_t Pi = (Pi0 - Pi)/tau
!> Components are stored as (xx, yy, zz, xy, xz, yz). Viscosities are in
!> reduced units of m n0 wp a^2, i.e. the code value is eta_bar * n0.
module hydft_stress
  use hydft_kinds
  use hydft_utils, only: lower
  use hydft_spectral
  implicit none
  private
  public :: stress_t, ncomp

  integer, parameter :: ncomp = 6

  type :: stress_t
    character(len=16) :: name = 'newtonian'
    logical  :: evolve = .false.      !< true for maxwell
    real(rp) :: tau = 0.0_rp
    real(rp) :: eta = 0.0_rp          !< code units (n0 * eta_bar)
    real(rp) :: xi = 0.0_rp
    real(rp), allocatable :: gu(:,:,:,:,:)   !< work: du_i/dx_j
  contains
    procedure :: setup => stress_setup
    procedure :: pi0 => stress_pi0
    procedure :: divergence => stress_divergence
  end type stress_t

contains

  subroutine stress_setup(s, name, eta_bar, xi_bar, tau, n0, spc)
    class(stress_t), intent(inout) :: s
    character(len=*), intent(in) :: name
    real(rp), intent(in) :: eta_bar, xi_bar, tau, n0
    type(spec_t), intent(in) :: spc
    s%name = lower(name)
    select case (trim(s%name))
    case ('newtonian')
      s%evolve = .false.; s%tau = 0.0_rp
    case ('maxwell')
      s%evolve = tau > 0.0_rp
      s%tau = tau
      if (.not. s%evolve) call warn('maxwell closure with tau = 0 reduces to newtonian')
    case default
      call fatal('unknown stress closure: '//trim(name))
    end select
    s%eta = eta_bar*n0
    s%xi = xi_bar*n0
    if (allocated(s%gu)) deallocate(s%gu)
    allocate(s%gu(spc%g%nx, spc%g%ny, spc%g%nz, 3, 3))
  end subroutine stress_setup

  !> Newtonian stress from the velocity field u(:,:,:,1:3).
  subroutine stress_pi0(s, spc, u, pi)
    class(stress_t), intent(inout) :: s
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: u(:,:,:,:)
    real(rp), intent(out) :: pi(:,:,:,:)
    integer :: i, j
    real(rp), allocatable :: divu(:,:,:)
    do i = 1, 3
      if (i <= spc%g%ndim) then
        call spc%grad(u(:,:,:,i), s%gu(:,:,:,i,:))
      else
        s%gu(:,:,:,i,:) = 0.0_rp
      end if
    end do
    allocate(divu(spc%g%nx, spc%g%ny, spc%g%nz))
    divu = s%gu(:,:,:,1,1) + s%gu(:,:,:,2,2) + s%gu(:,:,:,3,3)
    do i = 1, 3
      pi(:,:,:,i) = s%eta*(2.0_rp*s%gu(:,:,:,i,i) - (2.0_rp/3.0_rp)*divu) + s%xi*divu
    end do
    j = 4
    pi(:,:,:,4) = s%eta*(s%gu(:,:,:,1,2) + s%gu(:,:,:,2,1))
    pi(:,:,:,5) = s%eta*(s%gu(:,:,:,1,3) + s%gu(:,:,:,3,1))
    pi(:,:,:,6) = s%eta*(s%gu(:,:,:,2,3) + s%gu(:,:,:,3,2))
  end subroutine stress_pi0

  !> f_i = d_j Pi_ij
  subroutine stress_divergence(s, spc, pi, f)
    class(stress_t), intent(inout) :: s
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: pi(:,:,:,:)
    real(rp), intent(out) :: f(:,:,:,:)
    real(rp), allocatable :: v(:,:,:,:)
    integer, parameter :: idx(3,3) = reshape([1, 4, 5, 4, 2, 6, 5, 6, 3], [3, 3])
    integer :: i, j
    allocate(v(spc%g%nx, spc%g%ny, spc%g%nz, 3))
    do i = 1, 3
      if (i > spc%g%ndim) then
        f(:,:,:,i) = 0.0_rp
        cycle
      end if
      do j = 1, 3
        v(:,:,:,j) = pi(:,:,:,idx(i,j))
      end do
      call spc%div(v, f(:,:,:,i))
    end do
  end subroutine stress_divergence

end module hydft_stress

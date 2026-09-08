! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Abstract free-energy term. Each term provides
!>   add_mu   : adds its functional derivative mu_i[n] = dF_i/dn on the grid,
!>   kernel   : its second functional derivative K_i(k) = d^2F_i/dn dn |_{n0}
!>              about the uniform state, an isotropic function of |k|,
!>   energy   : F_i[n] (volume integral),
!>   setup    : one-time precomputation on the grid (kernel arrays).
!> The nonlinear solver uses add_mu; the linear-response tool uses kernel.
!> Test test_kernels checks that the two agree.
module hydft_free_energy
  use hydft_kinds
  use hydft_spectral
  implicit none
  private
  public :: fe_term_t, fe_box_t

  type, abstract :: fe_term_t
    character(len=16) :: name = ''
    real(rp) :: n0 = 3.0_rp/(4.0_rp*pi)
    real(rp) :: beta = 3.0_rp
  contains
    procedure(setup_iface),  deferred :: setup
    procedure(add_mu_iface), deferred :: add_mu
    procedure(kernel_iface), deferred :: kernel
    procedure(energy_iface), deferred :: energy
  end type fe_term_t

  abstract interface
    subroutine setup_iface(t, spc)
      import :: fe_term_t, spec_t
      class(fe_term_t), intent(inout) :: t
      type(spec_t), intent(inout) :: spc
    end subroutine setup_iface
    subroutine add_mu_iface(t, spc, n, nhat, mu)
      import :: fe_term_t, spec_t, rp
      class(fe_term_t), intent(inout) :: t
      type(spec_t), intent(inout) :: spc
      real(rp), intent(in) :: n(:,:,:)
      complex(rp), intent(in) :: nhat(:,:,:)
      real(rp), intent(inout) :: mu(:,:,:)
    end subroutine add_mu_iface
    real(rp) function kernel_iface(t, k)
      import :: fe_term_t, rp
      class(fe_term_t), intent(in) :: t
      real(rp), intent(in) :: k
    end function kernel_iface
    real(rp) function energy_iface(t, spc, n, nhat)
      import :: fe_term_t, spec_t, rp
      class(fe_term_t), intent(inout) :: t
      type(spec_t), intent(inout) :: spc
      real(rp), intent(in) :: n(:,:,:)
      complex(rp), intent(in) :: nhat(:,:,:)
    end function energy_iface
  end interface

  !> Container so that arrays of polymorphic terms can be built.
  type :: fe_box_t
    class(fe_term_t), allocatable :: t
  end type fe_box_t

end module hydft_free_energy

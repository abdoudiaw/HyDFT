! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Classical ideal-gas free energy F_id = beta^-1 int n (ln n - 1) (thermal
!> wavelength constant dropped): mu = beta^-1 ln n, K = 1/(beta n0).
module hydft_fe_ideal_classical
  use hydft_kinds
  use hydft_spectral
  use hydft_free_energy
  implicit none
  private
  public :: fe_ideal_t

  type, extends(fe_term_t) :: fe_ideal_t
  contains
    procedure :: setup => ideal_setup
    procedure :: add_mu => ideal_add_mu
    procedure :: kernel => ideal_kernel
    procedure :: energy => ideal_energy
  end type fe_ideal_t

contains

  subroutine ideal_setup(t, spc)
    class(fe_ideal_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    t%name = 'ideal'
  end subroutine ideal_setup

  subroutine ideal_add_mu(t, spc, n, nhat, mu)
    class(fe_ideal_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    real(rp), intent(inout) :: mu(:,:,:)
    if (minval(n) <= 0.0_rp) call fatal('ideal term: density is not positive')
    mu = mu + log(n)/t%beta
  end subroutine ideal_add_mu

  real(rp) function ideal_kernel(t, k) result(kk)
    class(fe_ideal_t), intent(in) :: t
    real(rp), intent(in) :: k
    kk = 1.0_rp/(t%beta*t%n0)
  end function ideal_kernel

  real(rp) function ideal_energy(t, spc, n, nhat) result(f)
    class(fe_ideal_t), intent(inout) :: t
    type(spec_t), intent(inout) :: spc
    real(rp), intent(in) :: n(:,:,:)
    complex(rp), intent(in) :: nhat(:,:,:)
    f = spc%mean(n*(log(n) - 1.0_rp))*spc%g%volume/t%beta
  end function ideal_energy

end module hydft_fe_ideal_classical

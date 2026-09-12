! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Budgets (mass, momentum, kinetic + free energy) and spectral probes.
module hydft_diagnostics
  use hydft_kinds
  use hydft_utils
  use hydft_model
  implicit none
  private
  public :: diag_t

  type :: diag_t
    integer :: udiag = -1, uprobe = -1
    integer :: probe_mode(3) = [1, 0, 0]
    real(rp) :: mass0 = 0.0_rp
  contains
    procedure :: open => diag_open
    procedure :: write => diag_write
    procedure :: probe => diag_probe
    procedure :: close => diag_close
    procedure :: budgets => diag_budgets
  end type diag_t

contains

  subroutine diag_open(d, prefix, probe_mode)
    class(diag_t), intent(inout) :: d
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: probe_mode(3)
    d%probe_mode = probe_mode
    d%udiag = open_file(trim(prefix)//'_diag.dat', 'write')
    write(d%udiag,'(a)') '# step   t   dt   mass   momx   momy   momz   E_kin   F   E_tot   n_min   n_max'
    d%uprobe = open_file(trim(prefix)//'_probe.dat', 'write')
    write(d%uprobe,'(a,3i4,a)') '# probe mode', probe_mode, ':   t   Re(dn_k/n0)   Im(dn_k/n0)   |dn_k/n0|   Re(mom_k)   Im(mom_k)'
  end subroutine diag_open

  subroutine diag_budgets(d, m, mass, mom, ekin, fe)
    class(diag_t), intent(inout) :: d
    type(model_t), intent(inout) :: m
    real(rp), intent(out) :: mass, mom(3), ekin, fe
    integer :: i
    real(rp) :: v
    v = m%spc%g%volume
    mass = m%spc%mean(m%sp%n)*v
    do i = 1, 3
      mom(i) = m%spc%mean(m%sp%mom(:,:,:,i))*v
    end do
    ekin = 0.0_rp
    do i = 1, m%spc%g%ndim
      ekin = ekin + m%spc%mean(m%sp%mom(:,:,:,i)**2/m%sp%n)
    end do
    ekin = 0.5_rp*ekin*v/m%sp%mass
    call m%spc%fwd(m%sp%n, m%nhat)
    fe = m%sp%fun%energy(m%spc, m%sp%n, m%nhat)
    if (m%ext%active) fe = fe + m%spc%mean(m%sp%n*m%ext%v)*v
  end subroutine diag_budgets

  subroutine diag_write(d, m)
    class(diag_t), intent(inout) :: d
    type(model_t), intent(inout) :: m
    real(rp) :: mass, mom(3), ekin, fe
    call d%budgets(m, mass, mom, ekin, fe)
    write(d%udiag,'(i9,11es18.9e3)') m%step, m%t, m%dt, mass, mom, ekin, fe, ekin + fe, minval(m%sp%n), maxval(m%sp%n)
  end subroutine diag_write

  subroutine diag_probe(d, m)
    class(diag_t), intent(inout) :: d
    type(model_t), intent(inout) :: m
    complex(rp) :: a, b
    a = m%spc%mode(m%sp%n, d%probe_mode(1), d%probe_mode(2), d%probe_mode(3))/m%sp%un%n0
    b = m%spc%mode(m%sp%mom(:,:,:,1), d%probe_mode(1), d%probe_mode(2), d%probe_mode(3))
    write(d%uprobe,'(6es18.9e3)') m%t, real(a), aimag(a), abs(a), real(b), aimag(b)
  end subroutine diag_probe

  subroutine diag_close(d)
    class(diag_t), intent(inout) :: d
    if (d%udiag >= 0) close(d%udiag)
    if (d%uprobe >= 0) close(d%uprobe)
    d%udiag = -1; d%uprobe = -1
  end subroutine diag_close

end module hydft_diagnostics

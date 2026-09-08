! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Mass and momentum are conserved to round-off; the total energy of a
!> Newtonian run does not increase.
program test_conservation
  use hydft_kinds
  use hydft_params
  use hydft_model
  use hydft_initial
  use hydft_diagnostics
  implicit none
  integer :: nfail
  nfail = 0
  call run_case('newtonian', 2.0_rp, 0.0_rp)
  call run_case('maxwell', 0.2_rp, 2.0_rp)
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_conservation: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_conservation: all passed'
contains
  subroutine run_case(closure, eta, tau)
    character(len=*), intent(in) :: closure
    real(rp), intent(in) :: eta, tau
    type(params_t) :: p
    type(model_t) :: m
    type(diag_t) :: d
    real(rp) :: mass0, mom0(3), ekin0, fe0, mass, mom(3), ekin, fe, err
    integer :: i
    p%kind = 'ion'; p%gamma = 2.0_rp; p%kappa = 0.1_rp
    p%nx = 32; p%ny = 4; p%nz = 1; p%lx = 20.0_rp; p%ly = 5.0_rp
    p%hnc_nr = 4096; p%hnc_rmax = 60.0_rp
    p%transport = 'constant'; p%eta = eta; p%closure = closure; p%tau_model = 'constant'; p%tau = tau
    p%cfl = 0.3_rp
    call m%init(p, verbose=.false.)
    call set_initial_condition(m%spc, m%sp, 'random', 0.02_rp, [1,0,0], 1.0_rp, 7)
    call d%budgets(m, mass0, mom0, ekin0, fe0)
    do i = 1, 100
      call m%compute_dt()
      call m%advance()
    end do
    call d%budgets(m, mass, mom, ekin, fe)
    write(*,'(a,a,a,f8.3,a,es10.2)') '  case ', closure, '  t =', m%t, '  dt =', m%dt
    err = abs(mass/mass0 - 1.0_rp); call check(err < 1e-12_rp, 'mass conserved', err)
    err = maxval(abs(mom - mom0))/(mass0*sqrt(ekin/mass0)); call check(err < 1e-12_rp, 'momentum conserved', err)
    call check(ekin > 0.0_rp, 'kinetic energy developed', ekin)
    if (closure == 'newtonian') then
      err = (ekin + fe - ekin0 - fe0)/abs(ekin0 + fe0)
      call check(err < 1e-8_rp, 'total energy non-increasing', err)
      write(*,'(a,es12.4)') '  info  relative energy change =', err
    end if
  end subroutine run_case
  subroutine check(ok, name, val)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: name
    real(rp), intent(in) :: val
    if (ok) then
      write(*,'(a,a,es12.3)') '    pass  ', name, val
    else
      write(*,'(a,a,es12.3)') '    FAIL  ', name, val
      nfail = nfail + 1
    end if
  end subroutine check
end program test_conservation

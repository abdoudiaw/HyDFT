! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Driven mode: the measured dn(k,w)/dv_ext(k,w) matches chi(k,w) from the
!> linear-response module (Newtonian and Maxwell).
program test_driven_mode
  use hydft_kinds
  use hydft_params
  use hydft_model
  use hydft_initial
  use hydft_linear_response
  use hydft_analysis
  implicit none
  integer :: nfail
  nfail = 0
  call run_case('newtonian', 0.0_rp, 0.7_rp)
  call run_case('maxwell', 1.0_rp, 0.7_rp)
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_driven_mode: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_driven_mode: all passed'
contains
  subroutine run_case(closure, tau, wd_over_w0)
    character(len=*), intent(in) :: closure
    real(rp), intent(in) :: tau, wd_over_w0
    type(params_t) :: p
    type(model_t) :: m
    type(linear_t) :: lin
    real(rp), allocatable :: ts(:)
    complex(rp), allocatable :: zs(:)
    real(rp) :: q, wd, err, scatter, amp
    complex(rp) :: chi_lin, chi_sim, a
    integer :: n, mode
    mode = 2
    amp = 1.0e-5_rp
    p%kind = 'ion'; p%gamma = 2.0_rp; p%kappa = 0.1_rp
    p%nx = 32; p%lx = 20.0_rp; p%hnc_nr = 4096; p%hnc_rmax = 60.0_rp
    p%transport = 'constant'; p%eta = 0.3_rp; p%closure = closure; p%tau_model = 'constant'; p%tau = tau
    p%cfl = 0.1_rp
    q = twopi*mode/p%lx
    ! frequency of the drive relative to the undamped mode frequency
    p%vext = 'driven'; p%vext_amplitude = amp; p%vext_mode = [mode,0,0]; p%vext_ramp = 10.0_rp
    call m%init(p, verbose=.false.)
    call lin%init(m%sp)
    wd = wd_over_w0*q*sqrt(lin%c2(q))
    m%ext%omega = wd
    call set_initial_condition(m%spc, m%sp, 'uniform', 0.0_rp, [1,0,0], 1.0_rp, 1)
    allocate(ts(200000), zs(200000))
    n = 0
    do while (m%t < 300.0_rp)
      call m%compute_dt()
      call m%advance()
      n = n + 1
      ts(n) = m%t
      zs(n) = m%spc%mode(m%sp%n, mode, 0, 0)
    end do
    ! drive A cos(kx - w t) has +k coefficient (A/2) e^{-i w t}
    call fit_driven_response(ts(1:n), zs(1:n), wd, 250.0_rp, a, scatter)
    chi_sim = 2.0_rp*a/amp
    chi_lin = lin%chi(q, wd)
    write(*,'(a,a,a,f7.4,a,es10.2)') '  case ', closure, '  w_drive =', wd, '  transient scatter =', scatter
    write(*,'(a,2es14.6,a,2es14.6)') '    chi linear =', chi_lin, '   simulation =', chi_sim
    err = abs(chi_sim - chi_lin)/abs(chi_lin)
    call check(err < 0.02_rp, 'chi(k,w) from driven run', err)
    call check(scatter < 0.05_rp, 'transient decayed', scatter)
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
end program test_driven_mode

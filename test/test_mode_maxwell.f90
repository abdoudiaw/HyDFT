! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Single-mode simulation vs linear response with the Maxwell closure:
!> elastic regime (w tau >> 1) and the tau -> 0 limit converging to Newtonian.
program test_mode_maxwell
  use hydft_kinds
  use hydft_params
  use hydft_model
  use hydft_initial
  use hydft_linear_response
  use hydft_analysis
  implicit none
  integer :: nfail
  real(rp) :: om_n, ga_n, om_s, ga_s
  nfail = 0
  call run_case('newtonian', 0.1_rp, 0.0_rp, om_n, ga_n)
  call run_case('maxwell', 0.1_rp, 3.0_rp, om_s, ga_s)
  call run_case('maxwell', 0.1_rp, 0.02_rp, om_s, ga_s)
  write(*,'(a,2f10.5,a,2f10.5)') '  info  small-tau maxwell (w,g) =', om_s, ga_s, '   newtonian =', om_n, ga_n
  call check(abs(om_s/om_n - 1.0_rp) < 0.01_rp .and. abs(ga_s/ga_n - 1.0_rp) < 0.05_rp, 'tau -> 0 converges to Newtonian', abs(ga_s/ga_n - 1.0_rp))
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_mode_maxwell: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_mode_maxwell: all passed'
contains
  subroutine run_case(closure, eta, tau, om, ga)
    character(len=*), intent(in) :: closure
    real(rp), intent(in) :: eta, tau
    real(rp), intent(out) :: om, ga
    type(params_t) :: p
    type(model_t) :: m
    type(linear_t) :: lin
    real(rp), allocatable :: ts(:), ys(:)
    real(rp) :: q, err
    complex(rp) :: w
    integer :: n, nz, mode
    mode = 2
    p%kind = 'ion'; p%gamma = 2.0_rp; p%kappa = 0.1_rp
    p%nx = 32; p%lx = 20.0_rp; p%hnc_nr = 4096; p%hnc_rmax = 60.0_rp
    p%transport = 'constant'; p%eta = eta; p%closure = closure; p%tau_model = 'constant'; p%tau = tau
    p%cfl = 0.1_rp
    call m%init(p, verbose=.false.)
    call set_initial_condition(m%spc, m%sp, 'mode', 1.0e-4_rp, [mode,0,0], 1.0_rp, 1)
    call lin%init(m%sp)
    q = twopi*mode/p%lx
    w = lin%dispersion(q)
    allocate(ts(200000), ys(200000))
    n = 0
    do while (m%t < 120.0_rp)
      call m%compute_dt()
      call m%advance()
      n = n + 1
      ts(n) = m%t
      ys(n) = real(m%spc%mode(m%sp%n, mode, 0, 0))
    end do
    call fit_damped_mode(ts(1:n), ys(1:n), 20.0_rp, om, ga, nz)
    write(*,'(a,a,a,f6.3,a,f6.2,a,i0,a,es10.2)') '  case ', closure, '  tau =', tau, '  w tau =', real(w)*tau, '  steps = ', n, '  dt =', m%dt
    write(*,'(a,2f12.6,a,2f12.6,a,i0)') '    linear (Re,Im) =', real(w), aimag(w), '   simulation =', om, ga, '   zeros = ', nz
    err = abs(om/real(w) - 1.0_rp); call check(err < 0.01_rp, 'frequency', err)
    err = abs(ga/aimag(w) - 1.0_rp); call check(err < 0.05_rp, 'damping', err)
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
end program test_mode_maxwell

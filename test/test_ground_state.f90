! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> DDFT closure: with a static external potential and viscous damping the
!> flow relaxes to the density-functional ground state mu[n] + v_ext = const,
!> and the linear part matches dn_k = -dv_k/K(k).
program test_ground_state
  use hydft_kinds
  use hydft_params
  use hydft_model
  use hydft_initial
  use hydft_linear_response
  implicit none
  integer :: nfail
  nfail = 0
  call run_case('ion')
  call run_case('electron')
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_ground_state: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_ground_state: all passed'
contains
  subroutine run_case(kind)
    character(len=*), intent(in) :: kind
    type(params_t) :: p
    type(model_t) :: m
    type(linear_t) :: lin
    real(rp), allocatable :: tot(:,:,:)
    real(rp) :: q, amp, err, ekin, spread
    complex(rp) :: dnk, pred
    integer :: mode
    mode = 1
    p%kind = kind; p%gamma = 2.0_rp; p%kappa = 0.1_rp; p%rs = 1.86_rp
    if (kind == 'electron') p%gamma = 1.0_rp
    p%nx = 32; p%lx = 20.0_rp; p%hnc_nr = 4096; p%hnc_rmax = 60.0_rp
    p%transport = 'constant'; p%eta = 1.0_rp; p%closure = 'newtonian'
    p%cfl = 0.2_rp
    amp = 1.0e-3_rp
    p%vext = 'mode'; p%vext_amplitude = amp; p%vext_mode = [mode,0,0]
    call m%init(p, verbose=.false.)
    call lin%init(m%sp)
    call set_initial_condition(m%spc, m%sp, 'uniform', 0.0_rp, [1,0,0], 1.0_rp, 1)
    do while (m%t < 150.0_rp)
      call m%compute_dt()
      call m%advance()
    end do
    ekin = m%spc%mean(sum(m%sp%mom**2, dim=4)/m%sp%n)
    call m%spc%fwd(m%sp%n, m%nhat)
    call m%sp%fun%mu(m%spc, m%sp%n, m%nhat, m%mu)
    allocate(tot(p%nx,1,1))
    tot = m%mu + m%ext%v
    spread = (maxval(tot) - minval(tot))/amp
    q = twopi*mode/p%lx
    dnk = m%spc%mode(m%sp%n, mode, 0, 0)
    pred = -0.5_rp*amp/m%sp%fun%kernel(q)
    write(*,'(a,a,a,es10.2,a,es10.2)') '  case ', trim(kind), '  residual kinetic energy =', ekin, '  spread of mu+v_ext / A =', spread
    call check(spread < 2.0e-3_rp, 'mu[n] + v_ext uniform', spread)
    err = abs(dnk - pred)/abs(pred)
    call check(err < 1.0e-2_rp, 'dn_k = -dv_k/K(k)', err)
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
end program test_ground_state

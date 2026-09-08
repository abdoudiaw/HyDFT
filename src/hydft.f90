! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Top-level driver: read the input file, build the model, run the time loop.
module hydft
  use hydft_kinds
  use hydft_params
  use hydft_model
  use hydft_initial
  use hydft_diagnostics
  use hydft_io
  implicit none
  private
  public :: hydft_run

contains

  subroutine hydft_run(infile)
    character(len=*), intent(in) :: infile
    type(params_t) :: p
    type(model_t) :: m
    type(diag_t) :: d
    integer :: irecompute
    real(rp) :: t_wall0, t_wall1

    call p%read(infile)
    call p%print()
    call m%init(p, verbose=.true.)
    call set_initial_condition(m%spc, m%sp, p%ic, p%ic_amplitude, p%ic_mode, p%ic_width, p%ic_seed)
    call d%open(p%prefix, p%probe_mode)
    call m%compute_dt()
    write(*,'(a,es12.4)') '  initial dt =', m%dt
    call d%write(m)
    call d%probe(m)
    if (p%iout > 0) call write_fields(m, p%prefix)
    call cpu_time(t_wall0)
    irecompute = 0
    do while (m%t < p%tend .and. m%step < p%nsteps_max)
      if (irecompute == 0) call m%compute_dt()
      irecompute = mod(irecompute + 1, 20)
      if (m%t + m%dt > p%tend) m%dt = p%tend - m%t
      call m%advance()
      if (mod(m%step, p%iprobe) == 0) call d%probe(m)
      if (mod(m%step, p%idiag) == 0) call d%write(m)
      if (p%iout > 0) then
        if (mod(m%step, p%iout) == 0) call write_fields(m, p%prefix)
      end if
    end do
    call cpu_time(t_wall1)
    if (mod(m%step, p%idiag) /= 0) call d%write(m)
    call d%close()
    write(*,'(a,i0,a,f12.4,a,f10.2,a)') 'hydft: finished ', m%step, ' steps, t =', m%t, '  (', t_wall1 - t_wall0, ' s)'
  end subroutine hydft_run

end module hydft

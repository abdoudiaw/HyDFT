! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Linear response from the same input file as the simulation:
!>   <prefix>_dispersion.dat : q, Re w, Im w, sqrt(q^2 c2) (undamped), S(q)
!>   <prefix>_dsf.dat        : w, S(q,w)/n0 for each q in q_dsf
program hydft_linear
  use hydft_kinds
  use hydft_utils
  use hydft_params
  use hydft_spectral
  use hydft_species
  use hydft_linear_response
  implicit none
  type(params_t) :: p
  type(spec_t) :: spc
  type(species_t), target :: sp
  type(linear_t) :: lin
  character(len=256) :: infile
  integer :: u, i, j, nqd
  real(rp) :: q, w, sq, dw
  complex(rp) :: wc

  if (command_argument_count() < 1) then
    write(*,'(a)') 'usage: hydft_linear input.in'
    stop
  end if
  call get_command_argument(1, infile)
  call p%read(infile)
  call p%print()
  call spc%init(p%nx, p%ny, p%nz, p%lx, p%ly, p%lz, p%dealias)
  call sp%build(p, spc)
  call lin%init(sp)

  u = open_file(trim(p%prefix)//'_dispersion.dat', 'write')
  write(u,'(a)') '# q    Re(w)/wp    Im(w)/wp    w0/wp=q*sqrt(n0 K/m)    S(q)'
  do i = 0, p%nq
    q = p%qmin + (p%qmax - p%qmin)*i/max(p%nq, 1)
    wc = lin%dispersion(q)
    sq = -1.0_rp
    if (allocated(sp%st%k)) sq = sp%st%s_of_k(q)
    write(u,'(5es18.8e3)') q, real(wc), aimag(wc), q*sqrt(max(lin%c2(q), 0.0_rp)), sq
  end do
  close(u)

  nqd = count(p%q_dsf > 0.0_rp)
  u = open_file(trim(p%prefix)//'_dsf.dat', 'write')
  write(u,'(a)', advance='no') '# w/wp'
  do j = 1, size(p%q_dsf)
    if (p%q_dsf(j) > 0.0_rp) write(u,'(a,f8.4)', advance='no') '   S(q,w)/n0 q=', p%q_dsf(j)
  end do
  write(u,*)
  dw = 2.0_rp*p%omega_max/p%nomega
  do i = 0, p%nomega
    w = -p%omega_max + dw*i
    write(u,'(es18.8e3)', advance='no') w
    do j = 1, size(p%q_dsf)
      if (p%q_dsf(j) > 0.0_rp) write(u,'(es18.8e3)', advance='no') lin%dsf(p%q_dsf(j), w)
    end do
    write(u,*)
  end do
  close(u)
  write(*,'(a)') 'hydft_linear: wrote '//trim(p%prefix)//'_dispersion.dat and '//trim(p%prefix)//'_dsf.dat'
end program hydft_linear

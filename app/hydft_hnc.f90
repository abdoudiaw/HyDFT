! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Standalone structure calculation: potential + parameters -> g(r), c(k), S(k).
!> Usage: hydft_hnc [input.in]   (namelist &hnc; defaults below)
program hydft_hnc
  use hydft_kinds
  use hydft_utils
  use hydft_pair_potential
  use hydft_oz_hnc
  use hydft_structure_source
  implicit none
  character(len=32)  :: potential = 'yukawa'
  real(rp) :: gamma = 10.0_rp, kappa = 1.0_rp, rs = 1.0_rp, lambda = -1.0_rp
  integer  :: nr = 8192, maxiter = 5000
  real(rp) :: rmax = 80.0_rp, tol = 1.0e-9_rp, mix = 0.2_rp
  character(len=256) :: prefix = 'hnc'
  logical :: verbose = .false.
  namelist /hnc/ potential, gamma, kappa, rs, lambda, nr, rmax, tol, mix, maxiter, prefix, verbose
  class(pair_potential_t), allocatable :: pot
  type(structure_t) :: st
  type(hnc_result_t) :: res
  character(len=256) :: infile
  integer :: u, i
  real(rp) :: n0

  if (command_argument_count() >= 1) then
    call get_command_argument(1, infile)
    u = open_file(infile, 'read')
    read(u, nml=hnc)
    close(u)
  end if
  n0 = 3.0_rp/(4.0_rp*pi)
  if (lambda < 0.0_rp) lambda = 2.0_rp*sqrt(pi*gamma/rs)   ! HM electron-pair de Broglie length / a
  call make_pair_potential(potential, gamma, kappa, lambda, pot)
  write(*,'(a,a,a,f8.3,a,f8.3)') 'hydft_hnc: potential ', trim(potential), '  Gamma =', gamma, '  kappa =', kappa
  call structure_from_hnc(pot, n0, nr, rmax, tol, mix, maxiter, st, with_derivative=.true., verbose=verbose, res_out=res)
  write(*,'(a,i0,a,es10.2)') '  converged in ', res%niter, ' iterations, residual ', res%residual
  write(*,'(a,f12.6,a,f12.6)') '  u_ex/NkT = ', st%u_ex, '   du/dGamma = ', st%du_dgamma
  write(*,'(a,f12.6)') '  S(k->0) = ', st%s(1)
  u = open_file(trim(prefix)//'_gr.dat', 'write')
  write(u,'(a)') '# r/a    g(r)    c(r)'
  do i = 1, res%nr
    write(u,'(3es20.10)') res%r(i), res%g(i), res%cr(i)
  end do
  close(u)
  call structure_write(st, trim(prefix)//'_sk.dat')
  write(*,'(a)') '  wrote '//trim(prefix)//'_gr.dat and '//trim(prefix)//'_sk.dat'
end program hydft_hnc

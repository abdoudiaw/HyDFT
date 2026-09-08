! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Where the static structure comes from: an in-code HNC calculation or a
!> file of S(k) or c(k) (e.g. from molecular dynamics). Either way the result
!> is a table c(k) on a radial k grid that the correlation functional and the
!> linear-response tool interpolate.
module hydft_structure_source
  use hydft_kinds
  use hydft_utils
  use hydft_pair_potential
  use hydft_oz_hnc
  implicit none
  private
  public :: structure_t, structure_from_hnc, structure_from_file, structure_write

  type :: structure_t
    character(len=32) :: source = 'none'
    real(rp) :: n0 = 3.0_rp/(4.0_rp*pi)
    real(rp), allocatable :: k(:)     !< radial wavenumbers (increasing, k > 0)
    real(rp), allocatable :: c(:)     !< total direct correlation function c(k)
    real(rp), allocatable :: s(:)     !< static structure factor S(k) = 1/(1 - n0 c(k))
    real(rp) :: u_ex = 0.0_rp         !< excess energy per particle / kT (HNC only; 0 otherwise)
    real(rp) :: du_dgamma = 0.0_rp    !< d u_ex / d Gamma (HNC only, by finite difference)
    logical :: has_energy = .false.
  contains
    procedure :: c_of_k => structure_c_of_k
    procedure :: s_of_k => structure_s_of_k
  end type structure_t

contains

  !> c(k) interpolated at k (linear); constant extrapolation beyond the table.
  pure real(rp) function structure_c_of_k(st, k) result(c)
    class(structure_t), intent(in) :: st
    real(rp), intent(in) :: k
    c = interp_linear(st%k, st%c, min(k, st%k(size(st%k))))
  end function structure_c_of_k

  pure real(rp) function structure_s_of_k(st, k) result(s)
    class(structure_t), intent(in) :: st
    real(rp), intent(in) :: k
    s = 1.0_rp/(1.0_rp - st%n0*st%c_of_k(k))
  end function structure_s_of_k

  !> Solve HNC for `pot`; optionally also the Gamma-derivative of the excess energy.
  subroutine structure_from_hnc(pot, n0, nr, rmax, tol, mix, maxiter, st, with_derivative, verbose, res_out)
    class(pair_potential_t), intent(inout) :: pot
    real(rp), intent(in) :: n0, rmax, tol, mix
    integer, intent(in) :: nr, maxiter
    type(structure_t), intent(out) :: st
    logical, intent(in), optional :: with_derivative, verbose
    type(hnc_result_t), intent(out), optional :: res_out
    type(hnc_result_t) :: res, rp1, rm1
    real(rp) :: g0, dg
    logical :: wd
    wd = .false.
    if (present(with_derivative)) wd = with_derivative
    call hnc_solve(pot, n0, nr, rmax, tol, mix, maxiter, res, verbose)
    if (.not. res%converged) call fatal('structure: HNC did not converge')
    st%source = 'hnc'
    st%n0 = n0
    st%k = res%k
    st%c = res%ck
    st%s = res%sk
    st%u_ex = res%u_ex
    st%has_energy = .true.
    if (wd) then
      g0 = pot%gamma
      dg = 0.02_rp*g0
      pot%gamma = g0 + dg
      call hnc_solve(pot, n0, nr, rmax, tol, mix, maxiter, rp1)
      pot%gamma = g0 - dg
      call hnc_solve(pot, n0, nr, rmax, tol, mix, maxiter, rm1)
      pot%gamma = g0
      if (.not. (rp1%converged .and. rm1%converged)) call fatal('structure: HNC derivative runs did not converge')
      st%du_dgamma = (rp1%u_ex - rm1%u_ex)/(2.0_rp*dg)
    end if
    if (present(res_out)) res_out = res
  end subroutine structure_from_hnc

  !> Read a two-column file (k, S(k)) or (k, c(k)); `kind` = 's' or 'c'. Lines
  !> starting with '#' are skipped. k must be in units of 1/a, increasing.
  subroutine structure_from_file(path, kind, n0, st)
    character(len=*), intent(in) :: path, kind
    real(rp), intent(in) :: n0
    type(structure_t), intent(out) :: st
    integer :: u, ios, n, i
    character(len=512) :: line
    real(rp) :: a, b
    real(rp), allocatable :: kk(:), vv(:)
    u = open_file(path, 'read')
    n = 0
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) a, b
      if (ios /= 0) cycle
      n = n + 1
    end do
    if (n < 2) call fatal('structure file has fewer than two data rows: '//trim(path))
    allocate(kk(n), vv(n))
    rewind(u)
    i = 0
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) a, b
      if (ios /= 0) cycle
      i = i + 1
      kk(i) = a; vv(i) = b
    end do
    close(u)
    st%source = 'file'
    st%n0 = n0
    st%k = kk
    select case (lower(kind(1:1)))
    case ('s')
      st%s = vv
      st%c = (1.0_rp - 1.0_rp/vv)/n0
    case ('c')
      st%c = vv
      st%s = 1.0_rp/(1.0_rp - n0*vv)
    case default
      call fatal('structure file kind must be s or c')
    end select
    do i = 2, n
      if (st%k(i) <= st%k(i-1)) call fatal('structure file: k must be strictly increasing')
    end do
    if (st%k(1) <= 0.0_rp) call fatal('structure file: k must be positive')
  end subroutine structure_from_file

  !> Write k, c(k), S(k) to a file.
  subroutine structure_write(st, path)
    type(structure_t), intent(in) :: st
    character(len=*), intent(in) :: path
    integer :: u, i
    u = open_file(path, 'write')
    write(u,'(a)') '# k*a    c(k)    S(k)'
    do i = 1, size(st%k)
      write(u,'(3es20.10)') st%k(i), st%c(i), st%s(i)
    end do
    close(u)
  end subroutine structure_write

end module hydft_structure_source

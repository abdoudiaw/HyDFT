! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Ornstein-Zernike equation with the hypernetted-chain closure on a uniform
!> radial grid. Fourier-Bessel transforms are done with FFTW sine transforms
!> (DST-I). Picard iteration with Ng acceleration on the short-range indirect
!> correlation function gamma_s = h - c_s, c_s = c + beta v_L. A bridge
!> function hook (bridge = 0, plain HNC) is kept for later strong-coupling work.
!>
!> Reduced units: r in a (Wigner-Seitz radius), n0 = 3/(4 pi), energies in kT.
module hydft_oz_hnc
  use, intrinsic :: iso_c_binding
  use hydft_kinds
  use hydft_utils, only: trapz
  use hydft_pair_potential
  implicit none
  private
  public :: hnc_result_t, hnc_solve, fourier_bessel_t

  include 'fftw3.f03'

  !> Radial Fourier-Bessel transform pair on r_j = j dr, k_m = m dk, dk = pi/((n+1) dr).
  type :: fourier_bessel_t
    integer :: n = 0
    real(rp) :: dr = 0.0_rp, dk = 0.0_rp
    real(rp), allocatable :: r(:), k(:)
    type(c_ptr) :: plan = c_null_ptr, buf_in = c_null_ptr, buf_out = c_null_ptr
    real(c_double), pointer :: win(:) => null(), wout(:) => null()
  contains
    procedure :: init => fb_init
    procedure :: r2k => fb_r2k
    procedure :: k2r => fb_k2r
    procedure :: destroy => fb_destroy
  end type fourier_bessel_t

  type :: hnc_result_t
    integer :: nr = 0
    real(rp) :: n0 = 3.0_rp/(4.0_rp*pi)
    real(rp) :: dr = 0.0_rp, dk = 0.0_rp
    real(rp), allocatable :: r(:), k(:)
    real(rp), allocatable :: g(:)       !< pair distribution g(r)
    real(rp), allocatable :: cr(:)      !< direct correlation c(r)
    real(rp), allocatable :: ck(:)      !< c(k) (total, includes -beta v(k))
    real(rp), allocatable :: hk(:)      !< h(k)
    real(rp), allocatable :: sk(:)      !< S(k) = 1 + n0 h(k)
    real(rp) :: u_ex = 0.0_rp           !< excess (correlation) energy per particle / kT, (n/2) int (g-1) beta v
    logical :: converged = .false.
    integer :: niter = 0
    real(rp) :: residual = 0.0_rp
  end type hnc_result_t

contains

  subroutine fb_init(fb, n, rmax)
    class(fourier_bessel_t), intent(inout) :: fb
    integer, intent(in) :: n
    real(rp), intent(in) :: rmax
    integer :: j
    call fb%destroy()
    fb%n = n
    fb%dr = rmax/(n + 1)
    fb%dk = pi/((n + 1)*fb%dr)
    allocate(fb%r(n), fb%k(n))
    do j = 1, n
      fb%r(j) = j*fb%dr
      fb%k(j) = j*fb%dk
    end do
    fb%buf_in = fftw_alloc_real(int(n, c_size_t))
    fb%buf_out = fftw_alloc_real(int(n, c_size_t))
    call c_f_pointer(fb%buf_in, fb%win, [n])
    call c_f_pointer(fb%buf_out, fb%wout, [n])
    fb%plan = fftw_plan_r2r_1d(n, fb%win, fb%wout, FFTW_RODFT00, FFTW_MEASURE)
  end subroutine fb_init

  !> f(k) = 4 pi / k int_0^inf r f(r) sin(k r) dr
  subroutine fb_r2k(fb, fr, fk)
    class(fourier_bessel_t), intent(inout) :: fb
    real(rp), intent(in) :: fr(:)
    real(rp), intent(out) :: fk(:)
    fb%win = fb%r*fr
    call fftw_execute_r2r(fb%plan, fb%win, fb%wout)
    fk = twopi*fb%dr*fb%wout/fb%k
  end subroutine fb_r2k

  !> f(r) = 1/(2 pi^2 r) int_0^inf k f(k) sin(k r) dk
  subroutine fb_k2r(fb, fk, fr)
    class(fourier_bessel_t), intent(inout) :: fb
    real(rp), intent(in) :: fk(:)
    real(rp), intent(out) :: fr(:)
    fb%win = fb%k*fk
    call fftw_execute_r2r(fb%plan, fb%win, fb%wout)
    fr = fb%dk*fb%wout/(4.0_rp*pi*pi*fb%r)
  end subroutine fb_k2r

  subroutine fb_destroy(fb)
    class(fourier_bessel_t), intent(inout) :: fb
    if (c_associated(fb%plan)) call fftw_destroy_plan(fb%plan)
    if (c_associated(fb%buf_in)) call fftw_free(fb%buf_in)
    if (c_associated(fb%buf_out)) call fftw_free(fb%buf_out)
    fb%plan = c_null_ptr; fb%buf_in = c_null_ptr; fb%buf_out = c_null_ptr
    nullify(fb%win, fb%wout)
    if (allocated(fb%r)) deallocate(fb%r, fb%k)
  end subroutine fb_destroy

  !> Solve OZ + HNC for potential `pot` at number density n0 (reduced units).
  !> A direct solve is attempted first; if it is unstable, the coupling is
  !> ramped up from a weakly coupled converged solution (continuation).
  subroutine hnc_solve(pot, n0, nr, rmax, tol, mix, maxiter, res, verbose)
    class(pair_potential_t), intent(in) :: pot
    real(rp), intent(in) :: n0, rmax, tol, mix
    integer, intent(in) :: nr, maxiter
    type(hnc_result_t), intent(out) :: res
    logical, intent(in), optional :: verbose
    class(pair_potential_t), allocatable :: p
    real(rp), allocatable :: gs(:), gs_init(:)
    real(rp) :: gtarget, g
    logical :: unstable, verb
    integer :: istep
    verb = .false.
    if (present(verbose)) verb = verbose
    allocate(gs(nr), gs_init(nr))
    call hnc_solve_direct(pot, n0, nr, rmax, tol, mix, maxiter, res, unstable, gs, verbose)
    if (.not. unstable) return
    allocate(p, source=pot)
    gtarget = pot%gamma
    g = min(gtarget, 2.0_rp)
    gs_init = 0.0_rp
    do istep = 1, 200
      p%gamma = g
      if (verb) write(*,'(a,f10.4)') '  hnc: continuation at Gamma = ', g
      call hnc_solve_direct(p, n0, nr, rmax, tol, mix, maxiter, res, unstable, gs, verbose, gamma_s_init=gs_init)
      if (unstable) call fatal('hnc: continuation in Gamma failed')
      gs_init = gs
      if (g >= gtarget) exit
      g = min(gtarget, g*1.3_rp)
    end do
  end subroutine hnc_solve

  subroutine hnc_solve_direct(pot, n0, nr, rmax, tol, mix, maxiter, res, unstable_out, gs_out, verbose, gamma_s_init)
    class(pair_potential_t), intent(in) :: pot
    real(rp), intent(in) :: n0, rmax, tol, mix
    integer, intent(in) :: nr, maxiter
    type(hnc_result_t), intent(out) :: res
    logical, intent(out) :: unstable_out
    real(rp), intent(out) :: gs_out(:)
    logical, intent(in), optional :: verbose
    real(rp), intent(in), optional :: gamma_s_init(:)   !< starting gamma_s(r) for continuation
    type(fourier_bessel_t) :: fb
    real(rp), allocatable :: bvs(:), bvl_r(:), bvl_k(:), gs(:), gnew(:), cs(:), csk(:), ck(:), gk(:)
    real(rp), allocatable :: hist_in(:,:), hist_out(:,:), d0(:), d01(:), d02(:)
    real(rp) :: resid, resid_prev, a11, a12, a22, b1, b2, det, c1, c2, ex
    integer :: it, j, nhist, irestart, it_total, ng_cooldown
    logical :: accepted
    real(rp), allocatable :: gtry(:), gk_mf(:)
    real(rp) :: mixc
    logical :: unstable
    integer, parameter :: npicard = 5
    logical :: verb
    real(rp), parameter :: expcap = 40.0_rp

    verb = .false.
    if (present(verbose)) verb = verbose
    call fb%init(nr, rmax)
    allocate(bvs(nr), bvl_r(nr), bvl_k(nr), gs(nr), gnew(nr), cs(nr), csk(nr), ck(nr), gk(nr))
    allocate(hist_in(nr,3), hist_out(nr,3), d0(nr), d01(nr), d02(nr), gtry(nr), gk_mf(nr))
    do j = 1, nr
      bvs(j)   = pot%beta_v_short(fb%r(j))
      bvl_r(j) = pot%beta_v_long_r(fb%r(j))
      bvl_k(j) = pot%beta_v_long_k(fb%k(j))
    end do

    ! initial guess: mean-field (Debye-Hueckel) c(k) = -beta v(k), as in the notebook HNC solver
    do j = 1, nr
      ck(j) = -pot%beta_v_k(fb%k(j))
      gk(j) = n0*ck(j)*ck(j)/(1.0_rp - n0*ck(j)) - bvl_k(j)
    end do
    call fb%k2r(gk, gs)
    if (present(gamma_s_init)) then
      if (size(gamma_s_init) == nr) gs = gamma_s_init
    end if
    mixc = mix
    res%converged = .false.
    do j = 1, nr
      ck(j) = -pot%beta_v_k(fb%k(j))
      gk_mf(j) = n0*ck(j)*ck(j)/(1.0_rp - n0*ck(j)) - bvl_k(j)
    end do
    restart: do irestart = 1, 4
      hist_in = 0.0_rp; hist_out = 0.0_rp
      nhist = 0
      ng_cooldown = 0
      resid_prev = 1.0e100_rp
      unstable = .false.
      it_total = 0
      do it = 1, maxiter
        it_total = it
        ! closure
        do j = 1, nr
          ex = -bvs(j) + gs(j)
          cs(j) = exp(min(ex, expcap)) - 1.0_rp - gs(j)
        end do
        call fb%r2k(cs, csk)
        ck = csk - bvl_k
        if (any(n0*ck >= 1.0_rp)) then
          unstable = .true.
          exit
        end if
        do j = 1, nr
          gk(j) = n0*ck(j)*ck(j)/(1.0_rp - n0*ck(j)) - bvl_k(j)
        end do
        call fb%k2r(gk, gnew)
        resid = maxval(abs(gnew - gs))
        if (verb) write(*,'(a,i6,a,es12.4)') '  hnc iter ', it, '  residual ', resid
        if (resid < tol) then
          gs = gnew
          res%converged = .true.
          exit
        end if
        if (resid > 1.0e6_rp) then
          unstable = .true.
          exit
        end if
        ! Ng acceleration with the last three (input, output) pairs; fall back to Picard mixing
        hist_in(:,3) = hist_in(:,2); hist_in(:,2) = hist_in(:,1); hist_in(:,1) = gs
        hist_out(:,3) = hist_out(:,2); hist_out(:,2) = hist_out(:,1); hist_out(:,1) = gnew
        nhist = min(nhist + 1, 3)
        if (resid > 10.0_rp*resid_prev) ng_cooldown = 10
        if (nhist == 3 .and. it > npicard .and. ng_cooldown == 0) then
          d0 = hist_out(:,1) - hist_in(:,1)
          d01 = d0 - (hist_out(:,2) - hist_in(:,2))
          d02 = d0 - (hist_out(:,3) - hist_in(:,3))
          a11 = dot_product(d01, d01); a12 = dot_product(d01, d02); a22 = dot_product(d02, d02)
          b1 = dot_product(d0, d01); b2 = dot_product(d0, d02)
          det = a11*a22 - a12*a12
          accepted = .false.
          if (abs(det) > 1.0e-14_rp*a11*a22 .and. a11 > 0.0_rp .and. a22 > 0.0_rp) then
            c1 = (b1*a22 - b2*a12)/det
            c2 = (a11*b2 - a12*b1)/det
            if (abs(c1) < 20.0_rp .and. abs(c2) < 20.0_rp) then
              gtry = (1.0_rp - c1 - c2)*hist_out(:,1) + c1*hist_out(:,2) + c2*hist_out(:,3)
              if (maxval(gtry - bvs) < expcap) then
                gs = gtry
                accepted = .true.
              end if
            end if
          end if
          if (.not. accepted) gs = (1.0_rp - mixc)*gs + mixc*gnew
        else
          gs = (1.0_rp - mixc)*gs + mixc*gnew
          if (ng_cooldown > 0) ng_cooldown = ng_cooldown - 1
        end if
        resid_prev = resid
      end do
      if (res%converged .or. .not. unstable) exit restart
      ! diverged: restart from scratch with a smaller mixing parameter
      mixc = 0.5_rp*mixc
      call fb%k2r(gk_mf, gs)
      if (present(gamma_s_init)) then
        if (size(gamma_s_init) == nr) gs = gamma_s_init
      end if
      if (verb) write(*,'(a,es10.2)') '  hnc: unstable, restarting with mix = ', mixc
    end do restart
    it = it_total
    res%niter = min(it, maxiter)
    res%residual = resid
    unstable_out = unstable
    gs_out = gs
    if (unstable) then
      call fb%destroy()
      return
    end if
    if (.not. res%converged) call warn('hnc: not converged; residual printed in result')

    ! final quantities
    res%nr = nr; res%n0 = n0; res%dr = fb%dr; res%dk = fb%dk
    allocate(res%r(nr), res%k(nr), res%g(nr), res%cr(nr), res%ck(nr), res%hk(nr), res%sk(nr))
    res%r = fb%r; res%k = fb%k
    do j = 1, nr
      ex = -bvs(j) + gs(j)
      res%g(j) = exp(min(ex, expcap))
      cs(j) = res%g(j) - 1.0_rp - gs(j)
    end do
    call fb%r2k(cs, csk)
    res%ck = csk - bvl_k
    res%hk = res%ck/(1.0_rp - n0*res%ck)
    res%sk = 1.0_rp + n0*res%hk
    res%cr = cs - bvl_r
    ! excess energy per particle: (n/2) int d^3r (g-1) beta v(r)
    do j = 1, nr
      gnew(j) = fb%r(j)**2*(res%g(j) - 1.0_rp)*pot%beta_v(fb%r(j))
    end do
    res%u_ex = 2.0_rp*pi*n0*(fb%dr*sum(gnew))   ! rectangle rule; integrand vanishes at both ends
    call fb%destroy()
  end subroutine hnc_solve_direct

end module hydft_oz_hnc

! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> HNC: transform normalisation, OCP energy vs the Slattery-Doolen-DeWitt fit,
!> Yukawa S(q) vs mean field at long wavelength, electron QSP case runs.
program test_hnc
  use hydft_kinds
  use hydft_pair_potential
  use hydft_oz_hnc
  use hydft_structure_source
  implicit none
  integer :: nfail, j, nr
  type(fourier_bessel_t) :: fb
  real(rp), allocatable :: fr(:), fk(:), fr2(:)
  class(pair_potential_t), allocatable :: pot
  type(hnc_result_t) :: res
  type(structure_t) :: st
  real(rp) :: err, n0, kap, gam, ufit, smf, q
  nfail = 0
  n0 = 3.0_rp/(4.0_rp*pi)

  ! Fourier-Bessel transform of e^{-r} -> 8 pi/(k^2+1)^2, and back
  nr = 4096
  call fb%init(nr, 60.0_rp)
  allocate(fr(nr), fk(nr), fr2(nr))
  fr = exp(-fb%r)
  call fb%r2k(fr, fk)
  err = maxval(abs(fk*(fb%k**2 + 1.0_rp)**2/(8.0_rp*pi) - 1.0_rp), mask=(fb%k < 20.0_rp))
  call check(err < 1e-4_rp, 'FB transform of e^-r', err)
  call fb%k2r(fk, fr2)
  err = maxval(abs(fr2 - fr))
  call check(err < 1e-12_rp, 'FB round trip', err)
  call fb%destroy()

  ! OCP at Gamma = 10 and 50: HNC energy within 3% of the SDD fit
  do j = 1, 2
    gam = merge(10.0_rp, 50.0_rp, j == 1)
    call make_pair_potential('coulomb', gam, 0.0_rp, 1.0_rp, pot)
    call hnc_solve(pot, n0, 8192, 80.0_rp, 1e-9_rp, 0.2_rp, 5000, res)
    ufit = -0.897744_rp*gam + 0.95043_rp*gam**0.25_rp + 0.18956_rp*gam**(-0.25_rp) - 0.81487_rp
    err = abs(res%u_ex/ufit - 1.0_rp)
    call check(res%converged, 'OCP HNC converged', real(res%niter, rp))
    call check(err < 0.03_rp, 'OCP energy vs SDD fit', err)
    call check(abs(res%g(1)) < 1e-6_rp .and. abs(res%g(res%nr) - 1.0_rp) < 1e-6_rp, 'g(0)=0, g(inf)=1', res%g(res%nr))
    ! perfect screening: S(k) -> k^2/(3 Gamma) as k -> 0
    err = abs(res%sk(1)/(res%k(1)**2/(3.0_rp*gam)) - 1.0_rp)
    call check(err < 0.05_rp, 'OCP S(k->0) ~ k^2/3Gamma', err)
  end do

  ! Yukawa Gamma = 5, kappa = 1 (PRE Fig. 3b): S(q) close to mean field (q^2+kappa^2)/(q^2+kappa^2+3Gamma) at q <= 0.5
  gam = 5.0_rp; kap = 1.0_rp
  call make_pair_potential('yukawa', gam, kap, 1.0_rp, pot)
  call structure_from_hnc(pot, n0, 8192, 80.0_rp, 1e-9_rp, 0.2_rp, 5000, st, with_derivative=.true.)
  q = 0.5_rp
  smf = (q*q + kap*kap)/(q*q + kap*kap + 3.0_rp*gam)
  err = abs(st%s_of_k(q)/smf - 1.0_rp)
  call check(err < 0.25_rp, 'Yukawa S(0.5) vs mean field', err)
  call check(st%u_ex < 0.0_rp .and. st%du_dgamma < 0.0_rp, 'Yukawa u_ex < 0, du/dGamma < 0', st%u_ex)
  write(*,'(a,f10.5,a,f10.5)') '  info  Yukawa(5,1) u_ex = ', st%u_ex, '  du/dG = ', st%du_dgamma

  ! Strong coupling Yukawa (PRE Fig. 4): Gamma = 150, kappa = 0.1 converges
  call make_pair_potential('yukawa', 150.0_rp, 0.1_rp, 1.0_rp, pot)
  call hnc_solve(pot, n0, 8192, 80.0_rp, 1e-8_rp, 0.1_rp, 20000, res)
  call check(res%converged, 'Yukawa(150,0.1) converged', real(res%niter, rp))
  call check(maxval(res%g) > 1.5_rp, 'Yukawa(150,0.1) first peak > 1.5', maxval(res%g))

  ! Electron QSP, rs = 1.86, Gamma = 1 (Sci. Rep. Fig. 2)
  call make_pair_potential('hansen_mcdonald', 1.0_rp, 0.0_rp, 2.0_rp*sqrt(pi*1.0_rp/1.86_rp), pot)
  call hnc_solve(pot, n0, 8192, 80.0_rp, 1e-9_rp, 0.2_rp, 5000, res)
  call check(res%converged, 'HM QSP converged', real(res%niter, rp))
  call check(res%g(1) > 0.0_rp .and. res%g(1) < 1.0_rp, 'HM QSP g(0) finite', res%g(1))
  err = abs(res%sk(1)/(res%k(1)**2/3.0_rp) - 1.0_rp)
  call check(err < 0.05_rp, 'HM QSP perfect screening', err)

  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_hnc: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_hnc: all passed'
contains
  subroutine check(ok, name, val)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: name
    real(rp), intent(in) :: val
    if (ok) then
      write(*,'(a,a,es12.3)') '  pass  ', name, val
    else
      write(*,'(a,a,es12.3)') '  FAIL  ', name, val
      nfail = nfail + 1
    end if
  end subroutine check
end program test_hnc

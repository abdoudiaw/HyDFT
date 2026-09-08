! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Linear response: classical Yukawa reproduces PRE Eq. 54, tau -> 0 gives the
!> Navier-Stokes root, the electron TF + Hartree case gives Bohm-Gross, the DSF
!> obeys detailed balance.
program test_dispersion
  use hydft_kinds
  use hydft_params
  use hydft_spectral
  use hydft_species
  use hydft_linear_response
  implicit none
  integer :: nfail, i, j
  type(params_t) :: p
  type(spec_t) :: spc
  type(species_t), target :: sp
  type(linear_t) :: lin
  real(rp) :: q, sq, err, nu, c2, wbg, x
  complex(rp) :: w, wq, wns, roots(3)
  nfail = 0

  ! ---- classical Yukawa ions, Maxwell closure with constant transport
  p%kind = 'ion'; p%gamma = 2.0_rp; p%kappa = 0.1_rp
  p%nx = 16; p%lx = 20.0_rp; p%hnc_nr = 4096; p%hnc_rmax = 60.0_rp
  p%transport = 'constant'; p%eta = 0.3_rp; p%xi = 0.05_rp
  p%closure = 'maxwell'; p%tau_model = 'constant'; p%tau = 2.0_rp
  call spc%init(p%nx, 1, 1, p%lx, 1.0_rp, 1.0_rp)
  call sp%build(p, spc, verbose=.false.)
  call lin%init(sp)
  err = abs(lin%eta_bar - (4.0_rp*0.3_rp/3.0_rp + 0.05_rp)); call check(err < 1e-14_rp, 'eta_l = 4 eta/3 + xi', err)
  do i = 1, 5
    q = 0.4_rp*i
    sq = sp%st%s_of_k(q)
    ! kernel assembly: n0 K = 1/(beta S(q))
    err = abs(lin%c2(q)*sp%un%beta*sq - 1.0_rp); call check(err < 1e-10_rp, 'n0 K(q) = 1/(beta S(q))', err)
    w = lin%dispersion(q, roots)
    ! PRE Eq. 54 residual
    err = abs(w*w - q*q/(3.0_rp*p%gamma*sq) + i_unit*lin%eta_bar*q*q*w/(1.0_rp - i_unit*w*p%tau))
    call check(err < 1e-12_rp, 'Eq. 54 satisfied', err)
    call check(real(w) > 0.0_rp .and. aimag(w) < 0.0_rp, 'propagating damped root', aimag(w))
    ! all three roots satisfy the cubic
    err = 0.0_rp
    do j = 1, 3
      err = max(err, abs(lin%denominator(q, roots(j))*(1.0_rp - i_unit*roots(j)*p%tau)))
    end do
    call check(err < 1e-10_rp, 'all roots of the cubic', err)
  end do
  ! tau -> 0 limit: Newtonian root -i nu/2 + sqrt(c^2 k^2 - nu^2/4)
  p%closure = 'newtonian'
  call sp%build(p, spc, verbose=.false.)
  call lin%init(sp)
  q = 1.0_rp
  nu = q*q*lin%eta_bar
  c2 = lin%c2(q)
  wns = -0.5_rp*i_unit*nu + sqrt(cmplx(q*q*c2 - 0.25_rp*nu*nu, 0.0_rp, rp))
  w = lin%dispersion(q)
  err = abs(w - wns); call check(err < 1e-12_rp, 'Newtonian root', err)
  ! static response chi(q,0) = -1/K
  err = abs(real(lin%chi(q, 0.0_rp)) + 1.0_rp/sp%fun%kernel(q)); call check(err < 1e-10_rp, 'chi(q,0) = -1/K', err)
  ! classical DSF: positive, symmetric
  err = abs(lin%dsf(q, 0.7_rp) - lin%dsf(q, -0.7_rp)); call check(err < 1e-14_rp, 'classical S(q,w) even', err)
  call check(lin%dsf(q, 0.7_rp) > 0.0_rp, 'classical S(q,w) > 0', lin%dsf(q, 0.7_rp))

  ! ---- electrons: TF + Hartree, no viscosity -> Bohm-Gross-like w^2 = wp^2 + q^2 n0 K_TF
  p%kind = 'electron'; p%rs = 1.86_rp; p%gamma = 1.0_rp
  p%terms = 'tf,hartree'; p%structure = 'none'; p%eta = 0.0_rp; p%xi = 0.0_rp
  call sp%build(p, spc, verbose=.false.)
  call lin%init(sp)
  do i = 1, 4
    q = 0.5_rp*i
    w = lin%dispersion(q)
    err = abs(w*w - q*q*lin%c2(q)); call check(err < 1e-12_rp, 'electron undamped root', err)
    err = abs(aimag(w)); call check(err < 1e-12_rp, 'no damping without viscosity', err)
    err = abs(real(w)**2 - (1.0_rp + q*q*sp%un%n0*sp%fun%terms(1)%t%kernel(q)))
    call check(err < 1e-12_rp, 'w^2 = wp^2 + q^2 n0 K_TF', err)
  end do
  ! classical limit of TF: K_TF -> 1/(beta n0) (Bohm-Gross w^2 = 1 + q^2/(3 Gamma))
  p%gamma = 0.01_rp
  call sp%build(p, spc, verbose=.false.)
  call lin%init(sp)
  q = 1.0_rp
  w = lin%dispersion(q)
  wbg = sqrt(1.0_rp + q*q/(3.0_rp*p%gamma))
  err = abs(real(w)/wbg - 1.0_rp); call check(err < 5e-3_rp, 'Bohm-Gross in the classical limit', err)
  write(*,'(a,f10.4)') '  info  theta =', sp%un%theta
  ! quantum DSF: detailed balance S(-w) = exp(-beta hbar w) S(w), with viscosity so Im chi /= 0
  p%gamma = 1.0_rp; p%eta = 0.2_rp
  call sp%build(p, spc, verbose=.false.)
  call lin%init(sp)
  q = 1.0_rp; x = 0.8_rp
  err = abs(lin%dsf(q, -x) - exp(-sp%un%beta*sp%un%hbar*x)*lin%dsf(q, x))/lin%dsf(q, x)
  call check(err < 1e-12_rp, 'detailed balance', err)
  call check(lin%dsf(q, x) > 0.0_rp, 'quantum S(q,w) > 0', lin%dsf(q, x))

  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_dispersion: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_dispersion: all passed'
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
end program test_dispersion

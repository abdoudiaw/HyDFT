! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Stanton-Murillo transport fits: continuity of the K_nm fits at g = 1, self-diffusion
!> against the MD values of PRE 93, 043203 Table I, the binary-mixture formulas reducing
!> to the one-component ones, and the electron viscosity interpolation.
program test_transport
  use hydft_kinds
  use hydft_stanton_murillo
  use hydft_units
  implicit none
  integer :: nfail, i, j, idx
  real(rp) :: err, kl, kr, d, dref, eta, kth, kt, d12, eta_l, eta_cv, eta_sm
  type(units_t) :: un
  ! Table I (kappa, Gamma, D*) for Gamma <= 5
  real(rp), parameter :: kap(4) = [0.5_rp, 1.0_rp, 1.5_rp, 2.0_rp]
  real(rp), parameter :: gam(7) = [0.1_rp, 0.3_rp, 0.5_rp, 0.7_rp, 1.0_rp, 2.0_rp, 5.0_rp]
  real(rp), parameter :: dmd(7, 4) = reshape([ &
    143.0_rp, 16.4_rp, 6.4_rp, 3.76_rp, 2.14_rp, 0.82_rp, 0.28_rp, &
    173.0_rp, 20.9_rp, 8.41_rp, 4.72_rp, 2.54_rp, 1.02_rp, 0.31_rp, &
    209.0_rp, 24.5_rp, 10.0_rp, 5.92_rp, 3.15_rp, 1.19_rp, 0.356_rp, &
    241.0_rp, 28.4_rp, 12.3_rp, 7.32_rp, 3.72_rp, 1.52_rp, 0.504_rp], [7, 4])
  nfail = 0

  ! C^1-continuous fits: both branches agree at g = 1 to the fit accuracy (~1e-2 relative)
  do idx = IDX_K11, IDX_K22
    kl = sm_knm(idx, 1.0_rp - 1e-12_rp); kr = sm_knm(idx, 1.0_rp + 1e-12_rp)
    err = abs(kl/kr - 1.0_rp); call check(err < 5e-3_rp, 'K_nm continuous at g = 1', err)
  end do
  ! K_22(1) from the strongly coupled branch: b0/(1 + b3 + b4) = 0.32237
  err = abs(sm_k22(1.0_rp) - 0.32237_rp); call check(err < 1e-4_rp, 'K_22(1)', err)
  ! weak-coupling limit K_11 ~ -(1/4) ln(a1 g)
  err = abs(sm_k11(1e-6_rp)/(-0.25_rp*log(1.4660e-6_rp)) - 1.0_rp); call check(err < 1e-4_rp, 'K_11 weak coupling', err)

  ! self-diffusion vs MD (effective-Boltzmann with the effective screening length, Fig. 10):
  ! within 25 % for Gamma <= 2 (worst: kappa = 2, Gamma = 1, +20 %), 35 % at Gamma = 5 where EB falls below MD
  do j = 1, 4
    do i = 1, 7
      d = sm_diffusion(gam(i), kap(j)); dref = dmd(i, j)
      err = abs(d/dref - 1.0_rp)
      if (gam(i) <= 2.0_rp) then
        call check(err < 0.25_rp, 'D* vs MD (Table I), Gamma <= 2', err)
      else
        call check(err < 0.35_rp, 'D* vs MD (Table I), Gamma = 5', err)
      end if
    end do
  end do
  ! viscosity and conductivity share K_22: K*/eta* = 15/4 (Prandtl-like ratio)
  err = abs(sm_conductivity(1.5_rp, 0.7_rp)/sm_viscosity(1.5_rp, 0.7_rp) - 3.75_rp); call check(err < 1e-12_rp, 'K*/eta* = 15/4', err)

  ! binary mixture with identical species reduces to the one-component results for any x1
  do i = 1, 3
    call sm_mixture(2.0_rp, 0.5_rp, 0.25_rp*i, 1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp, d12, eta, kth, kt)
    err = abs(d12/sm_diffusion(2.0_rp, 0.5_rp) - 1.0_rp);    call check(err < 1e-10_rp, 'mixture D_12 -> D', err)
    err = abs(eta/sm_viscosity(2.0_rp, 0.5_rp) - 1.0_rp);    call check(err < 1e-10_rp, 'mixture eta_tot -> eta', err)
    err = abs(kth/sm_conductivity(2.0_rp, 0.5_rp) - 1.0_rp); call check(err < 1e-10_rp, 'mixture K_tot -> K', err)
    call check(abs(kt) < 1e-12_rp, 'mixture k_T = 0 for identical species', abs(kt))
  end do
  ! H+/He2+ 50/50 BIM of Fig. 12 (Gamma_12 = 0.8, kappa = 0); values cross-checked with python/stanton_murillo.py
  call sm_mixture(0.8_rp, 0.0_rp, 0.5_rp, 1.0_rp, 2.0_rp, 1.0_rp, 4.0_rp, d12, eta, kth, kt)
  write(*,'(a,4es14.6)') '  info  H-He mixture d12*, eta*, K*, k_T =', d12, eta, kth, kt
  call sm_mixture(0.8_rp, 0.0_rp, 0.5_rp, 2.0_rp, 1.0_rp, 4.0_rp, 1.0_rp, d, eta_l, eta_cv, eta_sm)
  err = abs(d12/d - 1.0_rp) + abs(eta/eta_l - 1.0_rp) + abs(kth/eta_cv - 1.0_rp) + abs(kt + eta_sm)
  call check(err < 1e-12_rp, 'mixture symmetric under 1 <-> 2', err)

  ! electron viscosity: r_s = 1.86, Gamma = 1 (theta = 1.010, lambda_TF/a = 0.6472) -> eta_l = 1.312
  ! (independent scipy evaluation in examples/electron_dsf/electron_dsf.ipynb)
  call un%init('electron', 1.0_rp, 0.0_rp, 1.86_rp, 1.0_rp)
  err = abs(un%lambda_tf - 0.6472_rp); call check(err < 5e-4_rp, 'lambda_TF/a (r_s = 1.86, Gamma = 1)', err)
  call sm_electron_viscosity(un%rs, un%gamma, un%theta, un%lambda_tf, eta_l, eta_cv, eta_sm)
  err = abs(eta_l - 1.312_rp); call check(err < 3e-3_rp, 'electron eta_l (r_s = 1.86, Gamma = 1)', err)
  err = abs(eta_cv - 0.00467_rp); call check(err < 2e-5_rp, 'electron eta_CV', err)
  ! classical limit of lambda_TF: Debye length 1/sqrt(3 Gamma) at small Gamma (theta >> 1)
  call un%init('electron', 0.01_rp, 0.0_rp, 1.0_rp, 1.0_rp)
  err = abs(un%lambda_tf*sqrt(3.0_rp*0.01_rp) - 1.0_rp); call check(err < 5e-3_rp, 'lambda_TF -> Debye', err)

  if (nfail == 0) then
    write(*,'(a)') 'test_transport: PASS'
  else
    write(*,'(a,i0,a)') 'test_transport: FAIL (', nfail, ' checks)'
    error stop 1
  end if
contains
  subroutine check(ok, what, val)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: what
    real(rp), intent(in) :: val
    if (ok) then
      write(*,'(a,a,es10.2)') '  ok    ', what, val
    else
      write(*,'(a,a,es10.2)') '  FAIL  ', what, val
      nfail = nfail + 1
    end if
  end subroutine check
end program test_transport

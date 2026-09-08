! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Each free-energy term's linear kernel K_i(k) must equal the finite
!> difference of its nonlinear functional derivative mu_i[n] about n0.
!> This ties the solver and the linear-response code paths together and
!> catches double counting of the mean field between Hartree and RY terms.
program test_kernels
  use hydft_kinds
  use hydft_params
  use hydft_spectral
  use hydft_species
  implicit none
  integer :: nfail
  nfail = 0
  call run_case('ion', 'ideal,hartree,ry')
  call run_case('electron', 'tf,gradient,hartree,ry')
  call run_case('electron', 'tf,gradient,ry')
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_kernels: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_kernels: all passed'
contains
  subroutine run_case(kind, terms)
    character(len=*), intent(in) :: kind, terms
    type(params_t) :: p
    type(spec_t) :: spc
    type(species_t) :: sp
    real(rp), allocatable :: n(:,:,:), mu(:,:,:), mu0(:,:,:), n0a(:,:,:)
    complex(rp), allocatable :: nhat(:,:,:), n0hat(:,:,:)
    real(rp), parameter :: eps = 1.0e-6_rp
    real(rp) :: k, kfd, kan, err
    integer :: m, i, it
    p%kind = kind; p%terms = terms
    p%gamma = 2.0_rp; p%kappa = 0.1_rp; p%rs = 1.86_rp
    if (kind == 'electron') p%gamma = 1.0_rp
    p%nx = 64; p%ny = 1; p%nz = 1; p%lx = 20.0_rp
    p%hnc_nr = 4096; p%hnc_rmax = 60.0_rp
    p%transport = 'constant'; p%closure = 'newtonian'
    call spc%init(p%nx, p%ny, p%nz, p%lx, p%ly, p%lz)
    call sp%build(p, spc, verbose=.false.)
    allocate(n(p%nx,1,1), mu(p%nx,1,1), mu0(p%nx,1,1), n0a(p%nx,1,1), nhat(spc%g%nkx,1,1), n0hat(spc%g%nkx,1,1))
    n0a = sp%un%n0
    call spc%fwd(n0a, n0hat)
    write(*,'(a,a,a,a)') '  case ', trim(kind), ': ', trim(terms)
    do m = 1, 6
      k = twopi*m/p%lx
      do i = 1, p%nx
        n(i,1,1) = sp%un%n0*(1.0_rp + eps*cos(k*spc%g%x(i)))
      end do
      call spc%fwd(n, nhat)
      do it = 1, sp%fun%nterms
        mu = 0.0_rp; mu0 = 0.0_rp
        call sp%fun%terms(it)%t%add_mu(spc, n, nhat, mu)
        call sp%fun%terms(it)%t%add_mu(spc, n0a, n0hat, mu0)
        kfd = 2.0_rp*real(spc%mode(mu - mu0, m, 0, 0))/(eps*sp%un%n0)
        kan = sp%fun%terms(it)%t%kernel(k)
        err = abs(kfd - kan)/max(abs(kan), 1.0e-3_rp)
        call check(err < 2.0e-4_rp, trim(sp%fun%terms(it)%t%name)//' K(k) vs FD', err)
      end do
      call sp%fun%mu(spc, n, nhat, mu)
      call sp%fun%mu(spc, n0a, n0hat, mu0)
      kfd = 2.0_rp*real(spc%mode(mu - mu0, m, 0, 0))/(eps*sp%un%n0)
      kan = sp%fun%kernel(k)
      err = abs(kfd - kan)/abs(kan)
      call check(err < 2.0e-4_rp, 'total K(k) vs FD', err)
    end do
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
end program test_kernels

! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> The time-dependent model: assembles the right-hand sides of
!>   d_t n        = -div(m n u)/m
!>   d_t (m n u)  = -div(m n u u) - n grad mu[n] + n F_ext + div Pi
!>   d_t Pi       = -u.grad Pi + (Pi0 - Pi)/tau          (maxwell only)
!> and advances them with the low-storage third-order Runge-Kutta scheme of
!> Wray (as in CaNS). The stiff relaxation is applied per stage with the exact
!> solution of d_t Pi = (Pi0(t) - Pi)/tau for a source varying linearly over the
!> stage between Pi0(u_old) and Pi0(u_new):
!>   Pi <- Pi0_new - (Pi0_new - Pi0_old) (tau/h)(1 - e^{-h/tau}) + (Pi - Pi0_old) e^{-h/tau},
!> which is second-order accurate for tau >> h and exact (Pi = Pi0) as tau -> 0,
!> so tau is never a stability limit and the Newtonian closure is recovered.
module hydft_model
  use hydft_kinds
  use hydft_params
  use hydft_spectral
  use hydft_species
  use hydft_stress, only: ncomp
  use hydft_external_potential
  implicit none
  private
  public :: model_t

  type :: model_t
    type(params_t) :: p
    type(spec_t) :: spc
    type(species_t) :: sp
    type(vext_t) :: ext
    real(rp) :: t = 0.0_rp
    real(rp) :: dt = 0.0_rp
    integer  :: step = 0
    ! work arrays
    real(rp), allocatable :: u(:,:,:,:), mu(:,:,:), gmu(:,:,:,:), flux(:,:,:,:), divf(:,:,:)
    real(rp), allocatable :: pi(:,:,:,:), pi0(:,:,:,:), pi0_old(:,:,:,:), fpi(:,:,:,:), gpi(:,:,:,:)
    real(rp), allocatable :: dn(:,:,:), dmom(:,:,:,:), dpi(:,:,:,:)
    real(rp), allocatable :: dn_o(:,:,:), dmom_o(:,:,:,:), dpi_o(:,:,:,:)
    complex(rp), allocatable :: nhat(:,:,:)
  contains
    procedure :: init => model_init
    procedure :: rhs => model_rhs
    procedure :: compute_dt => model_compute_dt
    procedure :: advance => model_advance
    procedure :: velocity => model_velocity
  end type model_t

  real(rp), parameter :: rk_a(3) = [32.0_rp/60.0_rp, 25.0_rp/60.0_rp, 45.0_rp/60.0_rp]
  real(rp), parameter :: rk_b(3) = [0.0_rp, -17.0_rp/60.0_rp, -25.0_rp/60.0_rp]

contains

  subroutine model_init(m, p, verbose)
    class(model_t), intent(inout) :: m
    type(params_t), intent(in) :: p
    logical, intent(in), optional :: verbose
    integer :: nx, ny, nz
    m%p = p
    call m%spc%init(p%nx, p%ny, p%nz, p%lx, p%ly, p%lz, p%dealias)
    call m%sp%build(p, m%spc, verbose)
    call m%sp%allocate_state(m%spc)
    call m%ext%setup(m%spc, p%vext, p%vext_amplitude, p%vext_mode, p%vext_omega, p%vext_width, p%vext_ramp)
    nx = p%nx; ny = p%ny; nz = p%nz
    allocate(m%u(nx,ny,nz,3), m%mu(nx,ny,nz), m%gmu(nx,ny,nz,3), m%flux(nx,ny,nz,3), m%divf(nx,ny,nz))
    allocate(m%pi(nx,ny,nz,ncomp), m%pi0(nx,ny,nz,ncomp), m%pi0_old(nx,ny,nz,ncomp), m%fpi(nx,ny,nz,3), m%gpi(nx,ny,nz,3))
    allocate(m%dn(nx,ny,nz), m%dmom(nx,ny,nz,3), m%dpi(nx,ny,nz,ncomp))
    allocate(m%dn_o(nx,ny,nz), m%dmom_o(nx,ny,nz,3), m%dpi_o(nx,ny,nz,ncomp))
    allocate(m%nhat(m%spc%g%nkx,ny,nz))
    m%t = 0.0_rp; m%step = 0
    m%dt = p%dt
  end subroutine model_init

  subroutine model_velocity(m)
    class(model_t), intent(inout) :: m
    integer :: d
    do d = 1, 3
      m%u(:,:,:,d) = m%sp%mom(:,:,:,d)/(m%sp%mass*m%sp%n)
    end do
  end subroutine model_velocity

  !> Right-hand sides for the current state at time t.
  subroutine model_rhs(m, t, dn, dmom, dpi)
    class(model_t), intent(inout) :: m
    real(rp), intent(in) :: t
    real(rp), intent(out) :: dn(:,:,:), dmom(:,:,:,:), dpi(:,:,:,:)
    integer :: i, j, d, c
    associate (sp => m%sp, spc => m%spc)
      if (minval(sp%n) <= 0.0_rp) call fatal('density became non-positive; reduce the amplitude or the time step')
      call m%velocity()
      call spc%fwd(sp%n, m%nhat)
      ! continuity
      call spc%div(sp%mom, dn)
      dn = -dn/sp%mass
      ! chemical potential and its gradient
      call sp%fun%mu(spc, sp%n, m%nhat, m%mu)
      call spc%grad(m%mu, m%gmu)
      ! external force
      if (m%ext%active .and. m%ext%time_dependent) call m%ext%eval(spc, t)
      ! dissipative stress
      if (sp%stress%evolve) then
        m%pi = sp%pi
      else
        call sp%stress%pi0(spc, m%u, m%pi)
      end if
      call sp%stress%divergence(spc, m%pi, m%fpi)
      ! momentum
      do i = 1, 3
        if (i > spc%g%ndim) then
          dmom(:,:,:,i) = 0.0_rp
          cycle
        end if
        do j = 1, 3
          m%flux(:,:,:,j) = sp%mom(:,:,:,i)*m%u(:,:,:,j)
        end do
        call spc%div(m%flux, m%divf)
        ! the internal force -n grad mu has zero mean exactly (translation invariance
        ! of F[n]); remove the spectral-truncation residual so momentum is conserved
        m%divf = m%divf + sp%n*m%gmu(:,:,:,i)
        dmom(:,:,:,i) = -(m%divf - spc%mean(m%divf)) + m%fpi(:,:,:,i)
        if (m%ext%active) dmom(:,:,:,i) = dmom(:,:,:,i) + sp%n*m%ext%f(:,:,:,i)
      end do
      ! stress advection (relaxation is applied in the stage update)
      if (sp%stress%evolve) then
        do c = 1, ncomp
          call spc%grad(sp%pi(:,:,:,c), m%gpi)
          dpi(:,:,:,c) = 0.0_rp
          do d = 1, spc%g%ndim
            dpi(:,:,:,c) = dpi(:,:,:,c) - m%u(:,:,:,d)*m%gpi(:,:,:,d)
          end do
        end do
      end if
      ! dealias
      call spc%dealias(dn)
      do i = 1, spc%g%ndim
        call spc%dealias(dmom(:,:,:,i))
      end do
      if (sp%stress%evolve) then
        do c = 1, ncomp
          call spc%dealias(dpi(:,:,:,c))
        end do
      end if
    end associate
  end subroutine model_rhs

  !> Time step from acoustic/elastic (RK3 imaginary-axis limit sqrt3),
  !> viscous (real-axis limit 2.5) and advective constraints.
  subroutine model_compute_dt(m)
    class(model_t), intent(inout) :: m
    real(rp) :: wmax, k, c2, numax, umax, dtm, dxmin, kk
    integer :: i, j, l
    associate (sp => m%sp, spc => m%spc)
      if (m%p%dt > 0.0_rp) then
        m%dt = m%p%dt
        return
      end if
      wmax = 0.0_rp
      do l = 1, spc%g%nz; do j = 1, spc%g%ny; do i = 1, spc%g%nkx
        if (spc%g%mask(i,j,l) == 0.0_rp) cycle
        k = spc%g%kabs(i,j,l)
        if (k == 0.0_rp) cycle
        c2 = sp%un%n0*sp%fun%kernel(k)/sp%mass
        if (sp%stress%evolve) c2 = c2 + sp%tr%eta_l/(sp%mass*sp%stress%tau)
        kk = k*k*max(c2, 0.0_rp)
        wmax = max(wmax, sqrt(kk))
      end do; end do; end do
      numax = 0.0_rp
      if (.not. sp%stress%evolve) numax = sp%tr%eta_l/sp%mass*spc%g%kmax_dealiased**2
      call m%velocity()
      umax = maxval(abs(m%u))
      dxmin = spc%g%dx
      if (spc%g%ny > 1) dxmin = min(dxmin, spc%g%dy)
      if (spc%g%nz > 1) dxmin = min(dxmin, spc%g%dz)
      dtm = huge(1.0_rp)
      if (wmax > 0.0_rp) dtm = min(dtm, 1.7_rp/wmax)
      if (numax > 0.0_rp) dtm = min(dtm, 2.5_rp/numax)
      if (umax > 0.0_rp) dtm = min(dtm, dxmin/umax)
      if (dtm == huge(1.0_rp)) call fatal('cannot determine a time step; set dt in &time')
      m%dt = m%p%cfl*dtm
    end associate
  end subroutine model_compute_dt

  !> One RK3 step of size m%dt.
  subroutine model_advance(m)
    class(model_t), intent(inout) :: m
    integer :: s, i, c
    real(rp) :: dts, fac, phi, ts
    associate (sp => m%sp, spc => m%spc)
      ts = m%t
      do s = 1, 3
        call m%rhs(ts, m%dn, m%dmom, m%dpi)
        if (sp%stress%evolve) call sp%stress%pi0(spc, m%u, m%pi0_old)   ! m%u is the stage-start velocity
        if (s == 1) then
          sp%n = sp%n + m%dt*rk_a(1)*m%dn
          sp%mom = sp%mom + m%dt*rk_a(1)*m%dmom
          if (sp%stress%evolve) sp%pi = sp%pi + m%dt*rk_a(1)*m%dpi
        else
          sp%n = sp%n + m%dt*(rk_a(s)*m%dn + rk_b(s)*m%dn_o)
          sp%mom = sp%mom + m%dt*(rk_a(s)*m%dmom + rk_b(s)*m%dmom_o)
          if (sp%stress%evolve) sp%pi = sp%pi + m%dt*(rk_a(s)*m%dpi + rk_b(s)*m%dpi_o)
        end if
        m%dn_o = m%dn; m%dmom_o = m%dmom
        if (sp%stress%evolve) m%dpi_o = m%dpi
        dts = (rk_a(s) + rk_b(s))*m%dt
        ts = ts + dts
        ! exact relaxation towards the Newtonian stress, source linear over the stage
        if (sp%stress%evolve) then
          call m%velocity()
          call sp%stress%pi0(spc, m%u, m%pi0)
          fac = exp(-dts/sp%stress%tau)
          phi = (sp%stress%tau/dts)*(1.0_rp - fac)
          sp%pi = m%pi0 - (m%pi0 - m%pi0_old)*phi + (sp%pi - m%pi0_old)*fac
          do c = 1, ncomp
            call spc%dealias(sp%pi(:,:,:,c))
          end do
        end if
        call spc%dealias(sp%n)
        do i = 1, spc%g%ndim
          call spc%dealias(sp%mom(:,:,:,i))
        end do
      end do
      m%t = m%t + m%dt
      m%step = m%step + 1
    end associate
  end subroutine model_advance

end module hydft_model

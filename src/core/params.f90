! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Run-time parameters read from a single namelist input file (hydft.in),
!> in the spirit of CaNS' dns.in. Every physics choice is a string selected
!> here and turned into an object by the factories in the model layer.
module hydft_params
  use hydft_kinds
  use hydft_utils
  implicit none
  private
  public :: params_t

  type :: params_t
    ! &grid
    integer  :: nx = 64, ny = 1, nz = 1
    real(rp) :: lx = 20.0_rp, ly = 1.0_rp, lz = 1.0_rp
    logical  :: dealias = .true.
    ! &species
    character(len=16)  :: kind = 'ion'                    !< 'ion' (classical) | 'electron' (quantum)
    real(rp) :: gamma = 2.0_rp                            !< coupling Gamma
    real(rp) :: kappa = 0.1_rp                            !< screening kappa = a/lambda_D (ions)
    real(rp) :: rs = 1.86_rp                              !< density parameter (electrons)
    real(rp) :: theta = -1.0_rp                           !< degeneracy T/T_F; if > 0 overrides gamma (electrons)
    real(rp) :: mass = 1.0_rp
    character(len=64)  :: terms = ''                      !< comma list: ideal, tf, gradient, hartree, ry
    real(rp) :: tfk_gamma = 1.0_rp/9.0_rp                 !< gradient-correction prefactor (1/9 Kirzhnits, 1 vW)
    character(len=16)  :: closure = 'newtonian'           !< 'newtonian' | 'maxwell'
    character(len=16)  :: transport = 'constant'          !< 'constant' | 'yukawa_fit'
    real(rp) :: eta = 0.0_rp                              !< shear viscosity / (m n0 wp a^2)
    real(rp) :: xi = 0.0_rp                               !< bulk viscosity / (m n0 wp a^2)
    real(rp) :: tau = 0.0_rp                              !< Maxwell time * wp (constant model)
    character(len=16)  :: tau_model = 'ichimaru'          !< 'constant' | 'ichimaru' (Eq. 37)
    real(rp) :: adiabatic_index = 5.0_rp/3.0_rp
    character(len=16)  :: structure = 'hnc'               !< 'hnc' | 'file' | 'none'
    character(len=256) :: structure_file = ''
    character(len=1)   :: structure_kind = 's'            !< file holds 's' (S(k)) or 'c' (c(k))
    ! &hnc
    integer  :: hnc_nr = 8192, hnc_maxiter = 20000
    real(rp) :: hnc_rmax = 80.0_rp, hnc_tol = 1.0e-9_rp, hnc_mix = 0.2_rp
    ! &time
    real(rp) :: dt = 0.0_rp                               !< 0 -> from CFL
    real(rp) :: cfl = 0.4_rp
    real(rp) :: tend = 50.0_rp
    integer  :: nsteps_max = 10000000
    integer  :: iout = 0                                  !< field output interval (0: never)
    integer  :: iprobe = 1                                !< probe output interval
    integer  :: idiag = 10                                !< diagnostics interval
    ! &initial
    character(len=16) :: ic = 'mode'                      !< 'uniform' | 'mode' | 'random' | 'gaussian'
    real(rp) :: ic_amplitude = 1.0e-4_rp                  !< relative density amplitude
    integer  :: ic_mode(3) = [1, 0, 0]
    real(rp) :: ic_width = 1.0_rp
    integer  :: ic_seed = 1234
    ! &external
    character(len=16) :: vext = 'none'                    !< 'none' | 'mode' | 'driven' | 'gaussian'
    real(rp) :: vext_amplitude = 0.0_rp                   !< potential energy amplitude (units m wp^2 a^2)
    integer  :: vext_mode(3) = [1, 0, 0]
    real(rp) :: vext_omega = 0.0_rp
    real(rp) :: vext_width = 1.0_rp
    real(rp) :: vext_ramp = 0.0_rp                        !< switch-on time for the driven potential
    ! &probe
    integer  :: probe_mode(3) = [1, 0, 0]
    ! &linear
    real(rp) :: qmin = 0.05_rp, qmax = 2.5_rp
    integer  :: nq = 50
    real(rp) :: omega_max = 3.0_rp
    integer  :: nomega = 600
    real(rp) :: q_dsf(8) = [0.5_rp, 1.0_rp, 1.5_rp, -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp]
    ! &output
    character(len=256) :: prefix = 'out'
  contains
    procedure :: read => params_read
    procedure :: print => params_print
  end type params_t

contains

  subroutine params_read(p, path)
    class(params_t), intent(inout) :: p
    character(len=*), intent(in) :: path
    integer :: u, ios
    ! local mirrors for the namelists
    integer  :: nx, ny, nz
    real(rp) :: lx, ly, lz
    logical  :: dealias
    character(len=16)  :: kind, closure, transport, tau_model, structure
    real(rp) :: gamma, kappa, rs, theta, mass, tfk_gamma, eta, xi, tau, adiabatic_index
    character(len=64)  :: terms
    character(len=256) :: structure_file
    character(len=1)   :: structure_kind
    integer  :: nr, maxiter
    real(rp) :: rmax, tol, mix
    real(rp) :: dt, cfl, tend
    integer  :: nsteps_max, iout, iprobe, idiag
    character(len=16) :: type
    real(rp) :: amplitude, width, omega, ramp
    integer  :: mode(3), seed
    real(rp) :: qmin, qmax, omega_max, q_dsf(8)
    integer  :: nq, nomega
    character(len=256) :: prefix
    namelist /grid/ nx, ny, nz, lx, ly, lz, dealias
    namelist /species/ kind, gamma, kappa, rs, theta, mass, terms, tfk_gamma, closure, transport, eta, xi, tau, &
                       tau_model, adiabatic_index, structure, structure_file, structure_kind
    namelist /hnc/ nr, rmax, tol, mix, maxiter
    namelist /time/ dt, cfl, tend, nsteps_max, iout, iprobe, idiag
    namelist /initial/ type, amplitude, mode, width, seed
    namelist /external/ type, amplitude, mode, omega, width, ramp
    namelist /probe/ mode
    namelist /linear/ qmin, qmax, nq, omega_max, nomega, q_dsf
    namelist /output/ prefix

    u = open_file(path, 'read')

    nx = p%nx; ny = p%ny; nz = p%nz; lx = p%lx; ly = p%ly; lz = p%lz; dealias = p%dealias
    rewind(u); read(u, nml=grid, iostat=ios); call chk(ios, 'grid')
    p%nx = nx; p%ny = ny; p%nz = nz; p%lx = lx; p%ly = ly; p%lz = lz; p%dealias = dealias

    kind = p%kind; gamma = p%gamma; kappa = p%kappa; rs = p%rs; theta = p%theta; mass = p%mass
    terms = p%terms; tfk_gamma = p%tfk_gamma; closure = p%closure; transport = p%transport
    eta = p%eta; xi = p%xi; tau = p%tau; tau_model = p%tau_model; adiabatic_index = p%adiabatic_index
    structure = p%structure; structure_file = p%structure_file; structure_kind = p%structure_kind
    rewind(u); read(u, nml=species, iostat=ios); call chk(ios, 'species')
    p%kind = lower(kind); p%gamma = gamma; p%kappa = kappa; p%rs = rs; p%theta = theta; p%mass = mass
    p%terms = lower(terms); p%tfk_gamma = tfk_gamma; p%closure = lower(closure); p%transport = lower(transport)
    p%eta = eta; p%xi = xi; p%tau = tau; p%tau_model = lower(tau_model); p%adiabatic_index = adiabatic_index
    p%structure = lower(structure); p%structure_file = structure_file; p%structure_kind = lower(structure_kind)

    nr = p%hnc_nr; rmax = p%hnc_rmax; tol = p%hnc_tol; mix = p%hnc_mix; maxiter = p%hnc_maxiter
    rewind(u); read(u, nml=hnc, iostat=ios); call chk(ios, 'hnc')
    p%hnc_nr = nr; p%hnc_rmax = rmax; p%hnc_tol = tol; p%hnc_mix = mix; p%hnc_maxiter = maxiter

    dt = p%dt; cfl = p%cfl; tend = p%tend; nsteps_max = p%nsteps_max; iout = p%iout; iprobe = p%iprobe; idiag = p%idiag
    rewind(u); read(u, nml=time, iostat=ios); call chk(ios, 'time')
    p%dt = dt; p%cfl = cfl; p%tend = tend; p%nsteps_max = nsteps_max; p%iout = iout; p%iprobe = iprobe; p%idiag = idiag

    type = p%ic; amplitude = p%ic_amplitude; mode = p%ic_mode; width = p%ic_width; seed = p%ic_seed
    rewind(u); read(u, nml=initial, iostat=ios); call chk(ios, 'initial')
    p%ic = lower(type); p%ic_amplitude = amplitude; p%ic_mode = mode; p%ic_width = width; p%ic_seed = seed

    type = p%vext; amplitude = p%vext_amplitude; mode = p%vext_mode; omega = p%vext_omega; width = p%vext_width; ramp = p%vext_ramp
    rewind(u); read(u, nml=external, iostat=ios); call chk(ios, 'external')
    p%vext = lower(type); p%vext_amplitude = amplitude; p%vext_mode = mode; p%vext_omega = omega; p%vext_width = width; p%vext_ramp = ramp

    mode = p%probe_mode
    rewind(u); read(u, nml=probe, iostat=ios); call chk(ios, 'probe')
    p%probe_mode = mode

    qmin = p%qmin; qmax = p%qmax; nq = p%nq; omega_max = p%omega_max; nomega = p%nomega; q_dsf = p%q_dsf
    rewind(u); read(u, nml=linear, iostat=ios); call chk(ios, 'linear')
    p%qmin = qmin; p%qmax = qmax; p%nq = nq; p%omega_max = omega_max; p%nomega = nomega; p%q_dsf = q_dsf

    prefix = p%prefix
    rewind(u); read(u, nml=output, iostat=ios); call chk(ios, 'output')
    p%prefix = prefix
    close(u)

    if (len_trim(p%terms) == 0) then
      if (p%kind == 'electron') then
        p%terms = 'tf,gradient,hartree,ry'
      else
        p%terms = 'ideal,hartree,ry'
      end if
    end if
    if (p%kind == 'electron' .and. p%theta > 0.0_rp) then
      ! Gamma = 2 rs / (theta (9 pi/4)^(2/3))
      p%gamma = 2.0_rp*p%rs/(p%theta*(9.0_rp*pi/4.0_rp)**(2.0_rp/3.0_rp))
    end if
    if (p%kind == 'electron') p%kappa = 0.0_rp
  contains
    subroutine chk(ios, name)
      integer, intent(in) :: ios
      character(len=*), intent(in) :: name
      if (ios > 0) call fatal('error reading namelist &'//name//' in '//trim(path))
      ! ios < 0: namelist absent -> keep defaults
    end subroutine chk
  end subroutine params_read

  subroutine params_print(p, unit)
    class(params_t), intent(in) :: p
    integer, intent(in), optional :: unit
    integer :: u
    u = stdout
    if (present(unit)) u = unit
    write(u,'(a)') 'HyDFT parameters'
    write(u,'(a,3i6,a,3f10.4)') '  grid      nx,ny,nz =', p%nx, p%ny, p%nz, '   lx,ly,lz =', p%lx, p%ly, p%lz
    write(u,'(a,a,a,f9.4,a,f7.3,a,f7.3)') '  species   ', trim(p%kind), '  Gamma =', p%gamma, '  kappa =', p%kappa, '  rs =', p%rs
    write(u,'(a,a)') '  terms     ', trim(p%terms)
    write(u,'(a,a,a,a,a,a)') '  closure   ', trim(p%closure), '   transport ', trim(p%transport), '   structure ', trim(p%structure)
    write(u,'(a,3es11.3)') '  eta,xi,tau (input) =', p%eta, p%xi, p%tau
    write(u,'(a,es10.2,a,f8.3,a,f10.3)') '  time      dt =', p%dt, '  cfl =', p%cfl, '  tend =', p%tend
    write(u,'(a,a,a,es10.2,a,3i4)') '  initial   ', trim(p%ic), '  amplitude =', p%ic_amplitude, '  mode =', p%ic_mode
    write(u,'(a,a,a,es10.2,a,3i4,a,f8.4)') '  external  ', trim(p%vext), '  amplitude =', p%vext_amplitude, '  mode =', p%vext_mode, '  omega =', p%vext_omega
  end subroutine params_print

end module hydft_params

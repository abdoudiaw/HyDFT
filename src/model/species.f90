! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> A species: its units, structure, free-energy functional, transport,
!> stress closure and state arrays (n, momentum, dissipative stress).
!> `build` is the factory that turns the input parameters into objects; it is
!> shared by the simulation, the linear-response tool and the tests.
module hydft_species
  use hydft_kinds
  use hydft_utils
  use hydft_params
  use hydft_units
  use hydft_spectral
  use hydft_pair_potential
  use hydft_structure_source
  use hydft_free_energy
  use hydft_fe_ideal_classical
  use hydft_fe_kinetic_tfk
  use hydft_fe_hartree
  use hydft_fe_correlation_ry
  use hydft_functional
  use hydft_transport
  use hydft_stress
  implicit none
  private
  public :: species_t

  type :: species_t
    character(len=16) :: kind = 'ion'
    type(units_t) :: un
    type(structure_t) :: st
    type(functional_t) :: fun
    type(transport_t) :: tr
    type(stress_t) :: stress
    real(rp) :: mass = 1.0_rp
    real(rp), allocatable :: n(:,:,:)        !< density
    real(rp), allocatable :: mom(:,:,:,:)    !< momentum density m n u (3 components)
    real(rp), allocatable :: pi(:,:,:,:)     !< dissipative stress (6 components; allocated for maxwell)
  contains
    procedure :: build => species_build
    procedure :: allocate_state => species_allocate_state
  end type species_t

contains

  subroutine species_build(sp, p, spc, verbose)
    class(species_t), intent(inout) :: sp
    type(params_t), intent(in) :: p
    type(spec_t), intent(inout) :: spc
    logical, intent(in), optional :: verbose
    class(pair_potential_t), allocatable :: pot
    type(fe_ideal_t) :: ideal
    type(fe_tf_t) :: tf
    type(fe_gradient_t) :: grad
    type(fe_hartree_t) :: hartree
    type(fe_ry_t) :: ry
    character(len=64) :: list
    character(len=16) :: item
    integer :: ipos
    logical :: verb, need_deriv, has_hartree

    verb = .true.
    if (present(verbose)) verb = verbose
    sp%kind = p%kind
    sp%mass = p%mass
    if (allocated(sp%fun%terms)) deallocate(sp%fun%terms)
    sp%fun%nterms = 0
    call sp%un%init(p%kind, p%gamma, p%kappa, p%rs, p%mass)
    if (verb) call sp%un%print()

    ! ---- structure
    need_deriv = (p%closure == 'maxwell' .and. p%tau_model == 'ichimaru')
    select case (trim(p%structure))
    case ('hnc')
      if (sp%un%quantum) then
        call make_pair_potential('hansen_mcdonald', sp%un%gamma, 0.0_rp, sp%un%lambda_qsp, pot)
      else
        call make_pair_potential('yukawa', sp%un%gamma, sp%un%kappa, 1.0_rp, pot)
      end if
      if (verb) write(*,'(a,a,a)') '  structure HNC with ', trim(pot%name), ' potential ...'
      call structure_from_hnc(pot, sp%un%n0, p%hnc_nr, p%hnc_rmax, p%hnc_tol, p%hnc_mix, p%hnc_maxiter, sp%st, &
                              with_derivative=need_deriv)
      if (verb) write(*,'(a,f10.5,a,f10.5)') '            u_ex/NkT =', sp%st%u_ex, '   du/dGamma =', sp%st%du_dgamma
    case ('file')
      if (len_trim(p%structure_file) == 0) call fatal('structure = file needs structure_file')
      call structure_from_file(p%structure_file, p%structure_kind, sp%un%n0, sp%st)
      if (verb) write(*,'(a,a)') '  structure read from ', trim(p%structure_file)
    case ('none')
      sp%st%source = 'none'
    case default
      call fatal('unknown structure source: '//trim(p%structure))
    end select

    ! ---- functional
    list = p%terms
    if (len_trim(list) == 0) then
      if (sp%un%quantum) then
        list = 'tf,gradient,hartree,ry'
      else
        list = 'ideal,hartree,ry'
      end if
    end if
    has_hartree = index(list, 'hartree') > 0
    do
      if (len_trim(list) == 0) exit
      ipos = index(list, ',')
      if (ipos == 0) then
        item = adjustl(list); list = ''
      else
        item = adjustl(list(1:ipos-1)); list = list(ipos+1:)
      end if
      select case (trim(item))
      case ('ideal')
        ideal%n0 = sp%un%n0; ideal%beta = sp%un%beta; ideal%name = 'ideal'
        call sp%fun%add(ideal)
      case ('tf')
        if (.not. sp%un%quantum) call fatal('tf term needs kind = electron')
        tf%n0 = sp%un%n0; tf%beta = sp%un%beta; tf%cn = sp%un%cn; tf%alpha0 = sp%un%alpha0; tf%name = 'tf'
        call sp%fun%add(tf)
      case ('gradient')
        if (.not. sp%un%quantum) call fatal('gradient term needs kind = electron')
        grad%n0 = sp%un%n0; grad%beta = sp%un%beta; grad%cn = sp%un%cn; grad%alpha0 = sp%un%alpha0
        grad%gamma_g = p%tfk_gamma; grad%name = 'gradient'
        grad%cg = p%tfk_gamma*(3.0_rp*sqrt2*pi*pi/8.0_rp)*sp%un%hbar**5*sp%un%beta**1.5_rp/sp%un%mass**2.5_rp
        call sp%fun%add(grad)
      case ('hartree')
        hartree%n0 = sp%un%n0; hartree%beta = sp%un%beta; hartree%e2 = sp%un%e2; hartree%kappa = sp%un%kappa
        hartree%name = 'hartree'
        call sp%fun%add(hartree)
      case ('ry')
        if (sp%st%source == 'none') call fatal('ry term needs a structure source (hnc or file)')
        ry%n0 = sp%un%n0; ry%beta = sp%un%beta; ry%e2 = sp%un%e2; ry%kappa = sp%un%kappa
        ry%subtract_mean_field = has_hartree
        ry%st = sp%st
        ry%name = 'ry'
        call sp%fun%add(ry)
      case ('')
      case default
        call fatal('unknown free-energy term: '//trim(item))
      end select
    end do
    sp%fun%n0 = sp%un%n0
    call sp%fun%setup(spc)

    ! ---- transport and closure
    call sp%tr%setup(p%transport, p%tau_model, p%eta, p%xi, p%tau, sp%un%gamma, sp%un%kappa, p%adiabatic_index, &
                     sp%st, p%closure == 'maxwell')
    if (verb) call sp%tr%print()
    call sp%stress%setup(p%closure, sp%tr%eta, sp%tr%xi, sp%tr%tau, sp%un%n0, spc)
  end subroutine species_build

  subroutine species_allocate_state(sp, spc)
    class(species_t), intent(inout) :: sp
    type(spec_t), intent(in) :: spc
    integer :: nx, ny, nz
    nx = spc%g%nx; ny = spc%g%ny; nz = spc%g%nz
    if (allocated(sp%n)) deallocate(sp%n, sp%mom)
    allocate(sp%n(nx,ny,nz), sp%mom(nx,ny,nz,3))
    sp%n = sp%un%n0
    sp%mom = 0.0_rp
    if (allocated(sp%pi)) deallocate(sp%pi)
    if (sp%stress%evolve) then
      allocate(sp%pi(nx,ny,nz,ncomp))
      sp%pi = 0.0_rp
    end if
  end subroutine species_allocate_state

end module hydft_species

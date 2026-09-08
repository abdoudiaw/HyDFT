! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> The file structure source reproduces the HNC source when fed its own output.
program test_structure_file
  use hydft_kinds
  use hydft_pair_potential
  use hydft_structure_source
  use hydft_utils
  implicit none
  class(pair_potential_t), allocatable :: pot
  type(structure_t) :: st, sf, sc
  real(rp) :: n0, q, err, emax
  integer :: u, i, nfail
  nfail = 0
  n0 = 3.0_rp/(4.0_rp*pi)
  call make_pair_potential('yukawa', 10.0_rp, 0.5_rp, 1.0_rp, pot)
  call structure_from_hnc(pot, n0, 4096, 60.0_rp, 1e-9_rp, 0.2_rp, 5000, st)
  ! write S(k) with a comment line, read back as 's'
  u = open_file('test_structure_sk.dat', 'write')
  write(u,'(a)') '# k  S(k) from HNC'
  do i = 1, size(st%k)
    write(u,'(2es22.12)') st%k(i), st%s(i)
  end do
  close(u)
  call structure_from_file('test_structure_sk.dat', 's', n0, sf)
  u = open_file('test_structure_ck.dat', 'write')
  do i = 1, size(st%k)
    write(u,'(2es22.12)') st%k(i), st%c(i)
  end do
  close(u)
  call structure_from_file('test_structure_ck.dat', 'c', n0, sc)
  emax = 0.0_rp
  do i = 1, 200
    q = 0.1_rp + 0.05_rp*i
    err = abs(sf%c_of_k(q) - st%c_of_k(q)); emax = max(emax, err)
    err = abs(sc%c_of_k(q) - st%c_of_k(q)); emax = max(emax, err)
  end do
  call check(emax < 1e-6_rp, 'file source c(k) == hnc source c(k)', emax)
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_structure_file: ', nfail, ' failure(s)'; error stop 1
  end if
  write(*,'(a)') 'test_structure_file: all passed'
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
end program test_structure_file

! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Field output: raw little-endian stream files <prefix>_fld_NNNNNN.bin with a
!> small header (nx, ny, nz, ndim, ncomp_pi, t), then n, mom(1:3), and Pi if
!> present. python/hydft_io.py reads them.
module hydft_io
  use hydft_kinds
  use hydft_model
  implicit none
  private
  public :: write_fields

contains

  subroutine write_fields(m, prefix)
    type(model_t), intent(in) :: m
    character(len=*), intent(in) :: prefix
    character(len=512) :: fname
    integer :: u, ios, npi
    write(fname, '(a,a,i6.6,a)') trim(prefix), '_fld_', m%step, '.bin'
    open(newunit=u, file=fname, access='stream', form='unformatted', status='replace', action='write', iostat=ios)
    if (ios /= 0) call fatal('cannot write '//trim(fname))
    npi = 0
    if (allocated(m%sp%pi)) npi = size(m%sp%pi, 4)
    write(u) m%spc%g%nx, m%spc%g%ny, m%spc%g%nz, m%spc%g%ndim, npi
    write(u) m%t, m%spc%g%lx, m%spc%g%ly, m%spc%g%lz
    write(u) m%sp%n
    write(u) m%sp%mom
    if (npi > 0) write(u) m%sp%pi
    close(u)
  end subroutine write_fields

end module hydft_io

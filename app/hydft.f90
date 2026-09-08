! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Time-dependent HyDFT simulation:  hydft input.in
program hydft_main
  use hydft
  implicit none
  character(len=256) :: infile
  if (command_argument_count() < 1) then
    write(*,'(a)') 'usage: hydft input.in'
    stop
  end if
  call get_command_argument(1, infile)
  call hydft_run(infile)
end program hydft_main

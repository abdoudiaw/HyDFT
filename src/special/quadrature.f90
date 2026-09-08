! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Gauss-Legendre nodes and weights (Newton iteration on Legendre polynomials).
module hydft_quadrature
  use hydft_kinds
  implicit none
  private
  public :: gauss_legendre

contains

  !> Nodes x(1:n) in (-1,1) and weights w(1:n).
  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(rp), intent(out) :: x(n), w(n)
    integer :: i, j, it
    real(rp) :: z, z1, p1, p2, p3, pp
    do i = 1, (n+1)/2
      z = cos(pi*(i - 0.25_rp)/(n + 0.5_rp))
      do it = 1, 100
        p1 = 1.0_rp; p2 = 0.0_rp
        do j = 1, n
          p3 = p2; p2 = p1
          p1 = ((2.0_rp*j - 1.0_rp)*z*p2 - (j - 1.0_rp)*p3)/j
        end do
        pp = n*(z*p1 - p2)/(z*z - 1.0_rp)
        z1 = z
        z = z1 - p1/pp
        if (abs(z - z1) < 1.0e-15_rp) exit
      end do
      x(i) = -z; x(n+1-i) = z
      w(i) = 2.0_rp/((1.0_rp - z*z)*pp*pp); w(n+1-i) = w(i)
    end do
  end subroutine gauss_legendre

end module hydft_quadrature

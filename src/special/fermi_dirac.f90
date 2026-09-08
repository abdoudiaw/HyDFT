! -
!
! SPDX-FileCopyrightText: Abdou Diaw and the HyDFT contributors
! SPDX-License-Identifier: MIT
!
! -
!> Fermi-Dirac integrals I_p(alpha) = int_0^inf x^p/(1+exp(x-alpha)) dx for
!> p = -1/2, 1/2, 3/2, the alpha-derivatives of I_{-1/2}, and the inverse of
!> I_{1/2}. Values are tabulated once on a fine alpha grid by Gauss-Legendre
!> quadrature and interpolated with cubic Hermite splines that use the exact
!> derivative relations; outside the table the Sommerfeld (alpha large) and
!> classical (alpha very negative) expansions are used.
module hydft_fermi_dirac
  use hydft_kinds
  use hydft_quadrature
  implicit none
  private
  public :: fd_init, fd_im12, fd_ip12, fd_ip32, fd_dim12, fd_d2im12, fd_alpha_from_ip12, fd_xi, fd_dxi

  real(rp), parameter :: amin = -30.0_rp, amax = 300.0_rp, h = 0.05_rp
  integer,  parameter :: ntab = nint((amax - amin)/h) + 1
  integer,  parameter :: ngl = 20
  ! table columns: 1 I_{-1/2}, 2 I_{1/2}, 3 I_{3/2}, 4 I'_{-1/2}, 5 I''_{-1/2}, 6 I'''_{-1/2}
  real(rp), allocatable, save :: tab(:,:)
  logical, save :: ready = .false.

contains

  !> Build the table (idempotent).
  subroutine fd_init()
    real(rp) :: xg(ngl), wg(ngl), alpha, xmax, a, b, xm, xr, x, t, f, e, s(6), p1, p2, p3
    integer :: i, ip, npan, q
    if (ready) return
    allocate(tab(6, ntab))
    call gauss_legendre(ngl, xg, wg)
    do i = 1, ntab
      alpha = amin + (i-1)*h
      xmax = max(alpha, 0.0_rp) + 45.0_rp
      s = 0.0_rp
      ! first panel x in [0,1] with x = t^2 to remove the x^{-1/2} singularity
      do q = 1, ngl
        t = 0.5_rp*(xg(q) + 1.0_rp)
        x = t*t
        call fermi_derivs(x - alpha, f, p1, p2, p3)
        e = wg(q)*0.5_rp*2.0_rp           ! dx = 2 t dt, and x^p = t^{2p}
        s(1) = s(1) + e*f                 ! x^{-1/2} * 2t = 2
        s(2) = s(2) + e*t*t*f             ! x^{1/2} * 2t = 2 t^2
        s(3) = s(3) + e*t**4*f            ! x^{3/2} * 2t = 2 t^4
        s(4) = s(4) + e*p1
        s(5) = s(5) + e*p2
        s(6) = s(6) + e*p3
      end do
      ! remaining panels of unit width on [1, xmax]
      npan = ceiling(xmax - 1.0_rp)
      do ip = 1, npan
        a = 1.0_rp + (ip-1); b = min(a + 1.0_rp, xmax)
        xm = 0.5_rp*(a + b); xr = 0.5_rp*(b - a)
        do q = 1, ngl
          x = xm + xr*xg(q)
          call fermi_derivs(x - alpha, f, p1, p2, p3)
          e = wg(q)*xr
          t = 1.0_rp/sqrt(x)
          s(1) = s(1) + e*t*f
          s(2) = s(2) + e*sqrt(x)*f
          s(3) = s(3) + e*x*sqrt(x)*f
          s(4) = s(4) + e*t*p1
          s(5) = s(5) + e*t*p2
          s(6) = s(6) + e*t*p3
        end do
      end do
      tab(:, i) = s
    end do
    ready = .true.
  end subroutine fd_init

  !> Fermi function f = 1/(1+e^y) with y = x - alpha and its alpha-derivatives
  !> d/dalpha f = f(1-f), d2 = f(1-f)(1-2f), d3 = f(1-f)(1-6f+6f^2).
  pure subroutine fermi_derivs(y, f, p1, p2, p3)
    real(rp), intent(in) :: y
    real(rp), intent(out) :: f, p1, p2, p3
    real(rp) :: g
    if (y > 0.0_rp) then
      g = exp(-y); f = g/(1.0_rp + g)
    else
      g = exp(y); f = 1.0_rp/(1.0_rp + g)
    end if
    p1 = f*(1.0_rp - f)
    p2 = p1*(1.0_rp - 2.0_rp*f)
    p3 = p1*(1.0_rp - 6.0_rp*f + 6.0_rp*f*f)
  end subroutine fermi_derivs

  !> Cubic Hermite interpolation of column `col` whose alpha-derivative is column `dcol`
  !> (or, if dcol == 0, the derivative is given by relation p*I_{p-1}: handled by caller).
  real(rp) function hermite(col, alpha, deriv) result(v)
    integer, intent(in) :: col
    real(rp), intent(in) :: alpha
    real(rp), intent(in) :: deriv(2)   ! derivative at the two bracketing nodes
    integer :: i
    real(rp) :: t, t2, t3, y0, y1, h00, h10, h01, h11
    i = int((alpha - amin)/h) + 1
    i = min(max(i, 1), ntab - 1)
    t = (alpha - (amin + (i-1)*h))/h
    t2 = t*t; t3 = t2*t
    h00 = 2*t3 - 3*t2 + 1; h10 = t3 - 2*t2 + t; h01 = -2*t3 + 3*t2; h11 = t3 - t2
    y0 = tab(col, i); y1 = tab(col, i+1)
    v = h00*y0 + h10*h*deriv(1) + h01*y1 + h11*h*deriv(2)
  end function hermite

  integer function node(alpha)
    real(rp), intent(in) :: alpha
    node = int((alpha - amin)/h) + 1
    node = min(max(node, 1), ntab - 1)
  end function node

  ! ---- Sommerfeld expansion (alpha >> 1): m-th alpha-derivative of I_p
  pure real(rp) function sommerfeld(p, alpha, m) result(v)
    real(rp), intent(in) :: p, alpha
    integer, intent(in) :: m
    real(rp) :: c(0:2), ex, coef
    integer :: j, l
    c(0) = 1.0_rp
    c(1) = (pi**2/6.0_rp)*p*(p + 1.0_rp)
    c(2) = (7.0_rp*pi**4/360.0_rp)*(p + 1.0_rp)*p*(p - 1.0_rp)*(p - 2.0_rp)
    v = 0.0_rp
    do j = 0, 2
      ex = p + 1.0_rp - 2*j
      coef = c(j)/(p + 1.0_rp)
      do l = 0, m - 1
        coef = coef*(ex - l)
      end do
      v = v + coef*alpha**(ex - m)
    end do
  end function sommerfeld

  ! ---- classical expansion (alpha << -1): m-th derivative of I_p
  pure real(rp) function classical(p, alpha, m) result(v)
    real(rp), intent(in) :: p, alpha
    integer, intent(in) :: m
    integer :: j
    real(rp) :: sgn
    v = 0.0_rp; sgn = 1.0_rp
    do j = 1, 12
      v = v + sgn*exp(j*alpha)*real(j, rp)**(m - p - 1.0_rp)
      sgn = -sgn
    end do
    v = v*gamma(p + 1.0_rp)
  end function classical

  !> I_{-1/2}(alpha)
  real(rp) function fd_im12(alpha) result(v)
    real(rp), intent(in) :: alpha
    integer :: i
    call fd_init()
    if (alpha > amax) then
      v = sommerfeld(-0.5_rp, alpha, 0)
    else if (alpha < amin) then
      v = classical(-0.5_rp, alpha, 0)
    else
      i = node(alpha)
      v = hermite(1, alpha, [tab(4,i), tab(4,i+1)])
    end if
  end function fd_im12

  !> I_{1/2}(alpha)
  real(rp) function fd_ip12(alpha) result(v)
    real(rp), intent(in) :: alpha
    integer :: i
    call fd_init()
    if (alpha > amax) then
      v = sommerfeld(0.5_rp, alpha, 0)
    else if (alpha < amin) then
      v = classical(0.5_rp, alpha, 0)
    else
      i = node(alpha)
      v = hermite(2, alpha, 0.5_rp*[tab(1,i), tab(1,i+1)])
    end if
  end function fd_ip12

  !> I_{3/2}(alpha)
  real(rp) function fd_ip32(alpha) result(v)
    real(rp), intent(in) :: alpha
    integer :: i
    call fd_init()
    if (alpha > amax) then
      v = sommerfeld(1.5_rp, alpha, 0)
    else if (alpha < amin) then
      v = classical(1.5_rp, alpha, 0)
    else
      i = node(alpha)
      v = hermite(3, alpha, 1.5_rp*[tab(2,i), tab(2,i+1)])
    end if
  end function fd_ip32

  !> dI_{-1/2}/dalpha
  real(rp) function fd_dim12(alpha) result(v)
    real(rp), intent(in) :: alpha
    integer :: i
    call fd_init()
    if (alpha > amax) then
      v = sommerfeld(-0.5_rp, alpha, 1)
    else if (alpha < amin) then
      v = classical(-0.5_rp, alpha, 1)
    else
      i = node(alpha)
      v = hermite(4, alpha, [tab(5,i), tab(5,i+1)])
    end if
  end function fd_dim12

  !> d2I_{-1/2}/dalpha2
  real(rp) function fd_d2im12(alpha) result(v)
    real(rp), intent(in) :: alpha
    integer :: i
    call fd_init()
    if (alpha > amax) then
      v = sommerfeld(-0.5_rp, alpha, 2)
    else if (alpha < amin) then
      v = classical(-0.5_rp, alpha, 2)
    else
      i = node(alpha)
      v = hermite(5, alpha, [tab(6,i), tab(6,i+1)])
    end if
  end function fd_d2im12

  !> xi(alpha) = I'_{-1/2}/I_{-1/2}^2, the finite-temperature gradient-correction weight (Sci. Rep. Eq. 18).
  real(rp) function fd_xi(alpha) result(v)
    real(rp), intent(in) :: alpha
    real(rp) :: im
    im = fd_im12(alpha)
    v = fd_dim12(alpha)/(im*im)
  end function fd_xi

  !> d xi / d alpha
  real(rp) function fd_dxi(alpha) result(v)
    real(rp), intent(in) :: alpha
    real(rp) :: im, d1, d2
    im = fd_im12(alpha); d1 = fd_dim12(alpha); d2 = fd_d2im12(alpha)
    v = (d2*im - 2.0_rp*d1*d1)/im**3
  end function fd_dxi

  !> Inverse of I_{1/2}: alpha such that I_{1/2}(alpha) = y (y > 0). Safeguarded Newton.
  real(rp) function fd_alpha_from_ip12(y) result(alpha)
    real(rp), intent(in) :: y
    real(rp) :: f, df, lo, hi, step
    integer :: it
    if (y <= 0.0_rp) call fatal('fermi_dirac: I_{1/2} inverse needs y > 0')
    ! initial guess from the two limits
    alpha = max(log(2.0_rp*y/sqrtpi), (1.5_rp*y)**(2.0_rp/3.0_rp) - 1.0_rp)
    if (y < 1.0_rp) alpha = log(2.0_rp*y/sqrtpi)
    lo = -huge(1.0_rp); hi = huge(1.0_rp)
    do it = 1, 100
      f = fd_ip12(alpha) - y
      if (f > 0.0_rp) then
        hi = min(hi, alpha)
      else
        lo = max(lo, alpha)
      end if
      df = 0.5_rp*fd_im12(alpha)
      step = f/df
      if (abs(step) < 1.0e-14_rp*max(1.0_rp, abs(alpha))) exit
      alpha = alpha - step
      if (alpha <= lo .or. alpha >= hi) then
        if (lo > -huge(1.0_rp) .and. hi < huge(1.0_rp)) then
          alpha = 0.5_rp*(lo + hi)
        end if
      end if
    end do
  end function fd_alpha_from_ip12

end module hydft_fermi_dirac

program laguerre_int
  implicit none

  integer, parameter:: N = 4
  integer, parameter:: maxita = 10
  double precision, parameter:: pi = acos(-1.0d0)

  integer:: i
  double precision:: ii, Is, Ie, h, theta, fi
  double precision, allocatable:: x(:), w(:)

!init
  allocate(x(N), w(N))

  call compute_intpoint(Ln, N, maxita, x, w)

  do i=1, N
    ii = w(i)*f(x(i))
  end do

  write(*, *) "-------------------------"
  write(*, *) "I = ", ii

  stop

contains

  subroutine compute_intpoint(Ln, n, maxita, x, w)
    implicit none

    interface
      function Ln(x, nn) result(z)
        double precision, intent(in):: x
        integer, intent(in):: nn
        double precision:: z
      end function Ln
    end interface

    integer, intent(in):: n, maxita
    double precision, intent(inout):: x(:), w(:)

    integer:: i, ita
    double precision:: s, deltax, ni

  !init
    do i=1, n
      x(i) = (pi**2 * ((i-1) + 0.75)**2)/(4*n)
    end do
  !mainloop
    write(*, *) "point"
    do ita=1, maxita
      s = 0
      do i=1, n
        deltax = -(x(i)*Ln(x(i), n))/((x(i) - n - 1)*Ln(x(i), n) + Ln(x(i), n+1))
        x(i) = x(i) + deltax
        s = s + x(i)
      end do

      write(*, *) "-------------------------"
      write(*, *) "x", x
      write(*, *) "s", s
    end do

    ni = 1
    do i=1, n
      ni = ni*i
    end do

    do i=1, n
      w(i) = (ni**2 * x(i))/(Ln(x(i), n+1)**2)
    end do

    write(*, *) "-------------------------"
    write(*, *) "w", w

  end subroutine compute_intpoint

  function f(x) result(y)
    implicit none

    double precision, intent(in):: x
    double precision:: y

    y = sin(x**2)

    return

  end function f

  recursive function Ln(x, nn) result(z)
    implicit none
    double precision, intent(in):: x
    integer, intent(in):: nn
    double precision:: Ln_1, Ln_2
    double precision:: z

    if(nn == 0) then
      z = 1
    else if(nn == 1) then
      z = 1 - x
    else
      Ln_1 = Ln(x, nn-1)
      Ln_2 = Ln(x, nn-2)
      z = (2*nn - 1 - x)*Ln_1 - ((nn - 1)**2)*Ln_2
    end if

    return
  end function Ln

end program laguerre_int

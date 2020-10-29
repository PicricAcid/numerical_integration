program hermit_int
  implicit none

  integer, parameter:: N = 4
  integer, parameter:: maxita = 10
  double precision, parameter:: pi = acos(-1.0d0)

  integer:: i
  double precision:: ii, Is, Ie, h, theta, fi
  double precision, allocatable:: x(:), w(:)

!init
  allocate(x(N), w(N))

  call compute_intpoint(Hn, N, maxita, x, w)

  do i=1, N
    ii = w(i)*f(x(i))
  end do

  write(*, *) "-------------------------"
  write(*, *) "I = ", ii

  stop

contains

  subroutine compute_intpoint(Hn, n, maxita, x, w)
    implicit none

    interface
      function Hn(x, nn) result(z)
        double precision, intent(in):: x
        integer, intent(in):: nn
        double precision:: z
      end function Hn
    end interface

    integer, intent(in):: n, maxita
    double precision, intent(inout):: x(:), w(:)

    integer:: i, ita
    double precision:: s, deltax, ni, insq

  !init
    do i=1, n
      if(mod(i, 2) == 0) then
        insq = 4*n + 1
        x(i) = (1/sqrt(insq))*((2*i - 1)/2)*pi
      else
        insq = 4*n + 3
        x(i) = (1/sqrt(insq))*i*pi
      end if
    end do
  !mainloop
    write(*, *) "point"
    do ita=1, maxita
      s = 0
      do i=1, n
        deltax = -Hn(x(i), n)/(2*x(i)*Hn(x(i), n) - Hn(x(i), n+1))
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
      w(i) = (2**(n+1) * ni * sqrt(pi))/(Hn(x(i), n+1)**2)
    end do

    write(*, *) "-------------------------"
    write(*, *) "w", w

  end subroutine compute_intpoint

  function f(x) result(y)
    implicit none

    double precision, intent(in):: x
    double precision:: y

    y = 1

    return

  end function f

  recursive function Hn(x, nn) result(z)
    implicit none
    double precision, intent(in):: x
    integer, intent(in):: nn
    double precision:: Hn_1, Hn_2
    double precision:: z

    if(nn == 0) then
      z = 1
    else if(nn == 1) then
      z = 2*x
    else
      Hn_1 = Hn(x, nn-1)
      Hn_2 = Hn(x, nn-2)
      z = 2*x*Hn_1 - 2*(nn - 1)*Hn_2
    end if

    return
  end function Hn

end program hermit_int

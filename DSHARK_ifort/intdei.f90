!> Performs numerical integration of improper integrals using a double exponential quadrature method, written by Takuya Ooura, Reasearch Institute for Mathematical Science, Kyoto University (See http:\\www.kurims.kyoto-u.ac.jp/~ooura/intde.html).
!! \param f function to be integrated
!! \param n parameter which is supplied to f
!! \param h1 parameter which is supplied to f
!! \param h2 parameter which is supplied to f
!! \param h3 parameter which is supplied to f
!! \param case parameter which is supplied to f
!! \param a lower limit of integration
!! \param aw points and weights of the quadrature formula
!! \param i approximation to the integral
!! \param err estimate of the absolute error
subroutine intdei(f, n,h1,h2,h3,case, a, aw, i, err)
  ! intdei
  !     [description]
  !         I = integral of f(x) over (a,infinity), 
  !             f(x) has not oscillatory factor.
  !     [declaration]
  !         external f
  !     [usage]
  !         call intdeiini(lenaw, tiny, eps, aw)  ! initialization of aw
  !         ...
  !         call intdei(f, a, aw, i, err)
  !     [parameters]
  !         lenaw     : length of aw (integer)
  !         tiny      : minimum value that 1/tiny does not 
  !                     overflow (real*8)
  !         eps       : relative error requested (real*8)
  !         aw        : points and weights of the quadrature 
  !                     formula, aw(0...lenaw-1) (real*8)
  !         f         : integrand f(x) (real*8 function)
  !         a         : lower limit of integration (real*8)
  !         i         : approximation to the integral (real*8)
  !         err       : estimate of the absolute error (real*8)
  !     [remarks]
  !         initial parameters
  !             lenaw > 1000, 
  !             IEEE double :
  !                 lenaw = 8000
  !                 tiny = 1.0d-307
  !         function
  !             f(x) needs to be analytic over (a,infinity).
  !         relative error
  !             eps is relative error requested excluding 
  !             cancellation of significant digits.
  !             i.e. eps means : (absolute error) / 
  !                              (integral_a^infinity |f(x)| dx).
  !             eps does not mean : (absolute error) / I.
  !         error message
  !             err >= 0 : normal termination.
  !             err < 0  : abnormal termination.
  !                        i.e. convergent error is detected :
  !                            1. f(x) or (d/dx)^n f(x) has 
  !                               discontinuous points or sharp 
  !                               peaks over (a,infinity).
  !                               you must divide the interval 
  !                               (a,infinity) at this points.
  !                            2. relative error of f(x) is 
  !                               greater than eps.
  !                            3. f(x) has oscillatory factor 
  !                               and decay of f(x) is very slow 
  !                               as x -> infinity.
  implicit none
  complex :: f
  real :: a, aw(0 : *)
  complex :: i, ir, fp, fm, iback, irback
  real :: err
  integer :: n
  complex :: h1
  real :: h2
  integer :: h3
  integer :: case

  integer :: noff, lenawm, nk, k, j, jtmp, jm, m, klim
  real :: epsh, errt, errh, errd, h

  noff = 5
  lenawm = int(aw(0) + 0.5)
  nk = int(aw(1) + 0.5)
  epsh = aw(4)
  i = f(a + aw(noff),n,h1,h2,h3,case)
  ir = i * aw(noff + 1)
  i = i * aw(noff + 2)
  err = abs(i)
  k = nk + noff
  j = noff
10 continue
  j = j + 6
  fm = f(a + aw(j),n,h1,h2,h3,case)
  fp = f(a + aw(j + 1),n,h1,h2,h3,case)
  ir = ir + (fm * aw(j + 2) + fp * aw(j + 3))
  fm = fm * aw(j + 4)
  fp = fp * aw(j + 5)
  i = i + (fm + fp)
  err = err + (abs(fm) + abs(fp))
  if (aw(j) .gt. epsh .and. j .lt. k) goto 10
  errt = err * aw(3)
  errh = err * epsh
  errd = 1 + 2 * errh
  jtmp = j
  do while (abs(fm) .gt. errt .and. j .lt. k)
     j = j + 6
     fm = f(a + aw(j),n,h1,h2,h3,case)
     ir = ir + fm * aw(j + 2)
     fm = fm * aw(j + 4)
     i = i + fm
  end do
  jm = j
  j = jtmp
  do while (abs(fp) .gt. errt .and. j .lt. k)
     j = j + 6
     fp = f(a + aw(j + 1),n,h1,h2,h3,case)
     ir = ir + fp * aw(j + 3)
     fp = fp * aw(j + 5)
     i = i + fp
  end do
  if (j .lt. jm) jm = j
  jm = jm - (noff + 6)
  h = 1
  m = 1
  klim = k + nk
  do while (errd .gt. errh .and. klim .le. lenawm)
     iback = i
     irback = ir
20   continue
     jtmp = k + jm
     do j = k + 6, jtmp, 6
        fm = f(a + aw(j),n,h1,h2,h3,case)
        fp = f(a + aw(j + 1),n,h1,h2,h3,case)
        ir = ir + (fm * aw(j + 2) + fp * aw(j + 3))
        i = i + (fm * aw(j + 4) + fp * aw(j + 5))
     end do
     k = k + nk
     j = jtmp
30   continue
     j = j + 6
     fm = f(a + aw(j),n,h1,h2,h3,case)
     ir = ir + fm * aw(j + 2)
     fm = fm * aw(j + 4)
     i = i + fm
     if (abs(fm) .gt. errt .and. j .lt. k) goto 30
     j = jtmp
40   continue
     j = j + 6
     fp = f(a + aw(j + 1),n,h1,h2,h3,case)
     ir = ir + fp * aw(j + 3)
     fp = fp * aw(j + 5)
     i = i + fp
     if (abs(fp) .gt. errt .and. j .lt. k) goto 40
     if (k .lt. klim) goto 20
     errd = h * (abs(i - 2 * iback) + abs(ir - 2 * irback))
     h = h * 0.5
     m = m * 2
     klim = 2 * klim - noff
  end do
  i = i * h
  if (errd .gt. errh) then
     err = -errd * m
  else
     err = err * (aw(2) * m)
  end if
end subroutine intdei

!> Initializes the double exponential quadrature method used for numerical integration. Written by Takuya Ooura, Reasearch Institute for Mathematical Science, Kyoto University (See http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html)
!! \param lenaw length of array aw
!! \param tiny minimum value that 1/tiny does not overflow
!! \param eps relative error requested
!! \param aw points and weights of the quadrature formula
subroutine intdeiini(lenaw, tiny, eps, aw)
  implicit none
  integer :: lenaw
  real :: tiny, eps, aw(0 : lenaw - 1)
  real :: efs, hoff
  integer :: noff, nk, k, j
  real :: pi4, tinyln, epsln, h0, ehp, ehm
  real :: h, t, ep, em, xp, xm, wp, wm

  ! ---- adjustable parameter ----
  efs = 0.1
  hoff = 11.0
  ! ------------------------------
  pi4 = atan(1.0)
  tinyln = -log(tiny)
  epsln = 1 - log(efs * eps)
  h0 = hoff / epsln
  ehp = exp(h0)
  ehm = 1 / ehp
  aw(2) = eps
  aw(3) = exp(-ehm * epsln)
  aw(4) = sqrt(efs * eps)
  noff = 5
  aw(noff) = 1
  aw(noff + 1) = 4 * h0
  aw(noff + 2) = 2 * pi4 * h0
  h = 2
  nk = 0
  k = noff + 6
10 continue
  t = h * 0.5
20 continue
  em = exp(h0 * t)
  ep = pi4 * em
  em = pi4 / em
  j = k
30 continue
  xp = exp(ep - em)
  xm = 1 / xp
  wp = xp * ((ep + em) * h0)
  wm = xm * ((ep + em) * h0)
  aw(j) = xm
  aw(j + 1) = xp
  aw(j + 2) = xm * (4 * h0)
  aw(j + 3) = xp * (4 * h0)
  aw(j + 4) = wm
  aw(j + 5) = wp
  ep = ep * ehp
  em = em * ehm
  j = j + 6
  if (ep .lt. tinyln .and. j .le. lenaw - 6) goto 30
  t = t + h
  k = k + nk
  if (t .lt. 1) goto 20
  h = h * 0.5
  if (nk .eq. 0) then
     if (j .gt. lenaw - 12) j = j - 6
     nk = j - noff
     k = k + nk
     aw(1) = nk
  end if
  if (2 * k - noff - 6 .le. lenaw) goto 10
  aw(0) = k - 6
end subroutine intdeiini

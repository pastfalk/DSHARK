subroutine cerror( z, cer )
  use param_mod
  !*****************************************************************************80
  !
  !! CERROR computes the error function for a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, complex ( kind = 8 ) Z, the argument.
  !
  !    Output, complex ( kind = 8 ) CER, the function value.
  !
  implicit none

  real :: a0
  complex :: c0
  complex :: cer
  complex :: cl
  complex :: cr
  complex :: cs
  integer :: k2
  complex :: z
  complex :: z1

  a0 = abs ( z )
  c0 = exp ( - z * z )

  z1 = z

  if ( real (z) .lt. 0.0 ) then
     z1 = - z
  end if

  if ( a0 .le. 5.8 ) then    

     cs = z1
     cr = z1
     do k2 = 1, 120
        cr = cr * z1 * z1 / ( k2 + 0.5 )
        cs = cs + cr
        if ( abs ( cr / cs ) .lt. 10.0**(-15) ) then
           exit
        end if
     end do

     cer = 2.0 * c0 * cs / sqrt ( pi )

  else

     cl = 1.0 / z1              
     cr = cl
     do k2 = 1, 13
        cr = -cr * ( k2 - 0.5 ) / ( z1 * z1 )
        cl = cl + cr
        if ( abs ( cr / cl ) .lt. 10.0**(-15) ) then
           exit
        end if
     end do

     cer = 1.0 - c0 * cl / sqrt ( pi )

  end if

  if ( real ( z) .lt. 0.0 ) then
     cer = -cer
  end if

end subroutine cerror

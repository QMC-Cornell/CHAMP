      function gauss()
! Written by Cyrus Umrigar
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Generate a gaussian random number using Box-Mueller method.
! Rewrote inline function as a function so that there would be no ambiquity as to
! the order in which the 2 rannyu's are evaluated.
! Could generate 2 numbers for almost the same price, but for
! backward compatibility generate just 1.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use const_mod
      implicit real*8(a-h,o-z)

      gauss=dcos(2*pi*rannyu(0))
      gauss=gauss*sqrt(-2*dlog(rannyu(0)))

      return
      end

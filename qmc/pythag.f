      FUNCTION pythag1(a,b)
!---------------------------------------------------
!     pythag: version from Numerical Recipes
!---------------------------------------------------
      DOUBLE PRECISION a,b,pythag1
      DOUBLE PRECISION absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag1=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag1=0.
        else
          pythag1=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END

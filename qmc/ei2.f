      FUNCTION ei2(x)
!-------------------------------------------------------------------------
!     General Ei(x) valid for x>0 and x<0
!-------------------------------------------------------------------------
      REAL*8 ei2,x

!     for x>0, use ei function of numerical recipes
      if(x > 0.d0) then
       ei2 = ei(x)

!     for x<0, use the relation: Ei(x) = -E1(-x) for x<0
      elseif (x < 0.d0) then
       ei2 = -expint(1,-x)

!     for x=0, ei(x) --> -infinity
      else
       ei2 = -1.d50
      endif

      return
      END

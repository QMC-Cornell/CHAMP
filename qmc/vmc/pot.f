      subroutine pot_nn(cent,znuc,iwctype,ncent,pecent)
c Written by Cyrus Umrigar
c get nuclear potential energy
      use control_mod
      use dim_mod
      use pseudo_mod
      use contrl_per_mod
      implicit real*8(a-h,o-z)

      dimension znuc(*),cent(3,ncent),iwctype(ncent)

      if(nloc.ge.0) then
        if(iperiodic.eq.0) then
          pecent=0
          do 20 i=2,ncent
            j1=i-1
            do 20 j=1,j1
              r2=0
              do 10 k=1,ndim
   10           r2=r2+(cent(k,i)-cent(k,j))**2
              r=dsqrt(r2)
   20         pecent=pecent+znuc(iwctype(i))*znuc(iwctype(j))/r
         else
          call pot_nn_ewald(cent,znuc,iwctype,ncent,pecent)
        endif
      endif

      return
      end

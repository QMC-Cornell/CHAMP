      subroutine pot_nn(cent,znuc,iwctype,ncent,pecent)
! Written by Cyrus Umrigar
! get nuclear potential energy
      use control_mod
      use dim_mod
      use pseudo_mod
      use contrl_per_mod
      use fragments_mod
      implicit real*8(a-h,o-z)

      dimension znuc(*),cent(3,ncent),iwctype(ncent)

      if(nloc.ge.0) then
        if(iperiodic.eq.0) then
          pecent=0
          if (l_fragments) pecent_frag=0
          do 20 i=2,ncent
            j1=i-1
            do 20 j=1,j1
              r2=0
              do 10 k=1,ndim
   10           r2=r2+(cent(k,i)-cent(k,j))**2
              r=dsqrt(r2)
              pot=znuc(iwctype(i))*znuc(iwctype(j))/r
              if (l_fragments) then
!                pecent_frag(iwfragnucl(i))=pecent_frag(iwfragnucl(i))+pot/2
!                pecent_frag(iwfragnucl(j))=pecent_frag(iwfragnucl(j))+pot/2
                iwfrag=merge(iwfragnucl(i),nfrag+1,iwfragnucl(j).EQ.iwfragnucl(i))
                pecent_frag(iwfrag)=pecent_frag(iwfrag)+pot
              endif
   20         pecent=pecent+pot
!   20         pecent=pecent+znuc(iwctype(i))*znuc(iwctype(j))/r
         else
          call pot_nn_ewald(cent,znuc,iwctype,ncent,pecent)
        endif
      endif

      return
      end

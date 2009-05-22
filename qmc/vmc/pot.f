      subroutine pot_nn(cent,znuc,iwctype,ncent,pecent)
      use control_mod
c Written by Cyrus Umrigar
c get nuclear potential energy
      use dim_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'pseudo.h'
!JT      include 'force.h'

!JT      common /dim/ ndim
      common /contrl_per/ iperiodic,ibasis
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
      dimension znuc(MCTYPE),cent(3,MCENT),iwctype(MCENT)

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

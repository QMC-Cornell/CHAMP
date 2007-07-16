      function zeta_cusp(zeta_old,znuc,term) ! author: Cyrus Umrigar
c Written by Cyrus Umrigar
c Solve quartic cusp-condition equation iteratively for zeta
      implicit real*8(a-h,o-z)
      parameter(cmix=0.5d0,cmix_min=1d-6,cmix_max=1d0)

      do 10 iter=1,100
        zeta_cusp=znuc+term/zeta_old**3
        diff=zeta_cusp-zeta_old
        if(iter.eq.1 .or. diff_old.eq.diff) then
          cmixx=cmix
         else
          cmixx=min(cmix_max,max(cmix_min,diff_i/(diff_old-diff)))
        endif
c       write(6,*) zeta_old,zeta_cusp,zeta_old+cmixx*diff,znuc
        zeta_old=min(zeta_old+cmixx*diff,2*zeta_old)
        if(abs(zeta_old-zeta_cusp).lt.1.d-15*zeta_old) goto 20
        diff_old=diff
   10   diff_i=cmixx*diff
      stop 'zeta_cusp not converged'
   20 zeta_cusp=zeta_old
c     write(6,'(''iter,zeta='',i3,f15.9)') iter,zeta_old
      return
      end

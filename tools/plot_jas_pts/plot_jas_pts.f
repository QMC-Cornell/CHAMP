      program plot_jas_pts
c Written by Cyrus Umrigar
c Generate 2-electron mc_configs suitable for plotting Jastrow
      implicit real*8(a-h,o-z)
      parameter(ntheta=11)

      pi=4*datan(1.d0)
      read(5,*) nr
      do 30 ir=1,nr
        read(5,*) r
        x1=r
        do 30 i=1,2
          r2=i*r
          do 10 itheta=1,ntheta
            x2=r2*cos((itheta-1)*pi/(ntheta-1))
            y2=sqrt(r2*r2-x2*x2)
            if(abs(x2).lt.1.d-6) x2=1.d-6
            if(abs(y2).lt.1.d-6) y2=1.d-6
   10       write(6,'(f5.2,'' 1.d-9 1.d-9'',x,2f10.6,'' 1.d-8 1 1. 1.'')') x1,x2,y2
   30   write(6,*)

      stop
      end

      program psp_needs_tm
c Convert h psp from Needs format to TM format
c The grid is r_i=r0*(exp(h*(i-1))-1)
c If we only know r_i but not r0 or h we can get r0 and h from
c r_3/r_2 = (exp(2a)-1)/(exp(h)-1)
c Defining x= exp(h), we solve quadratic equation for x to get
c x=(r_3-r_2)/r_2 or x=1.  We want the former root.
c We also have the option of writing only every n_th pt, starting
c from the first

      implicit real*8(h-h,o-z)
      parameter(MR=4000)
      character*8 units
      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 ititle(7),iray(6)

      dimension r(MR),psp(MR)

      pi=4*datan(1.d0)

      n_th=1

      read(5,*)
      read(5,*)
      read(5,*) z,zion
      read(5,*)
      read(5,'(a8)') units
      if(units.ne.' rydberg') stop 'expecting rydberg units'
      read(5,*)
      read(5,*)
      read(5,*)
      read(5,*)
      read(5,*)
      read(5,*) nr
      if(nr.gt.MR) stop 'nr > MR'
      read(5,*)
      read(5,*) (r(i),i=1,nr)
      x=(r(3)-r(2))/r(2)
      h=dlog(x)
      r0=r(2)/(x-1)
      h=n_th*h
      nrr=(nr-1)/n_th+1

      write(6,'(''r0,h='',1pd10.4,0p9f20.14)') r0,h

      open(1,file='pseudo.dat.1',form='unformatted')
      open(2,file='pseudopot1',form='formatted')

c Set comment variables to null for the moment
      icorr='ca'
      nameat='Si'
      irel='nrl'
      nicore='nc  '
      do 4 i=1,7
    4   ititle(i)=''
      do 6 i=1,6
    6   iray(i)=''

c Set number of psp channels to 3 for the moment
      npotd=3

c Warning: tmp
      write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     & (ititle(i),i=1,7),npotd,npotu,nrr-1,r0,h,zion
      write(2,'(a2,x,a2,x,a3,x,a4,x,13(a10,x))')
     & nameat,icorr,irel,nicore,(iray(i),i=1,6),(ititle(i),i=1,7)
      write(2,'(2i2,i5,1p2d22.15,0pf7.3,'' npotd,npotu,nr-1,r0,h,zion'')
     &') npotd,npotu,nrr-1,r0,h,zion
c     write(1) (r(i),i=2,nr)
c     write(2,'(1p10d22.15,0p)') (r(i),i=2,nr)
      write(1) (r(i),i=1+n_th,nr,n_th)
      write(2,'(1p10d22.15,0p)') (r(i),i=1+n_th,nr,n_th)

      do 10 l=1,10
        read(5,*,end=999)
        read(5,*,end=999) (psp(i),i=1,nr)
c       write(1) l-1,(psp(i),i=2,nr)
c  10   write(2,'(i2,(1p10d23.15,0p))') l-1,(psp(i),i=2,nr)
        write(1) l-1,(psp(i),i=1+n_th,nr,n_th)
   10   write(2,'(i2,(1p10d23.15,0p))') l-1,(psp(i),i=1+n_th,nr,n_th)

c Write garbage to satisfy read of cdc, cdd in psp. code, routine vionic.f
c Factor of 1/2 in cdd to get down-spin density, though Martins' program does
c not seem to care.
  999 write(1) (0.d0,ir=1+n_th,nr,n_th)
      r_cut=2.25d0
      alpha=1/r_cut**2
      anorm=0.5d0*zion*(alpha/pi)**1.5d0
c     write(1) (anorm*r(ir)*exp(-alpha*r(ir)**2),ir=2,nr)
      write(1) (anorm*r(ir)*exp(-alpha*r(ir)**2),ir=1+n_th,nr,n_th)
      write(2,'(1p10d23.15,0p)') (anorm*r(ir)*exp(-alpha*r(ir)**2),ir=1+n_th,nr,n_th)

      write(6,'(''2nd, last radial grid pt'',2d12.4)') r(2),r(nr)
c     write(6,*) (r(ir),ir=1,nr)
c     write(6,*) 'fake density*r'
c     write(6,*) (anorm*r(ir)*exp(-alpha*r(ir)**2),ir=1,nr)

      stop
      end

      subroutine finwrt_mpi
c MPI version created by Claudia Filippi starting from serial version
c routine to print out final results

# if defined (MPI)

      use mpi_mod

      implicit real*8(a-h,o-z)
c     character*16 mode
c     character*80 title
c     character*24 date


c     parameter (one=1.d0,half=.5d0)

      common /estsum/ esum1,esum(MFORCE),pesum,peisum,tpbsum,tjfsum,r2sum,accsum
      common /estcum/ ecum1,ecum(MFORCE),pecum,peicum,tpbcum,tjfcum,r2cum,acccum,iblk
      common /est2cm/ ecm21,ecm2(MFORCE),pecm2,peicm2,tpbcm2,tjfcm2,r2cm2
      common /estsig/ wsum1s(MFORCE),esum1s(MFORCE),ecum1s(MFORCE),ecm21s(MFORCE)
      common /stepv/try(NRAD),suc(NRAD),trunfb(NRAD),rprob(NRAD), !JT
     &ekin(NRAD),ekin2(NRAD)
      common /forcepar/ deltot(MFORCE),nforce,istrech
c     common /forcest/ fcum(MFORCE),fcm2(MFORCE)
c     common /forcewt/ wsum(MFORCE),wcum(MFORCE)
c     common /forcjac/ ajacob
      common /denupdn/ rprobup(NRAD),rprobdn(NRAD)
      common /pairden/ xx0probut(0:NAX,-NAX:NAX,-NAX:NAX),xx0probuu(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probud(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdt(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probdu(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdd(0:NAX,-NAX:NAX,-NAX:NAX),
     &den2d_t(-NAX:NAX,-NAX:NAX),den2d_d(-NAX:NAX,-NAX:NAX),den2d_u(-NAX:NAX,-NAX:NAX),
     &delxi,xmax,xfix(3),ifixe

c     dimension trunfbt(NRAD),rprobt(NRAD),ekint(NRAD),ekin2t(NRAD)
      dimension xx0probt(0:NAX,-NAX:NAX,-NAX:NAX),den2dt(-NAX:NAX,-NAX:NAX)
      dimension rprobt(NRAD),tryt(NRAD),suct(NRAD)
     &,work(MFORCE)


c     err(x,x2,i)=dsqrt(abs(x2/wcum(i)-(x/wcum(i))**2)/iblk)
c     err1(x,x2)=dsqrt(dabs(x2/passes-(x/passes)**2)/passes)

      call mpi_allreduce(ecum1,ecum1t,1,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ecm21,ecm21t,1,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      ecum1=ecum1t
      ecm21=ecm21t

      call mpi_allreduce(wsum1s,work,nforce,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 2 ifr=1,nforce
    2   wsum1s(ifr)=work(ifr)

      call mpi_allreduce(esum1s,work,nforce,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 4 ifr=1,nforce
    4   esum1s(ifr)=work(ifr)

      call mpi_allreduce(ecum1s,work,nforce,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 6 ifr=1,nforce
    6   ecum1s(ifr)=work(ifr)

      call mpi_allreduce(ecm21s,work,nforce,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 8 ifr=1,nforce
    8   ecm21s(ifr)=work(ifr)

      write(6,*) 'rprob'
      call mpi_allreduce(rprob,rprobt,NRAD,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 10 i=1,NRAD
   10   rprob(i)=rprobt(i)
      call mpi_allreduce(rprobup,rprobt,NRAD,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 12 i=1,NRAD
   12   rprobup(i)=rprobt(i)
      call mpi_allreduce(rprobdn,rprobt,NRAD,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 14 i=1,NRAD
   14   rprobdn(i)=rprobt(i)

      call mpi_allreduce(suc,suct,NRAD,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 20 i=1,NRAD
   20   suc(i)=suct(i)
      call mpi_allreduce(try,tryt,NRAD,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      do 30 i=1,NRAD
   30   try(i)=tryt(i)

      if(ifixe.ne.0) then

        naxt=(NAX+1)*(2*NAX+1)*(2*NAX+1)
        call mpi_allreduce(xx0probut,xx0probt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 40 i1=0,NAX
          do 40 i2=-NAX,NAX
            do 40 i3=-NAX,NAX
   40         xx0probut(i1,i2,i3)=xx0probt(i1,i2,i3)
        call mpi_allreduce(xx0probud,xx0probt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 41 i1=0,NAX
          do 41 i2=-NAX,NAX
            do 41 i3=-NAX,NAX
   41         xx0probud(i1,i2,i3)=xx0probt(i1,i2,i3)
        call mpi_allreduce(xx0probuu,xx0probt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 42 i1=0,NAX
          do 42 i2=-NAX,NAX
            do 42 i3=-NAX,NAX
   42         xx0probuu(i1,i2,i3)=xx0probt(i1,i2,i3)
        call mpi_allreduce(xx0probdu,xx0probt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 43 i1=0,NAX
          do 43 i2=-NAX,NAX
            do 43 i3=-NAX,NAX
   43         xx0probdu(i1,i2,i3)=xx0probt(i1,i2,i3)
        call mpi_allreduce(xx0probdd,xx0probt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 44 i1=0,NAX
          do 44 i2=-NAX,NAX
            do 44 i3=-NAX,NAX
   44         xx0probdd(i1,i2,i3)=xx0probt(i1,i2,i3)
        call mpi_allreduce(xx0probdt,xx0probt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 45 i1=0,NAX
          do 45 i2=-NAX,NAX
            do 45 i3=-NAX,NAX
   45         xx0probdt(i1,i2,i3)=xx0probt(i1,i2,i3)
        naxt=(2*NAX+1)*(2*NAX+1)
        call mpi_allreduce(den2d_t,den2dt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 50 i1=-NAX,NAX
          do 50 i2=-NAX,NAX
   50       den2d_t(i1,i2)=den2dt(i1,i2)
        call mpi_allreduce(den2d_d,den2dt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 51 i1=-NAX,NAX
          do 51 i2=-NAX,NAX
   51       den2d_d(i1,i2)=den2dt(i1,i2)
        call mpi_allreduce(den2d_u,den2dt,naxt,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
        do 52 i1=-NAX,NAX
          do 52 i2=-NAX,NAX
   52       den2d_u(i1,i2)=den2dt(i1,i2)

      endif

c     call mpi_allreduce(trunfb,trunfbt,NRAD,mpi_double_precision
c    &,mpi_sum,MPI_COMM_WORLD,ierr)
c     do 40 i=1,NRAD
c  40   trunfb(i)=trunfbt(i)
c     call mpi_allreduce(ekin,ekint,NRAD,mpi_double_precision
c    &,mpi_sum,MPI_COMM_WORLD,ierr)
c     do 50 i=1,NRAD
c  50   ekin(i)=ekint(i)
c     call mpi_allreduce(ekin2,ekin2t,NRAD,mpi_double_precision
c    &,mpi_sum,MPI_COMM_WORLD,ierr)
c     do 60 i=1,NRAD
c  60   ekin2(i)=ekin2t(i)

      call mpi_allreduce(accsum,accept,1,mpi_double_precision
     &,mpi_sum,MPI_COMM_WORLD,ierr)
      accsum=accept

c call to grad_hess_jas_mpi moved to ../finwrt.f for the moment
c     write(6,*) 'before grad_hess_jas_mpi'
c     call grad_hess_jas_mpi
c     write(6,*) 'after grad_hess_jas_mpi'

# endif
      return
      end

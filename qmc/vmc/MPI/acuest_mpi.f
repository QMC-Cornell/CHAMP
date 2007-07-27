      subroutine acuest_mpi
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

# if defined (MPI)

      use mpi_mod

      implicit real*8(a-h,o-z)
      character*16 mode

!JT      parameter (half=.5d0)

      common /dim/ ndim
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fcum(MFORCE),fcm2(MFORCE)
      common /forcewt/ wsum(MFORCE),wcum(MFORCE)
      common /forcjac/ ajacob

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /contr3/ mode
      common /config/ xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)
     &,vnew(3,MELEC),psi2o(MFORCE),psi2n(MFORCE),eold(MFORCE),enew(MFORCE)
     &,peo,pen,tjfn,tjfo,psido,psijo
     &,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)
     &,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)
     &,nearesto(MELEC),nearestn(MELEC),delttn(MELEC)
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad
      common /estsum/ esum1,eaverage(MFORCE+6)
      common /estcum/ ecum1,ecum(MFORCE),pecum,peicum,tpbcum,tjfcum,r2cum,acccum,iblk
      common /est2cm/ ecm21,ecm2(MFORCE),pecm2,peicm2,tpbcm2,tjfcm2,r2cm2
      common /estsig/ wsum1s(MFORCE),esum1s(MFORCE),ecum1s(MFORCE),ecm21s(MFORCE)
      common /stepv/try(NRAD),suc(NRAD),trunfb(NRAD),rprob(NRAD), !JT
     &ekin(NRAD),ekin2(NRAD)
      common /denupdn/ rprobup(NRAD),rprobdn(NRAD)
      common /doefp/ nefp
      common /div_v/ div_vo(MELEC)
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)


      dimension xstrech(3,MELEC),ecollect(MFORCE+6),wcollect(MFORCE)

c statement function for error calculation
c     err(x,x2)=dsqrt(dabs(x2/iblk-(x/iblk)**2)/iblk)
      err(x,x2,i)=dsqrt(abs(x2/wcum(i)-(x/wcum(i))**2)/iblk)

c xsum = sum of values of x from metrop
c xnow = average of values of x from metrop
c xcum = accumulated sums of xnow
c xcm2 = accumulated sums of xnow**2
c xave = current average value of x
c xerr = current error of x

c eaverage(1-MFORCE) = esum(1-MFORCE)
c eaverage(MFORCE+1) = pesum
c eaverage(MFORCE+2) = peisum
c eaverage(MFORCE+3) = tpbsum
c eaverage(MFORCE+4) = tjfsum
c eaverage(MFORCE+5) = r2sum
c eaverage(MFORCE+6) = accsum


c Note we do not reduce the 1-electron move quantities here because
c they are not printed out from acuest.  So we just cumulate the
c quantities on individual processors, and reduce the cumulated
c quantities in finwrt_mpi

      call mpi_allreduce(eaverage,ecollect,MFORCE+6
     &,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      call mpi_allreduce(wsum,wcollect,nforce
     &,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

c Warning this flush and barrier should not be necessary
      call systemflush(6)
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      do 10 ifr=1,MFORCE+6
        eaverage(ifr)=ecollect(ifr)
   10   if(ifr.le.nforce) wsum(ifr)=wcollect(ifr)
c     pesum=ecollect(MFORCE+1)
c     peisum=ecollect(MFORCE+2)
c     tpbsum=ecollect(MFORCE+3)
c     tjfsum=ecollect(MFORCE+4)
c     r2sum=ecollect(MFORCE+5)
c     accsum=ecollect(MFORCE+6)

# endif

      return
      end

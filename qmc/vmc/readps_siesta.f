      subroutine readps_siesta
c Written by W. Al-Saidi 

c read Troullier-Martins pseudopotentials
c reads in r*v in ryd.
c does 3 conversions: a) ryd -> Har, b) r*v -> v and
c c) subtracts out local part from all except highest l component.
c Also eval pot. at 0 and initializes quadrature pts.
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'

      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 ititle(7),iray(6)
      character*20 filename,atomtyp
      character*10 ntitle(20)

      parameter (ncoef=5)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /pseudo_tm/ rmax_coul(MCTYPE),rmax_nloc(MCTYPE),exp_h_ps(MCTYPE),r0_ps(MCTYPE)
     &,vpseudo(MPS_GRID,MCTYPE,MPS_L),d2pot(MPS_GRID,MCTYPE,MPS_L),igrid_ps(MCTYPE),nr_ps(MCTYPE)
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc

      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad

      dimension r(MPS_GRID),y(ncoef),ce(ncoef),dmatr(ncoef*ncoef)
      dimension work(MPS_GRID)

c Warning temp:
      dimension vpot(2)

      do 200 ict=1,nctype

        igrid_ps(ict)=3

        if(ict.lt.10) then
          write(atomtyp,'(i1)') ict
         elseif(ict.lt.100) then
          write(atomtyp,'(i2)') ict
         else
          stop 'readps_tm, nctype>100'
        endif

        filename='pseudopot'//atomtyp(1:index(atomtyp,' ')-1)
        open(1,file=filename,status='old',form='formatted',err=999)
        write(6,'(''Reading pseudopotential file '',a20)') filename

c The TM psp. format has npotd and npotu for down and up, but we just use one of them
c They are the number of different l components of the psp.


        read(1,*) (ititle(i), i=1,4)
        write(6,*) (ititle(i), i=1,4)
        read(1,*) (ititle(i), i=1,7)
        write(6,*) (ititle(i), i=1,7)
        read(1,*) 
c        write(6,*) (ntitle(i), i=1,20)
        read(1,*) npotd(ict),npotu,nrm,r0_ps(ict),h_ps,zion
        write(6,'(''Center type'',i2,'' has'',i2,'' pseudopotential L components, and component''
     &  ,i2,'' is chosen to be local'')') ict,npotd(ict),lpotp1(ict)
        lpotp1(ict)=0
        if(lpotp1(ict).gt.npotd(ict)) then
           write(6,'(''lpotp1(ict),npotd(ict)='',2i3)') lpotp1(ict),npotd(ict)
           stop 'Cannot choose local psp. to be > number of l components, lpotp1(ict) > npotd(ict)'
        endif
        if(npotd(ict).gt.MPS_L) stop 'npotd(ict).gt.MPS_L'

c If the local pseudopot component is not set in input, set it here
        if(lpotp1(ict).le.0) then
          lpotp1(ict)=npotd(ict)
          write(6,'(''Center type'',i4,'' local pseudopot component is'',i3)') ict,lpotp1(ict)
        endif

        nr=nrm+1
        nr_ps(ict)=nr
        write(6,'(''ict,nr_ps(ict)'',9i5)') ict,nr_ps(ict)

        if(znuc(ict).ne.zion) then
          write(6,'(''znuc(ict) != zion in readps_tm'',2f6.1)') znuc(ict),zion
          stop 'znuc(ict) != zion in readps_tm'
        endif
        if(nr.gt.MPS_GRID) then
          write(6,'(''nr > MPS_GRID'',2i6)') nr,MPS_GRID
          stop 'nr > MPS_GRID'
        endif

        r(1)=0.d0
        read(1,*) 
        read(1,*) (r(i),i=2,nr)
c     Read in psp and divide by 2 to convert from Ryd. to H.
        do i=1,npotd(ict)
           read(1,*)
           read(1,*) ll
           read(1,*)(vpseudo(j,ict,i),j=2,nr)
           do ir=2,nr
              vpseudo(ir,ict,i)=0.5d0*vpseudo(ir,ict,i)/r(ir)
           enddo
        enddo


!!! ploting pot
        if(ipr.ge.0) then
           do ir=2, nr
              write (55+ict,'(5F14.6)') r(ir), (vpseudo(ir,ict,i), i=1, npotd(ict))
           enddo
        endif
!!!

!!!!
!!!!! fix the noise at the begining
!           do i=1,npotd(ict)
!              do ir=1, 5
c     vpseudo(ir,ict,i) = vpseudo(6,ict,i)
!              enddo
!           enddo
!!!!!!



c     Find the point beyond which r*v differs from zion by no more than .5*d-6.
c irmax_coul is used for the endpoint of the spline where it is assumed that
c the derivative of the local component is zion/r(irmax_coul)**2 and that of 
c the local component is 0.  Also, rmax_coul is used in splfit_ps when it is
c called in a calculation of a periodic system.
        rmax_coul(ict)=0.d0
        irmax_coul=0
        do 50 ir=nr,1,-1
          do 50 i=1,npotd(ict)
            if(dabs(r(ir)*vpseudo(ir,ict,i)+zion).gt..5d-6) then
              rmax_coul(ict)=max(rmax_coul(ict),r(ir))
              irmax_coul=ir
              goto 60
            endif
   50   continue
   60   irmax_coul=min(irmax_coul+5,nr)

c Find the point beyond which the various v components differ from each other by no more than .5*d-6
        rmax_nloc(ict)=0.d0
        irmax_nloc=0
        do 70 ir=nr,1,-1
          do 70 i=2,npotd(ict)
            do 70 j=1,i-1
              if(dabs(vpseudo(ir,ict,i)-vpseudo(ir,ict,j)).gt..5d-6) then
                rmax_nloc(ict)=max(rmax_nloc(ict),r(ir))
                irmax_nloc=ir
                goto 80
              endif
   70   continue
   80   irmax_nloc=irmax_nloc+1
        rmax_nloc(ict)=r(irmax_nloc)

c rmax_nloc is used in getvps_champ to decide whether to calculate calculate nonloc part of psp
c or to set it to zero.  irmax_coul is used to decide how far out to spline the psp. components
c so irmax_coul must be >= irmax_nloc.
        irmax_coul=max(irmax_coul,irmax_nloc)
        rmax_coul(ict)=r(irmax_coul)

        write(6,'(''center '',i3,'' pseudopot rmax_coul,irmax_coul,rmax_nloc,irmax_nloc= '',2(f6.2,i5))')
     &  ict,rmax_coul(ict),irmax_coul,rmax_nloc(ict),irmax_nloc

c       rmax_coul(ict)=0.d0
c       irmax_coul=0
c       do 100 i=1,npotd(ict)
c         if(nloc.eq.2) then
c           read(1) ii,(vpseudo(j,ict,i),j=2,nr)
c          elseif(nloc.eq.3) then
c           read(1,*) ii,(vpseudo(j,ict,i),j=2,nr)
c         endif
c         vpseudo(1,ict,i)=0.d0
c         do 50 j=nr,1,-1
cc          write(6,'(''vps'',f8.4,f12.6,d12.4)') r(j),vpseudo(j,ict,i)/2,vpseudo(j,ict,i)/2+zion
c           if(dabs(vpseudo(j,ict,i)+2.d0*zion).gt.1.d-6) then
c             rmax_coul(ict)=max(rmax_coul(ict),r(j))
c             irmax_coul=max(irmax_coul,j)
c             goto 60
c           endif
c  50     continue
c  60     do 100 j=2,nr
c           vpseudo(j,ict,i)=0.5d0*vpseudo(j,ict,i)/r(j)
c 100   continue
c       irmax_coul=irmax_coul+5
c       write(6,'(''center '',i3,'' pseudopot rmax,irmax_coul= '',f6.2,i5)')
c    &  ict,rmax_coul(ict),irmax_coul

        close(1)

c If the local pseudopot component is not set in input, set it here
c Already done before here
c       if(lpotp1(ict).le.0) then
c         lpotp1(ict)=npotd(ict)
c         write(6,'(''Center type'',i4,'' local pseudopot component is'',i3)') ict,lpotp1(ict)
c       endif

        exp_h_ps(ict)=dexp(h_ps)
        write(6,'(''potential grid param= '',2f10.4)') r0_ps(ict),exp_h_ps(ict)

        if(ipr.ge.5) then
          write(38,'(''r(j)  (vpseudo(j,ict,i),i=1,npotd(ict))  -znuc(ict)/r(j)'')')
          do 104 j=2,nr
            if(r(j).gt.0.d0) then
!             write(38,'(1pd12.6,9d14.6)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict)),-znuc(ict)/r(j)
              write(38,'(1pd12.6,0p,9g22.14)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict)),-znuc(ict)/r(j)
             else
!             write(38,'(1pd12.6,9d14.6)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict))
              write(38,'(1pd12.6,0p,9g22.14)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict))
            endif
  104     continue
        endif

        do 110 i=1,npotd(ict)
          if(i.ne.lpotp1(ict)) then
            do 105 j=2,nr
              if(ipr.ge.1) then
              endif
  105         vpseudo(j,ict,i)=vpseudo(j,ict,i)-vpseudo(j,ict,lpotp1(ict))
          endif

  110   continue

        if(rmax_coul(ict).eq.0.d0) goto 200
        do 190 i=1,npotd(ict)
c small radii pot(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
          ll=0
          do 120 jj=1,ncoef
            y(jj)=vpseudo(jj+1,ict,i)
            do 120 ii=1,ncoef
              ll=ll+1
  120         dmatr(ll)=r(ii+1)**(jj-1)
          call dgelg(y,dmatr,ncoef,1,1.d-8,ier)

          do 125 icoef=1,ncoef
  125     ce(icoef)=y(icoef)

c         write(6,'(''coefficients'',1p10e22.10)') (ce(iff),iff=1,ncoef)
          if(ipr.ge.1) then
            write(6,'(''check the small radius expansion for l= '',i3)')i-1
            write(6,'(''irad, rad, extrapolated, correct value, diff'')')
            do 130 ir=1,10
              val=ce(1)
              do 128 icoef=2,ncoef
  128           val=val+ce(icoef)*r(ir)**(icoef-1)
c               if(ir.eq.1) vpseudo(ir,ict,i)=val
  130           write(6,'(i2,1p3e16.8,1p1e11.2)')
     &          ir,r(ir),val,vpseudo(ir,ict,i),val-vpseudo(ir,ict,i)
          endif

          vpseudo(1,ict,i)=ce(1)

          dpot1=ce(2)
          do 135 icoef=3,ncoef
  135     dpot1=dpot1+(icoef-1)*ce(icoef)*r(1)**(icoef-2)

          dpotn=0.d0
          if(i.eq.lpotp1(ict)) dpotn=zion/r(irmax_coul)**2

          write(6,'(''dpot1,dpotn'',1p2e15.5)') dpot1,dpotn

c get second derivative for spline fit
          call spline2(r,vpseudo(1,ict,i),irmax_coul,dpot1,dpotn,d2pot(1,ict,i),work)

          do 190 j=1,nr
            if(ipr.ge.2) then
              if(i.eq.1) write(35,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
              if(i.eq.2) write(36,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
              if(i.eq.3) write(37,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
c             if(i.eq.1.and.i.ne.lpotp1(ict)
c    &        write(38,'(1p5d14.6)') r(j),vpseudo(j,ict,i)+vpseudo(j,ict,lpotp1(ict)),vpseudo(j,ict,lpotp1(ict)),-znuc(ict)/r(j)
c             if(i.eq.2.and.i.ne.lpotp1(ict))
c    &        write(39,'(1p5d14.6)') r(j),vpseudo(j,ict,i)+vpseudo(j,ict,lpotp1(ict)),vpseudo(j,ict,lpotp1(ict)),-znuc(ict)/r(j)
            endif
  190   continue

c Warning: Temporary print of psp
c       write(20,'(2f9.5,'' znuc,rpotc'')') znuc(1)
c       write(20,'(i6,f18.14,'' nr_ps(ict),exp_h_ps(ict)'')') nr,exp_h_ps(ict)
c       do 195 ir=1,nr
c 195     write(20,'(2d20.12)') r(ir),vpseudo(ir,ict,lpotp1(ict))

        
!!! ploting pot
        if(ipr.ge.0) then
           do ir=2, nr
              write (75+ict,'(5F14.6)') r(ir), (vpseudo(ir,ict,i), i=1, npotd(ict))
           enddo
        endif
!!!



 200  continue
      

c     call gesqua(nquad,xq0,yq0,zq0,wq)
c     call gesqua(nquad,xq,yq,zq,wq)

c      write(6,'(''quadrature points'')')
c      do 210 i=1,nquad
c  210   write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

c Temporary check of splining
c     write(6,'(''r vpot(ict=2) '')')
c     rr=1.d-6 
c     ict=2
c     do 230 i=1,4000
c       do 220 l=1,2
c 220     call splfit_ps(rr,l,ict,vpot(l))
c       write(6,'(9d22.14)') rr,(vpot(l),l=1,2)
c       if(rr.gt.5.d0) goto 240
c 230 rr=rr*1.01
c 240 write(6,'(''end of r vpot(ict=2) '')')

      return

  999 write(6,'(''Error: Pseudopot. file '',a20,'' is missing'')') filename
      stop 'Pseudopot. file is missing'

      end

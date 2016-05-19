      subroutine readps_tm
! Written by Claudia Filippi.  Modified by Cyrus Umrigar

! read Troullier-Martins pseudopotentials
! reads in r*v in ryd.
! does 3 conversions: a) ryd -> Har, b) r*v -> v and
! c) subtracts out local part from all except highest l component.
! Also eval pot. at r=0
      use all_tools_mod
      use atom_mod
      use const_mod
      use pseudo_mod
      use qua_mod
      use pseudo_tm_mod
      implicit real*8(a-h,o-z)

      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 ititle(7),iray(6)
      character*20 filename,atomtyp

      parameter (ncoef=5)

      dimension r(MPS_GRID),y(ncoef),ce(ncoef),dmatr(ncoef*ncoef)
      dimension work(MPS_GRID)

! Warning temp:
      dimension vpot(2)

      call alloc ('rmax_coul', rmax_coul, nctype)
      call alloc ('rmax_nloc', rmax_nloc, nctype)
      call alloc ('exp_h_ps', exp_h_ps, nctype)
      call alloc ('r0_ps', r0_ps, nctype)
      call alloc ('vpseudo', vpseudo, MPS_GRID, nctype, MPS_L)
      call alloc ('d2pot', d2pot, MPS_GRID, nctype, MPS_L)
      call alloc ('igrid_ps', igrid_ps, nctype)
      call alloc ('nr_ps', nr_ps, nctype)

      do 200 ict=1,nctype

        igrid_ps(ict)=3

        if(ict.lt.10) then
          write(atomtyp,'(i1)') ict
         elseif(ict.lt.100) then
          write(atomtyp,'(i2)') ict
         else
          stop 'readps_tm, nctype>100'
        endif

        if(nloc.eq.2) then
          filename='pseudo.dat.'//atomtyp(1:index(atomtyp,' ')-1)
          open(1,file=filename,status='old',form='unformatted',err=999)
         elseif(nloc.eq.3) then
          filename='pseudopot'//atomtyp(1:index(atomtyp,' ')-1)
          open(1,file=filename,status='old',form='formatted',err=999)
        endif
        write(6,'(''Reading pseudopotential file '',a20)') filename

! The TM psp. format has npotd and npotu for down and up, but we just use one of them
! They are the number of different l components of the psp.
        if(nloc.eq.2) then
          read(1) nameat,icorr,irel,nicore,(iray(i),i=1,6), &
     &    (ititle(i),i=1,7),npotd(ict),npotu,nrm,r0_ps(ict),h_ps,zion
         elseif(nloc.eq.3) then
          read(1,'(a2,x,a2,x,a3,x,a4,x,13(a10,x))') nameat,icorr,irel,nicore &
     &    ,(iray(i),i=1,6),(ititle(i),i=1,7)
          read(1,*) npotd(ict),npotu,nrm,r0_ps(ict),h_ps,zion
        endif
        write(6,'(''Center type'',i2,'' has'',i2,'' pseudopotential L components, and component'' &
     &  ,i2,'' is chosen to be local'')') ict,npotd(ict),lpotp1(ict)

        if(lpotp1(ict).gt.npotd(ict)) then
          write(6,'(''lpotp1(ict),npotd(ict)='',2i3)') lpotp1(ict),npotd(ict)
          stop 'Cannot choose local psp. to be > number of l components, lpotp1(ict) > npotd(ict)'
        endif
        if(npotd(ict).gt.MPS_L) stop 'npotd(ict).gt.MPS_L'

! If the local pseudopot component is not set in input, set it here
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
        if(nloc.eq.2) then
          read(1) (r(i),i=2,nr)
         elseif(nloc.eq.3) then
          read(1,*) (r(i),i=2,nr)
        endif

! Read in psp and divide by 2 to convert from Ryd. to H.
        do 40 i=1,npotd(ict)
          if(nloc.eq.2) then
            read(1) ii,(vpseudo(j,ict,i),j=2,nr)
           elseif(nloc.eq.3) then
            read(1,*) ii,(vpseudo(j,ict,i),j=2,nr)
          endif
          do 40 ir=2,nr
   40       vpseudo(ir,ict,i)=0.5d0*vpseudo(ir,ict,i)/r(ir)

! Find the point beyond which r*v differs from zion by no more than .5*d-6.
! irmax_coul is used for the endpoint of the spline where it is assumed that
! the derivative of the local component is zion/r(irmax_coul)**2 and that of
! the local component is 0.  Also, rmax_coul is used in splfit_ps when it is
! called in a calculation of a periodic system.
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

! Find the point beyond which the various v components differ from each other by no more than .5*d-6
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

! rmax_nloc is used in getvps_champ to decide whether to calculate calculate nonloc part of psp
! or to set it to zero.  irmax_coul is used to decide how far out to spline the psp. components
! so irmax_coul must be >= irmax_nloc.
        irmax_coul=max(irmax_coul,irmax_nloc)
        rmax_coul(ict)=r(irmax_coul)

        write(6,'(''center '',i3,'' pseudopot rmax_coul,irmax_coul,rmax_nloc,irmax_nloc= '',2(f6.2,i5))') &
     &  ict,rmax_coul(ict),irmax_coul,rmax_nloc(ict),irmax_nloc

!       rmax_coul(ict)=0.d0
!       irmax_coul=0
!       do 100 i=1,npotd(ict)
!         if(nloc.eq.2) then
!           read(1) ii,(vpseudo(j,ict,i),j=2,nr)
!          elseif(nloc.eq.3) then
!           read(1,*) ii,(vpseudo(j,ict,i),j=2,nr)
!         endif
!         vpseudo(1,ict,i)=0.d0
!         do 50 j=nr,1,-1
!c          write(6,'(''vps'',f8.4,f12.6,d12.4)') r(j),vpseudo(j,ict,i)/2,vpseudo(j,ict,i)/2+zion
!           if(dabs(vpseudo(j,ict,i)+2.d0*zion).gt.1.d-6) then
!             rmax_coul(ict)=max(rmax_coul(ict),r(j))
!             irmax_coul=max(irmax_coul,j)
!             goto 60
!           endif
!  50     continue
!  60     do 100 j=2,nr
!           vpseudo(j,ict,i)=0.5d0*vpseudo(j,ict,i)/r(j)
! 100   continue
!       irmax_coul=irmax_coul+5
!       write(6,'(''center '',i3,'' pseudopot rmax,irmax_coul= '',f6.2,i5)')
!    &  ict,rmax_coul(ict),irmax_coul

        close(1)

! If the local pseudopot component is not set in input, set it here
! Already done before here
!       if(lpotp1(ict).le.0) then
!         lpotp1(ict)=npotd(ict)
!         write(6,'(''Center type'',i4,'' local pseudopot component is'',i3)') ict,lpotp1(ict)
!       endif

        exp_h_ps(ict)=dexp(h_ps)
        write(6,'(''potential grid param= '',2f10.4)') r0_ps(ict),exp_h_ps(ict)

        if(ipr.ge.0) then
         write(38,'(''r(j)  (vpseudo(j,ict,i),i=1,npotd(ict))  -znuc(ict)/r(j)'')')
          do 104 j=2,nr
            if(r(j).gt.0.d0) then
!             write(38,'(1pd13.6,9d14.6)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict)),-znuc(ict)/r(j)
             write(38,'(1pd13.6,0p,9g22.14)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict)),-znuc(ict)/r(j)
             else
!             write(38,'(1pd13.6,9d14.6)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict))
             write(38,'(1pd13.6,0p,9g22.14)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict))
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
! small radii pot(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
          ll=0
          do 120 jj=1,ncoef
            y(jj)=vpseudo(jj+1,ict,i)
            do 120 ii=1,ncoef
              ll=ll+1
  120         dmatr(ll)=r(ii+1)**(jj-1)
          call dgelg(y,dmatr,ncoef,1,1.d-8,ier)

          do 125 icoef=1,ncoef
  125     ce(icoef)=y(icoef)

!         write(6,'(''coefficients'',1p10e22.10)') (ce(iff),iff=1,ncoef)
          if(ipr.ge.1) then
            write(6,'(''check the small radius expansion for l= '',i3)')i-1
            write(6,'(''irad, rad, extrapolated, correct value, diff'')')
            do 130 ir=1,10
              val=ce(1)
              do 128 icoef=2,ncoef
  128           val=val+ce(icoef)*r(ir)**(icoef-1)
!               if(ir.eq.1) vpseudo(ir,ict,i)=val
  130           write(6,'(i2,1p3e16.8,1p1e11.2)') &
     &          ir,r(ir),val,vpseudo(ir,ict,i),val-vpseudo(ir,ict,i)
          endif

          vpseudo(1,ict,i)=ce(1)

          dpot1=ce(2)
          do 135 icoef=3,ncoef
  135     dpot1=dpot1+(icoef-1)*ce(icoef)*r(1)**(icoef-2)

          dpotn=0.d0
          if(i.eq.lpotp1(ict)) dpotn=zion/r(irmax_coul)**2

          write(6,'(''dpot1,dpotn'',1p2e15.5)') dpot1,dpotn

! get second derivative for spline fit
          call spline2(r,vpseudo(1,ict,i),irmax_coul,dpot1,dpotn,d2pot(1,ict,i),work)

          do 190 j=1,nr
            if(ipr.ge.2) then
              if(i.eq.1) write(35,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
              if(i.eq.2) write(36,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
              if(i.eq.3) write(37,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
!             if(i.eq.1.and.i.ne.lpotp1(ict)
!    &        write(38,'(1p5d14.6)') r(j),vpseudo(j,ict,i)+vpseudo(j,ict,lpotp1(ict)),vpseudo(j,ict,lpotp1(ict)),-znuc(ict)/r(j)
!             if(i.eq.2.and.i.ne.lpotp1(ict))
!    &        write(39,'(1p5d14.6)') r(j),vpseudo(j,ict,i)+vpseudo(j,ict,lpotp1(ict)),vpseudo(j,ict,lpotp1(ict)),-znuc(ict)/r(j)
            endif
  190   continue

! Warning: Temporary print of psp
!       write(20,'(2f9.5,'' znuc,rpotc'')') znuc(1)
!       write(20,'(i6,f18.14,'' nr_ps(ict),exp_h_ps(ict)'')') nr,exp_h_ps(ict)
!       do 195 ir=1,nr
! 195     write(20,'(2d20.12)') r(ir),vpseudo(ir,ict,lpotp1(ict))

  200 continue

!     call gesqua(nquad,xq0,yq0,zq0,wq)
!     call gesqua(nquad,xq0,yq0,zq0,wq)
!     call gesqua(nquad,xq,yq,zq,wq)

!     write(6,'(''quadrature points'')')
!     do 210 i=1,nquad
! 210   write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

! Temporary check of splining
!     write(6,'(''r vpot(ict=2) '')')
!     rr=1.d-6
!     ict=2
!     do 230 i=1,4000
!       do 220 l=1,2
! 220     call splfit_ps(rr,l,ict,vpot(l))
!       write(6,'(9d22.14)') rr,(vpot(l),l=1,2)
!       if(rr.gt.5.d0) goto 240
! 230 rr=rr*1.01
! 240 write(6,'(''end of r vpot(ict=2) '')')

      return

  999 write(6,'(''Error: Pseudopot. file '',a20,'' is missing'')') filename
      stop 'Pseudopot. file is missing'

      end
!-----------------------------------------------------------------------

      subroutine getvps_champ(r_en,iel)
! compute pseudopotential for electron iel

      use atom_mod
      use const_mod
      use pseudo_mod
      use pseudo_tm_mod
      implicit real*8(a-h,o-z)

      dimension r_en(nelec,ncent)

      do 10 ic=1,ncent
        ict=iwctype(ic)

        r=r_en(iel,ic)
! local potential
        if(r.lt.rmax_coul(ict)) then
          call splfit_ps(r,lpotp1(ict),ict,vpot)
          vps(iel,ic,lpotp1(ict))=vpot
         else
          vps(iel,ic,lpotp1(ict))=-znuc(ict)/r
        endif
!       write(6,'(''ic,iel,r,vpot='',2i3,f6.3,f9.5)') ic,iel,r,vps(iel,ic,lpotp1(ict))
! non-local pseudopotential
        do 10 l=1,npotd(ict)
          if(l.ne.lpotp1(ict)) then
            if(r.lt.rmax_nloc(ict)) then
              call splfit_ps(r,l,ict,vpot)
              vps(iel,ic,l)=vpot
             else
              vps(iel,ic,l)=0.d0
            endif
          endif
   10 continue

      return
      end
!-----------------------------------------------------------------------

      subroutine splfit_ps(r,l,ict,vpot)
! get spline fit at r=r(iel,ic) of pseudopotential for center-type ict and
! angular momentum l stored on shifted exponential grid
! Note: I check if r < rmax_coul(ict) because this routine is called from
! ewald without going through getvps_tm.
! We assume that rmax_nloc(ict) <= rmax_coul(ict).

      use atom_mod
      use pseudo_mod
      use pseudo_tm_mod
      implicit real*8(a-h,o-z)

      if(r.lt.rmax_coul(ict)) then

        if(igrid_ps(ict).eq.1)then
! Warning: this may need fixing to h_ps
          stop 'May need fixing in splfit_ps'
          xr=(r-r0_ps(ict))/exp_h_ps(ict)+1
          jx=int(xr)
          ref0=r0_ps(ict)+exp_h_ps(ict)*(jx-1)
          ref1=ref0+exp_h_ps(ict)
          delh=exp_h_ps(ict)
         elseif(igrid_ps(ict).eq.2)then
          h_ps=dlog(exp_h_ps(ict))
          xr=dlog(r/r0_ps(ict))/h_ps+1
          jx=int(xr)
          ref0=r0_ps(ict)*exp_h_ps(ict)**(jx-1)
          ref1=ref0*exp_h_ps(ict)
          delh=ref1-ref0
         elseif(igrid_ps(ict).eq.3)then
          h_ps=dlog(exp_h_ps(ict))
          xr=(dlog((r+r0_ps(ict))/r0_ps(ict)))/h_ps+1
          jx=int(xr)
          ref0=r0_ps(ict)*exp_h_ps(ict)**(jx-1)-r0_ps(ict)
          ref1=(ref0+r0_ps(ict))*exp_h_ps(ict)-r0_ps(ict)
          delh=ref1-ref0
        endif

        if(jx.lt.1) then
          write(6,'(''ict,jx,xr,r,r0_ps(ict),exp_h_ps(ict),h_ps='',2i3,5d12.4)') &
     &    ict,jx,xr,r,r0_ps(ict),exp_h_ps(ict),h_ps
          write(6,'(''Warning: index < 1 in splfit_tm, r,xr='',2d12.4)') r,xr
          jx=1
        endif

        if(jx.gt.MPS_GRID) then
          write(6,'(''ict,jx,xr,r,r0_ps(ict),exp_h_ps(ict),h_ps='',2i3,5d12.4)') &
     &    ict,jx,xr,r,r0_ps(ict),exp_h_ps(ict),h_ps
          stop 'index > MPS_GRID in splfit_tm'
        endif

! cubic spline interpolation

        bb=(r-ref0)/delh
        aa=(ref1-r)/delh
        cc=aa*(aa**2-1.d0)*delh**2/6.d0
        dd=bb*(bb**2-1.d0)*delh**2/6.d0
! Take into account that for igrid=2 we added in an additional point at r=0.
        if(igrid_ps(ict).eq.2) jx=jx+1
        vpot=aa*vpseudo(jx,ict,l)+bb*vpseudo(jx+1,ict,l)+cc*d2pot(jx,ict,l)+dd*d2pot(jx+1,ict,l)
       else
        if(l.eq.lpotp1(ict)) then
          vpot=-znuc(ict)/r
         else
          vpot=0
        endif
      endif

      return
      end

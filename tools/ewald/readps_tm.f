      subroutine readps_tm
c Written by Claudia Filippi
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

      parameter (ncoef=5)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /pseudo_tm/ rmax_coul(MCTYPE),rmax_nloc(MCTYPE),arg_ps(MCTYPE),r0_ps(MCTYPE)
     &,vpseudo(MPS_GRID,MCTYPE,MPS_L),d2pot(MPS_GRID,MCTYPE,MPS_L),igrid_ps(MCTYPE),nr_ps(MCTYPE)
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc

      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad

      dimension r(MPS_GRID),y(ncoef),ce(ncoef),dmatr(ncoef*ncoef)
      dimension work(MPS_GRID)

      do 200 ic=1,nctype

        if(ic.lt.10) then
          write(atomtyp,'(i1)') ic
         elseif(ic.lt.100) then
          write(atomtyp,'(i2)') ic
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

c The TM psp. format has npotd and npotu for down and up, but we just use one of them
        if(nloc.eq.2) then
          read(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     &  (ititle(i),i=1,7),npotd(ic),npotu,nrm,r0_ps(ic),arg_ps(ic),zion
         elseif(nloc.eq.3) then
          read(1,'(a2,x,a2,x,a3,x,a4,x,13(a10,x))') nameat,icorr,irel,nicore
     &  ,(iray(i),i=1,6),(ititle(i),i=1,7)
          read(1,*) npotd(ic),npotu,nrm,r0_ps(ic),arg_ps(ic),zion
        endif
        write(6,'(''Center type'',i2,'' has'',i2,'' pseudopotential L components, and component''
     &  ,i2,'' is chosen to be local'')') ic,npotd(ic),lpotp1(ic)

        if(lpotp1(ic).gt.npotd(ic)) then
          write(6,'(''lpotp1(ic),npotd(ic)='',2i3)') lpotp1(ic).gt.npotd(ic)
          stop 'Cannot choose local psp. to be > number of l components, lpotp1(ic) > npotd(ic)'
        endif

        nr=nrm+1
        nr_ps(ic)=nr
        write(6,'(''ic,nr_ps(ic)'',9i5)') ic,nr_ps(ic)

        if(znuc(ic).ne.zion) then
          write(6,'(''znuc(ic) != zion in readps_tm'',2f6.1)') znuc(ic),zion
          stop 'znuc(ic) != zion in readps_tm'
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

        rmax_coul(ic)=0.d0
        jmax=0
        do 100 i=1,npotd(ic)
          if(nloc.eq.2) then
            read(1) ii,(vpseudo(j,ic,i),j=2,nr)
           elseif(nloc.eq.3) then
            read(1,*) ii,(vpseudo(j,ic,i),j=2,nr)
          endif
          vpseudo(1,ic,i)=0.d0
          do 50 j=nr,1,-1
c           write(6,'(''vps'',f8.4,f12.6,d12.4)') r(j),vpseudo(j,ic,i)/2,vpseudo(j,ic,i)/2+zion
            if(dabs(vpseudo(j,ic,i)+2.d0*zion).gt.1.d-6) then
              rmax_coul(ic)=max(rmax_coul(ic),r(j))
              jmax=max(jmax,j)
              goto 60
            endif
   50     continue
   60     do 100 j=2,nr
            vpseudo(j,ic,i)=0.5d0*vpseudo(j,ic,i)/r(j)
  100   continue
        jmax=jmax+5
        write(6,'(''center '',i3,'' pseudopot rmax,jmax= '',f6.2,i4)')
     &  ic,rmax_coul(ic),jmax

        close(1)

c If the local pseudopot component is not set in input, set it here
        if(lpotp1(ic).lt.0) then
          lpotp1(ic)=npotd(ic)
          write(6,'(''local pseudopot component is'',i3)') lpotp1(ic)
        endif

        arg_ps(ic)=dexp(arg_ps(ic))
        write(6,'(''potential grid param= '',2f10.4)') r0_ps(ic),arg_ps(ic)

        if(ipr.ge.1) then
          do 104 j=2,nr
  104       write(38,'(1pd12.6,9d14.6)') r(j),(vpseudo(j,ic,i),i=1,npotd(ic)),-znuc(ic)/r(j)
        endif

        do 110 i=1,npotd(ic)
          if(i.ne.lpotp1(ic)) then
            do 105 j=2,nr
              if(ipr.ge.1) then
              endif
  105         vpseudo(j,ic,i)=vpseudo(j,ic,i)-vpseudo(j,ic,lpotp1(ic))
          endif

  110   continue

        if(rmax_coul(ic).eq.0.d0) goto 200
        do 190 i=1,npotd(ic)
c small radii pot(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
          ll=0
          do 120 jj=1,ncoef
            y(jj)=vpseudo(jj+1,ic,i)
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
c               if(ir.eq.1) vpseudo(ir,ic,i)=val
  130           write(6,'(i2,1p3e16.8,1p1e11.2)')
     &          ir,r(ir),val,vpseudo(ir,ic,i),val-vpseudo(ir,ic,i)
          endif

          vpseudo(1,ic,i)=ce(1)

          dpot1=ce(2)
          do 135 icoef=3,ncoef
  135     dpot1=dpot1+(icoef-1)*ce(icoef)*r(1)**(icoef-2)

          dpotn=0.d0
          if(i.eq.lpotp1(ic)) dpotn=zion/r(jmax)**2

          write(6,'(''dpot1,dpotn'',1p2e15.5)') dpot1,dpotn

c get second derivative for spline fit
c Warning: this is done only upto jmax, not upto nr, so when using spline
c make sure it is not used beyond jmax.
          call spline2(r,vpseudo(1,ic,i),jmax,dpot1,dpotn,
     &    d2pot(1,ic,i),work)

          do 190 j=1,nr
            if(ipr.ge.1) then
              if(i.eq.1) write(35,'(1p5d14.6)') r(j),vpseudo(j,ic,i),d2pot(j,ic,i),-znuc(ic)/r(j),-2*znuc(ic)/r(j)**3
              if(i.eq.2) write(36,'(1p5d14.6)') r(j),vpseudo(j,ic,i),d2pot(j,ic,i),-znuc(ic)/r(j),-2*znuc(ic)/r(j)**3
              if(i.eq.3) write(37,'(1p5d14.6)') r(j),vpseudo(j,ic,i),d2pot(j,ic,i),-znuc(ic)/r(j),-2*znuc(ic)/r(j)**3
c             if(i.eq.1.and.i.ne.lpotp1(ic)
c    &        write(38,'(1p5d14.6)') r(j),vpseudo(j,ic,i)+vpseudo(j,ic,lpotp1(ic)),vpseudo(j,ic,lpotp1(ic)),-znuc(ic)/r(j)
c             if(i.eq.2.and.i.ne.lpotp1(ic))
c    &        write(39,'(1p5d14.6)') r(j),vpseudo(j,ic,i)+vpseudo(j,ic,lpotp1(ic)),vpseudo(j,ic,lpotp1(ic)),-znuc(ic)/r(j)
            endif
  190   continue

c Warning: Temporary print of psp
c       write(20,'(2f9.5,'' znuc,rpotc'')') znuc(1)
c       write(20,'(i6,f18.14,'' nr_ps(ict),arg_ps(ict)'')') nr,arg_ps(ic)
c       do 195 ir=1,nr
c 195     write(20,'(2d20.12)') r(ir),vpseudo(ir,ic,lpotp1(ic))

  200 continue

      call gesqua(nquad,xq0,yq0,zq0,wq)
c     call gesqua(nquad,xq,yq,zq,wq)

      write(6,'(''quadrature points'')')
      do 210 i=1,nquad
  210   write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

      do 240 j=1,nr
  240   if(r(j).gt.2..and.r(j).lt.5.)
     & write(6,'(f6.3,2f10.5,2d10.2)') r(j),4/r(j),vpseudo(j,1,lpotp1(1)),4/r(j)+vpseudo(j,1,lpotp1(1)),d2pot(j,1,lpotp1(1))

      return

  999 write(6,'(''Error: Pseudopot. file '',a20,'' is missing'')') filename
      stop 'Pseudopot. file is missing'

      end
c-----------------------------------------------------------------------

c compute tm-pseudopotential for electron iel
      subroutine getvps_tm(r_en,iel)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /pseudo_tm/ rmax_coul(MCTYPE),rmax_nloc(MCTYPE),arg_ps(MCTYPE),r0_ps(MCTYPE)
     &,vpseudo(MPS_GRID,MCTYPE,MPS_L),d2pot(MPS_GRID,MCTYPE,MPS_L),igrid_ps(MCTYPE),nr_ps(MCTYPE)
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc

      dimension r_en(MELEC,MCENT)

      do 10 ic=1,ncent
        ict=iwctype(ic)

        r=r_en(iel,ic)
c local potential
        if(r.lt.rmax_coul(ict)) then
          call splfit_tm(r,lpotp1(ict),ict,vpot)
          vps(iel,ic,lpotp1(ict))=vpot
         else
          vps(iel,ic,lpotp1(ict))=-znuc(ict)/r
        endif
c non-local pseudopotential
        do 10 l=1,npotd(ict)
          if(l.ne.lpotp1(ict)) then
            if(r.lt.rmax_coul(ict)) then
              call splfit_tm(r,l,ict,vpot)
              vps(iel,ic,l)=vpot
             else
              vps(iel,ic,l)=0.d0
            endif
          endif
   10 continue

      return
      end
c-----------------------------------------------------------------------

      subroutine splfit_tm(r,l,ict,vpot)
c get spline_fit at r=r(iel,ic) of TM potential for center-type ict and
c angular momentum l stored on shifted exponential grid
c Note: I check if r < rmax_coul(ict) because this routine is called from
c ewald without going through getvps_tm.

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /pseudo_tm/ rmax_coul(MCTYPE),rmax_nloc(MCTYPE),arg_ps(MCTYPE),r0_ps(MCTYPE)
     &,vpseudo(MPS_GRID,MCTYPE,MPS_L),d2pot(MPS_GRID,MCTYPE,MPS_L),igrid_ps(MCTYPE),nr_ps(MCTYPE)

c     if(r.lt.rmax_coul(ict)) then
        dlogag=dlog(arg_ps(ict))
        xr=(dlog((r+r0_ps(ict))/r0_ps(ict)))/dlogag+1.d0
        jx=int(xr)
        if(jx.lt.1) then
          write(6,'(''ict,jx,xr,r,r0_ps(ict),arg_ps(ict),dlogag='',2i3,5d12.4)')
     &    ict,jx,xr,r,r0_ps(ict),arg_ps(ict),dlogag
          write(6,'(''Warning: index < 1 in splfit_tm, r,xr='',2d12.4)') r,xr
          jx=1
        endif

        if(jx.gt.MPS_GRID) then
          write(6,'(''ict,jx,xr,r,r0_ps(ict),arg_ps(ict),dlogag='',2i3,5d12.4)')
     &    ict,jx,xr,r,r0_ps(ict),arg_ps(ict),dlogag
          stop 'index > MPS_GRID in splfit_tm'
        endif

        ref0=r0_ps(ict)*arg_ps(ict)**(jx-1)-r0_ps(ict)
        ref1=(ref0+r0_ps(ict))*arg_ps(ict)-r0_ps(ict)
        delh=ref1-ref0

c cubic spline interpolation

        bb=(r-ref0)/delh
        aa=(ref1-r)/delh
        cc=aa*(aa**2-1.d0)*delh**2/6.d0
        dd=bb*(bb**2-1.d0)*delh**2/6.d0
        vpot=aa*vpseudo(jx,ict,l)+bb*vpseudo(jx+1,ict,l)+
     &  cc*d2pot(jx,ict,l)+dd*d2pot(jx+1,ict,l)
c      else
c       vpot=-znuc(ict)/r
c     endif

      return
      end

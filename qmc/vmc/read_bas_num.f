      subroutine read_bas_num(iwf)
c Written by Claudia Filippi. Modified by Cyrus Umrigar
c Reads in localized orbitals on a radial grid
c If igrid(ic).eq.2 .and. (r0_bas(ic).le.0.d0 .or. exp_h_bas(ic).le.0.d0) then grid parameters are deduced later from r values

      use atom_mod
      use const_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'numbas.h'
!JT      include 'pseudo.h'
!JT      include 'force.h'

      character*20 filename,wforce,atomtyp
      character*80 title
c     character*20 lcent

      parameter(NCOEF=5)

      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
      common /numexp/ce(NCOEF,MRWF,MCTYPE,MFORCE),ae(2,MRWF,MCTYPE,MFORCE)
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension x(MRWF_PTS),work(MRWF_PTS),y(NCOEF),dmatr(NCOEF*NCOEF)
     &,l(MRWF),icusp(MCTYPE)

c nrbas = number of radial basis functions for each center

c igrid = 1 linear              r(i+1)=h_bas+r(i), r(1)=r0_bas
c         2 exponential         r(i+1)=exp_h_bas*r(i), r(1)=r0_bas
c         3 shifted exponential r(i+1)=r0_bas*(exp_h_bas**(i-1)-1)
c           r(n) is read in, r0_bas=r(n)/(exp_h_bas**(nr-1)-1)

      if(iwf.lt.10) then
        write(wforce,'(i1)') iwf
       elseif(iwf.lt.100) then
        write(wforce,'(i2)') iwf
       else
        stop 'read_bas_num, wforce >= 100'
      endif

      write(6,'(/,''Reading numerical radial basis functions'')')
      do 100 ic=1,nctype

        if(ic.lt.10) then
          write(atomtyp,'(i1)') ic
         elseif(ic.lt.100) then
          write(atomtyp,'(i2)') ic
         else
          stop 'read_bas_num, problem atomtyp'
        endif

        filename='basis.'//atomtyp(1:index(atomtyp,' ')-1)
        if(iwf.ge.2)
     &  filename=filename(1:index(filename,' ')-1)//'.'//wforce
        open(21,file=filename,status='old',err=999)

c position file to skip comments
        write(6,'(''Reading numerical basis function file'')')
        title(1:1)='#'
        do while(title(1:1).eq.'#')
          read(21,'(a80)') title
c         write(6,'(a80)') title
        enddo

c       read(21,*) nrbas(ic),igrid(ic),nr(ic),exp_h_bas(ic),r0_bas(ic),icusp(ic)
        read(title,*) nrbas(ic),igrid(ic),nr(ic),exp_h_bas(ic),r0_bas(ic),icusp(ic)
        write(6,'(''ic,nrbas,igrid,nr,exp_h_bas,r0_bas,icusp=''i3,3i5,2f10.6,i3)')
     &  ic,nrbas(ic),igrid(ic),nr(ic),exp_h_bas(ic),r0_bas(ic),icusp(ic)

        if(nrbas(ic).gt.MRWF) stop 'nrbas gt MRWF'
        if(nr(ic).gt.MRWF_PTS) stop 'nr gt MRWF_PTS'
        if(igrid(ic).ne.1.and.igrid(ic).ne.2.and.igrid(ic).ne.3)
     &  stop 'grid not implemented'
        if(igrid(ic).eq.3 .and. r0_bas(ic).gt.0.d0 .and. r0_bas(ic).lt.1.d0)
     &  stop 'For igrid=3 r0_bas should be the last point not the first point on the grid'
        if(igrid(ic).eq.2 .and. r0_bas(ic).gt.0.d0 .and. r0_bas(ic).gt.1.d0)
     &  stop 'For igrid=2 r0_bas should be the last point not the first point on the grid'

        if(nloc.eq.0) read(21,*) (l(irb),irb=1,nrbas(ic))

        do 10 ir=1,nr(ic)
          read(21,*) x(ir),(rwf(ir,irb,ic,iwf),irb=1,nrbas(ic))
   10     continue

        do 12 irb=1,nrbas(ic)
   12     write(6,'(''center'',i3,'' numerical radial basis'',i3,'' has l='',9f8.5)')
     &    ic,irb,log(rwf(2,irb,ic,iwf)/rwf(1,irb,ic,iwf))/log(x(2)/x(1)),
     &    log(rwf(3,irb,ic,iwf)/rwf(2,irb,ic,iwf))/log(x(3)/x(2))

        if(igrid(ic).eq.2.and.exp_h_bas(ic).le.1.d0) exp_h_bas(ic)=x(2)/x(1)
        if(igrid(ic).eq.3) r0_bas(ic)=r0_bas(ic)/(exp_h_bas(ic)**(nr(ic)-1)-1.d0)

        if(igrid(ic).eq.2 .and. (r0_bas(ic).le.0.d0 .or. exp_h_bas(ic).le.0.d0)) then
          r0_bas(ic)=x(1)
c         r0_bas(ic)=x(nr(ic))
          exp_h_bas(ic)=x(2)/x(1)
          write(6,'(''Grid parameters deduced from grid values are, r0_bas(ic),exp_h_bas(ic)='',9f10.5)')
     &    r0_bas(ic),exp_h_bas(ic)
        endif

        do 100 irb=1,nrbas(ic)

        if(nloc.eq.0.and.l(irb).eq.0.and.icusp(ic).eq.1) then

c small radii wf(r)=ce1-znuc*ce1*r+ce3*r**2+ce4*r**3+ce5*r**4
          do 15 ii=1,ncoef-1
   15       dmatr(ii)=1.d0-znuc(ic)*x(ii)
          y(1)=rwf(1,irb,ic,iwf)
          ll=ncoef-1
          do 16 jj=2,ncoef-1
            y(jj)=rwf(jj,irb,ic,iwf)
            do 16 ii=2,ncoef-1
              ll=ll+1
   16         dmatr(ll)=x(ii)**jj

          call dgelg(y,dmatr,ncoef-1,1,1.d-8,ier)
          ce(1,irb,ic,iwf)=y(1)
          ce(2,irb,ic,iwf)=-znuc(ic)*ce(1,irb,ic,iwf)
          ce(3,irb,ic,iwf)=y(2)
          ce(4,irb,ic,iwf)=y(3)
          ce(5,irb,ic,iwf)=y(4)

         else

c small radii wf(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
          ll=0
          do 25 jj=1,ncoef
            y(jj)=rwf(jj,irb,ic,iwf)
            do 25 ii=1,ncoef
              ll=ll+1
   25         dmatr(ll)=x(ii)**(jj-1)
          call dgelg(y,dmatr,ncoef,1,1.d-8,ier)

          do 26 icoef=1,ncoef
   26       ce(icoef,irb,ic,iwf)=y(icoef)

        endif

        if(ipr.ge.1) then
          write(6,'(''coefficients'',1p10d22.10)')
     &             (ce(iff,irb,ic,iwf),iff=1,ncoef)
          write(6,'(''check the small radius expansion of radial basis fn'',i3)') irb
          write(6,'(''irad, rad, extrapolated value, correct value'')')
        endif
        do 30 ir=1,10
          val=ce(1,irb,ic,iwf)
          do 28 icoef=2,ncoef
   28     val=val+ce(icoef,irb,ic,iwf)*x(ir)**(icoef-1)
          if(ipr.ge.1) write(6,'(i2,1p3d22.14,1pd8.0)')ir,x(ir),val,rwf(ir,irb,ic,iwf)
     &    ,val-rwf(ir,irb,ic,iwf)
          if(abs(val-rwf(ir,irb,ic,iwf))/rwf(ir,irb,ic,iwf).gt.1.d-2 .and. rwf(ir,irb,ic,iwf).ne.0.d0) then
            write(6,'(''irb,ir,val,rwf(ir,irb,ic,iwf)'',2i5,9d12.4)') irb,ir,val,rwf(ir,irb,ic,iwf)
            write(6,'(''Warning: fit of radial function at small radii not good'')')
c           stop 'fit of radial function at small radii not good'
          endif
   30   continue

        dwf1=0.d0
        do 32 icoef=2,ncoef
   32     dwf1=dwf1+(icoef-1)*ce(icoef,irb,ic,iwf)*x(1)**(icoef-2)

c large radii wf(r)=a0*exp(-ak*r) for ndim=3
c             wf(r)=a0*exp(-ak*r^2) for ndim=2, ak should give weff/2
c       xm=0.5d0*(x(nr(ic))+x(nr(ic)-1))
        wfm=0.5d0*(rwf(nr(ic),irb,ic,iwf)+rwf(nr(ic)-1,irb,ic,iwf))
        dwfm=(rwf(nr(ic),irb,ic,iwf)-rwf(nr(ic)-1,irb,ic,iwf))/
     &  (x(nr(ic))-x(nr(ic)-1))
        if(ndim.eq.3) then
          if(dabs(wfm).gt.1.d-99) then
            ae(2,irb,ic,iwf)=-dwfm/wfm
            ae(1,irb,ic,iwf)=rwf(nr(ic),irb,ic,iwf)*
     &                     dexp(ae(2,irb,ic,iwf)*x(nr(ic)))
            dwfn=-ae(2,irb,ic,iwf)*rwf(nr(ic),irb,ic,iwf)
          else
            ae(1,irb,ic,iwf)=0.d0
            ae(2,irb,ic,iwf)=0.d0
            dwfn=0.d0
          endif
        elseif(ndim.eq.2) then
          if(dabs(wfm).gt.1.d-99) then
            ae(2,irb,ic,iwf)=-0.5d0*dwfm/(wfm*x(nr(ic)))
c            ae(2,irb,ic,iwf)=we/2   !  correct expression for parabolic confinement
            ae(1,irb,ic,iwf)=rwf(nr(ic),irb,ic,iwf)*
     &                     dexp(ae(2,irb,ic,iwf)*x(nr(ic))*x(nr(ic)))
            dwfn=-2*x(nr(ic))*ae(2,irb,ic,iwf)*rwf(nr(ic),irb,ic,iwf)
          else
            ae(1,irb,ic,iwf)=0.d0
            ae(2,irb,ic,iwf)=0.d0
            dwfn=0.d0
          endif
        else
          stop 'ndim must be 2 or 3 in read_bas_num'
        endif

        if(ipr.ge.1) then
          write(6,'(''a0,ak'',1p2d22.10)')
     &                ae(1,irb,ic,iwf),ae(2,irb,ic,iwf)
          write(6,'(''check the large radius expansion'')')
          write(6,'(''irad, rad, extrapolated value, correct value'')')
        endif
        do 40 ir=1,10
          if(ndim.eq.3) then
            val=ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*x(nr(ic)-ir))
          elseif(ndim.eq.2) then
            val=ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*x(nr(ic)-ir)*x(nr(ic)-ir))
          else
            stop 'ndim must be 2 or 3 in read_bas_num'
          endif
          if(ipr.ge.1) write(6,'(i2,1p3d22.14,1pd8.0)')
     &    ir,x(nr(ic)-ir),val,rwf(nr(ic)-ir,irb,ic,iwf)
     &    ,val-rwf(nr(ic)-ir,irb,ic,iwf)
          if(abs(val-rwf(nr(ic)-ir,irb,ic,iwf))/rwf(nr(ic)-ir,irb,ic,iwf).gt.1.d0)
     &       write(*,*) 'Warning: fit of radial function at large radii not good'
c     &    stop 'fit of radial function at large radii not good'
   40   enddo
        if(ipr.ge.1) write(6,*) 'dwf1,dwfn',dwf1,dwfn

        call spline2(x,rwf(1,irb,ic,iwf),nr(ic),dwf1,dwfn,
     &  d2rwf(1,irb,ic,iwf),work)

      close(21)
  100 continue

c TEMPORARY debug
c     do 130 jwf=1,nforce
c       do 130 ic=1,nctype
c         if(ic.lt.10) then
c           write(lcent,'(i1)') ic
c          elseif(ic.lt.100) then
c           write(lcent,'(i2)') ic
c         else
c           stop 'problem with spline.test'
c         endif
c         if(jwf.lt.10) then
c           write(wforce,'(i1)') jwf
c          elseif(jwf.lt.100) then
c           write(wforce,'(i2)') jwf
c         endif
c         filename='spline.chk.'//lcent(1:index(lcent,' ')-1)//wforce
c         open(22,file=filename,status='unknown')
c         do 120 ir=1,nr(ic)-1
c           ii=1
c           do 110 irb=1,nrbas(ic)
c             call splfit_bas(x(ir),irb,ic,jwf,work(ii),1)
c110          ii=ii+3
c120        write(22,'(1p40e20.12)') x(ir),(work(j),j=1,ndim*nrbas(ic))
c130      close(22)

      return

  999 write(6,'(''Error: file '',a20,'' is missing'')') filename
      stop 'Numerical basis function file is missing'

      end

      subroutine read_bas_num(iwf)
c Written by Claudia Filippi. Modified by Cyrus Umrigar
c Reads in localized orbitals on a radial grid
c If igrid(ict).eq.2 .and. (r0_bas(ict).le.0.d0 .or. exp_h_bas(ict).le.0.d0) then grid parameters are deduced later from r values
      use all_tools_mod
      use atom_mod
      use const_mod
      use dim_mod
      use coefs_mod,  only : nbasis
      use basis1_mod, only : zex
      use basis2_mod, only : nbasis_ctype, iwrwf2
      use numbas_mod
      use pseudo_mod
      use forcepar_mod
      use numexp_mod
      implicit real*8(a-h,o-z)

      character*20 filename,wforce,atomtyp
      character*80 title
c     character*20 lcent

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      double precision, allocatable :: x(:),work(:),rwf_tmp(:,:)
      integer, allocatable :: iwrwf_tmp(:)
      integer, allocatable :: l(:)
      dimension y(NCOEF),dmatr(NCOEF*NCOEF),icusp(nctype)

c nrbas_numerical = number of radial basis functions for each centertype

c igrid = 1 linear              r(i+1)=h_bas+r(i), r(1)=r0_bas
c         2 exponential         r(i+1)=exp_h_bas*r(i), r(1)=r0_bas
c         3 shifted exponential r(i+1)=r0_bas*(exp_h_bas**(i-1)-1)
c           r(n) is read in, r0_bas=r(n)/(exp_h_bas**(nr-1)-1)

c If there are no numerical radial basis functions just initialize nrbas_numerical, nrbas and return
      call alloc ('nrbas_numerical', nrbas_numerical, nctype)
      call alloc ('nrbas', nrbas, nctype)
      if(minval(zex(:,1)).ne.0.d0) then
        nrbas_numerical=0
        nrbas=nrbas_analytical
        return
      endif
        
      MRWF =  0
      MRWF_PTS = 0

      if(iwf.lt.10) then
        write(wforce,'(i1)') iwf
       elseif(iwf.lt.100) then
        write(wforce,'(i2)') iwf
       else
        stop 'read_bas_num, wforce >= 100'
      endif

      call alloc ('exp_h_bas', exp_h_bas, nctype)
      call alloc ('r0_bas', r0_bas, nctype)
      call alloc ('igrid', igrid, nctype)
      call alloc ('nr', nr, nctype)

      write(6,'(''beginning of read_bas_num: iwrwf2='',40i3)') (iwrwf2(ib),ib=1,nbasis)

      write(6,'(''iwrwf2 before  compact='',40i3)') (iwrwf2(ib),ib=1,nbasis)

      ib=0
      do 110 ict=1,nctype

        write(6,'(''beginning of read_bas_num: ict, iwrwf='',i3,3x,40i3)') ict, (iwrwf(ibct,ict),ibct=1,nbasis_ctype(ict))

        if(ict.lt.10) then
          write(atomtyp,'(i1)') ict
         elseif(ict.lt.100) then
          write(atomtyp,'(i2)') ict
         else
          stop 'read_bas_num, problem atomtyp'
        endif

        filename='basis.'//atomtyp(1:index(atomtyp,' ')-1)
        if(iwf.ge.2) filename=filename(1:index(filename,' ')-1)//'.'//wforce
        open(21,file=filename,status='old',err=999)

c position file to skip comments
        write(6,'(/3a)') ' Reading numerical radial basis function file >',trim(filename),'<'
        title(1:1)='#'
        do while(title(1:1).eq.'#')
          read(21,'(a80)') title
c         write(6,'(a80)') title
        enddo

c       read(21,*) nrbas_numerical(ict),igrid(ict),nr(ict),exp_h_bas(ict),r0_bas(ict),icusp(ict)
        read(title,*) nrbas_numerical(ict),igrid(ict),nr(ict),exp_h_bas(ict),r0_bas(ict),icusp(ict)
        write(6,'('' ict,nrbas_numerical,igrid,nr,exp_h_bas,r0_bas,icusp=''i3,3i5,2f10.6,i3)')
     &  ict,nrbas_numerical(ict),igrid(ict),nr(ict),exp_h_bas(ict),r0_bas(ict),icusp(ict)

        MRWF = max (MRWF, nrbas_analytical(ict)+nrbas_numerical(ict))
        MRWF_PTS = max (MRWF_PTS, nr(ict))

        if(igrid(ict).ne.1.and.igrid(ict).ne.2.and.igrid(ict).ne.3) stop 'grid not implemented'
        if(igrid(ict).eq.3 .and. r0_bas(ict).gt.0.d0 .and. r0_bas(ict).lt.1.d0)
     &  stop 'For igrid=3 r0_bas should be the last point not the first point on the grid'
        if(igrid(ict).eq.2 .and. r0_bas(ict).gt.0.d0 .and. r0_bas(ict).gt.1.d0)
     &  stop 'For igrid=2 r0_bas should be the last point not the first point on the grid'

        call alloc ('l', l, MRWF)
        if(nloc.eq.0) read(21,*) (l(irb),irb=1,nrbas_numerical(ict))

        call alloc ('x', x, MRWF_PTS)
        call alloc ('rwf', rwf, MRWF_PTS, MRWF, nctype, nwf)
        do 10 ir=1,nr(ict)
          read(21,*) x(ir),(rwf(ir,irb,ict,iwf),irb=1,nrbas_numerical(ict))
   10     continue

c Compactify numerical basis read in from basis files, calculate nrbas_numerical and modify iwrwf correspondingly.
c Those basis functions with zex(ib)==0.d0 are numerical contracted functions
        write(6,'(''iwrwf before compact='',40i3)') (iwrwf(ibct,ict),ibct=1,nbasis_ctype(ict))
        call alloc ('rwf_tmp', rwf_tmp, nr(ict), nrbas_numerical(ict))
        call alloc ('iwrwf_tmp', iwrwf_tmp, nbasis_ctype(ict))

        do 12 ibct=1,nbasis_ctype(ict)
   12     iwrwf_tmp(ibct)=iwrwf(ibct,ict)
        nrbas_numerical(ict)=0
        do 20 ibct=1,nbasis_ctype(ict)
          ib=ib+1
          irb=iwrwf(ibct,ict)
          if(zex(ib,iwf).eq.0.d0) then
            do 14 jbct=1,ibct-1
              if(iwrwf(ibct,ict).eq.iwrwf(jbct,ict)) then
                jrb=iwrwf_tmp(jbct)-nrbas_analytical(ict)
                goto 18
              endif
   14       continue
            nrbas_numerical(ict)=nrbas_numerical(ict)+1
            jrb=nrbas_numerical(ict)
            do 16 ir=1,nr(ict)
   16         rwf_tmp(ir,jrb)=rwf(ir,irb,ict,iwf)
   18       iwrwf_tmp(ibct)=nrbas_analytical(ict)+jrb
            do 19 icent=1,ncent_ctype(ict)
   19         iwrwf2((icent-1)*nbasis_ctype(ict)+ib)=iwrwf_tmp(ibct)
          endif
   20   continue
        nrbas(ict)=nrbas_analytical(ict)+nrbas_numerical(ict)
        do 22 irb=1,nrbas_numerical(ict)
          do 22 ir=1,nr(ict)
   22       rwf(ir,irb,ict,iwf)=rwf_tmp(ir,irb)
        do 23 ibct=1,nbasis_ctype(ict)
          iwrwf(ibct,ict)=iwrwf_tmp(ibct)
          if(iwrwf(ibct,ict).gt.nrbas(ict)) then
            write(6,'(''ict,ibct,iwrwf(ib,ict),nrbas(ict)'',9i3)') ict,ibct,iwrwf(ibct,ict),nrbas(ict)
            stop 'iwrwf(ib,ict) > nrbas(ict)'
          endif
   23   continue
        deallocate(rwf_tmp,stat=istat)
        deallocate(iwrwf_tmp,stat=istat2)
        if(istat.ne.0 .or. istat2.ne.0) then
          write(6,'(''read_orb_bas: Error deallocating rwf_tmp or iwrwf_tmp, istat='',i2)') istat
          stop 'read_orb_bas: Error deallocating rwf_tmp or iwrwf_tmp'
        endif
        write(6,'(''iwrwf after  compact='',40i3)') (iwrwf(ibct,ict),ibct=1,nbasis_ctype(ict))
        write(6,'(''iwrwf2 after  compact='',40i3)') (iwrwf2(ib2),ib2=1,nbasis)

        do 25 irb=1,nrbas_numerical(ict)
   25     write(6,'('' center'',i3,'' numerical radial basis'',i3,'' has l='',9f8.5)') ict,irb,
c    &    log(rwf(2,irb,ict,iwf)/rwf(1,irb,ict,iwf))/log(x(2)/x(1)),
     &    log(rwf(3,irb,ict,iwf)/rwf(2,irb,ict,iwf))/log(x(3)/x(2))

        if(igrid(ict).eq.2.and.exp_h_bas(ict).le.1.d0) exp_h_bas(ict)=x(2)/x(1)
        if(igrid(ict).eq.3) r0_bas(ict)=r0_bas(ict)/(exp_h_bas(ict)**(nr(ict)-1)-1.d0)

        if(igrid(ict).eq.2 .and. (r0_bas(ict).le.0.d0 .or. exp_h_bas(ict).le.0.d0)) then
          r0_bas(ict)=x(1)
c         r0_bas(ict)=x(nr(ict))
          exp_h_bas(ict)=x(2)/x(1)
          write(6,'('' Grid parameters deduced from grid values are, r0_bas(ict),exp_h_bas(ict)='',9f10.5)')
     &    r0_bas(ict),exp_h_bas(ict)
        endif

        call alloc ('ce', ce, NCOEF, MRWF, nctype, nwf)
        call alloc ('ae', ae, 2, MRWF, nctype, nwf)

        do 100 irb=1,nrbas_numerical(ict)

          if(nloc.eq.0.and.l(irb).eq.0.and.icusp(ict).eq.1) then

c small radii wf(r)=ce1-znuc*ce1*r+ce3*r**2+ce4*r**3+ce5*r**4
            do 40 ii=1,ncoef-1
   40         dmatr(ii)=1.d0-znuc(ict)*x(ii)
            y(1)=rwf(1,irb,ict,iwf)
            ll=ncoef-1
            do 45 jj=2,ncoef-1
              y(jj)=rwf(jj,irb,ict,iwf)
              do 45 ii=2,ncoef-1
                ll=ll+1
   45           dmatr(ll)=x(ii)**jj

            call dgelg(y,dmatr,ncoef-1,1,1.d-8,ier)
            ce(1,irb,ict,iwf)=y(1)
            ce(2,irb,ict,iwf)=-znuc(ict)*ce(1,irb,ict,iwf)
            ce(3,irb,ict,iwf)=y(2)
            ce(4,irb,ict,iwf)=y(3)
            ce(5,irb,ict,iwf)=y(4)

           else

c small radii wf(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
            ll=0
            do 50 jj=1,ncoef
              y(jj)=rwf(jj,irb,ict,iwf)
              do 50 ii=1,ncoef
                ll=ll+1
   50           dmatr(ll)=x(ii)**(jj-1)
            call dgelg(y,dmatr,ncoef,1,1.d-8,ier)

            do 55 icoef=1,ncoef
   55         ce(icoef,irb,ict,iwf)=y(icoef)

          endif

          if(ipr.ge.1) then
            write(6,'('' coefficients'',1p10d22.10)') (ce(iff,irb,ict,iwf),iff=1,ncoef)
            write(6,'('' check the small radius expansion of radial basis fn'',i3)') irb
            write(6,'('' irad, rad, extrapolated value, correct value'')')
          endif
          do 70 ir=1,10
            val=ce(1,irb,ict,iwf)
            do 60 icoef=2,ncoef
   60       val=val+ce(icoef,irb,ict,iwf)*x(ir)**(icoef-1)
            if(ipr.ge.1) write(6,'(i2,1p3d22.14,1pd8.0)')ir,x(ir),val,rwf(ir,irb,ict,iwf),val-rwf(ir,irb,ict,iwf)
            if(abs(val-rwf(ir,irb,ict,iwf))/rwf(ir,irb,ict,iwf).gt.1.d-2 .and. rwf(ir,irb,ict,iwf).ne.0.d0) then
              write(6,'('' irb,ir,val,rwf(ir,irb,ict,iwf)'',2i5,9d12.4)') irb,ir,val,rwf(ir,irb,ict,iwf)
              write(6,'('' Warning: fit of radial function at small radii not good'')')
c             stop 'fit of radial function at small radii not good'
            endif
   70     continue

          dwf1=0.d0
          do 80 icoef=2,ncoef
   80     dwf1=dwf1+(icoef-1)*ce(icoef,irb,ict,iwf)*x(1)**(icoef-2)

c large radii wf(r)=a0*exp(-ak*r) for ndim=3
c         wf(r)=a0*exp(-ak*r^2) for ndim=2, ak should give weff/2
c         xm=0.5d0*(x(nr(ict))+x(nr(ict)-1))
          wfm=0.5d0*(rwf(nr(ict),irb,ict,iwf)+rwf(nr(ict)-1,irb,ict,iwf))
          dwfm=(rwf(nr(ict),irb,ict,iwf)-rwf(nr(ict)-1,irb,ict,iwf))/
     &    (x(nr(ict))-x(nr(ict)-1))
          if(ndim.eq.3) then
            if(dabs(wfm).gt.1.d-99) then
              ae(2,irb,ict,iwf)=-dwfm/wfm
              ae(1,irb,ict,iwf)=rwf(nr(ict),irb,ict,iwf)*
     &                       dexp(ae(2,irb,ict,iwf)*x(nr(ict)))
              dwfn=-ae(2,irb,ict,iwf)*rwf(nr(ict),irb,ict,iwf)
            else
              ae(1,irb,ict,iwf)=0.d0
              ae(2,irb,ict,iwf)=0.d0
              dwfn=0.d0
            endif
          elseif(ndim.eq.2) then
            if(dabs(wfm).gt.1.d-99) then
              ae(2,irb,ict,iwf)=-0.5d0*dwfm/(wfm*x(nr(ict)))
c              ae(2,irb,ict,iwf)=we/2   !  correct expression for parabolic confinement
              ae(1,irb,ict,iwf)=rwf(nr(ict),irb,ict,iwf)*
     &                       dexp(ae(2,irb,ict,iwf)*x(nr(ict))*x(nr(ict)))
              dwfn=-2*x(nr(ict))*ae(2,irb,ict,iwf)*rwf(nr(ict),irb,ict,iwf)
            else
              ae(1,irb,ict,iwf)=0.d0
              ae(2,irb,ict,iwf)=0.d0
              dwfn=0.d0
            endif
          else
            stop 'ndim must be 2 or 3 in read_bas_num'
          endif

          if(ipr.ge.1) then
            write(6,'('' a0,ak'',1p2d22.10)')
     &                  ae(1,irb,ict,iwf),ae(2,irb,ict,iwf)
            write(6,'('' check the large radius expansion'')')
            write(6,'('' irad, rad, extrapolated value, correct value'')')
          endif
          do 90 ir=1,10
            if(ndim.eq.3) then
              val=ae(1,irb,ict,iwf)*dexp(-ae(2,irb,ict,iwf)*x(nr(ict)-ir))
            elseif(ndim.eq.2) then
              val=ae(1,irb,ict,iwf)*dexp(-ae(2,irb,ict,iwf)*x(nr(ict)-ir)*x(nr(ict)-ir))
            else
              stop 'ndim must be 2 or 3 in read_bas_num'
            endif
            if(ipr.ge.1) write(6,'(i2,1p3d22.14,1pd8.0)')
     &      ir,x(nr(ict)-ir),val,rwf(nr(ict)-ir,irb,ict,iwf)
     &      ,val-rwf(nr(ict)-ir,irb,ict,iwf)
            if(abs(val-rwf(nr(ict)-ir,irb,ict,iwf))/rwf(nr(ict)-ir,irb,ict,iwf).gt.1.d0)
     &         write(*,*) 'Warning: fit of radial function at large radii not good'
c    &         stop 'fit of radial function at large radii not good'
   90     enddo
          if(ipr.ge.1) write(6,*) 'dwf1,dwfn',dwf1,dwfn

          call alloc ('d2rwf', d2rwf, MRWF_PTS, MRWF, nctype, nwf)
          call alloc ('work', work, MRWF_PTS)
  100     call spline2(x,rwf(1,irb,ict,iwf),nr(ict),dwf1,dwfn,d2rwf(1,irb,ict,iwf),work)

        ib=ib+(ncent_ctype(ict)-1)*nbasis_ctype(ict)
        close(21)
  110 continue
      write(6,'(''end       of read_bas_num: iwrwf2='',40i3)') (iwrwf2(ib),ib=1,nbasis)

c TEMPORARY debug
c     do 130 jwf=1,nforce
c       do 130 ict=1,nctype
c         if(ict.lt.10) then
c           write(lcent,'(i1)') ict
c          elseif(ict.lt.100) then
c           write(lcent,'(i2)') ict
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
c         do 120 ir=1,nr(ict)-1
c           ii=1
c           do 110 irb=1,nrbas_numerical(ict)
c             call splfit_bas(x(ir),irb,ict,jwf,work(ii),1)
c110          ii=ii+3
c120        write(22,'(1p40e20.12)') x(ir),(work(j),j=1,ndim*nrbas_numerical(ict))

      return

  999 write(6,'(''Error: file '',a20,'' is missing'')') filename
      stop 'Numerical basis function file is missing'

      end

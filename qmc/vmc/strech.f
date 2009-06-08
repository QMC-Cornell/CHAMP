      subroutine strech(x,xstrech,ajacob,ifr,istrech_el)
c Written by Cyrus Umrigar and Claudia Filippi
c Uses the coordinate transform described in:
c 1) Two Aspects of Quantum Monte Carlo: Determination of Accurate Wavefunctions and
c    Determination of Potential Energy Surfaces of Molecules, C.J. Umrigar,
c    Int. J. Quant. Chem. Symp., 23, 217 (1989).
c 2) Correlated sampling in quantum Monte Carlo: A route to forces,
c    Claudia Filippi and C. J. Umrigar, Phys. Rev. B., 61, R16291, (2000).

c stretch space so that electrons close to a nucleus move almost
c rigidly with that nucleus
c returns secondary geom. nuclear positions, stretched electron positions and
c jacobian of the transformation
      use constants_mod
      use nuclei_mod
      use atom_mod
      use const_mod
      use dim_mod
      use forcepar_mod
      use wfsec_mod
      use contr3_mod
      use force_dmc_mod
      implicit real*8(a-h,o-z)
      character*64 filename

      common /wfname/ filename

      dimension x(3,nelec),xstrech(3,nelec)
      dimension wt(ncent),dvol(3,3),dwt(3,ncent),dwtsm(3),cent_str(3,ncent)

!     JT: cent_ref, delc and pecentn used to be local, now there are in nuclei module

      save alfstr

      pecent=pecentn(ifr)
      do 1 icent=1,ncent
        do 1 k=1,ndim
    1     cent(k,icent)=cent_ref(k,icent)+delc(k,icent,ifr)

      ajacob=one

      if(istrech_el.eq.0) return

      do 2 i=1,nelec
        do 2 k=1,ndim
    2     xstrech(k,i)=x(k,i)

      if(istrech.eq.0) return


      do 50 i=1,nelec
        wtsm=zero
c initialize volume change matrix
        do 5 k=1,ndim
          dwtsm(k)=zero
          do 5 j=1,ndim
            dvol(j,k)=zero
            if(j.eq.k) dvol(j,k)=one
    5   continue

        do 20 icent=1,ncent
          dist2=zero
          do 10 k=1,ndim
   10       dist2=dist2+(x(k,i)-cent_ref(k,icent))**2
            dist=dsqrt(dist2)
            if(istrech.eq.1) then
              wt(icent)=dexp(-alfstr*dist)
            elseif(istrech.eq.2) then
              wt(icent)=one/dist**alfstr
            elseif(istrech.eq.3) then
              wt(icent)=dexp(alfstr/dist)
            else
              stop 'istrech must be 1, 2 or 3.'
            endif
            wtsm=wtsm+wt(icent)
            do 20 k=1,ndim
              if(istrech.eq.1) dwt(k,icent)=-alfstr*dexp(-alfstr*dist)*
     &        (x(k,i)-cent_ref(k,icent))/dist
              if(istrech.eq.2) dwt(k,icent)=-alfstr*
     &        (x(k,i)-cent_ref(k,icent))/dist**(alfstr+2)
              if(istrech.eq.3) dwt(k,icent)=-alfstr*dexp(alfstr/dist)*
     &        (x(k,i)-cent_ref(k,icent))/(dist2*dist)
              dwtsm(k)=dwtsm(k)+dwt(k,icent)
   20     continue
          wtsmi=one/wtsm

          do 40 icent=1,ncent
            do 40 k=1,ndim
              dwt(k,icent)=(wtsm*dwt(k,icent)-wt(icent)*dwtsm(k))/
     &        wtsm**2
              dvol(1,k)=dvol(1,k)+dwt(k,icent)*delc(1,icent,ifr)
              dvol(2,k)=dvol(2,k)+dwt(k,icent)*delc(2,icent,ifr)
   40         dvol(3,k)=dvol(3,k)+dwt(k,icent)*delc(3,icent,ifr)
          call matinv(dvol,3,det)
          ajacob=ajacob*det
          do 50 icent=1,ncent
            wt(icent)=wt(icent)*wtsmi
            do 50 k=1,ndim
              xstrech(k,i)=xstrech(k,i)+wt(icent)*delc(k,icent,ifr)
   50 continue

      return
c-----------------------------------------------------------------------

c read force parameters and set up n-n potential energy at displaced positions
      entry readforce

      open(3,file='input.force',status='old',form='formatted')

      read(3,*) istrech,alfstr
      write(6,'(''istrech,alfstr ='',i4,2f10.5)') istrech,alfstr

      call alloc ('delc', delc, 3, ncent, nwf)
      do 60 ifl=1,nforce
        do 60 ic=1,ncent
          read(3,*)  (delc(k,ic,ifl),k=1,ndim)
   60     write(6,'(''center '',i2,'' conf '',i2,'' displace '',
     &    3f10.2)') ic,ifl,(delc(k,ic,ifl),k=1,ndim)

      call alloc ('iwftype', iwftype, nforce)
      read(3,*) nwftype,(iwftype(i),i=1,nforce),filename
      
      if(nwftype.gt.nforce) stop 'nwftype gt nforce'
!JT      if(nwftype.gt.MWF) stop 'nwftype gt MWF'
      if(iwftype(1).ne.1) stop 'iwftype(1) ne 1'

      if(index(mode,'dmc').ne.0) then
        read(3,*) nwprod,itausec
        write(6,'(''nwprod,itausec='',2i4)') nwprod,itausec
        if(nwprod.lt.1) stop 'nwprod must be 1 or more'
        if(itausec.ne.0.and.itausec.ne.1) stop 'itausec must be 0 or 1'
        call object_modified ('nwprod')
      endif

      close(3)

      call alloc ('cent_ref', cent_ref, 3, ncent)
      do 65 k=1,ndim
        do 65 icent=1,ncent
   65     cent_ref(k,icent)=cent(k,icent)

c     do 80 ifl=1,nforce
c       pecentn(ifl)=zero
c       do 80 i=2,ncent
c         j1=i-1
c         do 80 j=1,j1
c           r2=zero
c           do 70 k=1,ndim
c  70         r2=r2+(cent(k,i)-cent(k,j)
c    &             +delc(k,i,ifl)-delc(k,j,ifl))**2
c         r=dsqrt(r2)
c  80     pecentn(ifl)=pecentn(ifl)+znuc(iwctype(i))*znuc(iwctype(j))/r

!      write(6,'(''pecentn(ifl)1'',9f9.4)') (pecentn(ifl),ifl=1,nforce)
      call alloc ('pecentn', pecentn, nwf)
      do 80 ifl=1,nforce
        do 70 i=1,ncent
            do 70 k=1,ndim
   70         cent_str(k,i)=cent(k,i)+delc(k,i,ifl)
   80     call pot_nn(cent_str,znuc,iwctype,ncent,pecentn(ifl))
      write(6,'(''pecentn(ifl)2'',9f9.4)') (pecentn(ifl),ifl=1,nforce)

      write(6,'(''n-n potential energies '',10f10.5)')(pecentn(ifl),ifl=1,nforce)

      call alloc ('deltot', deltot, nforce)
      do 110 ifl=1,nforce
        deltot(ifl)=zero
        rsq=zero
        do 100 jc=1,ncent
          do 100 k=1,ndim
            rcm=zero
            do 90 ic=1,ncent
              rcm=rcm+delc(k,ic,ifl)
   90         rsq=rsq+
     &        (cent(k,ic)+delc(k,ic,ifl)-cent(k,jc)-delc(k,jc,ifl))**2
            rcm=rcm/ncent
  100       deltot(ifl)=deltot(ifl)+(delc(k,jc,ifl)-rcm)**2
        if(ifl.eq.1) rsq1=rsq
c Warning: TEMPORARY: multiplication by ncent right for diatomics
        deltot(ifl)=sign(dsqrt(deltot(ifl)*ncent),rsq-rsq1)
c Warning: this is a temporary fix to put deltot=1 if one is doing excited
c states rather than forces
  110   if(deltot(ifl).eq.zero) deltot(ifl)=1

      return
c-----------------------------------------------------------------------

c For wf. optim. copy first geometry info to next two
      entry wf_copy2

      call alloc ('pecentn', pecentn, nwf)
      call alloc ('delc', delc, 3, ncent, nwf)
      do 160 ifl=1,nforce
        pecentn(ifl)=pecent
        do 160 ic=1,ncent
          do 160 k=1,ndim
  160       delc(k,ic,ifl)=0

      call alloc ('deltot', deltot, nforce)
      deltot (:) = 0.d0

       call alloc ('cent_ref', cent_ref, 3, ncent)
        do 165 ic=1,ncent
          do 165 k=1,ndim
  165       cent_ref(k,ic)=cent(k,ic)

c     do 170 iwf=2,nwftype
c 170   call basis_norm(iwf,1)


      return
      end

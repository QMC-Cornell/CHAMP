      subroutine cdeterminante(iel,x,rvec_en,r_en,ddet_det,determ)
c same subroutine as determinant() adapted to complex orbitals/determinants
c by A.D.Guclu Feb2004
c can deal with any complex determinant, provided that the determinantal
c coefficients are real.

      use control_mod
      use basic_tools_mod
      use cslater_mod

      implicit real*8(a-h,o-z)

!JT      parameter(one=1.d0)

c complex locals:
      complex*16 cddet_det(3,MELEC),cdeterm,cdetinv
      complex*16 corb(MORB)
      complex*16 cratio(MDET),csum,cterm

c complex commons:
c      complex*16 cslmui,cslmdi,cfpu,cfpd,cfppu,cfppd,cdetu,cdetd,cddeti_deti,cd2edeti_deti
c      complex*16 cdeti_det,cddeti_det,cd2deti_det,cd2det_det
      complex*16 cslmin,cdetn,cddeti_detin,cd2edeti_detin,cdorb,cddorb

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contrl_per/ iperiodic,ibasis
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
c      common /cslater/ cslmui(MMAT_DIM,MDET),cslmdi(MMAT_DIM,MDET)
c     &,cfpu(3,MMAT_DIM,MDET),cfpd(3,MMAT_DIM,MDET)
c     &,cfppu(MMAT_DIM,MDET),cfppd(MMAT_DIM,MDET)
c     &,cdetu(MDET),cdetd(MDET)
c     &,cddeti_deti(3,MELEC,MDET),cd2edeti_deti(MELEC,MDET),cdeti_det(MCSF),cddeti_det(3,MELEC,MCSF),cd2deti_det(MCSF),cd2det_det
      common /cslatn/ cslmin(MMAT_DIM,MDETUD),cdetn(MDETUD)
     &,cddeti_detin(3,MELEC,MDETUD),cd2edeti_detin(MELEC,MDETUD)
     &,cdorb(3,MORB),cddorb(MORB)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),ddet_det(3,*)

c allocate memory:
      n2=nelec*nelec
      call alloc('cslmui',cslmui,n2,ndetup)
      call alloc('cslmdi',cslmdi,n2,ndetdn)
      call alloc('cfpu',cfpu,ndim,n2,ndetup)
      call alloc('cfpd',cfpd,ndim,n2,ndetdn)
      call alloc('cdetu',cdetu,ndetup)
      call alloc('cdetd',cdetd,ndetdn)
      call alloc('cddeti_deti',cddeti_deti,ndim,nelec,ndet)
      call alloc('cd2edeti_deti',cd2edeti_deti,nelec,ndet)

      do 10 i=1,nelec
        do 10 idim=1,ndim
   10     cddet_det(idim,i)=dcmplx(0,0)
c      cdeterm=dcmplx(0,0)

c get orbitals and derivatives for all electron iel
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          call corbitals_loc_ana_grade(iel,rvec_en,r_en,corb,cdorb,cddorb)
         else
          stop 'complex calculations of numerical orbitals not implemented'
        endif

      else
          stop 'complex calculations of periodic systems not implemented'
      endif

      if(ipr.ge.4) then
        do 26 iorb=1,norb
          write(6,'(''iorb,corb='',i3,(30f9.5))') iorb,corb(iorb)
          write(6,'(''iorb,cdorb1='',i3,(30f9.5))') iorb,cdorb(1,iorb)
          write(6,'(''iorb,cdorb2='',i3,(30f9.5))') iorb,cdorb(2,iorb)
   26     write(6,'(''iorb,cddorb='',i3,(30f9.5))') iorb,cddorb(iorb)
      endif

      if(iel.le.nup) then

      ikel=nup*(iel-1)

      do 30 idet=1,ndetup
        cratio(idet)=dcmplx(0,0)
        do 30 j=1,nup
   30     cratio(idet)=cratio(idet)+cslmui(j+ikel,idet)*corb(iworbdup(j,idet))

      do 50 idet=1,ndetup
        cdetn(idet)=cdetu(idet)*cratio(idet)
        do 45 i=1,nup
          if(i.ne.iel) then
            ik=nup*(i-1)
            csum=dcmplx(0,0)
            do 35 j=1,nup
   35         csum=csum+cslmui(j+ik,idet)*corb(iworbdup(j,idet))
            csum=csum/cratio(idet)
            do 40 j=1,nup
   40         cslmin(j+ik,idet)=cslmui(j+ik,idet)-cslmui(j+ikel,idet)*csum
          endif
   45   continue
        do 50 j=1,nup
   50     cslmin(j+ikel,idet)=cslmui(j+ikel,idet)/cratio(idet)

      else

      ikel=ndn*(iel-nup-1)

      do 55 idet=1,ndetdn
        cratio(idet)=dcmplx(0,0)
        do 55 j=1,ndn
   55     cratio(idet)=cratio(idet)+cslmdi(j+ikel,idet)*corb(iworbddn(j,idet))

      do 75 idet=1,ndetdn
        cdetn(idet)=cdetd(idet)*cratio(idet)
        do 70 i=1,ndn
          if(i+nup.ne.iel) then
            ik=ndn*(i-1)
            csum=dcmplx(0,0)
            do 60 j=1,ndn
   60         csum=csum+cslmdi(j+ik,idet)*corb(iworbddn(j,idet))
            csum=csum/cratio(idet)
            do 65 j=1,ndn
   65        cslmin(j+ik,idet)=cslmdi(j+ik,idet)-cslmdi(j+ikel,idet)*csum
          endif
   70   continue
        do 75 j=1,ndn
   75     cslmin(j+ikel,idet)=cslmdi(j+ikel,idet)/cratio(idet)

      endif


        if(iel.le.nup) then
          do 85 idet=1,ndetup

c        cterm=cdetn(idet)*cdetd(idet)*cdet(idet,iwf)
          do 80 i=1,nup
            do 78 idim=1,ndim
   78         cddeti_detin(idim,i,idet)=dcmplx(0,0)
c what is the next line for?? :
   80         cd2edeti_detin(i,idet)=dcmplx(0,0)
          ik=-nup
          do 85 i=1,nup
            ik=ik+nup
            if(i.ne.iel) then
              do 82 j=1,nup
                do 82 idim=1,ndim
   82             cddeti_detin(idim,i,idet)=cddeti_detin(idim,i,idet)+cslmin(j+ik,idet)*cfpu(idim,j+ik,idet)
             else
               do 84 j=1,nup
                 do 84 idim=1,ndim
   84              cddeti_detin(idim,i,idet)=cddeti_detin(idim,i,idet)+cslmin(j+ik,idet)*cdorb(idim,iworbdup(j,idet))
            endif
   85       continue

         else

       do 95 idet=1,ndetdn
c        cterm=cdetu(idet)*cdetn(idet)*cdet(idet,iwf)
          do 90 i=nup+1,nelec
            do 90 idim=1,ndim
   90         cddeti_detin(idim,i,idet)=dcmplx(0,0)
          ik=-ndn
          do 95 i=nup+1,nelec
            ik=ik+ndn
            if(i.ne.iel) then
              do 92 j=1,ndn
                do 92 idim=1,ndim
   92             cddeti_detin(idim,i,idet)=cddeti_detin(idim,i,idet)+cslmin(j+ik,idet)*cfpd(idim,j+ik,idet)
             else
              do 94 j=1,ndn
                do 94 idim=1,ndim
   94             cddeti_detin(idim,i,idet)=cddeti_detin(idim,i,idet)+cslmin(j+ik,idet)*cdorb(idim,iworbddn(j,idet))
             endif
   95     continue

      endif

      cdeterm=dcmplx(0,0)
      do 115 icsf=1,ncsf
        do 115 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
          if(iel.le.nup) then
            if(ndn.ge.1) then
              cterm=cdetn(iwdetup(idet))*cdetd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
             else
              cterm=cdetn(iwdetup(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
            endif
            do 106 i=1,nup
              iwdet=iwdetup(idet)
              do 106 idim=1,ndim
  106           cddet_det(idim,i)=cddet_det(idim,i)+cddeti_detin(idim,i,iwdet)*cterm
            do 107 i=nup+1,nelec
              iwdet=iwdetdn(idet)
              do 107 idim=1,ndim
  107           cddet_det(idim,i)=cddet_det(idim,i)+cddeti_deti(idim,i,iwdet)*cterm
           else
            cterm=cdetu(iwdetup(idet))*cdetn(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
            do 108 i=1,nup
              iwdet=iwdetup(idet)
              do 108 idim=1,ndim
  108           cddet_det(idim,i)=cddet_det(idim,i)+cddeti_deti(idim,i,iwdet)*cterm
            do 109 i=nup+1,nelec
              iwdet=iwdetdn(idet)
              do 109 idim=1,ndim
  109           cddet_det(idim,i)=cddet_det(idim,i)+cddeti_detin(idim,i,iwdet)*cterm
          endif
  115     cdeterm=cdeterm+cterm

      cdetinv=one/cdeterm
      determ=cdabs(cdeterm)
      do 120 i=1,nelec
        do 120 idim=1,ndim
          cddet_det(idim,i)=cddet_det(idim,i)*cdetinv
  120     ddet_det(idim,i)=dreal(cddet_det(idim,i))
      return
      end

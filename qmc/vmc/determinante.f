      subroutine determinante(iel,x,rvec_en,r_en,ddet_det,determ)
c Written by Claudia Filippi by modifying determinant, modified by Cyrus Umrigar
      use control_mod
      use dorb_mod
      use orbitals_mod, only: orb_tot_nb
      use slatn_mod
      use orbe_mod
      use dets_mod
      use slater_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'force.h'

!JT      parameter(one=1.d0)

c     common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contrl_per/ iperiodic,ibasis
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
c     common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
c    &,d2phin(MBASIS,MELEC)
c     common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /kinet/ ekineo(MELEC),ekinen(MELEC)
!JT      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
!JT     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
!JT      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
!JT     &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
!JT     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
!JT     &,detu(MDETUD),detd(MDETUD)
!JT     &,ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
!JT     &,detij_det(MPARMD,MPARMD)
!JT      common /slatn/ slmin(MMAT_DIM,MDETUD),detn(MDETUD)
!JT     &,ddeti_detin(3,MELEC,MDETUD),d2edeti_detin(MELEC,MDETUD)
!JT     &,dorb(3,MORB),ddorb(MORB)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
!JT      dimension orbe(orb_tot_nb),ddet_det(3,*),ratio(MDET)
      dimension ddet_det(3,*),ratio(MDET)

      do 10 i=1,nelec
        ddet_det(1,i)=0
        ddet_det(2,i)=0
   10   ddet_det(3,i)=0
c     determ=0

c get orbitals and derivatives for all electron iel
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          call orbitals_loc_ana_grade(iel,rvec_en,r_en,orbe,dorbe,ddorbe)
         else
          call orbitals_loc_num_grade(iel,x,orbe,dorbe,ddorbe)
        endif

       else

        if(inum_orb.eq.0) then
          call orbitals_pw_grade(x(1,iel),orbe,dorbe,ddorbe)
         else
          call orbitals_period_num_grade(x(1,iel),orbe,dorbe,ddorbe)
        endif

      endif


      if(iel.le.nup) then

        ikel=nup*(iel-1)

        do 30 idet=1,ndetup
          ratio(idet)=0
          do 30 j=1,nup
   30       ratio(idet)=ratio(idet)+slmui(j+ikel,idet)*orbe(iworbdup(j,idet))

        do 50 idet=1,ndetup
          detn(idet)=detu(idet)*ratio(idet)
          do 45 i=1,nup
            if(i.ne.iel) then
              ik=nup*(i-1)
              sum=0
              do 35 j=1,nup
   35           sum=sum+slmui(j+ik,idet)*orbe(iworbdup(j,idet))
              sum=sum/ratio(idet)
              do 40 j=1,nup
   40           slmin(j+ik,idet)=slmui(j+ik,idet)-slmui(j+ikel,idet)*sum
            endif
   45     continue
          do 50 j=1,nup
   50       slmin(j+ikel,idet)=slmui(j+ikel,idet)/ratio(idet)

       else

        ikel=ndn*(iel-nup-1)

        do 55 idet=1,ndetdn
          ratio(idet)=0
          do 55 j=1,ndn
   55       ratio(idet)=ratio(idet)+slmdi(j+ikel,idet)*orbe(iworbddn(j,idet))

        do 75 idet=1,ndetdn
          detn(idet)=detd(idet)*ratio(idet)
          do 70 i=1,ndn
            if(i+nup.ne.iel) then
              ik=ndn*(i-1)
              sum=0
              do 60 j=1,ndn
   60           sum=sum+slmdi(j+ik,idet)*orbe(iworbddn(j,idet))
              sum=sum/ratio(idet)
              do 65 j=1,ndn
   65           slmin(j+ik,idet)=slmdi(j+ik,idet)-slmdi(j+ikel,idet)*sum
            endif
   70     continue
          do 75 j=1,ndn
   75       slmin(j+ikel,idet)=slmdi(j+ikel,idet)/ratio(idet)

      endif


      if(iel.le.nup) then
        do 85 idet=1,ndetup

c         term=detn(idet)*detd(idet)*cdet(idet,iwf)
          do 80 i=1,nup
            ddeti_detin(1,i,idet)=0
            ddeti_detin(2,i,idet)=0
            ddeti_detin(3,i,idet)=0
   80       d2edeti_detin(i,idet)=0
          ik=-nup
          do 85 i=1,nup
            ik=ik+nup
            if(i.ne.iel) then
              do 82 j=1,nup
                ddeti_detin(1,i,idet)=ddeti_detin(1,i,idet)+slmin(j+ik,idet)*fpu(1,j+ik,idet)
                ddeti_detin(2,i,idet)=ddeti_detin(2,i,idet)+slmin(j+ik,idet)*fpu(2,j+ik,idet)
   82           ddeti_detin(3,i,idet)=ddeti_detin(3,i,idet)+slmin(j+ik,idet)*fpu(3,j+ik,idet)
             else
              do 84 j=1,nup
                ddeti_detin(1,i,idet)=ddeti_detin(1,i,idet)+slmin(j+ik,idet)*dorbe(1,iworbdup(j,idet))
                ddeti_detin(2,i,idet)=ddeti_detin(2,i,idet)+slmin(j+ik,idet)*dorbe(2,iworbdup(j,idet))
   84           ddeti_detin(3,i,idet)=ddeti_detin(3,i,idet)+slmin(j+ik,idet)*dorbe(3,iworbdup(j,idet))
            endif
   85     continue

         else

       do 95 idet=1,ndetdn
c         term=detu(idet)*detn(idet)*cdet(idet,iwf)
          do 90 i=nup+1,nelec
            ddeti_detin(1,i,idet)=0
            ddeti_detin(2,i,idet)=0
   90       ddeti_detin(3,i,idet)=0
          ik=-ndn
          do 95 i=nup+1,nelec
            ik=ik+ndn
            if(i.ne.iel) then
              do 92 j=1,ndn
                ddeti_detin(1,i,idet)=ddeti_detin(1,i,idet)+slmin(j+ik,idet)*fpd(1,j+ik,idet)
                ddeti_detin(2,i,idet)=ddeti_detin(2,i,idet)+slmin(j+ik,idet)*fpd(2,j+ik,idet)
   92           ddeti_detin(3,i,idet)=ddeti_detin(3,i,idet)+slmin(j+ik,idet)*fpd(3,j+ik,idet)
             else
              do 94 j=1,ndn
                ddeti_detin(1,i,idet)=ddeti_detin(1,i,idet)+slmin(j+ik,idet)*dorbe(1,iworbddn(j,idet))
                ddeti_detin(2,i,idet)=ddeti_detin(2,i,idet)+slmin(j+ik,idet)*dorbe(2,iworbddn(j,idet))
   94           ddeti_detin(3,i,idet)=ddeti_detin(3,i,idet)+slmin(j+ik,idet)*dorbe(3,iworbddn(j,idet))
            endif
   95     continue

      endif

      determ=0
      do 115 icsf=1,ncsf
        do 115 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
          if(iel.le.nup) then
            if(ndn.ge.1) then
              term=detn(iwdetup(idet))*detd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
             else
              term=detn(iwdetup(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
            endif
            do 106 i=1,nup
              iwdet=iwdetup(idet)
              ddet_det(1,i)=ddet_det(1,i)+ddeti_detin(1,i,iwdet)*term
              ddet_det(2,i)=ddet_det(2,i)+ddeti_detin(2,i,iwdet)*term
  106         ddet_det(3,i)=ddet_det(3,i)+ddeti_detin(3,i,iwdet)*term
            do 107 i=nup+1,nelec
              iwdet=iwdetdn(idet)
              ddet_det(1,i)=ddet_det(1,i)+ddeti_deti(1,i,iwdet)*term
              ddet_det(2,i)=ddet_det(2,i)+ddeti_deti(2,i,iwdet)*term
  107         ddet_det(3,i)=ddet_det(3,i)+ddeti_deti(3,i,iwdet)*term
           else
            term=detu(iwdetup(idet))*detn(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
            do 108 i=1,nup
              iwdet=iwdetup(idet)
              ddet_det(1,i)=ddet_det(1,i)+ddeti_deti(1,i,iwdet)*term
              ddet_det(2,i)=ddet_det(2,i)+ddeti_deti(2,i,iwdet)*term
  108         ddet_det(3,i)=ddet_det(3,i)+ddeti_deti(3,i,iwdet)*term
            do 109 i=nup+1,nelec
              iwdet=iwdetdn(idet)
              ddet_det(1,i)=ddet_det(1,i)+ddeti_detin(1,i,iwdet)*term
              ddet_det(2,i)=ddet_det(2,i)+ddeti_detin(2,i,iwdet)*term
  109         ddet_det(3,i)=ddet_det(3,i)+ddeti_detin(3,i,iwdet)*term
          endif
  115     determ=determ+term

      detinv=one/determ
      do 120 i=1,nelec
        ddet_det(1,i)=ddet_det(1,i)*detinv
        ddet_det(2,i)=ddet_det(2,i)*detinv
  120   ddet_det(3,i)=ddet_det(3,i)*detinv

      return
      end

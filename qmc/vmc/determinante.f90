      subroutine determinante(iel,x,rvec_en,r_en,ddet_det,determ)
! Written by Claudia Filippi by modifying determinant, modified by Cyrus Umrigar
      use constants_mod
      use control_mod
      use dorb_mod
      use slatn_mod
      use orbe_mod
      use dets_mod
      use slater_mod
      use const_mod
      use contr2_mod
      use wfsec_mod
      use contrl_per_mod
      implicit real*8(a-h,o-z)


      dimension x(3,*),rvec_en(3,nelec,*),r_en(nelec,*)
      dimension ddet_det(3,*),ratio(ndet)

      do 10 i=1,nelec
        ddet_det(1,i)=0
        ddet_det(2,i)=0
   10   ddet_det(3,i)=0
!     determ=0

! get orbitals and derivatives for all electron iel
      if(iperiodic.eq.0 .or. iperiodic.eq.1) then

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


!     write(6,'(''iel,orbe='',i3,(30f9.5))') iel,(orbe(iorb),iorb=1,nelec)

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

!         term=detn(idet)*detd(idet)*cdet(idet,iwf)
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
!         term=detu(idet)*detn(idet)*cdet(idet,iwf)
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

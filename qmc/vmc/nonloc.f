      subroutine nonloc(x,rshift,rvec_en,r_en,detu,detd,slmui,slmdi,vpsp)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use control_mod
      use deriv_orb_mod
      use periodic_jastrow_mod !WAS
      use atom_mod
      use dets_mod

      use const_mod
      use dim_mod
      use pseudo_mod
      use contrl_per_mod
      use periodic_mod
      use qua_mod
      implicit real*8(a-h,o-z)
!JT      parameter (one=1.d0)


!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /contrl_per/ iperiodic,ibasis

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
!JT     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
!JT      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
!JT     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad

!JT      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
!JT     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
!JT     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
!JT     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
!JT     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
!JT     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
!JT     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
!JT     &,znuc_sum,znuc2_sum,vcell,vcell_sim
!JT     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
!JT     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
!JT     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

      dimension x(3,*),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)
      dimension rr_en(MELEC,MCENT),rr_en2(MELEC,MCENT),rr_en_sav(MCENT),rr_en2_sav(MCENT)
     &,xsav(3),rshift_sav(3,MCENT),rvec_en_sav(3,MCENT),r_en_sav(MCENT),vpot(MPS_L)


      do 10 i=1,nelec
        do 10 ic=1,ncent
          call scale_dist(r_en(i,ic),rr_en(i,ic),1)
   10     call scale_dist(r_en(i,ic),rr_en2(i,ic),3)

c     write(6,'(''x='',30f9.4)') ((x(k,i),k=1,ndim),i=1,nelec)
c     write(6,'(''r_en='',30f9.4)') ((r_en(i,ic),i=1,nelec),ic=1,ncent)
c     write(6,'(''rvec_en='',60f9.4)') (((rvec_en(k,i,ic),k=1,ndim),i=1,nelec),ic=1,ncent)


! JT beg
      if (l_opt_orb_energy) then
       call object_provide_by_index (param_orb_nb_index)
       call object_alloc ('vpot_ex', vpot_ex, MPS_L, param_orb_nb)
       call object_alloc ('vpsp_ex', vpsp_ex, param_orb_nb)
       vpsp_ex = 0.d0
      endif
! JT end

      vpsp=0
      do 100 ic=1,ncent

        ict=iwctype(ic)
        do 100 i=1,nelec

c vps was calculated by calling getvps_tm from nonloc_pot
          iskip=1
          do 15 l=1,npotd(ict)
   15       if(l.ne.lpotp1(ict) .and. dabs(vps(i,ic,l)).gt.1.d-4) iskip=0

          if(iskip.eq.0) then

            ri=one/r_en(i,ic)

            do 20 l=1,npotd(ict)
   20         if(l.ne.lpotp1(ict)) vpot(l)=0

            if (l_opt_orb_energy) then  !JT
               vpot_ex = 0.d0           !JT
            endif                       !JT

            do 28 k=1,ndim
   28         xsav(k)=x(k,i)
            do 30 jc=1,ncent
              r_en_sav(jc)=r_en(i,jc)
              rr_en_sav(jc)=rr_en(i,jc)
              rr_en2_sav(jc)=rr_en2(i,jc)
              do 30 k=1,ndim
                rshift_sav(k,jc)=rshift(k,i,jc)
   30           rvec_en_sav(k,jc)=rvec_en(k,i,jc)

            do 60 iq=1,nquad
              costh=rvec_en_sav(1,ic)*xq(iq)
     &             +rvec_en_sav(2,ic)*yq(iq)
     &             +rvec_en_sav(3,ic)*zq(iq)
              costh=costh*ri

              if(iperiodic.eq.0) then
                x(1,i)=r_en(i,ic)*xq(iq)+cent(1,ic)
                x(2,i)=r_en(i,ic)*yq(iq)+cent(2,ic)
                x(3,i)=r_en(i,ic)*zq(iq)+cent(3,ic)
               else
                x(1,i)=r_en(i,ic)*xq(iq)+cent(1,ic)+rshift(1,i,ic)
                x(2,i)=r_en(i,ic)*yq(iq)+cent(2,ic)+rshift(2,i,ic)
                x(3,i)=r_en(i,ic)*zq(iq)+cent(3,ic)+rshift(3,i,ic)
              endif

              do 40 jc=1,ncent
                do 38 k=1,ndim
   38             rvec_en(k,i,jc)=x(k,i)-cent(k,jc)

                if(jc.ne.ic) then
                  if(iperiodic.eq.0) then
                    r_en(i,jc)=0
                    do 39 k=1,ndim
   39                 r_en(i,jc)=r_en(i,jc)+rvec_en(k,i,jc)**2
                    r_en(i,jc)=dsqrt(r_en(i,jc))
                   else
                    call find_image4(rshift(1,i,jc),rvec_en(1,i,jc),r_en(i,jc),rlatt,rlatt_inv)
                  endif

                  call scale_dist(r_en(i,jc),rr_en(i,jc),1)
                  call scale_dist(r_en(i,jc),rr_en2(i,jc),3)
                endif
   40         continue

              iel=i

              electron = iel !JT
              call object_modified_by_index (electron_index) !JT

              call nonlocd(iel,x(1,i),rvec_en,r_en,detu,detd,slmui,slmdi,deter)
c             call nonlocd(iel,x,rvec_en,r_en,detu,detd,slmui,slmdi,deter)

cWAS              call nonlocj(iel,x,rshift,rr_en,rr_en2,value)

              call nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,value)

!WAS
              if (do_pjas) then
                 call nonloc_pjas (iel, x(:,1:melec), value)
              endif
!WAS

              if(ipr.ge.4) then
                write(6,'(''rr_en,rr_en2'',2d14.6)') rr_en(1,1),rr_en2(1,1)
                write(6,'(''ic,i,iq,deter,value'',3i3,2d14.6)') ic,i,iq,deter,value
              endif

              do 50 l=1,npotd(ict)
                if(l.ne.lpotp1(ict)) then
                  vpot(l)=vpot(l)+wq(iq)*yl0(l,costh)*deter*exp(value)

! JT beg
!             For singly-excited wave functions
              if (l_opt_orb_energy) then
                 call object_provide_by_index (psid_ex_in_x_index)
                 do iex = 1, param_orb_nb
                  vpot_ex(l,iex)=vpot_ex(l,iex)+wq(iq)*yl0(l,costh)*psid_ex_in_x(iex)*exp(value)
                 enddo
              endif
! JT end

                endif
   50         continue

   60       continue

            do 68 k=1,ndim
   68         x(k,i)=xsav(k)
            do 70 jc=1,ncent
              r_en(i,jc)=r_en_sav(jc)
              rr_en(i,jc)=rr_en_sav(jc)
              rr_en2(i,jc)=rr_en2_sav(jc)
              do 70 k=1,ndim
                rshift(k,i,jc)=rshift_sav(k,jc)
   70           rvec_en(k,i,jc)=rvec_en_sav(k,jc)

            do 80 l=1,npotd(ict)
              if(l.ne.lpotp1(ict)) then
                vpsp=vpsp+vps(i,ic,l)*vpot(l)
                if(ipr.ge.4) write(6,'(''nonloc: i,ic,l,vps(i,ic,l),vpot(l),vpsp'',3i5,9d12.4)')
     &          i,ic,l,vps(i,ic,l),vpot(l),vpsp

! JT beg
!             For singly-excited wave functions
              if (l_opt_orb_energy) then
                 do iex = 1, param_orb_nb
                  vpsp_ex(iex)=vpsp_ex(iex)+vps(i,ic,l)*vpot_ex(l,iex)
                 enddo
              endif
! JT end


              endif
   80       continue

          endif
  100 continue

      call object_modified_by_index (vpsp_ex_index) ! JT


c     write(6,'(''x='',30f9.4)') ((x(k,i),k=1,ndim),i=1,nelec)
c     write(6,'(''r_en='',30f9.4)') ((r_en(i,ic),i=1,nelec),ic=1,ncent)
c     write(6,'(''rvec_en='',60f9.4)') (((rvec_en(k,i,ic),k=1,ndim),i=1,nelec),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------

      function yl0(l,costh)
c (2L+1)*P_L(costh)
c This is not quite Y_L0 but sqrt(4pi/(2L+1)) Y_L0
c Note that the associated P_L^m and the unassociated P_L Legendre polynomials are the same for m=0.
c l is actually L+1.

      implicit real*8(a-h,o-z)

      if(l.eq.1) then
        yl0=1.d0
       elseif(l.eq.2) then
        yl0=3.d0*costh
       elseif(l.eq.3) then
        yl0=2.5d0*(3*costh*costh-1)
       elseif(l.eq.4) then
        yl0=3.5d0*costh*(5*costh*costh-3)
       elseif(l.eq.5) then
        yl0=1.125d0*(35*costh**4-30*costh**2+3)
       else
        stop 'yl0 implemented to l=4 only (Warning: l is l+1)'
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine nonlocd(iel,x,rvec_en,r_en,detu,detd,slmui,slmdi,determ)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use all_tools_mod
      use control_mod
      use eloc_mod
      use dorb_mod
      use slatn_mod
      use orbe_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use contr2_mod
      use contrl_opt2_mod
      use wfsec_mod
      use contrl_per_mod
      use contr3_mod
      use phifun_mod
      implicit real*8(a-h,o-z)


c     common /dim/ ndim
!JT      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
!JT     &,ifock,i3body,irewgt,iaver,istrch
!JT     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
!JT      common /contr3/ mode
!JT      common /contrl_per/ iperiodic,ibasis
!JT      common /contrl_opt2/ igradhess,iadd_diag_opt
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /slatn2/ deti_new(MPARMD)
!JT      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
!JT     &,d2phin(MBASIS,MELEC)
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
!JT      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
!JT     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
!JT      common /wfsec/ iwftype(MFORCE),iwf,nwftype
!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

!JT     common /orbe/ orbe(MORB),detn(MDETUD) !JT
      dimension x(3),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)
      dimension ratio(MDET)

c     determ=0

c get orbitals for all electron iel
c Note x is 3-dim but rvec_en,r_en have info about all electrons.
c So, iel is passed to select elements of rvec_en,r_en,phin and for IO.
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          call orbitals_loc_anae(iel,rvec_en,r_en,orbe)
         else
          call orbitals_loc_nume(x,orbe)
        endif

       else

        if(inum_orb.eq.0) then
          call orbitals_pwe(iel,x,orbe)
         else
          call orbitals_period_nume(x,orbe)
        endif

      endif

      call object_modified_by_index (orbe_index) !JT

      if(iel.le.nup) then

        ikel=nup*(iel-1)

        do 30 idet=1,ndetup
          ratio(idet)=0
          do 30 j=1,nup
c           write(6,'(''slmui(j+ikel,idet)*orbe(iworbdup(j,idet))'',9d12.4)') slmui(j+ikel,idet),orbe(iworbdup(j,idet))
   30       ratio(idet)=ratio(idet)+slmui(j+ikel,idet)*orbe(iworbdup(j,idet))

        do 50 idet=1,ndetup
c         write(6,'(''detu(idet),ratio(idet)'',9d12.4)') detu(idet),ratio(idet)
   50     detn(idet)=detu(idet)*ratio(idet)

       else

        ikel=ndn*(iel-nup-1)

        do 55 idet=1,ndetdn
          ratio(idet)=0
          do 55 j=1,ndn
   55       ratio(idet)=ratio(idet)+slmdi(j+ikel,idet)*orbe(iworbddn(j,idet))

        do 75 idet=1,ndetdn
   75     detn(idet)=detd(idet)*ratio(idet)

      endif

      call object_modified_by_index (detn_index) !JT

      determ=0
      do 115 icsf=1,ncsf
        do 115 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
          if(iel.le.nup) then
            term=detn(iwdetup(idet))*detd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
           else
            term=detu(iwdetup(idet))*detn(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
          endif
  115     determ=determ+term

c Derivatives wrt to csf_coefs for optimizing them
      if(index(mode,'fit').ne.0 .or. igradhess.ge.1 .or. l_opt_csf) then
        do 140 iparm=1,nparmcsf
          icsf=iwcsf(iparm)
          deti_new(iparm)=0
          do 140 idet_in_csf=1,ndet_in_csf(icsf)
            idet=iwdet_in_csf(idet_in_csf,icsf)
            if(iel.le.nup) then
              term=detn(iwdetup(idet))*detd(iwdetdn(idet))*cdet_in_csf(idet_in_csf,icsf)
             else
              term=detu(iwdetup(idet))*detn(iwdetdn(idet))*cdet_in_csf(idet_in_csf,icsf)
            endif
  140       deti_new(iparm)=deti_new(iparm)+term
      endif

      return
      end
c-----------------------------------------------------------------------

!      subroutine nonlocj(iel,x,rshift,rr_en,rr_en2,value)
!WAS
      subroutine nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,value)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use control_mod
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use contrl_per_mod
      use jaspar_mod
      use bparm_mod
      use periodic_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'ewald.h'
!JT      include 'force.h'

!JT      parameter (half=.5d0)

!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /contrl_per/ iperiodic,ibasis
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
!JT      common /jaspar/ nspin1,nspin2,sspin,sspinn,is

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /bparm/ nspin2b,nocuspb
!JT      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
!JT     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
!JT     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
!JT     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
!JT     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
!JT     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
!JT     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
!JT     &,znuc_sum,znuc2_sum,vcell,vcell_sim
!JT     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
!JT     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
!JT     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

      common /jaso/ fso(MELEC,MELEC),fij(3,MELEC,MELEC)
     &,d2ij(MELEC,MELEC),d2,fsum,fjo(3,MELEC)

      dimension x(3,*),rshift(3,MELEC,MCENT),rr_en(MELEC,MCENT),rr_en2(MELEC,MCENT)
     &,fsn(MELEC,MELEC),dx(3)

!WAS
      dimension r_en(MELEC,MCENT)
!!

      fsumn=0

      if(nelec.lt.2) goto 47

      do 45 jj=1,nelec

        if(jj.eq.iel) goto 45
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif

        sspinn=1
        ipar=0
        isb=1
        if(i.le.nup .or. j.gt.nup) then
          if(nspin2b.eq.2) then
            isb=2
           elseif(nocuspb.eq.0) then
            sspinn=half
          endif
          ipar=1
        endif

        do 10 k=1,ndim
   10     dx(k)=x(k,jj)-x(k,iel)

        if(iperiodic.eq.0) then
          rij=0
          do 20 k=1,ndim
   20       rij=rij+dx(k)**2
          rij=dsqrt(rij)
         else
          call find_image3(dx,rij,rlatt_sim,rlatt_sim_inv)
        endif

c e-e terms
        call scale_dist(rij,u,2)

        fsn(i,j)=psibnl(u,isb,ipar)

c e-e-n terms
c The scaling is switched in psinl, so do not do it here.
c     if(isc.ge.12) call scale_dist(rij,u,3)
      call scale_dist(rij,u,4)

        do 40 ic=1,ncent
          it=iwctype(ic)
   40     fsn(i,j)=fsn(i,j) +
!!
!!     &    psinl(u,rshift(1,i,ic),rshift(1,j,ic),rr_en2(i,ic),rr_en2(j,ic),it)
!WAS
     &    psinl(u,rshift(1,i,ic),rshift(1,j,ic),r_en(i,ic),r_en(j,ic),
     &         rr_en2(i,ic),rr_en2(j,ic),it)
!!

        fsumn=fsumn+fsn(i,j)-fso(i,j)

   45 continue

c e-n terms
   47 fsn(iel,iel)=0
      do 50 ic=1,ncent
        it=iwctype(ic)
   50   fsn(iel,iel)=fsn(iel,iel)+psianl(rr_en(iel,ic),it)

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)

      value=fsumn

      return
      end

      subroutine grad_hess_jas_sum(p,q,enew,eold,wt,wi_w)
c Note that if it is called from a 1-electron move routine then p=1 and only new values are relevant.
c Hence it is not necessary to call grad_hess_jas_save.
      use mpi_mod
      use atom_mod
      use slater_mod
      use optim_mod
      use const_mod
      implicit real*8(a-h,o-z)

      parameter(factor_max=1.d2,ratio_max=1.1d0)

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
!JT    &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
!JT     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
!JT     &,detu(MDETUD),detd(MDETUD)
!JT     &,ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
!JT     &,detij_det(MPARMD,MPARMD)
      common /bparm/ nspin2b,nocuspb

      common /delocc/ denergy(MPARM)
      common /derivjas/ gvalue(MPARMJ),g(3,MELEC,MPARMJ),d2g(MPARMJ)
     &,go(MELEC,MELEC,MPARMJ)

      common /gradhessder/ dj(MPARM),dj_e(MPARM),dj_de(MPARM,MPARM),dj_dj(MPARM,MPARM),dj_dj_e(MPARM,MPARM)
     &,de(MPARM),d2j(MPARM,MPARM),d2j_e(MPARM,MPARM),de_e(MPARM),e2(MPARM),dj_e2(MPARM),de_de(MPARM,MPARM)
     &,w_i(MPARM),w_i_e(MPARM)
      common /gradhessdero/ deti_det_old(MPARMD),gvalue_old(MPARMJ),denergy_old(MPARM)
     &,d1d2a_old(MCTYPE),d2d2a_old(MCTYPE),d1d2b_old(2),d2d2b_old(2),didk_old(MPARMJ)
     &,detij_det_old(MPARMD,MPARMD)

      common /gradhessj_nonlin/ d1d2a(MCTYPE),d2d2a(MCTYPE),d1d2b(2),d2d2b(2),didk(MPARMJ)

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype
      common /pointer/ npointa(MPARMJ*NCTYP3X)

      common /contrl_opt2/ igradhess,iadd_diag_opt

      dimension wi_w(MPARM)

      save wt_tot

      data wt_tot /0.d0/

      if(igradhess.eq.0) return

      wt_tot=wt_tot+wt

c Warning: temp
c     ratio=1
c     ratio_old=1
c     do 5 iparm=1,nparmcsf
cc      write(6,'(''iparm,deti_det(iparm),deti_det_old(iparm)'',i3,9d12.4)') iparm,deti_det(iparm),deti_det_old(iparm)
c       ratio=max(ratio,abs(deti_det(iparm)))
c   5   ratio_old=max(ratio_old,abs(deti_det_old(iparm)))
c     ratio_lim=(ratio_max*ratio-1)/(ratio_max+ratio-2)
c     ratio_old_lim=(ratio_max*ratio_old-1)/(ratio_max+ratio_old-2)
cc    ratio_inv=1/ratio
cc    ratio_old_inv=1/ratio_old
c     ratio_inv=ratio_lim/ratio
c     ratio_old_inv=ratio_old_lim/ratio_old
cc    write(6,'(''ratio,ratio_old'',9d12.4)') ratio,ratio_old

      ratio_inv=1
      ratio_old_inv=1

c For dj_de compute both the (iparm,jparm) and the (jparm,iparm) elements making sure not to do the diagonal twice.
      do 10 iparm=1,nparmcsf+nparmot+nparmj+nparms
        if(iparm.le.nparmcsf+nparmot) then
c Warning next 2 lines tmp
c         ratio_inv=min(1.d0,ratio_max/deti_det(iparm))
c         ratio_old_inv=min(1.d0,ratio_max/deti_det_old(iparm))

          psii_psi=deti_det(iparm)*ratio_inv
          psii_psi_old=deti_det_old(iparm)*ratio_old_inv
          denergyi=denergy(iparm)*ratio_inv
          denergyi_old=denergy_old(iparm)*ratio_old_inv
c Warning temp printout, can create huge output
c         if(abs(psii_psi).gt.factor_max .or. abs(denergy(iparm)).gt.10.d0*factor_max) then
cc          factori=factor_max
c           write(6,'(''iparm,psii_psi,denergy(iparm),enew'',i3,9d12.4)')
c    &      iparm,psii_psi,denergy(iparm),enew
c         endif
c         if(abs(psii_psi_old).gt.factor_max) then
c           write(6,'(''iparm,psii_psi_old,denergy_old(iparm),eold'',i3,9d12.4)')
c    &      iparm,psii_psi_old,denergy_old(iparm),eold
c         endif
         else
          i=iparm-nparmcsf-nparmot
          psii_psi=gvalue(i)
          psii_psi_old=gvalue_old(i)
          denergyi=denergy(iparm)
          denergyi_old=denergy_old(iparm)
        endif
        dj(iparm)=dj(iparm)      +wt*(p*psii_psi+q*psii_psi_old)
        dj_e(iparm)=dj_e(iparm)  +wt*(p*psii_psi*enew+q*psii_psi_old*eold)
        dj_e2(iparm)=dj_e2(iparm)+wt*(p*psii_psi*enew**2+q*psii_psi_old*eold**2)
        de(iparm)=de(iparm)      +wt*(p*denergyi+q*denergyi_old)
        de_e(iparm)=de_e(iparm)  +wt*(p*denergyi*enew+q*denergyi_old*eold)
        e2(iparm)=e2(iparm)      +wt*(p*enew**2+q*eold**2)
        w_i(iparm)=w_i(iparm)    +wt*wi_w(iparm)
        w_i_e(iparm)=w_i_e(iparm)+wt*wi_w(iparm)*(p*enew+q*eold)
        if(ipr.ge.1) write(6,'(''iparm,denergyi,enew,eold,='',i5,9d12.4)') iparm,denergyi,enew,eold
        if(ipr.ge.1) write(6,'(''iparm,dj(iparm),dj_e(iparm),dj_e2(iparm),de(iparm),de_e(iparm),e2(iparm)'',i5,9d12.4)')
     &  iparm,dj(iparm),dj_e(iparm),dj_e2(iparm),de(iparm),de_e(iparm),e2(iparm)
        do 10 jparm=1,iparm
          if(jparm.le.nparmcsf+nparmot) then
            psij_psi=deti_det(jparm)*ratio_inv
            psij_psi_old=deti_det_old(jparm)*ratio_old_inv
            denergyj=denergy(jparm)*ratio_inv
            denergyj_old=denergy_old(jparm)*ratio_old_inv
           else
            j=jparm-nparmcsf-nparmot
            psij_psi=gvalue(j)
            psij_psi_old=gvalue_old(j)
            denergyj=denergy(jparm)
            denergyj_old=denergy_old(jparm)
          endif
          dj_dj(iparm,jparm)=dj_dj(iparm,jparm)+wt*(p*psii_psi* psij_psi+q*psii_psi_old*psij_psi_old)
          dj_de(iparm,jparm)=dj_de(iparm,jparm)+wt*(p*psii_psi*denergyj+q*psii_psi_old*denergyj_old)
          if(jparm.lt.iparm) dj_de(jparm,iparm)=dj_de(jparm,iparm)+wt*(p*psij_psi*denergyi+q*psij_psi_old*denergyi_old)
          dj_dj_e(iparm,jparm)=dj_dj_e(iparm,jparm)+wt*(p*psii_psi* psij_psi*enew+q*psii_psi_old*psij_psi_old*eold)
          de_de(iparm,jparm)=de_de(iparm,jparm)+wt*(p*denergyi*denergyj+q*denergyi_old*denergyj_old)
      if(ipr.ge.1)
     &write(6,'(''i,j,dj_dj(iparm,jparm),dj_de(iparm,jparm),dj_de(jparm,iparm),dj_dj_e(iparm,jparm),de_de(iparm,jparm)'',2i5,9d12.4)
     &') i,j,dj_dj(iparm,jparm),dj_de(iparm,jparm),dj_de(jparm,iparm),dj_dj_e(iparm,jparm),de_de(iparm,jparm)
   10 continue

c second derivatives involving scalek parameter:
      if(nparms.eq.1) then
        iparm=nparmcsf+nparmot+1
        do 12 j=1,nparmj+nparms
          jparm=nparmcsf+nparmot+j
          d2j(jparm,iparm)=d2j(jparm,iparm)+wt*(p*didk(j)+q*didk_old(j))
          d2j_e(jparm,iparm)=d2j_e(jparm,iparm)+wt*(p*didk(j)*enew+q*didk_old(j)*eoldr)
c         d2j(jparm,iparm)=0.d0
c         d2j_e(jparm,iparm)=0.d0
   12   continue
      endif

c second derivatives involving orbitals:
      if(nparmot.gt.0) then
        if(nparmcsf.gt.0) stop 'cant optimize both orbitals and csf yet'
c (easily doable but for testing dont optimize both for now)
        do i=1,nparmot
          iparm=i+nparmcsf
          do j=1,i
            jparm=j+nparmcsf
            d2j(iparm,jparm)=d2j(iparm,jparm)+wt*(p*detij_det(iparm,jparm)
     &                                       +q*detij_det_old(iparm,jparm))
            d2j_e(iparm,jparm)=d2j_e(iparm,jparm)+wt*(p*detij_det(iparm,jparm)*enew
     &                                       +q*detij_det_old(iparm,jparm)*eold)
            d2j(jparm,iparm)=d2j(iparm,jparm)
            d2j_e(jparm,iparm)=d2j_e(iparm,jparm)
          enddo
        enddo
      endif

c      do i=1,nparmot
c        write(*,*) 'i,j,d2j(i,j)=',(d2j(i,j),j=1,i)
c        write(*,*) 'i,j,detij_det(i,j)=',(detij_det(i,j),j=1,i)
c      enddo


      nparm0=nparmcsf+nparmot+nparms
      do 20 it=1,nctype
        nparm0=nparm0+npointa(it)
        do 20 i=1,nparma(it)
c swap next 2 lines
          iparm=nparm0+i
          if(iwjasa(i,it).eq.2) then
            d2j(iparm,iparm)=d2j(iparm,iparm)+wt*(p*d2d2a(it)+q*d2d2a_old(it))
            d2j_e(iparm,iparm)=d2j_e(iparm,iparm)+wt*(p*d2d2a(it)*enew+q*d2d2a_old(it)*eold)
      if(ipr.ge.1) write(6,'(''iparm,it,d2d2a(it),d2d2a_old(it),d2j(iparm,iparm)
     &,enew,eold,d2j_e(iparm,iparm)'',2i3,9d12.4)')
     &iparm,it,d2d2a(it),d2d2a_old(it),d2j(iparm,iparm),enew,eold,d2j_e(iparm,iparm)
            do 15 j=1,nparma(it)
              if(iwjasa(j,it).eq.1) then
                jparm=nparm0+j
                sav1=wt*(p*d1d2a(it)+q*d1d2a_old(it))
                sav2=wt*(p*d1d2a(it)*enew+q*d1d2a_old(it)*eold)
                if(jparm.gt.iparm) then
                  d2j(jparm,iparm)=d2j(jparm,iparm)+sav1
                  d2j_e(jparm,iparm)=d2j_e(jparm,iparm)+sav2
                 else
                  d2j(iparm,jparm)=d2j(iparm,jparm)+sav1
                  d2j_e(iparm,jparm)=d2j_e(iparm,jparm)+sav2
                endif
              endif
   15       continue
          endif
   20 continue

      nparm0=nparmcsf+nparmot+npointa(nctype)+nparma(nctype)+nparms
      do 30 isb=1,nspin2b
        if(isb.eq.2) nparm0=nparm0+nparmb(1)
        do 30 i=1,nparmb(isb)
          iparm=nparm0+i
          if(iwjasb(i,isb).eq.2) then
            d2j(iparm,iparm)=d2j(iparm,iparm)+wt*(p*d2d2b(isb)+q*d2d2b_old(isb))
            d2j_e(iparm,iparm)=d2j_e(iparm,iparm)+wt*(p*d2d2b(isb)*enew+q*d2d2b_old(isb)*eold)
            do 25 j=1,nparmb(isb)
              if(iwjasb(j,isb).eq.1) then
                jparm=nparm0+j
                sav1=wt*(p*d1d2b(isb)+q*d1d2b_old(isb))
                sav2=wt*(p*d1d2b(isb)*enew+q*d1d2b_old(isb)*eold)
                if(jparm.gt.iparm) then
                  d2j(jparm,iparm)=d2j(jparm,iparm)+sav1
                  d2j_e(jparm,iparm)=d2j_e(jparm,iparm)+sav2
                 else
                  d2j(iparm,jparm)=d2j(iparm,jparm)+sav1
                  d2j_e(iparm,jparm)=d2j_e(iparm,jparm)+sav2
                endif
              endif
   25       continue
          endif
   30 continue

c      do i=1,nparmot
c        write(*,*) 'i,d2j(i,j,)=',i,(d2j(i,j),j=1,nparmot)
c      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine grad_hess_jas_cum(wsum,enow)

      use optim_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include '../fit/fit.h'

      common /gradhessder/ dj(MPARM),dj_e(MPARM),dj_de(MPARM,MPARM),dj_dj(MPARM,MPARM),dj_dj_e(MPARM,MPARM)
     &,de(MPARM),d2j(MPARM,MPARM),d2j_e(MPARM,MPARM),de_e(MPARM),e2(MPARM),dj_e2(MPARM),de_de(MPARM,MPARM)
     &,w_i(MPARM),w_i_e(MPARM)

      common /gradjerr/ grad_jas_bcum(MPARMJ),grad_jas_bcm2(MPARMJ),
     &dj_e_bsum(MPARMJ),dj_bsum(MPARMJ),dj_e_save(MPARMJ),dj_save(MPARMJ),e_bsum
      common /gradjerrb/ ngrad_jas_blocks,ngrad_jas_bcum,nb_current

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

      common /contrl_opt2/ igradhess,iadd_diag_opt

      dimension dj_e_b(MPARMJ),dj_b(MPARMJ)

      if(igradhess.eq.0.or.ngrad_jas_blocks.eq.0) return

      do 10 i=1,nparmj+nparms
        dj_e_b(i)=dj_e(i)-dj_e_save(i)
   10   dj_b(i)=dj(i)-dj_save(i)

      e_bsum=e_bsum+enow
      do 20 i=1,nparmj+nparms
        dj_e_bsum(i)=dj_e_bsum(i)+dj_e_b(i)/wsum
   20   dj_bsum(i)=dj_bsum(i)+dj_b(i)/wsum

      nb_current=nb_current+1
      if(nb_current.eq.ngrad_jas_blocks)then
        nb_current=0
        ngrad_jas_bcum=ngrad_jas_bcum+1
        eb=e_bsum/dble(ngrad_jas_blocks)
        e_bsum=0
        do 30 i=1,nparmj+nparms
          gnow=2*(dj_e_bsum(i)-dj_bsum(i)*eb)/dble(ngrad_jas_blocks)
          grad_jas_bcum(i)=grad_jas_bcum(i)+gnow
          grad_jas_bcm2(i)=grad_jas_bcm2(i)+gnow**2
          dj_e_bsum(i)=0
   30     dj_bsum(i)=0
      endif

      do 40 i=1,nparmj+nparms
        dj_e_save(i)=dj_e(i)
   40   dj_save(i)=dj(i)

      return
      end
c-----------------------------------------------------------------------
      subroutine grad_hess_jas_save

      use atom_mod
      use slater_mod
      use optim_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include '../fit/fit.h'

      common /contrl_opt/ nparm,nsig,ncalls,iopt,ipr_opt
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!JT      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
!JT     &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
!JT     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
!JT     &,detu(MDETUD),detd(MDETUD)
!JT     &,ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
!JT     &,detij_det(MPARMD,MPARMD)
      common /bparm/ nspin2b,nocuspb

      common /gradhessdero/ deti_det_old(MPARMD),gvalue_old(MPARMJ),denergy_old(MPARM)
     &,d1d2a_old(MCTYPE),d2d2a_old(MCTYPE),d1d2b_old(2),d2d2b_old(2),didk_old(MPARMJ)
     &,detij_det_old(MPARMD,MPARMD)
      common /delocc/ denergy(MPARM)
      common /derivjas/ gvalue(MPARMJ),g(3,MELEC,MPARMJ),d2g(MPARMJ)
     &,go(MELEC,MELEC,MPARMJ)
      common /gradhessj_nonlin/ d1d2a(MCTYPE),d2d2a(MCTYPE),d1d2b(2),d2d2b(2),didk(MPARMJ)

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype
      common /contrl_opt2/ igradhess,iadd_diag_opt

      if(igradhess.eq.0) return

      do 5 i=1,nparmcsf+nparmot
        deti_det_old(i)=deti_det(i)
        denergy_old(i)=denergy(i)
        do 5 j=1,nparmcsf+nparmot
    5     detij_det_old(i,j)=detij_det(i,j)

      do 10 i=1,nparmj+nparms
        iparm=nparmcsf+nparmot+i
        gvalue_old(i)=gvalue(i)
        didk_old(i)=didk(i)
   10   denergy_old(iparm)=denergy(iparm)

      do 20 it=1,nctype
        if(ipr_opt.ge.3) write(6,'(''it,d1d2a(it),d2d2a_old(it)'',i3,2d12.4)') it,d1d2a(it),d2d2a(it)
        d1d2a_old(it)=d1d2a(it)
   20   d2d2a_old(it)=d2d2a(it)

      do 30 isb=1,nspin2b
        d1d2b_old(isb)=d1d2b(isb)
   30   d2d2b_old(isb)=d2d2b(isb)

      return
      end
c-----------------------------------------------------------------------
      subroutine grad_hess_jas_init

      use optim_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include '../fit/fit.h'

      common /contrl_opt/ nparm,nsig,ncalls,iopt,ipr_opt
      common /contrl_opt2/ igradhess,iadd_diag_opt

      common /gradhessder/ dj(MPARM),dj_e(MPARM),dj_de(MPARM,MPARM),dj_dj(MPARM,MPARM),dj_dj_e(MPARM,MPARM)
     &,de(MPARM),d2j(MPARM,MPARM),d2j_e(MPARM,MPARM),de_e(MPARM),e2(MPARM),dj_e2(MPARM),de_de(MPARM,MPARM)
     &,w_i(MPARM),w_i_e(MPARM)

      common /gradjerr/ grad_jas_bcum(MPARMJ),grad_jas_bcm2(MPARMJ),
     &dj_e_bsum(MPARMJ),dj_bsum(MPARMJ),dj_e_save(MPARMJ),dj_save(MPARMJ),e_bsum

      common /gradjerrb/ ngrad_jas_blocks,ngrad_jas_bcum,nb_current

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype

      if(igradhess.eq.0) return

      do 10 i=1,nparm
        dj(i)=0
        de(i)=0
        dj_e(i)=0
        de_e(i)=0
        dj_e2(i)=0
        e2(i)=0
        w_i(i)=0
        w_i_e(i)=0
        do 10 j=1,i
          dj_de(i,j)=0
          dj_de(j,i)=0
          dj_dj(i,j)=0
          dj_dj_e(i,j)=0
          d2j(i,j)=0
          d2j_e(i,j)=0
   10     de_de(i,j)=0

      e_bsum=0
      do 20 i=1,nparmj+nparms
        grad_jas_bcum(i)=0
        grad_jas_bcm2(i)=0
        dj_e_bsum(i)=0
        dj_bsum(i)=0
        dj_e_save(i)=0
   20   dj_save(i)=0

      nb_current=0
      ngrad_jas_bcum=0
      ngrad_jas_blocks=0

      return
      end
c-----------------------------------------------------------------------
      subroutine grad_hess_jas_dump(iu)

      use optim_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include '../fit/fit.h'

      common /gradhessder/ dj(MPARM),dj_e(MPARM),dj_de(MPARM,MPARM),dj_dj(MPARM,MPARM),dj_dj_e(MPARM,MPARM)
     &,de(MPARM),d2j(MPARM,MPARM),d2j_e(MPARM,MPARM),de_e(MPARM),e2(MPARM),dj_e2(MPARM),de_de(MPARM,MPARM)
     &,w_i(MPARM),w_i_e(MPARM)

      common /gradjerr/ grad_jas_bcum(MPARMJ),grad_jas_bcm2(MPARMJ),
     &dj_e_bsum(MPARMJ),dj_bsum(MPARMJ),dj_e_save(MPARMJ),dj_save(MPARMJ),e_bsum

      common /gradjerrb/ ngrad_jas_blocks,ngrad_jas_bcum,nb_current

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

      common /contrl_opt2/ igradhess,iadd_diag_opt

      if(igradhess.eq.0) return
c to do: write out which parameters are being varied -> check for restart
c Warning: Except for dj_de the rest are sym. so we do not really need to write entire matrix
      write(iu) nparmj
      write(iu) (dj(i),de(i),dj_e(i),de_e(i),dj_e2(i),e2(i),i=1,nparmj)
      write(iu) ((dj_de(i,j),j=1,nparmj),i=1,nparmj)
c     write(iu) ((dj_dj(i,j),dj_dj_e(i,j),j=1,i),i=1,nparmj)
      write(iu) ((dj_dj(i,j),dj_dj_e(i,j),j=1,nparmj),i=1,nparmj)
      write(iu) ((d2j(i,j),d2j_e(i,j),j=1,nparmj),i=1,nparmj)
      write(iu) ((de_de(i,j),j=1,nparmj),i=1,nparmj)
      if(ngrad_jas_blocks.gt.0)
     & write(iu) (grad_jas_bcum(i),grad_jas_bcm2(i),i=1,nparmj),ngrad_jas_bcum

      return
      end
c-----------------------------------------------------------------------
      subroutine grad_hess_jas_rstrt(iu)

      use optim_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include '../fit/fit.h'

      common /gradhessder/ dj(MPARM),dj_e(MPARM),dj_de(MPARM,MPARM),dj_dj(MPARM,MPARM),dj_dj_e(MPARM,MPARM)
     &,de(MPARM),d2j(MPARM,MPARM),d2j_e(MPARM,MPARM),de_e(MPARM),e2(MPARM),dj_e2(MPARM),de_de(MPARM,MPARM)

      common /gradjerr/ grad_jas_bcum(MPARMJ),grad_jas_bcm2(MPARMJ),
     &dj_e_bsum(MPARMJ),dj_bsum(MPARMJ),dj_e_save(MPARMJ),dj_save(MPARMJ),e_bsum

      common /gradjerrb/ ngrad_jas_blocks,ngrad_jas_bcum,nb_current

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

      common /contrl_opt2/ igradhess,iadd_diag_opt

      if(igradhess.eq.0) return

      read(iu) nparmj
      read(iu) (dj(i),de(i),dj_e(i),de_e(i),dj_e2(i),e2(i),i=1,nparmj)
      read(iu) ((dj_de(i,j),j=1,nparmj),i=1,nparmj)
c     read(iu) ((dj_dj(i,j),dj_dj_e(i,j),j=1,i),i=1,nparmj)
      read(iu) ((dj_dj(i,j),dj_dj_e(i,j),j=1,nparmj),i=1,nparmj)
      read(iu) ((d2j(i,j),d2j_e(i,j),j=1,nparmj),i=1,nparmj)
      read(iu) ((de_de(i,j),j=1,nparmj),i=1,nparmj)
      if(ngrad_jas_blocks.gt.0)
     & read(iu) (grad_jas_bcum(i),grad_jas_bcm2(i),i=1,nparmj),ngrad_jas_bcum

      do 10 i=1,nparmj
        dj_e_save(i)=dj_e(i)
   10   dj_save(i)=dj(i)

      return
      end
c-----------------------------------------------------------------------
      subroutine grad_hess_jas_fin(passes,eave)

      use deriv_mod
      use opt_lin_mod
      use opt_nwt_mod
      use atom_mod
      use optim_mod
      implicit real*8(a-h,o-z)

      character*20 fmt

      common /contrl_opt/ nparm,nsig,ncalls,iopt,ipr_opt
      common /gradhessder/ dj(MPARM),dj_e(MPARM),dj_de(MPARM,MPARM),dj_dj(MPARM,MPARM),dj_dj_e(MPARM,MPARM)
     &,de(MPARM),d2j(MPARM,MPARM),d2j_e(MPARM,MPARM),de_e(MPARM),e2(MPARM),dj_e2(MPARM),de_de(MPARM,MPARM)
     &,w_i(MPARM),w_i_e(MPARM)

      common /gradjerr/ grad_jas_bcum(MPARMJ),grad_jas_bcm2(MPARMJ),
     &dj_e_bsum(MPARMJ),dj_bsum(MPARMJ),dj_e_save(MPARMJ),dj_save(MPARMJ),e_bsum

      common /gradjerrb/ ngrad_jas_blocks,ngrad_jas_bcum,nb_current

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent

      common /contrl_opt2/ igradhess,iadd_diag_opt

      common /gradhess/ grad(MPARM),grad_var(MPARM),hess(MPARM,MPARM),hess_var(MPARM,MPARM),gerr(MPARM),
     &add_diag(3),energy(3),energy_sigma(3),energy_err(3),force(3),force_err(3),
     &eig_min,eig_max,p_var,tol_energy,nopt_iter,nblk_max
      common /linear/ ham(MPARM,MPARM),ovlp(MPARM,MPARM),coef(MPARM,MPARM)


      errn(x,x2,n)=dsqrt(dabs(x2/dble(n)-(x/dble(n))**2)/dble(n))

      if(igradhess.eq.0) return

c Compute gradient.  For the gradient of the variance, grad_var, try both unweighted and weighted.
      grad_norm=0
      grad_var_norm=0
      do 20 i=1,nparm
c       write(6,'(''i,dj_e(i),eave,dj(i),passes'',i5,9d12.4)') i,dj_e(i),eave,dj(i),passes
        grad(i)=2*(dj_e(i)-eave*dj(i))/passes
c       grad_var(i)=2*(de_e(i)-eave*de(i))/passes
        grad_var(i)=2*(de_e(i)-eave*de(i) + dj_e2(i)-dj(i)*e2(i)/passes
     &  -2*eave*(dj_e(i)-eave*dj(i)))/passes  ! equivalent to next line
c    &  -eave*grad(i)*passes)/passes          ! equivalent to previous line
        grad_norm=grad_norm+grad(i)**2
   20   grad_var_norm=grad_var_norm+grad_var(i)**2
      grad_norm=sqrt(grad_norm)
      grad_var_norm=sqrt(grad_var_norm)

c Before computing Hessian, compute the ratio to be used in rescaling the low fluctuation terms:
c Use symmetrized dj_de
      topsum_j=0
      botsum_j=0
      topsum_csf_j=0
      botsum_csf_j=0
      do 25 i=1,nparm
        do 25 j=1,i
           if(j.gt.nparmcsf+nparmot) then
             bot_j=dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes
             top_j=bot_j+2*(2*(dj_dj_e(i,j)-eave*dj_dj(i,j))-grad(i)*dj(j)-grad(j)*dj(i))
             if(ipr_opt.ge.-2) write(6,'(''i,j,top_j,bot_j,ratio='',i2,i3,1p2d11.3,0pf6.2)') i,j,top_j,bot_j,top_j/bot_j
             botsum_j=botsum_j+bot_j
             topsum_j=topsum_j+top_j
             check=(2*(2*(dj_dj_e(i,j)-eave*dj_dj(i,j))
     &             -grad(i)*dj(j)-grad(j)*dj(i))
     &             +dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes)
             if(abs(top_j-check).gt.1.d-5*abs(top_j)) then
c            if(abs(top_j-check).gt.1.d-2*abs(top_j)) then
               write(6,'(''i,j,dj_de(i,j),dj_de(j,i),dj(i),de(j),dj_dj_e(i,j),eave*dj_dj(i,j),eave,grad(j)='',2i4,9d12.4)')
     &         i,j,dj_de(i,j),dj_de(j,i),dj(i),de(j),dj_dj_e(i,j),eave*dj_dj(i,j),eave,grad(j)
               write(6,'(''bot_j,top_j,check'',3d21.14,d12.4)') bot_j,top_j,check,top_j-check
               stop 'top_j check'
             endif
           endif
   25 continue
      ratio_j=topsum_j/botsum_j
c     ratio_csf_j=topsum_csf_j/botsum_csf_j
      if(ipr_opt.ge.-4) write(6,'(''ratio_j,grad_norm,grad_var_norm='',f9.5,9f11.5)') ratio_j,grad_norm,grad_var_norm

c Note the "2" in 2*dj_dj_e(i,j) and 2*eave*dj_dj(i,j) is because d2j is not the
c second derivative of psi divided by psi, but rather the second derivative of the log of psi.
c Note by A.D.Guclu: above statement is true for jastrow parameters. For determinantal
c parameters d2j IS taken to be second derivative of psi divided by psi.

      if(idtask.eq.0) then
        open(20,file='hess_pieces',status='unknown')
        open(21,file='grad_hess',status='unknown')
c       open(22,file='lin_ham_overlap',status='unknown')
       else
        open(20,file='/dev/null')
        open(21,file='/dev/null')
c       open(22,file='/dev/null')
      endif

c In hess_pieces write the expression for hess in our PRL but then evaluate hess for writing in
c grad_hess
      write(20,'(5i4,'' nparm,nparmcsf,nparmot,nparmj,nparms='')') nparm,nparmcsf,nparmot,nparmj,nparms
      write(fmt,'(''('',i3,''g14.6,a)'')') nparm
      write(20,fmt) (grad(i),i=1,nparm), ' (grad(i),i=1,nparm)'
      write(20,fmt) (grad_var(i),i=1,nparm), ' (grad_var(i),i=1,nparm)'
c     write(20,'(''grad    ='',50f8.4)') (grad(i),i=1,nparm)
c     write(20,'(''grad_var='',50f8.4)') (grad_var(i),i=1,nparm)
      write(20,'('' j  k     A        B1/2        B2        C_ij        C_ji       D-sym(C)     h(UF)'')')
      do 30 i=1,nparm
        do 30 j=1,i
casym   do 30 j=1,nparm
          if(i.le.nparmcsf+nparmot) then
            second_der=0
           elseif(j.gt.nparmcsf+nparmot) then
            second_der=d2j_e(i,j)-eave*d2j(i,j) + dj_dj_e(i,j)-eave*dj_dj(i,j)
           else
            second_der=dj_dj_e(i,j)-eave*dj_dj(i,j)
          endif
c         hess(i,j)=(2*( d2j_e(i,j)-eave*d2j(i,j) +2*(dj_dj_e(i,j)-eave*dj_dj(i,j))
          hess(i,j)=(2*( second_der + (dj_dj_e(i,j)-eave*dj_dj(i,j))
     &             -grad(i)*dj(j)-grad(j)*dj(i) )
     &             +dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes)/passes
c         write(20,'(2i3,2f13.4,9f15.4)') i,j, 2*(d2j_e(i,j)-eave*d2j(i,j))/passes,2*(dj_dj_e(i,j)-eave*dj_dj(i,j))/passes,
          write(20,'(i2,i3,f9.4,1p,20g12.5)') i,j, 2*(d2j_e(i,j)-eave*d2j(i,j))/passes,2*(dj_dj_e(i,j)-eave*dj_dj(i,j))/passes,
     &    2*(-grad(i)*dj(j)-grad(j)*dj(i))/passes, (dj_de(i,j)+dj_de(j,i))/passes,-(dj(i)*de(j)+dj(j)*de(i))/passes**2,
c    &    2*(-grad(i)*dj(j)-grad(j)*dj(i))/passes, dj_de(i,j)/passes,dj_de(j,i)/passes,-dj(i)*de(j)/passes**2,-dj(j)*de(i)/passes**2,
     &    hess(i,j),
     &    -dj(i)/passes,-dj(j)/passes, de(i)/passes,de(j)/passes
          if(i.le.nparmcsf+nparmot) then                                       ! i,j < nparmcsf
c Note: unlike for jastrow, for determinant params d2j is second deriv of psi divided by psi.
c The d2j are not needed for CSF parameters but are needed for the orbital parameters.
            hess(i,j)=(2*( d2j_e(i,j)-eave*d2j(i,j)+ dj_dj_e(i,j)-eave*dj_dj(i,j)
     &      -grad(i)*dj(j)-grad(j)*dj(i) )
     &      +dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes)/passes
           elseif(j.gt.nparmcsf+nparmot) then                                  ! i,j > nparmcsf
            if(mod(iopt/100,10).eq.1) then ! do not rescale Hessian
              hess(i,j)=(2*( second_der + dj_dj_e(i,j)-eave*dj_dj(i,j)
     &        -grad(i)*dj(j)-grad(j)*dj(i) )
     &        +dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes)/passes
             else                          ! use rescaled Hessian for Jastrow parameters only
              hess(i,j)=(2*( d2j_e(i,j)-eave*d2j(i,j))
     &        +ratio_j*(dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes))/passes
            endif
           else   ! i>nparmcsf, j<nparmcsf
c           hess(i,j)=ratio_csf_j*(2*( 2*(dj_dj_e(i,j)-eave*dj_dj(i,j))
            hess(i,j)=(2*( 2*(dj_dj_e(i,j)-eave*dj_dj(i,j))
     &      -grad(i)*dj(j)-grad(j)*dj(i) )
     &      +dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes)/passes
          endif

c         hess(i,j)=(2*( d2j_e(i,j)-eave*d2j(i,j))
c    &             +ratio_j*(dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes))/passes
   30     hess_var(i,j)=2*((de_de(i,j)-grad(i)*de(j)-grad(j)*de(i))/passes+grad(i)*grad(j))
casym     if(j.le.i) then
casym       hess(i,j)=(2*(2*dj_dj_e(i,j)-grad(i)*dj(j)-grad(j)*dj(i)-2*eave*dj_dj(i,j)
casym&               +(d2j_e(i,j)-eave*d2j(i,j))
casym&               +dj_de(i,j)-dj(i)*de(j)/passes))/passes
casym      else
casym       hess(i,j)=(2*(2*dj_dj_e(j,i)-grad(i)*dj(j)-grad(j)*dj(i)-2*eave*dj_dj(j,i)
casym&               +(d2j_e(j,i)-eave*d2j(j,i))
casym&               +dj_de(i,j)-dj(i)*de(j)/passes))/passes
casym     endif
casym   30 continue

      close(20)

      write(21,'(''1 ntimes'')')
      write(21,'(i3,'' nctype'')') nctype
      write(21,'(i5,'' .1 0. nparm,add_diag,p_var'')') nparm
      write(21,'(''grad    '',(1p5g22.14))') (grad(i),i=1,nparm)
      write(21,'(''grad_var'',(1p5g22.14))') (grad_var(i),i=1,nparm)
      write(21,'(''hess    '',(1p5g22.14))') ((hess(i,j),j=1,i),i=1,nparm)
      write(21,'(''hess_var'',(1p5g22.14))') ((hess_var(i,j),j=1,i),i=1,nparm)
casym write(21,'(''hess'',(1p5g22.14))') ((hess(i,j),j=1,nparm),i=1,nparm)
      if(ngrad_jas_blocks.gt.0) then
        do 40 i=1,nparm
   40     gerr(i)=errn(grad_jas_bcum(i),grad_jas_bcm2(i),ngrad_jas_bcum)
         write(21,'(''dgrad'',(1p5g22.14))') (gerr(i),i=1,nparm)
      endif

      close(21)

c compute <dj H dj> and <dj dj> for "linear method"
      open(22,file='lin_jas.dat',status='unknown')
      write(6,'(''opening lin_jas.dat'')')
c here, ham = <dj H dj>
      ham(1,1)=eave
      do 45 i=1,nparm
        if(mod(iopt/10,10).le.1) then
          ham(1,i+1)=(de(i)+dj_e(i))/passes
          ham(i+1,1)=dj_e(i)/passes
         elseif(mod(iopt/10,10).eq.2) then
          ham(1,i+1)=dj_e(i)/passes
          ham(i+1,1)=ham(1,i+1)
        endif
   45 continue

      do 47 i=1,nparm
        do 47 j=1,i
          if(mod(iopt/10,10).le.1) then
            ham(i+1,j+1)=(dj_de(i,j)+dj_dj_e(i,j))/passes
            ham(j+1,i+1)=(dj_de(j,i)+dj_dj_e(i,j))/passes
           elseif(mod(iopt/10,10).eq.2) then
            ham(i+1,j+1)=(0.5d0*(dj_de(i,j)+dj_de(j,i)-(dj(i)*de(j)+dj(j)*de(i))/passes)+dj_dj_e(i,j))/passes
c           ham(i+1,j+1)=(0.5d0*(dj_de(i,j)+dj_de(j,i))+dj_dj_e(i,j))/passes
            ham(j+1,i+1)=ham(i+1,j+1)
          endif
   47 continue

c Symmetrize H
c     if(mod(iopt/10,10).eq.2) then
c       do 48 i=1,nparm
c         do 48 j=1,i-1
c           ham(i,j)=0.5d0*(ham(i,j)+ham(j,i))
c 48        ham(j,i)=ham(i,j)
c     endif

c     write(22,'(''dj_dj_e(i,i)'',20f9.4)') (dj_dj_e(i,i)/passes,i=1,nparm)
c     write(22,'(''dj_de(i,i)'',20f9.4)') (dj_de(i,i)/passes,i=1,nparm)
c     write(22,'(''dj(i)*de(i)'',20f9.4)') (dj(i)*de(i)/passes**2,i=1,nparm)
c     write(22,'(''1 ntimes'')')
      write(22,'(i3,'' nctype'')') nctype
      write(22,'(i5,'' 0 nparm,isym'')') nparm
c     write(22,'(''dj H dj'',(1p5g22.14))') ((ham(i,j),j=1,nparm+1),i=1,nparm+1)
      write(22,'((1p5g22.14))') ((ham(i,j),j=1,nparm+1),i=1,nparm+1)

c Temp: write hamiltonian with hopefully lower variance.  It is not clear if
c this should hurt or help.  It could hurt because we are violating Peter's zero-variance condition.
c In practice it seems to neither hurt nor help.
c here, ham = <dj H dj>
c     ham(1,1)=eave
c     do 51 i=1,nparm
c       ham(1,i+1)=dj_e(i)/passes
c 51    ham(i+1,1)=dj_e(i)/passes
c     do 52 i=1,nparm
c       do 52 j=1,i
c         ham(i+1,j+1)=(dj_de(i,j)-dj(i)*de(j)/passes+dj_dj_e(i,j))/passes
c 52      ham(j+1,i+1)=(dj_de(j,i)-dj(j)*de(i)/passes+dj_dj_e(i,j))/passes

c here, ovlp = <dj dj>
      ovlp(1,1)=1
      do 55 i=1,nparm
  55    ovlp(i+1,1)=dj(i)/passes

      do 60 i=1,nparm
        do 60 j=1,i
  60      ovlp(i+1,j+1)=dj_dj(i,j)/passes
c     write(22,'(''dj dj'',(1p5g22.14))') ((ovlp(i,j),j=1,i),i=1,nparm+1)
      write(22,'((1p5g22.14))') ((ovlp(i,j),j=1,i),i=1,nparm+1)
      close(22)

      return
      end

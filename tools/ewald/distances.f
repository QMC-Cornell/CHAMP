      subroutine distances(x,pe)
c Written by Cyrus Umrigar
c calculate interparticle distances

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'pseudo.h'

      common /dim/ ndim
      common /contrl_per/ iperiodic,ibasis
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      common /dot/ w0,p1,p2,p3,p4

c Warning: temporary
      common /tempor2/ pe_en,pe_ee

c     dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),
c    &rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      dimension x(3,*)

c  pe from nucleus-nucleus repulsion
      pe=pecent

      if(iperiodic.eq.0) then

c Calculate e-N inter-particle distances
        do 26 ic=1,ncent
          do 26 i=1,nelec
            r_en(i,ic)=0
            do 25 k=1,ndim
              rvec_en(k,i,ic)=x(k,i)-cent(k,ic)
   25         r_en(i,ic)=r_en(i,ic)+rvec_en(k,i,ic)**2
            r_en(i,ic)=dsqrt(r_en(i,ic))
            if(nloc.eq.0) pe=pe-znuc(iwctype(ic))/r_en(i,ic)
            if(nloc.eq.-1) pe=pe+0.5d0*(w0*r_en(i,ic))**2
            if(nloc.eq.-2) pe=pe+p1*rvec_en(1,i,ic)**4+p2*rvec_en(2,i,ic)**4
     &      -2*p3*(rvec_en(1,i,ic)*rvec_en(2,i,ic))**2
     &      +p4*(rvec_en(1,i,ic)-rvec_en(2,i,ic))*rvec_en(1,i,ic)*rvec_en(2,i,ic)*r_en(i,ic)
   26   continue

c Calculate e-e inter-particle distances
        ij=0
        do 29 i=2,nelec
          do 29 j=1,i-1
            ij=ij+1
            r_ee(ij)=0
            do 28 k=1,ndim
              rvec_ee(k,ij)=x(k,i)-x(k,j)
   28         r_ee(ij)=r_ee(ij)+rvec_ee(k,ij)**2
            r_ee(ij)=dsqrt(r_ee(ij))
            pe=pe+1/r_ee(ij)
   29   continue

       else

        call pot_ee_ewald(x,pe_ee)
        call pot_en_ewald(x,pe_en)
        pe=pe+pe_en+pe_ee

c     write(6,*) 'in distances'
      write(6,'(''r_en(i,j)'',9f9.5)') ((r_en(i,j),i=1,nelec),j=1,ncent)
      write(6,'(''r_ee(ij)'',9f9.5)') (r_ee(ij),ij=1,nelec*(nelec-1)/2)
c     write(6,'(''rvec_ee(k,ij)'',9f12.4)') ((rvec_ee(k,ij),k=1,ndim),ij=1,nelec*(nelec-1)/2)
      if(ipr.ge.3) write(6,'(''pe,pe_en(loc),pe_ee'',11f9.5)') pe,pe_en,pe_ee

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine distancese(iel,x)
c Written by Cyrus Umrigar
c calculate distances of electron iel to all other particles

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'ewald.h'

      common /dim/ ndim
      common /contrl_per/ iperiodic,ibasis
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      common /distances_sav/ rshift_sav(3,MCENT),rvec_en_sav(3,MCENT),r_en_sav(MCENT),rvec_ee_sav(3,MELEC),r_ee_sav(MELEC)
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),rlatt_sim_inv(3,3),glatt_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

      dimension x(3,*)

c Calculate e-N inter-particle distances
      do 30 ic=1,ncent
        r_en_sav(ic)=r_en(iel,ic)
        do 10 k=1,ndim
          rshift_sav(k,ic)=rshift(k,iel,ic)
          rvec_en_sav(k,ic)=rvec_en(k,iel,ic)
   10     rvec_en(k,iel,ic)=x(k,iel)-cent(k,ic)
        if(iperiodic.eq.0) then
          r_en(iel,ic)=0
          do 20 k=1,ndim
   20       r_en(iel,ic)=r_en(iel,ic)+rvec_en(k,iel,ic)**2
          r_en(iel,ic)=dsqrt(r_en(iel,ic))
         else
          call find_image4(rshift(1,iel,ic),rvec_en(1,iel,ic),r_en(iel,ic),rlatt,rlatt_inv)
        endif
   30 continue

c Calculate e-e inter-particle distances
      do 60 jj=1,nelec

        if(jj.eq.iel) goto 60
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif
        ij=((i-1)*(i-2))/2+j

        r_ee_sav(jj)=r_ee(ij)
        do 40 k=1,ndim
          rvec_ee_sav(k,jj)=rvec_ee(k,ij)
   40     rvec_ee(k,ij)=x(k,i)-x(k,j)
        if(iperiodic.eq.0) then
          r_ee(ij)=0
          do 50 k=1,ndim
   50       r_ee(ij)=r_ee(ij)+rvec_ee(k,ij)**2
          r_ee(ij)=dsqrt(r_ee(ij))
         else
          call find_image3(rvec_ee(1,ij),r_ee(ij),rlatt_sim,rlatt_sim_inv)
        endif
   60 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine distancese_restore(iel)
c Written by Cyrus Umrigar
c restore interparticle distances (called if move rejected)

      implicit real*8(a-h,o-z)

      include 'vmc.h'

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      common /distances_sav/ rshift_sav(3,MCENT),rvec_en_sav(3,MCENT),r_en_sav(MCENT),rvec_ee_sav(3,MELEC),r_ee_sav(MELEC)

c Restore e-N inter-particle distances
      do 25 ic=1,ncent
        r_en(iel,ic)=r_en_sav(ic)
        do 25 k=1,ndim
          rshift(k,iel,ic)=rshift_sav(k,ic)
   25     rvec_en(k,iel,ic)=rvec_en_sav(k,ic)

c Restore e-e inter-particle distances
      do 29 jj=1,nelec

        if(jj.eq.iel) goto 29
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif
        ij=((i-1)*(i-2))/2+j

        r_ee(ij)=r_ee_sav(jj)
        do 28 k=1,ndim
   28     rvec_ee(k,ij)=rvec_ee_sav(k,jj)
   29 continue

      return
      end
c-----------------------------------------------------------------------

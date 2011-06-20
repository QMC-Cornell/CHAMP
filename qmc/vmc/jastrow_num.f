      subroutine jastrow_num(x,v,d2,div_vj,value)
c Written by Cyrus Umrigar
c **Warning** This routine needs to be upgraded to calculate distances
c correctly for periodic systems if we add in capability to use
c numerical Laplacian for periodic systems.
      use constants_mod
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use wfsec_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar6_mod
      use bparm_mod
      use distance_mod
      implicit real*8(a-h,o-z)

      parameter (eps=.5d-4,eps2=2.d0*eps,eps4=4.d0*eps,epssq=eps**2
     &,eps2sq=eps2**2)
      parameter (d1b12=8.333333333333333d-2,d2b3=0.666666666666667d0,
     &d4b3=1.333333333333333d0)

c subroutine to calculate jastrow factor,its derivatives
c and the potential
c Warning: div_vj not yet implememnted

      dimension x(3,*),v(3,*),div_vj(nelec), rtemp(3)
      dimension rp(3,nelec,ncent),rm(3,nelec,ncent)
     &,rp2(3,nelec,ncent),rm2(3,nelec,ncent)
     &,rrp(3,nelec_pair),rrm(3,nelec_pair),rrp2(3,nelec_pair),rrm2(3,nelec_pair)

      do 10 i=1,nelec
      v(1,i)=zero
      v(2,i)=zero
   10 v(3,i)=zero
      d2=zero
      usum=zero

      if(ijas.eq.1) then
        ncentt=1
       else
        ncentt=ncent
      endif

c Calculate e-N and e-e inter-particle distances
        do 29 ic=1,ncent
        ij=0
        do 29 i=1,nelec
          r_en(i,ic)=zero
          if(iperiodic.eq.1 .and. nloc.eq.-4 .and. ic.eq.1) then
            r_en(i,ic)=x(2,i)**2  ! distance from electron to axis
            rp(2,i,ic)=r_en(i,ic)+(x(2,i))*eps2+epssq
            rm(2,i,ic)=r_en(i,ic)-(x(2,i))*eps2+epssq
            rp2(2,i,ic)=r_en(i,ic)+(x(2,i))*eps4+eps2sq
            rm2(2,i,ic)=r_en(i,ic)-(x(2,i))*eps4+eps2sq
            rp(2,i,ic)=dsqrt(rp(2,i,ic))
            rm(2,i,ic)=dsqrt(rm(2,i,ic))
            rp2(2,i,ic)=dsqrt(rp2(2,i,ic))
            rm2(2,i,ic)=dsqrt(rm2(2,i,ic))
            r_en(i,ic)=dsqrt(r_en(i,ic))
            rp(1,i,ic) = r_en(i,ic)
            rm(1,i,ic) = r_en(i,ic)
            rp2(1,i,ic) = r_en(i,ic)
            rm2(1,i,ic) = r_en(i,ic)
          else
            do 25 m=1,ndim
   25         r_en(i,ic)=r_en(i,ic)+(x(m,i)-cent(m,ic))**2
            do 26 m=1,ndim
              rp(m,i,ic)=r_en(i,ic)+(x(m,i)-cent(m,ic))*eps2+epssq
              rm(m,i,ic)=r_en(i,ic)-(x(m,i)-cent(m,ic))*eps2+epssq
              rp2(m,i,ic)=r_en(i,ic)+(x(m,i)-cent(m,ic))*eps4+eps2sq
              rm2(m,i,ic)=r_en(i,ic)-(x(m,i)-cent(m,ic))*eps4+eps2sq
              rp(m,i,ic)=dsqrt(rp(m,i,ic))
              rm(m,i,ic)=dsqrt(rm(m,i,ic))
              rp2(m,i,ic)=dsqrt(rp2(m,i,ic))
   26         rm2(m,i,ic)=dsqrt(rm2(m,i,ic))
            r_en(i,ic)=dsqrt(r_en(i,ic))
          endif

          do 28 j=1,i-1
            ij=ij+1
            if(iperiodic.eq.0) then
              r_ee(ij)=(x(1,i)-x(1,j))**2 + (x(2,i)-x(2,j))**2
     &        + (x(3,i)-x(3,j))**2
              do 27 m=1,ndim
                rrp(m,ij)=r_ee(ij)+(x(m,i)-x(m,j))*eps2+epssq
                rrm(m,ij)=r_ee(ij)-(x(m,i)-x(m,j))*eps2+epssq
                rrp2(m,ij)=r_ee(ij)+(x(m,i)-x(m,j))*eps4+eps2sq
                rrm2(m,ij)=r_ee(ij)-(x(m,i)-x(m,j))*eps4+eps2sq
                rrp(m,ij)=dsqrt(rrp(m,ij))
                rrm(m,ij)=dsqrt(rrm(m,ij))
                rrp2(m,ij)=dsqrt(rrp2(m,ij))
   27           rrm2(m,ij)=dsqrt(rrm2(m,ij))
            else
              rtemp(:) = x(:,i)-x(:,j)
              call find_image_1d(rtemp, r_ee(ij))
              r_ee(ij) = r_ee(ij)*r_ee(ij)  ! since the below lines need
                                            !     r_ee**2
              rrp(:,ij)=dsqrt(r_ee(ij)+rtemp(:)*eps2+epssq)
              rrm(:,ij)=dsqrt(r_ee(ij)-rtemp(:)*eps2+epssq)
              rrp2(:,ij)=dsqrt(r_ee(ij)+rtemp(:)*eps4+eps2sq)
              rrm2(:,ij)=dsqrt(r_ee(ij)-rtemp(:)*eps4+eps2sq)
            endif
   28       r_ee(ij)=dsqrt(r_ee(ij))
   29   continue

c if nelec is < 2 then no pairs so exit
      if(nelec.lt.2) goto 60

      ij=0
      do 40 i=2,nelec
c3      jk=0
        im1=i-1
        do 40 j=1,im1
          ij=ij+1
          if(i.le.nup .or. j.gt.nup) then
c           parallel spins
            if(nspin2.ge.2) then
              is=2
              sspinn=one
              if(nspin2.eq.3 .and. j.gt.nup) is=3
             else
              is=1
              if(ndim.eq.3) then
                sspinn=half
               elseif(ndim.eq.2) then
                sspinn=third
              endif
            endif
           else
c           anti-parallel spins
            sspinn=one
            is=1
          endif

          if(ijas.ge.3.and.ijas.le.6) then
            sspinn=one
            ipar=0
            if(i.le.nup .or. j.gt.nup) then
              isb=is
              if(nspin2b.eq.2) then
                isb=2
               elseif(nocuspb.eq.0) then
                if(ndim.eq.3) then
                  sspinn=half
                 elseif(ndim.eq.2) then
                  sspinn=third
                endif
              endif
              ipar=1
             else
              isb=1
            endif

            psibc=psib(r_ee(ij),isb,ipar)
            usum=usum+psibc

            do 30 m=1,ndim
              psibpi=psib(rrp(m,ij),isb,ipar)
              psibmi=psib(rrm(m,ij),isb,ipar)
              psibpi2=psib(rrp2(m,ij),isb,ipar)
              psibmi2=psib(rrm2(m,ij),isb,ipar)

              ftmp=(-d1b12*(psibpi2-psibmi2)+d2b3*(psibpi-psibmi))/eps
              v(m,i)=v(m,i)+ftmp
              v(m,j)=v(m,j)-ftmp

   30         d2=d2+two*(-d1b12*((psibmi2-psibc)+(psibpi2-psibc))
     &                   +d4b3 *((psibmi -psibc)+(psibpi -psibc)))
          endif

          do 40 ic=1,ncentt
          it=iwctype(ic)

          wtj=one
          if(ijas.eq.2) wtj=one/ncentt

          sspin=sspinn*wtj
          psic=psi(r_ee(ij),r_en(i,ic),r_en(j,ic),it)
          usum=usum+psic

          do 35 m=1,ndim

            psipi=psi(rrp(m,ij),rp(m,i,ic),r_en(j,ic),it)
            psimi=psi(rrm(m,ij),rm(m,i,ic),r_en(j,ic),it)
            psipj=psi(rrm(m,ij),r_en(i,ic),rp(m,j,ic),it)
            psimj=psi(rrp(m,ij),r_en(i,ic),rm(m,j,ic),it)
            psipi2=psi(rrp2(m,ij),rp2(m,i,ic),r_en(j,ic),it)
            psimi2=psi(rrm2(m,ij),rm2(m,i,ic),r_en(j,ic),it)
            psipj2=psi(rrm2(m,ij),r_en(i,ic),rp2(m,j,ic),it)
            psimj2=psi(rrp2(m,ij),r_en(i,ic),rm2(m,j,ic),it)

            v(m,i)=v(m,i)+(-d1b12*(psipi2-psimi2)
     &      + d2b3*(psipi-psimi))/eps
            v(m,j)=v(m,j)+(-d1b12*(psipj2-psimj2)
     &      + d2b3*(psipj-psimj))/eps

   35       d2=d2-d1b12*((psimi2-psic)+(psipi2-psic))
     &           +d4b3 *((psimi -psic)+(psipi -psic))
     &           -d1b12*((psimj2-psic)+(psipj2-psic))
     &           +d4b3 *((psimj -psic)+(psipj -psic))

c Calculate e-e-e-n terms.  Presently only for parallel-spin electrons.
c Needs to be checked and more important 2 par and 1 anti par terms need
c to be added.
c3        if(i3body.ge.1) then
c3        do 37 k=1,j-1
c3          ik=((i-1)*(i-2))/2+k
c3          jk=jk+1
c3          f31c=f31(r_ee(ij),r_ee(ik),r_ee(jk),r_en(i,ic),
c3   &      r_en(j,ic),r_en(k,ic))
c3          usum=usum+f31c
c3          do 37 m=1,ndim
c3          f31pi=f31(rrp(m,ij),rrp(m,ik),r_ee(jk),rp(m,i,ic),
c3   &      r_en( j,ic),r_en( k,ic))
c3          f31mi=f31(rrm(m,ij),rrm(m,ik),r_ee(jk),rm(m,i,ic),
c3   &      r_en( j,ic),r_en( k,ic))
c3          f31pj=f31(rrm(m,ij),r_ee(ik),rrp(m,jk),r_en( i,ic),
c3   &      rp(m,j,ic),r_en( k,ic))
c3          f31mj=f31(rrp(m,ij),r_ee(ik),rrm(m,jk),r_en( i,ic),
c3   &      rm(m,j,ic),r_en( k,ic))
c3          f31pk=f31(r_ee(ij),rrm(m,ik),rrm(m,jk),r_en( i,ic),
c3   &      r_en( j,ic),rp(m,k,ic))
c3          f31mk=f31(r_ee(ij),rrp(m,ik),rrp(m,jk),r_en( i,ic),
c3   &      r_en( j,ic),rm(m,k,ic))
c3          f31pi2=f31(rrp2(m,ij),rrp2(m,ik),r_ee(jk),rp2(m,i,ic),
c3   &      r_en( j,ic),r_en( k,ic))
c3          f31mi2=f31(rrm2(m,ij),rrm2(m,ik),r_ee(jk),rm2(m,i,ic),
c3   &      r_en( j,ic),r_en( k,ic))
c3          f31pj2=f31(rrm2(m,ij),r_ee(ik),rrp2(m,jk),r_en( i,ic),
c3   &      rp2(m,j,ic),r_en( k,ic))
c3          f31mj2=f31(rrp2(m,ij),r_ee(ik),rrm2(m,jk),r_en( i,ic),
c3   &      rm2(m,j,ic),r_en( k,ic))
c3          f31pk2=f31(r_ee(ij),rrm2(m,ik),rrm2(m,jk),r_en( i,ic),
c3   &      r_en( j,ic),rp2(m,k,ic))
c3          f31mk2=f31(r_ee(ij),rrp2(m,ik),rrp2(m,jk),r_en( i,ic),
c3   &      r_en( j,ic),rm2(m,k,ic))
c3          v(m,i)=v(m,i)+(-d1b12*(f31pi2-f31mi2) + d2b3*(f31pi-f31mi))
c3   &      /eps
c3          v(m,j)=v(m,j)+(-d1b12*(f31pj2-f31mj2) + d2b3*(f31pj-f31mj))
c3   &      /eps
c3          v(m,k)=v(m,k)+(-d1b12*(f31pk2-f31mk2) + d2b3*(f31pk-f31mk))
c3   &      /eps
c3 37       d2=d2-d1b12*((f31mi2-f31c)+(f31pi2-f31c))
c3   &           +d4b3 *((f31mi -f31c)+(f31pi -f31c))
c3   &           -d1b12*((f31mj2-f31c)+(f31pj2-f31c))
c3   &           +d4b3 *((f31mj -f31c)+(f31pj -f31c))
c3   &           -d1b12*((f31mk2-f31c)+(f31pk2-f31c))
c3   &           +d4b3 *((f31mk -f31c)+(f31pk -f31c))
c3          endif

   40   continue

      if(ijas.ge.3.and.ijas.le.6) then
        do 50 i=1,nelec
          is=1

          do 50 ic=1,ncentt
            it=iwctype(ic)

            psiac=psia(r_en(i,ic),it)
            usum=usum+psiac

            do 45 m=1,ndim
              psiapi=psia(rp(m,i,ic),it)
              psiami=psia(rm(m,i,ic),it)
              psiapi2=psia(rp2(m,i,ic),it)
              psiami2=psia(rm2(m,i,ic),it)

              v(m,i)=v(m,i)+(-d1b12*(psiapi2-psiami2)
     &        + d2b3*(psiapi-psiami))/eps

   45         d2=d2-d1b12*((psiami2-psiac)+(psiapi2-psiac))
     &             +d4b3 *((psiami -psiac)+(psiapi -psiac))
   50   continue
      endif

      d2=d2/eps**2

   60 continue

c Warning: c1_jas6 below needs changing now that we have different
c ones for en and ee, but since ijas=6  is never used, do not bother.
c     if(ijas.eq.6) then
c       term=1/(c1_jas6*scalek(iwf))
c       usum=term*usum
c       d2=term*d2
c       do 70 i=1,nelec
c         do 70 k=1,ndim
c  70       v(k,i)=term*v(k,i)
c     endif

      value=usum

      return
      end

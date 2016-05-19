      subroutine orbitals_pw(x,orb,dorb,ddorb)
! Written by Cyrus Umrigar
! Calculate pw orbitals, gradient and laplacian.
! isortg could be used to map g-vectors from iv to ig and
! isortk could be used to map k-vectors.
! At present it is assumed that both g- and k-vectors are in the correct order.
! Note: computation time could be reduced by recognizing symmetry related k-pts
! and calculating the k-independent part of the orbital just once.

      use coefs_mod
      use const_mod
      use dim_mod
      use periodic_mod
      use pworbital_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),orb(nelec,*),dorb(3,nelec,*),ddorb(nelec,*)
      dimension dcos_rp(3),dsin_rm(3),dcos_ip(3),dsin_im(3) &
!    &,cos_g(MELEC,NGVECX),sin_g(MELEC,NGVECX),dcos_g(3,MELEC,NGVECX),dsin_g(3,MELEC,NGVECX)
!    &,ddcos_g(MELEC,NGVECX),ddsin_g(MELEC,NGVECX)
!    &,cos_k(MELEC,nkvec),sin_k(MELEC,nkvec),dcos_k(3,MELEC,nkvec),dsin_k(3,MELEC,nkvec)
!    &,ddcos_k(MELEC,nkvec),ddsin_k(MELEC,nkvec)
     &,cos_g(ngvec),sin_g(ngvec),dcos_g(3,ngvec),dsin_g(3,ngvec) &
     &,ddcos_g(ngvec),ddsin_g(ngvec) &
     &,cos_k(nkvec),sin_k(nkvec),dcos_k(3,nkvec),dsin_k(3,nkvec) &
     &,ddcos_k(nkvec),ddsin_k(nkvec),dterm1(3),dterm2(3)

!     do 5 iorb=1,norb
!       do 5 iel=1,nelec
!         orb(iel,iorb)=0
!         ddorb(iel,iorb)=0
!         do 5 k=1,ndim
!   5       dorb(k,iel,iorb)=0

      do 130 iel=1,nelec

! compute cos(g.r), sin(g.r) and derivatives
!     call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x,nelec,ng1d,cos_g,sin_g
      call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x(1,iel),ng1d,cos_g,sin_g &
     &,dcos_g,dsin_g,ddcos_g,ddsin_g,rkvec_shift)

!     write(6,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
!    &cos_g(1,1),sin_g(1,1),(dcos_g(k,1,1),k=1,ndim),(dsin_g(k,1,1),k=1,ndim),ddcos_g(1,1),ddsin_g(1,1)
!     write(6,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
!    &cos_g(1,2),sin_g(1,2),(dcos_g(k,1,2),k=1,ndim),(dsin_g(k,1,2),k=1,ndim),ddcos_g(1,2),ddsin_g(1,2)

! compute cos(k.r), sin(k.r) and derivatives
!     call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x,nelec,ng1d_sim,cos_k,sin_k
      call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x(1,iel),ng1d_sim,cos_k,sin_k &
     &,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)

      if(ipr.ge.4) then
        write(6,'(''rkvec'',9f9.5)') ((rkvec(k,ikvec),k=1,ndim),ikvec=1,nkvec)
        write(6,'(''cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k='',10f9.4)') &
     &  (cos_k(ik),sin_k(ik),(dcos_k(k,ik),k=1,ndim),(dsin_k(k,ik),k=1,ndim),ddcos_k(ik),ddsin_k(ik),ik=1,nkvec)
      endif

        iorb=0
        jorb=0
        do 130 ikvec=1,nkvec
          do 130 iband=1,nband(ikvec)
            jorb=jorb+1

            cos_rp=0
            sin_rm=0
            cos_ip=0
            sin_im=0
            do 60 k=1,ndim
              dcos_rp(k)=0
              dsin_rm(k)=0
              dcos_ip(k)=0
   60         dsin_im(k)=0
            ddcos_rp=0
            ddsin_rm=0
            ddcos_ip=0
            ddsin_im=0
!           do 80 iv=2,ngorb(ikvec)
!             ig=isortg(iv,ikvec)
            do 80 iv=2,ngvec_orb
              ig=iv
              cos_rp=cos_rp+cos_g(ig)*c_rp(iv,jorb)
              sin_rm=sin_rm+sin_g(ig)*c_rm(iv,jorb)
              cos_ip=cos_ip+cos_g(ig)*c_ip(iv,jorb)
              sin_im=sin_im+sin_g(ig)*c_im(iv,jorb)
!             write(6,'(''iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)'',3i5,9d12.4)')
!    & iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)
              do 70 k=1,ndim
                dcos_rp(k)=dcos_rp(k)+dcos_g(k,ig)*c_rp(iv,jorb)
                dsin_rm(k)=dsin_rm(k)+dsin_g(k,ig)*c_rm(iv,jorb)
                dcos_ip(k)=dcos_ip(k)+dcos_g(k,ig)*c_ip(iv,jorb)
   70           dsin_im(k)=dsin_im(k)+dsin_g(k,ig)*c_im(iv,jorb)
              ddcos_rp=ddcos_rp+ddcos_g(ig)*c_rp(iv,jorb)
              ddsin_rm=ddsin_rm+ddsin_g(ig)*c_rm(iv,jorb)
              ddcos_ip=ddcos_ip+ddcos_g(ig)*c_ip(iv,jorb)
   80         ddsin_im=ddsin_im+ddsin_g(ig)*c_im(iv,jorb)

!           write(6,'(''dcos_k(k,ikvec),dsin_k(k,ikvec),dcos_rp(k),dsin_rm(k),dsin_im(k),dcos_ip(k)'',30f9.5)')
!    &(dcos_k(k,ikvec),k=1,ndim),(dsin_k(k,ikvec),k=1,ndim),(dcos_rp(k),k=1,ndim),(dsin_rm(k),k=1,ndim),(dsin_im(k),k=1,ndim),(dcos_ip(k),k=1,ndim
!    &)

            term1=c_rp(1,jorb)+cos_rp-sin_im
            term2=c_ip(1,jorb)+cos_ip+sin_rm
            ddterm1=ddcos_rp-ddsin_im
            ddterm2=ddcos_ip+ddsin_rm
            do 85 k=1,ndim
              dterm1(k)=dcos_rp(k)-dsin_im(k)
   85         dterm2(k)=dcos_ip(k)+dsin_rm(k)

! Calculate psi_+ orbital if there are 2 indep states or if the + state is the one kept.
            if(k_inv(ikvec).eq.2. .or. ireal_imag(iorb+1).eq.1) then

            iorb=iorb+1

            if(ipr.ge.4) write(6,'(''1iorb,ireal_imag(iorb)'',9i5)') iorb,ireal_imag(iorb)

!           write(6,'(''ikvec,iel,iorb,term1,term2''3i5,9f10.5)')
!    &      ikvec,iel,iorb,term1,term2

            orb(iel,iorb)=cos_k(ikvec)*term1 &
     &                   -sin_k(ikvec)*term2
            if(ipr.ge.4) write(6,'(''1orb(iel,iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb) &
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)') &
     & iel,iorb,orb(iel,iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im

            do 90 k=1,ndim
   90         dorb(k,iel,iorb)=dcos_k(k,ikvec)*term1 &
     &                        -dsin_k(k,ikvec)*term2 &
     &                        +cos_k(ikvec)*dterm1(k) &
     &                        -sin_k(ikvec)*dterm2(k)

            ddorb(iel,iorb)=ddcos_k(ikvec)*term1 &
     &                     -ddsin_k(ikvec)*term2 &
     &                     +cos_k(ikvec)*ddterm1 &
     &                     -sin_k(ikvec)*ddterm2
!      write(6,'(''ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm='',9f9.4)')
!    &ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm

!           write(6,'(''orb'',2i5,9d12.4)') iel,iorb,orb(iel,iorb),(dorb(k,iel,iorb),k=1,ndim)
            do 100 k=1,ndim
  100         ddorb(iel,iorb)=ddorb(iel,iorb) &
     &                       +2*(dcos_k(k,ikvec)*dterm1(k) &
     &                          -dsin_k(k,ikvec)*dterm2(k))

            if(k_inv(ikvec).eq.1 .or. iorb.eq.norb) goto 130
            endif

! Calculate psi_- orbital if there are 2 indep states or if the - state is the one kept.
            iorb=iorb+1


            if(ipr.ge.4) write(6,'(''2iorb,ireal_imag(iorb)'',9i5)') iorb,ireal_imag(iorb)

!           write(6,'(''ikvec,iel,iorb,term1,term2''3i5,9f10.5)')
!    &      ikvec,iel,iorb,term1,term2

            orb(iel,iorb)=cos_k(ikvec)*term2 &
     &                   +sin_k(ikvec)*term1
            if(ipr.ge.4) write(6,'(''2orb(iel,iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb) &
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)') &
     & iel,iorb,orb(iel,iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im

            do 110 k=1,ndim
  110         dorb(k,iel,iorb)=dcos_k(k,ikvec)*term2 &
     &                        +dsin_k(k,ikvec)*term1 &
     &                        +cos_k(ikvec)*dterm2(k) &
     &                        +sin_k(ikvec)*dterm1(k)

            ddorb(iel,iorb)=ddcos_k(ikvec)*term2 &
     &                     +ddsin_k(ikvec)*term1 &
     &                     +cos_k(ikvec)*ddterm2 &
     &                     +sin_k(ikvec)*ddterm1
!           write(6,'(''orb2'',2i5,9d12.4)') iel,iorb,orb(iel,iorb),(dorb(k,iel,iorb),k=1,ndim)
            do 120 k=1,ndim
  120         ddorb(iel,iorb)=ddorb(iel,iorb) &
     &                       +2*(dcos_k(k,ikvec)*dterm2(k) &
     &                          +dsin_k(k,ikvec)*dterm1(k))
!           endif
!           endif

  130       continue

      if(ipr.ge.4) write(6,'(''orbitals_pw:'',i4,'' electrons placed in'',i4,'' orbitals'')') nelec,iorb

      return
      end
!-----------------------------------------------------------------------

      subroutine orbitals_pwe(iel,x,orb)
! Written by Cyrus Umrigar
! Calculate pw orbitals.
! isortg could be used to map g-vectors from iv to ig and
! isortk could be used to map k-vectors.
! At present it is assumed that both g- and k-vectors are in the correct order.

      use coefs_mod
      use const_mod
      use dim_mod
      use periodic_mod
      use pworbital_mod
      implicit real*8(a-h,o-z)

      dimension x(3),orb(*)
!     dimension dcos_rp(3),dsin_rm(3),dcos_ip(3),dsin_im(3)
!    &,cos_g(MELEC,NGVECX),sin_g(MELEC,NGVECX),dcos_g(3,MELEC,NGVECX),dsin_g(3,MELEC,NGVECX)
!    &,ddcos_g(MELEC,NGVECX),ddsin_g(MELEC,NGVECX)
!    &,cos_k(MELEC,nkvec),sin_k(MELEC,nkvec),dcos_k(3,MELEC,nkvec),dsin_k(3,MELEC,nkvec)
!    &,ddcos_k(MELEC,nkvec),ddsin_k(MELEC,nkvec)
      dimension cos_g(ngvec),sin_g(ngvec),dcos_g(3,ngvec),dsin_g(3,ngvec) &
     &,ddcos_g(ngvec),ddsin_g(ngvec) &
     &,cos_k(nkvec),sin_k(nkvec),dcos_k(3,nkvec),dsin_k(3,nkvec) &
     &,ddcos_k(nkvec),ddsin_k(nkvec)

!     do 5 iorb=1,norb
!   5     orb(iorb)=0

!     do 130 iel=1,nelec

! compute cos(g.r), sin(g.r) and derivatives
!     call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x,nelec,ng1d,cos_g,sin_g
!    &,dcos_g,dsin_g,ddcos_g,ddsin_g,rkvec_shift,0)
      call cossin_psi_g(glatt,gnorm,igmult,ngnorm_orb,gvec,igvec,ngvec_orb,x,ng1d,cos_g,sin_g &
     &,dcos_g,dsin_g,ddcos_g,ddsin_g,rkvec_shift)

!     write(6,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
!    &cos_g(1,1),sin_g(1,1),(dcos_g(k,1,1),k=1,ndim),(dsin_g(k,1,1),k=1,ndim),ddcos_g(1,1),ddsin_g(1,1)
!     write(6,'(''cos_g,sin_g,dcos_g,dsin_g,ddcos_g,ddsin_g='',30f9.4)')
!    &cos_g(1,2),sin_g(1,2),(dcos_g(k,1,2),k=1,ndim),(dsin_g(k,1,2),k=1,ndim),ddcos_g(1,2),ddsin_g(1,2)

! compute cos(k.r), sin(k.r) and derivatives
!     call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x,nelec,ng1d_sim,cos_k,sin_k
      call cossin_psi_k(glatt_sim,rknorm,rkvec,kvec,nkvec,x,ng1d_sim,cos_k,sin_k &
     &,dcos_k,dsin_k,ddcos_k,ddsin_k,rkvec_shift)

!     write(6,'(''cos_k,sin_k,dcos_k,dsin_k,ddcos_k,ddsin_k='',30f9.4)')
!    &cos_k(1,1),sin_k(1,1),(dcos_k(k,1,1),k=1,ndim),(dsin_k(k,1,1),k=1,ndim),ddcos_k(1,1),ddsin_k(1,1)


        iorb=0
        jorb=0
        do 130 ikvec=1,nkvec
          do 130 iband=1,nband(ikvec)
            jorb=jorb+1

            cos_rp=0
            sin_rm=0
            cos_ip=0
            sin_im=0
!           do 80 iv=2,ngorb(ikvec)
!             ig=isortg(iv,ikvec)
            do 80 iv=2,ngvec_orb
              ig=iv
              cos_rp=cos_rp+cos_g(ig)*c_rp(iv,jorb)
              sin_rm=sin_rm+sin_g(ig)*c_rm(iv,jorb)
              cos_ip=cos_ip+cos_g(ig)*c_ip(iv,jorb)
              sin_im=sin_im+sin_g(ig)*c_im(iv,jorb)
!             write(6,'(''iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)'',3i5,9d12.4)')
!    & iel,ig,jorb,cos_rp,cos_g(ig),c_rp(iv,jorb)
   80         continue

!           write(6,'(''dcos_k(k,iel,ikvec),dsin_k(k,iel,ikvec),dcos_rp(k),dsin_rm(k),dsin_im(k),dcos_ip(k)'',30f9.5)')
!    &(dcos_k(k,iel,ikvec),k=1,ndim),(dsin_k(k,iel,ikvec),k=1,ndim),(dcos_rp(k),k=1,ndim),(dsin_rm(k),k=1,ndim),(dsin_im(k),k=1,ndim),(dcos_ip(k),k=1,ndim
!    &)

            term1=c_rp(1,jorb)+cos_rp-sin_im
            term2=c_ip(1,jorb)+cos_ip+sin_rm

!           write(6,'(''iorb,term1,term2,cos_k(ikvec)*term1-sin_k(ikvec)*term2,cos_k(ikvec)*term2+sin_k(ikvec)*term1'',i5,9d12.4)')
!    & iorb,term1,term2,cos_k(ikvec)*term1-sin_k(ikvec)*term2,cos_k(ikvec)*term2+sin_k(ikvec)*term1


! Calculate psi_+ orbital if there are 2 indep states or if the + state is the one kept.
            if(k_inv(ikvec).eq.2. .or. ireal_imag(iorb+1).eq.1) then

            iorb=iorb+1

            orb(iorb)=cos_k(ikvec)*term1 &
     &               -sin_k(ikvec)*term2
            if(ipr.ge.4) then
            write(6,'(''1x='',3f9.5)') (x(k),k=1,ndim)
            write(6,'(''21orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb) &
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)') &
     & iel,iorb,orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im
            endif

!      write(6,'(''ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm='',9f9.4)')
!    &ddcos_k(ikvec),ddsin_k(ikvec),cos_k(ikvec),ddcos_rp,sin_k(ikvec),ddsin_rm

!           write(6,'(''orb'',2i5,9d12.4)') iel,iorb,orb(iorb),(dorb(k,iel,iorb),k=1,ndim)

            if(k_inv(ikvec).eq.1 .or. iorb.eq.norb) goto 130

            endif

! Calculate psi_- orbital if there are 2 indep states or if the - state is the one kept.
            iorb=iorb+1

            orb(iorb)=cos_k(ikvec)*term2 &
     &               +sin_k(ikvec)*term1
            if(ipr.ge.4) write(6,'(''22orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb) &
     &,cos_rp,cos_ip,sin_rm,sin_im='',2i5,20d12.4)') &
     & iel,iorb,orb(iorb),cos_k(ikvec),sin_k(ikvec),c_rp(1,jorb),c_ip(1,jorb),cos_rp,cos_ip,sin_rm,sin_im

!           endif
!           endif

  130       continue

      if(ipr.ge.4) write(6,'(''orbitals_pwe:'',i4,'' electrons placed in'',i4,'' orbitals'')') nelec,iorb

      return
      end

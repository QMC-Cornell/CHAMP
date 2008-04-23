Module pjas_setup_mod 
  use all_tools_mod 
  implicit none

  complex (dp)                           :: Z_0= (0._dp,0._dp)
  complex (dp)                           :: Z_i= (0._dp,1._dp)
  
  
  integer                                :: nstar_en=0
  integer                                :: nstar_ee =0

  logical                                :: not_include_zero = .true.
  integer                                :: nsym_crys
  real(dp), allocatable                  :: symrel (:,:,:)
  real(dp), allocatable                  :: tnons (:,:)

  integer                                :: nsym_crys_sim
  real(dp), allocatable                  :: symrel_sim (:,:,:)
  real(dp), allocatable                  :: tnons_sim (:,:)
  
  integer                                :: nstar , nbasis_pw
  integer                                :: kvec_pjasen (3)
  real(dp),allocatable                   :: rkv(:,:)
  integer,allocatable                    :: kv(:,:)
  real(dp), allocatable                  :: sk3(:)
  complex(dp),allocatable                :: phase(:)
  integer,allocatable                    :: mstar(:),istar(:),fstar(:)
  
  integer                                :: kvec_pjasee (3)
  integer                                :: nstar_sim , nbasis_pw_sim
  real(dp), allocatable                  :: rkv_sim (:,:), sk3_sim (:)
  integer,allocatable                    :: kv_sim(:,:)
  integer, allocatable                   :: istar_sim (:), mstar_sim (:), fstar_sim (:)
  complex(dp),allocatable                :: phase_sim(:)
 
  logical                                :: inversion = .true. ! defaul has inversion symm

contains


!!! prim_sim =1 for primitive cell and 2 for simulation cell
!!! input
!!            nsym_crys number of symmetries
!!!           ecut cuttoff energy for planewave generation
!!!           avec and bvec              :: crystall lattice and reciprocal lattice vectors
!!!           wsvol                      :: cell volume
!!!           symrel and tnons: rotational and translational symmetry operators
!!! output
!!!          nstar       number of stars 
!!!          nbasis_pw   number of planewaves in nstars
!!! also     generates
!!!  sk3: magnitude of each star, phase: phase for each vector, rkv lattice vector
!!!  mstar number of planewaves in each star, istar: index for the begining of each star
  
  subroutine create_pw_lattice(prim_sim,nstar,nbasis_pw,nsym_crys,ecut,avec,bvec,wsvol,symrel,tnons)
!---------------------------------------------------------------------------
! Description : generate the stars.                                                              
!                                                                           
! Created     : W. A. Al-Saidi, June 2007                          
!---------------------------------------------------------------------------
    implicit none 
    real(dp),parameter                   :: deltlg =1.0e-12_dp
    integer, parameter                   :: mlignd=2500
    integer                              :: prim_sim
    real(dp)                             :: symrel(3,3,nsym_crys), tnons(3,nsym_crys) !! symmetry input
    integer                              :: n3d_wis,nsym_crys, nstar, nbasis_pw, ntstrd
    real(dp)                             :: wsvol
    real(dp)                             :: avec(3,3), bvec (3,3), bij(3,3)
    real(dp)                             :: ecut
    complex(dp)                          :: cc
    real(dp)                             :: cx
    logical                              :: iprint,lg1,lg2,lg3
    integer                              :: hi,low
    real(dp)                             :: bi(3),bj(3),bk(3)
    integer                              :: mnbmax,i,j,ipas,jpas,kpas,mnm,k,l,mn,mntot,mlign,iprot,ii,mn2
    integer                              :: m1,m2,m3,i1,i2,mnbr,ntt,ir
    real(dp)                             :: g,gkmod,xli,q
    real(dp)                             :: bijk,bpi,bpj,bpk,bbi,bbj,gki,gkj,bbk,gkk
    integer                              :: ntstar,nop, istat, kk, nbasis_pw_t
    real(dp)                             :: temp , pi, twopi
    complex(dp), allocatable             :: etagv(:)
    real(dp), allocatable                :: ekin_t(:),rkv_t(:,:)
    integer, allocatable                 :: kv_t(:,:)
    real(dp), allocatable                :: resultt(:,:),res2(:)
    real(dp), allocatable                :: resmat(:,:)
    real(dp)                             :: gv(3),gij(3,3),gkt(3)
    real(dp)                             :: protyp(3)
    logical, allocatable                 :: impt(:,:)
    logical, allocatable                 :: used(:),impr(:)
    integer, allocatable                 :: indexx(:)
    real(dp), allocatable                :: sk3_t(:)
    complex(dp),allocatable              :: phase_t(:)
    integer     ,allocatable             :: mstar_t(:),istar_t(:)
    integer                              :: i3
    real(dp),allocatable                 :: sk3_tt(:), ekin_tt(:), rkv_tt (:,:)
    integer, allocatable                 :: mstar_tt(:), istar_tt(:), kv_tt(:,:)
    complex(dp), allocatable             :: phase_tt (:)
    integer                              :: j1, j2 , phase_zero, ist, r

    !
    !     calculate the metric in reciprocal space:
    !
    bij= 0 
    
    do j=1,3
       do i=1,3
          do  k=1,3
             bij(i,j)=bvec(k,i)*bvec(k,j)+bij(i,j)  
          enddo
       enddo
    enddo


    pi =4*atan(1.0_dp)
    twopi= 2*4*atan(1.0_dp)

    nbasis_pw=wsvol*ecut**1.5/(6*pi*pi)

    nbasis_pw= 8 * nbasis_pw 
    ntstrd=9*Nbasis_pw
    n3d_wis= ntstrd


    allocate(sk3_t(n3d_wis)) 
    allocate(phase_t(ntstrd)) 
    allocate(mstar_t(n3d_wis))
    allocate(istar_t(n3d_wis))
    
    allocate (ekin_t(ntstrd),rkv_t(3,ntstrd),kv_t(3,ntstrd))
    
    allocate(etagv(nsym_crys)) 
    allocate(resultt(3,ntstrd),res2(ntstrd))
    allocate(resmat(3,mlignd))
    allocate(impt(mlignd,nsym_crys))
    allocate(used(mlignd),impr(mlignd))
    allocate(indexx(ntstrd))


    phase_t=Z_0;
    mstar_t=0
    istar_t=0

    bi=bvec(:,1)
    bj=bvec(:,2)
    bk=bvec(:,3)
    
    mnbmax=ntstrd

    do j=1,3 
       gij(j,j)=one
       do  i=1,3
	  if (i.ne.j) gij(i,j)=zero
       enddo
    enddo
    !
    !     set up all g-vectors corresponding to a sphere of energy=ecut (in ry)
    !     this would be gmax=sqrt(ecut), but we also need q=g-g', so the
    !     the radius is twice this value:
    g=2.d0*sqrt(ecut)

    xli=1.0e-6_dp
    gkmod=g
   
    bijk=bi(1)*(bj(2)*bk(3)-bj(3)*bk(2)) +&
	 bi(2)*(bj(3)*bk(1)-bj(1)*bk(3)) +&
	 bi(3)*(bj(1)*bk(2)-bj(2)*bk(1))
    bijk=abs(bijk)
    bpi=gdot ( bi,bi,gij )
    bpj=gdot ( bj,bj,gij )
    bpk=gdot ( bk,bk,gij )
    bbi=bijk/sqrt( bpj*bpk - gdot( bj,bk,gij )**2 )
    bbj=bijk/sqrt( bpk*bpi - gdot( bk,bi,gij )**2 )
    bbk=bijk/sqrt( bpi*bpj - gdot( bi,bj,gij )**2 )
    gki=int ( gkmod / bbi + 1. )
    gkj=int ( gkmod / bbj + 1. )
    gkk=int ( gkmod / bbk + 1. )
    !
    !
    ipas=gki * 2. + 1.
    jpas=gkj * 2. + 1.
    kpas=gkk * 2. + 1.
    mnm=0
    do i=1,ipas
       gv(1)=-gki+i-1
       do j=1,jpas
	  gv(2)=-gkj+j-1
	  do k=1,kpas 
	     gv(3)=-gkk+k-1
	     gkmod=sqrt(gdot(gv,gv,bij))
             if (not_include_zero .and. gkmod == 0) cycle
	     if (gkmod.gt.g ) cycle
	     mnm=mnm + 1
	     !
	     if (mnm.gt.ntstrd) then
		write ( 6,*) 'create_pw_lattice: increase star dim.: mn,ntstrd=',mnm,ntstrd
		stop 'create_pw_lattice:  increase  star dimensions'
	     endif

             ekin_t(mnm)=gkmod+deltlg*mnm
             
	     do l=1,3
		rkv_t(l,mnm)=gv(l)
	     enddo
	  enddo 
       enddo
    enddo 
    !
    !
    !     construction of stars
    !
    call hsort(mnm,ekin_t(1:mnm),indexx(1:mnm))
    do  mn=1,mnm
       ekin_t(mn)=ekin_t(mn)-mn*deltlg
    enddo

    do  mn=1,mnm
       res2(mn)=ekin_t(indexx(mn))
       do  l=1,3
	  resultt(l,mn)=rkv_t(l,indexx(mn))
       enddo
    enddo

    !
    hi=0
    nstar=0
    mntot=0
    !
1001 continue

    low=hi+1
    !
    !     define star or superstar
    !
    gkmod=res2(low)
    hi=low
17  if(hi.eq.mnm)  goto 18
    if(res2(hi+1).le.gkmod+xli)then
       hi=hi+1
       if((hi-low).gt.mlignd)then
	  print*, 'hi=',hi
	  print*, 'low=',low
	  print*, mlignd
	  write( 6,'('' create_pw_lattice:mlignd'')')
          write(6,*) 'mlignd: check your symmetry operations and your lattice'
	  stop
       endif
       goto 17
    endif
18  continue
    mlign=hi-low+1

    do mn=1,mlign
       used(mn)=.false.
    enddo
    iprot=low
    !
    !     operate on prototype with group operations
203 continue
    !
    do j=1,3
       protyp(j) = resultt(j,iprot)
    enddo
    do mn=1,mlign
       impr(mn)=.true.
    enddo

    nop = nsym_crys 

    do ii=1,nop

       do j=1,3
	  gkt(j)=0.d0
	  do k=1,3
	     gkt(j)=gkt(j)+symrel(j,k,ii)*protyp(k)
	  enddo
       enddo
       lg1=.true.

       do mn= 1, mlign
	  mn2=mn+low-1
	  !
          temp =abs(gkt(1)-resultt(1,mn2))+abs(gkt(2)-resultt(2,mn2))+&
               abs(gkt(3)-resultt(3,mn2))

	  if(temp .gt. xli ) then
	     impt(mn,ii)=.true.
	  else


	     impt(mn,ii)=.false.
	     impr(mn)=.false.
	     lg1=.false.
	     cx= 0.d0
	     do j=1,3
		cx=cx-gkt(j)*tnons (j, ii) 
	     enddo
	     cx=cx*twopi
	     etagv(ii)=cmplx(cos(cx), sin(cx))
             !             if ((sin(cx)) > 0.00001) then 
             !                stop "complex phase with"
             !             endif
          endif

       enddo !mn align 
       if(lg1)then
	  write(6,*)'error in create_pw_lattice: unsymmetric superstar', nstar
          write(6,*) 'Symmetrey operations inconsistent with your lattice'
          write(6,*) 'create_pw_lattice'
          stop
       endif
       !
    enddo !! symmetry opers
    !
    mnbr=mlign
    do mn=1,mlign
       if ( impr(mn) )   mnbr=mnbr - 1
    enddo

    
    nstar=nstar + 1

    if (nstar.gt.n3d_wis)then
       write(6,'('' n3d_wis too small'')')
       write(6,*) 'n3d_wis too small'
       stop 
    endif
    mstar_t(nstar)=mnbr
    istar_t(nstar)=mntot+1
    ntt=nop/mnbr

    
    do mn=mlign,1,-1
       mn2 = mn + low - 1
       if (impr(mn).and. .not.used(mn)) iprot=mn2
    enddo 
    do mn=1,mlign
       mn2=mn+low-1
       if (impr(mn)) then 
          cycle
       endif
       !
       if(used(mn)) then
	  write(6,*)resultt(1,mn2),resultt(2,mn2),resultt(3,mn2)
	  write(*,*)'create_pw_lattice: attempt to reuse rec. lat. vector'
       endif
       used(mn)=.true.
       mntot=mntot + 1
       cc=z_0
       do ii=1,nop 
	  if(impt(mn,ii)) then 
             cycle 
          endif
	  cc=cc+etagv(ii)
       enddo 
       cc=cc/real(ntt,dp)
       phase_t(mntot)=cc
       ekin_t(mntot)=res2(mn2)**2
       sk3_t(nstar)=res2(mn2)**2
       do j=1,3
	  kv_t(j,mntot)=nint(resultt(j,mn2))
       enddo
       m1=kv_t(1,mntot)
       m2=kv_t(2,mntot)
       m3=kv_t(3,mntot)

       do j=1,3
	  rkv_t(j,mntot)=kv_t(1,mntot)*bvec(j,1)&
	       +kv_t(2,mntot)*bvec(j,2)&
	       +kv_t(3,mntot)*bvec(j,3)

       enddo

    enddo
    
    ntstar=istar_t(nstar)+mstar_t(nstar)-1

    if(ntstar.lt.hi)  goto 203
    
    if (hi.lt.mnm)   go to 1001

    do i=1,ntstar
       if (ekin_t(i).gt.ecut) then
          nbasis_pw_t=i-1
          go to 9999
       endif
    enddo
    write(6,*) 'create_pw_lattice: not enough plane waves for q-vector'
    stop 
9999 continue

!!! For simulation phase, remove  phase
    if (prim_sim==2) then 
       do j= 1, ntstar
          phase_t(j) = 1 
       enddo
    endif


!!! remove those stars which give zero contribution due to symmetry
    
    write(*,*) "nstar ",  nstar 

    allocate(mstar_tt(ntstar), istar_tt(ntstar),sk3_tt(ntstar))
    allocate(ekin_tt (ntstar), kv_tt(3,ntstar), phase_tt(ntstar),rkv_tt(3,ntstar))

    k=0
    istar_tt(1)=1
    do i=1,nstar
       i1=istar_t (i)
       i2=i1+mstar_t (i)-1
       phase_zero= 0

       do  j=i1,i2
          if (abs(phase_t(j)) > 0.00000001) then
             phase_zero = 1 
             exit
          endif
       enddo

!!$       !! debug
!!$       write(6,*) i, i1, i2, mstar_t (i), phase_zero 
!!$       do  j=i1,i2
!!$          write (6,5) j,ekin_t(j),kv_t(1,j),kv_t(2,j),kv_t(3,j),phase_t(j),&
!!$               rkv_t(1,j),rkv_t(2,j),rkv_t(3,j)
!!$       enddo
!!$       !! debug

       
       if (phase_zero == 1 ) then 
          k=k+1
          sk3_tt(k)= sk3_t (i)
          mstar_tt(k) = mstar_t (i)
          if (k==1) then 
             rkv_tt(:,1:i2-i1+1)= rkv_t (:,i1:i2)
             kv_tt (:,1:i2-i1+1)= kv_t (:,i1:i2)
             istar_tt(1)= 1 
             phase_tt(1:i2-i1+1)= phase_t (i1:i2)
             ekin_tt(1:i2-i1+1)= ekin_t (i1:i2)
          else
             j1= istar_tt(k-1)+ mstar_tt (k-1) 
             rkv_tt(:,j1:i2-i1+j1)= rkv_t (:,i1:i2)
             kv_tt (:,j1:i2-i1+j1)= kv_t (:,i1:i2)
             istar_tt(k)= j1 
             phase_tt(j1:i2-i1+j1) = phase_t (i1:i2)
             ekin_tt(j1:i2-i1+j1) =  ekin_t (i1:i2)
          endif
       else

          write(6,*) "star is zero by symmetry", i 
          !! debug
          write(6,*) i, i1, i2, mstar_t (i) 
          do  j=i1,i2
             write (6,5) j,ekin_t(j),kv_t(1,j),kv_t(2,j),kv_t(3,j),phase_t(j),&
                  rkv_t(1,j),rkv_t(2,j),rkv_t(3,j)
          enddo
          !! debug

       endif
    enddo
    nstar=k 
    
!!! pick only stars within this cutoff
    do i=1,ntstar
       if (sk3_tt(i).gt.ecut) then
          nstar=i-1
          exit
       endif
    enddo

    !! note checking is done on the orignal vectors before 
    !! reducing them by symmetry.
    !! check that the phases 
    do i=1,nstar
       i1=istar_tt (i)
       i2=i1+ (mstar_tt (i))/2 - 1 
       i3=i1+mstar_tt (i)-1
       kk=0
       do j=i1,i2
          if (maxval(real(kv_tt (:,j)+kv_tt(:,i3-kk)))  > 0.0000001) then 
             !             write(*,*) "kv_tt1 ",kv_tt(:,j)
             !             write(*,*) "kv_tt2 ",kv_tt(:,i3-kk)
             !             write(*,*) "diff is ", maxval(real(kv_tt (:,j)+kv_tt(:,i3-kk)))
             write(*,*) "you should use generalized stars which include G and -G "
             write(*,*) "add inversion symmetry to your list of symmetries"
             stop
          endif
          if (abs(phase_tt(j)-conjg(phase_tt(i3-kk))) > 0.0000001) then 
             write(*,*) "star kv_tt bar{kv_tt}=", i, kv_tt (:,j), kv_tt (:,i3-kk)
             write(*,*) "phase1 ",phase_tt(j)
             write(*,*) "phase2 ",phase_tt(i3-kk)
             write(*,*) "phases should be conjugate to each other"
             stop
          endif
          kk=kk+1
       enddo
    enddo


    if (prim_sim==1 .and. nstar_en > 0) then 
       if (nstar_en > nstar) nstar_en = nstar 
       nstar= nstar_en 
    endif

    
    if (prim_sim==2 .and. nstar_ee > 0) then 
       if (nstar_ee > nstar) nstar_ee = nstar 
       nstar= nstar_ee 
    endif


    i= nstar
    nbasis_pw_t =istar_tt (i) + mstar_tt (i)-1


    !! redimension   
    nbasis_pw = nbasis_pw_t

    !!include only one half 
    nbasis_pw = nbasis_pw/2 

    !!global variables for primitive cell
    if (prim_sim==1) then 
       
       allocate(rkv(3,nbasis_pw), stat=istat)
       !       allocate(kv(3,nbasis_pw), stat=istat)
       allocate (phase(nbasis_pw),mstar(nstar),istar(nstar))
       allocate (fstar(nstar))
       allocate (kv(3,nbasis_pw), stat=istat)

!!! remove the planewaves related by inversion symmetry
       r=1
       do i=1,nstar
          i1=istar_tt (i)
          i2=i1+ (mstar_tt (i)-1)/2
          do j=i1,i2
             rkv(:,r)= rkv_tt(:,j)
             kv(:,r)= kv_tt(:,j)
             phase (r)= phase_tt (j)
             ekin_tt (r)= ekin_tt (j) !! destroy ekin_tt 
             r=r+1
          enddo
          if (i==1) then 
             istar(1)= 1
             mstar(1)= mstar_tt(1)/2
          else
             istar(i)= istar(i-1) + mstar(i-1) 
             mstar(i)= mstar_tt(i)/2
          endif
       enddo
       do k=1, 3
          kvec_pjasen (k)= int (maxval (abs (real (kv (k,:)))))
       enddo
       write(6,*) " kvec_pjasen = ",  kvec_pjasen
       

       if (.not. inversion) then 
          !! double the size of sk3 
          allocate (sk3(2*nstar))
          k=1
          do i=1, nstar
             sk3 (k)= sk3_tt (i)
             sk3 (k+1)= sk3_tt (i)
             k=k+2 
          enddo
       else
          allocate (sk3(nstar))
          sk3= sk3_tt (1:nstar)
       endif

       
       do ist=1, nstar 
          i1 = istar (ist)
          i2 = i1+ mstar (ist) - 1 
          fstar(ist)= i2
       enddo

      write(6,'(a,I4,a,I4,a)') " Number of planewaves included in e-n periodic Jastrow= ", nbasis_pw, " in ", nstar, " stars"

      do i=1,nstar
          write(6,6) i,mstar(i),istar(i),sk3(i)
          i1=istar (i)
          i2=fstar (i)
          do j=i1,i2
             write (6,5) j,ekin_tt(j),kv(1,j),kv(2,j),kv(3,j),phase(j),&
                  rkv(1,j),rkv(2,j),rkv(3,j)
          enddo
       enddo

    else

       allocate (rkv_sim(3,nbasis_pw), stat=istat)
       allocate (sk3_sim(nstar),mstar_sim(nstar),istar_sim(nstar))
       allocate (fstar_sim (nstar))
       allocate (kv_sim(3,1:nbasis_pw))
       
!!! remove the planewaves related by inversion symmetry
       r=1
       do i=1,nstar
          i1=istar_tt (i)
          i2=i1+ (mstar_tt (i)-1)/2
          do j=i1,i2
             rkv_sim(:,r)= rkv_tt(:,j)
             kv_sim(:,r) = kv_tt (:,j)
             ekin_tt (r)= ekin_tt (j) !! destroy ekin_tt 
             phase_tt (r)= phase_tt (j) !! destroy phase_tt 
             r=r+1
          enddo
          if (i==1) then 
             istar_sim(1)= 1
             mstar_sim(1)= mstar_tt(1)/2
          else
             istar_sim(i)= istar_sim(i-1) + mstar_sim(i-1) 
             mstar_sim(i)= mstar_tt(i)/2
          endif
       enddo

       do k=1, 3
          kvec_pjasee (k)= int (maxval (abs (real (kv_sim (k,:)))))
       enddo
       write(6,*) " kvec_pjasee = ",  kvec_pjasee
       


       sk3_sim = sk3_tt (1:nstar)
       do ist=1, nstar_sim 
          i1=istar_sim (ist)
          i2=i1+ mstar_sim (ist) - 1 
          fstar_sim(ist)= i2
       enddo
       
       write(6,'(a,I4,a,I4,a)') " Number of planewaves included in e-e periodic Jastrow= ", nbasis_pw, " in ", nstar, " stars"

       do i=1,nstar
          write(6,6) i,mstar_sim(i),istar_sim(i),sk3_sim(i)
          i1=istar_sim (i)
          i2=fstar_sim (i)
          do  j=i1,i2
             write (6,5) j,ekin_tt(j),kv_sim(1,j),kv_sim(2,j),kv_sim(3,j),phase_tt(j),&
                  rkv_sim(1,j),rkv_sim(2,j),rkv_sim(3,j)
          enddo
       enddo
       
    endif
    
    
!!$    do i=1,nstar
!!$       write(6,6) i,mstar_tt(i),istar_tt(i),sk3_tt(i)
!!$       i1=istar_tt (i)
!!$       i2=i1+mstar_tt (i)-1
!!$       do  j=i1,i2
!!$          write (6,5) j,ekin_tt(j),kv_tt(1,j),kv_tt(2,j),kv_tt(3,j),phase_tt(j),&
!!$               rkv_tt(1,j),rkv_tt(2,j),rkv_tt(3,j)
!!$       enddo
!!$    enddo



!!$
!!$    write(101+prim_sim,*) nstar
!!$    write(101+prim_sim,'(120I5.2)') (mstar_tt(i),i=1,nstar)
!!$    write(101+prim_sim,'(120I5.2)') (istar_tt(i),i=1,nstar)
!!$    
!!$    do j=1,nbasis_pw 
!!$       write(101+prim_sim,'(12F14.8)') ekin_tt(j), rkv_tt(1,j),rkv_tt(2,j),rkv_tt(3,j), phase_tt(k)
!!$    enddo
    
    
    deallocate (etagv,resultt,res2,resmat,impt,used,indexx)
    deallocate (sk3_t,phase_t, mstar_t,istar_t,ekin_t,rkv_t,kv_t)
    

    deallocate(mstar_tt, istar_tt ,sk3_tt, ekin_tt , kv_tt , phase_tt,rkv_tt)
    

2   format(10x,i7,' plane waves in ',i6,' stars',/,&
	 10x,'gki gkj gkk = ',3f5.1)
5   format(15x, i5,f15.8,3i5,2f10.5,10x,3f10.5)
6   format(/3i5,5x,f15.8/)
7   format (3(/4x,3f10.6))
8   format (8f10.6)
    
  contains

    
    function  gdot(a,b,gij)
!---------------------------------------------------------------------------
! Description : do                                                             
!     *****
!     *****   dot product of two vectors with metric gij
!     *****
!     ***********************************************************
! Created     : W. A. Al-Saidi, June 2007                          
!---------------------------------------------------------------------------
      implicit none 
      real(dp)                           :: a(3),b(3),gij(3,3),gdot
      integer                            :: i,j
      gdot=zero
      do  j=1,3
	 do  i=1,3
	    gdot=gdot+a(i)*gij(i,j)*b(j)
	 enddo
      enddo
    end function gdot

    
    subroutine hsort(n,a,ind)
!---------------------------------------------------------------------------
! Description : sort an array                                                             
!                                                                           
! Created     : W. A. Al-Saidi, June 2007                          
!---------------------------------------------------------------------------
      implicit none 
      integer                            :: n
      real(dp)                           :: a(n),q
      integer                            :: ind(n)
      integer                            :: j,l,indxt,ir
      do j=1,n
	 ind(j)=j
      enddo
      if(n.le.1)return
      l=n/2+1
      ir=n
10    continue
      if(l.gt.1)then
	 l=l-1
	 indxt=ind(l)
	 q=a(indxt)
      else
	 indxt=ind(ir)
	 q=a(indxt)
	 ind(ir)=ind(1)
	 ir=ir-1
	 if(ir.eq.1)then
	    ind(1)=indxt
	    return
	 endif
      endif
      i=l
      j=l+l
20    if(j.le.ir)then
	 if(j.lt.ir)then
	    if(a(ind(j)).lt.a(ind(j+1)))j=j+1
	 endif
	 if(q.lt.a(ind(j)))then
	    ind(i)=ind(j)
	    i=j
	    j=j+j
	 else
	    j=ir+1
	 endif
	 goto 20
      endif
      ind(i)=indxt
      goto 10

    end subroutine hsort
   
  end subroutine create_pw_lattice
  

  
end Module pjas_setup_mod

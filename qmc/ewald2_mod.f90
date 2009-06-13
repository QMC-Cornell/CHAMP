module ewald2_mod 

! JT: put subroutines shells and sort in a module in order to be able to allocate dummy arguments

 contains

      subroutine shells(cutg,glatt,gdist,igvec,gvec,gnorm,igmult,ngvec_big,ngnorm_big,ng1d,icell)
! Written by Cyrus Umrigar

! icell = 0  primitive cell
!         1  simulation cell

      use basic_tools_mod, only: alloc
      use dim_mod
      implicit real*8(a-h,o-z)

      integer, allocatable :: igvec(:,:),igmult(:),ng1d(:)
      double precision, allocatable :: gvec(:,:),gnorm(:)
      dimension glatt(3,*),gdist(3)
      double precision, allocatable :: gnorm_tmp(:)

      call alloc ('ng1d', ng1d, 3)
      do 1 k=1,ndim
    1   ng1d(k)=int(cutg/gdist(k))

! JT: do first loop to compute ngvec_big
      cutg2=cutg**2
      ngvec_big=0
!     do 10 i1=-ng1d(1),ng1d(1)
      do 10 i1=0,ng1d(1)
        if(i1.ne.0) then
          i2min=-ng1d(2)
         else
          i2min=0
        endif
!       do 10 i2=-ng1d(2),ng1d(2)
        do 10 i2=i2min,ng1d(2)
          if(i2.ne.0.or.i1.ne.0) then
            i3min=-ng1d(3)
           else
            i3min=0
          endif
!         do 10 i3=-ng1d(3),ng1d(3)
          do 10 i3=i3min,ng1d(3)

            gx=i1*glatt(1,1)+i2*glatt(1,2)+i3*glatt(1,3)
            gy=i1*glatt(2,1)+i2*glatt(2,2)+i3*glatt(2,3)
            gz=i1*glatt(3,1)+i2*glatt(3,2)+i3*glatt(3,3)

            glen2=gx*gx+gy*gy+gz*gz

            if(glen2.le.cutg2) then
              ngvec_big=ngvec_big+1
!JT              if(icell.eq.0 .and. ngvec_big.gt.NGVEC_BIGX) then
!JT                stop 'ngvec_big > NGVEC_BIGX in shells'
!JT               elseif(icell.eq.1 .and. ngvec_big.gt.NGVEC_SIM_BIGX) then
!JT                stop 'ngvec_big > NGVEC_SIM_BIGX in shells'
!JT              endif

!JT              igvec(1,ngvec_big)=i1
!JT              igvec(2,ngvec_big)=i2
!JT              igvec(3,ngvec_big)=i3

!JT              gvec(1,ngvec_big)=gx
!JT              gvec(2,ngvec_big)=gy
!JT              gvec(3,ngvec_big)=gz

!JT              gnorm_tmp(ngvec_big)=dsqrt(glen2)
            endif
   10 continue


      call alloc ('igvec', igvec, 3, ngvec_big)
      call alloc ('gvec', gvec, 3, ngvec_big)
      call alloc ('gnorm_tmp', gnorm_tmp, ngvec_big)

! JT: redo same loop now that the arrays are allocated
      cutg2=cutg**2
      ngvec_big=0
!     do 11 i1=-ng1d(1),ng1d(1)
      do 11 i1=0,ng1d(1)
        if(i1.ne.0) then
          i2min=-ng1d(2)
         else
          i2min=0
        endif
!       do 11 i2=-ng1d(2),ng1d(2)
        do 11 i2=i2min,ng1d(2)
          if(i2.ne.0.or.i1.ne.0) then
            i3min=-ng1d(3)
           else
            i3min=0
          endif
!         do 11 i3=-ng1d(3),ng1d(3)
          do 11 i3=i3min,ng1d(3)

            gx=i1*glatt(1,1)+i2*glatt(1,2)+i3*glatt(1,3)
            gy=i1*glatt(2,1)+i2*glatt(2,2)+i3*glatt(2,3)
            gz=i1*glatt(3,1)+i2*glatt(3,2)+i3*glatt(3,3)

            glen2=gx*gx+gy*gy+gz*gz

            if(glen2.le.cutg2) then
              ngvec_big=ngvec_big+1
!JT              if(icell.eq.0 .and. ngvec_big.gt.NGVEC_BIGX) then
!JT                stop 'ngvec_big > NGVEC_BIGX in shells'
!JT               elseif(icell.eq.1 .and. ngvec_big.gt.NGVEC_SIM_BIGX) then
!JT                stop 'ngvec_big > NGVEC_SIM_BIGX in shells'
!JT              endif

              igvec(1,ngvec_big)=i1
              igvec(2,ngvec_big)=i2
              igvec(3,ngvec_big)=i3

              gvec(1,ngvec_big)=gx
              gvec(2,ngvec_big)=gy
              gvec(3,ngvec_big)=gz

              gnorm_tmp(ngvec_big)=dsqrt(glen2)
            endif
   11 continue

      call sort(igvec,gvec,gnorm_tmp,gnorm,igmult,ngvec_big,ngnorm_big,icell)

      return
      end subroutine
!-----------------------------------------------------------------------

      subroutine sort(igvec,gvec,gnorm_tmp,gnorm,igmult,ngvec_big,ngnorm_big,icell)
      use basic_tools_mod, only: alloc
      use dim_mod
      implicit real*8(a-h,o-z)
! Written by Cyrus Umrigar
! Use Shell-Metzger sort to put g-vectors in some standard order, so that
! the order they appear in is independent of cutg_sim_big.

      parameter(eps=1.d-12)

      integer, allocatable :: igmult(:)
      double precision, allocatable :: gnorm(:)
      dimension igvec(3,*),gvec(3,*),gnorm_tmp(*)

!     cost(igv1,igv2,igv3,gn)=igv3+10.d0**4*igv2+10.d0**8*igv1+10.d0**12*gn
!     cost(igv1,igv2,igv3,gn)=((igv3+10.d0**4*igv2)+10.d0**8*igv1)+10.d0**14*gn
      cost(igv1,igv2,igv3,gn)=((igv3+10.d0**2*igv2)+10.d0**4*igv1)+10.d0**14*gn

      lognb2=int(dlog(dfloat(ngvec_big))/dlog(2.d0)+1.d-14)
      m=ngvec_big
      do 20 nn=1,lognb2
        m=m/2
        k=ngvec_big-m
        do 20 j=1,k
          do 10 i=j,1,-m
            l=i+m
!           if(gnorm_tmp(l).gt.gnorm_tmp(i)-eps) goto 20
            if(cost(igvec(1,l),igvec(2,l),igvec(3,l),gnorm_tmp(l)).gt.cost(igvec(1,i),igvec(2,i),igvec(3,i),gnorm_tmp(i))) goto 20
            t=gnorm_tmp(i)
            gnorm_tmp(i)=gnorm_tmp(l)
            gnorm_tmp(l)=t
            do 10 k=1,ndim
              it=igvec(k,i)
              igvec(k,i)=igvec(k,l)
              igvec(k,l)=it
              t=gvec(k,i)
              gvec(k,i)=gvec(k,l)
   10         gvec(k,l)=t
   20     continue

! JT: do first loop to compute ngvec_big
! figure out the multiplicities and convert gnorm from being ngvec_big long to being ngnorm_big long
      ngnorm_big=1
      icount=0
      do 30 i=2,ngvec_big
        icount=icount+1
        if(gnorm_tmp(i)-gnorm_tmp(i-1).gt.eps) then
!JT          igmult(ngnorm_big)=icount
!JT          gnorm(ngnorm_big)=gnorm_tmp(i-1)
          ngnorm_big=ngnorm_big+1
!JT          if(icell.eq.0 .and. ngnorm_big.gt.NGNORM_BIGX) then
!JT            write(6,'(''ngnorm_big,NGNORM_BIGX='',2i8)') ngnorm_big,NGNORM_BIGX
!JT            stop 'ngnorm_big > NGNORM_BIGX in sort'
!JT           elseif(icell.eq.1 .and. ngnorm_big.gt.NGNORM_SIM_BIGX) then
!JT            write(6,'(''ngnorm_sim_big,NGNORM_SIM_BIGX='',2i8)') ngnorm_big,NGNORM_SIM_BIGX
!JT            stop 'ngnorm_sim_big > NGNORM_SIM_BIGX in sort'
!JT          endif
          icount=0
        endif
   30 continue
      call alloc ('igmult', igmult, ngnorm_big) 
      call alloc ('gnorm', gnorm, ngnorm_big) 
! JT: redo same loop now that the arrays are allocated
! figure out the multiplicities and convert gnorm from being ngvec_big long to being ngnorm_big long
      ngnorm_big=1
      icount=0
      do 31 i=2,ngvec_big
        icount=icount+1
        if(gnorm_tmp(i)-gnorm_tmp(i-1).gt.eps) then
          igmult(ngnorm_big)=icount
          gnorm(ngnorm_big)=gnorm_tmp(i-1)
          ngnorm_big=ngnorm_big+1
!JT          if(icell.eq.0 .and. ngnorm_big.gt.NGNORM_BIGX) then
!JT            write(6,'(''ngnorm_big,NGNORM_BIGX='',2i8)') ngnorm_big,NGNORM_BIGX
!JT            stop 'ngnorm_big > NGNORM_BIGX in sort'
!JT           elseif(icell.eq.1 .and. ngnorm_big.gt.NGNORM_SIM_BIGX) then
!JT            write(6,'(''ngnorm_sim_big,NGNORM_SIM_BIGX='',2i8)') ngnorm_big,NGNORM_SIM_BIGX
!JT            stop 'ngnorm_sim_big > NGNORM_SIM_BIGX in sort'
!JT          endif
          icount=0
        endif
   31 continue
      igmult(ngnorm_big)=icount+1
      gnorm(ngnorm_big)=gnorm_tmp(ngvec_big)

! Check that looping over norms and using multiplicities we get the right number
      icheck=0
      do 40 i=1,ngnorm_big
   40   icheck=icheck+igmult(i)
      if(icheck.ne.ngvec_big) stop 'problem in sort'

! Since vectors are sorted not just by norm but by 3 components also, make sure that
! norms are monotonically ascending within a tolerance.  I do not stop unless the violation
! is large because for lattice vectors that are very nearly equal, the norms may violate this.
      do 50 i=2,ngvec_big
        if(gnorm_tmp(i)+1.d-12.lt.gnorm_tmp(i-1)) then
          write(6,'(''Warning: gnorm_tmps are not monotonic, i,gnorm(i-1),gnorm(i),diff='',i5,3d12.4)')i,gnorm_tmp(i-1),gnorm_tmp(i),gnorm_tmp(i)-gnorm_tmp(i-1)
        if(gnorm_tmp(i)+1.d-8.lt.gnorm_tmp(i-1)) stop 'gnorm_tmps are not monotonic'
        endif
   50 continue
      do 60 i=2,ngnorm_big
        if(gnorm(i)+1.d-12.lt.gnorm(i-1)) then
          write(6,'(''Warning: gnorms are not monotonic, i,gnorm(i-1),gnorm(i),diff='',i5,3d12.4)')i,gnorm(i-1),gnorm(i),gnorm(i)-gnorm(i-1)
        if(gnorm(i)+1.d-8.lt.gnorm(i-1)) stop 'gnorms are not monotonic'
        endif
   60 continue

!     if(icell.eq.0) then
!       write(6,'(''shells and vectors in primitive cell'')')
!      else
!       write(6,'(''shells and vectors in simulation cell'')')
!     endif
!     j=0
!     do 100 i=1,ngnorm_big
!       do 100 im=1,igmult(i)
!         j=j+1
! 100     write(6,'(''CHECK '',i5,i4,i3,x,3i3,9f10.4)')
!    &    j,i,igmult(i),(igvec(k,j),k=1,ndim),gnorm(i),(gvec(k,j),k=1,ndim)
!     write(6,*)

      return
      end subroutine

end module ewald2_mod 

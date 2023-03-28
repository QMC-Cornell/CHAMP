
      subroutine cjastrow(cv,d2,div,absval)
! written by A.D.Guclu jul2004
! calculates composite fermion jastrow factor (vortices, e^k) related quantities.

      use const_mod
      use dim_mod
      use distance_mod
      implicit real*8(a-h,o-z)

      complex(dpc) cv

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot

      dimension cv(3,*),div(*)

      nn=2*nv
      if(idot.eq.2) nn=nn+1     ! for laughlin wave function

! first calculate e^K = PI_{j<i} (zji)^(2nv)
! with z=x-iy, zji=zj-zi (or zi-zj, doesnot matter..)

      ij=0
!      z=(1.d0,0.d0)
      absval=1.d0
      do 40 i=2,nelec
        do 40 j=1,i-1
          ij=ij+1
!         zij=(-rvec(1,ij),rvec(2,ij))
          absval=absval*r_ee(ij)
   40 continue
      absval=absval**nn

! now calculate the velocities. we only need the real parts. divergences are
! not needed for now...

      do 60 i1=1,nelec
        sumx=0.d0
        sumy=0.d0
        do 50 i2=1,nelec
! confusing part (reverify this later):
          if(i1.lt.i2) then
            ij=((i2-1)*(i2-2))/2+i1
            r2_ee=r_ee(ij)*r_ee(ij)
            sumx=sumx+rvec_ee(1,ij)/r2_ee
            sumy=sumy+rvec_ee(2,ij)/r2_ee
           elseif(i1.gt.i2) then
            ij=((i1-1)*(i1-2))/2+i2
            r2_ee=r_ee(ij)*r_ee(ij)
            sumx=sumx-rvec_ee(1,ij)/r2_ee
            sumy=sumy-rvec_ee(2,ij)/r2_ee
          endif
 50     continue
        cv(1,i1)=dcmplx(-nn*sumx,-nn*sumy)
        cv(2,i1)=dcmplx(-nn*sumy,nn*sumx)
 60   continue

      d2=0.d0
      do i=1,nelec
        div(i)=0.d0
      enddo

! there are pieces due to the gaussian tail of laughlin functions and proj. cfs:
      if(idot.ge.2) then
        sumr2=0.d0
        do 80 i=1,nelec
          sumr2=sumr2+r_en(i,1)*r_en(i,1)
          do 70 k=1,ndim
            cv(k,i)=cv(k,i)-dcmplx(we*rvec_en(k,i,1),0.d0)
 70       continue
          div(i)=-2*we
          d2=d2+div(i)
 80     continue
        absval=absval*dexp(-sumr2*we/2)
      endif

      return
      end


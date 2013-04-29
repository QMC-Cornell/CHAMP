      subroutine jastrow1(x,v,d2,div_vj,value)
! Written by Cyrus Umrigar

      use dets_mod
      use const_mod
      use dim_mod
      use wfsec_mod
      use distance_mod
      use jaspar1_mod
      implicit real*8(a-h,o-z)


!JT      parameter(zero=0.d0,one=1.d0,two=2.d0,third=1.d0/3.d0)




      dimension x(3,*),v(3,*),div_vj(*)

      fsum=zero
      if(nelec.lt.2) goto 60

      ij=0
      do 40 i=2,nelec
        im1=i-1
        do 40 j=1,im1
          ij=ij+1
          if(i.le.nup .or. j.gt.nup) then
            if(ndim.eq.3) then
              cjasa=cjas1(iwf)/2
             elseif(ndim.eq.2) then
              cjasa=cjas1(iwf)/3
            endif
           else
            cjasa=cjas1(iwf)
          endif

          rij=r_ee(ij)
          riji=one/rij

          bot=one/(one+cjas2(iwf)*rij)
          fsum=fsum+cjasa*rij*bot
          term=cjasa*riji*bot**2
          v(1,i)=v(1,i)+term*rvec_ee(1,ij)
          v(2,i)=v(2,i)+term*rvec_ee(2,ij)
          v(3,i)=v(3,i)+term*rvec_ee(3,ij)
          v(1,j)=v(1,j)-term*rvec_ee(1,ij)
          v(2,j)=v(2,j)-term*rvec_ee(2,ij)
          v(3,j)=v(3,j)-term*rvec_ee(3,ij)
          div_vj(i)=div_vj(i)+2*cjasa*riji*bot**3
          div_vj(j)=div_vj(j)+2*cjasa*riji*bot**3
          d2=d2+two*(ndim-1)*cjasa*riji*bot**3

   40     continue

   60 continue

      value=fsum

      return
      end

      subroutine jastrow1(x,v,d2,div_vj,value)
c Written by Cyrus Umrigar

      use dets_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'force.h'

!JT      parameter(zero=0.d0,one=1.d0,two=2.d0,third=1.d0/3.d0)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /jaspar1/ cjas1(MWF),cjas2(MWF)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

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

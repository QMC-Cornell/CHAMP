      subroutine splfit_bas(r,irb,ic,iwf,wfv,ider)
c Written by Claudia Filippi
c get spline_fit at r of basis fn irb of center ic and force iwf
c 1st and 2nd derivs also calculated if ider=1.
c ref0 is the nearest point smaller than r and ref1 the nearest one larger.

      use numbas_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'numbas.h'
!JT      include 'force.h'

      parameter(NCOEF=5)

!JT      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
!JT     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
!JT     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
      common /numexp/ce(NCOEF,MRWF,MCTYPE,MFORCE),ae(2,MRWF,MCTYPE,MFORCE)

      dimension wfv(3)

      if(igrid(ic).eq.1)then
        xr=(r-r0_bas(ic))/exp_h_bas(ic)+1
        jx=int(xr)
        ref0=r0_bas(ic)+exp_h_bas(ic)*(jx-1)
        ref1=ref0+exp_h_bas(ic)
        delh=exp_h_bas(ic)
       elseif(igrid(ic).eq.2)then
c       dlrad0=dlog(r0_bas(ic))
        h_bas=dlog(exp_h_bas(ic))
c       xr=(dlog(r)-dlrad0)/h_bas+1
        xr=dlog(r/r0_bas(ic))/h_bas+1
        jx=int(xr)
        ref0=r0_bas(ic)*exp_h_bas(ic)**(jx-1)
        ref1=ref0*exp_h_bas(ic)
        delh=ref1-ref0
       elseif(igrid(ic).eq.3)then
        h_bas=dlog(exp_h_bas(ic))
        xr=(dlog((r+r0_bas(ic))/r0_bas(ic)))/h_bas+1.d0
        jx=int(xr)
        ref0=r0_bas(ic)*exp_h_bas(ic)**(jx-1)-r0_bas(ic)
        ref1=(ref0+r0_bas(ic))*exp_h_bas(ic)-r0_bas(ic)
        delh=ref1-ref0
      endif

      if(jx.lt.1) then
c use small radius polynomial

        wfv(1)=ce(1,irb,ic,iwf)
        do 10 ii=2,ncoef
  10      wfv(1)=wfv(1)+ce(ii,irb,ic,iwf)*r**(ii-1)
        if(ider.eq.1) then
          wfv(2)=0.d0
          wfv(3)=0.d0
          do 20 ii=2,ncoef
            wfv(2)=wfv(2)+(ii-1)*ce(ii,irb,ic,iwf)*r**(ii-2)
  20        wfv(3)=wfv(3)+(ii-1)*(ii-2)*ce(ii,irb,ic,iwf)*r**(ii-3)
        endif

       elseif(jx.ge.nr(ic)) then
c use large radius exponential

        wfv(1)=ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*r)
        if(ider.eq.1) then
          wfv(2)=-ae(2,irb,ic,iwf)*wfv(1)
          wfv(3)=-ae(2,irb,ic,iwf)*wfv(2)
        endif

       else
c cubic spline interpolation

        bb=(r-ref0)/delh
        aa=(ref1-r)/delh
        cc=aa*(aa**2-1.d0)*delh**2/6.d0
        dd=bb*(bb**2-1.d0)*delh**2/6.d0
        wfv(1)=aa*rwf(jx,irb,ic,iwf)+bb*rwf(jx+1,irb,ic,iwf)+
     &         cc*d2rwf(jx,irb,ic,iwf)+dd*d2rwf(jx+1,irb,ic,iwf)

        if(ider.eq.1) then
          wfv(2)=(rwf(jx+1,irb,ic,iwf)-rwf(jx,irb,ic,iwf))/delh+
     &    (-(3.d0*aa**2-1.d0)*d2rwf(jx,irb,ic,iwf)
     &     +(3.d0*bb**2-1.d0)*d2rwf(jx+1,irb,ic,iwf))*delh/6.d0
          wfv(3)=aa*d2rwf(jx,irb,ic,iwf)+bb*d2rwf(jx+1,irb,ic,iwf)
        endif
      endif

      return
      end

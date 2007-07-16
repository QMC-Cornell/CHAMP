cDate of last modification: Tue May 14 09:21:24 EDT 1996
      function gaushn()
c gaushn: purpose: generate N(0,1) distributed random numbers
      implicit real*8(a-h,o-z)
      dimension ran(2)
      parameter (one=1,two=2)
      data iset/0/
      save gset,iset
      if (iset.eq.0) then
1       call ransr(ran,2)
        v1=two*ran(1)-one
        v2=two*ran(2)-one
        r=v1**2+v2**2
        if(r.ge.1..or.r.eq.0.) go to 1
        fac=sqrt(-2.*log(r)/r)
        gset=v1*fac
        gaushn=v2*fac
        iset=1
      else
        gaushn=gset
        iset=0
      endif
      return
      end

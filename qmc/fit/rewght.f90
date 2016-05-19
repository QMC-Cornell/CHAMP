      subroutine rewght(diff)
! Written by Cyrus Umrigar
! Calculate the weights due to the fact that we have sampled from
! one probability distribution, but want to calculate the expectation
! value with respect to another.
      use constants_mod
      use contr2_mod
      use confg_mod
      implicit real*8(a-h,o-z)

      parameter (biglog=500.d0)

      dimension diff(*)

      wtmax=iabs(mod(irewgt,100))

! If the normalizations of the old and new wavefns. is very different then
! when we take the exp it could be outside the range of the smallest and
! largest numbers that can be represented on the computer, so in that case
! renormalize the old wavefn. and recompute weights.
      av=0
      do 20 i=1,ndata
   20   av=av+psij(i)+dlog(dabs(psid(i)))-psio(i)
      av=av/ndata
      if(abs(av).gt.biglog) then
        write(6,'(''Warning: normalization of old wavefn. is being reset by '',d12.4)') exp(av)
        do 30 i=1,ndata
   30     psio(i)=psio(i)+av
      endif

      wghtsm=zero
      do 100 i=1,ndata
        exponent=2*(psij(i)+dlog(dabs(psid(i)))-psio(i))
        if(exponent.lt.biglog) then
          wght(i)=dvpdv(i)*exp(exponent)
         else
          wght(i)=exp(biglog)
        endif
        wghtsm=wghtsm+wght(i)
  100 continue

! Put max lim. on wght and calculate new wghtsm
      do 120 iter=1,5
        term=dfloat(ndata)/wghtsm
        eav=zero
        wghtsm=zero
        iflag=0
        do 110 i=1,ndata
          wght(i)=wght(i)*term
          if(wght(i).gt.wtmax+one) then
            wght(i)=wtmax
            iflag=1
          endif
          eav=eav+uwdiff(i)
  110     wghtsm=wghtsm+wght(i)
        if(iflag.eq.0) goto 130
  120 continue

! scale diffs by wght
  130 term=dfloat(ndata)/wghtsm
      eav=eguess+eav/wghtsm
      do 134 i=1,ndata
        wght(i)=wght(i)*term
        if(iaver.ge.1) diff(i)=diff(i)+eguess-eav
        diff(i)=diff(i)*dsqrt(wght(i))
  134 continue
!     wghtsm=ndata

      return
      end

      subroutine rewght(diff)
c Written by Cyrus Umrigar
c Calculate the weights due to the fact that we have sampled from
c one probability distribution, but want to calculate the expectation
c value with respect to another.
      implicit real*8(a-h,o-z)

      include '../vmc/vmc.h'

      parameter (zero=0.d0,one=1.d0,biglog=500.d0)

      include 'fit.h'

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata

      dimension diff(*)

      wtmax=iabs(mod(irewgt,100))

c If the normalizations of the old and new wavefns. is very different then
c when we take the exp it could be outside the range of the smallest and
c largest numbers that can be represented on the computer, so in that case
c renormalize the old wavefn. and recompute weights.
      av=0
      do 20 i=1,ndata
   20   av=av+psij(i)+dlog(dabs(psid(i)))-psio(i)
      av=av/ndata
      if(abs(av).gt.biglog) then
        write(6,'(''Warning: normalization of old wavefn. is being reset'')')
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

c Put max lim. on wght and calculate new wghtsm
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

c scale diffs by wght
  130 term=dfloat(ndata)/wghtsm
      eav=eguess+eav/wghtsm
      do 134 i=1,ndata
        wght(i)=wght(i)*term
        if(iaver.ge.1) diff(i)=diff(i)+eguess-eav
        diff(i)=diff(i)*dsqrt(wght(i))
  134 continue
c     wghtsm=ndata

      return
      end

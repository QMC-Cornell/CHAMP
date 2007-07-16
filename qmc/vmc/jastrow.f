      subroutine jastrow(x,v,d2,div_vj,value)
c Written by Cyrus Umrigar

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      include 'pseudo.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      dimension x(3,*),v(3,*)
     &,div_vj(MELEC)

      do 10 i=1,nelec
        div_vj(i)=0
        v(1,i)=0
        v(2,i)=0
   10   v(3,i)=0
      d2=0

      if(ijas.eq.1) then
	call jastrow1(x,v,d2,div_vj,value)
       elseif(ijas.eq.2) then
        call jastrow2(x,v,d2,div_vj,value)
       elseif(ijas.eq.3) then
        call jastrow3(x,v,d2,div_vj,value)
       elseif(ijas.ge.4.and.ijas.le.6) then
        call jastrow4(x,v,d2,div_vj,value)
      endif

      return
      end

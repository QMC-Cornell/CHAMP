      subroutine deriv_jastrow(x,v,d2,div_vj,value)
c Written by Claudia Filippi

      implicit real*8(a-h,o-z)
      parameter (zero=0.d0)

      include '../vmc/vmc.h'
      include '../vmc/force.h'
      include '../vmc/pseudo.h'

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      dimension x(3,*),v(3,*),div_vj(MELEC)

      do 10 i=1,nelec
        v(1,i)=zero
        v(2,i)=zero
   10   v(3,i)=zero
      d2=zero

      if(ijas.eq.1) then
        call jastrow1(x,v,d2,div_vj,value)
       elseif(ijas.eq.2) then
        call jastrow2(x,v,d2,div_vj,value)
       elseif(ijas.eq.3) then
        call deriv_jastrow3(x,v,d2,value)
       elseif(ijas.ge.4.and.ijas.le.6) then
        call deriv_jastrow4(x,v,d2,value)
      endif

      return
      end

      subroutine jastrowe(iel,x,v,value)
c Written by Claudia Filippi by modifying jastrow

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      dimension x(3,*),v(3,*)

      do 10 i=1,nelec
        v(1,i)=0
        v(2,i)=0
   10   v(3,i)=0

      if(ijas.eq.3) then
        call jastrow3e(iel,x,v,value)
       elseif(ijas.ge.4.and.ijas.le.6) then
        call jastrow4e(iel,x,v,value)
       else
        stop 'single-electron move only for jastrow3 and jastrow4'
      endif

      return
      end

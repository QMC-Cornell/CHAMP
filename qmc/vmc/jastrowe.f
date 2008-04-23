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

!! WAS
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
      common /pjase/ ido_pjasen, ido_pjasee, ido_pjas
!!!

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
!! WAS 
      if (ido_pjas .eq. 1) then 
         call pjas_jas_e_interface (iel, x, rvec_ee, v, value)
      endif
!!! 
      return
      end

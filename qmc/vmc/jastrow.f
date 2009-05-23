      subroutine jastrow(x,v,d2,div_vj,value)
c Written by Cyrus Umrigar

      use const_mod
      use contr2_mod
      use distance_mod
      use pjase_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'force.h'
!     ! WAS
!JT      common /pjase/ ido_pjasen, ido_pjasee, ido_pjas
!JT      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
!!!
!JT      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
!JT     &,ifock,i3body,irewgt,iaver,istrch
!JT     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
!JT      include 'pseudo.h'

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

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

!     ! WAS
!!!   add contrib due to long range jastrow
      if (ido_pjas .eq. 1) then
         call pjas_jas_interface (x, rvec_ee, v, d2, div_vj, value)
      endif


      return
      end

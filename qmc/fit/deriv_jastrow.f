      subroutine deriv_jastrow(x,v,d2,div_vj,value)
c Written by Claudia Filippi

      implicit real*8(a-h,o-z)
      parameter (zero=0.d0)

      include '../vmc/vmc.h'
      include '../vmc/force.h'
      include '../vmc/pseudo.h'
!     ! WAS
      common /pjase/ ido_pjasen, ido_pjasee, ido_pjas
      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)
!!!
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
!! WAS
!!! add contrib due to long range jastrow
      if (ido_pjas .eq. 1) then

c$$$         write(55,'(a,100F10.6)') "bv=",  v(:,1:nelec)
c$$$         write(55,'(a,10F10.6)') "bd2and value ",  d2, value
c$$$
         call pjas_deriv_jas_interface (x,rvec_ee, v, d2, div_vj, value)
c$$$
c$$$         write(55,'(a,100F10.6)') "v=",  v(:,1:nelec)
c$$$         write(55,'(a,10F10.6)') "d2and value ",  d2, value
      endif
!!!




      return
      end

      subroutine deriv_jastrow(x,v,d2,div_vj,value)
! Written by Claudia Filippi
      use constants_mod
      use const_mod
      use contr2_mod
      use distance_mod
      use pjase_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),v(3,*),div_vj(nelec)

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

!$$$         write(55,'(a,100F10.6)') "bv=",  v(:,1:nelec)
!$$$         write(55,'(a,10F10.6)') "bd2and value ",  d2, value
!$$$
         call pjas_deriv_jas_interface (x,rvec_ee, v, d2, div_vj, value)
!$$$
!$$$         write(55,'(a,100F10.6)') "v=",  v(:,1:nelec)
!$$$         write(55,'(a,10F10.6)') "d2and value ",  d2, value
      endif
!!!




      return
      end

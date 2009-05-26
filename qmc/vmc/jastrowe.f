      subroutine jastrowe(iel,x,v,value)
c Written by Claudia Filippi by modifying jastrow

      use const_mod
      use contr2_mod
      use distance_mod
      use pjase_mod
      implicit real*8(a-h,o-z)



!! WAS
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

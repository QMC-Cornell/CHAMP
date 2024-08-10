      subroutine jastrowe(iel,x,v,value)
! Written by Claudia Filippi by modifying jastrow

      use const_mod
      use contr2_mod
      use distance_mod
      use pjase_mod
      use contrldmc_mod, only: l_tau_diffusion, l_psi_approx
      implicit real*8(a-h,o-z)

      dimension x(3,*),v(3,*)

      do 10 i=1,nelec
        v(1,i)=0
        v(2,i)=0
   10   v(3,i)=0

      if(ijas.eq.3) then
        call jastrow3e(iel,x,v,value)
       elseif(ijas.ge.4.and.ijas.le.6) then
        if (l_tau_diffusion.OR.l_psi_approx) then
          call jastrow4e_lapl(iel,x,v,value)
        else
          call jastrow4e(iel,x,v,value)
        endif
       else
        stop 'single-electron move only for jastrow3 and jastrow4'
      endif

!! WAS
      if (ido_pjas .eq. 1) then
         call pjas_jas_e_interface (iel, x, rvec_ee, v, value)
      endif

      return
      end

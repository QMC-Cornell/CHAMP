      subroutine basis_norm(iwf,iflag)
! Written by Cyrus Umrigar
! Set normalization of basis fns.
! In 3d:
! Norm of radial part: ((2*zeta)^{2n+1}/(2n)!)^{1/2} for Slaters (n>0).
!                      (2*(2*zeta)^(n+1/2)/Gamma(n+1/2))^{1/2} for gaussians
! where                Gamma(n+1/2) = integral {t^{n+1/2} e^-t dt} = (2n-1)!! sqrt(pi)/2^n for gaussians
!                      Gamma(1/2)=sqrt(pi), Gamma(a+1)=a*Gamma(a), Gamma(a)=(a-1)!
! obtained by Integrate[r^2 (r^{n-1} Exp[-zeta r])^2,{r,0,Infinity}] for Slaters
! obtained by Integrate[r^2 (r^{n-1} Exp[-zeta r^2])^2,{r,0,Infinity}] for Gaussians
! Norm of angular part: ((2*l+1)/(4*pi))^{1/2}.
! In 2d:
! Norm of radial part: ((2*zeta)^{2n}/(2n-1)!)^{1/2}.
! Norm of angular part: (min(m+1,2)/(2*pi))^{1/2}.
! If numr =0 or -1 we are using analytic basis functions and we use normalization
!                  for angular and radial parts.  The -1 is just to tell it to order
!                  basis functions by all the s's first, then all the p's etc.
!                  instead of 1s,2s,2p,3s,3p,3d,...
!         =1       we are using numerical basis functions and we use normalization
!                  for angular part only
!                  The check has now been changed to iwrwf2(ib).le.nrbas_analytical(ict) 14 Oct 09
! Whether one is using Slater or gaussian basis fns. used to be inputted by having
! n1s,n2s etc. be either > 0 or < 0.  Now use which_analytical_basis
      use all_tools_mod
      use const_mod
      use control_mod
      use atom_mod, only : iwctype
      use orbitals_mod
      use coefs_mod
      use dim_mod
      use numbas_mod
      use basis1_mod
      use basis_mod, only : which_analytical_basis, norm_gauss_slat_exp_1 !fp
      use basis2_mod
      use basisnorm_mod
      use contr2_mod
      use contrl_per_mod
      implicit real*8(a-h,o-z)
! anorm stored for reuse in fit.  Since iwf=1 in fit, we omit iwf dependence.

      call alloc ('anorm', anorm, nbasis)
      call alloc ('n_bas', n_bas, nbasis)
      call alloc ('l_bas', l_bas, nbasis)
      call alloc ('m_bas', m_bas, nbasis)

      do 20 ib=1,nbasis
        n=n_bas(ib)
        if(ndim.eq.3) then
          ict=ictype_basis(ib)
          l=l_bas(ib)
          if(iwrwf2(ib).le.nrbas_analytical(ict)) then
            select case (trim(which_analytical_basis)) !fp
             case ('slater')    !fp
                anorm(ib)=sqrt((2*zex(ib,iwf))**(2*n+1)*(2*l+1)/(fact(2*n)*4*pi))
             case ('gaussian')  !fp
                n1=abs(n)
                anorm(ib)=sqrt(2*(2*zex(ib,iwf))**(n1+0.5d0)*(2*l+1)/(gamma1(n1)*4*pi))
             case ('gauss-slater') !fp
                anorm(ib) = norm_gauss_slat_exp_1(n) * sqrt((2*l+1)/(4*pi)) * zex(ib,iwf)**(n+0.5d0)   !fp
             case default
                write(6,*) 'basis_norm: Allowed basis types are slater gaussian gauss-slater!'
                stop 'basis_norm: Allowed basis types are slater gaussian gauss-slater!'
            end select         !fp
           else
            anorm(ib)=sqrt((2*l+1)/(4*pi))
          endif
         elseif(ndim.eq.2) then
          m=m_bas(ib)
          if(numr.le.0 .and. ibasis.lt.4) then
!         if(iwrwf2(ib).le.nrbas_analytical(ict) .and. ibasis.lt.4) then
            anorm(ib)=sqrt((2*zex(ib,iwf))**(2*n)*min(abs(m)+1,2)/(fact(2*n-1)*2*pi))
! The following change is not necessary at this time and has not been tested.
           elseif(numr.le.0 .and. (ibasis.ge.4 .and. ibasis.le.7)) then
!          elseif(iwrwf2(ib).le.nrbas_analytical(ict) .and. (ibasis.ge.4 .and. ibasis.le.7)) then
            anorm(ib)=dsqrt(1/pi)
           else
! Warning: temporarily commented out diff norm for m=0
!           anorm(ib)=sqrt(min(abs(m)+1,2)/(2*pi))
            anorm(ib)=sqrt(1/pi)
          endif
        endif
        if(iflag.eq.1) then
          do iorb=1,norb
            coef(ib,iorb,iwf)=coef(ib,iorb,iwf)*anorm(ib)
          enddo
        endif
   20 continue
      if(ipr.ge.0) write(6,'(''anorm='',20f10.6)') (anorm(ib),ib=1,nbasis)

      call object_modified ('anorm')   !JT
!      write(6,*) 'JT icusp=',icusp
!      if (icusp.eq.1) then
!       call ie
!       call equiv_bas
!       call impose_cusp_en_orb_occ
!      endif
      call object_modified ('coef')    !JT

      return
      end
!-----------------------------------------------------------------------

      function fact(n)
      implicit real*8(a-h,o-z)

      fact=1
      do 10 i=2,n
   10   fact=fact*i
      return
      end
!-----------------------------------------------------------------------

      function gamma1(n)
! W. Al-Saidi 4/11/07
! Used for norm of 3D Guassians, gamma1(n) = gamma(n+1/2) = integral {t^{n+1/2} e^-t dt} = (2n-1)!! sqrt(pi)/2^n
      implicit real*8(a-h,o-z)

      spi=1.77245385090552d0 ! sqrt{pi}

      if(n.eq.0) then
        gamma1= spi
       elseif(n.eq.1) then
         gamma1= spi/2
       elseif(n.eq.2) then
         gamma1=3*spi/4
       elseif(n.eq.3) then
         gamma1=15*spi/8
       elseif(n.eq.4) then
         gamma1=105*spi/16
       elseif(n.eq.5) then
         gamma1=945*spi/32
       elseif(n.eq.6) then
         gamma1=10395*spi/64
       else
         write(6,'(''gamma1 is not implemented  in basis_norm.f for n='',i3)') n
         stop "gamma1 is not implemented  in basis_norm.f for this n"
      endif
      return
      end

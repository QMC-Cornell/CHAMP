!**RM(6)
!*******************
      MODULE real_spherical_harmonics
!*******************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Module to the code extension by RM
!       for the recursive construction of higher
!       real spherical harmonics.
!---------------
      use all_tools_mod !JT
! Warning: If you change lmax, you need to recompile not only basis_fns.f but also basis_fns2e.f because
! module real_spherical_harmonics is used in both
      integer,parameter   :: lmax=4
!     integer,parameter   :: lmax=6
!MS Jellium sphere
!     integer,parameter   :: lmax=12
      integer             :: nmax
      integer,allocatable :: large_q(:,:)
!MS Jellium sphere
      real(dp),allocatable:: sq_coef_a_num(:),sq_coef_a_denum(:), &
     &                       sq_coef_r_denum(:,:)
      real(dp)            :: large_p(0:lmax,-lmax:lmax),d_large_p(1:3,0:lmax,-lmax:lmax), &
     &                       r_inv_power_minus(0:lmax),x_power_of(0:lmax), &
     &                       y_power_of(0:lmax),z_power_of(0:lmax), &
     &                       coef_ylm(0:lmax,0:lmax)=1.d0

      END MODULE real_spherical_harmonics
!**EndRM(6)

      subroutine basis_fns(iel,rvec_en,r_en)

! Written by Cyrus Umrigar.  L=4 added in by John Lawson(JWL).  Recursive Ylm for all L added in by Ryo Maezono
! Calculate 3-dim localised basis functions and their derivatives
! iel = 0, for all electrons
!    != 0, for electron iel

! In input:
! n1s,n2s,...     > 0 : Slater basis
! n1s,n2s,...     < 0 : Gaussian basis
! nsa,npa,nda         : asymptotic basis
! Here:
! n_bas2(irb,ict) > 0 : Slater basis
!                 < 0 : Gaussian basis
!                 = 0 : asymptotic basis

!     For notes on basis functions in QMC see fp notes (contain slater and gaussian 3-d)

      use all_tools_mod
      use constants_mod
      use const_mod
      use dim_mod
      use real_spherical_harmonics
      use atom_mod
      use basis1_mod
      use basis2_mod
      use basis_mod, only : which_analytical_basis !fp
      use numbas_mod
      use wfsec_mod
      use phifun_mod
      use contr_ylm_mod
      implicit real*8(a-h,o-z)
      real(dp) :: aux1
      integer :: itemp1

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

      dimension wfv(3,MRWF),xc(3),th(0:ML_BAS,0:ML_BAS),dth(3,0:ML_BAS,0:ML_BAS),ph(-ML_BAS:ML_BAS),dph(3,-ML_BAS:ML_BAS)

! Here we have additional normalization factors beyond those in basis_norm, viz., sqrt((2*l+1)/(4*pi))
! The additional normalization factors for d,f,g are Sqrt of 1/4, 3, 3/4; 1/4, 3/8, 15/4, 5/8; 1/64, 10/16, 5/16, 70/16, 35/64.
! JWL added cg0-4
      data cd0,cd1,cd2,cf0,cf1,cf2,cf3,cg0,cg1,cg2,cg3,cg4/ &
     &0.5d0,1.73205080756888d0,0.866025403784439d0,0.5d0,0.612372435695794d0,1.93649167310371d0,0.790569415042095d0, &
     &0.125d0,0.790569415042095d0,0.559016994374947d0,2.09165006633519d0,0.739509972887452d0/


      ider=1

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

      do 40 ie=nelec1,nelec2
        ib=0
        do 40 ic=1,ncent
          ict=iwctype(ic)

          if(r_en(ie,ic).eq.0.d0) then
            write(6,'(''Warning: basis_fns: r_en(ie,ic) = 0.d0'')')
            r_en(ie,ic)=1.d-99
            rvec_en(1,ie,ic)=1.d-99
          endif
          xc(1)=rvec_en(1,ie,ic)
          xc(2)=rvec_en(2,ie,ic)
          xc(3)=rvec_en(3,ie,ic)
          r=r_en(ie,ic)
          r2=r*r
          ri=1/r
          ri2=ri*ri
          ri3=ri2*ri
          ri4=ri3*ri
          ri5=ri4*ri

! Loops over analytical basis functions, then numerical.  One can optimize the exponents of the analytic ones.
          do 10 irb=1,nrbas_analytical(ict)
              n=n_bas2(irb,ict)
              rn=abs(n)
              rm3=r**(rn-3)
              rm2=rm3*r
              rm1=rm2*r
              if(n.ne.0) then
                 select case (trim(which_analytical_basis)) !fp
                 case ('slater') !fp
                    zr=zex2(irb,ict,iwf)*r
                    ex=dexp(-zr)
                    wfv(3,irb)=rm3*((rn-1)*(rn-2-zr)-(rn-1-zr)*zr)*ex
                 case ('gaussian') !fp
                    zr=2*zex2(irb,ict,iwf)*r2
                    ex=dexp(-0.5d0*zr)
                    wfv(3,irb)=rm3*((rn-1)*(rn-2)-(2*rn-1-zr)*zr)*ex
                 case ('gauss-slater') !fp
                    zr=(zex2(irb,ict,iwf)*r)**2/(1+zex2(irb,ict,iwf)*r)**2 * (2+zex2(irb,ict,iwf)*r)  !fp
!                    ex=dexp(-zr * (1+zex2(irb,ict,iwf)*r) / (2+zex2(irb,ict,iwf)*r) ) !fp
                    ex=dexp(-(zex2(irb,ict,iwf)*r)**2 / (1+zex2(irb,ict,iwf)*r) ) !fp
                    wfv(3,irb)=rm3*ex *((rn-1-zr)*(rn-2-zr)-zr*(1+2/(1+zex2(irb,ict,iwf)*r)-2/(2+zex2(irb,ict,iwf)*r))) !fp
                 case default
                    write(6,*) 'basis_fns: Allowed basis types are slater gaussian gauss-slater!'
                    stop 'basis_fns: Allowed basis types are slater gaussian gauss-slater!'
                 end select     !fp
              elseif(n.eq.0) then
!     Warning: Asymptotic and Gaussian not yet tested.
!     Asymptotic r^(rn-1)*Exp(-zeta*r), where rn=beta+1, beta=betaq/zeta-1, zeta=sqrt(-2*E_ion)?
                 write(6,*) 'basis_fns: ict=',ict
                 write(6,*) 'basis_fns: irb=',irb
                 write(6,*) 'basis_fns: n=',n
                 stop 'basis_fns: asymptotic not yet fully tested'
                 rn=betaq/zex2(irb,ict,iwf)
                 rm3=r**(rn-3)
                 rm2=rm3*r
                 rm1=rm2*r
                 wfv(3,irb)=rm3*((rn-1)*(rn-2-zr)-(rn-1-zr)*zr)*ex
              endif
              wfv(1,irb)=rm1*ex
   10         wfv(2,irb)=rm2*((rn-1)-zr)*ex

!       write(6,'(''ic,ict,nrbas(ict),wfv='',3i5,29d12.5)') ic,ict,nrbas(ict),(wfv(1,irb),irb=1,nrbas(ict))

          do 20 irb=1,nrbas_numerical(ict)
            rk=r
            call splfit_bas(rk,irb,ict,iwf,wfv(1,nrbas_analytical(ict)+irb),ider)
            if(wfv(1,nrbas_analytical(ict)+irb).eq.0.d0) wfv(1,nrbas_analytical(ict)+irb)=DBLMIN
   20     continue

       if(ipr.ge.3) write(6,'(''ic,ict,nrbas(ict),n_bas2='',3i5,29i12)') &
     &   ic,ict,nrbas(ict),(n_bas2(irb,1),irb=1,nrbas_analytical(ict))
       if(ipr.ge.3) write(6,'(''ic,ict,nrbas(ict),zex2  ='',3i5,29d12.5)') &
     &   ic,ict,nrbas(ict),(zex2(irb,1,1),irb=1,nrbas_analytical(ict))
       if(ipr.ge.3) write(6,'(''ic,ict,nrbas(ict),wfv   ='',3i5,29d12.5)') &
     &   ic,ict,nrbas(ict),(wfv(1,irb),irb=1,nrbas_analytical(ict))

!**RM(2)
          if(irecursion_ylm.eq.0)  then
            xx=xc(1)
            yy=xc(2)
            zz=xc(3)
            xx2=xx*xx
            yy2=yy*yy
            zz2=zz*zz

            xhat=xx*ri
            yhat=yy*ri
            zhat=zz*ri
          do ib2=1,nbasis_ctype(ict)
!          do 40 ib2=1,nbasis_ctype(ict)  !* Modified 'do 40' --> 'do-enddo' for
!                                         !  code extension.
!            xx=xc(1)
!            yy=xc(2)
!            zz=xc(3)
!            xx2=xx*xx
!            yy2=yy*yy
!            zz2=zz*zz
!            xhat=xx*ri
!            yhat=yy*ri
!            zhat=zz*ri
!**EndRM(2)

! Phi function

            ph(0)=1
            dph(1,0)=0
            dph(2,0)=0

            ph(1)=xx
            dph(1,1)=1
            dph(2,1)=0

            ph(-1)=yy
            dph(1,-1)=0
            dph(2,-1)=1

            ph(2)=xx2-yy2
            dph(1,2)= 2*xx
            dph(2,2)=-2*yy

            ph(-2)=2*xx*yy
            dph(1,-2)=2*yy
            dph(2,-2)=2*xx

!           ph(3)=(xx2-yy2)*xx-2*xx*yy2
            ph(3)=ph(2)*ph(1)-ph(-2)*ph(-1)
            dph(1,3)=3*ph(2)
            dph(2,3)=-3*ph(-2)

!           ph(-3)=3*xx2*y-yy2*y
            ph(-3)=ph(-2)*ph(1)+ph(-1)*ph(2)
            dph(1,-3)=3*ph(-2)
            dph(2,-3)=3*ph(2)

! JWL added l=4
            ph(4)=xx2*xx2-6*xx2*yy2+yy2*yy2
            dph(1,4)=4*xx*xx2-12*xx*yy2
            dph(2,4)=-12*xx2*yy+4*yy*yy2

            ph(-4)=4*xx*yy*(xx2-yy2)
            dph(1,-4)=12*xx2*yy-4*yy*yy2
            dph(2,-4)=4*xx*xx2-12*xx*yy2

! Theta function

            th(0,0)=1
            dth(1,0,0)=0
            dth(2,0,0)=0
            dth(3,0,0)=0

            th(1,0)=ri*zz
            dth(1,1,0)=-ri3*xx*zz
            dth(2,1,0)=-ri3*yy*zz
            dth(3,1,0)= ri3*(xx2+yy2)

            th(1,1)=ri
            dth(1,1,1)=-ri3*xx
            dth(2,1,1)=-ri3*yy
            dth(3,1,1)=-ri3*zz

!           th(2,0)=cd0*ri3*(3*zz**2-r2)
            th(2,0)=cd0*(3*zhat**2-1)
            dth(1,2,0)=-6*cd0*ri*zhat**2*xhat
            dth(2,2,0)=-6*cd0*ri*zhat**2*yhat
            dth(3,2,0)= 6*cd0*ri4*(xx2+yy2)*zz

            th(2,1)=cd1*ri2*zz
            dth(1,2,1)=-2*cd1*ri4*zz*xx
            dth(2,2,1)=-2*cd1*ri4*zz*yy
            dth(3,2,1)=   cd1*ri4*(xx2+yy2-zz2)

            th(2,2)=cd2*ri2
            dth(1,2,2)=-2*cd2*ri4*xx
            dth(2,2,2)=-2*cd2*ri4*yy
            dth(3,2,2)=-2*cd2*ri4*zz

            th(3,0)=cf0*ri3*zz*(2*zz2-3*(xx2+yy2))
            dth(1,3,0)= 3*cf0*ri5*(xx2+yy2-4*zz2)*zz*xx
            dth(2,3,0)= 3*cf0*ri5*(xx2+yy2-4*zz2)*zz*yy
            dth(3,3,0)=-3*cf0*ri5*((xx2+yy2)**2-4*zz2*(xx2+yy2))

            th(3,1)=cf1*ri3*(4*zz2-(xx2+yy2))
            dth(1,3,1)=cf1*ri5*(xx2+yy2-14*zz2)*xx
            dth(2,3,1)=cf1*ri5*(xx2+yy2-14*zz2)*yy
            dth(3,3,1)=cf1*ri5*(11*(xx2+yy2)-4*zz2)*zz

            th(3,2)=cf2*ri3*zz
            dth(1,3,2)=-3*cf2*ri5*zz*xx
            dth(2,3,2)=-3*cf2*ri5*zz*yy
            dth(3,3,2)=   cf2*ri5*(xx2+yy2-2*zz2)

            th(3,3)=cf3*ri3
            dth(1,3,3)=-3*cf3*ri5*xx
            dth(2,3,3)=-3*cf3*ri5*yy
            dth(3,3,3)=-3*cf3*ri5*zz

! JWL added l=4
            th(4,0)=cg0*(35*zz2*zz2-30*r2*zz2+3*r2*r2)*ri4
            dth(1,4,0)=cg0*(-60*xx*zz2+12*r2*xx)*ri4 + th(4,0)*(-4*ri2*xx)
            dth(2,4,0)=cg0*(-60*yy*zz2+12*r2*yy)*ri4 + th(4,0)*(-4*ri2*yy)
            dth(3,4,0)=cg0*(140*zz2*zz-(60*zz*zz2+60*r2*zz)+12*r2*zz)*ri4 + th(4,0)*(-4*ri2*zz)

            th(4,1)=cg1*zz*(7*zz2-3*r2)*ri4
            dth(1,4,1)=cg1*zz*(-6*xx)*ri4 + th(4,1)*(-4*ri2*xx)
            dth(2,4,1)=cg1*zz*(-6*yy)*ri4 + th(4,1)*(-4*ri2*yy)
            dth(3,4,1)=cg1*( (7*zz2-3*r2) + zz*(14*zz-6*zz) )*ri4 + th(4,1)*(-4*ri2*zz)

            th(4,2)=cg2*(7*zz2-r2)*ri4
            dth(1,4,2)=cg2*(-2*xx)*ri4 + th(4,2)*(-4*ri2*xx)
            dth(2,4,2)=cg2*(-2*yy)*ri4 + th(4,2)*(-4*ri2*yy)
            dth(3,4,2)=cg2*(14*zz-2*zz)*ri4 + th(4,2)*(-4*ri2*zz)

            th(4,3)=cg3*zz*ri4
            dth(1,4,3)=th(4,3)*(-4*ri2*xx)
            dth(2,4,3)=th(4,3)*(-4*ri2*yy)
            dth(3,4,3)=cg3*ri4 + th(4,3)*(-4*ri2*zz)

            th(4,4)=cg4*ri4
            dth(1,4,4)=th(4,4)*(-4*ri2*xx)
            dth(2,4,4)=th(4,4)*(-4*ri2*yy)
            dth(3,4,4)=th(4,4)*(-4*ri2*zz)

            ib=ib+1
            irb=iwrwf(ib2,ict)
            l=l_bas(ib)
            m=m_bas(ib)
            mabs=abs(m)
            ylm=th(l,mabs)*ph(m)
            phin(ib,ie)=ylm*wfv(1,irb)
            dphin(1,ib,ie)=ylm*xhat*wfv(2,irb)+(th(l,mabs)*dph(1,m)+dth(1,l,mabs)*ph(m))*wfv(1,irb)
            dphin(2,ib,ie)=ylm*yhat*wfv(2,irb)+(th(l,mabs)*dph(2,m)+dth(2,l,mabs)*ph(m))*wfv(1,irb)
            dphin(3,ib,ie)=ylm*zhat*wfv(2,irb)+(                    dth(3,l,mabs)*ph(m))*wfv(1,irb)
            d2phin(ib,ie)=ylm*(wfv(3,irb)+2*ri*wfv(2,irb)-l*(l+1)*ri2*wfv(1,irb))

            if(ipr.ge.1) then
              write(*,'(''ib2,l,m,mabs,ylm:'',4i6,4f16.8)') ib2,l,m,mabs,ylm,(dphin(k,ib,ie),k=1,ndim)
            endif

!     write(6,'(''irb,l,m,ylm,xhat,wfv(1,irb),wfv(2,irb),th(l,mabs),dph(1,m),dth(1,l,mabs),ph(m))'',3i4,9d16.8)')
!    &irb,l,m,ylm,xhat,wfv(1,irb),wfv(2,irb),th(l,mabs),dph(1,m),dth(1,l,mabs),ph(m)
      if(ipr.ge.4) write(6,'(''ie,ib,irb,l,m,ph(m),(dph(k,m),k=1,ndim),th(l,mabs),(dth(k,l,mabs),k=1,ndim) &
     &,wfv(1,irb),wfv(2,irb),wfv(3,irb)'',5i3,20f8.1)') &
     & ie,ib,irb,l,m,ph(m),(dph(k,m),k=1,ndim),th(l,mabs),(dth(k,l,mabs),k=1,ndim),wfv(1,irb),wfv(2,irb),wfv(3,irb)
!     write(6,'(''ie,ib,ib2,irb,phin(ib,ie),dphin,d2phin'',4i5,9d12.5)') ie,ib,ib2,irb,phin(ib,ie),(dphin(kk,ib,ie),kk=1,ndim),d2phin(ib,ie),wfv(1,irb)
!    &,wfv(1,irb)
!        write(6,'(''ib1,ylm'',i3,9f16.10)') ib,ylm
!        write(6,'(''ib1,etc'',i3,9f16.10)')
!    &   ib,ylm*wfv(1,irb),
!    &   ylm*xhat*wfv(2,irb)+(th(l,mabs)*dph(1,m)+dth(1,l,mabs)*ph(m))*wfv(1,irb),
!    &   ylm*yhat*wfv(2,irb)+(th(l,mabs)*dph(2,m)+dth(2,l,mabs)*ph(m))*wfv(1,irb),
!    &   ylm*zhat*wfv(2,irb)+(                    dth(3,l,mabs)*ph(m))*wfv(1,irb),
!    &   ylm*(wfv(3,irb)+2*ri*wfv(2,irb)-l*(l+1)*ri2*wfv(1,irb))
!**RM(5)
      enddo
      else
!*******************
!     New version for recursion construction of Ylm.
!*******************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Ylm is constructed in the form of
!       Ylm = (1/r^l)*[\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}]*coef_ylm(l,m)
!
!     - The recursion formula is based on
!          Y(l,l-1) ~ Y(l-l,l-1)*z       ! increment of 'l'
!          Y(l,m+1) ~ l_{+}*Y(l,m)  etc. ! increment of 'm'
!
!       where l_{+} = -i*(jx) + (jy), and
!             (jx)  = y*dz - z*dy     etc.
!
!     - Coefficents and exponents (a_s,l_s,m_s,n_s) is
!       already stored as module variables large_q(N_lms,1:4) as
!
!                large_q(N_LMs,1) = a_s
!                large_q(N_LMs,2) = l_s
!                large_q(N_LMs,3) = m_s
!                large_q(N_LMs,4) = n_s
!
!       at the beginning of the run by the subroutine
!         'read_input.f -->  call setup_spherical_harmonics'
!
!     - The normalization coefficient coef_ylm(l,m) is also
!       stored as module variable at the beginning of the run by
!       the subroutine
!         'read_input.f -->  call setup_coefficients_ylm'
!
!     - For given (x,y,z),
!       P_lm(x,y,z) = [\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}]/r^L,
!       and dP_lm/dx, dP_lm/dy, and dP_lm/dz are calculated
!       and stored in 'large_p(l,m)' and 'd_large_p(1:3,l,m)'
!       array by calling subroutine 'calculate_spherical_harmonics'.
!
!     - Before calling 'calculate_spherical_harmonics',
!       powers of x,y,z, and 1/r is stored as module variables
!       x_power_of(N) etc.
!---------------

!---------
!* construct power of (x,y,z,1/r)
!---------
        r_inv_power_minus(:) = 1.d0 !* clear first at every centre...
        x_power_of(:) = 0.d0        !* clear first at every centre...
        y_power_of(:) = 0.d0        !* clear first at every centre...
        z_power_of(:) = 0.d0        !* clear first at every centre...

        x_power_of(0) = 1.d0 ; y_power_of(0) = 1.d0 ; z_power_of(0) = 1.d0

        do itemp1 = 1, lmax
           r_inv_power_minus(itemp1) = r_inv_power_minus(itemp1-1)*ri
           x_power_of(itemp1) = x_power_of(itemp1-1)*xc(1)
           y_power_of(itemp1) = y_power_of(itemp1-1)*xc(2)
           z_power_of(itemp1) = z_power_of(itemp1-1)*xc(3)
        enddo

!----------
!*P_lm(x,y,z), dP_lm/dx, dP_lm/dy, and dP_lm/dz is calculated and stored.
!----------
        call calculate_spherical_harmonics

!----------
!*Orbital functions are constructed.
!----------
        do ib2 = 1,nbasis_ctype(ict)
          irb=iwrwf(ib2,ict)
          ib = ib + 1
          l = l_bas(ib)
          m = m_bas(ib)
          mabs=abs(m)
          if(l.gt.lmax) stop 'l>lmax in basis_fns.f'

          ylm           = large_p(l,m)*r_inv_power_minus(l)
          phin(ib,ie)   = ylm*wfv(1,irb)*coef_ylm(l,mabs)
          d2phin(ib,ie)  = coef_ylm(l,mabs)*ylm*( &
     &                      wfv(3,irb) + 2.d0*wfv(2,irb)*r_inv_power_minus(1) &
     &                      - dble(l*(l+1))*wfv(1,irb)*r_inv_power_minus(2))
          aux1           = (wfv(2,irb) - dble(l)*wfv(1,irb)*r_inv_power_minus(1)) &
     &                     *large_p(l,m)*r_inv_power_minus(1)
          dphin(1,ib,ie) = coef_ylm(l,mabs)*(aux1*xc(1) + d_large_p(1,l,m)*wfv(1,irb)) &
     &                     *r_inv_power_minus(l)
          dphin(2,ib,ie) = coef_ylm(l,mabs)*(aux1*xc(2) + d_large_p(2,l,m)*wfv(1,irb)) &
     &                     *r_inv_power_minus(l)
          dphin(3,ib,ie) = coef_ylm(l,mabs)*(aux1*xc(3) + d_large_p(3,l,m)*wfv(1,irb)) &
     &                     *r_inv_power_minus(l)
!       write(6,'(''ib2,ylm'',i3,9f16.10)') ib,ylm
!       write(6,'(''ib2,etc'',i3,9f16.10)')
!    &  ib,ylm*wfv(1,irb)*coef_ylm(l,mabs),
!    &  coef_ylm(l,mabs)*(aux1*xc(1) + d_large_p(1,l,m)*wfv(1,irb))*r_inv_power_minus(l),
!    &  coef_ylm(l,mabs)*(aux1*xc(2) + d_large_p(2,l,m)*wfv(1,irb)) *r_inv_power_minus(l),
!    &  coef_ylm(l,mabs)*(aux1*xc(3) + d_large_p(3,l,m)*wfv(1,irb)) *r_inv_power_minus(l),
!    &  coef_ylm(l,mabs)*ylm*( wfv(3,irb) + 2.d0*wfv(2,irb)*r_inv_power_minus(1) - dble(l*(l+1))*wfv(1,irb)*r_inv_power_minus(2))

!JWL,CJU
        if(ipr.ge.1) then
          write(*,'(''ib2,l,m,mabs,ylm:'',4i6,4f16.8)') ib2,l,m,mabs,coef_ylm(l,mabs)*ylm,(dphin(k,ib,ie),k=1,ndim)
        endif

        enddo
      endif
!**EndRM(5)

   40 continue

!     do 2000 l=1,nbasis
!2000   write(6,'(''basisfns='',20d12.5)') phin(l,nelec),(dphin(ie,l,nelec),ie=1,ndim),d2phin(l,nelec)

      call object_modified_by_index (phin_index) !JT
      call object_modified_by_index (dphin_index) !JT
      call object_modified_by_index (d2phin_index) !JT

      return
      end
!-----------------------------------------------------------------------

      subroutine basis_fns_2d(iel,rvec_en,r_en)
! Written by Cyrus Umrigar
! Calculate 2-dim localised basis functions and their derivatives for all electrons

! In input:
! n1s,n2s,...     > 0 : Slater basis
! nsa,npa,nda     = 0 : asymptotic basis
! n1s,nsa,...     < 0 : Gaussian basis
! Here:
! n_bas2(irb,ict) > 0 : Slater basis
!                 = 0 : asymptotic basis
!                 < 0 : Gaussian basis

      use atom_mod
      use basis1_mod
      use basis_mod, only : which_analytical_basis !fp
      use const_mod
      use numbas_mod
      use basis2_mod
      use wfsec_mod
      use phifun_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

      dimension wfv(3,MRWF),xc(3),ph(-ML_BAS:ML_BAS),dph(3,-ML_BAS:ML_BAS)

      ider=1

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

      do 40 ie=nelec1,nelec2
        ib=0
        do 40 ic=1,ncent
          ict=iwctype(ic)

          if(r_en(ie,ic).eq.0.d0) then
            r_en(ie,ic)=1.d-99
            rvec_en(1,ie,ic)=1.d-99
          endif
          xc(1)=rvec_en(1,ie,ic)
          xc(2)=rvec_en(2,ie,ic)
          r=r_en(ie,ic)
          r2=r*r
          ri=1/r
          ri2=ri*ri
          ri3=ri2*ri

! Loops over analytical basis functions, then numerical.  One can optimize the exponents of the analytic ones.
          do 10 irb=1,nrbas_analytical(ict)
              n=n_bas2(irb,ict)
              rn=abs(n)
              rm3=r**(rn-3)
              rm2=rm3*r
              rm1=rm2*r
              if(n.ne.0) then
                 select case (trim(which_analytical_basis)) !fp
                 case ('slater') !fp
                    zr=zex2(irb,ict,iwf)*r
                    ex=dexp(-zr)
                    wfv(3,irb)=rm3*((rn-1)*(rn-2-zr)-(rn-1-zr)*zr)*ex
                 case ('gaussian') !fp
                    zr=2*zex2(irb,ict,iwf)*r2
                    ex=dexp(-0.5d0*zr)
                    wfv(3,irb)=rm3*((rn-1)*(rn-2)-(2*rn-1-zr)*zr)*ex
                 case ('gauss-slater') !fp
                    zr=(zex2(irb,ict,iwf)*r)**2/(1+zex2(irb,ict,iwf)*r)**2 * (2+zex2(irb,ict,iwf)*r)  !fp
!                    ex=dexp(-zr * (1+zex2(irb,ict,iwf)*r) / (2+zex2(irb,ict,iwf)*r) ) !fp
                    ex=dexp(-(zex2(irb,ict,iwf)*r)**2 / (1+zex2(irb,ict,iwf)*r) ) !fp
                    wfv(3,irb)=rm3*ex *((rn-1-zr)*(rn-2-zr)-zr*(1+2/(1+zex2(irb,ict,iwf)*r)-2/(2+zex2(irb,ict,iwf)*r))) !fp
                 case default
                    write(6,*) 'basis_fns_2d: Allowed basis types are slater gaussian gauss-slater!'
                    stop 'basis_fns_2d: Allowed basis types are slater gaussian gauss-slater!'
                 end select     !fp
              elseif(n.eq.0) then
!     Warning: Asymptotic and Gaussian not yet tested.
!     Asymptotic r^(rn-1)*Exp(-zeta*r), where rn=beta+1, beta=betaq/zeta-1, zeta=sqrt(-2*E_ion)?
                 write(6,*) 'basis_fns: ict=',ict
                 write(6,*) 'basis_fns: irb=',irb
                 write(6,*) 'basis_fns: n=',n
                 stop 'basis_fns: asymptotic not yet fully tested'
                 rn=betaq/zex2(irb,ict,iwf)
                 rm3=r**(rn-3)
                 rm2=rm3*r
                 rm1=rm2*r
                 wfv(3,irb)=rm3*((rn-1)*(rn-2-zr)-(rn-1-zr)*zr)*ex
              endif
              wfv(1,irb)=rm1*ex
   10         wfv(2,irb)=rm2*((rn-1)-zr)*ex

!       write(6,'(''ic,ict,nrbas(ict),wfv='',3i5,29d12.5)') ic,ict,nrbas(ict),(wfv(1,irb),irb=1,nrbas(ict))

            do 20 irb=1,nrbas_numerical(ict)
              rk=r
              call splfit_bas(rk,irb,ict,iwf,wfv(1,nrbas_analytical(ict)+irb),ider)
              if(wfv(1,nrbas_analytical(ict)+irb).eq.0.d0) wfv(1,nrbas_analytical(ict)+irb)=DBLMIN
   20       continue

          xx=xc(1)
          yy=xc(2)
          xhat=xx*ri
          yhat=yy*ri

! Phi function

          ph(0)=1
          dph(1,0)=0
          dph(2,0)=0
          ph(1)=xhat
          ph(-1)=yhat
          dph(1,1)=yhat**2*ri
          dph(2,1)=-xhat*yhat*ri
          dph(1,-1)=-xhat*yhat*ri
          dph(2,-1)=xhat**2*ri
          do 30 i=2,ML_BAS
            ph(i)=ph(i-1)*ph(1)-ph(1-i)*ph(-1)
            ph(-i)=ph(1-i)*ph(1)+ph(i-1)*ph(-1)
            dph(1,i) = i*ph(-i)*ph(-1)*ri
            dph(2,i) =-i*ph(-i)*ph(1)*ri
            dph(1,-i)=-i*ph(i)*ph(-1)*ri
!           write(6,'(''xx,yy,ri,xhat,yhat,ph(i),ph(-i)='',19f9.5)')
!    &      xx,yy,ri,xhat,yhat,ph(i),ph(-i),dph(1,i),dph(2,i),dph(1,-i)
   30       dph(2,-i)= i*ph(i)*ph(1)*ri

          do 40 ib2=1,nbasis_ctype(ict)
            ib=ib+1
            irb=iwrwf(ib2,ict)
            m=m_bas(ib)
            ylm=ph(m)
            phin(ib,ie)=ylm*wfv(1,irb)
            dphin(1,ib,ie)=ylm*xhat*wfv(2,irb)+dph(1,m)*wfv(1,irb)
            dphin(2,ib,ie)=ylm*yhat*wfv(2,irb)+dph(2,m)*wfv(1,irb)
            d2phin(ib,ie)=ylm*(wfv(3,irb)+ri*wfv(2,irb)-m*m*ri2*wfv(1,irb))
!     write(6,'(''m_bas'',9i5)') m
!     write(6,'(''ie,ib,ph(m),(dph(k,m),ddph(m)'',2i5,9d12.5)')
!    & ie,ib,ph(m),(dph(k,m),k=1,ndim),-m*m*ri2*wfv(1,irb),wfv(1,irb),wfv(2,irb),wfv(3,irb),ylm*(wfv(3,irb)+ri*wfv(2,irb)),ylm*(-m*m*ri2*wfv(1,irb))
!     write(6,'(''ie,ib,irb,m,ph(m),(dph(k,m),k=1,ndim),wfv(1,irb),wfv(2,irb),wfv(3,irb)'',4i5,9d12.5)')
!    & ie,ib,irb,m,ph(m),(dph(k,m),k=1,ndim),wfv(1,irb),wfv(2,irb),wfv(3,irb)
!     write(6,'(''ie,ib,ib2,irb,phin(ib,ie),dphin,d2phin'',4i5,9d12.5)')
!    & ie,ib,ib2,irb,phin(ib,ie),(dphin(kk,ib,ie),kk=1,ndim),d2phin(ib,ie)
   40 continue

!     do 2000 l=1,nbasis
!2000   write(6,'(''basisfns='',20d12.5)') phin(l,nelec),(dphin(ie,l,nelec),ie=1,ndim),d2phin(l,nelec)

      return
      end
!-------------------------------------------------------------------------

      subroutine basis_fns_2dgauss(iel,rvec_en,r_en)

! Written by A.D.Guclu, Apr 2006
! 2-dimensional localized gaussian basis set.
! Main purpose is the study 2d wigner crystal.

! arguments: iel=0 -> all electron
!               >0 -> only electron iel
!            rvec_en=vector electron-nucleus
!                    (or electron-dot center in this context)

! output: phin,dphin, and d2phin are calculated

! Wave functions are given by (except the normalization sqrt(we/pi)):
! phi=dsqrt(xg3)*exp(-we*xg3/2*((x1-xg1)^2+(x2-xg2)^2))
! x1,x2 are the electronic positions, xg1 xg2 are the gaussian positions.
! normalization is taken care in (..)

      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      implicit real*8(a-h,o-z)


!     common /dim/ ndim
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
      else
        nelec1=iel
        nelec2=iel
      endif

      ic=1
!      write(6,*) 'in basis_fns'

      do ie=nelec1,nelec2
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)

!       write(6,*) 'x1,x2,we=',x1,x2,we

        do ib=1,nbasis
          wez=we*oparm(3,ib,iwf)
          wez2=wez*wez
          x1rel=x1-oparm(1,ib,iwf)
          x1rel2=x1rel*x1rel
          x2rel=x2-oparm(2,ib,iwf)
          x2rel2=x2rel*x2rel
          rrel2=x1rel2+x2rel2

          phin(ib,ie)=dsqrt(wez)*dexp(-0.5d0*wez*rrel2)

!         write(6,*) 'ib,ie,phin(ib,ie)=',ib,ie,phin(ib,ie)
!         write(6,*) 'oparm1,oparm2,oparm3=',oparm(1,ib,iwf),oparm(2,ib,iwf),oparm(3,ib,iwf)

          dphin(1,ib,ie)=-wez*x1rel*phin(ib,ie)
          dphin(2,ib,ie)=-wez*x2rel*phin(ib,ie)

          d2phin(ib,ie)=(wez2*rrel2-2*wez)*phin(ib,ie)

        enddo
      enddo

      return
      end
!--------------------------------------------------------------------------

      subroutine deriv_2dgauss(rvec_en,r_en)

! Written by A.D.Guclu, Apr 2006
! 2-dimensional localized gaussian basis set,
! and his derivatives wrt parameters.
! Main purpose is the study of 2d wigner crystal.

! arguments:
!            rvec_en=vector electron-nucleus
!                    (or electron-dot center in this context)

! output: phin,dphin,d2phin  = wfs and coo. derivatives
!         dparam, d2param  = parameter derivatives

! Wave functions are given by (except the normalization sqrt(1/pi)):
! phi=dsqrt(we*xg3)*exp(-we*xg3/2*((x1-xg1)^2+(x2-xg2)^2))
! x1,x2 are the electronic positions, xg1 xg2 are the gaussian positions.
! normalization is taken care in (..)

! parameter xg1,xg2,xg3 correspond to nparmo1,nparmo2,nparmo3 respectively.

      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      use deriv_phifun_mod
      implicit real*8(a-h,o-z)


!     common /dim/ ndim
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

      nelec1=1
      nelec2=nelec

      ic=1
!      write(6,*) 'in deriv_2dgauss'
      do ie=nelec1,nelec2
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)

! in the following we are losing some efficiency by calculating all the
! gaussians for each electron. Because up and down electrons do not share
! the same gaussian in crystals we could restrict the calculations...
! I will however keep it this way in case we are interested in other
! application than crystals.
        do ib=1,nbasis
          if(oparm(3,ib,iwf).lt.0.d0) then
            stop 'oparm(3,ib,iwf).lt.0.d0 in deriv_2dgauss. '
          endif
          wez=we*oparm(3,ib,iwf)
          wez2=wez*wez
          xg3i=1/oparm(3,ib,iwf)
          x1rel=x1-oparm(1,ib,iwf)
          x1rel2=x1rel*x1rel
          x2rel=x2-oparm(2,ib,iwf)
          x2rel2=x2rel*x2rel
          rrel2=x1rel2+x2rel2

! wfs and coo. derivatives:
          phin(ib,ie)=dsqrt(wez)*dexp(-0.5d0*wez*rrel2)

          dphin(1,ib,ie)=-wez*x1rel*phin(ib,ie)
          dphin(2,ib,ie)=-wez*x2rel*phin(ib,ie)

          d2phin(ib,ie)=(wez2*rrel2-2*wez)*phin(ib,ie)

! parameter derivatives:
          dparam(1,ib,ie)=-dphin(1,ib,ie)                                  ! wrt xg1
          dparam(2,ib,ie)=-dphin(2,ib,ie)                                  ! wrt xg2
          dparam(3,ib,ie)=0.5d0*(xg3i-we*rrel2)*phin(ib,ie)                ! wrt xg3

          d2param(1,1,ib,ie)=(wez2*x1rel2-wez)*phin(ib,ie)                 ! wrt xg1,xg1
          d2param(2,2,ib,ie)=(wez2*x2rel2-wez)*phin(ib,ie)                 ! wrt xg2,xg2
          d2param(3,3,ib,ie)=0.5d0*(0.5d0*(xg3i-we*rrel2)**2 &
     &                         -xg3i*xg3i)*phin(ib,ie)                     ! wrt xg3,xg3
          d2param(1,2,ib,ie)=wez2*x1rel*x2rel*phin(ib,ie)                  ! wrt xg1,xg2
          temp=we*((phin(ib,ie)+oparm(3,ib,iwf)*dparam(3,ib,ie)))
          d2param(1,3,ib,ie)=x1rel*temp                                    ! wrt xg1,xg3
          d2param(2,3,ib,ie)=x2rel*temp                                    ! wrt xg2,xg3
          d2param(2,1,ib,ie)=d2param(1,2,ib,ie)
          d2param(3,1,ib,ie)=d2param(1,3,ib,ie)
          d2param(3,2,ib,ie)=d2param(2,3,ib,ie)

! coo. derivatives of parameter derivatives:
          ddparam(1,1,ib,ie)=-d2param(1,1,ib,ie)                             ! wrt x1,xg1
          ddparam(2,1,ib,ie)=-d2param(2,1,ib,ie)                             ! wrt x2,xg1

          ddparam(1,2,ib,ie)=ddparam(2,1,ib,ie)                            ! wrt x1,xg2
          ddparam(2,2,ib,ie)=-d2param(2,2,ib,ie)                             ! wrt x2,xg2

          ddparam(1,3,ib,ie)=-d2param(1,3,ib,ie)                             ! wrt x1,xg3
          ddparam(2,3,ib,ie)=-d2param(2,3,ib,ie)                             ! wrt x2,xg3

          d2dparam(1,ib,ie)=wez*(2*dparam(1,ib,ie)+x1rel*d2phin(ib,ie))   ! laplacian of dparam(1,ib,ie)
          d2dparam(2,ib,ie)=wez*(2*dparam(2,ib,ie)+x2rel*d2phin(ib,ie))   ! laplacian of dparam(2,ib,ie)
          d2dparam(3,ib,ie)=xg3i*d2phin(ib,ie)-2*wez*dparam(3,ib,ie) &
     &+wez*we*(phin(ib,ie)+oparm(3,ib,iwf)*dparam(3,ib,ie))*rrel2     ! laplacian of dparam(3,ib,ie)
!         test=2*we*phin(ib,ie)+0.5d0*(5*xg3i-we*rrel2)*d2phin(ib,ie)
!         write(6,*) 'ib,ie,d2dparam,test=',ib,ie,d2dparam(3,ib,ie),test
!         d2dparam(3,ib,ie)=test

        enddo
      enddo

      return
      end

!--------------------------------------------------------------------------

      subroutine basis_fns_polargauss(iel,rvec_en,r_en)

! Written by A.D.Guclu, Jan 2007
! 2-dimensional localized gaussian basis set in polar coordinates
! Main purpose is the study quasi-1D wigner crystal. It can
! also be applied to 2d wigner crystals.

! arguments: iel=0 -> all electron
!               >0 -> only electron iel
!            rvec_en=vector electron-nucleus
!                    (or electron-dot center in this context)

! output: phin,dphin, and d2phin are calculated

! Wave functions are given by (except the normalization csnt)

! phi=dsqrt(xg3)*exp(-we*xg3/2*((x1-xg1)^2+(x2-xg2)^2))

! phi=dsqrt(xg3) * exp(-we*xg3/2*(xr-xr0)^2) * exp(xg4*(cos(xt-xt0)-1))

! where  xr=r and xt=\theta  are polar coordinates of electrons
!        xg1 and xg2 are polar coordinates of gaussians
!        xg3 and xg4 are gaussian width parameters.

! Here, we do not normalize the wfs



      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      implicit real*8(a-h,o-z)


!     common /dim/ ndim
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

!      write(6,*) 'oparm1s=',oparm(1,1,iwf),oparm(1,2,iwf),oparm(1,3,iwf),oparm(1,4,iwf)

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
      else
        nelec1=iel
        nelec2=iel
      endif

      ic=1
      expnorm=1.d0
!      write(6,*) 'in basis_fns'

      do ie=nelec1,nelec2
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)
        xr=r_en(ie,ic)
        xt=datan2(x2,x1)
        xri=1/xr
        xri2=xri*xri

!        write(6,*) 'x,y,xr,xt,xr*dcos(xt)=',x1,x2,xr,xt,xr*dcos(xt)

        do ib=1,nbasis

          xg1=oparm(1,ib,iwf)
          xg2=oparm(2,ib,iwf)
          xg3=oparm(3,ib,iwf)
          xg4=oparm(4,ib,iwf)
          wez=we*xg3
          wez2=wez*wez
          xrrel=xr-xg1
          xrrel2=xrrel*xrrel
          xtrel=xt-xg2
          cxtrel=dcos(xtrel)
          sxtrel=dsin(xtrel)

          phir=dsqrt(we)*dexp(-0.5d0*wez*xrrel2)
          phit=dexp(xg4*(cxtrel-expnorm))
          phin(ib,ie)=phir*phit

!         if(abs(phir).lt.1.d-300) then
!           write(6,'(''phir='',9d12.4)') phir
!           stop 'phir too small'
!         endif

          if(abs(phir).gt.1.d+300) then
            write(6,'(''phir='',9d12.4)') phir
            stop 'phir too large'
          endif

!         if(abs(phir).lt.1.d-300) then
!           write(6,'(''phit='',9d12.4)') phit
!           stop 'phit too small'
!         endif

          if(abs(phit).gt.1.d+300) then
            write(6,'(''phit='',9d12.4)') phit
            stop 'phit too large'
          endif

!          write(6,*) 'phin(ib,ie),phir,phit,ib,ie=',phin(ib,ie),phir,phit,ib,ie
!          write(6,*) 'xg1,xg2,xg3,xg4,wez,wez2,',
!     &'xrrel,xrrel2,xtrel,cxtrel,sxtrel=',
!     &xg1,xg2,xg3,xg4,wez,wez2,xrrel,xrrel2,xtrel,cxtrel,sxtrel
!          write(6,*) 'xg1,xg2,xg3,xg4=',xg1,xg2,xg3,xg4

          dpdxr=-wez*xrrel
          dpdxt=-xg4*sxtrel
          d2pdxr2=-wez+dpdxr*dpdxr
          d2pdxt2=-xg4*cxtrel+dpdxt*dpdxt

          dphin(1,ib,ie)=(dpdxr*x1*xri-dpdxt*x2*xri2)*phin(ib,ie)
          dphin(2,ib,ie)=(dpdxr*x2*xri+dpdxt*x1*xri2)*phin(ib,ie)

          d2phin(ib,ie)=(dpdxr*xri+d2pdxr2+d2pdxt2*xri2)*phin(ib,ie)

!         if(abs(phin(ib,ie)).lt.1.d-300) then
!           write(6,'(''phir,phit,phin(ib,ie),dphin(1,ib,ie),dphin(2,ib,ie),d2phin(ib,ie)='',9d12.4)')
!    &      phir,phit,phin(ib,ie),dphin(1,ib,ie),dphin(2,ib,ie),d2phin(ib,ie)
!           stop 'phir,phit,phin(ib,ie) too small'
!         endif

          if(abs(phin(ib,ie)).gt.1.d+300) then
            write(6,'(''phir,phit,phin(ib,ie),dphin(1,ib,ie),dphin(2,ib,ie),d2phin(ib,ie)='',9d12.4)') &
     &      phir,phit,phin(ib,ie),dphin(1,ib,ie),dphin(2,ib,ie),d2phin(ib,ie)
            stop 'phir,phit,phin(ib,ie) too large'
          endif

        enddo
      enddo

! Check that every electron has at least one nonzero basis function.
      phimax=maxval(abs(phin))
!     write(6,'(''phimax='',d12.4)') phimax

      do ie=nelec1,nelec2
        phicolmax=maxval(abs(phin(:,ie)))
        if(phicolmax.lt.1.d-100*phimax) write(6,'(''Warning: ie, phicolmax, phimax='',i5,9e12.4)')  ie, phicolmax, phimax
        if(phicolmax.eq.0.d0) then
          write(6,'(''Warning stop: ie, phicolmax, phimax='',i5,9e12.4)')  ie, phicolmax, phimax
          stop 'phicolmax = 0'
        endif
      enddo

!     do ib=1,nbasis
!       phirowmax=maxval(abs(phin(ib,:)))
!       if(phirowmax.lt.1.d-100*phimax) write(6,'(''Warning: ib, phirowmax, phimax='',i5,9e12.4)')  ib, phirowmax, phimax
!       if(phirowmax.eq.0.d0) then
!         write(6,'(''Warning stop: ib, phirowmax, phimax='',i5,9e12.4)')  ib, phicolmax, phimax
!         stop 'phirowmax = 0'
!       endif
!     enddo

      return
      end

!-------------------------------------------------------------------------
      subroutine deriv_polargauss(rvec_en,r_en)

! Written by A.D.Guclu, Jan 2007
! 2-dimensional localized gaussian basis set in polar coordinates
! Main purpose is the study quasi-1D wigner crystal. It can
! also be applied to 2d wigner crystals.

! arguments: iel=0 -> all electron
!               >0 -> only electron iel
!            rvec_en=vector electron-nucleus
!                    (or electron-dot center in this context)

! output: phin,dphin, and d2phin are calculated
!         dparam, d2param  = parameter derivatives

! Wave functions are given by (except the normalization csnt)

! phi=dsqrt(xg3)*exp(-we*xg3/2*((x1-xg1)^2+(x2-xg2)^2))

! phi=dsqrt(xg3) * exp(-we*xg3/2*(xr-xr0)^2) * exp(xg4*cos(xt-xt0))

! where  xr=r and xt=\theta  are polar coordinates of electrons
!        xg1 and xg2 are polar coordinates of gaussians
!        xg3 and xg4 are gaussian width parameters.

! Here, we do not normalize the wfs

      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      use deriv_phifun_mod
      implicit real*8(a-h,o-z)


!     common /dim/ ndim
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

      nelec1=1
      nelec2=nelec

      ic=1
      expnorm=1.d0

      do ie=nelec1,nelec2
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)
        xr=r_en(ie,ic)
        xt=datan2(x2,x1)
        xri=1/xr
        xri2=xri*xri

!       write(6,*) 'x,y,xr,xt,xr*dcos(xt)=',x1,x2,xr,xt,xr*dcos(xt)

        do ib=1,nbasis

          xg1=oparm(1,ib,iwf)
          xg2=oparm(2,ib,iwf)
          xg3=oparm(3,ib,iwf)
          xg4=oparm(4,ib,iwf)
          wez=we*xg3
          wez2=wez*wez
          xrrel=xr-xg1
          xrrel2=xrrel*xrrel
          xtrel=xt-xg2
          cxtrel=dcos(xtrel)
          sxtrel=dsin(xtrel)
! jacobian:
          dxrdx1=x1*xri
          dxrdx2=x2*xri
          dxtdx1=-x2*xri2
          dxtdx2=x1*xri2

! wfs and coo. derivatives:
          phir=dsqrt(we)*dexp(-0.5d0*wez*xrrel2)
          phit=dexp(xg4*(cxtrel-expnorm))
          phin(ib,ie)=phir*phit

!         write(6,*) 'phir,phit=',phir,phit

          tempr1=-wez*xrrel
          tempt1=-xg4*sxtrel
          tempr2=-wez+tempr1*tempr1
          tempt2=-xg4*cxtrel+tempt1*tempt1
          dpdxr=tempr1*phin(ib,ie)
          dpdxt=tempt1*phin(ib,ie)
          d2pdxr2=tempr2*phin(ib,ie)
          d2pdxt2=tempt2*phin(ib,ie)

          dphin(1,ib,ie)=dpdxr*dxrdx1+dpdxt*dxtdx1
          dphin(2,ib,ie)=dpdxr*dxrdx2+dpdxt*dxtdx2

          d2phin(ib,ie)=dpdxr*xri+d2pdxr2+d2pdxt2*xri2

! parameter derivatives:
          dparam(1,ib,ie)=-dpdxr                                 ! wrt xg1
          dparam(2,ib,ie)=-dpdxt                                 ! wrt xg2
          dparam(3,ib,ie)=-0.5d0*we*xrrel2*phin(ib,ie)           ! wrt xg3
          dparam(4,ib,ie)=(cxtrel-expnorm)*phin(ib,ie)           ! wrt xg4

          d2param(1,1,ib,ie)=d2pdxr2                             ! wrt xg1,xg1
          d2param(2,2,ib,ie)=d2pdxt2                             ! wrt xg2,xg2
          d2param(3,3,ib,ie)=-0.5d0*we*xrrel2*dparam(3,ib,ie)    ! wrt xg3,xg3
          d2param(4,4,ib,ie)=(cxtrel-expnorm)*dparam(4,ib,ie)    ! wrt xg4,xg4

          d2param(1,2,ib,ie)=-tempr1*dparam(2,ib,ie)
          d2param(1,3,ib,ie)=we*xrrel*(1.d0-0.5d0*wez*xrrel2)*phin(ib,ie)
          d2param(1,4,ib,ie)=-tempr1*dparam(4,ib,ie)
          d2param(2,3,ib,ie)=-tempt1*dparam(3,ib,ie)
          d2param(2,4,ib,ie)=sxtrel*(1.d0+xg4*(cxtrel-expnorm))*phin(ib,ie)
          d2param(3,4,ib,ie)=-0.5d0*we*xrrel2*dparam(4,ib,ie)

          d2param(2,1,ib,ie)=d2param(1,2,ib,ie)
          d2param(3,1,ib,ie)=d2param(1,3,ib,ie)
          d2param(4,1,ib,ie)=d2param(1,4,ib,ie)
          d2param(3,2,ib,ie)=d2param(2,3,ib,ie)
          d2param(4,2,ib,ie)=d2param(2,4,ib,ie)
          d2param(4,3,ib,ie)=d2param(3,4,ib,ie)

! coo. derivatives of parameter derivatives:

          ddparam(1,1,ib,ie)=                                 & ! wrt x1,xg1
     &-(d2param(1,1,ib,ie)*dxrdx1+d2param(1,2,ib,ie)*dxtdx1)
          ddparam(2,1,ib,ie)=                                 & ! wrt x2,xg1
     &-(d2param(1,1,ib,ie)*dxrdx2+d2param(1,2,ib,ie)*dxtdx2)

          ddparam(1,2,ib,ie)= &
     &-(d2param(2,1,ib,ie)*dxrdx1+d2param(2,2,ib,ie)*dxtdx1)
          ddparam(2,2,ib,ie)= &
     &-(d2param(2,1,ib,ie)*dxrdx2+d2param(2,2,ib,ie)*dxtdx2)

          ddparam(1,3,ib,ie)= &
     &-(d2param(3,1,ib,ie)*dxrdx1+d2param(3,2,ib,ie)*dxtdx1)
          ddparam(2,3,ib,ie)= &
     &-(d2param(3,1,ib,ie)*dxrdx2+d2param(3,2,ib,ie)*dxtdx2)

          ddparam(1,4,ib,ie)= &
     &-(d2param(4,1,ib,ie)*dxrdx1+d2param(4,2,ib,ie)*dxtdx1)
          ddparam(2,4,ib,ie)= &
     &-(d2param(4,1,ib,ie)*dxrdx2+d2param(4,2,ib,ie)*dxtdx2)


          d2dparam(1,ib,ie)=-d2param(1,1,ib,ie)*xri &
     &-2*wez2*xrrel*phin(ib,ie)+wez*xrrel*d2pdxr2 &
     &+xri2*tempt2*dparam(1,ib,ie)
          d2dparam(2,ib,ie)=-d2param(1,2,ib,ie)*xri &
     &+tempr2*dparam(2,ib,ie) &
     &+xri2*(-xg4*sxtrel*(1.d0+2*xg4*cxtrel)*phin(ib,ie) &
     &       +tempt2*dparam(2,ib,ie))
          d2dparam(3,ib,ie)=-d2param(1,3,ib,ie)*xri &
     &+we*(-1.d0+2*wez*xrrel2)*phin(ib,ie)-0.5d0*we*xrrel2*d2pdxr2 &
     &+xri2*tempt2*dparam(3,ib,ie)
          d2dparam(4,ib,ie)=-d2param(1,4,ib,ie)*xri &
     &+tempr2*dparam(4,ib,ie) &
     &+xri2*((-cxtrel+2*xg4*sxtrel*sxtrel)*phin(ib,ie) &
     &+tempt2*dparam(4,ib,ie))

        enddo
      enddo

      return
      end


!-------------------------------------------------------------------------

      subroutine basis_fns_2dgauss_noncirc(iel,rvec_en,r_en)

! Written by Abhijit C. Mehta, October 2008
! (modified version of basis_fns_2dgauss)
! 2-dimensional localized gaussian basis set,
!    gaussians can have diffferent widths in x- and y- directions
! Main purpose is the study 2d wigner crystal in quantum wires

! arguments: iel=0 -> all electron
!               >0 -> only electron iel
!            rvec_en=vector electron-nucleus
!                    (or electron-dot/"wire" center in this context)

! output: phin,dphin, and d2phin are calculated

! Wave functions are given by (except the normalization sqrt(we/pi)):

! phi=dsqrt(dsqrt(xg3*xg4)) * exp(-we*xg3/2*((x1-xg1)^2))
!                           * exp(-we*xg4/2*((x2-xg2)^2))

! where x1,x2 are the electronic x and y positions,
!       xg1, xg2 are the gaussian x and y positions.
!       xg3, xg4 are the gaussian x and y width parameters
!       (Note that xg3 and xg4 are the "omegas" (1/width)
!         in units of we (or wire_w in case of wire)))


      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      implicit real*8(a-h,o-z)


!     common /dim/ ndim
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
      else
        nelec1=iel
        nelec2=iel
      endif

      ic=1
!      write(6,*) 'in basis_fns'

      do ie=nelec1,nelec2
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)

!       write(6,*) 'x1,x2,we=',x1,x2,we

        do ib=1,nbasis
          wex=we*oparm(3,ib,iwf)
          wex2=wex*wex
          wey=we*oparm(4,ib,iwf)
          wey2=wey*wey
          x1rel=x1-oparm(1,ib,iwf)
          x1rel2=x1rel*x1rel
          x2rel=x2-oparm(2,ib,iwf)
          x2rel2=x2rel*x2rel

          phin(ib,ie)=dsqrt(dsqrt(wex*wey))*dexp(-0.5d0*wex*x1rel2)*dexp(-0.5d0*wey*x2rel2)

!         write(6,*) 'ib,ie,phin(ib,ie)=',ib,ie,phin(ib,ie)
!         write(6,*) 'oparm1,oparm2,oparm3,oparm4=',oparm(1,ib,iwf),oparm(2,ib,iwf),oparm(3,ib,iwf),oparm(4,ib,iwf)

          dphin(1,ib,ie)=-wex*x1rel*phin(ib,ie)
          dphin(2,ib,ie)=-wey*x2rel*phin(ib,ie)

          d2phin(ib,ie)=(wex2*x1rel2 + wey2*x2rel2 - wex - wey)*phin(ib,ie)

        enddo
      enddo

      return
      end

!--------------------------------------------------------------------------

      subroutine deriv_2dgauss_noncirc(rvec_en,r_en)

! Written by Abhijit C. Mehta, October 2008
! modified version of deriv_2dgauss
! 2-dimensional localized gaussian basis set, with different x- and y- widths,
! and the derivatives wrt parameters.
! Main purpose is the study of 2d wigner crystal in quantum wires

! arguments:
!            rvec_en=vector electron-nucleus
!                    (or electron-dot center in this context)

! output: phin,dphin,d2phin  = wfs and coo. derivatives
!         dparam, d2param, ddparam, d2dparam  = parameter derivatives

! Wave functions are given by (except the normalization sqrt(we/pi)):

! phi=dsqrt(dsqrt(xg3*xg4)) * exp(-we*xg3/2*((x1-xg1)^2))
!                           * exp(-we*xg4/2*((x2-xg2)^2))

! where x1,x2 are the electronic x and y positions,
!       xg1, xg2 are the gaussian x and y positions.
!       xg3, xg4 are the gaussian x and y width parameters
!       (Note that xg3 and xg4 are the "omegas" (1/width)
!         in units of we (or wire_w in case of wire)))

! parameters xg1,xg2,xg3,xg4 correspond to nparmo1,nparmo2,nparmo3,nparmo4
      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      use deriv_phifun_mod
      implicit real*8(a-h,o-z)


!     common /dim/ ndim
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

      nelec1=1
      nelec2=nelec

      ic=1
!      write(6,*) 'in deriv_2dgauss'
      do ie=nelec1,nelec2
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)

! in the following we are losing some efficiency by calculating all the
! gaussians for each electron. Because up and down electrons do not share
! the same gaussian in crystals we could restrict the calculations...
! I will however keep it this way in case we are interested in other
! application than crystals.
        do ib=1,nbasis
          if(oparm(3,ib,iwf).lt.0.d0) then
            stop 'oparm(3,ib,iwf).lt.0.d0 in deriv_2dgauss_noncirc. '
          elseif(oparm(4,ib,iwf).lt.0.d0) then
            stop 'oparm(4,ib,iwf).lt.0.d0 in deriv_2dgauss_noncirc. '
          endif

          wex=we*oparm(3,ib,iwf)
          wex2=wex*wex
          wey=we*oparm(4,ib,iwf)
          wey2=wey*wey
          xg3i=1/oparm(3,ib,iwf)
          xg4i=1/oparm(4,ib,iwf)
          x1rel=x1-oparm(1,ib,iwf)
          x1rel2=x1rel*x1rel
          x2rel=x2-oparm(2,ib,iwf)
          x2rel2=x2rel*x2rel

! wfs and coo. derivatives:
          phin(ib,ie)=dsqrt(dsqrt(wex*wey))*dexp(-0.5d0*wex*x1rel2)*dexp(-0.5d0*wey*x2rel2)

          dphin(1,ib,ie)=-wex*x1rel*phin(ib,ie)
          dphin(2,ib,ie)=-wey*x2rel*phin(ib,ie)
          tempd2 = (wex2*x1rel2 + wey2*x2rel2 - wex - wey)
          d2phin(ib,ie)=tempd2*phin(ib,ie)


! parameter derivatives:
          dparam(1,ib,ie)=-dphin(1,ib,ie)                                  ! wrt xg1
          dparam(2,ib,ie)=-dphin(2,ib,ie)                                  ! wrt xg2
          dparam(3,ib,ie)=0.5d0*(0.5d0*xg3i - we*x1rel2)* phin(ib,ie)      ! wrt xg3
          dparam(4,ib,ie)=0.5d0*(0.5d0*xg4i - we*x2rel2)* phin(ib,ie)      ! wrt xg4


          d2param(1,1,ib,ie)=(wex2*x1rel2-wex)*phin(ib,ie)                 ! wrt xg1,xg1
          d2param(2,2,ib,ie)=(wey2*x2rel2-wey)*phin(ib,ie)                 ! wrt xg2,xg2
          d2param(3,3,ib,ie)=(0.25d0*we*x1rel2*(we*x1rel2 - xg3i) &
     &                        - 0.1875d0*xg3i*xg3i)*phin(ib,ie)            ! wrt xg3,xg3
          d2param(4,4,ib,ie)=(0.25d0*we*x2rel2*(we*x2rel2 - xg4i) &
     &                        - 0.1875d0*xg4i*xg4i)*phin(ib,ie)            ! wrt xg4,xg4
          d2param(1,2,ib,ie)=wex*wey*x1rel*x2rel*phin(ib,ie)               ! wrt xg1,xg2

          d2param(1,3,ib,ie)=x1rel*(wex*dparam(3,ib,ie) + we*phin(ib,ie))  ! wrt xg1,xg3
          d2param(1,4,ib,ie)=x1rel*wex*dparam(4,ib,ie)                     ! wrt xg1,xg3
          d2param(2,3,ib,ie)=x2rel*wey*dparam(3,ib,ie)                     ! wrt xg2,xg3
          d2param(2,4,ib,ie)=x2rel*(wey*dparam(4,ib,ie) + we*phin(ib,ie))  ! wrt xg2,xg4
          d2param(3,4,ib,ie)=0.5d0*(0.5d0*xg4i - we*x2rel2)*dparam(3,ib,ie) ! wrt xg3,xg4

          d2param(2,1,ib,ie)=d2param(1,2,ib,ie)
          d2param(3,1,ib,ie)=d2param(1,3,ib,ie)
          d2param(4,1,ib,ie)=d2param(1,4,ib,ie)
          d2param(3,2,ib,ie)=d2param(2,3,ib,ie)
          d2param(4,2,ib,ie)=d2param(2,4,ib,ie)
          d2param(4,3,ib,ie)=d2param(3,4,ib,ie)

! coo. derivatives of parameter derivatives:
          ddparam(1,1,ib,ie)=-d2param(1,1,ib,ie)                            ! wrt x1,xg1
          ddparam(2,1,ib,ie)=-d2param(2,1,ib,ie)                            ! wrt x2,xg1

          ddparam(1,2,ib,ie)= ddparam(2,1,ib,ie)                            ! wrt x1,xg2
          ddparam(2,2,ib,ie)=-d2param(2,2,ib,ie)                             ! wrt x2,xg2

          ddparam(1,3,ib,ie)=-d2param(1,3,ib,ie)                             ! wrt x1,xg3
          ddparam(2,3,ib,ie)=-d2param(2,3,ib,ie)                             ! wrt x2,xg3

          ddparam(1,4,ib,ie)=-d2param(1,4,ib,ie)                             ! wrt x1,xg4
          ddparam(2,4,ib,ie)=-d2param(2,4,ib,ie)                             ! wrt x2,xg4

          d2dparam(1,ib,ie)=wex*(2*dparam(1,ib,ie)+x1rel*d2phin(ib,ie))  ! laplacian of dparam(1,ib,ie)
          d2dparam(2,ib,ie)=wey*(2*dparam(2,ib,ie)+x2rel*d2phin(ib,ie))   ! laplacian of dparam(2,ib,ie)
          d2dparam(3,ib,ie)=(2*wex2*xg3i*x1rel2 - we)*phin(ib,ie) &
     &                      + tempd2*dparam(3,ib,ie)  ! laplacian of dparam(3,ib,ie)
          d2dparam(4,ib,ie)=(2*wey2*xg4i*x2rel2 - we)*phin(ib,ie) &
     &                      + tempd2*dparam(4,ib,ie)  ! laplacian of dparam(4,ib,ie)

        enddo
      enddo

      return
      end



!-------------------------------------------------------------------------

      subroutine basis_fns_2dgauss_periodic(iel,rvec_en,r_en)

! Written by Abhijit C. Mehta, December 2009
! (modified version of basis_fns_2dgauss)
! 2-dimensional localized gaussian basis set,
!    gaussians can have diffferent widths in x- and y- directions
!    periodic in x-direction (designed for iperiodic = 1)
! Main purpose is the study 2d wigner crystal in quantum wires
!    (infinite periodic wire - i.e., 1D periodic boundary conditions)

! arguments: iel=0 -> all electron
!               >0 -> only electron iel
!            rvec_en=vector electron-nucleus
!                    (or electron-dot/"wire" center in this context)
!                    (i.e., distance from (0,0))

! output: phin,dphin, and d2phin are calculated

! Wave functions are given by (except the normalization sqrt(we/pi)):

! phi=dsqrt(dsqrt(xg3*xg4)) * exp(-we*xg3/2*((x1-xg1)^2))
!                           * exp(-we*xg4/2*((x2-xg2)^2))
!                + all periodic images
!  i.e. we take phi(x1) + phi(x1 + a) + phi(x1 + 2a) + phi(x1 - a) + ...

! where x1,x2 are the electronic x and y positions,
!       xg1, xg2 are the gaussian x and y positions.
!       xg3, xg4 are the gaussian x and y width parameters
!       (Note that xg3 and xg4 are the "omegas" (1/width)
!         in units of we (or wire_w in case of wire)))


      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      use periodic_1d_mod

      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)
      gausseps = 1.d-12  ! We sum gaussians with value greater than gausseps

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
      else
        nelec1=iel
        nelec2=iel
      endif

      ic=1
!      write(6,*) 'in basis_fns'

      do ie=nelec1,nelec2
!       modulo math was done in distances.f, so we don't need to do it here
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)

!       write(6,*) 'x1,x2,we=',x1,x2,we

!       assume that input parameters oparm are scaled to periodic BC's
        do ib=1,nbasis
          wex=we*oparm(3,ib,iwf)
          wex2=wex*wex
          wey=we*oparm(4,ib,iwf)
          wey2=wey*wey
          x1rel=x1-oparm(1,ib,iwf)
!         modulo math:
          if(x1rel.ge.(alattice/2.0)) then
             x1rel = x1rel - alattice
          elseif(x1rel.lt.(-alattice/2.0)) then
             x1rel = alattice + x1rel
          endif

          if((x1rel.gt.(alattice/2.0).or.(x1rel.lt.(-alattice/2.0)))) then
             write(6,*) 'Error in modulo math in basis_fns_2dgauss_periodic'
             stop 'Error in modulo math in basis_fns_2dgauss_periodic'
          endif

          x1rel2=x1rel*x1rel
          x2rel=x2-oparm(2,ib,iwf)
          x2rel2=x2rel*x2rel

          phinypart=dsqrt(dsqrt(wex*wey))*dexp(-0.5d0*wey*x2rel2)
          phinxpart=dexp(-0.5d0*wex*x1rel2)
!         The following line was the way we originally wanted to do this, but
!            there were strange compiler errors, and the code would often fail to
!            exit the following do loop if gausseps was too small.
!          For now, we just set gausseps = 1.d-12, but another thing to try
!            would be to try gausseps = abs(1.d-12 * phinxpart)
!          gausseps = 1.d-12 * phinxpart
          dphinxpart = x1rel*phinxpart
          d2phinxpart = x1rel2*phinxpart

          do icell = 1,3
             x1relleft = x1rel - alattice*icell
             x1relright = x1rel + alattice*icell
             x1relleft2 = x1relleft*x1relleft
             x1relright2 = x1relright*x1relright
             phileft = dexp(-0.5d0*(wex*x1relleft2))
             phiright = dexp(-0.5d0*(wex*x1relright2))
             if((phileft.le.gausseps).and.(phiright.le.gausseps)) exit
             phinxpart = phinxpart + phileft + phiright
             dphinxpart = dphinxpart + x1relleft*phileft + x1relright*phiright
             d2phinxpart = d2phinxpart + x1relleft2*phileft + x1relright2*phiright
!             if (icell.gt.2) then
!                write (6,*) 'Warning: in basis_fns_2dgauss_periodic: gaussians are wider than 2 unit cells.'
!             endif
          enddo

          phin(ib,ie)=phinypart*phinxpart


!         write(6,*) 'ib,ie,phin(ib,ie)=',ib,ie,phin(ib,ie)
!         write(6,*) 'oparm1,oparm2,oparm3,oparm4=',oparm(1,ib,iwf),oparm(2,ib,iwf),oparm(3,ib,iwf),oparm(4,ib,iwf)

          dphin(1,ib,ie)=-wex*dphinxpart*phinypart
          dphin(2,ib,ie)=-wey*x2rel*phin(ib,ie)

          d2phin(ib,ie)=wex2*d2phinxpart*phinypart + (wey2*x2rel2 - wex - wey)*phin(ib,ie)

        enddo
      enddo

      return
      end

!--------------------------------------------------------------------------

      subroutine deriv_2dgauss_periodic(rvec_en,r_en)

! Written by Abhijit C. Mehta, January 2010
! modified version of deriv_2dgauss
! 2-dimensional localized gaussian basis set, with different x- and y- widths,
! and the derivatives wrt parameters. for a periodic system
! Main purpose is the study of 2d wigner crystal in quantum wires
!  (Infinite wire, periodic BC's)

! arguments:
!            rvec_en=vector electron-nucleus
!                    (or electron-dot center in this context)

! output: phin,dphin,d2phin  = wfs and coo. derivatives
!         dparam, d2param, ddparam, d2dparam  = parameter derivatives

! Wave functions are given by (except the normalization sqrt(we/pi)):

! phi=dsqrt(dsqrt(xg3*xg4)) * exp(-we*xg3/2*((x1-xg1)^2))
!                           * exp(-we*xg4/2*((x2-xg2)^2))
!                + all periodic images
!  i.e. we take phi(x1) + phi(x1 + a) + phi(x1 + 2a) + phi(x1 - a) + ...

! where x1,x2 are the electronic x and y positions,
!       xg1, xg2 are the gaussian x and y positions.
!       xg3, xg4 are the gaussian x and y width parameters
!       (Note that xg3 and xg4 are the "omegas" (1/width)
!         in units of we (or wire_w in case of wire)))

! parameters xg1,xg2,xg3,xg4 correspond to nparmo1,nparmo2,nparmo3,nparmo4
      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      use deriv_phifun_mod
      use periodic_1d_mod
      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)
      gausseps = 1.d-12  ! We sum gaussians with value greater than gausseps

      nelec1=1
      nelec2=nelec

      ic=1
!      write(6,*) 'in deriv_2dgauss'
      do ie=nelec1,nelec2
!       modulo math already done in distances.f
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)

! in the following we are losing some efficiency by calculating all the
! gaussians for each electron. Because up and down electrons do not share
! the same gaussian in crystals we could restrict the calculations...
! I will however keep it this way in case we are interested in other
! application than crystals.
        do ib=1,nbasis
          if(oparm(3,ib,iwf).lt.0.d0) then
            stop 'oparm(3,ib,iwf).lt.0.d0 in deriv_2dgauss_periodic. '
          elseif(oparm(4,ib,iwf).lt.0.d0) then
            stop 'oparm(4,ib,iwf).lt.0.d0 in deriv_2dgauss_periodic. '
          endif

          wex=we*oparm(3,ib,iwf)
          wex2=wex*wex
          wey=we*oparm(4,ib,iwf)
          wey2=wey*wey
          xg3i=1/oparm(3,ib,iwf)
          xg4i=1/oparm(4,ib,iwf)
          x1rel=x1-oparm(1,ib,iwf)
!         modulo math:
          if(x1rel.ge.(alattice/2.0)) then
             x1rel = x1rel - alattice
          elseif(x1rel.lt.(-alattice/2.0)) then
             x1rel = alattice + x1rel
          endif
          if((x1rel.gt.(alattice/2.0).or.(x1rel.lt.(-alattice/2.0)))) then
             write(6,*) 'Error in modulo math in deriv_2dgauss_periodic'
             stop 'Error in modulo math in basis_fns_2dgauss_periodic'
          endif

          x1rel2=x1rel*x1rel
          x2rel=x2-oparm(2,ib,iwf)
          x2rel2=x2rel*x2rel

! wfs and coo. derivatives:

          phinypart=dsqrt(dsqrt(wex*wey))*dexp(-0.5d0*wey*x2rel2)
          phinxpart=dexp(-0.5d0*wex*x1rel2)
!         The following line was the way we originally wanted to do this, but
!            there were strange compiler errors, and the code would often fail to
!            exit the following do loop if gausseps was too small.
!          For now, we just set gausseps = 1.d-12, but another thing to try
!            would be to try gausseps = abs(1.d-12 * phinxpart)
!          gausseps = 1.0d-12*phinxpart
          dphinxpart = x1rel*phinxpart
          d2phinxpart = x1rel2*phinxpart
          d3phinxpart = x1rel*x1rel2*phinxpart
          d4phinxpart = x1rel2*x1rel2*phinxpart
          do icell = 1,3
             x1relleft = x1rel - alattice*icell
             x1relright = x1rel + alattice*icell
             x1relleft2 = x1relleft*x1relleft
             x1relright2 = x1relright*x1relright
             phileft = dexp(-0.5d0*(wex*x1relleft2))
             phiright = dexp(-0.5d0*(wex*x1relright2))
             if((phileft.le.gausseps).and.(phiright.le.gausseps)) exit
             phinxpart = phinxpart + phileft + phiright
             dphinxpart = dphinxpart + x1relleft*phileft + x1relright*phiright
             d2phinxpart = d2phinxpart + x1relleft2*phileft + x1relright2*phiright
             d3phinxpart = d3phinxpart + x1relleft*x1relleft2*phileft &
     &            + x1relright*x1relright2*phiright
             d4phinxpart = d4phinxpart + x1relleft2*x1relleft2*phileft &
     &            + x1relright2*x1relright2*phiright
!             if (icell.gt.2) then
!                write (6,*) 'Warning: in deriv_2dgauss_periodic gaussians are wider than 2 unit cells.'
!             endif
          enddo

          phin(ib,ie)=phinypart*phinxpart

          dphin(1,ib,ie)=-wex*dphinxpart*phinypart
          dphin(2,ib,ie)=-wey*x2rel*phin(ib,ie)

          d2phin(ib,ie)=wex2*d2phinxpart*phinypart + (wey2*x2rel2 - wex - wey)*phin(ib,ie)

!          tempd2 = (wex2*x1rel2 + wey2*x2rel2 - wex - wey)

          tempd4parm = 0.5d0*(0.5d0*xg4i - we*x2rel2)


! parameter derivatives:
          dparam(1,ib,ie)=-dphin(1,ib,ie)                                  ! wrt xg1
          dparam(2,ib,ie)=-dphin(2,ib,ie)                                  ! wrt xg2
          dparam(3,ib,ie)= phinypart*(0.5d0*(0.5d0*xg3i*phinxpart - we*d2phinxpart))     ! wrt xg3
          dparam(4,ib,ie)= tempd4parm * phin(ib,ie)      ! wrt xg4


          d2param(1,1,ib,ie)=wex2*d2phinxpart*phinypart - wex*phin(ib,ie)  ! wrt xg1,xg1
          d2param(2,2,ib,ie)=(wey2*x2rel2-wey)*phin(ib,ie)                 ! wrt xg2,xg2
          d2param(3,3,ib,ie)=(0.25d0*we*(we*d4phinxpart - xg3i*d2phinxpart) &
     &                        - 0.1875d0*xg3i*xg3i*phinxpart)*phinypart    ! wrt xg3,xg3
          d2param(4,4,ib,ie)=(0.25d0*we*x2rel2*(we*x2rel2 - xg4i) &
     &                        - 0.1875d0*xg4i*xg4i)*phin(ib,ie)            ! wrt xg4,xg4


          d2param(1,2,ib,ie)=wey*x2rel*dparam(1,ib,ie)                     ! wrt xg1,xg2
          d2param(1,3,ib,ie)=we*phinypart*(1.25d0*dphinxpart - 0.5d0*d3phinxpart)    ! wrt xg1,xg3
          d2param(1,4,ib,ie)=tempd4parm*dparam(1,ib,ie)     ! wrt xg1,xg3

          d2param(2,3,ib,ie)=x2rel*wey*dparam(3,ib,ie)                     ! wrt xg2,xg3
          d2param(2,4,ib,ie)=we*x2rel*(1.25d0 - 0.5d0*wey*x2rel2)*phin(ib,ie)  ! wrt xg2,xg4

          d2param(3,4,ib,ie)=tempd4parm*dparam(3,ib,ie) ! wrt xg3,xg4

          d2param(2,1,ib,ie)=d2param(1,2,ib,ie)
          d2param(3,1,ib,ie)=d2param(1,3,ib,ie)
          d2param(4,1,ib,ie)=d2param(1,4,ib,ie)
          d2param(3,2,ib,ie)=d2param(2,3,ib,ie)
          d2param(4,2,ib,ie)=d2param(2,4,ib,ie)
          d2param(4,3,ib,ie)=d2param(3,4,ib,ie)

! coo. derivatives of parameter derivatives:
          ddparam(1,1,ib,ie)=-d2param(1,1,ib,ie)                            ! wrt x1,xg1
          ddparam(2,1,ib,ie)=-d2param(2,1,ib,ie)                            ! wrt x2,xg1

          ddparam(1,2,ib,ie)= ddparam(2,1,ib,ie)                            ! wrt x1,xg2
          ddparam(2,2,ib,ie)=-d2param(2,2,ib,ie)                             ! wrt x2,xg2

          ddparam(1,3,ib,ie)=-d2param(1,3,ib,ie)                             ! wrt x1,xg3
          ddparam(2,3,ib,ie)=-d2param(2,3,ib,ie)                             ! wrt x2,xg3

          ddparam(1,4,ib,ie)=-d2param(1,4,ib,ie)                             ! wrt x1,xg4
          ddparam(2,4,ib,ie)=-d2param(2,4,ib,ie)                             ! wrt x2,xg4

          d2dparam(1,ib,ie)= wex*phinypart*(dphinxpart*(-3.d0*wex - wey - wey2*x2rel2) &
     &                              + wex2*d3phinxpart)              ! laplacian of dparam(1,ib,ie)
          d2dparam(2,ib,ie)=wey*(2*dparam(2,ib,ie)+x2rel*d2phin(ib,ie))   ! laplacian of dparam(2,ib,ie)
          d2dparam(3,ib,ie)=phinypart*( -(1.25d0*we + 0.25d0*wey*xg3i - 0.25d0*x2rel2*wey2)*phinxpart &
     &                          + we*(2.75d0*wex - 0.5d0*wey + 0.5d0*x2rel2*wey2)*d2phinxpart &
     &                          - 0.5d0*we*wex2)             ! laplacian of dparam(3,ib,ie)
          d2dparam(4,ib,ie)=tempd4parm*d2phin(ib,ie) + we*(2.0d0*wey*x2rel2 - 1.0d0)*phin(ib,ie)
                                                ! laplacian of dparam(4,ib,ie)

        enddo
      enddo

      return
      end
!-------------------------------------------------------------------------

      subroutine basis_fns_2dgauss_periodic2(iel,rvec_en,r_en)

! Written by Abhijit C. Mehta, July 2011
! (modified version of basis_fns_polargauss)
! 2-dimensional localized gaussian basis set,
!    gaussians can have diffferent widths in x- and y- directions
!    periodic in x-direction (designed for iperiodic = 1)
! Main purpose is the study 2d wigner crystal in quantum wires
!    (infinite periodic wire - i.e., 1D periodic boundary conditions)
!  This version uses an exp(cos(2*pi*x/alattice)) term instead of
!      starting with the exp(-x^2) term and summing images


! arguments: iel=0 -> all electron
!               >0 -> only electron iel
!            rvec_en=vector electron-nucleus
!                    (or electron-dot/"wire" center in this context)
!                    (i.e., distance from (0,0))

! output: phin,dphin, and d2phin are calculated

! Wave functions are given by (except the normalization sqrt(we/pi)):

! phi=dsqrt(we) * exp(-we*xg4/2*(x2-xg2)^2) * exp(xg3*cos(xt-xt0))
!    (xt - xt0) = 2*pi*(x1-xg1)/alattice

! where x1,x2 are the electronic x and y positions,
!       xg1, xg2 are the gaussian x and y positions.
!       xg3, xg4 are the gaussian x and y width parameters
!       (Note that xg3 and xg4 are the "omegas" (1/width)
!         in units of we (or wire_w in case of wire)))


      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      use periodic_1d_mod

      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
      else
        nelec1=iel
        nelec2=iel
      endif

      ic=1
      expnorm = 1.d0
!      write(6,*) 'in basis_fns'

      do ie=nelec1,nelec2
!       modulo math was done in distances.f, so we don't need to do it here
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)

!       write(6,*) 'x1,x2,we=',x1,x2,we

!       assume that input parameters oparm are scaled to periodic BC's
        do ib=1,nbasis
          xg3 = oparm(3,ib,iwf)
          wey=we*oparm(4,ib,iwf)
          wey2=wey*wey
          x1rel=x1-oparm(1,ib,iwf)  !no need to do mod math bc of cos
          cxtrel=dcos(2.*pi*x1rel/alattice)
          sxtrel=dsin(2.*pi*x1rel/alattice)
          x2rel=x2-oparm(2,ib,iwf)
          x2rel2=x2rel*x2rel

          phinypart=dsqrt(we)*dexp(-0.5d0*wey*x2rel2)
          phinxpart=dexp(xg3*(cxtrel - expnorm))

          phin(ib,ie)=phinypart*phinxpart


!         write(6,*) 'ib,ie,phin(ib,ie)=',ib,ie,phin(ib,ie)
!         write(6,*) 'oparm1,oparm2,oparm3,oparm4=',oparm(1,ib,iwf),oparm(2,ib,iwf),oparm(3,ib,iwf),oparm(4,ib,iwf)

          dphin(1,ib,ie)=-(xg3*2.*pi*sxtrel/alattice)*phin(ib,ie)
          dphin(2,ib,ie)=-wey*x2rel*phin(ib,ie)

          d2dx2term = xg3*(2.*pi/alattice)**2*(xg3*sxtrel**2-cxtrel)
          d2dy2term = wey2*x2rel2 - wey

          d2phin(ib,ie) = (d2dx2term + d2dy2term)*phin(ib,ie)

        enddo
      enddo

      return
      end

!--------------------------------------------------------------------------

      subroutine deriv_2dgauss_periodic2(rvec_en,r_en)

! Written by Abhijit C. Mehta, July 2011
! (modified version of deriv_polargauss)
! 2-dimensional localized gaussian basis set, and derivatives
!    gaussians can have diffferent widths in x- and y- directions
!    periodic in x-direction (designed for iperiodic = 1)
! Main purpose is the study 2d wigner crystal in quantum wires
!    (infinite periodic wire - i.e., 1D periodic boundary conditions)
!  This version uses an exp(cos(2*pi*x/alattice)) term instead of
!      starting with the exp(-x^2) term and summing images


! arguments: iel=0 -> all electron
!               >0 -> only electron iel
!            rvec_en=vector electron-nucleus
!                    (or electron-dot/"wire" center in this context)
!                    (i.e., distance from (0,0))

! output: phin,dphin, and d2phin are calculated

! Wave functions are given by (except the normalization sqrt(we/pi)):

! phi=dsqrt(we) * exp(-we*xg4/2*(x2-xg2)^2) * exp(xg3*cos(xt-xt0))
!    (xt - xt0) = 2*pi*(x1-xg1)/alattice

! where x1,x2 are the electronic x and y positions,
!       xg1, xg2 are the gaussian x and y positions.
!       xg3, xg4 are the gaussian x and y width parameters
!       (Note that xg3 and xg4 are the "omegas" (1/width)
!         in units of we (or wire_w in case of wire)))


      use coefs_mod
      use const_mod
      use wfsec_mod
      use phifun_mod
      use orbpar_mod
      use deriv_phifun_mod
      use periodic_1d_mod

      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      dimension rvec_en(3,nelec,*),r_en(nelec,*)

      nelec1=1
      nelec2=nelec

      ic=1
      expnorm = 1.d0
!      write(6,*) 'in deriv_2dgauss'

      do ie=nelec1,nelec2
!       modulo math was done in distances.f, so we don't need to do it here
        x1=rvec_en(1,ie,ic)
        x2=rvec_en(2,ie,ic)

!       write(6,*) 'x1,x2,we=',x1,x2,we

!       assume that input parameters oparm are scaled to periodic BC's
        do ib=1,nbasis
          xg3 = oparm(3,ib,iwf)
          wey=we*oparm(4,ib,iwf)
          wey2=wey*wey
          x1rel=x1-oparm(1,ib,iwf)  !no need to do mod math bc of cos
          cxtrel=dcos(2.*pi*x1rel/alattice)
          sxtrel=dsin(2.*pi*x1rel/alattice)
          x2rel=x2-oparm(2,ib,iwf)
          x2rel2=x2rel*x2rel

          temp2pida = 2.*pi/alattice

          phinypart=dsqrt(we)*dexp(-0.5d0*wey*x2rel2)
          phinxpart=dexp(xg3*(cxtrel - expnorm))

          phin(ib,ie)=phinypart*phinxpart


!         write(6,*) 'ib,ie,phin(ib,ie)=',ib,ie,phin(ib,ie)
!         write(6,*) 'oparm1,oparm2,oparm3,oparm4=',oparm(1,ib,iwf),oparm(2,ib,iwf),oparm(3,ib,iwf),oparm(4,ib,iwf)

          dphin(1,ib,ie)=-xg3*temp2pida*sxtrel*phin(ib,ie)
          dphin(2,ib,ie)=-wey*x2rel*phin(ib,ie)

          d2dx2term = xg3*(temp2pida*temp2pida)*(xg3*sxtrel**2-cxtrel)
          d2dy2term = wey2*x2rel2 - wey

          d2phin(ib,ie) = (d2dx2term + d2dy2term)*phin(ib,ie)

! parameter derivatives:
          dparam(1,ib,ie)=-dphin(1,ib,ie)                ! wrt xg1
          dparam(2,ib,ie)=-dphin(2,ib,ie)                ! wrt xg2
          dparam(3,ib,ie)=(cxtrel - expnorm)*phin(ib,ie) ! wrt xg3
          dparam(4,ib,ie)=-0.5d0*we*x2rel2*phin(ib,ie)   ! wrt xg4

          d2param(1,1,ib,ie)=d2dx2term*phin(ib,ie)  ! wrt xg1,xg1
          d2param(2,2,ib,ie)=d2dy2term*phin(ib,ie)  ! wrt xg2,xg2
          d2param(3,3,ib,ie)=(cxtrel - expnorm)*dparam(3,ib,ie)  ! wrt xg3,xg3
          d2param(4,4,ib,ie)=-0.5d0*we*x2rel2*dparam(4,ib,ie)    ! wrt xg4,xg4

          d2param(1,2,ib,ie)=wey*x2rel*dparam(1,ib,ie)  ! wrt xg1,xg2
          d2param(1,3,ib,ie)=-temp2pida*sxtrel*(phin(ib,ie) + xg3*dparam(3,ib,ie))   ! wrt xg1,xg3
          d2param(1,4,ib,ie)=-0.5d0*we*x2rel2*dparam(1,ib,ie)    ! wrt xg1,xg3

          d2param(2,3,ib,ie)=wey*x2rel*dparam(3,ib,ie)          ! wrt xg2,xg3
          d2param(2,4,ib,ie)=x2rel*(wey*dparam(4,ib,ie) + we*phin(ib,ie))  ! wrt xg2,xg4

          d2param(3,4,ib,ie)=-0.5d0*we*x2rel2*dparam(3,ib,ie) ! wrt xg3,xg4

          d2param(2,1,ib,ie)=d2param(1,2,ib,ie)
          d2param(3,1,ib,ie)=d2param(1,3,ib,ie)
          d2param(4,1,ib,ie)=d2param(1,4,ib,ie)
          d2param(3,2,ib,ie)=d2param(2,3,ib,ie)
          d2param(4,2,ib,ie)=d2param(2,4,ib,ie)
          d2param(4,3,ib,ie)=d2param(3,4,ib,ie)

! coo. derivatives of parameter derivatives:
          ddparam(1,1,ib,ie)=-d2param(1,1,ib,ie)        ! wrt x1,xg1
          ddparam(2,1,ib,ie)=-d2param(2,1,ib,ie)        ! wrt x2,xg1

          ddparam(1,2,ib,ie)= ddparam(2,1,ib,ie)        ! wrt x1,xg2
          ddparam(2,2,ib,ie)=-d2param(2,2,ib,ie)        ! wrt x2,xg2

          ddparam(1,3,ib,ie)=-d2param(1,3,ib,ie)        ! wrt x1,xg3
          ddparam(2,3,ib,ie)=-d2param(2,3,ib,ie)        ! wrt x2,xg3

          ddparam(1,4,ib,ie)=-d2param(1,4,ib,ie)        ! wrt x1,xg4
          ddparam(2,4,ib,ie)=-d2param(2,4,ib,ie)        ! wrt x2,xg4

!  parameter derivatives of laplacian:
          d2dparam(1,ib,ie)=(d2dx2term + d2dy2term)*dparam(1,ib,ie) &
     &     - xg3*(temp2pida**3)*sxtrel*(2.*xg3*cxtrel + 1.)*phin(ib,ie)
          d2dparam(2,ib,ie)=(d2dx2term + d2dy2term)*dparam(2,ib,ie) &
     &     - 2.*wey2*x2rel*phin(ib,ie)
          d2dparam(3,ib,ie)=(d2dx2term + d2dy2term)*dparam(3,ib,ie) &
     &     + temp2pida*temp2pida*(2.*xg3*sxtrel**2 - cxtrel)*phin(ib,ie)
          d2dparam(4,ib,ie)=(d2dx2term + d2dy2term)*dparam(4,ib,ie) &
     &     + we*(2.*wey*x2rel2 - 1.)*phin(ib,ie)

        enddo
      enddo

      return
      end


!--------------------------------------------------------------------------




!**RM(3) Following added by RM.
!*******************
      Subroutine calculate_spherical_harmonics
!*******************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - P_lm(x,y,z), dP_lm/dx, dP_lm/dy, and dP_lm/dz is calculated.
!       where, Y_lm = P_lm(x,y,z)/(r^l).
!
!     - Inputs (x,y,z) are already stored as module variables
!       x_power_of(N) etc. by subroutine basis_fns.
!
!     - P_lm etc are handled in the form of
!
!         P_lm(x,y,z) = [\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}]/r^L,
!         where the summation is taken under
!             l_s + m_s + n_s = L .
!         Then.
!             s = 1 : s_max,
!                     s_max =  (L+1)(L+2)/2
!
!         These quantities (coefficient and exponents) are already stored
!         as module variables large_q(N_lms,1:4)
!
!                large_q(N_LMs,1) = a_s
!                large_q(N_LMs,2) = l_s
!                large_q(N_LMs,3) = m_s
!                large_q(N_LMs,4) = n_s
!
!         These quantities are already stored as module variables
!         large_q(N_lms,1:4) by subroutine setup_spherical_harmonics.
!
!         The convention of the index N_lms is the followings:
!
!            N_LMs     L     M     s
!              1       0     0     1
!
!              2       1    -1     1
!              3       1    -1     2
!              4       1    -1     3
!
!              5       1     0     1
!              6       1     0     2
!              7       1     0     3
!
!              8       1    +1     1
!              9       1    +1     2
!             10       1    +1     3
!             ...    ...   ...   ...
!             Nmax  lmax  lmax  smax(lmax)
!
!             --> Nmax = (lmax + 1)*(3*lmax^3 + 17*lmax^2 + 28*lmax + 12)/12
!
!     - Results are stored as module variables
!
!           large_p  (l,m) = P_lm(x,y,z)
!         d_large_p(1,l,m) = dP_lm(x,y,z)/dx      etc.
!-------------------------------------------------
      use real_spherical_harmonics

      implicit none

      integer :: itemp1
      integer :: scount,mcount,lcount,m_tempmax,s_tempmax

!      call show_large_q(2,0)  !* Debug use.

      large_p(:,:) = 0.d0 ; d_large_p(:,:,:) = 0.d0

      lcount = 0
      m_tempmax = 2*lcount + 1
      s_tempmax = (lcount+1)*(lcount+2)/2
      scount = 1 ; mcount = 1

      do itemp1 = 1, nmax

         large_p(lcount,-lcount+mcount-1) = large_p(lcount,-lcount+mcount-1) &
     &        + dble(large_q(itemp1,1))*x_power_of(large_q(itemp1,2)) &
     &                                 *y_power_of(large_q(itemp1,3)) &
     &                                 *z_power_of(large_q(itemp1,4))

         d_large_p(1,lcount,-lcount+mcount-1) = d_large_p(1,lcount,-lcount+mcount-1) &
     &        + dble(large_q(itemp1,1)*large_q(itemp1,2)) &
     &                                *x_power_of(large_q(itemp1,2)-1) &
     &                                *y_power_of(large_q(itemp1,3)) &
     &                                *z_power_of(large_q(itemp1,4))

         d_large_p(2,lcount,-lcount+mcount-1) = d_large_p(2,lcount,-lcount+mcount-1) &
     &        + dble(large_q(itemp1,1)*large_q(itemp1,3)) &
     &                                *x_power_of(large_q(itemp1,2)) &
     &                                *y_power_of(large_q(itemp1,3)-1) &
     &                                *z_power_of(large_q(itemp1,4))

         d_large_p(3,lcount,-lcount+mcount-1) = d_large_p(3,lcount,-lcount+mcount-1) &
     &        + dble(large_q(itemp1,1)*large_q(itemp1,4)) &
     &                                *x_power_of(large_q(itemp1,2)) &
     &                                *y_power_of(large_q(itemp1,3)) &
     &                                *z_power_of(large_q(itemp1,4)-1)

         if(scount==s_tempmax)then
            scount = 0
            if(mcount==m_tempmax)then
               lcount = lcount + 1
               m_tempmax = 2*lcount + 1
               s_tempmax = (lcount+1)*(lcount+2)/2
               mcount = 0
            endif
            mcount = mcount + 1
         endif
         scount = scount + 1
      enddo

      end subroutine calculate_spherical_harmonics


!**************************************
      Subroutine setup_coefficients_ylm
!**************************************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Definition of the normalization coeficients is
!
!       Ylm = (1/r^l)*[\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}]*coef_ylm(l,m)
!
!       Here coef_ylm(l,m) is beyond the factor sqrt((2*l+1)/(4*pi)), then
!
!         coef_ylm(0,    0) = 1        ,
!         coef_ylm(1,\pm 1) = 1        ,
!         coef_ylm(1,    0) = 1        ,
!         coef_ylm(2,\pm 2) = sqrt(3)/2,
!         coef_ylm(2,\pm 1) = sqrt(3)  ,
!         coef_ylm(2,    0) = 1/2      ,... etc.
!
!       See, for exam.,
!       http://www.cachesoftware.com/mopac/Mopac2002manual/real_spherical_harmonics.html
!       for first several coeffs.
!
!     - coef_ylm(l,m) is given in the form as
!
!         = sqrt( [sq_coef_a_num(l)/sq_coef_a_denum(l)] * sq_coef_r_num(l,m) )
!
!       All 'sq_...' quantities are integer and already stored as module variables
!       in the preceding call
!         'read_input.f -->  call setup_spherical_harmonics'
!-------------

      use real_spherical_harmonics
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
      do l = 0,lmax
         do m = 0,l
            coef_ylm(l,m) = sqrt(dble(sq_coef_a_num(l)) &
     &             /(sq_coef_a_denum(l)*sq_coef_r_denum(l,m)))
!* For debug:
!           coef_ylm(l,m) = sqrt(1.d0/dble(sq_coef_r_denum(l,m)))
!           !* ^^ This should be 1.d0 for (l,m) = (l,l-1).
!*                Use this for debug.
         enddo
      enddo

      coef_ylm(1,0) = 1.d0       !* The only exception...
                                 !  This should be substituted
                                 !  after the recursion generation is performed.
      end subroutine setup_coefficients_ylm


!**************
      Subroutine setup_spherical_harmonics
!*************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Generating higher real spherical harmonics using recursion
!       relations:
!
!         Q(L,L-1) = Q(L-1,L-1)*z
!
!         Q(L,+(M+1)) = A(L)*R(L,M,+)*[(jx)*Q(L,-M) + (jy)*Q(L,+M)]
!         Q(L,+(M-1)) = A(L)*R(L,M,-)*[(jx)*Q(L,-M) - (jy)*Q(L,+M)]
!         Q(L,-(M-1)) = A(L)*R(L,M,-)*[(jx)*Q(L,+M) + (jy)*Q(L,-M)]
!         Q(L,-(M+1)) = A(L)*R(L,M,+)*[(jx)*Q(L,+M) - (jy)*Q(L,-M)]
!
!         where, Y_lm = P_lm(x,y,z)/(r^l),
!                       P_lm(x,y,z) =: Q(L,M) = 'large_q' array.
!
!     - See also the description of calculate_spherical_harmonics
!       for conventions.
!
!     - First several Q(L,M) upto p-wave is set manually. Then
!       higher ones are constructed recursively by calling
!       subroutine 'recursion_body' as
!
!       Do L=2, L_max
!         call recursion_body[ Q(L-1,\pm(L-1)) ]
!       Enddo
!
!     - Components for the normalization coefficients
!       sq_coef_a_num etc. are also stored here.
!       sq_coef_a_... are determined only by l
!       and constructed recursively here, whilst
!       sq_coef_r_... is depending on m as well
!       so this is constructed under the
!       'call recursion_body'
!----------------
      use real_spherical_harmonics

      implicit none
      integer :: ialloc
      integer :: l,smax
      integer :: l_alreadyset,s_alreadyset,n_alreadyset
      integer :: pos_fin,pos_init,neg_init,neg_fin
      integer :: itemp1

      allocate(sq_coef_a_num(0:lmax),sq_coef_a_denum(0:lmax),stat=ialloc)
      sq_coef_a_num(:)=0 ; sq_coef_a_denum(:)=1
      allocate(sq_coef_r_denum(0:lmax,0:lmax),stat=ialloc)
      sq_coef_r_denum(:,:)=1

      nmax = (lmax + 1)*(3*lmax**3 + 17*lmax**2 + 28*lmax + 12)/12

      allocate(large_q(1:nmax,1:4),stat=ialloc) ; large_q(:,:)=0

! * s- and p-wave is set manually.
      sq_coef_a_num(0) = 1 ; sq_coef_a_denum(0) = 1
      large_q(1,1) = 1 ; large_q(1,2) = 0 ; large_q(1,3) = 0 ;large_q(1,4) = 0 ! S( 0)

      sq_coef_a_num(1) = 2 ; sq_coef_a_denum(1) = 1
      sq_coef_r_denum(1,1) = 2
      large_q(2,1) = 1 ; large_q(2,2) = 0 ; large_q(2,3) = 1 ;large_q(2,4) = 0 ! P(-1)
      large_q(3:4,:) = 0

      large_q(5,1) = 1 ; large_q(5,2) = 0 ; large_q(5,3) = 0 ;large_q(5,4) = 1 ! P( 0)
      large_q(6:7,:) = 0

      large_q(8,1) = 1 ; large_q(8,2) = 1 ; large_q(8,3) = 0 ;large_q(8,4) = 0 ! P(+1)
      large_q(9:10,:) = 0

      l_alreadyset = 1

! * recursion generation
      do itemp1 = 2,lmax
         n_alreadyset = (l_alreadyset + 1)*(3*l_alreadyset**3 + 17*l_alreadyset**2 + &
     &                   28*l_alreadyset + 12)/12
         s_alreadyset = (l_alreadyset+1)*(l_alreadyset+2)/2

         pos_fin  = n_alreadyset
         pos_init = pos_fin  - s_alreadyset + 1
         neg_init = pos_init - s_alreadyset*(2*l_alreadyset)
         Neg_fin  = neg_init + s_alreadyset - 1

         !* Q_pos := Q(L-1,L-1) is included as large_q(s=pos_init:pos_fin)
         !* Q_neg := Q(L-1,1-L)                large_q(s=neg_init:neg_fin)

! * For given L, the size of Q is determined.
         l=l_alreadyset +1 ; smax = (l+1)*(l+2)/2

         sq_coef_a_num(l)   = sq_coef_a_num(l-1)*(2*l-1)
         sq_coef_a_denum(l) = sq_coef_a_denum(l-1)*(2*l-2)

         call recursion_body(large_q(pos_init:pos_fin,:),large_q(neg_init:neg_fin,:), &
     &                       l_alreadyset,s_alreadyset,n_alreadyset)
         l_alreadyset = l_alreadyset + 1
      enddo

!     call show_large_q(lmax,1) !* Debug use.

      END subroutine setup_spherical_harmonics


!**************
      SUBROUTINE recursion_body(theta_pos,theta_neg,l_alreadyset,s_alreadyset, &
     &n_alreadyset)
!*************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Generating higher real spherical harmonics using recursion
!       relations:
!
!         Q(L,L-1) = Q(L-1,L-1)*z
!
!         Q(L,+(M+1)) = A(L)*R(L,M,+)*[(jx)*Q(L,-M) + (jy)*Q(L,+M)]
!         Q(L,+(M-1)) = A(L)*R(L,M,-)*[(jx)*Q(L,-M) - (jy)*Q(L,+M)]
!         Q(L,-(M-1)) = A(L)*R(L,M,-)*[(jx)*Q(L,+M) + (jy)*Q(L,-M)]
!         Q(L,-(M+1)) = A(L)*R(L,M,+)*[(jx)*Q(L,+M) - (jy)*Q(L,-M)]
!
!         where, Y_lm = P_lm(x,y,z)/(r^l),
!                       P_lm(x,y,z) =: Q(L,M) = 'large_q' array.
!
!       So, for recursion construction (L-1) --> (L) call this
!       subroutine with specifying,
!
!             l_alreadyset = (L-1)
!             theta_pos    = Q(L-1,L-1)
!             theta_neg    = Q(L-1,1-L)
!
!     - See also the description of calculate_spherical_harmonics
!       for conventions.
!
!     - Within this subroutine the number L is fixed, so we use
!       allocatable 'q_array: q(M,s,1:4)' defined as
!
!                Convention of the array Q and q are
!                Q(N_LMs,1) = a_s =: q(M,s,1)
!                Q(N_LMs,2) = l_s =: q(M,s,2)
!                Q(N_LMs,3) = m_s =: q(M,s,3)
!                Q(N_LMs,4) = n_s =: q(M,s,4)
!       note,
!                Q(N_LMs,1:4): large_q(N_LMs,1:4) : module variable
!                q(  M,s,1:4):       q(  M,s,1:4) : local  variable
!-----------

      USE real_spherical_harmonics

      implicit none
      integer,allocatable :: q(:,:,:)
      integer,allocatable :: qtemp1(:,:),qtemp2(:,:),qtemp3(:,:), &
     &qtemp4(:,:)
      integer :: ialloc,icount1,icount2
      integer :: l,smax
      integer :: l_alreadyset,s_alreadyset,n_alreadyset
      integer :: theta_pos(s_alreadyset,1:4),theta_neg(s_alreadyset,1:4)
      integer :: itemp1,i
      integer :: ninit,nfin
      integer :: mcount

! * For given L,
!   Y(L,M) = [\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}]/r^L
!   q(M,1:smax,...) = [\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}]
!   smax is determined as smax =  (L+1)(L+2)/2

      l=l_alreadyset +1 ; smax = (l+1)*(l+2)/2

      allocate(q(-l:+l,1:smax,1:4),stat=ialloc) ; q(:,:,:)=0
      allocate(qtemp1(1:smax,1:4),qtemp2(1:smax,1:4), &
     &qtemp3(1:smax,1:4),qtemp4(1:smax,1:4),stat=ialloc)
      qtemp1(:,:)=0 ; qtemp2(:,:)=0 ; qtemp3(:,:)=0 ; qtemp4(:,:)=0


! * Setup the 'seed' for recursion Q(L-1,L-1) and Q(L-1,1-L)
!   and then Q(L,L-1) = Q(L-1,L-1)*z
!            Q(L,1-L) = Q(L-1,1-L)*z
      q(l-1,1:s_alreadyset,1:4) = theta_pos(1:s_alreadyset,1:4)
      q(l-1,1:s_alreadyset,4)   = q(l-1,1:s_alreadyset,4) + 1
      q(1-l,1:s_alreadyset,1:4) = theta_neg(1:s_alreadyset,1:4)
      q(1-l,1:s_alreadyset,4)   = q(1-l,1:s_alreadyset,4) + 1


! * Then up/down m, starting from \pm m, [m,-m] --> [(m+1),(m-1)]
!   Q(L,+(M+1)) = A(L)*R(L,M,+)*[(jx)*Q(L,-M) + (jy)*Q(L,+M)]
!   Q(L,+(M-1)) = A(L)*R(L,M,-)*[(jx)*Q(L,-M) - (jy)*Q(L,+M)]

      mcount = l-1
      call jxq(q(-mcount,:,:),smax,qtemp1)     != (jx)*Q(L,-M)
      call jyq(q(+mcount,:,:),smax,qtemp2)     != (jy)*Q(L,+M)
      call sum_q(qtemp1,qtemp2,+1,smax,qtemp3) != (jx)*Q(L,-M)+(jy)*Q(L,+M)
      call sum_q(qtemp1,qtemp2,-1,smax,qtemp4) != (jx)*Q(L,-M)-(jy)*Q(L,+M)
      icount1 = 1 ; icount2 = 1
      do i=1,smax
                                 !* Substituting into Q(L,+(M+1))
         if(qtemp3(i,1)/=0)then
            q(mcount+1,icount1,1)  =-qtemp3(i,1)
            q(mcount+1,icount1,2:4)= qtemp3(i,2:4)
            icount1 = icount1 + 1
         endif
                                 !* Substituting into Q(L,+(M-1))
         if(qtemp4(i,1)/=0)then
            q(mcount-1,icount2,1)  =-qtemp4(i,1)
            q(mcount-1,icount2,2:4)= qtemp4(i,2:4)
            icount2 = icount2 + 1
         endif
      enddo

! * Squared coefficients regarding to up/dn operation is stored.
      sq_coef_r_denum(l,mcount+1) = sq_coef_r_denum(l,mcount) &
     &                             *(l-mcount)*(l+mcount+1)
      sq_coef_r_denum(l,mcount-1) = sq_coef_r_denum(l,mcount) &
     &                             *(l+mcount)*(l-mcount+1)

! This factor is required for 'real' sperical harmonics taking into account the multiplicity +0 and -0.
      if(mcount==1) sq_coef_r_denum(l,mcount-1)=sq_coef_r_denum(l,mcount-1)*2


! * Then generate negative M part from \pm m, [m,-m] --> [-(m+1),-(m-1)]
!    Q(L,-(M-1)) = A(L)*R(L,M,-)*[(jx)*Q(L,+M) + (jy)*Q(L,-M)]
!    Q(L,-(M+1)) = A(L)*R(L,M,+)*[(jx)*Q(L,+M) - (jy)*Q(L,-M)]

      qtemp1(:,:)=0 ; qtemp2(:,:)=0 ; qtemp3(:,:)=0 ; qtemp4(:,:)=0
      call jxq(q(+mcount,:,:),smax,qtemp1)     != (jx)*Q(L,+M)
      call jyq(q(-mcount,:,:),smax,qtemp2)     != (jy)*Q(L,-M)
      call sum_q(qtemp1,qtemp2,+1,smax,qtemp3) != (jx)*Q(L,+M)+(jy)*Q(L,-M)
      call sum_q(qtemp1,qtemp2,-1,smax,qtemp4) != (jx)*Q(L,+M)-(jy)*Q(L,-M)

      icount1 = 1 ; icount2 = 1
      do i=1,smax
                                 !* Substituting into Q(L,-(M+1))
         if(qtemp4(i,1)/=0)then
            q(-(mcount+1),icount1,1)  =+qtemp4(i,1)
            q(-(mcount+1),icount1,2:4)= qtemp4(i,2:4)
            icount1 = icount1 + 1
         endif
                                 !* Substituting into Q(L,-(M-1))
         if((mcount-1)/=0)then
            if(qtemp3(i,1)/=0)then
               q(-(mcount-1),icount1,1)  =+qtemp3(i,1)
               q(-(mcount-1),icount1,2:4)= qtemp3(i,2:4)
               icount1 = icount1 + 1
            endif
         endif
      enddo


! * Down until to m=0 similarly.
      do itemp1 = l-2,1,-1
           !*    [m] --> [m-1]
           !      Q(L,+(M-1)) = A(L)*R(L,M,-)*[(jx)*Q(L,-M) - (jy)*Q(L,+M)]
         qtemp1(:,:)=0 ; qtemp2(:,:)=0 ; qtemp3(:,:)=0 ; qtemp4(:,:)=0
         call jxq(q(-itemp1,:,:),smax,qtemp1)
         call jyq(q(+itemp1,:,:),smax,qtemp2)
         call sum_q(qtemp1,qtemp2,-1,smax,qtemp3)

           !*    [m] --> [-(m-1)]
           !     Q(L,-(M-1)) = A(L)*R(L,M,-)*[(jx)*Q(L,+M) + (jy)*Q(L,-M)]
         call jxq(q(+itemp1,:,:),smax,qtemp1)
         call jyq(q(-itemp1,:,:),smax,qtemp2)
         call sum_q(qtemp1,qtemp2,+1,smax,qtemp4)
!        call sum_q(qtemp1,qtemp2,+1,smax,qtemp3) !<-- for Debug
                                                 ! to check if the identity
                                                ! (M=-0) = 0 is realized.
         icount1 = 1 ; icount2 = 1
         do i=1,smax
            if(qtemp3(i,1)/=0)then
               q(itemp1-1,icount2,1)  =-qtemp3(i,1)
               q(itemp1-1,icount2,2:4)= qtemp3(i,2:4)
               icount2 = icount2 + 1
            endif
            if((itemp1-1)/=0)then
               if(qtemp4(i,1)/=0)then
                  q(-(itemp1-1),icount1,1)  =+qtemp4(i,1)
                  q(-(itemp1-1),icount1,2:4)= qtemp4(i,2:4)
                  icount1 = icount1 + 1
               endif
            endif
         enddo

         sq_coef_r_denum(l,itemp1-1) = sq_coef_r_denum(l,itemp1) &
     &                                 *(l+itemp1)*(l-itemp1+1)

         if(itemp1==1) sq_coef_r_denum(l,itemp1-1)=sq_coef_r_denum(l,itemp1-1)*2  !* Taking care of M= +0 and -0.

      enddo

      ninit = n_alreadyset

!* Store results into mother Q-array
      do itemp1 = -l, l
         nfin = ninit + smax
         large_q(ninit+1:nfin,1)   = q(itemp1,:,1)
         large_q(ninit+1:nfin,2:4) = q(itemp1,:,2:4)
         ninit = nfin
      enddo

      deallocate(q,qtemp1,qtemp2,qtemp3,qtemp4)
      END SUBROUTINE recursion_body


!******************************************
      SUBROUTINE sum_q(q1,q2,sign,stemp,q3)
!******************************************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Do the summation or subtraction of q-array in the convention of
!       q(s,1) = a_s ; q(s,2) = l_s ; q(s,3) = m_s ; q(s,4) = n_s
!       representing q(:,:) == [\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}].
!
!     - When sign = +1 is specified, the output qtemp is given as
!       qtemp = q1 + q2, while with sign = -1 then
!       qtemp = q1 - q2.
!
!     - stemp should be smax of q's to be applied.
!-----------------------------------------
      implicit none

      integer :: sign,stemp
      integer :: q1(1:stemp,1:4),q2(1:stemp,1:4),q3(1:stemp,1:4),q2_copy(1:stemp,1:4)

      integer :: i,j,icount
      integer :: l

      l=2
      q3(:,:)=0
      q2_copy(:,1)   = sign*q2(:,1)
      q2_copy(:,2:4) = q2(:,2:4)
      icount = 0

      do i = 1,stemp
         if(q1(i,1)/=0)then
            icount = icount + 1
            Q3(icount,:) = q1(i,:)
            do j = 1,stemp
               if((q2_copy(j,2)==q1(i,2)).and.(q2_copy(j,3)==q1(i,3)).and. &
     &             (q2_copy(j,4)==q1(i,4)))then
               q3(icount,1) = q3(icount,1) + q2_copy(j,1)
               q2_copy(j,1) = 0
            endif
         enddo
      endif
      enddo

      do i = 1,stemp
         if(q2_copy(i,1)/=0)then
            icount = icount + 1
            q3(icount,:) = q2_copy(i,:)
         endif
      enddo
      END SUBROUTINE sum_q


!*******************************
      SUBROUTINE jyq(qin,s,qout)
!*******************************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Do the (jy)-operation of q-array in the convention of
!       q(s,1) = a_s ; q(s,2) = l_s ; q(s,3) = m_s ; q(s,4) = n_s
!       representing q(:,:) == [\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}].
!       where, (jy) := z*dx - x*dz etc. and then, in term of q-array structure
!       of (A,L,M,N), the operation corresponds to
!                (jy) --> (A*L,L-1,M,N+1) - (A*N,L+1,M,N-1)
!
!     - The output qout is given as
!       qout = (jy)*qin .
!
!     - s should be smax of q's to be applied.
!----------------------------
      implicit none

      integer :: s,i
      integer :: qin(1:s,1:4),qout(1:s,1:4)
      integer :: qtemp1(1:s,1:4),qtemp2(1:s,1:4)

      qtemp1(:,:)=0 ; qtemp2(:,:)=0
      do i = 1,s
         if(qin(i,1)/=0)then
            qtemp1(i,1) = qin(i,1)*qin(i,2) !* z*dx term
            qtemp1(i,2) = qin(i,2) - 1      !* z*dx term
            qtemp1(i,3) = qin(i,3)          !* z*dx term
            qtemp1(i,4) = qin(i,4) + 1      !* z*dx term
            qtemp2(i,1) = qin(i,1)*qin(i,4) !* x*dz term
            qtemp2(i,2) = qin(i,2) + 1      !* x*dz term
            qtemp2(i,3) = qin(i,3)          !* x*dz term
            qtemp2(i,4) = qin(i,4) - 1      !* x*dz term
         endif
      enddo
      call sum_q(qtemp1,qtemp2,-1,s,qout)
      END SUBROUTINE jyq

!***********************************
      SUBROUTINE jxq(qin,stemp,qout)
!***********************************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Do the (jx)-operation of q-array in the convention of
!       q(s,1) = a_s ; q(s,2) = l_s ; q(s,3) = m_s ; q(s,4) = n_s
!       representing q(:,:) == [\sum_s{a_s*(x^l_s)*(y^m_s)*(z^n_s)}].
!       where, (jx) := y*dz - z*dy etc.,
!       and then, in term of q-array structure of (A,L,M,N),
!       the operation corresponds to
!            (jx) --> (A*N,L,M+1,N-1) - (A*M,L,M-1,N+1)
!
!     - The output qout is given as
!       qout = (jx)*qin .
!
!     - s should be smax of q's to be applied.
!----------------------------
      implicit none

      integer :: stemp,i
      integer :: qin(1:stemp,1:4),qout(1:stemp,1:4)
      integer :: qtemp1(1:stemp,1:4),qtemp2(1:stemp,1:4)

      qtemp1(:,:)=0 ; qtemp2(:,:)=0
      do i = 1,stemp
         if(qin(i,1)/=0)then
            qtemp1(i,1) = qin(i,1)*qin(i,4) !* y*dz term
            qtemp1(i,2) = qin(i,2)          !* y*dz term
            qtemp1(i,3) = qin(i,3) + 1      !* y*dz term
            qtemp1(i,4) = qin(i,4) - 1      !* y*dz term
            qtemp2(i,1) = qin(i,1)*qin(i,3) !* z*dy term
            qtemp2(i,2) = qin(i,2)          !* z*dy term
            qtemp2(i,3) = qin(i,3) - 1      !* z*dy term
            qtemp2(i,4) = qin(i,4) + 1      !* z*dy term
         endif
      enddo
      call sum_q(qtemp1,qtemp2,-1,stemp,qout)
      END SUBROUTINE jxq


!***********************************************
      SUBROUTINE show_large_q(l_show_max,n_mask)
!***********************************************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 29/Dec./05 ; Start coding.
!     - 09/Jan./06 ; First completed.
!
!     //* Description *//
!     - Write out large_q(n,1:4) array upto l = l_show_max.
!       (debugging tool).

!     - n_mask = 0 : all components are written
!     - n_mask = 1 : only non zero components are written.
!-------------------------
      USE real_spherical_harmonics

      implicit none

      integer :: l_show_max,n_show_max,n_mask
      integer :: itemp1
      integer :: scount,mcount,lcount,m_tempmax,s_tempmax

      lcount = 0
      n_show_max = (l_show_max + 1)*(3*l_show_max**3 + 17*l_show_max**2 &
     &              + 28*l_show_max + 12)/12
      m_tempmax = 2*lcount + 1
      s_tempmax = (lcount+1)*(lcount+2)/2
      scount = 1 ; mcount = 1

      do itemp1 = 1, n_show_max
         if(scount==1) write(6,*) '** N,',itemp1,'L=',lcount,'M=',-lcount+mcount-1
         if(large_q(itemp1,1)/=0.or.(n_mask==0))then
            write(6,*) large_q(itemp1,:)
         endif

         if(scount==s_tempmax)then
            scount = 0
            if(mcount==m_tempmax)then
               lcount = lcount + 1
               m_tempmax = 2*lcount + 1
               s_tempmax = (lcount+1)*(lcount+2)/2
               mcount = 0
            endif
            mcount = mcount + 1
         endif
         scount = scount + 1
      enddo
      END SUBROUTINE show_large_q

!**EndRM(3)

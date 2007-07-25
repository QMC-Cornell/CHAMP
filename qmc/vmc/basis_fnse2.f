      subroutine basis_fnse2(ie,rvec_en,r_en)
c Written by Cyrus Umrigar
c Calculate 3-dim localised basis functions for electron ie

c In input:
c n1s,n2s,...     > 0 : Slater basis
c nsa,npa,nda     = 0 : asymptotic basis
c n1s,nsa,...     < 0 : Gaussian basis
c Here:
c n_bas2(irb,ict) > 0 : Slater basis
c                 = 0 : asymptotic basis
c                 < 0 : Gaussian basis

      use all_tools_mod !JT

!**RM(7)
      use real_spherical_harmonics
      implicit real*8(a-h,o-z)
      real(dp) :: aux1
      integer :: itemp1
      common /contr_ylm/ irecursion_ylm
!**EndRM(7)

c     common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE)
     &,n5f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ zex2(MRWF,MCTYPE,MWF),n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),iwrwf2(MBASIS)
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      dimension wfv(3,MRWF),xc(3)
     &,th(0:ML_BAS,0:ML_BAS),ph(-ML_BAS:ML_BAS)

c Here we have additional normalization factors beyond those in basis_norm, viz., sqrt((2*l+1)/(4*pi))
      data cd0,cd1,cd2,cf0,cf1,cf2,cf3/0.5d0,1.73205080756888d0,0.866025403784439d0,
     &0.5d0,0.612372435695794d0,1.93649167310371d0,0.790569415042095d0/

      ider=0
c     do 40 ie=1,nelec
        ib=0
        do 40 ic=1,ncent
          ict=iwctype(ic)

          if(r_en(ie,ic).eq.0.d0) then
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

          if(numr.le.0) then

            do 10 irb=1,nrbas(ict)
              n=n_bas2(irb,ict)
              rn=abs(n)
              rm1=r**(rn-1)
c Slater r^(n-1)*Exp(-zeta*r)
              if(n.gt.0) then
                zr=zex2(irb,ict,iwf)*r
                ex=dexp(-zr)
               elseif(n.eq.0) then
c Warning: Asymptotic and Gaussian not yet tested.
c Asymptotic r^(rn-1)*Exp(-zeta*r), where rn=beta+1, beta=betaq/zeta-1, zeta=sqrt(-2*E_ion)?
                stop 'asymptotic not yet fully tested'
                rn=betaq/zex2(irb,ict,iwf)
                rm1=r**(rn-1)
               elseif(n.lt.0) then
c Gaussian  r^(n-1)*Exp(-zeta*r^2) (to be checked)
                zr=2*zex2(irb,ict,iwf)*r2
                ex=dexp(-0.5d0*zr)
              endif
   10         wfv(1,irb)=rm1*ex

           else

            rk=r
            do 20 irb=1,nrbas(ict)
              call splfit_bas(rk,irb,ict,iwf,wfv(1,irb),ider)
   20         if(wfv(1,irb).eq.0.d0) wfv(1,irb)=DBLMIN

          endif

c         write(6,'(''ict,nbasis_ctype'',9i5)') ict,nbasis_ctype(ict)

!**RM(8)
          if(irecursion_ylm.eq.0)  then
           do ib2=1,nbasis_ctype(ict)
!          do 40 ib2=1,nbasis_ctype(ict)  !* Modified 'do 40' --> 'do-enddo' for
!                                         !  code extension.
!**EndRM(8)
            xx=xc(1)
            yy=xc(2)
            zz=xc(3)
            xx2=xx*xx
            yy2=yy*yy
            zz2=zz*zz

c           xhat=xx*ri
c           yhat=yy*ri
            zhat=zz*ri

c Phi function

            ph(0)=1

            ph(1)=xx

            ph(-1)=yy

            ph(2)=xx2-yy2

            ph(-2)=2*xx*yy

c           ph(3)=(xx2-yy2)*xx-2*xx*yy2
            ph(3)=ph(2)*ph(1)-ph(-2)*ph(-1)

            ph(-3)=ph(-2)*ph(1)+ph(-1)*ph(2)

c Theta function

            th(0,0)=1

            th(1,0)=ri*zz

            th(1,1)=ri

c           th(2,0)=cd0*ri3*(3*zz**2-r2)
            th(2,0)=cd0*(3*zhat**2-1)

            th(2,1)=cd1*ri2*zz

            th(2,2)=cd2*ri2

            th(3,0)=cf0*ri3*zz*(2*zz2-3*(xx2+yy2))

            th(3,1)=cf1*ri3*(4*zz2-(xx2+yy2))

            th(3,2)=cf2*ri3*zz

            th(3,3)=cf3*ri3

            ib=ib+1
            irb=iwrwf(ib2,ict)
            l=l_bas(ib)
            m=m_bas(ib)
            mabs=abs(m)
            ylm=th(l,mabs)*ph(m)
            phin(ib,ie)=ylm*wfv(1,irb)

!**RM(9)
       enddo
      else
!*******************
!     New version for recursion construction of Ylm.
!*******************
!
!     Ryo MAEZONO / C.J. Umrigar
!
!     //* Revisions   *//
!     - 01/Feb./06 ; First completed.
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
!         'main_mov1.f -->  call setup_spherical_harmonics'
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

!* construct power of (x,y,z,1/r)
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

!* P_lm(x,y,z), dP_lm/dx, dP_lm/dy, and dP_lm/dz is calculated and stored.
         call calculate_spherical_harmonics

!* Orbital functions are constructed.
         do ib2 = 1,nbasis_ctype(ict)
            irb=iwrwf(ib2,ict)
            ib = ib + 1
            l = l_bas(ib)
            m = m_bas(ib)
            mabs=abs(m)
            if(l.gt.lmax) stop 'l>lmax in basis_fns.f'

            ylm           = large_p(l,m)*r_inv_power_minus(l)
            phin(ib,ie)   = ylm*wfv(1,irb)*coef_ylm(l,mabs)
           d2phin(ib,ie)  = coef_ylm(l,mabs)*ylm*(
     &                       wfv(3,irb) + 2.d0*wfv(2,irb)*r_inv_power_minus(1)
     &                       - dble(l*(l+1))*wfv(1,irb)*r_inv_power_minus(2))
           aux1           = (wfv(2,irb) - dble(l)*wfv(1,irb)*r_inv_power_minus(1))
     &                      *large_p(l,m)*r_inv_power_minus(1)
           dphin(1,ib,ie) = coef_ylm(l,mabs)*(aux1*xc(1) + d_large_p(1,l,m)*wfv(1,irb))
     &                      *r_inv_power_minus(l)
           dphin(2,ib,ie) = coef_ylm(l,mabs)*(aux1*xc(2) + d_large_p(2,l,m)*wfv(1,irb))
     &                      *r_inv_power_minus(l)
           dphin(3,ib,ie) = coef_ylm(l,mabs)*(aux1*xc(3) + d_large_p(3,l,m)*wfv(1,irb))
     &                      *r_inv_power_minus(l)
         enddo
      endif
!**EndRM(9)
   40 continue

      call object_modified_by_index (phin_index) !JT
      call object_modified_by_index (dphin_index) !JT
      call object_modified_by_index (d2phin_index) !JT

      return
      end
c-----------------------------------------------------------------------

      subroutine basis_fns_2de2(ie,rvec_en,r_en)
c Written by Cyrus Umrigar
c Calculate 2-dim localised basis functions for electron ie

c In input:
c n1s,n2s,...     > 0 : Slater basis
c nsa,npa,nda     = 0 : asymptotic basis
c n1s,nsa,...     < 0 : Gaussian basis
c Here:
c n_bas2(irb,ict) > 0 : Slater basis
c                 = 0 : asymptotic basis
c                 < 0 : Gaussian basis

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'numbas.h'
      include 'force.h'

c     common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE)
     &,n5f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ zex2(MRWF,MCTYPE,MWF),n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),iwrwf2(MBASIS)
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      dimension wfv(3,MRWF),xc(3)
     &,ph(-ML_BAS:ML_BAS)

      ider=0

c     do 40 ie=nelec1,nelec2
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

          if(numr.le.0) then

            do 10 irb=1,nrbas(ict)
              n=n_bas2(irb,ict)
              rn=abs(n)
c             rnm2=max(0.d0,rn-2)
              rm1=r**(rn-1)
c Slater r^(n-1)*Exp(-zeta*r)
              if(n.gt.0) then
                zr=zex2(irb,ict,iwf)*r
                ex=dexp(-zr)
               elseif(n.eq.0) then
c Warning: Asymptotic and Gaussian not yet tested.
c Asymptotic r^(rn-1)*Exp(-zeta*r), where rn=beta+1, beta=betaq/zeta-1, zeta=sqrt(-2*E_ion)?
                stop 'asymptotic not yet fully tested'
                rn=betaq/zex2(irb,ict,iwf)
                rm1=r**(rn-1)
               elseif(n.lt.0) then
c Gaussian  r^(n-1)*Exp(-zeta*r^2) (to be checked)
c               rnm2=rn-2
                zr=2*zex2(irb,ict,iwf)*r2
                ex=dexp(-0.5d0*zr)
              endif
   10         wfv(1,irb)=rm1*ex

           else

            rk=r
            do 20 irb=1,nrbas(ict)
   20         call splfit_bas(rk,irb,ict,iwf,wfv(1,irb),ider)

          endif

          do 40 ib2=1,nbasis_ctype(ict)

            xx=xc(1)
            yy=xc(2)
            xhat=xx*ri
            yhat=yy*ri

c Phi function

            ph(0)=1
            ph(1)=xhat
            ph(-1)=yhat
            do 30 i=2,ML_BAS
              ph(i)=ph(i-1)*ph(1)-ph(1-i)*ph(-1)
   30         ph(-i)=ph(1-i)*ph(1)+ph(i-1)*ph(-1)

            ib=ib+1
            irb=iwrwf(ib2,ict)
            m=m_bas(ib)
            ylm=ph(m)
            phin(ib,ie)=ylm*wfv(1,irb)
   40 continue

      return
      end

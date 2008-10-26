      subroutine orbitals_loc_ana(iel,rvec_en,r_en,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c Calculate localized orbitals and derivatives for all or 1 electrons

      use all_tools_mod !JT
      implicit real*8(a-h,o-z)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl_per/ iperiodic,ibasis
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB)

c Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

c get basis functions
      if(ndim.eq.3) then
        call basis_fns(iel,rvec_en,r_en)
       elseif(ndim.eq.2) then
         if(ibasis.eq.1) then
           call basis_fns_2d(iel,rvec_en,r_en)
         elseif(ibasis.eq.4) then
           call basis_fns_2dgauss(iel,rvec_en,r_en)
         elseif(ibasis.eq.5) then
           call basis_fns_polargauss(iel,rvec_en,r_en)
         elseif(ibasis.eq.6) then
           call basis_fns_2dgauss_noncirc(iel,rvec_en,r_en)
         else
           stop 'orbitals_loc_ana: ibasis must be 1,4,5 or 6 for 2d systems'
         endif
      endif

      do 25 iorb=1,norb
        do 25 ie=nelec1,nelec2
          orb(ie,iorb)=0
          dorb(1,ie,iorb)=0
          dorb(2,ie,iorb)=0
          dorb(3,ie,iorb)=0
          ddorb(ie,iorb)=0
          do 25 m=1,nbasis
      if(ipr.ge.5) write(6,'(''iorb,ie,m,iwf,coef(m,iorb,iwf),phin(m,ie)'',3i3,9g13.6)') iorb,ie,m,iwf,coef(m,iorb,iwf),phin(m,ie)
c     &,coef(m,iorb,iwf)*phin(m,ie)
            orb(ie,iorb)=orb(ie,iorb)+coef(m,iorb,iwf)*phin(m,ie)
            dorb(1,ie,iorb)=dorb(1,ie,iorb)+coef(m,iorb,iwf)*dphin(1,m,ie)
            dorb(2,ie,iorb)=dorb(2,ie,iorb)+coef(m,iorb,iwf)*dphin(2,m,ie)
            dorb(3,ie,iorb)=dorb(3,ie,iorb)+coef(m,iorb,iwf)*dphin(3,m,ie)
   25       ddorb(ie,iorb)=ddorb(ie,iorb)+coef(m,iorb,iwf)*d2phin(m,ie)

c      do  iorb=1,norb
c        do  ie=nelec1,nelec2
c          write(*,*) 'orb(ie,iorb),ie,iorb=',orb(ie,iorb),ie,iorb
c        enddo
c      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine orbitals_loc_anae(iel,rvec_en,r_en,orb)
c Written by Cyrus Umrigar
c Calculate localized orbitals for electron iel

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /dim/ ndim
c     common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,orb(MORB)

c get basis functions
      if(ndim.eq.3) then
        call basis_fnse2(iel,rvec_en,r_en)
       elseif(ndim.eq.2) then
        call basis_fns_2de2(iel,rvec_en,r_en)
      endif

      do 25 iorb=1,norb
          orb(iorb)=0
          do 25 m=1,nbasis
   25       orb(iorb)=orb(iorb)+coef(m,iorb,iwf)*phin(m,iel)
      return
      end

c---------------------------------------------------------------------------

      subroutine deriv_orbitals(rvec_en,r_en,orb,dorb,ddorb,dporb,d2porb
     &          ,ddporb,d2dporb)
c Written by A.D.Guclu (Apr 2005) starting from orbitals_loc_ana.f
c Calculate localized orbitals, coo. and parameter derivatives for all electrons


      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl_per/ iperiodic,ibasis
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC)
      common /deriv_phifun/ dparam(MOTYPE,MBASIS,MELEC)
     &,d2param(MOTYPE,MOTYPE,MBASIS,MELEC),ddparam(3,MOTYPE,MBASIS,MELEC)
     &,d2dparam(MOTYPE,MBASIS,MELEC)
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
      dimension orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB)
      dimension dporb(MOTYPE,MELEC,MORB),d2porb(MOTYPE,MOTYPE,MELEC,MORB)
      dimension ddporb(3,MOTYPE,MELEC,MORB),d2dporb(MOTYPE,MELEC,MORB)

      nelec1=1
      nelec2=nelec

c get basis functions
      if(ndim.ne.2) stop 'deriv_orbitals: ndim must be 2'
      if(ibasis.eq.4) then
        call deriv_2dgauss(rvec_en,r_en)
      elseif(ibasis.eq.5) then
        call deriv_polargauss(rvec_en,r_en)
      elseif(ibasis.eq.6) then
        call deriv_2dgauss_noncirc(rvec_en,r_en)
      else
        stop 'deriv_orbitals: ibasis must be 4, 5, or 6'
      endif

      do iorb=1,norb
        do ie=nelec1,nelec2

          orb(ie,iorb)=0
          do idim=1,ndim
            dorb(idim,ie,iorb)=0
          enddo
          ddorb(ie,iorb)=0
          do ip=1,notype
            dporb(ip,ie,iorb)=0
            d2dporb(ip,ie,iorb)=0
            do idim=1,ndim
              ddporb(idim,ip,ie,iorb)=0
            enddo
          enddo
          do ip=1,notype
            do jp=1,notype
              d2porb(ip,jp,ie,iorb)=0
            enddo
          enddo

          do m=1,nbasis

            orb(ie,iorb)=orb(ie,iorb)+coef(m,iorb,iwf)*phin(m,ie)
c            write(*,*) 'ie,iorb,m,coef(m,iorb,iwf),phin(m,ie)='
c     &,ie,iorb,m,coef(m,iorb,iwf),phin(m,ie)
            do idim=1,ndim
              dorb(idim,ie,iorb)=dorb(idim,ie,iorb)+coef(m,iorb,iwf)*dphin(idim,m,ie)
            enddo
            ddorb(ie,iorb)=ddorb(ie,iorb)+coef(m,iorb,iwf)*d2phin(m,ie)

            do ip=1,notype
              dporb(ip,ie,iorb)=dporb(ip,ie,iorb)+coef(m,iorb,iwf)*dparam(ip,m,ie)
              d2dporb(ip,ie,iorb)=d2dporb(ip,ie,iorb)+coef(m,iorb,iwf)*d2dparam(ip,m,ie)
              do idim=1,ndim
                ddporb(idim,ip,ie,iorb)=ddporb(idim,ip,ie,iorb)
     &           +coef(m,iorb,iwf)*ddparam(idim,ip,m,ie)
              enddo
            enddo
            do ip=1,notype
              do jp=1,notype
                d2porb(ip,jp,ie,iorb)=d2porb(ip,jp,ie,iorb)+coef(m,iorb,iwf)*d2param(ip,jp,m,ie)
              enddo
            enddo
          enddo
        enddo

c        write(*,*) 'iorb,orb(ie,iorb)=',(orb(ie,iorb),ie=1,nelec)
c        write(*,*) 'phin(m,ie)=',(phin(iorb,ie),ie=1,nelec)

      enddo

      return
      end

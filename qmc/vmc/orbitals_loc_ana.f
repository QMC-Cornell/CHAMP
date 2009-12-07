      subroutine orbitals_loc_ana(iel,rvec_en,r_en,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c Calculate localized orbitals and derivatives for all or 1 electrons
      use all_tools_mod
      use control_mod
      use coefs_mod
      use const_mod
      use dim_mod
      use wfsec_mod
      use contrl_per_mod
      use phifun_mod
      use atom_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent)
     &,orb(nelec,norb),dorb(3,nelec,norb),ddorb(nelec,norb)

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
c          write(6,*) 'orb(ie,iorb),ie,iorb=',orb(ie,iorb),ie,iorb
c        enddo
c      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine orbitals_loc_anae(iel,rvec_en,r_en,orb)
c Written by Cyrus Umrigar
c Calculate localized orbitals for electron iel
      use coefs_mod
      use dim_mod
      use wfsec_mod
      use phifun_mod
      use atom_mod
      use const_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent)
     &,orb(norb)

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
      use control_mod
      use coefs_mod
      use const_mod
      use dim_mod
      use wfsec_mod
      use contrl_per_mod
      use phifun_mod
      use optimo_mod
      use atom_mod
      use deriv_phifun_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent)
      dimension orb(nelec,norb),dorb(3,nelec,norb),ddorb(nelec,norb)
      dimension dporb(notype,nelec,norb),d2porb(notype,notype,nelec,norb)
      dimension ddporb(3,notype,nelec,norb),d2dporb(notype,nelec,norb)

      nelec1=1
      nelec2=nelec

!     JT: should move these allocations outside the subroutine for efficiency?
      call alloc ('dparam', dparam, notype, nbasis, nelec)
      call alloc ('d2param', d2param, notype, notype, nbasis, nelec)
      call alloc ('ddparam', ddparam, 3, notype, nbasis, nelec)
      call alloc ('d2dparam', d2dparam, notype, nbasis, nelec)

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
c            write(6,*) 'ie,iorb,m,coef(m,iorb,iwf),phin(m,ie)='
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

c        write(6,*) 'iorb,orb(ie,iorb)=',(orb(ie,iorb),ie=1,nelec)
c        write(6,*) 'phin(m,ie)=',(phin(iorb,ie),ie=1,nelec)

      enddo

      return
      end

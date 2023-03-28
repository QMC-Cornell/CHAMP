      subroutine corbitals_loc_ana_grade(iel,rvec_en,r_en,corb,cdorb,cddorb)
! orbitals_loc_ana_grade adapted to complex orbitals by A.D.Guclu Feb2004
! Calculate localized orbitals and derivatives for all or 1 electrons
      use all_tools_mod
      use coefs_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use wfsec_mod
      use atom_mod
      use cphifun_mod
      implicit real*8(a-h,o-z)

      complex(dpc) corb,cdorb,cddorb

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent) &
     &,corb(norb),cdorb(3,norb),cddorb(norb)

!     JT: should move these allocations outside the subroutine for efficiency?
      call alloc ('cphin', cphin, nbasis, nelec)
      call alloc ('cdphin', cdphin, 3, nbasis, nelec, nelec)
      call alloc ('cd2phin', cd2phin, nbasis, nelec)

! get basis functions
      if(numr.eq.0) then
        call cbasis_fns_fd(iel,rvec_en,r_en)
      else if(numr.eq.1) then
        call cbasis_fns_num(iel,rvec_en,r_en)
      endif

      do 25 iorb=1,norb
          corb(iorb)=dcmplx(0,0)
          do 10 idim=1,ndim
   10       cdorb(idim,iorb)=dcmplx(0,0)
          cddorb(iorb)=dcmplx(0,0)
          do 25 m=1,nbasis
           if(ipr.ge.2) write(6,*) 'cphin,cdphin1,cdphin2,cd2phin=' &
     &              ,cphin(m,iel),cdphin(1,m,iel,1),cdphin(2,m,iel,1),cd2phin(m,iel)
           corb(iorb)=corb(iorb)+coef(m,iorb,iwf)*cphin(m,iel)
           do 15 idim=1,ndim
   15        cdorb(idim,iorb)=cdorb(idim,iorb)+coef(m,iorb,iwf)*cdphin(idim,m,iel,1)
   25      cddorb(iorb)=cddorb(iorb)+coef(m,iorb,iwf)*cd2phin(m,iel)

      if(ipr.ge.4) then
        write(6,'(''corb='',100es12.4)') (corb(iorb),iorb=1,norb)
        write(6,'(''cdorb='',100es12.4)') ((cdorb(idim,iorb),idim=1,ndim),iorb=1,norb)
        write(6,'(''cddorb='',100es12.4)') (cddorb(iorb),iorb=1,norb)
      endif

      return
      end

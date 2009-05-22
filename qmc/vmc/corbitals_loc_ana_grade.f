      subroutine corbitals_loc_ana_grade(iel,rvec_en,r_en,corb,cdorb,cddorb)
c orbitals_loc_ana_grade adapted to complex orbitals by A.D.Guclu Feb2004
c Calculate localized orbitals and derivatives for all or 1 electrons

      use orbitals_mod, only: orb_tot_nb
      use coefs_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use wfsec_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'force.h'
!JT      include 'numbas.h'

      complex*16 corb,cdorb,cddorb
      complex*16 cphin,cdphin,cd2phin

!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /cphifun/ cphin(MBASIS,MELEC),cdphin(3,MBASIS,MELEC,MELEC)
     &,cd2phin(MBASIS,MELEC)
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
!JT      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
!JT     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
!JT     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
!JT      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,corb(orb_tot_nb),cdorb(3,orb_tot_nb),cddorb(orb_tot_nb)

c get basis functions
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
c           write(*,*) 'cphin,cdphin1,cdphin2,cd2phin='
c     &              ,cphin(m,iel),cdphin(1,m,iel),cdphin(2,m,iel),cd2phin(m,iel)
           corb(iorb)=corb(iorb)+coef(m,iorb,iwf)*cphin(m,iel)
           do 15 idim=1,ndim
   15        cdorb(idim,iorb)=cdorb(idim,iorb)+coef(m,iorb,iwf)*cdphin(idim,m,iel,1)
   25      cddorb(iorb)=cddorb(iorb)+coef(m,iorb,iwf)*cd2phin(m,iel)

      return
      end

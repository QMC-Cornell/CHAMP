      subroutine corbitals_loc_ana(iel,rvec_en,r_en,corb,cdorb,cddorb)
! orbitals_loc_ana adapted to complex orbitals by A.D.Guclu, Feb2004
! Calculate localized orbitals and derivatives for all or 1 electrons
      use all_tools_mod
      use coefs_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use wfsec_mod
      use cphifun_mod
      implicit real*8(a-h,o-z)

      complex(dpc) corb,cdorb,cddorb

      common /compferm/ emagv,nv,idot

      dimension rvec_en(3,nelec,*),r_en(nelec,*) &
     &,corb(nelec,norb),cdorb(3,nelec,nelec,norb),cddorb(nelec,norb)

!     JT: should move these allocations outside the subroutine for efficiency?
      call alloc ('cphin', cphin, nbasis, nelec)
      call alloc ('cdphin', cdphin, 3, nbasis, nelec, nelec)
      call alloc ('cd2phin', cd2phin, nbasis, nelec)

! Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

! get basis functions
! nej controls whether if we have correlated basis set
      if(numr.eq.1) then
        call cbasis_fns_num(iel,rvec_en,r_en)
        nej=1
      elseif(idot.eq.3) then
        if(iel.ne.0) stop '1 electron move not possible with projected comp.ferm.'
        call cbasis_fns_cf(rvec_en)
        nej=nelec
      else
        call cbasis_fns_fd(iel,rvec_en,r_en)
        nej=1
      endif

!      do 5 ib=1,nbasis
!        write(6,'(''ib,cphin='',i3,(30f9.5))') ib,(cphin(i,ib),i=1,nelec)
!        write(6,'(''ib,cdphin1='',i3,(30f9.5))') ib,(cdphin(1,i,ib),i=1,nelec)
!        write(6,'(''ib,cdphin2='',i3,(30f9.5))') ib,(cdphin(1,i,ib),i=1,nelec)
! 5      write(6,'(''ib,cd2phin='',i3,(30f9.5))') ib,(cd2phin(i,ib),i=1,nelec)

      do 25 iorb=1,norb
        do 25 ie=nelec1,nelec2
          corb(ie,iorb)=dcmplx(0,0)
            do 10 idim=1,ndim
              do 10 je=1,nej
   10         cdorb(idim,ie,je,iorb)=dcmplx(0,0)
          cddorb(ie,iorb)=dcmplx(0,0)

          do 25 m=1,nbasis
            corb(ie,iorb)=corb(ie,iorb)+coef(m,iorb,iwf)*cphin(m,ie)
            do 15 idim=1,ndim
              do 15 je=1,nej
   15           cdorb(idim,ie,je,iorb)=cdorb(idim,ie,je,iorb)+coef(m,iorb,iwf)*cdphin(idim,m,ie,je)
   25       cddorb(ie,iorb)=cddorb(ie,iorb)+coef(m,iorb,iwf)*cd2phin(m,ie)

      return
      end
!-----------------------------------------------------------------------


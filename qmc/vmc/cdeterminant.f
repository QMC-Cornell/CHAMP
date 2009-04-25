      subroutine cdeterminant(x,rvec_en,r_en,cddet_det,d2lndet,div_vd,determ)
c same subroutine as determinant() adapted to complex orbitals/determinants
c by A.D.Guclu Feb2004
c can deal with any complex determinant, provided that the determinantal
c coefficients are real.

      use control_mod
      use basic_tools_mod
      use cslater_mod

      implicit real*8(a-h,o-z)
!      include 'vmc.h'
!      include 'force.h'
!      include 'fit.h'

!      parameter (one=1.d0,half=0.5d0)
      character*16 mode

c routine to calculate the value, gradient and Laplacian of the
c determinantal part of the wavefunction
c Also derivatives wrt csf_coefs for optimizing them.

c cdeterm =  complex value of the determinants
c determ  =  magnitude of cdeterm
c cddet_det   =  d/dx(ln(det)))
c d2lndet  =  real(d2/dx2(det)/det) - abs(cddet_det)**2

c Note that the first dimension of the slater matrices is MMAT_DIM = (MELEC/2)**2.
c The first dimension of the Slater matrices must be at least max(nup**2,ndn**2)
c So, we check in read_input that nup and ndn are each <= MELEC/2.

c complex argument:
      complex*16 cddet_det(3,*)

c complex local:
      complex*16 cd2lndet,cdeterm,cterm,cdetinv
      complex*16 corb(MELEC,MORB),cdorb(3,MELEC,MELEC,MORB),cddorb(MELEC,MORB)
      complex*16 cekinen(MELEC)
      complex*16 cauxx(MELEC,MELEC)

c following complex are now defined in module cslater_mod
c     complex*16 cslmui,cslmdi,cfpu,cfpd,cfppu,cfppd,cdetu,cdetd,cddeti_deti,cd2edeti_deti
c     complex*16 cdeti_det,cddeti_det,cd2deti_det,cd2det_det

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode
      common /contrl_opt2/ igradhess,iadd_diag_opt
      common /contrl_per/ iperiodic,ibasis
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /kinet/ ekineo(MELEC),ekinen(MELEC)
      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
c      common /cslater/ cslmui(MMAT_DIM,MDET),cslmdi(MMAT_DIM,MDET)
c     &,cfpu(3,MMAT_DIM,MDET),cfpd(3,MMAT_DIM,MDET)
c     &,cfppu(MMAT_DIM,MDET),cfppd(MMAT_DIM,MDET)
c     &,cdetu(MDET),cdetd(MDET)
c     &,cddeti_deti(3,MELEC,MDET),cd2edeti_deti(MELEC,MDET)
c     &,cdeti_det(MCSF),cddeti_det(3,MELEC,MCSF),cd2deti_det(MCSF),cd2det_det

      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /dojasderiv/ ijasderiv
      common /optim/ lo(MORB),npoint(MORB),
     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
     &imnbas(MCENT),
     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
     &necn,nebase

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),div_vd(MELEC)

c allocate memory (maybe better to allocate in read_input ?) :
      n2=nelec*nelec
      call alloc('cslmui',cslmui,n2,ndetup)
      call alloc('cslmdi',cslmdi,n2,ndetdn)
      call alloc('cfpu',cfpu,ndim,n2,ndetup)
      call alloc('cfpd',cfpd,ndim,n2,ndetdn)
      call alloc('cfppu',cfppu,n2,ndetup)
      call alloc('cfppd',cfppd,n2,ndetdn)
      call alloc('cdetu',cdetu,ndetup)
      call alloc('cdetd',cdetd,ndetdn)
      call alloc('cddeti_deti',cddeti_deti,ndim,nelec,ndet)
      call alloc('cd2edeti_deti',cd2edeti_deti,nelec,ndet)

c initialize the derivative arrays to zero
      do 10 i=1,nelec
        cekinen(i)=dcmplx(0,0)
        do 10 idim=1,ndim
   10     cddet_det(idim,i)=dcmplx(0,0)
      cd2lndet=dcmplx(0,0)

c initialize the determinant arrays to one
      do 20 idet=1,ndetup
   20   cdetu(idet)=dcmplx(one,0)
      do 22 idet=1,ndetdn
   22   cdetd(idet)=dcmplx(one,0)
c      cdeterm=dcmplx(0,0)

c get orbitals and derivatives for all electrons
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          call corbitals_loc_ana(0,rvec_en,r_en,corb,cdorb,cddorb)
         else
          stop 'complex calculations of num. orbs not implemented yet'
        endif

      else
          stop 'complex calculations of periodic system not implemented'
      endif

      if(ipr.ge.4) then
        do 26 iorb=1,norb
          write(6,'(''iorb,corb='',i3,(30f9.5))') iorb,(corb(i,iorb),i=1,nelec)
          write(6,'(''iorb,cdorb1='',i3,(30f9.5))') iorb,(cdorb(1,i,1,iorb),i=1,nelec)
          write(6,'(''iorb,cdorb2='',i3,(30f9.5))') iorb,(cdorb(2,i,1,iorb),i=1,nelec)
   26     write(6,'(''iorb,cddorb='',i3,(30f9.5))') iorb,(cddorb(i,iorb),i=1,nelec)
      endif


c The 3 nested loops are over
c 1) up electrons,
c 2) determinants
c 3) basis states setting up transpose of the Slater
c matrix in slmui to get inverse transpose.
c Also put derivatives in fpu and fppu.
      ik=-nup
      do 30 i=1,nup
        ik=ik+nup
        do 30 idet=1,ndetup
          jk=-nup
          do 30 j=1,nup
            jk=jk+nup
            cslmui(i+jk,idet)=corb(i,iworbdup(j,idet))
            do 28 idim=1,ndim
   28         cfpu(idim,j+ik,idet)=cdorb(idim,i,1,iworbdup(j,idet))
   30       cfppu(j+ik,idet)=cddorb(i,iworbdup(j,idet))

c loop through number of determinants calculating the inverse
c transpose matrices and their determinants
c (cmatinv from Wolfgang's subroutine. is the matrix transfer really necessary?)
      do 40 idet=1,ndetup
        do in=0,nup*nup-1
            i=in/nup+1
            j=mod(in,nup)+1
            cauxx(i,j)=cslmui(in+1,idet)
        enddo
        call cmatinv(cauxx,nup,cdetu(idet))
        do in=0,nup*nup-1
            i=in/nup+1
            j=mod(in,nup)+1
            cslmui(in+1,idet)=cauxx(i,j)
        enddo
   40 continue

c repeat above for down spins
      ik=-ndn
      do 50 i=1,ndn
        ik=ik+ndn
        do 50 idet=1,ndetdn
          jk=-ndn
          do 50 j=1,ndn
            jk=jk+ndn
            cslmdi(i+jk,idet)=corb(i+nup,iworbddn(j,idet))
            do 45 idim=1,ndim
   45         cfpd(idim,j+ik,idet)=cdorb(idim,i+nup,1,iworbddn(j,idet))
   50       cfppd(j+ik,idet)=cddorb(i+nup,iworbddn(j,idet))

      do 60 idet=1,ndetdn
        do in=0,ndn*ndn-1
            i=in/ndn+1
            j=mod(in,ndn)+1
            cauxx(i,j)=cslmdi(in+1,idet)
        enddo
        call cmatinv(cauxx,ndn,cdetd(idet))
        do in=0,ndn*ndn-1
            i=in/ndn+1
            j=mod(in,ndn)+1
            cslmdi(in+1,idet)=cauxx(i,j)
        enddo
   60 continue

c      if(ipr.ge.4) write(6,'(''cdetu,cdetd,cdet'',9d12.5)') cdetu(1),cdetd(1),cdet(1,1)

c set up sum of slater determinants along with their
c coefficients that were included in input data
      do 70 idet=1,max(ndetup,ndetdn)
c        cdeterm=cdeterm+cdetu(idet)*cdetd(idet)*cdet(idet,iwf)

c zero out temporary derivative arrays
        do 70 i=1,nelec
          do 65 idim=1,ndim
   65       cddeti_deti(idim,i,idet)=dcmplx(0,0)
   70     cd2edeti_deti(i,idet)=dcmplx(0,0)

c loop through up spin electrons
c take inner product of transpose inverse with derivative
c vectors to get (1/detup)*d(detup)/dx and (1/detup)*d2(detup)/dx**2

      do 80 idet=1,ndetup
        ik=-nup
        do 80 i=1,nup
          ik=ik+nup
          do 80 j=1,nup
            do 75 idim=1,ndim
   75         cddeti_deti(idim,i,idet)=cddeti_deti(idim,i,idet)+cslmui(j+ik,idet)*cfpu(idim,j+ik,idet)
   80       cd2edeti_deti(i,idet)=cd2edeti_deti(i,idet)+cslmui(j+ik,idet)*cfppu(j+ik,idet)

c repeat above for down spins
      do 90 idet=1,ndetdn
        ik=-ndn
        do 90 i=nup+1,nelec
          ik=ik+ndn
          do 90 j=1,ndn
            do 85 idim=1,ndim
   85         cddeti_deti(idim,i,idet)=cddeti_deti(idim,i,idet)+cslmdi(j+ik,idet)*cfpd(idim,j+ik,idet)
   90         cd2edeti_deti(i,idet)=cd2edeti_deti(i,idet)+cslmdi(j+ik,idet)*cfppd(j+ik,idet)

c combine results for up and down spins to get d(det)/dx
c and d2(det)/dx in ddet_det and d2lndet respectively

      cdeterm=0
      do 117 icsf=1,ncsf
        do 117 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
c         cterm=cdetu(idet)*cdetd(idet)*cdet(idet,iwf)
          if(ndn.ge.1) then
            cterm=cdetu(iwdetup(idet))*cdetd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
           else
            cterm=cdetu(iwdetup(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
          endif
          cdeterm=cdeterm+cterm
          do 115 i=1,nup
            iwdet=iwdetup(idet)
            do 113 idim=1,ndim
  113         cddet_det(idim,i)=cddet_det(idim,i)+cddeti_deti(idim,i,iwdet)*cterm
            cekinen(i)=cekinen(i)+cd2edeti_deti(i,iwdet)*cterm
  115       cd2lndet=cd2lndet+cd2edeti_deti(i,iwdet)*cterm
          do 117 i=nup+1,nelec
            iwdet=iwdetdn(idet)
            do 116 idim=1,ndim
  116         cddet_det(idim,i)=cddet_det(idim,i)+cddeti_deti(idim,i,iwdet)*cterm
            cekinen(i)=cekinen(i)+cd2edeti_deti(i,iwdet)*cterm
  117       cd2lndet=cd2lndet+cd2edeti_deti(i,iwdet)*cterm

c form inverse of full sum of determinants
      cdetinv=one/cdeterm
      determ=cdabs(cdeterm)

c multiply through to set up logarithmic first and second derivatives.

c      cd2lndet=cd2lndet*cdetinv
c d2lndet=dreal(cd2lndet), but ignore the imaginary part
      d2lndet=dreal(cd2lndet)*dreal(cdetinv)-dimag(cd2lndet)*dimag(cdetinv)
      do 120 i=1,nelec
c div_vd(i)=dreal(cdiv_vd(i))
        div_vd(i)=dreal(cekinen(i))*dreal(cdetinv)-dimag(cekinen(i))*dimag(cdetinv)
        ekinen(i)=-half*div_vd(i)
        do 120 k=1,ndim
          cddet_det(k,i)=cddet_det(k,i)*cdetinv
          temp=cdabs(cddet_det(k,i))**2
          div_vd(i)=div_vd(i)-temp
  120     d2lndet=d2lndet-temp

c Derivatives wrt to csf_coefs for optimizing them
c Note that the arrays that are needed for vmc and dmc are over ndet but
c those that are needed for optimization only are over nparmd.
      if(index(mode,'fit').ne.0 .or. igradhess.gt.0) then
        d2det_det=0
        do 125 i=1,nelec
  125     d2det_det=d2det_det-2*ekinen(i)
        do 140 iparm=1,nparmcsf
          stop 'optim. of cfs coeffs not possible for complex wfs yet'
  140   continue
      endif

      return
      end

c----------------------------------------------------------------------------------

      subroutine cdeterminant_cf(x,rvec_en,r_en,cddet_det,d2lndet,div_vd,determ)
c same subroutine as cdeterminant() adapted to composite fermions
c laplacian (LLL wfs are analytical) is zero and ignored here.
c by A.D.Guclu sep2004
c can deal only with spin polarized systems for the moment

      use control_mod
      use basic_tools_mod
      use cslater_cf_mod

      implicit real*8(a-h,o-z)
!      include 'vmc.h'
!      include 'force.h'
!      include 'fit.h'

!      parameter (one=1.d0,half=0.5d0)
      character*16 mode

c cdeterm =  complex value of the determinants
c determ  =  magnitude of cdeterm
c cddet_det   =  d/dx(ln(det)))
c d2lndet  =  real(d2/dx2(det)/det) - abs(cddet_det)**2

c note that the dimension of the slater matrices is assumed
c to be given by MMAT_DIM = (MELEC/2)**2, that is there are
c as many ups as downs. If this is not true then be careful if
c nelec is close to MELEC. The Slater matrices must be
c dimensioned at least max(nup**2,ndn**2)

c complex argument:
      complex*16 cddet_det(3,*)

c complex local:
      complex*16 cd2lndet,cdeterm
      complex*16 corb(MELEC,MORB),cdorb(3,MELEC,MELEC,MORB),cddorb(MELEC,MORB)
      complex*16 cekinen(MELEC)
      complex*16 cterm,cdetinv
      complex*16 cauxx(MELEC,MELEC)

c following are defined now defined in cslater_cf_mod for dynamic allo
c      complex*16 cslmui(MMAT_DIM,MDET),cslmdi(MMAT_DIM,MDET)
c     &,cfpu(3,MMAT_DIM,MELEC,MDET),cfpd(3,MMAT_DIM,MELEC,MDET)
c     &,cdetu(MDET),cdetd(MDET)
c     &,cddeti_deti(3,MELEC,MDET),cd2edeti_deti(MELEC,MDET)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode
      common /contrl_opt2/ igradhess,iadd_diag_opt
      common /contrl_per/ iperiodic,ibasis
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /kinet/ ekineo(MELEC),ekinen(MELEC)
      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /dojasderiv/ ijasderiv
      common /optim/ lo(MORB),npoint(MORB),
     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
     &imnbas(MCENT),
     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
     &necn,nebase

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),div_vd(MELEC)

c allocate memory:
      n2=nelec*nelec
      call alloc('cslmui',cslmui,n2,ndetup)
      call alloc('cslmdi',cslmdi,n2,ndetdn)
      call alloc('cfpu',cfpu,ndim,n2,nelec,ndetup)
      call alloc('cfpd',cfpd,ndim,n2,nelec,ndetdn)
      call alloc('cdetu',cdetu,ndetup)
      call alloc('cdetd',cdetd,ndetdn)
      call alloc('cddeti_deti',cddeti_deti,ndim,nelec,ndet)
      call alloc('cd2edeti_deti',cd2edeti_deti,nelec,ndet)

c initialize the derivative arrays to zero
      do 10 i=1,nelec
        cekinen(i)=dcmplx(0,0)   ! not used
        do 10 idim=1,ndim
   10     cddet_det(idim,i)=dcmplx(0,0)
      cd2lndet=dcmplx(0,0)       ! not used

c initialize the determinant arrays to one
      do 20 idet=1,ndetup
   20   cdetu(idet)=dcmplx(one,0)
      do 22 idet=1,ndetdn
   22   cdetd(idet)=dcmplx(one,0)
      cdeterm=dcmplx(0,0)

c      write(*,*) 'flag1'

c      if(nelec.ne.nup) stop 'system must be fully polarized for projected composite fermions!'

c get orbitals and derivatives for all electrons
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          call corbitals_loc_ana(0,rvec_en,r_en,corb,cdorb,cddorb)
         else
          stop 'complex calculations of num. orbs not implemented yet'
        endif

      else
          stop 'complex calculations of periodic system not implemented'
      endif

c      if(ipr.ge.4) then
c        do 26 iorb=1,norb
c          write(6,'(''iorb,corb='',i3,(30f9.5))') iorb,(corb(i,iorb),i=1,nelec)
c          write(6,'(''iorb,cdorb1='',i3,(30f9.5))') iorb,(cdorb(1,i,iorb),i=1,nelec)
c          write(6,'(''iorb,cdorb2='',i3,(30f9.5))') iorb,(cdorb(2,i,iorb),i=1,nelec)
c   26     write(6,'(''iorb,cddorb='',i3,(30f9.5))') iorb,(cddorb(i,iorb),i=1,nelec)
c      endif

c The 3 nested loops are over
c 1) up electrons,
c 2) determinants
c 3) basis states setting up transpose of the Slater
c matrix in slmui to get inverse transpose.

      ik=-nup
      do 30 i=1,nup
        ik=ik+nup
        do 30 idet=1,ndetup
          jk=-nup
          do 30 j=1,nup
            jk=jk+nup
            cslmui(i+jk,idet)=corb(i,iworbdup(j,idet))  ! transpose
            do 30 idim=1,ndim
c WARNING, should the following loop go up to nelec for partial
c polarization? need to think about it. composite fermions have never been
c tested for partially polarized systems.
              do 30 k=1,nup
                cfpu(idim,j+ik,k,idet)=cdorb(idim,i,k,iworbdup(j,idet))
   30 continue

c loop through number of determinants calculating the inverse
c transpose matrices and their determinants
c (cmatinv from Wolfgang's subroutine. is the matrix transfer really necessary?)
      do 40 idet=1,ndetup
        do in=0,nup*nup-1
            i=in/nup+1
            j=mod(in,nup)+1
            cauxx(i,j)=cslmui(in+1,idet)
        enddo
        call cmatinv(cauxx,nup,cdetu(idet))
c        write(*,*) 'det=',cdetu(idet)
        do in=0,nup*nup-1
            i=in/nup+1
            j=mod(in,nup)+1
            cslmui(in+1,idet)=cauxx(i,j)
        enddo
   40   continue

c repeat above for down spins
      ik=-ndn
      do 50 i=1,ndn
        stop 'down electron not tested for composite fermions'
        ik=ik+ndn
        do 50 idet=1,ndetdn
          jk=-ndn
          do 50 j=1,ndn
            jk=jk+ndn
            cslmdi(i+jk,idet)=corb(i+nup,iworbddn(j,idet))
            do 45 idim=1,ndim
              do 45 k=1,ndn
   45           cfpd(idim,j+ik,k,idet)=cdorb(idim,i+nup,k,iworbddn(j,idet))
   50       continue

      do 60 idet=1,ndetdn
        do in=0,ndn*ndn-1
            i=in/ndn+1
            j=mod(in,ndn)+1
            cauxx(i,j)=cslmdi(in+1,idet)
        enddo
        call cmatinv(cauxx,ndn,cdetd(idet))
        do in=0,ndn*ndn-1
            i=in/ndn+1
            j=mod(in,ndn)+1
            cslmdi(in+1,idet)=cauxx(i,j)
        enddo
   60 continue

c      write(*,*) 'flag3'

c      if(ipr.ge.4) write(6,'(''cdetu,cdetd,cdet'',9d12.5)') cdetu(1),cdetd(1),cdet(1,1)

c set up sum of slater determinants along with their
c coefficients that were included in input data
      do 70 idet=1,max(ndetup,ndetdn)
c        cdeterm=cdeterm+cdetu(idet)*cdetd(idet)*cdet(idet,iwf)

c zero out temporary derivative arrays
        do 70 i=1,nelec
          do 65 idim=1,ndim
   65       cddeti_deti(idim,i,idet)=dcmplx(0,0)
   70     cd2edeti_deti(i,idet)=dcmplx(0,0)   !  not used

c loop through up spin electrons
c take inner product of transpose inverse with derivative
c vectors to get (1/detup)*d(detup)/dx and (1/detup)*d2(detup)/dx**2

      do 75 idet=1,ndetup
        do 75 i=1,nup
          kk=-nup
          do 75 k=1,nup
            kk=kk+nup
            do 75 j=1,nup
              do 75 idim=1,ndim
   75           cddeti_deti(idim,i,idet)=cddeti_deti(idim,i,idet)+cslmui(j+kk,idet)*cfpu(idim,j+kk,i,idet)

c repeat above for down spins
      do 85 idet=1,ndetdn
        do 85 i=1,ndn
          kk=-ndn
          do 85 k=1,ndn
            kk=kk+ndn
            do 85 j=1,ndn
              do 85 idim=1,ndim
   85           cddeti_deti(idim,i+nup,idet)=cddeti_deti(idim,i+nup,idet)+cslmdi(j+kk,idet)*cfpd(idim,j+kk,i,idet)


c        cterm=cdetu(idet)*cdetd(idet)*cdet(idet,iwf)
c        do 95 i=1,nelec
c          do 95 idim=1,ndim
c   95       cddet_det(idim,i)=cddet_det(idim,i)+cddeti_deti(idim,i,idet)*cterm
  110 continue

c      write(*,*) 'flag4'

      cdeterm=0
      do 115 icsf=1,ncsf
        do 115 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
c         cterm=cdetu(idet)*cdetd(idet)*cdet(idet,iwf)
          if(ndn.ge.1) then
            cterm=cdetu(iwdetup(idet))*cdetd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
           else
            cterm=cdetu(iwdetup(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
          endif
          cdeterm=cdeterm+cterm
          do 115 i=1,nelec
            do 115 idim=1,ndim
  115         cddet_det(idim,i)=cddet_det(idim,i)+cddeti_deti(idim,i,idet)*cterm

c form inverse of full sum of determinants
      cdetinv=one/cdeterm
      determ=cdabs(cdeterm)

      d2lndet=0.d0
      do 120 i=1,nelec
c div_vd(i)=dreal(cdiv_vd(i))
        div_vd(i)=0.d0
        ekinen(i)=0.d0
        do 120 k=1,ndim
          cddet_det(k,i)=cddet_det(k,i)*cdetinv
          temp=cdabs(cddet_det(k,i))**2
          div_vd(i)=div_vd(i)-temp
  120     d2lndet=d2lndet-temp

c      write(*,*) 'flag5'

c Derivatives wrt to csf_coefs for optimizing them
c Note that the arrays that are needed for vmc and dmc are over ndet but
c those that are needed for optimization only are over nparmd.
      if(index(mode,'fit').ne.0 .or. igradhess.gt.0) then
        d2det_det=0
        do 125 i=1,nelec
  125     d2det_det=d2det_det-2*ekinen(i)
        do 140 iparm=1,nparmcsf
          stop 'optimization of cfs coefficients not possible for complex wfs yet'

c          icsf=iwcsf(iparm)
c          d2deti_det(iparm)=0
c          deti_det(iparm)=0
c          do 130 i=1,nelec
c            do 130 k=1,ndim
c  130         ddeti_det(k,i,iparm)=0
c          do 140 idet_in_csf=1,ndet_in_csf(icsf)
c            idet=iwdet_in_csf(idet_in_csf,icsf)
c            term=detu(idet)*detd(idet)*cdet_in_csf(idet_in_csf,icsf)*detinv
c            deti_det(iparm)=deti_det(iparm)+term
c            do 140 i=1,nelec
c              d2deti_det(iparm)=d2deti_det(iparm)+d2edeti_deti(i,idet)*term
c              do 140 k=1,ndim
c  140           ddeti_det(k,i,iparm)=ddeti_det(k,i,iparm)+ddeti_deti(k,i,idet)*term

  140   continue

c      if(ipr.ge.4) write(6,'(''deti_det(iparm) in determinant'',40d12.4)') (deti_det(iparm),iparm=1,nparmcsf)

      endif

c      write(*,*) 'flag6'

c release memory. we do not need to keep these because no vmc_mov1 for cfermions.
c on the other hand, releasing not very meaningful since the rest of the code
c do not use memory allocation? keep it for now.
      call release('cslmui',cslmui)
      call release('cslmdi',cslmdi)
      call release('cfpu',cfpu)
      call release('cfpd',cfpd)
      call release('cdetu',cdetu)
      call release('cdetd',cdetd)
      call release('cddeti_deti',cddeti_deti)
      call release('cd2edeti_deti',cd2edeti_deti)

c      write(*,*) 'flag7'

      return
      end

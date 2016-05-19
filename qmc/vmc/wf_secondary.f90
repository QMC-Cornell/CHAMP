      subroutine wf_secondary
! Written by Claudia Filippi and Cyrus Umrigar
! Read parameters for secondary wavefns.
      use constants_mod
      use control_mod
      use orbitals_mod
      use mpi_mod
      use atom_mod
      use dorb_mod
      use coefs_mod
      use dets_mod
      use basis1_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use contr2_mod
      use wfsec_mod
      use contrl_per_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use bparm_mod
      use jaspar1_mod
      use jaspar2_mod
      implicit real*8(a-h,o-z)

      character*20 filename,wfile,fmt
      character*30 section

      common /wfname/ filename

      do 400   iwt=2,nwftype

      if(iwt.lt.10) then
        write(wfile,'(i1)') iwt
       elseif(iwt.lt.100) then
        write(wfile,'(i2)') iwt
       else
        stop 'wf_secondary, increase MFORCE'
      endif

      wfile=filename(1:index(filename,' ')-1)//'.'//wfile(1:index(wfile,' ')-1)
      open(3,file=wfile,status='old')

! Determinantal section
      read(3,*) section
      write(6,'(/,a30,''(secondary)'',/)') section

      read(3,*)
      read(3,*)

      do 10 i=1,nctype
        read(3,*)
        read(3,*)
   10 continue

      if(ibasis.eq.1 .or. ibasis.eq.3) then
        write(6,'(''secondary orbital coefficients'')')
        do 20 i=1,orb_tot_nb
          read(3,*) (coef(j,i,iwt),j=1,nbasis)
   20     write(6,'(12f10.6)') (coef(j,i,iwt),j=1,nbasis)
        read(3,*) (zex(i,iwt),i=1,nbasis)
        write(6,'(''secondary screening constants'')')
        write(6,'(12f10.6)') (zex(i,iwt),i=1,nbasis)
      endif

!     read(3,*) (cdet(i,iwt),i=1,ndet)
!     write(6,'(/,''determinant coefficients'')')
!     write(6,'(20f10.6)') (cdet(idet,iwt),idet=1,ndet)

!     do 30 i=1,ndet
!       read(3,*) (iworbd(j,i),j=1,nelec)
!       write(6,*) (iworbd(j,i),j=1,nelec)
!       do 25 j=1,nelec
!         if(iworbd(j,i).gt.norb) then
!           write(6,'(''j,i,iwt,iworbd(j,i),norb='',9i5)') j,i,iwt,iworbd(j,i),norb
!           stop 'iworbd(j,i) > norb'
!         endif
!  25   continue
!       write(fmt,'(''('',i3,''i4,3x,'',i3,''i4)'')') nup,ndn
!  30   write(6,fmt) (iworbd(j,i),j=1,nup),(iworbd(j+nup,i),j=1,ndn)

      write(6,'(/,''secondary orbitals in determinants'')')
      do 72 i=1,ndet
        read(3,*) (iworbd(j,i),j=1,nelec)
        do 70 j=1,nelec
          if(iworbd(j,i).gt.orb_tot_nb) then
            write(6,'(''j,i,iwt,iworbd(j,i),orb_tot_nb='',9i5)') j,i,iwt,iworbd(j,i),orb_tot_nb
            stop 'iworbd(j,i) > orb_tot_nb'
          endif
   70     continue
        if(nup+ndn.lt.60) then
          write(fmt,'(''('',i2,''i3,3x,'',i2,''i3)'')') nup,ndn
          write(6,fmt) (iworbd(j,i),j=1,nup),(iworbd(j+nup,i),j=1,ndn)
         else
          write(6,'(30i4)') (iworbd(j,i),j=1,nup)
          write(6,'(30i4)') (iworbd(j+nup,i),j=1,ndn)
        endif
        do 71 j=2,nup
          do 71 jj=1,j-1
   71       if(iworbd(jj,i).eq.iworbd(j,i)) stop 'An up-spin determinant has 2 identical orbitals'
        do 72 j=2,ndn
          do 72 jj=1,j-1
   72       if(iworbd(jj+nup,i).eq.iworbd(j+nup,i)) stop 'A down-spin determinant has 2 identical orbitals'

      read(3,*)
      read(3,*) (csf_coef(icsf,iwt),icsf=1,ncsf)
      write(6,'(/,''secondary CSF coefs='',20f10.6)') (csf_coef(icsf,iwt),icsf=1,ncsf)
!     read(3,*) (ndet_in_csf(icsf),icsf=1,ncsf)
      read(3,*)
!     write(6,'(''ndet_in_csf='',20i4)') (ndet_in_csf(icsf),icsf=1,ncsf)
!     do 75 idet=1,ndet
!  75   iflag(idet)=0
      do 85 icsf=1,ncsf
!       read(3,*) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
        read(3,*)
!       write(6,'(''CSF'',i4,'' iwdet_in_csf='',10i4)') icsf,(iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
!       read(3,*) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
        read(3,*)
!       write(6,'(''CSF'',i4,'' cdet_in_csf='',10f8.5)') icsf,(cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
!       do 85 idet_in_csf=1,ndet_in_csf(icsf)
!         iflag(iwdet_in_csf(idet_in_csf,icsf))=1
!  85     if(iwdet_in_csf(idet_in_csf,icsf).gt.ndet) stop 'iwdet_in_csf(idet_in_csf,icsf) > ndet'
   85 continue
!cCheck if all the determinants are used in CSFs
!     do 90 idet=1,ndet
!  90   if(iflag(idet).eq.0) write(6,'(''Warning: determinant'',i3,'' is unused'')') idet

!     if((ibasis.eq.1.or.ibasis.eq.3).and.numr.gt.0) call read_bas_num(iwt)
!     if((ibasis.eq.1.or.ibasis.eq.3).and.minval(zex(:,1)).eq.0.d0) call read_bas_num(iwt)
      if(ibasis.eq.1.or.ibasis.eq.3) call read_bas_num(iwt)

! Jastrow section
      read(3,*) section
      write(6,'(/,a30,''(secondary)'',/)') section

      read(3,*)
      read(3,*)

      write(6,'(''secondary ijas,isc='',9i2)') ijas,isc

      if(ijas.eq.1) then
        read(3,*) cjas1(iwt),cjas2(iwt)
        write(6,'(''jastrow numerator,denominator ='',2f10.5)') &
     &  cjas1(iwt),cjas2(iwt)
       elseif(ijas.eq.2) then
        nparm_read=69
        if(isc.ge.2) read(3,*) scalek(iwt),a21
        write(6,'(''scalek,a21='',t31,9f10.5)') scalek(iwt),a21
        do 270 isp=nspin1,nspin2
          read(3,*) (a1(iparm,isp,iwt),iparm=1,nparm_read)
          if(ncent.gt.1.and.a1(2,isp,iwt).ne.zero) &
     &    write(6,'(''WARNING e-n cusp condition cannot be imposed'', &
     &    '' for molecules'',/,''with present weighted form of'', &
     &    '' Jastrow'')')
          write(6,'(''a='',x,7f10.6,(8f10.6))') &
     &                 (a1(iparm,isp,iwt),iparm=1,nparm_read)
  270   continue
        do 275 isp=nspin1,nspin2
          read(3,*) (a2(iparm,isp,iwt),iparm=1,nparm_read)
  275     write(6,'(''b='',x,7f10.6,(8f10.6))') &
     &                 (a2(iparm,isp,iwt),iparm=1,nparm_read)
       elseif(ijas.eq.3) then
        nparm_read=2
        nparmc_read=(nord**3+5*nord)/6+nord**2+nord
        if(isc.ge.2) then
          read(3,*) scalek(iwt),a21
          write(6,'(''scalek(iwt),a21='',2f10.5)') scalek(iwt),a21
        endif
        read(3,*) (a(iparm,iwt),iparm=1,nparm_read)
        write(6,'(''a='',x,7f10.6,(8f10.6))')(a(iparm,iwt),iparm=1,nparm_read)
        do 280 isp=nspin1,nspin2b
          read(3,*) (b(iparm,isp,iwt),iparm=1,nparm_read)
  280     write(6,'(''b='',x,7f10.6,(8f10.6))') &
     &                (b(iparm,isp,iwt),iparm=1,nparm_read)
        do 290 it=1,nctype
          read(3,*) (c(iparm,it,iwt),iparm=1,nparmc_read)
  290     write(6,'(''c='',x,7f10.6,(8f10.6))') (c(iparm,it,iwt), &
     &    iparm=1,nparmc_read)
        if(ifock.gt.0) then
          nfock=9
          if(ifock.eq.2) nfock=15
          do 300 it=1,nctype
            read(3,*) (fck(iparm,it,iwt),iparm=1,nfock)
            if(ifock.gt.2) then
              call scale3(iwt,it)
            endif
            write(6,'(''f='',x,7f10.6,(8f10.6))') (fck(iparm,it,iwt), &
     &      iparm=1,nfock)
  300     continue
        endif
       elseif(ijas.ge.4.and.ijas.le.6) then
        if(ifock.gt.0) stop 'fock not yet implemented for ijas=4,5,6'
        read(3,*)
        nparma_read=2+max(0,norda-1)
        nparmb_read=2+max(0,nordb-1)
        nparmc_read=nterms4(nordc)

        if(isc.ge.2) then
          read(3,*) scalek(iwt),a21
          write(6,'(''scalek(iwt),a21='',2f10.5)') scalek(iwt),a21
        endif
        do 301 it=1,nctype
           read(3,*) (a4(iparm,it,iwt),iparm=1,nparma_read)
  301      write(6,'(''a='',x,7f10.6,(8f10.6))') &
     &                     (a4(iparm,it,iwt),iparm=1,nparma_read)
        do 302 isp=nspin1,nspin2b
          read(3,*) (b(iparm,isp,iwt),iparm=1,nparmb_read)
  302     write(6,'(''b='',x,7f10.6,(8f10.6))') &
     &                 (b(iparm,isp,iwt),iparm=1,nparmb_read)
        do 303 it=1,nctype
          read(3,*) (c(iparm,it,iwt),iparm=1,nparmc_read)
  303     write(6,'(''c='',x,7f10.6,(8f10.6))') (c(iparm,it,iwt), &
     &    iparm=1,nparmc_read)
! Note: Fock terms yet to be put in ijas=4,5,6.
      endif

      close(3)

      call basis_norm(iwt,1)
  400 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine wf_copy
! Written by Cyrus Umrigar
! Copy all the parameters that are read in, from iadd_diag=1 to iadd_diag=2,3 for use in correlated sampling
      use control_mod
      use orbitals_mod
      use atom_mod
      use coefs_mod
      use dets_mod
      use basis1_mod
      use dim_mod
      use numbas_mod
      use basis2_mod
      use contr2_mod
      use forcepar_mod
      use wfsec_mod
      use contrl_per_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use bparm_mod
      use optimo_mod
      use orbpar_mod
      use numexp_mod
      implicit real*8(a-h,o-z)

      nparma_read=2+max(0,norda-1)
      nparmb_read=2+max(0,nordb-1)
      nparmc_read=nterms4(nordc)

      nwftype=3
      istrech=0

! Copy all the parameters that are read in, from iadd_diag=1 to iadd_diag=2,3 for use in correlated sampling

! Determinants
! Warning: I need to do this for other types of basis functions, both analytic and numerical
! Warning: I have not been careful about doing it only when needed, because it does not
! matter to do extra anyway.
! Analytical basis functions
      do 42 iadd_diag=2,3
        if(inum_orb.eq.0) then
          if(ibasis.eq.1 .or. ibasis.eq.3) then
            do 20 i=1,orb_tot_nb
              do 20 j=1,nbasis
   20         coef(j,i,iadd_diag)=coef(j,i,1)
            do 30 j=1,nbasis
   30         zex(j,iadd_diag)=zex(j,1)
            do 35 ictype=1,nctype
              do 35 irwf=1,nrbas_analytical(ictype)
   35           zex2(irwf,ictype,iadd_diag)=zex2(irwf,ictype,1)
          endif
          if(ibasis.ge.4 .and. ibasis.le.7) then
            do 37 i=1,orb_tot_nb
              do 37 j=1,nbasis
   37           coef(j,i,iadd_diag)=coef(j,i,1)
          endif
        endif

        do 40 i=1,ncsf
   40     csf_coef(i,iadd_diag)=csf_coef(i,1)

        do 42 it=1,notype
          do 42 ip=1,nbasis
   42       oparm(it,ip,iadd_diag)=oparm(it,ip,1)

! Analytical orbitals and numerical radial basis functions
      if(inum_orb.eq.0 .and. (ibasis.eq.1.or.ibasis.eq.3)) then
      do 45 iadd_diag=2,nwftype
        do 45 ict=1,nctype
          do 45 irb=1,nrbas_numerical(ict)
            do 45 ir=1,nr(ict)
              rwf(ir,irb,ict,iadd_diag)=rwf(ir,irb,ict,1)
   45         d2rwf(ir,irb,ict,iadd_diag)=d2rwf(ir,irb,ict,1)

! Warning: Although ae and ce have MFORCE rather than MWF in the dimensions, I think they should be MWF
      do 50 ifr=2,nforce
        do 50 ict=1,nctype
          do 50 irb=1,nrbas_numerical(ict)
            do 50 icoef=1,NCOEF
              if(icoef.le.2) ae(icoef,irb,ict,ifr)=ae(icoef,irb,ict,1)
   50         ce(icoef,irb,ict,ifr)=ce(icoef,irb,ict,1)

      endif

! Jastrow
      call alloc ('iwftype', iwftype, max(3,nforce))
      do 90 iadd_diag=2,3
      iwftype(iadd_diag)=iadd_diag
      scalek(iadd_diag)=scalek(1)
      do 52 ict=1,nctype
        do 52 i=1,nparma_read
   52     a4(i,ict,iadd_diag)=a4(i,ict,1)
      do 54 isp=nspin1,nspin2b
        do 54 i=1,nparmb_read
   54     b(i,isp,iadd_diag)=b(i,isp,1)
      do 56 ict=1,nctype
        do 56 i=1,nparmc_read
   56     c(i,ict,iadd_diag)=c(i,ict,1)

! Warning: Is this call to basis_norm needed?
!  90 call basis_norm(iadd_diag,1)
   90 continue

!     Copy orbital coefficients
      if (inum_orb == 0) then
        do iadd_diag = 2, 3
          coef_orb_on_norm_basis (:,:,iadd_diag) = coef_orb_on_norm_basis (:,:,1)
        enddo

        if (l_opt_exp .and. trim(basis_functions_varied) == 'orthonormalized') then
          do iadd_diag = 2, 3
            coef_orb_on_ortho_basis (:,:,iadd_diag) = coef_orb_on_ortho_basis (:,:,1)
          enddo
        endif
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine wf_save
! Written by Cyrus Umrigar
      use all_tools_mod
      use nuclei_mod
      use basis_mod
      use orbitals_mod
      use atom_mod
      use coefs_mod
      use dets_mod
      use basis1_mod
      use dim_mod
      use basis2_mod
      use contr2_mod
      use contrl_per_mod
      use forcepar_mod
      use wfsec_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use bparm_mod
      use optimo_mod
      use orbpar_mod
      use numbas_mod
      implicit real*8(a-h,o-z)
!     common /contrl_per/ iperiodic,ibasis

      save scalek_sav,nparma_read,nparmb_read,nparmc_read

      nparma_read=2+max(0,norda-1)
      nparmb_read=2+max(0,nordb-1)
      nparmc_read=nterms4(nordc)

!cDeterminants
!cWarning: I need to do this for other types of basis functions, both analytic and numerical
! Warning: Also I need to edit and uncomment foll. lines
!     if(ibasis.eq.1 .or. ibasis.eq.3) then
!       do 20 i=1,norb
!         do 20 j=1,nbasis
!  20     coef(j,i,iadd_diag)=coef(j,i,1)
!       do 30 j=1,nbasis
!  30     zex(j,iadd_diag)=zex(j,1)
!       do 35 ictype=1,nctype
!         do 35 irwf=1,MRWF
!  35       zex2(irwf,ictype,iadd_diag)=zex2(irwf,ictype,1)
!     endif
!     do 40 i=1,ndet
!  40 cdet(i,iadd_diag)=cdet(i,1)

      call alloc ('csf_coef_sav', csf_coef_sav, ncsf)
      do 50 i=1,ncsf
   50   csf_coef_sav(i)=csf_coef(i,1)

      call alloc ('oparm_sav', oparm_sav, notype, nbasis)
      do 51 it=1,notype
        do 51 ip=1,nbasis
   51     oparm_sav(it,ip)=oparm(it,ip,1)

! Jastrow
      call alloc ('a4_sav', a4_sav, nparma_read, nctype)
      call alloc ('b_sav', b_sav, nparmb_read, nspin2b-nspin1+1)
      call alloc ('c_sav', c_sav, nparmc_read, nctype)
      scalek_sav=scalek(1)
      do 52 ict=1,nctype
        do 52 i=1,nparma_read
   52     a4_sav(i,ict)=a4(i,ict,1)
      do 54 isp=nspin1,nspin2b
        do 54 i=1,nparmb_read
   54     b_sav(i,isp)=b(i,isp,1)
      do 56 ict=1,nctype
        do 56 i=1,nparmc_read
   56     c_sav(i,ict)=c(i,ict,1)

!     save orbital coefficients
      if (inum_orb == 0) then
       call object_provide ('nbasis')
       call object_provide ('orb_tot_nb')

       call object_provide ('coef')
       call object_alloc ('coef_sav', coef_sav, nbasis, orb_tot_nb)
       coef_sav (1:nbasis, 1:orb_tot_nb) = coef (1:nbasis, 1:orb_tot_nb, 1)
       call object_modified ('coef_sav')

       call object_provide ('coef_orb_on_norm_basis')
       call object_alloc ('coef_orb_on_norm_basis_sav', coef_orb_on_norm_basis_sav, nbasis, orb_tot_nb)
       coef_orb_on_norm_basis_sav (1:nbasis, 1:orb_tot_nb) = coef_orb_on_norm_basis (1:nbasis, 1:orb_tot_nb, 1)
       call object_modified ('coef_orb_on_norm_basis_sav')

       if (l_opt_exp .and. trim(basis_functions_varied) == 'orthonormalized') then
        call object_provide ('coef_orb_on_ortho_basis')
        call object_alloc ('coef_orb_on_ortho_basis_sav', coef_orb_on_ortho_basis_sav, nbasis, orb_tot_nb)
        coef_orb_on_ortho_basis_sav (1:nbasis, 1:orb_tot_nb) = coef_orb_on_ortho_basis (1:nbasis, 1:orb_tot_nb, 1)
        call object_modified ('coef_orb_on_ortho_basis_sav')
       endif

!      save exponents
!      if (numr.le.0 .and. (ibasis.eq.1 .or. ibasis.eq.3)) then
       if (maxval(zex).ne.0.d0 .and. (ibasis.eq.1 .or. ibasis.eq.3)) then
         call object_provide ('nctype')
         call object_alloc ('zex_sav', zex_sav, nbasis)
         call object_alloc ('zex2_sav', zex2_sav, MRWF, nctype)
         zex_sav (1:nbasis) = zex (1:nbasis,1)
!        zex2_sav (1:MRWF,1:nctype) = zex2 (1:MRWF,1:nctype,1)
         do ict=1,nctype
           zex2_sav (1:nrbas_analytical(ict),1:nctype) = zex2 (1:nrbas_analytical(ict),1:nctype,1)
         enddo
         call object_modified ('zex_sav')
         call object_modified ('zex2_sav')
       endif
      endif

!     save nuclear coordinates
      call object_provide ('ndim')
      call object_provide ('ncent')
      call object_provide ('cent')
      call object_alloc ('cent_sav', cent_sav, ndim, ncent)
      cent_sav (1:ndim,1:ncent) = cent (1:ndim,1:ncent)
      call object_modified ('cent_sav')

      return
!-----------------------------------------------------------------------
      entry wf_restore

      write(6,*)
      write(6,'(a)') 'Restoring wavefunction.'

      do 150 i=1,ncsf
  150   csf_coef(i,1)=csf_coef_sav(i)

      do 151 it=1,notype
        do 151 ip=1,nbasis
  151     oparm(it,ip,1)=oparm_sav(it,ip)

      scalek(1)=scalek_sav
      do 152 ict=1,nctype
        do 152 i=1,nparma_read
  152     a4(i,ict,1)=a4_sav(i,ict)
      do 154 isp=nspin1,nspin2b
        do 154 i=1,nparmb_read
  154     b(i,isp,1)=b_sav(i,isp)
      do 156 ict=1,nctype
        do 156 i=1,nparmc_read
  156     c(i,ict,1)=c_sav(i,ict)

!     restore orbital coefficients (must restore only for iwf=1 to avoid problems)
      if (inum_orb == 0) then
        call object_valid_or_die ('coef_sav')
        coef (1:nbasis, 1:orb_tot_nb, 1) = coef_sav (1:nbasis, 1:orb_tot_nb)
        call object_modified ('coef')

        call object_valid_or_die ('coef_orb_on_norm_basis_sav')
        coef_orb_on_norm_basis (1:nbasis, 1:orb_tot_nb, 1) = coef_orb_on_norm_basis_sav (1:nbasis, 1:orb_tot_nb)
        call object_modified ('coef_orb_on_norm_basis')

        if (l_opt_exp .and. trim(basis_functions_varied) == 'orthonormalized') then
         call object_valid_or_die ('coef_orb_on_ortho_basis_sav')
         coef_orb_on_ortho_basis (1:nbasis, 1:orb_tot_nb, 1) = coef_orb_on_ortho_basis_sav (1:nbasis, 1:orb_tot_nb)
         call object_modified ('coef_orb_on_ortho_basis')
        endif

!       restore exponents (must restore only for iwf=1 to avoid problems)
!       if (numr.le.0 .and. (ibasis.eq.1 .or. ibasis.eq.3)) then
        if (maxval(zex).ne.0.d0 .and. (ibasis.eq.1 .or. ibasis.eq.3)) then
          call object_valid_or_die ('zex_sav')
          call object_valid_or_die ('zex2_sav')
          zex (1:nbasis,1) = zex_sav (1:nbasis)
!         zex2 (1:MRWF,1:nctype,1) = zex2_sav (1:MRWF,1:nctype)
          do ict=1,nctype
            zex2 (1:nrbas_analytical(ict),1:nctype,1) = zex2_sav (1:nrbas_analytical(ict),1:nctype)
          enddo
          call object_modified ('zex')
          call object_modified ('zex2')
        endif
      endif

!     restore nuclear coordinates
      call object_valid_or_die ('cent_sav')
      cent (1:ndim,1:ncent) = cent_sav (1:ndim,1:ncent)
      call object_modified ('cent')

      return
      end
!-----------------------------------------------------------------------
      subroutine wf_best_save

! Written by Cyrus Umrigar
! Saves the best wf yet and writes it out at end of run
      use all_tools_mod
      use orbitals_mod
      use basis_mod
      use optimization_mod
      use optim_mod, only : nparml,nparme,nparmcsf,nparms,nparmg,nparmj
      use atom_mod
      use coefs_mod
      use dets_mod
      use basis1_mod
      use dim_mod
      use basis2_mod
      use contr2_mod
      use forcepar_mod
      use wfsec_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use bparm_mod
      use pjas_mod, only: param_pjas_nb, pjas_parms, pjas_parms_best
      use optimo_mod
      use orbpar_mod
      use contrl_per_mod
      use numbas_mod
      implicit real*8(a-h,o-z)

      character*30 fmt

!     common /contrl_per/ iperiodic,ibasis

      save nparma_read,nparmb_read,nparmc_read

      nparma_read=2+max(0,norda-1)
      nparmb_read=2+max(0,nordb-1)
      nparmc_read=nterms4(nordc)

!cDeterminants
!cWarning: I need to do this for other types of basis functions, both analytic and numerical
! Warning: Also I need to edit and uncomment foll. lines and then I need common /coefs/ too.
!     if(ibasis.eq.1 .or. ibasis.eq.3) then
!       do 20 i=1,norb
!         do 20 j=1,nbasis
!  20     coef(j,i,iadd_diag)=coef(j,i,1)
!       do 30 j=1,nbasis
!  30     zex(j,iadd_diag)=zex(j,1)
!       do 35 ictype=1,nctype
!         do 35 irwf=1,MRWF
!  35       zex2(irwf,ictype,iadd_diag)=zex2(irwf,ictype,1)
!     endif
!     do 40 i=1,ndet
!  40 cdet(i,iadd_diag)=cdet(i,1)

      call alloc ('csf_coef_best', csf_coef_best, ncsf)
      do 50 i=1,ncsf
   50   csf_coef_best(i)=csf_coef(i,1)

      call object_modified ('csf_coef_best')

      call alloc ('oparm_best', oparm_best, notype, nbasis)
      do 51 it=1,notype
        do 51 ip=1,nbasis
   51     oparm_best(it,ip)=oparm(it,ip,1)

! Jastrow
      call alloc ('a4_best', a4_best, nparma_read, nctype)
      call alloc ('b_best', b_best, nparmb_read, nspin2b-nspin1+1)
      call alloc ('c_best', c_best, nparmc_read, nctype)
      scalek_best=scalek(1)
      do 52 ict=1,nctype
        do 52 i=1,nparma_read
   52     a4_best(i,ict)=a4(i,ict,1)
      do 54 isp=nspin1,nspin2b
        do 54 i=1,nparmb_read
   54     b_best(i,isp)=b(i,isp,1)
      do 56 ict=1,nctype
        do 56 i=1,nparmc_read
   56     c_best(i,ict)=c(i,ict,1)

      call object_modified ('a4_best')
      call object_modified ('b_best')
      call object_modified ('c_best')

! Periodic Jastrow
      if(param_pjas_nb.gt.0) then
        call alloc ('pjas_parms_best', pjas_parms_best, param_pjas_nb)
        pjas_parms_best(1:param_pjas_nb)=pjas_parms(1:param_pjas_nb,1)
      endif

      call object_modified ('pjas_parms_best')

! Orbital coefficients
      if (inum_orb == 0) then
        call object_provide ('nbasis')
        call object_provide ('orb_tot_nb')

        call object_provide ('coef')
        call object_alloc ('coef_best', coef_best, nbasis, orb_tot_nb)
        coef_best (1:nbasis, 1:orb_tot_nb) = coef (1:nbasis, 1:orb_tot_nb, 1)
        call object_modified ('coef_best')


        call object_provide ('coef_orb_on_norm_basis')
        call object_alloc ('coef_orb_on_norm_basis_best', coef_orb_on_norm_basis_best, nbasis, orb_tot_nb)
        coef_orb_on_norm_basis_best (1:nbasis, 1:orb_tot_nb) = coef_orb_on_norm_basis (1:nbasis, 1:orb_tot_nb, 1)
        call object_modified ('coef_orb_on_norm_basis_best')

        if (l_opt_exp .and. trim(basis_functions_varied) == 'orthonormalized') then
         call object_provide ('coef_orb_on_ortho_basis')
         call object_alloc ('coef_orb_on_ortho_basis_best', coef_orb_on_ortho_basis_best, nbasis, orb_tot_nb)
         coef_orb_on_ortho_basis_best (1:nbasis, 1:orb_tot_nb) = coef_orb_on_ortho_basis (1:nbasis, 1:orb_tot_nb, 1)
         call object_modified ('coef_orb_on_ortho_basis_best')
        endif

!       save exponents
!       if (numr.le.0) then
!       if (maxval(zex).ne.0.d0) then
        if (maxval(zex).ne.0.d0 .and. (ibasis.eq.1 .or. ibasis.eq.3)) then
          call object_provide ('nctype')
          call object_alloc ('zex_best', zex_best, nbasis)
          call object_alloc ('zex2_best', zex2_best, MRWF, nctype)
          zex_best (1:nbasis) = zex (1:nbasis,1)
          do ict=1,nctype
           zex2_best (1:nrbas_analytical(ict),1:nctype) = zex2 (1:nrbas_analytical(ict),1:nctype,1)
          enddo
          call object_modified ('zex_best')
          call object_modified ('zex2_best')
        endif
      endif

!     save nuclear coordinates
      call object_provide ('ndim')
      call object_provide ('ncent')
      call object_provide ('cent')
      call object_alloc ('cent_best', cent_best, ndim, ncent)
      cent_best (1:ndim,1:ncent) = cent (1:ndim,1:ncent)
      call object_modified ('cent_best')

      return
!-----------------------------------------------------------------------
      entry wf_best_restore

      write(6,*)
      write(6,'(a)') 'Restoring best wavefunction.'

      do 150 i=1,ncsf
  150   csf_coef(i,1)=csf_coef_best(i)

      do 151 it=1,notype
        do 151 ip=1,nbasis
  151     oparm(it,ip,1)=oparm_best(it,ip)

      scalek(1)=scalek_best
      do 152 ict=1,nctype
        do 152 i=1,nparma_read
  152     a4(i,ict,1)=a4_best(i,ict)
      do 154 isp=nspin1,nspin2b
        do 154 i=1,nparmb_read
  154     b(i,isp,1)=b_best(i,isp)
      do 156 ict=1,nctype
        do 156 i=1,nparmc_read
  156     c(i,ict,1)=c_best(i,ict)

! Periodic Jastrow
      if(param_pjas_nb.gt.0) then
        pjas_parms(1:param_pjas_nb,1)=pjas_parms_best(1:param_pjas_nb)
      endif

!     restore orbital coefficients (must restore only for iwf=1 to avoid problems)
      if (inum_orb == 0) then
        call object_valid_or_die ('coef_best')
        coef (1:nbasis, 1:orb_tot_nb, 1) = coef_best (1:nbasis, 1:orb_tot_nb)
        call object_modified ('coef')

        call object_valid_or_die ('coef_orb_on_norm_basis_best')
        coef_orb_on_norm_basis (1:nbasis, 1:orb_tot_nb, 1) = coef_orb_on_norm_basis_best (1:nbasis, 1:orb_tot_nb)
        call object_modified ('coef_orb_on_norm_basis')

        if (l_opt_exp .and. trim(basis_functions_varied) == 'orthonormalized') then
         call object_valid_or_die ('coef_orb_on_ortho_basis_best')
         coef_orb_on_ortho_basis (1:nbasis, 1:orb_tot_nb, 1) = coef_orb_on_ortho_basis_best (1:nbasis, 1:orb_tot_nb)
         call object_modified ('coef_orb_on_ortho_basis')
        endif


!       restore exponents (must restore only for iwf=1 to avoid problems)
!       if (numr.le.0 .and. (ibasis.eq.1 .or. ibasis.eq.3)) then
        if (maxval(zex).ne.0.d0 .and. (ibasis.eq.1 .or. ibasis.eq.3)) then
          call object_valid_or_die ('zex_best')
          call object_valid_or_die ('zex2_best')
          zex (1:nbasis,1) = zex_best (1:nbasis)
          do ict=1,nctype
            zex2 (1:nrbas_analytical(ict),1:nctype,1) = zex2_best (1:nrbas_analytical(ict),1:nctype)
          enddo
          call object_modified ('zex')
          call object_modified ('zex2')
        endif
      endif

!     restore nuclear coordinates
      call object_valid_or_die ('cent_best')
      cent (1:ndim,1:ncent) = cent_best (1:ndim,1:ncent)
      call object_modified ('cent')

      return
!-----------------------------------------------------------------------
      entry wf_best_write

      close(2)
      if(idtask.eq.0) then
        open(2,file='wavefn_best')
       else
        open(2,file='/dev/null')
      endif

      write(6,'(/,''Best wave function:'')')

      if(nparmcsf.gt.0) then
        write(fmt,'(''(''i2,''f15.8,a)'')') ncsf
        write(6,fmt) (csf_coef_best(i),i=1,ncsf),' (csf_coef_best(icsf),icsf=1,ncsf)'
        write(2,fmt) (csf_coef_best(i),i=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      endif

      if(nparmot.gt.0) then
        do it=1,notype
          write(fmt,'(''(''i2,''f15.8,a)'')') nbasis
          if(ibasis.eq.5) then
            if(it.eq.1) then
              write(6,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_rad_pos_best(it,i),i=1,nbasis)'
              write(2,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_rad_pos_best(it,i),i=1,nbasis)'
            elseif(it.eq.2) then
              write(6,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_ang_pos_best(it,i),i=1,nbasis)'
              write(2,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_ang_pos_best(it,i),i=1,nbasis)'
            elseif(it.eq.3) then
              write(6,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_rad_width_best(it,i),i=1,nbasis)'
              write(2,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_rad_width_best(it,i),i=1,nbasis)'
            elseif(it.eq.4) then
              write(6,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_ang_width_best(it,i),i=1,nbasis)'
              write(2,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_ang_width_best(it,i),i=1,nbasis)'
            endif
          else
            if(it.eq.1) then
              write(6,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_x_pos_best(it,i),i=1,nbasis)'
              write(2,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_x_pos_best(it,i),i=1,nbasis)'
            elseif(it.eq.2) then
              write(6,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_y_pos_best(it,i),i=1,nbasis)'
              write(2,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_y_pos_best(it,i),i=1,nbasis)'
            elseif(it.eq.3) then
              write(6,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_x_width_best(it,i),i=1,nbasis)'
              write(2,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_x_width_best(it,i),i=1,nbasis)'
            elseif(it.eq.4) then
              write(6,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_y_width_best(it,i),i=1,nbasis)'
              write(2,fmt) (oparm_best(it,i),i=1,nbasis),' (floating_gauss_y_width_best(it,i),i=1,nbasis)'
            endif
          endif
        enddo
      endif

      if(nparms.gt.0) then
        write(6,'(/,''scalek parameter:'')')
        write(6,'(f9.6,'' 0. scalek_best,a21'')') scalek_best
        write(2,'(f9.6,'' 0. scalek_best,a21'')') scalek_best
      endif

      if(nparmj.gt.0) then
        write(6,'(a)') 'Jastrow parameters:'
        write(2,'(a)') 'Jastrow parameters:'
        nparma_read=2+max(0,norda-1)
        nparmb_read=2+max(0,nordb-1)
        nparmc_read=nterms4(nordc)

        if(nparma_read.gt.0) then
!         write(fmt,'(''(''i2,''f15.8,\'\' (a(iparmj),iparmj=1,nparma)\'\')'')') nparma_read
          write(fmt,'(''(1p'',i2,''g22.14,a)'')') nparma_read
         else
!         write(fmt,'(''(\'\' (a(iparmj),iparmj=1,nparma)\'\')'')')
          write(fmt,'(''(a)'')')
        endif
        do 80 ict=1,nctype
          write(6,fmt) (a4_best(i,ict),i=1,nparma_read),' (a_best(iparmj),iparmj=1,nparma)'
   80     write(2,fmt) (a4_best(i,ict),i=1,nparma_read),' (a(iparmj),iparmj=1,nparma)'

        if(nparmb_read.gt.0) then
!         write(fmt,'(''(''i2,''f15.8,\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')') nparmb_read
          write(fmt,'(''(1p'',i2,''g22.14,a)'')') nparmb_read
         else
!         write(fmt,'(''(\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')')
          write(fmt,'(''(a)'')')
        endif
        do 85 isp=nspin1,nspin2b
          write(6,fmt) (b_best(i,isp),i=1,nparmb_read),' (b_best(iparmj),iparmj=1,nparmb)'
   85     write(2,fmt) (b_best(i,isp),i=1,nparmb_read),' (b(iparmj),iparmj=1,nparmb)'

        if(nparmc_read.gt.0) then
!         write(fmt,'(''(''i2,''f15.8,\'\' (c(iparmj),iparmj=1,nparmc)\'\')'')') nparmc_read
          write(fmt,'(''(1p'',i2,''g22.14,a)'')') nparmc_read
         else
!         write(fmt,'(''(\'\' (c(iparmj),iparmj=1,nparmc)\'\')'')')
          write(fmt,'(''(a)'')')
        endif
        do 90 ict=1,nctype
          write(6,fmt) (c_best(i,ict),i=1,nparmc_read),' (c_best(iparmj),iparmj=1,nparmc)'
   90     write(2,fmt) (c_best(i,ict),i=1,nparmc_read),' (c(iparmj),iparmj=1,nparmc)'
      endif

      return
      end

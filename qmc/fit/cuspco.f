      subroutine cuspco(diff,iprin)
c Written by Cyrus Umrigar
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:::  If icusp.ge.0 impose the cusp condition for each s orbital.    :::
c:::  Now we do this exactly, previously it was done exactly for     :::
c:::  atoms but iteratively for molecules.                           :::
c:::  In any case calculate cusp penalty for 0th order e-N cusp.     :::
c:::  If icusp.ge.0 then this just serves to confirm that we imposed :::
c:::  the cusp exactly.                                              :::
c:::  The number of these cusp conditions is ncent*norb              :::
!
! J. Toulouse - 06 Jan 05:
! l_impose_cusp_en_occ = true : impose cusp conditions only on occupied orbitals
! J. Toulouse - 08 Jan 05: change coef(i,j,1) -> coef(i,j,iwf)
!                        : change a(i,j,1) -> a(i,j,iwf)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use atom_mod
      use orbitals_mod ! JT
      implicit real*8(a-h,o-z)


      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      common /optim/ lo(MORB),npoint(MORB),
     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
     &imnbas(MCENT),
     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
     &necn,nebase

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /coefs2/ coef2(MBASIS,MORB,MCENT)

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE)
     &,n2s(MCTYPE),n2p(-1:1,MCTYPE)
     &,n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
     &,n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE),n4f(-3:3,MCTYPE)
     &,n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
     &,n5g(-4:4,MCTYPE)
     &,n6d(-2:2,MCTYPE),n6f(-3:3,MCTYPE),n6g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,n7g(-4:4,MCTYPE),n7h(-5:5,MCTYPE),n7i(-6:6,MCTYPE)
     &,n8i(-6:6,MCTYPE),n8j(-7:7,MCTYPE)
     &,n9k(-8:8,MCTYPE)
     &,n10l(-9:9,MCTYPE)
     &,n11m(-10:10,MCTYPE)
     &,n12n(-11:11,MCTYPE)
     &,n13o(-12:12,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ zex2(MRWF,MCTYPE,MWF),n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),iwrwf2(MBASIS)
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar2/ a1(MPARMJ,3,MWF),a2(MPARMJ,3,MWF)
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn

      common /wfsec/ iwftype(MFORCE),iwf,nwftype !JT

      dimension diff(*),orb(MORB),orb2(MORB)

!      data zero,one,three/0.d0,1.d0,3.d0/
!      data half,d1b4/0.5d0,.25d0/
      data d1b4/.25d0/

c     write(6,'(''n_bas,l_bas='',30(2i2,2x))') (n_bas(i),l_bas(i),i=1,nbasis)
c The shifting is done in the calling routine, so do not shift here
      ishft=0

c Figure out part of nuclear cusp from the Jastrow factor.
c We assume that the coef. of t, a1(3,i)=0 as required by symmetry.
c The 0th order e-N cusp-cond is
c ndn*a1(2,1)+(nup-1)*a1(2,2) = nup*a1(2,1)+(ndn-1)*a1(2,3)
c So, when nspin2=3, only 2 of the 3 a1(2,isp) coefs are independent.
c If nspin2<=2 so that a1(i,3)=a1(i,2) and ndn.ne.nup then this requires
c that a1(2,2)=a1(2,1)
c This can only be consistent with a1(i,2)=a1(i,3)=a1(i,1)/2 (nspin2=1) if
c a1(2,1)=a1(2,2)=a1(2,3)=0
c So, if ndn.ne.nup and nspin2=1, we must fix the a1(2,isp)=0.
c If ndn=nup, then there is no restriction from the 0th order e-N cusp
c and one is free to use a non-zero a1(i,2)=a1(i,3)=a1(i,1)/2 when nspin2=1.

c If nbasis=1 and nsa=1, aa1 is expressed in terms of zex. For ijas=2,
c impose in the input a1(2,1) (singlet) or a1(2,2) (triplet) equal to aa1
c and do not vary it.
      if(ijas.eq.2) then
        if(nspin1.eq.1 .and. nspin2.eq.1) then
          if(nelec.eq.2) then
            aa1=a1(2,1,iwf)
           else
            aa1=d1b4*(three*(nup+ndn)-2)*a1(2,1,iwf)
          endif
         elseif(nspin1.eq.2 .and. nspin2.eq.2) then
          aa1=(nup-1)*a1(2,2,iwf)
         elseif(nspin1.eq.1 .and. nspin2.eq.2) then
          aa1=half*((nup+ndn)*a1(2,1,iwf)+(nup+ndn-2)*a1(2,2,iwf))
         elseif(nspin1.eq.1 .and. nspin2.eq.3) then
          aa1=nup*a1(2,1,iwf)+(ndn-1)*a1(2,3,iwf)
        endif
       elseif(ijas.eq.3) then
        aa1=a(1,1)
       elseif(ijas.ge.4.and.ijas.le.6) then
c       write(6,*) a4(1,1,1)
        if(nctype.gt.1.and.a4(1,1,iwf).ne.0.d0) stop 'cuspco for ijas=4,5,6 presently for nctype=1 only'
        aa1=a4(1,1,iwf)
       else
        aa1=zero
      endif

c Impose e-N cusp conds for each atom and each orbital.
c This was previously being done iteratively but is now being
c done exactly.  However, it assumes that iebasi and ieorb
c are such that the coefs of the higher numbered atoms are being set
c equal to those of the lower numbered ones, rather than the reverse.
c All our inputs are in this form anyway.
c First, equiv_bas is called from func to set coef2 to coef for basis functions
c that are symmetry related to the one being adjusted and 0 for others. Then
c cusorb calculates the value of the orbital at each nuclear position and
c cusorb_equiv calculates the value of the part of the orbital coming
c from the coefficient that is being adjusted and the symmetry related coeficients
c at each nuclear position.

c         c_i*xi_i + sum_{j=same atom other 1s bases} c_j*xi_j + sum_{j=same atom 2s bases} c_j
c (Z+a) = ---------------------------------------------------------------------------------------------------------------------
c         c_i(1+psi_{other_atom,equiv_bases}/c_i) + sum_{j=other_atom,inequiv_bases} psi + sum_{j=same atom other 1s bases} c_j

c         t + (Z+a)*psi_{other_atom,inequiv_bases}
c   c_i = ----------------------------------------
c         xi_i - (Z+a)*(1+t2)
c where
c     t = sum_{j=same atom other 1s bases} (Z+a-xi_j)*c_j + sum_{j=same atom 2s bases} c_j
c    t2 = psi_{other_atom,equiv_bases}/c_i

      if(icusp.ge.0) then
        do 5 i=1,necn
    5     coef(iebasi(1,i),ieorb(1,i),iwf)=sign(one,dfloat(ieorb(2,i)))
     &    *coef(iebasi(2,i),iabs(ieorb(2,i)),iwf)
        do 10 i=1,nebase
   10     zex(iebase(1,i),iwf)=zex(iebase(2,i),iwf)
        do 35 icent=1,ncent
         if(icent.eq.ncent) then
           imxbas=nbasis
          else
           imxbas=imnbas(icent+1)-1
         endif
         ibas=imnbas(icent)
c     write(6,'(''ibas,imnbas(icent),imxbas'',9i5)') ibas,imnbas(icent),imxbas
         call cusorb(icent,orb)
         call cusorb_equiv(icent,orb2)
         do 35 iorb=1,norb
! JT beg: skip non-occupied orbitals
           if(l_impose_cusp_en_occ) then
            call object_provide ('orb_occ_in_wf')
            if(.not. orb_occ_in_wf(iorb)) cycle
           endif
! JT end
           if(lo(iorb).eq.0) then
             if(imxbas.eq.ibas) then
               if(n_bas(ibas).eq.0) then
                 beta=betaq/zex(ibas,iwf)-one
                 aa1=zex(ibas,iwf)-znuc(iwctype(icent))-beta
                else
                 zex(ibas,iwf)=znuc(iwctype(icent))+aa1
               endif
              else
               term=zero
               orbsam=coef(ibas,iorb,iwf)
               do 20 ib=ibas+1,imxbas
                 if(l_bas(ib).eq.0) then
                  if(n_bas(ib).eq.0) then
                    beta=betaq/zex(ib,iwf)-one
                    term=term+(znuc(iwctype(icent))+aa1
     &              -(zex(ib,iwf)-beta))*coef(ib,iorb,iwf)
                    orbsam=orbsam+coef(ib,iorb,iwf)
                   elseif(n_bas(ib).eq.1) then
                    term=term+(znuc(iwctype(icent))+aa1-zex(ib,iwf))*
     &              coef(ib,iorb,iwf)
                    orbsam=orbsam+coef(ib,iorb,iwf)
                   elseif(n_bas(ib).eq.2) then
                    term=term+coef(ib,iorb,iwf)
                  endif
                 endif
   20            continue
                 other_atom=orb(iorb)-orbsam
                 other_atom_equiv_bas=orb2(iorb)-coef2(ibas,iorb,icent)
c          write(6,'(''iorb,other_atom,other_atom_equiv_bas'',i5,9f9.5)') iorb,other_atom,other_atom_equiv_bas,
c    &     orb(iorb),orbsam,orb2(iorb),coef2(ibas,iorb,icent)
                 other_atom_ineqv_bas=other_atom-other_atom_equiv_bas
                 if(coef(ibas,iorb,iwf).ne.0) then
c                  write(6,'(''ibas,iorb,coef(ibas,iorb,iwf)'',2i5,9f9.5)') ibas,iorb,coef(ibas,iorb,iwf)
                   term2=other_atom_equiv_bas/coef(ibas,iorb,iwf)
                  else
                   term2=0
                 endif
c           write(6,'(''coef(ibas,iorb,iwf)1'',9f9.5)') coef(ibas,iorb,iwf)
c           write(6,'(''zex(ibas,iwf),znuc(iwctype(icent)),aa1,term2'',9f9.5)') zex(ibas,iwf),znuc(iwctype(icent)),aa1,term2
                 coef(ibas,iorb,iwf)=
     &           (term+(znuc(iwctype(icent))+aa1)*other_atom_ineqv_bas)
     &           /(zex(ibas,iwf)-(znuc(iwctype(icent))+aa1)*(1+term2))
c           write(6,'(''coef(ibas,iorb,iwf)2'',9f9.5)') coef(ibas,iorb,iwf)
             endif
c            write(6,'(''other_atom'',3i2,9d12.5)') icent,iorb,ibas,
c    &       coef(ibas,iorb,iwf),other_atom,orb(iorb),orb2(iorb),
c    &       other_atom_equiv_bas,other_atom_ineqv_bas
           endif
           do 22 i=1,necn
   22        coef(iebasi(1,i),ieorb(1,i),iwf)=sign(one,dfloat(ieorb(2,i)))
     &       *coef(iebasi(2,i),iabs(ieorb(2,i)),iwf)
   35     continue
        endif

c Check 0th order e-N cusp conds. for each atom and each orbital
c If icusp.ge.0 then all the diffs should be 0.
        do 55 icent=1,ncent
         if(icent.eq.ncent) then
           imxbas=nbasis
          else
           imxbas=imnbas(icent+1)-1
         endif
         ibas=imnbas(icent)
         call cusorb(icent,orb)
         do 50 iorb=1,norb
! JT beg: skip non-occupied orbitals
           if(l_impose_cusp_en_occ) then
            call object_provide ('orb_occ_in_wf')
            if(.not. orb_occ_in_wf(iorb)) cycle
           endif
! JT end
           if(lo(iorb).eq.0) then
             top=zero
             bot=zero
             do 40 ib=ibas,imxbas
               if(l_bas(ib).eq.0) then
                  if(n_bas(ib).eq.0) then
                   beta=betaq/zex(ib,iwf)-one
                   top=top-(zex(ib,iwf)-beta)*coef(ib,iorb,iwf)
                   bot=bot+coef(ib,iorb,iwf)
                  elseif(n_bas(ib).eq.1) then
                   top=top-zex(ib,iwf)*coef(ib,iorb,iwf)
                   bot=bot+coef(ib,iorb,iwf)
                  elseif(n_bas(ib).eq.2) then
                   top=top+coef(ib,iorb,iwf)
                 endif
               endif
 40          continue
             differ=top+orb(iorb)*(znuc(iwctype(icent))+aa1)
c      write(6,'(''top,bot,orb(iorb),znuc(iwctype(icent)),aa1'',9d12.4)')
c    & top,bot,orb(iorb),znuc(iwctype(icent)),aa1

!             if (iprin.ge.1) then
!              write(6,'(''coef  of cnst - znuc'',5x,f12.6
!     &       ,4x,''orb'',i2,3x,''atom'',i2)')
!     &             differ, iorb,icent
!             endif
!             if(abs(differ).ge.1.d-13) write(6,'(''differ'',2i2,d12.4)')
!     &       icent,iorb,differ
             write(6,'(a,i2,a,i2,a,f12.6)') 'atom # ',icent,' orbital # ',iorb,' : log_deriv + znuc =',differ
             diff(ishft+1) = differ

             ishft=ishft+1
           endif
 50      continue
 55     continue

      return
      end

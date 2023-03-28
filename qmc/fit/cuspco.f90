      subroutine cuspco(diff,iprin)
! Written by Cyrus Umrigar
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::  If icusp.ge.0 impose the cusp condition for each s orbital.    :::
!:::  Now we do this exactly, previously it was done exactly for     :::
!:::  atoms but iteratively for molecules.                           :::
!:::  In any case calculate cusp penalty for 0th order e-N cusp.     :::
!:::  If icusp.ge.0 then this just serves to confirm that we imposed :::
!:::  the cusp exactly.                                              :::
!:::  The number of these cusp conditions is ncent*norbc             :::
!:::  where norbc is the number of s-like orbitals (occup or total   :::
!:::  depending on whether l_impose_cusp_en_occ = true or not.       :::
!
! J. Toulouse - 06 Jan 05:
! l_impose_cusp_en_occ = true : impose cusp conditions only on occupied orbitals
! J. Toulouse - 08 Jan 05: change coef(i,j,1) -> coef(i,j,iwf)
!                        : change a(i,j,1) -> a(i,j,iwf)
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use constants_mod
      use atom_mod
      use orbitals_mod
      use cusp_mod
      use coefs_mod
      use coefs2_mod
      use dets_mod
      use optim_mod
      use basis1_mod
      use const_mod
      use basis2_mod
      use contr2_mod
      use wfsec_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar2_mod
      implicit real*8(a-h,o-z)

      dimension diff(*),orb(norb),orb2(norb)

      data d1b4/.25d0/

!     write(6,'(''n_bas,l_bas='',30(2i2,2x))') (n_bas(i),l_bas(i),i=1,nbasis)
! The shifting is done in the calling routine, so do not shift here
      ishft=0

! Figure out part of nuclear cusp from the Jastrow factor.
! We assume that the coef. of t, a1(3,i)=0 as required by symmetry.
! The 0th order e-N cusp-cond is
! ndn*a1(2,1)+(nup-1)*a1(2,2) = nup*a1(2,1)+(ndn-1)*a1(2,3)
! So, when nspin2=3, only 2 of the 3 a1(2,isp) coefs are independent.
! If nspin2<=2 so that a1(i,3)=a1(i,2) and ndn.ne.nup then this requires
! that a1(2,2)=a1(2,1)
! This can only be consistent with a1(i,2)=a1(i,3)=a1(i,1)/2 (nspin2=1) if
! a1(2,1)=a1(2,2)=a1(2,3)=0
! So, if ndn.ne.nup and nspin2=1, we must fix the a1(2,isp)=0.
! If ndn=nup, then there is no restriction from the 0th order e-N cusp
! and one is free to use a non-zero a1(i,2)=a1(i,3)=a1(i,1)/2 when nspin2=1.

! If nbasis=1 and nsa=1, aa1 is expressed in terms of zex. For ijas=2,
! impose in the input a1(2,1) (singlet) or a1(2,2) (triplet) equal to aa1
! and do not vary it.
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
!       write(6,*) a4(1,1,1)
        if(nctype.gt.1.and.a4(1,1,iwf).ne.0.d0) stop 'cuspco for ijas=4,5,6 presently for nctype=1 only'
        aa1=a4(1,1,iwf)
       else
        aa1=zero
      endif

! Impose e-N cusp conds for each atom and each orbital.
! This was previously being done iteratively but is now being
! done exactly.  However, it assumes that iebasi and ieorb
! are such that the coefs of the higher numbered atoms are being set
! equal to those of the lower numbered ones, rather than the reverse.
! All our inputs are in this form anyway.
! First, equiv_bas is called from func to set coef2 to coef for basis functions
! that are symmetry related to the one being adjusted and 0 for others. Then
! cusorb calculates the value of the orbital at each nuclear position and
! cusorb_equiv calculates the value of the part of the orbital coming
! from the coefficient that is being adjusted and the symmetry related coeficients
! at each nuclear position.

!         c_i*xi_i + sum_{j=same atom other 1s bases} c_j*xi_j + sum_{j=same atom 2s bases} c_j
! (Z+a) = ---------------------------------------------------------------------------------------------------------------------
!         c_i(1+psi_{other_atom,equiv_bases}/c_i) + sum_{j=other_atom,inequiv_bases} psi + sum_{j=same atom other 1s bases} c_j

!         t + (Z+a)*psi_{other_atom,inequiv_bases}
!   c_i = ----------------------------------------
!         xi_i - (Z+a)*(1+t2)
! where
!     t = sum_{j=same atom other 1s bases} (Z+a-xi_j)*c_j + sum_{j=same atom 2s bases} c_j
!    t2 = psi_{other_atom,equiv_bases}/c_i

      if(icusp.ge.0) then
        do 5 i=1,necn
    5     coef(iebasi(1,i),ieorb(1,i),iwf)=sign(one,dfloat(ieorb(2,i))) &
     &    *coef(iebasi(2,i),iabs(ieorb(2,i)),iwf)
        do 10 i=1,nebase
   10     zex(iebase(1,i),iwf)=zex(iebase(2,i),iwf)
        do 36 icent=1,ncent
          if(icent.eq.ncent) then
            imxbas=nbasis
          else
            imxbas=imnbas(icent+1)-1
          endif
          ibas=imnbas(icent)
!         write(6,'(''ibas,imnbas(icent),imxbas'',9i5)') ibas,imnbas(icent),imxbas
          call cusorb(icent,orb)
          call cusorb_equiv(icent,orb2)

          do 35 iorb=1,norb

            if(abs(orb(iorb)).le.1e-9) cycle ! Impose en cusp conditions only on orbitals that are non-zero at that nucleus

            if(l_impose_cusp_en_occ) then ! skip unoccupied orbitals
              call object_provide ('orb_occ_in_wf')
              if(.not. orb_occ_in_wf(iorb)) cycle
            endif

            if(imxbas.eq.ibas) then ! If there is only 1 basis function on atom icent
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
                    term=term+(znuc(iwctype(icent))+aa1-(zex(ib,iwf)-beta))*coef(ib,iorb,iwf)
                    orbsam=orbsam+coef(ib,iorb,iwf)
                  elseif(n_bas(ib).eq.1) then
                    term=term+(znuc(iwctype(icent))+aa1-zex(ib,iwf))*coef(ib,iorb,iwf)
                    orbsam=orbsam+coef(ib,iorb,iwf)
                  elseif(n_bas(ib).eq.2) then
                    term=term+coef(ib,iorb,iwf)
                  endif
                endif
   20           continue
                other_atom=orb(iorb)-orbsam
                other_atom_equiv_bas=orb2(iorb)-coef2(ibas,iorb,icent)
!               write(6,'(''iorb,other_atom,other_atom_equiv_bas'',i5,9f9.5)') iorb,other_atom,other_atom_equiv_bas,orb(iorb),orbsam,orb2(iorb),coef2(ibas,iorb,icent)
                other_atom_ineqv_bas=other_atom-other_atom_equiv_bas
                if(coef(ibas,iorb,iwf).ne.0) then
!                 write(6,'(''ibas,iorb,coef(ibas,iorb,iwf)'',2i5,9f9.5)') ibas,iorb,coef(ibas,iorb,iwf)
                  term2=other_atom_equiv_bas/coef(ibas,iorb,iwf)
                else
                  term2=0
                endif
!               write(6,'(''coef(ibas,iorb,iwf)1'',9f9.5)') coef(ibas,iorb,iwf)
!               write(6,'(''zex(ibas,iwf),znuc(iwctype(icent)),aa1,term2'',9f9.5)') zex(ibas,iwf),znuc(iwctype(icent)),aa1,term2
                coef(ibas,iorb,iwf)=(term+(znuc(iwctype(icent))+aa1)*other_atom_ineqv_bas) / (zex(ibas,iwf)-(znuc(iwctype(icent))+aa1)*(1+term2))
!               write(6,'(''coef(ibas,iorb,iwf)2'',9f9.5)') coef(ibas,iorb,iwf)
            endif

!           write(6,'(''other_atom'',3i2,9d12.5)') icent,iorb,ibas,
!  &        coef(ibas,iorb,iwf),other_atom,orb(iorb),orb2(iorb),
!  &        other_atom_equiv_bas,other_atom_ineqv_bas
            do 22 i=1,necn
   22         coef(iebasi(1,i),ieorb(1,i),iwf)=sign(one,dfloat(ieorb(2,i)))*coef(iebasi(2,i),iabs(ieorb(2,i)),iwf)
   35     continue
   36   continue
      endif

! Check 0th order e-N cusp conds. for each atom and each orbital
! If icusp.ge.0 then all the diffs should be 0.
      do 55 icent=1,ncent
        if(icent.eq.ncent) then
          imxbas=nbasis
         else
          imxbas=imnbas(icent+1)-1
        endif
        ibas=imnbas(icent)
        call cusorb(icent,orb)

        do 50 iorb=1,norb

          if(abs(orb(iorb)).le.1e-9) cycle ! Impose en cusp conditions only on orbitals that are non-zero at that nucleus

          if(l_impose_cusp_en_occ) then ! skip unoccupied orbitals
            call object_provide ('orb_occ_in_wf')
            if(.not. orb_occ_in_wf(iorb)) cycle
          endif

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
   40     continue
          differ=top+orb(iorb)*(znuc(iwctype(icent))+aa1)

          if (iprin.ge.1) then
            write(6,'(''atom '',i3, '' orbital'',i3,'': log_deriv + znuc ='',f12.6)') icent,iorb,differ
          endif
          diff(ishft+1) = differ

          ishft=ishft+1
   50   continue
   55 continue

      return
      end

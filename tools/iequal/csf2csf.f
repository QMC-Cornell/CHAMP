      program csf2csf
c Written by Cyrus Umrigar
c John Lawson modified det_to_csf to deal with GENCI rather than GUGA output.
c This code has been carried over to this program, but probably it does not work here.

c Needs cleaning!!

c Since GAMESS can only deal with some point groups this program takes the
c CSFs from GAMESS and merges them to make a smaller number of CSFs that
c correspond to the higher point group symmetry.
c At first I had written det_to_csf to start with dets from GAMESS and combine
c them into CSFs, then I wrote this program to start with CSFs from GAMESS.
c It identifies the CSFs by just looking to see which pairs are
c related by certain multiplicative factors or which triplets are
c related with multiplicative factors of 1.
c It is not guaranteed to find all relationships, but it finds almost all of them.
c Also, it could find a relationship because of a chance equality of the CSF's
c that does not really exist.  I reduce the chance of that by using orbital
c eigenvalues or natural orbital occupation numbers to label determinants that could
c enter in the same CSF.
c Output written in fashion to be useful to serve as input to both old and
c new versions of optimization program, but presently only the one for the
c new version is written out after checking the triplets.
c The input can come either from the output of a GUGA or a GENCI run using
c GAMESS.  If it is a GENCI run, then gamess2qmc does not write out the
c determinants in CHAMP format, so this program does that too, but that is
c assuming that the GENCI part has been transfered correctly from det_to_csf.

c iused_orbs = 0  write out all orbs
c            = 1  write out used orbs only

      implicit real*8(a-h,o-z)
      character*80 fmt,fmt2,mode
      parameter(MDET=20000,MCSF=8000,MDET_IN_CSF=200,MELEC=100,MORB=100,MBASIS=500,MLABEL=MDET,MACT=100)
      dimension cdet(MDET),iwdet(MDET),iedet(2,MDET),frac(2,MDET)
     &,iflag_det(MDET),iflag_label_det(MLABEL),iflag_csf(MCSF), csf_coef(MCSF),ndet_in_csf(MCSF)
     &,iwdet_in_csf(MDET_IN_CSF,MCSF),cdet_in_csf(MDET_IN_CSF,MCSF)
     &,iworbd(MELEC,MDET),eigs(MORB),eigs_in_det(MELEC,MDET),label_det(MDET)
     &,coef(MBASIS,MORB),iflag_orb(MORB),map_orb(MORB)
     &,abs_csf_coef(MCSF),csf_coef_tmp(MCSF),eigs_tmp(MORB)
      character(MORB) alphain(MDET),betain(MDET),ctmp,ctmp2
      integer alpha(MDET,MACT),beta(MDET,MACT)

      read(5,*) mode
      read(5,*) iused_orbs
      write(6,'(''Reading '',a,'' mode input'')') mode
      read(5,*) cutoff_g2q,cutoff_d2c
      write(6,'(''Cutoffs in gamess2qmc and dets2csf are'',2f8.4)') cutoff_g2q,cutoff_d2c
      if(cutoff_d2c.lt.cutoff_g2q) stop 'cutoff_d2c must be >= cutoff_g2q'

      if(index(mode,'guga').ne.0) then

c eps should be something like 4.d-6 to 4.d-5
        read(5,*) ndet,nelec,nup,norb,nbasis,eps
        write(6,'(''This program is far from foolproof since it looks for near equalities'')')
        write(6,'(''You can play with eps in a range something like 1.d-7 to 1.d-5'')')
        write(6,'(''ndet,nelec,nup,norb,nbasis,eps='',5i4,1pd9.2)') ndet,nelec,nup,norb,nbasis,eps
        if(ndet.gt.MDET) stop 'ndet>MDET'
        if(nelec.gt.MELEC) stop 'nelec>MELEC'
        if(norb.gt.MORB) stop 'norb>MORB'
        if(nbasis.gt.MBASIS) stop 'nbasis>MBASIS'

        read(5,*) (cdet(i),i=1,ndet)
        write(6,'(''cdet='',20f9.6)') (cdet(i),i=1,ndet)

        read(5,*) (eigs(i),i=1,norb)
        write(6,'(''eigs(i)='',20f10.6)') (eigs(i),i=1,norb)

        do 5 iorb=1,norb
    5     read(5,*) (coef(ibasis,iorb),ibasis=1,nbasis)

        do 10 idet=1,ndet
          read(5,*) (iworbd(iel,idet),iel=1,nelec)
   10     write(6,'(90i4)') (iworbd(iel,idet),iel=1,nelec)

       elseif(index(mode,'genci').ne.0) then

c eps should be something like 4.d-6
        read(5,*) ndet,nelec,norb,eps,ncore,nalpha,nbeta,nact
        write(6,'(''ndet,nelec,norb,eps,nalpha,nbeta,nact='',3i4,1pd9.2,0p,9i3)') ndet,nelec,norb,eps,nalpha,nbeta,nact
        write(6,'(''ndet='',i4,'',ncore='',i2,'',nalpha='',i2,'',nbeta='',i2,'',nact='',i2)') ndet,ncore,nalpha,nbeta,nact
        if(ndet.gt.MDET) stop 'ndet>MDET'
        if(nelec.gt.MELEC) stop 'nelec>MELEC'
        if(norb.gt.MORB) stop 'norb>MORB'
        if(nelec.ne.2*ncore+nalpha+nbeta) stop 'nelec != 2*ncore+nalpha+nbeta'
        if(nact.gt.MACT) stop 'nact>MACT'

        read(5,*) (eigs(i),i=1,norb)
        write(6,'(''eigs(i)='',20f10.6)') (eigs(i),i=1,norb)

        write(6,'(''Reading in determinants and coefficients'')')
        write(fmt,'(''(a,A''i3,'',a,A''i3,'',a,f9.6)'')') nact,nact
        do 20 i=1, ndet
          read(5,*) alphain(i),ctmp,betain(i),ctmp2,cdet(i)
   20     write(6,fmt) 'alpha=',alphain(i),', beta=',betain(i),', cdet=',cdet(i)
        write(6,*)

c put genci determinants into CHAMP format
        norb_max=0
        do 60 idet=1,ndet
          do 30 j=1,ncore
            alpha(idet,j)=j
            beta(idet,j)=j
            iworbd(j,idet)=j
   30       iworbd(j+ncore+nalpha,idet)=j

          iel=ncore
          ctmp=alphain(idet)
          do 40 j=1,nact
            if(ctmp(j:j).eq.'1') then
              iel=iel+1
              if(iel.gt.MELEC) stop 'iel>MELEC'
              if(iel.gt.(nalpha+ncore)) stop 'iel>nalpha+ncore'
              if(ncore+j.gt.MORB) stop 'ncore+j>MORB'
              norb_max=max(ncore+j,norb_max)
              alpha(idet,iel)=ncore+j
              iworbd(iel,idet)=ncore+j
            endif
   40     continue

          iel=ncore
          ctmp=betain(idet)
          do 50 j=1,nact
            if(ctmp(j:j).eq.'1') then
              iel=iel+1
              if(iel+ncore+nalpha.gt.MELEC) stop 'iel+ncore+nalpha>MELEC'
              if(iel.gt.(nbeta+ncore)) stop 'iel>nbeta+ncore'
              if(ncore+j.gt.MORB) stop 'ncore+j>MORB'
              norb_max=max(ncore+j,norb_max)
              beta(idet,iel)=ncore+j
              iworbd(iel+ncore+nalpha,idet)=ncore+j
            endif
   50     continue
          write(6,'(''idet,iworbd(iel,idet)'',i3,2x,100i4)') idet,(iworbd(iel,idet),iel=1,nelec)
   60   continue

c It could happen that one or more of the lower orbs does not appear in any
c of the determinants, but for the moment we assume that they all appear.
        if(norb.ne.norb_max) then
          write(6,'(''norb,norb_max='',2i5)') norb,norb_max
          stop 'The value of norb read in does not match that deduced from the genci determinants'
        endif

        write(6,'(''Determinants converted to CHAMP format'')')
c       write(fmt2,'(''(''i3,''i3,a)'')') 2*ncore+nalpha+nbeta
        write(fmt2,'(''(''i3,''i3,x,''i3,''i3,a)'')') ncore+nalpha,ncore+nbeta
        do 70 i=1, ndet-1
   70     write(6,fmt2) (alpha(i,j),j=1,ncore+nalpha),(beta(i,j),j=1,ncore+nbeta)
        write(6,fmt2) (alpha(ndet,j),j=1,ncore+nalpha),(beta(ndet,j),j=1,ncore+nbeta),' (iworbd(j,idet),j=1,nelec)'
        write(6,*)

       else

        write(6,'(''mode should be either guga or genci'')')
        stop 'mode should be either guga or genci'

      endif

c Write out eigs just to see in practice how close the ones that should be degenerate are, since only those within eps can be put in a CSF.
      do 80 i=1,norb
   80   eigs_tmp(i)=eigs(i)
      call shell(eigs_tmp,norb)
      write(6,'(/,''Only determinants with orbitals whose eigs are within 0.1*eps of each other can be combined into CSFs'')')
      write(6,'(''This helps in choosing eps'')')
      write(6,'(''sorted eigs'',2000f13.9)') (eigs_tmp(i),i=1,norb)

c Write out the small eigenvalue diffs.
c In loop over 130 we are giving eigenvals. that match within 0.1d0*eps the same label.
c If the eigenvalues are close, but do not fit this criterion, write out warning.
      ndiffs=0
      do 90 i=2,norb
        if(eigs_tmp(i)-eigs_tmp(i-1).gt.0.1d0*eps .and. eigs_tmp(i)-eigs_tmp(i-1).lt.100*eps) then
          write(6,'(''Warning: eigs. close but not close enough to get same label, iorb,eigs='',2i6,2f13.9)')
     &    i-1,i,eigs_tmp(i-1),eigs_tmp(i)
          ndiffs=ndiffs+1
          eigs_tmp(ndiffs)=eigs_tmp(i)-eigs_tmp(i-1)
        endif
   90 continue
c     write(6,'(i4,'' Eigenvalue diffs are'')') ndiffs
c     write(6,'(''Eig diffs'',2000f10.6)') (eigs_tmp(i),i=1,ndiffs)
      call shell(eigs_tmp,ndiffs)
      if(ndiffs.ge.1) write(6,'(''Warning: Eig diffs'',2000f12.9)') (eigs_tmp(i),i=1,ndiffs)

c Sort eigenvalues and match eigenvalue sets to label dets. so that those dets. that
c can be in the same CSF have the same label.
c For the eigenvlues we use 0.1*eps rather than eps as the matching criterion.
      write(6,'(/,''Only determinants with orbitals whose eigs are within 0.1*eps of each other can be combined into CSFs'')')
      do 120 idet=1,ndet
        do 110 iel=1,nelec
  110     eigs_in_det(iel,idet)=eigs(iworbd(iel,idet))
        call shell(eigs_in_det(1,idet),nelec)
  120   write(6,'(''sorted eigs'',i4,2000f13.9)') idet,(eigs_in_det(iel,idet),iel=1,nelec)

      do 125 idet=2,ndet
  125   label_det(idet)=0

      nlabel=1
      label_det(1)=1
      do 150 idet=2,ndet
        do 140 jdet=1,idet-1
          do 130 iel=1,nelec
            if(abs(eigs_in_det(iel,idet)-eigs_in_det(iel,jdet)).gt.0.1d0*eps) goto 140
  130     continue
          label_det(idet)=label_det(jdet)
          goto 150
  140   continue
        nlabel=nlabel+1
        label_det(idet)=nlabel
  150 continue
      write(6,'(''nlabel,label_det(idet)='',i3,2x,100i3,(100i3))') nlabel,(label_det(idet),idet=1,ndet)

      do 160 i=1,ndet
        frac(1,i)=1
        iflag_det(i)=0
  160   iflag_label_det(i)=0

c Read CSFs from gamess2qmc and flag those determinants that are in CSFs with csf_coef >= cutoff_d2c
c Check that all dets within a CSF have the same label_det
c Flag the label_det's that have csf_coef(icsf) >= cutoff_d2c
c Check that all dets in a given CSF have the same label_det
      read(5,*) ncsf
      write(6,'(''No of CSFs read in='',i4)') ncsf
      if(ncsf.gt.MCSF) stop 'ncsf>MCSF'
      read(5,*) (csf_coef(icsf),icsf=1,ncsf)
      write(6,'(''Input csf_coef='',100f8.4)') (csf_coef(icsf),icsf=1,ncsf)
      read(5,*) (ndet_in_csf(icsf),icsf=1,ncsf)
      write(6,'(''Input  ndet_in_csf'',100i4)') (ndet_in_csf(icsf),icsf=1,ncsf)
      do 165 icsf=1,ncsf
        read(5,*) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
        read(5,*) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))
        do 162 idet_in_csf=1,ndet_in_csf(icsf)
  162     if(label_det(iwdet_in_csf(idet_in_csf,icsf)).ne.label_det(iwdet_in_csf(1,icsf)))
     &    stop 'Not all dets in CSF have the same label'
        if(abs(csf_coef(icsf)).ge.cutoff_d2c) then
          do 164 idet_in_csf=1,ndet_in_csf(icsf)
            iflag_det(iwdet_in_csf(idet_in_csf,icsf))=1
  164       iflag_label_det(label_det(iwdet_in_csf(idet_in_csf,icsf)))=1
        endif
  165 continue
c     write(6,'(''label_det(~icsf)='',100i8)') (label_det(iwdet_in_csf(1,icsf)),icsf=1,ncsf)
      write(6,'(''nlabel,iflag_det(idet)='',i3,2x,90i3)') nlabel,(iflag_det(idet),idet=1,ndet)
      write(6,'(''nlabel,iflag_label_det(ilabel)='',i3,2x,90i3)') nlabel,(iflag_label_det(ilabel),ilabel=1,nlabel)

c Write out csf_coefs just to see in practice how close the ones that should be identical are, since only those within eps can be put in a CSF.
      do i=1,ncsf
        csf_coef_tmp(i)=abs(csf_coef(i))
      enddo
      call shell(csf_coef_tmp,ncsf)
      write(6,'(/,''Only determinants in CSFs whose abs(csf_coef) are within eps*sqrt(abs(csf_coef(icsf))) of each other can be'',
     &'' combined into CSFs'')')
      write(6,'(''This helps in choosing eps'')')
      write(6,'(''sorted abs(csf_coef)'',i4,100000f10.6)') idet,(csf_coef_tmp(i),i=1,ncsf)

c Write out the small csf_coef diffs.
      ndiffs=0
      do i=2,ncsf
        if(csf_coef_tmp(i).gt.csf_coef_tmp(i-1) .and. csf_coef_tmp(i)-csf_coef_tmp(i-1).lt.100*eps) then
          ndiffs=ndiffs+1
          abs_csf_coef(ndiffs)=csf_coef_tmp(i)
          csf_coef_tmp(ndiffs)=csf_coef_tmp(i)-csf_coef_tmp(i-1)
        endif
      enddo
      write(6,'(i4,'' abs(csf_coef) diffs are'')') ndiffs
      write(6,'(/,''abs(csf_coef) diffs'',2000f10.6)') (csf_coef_tmp(i),i=1,ndiffs)
      write(6,'(''And the corresponding abs(csf_coef) are'')')
      write(6,'(''abs(csf_coef)      '',2000f10.6)') (abs_csf_coef(i),i=1,ndiffs)
      call shell(csf_coef_tmp,ndiffs)
c     write(6,'(''abs(csf_coef) diffs'',2000f10.6)') (csf_coef_tmp(i),i=1,ndiffs)
      ndiffs2=1
      do i=2,ndiffs
        if(csf_coef_tmp(i).gt.csf_coef_tmp(i-1)+1.d-9) then
          ndiffs2=ndiffs2+1
          abs_csf_coef(ndiffs2)=abs_csf_coef(i)
          csf_coef_tmp(ndiffs2)=csf_coef_tmp(i)
        endif
      enddo
      write(6,'(/,i4,'' abs(csf_coef) diffs are'')') ndiffs2
      write(6,'(''abs(csf_coef) diffs'',2000f10.6)') (csf_coef_tmp(i),i=1,ndiffs2)
c abs_csf_coef was not sorted when csf_coef_tmp was sorted, so the next write is not correct.
c     write(6,'(''abs(csf_coef)      '',2000f10.6)') (abs_csf_coef(i),i=1,ndiffs2)

c Make cdet_in_csf more accurate since they are printed out only to an accuracy of 10^-6
      do 167 icsf=1,ncsf
        do 167 idet_in_csf=1,ndet_in_csf(icsf)
          do 166 inum=1,50
            do 166 iden=1,20000
              if(abs(abs(cdet_in_csf(idet_in_csf,icsf))-sqrt(dfloat(inum)/dfloat(iden))).le.1.d-6) then
c               write(6,*) 'Match',abs(cdet_in_csf(idet_in_csf,icsf)),dfloat(inum)/dfloat(iden),inum,iden
                write(6,'(''Match: inum,iden'',2i6)') inum,iden
                cdet_in_csf(idet_in_csf,icsf)=sign(sqrt(dfloat(inum)/dfloat(iden)),cdet_in_csf(idet_in_csf,icsf))
                goto 167
              endif
  166   continue
        write(6,'(''Warning: failed to match cdet_in_csf(idet_in_csf,icsf)='',f8.6)') cdet_in_csf(idet_in_csf,icsf)
        stop 'failed to match up cdet_in_csf(idet_in_csf,icsf)'
  167 continue

c flag CSFs
c a) that have csf_coef(icsf).ge.cutoff_d2c by iflag_csf(icsf)=2
c b) that have dets that have been flagged by  iflag_csf(icsf)=1
      do 168 icsf=1,ncsf
  168   iflag_csf(icsf)=0
      do 171 icsf=1,ncsf
        if(abs(csf_coef(icsf)).ge.cutoff_d2c) then
          iflag_csf(icsf)=2
         else
c         do 170 idet_in_csf=1,ndet_in_csf(icsf)
c           write(6,'(''csf_coef(icsf),icsf,idet_in_csf,label_det(iwdet_in_csf(idet_in_csf,icsf))'',f8.4,9i5)')
c    &      csf_coef(icsf),icsf,idet_in_csf,label_det(iwdet_in_csf(idet_in_csf,icsf))
            do 170 ilabel=1,nlabel
c 170         if(label_det(iwdet_in_csf(idet_in_csf,icsf)).eq.ilabel .and. iflag_label_det(ilabel).eq.1) iflag_csf(icsf)=1
  170         if(label_det(iwdet_in_csf(1,icsf)).eq.ilabel .and. iflag_label_det(ilabel).eq.1) iflag_csf(icsf)=1
        endif
  171 continue

c Temp printout
      write(6,'(/,''1Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i3,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 1340 icsf=1,ncsf
        if(ndet_in_csf(icsf).le.0) write(6,'(''Warning: icsf,ndet_in_csf(icsf)='',9i6)') icsf,ndet_in_csf(icsf)
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 1340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c Zero out cdets since we will recompute them for the CSFs that are kept
c In order to do a stable sort (preserve original order when coefs. are equal)
c we subtract a tiny amount depending on original order
      do 172 idet=1,ndet
  172   cdet(idet)=-idet*1.d-12

c Recompute cdet including those with csf_coef>=cutoff_d2c (> cutoff_g2q) or iflag_csf(icsf)=1
      do 176 icsf=1,ncsf
        write(6,'(''icsf,iflag_csf(icsf),label_det(iwdet_in_csf(1,icsf)),csf_coef(icsf)'',3i5,f10.6)')
     &  icsf,iflag_csf(icsf),label_det(iwdet_in_csf(1,icsf)),csf_coef(icsf)
        if(abs(csf_coef(icsf)).ge.cutoff_d2c .or. iflag_csf(icsf).eq.1) then
          do 174 idet_in_csf=1,ndet_in_csf(icsf)
  174       cdet(iwdet_in_csf(idet_in_csf,icsf))=cdet(iwdet_in_csf(idet_in_csf,icsf))+csf_coef(icsf)*cdet_in_csf(idet_in_csf,icsf)
        endif
  176 continue
      write(6,'(''cdets before sorting'',100f8.4,(100f8.4))') (cdet(idet),idet=1,ndet)

c Temp printout
      write(6,'(/,''2Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i3,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 2340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 2340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c Count how many determinants have |cdet|>1.d-6, sort the cdets and the corresponding dets and iwdet_in_csf
c     call shell_abs_cdet(cdet,iworbd,ndet_in_csf,iwdet_in_csf,nelec,ndet,ncsf,MELEC,MDET,MDET_IN_CSF,MCSF)
      write(6,'(''label_det'',i4,99i8,(100i8))') (label_det(idet),idet=1,ndet)

c Temp printout
      write(6,'(/,''3Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i3,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 3340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 3340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c Flag the orbitals that are used
      if(iused_orbs.ne.0) then
        do 180 iorb=1,norb
  180   iflag_orb(iorb)=0

        do 185 idet=1,ndet
          do 185 iel=1,nelec
  185       iflag_orb(iworbd(iel,idet))=1
        write(6,'(''iflag_orb='',200i2)') (iflag_orb(iorb),iorb=1,norb)

c Create a map of the used orbitals that can be used to eliminate those that are not used
        map_orb(1)=1
        norb_used=1
        do 190 iorb=2,norb
          if(iflag_orb(iorb).eq.1) then
            map_orb(iorb)=map_orb(iorb-1)+1
            norb_used=norb_used+1
           else
            map_orb(iorb)=map_orb(iorb-1)
          endif
  190   continue
        write(6,'(i4,'' out of'',i5,'' orbitals are used'')') norb_used,norb
        write(6,'(''map_orb='',200i4)') (map_orb(iorb),iorb=1,norb)

c Eliminate the unused orbitals and update the corresponding eigs
        do 195 iorb=1,norb
          if(iflag_orb(iorb).eq.1) then
            do 192 ibasis=1,nbasis
  192         coef(ibasis,map_orb(iorb))=coef(ibasis,iorb)
            eigs(map_orb(iorb))=eigs(iorb)
          endif
  195   continue
        if(norb_used.ne.norb) then
          write(6,'(''norb_used is < norb, norb_used, norb='',9i5)') norb_used,norb
         else
          write(6,'(''norb_used = norb ='',9i5)') norb
        endif
        norb=norb_used

c Update iworbd
        do 198 idet=1,ndet
          do 198 iel=1,nelec
  198       iworbd(iel,idet)=map_orb(iworbd(iel,idet))
      endif

c Sort eigenvalues and match eigenvalue sets to label dets. so that those dets. that
c can be in the same CSF have the same label.
c For the eigenvlues we use 0.1*eps rather than eps as the matching criterion.
      do 220 idet=1,ndet
        do 210 iel=1,nelec
  210     eigs_in_det(iel,idet)=eigs(iworbd(iel,idet))
        call shell(eigs_in_det(1,idet),nelec)
  220   write(6,'(''sorted eigs'',i4,2000f13.9)') idet,(eigs_in_det(iel,idet),iel=1,nelec)

      do 225 idet=2,ndet
  225   label_det(idet)=0

      nlabel=1
      label_det(1)=1
      do 250 idet=2,ndet
        do 240 jdet=1,idet-1
          do 230 iel=1,nelec
            if(abs(eigs_in_det(iel,idet)-eigs_in_det(iel,jdet)).gt.0.1d0*eps) goto 240
  230     continue
          label_det(idet)=label_det(jdet)
          goto 250
  240   continue
        nlabel=nlabel+1
        label_det(idet)=nlabel
  250 continue
      write(6,'(''nlabel,label_det(idet)='',i3,2x,1000i3,(100i3))') nlabel,(label_det(idet),idet=1,ndet)


c Zero out frac and iflag_det for reuse
      do 270 i=1,ndet
        frac(1,i)=1
  270   iflag_det(i)=0

c Reform new CSFs from determinants
c Sometimes cdet can be very small even though csf_coef are not because different csf_coefs can add
c to produce a very small cdet.  So, while checking, check not only
c if(abs(cdet(j)-isign*times*cdet(i)).lt.eps) but rather
c if(abs(cdet(j)-isign*times*cdet(i)).lt.min(eps,0.2d0*abs(cdet(j))))
c iflag_det=1 marks the determinants that have been matched to an earlier determinant.
c So, the number of dets with iflag_det=0 is the number of determinantal free parameters.
c Later we will combine CSFs.
c     nedet=0
c     ncsf=0
c     do 300 i=1,ndet
c       if(iflag_det(i).eq.1) goto 300
c       ncsf=ncsf+1
c       if(ncsf.gt.MCSF) stop 'ncsf>MCSF'
c       ndet_in_csf(ncsf)=1
c       iwdet_in_csf(1,ncsf)=i
c       csf_coef(ncsf)=cdet(i)
c       cdet_in_csf(1,ncsf)=1
c       do 290 j=i+1,ndet
c         if(iflag_det(j).eq.1 .or. label_det(i).ne.label_det(j)) goto 290
c         do 280 itimes=2,8
c           times=0.5d0*itimes
c           do 280 isign=-1,1,2
c             if(abs(cdet(j)-isign*times*cdet(i)).lt.min(eps,0.2d0*abs(cdet(j)))) then
c               iflag_det(j)=1
c               nedet=nedet+1
c               iedet(1,nedet)=j
c               iedet(2,nedet)=i
c               frac(2,nedet)=isign*times
c               ndet_in_csf(ncsf)=ndet_in_csf(ncsf)+1
c               if(ndet_in_csf(ncsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(ncsf) > MDET_IN_CSF'
c               iwdet_in_csf(ndet_in_csf(ncsf),ncsf)=j
c               cdet_in_csf(ndet_in_csf(ncsf),ncsf)=isign*times
c             endif
c             if(abs(isign*times*cdet(j)-cdet(i)).lt.min(eps,0.2d0*abs(cdet(j))) .and. itimes.ne.2) then
c               iflag_det(j)=1
c               nedet=nedet+1
c               iedet(1,nedet)=j
c               iedet(2,nedet)=i
c               frac(2,nedet)=1.d0/(isign*times)
c               ndet_in_csf(ncsf)=ndet_in_csf(ncsf)+1
c               if(ndet_in_csf(ncsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(ncsf) > MDET_IN_CSF'
c               iwdet_in_csf(ndet_in_csf(ncsf),ncsf)=j
c               cdet_in_csf(ndet_in_csf(ncsf),ncsf)=1.d0/(isign*times)
c             endif
c 280     continue
c 290   continue
c 300 continue

c Make largest cdet_in_csf of each CSF be 1.
      do 300 icsf=1,ncsf
        cdet_in_csf_max=0
        do 280 idet_in_csf=1,ndet_in_csf(icsf)
          if(abs(cdet_in_csf(idet_in_csf,icsf)).gt.abs(cdet_in_csf_max)) then
            cdet_in_csf_max=cdet_in_csf(idet_in_csf,icsf)
            iwdet_in_csf_max=idet_in_csf
          endif
  280   continue
        do 290 idet_in_csf=1,ndet_in_csf(icsf)
          cdet_in_csf(idet_in_csf,icsf)=cdet_in_csf(idet_in_csf,icsf)/cdet_in_csf_max
  290   continue
  300   csf_coef(icsf)=csf_coef(icsf)*cdet_in_csf_max

c Create old version inputs (cdet's rather than csf_coef's)
c     write(6,'(/,''4Inputs for old version of CHAMP'')')
c     write(6,'(i3,'' nedet'')') nedet
c     if(nedet.gt.0) then
c       write(fmt,'(''(''i3,''(i5,i4),a)'')') nedet
c       write(6,fmt) ((iedet(k,i),k=1,2),i=1,nedet),' ((iedet(k,i),k=1,2),i=1,nedet)'
c       write(fmt,'(''(''i3,''(f6.1,f5.1),a)'')') nedet
c       write(6,fmt) ((frac(k,i),k=1,2),i=1,nedet),' ((frac(k,i),k=1,2),i=1,nedet)'
c      else
c       write(6,'(a)') '((iedet(k,i),k=1,2),i=1,nedet)'
c       write(6,'(a)') '((frac(k,i),k=1,2),i=1,nedet)'
c     endif

c     nparmd=0
c     do 330 i=2,ndet
c       if(iflag_det(i).eq.0) then
c         nparmd=nparmd+1
c         iwdet(nparmd)=i
c       endif
c 330 continue
c     write(6,'(i4,'' nparmd'')') nparmd
c     if(nparmd.gt.0) then
c       write(fmt,'(''(''i5,''i5,a)'')') nparmd
c       write(6,fmt) (iwdet(ipar),ipar=1,nparmd),' (iwdet(ipar),ipar=1,nparmd)'
c      else
c       write(6,'(a)') '(iwdet(ipar),ipar=1,nparmd)'
c     endif

c Temp printout
      write(6,'(/,''4Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i4,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 4340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 4340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'


c Sort the CSFs by the absolute value of csf_coef
      write(6,'(''before shell_abs_csf_coef'')')
      call shell_abs_csf_coef(csf_coef,ncsf,ndet_in_csf,iwdet_in_csf,cdet_in_csf,MDET_IN_CSF,MCSF,1.d-6)
      write(6,'(''after shell_abs_csf_coef'')')

c Temp printout
      write(6,'(/,''5Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i4,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 5340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
c       do 5710 idet_in_csf=1,ndet_in_csf(icsf)
c5710     if(abs(0.5d0*nint(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-3)
c    &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f9.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 5340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c Combine CSFs that have related coefs
c The following puts the icsf CSF into the jcsf CSF and eliminates the icsf CSF.
c It then moves csfs from jcp on down a step and goes back to the very beginning.  Dumb way, but it works.
c Note: it is probably better to have icsf be the smallest index (largest coef) rather
c than the largest index (smallest coef) of the two because then one ends up with
c smaller CSF coefs. and if they are small enough then we can drop them.
  350 ncsf_tmp=ncsf
c     do 400 icsf=2,ncsf_tmp
      do 400 icsf=1,ncsf_tmp-1
c       do 400 jcsf=1,icsf-1
        do 400 jcsf=icsf+1,ncsf_tmp
          if(label_det(iwdet_in_csf(1,icsf)).eq.label_det(iwdet_in_csf(1,jcsf))) then
            do 390 inum=2,20
            do 390 iden=2,10
              times=dfloat(inum)/dfloat(iden)
              do 390 isignj=-1,1,2
                if(abs(csf_coef(icsf)-isignj*times*csf_coef(jcsf)).lt.eps .and.
     &             abs(csf_coef(icsf)-isignj*times*csf_coef(jcsf)).lt.1000*eps*sqrt(abs(csf_coef(icsf)))) then
                  write(6,'(''inum,iden,icsf,jcsf,csf_coef(icsf),csf_coef(jcsf)='',2i3,2i4,2f6.3)')
     &            inum,iden,icsf,jcsf,csf_coef(icsf),csf_coef(jcsf)
                  if(ndet_in_csf(jcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(jcsf) > MDET_IN_CSF'
                  do 370 idet_in_csf=1,ndet_in_csf(jcsf)
                    iwdet_in_csf(ndet_in_csf(icsf)+idet_in_csf,icsf)=iwdet_in_csf(idet_in_csf,jcsf)
  370               cdet_in_csf(ndet_in_csf(icsf)+idet_in_csf,icsf)=cdet_in_csf(idet_in_csf,jcsf)*csf_coef(jcsf)/csf_coef(icsf)
                  ndet_in_csf(icsf)=ndet_in_csf(icsf)+ndet_in_csf(jcsf)

                  do 380 lcsf=jcsf+1,ncsf_tmp
                    csf_coef(lcsf-1)=csf_coef(lcsf)
                    ndet_in_csf(lcsf-1)=ndet_in_csf(lcsf)
                    do 380 idet_in_csf=1,ndet_in_csf(lcsf)
                      iwdet_in_csf(idet_in_csf,lcsf-1)=iwdet_in_csf(idet_in_csf,lcsf)
  380                 cdet_in_csf(idet_in_csf,lcsf-1)=cdet_in_csf(idet_in_csf,lcsf)
                  ncsf=ncsf-1
                  goto 350
                endif
  390       continue
          endif
  400 continue

c And do it again
  401 ncsf_tmp=ncsf
      do 405 icsf=1,ncsf_tmp-1
        do 405 jcsf=icsf+1,ncsf_tmp
          if(abs(csf_coef(icsf)-isignj*times*csf_coef(jcsf)).lt.eps .and.
     &       abs(csf_coef(icsf)-isignj*times*csf_coef(jcsf)).lt.1000*eps*sqrt(abs(csf_coef(icsf)))) then
            do 404 inum=2,20
            do 404 iden=2,10
              times=dfloat(inum)/dfloat(iden)
              do 404 isignj=-1,1,2
                if(abs(csf_coef(icsf)-isignj*times*csf_coef(jcsf)).lt.eps*sqrt(abs(csf_coef(icsf)))) then
                  write(6,'(''2nd round inum,iden,icsf,jcsf,csf_coef(icsf),csf_coef(jcsf)='',2i3,2i4,2f6.3)')
     &            inum,iden,icsf,jcsf,csf_coef(icsf),csf_coef(jcsf)
                  if(ndet_in_csf(jcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(jcsf) > MDET_IN_CSF'
                  do 402 idet_in_csf=1,ndet_in_csf(jcsf)
                    iwdet_in_csf(ndet_in_csf(icsf)+idet_in_csf,icsf)=iwdet_in_csf(idet_in_csf,jcsf)
  402               cdet_in_csf(ndet_in_csf(icsf)+idet_in_csf,icsf)=cdet_in_csf(idet_in_csf,jcsf)*csf_coef(jcsf)/csf_coef(icsf)
                  ndet_in_csf(icsf)=ndet_in_csf(icsf)+ndet_in_csf(jcsf)

                  do 403 lcsf=jcsf+1,ncsf_tmp
                    csf_coef(lcsf-1)=csf_coef(lcsf)
                    ndet_in_csf(lcsf-1)=ndet_in_csf(lcsf)
                    do 403 idet_in_csf=1,ndet_in_csf(lcsf)
                      iwdet_in_csf(idet_in_csf,lcsf-1)=iwdet_in_csf(idet_in_csf,lcsf)
  403                 cdet_in_csf(idet_in_csf,lcsf-1)=cdet_in_csf(idet_in_csf,lcsf)
                  ncsf=ncsf-1
                  goto 401
                endif
  404       continue
          endif
  405 continue

c Temp printout
      write(6,'(/,''6Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i4,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 6340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
c       do 6710 idet_in_csf=1,ndet_in_csf(icsf)
c6710     if(abs(0.5d0*nint(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-3)
c    &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f9.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 6340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c If the same det appears twice in the same CSF because of combining the CSF's above, then combine the 2 determinants
      do 420 icsf=1,ncsf
        do 420 idet_in_csf=2,ndet_in_csf(icsf)
          do 410 jdet_in_csf=1,idet_in_csf-1
            write(6,'(''iwdet_in_csf(jdet_in_csf,icsf),'',3i5,f10.6)')
     &      icsf,iwdet_in_csf(jdet_in_csf,icsf),iwdet_in_csf(idet_in_csf,icsf),cdet_in_csf(idet_in_csf,icsf)
            if(iwdet_in_csf(jdet_in_csf,icsf).eq.iwdet_in_csf(idet_in_csf,icsf) .and. cdet_in_csf(idet_in_csf,icsf).ne.0.d0) then
              cdet_in_csf(jdet_in_csf,icsf)=cdet_in_csf(jdet_in_csf,icsf)+cdet_in_csf(idet_in_csf,icsf)
              cdet_in_csf(idet_in_csf,icsf)=0.d0
              goto 420
            endif
  410   continue
  420 continue

c Make largest cdet_in_csf be 1
      do 424 icsf=1,ncsf
        cdet_in_csf_max=0
        do 422 idet_in_csf=1,ndet_in_csf(icsf)
          if(abs(cdet_in_csf(idet_in_csf,icsf)).gt.abs(cdet_in_csf_max)) then
            cdet_in_csf_max=cdet_in_csf(idet_in_csf,icsf)
            iwdet_in_csf_max=idet_in_csf
          endif
  422   continue
        do 423 idet_in_csf=1,ndet_in_csf(icsf)
  423     cdet_in_csf(idet_in_csf,icsf)=cdet_in_csf(idet_in_csf,icsf)/cdet_in_csf_max
  424     csf_coef(icsf)=csf_coef(icsf)*cdet_in_csf_max

c Temp printout
      write(6,'(/,''7Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i4,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 7340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
c       do 7710 idet_in_csf=1,ndet_in_csf(icsf)
c7710     if(abs(0.5d0*nint(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-3)
c    &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f9.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 7340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c Remove duplicate iwdet_in_csf's and reduce ndet_in_csf(icsf) accordingly
      do 430 icsf=1,ncsf
        ndet_in_csf_new=0
        do 425 idet_in_csf=1,ndet_in_csf(icsf)
          if(abs(cdet_in_csf(idet_in_csf,icsf)).gt.1.d-4) then
            ndet_in_csf_new=ndet_in_csf_new+1
            cdet_in_csf(ndet_in_csf_new,icsf)=cdet_in_csf(idet_in_csf,icsf)
            iwdet_in_csf(ndet_in_csf_new,icsf)=iwdet_in_csf(idet_in_csf,icsf)
          endif
  425   continue
        ndet_in_csf(icsf)=ndet_in_csf_new
  430 continue


c Create new version inputs
      write(6,'(/,''7aInputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i4,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 440 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
c       do 435 idet_in_csf=1,ndet_in_csf(icsf)
c 435   if(abs(cdet_in_csf(idet_in_csf,icsf)).le.1.d-4)
c    &  write(6,'(''Warning: cdet_in_csf(idet_in_csf,icsf)='',d12.4)') cdet_in_csf(idet_in_csf,icsf)
c       do 436 idet_in_csf=1,ndet_in_csf(icsf)
c 436     if(abs(0.5d0*nint(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-3)
c    &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f9.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
  440   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c When the same determinant appears in more than one CSF then there can be linear
c relations between the coefs. of the determinants.  The following detects these and
c puts the icsf CSF into the jcsf and kcsf CSFs and eliminates the icsf CSF.
c Note: it is probably better to have icsf be the smallest index (largest coef) rather
c than the largest index (smallest coef) of the three because then one ends up with
c smaller CSF coefs. and if they are small enough then we can drop them.
  450 ncsf_tmp=ncsf
c     do 500 icsf=2,ncsf_tmp
      do 500 icsf=1,ncsf_tmp-2
c       do 500 jcsf=1,icsf-1
        do 500 jcsf=icsf+1,ncsf_tmp-1
c         do 500 kcsf=1,jcsf-1
          do 500 kcsf=jcsf+1,ncsf_tmp
            if(label_det(iwdet_in_csf(1,icsf)).eq.label_det(iwdet_in_csf(1,jcsf)) .and.
     &         label_det(iwdet_in_csf(1,icsf)).eq.label_det(iwdet_in_csf(1,kcsf))) then
            do 490 isignj=-1,1,2
              do 490 isignk=-1,1,2
                if(abs(csf_coef(icsf)-isignj*csf_coef(jcsf)-isignk*csf_coef(kcsf)).lt.eps*sqrt(abs(csf_coef(icsf)))) then
                  write(6,'(''icsf,jcsf,kcsf,csf_coef(icsf),csf_coef(jcsf),csf_coef(kcsf)='',3i3,3f6.3)')
     &            icsf,jcsf,kcsf,csf_coef(icsf),csf_coef(jcsf),csf_coef(kcsf)
                  if(ndet_in_csf(jcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(jcsf) > MDET_IN_CSF'
                  if(ndet_in_csf(kcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(kcsf) > MDET_IN_CSF'
                  do 470 idet_in_csf=1,ndet_in_csf(icsf)
                    iwdet_in_csf(ndet_in_csf(jcsf)+idet_in_csf,jcsf)=iwdet_in_csf(idet_in_csf,icsf)
                    iwdet_in_csf(ndet_in_csf(kcsf)+idet_in_csf,kcsf)=iwdet_in_csf(idet_in_csf,icsf)
                    cdet_in_csf(ndet_in_csf(jcsf)+idet_in_csf,jcsf)=isignj*cdet_in_csf(idet_in_csf,icsf)
  470               cdet_in_csf(ndet_in_csf(kcsf)+idet_in_csf,kcsf)=isignk*cdet_in_csf(idet_in_csf,icsf)
                  ndet_in_csf(jcsf)=ndet_in_csf(jcsf)+ndet_in_csf(icsf)
                  ndet_in_csf(kcsf)=ndet_in_csf(kcsf)+ndet_in_csf(icsf)
                  do 480 lcsf=icsf+1,ncsf_tmp
                    csf_coef(lcsf-1)=csf_coef(lcsf)
                    ndet_in_csf(lcsf-1)=ndet_in_csf(lcsf)
                    do 480 idet_in_csf=1,ndet_in_csf(lcsf)
                      iwdet_in_csf(idet_in_csf,lcsf-1)=iwdet_in_csf(idet_in_csf,lcsf)
  480                 cdet_in_csf(idet_in_csf,lcsf-1)=cdet_in_csf(idet_in_csf,lcsf)
                  ncsf=ncsf-1
                  goto 450
                endif
  490       continue
            endif
  500 continue

      write(6,'(/,''If any of the foll. numbers look rational then there was a failure in recognizing a relationship'')')
      if(ncsf.le.100) then
        do 503 icsf=1,ncsf-1
          do 503 jcsf=icsf+1,ncsf
  503       write(6,'(''icsf,jcsf,ratio='',2i4,9f12.6)')
     &      icsf,jcsf,(csf_coef(icsf)/csf_coef(jcsf))**2,(csf_coef(jcsf)/csf_coef(icsf))**2
      endif

c Sort the CSFs by the absolute value of csf_coef
      write(6,'(''before shell_abs_csf_coef'')')
      call shell_abs_csf_coef(csf_coef,ncsf,ndet_in_csf,iwdet_in_csf,cdet_in_csf,MDET_IN_CSF,MCSF,cutoff_d2c)
      write(6,'(''after shell_abs_csf_coef'')')

c Keep only CSF's with csf_coef larger than cutoff_d2c and corresponding dets
c This is not done quite correctly because of the normalization of the CSFs
c     ncsf_new=0
c     ndet_new=0
c     do 510 icsf=1,ncsf
c       csf_norm=0
c       do 505 idet_in_csf=1,ndet_in_csf(icsf)
c 505     csf_norm=csf_norm+cdet_in_csf(idet_in_csf,icsf)**2
c       csf_norm=sqrt(csf_norm)
c       if(abs(csf_norm*csf_coef(icsf)).ge.cutoff_d2c) then
c         ncsf_new=icsf
c         do 507 idet_in_csf=1,ndet_in_csf(icsf)
c 507       ndet_new=max(ndet_new,iwdet_in_csf(idet_in_csf,icsf))
c       endif
c 510 continue
c     write(6,'(''ncsf,ncsf_new,ndet,ndet_new='',9i5)') ncsf,ncsf_new,ndet,ndet_new
c     ncsf=ncsf_new
c     ndet=ndet_new

      if(csf_coef(1).lt.0) then
        do 520 icsf=1,ncsf
  520     csf_coef(icsf)=-csf_coef(icsf)
      endif

c Zero out cdets since we will recompute them for the CSFs that are kept
c In order to do a stable sort (preserve original order when coefs. are equal)
c we subtract a tiny amount depending on original order
      do 572 idet=1,ndet
  572   cdet(idet)=-idet*1.d-12

ccRecompute cdet including those with csf_coef>=cutoff_d2c (> cutoff_g2q) or iflag_csf(icsf)=1
c Recompute cdet
      do 576 icsf=1,ncsf
        write(6,'(''icsf,iflag_csf(icsf),label_det(iwdet_in_csf(1,icsf)),csf_coef(icsf)'',3i5,f10.6)')
     &  icsf,iflag_csf(icsf),label_det(iwdet_in_csf(1,icsf)),csf_coef(icsf)
c       if(abs(csf_coef(icsf)).ge.cutoff_d2c .or. iflag_csf(icsf).eq.1) then
          do 574 idet_in_csf=1,ndet_in_csf(icsf)
  574       cdet(iwdet_in_csf(idet_in_csf,icsf))=cdet(iwdet_in_csf(idet_in_csf,icsf))+csf_coef(icsf)*cdet_in_csf(idet_in_csf,icsf)
c       endif
  576 continue
      write(6,'(''cdets before sorting'',100f8.4,(100f8.4))') (cdet(idet),idet=1,ndet)

c Count how many determinants have |cdet|>1.d-8 (done just to reduce print out)
c Sort the cdets and the corresponding dets and iwdet_in_csf
      call shell_abs_cdet(cdet,iworbd,ndet_in_csf,iwdet_in_csf,nelec,ndet,ncsf,MELEC,MDET,MDET_IN_CSF,MCSF)
      write(6,'(''cdets after sorting'',100f8.4,/,(100f8.4))') (cdet(idet),idet=1,ndet)
      write(6,'(''label_det'',i4,99i8,/,(100i8))') (label_det(idet),idet=1,ndet)

c Determine ndet by checking which appear in the kept csfs.
c Some csfs may combine to give very small cdet's and those could be
c eliminated, but at present we just print informative warning.
      ndet=0
      do 578 icsf=1,ncsf
        do 578 idet_in_csf=1,ndet_in_csf(icsf)
c         if(abs(cdet(iwdet_in_csf(idet_in_csf,icsf))).lt.1.d-6)
c    &    write(6,'(''Warning: det'',i5,'' could probably be eliminated but it just costs a bit of extra time to keep it'')')
c    &    iwdet_in_csf(idet_in_csf,icsf)
  578     ndet=max(ndet,iwdet_in_csf(idet_in_csf,icsf))
      write(6,'(''used ndet='',i5)') ndet

c Temp printout
      write(6,'(/,''8Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i4,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 8340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
c       do 8710 idet_in_csf=1,ndet_in_csf(icsf)
c8710     if(abs(0.5d0*nint(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-3)
c    &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f9.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 8340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

      if(iused_orbs.ne.0) then
c Flag the orbitals that are used
        do 580 iorb=1,norb
  580   iflag_orb(iorb)=0

        do 585 idet=1,ndet
          do 585 iel=1,nelec
  585       iflag_orb(iworbd(iel,idet))=1
        write(6,'(''iflag_orb='',200i2)') (iflag_orb(iorb),iorb=1,norb)

c Create a map of the used orbitals that can be used to eliminate those that are not used
        map_orb(1)=1
        norb_used=1
        do 590 iorb=2,norb
          if(iflag_orb(iorb).eq.1) then
            map_orb(iorb)=map_orb(iorb-1)+1
            norb_used=norb_used+1
           else
            map_orb(iorb)=map_orb(iorb-1)
          endif
  590   continue
        write(6,'(i4,'' out of'',i5,'' orbitals are used'')') norb_used,norb
        write(6,'(''map_orb='',200i4)') (map_orb(iorb),iorb=1,norb)

c Eliminate the unused orbitals and update the corresponding eigs
        do 595 iorb=1,norb
          if(iflag_orb(iorb).eq.1) then
            do 592 ibasis=1,nbasis
  592         coef(ibasis,map_orb(iorb))=coef(ibasis,iorb)
            eigs(map_orb(iorb))=eigs(iorb)
          endif
  595 continue
        norb=norb_used

c Update iworbd
        do 598 idet=1,ndet
          do 598 iel=1,nelec
  598       iworbd(iel,idet)=map_orb(iworbd(iel,idet))
      endif


c Sort eigenvalues and match eigenvalue sets to label dets. so that those dets. that
c can be in the same CSF have the same label.
c For the eigenvlues we use 0.1*eps rather than eps as the matching criterion.
      do 620 idet=1,ndet
        do 610 iel=1,nelec
  610     eigs_in_det(iel,idet)=eigs(iworbd(iel,idet))
        call shell(eigs_in_det(1,idet),nelec)
  620   write(6,'(''sorted eigs'',i4,2000f13.9)') idet,(eigs_in_det(iel,idet),iel=1,nelec)

      do 625 idet=2,ndet
  625   label_det(idet)=0

      nlabel=1
      label_det(1)=1
      do 650 idet=2,ndet
        do 640 jdet=1,idet-1
          do 630 iel=1,nelec
            if(abs(eigs_in_det(iel,idet)-eigs_in_det(iel,jdet)).gt.0.1d0*eps) goto 640
  630     continue
          label_det(idet)=label_det(jdet)
          goto 650
  640   continue
        nlabel=nlabel+1
        label_det(idet)=nlabel
  650 continue
      write(6,'(''nlabel,label_det(idet)='',i3,2x,1000i3,(100i3))') nlabel,(label_det(idet),idet=1,ndet)

c Temp printout
      write(6,'(/,''9Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i4,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 9340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
c       do 9710 idet_in_csf=1,ndet_in_csf(icsf)
c9710     if(abs(0.5d0*nint(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-3)
c    &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f9.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
 9340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c Make cdet_in_csf be integers.
c Change csf_coef so that they are appropriate not for either the existing cdet_in_csf
c or the new cdet_in_csf but for normalized csfs.
c First calculate normalization and adjust csf_coef to correspond to that.
      do 670 icsf=1,ncsf
        csf_norm=0
        do 660 idet_in_csf=1,ndet_in_csf(icsf)
  660     csf_norm=csf_norm+cdet_in_csf(idet_in_csf,icsf)**2
        csf_norm=sqrt(csf_norm)
        csf_coef(icsf)=csf_coef(icsf)*csf_norm
        do 670 idet_in_csf=1,ndet_in_csf(icsf)
  670     cdet_in_csf(idet_in_csf,icsf)=cdet_in_csf(idet_in_csf,icsf)/csf_norm

c Now make cdet_in_csf integers to make the input file smaller.  The qmc program will
c turn these back into normalized csfs.
      do 690 icsf=1,ncsf
        cdet_in_csf_min=9.d99
        do 680 idet_in_csf=1,ndet_in_csf(icsf)
  680     cdet_in_csf_min=min(cdet_in_csf_min,abs(cdet_in_csf(idet_in_csf,icsf)))
        do 685 idet_in_csf=1,ndet_in_csf(icsf)
  685     cdet_in_csf(idet_in_csf,icsf)=cdet_in_csf(idet_in_csf,icsf)/cdet_in_csf_min

  686   continue

        do 690 idet_in_csf=1,ndet_in_csf(icsf)
          do 688 iden=2,50
            do 688 inum=iden+1,50
              if(dfloat(inum)/dfloat(iden)-dfloat(inum/iden). gt. 1.d-9 .and.
     &          abs(abs(cdet_in_csf(idet_in_csf,icsf))-(dfloat(inum)/dfloat(iden))).lt.1.d-3) then
                write(6,'(''inum,iden,icsf='',2i3,i4)') inum,iden,icsf
                do 687 jdet_in_csf=1,ndet_in_csf(icsf)
  687             cdet_in_csf(jdet_in_csf,icsf)=cdet_in_csf(jdet_in_csf,icsf)*iden
                goto 686
              endif
  688   continue
  690 continue
      write(6,'(''ndet='',i5)') ndet

c Calculate total norm of the csf's kept.
      csf_sum=0
      do 695 icsf=1,ncsf
  695   csf_sum=csf_sum+csf_coef(icsf)**2
      csf_sum=sqrt(csf_sum)
      write(6,'(''csf_sum='',f8.6)') csf_sum

c Check if there are close csf coefs.  At this point they are in order of decreasing magnitude, so upper limit of iden is inum.
      do 699 icsf=1,ncsf-1
        do 699 jcsf=icsf+1,ncsf
          if(label_det(iwdet_in_csf(1,icsf)).eq.label_det(iwdet_in_csf(1,jcsf))) then
            do 698 inum=1,20
            do 698 iden=1,10
              times=dfloat(inum)/dfloat(iden)
              do 698 isignj=-1,1,2
                if(abs(csf_coef(icsf)-isignj*times*csf_coef(jcsf)).lt.100*eps*sqrt(abs(csf_coef(icsf)))) then
                  write(6,'(''Warning close coefs.: inum,iden,icsf,jcsf,csf_coef(icsf),csf_coef(jcsf)='',2i3,2i5,2f13.9)')
     &            inum,iden,icsf,jcsf,csf_coef(icsf),csf_coef(jcsf)
                  goto 699
                endif
  698       continue
          endif
  699 continue

c Create new version inputs
c While doing that, check if any cdet_in_csf is not an integer or half-integer.
c I need to do the checking completely differently using orbital energies or occupations.
      write(6,'(/,''Inputs for new version of CHAMP with linear relations recognized'')')

c First write the easy lines that are to be swapped into the input file
c     write(6,'(a1,''ncsf='',i2,'' ndet='',i3,'' norb='',i3,'' csf_sum='',f8.6,a1,a)')
c    &"'",ncsf,ndet,norb,csf_sum,"'",'  title'
      write(6,'(a1,''ncsf='',i3,'' ndet='',i4,'' norb='',i3,'' csf_sum='',f8.6,'' cutoff_g2q='',f6.4,'' cutoff_d2c='',f5.3,
     &'' eps='',es7.1,a1,a)')
     &"'",ncsf,ndet,norb,csf_sum,cutoff_g2q,cutoff_d2c,eps,"'",'  title'
      write(6,'(i4,2i4,t42,a)') ndet,nbasis,norb, 'ndet,nbasis,norb'
      nparm=24+ncsf-1
      write(6,'(''1000 '',i3,'' -1 1 5 1000 21101 1 NDATA,NPARM,icusp,icusp2,NSIG,NCALLS,iopt,ipr'')') nparm
      write(6,'(''0  4  5  15  0 '',i3,'' 0 0  nparml,nparma,nparmb,nparmc,nparmf,nparmcsf,nparms,nparmg'')') ncsf-1

c Write the line for making 0Info
      write(6,'(''MCSCF->CI natural orbs with cutoff_g2q='',f5.3,'', cutoff_d2c='',f5.3)') cutoff_g2q,cutoff_d2c

      write(6,'($,8i2)') (icsf,icsf=2,min(ncsf,9))
      write(6,'($,99i3)') (icsf,icsf=10,min(ncsf,99))
      write(6,'($,99i4)') (icsf,icsf=100,min(ncsf,999))
      write(6,'($,99i5)') (icsf,icsf=1000,min(ncsf,9999))
      write(6,'('' (iwcsf(iparm),iparm=1,nparmcsf)'')')

c Write out orbs if we are keeping only the orbitals used in the CSFs (iused_orbs=1).  If we are not doing orbital optimization by
c rotations, then that is all we need, but if we do orbital optimization by rotations then we need all orbitals.
      if(iused_orbs.eq.1) then
        do 706 iorb=1,norb
          write(fmt,'(''(1p,'',i3''d16.8,a)'')') nbasis
          if(iorb.eq.1) then
            write(6,fmt) (coef(ibasis,iorb),ibasis=1,nbasis),' ((coef(ibasis,iorb),ibasis=1,nbasis),iorb=1,norb)'
           else
            write(6,fmt) (coef(ibasis,iorb),ibasis=1,nbasis)
          endif
  706   continue
      endif

      ndn=nelec-nup
      do 707 idet=1,ndet
        if(nup-1.ge.1 .and. ndn.ge.1) then
          write(fmt,'(''(i3,'',i2,''i4,3x,'',i2,''i4,i5,a)'')') nup-1,ndn
         elseif(nup-1.ge.1 .and. ndn.eq.0) then
          write(fmt,'(''(i3,'',i2,''i4,3x,''i5,a)'')') nup-1
         elseif(nup-1.eq.0 .and. ndn.ge.1) then
          write(fmt,'(''(i3,3x,'',i2,''i4,i5,a)'')') ndn
        endif
        if(idet.eq.ndet) then
          write(6,fmt) (iworbd(iel,idet),iel=1,nelec), label_det(idet), ' (iworbd(iel,idet),iel=1,nelec), label_det(idet)'
         else
          write(6,fmt) (iworbd(iel,idet),iel=1,nelec), label_det(idet)
        endif
  707 continue

      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i4,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i4,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 720 icsf=1,ncsf
        write(fmt,'(''(''i3,''i4,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        do 710 idet_in_csf=1,ndet_in_csf(icsf)
          if(iwdet_in_csf(idet_in_csf,icsf).gt.ndet) then
            write(6,'(''icsf='',i4,'' iwdet_in_csf,ndet='',80i4)') icsf,iwdet_in_csf(idet_in_csf,icsf),ndet
            write(6,'(''cdet_in_csf(idet_in_csf,icsf)='',d12.4)') cdet_in_csf(idet_in_csf,icsf)
            write(6,'(''cdet(iwdet_in_csf(idet_in_csf,icsf))='',d12.4)') cdet(iwdet_in_csf(idet_in_csf,icsf))
            stop 'iwdet_in_csf(idet_in_csf,icsf) > ndet'
          endif
c         if(abs(cdet(iwdet_in_csf(idet_in_csf,icsf))).lt.1.d-6)
c    &    write(6,'(''Warning: det'',i5,'' could probably be eliminated but it just costs a bit of extra time to keep it'')')
c    &    iwdet_in_csf(idet_in_csf,icsf)
c 710     if(abs(0.5d0*nint(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-3)
c    &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f9.6)') cdet_in_csf(idet_in_csf,icsf)
  710   continue
        write(fmt,'(''(''i3,''i4,a)'')') ndet_in_csf(icsf)
  720   write(6,fmt) (nint(cdet_in_csf(idet_in_csf,icsf)),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

      stop
      end
c-----------------------------------------------------------------------
      subroutine shell(d,n)
      implicit real*8(a-h,o-z)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:::  shell-metzger sort in ascending order.          ...Cyrus 1979  :::
c:::  modified slightly for readibility.          ...Cyrus 7 dec 83  :::
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      dimension d(n)
      lognb2=int(dlog(dfloat(n))/dlog(2.d0)+1.d-14)
      m=n
      do 20 nn=1,lognb2
        m=m/2
        k=n-m
        do 20 j=1,k
          do 10 i=j,1,-m
            l=i+m
            if(d(l).ge.d(i)) goto 20
            t=d(i)
            d(i)=d(l)
            d(l)=t
   10       continue
   20     continue
      return
      end
c-----------------------------------------------------------------------
      subroutine shell_abs_cdet(cdet,iworbd,ndet_in_csf,iwdet_in_csf,nelec,ndet,ncsf,MELEC,MDET,MDET_IN_CSF,MCSF)
      implicit real*8(a-h,o-z)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:::  shell-metzger sort in descending order.         ...Cyrus 1979  :::
c:::  modified slightly for readibility.          ...Cyrus 7 dec 83  :::
c:::  modified for cdet, iworbd sorting.         ...Cyrus 30 apr 06  :::
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      dimension cdet(ndet),iworbd(MELEC,MDET),ndet_in_csf(MCSF),iwdet_in_csf(MDET_IN_CSF,MCSF)
      lognb2=int(dlog(dfloat(ndet))/dlog(2.d0)+1.d-14)
      m=ndet
      do 20 nn=1,lognb2
        m=m/2
        k=ndet-m
        do 20 j=1,k
          do 10 i=j,1,-m
            l=i+m
            if(abs(cdet(l)).le.abs(cdet(i))) goto 20
            t=cdet(i)
            cdet(i)=cdet(l)
            cdet(l)=t
            do 5 ielec=1,nelec
              it=iworbd(ielec,i)
              iworbd(ielec,i)=iworbd(ielec,l)
    5         iworbd(ielec,l)=it
            do 7 icsf=1,ncsf
              do 7 idet_in_csf=1,ndet_in_csf(icsf)
                if(iwdet_in_csf(idet_in_csf,icsf).eq.i) then
                  iwdet_in_csf(idet_in_csf,icsf)=l
                  write(6,'(''switched i to l'',9i5)') i,l
                 elseif(iwdet_in_csf(idet_in_csf,icsf).eq.l) then
                  iwdet_in_csf(idet_in_csf,icsf)=i
                  write(6,'(''switched l to i'',9i5)') l,i
                endif
c               it=iwdet_in_csf(idet_in_csf,i)
c               iwdet_in_csf(idet_in_csf,i)=iwdet_in_csf(idet_in_csf,l)
c   7           iwdet_in_csf(idet_in_csf,l)=it
    7       continue
   10     continue
   20 continue

c Calculate number of determinants with abs(cdet) >= 1.d-8
c This is done just to reduce the printout.  Later we will calculate ndet by checking
c which dets appear in the csfs that are kept.
      ndet_new=0
      do 30 idet=1,ndet
c       if(abs(cdet(idet)).le.1.d-5) write(6,'(''Warning small cdet('',i3,'' )='',d12.4)') idet,cdet(idet)
   30   if(abs(cdet(idet)).ge.1.d-8) ndet_new=ndet_new+1
      write(6,'(''ndet,ndet with abs coefs. > 1.d-8'',9i5)') ndet,ndet_new
c     write(6,'(''cdet='',100f8.4)') (cdet(idet),idet=1,ndet_new)
      ndet=ndet_new

      return
      end
c-----------------------------------------------------------------------
      subroutine shell_abs_csf_coef(csf_coef,ncsf,ndet_in_csf,iwdet_in_csf,cdet_in_csf,MDET_IN_CSF,MCSF,cutoff)
      implicit real*8(a-h,o-z)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:::  shell-metzger sort in descending order.         ...Cyrus 1979  :::
c:::  modified slightly for readibility.          ...Cyrus 7 dec 83  :::
c:::  modified for csf_coef, iwdet_in_csft,cdet_in_csf sorting       :::
c:::                                             ...Cyrus 30 apr 06  :::
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      dimension csf_coef(MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_IN_CSF,MCSF),cdet_in_csf(MDET_IN_CSF,MCSF)
     &,iwdet_in_csft(MDET_IN_CSF),cdet_in_csft(MDET_IN_CSF),csf_norm(MCSF)

c Find normalization of cdet_in_csf's
      do 3 icsf=1,ncsf
        csf_norm(icsf)=0
        do 2 idet_in_csf=1,ndet_in_csf(icsf)
    2     csf_norm(icsf)=csf_norm(icsf)+cdet_in_csf(idet_in_csf,icsf)**2
    3   csf_norm(icsf)=sqrt(csf_norm(icsf))

      lognb2=int(dlog(dfloat(ncsf))/dlog(2.d0)+1.d-14)
      m=ncsf
      do 20 nn=1,lognb2
        m=m/2
        k=ncsf-m
        do 20 j=1,k
          do 10 i=j,1,-m
            l=i+m
            if(abs(csf_coef(l)*csf_norm(l)).le.abs(csf_coef(i)*csf_norm(i))) goto 20
            do 5 idet_in_csf=1,ndet_in_csf(i)
              iwdet_in_csft(idet_in_csf)=iwdet_in_csf(idet_in_csf,i)
    5         cdet_in_csft(idet_in_csf)=cdet_in_csf(idet_in_csf,i)
            do 6 idet_in_csf=1,ndet_in_csf(l)
              iwdet_in_csf(idet_in_csf,i)=iwdet_in_csf(idet_in_csf,l)
    6         cdet_in_csf(idet_in_csf,i)=cdet_in_csf(idet_in_csf,l)
            do 7 idet_in_csf=1,ndet_in_csf(i)
              iwdet_in_csf(idet_in_csf,l)=iwdet_in_csft(idet_in_csf)
    7         cdet_in_csf(idet_in_csf,l)=cdet_in_csft(idet_in_csf)
            it=ndet_in_csf(i)
            ndet_in_csf(i)=ndet_in_csf(l)
            ndet_in_csf(l)=it
            t=csf_coef(i)
            csf_coef(i)=csf_coef(l)
            csf_coef(l)=t
            t=csf_norm(i)
            csf_norm(i)=csf_norm(l)
            csf_norm(l)=t
   10       continue
   20     continue

c Calculate number of CSFs with abs(csf_coef) >= cutoff
      ncsf_new=0
      do 30 icsf=1,ncsf
   30   if(abs(csf_coef(icsf)*csf_norm(icsf)).ge.cutoff) ncsf_new=ncsf_new+1
      write(6,'(''ncsf,ncsf_new='',9i5)') ncsf,ncsf_new
      write(6,'(''csf_coef='',100f8.4)') (csf_coef(icsf),icsf=1,ncsf_new)
      ncsf=ncsf_new

      return
      end

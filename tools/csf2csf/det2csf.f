      program det2csf
c Written by Cyrus Umrigar
c Modified by John Lawson to deal with GENCI output.

c Construct CSFs and CSF coefs. from determinantal coefs.
c Reads in determinantal coefs., cdet, and checks to see which pairs are
c related by certain multiplicative factors or which triplets are
c related with multiplicative factors of 1.
c It is not guaranteed to find all relationships, but it finds most of them.
c Also, it could find a relationship because of a chance equality of the cdet's
c that does not really exist.  I reduce the chance of that by using orbital
c eigenvalues or occupation numbers to label determinants that could enter in
c the same CSF.
c Output written in fashion to be useful to serve as input to both old and
c new versions of optimization program, but presently only the one for the
c new version is written out after checking the triplets.
c The input can come either from the output of a GUGA or a GENCI run using
c GAMESS.  If it is a GENCI run, then gamess2qmc does not write out the
c determinants in CHAMP format, so this program does that too.

      implicit real*8(a-h,o-z)
      character*80 fmt,fmt2,mode
      parameter(MDET=1000,MCSF=500,MDET_IN_CSF=40,MELEC=100,MORB=100,MBASIS=500,MLABEL=MDET,MACT=100)
      dimension cdet(MDET),iwdet(MDET),iedet(2,MDET),frac(2,MDET)
     &,iflag_det(MDET),iflag_label_det(MLABEL),iflag_csf(MCSF), csf_coef(MCSF),ndet_in_csf(MCSF)
     &,iwdet_in_csf(MDET_IN_CSF,MCSF),cdet_in_csf(MDET_IN_CSF,MCSF)
     &,iworbd(MELEC,MDET),eigs(MORB),eigs_in_det(MELEC,MDET),label_det(MDET)
     &,coef(MBASIS,MORB),iflag_orb(MORB),map_orb(MORB)
      character(MORB) alphain(MDET),betain(MDET),ctmp,ctmp2
      integer alpha(MDET,MACT),beta(MDET,MACT)

      read(5,*) mode
      write(6,'(''Reading '',a,'' mode input'')') mode
      read(5,*) cutoff_g2q,cutoff_d2c
      write(6,'(''Cutoffs in gamess2qmc and dets2csf are'',2f8.4)') cutoff_g2q,cutoff_d2c
      if(cutoff_d2c.lt.cutoff_g2q) stop 'cutoff_d2c must be >= cutoff_g2q'

      if(index(mode,'guga').ne.0) then

c eps should be something like 4.d-6
        read(5,*) ndet,nelec,nup,norb,nbasis,eps
        write(6,'(''ndet,nelec,nup,norb,nbasis,eps='',5i4,1pd9.2)') ndet,nelec,nup,norb,nbasis,eps
        if(ndet.gt.MDET) stop 'ndet>MDET'
        if(nelec.gt.MELEC) stop 'nelec>MELEC'
        if(norb.gt.MORB) stop 'norb>MORB'
        if(nbasis.gt.MBASIS) stop 'nbasis>MBASIS'

        read(5,*) (cdet(i),i=1,ndet)
        write(6,'(''cdet='',20f9.6)') (cdet(i),i=1,ndet)

        read(5,*) (eigs(i),i=1,norb)
        write(6,'(''eigs(i)='',20f9.6)') (eigs(i),i=1,norb)

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


c Sort eigenvalues and match eigenvalue sets to label dets. so that those dets. that
c can be in the same CSF have the same label.
      do 120 idet=1,ndet
        do 110 iel=1,nelec
  110     eigs_in_det(iel,idet)=eigs(iworbd(iel,idet))
        call shell(eigs_in_det(1,idet),nelec)
  120   write(6,'(''sorted eigs'',i3,90f10.6)') idet,(eigs_in_det(iel,idet),iel=1,nelec)

      do 125 idet=2,ndet
  125   label_det(idet)=0
 
      nlabel=1
      label_det(1)=1
      do 150 idet=2,ndet
        do 140 jdet=1,idet-1
          do 130 iel=1,nelec
            if(abs(eigs_in_det(iel,idet)-eigs_in_det(iel,jdet)).gt.eps) goto 140
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
      read(5,*) (csf_coef(icsf),icsf=1,ncsf)
      write(6,'(''Input csf_coef='',100f8.4)') (csf_coef(icsf),icsf=1,ncsf)
      read(5,*) (ndet_in_csf(icsf),icsf=1,ncsf)
      write(6,'(''Input  ndet_in_csf'',100i4)') (ndet_in_csf(icsf),icsf=1,ncsf)
      do 166 icsf=1,ncsf
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
  166 continue
c     write(6,'(''label_det(~icsf)='',100i8)') (label_det(iwdet_in_csf(1,icsf)),icsf=1,ncsf)
      write(6,'(''nlabel,iflag_det(idet)='',i3,2x,90i3)') nlabel,(iflag_det(idet),idet=1,ndet)
      write(6,'(''nlabel,iflag_label_det(ilabel)='',i3,2x,90i3)') nlabel,(iflag_label_det(ilabel),ilabel=1,nlabel)

c flag CSFs
c a) that have csf_coef(icsf).ge.cutoff_d2c by 2
c b) that have dets that have been flagged by 1
      do 168 icsf=1,ncsf
  168   iflag_csf(i)=0
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

c Count how many determinants have |cdet|>1.d-6, sort the cdets and the corresponding dets and iwdet_in_csf
c     call shell_abs_cdet(cdet,iworbd,iwdet_in_csf,nelec,ndet,MELEC,MDET)
      call shell_abs_cdet(cdet,iworbd,ndet_in_csf,iwdet_in_csf,nelec,ndet,ncsf,MELEC,MDET,MDET_IN_CSF,MCSF)
      write(6,'(''label_det'',i4,99i8,(100i8))') (label_det(idet),idet=1,ndet)

c Flag the orbitals that are used
      do 180 iorb=1,norb
!JT  180   iflag_orb(iorb)=0
  180   iflag_orb(iorb)=1 !JT: I need all the orbitals

      do 185 idet=1,ndet
        do 185 iel=1,nelec
  185     iflag_orb(iworbd(iel,idet))=1
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
  190 continue
      write(6,'(i4,'' out of'',i5,'' orbitals are used'')') norb_used,norb
      write(6,'(''map_orb='',200i4)') (map_orb(iorb),iorb=1,norb)

c Eliminate the unused orbitals and update the corresponding eigs
      do 195 iorb=1,norb
        if(iflag_orb(iorb).eq.1) then
          do 192 ibasis=1,nbasis
  192       coef(ibasis,map_orb(iorb))=coef(ibasis,iorb)
          eigs(map_orb(iorb))=eigs(iorb)
        endif
  195 continue
      norb=norb_used

c Update iworbd
      do 198 idet=1,ndet
        do 198 iel=1,nelec
  198     iworbd(iel,idet)=map_orb(iworbd(iel,idet))

c Sort eigenvalues and match eigenvalue sets to label dets. so that those dets. that
c can be in the same CSF have the same label.
      do 220 idet=1,ndet
        do 210 iel=1,nelec
  210     eigs_in_det(iel,idet)=eigs(iworbd(iel,idet))
        call shell(eigs_in_det(1,idet),nelec)
  220   write(6,'(''sorted eigs'',i3,90f10.6)') idet,(eigs_in_det(iel,idet),iel=1,nelec)

      do 225 idet=2,ndet
  225   label_det(idet)=0

      nlabel=1
      label_det(1)=1
      do 250 idet=2,ndet
        do 240 jdet=1,idet-1
          do 230 iel=1,nelec
            if(abs(eigs_in_det(iel,idet)-eigs_in_det(iel,jdet)).gt.eps) goto 240
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
      nedet=0
      ncsf=0
      do 300 i=1,ndet
        if(iflag_det(i).eq.1) goto 300
        ncsf=ncsf+1
        if(ncsf.gt.MCSF) stop 'ncsf>MCSF'
        ndet_in_csf(ncsf)=1
        iwdet_in_csf(1,ncsf)=i
        csf_coef(ncsf)=cdet(i)
        cdet_in_csf(1,ncsf)=1
        do 290 j=i+1,ndet
          if(iflag_det(j).eq.1 .or. label_det(i).ne.label_det(j)) goto 290
          do 280 itimes=2,8
            times=0.5d0*itimes
            do 280 isign=-1,1,2
              if(abs(cdet(j)-isign*times*cdet(i)).lt.min(eps,0.2d0*abs(cdet(j)))) then
                iflag_det(j)=1
                nedet=nedet+1
                iedet(1,nedet)=j
                iedet(2,nedet)=i
                frac(2,nedet)=isign*times
                ndet_in_csf(ncsf)=ndet_in_csf(ncsf)+1
                if(ndet_in_csf(ncsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(ncsf) > MDET_IN_CSF'
                iwdet_in_csf(ndet_in_csf(ncsf),ncsf)=j
                cdet_in_csf(ndet_in_csf(ncsf),ncsf)=isign*times
              endif
              if(abs(isign*times*cdet(j)-cdet(i)).lt.min(eps,0.2d0*abs(cdet(j))) .and. itimes.ne.2) then
                iflag_det(j)=1
                nedet=nedet+1
                iedet(1,nedet)=j
                iedet(2,nedet)=i
                frac(2,nedet)=1.d0/(isign*times)
                ndet_in_csf(ncsf)=ndet_in_csf(ncsf)+1
                if(ndet_in_csf(ncsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(ncsf) > MDET_IN_CSF'
                iwdet_in_csf(ndet_in_csf(ncsf),ncsf)=j
                cdet_in_csf(ndet_in_csf(ncsf),ncsf)=1.d0/(isign*times)
              endif
  280     continue
  290   continue
  300 continue

c Create old version inputs
      write(6,'(/,''Inputs for old version of CHAMP'')')
      write(6,'(i3,'' nedet'')') nedet
      if(nedet.gt.0) then
        write(fmt,'(''(''i3,''(i5,i4),a)'')') nedet
        write(6,fmt) ((iedet(k,i),k=1,2),i=1,nedet),' ((iedet(k,i),k=1,2),i=1,nedet)'
        write(fmt,'(''(''i3,''(f6.1,f5.1),a)'')') nedet
        write(6,fmt) ((frac(k,i),k=1,2),i=1,nedet),' ((frac(k,i),k=1,2),i=1,nedet)'
       else
        write(6,'(a)') '((iedet(k,i),k=1,2),i=1,nedet)'
        write(6,'(a)') '((frac(k,i),k=1,2),i=1,nedet)'
      endif

      nparmd=0
      do 330 i=2,ndet
        if(iflag_det(i).eq.0) then
          nparmd=nparmd+1
          iwdet(nparmd)=i
        endif
  330 continue
      write(6,'(i3,'' nparmd'')') nparmd
      if(nparmd.gt.0) then
        write(fmt,'(''(''i3,''i3,a)'')') nparmd
        write(6,fmt) (iwdet(ipar),ipar=1,nparmd),' (iwdet(ipar),ipar=1,nparmd)'
       else
        write(6,'(a)') '(iwdet(ipar),ipar=1,nparmd)'
      endif

c Create new version inputs
      write(6,'(/,''Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i3,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i3,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i3,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 340 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
  340   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c When the same determinant appears in more than one CSF then there can be linear
c relations between the coefs. of the determinants.  The following detects these and
c puts the icsf CSF into the jcsf and kcsf CSFs and eliminates the icsf CSF.
c Note: it is probably better to have icsf be the smallest index (largest coef) rather
c than the largest index (smallest coef) of the three because then one ends up with
c smaller CSF coefs. and if they are small enough then we can drop them.
  350 ncsf_tmp=ncsf
c     do 400 icsf=2,ncsf_tmp
      do 400 icsf=1,ncsf_tmp-2
c       do 400 jcsf=1,icsf-1
        do 400 jcsf=icsf+1,ncsf_tmp-1
c         do 400 kcsf=1,jcsf-1
          do 400 kcsf=jcsf+1,ncsf_tmp
            if(label_det(iwdet_in_csf(1,icsf)).eq.label_det(iwdet_in_csf(1,jcsf)) .and.
     &         label_det(iwdet_in_csf(1,icsf)).eq.label_det(iwdet_in_csf(1,kcsf))) then
            do 390 isignj=-1,1,2
              do 390 isignk=-1,1,2
                if(abs(csf_coef(icsf)-isignj*csf_coef(jcsf)-isignk*csf_coef(kcsf)).lt.eps) then
                  write(6,'(''icsf,jcsf,kcsf,csf_coef(icsf),csf_coef(jcsf),csf_coef(kcsf)='',3i3,3f6.3)')
     &            icsf,jcsf,kcsf,csf_coef(icsf),csf_coef(jcsf),csf_coef(kcsf)
                  if(ndet_in_csf(jcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(jcsf) > MDET_IN_CSF'
                  if(ndet_in_csf(kcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(kcsf) > MDET_IN_CSF'
                  do 370 idet_in_csf=1,ndet_in_csf(icsf)
                    iwdet_in_csf(ndet_in_csf(jcsf)+idet_in_csf,jcsf)=iwdet_in_csf(idet_in_csf,icsf)
                    iwdet_in_csf(ndet_in_csf(kcsf)+idet_in_csf,kcsf)=iwdet_in_csf(idet_in_csf,icsf)
                    cdet_in_csf(ndet_in_csf(jcsf)+idet_in_csf,jcsf)=isignj*cdet_in_csf(idet_in_csf,icsf)
  370               cdet_in_csf(ndet_in_csf(kcsf)+idet_in_csf,kcsf)=isignk*cdet_in_csf(idet_in_csf,icsf)
                  ndet_in_csf(jcsf)=ndet_in_csf(jcsf)+ndet_in_csf(icsf)
                  ndet_in_csf(kcsf)=ndet_in_csf(kcsf)+ndet_in_csf(icsf)
                  do 380 lcsf=icsf+1,ncsf_tmp
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

      write(6,'(/,''If any of the foll. numbers look rational then there was a failure in recognizing a relationship'')')
      do 405 icsf=1,ncsf-1
        do 405 jcsf=icsf+1,ncsf
  405     write(6,'(''icsf,jcsf,ratio='',2i4,9f12.6)')
     &    icsf,jcsf,(csf_coef(icsf)/csf_coef(jcsf))**2,(csf_coef(jcsf)/csf_coef(icsf))**2

c Keep only CSF's with csf_coef larger than cutoff_d2c and corresponding dets
c This is not done quite correctly because of the normalization of the CSFs
      ncsf_new=0
      ndet_new=0
      do 310 icsf=1,ncsf
        csf_norm=0
        do 305 idet_in_csf=1,ndet_in_csf(icsf)
  305     csf_norm=csf_norm+cdet_in_csf(idet_in_csf,icsf)**2
        csf_norm=sqrt(csf_norm)
        if(abs(csf_norm*csf_coef(icsf)).ge.cutoff_d2c) then
          ncsf_new=icsf
          do 307 idet_in_csf=1,ndet_in_csf(icsf)
  307       ndet_new=max(ndet_new,iwdet_in_csf(idet_in_csf,icsf))
        endif
  310 continue
      write(6,'(''ncsf,ncsf_new,ndet,ndet_new='',9i5)') ncsf,ncsf_new,ndet,ndet_new
      ncsf=ncsf_new
      ndet=ndet_new

      if(csf_coef(1).lt.0) then
        do 320 icsf=1,ncsf
  320     csf_coef(icsf)=-csf_coef(icsf)
      endif

c Create new version inputs
c While doing that, check if any cdet_in_csf is not an integer or half-integer.
c I need to do the checking completely differently using orbital energies or occupations.
      write(6,'(/,''Inputs for new version of CHAMP with linear relations recognized'')')

      write(6,'(3i4,t42,a)') ndet,nbasis,norb, 'ndet,nbasis,norb'
      write(6,'(2i4,a)') ndet,ncsf, ' ndet,ncsf'

      do 406 iorb=1,norb
        write(fmt,'(''(1p,'',i3''d16.8,a)'')') nbasis
        if(iorb.eq.1) then
          write(6,fmt) (coef(ibasis,iorb),ibasis=1,nbasis),' ((coef(ibasis,iorb),ibasis=1,nbasis),iorb=1,norb)'
         else
          write(6,fmt) (coef(ibasis,iorb),ibasis=1,nbasis)
        endif
  406 continue

      ndn=nelec-nup
      do 407 idet=1,ndet
        write(fmt,'(''(i3,'',i2,''i4,3x,'',i2,''i4,i5,a)'')') nup-1,ndn
        if(idet.eq.ndet) then
          write(6,fmt) (iworbd(iel,idet),iel=1,nelec), label_det(idet), ' (iworbd(iel,idet),iel=1,nelec), label_det(idet)'
         else
          write(6,fmt) (iworbd(iel,idet),iel=1,nelec), label_det(idet)
        endif
  407 continue
 
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i3,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i3,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i3,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 420 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        do 410 idet_in_csf=1,ndet_in_csf(icsf)
  410     if(abs(0.5d0*int(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-10)
     &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f8.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
  420   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
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
                it=iwdet_in_csf(idet_in_csf,i)
                iwdet_in_csf(idet_in_csf,i)=iwdet_in_csf(idet_in_csf,l)
    7           iwdet_in_csf(idet_in_csf,l)=it
   10       continue
   20     continue

c Calculate number of determinants with abs(cdet) >= 1.d-6
      ndet_new=0
      do 30 idet=1,ndet
   30   if(abs(cdet(idet)).ge.1.d-6) ndet_new=ndet_new+1
      write(6,'(''ndet,ndet_new='',9i5)') ndet,ndet_new
      write(6,'(''cdet='',100f8.4)') (cdet(idet),idet=1,ndet_new)
      ndet=ndet_new

      return
      end

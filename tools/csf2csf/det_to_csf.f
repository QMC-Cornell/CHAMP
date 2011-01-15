      program det_to_csf
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
      parameter(MDET=1000,MCSF=500,MDET_IN_CSF=40,MELEC=100,MORB=100,MACT=100)
      dimension cdet(MDET),iwdet(MDET),iedet(2,MDET),frac(2,MDET)
     &,iflag(MDET), csf_coef(MCSF),ndet_in_csf(MCSF)
     &,iwdet_in_csf(MDET_IN_CSF,MCSF),cdet_in_csf(MDET_IN_CSF,MCSF)
     &,iworbd(MELEC,MDET),eigs(MORB),eigs_in_det(MELEC,MDET),label_det(MDET)
      character(MORB) alphain(MDET),betain(MDET),ctmp,ctmp2
      integer alpha(MDET,MACT),beta(MDET,MACT)

c     read(5,*) mode
c     read(5,*) cutoff_g2q,cutoff_d2c
c     if(cutoff_d2c.lt.cutoff_g2q) stop 'cutoff_d2c must be >= cutoff_g2q'

      if(index(mode,'guga').ne.0) then

c eps should be something like 4.d-6
        read(5,*) ndet,nelec,norb,eps
        write(6,'(''ndet,nelec,norb,eps='',3i4,1pd9.2)') ndet,nelec,norb,eps
        if(ndet.gt.MDET) stop 'ndet>MDET'
        if(nelec.gt.MELEC) stop 'nelec>MELEC'
        if(norb.gt.MORB) stop 'norb>MORB'

        read(5,*) (cdet(i),i=1,ndet)
        write(6,'(''cdet='',20f9.6)') (cdet(i),i=1,ndet)

        read(5,*) (eigs(i),i=1,norb)
        write(6,'(''eigs(i)='',20f9.6)') (eigs(i),i=1,norb)

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
      write(6,'(''nlabel,label_det(idet)='',i3,2x,90i3)') nlabel,(label_det(idet),idet=1,ndet)

      do 160 i=1,ndet
        frac(1,i)=1
  160   iflag(i)=0

c Sometimes cdet can be very small even though csf_coef are not because different csf_coefs can add
c to produce a very small cdet.  So, while checking, check not only
c if(abs(cdet(j)-isign*times*cdet(i)).lt.eps) but rather
c if(abs(cdet(j)-isign*times*cdet(i)).lt.min(eps,0.2d0*abs(cdet(j))))
      nedet=0
      ncsf=0
      do 200 i=1,ndet
        if(iflag(i).eq.1) goto 200
        ncsf=ncsf+1
        if(ncsf.gt.MCSF) stop 'ncsf>MCSF'
        ndet_in_csf(ncsf)=1
        iwdet_in_csf(1,ncsf)=i
        csf_coef(ncsf)=cdet(i)
        cdet_in_csf(1,ncsf)=1
        do 190 j=i+1,ndet
          if(iflag(j).eq.1 .or. label_det(i).ne.label_det(j)) goto 190
          do 180 itimes=2,8
            times=0.5d0*itimes
            do 180 isign=-1,1,2
              if(abs(cdet(j)-isign*times*cdet(i)).lt.min(eps,0.2d0*abs(cdet(j)))) then
                iflag(j)=1
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
                iflag(j)=1
                nedet=nedet+1
                iedet(1,nedet)=j
                iedet(2,nedet)=i
                frac(2,nedet)=1.d0/(isign*times)
                ndet_in_csf(ncsf)=ndet_in_csf(ncsf)+1
                if(ndet_in_csf(ncsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(ncsf) > MDET_IN_CSF'
                iwdet_in_csf(ndet_in_csf(ncsf),ncsf)=j
                cdet_in_csf(ndet_in_csf(ncsf),ncsf)=1.d0/(isign*times)
              endif
  180     continue
  190   continue
  200 continue

      do 210 i=1,nedet
  210   write(6,'(2i4,2f7.2)') (iedet(k,i),k=1,2),(frac(k,i),k=1,2)

      if(csf_coef(1).lt.0) then
        do 120 icsf=1,ncsf
  120     csf_coef(icsf)=-csf_coef(icsf)
      endif

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
      do 130 i=2,ndet
        if(iflag(i).eq.0) then
          nparmd=nparmd+1
          iwdet(nparmd)=i
        endif
  130 continue
      write(6,'(i3,'' nparmd'')') nparmd
      write(fmt,'(''(''i3,''i3,a)'')') nparmd
      write(6,fmt) (iwdet(ipar),ipar=1,nparmd),' (iwdet(ipar),ipar=1,nparmd)'

c Create new version inputs
      write(6,'(/,''Inputs for new version of CHAMP'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i3,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i3,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i3,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 140 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
  140   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c When the same determinant appears in more than one CSF then there can be linear
c relations between the coefs. of the determinants.  The following detects these and
c puts the icsf CSF into the jcsf and kcsf CSFs and eliminates the icsf CSF.
c Note: it is probably better to have icsf be the smallest index (largest coef) rather
c than the largest index (smallest coef) of the three because then one ends up with
c smaller CSF coefs. and if they are small enough then we can drop them.
  150 ncsf_tmp=ncsf
c     do 200 icsf=2,ncsf_tmp
      do 200 icsf=1,ncsf_tmp-2
c       do 200 jcsf=1,icsf-1
        do 200 jcsf=icsf+1,ncsf_tmp-1
c         do 200 kcsf=1,jcsf-1
          do 200 kcsf=jcsf+1,ncsf_tmp
            if(label_det(iwdet_in_csf(1,icsf)).eq.label_det(iwdet_in_csf(1,jcsf)) .and.
     &         label_det(iwdet_in_csf(1,icsf)).eq.label_det(iwdet_in_csf(1,kcsf))) then
            do 190 isignj=-1,1,2
              do 190 isignk=-1,1,2
                if(abs(csf_coef(icsf)-isignj*csf_coef(jcsf)-isignk*csf_coef(kcsf)).lt.eps) then
                  write(6,'(''icsf,jcsf,kcsf,csf_coef(icsf),csf_coef(jcsf),csf_coef(kcsf)='',3i3,3f6.3)')
     &            icsf,jcsf,kcsf,csf_coef(icsf),csf_coef(jcsf),csf_coef(kcsf)
                  if(ndet_in_csf(jcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(jcsf) > MDET_IN_CSF'
                  if(ndet_in_csf(kcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(kcsf) > MDET_IN_CSF'
                  do 170 idet_in_csf=1,ndet_in_csf(icsf)
                    iwdet_in_csf(ndet_in_csf(jcsf)+idet_in_csf,jcsf)=iwdet_in_csf(idet_in_csf,icsf)
                    iwdet_in_csf(ndet_in_csf(kcsf)+idet_in_csf,kcsf)=iwdet_in_csf(idet_in_csf,icsf)
                    cdet_in_csf(ndet_in_csf(jcsf)+idet_in_csf,jcsf)=isignj*cdet_in_csf(idet_in_csf,icsf)
  170               cdet_in_csf(ndet_in_csf(kcsf)+idet_in_csf,kcsf)=isignk*cdet_in_csf(idet_in_csf,icsf)
                  ndet_in_csf(jcsf)=ndet_in_csf(jcsf)+ndet_in_csf(icsf)
                  ndet_in_csf(kcsf)=ndet_in_csf(kcsf)+ndet_in_csf(icsf)
                  do 180 lcsf=icsf+1,ncsf_tmp
                    csf_coef(lcsf-1)=csf_coef(lcsf)
                    ndet_in_csf(lcsf-1)=ndet_in_csf(lcsf)
                    do 180 idet_in_csf=1,ndet_in_csf(lcsf)
                      iwdet_in_csf(idet_in_csf,lcsf-1)=iwdet_in_csf(idet_in_csf,lcsf)
  180                 cdet_in_csf(idet_in_csf,lcsf-1)=cdet_in_csf(idet_in_csf,lcsf)
                  ncsf=ncsf-1
                  goto 150
                endif
  190       continue
            endif
  200 continue

      write(6,'(/,''If any of the foll. numbers look rational then there was a failure in recognizing a relationship'')')
      do 205 icsf=1,ncsf-1
        do 205 jcsf=icsf+1,ncsf
  205     write(6,'(''icsf,jcsf,ratio='',2i4,9f12.6)')
     &    icsf,jcsf,(csf_coef(icsf)/csf_coef(jcsf))**2,(csf_coef(jcsf)/csf_coef(icsf))**2

c Create new version inputs
c While doing that, check if any cdet_in_csf is not an integer or half-integer.
c I need to do the checking completely differently using orbital energies or occupations.
      write(6,'(/,''Inputs for new version of CHAMP with linear relations recognized'')')
      write(6,'(i5,'' ncsf'')') ncsf
c     write(fmt,'(''(''i3,''f12.8,\'\' (csf_coef(icsf),icsf=1,ncsf)\'\')'')') ncsf
      write(fmt,'(''(''i3,''f12.8,a)'')') ncsf
      write(6,fmt) (csf_coef(icsf),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      write(fmt,'(''(''i3,''i3,a)'')') ncsf
      write(6,fmt) (ndet_in_csf(icsf),icsf=1,ncsf),' (ndet_in_csf(icsf),icsf=1,ncsf)'
      do 220 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        do 210 idet_in_csf=1,ndet_in_csf(icsf)
  210     if(abs(0.5d0*int(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-10)
     &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f8.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
  220   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

      stop
      end
c-----------------------------------------------------------------------
      SUBROUTINE SHELL(D,N)
      IMPLICIT REAL*8(A-H,O-Z)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:::  SHELL-METZGER SORT IN ASCENDING ORDER.          ...CYRUS 1979  :::
C:::  MODIFIED SLIGHTLY FOR READIBILITY.          ...CYRUS 7 DEC 83  :::
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      DIMENSION D(N)
      LOGNB2=INT(DLOG(DFLOAT(N))/DLOG(2.D0)+1.D-14)
      M=N
      DO 20 NN=1,LOGNB2
        M=M/2
        K=N-M
        DO 20 J=1,K
          DO 10 I=J,1,-M
            L=I+M
            IF (D(L).GT.D(I)) GOTO 20
            T=D(I)
            D(I)=D(L)
            D(L)=T
   10       CONTINUE
   20     CONTINUE
      RETURN
      END

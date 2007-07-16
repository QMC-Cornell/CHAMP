c Reads in determinantal coefs. cdet, and checks to see which are
c equal or equal and opposite possibly with a multiplier of 2.
c Does not find all dependencies, when there are more than pairs
c of related coefs, that are not related in the above simple way.
c Output written in fashion to be useful to serve as input to
c both old and new versions of optimization program.
c Not clear what is best value of eps

c Written by Cyrus Umrigar

      implicit real*8(a-h,o-z)
      character*80 fmt
      parameter (MDET=1000,MCSF=500,MDET_IN_CSF=100)
      dimension cdet(MDET),iwdet(MDET),iedet(2,MDET),frac(2,MDET)
     &,iflag(MDET), csf_coef(MCSF),ndet_in_csf(MCSF)
     &,iwdet_in_csf(MDET_IN_CSF,MCSF),cdet_in_csf(MDET_IN_CSF,MCSF)

c eps should be something like 4.d-6
      read(5,*) ndet, eps
      if(ndet.gt.MDET) stop 'ndet>MDET'
      write(6,'(''ndet='',i4)') ndet
      read(5,*) (cdet(i),i=1,ndet)
      write(6,'(''cdet='',20f9.6)') (cdet(i),i=1,ndet)
      nedet=0

      do 10 i=1,ndet
        frac(1,i)=1
   10   iflag(i)=0

      ncsf=0
      do 40 i=1,ndet
        if(iflag(i).eq.1) goto 40
        ncsf=ncsf+1
        if(ncsf.gt.MCSF) stop 'ncsf>MCSF'
        ndet_in_csf(ncsf)=1
        iwdet_in_csf(1,ncsf)=i
        csf_coef(ncsf)=cdet(i)
        cdet_in_csf(1,ncsf)=1
        do 20 j=i+1,ndet
          if(iflag(j).eq.1) goto 20
c         if(abs(cdet(j)-cdet(i)).lt.eps) then
c           iflag(j)=1
c           nedet=nedet+1
c           iedet(1,nedet)=j
c           iedet(2,nedet)=i
c           frac(2,nedet)=1
c           ndet_in_csf(ncsf)=ndet_in_csf(ncsf)+1
c           if(ndet_in_csf(ncsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(ncsf) > MDET_IN_CSF'
c           iwdet_in_csf(ndet_in_csf(ncsf),ncsf)=j
c           cdet_in_csf(ndet_in_csf(ncsf),ncsf)=1
c         endif
c         if(abs(cdet(j)+cdet(i)).lt.eps) then
c           iflag(j)=1
c           nedet=nedet+1
c           iedet(1,nedet)=j
c           iedet(2,nedet)=-i
c           frac(2,nedet)=1
c           ndet_in_csf(ncsf)=ndet_in_csf(ncsf)+1
c           if(ndet_in_csf(ncsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(ncsf) > MDET_IN_CSF'
c           iwdet_in_csf(ndet_in_csf(ncsf),ncsf)=j
c           cdet_in_csf(ndet_in_csf(ncsf),ncsf)=-1
c         endif
          do 15 itimes=2,8
          times=0.5d0*itimes
          do 15 isign=-1,1,2
c           if(itimes.eq.3) write(6,'(''i,j,cdet(i),cdet(j)'',2i3,2f6.2)') i,j,cdet(i),cdet(j)
c           write(6,'(''i,j,cdet(i),cdet(j)'',2i3,2f6.2)') i,j,cdet(i),cdet(j)
c           if(itimes.eq.3 .and. i.gt.8 .and. j.gt.8) write(6,'(''diff'',d12.4)') abs(cdet(j)-isign*times*cdet(i))
c           if(itimes.eq.3 .and. i.gt.8 .and. j.gt.8) write(6,'(''diff'',d12.4)') abs(cdet(j)+isign*times*cdet(i))
c           if(itimes.eq.3 .and. i.gt.8 .and. j.gt.8) write(6,'(''diff'',d12.4)') abs(isign*times*cdet(j)+cdet(i))
c           if(itimes.eq.3 .and. i.gt.8 .and. j.gt.8) write(6,'(''diff'',d12.4)') abs(isign*times*cdet(j)-cdet(i))
            if(abs(cdet(j)-isign*times*cdet(i)).lt.eps) then
c             write(6,'(''e1'')')
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
            if(abs(isign*times*cdet(j)-cdet(i)).lt.eps .and. itimes.ne.2) then
c             write(6,'(''e3'')')
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
   15     continue
   20   continue
   40 continue

      do 50 i=1,nedet
   50   write(6,'(2i4,2f6.1)') (iedet(k,i),k=1,2),(frac(k,i),k=1,2)

      if(csf_coef(1).lt.0) then
        do 60 icsf=1,ncsf
   60     csf_coef(icsf)=-csf_coef(icsf)
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
      do 70 i=2,ndet
        if(iflag(i).eq.0) then
          nparmd=nparmd+1
          iwdet(nparmd)=i
        endif
   70 continue
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
      do 80 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
   80   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

c When the same determinant appears in more than one CSF then there can be linear
c relations between the coefs. of the determinants.  The following detects these and
c puts the icsf CSF into the jcsf and kcsf CSFs and eliminates the icsf CSF.
c Note: it is probably better to have icsf be the smallest index (largest coef) rather
c than the largest index (smallest coef) of the three because then one ends up with
c smaller CSF coefs. and if they are small enough then we can drop them.
   85 ncsf_tmp=ncsf
c     do 100 icsf=2,ncsf_tmp
      do 100 icsf=1,ncsf_tmp-2
c       do 100 jcsf=1,icsf-1
        do 100 jcsf=icsf+1,ncsf_tmp-1
c         do 100 kcsf=1,jcsf-1
          do 100 kcsf=jcsf+1,ncsf_tmp
c           if(icsf.ne.jcsf .and. icsf.ne.kcsf .and. jcsf.ne.kcsf) then
c      write(6,'(''i,j,k,csf_coef(icsf)-csf_coef(jcsf)-csf_coef(kcsf)'',3i2,d12.4)')
c    & icsf,jcsf,kcsf,csf_coef(icsf)-csf_coef(jcsf)-csf_coef(kcsf)
            do 98 isignj=-1,1,2
            do 98 isignk=-1,1,2
            if(abs(csf_coef(icsf)-isignj*csf_coef(jcsf)-isignk*csf_coef(kcsf)).lt.eps) then
              write(6,'(''icsf,jcsf,kcsf,csf_coef(icsf),csf_coef(jcsf),csf_coef(kcsf)='',3i3,3f6.3)')
     &        icsf,jcsf,kcsf,csf_coef(icsf),csf_coef(jcsf),csf_coef(kcsf)
              if(ndet_in_csf(jcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(jcsf) > MDET_IN_CSF'
              if(ndet_in_csf(kcsf)+ndet_in_csf(icsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(kcsf) > MDET_IN_CSF'
              do 90 idet_in_csf=1,ndet_in_csf(icsf)
                iwdet_in_csf(ndet_in_csf(jcsf)+idet_in_csf,jcsf)=iwdet_in_csf(idet_in_csf,icsf)
                iwdet_in_csf(ndet_in_csf(kcsf)+idet_in_csf,kcsf)=iwdet_in_csf(idet_in_csf,icsf)
                cdet_in_csf(ndet_in_csf(jcsf)+idet_in_csf,jcsf)=isignj*cdet_in_csf(idet_in_csf,icsf)
   90           cdet_in_csf(ndet_in_csf(kcsf)+idet_in_csf,kcsf)=isignk*cdet_in_csf(idet_in_csf,icsf)
              ndet_in_csf(jcsf)=ndet_in_csf(jcsf)+ndet_in_csf(icsf)
              ndet_in_csf(kcsf)=ndet_in_csf(kcsf)+ndet_in_csf(icsf)
c             if(ndet_in_csf(jcsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(jcsf) > MDET_IN_CSF'
c             if(ndet_in_csf(kcsf).gt.MDET_IN_CSF) stop 'ndet_in_csf(kcsf) > MDET_IN_CSF'
              do 95 lcsf=icsf+1,ncsf_tmp
                csf_coef(lcsf-1)=csf_coef(lcsf)
                ndet_in_csf(lcsf-1)=ndet_in_csf(lcsf)
                do 95 idet_in_csf=1,ndet_in_csf(lcsf)
                  iwdet_in_csf(idet_in_csf,lcsf-1)=iwdet_in_csf(idet_in_csf,lcsf)
   95             cdet_in_csf(idet_in_csf,lcsf-1)=cdet_in_csf(idet_in_csf,lcsf)
              ncsf=ncsf-1
              goto 85
            endif
   98       continue
c           endif
  100 continue

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
      do 120 icsf=1,ncsf
        write(fmt,'(''(''i3,''i5,a)'')') ndet_in_csf(icsf)
        write(6,fmt) (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
        do 110 idet_in_csf=1,ndet_in_csf(icsf)
  110     if(abs(0.5d0*int(2*cdet_in_csf(idet_in_csf,icsf))-cdet_in_csf(idet_in_csf,icsf)).gt.1.d-10)
     &    write(6,'(''Warning:cdet_in_csf(idet_in_csf,icsf)='',f8.6)') cdet_in_csf(idet_in_csf,icsf)
        write(fmt,'(''(''i3,''f5.1,a)'')') ndet_in_csf(icsf)
  120   write(6,fmt) (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf)),
     &  ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'

      stop
      end

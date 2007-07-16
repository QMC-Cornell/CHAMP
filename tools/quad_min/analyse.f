      program analyse
      implicit real*8(a-h,o-z)
c Written by Cyrus Umrigar
c Analyse fluctuations in various terms of hessian
      parameter(MPARM=100,MPARM2=MPARM*(MPARM+1)/2,MTERM=20)

      character*20 filename
      character*64 title

      dimension hess_pieces(MTERM),hess(MPARM2,MTERM),hess2(MPARM2,MTERM),jparm(MPARM2),kparm(MPARM2)

c Read in 9 terms.  Calculate 2 more as linear comninations
      nterm=6
      nterm2=nterm+3
c     nterm2=nterm+1

      read(5,*) nparm,irun_beg,nrun
      write(6,'(''nparm,irun_beg,nrun='',9i4)') nparm,irun_beg,nrun
      nparm2=nparm*(nparm+1)/2
      if(nparm.gt.MPARM) stop 'parm>MPARM'
      if(nterm2.gt.MTERM) stop 'nterm2>MTERM'

      do 10 iparm2=1,nparm2
        do 10 iterm=1,nterm2
          hess(iparm2,iterm)=0
   10     hess2(iparm2,iterm)=0

      do 50 irun=irun_beg,irun_beg+nrun-1
        if(irun.le.9) then
          write(filename,'(''hess_pieces.'',i1)') irun
         elseif(irun.le.99) then
          write(filename,'(''hess_pieces.'',i2)') irun
         elseif(irun.le.999) then
          write(filename,'(''hess_pieces.'',i3)') irun
        endif
        open(22,file=filename,status='unknown')
        do 15 i=1,3
   15     read(22,*) title
        do 30 iparm2=1,nparm2
          read(22,*) jparm(iparm2),kparm(iparm2),(hess_pieces(iterm),iterm=1,nterm)
c         write(6,*) jparm(iparm2),kparm(iparm2),(hess_pieces(iterm),iterm=1,nterm)
          do 20 iterm=1,nterm
            hess(iparm2,iterm)=hess(iparm2,iterm)+hess_pieces(iterm)
   20       hess2(iparm2,iterm)=hess2(iparm2,iterm)+hess_pieces(iterm)**2
            term=hess_pieces(2)+hess_pieces(3)+hess_pieces(4)+hess_pieces(5)
            hess(iparm2,nterm+1)=hess(iparm2,nterm+1)+term
            hess2(iparm2,nterm+1)=hess2(iparm2,nterm+1)+term**2
c           term=hess_pieces(1)+2*hess_pieces(2)+hess_pieces(3)+hess_pieces(4)+hess_pieces(5)
            term=2*hess_pieces(2)+hess_pieces(3)+hess_pieces(4)+hess_pieces(5)
            hess(iparm2,nterm+2)=hess(iparm2,nterm+2)+term
            hess2(iparm2,nterm+2)=hess2(iparm2,nterm+2)+term**2
            term=hess_pieces(4)+hess_pieces(5)
            hess(iparm2,nterm+3)=hess(iparm2,nterm+3)+term
   30       hess2(iparm2,nterm+3)=hess2(iparm2,nterm+3)+term**2
   50 close(22)

      do 60 iparm2=1,nparm2
        do 60 iterm=1,nterm2
          hess(iparm2,iterm)=hess(iparm2,iterm)/nrun
          hess2(iparm2,iterm)=hess2(iparm2,iterm)/nrun
   60     hess2(iparm2,iterm)=sqrt(hess2(iparm2,iterm)-hess(iparm2,iterm)**2)

      write(6,'('' j  k   1         2         3         4         5         6        7         8          9     6/7   6/8   8/9   
     &4/9'')')
      write(6,'('' j  k   1         2         3         4         5         6      2+3+4+5  2*2+3+4+5    4+5    6/7   6/8   8/9   
     &4/9'')')

      do 70 iparm2=1,nparm2
c  70   write(6,'(i2,i3,1p9g11.3,0p9f5.1)') jparm(iparm2),kparm(iparm2),(hess2(iparm2,iterm),iterm=1,nterm2)
   70   write(6,'(i2,i3,f7.3,8f10.4,9f6.2)') jparm(iparm2),kparm(iparm2),(hess2(iparm2,iterm),iterm=1,nterm2)
     &   ,hess2(iparm2,nterm)/hess2(iparm2,nterm+1)
     &   ,hess2(iparm2,nterm)/hess2(iparm2,nterm+2)
     &   ,hess2(iparm2,nterm+2)/hess2(iparm2,nterm+3)
     &   ,hess2(iparm2,4)/hess2(iparm2,nterm+3)

      write(6,*)
      write(6,'('' j  k   1         2         3         4         5         6        7         8          9     6/7   6/8   8/9   
     &4/9'')')
      write(6,'('' j  k   1         2         3         4         5         6      2+3+4+5  2*2+3+4+5    4+5    6/7   6/8   8/9   
     &4/9'')')

      do 90 iparm2=1,nparm2
c  90   write(6,'(i2,i3,1p9g11.3,0p9f6.2)') jparm(iparm2),kparm(iparm2),(hess(iparm2,iterm),iterm=1,nterm2)
   90   write(6,'(i2,i3,f7.3,8f10.4,9f6.2)') jparm(iparm2),kparm(iparm2),(hess(iparm2,iterm),iterm=1,nterm2)
     &  ,hess(iparm2,nterm)/hess(iparm2,nterm+1)
     &  ,hess(iparm2,nterm)/hess(iparm2,nterm+2)
     &  ,hess(iparm2,nterm+2)/hess(iparm2,nterm+3)
     &  ,hess(iparm2,4)/hess(iparm2,nterm+3)
 
      stop
      end

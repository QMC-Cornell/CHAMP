      program linear
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'fit.h'
      parameter(MDIM=100,eps=1.d-12)
      character*50 fmt
      common /atom/ nctype
      common /jaspar4/ c(83,MCTYPE,MWF),norda,nordb,nordc
      common /orbit/ lo(MORB),n(MBASIS),l(MBASIS),npoint(MORB),
     &iwjasa(83,NCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),
     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
     &iwdet(MPARM),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
     &iedet(ICX,MDET),imnbas(MCENT),
     &nparml,nparme,nparmd,nparms,nparmg,nparm_read,nparmj,
     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
     &nebase,necn,nedet

      dimension hmat(MDIM,MDIM),overlap(MDIM,MDIM),product(MDIM,MDIM),
     &overlap_sav(MDIM,MDIM),coef(MDIM,MDIM)
      dimension b(MORDJ1,2),a(MORDJ1,MCTYPE)

      parameter(MWORK=MDIM*MDIM+1024)
      dimension scratch(MDIM,MDIM)
      dimension work(MWORK),eigv(MDIM),eigvi(MDIM),eig_denom(MDIM)

c     read(5,*) ntimes
c     read(5,*) nctype
      read(5,*) nparmj,isym
c     read(5,*) nparmj,(nparma(ict),ict=1,nctype),nparmb(1),(nparmc(ict),ict=1,nctype)
      read(5,*) ((hmat(i,j),j=1,nparmj+1),i=1,nparmj+1)
      read(5,*) ((overlap(i,j),j=1,i),i=1,nparmj+1)

c for high printout, set idebug=1
      idebug=1
c for symmetrization, set isym=1
c     isym=0
c to get eigv different than the lowest, set e_target >= dabs(eigv)
      e_target=0

      do 1 i=1,MDIM
        do 1 j=1,MDIM
    1     coef(i,j)=0

      nparmj=nparmj+1

      do 2 i=1,nparmj
        overlap_sav(i,i)=overlap(i,i)
        do 2 j=1,i-1
          overlap_sav(i,j)=overlap(i,j)
          overlap_sav(j,i)=overlap(i,j)
    2     overlap(j,i)=overlap(i,j)

      if(isym.eq.1) then
        do 3 i=1,nparmj
          do 3 j=1,i-1
            hmat(i,j)=0.5d0*(hmat(i,j)+hmat(j,i))
    3       hmat(j,i)=hmat(i,j)
      endif

      if(idebug.gt.0) then
        do 6 i=1,nparmj
    6     write(6,'(''o      '',83e16.8)') (overlap(i,j),j=1,i)
        write(6,*)

        do 8 i=1,nparmj
    8     write(6,'(''hmat   '',83e16.8)') (hmat(i,j),j=1,nparmj)
        write(6,*)
      endif

      write(6,'(''diagonal H/O'',9f9.4)') (hmat(i,i)/overlap(i,i),i=1,nparmj)


      write(*,*) 'nparmj=',nparmj
      do j = 1, nparmj
      write(*,*) 'hmat=',(hmat(i,j),i=1,nparmj)
      enddo
      do j = 1, nparmj
      write(*,*) 'overlap=',(overlap(i,j),i=1,nparmj)
      enddo
      call dggev('N','V',nparmj, hmat, MDIM, overlap, MDIM, eigv, eigvi,
     &                  eig_denom, eigl, MDIM, coef, MDIM, work, -1, info )
      if(work(1).gt.MWORK) stop 'work(1).gt.MWORK'
      if(info.ne.0) stop 'info from dggev != 0'
      lwork=work(1)
      write(6,'(''optimal lwork='',i6)') lwork
      call dggev('N','V',nparmj, hmat, MDIM, overlap, MDIM, eigv, eigvi,
     &                  eig_denom, eigl, MDIM, coef, MDIM, work, lwork, info )
      if(info.ne.0) write(6,*) 'Warning dgeev: info =', info
      if(info.ne.0) stop 'info from dggev != 0'

c First check that MWORK is large enough
c     call dgeev('N','V',nparmj,product,MDIM,eigv,eigvi,scratch, MDIM,
c    &            coef, MDIM, work, -1,info) 
c     write(6,'(''optimal MWORK='',f8.1)') work(1)
c     if(work(1).gt.MWORK) stop 'work(1).gt.MWORK'
c     call dgeev('N','V',nparmj,product,MDIM,eigv,eigvi,scratch, MDIM,
c    &            coef, MDIM, work, MWORK,info) 
c     if(info.ne.0) write(6,*) 'Warning dgeev: info =', info

      write(6,'(''eigenvalue_denominators in dggev='',20d10.2)') (eig_denom(i),i=1,nparmj)
      do 15 i=1,nparmj
        eigv(i)=eigv(i)/eig_denom(i)
        eigvi(i)=eigvi(i)/eig_denom(i)
   15   write(6,'(''eigenvalue '',f12.6,'' + i* '',f10.6)')
     &  eigv(i),eigvi(i)
      write(6,*)

c The foll. is confusing
c     i0=0
c     e0=e_target
c     do 20 i=1,nparmj
c       if(e_target+eigv(i).lt.e0) then
c         e0=e_target+eigv(i)
c         i0=i
c       endif
c  20 continue
c     if(i0.eq.0) stop 'No eigenvalues below 0'

      coef_max=0
      do 20 i=1,nparmj
        if(abs(coef(1,i)).gt.coef_max) then
          coef_max=abs(coef(1,i))
          i0=i
        endif
   20 continue

      do 22 i=1,min(nparmj,5)
   22   write(6,'(''coef'',9f9.6)') (coef(j,i),j=1,nparmj)

      write(6,'(''eigv  '',f10.3)') eigv(i0)
      write(6,'(''coef0 '',100f10.3)') coef(1,i0)
      dnorm=1.d0/coef(1,i0)
      write(6,'(''dparm '',100f10.6)') (coef(i,i0)*dnorm,i=2,nparmj)
      write(6,*)

c Add change to old parameters
c     open(1,file='jas_old')
c     read(1,*) norda,nordb,nordc
c     read(1,*) (nparma(ict),ict=1,nctype),nparmb(1),(nparmc(ict),ict=1,nctype)
c     nparma_read=2+max(0,norda-1)
c     nparmb_read=2+max(0,nordb-1)
c     nparmc_read=nterms4(nordc)

c     do 46 ict=1,nctype
c  46   read(1,*) (a(i,ict),i=1,nparma_read)
c     read(1,*) (b(i,1),i=1,nparmb_read)
c     do 47 ict=1,nctype
c  47   read(1,*) (c(i,ict,1),i=1,nparmc_read)
c     do 48 ict=1,nctype
c  48   read(1,*) (iwjasa(iparm,ict),iparm=1,nparma(ict))
c     read(1,*) (iwjasb(iparm,1),iparm=1,nparmb(1))
c     do 49 ict=1,nctype
c  49   read(1,*) (iwjasc(iparm,ict),iparm=1,nparmc(ict))
c     call cuspinit4(1)

cc    nparm=0
c     iparm=0
c     do 50 ict=1,nctype
c       do 50 i=1,nparma(ict)
c         iparm=iparm+1
c  50     a(iwjasa(i,ict),ict)=a(iwjasa(i,ict),ict)+coef(iparm+1,i0)*dnorm
c     do 60 i=1,nparmb(1)
c       iparm=iparm+1
c  60   b(iwjasb(i,1),1)=b(iwjasb(i,1),1)+coef(iparm+1,i0)*dnorm
c     do 70 ict=1,nctype
c       do 70 i=1,nparmc(ict)
c         iparm=iparm+1
c  70     c(iwjasc(i,ict),ict,1)=c(iwjasc(i,ict),ict,1)+coef(iparm+1,i0)*dnorm
c     call cuspexact4(1)

c     open(2,file='jas_new')

c     if(nparma_read.gt.0) then
cc      write(fmt,'(''(''i2,''f13.8,\'\' (a(iparmj),iparmj=1,nparma)\'\')'')') nparma_read
c       write(fmt,'(''(''i2,''f13.8,a28)'')') nparma_read
c      else
cc      write(fmt,'(''(\'\' (a(iparmj),iparmj=1,nparma)\'\')'')')
c       write(fmt,'(''(a28)'')')
c     endif
c     do 80 ict=1,nctype
c  80   write(2,fmt) (a(i,ict),i=1,nparma_read),' (a(iparmj),iparmj=1,nparma)'

c     if(nparmb_read.gt.0) then
cc      write(fmt,'(''(''i2,''f13.8,\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')') nparmb_read
c       write(fmt,'(''(''i2,''f13.8,a28)'')') nparmb_read
c      else
cc      write(fmt,'(''(\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')')
c       write(fmt,'(''(a28)'')')
c     endif
c     write(2,fmt) (b(i,1),i=1,nparmb_read),' (b(iparmj),iparmj=1,nparmb)'

c     if(nparmc_read.gt.0) then
cc      write(fmt,'(''(''i2,''f13.8,\'\' (c(iparmj),iparmj=1,nparmc)\'\')'')') nparmc_read
c       write(fmt,'(''(''i2,''f13.8,a28)'')') nparmc_read
c      else
cc      write(fmt,'(''(\'\' (c(iparmj),iparmj=1,nparmc)\'\')'')')
c       write(fmt,'(''(a28)'')')
c     endif
c     do 90 ict=1,nctype
c  90   write(2,fmt) (c(i,ict,1),i=1,nparmc_read),' (c(iparmj),iparmj=1,nparmc)'

      stop
      end
c-----------------------------------------------------------------------
      function nterms4(nord)
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)

      i=0
      do 20 n=2,nord
        do 20 k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do 20 l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              i=i+1
            endif
   20 continue
      nterms4=i
c     write(6,'(''nterms4='',i5)') nterms4
      return
      end
c-----------------------------------------------------------------------
      subroutine cuspinit4(iprin)
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

c     parameter(NEQSX=2*(MORDJ-1),MTERMS=55)
      parameter(NEQSX=6*MORDJ,MTERMS=55)

      common /jaspar4/ c(83,MCTYPE,MWF),norda,nordb,nordc
      common /cuspmat4/ d(NEQSX,MTERMS),iwc4(NEQSX),nterms

      if(nordc.eq.0) return

      do 10 n=1,2*(nordc-1)
        iwc4(n)=0
        do 10 i=1,MTERMS
   10   d(n,i)=0

      i=0
      do 20 n=2,nordc
        do 20 k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do 20 l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              i=i+1
              if(i.gt.MTERMS) stop 'nterms>MTERMS in cuspinit4'
              if(k.eq.1.and.iwc4(n-1).eq.0) iwc4(n-1)=i
              if(k.eq.0.and.iwc4(n+nordc-2).eq.0) iwc4(n+nordc-2)=i
c             write(6,'(9i4)') i,n,k,l,m
              iord=l+2*m
              d(iord,i)=d(iord,i)-2*k
              if(l+m.gt.0) then
                iord=nordc-1+k+m
                d(iord,i)=d(iord,i)-(l+m)
              endif
              if(m.gt.0) then
                iord=nordc-1+k+l+m
                d(iord,i)=d(iord,i)-m
              endif
            endif
   20 continue
      nterms=i
c     write(6,'(''# of e-e-n terms, nterms='',i5)') nterms
c     write(6,'(''d matrix:'')')
c     write(6,'(55i2)') (i,i=1,nterms)
c     do 30 n=1,2*(nordc-1)
c  30   write(6,'(55i2)') (nint(d(n,i)),i=1,nterms)
c     write(6,'(''coefs. fixed by cusp conditions are'')')
c     write(6,'(55i3)') (iwc4(i),i=1,2*(nordc-1))

      call checkdepend4

      return
      end
c-----------------------------------------------------------------------
      subroutine checkdepend4

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'fit.h'
      include 'force.h'

c     parameter(NEQSX=2*(MORDJ-1),MTERMS=55)
      parameter(NEQSX=6*MORDJ,MTERMS=55)

      common /jaspar4/ c(83,MCTYPE,MWF),norda,nordb,nordc
      common /cuspmat4/ d(NEQSX,MTERMS),iwc4(NEQSX),nterms

      common /orbit/ lo(MORB),n(MBASIS),lb(MBASIS),npoint(MORB),
     &iwjasa(83,NCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),
     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
     &iwdet(MPARM),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
     &iedet(ICX,MDET),imnbas(MCENT),
     &nparml,nparme,nparmd,nparms,nparmg,nparm_read,nparmj,
     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
     &nebase,necn,nedet

      common /atom/ nctype

      common /vardep/ nvdepend(NEQSX,MCTYPE),iwdepend(NEQSX,83,MCTYPE)
     &,cdep(NEQSX,83,MCTYPE)

      neqs=2*(nordc-1)

      do 2 i=1,neqs
        do 2 it=1,nctype
          do 2 l=1,nparmc(it)
            if(iwjasc(l,it) .eq. iwc4(i)) then
              write(6,'(''Error: J_een parameter'',i3,
     &        '' is dependent and should not be varied'')') iwjasc(l,it)
              stop 'You are trying to vary a dependent J_een parameter'
            endif
    2 continue

      do 5 i=1,neqs
        do 5 it=1,nctype
          nvdepend(i,it)=0
          do 5 l=1,nparmc(it)
    5       iwdepend(i,l,it)=0

c Figure out dependence of all dependent variables except that from
c the 2nd order e-n cusp cond.
      do 40 i=1,neqs
        if(i.eq.nordc) goto 40
        do 20 j=1,nterms
          if(dabs(d(i,j)).gt.1.d-10) then
            do 10 it=1,nctype
              do 10 l=1,nparmc(it)
                if(j.eq.iwjasc(l,it)) then
                  if(j.eq.iwc4(i)) stop 'do not vary dep. variable'
                  nvdepend(i,it)=nvdepend(i,it)+1
                  iwdepend(i,nvdepend(i,it),it)=l
                  cdep(i,nvdepend(i,it),it)=-d(i,j)/d(i,iwc4(i))
                endif
   10       continue
          endif
   20   continue
   40 continue

c Now do the one from the 2nd order e-n cusp cond.
c The dep. variable from the 2nd-order e-n cc depends directly on all the
c other dependent variables and only on the other dependent variables.
c Figure out what dependence on the independent variables is implied by that.

c Since it depends directly on all other dependent variables, it depends
c indirectly on all the independent variables that are being varied.
      do 50 it=1,nctype
      do 50 l=1,nparmc(it)
        nvdepend(nordc,it)=nvdepend(nordc,it)+1
   50   iwdepend(nordc,nvdepend(nordc,it),it)=l

      do 70 i=1,neqs
        if(i.eq.nordc) goto 70
        factor=-d(nordc,iwc4(i))/d(nordc,iwc4(nordc))
        do 60 j=1,nterms
          do 60 it=1,nctype
            do 60 l=1,nparmc(it)
              if(j.eq.iwjasc(l,it)) then
                cdep(nordc,l,it)= cdep(nordc,l,it)-
     &          factor*(d(i,j)/d(i,iwc4(i)))
              endif
   60   continue
   70 continue

c     write(6,'(''i  it nvdepend iwdepend or cdep'')')
c     do 80 i=1,neqs
c       do 80 it=1,nctype
c         write(6,'(i2,2i4,2x,15i4)') i,it,nvdepend(i,it),
c    &   (iwdepend(i,j,it),j=1,nvdepend(i,it))
c  80     write(6,'(i2,2i4,2x,15f4.0)') i,it,nvdepend(i,it),
c    &   (cdep(i,j,it),j=1,nvdepend(i,it))

      return
      end
c-----------------------------------------------------------------------
      subroutine cuspexact4(iprin)
c Written by Cyrus Umrigar

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

c For Jastrow4 NEQSX=2*(MORDJ-1) is sufficient.
c For Jastrow3 NEQSX=2*MORDJ should be sufficient.
c I am setting NEQSX=6*MORDJ simply because that is how it was for
c Jastrow3 for reasons I do not understand.
c     parameter(NEQSX=2*(MORDJ-1),MTERMS=55)
      parameter(NEQSX=6*MORDJ,MTERMS=55)

c The last 2 columns are what we care about in the foll. table
c------------------------------------------------------------------------------
c ord  # of   cumul.    # terms  # terms   # 3-body  Cumul. #      Cumul # indep
c      terms  # terms   even t   odd t      terms    3-body terms  3-body terms
c  n  (n+1)* (n^3+5n)/6         int((n+1)/2            nterms
c    (n+2)/2  +n^2+n           *int((n+2)/2
c------------------------------------------------------------------------------
c  1     3       3        2         1          0         0              0
c  2     6       9        4         2          2         2              0
c  3    10      19        6         4          4         6              2
c  4    15      34        9         6          7        13              7
c  5    21      55       12         9         10        23             15
c  6    28      83       16        12         14        37             27
c  7    36     119       20        16         18        55             43
c------------------------------------------------------------------------------

c Dependent coefs. fixed by e-e and e-n cusp conditions resp. are;
c order:   2  3  4  5  6  7  2  3  4  5  6  7
c coefs:   1  4 10 19 32 49  2  6 12 22 35 53

c So the terms varied for a 5th, 6th order polynomial are:
c    3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 (iwjasc(iparm),iparm=1,nparmc)
c    3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 24 25 26 27 28 29 30 31    33 34    36 37 (iwjasc(iparm),iparm=1,nparmc)


c All the dependent variables, except one (the one from the 2nd order
c e-n cusp) depend only on independent variables.  On the other hand
c the one from the 2nd order e-n cusp depends only on other dependent
c variables.

      common /atom/ nctype
      common /jaspar4/ c(83,MCTYPE,MWF),norda,nordb,nordc
      common /cuspmat4/ d(NEQSX,MTERMS),iwc4(NEQSX),nterms

      do 100 it=1,nctype

c Set dep. variables from e-e cusp
        do 20 i=1,nordc-1
          sum=0
          do 10 j=1,nterms
   10       if(j.ne.iwc4(i)) sum=sum+d(i,j)*c(j,it,1)
   20     c(iwc4(i),it,1)=-sum/d(i,iwc4(i))

c Set dep. variables from 3rd and higher order e-n cusp
        do 40 i=nordc+1,2*(nordc-1)
          sum=0
          do 30 j=1,nterms
   30       if(j.ne.iwc4(i)) sum=sum+d(i,j)*c(j,it,1)
   40     c(iwc4(i),it,1)=-sum/d(i,iwc4(i))

c Set dep. variables from 2nd order e-n cusp
        if(nordc.gt.1) then
          i=nordc
          sum=0
          do 50 j=1,nterms
   50       if(j.ne.iwc4(i)) sum=sum+d(i,j)*c(j,it,1)
          c(iwc4(i),it,1)=-sum/d(i,iwc4(i))
        endif
c     write(6,'(''nordc,nctype,nterms='',9i5)') nordc,nctype,nterms

 100  continue

      return
      end

c-----------------------------------------------------------------------
      subroutine matinv(a,mx,n)
      implicit real*8(a-h,o-z)

c routine to calculate inverse of matrix a
c assumed to be dimensioned a(mx,mx).
c the matrix a is replaced by its inverse.

      dimension a(mx,n)

      call dpotrf('L',n,a,mx,info)
      if (info.ne.0) then
       write(6,*) ' dpotrf: info =', info
      endif
      call dpotri('L',n,a,mx,info)
      if (info.ne.0) then
       write(6,*)  ' dpotri: info =', info
      endif

      do 10 i=1,n
        do 10 j=i+1,n
   10   a(i,j)=a(j,i)

c     do 20 i=1,n
c  20   write(6,'(10e10.3)') (a(i,j),j=1,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine prod(a,b,c,mx,n)
      implicit real*8(a-h,o-z)

      dimension a(mx,n),b(mx,n),c(mx,n)

      do 10 j=1,n
        do 10 i=1,n
          c(i,j)=0.d0
          do 10 k=1,n
  10        c(i,j)=c(i,j)+a(i,k)*b(k,j)

c     call dgemm('N','N',n,n,n,1.d0,a,mx,
c    &           b, mx, 0.d0, c, mx)

      return
      end

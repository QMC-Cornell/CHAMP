      program newton_chol
c Used to find the change in the wavefunction parameters.
c Given gradient, grad, and Hessian, hess, solve hess*x=grad for x.
c The required change in the parameters is -x.
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'fit.h'
      parameter(MDIM=100)
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

      dimension hess(MDIM,MDIM),grad(MDIM),dgrad(MDIM),
     &hess_sav(MDIM,MDIM),grad_sav(MDIM),grad_cal(MDIM)
      dimension index(MDIM),indexi(MDIM)
      dimension b(MORDJ1,2),a(MORDJ1,MCTYPE)

      read(5,*) ntimes
      do 40 itimes=1,ntimes
      read(5,*) nctype,istdev
      read(5,*) nparmj,add_diag
      read(5,*) (grad(i),i=1,nparmj)
      read(5,*) ((hess(i,j),j=1,i),i=1,nparmj)
      if(istdev.gt.0) read(5,*) (dgrad(i),i=1,nparmj)

      nparmj_sav=nparmj
      igood=0
      do 1 i=1,nparmj
        indexi(i)=0
        if(dabs(grad(i)).gt.istdev*dgrad(i)) then
          igood=igood+1
          index(igood)=i
          indexi(i)=igood
        endif
    1 continue
      nparmj=igood

      do 5 i=1,nparmj
        hess_sav(i,i)=hess(index(i),index(i))
        grad_sav(i)=grad(index(i))
        do 5 j=1,i-1
          hess_sav(i,j)=hess(index(i),index(j))
          hess_sav(j,i)=hess(index(i),index(j))
    5     hess(index(j),index(i))=hess(index(i),index(j))
      write(6,*)
      write(6,'(''nparmj='',i5)') nparmj
      write(6,'(''grad='',9g12.4)') (grad_sav(i),i=1,nparmj)
      do 6 i=1,nparmj
    6   write(6,'(''hess='',9g12.4)') (hess_sav(i,j),j=1,i)

c Do cholesky decomposition
c If A is not positive definite, iterate on it, each time
c adding a constant, add_diag, to diagonal till it is positive definite.
c We can come up with a better method later.
c Peter suggests first diagonalizing A by a unitary transform
c taking the inverse of the diagonalized matrix to construct the
c inverse of A, but resetting the diagonal elements of the diagonalized A
c to a minimum positive value if they are less than this value.
c Warning: temp
      ipos_def=0
      do 40 iter=1,100000
c     do 10 iter=1,50
    7   do 8 i=1,nparmj
          grad(i)=grad_sav(i)
          do 8 j=1,nparmj
    8     hess(i,j)=hess_sav(i,j)
        call chlsky(hess,nparmj,MDIM,ierr)
c Warning: temp
c       if(ierr.eq.0) goto 15
c         do 10 i=1,nparmj
c  10       hess_sav(i,i)=hess_sav(i,i)+add_diag
        if(ierr.ne.0) goto 35

      ipos_def=ipos_def+1
c Symmetrize decomposed matrix (needs to be done before calling uxb
c or need to modify uxb)
   15 do 20 i=1,nparmj
         do 20 j=i+1,nparmj
   20       hess(i,j)=hess(j,i)

c Solve equations
      call lxb(hess,nparmj,MDIM,grad)
      call uxb(hess,nparmj,MDIM,grad)

c test solution
      do 30 i=1,nparmj
        grad_cal(i)=0
        do 30 j=1,nparmj
   30     grad_cal(i)=grad_cal(i)+hess_sav(i,j)*grad(j)
      write(6,'(''grad_cal='',9g12.4)') (grad_cal(i),i=1,nparmj)
      write(6,'(''-x='',9f15.9)') (-grad(i),i=1,nparmj)

   35 do 38 i=1,nparmj
   38   hess_sav(i,i)=hess_sav(i,i)+add_diag
c  40 if(iter.ge.2 .and. ierr.eq.0) goto 45
   40 if(ipos_def.ge.2) goto 45

   45 write(6,*) 'converged in iterations ',iter

c Add change to old parameters
      open(1,file='jas_old')
      read(1,*) norda,nordb,nordc
      read(1,*) (nparma(ict),ict=1,nctype),nparmb(1),(nparmc(ict),ict=1,nctype)
c     nparma_read=norda+1
c     nparmb_read=nordb+1
      nparma_read=2+max(0,norda-1)
      nparmb_read=2+max(0,nordb-1)
      nparmc_read=nterms4(nordc)

      do 46 ict=1,nctype
   46   read(1,*) (a(i,ict),i=1,nparma_read)
      read(1,*) (b(i,1),i=1,nparmb_read)
      do 47 ict=1,nctype
   47   read(1,*) (c(i,ict,1),i=1,nparmc_read)
      do 48 ict=1,nctype
   48   read(1,*) (iwjasa(iparm,ict),iparm=1,nparma(ict))
      read(1,*) (iwjasb(iparm,1),iparm=1,nparmb(1))
      do 49 ict=1,nctype
   49   read(1,*) (iwjasc(iparm,ict),iparm=1,nparmc(ict))
      call cuspinit4(1)

      iparm=0
      do 50 ict=1,nctype
        do 50 i=1,nparma(ict)
          iparm=iparm+1
   50     a(iwjasa(i,ict),ict)=a(iwjasa(i,ict),ict)-grad(indexi(iparm))
      do 60 i=1,nparmb(1)
        iparm=iparm+1
   60   b(iwjasb(i,1),1)=b(iwjasb(i,1),1)-grad(indexi(iparm))
      do 70 ict=1,nctype
        do 70 i=1,nparmc(ict)
          iparm=iparm+1
   70     c(iwjasc(i,ict),ict,1)=c(iwjasc(i,ict),ict,1)-grad(indexi(iparm))
      call cuspexact4(1)

      open(2,file='jas_new')
      if(nparma_read.gt.0) then
c       write(fmt,'(''(''i2,''f13.8,\'\' (a(iparmj),iparmj=1,nparma)\'\')'')') nparma_read
        write(fmt,'(''(''i2,''f13.8,a28)'')') nparma_read
       else
c       write(fmt,'(''(\'\' (a(iparmj),iparmj=1,nparma)\'\')'')')
        write(fmt,'(''(a28)'')')
      endif
      do 80 ict=1,nctype
   80   write(2,fmt) (a(i,ict),i=1,nparma_read),' (a(iparmj),iparmj=1,nparma)'

      if(nparmb_read.gt.0) then
c       write(fmt,'(''(''i2,''f13.8,\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')') nparmb_read
        write(fmt,'(''(''i2,''f13.8,a28)'')') nparmb_read
       else
c       write(fmt,'(''(\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')')
        write(fmt,'(''(a28)'')')
      endif
      write(2,fmt) (b(i,1),i=1,nparmb_read),' (b(iparmj),iparmj=1,nparmb)'

      if(nparmc_read.gt.0) then
c       write(fmt,'(''(''i2,''f13.8,\'\' (c(iparmj),iparmj=1,nparmc)\'\')'')') nparmc_read
        write(fmt,'(''(''i2,''f13.8,a28)'')') nparmc_read
       else
c       write(fmt,'(''(\'\' (c(iparmj),iparmj=1,nparmc)\'\')'')')
        write(fmt,'(''(a28)'')')
      endif
      do 90 ict=1,nctype
   90   write(2,fmt) (c(i,ict,1),i=1,nparmc_read),' (c(iparmj),iparmj=1,nparmc)'

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
      subroutine chlsky(a,n,np,ierr)
c chlsky: purpose: cholesky decomposition of a and determinant
c in: matrix a of order n stored with physical dimension np
c out: lower triangular matrix stored in lower portion of a
c note: lower triangular portion of original a is overwritten
      implicit real*8(a-h,o-z)
      dimension a(np,np)
      parameter (ZERO=0,ONE=1)

      diag_prod=1
      do j=1,n
         diag_prod=diag_prod*a(j,j)
      enddo

      det=1
      ierr=0
      do j=1,n
c        diag_prod=diag_prod*a(j,j)
         if(j.gt.1) then
            jm1=j-1
            do k=j,n
               sum=0
               do ip=1,jm1
                  sum=sum+a(k,ip)*a(j,ip)
               enddo
               a(k,j)=a(k,j)-sum
            enddo
         endif

         det=det*a(j,j)
         if(a(j,j).le.ZERO) then
           write(6,'(''Warning: '',i2,'' element of a is <0'',d9.2)') j,a(j,j)
           ierr=j
           return
         endif

         s=ONE/sqrt(a(j,j))
         do k=j,n
            a(k,j)=a(k,j)*s
         enddo
      enddo

      deta=1
      do j=1,n
         deta=deta*a(j,j)**2
      enddo

      write(6,'(''diag_prod,det,deta='',9d12.4)') diag_prod,det,deta

      return
      end
c-----------------------------------------------------------------------
      subroutine lxb(a,n,np,b)
c lxb: purpose: solve equation Lx=b, for lower triangular matrix L=a
c Golub and van Loan: algorithm 3.1.3
c in:  a = matrix of order n with physical dimension np
c      b = vector of order n
c out: b = solution of eqs.; overwites original b
      implicit real*8(a-h,o-z)
      dimension a(np,np),b(np)
      do j=1,n-1
         b(j)=b(j)/a(j,j)
         do k=j+1,n
            b(k)=b(k)-b(j)*a(k,j)
         enddo
      enddo
      b(n)=b(n)/a(n,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine uxb(a,n,np,b)
c uxb: purpose: solve equation Ux=b, for upper triangular matrix U=a
c Golub and van Loan: algorithm 3.1.4
c in:  a = matrix of order n with physical dimension np
c      b = vector of order n
c out: b = solution of eqs.; overwites original b
      implicit real*8(a-h,o-z)
      dimension a(np,np),b(np)
      do j=n,2,-1
         b(j)=b(j)/a(j,j)
         do k=1,j-1
            b(k)=b(k)-b(j)*a(k,j)
         enddo
      enddo
      b(1)=b(1)/a(1,1)
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
      write(6,'(''nordc,nctype,nterms='',9i5)') nordc,nctype,nterms

 100  continue

      return
      end

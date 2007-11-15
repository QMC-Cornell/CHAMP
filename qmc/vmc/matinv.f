      subroutine matinv(a,nsub,det)
c Written by Kevin Schmidt or from some standard library?

      implicit real*8(a-h,o-z)
      parameter (zero=0.d0,one=1.d0)

      include 'vmc.h'

c routine to calculate inverse and determinant of matrix a
c assumed to be dimensioned a(nsub,nsub).
c the matrix a is replaced by its inverse.

      dimension a(*)
      dimension ipivot(MELEC),atemp(MELEC)

      n=nsub
      do 5 i=1,n
    5   ipivot(i)=i

c initialize determinant
      determ=one

c loop through columns
      iclm=-n
      do 25 i=1,n
        iclm=iclm+n

c loop through rows and select row with largest element
        adiag=a(ipivot(i)+iclm)
        idiag=i
        do 15 k=i,n
          if(dabs(a(ipivot(k)+iclm)).gt.dabs(adiag)) then
            adiag=a(ipivot(k)+iclm)
            idiag=k
          endif
   15   continue

c interchange pointers if row different from
c original is selected and change sign of determinant because
c of interchange
        if(idiag.ne.i) then
          determ=-determ
          itemp=ipivot(i)
          ipivot(i)=ipivot(idiag)
          ipivot(idiag)=itemp
        endif

c update determinant
        determ=adiag*determ
        a(ipivot(i)+iclm)=one
        if(adiag.eq.0.d0) then
          write(6,'(''adiag=0 in matinv, slater matrix is singular., dimension='',i5,/,
     &    ''Possible fix: change renormaliz. rnorm in read_orb_pw.f if doing periodic system'',/,
     &    ''Another possible reason if this occurs during orbital optimization is that orbitals are linearly dep.'',/,
     &    ''which happens if norb=nbasis and cusps are imposed on the s orbitals'',/,
     &    ''The dimension gives you a clue as to where matinv was called from.'')') nsub
          stop 'adiag=0 in matinv'
        endif
        adiagi=one/adiag

c scale row by inverse diagonal
        call scal(n,adiagi,a(ipivot(i)),n)

c loop through other rows
c if not current row, then row reduce
        do 20 j=1,n
          if(j.ne.ipivot(i)) then
            t=-a(j+iclm)
            a(j+iclm)=zero
            call axpy(n,t,a(ipivot(i)),n,a(j),n)
          endif
   20   continue
   25 continue

c interchange elements to unpivot inverse matrix
c the following is equivalent to:
c      anew(i,ipivot(j))=aold(ipivot(i),j)
      jn=-n
      do 52 j=1,n
        jn=jn+n
        do 40 i=1,n
   40     atemp(i)=a(i+jn)
        do 50 i=1,n
   50     a(i+jn)=atemp(ipivot(i))
   52 continue

      do 55 j=1,n
   55   ipivot(j)=(ipivot(j)-1)*n

      do 90 i=1,n
        jn=-n
        do 70 j=1,n
          jn=jn+n
   70     atemp(j)=a(i+jn)
        do 80 j=1,n
   80     a(i+ipivot(j))=atemp(j)
   90 continue

      det=determ

      return
      end
c-----------------------------------------------------------------------
      subroutine axpy(n,a,x,nxskip,y,nyskip)
      implicit real*8(a-h,o-z)

c simplified blas routine to calculate y=a*x+y where
c x and y are arrays, and a is a scalar
c n      = number of elements in x and y to calculate
c nxskip = x stride
c nyskip = y stride

      dimension x(*),y(*)

      ix=1-nxskip
      iy=1-nyskip
      do 10 i=1,n
        ix=ix+nxskip
        iy=iy+nyskip
   10   y(iy)=a*x(ix)+y(iy)
      return
      end
c-----------------------------------------------------------------------
      subroutine scal(n,scalor,x,nskip)
      implicit real*8(a-h,o-z)

c simplified blas routine to scale x by scalor
c n     = number of elements in x to be scaled
c nskip = stride

      dimension x(*)

      ix=1-nskip
      do 10 i=1,n
        ix=ix+nskip
   10   x(ix)=scalor*x(ix)
      return
      end

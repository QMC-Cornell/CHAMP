      subroutine matinv(a,nsub,det)
! Written by Kevin Schmidt or from some standard library?

      implicit real*8(a-h,o-z)
      parameter (zero=0.d0,one=1.d0)

! routine to calculate inverse and determinant of matrix a
! assumed to be dimensioned a(nsub,nsub).
! the matrix a is replaced by its inverse.

      dimension a(nsub**2)
      dimension ipivot(nsub),atemp(nsub)

!     rmaxval=maxval(a,nsub**2)
!     write(6,'(''rmaxval='',es12.4)') rmaxval

      n=nsub
!      write(6,*) 'Inverting this matrix:'
!      do irow=0,(n-1)
!         write(6,'(100f6.2)') (a(n*irow+k),k=1,n)
!      enddo
      do 5 i=1,n
    5   ipivot(i)=i

! initialize determinant
      determ=one

! loop through columns
      iclm=-n
      do 25 i=1,n
        iclm=iclm+n

! loop through rows and select row with largest element
        adiag=a(ipivot(i)+iclm)
        idiag=i
        do 15 k=i,n
          if(dabs(a(ipivot(k)+iclm)).gt.dabs(adiag)) then
            adiag=a(ipivot(k)+iclm)
            idiag=k
          endif
   15   continue

! interchange pointers if row different from
! original is selected and change sign of determinant because
! of interchange
        if(idiag.ne.i) then
          determ=-determ
          itemp=ipivot(i)
          ipivot(i)=ipivot(idiag)
          ipivot(idiag)=itemp
        endif

! update determinant
        determ=adiag*determ
        a(ipivot(i)+iclm)=one
!       write(6,*)
!       do irow=0,(n-1)
!           write(6,'(100es9.1)') (a(irow*n+k),k=1,n)
!       enddo
        if(adiag.eq.0.d0) then
          write(6,'(''in matinv, n='',i5)') n
          call systemflush(6)
          write(6,'(''adiag=0 in matinv, slater matrix is singular., dimension='',i5,/,
     &    ''Possible fix: change renormaliz. rnorm in read_orb_pw.f if doing periodic system.'',/,
     &    ''Another possible reason if this occurs during orbital optimization is that orbitals are linearly dep.'',/,
     &    ''which happens if norb=nbasis and cusps are imposed on the s orbitals.'',/,
     &    ''Yet another possible reason is failure to do a "make clean" between "make" and "make mpi" or vice versa.'',/,
     &    ''The dimension gives you a clue as to where matinv was called from.'')') nsub
          do irow=0,(n-1)
              write(6,'(100es9.1)') (a(irow*n+k),k=1,n)
          enddo
          call systemflush(6)
          stop 'adiag=0 in matinv'
        endif
        adiagi=one/adiag

! scale row by inverse diagonal
        call scal(n,adiagi,a(ipivot(i)),n)

! loop through other rows
! if not current row, then row reduce
        do 20 j=1,n
          if(j.ne.ipivot(i)) then
            t=-a(j+iclm)
            a(j+iclm)=zero
            call axpy(n,t,a(ipivot(i)),n,a(j),n)
          endif
   20   continue
   25 continue

! interchange elements to unpivot inverse matrix
! the following is equivalent to:
!      anew(i,ipivot(j))=aold(ipivot(i),j)
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
!-----------------------------------------------------------------------
      subroutine axpy(n,a,x,nxskip,y,nyskip)
      implicit real*8(a-h,o-z)

! simplified blas routine to calculate y=a*x+y where
! x and y are arrays, and a is a scalar
! n      = number of elements in x and y to calculate
! nxskip = x stride
! nyskip = y stride

      dimension x(*),y(*)

      ix=1-nxskip
      iy=1-nyskip
      do 10 i=1,n
        ix=ix+nxskip
        iy=iy+nyskip
   10   y(iy)=a*x(ix)+y(iy)
      return
      end
!-----------------------------------------------------------------------
      subroutine scal(n,scalar,x,nskip)
      implicit real*8(a-h,o-z)

! simplified blas routine to scale x by scalar
! n     = number of elements in x to be scaled
! nskip = stride

      dimension x(*)

      ix=1-nskip
      do 10 i=1,n
        ix=ix+nskip
   10   x(ix)=scalar*x(ix)
      return
      end

c Simple test routine to check quad_min
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)
      parameter(MPTS=100)
      dimension e(MPTS),add_diag(MPTS)

      read(5,*) npts,eig_min,eig_max
      if(npts.gt.MPTS) stop 'npts>MPTS'
      do 10 k=1,npts
   10   read(5,*) add_diag(k),e(k)

      call quad_min(e,add_diag,npts,eig_min,eig_max)
      stop
      end
c-----------------------------------------------------------------------
      subroutine quad_min(e,add_diag,npts,eig_min,eig_max)
c Find the value of add_diag at which e is minimum by fitting to a 3-parameter function.
c First tried a parabola in log10(add_diag) and that worked well, but a modification
c of this worked even better.
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)
      parameter(MFUNC=3)
      dimension e(npts),add_diag(npts),add_diag_log(npts),a(MFUNC,MFUNC),b(MFUNC)

      nfunc=MFUNC

      if(eig_max.lt.1.d99) then
        dlog10_eig_max=dlog10(eig_max)
       else
        write(6,'(''dlog10_eig_max reset to 1.d300'')')
        dlog10_eig_max=1.d300
      endif

      do 5 k=1,npts
c   5   add_diag_log(k)=dlog10(add_diag(k))
        if(add_diag(k).lt.1.d99) then
          tmp=dlog10(add_diag(k)+max(eig_min,0.d0))
          add_diag_log(k)=tmp/(1+tmp/dlog10_eig_max)
         else
          add_diag_log(k)=dlog10(eig_max)
        endif
    5 continue

      do 30 i=1,nfunc
        b(i)=0
        do 10 k=1,npts
   10     b(i)=b(i)+e(k)*add_diag_log(k)**(i-1)
        do 30 j=1,i
          a(i,j)=0
          do 20 k=1,npts
   20       a(i,j)=a(i,j)+add_diag_log(k)**(i+j-2)
   30     a(j,i)=a(i,j)

      write(6,*) b
      write(6,*)
      write(6,'(3f15.6)') a

c Do cholesky decomposition
      call chlsky(a,nfunc,MFUNC,ierr)
      if(ierr.ne.0) stop 'ierr ne 0 in chlsky'

c Symmetrize decomposed matrix (needs to be done before calling uxb
c or need to modify uxb)
      do 40 i=1,nfunc
        do 40 j=i+1,nfunc
   40     a(i,j)=a(j,i)

c Solve linear equations
      call lxb(a,nfunc,MFUNC,b)
      call uxb(a,nfunc,MFUNC,b)

      write(6,'(''adiag_log,e'',2f10.5)') (add_diag_log(k),e(k),k=1,npts)
      write(6,'(''b='',3d12.4)') b

      e_min=9.d99
      rms=0
      do 50 k=1,npts
        ee=b(1)+b(2)*add_diag_log(k)+b(3)*add_diag_log(k)**2
        if(e(k).lt.e_min) then
          k_min=k
          e_min=e(k)
        endif
        rms=rms+(ee-e(k))**2
   50   write(6,'(9g14.6)') add_diag(k),ee,e(k),ee-e(k)
      rms=dsqrt(rms/npts)
      write(6,'(''rms error is'',d12.4)') rms

      if(b(3).gt.0) then
        add_diag_log_min=-0.5d0*b(2)/b(3)
       else
c       write(6,'(''b='',3d12.4)') b
        if(k_min.eq.1) then
          add_diag_log_min=2*add_diag_log(1)-add_diag_log(3)
         elseif(k_min.eq.3) then
          add_diag_log_min=5*add_diag_log(3)-4*add_diag_log(1)
        endif
      endif
c     add_diag_min=10**add_diag_log_min
      add_diag_min=10**(add_diag_log_min/(1-add_diag_log_min/dlog10_eig_max))-max(eig_min,0.d0)
      e_min=b(1)+b(2)*add_diag_log_min+b(3)*add_diag_log_min**2
      write(6,'(''add_diag_min,e_min='',2g14.6)') add_diag_min,e_min

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

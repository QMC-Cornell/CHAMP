cDate of last modification: Tue Apr 19 09:54:33 EDT 1994
      subroutine cholesky_piv(A,ipiv,iperm,n,np,ir,eps)
c cholesky_piv: purpose Cholesky with symmetric pivoting
c A	= matrix of order n of physical dimension np
c	  on input: matrix to be decomposed,
c	  on output: upper and lower triangular 1:ir by 1:ir Cholesky factors
c ipiv	= the Cholesky decomposition permuted states i and ipiv(i) at step i
c iperm	= the input matrix A_{in}(i,j) can be reconstituted from the output
c	  matrix A_{out} as follows:
c	    A'(i,j)=sum_{k=1}^{min(i,j,ir)} A(i,k)*A(k,j)
c           A(i,j)=A'(iperm(i),iperm(j))
c         To use for the solution of A.x=b
c         c'(iperm(i))=c(i) for all i
c         call lxb and then uxb with c'
c         solution: c(i)=c'(iperm(i))
c         Note this solution reproduces only the non-singular components of b
c         if the matrix A is rank defincient.
c ir	= apparent rank of the matrix
c eps	= if at step i of the algorithm A(i,i)/|A| < eps, where A(i,i) is the current
c	  not the original matrix, the matrix is considered
c	  rank deficient; |A|=max_i A(i,i) for the input matrix.
      implicit real*8(a-h,o-z)
      dimension A(np,np),ipiv(n),iperm(n)
      ir=0
      amax=0
      do i=1,n
        iperm(i)=i
        amax=max(amax,A(i,i))
      enddo
      aeps=sqrt(amax)*eps
      do k=1,n
        imax=k
        amax=A(k,k)
        do l=k+1,n
          if(A(l,l).gt.amax) then
            imax=l
            amax=A(l,l)
          endif
        enddo
        if(amax.gt.aeps) then
          ir=ir+1
          ipiv(k)=imax
          itemp=iperm(k)
          iperm(k)=iperm(imax)
          iperm(imax)=itemp
          if(imax.ne.k) then
            do i=1,n
              temp=A(i,k)
              A(i,k)=A(i,imax)
              A(i,imax)=temp
            enddo
            do i=1,n
              temp=A(k,i)
              A(k,i)=A(imax,i)
              A(imax,i)=temp
            enddo
          endif
          A(k,k)=sqrt(A(k,k))
          temp=A(k,k)
          do i=k+1,n
            A(i,k)=A(i,k)/temp
            A(k,i)=A(i,k)
          enddo
          do j=k+1,n
            do i=j,n
              A(i,j)=A(i,j)-A(i,k)*A(j,k)
              A(j,i)=A(i,j)
            enddo
          enddo
        endif
      enddo
      do i=1,n
      ipiv(iperm(i))=i
      enddo
      return
      end

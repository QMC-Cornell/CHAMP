c $Log: svd_gaus_test.F,v $
c Revision 1.1.1.1  2006/10/08 19:38:17  toulouse
c champ
c
c Revision 1.4  2003/07/14 17:41:49  nigh
c *** empty log message ***
c
c Revision 1.3  2002/12/21 03:58:57  nigh
c replaced _wrap by _qa
c
c Revision 1.2  2002/12/20 19:00:32  nigh
c Cyrus's changes: some exact derivatives
c
c Revision 1.1  2002/03/04 16:11:12  nigh
c added LAPACK precompiler option
c
c Revision 1.7  2002/01/25 16:19:42  nigh
c output changes
c If next line is not commented use LAPACK else use LINPACK
c #define LAPACK
      subroutine svd_gaus_test(A,B,n,np,ajac,m,mp,eps,dseed,newA,r,
     &  problog,
     &  pronly,T,parmar,s,u,v)
c mvgaus: purpose: sample from multivariate gaussian/calculate probablility
c         density prop. to exp[-(r.A.r/2+B.r)/T], where A'=sqrt(A**2)
c in:
c A	= symmetric matrix of dimension n and physical dimension np
c         A is used only to check consistency of Hessian and Jacobian
c B	= vector of order n
c ajac	= m x n jacobian of leading physical dimension mp; ajac is
c	  destroyed
c eps	= singular value s_i is cutoff from below by max(s_i,max_k(s_k)*eps)
c dseed = seed from random number generator
c newA	= true for the first call with matrix A
c         false for repeat call with same A
c r	= n-vector or which probability will be calculated, if
c          pronly = .true.
c parmar= shift to be applied so that transpose(ajac).ajac + parmar = A
c out:
c B	= inverse(A).B
c r	= random n-vector sampled from gaussian
c problog=log( sqrt{det[A/T]/(2 pi)^n} exp(-B.A^(-1).B/2T) exp[-(r.A.r/2+B.r)/T)
c In this normalization factor:
c (a) the factors of 1/sqrt(2 pi) are the usual normalization
c     factors of the N(0,1) distribution.
c (b) Sqrt(det[A]) is the Jacobian due to the transformation to
c     principal axes.
c (c) The exponential factor comes about by completing the square and
c     shifting the coordinates.
c parmar= input value is changed to satisfy
c	  parmar * DBL_EPSILON < max_i(s_i)^2
c	  and parmar > DBL_EPSILON * min_i(s_i)^2,
c	  where the s_i are the singular values of the Hamiltonian.
c decompose ajac(m x n) = u(m x n).diag[s(1),...,s(n)].vt(n x n)
c In general ajac(m x n) = u(m x p).diag[s(1),...,s(p)].vt(p x n), p=min(m,n).
c If m < n, u is square, while ajac and vt have the same shape.
c In this program vt is square, while u and ajac have the same shape, which
c works so long as n<=m (# of parameters <= # of data).

      implicit real*8(a-h,o-z)

      include 'parameters.h'
      parameter (TWO_PI_INV=0.159154943091895336d0)
      dimension A(np,n),B(n),r(n),ajac(mp,n)
#ifdef LAPACK
      dimension s(n),e(n),vt(n,n),y(n)
#else
      dimension s(m+1),e(n),u(n,m),v(n,n),y(m)
#endif
      logical newA,pronly,zeroT
      save amins,amaxs,n_act
c ipr_com determines amount of output
      logical called_by_qa
      common /quenchsim_pr/ ipr_com,called_by_qa
      logical first_call
      save first_call
      data first_call/.true./

      if(m .lt. n) then
c need to fix routine to be able to deal with fewer data than parameters
        write(6,*) 'svd_gaus_test: m .lt. n',m,n
        stop 'svd_gaus_test: m .lt. n'
      endif

      if(T.eq.ZERO) then
c T0 has no meaning here it's defined to avoid zero devides
        zeroT=.true.
        T0=ONE
        T0log=ZERO
      else
        zeroT=.false.
        T0=T
        T0log=log(T0)
        T0sqrt=sqrt(T0)
      endif

      if(pronly) then
c check validity of Hessian
        if(newA) then
          aa=ZERO
          do k=1,m
            aa=aa+ajac(k,1)*ajac(k,1)
          enddo
c These lines commented out after talking to Peter because it is OK to use approx. Hessian
c         if(abs(aa+parmar-A(1,1)).gt.abs(A(1,1))*1e-12*m) then
c           stop 'svd_gaus: A<>Transpose(Jacobian).Jacobian'
c         endif
        endif

        s2=ZERO
        do i=1,n
          ri=r(i)
          s2=s2-B(i)*ri
          do j=1,i-1
            s2=s2-ri*A(i,j)*r(j)
          enddo
          s2=s2-HALF*(A(i,i)+parmar)*ri**2
        enddo
      endif

      if(newA) then
#ifdef LAPACK
        if(first_call) write(6,*) 'svd_gaus: LAPACK'
        lwork=-1
        call dgesvd('n','a',m,n,ajac,mp,s,u,m,vt,n,will_work,lwork,info)
        lwork=will_work
        call dgesvd_qa('n','a',m,n,ajac,mp,s,u,m,vt,n,lwork,info)
#else
        if(first_call) write(6,*) 'svd_gaus: LINPACK'
        job=21
        job=01
        call dsvdc(ajac,mp,m,n,s,e,u,m,v,n,y,job,info)
        if(info.ne.0) then
          write(6,'(''svd_gaus_test: Warning info.ne.0 in dsvdc'',i5)')
     &    info
          stop 'svd_gaus_test: Warning info.ne.0 in dsvdc'
        endif
#endif
        first_call=.false.
        amins=s(1)
        amaxs=s(1)
        do i=2,n
          amins=min(amins,s(i))
          amaxs=max(amaxs,s(i))
        enddo
        if(ipr_com.ge.1) then
          write(6,'(''svd_gaus: amins,amaxs '',2d9.2)') amins,amaxs
          write(6,'(''sing. val. ratios:'',8g9.2)')
     &      (s(i)/amaxs,i=1,n)
        endif
        cutoff=amaxs*dbl_epsilon
        n_act=n
        do i=1,n
c         s(i)=max(s(i),cutoff)
          if(s(i).lt.cutoff) then
            s(i)=ZERO
            n_act=n_act-1
          endif
        enddo
      endif
c keep parmar within bounds
      parmar=min(parmar,TENTH*amaxs**2/dbl_epsilon,sqrt(dbl_max))
      parmar=max(parmar,TEN*dbl_epsilon*amins**2,sqrt(dbl_min))
c compute log of determinant
      d1=ZERO
      do i=1,n
        if(s(i).gt.ZERO) d1=d1+log(TWO_PI_INV*(s(i)**2+parmar))
      enddo

      if(pronly) then
        b_ainv_b=ZERO
        do i=1,n
          sum=0
          do l=1,n
#ifdef LAPACK
            sum=sum+vt(i,l)*B(l)
#else
            sum=sum+v(l,i)*B(l)
#endif
          enddo
          if(s(i).gt.ZERO) b_ainv_b=b_ainv_b+sum**2/(s(i)**2+parmar)
        enddo
        problog=(s2-HALF*b_ainv_b)/T0+HALF*(d1-n_act*T0log)
        return
      else
c sample from Gaussian
        s2=ZERO
        do k=1,n
          if(s(k).gt.ZERO) then
            if(zeroT) then
              sg=ZERO
            else
              sg=T0sqrt*gaushn()
            endif
            s2=s2-sg**2
            sg=sg/sqrt(s(k)**2+parmar)
          else
            sg=ZERO
          endif
          tmp=ZERO
          do l=1,n
#ifdef LAPACK
            tmp=tmp+vt(k,l)*B(l)
#else
            tmp=tmp+v(l,k)*B(l)
#endif
          enddo
          if(s(k).gt.ZERO) then
            e(k)=tmp/(s(k)**2+parmar)
          else
            e(k)=ZERO
          endif
          y(k)=sg-e(k)
        enddo
        do i=1,n
          sum1=0
          sum2=0
          do k=1,n
#ifdef LAPACK
            sum1=sum1+vt(k,i)*y(k)
            sum2=sum2+vt(k,i)*e(k)
#else
            sum1=sum1+v(i,k)*y(k)
            sum2=sum2+v(i,k)*e(k)
#endif
          enddo
          r(i)=sum1
          B(i)=sum2
        enddo
        problog=HALF*(s2/T0-n_act*T0log+d1)
        return
      endif
      end

#ifdef LAPACK
      subroutine dgesvd_qa(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,lwork,
     &  info)
c routine serves to make array "work" dynamic
      implicit real*8(a-h,o-z)
      character*1 jobu,jobvt
      dimension a(lda,n),s(n),vt(n,n)
      dimension work(lwork)
      call dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
      return
      end
#endif

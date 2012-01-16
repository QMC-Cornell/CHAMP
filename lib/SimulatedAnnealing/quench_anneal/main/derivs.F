c $Log: derivs.F,v $
c Revision 1.1.1.1  2006/10/08 19:38:17  toulouse
c champ
c
c Revision 1.1  2003/07/14 17:41:37  nigh
c *** empty log message ***
c
c Revision 1.6  2002/12/20 19:00:32  nigh
c Cyrus's changes: some exact derivatives
c
c Revision 1.5  2002/01/18 18:39:17  nigh
c introduced two-dimensional diff in func_mpi and derivs
c
c
cDate of last modification: Tue Jan  5 12:57:02 EST 1999
      subroutine derivs(A,B,C,parm,ndata,nparm,nanalytic,diff,func,
     &  jacobian,aiimax,pmarquardt,central,use_sav,update_B,pdif,ajac,
     &  eps_diff)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c of an objective function which is a sum of squares.
c A     = approximate hessian of chi^2 a la Marquardt-Levenberg
c B     = 1st derivative
c C 	= scaling vector
c parm  = parameters
c ndata = number of data pts
c nparm = number of parameters
c diff  = residues of the model at each of the data pnts; must correspond EXACTLY
c to parm(i), unless update_B=.not.use_sav=.true.)
c inconsistencies due to progamming errors can be checked with DEBUG_DERIVS
c#define DEBUG_DERIVS
c func  = sum of squares (chisq)
c jacobian = if nanalytic>0 this routine is called to evaluate the
c            jacobian
c aiimax= max diagonal element of A (may be used for scaling)
c pmarquardt= marquardt parameter
c central= .true. for central differences, .false. for forward
c use_sav = .true. no new derivatives are calculated or rank-1 update is used
c                  scaling vector is recalculated
c	    .false. new derivative are calculated
c update_B= use rank-1 update to compute B
c pdif	= reverse difference in parameters
c ajac	= Jacobian: ajac(k,l) = df_k/dp_l (f are residues; p are
c	  parameters)
c eps_diff= relative numerical accuracy of residues, used if greater
c           than DBL_EPSILON for increments in numerical differentiation
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include '../include/parameters.h'
      external func,jacobian
      common /ichange/ ichange
      dimension a(nparm,nparm),b(nparm),C(nparm)
      common/quenchsim_diff/ndim1_diff,ndim2_diff
      dimension parm(nparm),diff(ndata,ndim2_diff),
     &  diffm1(ndata,ndim2_diff),diffp1(ndata,ndim2_diff),
     &  diffm3(ndata,ndim2_diff),diffp3(ndata,ndim2_diff),
     &  ajac(ndata,nparm),pdif(nparm)
      logical central,use_sav,update_B,cntr
c common /quenchsim_pr/ to communicate amount of output
      logical called_by_qa
      common /quenchsim_pr/ ipr_com,called_by_qa
      parameter (FORTY_EIGHTH=ONE/48)
      logical central_derivs_qs
c allow user to switch to central 3-pnt or 5-pnt finite diffs (central point
c included in count, but only used for forward diffs.)
      common /quenchsim_derivs/central_derivs_qs,num_deriv_points_qs
#ifdef DEBUG_DERIVS
c consistency check for debug only; may double func execution time 
c checks to see if parameters were not changed after diffs were calculated
      dimension diff_debug(ndata,ndim2_diff)
      if(update_B .eq. .false. .and. use_sav .eq. .false.) then
        chi2_debug=func_mpi(func,ndata,nparm,parm,diff_debug,iflag)
        do i=1,ndata
          if(diff_debug(i,1) .ne. diff(i,1)) then
            write(6,*) 'derivs: diff_debug(i) .ne. diff(i)',i,diff_debug(i,1),diff(i,1)
            stop 'derivs: diff_debug .ne. diff'
          endif
        enddo
      endif
#endif
c derivatives
      if(.not.use_sav) then
        nnumerical=nparm-nanalytic
        if(nanalytic.gt.0) then
          call jacobian(ndata,nparm,nanalytic,parm,ajac)
        endif
        if(nnumerical.gt.0) then
          if(num_deriv_points_qs .ne. 3 .and.
     &      num_deriv_points_qs .ne. 5) then
c           write(6,*) 'derivs: incorrect value of num_deriv_points_qs',
c    &        num_deriv_points_qs
            write(6,'(''derivs: incorrect value of num_deriv_points_qs''
     &      ,i30)') num_deriv_points_qs
            stop 'derivs: incorrect value of num_deriv_points_qs'
          endif
          cntr=central .or. central_derivs_qs
          epsder1f=sqrt(max(eps_diff,DBL_EPSILON))
          epsder1=max(eps_diff,DBL_EPSILON)**THIRD
          epsder2=max(eps_diff,DBL_EPSILON)**FIFTH
          do i=1,nnumerical
            parmsav=parm(i)
            if(cntr) then
              if(num_deriv_points_qs .eq. 3) then
                eps=max(epsder1,epsder1*abs(parm(i)))
              else ! assume that 5  points are required
                eps=max(epsder2,epsder2*abs(parm(i)))
              endif
            else
              eps=max(epsder1f,epsder1f*abs(parm(i)))
            endif
            parm(i)=parmsav+eps
c compute the exact difference eps (see Numerical Recipes, Chap 5.7)
            eps=difference(parm(i),parmsav)
            factor=ONE/eps
            half_factor=HALF*factor
            forty8th_factor=FORTY_EIGHTH*factor
            f=func_mpi(func,ndata,nparm,parm,diffp1,iflag)
            ichange=0
c           write(6,*) 'in derivs',ichange
            if(cntr) then
              parm(i)=parmsav-eps
              f=func_mpi(func,ndata,nparm,parm,diffm1,iflag)
              do k=1,ndata
                ajac(k,i)=half_factor*(diffp1(k,1)-diffm1(k,1))
              enddo
              if(num_deriv_points_qs .gt. 3) then
                parm(i)=parmsav+THREE*eps
                f=func_mpi(func,ndata,nparm,parm,diffp3,iflag)
                parm(i)=parmsav-THREE*eps
                f=func_mpi(func,ndata,nparm,parm,diffm3,iflag)
                do k=1,ndata
                  tmp=ajac(k,i)
                  ajac(k,i)=
     &              forty8th_factor*(diffm3(k,1)-diffp3(k,1)+
     &              THREE*(diffp1(k,1)-diffm1(k,1)))+ajac(k,i)
                enddo
              endif
            else
              do k=1,ndata
                ajac(k,i)=factor*(diffp1(k,1)-diff(k,1))
              enddo
            endif
            parm(i)=parmsav
          enddo
        endif

c gradient
        do i=1,nparm
          b(i)=ZERO
          do k=1,ndata
              b(i)=b(i)+diff(k,1)*ajac(k,i)
          enddo
          b(i)=TWO*b(i)
        enddo
c approximate Hessian
        do i=1,nparm
          do j=1,nparm
            a(i,j)=ZERO
            do k=1,ndata
              a(i,j)=a(i,j)+ajac(k,i)*ajac(k,j)
            enddo
            a(i,j)=TWO*a(i,j)
          enddo
        enddo
c Jacobian
        do iparm=1,nparm
          do idata=1,ndata
            ajac(idata,iparm)=sqrt(TWO)*ajac(idata,iparm)
          enddo
        enddo
      endif

      if(ipr_com.eq.999) then
        write(6,'(''derivs: Hessian'')')
        do ipr=1,nparm
          write(6,'((5(2(i3,1x),g9.2)))')
     &      (ipr,jpr,a(ipr,jpr),jpr=1,nparm)
        enddo
      endif

c rank-1 update of B
      if(update_B) then
        do i=1,nparm
          tmp=B(i)
          do j=1,nparm
            tmp=tmp-A(i,j)*pdif(j)
          enddo
          B(i)=tmp
        enddo
      endif

c rescaling vector
      aiimax=a(1,1)
      aiimin=a(1,1)
      do iparm=2,nparm
        aiimax=max(aiimax,a(iparm,iparm))
        aiimin=min(aiimin,a(iparm,iparm))
      enddo

      f=sqrt(aiimax*pmarquardt)
      if(ipr_com.ge.2) then
        write(6,'(''derivs: aiimax,aiimin,pmarquardt'',9(1x,g12.6))')aiimax,aiimin,pmarquardt
c       write(6,'(''derivs: det for 2 parm prob'',9(1x,g12.6))') a(1,1)*a(2,2)-a(1,2)*a(2,1)
      endif
      do iparm=1,nparm
c ?? This will scale the effective hessian so that all diagonal elements are
c ?? one.  This may be not such a great idea: it may interfere with a pivoting
c ?? strategy in the Cholesky decomposition.
c ?? choose either next line or the five following it
c       c(iparm)=sqrt(aiimax*pmarquardt+a(iparm,iparm))
        c(iparm)=ZERO
        do idata=1,ndata
          c(iparm)=c(iparm)+abs(ajac(idata,iparm))
        enddo
        c(iparm)=c(iparm)+f
      enddo

c     call fullrankQ(A,nparm)
      return
      end
c-----------------------------------------------------------------------

      block data quenchsim_derivs_block
      logical central_derivs_qs
      common / quenchsim_derivs / central_derivs_qs,
     &   num_deriv_points_qs
      data central_derivs_qs/.false./, num_deriv_points_qs/3/
      end

c-----------------------------------------------------------------------
      function difference(a,b)
c force the compiler to write results to memory.  If it has higher
c internal CPU accuracy we want eps_diff to be such that
c (par+eps_diff)-par=eps_diff EXACTLY without roundoff error.
c See for Numerical Recipes, Chap 5.7
      implicit real*8(a-h,o-z)
      difference=a-b
      return
      end

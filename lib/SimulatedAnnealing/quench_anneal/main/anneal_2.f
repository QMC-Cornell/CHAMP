c $Log: anneal_2.f,v $
c Revision 1.10  2002/12/20 19:00:32  nigh
c Cyrus's changes: some exact derivatives
c
c Revision 1.9  2002/03/04 14:21:32  nigh
c new version with LAPACK and LINPACK; see svd_gaus_test
c
c Revision 1.8  2002/01/06 20:24:46  nigh
c Cyrus's new version with unused variables removed
c
c Revision 1.7  2002/01/04 15:03:54  nigh
c fixed log lines
c
c Revision 1.6  2002/01/04 15:00:36  nigh
c made character variable mesg consistent with calls
c
c Revision 1.5  2002/01/02 15:50:01  nigh
c added after making changes in svd_gaus_test etc. to save s,u,v
c
      subroutine step(A,B,Asav,Bsav,ndata,nparm,nanalytic,p,pdif,diff
     &,pmarquardt,dseed,newA,pronly,prob,func,jacobian,T,epsg,epsp
     &,converg,mesg,central,use_sav,update_B,cholesky,pratio,rot_wt,
     &eps_diff,aiimax,ajac,s_svd,u_svd,v_svd)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ?? this routine has not been tested for newA=.false. for which it's not
c ?? being used currrently
c A	= input: approximate hessian if newA=.true.
c	  output: modified for use on next with newA=.false.
c B	= input: gradient if newA=.true.
c	  output: modified for use on next with newA=.false.
c Bsave	= output: gradient if newA=.false., not used otherwise
c p	= parameters passed to func
c pdif	= reverse change in parameters
c diff	= residues computed by func
c pmarquardt = marquardt parameter is changed if too big or too small
c dseed	= seed used by random number generator
c newA	= Cholesky decomposition is recomputed only for newA=.true.
c pronly= compute log. probability for given pdif
c prob	= above probability
c func  = function to evaluate residues
c jacobian = if nanalytic>0 this routine is called to evaluate the
c            jacobian
c T	= temperature determines variance of noise (see lugaus_new)
c epsg	= convergence if gnorm < epsg,
c 	  where gnorm = 1-norm of scaled gradient with components B(i)/A(i,i)
c epsp	= convergence if for all i:
c	  pnorm <= epsp abs(1-max(1,pnorm/pnorm'))
c	  where pnorm = 1-norm of the vector with components
c	  pdif(i)/max(1,|p(i)|) and pnorm' the value of this quantity
c	  at the previous iteration
c converg = .true. if convergence has occurred
c mesg        = "gr" and/or "pr" depending on convergence condition
c central = central or forward differences (see derivs). To be inialized by calling
c	    routine; computed by this one ever after.
c use_sav = true if old derivatives etc. are to be used as saved in variables
c	    with names ending in "sav"
c update_B = use A to evaluate B at a nearby point
c pnorm1 = previous value of pnorm
c cholesky      = .true. use LU decomposition; .false. use singular
c                 value decomposition
c pratio	= ratio of norms: pdif's with zero Marquardt parameter
c		  over non-zero value (beware of cutoff in svd!)
c	 	  This ratio is only calculated for cholesky=.false.
c eps_diff= relative numerical accuracy of residues, used if greater
c           than DBL_EPSILON for increments in numerical differentiation
c aiimax	= maximum diagonal element of the Hessian
c ajac          = jacobian matrix (had to be saved between calls)
c s_svd,u_svd,v_svd = arrays saved for svd_gaus_test
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include 'parameters.h'
      logical newA,pronly,converg,central,use_sav,update_B,
     &  cholesky
c ?? This routine needs to be cleaned up: the *sav and *_s are
c ?? not required simultaneously etc.
      common/quenchsim_diff/ndim1_diff,ndim2_diff
      dimension A(nparm,nparm),Asav(nparm,nparm),B(nparm),Bsav(nparm)
     &  ,C(nparm),p(nparm),pdif(nparm),diff(ndata,ndim2_diff)
     &  ,ipermu(nparm),A_s(nparm,nparm),B_s(nparm),C_s(nparm)
     &  ,pdif_s(nparm),ajac(ndata,nparm),ajac_s(ndata,nparm)
     &  ,tmp(nparm),psav(nparm)
      dimension s_svd(ndata+1),u_svd(nparm,ndata),v_svd(nparm,nparm)
      external func,jacobian
      character*(*) mesg
      character*3 grad_msg,parm_msg

c common /quenchsim_pr/ to communicate amount of output
      logical called_by_qa
      common /quenchsim_pr/ ipr_com,called_by_qa

      converg=.false.
c     use_sav=.false.
      if(newA) then
        if(.not.use_sav) then
c compute and save derivatives
          if(update_B) stop 'step: update_B=.not.use_sav=.true. ?'
          call derivs(A,B,C,p,ndata,nparm,nanalytic,diff,func,jacobian,
     &    aiimax,pmarquardt,central,use_sav,update_B,pdif,ajac,eps_diff)
          call rank(A,Asav,nparm,ipermu,nparm_s)
          do iparm=1,nparm
            Bsav(iparm)=B(iparm)
            do jparm=1,iparm
              Asav(iparm,jparm)=A(iparm,jparm)
              Asav(jparm,iparm)=A(iparm,jparm)
            enddo
          enddo
        else
c restore derivatives aiimax and C
          do iparm=1,nparm
            B(iparm)=Bsav(iparm)
            do jparm=1,iparm
              A(iparm,jparm)=Asav(iparm,jparm)
              A(jparm,iparm)=Asav(iparm,jparm)
            enddo
          enddo
          call derivs(A,B,C,p,ndata,nparm,nanalytic,diff,func,jacobian,
     &    aiimax,pmarquardt,central,use_sav,update_B,pdif,ajac,eps_diff)
        endif
        call select(A,B,ajac,ndata,C,nparm,pdif,A_s,B_s,ajac_s,C_s,
     &    pdif_s,ipermu,nparm_s)
c ?? to have any effect on the algorithm the call to scale should be in
c ?? the next line, but our way of rescaling is probably not good.
c ?? Maybe just a rescaling based on the original point would suffice.
c ?? See the discussion More (The LM algorithm: implemetation and
c ?? theory)
c       call scale(A_s,B_s,ajac_s,ndata,C_s,nparm_s,pdif_s,pronly)
        if(cholesky)
     &  call marquardt(A_s,nparm_s,nparm,pmarquardt,aiimax,C_s)
c       write(6,*) 'A_s'
c       do i=1,nparm_s
c         write(6,*) (A_s(i,j),j=1,nparm_s)
c       enddo
c        call marquardt(A,nparm,nparm,pmarquardt,aiimax)
c        write(6,*) 'A'
c        do i=1,nparm
c          write(6,*) (A_s(i,j),j=1,nparm)
c        enddo
      else
        stop 'step: entering uncharted territory'
      endif

      if(pronly) then
        do iparm=1,nparm
          psav(iparm)=pdif(iparm)
        enddo
      endif

c check convergence for gradient
c ?? Add something to denom in case a diagonal element of A is 0,
c as is the case for one of minpack test cases.
      gnorm=ZERO
      do iparm=1,nparm
        gnorm=gnorm+abs(B(iparm)
c    &    /A(iparm,iparm))
     &    /(A(iparm,iparm)+min(ONE,pmarquardt*aiimax)))
      enddo
c      write(6,*) 'gnorm:',gnorm
      if(gnorm.le.epsg*nparm) then
        converg=.true.
        grad_msg='gr '
      else
        grad_msg=' '
      endif

c ?? It remains to be seen that scaling helps the numerical stability ...
c ?? In particular if we switch to Cholesky with pivoting, scaling may
c ?? be counter productive.
c ?? the calls to scale and unscale may be removed simultaneously
c ?? Scaling in this place makes shifting of the singular values in
c ?? svd_gaus inequivalent to the effect of adding a multiple of
c ?? unity to the Hessian
c      call scale(A_s,B_s,ajac_s,ndata,C_s,nparm_s,pdif_s,pronly)
c      call scale(A,B,C,nparm,pdif,pronly)

c Note that if the call to marquardt is here rather than ~10 lines up,
c the tendency to go to a saddle point at high temperatures
c (i.e. temperatures where diffusion > drift) is greatly enhanced.
c     call marquardt(A,nparm,nparm,pmarquardt*aiimax)

c ?? The call to fullranQ is for testing for singularities in the Hessian.
c ?? The number of those is counted by ising in common singularities.
      if(cholesky) then
        call fullrankQ(A_s,nparm_s)
        call lugaus_new(A_s,B_s,nparm_s,nparm,ajac_s,ndata,ndata,eps,
     &    dseed,newA,pdif_s,prob,pronly,T,pmarquardt*aiimax)
c       call lugaus_new(A,B,nparm,nparm,ajac,ndata,ndata,eps,dseed,
c     &   newA,pdif,prob,pronly,T)
c       write(6,*) 'pdif'
c       write(6,*) (pdif(i),i=1,nparm)
c ?? to make svd_gaus and lugaus_new more analogous the next call to
c ?? uxb should probably be incorporated in lugaus_new
c compute inverse(A).B
        call uxb(A_s,nparm_s,nparm,B_s)
c       call uxb(A,nparm,nparm,B)
      else
c after this call B_s will contain
c inverse(transpose(ajac_s).ajac_s+pmarquardt*aiimax).B
c at T=0 this is minus the drift contribution to pdif
c levenberg uses the result stored in B_s

c rot_wt allows interpolation between shifting the hessian (rot_wt=1) and
c rescaling pdiff and B_s to give a multiple of Newton step + noise
c (rot_wt=0).  rot_wt=0 has to excluded to avoid zero divide.
        if(rot_wt.gt.ONE.or.rot_wt.le.ZERO) stop 'step: invalid rot_wt'
        pmara=ZERO
        do iparm=1,nparm_s
          tmp(iparm)=B_s(iparm)
        enddo
        call svd_gaus_test(A_s,B_s,nparm_s,nparm,ajac_s,ndata,ndata,
     &    ZERO,
     &    dseed,newA,pdif_s,prob,pronly,T,pmara,s_svd,u_svd,v_svd)
        pnorm0_2=ZERO
        do iparm=1,nparm_s
          B_s(iparm)=tmp(iparm)
          pnorm0_2=pnorm0_2+pdif_s(iparm)**2
        enddo
        pmara=pmarquardt*aiimax*rot_wt
        call svd_gaus_test(A_s,B_s,nparm_s,nparm,ajac_s,ndata,ndata,
     &    ZERO,
     &    dseed,.false.,pdif_s,prob,pronly,T,pmara,s_svd,u_svd,v_svd)
        pmarquardt=pmara/(aiimax*rot_wt)
        factor=one/(one+(one-rot_wt)*pmarquardt*aiimax)
        do iparm=1,nparm_s
          pdif_s(iparm)=factor*pdif_s(iparm)
          B_s(iparm)=factor*B_s(iparm)
        enddo
        pnorm1_2=ZERO
        do iparm=1,nparm_s
          pnorm1_2=pnorm1_2+pdif_s(iparm)**2
        enddo
c The next line just avoids zero divides in rare cases. (Peter, Wed May 19 1999)
        if(pnorm1_2 .eq. ZERO) pnorm1_2=ONE
        pratio=sqrt(pnorm0_2/pnorm1_2)
      endif
c      call unscale(B_s,pdif_s,ajac_s,ndata,C_s,nparm_s,prob)
c      call unscale(B,pdif,C,nparm,prob)
      call unselect(pdif_s,B_s,C_s,nparm_s,pdif,B,C,nparm,ipermu)
      if(nparm_s.lt.nparm) call
     &  addgrad(pdif,B,Asav,Bsav,nparm_s,nparm,ipermu,pmarquardt)
      if(pronly) then
        do iparm=1,nparm
          pdif(iparm)=psav(iparm)
        enddo
      endif

c check convergence for parameters
      jparm=0
      istring=0
      do iparm=1,nparm
        if(abs(pdif(iparm)).gt.epsp*max(ONE,abs(p(iparm)))) goto 10
        jparm=iparm
      enddo
10    continue
      if(jparm.lt.nparm) then
        parm_msg=' '
      else
        converg=.true.
        parm_msg='pr '
      endif

      mesg=parm_msg//grad_msg

c check for swiching to central differences
      jparm=0
      epsder1=max(eps_diff,DBL_EPSILON)**THIRD
      do iparm=1,nparm
        if(abs(pdif(iparm)).gt.epsder1*max(ONE,abs(p(iparm)))) goto 20
        jparm=iparm
      enddo
20    continue
      if(jparm.eq.nparm) then
        central=.true.
        if(ipr_com.ge.2) write(6,*) " step: central diffs"
      else
        central=.false.
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine rank(A,Asav,nparm,ipermu,irank)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c select states to form a non-singular sub-matrix from the Hessian A
c this matrix is given by A(ipermu(1:irank),ipermu(1:irank))
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include 'parameters.h'
      dimension A(nparm,nparm),Asav(nparm,nparm),
     &  ipermu(nparm),ipermu_inv(nparm)
c epsnum determines what are called singularities: a conservative choice is
c the sqrt of the relative floating point accuracy
      do i=1,nparm
        do j=1,nparm
          Asav(i,j)=A(i,j)
        enddo
      enddo
c     epsnum=ZERO
      epsnum=1d-14
      call cholesky_piv(Asav,ipermu_inv,ipermu,nparm,nparm,irank,epsnum)
c      print '((''rank = '',i3,'',''20(i3,1x)))',
c     &  irank,(ipermu(i),i=1,nparm)
      return
      end
c-----------------------------------------------------------------------

      subroutine select(A,B,ajac,ndata,C,nparm,pdif,A_s,B_s,ajac_s,
     &  C_s,pdif_s,ipermu,nparm_s)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include 'parameters.h'
      dimension A(nparm,nparm),B(nparm),C(nparm),pdif(nparm),
     &  A_s(nparm,nparm),B_s(nparm),C_s(nparm),pdif_s(nparm),
     &  ipermu(nparm),ajac(ndata,nparm),ajac_s(ndata,nparm)
c ?? to enable pivoting comment out four next for lines
      nparm_s=nparm
      do i=1,nparm
        ipermu(i)=i
      enddo
      do i=1,nparm_s
        B_s(i)=B(ipermu(i))
        C_s(i)=C(ipermu(i))
        pdif_s(i)=pdif(ipermu(i))
        do j=1,nparm_s
          A_s(i,j)=A(ipermu(i),ipermu(j))
        enddo
        do j=1,ndata
          ajac_s(j,i)=ajac(j,ipermu(i))
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine unselect(pdif_s,B_s,C_s,nparm_s,pdif,B,C,nparm,ipermu)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include 'parameters.h'
      dimension pdif_s(nparm_s),B_s(nparm_s),C_s(nparm_s),
     $  pdif(nparm),B(nparm),C(nparm),ipermu(nparm)
      do i=1,nparm_s
        pdif(ipermu(i))=pdif_s(i)
        B(ipermu(i))=B_s(i)
c        C(ipermu(i))=C_s(i)
      enddo
      do i=nparm_s+1,nparm
        pdif(ipermu(i))=ZERO
        B(ipermu(i))=ZERO
c        C(ipermu(i))=ZERO
      enddo
c      write(6,*) 'unselect: pdif:',(pdif(i),i=1,nparm)
      return
      end
c-----------------------------------------------------------------------

      subroutine addgrad(pdif,B,Asav,Bsav,nparm_s,nparm,ipermu,
     &  pmarquardt)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include 'parameters.h'
      dimension pdif(nparm),B(nparm),Asav(nparm,nparm),
     &Bsav(nparm),ipermu(nparm)
      gnorm_r=ZERO
      prgr=ZERO
      do i=1,nparm_s
        gnorm_r=gnorm_r+Bsav(ipermu(i))**2
        prgr=prgr+Bsav(ipermu(i))*B(ipermu(i))
      enddo
      gnorm_s=ZERO
      aiimax_s=ZERO
      do i=nparm_s+1,nparm
        gnorm_s=gnorm_s+Bsav(ipermu(i))**2
        aiimax_s=max(aiimax_s,Asav(ipermu(i),ipermu(i)))
      enddo
      write(6,*) 'addgrad: aiimax_s',aiimax_s,(Asav(i,i),i=1,nparm)
      factor=gnorm_s/(gnorm_r*(ONE+pmarquardt))
      factor=prgr/gnorm_r**2
      do i=nparm_s+1,nparm
        B(ipermu(i))=factor*Bsav(ipermu(i))/Asav(ipermu(i),ipermu(i))
        B(ipermu(i))=factor*Bsav(ipermu(i))
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine marquardt(a,nparm,nparm_phys,pmarquardt,aiimax,C)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c add pmarquardt*aiimax times the unit matrix to the LM hessian
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit real*8(a-h,o-z)
      include 'parameters.h'
      common /temp/ sqrttemp,sqrttemp_inv
      dimension a(nparm_phys,nparm_phys),C(nparm_phys)
      amaxa=a(1,1)
      amina=a(1,1)
      do i=1,nparm
        amaxa=max(amaxa,a(i,i))
        amina=min(amaxa,a(i,i))
      enddo
      amaxa=amaxa/aiimax
      amina=amina/aiimax
c keep parmar within bounds
      pmarquardt=min(pmarquardt,amaxa/DBL_EPSILON)
c?? The following choice may well be too small. We've only used
c?? the next line with amina->amaxa, which should be 1.
      pmarquardt=max(pmarquardt,DBL_EPSILON*amina,sqrt(DBL_MIN))
c     pmarquardt=max(pmarquardt,DBL_EPSILON*amaxa)

      pa=pmarquardt*aiimax
      do i=1,nparm
c use of pure gradient directions with unscaled parameters
         a(i,i)=a(i,i)+pa
c ?? to enable scaling of the marquardt parameter use one of the next lines
c        a(i,i)=a(i,i)+pmarquardt*C(i)**2
c        a(i,i)=(ONE+pmarquardt)*a(i,i)
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine scale(A,B,ajac,ndata,C,nparm,pdif,pronly)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      dimension A(nparm,nparm),B(nparm),C(nparm),pdif(nparm),
     &  ajac(ndata,nparm)
      logical pronly

      do iparm=1,nparm
        B(iparm)=B(iparm)/C(iparm)
        if(pronly) pdif(iparm)=pdif(iparm)*C(iparm)
        do jparm=1,nparm
          A(iparm,jparm)=A(iparm,jparm)/(C(iparm)*C(jparm))
        enddo
        do idata=1,ndata
          ajac(idata,iparm)=ajac(idata,iparm)/C(iparm)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine unscale(B,pdif,ajac,ndata,C,nparm,prob)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      dimension B(nparm),pdif(nparm),C(nparm),ajac(ndata,nparm)
c Note that b at this point is the drift, not the gradient,
c (and therefore covariant rather than contravariant)
c and so needs to be scaled rather than unscaled.

      do iparm=1,nparm
        B(iparm)=B(iparm)/C(iparm)
        pdif(iparm)=pdif(iparm)/C(iparm)
        prob=prob+log(c(iparm))
c ?? is it necessary to unscale the Jacobian ?
        do idata=1,ndata
          ajac(idata,iparm)=ajac(idata,iparm)*C(iparm)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------

      logical function stop_chi2(chi2_new,chi2_old,eps_chi2,max_number,number_cur)
c determines convergence of chi2
c chi2_new          = latest chi2 value
c chi2_old          = smallest chi2 values so far, in ascending order (upto max_number of them)
c eps_chi2          = converged if smallest chi2 and the max_number^th smallest chisq differ by
c                     less than a factor eps_chi2 relative to the smallest chi2
c max_number        = number of chi2's saved and compared
c number_cur        = number of chi2's currently saved
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, September 1996.
c Fixed bug (le in dowhile was lt) and improved stopping criterion, Cyrus, 18 Nov 2011
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      dimension chi2_old(max_number+1)
      i=1
      if(number_cur.gt.0) then
        do while(i.le.number_cur.and.chi2_new.gt.chi2_old(i))
          i=i+1
        enddo
      endif
      if(number_cur.lt.max_number) number_cur=number_cur+1
      do j=number_cur,i+1,-1
        chi2_old(j)=chi2_old(j-1)
      enddo
      chi2_old(i)=chi2_new
c     write(6,'(''stop_chi2:'',(100es12.3))') (chi2_old(ip),ip=1,number_cur)
      if((number_cur.eq.max_number .and. chi2_old(max_number)-chi2_old(1).le.chi2_old(1)*eps_chi2) .or. chi2_old(1).eq.0.d0) then
        stop_chi2=.true.
      else
        stop_chi2=.false.
      endif
      return
      end

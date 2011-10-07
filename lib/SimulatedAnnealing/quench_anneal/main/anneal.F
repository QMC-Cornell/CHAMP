c $Log: anneal.F,v $
c Revision 1.2  2007/07/12 15:39:21  toulouse
c updating Makefile for Argonne bluegene
c J. Toulouse
c
c Revision 1.1.1.1  2006/10/08 19:38:17  toulouse
c champ
c
c Revision 1.17  2005/10/18 14:40:42  nigh
c removed redundant declarations
c
c Revision 1.16  2004/02/21 04:02:35  nigh
c trivial changes
c
c Revision 1.15  2004/02/06 14:11:20  nigh
c *** empty log message ***
c
c Revision 1.14  2003/07/14 17:41:36  nigh
c *** empty log message ***
c
c Revision 1.13  2003/01/02 20:01:43  nigh
c changed commentary for jacobian
c
c Revision 1.12  2002/12/20 19:00:32  nigh
c Cyrus's changes: some exact derivatives
c
c Revision 1.11  2002/03/04 14:21:30  nigh
c new version with LAPACK and LINPACK; see svd_gaus_test
c
c Revision 1.10  2002/01/17 22:27:05  nigh
c *** empty log message ***
c
c Look for "c ??" to find questionable code
c Lines to return to best values found at each temperature can restored by
c setting the following switch to 1
#define ReturnToBest 0
      subroutine anneal(func,jacobian,nanalytic,pold,pmarquardt,tau,
     &  noutput,nstep,ndata,nparm,ipr,diff,chisqmin,epsg,epsp,epsch2,
     &  converg,mesg,dseed,temp,ntemp,cholesky,rot_wt,eps_diff)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c func		= function appearing in sum of squares
c                 function func(ndata,nparm,parm,diff,iflag)
c jacobian	= if nanalytic>0 jacobian evaluates jacobian
c                 subroutine jacobian(ndata,nparm,nanalytic,parm,ajac)
c                 (see subroutine derivs for details)
c pold		= initial parameter estimates
c pmarquardt	= initial value of Marquardt parameter
c tau		= sets time scale of changes of the Marquardt parameter
c noutput	= produce output every this many steps
c nstep		= number of anneal steps performed per temperature
c ndata		= number of data to be fitted
c nparm		= number of parameters
c ipr           = print flag: ipr.lt.0 no output, except fatal errors
c                             ipr.eq.0 almost no output
c                             ipr.eq.1 parameters and object function for
c                                      interation steps determined by noutput
c                             ipr.eq.2 full diagnostic output at those times
c epsg,epsp,epsch2,converg,mesg = see subr. step
c epsg		= convergence if norm of gradient, gnorm < epsg,
c	          where gnorm = 1-norm of scaled gradient with
c		  components B(i)/A(i,i)
c epsp		= convergence if for all i:
c	          pnorm <= epsp abs(1-max(1,pnorm/pnorm'))
c		  where pnorm = 1-norm of the vector with components
c	          pdif(i)/max(1,|p(i)|) and pnorm' the value of this
c                 quantity at the previous iteration
c epsch2        = convergence if max_number of sum-of-square values
c                 have relative differences less than epsch2.
c                 max_number=4 by default. see common/chi2_min/
c temp		= temp(1:ntemp) temperatures at which annealing will done
c cholesky	= .true. use LU decomposition; .false. use singular
c		  value decomposition
c eps_diff      = relative numerical accuracy of residues, used if greater
c                 than dbl_epsilon for increments in numerical
c                 differentiation
c dseed		= probably obsolete thoughout, however ransr must be
c		  initialized with a call to ransi

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include '../include/parameters.h'
      common /temp/ sqrttemp,sqrttemp_inv
      dimension A(nparm,nparm),Asav(nparm,nparm),B(nparm),
     &  Bsav(nparm),pold_sav(nparm),pnew(nparm),
     &  pbest(nparm),pdif(nparm)
      common /singularies/ ising
      common/quenchsim_diff/ndim1_diff,ndim2_diff
      dimension pold(nparm),diff(ndata,ndim2_diff),temp(ntemp)
      logical newA,pronly,converg,central,use_sav,use_sav1,
     $  update_B,detailed_balance,cholesky
      integer accepted
      external func,jacobian
      character*10 mesg
c common /quenchsim_pr/ to communicate amount of output
      logical called_by_qa
      common /quenchsim_pr/ ipr_com,called_by_qa
c common /chi2_min/ stores parameters for stop_chi2
      parameter (MAX_NUMBER_DIM=4)
      common /chi2_min/ chi2_old(MAX_NUMBER_DIM),max_number
      dimension ajac(ndata,nparm)
      dimension s_svd(ndata+1),u_svd(nparm,ndata),v_svd(nparm,nparm)
      ipr_com=ipr
      if(ipr.ge.0)
     &  write
     &    (6,'(/,''Anneal: global least-squares optimizer. Authors:'',
     &    '' Peter Nightingale and Cyrus Umrigar'')')

      if(nparm.lt.0) then
        write(6,'(''Anneal: nparm <= 0'',i8)')nparm
        stop 'Anneal: nparm < 0'
      endif

      if(max_number.gt.MAX_NUMBER_DIM) then
        write(6,'(''max_number.gt.MAX_NUMBER_DIM'',2i20)')
     &    max_number,MAX_NUMBER_DIM
        stop 'Anneal: max_number.gt.MAX_NUMBER_DIM'
      endif

      ising=0
      central=.false.
c     twopi=EIGHT*atan(ONE)

c     read (5,*) ntemp,(temp(i),i=1,max(ntemp,1))
c     write (6,'(''ntemp,temp='',i5,/,(9f8.4))') ntemp,
c    1(temp(i),i=1,max(ntemp,0))
c ??  if(ntemp.le.0) temp(1)=1

c     read *,dseed
c     write (6,'(''dseed'',f20.10,/)') dseed

c     call ransi(int(dseed))

      pmarquardt_sav=pmarquardt
      epsg_sav=epsg
      epsp_sav=epsp
      noutput_sav=noutput
      if(noutput.le.0) noutput=1

      amemory=exp(-ONE/tau)
      av_mov_acc=HALF
c Find the initial value of efo to set initial temperature
      sqrttemp=sqrt(temp(1))
      sqrttemp_inv=ONE/sqrttemp
      called_by_qa=.true.
      efo=func_mpi(func,ndata,nparm,pold,diff,iflag)
      called_by_qa=.false.
      tempsav=temp(1)
      chisqmin=efo

      do iparm=1,nparm
        pbest(iparm)=pold(iparm)
      enddo

c Initialize efo
      Tmin=ONE/1 000 000 000
      newA=.true.
      use_sav=.false.
      do itemp=1,iabs(ntemp)
        chisqcur=efo
c Update the temperature using the cooling scheme
        if(ntemp.le.0) temp(itemp)=Tnext(T0,Tmin,itemp,-ntemp,
     &   chisqcur,chisqmin)
c ??
        if(itemp.eq.1) temp(1)=min(tempsav,temp(1))
        write(6,'(''itemp,temp'',i5,2d12.4)') itemp,temp(1),temp(itemp)
c       write(6,'(''temp,chisqmin'',2d14.6)') temp(itemp),chisqmin
ctemp if(temp(itemp).lt.chisqmin.and.itemp.eq.1) temp(itemp)=chisqmin
        sqrttemp=sqrt(temp(itemp))
        sqrttemp_inv=ONE/sqrttemp
c       called_by_qa=.true.
c       efo=func_mpi(func,ndata,nparm,pold,diff,iflag)
c       called_by_qa=.false.
        tempsav=temp(itemp)
        av=0
        av2=0
        iacceptance=0
        ipassed=0
        do istep=1,nstep
           pmarquardt=parmars(pmarquardt,av_mov_acc,tau)

           pronly=.false.
           update_B=.false.
           call step(A,B,Asav,Bsav,ndata,nparm,nanalytic,pold,pdif,diff,
     &       pmarquardt,dseed,newA,pronly,probo_to_n,func,jacobian,
     &       temp(itemp),epsg,epsp,converg,mesg,central,
     &       use_sav,update_B,cholesky,pratio,rot_wt,eps_diff,aiimax,
     &       ajac,s_svd,u_svd,v_svd)

           if(converg) then
             ipassed=ipassed+1
           else
             ipassed=0
           endif

           drifto=0
           pdifo=0
           dotprodo=0
           changeo=0
c in quench B=-pdif, because B = inverse(Marquardt-shiftes quasi-Hessian) times gradient
c in anneal the latter also holds, but the former does not because of the gaussian noise
c added to the gradient.
           do iparm=1,nparm
             drifto=drifto+B(iparm)**2
             pdifo=pdifo+pdif(iparm)**2
             dotprodo=dotprodo+B(iparm)*pdif(iparm)
             changeo=changeo+Bsav(iparm)*pdif(iparm)
           enddo

           do j=1,nparm
             pnew(j)=pold(j)+pdif(j)
           enddo
c Calculate reverse probs
           do iparm=1,nparm
             pdif(iparm)=-pdif(iparm)
           enddo

           pronly=.true.
           detailed_balance=.false.
           if(detailed_balance) then
c ?? The next code is what has to be done to be exact for detailed balance.
c ?? However, this was observed to sometimes give very low acceptances
c ?? unless a very large Mar. par. was used.
             use_sav1=.false.
             update_B=.false.
c the following call saves A and B in A and B which means that they are not
c saved
             call step(A,B,A,B,ndata,nparm,nanalytic,pnew,pdif,diff,
     &         pmarquardt,dseed,newA,pronly,probn_to_o,func,jacobian,
     &         temp(itemp), epsg,epsp,converg,mesg,central,
     &         use_sav1,update_B,cholesky,pratio,rot_wt,eps_diff,aiimax,
     &         ajac,s_svd,u_svd,v_svd)
           else
c ?? Maybe this problem goes away if the same A and an updated B is used for
c ?? the reverse move.  This is faster too!  Only exact for a quadratic chi^2!
             use_sav1=.true.
             update_B=.true.
             call step(A,B,Asav,Bsav,ndata,nparm,nanalytic,pnew,pdif,
     &         diff,pmarquardt,dseed,newA,pronly,probn_to_o,func,
     &         jacobian,temp(itemp),epsg,epsp,converg,mesg,central,
     &         use_sav1,update_B,cholesky,pratio,rot_wt,eps_diff,aiimax,
     &         ajac,s_svd,u_svd,v_svd)
           endif

           driftn=0
           dotprodn=0
           changen=0
           do iparm=1,nparm
             driftn=driftn+B(iparm)**2
             dotprodn=dotprodn+B(iparm)*pdif(iparm)
             changen=changen+Bsav(iparm)*pdif(iparm)
           enddo

           called_by_qa=.true.
           efn=func_mpi(func,ndata,nparm,pnew,diff,iflag)
           called_by_qa=.false.
           if(efn.lt.chisqmin) then
             temp_chisqmin=temp(itemp)
             chisqmin=efn
             do iparm=1,nparm
               pbest(iparm)=pnew(iparm)
             enddo
           endif
           attemptratio=(probn_to_o-probo_to_n)
c Instead of a probability p the following quantity is T log(p)
           r=(efo-efn)+temp(itemp)*attemptratio
           efo_sav=efo
           do j=1,nparm
             pold_sav(j)=pold(j)
           enddo
           call ransr(ran,1)
c Warning temporary change
c ??       if(r.gt.temp(itemp)*log(ran)) then
c ??       if(efo-efn.gt.2*temp(itemp) .or. r.gt.temp(itemp)*log(ran)) then
c ??       if(efo.gt.efn .or. r.gt.    (itemp)*log(ran)) then
c ??       if(efn-efo.lt.temp(itemp) .or. r.gt.temp(itemp)*log(ran))then
           if(efn-efo.lt.2*temp(itemp) .or. r.gt.temp(itemp)*log(ran))
     &       then
             accepted=1
             use_sav=.false.
             do j=1,nparm
               pold(j)=pnew(j)
             enddo
             efo=efn
             iacceptance=iacceptance+1
             av_mov_acc=amemory*av_mov_acc+(ONE-amemory)
           else
             accepted=-1
             use_sav=.true.
             av_mov_acc=amemory*av_mov_acc
           endif
c If we ever need the next line remember that r used to be the log of
c the Boltzmann distribution; now it's T times that
c          av_mov_acc=amemory*av_mov_acc+(ONE-amemory)*exp(min(r,ZERO))
           av=av+efo
           av2=av2+efo**2

           acceptance=iacceptance
           acceptance=acceptance/istep
           istep_sav=istep

c          write(1,*) efo,(pold(j),j=1,nparm)
c          write(3,*) 2
c          write(3,*) ' '
c          write(3,'(1a,4(f10.5,1x))') 'H',0,0,0
c          write(3,'(1a,4(f10.5,1x))') 'C',(pold(j),j=1,nparm),
c    &       log(efo)
           if(mod(istep-1,noutput).eq.0)
     &     call outputs(istep,nparm,pold_sav,pdif,efo_sav,efn,probn_to_o
     &     ,probo_to_n,
     &     r,pmarquardt,acceptance,accepted,temp(itemp),drifto,pdifo,
     &     dotprodo,driftn,dotprodn,changeo,changen)
           if(ipassed.ge.5) goto 10

        enddo
#if ReturnToBest == 1
c Return to the best point found if
c a) temperature is > 2 times lower than the temperature at which best
c    occurred or at which last return occurred, or
c b) if this is the last temperature.
        if(2*temp(itemp).le.temp_chisqmin .or. itemp.eq.ntemp) then
          temp_chisqmin=temp(itemp)
          efo=chisqmin
          do iparm=1,nparm
            pold(iparm)=pbest(iparm)
          enddo
          use_sav=.false.
        endif
#endif

10      continue
        acceptance=iacceptance
        acceptance=acceptance/istep_sav
        av=av/istep_sav
        av2=av2/istep_sav

        Tprint=max(temp(itemp),REALMIN)
        write(6,'(/,''T,Emin,E,E/T,C,acc.,marq.:'',i3,1x,g8.2,1x,
     &    g12.6,2(1x,g11.5),1x,f9.5,1x,f7.4,1x,g9.2)')
     &    itemp,temp(itemp),chisqmin,av,av/Tprint,
     &    (av2-av**2)/Tprint**2,acceptance,pmarquardt
        write(6,'(''chisqmin '',d12.6,/,(''parm'',6(1x,g12.6)))')
     &    chisqmin,(pbest(ip),ip=1,nparm)
        write(6,*)
      enddo

c Call func so that diffs correspond to best parms ever found
      if(nstep.gt.0) then
        called_by_qa=.true.
        efo=func_mpi(func,ndata,nparm,pold,diff,iflag)
        called_by_qa=.false.
      endif

      if(ising.ne.0) write(6,'(''Anneal: detected '',i6,
     &  '' numerical singularities; see fullrankQ'')') ising

      pmarquardt=pmarquardt_sav
      epsg=epsg_sav
      epsp=epsp_sav
      noutput=noutput_sav

      return
      end
c-----------------------------------------------------------------------

      subroutine quench1(func,jacobian,nanalytic,pold,nstep,
     &ndata,nparm,diff,efo,epsg,epsp,converg,mesg,ibold)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c This subroutine is just to have a simplified calling sequence
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include '../include/parameters.h'
      common/quenchsim_diff/ndim1_diff,ndim2_diff
      dimension pold(nparm),diff(ndata,ndim2_diff)
      external func,jacobian
      character*10 mesg
      logical converg,cholesky

c     pmarquardt=PRECISION
      pmarquardt=1.d-3
      tau=2
      noutput=1
      ipr=0
      cholesky=.true.
      rot_wt=1
      eps_diff=100*DBL_EPSILON
      epschi2=0
      call quench(func,jacobian,nanalytic,pold,pmarquardt,tau,
     &  noutput,nstep,ndata,nparm,ipr,diff,efo,epsg,epsp,epschi2,
     &  converg,mesg,ibold,cholesky,rot_wt,eps_diff)

      return
      end
c-----------------------------------------------------------------------

      subroutine quench(func,jacobian,nanalytic,pold,pmarquardt,tau,
     &  noutput,nstep,ndata,nparm,ipr,diff,efo,epsg,epsp,epsch2,converg,
     &  mesg,ibold,cholesky,rot_wt,eps_diff)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c func		= function appearing in sum of squares
c                 function func(ndata,nparm,parm,diff,iflag)
c jacobian	= if nanalytic>0 jacobian evaluates jacobian
c                 subroutine jacobian(ndata,nparm,nanalytic,parm,ajac)
c                 (see subroutine derivs for details)
c nanalytic     = Number of analytic derivatives provided by jacobian
c pold		= initial parameter estimates
c pmarquardt	= initial value of the Marquardt parameter
c tau		= sets time scale of changes of the Marquardt parameter
c noutput	= produce output every this many steps
c nstep		= maximum number of minimization steps performed
c ndata		= number of data to be fitted
c nparm		= number of parameters
c ipr		= print flag: ipr.lt.0 no output, except fatal errors
c			      ipr.eq.0 almost no output
c			      ipr.eq.1 parameters and object function for
c				       interation steps determined by noutput
c                             ipr.eq.2 full diagnostic output at those times
c diff          = errors at each data point (residuals)
c efo           = sum of squares after minimization
c epsg,epsp,epsch2,converg,mesg = see subr. step
c epsg		= convergence if norm of gradient, gnorm < epsg,
c	          where gnorm = 1-norm of scaled gradient with
c		  components B(i)/A(i,i)
c epsp		= convergence if for all i:
c	          pnorm <= epsp abs(1-max(1,pnorm/pnorm'))
c		  where pnorm = 1-norm of the vector with components
c	          pdif(i)/max(1,|p(i)|) and pnorm' the value of this
c                 quantity at the previous iteration
c epsch2        = convergence if max_number of sum-of-square values
c                 have relative differences less than epsch2.
c                 max_number=4 by default. see common/chi2_min/
c converg 	= .true. if convergence has occurred
c mesg  	= "gr", "pr", "chi" or combinations depending on convergence
c		  criteria met
c ibold 	= boldness level of moves.
c                 0 = no uphill moves accepted
c                 abs(ibold)=1-4 uphill moves accepted with increasing boldness.
c                 Usually a non-zero value of ibold gives faster convergence
c                 but it can occasionally lead to oscillations.
c                 if ibold>0, ibold is adjusted dynamically.
c                 although it seems like a good idea, tests show that dynamic ibold
c                 makes convergence slower, so it is better to use negative ibolds
c cholesky	= .true. use LU decomposition; .false. use singular
c		  value decomposition
c rot_wt 	= interpolate between shifted Gauss-Newton (1) and
c 		  unshifted Gauss-Newton (0).
c eps_diff	= relative numerical accuracy of residues; used if greater
c 		  than dbl_epsilon for increments in numerical differentiation

c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit real*8(a-h,o-z)
      include '../include/parameters.h'
      common /temp/ sqrttemp,sqrttemp_inv
      common/quenchsim_diff/ndim1_diff,ndim2_diff
      dimension A(nparm,nparm),Asav(nparm,nparm),B(nparm),
     &  Bsav(nparm),pold_sav(nparm),pnew(nparm),
     &  pbest(nparm),pdif(nparm),diff_best(ndata),
     &  diff_new(ndata,ndim2_diff)
      dimension pold(nparm),diff(ndata,ndim2_diff),pdif_sav(nparm)
      common /singularies/ ising
c common /quench_info/ used for saving information not returned via
c subroutine parameters
      common /quench_info/ pmarquardt_info,stderr_info,istep_info
      common /ichange/ ichange
c common /quenchsim_pr/ to communicate amount of output
      logical called_by_qa
      common /quenchsim_pr/ ipr_com,called_by_qa
c common /chi2_min/ stores parameters for stop_chi2
      parameter (MAX_NUMBER_DIM=4)
      common /chi2_min/ chi2_old(MAX_NUMBER_DIM),max_number
      logical newA,pronly,central,use_sav,update_B,cholesky
      integer accepted
      external func,jacobian
      character*10 mesg
      character*6 msg
      logical converg,stop_chi2
      dimension ajac(ndata,nparm)
      dimension s_svd(ndata+1),u_svd(nparm,ndata),v_svd(nparm,nparm)

      eps_chi2=epsch2
      ipr_com=ipr
      if(ipr.ge.0)
     &  write
     &  (6,'(/,''Quench: local least-squares optimizer. Authors:'',
     &  '' Peter Nightingale & Cyrus Umrigar'')')

      if(nparm.lt.0) then
        write(6,'(''Quench: nparm <= 0'',i8)')nparm
        stop 'Quench: nparm < 0'
      endif

      if(max_number.gt.MAX_NUMBER_DIM) then
        write(6,'(''max_number.gt.MAX_NUMBER_DIM'',2i20)')
     &    max_number,MAX_NUMBER_DIM
        stop 'Quench: max_number.gt.MAX_NUMBER_DIM'
      endif


      number_cur=0
      ibold_sav=ibold
      pmarquardt_sav=pmarquardt

      epsg_sav=epsg
      epsp_sav=epsp
      noutput_sav=noutput
      if(noutput.le.0) noutput=1
      if(epsg.le.ZERO .and. epsp.le.ZERO) then
        if(ipr.ge.0) write(6,'(''epsg,epsp reset to 1.d-6'')')
        epsg=1.d-6
        epsp=1.d-6
      endif

      ising=0
      central=.false.

c temp is meaningless here but needed for func
      sqrttemp=1
      sqrttemp_inv=sqrttemp
c dseed is meaningless here but needed for rannsu
      dseed=1234

      amemory=exp(-ONE/tau)
c av_mov_acc=ZERO
      iacceptance=0
      ichange=2

c Initialize efo,A,B,chisqmin,pbest,diff_best
      called_by_qa=.true.
      efo=func_mpi(func,ndata,nparm,pold,diff,iflag)
      called_by_qa=.false.
      if(nparm .le. 0 .or. ndata .le. 0) then
c       if(ipr.ge.0) write(6,'(''quench: nparm .le. 0 .or. ndata .le. 0'',2i9)') nparm,ndata
c       return
        write(6,'(''quench: nparm .le. 0 .or. ndata .le. 0'',2i9)') nparm,ndata
        stop 'Does not make sense to optimize with nparm or ndata <=0'
      endif
      chisqmin=efo
      chisqminold=efo
      do iparm=1,nparm
        pbest(iparm)=pold(iparm)
      enddo
      do idata=1,ndata
        diff_best(idata)=diff(idata,1)
      enddo

      newA=.true.
      use_sav=.false.
      av_mov_ratio=1
      accepted=0
      if(pmarquardt.lt.1.d-6) then
        av_mov_acc=0
      else
        av_mov_acc=1
      endif
      pmarquardt_orig=ZERO
      drift_orig=ZERO
      do iparm=1,nparm
        pdif_sav(iparm)=ZERO
      enddo
      factor_sav=ONE
      ichange=1
c make sure the following variables are initialized
      drift_rat=ZERO
      drift_rat_sav=ZERO
c main loop
      istep_best=1
      do istep=1,nstep
        pmarquardt_old=pmarquardt_orig
        drift_old=drift_orig
        if(istep.ge.2) pmarquardt= parmarm(pmarquardt,av_mov_acc,
     &    accepted,cos_on,cholesky)

        pronly=.false.
        update_B=.false.
        call step(A,B,Asav,Bsav,ndata,nparm,nanalytic,pold,pdif,diff,
     &    pmarquardt,dseed,newA,pronly,probo_to_n,func,jacobian,
     &    ZERO,epsg,epsp,converg,msg,central,use_sav,
     &    update_B,cholesky,pratio,rot_wt,eps_diff,aiimax,ajac,
     &    s_svd,u_svd,v_svd)
        mesg=msg

        drift=0
        grad=0
        change=0
        dotprod_on=0
c in quench B=-pdif, because B = inverse(Marquardt-shiftes quasi-Hessian) times gradient
c in anneal the latter also holds, but the former does not because of the gaussian noise
c added to the gradient.
        do iparm=1,nparm
          drift=drift+B(iparm)**2
          grad=grad+Bsav(iparm)**2
          change=change-Bsav(iparm)*B(iparm)
          dotprod_on=dotprod_on-pdif_sav(iparm)*B(iparm)
        enddo
        drift=sqrt(drift)
        if(istep.eq.1) drift_sav=drift
        grad=sqrt(grad)

c save pmarquardt and drift because we may change them if cos_on<0
        pmarquardt_orig=pmarquardt
        drift_orig=drift

c If pmarquardt is having a negligible effect, boost it to have about
c a 10% impact on the drift.  Two alternative formulas are used.  The 2nd
c (which is commented out) is more justifiable, but it make little
c difference anyway.
        if(accepted.eq.-1) then
          pmarquardt_dif=pmarquardt-pmarquardt_old
          drift_dif=drift_old-drift
          if(TEN*abs(drift_dif).lt.drift
c    &      .and.abs(drift_dif).gt.sqrt(PRECISION)*drift) then
     &      .and.drift_dif.ne.ZERO) then
          if(mod(istep-1,noutput).eq.0 .and. ipr.ge.2)
     &      write(6,'(''p_dif,d_dif,mar,drift'',9g12.4)')
     &      pmarquardt_dif,drift_dif,pmarquardt,drift
c make sure that changing the Marquardt parameter changes pdiff
          if(abs(drift_dif).gt.DBL_EPSILON*drift) then
            pmarquardt=
     &        max(pmarquardt+.1*drift*pmarquardt_dif/drift_dif,
     &        pmarquardt)
            else
              pmarquardt=aiimax*DBL_EPSILON
            endif
c           drift_rat=.9d0
c           eig_rat=(drift_old-drift)/
c    &      (drift*pmarquardt-drift_old*pmarquardt_old)
c           pmarquardt=max((ONE/drift_rat-ONE)/eig_rat+
c    &      pmarquardt/drift_rat,pmarquardt)
            call step(A,B,Asav,Bsav,ndata,nparm,nanalytic,pold,pdif,
     &        diff,pmarquardt,dseed,newA,pronly,probo_to_n,func,
     &        jacobian,ZERO,epsg,epsp,converg,msg,central,
     &        use_sav,update_B,cholesky,pratio,rot_wt,eps_diff,aiimax,
     &        ajac,s_svd,u_svd,v_svd)
            mesg=msg

            drift=0
            grad=0
            change=0
            dotprod_on=0
            do iparm=1,nparm
              drift=drift+B(iparm)**2
              grad=grad+Bsav(iparm)**2
              change=change-Bsav(iparm)*B(iparm)
              dotprod_on=dotprod_on-pdif_sav(iparm)*B(iparm)
            enddo
            drift=sqrt(drift)
            grad=sqrt(grad)
            if(mod(istep-1,noutput).eq.0 .and. ipr.ge.2)
     &        write(6,'(''p_new,drift_new,dr_rat'',9g12.4)')
     &        pmarquardt,drift,(drift_old-drift)/drift_old
          endif
        endif

        if(iacceptance.ge.1) then
          cos_on=dotprod_on/(drift_sav*drift)
        else
          cos_on=0
        endif

        factor=ONE
c Deal with non-quadratic zero minima and all non-zero minima.
c This factor can be larger or smaller than 1
        if(drift_sav .gt. 0) then
          drift_rat_raw=cos_on*drift/drift_sav
        else
c The next line just avoids zero divides in rare cases. Division can probably
c avoided altogether. (Peter, Wed May 19 1999)
          drift_rat_raw=cos_on
        endif
c       drift_rat=(drift_rat_raw+factor_sav-ONE)/factor_sav
        drift_rat=drift_rat_raw+ONE-ONE/factor_sav
        if(accepted.ge.0 .and. ONE-abs(cos_on) .lt. .01d0 .and.
     &    abs(drift_rat-drift_rat_sav).lt..1 .and.
     &    abs(drift_rat).lt.ONE) then

c         factor=max(.1d0,min(TEN,ONE/(ONE-drift_rat)))
          factor=max(min(.1d0,.1d0*factor_sav),
     &      min(max(TEN,TEN*factor_sav),ONE/(ONE-drift_rat)))
          if(mod(istep-1,noutput).eq.0 .and. ipr.ge.2)
     &      write(6,'(''non-quad factor'',9d12.4)') factor
     &      ,drift_rat,drift_rat_raw
        endif

c The foll. lines are supposed to take one directly to the quadratic
c minimum along the line.  (They are activated only if the moves are
c approximately on a line.)
c Although this seems like a good idea, it leads to more function evals.
cc      if(accepted.ge.0 .and. ONE-abs(cos_on).lt..01.and.factor.eq.ONE)
c       if(accepted.ge.0 .and. ONE+cos_on .lt. .01 .and. factor.eq.ONE)
c    &    then
c         derlin=change/drift
c         derlin2=(efo_sav-efo-derlin*sign(ONE,-cos_on)*drift_sav)
c    &    *TWO/drift_sav**2
c         if(derlin2.gt.ZERO) factor=-derlin/(derlin2*drift)
c         write(6,'(''derlin,derlin2,factor'',9d12.4)')
c    &    derlin,derlin2,factor
c       endif

c This factor must be <= 1.
c It may be better to leave it unchanged if accepted < 0.
        if(cos_on.lt.0 .and. factor.eq.ONE) then
          if(accepted.ge.0) then
            factor=max(ONE+cos_on,min(ONE,HALF*drift_sav/drift))
c           factor=min(ONE,drift_sav/drift)
c          else
c           factor=HALF
          endif
          if(mod(istep-1,noutput).eq.0 .and. ipr.ge.2)
     &      write(6,'(''cos_on neg factor'',d12.4)') factor
        endif

c Do not allow current drift to be more than 10 times previous one.
c Although this seems like a good idea, it leads to more function evals.
c       if(drift.gt.10*drift_sav .and. iacceptance.gt.0 .and.
c    &    factor.eq.ONE) then
c         factor=10*drift_sav/drift
c         write(6,'(''drift big factor'',d12.4)') factor
c       endif

        factor_sav=factor
        drift_rat_sav=drift_rat

c Change pdif's etc by the factor calculated in one of above blocks.
c Note that this move corresponds only approximately to that
c with the changed pmarquardt.  This could be taken care of by calling
c step again with the new pmarquardt, but then some moves that were nearly
c along a line may no longer be.
        if(factor.ne.ONE) then
          do j=1,nparm
            B(j)=B(j)*factor
          enddo
          drift=drift*factor
          change=change*factor
          dotprod_on=dotprod_on*factor
          pmarquardt=pmarquardt/factor
        endif

        do j=1,nparm
          pnew(j)=pold(j)-B(j)
        enddo

        efo_sav=efo
        do j=1,nparm
          pold_sav(j)=pold(j)
        enddo

        called_by_qa=.true.
        efn=func_mpi(func,ndata,nparm,pnew,diff_new,iflag)
        called_by_qa=.false.

        chisqminold=chisqmin
        if(efn.lt.chisqmin) then
          chisqmin=efn
          do iparm=1,nparm
            pbest(iparm)=pnew(iparm)
          enddo
          do idata=1,ndata
            diff_best(idata)=diff_new(idata,1)
          enddo
          istep_best=istep

          if(ibold_sav.gt.0 .and. ibold.lt.ibold_sav) then
            if(ipr.ge.2) write(6,'(''switch up to ibold='',i2)') ibold+1
            ibold=ibold+1
          endif

        endif

c ?? Warning accepting some uphill moves
        if(efn.lt.efo) then
          if(stop_chi2(efn,chi2_old,eps_chi2,max_number,number_cur))
     &      then
            converg=.true.
            mesg=msg//'chi'
          endif
          do j=1,nparm
            pold(j)=pnew(j)
          enddo
          do idata=1,ndata
            diff(idata,1)=diff_new(idata,1)
          enddo
          av_mov_ratio=amemory*av_mov_ratio+(ONE-amemory)*efn/efo
          efo=efn
          use_sav=.false.
          accepted=max(accepted+1,1)
          iacceptance=iacceptance+1
c         av_mov_acc=amemory*av_mov_acc+(ONE-amemory)
          if(cos_on.gt.ZERO) then
            av_mov_acc=amemory*av_mov_acc+(ONE-amemory)
          else
            av_mov_acc=amemory*av_mov_acc+(ONE-amemory)*HALF
          endif
c        elseif(efn.lt.2*efo.and.av_mov_ratio.lt.ONE) then
        elseif((abs(ibold).ge.1 .and.
     &      iacceptance.eq.0 .and. efn.lt.2*chisqmin) .or.
     &      (abs(ibold).eq.1 .and. efn*(ONE-cos_on).lt.chisqmin) .or.
     &      (abs(ibold).eq.2 .and. efn*(ONE-cos_on)**2.lt.chisqmin) .or.
     &      (abs(ibold).eq.3 .and. efn*(ONE-cos_on).lt.efo) .or.
     &      (abs(ibold).eq.4 .and. efn*(ONE-cos_on)**2.lt.efo)) then
          do j=1,nparm
            pold(j)=pnew(j)
          enddo
          do idata=1,ndata
            diff(idata,1)=diff_new(idata,1)
          enddo
          av_mov_ratio=amemory*av_mov_ratio+(ONE-amemory)*efn/efo
          efo=efn
          use_sav=.false.
          accepted=0
          iacceptance=iacceptance+1
          av_mov_acc=amemory*av_mov_acc+(ONE-amemory)*HALF
        else
          use_sav=.true.
          accepted=min(accepted-1,-1)
          av_mov_acc=amemory*av_mov_acc
        endif

        if(accepted.ge.0) ichange=1
        if(accepted.ge.0.or.istep.eq.1) then
          do j=1,nparm
            pdif_sav(j)=-B(j)
          enddo
          drift_sav=drift
        endif

c       write(2,*) efo,(pold(j),j=1,nparm)
c       write(4,*) 1
c       write(4,*) ' '
c       write(4,'(1a,4(f10.5,1x))') 'H',0,0,0
c       write(4,'(1a,4(f10.5,1x))') 'C',(pold(j),j=1,nparm),log(efo)
        acceptance=iacceptance
        acceptance=acceptance/istep
        istep_sav=istep
        if(mod(istep-1,noutput).eq.0 .and. ipr.gt.0)
     &    call outputm(istep,nparm,pold_sav,B,efo_sav,efn,chisqminold,
     &    pmarquardt,acceptance,accepted,drift,grad,change,cos_on,
     &    av_mov_acc,av_mov_ratio)

c Temporarily change to ibold=0 near end of run
        if(ibold.ne.0.and.istep.gt.max(nstep-10,4*nstep/5)) then
          if(ipr.ge.2) write(6,'(''switch to ibold=0'')')
          ibold=0
        endif

c If no better minimum is found for some time temporarily reduce ibold
        if(ibold.gt.0 .and.
     &  (10*(istep-istep_best).gt.nstep .or.istep-istep_best.gt.10))then
          if(ipr.ge.2) write(6,'(''switch dn to ibold='',i2)') ibold-1
          ibold=ibold-1
        endif

        if(converg) goto 10
      enddo

 10   continue

c Return to best parms found
      efo=chisqmin
      do iparm=1,nparm
        pold(iparm)=pbest(iparm)
      enddo
      do idata=1,ndata
        diff(idata,1)=diff_best(idata)
      enddo

c Call func so that diffs correspond to best parms ever found
      if(nstep.gt.0) then
        called_by_qa=.true.
        efo=func_mpi(func,ndata,nparm,pold,diff,iflag)
        called_by_qa=.false.
      endif

      err=sqrt(efo/ndata)
c     write(6,'(''nstep,acceptance,pmarquardt,chisqmin,err:'',
c    & i6,f9.6,3(1x,g12.6))') istep_sav,acceptance,pmarquardt,efo,err
      if(ipr.ge.0) then
        write(6,'(''istep,chisqmin,err:'',i6,2(1x,g15.9))')
     &    istep_sav,efo,err
        write(6,'(''acceptance,marquardt parm'',f9.6,1x,g10.2)')
     &    acceptance,pmarquardt
        write(6,'((''parm'',6(1x,g12.6)))') (pold(ip),ip=1,nparm)
        write(6,*)
c       if(ising.ne.0) write(6,'(''Quench: detected '',i6,
c    &  '' numerical singularities; see fullrankQ'')') ising
        call systemflush(6)
      endif

c save relevant information in common block quench_info
      pmarquardt_info=pmarquardt
      stderr_info=err
      istep_info=istep_sav
c restore some input parameters
      ibold=ibold_sav
      pmarquardt=pmarquardt_sav
      epsg=epsg_sav
      epsp=epsp_sav
      noutput=noutput_sav
      write(6,'(''Returning from quench'')')
      call systemflush(6)
      return
      end
c-----------------------------------------------------------------------

      function Tnext(T0,Tmin,n,ntot,chisqcur,chisqmin)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c compute next temperature in cooling scheme:
c Tnext=T0*(1-n/ntot)**alfacool+Tmin or
c Tnext=Tmin*(T0/Tmin)**((n-ntot)/ntot)
c depending on the following switch
#define exponential 2
c T0 is updated using:
c chisqcur,chisqmin: current and current minimum values of chi-square
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include '../include/parameters.h'
c ??  factor=20
      factor=1
#if exponential == 1
      if(n.lt.ntot) then
        T0=max(factor*chisqmin,Tmin)
        Tnext=Tmin*(T0/Tmin)**((ntot-n)/dfloat(ntot-1))
      else
        Tnext=ZERO
      endif
      write(6,'(''Tnext: Tnext chisqmin chisqcur '',3(e10.3,1x))')
     &  Tnext,chisqmin,chisqcur
#elif exponential == 2
      if(n.lt.ntot) then
        T0=factor*chisqmin
        Tnext=T0/(ONE+exp(TEN*(n-.8*ntot)/dfloat(ntot)))
      else
        Tnext=ZERO
      endif
      write(6,'(''Tnext: Tnext chisqmin chisqcur '',3(e10.3,1x))')
     &  Tnext,chisqmin,chisqcur
#else
      T0=max(factor*chisqmin,Tmin)
      alfacool=4
      Tnext=T0*(1-n/dfloat(ntot))**alfacool+Tmin
      write(6,*) '(''Tnext: alfacool Tnext chisqmin chisqcur '',
     &  4(e10.3,1x))', alfacool,Tnext,chisqmin,chisqcur
#endif
#undef exponential
      return
      end
c-----------------------------------------------------------------------

      function parmars(parmar0,a,tau)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c update Marquardt parameter for simulated annealing
c parmars	= new value of M. parameter
c parmar0	= old value of M. parameter
c tau		= relaxation time: the M. parameter may grow or decrease
c		  by a factor CONSTMAR over this time
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include '../include/parameters.h'
c     parmar=a**EXPMAR-CONSTMAR
      parmars=parmar0*CONSTMAR**((ONE-TWO*a)/tau)
c     write(6,'(''parmar,a'',d12.5,f6.3)') parmar,a
      return
      end
c-----------------------------------------------------------------------

      function parmarm(parmar0,a,accepted,cos_on,cholesky)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c update Marquardt parameter for quench
c the factor by which the parameter is changed is itself adjusted
c such that it is close to 1 when the av_mov_acc is close to .5 and
c large when av-mov_acc is 0 or 1.
c parmarm       = new value of M. parameter
c parmar0       = old value of M. parameter
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      include '../include/parameters.h'
      integer accepted
      logical cholesky

      factor=min(100.d0,
     &  max(1.1d0,ONE/(DBL_EPSILON+(ONE-abs(2*a-ONE)))**2))

      if(accepted.ge.1) then
        if(cos_on.ge.ZERO) then
          parmarm=parmar0/factor
         else
          parmarm=parmar0
        endif
       elseif(accepted.le.-1) then
        if(cos_on.gt.ZERO) then
          parmarm=parmar0*sqrt(factor)
         else
          parmarm=parmar0*factor
        endif
       else
        parmarm=parmar0*sqrt(factor)
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine outputs(istep,nparm,pold,pdifm,efo,efn,probn_to_o,
     &  probo_to_n,
     &  r,pmarquardt,acceptance,accepted,temp,drifto,pdifo,
     &  dotprodo,driftn,dotprodn,changeo,changen)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      integer accepted
      dimension pold(*),pdifm(*)
c common /quenchsim_pr/ to communicate amount of output
      logical called_by_qa
      common /quenchsim_pr/ ipr_com,called_by_qa

      write(6,'(''temp,istep='',d12.4,i9)') temp,istep
      write(6,'((''pold'',6(1x,g12.6)))') (pold(i),i=1,nparm)
      write(6,'((''pdif'',6(1x,g12.6)))') (-pdifm(i),i=1,nparm)
      write(6,'(''efo,efn,n_to_o,o_to_n,r'',9d9.2)') efo,
     &   efn,probn_to_o,probo_to_n,r
      write(6,'(''acceptance mar.par. accptd '',f10.4,1x,d9.2,i3)')
     & acceptance,pmarquardt,accepted
      if(ipr_com.ge.2) then
        cosold=-dotprodo/sqrt(drifto*pdifo)
        cosnew= dotprodn/sqrt(driftn*pdifo)
        write(6,'(''drifto,pdifo,driftn'',9d12.4)')
     &    drifto,pdifo,driftn
        write(6,'(''cosold,cosnew'',9d12.4)') cosold,cosnew
        write(6,'(''changeo,changen'',9d12.4)') changeo,changen
      endif
      write(6,*)
      call systemflush(6) !JT
      return
      end
c-----------------------------------------------------------------------

      subroutine outputm(istep,nparm,pold,pdifm,efo,efn,chisqminold,
     &  pmarquardt,acceptance,accepted,drift,grad,change,cos_on,
     &  av_mov_acc,av_mov_ratio)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      integer accepted
      dimension pold(*),pdifm(*)
c common /quenchsim_pr/ to communicate amount of output
      logical called_by_qa
      common /quenchsim_pr/ ipr_com,called_by_qa

      write(6,'(''istep='',i9)') istep
      if(accepted.ge.0)
     &write(6,'((''pold'',6(1x,g12.6)))') (pold(i),i=1,nparm)
      write(6,'((''pdif'',6(1x,g12.6)))') (-pdifm(i),i=1,nparm)
      if(efn.lt.chisqminold) then
        write(6,'((''pnew'',6(1x,g12.6)))') (pold(i)-pdifm(i),i=1,nparm)
        write(6,'(''efo,efn,efo-efn'',g11.4,9g9.2)') efo,efn,efo-efn
     &  ,chisqminold-efn
       else
        write(6,'(''efo,efn,efo-efn'',g11.4,9g9.2)') efo,efn,efo-efn
      endif
      if(ipr_com.ge.0) then
        write(6,'(''mar.par. acceptance av_mov_acc accptd '',
     &  d9.2,2f8.4,i3)') pmarquardt,acceptance,av_mov_acc,accepted
        if(drift*grad.ne.0) then
          cosold=-change/(drift*grad)
         else
          cosold=10d99
        endif
        write(6,'(''drift,grad'',9d12.4)') drift,grad
        write(6,'(''change'',9d12.4)') change
        write(6,'(''cosold cos_on'',9d12.4)') cosold,cos_on
      endif
      write(6,*)
      call systemflush(6) !JT
      return
      end
c-----------------------------------------------------------------------

      subroutine fullrankQ(A,nparm)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Checks to see if matrix A is numerically of full rank.
c If not a message is generated but no other action is undertaken.
c Early version of non-linear least-squares optimizer
c Copyright: M. Peter Nightingale and Cyrus J. Umrigar, May 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit real*8(a-h,o-z)
      dimension a(nparm,nparm),testa(nparm,nparm),ipiv(nparm),
     &  iperm(nparm)
      common /singularies/ ising
      do i=1,nparm
        do j=1,nparm
          testa(i,j)=a(i,j)
        enddo
      enddo
      eps=1e-8
      call cholesky_piv(testa,ipiv,iperm,nparm,nparm,ir,eps)
      if(ir.ne.nparm) then
c ??    write(6,*) "fullrankQ: singular Hessian: rank deficiency = ",
c ?? &  nparm-ir
        ising=ising+1
c ??    write(6,*) (a(i,i),i=1,nparm)
c ??    write(6,*) (iperm(i),i=1,nparm)
c ??    if(ising.eq.1) then
c ??      do i=1,nparm
c ??        write(1,'(10f25.16,1x)')(a(i,j),j=1,nparm)
c ??      enddo
c ??    endif
      endif
      return
      end
c-----------------------------------------------------------------------

      block data chi2_min_block
      implicit real*8(a-h,o-z)
      parameter (MAX_NUMBER_DIM=4)
      common /chi2_min/ chi2_old(MAX_NUMBER_DIM),max_number
      data max_number/4/
      end

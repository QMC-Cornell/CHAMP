      subroutine fit
! Written by Cyrus Umrigar
! Warning: Moved some of the input to read_input.f and
!          Those variables get passed back to here via use statements.
! Uses the variance minimization algorithm described in:
! 1) Optimized Trial Wave Functions for Quantum Monte Carlo Calculations,
!    C.J. Umrigar, K.G. Wilson and J.W. Wilkins, Phys. Rev. Lett.,60, 1719 (1988).
! 2) A Method for Determining Many-Body Wavefunctions, C.J. Umrigar, K.G. Wilson and J.W. Wilkins,
!    {\it Computer Simulation Studies in Condensed Matter Physics: Recent Developments,}
!    ed. by D.P. Landau, K.K. Mon and H.B. Schuttler (Springer-Verlag 1988).
! 3) Two Aspects of Quantum Monte Carlo: Determination of Accurate Wavefunctions and
!    Determination of Potential Energy Surfaces of Molecules, C.J. Umrigar,
!    Int. J. Quant. Chem. Symp., 23, 217 (1989).
      use all_tools_mod
      use mpi_mod
      use constants_mod
      use control_mod
      use allocations_mod
      use atom_mod
      use orbitals_mod
      use dorb_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use basis1_mod
      use basis2_mod
      use basisnorm_mod
      use contr2_mod
      use pseudo_mod
      use contrl_per_mod
      use contrl_opt_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use jaspar6_mod
      use bparm_mod
      use pointer_mod
      use pars_mod
      use jaspar1_mod
      use jaspar2_mod
      use ncusp_mod
      use confg_mod
      use contrl_opt_mod
      use mpioffset_mod
      implicit real*8(a-h,o-z)
      character*80 fmt
      character*10 mesg
      logical converg
      external func,jacobian

      common /focsav/ a4sav,a5sav,a6sav,a7sav

      common /fcn_calls/icalls
      common /update/ ichange

      real(dp), allocatable :: diff(:), parm(:)

! Inputs are:
! 1) ndata,nparm,ijas,icusp,icusp2,isc,nsig,ncalls,iopt,ipr_opt
!    ndata  = no. of configs to optimize over
!    nparm  = no. of parameters to optimize
!    ijas   = form of wavefunction
!    icusp  = <= -1  zeroth order e-N cusp not imposed
!             >=  0  zeroth order e-N cusp imposed "exactly"
!    icusp2   is for imposing e-N and e-e cusps
!             >=  1 impose cusps exactly
!             <= -2 impose/check cusps via penalty (used mostly for ijas=2 for which exact imposition is difficult)
!                The penalty wt. cuspwt depends on whether icusp2 is divisible by 3 or not.
!                e.g. icusp2=-101 gives a wt of 100, icusp2=-102 gives a wt of .01
!                     icusp2=-1001 gives a wt of 1000, icusp2=-1002 gives a wt of .001
!    isc      is used for 3 purposes
!             a) for using spq rather than srt variables
!             b) for calling numerical (jastrow_num,psi) rather than analytical (jastrow)
!                gradient and Laplacian. Analytic if >= -5
!             c) for scaling ri,rj,rij in Pade.
!                -2,-7  (1-exp(scalek(1)*r)/scalek(1)
!                -3,-8  [1-exp{-scalek(1)*r-(scalek(1)*r)**2/2}]/scalek(1)
!                -4,-9  r/(1+scalek(1)*r)
!                -5,-10 r/{1+(scalek(1)*r)**2}**.5
!    nsig     = no. of significant digits required in parameters
!    ncalls   = Max. no. of function calls to make
!    iopt     In fit:
!               Choose between zxssq and quench
!             = 0,1 call zxssq from IMSL (obsolete)
!             = 2   call quench written by Peter Nightingale and Cyrus Umrigar
!             = 0 Do not check for strict downward descent in zxssq
!             = 1 Strict downward descent in zxssq
!               zxssq is obsolete so, if in fit mode, we reset iopt to 2 in read_input
!             In vmc:
!               Choose between Newton, linear and perturbation theory
!             = 1 linear method
!             = 2 modified Newton method
!             = 3 perturbation theory
!    ipr_opt  <= -2  Minimal print out
!             >= -1  Print cusp monotonicity conditions
!             >=  0  Print configs and errors on 6
!             >=  2  Print out configs and wavefunction on 2
! 2) i3body,irewgt,iaver,istrch
! 3) ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
! nparml   no. of linear    parameters
! nparme   no. of exponent  parameters
! nparmd   no. of determin  parameters (obsolete)
! nparmcsf no. of CSF       parameters
! nparms   no. of scale     parameters for good Jastrow (scalek(1))
! nparmg   no. of GAM       parameters (a21)
! nparmj   no. of Jastrow   parameters
! nparma   no. of Jastrow a parameters in complicated Jastrow
! nparmb   no. of Jastrow b parameters in complicated Jastrow
! nparmc   no. of Jastrow c parameters in complicated Jastrow
! nparmf   no. of Jastrow fck parameters in complicated Jastrow

      call my_second(0,'fit   ')

      call common_allocations

      call alloc ('psid', psid, ndata)
      call alloc ('psij', psij, ndata)
      call alloc ('psio', psio, ndata)
      call alloc ('eold', eold, ndata)
      call alloc ('uwdiff', uwdiff, ndata)
      call alloc ('wght', wght, ndata)
      call alloc ('dvpdv', dvpdv, ndata)

      allocate(ircounts(0:nproc))
      allocate(idispls(0:nproc))

! Calculate distances of atoms from center for use in cusorb
!JT comment out because dist_cent is not used
!JT      call alloc ('dist_cent', dist_cent, ncent)
!JT      do 2 i=1,ncent
!JT        dist_cent(i)=0
!JT        do 1 k=1,ndim
!JT    1     dist_cent(i)=dist_cent(i)+cent(k,i)**2
!JT    2   dist_cent(i)=sqrt(dist_cent(i))

! Save values of a(4) thru a(7), in order to adiabatically go from
! 4th order Pade without Fock terms to one with
      if(ifock.gt.0.and.ijas.eq.2) then
        a4sav=a1(4,1,1)
        a5sav=a1(5,1,1)
        a6sav=a1(6,1,1)
        a7sav=a1(7,1,1)
      endif
      if(ifock.eq.4.and.ijas.eq.3) then
        do it=1,nctype
          a4sav=c(4,it,1)
          a5sav=c(5,it,1)
          a6sav=c(7,it,1)
          a7sav=c(9,it,1)
          call scale20(1,it)
        enddo
      endif

! Optimization section
!     read(5,*) section
!     write(6,'(a)') section

!     read(5,*)
!     read(5,*) ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt
!     read(5,*) i3body,irewgt,iaver,istrch
!     read(5,*) ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,
!    &idbdu,idbdt

!     if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt
!    &.gt.0.and.(ijas.ne.2))
!    &stop 'ipos+...>0 checkjas2 exists only be used with Jastrow2'
!     if(mod(irewgt,100).eq.1) then
!       write(6,*) '**Warning irewgt=1 reset to irewgt=10'
!       irewgt=irewgt+9
!     endif
!     do 4 i=1,ndata
!  4    wght(i)=one

!     write(6,'(/,''ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt=''
!    &,20i4)') ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt
!     write(6,'(''i3body,irewgt,iaver,istrch'',9i5)')
!    &i3body,irewgt,iaver,istrch
!     write(6,'(''ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt'',10i8)')
!    & ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdt,idbdu

!     if(isc.eq.2) write(6,'(
!    &''dist scaled r=(1-exp(-scalek*r))/scalek'')')
!     if(isc.eq.3) write(6,'(
!    &''dist scaled r=(1-exp(-scalek*r-(scalek*r)**2/2))/scalek'')')
!     if(isc.eq.4) write(6,'(
!    &''dist scaled r=r/(1+scalek*r)'')')
!     if(isc.eq.5) write(6,'(
!    &''dist scaled r=r/(1+(scalek*r)**2)**.5'')')

!     if(ijas.eq.1) write(6,'(''Conventional Jastrow'')')
!     if(ijas.eq.2) write(6,'(''Exp. Pade + non-anal terms'')')
!     if(ijas.eq.3) write(6,'(''Standard form'')')
!     if(ijas.eq.4) write(6,'(''New transferable standard form 4'')')
!     if(ijas.eq.5) write(6,'(''New transferable standard form 5'')')
!     if(ijas.eq.6) write(6,'(''New transferable standard form 6'')')

!     if(icusp.ge.0) write(6,'(''Nuclear cusp constraint is imposed'')')

!     read(5,*) (lo(iorb),iorb=1,norb)
!     write(6,'(''lo='',20i3)') (lo(iorb),iorb=1,norb)
!c    read(5,*) (n(ib),l(ib),ib=1,nbasis)
!c    write(6,'(''n,l='',20(2i3,1x))') (n(ib),l(ib),ib=1,nbasis)

!     if(ijas.le.3) then
!       na1=nspin1
!       na2=nspin2
!      else
!       na1=1
!       na2=nctype
!     endif
!     read(5,*) nparml,(nparma(ia),ia=na1,na2),
!    &(nparmb(isp),isp=nspin1,nspin2b),(nparmc(it),it=1,nctype),
!    &(nparmf(it),it=1,nctype),nparmd,nparms,nparmg

!     if(ijas.ge.4.and.ijas.le.6) then
!       do 5 it=1,nctype
!         if(numr.le.0) then
!cAll-electron with analytic slater basis
!           if((norda.eq.0.and.nparma(it).gt.0)
!    &      .or.(norda.gt.0 .and. nparma(it).gt.norda+1)) then
!             write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
!             stop 'nparma too large for norda in all-electron calculation'
!           endif
!          else
!cPseudopotential with numerical basis (cannot vary a(1) or a(2)
!           if(norda.eq.1) stop 'makes no sense to have norda=1 for numr>0'
!           if((norda.eq.0.and.nparma(it).gt.0)
!    &      .or.(norda.gt.0 .and. nparma(it).gt.norda-1)) then
!             write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
!             stop 'nparma too large for norda in pseudopot calculation'
!           endif
!         endif
!         if(isc.le.7 .and.
!    &       ((nordc.le.2.and.nparmc(it).gt.0)
!    &    .or.(nordc.eq.3.and.nparmc(it).gt.2)
!    &    .or.(nordc.eq.4.and.nparmc(it).gt.7)
!    &    .or.(nordc.eq.5.and.nparmc(it).gt.15)
!    &    .or.(nordc.eq.6.and.nparmc(it).gt.27)
!    &    .or.(nordc.eq.7.and.nparmc(it).gt.43))) then
!           write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
!           stop 'nparmc too large for nordc in J_een with cusp conds'
!         endif
!         if(isc.gt.7 .and.
!    &       ((nordc.le.1.and.nparmc(it).gt.0)
!    &    .or.(nordc.eq.2.and.nparmc(it).gt.2)
!    &    .or.(nordc.eq.3.and.nparmc(it).gt.6)
!    &    .or.(nordc.eq.4.and.nparmc(it).gt.13)
!    &    .or.(nordc.eq.5.and.nparmc(it).gt.23)
!    &    .or.(nordc.eq.6.and.nparmc(it).gt.37)
!    &    .or.(nordc.eq.7.and.nparmc(it).gt.55))) then
!           write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
!           stop 'nparmc too large for nordc without cusp conds'
!         endif
!   5   continue
!cFor the b coefs. we assume that b(1) is fixed by the cusp-cond.
!       do 6 isp=1,nspin1,nspin2b
!           if((nordb.eq.0.and.nparmb(isp).gt.0).or.(nordb.gt.0 .and. nparmb(isp).gt.nordb)) then
!             write(6,'(''isp,nordb,nparmb(isp)'',3i5)') isp,nordb,nparmb(isp)
!             stop 'nparmb too large for nordb'
!           endif
!   6   continue
!     endif

!ccompute nparmj and nparme
!     nparmj=0
!     npointa(1)=nparmj
!     do 7 ia=na1,na2
!       if(ia.gt.1) npointa(ia)=npointa(ia-1)+nparma(ia-1)
!   7   nparmj=nparmj+nparma(ia)
!     do 8 isp=nspin1,nspin2b
!   8   nparmj=nparmj+nparmb(isp)
!     npoint(1)=nparmj
!     do 9 it=1,nctype
!       if(it.gt.1) npoint(it)=npoint(it-1)+nparmc(it-1)
!   9   nparmj=nparmj+nparmc(it)+nparmf(it)
!     nparme=nparm-nparml-nparmj-nparmd-nparms-nparmg
!     write(6,'(''No of linear coefs, exponents, Jastrow, det, scale parms varied='',9i5)')
!    &nparml, nparme, nparmj, nparmd, nparms
!     if(nparme.lt.0) stop 'nparme < 0'
!     if(nparme.gt.nbasis) stop 'nparme > nbasis'
!     if(nparme.gt.0 .and. numr.gt.0) stop 'nparme > 0 and numr > 0'
!     if(nparml.lt.0 .or. nparmj.lt.0 .or. nparmd.lt.0 .or. nparms.lt.0 .or.nparmg.lt.0)
!    &stop 'nparm? must be >= 0'
!     if(nparms.gt.1) stop 'nparms must be 0 or 1'

!     read(5,*) (iworb(iparm),iwbasi(iparm),iparm=1,nparml)
!     write(6,'(''lin. coefs. of orbs varied='',10(2i3,2x))')
!    &(iworb(iparm),iwbasi(iparm),iparm=1,nparml)

!     read(5,*) (iwbase(iparm),iparm=1,nparme)
!     write(6,'(''exponents varied='',20i3)') (iwbase(iparm),iparm=1,
!    &nparme)

!     read(5,*) (iwdet(iparm),iparm=1,nparmd)
!     write(6,'(''determinantal coefs varied='',20i3)')
!    &(iwdet(iparm),iparm=1,nparmd)

!     if(ijas.eq.2.or.ijas.eq.3) then
!       write(6,'(''Correl. params. that are varied are:'')')
!       do 14 isp=nspin1,nspin2
!         read(5,*) (iwjasa(iparm,isp),iparm=1,nparma(isp))
!  14     write(6,'(''a: '',30i3)') (iwjasa(iparm,isp),iparm=1,
!    &    nparma(isp))
!       do 16 isp=nspin1,nspin2b
!         read(5,*) (iwjasb(iparm,isp),iparm=1,nparmb(isp))
!  16     write(6,'(''b: '',30i3)') (iwjasb(iparm,isp),iparm=1,
!    &    nparmb(isp))
!      elseif(ijas.ge.4.and.ijas.le.6) then
!       do 18 it=1,nctype
!         read(5,*) (iwjasa(iparm,it),iparm=1,nparma(it))
!  18     write(6,'(''a: '',30i3)') (iwjasa(iparm,it),iparm=1,
!    &    nparma(it))
!       do 20 isp=nspin1,nspin2b
!         read(5,*) (iwjasb(iparm,isp),iparm=1,nparmb(isp))
!  20     write(6,'(''b: '',30i3)') (iwjasb(iparm,isp),iparm=1,
!    &    nparmb(isp))
!     endif
!     if(ijas.ge.3.and.ijas.le.6) then
!       do 25 it=1,nctype
!         read(5,*) (iwjasc(iparm,it),iparm=1,nparmc(it))
!  25     write(6,'(''c: '',60i3)') (iwjasc(iparm,it),iparm=1,
!    &    nparmc(it))
!       if(ifock.gt.0) then
!         do 26 it=1,nctype
!           read(5,*) (iwjasf(iparm,it),iparm=1,nparmf(it))
!  26       write(6,'(''f: '',30i3)') (iwjasf(iparm,it),iparm=1,
!    &      nparmf(it))
!       endif
!     endif

!! set quadrature points
!      if(nloc.gt.0) call rotqua
!
!      write(6,'(''nparm,nparml,nparmj,nparmcsf,nparms,nparmg,nparme='',9i5)') nparm,nparml,nparmj,nparmcsf,nparms,nparmg,nparme
!      write(6,'(''total number of parameters, nparm='',i5)') nparm
!
!      if(mbasis_ctype.gt.0) then   ! Needed since nbasis_ctype(it) is not initialized for quantum rings
!        call alloc ('imnbas', imnbas, ncent)
!        imnbas(1)=1
!        do i=1,ncent-1
!          it=iwctype(i)
!          imnbas(i+1)=imnbas(i)+nbasis_ctype(it)
!        enddo
!      endif
!
!      read(5,*) necn,nebase
!
!! If necn>=0, figure out which orbital (LCAO) coefs are varied and which are constrained to be equal
!      if(necn.ge.0) then
!        write(6,'(''No of linear coefs, exponents set equal='',3i5)') necn,nebase
!
!        call alloc ('ieorb', ieorb, 2, necn)
!        call alloc ('iebasi', iebasi, 2, necn)
!        read(5,*) ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)
!
!      else
!
!        call orb_params
!
!      endif
!
!      do j=1,necn
!        do i=1,2
!          if(ieorb(i,j).gt.norb) stop 'ieorb(i,j).gt.norb'
!          if(iebasi(i,j).gt.nbasis) stop 'iebasi(i,j).gt.nbasis'
!        enddo
!      enddo
!
!      write(6,'(''lin. coefs of orbs (ieorb,iebasi) set equal='',1(2(2i3,2x),2x))') ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)
!      write(6,'(''Here1'')')
!      write(6,'(''lin. coefs of orbs (ieorb,iebasi) set equal='',5(2(2i3,2x),2x))') ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)
!
!! If nebase>=0, figure out which basis exponents varied and which are constrained to be equal
!      if(nebase.ge.0) then
!        write(6,'(''No of linear coefs, exponents set equal='',3i5)') necn,nebase
!
!        call alloc ('iebase', iebase, 2, nbasis)
!        read(5,*) ((iebase(i,j),i=1,2),j=1,nebase)
!        write(6,'(''expon. (iebase) set equal='',10(2i3,2x))') ((iebase(i,j),i=1,2),j=1,nebase)
!
!        do j=1,nebase
!          do i=1,2
!            if(iebase(i,j).gt.nbasis) stop 'iebase(i,j).gt.nbasis'
!          enddo
!        enddo
!
!      else
!
!        call basis_exp_params
!
!      endif
!
!
!      call systemflush(6)
!
!      read(5,*) (ipivot(j),j=1,norb)
!      write(6,'(''ipivot='',10i4)') (ipivot(j),j=1,norb)
!
!      read(5,*) eguess
!      write(6,'(''eguess='',f12.6)') eguess
!
!      read(5,*) pmarquardt,tau_fit,noutput,nstep_fit,ibold
!      write(6,'(''pmarquardt,tau_fit,ibold '',2g9.2,i3)') pmarquardt,tau_fit,ibold
!
!      read(5,'(2l2)') analytic_jacobian,cholesky
!      write(6,'(''analytic_jacobian,cholesky'',2l2)') analytic_jacobian,cholesky
!
!      if(analytic_jacobian .and. irewgt.ne.0) write(6,'(''*** Warning *** analytic derivs not correctly implemented when points are reweighted'')')
!      if(analytic_jacobian .and. nparmf(1).gt.0) stop 'analytic derivatives not implemented yet for Fock parameters'

! If necn<0, call orb_params figure out which orbital (LCAO) coefs are varied and which are constrained to be equal
      if(necn < 0 .or. nebase < 0) then
        call orb_params
        write(6,'(/,''After orb_params: nparm, nparml, nparme, nparmcsf, nparmj, nparms='',6i5,/)') nparm, nparml, nparme, nparmcsf, nparmj, nparms
!       write(6,'(''lin. coefs of orbs (ieorb,iebasi) set equal='',/,5(2(2i3,2x),2x))') ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)
!       write(6,'(''expon. (iebase) set equal='',/,10(2i3,2x))') ((iebase(i,j),i=1,2),j=1,max(0,nebase))
      endif

      if(ndata.lt.nparm) then
        write(6,'(''ndata, nparm='',2i6)') ndata, nparm ; call systemflush(6)
        stop 'ndata < nparm'
      endif

! If we are doing analytic derivatives wrt the optimization parameters
! then we should not optimize the 1st a parameter because then coef
! changes to maintain the cusp, but the change of the residues wrt
! the change in coef is not calculated analytically.
! If there is only one atom type, then we could easily modify the program
! to include the 1st a parameter in the numerical rather than the analytic
! set, but it is hard to do this when there is more than one atom type
! because in quench it is assumed that the numerical parameters appear before
! the analytic ones.
      if(analytic_jacobian .and. icusp.eq.1) then
        do it=1,nctype
          do iparm=1,nparma(it)
            if(iwjasa(iparm,it).eq.1) stop 'Does not make sense to optimize linear e-n parameter if one is using analytic jacobian and there is more than one atom type.  For one atom type program needs modification to do this.'
          enddo
        enddo
      endif

! norbc = number of e-n cusp conditions.  At each nucleus there is one for every nonzero occupied orbital.
! But some orbitals can have an s component at some nuclei and not others.
! Note at present norb has been reset to norb_occ, so only occupied orbitals are considered here.
! Warning: We are presently incorrectly assuming that an orbital is nonzero at all or none of the nuclei, so number of penalties is ncent*norbc.
! We are calling cuspco from func not jacobian, but func is called whether analytic_jacobian is true or not
! Now, lo from input is no longer used and can be removed.
!     write(6,'(''lo='',1000i2)') lo
      norbc=0
!     do iorb=1,norb
!       if(lo(iorb).eq.0 .and. numr.le.0) norbc=norbc+1
!     enddo
!     write(6,'(''norbc='',i5)') norbc

!     write(6,'(/,''norb='',i5)') norb
      if(icusp.ge.0) then
        do iorb=1,norb
          do icent=1,ncent
            ibas=imnbas(icent)
            if(abs(coef(ibas,iorb,1)).ge.1d-9) then
              norbc=norbc+1
              lo(iorb)=0
              exit
            endif
          enddo
        enddo
      endif
      write(6,'(/,''norb='',i5,'', Number of orbs nonzero at any nucleus, norbc='',i5)') norb, norbc
      write(6,'(''lo='',1000i2)') lo

! The penalties can come from
! a) imposing en cusp on s orbitals,
! b) imposing cusp conds. on analyt. een terms in Jastrow (ijas=2,3)
! c) imposing cusp conds. on Fock een terms in Jastrow (ijas=3)
! d) imposing shape conditions on Jastrow (ijas=2)
! ncent*norbc = number of en cusp conds. on orbitals that are nonzero at nuclei (presently assume if nonzero at one nucleus, it is nonzero at all, though this could easily be changed)
! ncuspc = number of cusp conditions from Jastrow (actual # of conds. is ncuspc*(nspin2-nspin1+1)) for ijas=2,3
! nfockc = number of cusp conditions from Fock terms for ijas=3
! ncnstr = number of shape conditions, (actual # of conds. is ncnstr*(nspin2-nspin1+1)) for ijas=2
      nfock=0
      nfockc=0
      ncuspc=0
      cuspwt=one
      if(iabs(icusp2).ge.2) then
        if(ijas.eq.2) then
          ncuspc=38
        elseif(ijas.eq.3) then
          ncuspc=2*nord*nctype
          if(ifock.eq.2) nfock=2
          nfockc=nfock*nctype
          write(6,'(''nfock,nfockc'',3i5)') nfock,nfockc
        endif
        ndata2=ndata+ncent*norbc+ncuspc*(nspin2-nspin1+1)+nfockc
        if(mod(icusp2,3).ne.0) then
          cuspwt=    dfloat(iabs(icusp2)-1)
        else
          cuspwt=one/dfloat(iabs(icusp2)-2)
        endif
      else
        ndata2=ndata+ncent*norbc
      endif

! The call to checkjas2 here is only to calculate ncnstr; the diffs are not calculated and so do not need correct dimension
! After the call to checkjas2, diff is reallocated correctly.
      ncnstr=0
      ishft=ndata2+1
      call alloc ('diff', diff, ishft) !JT
      if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt.gt.0 .and. ijas.eq.2) then
        icalcul_diff=0
        call checkjas2(scalek(1),1,ncnstr,diff(ishft),ipr_opt,0,icalcul_diff)
      endif
      ndata2=ndata2+ncnstr*(nspin2-nspin1+1)
      call alloc ('diff', diff, ndata2)

      write(6,'(''ndata, ndata2, # of s-orb. cusp cond., # of analyt. cusp cond., # of Fock cond., # of shape cond. ='',6i5)') &
     &   ndata,ndata2,ncent*norbc,ncuspc*(nspin2-nspin1+1),nfockc,ncnstr*(nspin2-nspin1+1)

!     if(mbasis_ctype.gt.0) then   ! Needed since nbasis_ctype(it) is not initialized for quantum rings
!       call alloc ('imnbas', imnbas, ncent)
!       imnbas(1)=1
!       do 36 i=1,ncent-1
!         it=iwctype(i)
!  36     imnbas(i+1)=imnbas(i)+nbasis_ctype(it)
!     endif

      if(ijas.eq.2 .and. iabs(icusp2).ge.2) then
        do isp=nspin1,nspin2
          ishft=ncuspc*(isp-nspin1)
          if(nspin2.eq.2 .and. (nup-ndn).ne.0) a1(2,2,1)=a1(2,1,1)
          if(nspin2.eq.3 .and. isp.eq.nspin2) a1(2,3,1)=((ndn-nup)*a1(2,1,1)+(nup-1)*a1(2,2,1))/(ndn-1)
          call cuspcheck2(scalek(1),a1(1,isp,1),a2(1,isp,1), diff(ndata+ishft+1),isp,nspin1,ncuspc,1)
        enddo
      elseif(ijas.eq.3) then
        call cuspcheck3(diff(ndata+1),1)
      endif

! Moved to read_input.f
!     if(icusp2.ge.1.and.ijas.eq.3.and.isc.le.7) call cuspinit3(1)
!     if(icusp2.ge.1.and.ijas.eq.4.and.isc.le.7) call cuspinit4(1)

      ishft=ncuspc*(nspin2-nspin1+1)+nfockc

! If we are doing an all-electron calculation or using a chemistry pseudopotential, some of which
! have a negative divergence, and using analytical basis functions then compute
! the cusp violation for each occupied s-like orbital at each nucleus.
! If icusp.ge.0 then impose the cusp condition.  In that case the computed cusp violation must be 0.
!     if((nloc.eq.0. .or. nloc.eq.6) .and. numr.le.0) then
      if(icusp.ge.0) then
! Calculate coefs to construct the piece of the orbital that comes
! from basis fns that are related by symmetry.  This is needed to
! impose cusp conditions when there is more than one atom.
        call equiv_bas

        call cuspco(diff(ndata+ishft+1),1)
      endif

      call alloc('parm', parm, nparm)

!     write(6,'(''coef='')')
!     do iorb=1,norb
!       write(6,'(1000f9.6)') (coef(ibasis,iorb,1),ibasis=1,nbasis)
!     enddo

      do iparm=1,nparml
        parm(iparm)=coef(iwbasi(iparm),iworb(iparm),1)
      enddo
      do iparm=1,nparme
        parm(nparml+iparm)=zex(iwbase(iparm),1)
      enddo
      do iparm=1,nparmcsf
        parm(nparml+nparme+iparm)=csf_coef(iwcsf(iparm),1)
      enddo
      if(nparms.eq.1) parm(nparml+nparme+nparmcsf+1)=scalek(1)
      if(nparmg.eq.1) parm(nparml+nparme+nparmcsf+nparms+1)=a21
      if(ijas.eq.1) then
        if(nparmj.ge.1) parm(nparm)=cjas2(1)
      elseif(ijas.eq.2) then
        ntmp=nparmj
        do isp=nspin1,nspin2
          do iparm=1,nparma(isp)
            parm(nparm-ntmp+iparm)=a1(iwjasa(iparm,isp),isp,1)
          enddo
          ntmp=ntmp-nparma(isp)
        enddo
        do isp=nspin1,nspin2
          do iparm=1,nparmb(isp)
            parm(nparm-ntmp+iparm)=a2(iwjasb(iparm,isp),isp,1)
          enddo
          ntmp=ntmp-nparmb(isp)
        enddo
      elseif(ijas.eq.3) then
        ntmp=nparmj
        do iparm=1,nparma(1)
          parm(nparm-ntmp+iparm)=a(iwjasa(iparm,1),1)
        enddo
        ntmp=ntmp-nparma(1)
        do isp=nspin1,nspin2b
          do iparm=1,nparmb(isp)
            parm(nparm-ntmp+iparm)=b(iwjasb(iparm,isp),isp,1)
          enddo
          ntmp=ntmp-nparmb(isp)
        enddo
        do it=1,nctype
          do iparm=1,nparmc(it)
            parm(nparm-ntmp+iparm)=c(iwjasc(iparm,it),it,1)
          enddo
          ntmp=ntmp-nparmc(it)
        enddo
        if(ifock.gt.0) then
          do it=1,nctype
            do iparm=1,nparmf(it)
              parm(nparm-ntmp+iparm)=fck(iwjasf(iparm,it),it,1)
            enddo
            ntmp=ntmp-nparmf(it)
          enddo
        endif
      elseif(ijas.ge.4.and.ijas.le.6) then
        ntmp=nparmj
        do it=1,nctype
          do iparm=1,nparma(it)
            parm(nparm-ntmp+iparm)=a4(iwjasa(iparm,it),it,1)
          enddo
          ntmp=ntmp-nparma(it)
        enddo
        do isp=nspin1,nspin2b
          do iparm=1,nparmb(isp)
            parm(nparm-ntmp+iparm)=b(iwjasb(iparm,isp),isp,1)
          enddo
          ntmp=ntmp-nparmb(isp)
        enddo
        do it=1,nctype
          do iparm=1,nparmc(it)
            parm(nparm-ntmp+iparm)=c(iwjasc(iparm,it),it,1)
          enddo
          ntmp=ntmp-nparmc(it)
        enddo
      endif

      open(1,file='mc_configs',status='old',form='formatted')
      rewind 1
      call alloc ('x', x, 3, nelec, ndata)
      do k=1,ndata
        if(irewgt.ne.0) then
          read(1,*,err=999) ((x(i,j,k),i=1,ndim),j=1,nelec),isgn,psio(k),eold(k)
        else
          read(1,*,err=999) ((x(i,j,k),i=1,ndim),j=1,nelec)
        endif
      enddo
      close(1)

      if(istrch.ge.1.or.irewgt.ne.0) then
        do idata=1,ndata
          dvpdv(idata)=one
        enddo
      endif
!  If nucleii have been moved, move electrons
      if(istrch.ge.1) call strech_fit !JT

!  Set the quadrature points for nonlocal pseudopotential evaluation
      if(nloc.gt.0) call rotqua

      if(iopt.le.1) then

!       if(ncalls.gt.0) call zxssq2(func,ndata2,nparm,nsig,zero,zero,ncalls,iopt,popt,parm,err2,diff,xjac,ndata,xjtj,work,infer,ier)

!       write(6,'(/,''nsig, infer, ier='',3i5)') nsig,infer,ier
!       write(6,'(''if infer=1 then the param. estim. agree on 2 succesive iter. to'',i3,'' sig. digits'',/, &
!    &  ''if infer=2 the residual sum of squares estimates have rel diff lt'',d12.4,/,''if infer=3 the norm of the approx. gradient is lt'',d12.4)') nsig,eps1,del
!       write(6,'(/,''norm of grad, marquardt scaling param. ='',2d12.4)') work(1),work(4)
!       write(6,'(''no of func eval., no of iter, no of sig. digits in param. ='',2f6.0,f5.1)') work(2),work(5),work(3)

       else

! Of the nparm parameters, quench computes the elements of the Jacobian matrix
! for the first nparm-nanalytic, numerically in a call to func
! for the last nanalytic, analytically in a call to jacobian.
        epsp=0.01d0*10.d0**(-nsig)
        epsg=epsp
        epsch2=epsp
        ichange=1
        rot_wt=1.d0
        eps_diff=1.d-15
        if(ianalyt_lap.eq.0) eps_diff=1.d-8
        if(analytic_jacobian) then
          nanalytic=nparmcsf+nparmj+nparms
          write(6,'(''Number of parameters optimized analytically='',i4)') nanalytic
        else
          nanalytic=0
        endif
        call systemflush(6)

#ifndef NOQUENCH
! Warning tmp
        write(6,'(''ndata,ndata2='',9i5)') ndata,ndata2
!       ndata2=ndata
        call quench(func, jacobian, nanalytic, parm, pmarquardt, tau_fit, noutput, &
     &  nstep_fit, ndata2, nparm, ipr_opt, diff, err2, epsg, epsp, epsch2, converg, mesg, ibold, cholesky, rot_wt, eps_diff)
        call systemflush(6)
#else
        stop 'needs quench library'
#endif

        if(converg) then
          write(6,'(''chisq='',es13.6,i5,'' func evals, convergence: '',a10)') err2,icalls,mesg
        else
          write(6,'(''chisq='',es13.6,i5,'' func evals, no convergence'')') err2,icalls
        endif
        call systemflush(6)

      endif

! The foll. call to func is no longer needed since quench takes care of it.
! Warning: I am almost but not completely sure about this. It may be needed for dependent parms.
!     dum_func=func(ndata2,nparm,parm,diff,0)
!     dum_func=func(ndata2,nparm,parm,diff,1)

      if(icusp2.ge.1 .and. isc.le.7) then
        if(ijas.eq.2) then
          do isp=nspin1,nspin2
            call cuspexact2(a1(1,isp,1),a2(1,isp,1))
          enddo
        endif
        if(ijas.eq.3) call cuspexact3(1)
        if(ijas.ge.4.and.ijas.le.6) call cuspexact4(1,1)
      endif

      if(iabs(icusp2).ge.2) then
        if(ijas.eq.2) then
          do isp=nspin1,nspin2
            ishft=ncuspc*(isp-nspin1)
            if(nspin2.eq.2 .and. (nup-ndn).ne.0) a1(2,2,1)=a1(2,1,1)
            if(nspin2.eq.3 .and. isp.eq.nspin2) a1(2,3,1)=((ndn-nup)*a1(2,1,1)+(nup-1)*a1(2,2,1))/(ndn-1)
            call cuspcheck2(scalek(1),a1(1,isp,1),a2(1,isp,1),diff(ndata+ishft+1),isp,nspin1,ncuspc,1)
          enddo
        elseif(ijas.eq.3) then
          call cuspcheck3(diff(ndata+1),1)
        endif
      endif
      ishft=ncuspc*(nspin2-nspin1+1)+nfockc
!     if((nloc.eq.0. .or. nloc.eq.6) .and. numr.le.0) call cuspco(diff(ndata+ishft+1),1)
      if(icusp.ge.0) call cuspco(diff(ndata+ishft+1),1)

      do i=1,necn
        coef(iebasi(1,i),ieorb(1,i),1)=sign(one,dfloat(ieorb(2,i))) * coef(iebasi(2,i),iabs(ieorb(2,i)),1)
      enddo
      do i=1,nebase
        zex(iebase(1,i),1)=zex(iebase(2,i),1)
      enddo

! If the exponents of the basis functions are optimized, recompute normalizations
      if(nparme.gt.0) then
        if(ibasis.eq.3) then
          call basis_norm_dot(0)
        else
          call basis_norm(1,0)
        endif
      endif

! Reset norb to the total number of orbs rather than the number of orbs in at least one determinant
      call object_provide ('orb_tot_nb')
      norb=orb_tot_nb
      write(6,'(''norb reset to total number of orbs, ,norb='',9i5)') norb

      if(inum_orb.eq.0) then
        do ib=1,nbasis
          do iorb=1,norb
            coef(ib,iorb,1)=coef(ib,iorb,1)/anorm(ib)
          enddo
        enddo
      endif

! Pivot among orbitals of the same symmetry.
! Orbitals used to pivot other orbitals must appear in all the determinants (core states).
! Orbs of the same symmetry have the same abs(ipivot).
! Positive means use it to pivot other orbs of same sym; negative means it gets pivoted by other orbs of that sym.

      if(nparml.gt.0.or.nparme.gt.0) then ! lcao or exponential params are optimized
! Pivot to make some coefs. zero
        do iorb=1,norb
          indexi=ipivot(iorb)
          if(indexi.ne.0) then
            comax=zero
            do ib=1,nbasis
              if(dabs(coef(ib,iorb,1)).gt.comax) then
                comax=dabs(coef(ib,iorb,1))
                jb=ib
              endif
            enddo ! ib
            do jorb=1,norb
             indexj=ipivot(jorb)
             ipij=iabs(indexi)-iabs(indexj)
             if((ipij.eq.0.and.indexi.gt.0).and.iorb.ne.jorb) then
               rat=coef(jb,jorb,1)/coef(jb,iorb,1)
               do kbasis=1,nbasis
                 coef(kbasis,jorb,1)=coef(kbasis,jorb,1) - rat*coef(kbasis,iorb,1)
               enddo
             endif
            enddo ! jorb
          endif
        enddo ! iorb

! Make the largest in absolute magnitude coef be 1.  This renormalizes the
! determinant, so adjust cdet (in old code) or csf_coef (in new code) to compensate.
! Warning: The foll. as implemented, is not always correct.  
! It does not work for N2/jas4/cvb1/shci_real_rotated_csf061.  
! Also, from before I have a note that say that it does does not work for N, but I have not checked that.
! Earlier note:  It assumes that if there are several determinants that make up a CSF, then each determinant has the
! same orbitals (though of course divided differently among up and down determinants).
! For example if we consider for the N atom the CSF that comes from exciting
! (a p->d for up-spin and a s->p for dn-spin) then the determinants that make up the
! CSF do not all have the same orbitals.  However, they have the same number of p and d
! orbitals and so if the same renormalization is applied to all the p's and all the d's
! then it is correct I think.
!       do iorb=1,norb
!         comax=zero
!         do ibas=1,nbasis
!           if(dabs(coef(ibas,iorb,1)).gt.dabs(comax)) then
!             comax=coef(ibas,iorb,1)
!           endif
!         enddo ! ibas
!         do icsf=1,ncsf
!           idet=iwdet_in_csf(1,icsf)
!           do j=1,nelec
!             if(iworbd(j,idet).eq.iorb) csf_coef(icsf,1)=csf_coef(icsf,1)*comax
!           enddo
!         enddo
!         do ibas=1,nbasis
!           coef(ibas,iorb,1)=coef(ibas,iorb,1)/comax
!         enddo
!       enddo ! iorb

      endif ! lcao or exponential params are optimized

! Renormalize so first CSF coef. is 1.  No good reason to do this.
!     if((nparml.gt.0.or.nparme.gt.0.or.nparmcsf.gt.0).and.csf_coef(1,1).ne.0.d0) then
!       term=one/csf_coef(1,1)
!       do 240 icsf=1,ncsf
! 240     csf_coef(icsf,1)=csf_coef(icsf,1)*term
!     endif

      if(ijas.ge.4) then
        call set_scale_dist(0,1)
        call plot_jas
      endif

      write(6,'(/,''Final parameters:'')')

      if(iperiodic.eq.0) then
        if(nbasis.gt.9999) stop 'nbasis > 9999 in fit: increase i4 below'
        do iorb=1,norb
          if(iorb.eq.1) then
!           write(fmt,'(''(''i3,''f14.8,\'\' ((coef(j,i),j=1,nbasis),i=1,norb)\'\')'')') nbasis
            write(fmt,'(''(''i3,''f14.8,a38)'')') nbasis
            write(6,fmt) (coef(ib,iorb,1),ib=1,nbasis),' ((coef(j,i),j=1,nbasis),i=1,norb)'
          else
!           write(fmt,'(a1,i3,a6)') '(',nbasis,'f14.8)'
            write(fmt,'(''(''i3,''f14.8)'')') nbasis
            write(6,fmt) (coef(ib,iorb,1),ib=1,nbasis)
          endif
        enddo

!       write(fmt,'(''(''i3,''f14.8,\'\' (zex(i),i=1,nbasis)\'\')'')') nbasis
!       write(6,fmt) (zex(ib,1),ib=1,nbasis)
        write(fmt,'(''(''i3,''f14.8,a20)'')') nbasis
        write(6,fmt) (zex(ib,1),ib=1,nbasis),' (zex(i),i=1,nbasis)'

        if(ndet.gt.9999) stop 'ndet > 9999 in fit: increase i4 below'
        write(fmt,'(''(''i4,''f12.8,a)'')') ncsf
        write(6,fmt) (csf_coef(icsf,1),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      endif

!     write(6,'(/,''Final Jastrow parameters:'')')

      write(6,'(2f13.8,'' scalek,a21'')') scalek(1),a21
      call systemflush(6)

      if(ijas.eq.1) write(6,'(''cjas2='',f12.8)') cjas2(1)

      if(ijas.eq.2) then

        do isp=nspin1,nspin2
          if(nparm_read.gt.0) then
            write(fmt,'(''(''i2,''f16.8,a27'')') nparm_read
           else
            write(fmt,'(''(a27)'')') nparm_read
          endif
          write(6,fmt) (a1(i,isp,1),i=1,nparm_read),' (a(iparmj),iparmj=1,nparma)'
        enddo
        do isp=nspin1,nspin2
          if(nparm_read.gt.0) then
            write(fmt,'(''(''i2,''f19.11,a27'')') nparm_read
          else
            write(fmt,'(''(a28)'')')
          endif
          write(6,fmt) (a2(i,isp,1),i=1,nparm_read),' (b(iparmj),iparmj=1,nparmb)'
        enddo

      elseif(ijas.ge.3.and.ijas.le.6) then ! ijas /= 2
        if(ijas.eq.3) then
          if(nparm_read.gt.0) then
            write(fmt,'(''(''i2,''f19.11,a28)'')') nparm_read
          else
            write(fmt,'(''(a28)'')')
          endif
          write(6,fmt) (a(i,1),i=1,nparm_read),' (a(iparmj),iparmj=1,nparma)'
          nparmc_read=(nord**3+5*nord)/6+nord**2+nord
        elseif(ijas.ge.4.and.ijas.le.6) then
          nparma_read=2+max(0,norda-1)
          do it=1,nctype
            if(nparma_read.gt.0) then
              write(fmt,'(''(''i2,''f19.11,a28)'')') nparma_read
             else
              write(fmt,'(''(a28)'')')
            endif
            write(6,fmt) (a4(i,it,1),i=1,nparma_read),' (a(iparmj),iparmj=1,nparma)'
          enddo
          nparm_read=2+max(0,nordb-1)
          nparmc_read=nterms4(nordc)
        endif

        do isp=nspin1,nspin2b
          if(nparm_read.gt.0) then
            write(fmt,'(''(''i2,''f19.11,a28)'')') nparm_read
          else
            write(fmt,'(''(a28)'')')
          endif
          write(6,fmt) (b(i,isp,1),i=1,nparm_read),' (b(iparmj),iparmj=1,nparmb)'
        enddo

        do it=1,nctype
          if(nparmc_read.gt.0) then
            write(fmt,'(''(''i2,''f19.11,a28)'')') nparmc_read
          else
            write(fmt,'(''(a28)'')')
          endif
          write(6,fmt) (c(i,it,1),i=1,nparmc_read),' (c(iparmj),iparmj=1,nparmc)'
        enddo

        if(ifock.gt.0) then
          do it=1,nctype
            write(6,'(15f19.11,'' (f(iparmj),iparmj=1,nparmf)'')') (fck(i,it,1),i=1,15)
          enddo
        endif
      endif
      call systemflush(6)

      if(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17) write(6,'(2f12.8,'' cutjas_en,cutjas_ee'')') cutjas_en,cutjas_ee

      write(6,'(''cusp penalty weighted by'',f10.5)') cuspwt
      call systemflush(6)

      if(ncalls.gt.0 .or. iopt.gt.1) errc=dsqrt(err2/ndata2)

      if((iopt.le.1 .and. ncalls.le.0) .or. ndata.ne.ndata2) then
        err2=zero
        do i=1,ndata
          err2=err2+diff(i)**2
        enddo
      endif
      err=dsqrt(err2/ndata)

      if(iopt.le.1 .and. ncalls.le.0 .and. ndata.ne.ndata2) then
        errc=err2
        nx=ndata+ncuspc*(nspin2-nspin1+1)+nfockc+ncent*norbc
        do i=ndata+1,nx
          errc=errc+(diff(i)*cuspwt)**2
        enddo
        do i=nx+1,ndata2
          errc=errc+diff(i)**2
        enddo
        errc=dsqrt(errc/ndata2)
      endif

      if(iabs(icusp2).ge.2) then
        if(mod(icusp2,3).eq.0) then
          write(6,'(''rms fit error,errc/'',i5,2es12.4,'' ndata,ndata2'',2i5)') iabs(icusp2)-2,err,errc,ndata,ndata2
        else
          write(6,'(''rms fit error,errc*'',i5,2es12.4,'' ndata,ndata2'',2i5)') nint(cuspwt),err,errc,ndata,ndata2
        endif
       else
        write(6,'(''rms fit error,errc='',2es12.4,'' ndata,ndata2='',2i5)') err,errc,ndata,ndata2
      endif

      if(ipr_opt.ge.2) write(6,'(/,''config #,   etrial,       efit,      diff,       relerr,    wght,     rmin,      rmax,      r12 min'')')
      call systemflush(6)
      if(ipr_opt.ge.2) rewind 2

      rminav=zero
      rmaxav=zero
      r12minav=zero
      eav=zero
      reldif=zero
      relerr=zero
      wtone=zero
      ovltop=zero
      ovlbot=zero
      do i=1,ndata
        yfit=eguess+uwdiff(i)
        rmin=1.d99
        rmax=0
        r12min=1.d99
        do j=1,nelec
          r1=0
          do k=1,ndim
            r1=r1+x(k,j,i)**2
          enddo
          rmin=dmin1(r1,rmin)
          rmax=dmax1(r1,rmax)
          do jj=1,j-1
            r12=0
            do k=1,ndim
              r12=r12+(x(k,j,i)-x(k,jj,i))**2
            enddo
            r12min=dmin1(r12,r12min)
          enddo
        enddo ! nelec
        rmin=dsqrt(rmin)
        rmax=dsqrt(rmax)
        r12min=dsqrt(r12min)
        rminav=rminav+rmin
        rmaxav=rmaxav+rmax
        r12minav=r12minav+r12min
        if(irewgt.ne.0) then
          eav=eav+yfit*wght(i)
          reldif=reldif+yfit*wght(i)-eold(i)
          relerr=relerr+(yfit-eold(i))**2*wght(i)
          if(wght(i).le.one) then
            wtone=wtone+wght(i)**2
          else
            wtone=wtone+one/wght(i)**2
          endif
          ovltop=ovltop+dsqrt(wght(i))
          ovlbot=ovlbot+wght(i)
        else
          eav=eav+yfit
          reldif=reldif+yfit-eold(i)
          relerr=relerr+(yfit-eold(i))**2
        endif
        if(ipr_opt.ge.2) write(6,'(i5,2f15.5,es12.4,2f9.3,2f12.6,f9.3)') &
     &  i,eguess,yfit, uwdiff(i),yfit-eold(i),wght(i),rmin,rmax,r12min
         if(ipr_opt.ge.2) then
           write(fmt,'(''('',i6,''f13.8,i3,es12.4,f12.5,2es12.4)'')')ndim*nelec
!          write(2,fmt)((x(k,j,i),k=1,ndim),j=1,nelec),psid(i),psij(i),psid(i)*exp(psij(i)),yfit
           write(2,fmt)((x(k,j,i),k=1,ndim),j=1,nelec),int(sign(1.d0,psid(i))),log(dabs(psid(i)))+psij(i),yfit,psid(i),psij(i)
         endif
      enddo ! ndata

      if(iabs(icusp2).ge.2 .and. ipr_opt.ge.0) then
        do i=ndata+1,ndata2
          write(6,'(i5,20x,es13.5)') i,diff(i)
        enddo
      endif

      rminav=rminav/ndata
      rmaxav=rmaxav/ndata
      r12minav=r12minav/ndata
      eav=eav/ndata
      reldif=reldif/ndata
      relerr=dsqrt(max(0.d0,relerr/ndata-reldif*reldif))
      wtone=dsqrt(wtone/ndata)
      if(ovlbot.ne.zero) ovlap=ovltop/dsqrt(ndata*ovlbot)
      write(6, '(''ndata, eguess, eav, rms fit error, relerr, wtone, ovlap, rminav, rmaxav, r12minav'', &
     &i5,2f12.6,2f10.6,2x,f6.3,f7.4,2x,9f6.3)') &
     &ndata,eguess,eav,err,relerr, wtone,ovlap, rminav,rmaxav,r12minav
      if(analytic_jacobian.and.irewgt.ne.0) write(6,'(''*** Warning: analytic derivs not correctly implemented when points are reweighted'')')
      if(eguess.gt.eav+0.1*err) write(6,'(''*** Warning: eguess is greater than average energy on fit pts + 0.1*std. Decrease eguess'')')
      if(eguess.lt.eav-0.3*err) write(6,'(''*** Warning: eguess is less than average energy on fit pts - 0.3*std. Increase eguess'')')

      if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt.gt.0 .and. ijas.eq.2) then
        ishft=ndata+ncuspc*(nspin2-nspin1+1)+nfockc+ncent*norbc+1
        icalcul_diff=1
        do isp=nspin1,nspin2
          call checkjas2(scalek(1),isp,ncnstr,diff(ishft),ipr_opt,1,icalcul_diff)
          ishft=ishft+ncnstr
        enddo
      endif

      call my_second(2,'all   ')
      call systemflush(6)

      return

  999 write(6,'(''Error reading mc_configs, config.'',i5)') k
      stop 'Error reading mc_configs, possibly not enough configs.'

      end

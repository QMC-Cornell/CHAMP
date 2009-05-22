      subroutine fit
c Written by Cyrus Umrigar
c Warning: Moved some of the input to read_input.f and commented it out here
c          Need to make sure those variables get passed back to here via commons.
c Uses the variance minimization algorithm described in:
c 1) Optimized Trial Wave Functions for Quantum Monte Carlo Calculations,
c    C.J. Umrigar, K.G. Wilson and J.W. Wilkins, Phys. Rev. Lett.,60, 1719 (1988).
c 2) A Method for Determining Many-Body Wavefunctions, C.J. Umrigar, K.G. Wilson and J.W. Wilkins,
c    {\it Computer Simulation Studies in Condensed Matter Physics: Recent Developments,}
c    ed. by D.P. Landau, K.K. Mon and H.B. Schuttler (Springer-Verlag 1988).
c 3) Two Aspects of Quantum Monte Carlo: Determination of Accurate Wavefunctions and
c    Determination of Potential Energy Surfaces of Molecules, C.J. Umrigar,
c    Int. J. Quant. Chem. Symp., 23, 217 (1989).
      use control_mod
      use allocations_mod
      use atom_mod
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
      implicit real*8(a-h,o-z)
c     character*16 mode
      character*80 fmt
c      character*30 section
      character*10 mesg
      logical converg,analytic,cholesky
      external func,jacobian

!JT      include '../vmc/vmc.h'
!JT      include 'fit.h'
!JT      include '../vmc/pseudo.h'
!JT      include '../vmc/numbas.h'
!JT      include '../vmc/force.h'
c     parameter(MXJTJ=(MPARM*(MPARM+1))/2,MWORK=4*MDATA+5*MPARM+MXJTJ)
!JT      parameter(zero=0.d0,one=1.d0,two=2.d0,four=4.d0)

!JT      common /dim/ ndim
      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Zfock
!JT      common /contrl_per/ iperiodic,ibasis
!JT      common /contrl_opt/ nparm,nsig,ncalls,iopt,ipr_opt
!JT      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
!JT     &,ifock,i3body,irewgt,iaver,istrch
!JT     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
c     common /contr3/ mode
      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /atom2/ dist_cent(MCENT)

!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
!JT      common /basis/ zex(MBASIS,MWF),betaq
!JT     &,n1s(MCTYPE)
!JT     &,n2s(MCTYPE),n2p(-1:1,MCTYPE)
!JT     &,n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
!JT     &,n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE),n4f(-3:3,MCTYPE)
!JT     &,n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
!JT     &,n5g(-4:4,MCTYPE)
!JT     &,n6d(-2:2,MCTYPE),n6f(-3:3,MCTYPE),n6g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
!JT     &,n7g(-4:4,MCTYPE),n7h(-5:5,MCTYPE),n7i(-6:6,MCTYPE)
!JT     &,n8i(-6:6,MCTYPE),n8j(-7:7,MCTYPE)
!JT     &,n9k(-8:8,MCTYPE)
!JT     &,n10l(-9:9,MCTYPE)
!JT     &,n11m(-10:10,MCTYPE)
!JT     &,n12n(-11:11,MCTYPE)
!JT     &,n13o(-12:12,MCTYPE)
!JT     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
!JT      common /basis2/ zex2(MRWF,MCTYPE,MWF),n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
!JT     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
!JT     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),iwrwf2(MBASIS)
!JT      common /basisnorm/ anorm(MBASIS)
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
!JT      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
!JT     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
c     common /wcsf/ frac(ICX,MDET),icsf(ICSFX)

!JT      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar1/ cjas1(MWF),cjas2(MWF)
      common /jaspar2/ a1(MPARMJ,3,MWF),a2(MPARMJ,3,MWF)
!JT      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
!JT     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
!JT      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
!JT      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
!JT     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
!JT     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
!JT     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
!JT     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
!JT     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
!JT      common /bparm/ nspin2b,nocuspb
      common /ncusp/ norbc,ncuspc,nfockc,nfock,ncnstr

!JT      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
!JT     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
!JT      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
!JT     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
!JT     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
!JT      common /pointer/ npointa(MPARMJ*NCTYP3X)
      common /focsav/ a4sav,a5sav,a6sav,a7sav

      common /fcn_calls/icalls
      common /update/ ichange

      dimension parm(MPARM),diff(MDATA),ipivot(MORB)
c For zxssq
c     dimension,popt(4),xjac(MDATA,MPARM),xjtj(MXJTJ),work(MWORK)

c Inputs are:
c 1) ndata,nparm,ijas,icusp,icusp2,isc,nsig,ncalls,iopt,ipr_opt
c    ndata  = no. of configs to optimize over
c    nparm  = no. of parameters to optimize
c    ijas   = form of wavefunction
c    icusp  = <= -1  zeroth order e-N cusp not imposed
c             >=  0  zeroth order e-N cusp imposed "exactly"
c    icusp2   is for imposing e-N and e-e cusps
c             <= -2 impose cusps via penalty
c    isc      is used for 3 purposes
c             a) for using spq rather than srt variables
c             b) for calling numerical (jastrow_num,psi) rather than analytical (jastrow)
c                gradient and Laplacian. Analytic if >= -5
c             c) for scaling ri,rj,rij in Pade.
c                -2,-7  (1-exp(scalek(1)*r)/scalek(1)
c                -3,-8  [1-exp{-scalek(1)*r-(scalek(1)*r)**2/2}]/scalek(1)
c                -4,-9  r/(1+scalek(1)*r)
c                -5,-10 r/{1+(scalek(1)*r)**2}**.5
c    nsig     = no. of significant digits required in parameters
c    ncalls   = Max. no. of function calls to make
c    iopt     In fit:
c               Choose between zxssq and quench
c             = 0,1 call zxssq from IMSL (obsolete)
c             = 2   call quench written by Peter Nightingale and Cyrus Umrigar
c             = 0 Do not check for strict downward descent in zxssq
c             = 1 Strict downward descent in zxssq
c               zxssq is obsolete so, if in fit mode, we reset iopt to 2 in read_input
c             In vmc:
c               Choose between Newton, linear and perturbation theory
c             = 1 linear method
c             = 2 modified Newton method
c             = 3 perturbation theory
c    ipr_opt  <= -2  Minimal print out
c             >= -1  Print cusp monotonicity conditions
c             >=  0  Print configs and errors on 6
c             >=  2  Print out configs and wavefunction on 2
c 2) i3body,irewgt,iaver,istrch
c 3) ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
c nparml   no. of linear    parameters
c nparme   no. of exponent  parameters
c nparmd   no. of determin  parameters (obsolete)
c nparmcsf no. of CSF       parameters
c nparms   no. of scale     parameters for good Jastrow (scalek(1))
c nparmg   no. of GAM       parameters (a21)
c nparmj   no. of Jastrow   parameters
c nparma   no. of Jastrow a parameters in complicated Jastrow
c nparmb   no. of Jastrow b parameters in complicated Jastrow
c nparmc   no. of Jastrow c parameters in complicated Jastrow
c nparmf   no. of Jastrow fck parameters in complicated Jastrow

      call my_second(0,'begin ')

c     mode='fit         '
!JT      call read_input
      write(6,'(''returned to fit from read_input'')')

      call common_allocations

c Calculate distances of atoms from center for use in cusorb
      do 2 i=1,ncent
        dist_cent(i)=0
        do 1 k=1,ndim
    1     dist_cent(i)=dist_cent(i)+cent(k,i)**2
    2   dist_cent(i)=sqrt(dist_cent(i))

c Save values of a(4) thru a(7), in order to adiabatically go from
c 4th order Pade without Fock terms to one with
      if(ifock.gt.0.and.ijas.eq.2) then
        a4sav=a1(4,1,1)
        a5sav=a1(5,1,1)
        a6sav=a1(6,1,1)
        a7sav=a1(7,1,1)
      endif
      if(ifock.eq.4.and.ijas.eq.3) then
        do 3 it=1,nctype
          a4sav=c(4,it,1)
          a5sav=c(5,it,1)
          a6sav=c(7,it,1)
          a7sav=c(9,it,1)
    3     call scale20(1,it)
      endif

c Optimization section
c     read(5,*) section
c     write(6,'(a)') section

c     read(5,*)
c     read(5,*) ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt
c     read(5,*) i3body,irewgt,iaver,istrch
c     read(5,*) ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,
c    &idbdu,idbdt

c     if(ndata.gt.MDATA) stop 'MDATA exceeded'
c     if(nparm.gt.MPARM) stop 'MPARM exceeded'

c     if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt
c    &.gt.0.and.(ijas.ne.2))
c    &stop 'ipos+...>0 checkjas2 exists only be used with Jastrow2'
c     if(mod(irewgt,100).eq.1) then
c       write(6,*) '**Warning irewgt=1 reset to irewgt=10'
c       irewgt=irewgt+9
c     endif
c     do 4 i=1,ndata
c  4    wght(i)=one

c     write(6,'(/,''ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt=''
c    &,20i4)') ndata,nparm,icusp,icusp2,nsig,ncalls,iopt,ipr_opt
c     write(6,'(''i3body,irewgt,iaver,istrch'',9i5)')
c    &i3body,irewgt,iaver,istrch
c     write(6,'(''ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt'',10i8)')
c    & ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdt,idbdu

c     if(isc.eq.2) write(6,'(
c    &''dist scaled r=(1-exp(-scalek*r))/scalek'')')
c     if(isc.eq.3) write(6,'(
c    &''dist scaled r=(1-exp(-scalek*r-(scalek*r)**2/2))/scalek'')')
c     if(isc.eq.4) write(6,'(
c    &''dist scaled r=r/(1+scalek*r)'')')
c     if(isc.eq.5) write(6,'(
c    &''dist scaled r=r/(1+(scalek*r)**2)**.5'')')

c     if(ijas.eq.1) write(6,'(''Conventional Jastrow'')')
c     if(ijas.eq.2) write(6,'(''Exp. Pade + non-anal terms'')')
c     if(ijas.eq.3) write(6,'(''Standard form'')')
c     if(ijas.eq.4) write(6,'(''New transferable standard form 4'')')
c     if(ijas.eq.5) write(6,'(''New transferable standard form 5'')')
c     if(ijas.eq.6) write(6,'(''New transferable standard form 6'')')

c     if(icusp.ge.0) write(6,'(''Nuclear cusp constraint is imposed'')')

c     read(5,*) (lo(iorb),iorb=1,norb)
c     write(6,'(''lo='',20i3)') (lo(iorb),iorb=1,norb)
cc    read(5,*) (n(ib),l(ib),ib=1,nbasis)
cc    write(6,'(''n,l='',20(2i3,1x))') (n(ib),l(ib),ib=1,nbasis)

c     if(ijas.le.3) then
c       na1=nspin1
c       na2=nspin2
c      else
c       na1=1
c       na2=nctype
c     endif
c     read(5,*) nparml,(nparma(ia),ia=na1,na2),
c    &(nparmb(isp),isp=nspin1,nspin2b),(nparmc(it),it=1,nctype),
c    &(nparmf(it),it=1,nctype),nparmd,nparms,nparmg

c     if(ijas.ge.4.and.ijas.le.6) then
c       do 5 it=1,nctype
c         if(numr.le.0) then
ccAll-electron with analytic slater basis
c           if((norda.eq.0.and.nparma(it).gt.0)
c    &      .or.(norda.gt.0 .and. nparma(it).gt.norda+1)) then
c             write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
c             stop 'nparma too large for norda in all-electron calculation'
c           endif
c          else
ccPseudopotential with numerical basis (cannot vary a(1) or a(2)
c           if(norda.eq.1) stop 'makes no sense to have norda=1 for numr>0'
c           if((norda.eq.0.and.nparma(it).gt.0)
c    &      .or.(norda.gt.0 .and. nparma(it).gt.norda-1)) then
c             write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
c             stop 'nparma too large for norda in pseudopot calculation'
c           endif
c         endif
c         if(isc.le.7 .and.
c    &       ((nordc.le.2.and.nparmc(it).gt.0)
c    &    .or.(nordc.eq.3.and.nparmc(it).gt.2)
c    &    .or.(nordc.eq.4.and.nparmc(it).gt.7)
c    &    .or.(nordc.eq.5.and.nparmc(it).gt.15)
c    &    .or.(nordc.eq.6.and.nparmc(it).gt.27)
c    &    .or.(nordc.eq.7.and.nparmc(it).gt.43))) then
c           write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
c           stop 'nparmc too large for nordc in J_een with cusp conds'
c         endif
c         if(isc.gt.7 .and.
c    &       ((nordc.le.1.and.nparmc(it).gt.0)
c    &    .or.(nordc.eq.2.and.nparmc(it).gt.2)
c    &    .or.(nordc.eq.3.and.nparmc(it).gt.6)
c    &    .or.(nordc.eq.4.and.nparmc(it).gt.13)
c    &    .or.(nordc.eq.5.and.nparmc(it).gt.23)
c    &    .or.(nordc.eq.6.and.nparmc(it).gt.37)
c    &    .or.(nordc.eq.7.and.nparmc(it).gt.55))) then
c           write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
c           stop 'nparmc too large for nordc without cusp conds'
c         endif
c   5   continue
ccFor the b coefs. we assume that b(1) is fixed by the cusp-cond.
c       do 6 isp=1,nspin1,nspin2b
c           if((nordb.eq.0.and.nparmb(isp).gt.0).or.(nordb.gt.0 .and. nparmb(isp).gt.nordb)) then
c             write(6,'(''isp,nordb,nparmb(isp)'',3i5)') isp,nordb,nparmb(isp)
c             stop 'nparmb too large for nordb'
c           endif
c   6   continue
c     endif

cccompute nparmj and nparme
c     nparmj=0
c     npointa(1)=nparmj
c     do 7 ia=na1,na2
c       if(ia.gt.1) npointa(ia)=npointa(ia-1)+nparma(ia-1)
c   7   nparmj=nparmj+nparma(ia)
c     do 8 isp=nspin1,nspin2b
c   8   nparmj=nparmj+nparmb(isp)
c     npoint(1)=nparmj
c     do 9 it=1,nctype
c       if(it.gt.1) npoint(it)=npoint(it-1)+nparmc(it-1)
c   9   nparmj=nparmj+nparmc(it)+nparmf(it)
c     nparme=nparm-nparml-nparmj-nparmd-nparms-nparmg
c     write(6,'(''No of linear coefs, exponents, Jastrow, det, scale parms varied='',9i5)')
c    &nparml, nparme, nparmj, nparmd, nparms
c     if(nparme.lt.0) stop 'nparme < 0'
c     if(nparme.gt.nbasis) stop 'nparme > nbasis'
c     if(nparme.gt.0 .and. numr.gt.0) stop 'nparme > 0 and numr > 0'
c     if(nparml.lt.0 .or. nparmj.lt.0 .or. nparmd.lt.0 .or. nparms.lt.0 .or.nparmg.lt.0)
c    &stop 'nparm? must be >= 0'
c     if(nparms.gt.1) stop 'nparms must be 0 or 1'

c     read(5,*) (iworb(iparm),iwbasi(iparm),iparm=1,nparml)
c     write(6,'(''lin. coefs. of orbs varied='',10(2i3,2x))')
c    &(iworb(iparm),iwbasi(iparm),iparm=1,nparml)

c     read(5,*) (iwbase(iparm),iparm=1,nparme)
c     write(6,'(''exponents varied='',20i3)') (iwbase(iparm),iparm=1,
c    &nparme)

c     read(5,*) (iwdet(iparm),iparm=1,nparmd)
c     write(6,'(''determinantal coefs varied='',20i3)')
c    &(iwdet(iparm),iparm=1,nparmd)

c     if(ijas.eq.2.or.ijas.eq.3) then
c       write(6,'(''Correl. params. that are varied are:'')')
c       do 14 isp=nspin1,nspin2
c         read(5,*) (iwjasa(iparm,isp),iparm=1,nparma(isp))
c  14     write(6,'(''a: '',30i3)') (iwjasa(iparm,isp),iparm=1,
c    &    nparma(isp))
c       do 16 isp=nspin1,nspin2b
c         read(5,*) (iwjasb(iparm,isp),iparm=1,nparmb(isp))
c  16     write(6,'(''b: '',30i3)') (iwjasb(iparm,isp),iparm=1,
c    &    nparmb(isp))
c      elseif(ijas.ge.4.and.ijas.le.6) then
c       do 18 it=1,nctype
c         read(5,*) (iwjasa(iparm,it),iparm=1,nparma(it))
c  18     write(6,'(''a: '',30i3)') (iwjasa(iparm,it),iparm=1,
c    &    nparma(it))
c       do 20 isp=nspin1,nspin2b
c         read(5,*) (iwjasb(iparm,isp),iparm=1,nparmb(isp))
c  20     write(6,'(''b: '',30i3)') (iwjasb(iparm,isp),iparm=1,
c    &    nparmb(isp))
c     endif
c     if(ijas.ge.3.and.ijas.le.6) then
c       do 25 it=1,nctype
c         read(5,*) (iwjasc(iparm,it),iparm=1,nparmc(it))
c  25     write(6,'(''c: '',60i3)') (iwjasc(iparm,it),iparm=1,
c    &    nparmc(it))
c       if(ifock.gt.0) then
c         do 26 it=1,nctype
c           read(5,*) (iwjasf(iparm,it),iparm=1,nparmf(it))
c  26       write(6,'(''f: '',30i3)') (iwjasf(iparm,it),iparm=1,
c    &      nparmf(it))
c       endif
c     endif

c nparm is read in in read_input.f but is not in common.  So recompute it.  No need.  Now it is in common.
c     nparm=nparml+nparmj+nparmd+nparms+nparmg+nparme
c     write(6,'(''nparm,nparml,nparmj,nparmd,nparms,nparmg,nparme='',9i5)') nparm,nparml,nparmj,nparmd,nparms,nparmg,nparme
      write(6,'(''nparm,nparml,nparmj,nparmcsf,nparms,nparmg,nparme='',9i5)') nparm,nparml,nparmj,nparmcsf,nparms,nparmg,nparme
      write(6,'(''total number of parameters, nparm='',i5)') nparm

      read(5,*) necn,nebase
      write(6,'(''No of linear coefs, exponents set equal='',3i5)') necn,nebase
      if(necn.gt.MPARM) stop 'fit: necn>MPARM'
      if(nebase.gt.MBASIS) stop 'fit: nebase>MBASIS'

      read(5,*) ((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)
      write(6,'(''lin. coefs of orbs set equal='',5(2(2i3,2x),2x))')
     &((ieorb(i,j),iebasi(i,j),i=1,2),j=1,necn)

      do 27 j=1,necn
        do 27 i=1,2
          if(ieorb(i,j).gt.norb) stop 'ieorb(i,j).gt.norb'
   27     if(iebasi(i,j).gt.nbasis) stop 'iebasi(i,j).gt.nbasis'

      read(5,*) ((iebase(i,j),i=1,2),j=1,nebase)
      write(6,'(''expon. set equal='',10(2i3,2x))')
     &((iebase(i,j),i=1,2),j=1,nebase)

      do 28 j=1,nebase
        do 28 i=1,2
   28     if(iebase(i,j).gt.nbasis) stop 'iebase(i,j).gt.nbasis'

c     if(nedet.lt.0) then
c       read(5,*) (icsf(j),j=1,iabs(nedet))
c       write(6,'(''icsf='',50i2)') (icsf(j),j=1,iabs(nedet))
c      else
c       do 31 j=1,nedet
c  31     icsf(j)=2
c     endif

c     read(5,*) ((iedet(i,j),i=1,icsf(j)),j=1,iabs(nedet))
c     write(6,'(''coef. of det. set equal='',10(2i3,2x))')
c    &((iedet(i,j),i=1,icsf(j)),j=1,iabs(nedet))

c     if(nedet.lt.0) then
c       read(5,*) ((frac(i,j),i=1,icsf(j)),j=1,iabs(nedet))
c       write(6,'(''frac='',10f6.1)') ((frac(i,j),i=1,icsf(j)),j=1,iabs(nedet))
c      else
c       do 32 j=1,nedet
c  32     frac(2,j)=1.d0
c     endif
      call systemflush(6)

      read(5,*) (ipivot(j),j=1,norb)
      write(6,'(''ipivot='',10i4)') (ipivot(j),j=1,norb)

      read(5,*) eguess
      write(6,'(''eguess='',f12.6)') eguess

      read(5,*) pmarquardt,tau,noutput,nstep,ibold
      write(6,'(''pmarquardt,tau,ibold '',2g8.2,i3)')
     &pmarquardt,tau,ibold
      read(5,'(2l2)') analytic,cholesky
      write(6,'(''analytic,cholesky'',2l2)') analytic,cholesky

      if(analytic .and. irewgt.ne.0)
     &write(6,'(''*** Warning *** analytic derivs not correctly implemented when points are reweighted'')')
      if(analytic .and. nparmf(1).gt.0) stop 'analytic derivatives not implemented yet for Fock parameters'

c If we are doing analytic derivatives wrt the optimization parameters
c then we should not optimize the 1st a parameter because then coef
c changes to maintain the cusp, but the change of the residues wrt
c the change in coef is not calculated analytically.
c If there is only one atom type, then we could easily modify the program
c to include the 1st a parameter in the numerical rather than the analytic
c set, but it is hard to do this when there is more than one atom type
c because in quench it is assumed that the numerical parameters appear before
c the analytic ones.
      if(analytic .and. icusp.eq.1) then
        do 33 it=1,nctype
          do 33 iparm=1,nparma(it)
   33       if(iwjasa(iparm,it).eq.1) stop 'Does not make sense to optimize linear e-n parameter if one is using analytic jacobian
     & and there is more than one atom type.  For one atom type program needs modification to do this.'
      endif

c norbc = nuclear cusp conditions
      norbc=0
      do 34 iorb=1,norb
c  34   if(lo(iorb).eq.0) norbc=norbc+1
   34   if(lo(iorb).eq.0 .and. .not.analytic) norbc=norbc+1
      if(numr.gt.0) norbc=0

c ncuspc = number of cusp conditions from Jastrow
c nfockc = number of cusp conditions from Fock terms (ijas=3)
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
        ndata2=ndata+ncuspc*(nspin2-nspin1+1)+ncent*norbc+nfockc
        if(mod(icusp2,3).ne.0) then
          cuspwt=    dfloat(iabs(icusp2)-1)
         else
          cuspwt=one/dfloat(iabs(icusp2)-2)
        endif
       else
c We are presently calculating a penalty in func but not in jacobian
c so if analytic is true set ndata2=ndata
        if(analytic) then
          ndata2=ndata
         else
          ndata2=ndata+ncent*norbc
        endif
      endif

c isp=1, it is only to determine ncnstr
      ncnstr=0
      ishft=ndata2+1
      if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt
     &.gt.0 .and. ijas.eq.2)
     &call checkjas2(scalek(1),1,ncnstr,diff(ishft),ipr_opt,0)
      ndata2=ndata2+ncnstr*(nspin2-nspin1+1)

      write(6,'(''No of data points fitted to, # of cusp cond='',3i5)')
     &   ndata,ncuspc*(nspin2-nspin1+1)+ncent*norbc+nfockc,
     &   ncnstr*(nspin2-nspin1+1)
      write(6,'(''ndata2'',i5)') ndata2
      if(ndata2.gt.mdata) stop 'ndata2 exceeds mdata'

      imnbas(1)=1
      do 36 i=1,ncent-1
        it=iwctype(i)
   36   imnbas(i+1)=imnbas(i)+nbasis_ctype(it)

      if(iabs(icusp2).ge.2) then
        if(ijas.eq.2) then
          do 37 isp=nspin1,nspin2
            ishft=ncuspc*(isp-nspin1)
            if(nspin2.eq.2 .and. (nup-ndn).ne.0) a1(2,2,1)=a1(2,1,1)
            if(nspin2.eq.3 .and. isp.eq.nspin2)
     &      a1(2,3,1)=((ndn-nup)*a1(2,1,1)+(nup-1)*a1(2,2,1))/(ndn-1)
   37       call cuspcheck2(scalek(1),a1(1,isp,1),a2(1,isp,1),
     &      diff(ndata+ishft+1),isp,nspin1,ncuspc,1)
         elseif(ijas.eq.3) then
          call cuspcheck3(diff(ndata+1),1)
        endif
      endif

c Moved to read_input.f
c     if(icusp2.ge.1.and.ijas.eq.3.and.isc.le.7) call cuspinit3(1)
c     if(icusp2.ge.1.and.ijas.eq.4.and.isc.le.7) call cuspinit4(1)

      ishft=ncuspc*(nspin2-nspin1+1)+nfockc

      if((nloc.eq.0. .or. nloc.eq.6) .and. numr.le.0) then
c Calculate coefs to construct the piece of the orbital that comes
c from basis fns that are related by symmetry.  This is needed to
c impose cusp conditions when there is more than one atom.
        call equiv_bas

        call cuspco(diff(ndata+ishft+1),1)
      endif

      do 40 iparm=1,nparml
   40   parm(iparm)=coef(iwbasi(iparm),iworb(iparm),1)
      do 42 iparm=1,nparme
   42   parm(nparml+iparm)=zex(iwbase(iparm),1)
c     do 45 iparm=1,nparmd
c  45   parm(nparml+nparme+iparm)=cdet(iwdet(iparm),1)
      do 45 iparm=1,nparmcsf
   45   parm(nparml+nparme+iparm)=csf_coef(iwcsf(iparm),1)
c     if(nparms.eq.1) parm(nparml+nparme+nparmd+1)=scalek(1)
c     if(nparmg.eq.1) parm(nparml+nparme+nparmd+nparms+1)=a21
      if(nparms.eq.1) parm(nparml+nparme+nparmcsf+1)=scalek(1)
      if(nparmg.eq.1) parm(nparml+nparme+nparmcsf+nparms+1)=a21
      if(ijas.eq.1) then
        if(nparmj.ge.1) parm(nparm)=cjas2(1)
       elseif(ijas.eq.2) then
        ntmp=nparmj
        do 49 isp=nspin1,nspin2
          do 48 iparm=1,nparma(isp)
   48       parm(nparm-ntmp+iparm)=a1(iwjasa(iparm,isp),isp,1)
   49     ntmp=ntmp-nparma(isp)
        do 51 isp=nspin1,nspin2
          do 50 iparm=1,nparmb(isp)
   50       parm(nparm-ntmp+iparm)=a2(iwjasb(iparm,isp),isp,1)
   51     ntmp=ntmp-nparmb(isp)
       elseif(ijas.eq.3) then
        ntmp=nparmj
        do 52 iparm=1,nparma(1)
   52     parm(nparm-ntmp+iparm)=a(iwjasa(iparm,1),1)
        ntmp=ntmp-nparma(1)
        do 54 isp=nspin1,nspin2b
          do 53 iparm=1,nparmb(isp)
   53       parm(nparm-ntmp+iparm)=b(iwjasb(iparm,isp),isp,1)
   54     ntmp=ntmp-nparmb(isp)
        do 58 it=1,nctype
          do 57 iparm=1,nparmc(it)
   57       parm(nparm-ntmp+iparm)=c(iwjasc(iparm,it),it,1)
   58     ntmp=ntmp-nparmc(it)
        if(ifock.gt.0) then
          do 60 it=1,nctype
            do 59 iparm=1,nparmf(it)
   59         parm(nparm-ntmp+iparm)=fck(iwjasf(iparm,it),it,1)
   60       ntmp=ntmp-nparmf(it)
        endif
       elseif(ijas.ge.4.and.ijas.le.6) then
        ntmp=nparmj
        do 62 it=1,nctype
          do 61 iparm=1,nparma(it)
   61       parm(nparm-ntmp+iparm)=a4(iwjasa(iparm,it),it,1)
   62     ntmp=ntmp-nparma(it)
        do 64 isp=nspin1,nspin2b
          do 63 iparm=1,nparmb(isp)
   63       parm(nparm-ntmp+iparm)=b(iwjasb(iparm,isp),isp,1)
   64     ntmp=ntmp-nparmb(isp)
        do 66 it=1,nctype
          do 65 iparm=1,nparmc(it)
   65       parm(nparm-ntmp+iparm)=c(iwjasc(iparm,it),it,1)
   66     ntmp=ntmp-nparmc(it)
      endif

      open(1,file='mc_configs',status='old',form='formatted')
      rewind 1
      do 67 k=1,ndata
        if(irewgt.ne.0) then
          read(1,*,err=999) ((x(i,j,k),i=1,ndim),j=1,nelec),isgn,psio(k),eold(k)
         else
          read(1,*,err=999) ((x(i,j,k),i=1,ndim),j=1,nelec)
        endif
   67 continue
      close(1)

      if(istrch.ge.1.or.irewgt.ne.0) then
        do 68 idata=1,ndata
   68     dvpdv(idata)=one
      endif
c  If nucleii have been moved, move electrons
      if(istrch.ge.1) call strech_fit !JT

      write(6,'(''iopt='',i3)') iopt
      if(iopt.le.1) then

c       if(ncalls.gt.0) call zxssq2(func,ndata2,nparm,nsig,zero,zero,
c    &  ncalls,iopt,popt,parm,err2,diff,xjac,MDATA,xjtj,work,infer,ier)

c       write(6,'(/,''nsig, infer, ier='',3i5)') nsig,infer,ier
c       write(6,'(''if infer=1 then the param. estim. agree on 2 '',
c    &  ''succesive iter. to'',i3,'' sig. digits'',/,
c    &  ''if infer=2 the residual sum of squares estimates have '',
c    &  ''rel diff lt'',d12.4,/,''if infer=3 the norm of the approx. '',
c    &  ''gradient is lt'',d12.4)') nsig,eps1,del
c       write(6,'(/,''norm of grad, marquardt scaling param. ='',
c    &  2d12.4)') work(1),work(4)
c       write(6,'(''no of func eval., no of iter, no of sig. digits in
c    &  param. ='',2f6.0,f5.1)') work(2),work(5),work(3)

       else

        epsp=0.01d0*10.d0**(-nsig)
        epsg=epsp
        epsch2=epsp
        ichange=1
        rot_wt=1.d0
c       parmarmin=1.d-16
        eps_diff=1.d-15
        if(ianalyt_lap.eq.0) eps_diff=1.d-8
        if(analytic) then
          nanalytic=nparmcsf+nparmj+nparms
          write(6,'(''Number of parameters optimized analytically='',i4)') nanalytic
         else
          nanalytic=0
        endif
        call systemflush(6)
c       call quench(func,jacobian,analytic,parm,pmarquardt,tau,noutput,
c    &  nstep,ndata2,nparm,ipr_opt,diff,err2,epsg,epsp,converg,mesg,
c    &  ibold,cholesky,rot_wt,parmarmin,eps_diff)
c       call quench(func,jacobian,analytic,parm,pmarquardt,tau,noutput,
c    &  nstep,ndata2,nparm,ipr_opt,diff,err2,epsg,epsp,epsch2,converg,mesg,
c    &  ibold,cholesky,rot_wt,eps_diff)

#ifndef NOQUENCH
        call quench(func,jacobian,nanalytic,parm,pmarquardt,tau,noutput,
     &  nstep,ndata2,nparm,ipr_opt,diff,err2,epsg,epsp,epsch2,converg,mesg,
     &  ibold,cholesky,rot_wt,eps_diff)
#else
        stop 'needs quench library'
#endif

        if(converg) then
          write(6,'(''chisq='',d12.6,i5,'' func evals, convergence: '',
     &    a10)') err2,icalls,mesg
         else
          write(6,'(''chisq='',d12.6,i5,'' func evals, no convergence''
     &    )') err2,icalls
        endif

      endif

c The foll. call to func is no longer needed since quench takes care of it.
c Warning: I am almost but not completely sure about this. It may be needed for dependent parms.
c     dum_func=func(ndata2,nparm,parm,diff,0)
c     dum_func=func(ndata2,nparm,parm,diff,1)

      if(icusp2.ge.1 .and. isc.le.7) then
        do 95 isp=nspin1,nspin2
  95      if(ijas.eq.2) call cuspexact2(a1(1,isp,1),a2(1,isp,1))
        if(ijas.eq.3) call cuspexact3(1)
        if(ijas.ge.4.and.ijas.le.6) call cuspexact4(1,1)
      endif

      if(iabs(icusp2).ge.2) then
        if(ijas.eq.2) then
          do 100 isp=nspin1,nspin2
            ishft=ncuspc*(isp-nspin1)
            if(nspin2.eq.2 .and. (nup-ndn).ne.0) a1(2,2,1)=a1(2,1,1)
            if(nspin2.eq.3 .and. isp.eq.nspin2)
     &      a1(2,3,1)=((ndn-nup)*a1(2,1,1)+(nup-1)*a1(2,2,1))/(ndn-1)
  100       call cuspcheck2(scalek(1),a1(1,isp,1),a2(1,isp,1),
     &      diff(ndata+ishft+1),isp,nspin1,ncuspc,1)
         elseif(ijas.eq.3) then
          call cuspcheck3(diff(ndata+1),1)
        endif
      endif
      ishft=ncuspc*(nspin2-nspin1+1)+nfockc
      if((nloc.eq.0. .or. nloc.eq.6) .and. numr.le.0) call cuspco(diff(ndata+ishft+1),1)

      do 120 i=1,necn
  120   coef(iebasi(1,i),ieorb(1,i),1)=sign(one,dfloat(ieorb(2,i)))*
     &  coef(iebasi(2,i),iabs(ieorb(2,i)),1)
      do 130 i=1,nebase
  130   zex(iebase(1,i),1)=zex(iebase(2,i),1)
c     do 140 i=1,iabs(nedet)
c       cdet(iedet(1,i),1)=zero
c       do 140 j=2,icsf(i)
c 140     cdet(iedet(1,i),1)=cdet(iedet(1,i),1)+frac(j,i)*
c    &    sign(one,dfloat(iedet(j,i)))*cdet(iabs(iedet(j,i)),1)

c If the exponents of the basis functions are optimized, recompute normalizations
       if(nparme.gt.0) then
         if(ibasis.eq.3) then
           call basis_norm_dot(0)
         else
           call basis_norm(1,0)
         endif
       endif


c     fourpi=four*pi
      do 150 ib=1,nbasis
        do 150 iorb=1,norb
  150     coef(ib,iorb,1)=coef(ib,iorb,1)/anorm(ib)

c Pivot among orbitals of the same symmetry.
c Orbitals used to pivot other orbitals must appear
c in all the determinants (core states).

      if(nparml.gt.0.or.nparme.gt.0) then
c Pivot to make some coefs. zero
        do 200 iorb=1,norb
          indexi=ipivot(iorb)
          if(indexi.ne.0) then
            comax=zero
            do 170 ib=1,nbasis
              if(dabs(coef(ib,iorb,1)).gt.comax) then
                comax=dabs(coef(ib,iorb,1))
                jb=ib
              endif
  170       continue
            do 190 jorb=1,norb
             indexj=ipivot(jorb)
             ipij=iabs(indexi)-iabs(indexj)
             if((ipij.eq.0.and.indexi.gt.0).and.iorb.ne.jorb) then
               rat=coef(jb,jorb,1)/coef(jb,iorb,1)
               do 180 kbasis=1,nbasis
  180            coef(kbasis,jorb,1)=coef(kbasis,jorb,1)
     &                          -rat*coef(kbasis,iorb,1)
             endif
  190       continue
          endif
  200   continue

c Make the largest in absolute magnitude coef be 1.  This renormalizes the
c determinant, so adjust cdet (in old code) or csf_coef (in new code) to compensate.
        do 230 iorb=1,norb
          comax=zero
          do 210 ibas=1,nbasis
            if(dabs(coef(ibas,iorb,1)).gt.dabs(comax)) then
              comax=coef(ibas,iorb,1)
            endif
  210   continue
c         do 220 idet=1,ndet
c           do 220 j=1,nelec
c 220         if(iworbd(j,idet).eq.iorb) cdet(idet,1)=cdet(idet,1)*comax
c Warning: I am not sure if the foll. is always correct.  It assumes that if there
c are several determinants that make up a CSF, then each determinant has the
c same orbitals (though of course divided differently among up and down determinants).
c For example if we consider for the N atom the CSF that comes from exciting
c (a p->d for up-spin and a s->p for dn-spin) then the determinants that make up the
c CSF do not all have the same orbitals.  However, they have the same number of p and d
c orbitals and so if the same renormalization is applied to all the p's and all the d's
c then it is correct I think.
          do 220 icsf=1,ncsf
            idet=iwdet_in_csf(1,icsf)
            do 220 j=1,nelec
  220         if(iworbd(j,idet).eq.iorb) csf_coef(icsf,1)=csf_coef(icsf,1)*comax
            do 230 ibas=1,nbasis
  230         coef(ibas,iorb,1)=coef(ibas,iorb,1)/comax
      endif

c     if((nparml.gt.0.or.nparme.gt.0.or.nparmd.gt.0).and.cdet(1,1).ne.0.d0) then
c       term=one/cdet(1,1)
c       do 240 idet=1,ndet
c 240     cdet(idet,1)=cdet(idet,1)*term
      if((nparml.gt.0.or.nparme.gt.0.or.nparmcsf.gt.0).and.csf_coef(1,1).ne.0.d0) then
        term=one/csf_coef(1,1)
        do 240 icsf=1,ncsf
  240     csf_coef(icsf,1)=csf_coef(icsf,1)*term
      endif

      if(ijas.ge.4) then
        call set_scale_dist(0,1)
        call plot_jas
      endif

      write(6,'(/,''Final parameters:'')')

      if(iperiodic.eq.0) then
        if(nbasis.gt.999) stop 'nbasis > 999 in fit: increase i3 below'
        do 250 iorb=1,norb
          if(iorb.eq.1) then
c           write(fmt,'(''(''i3,''f14.8,\'\' ((coef(j,i),j=1,nbasis),i=1,norb)\'\')'')') nbasis
            write(fmt,'(''(''i3,''f14.8,a38)'')') nbasis
            write(6,fmt) (coef(ib,iorb,1),ib=1,nbasis),' ((coef(j,i),j=1,nbasis),i=1,norb)'
           else
c           write(fmt,'(a1,i3,a6)') '(',nbasis,'f14.8)'
            write(fmt,'(''(''i3,''f14.8)'')') nbasis
            write(6,fmt) (coef(ib,iorb,1),ib=1,nbasis)
          endif
  250   continue

c       write(fmt,'(''(''i3,''f14.8,\'\' (zex(i),i=1,nbasis)\'\')'')') nbasis
c       write(6,fmt) (zex(ib,1),ib=1,nbasis)
        write(fmt,'(''(''i3,''f14.8,a20)'')') nbasis
        write(6,fmt) (zex(ib,1),ib=1,nbasis),' (zex(i),i=1,nbasis)'

c       if(ndet.lt.10) then
c         write(fmt,'(a1,i1,a6)') '(',ndet,'f12.8)'
c        elseif(ndet.lt.100) then
c         write(fmt,'(a1,i2,a6)') '(',ndet,'f12.8)'
c        elseif(ndet.lt.1000) then
c         write(fmt,'(a1,i3,a6)') '(',ndet,'f12.8)'
c       endif

        if(ndet.gt.999) stop 'ndet > 999 in fit: increase i3 below'
c       write(6,*) 'Coef of dets.='
c       write(fmt,'(''(''i3,''f12.8,\'\' (cdet(i),i=1,ndet)\'\')'')') ndet
c       write(6,fmt) (cdet(idet,1),idet=1,ndet)
c       write(fmt,'(''(''i3,''f12.8,a19)'')') ndet
c       write(6,fmt) (cdet(idet,1),idet=1,ndet),' (cdet(i),i=1,ndet)'
        write(fmt,'(''(''i3,''f12.8,a)'')') ncsf
        write(6,fmt) (csf_coef(icsf,1),icsf=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
      endif

c     write(6,'(/,''Final Jastrow parameters:'')')

c     write(6,'(''scalek,a21='',9f13.8)') scalek(1),a21
      write(6,'(2f13.8,'' scalek,a21'')') scalek(1),a21

      if(ijas.eq.1) write(6,'(''cjas2='',f12.8)') cjas2(1)

      if(ijas.eq.2) then

        do 260 isp=nspin1,nspin2
c 260     write(6,'(''Jastrow a='',7f16.8,/,(8f16.8))') (a1(i,isp,1),i=1,nparm_read)
          if(nparm_read.gt.0) then
c           write(fmt,'(''(''i2,''f16.8,\'\' (a(iparmj),iparmj=1,nparma)\'\')'')') nparm_read
            write(fmt,'(''(''i2,''f16.8,a27'')') nparm_read
           else
            write(fmt,'(''(a27)'')') nparm_read
          endif
  260     write(6,fmt) (a1(i,isp,1),i=1,nparm_read),' (a(iparmj),iparmj=1,nparma)'
        do 265 isp=nspin1,nspin2
          if(nparm_read.gt.0) then
            write(fmt,'(''(''i2,''f16.8,a27'')') nparm_read
           else
c           write(fmt,'(''(\'\' (a(iparmj),iparmj=1,nparma)\'\')'')')
            write(fmt,'(''(a28)'')')
          endif
  265     write(6,fmt) (a2(i,isp,1),i=1,nparm_read),' (b(iparmj),iparmj=1,nparmb)'

       elseif(ijas.ge.3.and.ijas.le.6) then
        if(ijas.eq.3) then
c         write(6,'(''Jastrow a='',7f16.8,/,(8f16.8))') (a(i,1),i=1,nparm_read)
          if(nparm_read.gt.0) then
c           write(fmt,'(''(''i2,''f16.8,\'\' (a(iparmj),iparmj=1,nparma)\'\')'')') nparm_read
            write(fmt,'(''(''i2,''f16.8,a28)'')') nparm_read
           else
c           write(fmt,'(''(\'\' (a(iparmj),iparmj=1,nparma)\'\')'')')
            write(fmt,'(''(a28)'')')
          endif
          write(6,fmt) (a(i,1),i=1,nparm_read),' (a(iparmj),iparmj=1,nparma)'
          nparmc_read=(nord**3+5*nord)/6+nord**2+nord
         elseif(ijas.ge.4.and.ijas.le.6) then
          nparma_read=2+max(0,norda-1)
          do 267 it=1,nctype
c 267       write(6,'(''Jastrow a='',7f16.8,/,(8f16.8))') (a4(i,it,1),i=1,nparma_read)
            if(nparma_read.gt.0) then
c             write(fmt,'(''(''i2,''f16.8,\'\' (a(iparmj),iparmj=1,nparma)\'\')'')') nparma_read
              write(fmt,'(''(''i2,''f16.8,a28)'')') nparma_read
             else
c             write(fmt,'(''(\'\' (a(iparmj),iparmj=1,nparma)\'\')'')')
              write(fmt,'(''(a28)'')')
            endif
  267       write(6,fmt) (a4(i,it,1),i=1,nparma_read),' (a(iparmj),iparmj=1,nparma)'
          nparm_read=2+max(0,nordb-1)
          nparmc_read=nterms4(nordc)
        endif

        do 270 isp=nspin1,nspin2b
c 270     write(6,'(''Jastrow b='',7f16.8,/,(8f16.8))') (b(i,isp,1),i=1,nparm_read)
          if(nparm_read.gt.0) then
c           write(fmt,'(''(''i2,''f16.8,\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')') nparm_read
            write(fmt,'(''(''i2,''f16.8,a28)'')') nparm_read
           else
c           write(fmt,'(''(\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')')
            write(fmt,'(''(a28)'')')
          endif
  270     write(6,fmt) (b(i,isp,1),i=1,nparm_read),' (b(iparmj),iparmj=1,nparmb)'

        do 280 it=1,nctype
c 280     write(6,'(''Jastrow c='',7f16.8,/,(8f16.8))') (c(i,it,1),i=1,nparmc_read)
          if(nparmc_read.gt.0) then
            write(fmt,'(''(''i2,''f16.8,a28)'')') nparmc_read
           else
            write(fmt,'(''(a28)'')')
          endif
  280     write(6,fmt) (c(i,it,1),i=1,nparmc_read),' (c(iparmj),iparmj=1,nparmc)'

        if(ifock.gt.0) then
          do 285 it=1,nctype
c           write(6,'(''Jastrow f='',7f16.8,/,(8f16.8))')(fck(i,it,1),i=1,15)
  285       write(6,'(15f16.8,'' (f(iparmj),iparmj=1,nparmf)'')') (fck(i,it,1),i=1,15)
        endif
      endif

      if(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17) write(6,'(2f12.8,'' cutjas_en,cutjas_ee'')') cutjas_en,cutjas_ee

      write(6,'(''cusp penalty weighted by'',f10.5)') cuspwt

      if(ncalls.gt.0 .or. iopt.gt.1) errc=dsqrt(err2/ndata2)

      if((iopt.le.1 .and. ncalls.le.0) .or. ndata.ne.ndata2) then
        err2=zero
        do 290 i=1,ndata
  290     err2=err2+diff(i)**2
      endif
      err=dsqrt(err2/ndata)

      if(iopt.le.1 .and. ncalls.le.0 .and. ndata.ne.ndata2) then
        errc=err2
        nx=ndata+ncuspc*(nspin2-nspin1+1)+nfockc+ncent*norbc
        do 300 i=ndata+1,nx
  300     errc=errc+(diff(i)*cuspwt)**2
        do 310 i=nx+1,ndata2
  310     errc=errc+diff(i)**2
        errc=dsqrt(errc/ndata2)
      endif

      if(iabs(icusp2).ge.2) then
        if(mod(icusp2,3).eq.0) then
          write(6,'(''rms fit error,errc/'',i5,2d12.4,'' ndata,ndata2''
     &    ,2i5)') iabs(icusp2)-2,err,errc,ndata,ndata2
         else
          write(6,'(''rms fit error,errc*'',i5,2d12.4,'' ndata,ndata2''
     &    ,2i5)') nint(cuspwt),err,errc,ndata,ndata2
        endif
       else
        write(6,'(''rms fit error,errc='',2d12.4,'' ndata,ndata2=''
     &  ,2i5)') err,errc,ndata,ndata2
      endif

      if(ipr_opt.ge.2) write(6,'(/,''config #,   etrial,   efit,      diff
     &  relerr,   wght,     rmin,  rmax,  r12min'')')
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
      do 320 i=1,ndata
        yfit=eguess+uwdiff(i)
        rmin=1.d99
        rmax=0
        r12min=1.d99
        do 318 j=1,nelec
          r1=0
          do 312 k=1,ndim
  312       r1=r1+x(k,j,i)**2
          rmin=dmin1(r1,rmin)
          rmax=dmax1(r1,rmax)
          do 318 jj=1,j-1
            r12=0
            do 316 k=1,ndim
  316       r12=r12+(x(k,j,i)-x(k,jj,i))**2
  318       r12min=dmin1(r12,r12min)
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
c       if(ipr_opt.ge.0) write(6,'(i5,2f10.5,d13.5,2f9.3,2f12.6,f9.3)')
        if(ipr_opt.ge.2) write(6,'(i5,2f15.5,d13.5,2f9.3,2f12.6,f9.3)')
     &  i,eguess,yfit, uwdiff(i),yfit-eold(i),wght(i),rmin,rmax,r12min
         if(ipr_opt.ge.2) then
           if(ndim*nelec.lt.100) then
             write(fmt,'(a1,i2,a18)') '(',ndim*nelec,'f8.4,3d14.6,f12.5)'
            elseif(ndim*nelec.lt.1000) then
             write(fmt,'(a1,i3,a18)') '(',ndim*nelec,'f8.4,3d14.6,f12.5)'
           endif
           write(2,fmt)((x(k,j,i),k=1,ndim),j=1,nelec),psid(i),psij(i),
     &     psid(i)*exp(psij(i)),yfit
         endif
  320 continue

      if(iabs(icusp2).ge.2 .and. ipr_opt.ge.0) then
        do 330 i=ndata+1,ndata2
  330     write(6,'(i5,20x,d13.5)') i,diff(i)
      endif

      rminav=rminav/ndata
      rmaxav=rmaxav/ndata
      r12minav=r12minav/ndata
      eav=eav/ndata
      reldif=reldif/ndata
      relerr=dsqrt(max(0.d0,relerr/ndata-reldif*reldif))
      wtone=dsqrt(wtone/ndata)
      if(ovlbot.ne.zero) ovlap=ovltop/dsqrt(ndata*ovlbot)
      write(6,'(''ndata,rms fit error,eguess,eav,relerr,wtone,ovlap,rminav,rmaxav,r12minav'',
     &i5,f10.6,2f11.6,f10.6,2x,f6.3,f7.4,2x,9f6.3)')
     &ndata,err,eguess,eav,relerr, wtone,ovlap, rminav,rmaxav,r12minav
      if(analytic.and.irewgt.ne.0)
     &write(6,'(''*** Warning: analytic derivs not correctly implemented when points are reweighted'')')
      if(eguess.gt.eav+0.1*err)
     &write(6,'(''*** Warning: eguess is greater than average energy on fit pts + 0.1*std. Decrease eguess'')')
      if(eguess.lt.eav-0.3*err)
     &write(6,'(''*** Warning: eguess is less than average energy on fit pts - 0.3*std. Increase eguess'')')

      if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt
     &.gt.0) then
        ishft=ndata+ncuspc*(nspin2-nspin1+1)+nfockc+ncent*norbc+1
        do 340 isp=nspin1,nspin2
          call checkjas2(scalek(1),isp,ncnstr,diff(ishft),ipr_opt,1)
  340     ishft=ishft+ncnstr
      endif

      call my_second(2,'all   ')

      return

  999 write(6,'(''Error reading mc_configs, config.'',i5)') k
      stop 'Error reading mc_configs, possibly not enough configs.'

      end

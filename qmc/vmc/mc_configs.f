!=============================================================================
      subroutine mc_configs_read !JT
!-----------------------------------------------------------------------------
!=============================================================================

      implicit none

      character*16 mode
      common /contr3/ mode


      if(index(mode,'mpi').ne.0) then
         call mc_configs_read_mpi
      else
         call mc_configs_read_notmpi
      endif

      return !JT
c-----------------------------------------------------------------------
      entry mc_configs_write !JT

      if(index(mode,'mpi').ne.0) then
         call mc_configs_write_mpi
      else
         call mc_configs_write_notmpi
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine mc_configs_read_notmpi !JT
c If isite<=0 reads initial MC configuration from mc_configs_start
c         >=1 gets initial MC configuration by calling subroutine sites
c Write mc_configs_new at end of run to provide configurations for fit optimization

      use all_tools_mod
      use walkers_mod
      use atom_mod
      use config_mod
      use contrl_mod
      use const_mod
      implicit real*8(a-h,o-z)

      common /dim/ ndim

!JT      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
!JT      common /config/ xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)
!JT     &,vnew(3,MELEC),psi2o(MFORCE),psi2n(MFORCE),eold(MFORCE),enew(MFORCE)
!JT     &,peo,pen,peio,pein,tjfn,tjfo,psido,psijo
!JT     &,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)
!JT     &,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)
!JT     &,nearesto(MELEC),nearestn(MELEC),delttn(MELEC)
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /pairden/ xx0probut(0:NAX,-NAX:NAX,-NAX:NAX),xx0probuu(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probud(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdt(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probdu(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdd(0:NAX,-NAX:NAX,-NAX:NAX),
     &den2d_t(-NAX:NAX,-NAX:NAX),den2d_d(-NAX:NAX,-NAX:NAX),den2d_u(-NAX:NAX,-NAX:NAX),
     &delxi,xmax,xfix(3),ifixe

      dimension nsite(MCENT)

c Truncate fbias so that it is never negative, and the quantity
c sampled is never negative
c     fbias=dmin1(two,dmax1(zero,fbias))

      call alloc ('xold', xold, 3, nelec)

      if(irstar.ne.1) then

c if isite<=0 then get initial configuration from mc_configs_start
c if isite>=1 then get initial configuration from sites routine
        if(isite.ge.1) goto 393
        open(9,err=393,file='mc_configs_start')
        rewind 9
        read(9,*,end=393,err=393) ((xold(k,i),k=1,ndim),i=1,nelec)
        write(6,'(/,''initial configuration read from unit 9'')')
        close(9)
        goto 395

  393   l=0
        do 394 i=1,ncent
          nsite(i)=int(znuc(iwctype(i))+0.5d0)
          l=l+nsite(i)
          if(l.gt.nelec) then
            nsite(i)=nsite(i)-(l-nelec)
            l=nelec
          endif
  394     continue
        if(l.lt.nelec) nsite(1)=nsite(1)+(nelec-l)
        call sites(xold,nelec,nsite)
        call object_modified ('xold') ! JT
!JT        write(6,'(/,''initial configuration from sites'')')

  395   continue

c fix the position of electron i=ifixe for pair-density calculation:
        if(ifixe.gt.0) then
          do 398 k=1,ndim
  398       xold(k,ifixe)=xfix(k)
        endif


c If we are moving one electron at a time, then we need to initialize
c xnew, since only the first electron gets initialized in metrop
c       do 400 i=1,nelec
c         do 400 k=1,ndim
c 400       xnew(k,i)=xold(k,i)
      endif

c If nconf_new > 0 then we want to write nconf_new configurations from each processor for a future
c optimization or dmc calculation. So figure out how often we need to write a
c configuration to produce nconf_new configurations. If nconf_new = 0
c then set up so no configurations are written.
      if(nconf_new.ne.0) then
        open(7,form='formatted',file=trim(file_mc_configs_out))
        rewind 7
c       write(7,'(i5)') nconf_new
      endif
      return
c-----------------------------------------------------------------------
      entry mc_configs_write_notmpi !JT

c write out last configuration to unit mc_configs_start
      open(9,status='unknown',file='mc_configs_start')
      write(9,*) ((xold(k,i),k=1,ndim),i=1,nelec)
      close(9)

      return
      end

      subroutine mc_configs_read_mpi
c If isite<=0 reads initial MC configuration from mc_configs_start
c         >=1 gets initial MC configuration by calling subroutine sites
c Write mc_configs_new.<iproc> at end of run to provide configurations for fit optimization


# if defined (MPI)

      use walkers_mod
      use mpi_mod
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

c     dimension irn(4)
      dimension nsite(MCENT)
      dimension istatus(MPI_STATUS_SIZE)

      character*30 filename

c set the random number seed differently on each processor
c call to setrn is in read_input since irn is local there
c and in startr (entry in dumper.f) after reading in irn from unit 10.
c rnd itself is unused.
      if(irstar.ne.1) then
        do 95 id=1,(3*nelec)*idtask
   95     rnd=rannyu(0)

      call alloc ('xold', xold, 3, nelec)

c if isite=1 then get initial configuration from sites routine
        if(isite.ge.1) goto 393
          open(9,err=393,file='mc_configs_start')
          rewind 9
          do 392 id=0,idtask
  392       read(9,*,end=393,err=393) ((xold(k,i),k=1,ndim),i=1,nelec)
          write(6,'(/,''initial configuration from unit 9'')')
         goto 395
  393     l=0
          do 394 i=1,ncent
            nsite(i)=int(znuc(iwctype(i))+0.5d0)
            l=l+nsite(i)
            if(l.gt.nelec) then
              nsite(i)=nsite(i)-(l-nelec)
              l=nelec
            endif
  394       continue
          if(l.lt.nelec) nsite(1)=nsite(1)+(nelec-l)
          call sites(xold,nelec,nsite)
!JT          write(6,'(/,''initial configuration from sites'')')
  395 continue

c fix the position of electron i=ifixe for pair-density calculation:
      if(ifixe.gt.0) then
        do 398 k=1,ndim
  398     xold(k,ifixe)=xfix(k)
      endif

c If we are moving one electron at a time, then we need to initialize
c xnew, since only the first electron gets initialized in metrop
c       do 400 i=1,nelec
c         do 400 k=1,ndim
c 400       xnew(k,i)=xold(k,i)
      endif

c If nconf_new > 0 then we want to dump configurations for a future
c optimization or dmc calculation. So figure out how often we need to write a
c configuration to produce nconf_new configurations. If nconf_new = 0
c then set up so no configurations are written.
      if(nconf_new.ne.0) then
        write(6,'(''idtask,nconf_new'',9i5)') idtask,nconf_new
        if(idtask.lt.10) then
          write(filename,'(i1)') idtask
         elseif(idtask.lt.100) then
          write(filename,'(i2)') idtask
         elseif(idtask.lt.1000) then
          write(filename,'(i3)') idtask
         else
          write(filename,'(i4)') idtask
        endif
c The next line put in by Julien fails to put the numeric subscript on the filename.
c So temporarily go back to the hard-coded name we had before.
c       filename=file_mc_configs_out//filename(1:index(filename,' ')-1)
        filename='mc_configs_new'//filename(1:index(filename,' ')-1)
        open(7,form='formatted',file=filename)
!CU     open(7,form='formatted',file=file_mc_configs_out)
        rewind 7
c       write(7,'(i5)') nconf_new
      endif

# endif
      return
c-----------------------------------------------------------------------
      entry mc_configs_write_mpi

# if defined (MPI)

c write out last configuration to unit mc_configs_start
      if(idtask.ne.0) then
        call mpi_isend(xold,3*nelec,mpi_double_precision,0
     &  ,1,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
       else
        close(9)
        open(9,status='unknown',file='mc_configs_start')
        write(9,*) ((xold(k,i),k=1,ndim),i=1,nelec)
        do 450 id=1,nproc-1
          call mpi_recv(xnew,3*nelec,mpi_double_precision,id
     &    ,1,MPI_COMM_WORLD,istatus,ierr)
  450     write(9,*) ((xnew(k,i),k=1,ndim),i=1,nelec)
        close(9)
      endif

# endif

      return
      end

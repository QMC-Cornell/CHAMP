      subroutine mc_configs_read_mpi
! If isite<=0 reads initial MC configuration from mc_configs_start
!         >=1 gets initial MC configuration by calling subroutine sites
! Write mc_configs_new.<iproc> at end of run to provide configurations for fit optimization


# if defined (MPI)

      use walkers_mod
      use mpi_mod
      use atom_mod
      use config_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use pairden_mod
      implicit real*8(a-h,o-z)

      dimension nsite(ncent)
      dimension istatus(MPI_STATUS_SIZE)

      character*30 filename

! Set the random number seed differently on each processor.
! Seed is set for serial run by calling setrn in read_input and is reset for all except the first process
! by calling setrn again in vmc/MPI/mc_configs_read_mpi for vmc mpi runs and in
! dmc/dmc_elec/MPI/open_files_mpi for dmc mpi runs.
! It is also set in startr (entry in dumper.f) after reading in irand_seed from unit 10.
! rnd itself is unused.
! Frank's temporary fix for a better choice of random number seeds for parallel run.
      if(irstar.ne.1) then
        do 95 id=1,ndim*nelec*idtask
   95     rnd=rannyu(0)
! The next call to savern is not really needed but since Claudia put it in, I am too for consistency
        call savern(irand_seed)
        do i =1,4
          irand_seed(i)=mod(irand_seed(i)+int(rannyu(0)*idtask*9999),9999)
        enddo
        call setrn(irand_seed)
        write(6,'(''irand_seed='',4i5)') irand_seed

        call alloc ('xold', xold, 3, nelec)

! if isite=1 then get initial configuration from sites routine
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

! fix the position of electron i=ifixe for pair-density calculation:
      if(ifixe.gt.0) then
        do 398 k=1,ndim
  398     xold(k,ifixe)=xfix(k)
      endif

! If we are moving one electron at a time, then we need to initialize
! xnew, since only the first electron gets initialized in metrop
!       do 400 i=1,nelec
!         do 400 k=1,ndim
! 400       xnew(k,i)=xold(k,i)
      endif

! If nconf_new > 0 then we want to dump configurations for a future
! optimization or dmc calculation. So figure out how often we need to write a
! configuration to produce nconf_new configurations. If nconf_new = 0
! then set up so no configurations are written.
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
! The next line put in by Julien fails to put the numeric subscript on the filename.
! So temporarily go back to the hard-coded name we had before.
!       filename=file_mc_configs_out//filename(1:index(filename,' ')-1)
        filename='mc_configs_new'//filename(1:index(filename,' ')-1)
        open(7,form='formatted',file=filename)
!CU     open(7,form='formatted',file=file_mc_configs_out)
        rewind 7
!       write(7,'(i5)') nconf_new
      endif

# endif
      return
!-----------------------------------------------------------------------
      entry mc_configs_write_mpi

# if defined (MPI)

! write out last configuration to unit mc_configs_start
      if(idtask.ne.0) then
        call mpi_isend(xold,3*nelec,mpi_double_precision,0 &
     &  ,1,MPI_COMM_WORLD,irequest,ierr)
        call MPI_Wait(irequest,istatus,ierr)
       else
        close(9)
        open(9,status='unknown',file='mc_configs_start')
        write(9,*) ((xold(k,i),k=1,ndim),i=1,nelec)
        do 450 id=1,nproc-1
          call mpi_recv(xnew,3*nelec,mpi_double_precision,id &
     &    ,1,MPI_COMM_WORLD,istatus,ierr)
  450     write(9,*) ((xnew(k,i),k=1,ndim),i=1,nelec)
        close(9)
      endif

# endif

      return
      end

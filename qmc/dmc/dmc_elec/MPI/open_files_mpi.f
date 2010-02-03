      subroutine open_files_mpi
c Written by Cyrus Umrigar
c Open files that are different for serial and parallel runs
c Also do other things that differ for serial and parallel runs, such as setting random number seeds

# if defined (MPI)
      use all_tools_mod
      use dmc_mod
      use walkers_mod
      use dets_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use contr3_mod, only : mode
      use periodic_mod
      use config_dmc_mod
      use contrldmc_mod
      implicit real*8(a-h,o-z)

      character*27 fmt
      character*20 filename

c set the random number seed differently on each processor
c Seed is set for serial run by calling setrn in read_input and is reset by calling setrn again in
c vmc/MPI/mc_configs_read_mpi for vmc mpi runs and in
c dmc/dmc_elec/MPI/open_files_mpi for dmc mpi runs.
c It is also set in startr (entry in dumper.f) after reading in irand_seed from unit 10.
c rnd itself is unused.
c Frank's temporary fix for a better choice of random number seeds for parallelel run.
      if(irstar.ne.1) then
        do 95 id=1,ndim*nelec*idtask
   95     rnd=rannyu(0)
        do i =1,4
          irand_seed(i)=mod(int(irand_seed(i)+rannyu(0)*nelec*nstep*nblk*idtask),7777)
        enddo
        call setrn(irand_seed)
        write(6,'(''irand_seed='',4i5)') irand_seed
      endif

      if(ipr.gt.-2 .and. irstar.eq.0) then
         if (mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') then
            if (idtask == 0) then 
               filename='walkalize'
               open(11,file=filename)
               rewind 11
               write(11,'(i3,'' nblkeq to be added to nblock at file end'')') nblkeq
            end if
         else
            if(idtask.le.9) then
               write(filename,'(''walkalize.'',i1)') idtask
            elseif(idtask.le.99) then
               write(filename,'(''walkalize.'',i2)') idtask
            elseif(idtask.le.999) then
               write(filename,'(''walkalize.'',i3)') idtask
            elseif(idtask.le.9999) then
               write(filename,'(''walkalize.'',i4)') idtask
            elseif(idtask.le.99999) then
               write(filename,'(''walkalize.'',i5)') idtask
            elseif(idtask.le.999999) then
               write(filename,'(''walkalize.'',i6)') idtask
            else
               stop 'idtask > 999999'
            endif
            open(11,file=filename)
            rewind 11
c     write(11,*) 'Move line nstep*(2*nblkeq+nblk)+1 here and delete this line'
            write(11,'(i3,'' nblkeq to be added to nblock at file end'')') nblkeq
         end if

      endif

      call get_initial_walkers

c initialize sums and averages and reset nconf_global if there is one global population on all processors
c I do not see why another call to zeres0_dmc is needed for a global population.  I am commenting it out so mpi1 and mpi2,mpi3 give same energies for 1 processor.
c     if(irstar.ne.1 .and. (mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3')) then
c       call my_second(1,'zeres0')
c       call zeres0_dmc
ccThis is now done in read_input
cc      if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') then
cc        nconf_global=nconf_global*nproc
cc      endif
c     endif

c If nconf_new > 0 then we want to dump configurations for a future
c optimization or dmc calculation. So figure out how often we need to write a
c configuration to produce nconf_new configurations. If nconf_new = 0
c then set up so no configurations are written.
      if(nconf_new.eq.0) then
c       ngfmc=2*nstep*nblk
       else
c       ngfmc=max(1,(nstep*nblk)*nconf/nconf_new)
        if(idtask.lt.10) then
          write(filename,'(i1)') idtask
         elseif(idtask.lt.100) then
          write(filename,'(i2)') idtask
         elseif(idtask.lt.1000) then
          write(filename,'(i3)') idtask
         elseif(idtask.lt.10000) then
          write(filename,'(i4)') idtask
         else
          stop 'idtask>=10000'
        endif
        filename='mc_configs_new'//filename(1:index(filename,' ')-1)
        open(7,form='formatted',file=filename)
        rewind 7
c       write(7,'(i5)') nconf_new
      endif

# endif
      return
      end

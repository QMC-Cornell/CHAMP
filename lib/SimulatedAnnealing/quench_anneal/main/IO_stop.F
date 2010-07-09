c $Log: IO_stop.F,v $
c Revision 1.1.1.1  2006/10/08 19:38:16  toulouse
c champ
c
c Revision 1.3  2001/09/28 20:22:24  nigh
c *** empty log message ***
c
c Revision 1.2  2001/09/28 20:04:47  nigh
c fixed mpi_abort: 2->3 arguments
c
cDate of last modification: Fri Feb  6 17:33:27 EST 1998
      subroutine stop_process(jproc,message)
c stop a process and print a message
c jproc	= process number
c message	= message to be printed file stop_file.<process number>
      include '../include/IO_stop.h'
      include 'mpif.h'
      character*(*) message
      call open1(IUNIT_ERR,STOP_FILE_NAME,ISFN_LENGTH,'append',6,
     &  jproc,MXPROCS)
      write(IUNIT_ERR,*) ' process ',jproc,' stopped'
      write(IUNIT_ERR,'((a),(a))') ' reason ',message
      close(IUNIT_ERR)
      call mpi_abort(MPI_COMM_WORLD,ierr_final,ierr)
      end
c-----------------------------------------------------------------------

      subroutine open1(iunit,file_name,name_length,status,length_status,
     &  jproc,idispose)
c open IO unit in sequential append mode
c iunit		= unit number
c file_name	= file name will be file_name.<process number>
c name_length	= number of characters in file_name
c status	= IO status in open
c length_status	= number of characters in status
c jproc		= process number
c idispose	= standard output will be disposed of if jproc > idispose
      include '../include/IO_stop.h'
      character*1 char_tmp1
      character*2 char_tmp2
      character*3 proc_ident
      character*6 unit_num
      character*(*) file_name
      character*(*) status
      character*20 scratch_name
      character*84 full_file_name
      logical DEBUG
      parameter (DEBUG=.false.)
      if(DEBUG)
     &   write(6,*) 'open1: iunit,file_name,name_length,status,jproc ',
     &  iunit,file_name(1:name_length),name_length,status,jproc
      if(jproc.gt.MXPROCS.or.jproc.lt.0) then
        write(6,*)   'open1: not 0<=jproc<=MXPROCS',jproc,MXPROCS
        call stop_process(0,'open1: not 0<=jproc<=MXPROCS')
      endif
      if ( name_length + 4 .gt. 84 ) then
        call stop_process(jproc,
     &    'open1: character variable full_file_name too short')
      endif
      if(jproc.gt.idispose) then
        close(iunit)
        open(unit=iunit,status='scratch')
        inquire(unit=iunit,name=scratch_name)
c The next line works only for UNIX: whatever is written on unit=iunit
c goes immediately down the bit drain. (Does not work at Sandia.)
        call system('ln -sf /dev/null '//scratch_name)
      else
        if(jproc.ge.100) then
          write(proc_ident,'(i3)') jproc
        elseif(jproc.ge.10) then
          write(char_tmp2 ,'(i2)') jproc
          proc_ident='0'//char_tmp2
        else
          write(char_tmp1 ,'(i1)') jproc
          proc_ident='00'//char_tmp1
        endif
        if(status(1:length_status).eq.'append') then
	  full_file_name=file_name(1:name_length)//'.'//proc_ident
          open (unit=iunit,
     &      file=full_file_name,
     &      status='unknown',err=10)
        else
	  full_file_name=file_name(1:name_length)//'.'//proc_ident
          open (unit=iunit,
     &      file=full_file_name,
     &      status=status(1:length_status),err=10)
        endif
30      read (iunit,"(a)",err=10,end=20)
        goto 30
20      continue
      endif
c     backspace(iunit)
      return

10    continue
      if(iunit.eq.IUNIT_ERR) then
        write(6,*) 'open1: error opening unit = ',iunit
c        call mpi_finalize(ierr_final)
c        stop '''emergency stop in open'''
        call mpi_abort(MPI_COMM_WORLD,ierr_final,ierr)
      else
        write(unit_num,'(i6)') iunit
        call stop_process(jproc,'open1: error opening unit '//unit_num)
      endif
      end
c-----------------------------------------------------------------------

      subroutine open_stio(jproc)
c jproc		= number of calling process

c open FORTRAN standard IO units IUNIT_STI (usually 5) and IUNIT_STO
c (usually 6). See IO_stop.h.
c the file names are taken from a file with the name FILE_NAMES,
c as defined in IO_stop.h
c unit 5 is opened in read only mode.
c unit 6 is given an extra tag identifying the process
c if output file name is - standard output is used

c jproc	= process number

c FILE_NAMES looks as follows:
c <name of input file> or - if standard input
c <name of output file> or - if standard output
c idispose	= dispose of standard output if jproc > idispose

      include '../include/IO_stop.h'
      character*80 file_name,input_file,output_file

#if defined (IRIX)
      open(unit=IUNIT_NAMES,file=FILE_NAMES,dispose='readonly',err=10)
#else
#if defined (AIX)
      open(unit=IUNIT_NAMES,file=FILE_NAMES,action='read',err=10)
#else
#if defined (CLOBBER)
      open(unit=IUNIT_NAMES,file=FILE_NAMES,status='old',err=10)
#else
      stop 'system is not IRIX or AIX and CLOBBER is not defined'
#endif
#endif
#endif
      read(IUNIT_NAMES,'(a80)') file_name

      input_file=file_name(1:index(file_name,' ')-1)
      if(input_file(1:name_length).ne.'-') then
#if defined (IRIX)
        open(unit=IUNIT_STI,file=input_file,dispose='readonly',err=20)
#else
#if defined (AIX)
        open(unit=IUNIT_STI,file=input_file,action='read',err=20)
#else
#if defined (CLOBBER)
        open(unit=IUNIT_STI,file=input_file,status='unknown',err=20)
#else
      stop 'system is not IRIX or AIX and CLOBBER is not defined'
#endif
#endif
#endif
      endif

      read(IUNIT_NAMES,'(a80)') file_name
      read(IUNIT_NAMES,*) idispose
      if(jproc.le.idispose) then
        name_length=index(file_name,' ')-1
        output_file=file_name(1:name_length)
      else
        name_length=9
        output_file='/dev/null'
      endif
c open tagged output files
      if(output_file(1:name_length).ne.'-') then
          call open1(IUNIT_STO,output_file,name_length,'new',3,jproc,
     &      idispose)
      else
        if(jproc.gt.idispose) open(IUNIT_STO,file='/dev/null')
      endif
      write(6,*) 'open_stio: output= ',output_file
      write(6,*) 'open_stio: input = ',input_file
      return

10    continue
      call stop_process(jproc,'open_stio: error opening unit=IUNIT_NAMES')
20    continue
      call stop_process(jproc,'open_stio: error opening unit=IUNIT_STI  ')

      end

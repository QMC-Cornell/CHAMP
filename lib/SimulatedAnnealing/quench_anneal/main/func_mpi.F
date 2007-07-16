c $Log: func_mpi.F,v $
c Revision 1.1.1.1  2006/10/08 19:38:17  toulouse
c champ
c
c Revision 1.11  2005/10/18 14:40:42  nigh
c removed redundant declarations
c
c Revision 1.10  2004/02/06 14:11:20  nigh
c *** empty log message ***
c
c Revision 1.9  2002/12/20 19:00:32  nigh
c Cyrus's changes: some exact derivatives
c
c Revision 1.8  2002/03/04 14:21:32  nigh
c new version with LAPACK and LINPACK; see svd_gaus_test
c
c Revision 1.6  2002/01/18 18:39:17  nigh
c introduced two-dimensional diff in func_mpi and derivs
c
c
c Modification Mon Jan 25 1999:
c ndata is no longer required to be divisible by nproc.
c There is a new variable, more_data_qa; now ndata-more_data_qa should
c be divisible by nproc.  The first more_data_qa data are processes
c only by process 0.
      function func_mpi(func,ndata,nparm,parm,diff,iflag)
c wrapper to parallelize function calls
c func	= function to be called
c ndata         = number of data to be fitted
c nparm         = number of parameters
c diff		= residues of which func computes sum-of-squares
c iflag		= not used
      implicit real*8(a-h,o-z)
      external func
      common/quenchsim_diff/ndim1_diff,ndim2_diff
      dimension diff(ndata,ndim2_diff),parm(nparm)
#ifdef MPI
      include '../include/IO_stop.h'
      save time_g
      include 'mpif.h'
      include '../include/parameters.h'

      dimension diffl(ndata,ndim2_diff)
      dimension time_g(MXPROCS)
      dimension time_new_g(MXPROCS),ibegins(MXPROCS)
      dimension icnts(MXPROCS)
      do i=1,nproc
	icnts(i)=1
      enddo

      if(ndata.lt.nproc) then
        write(6,'((a),i4,(a))')
     &    'proc=',iproc,': func_mpi: ndata.lt.nproc'
        call stop(iproc,'func_mpi: ndata.lt.nproc')
      endif

c all processes get chunks that take equal amount of wall time
      call chunks_mpi(ibegins,time_g,ndata-more_data_qa)
      ibegin=ibegins(iproc+1)
      if(iproc.eq.nproc-1) then
        iend=ndata
      else
        iend=ibegins(iproc+2)-1
      endif

c each process computes a chunk of the diff's and stores the result in diffl
c the chunk has length ichunk and belongs to the data starting at
c data(ioffset+1)
      ichunk=iend-ibegin+1
      ioffset=ibegin-1
      do i=1,ndata
        diffl(i,1)=0
      enddo
      time_l=MPI_WTIME()
c first call to func
      ndim1_diff=ndata
      f=func(ichunk,nparm,parm,diffl,ioffset)
      time_l=MPI_WTIME()-time_l
c the next call seemed to give trouble on exeter; god only knows why.
c     call MPI_ALLGATHER(time_l,1,MPI_DOUBLE_PRECISION,time_new_g,
c    &  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
      call allgatherv(time_l,time_new_g,icnts)
c     do i=1,nproc
c       time_new_g(i)=time_l
c     enddo
      if(show_load_chunks) then
        write(6,'((a),(6(g10.3,x)))') 'func_mpi: times: ',
     &    (time_new_g(i),i=1,nproc)
        write(6,'((a),(6(i5,x)))') 'func_mpi: chunks: ',
     &    (ichunks(i),i=1,nproc)
      endif
      do i=1,nproc
        time_g_tot(i)=time_g_tot(i)+time_new_g(i)
        time_g(i)=THIRD*time_g(i)+
     &    TWO_THIRDS*time_new_g(i)/ichunks(i)
      enddo

c combine all chunks in diff
c     do k=1,ndim2_diff
c       call allgatherv(diffl(1,k),diff(1,k),ichunks)
c     enddo
      call allgatherv_2d1(diffl,diff,ndata,ndim2_diff,ichunks)
c final call to func
      ioffset=-1
      func_mpi=func(ndata,nparm,parm,diff,ioffset)

#else
      ioffset=0
      ndim1_diff=ndata
      f=func(ndata,nparm,parm,diff,ioffset)
      ioffset=-1
      func_mpi=func(ndata,nparm,parm,diff,ioffset)
#endif

      return
      end
c-----------------------------------------------------------------------

#ifdef MPI
      subroutine chunks_mpi(ibegins,time_g,ndata)
      implicit real*8(a-h,o-z)
c generate chunks of size inversely proportional to computation time
c ibegins       = starting points of the chunks
c ichunks       = size of the chunks (in common/chunks_mpi_c/)
c time_g        = wall clock time of the processes/unit work
c ndata         = number to be partitioned in chunks
c iatom         = chunks should be multiples of iatom
      include '../include/IO_stop.h'
      dimension ibegins(MXPROCS),time_g(MXPROCS)
      ndata_chunks=ndata/iatom
      if(first_call) then
        first_call=.false.
        do i=1,nproc
          time_g(i)=1d-30
        enddo
        if(iatom*ndata_chunks.ne.ndata) then
          write(6,'((a),i10,(a),i10)') 'chunks_mpi: ndata=',ndata,
     &      ' not divisible by iatom=',iatom
          call stop(iproc,'chunks_mpi: ndata not divisible by iatom')
        endif
        if(iatom*nproc.gt.ndata) then
          write(6,'((a),i3,1x,i3,i6)')
     &      'chunks_mpi: iatom*nproc > ndata',iatom,nproc,ndata
          call stop(iproc,'chunks_mpi: iatom*nproc > ndata')
        endif
      endif
      total=0
      do i=1,nproc
        if(static_chunks) time_g(i)=time_g(1)
        total=total+1/time_g(i)
      enddo
      mx_chunk=0
      ndata_chunks_rest=ndata_chunks
      rest=total
      do i=1,nproc
        ichunks(i)=min(ndata_chunks_rest-nproc+i,
     &    max(1,nint(ndata_chunks_rest/(time_g(i)*rest))))
        rest=rest-1/time_g(i)
        ndata_chunks_rest=ndata_chunks_rest-ichunks(i)
        if(ichunks(i).gt.mx_chunk) then
          mx_chunk=ichunks(i)
          mx_num=i
        endif
      enddo
      do i=1,nproc
        itot=itot+ichunks(i)
        ichunks(i)=ichunks(i)*iatom
      enddo
      ibegins(1)=1
      ichunks(1)=ichunks(1)+more_data_qa
      do i=2,nproc
        ibegins(i)=ibegins(i-1)+ichunks(i-1)
      enddo
      ibegins(1)=1
      return
      end
c-----------------------------------------------------------------------

      subroutine init_atom_mpi(jatom)
c jatom          = chunks should be multiples of iatom
      implicit real*8(a-h,o-z)
      include '../include/IO_stop.h'
      iatom=jatom
      first_call=.true.
      return
      end
c-----------------------------------------------------------------------

      subroutine init_more_data(more_data_qa1)
c more_data_qa1 = number of extra data to be processed only by process 0
      implicit real*8(a-h,o-z)
      include '../include/IO_stop.h'
      more_data_qa=more_data_qa1
      return
      end

c-----------------------------------------------------------------------

      subroutine allgatherv(from,to,icounts)
c share among all processes contiguously chunks of individual processes
c from(1:icounts(i)), i=0,...,nproc-1, concatenated in array "to"
c iproc = 0,...,nproc-1 is current process
c nproc = total number of processes
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      include '../include/IO_stop.h'
      dimension icounts(0:MXPROCS-1),idispls(0:MXPROCS-1)
      dimension from(1:*),to(1:*)
      if(nproc.gt.MXPROCS) then
        write(6,*) 'allgatherv: nproc>MXPROCS',nproc,MXPROCS
        call stop(iproc,'allgatherv: nproc>MXPROCS')
      endif
      idispls(0)=0
      do i=1,nproc-1
        idispls(i)=idispls(i-1)+icounts(i-1)
      enddo
      call mpi_allgatherv
     &  (from,icounts(iproc),MPI_DOUBLE_PRECISION,to,icounts,idispls
     &  ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
      return
      end
c-----------------------------------------------------------------------

#ifdef NON_STANDARD
      subroutine allgatherv_2d2(from,to,idim_phy1,idim_act1,icounts)
c before intelligent design created allgatherv_2d1, this was allgatherv_2d
c share among all processes contiguously with respect to SECOND array index
c chunks of individual processes
c concatenate from(1:idim_act1,1:icounts(i)), in "to" for increasing i
c for all nproc processes i=0,...,nproc-1
c from      = origin; the total number of elements of "from" may be no less than the actual number
c             of elements of "to"
c to        = destination
c idim_phy1 = physical dimension of first index of "from" and "to"
c idim_act1 = actual dimension of first index of "from"
c icounts   = icounts(iproc) actual dimension of second index of "from"
c             note: has to be defined for iproc=0,...,nproc-1
c From 'IO_stop.h' use:
c iproc     = number of current process
c nproc     = total number of processes
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      include '../include/IO_stop.h'
      dimension from(idim_phy1,1:*)
      dimension   to(idim_phy1,1:*)
      dimension icounts(0:MXPROCS-1),jcounts(0:MXPROCS-1),
     &  idispls(0:MXPROCS-1),jdispls(0:MXPROCS-1)
c dynamic equivalence (from,tmp2),(to,tmp1)
      pointer (itmp1_ptr,tmp1(1:*)),(itmp2_ptr,tmp2(1:*))
      itmp1_ptr=loc(to)
      itmp2_ptr=loc(from)

      if(nproc.gt.MXPROCS) then
        write(6,*) 'allgatherv_2d2: nproc>MXPROCS',nproc,MXPROCS
        call stop(iproc,'allgatherv_2d2: nproc>MXPROCS')
      endif

      do i=0,nproc-1
        jcounts(i)=icounts(i)*idim_act1
      enddo
      jdispls(0)=0
      icounts_tot=icounts(0)
      do i=1,nproc-1
        jdispls(i)=jdispls(i-1)+jcounts(i-1)
        icounts_tot=icounts_tot+icounts(i)
      enddo
      i3=0
      do i2=1,icounts(iproc)
        do i1=1,idim_act1
          i3=i3+1
          tmp1(i3)=from(i1,i2)
        enddo
      enddo

      call mpi_allgatherv
     &  (tmp1,jcounts(iproc),MPI_DOUBLE_PRECISION,tmp2,jcounts,jdispls
     &  ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)

      i3=0
      do i2=1,icounts_tot
        do i1=1,idim_act1
          i3=i3+1
          to(i1,i2)=tmp2(i3)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine allgatherv_2d1(from,to,idim_phy1,idim2,icounts)
c before intelligent design created this routine, allgatherv_2d2 was allgatherv_2d
c share among all processes contiguously with respect to FIRST array index
c chunks of individual processes
c concatenate from(1:icounts(i),1:idim2), in "to" for increasing i
c for all nproc processes i=0,...,nproc-1
c from      = origin; the total number of elements of "from" may be no less than the actual number
c             of elements of "to"
c to        = destination
c idim_phy1 = physical dimension of first index of "from" and "to"
c idim2     = second dimension of of "from"
c icounts   = icounts(iproc) actual dimension of second index of "from"
c             note: has to be defined for iproc=0,...,nproc-1
c From 'IO_stop.h' use:
c iproc     = number of current process
c nproc     = total number of processes
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      include '../include/IO_stop.h'
      dimension from(idim_phy1,1:*)
      dimension   to(idim_phy1,1:*)
      dimension icounts(0:MXPROCS-1),jcounts(0:MXPROCS-1),
     &  idispls(0:MXPROCS-1),jdispls(0:MXPROCS-1)
c dynamic equivalence (from,tmp2),(to,tmp1)
      pointer (itmp1_ptr,tmp1(1:*)),(itmp2_ptr,tmp2(1:*))
      itmp1_ptr=loc(to)
      itmp2_ptr=loc(from)

      if(nproc.gt.MXPROCS) then
        write(6,*) 'allgatherv_2d1: nproc>MXPROCS',nproc,MXPROCS
        call stop(iproc,'allgatherv_2d1: nproc>MXPROCS')
      endif

      do i=0,nproc-1
        jcounts(i)=icounts(i)*idim2
      enddo
      jdispls(0)=0
      icounts_tot=icounts(0)
      do i=1,nproc-1
        jdispls(i)=jdispls(i-1)+jcounts(i-1)
        icounts_tot=icounts_tot+icounts(i)
      enddo
      i3=0
      do i2=1,icounts(iproc)
        do i1=1,idim2
          i3=i3+1
          tmp1(i3)=from(i2,i1)
        enddo
      enddo

      call mpi_allgatherv
     &  (tmp1,jcounts(iproc),MPI_DOUBLE_PRECISION,tmp2,jcounts,jdispls
     &  ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)

      i3=0
      do i2=1,icounts_tot
        do i1=1,idim2
          i3=i3+1
          to(i2,i1)=tmp2(i3)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------

#ifdef TEST_ALLGATHERV_2D1
      program test_allgatherv_2d1
      implicit real*8(a-h,o-z)
      logical initialize,show_load_chunks
      parameter (IDIM1=100,IDIM2=4,MXPROCS=10)
      dimension x(IDIM1,IDIM2),y(IDIM1,IDIM2),icounts(MXPROCS)
      initialize=.true.
      show_load_chunks=.false.
      call start_mpi(initialize,iproc,nproc,show_load_chunks)
      do i=1,IDIM2
        do j=1,iproc+1
          x(j,i)=iproc+1+10*i+100*j
        enddo
      enddo
      itot=0
      do i=1,nproc
        icounts(i)=i
        itot=itot+icounts(i)
      enddo
      call allgatherv_2d1(x,y,IDIM1,IDIM2,icounts)
      do i=1,IDIM2
        print 1000,i,(y(j,i),j=1,itot)
 1000   format(i3,<itot>(1x,f4.0))
      enddo
      call stop_mpi(.true.)
      end
#endif

#endif
#endif
c-----------------------------------------------------------------------

      block data quenchsim_diff_data
c set default so that diff in function func is a one-dimensional array
      common/quenchsim_diff/ndim1_diff,ndim2_diff
      data ndim2_diff/1/
      end

      block data chunk_data
c set default so that chunks are static
      include '../include/IO_stop.h'
      common/chunk/static_chunks
      logical static_chunks
      data static_chunks/.true./,more_data_qa/0/
      end

      subroutine set_static_chunks(logical)
      include '../include/IO_stop.h'
      logical logical
      static_chunks=logical
      return
      end

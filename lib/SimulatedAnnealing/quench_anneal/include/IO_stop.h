c Date of last modification: Mon Dec 11 18:09:58 EST 2000

      integer IUNIT_ERR
      parameter( IUNIT_ERR = 99 )

      character*10 STOP_FILE_NAME
      integer ISFN_LENGTH
      parameter( STOP_FILE_NAME = 'stop_file0', ISFN_LENGTH = 10)

      integer MXPROCS
      parameter (MXPROCS=64)

      integer IUNIT_NAMES, IUNIT_STI, IUNIT_STO
      parameter (IUNIT_NAMES=10,IUNIT_STI=5,IUNIT_STO=6)

      character*11 FILE_NAMES
      parameter (FILE_NAMES='file_names0')

      logical first_call, show_load_chunks
      integer iatom, ichunks, iproc, nproc, more_data_qa
      double precision time_g_tot
      common / chunks_mpi_c / time_g_tot(MXPROCS), iatom, first_call,
     &  ichunks(MXPROCS), iproc, nproc, show_load_chunks, more_data_qa

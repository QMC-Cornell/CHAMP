c SimulatedAnnealing/quench_anneal/tests/test.f
      program test

      implicit real*8(a-h,o-z)

      include 'dimen.h'

      parameter( ZERO = 0.0d0 )

      character*10 mesg
      logical cholesky, converg
      integer nparm, ndata, nsig, noutput, ipr, maxstep,
     &     ibold, i, j, iseed, ntemp, nstep
      real*8 pold(MPARM), ptrue(MPARM), diff(MDATA),
     &     temp(MTEMP), epsp, pmarquardt, epsg, epschi, tau,
     &     rot_wt, eps_diff, efo, chisq, epsch2, dseed

      real*8 func
      external func, jacobian

#ifdef MPI
      include 'mpif.h'
      logical initialize,show_load_chunks,finalize
      initialize=.true.
      show_load_chunks=.false.
      call start_mpi(initialize,iproc,nproc,show_load_chunks)
      call init_atom_mpi(1)
#endif

c     number of parameters
      nparm = 2
      if ( nparm .gt. MPARM ) then
         stop 'test: nparm > MPARM'
      endif

c     number of data points
      ndata = 6
      if ( ndata .gt. MDATA ) then
         stop 'test: ndata > MDATA'
      endif

c     number of different temperatures in annealing schedule
      ntemp = 2
      if ( ntemp .gt. MTEMP ) then
         stop 'test: ntemp > MTEMP'
      endif

c     We will be fitting a function to a set of ndata (x,y) data points
c     (only init_func and func need access to these data points, which
c     are stored in a common block invisible to this program and to
c     quench-anneal).  The input data points are in fact obtained by
c     sampling the fitting function at ndata different x values for a
c     particular set of parameters ptrue.  (There is also a facility to
c     add noise to the sampled values --- see subroutine init_func.)
c     The input parameters ptrue are chosen here, and init_func is then
c     called to generate the data.  After running quench-anneal, the
c     output parameters returned in pold should be close to those in
c     ptrue.
      do i = 1, nparm
         ptrue(i) = 0.0d0
      end do
      call init_func( ptrue, nparm, ndata )

      nsig = 12
      epsp = 10.0d0**(-nsig)
      pmarquardt = 1.0d-5
      epsg = epsp
      epschi = 0.001d0
      noutput = 2
      ipr = 1
      maxstep = 1000
      tau = 4.0d0
      ibold = 4
      rot_wt = 1.0d0
      eps_diff = ZERO
      cholesky = .true.
      nanalytic = 0

      write(6,*) ' test: quench1'
      do i = 1, nparm
         pold(i) = dble(i)
      end do
      call quench1( func, jacobian, nanalytic, pold, maxstep,
     &     ndata, nparm, diff, efo, epsg, epsp, converg,
     &     mesg, ibold )

      do i = 1, 2

         do j = 1, nparm
            pold(j) = dble(j)
         end do
         write(6,*)
         write(6,*) ' test: quench'

         call quench( func, jacobian, nanalytic, pold, pmarquardt,
     &        tau, noutput, maxstep, ndata, nparm, ipr, diff,
     &        chisq, epsg, epsp, epschi, converg, mesg, ibold,
     &        cholesky, rot_wt, eps_diff )

         if ( cholesky ) then
            write(6, '(''test:  LU-mode: minimum at: '')' )
            write(6, '(5f15.10)' ) (pold(j), j=1,nparm)
         else
            write(6, '(''test:  SVD-mode: minimum at: '')' )
            write(6, '(5f15.10)' ) (pold(j), j=1,nparm)
         endif

         cholesky = .not.cholesky

      enddo

      iseed = 32123
      call ransi( iseed )
      epsp = ZERO
      epsg = ZERO
      epsch2 = ZERO
      do i = 1, ntemp
         temp(i) = (10.0**(ntemp+1-i)) * chisq
      end do
      nstep = 10

c Run test in cholesky mode, then in svd mode
      do i = 1, 2

c        do j = 1, nparm
c           pold(j) = dble(j)
c        end do
         write(6,*)
         write(6,*) ' test: anneal'

         call anneal( func, jacobian, nanalytic, pold, pmarquardt,
     &        tau, noutput, nstep, ndata, nparm, ipr, diff,
     &        chisq, epsg, epsp, epsch2, converg, mesg, dseed,
     &        temp, ntemp, cholesky, rot_wt, eps_diff )

         if ( cholesky ) then
            write(6, '(''test:  LU-mode: minimum at: '')' )
            write(6, '(5f15.10)' ) (pold(j), j=1,nparm)
         else
            write(6, '(''test:  SVD-mode: minimum at: '')' )
            write(6, '(5f15.10)' ) (pold(j), j=1,nparm)
         endif

         cholesky = .not.cholesky

      enddo

#ifdef MPI
      finalize=.true.
      call stop_mpi(finalize)
#endif

      end
c-----------------------------------------------------------------------

      subroutine jacobian
c     jacobian: dummy routine to keep the loader happy
      stop 'jacobian should not have been called'
      end

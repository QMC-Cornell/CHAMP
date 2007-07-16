***********************************************************************

      subroutine get_job_name(job_name)

c Written by Devrim Guclu for Cyrus Umrigar
*----------------------------------------------------------------------
* description: reads the argument which is the job name which is used
*              to define the input file name and output
*              files
* arguments:   job_name = argv
*----------------------------------------------------------------------

      character*35      job_name
      character*35      argv

      write(*,*) '* Getting job argument'

      call getarg(1,argv)
      if (argv(1:1).ne.' ') then
         job_name = argv
      else
         write(*,*) ' Use: prog_name job_name '
         stop
      end if
      write(*,*) '  Job name = ',job_name
      write(*,*) '  OK.'
      write(*,*)

      return
      end

**********************************************************************

      subroutine read_input(job_name)

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: reads the input file/ check for errors
* arguments: job_name
* return values: initial values for nelec,ndet,det, and cdet are
*                read from the file "job_name.in"
*---------------------------------------------------------------------
      include         'maxdim.h'
      include         'qnumbers.h'

c arguments:
      character*35    job_name

c locals:
      character*35    inputfile
      integer         idet,in,iq,indet
      integer         iqns,iqnn,iqnm,iqnl
      real*8          amp

c commons:
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)

      common/psi_n/   cdet,ndim,nelec,ndet,det

      write(*,*) '* Reading the input file'

c define input file name
      ii=1
      do 10 while (job_name(ii:ii).ne.' ')
         ii=ii+1
 10   enddo
      lname=ii-1
      inputfile=job_name
      inputfile(lname+1:lname+3)='.in'
      write(*,*) '  Input file name = ',inputfile

c open the file
      open(33,file=inputfile,status='old')

c read the file
      read(33,*) ndim,nelec,ndet
      write(*,*) '  Number of spatial dimensions = ', ndim
      write(*,*) '  Number of electrons          = ', nelec
      write(*,*) '  Number of determinants       = ', ndet

      if (ndim.lt.2 .or. ndim.gt.3) then
         write(*,*) 'PROBLEM! ndim must be 2 or 3!'
         stop
      endif
      if (nelec.gt.MAXNELEC) then
         write(*,*) 'PROBLEM! nelec too big!'
         stop
      endif
      indet=nint(dfloat(MAXNDET)/2.d0**nelec)
      if (ndet.gt.indet) then
         write(*,*) 'PROBLEM! ndet too big!'
         stop
      endif

      do 20 idet=1,ndet
         read(33,*) amp
         write(*,*)
         write(*,*) '  Determinant: ',idet,', Coeff.: ',amp
         cdet(idet)=dcmplx(amp,0.d0)
         do 15 in=1,nelec
            read(33,*) (det(idet,in,iq),iq=qns,qnm)
            write(*,'(''     Electron'',i3,'' : | s = '',i2,
     &        '' , n = '',i2,'' , l = '',i2,'' , m = '',i2,'' )'')')
     &              in,det(idet,in,qns),det(idet,in,qnn),
     &                     det(idet,in,qnl),det(idet,in,qnm)

c check if quantum numbers are in the minimum/maximum limits:
            iqns=det(idet,in,qns)
            iqnn=det(idet,in,qnn)
            iqnl=det(idet,in,qnl)
            iqnm=det(idet,in,qnm)
            if (abs(iqns).ne.MAXS) then
               write(*,*) 'PROBLEM! s must be +/-1!'
               stop
            endif
            if (iqnn.lt.0.or.iqnn.gt.MAXN) then
               write(*,*) 'PROBLEM! n must be between 0,MAXN!'
               stop
            endif
            if (iqnl.lt.0.or.iqnl.gt.MAXL) then
c this limit for l must be respected for 2D systems also, because the
c ordering subtroutine def_ses does not know about the the
c dimensionality of the problem
               write(*,*) 'PROBLEM! qnl must be between 0,MAXL!'
               stop
            endif
            if (abs(iqnm).gt.MAXM) then
               write(*,*) 'PROBLEM! qnm must be between -MAXM,MAXM!'
               stop
            endif

c check if the quantum numbers satisfy "atomic orbital rules" in 3D.
            if (iqnl.ge.iqnn .and. ndim.eq.3) then
               write(*,*) 'PROBLEM! For 3D systems l must be smaller than n-1'
               stop
            endif
            if (abs(iqnm).gt.iqnl .and. ndim.eq.3) then
               write(*,*) 'PROBLEM! For 3D systems l must be higher than or equal to abs(m)'
               stop
            endif
 15      enddo
 20   enddo

      close(33)

      write(*,*) '  OK.'
      write(*,*)
      return
      end

**********************************************************************

      subroutine write_output(job_name)

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: defines output file name, writes the results.
* arguments: job_name=argv
* return values: none
*---------------------------------------------------------------------
      include         'maxdim.h'
      include         'qnumbers.h'

c arguments:
      character*35    job_name

c locals:
      real*8          norm2
      character*35    outputfile
      integer         idet,in,rcount,icount

c commons:
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)

      common/psi_n/   cdet,ndim,nelec,ndet,det

      write(*,*) '* Writing input file'

c define input file name
      ii=1
      do 10 while (job_name(ii:ii).ne.' ')
        ii=ii+1
 10   enddo
      lname=ii-1
      outputfile=job_name
      outputfile(lname+1:lname+4)='.out'
      write(*,*) '  output file name = ',outputfile

c open the file
      open(33,file=outputfile)

      norm2=0.d0
      rcount=0
      icount=0
      write(33,*) ' REAL PART:'
      write(33,*)
      do 20 idet=1,ndet
        if (dreal(cdet(idet)).ne.0.d0) then
          rcount=rcount+1
          write(33,*) ' Coeff.: ', dreal(cdet(idet))
          do 17 in=1,nelec
             write(33,'(''     Electron'',i3,'' : | s = '',i2,
     &        '' , n = '',i2,'' , l = '',i2,'' , m = '',i2,'' )'')')
     &              in,det(idet,in,qns),det(idet,in,qnn),
     &                     det(idet,in,qnl),det(idet,in,qnm)
 17       enddo
          norm2=norm2+dreal(cdet(idet))*dreal(cdet(idet))
        endif
 20   enddo
      if (rcount.eq.0) write(33,*) ' No real part'

      write(33,*)
      write(33,*) ' IMAGINARY PART:'
      write(33,*)
      do 30 idet=1,ndet
        if (dimag(cdet(idet)).ne.0.d0) then
          icount=icount+1
          write(33,*) ' Coeff.:', dimag(cdet(idet))
          do 25 in=1,nelec
            write(33,'(''     Electron'',i3,'' : | s = '',i2,
     &        '' , n = '',i2,'' , l = '',i2,'' , m = '',i2,'' )'')')
     &              in,det(idet,in,qns),det(idet,in,qnn),
     &                     det(idet,in,qnl),det(idet,in,qnm)
 25       enddo
          norm2=norm2+dimag(cdet(idet))*dimag(cdet(idet))
        endif
 30   enddo
      if (icount.eq.0) then
         write(33,*) ' No imaginary part'
         write(33,*)
         if (rcount.eq.0) then
            write(33,*) ' WARNING! Null wavefunction. Verify the initial',
     &                  ' wavefunction'
         endif
      endif
      write(33,*)
      write(33,*) ' NORMALIZATION Norm^2 : ', norm2


      close(33)

      write(*,*) '  OK.'

      return
      end

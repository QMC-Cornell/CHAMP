***********************************************************************

      subroutine get_job_name

c Written by Devrim Guclu for Cyrus Umrigar
*----------------------------------------------------------------------
* description: reads the argument which is the job name which is used
*              to define the input file name and output
*              files
* arguments:   job_name = argv
*----------------------------------------------------------------------

      character*35  job_name
      character*35  dstrfile,inputfile,outfile,s2file,detfile,csffile,conffile,dpfile
	common/bfield/  ibfield,dstrfile

      write(*,*) '* Getting job argument'


	read(*,*) job_name
      write(*,*) '  Job name = ',job_name
      write(*,*) '  OK.'

c define file names
      ii=1
      do 10 while (job_name(ii:ii).ne.' ')
         ii=ii+1
 10   enddo
      lname=ii-1
      inputfile=job_name
      inputfile(lname+1:lname+3)='.in'
      write(*,*) '  Input file name = ',inputfile

      outfile=job_name
      outfile(lname+1:lname+4)='.out'
      write(*,*) '  output file name = ',outfile

      s2file=job_name
      s2file(lname+1:lname+3)='.s2'
      write(*,*) '  S^2 eigenvector file name = ',s2file

      detfile=job_name
      detfile(lname+1:lname+4)='.det'
      write(*,*) '  determinant file name = ',detfile

      csffile=job_name
      csffile(lname+1:lname+4)='.csf'
      write(*,*) '  CSF file name = ',csffile

      conffile=job_name
      conffile(lname+1:lname+5)='.conf'
      write(*,*) '  configuration file name = ',conffile

      dstrfile=job_name
      dstrfile(lname+1:lname+5)='.dstr'
      write(*,*) '  distribution file name = ',dstrfile

c      dpfile=job_name
c      dpfile(lname+1:lname+3)='.dp'
c      write(*,*) '  linear dependance file name = ',dpfile

c open the file
      open(18,file=inputfile,status='UNKNOWN')
	open(19,file=conffile,status='UNKNOWN')
	open(20,file=s2file,status='UNKNOWN')

c	if (ibfield .eq. 1) open(11,file=dstrfile)

	open(15,file=outfile,status='UNKNOWN',access='sequential',form='formatted')
c	open(16,file=detfile,status='UNKNOWN',access='sequential',form='formatted')
c	open(17,file=csffile,status='UNKNOWN',access='sequential',form='formatted')
	open(16,file=detfile,status='REPLACE',access='sequential',form='formatted')
	open(17,file=csffile,status='REPLACE',access='sequential',form='formatted')

c	open(27,file=dpfile,status='REPLACE',access='sequential',form='formatted')

      return
      end

**********************************************************************

      subroutine read_input

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





c read the file
      read(18,*) ndim,nelec,ndet,icsf,iconf
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
         read(18,*) amp
         write(*,*)
         write(*,*) '  Determinant: ',idet,', Coeff.: ',amp
         cdet(idet)=dcmplx(amp,0.d0)
         do 15 in=1,nelec
            read(18,*) (det(idet,in,iq),iq=qns,qnm)
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


      write(*,*) '  OK.'
      write(*,*)
      return
      end

**********************************************************************

      subroutine write_output(mmax)

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
      real*8          norm2, coeff(MAXNDET), coeffi(MAXNDET)
      character*35    outputfile
      integer         idet,in,rcount,icount, mmax(0:MAXN), mplus(MAXNELEC)
	integer         morb(2*MAXM+1)


c commons:
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)

      common/psi_n/   cdet,ndim,nelec,ndet,det

      write(*,*) '* Writing input file'



c assume 1st Landau Level(LL) orbitals with m=0, 1, -1, 2, -2, ... morbmax, -morbmax
c in order, then 2nd LL orbitals with m=0, 1, -1, 2, -2, ..., then 3rd LL
c      write(*,*) ' input maximum of m for 1st Landau Level'
c	read(*,*)  morbmax
c      write(*,*) ' input maximum of m for 2nd Landau Level'
c	read(*,*)  morbmax1
c	if (morbmax .gt.MAXM) then
c         write(*,*) '  PROBLEM! cannot handle such a large m value'
c         stop
c      endif


      rcount=0
      icount=0
	write(15,'(//,i3,6x,"electrons",/)') nelec
      write(15,*) ' REAL PART:'
	norm2=0.0
      do idet=1,ndet
        if (abs(dreal(cdet(idet))).gt.1.d-6) rcount=rcount+1
	  norm2=norm2+dreal(cdet(idet))**2	
      enddo
      write(15,'(i3,6x,"determinantal coefficients")')  rcount
      rcount=0
      do 20 idet=1,ndet
        if (abs(dreal(cdet(idet))).gt.1.d-6) then
	    rcount=rcount+1
c	    if (rcount .eq. 1)	norm2=abs(dreal(cdet(idet)))
          coeff(rcount)=dreal(cdet(idet))/dsqrt(norm2)
	    do ielec=1,nelec
		 mplus(ielec)=0
	     do ii=0,det(idet,ielec,qnn)-1
	      mplus(ielec)=mplus(ielec)+2*mmax(ii)+1
	     enddo
	    enddo
          write(15,'(<nelec>i3)') (2*abs(det(idet,in,qnm))
     &	+(1-sign(1,det(idet,in,qnm)-1))/2+mplus(in),in=1,nelec)
        endif
20	enddo

      if (rcount.ne.0) write(15,'(<rcount>f12.6)') (coeff(i),i=1,rcount)

      write(15,*)
      write(15,*) ' IMAGINARY PART:'
	norm2=0.0
      do idet=1,ndet
        if (abs(dimag(cdet(idet))).gt.1.d-6) icount=icount+1
	  norm2=norm2+dimag(cdet(idet))**2	
      enddo
      write(15,'(i3,6x,"determinantal coefficients")')  icount
      icount=0
      do 30 idet=1,ndet
        if (abs(dimag(cdet(idet))).gt.1.d-6) then
          icount=icount+1
c	    if (icount .eq. 1)	norm2=abs(dimag(cdet(idet)))
          coeffi(icount)=dimag(cdet(idet))/dsqrt(norm2)
	    do ielec=1,nelec
		 mplus(ielec)=0
	     do ii=0,det(idet,ielec,qnn)-1
	      mplus(ielec)=mplus(ielec)+2*mmax(ii)+1
	     enddo
	    enddo
          write(15,'(<nelec>i3)') (2*abs(det(idet,in,qnm))
     &	+(1-sign(1,det(idet,in,qnm)-1))/2+mplus(in),in=1,nelec)
        endif
 30   enddo
      if (icount.ne.0) write(15,'(<icount>f12.6)') (coeffi(i),i=1,icount)
      if (icount.eq.0) then
         if (rcount.eq.0) then
            write(15,*) ' WARNING! Null wavefunction. Verify the initial',
     &                  ' wavefunction'
         endif
      endif
	



      write(*,*) '  OK.'

      return
      end

      subroutine write_detcsf
      implicit real*8(a-h,o-z)
	character*35    job_name
	logical eq_csf
	parameter(max=15)
	dimension idet(200,max),iidet(200,max),ioldet(9999,max),inewdet(200,max)
	dimension indet(200),coeff(200),coeffi(200)
	dimension icsfd(1000,200),csfc(1000,200),ndet_csf(1000)
	
	idcount=0
	ncsfsum=0
c read in determinants and their co-efficients	
	read(15,*) nelec
	read(15,*)
	read(15,*)
	read(15,*) nrdet
	if (nrdet .gt. 200) then
	 write(*,*) 'increase array dimension!'
	 stop
	endif
	do i=1,nrdet
	 read(15,*) (idet(i,j),j=1,nelec)
c	 write(*,*) (idet(j),j=1,nelec)
	enddo
	if (nrdet .ne. 0) read(15,*) (coeff(i),i=1,nrdet)
	read(15,*)
	read(15,*)
	read(15,*) nidet
	if (nidet .gt. 200) then
	 write(*,*) 'increase array dimension!'
	 stop
	endif
	do i=1,nidet
	 read(15,*) (iidet(i,j),j=1,nelec)
	enddo
	if (nidet .ne. 0) read(15,*) (coeffi(i),i=1,nidet)


c read in old determinants
	rewind(16)
	read(16,*,end=10) idcount
	do i=1,idcount
	 read(16,*) (ioldet(i,j),j=1,nelec)
	enddo
10	write(*,*) 'Begin writing', idcount
c	 stop

c	if (nrdet .eq. 0) goto 21
c compare real part with old determinants
	call detcomp(idet,nrdet,ioldet,idcount,inewdet,newdet,indet,nelec)

c append new determinants	

	if (newdet.ne.0) then
	do i=1,newdet
	 do j=1,nelec
	  ioldet(idcount+i,j)=inewdet(i,j)
	 enddo
	enddo
	idcount=idcount+newdet
	rewind(16)
	write(16,'(i5,3x,"ndet")') idcount
	do i=1,idcount
	 idig=1
	 if(ioldet(i,1).gt.9) idig=2
c	 if(ioldet(i,1).gt.99) idig=3
	 write(16, '(i<idig>,<nelec-1>i3)') (ioldet(i,j),j=1,nelec)
	enddo
	endif


c read in old CSFs. If there is repeat, use imaginary part instead.


	eq_csf=.false.
	rewind(17)
	read(17,*,end=21) ncsfsum
	read(17,*) (ndet_csf(i),i=1,ncsfsum)
	do i=1,ncsfsum
	 read(17,*) (icsfd(i,j),j=1,ndet_csf(i))
	 read(17,*) (csfc(i,j),j=1,ndet_csf(i))
	enddo
	if (newdet .ne. 0) goto 22
	do i=1,ncsfsum
	 eq_csf=.true.
	 if(nrdet .ne. ndet_csf(i)) then
	  eq_csf=.false.
	  goto 20
	 endif
	 do j=1,nrdet
	  if(indet(j) .ne. icsfd(i,j)) then
	   eq_csf=.false.
	   goto 20
	  endif
	  if(j .eq. 1) ratio=csfc(i,j)/coeff(j)
	  if(csfc(i,j)/coeff(j) .ne. ratio) then
	   eq_csf=.false.
c	   write(*,*) 'CSFs use same determinants, check for linear dependence'
c	   pause
	   goto 20
	  endif
	 enddo
20	 if (eq_csf) go to 21
	enddo
c	if (not(eq_csf)) go to 22

c use imaginary part as CSF, if real part is already used.
c Is it possibble both real part and imaginary part are used?
c No, unless the complex wavefunctions are linear dependent.
21	if (eq_csf .or. (nrdet .eq. 0)) then

	call detcomp(iidet,nidet,ioldet,idcount,inewdet,newdet,indet,nelec)

c append new determinants
	do i=1,newdet
	 do j=1,nelec
	  ioldet(idcount+i,j)=inewdet(i,j)
	 enddo
	enddo
	idcount=idcount+newdet
	rewind(16)
	write(16,'(i5,3x,"ndet")') idcount
	do i=1,idcount
	 idig=1
	 if(ioldet(i,1).gt.9) idig=2
	 write(16, '(i<idig>,<nelec-1>i3)') (ioldet(i,j),j=1,nelec)
	enddo
	
	nrdet=nidet
	do i=1,nidet
	 coeff(i)=coeffi(i)
	enddo
	if (nrdet .eq. 0) then
	 write(*,*) 'Something wrong, check the input files'
	 stop
	endif

	endif

c rewrite ncsfsum, ndet_csf(i) and CSFs

22	ncsfsum=ncsfsum+1
	ndet_csf(ncsfsum)=nrdet
	do i=1,nrdet
	 icsfd(ncsfsum,i)=indet(i)
	 csfc(ncsfsum,i)=coeff(i)
	enddo 	
	rewind(17)
	write (17,'(i4,3x,"ncsfsum")') ncsfsum
	write(17,'(<ncsfsum>i4,3x,"(ndet_in_csf(icsf),icsf=1,ncsfsum)")') (ndet_csf(i), i=1,ncsfsum)
	do i=1,ncsfsum
	 idig=1
	 if(icsfd(i,1).gt.9) idig=2
	 if(icsfd(i,1).gt.99) idig=3
	 if(icsfd(i,1).gt.999) idig=4
	 write(17,'(i<idig>,500i5)') (icsfd(i,j), j=1,ndet_csf(i))
	 idig1=8
	 if(csfc(i,1).lt.0.0) idig1=9
	 write(17,'(f<idig1>.6,500f12.6)') (csfc(i,j), j=1,ndet_csf(i))
	enddo

	return
	end

*******************************************************************************
	subroutine detcomp(idet,ndet,ioldet,idcount,inewdet,newdet,indet,nelec)
      implicit real*8(a-h,o-z)
	logical eq_det
	parameter(max=15)
	dimension idet(200,max),ioldet(9999,max),inewdet(200,max),indet(200)

c compare determinants with old ones, obtain the order numbers for repeated one,
c add previously unused ones and assign order numbers.
	newdet=0
	if (idcount .eq. 0) then
	 do i=1,ndet
	  newdet=newdet+1
	  indet(i)=newdet
	  do j=1,nelec
	   inewdet(i,j)=idet(i,j)
	  enddo
	 enddo
	else
	 do i=1,ndet
	  do k=1,idcount
	   eq_det=.true.
	   do j=1,nelec
	    if(idet(i,j) .ne. ioldet(k,j)) then
	     eq_det=.false.
	     goto 11
	    endif
	   enddo
11	   if(eq_det) goto 12
	  enddo
12	  if(eq_det) then
	   indet(i)=k
	  else
	   newdet=newdet+1
	   indet(i)=idcount+newdet
	   do j=1,nelec
	    inewdet(newdet,j)=idet(i,j)
	   enddo
	  endif
	 enddo
	endif

	return
	end	

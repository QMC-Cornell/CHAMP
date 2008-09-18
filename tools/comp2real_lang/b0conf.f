*******************************************************************************
	subroutine b0configs(nu, nd, ltot)
	parameter(max=15)
	integer norb(2000,max),morb(2000,max),norb0(3000,max),morb0(3000,max)
	integer neshl(2,max),nesh(2,max),nds(0:10),idscnt(10)
	integer nexct
	character*35    dstrfile
	logical eq_conf,newdstr
	common/bfield/  ibfield,dstrfile

c nushlf/ndshlf: the shell before the first non-full shell for up/down-spin electrons
c nushlt/ndshlt: number of total shells populated for up/down-spin electron
c neshl: number of up/down-spin electron in each shell

	nushlf=int(sqrt(2.0*nu))
	if (nushlf*(nushlf+1)/2 .gt. nu) nushlf=nushlf-1
	nu1=nu-nushlf*(nushlf+1)/2
	if (nu1 .eq. 0) then
	 nushlt=nushlf
	else
	 nushlt=nushlf+1
	endif
	do i=1,nushlf
	neshl(1,i)=i
	enddo
	do i=nushlf+1,max
	neshl(1,i)=0
	enddo
	if (nu1 .ne. 0) neshl(1,nushlt)=nu1

	ndshlf=int(sqrt(2.0*nd))
	if (ndshlf*(ndshlf+1)/2 .gt. nd) ndshlf=ndshlf-1
	nd1=nd-ndshlf*(ndshlf+1)/2
	if (nd1 .eq. 0) then
 	 ndshlt=ndshlf
	else
  	 ndshlt=ndshlf+1
	endif
	do i=1,ndshlf
	neshl(2,i)=i
	enddo
	do i=ndshlf+1,max
	neshl(2,i)=0
	enddo
	if (nd1 .ne. 0) neshl(2,ndshlt)=nd1

	lutot=0
	ldtot=0
	do i=1,nu1
	 lutot=lutot+(nushlt-1)-2*(i-1)
	enddo
	do i=1,nd1
	 ldtot=ldtot+(ndshlt-1)-2*(i-1)
	enddo
	ltotm=lutot+ldtot
	write (*,*) 'max. of total L for non-interacting GS configurations: Ltotmax=', ltotm
	write (*,*) 'choose Ltot>=0: Ltotmax, Ltotmax-2, Ltotmax-4...'
	read (*,*) ltot
	if (ltot.lt.0 .or. ltot .gt. ltotm .or. mod(ltotm-ltot,2).eq.1) then
	 write(*,*) 'wrong Ltot value'
	 stop
	endif

	call get_job_name
	
	open(unit=9,file='lorb.tmp')
	open(unit=10,file='lshl.tmp')
	open(unit=11,file=dstrfile)
	open(unit=12,file='orbt.tmp')
c	open(unit=19,file=conffile)

	call writedstr(neshl,nushlf,ndshlf,nushlt,ndshlt,ltot)
	nds(0)=1

	write(*,*) '# of excitements included, choose from 0, 2, 4,...'
	write(*,*) '0 for GS confg., 2 is usually sufficeint, play with 4 and up to find bugs :)'
	read(*,*) nexct
	if (mod(nexct,2).eq.1) then
	 write(*,*) 'should be multiple of 2'
	 stop
c	else if (nexct .eq. 0) then
c	 goto 8
	else
c 	 call move2(neshl,nushlf,ndshlf,nushlt,ndshlt,ltot)
c	 if (nexct .eq. 4) then
	do iexct=1,nexct/2
	  close(11)
	  open(unit=21,file=dstrfile)
	  open(unit=11,file='dstrb.tmp')
	  if (iexct.gt.1) then
c don't repeat dist. with 0,2,...,nexct-2 excitements
	  do idstr=1,nds(iexct-2)
	   read(21,'(5i3)',end=4) 	ndscnt,nushlf,ndshlf,nushlt,ndshlt
	   read(21,'(20i3)') (neshl(1,i), i=1,nushlt)
	   read(21,'(20i3)') (neshl(2,i), i=1,ndshlt)
	  enddo
	  endif
	  do idstr=1,500
	   read(21,'(5i3)',end=4) 	ndscnt,nushlf,ndshlf,nushlt,ndshlt
c flash out the old distribution
	   do i=1,max
	    neshl(1,i)=0
	    neshl(2,i)=0
	   enddo
	   read(21,'(20i3)') (neshl(1,i), i=1,nushlt)
	   read(21,'(20i3)') (neshl(2,i), i=1,ndshlt)
	   call move2(neshl,nushlf,ndshlf,nushlt,ndshlt,ltot)
	  enddo
c record # of distributions with up to nexct-2 excitements
4	  nds(iexct-1)=ndscnt
	  idscnt(iexct)=0
	  rewind(11)
	  do idstr=1,500
	   read(11,'(5i3)',end=7) ndscnt,nushlf,ndshlf,nushlt,ndshlt
	   read(11,'(20i3)') (neshl(1,i), i=1,nushlt)
	   read(11,'(20i3)') (neshl(2,i), i=1,ndshlt)
	   newdstr=.true.
	   rewind(11)
	   do jdstr=1,idstr-1
	    read(11,'(5i3)') ndscnt,nushf,ndshf,nusht,ndsht
	    read(11,'(20i3)') (nesh(1,j), j=1,nusht)
	    read(11,'(20i3)') (nesh(2,j), j=1,ndsht)
	    if (nushlf.eq.nushf .and. ndshlf.eq.ndshf .and. nushlt.eq.nusht .and. ndshlt.eq.ndsht) then
		 do i=1,nushlt
		  if (neshl(1,i).ne.nesh(1,i)) then
		   goto 5
		  endif
		 enddo
		 do i=1,ndshlt
		  if (neshl(2,i).ne.nesh(2,i)) then
		   goto 5
		  endif
		 enddo
		 newdstr=.false.
	    endif
5	   enddo
	   read(11,'(5i3)') ndscnt,nushf,ndshf,nusht,ndsht
	   read(11,'(20i3)') (neshl(1,j), j=1,nushlt)
	   read(11,'(20i3)') (neshl(2,j), j=1,ndshlt)
	   if(newdstr) then
	    idscnt(iexct)=idscnt(iexct)+1
	    write(21,'(5i3)') nds(iexct-1)+idscnt(iexct),nushlf,ndshlf,nushlt,ndshlt
	    write(21,'(20i3)') (neshl(1,i), i=1,nushlt)
		if (ndshlt.ne.0) then
			write(21,'(20i3)') (neshl(2,i), i=1,ndshlt)
		else
			write(21,'("  0")')
		endif
	   endif

	  enddo
7	  close(11)
	  close(21)
	  open(unit=11,file=dstrfile)
	 enddo
	endif

8	rewind(11)
	do idstr=1,500
	 read(11,'(5i3)',end=9) 	ndscnt,nushlf,ndshlf,nushlt,ndshlt
c flash out the old distribution
	 do i=1,max
	  neshl(1,i)=0
	  neshl(2,i)=0
	 enddo
	 read(11,'(20i3)') (neshl(1,i), i=1,nushlt)
	 read(11,'(20i3)') (neshl(2,i), i=1,ndshlt)
	 write(*,*) ndscnt
	 call occup(neshl,nushlf,ndshlf,nushlt,ndshlt,ltot)
	enddo
      write(*,*) 'more than 500 distributions, increase do loop limit'
	stop
9	continue
	
c read in all raw configurations
	rewind(12)
	ntot=nu+nd
	nconf0=0
	do i=1,2999
	 laccum=0
	 do j=1,ntot
	  read(12,*,end=10) norb0(i,j),morb0(i,j)
	  laccum=laccum+morb0(i,j)
	 enddo
	 if (laccum .ne. ltot) then
	  write(*,*) 'Total L is wrong, check orbitals:', i
	  do j=1,ntot
	   write(*,*) norb0(i,j),morb0(i,j)
	  enddo
	  do j=1,ntot
	   write(*,*) norb0(i-1,j),morb0(i-1,j)
	  enddo
	  stop
	 endif
	 nconf0=nconf0+1
	enddo
      write(*,*) 'more than 2000 configurations, increase do loop limit'
	stop
10	continue
	write(*,*) 'totoal # of configuration (with repeat)',nconf0

c get rid of redundant configurations
	do j=1,ntot
	 norb(1,j)=norb0(1,j)
	 morb(1,j)=morb0(1,j)
	enddo
	nconf=1
	do i=2, nconf0
	 do iconf=1,nconf
	  eq_conf=.true.
	  do j=1,ntot
	   if((morb0(i,j).ne.morb(iconf,j)) .or. (norb0(i,j).ne.norb(iconf,j))) then
	    eq_conf=.false.
	    goto 11
	   endif
	  enddo
11	  if(eq_conf) goto 12
	 enddo
12	 if(.not.eq_conf) then
c uncomment the next 3 lines to get rid of n>=2 orbitals
c	  do j=1,ntot
c	   if (norb0(i,j).ge.2) goto 20
c	  enddo
	  nconf=nconf+1
	  do j=1,ntot
	   norb(nconf,j)=norb0(i,j)
	   morb(nconf,j)=morb0(i,j)
	  enddo
	 endif
20	enddo
	write(*,*) 'new # of configuration',nconf

	write(19,'(i5,3x,"nconf")')nconf
	do i=1,nconf
	 do j=1,ntot
	  write(19,'(2i5)') norb(i,j), morb(i,j)
	 enddo
	 write(19,'(/)')
	enddo
	
	close(9)
	close(10)
	close(11)
	close(12)
c	close(19)
100	return
	end
*******************************************************************************
	subroutine move2(neshl,nushlf,ndshlf,nushlt,ndshlt,ltot)
	parameter(max=20)
	integer neshl(2,max),neshln(2,max),neshln1(2,max)

c add configurations with 2 energy excitements original one

c move 2 up-spin electrons past single energy gap, by repeating moving 1 up-spin electron
c up 1 shell, the situation moving 1 up-spin electron past two energy gaps is also covered.

	do i=nushlt,1,-1
	if (neshl(1,i) .ne. 0 .and. neshl(1,i+1) .ne. i+1) then
       call shlcp(neshl,nushlf,ndshlf,nushlt,ndshlt,
     &           neshln1,nushlfn1,ndshlfn1,nushltn1,ndshltn1)

	 neshln1(1,i)=neshln1(1,i)-1
	 neshln1(1,i+1)=neshln1(1,i+1)+1
       if (i+1 .gt. nushltn1) nushltn1=i+1

	 do j=i+1,1,-1
    	 if (neshln1(1,j).ne.0 .and. neshln1(1,j+1).ne.j+1) then
c get rid of redudancy: i->i+1->i+2 and (i+1->i+2, i->i+1)
	  if (j.eq.i+1 .and. neshl(1,j).ne.0) goto 8
        call shlcp(neshl,nushlf,ndshlf,nushlt,ndshlt,
     &            neshln,nushlfn,ndshlfn,nushltn,ndshltn)

	  neshln(1,i)=neshln(1,i)-1
	  neshln(1,i+1)=neshln(1,i+1)+1
	  if (i+1 .gt. nushltn) nushltn=i+1
	  if (i .le. nushlf) nushlfn=i-1
	  neshln(1,j)=neshln(1,j)-1
	  neshln(1,j+1)=neshln(1,j+1)+1
	  if (j+1 .gt. nushltn) nushltn=j+1
	  if (j .le. nushlf) nushlfn=j-1

	  call writedstr(neshln,nushlfn,ndshlfn,nushltn,ndshltn,ltot)
	 endif
8	 enddo
	endif
	enddo

c move 1 up-spin electron and 1 down-spin electron past single energy gap
	do i=nushlt,1,-1
	do j=ndshlt,1,-1
	if (neshl(1,i) .ne. 0 .and. neshl(1,i+1) .ne. i+1 .and.
     &	neshl(2,j) .ne. 0 .and. neshl(2,j+1) .ne. j+1) then
	 call shlcp(neshl,nushlf,ndshlf,nushlt,ndshlt,
     &           neshln,nushlfn,ndshlfn,nushltn,ndshltn)

	 neshln(1,i)=neshln(1,i)-1
	 neshln(1,i+1)=neshln(1,i+1)+1
	 if (i+1 .gt. nushltn) nushltn=i+1
	 if (i .le. nushlf) nushlfn=i-1
	 neshln(2,j)=neshln(2,j)-1
	 neshln(2,j+1)=neshln(2,j+1)+1
	 if (j+1 .gt. ndshltn) ndshltn=j+1
	 if (j .le. ndshlf) ndshlfn=j-1

	 call writedstr(neshln,nushlfn,ndshlfn,nushltn,ndshltn,ltot)
	endif
	enddo
	enddo

c move 2 down-spin electrons past single energy gap. And by consecutively moving 1 down-spin electron
c up 1 shell (j=i+1), the situation moving 1 down-spin electron past two energy gaps is also covered.

	do i=ndshlt,1,-1
	if (neshl(2,i) .ne. 0 .and. neshl(2,i+1) .ne. i+1) then
       call shlcp(neshl,nushlf,ndshlf,nushlt,ndshlt,
     &           neshln1,nushlfn1,ndshlfn1,nushltn1,ndshltn1)

	 neshln1(2,i)=neshln1(2,i)-1
	 neshln1(2,i+1)=neshln1(2,i+1)+1
       if (i+1 .gt. ndshltn1) ndshltn1=i+1

	 do j=i+1,1,-1
    	 if (neshln1(2,j).ne.0 .and. neshln1(2,j+1).ne.j+1) then
c get rid of redudancy: i->i+1->i+2 and (i+1->i+2, i->i+1)
	  if (j.eq.i+1 .and. neshl(2,j).ne.0) goto 9
        call shlcp(neshl,nushlf,ndshlf,nushlt,ndshlt,
     &            neshln,nushlfn,ndshlfn,nushltn,ndshltn)

	  neshln(2,i)=neshln(2,i)-1
	  neshln(2,i+1)=neshln(2,i+1)+1
	  if (i+1 .gt. ndshltn) ndshltn=i+1
	  if (i .le. ndshlf) ndshlfn=i-1
	  neshln(2,j)=neshln(2,j)-1
	  neshln(2,j+1)=neshln(2,j+1)+1
	  if (j+1 .gt. ndshltn) ndshltn=j+1
	  if (j .le. ndshlf) ndshlfn=j-1

	  call writedstr(neshln,nushlfn,ndshlfn,nushltn,ndshltn,ltot)
	 endif
9	 enddo
	endif
	enddo

	return
	end
*******************************************************************************
	subroutine shlcp(nesh,n1f,n2f,n1,n2,neshnew,n1fnew,n2fnew,n1new,n2new)
	parameter(max=20)
	integer nesh(2,max),neshnew(2,max)

c copy shell distribution

	do i=1,max
	neshnew(1,i)=nesh(1,i)
	neshnew(2,i)=nesh(2,i)
	enddo

	n1new=n1
	n2new=n2
	n1fnew=n1f
	n2fnew=n2f

	return
	end
*******************************************************************************
	subroutine writedstr(nesh,n1f,n2f,n1,n2,lt)
	integer nesh(2,20)

c write down shell distribution
	ndscnt=ndscnt+1
	write(11,'(5i3,3x,"# of dist., up/dn-spin full shell, up/dn-spin total shell")')
     &	ndscnt,n1f,n2f,n1,n2	
	write(11,'(<n1>i3,3x,"polulation of up-spin ele. in shells")') (nesh(1,i), i=1,n1)
	if (n2.ne.0) then
	 write(11,'(<n2>i3,3x,"polulation of down-spin ele. in shells")') (nesh(2,i), i=1,n2)
	else
	 write(11,'("  0   polulation of down-spin ele. in shells")')
	endif

	return
	end
*******************************************************************************
	subroutine occup(nesh,n1f,n2f,n1,n2,lt)
	integer nesh(2,20),ne(20),lm(20),lshl(20)
	common/orb_zero/no0cnt

c count open shells (for up/down-spin lectrons)
c lm: maximum of L in an open shell with nesh(1/2,i) electrons
	no0cnt=0
	do i=1,max(n1,n2)
	 if ((i-nesh(1,i))*nesh(1,i).ne.0) then
	 no0cnt=no0cnt+1
	 lm(no0cnt)=(i-nesh(1,i))*nesh(1,i)
	 ne(no0cnt)=nesh(1,i)
	 lshl(no0cnt)=i
	 endif
	 if ((i-nesh(2,i))*nesh(2,i).ne.0) then
	 no0cnt=no0cnt+1
	 lm(no0cnt)=(i-nesh(2,i))*nesh(2,i)
	 ne(no0cnt)=nesh(2,i)
	 lshl(no0cnt)=i
	 endif
	enddo
	ncnt=no0cnt

	nshl=max(n1,n2)
	ltot=0
	if (ncnt.ne.0)	then
	 call accum(nesh,nshl,ncnt,lm,ne,lshl,ltot,lt)
	else if (lt.eq.0) then
	 do i=1,n1
	     call writefsh(i,nesh(1,i),nesh(2,i))
	 enddo
	 write(12,'(/)')
	endif

	return
	end

*******************************************************************************
	RECURSIVE subroutine accum(nesh,nsh,ncnt,lm,ne,lshl,ltot,lt)
	integer lm(20),lshl(20),ne(20),lgood(20),nesh(2,20)
	integer norb(20,10,20),lorb(20,10,20),nset(20)
	common/orb_zero/no0cnt,lgood

c find out all possible lgood (sum of l in an oepn shell) cominations to add up to
c total angular momentum lt
	do i=-lm(ncnt),lm(ncnt),2
	 ltot=ltot+i
	 lgood(ncnt)=i
	 if (ncnt.ne.1) then
	  call accum(nesh,nsh,ncnt-1,lm,ne,lshl,ltot,lt)
	  ltot=ltot-i
	 else
	  if (ltot.eq.lt) then
	   do j=1,no0cnt
	    write(10,'(3i5)') lshl(j),ne(j),lgood(j)
	    lbegin=0
	    call orb(1-lshl(j),lshl(j),ne(j),lbegin,lgood(j))
	    nset(j)=0
	    rewind(9)
	    do k=1,10
		 do n=1,ne(j)	
	      read(9,*,end=10) norb(j,k,n), lorb(j,k,n)
	     enddo
	     nset(j)=nset(j)+1
	    enddo
	    write(*,*) 'more than 10 sets, increase do loop limit'
	    stop
10	    rewind(9)
c          write(*,*) nset(j)
	   enddo
	   write(10,'(/)')

	   nopensh=no0cnt
	   call writeorb(nesh,nsh,nopensh,lshl,ne,nset,norb,lorb)

	  endif
	  ltot=ltot-i
	 endif
	enddo

	return
	end

*******************************************************************************
	RECURSIVE subroutine orb(lmin,ish,nele,ltot,lt)
	integer norb(10),lorb(10)
	common/nele_shl/n0,norb,lorb

c find out electron occupation in an open shell to add up to
c a sum. of lt
	if (lmin .eq. 1-ish) n0=nele
	do i=lmin,ish-1,2
	 ltot=ltot+i
	 lorb(nele)=i
	 norb(nele)=(ish-1-abs(i))/2
	 if (nele.ne.1) then
	  call orb(lorb(nele)+2,ish,nele-1,ltot,lt)
	  ltot=ltot-i
	 else
	  if(ltot.eq.lt) then
c	   write(9,'(i3,3x,"shell")') ish
c	   write(9,'(i3,3x,"electrons@shell")') n0
	   do n=1,n0
	    write(9,'(2i5)') norb(n),lorb(n)
	   enddo
c	   write(9,'(/)')	
	  endif
	  ltot=ltot-i
	 endif
	enddo

	return
	end
*******************************************************************************
	RECURSIVE subroutine writeorb(nesh,nsh,nopensh,lshl,ne,nset,norb,lorb)
	integer norb(20,10,20),lorb(20,10,20),nset(20),lshl(20),ne(20),nesh(2,20)
	integer norb1(20),lorb1(20),norb2(20),lorb2(20),iset(20)
	common/orb_zero/no0cnt
	common/orb_set/iset

c write down electron occupation
	do i=1,nset(no0cnt+1-nopensh)
	 iset(no0cnt+1-nopensh)=i
	 if (nopensh .ne. 1) then
	  call writeorb(nesh,nsh,nopensh-1,lshl,ne,nset,norb,lorb)
	 else
	  do ish=1,nsh
	    npartial=0
	    do j=1,no0cnt
	     if (ish .eq. lshl(j)) then
		  npartial=npartial+1
	      if(npartial.eq.1) then
	       do n=1,ne(j)
	       norb1(n)=norb(j,iset(j),n)
		   lorb1(n)=lorb(j,iset(j),n)
		   enddo
		  else if (npartial.eq.2) then
	       do n=1,ne(j)
	       norb2(n)=norb(j,iset(j),n)
		   lorb2(n)=lorb(j,iset(j),n)
		   enddo
		  endif
	     endif
		enddo
	    if(npartial .ge. 3) then
		 write(*,*) 'check the open shells'
		 stop
	    endif
		if(npartial .eq. 0) then
	     call writefsh(ish,nesh(1,ish),nesh(2,ish))
	    elseif (npartial .eq. 1) then
	     call writesh1(ish,nesh(1,ish),nesh(2,ish),norb1,lorb1)
	    elseif (npartial .eq. 2) then
	     call writesh2(ish,nesh(1,ish),nesh(2,ish),norb1,lorb1,norb2,lorb2)
	    endif
	    if (ish.eq.nsh) write(12,'(/)')
10	  enddo
	 endif
	enddo

	return
	end
*******************************************************************************
	subroutine writefsh(ish,nu,nd)

c write down a full/half-full shell
	if (nu.eq.ish .and. nd.eq.ish) then
	   do j=1,(ish+1)/2
	   write(12,'(2i6)') j-1,ish-1-2*(j-1)
	   write(12,'(2i6)') j-1,ish-1-2*(j-1)
	   if (ish-1-2*(j-1).ne.0) then
	    write(12,'(2i6)') j-1,-(ish-1)+2*(j-1)
	    write(12,'(2i6)') j-1,-(ish-1)+2*(j-1)
	   endif
	   enddo
	else if (nu.eq.ish .or. nd.eq.ish) then
	   do j=1,ish/2
	   write(12,'(2i6)') j-1,ish-1-2*(j-1)
	   if (ish-1-2*(j-1).ne.0) write(12,'(2i6)') j-1,-(ish-1)+2*(j-1)
	   enddo
	   if (mod(ish,2).eq.1) write(12,'(2i6)') ish/2,0
	endif

	return
	end
*******************************************************************************
	subroutine writesh1(ish,nu,nd,norb,lorb)
	integer lorb(10),norb(10)

c write down one open shell(up/down spin)
	if(nu.eq.0 .or. nd.eq.0) then
	 ne=nu+nd
	 do j=1,(ish+1)/2
	  do i=1,ne
	   if(norb(i).eq.j-1 .and. abs(lorb(i)).eq.ish-1-2*(j-1)) write(12,'(2i6)') norb(i),lorb(i)
	  enddo
	 enddo
	elseif(nu.eq.ish .or. nd.eq.ish) then
	 ne=nu+nd-ish
	 do j=1,(ish+1)/2
	  write(12,'(2i6)') j-1,ish-1-2*(j-1)
	  do i=1,ne
	   if(norb(i).eq.j-1 .and. abs(lorb(i)).eq.ish-1-2*(j-1)) write(12,'(2i6)') norb(i),lorb(i)
	  enddo
	  if (ish-1-2*(j-1).ne.0) write(12,'(2i6)') j-1,-(ish-1)+2*(j-1)
	 enddo
	endif

	return
	end
*******************************************************************************
	subroutine writesh2(ish,nu,nd,norb1,lorb1,norb2,lorb2)
	integer norb1(10),lorb1(10),norb2(10),lorb2(10)

c write down two open shells(up & down spin)
	do j=1,(ish+1)/2
	  do i=1,nu
	   if(norb1(i).eq.j-1 .and. lorb1(i).eq.ish-1-2*(j-1)) write(12,'(2i6)') norb1(i),lorb1(i)
	  enddo
	  do i=1,nd
	   if(norb2(i).eq.j-1 .and. lorb2(i).eq.ish-1-2*(j-1)) write(12,'(2i6)') norb2(i),lorb2(i)
	  enddo
	 if (ish-1-2*(j-1).ne.0) then
	  do i=1,nu
	   if(norb1(i).eq.j-1 .and. -lorb1(i).eq.ish-1-2*(j-1)) write(12,'(2i6)') norb1(i),lorb1(i)
	  enddo
	  do i=1,nd
	   if(norb2(i).eq.j-1 .and. -lorb2(i).eq.ish-1-2*(j-1)) write(12,'(2i6)') norb2(i),lorb2(i)
	  enddo
	 endif
	enddo

	return
	end
*******************************************************************************

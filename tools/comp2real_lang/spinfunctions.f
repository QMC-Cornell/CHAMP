      subroutine Spin_eigenfunctions(NA,NB,IS)
      implicit real*8(a-h,o-z)
	parameter(max=15)
	dimension IPR(2000,max)
	dimension Nib(1000,max)


c	PRINT *,'number of up electrons='
c	READ *, NA
c	PRINT *,'number of down electrons='
c	READ *, NB
c	PRINT *,'2*Spin ='
c	READ *, IS
	N=NA+NB
c	IS=NA-NB

	write(*,*)'Total Number of electrons=',N
      write(*,*)'Total Spin S=             ',dfloat(IS)*0.5d0
      write(*,*)'Total Spin projection S_z=',dfloat(NA-NB)*0.5d0

      call primitiv_spinfkts(NA,NB,iprim,ipr)

      call branching_symbols(N,IS,idim,Nib)

      call S2_coeff(NA,NB,IS,iprim,idim,Nib,ipr)


      return
	end


**********************************************************************************************
**********************************************************************************************
*** From"The construction of Spin eigenfunctions",Ruben Pauncz,Kluwer Academic/Plenum Publishers
      subroutine  branching_symbols(Nelec,ISpin,ndim,Nib)
**** A Spin S eigenfunction Psi (S^2 Psi=S(S+1)Psi)
**** can be uniquely represented by a ordered series of branching-symbols
**** B(N,S,I,J) N=number of electrons
****            S = total Spin number
****            I = 1...f(N,S) = degeneracy of state with Spin S
****            J = 1...N

**** Example: for N=3, S=1/2, one can have a state
**** from  1/2 +1/2 -1/2 => 1 1 2
**** and   1/2 -1/2 +1/2 => 1 2 1
**** which follwos from the addition and subtraction equations (1),(2):

**** X(N,S+1/2,M+1/2) = [(S+M+1)^0.5 X(N-1,S,M)alpha(N)
****                   -(S-M)^0.5 X(N-1,S,M+1)beta(N)     (1)

**** X(N,S-1/2,M+1/2) = [-(S-M)^0.5 X(N-1,S,M)alpha(N)
****                   +(S+M+1)^0.5 X(N-1,S,M+1)beta(N)   (2)

**** We therefore write IB(3,2*S=1,1,1)=1
****                    IB(3,2*S=1,1,2)=1
****                    IB(3,2*S=1,1,3)=2
****                and
****                    IB(3,2*S=1,2,1)=1
****                    IB(3,2*S=1,2,2)=2
****                    IB(3,2*S=1,2,3)=1

      implicit real*8(a-h,o-z)
	parameter(max=15)
	dimension Nib(1000,max)
	dimension jf(max,0:max),ib(max,0:max,1000,max)
	open(6,file='branch1.lis',status='replace',access='sequential')
	open(5,file='branch.dat',status='replace',access='sequential')
	open(7,file='branch2.lis',status='replace',access='sequential')

      nm=nelec
****** jf(N,S) is the degeneracy of Spin-eigenfunctions of N particles and Spin S
*******             N          N
****** jf(N,S)=  1/2N-S  -  1/2N-S-1
      jf(2,2)=1
	jf(2,0)=1
      do i=3,nm
	jf(i,i)=1
	enddo
	do 10 n=3,nm,2
	do i=1,n-1,2
	jf(n,i)=jf(n-1,i-1)+jf(n-1,i+1)
	enddo
	n1=n+1
	jf(n1,0)=jf(n,1)
	do i=1,n-1
	jf(n1,i)=jf(n,i-1)+jf(n,i+1)
	enddo
10	continue

      do i=2,nm,2
	write(6,100)i,(jf(i,j),j=0,i,2)
100    format(6x,i2,<i>i7)
      i1=i+1
	write(6,101)i,(jf(i1,j),j=0,i1,2)
101    format(6x,i2,<i1>i7)
      enddo
	
	write(6,*)'***************************************'
	write(6,*)'The dimensions are printed out'
	write(6,*)'***************************************'

**** Calculating the branching diagram symbols

***** Two electron triplet state
      ib(2,2,1,1)=1
      ib(2,2,1,2)=1

***** Two electron singlet state
      ib(2,0,1,1)=1
      ib(2,0,1,2)=2
***** n=3 is=1
      ib(3,1,1,1)=1
      ib(3,1,1,2)=1
      ib(3,1,1,3)=2

      ib(3,1,2,1)=1
      ib(3,1,2,2)=2
      ib(3,1,2,3)=1
***** n=3 is=3
      ib(3,3,1,1)=1
      ib(3,3,1,2)=1
      ib(3,3,1,3)=1

      do n=3,nm
	do j=1,n
	ib(n,n,1,j)=1
	enddo
	enddo
***** even number of electrons
      do 11 n=4,nm,2
	do 12 is=0,n,2
	write(6,110)n,is
110	format(/,12x,'n=',i2,6x,'2*S=',i2,/)
	if1=jf(n-1,is+1)
	do k=1,if1
	do j=1,n-1
	ib(n,is,k,j)=ib(n-1,is+1,k,j)
	enddo
	ib(n,is,k,n)=2
	write(6,120)k,(ib(n,is,k,j),j=1,n)
120	format(6x,i4,6x,<n>i2)
	enddo
      if(is.eq.0)goto 12
	if2=jf(n-1,is-1)
	do k=1,if2
	kk=if1+k
	do j=1,n-1
	ib(n,is,kk,j)=ib(n-1,is-1,k,j)
	enddo
	ib(n,is,kk,n)=1
	write(6,120)kk,(ib(n,is,kk,j),j=1,n)
	enddo
12    continue

***** odd number of electrons
      n1=n+1
	do 13 is=1,n1,2
	write(6,110)n1,is
	if1=jf(n,is+1)
	do k=1,if1
	do j=1,n
	ib(n1,is,k,j)=ib(n,is+1,k,j)
	enddo
	ib(n1,is,k,n1)=2
	write(6,121)k,(ib(n1,is,k,j),j=1,n1)
	enddo
121	format(6x,i4,6x,<n1>i2)
	if2=jf(n,is-1)
	do k=1,if2
	kk=if1+k
	do j=1,n
	ib(n1,is,kk,j)=ib(n,is-1,k,j)
	enddo
	ib(n1,is,kk,n1)=1
	write(6,121)kk,(ib(n1,is,kk,j),j=1,n1)
	enddo
13	continue
11	continue

      close(5)
      close(6)

      ndim=jf(nm,ISpin)

      do i=1,nm
	do j=1,ndim
	  Nib(j,i)=ib(nm,Ispin,j,i)
	enddo
	enddo
	write(7,110)n,ispin

      do k=1,jf(nm,ISpin)
	write(7,120)k,(Nib(k,j),j=1,nm)
	enddo
122	format(7x,i4,6x,<nm>i2)

	write(*,*)'end branchin_symbols'

	return
	end		

*********************************************************************
*********************************************************************
      subroutine  primitiv_spinfkts(NA,NB,idim,ipath)
*** Primitive functions for a certain S_z= (NA-NB)*1/2
*** are given by the arramgement of spin up and downs.
*** There are N!/(NA! NB!) primitive spin functions, N=NA+NB
*** The primitve functions are  eigenfunctions of S_z but in general
*** not eigenfunction of S
*** The primitive fucntions are used to construct eigenfunctions of S
      implicit real*8(a-h,o-z)
	parameter(max=15)
	dimension j(0:max),ibranch(1000,max),ipar(max)
	dimension ibeta(2000,max),index(2000,max),ipath(2000,max)	
	open(6,file='prim.lis',status='replace',access='sequential',form='formatted')

c	PRINT *,'number of up electrons='
c	READ *, NA
c	PRINT *,'number of down electrons='
c	READ *, NB

	n=na+nb
	write(6,100)na,nb,n
100	format(6x,'na=',i2,6x,'nb=',i2,6x,'n=',i2)
      icount=0

	j(0)=0
	do i=1,nb
	j(i)=na
	enddo

20	icount=icount+1
	do i=1,nb
	index(icount,i)=j(i)
	ibeta(icount,i)=j(i)+i
	enddo

	do 30 i=1,nb
	if(j(i).gt.j(i-1))then
	   j(i)=j(i)-1
	   do k=1,i-1
         j(k)=j(i)
	   enddo
	   goto 20
	endif
30    continue

      idim=icount !idim is the number of primitive functions
101   format(6x,60('='))
      write(6,101)
      write(6,104)
104   format(16x,'index',6x,'beta places',6x,'path diagram symbols')
      write(6,101)

	do i=1,idim
	do k=1,n
	ipath(i,k)=1
	enddo
	enddo

	do i=1,idim
	do k=1,nb
	ik=ibeta(i,k)
	ipath(i,ik)=2
	enddo
	enddo
	icount=0
	do 300 i=1,idim
	mm=0
	do 310 k=1,n
	if(ipath(i,k).eq.1)goto 302
	mm=mm-1
	if(mm.lt.0)goto 300
	goto 304

302   mm=mm+1
304   ipar(k)=mm
310   continue

      icount=icount+1
	do k=1,n
	ibranch(icount,k)=ipath(i,k)
	enddo
300   continue

      do i=1,idim
	write(6,103)i,(index(i,k),k=1,nb),(ibeta(i,k),k=1,nb),(ipath(i,k),k=1,n)
103   format(4x,i4,6x,2i3,6x,2i3,6x,6i3)
      enddo
      close(6)
	return
	end
	
*********************************************************************
*********************************************************************
      subroutine  S2_coeff(NA,NB,IS,iprim,idim,Nib,ipr)
*** Find the S eigenfunctions for a system with given NA,NB,S.
*** Supply all irpim=N!/NA!NB! primitive spin functions,
*** and all idim=F(N,S) branching symbols


*** Some background:
*** Write

**** X(N,S,M,Bi) = 1/(2S)^0.5  [(S+M)^0.5   X(N-1,S-1/2,M-1/2;Bi')alpha(N)
****                            -(S-M)^0.5  X(N-1,S-1/2,M+1/2;Bi')beta(N)
****             = C(1,1,S,M) X(N-1,S-1/2,M-1/2,Bi')alpha(N)
****             + C(1,2,S,M) X(N-1,S-1/2,M+1/2,Bi')beta(N)  (1)

**** X(N,S,M,Bi) = 1/(2S+2)^0.5[-(S-M+1)^0.5 X(N-1,S+1/2,M-1/2;Bi') alpha(N)
****              +              (S+M+1)^0.5 X(N-1,S+1/2,M+1/2;Bi') beta(N)
****             = C(2,1,S,M) X(N-1,S+1/2,M-1/2,Bi')alpha(N)
****             + C(2,2,S,M) X(N-1,S+1/2,M+1/2,Bi')beta(N)  (1)

**** In order to uniquely specify S eigenfunctions for N particles
**** without keeping track of all the corresponding N-1 S+- 1/2 eigenfunctions,
**** one can use the indexing of the C's: 1st index: 1 for addition, 2 for subtraction
****                                      2nd index: 1 for addition alpha, 2 for beta

      implicit real*8(a-h,o-z)
	parameter(max=15)
	dimension Nib(1000,max),IB(max),IP(max),IPR(2000,max),coeff(2000),icsfconf(2000)
	real*8 in,id
	real coefa(1000), coefb(1000)
	dimension norb(max),lorb(max),morb(max),ieorb(max),iexch(2000),nup(2000,max),ndn(2000,max)
	dimension ndeto(4000),coeffo(4000,100),norbo(4000,100,max),morbo(4000,100,max),icsfconfo(4000)	
	logical eq_prim, inontrival, eq_csf
	character*35    dstrfile
	common/bfield/  ibfield,dstrfile
	common/csf_conf/ icsfconf

c	open(19,file='orb_in_det',status='old',access='sequential',form='formatted')
c	open(16,file='ssquare.in',status='replace',access='sequential',form='formatted')
	

	
c	open(5,file='coeff.dat',status='old',access='sequential',form='formatted')

      N=NA+NB
	one=1.d0
	coefflimit=1.d-6
	ncsf=0

c read in orbitals occupied from *.conf file
	iemax=1
	rewind(19)
	read(19,*) nconf

	DO iconf=1,nconf
	 if (ibfield .eq. 1) then
	  do j=1,n
	  read(19,*) norb(j), morb(j)
	  lorb(j)=abs(morb(j))
	  write(*,*) ' n, l, m'
	  write(*,*) norb(j), lorb(j), morb(j)
	  enddo
	 else
	  read(19,*) (morb(i),i=1,n)
	  do j=1,n
	  lorb(j)=abs(morb(j))
	  norb(j)=0
	  write(*,*) ' n, l, m'
	  write(*,*) norb(j), lorb(j), morb(j)
	  ieorb(j)=2*norb(j)+lorb(j)+1
	  if(ieorb(j) .gt. iemax) iemax=ieorb(j)
	  enddo
	 endif

	open(6,file='coeff.lis',status='replace',access='sequential',form='formatted')

	open(28,file='csf.tmp',status='replace',access='sequential',form='formatted')

c	read(5,*)n,iprim,idim
	write(6,100)n,IS,NA-NB,iprim,idim
100   format(6x,'N=',i3,6x,'2*S=',i3,6x,'2*Sz=',i3,6x,'IPRIM=',i4,6x,'IDIM=',i4)
      write(6,*)'Primitiv functions as basis for S eigenfunctions (1=up,2=down)'
**** Read all possible primitive functions iprim=N!/(NA! NB!)
      do i=1,iprim
c	read(5,*)(ipr(i,j),j=1,n)
	write(6,120)i,(ipr(i,j),j=1,n)
120   format(5x,i4,6x,50i2)
	enddo

c count number of CSF for each configuration 	
	ncsf1=0

**** Read all possible branching symbols idim=f(N,S)=degenracy of Spin eigenstate S
	write(6,*)'=============================================='

	do 300 ibranch=1,idim

c	read(5,*)(ib(i),i=1,n)
	do j=1,N
	 ib(j)=Nib(ibranch,j)
	enddo

101   format(6x,'Branching diagram: ',12i2)
      write(6,*)'**********************************************'

	do 200 icount=1,iprim
	do k=1,n
	ip(k)=ipr(icount,k)
	enddo
	write(6,102)icount,(ip(i),i=1,n)
102   format(5x,'inum=',i4,6x,'Path diagram: ',15i2)

	in=1
	id=1
	ks=0
	km=0
	do k=1,n
	ipk=ip(k)
	ibk=ib(k)

	if(ibk.eq.2)goto 2
	ks=ks+1
	goto 3
2     ks=ks-1
3     continue

      if(ipk.eq.2)goto 22
	km=km+1
	goto 33
22    km=km-1
33    continue

      call ic(ibk,ipk,ks,km,inum,iden)
c	in=in*inum
c	id=id*iden
	in=in*inum/2.0
	id=id*iden/2.0
	enddo

c	n2=2**n
c	in=in/n2
c	id=id/n2
	a=abs(in)
	if(a.eq.0)goto 44
	su=in/a
	s=su*sqrt(a)
	goto 45
44    s=0.d0
45    continue
      write(6,107)icount,int(in),s
	coeff(icount)=s
107   format(15x,i4,6x,i8,3x,f12.6,6x,i6)
200   continue
      write(6,*)ibranch,'. Coefficients for S eigenfunction'
      write(6,109)(coeff(i),i=1,iprim)
109	format(6x,<iprim>f12.6)
      write(6,108) int(id)
108   format(6x,'Square of the normalization constant=',i10)




c set coeff to 0 for null determinants according to Pauli principle
c and calculate number of exchanges of orbitals needed to order determinant
c to D_up D-dn
c NOTE S^2 eigenfunctions won't remain orthogonal to each other after the deduction
c hence some redundant ones in final result.
	do i=1,iprim
	do j=1,n-1
	 do j1=j+1,n
	  if ((morb(j1) .eq. morb (j)) .and. (ipr(i,j1) .eq. ipr(i, j)) .and.
     &   (norb(j1) .eq. norb (j))) then
	   coeff(i)=0.0
         goto 121
	  endif
	 enddo
      enddo
121	iexch(i)=0
	iup=0
	idn=0
	if (abs(coeff(i)) .gt. coefflimit) then
	 do j=1,n
	  if (ipr(i,j) .eq. 1) then
	   iexch(i)=iexch(i)+j
	   iup=iup+1
	   nup(i,iup)=j
	  else
	   idn=idn+1
	   ndn(i,idn)=j
	  endif
	 enddo
	 iexch(i)=iexch(i)-NA*(NA+1)/2
	endif
	enddo


c combine same determinants(after exchange of orbitals)	
	do i=1, iprim-1
	if (abs(coeff(i)) .gt. coefflimit) then
	 do i1=i+1, iprim
	  if (abs(coeff(i1)) .gt. coefflimit) then
	   eq_prim=.true.
	   do k=1,NA
	    if (morb(nup(i,k)) .ne. morb(nup(i1,k)).or.norb(nup(i,k)) .ne. norb(nup(i1,k))) then
		 eq_prim=.false.
	     goto 122
	    endif
	   enddo
	   do k=1,NB
	    if (morb(ndn(i,k)) .ne. morb(ndn(i1,k)).or.norb(nup(i,k)) .ne. norb(nup(i1,k))) then
		 eq_prim=.false.
	     goto 122
	    endif
	   enddo
	  endif
122     if(eq_prim) then
	   coeff(i)=coeff(i)+coeff(i1)*(-1)**(iexch(i)+iexch(i1))
	   coeff(i1)=0.0
	  endif	
       enddo
	endif
	enddo


	
c	write(*,*) (coeff(i),i=1,iprim)
c	pause

c normalize det. coeff. and see if this CSF is trival
c inonzero: number of non-zero determinants 	
	inontrival=.false.
	coeff_norm2=0.0
	do i=1,iprim
	 coeff_norm2=coeff_norm2+coeff(i)*coeff(i)
	enddo
	if(coeff_norm2 .gt. coefflimit) then
	  inontrival=.true.
	  inonzero=0
	  coeff_norm=sqrt(coeff_norm2)
c	write(*,*) coeff(i)
	  do i=1,iprim
	   coeff(i)=coeff(i)*(-1)**iexch(i)/coeff_norm
	   if(abs(coeff(i)) .gt. coefflimit) then
	    inonzero=inonzero+1
c make sure the 1st non-zero determ. coeff. be positive, so no duplicate CSF
c with just opposite signs on determ. coeff. remains after comparison later.
	    coeff_norm=coeff_norm*(-1)**iexch(i)
	    coeff(i)=abs(coeff(i))
	    do i1=i+1,iprim
	     coeff(i1)=coeff(i1)*(-1)**iexch(i1)/coeff_norm
	     if(abs(coeff(i1)) .gt. coefflimit) inonzero=inonzero+1
	    enddo
		goto 123	
	   endif
	  enddo
	endif
	

123	ndim=2
	if(inontrival) then
	 ncsf1=ncsf1+1
	 icsfconfo(ncsf1)=iconf
	 write(28,'(3i3,6x,"ndim, nelec, ndet")') ndim,n,inonzero
	 do i=1, iprim
	  if(abs(coeff(i)) .gt. coefflimit) then
	   write(28,'(/,f12.8,6x,"cdet")') coeff(i)
         write(28,'(4i3,6x,"2s,n,l,m")') 1,norb(nup(i,1)),lorb(nup(i,1)),morb(nup(i,1))
	   do k=2,NA
	    write(28,'(4i3)') 1,norb(nup(i,k)),lorb(nup(i,k)),morb(nup(i,k))
	   enddo
	   do k=1,NB
	    write(28,'(4i3)') -1,norb(ndn(i,k)),lorb(ndn(i,k)),morb(ndn(i,k))
	   enddo
	  endif
	 enddo
	 write(28,'(//)')
	endif

300   continue
      close(6)
	rewind(28)
	do icsf=1,ncsf1
	 read(28,*) ndimo,neleco,ndeto(icsf)
	 do idet=1,ndeto(icsf)
	  read(28,*) coeffo(icsf,idet)
	  do i=1,n
	   read(28,*) ispin,norbo(icsf,idet,i),lorbo,morbo(icsf,idet,i)
	  enddo
	 enddo
	enddo
	ncsfn=ncsf1
	close(28)

c flash out duplicate CSF
	do icsf=1,ncsf1-1
	if (ndeto(icsf) .ne. 0) then
	 do icsf1=icsf+1,ncsf1
	  eq_csf=.true.
	  if(ndeto(icsf1) .ne. ndeto(icsf)) then
	   eq_csf=.false.
	   goto 310
	  endif
	   do idet=1,ndeto(icsf)
	    do i=1,n
	     if (norbo(icsf1,idet,i) .ne. norbo(icsf,idet,i) .or.
     &		 morbo(icsf1,idet,i) .ne. morbo(icsf,idet,i)) then
		  eq_csf=.false.
		  goto 310
		 endif
	    enddo
	    if (coeffo(icsf1,idet) .ne. coeffo(icsf,idet)) then
	     eq_csf=.false.
	     goto 310
		endif
	   enddo

310	  if(eq_csf) then
	   ndeto(icsf1)=0
	   ncsfn=ncsfn-1
	  endif
	 enddo
	endif
	enddo


	icsfcnt=0
	write(20,'(i3,3x,"CSF(s) constructed from configuration",i4,/)') ncsfn, iconf
	do icsf=1,ncsf1
	if (ndeto(icsf) .ne. 0) then
	 icsfcnt=icsfcnt+1
	 icsfconf(icsfcnt)=icsfconfo(icsf)
	 write(20,'(3i3,i5,6x,"ndim, nelec, ndet, icsf")') ndim,n,ndeto(icsf),icsfcnt+ncsf
	 do idet=1,ndeto(icsf)
	   write(20,'(/,f12.8,6x,"cdet")') coeffo(icsf,idet)
         write(20,'(4i3,6x,"2s,n,l,m")') 1,norbo(icsf,idet,1),abs(morbo(icsf,idet,1)),morbo(icsf,idet,1)
	   do i=2,NA
	    write(20,'(4i3)') 1,norbo(icsf,idet,i),abs(morbo(icsf,idet,i)),morbo(icsf,idet,i)
	   enddo
	   do i=1,NB
	    write(20,'(4i3)') -1,norbo(icsf,idet,i+NA),abs(morbo(icsf,idet,i+NA)),morbo(icsf,idet,i+NA)
	   enddo
	 enddo
	 write(20,'(//)')
	endif
	enddo
	ncsf=ncsf+ncsfn

	ENDDO

	write(20,'("Total number of CSFs:",i4,//)') ncsf

	rewind(20)
	write(18,'(i5,3x,"ncsf: Total number of CSFs",//)') ncsf
	do iconf=1,nconf
	 read(20,*) ncsfinconf
	 do i=1,ncsfinconf
	  read(20,*) ndim, nelec, ndet, icsf
	  write(18,'(3i3,2i5,6x,"ndim, nelec, ndet, icsf, iconf")') ndim,n,ndet,icsf,iconf
	  do idet=1,ndet
	   read(20,*) coeffi
	   write(18,'(/,f12.8,6x,"cdet")') coeffi
	   do ielec=1,n
	    read(20,*) is2,norbb,lorbb,morbb
		write(18,'(4i3)') is2,norbb,lorbb,morbb
	   enddo
	  enddo
	  write(18,'(//)')
	 enddo
	enddo
	close(20)

      return
	end
*********************************************************************
*********************************************************************
      subroutine  IC(ibb,ipp,is,im,inum,iden)
      implicit double precision(a-h,o-z)
	if(ibb.eq.1) goto 1
	goto 2
1	if(ipp.eq.1) goto 11
	goto 12
2     if(ipp.eq.1)goto 21
      goto 22
11    continue
      inum=is+im
	iden=is*2
	goto 3
12    continue
      inum=is-im
	iden=is*2
	goto 3
21    continue
      inum=-(is-im+2)
	iden=(is+2)*2
	goto 3
22    continue
      inum=is+im+2
	iden=(is+2)*2
3	return
	end

**********************************************************************************************
**********************************************************************************************
*** From"The construction of Spin eigenfunctions",Ruben Pauncz,Kluwer Academic/Plenum Publishers
      subroutine  branching_symbols_old(Nelec,ISpin,ndim,Nib)
**** A Spin S eigenfunction Psi (S^2 Psi=S(S+1)Psi)
**** can be uniquely represented by a ordered series of branching-symbols
**** B(N,S,I,J) N=number of electrons
****            S = total Spin number
****            I = 1...f(N,S) = degeneracy of state with Spin S
****            J = 1...N

**** Example: for N=3, S=1/2, one can have a state
**** from  1/2 +1/2 -1/2 => 1 1 2
**** and   1/2 -1/2 +1/2 => 1 2 1
**** which follwos from the addition and subtraction equations (1),(2):

**** X(N,S+1/2,M+1/2) = [(S+M+1)^0.5 X(N-1,S,M)alpha(N)
****                   -(S-M)^0.5 X(N-1,S,M+1)beta(N)     (1)

**** X(N,S-1/2,M+1/2) = [-(S-M)^0.5 X(N-1,S,M)alpha(N)
****                   +(S+M+1)^0.5 X(N-1,S,M+1)beta(N)   (2)

**** We therefore write IB(3,2*S=1,1,1)=1
****                    IB(3,2*S=1,1,2)=1
****                    IB(3,2*S=1,1,3)=2
****                and
****                    IB(3,2*S=1,2,1)=1
****                    IB(3,2*S=1,2,2)=2
****                    IB(3,2*S=1,2,3)=1

      implicit real*8(a-h,o-z)
	parameter(max=15)
	dimension Nib(1000,max)
	dimension jf(max,0:max),ib(max,0:max,1000,max)
	open(6,file='branch1.lis',status='replace',access='sequential')
	open(5,file='branch.dat',status='replace',access='sequential')
	open(7,file='branch2.lis',status='replace',access='sequential')

      nm=nelec
****** jf(N,S) is the degeneracy of Spin-eigenfunctions of N particles and Spin S
*******             N          N
****** jf(N,S)=  1/2N-S  +  1/2N-S-1
      jf(2,2)=1
	jf(2,0)=1
      do i=3,nm
	jf(i,i)=1
	enddo
	do 10 n=3,nm,2
	do i=1,n-1,2
	jf(n,i)=jf(n-1,i-1)+jf(n-1,i+1)
	enddo
	n1=n+1
	jf(n1,0)=jf(n,1)
	do i=1,n-1
	jf(n1,i)=jf(n,i-1)+jf(n,i+1)
	enddo
10	continue

      do i=2,nm,2
	write(6,100)i,(jf(i,j),j=0,i,2)
100    format(6x,i2,<i>i7)
      i1=i+1
	write(6,101)i,(jf(i1,j),j=0,i1,2)
101    format(6x,i2,<i1>i7)
      enddo
	
	write(6,*)'***************************************'
	write(6,*)'The dimensions are printed out'
	write(6,*)'***************************************'

**** Calculating the branching diagram symbols

***** Two electron triplet state
      ib(2,2,1,1)=1
      ib(2,2,1,2)=1

***** Two electron singlet state
      ib(2,0,1,1)=1
      ib(2,0,1,2)=2
***** n=3 is=1
      ib(3,1,1,1)=1
      ib(3,1,1,2)=1
      ib(3,1,1,3)=2

      ib(3,1,2,1)=1
      ib(3,1,2,2)=2
      ib(3,1,2,3)=1
***** n=3 is=3
      ib(3,3,1,1)=1
      ib(3,3,1,2)=1
      ib(3,3,1,3)=1

      do n=3,nm
	do j=1,n
	ib(n,n,1,j)=1
	enddo
	enddo
***** even number of electrons
      do 11 n=4,nm,2
	do 12 is=0,n,2
	write(6,110)n,is
110	format(/,12x,'n=',i2,6x,'2*S=',i2,/)
	if1=jf(n-1,is+1)
	do k=1,if1
	do j=1,n-1
	ib(n,is,k,j)=ib(n-1,is+1,k,j)
	enddo
	ib(n,is,k,n)=2
	write(6,120)k,(ib(n,is,k,j),j=1,n)
120	format(6x,i4,6x,<n>i2)
	enddo
      if(is.eq.0)goto 12
	if2=jf(n-1,is-1)
	do k=1,if2
	kk=if1+k
	do j=1,n-1
	ib(n,is,kk,j)=ib(n-1,is-1,k,j)
	enddo
	ib(n,is,kk,n)=1
	write(6,120)kk,(ib(n,is,kk,j),j=1,n)
	enddo
12    continue
***** odd number of electrons
      n1=n+1
	do 13 is=1,n1,2
	write(6,110)n1,is
	if1=jf(n,is+1)
	do k=1,if1
	do j=1,n
	ib(n1,is,k,j)=ib(n,is+1,k,j)
	enddo
	ib(n1,is,k,n1)=2
	write(6,121)k,(ib(n1,is,k,j),j=1,n1)
	enddo
121	format(6x,i4,6x,<n1>i2)
	if2=jf(n,is-1)
	do k=1,if2
	kk=if1+k
	do j=1,n
	ib(n1,is,kk,j)=ib(n,is-1,k,j)
	enddo
	ib(n1,is,kk,n1)=1
	write(6,121)kk,(ib(n1,is,kk,j),j=1,n1)
	enddo
13	continue
11	continue

      close(5)
      close(6)

      ndim=jf(nm,ISpin)

      do i=1,nm
	do j=1,ndim
	  Nib(j,i)=ib(nm,Ispin,j,i)
	enddo
	enddo
	write(7,110)n,ispin

      do k=1,jf(nm,ISpin)
	write(7,120)k,(Nib(k,j),j=1,nm)
	enddo
122	format(7x,i4,6x,<nm>i2)

	write(*,*)'end branchin_symbols'

	return
	end		


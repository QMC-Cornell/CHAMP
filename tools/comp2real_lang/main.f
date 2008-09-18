***********************************************************************
      program qdot_csf
c Written by Lang Zeng, based on program by Devrim Guclu and Wolfgang Geist.
*----------------------------------------------------------------------
* title: CSF for Quantum Dot
* file: main.f filem.f order.f comp2real.f b0conf.f spinfunctions.f detcsf.f
* author: L Zeng
* date: Aug 2007
* modification:
* description: Find out possible configurations for given # of up/down-spin
*              electrons and total angular momentum L. Then produce CSFs from
*              those configurations and write out the real and imag. parts of
*              CSFs (eigenfunctions of S^2, in the form of a linear comb. of
*              Slater determinants), by converting the complex orbitals into
*              real orbitals. The results can be used for QMC calculations.
*----------------------------------------------------------------------
* LIST OF VARIABLES:
*
* job_name = job_name used to define output files
* ndim     = number of spatial dimensions
* nelec    = number of electrons
* ndet     = number of determinants (which increases during expansion)
* nconf    = number of configurations
* ncsf     = number of CSFs
* det(idet,in,qnx) = quantum number "qnx" of in^{th} electron in the
*                    idet^{th} determinant
* cdet(idet)    = coefficient of the idet^{th} determinant.
*----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      include         'maxdim.h'
      include         'qnumbers.h'
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)

	integer         norb(MAXNELEC),morb(MAXNELEC)
	integer     	morbmax(0:MAXN)
      character*35 dstrfile

      common/psi_n/   cdet,ndim,nelec,ndet,det
	common/bfield/  ibfield, dstrfile
	common/basis/  ibasis
	
c	call get_job_name
c	call lineardp
c	goto 10

	write(*,*) 'Number of upspin electrons'
	read(*,*) nup
	write(*,*) 'Number of downspin electrons'
	read(*,*) ndn
	nelec=nup+ndn
	write(*,*) '2*Spin ='
	read(*,*) ns2
	if (nup .lt. ndn) then
	 write(*,*) 'N_up<N_dn?'
	 stop
	endif

	call configs(nup,ndn,ltot)
	call Spin_eigenfunctions(nup,ndn,ns2)

	rewind(18)
	read(18,*) ncsf
	write(15, '(i5,6x,"ncsf")') ncsf

	rewind(19)
	read(19,*) nconf
	do i=0,MAXN
	 morbmax(i)=0
	enddo
	nmax=0	
	do iconf=1,nconf
	 if (ibfield .eq. 1) then
	  do j=1,nelec
	  read(19,*) norb(j), morb(j)
	  if(norb(j).gt.MAXN .or. abs(morb(j)).gt.MAXM) then
	   write(*,*) 'Quantum number(s) too large!'
	   stop
	  endif
	  if(abs(morb(j)) .gt. morbmax(norb(j))) morbmax(norb(j))=abs(morb(j))
	  if(norb(j).gt.nmax) nmax=norb(j)
	  enddo
	 elseif (ibfield .eq. 2) then
	  read(19,*) (morb(i),i=1,nelec)
	  do j=1,nelec
	  if(abs(morb(j)) .gt. morbmax(1)) morbmax(1)=abs(morb(j))
	  enddo
	 endif
	 read(19,*)
	enddo


	do i=1,ncsf
      call read_input

c update the sign of coefficents for 3D system
      call sign_3D

c define an ordered list of single electron states (orbitals)
      call def_ses

c the big expansion:
      call expand_det

c order the orbitals within each determinant, and get rid of
c null determinants due to Pauli pr.
      call order_det

c get rid of doubled determinants:
      call reduce_det

      call write_output(morbmax)

	enddo

	rewind(15)
	read(15,*) ncsf
	do i=1,ncsf
	call write_detcsf
	enddo

	write(*,*) 'check the linear dependance among truncated CSFs'
5	call lineardp

	write(*,*) 'choose type of basis to use: 1 for real lda orbitals'
	write(*,*) '                             5 for real part of Fork-Darwin state'
	read(*,*) ibasis

c	if (ibasis.eq.5) then
	 write(*,*) 'lambda?'
	 read(*,*) lambda
	 w0=1.0/lambda**2
	 etry=1.0
	 call write_vmc1(nelec,nup,w0,ns2,ltot,etry)
	 call write_vmc(nelec,ncsf,norb,nmax,morbmax)
	 call write_vmc2(ncsf,norb)
c	endif

	close(15)
	close(16)
	close(17)
	close(18)
	close(19)
	close(27)
10    stop
      end
*******************************************************************************
	subroutine configs(nup,ndn,ltot)
	implicit real*8(a-h,o-z)
	logical badconf
	integer morb(500,20)
      character*35 dstrfile
	common/orb_conf/ nele,npair,nconf,morb,badconf
	common/bfield/  ibfield, dstrfile

	
	nele=nup+ndn

	npair=min(nup,ndn)
	lmdd=nup*(nup-1)/2+ndn*(ndn-1)/2

	write(*,*) 'choose 1. zero B field '
	write(*,*) '       2. mid/large B field'
	read(*,*) ibfield
	if (ibfield .eq. 1) then
	 call b0configs(nup, ndn, ltot)
	else
	 write(*,*) 'input L total'
	 read(*,*) ltot
	 if (ltot .ge. lmdd) then
	  call get_job_name
	  nconf=1
	  call mddconfigs(nup, ndn, ltot)
	  nconf=nconf-1
	  write(19,'(i3,6x,"nconf")') nconf
	  do j=1,nconf
	   write(19,'(<nele>i3)') (morb(nconf+1-j,i), i=1,nele)
	  enddo
	 else
	  write(*,*) 'only work for L >=', lmdd
	  stop
	 endif
	endif




	return
	end
*******************************************************************************
	RECURSIVE subroutine mddconfigs(nu, nd, ltot)
	implicit real*8(a-h,o-z)
	logical badconf
	integer morb(500,20)

	common/orb_conf/ nele,npair,nconf,morb,badconf
c For mid/high B field, where the GS conifguration populates only in LLL(lowest
c landau level with non-negative l values so that total angular momemtum
c L>=L0=nup*(nup-1)/2+ndn*(ndn-1)/2
c where for L=L0, electrons occupy (n,l)=(0,0), (0,1)...(0,nup-1) for up-spin
c electrons and (0,0), (0,1)...(0,ndn-1) for down-spin electrons
c First find out the largest l value possible for electrons lmax=nup+1+L-L0
c deduct lmax out from the total L, decuce to the problem of L-lmax for n-1 QD,
c iterate down to QD with 1 electron. Then lower to lmax-1, repeat the iteration.
	lmdd=nu*(nu-1)/2+nd*(nd-1)/2
	nel=nu+nd
	If (nel .eq. 1) then
	 morb(nconf,1)=ltot
	 badconf=.false.
	 if(morb(nconf,1).gt.morb(nconf,2)) then
	  badconf=.true.
	 elseif(morb(nconf,1).eq.morb(nconf,2).and.morb(nconf,1).eq.morb(nconf,3)) then
        badconf=.true.
	 endif
	 if(badconf) goto 11
	 np=0
	 do i=1,nele-1
	  if(morb(nconf,i).eq.morb(nconf,i+1)) np=np+1
	 enddo
	 if(np.gt.npair) badconf=.true.
	 if(badconf) goto 11
c	 do i=nele,1,-1
c	  if(morb(nconf,i) .eq. 0 .and. nconf .ne. 1) then
c	   morb(nconf,i)=morb(nconf-1,i)
c	  else
c	   goto 10
c	  endif
c	 enddo	
	
10	 write(*,*) 'nconf=',nconf
	 write(*,*) (morb(nconf,i), i=1,nele)
	 nconf=nconf+1
11	 return
	endif

	if(nel .eq. nele) then
	 lmax=ltot-lmdd+nu-1
	elseif(nel .eq. nele-1) then
c to order orbitals in non-decreasing l value
	 lmax=min(ltot-lmdd+nu-1,morb(nconf,nel+1))
	 else
c to eliminate 3 consecutive same l
	  lmax=min(ltot-lmdd+nu-1,morb(nconf,nel+1),morb(nconf,nel+2)-1)
	endif
	
	do i=lmax, nu-1, -1
c	 if((morb(nconf,nel+1) .ne. 0) .and. (i .gt. morb(nconf,nel+1))) then
c	  nconf=nconf-1
c	  return
c	 endif
	 morb(nconf,nel)=i
c	 if(nconf.ne.1 .and. morb(nconf,nel+1).eq.0 .and. morb(nconf-1,nel+1).ne.0) morb(nconf,nel+1)=morb(nconf-1,nel+1)
	 do j=nele,nel+1,-1
	  if(morb(nconf,j) .eq. 0 .and. nconf .ne. 1) then
	   morb(nconf,j)=morb(nconf-1,j)
	  else
	   goto 20
	  endif
	 enddo	
20	 if (nu .eq. nd) then
	  call mddconfigs(nu, nd-1, ltot-i)
	 else
	  call mddconfigs(nu-1, nd, ltot-i)
	 endif
	enddo

	return
	end
*******************************************************************************
	subroutine lineardp

	parameter(maxd=500)	
	real dcoeff(2000,200), dcoeffdp(maxd,9999)
	integer ndett, ncsft
	dimension icsfdp(0:maxd),irdet(20),iidet(20),coeff(200)
	integer ndetcsf(2000), ndetcsfdp(1000), idet(2000,200), idetdp(1000,200)
	logical icsfidp(2000)
	real*8 detdp, epsdp, dpmtrx(maxd,maxd)

	epsdp=1.d-6


	rewind(16)
	rewind(17)

	read(16,*) ndett
c	close(16)

	read(17,*) ncsft
	read(17,*) (ndetcsf(i),i=1,ncsft)

	do i=1,ncsft
	 read(17,*) (idet(i,j),j=1,ndetcsf(i))
	 read(17,*) (dcoeff(i,j),j=1,ndetcsf(i))
	enddo

	do i=1,ncsft
	 icsfidp(i)=.false.
	enddo

c find candidates of dependant CSFs
	call findindp(ndett, ncsft, ndetcsf, idet, icsfidp)

	ncsfdp=ncsft
	write (*,*) ncsfdp
	do i=1,ncsft
	 if(icsfidp(i)) ncsfdp=ncsfdp-1
	enddo
	write (*,*) ncsfdp

	open(27,file='csfdp',status='UNKNOWN',access='sequential',form='formatted')
	write(27,'(i4,6x,"number of linear dp candidates")') ncsfdp
	idpcount=0
	do i=1,ncsft
	if(.not.icsfidp(i)) then
	 idpcount=idpcount+1
	 write(27,'(i3,":",i4,3x,"CSF_dp:CSF")') idpcount,i
	 icsfdp(idpcount)=i 	
	 write(27,'(200i5)') (idet(i,j),j=1,ndetcsf(i))
	 write(27,'(200f12.6)') (dcoeff(i,j),j=1,ndetcsf(i))
	 do id=1,ndett
	  dcoeffdp(idpcount,id)=0.0
	  do idd=1,ndetcsf(i)
	   if (idet(i,idd) .eq. id) dcoeffdp(idpcount,id)=dcoeff(i,idd)
	  enddo
	 enddo
c	 write(27,'(500f12.6)') (dcoeffdp(idpcount,j),j=1,ndett)
	endif
	enddo
	close(27)

c normalize coefficients, to avoid false non-zero determinant.
	do i=1,ncsfdp
	 dnorm=0.0
	 do id=1,ndett
	  dnorm=dnorm+dcoeffdp(i,id)**2
	 enddo
	 dnorm=sqrt(dnorm)
	 do id=1,ndett
	  dcoeffdp(i,id)=dcoeffdp(i,id)/dnorm
	 enddo
	enddo
	
	do i=1,ncsfdp
	 do j=i,ncsfdp
	  dpmtrx(i,j)=0.0
	  do k=1,ndett
	   dpmtrx(i,j)=dpmtrx(i,j)+dcoeffdp(i,k)*dcoeffdp(j,k)
	  enddo
	 enddo
	enddo

	do i=1,ncsfdp
	 do j=1,i-1
	  dpmtrx(i,j)=dpmtrx(j,i)
	 enddo
	enddo

	write(*,*) 'possible linear dependent CSF #', ncsfdp
	if(ncsfdp.eq.0) goto 110
c	open(28,file='dpmtrx',status='UNKNOWN',access='sequential',form='formatted')
c	do i=1,ncsfdp
c	 write(28,'((500f12.6))') (dpmtrx(i,j),j=1,ncsfdp)
c	enddo
c	close(28)
c check if the determinant is non-zero, hence no linear dependency
	call Calculate_deter(ncsfdp,dpmtrx,epsdp,detdp)
	if(detdp.gt.epsdp*100) then
	 write(*,*) 'No linear dependency'
	 goto 110
	endif

c	pause
c find out linear dependency and rewrite .det and .csf files
	call dprewrite(ncsfdp,icsfdp,dpmtrx,epsdp,detdp)
110	return
	end
*******************************************************************************
	subroutine dprewrite(ncsfdp,icsfdp,dpmtrx,epsdp,detdp)

	parameter(maxd=500)	
	real dcoeff(2000,200),coeff(200)
	dimension irdet(20),iidet(20),icsfdp(0:maxd), ndetcsf(2000), idet(2000,200)
	real*8 detdp, epsdp, dpmtrx(maxd,maxd)

c idp: the order of ld(linear dependency) candidates
c icsfdp(idp): the actual CSF

	idp=2
206	call Calculate_deter(idp,dpmtrx,epsdp,detdp)
	
	if(detdp.lt.epsdp*100) then
		rewind(15)
		read(15,*) ncsf
		do icsf=1,icsfdp(idp)-1
		 read(15,*) ne
		 read(15,*)
		 read(15,*)
		 read(15,*) nrdet
		 do i=1,nrdet
		  read(15,*) (irdet(j),j=1,ne)
		 enddo
		 if (nrdet .ne. 0) read(15,*) (coeff(i),i=1,nrdet)
		 read(15,*)
		 read(15,*)
		 read(15,*) nidet
		 do i=1,nidet
		  read(15,*) (iidet(j),j=1,ne)
		 enddo
		 if (nidet .ne. 0) read(15,*) (coeff(i),i=1,nidet)
		enddo
		write(*,*) 'adding more determinants'	
		call write_detcsf

		rewind(17)
		read(17,*) ncsft
		read(17,*) (ndetcsf(i),i=1,ncsft)
		do i=1,ncsft
		 read(17,*) (idet(i,j),j=1,ndetcsf(i))
		 read(17,*) (dcoeff(i,j),j=1,ndetcsf(i))
		enddo
c the new csf is at the end of .csf file, use it to replace the old one with ld.
		ncsfsum=ncsft-1
		rewind(17)
		write (17,'(i4,3x,"ncsfsum")') ncsfsum
		ndetcsf(icsfdp(idp))=ndetcsf(ncsft)
		write(17,'(<ncsfsum>i4,3x,"(ndet_in_csf(icsf),icsf=1,ncsfsum)")') (ndetcsf(i), i=1,ncsfsum)
	 do icsf=1,icsfdp(idp)-1
	  idig=1
	  if(idet(icsf,1).gt.9) idig=2
	  if(idet(icsf,1).gt.99) idig=3
	  if(idet(icsf,1).gt.999) idig=4
	  write(17,'(i<idig>,500i5)') (idet(icsf,j), j=1,ndetcsf(icsf))
	  idig1=8
	  if(dcoeff(icsf,1).lt.0.0) idig1=9
	  write(17,'(f<idig1>.6,500f12.6)') (dcoeff(icsf,j), j=1,ndetcsf(icsf))
	 enddo
	 inew=ncsft
	 idig=1
	 if(idet(inew,1).gt.9) idig=2
	 if(idet(inew,1).gt.99) idig=3
	 if(idet(inew,1).gt.999) idig=4
	 write(17,'(i<idig>,500i5)') (idet(inew,j), j=1,ndetcsf(inew))
	 idig1=8
	 if(dcoeff(inew,1).lt.0.0) idig1=9
	 write(17,'(f<idig1>.6,500f12.6)') (dcoeff(inew,j), j=1,ndetcsf(inew))
 	 do icsf=icsfdp(idp)+1,ncsfsum
	  idig=1
	  if(idet(icsf,1).gt.9) idig=2
	  if(idet(icsf,1).gt.99) idig=3
	  if(idet(icsf,1).gt.999) idig=4
	  write(17,'(i<idig>,500i5)') (idet(icsf,j), j=1,ndetcsf(icsf))
	  idig1=8
	  if(dcoeff(icsf,1).lt.0.0) idig1=9
	  write(17,'(f<idig1>.6,500f12.6)') (dcoeff(icsf,j), j=1,ndetcsf(icsf))
	 enddo
	 goto 208
	else
	 idp=idp+1
	 if(idp.le.ncsfdp) goto 206
	endif

208	call lineardp


210	return
	end
*******************************************************************************
	RECURSIVE subroutine findindp(ndet, ncsf, ndetcsf, idet, icsfidp)
	integer ndett, ncsft, ncount(9999), ndetcsf(2000), idet(2000,200)
	logical icsfidp(2000), idp
c ncount(i): times det. i appears
c
c pick out CSFs with only repeated determiants
	do i=1,ndet
	 ncount(i)=0
	 do j=1,ncsf
	 if (.not.icsfidp(j)) then
	  do k=1,ndetcsf(j)
	   if (idet(j,k) .eq. i) ncount(i)=ncount(i)+1
	   if (ncount(i) .gt. 1) goto 20
	  enddo
	 endif
	 enddo
20	enddo

	nsingle=0
	do i=1,ndet
	 if (ncount(i).eq.1) nsingle=nsingle+1
	enddo
	if (nsingle.eq.0) goto 30


	do i=1,ndet
	if (ncount(i).eq.1) then
	 do j=1,ncsf
	 if (.not.icsfidp(j)) then
	  do k=1,ndetcsf(j)
	   if (idet(j,k) .eq. i) icsfidp(j)=.true.
	  enddo
	 endif
	 enddo
	endif
	enddo
	
	call findindp(ndet, ncsf, ndetcsf, idet, icsfidp)

30	return
	end
*******************************************************************************

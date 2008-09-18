      subroutine write_vmc1(nel,nup,w0,Is,ltot,etot)
      implicit real*8(a-h,o-z)
      character*20 filename
	character*4 fileext
	common/basis/  ibasis

c	common/ivmc/fileext,tau,delta

c      write(*,*)'ivmc extension?'
c	read(*,*)fileext
c	filename='ivmc.5'
c	filename(6:8)=fileext
	if (ibasis.eq.5) then
       open(unit=800,file='ivmc.5',status='replace',form='formatted')
	elseif (ibasis.eq.1) then
       open(unit=800,file='i555',status='replace',form='formatted')
	endif

	tau=6.0
	delta=25.0

	if (ibasis.eq.5)  write(800,*)5,1
	if(nel.lt.10)then
	write(800,'(''DotN='',i1,''l='',i1,''Is='',i1''        title'')'),nel,ltot,is
	else
	write(800,'(''DotN='',i2,''l='',i1,''Is='',i1''        title'')'),nel,ltot,is
      endif
      write(800,'(''1837465927472523                 irn'')')
      write(800,'(''0'',i2,''                              iperiodic,ibasis'')') ibasis
      write(800,'(''0.5''  f8.4''       Hartrees     hb,etrial,eunit'')'),etot
      write(800,'(''100  100  1 100   100           nstep,nblk,nblkeq,nconf,nconf_new,nconf_new'')')
      write(800,'(''0    0    1    -1                idump,irstar,isite,ipr '')')
      write(800,'(i2, f6.2, f4.1, f4.1, f4.1 ''          imetro delta,deltar,deltat,fbias '')')1,delta,4.,1.,1.
      write(800,'(''2 1 1 1 1 0 0 0 0                idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e'')')
      write(800,'(i3, 3x, f6.3''                          nfprod,tau'')'),50,tau
      write(800,'(''-1  0   1  0                      nloc,numr,nforce,nefp'')')
c      write(800,'(''.25 0. 0.                        w0,bext,glande'')')
      write(800,'(f12.8,f6.2,f6.2''  w0,bext,glande'')'),w0,0.,0.


c      write(800,'(''5 3                              nelec,nup'')')
      write(800,'(i2,i2''     nelec,nup'')'),nel,nup
	write(800,*)
      write(800,'(''* Geometry section'')')
      write(800,'(''2                                ndim'')')
      write(800,'(''1 1                              nctype,ncent'')')
      write(800,'(''1                                (iwctype(i),i=1,ncent)'')')
      write(800,'(f4.1,''                               (znuc(i),i=1,nctype)'')'),float(nel)
      write(800,'(''0. 0. 0.                         ((cent(ii,i),ii=1,3),i=1,ncent)'')')
	write(800,*)
      write(800,'(''* Determinantal section'')')
      write(800,'(''0                                inum_orb'')')


      return
	end
***************************************************************
***************************************************************
      subroutine write_vmc2(nlsconf,norb)
      implicit real*8(a-h,o-z)
      real*8 A(6),B(6),C(23)
	common/readinput/nreadparam
	character*4 fileext
c	common/ivmc/fileext,tau

	character*20 fmt

      write(800,*)
      write(800,'(''* Jastrow section'')')
      write(800,'(''1             ianalyt_lap'')')
      write(800,'(''4 4 1 1 5 0   ijas,isc, nspin1,nspin2, nord, ifock'')')
      write(800,'(''5 5 5         norda,nordb,nordc'')')


      if(nreadparam.eq.1)then

	open (unit=55,status='old',form='formatted',file='parameters.for')      	
	read(55,*)scalek,a21
      read(55,*)(a(i),i=1,6)
	read(55,*)(b(i),i=1,6)
	read(55,*)(c(i),i=1,23)
      close(55)
      write(800,'(2(E14.8,2x) '' scalek,a21'')'),scalek,a21
      write(800,'(6(E14.8,2x) '' (a(iparmj),iparmj=1,nparma)'')'),(a(i),i=1,6)
      write(800,'(6(E14.8,2x)'' (b(iparmj),iparmj=1,nparmb)'')'),(b(i),i=1,6)
      write(800,'(23(E14.8,2x) '' (c(iparmj),iparmj=1,nparmc)'')'),(c(i),i=1,23)

      else

      write(800,'(''0.1 0. scalek,a21'')')
      write(800,'(''0. 0. 0. 0. 0. 0. (a(iparmj),iparmj=1,nparma)'')')
      write(800,'(''1. 0. 0. 0. 0. 0. (b(iparmj),iparmj=1,nparmb)'')')
      write(800,'(''0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. (c(iparmj),iparmj=1,nparmc)'')')
      endif
      write(800,*)
      write(800,'(''* Optional features'')')
      write(800,'(''&opt_list '')')
      write(800,'(''/'')')
      write(800,*)

      write(fmt,'(''(''i3,''i3)'')')norb

      write(800,'(''* Optimization section'')')
      write(800,'(''4 1000 1.d-8 0. 1.d-6     nopt_iter,nblk_max,add_diag(1),p_var,tol_energy'')')
      write(800,'(''2000 '',i3, '' -1 1 6 1000 21101 1 NDATA, NPARM, icusp,icusp2, NSIG,NCALLS, iopt,ipr'')') nlsconf+24
      write(800,'(''0 0 0 0 i3body,irewgt,iaver,istrech'')')
      write(800,'(''0 0 0 0 0 0 0 0 0 0 ipos,idcds,idcdr,idcdt,id2cds,id2cdr,id2cdt,idbds,idbdr,idbdt'')')

      write(800, fmt),(-1,i=1,norb)

      write(800,'(''0  4 5 15 0  '',i3,'' 1 0  nparml, nparma,nparmb,nparmc,nparmf, nparmcsf,nparms,nparmg'')') nlsconf-1
      write(800,'(''(iworb(iparm),iwbasi(iparm),iparm=1,nlarml)'')')
      write(800,'(''    (iwbase(iparm),iparm=1,nparm-nparml)'')')
	if (nlsconf.gt.1) then
      write(800,'(<nlsconf-1>i4, '' (iwcsf(iparm),iparm=1,nparmcsf)'')') (i,i=2,nlsconf)
	else
      write(800,'('' (iwcsf(iparm),iparm=1,nparmcsf)'')')
	endif
      write(800,'(''    3 4 5 6 (iwjasa(iparm),iparm=1,nparma)'')')
      write(800,'(''2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)'')')
      write(800,'(''    3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 (iwjasc(iparm),iparm=1,nparmc)'')')
      write(800,'(''0 0        necn,nebase'')')
      write(800,'(''           ((ieorb(j,i),iebasi(j,i),j=1,2),i=1,necn)'')')
      write(800,'(''           ((iebase(j,i),j=1,2),i=1,nebase)'')')

      write(800,fmt),(0,i=1,norb)

      write(800,'(''7.080   eave'')')
      write(800,'(''1.d-15 5. 1 20 4 pmarquardt,tau,noutput,nstep,ibold'')')
      write(800,'(''T F analytic,cholesky'')')

      close(800)
      return
	end
***************************************************************
      subroutine write_vmc(ne,ncsf,norb,nmax,mmax)
	integer mmax(0:20),iwrwf(0:20,50),iwdet(8000,20),idet(2000,200),ndetcsf(2000)
	real coe(50),dcoeff(2000,200)
      character*20 fmt
	common/basis/  ibasis

	norb=0
	do i=0,nmax
	 norb=norb+2*mmax(i)+1
	enddo

	rewind(16)
	read(16,*) ndet
	do i=1,ndet
	 read(16,*) (iwdet(i,j),j=1,ne)
	enddo
	
	write(800,'(i4,i4,i4''  ndet,nbasis,norb'')') ndet,norb,norb
	if (ibasis.eq.5) then
       write(800,'(i4,''  Number of ll-levels'')') nmax+1
	 do i=0,nmax
	  write(800,'(100i5)') i, 2*mmax(i)+1, 0, (ll,-ll,ll=1,mmax(i))
	 enddo
	endif

	if (ibasis.eq.1) then
	 do i=0,nmax
	  write(*,'(''for main quantum number n='', i2)') i
	  write(*,'(''we use orbitals with l='',100i5)') 0, (ll,-ll,ll=1,mmax(i))
	  write(*,*) 'find out corresponding orbitals in basis.n file:'
	  read(*,*) (iwrwf(i,ib),ib=1,2*mmax(i)+1)
	 enddo	
       write(800,'(<norb>i5,'' (iwrwf(ib),ib=1,nbasis)'')') ((iwrwf(in,ib),ib=1,2*mmax(in)+1),in=0,nmax)
       write(800,'(<norb>i5,'' (m_bas(i),i=1,nbasis)'')') (0, (ll,-ll,ll=1,mmax(in)),in=0,nmax)
	endif


      write(fmt,'(''(''i3,''f4.1)'')') norb
	do i=1,norb
      do j=1,norb
	   if(i.eq.j)then
	     coe(j)=1.
	   else
	     coe(j)=0.
	   endif
	enddo
      write(800,fmt)(coe(j),j=1,norb)
	enddo
      write(800,fmt)(1.,j=1,norb)

      write(800,'(<ne>i3,'' (iworbd(j,idet),j=1,nelec)'')') (iwdet(1,j),j=1,ne)
	do i=2,ndet
	 write(800,'(<ne>i3)') (iwdet(i,j),j=1,ne)
	enddo

	rewind(17)
	read(17,*) ncsf
	read(17,*) (ndetcsf(i),i=1,ncsf)

	do i=1,ncsf
	 read(17,*) (idet(i,j),j=1,ndetcsf(i))
	 read(17,*) (dcoeff(i,j),j=1,ndetcsf(i))
	enddo

      write(800,'(i4," ncsf")'),ncsf
	write(800,'(<ncsf>f6.3,'' (csf_coef(icsf),icsf=1,ncsf)'')') 1.0, (0.0,i=1,ncsf-1)
	write(800,'(<ncsf>i4,'' (ndet_in_csf(icsf),icsf=1,ncsf)'')') (ndetcsf(i),i=1,ncsf)

	write(800,'(200i5,'' (iwcsf(idet_in_csf,1),idet_in_csf=1,ndet_in_csf(1))'')') (idet(1,j),j=1,ndetcsf(1))
	write(800,'(200f12.6,'' (cdet_csf(idet_in_csf,1),idet_in_csf=1,ndet_in_csf(1))'')') (dcoeff(1,j),j=1,ndetcsf(1))
	do i=2,ncsf
	 write(800,'(200i5)') (idet(i,j),j=1,ndetcsf(i))
	 write(800,'(200f12.6)') (dcoeff(i,j),j=1,ndetcsf(i))
	enddo

	return
	end


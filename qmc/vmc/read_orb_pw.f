      subroutine read_orb_pw_real
c Written by Cyrus Umrigar
c Reads in pw basis orbitals that have already been converted to be real.
c Presently not used.

      use all_tools_mod

      implicit real*8(a-h,o-z)


      common /dim/ ndim
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /pworbital/c_rp(NGVECX,MORB_OCC),c_rm(NGVECX,MORB_OCC),c_ip(NGVECX,MORB_OCC)
     &,c_im(NGVECX,MORB_OCC),ngorb(MORB),isortg(NGVECX,MORB),isortk(MKPTS),icmplx
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

      dimension rkvec_tmp(3)

      open(3,file='orbitals_pw')
      read(3,*) nkvec,ngvec
      if(nkvec.gt.MKPTS) stop 'nkvec>MKPTS in read_orb_pw'

      jorb=0
      jorba=0
      do 50 ikv=1,nkvec
        read(3,*) ikvec,nband(ikv),(rkvec_tmp(k),k=1,ndim)
        if(ikvec.ne.ikv) stop 'ikvec.ne.ikv in read_orb_pw_real'
        do 10 k=1,ndim
          if(abs(rkvec_tmp(k)-rkvec(k,ikv)).gt.1.d-5) then
            write(6,'(''rkvec_tmp!=rkvec in read_orb_pw_tm.  Since kvecs are now generated in standard order indep of cutg_sim_big''
     &      ,'' likely reasons are that orb_pw_tm needs to be regenerated with kvecs in new order or the k-vec. shift is not right''
     &      )')
            stop 'rkvec_tmp!=rkvec in read_orb_pw_tm.  Likely reasons are orb_pw_tm needs to be regenerated with kvecs in new order
     & or k-vec shift is wrong'
          endif
   10   continue
        do 50 iband=1,nband(ikv)
          jorb=jorb+1
          jorba=jorba+k_inv(ikv)
          read(3,*) ib,eig
          if(ib.ne.iband) stop 'ib.ne.iband in read_orb_pw_real'
          write(6,'(''iband,eig='',i3,f9.5)') ib,eig
   50     read(3,*) (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec)

c Note that jorba can be > nord if ndet>1
      write(6,'(i3,'' orbitals read in from file orbitals_pw'')') jorb
      if(jorba.lt.norb) then
        write(6,'(''jorba,norb='',2i5)') jorba,norb
        stop 'jorba < norb in read_orb_pw'
      endif
      close(3)

      return
      end
c-----------------------------------------------------------------------

      subroutine read_orb_pw
c Written by Cyrus Umrigar
c Reads in pw basis orbitals and convert them to real ones suitable for qmc.
c Warning: At present NGVECX is used to dimension not only quantities that are
c used for the Ewald sum, but also quantities used for the wavefunction pw coefficients.
c There is no reason for the latter to be tied to the former.

c I write out orbitals_pw_lagrange so that in future runs I can just read that file in
c rather than orbitals_pw_tm, but at the moment I am not using that feature.
c Also, I first write out a temporary fort.3 and then delete it just because
c it is only after one has processed all the k-pts that one knows how big ngvec_orb is.
c However, that causes problems when running with mpi, so comment out that part.

      use bwfdet_mod
      use all_tools_mod
      implicit real*8(a-h,o-z)
      character*20 fmt
      character*16 mode,iorb_format

      parameter(eps=1.d-3)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode
      common /contr_names/ iorb_format
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
      common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,ngrid_orbx,ngrid_orby,ngrid_orbz
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange
c    &,orb_splines(8,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:
c    &MGRID_ORB_PER-1,MORB_OCC)
c    &,grid_orbx(0:MGRID_ORB_PER-1),grid_orby(
c    &0:MGRID_ORB_PER-1),grid_orbz(0:MGRID_ORB_PER-1)
c    &,orb_splines_explicit(4,4,4,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:
c    &MGRID_ORB_PER-1,MORB_OCC)
c atom included just to be able to shift test pt. for checking orbs by cent(k,1).
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /tempor_test/ igvec_dft(3,NGVEC_BIGX),iwgvec(NGVEC2X),c_real(NGVEC2X),c_imag(NGVEC2X)
     &,map_gvecdft_gvec(NGVEC2X),isign_gvecdft_gvec(NGVEC2X)
     &,orb(MORB),dorb(3,MORB),ddorb(MORB)
     &,orb_si(MORB),dorb_si(3,MORB),ddorb_si(MORB),iflag(MORB),rnorm,r(3)
     &,rkvec_tmp(3),rkvec_tmp2(3),ngg(MKPTS),ngvec_dft
c    &,orb_splines_tmp(10)
c    &,wk(80*MGRID_ORB_PER*MGRID_ORB_PER*MGRID_ORB_PER)
c    &,orb_tmp(MORB_OCC),dorb_tmp(3,MORB_OCC),ddorb_tmp(MORB_OCC)
c    &,ict(10),ddorb_splines_tmp(3),ibcxmin,ibcymin,ibczmin,ibcxmax,
c    &ibcymax,ibczmax

      dimension ipoint(MORB)
      real*8 r_basis(3),xi,yi,zi
     &      ,bcxmin(0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,MORB),
     &       bcxmax(0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,MORB),
     &       bcymin(0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,MORB),
     &       bcymax(0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,MORB),
     &       bczmin(0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,MORB),
     &       bczmax(0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,MORB)

c     write(6,*)'rlatt_sim=',((rlatt_sim(k,i),i=1,3),k=1,3)
c     stop
      write(6,'(/,''Reading in orbitals for periodic system'',/)')

c Factor for renormalizing orbitals to avoid overflow when calculating determinant.
c If number of electrons is small do not renormalize, so old outputs do not change.
      if(nelec.le.256) then
        rnorm=1.d0
       else
        rnorm=0.1d0
        write(6,'(''Info: orbitals renormalized by 0.1 to prevent overflow'')')
      endif

c Flag the orbitals in orbitals_pw_tm that are actually used
c norb in the input file is >= norb actually used, since some may be skipped.
c Reset norb to the number actually used
c JT: if iorb_used=0, then use all orbitals
      do 4 i=1,norb
        if (iorb_used.eq.0) then
         iflag(i)=1
        else
         iflag(i)=0
        endif
    4  continue
      do 5 idet=1,ndet
        do 5 ielec=1,nelec
    5     iflag(iworbd(ielec,idet))=1
      norb_used=0
      do 6 iorb=1,norb
        if(iflag(iorb).eq.1) then
          norb_used=norb_used+1
          ipoint(iorb)=norb_used
        endif
    6 continue

      if(norb_used.lt.norb) then
        write(6,'(/,''norb(read in), norb(used)='',2i5)') norb,norb_used
        write(6,'(/,''new mapping of orbitals in determinants'')')
        do 8 idet=1,ndet
          do 7 ielec=1,nelec
    7       iworbd(ielec,idet)=ipoint(iworbd(ielec,idet))
          if(nup+ndn.lt.60) then
            write(fmt,'(''('',i2,''i3,3x,'',i2,''i3)'')') nup,ndn
            write(6,fmt) (iworbd(j,idet),j=1,nup),(iworbd(j+nup,idet),j=1,ndn)
           else
            write(6,'(30i4)') (iworbd(j,idet),j=1,nup)
            write(6,'(30i4)') (iworbd(j+nup,idet),j=1,ndn)
          endif
    8   continue
        norb=norb_used
      endif


c Set coordinates of test point
      r(1)=.1d0+cent(1,1)
      r(2)=.2d0+cent(2,1)
      r(3)=.3d0+cent(3,1)

c Read in orbitals on grid if file exists and has correct mesh
      if(inum_orb.ne.0) then
c ** Warning: set ngrid_orbx better
c       ngrid_orbx=ngvec**(1.d0/3.d0)
c       ngrid_orby=ngrid_orbx
c       ngrid_orbz=ngrid_orbx
        if(inum_orb.ne.6) then
           if(inum_orb.eq.4) then
             write(6,'(''Using Lagrange polynomial interpolation'')')
c          elseif(inum_orb.eq.5) then
c            write(6,'(''Using piecewise-polynomial spline interpolation
c    &'')')
           endif
           if(MGRID_ORB_PER.lt.5) stop 'MGRID_ORB_PER < 5'
           if(MGRID_ORB_PER.lt.21) write(6,'(''Warning: MGRID_ORB_PER < 21 ='',i5)') MGRID_ORB_PER
           ngrid_orbx=MGRID_ORB_PER
           ngrid_orby=MGRID_ORB_PER
           ngrid_orbz=MGRID_ORB_PER
           write(6,'(''ngrid_orbx,ngrid_orby,ngrid_orbz='',9i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
        else
	   write(6,'(''Using B-spline approximation'')')
        endif

        num_orb_exist=0

cwparker Read in Lagrange interpolated numerical orbitals from file

        if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
          open(4,file='orbitals_num_lagrange',form='unformatted',
     &status='old',err=10)
          num_orb_exist=1
c         write(6,'(''Lagrange interpolation grid is'',3i5)') ngrid_orbx,ngrid_orby,ngrid_orbz
          read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf
          if(ngrid_orbxf.ne.ngrid_orbx .or. ngrid_orbyf.ne.ngrid_orby
     &       .or.  ngrid_orbzf.ne.ngrid_orbz) then
            write(6,'(''orbital grids do not match, program, orbitals_nu
     &m:'',3i4,x,3i4)')
     &      ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,
     &      ngrid_orbzf
            stop 'orbital grids do not match'
          endif
          read(4) ((((orb_num(iorb,ix,iy,iz),iorb=1,norb),ix=0,
     &             ngrid_orbx-1),iy=0,ngrid_orby-1),iz=0,ngrid_orbz-1)
          read(4) ((((ddorb_num(iorb,ix,iy,iz),iorb=1,norb),ix=0,
     &             ngrid_orbx-1),iy=0,ngrid_orby-1),iz=0,ngrid_orbz-1)
c         read(4) (((((dorb_num(k,iorb,ix,iy,iz),k=1,ndim),iorb=1,norb),
c    &          ix=0,ngrid_orbx-1),iy=0,ngrid_orby-1),iz=0,ngrid_orbz-1)
c I break up the loop over k because otherwise the ifort compiler
c complains that record is too big on reading it back in
          do 9 k=1,3
    9       read(4) ((((dorb_num(k,iorb,ix,iy,iz),iorb=1,norb),ix=0,
     &               ngrid_orbx-1),iy=0,ngrid_orby-1),iz=0,ngrid_orbz-1)

          write(6,'(''Done reading orbitals_num_lagrange'')')
          close(4)
        endif

cwparker Read in spline interpolated numerical orbitals from file

c     if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
c         open(4,file='orbitals_num_splines',form='unformatted',
c    &status='old',err=10)
c         num_orb_exist=1
c         read(4) ngrid_orbxf,ngrid_orbyf,ngrid_orbzf,norbf
c         write(6,*) "read in gridsizes,orbitals=",ngrid_orbxf,
c    &ngrid_orbyf,ngrid_orbzf,norbf
c         if(ngrid_orbxf.ne.ngrid_orbx .or. ngrid_orbyf.ne.ngrid_orby
c    &       .or.  ngrid_orbzf.ne.ngrid_orbz) then
c           write(6,'(''orbital grid sizes do not match, program, orbita
c    &ls_num:'',3i4,x,3i4)')
c    &      ngrid_orbx,ngrid_orby,ngrid_orbz,ngrid_orbxf,ngrid_orbyf,
c    &      ngrid_orbzf
c           stop 'orbital grid sizes do not match'
c         endif
c         if(norbf.ne.norb) then
c           write(6,'(''number of orbitals does not match, program, orbi
c    &tals_num:'',i4,x,i4)')
c    &      norb,norbf
c           stop 'number of orbitals does not match'
c         endif
!cwparker Read in the grids
c         do ix=0,ngrid_orbx-1
c            read(4)grid_orbx(ix)
c         enddo
c         do iy=0,ngrid_orby-1
c            read(4)grid_orby(iy)
c         enddo
c         do iz=0,ngrid_orbz-1
c            read(4)grid_orbz(iz)
c         enddo

!cwparker Read in the values of the grid
c          do ix=0,ngrid_orbx-1
c             do iy=0,ngrid_orby-1
c                do iz= 0,ngrid_orbz-1
c                   do iorb=1,norb
cc Compact splines
c                      do i=1,8
c                         read(4) orb_splines(i,ix,iy,iz,iorb)
c                      enddo
cc Explicit splines
cc                     do i=1,4
cc                        do j=1,4
cc                           do k=1,4
cc                              read(4) orb_splines_explicit(i,j,k,ix,
cc    &                                                      iy,iz,iorb)
cc                           enddo
cc                        enddo
cc                     enddo

c                   enddo
c                enddo
c             enddo
c          enddo

c        endif

cwparker Read in B-spline coefficients from file
        if(inum_orb.eq.6 .or .inum_orb.eq.-6) then
           write(6,*)
           write(6,*)'Reading blip wave function and associated data'
           write(6,*)'=============================================='
cwparker pecent is reset later in the read_input subroutine
           call readbwf(pecent)
           call bwfdet_setup(rnorm)
!     ! WAS
           ngvec_orb = nwvec !! defined in module bwf
           call object_modified ("ngvec_orb")
!     !

           write(6,'(''Done with B-spline setup'')')

        endif

        return
      endif

   10 if(inum_orb.ne.0 .and. num_orb_exist.eq.0)
     &write(6,'(''Warning: orbitals_num_lagrange not read correctly'')')

c Read in orbitals in Jose-Luis Martins or PWSCF format
      if(index(iorb_format,'tm').ne.0) then
        write(6,'(''Reading in orbitals in Jose-Luis Martins format'')')
        call read_orb_pw_tm
       elseif(index(iorb_format,'pwscf').ne.0) then
        write(6,'(''Reading in orbitals in PWSCF format'')')
        call read_orb_pw_pwscf
       else
        stop 'iorb_format should be tm or pwscf'
      endif

      nsum=0
      do 60 i=1,ngnorm
        nsum=nsum+igmult(i)
        if(nsum.ge.ngvec_orb) then
          ngnorm_orb=i
          goto 70
        endif
   60 continue
   70 write(6,'(''ngnorm_orb,ngvec_orb='',9i6)') ngnorm_orb,ngvec_orb

c Test to see if orbitals and derivs. calculated smart way are correct by comparing to dumb way from above
      call my_second(1,'orb_co')
c1    nelec_sav=nelec
c1    nelec=1
c1    call orbitals_pw(r,orb,dorb,ddorb)
      call orbitals_pw_grade(r,orb,dorb,ddorb)
      write(6,'(''orb_complic='',2i3,5d15.8)') (iorb,ireal_imag(iorb),orb(iorb),ddorb(iorb),
     &(dorb(k,iorb),k=1,ndim),iorb=1,norb)
      call my_second(2,'orb_co')

c Check for linear dependencies among orbitals
      write(6,'(/,''Checking for linear dependency among orbitals'',/)')
      istop=0
      do 72 i=2,norb
        do 72 j=1,i-1
          idep=1
          ratio=orb(i)/orb(j)
          if(abs(ratio*(ddorb(j)/ddorb(i))-1).gt.eps) then
            idep=0
          endif
          do 71 k=1,ndim
            if(abs(ratio*dorb(k,j)/dorb(k,i)-1).gt.eps) then
              idep=0
            endif
   71     continue
          if(idep.eq.1) then
            write(6,'(''Warning: orbitals'',i4,'' and'',i4,'' are linearly dependent'')') i,j
            istop=1
          endif
   72   continue

      if(istop.eq.1) stop 'looks like orbitals are linearly dependent'

c Check if orbitals and derivs. calculated by slow and fast method are the same.
      do 73 i=1,norb
        if(abs(orb(i)-orb_si(i)).gt.1.d-10) then
            write(6,'(''i,orb(i),orb_si(i),orb(i)-orb_si(i)='',i5,9d21.14)')
     &      i,orb(i),orb_si(i),orb(i)-orb_si(i)
          stop 'orb(i).ne.orb_si(i) in read_orb_pw'
        endif
        if(abs(ddorb(i)-ddorb_si(i)).gt.1.d-10) then
            write(6,'(''Warning: i,ddorb(i),ddorb_si(i),ddorb(i)-ddorb_si(i)='',i5,9d21.14)')
     &      i,ddorb(i),ddorb_si(i),ddorb(i)-ddorb_si(i)
          if(abs(ddorb(i)-ddorb_si(i)).gt.1.d-3) stop 'ddorb(i).ne.ddorb_si(i) in read_ddorb_pw'
        endif
        do 73 k=1,ndim
          if(abs(dorb(k,i)-dorb_si(k,i)).gt.1.d-10) then
            write(6,'(''k,i,dorb(k,i),dorb_si(k,i),dorb(k,i)-dorb_si(k,i)='',2i5,9d21.14)')
     &      k,i,dorb(k,i),dorb_si(k,i),dorb(k,i)-dorb_si(k,i)
            stop 'dorb(k,i).ne.dorb_si(k,i) in read_orb_pw'
          endif
   73 continue

c1    nelec=nelec_sav

c Calculate orbitals on grid and save them if orbitals_num does not exist
      if(inum_orb.ne.0 .and. num_orb_exist.eq.0) then

c1      nelec_sav=nelec
c1      nelec=1

c       call my_second(0,'begin ')
        call my_second(1,'orb_pw')

cwparker Calculate orbitals at grid points for Lagrange interpolation

        if(inum_orb.eq.4 .or. inum_orb.eq.-4) then

          do 76 ix=0,ngrid_orbx-1
            r_basis(1)=ix/dfloat(ngrid_orbx)
            do 76 iy=0,ngrid_orby-1
              r_basis(2)=iy/dfloat(ngrid_orby)
              do 76 iz=0,ngrid_orbz-1
                r_basis(3)=iz/dfloat(ngrid_orbz)
c               write(6,'(''ix,iy,iz,r_basis'',3i3,9f8.4)') ix,iy,iz,r_basis(1),r_basis(2),r_basis(3)

c Convert to cartesian coodinates
              do 74 k=1,ndim
                r(k)=0
                do 74 i=1,ndim
   74             r(k)=r(k)+rlatt_sim(k,i)*r_basis(i)

              call orbitals_pw_grade(r,orb,dorb,ddorb)
c             write(6,'(''orb_complic='',i3,5f12.8)') (iorb,orb(iorb),ddorb(iorb),(dorb(k,iorb),k=1,ndim),iorb=1,norb)
              do 76 iorb=1,norb
c               write(6,'(''ix'',4i2,9f8.4)')ix,iy,iz,iorb,(r_basis(k),k=1,ndim),(r(k),k=1,ndim),orb(iorb)
                orb_num(iorb,ix,iy,iz)=orb(iorb)
                ddorb_num(iorb,ix,iy,iz)=ddorb(iorb)
                do 76 k=1,ndim
   76             dorb_num(k,iorb,ix,iy,iz)=dorb(k,iorb)

cwparker End if block for Lagrange interpolation setup and printout
        endif

cwparker Calculate orbitals on grid points for spline interpolation

c       if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
c        write(6,*)'Calculating spline coefficients'
!wparker Set up an evenly spaced grid for the control points in each direction
!	 with values from 0 to 1

c         do ix=0,ngrid_orbx-1
c            r_basis(1)=ix/dfloat(ngrid_orbx-1)
c            grid_orbx(ix)=r_basis(1)

c            do iy=0,ngrid_orby-1
c               r_basis(2)=iy/dfloat(ngrid_orby-1)
c               grid_orby(iy)=r_basis(2)

c               do iz=0,ngrid_orbz-1
c                  r_basis(3)=iz/dfloat(ngrid_orbz-1)
c                  grid_orbz(iz)=r_basis(3)

!cwparker Convert the current control point, r_basis(i), into Cartesian
!c        coordinates because the orbitals evaluation routine needs a
!c        point in Cartesian coordinates
!
c                    do k=1,ndim
c                       r(k)=0
c                       do i=1,ndim
c                          r(k)=r(k)+rlatt_sim(k,i)*r_basis(i)
c                       enddo
c                    enddo
!
!cwparker Get the analytic value of the orbitals and their derivatives
!c        at the current control point
!
c                    call orbitals_pw_grade(r,orb,dorb,ddorb)

!cwparker Set the value of the splines function to the values of the
!c        orbitals at the control point
!
c                    do iorb=1,norb
!cwparker Compact splines
c                       orb_splines(1,ix,iy,iz,iorb)=orb(iorb)
!cwparker Explicit splines
!c                      orb_splines_explicit(1,1,1,ix,iy,iz,iorb)=
!c    &                                                         orb(iorb)

!cwparker Set boundary conditions if at a grid point on the boundary
!c		(i.e. ix = 0 or ngrid_orbx-1, etc.)
c                      if (ix.eq.0) then
c                         bcxmin(iy,iz,iorb)=rlatt_sim(1,1)*dorb(1,iorb)
c    &                                      +rlatt_sim(1,2)*dorb(2,iorb)
c    &                                      +rlatt_sim(1,3)*dorb(3,iorb)
c                      endif
c                      if (iy.eq.0) then
c                         bcymin(iz,ix,iorb)=rlatt_sim(2,1)*dorb(1,iorb)
c    &                                      +rlatt_sim(2,2)*dorb(2,iorb)
c    &                                      +rlatt_sim(2,3)*dorb(3,iorb)
c                      endif
c                      if (iz.eq.0) then
c                         bczmin(ix,iy,iorb)=rlatt_sim(3,1)*dorb(1,iorb)
c    &                                      +rlatt_sim(3,2)*dorb(2,iorb)
c    &                                      +rlatt_sim(3,3)*dorb(3,iorb)
c                      endif
c                      if (ix.eq.ngrid_orbx-1) then
c                         bcxmax(iy,iz,iorb)=rlatt_sim(1,1)*dorb(1,iorb)
c    &                                      +rlatt_sim(1,2)*dorb(2,iorb)
c    &                                      +rlatt_sim(1,3)*dorb(3,iorb)
c                      endif
c                      if (iy.eq.ngrid_orby-1) then
c                         bcymax(iz,ix,iorb)=rlatt_sim(2,1)*dorb(1,iorb)
c    &                                      +rlatt_sim(2,2)*dorb(2,iorb)
c    &                                      +rlatt_sim(2,3)*dorb(3,iorb)
c                      endif
c                      if (iz.eq.ngrid_orbz-1) then
c                         bczmax(ix,iy,iorb)=rlatt_sim(3,1)*dorb(1,iorb)
c    &                                      +rlatt_sim(3,2)*dorb(2,iorb)
c    &                                      +rlatt_sim(3,3)*dorb(3,iorb)
c                      endif
c                   enddo

!cwparker End do loops for the grid of control points for the spline
c               enddo
c            enddo
c         enddo
!
!cwparker Output x-boundary condition arrays in first orbital (for testing)
!c        iorb=1
!c        do iy=0,ngrid_orby-1
!c           do iz=0,ngrid_orbz-1
!c              write(6,*)'bcxmin(',iy,',',iz,')=',bcxmin(iy,iz,iorb)
!c           enddo
!c        enddo
!c        write(6,*)
!c        do iy=0,ngrid_orby-1
!c           do iz=0,ngrid_orbz-1
!c              write(6,*)'bcxmax(',iy,',',iz,')=',bcxmax(iy,iz,iorb)
!c           enddo
!c        enddo
!c        stop 'Look at x-boundary condition arrays above'
!
!cwparker End if block for spline grid setup
c        endif

        if(inum_orb.eq.6 .or. inum_orb.eq.-6) then
          stop 'Blip data must already exist in bwfn.data.'
        endif

        call my_second(1,'orb_nu')

cwparker Set up the boundary conditions for the splines
c        ibcimin = -1 (i=x,y,z) are for periodic boundary conditions
c        ibcimin = 1 are for specified slopes at the boundaries

        if(rkvec_shift(1).eq.0.d0) then
           ibcxmin=-1
           ibcxmax=0
        elseif(rkvec_shift(1).eq.0.5d0) then
           ibcxmin=1
           ibcxmax=1
        endif
        if(rkvec_shift(2).eq.0.d0) then
           ibcymin=-1
           ibcymax=0
        elseif(rkvec_shift(2).eq.0.5d0) then
           ibcymin=1
           ibcymax=1
        endif
        if(rkvec_shift(3).eq.0.d0) then
           ibczmin=-1
           ibczmax=0
        elseif(rkvec_shift(3).eq.0.5d0) then
           ibczmin=1
           ibczmax=1
        endif

cwparker Set up the size of the work space needed for spline calculations
c        This is as per instructions in the PSPLINE library

c Compact representation of splines
c       nwk=80*ngrid_orbx*ngrid_orby*ngrid_orbz
c       do i=1,nwk
c          wk(i)=0.d0
c       enddo
c Explicit representation of splines
c       if(ibczmin.eq.1 .or. ibczmin.eq.2 .or. ibczmax.eq.1 .or.
c    &     ibczmax.eq.2) then
c           nwk=16*ngrid_orbx*ngrid_orby*ngrid_orbz
c       else
c           nwk=20*ngrid_orbx*ngrid_orby+10*max(ngrid_orbx,ngrid_orby,
c    &                                          ngrid_orbz)
c       endif

c       if(inum_orb.eq.5 .or. inum_orb.eq.-5) then

!cwparker Make the splines for each orbital, taking the orb_splines array
!c        made above as the input

c         do iorb=1,norb
!cwparker Call for compact version of splines
c               write(6,*)'Calling r8mktricubw with:'
c               write(6,*)'iorb=',iorb
c               write(6,*)'grid_orbx(0),grid_orbx(n-1)=',
c    &grid_orbx(0),grid_orbx(ngrid_orbx-1)
c               write(6,*)'grid_orby(0),grid_orby(n-1)=',
c    &grid_orby(0),grid_orby(ngrid_orby-1)
c               write(6,*)'grid_orbz(0),grid_orbz(n-1)=',
c    &grid_orbz(0),grid_orbz(ngrid_orbz-1)
c               write(6,*)'ngrid_orbx,ngrid_orby,ngrid_orbz=',ngrid_orbx
c    &,ngrid_orby,ngrid_orbz
c               write(6,*)'orb_splines(1,0,0,0,',iorb,')=',
c    &orb_splines(1,0,0,0,iorb)
c               write(6,*)'orb_splines(1,n-1,0,0,',iorb,')=',
c    &orb_splines(1,ngrid_orbx-1,0,0,iorb)
c               write(6,*)'orb_splines(1,0,n-1,0,',iorb,')=',
c    &orb_splines(1,0,ngrid_orby-1,0,iorb)
c               write(6,*)'orb_splines(1,0,0,n-1,',iorb,')=',
c    &orb_splines(1,0,0,ngrid_orbz-1,iorb)
c               write(6,*)'ibcxmin,ibcymin,ibczmin=',ibcxmin,ibcymin,
c    &ibczmin
c               if(ibcxmin.eq.1) then
c                  write(6,*)'bcxmin(0,0,',iorb,'),bcxmin(0,n-1,',iorb,'
c    &)=',bcxmin(0,0,iorb),bcxmin(0,ngrid_orbx-1,iorb)
c                  write(6,*)'bcxmax(0,0,',iorb,'),bcxmax(0,n-1,',iorb,'
c    &)=',bcxmax(0,0,iorb),bcxmax(0,ngrid_orbx-1,iorb)
c               endif
c               if(ibcymin.eq.1) then
c                  write(6,*)'bcymin(0,0,',iorb,'),bcymin(0,n-1,',iorb,'
c    &)=',bcymin(0,0,iorb),bcymin(0,ngrid_orby-1,iorb)
c                  write(6,*)'bcymax(0,0,',iorb,'),bcymax(0,n-1,',iorb,'
c    &)=',bcymax(0,0,iorb),bcymax(0,ngrid_orby-1,iorb)
c               endif
c               if(ibczmin.eq.1) then
c                  write(6,*)'bczmin(0,0,',iorb,'),bczmin(0,n-1,',iorb,'
c    &)=',bczmin(0,0,iorb),bczmin(0,ngrid_orbz-1,iorb)
c                  write(6,*)'bczmax(0,0,',iorb,'),bczmax(0,n-1,',iorb,'
c    &)=',bczmax(0,0,iorb),bczmax(0,ngrid_orbz-1,iorb)
c               endif
c               write(6,*)'nwk=',nwk

c              call r8mktricubw(grid_orbx,ngrid_orbx,grid_orby,
c    &                         ngrid_orby,grid_orbz,ngrid_orbz,
c    &                         orb_splines(1,0,0,0,iorb),
c    &                         ngrid_orby,ngrid_orbz,
c    &                         ibcxmin,bcxmin(0,0,iorb),
c    &                         ibcxmax,bcxmax(0,0,iorb),ngrid_orbx,
c    &                         ibcymin,bcymin(0,0,iorb),
c    &                         ibcymax,bcymax(0,0,iorb),ngrid_orby,
c    &                         ibczmin,bczmin(0,0,iorb),
c    &                         ibczmax,bczmax(0,0,iorb),ngrid_orbz,
c    &                         wk,nwk,1,1,1,ier)

c              if(ier.ne.0) stop 'error in r8mktricubw'
c               write(6,*)'orb_splines(1,0,0,0,',iorb,')=',
c    &orb_splines(1,0,0,0,iorb)
c               write(6,*)

!cwparker Call for explicit version of splines
!c             call r8tcspline(grid_orbx,ngrid_orbx,grid_orby,
!c    &                        ngrid_orby,grid_orbz,ngrid_orbz,
!c    &                        orb_splines_explicit(1,1,1,0,0,0,iorb),
!c    &                        ngrid_orbx,ngrid_orby,
!c    &                        ibcxmin,bcxmin(0,0,iorb),
!c    &                        ibcxmax,bcxmax(0,0,iorb),ngrid_orbx,
!c    &                        ibcymin,bcymin(0,0,iorb),
!c    &                        ibcymax,bcymax(0,0,iorb),ngrid_orby,
!c    &                        ibczmin,bczmin(0,0,iorb),
!c    &                        ibczmax,bczmax(0,0,iorb),ngrid_orbz,
!c    &                        wk,nwk,1,1,1,ier)
c         enddo

!cwparker End if block for splines
c       endif

c Check values on grid pts.
c       do 78 ix=0,ngrid_orbx-1
c       r_basis(1)=ix/dfloat(ngrid_orbx)
c       do 78 iy=0,ngrid_orby-1
c       r_basis(2)=iy/dfloat(ngrid_orby)
c       do 78 iz=0,ngrid_orbz-1
c       r_basis(3)=iz/dfloat(ngrid_orbz)

ccConvert to cartesian coodinates
c       do 77 k=1,ndim
c         r(k)=0
c         do 77 i=1,ndim
c  77       r(k)=r(k)+rlatt_sim(k,i)*r_basis(i)

cc      call orbitals_pw_grade(r,orb,dorb,ddorb)
cc      write(6,'(''orb_complic2='',i3,5f12.8)') (iorb,orb(iorb),ddorb(iorb),(dorb(k,iorb),k=1,ndim),iorb=1,norb)
c  78   call orbitals_period_num_grade(r,orb,dorb,ddorb)
cc 78   write(6,'(''orb_complic3='',i3,5f12.8)') (iorb,orb(iorb),ddorb(iorb),(dorb(k,iorb),k=1,ndim),iorb=1,norb)

        if(index(mode,'mpi').eq.0 .or. idtask.eq.0) then

c       if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
c         open(4,file='orbitals_num.splines',form='unformatted',
c    &status='new')
c         write(4) ngrid_orbx,ngrid_orby,ngrid_orbz,norb

!cwparker Write out the grids
c         do ix=0,ngrid_orbx-1
c            write(4)grid_orbx(ix)
c         enddo
c         do iy=0,ngrid_orby-1
c            write(4)grid_orby(iy)
c         enddo
c         do iz=0,ngrid_orbz-1
c            write(4)grid_orbz(iz)
c         enddo
!
!cwparker Write out the grid point values
c         do ix=0,ngrid_orbx-1
c            do iy=0,ngrid_orby-1
c               do iz=0,ngrid_orbz-1
c                  do iorb=1,norb
ccwparker Compact splines
c                     do i=1,8
c                        write(4)orb_splines(i,ix,iy,iz,iorb)
c                     enddo
ccwparker Explicit splines
cc                     do i=1,4
cc                        do j=1,4
cc                           do k=1,4
cc                              write(4) orb_splines_explicit(i,j,k,ix,
cc    &                                                       iy,iz,iorb)
cc                           enddo
cc                        enddo
cc                     enddo

c                  enddo
c               enddo
c            enddo
c         enddo

c         endif

          if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
          open(4,file='orbitals_num_lagrange',form='unformatted',
     &status='new')
cwparker	Unformatted output
          write(4) ngrid_orbx,ngrid_orby,ngrid_orbz
          write(4) ((((orb_num(iorb,ix,iy,iz),iorb=1,norb),ix=0,ngrid_orbx-1),iy=0,ngrid_orby-1),iz=0,ngrid_orbz-1)
          write(4) ((((ddorb_num(iorb,ix,iy,iz),iorb=1,norb),ix=0,ngrid_orbx-1),iy=0,ngrid_orby-1),iz=0,ngrid_orbz-1)

c         write(4) (((((dorb_num(k,iorb,ix,iy,iz),k=1,ndim),iorb=1,norb),ix=0,ngrid_orbx-1),iy=0,ngrid_orby-1),iz=0,ngrid_orbz-1)
c I break up the loop over k because otherwise the ifort compiler complains that record is too big on reading it back in
          do 79 k=1,3
   79       write(4) ((((dorb_num(k,iorb,ix,iy,iz),iorb=1,norb),ix=0,ngrid_orbx-1),iy=0,ngrid_orby-1),iz=0,ngrid_orbz-1)
          endif
          close(4)
        endif
c       write(14,'(51f12.6)') (orb_num(1,ix,1,1),ix=0,ngrid_orbx-1)
c       write(14,'(51f12.6)') (orb_num(2,ix,2,1),ix=0,ngrid_orbx-1)
c       write(14,'(51f12.6)') (orb_num(3,ix,1,1),ix=0,ngrid_orbx-1)
c       write(14,'(51f12.6)') (orb_num(4,ix,2,1),ix=0,ngrid_orbx-1)
c       write(14,'(51f12.6)') (orb_num(1,1,iy,1),iy=0,ngrid_orby-1)
c       write(14,'(51f12.6)') (orb_num(2,2,iy,1),iy=0,ngrid_orby-1)
c       write(14,'(51f12.6)') (orb_num(3,1,iy,1),iy=0,ngrid_orby-1)
c       write(14,'(51f12.6)') (orb_num(4,2,iy,1),iy=0,ngrid_orby-1)

        call my_second(2,'orb_nu')

c1      nelec=nelec_sav

      endif

c Write file for subsequent read-in, though at present I am not using it.
c Comment it out since it causes problems with mpi.
c     open(4,file='orbitals_pw')
c     write(4,'(i1,3i5,'' icmplx,nkvec,ngnorm_orb,ngvec_orb'')') icmplx,nkvec,ngnorm_orb,ngvec_orb
c     rewind 3
c     read(3,*) nkvec_tmp,ngvec_tmp
c     do 80 ikv=1,nkvec
c       read(3,'(2i4,3f9.5)') ikv_tmp,nband(ikv),(rkvec_tmp(k),k=1,ndim)
c       write(4,'(2i4,3f9.5'' ikvec, nband, rkvec(in recip. lat. units)'')') ikv,nband(ikv),(rkvec_tmp(k),k=1,ndim)
c       do 80 iband=1,nband(ikv)
c         read(3,'(2i5,f10.6)') ikv_tmp,iband_tmp,eig
c         write(4,'(2i5,f10.6,'' ikvec, iband, eig (Ha)'')') ikv,iband,eig
c         if(k_inv(ikv).eq.1) then
c           read(3,'(1p2d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),igv=1,ngvec)
c           write(4,'(1p2d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),igv=1,ngvec_orb)
c          else
c           read(3,'(1p4d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec)
c           write(4,'(1p4d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec_orb)
c         endif
c  80     continue

c     close(3,status='delete')
c     close(4)
      return

      end
c-----------------------------------------------------------------------

      subroutine read_orb_pw_tm
c Written by Cyrus Umrigar
c Reads in pw basis orbitals in Jose-Luis Martins' format and convert them to real ones suitable for qmc.
c Warning: At present NGVECX is used to dimension not only quantities that are
c used for the Ewald sum, but also quantities used for the wavefunction pw coefficients.
c There is no reason for the latter to be tied to the former.

c I write out orbitals_pw so that in future runs I can just read that file in
c rather than orbitals_pw_tm, but at the moment I am not using that feature.
c Also, I first write out a temporary fort.3 and then delete it just because
c it is only after one has processed all the k-pts that one knows how big ngvec_orb is.
c However, that causes problems when running with mpi, so comment out that part.

      use all_tools_mod
      use orbitals_mod

      implicit real*8(a-h,o-z)
      parameter(eps=1.d-6)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /pworbital/c_rp(NGVECX,MORB_OCC),c_rm(NGVECX,MORB_OCC),c_ip(NGVECX,MORB_OCC)
     &,c_im(NGVECX,MORB_OCC),ngorb(MORB),isortg(NGVECX,MORB),isortk(MKPTS),icmplx
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

c igvec_dft needs to be dimensioned with NGVEC_BIGX since in file orbitals_pw_tm the
c list of g vectors at the top is longer than what is actually used.
c The other arrays are dimensioned NGVEC2X rather than NGVECX because planewave code does not
c combine coefs. of G and -G, whereas QMC code does.
      common /tempor_test/ igvec_dft(3,NGVEC_BIGX),iwgvec(NGVEC2X),c_real(NGVEC2X),c_imag(NGVEC2X)
     &,map_gvecdft_gvec(NGVEC2X),isign_gvecdft_gvec(NGVEC2X)
     &,orb(MORB),dorb(3,MORB),ddorb(MORB)
     &,orb_si(MORB),dorb_si(3,MORB),ddorb_si(MORB),iflag(MORB),rnorm,r(3)
     &,rkvec_tmp(3),rkvec_tmp2(3),ngg(MKPTS),ngvec_dft

      call file_exist_or_die (file_orbitals_pw_tm_in)
      open(30,file=file_orbitals_pw_tm_in,err=999)
      read(30,*) icmplx,ngvec_dft
      if(ngvec_dft.gt.NGVEC_BIGX) then
        write(6,'(''ngvec_dft,NGVEC_BIGX,NGVECX,NGNORMX'',9i10)') ngvec_dft,NGVEC_BIGX,NGVECX,NGNORMX
        stop 'ngvec_dft > NGVEC_BIGX in read_orb_pw_tm'
      endif
      do 20 i=1,ngvec_dft
   20   read(30,*) (igvec_dft(k,i),k=1,ndim)
c  20   write(6,*) (igvec_dft(k,i),k=1,ndim)

      read(30,*) nkvec_tmp
      if(nkvec_tmp.ne.nkvec) then
        write(6,'(''nkvec_tmp,nkvec='',9i5)') nkvec_tmp,nkvec
        stop 'nkvec_tmp != nkvec in read_orb_pw_tm'
      endif

      if(ipr.ge.1) write(3,'(2i6,'' nkvec,ngvec'')') nkvec,ngvec
      if(ipr.ge.1) write(6,'(2i6,'' nkvec,ngvec'')') nkvec,ngvec
      if(nkvec.gt.MKPTS) stop 'nkvec>MKPTS in read_orb_pw_tm'

      call my_second(1,'orb_si')

c iorb  goes over all orbitals
c iorba goes over all orbs and takes into account that for k.ne.G/2 we get 2 indep.
c       orbs by combining Psi_k with Psi_-k.
c jorb  goes over all occup. orbitals
c jorba goes over all occup. orbitals and takes into account that for k.ne.G/2 we get 2 indep.
c       orbs by combining Psi_k with Psi_-k.
      iorba=0
      jorb=0
      jorba=0
      ngvec_orb=0
      eigmax=-1.d99
      eigmin=1.d99
      do 90 ikv=1,nkvec
        read(30,*) ikvec,nband(ikv),ngvec_dftorb,(rkvec_tmp(k),k=1,ndim)
        call object_modified ('nband')
        if(ikvec.ne.ikv) stop 'ikvec.ne.ikv in read_orb_pw_tm'
        if(ngvec_dftorb.gt.NGVEC2X) stop 'ngvec_dftorb>NGVEC2X in read_orb_pw_tm'
        ngg(ikv)=ngvec_dftorb
        if(icmplx.eq.0) then
          do 21 ig=1,ngvec_dftorb
   21       c_imag(ig)=0
        endif

c       if(ipr.ge.1) write(3,'(2i4,3f9.5'' ikvec, nband, rkvec(in cartesian units)'')') ikv,nband(ikv),(rkvec(k,ikv),k=1,ndim)
        if(ipr.ge.1) write(3,'(2i4,3f9.5'' ikvec, nband, rkvec(in recip. lat. units)'')') ikv,nband(ikv),(rkvec_tmp(k),k=1,ndim)
        write(6,'(''ikvec, nband='',2i4)') ikv,nband(ikv)
        do 22 k=1,ndim
   22     rkvec_tmp2(k)=rkvec_tmp(k)
        do 23 k=1,ndim
          rkvec_tmp(k)=0
          do 23 i=1,ndim
   23       rkvec_tmp(k)=rkvec_tmp(k)+rkvec_tmp2(i)*glatt(k,i)
        write(6,'(/,''rkvec_tmp in recip. lat. vec. units'',9f9.4)') (rkvec_tmp2(k),k=1,ndim)
        write(6,'(''rkvec_tmp in cartesian coodinates'',9f9.4)') (rkvec_tmp(k),k=1,ndim)
        do 24 k=1,ndim
          if(abs(rkvec_tmp(k)-rkvec(k,ikv)).gt.1.d-5) then
            write(6,'(''kvec='',9f9.6)') (kvec(kk,ikv)+rkvec_shift(kk),kk=1,ndim)
            write(6,'(''rkvec_tmp,rkvec='',9f9.6)') (rkvec_tmp(kk),kk=1,ndim),(rkvec(kk,ikv),kk=1,ndim)
            write(6,'(''rkvec_tmp!=rkvec in read_orb_pw_tm.  Since kvecs are now generated in standard order indep of cutg_sim_big''
     &      ,'' Likely reasons are that orb_pw_tm needs to be regenerated with kvecs in new order or the k-vec. shift is not right''
     &      )')
            stop 'rkvec_tmp!=rkvec in read_orb_pw_tm.  Likely reasons are orb_pw_tm needs to be regenerated with kvecs in new order
     & or k-vec shift is wrong'
          endif
   24   continue
c       if(ngvec_dftorb.gt.NGVEC2X) stop 'ngvec_dftorb>NGVEC2X in read_orb_pw_tm'
        read(30,*) (iwgvec(j),j=1,ngvec_dftorb)

c Establish the mapping between DFT and QMC g-vectors.
c Note that the DFT g-vectors for each k-pt are different and are themselves mapped
c to the list of gvectors at beginning of orbitals_pw_tm file.
        ngvec_found=0
        do 35 igv=1,ngvec
          if(igvec(1,igv).eq.0 .and. igvec(2,igv).eq.0 .and. igvec(3,igv).eq.0) then
            isign_min=1
           else
            isign_min=-1
          endif
          do 30 ig=1,ngvec_dftorb
            do 30 isign=1,isign_min,-2
              norm=0
c             write(6,*) (igvec_dft(k,iwgvec(ig)),k=1,ndim),(igvec(k,igv),k=1,ndim)
              do 25 k=1,ndim
   25           norm=norm+(igvec_dft(k,iwgvec(ig))-isign*igvec(k,igv))**2
              if(norm.eq.0) then
                ngvec_found=ngvec_found+1
                ngvec_orb=max(ngvec_orb,igv)
                map_gvecdft_gvec(ig)=igv
                isign_gvecdft_gvec(ig)=isign
c               write(6,'(''igv,ig,isign,ngvec_found'',9i5)') igv,ig,isign,ngvec_found
              endif
   30     continue
c         write(6,'(''igv,jorb,c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb)'',2i5,9d12.4)')
c    & igv,jorb,c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb)
   35   continue

        call object_modified ('ngvec_orb')

        if(ngvec_found.ne.ngvec_dftorb) then
          if(ipr.ge.1) write(3,'(''ngvec_dftorb,ngvec_found='',2i5)') ngvec_dftorb,ngvec_found
          write(6,'(''ngvec_dftorb,ngvec_found='',2i5)') ngvec_dftorb,ngvec_found
          if(ngvec_found.lt.ngvec_dftorb) then
             write(6,'(''Increase cutg to about'',f6.2,'' or more'')') cutg*(dfloat(ngvec_dftorb)/ngvec_found)**(1.d0/3.d0)
             stop 'probably need to increase cutg to generate more g-vectors'
          endif
          stop 'ngvec_found != ngvec_dftorb in read_orb_pw_tm'
        endif

        nband_tmp=0
        do 80 iband=1,nband(ikv)
          iorba=iorba+k_inv(ikv)
          if(iorba.gt.MORB) stop 'iorba>MORB in read_orb_pw_tm'
          read(30,*) ibandx,eig
          if(ibandx.ne.iband) stop 'ibandx.ne.iband in read_orb_pw_tm'
          if(icmplx.ne.0) then
            read(30,*) (c_real(ig),c_imag(ig),ig=1,ngvec_dftorb)
           else
            read(30,*) (c_real(ig),ig=1,ngvec_dftorb)
          endif

c eigmin is the energy of the lowest unocc. orb, eigmax that of the highest occ. one.
c I have replaced lines below with double 'if' lines because if we do not then
c (iorba-1) may be 0 when k_inv(ikv).eq.1 and accessing iflag(0) is illegal.
c         if(k_inv(ikv).eq.2.and.iflag(iorba-1).ne.iflag(iorba)) stop
c    &    'if k_inv=2 then the 2 indep. states should have same occupancy
c    &    (except for some excited states)'
c         if((k_inv(ikv).eq.1.and.iflag(iorba).eq.0) .or.
c    &    (k_inv(ikv).eq.2.and.iflag(iorba-1).eq.0.and.iflag(iorba).eq.0)) then
c           eigmin=min(eigmin,eig)
c           goto 80
c         endif
          if(k_inv(ikv).eq.2) then
            if(iflag(iorba-1).ne.iflag(iorba)) then
              write(6,'(''ikv,k_inv(ikv),iband,iorba,eig='',4i5,9f9.5)') ikv,k_inv(ikv),iband,iorba,eig
              stop 'if k_inv=2 then the 2 indep. states should have same occupancy (except for some excited states).
     & Likely reason is iworbd is incorrect'
            endif
          endif
          if(k_inv(ikv).eq.1.and.iflag(iorba).eq.0) then
            eigmin=min(eigmin,eig)
            goto 80
          endif
          if(k_inv(ikv).eq.2) then
            if(iflag(iorba-1).eq.0.and.iflag(iorba).eq.0) then
              eigmin=min(eigmin,eig)
              goto 80
            endif
          endif

          eigmax=max(eigmax,eig)
          jorb=jorb+1
          jorba=jorba+k_inv(ikv)
          if(jorba.gt.MORB_OCC) stop 'jorba > MORB_OCC in read_orb_pw_tm'
          nband_tmp=nband_tmp+1

c If there is only one linearly indep. state formed from psi_k and psi_-k then
c determine if that state is real or imaginary.  This is possibly not reliable
c so redo it later using values of orbitals at some point.
          if(k_inv(ikv).eq.1) then
            if(rknorm(ikv).eq.0.d0) then
              ig_min=2
             else
              ig_min=1
            endif
            sum=0
            sum_abs=0
            do 40 ig=ig_min,ngvec_dftorb
              sum=sum+c_real(ig)
   40         sum_abs=sum_abs+abs(c_real(ig))
            if(abs(sum/sum_abs).gt.1.d-6) then
              ireal_imag(jorba)=1
             else
              ireal_imag(jorba)=2
            endif
            if(ipr.ge.1) write(6,'(''ikv,iband,ireal_imag,sum,sum_abs='',3i4,9d12.4)')
     &      ikv,iband,ireal_imag(jorba),sum,sum_abs
           else
            ireal_imag(jorba-1)=0
            ireal_imag(jorba)=0
          endif

c Set ngvec of them to 0 because we do not know until we have processed all k-pts
c what the final value of ngvec_orb will be, and ngvec is an upper bound to ngvec_orb.
          do 45 igv=1,ngvec
            c_rp(igv,jorb)=0
            c_rm(igv,jorb)=0
            c_ip(igv,jorb)=0
   45       c_im(igv,jorb)=0

          do 50 ig=1,ngvec_dftorb
            igv=map_gvecdft_gvec(ig)
            isign=isign_gvecdft_gvec(ig)
            c_rp(igv,jorb)=c_rp(igv,jorb)+c_real(ig)
            c_rm(igv,jorb)=c_rm(igv,jorb)+c_real(ig)*isign
            c_ip(igv,jorb)=c_ip(igv,jorb)+c_imag(ig)
   50       c_im(igv,jorb)=c_im(igv,jorb)+c_imag(ig)*isign
c           write(6,'(''igv,jorb,c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb)'',2i5,9d12.4)')
c    & igv,jorb,c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb)

          do 60 igv=1,ngvec
            c_rp(igv,jorb)=rnorm*c_rp(igv,jorb)
            c_rm(igv,jorb)=rnorm*c_rm(igv,jorb)
            c_ip(igv,jorb)=rnorm*c_ip(igv,jorb)
   60       c_im(igv,jorb)=rnorm*c_im(igv,jorb)
c         call my_second(2,'42    ')

          call object_modified ('c_rp')
          call object_modified ('c_rm')
          call object_modified ('c_ip')
          call object_modified ('c_im')

          if(ipr.ge.1) write(3,'(i5,f10.6,'' iband, eig (Ha)'')') iband,eig
          if(ipr.ge.0) write(6,'(i5,f10.6,'' iband, eig (Ha)'')') iband,eig
          if(k_inv(ikv).eq.1) then
            if(ipr.ge.1) write(3,'(1p2d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),igv=1,ngvec)
           else
            if(ipr.ge.1) write(3,'(1p4d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec)
          endif

c We calculate orbitals dumb way and then smart way so that they can be compared
c by eye to make sure smart way is right.
c Calculation of orbitals dumb way with orbitals_pw2 is time consuming.
c Also, choose which state to keep when there is only one independent state (k=G/2)
c     call my_second(1,'orb_s1')
c Test to see if orbitals and derivs. calculated correctly
          call orbitals_pw2(ikv,iband,jorba)
c If there is only 1 independent state, choose it to be the one with the largest absolute value
c unless it is linearly dependent on an already chosen state
          if(k_inv(ikv).eq.1) then

            if(jorba.ge.2) then
              if(orb_si(jorba+1)*ddorb_si(jorba+1).eq.0.d0 .or.
c    &        abs(orb_si(jorba+1)*ddorb_si(jorba-1)/(orb_si(jorba-1)*ddorb_si(jorba+1))-1).lt.eps) then
     &        abs(orb_si(jorba+1)*ddorb_si(jorba-1)-orb_si(jorba-1)*ddorb_si(jorba+1)).lt.
     &        abs(eps*orb_si(jorba-1)*ddorb_si(jorba+1))) then
                ireal_imag(jorba)=1
                goto 70
              endif
              if(orb_si(jorba)*ddorb_si(jorba).eq.0.d0 .or.
c    &        abs(orb_si(jorba)*ddorb_si(jorba-1)/(orb_si(jorba-1)*ddorb_si(jorba))-1).lt.eps) then
     &        abs(orb_si(jorba)*ddorb_si(jorba-1)-orb_si(jorba-1)*ddorb_si(jorba)).lt.
     &        abs(eps*orb_si(jorba-1)*ddorb_si(jorba))) then
                ireal_imag(jorba)=2
                goto 70
              endif
            endif
            if(jorba.ge.3) then
              if(orb_si(jorba+1)*ddorb_si(jorba+1).eq.0.d0 .or.
c    &        abs(orb_si(jorba+1)*ddorb_si(jorba-2)/(orb_si(jorba-2)*ddorb_si(jorba+1))-1).lt.eps) then
     &        abs(orb_si(jorba+1)*ddorb_si(jorba-2)-orb_si(jorba-2)*ddorb_si(jorba+1)).lt.
     &        abs(eps*orb_si(jorba-2)*ddorb_si(jorba+1))) then
                ireal_imag(jorba)=1
                goto 70
              endif
              if(orb_si(jorba)*ddorb_si(jorba).eq.0.d0 .or.
c    &        abs(orb_si(jorba)*ddorb_si(jorba-2)/(orb_si(jorba-2)*ddorb_si(jorba))-1).lt.eps) then
     &        abs(orb_si(jorba)*ddorb_si(jorba-2)-orb_si(jorba-2)*ddorb_si(jorba)).lt.
     &        abs(eps*orb_si(jorba-2)*ddorb_si(jorba))) then
                ireal_imag(jorba)=2
                goto 70
              endif
            endif
            if(abs(orb_si(jorba)/orb_si(jorba+1)).gt.1.d0) then
              ireal_imag(jorba)=1
             else
              ireal_imag(jorba)=2
            endif

   70       if(ireal_imag(jorba).eq.1) then
              write(6,'(''iband,ireal_imag,orb_simple='',i3,i2,5d15.8)')
     &        iband,ireal_imag(jorba),orb_si(jorba),ddorb_si(jorba),(dorb_si(k,jorba),k=1,ndim)
c Warning: tmp
              write(6,'(''iband,           orb_simple='',i3,2x, 5d15.8)')
     &        iband,orb_si(jorba+1),ddorb_si(jorba+1),(dorb_si(k,jorba+1),k=1,ndim)
             else
c Warning: tmp
              write(6,'(''iband,           orb_simple='',i3,2x, 5d15.8)')
     &        iband,orb_si(jorba),ddorb_si(jorba),(dorb_si(k,jorba),k=1,ndim)

              write(6,'(''iband,ireal_imag,orb_simple='',i3,i2,5d15.8)')
     &        iband,ireal_imag(jorba),orb_si(jorba+1),ddorb_si(jorba+1),(dorb_si(k,jorba+1),k=1,ndim)
c Save the orb_si and 2nd deriv in order to check for linear dependency
              orb_si(jorba)=orb_si(jorba+1)
              ddorb_si(jorba)=ddorb_si(jorba+1)
              do  k=1,ndim
                dorb_si(k,jorba)=dorb_si(k,jorba+1)
              enddo
            endif
           else
            ireal_imag(jorba-1)=0
            ireal_imag(jorba)=0
            write(6,'(''iband,ireal_imag,orb_simple='',i3,i2,5d15.8)')
     &      iband,ireal_imag(jorba-1),orb_si(jorba-1),ddorb_si(jorba-1),(dorb_si(k,jorba-1),k=1,ndim)
            write(6,'(''iband,ireal_imag,orb_simple='',i3,i2,5d15.8)')
     &      iband,ireal_imag(jorba),orb_si(jorba),ddorb_si(jorba),(dorb_si(k,jorba),k=1,ndim)
          endif

   80   continue
   90   nband(ikv)=nband_tmp
      call my_second(2,'orb_si')

      write(6,'(2i4,'' orbitals read in from file orbitals_pw_tm'')') jorb,jorba
      if(jorba.lt.norb) then
        write(6,'(''jorba,norb='',2i5)') jorba,norb
        stop 'jorba < norb in read_orb_pw_tm'
      endif

      write(6,'(/,''energy of highest  occupied orbital='',f10.6)') eigmax
      if(eigmin.ne.1.d99) then
        write(6,'(''energy of lowest unoccupied orbital='',f10.6)') eigmin
       else
        write(6,'(''no unoccupied orbitals'')')
      endif
      if(eigmin.lt.eigmax) write(6,'(''Warning: Energy of lowest unoccupied
     &orbital is lower than that of the highest occupied one.'',/,
     &''This may be because you are calculating an excited state, or it may be
     & that you are calculating a ground state but the single particle'',/,
     &''band structure has a band overlap, though in fact there is none.'')')

      close(30)
      return

  999 write(6,'(''Error: file orbitals_pw_tm is missing'')')
      stop 'file orbitals_pw_tm is missing'
      end
c-----------------------------------------------------------------------

      subroutine read_orb_pw_pwscf
c Written by Cyrus Umrigar
c Reads in pw basis orbitals in PWSCF format and convert them to real ones suitable for qmc.
c Warning: At present NGVECX is used to dimension not only quantities that are
c used for the Ewald sum, but also quantities used for the wavefunction pw coefficients.
c There is no reason for the latter to be tied to the former.

c I write out orbitals_pw so that in future runs I can just read that file in
c rather than orbitals_pw_pwscf, but at the moment I am not using that feature.
c Also, I first write out a temporary fort.3 and then delete it just because
c it is only after one has processed all the k-pts that one knows how big ngvec_orb is.
c However, that causes problems when running with mpi, so comment out that part.

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'ewald.h'
      include 'numorb.h'
      parameter(eps=1.d-4)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /pworbital/c_rp(NGVECX,MORB_OCC),c_rm(NGVECX,MORB_OCC),c_ip(NGVECX,MORB_OCC)
     &,c_im(NGVECX,MORB_OCC),ngorb(MORB),isortg(NGVECX,MORB),isortk(MKPTS),icmplx
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

c igvec_dft needs to be dimensioned with NGVEC_BIGX since in file orbitals_pw_pwscf the
c list of g vectors at the top is longer than what is actually used.
c The other arrays are dimensioned NGVEC2X rather than NGVECX because planewave code does not
c combine coefs. of G and -G, whereas QMC code does.
      common /tempor_test/ igvec_dft(3,NGVEC_BIGX),iwgvec(NGVEC2X),c_real(NGVEC2X),c_imag(NGVEC2X)
     &,map_gvecdft_gvec(NGVEC2X),isign_gvecdft_gvec(NGVEC2X)
     &,orb(MORB),dorb(3,MORB),ddorb(MORB)
     &,orb_si(MORB),dorb_si(3,MORB),ddorb_si(MORB),iflag(MORB),rnorm,r(3)
     &,rkvec_tmp(3),rkvec_tmp2(3),ngg(MKPTS),ngvec_dft
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      dimension centx(3),gvec_dft(3),gvec_latt(3),rkvec_latt(3)

      complex*16 c_complex_tmp

c The current version of PWSCF does not exploit inversion symmetry to make the
c PW coefs. real.  So, set icmplx=1.  Dario Alfe tells me that the older versions
c did exploit inversion symmetry.  If they restore this feature then I should
c check if there is inversion symmetry about the origin and if there is then real
c plane wave coefs. are read in, otherwise complex.
      icmplx=1

      open(30,file='orbitals_pw_pwscf',err=999)

      do 12 i=1,34
   12   read(30,*)
      read(30,*) ncentx
      if(ncentx.ne.ncent) stop 'ncentx.ne.ncent in read_orb_pw_casino'
      read(30,*)
      do 14 i=1,ncent
        call systemflush(6)
        read(30,*) iznuc_allelec,(centx(k),k=1,ndim)
        do 14 k=1,ndim
          if(abs(centx(k)-cent(k,i)).gt.1.d-6) then
            write(6,'(''k,i,centx(k),cent(k,i)='',2i5,9f9.5)') k,i,centx(k),cent(k,i)
            stop 'centx(k).ne.cent(k,i) in read_orb_pw_casino'
          endif
   14   continue
      do 16 i=1,8
   16   read(30,*)
      read(30,*) ngvec_dft
      if(ngvec_dft.gt.NGVEC_BIGX) stop 'ngvec_dft > NGVEC_BIGX in read_orb_pw_tm'
c In orbitals_pw_tm ngvec_dftorb depends on k vector and can be smaller than ngvec_dft.
c Here it is the same.
      ngvec_dftorb=ngvec_dft

      read(30,*)
      do 19 igv=1,ngvec_dft
        read(30,*) (gvec_dft(k),k=1,ndim)
c  19   write(6,*) (igvec_dft(k,igv),k=1,ndim)
        do 19 k=1,ndim
          gvec_latt(k)=0
          do 17 i=1,ndim
   17       gvec_latt(k)=gvec_latt(k)+glatt_inv(k,i)*gvec_dft(i)
          if(abs(gvec_latt(k)-nint(gvec_latt(k))).gt.1.d-6) then
            write(6,'(''gvec_latt(k),nint(gvec_latt(k))='',f9.5,i5)') gvec_latt(k),nint(gvec_latt(k))
            stop 'gvec_latt(k)-nint(gvec_latt(k)).gt.1.d-10'
          endif
   19     igvec_dft(k,igv)=nint(gvec_latt(k))


      do 20 i=1,4
   20   read(30,*)

      read(30,*) nkvec_tmp
      write(6,*) nkvec_tmp
      if(nkvec_tmp.ne.nkvec) then
        write(6,'(''nkvec_tmp,nkvec='',9i5)') nkvec_tmp,nkvec
        stop 'nkvec_tmp != nkvec in read_orb_pw_pwscf'
      endif

      if(ipr.ge.1) write(3,'(2i6,'' nkvec,ngvec'')') nkvec,ngvec
      if(ipr.ge.1) write(6,'(2i6,'' nkvec,ngvec'')') nkvec,ngvec
      if(nkvec.gt.MKPTS) stop 'nkvec>MKPTS in read_orb_pw_pwscf'

      call my_second(1,'orb_si')

c iorb  goes over all orbitals
c iorba goes over all orbs and takes into account that for k.ne.G/2 we get 2 indep.
c       orbs by combining Psi_k with Psi_-k.
c jorb  goes over all occup. orbitals
c jorba goes over all occup. orbitals and takes into account that for k.ne.G/2 we get 2 indep.
c       orbs by combining Psi_k with Psi_-k.
      iorba=0
      jorb=0
      jorba=0
      ngvec_orb=0
      eigmax=-1.d99
      eigmin=1.d99
      do 90 ikv=1,nkvec
        read(30,*)
        read(30,*) ikvec,nband(ikv),nband_dn,(rkvec(k,ikv),k=1,ndim)
        if(nband_dn.ne.0) stop 'for the moment nband_dn must be 0'
        if(ikvec.ne.ikv) stop 'ikvec.ne.ikv in read_orb_pw_pwscf'
        ngg(ikv)=ngvec_dftorb
        if(icmplx.eq.0) then
          do 21 ig=1,ngvec_dftorb
   21       c_imag(ig)=0
        endif

c       if(ipr.ge.1) write(3,'(2i4,3f9.5'' ikvec, nband, rkvec(in cartesian units)'')') ikv,nband(ikv),(rkvec(k,ikv),k=1,ndim)
        write(6,'(''ikvec, nband='',2i4)') ikv,nband(ikv)

c Figure out k-vector in primitive-cell lattice vector units
        write(6,'(/,i2,'' k-vectors (shifted) in recip. latt. units, and wts, for input to pw program'')') nkvec
        do 22 k=1,ndim
          rkvec_latt(k)=0
          do 22 i=1,ndim
   22       rkvec_latt(k)=rkvec_latt(k)+glatt_inv(k,i)*rkvec(i,ikv)
        write(6,'(''k-vec('',i2,'')='',3f14.10,f4.0)') ikv,(rkvec_latt(k),k=1,ndim),dfloat(k_inv(ikv))

c If k is half a primitive-cell lattice vector then k_inv=1, else it is 2.
        rknorm(ikv)=0
        k_inv(ikv)=1
        do 23 k=1,ndim
          if(abs(rkvec_latt(k)-.5d0*nint(2*rkvec_latt(k))).gt.1.d-12) k_inv(ikv)=2
   23     rknorm(ikv)=rknorm(ikv)+rkvec(k,ikv)**2
        rknorm(ikv)=sqrt(rknorm(ikv))

c Figure out k-vector in reciprocal simulation-cell lattice vector units with k-shift removed
        write(6,'(/,i2,'' k-vectors (unshifted) in recip. simul. cell latt. units'')') nkvec
        do 25 k=1,ndim
   25     rkvec_tmp(k)=rkvec(k,ikv)-rkvec_shift(k)
        do 28 k=1,ndim
          rkvec_latt(k)=0
          do 26 i=1,ndim
   26       rkvec_latt(k)=rkvec_latt(k)+glatt_sim_inv(k,i)*rkvec_tmp(i)
          if(abs(rkvec_latt(k)-nint(rkvec_latt(k))).gt.1.d-12)
     &    stop 'unshifted k-vector is not a integer in recip. simul. cell latt. units'
   28   kvec(k,ikv)=nint(rkvec_latt(k))
        write(6,'(''k-vec('',i2,'')='',3f14.10)') ikv,(rkvec_latt(k),k=1,ndim)

c I map iwgvec(j)=j in order to keep most of the code the same for orbitals read from TM and other PW codes.
c       read(30,*) (iwgvec(j),j=1,ngvec_dftorb)
        do j=1,ngvec_dftorb
          iwgvec(j)=j
        enddo

c Establish the mapping between DFT and QMC g-vectors.
c Note that the DFT g-vectors for each k-pt are different and are themselves mapped
c to the list of gvectors at beginning of orbitals_pw_tm file but not for orbitals_pw_pwscf
        ngvec_found=0
        do 35 igv=1,ngvec
          if(igvec(1,igv).eq.0 .and. igvec(2,igv).eq.0 .and. igvec(3,igv).eq.0) then
            isign_min=1
           else
            isign_min=-1
          endif
          do 32 ig=1,ngvec_dftorb
            do 32 isign=1,isign_min,-2
              norm=0
c             write(6,*) (igvec_dft(k,iwgvec(ig)),k=1,ndim),(igvec(k,igv),k=1,ndim)
              do 30 k=1,ndim
   30           norm=norm+(igvec_dft(k,iwgvec(ig))-isign*igvec(k,igv))**2

              if(norm.eq.0) then
                ngvec_found=ngvec_found+1
                ngvec_orb=max(ngvec_orb,igv)
                map_gvecdft_gvec(ig)=igv
                isign_gvecdft_gvec(ig)=isign
c               write(6,'(''igv,ig,isign,ngvec_found'',9i5)') igv,ig,isign,ngvec_found
              endif
   32     continue
c         write(6,'(''igv,jorb,c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb)'',2i5,9d12.4)')
c    & igv,jorb,c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb)
   35   continue

        if(ngvec_found.ne.ngvec_dftorb) then
          if(ipr.ge.1) write(3,'(''ngvec_dftorb,ngvec_found='',2i5)') ngvec_dftorb,ngvec_found
          write(6,'(''ngvec_dftorb,ngvec_found='',2i5)') ngvec_dftorb,ngvec_found
          if(ngvec_found.lt.ngvec_dftorb) then
             write(6,'(''Increase cutg to about'',f6.2,'' or more'')') cutg*(dfloat(ngvec_dftorb)/ngvec_found)**(1.d0/3.d0)
             stop 'probably need to increase cutg to generate more g-vectors'
          endif
          stop 'ngvec_found != ngvec_dftorb in read_orb_pw_pwscf'
        endif

        nband_tmp=0
        write(6,*)'ngvec_dftorb=',ngvec_dftorb
        do 80 iband=1,nband(ikv)
          iorba=iorba+k_inv(ikv)
          if(iorba.gt.MORB) stop 'iorba>MORB in read_orb_pw_tm'
          read(30,*)
c         write(6,*) 'Skipped "Band, spin, eigenvalue (au)" line'
          read(30,*) ibandx,ispin,eig
c         write(6,*) ibandx,ispin,eig
          read(30,*)
c         write(6,*) 'Skipped "Eigenvector coeffecients" line'
          if(ibandx.ne.iband) stop 'ibandx.ne.iband in read_orb_pw_pwscf'
          if(ispin.ne.1) stop 'ispin.ne.1 in read_orb_pw_pwscf'
          if(icmplx.ne.0) then
cwparker Orignal statement
c           read(30,*) (c_real(ig),c_imag(ig),ig=1,ngvec_dftorb)
cwparker Modified statement
            do ig=1,ngvec_dftorb
               read(30,*) c_complex_tmp
c              write(6,*) ig,c_complex_tmp
               c_real(ig)=real(c_complex_tmp)
               c_imag(ig)=aimag(c_complex_tmp)
            enddo
           else
            read(30,*) (c_real(ig),ig=1,ngvec_dftorb)
          endif

c eigmin is the energy of the lowest unocc. orb, eigmax that of the highest occ. one.
c I have replaced lines below with double 'if' lines because if we do not then
c (iorba-1) may be 0 when k_inv(ikv).eq.1 and accessing iflag(0) is illegal.
c         if(k_inv(ikv).eq.2.and.iflag(iorba-1).ne.iflag(iorba)) stop
c    &    'if k_inv=2 then the 2 indep. states should have same occupancy
c    &    (except for some excited states)'
c         if((k_inv(ikv).eq.1.and.iflag(iorba).eq.0) .or.
c    &    (k_inv(ikv).eq.2.and.iflag(iorba-1).eq.0.and.iflag(iorba).eq.0)) then
c           eigmin=min(eigmin,eig)
c           goto 80
c         endif
          if(k_inv(ikv).eq.2) then
            if(iflag(iorba-1).ne.iflag(iorba)) then
              write(6,'(''ikv,k_inv(ikv),iband,iorba,eig='',4i5,9f9.5)') ikv,k_inv(ikv),iband,iorba,eig
              stop 'if k_inv=2 then the 2 indep. states should have same occupancy (except for some excited states).
     & Likely reason is iworbd is incorrect'
            endif
          endif
          if(k_inv(ikv).eq.1.and.iflag(iorba).eq.0) then
            eigmin=min(eigmin,eig)
            goto 80
          endif
          if(k_inv(ikv).eq.2) then
            if(iflag(iorba-1).eq.0.and.iflag(iorba).eq.0) then
              eigmin=min(eigmin,eig)
              goto 80
            endif
          endif

          eigmax=max(eigmax,eig)
          jorb=jorb+1
          jorba=jorba+k_inv(ikv)
          if(jorba.gt.MORB_OCC) stop 'jorba > MORB_OCC in read_orb_pw_pwscf'
          nband_tmp=nband_tmp+1

c If there is only one linearly indep. state formed from psi_k and psi_-k then
c determine if that state is real or imaginary.  This is possibly not reliable
c so redo it later using values of orbitals at some point.
          if(k_inv(ikv).eq.1) then
            if(rknorm(ikv).eq.0.d0) then
              ig_min=2
             else
              ig_min=1
            endif
            sum=0
            sum_abs=0
            do 40 ig=ig_min,ngvec_dftorb
              sum=sum+c_real(ig)
   40         sum_abs=sum_abs+abs(c_real(ig))
            if(abs(sum/sum_abs).gt.1.d-6) then
              ireal_imag(jorba)=1
             else
              ireal_imag(jorba)=2
            endif
            if(ipr.ge.1) write(6,'(''ikv,iband,ireal_imag,sum,sum_abs='',3i4,9d12.4)')
     &      ikv,iband,ireal_imag(jorba),sum,sum_abs
           else
            ireal_imag(jorba-1)=0
            ireal_imag(jorba)=0
          endif

c Set ngvec of them to 0 because we do not know until we have processed all k-pts
c what the final value of ngvec_orb will be, and ngvec is an upper bound to ngvec_orb.
          do 45 igv=1,ngvec
            c_rp(igv,jorb)=0
            c_rm(igv,jorb)=0
            c_ip(igv,jorb)=0
   45       c_im(igv,jorb)=0

          do 50 ig=1,ngvec_dftorb
            igv=map_gvecdft_gvec(ig)
            isign=isign_gvecdft_gvec(ig)
            c_rp(igv,jorb)=c_rp(igv,jorb)+c_real(ig)
            c_rm(igv,jorb)=c_rm(igv,jorb)+c_real(ig)*isign
            c_ip(igv,jorb)=c_ip(igv,jorb)+c_imag(ig)
   50       c_im(igv,jorb)=c_im(igv,jorb)+c_imag(ig)*isign
c           write(6,'(''igv,jorb,c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb)'',2i5,9d12.4)')
c    & igv,jorb,c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb)

          do 60 igv=1,ngvec_dft
            c_rp(igv,jorb)=rnorm*c_rp(igv,jorb)
            c_rm(igv,jorb)=rnorm*c_rm(igv,jorb)
            c_ip(igv,jorb)=rnorm*c_ip(igv,jorb)
   60       c_im(igv,jorb)=rnorm*c_im(igv,jorb)
c         call my_second(2,'42    ')

          if(ipr.ge.1) write(3,'(i5,f10.6,'' iband, eig (Ha)'')') iband,eig
          if(ipr.ge.0) write(6,'(i5,f10.6,'' iband, eig (Ha)'')') iband,eig
          if(k_inv(ikv).eq.1) then
            if(ipr.ge.1) write(3,'(1p2d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),igv=1,ngvec)
           else
            if(ipr.ge.1) write(3,'(1p4d22.14)') (c_rp(igv,jorb),c_rm(igv,jorb),c_ip(igv,jorb),c_im(igv,jorb),igv=1,ngvec)
          endif

c We calculate orbitals dumb way and then smart way so that they can be compared
c by eye to make sure smart way is right.
c Calculation of orbitals dumb way with orbitals_pw2 is time consuming.
c Also, choose which state to keep when there is only one independent state (k=G/2)
c     call my_second(1,'orb_s1')
c Test to see if orbitals and derivs. calculated correctly
          call orbitals_pw2(ikv,iband,jorba)
c If there is only 1 independent state, choose it to be the one with the largest absolute value
c unless it is linearly dependent on an already chosen state
          if(k_inv(ikv).eq.1) then
            if(jorba.ge.2) then
              if(abs(orb_si(jorba+1)*ddorb_si(jorba-1)/(orb_si(jorba-1)*ddorb_si(jorba+1))-1).lt.eps) then
                ireal_imag(jorba)=1
                goto 70
              endif
              if(abs(orb_si(jorba)*ddorb_si(jorba-1)/(orb_si(jorba-1)*ddorb_si(jorba))-1).lt.eps) then
                ireal_imag(jorba)=2
                goto 70
              endif
            endif
            if(jorba.ge.3) then
              if(abs(orb_si(jorba+1)*ddorb_si(jorba-2)/(orb_si(jorba-2)*ddorb_si(jorba+1))-1).lt.eps) then
                ireal_imag(jorba)=1
                goto 70
              endif
              if(abs(orb_si(jorba)*ddorb_si(jorba-2)/(orb_si(jorba-2)*ddorb_si(jorba))-1).lt.eps) then
                ireal_imag(jorba)=2
                goto 70
              endif
            endif
            if(abs(orb_si(jorba)/orb_si(jorba+1)).gt.1.d0) then
              ireal_imag(jorba)=1
             else
              ireal_imag(jorba)=2
            endif

   70       if(ireal_imag(jorba).eq.1) then
              write(6,'(''iband,ireal_imag,orb_simple='',i3,i2,5d15.8)')
     &        iband,ireal_imag(jorba),orb_si(jorba),ddorb_si(jorba),(dorb_si(k,jorba),k=1,ndim)
c Warning: tmp
              write(6,'(''iband,           orb_simple='',i3,2x, 5d15.8)')
     &        iband,orb_si(jorba+1),ddorb_si(jorba+1),(dorb_si(k,jorba+1),k=1,ndim)
             else
c Warning: tmp
              write(6,'(''iband,           orb_simple='',i3,2x, 5d15.8)')
     &        iband,orb_si(jorba),ddorb_si(jorba),(dorb_si(k,jorba),k=1,ndim)

              write(6,'(''iband,ireal_imag,orb_simple='',i3,i2,5d15.8)')
     &        iband,ireal_imag(jorba),orb_si(jorba+1),ddorb_si(jorba+1),(dorb_si(k,jorba+1),k=1,ndim)
c Save the orb_si and 2nd deriv in order to check for linear dependency
              orb_si(jorba)=orb_si(jorba+1)
              ddorb_si(jorba)=ddorb_si(jorba+1)
              do  k=1,ndim
                dorb_si(k,jorba)=dorb_si(k,jorba+1)
              enddo
            endif
           else
            ireal_imag(jorba-1)=0
            ireal_imag(jorba)=0
            write(6,'(''iband,ireal_imag,orb_simple='',i3,i2,5d15.8)')
     &      iband,ireal_imag(jorba-1),orb_si(jorba-1),ddorb_si(jorba-1),(dorb_si(k,jorba-1),k=1,ndim)
            write(6,'(''iband,ireal_imag,orb_simple='',i3,i2,5d15.8)')
     &      iband,ireal_imag(jorba),orb_si(jorba),ddorb_si(jorba),(dorb_si(k,jorba),k=1,ndim)
          endif

   80   continue
   90   nband(ikv)=nband_tmp
      call my_second(2,'orb_si')

      write(6,'(2i4,'' orbitals read in from file orbitals_pw_pwscf'')') jorb,jorba
      if(jorba.lt.norb) then
        write(6,'(''jorba,norb='',2i5)') jorba,norb
        stop 'jorba < norb in read_orb_pw_pwscf'
      endif

      write(6,'(/,''energy of highest  occupied orbital='',f10.6)') eigmax
      if(eigmin.ne.1.d99) then
        write(6,'(''energy of lowest unoccupied orbital='',f10.6)') eigmin
       else
        write(6,'(''no unoccupied orbitals'')')
      endif
      if(eigmin.lt.eigmax) write(6,'(''Warning: Energy of lowest unoccupied
     &orbital is lower than that of the highest occupied one.'',/,
     &''This may be because you are calculating an excited state, or it may be
     & that you are calculating a ground state but the single particle'',/,
     &''band structure has a band overlap, though in fact there is none.'')')

      close(30)
      return

  999 write(6,'(''Error: file orbitals_pw_pwscf is missing'')')
      stop 'file orbitals_pw_pwscf is missing'
      end
c-----------------------------------------------------------------------
      subroutine orbitals_pw2(ikv,iband,jorba)
c Written by Cyrus Umrigar
c Calculate pw orbitals.
c isortg could be (but is not) used to map g-vectors from iv to ig.
c At present it is assumed that k-vectors are in the correct order, but
c if not one could use isortk to map iorb.
c Note: This is the straightforward, expensive evaluation for checking purposes only.

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'ewald.h'

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

      common /tempor_test/ igvec_dft(3,NGVEC_BIGX),iwgvec(NGVEC2X),c_real(NGVEC2X),c_imag(NGVEC2X)
     &,map_gvecdft_gvec(NGVEC2X),isign_gvecdft_gvec(NGVEC2X)
     &,orb(MORB),dorb(3,MORB),ddorb(MORB)
     &,orb_si(MORB),dorb_si(3,MORB),ddorb_si(MORB),iflag(MORB),rnorm,r(3)
     &,rkvec_tmp(3),rkvec_tmp2(3),ngg(MKPTS),ngvec_dft

c1    dimension r(3),orb(MELEC,*),dorb(3,MELEC,*),ddorb(MELEC,*)
c     dimension dcos_rp(3),dsin_rm(3),dcos_ip(3),dsin_im(3)
c    &,cos_g(MELEC,NGVECX),sin_g(MELEC,NGVECX),dcos_g(3,MELEC,NGVECX),dsin_g(3,MELEC,NGVECX)
c    &,ddcos_g(MELEC,NGVECX),ddsin_g(MELEC,NGVECX)
c    &,cos_k(MELEC,MKPTS),sin_k(MELEC,MKPTS),dcos_k(3,MELEC,MKPTS),dsin_k(3,MELEC,MKPTS)
c    &,ddcos_k(MELEC,MKPTS),ddsin_k(MELEC,MKPTS)
      dimension gvec_dft(3,NGVEC_BIGX),gnorm_dft(NGVEC_BIGX)

c     write(6,'(''nelec,norb,nkvec in orbitals_pw2'',9i5)') nelec,norb,nkvec

      do 26 ig=1,ngvec_dft
        do 23 k=1,ndim
          gvec_dft(k,ig)=0
          do 23 i=1,ndim
   23       gvec_dft(k,ig)=gvec_dft(k,ig)+igvec_dft(i,ig)*glatt(k,i)
c       write(6,'(/,''igvec_dft in recip. lat. vec. units'',9i4)') (igvec_dft(i,ig),i=1,ndim)
c  26   write(6,'(''gvec_dft in cartesian coodinates'',9f9.4)') (gvec_dft(k,ig),k=1,ndim)
   26 continue

      if(ipr.ge.1) then
        if(ikv.eq.1.and.iband.eq.1) then
          open(4,file='wf_alfe')
          write(4,'(i5,'' g-vectors are:'')') ngg(ikv)
c Warning this may not work if there is more than one k-vec because ngg(ikv) may be larger than ngg(1).
          write(4,'(3f19.15)') ((gvec_dft(k,iwgvec(ig)),k=1,ndim),ig=1,ngg(ikv))
        endif
        if(ikv.eq.1) write(4,'(''k-vector is'',3f12.8)') (rkvec(k,ikv),k=1,ndim)
        write(4,'(''iband='',i5)') iband
        write(4,'(1p2g23.15)') (c_real(ig),c_imag(ig),ig=1,ngg(ikv))
      endif

c Warning: c_real and c_imag are only stored for particular k-vector and band
c Warning: for the moment just use this for testing 1 electron
c1    nelecc=1
c1    do 80 i=1,nelecc
c     iorb=0
c     do 80 ikv=1,nkvec
c       do 80 iband=1,nband(ikv)

        iorb=jorba-k_inv(ikv)
c Calculate psi_+ orbital if that is the one kept or if both are kept
c       if(ireal_imag(iorb).eq.0 .or. ireal_imag(iorb).eq.1) then

        iorb=iorb+1
        orb_si(iorb)=0
        ddorb_si(iorb)=0
        do 29 k=1,ndim
   29     dorb_si(k,iorb)=0
        do 40 ig=1,ngg(ikv)
          ig2=iwgvec(ig)
          dot=0
          gnorm_dft(ig2)=0
          do 30 k=1,ndim
            gnorm_dft(ig2)=gnorm_dft(ig2)+(rkvec(k,ikv)+gvec_dft(k,ig2))**2
   30       dot=dot+(rkvec(k,ikv)+gvec_dft(k,ig2))*r(k)
          orb_si(iorb)=orb_si(iorb)+c_real(ig)*cos(dot)-c_imag(ig)*sin(dot)
          if(ipr.ge.4 .and. ig.le.22) write(6,'(''rkvec+gvec'',2i4,7f9.4,f18.12)')
     &    ig,ig2,(rkvec(k,ikv)+gvec_dft(k,ig2),k=1,ndim),c_real(ig),dot,cos(dot),sin(dot),orb_si(iorb)
          do 35 k=1,ndim
   35       dorb_si(k,iorb)=dorb_si(k,iorb)+(rkvec(k,ikv)+gvec_dft(k,ig2))*(-c_real(ig)*sin(dot)-c_imag(ig)*cos(dot))
   40     ddorb_si(iorb)=ddorb_si(iorb)-gnorm_dft(ig2)*(c_real(ig)*cos(dot)-c_imag(ig)*sin(dot))

        orb_si(iorb)=orb_si(iorb)*rnorm
        ddorb_si(iorb)=ddorb_si(iorb)*rnorm
        do 55 k=1,ndim
   55      dorb_si(k,iorb)=dorb_si(k,iorb)*rnorm

c       write(6,'(''ikv,iband,nband(ikv),k_inv(ikv)'',9i5)') ikv,iband,nband(ikv),k_inv(ikv)
c       write(6,'(''real orb_si='',i5,9d12.4)') iorb,orb_si(iorb)
c       endif

c Calculate psi_- orbital if that is the one kept or if both are kept
c       if(ireal_imag(iorb).eq.0 .or. ireal_imag(iorb).eq.2) then

        iorb=iorb+1
c This routine is just used for testing and the dimensioning is such that
c if norb=MORB then iorb could be MORB+1, because even if one is kept, both are
c calculated in order to decide which to keep.  So do the foll. check:
        if(iorb.gt.MORB) then
          write(6,*)'iorb,MORB=',iorb,morb
          stop 'iorb>MORB in orbitals_pw2'
        endif

        orb_si(iorb)=0
        ddorb_si(iorb)=0
        do 59 k=1,ndim
   59     dorb_si(k,iorb)=0
        do 70 ig=1,ngg(ikv)
          ig2=iwgvec(ig)
          dot=0
          gnorm_dft(ig2)=0
          do 60 k=1,ndim
            gnorm_dft(ig2)=gnorm_dft(ig2)+(rkvec(k,ikv)+gvec_dft(k,ig2))**2
   60       dot=dot+(rkvec(k,ikv)+gvec_dft(k,ig2))*r(k)
          orb_si(iorb)=orb_si(iorb)+c_real(ig)*sin(dot)+c_imag(ig)*cos(dot)
c         if(ipr.ge.4 .and. ig.le.22) write(6,'(''rkvec+gvec'',2i4,7f9.4,f18.12)')
c    &    ig,ig2,(rkvec(k,ikv)+gvec_dft(k,ig2),k=1,ndim),c_real(ig),dot,cos(dot),sin(dot),orb_si(iorb)
          if(ipr.ge.4 .and. ig.le.22) write(6,'(''rkvec,gvec'',8f9.4)') (rkvec(k,ikv),gvec_dft(k,ig2),k=1,ndim)
          if(ipr.ge.4 .and. ig.le.22) write(6,'(''rkvec+gvec'',2i4,8f9.4,f18.12)')
     &    ig,ig2,(rkvec(k,ikv)+gvec_dft(k,ig2),k=1,ndim),c_real(ig),r(1),dot,cos(dot),sin(dot),orb_si(iorb)
          do 65 k=1,ndim
   65       dorb_si(k,iorb)=dorb_si(k,iorb)+(rkvec(k,ikv)+gvec_dft(k,ig2))*(c_real(ig)*cos(dot)-c_imag(ig)*sin(dot))
   70     ddorb_si(iorb)=ddorb_si(iorb)-gnorm_dft(ig2)*(c_real(ig)*sin(dot)+c_imag(ig)*cos(dot))

        orb_si(iorb)=orb_si(iorb)*rnorm
        ddorb_si(iorb)=ddorb_si(iorb)*rnorm
        do 75 k=1,ndim
   75      dorb_si(k,iorb)=dorb_si(k,iorb)*rnorm

c       write(6,'(''ikv,iband,nband(ikv),k_inv(ikv)'',9i5)') ikv,iband,nband(ikv),k_inv(ikv)
c       write(6,'(''imag orb_si='',i5,9d12.4)') iorb,orb_si(iorb)

c       endif
c       endif

c  80   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine print_orbitals_pw

      use bwfdet_mod
      implicit real*8(a-h,o-z)
      character*20 fmt
      character*16 mode,iorb_format

      include 'vmc.h'
      include 'force.h'
      include 'ewald.h'
      include 'numorb.h'
      parameter(eps=1.d-3)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode
      common /contr_names/ iorb_format
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /dorb/ iworbd(MELEC,MDET)

      common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,ngrid_orbx,ngrid_orby,ngrid_orbz
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange
c    &,orb_splines(8,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:
c    &MGRID_ORB_PER-1,MORB_OCC)
c    &,grid_orbx(0:MGRID_ORB_PER-1),grid_orby(
c    &0:MGRID_ORB_PER-1),grid_orbz(0:MGRID_ORB_PER-1)
c    &,orb_splines_explicit(4,4,4,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:
c    &MGRID_ORB_PER-1,MORB_OCC)

      common /periodic2/ rkvec_shift_latt(3)

      dimension r(3),r_basis(3),r_test(3,MGRID_ORB_PER**3)
      dimension orb(MORB_OCC),dorb(3,MORB_OCC),ddorb(MORB_OCC)
c     dimension orb_splines_tmp(10),ict(10),ddorb_splines_tmp(3)
      dimension orb_tmp(MORB_OCC),dorb_tmp(3,MORB_OCC),
     &          ddorb_tmp(MORB_OCC)
      dimension orb_blip_tmp(MORB_OCC,MDET),
     &          dorb_blip_tmp(3,MORB_OCC,MDET),
     &          ddorb_blip_tmp(MORB_OCC,MDET)


      integer i,k,ix,iy,iz,iorb,ier,isgn,npts,npts_max,xfac,yfac,zfac


      ngrid_orbx=MGRID_ORB_PER
      ngrid_orby=MGRID_ORB_PER
      ngrid_orbz=MGRID_ORB_PER

      open(14,file='print_orbitals.data')
      read(14,*)npts,xfac,yfac,zfac
      close(14)

      npts_max=ngrid_orbx*ngrid_orby*ngrid_orbz

      if(npts.gt.npts_max) stop 'npts must be less than MGRID_ORB_PER^3'

      if(xfac.ne.0) then
         do i=0,npts-1
            r_test(1,i)=i/dfloat(npts-1)
            r_test(2,i)=r_test(1,i)*dfloat(yfac)/dfloat(xfac)
            r_test(3,i)=r_test(1,i)*dfloat(zfac)/dfloat(xfac)
         enddo
      elseif(yfac.ne.0) then
         do iy=0,npts-1
            r_test(1,i)=0.d0
            r_test(2,i)=iy/dfloat(npts-1)
            r_test(3,i)=r_test(2,i)*dfloat(zfac)/dfloat(yfac)
         enddo
      elseif(zfac.ne.0) then
         do iz=0,npts-1
            r_test(1,i)=0.d0
            r_test(2,i)=0.d0
            r_test(3,i)=i/dfloat(npts-1)
         enddo
      else
         stop '000 is not a valid direction'
      endif

c File open
      if(inum_orb.eq.0) then
          open(15,file='pw_test.data')
          write(15,*)'#',npts,'points'
          write(15,*)'#',xfac,yfac,zfac,' direction'
      endif
      if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
          open(15,file='lagrange_test.data')
          write(15,*)'#',npts,'points'
          write(15,*)'#',xfac,yfac,zfac,' direction'
      endif
!     if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
!         open(15,file='spline_test.data')
!         write(15,*)'#',npts,'points'
!         write(15,*)'#',xfac,yfac,zfac,' direction'
!        do i=1,10
!           ict(i)=1
!        enddo
!     endif
      if(inum_orb.eq.6 .or. inum_orb.eq.-6) then
          open(15,file='blip_test.data')
          write(15,*)'#',npts,'points'
          write(15,*)'#',xfac,yfac,zfac,' direction'
      endif



      do ix=0,npts-1
         r_basis(1)=r_test(1,ix)
         r_basis(2)=r_test(2,ix)
         r_basis(3)=r_test(3,ix)

c        do ix=0,2*(ngrid_orbx)*(ngrid_orbx-1)-2
c           r_basis(1)=ix/dfloat(2*(ngrid_orbx)*(ngrid_orbx-1)-2)

c           do iy=0,2*(ngrid_orby)*(ngrid_orby-1)-2
c             r_basis(2)=iy/dfloat(2*(ngrid_orby)*(ngrid_orby-1)-2)

c              do iz=0,2*(ngrid_orbz)*(ngrid_orbz-1)-2
c                r_basis(3)=iz/dfloat(2*(ngrid_orbz)*(ngrid_orbz-1)-2)

                 isgn=1
                 do k=1,ndim
                   if(rkvec_shift_latt(k).ne.0.d0) then
                     if(r_basis(k).ge.0.d0) then
                       isgn=isgn*(-1)**int(r_basis(k))
                      else
                       isgn=isgn*(-1)**(int(r_basis(k))+1)
                     endif
                   endif
                 enddo

cwparker Go along the <100> direction
c                if (iy.eq.0 .and. iz.eq.0) then
cwparker Go along the <010> direction
c                if (ix.eq.0 .and. iz.eq.0) then
cwparker Go along the <001> direction
c                if (ix.eq.0 .and. iy.eq.0) then
cwparker Go along the <110> direction
c                if (ix.eq.iy .and. iz.eq.0) then
cwparker Go along the <111> direction
c                if (ix.eq.iy .and. iy.eq.iz) then

                    do k=1,ndim
                       r(k)=0
                       do i=1,ndim
                          r(k)=r(k)+rlatt_sim(k,i)*r_basis(i)
                       enddo
                    enddo

cwparker Plane wave printout
                 if(inum_orb.eq.0) then

                   call orbitals_pw_grade(r,orb,dorb,ddorb)
                   do iorb=1,norb
                      if(xfac.ne.0) then
cwparker write statement for xyz directions
                       write(15,'(''r(1),orb('',i2,''),dorb('',i2,''),''
     &                   ,''ddorb('',i2,'')='',6f14.8)')iorb,iorb,iorb,
     &                   r_basis(1),orb(iorb),dorb(1,iorb),dorb(2,iorb),
     &                   dorb(3,iorb),ddorb(iorb)
                      elseif(yfac.ne.0) then
cwparker write statement for 0yz directions
                       write(15,'(''r(2),orb('',i2,''),dorb('',i2,''),''
     &                   ,''ddorb('',i2,'')='',6f14.8)')iorb,iorb,iorb,
     &                   r_basis(2),orb(iorb),dorb(1,iorb),dorb(2,iorb),
     &                   dorb(3,iorb),ddorb(iorb)
                      elseif(zfac.ne.0) then
cwparker write statement for 00z directions
                       write(15,'(''r(3),orb('',i2,''),dorb('',i2,''),''
     ^                   ,''ddorb('',i2,'')='',6f14.8)')iorb,iorb,iorb,
     &                   r_basis(3),orb(iorb),dorb(1,iorb),dorb(2,iorb),
     &                   dorb(3,iorb),ddorb(iorb)
                      endif
                   enddo
                 endif
cwparker End of plane wave printout

cwparker Lagrange polynomial printout
                 if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
                   xi=r_basis(1)*ngrid_orbx
                   yi=r_basis(2)*ngrid_orby
                   zi=r_basis(3)*ngrid_orbz

                   call interpol_orb(ngrid_orbx,ngrid_orby,ngrid_orbz,
     &               xi,yi,zi,orb_tmp,dorb_tmp,ddorb_tmp)

                   do iorb=1,norb
                      if(xfac.ne.0) then
cwparker write statement for xyz directions
                        write(15,'(''r(1),orb('',i2,''),dorb('',i2,''),
     &ddorb('',i2,'')='',6f14.8)')iorb,iorb,iorb,r_basis(1),
     &isgn*orb_tmp(iorb),isgn*dorb_tmp(1,iorb),isgn*dorb_tmp(2,iorb),
     &isgn*dorb_tmp(3,iorb),isgn*ddorb_tmp(iorb)
cwparker write statement for 0yz directions
                      elseif(yfac.ne.0) then
                        write(15,'(''r(2),orb('',i2,''),dorb('',i2,''),
     &ddorb('',i2,'')='',6f14.8)')iorb,iorb,iorb,r_basis(2),
     &isgn*orb_tmp(iorb),isgn*dorb_tmp(1,iorb),isgn*dorb_tmp(2,iorb),
     &isgn*dorb_tmp(3,iorb),isgn*ddorb_tmp(iorb)
cwparker write statement for 00z directions
                      elseif(zfac.ne.0) then
                        write(15,'(''r(3),orb('',i2,''),dorb('',i2,''),
     &ddorb('',i2,'')='',6f14.8)')iorb,iorb,iorb,r_basis(3),
     &isgn*orb_tmp(iorb),isgn*dorb_tmp(1,iorb),isgn*dorb_tmp(2,iorb),
     &isgn*dorb_tmp(3,iorb),isgn*ddorb_tmp(iorb)
                      endif
                   enddo
                 endif
cwparker End of Lagrange polynomial printout

cwparker Spline printout
!                 if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
!
!                   do iorb=1,norb
!
!cwparker Compact splines
!                      call r8evtricub(r_basis(1),r_basis(2),r_basis(3),
!     &                                grid_orbx,ngrid_orbx,grid_orby,
!     &                                ngrid_orby,grid_orbz,ngrid_orbz
!     &                                ,1,1,1,orb_splines(1,0,0,0,iorb),
!     &                                ngrid_orby,ngrid_orbz,
!     &                                ict,orb_splines_tmp,ier)
!
!             if(ier.ne.0) stop 'error in r8evtricub'
!cwparker Explicit splines
!c                     call r8tcspeval(r_basis(1),r_basis(2),r_basis(3),
!c    &                                ict,orb_splines_tmp,grid_orbx,
!c    &                                ngrid_orbx,grid_orby,ngrid_orby,
!c    &                                grid_orbz,ngrid_orbz,1,1,1,
!c    &                                orb_splines_explicit(1,1,1,0,0,0,
!c    &                                                            iorb),
!c    &                                ngrid_orbx,ngrid_orby,ier)
!
!             if(ier.ne.0) stop 'error in r8tcspeval'
!
!             orb_tmp(iorb)    = orb_splines_tmp(1)
!
!             dorb_tmp(1,iorb) = rlatt_sim_inv(1,1)*orb_splines_tmp(2)+
!     &                          rlatt_sim_inv(1,2)*orb_splines_tmp(3)+
!     &                          rlatt_sim_inv(1,3)*orb_splines_tmp(4)
!
!             dorb_tmp(2,iorb) = rlatt_sim_inv(2,1)*orb_splines_tmp(2)+
!     &                          rlatt_sim_inv(2,2)*orb_splines_tmp(3)+
!     &                          rlatt_sim_inv(2,3)*orb_splines_tmp(4)
!
!             dorb_tmp(3,iorb) = rlatt_sim_inv(3,1)*orb_splines_tmp(2)+
!     &                          rlatt_sim_inv(3,2)*orb_splines_tmp(3)+
!     &                          rlatt_sim_inv(3,3)*orb_splines_tmp(4)
!
!              ddorb_splines_tmp(1)=rlatt_sim_inv(1,1)**2*
!     &                      orb_splines_tmp(5)+
!     &                      2*rlatt_sim_inv(1,1)*rlatt_sim_inv(1,2)*
!     &                      orb_splines_tmp(8)+rlatt_sim_inv(1,2)**2*
!     &                      orb_splines_tmp(6)+2*rlatt_sim_inv(1,2)*
!     &                      rlatt_sim_inv(1,3)*orb_splines_tmp(10)+
!     &                      rlatt_sim_inv(1,3)**2*orb_splines_tmp(7)+
!     &                      2*rlatt_sim_inv(1,3)*rlatt_sim_inv(1,1)*
!     &                      orb_splines_tmp(9)
!
!             ddorb_splines_tmp(2)=rlatt_sim_inv(2,1)**2*
!     &                      orb_splines_tmp(5)+
!     &                      2*rlatt_sim_inv(2,1)*rlatt_sim_inv(2,2)*
!     &                      orb_splines_tmp(8)+rlatt_sim_inv(2,2)**2*
!     &                      orb_splines_tmp(6)+2*rlatt_sim_inv(2,2)*
!     &                      rlatt_sim_inv(2,3)*orb_splines_tmp(10)+
!     &                      rlatt_sim_inv(2,3)**2*orb_splines_tmp(7)+
!     &                      2*rlatt_sim_inv(2,3)*rlatt_sim_inv(2,1)*
!     &                      orb_splines_tmp(9)
!
!             ddorb_splines_tmp(3)=rlatt_sim_inv(3,1)**2*
!     &                      orb_splines_tmp(5)+
!     &                      2*rlatt_sim_inv(3,1)*rlatt_sim_inv(3,2)*
!     &                      orb_splines_tmp(8)+rlatt_sim_inv(3,2)**2*
!     &                      orb_splines_tmp(6)+2*rlatt_sim_inv(3,2)*
!     &                      rlatt_sim_inv(3,3)*orb_splines_tmp(10)+
!     &                      rlatt_sim_inv(3,3)**2*orb_splines_tmp(7)+
!     &                      2*rlatt_sim_inv(3,3)*rlatt_sim_inv(3,1)*
!     &                      orb_splines_tmp(9)
!
!             ddorb_tmp(iorb) =  ddorb_splines_tmp(1)+
!     &                          ddorb_splines_tmp(2)+
!     &                          ddorb_splines_tmp(3)
!
!             orb(iorb)=orb_tmp(iorb)*isgn
!             do k=1,ndim
!                dorb(k,iorb)=dorb_tmp(k,iorb)*isgn
!             enddo
!             ddorb(iorb)=ddorb_tmp(iorb)*isgn
!
!                      if(xfac.ne.0) then
!cwparker write statement for xyz directions
!                        write(15,'(''r(1),orb('',i2,''),dorb('',i2,''),
!     &ddorb('',i2,'')= '',6f14.8)')iorb,iorb,iorb,r_basis(1),
!     &orb_tmp(iorb),dorb_tmp(1,iorb),dorb_tmp(2,iorb),
!     &dorb_tmp(3,iorb),ddorb_tmp(iorb)
!                      elseif(yfac.ne.0) then
!cwparker write statement for 0yz directions
!                        write(15,'(''r(2),orb('',i2,''),dorb('',i2,''),
!     &ddorb('',i2,'')= '',6f14.8)')iorb,iorb,iorb,r_basis(2),
!     &orb_tmp(iorb),dorb_tmp(1,iorb),dorb_tmp(2,iorb),
!     &dorb_tmp(3,iorb),ddorb_tmp(iorb)
!                      elseif(zfac.ne.0) then
!cwparker write statement for 00z directions
!                        write(15,'(''r(3),orb('',i2,''),dorb('',i2,''),
!     &ddorb('',i2,'')= '',6f14.8)')iorb,iorb,iorb,r_basis(3),
!     &orb_tmp(iorb),dorb_tmp(1,iorb),dorb_tmp(2,iorb),
!     &dorb_tmp(3,iorb),ddorb_tmp(iorb)
!                      endif
!                   enddo
!
!                 endif
!cwparker End of Spline printout

cwparker Blip printout
                 if(inum_orb.eq.6 .or. inum_orb.eq.-6) then
                    call bwfdet_main(r,1,1,1,orb_blip_tmp,
     &                               dorb_blip_tmp,ddorb_blip_tmp)

                   do iorb=1,norb
                      if(xfac.ne.0) then
cwparker write statement for 1yz directions
                        write(15,'(''r(1),orb('',i2,''),dorb('',i2,''),
     &ddorb('',i2,''),'',6f14.8)') iorb,iorb,iorb,r_basis(1),
     &orb_blip_tmp(iorb,1),dorb_blip_tmp(1,iorb,1),
     &dorb_blip_tmp(2,iorb,1),dorb_blip_tmp(3,iorb,1),
     &ddorb_blip_tmp(iorb,1)
                      elseif(yfac.ne.0) then
cwparker write statement for x1z directions
                        write(15,'(''r(2),orb('',i2,''),dorb('',i2,''),
     &ddorb('',i2,''),'',6f14.8)') iorb,iorb,iorb,r_basis(2),
     &orb_blip_tmp(iorb,1),dorb_blip_tmp(1,iorb,1),
     &dorb_blip_tmp(2,iorb,1),dorb_blip_tmp(3,iorb,1),
     &ddorb_blip_tmp(iorb,1)
                      elseif(zfac.ne.0) then
cwparker write statement for xy1 directions
                        write(15,'(''r(3),orb('',i2,''),dorb('',i2,''),
     &ddorb('',i2,''),'',6f14.8)') iorb,iorb,iorb,r_basis(3),
     &orb_blip_tmp(iorb,1),dorb_blip_tmp(1,iorb,1),
     &dorb_blip_tmp(2,iorb,1),dorb_blip_tmp(3,iorb,1),
     &ddorb_blip_tmp(iorb,1)
                      endif
                   enddo
                 endif
cwparker End of blip printout

c                endif

cwparker End of r_test do loop
         enddo

cwparker End of directional choice (<100>,<010>,<001> or  <110>)
c              enddo
c           enddo
c        enddo
cwparker End of grid do loops

cwparker Close test data file
       close(15)

cwparker Stop with appropriate message
       if(inum_orb.eq.0) then
         stop 'Done with plane wave printout; see pw_test.data'
       endif
       if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
         stop 'Done with Lagrange printout; see lagrange_test.data'
       endif
c      if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
c        stop 'Done with pp-spline prinout; see spline_test.data'
c      endif
       if(inum_orb.eq.6 .or. inum_orb.eq.-6) then
         stop 'Done with B-spline printout; see blip_test.data'
       endif

       stop 'No orbital printout made; inum_orb invalid'

       return
       end

c-----------------------------------------------------------------------
c     subroutine energy_test

c     implicit real*8(a-h,o-z)
c     character*20 fmt
c     character*16 mode,iorb_format

c     include 'vmc.h'
c     include 'force.h'
c     include 'ewald.h'
c     include 'numorb.h'
c     parameter(eps=1.d-3)

c     common /dim/ ndim
c     common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
c     common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
c    &,ifock,i3body,irewgt,iaver,istrch
c    &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
c     common /contr3/ mode
c     common /contr_names/ iorb_format
c     common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
c     common /dorb/ iworbd(MELEC,MDET)

c     common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
c    &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
c    &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
c    &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
c    &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
c    &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
c    &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
c    &,znuc_sum,znuc2_sum,vcell,vcell_sim
c    &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
c    &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
c    &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange
c    &,orb_splines(8,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:
c    &MGRID_ORB_PER-1,MORB_OCC),grid_orbx(0:MGRID_ORB_PER-1),grid_orby(
c    &0:MGRID_ORB_PER-1),grid_orbz(0:MGRID_ORB_PER-1)

c     common /periodic2/ rkvec_shift_latt(3)

c     dimension r(3),r_basis(3)
c     integer i,k,ix,iy,iz,iorb


c     ngrid_orbx=MGRID_ORB_PER
c     ngrid_orby=MGRID_ORB_PER
c     ngrid_orbz=MGRID_ORB_PER

c     if(inum_orb.eq.0) then
c         open(15,file='pw_energy_test.data')
c     endif
c     if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
c         open(15,file='lagrange_energy_test.data')
c     endif
c     if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
c         open(15,file='spline_energy_test.data')
c     endif
c     if(inum_orb.eq.6 .or. inum_orb.eq.-6) then
c         open(15,file='blip_energy_test.data')
c     endif

c        do ix=0,2*(ngrid_orbx)*(ngrid_orbx-1)-2
c           r_basis(1)=ix/dfloat(2*(ngrid_orbx)*(ngrid_orbx-1)-2)

c           do iy=0,2*(ngrid_orby)*(ngrid_orby-1)-2
c             r_basis(2)=iy/dfloat(2*(ngrid_orby)*(ngrid_orbx-1)-2)

c              do iz=0,2*(ngrid_orbz)*(ngrid_orbz-1)-2
c                r_basis(3)=iz/dfloat(2*(ngrid_orbz)*(ngrid_orbz-1)-2)

cwparker Go along the <100> direction
c                if (iy.eq.0 .and. iz.eq.0) then
cwparker Go along the <001> direction
c                if (ix.eq.0 .and. iy.eq.0) then
cwparker Go along the <110> direction
c                if (ix.eq.iy .and. iz.eq.0) then
cwparker Go along the <111> direction
c                if (ix.eq.iy .and. iy.eq.iz) then

c                   do k=1,ndim
c                      r(k)=0
c                      do i=1,ndim
c                         r(k)=r(k)+rlatt_sim(k,i)*r_basis(i)
c                      enddo
c                   enddo

c                   call hpsi(r,determinant,jastrow,velocity,div_v,
c    &                        laplacian,potential,eeinteraction,energy,
c    &                        denergy,ifr)
c                   write(15,'(''r(1),E,V='',3f9.4)')
c    &                   r_basis(1),energy,potential

c                endif
cwparker End of directional choice (<100>,<110> or <100>)

c              enddo
c           enddo
c        enddo
cwparker End of grid do loops

cwparker Close test data file
c      close(15)

cwparker Stop with appropriate message
c      if(inum_orb.eq.0) then
c        stop 'done with plane wave test; see pw_energy_test.data'
c      endif
c      if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
c        stop 'done with Lagrange test; see lagrange_energy_test.data'
c      endif
c      if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
c        stop 'done with spline test; see spline_energy_test.data'
c      endif
c      if(inum_orb.eq.6 .or. inum_orb.eq.-6) then
c        stop 'done with blip test; see blip_energy_test.data'
c      endif

c      stop 'no test done; inum_orb not valid'

c      return
c      end


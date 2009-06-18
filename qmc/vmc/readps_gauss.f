      subroutine readps_gauss
c Written by Claudia Filippi (or Friedemann Schautz?), modified by Cyrus Umrigar
c read 'Quantum-chemist' gauss pseudopotentials
c file format: one text file with basename gauss_ecp.dat
c              for each atom type
c first line : arbitrary label (written to log-file)
c second line: number of projectors + 1 (i.e. total number of components)
c remaining lines: components in the order (local,L=0,L=1 ...)
c     repeated for each component
c        number terms
c        repeated for each term in this component
c          coefficient, power, exponent

c NOTE: as usual power n means r**(n-2)

      use all_tools_mod
      use atom_mod
      use pseudo_mod
      use qua_mod
      use gauss_ecp_mod
      implicit real*8(a-h,o-z)

      character*80 label
      character*20 filename,atomtyp
      integer MGAUSS

      call alloc ('necp_term', necp_term, MPS_L, nctype)

      MGAUSS=0

      do 100 ict=1,nctype

        if(ict.lt.10) then
          write(atomtyp,'(i1)') ict
         elseif(ict.lt.100) then
          write(atomtyp,'(i2)') ict
         else
          stop 'readps_tm, nctype>100'
        endif

        filename='pseudopot'//atomtyp(1:index(atomtyp,' ')-1)
        open(1,file=filename,status='old',form='formatted',err=999)
        write(6,'(/3a)') ' Reading gaussian pseudopotential file >', trim(filename),'<'

        read(1,'(a100)',err=1000,end=1001) label
        write(6,'('' ECP for center type'',i4,'' label= '',a80)') ict,label

c Number of l components
        read(1,*,err=1000,end=1001) npotd(ict)

c If the local pseudopot component is not set in input, set it here
        if(lpotp1(ict).lt.0) then
          lpotp1(ict)=npotd(ict)
          write(6,'('' Center type'',i4,'' local pseudopot component is'',i3)') ict,lpotp1(ict)
        endif

        write(6,'('' Center type'',i2,'' has'',i2,'' pseudopotential L components, and component''
     &  ,i2,'' is chosen to be local'')') ict,npotd(ict),lpotp1(ict)

        if(lpotp1(ict).gt.npotd(ict)) then
          write(6,'(''lpotp1(ict),npotd(ict)='',2i3)') lpotp1(ict),npotd(ict)
          stop 'Cannot choose local psp. to be > number of l components, lpotp1(ict) > npotd(ict)'
        endif
        if(npotd(ict).gt.MPS_L) stop 'npotd(ict).gt.MPS_L'

c read terms of local part and all non-local parts
c local part first in file, but stored at index lpotp1
c since psps. are stored in order of ascending l.

        do 50 l=1,npotd(ict)
c         if(l.eq.lpotp1(ict))then
          if(l.eq.1)then
            idx=lpotp1(ict)
           else
            idx=l-1
          endif
          read(1,*,err=999,end=1000) necp_term(idx,ict)
          MGAUSS = max(MGAUSS, necp_term(idx,ict)) !JT

!JT          if(necp_term(idx,ict).gt.MGAUSS) stop 'necp_term(idx,ict) > MGAUSS'
          write(6,'('' L component, #terms='',2i4)') l,necp_term(idx,ict)
          call alloc ('ecp_coef', ecp_coef, MGAUSS, MPS_L, nctype)
          call alloc ('necp_power', necp_power, MGAUSS, MPS_L, nctype)
          call alloc ('ecp_exponent', ecp_exponent, MGAUSS, MPS_L, nctype)
          do 50 i=1,necp_term(idx,ict)
            read(1,*,err=999,end=1000) ecp_coef(i,idx,ict),
     &        necp_power(i,idx,ict),ecp_exponent(i,idx,ict)
            write(6,'('' Psp coef, power, expo'',f16.8,i2,f16.8)') ecp_coef(i,idx,ict),necp_power(i,idx,ict)
     &        ,ecp_exponent(i,idx,ict)

   50   continue

c Find the point beyond which the various v components differ from each other by no more than .5*d-6
c The angular integration needs to be done only within rmax_nloc, so the smaller this number is the
c less the computational expense per MC step.
        r=5.d0
        do 70 ir=1,500
          do 70 l=1,npotd(ict)
            if(l.ne.lpotp1(ict)) then
              call gauss_pot(r,l,ict,vpot)
              if(dabs(vpot).gt..5d-6) goto 80
            endif
   70   r=r-0.01d0
   80   rmax_nloc=r
  100 write(6,'('' center '',i3,'' pseudopot rmax_nloc= '',f6.2)') ict,rmax_nloc

c     call gesqua(nquad,xq0,yq0,zq0,wq)

      return

  999 write(6,'(''Error: Pseudopot. file '',a20,'' is missing'')') filename
      stop 'Pseudopot. file is missing'

 1000 write(6,'(''Error when reading gaussian psp. of center type'',i4)') ict
      stop 'readps_gauss: error while reading gaussian pseudopotential'

 1001 write(6,'(''End of file when reading gaussian psp. of center type'',i4)') ict
      stop 'readps_gauss: end of file while reading gaussian pseudopotential'

      end
c-----------------------------------------------------------------------

c compute gauss-pseudopotential for electron iel
      subroutine getvps_gauss(r_en,iel)

      use atom_mod
      use const_mod
      use pseudo_mod
      implicit real*8(a-h,o-z)

      dimension r_en(nelec,ncent)

      do 10 ic=1,ncent
        ict=iwctype(ic)
        r=max(1.d-10,r_en(iel,ic))
        do 10 l=1,npotd(ict)
          call gauss_pot(r,l,ict,vpot)
          if(l.eq.lpotp1(ict)) vpot=vpot-znuc(ict)/r
   10     vps(iel,ic,l)=vpot

      return
      end
c-----------------------------------------------------------------------
      subroutine gauss_pot(r,l,ict,vpot)
      use gauss_ecp_mod
      implicit real*8(a-h,o-z)

      v=0
      rsq=r**2

      do i=1,necp_term(l,ict)
        if(necp_power(i,l,ict).ne.2)then
          p=r**(necp_power(i,l,ict)-2)
         else
          p=1
        endif
        e=ecp_coef(i,l,ict)*exp(-ecp_exponent(i,l,ict)*rsq)
        v=v+p*e
      enddo
      vpot=v

      return
      end

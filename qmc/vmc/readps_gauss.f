      subroutine readps_gauss
c Written by Claudia Filippi (or Friedemann Schautz?)
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
c
c NOTE: as usual power n means r**(n-2)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'

      character*80 label
      character*20 filename,atomtyp

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc

      common /gauss_ecp/ necp_term(MPS_L,MCTYPE),necp_power(MGAUSS
     &     ,MPS_L,MCTYPE),ecp_coef(MGAUSS,MPS_L,MCTYPE)
     &     ,ecp_exponent(MGAUSS,MPS_L,MCTYPE)

      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad

      do 200 ic=1,nctype

        if(ic.lt.10) then
          write(atomtyp,'(i1)') ic
         elseif(ic.lt.100) then
          write(atomtyp,'(i2)') ic
         else
          stop 'readps_tm, nctype>100'
        endif

        filename='pseudopot'//atomtyp(1:index(atomtyp,' ')-1)
        open(1,file=filename,status='old',form='formatted',err=999)
        write(6,'(''Reading gaussian pseudopotential file '',a)') filename

c label 
        read(1,'(a100)',err=1000,end=1001) label
        write(6,'(''ECP for center type'',i4,'' label= '',a80)') ic,label

c Number of l components
        read(1,*,err=1000,end=1001) npotd(ic)

c If the local pseudopot component is not set in input, set it here
        if(lpotp1(ic).lt.0) then
          lpotp1(ic)=npotd(ic)
          write(6,'(''Center type'',i4,'' local pseudopot component is'',i3)') ic,lpotp1(ic)
        endif

        write(6,'(''Center type'',i2,'' has'',i2,'' pseudopotential L components, and component''
     &  ,i2,'' is chosen to be local'')') ic,npotd(ic),lpotp1(ic)
 
        if(lpotp1(ic).gt.npotd(ic)) then
          write(6,'(''lpotp1(ic),npotd(ic)='',2i3)') lpotp1(ic),npotd(ic)
          stop 'Cannot choose local psp. to be > number of l components, lpotp1(ic) > npotd(ic)'
        endif
        if(npotd(ic).gt.MPS_L) stop 'npotd(ic).gt.MPS_L'

c read terms of local part and all non-local parts
c local part first in file, but stored at index lpotp1
c since psps. are stored in order of ascending l.

        do 200 l=1,npotd(ic)
c         if(l.eq.lpotp1(ic))then
          if(l.eq.1)then
            idx=lpotp1(ic)
           else
            idx=l-1
          endif
          read(1,*,err=999,end=1000) necp_term(idx,ic)

          if(necp_term(idx,ic).gt.MGAUSS) stop 'necp_term(idx,ic) > MGAUSS'
          write(6,'(''L component, #terms='',2i4)') l,necp_term(idx,ic)
          do 200 i=1,necp_term(idx,ic)
            read(1,*,err=999,end=1000) ecp_coef(i,idx,ic),
     &        necp_power(i,idx,ic),ecp_exponent(i,idx,ic)
            write(6,'(''Psp coef, power, expo'',f16.8,i2,f16.8)') ecp_coef(i,idx,ic),necp_power(i,idx,ic)
     &        ,ecp_exponent(i,idx,ic)
  200   continue

c     call gesqua(nquad,xq0,yq0,zq0,wq)
c     call gesqua(nquad,xq0,yq0,zq0,wq)

      return

  999 write(6,'(''Error: Pseudopot. file '',a20,'' is missing'')') filename
      stop 'Pseudopot. file is missing'

 1000 write(6,'(''Error when reading gaussian psp. of center type'',i4)') ic
      stop 'readps_gauss: error while reading gaussian pseudopotential'

 1001 write(6,'(''End of file when reading gaussian psp. of center type'',i4)') ic
      stop 'readps_gauss: end of file while reading gaussian pseudopotential'

      end
c-----------------------------------------------------------------------

c compute gauss-pseudopotential for electron iel
      subroutine getvps_gauss(r_en,iel)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc

      dimension r_en(MELEC,MCENT)

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
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      common /gauss_ecp/ necp_term(MPS_L,MCTYPE),necp_power(MGAUSS
     &     ,MPS_L,MCTYPE),ecp_coef(MGAUSS,MPS_L,MCTYPE)
     &     ,ecp_exponent(MGAUSS,MPS_L,MCTYPE)

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

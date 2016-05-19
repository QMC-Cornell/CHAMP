      subroutine readps
! Written by Claudia Filippi
! Read pseudopotential in Fahy format
      use all_tools_mod
      use atom_mod
      use pseudo_mod
      use qua_mod
      use pseudo_fahy_mod
      implicit real*8(a-h,o-z)

      character*20 filename,atomtyp

! nquad = number of quadrature points
! nlang = number of non-local potentials
! rcmax = cutoff radius for non-local potential
! npotl = number of mesh point for local potential
! dradl = spacing of uniform mesh for local potential

      call alloc ('potl', potl, MPS_GRID, nctype)
      call alloc ('ptnlc', ptnlc, MPS_GRID, nctype, MPS_L)
      call alloc ('ptnlc', ptnlc, MPS_GRID, nctype, MPS_L)
      call alloc ('dradl', dradl, nctype)
      call alloc ('drad', drad, nctype)
      call alloc ('rcmax', rcmax, nctype)
      call alloc ('npotl', npotl, nctype)
      call alloc ('nlrad', nlrad, nctype)

      do 20 ict=1,nctype

      if(ict.lt.10) then
        write(atomtyp,'(i1)') ict
       elseif(ict.lt.100) then
        write(atomtyp,'(i2)') ict
       else
        stop 'readps_tm, problem atomtyp'
      endif

      filename='ps.data.'//atomtyp(1:index(atomtyp,' ')-1)
      open(3,file=filename,status='old',form='formatted')

      read(3,*) nquad
      write(6,'(''quadrature points'',i4)') nquad

      read(3,*) nlang,rcmax(ict)
! If the local pseudopot component is not set in input, set it here
      if(lpotp1(ict).lt.0) then
        lpotp1(ict)=nlang+1
        write(6,'(''local pseudopot component is'',i3)') lpotp1(ict)
      endif

! local potential
      read(3,*)
      read(3,*) npotl(ict),nzion,dradl(ict)
      if(npotl(ict).gt.MPS_GRID) stop 'npotl gt MPS_GRID'
      if(nzion.ne.int(znuc(iwctype(ict)))) stop 'nzion ne znuc'

      read(3,*) (potl(i,ict),i=1,npotl(ict))

      do 5 i=1,npotl(ict)
  5     write(33,*) (i-1)*dradl(ict),potl(i,ict)

! non-local potential
      read(3,*)
      read(3,*) nlrad(ict),drad(ict)
      if(nlrad(ict).gt.MPS_GRID) stop 'nlrad gt MPS_GRID'

      if(drad(ict)*(nlrad(ict)-1).le.rcmax(ict)) then
        write(6,'(''non-local table max radius = '', &
     &  f10.5,'' too small for cut-off = '',f10.5)') &
     &  drad(ict)*(nlrad(ict)-1),rcmax(ict)
        stop
      endif

      do 10 l=1,nlang
        read(3,*)
        read(3,*) (ptnlc(i,ict,l),i=1,nlrad(ict))
      do 10 i=1,nlrad(ict)
 10     write(34,*) (i-1)*drad(ict),ptnlc(i,ict,l)

      close(3)
 20   continue


!     call gesqua (nquad,xq0,yq0,zq0,wq)
!     call gesqua (nquad,xq0,yq0,zq0,wq)
!     call gesqua (nquad,xq,yq,zq,wq)

!     write(6,'(''quadrature points'')')
!     do 30 i=1,nquad
!30     write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

      return
      end
!-----------------------------------------------------------------------
      subroutine getvps_fahy(rad,iel)
! Written by Claudia Filippi
! compute Fahy-pseudopotential for electron iel
      use atom_mod
      use const_mod
      use pseudo_mod
      use pseudo_fahy_mod
      implicit real*8(a-h,o-z)

      dimension rad(nelec,ncent)

      do 10 ic=1,ncent
        ict=iwctype(ic)
        r=rad(iel,ic)
! local potential
        if(r.lt.(npotl(ict)-1)*dradl(ict)) then
          ri=r/dradl(ict)
          ir=int(ri)
          ri=ri-dfloat(ir)
          ir=ir+1
          vps(iel,ic,lpotp1(ict))=potl(ir+1,ict)*ri+(1.d0-ri)*potl(ir,ict)
         else
          vps(iel,ic,lpotp1(ict))=-znuc(ict)/r
        endif
! non-local pseudopotential
        do 10 l=1,npotd(ict)
          if(l.ne.lpotp1(ict)) then
            if(r.lt.rcmax(ict)) then
              ri=r/drad(ict)
              ir=int(ri)
              ri=ri-dfloat(ir)
              ir=ir+1
              vps(iel,ic,l)=ptnlc(ir+1,ict,l)*ri+(1.d0-ri)*ptnlc(ir,ict,l)
             else
              vps(iel,ic,l)=0.0d0
            endif
          endif
   10 continue

      return
      end

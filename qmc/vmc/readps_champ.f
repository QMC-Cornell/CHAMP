      subroutine readps_champ
! Written by Cyrus Umrigar

! Read CHAMP-formatted nonlocal pseudopotentials.
! Reads also FHI-formatted nonlocal pseudopotentials, though I am not
! quite sure what some of the entries in the file are.
! Reads in v in Hartrees and subtracts out local part from all except
! the lpotp1 component.  If lpotp1<=0 then set lpotp1 so that the highest l
! channel is the local one.
! Also initializes quadrature pts.
! rmax_coul is the point at which the psp. becomes -Z/r to within eps.
! rmax_nloc is the point at which the psp. becomes local to within eps.
! For TM psps. rmax_nloc is considerably smaller than rmax_coul and so
! considerable computer time can be saved by evaluating the nonlocal
! components only for r < rmax_nloc rather than r < rmax_coul.

! Can use 3 different grids:
! igrid_ps=1, linear,              r(i)=r0_ps+(i-1)*h_ps
!         =2, exponential,         r(i)=r0_ps*exp((i-1)*h_ps)
!         =3, shifted exponential, r(i)=r0_ps*(exp((i-1)*h_ps)-1)
! The prefered grid is 3.
      use all_tools_mod
      use atom_mod
      use const_mod
      use pseudo_mod
      use qua_mod
      use pseudo_tm_mod
      implicit real*8(a-h,o-z)

      character*20 filename,atomtyp
      character*80 title

      dimension r(MPS_GRID),work(MPS_GRID)

      call alloc ('rmax_coul', rmax_coul, nctype)
      call alloc ('rmax_nloc', rmax_nloc, nctype)
      call alloc ('exp_h_ps', exp_h_ps, nctype)
      call alloc ('r0_ps', r0_ps, nctype)
      call alloc ('vpseudo', vpseudo, MPS_GRID, nctype, MPS_L)
      call alloc ('d2pot', d2pot, MPS_GRID, nctype, MPS_L)
      call alloc ('igrid_ps', igrid_ps, nctype)
      call alloc ('nr_ps', nr_ps, nctype)

      do 200 ict=1,nctype

        if(ict.lt.10) then
          write(atomtyp,'(i1)') ict
         elseif(ict.lt.100) then
          write(atomtyp,'(i2)') ict
         else
          stop 'readps_tm, nctype>100'
        endif

        filename='pseudopot'//atomtyp(1:index(atomtyp,' ')-1)
        open(1,file=filename,status='old',form='formatted',err=999)

        if(nloc.eq.4) then
          write(6,'(/3a)') ' Reading CHAMP format pseudopotential file >',trim(filename),'<'

! position file to skip an arbitrary number of comment lines but write out the first one
! if it exists
          title(1:1)='#'
          i=0
          do while(title(1:1).eq.'#')
            if(i.eq.1) write(6,'(x,a)') trim(title)
            i=i+1
            read(1,'(a80)') title
          enddo

! The TM psp. format has npotd and npotu for down and up, but we just use one of them
! They are the number of different l components of the psp.

          read(title,*) npotd(ict),zion,r_asymp
          write(6,'(a,i2,a,i2,a,f4.0,a,f8.3)') ' ict=',ict,', npotd(ict)=',npotd(ict),', zion=',zion, ', r_asymp=',r_asymp
          if(npotd(ict).le.0 .or. npotd(ict).gt.MPS_L) stop 'npotd must be > 0 and <= MPS_L'

          if(lpotp1(ict).gt.npotd(ict)) then
            write(6,'(''lpotp1(ict),npotd(ict)='',2i3)') lpotp1(ict),npotd(ict)
            stop 'Cannot choose local psp. to be > number of l components, lpotp1(ict) > npotd(ict)'
          endif
          if(npotd(ict).gt.MPS_L) stop 'npotd(ict).gt.MPS_L'

! If the local pseudopot component is not set in input, set it here
          if(lpotp1(ict).le.0) then
            lpotp1(ict)=npotd(ict)
            write(6,'('' center type'',i4,'' local pseudopot component reset to'',i3)') ict,lpotp1(ict)
          endif

          write(6,'('' center type'',i2,'' has'',i2,'' pseudopotential L components, and component''
     &    ,i2,'' is chosen to be local'')') ict,npotd(ict),lpotp1(ict)

          if(znuc(ict).ne.zion) then
            write(6,'(''znuc(ict) != zion in readps_tm'',2f6.1)') znuc(ict),zion
            stop 'znuc(ict) != zion in readps_tm'
          endif

          read(1,*) igrid_ps(ict),nr_ps(ict),r0_ps(ict),h_ps
          nr=nr_ps(ict)
          write(6,'(a,i2,a,i5,a,1pd22.15,a,0pf8.5)') ' igrid_ps(ict)=',igrid_ps(ict),', nr_ps(ict)='
     &    ,nr_ps(ict),', r0_ps(ict)=',r0_ps(ict),', h_ps=',h_ps
          exp_h_ps(ict)=exp(h_ps)

          if(igrid_ps(ict).lt.1 .or. igrid_ps(ict).gt.3) stop 'igrid_ps(ict) must be 1 or 2 or 3'
          if(igrid_ps(ict).lt.1 .and. r0_ps(ict).ne.0.d0) stop 'if igrid_ps(ict)=1 r0_ps(ict) must be 0'

          if(nr.lt.100) then
            write(6,'(''nr in psp grid too small'',2i6)') nr
            stop 'nr in psp grid too small'
          endif
          if(nr.gt.MPS_GRID .or. igrid_ps(ict).eq.2.and.nr.gt.MPS_GRID-1) then
            write(6,'(''nr > MPS_GRID'',2i6)') nr,MPS_GRID
            stop 'nr > MPS_GRID'
          endif

          if(igrid_ps(ict).eq.1 .or. igrid_ps(ict).eq.3) then
            nr=nr_ps(ict)
            do 10 ir=1,nr_ps(ict)
   10         read(1,*) r(ir),(vpseudo(ir,ict,i),i=1,npotd(ict))
            if(r(1).ne.0.d0) stop 'if igrid_ps is 1 or 3, r(1) must be 0'
           else
            nr_ps(ict)=nr_ps(ict)+1
            nr=nr_ps(ict)
            nrm1=nr-1
            do 20 ir=2,nr_ps(ict)
   20         read(1,*) r(ir),(vpseudo(ir,ict,i),i=1,npotd(ict))
            r(1)=0
            if(r0_ps(ict).le.0.d0 .or. h_ps.le.0.d0) then
              r0_ps(ict)=r(2)
              exp_h_ps(ict)=r(3)/r(2)
              h_ps=dlog(exp_h_ps(ict))
              write(6,'('' Grid parameters deduced from grid values are, r0_ps(ict),h_ps,exp_h_ps(ict)='',9f10.5)')
     &        r0_ps(ict),h_ps,exp_h_ps(ict)
            endif
            do 30 i=1,npotd(ict)
              call intpol(r(2),vpseudo(2,ict,i),nrm1,r(1),vpseudo(1,ict,i),1,3)
   30         write(6,'('' Interpolated psp'',9f16.12)') (vpseudo(ir,ict,i),ir=1,5)
          endif

         elseif(nloc.eq.5) then

          write(6,'(''Reading FHI format pseudopotential file '',a20)') filename
          read(1,'(a80)') title
          write(6,'(a)') trim(title)
          if(index(title,' GEN ').ne.0 .or. index(title,' Gen ').ne.0 .or. index(title,' gen ').ne.0) then
            write(6,'(''Use nloc=6 for GAMESS formatted pseudopotentials'')')
            stop 'Use nloc=6 for GAMESS formatted pseudopotentials'
          endif
! position file to skip items we do not use and do not know what they are
          read(1,*) crap,zion
          read(1,*) (junk,i=1,4),nr_ps(ict)
!         write(6,'(''nr_ps for atomtype'',i3,'' is'',i5)') ict,nr_ps(ict)
          do i=1,4
            read(1,*)
          enddo
          read(1,*) crap,npotd(ict)
          do i=1,10
            read(1,*)
          enddo
          igrid_ps(ict)=2
          r0_ps(ict)=0.d0
          h_ps=0.d0

! The TM psp. format has npotd and npotu for down and up, but we just use one of them
! They are the number of different l components of the psp.
          write(6,'(''ict,npotd(ict),zion'',2i2,f4.0)') ict,npotd(ict),zion
          if(npotd(ict).le.0 .or. npotd(ict).gt.MPS_L) stop 'npotd must be > 0 and <= MPS_L'

          if(lpotp1(ict).gt.npotd(ict)) then
            write(6,'(''lpotp1(ict),npotd(ict)='',2i3)') lpotp1(ict),npotd(ict)
            stop 'Cannot choose local psp. to be > number of l components, lpotp1(ict) > npotd(ict)'
          endif
          if(npotd(ict).gt.MPS_L) stop 'npotd(ict).gt.MPS_L'

! If the local pseudopot component is not set in input, set it here
          if(lpotp1(ict).le.0) then
            lpotp1(ict)=npotd(ict)
            write(6,'(''Center type'',i4,'' local pseudopot component reset to'',i3)') ict,lpotp1(ict)
          endif

          write(6,'(''Center type'',i2,'' has'',i2,'' pseudopotential L components, and component''
     &    ,i2,'' is chosen to be local'')') ict,npotd(ict),lpotp1(ict)

          if(znuc(ict).ne.zion) then
            write(6,'(''znuc(ict) != zion in readps_tm'',2f6.1)') znuc(ict),zion
            stop 'znuc(ict) != zion in readps_tm'
          endif

          nr=nr_ps(ict)
          write(6,'(''igrid_ps(ict),nr_ps(ict),r0_ps(ict),h_ps='',i2,i5,1pd22.15,0pf8.5)')
     &    igrid_ps(ict),nr_ps(ict),r0_ps(ict),h_ps
          exp_h_ps(ict)=exp(h_ps)

          if(igrid_ps(ict).lt.1 .or. igrid_ps(ict).gt.3) stop 'igrid_ps(ict) must be 1 or 2 or 3'
          if(igrid_ps(ict).lt.1 .and. r0_ps(ict).ne.0.d0) stop 'if igrid_ps(ict)=1 r0_ps(ict) must be 0'

          if(nr.lt.100) then
            write(6,'(''nr in psp grid too small'',2i6)') nr
            stop 'nr in psp grid too small'
          endif
          if(nr.gt.MPS_GRID .or. igrid_ps(ict).eq.2.and.nr.gt.MPS_GRID-1) then
            write(6,'(''nr > MPS_GRID'',2i6)') nr,MPS_GRID
            stop 'nr > MPS_GRID'
          endif

          nr_ps(ict)=nr_ps(ict)+1
          nr=nr_ps(ict)
          nrm1=nr-1
          do 40 i=1,npotd(ict)
            read(1,*)
            do 40 ir=2,nr_ps(ict)
   40         read(1,*) junk,r(ir),crap,vpseudo(ir,ict,i)
          r(1)=0
          if(r0_ps(ict).le.0.d0 .or. h_ps.le.0.d0) then
            r0_ps(ict)=r(2)
            exp_h_ps(ict)=r(3)/r(2)
            h_ps=dlog(exp_h_ps(ict))
            write(6,'(''Grid parameters deduced from grid values are, r0_ps(ict),h_ps,exp_h_ps(ict)='',9f10.5)')
     &      r0_ps(ict),h_ps,exp_h_ps(ict)
          endif
          do 45 i=1,npotd(ict)
            call intpol(r(2),vpseudo(2,ict,i),nrm1,r(1),vpseudo(1,ict,i),1,3)
   45     write(6,'(''Interpolated psp'',9f16.12)') (vpseudo(ir,ict,i),ir=1,5)

         else

          write(6,'(''readps_champ should only be called if nloc = 4 or 5'')')
          stop 'readps_champ should only be called if nloc = 4 or 5'

        endif

        if(r0_ps(ict).lt.0.d0 .or. r0_ps(ict).gt.1.d-2) stop 'r0_ps in psp grid is not reasonable'
        if(h_ps.le.0.d0 .or. h_ps.gt.1.d-1) stop 'h_ps in psp grid is not reasonable'

        close(1)

! Find the point beyond which r*v differs from zion by no more than .5*d-6.
! irmax_coul is used for the endpoint of the spline where it is assumed that
! the derivative of the local component is zion/r(irmax_coul)**2 and that of
! the local component is 0.  Also, rmax_coul is used in splfit_ps when it is
! called in a calculation of a periodic system.
        rmax_coul(ict)=0.d0
        irmax_coul=0
        do 50 ir=nr,1,-1
          do 50 i=1,npotd(ict)
            if(dabs(r(ir)*vpseudo(ir,ict,i)+zion).gt..5d-6) then
              rmax_coul(ict)=max(rmax_coul(ict),r(ir))
              irmax_coul=ir
              goto 60
            endif
   50   continue
   60   irmax_coul=min(irmax_coul+5,nr)

! Find the point beyond which the various v components differ from each other by no more than .5*d-6
        rmax_nloc(ict)=0.d0
        irmax_nloc=0
        do 70 ir=nr,1,-1
          do 70 i=2,npotd(ict)
            do 70 j=1,i-1
              if(dabs(vpseudo(ir,ict,i)-vpseudo(ir,ict,j)).gt..5d-6) then
                rmax_nloc(ict)=max(rmax_nloc(ict),r(ir))
                irmax_nloc=ir
                goto 80
              endif
   70   continue
   80   irmax_nloc=irmax_nloc+1
        rmax_nloc(ict)=r(irmax_nloc)

! rmax_nloc is used in getvps_champ to decide whether to calculate calculate nonloc part of psp
! or to set it to zero.  irmax_coul is used to decide how far out to spline the psp. components
! so irmax_coul must be >= irmax_nloc.
        irmax_coul=max(irmax_coul,irmax_nloc)
        rmax_coul(ict)=r(irmax_coul)

        write(6,'('' center '',i3,'' pseudopot rmax_coul,irmax_coul,rmax_nloc,irmax_nloc= '',2(f6.2,i5))')
     &  ict,rmax_coul(ict),irmax_coul,rmax_nloc(ict),irmax_nloc

        if(ipr.ge.1) then
          write(38,'(''r(j)  (vpseudo(j,ict,i),i=1,npotd(ict))  -znuc(ict)/r(j)'')')
          do 104 j=2,nr
            if(r(j).gt.0.d0) then
              write(38,'(1pd12.6,9d14.6)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict)),-znuc(ict)/r(j)
             else
              write(38,'(1pd12.6,9d14.6)') r(j),(vpseudo(j,ict,i),i=1,npotd(ict))
            endif
  104     continue
        endif

        do 110 i=1,npotd(ict)
          if(i.ne.lpotp1(ict)) then
            do 105 j=1,nr
  105         vpseudo(j,ict,i)=vpseudo(j,ict,i)-vpseudo(j,ict,lpotp1(ict))
          endif

  110   continue

        if(rmax_coul(ict).eq.0.d0) goto 200

        do 190 i=1,npotd(ict)

! Construct the spline

! Warning: the next line is correct only for pseudopotentials that have zero derivative at the origin.
! At present this routine is only used for such pseudopotentials, since for the Dolg pseudopotentials
! we call readps_gauss.f
          dpot1=0.d0

! Set derivative at end point equal to 0 for nonlocal components and Z/r^2 for local
          dpotn=0.d0
          if(i.eq.lpotp1(ict)) dpotn=zion/r(irmax_coul)**2
!         if(i.eq.lpotp1(ict)) call deriv_intpol(r,vpseudo(1,ict,i),nr,r(irmax_coul),dpotn,irmax_coul,3)

          write(6,'(a,1p1e15.5,a,1p1e15.5)') ' dpot1=',dpot1,',dpotn=',dpotn

! get second derivative for spline fit
          call spline2(r,vpseudo(1,ict,i),irmax_coul,dpot1,dpotn,d2pot(1,ict,i),work)

          do 190 j=1,nr
            if(ipr.ge.2) then
              if(i.eq.1) write(35,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
              if(i.eq.2) write(36,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
              if(i.eq.3) write(37,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
            endif
  190   continue

  200 continue

!     call gesqua(nquad,xq0,yq0,zq0,wq)
!     call gesqua(nquad,xq0,yq0,zq0,wq)
!     call gesqua(nquad,xq,yq,zq,wq)

!     write(6,'(''quadrature points'')')
!     do 210 i=1,nquad
! 210   write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

      return

  999 write(6,'(''Error: Pseudopot. file '',a20,'' is missing'')') filename
      stop 'Pseudopot. file is missing'

      end

      subroutine ie
c Written by Cyrus Umrigar
c Figure out iwbase, iwbasi, iworb, iebase, iebasi, ieorb
c for imposing symmetries while optimizing orbital parameters
c Make nparml too large and so has too many iworb's because it does not
c know about pivoting and cusp conditions yet.
c Also, needs to be corrected for more than one atom or atom type or for
c parameters that are negatives of each other
c necn -> necoef
      use all_tools_mod
      use atom_mod
      use coefs_mod
      use optim_mod
      use basis1_mod
      use numbas_mod
      use basis2_mod
      use lbas_mod
      implicit real*8(a-h,o-z)
      parameter(eps=1.d-5)

c      character*80 fmt


!      read(5,*) nctype,ncent
!JT      write(6,'(''nctype,ncent='',2i3)') nctype,ncent
!      read(5,*) (iwctype(icent),icent=1,ncent)
!JT      write(6,'(''iwctype='',20i3)') (iwctype(icent),icent=1,ncent)
!      read(5,*) numr
!JT      write(6,'(''numr='',i2)') numr
      if(numr.eq.1) then
!        read(5,*) (nrbas(ic),ic=1,nctype)
!JT        write(6,'(''nrbas='',20i3)') (nrbas(ic),ic=1,nctype)
      endif
!      read(5,*) norb,nbasis
!JT      write(6,'(''norb,nbasis='',9i5)') norb,nbasis
!      if(norb.gt.MORB) stop 'norb>MORB'

!      if(numr.eq.0 .or. numr.eq.1) then
!        call read_orb_loc_ana(lbasis)
!       elseif(numr.eq.-1) then
!        call read_orb_loc_ana2(lbasis)
!       else
!        stop 'numr must be between -1 and 1'
!      endif

!      write(6,'(''read in analytical orbitals'')')

!JT      write(6,'(''lbasis'',20a)') (lbasis(i),i=1,nbasis)

c I have to put in more features.
c At present I assume that if a coef is zero, then it is not to be varied and
c it need not be equal to other zero coefs.  Neither assumption is correct.
c 1st basis function is not varied because normalization not relevant
      nparml=0
      necn=0
      do 20 k=1,norb
        do 20 i=1,nbasis
          do 15 l=1,norb
            do 10 j=1,nbasis
              if((l-1)*nbasis+j .ge. (k-1)*nbasis+i) goto 15
              indexi=index(lbasis(i),' ',.true.)
              indexj=index(lbasis(j),' ',.true.)
              if(lbasis(i)(indexi+1:indexi+2).eq.lbasis(j)(indexj+1:indexj+2) .and. zex(i,1).eq.zex(j,1)) then
                if( ((i.ne.j).or.(k.ne.l)) .and. abs(coef(i,k,1)-coef(j,l,1)).le.eps .and. abs(coef(i,k,1)).gt.eps ) then
                  necn=necn+1
                  call alloc ('iebasi', iebasi, 2, necn)
                  call alloc ('ieorb', ieorb, 2, necn)
                  iebasi(1,necn)=i
                  iebasi(2,necn)=j
                  ieorb(1,necn)=k
                  ieorb(2,necn)=l
                  goto 20
                elseif( ((i.ne.j).or.(k.ne.l)) .and. abs(coef(i,k,1)+coef(j,l,1)).le.eps .and. abs(coef(i,k,1)).gt.eps ) then
                  necn=necn+1
                  call alloc ('iebasi', iebasi, 2, necn)
                  call alloc ('ieorb', ieorb, 2, necn)
                  iebasi(1,necn)=i
                  iebasi(2,necn)=j
                  ieorb(1,necn)=k
                  ieorb(2,necn)=-l
                  goto 20
                endif
              endif
   10       continue
   15     continue
!JT        if(abs(coef(i,k,1)).ge.eps) then
!JT          nparml=nparml+1
!JT          if(nparml.gt.MPARM) stop 'nparml > MPARM'
!JT          iwbasi(nparml)=i
!JT          iworb(nparml)=k
!JT        endif
   20 continue

!JT      write(6,'(2i4,'' nparml,necn'')') nparml,necn
!JT      write(fmt,'(''(''i3,''(2i3,x),a)'')') nparml
!JT      write(6,fmt) (iworb(i),iwbasi(i),i=1,nparml),' (iworb(i),iwbasi(i),i=1,nparml)'
!JT      write(fmt,'(''(''i3,''(2(2i3,x),x),a)'')') necn
!JT      write(6,fmt) ((ieorb(k,i),iebasi(k,i),k=1,2),i=1,necn),' ((ieorb(k,i),iebasi(k,i),k=1,2),i=1,necn)'

      nparme=1
      call alloc ('iwbase', iwbase, nparme)
      iwbase(1)=1
      nebase=0
      do 40 i=2,nbasis
        do 30 j=1,i-1
          indexi=index(lbasis(i),' ',.true.)
          indexj=index(lbasis(j),' ',.true.)
          if(lbasis(i)(indexi+1:indexi+2).eq.lbasis(j)(indexj+1:indexj+2) .and. zex(i,1).eq.zex(j,1)) then
            nebase=nebase+1
            call alloc ('iebase', iebase, 2, nebase)
            iebase(1,nebase)=i
            iebase(2,nebase)=j
            goto 40
          endif
   30   continue
        nparme=nparme+1
        call alloc ('iwbase', iwbase, nparme)
        iwbase(nparme)=i
   40 continue

!JT      write(6,'(2i4,'' nparme,nebase'')') nparme,nebase
!JT      write(fmt,'(''(''i3,''i3,a)'')') nparme
!JT      write(6,fmt) (iwbase(i),i=1,nparme),' (iwbase(i),i=1,nparme)'
!JT      write(fmt,'(''(''i3,''(2i3,x),a)'')') nebase
!JT      write(6,fmt) ((iebase(k,i),k=1,2),i=1,nebase),' ((iebase(k,i),k=1,2),i=1,nebase)'

c For imposing the cusp conditions figure out which orbs are s-like
c and encode that in the lo array.
c Warning: this is presently programmed only for certain molecules.
!JT    do 50 iorb=1,norb
!JT      if(coef(1,iorb,1).ne.0.d0) then
!JT        lo(iorb)=0
!JT     else
!JT       lo(iorb)=1
!JT     endif
!JT 50 continue
!JT   write(fmt,'(''(i1,''i3''i2,x,a)'')') norb-1
!JT   write(6,fmt) (lo(iorb),iorb=1,norb),'(lo(iorb),iorb=1,norb)'

      end

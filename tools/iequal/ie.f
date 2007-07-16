      program ie
c Written by Cyrus Umrigar
c Figure out iwbase, iwbasi, iworb, iebase, iebasi, ieorb
c for imposing symmetries while optimizing orbital parameters
c Make nparml too large and so has too many iworb's because it does not
c know about pivoting and cusp conditions yet.
c Also, needs to be corrected for more than one atom or atom type or for
c parameters that are negatives of each other
c necn -> necoef
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'fit.h'
      include 'numbas.h'
      include 'force.h'
      parameter(eps=1.d-5)

      character*80 fmt
      character*10 lbasis(MBASIS)

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),zex2(MRWF,MCTYPE,MWF),iwrwf2(MBASIS)
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

      common /optim/ lo(MORB),npoint(MORB),
     &iwjasa(83,NCTYP3X),iwjasb(83,3),iwjasc(83,MCTYPE),
     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
     &iedet(ICX,MDET),imnbas(MCENT),
     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
     &necn,nebase

      read(5,*) nctype,ncent
      write(6,'(''nctype,ncent='',2i3)') nctype,ncent
      if(nctype.gt.MCTYPE) stop 'ctype>MCTYPE'
      if(ncent.gt.MCENT) stop 'ncent>MCENT'
      read(5,*) (iwctype(icent),icent=1,ncent)
      write(6,'(''iwctype='',20i3)') (iwctype(icent),icent=1,ncent)
      read(5,*) numr
      write(6,'(''numr='',i2)') numr
      if(numr.eq.1) then
        read(5,*) (nrbas(ic),ic=1,nctype)
        write(6,'(''nrbas='',20i3)') (nrbas(ic),ic=1,nctype)
      endif
      read(5,*) norb,nbasis
      write(6,'(''norb,nbasis='',9i5)') norb,nbasis
      if(norb.gt.MORB) stop 'norb>MORB'
      if(nbasis.gt.MBASIS) stop 'nbasis>MBASIS'

      if(numr.eq.0 .or. numr.eq.1) then
        call read_orb_loc_ana(lbasis)
       elseif(numr.eq.-1) then
        call read_orb_loc_ana2(lbasis)
       else
        stop 'numr must be between -1 and 1'
      endif

      write(6,'(''read in analytical orbitals'')')

      write(6,'(''lbasis'',20a)') (lbasis(i),i=1,nbasis)
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
                  iebasi(1,necn)=i
                  iebasi(2,necn)=j
                  ieorb(1,necn)=k
                  ieorb(2,necn)=l
                  goto 20
                elseif( ((i.ne.j).or.(k.ne.l)) .and. abs(coef(i,k,1)+coef(j,l,1)).le.eps .and. abs(coef(i,k,1)).gt.eps ) then
                  necn=necn+1
                  iebasi(1,necn)=i
                  iebasi(2,necn)=j
                  ieorb(1,necn)=k
                  ieorb(2,necn)=-l
                  goto 20
                endif
              endif
   10       continue
   15     continue
          if(abs(coef(i,k,1)).ge.eps) then
            nparml=nparml+1
            if(nparml.gt.MPARM) stop 'nparml > MPARM'
            iwbasi(nparml)=i
            iworb(nparml)=k
          endif
   20 continue

      write(6,'(2i4,'' nparml,necn'')') nparml,necn
      write(fmt,'(''(''i3,''(2i3,x),a)'')') nparml
      write(6,fmt) (iworb(i),iwbasi(i),i=1,nparml),' (iworb(i),iwbasi(i),i=1,nparml)'
      write(fmt,'(''(''i3,''(2(2i3,x),x),a)'')') necn
      write(6,fmt) ((ieorb(k,i),iebasi(k,i),k=1,2),i=1,necn),' ((ieorb(k,i),iebasi(k,i),k=1,2),i=1,necn)'

      nparme=1
      iwbase(1)=1
      nebase=0
      do 40 i=2,nbasis
        do 30 j=1,i-1
          indexi=index(lbasis(i),' ',.true.)
          indexj=index(lbasis(j),' ',.true.)
          if(lbasis(i)(indexi+1:indexi+2).eq.lbasis(j)(indexj+1:indexj+2) .and. zex(i,1).eq.zex(j,1)) then
            nebase=nebase+1
            iebase(1,nebase)=i
            iebase(2,nebase)=j
            goto 40
          endif
   30   continue
        nparme=nparme+1
        iwbase(nparme)=i
   40 continue

      write(6,'(2i4,'' nparme,nebase'')') nparme,nebase
      write(fmt,'(''(''i3,''i3,a)'')') nparme
      write(6,fmt) (iwbase(i),i=1,nparme),' (iwbase(i),i=1,nparme)'
      write(fmt,'(''(''i3,''(2i3,x),a)'')') nebase
      write(6,fmt) ((iebase(k,i),k=1,2),i=1,nebase),' ((iebase(k,i),k=1,2),i=1,nebase)'

c For imposing the cusp conditions figure out which orbs are s-like
c and encode that in the lo array.
c Warning: this is presently programmed only for certain molecules.
      do 50 iorb=1,norb
        if(coef(1,iorb,1).ne.0.d0) then
          lo(iorb)=0
         else
          lo(iorb)=1
        endif
   50 continue
      write(fmt,'(''(i1,''i3''i2,x,a)'')') norb-1
      write(6,fmt) (lo(iorb),iorb=1,norb),'(lo(iorb),iorb=1,norb)'

      stop
      end
c-----------------------------------------------------------------------

c read_orb_loc_ana and read_orb_loc_ana2 are the same as in vmc/read_orb_ana.f
c except that the subroutine argument (lbasis) has been added in.

      subroutine read_orb_loc_ana(lbasis)
c Written by Cyrus Umrigar
c Reads in localized orbitals, in one or more of
c 1) a slater basis
c 2) an asymptotic basis
c 3) a gaussian basis
c For the p and d basis functions, I maintain the old order to avoid
c having to change old inputs.
c The order I read in the p functions is x,y,z, which is m=1,-1,0 and
c the order that I read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, which is m=0,2,-2,1,-1,
c in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.

      implicit real*8(a-h,o-z)

      parameter(nprime=5)

      include 'vmc.h'
c     include 'ewald.h'
c     include 'pseudo.h'
      include 'numbas.h'
      include 'force.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),zex2(MRWF,MCTYPE,MWF),iwrwf2(MBASIS)

      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

      character*6 l1s(nprime),l2s(nprime),l2p(nprime,-1:1)
     &,l3s(nprime),l3p(nprime,-1:1),l3d(nprime,-2:2)
c    &,l4s(nprime),l4p(nprime,-1:1),l4d(nprime,-2:2)
     &,l4s(nprime),l4p(nprime,-1:1)
     &,lsa,lpa(-1:1),lda(-2:2)
      character*2 lmbas
      character*4 lcent,l4f(nprime),l5g(nprime),l6h(nprime)
      character*10 lbasis(MBASIS)

      data l1s /'    1s','   1s''','   1s"','  1s"''','  1s""'/
      data l2s /'    2s','   2s''','   2s"','  2s"''','  2s""'/
      data l3s /'    3s','   3s''','   3s"','  3s"''','  3s""'/
      data l4s /'    4s','   4s''','   4s"','  4s"''','  4s""'/
      data l2p /'   2py','  2py''','  2py"',' 2py"''',' 2py""'
     &         ,'   2pz','  2pz''','  2pz"',' 2pz"''',' 2pz""'
     &         ,'   2px','  2px''','  2px"',' 2px"''',' 2px""'/
      data l3p /'   3py','  3py''','  3py"',' 3py"''',' 3py""'
     &         ,'   3pz','  3pz''','  3pz"',' 3pz"''',' 3pz""'
     &         ,'   3px','  3px''','  3px"',' 3px"''',' 3px""'/
      data l4p /'   4py','  4py''','  4py"',' 4py"''',' 4py""'
     &         ,'   4pz','  4pz''','  4pz"',' 4pz"''',' 4pz""'
     &         ,'   4px','  4px''','  4px"',' 4px"''',' 4px""'/
      data l3d /'  3dxy',' 3dxy''',' 3dxy"','3dxy"''','3dxy""'
     &         ,'  3dyz',' 3dyz''',' 3dyz"','3dyz"''','3dyz""'
     &         ,'  3dzr',' 3dzr''',' 3dzr"','3dzr"''','3dzr""'
     &         ,'  3dxz',' 3dxz''',' 3dxz"','3dxz"''','3dxz""'
     &         ,'  3dx2',' 3dx2''',' 3dx2"','3dx2"''','3dx2""'/
c     data l4d /'  4dxy',' 4dxy''',' 4dxy"','4dxy"''','4dxy""'
c    &         ,'  4dyz',' 4dyz''',' 4dyz"','4dyz"''','4dyz""'
c    &         ,'  4dzr',' 4dzr''',' 4dzr"','4dzr"''','4dzr""'
c    &         ,'  4dxz',' 4dxz''',' 4dxz"','4dxz"''','4dxz""'
c    &         ,'  4dx2',' 4dx2''',' 4dx2"','4dx2"''','4dx2""'/
      data l4f /' f','f''','f"','f"''','f""'/
      data l5g /' g','g''','g"','g"''','g""'/
      data l6h /' h','h''','h"','h"''','h""'/
      data lsa /'    sa'/
      data lpa /'   pya','   pza','   pxa'/
      data lda /'   dxy','   dyz','   dzr','   dxz','   dx2'/

c Note: the order that we read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.
      if(numr.le.0) then
        do 10 i=1,nctype
   10     read(5,*) n1s(i)
     &    ,n2s(i),n2p(1,i),n2p(-1,i),n2p(0,i)
     &    ,n3s(i),n3p(1,i),n3p(-1,i),n3p(0,i)
     &    ,n3d(0,i),(n3d(m,i),n3d(-m,i),m=2,1,-1)
     &    ,n4s(i),n4p(1,i),n4p(-1,i),n4p(0,i)
c    &    ,n4d(0,i),(n4d(m,i),n4d(-m,i),m=2,1,-1)
     &    ,nsa(i),npa(1,i),npa(-1,i),npa(0,i)
     &    ,nda(0,i),(nda(m,i),nda(-m,i),m=2,1,-1)

        nbas_tot=0
        do 24 ic=1,ncent
          ict=iwctype(ic)
          nbas_typ=           iabs(n1s(ict))+iabs(n2s(ict))+iabs(n3s(ict))+iabs(n4s(ict))
          nbas_tot=nbas_tot + iabs(n1s(ict))+iabs(n2s(ict))+iabs(n3s(ict))+iabs(n4s(ict))
          do 21 m=-1,1
            nbas_typ=nbas_typ + iabs(n2p(m,ict)) + iabs(n3p(m,ict)) + iabs(n4p(m,ict))
   21       nbas_tot=nbas_tot + iabs(n2p(m,ict)) + iabs(n3p(m,ict)) + iabs(n4p(m,ict))
          do 22 m=-2,2
            nbas_typ=nbas_typ + iabs(n3d(m,ict)) + iabs(n4d(m,ict))
   22       nbas_tot=nbas_tot + iabs(n3d(m,ict)) + iabs(n4d(m,ict))
          do 23 m=-3,3
            nbas_typ=nbas_typ + iabs(n4f(m,ict))
   23       nbas_tot=nbas_tot + iabs(n4f(m,ict))
   24     nbasis_ctype(ict)=nbas_typ
        if(nbas_tot.ne.nbasis) then
          write(6,'(''nbas_tot,nbasis='',9i4)') nbas_tot,nbasis
          stop 'nbas_tot not equal to nbasis'
        endif

        betaq=0
        do 25 ic=1,ncent
   25     betaq=betaq+znuc(iwctype(ic))
        betaq=betaq-nelec+1

       else

        do 49 ict=1,nctype
          write(6,'(/,''center type'',i3)') ict
          read(5,*) n1s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &    ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &    ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &    ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &    ,(n6h(m,ict),n6h(-m,ict),m=5,1,-1),n6h(0,ict)
          write(6,'(''number of basis fns. for each l:'',100i2)')
     &    n1s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &    ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &    ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &    ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &    ,(n6h(m,ict),n6h(-m,ict),m=5,1,-1),n6h(0,ict)

          nbas=iabs(n1s(ict))
          n2s(ict)=0
          n3s(ict)=0
          n4s(ict)=0
          do 41 m=-1,1
            n3p(m,ict)=0
            n4p(m,ict)=0
            if(n2p(m,ict).ne.0 .and. ML_BAS.lt.1) stop 'ML_BAS in vmc.h too small'
   41       nbas=nbas+ iabs(n2p(m,ict))
          do 42 m=-2,2
            n4d(m,ict)=0
            if(n3d(m,ict).ne.0 .and. ML_BAS.lt.2) stop 'ML_BAS in vmc.h too small'
   42       nbas=nbas+ iabs(n3d(m,ict))
          do 43 m=-3,3
            if(n4f(m,ict).ne.0 .and. ML_BAS.lt.3) stop 'ML_BAS in vmc.h too small'
   43       nbas=nbas+ iabs(n4f(m,ict))
          do 44 m=-4,4
            if(n5g(m,ict).ne.0 .and. ML_BAS.lt.4) stop 'ML_BAS in vmc.h too small'
   44       nbas=nbas+ iabs(n5g(m,ict))
          do 45 m=-5,5
            if(n6h(m,ict).ne.0 .and. ML_BAS.lt.5) stop 'ML_BAS in vmc.h too small'
   45       nbas=nbas+ iabs(n6h(m,ict))
          if(nbas.gt.MBASIS_CTYPE) stop 'nbas > MBASIS_CTYPE'
          nbasis_ctype(ict)=nbas

          write(6,'(''Reading iwrwf'')')
          read(5,*) (iwrwf(ib,ict),ib=1,nbas)
          write(6,'(''Center'',i5,'' uses radial bas. fns:'',(100i3))')
     &    ict,(iwrwf(ib,ict),ib=1,nbas)

          do 46 ib=1,nbas
            if(iwrwf(ib,ict).gt.nrbas(ict)) then
              write(6,'(''ict,ib,iwrwf(ib,ict),nrbas(ict)'',9i3)')
     &        ict,ib,iwrwf(ib,ict),nrbas(ict)
              stop 'iwrwf(ib,ict) > nrbas(ict)'
            endif
   46     continue

          nrwf=1
          do 48 ib=2,nbas
            do 47 ib2=1,ib
   47         if(iwrwf(ib2,ict).eq.iwrwf(ib,ict)) goto 48
            nrwf=nrwf+1
   48     continue
   49     if(nrwf.gt.MRWF) stop 'nrwf > MRWF'

        nbas_tot=0
        do 55 ic=1,ncent
          ict=iwctype(ic)
          nbas_tot=nbas_tot+ iabs(n1s(ict))
          do 51 m=-1,1
   51       nbas_tot=nbas_tot+ iabs(n2p(m,ict))
          do 52 m=-2,2
   52       nbas_tot=nbas_tot+ iabs(n3d(m,ict))
          do 53 m=-3,3
   53       nbas_tot=nbas_tot+ iabs(n4f(m,ict))
          do 54 m=-4,4
   54       nbas_tot=nbas_tot+ iabs(n5g(m,ict))
          do 55 m=-5,5
   55       nbas_tot=nbas_tot+ iabs(n6h(m,ict))
        if(nbas_tot.ne.nbasis) stop 'nbas_tot not equal to nbasis'

      endif

      write(6,'(/,''center type'',(12i4))') (i,i=1,nctype)
      write(6,'(/,''1s'',t11,(12i4))') (n1s(i),i=1,nctype)
      if(numr.le.0) write(6,'(''2s'',t11,(12i4))') (n2s(i),i=1,nctype)
      write(6,'(''2px'',t11,(12i4))') (n2p(1,i),i=1,nctype)
      write(6,'(''2py'',t11,(12i4))') (n2p(-1,i),i=1,nctype)
      write(6,'(''2pz'',t11,(12i4))') (n2p(0,i),i=1,nctype)
      if(numr.le.0) then
        write(6,'(''3s'',t11,(12i4))') (n3s(i),i=1,nctype)
        write(6,'(''3px'',t11,(12i4))') (n3p(1,i),i=1,nctype)
        write(6,'(''3py'',t11,(12i4))') (n3p(-1,i),i=1,nctype)
        write(6,'(''3pz'',t11,(12i4))') (n3p(0,i),i=1,nctype)
      endif
      write(6,'(''3dzr'',t11,(12i4))') (n3d(0,i),i=1,nctype)
      write(6,'(''3dx2'',t11,(12i4))') (n3d(2,i),i=1,nctype)
      write(6,'(''3dxy'',t11,(12i4))') (n3d(-2,i),i=1,nctype)
      write(6,'(''3dxz'',t11,(12i4))') (n3d(1,i),i=1,nctype)
      write(6,'(''3dyz'',t11,(12i4))') (n3d(-1,i),i=1,nctype)
      if(numr.le.0) then
        write(6,'(''4s'',t11,(12i4))') (n4s(i),i=1,nctype)
        write(6,'(''4px'',t11,(12i4))') (n4p(1,i),i=1,nctype)
        write(6,'(''4py'',t11,(12i4))') (n4p(-1,i),i=1,nctype)
        write(6,'(''4pz'',t11,(12i4))') (n4p(0,i),i=1,nctype)
c       write(6,'(''4dzr'',t11,(12i4))') (n4d(0,i),i=1,nctype)
c       write(6,'(''4dx2'',t11,(12i4))') (n4d(2,i),i=1,nctype)
c       write(6,'(''4dxy'',t11,(12i4))') (n4d(-2,i),i=1,nctype)
c       write(6,'(''4dxz'',t11,(12i4))') (n4d(1,i),i=1,nctype)
c       write(6,'(''4dyz'',t11,(12i4))') (n4d(-1,i),i=1,nctype)
        write(6,'(''sa'',t11,(12i4))') (nsa(i),i=1,nctype)
        write(6,'(''pxa'',t11,(12i4))') (npa(1,i),i=1,nctype)
        write(6,'(''pya'',t11,(12i4))') (npa(-1,i),i=1,nctype)
        write(6,'(''pza'',t11,(12i4))') (npa(0,i),i=1,nctype)
        write(6,'(''dzra'',t11,(12i4))') (nda(0,i),i=1,nctype)
        write(6,'(''dx2a'',t11,(12i4))') (nda(2,i),i=1,nctype)
        write(6,'(''dxya'',t11,(12i4))') (nda(-2,i),i=1,nctype)
        write(6,'(''dxza'',t11,(12i4))') (nda(1,i),i=1,nctype)
        write(6,'(''dyza'',t11,(12i4))') (nda(-1,i),i=1,nctype)
       else
        do 60 m=-3,3
   60     write(6,'(''4f('',i2,'')'',t11,(12i4))') m,(n4f(m,i),i=1,nctype)
        do 70 m=-4,4
   70     write(6,'(''5g('',i2,'')'',t11,(12i4))') m,(n5g(m,i),i=1,nctype)
        do 80 m=-5,5
   80     write(6,'(''6h('',i2,'')'',t11,(12i4))') m,(n6h(m,i),i=1,nctype)
      endif
      write(6,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=1,nctype)
      write(6,*)

c the following sets up character strings that identify basis
c functions on centers for the printout of coefficients and
c screening constants. Character strings are of form
c 2px'(3) for the second 2px function on the third center

      ib=0
      do 250 ic=1,ncent
        i=iwctype(ic)

        write(lcent,'(''('',i2,'')'')') ic
        if(ic.lt.10)  write(lcent,'(''(0'',i1,'')'')') ic

        do 90 j=1,iabs(n1s(i))
c         if(iabs(n1s(i)).gt.nprime) stop 'number of 1s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(1,n1s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
   90     lbasis(ib)=l1s(min(j,nprime))//lcent

        do 100 j=1,iabs(n2s(i))
c         if(iabs(n2s(i)).gt.nprime) stop 'number of 2s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(2,n2s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  100     lbasis(ib)=l2s(min(j,nprime))//lcent

        do 110 m=1,0,-1
          do 105 j=1,iabs(n2p(m,i))
c           if(iabs(n2p(m,i)).gt.nprime) stop 'number of 2p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(2,n2p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  105       lbasis(ib)=l2p(min(j,nprime),m)//lcent
          if(m.eq.0) goto 110
          do 108 j=1,iabs(n2p(-m,i))
c           if(iabs(n2p(-m,i)).gt.nprime) stop 'number of 2p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(2,n2p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  108       lbasis(ib)=l2p(min(j,nprime),-m)//lcent
  110   continue

        do 120 j=1,iabs(n3s(i))
c         if(iabs(n3s(i)).gt.nprime) stop 'number of 3s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(3,n3s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  120     lbasis(ib)=l3s(min(j,nprime))//lcent

        do 130 m=1,0,-1
          do 125 j=1,iabs(n3p(m,i))
c           if(iabs(n3p(m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  125       lbasis(ib)=l3p(min(j,nprime),m)//lcent
          if(m.eq.0) goto 130
          do 128 j=1,iabs(n3p(-m,i))
c           if(iabs(n3p(-m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  128       lbasis(ib)=l3p(min(j,nprime),-m)//lcent
  130   continue

          m=0
          do 132 j=1,iabs(n3d(m,i))
c           if(iabs(n3d(m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  132       lbasis(ib)=l3d(min(j,nprime),m)//lcent
        do 140 m=2,1,-1
          do 135 j=1,iabs(n3d(m,i))
c           if(iabs(n3d(m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  135       lbasis(ib)=l3d(min(j,nprime),m)//lcent
          do 140 j=1,iabs(n3d(-m,i))
c           if(iabs(n3d(-m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(-m,i))
            l_bas(ib)=2
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  140       lbasis(ib)=l3d(min(j,nprime),-m)//lcent

        do 150 j=1,iabs(n4s(i))
c         if(iabs(n4s(i)).gt.nprime) stop 'number of 4s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(4,n4s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          ib=ib+1
  150     lbasis(ib)=l4s(min(j,nprime))//lcent

        do 160 m=1,0,-1
          do 155 j=1,iabs(n4p(m,i))
c           if(iabs(n4p(m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  155       lbasis(ib)=l4p(min(j,nprime),m)//lcent
          if(m.eq.0) goto 160
          do 158 j=1,iabs(n4p(-m,i))
c           if(iabs(n4p(-m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  158       lbasis(ib)=l4p(min(j,nprime),-m)//lcent
  160   continue

c         m=0
c         do 162 j=1,iabs(n4d(m,i))
c           if(iabs(n4d(m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
c           ib=ib+1
c           n_bas(ib)=sign(4,n4d(m,i))
c           l_bas(ib)=2
c           m_bas(ib)=m
c           icenter_basis(ib)=ic
c           ictype_basis(ib)=i
c 162       lbasis(ib)=l4d(min(j,nprime),m)//lcent
c       do 170 m=2,1,-1
c         do 165 j=1,iabs(n4d(m,i))
c           if(iabs(n4d(m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
c           ib=ib+1
c           n_bas(ib)=sign(4,n4d(m,i))
c           l_bas(ib)=2
c           m_bas(ib)=m
c           icenter_basis(ib)=ic
c           ictype_basis(ib)=i
c 165       lbasis(ib)=l4d(min(j,nprime),m)//lcent
c         do 170 j=1,iabs(n4d(-m,i))
c           if(iabs(n4d(-m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
c           ib=ib+1
c           n_bas(ib)=sign(4,n4d(-m,i))
c           l_bas(ib)=2
c           m_bas(ib)=-m
c           icenter_basis(ib)=ic
c           ictype_basis(ib)=i
c 170       lbasis(ib)=l4d(min(j,nprime),-m)//lcent

c Ordering to be made consistent
        do 180 m=-3,3
          do 180 j=1,iabs(n4f(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 4f basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4f(m,i))
            l_bas(ib)=3
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
c           write(lmbas,'(''('',i2,'')'')') m
            write(lmbas,'(i2)') m
  180       lbasis(ib)=l4f(min(j,nprime))//lmbas//lcent

        do 190 m=-4,4
          do 190 j=1,iabs(n5g(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 5g basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5g(m,i))
            l_bas(ib)=4
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') m
  190       lbasis(ib)=l5g(min(j,nprime))//lmbas//lcent

        do 200 m=-5,5
          do 200 j=1,iabs(n6h(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 6h basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(6,n6h(m,i))
            l_bas(ib)=5
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') m
  200       lbasis(ib)=l6h(min(j,nprime))//lmbas//lcent

      if(iabs(nsa(i)).ge.1) then
        ib=ib+1
        n_bas(ib)=0
        l_bas(ib)=0
        m_bas(ib)=0
        icenter_basis(ib)=ic
        ictype_basis(ib)=i
        lbasis(ib)=lsa//lcent
      endif

      do 220 m=-1,1
        if(iabs(npa(m,i)).ge.1) then
          ib=ib+1
          n_bas(ib)=0
          l_bas(ib)=1
          m_bas(ib)=m
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          lbasis(ib)=lpa(m)//lcent
        endif
  220 continue

      do 230 m=-2,2
        if(iabs(nda(m,i)).ge.1) then
          ib=ib+1
          n_bas(ib)=0
          l_bas(ib)=2
          m_bas(ib)=m
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          lbasis(ib)=lda(m)//lcent
        endif
  230 continue

  250 continue

      write(6,'(/,(12a10))') (lbasis(j),j=1,nbasis)
      write(6,'(''orbital coefficients'')')
      do 260 i=1,norb
        read(5,*) (coef(j,i,1),j=1,nbasis)
  260   write(6,'(12f10.6)') (coef(j,i,1),j=1,nbasis)

c Warning: Foll. "if" commented out till I modify change_input to work without the zex
c     if(inumr.ge.1) then
        write(6,'(''screening constants'')')
        read(5,*) (zex(i,1),i=1,nbasis)
        write(6,'(12f10.6)') (zex(i,1),i=1,nbasis)
        do 262 i=1,nbasis
  262     if(zex(i,1).le.0.d0) stop 'exponent zex should be > 0'
c     endif

c If an analytical basis is used find out how many distinct radial orbitals there are
c and store n and zex for each of these.
c This is done so that only the unique analytic radial wavefns. are evaluated in basis_fns.f and basis_fnse.f
c zex2 is used in func when optimizing the exponents.
c ib  goes up to the total number of basis functions
c ib2 goes up to the number of basis functions on each center type
c irb goes up to the number of radial basis functions on each center type
      if(numr.le.0) then
        write(6,'(''n_bas(ib),l_bas(ib)'',50(i2,x))') (n_bas(ib),l_bas(ib),ib=1,nbasis)
        ib=0
        ib4=0
        do 290 ic=1,ncent
          ict=iwctype(ic)
          write(6,'(''n_bas(ib),l_bas(ib)'',50(i2,x))') (n_bas(i),l_bas(i),i=1,nbasis_ctype(ict))
          nrbas(ict)=0
          do 280 ib2=1,nbasis_ctype(ict)
            ib=ib+1
            if(ib.ne.ib2+ib4) stop 'ib .ne. ib2+ib4'
            do 270 ib3=1,ib2-1
              if(n_bas(ib).eq.n_bas(ib3+ib4)
     &        .and. zex(ib,1).eq.zex(ib3+ib4,1)) then
                irb=iwrwf2(ib3+ib4)
                goto 275
              endif
  270       continue
            nrbas(ict)=nrbas(ict)+1
            irb=nrbas(ict)
            if(irb.gt.MRWF) stop 'nbas > MRWF'
            n_bas2(irb,ict)=n_bas(ib)
            zex2(irb,ict,1)=zex(ib,1)
  275       iwrwf2(ib)=irb
  280       iwrwf(ib2,ict)=irb
  290       ib4=ib4+nbasis_ctype(ict)
        write(6,'(''nrbas(ict)='',40i4)') (nrbas(ict),ict=1,nctype)
        write(6,'(''iwrwf2(ib)='',40i3)') (iwrwf2(ib),ib=1,nbasis)
        write(6,'(''iwrwf(ib2,ict)='',40i3)') ((iwrwf(ib2,ict),ib2=1,nbasis_ctype(ict)),ict=1,nctype)
      endif

c     if(numr.le.0) then
c       ib=0
c       ib4=0
c       do 290 ict=1,nctype
c         irb=0
c         do 290 ib2=1,nbasis_ctype(ict)
c           ib=ib+1
c           do 270 ib3=1,ib2-1
c             write(6,'(''ict,ib2,ib3,ib4,ib='',9i5)') ict,ib2,ib3,ib4,ib
c 270         if(n_bas(ib)-l_bas(ib).eq.n_bas(ib3+ib4)-l_bas(ib3+ib4)
c    &        .and. zex(ib,1).eq.zex(ib3+ib4,1)) goto 280
c           irb=irb+1
c           n_bas2(irb,ict)=n_bas(ib)
c           zex2(irb,ict,1)=zex(ib,1)
c           iwrwf2(ib)=irb
c 280       iwrwf(ib2,ict)=irb
c           if(irb.gt.MRWF) stop 'nbas > MRWF'
c           nrbas(ict)=irb
c 290       ib4=ib4+nbasis_ctype(ict)
c     endif

      return
      end
c-----------------------------------------------------------------------

      subroutine read_orb_loc_ana2(lbasis)
c Written by Cyrus Umrigar
c Reads in localized orbitals, in one or more of
c 1) a slater basis
c 2) an asymptotic basis
c 3) a gaussian basis
c For the p and d basis functions, I maintain the old order to avoid
c having to change old inputs.
c The order I read in the p functions is x,y,z, which is m=1,-1,0 and
c the order that I read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, which is m=0,2,-2,1,-1,
c in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.
c This routine differs from read_orb_loc_ana in that the order of the fns.
c is changed from
c 1s, 2s, 2px, 2py, ...          to
c 1s, 2s, 3s, 4s,  2px, 2py, 2pz, 3px, ... ., 4pz,  3d, 4f, 5g, 6h

      implicit real*8(a-h,o-z)

      parameter(nprime=5)

      include 'vmc.h'
c     include 'ewald.h'
c     include 'pseudo.h'
      include 'numbas.h'
      include 'force.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE),n3s(MCTYPE),n3p(-1:1,MCTYPE)
     &,n3d(-2:2,MCTYPE),n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
     &,n4f(-3:3,MCTYPE),n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /basis2/ n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),zex2(MRWF,MCTYPE,MWF),iwrwf2(MBASIS)

      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)

      character*6 l1s(nprime),l2s(nprime),l2p(nprime,-1:1)
     &,l3s(nprime),l3p(nprime,-1:1),l3d(nprime,-2:2)
c    &,l4s(nprime),l4p(nprime,-1:1),l4d(nprime,-2:2)
     &,l4s(nprime),l4p(nprime,-1:1)
     &,lsa,lpa(-1:1),lda(-2:2)
      character*2 lmbas
      character*4 lcent,l4f(nprime),l5g(nprime),l6h(nprime)
      character*10 lbasis(MBASIS)

      data l1s /'    1s','   1s''','   1s"','  1s"''','  1s""'/
      data l2s /'    2s','   2s''','   2s"','  2s"''','  2s""'/
      data l3s /'    3s','   3s''','   3s"','  3s"''','  3s""'/
      data l4s /'    4s','   4s''','   4s"','  4s"''','  4s""'/
      data l2p /'   2py','  2py''','  2py"',' 2py"''',' 2py""'
     &         ,'   2pz','  2pz''','  2pz"',' 2pz"''',' 2pz""'
     &         ,'   2px','  2px''','  2px"',' 2px"''',' 2px""'/
      data l3p /'   3py','  3py''','  3py"',' 3py"''',' 3py""'
     &         ,'   3pz','  3pz''','  3pz"',' 3pz"''',' 3pz""'
     &         ,'   3px','  3px''','  3px"',' 3px"''',' 3px""'/
      data l4p /'   4py','  4py''','  4py"',' 4py"''',' 4py""'
     &         ,'   4pz','  4pz''','  4pz"',' 4pz"''',' 4pz""'
     &         ,'   4px','  4px''','  4px"',' 4px"''',' 4px""'/
      data l3d /'  3dxy',' 3dxy''',' 3dxy"','3dxy"''','3dxy""'
     &         ,'  3dyz',' 3dyz''',' 3dyz"','3dyz"''','3dyz""'
     &         ,'  3dzr',' 3dzr''',' 3dzr"','3dzr"''','3dzr""'
     &         ,'  3dxz',' 3dxz''',' 3dxz"','3dxz"''','3dxz""'
     &         ,'  3dx2',' 3dx2''',' 3dx2"','3dx2"''','3dx2""'/
c     data l4d /'  4dxy',' 4dxy''',' 4dxy"','4dxy"''','4dxy""'
c    &         ,'  4dyz',' 4dyz''',' 4dyz"','4dyz"''','4dyz""'
c    &         ,'  4dzr',' 4dzr''',' 4dzr"','4dzr"''','4dzr""'
c    &         ,'  4dxz',' 4dxz''',' 4dxz"','4dxz"''','4dxz""'
c    &         ,'  4dx2',' 4dx2''',' 4dx2"','4dx2"''','4dx2""'/
      data l4f /' f','f''','f"','f"''','f""'/
      data l5g /' g','g''','g"','g"''','g""'/
      data l6h /' h','h''','h"','h"''','h""'/
      data lsa /'    sa'/
      data lpa /'   pya','   pza','   pxa'/
      data lda /'   dxy','   dyz','   dzr','   dxz','   dx2'/

c Note: the order that we read in the d functions is
c 3z^-r^2, x^2-y^2, xy, xz, yz, in order to be able to use old inputs.
c All others are read in the foll. order: l, -l, l-1, -(l-1), ..., 0,
c i.e. the order in which we were reading in the p functions.
      if(numr.le.0) then
        do 10 i=1,nctype
   10     read(5,*) n1s(i)
     &    ,n2s(i),n2p(1,i),n2p(-1,i),n2p(0,i)
     &    ,n3s(i),n3p(1,i),n3p(-1,i),n3p(0,i)
     &    ,n3d(0,i),(n3d(m,i),n3d(-m,i),m=2,1,-1)
     &    ,n4s(i),n4p(1,i),n4p(-1,i),n4p(0,i)
c    &    ,n4d(0,i),(n4d(m,i),n4d(-m,i),m=2,1,-1)
     &    ,nsa(i),npa(1,i),npa(-1,i),npa(0,i)
     &    ,nda(0,i),(nda(m,i),nda(-m,i),m=2,1,-1)

        nbas_tot=0
        do 24 ic=1,ncent
          ict=iwctype(ic)
          nbas_typ=           iabs(n1s(ict))+iabs(n2s(ict))+iabs(n3s(ict))+iabs(n4s(ict))
          nbas_tot=nbas_tot + iabs(n1s(ict))+iabs(n2s(ict))+iabs(n3s(ict))+iabs(n4s(ict))
          do 21 m=-1,1
            nbas_typ=nbas_typ + iabs(n2p(m,ict)) + iabs(n3p(m,ict)) + iabs(n4p(m,ict))
   21       nbas_tot=nbas_tot + iabs(n2p(m,ict)) + iabs(n3p(m,ict)) + iabs(n4p(m,ict))
          do 22 m=-2,2
            nbas_typ=nbas_typ + iabs(n3d(m,ict)) + iabs(n4d(m,ict))
   22       nbas_tot=nbas_tot + iabs(n3d(m,ict)) + iabs(n4d(m,ict))
          do 23 m=-3,3
            nbas_typ=nbas_typ + iabs(n4f(m,ict))
   23       nbas_tot=nbas_tot + iabs(n4f(m,ict))
   24     nbasis_ctype(ict)=nbas_typ
        if(nbas_tot.ne.nbasis) then
          write(6,'(''nbas_tot,nbasis='',9i4)') nbas_tot,nbasis
          stop 'nbas_tot not equal to nbasis'
        endif

        betaq=0
        do 25 ic=1,ncent
   25     betaq=betaq+znuc(iwctype(ic))
        betaq=betaq-nelec+1

       else

        do 49 ict=1,nctype
          write(6,'(/,''center type'',i3)') ict
          read(5,*) n1s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &    ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &    ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &    ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &    ,(n6h(m,ict),n6h(-m,ict),m=5,1,-1),n6h(0,ict)
          write(6,'(''number of basis fns. for each l:'',100i2)')
     &    n1s(ict),n2p(1,ict),n2p(-1,ict),n2p(0,ict)
     &    ,n3d(0,ict),(n3d(m,ict),n3d(-m,ict),m=2,1,-1)
     &    ,(n4f(m,ict),n4f(-m,ict),m=3,1,-1),n4f(0,ict)
     &    ,(n5g(m,ict),n5g(-m,ict),m=4,1,-1),n5g(0,ict)
     &    ,(n6h(m,ict),n6h(-m,ict),m=5,1,-1),n6h(0,ict)

          nbas=iabs(n1s(ict))
          n2s(ict)=0
          n3s(ict)=0
          n4s(ict)=0
          do 41 m=-1,1
            n3p(m,ict)=0
            n4p(m,ict)=0
            if(n2p(m,ict).ne.0 .and. ML_BAS.lt.1) stop 'ML_BAS in vmc.h too small'
   41       nbas=nbas+ iabs(n2p(m,ict))
          do 42 m=-2,2
            n4d(m,ict)=0
            if(n3d(m,ict).ne.0 .and. ML_BAS.lt.2) stop 'ML_BAS in vmc.h too small'
   42       nbas=nbas+ iabs(n3d(m,ict))
          do 43 m=-3,3
            if(n4f(m,ict).ne.0 .and. ML_BAS.lt.3) stop 'ML_BAS in vmc.h too small'
   43       nbas=nbas+ iabs(n4f(m,ict))
          do 44 m=-4,4
            if(n5g(m,ict).ne.0 .and. ML_BAS.lt.4) stop 'ML_BAS in vmc.h too small'
   44       nbas=nbas+ iabs(n5g(m,ict))
          do 45 m=-5,5
            if(n6h(m,ict).ne.0 .and. ML_BAS.lt.5) stop 'ML_BAS in vmc.h too small'
   45       nbas=nbas+ iabs(n6h(m,ict))
          if(nbas.gt.MBASIS_CTYPE) stop 'nbas > MBASIS_CTYPE'
          nbasis_ctype(ict)=nbas

          write(6,'(''Reading iwrwf'')')
          read(5,*) (iwrwf(ib,ict),ib=1,nbas)
          write(6,'(''Center'',i5,'' uses radial bas. fns:'',(100i3))')
     &    ict,(iwrwf(ib,ict),ib=1,nbas)

          do 46 ib=1,nbas
            if(iwrwf(ib,ict).gt.nrbas(ict)) then
              write(6,'(''ict,ib,iwrwf(ib,ict),nrbas(ict)'',9i3)')
     &        ict,ib,iwrwf(ib,ict),nrbas(ict)
              stop 'iwrwf(ib,ict) > nrbas(ict)'
            endif
   46     continue

          nrwf=1
          do 48 ib=2,nbas
            do 47 ib2=1,ib
   47         if(iwrwf(ib2,ict).eq.iwrwf(ib,ict)) goto 48
            nrwf=nrwf+1
   48     continue
   49     if(nrwf.gt.MRWF) stop 'nrwf > MRWF'

        nbas_tot=0
        do 55 ic=1,ncent
          ict=iwctype(ic)
          nbas_tot=nbas_tot+ iabs(n1s(ict))
          do 51 m=-1,1
   51       nbas_tot=nbas_tot+ iabs(n2p(m,ict))
          do 52 m=-2,2
   52       nbas_tot=nbas_tot+ iabs(n3d(m,ict))
          do 53 m=-3,3
   53       nbas_tot=nbas_tot+ iabs(n4f(m,ict))
          do 54 m=-4,4
   54       nbas_tot=nbas_tot+ iabs(n5g(m,ict))
          do 55 m=-5,5
   55       nbas_tot=nbas_tot+ iabs(n6h(m,ict))
        if(nbas_tot.ne.nbasis) stop 'nbas_tot not equal to nbasis'

      endif

      write(6,'(/,''center type'',(12i4))') (i,i=1,nctype)
      write(6,'(/,''1s'',t11,(12i4))') (n1s(i),i=1,nctype)
      if(numr.le.0) write(6,'(''2s'',t11,(12i4))') (n2s(i),i=1,nctype)
      write(6,'(''2px'',t11,(12i4))') (n2p(1,i),i=1,nctype)
      write(6,'(''2py'',t11,(12i4))') (n2p(-1,i),i=1,nctype)
      write(6,'(''2pz'',t11,(12i4))') (n2p(0,i),i=1,nctype)
      if(numr.le.0) then
        write(6,'(''3s'',t11,(12i4))') (n3s(i),i=1,nctype)
        write(6,'(''3px'',t11,(12i4))') (n3p(1,i),i=1,nctype)
        write(6,'(''3py'',t11,(12i4))') (n3p(-1,i),i=1,nctype)
        write(6,'(''3pz'',t11,(12i4))') (n3p(0,i),i=1,nctype)
      endif
      write(6,'(''3dzr'',t11,(12i4))') (n3d(0,i),i=1,nctype)
      write(6,'(''3dx2'',t11,(12i4))') (n3d(2,i),i=1,nctype)
      write(6,'(''3dxy'',t11,(12i4))') (n3d(-2,i),i=1,nctype)
      write(6,'(''3dxz'',t11,(12i4))') (n3d(1,i),i=1,nctype)
      write(6,'(''3dyz'',t11,(12i4))') (n3d(-1,i),i=1,nctype)
      if(numr.le.0) then
        write(6,'(''4s'',t11,(12i4))') (n4s(i),i=1,nctype)
        write(6,'(''4px'',t11,(12i4))') (n4p(1,i),i=1,nctype)
        write(6,'(''4py'',t11,(12i4))') (n4p(-1,i),i=1,nctype)
        write(6,'(''4pz'',t11,(12i4))') (n4p(0,i),i=1,nctype)
c       write(6,'(''4dzr'',t11,(12i4))') (n4d(0,i),i=1,nctype)
c       write(6,'(''4dx2'',t11,(12i4))') (n4d(2,i),i=1,nctype)
c       write(6,'(''4dxy'',t11,(12i4))') (n4d(-2,i),i=1,nctype)
c       write(6,'(''4dxz'',t11,(12i4))') (n4d(1,i),i=1,nctype)
c       write(6,'(''4dyz'',t11,(12i4))') (n4d(-1,i),i=1,nctype)
        write(6,'(''sa'',t11,(12i4))') (nsa(i),i=1,nctype)
        write(6,'(''pxa'',t11,(12i4))') (npa(1,i),i=1,nctype)
        write(6,'(''pya'',t11,(12i4))') (npa(-1,i),i=1,nctype)
        write(6,'(''pza'',t11,(12i4))') (npa(0,i),i=1,nctype)
        write(6,'(''dzra'',t11,(12i4))') (nda(0,i),i=1,nctype)
        write(6,'(''dx2a'',t11,(12i4))') (nda(2,i),i=1,nctype)
        write(6,'(''dxya'',t11,(12i4))') (nda(-2,i),i=1,nctype)
        write(6,'(''dxza'',t11,(12i4))') (nda(1,i),i=1,nctype)
        write(6,'(''dyza'',t11,(12i4))') (nda(-1,i),i=1,nctype)
       else
        do 60 m=-3,3
   60     write(6,'(''4f('',i2,'')'',t11,(12i4))') m,(n4f(m,i),i=1,nctype)
        do 70 m=-4,4
   70     write(6,'(''5g('',i2,'')'',t11,(12i4))') m,(n5g(m,i),i=1,nctype)
        do 80 m=-5,5
   80     write(6,'(''6h('',i2,'')'',t11,(12i4))') m,(n6h(m,i),i=1,nctype)
      endif
      write(6,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=1,nctype)
      write(6,*)

c the following sets up character strings that identify basis
c functions on centers for the printout of coefficients and
c screening constants. Character strings are of form
c 2px'(3) for the second 2px function on the third center

      ib=0
      do 250 ic=1,ncent
        i=iwctype(ic)

        write(lcent,'(''('',i2,'')'')') ic
        if(ic.lt.10)  write(lcent,'(''(0'',i1,'')'')') ic

        do 90 j=1,iabs(n1s(i))
c         if(iabs(n1s(i)).gt.nprime) stop 'number of 1s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(1,n1s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
   90     lbasis(ib)=l1s(min(j,nprime))//lcent

        do 100 j=1,iabs(n2s(i))
c         if(iabs(n2s(i)).gt.nprime) stop 'number of 2s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(2,n2s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  100     lbasis(ib)=l2s(min(j,nprime))//lcent

        do 120 j=1,iabs(n3s(i))
c         if(iabs(n3s(i)).gt.nprime) stop 'number of 3s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(3,n3s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
  120     lbasis(ib)=l3s(min(j,nprime))//lcent

        do 150 j=1,iabs(n4s(i))
c         if(iabs(n4s(i)).gt.nprime) stop 'number of 4s basis fns > nprime'
          ib=ib+1
          n_bas(ib)=sign(4,n4s(i))
          l_bas(ib)=0
          m_bas(ib)=0
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          ib=ib+1
  150     lbasis(ib)=l4s(min(j,nprime))//lcent

        do 110 m=1,0,-1
          do 105 j=1,iabs(n2p(m,i))
c           if(iabs(n2p(m,i)).gt.nprime) stop 'number of 2p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(2,n2p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  105       lbasis(ib)=l2p(min(j,nprime),m)//lcent
          do 125 j=1,iabs(n3p(m,i))
c           if(iabs(n3p(m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  125       lbasis(ib)=l3p(min(j,nprime),m)//lcent
          do 155 j=1,iabs(n4p(m,i))
c           if(iabs(n4p(m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4p(m,i))
            l_bas(ib)=1
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  155       lbasis(ib)=l4p(min(j,nprime),m)//lcent
          if(m.eq.0) goto 110
          do 108 j=1,iabs(n2p(-m,i))
c           if(iabs(n2p(-m,i)).gt.nprime) stop 'number of 2p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(2,n2p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  108       lbasis(ib)=l2p(min(j,nprime),-m)//lcent
          do 128 j=1,iabs(n3p(-m,i))
c           if(iabs(n3p(-m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  128       lbasis(ib)=l3p(min(j,nprime),-m)//lcent
          do 158 j=1,iabs(n4p(-m,i))
c           if(iabs(n4p(-m,i)).gt.nprime) stop 'number of 3p basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4p(-m,i))
            l_bas(ib)=1
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  158       lbasis(ib)=l4p(min(j,nprime),-m)//lcent
  110   continue

          m=0
          do 132 j=1,iabs(n3d(m,i))
c           if(iabs(n3d(m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  132       lbasis(ib)=l3d(min(j,nprime),m)//lcent
c         do 162 j=1,iabs(n4d(m,i))
c           if(iabs(n4d(m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
c           ib=ib+1
c           n_bas(ib)=sign(4,n4d(m,i))
c           l_bas(ib)=2
c           m_bas(ib)=m
c           icenter_basis(ib)=ic
c           ictype_basis(ib)=i
c 162       lbasis(ib)=l4d(min(j,nprime),m)//lcent
        do 175 m=2,1,-1
          do 135 j=1,iabs(n3d(m,i))
c           if(iabs(n3d(m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(m,i))
            l_bas(ib)=2
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  135       lbasis(ib)=l3d(min(j,nprime),m)//lcent
c         do 165 j=1,iabs(n4d(m,i))
c           if(iabs(n4d(m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
c           ib=ib+1
c           n_bas(ib)=sign(4,n4d(m,i))
c           l_bas(ib)=2
c           m_bas(ib)=m
c           icenter_basis(ib)=ic
c           ictype_basis(ib)=i
c 165       lbasis(ib)=l4d(min(j,nprime),m)//lcent
          do 140 j=1,iabs(n3d(-m,i))
c           if(iabs(n3d(-m,i)).gt.nprime) stop 'number of 3d basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(3,n3d(-m,i))
            l_bas(ib)=2
            m_bas(ib)=-m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
  140       lbasis(ib)=l3d(min(j,nprime),-m)//lcent
c         do 170 j=1,iabs(n4d(-m,i))
c           if(iabs(n4d(-m,i)).gt.nprime) stop 'number of 4d basis fns > nprime'
c           ib=ib+1
c           n_bas(ib)=sign(4,n4d(-m,i))
c           l_bas(ib)=2
c           m_bas(ib)=-m
c           icenter_basis(ib)=ic
c           ictype_basis(ib)=i
c 170       lbasis(ib)=l4d(min(j,nprime),-m)//lcent
  175     continue

c Ordering to be made consistent
        do 180 m=-3,3
          do 180 j=1,iabs(n4f(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 4f basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(4,n4f(m,i))
            l_bas(ib)=3
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
c           write(lmbas,'(''('',i2,'')'')') m
            write(lmbas,'(i2)') m
  180       lbasis(ib)=l4f(min(j,nprime))//lmbas//lcent

        do 190 m=-4,4
          do 190 j=1,iabs(n5g(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 5g basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(5,n5g(m,i))
            l_bas(ib)=4
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') m
  190       lbasis(ib)=l5g(min(j,nprime))//lmbas//lcent

        do 200 m=-5,5
          do 200 j=1,iabs(n6h(m,i))
c           if(iabs(n4f(m,i)).gt.nprime) stop 'number of 6h basis fns > nprime'
            ib=ib+1
            n_bas(ib)=sign(6,n6h(m,i))
            l_bas(ib)=5
            m_bas(ib)=m
            icenter_basis(ib)=ic
            ictype_basis(ib)=i
            write(lmbas,'(i2)') m
  200       lbasis(ib)=l6h(min(j,nprime))//lmbas//lcent

      if(iabs(nsa(i)).ge.1) then
        ib=ib+1
        n_bas(ib)=0
        l_bas(ib)=0
        m_bas(ib)=0
        icenter_basis(ib)=ic
        ictype_basis(ib)=i
        lbasis(ib)=lsa//lcent
      endif

      do 220 m=-1,1
        if(iabs(npa(m,i)).ge.1) then
          ib=ib+1
          n_bas(ib)=0
          l_bas(ib)=1
          m_bas(ib)=m
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          lbasis(ib)=lpa(m)//lcent
        endif
  220 continue

      do 230 m=-2,2
        if(iabs(nda(m,i)).ge.1) then
          ib=ib+1
          n_bas(ib)=0
          l_bas(ib)=2
          m_bas(ib)=m
          icenter_basis(ib)=ic
          ictype_basis(ib)=i
          lbasis(ib)=lda(m)//lcent
        endif
  230 continue

  250 continue

      write(6,'(/,(12a10))') (lbasis(j),j=1,nbasis)
      write(6,'(''orbital coefficients'')')
      do 260 i=1,norb
        read(5,*) (coef(j,i,1),j=1,nbasis)
  260   write(6,'(12f10.6)') (coef(j,i,1),j=1,nbasis)

c Warning: Foll. "if" commented out till I modify change_input to work without the zex
c     if(inumr.ge.1) then
        write(6,'(''screening constants'')')
        read(5,*) (zex(i,1),i=1,nbasis)
        write(6,'(12f10.6)') (zex(i,1),i=1,nbasis)
        do 262 i=1,nbasis
  262     if(zex(i,1).le.0.d0) stop 'exponent zex should be > 0'
c     endif

c If an analytical basis is used find out how many distinct radial orbitals there are
c and store n and zex for each of these.
c This is done so that only the unique analytic radial wavefns. are evaluated in basis_fns.f and basis_fnse.f
c zex2 is used in func when optimizing the exponents.
c ib  goes up to the total number of basis functions
c ib2 goes up to the number of basis functions on each center type
c irb goes up to the number of radial basis functions on each center type
      if(numr.le.0) then
        write(6,'(''n_bas(ib),l_bas(ib)'',50(i2,x))') (n_bas(ib),l_bas(ib),ib=1,nbasis)
        ib=0
        ib4=0
        do 290 ic=1,ncent
          ict=iwctype(ic)
          write(6,'(''n_bas(ib),l_bas(ib)'',50(i2,x))') (n_bas(i),l_bas(i),i=1,nbasis_ctype(ict))
          nrbas(ict)=0
          do 280 ib2=1,nbasis_ctype(ict)
            ib=ib+1
            if(ib.ne.ib2+ib4) stop 'ib .ne. ib2+ib4'
            do 270 ib3=1,ib2-1
              if(n_bas(ib).eq.n_bas(ib3+ib4)
     &        .and. zex(ib,1).eq.zex(ib3+ib4,1)) then
                irb=iwrwf2(ib3+ib4)
                goto 275
              endif
  270       continue
            nrbas(ict)=nrbas(ict)+1
            irb=nrbas(ict)
            if(irb.gt.MRWF) stop 'nbas > MRWF'
            n_bas2(irb,ict)=n_bas(ib)
            zex2(irb,ict,1)=zex(ib,1)
  275       iwrwf2(ib)=irb
  280       iwrwf(ib2,ict)=irb
  290       ib4=ib4+nbasis_ctype(ict)
        write(6,'(''nrbas(ict)='',40i4)') (nrbas(ict),ict=1,nctype)
        write(6,'(''iwrwf2(ib)='',40i3)') (iwrwf2(ib),ib=1,nbasis)
        write(6,'(''iwrwf(ib2,ict)='',40i3)') ((iwrwf(ib2,ict),ib2=1,nbasis_ctype(ict)),ict=1,nctype)
      endif

c     if(numr.le.0) then
c       ib=0
c       ib4=0
c       do 290 ict=1,nctype
c         irb=0
c         do 290 ib2=1,nbasis_ctype(ict)
c           ib=ib+1
c           do 270 ib3=1,ib2-1
c             write(6,'(''ict,ib2,ib3,ib4,ib='',9i5)') ict,ib2,ib3,ib4,ib
c 270         if(n_bas(ib)-l_bas(ib).eq.n_bas(ib3+ib4)-l_bas(ib3+ib4)
c    &        .and. zex(ib,1).eq.zex(ib3+ib4,1)) goto 280
c           irb=irb+1
c           n_bas2(irb,ict)=n_bas(ib)
c           zex2(irb,ict,1)=zex(ib,1)
c           iwrwf2(ib)=irb
c 280       iwrwf(ib2,ict)=irb
c           if(irb.gt.MRWF) stop 'nbas > MRWF'
c           nrbas(ict)=irb
c 290       ib4=ib4+nbasis_ctype(ict)
c     endif

      return
      end

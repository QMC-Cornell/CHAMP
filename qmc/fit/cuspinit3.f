      subroutine cuspinit3(iprin)
c Written by Claudia Filippi

      implicit real*8(a-h,o-z)

      parameter(zero=0.d0,one=1.d0)

      include '../vmc/vmc.h'
      include '../vmc/force.h'

      parameter(NEQSX=6*MORDJ)
      parameter(job=11)

      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /cuspmat/ cm(NEQSX,NEQSX),iwc(NEQSX),neqs,ishe

c     dimension q(2*MORDJ,NEQSX),qp(0:2*MORDJ,NEQSX)
      dimension qp(0:2*MORDJ,NEQSX)
      dimension ipvt(NEQSX),work(NEQSX),det(2),cs(NEQSX,NEQSX)

      ishe=nord
      neqs=2*nord

c Put the i^th dependent variable in iwc(i)
      jj=1
      ll=0
      do 3 jp=1,nord
        iu=iorder(jp-1)+(jp-1)*jp/2+1
        is=iu+jp
        iwc(jj)=iu
        jj=jj+1
        iwc(jj)=is
        jj=jj+1
   3    continue

      nvar=neqs
      call piksrt(nvar,iwc)

      do 25 ii=1,nvar
        do 25 jp=1,nord
c         q(jp+nord,ii)=zero
          qp(jp+nord-1,ii)=zero
c         q(jp,ii)=zero
  25      qp(jp-1,ii)=zero

      jj=1
      ll=0
      do 30 jp=1,nord
        jpsh=jp+nord
        do 30 ju=jp,0,-1
          jsx=jp-ju
          do 30 js=jsx,0,-1
            ll=ll+1
            jt=jsx-js
            jsg=1

            if(jj.le.nvar.and.ll.eq.iwc(jj)) then
              if(jt.eq.0) then
c               if(ju.eq.0) q(jpsh,jj)=one
                if(ju.eq.1) qp(jpsh-1,jj)=one
              endif
              if(mod(jt,2).ne.0) jsg=-1
c             q(jp,jj)=jsg
              qp(jp-1,jj)=jsg*(js-jt)
              jj=jj+1
            endif
  30        continue

      do 40 jp=1,nord
        do 40 jj=1,nvar
          cs(jp+ishe,jj)=qp(jp+nord-1,jj)
  40      cs(jp,jj)=qp(jp-1,jj)

      if(iprin.eq.1) then
        write(6,'(''dependent variables'',40i3)')
     &  (iwc(jj),jj=1,nvar)
        do 65 ii=1,nvar
  65      write(6,'(i3,'')'',30f8.2)') ii,(cs(ii,jj),jj=1,nvar)
      endif

      call dgefa(cs,NEQSX,nvar,ipvt,info)
      write(6,'(''dgefa status='',i5)') info
      call dgedi(cs,NEQSX,nvar,ipvt,det,work,job)
      determinant=det(1)*10.d0**det(2)
      if(iprin.eq.1) write(6,'(''determinat='',3f12.4)')
     &determinant,(det(k),k=1,2)

      do 70 i=1,nvar
        do 70 j=1,nvar
  70      cm(i,j)=cs(i,j)

      call checkdepend3

      return
      end

c-----------------------------------------------------------------------
      function iorder(n)

      implicit real*8(a-h,o-z)

      iorder=(n**3+5*n)/6+n**2+n
      if(n.eq.0) iorder=0

      return
      end

c-----------------------------------------------------------------------
      subroutine checkdepend3

      use atom_mod
      use optim_mod
      implicit real*8(a-h,o-z)

!JT      include '../vmc/vmc.h'
!JT      include 'fit.h'

      parameter(NEQSX=6*MORDJ)

      common /cuspmat/ cm(NEQSX,NEQSX),iwc(NEQSX),neqs,ishe

!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent

      common /vardep/ nvdepend(NEQSX,MCTYPE),iwdepend(NEQSX,MPARMJ,MCTYPE)
     &,cdep(NEQSX,MPARMJ,MCTYPE)

      nvar=neqs

      do 5 i=1,nvar
        do 5 it=1,nctype
          nvdepend(i,it)=0
          do 5 l=1,nparmc(it)
  5         iwdepend(i,l,it)=0


      do 40 i=1,nvar
        do 40 j=1,nvar
          if(dabs(cm(i,j)).gt.1.d-10) then
            jorder=mod(j,ishe)
            if(jorder.eq.0) jorder=ishe
            jj=2*(jorder-1)+1
            ll=iorder(jorder-1)
            do 30 ju=jorder,0,-1
              jsx=jorder-ju
              do 30 js=jsx,0,-1
                ll=ll+1
                jt=jsx-js
                jsg=1
                if(jj.le.nvar.and.ll.eq.iwc(jj)) then
                  jj=jj+1
                 else
                  if(j.gt.ishe.and.jt.eq.0.and.ju.eq.1) then
                    do 10 it=1,nctype
                      do 10 l=1,nparmc(it)
                        if(ll.eq.iwjasc(l,it)) then
                          nvdepend(i,it)=nvdepend(i,it)+1
                          iwdepend(i,nvdepend(i,it),it)=l
                          cdep(i,nvdepend(i,it),it)=-cm(i,j)
                        endif
  10                continue
c                   write(6,'(''ivar'',i3,'' coef,idep'',f10.5,2i3)')
c    &              i,-cm(i,j),ll
                  endif
                  if(j.le.ishe.and.js-jt.ne.0) then
                    if(mod(jt,2).ne.0) jsg=-1
                    do 20 it=1,nctype
                      do 20 l=1,nparmc(it)
                        if(ll.eq.iwjasc(l,it)) then
                          nvdepend(i,it)=nvdepend(i,it)+1
                          iwdepend(i,nvdepend(i,it),it)=l
                          cdep(i,nvdepend(i,it),it)=-jsg*(js-jt)*cm(i,j)
                        endif
  20                continue
c                   write(6,'(''ivar'',i3,'' coef,idep'',f10.5,2i3)')
c    &              i,-jsg*(js-jt)*cm(i,j),ll
                  endif
                endif
  30        continue
          endif
  40  continue

c     do 60 it=1,nctype
c       do 60 i=1,nvar
c         write(6,'(''ivar'',i3,'' ndep'',i3)') i,nvdepend(i,it)
c         write(6,'(10i3)') (iwdepend(i,j,it),j=1,nvdepend(i,it))
c 60      write(6,'(10f10.3)') (cdep(i,j,it),j=1,nvdepend(i,it))

      return
      end

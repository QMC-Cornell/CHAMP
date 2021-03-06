      subroutine cuspinit3(iprin)
! Written by Claudia Filippi
      use all_tools_mod
      use constants_mod
      use jaspar3_mod
      use cuspmat_mod
      implicit real*8(a-h,o-z)

      parameter(job=11)

      dimension qp(0:2*nord,neqs)
      dimension ipvt(neqs),work(neqs),det(2),cs(neqs,neqs)

      ishe=nord
      neqs=2*nord

      call alloc ('cm', cm, neqs, neqs)
      call alloc ('iwc', iwc, neqs)

! Put the i^th dependent variable in iwc(i)
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
!         q(jp+nord,ii)=zero
          qp(jp+nord-1,ii)=zero
!         q(jp,ii)=zero
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
!               if(ju.eq.0) q(jpsh,jj)=one
                if(ju.eq.1) qp(jpsh-1,jj)=one
              endif
              if(mod(jt,2).ne.0) jsg=-1
!             q(jp,jj)=jsg
              qp(jp-1,jj)=jsg*(js-jt)
              jj=jj+1
            endif
  30        continue

      do 40 jp=1,nord
        do 40 jj=1,nvar
          cs(jp+ishe,jj)=qp(jp+nord-1,jj)
  40      cs(jp,jj)=qp(jp-1,jj)

      if(iprin.eq.1) then
        write(6,'(''dependent variables'',40i3)') (iwc(jj),jj=1,nvar)
        do 65 ii=1,nvar
  65      write(6,'(i3,'')'',30f8.2)') ii,(cs(ii,jj),jj=1,nvar)
      endif

      call dgefa(cs,nvar,nvar,ipvt,info)
      write(6,'(''dgefa status='',i5)') info
      call dgedi(cs,nvar,nvar,ipvt,det,work,job)
      determinant=det(1)*10.d0**det(2)
      if(iprin.eq.1) write(6,'(''determinat='',3f12.4)') determinant,(det(k),k=1,2)

      do 70 i=1,nvar
        do 70 j=1,nvar
  70      cm(i,j)=cs(i,j)

      call checkdepend3

      return
      end

!-----------------------------------------------------------------------
      function iorder(n)

      implicit real*8(a-h,o-z)

      iorder=(n**3+5*n)/6+n**2+n
      if(n.eq.0) iorder=0

      return
      end

!-----------------------------------------------------------------------
      subroutine checkdepend3

      use all_tools_mod
      use atom_mod
      use optim_mod
      use vardep_mod
      use cuspmat_mod
      implicit real*8(a-h,o-z)

      call alloc ('nvdepend', nvdepend, neqs, nctype)
      call alloc ('iwdepend', iwdepend, neqs, nparmj, nctype)
      call alloc ('cdep', cdep, neqs, nparmj, nctype)

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
!                   write(6,'(''ivar'',i3,'' coef,idep'',f10.5,2i3)') i,-cm(i,j),ll
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
!                   write(6,'(''ivar'',i3,'' coef,idep'',f10.5,2i3)') i,-jsg*(js-jt)*cm(i,j),ll
                  endif
                endif
  30        continue
          endif
  40  continue

!     do 60 it=1,nctype
!       do 60 i=1,nvar
!         write(6,'(''ivar'',i3,'' ndep'',i3)') i,nvdepend(i,it)
!         write(6,'(10i3)') (iwdepend(i,j,it),j=1,nvdepend(i,it))
! 60      write(6,'(10f10.3)') (cdep(i,j,it),j=1,nvdepend(i,it))

      return
      end

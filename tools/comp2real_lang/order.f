      subroutine def_ses

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: assigns a unique number to each set of (s,n,l,m) which
*              helps to define ordered single electron orbitals
* arguments: none
* return values: ord_ses(is,in,il,im)=i^{th} overall quantum number
* note: the array ord_ses may contain unphysical 3D atomic orbitals,
*       but this has no importance in this context. (later if we
*       need a specific order for qmc input we may need to change it)
*---------------------------------------------------------------------

      include         'maxdim.h'

c locals:

      integer         is,in,im,il,icount

c commons:

      integer         ord_ses(-MAXS:MAXS,0:MAXN,0:MAXL,-MAXM:MAXM)
      common/order/   ord_ses

      write(*,*) '* Creating ordered single electron states'
      icount=0
      do 10 is=MAXS,-MAXS,-2
        do 10 im=-MAXM,MAXM
          do 10 il=0,MAXL
            do 10 in=0,MAXN
              icount=icount+1
              ord_ses(is,in,il,im)=icount
c              write(*,*) is,in,il,im,icount
 10   continue

      write(*,*) '  OK.'
      write(*,*)
      return
      end

**********************************************************************
      subroutine order_det

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: bubblesort the orbitals in each determinant. Also
*              look for Pauli principle to eliminate null determinants
*              (we don't "kill" determinants but just set their coef. to
*               zero)
* arguments: none
* return values: det and cdet are modified
*---------------------------------------------------------------------

      include         'maxdim.h'
      include         'qnumbers.h'

c locals:
      integer         idet
      integer         is,in,il,im,istate,ii
      integer         js,jn,jl,jm,jstate,jj
      logical         null

c common variables:
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)
      integer         ord_ses(-MAXS:MAXS,0:MAXN,0:MAXL,-MAXM:MAXM)

      common/psi_n/   cdet,ndim,nelec,ndet,det
      common/order/   ord_ses

      write(*,*) '* Ordering each determinant'

      do 10 idet=1,ndet
        null=.false.
        ii=0
        do 10 while(.not.null .and. ii.lt.nelec-1)
          ii=ii+1
          is=det(idet,ii,qns)
          in=det(idet,ii,qnn)
          im=det(idet,ii,qnm)
          il=det(idet,ii,qnl)
          istate=ord_ses(is,in,il,im)
          jj=ii
          do 10 while(.not.null .and. jj.lt.nelec)
            jj=jj+1
            js=det(idet,jj,qns)
            jn=det(idet,jj,qnn)
            jm=det(idet,jj,qnm)
            jl=det(idet,jj,qnl)
            jstate=ord_ses(js,jn,jl,jm)
c            write(*,*) idet,ii,istate,jj,jstate
c check Pauli principle:
            if (jstate.eq.istate) then
c              write(*,*) 'null det:',idet,istate,jstate
              null=.true.
              cdet(idet)=dcmplx(0.d0,0.d0)
c switch electrons if not ordered:
            elseif (jstate.lt.istate) then

              det(idet,ii,qns)=js
              det(idet,ii,qnn)=jn
              det(idet,ii,qnl)=jl
              det(idet,ii,qnm)=jm

              det(idet,jj,qns)=is
              det(idet,jj,qnn)=in
              det(idet,jj,qnl)=il
              det(idet,jj,qnm)=im

              is=det(idet,ii,qns)
              in=det(idet,ii,qnn)
              im=det(idet,ii,qnm)
              il=det(idet,ii,qnl)
              istate=ord_ses(is,in,il,im)
c update amplitude due to switch:
              cdet(idet)=-cdet(idet)
            endif
 10   enddo

      write(*,*) '  OK.'
      write(*,*)
      return
      end

********************************************************************

      subroutine reduce_det

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: scans all the determinants and compares non-zero ones.
*              similar dets. are combined togather. Note that we do
*              not "kill" any determinants. Just set their amplitude
*              to zero.
* arguments: none
* return values: det and cdet are modified
*---------------------------------------------------------------------

      include         'maxdim.h'

c locals:
      integer         idet,jdet
      complex*16      czero
      logical         equal
c common variables:
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)
      common/psi_n/   cdet,ndim,nelec,ndet,det

c functions
      logical         compare

      write(*,*) '* Reducing wave-function'

      czero=dcmplx(0.d0,0.d0)

c compare only non-null determinants:
      do 10 idet=1,ndet-1
        if (cdet(idet).ne.czero) then
          do 5 jdet=idet+1,ndet
            if (cdet(jdet).ne.czero) then
              equal=compare(idet,jdet)
c if idet=jdet, add jdet to idet, and set jdet to zero
              if (equal) then
                cdet(idet)=cdet(idet)+cdet(jdet)
                cdet(jdet)=czero
              endif
            endif
 5        enddo
        endif
 10   enddo

      write(*,*) '  OK.'
      write(*,*)

      return
      end

**********************************************************************
      function compare(idet,jdet)

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: compares two ORDERED determinants
* arguments: idet,jdet : index of determinants to compare
* return values:compare=true if determinants are equal
* note: for many other calculations such as exact diag., this
*       algorthim is very intensive. If you know any better
*       determinant comparison algorithm, please let me know!!
*---------------------------------------------------------------------

      include         'maxdim.h'
      include         'qnumbers.h'

      logical         compare
c arguments:
      integer         idet,jdet

c locals
      integer         nn
      integer         is,in,im,il,istate
      integer         js,jn,jm,jl,jstate

c common variables:
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)
      integer         ord_ses(-MAXS:MAXS,0:MAXN,0:MAXL,-MAXM:MAXM)

      common/psi_n/   cdet,ndim,nelec,ndet,det
      common/order/   ord_ses

      compare=.true.
      nn=1
c keep comparing orbitals until 1 different is found:
      do 10 while (compare .and. nn.le.nelec)

        is=det(idet,nn,qns)
        in=det(idet,nn,qnn)
        im=det(idet,nn,qnm)
        il=det(idet,nn,qnl)
        istate=ord_ses(is,in,il,im)

        js=det(jdet,nn,qns)
        jn=det(jdet,nn,qnn)
        jm=det(jdet,nn,qnm)
        jl=det(jdet,nn,qnl)
        jstate=ord_ses(js,jn,jl,jm)

        if (istate .ne. jstate) compare=.false.
c is that also correct?: compare=istate.eq.jstate

        nn=nn+1
 10   enddo

      return
      end

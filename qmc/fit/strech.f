      subroutine strech_fit
! Written by Cyrus Umrigar and Claudia Filippi
! Uses the coordinate transform described in:
! 1) Two Aspects of Quantum Monte Carlo: Determination of Accurate Wavefunctions and
!    Determination of Potential Energy Surfaces of Molecules, C.J. Umrigar,
!    Int. J. Quant. Chem. Symp., 23, 217 (1989).
! 2) Correlated sampling in quantum Monte Carlo: A route to forces,
!    Claudia Filippi and C. J. Umrigar, Phys. Rev. B., 61, R16291, (2000).

! Stretch space so that electrons close to a nucleus move almost
! rigidly with that nucleus
      use atom_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use confg_mod
      implicit real*8(a-h,o-z)

      dimension delc(3,ncent),wt(ncent),dvol(3,3),dwt(3,ncent),dwtsm(3)
     &,oldcent(3,ncent)

!JT      save alfstr,delc
      save alfstr

      read(5,*) alfstr
      write(6,'(''istrch,alfstr'',i5,f7.2)') istrch,alfstr
      if(alfstr.eq.0.d0) return

! Read in old atom positions and find shift in positions
      do 5 icent=1,ncent
        read(5,*) (oldcent(k,icent),k=1,ndim)
        do 5 k=1,ndim
    5     delc(k,icent)=cent(k,icent)-oldcent(k,icent)

      do 50 idata=1,ndata
        do 50 i=1,nelec
          wtsm=0

! Initialize volume change matrix
          do 7 k=1,ndim
            dwtsm(k)=0
            do 7 j=1,ndim
              dvol(j,k)=0
              if(j.eq.k) dvol(j,k)=1
    7     continue

          do 20 icent=1,ncent
            dist2=0
            do 10 k=1,ndim
   10         dist2=dist2+(x(k,i,idata)-oldcent(k,icent))**2
            dist=dsqrt(dist2)
            if(istrch.eq.1) wt(icent)=dexp(-alfstr*dist)
            if(istrch.eq.2) wt(icent)=1/dist**alfstr
            if(istrch.eq.3) wt(icent)=dexp(alfstr/dist)
            wtsm=wtsm+wt(icent)
            do 20 k=1,ndim
              if(istrch.eq.1) dwt(k,icent)=-alfstr*dexp(-alfstr*dist)*
     &        (x(k,i,idata)-oldcent(k,icent))/dist
              if(istrch.eq.2) dwt(k,icent)=-alfstr*
     &        (x(k,i,idata)-oldcent(k,icent))/dist**(alfstr+2)
              if(istrch.eq.3) dwt(k,icent)=-alfstr*dexp(alfstr/dist)*
     &        (x(k,i,idata)-oldcent(k,icent))/(dist2*dist)
              dwtsm(k)=dwtsm(k)+dwt(k,icent)
   20     continue
          wtsmi=1/wtsm

          do 40 icent=1,ncent
            do 40 k=1,ndim

              dwt(k,icent)=(wtsm*dwt(k,icent)-wt(icent)*dwtsm(k))/
     &        wtsm**2
              dvol(1,k)=dvol(1,k)+dwt(k,icent)*delc(1,icent)
              dvol(2,k)=dvol(2,k)+dwt(k,icent)*delc(2,icent)
   40         dvol(3,k)=dvol(3,k)+dwt(k,icent)*delc(3,icent)
          call matinv(dvol,3,det)
          dvpdv(idata)=dvpdv(idata)*det
          do 50 icent=1,ncent
            wt(icent)=wt(icent)*wtsmi
            do 50 k=1,ndim
   50         x(k,i,idata)=x(k,i,idata)+wt(icent)*delc(k,icent)

      return
      end

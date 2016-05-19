      subroutine pairden2d(p,q,xold,xnew)

! Written by A.D.Guclu jun2005.
! Calculates the full pair-densities reducing the dimensionality
! by 1 due to circular symmetry (2+1d instead of 2+2d).
! For the moment does not distinguish between all and 1 electron calculation.
! The reason is that even when only 1 electron is moved, several of
! the relative distances changes, making it diffcult to
! keep track of all the rotated-relative distances.
! (not impossible, can be optimized)

      use dets_mod
      use const_mod
      use dim_mod
      use pairden_mod
      implicit real*8(a-h,o-z)

      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys
      dimension xold(3,nelec),xnew(3,nelec)

      do 30 ier=1,nelec      ! reference electron

        rold=0.d0
        rnew=0.d0
        do 10 idim=1,ndim
          rold=rold+xold(idim,ier)**2
          rnew=rnew+xnew(idim,ier)**2
   10   enddo
        rold=dsqrt(rold)
        rnew=dsqrt(rnew)
        thetao=datan2(xold(2,ier),xold(1,ier))
        thetan=datan2(xnew(2,ier),xnew(1,ier))
        if(icoosys.eq.1) then
          iro=nint(delxi(2)*rold)
          irn=nint(delxi(2)*rnew)
        else
          iro=nint(delradi*(rold - xfix(1)))
          irn=nint(delradi*(rnew - xfix(1)))
        endif
        if((iro.lt.0 .or. iro.gt.NAX) .and. (irn.lt.0 .or. irn.gt.NAX)) cycle
! electron relative to the reference electron
        do 20 ie2=1,nelec
          if(ie2.ne.ier) then
! rotate old and new coordinates
            call rotate(thetao,xold(1,ie2),xold(2,ie2),x1roto,x2roto)
            call rotate(thetan,xnew(1,ie2),xnew(2,ie2),x1rotn,x2rotn)
! put on the grid:
            if(icoosys.eq.1) then
              ix1roto=nint(delxi(1)*x1roto)
              ix2roto=nint(delxi(2)*x2roto)
              ix1rotn=nint(delxi(1)*x1rotn)
              ix2rotn=nint(delxi(2)*x2rotn)
            else
! same trick adapted to circular coordinates
              ix1roto=nint(delradi*(sqrt(x1roto**2 + x2roto**2)-rmean))
              ix1rotn=nint(delradi*(sqrt(x1rotn**2 + x2rotn**2)-rmean))
              ix2roto=nint(delti*(datan2(x2roto,x1roto)))
              ix2rotn=nint(delti*(datan2(x2rotn,x1rotn)))
            endif


! check if we are within grid limits, check spins, and collect data
!  -old config
            if(iro.le.NAX .and. iro.ge.0 .and. abs(ix1roto).le.NAX .and. abs(ix2roto).le.NAX) then
              if(ier.le.nup) then
                xx0probut(iro,ix1roto,ix2roto)=xx0probut(iro,ix1roto,ix2roto)+q
                if(ie2.le.nup) then
                  xx0probuu(iro,ix1roto,ix2roto)=xx0probuu(iro,ix1roto,ix2roto)+q
                else
                  xx0probud(iro,ix1roto,ix2roto)=xx0probud(iro,ix1roto,ix2roto)+q
                endif
              else
                xx0probdt(iro,ix1roto,ix2roto)=xx0probdt(iro,ix1roto,ix2roto)+q
                if(ie2.le.nup) then
                  xx0probdu(iro,ix1roto,ix2roto)=xx0probdu(iro,ix1roto,ix2roto)+q
                else
                  xx0probdd(iro,ix1roto,ix2roto)=xx0probdd(iro,ix1roto,ix2roto)+q
                endif
              endif
            endif
! -new config
            if(irn.le.NAX .and. irn.ge.0 .and. abs(ix1rotn).le.NAX .and. abs(ix2rotn).le.NAX) then
              if(ier.le.nup) then
                xx0probut(irn,ix1rotn,ix2rotn)=xx0probut(irn,ix1rotn,ix2rotn)+p
                if(ie2.le.nup) then
                  xx0probuu(irn,ix1rotn,ix2rotn)=xx0probuu(irn,ix1rotn,ix2rotn)+p
                else
                  xx0probud(irn,ix1rotn,ix2rotn)=xx0probud(irn,ix1rotn,ix2rotn)+p
                endif
              else
                xx0probdt(irn,ix1rotn,ix2rotn)=xx0probdt(irn,ix1rotn,ix2rotn)+p
                if(ie2.le.nup) then
                  xx0probdu(irn,ix1rotn,ix2rotn)=xx0probdu(irn,ix1rotn,ix2rotn)+p
                else
                  xx0probdd(irn,ix1rotn,ix2rotn)=xx0probdd(irn,ix1rotn,ix2rotn)+p
                endif
              endif
            endif
          endif
   20   enddo

   30 enddo

      return
      end

!------------------------------------------------------------------------------------

      subroutine rotate(theta,x1,x2,xrot1,xrot2)

! rotates (x1,x2) by theta. Result is (xrot1,xrot2)

      implicit real*8(a-h,o-z)

      thetarot=datan2(x2,x1)-theta
      r=dsqrt(x1*x1+x2*x2)
      xrot1=r*dcos(thetarot)
      xrot2=r*dsin(thetarot)

      return
      end

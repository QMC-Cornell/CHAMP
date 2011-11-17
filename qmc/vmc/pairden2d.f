      subroutine pairden2d(p,q,xold,xnew)

c Written by A.D.Guclu jun2005.
c Calculates the full pair-densities reducing the dimensionality
c by 1 due to circular symmetry (2+1d instead of 2+2d).
c For the moment does not distinguish between all and 1 electron calculation.
c The reason is that even when only 1 electron is moved, several of
c the relative distances changes, making it diffcult to
c keep track of all the rotated-relative distances.
c (not impossible, can be optimized)

      use dets_mod
      use const_mod
      use dim_mod
      use pairden_mod
      implicit real*8(a-h,o-z)

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
        iro=nint(delxi*rold)
        irn=nint(delxi*rnew)

c electron relative to the reference electron
        do 20 ie2=1,nelec
          if(ie2.ne.ier) then
c rotate old and new coordinates
            call rotate(thetao,xold(1,ie2),xold(2,ie2),x1roto,x2roto)
            call rotate(thetan,xnew(1,ie2),xnew(2,ie2),x1rotn,x2rotn)
c put on the grid:
c           if(icoosys.eq.1) then 
              ix1roto=nint(delxi*x1roto)
              ix2roto=nint(delxi*x2roto)
              ix1rotn=nint(delxi*x1rotn)
              ix2rotn=nint(delxi*x2rotn)
c           else
c same trick adapted to circular coordinates
c              ixo(1)=nint(delradi*(rold-rmean))
c              ixn(1)=nint(delradi*(rnew-rmean))
c              ixo(2)=nint(delti*(datan2(xold(2,i),xold(1,i))))
c              ixn(2)=nint(delti*(datan2(xnew(2,i),xnew(1,i))))
c           endif


c check if we are within grid limits, check spins, and collect data
c  -old config
            if(iro.le.NAX .and. abs(ix1roto).le.NAX .and. abs(ix2roto).le.NAX) then
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
c -new config
            if(irn.le.NAX .and. abs(ix1rotn).le.NAX .and. abs(ix2rotn).le.NAX) then
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

c------------------------------------------------------------------------------------

      subroutine rotate(theta,x1,x2,xrot1,xrot2)

c rotates (x1,x2) by theta. Result is (xrot1,xrot2)

      implicit real*8(a-h,o-z)

      thetarot=datan2(x2,x1)-theta
      r=dsqrt(x1*x1+x2*x2)
      xrot1=r*dcos(thetarot)
      xrot2=r*dsin(thetarot)

      return
      end

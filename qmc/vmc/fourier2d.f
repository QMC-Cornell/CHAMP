      subroutine fourierrk(p,q,xold,xnew)

c Written by A.D.Guclu aug2006.
c Calculates the power spectrum estimator in polar coo., defined as:
c     FT=abs( sum_i Fourier[Delta^2(r_i-r)])^2
c taking the absolute value before the data accumulation is essential to get the
c "internal structure" of the density.

c Notes:
c Delta^2(ri-r)=Delta(ri_r) Delta(thetai-theta) / r
c We consider only angular FT here.
c Fourier[Delta]=cos + i sin

c In practice, there is a uninteresting constant term that we substract:
c FT(renormalized)= r^2 FT - N(r)   where N(r) is number of electrons at r.

c During the print out we divide everything by r^2 to get the final result.
c FT(r,k=0)= electronic density


      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /pairden/ xx0probut(0:NAX,-NAX:NAX,-NAX:NAX),xx0probuu(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probud(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdt(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probdu(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdd(0:NAX,-NAX:NAX,-NAX:NAX),
     &den2d_t(-NAX:NAX,-NAX:NAX),den2d_d(-NAX:NAX,-NAX:NAX),den2d_u(-NAX:NAX,-NAX:NAX),
     &delxi,xmax,xfix(3),ifixe
      common /fourier/ fourierrk_u(0:NAX,0:NAK1),fourierrk_d(0:NAX,0:NAK1)
     &,fourierrk_t(0:NAX,0:NAK1),fourierkk_u(-NAK2:NAK2,-NAK2:NAK2),fourierkk_d(-NAK2:NAK2,-NAK2:NAK2)
     &,fourierkk_t(-NAK2:NAK2,-NAK2:NAK2),delk1,delk2,fmax1,fmax2,ifourier

      dimension xold(3,MELEC),xnew(3,MELEC)
      dimension fcos_uo(NAX),fcos_do(NAX)
      dimension fsin_uo(NAX),fsin_do(NAX)
      dimension fcos_un(NAX),fcos_dn(NAX)
      dimension fsin_un(NAX),fsin_dn(NAX)
      dimension nr_uo(NAX),nr_do(NAX),nr_un(NAX),nr_dn(NAX)

      do ik=0,NAK1
        fk=delk1*ik
c reset temporary arrays:
        do ix=1,NAX
          fcos_uo(ix)=0.d0
          fcos_do(ix)=0.d0
          fsin_uo(ix)=0.d0
          fsin_do(ix)=0.d0
          fcos_un(ix)=0.d0
          fcos_dn(ix)=0.d0
          fsin_un(ix)=0.d0
          fsin_dn(ix)=0.d0
          nr_uo(ix)=0
          nr_do(ix)=0
          nr_un(ix)=0
          nr_dn(ix)=0
        enddo

        do ie=1,nelec

          rold=0.d0
          rnew=0.d0
          do  idim=1,ndim
            rold=rold+xold(idim,ie)**2
            rnew=rnew+xnew(idim,ie)**2
          enddo
          rold=dsqrt(rold)
          rnew=dsqrt(rnew)
          thetao=datan2(xold(2,ie),xold(1,ie))
          thetan=datan2(xnew(2,ie),xnew(1,ie))
          iro=min(int(delxi*rold)+1,NAX)
          irn=min(int(delxi*rnew)+1,NAX)

          if(ie.le.nup) then
            fcos_uo(iro)=fcos_uo(iro)+dcos(fk*thetao)
            fsin_uo(iro)=fsin_uo(iro)+dsin(fk*thetao)
            fcos_un(irn)=fcos_un(irn)+dcos(fk*thetan)
            fsin_un(irn)=fsin_un(irn)+dsin(fk*thetan)
            nr_uo(iro)=nr_uo(iro)+1
            nr_un(irn)=nr_un(irn)+1
          else
            fcos_do(iro)=fcos_do(iro)+dcos(fk*thetao)
            fsin_do(iro)=fsin_do(iro)+dsin(fk*thetao)
            fcos_dn(irn)=fcos_dn(irn)+dcos(fk*thetan)
            fsin_dn(irn)=fsin_dn(irn)+dsin(fk*thetan)
            nr_do(iro)=nr_do(iro)+1
            nr_dn(irn)=nr_dn(irn)+1
          endif

        enddo

        do ir=1,NAX
          fourierrk_u(ir,ik)=fourierrk_u(ir,ik)+(fcos_uo(ir)*fcos_uo(ir)+fsin_uo(ir)*fsin_uo(ir))*q
          fourierrk_d(ir,ik)=fourierrk_d(ir,ik)+(fcos_do(ir)*fcos_do(ir)+fsin_do(ir)*fsin_do(ir))*q
          fourierrk_t(ir,ik)=fourierrk_t(ir,ik)+((fcos_uo(ir)+fcos_do(ir))**2
     &                                               +(fsin_uo(ir)+fsin_do(ir))**2)*q
          fourierrk_u(ir,ik)=fourierrk_u(ir,ik)+(fcos_un(ir)*fcos_un(ir)+fsin_un(ir)*fsin_un(ir))*p
          fourierrk_d(ir,ik)=fourierrk_d(ir,ik)+(fcos_dn(ir)*fcos_dn(ir)+fsin_dn(ir)*fsin_dn(ir))*p
          fourierrk_t(ir,ik)=fourierrk_t(ir,ik)+((fcos_un(ir)+fcos_dn(ir))**2
     &                                               +(fsin_un(ir)+fsin_dn(ir))**2)*p
c substract the infinite k limit:
          fourierrk_u(ir,ik)=fourierrk_u(ir,ik)-nr_uo(ir)*q
          fourierrk_d(ir,ik)=fourierrk_d(ir,ik)-nr_do(ir)*q
          fourierrk_t(ir,ik)=fourierrk_t(ir,ik)-(nr_uo(ir)+nr_do(ir))*q

          fourierrk_u(ir,ik)=fourierrk_u(ir,ik)-nr_un(ir)*p
          fourierrk_d(ir,ik)=fourierrk_d(ir,ik)-nr_dn(ir)*p
          fourierrk_t(ir,ik)=fourierrk_t(ir,ik)-(nr_un(ir)+nr_dn(ir))*p

        enddo
c        write(*,*) 'fourierrk_t=',fourierrk_t(ir,ik)
      enddo

      return
      end

c-------------------------------------------------------------------------------------
      subroutine fourierrk_old(p,q,xold,xnew)

c Written by A.D.Guclu aug2006.
c Calculates the "internal fourier transform" estimator in polar coo., defined as:
c     FT=abs( sum_i Fourier[Delta^2(r_i-r)])
c taking the absolute value before the data accumlation is essential to get the
c "internal structure" of the density. We are indeed calculation the following
c quantity:
c     <FT>=<Psi| abs( Fourier[rho(r)] ) |Psi>
c Notes:
c Delta^2(ri-r)=Delta(ri_r) Delta(thetai-theta) / r
c We consider only angular FT here.
c Fourier[Delta]=cos + i sin

c FT(r,k=0)= electronic density


      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /pairden/ xx0probut(0:NAX,-NAX:NAX,-NAX:NAX),xx0probuu(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probud(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdt(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probdu(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdd(0:NAX,-NAX:NAX,-NAX:NAX),
     &den2d_t(-NAX:NAX,-NAX:NAX),den2d_d(-NAX:NAX,-NAX:NAX),den2d_u(-NAX:NAX,-NAX:NAX),
     &delxi,xmax,xfix(3),ifixe
      common /fourier/ fourierrk_u(0:NAX,0:NAK1),fourierrk_d(0:NAX,0:NAK1)
     &,fourierrk_t(0:NAX,0:NAK1),fourierkk_u(-NAK2:NAK2,-NAK2:NAK2),fourierkk_d(-NAK2:NAK2,-NAK2:NAK2)
     &,fourierkk_t(-NAK2:NAK2,-NAK2:NAK2),delk1,delk2,fmax1,fmax2,ifourier

      dimension xold(3,MELEC),xnew(3,MELEC)
      dimension fcos_uo(NAX),fcos_do(NAX)
      dimension fsin_uo(NAX),fsin_do(NAX)
      dimension fcos_un(NAX),fcos_dn(NAX)
      dimension fsin_un(NAX),fsin_dn(NAX)

      big=1000.d0

      do ik=0,NAK1
        if(ik.eq.NAK1) then
          fk=big    !last point represents infinity. inconsistent with the new den2dwrt.f
        else
          fk=delk1*ik
        endif
c reset temporary arrays:
        do ix=1,NAX
          fcos_uo(ix)=0.d0
          fcos_do(ix)=0.d0
          fsin_uo(ix)=0.d0
          fsin_do(ix)=0.d0
          fcos_un(ix)=0.d0
          fcos_dn(ix)=0.d0
          fsin_un(ix)=0.d0
          fsin_dn(ix)=0.d0
        enddo

        do ie=1,nelec

          rold=0.d0
          rnew=0.d0
          do  idim=1,ndim
            rold=rold+xold(idim,ie)**2
            rnew=rnew+xnew(idim,ie)**2
          enddo
          rold=dsqrt(rold)
          rnew=dsqrt(rnew)
          thetao=datan2(xold(2,ie),xold(1,ie))
          thetan=datan2(xnew(2,ie),xnew(1,ie))
          iro=min(int(delxi*rold)+1,NAX)
          irn=min(int(delxi*rnew)+1,NAX)

          if(ie.le.nup) then
            fcos_uo(iro)=fcos_uo(iro)+dcos(fk*thetao)
            fsin_uo(iro)=fsin_uo(iro)+dsin(fk*thetao)
            fcos_un(irn)=fcos_un(irn)+dcos(fk*thetan)
            fsin_un(irn)=fsin_un(irn)+dsin(fk*thetan)
          else
            fcos_do(iro)=fcos_do(iro)+dcos(fk*thetao)
            fsin_do(iro)=fsin_do(iro)+dsin(fk*thetao)
            fcos_dn(irn)=fcos_dn(irn)+dcos(fk*thetan)
            fsin_dn(irn)=fsin_dn(irn)+dsin(fk*thetan)
          endif

        enddo

        do ir=1,NAX
          fourierrk_u(ir,ik)=fourierrk_u(ir,ik)+dsqrt(fcos_uo(ir)*fcos_uo(ir)+fsin_uo(ir)*fsin_uo(ir))*q
          fourierrk_d(ir,ik)=fourierrk_d(ir,ik)+dsqrt(fcos_do(ir)*fcos_do(ir)+fsin_do(ir)*fsin_do(ir))*q
          fourierrk_t(ir,ik)=fourierrk_t(ir,ik)+dsqrt((fcos_uo(ir)+fcos_do(ir))**2
     &                                               +(fsin_uo(ir)+fsin_do(ir))**2)*q
          fourierrk_u(ir,ik)=fourierrk_u(ir,ik)+dsqrt(fcos_un(ir)*fcos_un(ir)+fsin_un(ir)*fsin_un(ir))*p
          fourierrk_d(ir,ik)=fourierrk_d(ir,ik)+dsqrt(fcos_dn(ir)*fcos_dn(ir)+fsin_dn(ir)*fsin_dn(ir))*p
          fourierrk_t(ir,ik)=fourierrk_t(ir,ik)+dsqrt((fcos_un(ir)+fcos_dn(ir))**2
     &                                               +(fsin_un(ir)+fsin_dn(ir))**2)*p
        enddo
c        write(*,*) 'fourierrk_t=',fourierrk_t(ir,ik)
      enddo

      return
      end

c-------------------------------------------------------------------------------------

      subroutine fourierkk(p,q,xold,xnew)

c Written by A.D.Guclu aug2006.
c Calculates the "internal fourier transform" estimator in cartesian coo., defined as:
c     FT=abs( sum_i Fourier[Delta^2(r_i-r)])
c taking the absolute value before the data accumlation is essential to get the
c "internal structure" of the density. We are indeed calculation the following
c quantity:
c     <FT>=<Psi| abs( Fourier[rho(r)] ) |Psi>
c Notes:
c Delta^2(ri-r)=Delta(xi-x) Delta(yi-y)
c Fourier[Delta]=cos + i sin

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /fourier/ fourierrk_u(0:NAX,0:NAK1),fourierrk_d(0:NAX,0:NAK1)
     &,fourierrk_t(0:NAX,0:NAK1),fourierkk_u(-NAK2:NAK2,-NAK2:NAK2),fourierkk_d(-NAK2:NAK2,-NAK2:NAK2)
     &,fourierkk_t(-NAK2:NAK2,-NAK2:NAK2),delk1,delk2,fmax1,fmax2,ifourier

      dimension xold(3,MELEC),xnew(3,MELEC)

      twopi=2.d0*pi

      do ikx=-NAK2,NAK2
        fkx=delk2*ikx
        do iky=-NAK2,NAK2
          fky=delk2*iky
c first, for old x values:
          sumr_u=0.d0
          sumi_u=0.d0
          sumr_d=0.d0
          sumi_d=0.d0
          sumr_t=0.d0
          sumi_t=0.d0
          do ie=1,nup
            sumr_u=sumr_u+dcos(twopi*(fkx*xold(1,ie)+fky*xold(2,ie)))
            sumi_u=sumi_u+dsin(twopi*(fkx*xold(1,ie)+fky*xold(2,ie)))
          enddo
          do ie=nup+1,nelec
            sumr_d=sumr_d+dcos(twopi*(fkx*xold(1,ie)+fky*xold(2,ie)))
            sumi_d=sumi_d+dsin(twopi*(fkx*xold(1,ie)+fky*xold(2,ie)))
          enddo
          sumr_t=sumr_d+sumr_u
          sumi_t=sumi_d+sumi_u
          fourierkk_u(ikx,iky)=fourierkk_u(ikx,iky)+q*(sumr_u*sumr_u+sumi_u*sumi_u)
          fourierkk_d(ikx,iky)=fourierkk_d(ikx,iky)+q*(sumr_d*sumr_d+sumi_d*sumi_d)
          fourierkk_t(ikx,iky)=fourierkk_t(ikx,iky)+q*(sumr_t*sumr_t+sumi_t*sumi_t)

c now repeat above for new x values:
          sumr_u=0.d0
          sumi_u=0.d0
          sumr_d=0.d0
          sumi_d=0.d0
          sumr_t=0.d0
          sumi_t=0.d0
          do ie=1,nup
            sumr_u=sumr_u+dcos(twopi*(fkx*xnew(1,ie)+fky*xnew(2,ie)))
            sumi_u=sumi_u+dsin(twopi*(fkx*xnew(1,ie)+fky*xnew(2,ie)))
          enddo
          do ie=nup+1,nelec
            sumr_d=sumr_d+dcos(twopi*(fkx*xnew(1,ie)+fky*xnew(2,ie)))
            sumi_d=sumi_d+dsin(twopi*(fkx*xnew(1,ie)+fky*xnew(2,ie)))
          enddo
          sumr_t=sumr_d+sumr_u
          sumi_t=sumi_d+sumi_u
          fourierkk_u(ikx,iky)=fourierkk_u(ikx,iky)+p*(sumr_u*sumr_u+sumi_u*sumi_u)
          fourierkk_d(ikx,iky)=fourierkk_d(ikx,iky)+p*(sumr_d*sumr_d+sumi_d*sumi_d)
          fourierkk_t(ikx,iky)=fourierkk_t(ikx,iky)+p*(sumr_t*sumr_t+sumi_t*sumi_t)

        enddo
      enddo

      return
      end

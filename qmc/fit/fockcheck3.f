      subroutine fockcheck3(o,op,diff,isht,nfsh,it,ipr)
c Written by Cyrus Umrigar and Claudia Filippi

      use atom_mod
      implicit real*8(a-h,o-z)

!JT      include '../vmc/vmc.h'
!JT      include 'fit.h'
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Zfock
      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata
!JT      include '../vmc/force.h'

!JT      parameter(half=0.5d0,two=2.d0,three=3.d0
      parameter(d1b24=1.d0/24.d0,d1b4=0.25d0,d1b6=1.d0/6.d0,d1b12=1.d0/12.d0
     &,rt2=1.414213562373095d0,dln2=0.6931471805599453d0
     &,const1=(1.d0-dln2)/12.d0)
c    &,const1=(1.d0-dln2)/12.d0,const2=-(pi-2.d0)/(6.d0*pi))

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /ncusp/ norbc,ncuspc,nfockc,nfock,ncnstr

      dimension o(2*MORDJ),op(0:2*MORDJ),diff(*)

c f2e o(s^2) from phi20(r12=0)
c f2n o(r^2) from phi20(r1=0)
c f2elog o(s^2 log(s)) from phi21(r12=0)

      if(ifock.eq.4) then
        zfock=znuc(iwctype(it))
        f2n=-d1b6*(eguess+zfock*(zfock-dln2))-d1b24
        f2e=-d1b12*eguess + (const1-d1b4*zfock)*zfock

        o(2)=o(2)+fck(2,it,1)*f2n
        o(2+nord)=o(2+nord)+fck(2,it,1)*f2e+half*a21
       else
        o(2) =o(2)+fck(4,it,1)+fck(5,it,1)+fck(6,it,1)
     &            +fck(7,it,1)+fck(8,it,1)+fck(9,it,1)
        op(1)=op(1)+fck(5,it,1)+three*fck(7,it,1)+fck(8,it,1)
     &             +two*fck(9,it,1)
        o(2+nord) =o(2+nord) +rt2*(half*fck(5,it,1)+fck(7,it,1))
        op(1+nord)=op(1+nord)+rt2*(half*fck(4,it,1)+fck(9,it,1))
      endif

      if(ifock.ne.2) return

      nfsh=nfock/2
      ishn=nord
      ishe=2*nord+nfsh

      alp2n=2*(3*fck(11,it,1)+fck(12,it,1)+2*fck(13,it,1)
     &        -2*fck(14,it,1)-fck(15,it,1))
      alp2e=two*fck(13,it,1)

      diff(isht+ishn+1) = alp2n
      diff(isht+ishe+1) = alp2e

      if(ipr.ge.1) then
        write(6,'(''coefs of r2lg'',f12.6)') diff(isht+ishn+1)
        write(6,'(''coefs of s2lg'',f12.6)') diff(isht+ishe+1)
      endif

      return
      end

      subroutine fockexact3(o,op,it,ipr)
c Written by Cyrus Umrigar and Claudia Filippi

      use atom_mod
      use contr2_mod
      implicit real*8(a-h,o-z)

!JT      include '../vmc/vmc.h'
!JT      include 'fit.h'
!JT      include '../vmc/force.h'

!JT      parameter(zero=0.d0,half=0.5d0,two=2.d0,three=3.d0
      parameter(d1b24=1.d0/24.d0,d1b4=0.25d0,d1b6=1.d0/6.d0,d1b12=1.d0/12.d0
     &,rt2=1.414213562373095d0,dln2=0.6931471805599453d0
     &,pi=3.141592653589793d0
     &,const1=(1.d0-dln2)/12.d0,const2=-(pi-2.d0)/(6.d0*pi))

      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata

!JT      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
!JT     &,ifock,i3body,irewgt,iaver,istrch
!JT     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Zfock

!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord

      dimension o(2*MORDJ),op(0:2*MORDJ)

      if(ifock.eq.4) then
        zfock=znuc(iwctype(it))
        f2n=-d1b6*(eguess+zfock*(zfock-dln2))-d1b24
        f2e=-d1b12*eguess + (const1-d1b4*zfock)*zfock
        o(2)=o(2)+fck(2,it,1)*f2n
        o(2+nord)=o(2+nord)+fck(2,it,1)*f2e
       else
        o(2+nord) =o(2+nord) +rt2*(half*fck(5,it,1)+fck(7,it,1))
        op(1+nord)=op(1+nord)+rt2*(half*fck(4,it,1)+fck(9,it,1))
        o(2) =o(2)+fck(4,it,1)+fck(5,it,1)+fck(6,it,1)
     &            +fck(7,it,1)+fck(8,it,1)+fck(9,it,1)
        op(1)=op(1)+fck(5,it,1)+three*fck(7,it,1)
     &             +fck(8,it,1)+two*fck(9,it,1)
      endif

      fck(11,it,1)=(-fck(12,it,1)-2*fck(13,it,1)+2*fck(14,it,1)
     &    +fck(15,it,1))/three
      fck(13,it,1)=zero

      return
      end

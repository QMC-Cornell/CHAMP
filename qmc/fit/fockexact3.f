      subroutine fockexact3(o,op,it,ipr)
c Written by Cyrus Umrigar and Claudia Filippi
      use constants_mod
      use atom_mod
      use contr2_mod
      use jaspar3_mod
      use pars_mod
      use confg_mod
      implicit real*8(a-h,o-z)

      parameter(d1b24=1.d0/24.d0,d1b4=0.25d0,d1b6=1.d0/6.d0,d1b12=1.d0/12.d0
     &,rt2=1.414213562373095d0,dln2=0.6931471805599453d0
     &,pi=3.141592653589793d0
     &,const1=(1.d0-dln2)/12.d0,const2=-(pi-2.d0)/(6.d0*pi))

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

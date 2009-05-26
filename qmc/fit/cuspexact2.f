      subroutine cuspexact2(a,b)
c Written by Cyrus Umrigar
      use atom_mod
      use dim_mod
      implicit real*8(a-h,o-z)


      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata

      dimension a(*),b(*)

!JT      parameter(half=0.5d0)
      parameter(d1b4=0.25d0,d1b6=1.d0/6.d0,d1b12=1.d0/12.d0
     &,d1b24=1.d0/24.d0,rt2i=0.707106781186547d0
     &,dln2=0.6931471805599453d0, pi=3.141592653589793d0
     &,const1=(1.d0-dln2)/12.d0,const2=-(pi-2.d0)/(6.d0*pi))

c if(nparmb(1)+nparmb(2)+nparmb(3).ne.0) stop 'cuspab implemented so
c &far only ijas=2 and denominator terms=0'

c Warning: temporarily assume only one nucleus type

c f2e    o(s^2) from phi20(r12=0)
c f2elog o(s^2 log(s)) from phi21(r12=0)
c f2n    o(r^2) from phi20(r1=0)
      f2e=-d1b12*eguess +
     &(const1-d1b4*znuc(iwctype(1)))*znuc(iwctype(1))
      f2elog=const2*znuc(iwctype(1))
      f2n=-d1b6*(eguess+znuc(iwctype(1))*(znuc(iwctype(1))-dln2))
     &- d1b24

c b(21), b(20) are determined by b4=0 for e-e and e-n resp.
c The rest are determined by bp=0 for e-e and e-n.
      b(7)=-rt2i*b(35)
      b(14)=0
      b(24)=0
      b(39)=-2*b(50)
      b(5)=b(6) - (b(7)-b(9) + b(36) - b(53)-3*b(54)-2*b(55)
     &-4*b(56)-6*b(57) -3*b(58)-5*b(59)-7*b(60) -b(61)-3*b(62)-5*b(63)
     &-2*b(64)-4*b(65)-6*b(66) -b(67)-3*b(68)-5*b(69))/2
      b(11)=b(12) - (b(13)-b(15) + 2*(b(14)-b(16)) +b(17)-b(18)
     &-b(61)-3*b(62)-5*b(63)-2*b(64)-4*b(65)-6*b(66))/3
      b(21)=0
      b(22)= (3*(b(24)-b(26))
     &+2*(b(27)-b(28) +b(32)-b(33)) +b(23)-b(25) +b(30)-b(31))/4
      b(48)=-(b(49)+2*b(50)-2*b(51)-b(52))/3
      b(20)=0
      do 10 i=21,34
   10   b(20)=b(20)-b(i)
      do 20 i=67,69
   20   b(20)=b(20)-b(i)

      b1=b(2)
      b2=b(5)+rt2i*b(36)+f2e*b(44)
      b3=b(11)
      b4=b(21)
      bl2=b(38)+f2elog*b(40)
      bl3=2*b(48)+half*f2elog*znuc(iwctype(1))*b(45)


      ap1=a(1)

      a(7)=-rt2i*a(35) +ap1*b1
      a(14)=ap1*b2
      a(24)=ap1*b3
      a(39)=-2*a(50)+ap1*bl2

      b1=0
      b2=f2n*b(44)
      b3=0
      b4=0
      bl2=0
      bl3=d1b6*f2elog*b(45)

      do 32 i=1,ndim
   32   b1=b1+b(i)
      do 33 i=4,9
   33   b2=b2+b(i)
      do 34 i=10,19
   34   b3=b3+b(i)
      do 35 i=20,34
   35   b4=b4+b(i)
      do 36 i=35,37
   36   b2=b2+b(i)
      do 37 i=53,60
   37   b2=b2+b(i)
      do 38 i=61,66
   38   b3=b3+b(i)
      do 39 i=67,69
   39   b4=b4+b(i)
      do 40 i=47,52
   40   bl3=bl3+2*b(i)

      ap1=a(2)

      a(5)=a(6) - (a(7)-a(9) + a(36) - a(53)-3*a(54)-2*a(55)
     &-4*a(56)-6*a(57) -3*a(58)-5*a(59)-7*a(60) -a(61)-3*a(62)-5*a(63)
     &-2*a(64)-4*a(65)-6*a(66) -a(67)-3*a(68)-5*a(69)-ap1*b1)/2
      a(11)=a(12) - (a(13)-a(15) + 2*(a(14)-a(16)) +a(17)-a(18)
     &-a(61)-3*a(62)-5*a(63)-2*a(64)-4*a(65)-6*a(66)-ap1*b2)/3
      a(21)=a(22) - (3*(a(24)-a(26))
     &+2*(a(27)-a(28) +a(32)-a(33)) +a(23)-a(25) +a(30)-a(31)
     &-ap1*b3)/4
      a(48)=-(a(49)+2*a(50)-2*a(51)-a(52))/3

      return
      end

      double precision zex(MBASIS,MWF),betaq
      integer n1s(MCTYPE),n2s(MCTYPE),n2p(-1:1,MCTYPE)
      integer n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
      integer n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE)
      integer n4f(-3:3,MCTYPE),n5s(MCTYPE),n5p(-1:1,MCTYPE)
      integer n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
      integer n5g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
      integer nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)

      common /basis/ zex,betaq,n1s,n2s,n2p,n3s,n3p,n3d,n4s,n4p,n4d,n4f,n5s,n5p,n5d,n5f,n5g,n6h,nsa,npa,nda


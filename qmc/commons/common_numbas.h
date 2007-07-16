      double precision arg(MCTYPE),r0(MCTYPE),rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
      double precision d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
      integer numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE)
      integer iwrwf(MBASIS_CTYPE,MCTYPE)

      common /numbas/ arg,r0,rwf,d2rwf,numr,nrbas,igrid,nr,iwrwf

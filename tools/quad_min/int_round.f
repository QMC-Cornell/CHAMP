      do 10 k=1,10
      read(5,*) i
      j=int_round(i)
   10 write(6,*) i,j
      stop
      end
c------------------------------------------
      function int_round(i)
c Find the nearest larger "round" integer
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)

      ipow=10**int((dlog10(dfloat(i))))
      rest=dfloat(i)/ipow
      if(rest.le.2.d0) then
        j=2
       elseif(rest.le.5.d0) then
        j=5
       else
        j=10
      endif

      int_round=j*ipow
   10 write(6,*) i,rest,j,ipow,int_round

      return
      end

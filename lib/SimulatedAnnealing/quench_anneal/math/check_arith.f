      subroutine check_arith
c check_arith: purpose: check overflow handling of integer arithmetic
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::c
c Copyright: H.W.J. Bloete and M.P. Nightingale, June 1995  c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::c
      parameter(MULT=1664525,IADD=1013904223)
      parameter (I0=0,I1=1013904223,I2=1196435762,I3=-775096599,
     &  I4=-1426500812,I5=1649599747)
      iseed=I0
      iseed=iseed*MULT+IADD
      if(iseed.ne.I1) then
        write(6,'(2(1x,z8))')'check_arith: ',iseed,I1
        stop 'check_arith'
      endif
      iseed=iseed*MULT+IADD
      if(iseed.ne.I2) then
        write(6,'(2(1x,z8))')'check_arith: ',iseed,I2
        stop 'check_arith'
      endif
      iseed=iseed*MULT+IADD
      if(iseed.ne.I3) then
        write(6,'(2(1x,z8))')'check_arith: ',iseed,I3
        stop 'check_arith'
      endif
      iseed=iseed*MULT+IADD
      if(iseed.ne.I4) then
        write(6,'(2(1x,z8))')'check_arith: ',iseed,I4
        stop 'check_arith'
      endif
      iseed=iseed*MULT+IADD
      if(iseed.ne.I5) then
        write(6,'(2(1x,z8))')'check_arith: ',iseed,I5
        stop 'check_arith'
      endif
      return
      end

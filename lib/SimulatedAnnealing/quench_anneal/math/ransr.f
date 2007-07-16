cDate of last modification: Thu Aug 29 16:53:45 EDT 1996
      subroutine ransi(iseed)
c ransr: purpose: generate uniform random numbers (fast version)
c shift register random generator with very long period
c combined with linear congruential rn generator
c sequential version
c ransr and ransr_2 differ:
c the first is a shift register xor-ed with a linear congruential
c the second xor-es two shift registers
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::c
c Copyright: H.W.J. Bloete and M.P. Nightingale, June 1995  c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::c
      implicit real*8(a-h,o-z)
c The following parameters for a linear congruential rn generator are from
c Numerical Recipes "An Even Quicker and Dirtier Generator."
c It relies on integer arithmetic that ignores overlows.
c The second set is the same rn generator iterated twice.
      parameter (MULT1=1664525  ,IADD1=1013904223)
      parameter (MULT2=389569705,IADD2=1196435762)
c In k'=mod(k*mult+iadd+2**31,2**32)-2**31
c mod and shift are assumed to be done implicitly by hardware
c the next multiplier is no loner used
c      parameter (mult=32781)

      parameter (TWO=2,TM32=two**(-32),HALF=1/TWO)
      dimension rn(*)
      parameter (LENR=9689,IFDB=471)
      common/ransrb/ irs(LENR),inxt(LENR)
      save ipoint,ipoinf,k
      k=3**18+iseed
      l=MULT1*k+IADD1
c check assumption about integer arithmetic
      call check_arith
c     write(6,*) 'ransi: overflow handling of integer arithmetic is OK'

c initialize shift registers
      do i=1,LENR
        k=k*MULT2+IADD2
        l=l*MULT2+IADD2
        irs(i)=ieor(k,ishft(l,-16))
        inxt(i)=i+1
      enddo

      inxt(LENR)=1
      ipoint=1
      ipoinf=ifdb+1
      return

      entry ransr(rn,n)
c     calculate n random numbers
      do i=1,n
        l=ieor(irs(ipoint),irs(ipoinf))
        k=k*MULT1+IADD1
        rn(i)=ieor(k,l)*TM32+HALF
        irs(ipoint)=l
        ipoint=inxt(ipoint)
        ipoinf=inxt(ipoinf)
      enddo
      return
      end

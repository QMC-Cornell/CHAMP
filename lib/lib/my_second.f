      subroutine my_second (n,title)
c Written by Cyrus Umrigar
c Prints out cpu and wall-clock time,
c both since beginning of run and since last call.

c I used to use etime for cpu time and time for wall-clock time.
c I now use cpu_time and system_clock because they are standard.
c I tested out time and cpu_time for cpu time and
c system_clock and secnds for wall-clock time (system_clock is standard)
c to see if they are more accurate.  All four gave variations of
c a few tenths of a second for a job that ran less than 10 secs.
c cpu_time and time were the same using g77 in Linux
c     real etim
c     call system_clock(icount,icount_rate,icountmax)
c     call cpu_time(etim)
c     t1=secnds(0.)

      implicit real*8 (a-h,o-z)
      real*4 etim
c     real*4 tarray
c     dimension tarray(2)
      character(len=*) title
      integer time
      save icall,itim1,itim2,etim1,etim2
      data icall/0/

c     etim=etime(i)        ! unix
c     etim=etime(tarray)   ! linux
c     etim=mclock()*1.d-2  ! aix
c     etim=mclock()*1.d-6
      call cpu_time(etim)  ! standard but use pgf90, not pgf77

c     itim=time()
c     itim=time(0)         ! aix
      call system_clock(itim,icount_rate,icountmax) ! standard but but use pgf90, not pgf77

      if(icall.eq.0) then
        icall=1
        itim1=itim
        itim2=itim
        etim1=etim
        etim2=etim
      endif
c     itimtot=itim-itim1
c     itimlast=itim-itim2
      itimtot=nint(dfloat(itim-itim1)/icount_rate)
      itimlast=nint(dfloat(itim-itim2)/icount_rate)
      itim2=itim
      etimtot=etim-etim1
      etimlast=etim-etim2
      etim2=etim
!JT      if(n.eq.1) write (6,'(''BEGINNING OF '',a6,'' CP, REAL TIME IS '',
!JT     &2f11.2,2i7)') title,etimtot,etimlast,itimtot,itimlast
!JT      if(n.eq.2) write (6,'(''END       OF '',a6,'' CP, REAL TIME IS '',
!JT     &2f11.2,2i7)') title,etimtot,etimlast,itimtot,itimlast
      if(n.eq.1) write (6,'(a,a,a,f11.2,a,f11.2,a)') 'Beginning of ',title,
     &': total CPU time is',etimtot,' s, CPU time since last check is',etimlast,' s'
      if(n.eq.2) write (6,'(a,a,a,f11.2,a,f11.2,a)') 'End       of ',title,
     &': total CPU time is',etimtot,' s, CPU time since last check is',etimlast,' s'
      return
      end

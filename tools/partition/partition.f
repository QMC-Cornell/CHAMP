c Partition a given angular momentum into N non-negative integers
c    This is needed for constructing determinants for fractional quantum Hall regime.
c    This is the dumb version for 6 particles.  Should be done recursively for any number.

c Print out also:
c 1) the number of partitions, "icount",
c 2) the number of partitions where the largest l_i is <= 2L/N, "icount2"
c    the rough argument for this is that:
c    the average angular momentum per pair is L/((N(N-1)/2),
c    so the angular momentum of a given electron about all N-1 others is 2*L/N.
c    Why this average angular momentum about all others should be related to the
c    maximum angular momentum about the center of the dot is not clear to me.
c    The simpler way of viewing it is just that we restrict the max. 1-body angular
c    momentum to be less than twice the average 1-body angular momentum.
c 3) the number of disconnected l regions "idiscon".
c    The L values for which zero disconnected l regions are possible correspond,
c    at least for N=6, to 1st set of magic numbers of Seki, Kuramoto and Nishino,
c    J. Phys. Soc. Jap. 65, 3945 (1996).
c    My hope is that by using 2 criteria -- max 1-body ang. momentum and max number
c    of disconnected regions, the number of determinants can be kept small.
c 4) The number of partitions for which l_i is <= 2L/N, and idiscon <= 1.

      implicit real*8(a-h,o-z)
      parameter(N=6)
      dimension li(N)

      read(5,*) Lmax
      Lmin=(N*(N-1))/2
      do 20 L=Lmin,Lmax
        icount=0
        icount2=0
        icount3=0
        do 10 l1=0,L
          li(1)=l1
          do 10 l2=l1+1,L
            li(2)=l2
            do 10 l3=l2+1,L
              li(3)=l3
              do 10 l4=l3+1,L
                li(4)=l4
                do 10 l5=l4+1,L
                  li(5)=l5
                  do 10 l6=l5+1,L
                    li(6)=l6
                    ltot=0
                    do 5 i=1,N
    5                  ltot=ltot+li(i)
                    if(ltot.eq.L) then
                      icount=icount+1
                      if(li(N).le.(2*L)/N) icount2=icount2+1
                      idiscon=0
                      do 7 i=1,N-1
    7                   if(li(i+1).ne.li(i)+1) idiscon=idiscon+1
                      if(li(N).le.(2*L)/N .and. idiscon.le.1) icount3=icount3+1
                      write(6,'(''L, l_i'',i4,x,6i3,i5)')
     &                L,(li(i),i=1,N),idiscon
                    endif
   10    continue
   20    write(6,'(''L,icount,icount2,icount3='',i2,9i5,/)') L,icount,icount2,icount3
       stop
       end

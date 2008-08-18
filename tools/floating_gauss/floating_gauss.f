c Used for 2 purposes:
c Initial positions of floating gaussians for Devrim's quantum rings
c Reduce existing positions to existing positions [0,2pi)
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)
      parameter(MELEC=1000)
      dimension angles(MELEC)

      twopi=8*datan(1.d0)

      read(5,*) nelec,nup
      write(6,'(''# of total, up-spin and dn-spin electrons='',3i4)') nelec,nup,ndn
      ndn=nelec-nup
      
c Reduce existing positions to existing positions [0,2pi)
      read(5,*) (angles(i),i=1,nelec)
      write(6,'(/,''Original positions and ones modulo 2pi'')')
      do 10 i=1,nelec
        write(6,'(i2,9f9.4)') i,angles(i),angles(i)-(floor(angles(i)/twopi)*twopi)
   10   angles(i)=angles(i)-(floor(angles(i)/twopi)*twopi)
      write(6,'(100f11.8)') (angles(i)-(floor(angles(i)/twopi)*twopi),i=1,nelec)

c Sort the up and the down separately assuming that the first nup are up.
      write(6,'(/,''Positions sorted separately for up and down'')')
      call shell(angles(1),nup)
      call shell(angles(nup+1),ndn)
      do 15 i=1,nelec
   15   write(6,'(i2,9f9.4)') i,angles(i)
      write(6,'(100f11.8)') (angles(i),i=1,nelec)

c Sort all nelec electrons
      write(6,'(/,''Positions sorted for all electrons'')')
      call shell(angles,nelec)
      do 17 i=1,nelec
   17   write(6,'(i2,9f9.4)') i,angles(i)
      write(6,'(100f11.8)') (angles(i),i=1,nelec)

c Now write out odd positions and then even positions so that the up and dn electrons are alternating
      write(6,'(/,''Odd sorted positions followed by even sorted positions'')')
      write(6,'(100f11.8)') (angles(2*i-1),i=1,nup),(angles(2*i),i=1,ndn)

c Symmetrize them
c Symmetrized sorted positions
      write(6,'(/,''Symetrized positions sorted for all electrons'')')
      angles(1)=0
      if(mod(nelec,2).eq.0) angles(nelec/2+1)=.5d0*twopi
      do 30 i=2,(nelec+1)/2
        angles(i)=.5d0*(angles(i)+(twopi-angles(nelec+2-i)))
   30   angles(nelec+2-i)=twopi-angles(i)
      write(6,'(100f11.8)') (angles(i),i=1,nelec)

c Now write out odd positions and then even positions so that the up and dn electrons are alternating
      write(6,'(/,''Symmetrized odd sorted positions followed by even sorted positions'')')
      write(6,'(100f11.8)') (angles(2*i-1),i=1,nup),(angles(2*i),i=1,ndn)

c Now create evenly spaced electrons instead of those read in
c If we say the ~N/2 up-spin electrons are in the first N/2 gaussians then we need the foll:
      write(6,'(/,''Antiferromag. even spacing if the nup up-spin electrons are in the first nup orbs'')')
      do 50 i=1,nelec
        if(i.le.(nelec+1)/2) then
          angles(i)=2*(i-1)*twopi/nelec
         else
          angles(i)=2*(i-(nelec+1)/2-0.5d0)*twopi/nelec
        endif
   50 write(6,'(i2,f9.4)') i,angles(i)
      write(6,'(100f11.8)') (angles(i),i=1,nelec)

c If we say the ~N/2 up-spin electrons are in the odd gaussians
c and the ~N/2 dn-spin in the even gaussians then we need the foll:
      write(6,'(/,''Antiferromag. even spacing if the nup up-spin electrons are in the odd orbs'')')
      do 60 i=1,nelec
        angles(i)=(i-1)*twopi/nelec
   60 write(6,'(i2,f9.4)') i,angles(i)
      write(6,'(100f11.8)') (angles(i),i=1,nelec)

      stop
      end
c-----------------------------------------------------------------------------------------
      SUBROUTINE SHELL(D,N)
      IMPLICIT REAL*8(A-H,O-Z)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:::  SHELL-METZGER SORT IN ASCENDING ORDER.          ...CYRUS 1979  :::
C:::  MODIFIED SLIGHTLY FOR READIBILITY.          ...CYRUS 7 DEC 83  :::
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      DIMENSION D(N)
      LOGNB2=INT(DLOG(DFLOAT(N))/DLOG(2.D0)+1.D-14)
      M=N
      DO 20 NN=1,LOGNB2
        M=M/2
        K=N-M
        DO 20 J=1,K
          DO 10 I=J,1,-M
            L=I+M
            IF (D(L).GT.D(I)) GOTO 20
            T=D(I)
            D(I)=D(L)
            D(L)=T
   10       CONTINUE
   20     CONTINUE
      RETURN
      END

      program change_basis
c Written by Cyrus Umrigar
c Given atom positions for one set of lattice generators, find the positions in another
c set that generates identically the same lattice (the same lattice rotated is not
c good enough.  This is needed because qmc program expects lattice generators to be
c the smallest possible lattice vectors.
c R = CA = BD, so B=CAB^-1
c Also, given reciprocal lattice vectors in the first basis, find those in
c the other basis.  Note that if H is the reciprocal lattice of B, then
c H^-1 is B^transpose/(2*pi).
c Finally print out the angles between the lattice vectors.
c Warning: I check to make sure that the 2 sets of lattice generators yield
c the same volume for the cell, but I should really check to see that they
c yield the same lattice.  The same lattice rotated is NOT good enough since
c the atom positions and k-vectors will then be wrong.  So, I should modify
c program to deduce B lattice generators, from the A ones, such that the B are
c the smallest possible, using code in QMC program.

      implicit real*8(a-h,o-z)

      dimension a(3,3),a_inv(3,3),b(3,3),b_inv(3,3),c(3),d(3),r(3)
     &,g(3,3),h(3,3)

      twopi=8*datan(1.d0)

      read(5,*) alattice
      read(5,*) ((a(i,k),k=1,3),i=1,3)
      read(5,*) blattice
      read(5,*) ((b(i,k),k=1,3),i=1,3)

      do 5 i=1,3
        do 5 k=1,3
          a(k,i)=a(k,i)*alattice
          b(k,i)=b(k,i)*blattice
          a_inv(k,i)=a(k,i)
    5     b_inv(k,i)=b(k,i)

      call matinv(a_inv,3,deta)
      call matinv(b_inv,3,detb)

      write(6,'(''volume or a,b='',2f12.6)') deta,detb
c     if(abs(deta-detb).gt.1.d-4) stop 'volumes of 2 lattices should be the same'
      if(abs(deta-detb).gt.1.d-4) write(6,'(''Warning: volumes of 2 lattices are not the same'')')

      read(5,*) ncent
      do 90 iatom=1,ncent
        read(5,*) c
        write(6,'(/,''atom position in old basis='',(3f14.9))') (c(k),k=1,3)
        do 10 k=1,3
          r(k)=0
          do 10 i=1,3
   10       r(k)=r(k)+c(i)*a(i,k)
        write(6,'(''atom position in cartesian='',(3f14.9))') (r(k),k=1,3)
        do 20 k=1,3
          d(k)=0
          do 20 i=1,3
   20       d(k)=d(k)+r(i)*b_inv(i,k)

        write(6,'(''atom position in new basis='',(3f14.9))') (d(k),k=1,3)

c Check we did it right
        do 30 k=1,3
          r(k)=0
          do 30 i=1,3
   30       r(k)=r(k)+d(i)*b(i,k)
   90   write(6,'(''atom position in cartesian='',(3f14.9))') (r(k),k=1,3)

      deta1=twopi/deta
      g(1,1)=deta1*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      g(2,1)=deta1*(a(3,2)*a(1,3)-a(3,3)*a(1,2))
      g(3,1)=deta1*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      g(1,2)=deta1*(a(2,3)*a(3,1)-a(2,1)*a(3,3))
      g(2,2)=deta1*(a(3,3)*a(1,1)-a(3,1)*a(1,3))
      g(3,2)=deta1*(a(1,3)*a(2,1)-a(1,1)*a(2,3))
      g(1,3)=deta1*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      g(2,3)=deta1*(a(3,1)*a(1,2)-a(3,2)*a(1,1))
      g(3,3)=deta1*(a(1,1)*a(2,2)-a(1,2)*a(2,1))

      detb1=twopi/detb
      h(1,1)=detb1*(b(2,2)*b(3,3)-b(2,3)*b(3,2))
      h(2,1)=detb1*(b(3,2)*b(1,3)-b(3,3)*b(1,2))
      h(3,1)=detb1*(b(1,2)*b(2,3)-b(1,3)*b(2,2))
      h(1,2)=detb1*(b(2,3)*b(3,1)-b(2,1)*b(3,3))
      h(2,2)=detb1*(b(3,3)*b(1,1)-b(3,1)*b(1,3))
      h(3,2)=detb1*(b(1,3)*b(2,1)-b(1,1)*b(2,3))
      h(1,3)=detb1*(b(2,1)*b(3,2)-b(2,2)*b(3,1))
      h(2,3)=detb1*(b(3,1)*b(1,2)-b(3,2)*b(1,1))
      h(3,3)=detb1*(b(1,1)*b(2,2)-b(1,2)*b(2,1))

      read(5,*) nkvec
      do 190 ikvec=1,nkvec
        read(5,*) c
        write(6,'(/,''k-vector in old basis='',(3f7.3))') (c(k),k=1,3)
        do 110 k=1,3
          r(k)=0
          do 110 i=1,3
  110       r(k)=r(k)+c(i)*g(i,k)
        write(6,'(''k-vector in cartesian='',(3f15.9))') (r(k),k=1,3)
        do 120 k=1,3
          d(k)=0
          do 120 i=1,3
  120       d(k)=d(k)+r(i)*b(k,i)/twopi

        write(6,'(''k-vector in new basis='',(3f7.3))') (d(k),k=1,3)

c Check we did it right
        do 130 k=1,3
          r(k)=0
          do 130 i=1,3
  130       r(k)=r(k)+d(i)*h(i,k)
  190   write(6,'(''k-vector in cartesian='',(3f15.9))') (r(k),k=1,3)

      do 210 i=2,3
        do 210 j=1,i-1
          dota=0
          absai=0
          absaj=0
          dotb=0
          absbi=0
          absbj=0
          do 200 k=1,3
            absai=absai+a(i,k)**2
            absaj=absaj+a(j,k)**2
            dota=dota+a(i,k)*a(j,k)
            absbi=absbi+b(i,k)**2
            absbj=absbj+b(j,k)**2
  200       dotb=dotb+b(i,k)*b(j,k)
  210     write(6,'(''angle between vectors'',2i2,'' of a,b lattice='',2f8.3)') j,i
     &    ,(360/twopi)*dacos(dota/sqrt(absai*absaj))
     &    ,(360/twopi)*dacos(dotb/sqrt(absbi*absbj))

      do 310 i=1,3
        absa=0
        absb=0
        do 300 k=1,3
          absa=absa+a(i,k)**2
  300     absb=absb+b(i,k)**2
  310   write(6,'(''length of'',i2,'' a,b lattice vector='',2f9.6)') i,sqrt(absa),sqrt(absb)

      stop
      end

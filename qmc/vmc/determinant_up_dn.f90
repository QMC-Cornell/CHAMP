      subroutine determinant_up_dn
! Written by Cyrus Umrigar
      use all_tools_mod
      use dorb_mod
      use dets_mod
      implicit real*8(a-h,o-z)

      call alloc ('iwdetup', iwdetup, ndet)
      call alloc ('iwdetdn', iwdetdn, ndet)

      call alloc ('iworbdup', iworbdup, nup, 1)
      call alloc ('iworbddn', iworbddn, ndn, 1)

      ndetup=0
      ndetdn=0
      if(nup.ge.1) then
        ndetup=1
        iwdetup(1)=1
        call alloc ('iworbdup', iworbdup, nup, 1)
        do 10  iup=1,nup
   10     iworbdup(iup,1)=iworbd(iup,1)
      endif

      if(ndn.ge.1) then
        ndetdn=1
        iwdetdn(1)=1
        call alloc ('iworbddn', iworbddn, ndn, 1)
        do 20  idn=1,ndn
   20     iworbddn(idn,1)=iworbd(nup+idn,1)
      endif

      do 80 idet=2,ndet
        do 40 idetup=1,ndetup
          do 30 iup=1,nup
   30       if(iworbd(iup,idet).ne.iworbdup(iup,idetup)) goto 40
        iwdetup(idet)=idetup
        do 35 iup=1,nup
   35     iworbdup(iup,idetup)=iworbd(iup,idet)
        goto 50
   40   continue
        ndetup=ndetup+1
        iwdetup(idet)=ndetup
        call alloc ('iworbdup', iworbdup, nup, ndetup)
        do 45 iup=1,nup
   45     iworbdup(iup,ndetup)=iworbd(iup,idet)

   50   do 70 idetdn=1,ndetdn
          do 60 idn=1,ndn
   60       if(iworbd(nup+idn,idet).ne.iworbddn(idn,idetdn)) goto 70
        iwdetdn(idet)=idetdn
        do 65 idn=1,ndn
   65     iworbddn(idn,idetdn)=iworbd(nup+idn,idet)

        goto 80
   70   continue
        ndetdn=ndetdn+1
        iwdetdn(idet)=ndetdn
        call alloc ('iworbddn', iworbddn, ndn, ndetdn)
        do 75 idn=1,ndn
   75     iworbddn(idn,ndetdn)=iworbd(nup+idn,idet)


   80   continue
!JT   80 write(6,'(''idet,ndetup,iwdetup(idet),ndetdn,iwdetdn(idet)'',9i5)') idet,ndetup,iwdetup(idet),ndetdn,iwdetdn(idet)


      write(6,'(a,i5)') ' number of unique spin-up   determinants = ', ndetup
      write(6,'(a,i5)') ' number of unique spin-down determinants = ', ndetdn
      write(6,'(a)') ' unique spin-up determinants have orbitals:'
      do 90 idetup=1,ndetup
   90   write(6,'(a,i5,a,9000i4)') ' det # ',idetup, ': ',(iworbdup(iup,idetup),iup=1,nup)
      write(6,'(a)') ' unique spin-down determinants have orbitals:'
      do 95 idetdn=1,ndetdn
   95   write(6,'(a,i5,a,9000i4)') ' det # ',idetdn, ': ',(iworbddn(idn,idetdn),idn=1,ndn)

!     JT: warning: quick (and dirty) fix for dealing with the case ndn = 0
!     ndetdn is reset to 1, and the corresponding determinants will just have the value 1
!     this is useful for orbital optimization
      if (ndetdn .eq. 0) then
        ndetdn = 1
        call alloc ('iworbddn', iworbddn, ndn, ndetdn)
        do idet=1,ndet
         iwdetdn(idet)=1
        enddo
        write(6,'(a,i1)') ' Warning: no spin-down determinants, but ndetdn is reset to ',ndetdn
      endif

! JT: introduce variable ndetupdn replacing MDETUD
! JT: I am not sure if we really need ndetupdn
      ndetupdn = max(ndetup, ndetdn)
      call object_modified ('ndetup')
      call object_modified ('ndetdn')
      call object_modified ('ndetupdn')
      call object_modified ('iworbdup')
      call object_modified ('iworbddn')
      call object_modified ('iwdetup')
      call object_modified ('iwdetdn')

      return
      end

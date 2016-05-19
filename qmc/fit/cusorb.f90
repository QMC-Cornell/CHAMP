      subroutine cusorb(icent,orb)
! Written by Cyrus Umrigar
! Calculate all orbitals, orb, at position of nucleus icent
! J. Toulouse - 08 Jan 05: change coef(i,j,1) -> coef(i,j,iwf)
      use all_tools_mod
      use atom_mod
      use coefs_mod
      use dim_mod
      use wfsec_mod
      use phifun_mod
      use const_mod
      implicit real*8(a-h,o-z)

      parameter(eps=1.d-99)

      dimension orb(*),rvec_en(3,nelec,ncent),r_en(nelec,ncent)

!JT    iwf=1

      do 6 ic=1,ncent
        r_en(1,ic)=0
        do 5 k=1,ndim
          rvec_en(k,1,ic)=cent(k,icent)-cent(k,ic)
    5     r_en(1,ic)=r_en(1,ic)+rvec_en(k,1,ic)**2
        r_en(1,ic)=sqrt(r_en(1,ic))
        if(r_en(1,ic).eq.0.d0) then
          r_en(1,ic)=eps
          rvec_en(1,1,ic)=eps
        endif
    6 continue

! allocate basis functions
      call alloc ('phin', phin, nbasis, nelec)
      call alloc ('dphin', dphin, 3, nbasis, nelec)
      call alloc ('d2phin', d2phin, nbasis, nelec)

! Calculate basis functions, phin, at position of nucleus icent
      if(ndim.eq.3) then
        call basis_fnse2(1,rvec_en,r_en)
       elseif(ndim.eq.2) then
        call basis_fns_2de2(1,rvec_en,r_en)
      endif

!     write(6,'(''r_en='',9f9.5)') (r_en(1,ic),ic=1,ncent)
!     write(6,'(''phin='',9f19.15)') (phin(ib,1),ib=1,nbasis)

! Calculate orbitals, orb, at position of nucleus icent
      do 10 iorb=1,norb
        orb(iorb)=0
        do 10 m=1,nbasis
   10     orb(iorb)=orb(iorb)+coef(m,iorb,iwf)*phin(m,1)

      return
      end

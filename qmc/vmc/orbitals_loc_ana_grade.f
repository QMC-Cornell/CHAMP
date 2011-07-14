      subroutine orbitals_loc_ana_grade(iel,rvec_en,r_en,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c Calculate localized orbitals and derivatives for all or 1 electrons
      use control_mod
      use coefs_mod
      use dim_mod
      use wfsec_mod
      use contrl_per_mod
      use phifun_mod
      use const_mod
      use atom_mod
      implicit real*8(a-h,o-z)

      dimension rvec_en(3,nelec,ncent),r_en(nelec,ncent)
     &,orb(norb),dorb(3,norb),ddorb(norb)

c get basis functions
      if(ndim.eq.3) then
        call basis_fns(iel,rvec_en,r_en)
       elseif(ndim.eq.2) then
         if(ibasis.eq.1) then
           call basis_fns_2d(iel,rvec_en,r_en)
         elseif(ibasis.eq.4) then
           call basis_fns_2dgauss(iel,rvec_en,r_en)
         elseif(ibasis.eq.5) then
           call basis_fns_polargauss(iel,rvec_en,r_en)
         elseif(ibasis.eq.6) then
           call basis_fns_2dgauss_noncirc(iel, rvec_en,r_en)
         elseif(ibasis.eq.7) then
           if(iper_gaussian_type.eq.1) then 
             call basis_fns_2dgauss_periodic(iel,rvec_en,r_en)
           else 
             call basis_fns_2dgauss_periodic2(iel,rvec_en,r_en)
           endif
         else
           stop 'orbitals_loc_ana: ibasis must be 1,4,5, 6, or 7 for 2d systems'
         endif
      endif

      do 25 iorb=1,norb
          orb(iorb)=0
          dorb(1,iorb)=0
          dorb(2,iorb)=0
          dorb(3,iorb)=0
          ddorb(iorb)=0
          do 25 m=1,nbasis
           orb(iorb)=orb(iorb)+coef(m,iorb,iwf)*phin(m,iel)
           dorb(1,iorb)=dorb(1,iorb)+coef(m,iorb,iwf)*dphin(1,m,iel)
           dorb(2,iorb)=dorb(2,iorb)+coef(m,iorb,iwf)*dphin(2,m,iel)
           dorb(3,iorb)=dorb(3,iorb)+coef(m,iorb,iwf)*dphin(3,m,iel)
   25      ddorb(iorb)=ddorb(iorb)+coef(m,iorb,iwf)*d2phin(m,iel)

      return
      end

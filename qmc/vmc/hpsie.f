      subroutine hpsie(iel,coord,psid,psij,velocity)
      use control_mod
c Written by Claudia Filippi by modifying hpsi
c Calculates wavefunction and velocity contributions for electron iel

      use const_mod
      use dim_mod
      use contr2_mod
      use wfsec_mod
      use contrl_per_mod
      use jaspar4_mod
      use contr3_mod
      use distance_mod
      implicit real*8(a-h,o-z)

      dimension coord(3,*),velocity(3,*),vj(3,nelec),vd(3,nelec)

      iwf=iwftype(1)

c Warning: When one has more than 1 walker and one is including een terms
c in J, then we have to call distances rather than distancese because
c r_en of electrons, other than the moved electron are not correctly
c restored because the array in which they are saved does not have MWALK.
c These are needed only if nordc.gt.0.
c Note that calling distances from the periodic code is a bad idea since it
c changes the scaling of the code with number of electrons.  However, at present
c we are not using een terms for the periodic code, so this is not a problem.
c Also, once I fix it, possibly in hpsi there is no need to call distances
c (just calculate pe instead) if nforce=1.
      if(index(mode,'dmc').ne.0 .and. nordc.gt.0) then
        call distances(coord,pe,pei)   ! this is needed unless distances is called in the dmc routine
      else
        call distancese(iel,coord)
      endif
c     write(6,'(''r_en1='',9d12.5)') ((r_en(i,j),i=1,nelec),j=1,1,ncent)
c     write(6,'(''r_ee1='',9d12.5)') (r_ee(ij),ij=1,nelec*(nelec-1)/2)

c get contribution from jastrow
      if(ianalyt_lap.eq.1) then
        call jastrowe(iel,coord,vj,psij)
       else
        stop 'numerical one-electron move not implemented'
      endif

c get contribution from determinants
      if(ibasis.eq.3) then
c complex case:
        call cdeterminante(iel,coord,rvec_en,r_en,vd,psid)
      else
        call determinante(iel,coord,rvec_en,r_en,vd,psid)
      endif
c combine to get derivatives of full wavefunction
c Gradient of ln(JD) is the sum of the gradient of ln(J) and ln(D).
      do 10 i=1,nelec
        do 10 k=1,ndim
          velocity(k,i)=vj(k,i)+vd(k,i)
c         write (6,*) 'in hpsie'
c         write (6,*) i, velocity(k,i), vj(k,i), vd(k,i)
   10     continue
      return
      end

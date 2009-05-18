      subroutine hpsie(iel,coord,psid,psij,velocity)
      use control_mod
c Written by Claudia Filippi by modifying hpsi
c Calculates wavefunction and velocity contributions for electron iel

      use const_mod
      implicit real*8(a-h,o-z)
      character*16 mode

!JT      include 'vmc.h'
!JT      include 'pseudo.h'
!JT      include 'force.h'

      common /dim/ ndim
      common /contrl_per/ iperiodic,ibasis
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode
c Warning: jaspar4 put in just to have nordc for temporary fix
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      dimension coord(3,*),velocity(3,*),vj(3,MELEC),vd(3,MELEC)

      iwf=iwftype(1)

c Warning: When one has more than 1 walker and one is including een terms
c in J, then we have to call distances rather than distancese because
c r_en of electrons, other than the moved electron are not correctly
c restored because the array in which they are saved does not have MWALK.
c Note that calling distances from the periodic code is a bad idea since it
c changes the scaling of the code with number of electrons.  However, at present
c we are not using een terms for the periodic code, so this is not a problem.
c Also, once I fix it, possibly in hpsi there is no need to call distances
c (just calculate pe instead) if nforce=1.
      if(index(mode,'dmc').ne.0 .and. nordc.gt.0) then
        call distances(coord,pe,pei)
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
   10     velocity(k,i)=vj(k,i)+vd(k,i)

      return
      end

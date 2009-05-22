      subroutine determinant(x,rvec_en,r_en,ddet_det,d2lndet,div_vd,determ)
c Written by Cyrus Umrigar starting from Kevin Schmidt's routine
      use all_tools_mod
      use control_mod
      use deriv_csf_mod
      use orb_mod
      use dorb_mod
      use coefs_mod
      use dets_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      implicit real*8(a-h,o-z)

!JT      parameter (one=1.d0,half=0.5d0)
      character*16 mode

c Routine to calculate the value, gradient and Laplacian of the
c determinantal part of the wavefunction.
c Also derivatives wrt csf_coefs for optimizing them.

c Naming conventions are:
c "d" at beginning of name indicates gradient
c "d2" at beginning of name indicates Laplacian
c "i" at end of portion of name indicates derivative wrt. csf_coef(i).
c "e" at end of portion of name indicates that it is not summed over all electrons
c "ln" indicates natural logarithm
c k=cartesian component, i=electron index
c determ        =  sum of determinants
c ddet_det(k,i) =  gradient(ln(det))
c d2lndet       =  Laplacian(ln(det)) = Laplacian(det)/det - (gradient(det)/det)**2
c div_vd(i)     =  Laplacian(ln(det)) = Laplacian(det)/det - (gradient(det)/det)**2 for each electron
c ekinen(i)     =  -0.5 Laplacian(det)/det

c Note that the first dimension of the slater matrices is MMAT_DIM = (MELEC/2)**2.
c The first dimension of the Slater matrices must be at least max(nup**2,ndn**2)
c So, we check in read_input that nup and ndn are each <= MELEC/2.

!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode
      common /contrl_opt2/ igradhess,iadd_diag_opt
      common /contrl_per/ iperiodic,ibasis
      common /contrl_opt/ nparm,nsig,ncalls,iopt,ipr_opt
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
c     common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
c    &,d2phin(MBASIS,MELEC)
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /kinet/ ekineo(MELEC),ekinen(MELEC)
!JT      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
!JT     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
c Note: d2edeti_deti(MELEC,MDET) need not be in common
!JT      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
!JT     &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
!JT     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
!JT     &,detu(MDETUD),detd(MDETUD)
!JT     &,ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
!JT     &,detij_det(MPARMD,MPARMD)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /dojasderiv/ ijasderiv
!JT     common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype

!JT      common /orb/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB)  !JT

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),ddet_det(3,*),div_vd(MELEC)
      dimension dporb(MOTYPE,MELEC,MORB),d2porb(MOTYPE,MOTYPE,MELEC,MORB)
      dimension ddporb(3,MOTYPE,MELEC,MORB),d2dporb(MOTYPE,MELEC,MORB)

c initialize the derivative arrays to zero
      do 10 i=1,nelec
        ekinen(i)=0
        ddet_det(1,i)=0
        ddet_det(2,i)=0
   10   ddet_det(3,i)=0
      d2lndet=0

c initialize the determinant arrays to one
      do 20 idet=1,ndetup
   20   detu(idet)=one
      do 22 idet=1,ndetdn
   22   detd(idet)=one
c     determ=0

c get orbitals and derivatives for all electrons
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          if(nparmot.eq.0 .or. igradhess.eq.0) then
            call orbitals_loc_ana(0,rvec_en,r_en,orb,dorb,ddorb)
           else
            call deriv_orbitals(rvec_en,r_en,orb,dorb,ddorb,dporb,d2porb,ddporb,d2dporb)
          endif
         else
          call orbitals_loc_num(0,x,orb,dorb,ddorb)
        endif

       else

        if(inum_orb.eq.0) then
          call orbitals_pw(x,orb,dorb,ddorb)
c         write(6,'(''x(1,1),orb(1,1) from pw'',9f12.8)') x(1,1),orb(1,1),dorb(1,1,1),ddorb(1,1)
         else
          call orbitals_period_num(x,orb,dorb,ddorb)
c         write(6,'(''x(1,1),orb(1,1) from nu'',9f12.8)') x(1,1),orb(1,1),dorb(1,1,1),ddorb(1,1)
        endif

      endif

      call object_modified_by_index (orb_index)   !JT
      call object_modified_by_index (dorb_index)   !JT
      call object_modified_by_index (ddorb_index)   !JT

      if(ipr.ge.4) then
        do 26 iorb=1,norb
          write(6,'(''iorb,orb='',i3,(30f9.5))') iorb,(orb(i,iorb),i=1,nelec)
          write(6,'(''iorb,dorb1='',i3,(30f9.5))') iorb,(dorb(1,i,iorb),i=1,nelec)
          write(6,'(''iorb,dorb2='',i3,(30f9.5))') iorb,(dorb(2,i,iorb),i=1,nelec)
          write(6,'(''iorb,dorb3='',i3,(30f9.5))') iorb,(dorb(3,i,iorb),i=1,nelec)
   26     write(6,'(''iorb,ddorb='',i3,(30f9.5))') iorb,(ddorb(i,iorb),i=1,nelec)
      endif

c !fp
c Not only are we looking at the transpose of the Slater,
c but also, we fold the the matrix into a single column
c example: the following slater matrix
c ( phi_1(r1) phi_1(r2) )
c |                     |
c ( phi_2(r1) phi_2(r2) )
c is represented as
c ( phi_1(r1) )
c | phi_1(r2) |
c | phi_2(r1) |
c ( phi_2(r2) )
c !fp
c The 3 nested loops are over
c 1) up electrons,
c 2) determinants
c 3) basis states setting up transpose of the Slater
c matrix in slmui to get inverse transpose.
c Also put derivatives in fpu and fppu.
c iworbdup which orbitals enter which determinants
c iworbdup(i,j) is ith orbitals of jth determinant

c fpu example: (for a idet=1)
c ( dphi_1(r1) / dx  dphi_2(r1) / dx  dphi_1(r2) / dx  dphi_2(r2) / dx )
c | dphi_1(r1) / dy  dphi_2(r1) / dy  dphi_1(r2) / dy  dphi_2(r2) / dy |
c ( dphi_1(r1) / dz  dphi_2(r1) / dz  dphi_1(r2) / dz  dphi_2(r2) / dz )

c fppu example: (for idet=1)
c ( d^2 phi_1(r1) / dx^2 + d^2 phi_1(r1) / dy^2 + d^2 phi_1(r1) / dz^2 )
c ( d^2 phi_2(r1) / dx^2 + d^2 phi_2(r1) / dy^2 + d^2 phi_2(r1) / dz^2 )
c ( d^2 phi_1(r2) / dx^2 + d^2 phi_1(r2) / dy^2 + d^2 phi_1(r2) / dz^2 )
c ( d^2 phi_2(r2) / dx^2 + d^2 phi_2(r2) / dy^2 + d^2 phi_2(r2) / dz^2 )
c !fp
      ik=-nup
      do 30 i=1,nup
        ik=ik+nup
        do 30 idet=1,ndetup
          jk=-nup
          do 30 j=1,nup
            jk=jk+nup
            slmui(i+jk,idet)=orb(i,iworbdup(j,idet))
            fpu(1,j+ik,idet)=dorb(1,i,iworbdup(j,idet))
            fpu(2,j+ik,idet)=dorb(2,i,iworbdup(j,idet))
            fpu(3,j+ik,idet)=dorb(3,i,iworbdup(j,idet))
   30       fppu(j+ik,idet)=ddorb(i,iworbdup(j,idet))

c loop through number of determinants calculating the inverse
c transpose matrices and their determinants
      do 40 idet=1,ndetup
c       write(6,'(''before matinv 40'',i2,9d14.6)') idet,slmui(1,idet),slmui(2,idet),slmui(3,idet),slmui(4,idet)
   40   call matinv(slmui(1,idet),nup,detu(idet))
c  40   write(6,'(''after matinv 40'',i2,9d14.6)') idet,slmui(1,idet),slmui(2,idet),slmui(3,idet),slmui(4,idet)

      call object_modified_by_index (slmui_index) !JT
      call object_modified_by_index (detu_index)  !JT

c repeat above for down spins
      ik=-ndn
      do 50 i=1,ndn
        ik=ik+ndn
        do 50 idet=1,ndetdn
          jk=-ndn
          do 50 j=1,ndn
            jk=jk+ndn
            slmdi(i+jk,idet)=orb(i+nup,iworbddn(j,idet))
            fpd(1,j+ik,idet)=dorb(1,i+nup,iworbddn(j,idet))
            fpd(2,j+ik,idet)=dorb(2,i+nup,iworbddn(j,idet))
            fpd(3,j+ik,idet)=dorb(3,i+nup,iworbddn(j,idet))
   50       fppd(j+ik,idet)=ddorb(i+nup,iworbddn(j,idet))

      do 60 idet=1,ndetdn
   60   call matinv(slmdi(1,idet),ndn,detd(idet))

      call object_modified_by_index (slmdi_index)  !JT
      call object_modified_by_index (detd_index)  !JT

      if(ipr.ge.4) write(6,'(''detu,detd'',9d12.5)') detu(1),detd(1)

c set up sum of slater determinants along with their
c coefficients that were included in input data
      do 70 idet=1,max(ndetup,ndetdn)
c       determ=determ+detu(idet)*detd(idet)*cdet(idet,iwf)

c zero out temporary derivative arrays
        do 70 i=1,nelec
          ddeti_deti(1,i,idet)=0
          ddeti_deti(2,i,idet)=0
          ddeti_deti(3,i,idet)=0
  70      d2edeti_deti(i,idet)=0

c loop through up spin electrons
c take inner product of transpose inverse with derivative
c vectors to get (1/detup)*d(detup)/dx and (1/detup)*d2(detup)/dx**2
      do 80 idet=1,ndetup
        ik=-nup
        do 80 i=1,nup
          ik=ik+nup
          do 80 j=1,nup
            ddeti_deti(1,i,idet)=ddeti_deti(1,i,idet)+slmui(j+ik,idet)*fpu(1,j+ik,idet)
            ddeti_deti(2,i,idet)=ddeti_deti(2,i,idet)+slmui(j+ik,idet)*fpu(2,j+ik,idet)
            ddeti_deti(3,i,idet)=ddeti_deti(3,i,idet)+slmui(j+ik,idet)*fpu(3,j+ik,idet)
   80       d2edeti_deti(i,idet)=d2edeti_deti(i,idet)+slmui(j+ik,idet)*fppu(j+ik,idet)

c repeat above for down spins
      do 90 idet=1,ndetdn
        ik=-ndn
        do 90 i=nup+1,nelec
          ik=ik+ndn
          do 90 j=1,ndn
            ddeti_deti(1,i,idet)=ddeti_deti(1,i,idet)+slmdi(j+ik,idet)*fpd(1,j+ik,idet)
            ddeti_deti(2,i,idet)=ddeti_deti(2,i,idet)+slmdi(j+ik,idet)*fpd(2,j+ik,idet)
            ddeti_deti(3,i,idet)=ddeti_deti(3,i,idet)+slmdi(j+ik,idet)*fpd(3,j+ik,idet)
   90       d2edeti_deti(i,idet)=d2edeti_deti(i,idet)+slmdi(j+ik,idet)*fppd(j+ik,idet)

c combine results for up and down spins to get d(det)/dx
c and d2(det)/dx in ddet_det and d2lndet respectively
      determ=0
      do 117 icsf=1,ncsf
        do 117 idet_in_csf=1,ndet_in_csf(icsf)
          idet=iwdet_in_csf(idet_in_csf,icsf)
c         term=detu(idet)*detd(idet)*cdet(idet,iwf)
          if(ndn.ge.1) then
            term=detu(iwdetup(idet))*detd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
           else
            term=detu(iwdetup(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
          endif
          determ=determ+term
          do 115 i=1,nup
            iwdet=iwdetup(idet)
            ddet_det(1,i)=ddet_det(1,i)+ddeti_deti(1,i,iwdet)*term
            ddet_det(2,i)=ddet_det(2,i)+ddeti_deti(2,i,iwdet)*term
            ddet_det(3,i)=ddet_det(3,i)+ddeti_deti(3,i,iwdet)*term
            ekinen(i)=ekinen(i)+d2edeti_deti(i,iwdet)*term
  115       d2lndet=d2lndet+d2edeti_deti(i,iwdet)*term
          do 117 i=nup+1,nelec
            iwdet=iwdetdn(idet)
            ddet_det(1,i)=ddet_det(1,i)+ddeti_deti(1,i,iwdet)*term
            ddet_det(2,i)=ddet_det(2,i)+ddeti_deti(2,i,iwdet)*term
            ddet_det(3,i)=ddet_det(3,i)+ddeti_deti(3,i,iwdet)*term
            ekinen(i)=ekinen(i)+d2edeti_deti(i,iwdet)*term
  117       d2lndet=d2lndet+d2edeti_deti(i,iwdet)*term

c Inverse of full sum of determinants
      detinv=one/determ

c multiply through to set up logarithmic first and second derivatives
      d2lndet=d2lndet*detinv
      do 120 i=1,nelec
        div_vd(i)=ekinen(i)*detinv
        ekinen(i)=-half*ekinen(i)*detinv
        do 120 k=1,ndim
          ddet_det(k,i)=ddet_det(k,i)*detinv
          div_vd(i)=div_vd(i)-ddet_det(k,i)**2
  120     d2lndet=d2lndet-ddet_det(k,i)**2

c Derivatives wrt to csf_coefs for optimizing them
c Note that the arrays that are needed for vmc and dmc are over ndet but
c those that are needed for optimization only are over nparmcsf.
      if(index(mode,'fit').ne.0 .or. igradhess.gt.0 .or. l_opt_csf) then
        d2det_det=0
        do 125 i=1,nelec
  125     d2det_det=d2det_det-2*ekinen(i)
        do 150 iparm=1,nparmcsf
          icsf=iwcsf(iparm)
          d2deti_det(iparm)=0
          deti_det(iparm)=0
          do 130 i=1,nelec
            do 130 k=1,ndim
  130         ddeti_det(k,i,iparm)=0
          do 150 idet_in_csf=1,ndet_in_csf(icsf)
            idet=iwdet_in_csf(idet_in_csf,icsf)
            if(ndn.ge.1) then
              term=detu(iwdetup(idet))*detd(iwdetdn(idet))*cdet_in_csf(idet_in_csf,icsf)*detinv
             else
              term=detu(iwdetup(idet))*cdet_in_csf(idet_in_csf,icsf)*detinv
            endif
            deti_det(iparm)=deti_det(iparm)+term
            do 140 i=1,nup
              iwdet=iwdetup(idet)
              d2deti_det(iparm)=d2deti_det(iparm)+d2edeti_deti(i,iwdet)*term
              do 140 k=1,ndim
  140           ddeti_det(k,i,iparm)=ddeti_det(k,i,iparm)+ddeti_deti(k,i,iwdet)*term
            do 150 i=nup+1,nelec
              iwdet=iwdetdn(idet)
              d2deti_det(iparm)=d2deti_det(iparm)+d2edeti_deti(i,iwdet)*term
              do 150 k=1,ndim
  150           ddeti_det(k,i,iparm)=ddeti_det(k,i,iparm)+ddeti_deti(k,i,iwdet)*term
        if(ipr.ge.4) write(6,'(''deti_det(iparm) in determinant'',40d12.4)') (deti_det(iparm),iparm=1,nparmcsf)

c Derivatives with respect to orbital parameters (not orbital coefficients!).
        if(iopt.eq.2) then
          do iparm=1,nparmcsf+nparmot
            do jparm=1,nparmcsf+nparmot
              detij_det(iparm,jparm)=0
            enddo
          enddo
        endif

! JT beg : deti_det for first csf (needed for unitary parametrization)
!          icsf=1
!          det1_det=0
!          do idet_in_csf=1,ndet_in_csf(icsf)
!            idet=iwdet_in_csf(idet_in_csf,icsf)
!            term=detu(idet)*detd(idet)*cdet_in_csf(idet_in_csf,icsf)*detinv
!            det1_det=det1_det+term
!          enddo
!          call object_modified_by_index (det1_det_index)
! JT end

        if(nparmot.gt.0) call deriv_det_orb(orb,dorb,ddorb,dporb,d2porb,ddporb,d2dporb,detinv)
      endif

      return
      end

c--------------------------------------------------------------------------------------
      subroutine deriv_det_orb(orb,dorb,ddorb,dporb,d2porb,ddporb,d2dporb,detinv)
c Written by A.D.Guclu, Apr 2006
c Calculates derivatives wrt orbital parameters for optimization

c explanations on local variables:

c detui(iparm,idet)  = derivative wrt the parameter iparm in up determ. idet
c detuij(jparm,idet) = derivative wrt the parameter jparm, of the current detui, in up determ. idet
c                       in other words, second derivatives wrt iparm&jparm, but iparm is not stored.
c ddetui(1:3,ie,idet)= velocity of the current detui, iparm is not stored.
c d2detui(idet)      = laplacian of the current detui, iparm is not stored.

      use dorb_mod
      use orbitals_mod, only: orb_tot_nb
      use dets_mod
      use slater_mod
      use optim_mod
      use const_mod
      use dim_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'force.h'
!JT      include '../fit/fit.h'

c commons
!JT      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl_opt/ nparm,nsig,ncalls,iopt,ipr_opt
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
!JT      common /dorb/ iworbd(MELEC,MDET),iworbdup(MELECUD,MDETUD),iworbddn(MELECUD,MDETUD)
!JT     &,iwdetup(MDET),iwdetdn(MDET),ndetup,ndetdn
!JT      common /slater/ slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD)
!JT     &,fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
!JT     &,fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
!JT     &,detu(MDETUD),detd(MDETUD)
!JT     &,ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD),deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD),d2deti_det(MPARMD),d2det_det
!JT     &,detij_det(MPARMD,MPARMD)
!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase
      common /optimo/ iwo(MORB,MOTYPE),nparmo(MOTYPE),nparmot,notype

c arguments:
      dimension orb(nelec,orb_tot_nb),dorb(3,nelec,orb_tot_nb),ddorb(nelec,orb_tot_nb)
      dimension dporb(MOTYPE,MELEC,MORB),d2porb(MOTYPE,MOTYPE,MELEC,MORB)
      dimension ddporb(3,MOTYPE,MELEC,MORB),d2dporb(MOTYPE,MELEC,MORB)

c local arrays:
      dimension detui(MPARMD,MDETUD),detdi(MPARMD,MDETUD)
      dimension detuij(MPARMD,MDETUD),detdij(MPARMD,MDETUD)
      dimension ddetui(3,MELEC,MDETUD),ddetdi(3,MELEC,MDETUD)
      dimension d2detui(MDETUD),d2detdi(MDETUD)
      dimension zvec(MELEC),anewi(MELEC,MELEC)

      logical found

c initializations of local and output variables:
c complete initialization of detij_det is done in determinant()
      do iparm0=1,nparmot
        iparm=iparm0+nparmcsf
        d2deti_det(iparm)=0
        deti_det(iparm)=0
        do idet=1,ndetup
          detui(iparm0,idet)=0
        enddo
        do idet=1,ndetdn
          detdi(iparm0,idet)=0
        enddo
        do i=1,nelec
          do k=1,ndim
            ddeti_det(k,i,iparm)=0
          enddo
        enddo
      enddo

c      write(*,*) 'detu(idet)',detu(1)
      iparm0=0
      do it=1,notype
        do ip=1,nparmo(it)
          iparm0=iparm0+1
c          write(*,*) 'it,ip,iparm0=',it,ip,iparm0

c more initializations
          do idet=1,ndetup
            do i=1,nelec
              do k=1,ndim
                ddetui(k,i,idet)=0
              enddo
            enddo
            d2detui(idet)=0
            do jparm=1,nparmot
              detuij(jparm,idet)=0
            enddo
          enddo
          do idet=1,ndetdn
            do i=1,nelec
              do k=1,ndim
                ddetdi(k,i,idet)=0
              enddo
            enddo
            d2detdi(idet)=0
            do jparm=1,nparmot
              detdij(jparm,idet)=0
            enddo
          enddo

c start the big loop on determinants
c look in the det if it contains the orbitals (is there a better way?)::
c UP ELECTRONS:
          do idet=1,ndetup
c            if(nup.gt.0) then
            found=.false.
            jopt=0
            do while(.not.found .and. jopt.lt.nup)
              jopt=jopt+1
              if(iworbdup(jopt,idet).eq.iwo(ip,it)) found=.true.
            enddo
c jopt is the orbital to optimize.
c calculate the inverse of the parameter derivative matrix using sherman-morrison:
c (using formulea 2.7.5 in numerical recipes)
c first calculate vector z
            if(found) then
              do iz=1,nup
                zvec(iz)=0
                do ie=1,nup
                  ike=(ie-1)*nup
                  zvec(iz)=zvec(iz)+dporb(it,ie,iworbdup(jopt,idet))*slmui(iz+ike,idet)
                enddo
              enddo
              beta=zvec(jopt)   !    beta=lambda + 1
              zvec(jopt)=zvec(jopt)-1.d0
              detui(iparm0,idet)=beta
              betai=1/beta
c             write(*,*) 'beta=',beta
c get the new inverse matrix:
              do i=1,nup
                jk=-nup
                do j=1,nup
                  jk=jk+nup
                  anewi(i,j)=slmui(i+jk,idet)-zvec(i)*slmui(jopt+jk,idet)*betai
                enddo
              enddo

c use sherman-morrison second time (but this time applied to determinants only)
c to get the coord. derivatives:
              do ie=1,nup
                do j=1,nup
                  if(j.eq.jopt) then
                    do idim=1,ndim
                      ddetui(idim,ie,idet)=ddetui(idim,ie,idet)
     &                     +ddporb(idim,it,ie,iworbdup(j,idet))*anewi(j,ie)*detui(iparm0,idet)
                    enddo
                    d2detui(idet)=d2detui(idet)
     &                   +d2dporb(it,ie,iworbdup(j,idet))*anewi(j,ie)*detui(iparm0,idet)
                  else
                    do idim=1,ndim
                      ddetui(idim,ie,idet)=ddetui(idim,ie,idet)
     &                     +dorb(idim,ie,iworbdup(j,idet))*anewi(j,ie)*detui(iparm0,idet)
                    enddo
                    d2detui(idet)=d2detui(idet)
     &                   +ddorb(ie,iworbdup(j,idet))*anewi(j,ie)*detui(iparm0,idet)
                  endif

                enddo
              enddo

c now get the second derivatives wrt optimization parameters
c if we are doing newton optimization
              if(iopt.eq.2) then
                jparm0=0
                do jt=1,notype
                  do jp=1,nparmo(jt)
                    jparm0=jparm0+1
                    if(jparm0.le.iparm0) then ! due to symmetry
                      if(jp.eq.ip) then ! parameters in same orbitals
                        do ie=1,nup
                          ike=(ie-1)*nup
                          detuij(jparm0,idet)=detuij(jparm0,idet)+
     &                       slmui(jopt+ike,idet)*d2porb(it,jt,ie,iworbdup(jopt,idet))
                        enddo
                      else      ! parameters in different orbitals
                        found=.false. ! look out for the orbital
                        jopt2=0
                        do while(.not.found .and. jopt2.lt.nup)
                          jopt2=jopt2+1
                          if(iworbdup(jopt2,idet).eq.iwo(jp,jt)) found=.true.
                        enddo
                        if(found) then
                          do ie=1,nup
                            detuij(jparm0,idet)=detuij(jparm0,idet)+
     &                           detui(iparm0,idet)*anewi(jopt2,ie)*
     &                           dporb(jt,ie,iworbdup(jopt2,idet))
                          enddo
                        endif
                      endif
                    endif
                  enddo         ! jp
                enddo           ! jt
              endif             ! end second derivatives

            endif               ! found
          enddo                 ! end up-electrons

c DOWN ELECTRONS
          do idet=1,ndetdn
            found=.false.
            jopt=0
            do while(.not.found .and. jopt.lt.ndn)
              jopt=jopt+1
              if(iworbddn(jopt,idet).eq.iwo(ip,it)) found=.true.
            enddo
c jopt is jopt'th down orbital to optimize.
c calculate the inverse of the parameter derivative matrix using sherman-morrison:
c (using formulea 2.7.5 in numerical recipes)
c first calculate vector z
            if(found) then
              do iz=1,ndn
                zvec(iz)=0
                do ie=1,ndn
                  ike=(ie-1)*ndn
                  zvec(iz)=zvec(iz)
     &                   +dporb(it,ie+nup,iworbddn(jopt,idet))*slmdi(iz+ike,idet)
                enddo
              enddo
              beta=zvec(jopt)
              zvec(jopt)=zvec(jopt)-1.d0
              detdi(iparm0,idet)=beta
              betai=1/beta
c                write(*,*) 'beta=',beta
c get the new inverse matrix:
              do i=1,ndn
                jk=-ndn
                do j=1,ndn
                  jk=jk+ndn
                  anewi(i,j)=slmdi(i+jk,idet)-zvec(i)*slmdi(jopt+jk,idet)*betai
                enddo
              enddo
c use sherman-morrison second time (but this time applied to determinants only)
c to get the coord. derivatives:
              do ie=1,ndn
                do j=1,ndn
                  if(j.eq.jopt) then
                    do idim=1,ndim
                      ddetdi(idim,ie,idet)=ddetdi(idim,ie,idet)
     &              +ddporb(idim,it,ie+nup,iworbddn(j,idet))*anewi(j,ie)*detdi(iparm0,idet)
                    enddo
                    d2detdi(idet)=d2detdi(idet)
     &                   +d2dporb(it,ie+nup,iworbddn(j,idet))*anewi(j,ie)*detdi(iparm0,idet)
                  else
                    do idim=1,ndim
                      ddetdi(idim,ie,idet)=ddetdi(idim,ie,idet)
     &              +dorb(idim,ie+nup,iworbddn(j,idet))*anewi(j,ie)*detdi(iparm0,idet)
                    enddo
                    d2detdi(idet)=d2detdi(idet)
     &              +ddorb(ie+nup,iworbddn(j,idet))*anewi(j,ie)*detdi(iparm0,idet)
                  endif
                enddo
              enddo
c              write(*,*) 'd2detdi(iparm0,idet)=',d2detdi(iparm0,idet)
c              write(*,*) 'detdi(iparm0,idet)=',detdi(iparm0,idet)

c now get the second derivatives wrt optimization parameters
c if we are doing newton optimization
              if(iopt.eq.2) then
                jparm0=0
                do jt=1,notype
                  do jp=1,nparmo(jt)
                    jparm0=jparm0+1
                    if(jparm0.le.iparm0) then ! due to symmetry
                      if(jp.eq.ip) then ! parameters in same orbitals
                        do ie=1,ndn
                          ike=(ie-1)*ndn
                          detdij(jparm0,idet)=detdij(jparm0,idet)+
     &                         slmdi(jopt+ike,idet)*d2porb(it,jt,ie+nup,iworbddn(jopt,idet))
                        enddo
                      else      ! parameters in different orbitals
                        found=.false. ! look out for the orbital
                        jopt2=0
                        do while(.not.found .and. jopt2.lt.ndn)
                          jopt2=jopt2+1
                          if(iworbddn(jopt2,idet).eq.iwo(jp,jt)) found=.true.
                        enddo
                        if(found) then
                          do ie=1,ndn
                            detdij(jparm0,idet)=detdij(jparm0,idet)+
     &                           detdi(iparm0,idet)*anewi(jopt2,ie)*
     &                           dporb(jt,ie+nup,iworbddn(jopt2,idet))
                          enddo
                        endif
                      endif
                    endif
                  enddo         ! jp
                enddo           ! jt
              endif             ! end second derivatives



            endif               ! found
          enddo                 ! down electrons


c now put determinants together to get the final results:
          iparm=iparm0+nparmcsf
          do icsf=1,ncsf
            do idet_in_csf=1,ndet_in_csf(icsf)

              idet=iwdet_in_csf(idet_in_csf,icsf)
              if(ndn.ge.1) then
                term=detu(iwdetup(idet))*detd(iwdetdn(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
              else
                term=detu(iwdetup(idet))*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
              endif
c              term=detu(idet)*detd(idet)*csf_coef(icsf,iwf)*cdet_in_csf(idet_in_csf,icsf)
              term=term*detinv
              deti_det(iparm)=deti_det(iparm)+
     &                 (detui(iparm0,iwdetup(idet))+detdi(iparm0,iwdetdn(idet)))*term
              d2deti_det(iparm)=d2deti_det(iparm)+
     &                 (d2detui(iwdetup(idet))+d2detdi(iwdetup(idet)))*term
              do i=1,nup
                iwdet=iwdetup(idet)
                d2deti_det(iparm)=d2deti_det(iparm)+
     &                 d2edeti_deti(i,iwdet)*detdi(iparm0,iwdet)*term
                do k=1,ndim
                  ddeti_det(k,i,iparm)=ddeti_det(k,i,iparm)+(ddetui(k,i,iwdet)+
     &            ddeti_deti(k,i,iwdet)*detdi(iparm0,iwdet))*term
                enddo
              enddo
              do i=nup+1,nelec
                iwdet=iwdetdn(idet)
                d2deti_det(iparm)=d2deti_det(iparm)+
     &                 d2edeti_deti(i,idet)*detui(iparm0,iwdet)*term
                do k=1,ndim
                  ddeti_det(k,i,iparm)=ddeti_det(k,i,iparm)+(ddetdi(k,i,iwdet)+
     &            ddeti_deti(k,i,iwdet)*detui(iparm0,iwdet))*term
                enddo
              enddo

              if(iopt.eq.2) then
                do jparm0=1,iparm0
                  jparm=jparm0+nparmcsf
                  if(ndn.ge.1) then
                    detij_det(iparm,jparm)=detij_det(iparm,jparm)
     &                   +term*(detuij(jparm0,iwdetup(idet))+detdij(jparm0,iwdetdn(idet))
     &                   +detui(jparm0,iwdetup(idet))*detdi(iparm0,iwdetdn(idet))
     &                   +detdi(jparm0,iwdetdn(idet))*detui(iparm0,iwdetup(idet)))
                  else
                    detij_det(iparm,jparm)=detij_det(iparm,jparm)
     &                   +term*(detuij(jparm0,iwdetup(idet)))
                  endif
                  detij_det(jparm,iparm)=detij_det(iparm,jparm)
                enddo
              endif
            enddo
          enddo


        enddo                   ! loop on ip
      enddo                     ! loop on it


      return
      end

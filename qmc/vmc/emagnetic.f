      subroutine emagnetic(ltot)
      use control_mod
c Written by A.D.Guclu, Feb 2004.  Modified by Cyrus Umrigar
c Magnetic energies due to excess of angular momentum and spin are calculated.
c We also verify that all the determinants have same angular momentum

      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'force.h'
!JT      include 'numbas.h'


      common /contrl_per/ iperiodic,ibasis
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /dorb/   iworbd(MELEC,MDET)
      common /coefs/  coef(MBASIS,MORB,MWF),nbasis,norb
      common /basis2/ zex2(MRWF,MCTYPE,MWF),n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
     &,icenter_basis(MBASIS),ictype_basis(MBASIS)
     &,nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),iwrwf2(MBASIS)
      common /basis3/ n_fd(MBASIS),m_fd(MBASIS),n_cf(MBASIS),ncfmax
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
      common /dot/    w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring

      write(6,'(''l_bas'',20i4)') (l_bas(ibas),ibas=1,nbasis)

c If complex orbitals are used it calculates total L here.  Otherwise it is inputted.
      if(ibasis.eq.3) then
        do 60 idet=1,ndet
          ltoti=0
          do 50 iel=1,nup+ndn
c Find angular mom. of the first non-zero basis function.
c Warning: we assume without verification that all the basis functions of a
c given orbital have same l, as should be the case for a "restricted" calculation
              nonzero_coef=0
              ibas=0
              do 40 while(nonzero_coef.eq.0)
                ibas=ibas+1
                if(coef(ibas,iworbd(iel,idet),1).ne.0.d0) then
                  nonzero_coef=1
                  if(numr.eq.0) then
c Fock-Darwin basis
                    ltoti=ltoti+m_fd(ibas)
                   else
c Radial fn. times complex spherical harmonic
                    ltoti=ltoti+l_bas(ibas)
                  endif
                endif
                if(ibas.gt.nbasis) stop 'all the coefficients are zero in emagnetic.f'
   40         enddo
   50       enddo
c       save the total angular momentum of the first determinant and
c       keep evaluating ltot for remaining determinants for debuging:
            if(idet.eq.1) then
              ltot=ltoti
             elseif(ltot.ne.ltoti) then
              stop 'determinants have different total angular momentum'
            endif
   60   enddo
      endif

c calculate magnetic energy due to angular momentum:
      emaglz=-0.5d0*bext*ltot
c calculate magnetic energy due to spin (zeeman term):
      emagsz=-0.25d0*glande*bext*(nup-ndn)
      write(6,*)
      write(6,'(''determinantal angular momentum, ltot ='',t31,i10)') ltot
      write(6,'(''angular mom. magnetic energy ='',t31,f10.5)') emaglz
      write(6,'(''spin magnetic energy ='',t31,f10.5)') emagsz

      return
      end

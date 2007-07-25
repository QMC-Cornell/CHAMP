c     subroutine quad_min(energy,add_diag,force,force_err,npts,eig_min,eig_max)
c     subroutine quad_min(energy_sav,energy_err_sav,energy,energy_err,force,force_err,ene_var,add_diag,npts,eig_min,eig_max,p_var)
c     subroutine quad_min(energy_sav,energy_err_sav,ene_var,npts)
      subroutine quad_min(ene_var,npts)
c Find the value of add_diag at which energy is minimum by fitting to a 3-parameter function.
c At present we are fitting a parabola in log10(add_diag) and that works well.
c In some toy tests a modification of this worked even better.
c Warning: at present we are fitting to the energy, even when we are optimizing a linear
c combination of energy and variance.  This results in add_diag gettting large if we
c start from a energy-minimized psi and then optimize the variance.  So, I should fit
c to the right linear combination.
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include '../fit/fit.h'

      parameter(MFUNC=3)
      common /contrl_opt/ nparm,nsig,ncalls,iopt,ipr_opt
      common /gradhess/ grad(MPARM),grad_var(MPARM),hess(MPARM,MPARM),hess_var(MPARM,MPARM),gerr(MPARM),
     &add_diag(3),energy(3),energy_sigma(3),energy_err(3),force(3),force_err(3),
     &eig_min,eig_max,p_var,tol_energy,nopt_iter,nblk_max

c     dimension energy(npts),energy_err(npts),force(npts),force_err(npts),ene_var(npts),
c    &add_diag(npts),add_diag_log(npts),a(MFUNC,MFUNC),b(MFUNC)
      dimension ene_var(npts),add_diag_log(MFUNC),a(MFUNC,MFUNC),b(MFUNC)

      nfunc=MFUNC

      if(eig_max.lt.1.d99) then
        dlog10_eig_max=dlog10(eig_max)
       else
        write(6,'(''dlog10_eig_max reset to 1.d300'')')
        dlog10_eig_max=1.d300
      endif

      do 5 k=1,npts
    5   add_diag_log(k)=dlog10(add_diag(k))
ctmp    if(add_diag(k).lt.1.d99) then
ctmp      tmp=dlog10(add_diag(k)+max(eig_min,0.d0))
ctmp      add_diag_log(k)=tmp/(1+tmp/dlog10_eig_max)
ctmp     else
ctmp      add_diag_log(k)=dlog10(eig_max)
ctmp    endif
ctmp5 continue

      do 30 i=1,nfunc
        b(i)=0
        do 10 k=1,npts
   10     b(i)=b(i)+ene_var(k)*add_diag_log(k)**(i-1)
        do 30 j=1,i
          a(i,j)=0
          do 20 k=1,npts
   20       a(i,j)=a(i,j)+add_diag_log(k)**(i+j-2)
   30     a(j,i)=a(i,j)

c     write(6,*) b
c     write(6,*)
c     write(6,'(3f15.6)') a

c Do cholesky decomposition
      call chlsky(a,nfunc,MFUNC,ierr)
      if(ierr.ne.0) stop 'ierr ne 0 in chlsky'

c Symmetrize decomposed matrix (needs to be done before calling uxb
c or need to modify uxb)
      do 40 i=1,nfunc
        do 40 j=i+1,nfunc
   40     a(i,j)=a(j,i)

c Solve linear equations
      call lxb(a,nfunc,MFUNC,b)
      call uxb(a,nfunc,MFUNC,b)

      write(6,*)
      write(6,'(a)') 'Find optimal add_diag by quadratic fit of 3 points:'
      write(6,'(a,3f10.5)') 'add_diag        = ', (add_diag(k),k=1,npts)
      write(6,'(a,3f10.5)') 'log(add_diag)   = ', (add_diag_log(k),k=1,npts)
      write(6,'(a,3f10.5)') 'energy_variance = ', (ene_var(k),k=1,npts)
      write(6,'(a,3f10.5)') 'energy          = ', (energy(k),k=1,npts)
      if(ipr_opt.ge.1) write(6,'(''parameters for fit of energy/variance, b='',3d12.4)') b

      ene_var_min=9.d99
      ene_var_max=-9.d99
c     energy_min=9.d99
c     energy_max=-9.d99
      rms=0
      k_min=0
      do 50 k=1,npts
        ee=b(1)+b(2)*add_diag_log(k)+b(3)*add_diag_log(k)**2
        if(ene_var(k).lt.ene_var_min) then
          k_min=k
          ene_var_min=ene_var(k)
        endif
        ene_var_max=max(ene_var_max,ene_var(k))
c       energy_min=min(energy_min,energy(k))
c       energy_max=max(energy_max,energy(k))
        rms=rms+(ee-ene_var(k))**2
   50   if(ipr_opt.ge.3) write(6,'(9g14.6)') add_diag(k),ee,ene_var(k),ee-ene_var(k)
      rms=dsqrt(rms/npts)
      if(ipr_opt.ge.2) write(6,'(''rms error in fit of energy or variance to get optimal add_diag is'',d12.4)') rms
      if(k_min.eq.0) then
        write(6,'(''None of the energies in quad_min are reasonable'')')
        stop 'None of the energies in quad_min are reasonable'
      endif

c roughly determine if difference between energies is statistically significant
c     if((energy_max-energy_min)**2.gt.energy_err(1)**2+energy_err(2)**2+energy_err(3)**2) then
c       ienergy_diff_signif=1
c      else
c       ienergy_diff_signif=0
c     endif

c We are assuming in the foll. that the three add_diag values are in the ratio 1,0.1,10
c or at least that the 2nd is the smallest and the 3rd is the largest.
c 1) if p_var<0.1 and force_err is small enough that force is significant, then make quadratic approximation
c    but limit it to changing add_diag(1) by no more than a factor of 100 in either direction
cx2) if energy(2) is not more than a std. dev. higher than in previous iteration use add_diag(2)
cx3) if energy(1) is not more than a std. dev. higher than in previous iteration use add_diag(1)
cx2) if energy(2) is not more than a std. dev. higher than energy(1) and the ene_var are close use add_diag(2)
cx3) if energy(1) is not more than a std. dev. higher than energy(3) and the ene_var are close use add_diag(1)
cx4) otherwise use add_diag corresponding to whichever of the 3 ene_var is lowest.
c 2) if p_var*(energy(2)+3*energy_err(2)-force_err(2))+q_var*energy_sigma(2)**2 <
c       p_var*(energy(1)+3*energy_err(1)             )+q_var*energy_sigma(1)**2 then use add_diag(2)
c 3) if p_var*(energy(1)+3*energy_err(1)             )+q_var*energy_sigma(1)**2 <
c       p_var*(energy(3)+3*energy_err(3)+force_err(3))+q_var*energy_sigma(3)**2 then use add_diag(1)
c    By having -force_err(2) in 2 and +force_err(3) in 3 we are favoring lower add_diag values when there
c    is no statistically significant difference.
c 3) otherwise use add_diag(3)
cx3) otherwise use add_diag corresponding to whichever of the 3 ene_var is lowest.
cx5) if energy(2) < energy(1) and add_diag is 1 or larger, then reduce add_diag by another factor of 10
c The reason for 4) and 5) is mostly to prevent add_diag_min from getting large when optimizing
c scalek with variance minimization.
c We do not have the error bar on the variance so use the first part of foll.
c if statement only if we are doing mostly energy minimization
      if(p_var.lt.0.1d0 .and. b(3).gt.0 .and. abs(force(2)).gt.3*force_err(2) .and. abs(force(3)).gt.3*force_err(3)) then
        iwadd_diag=0
        add_diag_log_min=-0.5d0*b(2)/b(3)
        add_diag_log_min=min(max(add_diag_log_min,add_diag_log(1)-2.d0),add_diag_log(1)+2.d0)
cc     elseif(energy(2).lt.energy_sav+sqrt(energy_err_sav**2+energy_err(2)**2)) then
cc      add_diag_log_min=add_diag_log(2)
cc     elseif(energy(1).lt.energy_sav+sqrt(energy_err_sav**2+energy_err(1)**2)) then
cc      add_diag_log_min=add_diag_log(1)
c      elseif(energy(2).lt.energy(1)+force_err(2) .and. ene_var_max-ene_var_min .lt. 1.d-3*abs(ene_var_max)) then
c       iwadd_diag=2
c       add_diag_log_min=add_diag_log(2)
c       if(energy(2).lt.energy(1) .and. add_diag_log_min.ge.0.d0) add_diag_log_min=add_diag_log_min-1.d0
c      elseif(energy(1).lt.energy(3)+force_err(3) .and. ene_var_max-ene_var_min .lt. 1.d-2*abs(ene_var_max) .and. k_min.eq.3) then
c       iwadd_diag=1
c       add_diag_log_min=add_diag_log(1)
c      else
c       iwadd_diag=k_min
c       add_diag_log_min=add_diag_log(k_min)
c       if(energy(2).lt.energy(1) .and. add_diag_log_min.ge.0.d0) add_diag_log_min=add_diag_log_min-1.d0
       elseif(ene_var(2)+(1-p_var)*(3*energy_err(2)-force_err(2)).lt.ene_var(1)+(1-p_var)*(3*energy_err(1))) then
        iwadd_diag=2
        add_diag_log_min=add_diag_log(2)
       elseif(ene_var(1)+(1-p_var)*3*energy_err(1).lt.ene_var(3)+(1-p_var)*(3*energy_err(3)+force_err(3))) then
        iwadd_diag=1
        add_diag_log_min=add_diag_log(1)
       else
        iwadd_diag=3
        add_diag_log_min=add_diag_log(3)
      endif
c     if(energy(2).lt.energy(1) .and. add_diag_log_min.ge.0.d0) then
c       add_diag_log_min=add_diag_log_min-1.d0
c       write(6,'(''aaaa'',9f12.5)') energy(2),energy(1),add_diag_log_min
c     endif
      write(6,'(''add_diag_min before bounding is'',1p,d12.4,0p)') 10**add_diag_log_min
c     add_diag_log_min=min(max(add_diag_log_min,-8.d0),6.d0)
      add_diag_log_min=max(add_diag_log_min,-8.d0)
      add_diag_min=10**add_diag_log_min
ctmp  add_diag_min=10**(add_diag_log_min/(1-add_diag_log_min/dlog10_eig_max))-max(eig_min,0.d0)
      ene_var_min=b(1)+b(2)*add_diag_log_min+b(3)*add_diag_log_min**2
      write(6,'(''iwadd_diag,add_diag_log_min,add_diag_min,ene_var_min='',i2,f7.3,1p,d9.2,0p,f12.6)')
     &iwadd_diag,add_diag_log_min,add_diag_min,ene_var_min

      add_diag(1)=add_diag_min

      return
      end

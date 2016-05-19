      subroutine quad_min(ene_var,npts)
! Find the value of add_diag at which energy is minimum by fitting to a 3-parameter function.
! At present we are fitting a parabola in log10(add_diag) and that works well.
! If the parabola has a maximum, use the best of the 3 sampled add_diag values.
! If the fit is bad (because one of the values is way different), make a conservative choice.
! Written by Cyrus Umrigar

      use gradhess_mod
      use contrl_opt_mod
      implicit real*8(a-h,o-z)
      parameter(MFUNC=3)

      dimension ene_var(npts),add_diag_log(MFUNC),a(MFUNC,MFUNC),b(MFUNC)

      nfunc=MFUNC

      do 5 k=1,npts
    5   add_diag_log(k)=dlog10(add_diag(k))

      do 30 i=1,nfunc
        b(i)=0
        do 10 k=1,npts
   10     b(i)=b(i)+ene_var(k)*add_diag_log(k)**(i-1)
        do 30 j=1,i
          a(i,j)=0
          do 20 k=1,npts
   20       a(i,j)=a(i,j)+add_diag_log(k)**(i+j-2)
   30     a(j,i)=a(i,j)

! Do cholesky decomposition
      call chlsky(a,nfunc,MFUNC,ierr)
      if(ierr.ne.0) stop 'ierr ne 0 in chlsky'

! Symmetrize decomposed matrix (needs to be done before calling uxb or need to modify uxb)
      do 40 i=1,nfunc
        do 40 j=i+1,nfunc
   40     a(i,j)=a(j,i)

! Solve linear equations
      call lxb(a,nfunc,MFUNC,b)
      call uxb(a,nfunc,MFUNC,b)

      write(6,*)
      write(6,'(a)') 'Find optimal add_diag by quadratic fit of 3 points:'
      write(6,'(a,3es10.3)') 'add_diag        = ', (add_diag(k),k=1,npts)
      write(6,'(a,3f10.5)') 'log(add_diag)   = ', (add_diag_log(k),k=1,npts)
      write(6,'(a,3f10.5)') 'energy_variance = ', (ene_var(k),k=1,npts)
      write(6,'(a,3f10.5)') 'energy          = ', (energy(k),k=1,npts)
      if(ipr_opt.ge.1) write(6,'(''parameters for fit of energy/variance, b='',3es12.4)') b

      ene_var_min=9.d99
      ene_var_max=-9.d99
      rms=0
      k_min=0
      do 50 k=1,npts
        ee=b(1)+b(2)*add_diag_log(k)+b(3)*add_diag_log(k)**2
        if(ene_var(k).lt.ene_var_min) then
          k_min=k
          ene_var_min=ene_var(k)
        endif
        ene_var_max=max(ene_var_max,ene_var(k))
        rms=rms+(ee-ene_var(k))**2
   50   if(ipr_opt.ge.3) write(6,'(9g14.6)') add_diag(k),ee,ene_var(k),ee-ene_var(k)
      rms=dsqrt(rms/npts)
      write(6,'(''rms error in fit of energy or variance to get optimal add_diag is'',es9.2)') rms
      if(k_min.eq.0) then
        write(6,'(''None of the energies in quad_min are reasonable'')')
        stop 'None of the energies in quad_min are reasonable'
      endif

! We are assuming in the foll. that the three add_diag values are in the ratio 1,0.1,10
! or at least that the 2nd is the smallest and the 3rd is the largest.
!x1) if b(3)> 0 and p_var<0.1 and force_err is small enough that force is significant,
! 1) if b(3)> 0
!    make quadratic approximation but limit it to changing add_diag(1) by no more than a factor of 100 in either direction
! 2) otherwise use add_diag corresponding to whichever of the 3 ene_var+(1-p_var)*3*energy_err is lowest.
! We do not have the error bar on the variance so add only (1-p_var)*3*energy_err
!x2) if energy(2) is not more than a std. dev. higher than in previous iteration use add_diag(2)
!x3) if energy(1) is not more than a std. dev. higher than in previous iteration use add_diag(1)
!x2) if energy(2) is not more than a std. dev. higher than energy(1) and the ene_var are close use add_diag(2)
!x3) if energy(1) is not more than a std. dev. higher than energy(3) and the ene_var are close use add_diag(1)
!x4) otherwise use add_diag corresponding to whichever of the 3 ene_var is lowest.
!x2) if p_var*(energy(2)+3*energy_err(2)-force_err(2))+q_var*energy_sigma(2)**2 <
!x      p_var*(energy(1)+3*energy_err(1)             )+q_var*energy_sigma(1)**2 then use add_diag(2)
!x3) if p_var*(energy(1)+3*energy_err(1)             )+q_var*energy_sigma(1)**2 <
!x      p_var*(energy(3)+3*energy_err(3)+force_err(3))+q_var*energy_sigma(3)**2 then use add_diag(1)
!x   By having -force_err(2) in 2 and +force_err(3) in 3 we are favoring lower add_diag values when there
!x   is no statistically significant difference.
!x3) otherwise use add_diag(3)
!x3) otherwise use add_diag corresponding to whichever of the 3 ene_var is lowest.
!x5) if energy(2) < energy(1) and add_diag is 1 or larger, then reduce add_diag by another factor of 10
!x The reason for 4) and 5) is mostly to prevent add_diag_min from getting large when optimizing
!x scalek with variance minimization.

!     if(p_var.lt.0.1d0 .and. b(3).gt.0 .and. abs(force(2)).gt.3*force_err(2) .and. abs(force(3)).gt.3*force_err(3)) then
      if(b(3).gt.0 .and. rms.lt.1.d-6) then ! quadratic fit has a minimum
        iwadd_diag=0
        add_diag_log_min=-0.5d0*b(2)/b(3)
        add_diag_log_min=min(max(add_diag_log_min,add_diag_log(1)-2.d0),add_diag_log(1)+2.d0)
        ene_var_min=b(1)+b(2)*add_diag_log_min+b(3)*add_diag_log_min**2
      elseif(rms.lt.1.d-6) then ! if quadratic fit is good but has maximum, choose the best of the 3 sampled add_diag values
        iwadd_diag=minloc(ene_var+(1-p_var)*3*energy_err,1)
        add_diag_log_min=add_diag_log(iwadd_diag)
        ene_var_min=ene_var(iwadd_diag)+(1-p_var)*3*energy_err(iwadd_diag)
      else ! rms is big usually when the value for (2) is very high and so quadratic would give optimal about halfway between (1) and (3), whereas it really should be at (3) or higher.  So, set it 100 times add_diag(1)
        add_diag_log_min=add_diag_log(1)+2
        ene_var_min=1.d9
      endif

! Make conservative choice based on energy_err because it can happen that (2) won only because we limited ratio of wts in metrop
! Remember add_diag increases from (2) to (1) to (3), so (3) is most conservative choice.
! If the error of (2) is considerably worse than error of (3), choose (1)
      if(add_diag_log_min.lt.add_diag_log(1) .and. energy_err(2) .gt. 1.1d0*energy_err(1)) then
        iwadd_diag=1
        add_diag_log_min=add_diag_log(1)
        ene_var_min=ene_var(1)+(1-p_var)*3*energy_err(1)
      endif
! If the errors of (1) and (2) are considerably worse than those of (3), choose (3)
      if(add_diag_log_min.lt.add_diag_log(3) .and. energy_err(1) .gt. 1.1d0*energy_err(3) .and. &
       &energy_err(2) .gt. 1.1d0*energy_err(3)) then
        iwadd_diag=3
        add_diag_log_min=add_diag_log(3)
        ene_var_min=ene_var(3)+(1-p_var)*3*energy_err(3)
      endif

      write(6,'(''add_diag_min before bounding is'',es12.4,0p)') 10**add_diag_log_min
      add_diag_log_min=max(add_diag_log_min,-12.d0)
      add_diag_min=10**add_diag_log_min
      write(6,'(''iwadd_diag,add_diag_log_min,add_diag_min,ene_var_min='',i2,f7.3,es10.3,1pg12.4)') &
     &iwadd_diag,add_diag_log_min,add_diag_min,ene_var_min

      add_diag(1)=add_diag_min

      return
      end

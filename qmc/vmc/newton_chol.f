      subroutine optim_options
c Written by Cyrus Umrigar
      use contrl_opt_mod
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include '../fit/fit.h'

      common /linear/ ham(MPARM,MPARM),ovlp(MPARM,MPARM),coef(MPARM,MPARM)
c     common /linear_sav/ ham_sav(MPARM,MPARM),ovlp_sav(MPARM,MPARM),renorm_ovlp(MPARM)
      common /linear_orig/ ovlp_orig(MPARM,MPARM)

      if(mod(iopt,10).eq.1) then
        nparmp1=nparm+1
c initialize superdiagonal elements of ovlp
        do 10 i=1,nparmp1
          do 10 j=1,i-1
   10       ovlp(j,i)=ovlp(i,j)

        do 15 i=1,nparmp1
          do 15 j=1,nparmp1
   15       ovlp_orig(i,j)=ovlp(i,j)

c Symmetrize H
        if(mod(iopt/10,10).eq.1) then
          do 20 i=1,nparmp1
            do 20 j=1,i-1
              ham(i,j)=0.5d0*(ham(i,j)+ham(j,i))
   20         ham(j,i)=ham(i,j)
        endif
c Use basis fns. orthogonal to current ground state
        if(mod(iopt/100,10).eq.1) then
          eave=ham(1,1)/ovlp(1,1)
          do 30 i=2,nparmp1
            do 30 j=2,nparmp1
   30         ham(i,j)=ham(i,j)-ovlp(1,i)*ham(1,j)-ovlp(1,j)*ham(i,1)+eave*ovlp(1,i)*ovlp(1,j)
          do 40 i=2,nparmp1
            ham(1,i)=ham(1,i)-eave*ovlp(1,i)
   40       ham(i,1)=ham(i,1)-eave*ovlp(i,1)
          do 50 i=2,nparmp1
            do 50 j=2,nparmp1
   50         ovlp(i,j)=ovlp(i,j)-ovlp(1,i)*ovlp(1,j)
          do 60 i=2,nparmp1
              ovlp(1,i)=0
   60         ovlp(i,1)=0
        endif
       elseif(mod(iopt,10).eq.2) then
c      elseif(mod(iopt,10).eq.3) then
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ham_ovlp_grad_hess
c Written by Cyrus Umrigar

c Setup and save Hamiltonian and overlap matrices for linear method, and, gradient and Hessian for Newton.
c Normalize them and save the normalization for transforming the parameter variations.
c Save grad_sav,hess_sav,ham_sav,ovlp_sav,renorm_ovlp
c In grad_sav and hess_sav, use appropriate linear combination of energy and variance.
c Evaluate the eigenvalues of the Hessian of the objective function (linear comb of energy and variance).
c In fact we should also renormalize grad_sav and hess_sav.
      use gradhess_mod
      use contrl_opt2_mod
      use contrl_opt_mod
      implicit real*8(a-h,o-z)
      parameter(MBUF=1024,MWORK=MBUF+MPARM*MPARM)
      common /gradhess_sav/ hess_sav(MPARM,MPARM),grad_sav(MPARM)
      common /linear/ ham(MPARM,MPARM),ovlp(MPARM,MPARM),coef(MPARM,MPARM)
      common /linear_sav/ ham_sav(MPARM,MPARM),ovlp_sav(MPARM,MPARM),renorm_ovlp(MPARM)
      dimension work(MWORK),eigv(MPARM)

c Needed for both linear and Newton methods

      nparmp1=nparm+1
      do 10 i=1,nparmp1
   10   renorm_ovlp(i)=sqrt(ovlp(i,i))
      if(ipr_opt.ge.1)  write(6,'(/,''renorm_ovlp='',100d12.4)') (renorm_ovlp(i),i=1,nparmp1)


c Linear method quantities
      if(mod(iopt,10).eq.1) then

      if(ipr_opt.ge.1) then
        write(6,'(/,''Hamiltonian and overlap matrices before renormalization'')')
        do 50 i=1,nparmp1
   50     write(6,'(''ovlp='',20g12.4)') (ovlp(i,j),j=1,i)
        write(6,*)

        do 60 i=1,nparmp1
   60     write(6,'(''ham='',20g12.4)') (ham(i,j),j=1,nparmp1)
        write(6,*)
      endif

      do 80 i=1,nparmp1
        do 80 j=1,i
   80     ovlp(i,j)=ovlp(i,j)/(renorm_ovlp(i)*renorm_ovlp(j))
      do 90 i=1,nparmp1
        do 90 j=1,nparmp1
c         write(6,'(''i,j,renorm_ovlp(i)*renorm_ovlp(j)'',2i3,20g14.4)') i,j,renorm_ovlp(i)*renorm_ovlp(j)
c         write(6,'(''i,j,ham(i,j)'',2i3,20g14.4)') i,j,ham(i,j)
   90     ham(i,j)=ham(i,j)/(renorm_ovlp(i)*renorm_ovlp(j))
c     endif

c set the superdiag elements of ovlp
      do 100 i=1,nparmp1
        do 100 j=1,i-1
  100     ovlp(j,i)=ovlp(i,j)

      do 110 i=1,nparmp1
        do 110 j=1,nparmp1
          ham_sav(i,j)=ham(i,j)
  110     ovlp_sav(i,j)=ovlp(i,j)

      if(ipr_opt.ge.1) then
        write(6,'(/,''Hamiltonian and overlap matrices after renormalization'')')
        do 150 i=1,nparmp1
c 150     write(6,'(''ovlp='',20g12.4)') (ovlp(i,j),j=1,i)
  150     write(6,'(''ovlp='',20g15.8)') (ovlp(i,j),j=1,i)
        write(6,*)

        do 160 i=1,nparmp1
c 160     write(6,'(''ham='',20g12.4)') (ham(i,j),j=1,nparmp1)
  160     write(6,'(''ham='',20g15.8)') (ham(i,j),j=1,nparmp1)
        write(6,*)
      endif
      write(6,'(/,''diagonal H/O'',9f9.4)') (ham(i,i)/ovlp(i,i),i=1,nparmp1)

      endif ! end of linear method quantities


c Newton method quantities
c Save the linear combination of gradient and Hessian we wish to optimize.
c It is not clear if it is best to take linear combination of hessian too or not.
c We could skip the rest if we are not doing Newton, but for the moment we keep it
c in order to see the eigenvalues of the Hessian.

      q_var=1-p_var
      grad_opt_norm=0
      do 210 i=1,nparm
        grad_sav(i)=(p_var*grad_var(i)+q_var*grad(i))/renorm_ovlp(i+1)
        grad_opt_norm=grad_opt_norm+grad_sav(i)**2
        do 210 j=1,i
          hess_sav(i,j)=(p_var*hess_var(i,j)+q_var*hess(i,j))/(renorm_ovlp(i+1)*renorm_ovlp(j+1))
  210     hess_sav(j,i)=hess_sav(i,j)
      grad_opt_norm=sqrt(grad_opt_norm)
      if(ipr_opt.ge.-4) write(6,'(''grad_opt_norm='',9f10.5)') grad_opt_norm

      do 220 i=1,nparm
        grad(i)=grad_sav(i)
        do 220 j=1,nparm
  220     hess(i,j)=hess_sav(i,j)

      if(ipr_opt.ge.0) then
        write(6,'(''grad='',10g14.6)') (grad(i),i=1,nparm)
        do 230 i=1,nparm
  230     write(6,'(''hess='',10g14.6)') (hess(i,j),j=1,i)
      endif
c     stop 'Warning: work on newton to be completed'

c Calculate eigenvalues of Hessian
      lwork=nparm**2+MBUF
      if(nparm.gt.MPARM) stop 'nparm > MPARM'
      call dsyev('V','U',nparm,hess,MPARM,eigv,work,lwork,info)
      if(ipr_opt.ge.-1) write(6,'(''eigs='',(1p90d9.1))') (eigv(i),i=1,nparm)

      eig_min=1.d99
      eig_max=-1.d99
      do 240 i=1,nparm
        eig_max=max(eig_max,eigv(i))
  240   eig_min=min(eig_min,eigv(i))
      if(ipr_opt.ge.-2) write(6,'(''eig_min,eig_max='',(1p20d8.1))') eig_min,eig_max

c Make sure that the add_diag values are not tiny compared to eig_min if we are using Newton method.
c I think it is not needed for the linear method, but it may be worth checking some more.
      if(iadd_diag_opt.eq.2 .and. add_diag(1).lt.1.d-1*eig_min) then
        write(6,'(''add_diag(1) increased to 1.d-1*eig_min='',1pd11.4)') 1.d-1*eig_min
        add_diag(1)=max(add_diag(1),1.d-1*eig_min)
        add_diag(2)=0.1d0*add_diag(1)
        add_diag(3)=10.d0*add_diag(1)
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine new_param(iadd_diag,ipr_eigs,ipr_new,iflag)
c Written by Cyrus Umrigar

c Used to find the change in the wavefunction parameters.
c Given gradient, grad, and Hessian, hess, solve hess*x=grad for x.
c The required change in the parameters is -x.
c If ipr_new = 0 then we print new parameters without _new subscript
c            = 1 then we print new parameters with _new subscript if iflag=0
c            = 2 then we print new parameters with _new subscript
      use contrl_opt_mod
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include '../fit/fit.h'

      if(mod(iopt,10).eq.1) then
        call linear_method(iadd_diag,ipr_eigs)
       elseif(mod(iopt,10).eq.2) then
        call newton_chol(iadd_diag,ipr_eigs)
       elseif(mod(iopt,10).eq.3) then
        call perturbation(iadd_diag,ipr_eigs)
      endif

      call update_params(iadd_diag,ipr_new,iflag)

      return
      end
c-----------------------------------------------------------------------

      subroutine newton_chol(iadd_diag,ipr_eigs)
c Written by Cyrus Umrigar

c Used to find the change in the wavefunction parameters.
c Given gradient, grad, and Hessian, hess, solve hess*x=grad for x.
c The required change in the parameters is -x.
c If ipr_new = 0 then we print new parameters without _new subscript
c            = 1 then we print new parameters with _new subscript if iflag=0
c            = 2 then we print new parameters with _new subscript
      use gradhess_mod
      use contrl_opt_mod
      implicit real*8(a-h,o-z)

      common /optim2/ dparm(MPARM)
      common /gradhess_sav/ hess_sav(MPARM,MPARM),grad_sav(MPARM)
      common /linear_sav/ ham_sav(MPARM,MPARM),ovlp_sav(MPARM,MPARM),renorm_ovlp(MPARM)
      dimension grad_cal(MPARM)

      if(ipr_opt.ge.-4) write(6,'(/,''add_diag('',i1,'')='',1pd11.4)') iadd_diag,add_diag(iadd_diag)

c Add max(-lowest eigenvalue,0) + add_diag to diagonal of hessian
c Another possibility is first diagonalizing hess by a unitary transform
c taking the inverse of the diagonalized matrix to construct the
c inverse of A, but resetting the diagonal elements of the diagonalized A
c to a minimum positive value if they are less than this value.
c We are resetting hess because it was destroyed by call to dsyev or chlsky.
      do 8 i=1,nparm
        grad(i)=grad_sav(i)
        do 8 j=1,nparm
          hess(i,j)=hess_sav(i,j)
          if(i.eq.j) then
c           hess(i,j)=hess(i,j)+2*max(-eig_min,0.d0)+add_diag(iadd_diag)
            hess(i,j)=hess(i,j)+max(-eig_min,0.d0)+add_diag(iadd_diag)
          endif
    8 continue
      if(ipr_eigs.ge.1 .and. ipr_opt.ge.1) then
        do 9 i=1,nparm
    9     write(6,'(''hess='',9g12.4)') (hess(i,j),j=1,i)
      endif

c Foll. lines are just to check change of eigenvals due to adding to diag
c     added=max(-eig_min,0.d0)+add_diag(iadd_diag)
c     eig_min_sav=eig_min
c     eig_max_sav=eig_max
c     call dsyev('V','U',nparm,hess,MPARM,eigv,work,lwork,info)
c     write(6,'(''eigs='',(1p20d9.1))') (eigv(i),i=1,nparm)

c     eig_min=1.d99
c     eig_max=-1.d99
c     do i=1,nparm
c       eig_max=max(eig_max,eigv(i))
c       eig_min=min(eig_min,eigv(i))
c     enddo
c     write(6,'(''eig_min_new,eig_max_new='',(1p20g8.1))') eig_min,eig_max
c     write(6,'(''d_eig_min,added='',(1p20g12.5))') eig_min-eig_min_sav,added,eig_min-eig_min_sav-added
c     write(6,'(''d_eig_max,added='',(1p20g12.5))') eig_max-eig_max_sav,added,eig_max-eig_max_sav-added

c     do i=1,nparm
c       do j=1,i
c         hess(i,j)=hess_sav(i,j)
c         hess(j,i)=hess_sav(j,i)
c       enddo
c     enddo


c Do cholesky decomposition
      call chlsky(hess,nparm,MPARM,ierr)
      if(ierr.ne.0) goto 99

c Symmetrize decomposed matrix (needs to be done before calling uxb
c or need to modify uxb)
      do 20 i=1,nparm
        do 20 j=i+1,nparm
   20     hess(i,j)=hess(j,i)

c Solve linear equations
      call lxb(hess,nparm,MPARM,grad)
      call uxb(hess,nparm,MPARM,grad)

c calculate rms change and test solution by checking that the Hessian multiplying the
c solution (which overwrites the gradient in array grad) gives the gradient
c To test solution first restore hess to value before calling chlsky

      do 28 i=1,nparm
        do 28 j=1,nparm
          hess(i,j)=hess_sav(i,j)
          if(i.eq.j) then
c           hess(i,j)=hess_sav(i,j)+2*max(-eig_min,0.d0)+add_diag(iadd_diag)
            hess(i,j)=hess_sav(i,j)+max(-eig_min,0.d0)+add_diag(iadd_diag)
          endif
   28     continue

      dparm_norm=0
      do 30 i=1,nparm
        dparm(i)=-grad(i)
        dparm_norm=dparm_norm+dparm(i)**2
        grad_cal(i)=0
        do 30 j=1,nparm
   30     grad_cal(i)=grad_cal(i)+hess(i,j)*grad(j)
      dparm_norm=sqrt(dparm_norm/nparm)
c     if(ipr_eigs.ge.1 .and. ipr_opt.ge.-4) write(6,'(''newton_chol: iadd_diag,dparm_norm='',i2,9f10.5)')
      if(ipr_opt.ge.-4) write(6,'(''newton_chol: iadd_diag,dparm_norm='',i2,9f10.5)')
     &iadd_diag,dparm_norm
      if(ipr_eigs.ge.1 .and. ipr_opt.ge.1) write(6,'(''grad_cal='',9g12.4)') (grad_cal(i),i=1,nparm)
c     if(ipr_eigs.ge.1 .and. ipr_opt.ge.0) write(6,'(''dparm='',100f10.6)') (dparm(i),i=1,nparm)

      write(6,'(''dparm='',100f10.6)') (dparm(i),i=1,nparm)

      do 40 i=1,nparm
   40  dparm(i)=dparm(i)/renorm_ovlp(i+1)

      write(6,'(''dparm='',100f10.6)') (dparm(i),i=1,nparm)

c Restore grad and hess to original value
      do 95 i=1,nparm
        grad(i)=grad_sav(i)
        do 95 j=1,nparm
   95     hess(i,j)=hess_sav(i,j)

      return
   99 stop 'ierr.ne.0 in chlsky'
      end
c-----------------------------------------------------------------------

      subroutine perturbation(iadd_diag,ipr_eigs)
c Written by Cyrus Umrigar

c Used to find the change in the wavefunction parameters.
c Given gradient, grad, and Hessian, hess, solve hess*x=grad for x.
c The required change in the parameters is -x.
c If ipr_new = 0 then we print new parameters without _new subscript
c            = 1 then we print new parameters with _new subscript if iflag=0
c            = 2 then we print new parameters with _new subscript
      use gradhess_mod
      use contrl_opt_mod
      implicit real*8(a-h,o-z)

      common /optim2/ dparm(MPARM)
c     common /gradhess_sav/ hess_sav(MPARM,MPARM),grad_sav(MPARM)
      common /linear/ ham(MPARM,MPARM),ovlp(MPARM,MPARM),coef(MPARM,MPARM)
      common /linear_sav/ ham_sav(MPARM,MPARM),ovlp_sav(MPARM,MPARM),renorm_ovlp(MPARM)
      dimension e_diag(MPARM)
c     dimension aa(MPARM,MPARM)

      write(6,'(''using perturbation theory'')')
      if(ipr_opt.ge.-4) write(6,'(/,''add_diag('',i1,'')='',1pd11.4)') iadd_diag,add_diag(iadd_diag)

      nparmp1=nparm+1

c initialize the superdiag elements of ovlp
      do 2 i=1,nparmp1
        do 2 j=1,i-1
    2     ovlp(j,i)=ovlp(i,j)

c Add max(-lowest eigenvalue,0) + add_diag to diagonal of Hessian
c which is equivalent to adding it to all except the first element of the diagonal
c of the Hamiltonian.
c We are resetting ovlp because it was destroyed by call to matinv
c We are resetting ham diagonal because it had add_diag(iadd_diag) added to it.
c We are resetting only the diagonal part of ham because that is all that is needed.
      do 6 i=1,nparmp1
        ham_sav(i,i)=ham(i,i)
        if(i.gt.1) then
c         ham(i,i)=ham(i,i)+max(-eig_min,0.d0)+add_diag(iadd_diag)
          ham(i,i)=ham(i,i)+add_diag(iadd_diag)
        endif
        do 6 j=1,nparmp1
          ovlp_sav(i,j)=ovlp(i,j)
    6 continue

c save the diagonal H/O
      do 10 i=1,nparmp1
   10   e_diag(i)=ham(i,i)/ovlp(i,i)
      write(6,'(''diagonal H/O='',50f9.3)') (e_diag(i),i=1,nparmp1)

c invert overlap matrix
c Warning: cannot use matinv because it assumes matrix is dimensioned same as filled part
      do 12 i=1,nparmp1
   12   write(6,'(''ovlp='',50f9.1)') (ovlp(i,j),j=1,nparmp1)

c     call matinv(ovlp,nparmp1,det)

ccDo cholesky decomposition
c     call chlsky(ovlp,nparmp1,MPARM,ierr)
c     if(ierr.ne.0) goto 99

ccSymmetrize decomposed matrix (needs to be done before calling uxb
ccor need to modify uxb)
c     do 20 i=1,nparmp1
c       do 20 j=i+1,nparmp1
c  20     ovlp(i,j)=ovlp(j,i)

ccSolve linear equations
c     do 30 i=1,nparmp1
c       call lxb(ovlp,nparm,MPARM,ovlp_inv()
c  30   call uxb(ovlp,nparm,MPARM,grad)

c Do cholesky decomposition using lower diagonal matrix
      call dpotrf('l',nparmp1,ovlp,MPARM,info)
c Invert using lower diagonal matrix.  ovlp is overwritten by its inverse.
      call dpotri('l',nparmp1,ovlp,MPARM,info)

c Symmetrize inverted matrix by copying lower part to upper
      do 15 i=1,nparmp1
        do 15 j=i+1,nparmp1
   15     ovlp(i,j)=ovlp(j,i)

c Temporarily check inversion
c     do 60 i=1,nparmp1
c       do 55 j=1,nparmp1
c         aa(i,j)=0
c         do 55 k=1,nparmp1
c  55     aa(i,j)=aa(i,j)+ovlp_sav(i,k)*ovlp(k,j)
c  60   write(6,'(''aa='',50f12.4)') (aa(i,j),j=1,nparmp1)

c multiply grad by inverse overlap
      do 25 i=1,nparm
        dparm(i)=0
        do 20 j=1,nparm
   20     dparm(i)=dparm(i)-ovlp(i+1,j+1)*grad(j)
   25   dparm(i)=dparm(i)/(e_diag(i+1)-e_diag(1))

      dparm_norm=0
      do 30 i=1,nparm
   30   dparm_norm=dparm_norm+dparm(i)**2
      dparm_norm=sqrt(dparm_norm/nparm)
      if(ipr_eigs.ge.1 .and. ipr_opt.ge.-4) write(6,'(''perturbation:iadd_diag,dparm_norm='',i2,9f10.5)')
     &iadd_diag,dparm_norm
c     if(ipr_eigs.ge.1 .and. ipr_opt.ge.0) write(6,'(''dparm='',100f10.6)') (dparm(i),i=1,nparm)
      write(6,'(''dparm='',100f10.6)') (dparm(i),i=1,nparm)

c Restore grad to original value
c     do 92 i=1,nparm
c  92   grad(i)=grad_sav(i)

c Restore ovlp and diagonal of ham to original value
      do 95 i=1,nparmp1
        ham(i,i)=ham_sav(i,i)
        do 95 j=1,nparmp1
   95     ovlp(i,j)=ovlp_sav(i,j)

      return
c  99 stop 'ierr.ne.0 in chlsky'
      end
c-----------------------------------------------------------------------
      subroutine linear_method(iadd_diag,ipr_eigs)
      use optim_mod
      use const_mod
      use gradhess_mod
      use contrl_opt_mod
      use optimo_mod
      implicit real*8(a-h,o-z)
      parameter(eps=1.d-9)

      common /optim2/ dparm(MPARM)
      common /linear/ ham(MPARM,MPARM),ovlp(MPARM,MPARM),coef(MPARM,MPARM)
      common /linear_sav/ ham_sav(MPARM,MPARM),ovlp_sav(MPARM,MPARM),renorm_ovlp(MPARM)
      common /linear_orig/ ovlp_orig(MPARM,MPARM)

c     dimension coef(MPARM,MPARM)
c     dimension ovlp_sav(MPARM,MPARM),renorm_ovlp(MPARM),coef(MPARM,MPARM)

      parameter(MWORK=MPARM*MPARM+1024)
c     dimension scratch(MPARM,MPARM)
      dimension work(MWORK),eigv(MPARM),eigvi(MPARM),eig_denom(MPARM)

c Since the first fn. is the current function nparm should not exceed MPARM-1
      if(nparm.gt.MPARM-1) stop 'nparm>MPARM-1'

c for high printout, set idebug=1
c     idebug=1
c for symmetrization, set isym=1
c     isym=0
c to get eigv different than the lowest, set e_target >= dabs(eigv)
c     e_target=0

      if(ipr_opt.ge.-4) write(6,'(/,''add_diag('',i1,'')='',1pd11.4)') iadd_diag,add_diag(iadd_diag)

      nparmp1=nparm+1

      do 1 i=1,nparmp1
        do 1 j=1,nparmp1
    1     coef(i,j)=0

c     if(isym.eq.1) then
c       do 3 i=1,nparmp1
c         do 3 j=1,i-1
c           ham(i,j)=0.5d0*(ham(i,j)+ham(j,i))
c   3       ham(j,i)=ham(i,j)
c     endif

c Add max(-lowest eigenvalue,0) + add_diag to diagonal of Hessian
c which is equivalent to adding it to all except the first element of the diagonal
c of the Hamiltonian.
c We are resetting ham and ovlp because they were destroyed by call to dggev.
      do 6 i=1,nparmp1
        do 6 j=1,nparmp1
          ham(i,j)=ham_sav(i,j)
          ovlp(i,j)=ovlp_sav(i,j)
          if(i.eq.j .and. i.gt.1) then
c           ham(i,j)=ham(i,j)+max(-eig_min,0.d0)+add_diag(iadd_diag)
c           ham(i,j)=ham(i,j)+max(-0.5d0*eig_min,0.d0)+add_diag(iadd_diag)
            ham(i,j)=ham(i,j)+add_diag(iadd_diag)
          endif
    6 continue

cc    if(ipr_eigs.ge.1 .and. ipr_opt.ge.1) then
c     if(ipr_opt.ge.1) then
c       do 8 i=1,nparmp1
c   8     write(6,'(''ovlp='',20g12.4)') (ovlp(i,j),j=1,i)
c       write(6,*)

c       do 9 i=1,nparmp1
c   9     write(6,'(''ham='',20g12.4)') (ham(i,j),j=1,nparmp1)
c       write(6,*)
c     endif


      call dggev('N','V',nparmp1, ham, MPARM, ovlp, MPARM, eigv, eigvi,
     &                  eig_denom, eigl, MPARM, coef, MPARM, work, -1, info )
      if(work(1).gt.MWORK) stop 'work(1).gt.MWORK'
      if(info.ne.0) stop 'info from dggev != 0'
      lwork=int(work(1))
      write(6,'(''optimal lwork='',i6)') lwork
      call dggev('N','V',nparmp1, ham, MPARM, ovlp, MPARM, eigv, eigvi,
     &                  eig_denom, eigl, MPARM, coef, MPARM, work, lwork, info )
      if(info.ne.0) write(6,*) 'Warning dggev: info =', info
      if(info.ne.0) stop 'info from dggev != 0'

c First check that MWORK is large enough
c     call dgeev('N','V',nparmp1,product,MPARM,eigv,eigvi,scratch, MPARM,
c    &            coef, MPARM, work, -1,info)
c     write(6,'(''optimal MWORK='',f8.1)') work(1)
c     if(work(1).gt.MWORK) stop 'work(1).gt.MWORK'
c     call dgeev('N','V',nparmp1,product,MPARM,eigv,eigvi,scratch, MPARM,
c    &            coef, MPARM, work, MWORK,info)
c     if(info.ne.0) write(6,*) 'Warning dgeev: info =', info

      write(6,'(''eigenvalue_denominators in dggev='',20d10.2)') (eig_denom(i),i=1,nparmp1)
      do 15 i=1,nparmp1
        eigv(i)=eigv(i)/eig_denom(i)
        eigvi(i)=eigvi(i)/eig_denom(i)
   15   write(6,'(''eigenvalue '',f16.6,'' + i* '',f10.6)')
     &  eigv(i),eigvi(i)
      write(6,*)

c The foll. is confusing
c     i0=0
c     e0=e_target
c     do 20 i=1,nparmp1
c       if(e_target+eigv(i).lt.e0) then
c         e0=e_target+eigv(i)
c         i0=i
c       endif
c  20 continue
c     if(i0.eq.0) stop 'No eigenvalues below 0'

c Find eigenvec with largest first coef.
c If 2 eigenvec have first coefs that are nearly equal in magnitude then
c choose the one for which the remaining coefs. are smaller in magnitude.
      coef_max=0
      i_coef_max=1
      do 19 i=1,nparmp1
        if(abs(coef(1,i))-coef_max.gt.-eps) then
          if(abs(coef(1,i))-coef_max.gt.eps) then
            coef_max=abs(coef(1,i))
            i_coef_max=i
           else
            sum1=0
            sum2=0
            do 18 j=2,nparmp1
              sum1=sum1+abs(coef(j,i_coef_max))
   18         sum2=sum2+abs(coef(j,i))
            if(sum2.lt.sum1) then
              coef_max=abs(coef(1,i))
              i_coef_max=i
            endif
          endif
        endif
   19 continue
      write(6,'(''The eigenvector with largest coefficient ratio is '',i4,'' eigenvalue='',f12.5,'' coef(1)='',f8.5)')
     &i_coef_max,eigv(i_coef_max),coef(1,i_coef_max)

c Find eigenvec with lowest eigenvalue that is not crazy
c If eigenvec with lowest eigenvalue is nearly degenerate with eigenvec with largest first coef
c then choose the one with largest first coef
      i_emin=0
      emin=9.d99
      do 20 i=1,nparmp1
        if(eigv(i).lt.emin .and. abs(eigv(i)-etrial).lt.10.d0) then
          if(eigv(i).lt.eigv(i_coef_max)-eps) then
            emin=eigv(i)
            i_emin=i
           else
            emin=eigv(i_coef_max)
            i_emin=i_coef_max
          endif
        endif
   20 continue
      if(i_emin.eq.0) then
        write(6,'(''Warning: No reasonable eigenvalues found; using eigenvec with largest 1st coef'')')
        i_emin=i_coef_max
      endif
      write(6,'(''The eigenvec with lowest reasonable eigenvalue is '',i4,'' eigenvalue='',f12.5,'' coef(1)='',f8.5)')
     &i_emin,eigv(i_emin),coef(1,i_emin)

c Make sure that this eigenvector has first coef of absolute magnitude 1
c     iflag=0
c     if(abs(coef_max-1.d0).gt.1.d-3) then
c       iflag=1
c       write(6,'(''iflag set to 1 in linear'')')
c     endif

c     do 22 i=1,min(nparmp1,5)
      do 22 i=1,nparmp1
   22   write(6,'(''i,coef'',i2,50f7.4)') i,(coef(j,i),j=1,nparmp1)

      if(mod(iopt/1000,10).eq.1) then
        i0=i_emin
        write(6,'(''Choosing eigenvector with lowest reasonable eigenvalue'')')
       else
        i0=i_coef_max
        write(6,'(''Choosing eigenvector with largest 1st coeff rel. to rest'')')
      endif

c Undo rescaling used to normalize H and O matrices
c Warning: For the moment we are doing it only for the ground state.
        do 23 i=1,nparmp1
   23     coef(i,i0)=coef(i,i0)/renorm_ovlp(i)

c calculate normalization
      sum=coef(1,i0)
      if(mod(iopt/10000,10).eq.0) then
        dnorm1=1.d0/coef(1,i0)
       elseif(mod(iopt/10000,10).eq.1) then
        do 24 i=2,nparmp1
          write(6,'(''ovlp_orig(1,i),coef(i,i0),ovlp_orig(1,i)*coef(i,i0)'',9f9.5)')
     &    ovlp_orig(1,i),coef(i,i0),ovlp_orig(1,i)*coef(i,i0)
   24     sum=sum-ovlp_orig(1,i)*coef(i,i0)
        dnorm1=1.d0/sum
       elseif(mod(iopt/10000,10).eq.2 .or. mod(iopt/10000,10).eq.3) then
        do 26 i=1,nparmcsf
          write(6,'(''ovlp_orig(1,i+1),coef(i+1,i0),ovlp_orig(1,i+1)*coef(i+1,i0)'',9f9.5)')
     &    ovlp_orig(1,i+1),coef(i+1,i0),ovlp_orig(1,i+1)*coef(i+1,i0)
   26     sum=sum-ovlp_orig(1,i+1)*coef(i+1,i0)
        dnorm1=1.d0/sum
        if(mod(iopt/10000,10).eq.3) dnorm2=1.d0/coef(1,i0)
       elseif(mod(iopt/10000,10).eq.4 .or. mod(iopt/10000,10).eq.5) then
        do 28 i=nparmcsf+2,nparmp1
          term1=coef(1,i0)*ovlp_orig(1,i)
          term2=coef(1,i0)*ovlp_orig(1,1)
          write(6,'(''term1,term2,ovlp_orig(1,i),ovlp_orig(i,1)'',9f10.5)') term1,term2,ovlp_orig(1,i),ovlp_orig(i,1)
c Warning: First I had 2, then nparmcsf+2, now there is a switch
          if(mod(iopt/10000,10).eq.4) then
            iparm_min=nparmcsf+2
           elseif(mod(iopt/10000,10).eq.5) then
            iparm_min=2
          endif
c         do 27 j=2,nparmp1
c         do 27 j=nparmcsf+2,nparmp1
          do 27 j=iparm_min,nparmp1
c           write(6,'(''coef(j,i0),ovlp_orig(i,j),ovlp_orig(1,j),coef(j,i0)*ovlp_orig(i,j),term1,term2'',9f10.5)')
c    &      coef(j,i0),ovlp_orig(i,j),ovlp_orig(1,j),coef(j,i0)*ovlp_orig(i,j),term1,term2
            term1=term1+coef(j,i0)*ovlp_orig(i,j)
   27       term2=term2+coef(j,i0)*ovlp_orig(1,j)
          write(6,'(''ovlp_orig(1,i),coef(i,i0),ovlp_orig(1,i)*coef(i,i0),term1,term2'',9f9.5)')
     &    ovlp_orig(1,i),coef(i,i0),ovlp_orig(1,i)*coef(i,i0),term1,term2
   28     sum=sum+(term1/term2)*coef(i,i0)
        dnorm1=1.d0/sum
      endif
      write(6,'(''eigv  '',f10.3)') eigv(i0)
      write(6,'(''coef0,dnorm1,dnorm2 '',9f16.11)') coef(1,i0),dnorm1,dnorm2
c     write(6,'(''dparm '',100f10.6)') (coef(i,i0)*dnorm1,i=2,nparmp1)
      write(6,*)

c     dparm_norm=0
      if(mod(iopt/10000,10).le.2 .or. mod(iopt/10000,10).eq.4 .or. mod(iopt/10000,10).eq.5) then
        do 30 iparm=2,nparmp1
   30     dparm(iparm-1)=coef(iparm,i0)*dnorm1
c  30     dparm_norm=dparm_norm+(coef(iparm,i0)*dnorm1)**2
       elseif(mod(iopt/10000,10).eq.3) then
        do 40 iparm=1,nparmcsf
   40     dparm(iparm)=coef(iparm+1,i0)*dnorm1
c  40     dparm_norm=dparm_norm+(coef(iparm+1,i0)*dnorm1)**2
        do 50 iparm=nparmcsf+1,nparmcsf+nparmj+nparmot+nparms
   50     dparm(iparm)=coef(iparm+1,i0)*dnorm2
c  50     dparm_norm=dparm_norm+(coef(iparm+1,i0)*dnorm2)**2
      endif

      dparm_norm=0
      do 60 iparm=1,nparm
   60   dparm_norm=dparm_norm+dparm(iparm)**2
      dparm_norm=sqrt(dparm_norm/nparm)
c     if(ipr_eigs.ge.1 .and. ipr_opt.ge.-4) write(6,'(''linear:iadd_diag,dparm_norm='',i2,9f10.5)')
      if(ipr_opt.ge.-4) write(6,'(''linear:iadd_diag,dparm_norm='',i2,9f10.5)')
     &iadd_diag,dparm_norm
      write(6,'(''dparm='',100f10.6)') (dparm(i),i=1,nparm)

c Restore ham and ovlp to original value
      do 95 i=1,nparmp1
        do 95 j=1,nparmp1
          ham(i,j)=ham_sav(i,j)
   95     ovlp(i,j)=ovlp_sav(i,j)

c Compute the norm of the move
      dnorm_move=0
      do 98 i=2,nparmp1
        do 98 j=2,nparmp1
   98     dnorm_move=dnorm_move+ovlp(i,j)*dparm(i-1)*dparm(j-1)
      dnorm_move=sqrt(dnorm_move)
      write(6,'(''Norm of the move is'',d12.4)') dnorm_move

      return
      end
c-----------------------------------------------------------------------

      subroutine update_params(iadd_diag,ipr_new,iflag)
c Written by Cyrus Umrigar

c Used to find the change in the wavefunction parameters.
c Given gradient, grad, and Hessian, hess, solve hess*x=grad for x.
c The required change in the parameters is -x.
c If ipr_new = 0 then we print new parameters without _new subscript
c            = 1 then we print new parameters with _new subscript if iflag=0
c            = 2 then we print new parameters with _new subscript
      use control_mod
      use atom_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use basis1_mod
      use contr2_mod
      use gradhess_mod
      use contrl_per_mod
      use contrl_opt_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use bparm_mod
      use optimo_mod
      implicit real*8(a-h,o-z)
      parameter(AMAX_NONLIN=100.d0)
      character*50 fmt

      common /orbpar/ oparm(MOTYPE,MBASIS,MWF)
      common /optim2/ dparm(MPARM)

c     dparm_norm=0
c     do 30 i=1,nparm
c       dparm(i)=-grad(i)
c       dparm_norm=dparm_norm+dparm(i)**2
c       grad_cal(i)=0
c       do 30 j=1,nparm
c  30     grad_cal(i)=grad_cal(i)+hess(i,j)*grad(j)
c     dparm_norm=sqrt(dparm_norm/nparm)
c     if(ipr_eigs.ge.1 .and. ipr_opt.ge.-4) write(6,'(''iadd_diag,dparm_norm='',i2,9f10.5)')
c    &iadd_diag,dparm_norm
c     if(ipr_eigs.ge.1 .and. ipr_opt.ge.1) write(6,'(''grad_cal='',9g12.4)') (grad_cal(i),i=1,nparm)
c     if(ipr_eigs.ge.1 .and. ipr_opt.ge.0) write(6,'(''dparm='',9f15.9)') (dparm(i),i=1,nparm)

c Add change to old parameters
c     open(1,file='jas_old')
c     read(1,*) ijas,isc,nspin1,nspin2,nord,ifock
c     read(1,*) norda,nordb,nordc
c     read(1,*) nparml,(nparma(ict),ict=1,nctype),(nparmb(isp),isp=nspin1,nspin2b),(nparmc(ict),ict=1,nctype)
      nparma_read=2+max(0,norda-1)
      nparmb_read=2+max(0,nordb-1)
      nparmc_read=nterms4(nordc)

c     do 45 ict=1,nctype
c  45   read(1,*) (a4(i,ict,1),i=1,nparma_read)
c     read(1,*) (b(i,1,1),i=1,nparmb_read)
c     do 46 ict=1,nctype
c  46   read(1,*) (c(i,ict,1),i=1,nparmc_read)
c     do 47 ict=1,nctype
c  47   read(1,*) (iwjasa(iparm,ict),iparm=1,nparma(ict))
c     do 48 isp=nspin1,nspin2b
c 48    read(1,*) (iwjasb(iparm,isp),iparm=1,nparmb(isp))
c     do 49 ict=1,nctype
c  49   read(1,*) (iwjasc(iparm,ict),iparm=1,nparmc(ict))
c     if(ijas.eq.4 .and. isc.le.9) call cuspinit4(1)

c Create the new parameters for this value of add_diag.
c iadd_diag=1 is done last so as not to overwrite it.
c First copy all the parameters that are read in, and then update the ones being changed.

c copy parameters
      do 50 i=1,ncsf
   50   csf_coef(i,iadd_diag)=csf_coef(i,1)

      do 51 it=1,notype
        do 51 ip=1,nbasis
   51     oparm(it,ip,iadd_diag)=oparm(it,ip,1)

      scalek(iadd_diag)=scalek(1)

      do 52 ict=1,nctype
        do 52 i=1,nparma_read
   52     a4(i,ict,iadd_diag)=a4(i,ict,1)
      do 54 isp=nspin1,nspin2b
        do 54 i=1,nparmb_read
   54     b(i,isp,iadd_diag)=b(i,isp,1)
      do 56 ict=1,nctype
        do 56 i=1,nparmc_read
   56     c(i,ict,iadd_diag)=c(i,ict,1)

c update parameters
      iparm=0
      do 58 i=1,nparmcsf
        iparm=iparm+1
   58   csf_coef(iwcsf(i),iadd_diag)=csf_coef(iwcsf(i),1)+dparm(iparm)


      do 59 it=1,notype
        do 59 ip=1,nparmo(it)
          iparm=iparm+1
   59     oparm(it,iwo(ip,it),iadd_diag)=oparm(it,iwo(ip,it),1)+dparm(iparm)

      if(nparms.eq.1) then
        iparm=iparm+1
        scalek(iadd_diag)=scalek(1)+dparm(iparm)
      endif

      do 60 ict=1,nctype
        do 60 i=1,nparma(ict)
          iparm=iparm+1
!JT          if(iwjasa(i,ict).eq.2 .and. a4(iwjasa(i,ict),ict,iadd_diag).gt.AMAX_NONLIN) stop 'probably do not want a(2) > AMAX_NONLIN'
   60     a4(iwjasa(i,ict),ict,iadd_diag)=a4(iwjasa(i,ict),ict,1)+dparm(iparm)
      do 65 isp=nspin1,nspin2b
        do 65 i=1,nparmb(isp)
          iparm=iparm+1
!JT          if(iwjasb(i,1).eq.2 .and. b(iwjasb(i,1),isp,1).gt.AMAX_NONLIN) stop 'probably do not want b(2) > AMAX_NONLIN'
   65     b(iwjasb(i,1),isp,iadd_diag)=b(iwjasb(i,1),isp,1)+dparm(iparm)
      do 70 ict=1,nctype
        do 70 i=1,nparmc(ict)
          iparm=iparm+1
   70     c(iwjasc(i,ict),ict,iadd_diag)=c(iwjasc(i,ict),ict,1)+dparm(iparm)

c set the dependent parameters by imposing cusp conditions
      if(ijas.eq.4 .and. isc.le.10) call cuspexact4(1,iadd_diag)

c Do minimal check that the move is not terrible by checking the average size of the move
c and that the nonlinear parameters in the exponent do not become < -1/scalek or too large.

c csf parameters:
      iparm=0
      iflag=0
      if(nparmcsf.gt.0) then
        dparm_norm=0
        do i=1,nparmcsf
          iparm=iparm+1
          dparm_norm=dparm_norm+dparm(iparm)**2
        enddo
        dparm_norm=sqrt(dparm_norm/nparmcsf)
        if(dparm_norm.gt.1.d0) then
          iflag=1
          write(6,'(''iadd_diag,dparm_norm=''
     &  ,i1,f9.2,'' iflag=1 because dparm_norm>1 for csf params'')') iadd_diag,dparm_norm
        endif
      endif
c orbital parameters (type 1,2,3 and 4)
      if(nparmo(1).gt.0) then
        dparm_norm=0
        do i=1,nparmo(1)
          iparm=iparm+1
          dparm_norm=dparm_norm+dparm(iparm)**2
        enddo
        dparm_norm=sqrt(dparm_norm/nparmo(1))
        if(dparm_norm.gt.1/(3*scalek(iadd_diag))) then
          iflag=1
          write(6,'(''iadd_diag,dparm_norm='',i1,f9.2,
     &'' iflag=1 because dparm_norm>1/3scalek for otype 1 params'')')
     &  iadd_diag,dparm_norm
        endif
      endif
      if(nparmo(2).gt.0) then
        dparm_norm=0
        do i=1,nparmo(2)
          iparm=iparm+1
          dparm_norm=dparm_norm+dparm(iparm)**2
        enddo
        dparm_norm=sqrt(dparm_norm/nparmo(2))
        if(ibasis.eq.4) then
          if(dparm_norm.gt.1/(3*scalek(iadd_diag))) then
            iflag=1
            write(6,'(''iadd_diag,dparm_norm='',i1,f9.2,
     &'' iflag=1 because dparm_norm>1/3scalek for otype 2 params'')') iadd_diag,dparm_norm
          endif
        elseif(ibasis.eq.6) then
          if(dparm_norm.gt.1/(3*scalek(iadd_diag))) then
            iflag=1
            write(6,'(''iadd_diag,dparm_norm='',i1,f9.2,
     &'' iflag=1 because dparm_norm>1/3scalek for otype 2 params'')') iadd_diag,dparm_norm
          endif
        elseif(ibasis.eq.5) then
          if(dparm_norm.gt..2d0) then
            iflag=1
            write(6,'(''iadd_diag,dparm_norm='',i1,f9.2,
     &'' iflag=1 because dparm_norm>0.2 for otype 2 params'')') iadd_diag,dparm_norm
          endif
        endif
      endif
      if(nparmo(3).gt.0) then
        dparm_norm=0
        do i=1,nparmo(3)
          iparm=iparm+1
          dparm_norm=dparm_norm+dparm(iparm)**2
        enddo
        dparm_norm=sqrt(dparm_norm/nparmo(3))
        if(dparm_norm.gt.2.d0) then
          iflag=1
          write(6,'(''iadd_diag,dparm_norm=''
     &,i1,f9.2,'' iflag=1 because dparm_norm>2 for otype 3 params'')') iadd_diag,dparm_norm
        endif
      endif
      if(nparmo(4).gt.0) then
        dparm_norm=0
        do i=1,nparmo(4)
          iparm=iparm+1
          dparm_norm=dparm_norm+dparm(iparm)**2
c          write(*,*) 'test: dparm(iparm)=',dparm(iparm)
        enddo
        dparm_norm=sqrt(dparm_norm/nparmo(4))
        if(dparm_norm.gt.3.d0) then
          iflag=1
          write(6,'(''iadd_diag,dparm_norm=''
     &,i1,f9.2,'' iflag=1 because dparm_norm>1 for otype 4 params'')') iadd_diag,dparm_norm
        endif
      endif
c scalek
      dparm_norm=0
      if(nparms.gt.0) then
        iparm=iparm+1
        dparm_norm=dabs(dparm(iparm))
        if(dparm_norm.gt..1d0) then
          iflag=1
          write(6,'(''iadd_diag,dparm_norm=''
     &  ,i1,f9.2,'' iflag=1 because dparm_norm>0.1 for scalek '')') iadd_diag,dparm_norm
        endif
      endif
c jastrow parameters
      if(nparmj.gt.0) then
        dparm_norm=0
        do i=1,nparmj
          iparm=iparm+1
          dparm_norm=dparm_norm+dparm(iparm)**2
        enddo
        dparm_norm=sqrt(dparm_norm/nparmj)
        if(dparm_norm.gt.max(10.d0,1.d0/(5*scalek(iadd_diag)))) then
          iflag=1
          write(6,'(''iadd_diag,dparm_norm='',i1,f9.2,
     &'' iflag=1 because dparm_norm > 10 or 1/5*scalek for jastrow params'')') iadd_diag,dparm_norm
        endif
      endif

c Do minimal check that the move is not terrible by checking the average size of the move
c and that the nonlinear parameters in the exponent do not become < -1/scalek or too large.
      iflag=0
      dparm_norm=0
      do 71 i=1,nparm
   71   dparm_norm=dparm_norm+dparm(i)**2
      dparm_norm=sqrt(dparm_norm/nparm)
      if(dparm_norm.gt.10.d0) then
        iflag=1
        write(6,'(''update_params:iadd_diag,dparm_norm='',i1,f9.2,'' iflag=1 because dparm_norm>10'')') iadd_diag,dparm_norm
      endif

      if(scalek(iadd_diag).le.0.d0) then
        write(6,'(''iadd_diag='',i1,'' iflag=1 because scalek<0'')') iadd_diag
        iflag=1
      endif

      if(nparmo(3).gt.0) then
        do ib=1,nbasis
          if(oparm(3,ib,iadd_diag).le.0.d0) then
            write(6,'(''iadd_diag='',i1,
     &         '' iflag=1 because oparm(3,..) < 0'')') iadd_diag
            iflag=1
          endif
        enddo
      endif

      if(nparmo(4).gt.0) then
        do ib=1,nbasis
          if(oparm(4,ib,iadd_diag).le.0.d0) then
            write(6,'(''iadd_diag='',i1,
     &         '' iflag=1 because oparm(4,..) < 0'')') iadd_diag
            iflag=1
          endif
        enddo
      endif

      if(isc.ne.8 .and. isc.ne.10) then
        parm2min=-scalek(1)
      else
        parm2min=-1.d0
      endif

      do 72 ict=1,nctype
        if(a4(2,ict,iadd_diag).lt.parm2min) write(6,'(''iadd_diag='',i1,'' iflag=1 because a4<-scalek'')') iadd_diag
        if(a4(2,ict,iadd_diag).gt.AMAX_NONLIN) write(6,'(''iadd_diag='',i1,'' iflag=1 because a4>AMAX_NONLIN'')') iadd_diag
   72   if(a4(2,ict,iadd_diag).lt.parm2min .or. a4(2,ict,iadd_diag).gt.AMAX_NONLIN) iflag=1
      do 74 isp=nspin1,nspin2b
        if(b(2,isp,iadd_diag).lt.parm2min) write(6,'(''iadd_diag='',i1,'' iflag=1 because b<-scalek'')') iadd_diag
        if(b(2,isp,iadd_diag).gt.AMAX_NONLIN) write(6,'(''iadd_diag='',i1,'' iflag=1 because b>AMAX_NONLIN'')') iadd_diag
   74   if(b(2,isp,iadd_diag).lt.parm2min .or. b(2,isp,iadd_diag).gt.AMAX_NONLIN) iflag=1

c Write out wavefn

      write(2,'(''iadd_diag,add_diag='',i1,1pd12.4)') iadd_diag,add_diag(iadd_diag)

      write(6,'(/,''New wave function:'')')

      if(ncsf.gt.0) then
        write(fmt,'(''('',i3,''f15.8,a)'')') ncsf
       else
        write(fmt,'(''(a)'')')
      endif
      if(ipr_new.eq.0 .or. (ipr_new.eq.1 .and. iflag.ne.0)) then
        write(6,fmt) (csf_coef(i,iadd_diag),i=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
        write(2,fmt) (csf_coef(i,iadd_diag),i=1,ncsf),' (csf_coef(icsf),icsf=1,ncsf)'
       else
        write(6,fmt) (csf_coef(i,iadd_diag),i=1,ncsf),' (csf_coef_new(icsf),icsf=1,ncsf)'
        write(2,fmt) (csf_coef(i,iadd_diag),i=1,ncsf),' (csf_coef_new(icsf),icsf=1,ncsf)'
      endif

      do it=1,notype
        write(fmt,'(''('',i3,''f15.8,a)'')') nbasis
        if(ipr_new.eq.0 .or. (ipr_new.eq.1 .and. iflag.ne.0)) then
          if(it.eq.1) then
            write(6,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_rad_pos(it,i),i=1,nbasis)'
            write(2,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_rad_pos(it,i),i=1,nbasis)'
           elseif(it.eq.2) then
            write(6,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_ang_pos(it,i),i=1,nbasis)'
            write(2,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_ang_pos(it,i),i=1,nbasis)'
           elseif(it.eq.3) then
            write(6,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_rad_width(it,i),i=1,nbasis)'
            write(2,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_rad_width(it,i),i=1,nbasis)'
           elseif(it.eq.4) then
            write(6,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_ang_width(it,i),i=1,nbasis)'
            write(2,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_ang_width(it,i),i=1,nbasis)'
          endif
         else
          if(it.eq.1) then
            write(6,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_rad_pos_new(it,i),i=1,nbasis)'
            write(2,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_rad_pos_new(it,i),i=1,nbasis)'
           elseif(it.eq.2) then
            write(6,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_ang_pos_new(it,i),i=1,nbasis)'
            write(2,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_ang_pos_new(it,i),i=1,nbasis)'
           elseif(it.eq.3) then
            write(6,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_rad_width_new(it,i),i=1,nbasis)'
            write(2,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_rad_width_new(it,i),i=1,nbasis)'
           elseif(it.eq.4) then
            write(6,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_ang_width_new(it,i),i=1,nbasis)'
            write(2,fmt) (oparm(it,i,iadd_diag),i=1,nbasis),' (floating_gauss_ang_width_new(it,i),i=1,nbasis)'
          endif
        endif
      enddo


      if(ipr_new.eq.0 .or. (ipr_new.eq.1 .and. iflag.ne.0)) then
        write(6,'(f9.6,'' 0. scalek,a21'')') scalek(iadd_diag)
        write(2,'(f9.6,'' 0. scalek,a21'')') scalek(iadd_diag)
       else
        write(6,'(f9.6,'' 0. scalek_new,a21'')') scalek(iadd_diag)
        write(2,'(f9.6,'' 0. scalek_new,a21'')') scalek(iadd_diag)
      endif

      if(nparma_read.gt.0) then
c       write(fmt,'(''('',i2,''f15.8,\'\' (a(iparmj),iparmj=1,nparma)\'\')'')') nparma_read
        write(fmt,'(''(1p'',i2,''g22.14,a)'')') nparma_read
       else
c       write(fmt,'(''(\'\' (a(iparmj),iparmj=1,nparma)\'\')'')')
        write(fmt,'(''(a)'')')
      endif
      do 80 ict=1,nctype
        if(ipr_new.eq.0 .or. (ipr_new.eq.1 .and. iflag.ne.0)) then
          write(6,fmt) (a4(i,ict,iadd_diag),i=1,nparma_read),' (a(iparmj),iparmj=1,nparma)'
          write(2,fmt) (a4(i,ict,iadd_diag),i=1,nparma_read),' (a(iparmj),iparmj=1,nparma)'
         else
          write(6,fmt) (a4(i,ict,iadd_diag),i=1,nparma_read),' (a_new(iparmj),iparmj=1,nparma)'
          write(2,fmt) (a4(i,ict,iadd_diag),i=1,nparma_read),' (a_new(iparmj),iparmj=1,nparma)'
        endif
   80 continue

      if(nparmb_read.gt.0) then
c       write(fmt,'(''('',i2,''f15.8,\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')') nparmb_read
        write(fmt,'(''(1p'',i2,''g22.14,a)'')') nparmb_read
       else
c       write(fmt,'(''(\'\' (b(iparmj),iparmj=1,nparmb)\'\')'')')
        write(fmt,'(''(a)'')')
      endif
      do 85 isp=nspin1,nspin2b
        if(ipr_new.eq.0 .or. (ipr_new.eq.1 .and. iflag.ne.0)) then
          write(6,fmt) (b(i,isp,iadd_diag),i=1,nparmb_read),' (b(iparmj),iparmj=1,nparmb)'
          write(2,fmt) (b(i,isp,iadd_diag),i=1,nparmb_read),' (b(iparmj),iparmj=1,nparmb)'
         else
          write(6,fmt) (b(i,isp,iadd_diag),i=1,nparmb_read),' (b_new(iparmj),iparmj=1,nparmb)'
          write(2,fmt) (b(i,isp,iadd_diag),i=1,nparmb_read),' (b_new(iparmj),iparmj=1,nparmb)'
        endif
   85 continue

      if(nparmc_read.gt.0) then
c       write(fmt,'(''('',i2,''f15.8,\'\' (c(iparmj),iparmj=1,nparmc)\'\')'')') nparmc_read
        write(fmt,'(''(1p'',i2,''g22.14,a)'')') nparmc_read
       else
c       write(fmt,'(''(\'\' (c(iparmj),iparmj=1,nparmc)\'\')'')')
        write(fmt,'(''(a)'')')
      endif
      do 90 ict=1,nctype
        if(ipr_new.eq.0 .or. (ipr_new.eq.1 .and. iflag.ne.0)) then
          write(6,fmt) (c(i,ict,iadd_diag),i=1,nparmc_read),' (c(iparmj),iparmj=1,nparmc)'
          write(2,fmt) (c(i,ict,iadd_diag),i=1,nparmc_read),' (c(iparmj),iparmj=1,nparmc)'
         else
          write(6,fmt) (c(i,ict,iadd_diag),i=1,nparmc_read),' (c_new(iparmj),iparmj=1,nparmc)'
          write(2,fmt) (c(i,ict,iadd_diag),i=1,nparmc_read),' (c_new(iparmj),iparmj=1,nparmc)'
        endif
   90 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine chlsky(a,n,np,ierr)
c chlsky: purpose: cholesky decomposition of a and determinant
c in: matrix a of order n stored with physical dimension np
c out: lower triangular matrix stored in lower portion of a
c note: lower triangular portion of original a is overwritten
      implicit real*8(a-h,o-z)
      dimension a(np,np)
      parameter (ZERO=0,ONE=1)

c     diag_prod=1
c     do j=1,n
c        diag_prod=diag_prod*a(j,j)
c     enddo

      det=1
      ierr=0
      do j=1,n
         if(j.gt.1) then
            jm1=j-1
            do k=j,n
               sum=0
               do ip=1,jm1
                  sum=sum+a(k,ip)*a(j,ip)
               enddo
               a(k,j)=a(k,j)-sum
            enddo
         endif

         det=det*a(j,j)
         if(a(j,j).le.ZERO) then
           write(6,'(''Warning: '',i2,'' element of a is <0'',d9.2)') j,a(j,j)
           ierr=j
           return
         endif

         s=ONE/sqrt(a(j,j))
         do k=j,n
            a(k,j)=a(k,j)*s
         enddo
      enddo

c     deta=1
c     do j=1,n
c        deta=deta*a(j,j)**2
c     enddo
c     write(6,'(''diag_prod,det,deta='',9d12.4)') diag_prod,det,deta

      return
      end
c-----------------------------------------------------------------------
      subroutine lxb(a,n,np,b)
c lxb: purpose: solve equation Lx=b, for lower triangular matrix L=a
c Golub and van Loan: algorithm 3.1.3
c in:  a = matrix of order n with physical dimension np
c      b = vector of order n
c out: b = solution of eqs.; overwites original b
      implicit real*8(a-h,o-z)
      dimension a(np,np),b(np)
      do j=1,n-1
         b(j)=b(j)/a(j,j)
         do k=j+1,n
            b(k)=b(k)-b(j)*a(k,j)
         enddo
      enddo
      b(n)=b(n)/a(n,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine uxb(a,n,np,b)
c uxb: purpose: solve equation Ux=b, for upper triangular matrix U=a
c Golub and van Loan: algorithm 3.1.4
c in:  a = matrix of order n with physical dimension np
c      b = vector of order n
c out: b = solution of eqs.; overwites original b
      implicit real*8(a-h,o-z)
      dimension a(np,np),b(np)
      do j=n,2,-1
         b(j)=b(j)/a(j,j)
         do k=1,j-1
            b(k)=b(k)-b(j)*a(k,j)
         enddo
      enddo
      b(1)=b(1)/a(1,1)
      return
      end

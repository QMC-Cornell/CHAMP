program symmetrize_den
! Symmetrize the Fourier coefficients of the density.  It assumes that the calculated Fourier components
! are accurate enough to get the sign right, but that need not be the case for the smaller Fourier aplitudes.
! Author: Cyrus Umrigar

  use types, only: i1b, ik, ik_vec, rk
  character*80 string
  integer, parameter :: MAX_SYM_OP=48
  real(rk), parameter :: eps=1.e-3
  logical rho_g_small

  real(rk) gvec(3,MAX_SYM_OP+1), gvec_abs_sorted(3,2), rho_g(MAX_SYM_OP+1), rho_g_err(MAX_SYM_OP+1), rho_g_abs_av, rho_g_err_av
  integer n_gvec, i_gvec, n_sym_op, i

  do i=1,8
    read(5,'(a)') string
    write(6,'(a)') string
    if(i.eq.6) then
      read(string,'(27x,i10)') n_gvec
    endif
  enddo

  n_sym_op=1
  read(5,*) i, gvec(:,n_sym_op), rho_g(n_sym_op), rho_g_err(n_sym_op)
  rho_g_abs_av=rho_g(n_sym_op)
  rho_g_err_av=rho_g_err(n_sym_op)
  gvec_abs_sorted(:,2)=abs(gvec(:,n_sym_op))
  call shell(gvec_abs_sorted(1,2),3)
  rho_g_small=.false.
  do i_gvec=2,n_gvec
    read(5,*) i, gvec(:,n_sym_op+1), rho_g(n_sym_op+1), rho_g_err(n_sym_op+1)
    gvec_abs_sorted(:,1)=abs(gvec(:,n_sym_op+1))
    call shell(gvec_abs_sorted(1,1),3)
    if(i_gvec.ge.3 .and. norm2(gvec_abs_sorted(:,1)-gvec_abs_sorted(:,2)).lt.eps) then
      n_sym_op=n_sym_op+1
      rho_g_abs_av=rho_g_abs_av+abs(rho_g(n_sym_op))
      rho_g_err_av=rho_g_abs_av+rho_g_err(n_sym_op)
      if(abs(rho_g(n_sym_op)).lt.2*rho_g_err(n_sym_op)) rho_g_small=.true.
    else
      rho_g_abs_av=rho_g_abs_av/n_sym_op
      rho_g_err_av=rho_g_err_av/n_sym_op**1.5_rk
      write(6,'(''n_sym_op='',i5)') n_sym_op
      do i_symp_op=1,n_sym_op
        if(rho_g_small.eqv..false.) then
          !write(6,'(i6,3f12.6,es17.8,es9.2)') i_gvec-n_sym_op+i_symp_op-1, gvec(:,i_symp_op), sign(rho_g_abs_av,rho_g(i_symp_op)), rho_g_err_av
          write(6,'(i6,3f12.6,es17.8,es9.2,i3)') i_gvec-n_sym_op+i_symp_op-1, gvec(:,i_symp_op), sign(rho_g_abs_av,rho_g(i_symp_op)), rho_g_err_av, nint(sign(1._rk,rho_g(i_symp_op)))
        else
          !write(6,'(i6,3f12.6,es17.8,es9.2)') i_gvec-n_sym_op+i_symp_op-1, gvec(:,i_symp_op), 0._rk, rho_g_err_av
          write(6,'(i6,3f12.6,es17.8,es9.2,i3)') i_gvec-n_sym_op+i_symp_op-1, gvec(:,i_symp_op), 0._rk, rho_g_err_av, nint(sign(1._rk,rho_g(i_symp_op)))
        endif
        call flush(6)
      enddo
      rho_g_abs_av=rho_g(n_sym_op+1)
      rho_g_err_av=rho_g_err(n_sym_op+1)
      gvec_abs_sorted(:,2)=gvec_abs_sorted(:,1)
      gvec(:,1)=gvec(:,n_sym_op+1)
      n_sym_op=1
      rho_g_small=.false.
    endif
  enddo
  rho_g_abs_av=rho_g_abs_av/n_sym_op
  rho_g_err_av=rho_g_err_av/n_sym_op**1.5_rk
  write(6,'(''n_sym_op='',i5)') n_sym_op
  do i_symp_op=1,n_sym_op
    if(rho_g_small.eqv..false.) then
      !write(6,'(i6,3f12.6,es17.8,es9.2)') n_gvec-n_sym_op+i_symp_op, gvec(:,i_symp_op), sign(rho_g_abs_av,rho_g(i_symp_op)), rho_g_err_av
      write(6,'(i6,3f12.6,es17.8,es9.2,i3)') n_gvec-n_sym_op+i_symp_op, gvec(:,i_symp_op), sign(rho_g_abs_av,rho_g(i_symp_op)), rho_g_err_av, nint(sign(1._rk,rho_g(i_symp_op)))
    else
      !write(6,'(i6,3f12.6,es17.8,es9.2)') n_gvec-n_sym_op+i_symp_op, gvec(:,i_symp_op), 0._rk, rho_g_err_av
      write(6,'(i6,3f12.6,es17.8,es9.2,i3)') n_gvec-n_sym_op+i_symp_op, gvec(:,i_symp_op), 0._rk, rho_g_err_av, nint(sign(1._rk,rho_g(i_symp_op)))
    endif
    call flush(6)
  enddo

  stop
end

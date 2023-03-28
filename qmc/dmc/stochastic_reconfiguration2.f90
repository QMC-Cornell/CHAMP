!subroutine stochastic_reconfiguration2(nwalk,wt,e_num,e_den)
subroutine stochastic_reconfiguration2
! Do Sandro Sorella's stochastic reconfiguration keeping just the reconfigured energy denominator and numerator equal to the true value.
! This does not require solving a set of linear equations.  If we keep any more observables constant, then we do.
! In this version we are assuming that PsiG=1.  In order to not have a sign problem, we need to impose that
! wt_ref*e_num <= 0 (assuming that the expectation value of the energy is negative).
! We have 2 options here for setting wt_ref, controlled by the abs_wt logical variable.
! 1) set wt_ref=-wt*sign(e_num),                                               , if abs_wt=true
! 2) set wt_ref=-wt*sign(e_num), if it is not sign violating and zero otherwise, if abs_wt=false (better choice)
! e_den = PsiT/PsiG
! e_num = H PsiT/PsiG
! wt    = Psi0*PsiG
! wtt   = wt * e_den = Psi0*PsiT.  For a dense wavefn, if this is negative, walker is sign violating.
! This routine can be used either for
! 1) PsiG=PsiT, in which case e_den=1 (this happens for real-space Hubbard with importance sampling)
! 2) PsiG=1, in which case e_den=PsiT (this happens for orbital-space: momentum-space Hubbard, chemistry and HEG without importance sampling)
! More generally e_den=PsiT/PsiG and e_num=E_loc*PsiT/PsiG.  wt represents Psi0*PsiG and wtt represents Psi0*PsiT.
! For sparse PsiT when e_den=0, we replace it with -eps*e_num, so that E_loc is a large negative number. (Warning: Check if this is best)

!use types, only : rk
 use, intrinsic :: iso_fortran_env, only: rk => real64 !JT: comment out for compilation on MESU
 use branch_mod, only : nwalk, eoldw, wt
 use const_mod, only : ipr

 implicit none
!integer, parameter :: rk = kind(1.0d0) ! JT: added for compilation on MESU

 integer iwalk
 real(rk) wt_ref(nwalk), wt_new(nwalk)
 real(rk) wt_sum, wt_ref_sum, wt_new_sum, e_sum, e_ref_sum, e2_ref_sum, e_new_sum, wt_av, wt_ref_av, wt_new_av, e_av, e_ref_av, e2_ref_av, e_new_av, & ! e2_sum, e2_av, &
& c, alpha, sign_cond, sign_cond_ref, sign_cond_new, min_ratio, max_ratio, p
!real(rk) :: eps=1.e-300_rk, eps2=1.e-6_rk
! false is better choice for abs_wt
 logical abs_wt /.false./
!logical use_new_wts

 if(minval(wt(1:nwalk)).ge.0.d0) return

 if(ipr.ge.1 .and. nwalk.le.60) then
   write(6,'(''  wt   e'')')
   do iwalk=1,nwalk
     write(6,'(9f10.6)') wt(iwalk), eoldw(iwalk,1)
   enddo
 endif

 wt_sum=sum(wt(1:nwalk))
 if(wt_sum.le.0._rk) then ! Switch all the wts because wt is PsiG*Psi0 and wt is PsiT*Psi0 and we want Psi0 to have +ve overlap with PsiT
   wt(1:nwalk)=-wt(1:nwalk)
   wt_sum=-wt_sum
 endif

 e_sum=dot_product(wt(1:nwalk),eoldw(1:nwalk,1))
 e_av=e_sum/wt_sum

 do iwalk=1,nwalk
     if(wt(iwalk).gt.0) then                                           ! Not sign violating
       wt_ref(iwalk)=wt(iwalk)
     else                                                               ! Sign violating
       if(abs_wt) then
         wt_ref(iwalk)=-wt(iwalk)
       else
         wt_ref(iwalk)=0
       endif
     endif
 enddo

 wt_ref_sum=sum(wt_ref(1:nwalk))
 e_ref_sum=dot_product(wt_ref(1:nwalk),eoldw(1:nwalk,1))
 e2_ref_sum=dot_product(wt_ref(1:nwalk),eoldw(1:nwalk,1)**2)

 wt_av=wt_sum/nwalk
 wt_ref_av=wt_ref_sum/nwalk
 e_ref_av=e_ref_sum/wt_ref_sum
 e2_ref_av=e2_ref_sum/wt_ref_sum

 !write(6,'(/,''wt_sum,wt_ref_sum='',9es12.4)') wt_sum,wt_ref_sum
 c=wt_sum/wt_ref_sum
 if(e2_ref_av-e_ref_av**2.ne.0._rk) then
   alpha=(e_av-e_ref_av)/(e2_ref_av-e_ref_av**2)
 else
   alpha=0 ! This can be violated if there are few walkers, so in calculation of alpha, denominator is zero but numerator is nonzero.  Just use more starting walkers to cure this.
 endif
 if(ipr.ge.1) write(6,'(/,''c, alpha='',f10.6,es14.6)') c, alpha

!write(6,'(''   wt        wt_ref     wt_new       eoldw  wt_new/wt  wt_new/wt_ref'')')
 do iwalk=1,nwalk
   wt_new(iwalk)=c*wt_ref(iwalk)*(1+alpha*(eoldw(iwalk,1)-e_ref_av))
   !!write(6,'(5f10.6,f11.6,9f10.6)') wt(iwalk),wt_ref(iwalk),wt_new(iwalk),eoldw(iwalk,1), wt_new(iwalk)/wt(iwalk), wt_new(iwalk)/wt_ref(iwalk)
   !write(6,'(3f10.6,2es12.4,es12.4,9f10.6)') wt(iwalk),wt_ref(iwalk),wt_new(iwalk),eoldw(iwalk,1), wt_new(iwalk)/wt(iwalk), wt_new(iwalk)/wt_ref(iwalk)
 enddo

 wt_new_sum=sum(wt_new(1:nwalk))
 e_new_sum=dot_product(wt_new(1:nwalk),eoldw(1:nwalk,1))
 wt_new_av=wt_new_sum/nwalk
 e_new_av=e_new_sum/wt_new_sum
!write(6,'(''wt_sum, wt_ref_sum, wt_new_sum, e_sum, e_ref_sum, e_new_sum ='',9f12.6)') wt_sum, wt_ref_sum, wt_new_sum, e_sum, e_ref_sum, e_new_sum
 if(ipr.ge.1) write(6,'(''wt_av, wt_ref_av, wt_new_av, e_av, e_ref_av, e_new_av, e-av-e_ref_av='',9f12.6)') wt_av, wt_ref_av, wt_new_av, e_av, e_ref_av, e_new_av, e_av-e_ref_av
!write(6,'(''e2_av-e_av**2, e2_ref_av-e_ref_av**2'',9es12.4)') e2_av-e_av**2, e2_ref_av-e_ref_av**2

!Check that total wt and average energy are unchanged
 if(abs((wt_sum-wt_new_sum)/wt_new_sum).gt.1.e-6_rk) then
   write(6,'(''Warning: wt_sum != wt_new_sum'',2f15.6)') wt_sum, wt_new_sum
   stop 'abs((wt_sum-wt_new_sum)/wt_new_sum).gt.1.e-6_rk'
 endif
 if(abs((e_sum-e_new_sum)/e_new_sum).gt.1.e-6_rk) then
   write(6,'(''Warning: e_sum != e_new_sum'',2f15.6)') e_sum, e_new_sum
   write(6,'(''e_av,e_ref_av,e2_ref_av,e_ref_av**2),(e_av-e_ref_av),(e2_ref_av-e_ref_av**2)'',9es14.6)') e_av,e_ref_av,e2_ref_av,e_ref_av**2,(e_av-e_ref_av),(e2_ref_av-e_ref_av**2)
   write(6,'(''abs((e_sum-e_new_sum)/e_new_sum).gt.1.e-6_rk) In calculation of alpha, denominator is zero but numerator is nonzero.  Just use more starting walkers to cure this.'')')
   stop 'abs((e_sum-e_new_sum)/e_new_sum).gt.1.e-6_rk) In calculation of alpha, denominator is zero but numerator is nonzero. Just use more starting walkers to cure this.'
 endif

! Note sign_cond_new can be < 1 only if min_ratio is < 0.  So, I use p*wt_new+(1-p)*wt, which should always have sign_cond_new=1.
! However, what is being written here is the sign_cond_new if we used wt_new rather than p*wt_new+(1-p)*wt.
 sign_cond=sum(wt(1:nwalk))/sum(abs(wt(1:nwalk)))
 sign_cond_ref=sum(wt_ref(1:nwalk))/sum(abs(wt_ref(1:nwalk)))
 sign_cond_new=sum(wt_new(1:nwalk))/sum(abs(wt_new(1:nwalk)))
 if(ipr.ge.1) write(6,'(''sign_cond, sign_cond_new='',9f16.12)') sign_cond, sign_cond_ref, sign_cond_new, sign_cond_ref-sign_cond, sign_cond_new-sign_cond
 if(sign_cond_ref.lt.0.999d0) then
   write(6,'(''wt,    wt,    wt_ref,    wt_new,   eden_tmp,   enum_tmp'')')
   do iwalk=1,nwalk
     write(6,'(9f10.6)') wt(iwalk), wt(iwalk), wt_ref(iwalk), wt_new(iwalk), eoldw(iwalk,1)
   enddo
 endif

 p=1
 min_ratio=minval(wt_new(1:nwalk)/wt(1:nwalk))
 max_ratio=maxval(wt_new(1:nwalk)/wt(1:nwalk))
 if(min_ratio.lt.0) p=1/(1-min_ratio)
 if(max_ratio.gt.2) p=min(p,1/(max_ratio-1))
 if(ipr.ge.1) write(6,'(''min, max wt_new/wt, wt_new/wt_ref'',9f10.6)') min_ratio, max_ratio, minval(wt_new(1:nwalk)/wt_ref(1:nwalk)),maxval(wt_new(1:nwalk)/wt_ref(1:nwalk))

 if(ipr.ge.1 .and. (min_ratio.lt.0 .or. max_ratio.gt.2)) then
   if(nwalk.le.2000) then
!    write(6,'(''Warning: minval(wt_new(1:nwalk)/wt(1:nwalk)), maxval(wt_new(1:nwalk)/wt(1:nwalk)), p='',2es12.4, f6.3)') minval(wt_new(1:nwalk)/wt(1:nwalk)), maxval(wt_new(1:nwalk)/wt(1:nwalk)), p
     write(6,'(''Warning: min_ratio, max_ratio, p='',2es12.4, f6.3)') min_ratio, max_ratio, p
     write(6,'(''   wt        wt_ref     wt_new         eoldw  wt_new/wt  wt_new/wt_ref'')')
     do iwalk=1,nwalk
       write(6,'(3f10.6,f11.6,9f10.6)') wt(iwalk),wt_ref(iwalk),wt_new(iwalk),eoldw(iwalk,1), wt_new(iwalk)/wt(iwalk), wt_new(iwalk)/wt_ref(iwalk)
     enddo
   endif
 endif

 wt(1:nwalk)=p*wt_new(1:nwalk)+(1-p)*wt(1:nwalk)
 
end subroutine stochastic_reconfiguration2

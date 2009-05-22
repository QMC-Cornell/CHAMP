module gradhess_mod

 use constants_mod
 implicit none
 save

 double precision grad(MPARM),grad_var(MPARM),hess(MPARM,MPARM),hess_var(MPARM,MPARM),gerr(MPARM)
 double precision add_diag(3),energy(3),energy_sigma(3),energy_err(3),force(3),force_err(3)
 double precision eig_min,eig_max,p_var,tol_energy
 integer nopt_iter,nblk_max
 
end module gradhess_mod

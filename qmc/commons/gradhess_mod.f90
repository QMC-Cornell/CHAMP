module gradhess_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: grad(:),grad_var(:),hess(:,:),hess_var(:,:),gerr(:)
 double precision add_diag(3),energy(3),energy_sigma(3),energy_err(3),force(3),force_err(3)
 double precision eig_min,eig_max,p_var,tol_energy
 integer nopt_iter,nblk_max
 
end module gradhess_mod

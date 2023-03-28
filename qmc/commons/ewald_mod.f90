module ewald_mod

 implicit none
 save 

 double precision, allocatable :: b_coul(:),y_coul(:),b_coul_sim(:),y_coul_sim(:),b_psp(:,:),y_psp(:,:),b_jas(:),y_jas(:)
 double precision, allocatable :: cos_n_sum(:),sin_n_sum(:),cos_e_sum(:),sin_e_sum(:)
 double precision, allocatable :: cos_e_sum_sim(:),sin_e_sum_sim(:),cos_p_sum(:),sin_p_sum(:)
 integer :: NG1DX

end module ewald_mod


module gradhessder_mod

 implicit none
 save

 double precision, allocatable :: dj(:),dj_e(:),dj_de(:,:),dj_dj(:,:),dj_dj_e(:,:)
 double precision, allocatable :: de(:),d2j(:,:),d2j_e(:,:),de_e(:),e2(:),dj_e2(:),de_de(:,:)
 double precision, allocatable :: w_i(:),w_i_e(:)

end module gradhessder_mod


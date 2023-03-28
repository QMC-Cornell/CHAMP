module orbital_num_mod

 implicit none
 save

 double precision, allocatable :: orb_num(:,:,:,:),xorb_grid(:),yorb_grid(:)
 double precision :: sizex,sizey,hx,hy,hxi,hyi
 integer ngrid_orbx,ngrid_orby
 integer :: ict(6)

end module orbital_num_mod

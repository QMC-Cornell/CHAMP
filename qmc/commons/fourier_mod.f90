module fourier_mod

 implicit none
 save

 double precision, allocatable :: fourierrk_u(:,:),fourierrk_d(:,:)
 double precision, allocatable :: fourierrk_t(:,:),fourierkk_u(:,:),fourierkk_d(:,:),fourierkk_t(:,:)
 double precision              :: delk1,delk2,fmax1,fmax2
 integer                       :: ifourier
 integer, parameter            :: NAK1=40
 integer, parameter            :: NAK2=1

end module fourier_mod

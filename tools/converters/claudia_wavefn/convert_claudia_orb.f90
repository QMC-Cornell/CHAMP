program convert_claudia_orb
! Remove extra 3s orbitals from Claudia's wavefn. for hexatriene
  implicit real*8(a-h,o-z)
  parameter(nbasis=142,norb=136)
  dimension coef(nbasis)

  do j=1,norb+1
    read(5,*) coef
    write(6,'(136es15.7)') (coef(i),i=1,2),(coef(i),i=4,19), (coef(i),i=21,36),(coef(i),i=38,53),(coef(i),i=55,70), (coef(i),i=72,87),(coef(i),i=89,142)
  enddo

  stop
end program convert_claudia_orb

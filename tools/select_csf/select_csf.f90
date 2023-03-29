program select_csf
implicit real*8(a-h,o-z)
integer icsf, icsf_new, det_in_icsf
integer, allocatable :: iworbd(:,:), ndet_in_csf(:), ndet_in_csf_new(:), iwdet_in_csf(:,:), iflag(:), idetnew_from_idetold(:), idet_index(:), idet_energy(:)
real*8, allocatable :: csf_coef(:), csf_coef_new(:), cdet_in_csf(:,:)
!real*8, parameter :: csf_coef_min=0.01d0
character*32 fmt

read(5,*) nup, ndn, ndet
read(5,*) csf_coef_min
nelec=nup+ndn

allocate(iworbd(nelec,ndet))
allocate(iflag(ndet))
allocate(idetnew_from_idetold(ndet))
allocate(idet_index(ndet))
allocate(idet_energy(ndet))

do idet=1,ndet
  read(5,*) (iworbd(j,idet),j=1,nelec), idet_index(idet), idet_energy(idet)
enddo

read(5,*) ncsf
allocate(csf_coef(ncsf))
allocate(ndet_in_csf(ncsf))
read(5,*) csf_coef(1:ncsf)
read(5,*) ndet_in_csf(1:ncsf)
allocate (iwdet_in_csf(maxval(ndet_in_csf),ncsf))
allocate (cdet_in_csf(maxval(ndet_in_csf),ncsf))
do icsf = 1, ncsf
  read(5,*) (iwdet_in_csf(det_in_icsf,icsf),det_in_icsf=1,ndet_in_csf(icsf))
  read(5,*) (cdet_in_csf(det_in_icsf,icsf),det_in_icsf=1,ndet_in_csf(icsf))
enddo

iflag(1:ndet)=0
do icsf = 1, ncsf
  if(abs(csf_coef(icsf)).gt.csf_coef_min) then
    do idet_in_csf=1,ndet_in_csf(icsf)
      iflag(iwdet_in_csf(idet_in_csf,icsf))=1
      if(iwdet_in_csf(idet_in_csf,icsf).gt.ndet) stop 'iwdet_in_csf(idet_in_csf,icsf) > ndet'
    enddo
  endif
enddo

idetnew_from_idetold(1:ndet)=0
ndet_new=0
write(fmt,'("(",i4,"i5,3x,",i4,"i5,3x,2i5)")') nup, ndn
!write(6,'(''fmt= '',a)') fmt ; call flush(6)
do idet=1,ndet
  if(iflag(idet).eq.1) then
    ndet_new=ndet_new+1
    idetnew_from_idetold(idet)=ndet_new
    write(6,fmt) (iworbd(j,idet),j=1,nelec), idet_index(idet), idet_energy(idet)
  endif
enddo

!write(6,*)
!write(6,'(2i5)') (idet, idetnew_from_idetold(idet), idet=1,ndet)
!write(6,*)

ncsf_new=count(abs(csf_coef).gt.csf_coef_min)
write(6,'(i5,a)') ncsf_new, ' ncsf'
!ncsf_new=count(abs(csf_coef).gt..1*csf_coef_min)
!write(6,'(''ncsf_new='',i5)') ncsf_new

allocate(csf_coef_new(ncsf_new))
allocate(ndet_in_csf_new(ncsf_new))

icsf_new=0
do icsf = 1, ncsf
  if(abs(csf_coef(icsf)).gt.csf_coef_min) then
    icsf_new=icsf_new+1
    csf_coef_new(icsf_new)=csf_coef(icsf)
    ndet_in_csf_new(icsf_new)=ndet_in_csf(icsf)
    write(fmt,'("(",i4,"i5,a)")') ndet_in_csf(icsf)
!   write(6,*) fmt ; call flush(6)
    write(6,fmt) (idetnew_from_idetold(iwdet_in_csf(det_in_icsf,icsf)),det_in_icsf=1,ndet_in_csf(icsf)), ' (iwdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
    write(fmt,'("(",i4,"f12.8,a)")') ndet_in_csf(icsf)
    write(6,fmt) (cdet_in_csf(det_in_icsf,icsf),det_in_icsf=1,ndet_in_csf(icsf)), ' (cdet_in_csf(idet_in_csf,icsf),idet_in_csf=1,ndet_in_csf(icsf))'
! write(6,'(a,i5)') ' CSF # ',icsf
! write(6,'(a,200i4)') ' determinants in CSF:',(iwdet_in_csf(det_in_icsf,icsf),det_in_icsf=1,ndet_in_csf(icsf))
! write(6,'(a,200f8.5)') ' coefficients:',(cdet_in_csf(det_in_icsf,icsf),det_in_icsf=1,ndet_in_csf(icsf))
  endif
enddo

write(6,*)
write(fmt,'("(",i5,"f12.8,a)")') ncsf_new
write(6,fmt) csf_coef_new(1:ncsf_new), ' (csf_coef(icsf),icsf=1,ncsf)'
write(fmt,'("(",i5,"i12,a)")') ncsf_new
write(6,fmt) ndet_in_csf_new(1:ncsf_new), ' (ndet_in_csf(icsf),icsf=1,ncsf)'
write(6,'(''Move the previous 2 lines after the 1st line'')')

write(6,'(/,''ndet_new='',i6)') ndet_new
write(6,'(''Use value in previous line for ndet'')')

stop
end

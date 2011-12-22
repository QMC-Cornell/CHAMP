program convert_claudia_csf
  implicit real*8(a-h,o-z)
  parameter(MDET=2000,MCSF=500,MDET_in_CSF=20)
  dimension ndet_in_csf(MCSF),iwdet_in_csf(MDET_in_CSF),cdet_in_csf(MDET_in_CSF)

  read(5,*) ncsf
  if(ncsf.gt.mcsf) stop 'ncsf>mcsf'

  do icsf=1,ncsf
    read(5,*) ndet_in_csf(icsf)
    if(ndet_in_csf(icsf).gt.MDET_in_CSF) stop 'ndet_in_csf(icsf) > MDET_in_CSF'
    do idet_in_csf=1,ndet_in_csf(icsf)
      read(5,*) iwdet_in_csf(idet_in_csf), cdet_in_csf(idet_in_csf)
    enddo
    write(6,'(20i10)') (iwdet_in_csf(idet_in_csf),idet_in_csf=1,ndet_in_csf(icsf))
    write(6,'(20f10.6)') (cdet_in_csf(idet_in_csf),idet_in_csf=1,ndet_in_csf(icsf))
  enddo
  write(6,'(9000i3)') (ndet_in_csf(icsf),icsf=1,ncsf)

  stop
end program convert_claudia_csf

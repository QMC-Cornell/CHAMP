!------------------------------------------------------------------!
! Plane waves -> blip conversion utility                           !
! Dario Alfe` 3.2001                                               !
!                                                                  !       
! Changes                                                          !       
! -------                                                          !
! 5.2002  MDT - Tidied up. Made conform to f90 standard.           ! 
! 11.2003 DA  - Added fast fourier transforms                      !       
! 10.2004 NDD - Tried to make the program a bit more user-friendly.!
! 10.2004 NDD - Some modifications for spin-polarized systems.     !
!------------------------------------------------------------------!


MODULE globals
! Some global variables.
 IMPLICIT NONE
 DOUBLE PRECISION xmul
 INTEGER igrad_lap
 LOGICAL ltest
END MODULE


 PROGRAM blip_wrapper
  USE globals
  IMPLICIT NONE
  INTEGER nbasis,nwvec,nkvec,idum,ierr
  INTEGER :: io=10
  DOUBLE PRECISION xk(3)
  DOUBLE PRECISION,PARAMETER :: eps=1.d-6
  LOGICAL pwreal
  CHARACTER(80) :: sline

  write(6,'(/''************************************************************'')')
  write(6,'(''        PLANE WAVES  ---->   BLIPS   TRANSFORMATION         '',&
   &/)')

  open(unit=io,file='pwfn.data',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening pwfn.data.  Stopping.'
   stop
  endif

! Rapid scan of file to get relevant array dimensions. Allocate arrays.

  call skipio(io,34)
  read(io,*)nbasis ; call skipio(io,nbasis+9)
  read(io,*)nwvec  ; call skipio(io,nwvec+5)
  read(io,*)nkvec          
  call skipio(io,1)
  read(io,*)idum,idum,idum,xk(:)
  call skipio(io,3)
  read(io,fmt='(a)')sline
  if(scan(sline,",")/=0)then
   pwreal=.false. ! No inversion symmetry --> complex PW coefficients
  else
   pwreal=.true.  ! Inversion symmetry --> real PW coefficients
  endif

  close(io)

  write(6,*)'Please enter XMUL (size of the blips grid; should be >= 1.0),'
  write(6,*)'IGRAD_LAP (0=approximate orbitals only; 1=orbitals+Laplacian),'
  write(6,*)'          (2=orbitals+gradient+Laplacian)'
  write(6,*)'and LTEST (.t. if you want to test the quality of the blips basis,'
  write(6,*)'.f. otherwise; see manual for information on test).'
  write(6,*)'XMUL, IGRAD_LAP, LTEST ?'
  read(5,*)xmul,igrad_lap,ltest
  write(6,*)

  if(pwreal)then
   if(nkvec==1.and.abs(xk(1))<eps.and.abs(xk(2))<eps.and.abs(xk(3))<eps)then
    call blipr
   else
    call blipkr
   endif
  else
   if(nkvec==1.and.abs(xk(1))<eps.and.abs(xk(2))<eps.and.abs(xk(3))<eps)then
    call blipp
   else
    call blipk
   endif
  endif

  write(6,*)'Program finished.  bwfn.data has been generated.'
  write(6,*)
     
 END PROGRAM blip_wrapper


 SUBROUTINE blipp
  USE globals
  USE singleton
  IMPLICIT NONE
  INTEGER i,j,l,nr(3),ibnd,ig,ig1,npw_not_zero,npwx,lx,ly,lz,nks,ik,&
   &nat,na,atomno,nbnd(2),ispin,nspin,nelec,isign,icount,indice,ierr
  INTEGER :: io=10,io2=11
  INTEGER,PARAMETER :: n_points_for_test=100
  INTEGER,ALLOCATABLE :: minus_g(:)
  DOUBLE PRECISION d,d2,dg(3),x,y,z,k,blip(5),gxmax,gymax,gzmax,dax,day,&
   &daz,r(3),xb(5),xd(5),xbd(5),at(3,3),bg(3,3),tv,tc,tv0,tc0,tmp,ranfx,tau(3),&
   &etot,ewld,ecut,xk(3),et,arg,norm_real,norm_imag
  DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0,tpi=2.d0*pi,eps=1.d-6
  DOUBLE PRECISION,ALLOCATABLE :: g(:,:), gcart(:,:), g2(:), avc(:,:,:), gamma(:)
  DOUBLE PRECISION,ALLOCATABLE :: avc2(:,:,:) !coefficient array for Laplacian
  DOUBLE PRECISION,ALLOCATABLE :: avc3(:,:,:),avc4(:,:,:),avc5(:,:,:) !gradient
  COMPLEX(KIND=KIND(0.d0)) :: eigr
  COMPLEX(KIND=KIND(0.d0)),PARAMETER :: imag=(0.d0,1.d0)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: evc(:), data(:)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: data2(:) !FFT array for Laplacian
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: data3(:),data4(:),data5(:) !gradient
  CHARACTER(80) title
  LOGICAL lsda,lreal,limag
  LOGICAL,ALLOCATABLE :: use_this_g(:)

  write(6,*)'Complex PW coefficients; GAMMA point'
  write(6,*)

  open(io2,file='pwfn.data',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening pwfn.data. Stopping.'
   stop
  endif
  open(io,file='bwfn.data',status='unknown',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening bwfn.data. Stopping.'
   stop
  endif
 
  call bwfn_info(io,io2,lsda)

  call skipio(io2,1) 
  write(io,'(a)')'Primitive lattice vectors (au) '
  read(io2,*)at(:,1)
  read(io2,*)at(:,2)
  read(io2,*)at(:,3)
  write(io,'(3f20.12)')at(:,1)
  write(io,'(3f20.12)')at(:,2)
  write(io,'(3f20.12)')at(:,3)

  call inve(at,bg)
  bg=transpose(bg)

  call skipio(io2,4) 
  write(io,'(a)')' '
  write(io,'(a)')'G VECTORS'
  write(io,'(a)')'---------'
  write(io,'(a)')'Number of G-vectors'
  read(io2,*) npwx
  write(io,*) npwx
  call skipio(io2,1) 
  write(io,'(a)')'Gx Gy Gz (au)'

  allocate(g(3,npwx),gcart(3,npwx),g2(npwx))
  gxmax=0 ; gymax=0 ; gzmax=0
  do ig=1,npwx
   read(io2,*)g(:,ig)
   write(io,'(3f20.12)') g(:,ig)
   g2(ig)=dot_product(g(:,ig),g(:,ig))
   gcart(:,ig)=g(:,ig)
   g(:,ig)=matmul(transpose(at),g(:,ig))/tpi
   if(abs(g(1,ig))>gxmax)gxmax=abs(g(1,ig))
   if(abs(g(2,ig))>gymax)gymax=abs(g(2,ig))
   if(abs(g(3,ig))>gzmax)gzmax=abs(g(3,ig))
  enddo
!  write(6,'('' Gmax:'',3f20.10)')gxmax, gymax, gzmax
  nr(1)=2*nint(gxmax)+2
  nr(2)=2*nint(gymax)+2
  nr(3)=2*nint(gzmax)+2
  nr=nr*xmul
! The following is to make nr an even number
  if(mod(nr(1),2)/=0)nr(1)=nr(1)+1
  if(mod(nr(2),2)/=0)nr(2)=nr(2)+1
  if(mod(nr(3),2)/=0)nr(3)=nr(3)+1
  write(6,'('' Grid : '',3i5)')nr
  write(6,*)
  dax=1.d0/dble(nr(1))
  day=1.d0/dble(nr(2))
  daz=1.d0/dble(nr(3))
  write(io,'(a)')'Blip grid'
  write(io,'(3i4)')nr

  call skipio(io2,4) 
  write(io,'(a)')' '
  write(io,'(a)')'WAVE FUNCTION'
  write(io,'(a)')'-------------'
  write(io,'(a)')'Number of k-points'
  read(io2,*)nks
  if(nks>1)then
   write(6,*)' !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
   write(6,*)' !!! WARNING, K-points, we should call blipk instead !!!' 
   write(6,*)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
  endif
  write(io,*)nks
  if(lsda)then
   nspin=2
  else
   nspin=1
  endif

  allocate(evc(npwx),use_this_g(npwx),minus_g(npwx),gamma(npwx))
  select case (igrad_lap)
    case(0)
        allocate( data(nr(1)*nr(2)*nr(3)) )
        allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1)); avc = 0.d0
    case(1)
        allocate( data(nr(1)*nr(2)*nr(3)) )
        allocate( data2(nr(1)*nr(2)*nr(3)) )
        allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1)); avc = 0.d0
        allocate(avc2(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1)); avc2 = 0.d0
    case(2)
        allocate( data(nr(1)*nr(2)*nr(3)) )
        allocate( data2(nr(1)*nr(2)*nr(3)) )
        allocate( data3(nr(1)*nr(2)*nr(3)) )
        allocate( data4(nr(1)*nr(2)*nr(3)) )
        allocate( data5(nr(1)*nr(2)*nr(3)) )
        allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1)); avc = 0.d0
        allocate(avc2(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1)); avc2 = 0.d0
        allocate(avc3(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1)); avc3 = 0.d0
        allocate(avc4(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1)); avc4 = 0.d0
        allocate(avc5(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1)); avc5 = 0.d0
  end select



! Calculating gamma
  do ig=1,npwx
   k=dax*g(1,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=1.5d0
   endif
   k=day*g(2,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=gamma(ig)*3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=gamma(ig)*1.5d0
   endif
   k=daz*g(3,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=gamma(ig)*3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=gamma(ig)*1.5d0
   endif
   gamma(ig)=1.d0/gamma(ig)
  enddo

! Set up to use only +G
  use_this_g=.true.
  minus_g=0
  minus_g(1)=1
  i=0
  if(abs(g(1,1))>eps.or.abs(g(1,1))>eps.or.abs(g(1,1))>eps) &
   stop 'First vector has to be the zero vector, reorder vectors in pwfn.data&
   & please'
  do ig=1,npwx 
   if(use_this_g(ig))then
! Find -G
    do_ig: do ig1=ig+1,npwx
     if(abs(g(1,ig)+g(1,ig1))<eps.and.abs(g(2,ig)+g(2,ig1))<eps.and. &
       &abs(g(3,ig)+g(3,ig1))<eps)then
      use_this_g(ig1)=.false.
      minus_g(ig)=ig1
      i=i+1
      exit do_ig
     endif
    enddo do_ig
   endif
  enddo
!  if(2*i+1/=npwx)stop 'I could not find all +G -G vectors.'

  do ik=1,nks
   call skipio(io2,1)
   write(io,'(a)')'k-point # ; # of bands (up spin/down spin) ; k-point coords&
    & (au)'
   read(io2,*)i,nbnd(1),nbnd(2),xk(:)
   if(abs(xk(1))>eps.or.abs(xk(2))>eps.or.abs(xk(3))>eps)then
    stop 'We should not use this routine but blipk.'
   endif
   write(io,'(3i4,3f20.16)')ik,nbnd(1),nbnd(2),xk(:)
   do ispin=1,nspin 
    if(ltest)then
     if(lsda)write(6,*)'Spin : ',ispin
     write(6,*)'Band ; alpha ; alpha (Laplacian) ; alpha (Gradient)'
    endif
    do ibnd=1,nbnd(ispin)
     call skipio(io2,1)
     write(io,'(a)')'Band, spin, eigenvalue (au)'
     read(io2,*)i,i,et 
     write(io,'(2i5,f20.12)')ibnd,ispin,et
     call skipio(io2,1)
! Read in the wavefunctions
     do ig=1,npwx
      read(io2,*)evc(ig)
     enddo
! Decide if take the real or the imaginary part
     norm_real=dble(evc(1))**2
     norm_imag=aimag(evc(1))**2
     do ig=2,npwx
      if(use_this_g(ig))then
       norm_real=norm_real+2*(abs(evc(ig)+conjg(evc(minus_g(ig))))/2)**2
       norm_imag=norm_imag+2*(abs(evc(ig)-conjg(evc(minus_g(ig))))/2)**2
      endif
     enddo
     lreal=.false. ; limag=.false.
     if(abs(norm_real-1.d0)<eps)lreal=.true.
     if(abs(norm_imag-1.d0)<eps)limag=.true.
     if(lreal.and.limag)stop 'WF both pure real and pure imaginary ???'
 
     if(lreal)write(6,*)'Real wave function, band : ',ibnd
     if(limag)write(6,*)'Imaginary wave function, band :',ibnd
!       write(6,*)norm_real,norm_imag
 
     if(.not.lreal.and..not.limag)then 
      if(norm_real>norm_imag)then
       lreal=.true.
      else
       limag=.true.
      endif
     endif

     write(io,'(a)') 'Blip coefficients'
     select case (igrad_lap)
     case(0)
        data=0.d0
        npw_not_zero=0
        if(lreal)then
         data(1)=dble(evc(1))/3.375d0
         do ig=2,npwx
          if(use_this_g(ig))then
           npw_not_zero=npw_not_zero+1
           icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
            &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
            &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
           data(icount)=(evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig) 
          endif
         enddo
! avc = \sum_g gamma(g) * Re[ evc(g) * exp( i.g.r ) ]
! If the wave function is real then evc(g) = conjg(evc(-g)) and the sum
! can be done only on +g because gamma(g) = gamma(-g).
! If the WF is not real but we want to take the real part then we take
! evc(ig) + conjg( evc(-g) ).
!      write(6,*)'Taking the REAL part of WF # ',ibnd,'; # of PW coefficients: ',&
!       &npw_not_zero
        elseif(limag)then
         data(1)=cmplx(0.d0,aimag(evc(1))/3.375d0)
         do ig=2,npwx
          if(use_this_g(ig))then
           npw_not_zero=npw_not_zero+1
           icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
            (nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
            (nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
           data(icount)=(evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig) 
          endif
         enddo
!    If the WF is not imaginary but we want to take the imaginary part then we take
!    evc(ig)-conjg(evc(-g)) 
!      write(6,*)'Taking the IMAGINARY part of WF # ',ibnd,'; # of PW &
!      &coefficients: ',npw_not_zero
        else
         stop 'Wavefunction is not real nor imaginary, what should I do?'
        endif
     
        call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
        data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
        icount=1
        do lz=0,nr(3)-1
         do ly=0,nr(2)-1
          do lx=0,nr(1)-1
           if(lreal)then
            avc(lx,ly,lz)=data(icount)
           else
            avc(lx,ly,lz)=aimag(data(icount))
           endif
           icount=icount+1
          enddo
         enddo
        enddo

        do lx=0,nr(1)-1
         do ly=0,nr(2)-1
          do lz=0,nr(3)-1 
           write(io,*)avc(lx,ly,lz)
          enddo
         enddo
        enddo
           
       if(ltest)then
        xb=0 ; xd=0 ; xbd=0
        do i=1,n_points_for_test 
         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)
         call blip3d(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d2=0
         dg=0
         if(lreal)then
          d=dble(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)+conjg(evc(minus_g(ig))))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+dble(eigr)
            d2=d2-dble(eigr)*g2(ig)
            dg(1)=dg(1)-tpi*g(1,ig)*aimag(eigr)
            dg(2)=dg(2)-tpi*g(2,ig)*aimag(eigr)
            dg(3)=dg(3)-tpi*g(3,ig)*aimag(eigr)
           endif
          enddo
         else
          d=aimag(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)-conjg(evc(minus_g(ig))))* &
            &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+aimag(eigr)
            d2=d2-aimag(eigr)*g2(ig) 
            dg(1)=dg(1)+tpi*g(1,ig)*dble(eigr)
            dg(2)=dg(2)+tpi*g(2,ig)*dble(eigr)
            dg(3)=dg(3)+tpi*g(3,ig)*dble(eigr)
           endif
          enddo
         endif
         dg=matmul(bg,dg)
         xb(1)=xb(1)+blip(1)**2
         xd(1)=xd(1)+d**2
         xbd(1)=xbd(1)+blip(1)*d
         xb(2)=xb(2)+blip(2)**2
         xd(2)=xd(2)+d2**2
         xbd(2)=xbd(2)+blip(2)*d2
         xb(3)=xb(3)+blip(3)**2
         xd(3)=xd(3)+dg(1)**2
         xbd(3)=xbd(3)+blip(3)*dg(1)
         xb(4)=xb(4)+blip(4)**2
         xd(4)=xd(4)+dg(2)**2
         xbd(4)=xbd(4)+blip(4)*dg(2)
         xb(5)=xb(5)+blip(5)**2
         xd(5)=xd(5)+dg(3)**2
         xbd(5)=xbd(5)+blip(5)*dg(3)
        enddo
       endif
       if(ltest)then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2,(xbd(3)/sqrt(xb(3)*xd(3)))**2,     & 
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2,(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif
     case(1) !orbitals + Laplacian
        data=0.d0
        data2=0.d0
        npw_not_zero=0
        if(lreal)then
         data(1)=dble(evc(1))/3.375d0
         do ig=2,npwx
          if(use_this_g(ig))then
           npw_not_zero=npw_not_zero+1
           icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
            &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
            &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
           data(icount)=(evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig) 
           data2(icount)=-g2(ig)*(evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig) 
          endif
         enddo
! avc = \sum_g gamma(g) * Re[ evc(g) * exp( i.g.r ) ]
! If the wave function is real then evc(g) = conjg(evc(-g)) and the sum
! can be done only on +g because gamma(g) = gamma(-g).
! If the WF is not real but we want to take the real part then we take
! evc(ig) + conjg( evc(-g) ).
!      write(6,*)'Taking the REAL part of WF # ',ibnd,'; # of PW coefficients: ',&
!       &npw_not_zero
        elseif(limag)then
         data(1)=cmplx(0.d0,aimag(evc(1))/3.375d0)
         do ig=2,npwx
          if(use_this_g(ig))then
           npw_not_zero=npw_not_zero+1
           icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
            (nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
            (nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
           data(icount)=(evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig) 
           data2(icount)=-g2(ig)*(evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig) 
          endif
         enddo
!    If the WF is not imaginary but we want to take the imaginary part then we take
!    evc(ig)-conjg(evc(-g)) 
!      write(6,*)'Taking the IMAGINARY part of WF # ',ibnd,'; # of PW &
!      &coefficients: ',npw_not_zero
        else
         stop 'Wavefunction is not real nor imaginary, what should I do?'
        endif
     
        call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
        call fftn(data2,(/nr(1),nr(2),nr(3)/),inv=.true.)
        data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
        data2=data2*sqrt(dble(nr(1)*nr(2)*nr(3)))
        icount=1
        do lz=0,nr(3)-1
         do ly=0,nr(2)-1
          do lx=0,nr(1)-1
           if(lreal)then
            avc(lx,ly,lz)=data(icount)
            avc2(lx,ly,lz)=data2(icount)
           else
            avc(lx,ly,lz)=aimag(data(icount))
            avc2(lx,ly,lz)=aimag(data2(icount))
           endif
           icount=icount+1
          enddo
         enddo
        enddo

        write(io,'(a)')'Orbital'
        do lx=0,nr(1)-1
         do ly=0,nr(2)-1
          do lz=0,nr(3)-1 
           write(io,*)avc(lx,ly,lz)
          enddo
         enddo
        enddo
        write(io,'(a)')'Laplacian'
        do lx=0,nr(1)-1
         do ly=0,nr(2)-1
          do lz=0,nr(3)-1 
           write(io,*)avc2(lx,ly,lz)
          enddo
         enddo
        enddo
           
       if(ltest)then
        xb=0 ; xd=0 ; xbd=0
        do i=1,n_points_for_test 
         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)
         call blip3d(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d2=0
         dg=0
         if(lreal)then
          d=dble(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)+conjg(evc(minus_g(ig))))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+dble(eigr)
            d2=d2-dble(eigr)*g2(ig)
            dg(1)=dg(1)-tpi*g(1,ig)*aimag(eigr)
            dg(2)=dg(2)-tpi*g(2,ig)*aimag(eigr)
            dg(3)=dg(3)-tpi*g(3,ig)*aimag(eigr)
           endif
          enddo
         else
          d=aimag(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)-conjg(evc(minus_g(ig))))* &
            &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+aimag(eigr)
            d2=d2-aimag(eigr)*g2(ig) 
            dg(1)=dg(1)+tpi*g(1,ig)*dble(eigr)
            dg(2)=dg(2)+tpi*g(2,ig)*dble(eigr)
            dg(3)=dg(3)+tpi*g(3,ig)*dble(eigr)
           endif
          enddo
         endif
         dg=matmul(bg,dg)
         xb(1)=xb(1)+blip(1)**2
         xd(1)=xd(1)+d**2
         xbd(1)=xbd(1)+blip(1)*d
         xb(2)=xb(2)+blip(2)**2
         xd(2)=xd(2)+d2**2
         xbd(2)=xbd(2)+blip(2)*d2
         xb(3)=xb(3)+blip(3)**2
         xd(3)=xd(3)+dg(1)**2
         xbd(3)=xbd(3)+blip(3)*dg(1)
         xb(4)=xb(4)+blip(4)**2
         xd(4)=xd(4)+dg(2)**2
         xbd(4)=xbd(4)+blip(4)*dg(2)
         xb(5)=xb(5)+blip(5)**2
         xd(5)=xd(5)+dg(3)**2
         xbd(5)=xbd(5)+blip(5)*dg(3)
        enddo
       endif
       if(ltest)then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2,(xbd(3)/sqrt(xb(3)*xd(3)))**2,     & 
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2,(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif
     case(2) !orbital, Laplacian and gradient
        data=0.d0
        data2=0.d0
        data3=0.d0
        data4=0.d0
        data5=0.d0
        npw_not_zero=0
        if(lreal)then
         data(1)=dble(evc(1))/3.375d0
         do ig=2,npwx
          if(use_this_g(ig))then
           npw_not_zero=npw_not_zero+1
           icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
            &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
            &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
           data(icount) =                 (evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig) 
           data2(icount)=         -g2(ig)*(evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig) 
           data3(icount)=imag*gcart(1,ig)*(evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig) 
           data4(icount)=imag*gcart(2,ig)*(evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig) 
           data5(icount)=imag*gcart(3,ig)*(evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig) 
          endif
         enddo
! avc = \sum_g gamma(g) * Re[ evc(g) * exp( i.g.r ) ]
! If the wave function is real then evc(g) = conjg(evc(-g)) and the sum
! can be done only on +g because gamma(g) = gamma(-g).
! If the WF is not real but we want to take the real part then we take
! evc(ig) + conjg( evc(-g) ).
!      write(6,*)'Taking the REAL part of WF # ',ibnd,'; # of PW coefficients: ',&
!       &npw_not_zero
        elseif(limag)then
         data(1)=cmplx(0.d0,aimag(evc(1))/3.375d0)
         do ig=2,npwx
          if(use_this_g(ig))then
           npw_not_zero=npw_not_zero+1
           icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
            (nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
            (nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
           data(icount) =                 (evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig) 
           data2(icount)=         -g2(ig)*(evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig) 
           data3(icount)=imag*gcart(1,ig)*(evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig) 
           data4(icount)=imag*gcart(2,ig)*(evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig) 
           data5(icount)=imag*gcart(3,ig)*(evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig) 
          endif
         enddo
!    If the WF is not imaginary but we want to take the imaginary part then we take
!    evc(ig)-conjg(evc(-g)) 
!      write(6,*)'Taking the IMAGINARY part of WF # ',ibnd,'; # of PW &
!      &coefficients: ',npw_not_zero
        else
         stop 'Wavefunction is not real nor imaginary, what should I do?'
        endif
     
        call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
        call fftn(data2,(/nr(1),nr(2),nr(3)/),inv=.true.)
        call fftn(data3,(/nr(1),nr(2),nr(3)/),inv=.true.)
        call fftn(data4,(/nr(1),nr(2),nr(3)/),inv=.true.)
        call fftn(data5,(/nr(1),nr(2),nr(3)/),inv=.true.)
        data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
        data2=data2*sqrt(dble(nr(1)*nr(2)*nr(3)))
        data3=data3*sqrt(dble(nr(1)*nr(2)*nr(3)))
        data4=data4*sqrt(dble(nr(1)*nr(2)*nr(3)))
        data5=data5*sqrt(dble(nr(1)*nr(2)*nr(3)))
        icount=1
        do lz=0,nr(3)-1
         do ly=0,nr(2)-1
          do lx=0,nr(1)-1
           if(lreal)then
            avc(lx,ly,lz)=data(icount)
            avc2(lx,ly,lz)=data2(icount)
            avc3(lx,ly,lz)=data3(icount)
            avc4(lx,ly,lz)=data4(icount)
            avc5(lx,ly,lz)=data5(icount)
           else
            avc(lx,ly,lz)=aimag(data(icount))
            avc2(lx,ly,lz)=aimag(data2(icount))
            avc3(lx,ly,lz)=aimag(data3(icount))
            avc4(lx,ly,lz)=aimag(data4(icount))
            avc5(lx,ly,lz)=aimag(data5(icount))
           endif
           icount=icount+1
          enddo
         enddo
        enddo

        write(io,'(a)')'Orbital'
        do lx=0,nr(1)-1
         do ly=0,nr(2)-1
          do lz=0,nr(3)-1 
           write(io,*)avc(lx,ly,lz)
          enddo
         enddo
        enddo
        write(io,'(a)')'Laplacian'
        do lx=0,nr(1)-1
         do ly=0,nr(2)-1
          do lz=0,nr(3)-1 
           write(io,*)avc2(lx,ly,lz)
          enddo
         enddo
        enddo
        write(io,'(a)')'Gradient - a1 direction'
        do lx=0,nr(1)-1
         do ly=0,nr(2)-1
          do lz=0,nr(3)-1 
           write(io,*)avc3(lx,ly,lz)
          enddo
         enddo
        enddo
        write(io,'(a)')'Gradient - a2 direction'
        do lx=0,nr(1)-1
         do ly=0,nr(2)-1
          do lz=0,nr(3)-1 
           write(io,*)avc4(lx,ly,lz)
          enddo
         enddo
        enddo
        write(io,'(a)')'Gradient - a3 direction'
        do lx=0,nr(1)-1
         do ly=0,nr(2)-1
          do lz=0,nr(3)-1 
           write(io,*)avc5(lx,ly,lz)
          enddo
         enddo
        enddo
           
       if(ltest)then
        xb=0 ; xd=0 ; xbd=0
        do i=1,n_points_for_test 
         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)
         call blip3d(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d2=0
         dg=0
         if(lreal)then
          d=dble(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)+conjg(evc(minus_g(ig))))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+dble(eigr)
            d2=d2-dble(eigr)*g2(ig)
            dg(1)=dg(1)-tpi*g(1,ig)*aimag(eigr)
            dg(2)=dg(2)-tpi*g(2,ig)*aimag(eigr)
            dg(3)=dg(3)-tpi*g(3,ig)*aimag(eigr)
           endif
          enddo
         else
          d=aimag(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)-conjg(evc(minus_g(ig))))* &
            &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+aimag(eigr)
            d2=d2-aimag(eigr)*g2(ig) 
            dg(1)=dg(1)+tpi*g(1,ig)*dble(eigr)
            dg(2)=dg(2)+tpi*g(2,ig)*dble(eigr)
            dg(3)=dg(3)+tpi*g(3,ig)*dble(eigr)
           endif
          enddo
         endif
         dg=matmul(bg,dg)
         xb(1)=xb(1)+blip(1)**2
         xd(1)=xd(1)+d**2
         xbd(1)=xbd(1)+blip(1)*d
         xb(2)=xb(2)+blip(2)**2
         xd(2)=xd(2)+d2**2
         xbd(2)=xbd(2)+blip(2)*d2
         xb(3)=xb(3)+blip(3)**2
         xd(3)=xd(3)+dg(1)**2
         xbd(3)=xbd(3)+blip(3)*dg(1)
         xb(4)=xb(4)+blip(4)**2
         xd(4)=xd(4)+dg(2)**2
         xbd(4)=xbd(4)+blip(4)*dg(2)
         xb(5)=xb(5)+blip(5)**2
         xd(5)=xd(5)+dg(3)**2
         xbd(5)=xbd(5)+blip(5)*dg(3)
        enddo
       endif
       if(ltest)then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2,(xbd(3)/sqrt(xb(3)*xd(3)))**2,     & 
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2,(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif

     end select
   enddo
   if(ltest)write(6,*)

  enddo
 enddo
 
 END SUBROUTINE blipp


 SUBROUTINE blipr
  USE globals
  USE singleton
  IMPLICIT NONE
  INTEGER i,j,l,nr(3),ibnd,ig,ig1,npw_not_zero,npwx,lx,ly,lz,nks,ik,nat,na, &
   &atomno,nbnd(2),ispin,nspin,nelec,isign,icount,indice,ierr
  INTEGER :: io=10,io2=11
  INTEGER,PARAMETER :: n_points_for_test=1000
  INTEGER,ALLOCATABLE :: minus_g(:)
  DOUBLE PRECISION :: d,d2,dg(3),x,y,z,k,blip(5),gxmax,gymax,gzmax, &
   &dax,day,daz,r(3),xb(5),xd(5),xbd(5),at(3,3),bg(3,3),tv,tc,tv0,tc0,tmp, &
   &ranfx,tau(3),etot,ewld,ecut,xk(3),et,arg,norm_real,norm_imag
  DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0,tpi=2.d0*pi,eps=1.d-6
  DOUBLE PRECISION,ALLOCATABLE :: g(:,:),gcart(:,:),g2(:),avc(:,:,:),evc(:),gamma(:)
  DOUBLE PRECISION,ALLOCATABLE :: avc2(:,:,:) !coefficients of Laplacian
  DOUBLE PRECISION,ALLOCATABLE :: avc3(:,:,:),avc4(:,:,:),avc5(:,:,:) !gradient
  COMPLEX(KIND=KIND(0.d0)) :: eigr
  COMPLEX(KIND=KIND(0.d0)),PARAMETER :: imag=(0.d0,1.d0)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: data(:)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: data2(:)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: data3(:),data4(:),data5(:)
  LOGICAL lsda,lreal,limag
  LOGICAL,ALLOCATABLE :: use_this_g(:)
  CHARACTER(80) title

  write(6,*)'Real PW coefficients; GAMMA point'
  write(6,*)

  open(io2,file='pwfn.data',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening pwfn.data. Stopping.'
   stop
  endif
  open(io,file='bwfn.data',status='unknown',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening bwfn.data. Stopping.'
   stop
  endif

  call bwfn_info(io,io2,lsda)

  call skipio(io2,1) 
  write(io,'(a)')'Primitive lattice vectors (au) '
  read(io2,*)at(:,1)
  read(io2,*)at(:,2)
  read(io2,*)at(:,3)
  write(io,'(3f20.12)')at(:,1)
  write(io,'(3f20.12)')at(:,2)
  write(io,'(3f20.12)')at(:,3)

  call inve(at,bg)
  bg=transpose(bg)

  call skipio(io2,4) 
  write(io,'(a)')' '
  write(io,'(a)')'G VECTORS'
  write(io,'(a)')'---------'
  write(io,'(a)')'Number of G-vectors'
  read(io2,*) npwx
  write(io,*) npwx
  call skipio(io2,1) 
  write(io,'(a)')'Gx Gy Gz (au)'

  allocate(g(3,npwx),gcart(3,npwx),g2(npwx))
  gxmax=0 ; gymax=0 ; gzmax=0
  do ig=1,npwx
   read(io2,*)g(:,ig)
   write(io,'(3f20.12)')g(:,ig)
   gcart(:,ig)=g(:,ig)
   g2(ig)=dot_product(g(:,ig),g(:,ig))
   g(:,ig)=matmul(transpose(at),g(:,ig))/tpi
   if(abs(g(1,ig))>gxmax)gxmax=abs(g(1,ig))
   if(abs(g(2,ig))>gymax)gymax=abs(g(2,ig))
   if(abs(g(3,ig))>gzmax)gzmax=abs(g(3,ig))
  enddo
!  write(6,'('' Gmax:'',3f20.10)')gxmax,gymax,gzmax
  nr(1)=2*nint(gxmax)+2
  nr(2)=2*nint(gymax)+2
  nr(3)=2*nint(gzmax)+2
  nr=nr*xmul
! The following is to make nr an even number
  if(mod(nr(1),2)/=0)nr(1)=nr(1)+1
  if(mod(nr(2),2)/=0)nr(2)=nr(2)+1
  if(mod(nr(3),2)/=0)nr(3)=nr(3)+1
  write(6,'('' Grid : '',3i5)')nr
  write(6,*)
  dax=1.d0/dble(nr(1))
  day=1.d0/dble(nr(2))
  daz=1.d0/nr(3)
  write(io,'(a)')'Blip grid'
  write(io,'(3i4)')nr

  call skipio(io2,4) 
  write(io,'(a)')' '
  write(io,'(a)')'WAVE FUNCTION'
  write(io,'(a)')'-------------'
  write(io,'(a)')'Number of k-points'
  read(io2,*)nks
  if(nks>1)then
   write(6,*)' !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
   write(6,*)' !!! WARNING, K-points scheme not yet implemented !!!' 
   write(6,*)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
  endif
  write(io,*)nks
  if(lsda)then
   nspin=2
  else
   nspin=1
  endif

  allocate(evc(npwx),use_this_g(npwx),minus_g(npwx),gamma(npwx))
  select case (igrad_lap)
  case (0)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate( data(nr(1)*nr(2)*nr(3)) )
  case (1)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc2(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate( data(nr(1)*nr(2)*nr(3)) )
    allocate( data2(nr(1)*nr(2)*nr(3)) )
  case (2)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc2(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc3(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc4(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc5(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate( data(nr(1)*nr(2)*nr(3)) )
    allocate( data2(nr(1)*nr(2)*nr(3)) )
    allocate( data3(nr(1)*nr(2)*nr(3)) )
    allocate( data4(nr(1)*nr(2)*nr(3)) )
    allocate( data5(nr(1)*nr(2)*nr(3)) )
  end select

! Calculating gamma
  do ig=1,npwx
   k=dax*g(1,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=1.5d0
   endif
   k=day*g(2,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=gamma(ig)*3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=gamma(ig)*1.5d0
   endif
   k=daz*g(3,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=gamma(ig)*3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=gamma(ig)*1.5d0
   endif
   gamma(ig)=1.d0/gamma(ig)
  enddo

! Set up to use only +G
  use_this_g=.true.
  minus_g=0
  minus_g(1)=1
  i=0

  if(abs(g(1,1))>eps.or.abs(g(1,1))>eps.or.abs(g(1,1))>eps) &
   stop 'first vector has to be the zero vector, reorder vectors in &
   &pwfn.data please'

  do ig=1,npwx 
   if(use_this_g(ig))then
! Set exp( i g.r )
! Find -G
    do_ig: do ig1=ig+1,npwx
     if(abs(g(1,ig)+g(1,ig1))<eps.and.abs(g(2,ig)+g(2,ig1))<eps.and. &
      &abs(g(3,ig)+g(3,ig1))<eps)then
      use_this_g(ig1)=.false.
      minus_g(ig)=ig1
      i=i+1
      exit do_ig
     endif
    enddo do_ig
   endif
  enddo
!  if(2*i+1/=npwx)stop 'I could not find all +G -G vectors'

  do ik=1,nks

   call skipio(io2,1)
   write(io,'(a)')'k-point # ; # of bands (up spin/down spin) ; k-point &
    &coords (au)'
   read(io2,*)i,nbnd(1),nbnd(2),xk(:)

   if(abs(xk(1))>eps.or.abs(xk(2))>eps.or.abs(xk(3))>eps)then
    write(6,*)' !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
    write(6,*)' !!! WARNING, non gamma scheme, we should call blipk !!!' 
    write(6,*)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
   endif

   write(io,'(3i4,3f20.16)')ik,nbnd(1),nbnd(2),xk(:)

   do ispin=1,nspin 

    if(ltest)then
     if(lsda)write(6,*)'Spin : ',ispin
     write(6,*)'Band ; alpha ; alpha (Laplacian) ; alpha (Gradient)'
    endif

    do ibnd=1,nbnd(ispin)

     call skipio(io2,1)
     write(io,'(a)')'Band, spin, eigenvalue (au)'
     read(io2,*)i,i,et 
     write(io,'(2i5,f20.12)')ibnd,ispin,et
     call skipio(io2,1)

! Read in the wavefunctions
     do ig=1,npwx
      read(io2,*) evc(ig)
     enddo

! Decide if take the real or the imaginary part
     norm_real=evc(1)**2
     norm_imag=0
     do ig=2,npwx
      if(use_this_g(ig))then
       norm_real=norm_real+2*(abs(evc(ig)+evc(minus_g(ig)))/2)**2
       norm_imag=norm_imag+2*(abs(evc(ig)-evc(minus_g(ig)))/2)**2
      endif
     enddo

     lreal=.false. ; limag=.false.
     if(abs(norm_real-1.d0)<eps)lreal=.true.
     if(abs(norm_imag-1.d0)<eps)limag=.true.
     if(lreal.and.limag)stop 'WF both pure real and pure imaginary ???'

     if(lreal)write(6,*)'Real wave function, band : ',ibnd
     if(limag)write(6,*)'Imaginary wave function, band : ',ibnd

     if(.not.lreal.and..not.limag)then 
      if(norm_real>norm_imag)then
       lreal=.true.
      else
       limag=.true.
      endif
     endif

     write(io,'(a)')'Blip coefficients'
     select case (igrad_lap)
     case (0)
       avc=0
       npw_not_zero=0
       
       data=0.d0
       if(lreal)then
        data(1)=dble(evc(1))/3.375d0
        do ig=2,npwx
         if(use_this_g(ig))then
          npw_not_zero=npw_not_zero+1
          icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
           &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
           &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
          data(icount)=(evc(ig)+evc(minus_g(ig)))*gamma(ig)              
         endif
        enddo
! avc = \sum_g gamma(g) * Re[ evc(g) * exp( i.g.r ) ]
! if the wavefunction is real then evc(g) = conjg(evc(-g)) and the sum
! can be done only on +g because gamma(g) = gamma(-g).
! if the WF is not real but we want to take the real part then we take
! evc(ig) + conjg( evc(-g) ) .....
        write(6,*)'Taking the real part of WF # ',ibnd,'; # of PW coefficients: ',&
         &npw_not_zero
       elseif(limag)then
        do ig=2,npwx
         if(use_this_g(ig))then
          npw_not_zero=npw_not_zero+1
          icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
           &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
           &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
          data(icount)=(evc(ig)-evc(minus_g(ig)))*gamma(ig)              
         endif
        enddo
        write(6,*)'Taking the imaginary part of WF # ',ibnd,'; # of PW &
         &coefficients: ',npw_not_zero
       else
        stop 'Wavefunction is not real nor imaginary, what should I do ???'
       endif

       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          if(lreal)then
           avc(lx,ly,lz)=data(icount)
          else
           avc(lx,ly,lz)=aimag(data(icount))
          endif
          icount=icount+1
         enddo
        enddo
       enddo

       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc(lx,ly,lz)
         enddo
        enddo
       enddo
                         
       if(ltest)then

        xb=0 ; xd=0 ; xbd=0

        do i=1,n_points_for_test 

         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)
         call blip3d(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d2=0
         dg=0
         d=0

         if(lreal)then

          d=dble(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)+evc(minus_g(ig)))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+dble(eigr)
            d2=d2-dble(eigr)*g2(ig)
            dg(1)=dg(1)-tpi*g(1,ig)*aimag(eigr)
            dg(2)=dg(2)-tpi*g(2,ig)*aimag(eigr)
            dg(3)=dg(3)-tpi*g(3,ig)*aimag(eigr)
           endif
          enddo

         else

          do ig=1,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)-evc(minus_g(ig)))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+aimag(eigr)
            d2=d2-aimag(eigr)*g2(ig) 
            dg(1)=dg(1)+tpi*g(1,ig)*dble(eigr)
            dg(2)=dg(2)+tpi*g(2,ig)*dble(eigr)
            dg(3)=dg(3)+tpi*g(3,ig)*dble(eigr)
           endif
          enddo

         endif

         dg=matmul(bg,dg)
         xb(1)=xb(1)+blip(1)**2
         xd(1)=xd(1)+d**2
         xbd(1)=xbd(1)+blip(1)*d
         xb(2)=xb(2)+blip(2)**2
         xd(2)=xd(2)+d2**2
         xbd(2)=xbd(2)+blip(2)*d2
         xb(3)=xb(3)+blip(3)**2
         xd(3)=xd(3)+dg(1)**2
         xbd(3)=xbd(3)+blip(3)*dg(1)
         xb(4)=xb(4)+blip(4)**2
         xd(4)=xd(4)+dg(2)**2
         xbd(4)=xbd(4)+blip(4)*dg(2)
         xb(5)=xb(5)+blip(5)**2
         xd(5)=xd(5)+dg(3)**2
         xbd(5)=xbd(5)+blip(5)*dg(3)

        enddo

       endif
   
       if(ltest) then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2,(xbd(3)/sqrt(xb(3)*xd(3)))**2,     &
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2,(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif
     case (1)
       avc=0
       avc2=0
       npw_not_zero=0
       
       data=0.d0
       data2=0.d0
       if(lreal)then
        data(1)=dble(evc(1))/3.375d0
        do ig=2,npwx
         if(use_this_g(ig))then
          npw_not_zero=npw_not_zero+1
          icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
           &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
           &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
          data(icount)=(evc(ig)+evc(minus_g(ig)))*gamma(ig)              
          data2(icount)=-g2(ig)*(evc(ig)+evc(minus_g(ig)))*gamma(ig)              
         endif
        enddo
! avc = \sum_g gamma(g) * Re[ evc(g) * exp( i.g.r ) ]
! if the wavefunction is real then evc(g) = conjg(evc(-g)) and the sum
! can be done only on +g because gamma(g) = gamma(-g).
! if the WF is not real but we want to take the real part then we take
! evc(ig) + conjg( evc(-g) ) .....
        write(6,*)'Taking the real part of WF # ',ibnd,'; # of PW coefficients: ',&
         &npw_not_zero
       elseif(limag)then
        do ig=2,npwx
         if(use_this_g(ig))then
          npw_not_zero=npw_not_zero+1
          icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
           &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
           &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
          data(icount)=(evc(ig)-evc(minus_g(ig)))*gamma(ig)              
          data2(icount)=-g2(ig)*(evc(ig)-evc(minus_g(ig)))*gamma(ig)              
         endif
        enddo
        write(6,*)'Taking the imaginary part of WF # ',ibnd,'; # of PW &
         &coefficients: ',npw_not_zero
       else
        stop 'Wavefunction is not real nor imaginary, what should I do ???'
       endif

       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data2,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data2=data2*sqrt(dble(nr(1)*nr(2)*nr(3)))
       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          if(lreal)then
           avc(lx,ly,lz)=data(icount)
           avc2(lx,ly,lz)=data2(icount)
          else
           avc(lx,ly,lz)=aimag(data(icount))
           avc2(lx,ly,lz)=aimag(data2(icount))
          endif
          icount=icount+1
         enddo
        enddo
       enddo

       write(io,'(a)')'Orbital'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Laplacian'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc2(lx,ly,lz)
         enddo
        enddo
       enddo
                         
       if(ltest)then

        xb=0 ; xd=0 ; xbd=0

        do i=1,n_points_for_test 

         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)
         call blip3d(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d2=0
         dg=0
         d=0

         if(lreal)then

          d=dble(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)+evc(minus_g(ig)))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+dble(eigr)
            d2=d2-dble(eigr)*g2(ig)
            dg(1)=dg(1)-tpi*g(1,ig)*aimag(eigr)
            dg(2)=dg(2)-tpi*g(2,ig)*aimag(eigr)
            dg(3)=dg(3)-tpi*g(3,ig)*aimag(eigr)
           endif
          enddo

         else

          do ig=1,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)-evc(minus_g(ig)))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+aimag(eigr)
            d2=d2-aimag(eigr)*g2(ig) 
            dg(1)=dg(1)+tpi*g(1,ig)*dble(eigr)
            dg(2)=dg(2)+tpi*g(2,ig)*dble(eigr)
            dg(3)=dg(3)+tpi*g(3,ig)*dble(eigr)
           endif
          enddo

         endif

         dg=matmul(bg,dg)
         xb(1)=xb(1)+blip(1)**2
         xd(1)=xd(1)+d**2
         xbd(1)=xbd(1)+blip(1)*d
         xb(2)=xb(2)+blip(2)**2
         xd(2)=xd(2)+d2**2
         xbd(2)=xbd(2)+blip(2)*d2
         xb(3)=xb(3)+blip(3)**2
         xd(3)=xd(3)+dg(1)**2
         xbd(3)=xbd(3)+blip(3)*dg(1)
         xb(4)=xb(4)+blip(4)**2
         xd(4)=xd(4)+dg(2)**2
         xbd(4)=xbd(4)+blip(4)*dg(2)
         xb(5)=xb(5)+blip(5)**2
         xd(5)=xd(5)+dg(3)**2
         xbd(5)=xbd(5)+blip(5)*dg(3)

        enddo

       endif
   
       if(ltest) then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2,(xbd(3)/sqrt(xb(3)*xd(3)))**2,     &
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2,(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif
     case (2)
       avc=0
       avc2=0
       avc3=0
       avc4=0
       avc5=0
       npw_not_zero=0
       
       data=0.d0
       data2=0.d0
       data3=0.d0
       data4=0.d0
       data5=0.d0
       if(lreal)then
        data(1)=dble(evc(1))/3.375d0
        do ig=2,npwx
         if(use_this_g(ig))then
          npw_not_zero=npw_not_zero+1
          icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
           &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
           &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
          data(icount) =                 (evc(ig)+evc(minus_g(ig)))*gamma(ig)              
          data2(icount)=         -g2(ig)*(evc(ig)+evc(minus_g(ig)))*gamma(ig)              
          data3(icount)=imag*gcart(1,ig)*(evc(ig)+evc(minus_g(ig)))*gamma(ig)              
          data4(icount)=imag*gcart(2,ig)*(evc(ig)+evc(minus_g(ig)))*gamma(ig)              
          data5(icount)=imag*gcart(3,ig)*(evc(ig)+evc(minus_g(ig)))*gamma(ig)              
         endif
        enddo
! avc = \sum_g gamma(g) * Re[ evc(g) * exp( i.g.r ) ]
! if the wavefunction is real then evc(g) = conjg(evc(-g)) and the sum
! can be done only on +g because gamma(g) = gamma(-g).
! if the WF is not real but we want to take the real part then we take
! evc(ig) + conjg( evc(-g) ) .....
        write(6,*)'Taking the real part of WF # ',ibnd,'; # of PW coefficients: ',&
         &npw_not_zero
       elseif(limag)then
        do ig=2,npwx
         if(use_this_g(ig))then
          npw_not_zero=npw_not_zero+1
          icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
           &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
           &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
          data(icount)=                  (evc(ig)-evc(minus_g(ig)))*gamma(ig)              
          data2(icount)=         -g2(ig)*(evc(ig)-evc(minus_g(ig)))*gamma(ig)              
          data3(icount)=imag*gcart(1,ig)*(evc(ig)-evc(minus_g(ig)))*gamma(ig)              
          data4(icount)=imag*gcart(2,ig)*(evc(ig)-evc(minus_g(ig)))*gamma(ig)              
          data5(icount)=imag*gcart(3,ig)*(evc(ig)-evc(minus_g(ig)))*gamma(ig)              
         endif
        enddo
        write(6,*)'Taking the imaginary part of WF # ',ibnd,'; # of PW &
         &coefficients: ',npw_not_zero
       else
        stop 'Wavefunction is not real nor imaginary, what should I do ???'
       endif

       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data2,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data3,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data4,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data5,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data2=data2*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data3=data3*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data4=data4*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data5=data5*sqrt(dble(nr(1)*nr(2)*nr(3)))
       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          if(lreal)then
           avc(lx,ly,lz)=data(icount)
           avc2(lx,ly,lz)=data2(icount)
           avc3(lx,ly,lz)=data3(icount)
           avc4(lx,ly,lz)=data4(icount)
           avc5(lx,ly,lz)=data5(icount)
          else
           avc(lx,ly,lz)=aimag(data(icount))
           avc2(lx,ly,lz)=aimag(data2(icount))
           avc3(lx,ly,lz)=aimag(data3(icount))
           avc4(lx,ly,lz)=aimag(data4(icount))
           avc5(lx,ly,lz)=aimag(data5(icount))
          endif
          icount=icount+1
         enddo
        enddo
       enddo

       write(io,'(a)')'Orbital'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Laplacian'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc2(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Gradient - a1 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc3(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Gradient - a2 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc4(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Gradient - a3 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc5(lx,ly,lz)
         enddo
        enddo
       enddo
                         
       if(ltest)then

        xb=0 ; xd=0 ; xbd=0

        do i=1,n_points_for_test 

         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)
         call blip3d(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d2=0
         dg=0
         d=0

         if(lreal)then

          d=dble(evc(1))
          do ig=2,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)+evc(minus_g(ig)))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+dble(eigr)
            d2=d2-dble(eigr)*g2(ig)
            dg(1)=dg(1)-tpi*g(1,ig)*aimag(eigr)
            dg(2)=dg(2)-tpi*g(2,ig)*aimag(eigr)
            dg(3)=dg(3)-tpi*g(3,ig)*aimag(eigr)
           endif
          enddo

         else

          do ig=1,npwx
           if(use_this_g(ig))then
            eigr=(evc(ig)-evc(minus_g(ig)))* &
             &exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
            d=d+aimag(eigr)
            d2=d2-aimag(eigr)*g2(ig) 
            dg(1)=dg(1)+tpi*g(1,ig)*dble(eigr)
            dg(2)=dg(2)+tpi*g(2,ig)*dble(eigr)
            dg(3)=dg(3)+tpi*g(3,ig)*dble(eigr)
           endif
          enddo

         endif

         dg=matmul(bg,dg)
         xb(1)=xb(1)+blip(1)**2
         xd(1)=xd(1)+d**2
         xbd(1)=xbd(1)+blip(1)*d
         xb(2)=xb(2)+blip(2)**2
         xd(2)=xd(2)+d2**2
         xbd(2)=xbd(2)+blip(2)*d2
         xb(3)=xb(3)+blip(3)**2
         xd(3)=xd(3)+dg(1)**2
         xbd(3)=xbd(3)+blip(3)*dg(1)
         xb(4)=xb(4)+blip(4)**2
         xd(4)=xd(4)+dg(2)**2
         xbd(4)=xbd(4)+blip(4)*dg(2)
         xb(5)=xb(5)+blip(5)**2
         xd(5)=xd(5)+dg(3)**2
         xbd(5)=xbd(5)+blip(5)*dg(3)

        enddo

       endif
   
       if(ltest) then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2,(xbd(3)/sqrt(xb(3)*xd(3)))**2,     &
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2,(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif
     end select
 
    enddo
    if(ltest)write(6,*)

   enddo
  enddo
 
 END SUBROUTINE blipr


 SUBROUTINE blipk
  USE globals
  USE singleton
  IMPLICIT NONE
  integer, parameter :: n_points_for_test = 100
  INTEGER i,j,l,nr(3),ibnd,ig,ig1,npwx,lx,ly,lz,nks,ik,nat,na,atomno,nbnd(2), &
   &ispin,nspin,nelec,isign,icount,indice,ierr
  INTEGER :: io=10,io2=11
  DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0,tpi=2.d0*pi,eps=1.d-8
  DOUBLE PRECISION x,y,z,k,alpha,gxmax,gymax,gzmax,alat,dax,day,daz, &
   &r(3),da(3),xb(5),xd(5),xbd(5),at(3,3),bg(3,3),tv,tc,tv0,tc0,tmp,tau(3),&
   &etot,ewld,ecut,xk(3),et,arg,ranfx
  DOUBLE PRECISION,ALLOCATABLE :: g(:,:),gcart(:,:),g2(:)
  COMPLEX(KIND=KIND(0.d0)) :: d,d2,dg(3),blip(5),eigr,psi
  COMPLEX(KIND=KIND(0.d0)),PARAMETER :: imag=(0.d0,1.d0)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: evc(:),gamma(:),avc(:,:,:),data(:)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: avc2(:,:,:),data2(:) !Laplacian
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: avc3(:,:,:),data3(:),avc4(:,:,:),&
   &data4(:),avc5(:,:,:),data5(:) !gradient
  LOGICAL lsda,lreal,limag
  CHARACTER(80) title

  write(6,*)'Complex PW coefficients; K-POINTS'
  write(6,*)

  open(io2,file='pwfn.data',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Cannot open pwfn.data.  Stopping.'
   stop
  endif
  open(io,file='bwfn.data',status='unknown',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Cannot open bwfn.data.  Stopping.'
   stop
  endif

  call bwfn_info(io,io2,lsda)

  call skipio(io2,1) 
  write(io,'(a)')'Primitive lattice vectors (au) '
  read(io2,*)at(:,1)
  read(io2,*)at(:,2)
  read(io2,*)at(:,3)
  write(io,'(3f20.12)')at(:,1)
  write(io,'(3f20.12)')at(:,2)
  write(io,'(3f20.12)')at(:,3)

  call inve(at,bg)
  bg=transpose(bg)

  call skipio(io2,4) 
  write(io,'(a)')' '
  write(io,'(a)')'G VECTORS'
  write(io,'(a)')'---------'
  write(io,'(a)')'Number of G-vectors'
  read(io2,*)npwx
  write(io,*)npwx
  call skipio(io2,1) 
  write(io,'(a)')'Gx Gy Gz (au)'

  allocate(g(3,npwx),gcart(3,npwx),g2(npwx))
  gxmax=0 ; gymax=0 ; gzmax=0

  do ig=1,npwx
   read(io2,*)g(:,ig)
   write(io,'(3f20.12)')g(:,ig)
   g2(ig)=dot_product(g(:,ig),g(:,ig))
   gcart(:,ig)=g(:,ig)
   g(:,ig)=matmul(transpose(at),g(:,ig))/tpi
   if(abs(g(1,ig))>gxmax)gxmax=abs(g(1,ig))
   if(abs(g(2,ig))>gymax)gymax=abs(g(2,ig))
   if(abs(g(3,ig))>gzmax)gzmax=abs(g(3,ig))
  enddo

!  write(6,'('' Gmax:'',3f20.10)')gxmax,gymax,gzmax
  nr(1)=2*nint(gxmax)+2
  nr(2)=2*nint(gymax)+2
  nr(3)=2*nint(gzmax)+2
  nr=nr*xmul
! The following is to make nr an even number
  if(mod(nr(1),2)/=0)nr(1)=nr(1)+1
  if(mod(nr(2),2)/=0)nr(2)=nr(2)+1
  if(mod(nr(3),2)/=0)nr(3)=nr(3)+1
  write(6,'('' Grid : '',3i5)')nr
  write(6,*)
  dax=1.d0/dble(nr(1))
  day=1.d0/dble(nr(2))
  daz=1.d0/dble(nr(3))
  write(io,'(a)')'Blip grid'
  write(io,'(3i4)')nr

  call skipio(io2,4) 
  write(io,'(a)')' '
  write(io,'(a)')'WAVE FUNCTION'
  write(io,'(a)')'-------------'
  write(io,'(a)')'Number of k-points'
  read(io2,*)nks
  write(io,*)nks
  if(lsda)then
   nspin=2
  else
   nspin=1
  endif

  select case (igrad_lap)
  case(0)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(data(nr(1)*nr(2)*nr(3)))
  case(1)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc2(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(data(nr(1)*nr(2)*nr(3)))
    allocate(data2(nr(1)*nr(2)*nr(3)))
  case(2)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc2(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc3(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc4(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc5(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(data(nr(1)*nr(2)*nr(3)))
    allocate(data2(nr(1)*nr(2)*nr(3)))
    allocate(data3(nr(1)*nr(2)*nr(3)))
    allocate(data4(nr(1)*nr(2)*nr(3)))
    allocate(data5(nr(1)*nr(2)*nr(3)))
  end select
  allocate(evc(npwx),gamma(npwx))

! Calculating gamma
  do ig=1,npwx
   k=dax*g(1,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=1.5d0
   endif
   k=day*g(2,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=gamma(ig)*3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=gamma(ig)*1.5d0
   endif
   k=daz*g(3,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=gamma(ig)*3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=gamma(ig)*1.5d0
   endif
   gamma(ig)=1.d0/gamma(ig)
  enddo

  if(abs(g(1,1))>eps.or.abs(g(1,1))>eps.or.abs(g(1,1))>eps) &
   stop 'First vector has to be the zero vector, reorder vectors &
    &in pwfn.data please'

  do ik=1,nks

   call skipio(io2,1)
   write(io,'(a)')'k-point # ; # of bands (up spin/down spin) ; k-point &
    &coords (au)'
   read(io2,*)i,nbnd(1),nbnd(2),xk(:)

   write(io,'(3i4,3f20.16)')ik,nbnd(1),nbnd(2),xk(:)

   do ispin=1,nspin 

!    xk(:)=matmul( transpose(at),xk(:))/tpi
    write(6,'('' k-point '',i5,"  :  (",f10.6,",",f10.6,",",f10.6,")")')ik,xk

    if(ltest)then
     if(lsda)write(6,*)'Spin : ',ispin
     write(6,*)'Band ; alpha ; alpha (Laplacian) ; alpha (Gradient)'
    endif

    do ibnd=1,nbnd(ispin)
     call skipio(io2,1)
     write(io,'(a)')'Band, spin, eigenvalue (au)'
     read(io2,*)i,i,et 
     write(io,'(2i5,f20.12)')ibnd,ispin,et
     call skipio(io2,1)
     
     select case (igrad_lap)
     case(0)
! Read in the wavefunctions
       data = 0.d0
       do ig=1,npwx
        read(io2,*)evc(ig)
        icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
         &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
         &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
        data(icount)=evc(ig)*gamma(ig)              
       enddo
       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))

       write(io,'(a)') 'Blip coefficients'
       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          avc(lx,ly,lz)=data(icount)
          icount=icount+1
         enddo
        enddo
       enddo

       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*) avc(lx,ly,lz)
         enddo
        enddo
       enddo
             
       if(ltest)then

        xb=0 ; xd=0 ; xbd =0
        do i=1,n_points_for_test 
         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)

         call blip3dk(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d=0
         d2=0
         dg=0

         do ig=1,npwx
           eigr=evc(ig)*exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
           d=d+eigr
           d2=d2-eigr*g2(ig)
           dg(1)=dg(1)-tpi*g(1,ig)*eigr*imag
           dg(2)=dg(2)-tpi*g(2,ig)*eigr*imag
           dg(3)=dg(3)-tpi*g(3,ig)*eigr*imag
         enddo

         dg=matmul(bg,dg)
         xb(1)=xb(1)+dble(blip(1))**2
         xd(1)=xd(1)+dble(d)**2
         xbd(1)=xbd(1)+dble(blip(1))*dble(d)
         xb(2)=xb(2)+dble(blip(2))**2
         xd(2)=xd(2)+dble(d2)**2
         xbd(2)=xbd(2)+dble(blip(2))*dble(d2)
         xb(3)=xb(3)+dble(blip(3))**2
         xd(3)=xd(3)+dble(dg(1))**2
         xbd(3)=xbd(3)+dble(blip(3))*dble(dg(1))
         xb(4)=xb(4)+dble(blip(4))**2
         xd(4)=xd(4)+dble(dg(2))**2
         xbd(4)=xbd(4)+dble(blip(4))*dble(dg(2))
         xb(5)=xb(5)+dble(blip(5))**2
         xd(5)=xd(5)+dble(dg(3))**2
         xbd(5)=xbd(5)+dble(blip(5))*dble(dg(3))

        enddo

       endif

       if(ltest) then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2) / sqrt(xb(2)*xd(2)))**2, &
         &(xbd(3) / sqrt(xb(3)*xd(3)))**2, (xbd(4) / sqrt(xb(4)*xd(4)))**2,      &
         &(xbd(5) / sqrt(xb(5)*xd(5)))**2
       endif
     case(1)
! Read in the wavefunctions
       data = 0.d0
       data2 = 0.d0
       do ig=1,npwx
        read(io2,*)evc(ig)
        icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
         &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
         &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
        data(icount)=evc(ig)*gamma(ig)              
        data2(icount)=-g2(ig)*evc(ig)*gamma(ig)              
       enddo
       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data2,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data2=data2*sqrt(dble(nr(1)*nr(2)*nr(3)))

       write(io,'(a)') 'Blip coefficients'
       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          avc(lx,ly,lz)=data(icount)
          avc2(lx,ly,lz)=data2(icount)
          icount=icount+1
         enddo
        enddo
       enddo

       write(io,'(a)') 'Orbital'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*) avc(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)') 'Laplacian'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*) avc2(lx,ly,lz)
         enddo
        enddo
       enddo
             
       if(ltest)then

        xb=0 ; xd=0 ; xbd =0
        do i=1,n_points_for_test 
         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)

         call blip3dk(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d=0
         d2=0
         dg=0

         do ig=1,npwx
           eigr=evc(ig)*exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
           d=d+eigr
           d2=d2-eigr*g2(ig)
           dg(1)=dg(1)-tpi*g(1,ig)*eigr*imag
           dg(2)=dg(2)-tpi*g(2,ig)*eigr*imag
           dg(3)=dg(3)-tpi*g(3,ig)*eigr*imag
         enddo

         dg=matmul(bg,dg)
         xb(1)=xb(1)+dble(blip(1))**2
         xd(1)=xd(1)+dble(d)**2
         xbd(1)=xbd(1)+dble(blip(1))*dble(d)
         xb(2)=xb(2)+dble(blip(2))**2
         xd(2)=xd(2)+dble(d2)**2
         xbd(2)=xbd(2)+dble(blip(2))*dble(d2)
         xb(3)=xb(3)+dble(blip(3))**2
         xd(3)=xd(3)+dble(dg(1))**2
         xbd(3)=xbd(3)+dble(blip(3))*dble(dg(1))
         xb(4)=xb(4)+dble(blip(4))**2
         xd(4)=xd(4)+dble(dg(2))**2
         xbd(4)=xbd(4)+dble(blip(4))*dble(dg(2))
         xb(5)=xb(5)+dble(blip(5))**2
         xd(5)=xd(5)+dble(dg(3))**2
         xbd(5)=xbd(5)+dble(blip(5))*dble(dg(3))

        enddo

       endif

       if(ltest) then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2) / sqrt(xb(2)*xd(2)))**2, &
         &(xbd(3) / sqrt(xb(3)*xd(3)))**2, (xbd(4) / sqrt(xb(4)*xd(4)))**2,      &
         &(xbd(5) / sqrt(xb(5)*xd(5)))**2
       endif
     case(2)
! Read in the wavefunctions
       data = 0.d0
       data2 = 0.d0
       data3 = 0.d0
       data4 = 0.d0
       data5 = 0.d0
       do ig=1,npwx
        read(io2,*)evc(ig)
        icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
         &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
         &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
        data(icount) =                 evc(ig)*gamma(ig)
        data2(icount)=         -g2(ig)*evc(ig)*gamma(ig)
        data3(icount)=imag*gcart(1,ig)*evc(ig)*gamma(ig)
        data4(icount)=imag*gcart(2,ig)*evc(ig)*gamma(ig)
        data5(icount)=imag*gcart(3,ig)*evc(ig)*gamma(ig)
       enddo
       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data2,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data3,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data4,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data5,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data2=data2*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data3=data3*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data4=data4*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data5=data5*sqrt(dble(nr(1)*nr(2)*nr(3)))

       write(io,'(a)') 'Blip coefficients'
       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          avc(lx,ly,lz)=data(icount)
          avc2(lx,ly,lz)=data2(icount)
          avc3(lx,ly,lz)=data3(icount)
          avc4(lx,ly,lz)=data4(icount)
          avc5(lx,ly,lz)=data5(icount)
          icount=icount+1
         enddo
        enddo
       enddo

       write(io,'(a)') 'Orbital'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*) avc(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)') 'Laplacian'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*) avc2(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)') 'Gradient - a1 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*) avc3(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)') 'Gradient - a2 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*) avc4(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)') 'Gradient - a3 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*) avc5(lx,ly,lz)
         enddo
        enddo
       enddo
             
       if(ltest)then

        xb=0 ; xd=0 ; xbd =0
        do i=1,n_points_for_test 
         r(1)=ranfx() ; r(2)=ranfx() ; r(3)=ranfx()
         x=r(1) ; y=r(2) ; z=r(3)

         call blip3dk(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d=0
         d2=0
         dg=0

         do ig=1,npwx
           eigr=evc(ig)*exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
           d=d+eigr
           d2=d2-eigr*g2(ig)
           dg(1)=dg(1)-tpi*g(1,ig)*eigr*imag
           dg(2)=dg(2)-tpi*g(2,ig)*eigr*imag
           dg(3)=dg(3)-tpi*g(3,ig)*eigr*imag
         enddo

         dg=matmul(bg,dg)
         xb(1)=xb(1)+dble(blip(1))**2
         xd(1)=xd(1)+dble(d)**2
         xbd(1)=xbd(1)+dble(blip(1))*dble(d)
         xb(2)=xb(2)+dble(blip(2))**2
         xd(2)=xd(2)+dble(d2)**2
         xbd(2)=xbd(2)+dble(blip(2))*dble(d2)
         xb(3)=xb(3)+dble(blip(3))**2
         xd(3)=xd(3)+dble(dg(1))**2
         xbd(3)=xbd(3)+dble(blip(3))*dble(dg(1))
         xb(4)=xb(4)+dble(blip(4))**2
         xd(4)=xd(4)+dble(dg(2))**2
         xbd(4)=xbd(4)+dble(blip(4))*dble(dg(2))
         xb(5)=xb(5)+dble(blip(5))**2
         xd(5)=xd(5)+dble(dg(3))**2
         xbd(5)=xbd(5)+dble(blip(5))*dble(dg(3))

        enddo

       endif

       if(ltest) then
        write(6,'(" ",i4,5f14.10)') ibnd, (xbd(1) / sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2) / sqrt(xb(2)*xd(2)))**2, &
         &(xbd(3) / sqrt(xb(3)*xd(3)))**2, (xbd(4) / sqrt(xb(4)*xd(4)))**2,      &
         &(xbd(5) / sqrt(xb(5)*xd(5)))**2
       endif
     end select
    enddo

    if(ltest)write(6,*)

   enddo
  enddo

  if(.not.ltest)write(6,*)

 END SUBROUTINE blipk


 SUBROUTINE blipkr 
  USE globals
  USE singleton
  IMPLICIT NONE
  INTEGER,PARAMETER :: n_points_for_test = 100
  INTEGER i,j,l,nr(3),ibnd,ig,ig1,npwx,lx,ly,lz,nks,ik,nat,na,atomno, &
   &nbnd(2),ispin,nspin,nelec,isign,icount,indice,ierr
  INTEGER :: io=10,io2=11
  DOUBLE PRECISION ranfx
  DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0,tpi=2.d0*pi,eps=1.d-8
  COMPLEX(KIND=KIND(0.d0)) :: eigr,psi,d,d2,dg(3),blip(5)
  COMPLEX(KIND=KIND(0.d0)),PARAMETER :: imag=(0.d0,1.d0)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: avc(:,:,:), data(:)
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: avc2(:,:,:), data2(:) !Laplacian
  COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: avc3(:,:,:), data3(:),avc4(:,:,:), &
   & data4(:),avc5(:,:,:),data5(:) !gradient

  DOUBLE PRECISION x,y,z,k,alpha,gxmax,gymax,gzmax,alat,dax,day,daz, &
   &r(3),da(3),xb(5),xd(5),xbd(5),at(3,3),bg(3,3),tv,tc,tv0,tc0,tmp,       &
   &tau(3),etot,ewld,ecut,xk(3),et,arg,norm_real,norm_imag
  DOUBLE PRECISION,ALLOCATABLE :: evc(:),g(:,:), gcart(:,:), g2(:), gamma(:)
  CHARACTER(80) title
  LOGICAL lsda,lreal,limag

  write(6,*)'Real PW coefficients, K-POINTS'
  write(6,*)

  open(io2,file='pwfn.data',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening pwfn.data. Stopping.'
   stop
  endif
  open(io,file='bwfn.data',status='unknown',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening bwfn.data. Stopping.'
   stop
  endif

  call bwfn_info(io,io2,lsda)

  call skipio(io2,1) 
  write(io,'(a)')'Primitive lattice vectors (au) '
  read(io2,*)at(:,1)
  read(io2,*)at(:,2)
  read(io2,*)at(:,3)
  write(io,'(3f20.12)')at(:,1)
  write(io,'(3f20.12)')at(:,2)
  write(io,'(3f20.12)')at(:,3)

  call inve(at,bg)
  bg=transpose(bg)

  call skipio(io2,4) 
  write(io,'(a)')' '
  write(io,'(a)')'G VECTORS'
  write(io,'(a)')'---------'
  write(io,'(a)')'Number of G-vectors'
  read(io2,*)npwx
  write(io,*)npwx
  call skipio(io2,1) 
  write(io,'(a)')'Gx Gy Gz (au)'

  allocate(g(3,npwx),gcart(3,npwx),g2(npwx))
  gxmax=0 ; gymax=0 ; gzmax=0

  do ig=1,npwx
   read(io2,*) g(:,ig)
   write(io,'(3f20.12)')g(:,ig)
   g2(ig)=dot_product(g(:,ig),g(:,ig))
   gcart(:,ig)=g(:,ig)
   g(:,ig)=matmul(transpose(at),g(:,ig))/tpi
   if(abs(g(1,ig))>gxmax)gxmax=abs(g(1,ig))
   if(abs(g(2,ig))>gymax)gymax=abs(g(2,ig))
   if(abs(g(3,ig))>gzmax)gzmax=abs(g(3,ig))
  enddo

!  write(6,'('' Gmax:'',3f20.10)')gxmax, gymax, gzmax
  nr(1)=2*nint(gxmax)+2
  nr(2)=2*nint(gymax)+2
  nr(3)=2*nint(gzmax)+2
  nr=nr*xmul
! the following is to make nr an even number
  if(mod(nr(1),2)/=0)nr(1)=nr(1)+1
  if(mod(nr(2),2)/=0)nr(2)=nr(2)+1
  if(mod(nr(3),2)/=0)nr(3)=nr(3)+1
  write(6,'('' Grid : '',3i5)')nr
  write(6,*)
  dax=1.d0/dble(nr(1))
  day=1.d0/dble(nr(2))
  daz=1.d0/dble(nr(3))
  write(io,'(a)')'Blip grid'
  write(io,'(3i4)')nr

  call skipio(io2,4) 
  write(io,'(a)')' '
  write(io,'(a)')'WAVE FUNCTION'
  write(io,'(a)')'-------------'
  write(io,'(a)')'Number of k-points'
  read(io2,*)nks
  write(io,*)nks
  if(lsda)then
   nspin=2
  else
   nspin=1
  endif

  allocate(evc(npwx),gamma(npwx))
  select case (igrad_lap)
  case (0)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(data(nr(1)*nr(2)*nr(3)))
  case (1)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc2(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(data(nr(1)*nr(2)*nr(3)))
    allocate(data2(nr(1)*nr(2)*nr(3)))
  case (2)
    allocate(avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc2(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc3(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc4(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(avc5(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1))
    allocate(data(nr(1)*nr(2)*nr(3)))
    allocate(data2(nr(1)*nr(2)*nr(3)))
    allocate(data3(nr(1)*nr(2)*nr(3)))
    allocate(data4(nr(1)*nr(2)*nr(3)))
    allocate(data5(nr(1)*nr(2)*nr(3)))
  end select

! Calculating gamma
  do ig=1,npwx
   k=dax*g(1,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=1.5d0
   endif
   k=day*g(2,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=gamma(ig)*3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=gamma(ig)*1.5d0
   endif
   k=daz*g(3,ig)*tpi
   if(abs(k)>eps)then
    gamma(ig)=gamma(ig)*3.d0/(k)**4*(3.d0-4.d0*cos(k)+cos(2*k)) 
   else
    gamma(ig)=gamma(ig)*1.5d0
   endif
   gamma(ig)=1.d0/gamma(ig)
  enddo

  if(abs(g(1,1))>eps.or.abs(g(1,1))>eps.or.abs(g(1,1))>eps) &
   &stop 'first vector has to be the zero vector, reorder vectors in pwfn.data&
   &please'

  do ik=1,nks
   call skipio(io2,1)
   write(io,'(a)')'k-point # ; # of bands (up spin/down spin) ; k-point &
    &coords (au)'
   read(io2,*)i,nbnd(1),nbnd(2),xk(:)
   write(io,'(3i4,3f20.16)') ik, nbnd(1), nbnd(2), xk(:)
   do ispin=1,nspin 
!    xk(:) =  matmul( transpose(at), xk(:) ) / tpi
    write(6,'('' k-point '',i5,"  :  (",f10.6,",",f10.6,",",f10.6,")")')ik,xk

    if(ltest)then
     if(lsda)write(6,*)'Spin : ',ispin
     write(6,*)'Band ; alpha ; alpha (Laplacian) ; alpha (Gradient)'
    endif
    do ibnd=1,nbnd(ispin)
     call skipio(io2,1)
     write(io,'(a)')'Band, spin, eigenvalue (au)'
     read(io2,*)i,i,et 
     write(io,'(2i5,f20.12)')ibnd,ispin,et
     call skipio(io2,1)
! Read in the wavefunctions
     select case (igrad_lap)
     case(0)
       data=0.d0
       do ig=1,npwx
        read(io2,*)evc(ig)
        icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
         &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
         &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
        data(icount)=evc(ig)*gamma(ig)              
       enddo

       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))

       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          avc(lx,ly,lz)=data(icount)
          icount=icount+1
         enddo
        enddo
       enddo

       write(io,'(a)')'Blip coefficients'
   
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc(lx,ly,lz)
         enddo
        enddo
       enddo
             
       if(ltest)then

        xb=0 ; xd=0 ; xbd=0

        do i = 1, n_points_for_test 
            
         r(1) = ranfx(); r(2) = ranfx();  r(3) = ranfx()
         x = r(1); y = r(2); z = r(3)
             
         call blip3dk(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d=0 ; d2=0 ; dg=0
             
         do ig=1,npwx
          eigr=evc(ig)*exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
          d=d+eigr
          d2=d2-eigr*g2(ig)
          dg(1)=dg(1)-tpi*g(1,ig)*eigr*imag
          dg(2)=dg(2)-tpi*g(2,ig)*eigr*imag
          dg(3)=dg(3)-tpi*g(3,ig)*eigr*imag
         enddo
             
         dg=matmul(bg,dg)
         xb(1)=xb(1)+dble(blip(1))**2
         xd(1)=xd(1)+dble(d)**2
         xbd(1)=xbd(1)+dble(blip(1))*dble(d)
         xb(2)=xb(2)+dble(blip(2))**2
         xd(2)=xd(2)+dble(d2)**2
         xbd(2)=xbd(2)+dble(blip(2))*dble(d2)
         xb(3)=xb(3)+dble(blip(3))**2
         xd(3)=xd(3)+dble(dg(1))**2
         xbd(3)=xbd(3)+dble(blip(3))*dble(dg(1))
         xb(4)=xb(4)+dble(blip(4))**2
         xd(4)=xd(4)+dble(dg(2))**2
         xbd(4)=xbd(4)+dble(blip(4))*dble(dg(2))
         xb(5)=xb(5)+dble(blip(5))**2
         xd(5)=xd(5)+dble(dg(3))**2
         xbd(5)=xbd(5)+dble(blip(5))*dble(dg(3))

        enddo
        
       endif
     
       if(ltest)then
        write(6,'(" ",i4,5f14.10)')ibnd,&
         &(xbd(1)/sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2, &
         &(xbd(3)/sqrt(xb(3)*xd(3)))**2, &
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2, &
         &(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif
     case(1)
       data=0.d0
       data2=0.d0
       do ig=1,npwx
        read(io2,*)evc(ig)
        icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
         &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
         &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
        data(icount) =        evc(ig)*gamma(ig)              
        data2(icount)=-g2(ig)*evc(ig)*gamma(ig)
       enddo

       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data2,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data2=data2*sqrt(dble(nr(1)*nr(2)*nr(3)))

       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          avc(lx,ly,lz)=data(icount)
          avc2(lx,ly,lz)=data2(icount)
          icount=icount+1
         enddo
        enddo
       enddo

       write(io,'(a)')'Blip coefficients'
   
       write(io,'(a)')'Orbital'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Laplacian'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc2(lx,ly,lz)
         enddo
        enddo
       enddo
             
       if(ltest)then

        xb=0 ; xd=0 ; xbd=0

        do i = 1, n_points_for_test 
            
         r(1) = ranfx(); r(2) = ranfx();  r(3) = ranfx()
         x = r(1); y = r(2); z = r(3)
             
         call blip3dk(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d=0 ; d2=0 ; dg=0
             
         do ig=1,npwx
          eigr=evc(ig)*exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
          d=d+eigr
          d2=d2-eigr*g2(ig)
          dg(1)=dg(1)-tpi*g(1,ig)*eigr*imag
          dg(2)=dg(2)-tpi*g(2,ig)*eigr*imag
          dg(3)=dg(3)-tpi*g(3,ig)*eigr*imag
         enddo
             
         dg=matmul(bg,dg)
         xb(1)=xb(1)+dble(blip(1))**2
         xd(1)=xd(1)+dble(d)**2
         xbd(1)=xbd(1)+dble(blip(1))*dble(d)
         xb(2)=xb(2)+dble(blip(2))**2
         xd(2)=xd(2)+dble(d2)**2
         xbd(2)=xbd(2)+dble(blip(2))*dble(d2)
         xb(3)=xb(3)+dble(blip(3))**2
         xd(3)=xd(3)+dble(dg(1))**2
         xbd(3)=xbd(3)+dble(blip(3))*dble(dg(1))
         xb(4)=xb(4)+dble(blip(4))**2
         xd(4)=xd(4)+dble(dg(2))**2
         xbd(4)=xbd(4)+dble(blip(4))*dble(dg(2))
         xb(5)=xb(5)+dble(blip(5))**2
         xd(5)=xd(5)+dble(dg(3))**2
         xbd(5)=xbd(5)+dble(blip(5))*dble(dg(3))

        enddo
        
       endif
     
       if(ltest)then
        write(6,'(" ",i4,5f14.10)')ibnd,&
         &(xbd(1)/sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2, &
         &(xbd(3)/sqrt(xb(3)*xd(3)))**2, &
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2, &
         &(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif
     case(2)
       data=0.d0
       data2=0.d0
       data3=0.d0
       data4=0.d0
       data5=0.d0
       do ig=1,npwx
        read(io2,*)evc(ig)
        icount=1+nint(g(1,ig))+nr(1)*indice(nint(g(1,ig)))+ &
         &(nr(2)*indice(nint(g(2,ig)))+nint(g(2,ig)))*nr(1)+ &
         &(nr(3)*indice(nint(g(3,ig)))+nint(g(3,ig)))*nr(1)*nr(2)  
        data(icount) =                 evc(ig)*gamma(ig)
        data2(icount)=         -g2(ig)*evc(ig)*gamma(ig)
        data3(icount)=imag*gcart(1,ig)*evc(ig)*gamma(ig)
        data4(icount)=imag*gcart(2,ig)*evc(ig)*gamma(ig)
        data5(icount)=imag*gcart(3,ig)*evc(ig)*gamma(ig)
       enddo

       call fftn(data,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data2,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data3,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data4,(/nr(1),nr(2),nr(3)/),inv=.true.)
       call fftn(data5,(/nr(1),nr(2),nr(3)/),inv=.true.)
       data=data*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data2=data2*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data3=data3*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data4=data4*sqrt(dble(nr(1)*nr(2)*nr(3)))
       data5=data5*sqrt(dble(nr(1)*nr(2)*nr(3)))

       icount=1
       do lz=0,nr(3)-1
        do ly=0,nr(2)-1
         do lx=0,nr(1)-1
          avc(lx,ly,lz)=data(icount)
          avc2(lx,ly,lz)=data2(icount)
          avc3(lx,ly,lz)=data3(icount)
          avc4(lx,ly,lz)=data4(icount)
          avc5(lx,ly,lz)=data5(icount)
          icount=icount+1
         enddo
        enddo
       enddo

       write(io,'(a)')'Blip coefficients'
   
       write(io,'(a)')'Orbital'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Laplacian'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc2(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Gradient - a1 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc3(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Gradient - a2 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc4(lx,ly,lz)
         enddo
        enddo
       enddo
       write(io,'(a)')'Gradient - a3 direction'
       do lx=0,nr(1)-1
        do ly=0,nr(2)-1
         do lz=0,nr(3)-1 
          write(io,*)avc5(lx,ly,lz)
         enddo
        enddo
       enddo
             
       if(ltest)then

        xb=0 ; xd=0 ; xbd=0

        do i = 1, n_points_for_test 
            
         r(1) = ranfx(); r(2) = ranfx();  r(3) = ranfx()
         x = r(1); y = r(2); z = r(3)
             
         call blip3dk(blip(1),blip(2),blip(3:5),r,avc,nr,bg,1,1)
         d=0 ; d2=0 ; dg=0
             
         do ig=1,npwx
          eigr=evc(ig)*exp(tpi*imag*(g(1,ig)*x+g(2,ig)*y+g(3,ig)*z))
          d=d+eigr
          d2=d2-eigr*g2(ig)
          dg(1)=dg(1)-tpi*g(1,ig)*eigr*imag
          dg(2)=dg(2)-tpi*g(2,ig)*eigr*imag
          dg(3)=dg(3)-tpi*g(3,ig)*eigr*imag
         enddo
             
         dg=matmul(bg,dg)
         xb(1)=xb(1)+dble(blip(1))**2
         xd(1)=xd(1)+dble(d)**2
         xbd(1)=xbd(1)+dble(blip(1))*dble(d)
         xb(2)=xb(2)+dble(blip(2))**2
         xd(2)=xd(2)+dble(d2)**2
         xbd(2)=xbd(2)+dble(blip(2))*dble(d2)
         xb(3)=xb(3)+dble(blip(3))**2
         xd(3)=xd(3)+dble(dg(1))**2
         xbd(3)=xbd(3)+dble(blip(3))*dble(dg(1))
         xb(4)=xb(4)+dble(blip(4))**2
         xd(4)=xd(4)+dble(dg(2))**2
         xbd(4)=xbd(4)+dble(blip(4))*dble(dg(2))
         xb(5)=xb(5)+dble(blip(5))**2
         xd(5)=xd(5)+dble(dg(3))**2
         xbd(5)=xbd(5)+dble(blip(5))*dble(dg(3))

        enddo
        
       endif
     
       if(ltest)then
        write(6,'(" ",i4,5f14.10)')ibnd,&
         &(xbd(1)/sqrt(xb(1)*xd(1)))**2, &
         &(xbd(2)/sqrt(xb(2)*xd(2)))**2, &
         &(xbd(3)/sqrt(xb(3)*xd(3)))**2, &
         &(xbd(4)/sqrt(xb(4)*xd(4)))**2, &
         &(xbd(5)/sqrt(xb(5)*xd(5)))**2
       endif
     end select
   enddo

   if(ltest)write(6,*)
     
  enddo
 enddo

 if(.not.ltest)write(6,*)

 END SUBROUTINE blipkr


 SUBROUTINE blip3dk(rpsi,lap,grad,r,avc,nr,bg,iw,igl)
!------------------------------------------------------------------------------
! This subroutine evaluates the value of a function, its gradient and laplacian
! at a vector point r, using the overlapping of blip functions. The blip grid is
! defined on a cubic cell, so r should always be given in units of the crystal
! lattice vectors.
!
! Input:  
!        r(3)                 position in units of lattice vectors
!        avc(0:nr,0:nr,0:nr)  blips coefficients
!        nr(3)                number of divisions for each side of the box 
!                             (defines the blip grid)
!        bg(3,3)              the reciprocal lattice vectors (in a.u./tpi)
!------------------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER iw,igl,i,nr(3),ix,iy,iz,ixp,ixpp,ixm,iyp,iypp,iym,izp,izpp,izm
  DOUBLE PRECISION r(3),theta,dtheta,d2theta,t(3),txm,tx,txp,txpp,tym,ty, &
   &typ,typp,tzm,tz,tzp,tzpp,dtxm,dtx,dtxp,dtxpp,dtym,dty,dtyp,dtypp,dtzm,&
   &dtz,dtzp,dtzpp,d2txm,d2tx,d2txp,d2txpp,d2tym,d2ty,d2typ,d2typp,d2tzm, &
   &d2tz,d2tzp,d2tzpp,bg(3,3)
  COMPLEX(KIND=KIND(0.d0)) avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1),rpsi,lap,   &
   &grad(3),d1(16)
  EXTERNAL theta,dtheta,d2theta

  rpsi=0.d0
  lap=0.d0
  grad=0.d0
  ix=int(mod(r(1)+abs(int(r(1)))+1,1.d0)*nr(1))
  iy=int(mod(r(2)+abs(int(r(2)))+1,1.d0)*nr(2))
  iz=int(mod(r(3)+abs(int(r(3)))+1,1.d0)*nr(3))
  
  if(ix<0.or.iy<0.or.iz<0)stop 'Negative index in blip3d.'

! The blips are defined as the product of onedimensional cubic splines, 
! these are different from zero only on four adjacent grid points.  It 
! follows that for any vector r the are only 64 overlapping three-dimensional
! blips, which are the product of all possible combinations of the three 
! one-dimensional splines.

! These are the extra 3 coefficients for each dimension needed
  ix=mod(ix,nr(1))
  iy=mod(iy,nr(2))
  iz=mod(iz,nr(3))
  ixp=mod(ix+1,nr(1))
  ixpp=mod(ix+2,nr(1))
  ixm=mod(ix-1+nr(1),nr(1))
  iyp=mod(iy+1,nr(2))
  iypp=mod(iy+2,nr(2))
  iym=mod(iy-1+nr(2),nr(2))
  izp=mod(iz+1,nr(3))
  izpp=mod(iz+2,nr(3))
  izm=mod(iz-1+nr(3),nr(3))

! Now calculate the 12 monodimensional blip functions
  t=mod(r+abs(int(r))+1,1.d0)*nr
  t=mod(t,dble(nr))

  txm=theta(t(1)-ix+1.d0)
  tx=theta(t(1)-ix)
  txp=theta(t(1)-ix-1.d0)
  txpp=theta(t(1)-ix-2.d0)
  
  tym=theta(t(2)-iy+1.d0)
  ty=theta(t(2)-iy)
  typ=theta(t(2)-iy-1.d0)
  typp=theta(t(2)-iy-2.d0)
  
  tzm=theta(t(3)-iz+1.d0)
  tz=theta(t(3)-iz)
  tzp=theta(t(3)-iz-1.d0)
  tzpp=theta(t(3)-iz-2.d0)
       
  d1(1)=avc(ix,iy,iz)*tz+avc(ix,iy,izp)*tzp+avc(ix,iy,izpp)*tzpp+ &
   &avc(ix,iy,izm)*tzm 
  d1(2)=avc(ix,iyp,iz)*tz+avc(ix,iyp,izp)*tzp+avc(ix,iyp,izpp)*tzpp+ &
   &avc(ix,iyp,izm)*tzm 
  d1(3)=avc(ix,iypp,iz)*tz+avc(ix,iypp,izp)*tzp+avc(ix,iypp,izpp)*tzpp+ &
   &avc(ix,iypp,izm)*tzm 
  d1(4)=avc(ix,iym,iz)*tz+avc(ix,iym,izp)*tzp+avc(ix,iym,izpp)*tzpp+ &
   &avc(ix,iym,izm)*tzm 

  d1(5)=avc(ixp,iy,iz)*tz+avc(ixp,iy,izp)*tzp+avc(ixp,iy,izpp)*tzpp+ &
   &avc(ixp,iy,izm)*tzm 
  d1(6)=avc(ixp,iyp,iz)*tz+avc(ixp,iyp,izp)*tzp+avc(ixp,iyp,izpp)*tzpp+ &
   &avc(ixp,iyp,izm)*tzm 
  d1(7)=avc(ixp,iypp,iz)*tz+avc(ixp,iypp,izp)*tzp+avc(ixp,iypp,izpp)*tzpp+ &
   &avc(ixp,iypp,izm)*tzm 
  d1(8)=avc(ixp,iym,iz)*tz+avc(ixp,iym,izp)*tzp+avc(ixp,iym,izpp)*tzpp+ &
   &avc(ixp,iym,izm)*tzm 

  d1(9)=avc(ixpp,iy,iz)*tz+avc(ixpp,iy,izp)*tzp+avc(ixpp,iy,izpp)*tzpp+ &
   &avc(ixpp,iy,izm)*tzm 
  d1(10)=avc(ixpp,iyp,iz)*tz+avc(ixpp,iyp,izp)*tzp+avc(ixpp,iyp,izpp)*tzpp+ &
   &avc(ixpp,iyp,izm)*tzm 
  d1(11)=avc(ixpp,iypp,iz)*tz+avc(ixpp,iypp,izp)*tzp+avc(ixpp,iypp,izpp)*tzpp+ &
   &avc(ixpp,iypp,izm)*tzm 
  d1(12)=avc(ixpp,iym,iz)*tz+avc(ixpp,iym,izp)*tzp+avc(ixpp,iym,izpp)*tzpp+ &
   &avc(ixpp,iym,izm)*tzm 

  d1(13)=avc(ixm,iy,iz)*tz+avc(ixm,iy,izp)*tzp+avc(ixm,iy,izpp)*tzpp+ &
   &avc(ixm,iy,izm)*tzm 
  d1(14)=avc(ixm,iyp,iz)*tz+avc(ixm,iyp,izp)*tzp+avc(ixm,iyp,izpp)*tzpp+ &
   &avc(ixm,iyp,izm)*tzm 
  d1(15)=avc(ixm,iypp,iz)*tz+avc(ixm,iypp,izp)*tzp+avc(ixm,iypp,izpp)*tzpp+ &
   &avc(ixm,iypp,izm)*tzm 
  d1(16)=avc(ixm,iym,iz)*tz+avc(ixm,iym,izp)*tzp+avc(ixm,iym,izpp)*tzpp+ &
   &avc(ixm,iym,izm)*tzm 

! The function
  rpsi=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*txp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*txpp+&
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

  if(igl/=1)return

  dtxm =dtheta(t(1)-ix+1.d0)*nr(1)
  dtx  =dtheta(t(1)-ix     )*nr(1)
  dtxp =dtheta(t(1)-ix-1.d0)*nr(1)
  dtxpp=dtheta(t(1)-ix-2.d0)*nr(1)
  
  dtym =dtheta(t(2)-iy+1.d0)*nr(2)
  dty  =dtheta(t(2)-iy     )*nr(2)
  dtyp =dtheta(t(2)-iy-1.d0)*nr(2)
  dtypp=dtheta(t(2)-iy-2.d0)*nr(2)
  
  dtzm =dtheta(t(3)-iz+1.d0)*nr(3)
  dtz  =dtheta(t(3)-iz     )*nr(3)
  dtzp =dtheta(t(3)-iz-1.d0)*nr(3)
  dtzpp=dtheta(t(3)-iz-2.d0)*nr(3)
  
  d2txm =d2theta(t(1)-ix+1.d0)*nr(1)**2
  d2tx  =d2theta(t(1)-ix     )*nr(1)**2
  d2txp =d2theta(t(1)-ix-1.d0)*nr(1)**2
  d2txpp=d2theta(t(1)-ix-2.d0)*nr(1)**2
  
  d2tym =d2theta(t(2)-iy+1.d0)*nr(2)**2
  d2ty  =d2theta(t(2)-iy     )*nr(2)**2
  d2typ =d2theta(t(2)-iy-1.d0)*nr(2)**2
  d2typp=d2theta(t(2)-iy-2.d0)*nr(2)**2
  
  d2tzm =d2theta(t(3)-iz+1.d0)*nr(3)**2
  d2tz  =d2theta(t(3)-iz     )*nr(3)**2
  d2tzp =d2theta(t(3)-iz-1.d0)*nr(3)**2
  d2tzpp=d2theta(t(3)-iz-2.d0)*nr(3)**2
  
! the laplacian: first term involving Theta''(x)
  lap=((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*d2tx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*d2txp+&
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*d2txpp+&
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*d2txm)* &
   &(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2)
! this is to go from the lattice to the cartesian grid

! the laplacian: term involving Theta'(x)Theta'(y)
  lap=lap+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*dtx+ &
   &(d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*dtxp+ &
   &(d1(9)*dty+d1(10)*dtyp+d1(11)*dtypp+d1(12)*dtym)*dtxpp+ &
   &(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*dtym)*dtxm)* &
   &2*(bg(1,1)*bg(1,2)+bg(2,1)*bg(2,2)+bg(3,1)*bg(3,2)) 
! this is to go from the lattice to the cartesian grid

! the gradient, first term involving Theta'(x)
  grad(1)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*dtxp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*dtxpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm

! second term of the laplacian involving Theta''(y)
  lap=lap+((d1(1)*d2ty+d1(2)*d2typ+d1(3)*d2typp+d1(4)*d2tym)*tx+ &
   &(d1(5)*d2ty+d1(6)*d2typ+d1(7)*d2typp+d1(8)*d2tym)*txp+ &
   &(d1(9)*d2ty+d1(10)*d2typ+d1(11)*d2typp+d1(12)*d2tym)*txpp+ &
   &(d1(13)*d2ty+d1(14)*d2typ+d1(15)*d2typp+d1(16)*d2tym)*txm)* &
   &(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2 )    
! this is to go from the lattice to the cartesian grid

! second term of the gradient involving Theta'(y)
  grad(2)=(d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+ &
   (d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+ &
   (d1(9)*dty+d1(10)*dtyp+d1(11)*dtypp+d1(12)*dtym)*txpp+ &
   (d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*dtym)*txm

! and now the third term of the laplacian involving Theta''(z)
  d1(1)=avc(ix,iy,iz)*d2tz+avc(ix,iy,izp)*d2tzp+avc(ix,iy,izpp)*d2tzpp+ &
   &avc(ix,iy,izm)*d2tzm 
  d1(2)=avc(ix,iyp,iz)*d2tz+avc(ix,iyp,izp)*d2tzp+avc(ix,iyp,izpp)*d2tzpp+ &
   &avc(ix,iyp,izm)*d2tzm 
  d1(3)=avc(ix,iypp,iz)*d2tz+avc(ix,iypp,izp)*d2tzp+avc(ix,iypp,izpp)*d2tzpp+ &
   &avc(ix,iypp,izm)*d2tzm 
  d1(4)=avc(ix,iym,iz)*d2tz+avc(ix,iym,izp)*d2tzp+avc(ix,iym,izpp)*d2tzpp+ &
   &avc(ix,iym,izm)*d2tzm 

  d1(5)=avc(ixp,iy,iz)*d2tz+avc(ixp,iy,izp)*d2tzp+avc(ixp,iy,izpp)*d2tzpp+ &
   &avc(ixp,iy,izm)*d2tzm 
  d1(6)=avc(ixp,iyp,iz)*d2tz+avc(ixp,iyp,izp)*d2tzp+avc(ixp,iyp,izpp)*d2tzpp+ &
   &avc(ixp,iyp,izm)*d2tzm 
  d1(7)=avc(ixp,iypp,iz)*d2tz+avc(ixp,iypp,izp)*d2tzp+avc(ixp,iypp,izpp)*d2tzpp&
   &+avc(ixp,iypp,izm)*d2tzm 
  d1(8)=avc(ixp,iym,iz)*d2tz+avc(ixp,iym,izp)*d2tzp+avc(ixp,iym,izpp)*d2tzpp+ &
   &avc(ixp,iym,izm)*d2tzm 

  d1(9)=avc(ixpp,iy,iz)*d2tz+avc(ixpp,iy,izp)*d2tzp+avc(ixpp,iy,izpp)*d2tzpp+ &
   &avc(ixpp,iy,izm)*d2tzm 
  d1(10)=avc(ixpp,iyp,iz)*d2tz+avc(ixpp,iyp,izp)*d2tzp+avc(ixpp,iyp,izpp)* &
   &d2tzpp+avc(ixpp,iyp,izm)*d2tzm 
  d1(11)=avc(ixpp,iypp,iz)*d2tz+avc(ixpp,iypp,izp)*d2tzp+avc(ixpp,iypp,izpp)* &
   &d2tzpp+avc(ixpp,iypp,izm)*d2tzm 
  d1(12)=avc(ixpp,iym,iz)*d2tz+avc(ixpp,iym,izp)*d2tzp+avc(ixpp,iym,izpp)* &
   &d2tzpp+avc(ixpp,iym,izm)*d2tzm 

  d1(13)=avc(ixm,iy,iz)*d2tz+avc(ixm,iy,izp)*d2tzp+avc(ixm,iy,izpp)*d2tzpp+ &
   &avc(ixm,iy,izm)*d2tzm 
  d1(14)=avc(ixm,iyp,iz)*d2tz+avc(ixm,iyp,izp)*d2tzp+avc(ixm,iyp,izpp)*d2tzpp+ &
   &avc(ixm,iyp,izm)*d2tzm 
  d1(15)=avc(ixm,iypp,iz)*d2tz+avc(ixm,iypp,izp)*d2tzp+avc(ixm,iypp,izpp)* &
   &d2tzpp+avc(ixm,iypp,izm)*d2tzm 
  d1(16)=avc(ixm,iym,iz)*d2tz+avc(ixm,iym,izp)*d2tzp+avc(ixm,iym,izpp)*d2tzpp+ &
   &avc(ixm,iym,izm)*d2tzm 

! Theta''(z)
  lap=lap+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym )*txp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*txpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm)* &
   &(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2)    
! this is to go from the lattice to the cartesian grid

! and the third term of the gradient involving Theta'(z)
  d1(1)=avc(ix,iy,iz)*dtz+avc(ix,iy,izp)*dtzp+avc(ix,iy,izpp)*dtzpp+ &
   &avc(ix,iy,izm)*dtzm 
  d1(2)=avc(ix,iyp,iz)*dtz+avc(ix,iyp,izp)*dtzp+avc(ix,iyp,izpp)*dtzpp+ &
   &avc(ix,iyp,izm)*dtzm 
  d1(3)=avc(ix,iypp,iz)*dtz+avc(ix,iypp,izp)*dtzp+avc(ix,iypp,izpp)*dtzpp+ &
   &avc(ix,iypp,izm)*dtzm 
  d1(4)=avc(ix,iym,iz)*dtz+avc(ix,iym,izp)*dtzp+avc(ix,iym,izpp)*dtzpp+ &
   &avc(ix,iym,izm)*dtzm 

  d1(5)=avc(ixp,iy,iz)*dtz+avc(ixp,iy,izp)*dtzp+avc(ixp,iy,izpp)*dtzpp+ &
   &avc(ixp,iy,izm)*dtzm 
  d1(6)=avc(ixp,iyp,iz)*dtz+avc(ixp,iyp,izp)*dtzp+avc(ixp,iyp,izpp)*dtzpp+ &
   &avc(ixp,iyp,izm)*dtzm 
  d1(7)=avc(ixp,iypp,iz)*dtz+avc(ixp,iypp,izp)*dtzp+avc(ixp,iypp,izpp)*dtzpp+ &
   &avc(ixp,iypp,izm)*dtzm 
  d1(8)=avc(ixp,iym,iz)*dtz+avc(ixp,iym,izp)*dtzp+avc(ixp,iym,izpp)*dtzpp+ &
   &avc(ixp,iym,izm)*dtzm 

  d1(9)=avc(ixpp,iy,iz)*dtz+avc(ixpp,iy,izp)*dtzp+avc(ixpp,iy,izpp)*dtzpp+ &
   &avc(ixpp,iy,izm)*dtzm 
  d1(10)=avc(ixpp,iyp,iz)*dtz+avc(ixpp,iyp,izp)*dtzp+avc(ixpp,iyp,izpp)*dtzpp+ &
   &avc(ixpp,iyp,izm)*dtzm 
  d1(11)=avc(ixpp,iypp,iz)*dtz+avc(ixpp,iypp,izp)*dtzp+avc(ixpp,iypp,izpp)* &
   &dtzpp+avc(ixpp,iypp,izm)*dtzm 
  d1(12)=avc(ixpp,iym,iz)*dtz+avc(ixpp,iym,izp)*dtzp+avc(ixpp,iym,izpp)*dtzpp+ &
   &avc(ixpp,iym,izm)*dtzm 

  d1(13)=avc(ixm,iy,iz)*dtz+avc(ixm,iy,izp)*dtzp+avc(ixm,iy,izpp) &
   &*dtzpp+avc(ixm,iy,izm)*dtzm 
  d1(14)=avc(ixm,iyp,iz)*dtz+avc(ixm,iyp,izp)*dtzp+avc(ixm,iyp,izpp)*dtzpp+&
   &avc(ixm,iyp,izm)*dtzm 
  d1(15)=avc(ixm,iypp,iz)*dtz+avc(ixm,iypp,izp)*dtzp+avc(ixm,iypp,izpp)*dtzpp+&
   & avc(ixm,iypp,izm)*dtzm 
  d1(16)=avc(ixm,iym,iz)*dtz+avc(ixm,iym,izp)*dtzp+avc(ixm,iym,izpp)*dtzpp+ &
   &avc(ixm,iym,izm)*dtzm 

! The laplacian: term involving Theta'(x)Theta'(z)
  lap=lap+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*dtxp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*dtxpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm)* &
   &2*(bg(1,1)*bg(1,3)+bg(2,1)*bg(2,3)+bg(3,1)*bg(3,3)) 
! This is to go from the lattice to the cartesian grid

! The laplacian: term involving Theta'(y)Theta'(z)
  lap=lap+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+ &
   &(d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+ &
   &(d1(9)*dty+d1(10)*dtyp+d1(11)*dtypp+d1(12)*dtym)*txpp+ &
   &(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*dtym)*txm)* &
   &2*(bg(1,2)*bg(1,3)+bg(2,2)*bg(2,3)+bg(3,2)*bg(3,3)) 
! This is to go from the lattice to the cartesian grid

  grad(3)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*txp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*txpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm

! Go to the cartesian grid  
  grad=matmul(bg,grad)
       
 END SUBROUTINE blip3dk


 SUBROUTINE blip3d(rpsi,lap,grad,r,avc,nr,bg,iw,igl)
!-----------------------------------------------------------------------------
! This subroutine evaluates the value of a function, its gradient and its 
! laplacian at a vector point r, using the overlapping of blip functions. 
! The blip grid is defined on a cubic cell, so r should always be given in 
! units of the crystal lattice vectors.
!
! Input:  
!        r(3)                 position in units of lattice vectors
!        avc(0:nr,0:nr,0:nr)  blips coefficients
!        nr(3)                number of divisions for each side of  
!                             the box (defines the blip grid)
!        bg(3,3)              the reciprocal lattice vectors (in a.u./tpi)
!-----------------------------------------------------------------------------
  
  IMPLICIT NONE
  INTEGER iw,igl,i,nr(3),ix,iy,iz,ixp,ixpp,ixm,iyp,iypp,iym,izp,izpp,izm
  DOUBLE PRECISION avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1),r(3),theta,dtheta,    &
   &d2theta,t(3),txm,tx,txp,txpp,tym,ty,typ,typp,tzm,tz,tzp,tzpp,           &
   &dtxm,dtx,dtxp,dtxpp,dtym,dty,dtyp,dtypp,dtzm,dtz,dtzp,dtzpp,            &
   &d2txm,d2tx,d2txp,d2txpp,d2tym,d2ty,d2typ,d2typp,d2tzm,d2tz,d2tzp,d2tzpp,&
   &d1(16),bg(3,3),rpsi,lap,grad(3)
  EXTERNAL theta,dtheta,d2theta
  
  rpsi=0.d0
  lap=0.d0
  grad=0.d0
  ix=int(mod(r(1)+abs(int(r(1)))+1,1.d0)*nr(1))
  iy=int(mod(r(2)+abs(int(r(2)))+1,1.d0)*nr(2))
  iz=int(mod(r(3)+abs(int(r(3)))+1,1.d0)*nr(3))
  
  if(ix<0.or.iy<0.or.iz<0)stop 'negative index in blip3d'
  
! The blips are defined as the product of one-dimensional cubic splines, these 
! are different from zero only on four adjacent grid points. It follows that 
! for any vector r the are only 64 overlapping threedimensional blips, which 
! are the product of all possible combinations of the three one-dimensional 
! splines.
  
! These are the extra 3 coefficients for each dimension needed
  ix=mod(ix,nr(1))
  iy=mod(iy,nr(2))
  iz=mod(iz,nr(3))
  ixp=mod(ix+1,nr(1))
  ixpp=mod(ix+2,nr(1))
  ixm=mod(ix-1+nr(1),nr(1))
  iyp=mod(iy+1,nr(2))
  iypp=mod(iy+2,nr(2))
  iym=mod(iy-1+nr(2),nr(2) )
  izp=mod(iz+1,nr(3))
  izpp=mod(iz+2,nr(3))
  izm=mod(iz-1+nr(3),nr(3))
  
! Now calculate the 12 monodimensional blip functions
  t=mod(r+abs(int(r))+1,1.d0)*nr
  t=mod(t,dble(nr))
  
  txm =theta(t(1)-ix+1.d0)
  tx  =theta(t(1)-ix     )
  txp =theta(t(1)-ix-1.d0)
  txpp=theta(t(1)-ix-2.d0)
  
  tym =theta(t(2)-iy+1.d0)
  ty  =theta(t(2)-iy     )
  typ =theta(t(2)-iy-1.d0)
  typp=theta(t(2)-iy-2.d0)
  
  tzm =theta(t(3)-iz+1.d0)
  tz  =theta(t(3)-iz     )
  tzp =theta(t(3)-iz-1.d0)
  tzpp=theta(t(3)-iz-2.d0)
  
  d1(1)=avc(ix,iy,iz)*tz+avc(ix,iy,izp)*tzp+avc(ix,iy,izpp)*tzpp+avc(ix,iy,izm)&
   &*tzm 
  d1(2)=avc(ix,iyp,iz)*tz+avc(ix,iyp,izp)*tzp+avc(ix,iyp,izpp)*tzpp+avc(ix,iyp,&
   &izm)*tzm 
  d1(3)=avc(ix,iypp,iz)*tz+avc(ix,iypp,izp)*tzp+avc(ix,iypp,izpp)*tzpp+ &
   &avc(ix,iypp, izm)*tzm 
  d1(4)=avc(ix,iym,iz)*tz+avc(ix,iym,izp)*tzp+avc(ix,iym,izpp)*tzpp+ &
   &avc(ix,iym,izm)*tzm 
  
  d1(5)=avc(ixp,iy,iz)*tz+avc(ixp,iy,izp)*tzp+avc(ixp,iy,izpp)*tzpp+ &
   &avc(ixp,iy,izm)*tzm 
  d1(6)=avc(ixp,iyp,iz)*tz+avc(ixp,iyp,izp)*tzp+avc(ixp,iyp,izpp)*tzpp+ &
   &avc(ixp, iyp, izm)*tzm 
  d1(7)=avc(ixp,iypp,iz)*tz+avc(ixp,iypp,izp)*tzp+avc(ixp,iypp,izpp)*tzpp+ &
   &avc(ixp,iypp,izm)*tzm 
  d1(8)=avc(ixp,iym,iz)*tz+avc(ixp,iym,izp)*tzp+avc(ixp,iym,izpp)*tzpp+ &
   &avc(ixp,iym,izm)*tzm 
  
  d1(9)=avc(ixpp,iy,iz)*tz+avc(ixpp,iy,izp)*tzp+avc(ixpp,iy,izpp)*tzpp+ &
   &avc(ixpp,iy,izm)*tzm 
  d1(10)=avc(ixpp,iyp,iz)*tz+avc(ixpp,iyp,izp)*tzp+avc(ixpp,iyp,izpp)*tzpp+ &
   &avc(ixpp,iyp,izm)*tzm 
  d1(11)=avc(ixpp,iypp,iz)*tz+avc(ixpp,iypp,izp)*tzp+avc(ixpp,iypp,izpp)*tzpp+ &
   &avc(ixpp,iypp,izm)*tzm 
  d1(12)=avc(ixpp,iym,iz)*tz+avc(ixpp,iym,izp)*tzp+avc(ixpp,iym,izpp)*tzpp+ &
   &avc(ixpp,iym,izm)*tzm 
  
  d1(13)=avc(ixm,iy,iz)*tz+avc(ixm,iy,izp)*tzp+avc(ixm,iy,izpp)*tzpp+ &
   &avc(ixm,iy,izm)*tzm 
  d1(14)=avc(ixm,iyp,iz)*tz+avc(ixm,iyp,izp)*tzp+avc(ixm,iyp,izpp)*tzpp+ &
   &avc(ixm,iyp,izm)*tzm 
  d1(15)=avc(ixm,iypp,iz)*tz+avc(ixm,iypp,izp)*tzp+avc(ixm,iypp,izpp)*tzpp+ &
   &avc(ixm,iypp,izm)*tzm 
  d1(16)=avc(ixm,iym,iz)*tz+avc(ixm,iym,izp)*tzp+avc(ixm,iym,izpp)*tzpp+ &
   &avc(ixm,iym,izm)*tzm 
  
! the function
  rpsi=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*txp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*txpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm
  
  if(igl/=1) return
  
  dtxm =dtheta(t(1)-ix+1.d0)*nr(1)
  dtx  =dtheta(t(1)-ix     )*nr(1)
  dtxp =dtheta(t(1)-ix-1.d0)*nr(1)
  dtxpp=dtheta(t(1)-ix-2.d0)*nr(1)
  
  dtym =dtheta(t(2)-iy+1.d0)*nr(2)
  dty  =dtheta(t(2)-iy     )*nr(2)
  dtyp =dtheta(t(2)-iy-1.d0)*nr(2)
  dtypp=dtheta(t(2)-iy-2.d0)*nr(2)
  
  dtzm =dtheta(t(3)-iz+1.d0)*nr(3)
  dtz  =dtheta(t(3)-iz     )*nr(3)
  dtzp =dtheta(t(3)-iz-1.d0)*nr(3)
  dtzpp=dtheta(t(3)-iz-2.d0)*nr(3)
  
  d2txm =d2theta(t(1)-ix+1.d0)*nr(1)**2
  d2tx  =d2theta(t(1)-ix     )*nr(1)**2
  d2txp =d2theta(t(1)-ix-1.d0)*nr(1)**2
  d2txpp=d2theta(t(1)-ix-2.d0)*nr(1)**2
  
  d2tym =d2theta(t(2)-iy+1.d0)*nr(2)**2
  d2ty  =d2theta(t(2)-iy     )*nr(2)**2
  d2typ =d2theta(t(2)-iy-1.d0)*nr(2)**2
  d2typp=d2theta(t(2)-iy-2.d0)*nr(2)**2
  
  d2tzm =d2theta(t(3)-iz+1.d0)*nr(3)**2
  d2tz  =d2theta(t(3)-iz     )*nr(3)**2
  d2tzp =d2theta(t(3)-iz-1.d0)*nr(3)**2
  d2tzpp=d2theta(t(3)-iz-2.d0)*nr(3)**2
  
! the laplacian: first term involving Theta''(x)
  lap=((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*d2tx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*d2txp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*d2txpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*d2txm)* &
   &(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2)    
! This is to go from the lattice to the cartesian grid
  
! The laplacian: term involving Theta'(x)Theta'(y)
  lap=lap+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*dtx+ &
   &(d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*dtxp+ &
   &(d1(9)*dty+d1(10)*dtyp+d1(11)*dtypp+d1(12)*dtym)*dtxpp+ &
   &(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*dtym)*dtxm)* &
   &2*(bg(1,1)*bg(1,2)+bg(2,1)*bg(2,2)+bg(3,1)*bg(3,2)) 
! This is to go from the lattice to the cartesian grid
  
  
! The gradient, first term involving Theta'(x)
  grad(1)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*dtxp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*dtxpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm
  
  
! Second term of the laplacian involving Theta''(y)
  lap=lap+((d1(1)*d2ty+d1(2)*d2typ+d1(3)*d2typp+d1(4)*d2tym)*tx+ &
   &(d1(5)*d2ty+d1(6)*d2typ+d1(7)*d2typp+d1(8)*d2tym)*txp+ &
   &(d1(9)*d2ty+d1(10)*d2typ+d1(11)*d2typp+d1(12)*d2tym)*txpp+ &
   &(d1(13)*d2ty+d1(14)*d2typ+d1(15)*d2typp+d1(16)*d2tym)*txm)* &
   &(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2)    
! This is to go from the lattice to the cartesian grid
  
! Second term of the gradient involving Theta'(y)
  grad(2)=(d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+ &
   &(d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+ &
   &(d1(9)*dty+d1(10)*dtyp+d1(11)*dtypp+d1(12)*dtym)*txpp+ &
   &(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*dtym)*txm
  
! And now the third term of the laplacian involving Theta''(z)
  d1(1)=avc(ix,iy,iz)*d2tz+avc(ix,iy,izp)*d2tzp+avc(ix,iy,izpp)*d2tzpp+ &
   &avc(ix,iy,izm)*d2tzm 
  d1(2)=avc(ix,iyp,iz)*d2tz+avc(ix,iyp,izp)*d2tzp+avc(ix,iyp,izpp)*d2tzpp+ &
   &avc(ix,iyp,izm)*d2tzm 
  d1(3)=avc(ix,iypp,iz)*d2tz+avc(ix,iypp,izp)*d2tzp+avc(ix,iypp,izpp)*d2tzpp+ &
   &avc(ix,iypp,izm)*d2tzm 
  d1(4)=avc(ix,iym,iz)*d2tz+avc(ix,iym,izp)*d2tzp+avc(ix,iym,izpp)*d2tzpp+ &
   &avc(ix,iym,izm)*d2tzm 
  
  d1(5)=avc(ixp,iy,iz)*d2tz+avc(ixp,iy,izp)*d2tzp+avc(ixp,iy,izpp)*d2tzpp+ &
   &avc(ixp,iy,izm)*d2tzm 
  d1(6)=avc(ixp,iyp,iz)*d2tz+avc(ixp,iyp,izp)*d2tzp+avc(ixp,iyp,izpp)*d2tzpp+ &
   &avc(ixp,iyp,izm)*d2tzm 
  d1(7)=avc(ixp,iypp,iz)*d2tz+avc(ixp,iypp,izp)*d2tzp+avc(ixp,iypp,izpp)* &
   &d2tzpp+avc(ixp,iypp,izm)*d2tzm 
  d1(8)=avc(ixp,iym,iz)*d2tz+avc(ixp,iym,izp)*d2tzp+avc(ixp,iym,izpp)*d2tzpp+ &
   &avc(ixp, iym, izm)*d2tzm 
  
  d1(9)=avc(ixpp,iy,iz)*d2tz+avc(ixpp,iy,izp)*d2tzp+avc(ixpp,iy,izpp)*d2tzpp+ &
   &avc(ixpp,iy,izm)*d2tzm 
  d1(10)=avc(ixpp,iyp,iz)*d2tz+avc(ixpp,iyp,izp)*d2tzp+avc(ixpp,iyp,izpp)* &
   &d2tzpp+avc(ixpp,iyp,izm)*d2tzm 
  d1(11)=avc(ixpp,iypp,iz)*d2tz+avc(ixpp,iypp,izp)*d2tzp+avc(ixpp,iypp,izpp)* &
   &d2tzpp+avc(ixpp,iypp,izm)*d2tzm 
  d1(12)=avc(ixpp,iym,iz)*d2tz+avc(ixpp,iym,izp)*d2tzp+avc(ixpp,iym,izpp)* &
   &d2tzpp+avc(ixpp,iym,izm)*d2tzm 
  
  d1(13)=avc(ixm,iy,iz)*d2tz+avc(ixm,iy,izp)*d2tzp+avc(ixm,iy,izpp)*d2tzpp+ &
   &avc(ixm,iy,izm)*d2tzm 
  d1(14)=avc(ixm,iyp,iz)*d2tz+avc(ixm,iyp,izp)*d2tzp+avc(ixm,iyp,izpp)*d2tzpp+ &
   &avc(ixm,iyp,izm)*d2tzm 
  d1(15)=avc(ixm,iypp,iz)*d2tz+avc(ixm,iypp,izp)*d2tzp+avc(ixm,iypp,izpp)* &
   &d2tzpp+avc(ixm,iypp,izm)*d2tzm 
  d1(16)=avc(ixm,iym,iz)*d2tz+avc(ixm,iym,izp)*d2tzp+avc(ixm,iym,izpp)*d2tzpp+ &
   &avc(ixm,iym,izm)*d2tzm 
  
! Theta''(z)
  lap=lap+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*txp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*txpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm)* &
   &(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2)    
! This is to go from the lattice to the cartesian grid
  
! And the third term of the gradient involving Theta'(z)
  d1(1)=avc(ix,iy,iz)*dtz+avc(ix,iy,izp)*dtzp+avc(ix,iy,izpp)*dtzpp+ &
   &avc(ix,iy,izm)*dtzm 
  d1(2)=avc(ix,iyp,iz)*dtz+avc(ix,iyp,izp)*dtzp+avc(ix,iyp,izpp)*dtzpp+ &
   &avc(ix,iyp,izm)*dtzm 
  d1(3)=avc(ix,iypp,iz)*dtz+avc(ix,iypp,izp)*dtzp+avc(ix,iypp,izpp)*dtzpp+ &
   &avc(ix,iypp,izm)*dtzm 
  d1(4)=avc(ix,iym,iz)*dtz+avc(ix,iym,izp)*dtzp+avc(ix,iym,izpp)*dtzpp+ &
   &avc(ix, iym, izm)*dtzm 
  
  d1(5)=avc(ixp,iy,iz)*dtz+avc(ixp,iy,izp)*dtzp+avc(ixp,iy,izpp)*dtzpp+ &
   &avc(ixp,iy,izm)*dtzm 
  d1(6)=avc(ixp,iyp,iz)*dtz+avc(ixp,iyp,izp)*dtzp+avc(ixp,iyp,izpp)*dtzpp+ &
   &avc(ixp,iyp,izm)*dtzm 
  d1(7)=avc(ixp,iypp,iz)*dtz+avc(ixp,iypp,izp)*dtzp+avc(ixp,iypp,izpp)*dtzpp+ &
   &avc(ixp,iypp,izm)*dtzm 
  d1(8)=avc(ixp,iym,iz)*dtz+avc(ixp,iym,izp)*dtzp+avc(ixp,iym,izpp)*dtzpp+ &
   &avc(ixp, iym, izm)*dtzm 
  
  d1(9)=avc(ixpp,iy,iz)*dtz+avc(ixpp,iy,izp)*dtzp+avc(ixpp,iy,izpp)*dtzpp+ &
   &avc(ixpp,iy,izm)*dtzm 
  d1(10)=avc(ixpp,iyp,iz)*dtz+avc(ixpp,iyp,izp)*dtzp+avc(ixpp,iyp,izpp)*dtzpp+ &
   &avc(ixpp,iyp,izm)*dtzm 
  d1(11)=avc(ixpp,iypp,iz)*dtz+avc(ixpp,iypp,izp)*dtzp+avc(ixpp,iypp,izpp)* &
   &dtzpp+avc(ixpp,iypp,izm)*dtzm 
  d1(12)=avc(ixpp,iym,iz)*dtz+avc(ixpp,iym,izp)*dtzp+avc(ixpp,iym,izpp)*dtzpp+ &
   &avc(ixpp,iym,izm)*dtzm 
  
  d1(13)=avc(ixm,iy,iz)*dtz+avc(ixm,iy,izp)*dtzp+avc(ixm,iy,izpp)*dtzpp+ &
   &avc(ixm,iy,izm)*dtzm 
  d1(14)=avc(ixm,iyp,iz)*dtz+avc(ixm,iyp,izp)*dtzp+avc(ixm,iyp,izpp)*dtzpp+ &
   &avc(ixm,iyp,izm)*dtzm 
  d1(15)=avc(ixm,iypp,iz)*dtz+avc(ixm,iypp,izp)*dtzp+avc(ixm,iypp,izpp)*dtzpp+ &
   &avc(ixm,iypp,izm)*dtzm 
  d1(16)=avc(ixm,iym,iz)*dtz+avc(ixm,iym,izp)*dtzp+avc(ixm,iym,izpp)*dtzpp+ &
   &avc(ixm,iym,izm)*dtzm 
  
! The laplacian: term involving Theta'(x)Theta'(z)
  lap=lap+((d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*dtx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*dtxp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*dtxpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*dtxm)* &
   &2*(bg(1,1)*bg(1,3)+bg(2,1)*bg(2,3)+bg(3,1)*bg(3,3)) 
! This is to go from the lattice to the cartesian grid
  
! The laplacian: term involving Theta'(y)Theta'(z)
  lap=lap+((d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym)*tx+ &
   &(d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym)*txp+ &
   &(d1(9)*dty+d1(10)*dtyp+d1(11)*dtypp+d1(12)*dtym)*txpp+ &
   &(d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*dtym)*txm)* &
   &2*(bg(1,2)*bg(1,3)+bg(2,2)*bg(2,3)+bg(3,2)*bg(3,3)) 
! This is to go from the lattice to the cartesian grid
  
  grad(3)=(d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym)*tx+ &
   &(d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym)*txp+ &
   &(d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym)*txpp+ &
   &(d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym)*txm
  
! Go to the cartesian grid  
  grad=matmul(bg,grad)
  
 END SUBROUTINE blip3d


 FUNCTION theta(x)
  IMPLICIT NONE
  DOUBLE PRECISION x,theta
  
  theta=0.d0
  if(abs(x)>2)then
   return
  elseif(abs(x)<1)then
   theta=1.d0-1.5d0*x**2+0.75d0*abs(x)**3
   return
  else
   theta=0.25d0*(2.d0-abs(x))**3  
  endif
  
 END FUNCTION theta


 FUNCTION dtheta(x)
  IMPLICIT NONE
  DOUBLE PRECISION x,dtheta
  
  dtheta=0.d0
  if(abs(x)>2)then
   return
  elseif(abs(x)<1)then
   dtheta=-3.d0*x+3*0.75d0*abs(x)*x
   return
  else
   dtheta=-0.75d0*(2.d0-abs(x))**2*x/abs(x)
  endif
  
 END FUNCTION dtheta


 FUNCTION d2theta(x)
  IMPLICIT NONE
  DOUBLE PRECISION x,d2theta
  
  d2theta=0.d0
  if(abs(x)>2)then
   return
  elseif(abs(x)<1)then
   d2theta=-3.d0+6*0.75d0*abs(x) 
   return
  else
   d2theta=1.5d0*(2.d0-abs(x))
  endif
  
 END FUNCTION d2theta


 SUBROUTINE skipio(io,n)
  do i=1,n
   read(io,*)
  enddo
 END SUBROUTINE skipio


 SUBROUTINE inve(v,inv)
!-----------------------!
! Inverts 3x3 matrices. !
!-----------------------!
  IMPLICIT NONE
  DOUBLE PRECISION v(3,3),inv(3,3),d
  d=v(1,1)*(v(2,2)*v(3,3)-v(2,3)*v(3,2))+ &
   &v(2,1)*(v(3,2)*v(1,3)-v(1,2)*v(3,3))+ &
   &v(3,1)*(v(1,2)*v(2,3)-v(1,3)*v(2,2))
  inv(1,1)=(v(2,2)*v(3,3)-v(2,3)*v(3,2))/d
  inv(1,2)=(v(3,2)*v(1,3)-v(1,2)*v(3,3))/d
  inv(1,3)=(v(1,2)*v(2,3)-v(1,3)*v(2,2))/d
  inv(2,1)=(v(3,1)*v(2,3)-v(2,1)*v(3,3))/d
  inv(2,2)=-(v(3,1)*v(1,3)-v(1,1)*v(3,3))/d
  inv(2,3)=(v(2,1)*v(1,3)-v(1,1)*v(2,3))/d
  inv(3,1)=(v(2,1)*v(3,2)-v(2,2)*v(3,1))/d
  inv(3,2)=(v(3,1)*v(1,2)-v(1,1)*v(3,2))/d
  inv(3,3)=(v(1,1)*v(2,2)-v(1,2)*v(2,1))/d
 END SUBROUTINE inve


 FUNCTION ranfx()
 !--------------------------!
 ! Random number generator. !
 !--------------------------!
  IMPLICIT DOUBLE PRECISION (a-h,o-z)
  DATA m12/4096/
  DATA f1/2.44140625d-04/,f2/5.96046448d-08/
  DATA f3/1.45519152d-11/
  DATA j1/3823/,j2/4006/,j3/2903/
  DATA i1/3823/,i2/4006/,i3/2903/
  k3=i3*j3
  l3=k3/m12
  k2=i2*j3+i3*j2+l3
  l2=k2/m12
  k1=i1*j3+i2*j2+i3*j1+l2
  l1=k1/m12
  i1=k1-l1*m12
  i2=k2-l2*m12
  i3=k3-l3*m12
  ranfx=f1*dble(i1)+f2*dble(i2)+f3*dble(i3)
 END FUNCTION ranfx


 FUNCTION indice(i)
  IMPLICIT NONE
  INTEGER i,indice
  if(i>=0)then
   indice=0
  else
   indice=1
  endif
 END FUNCTION indice


 SUBROUTINE bwfn_info(io,io2,lsda)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: io,io2
  LOGICAL,INTENT(out) :: lsda
  INTEGER atomno,nat,na,nelec
  DOUBLE PRECISION energy,tau(3)
  CHARACTER(80) char_80
  read(io2,'(a)')char_80
  write(io,'(a)')trim(char_80)
  call skipio(io2,4)
  write(io,'(a)')
  write(io,'(a)')'BASIC INFO'
  write(io,'(a)')'----------'
  write(io,'(a)')'Generated by:'
  read(io2,'(a)')char_80
  write(io,'(a)')trim(char_80)
  call skipio(io2,1)
  write(io,'(a)')'Method:'
  read(io2,'(a)')char_80
  write(io,'(a)')trim(char_80)
  call skipio(io2,1)
  write(io,'(a)')'DFT Functional:'         
  read(io2,'(a)')char_80
  write(io,'(a)')trim(char_80)
  call skipio(io2,1)
  write(io,'(a)')'Pseudopotential'
  read(io2,'(a)')char_80
  write(io,'(a)')trim(char_80)
  call skipio(io2,1)
  write(io,'(a)')'Plane wave cutoff (au)'
  read(io2,*)energy
  write(io,*)energy
  call skipio(io2,1)
  write(io,'(a)')'Spin polarized:'
  read(io2,*)lsda
  write(io,*)lsda 
  call skipio(io2,1)
  write(io,'(a)')'Total energy (au per primitive cell)' 
  read(io2,*)energy
  write(io,*)energy         
  call skipio(io2,1)
  write(io,'(a)')'Kinetic energy (au per primitive cell)' 
  read(io2,*)energy
  write(io,*)energy              
  call skipio(io2,1)
  write(io,'(a)')'Local potential energy (au per primitive cell)' 
  read(io2,*)energy
  write(io,*)energy 
  call skipio(io2,1)
  write(io,'(a)')'Non local potential energy(au per primitive cell)'
  read(io2,*)energy
  write(io,*)energy
  call skipio(io2,1)
  write(io,'(a)')'Electron electron energy (au per primitive cell)' 
  read(io2,*)energy
  write(io,*)energy
  call skipio(io2,1)
  write(io,'(a)')'Ion-ion energy (au per primitive cell)' 
  read(io2,*)energy
  write(io,*)energy
  call skipio(io2,1)
  write(io,'(a)')'Number of electrons per primitive cell'                 
  read(io2,*)nelec
  write(io,*)nelec
  call skipio(io2,4)
  write(io,'(a)')' '                 
  write(io,'(a)')'GEOMETRY'
  write(io,'(a)')'-------- '
  write(io,'(a)')'Number of atoms per primitive cell '
  read(io2,*)nat
  write(io,*)nat
  call skipio(io2,1) 
  write(io,'(a)')'Atomic number and position of the atoms(au) '
  do na=1,nat
   read(io2,*)atomno,tau(:)
   write(io,'(i6,3f20.12)') atomno, tau(:)
  enddo
 END SUBROUTINE bwfn_info

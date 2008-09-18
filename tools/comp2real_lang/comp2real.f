      subroutine expand_det

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: expands each orbital in each determinant
*              into real and complex parts. This done by doubling
*              the number of determinants each time we expand an orbital
* arguments: none
* return values: det, ndet, cdet are modified.
*---------------------------------------------------------------------

      include         'maxdim.h'
      include         'qnumbers.h'

c locals:
      integer         in,idet,abs_m
      integer         ndet_new
      integer         det_new(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet_new(MAXNDET)
      real*8          sqrt2

c common variables:
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)

      common/psi_n/   cdet,ndim,nelec,ndet,det

      write(*,*) '* Expanding determinants'

      sqrt2=dsqrt(2.d0)

      do 10 in=1,nelec
        ndet_new=0
        do 7 idet=1,ndet
          ndet_new=ndet_new+1

c copy determinant
          call copy_det(nelec,idet,det,ndet_new,det_new)
          cdet_new(ndet_new)=cdet(idet)

c determinant is doubled only if abs(m)>0:
          abs_m=abs(det(idet,in,qnm))
          if (abs(det(idet,in,qnm)).gt.0) then

c first, deal with +m state
            det_new(ndet_new,in,qnm)=abs_m
            cdet_new(ndet_new)=cdet(idet)/sqrt2

c now create -m state: copy determinant, change sign and amplitude.
            ndet_new=ndet_new+1
            call copy_det(nelec,idet,det,ndet_new,det_new)
            det_new(ndet_new,in,qnm)=-abs_m
            if (det(idet,in,qnm).gt.0) then
              cdet_new(ndet_new)=cdet(idet)
     &                                *dcmplx(0.d0,1.d0)/sqrt2
            else
              cdet_new(ndet_new)=cdet(idet)
     &                                *dcmplx(0.d0,-1.d0)/sqrt2
            endif
          endif
 7      enddo

        ndet=ndet_new
        if (ndet.gt.MAXNDET) then
          write(*,*) '  PROBLEM! Too much determinants created'
          stop
        endif


c overwrite det()
        do 8 idet=1,ndet
          call copy_det(nelec,idet,det_new,idet,det)
          cdet(idet)=cdet_new(idet)
 8      enddo

 10   enddo

      write(*,*) '  Number of determinants after expansion = ',ndet
      write(*,*) '  OK.'
      write(*,*)
      return
      end

**********************************************************************

      subroutine copy_det(nelec,idet1,det1,idet2,det2)

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: copies idet1^{th} determinant of det1 to idet2^{th}
*              determinant of det2
*---------------------------------------------------------------------

      include         'maxdim.h'

c arguments
      integer         nelec,idet1,idet2
      integer         det1(MAXNDET,MAXNELEC,MAXQN)
      integer         det2(MAXNDET,MAXNELEC,MAXQN)

c locals
      integer         in,iq

      do iq=1,MAXQN
        do in=1,nelec
          det2(idet2,in,iq)=det1(idet1,in,iq)
        enddo
      enddo

      return
      end

**********************************************************************

      subroutine sign_3D

c Written by Devrim Guclu for Cyrus Umrigar
*---------------------------------------------------------------------
* description: Applies sign convention of Spherical Harmonics in 3D
*              sys. assuming that it has not been applied in the input
*              data.
* arguments: none
* return values: cdet is modified for 3D systems.
*---------------------------------------------------------------------

      include         'maxdim.h'
      include         'qnumbers.h'

c locals:
      integer         in,idet,im
      integer         exp
      real*8          sign
c common variables:
      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)

      common/psi_n/   cdet,ndim,nelec,ndet,det

      if (ndim.eq.3) then
        write(*,*) '* WARNING: applying sign convention for ',
     &              'spherical harmonics'
        write(*,*) '  New determinant coefficients are = '
        do 20 idet=1,ndet
          exp=0
          do 15 in=1,nelec
            im=det(idet,in,qnm)
            exp=exp+im+abs(im)
 15       enddo
          sign=(-1)**(dfloat(exp)/2.d0)
          cdet(idet)=cdet(idet)*sign
          write(*,*) '   ',dreal(cdet(idet))
 20     enddo
        write(*,*) '  OK.'
        write(*,*)
      endif

      return
      end

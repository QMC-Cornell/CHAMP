! ==================================================================== *
      subroutine CNVDBL (String, value)
! -------------------------------------------------------------------- *
!
! Comments: Converts a string into a real*8 value
!
! author:F.Colonna
! date  :20 Dec 89
!--------------------------------------------------------------------- *
!
      implicit none
! i/o scalars:
      character*(*) String
      character*(30) Stemp
      double precision value
! local scalars:
      integer lenstr

! begin:
      call strim (String,lenstr)

      if(lenstr.gt.0) then
         Stemp = '                               '
         Stemp(31-lenstr:30) = String(1:lenstr)
         read(Stemp,'(D30.6)',err=999) value

      elseif(lenstr.lt.0) then
         write(6,'(a,i10)') 'error_CNVDBL> string length =',lenstr
         stop ' error_CNVDBL> string length < 0'

      else
         write(6,'(a)') ' warning_CNVDBL> zero length string. converted to 0.d0'
         value = 0.d0
      endif

      return

 999  continue
      write(6,'(3a)') 'warning_CNVDBL> conversion error for string "', String(1:lenstr),'" converted to 1.d+99'
      stop 'CNVDBL: conversion error'

! end of routine CNVDBL
      end

! ==================================================================== *
      subroutine CNVINT (String, value)
! -------------------------------------------------------------------- *
!
! Comments: Converts a String into an integer value
!
! author:F.Colonna
! date  :02 Jan 90
! -------------------------------------------------------------------- *
      implicit none
      character*(*) string
      character*30  text
      integer value, lenstr

      call strim (string, lenstr)

      if(lenstr.ne.0) then
         text = '                                        '
         text(31-lenstr:30) = string(1:lenstr)
         read(text,'(i30)') value
      else
         write(6,'(a)') 'warning_CNVINT> zero length string. converted to 0'
         value = 0
      endif

      return
      end

! ==================================================================== *
      subroutine CNVTLC (chr)
! -------------------------------------------------------------------- *
!
! Comments: Converts one character from Upper-case to Lower-case
!
! author:F.Colonna
! date  :20 Fev 91
! -------------------------------------------------------------------- *
      implicit none

! input/output scalars:
      character*1 chr

! local scalars:
      integer ichr
      character*1 lita, litz, biga, bigz, space
      parameter (lita='a',litz='z',biga='A',bigz='Z',space=' ')

! begin:
      ichr = ichar(chr)
      if(ichr.ge.ichar(biga).and.ichr.le.ichar(bigz)) then
        chr = char(ichr-(ichar(biga)-ichar(lita)))
      endif
      if(ichr.lt.ichar(space)) chr = space
      return
      end


! ==================================================================== *
      subroutine UPPLOW (string)
! -------------------------------------------------------------------- *
!
! Comments: Converts a string from Upper-case to Lower-case
!
! author:F.Colonna
! date  :20 Fev 91
!--------------------------------------------------------------------- *
      implicit none
! input scalars:
      character*(*) string
      integer i
      integer first, last

      call sizstr (string, first, last)

      if(first.eq.0) first=1
      if(last.gt.0) then
         do i = first,last
         call CNVTLC (string(i:i))
         enddo
      else
      !  write(6,'(a)')' upplow-w: zero lenstr String '
      endif
      return
      end

! ==================================================================== *
      subroutine STRIM (String, strlen)
! ----------------------------------------------------------------------
!
! Comments: cuts out ending blanks, returns new length
!
! ----------------------------------------------------------------------
      implicit none
      character*(*) string
      integer strlen
! begin:
      strlen=len(string)
      if(strlen.gt.0) then
   10   continue
        if(strlen.ge.1.and.string(strlen:strlen).eq.' ') then
         strlen=strlen-1
         go to 10
        endif
      else
      write(6,'(9x,a)')' error_STRIM> zero length string '
      endif
      return
      end

      subroutine SIZSTR (string, first, last)
! -------------------------------------------------------------------- *
!
! Comments: returns first and last non blank character
!
! author:F.Colonna
! date:Tue Sep 22 15:25:28 GMT 1992
! -------------------------------------------------------------------- *
!
      implicit none

! i/o scalars  :
      character*(*) string
      integer first, last
! local   :
      integer is, lenstr
      logical found, max

! begin   :
      found = .false.

      lenstr = len(string)

      if(lenstr.eq.0) then
          write(6,'(a)') 'error_sizstr> zero length string'
          stop '   error_sizstr> zero length  string'
      endif

      max = .false.
      is = lenstr

! search not blank from end.
      do while (.not.found.and..not.max)
           last = is
           found = string(is:is) .ne. ' '
           is = is - 1
           max = is .eq. 0
      enddo

      if(max) then
          if(.not.found) then
!           write(6,'(a)') 'warning_sizstr> empty string'
            first = 1
            last  = 0
            string = ' '
          else
            first = 1
            last  = 1
          endif
          return
      endif

! search not blank from beginning.
      found = .false.
      max = .false.
      is = 1
      do while (.not.found.and..not.max)
           first = is
           found = string(is:is) .ne. ' '
           is = is +1
           max = is .eq. lenstr+1
      enddo

      if(max) then
          write(6,'(a)') 'error_sizstr> strange  string'
          stop ' error_sizstr> strange  string'
      endif

      if(last.lt.first) then
          write(6,'(a,i5)') 'error_sizstr> first non-blank position=',first, ' error_sizstr> last  non-blank position=',last
          stop ' error_sizstr> last < first'
      endif

! end of routine SIZSTR
      return
      end

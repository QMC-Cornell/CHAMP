c an input file is copied to input and an output buffers, which may be the same
c subsequently the output buffer is updated at the end the output buffer has to be
c disposed of properly

      subroutine int_read0(integer,name,iunit_in,iunit_out)
      include 'general.h'
      include 'io.h'
      pointer (input_buffer_pointer,input_lines)
      character*132 input_lines(*)
      character*132 substrings(MX_SUBSTRINGS),string,head_name,head_sub1
      character name*(*)
      parameter(NUMBER_FIELDS=2,NUMBER_AGREE=1)

      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif

      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'int_read0: nstrings. eq. 0'
        stop 'int_read0: nstrings. eq. 0'
      endif
      head_name=substrings(1)
      do l=1,num_input_lines
        call parse_string(input_lines(l),substrings,nstrings,MX_SUBSTRINGS)
        if(nstrings .ge. 2) then
          head_sub1=substrings(1)
          if(head_sub1 .eq. head_name) then
            read(substrings(2),*) integer
            call int_write_line(name,integer,1,string)
            call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
            return
          endif
        endif
      enddo
      call input_not_found(head_name)
      return
      end
c-----------------------------------------------------------------------

      subroutine real_read0(real,name,iunit_in,iunit_out)
      include 'general.h'
      include 'io.h'
      pointer (input_buffer_pointer,input_lines)
      character*132 input_lines(*)
      character*132 substrings(MX_SUBSTRINGS),string,head_name,head_sub1
      character name*(*)
      parameter(NUMBER_FIELDS=2,NUMBER_AGREE=1)

      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif

      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'real_read0: nstrings. eq. 0'
        stop 'real_read0: nstrings. eq. 0'
      endif
      head_name=substrings(1)
      do l=1,num_input_lines
        call parse_string(input_lines(l),substrings,nstrings,MX_SUBSTRINGS)
        if(nstrings .ge. NUMBER_FIELDS) then
          head_sub1=substrings(1)
          if(head_sub1 .eq. head_name) then
            read(substrings(2),*) real
            call real_write_line(name,NUMBER_AGREE-1,real,NUMBER_AGREE-1,string)
            call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
            return
          endif
        endif
      enddo
      call input_not_found(head_name)
      return
      end
c-----------------------------------------------------------------------

      subroutine int_read1(integer,isize,name,iunit_in,iunit_out)
      include 'general.h'
      include 'io.h'
      pointer (input_buffer_pointer,input_lines)
      character*132 input_lines(*)
      dimension integer(isize)
      character*132 substrings(MX_SUBSTRINGS),string,head_name,head_sub1
      character name*(*)
      parameter(NUMBER_FIELDS=3,NUMBER_AGREE=2)
      dimension iout(NUMBER_AGREE)
      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif
      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'int_read1: nstrings. eq. 0'
        stop 'int_read1: nstrings. eq. 0'
      endif
      head_name=substrings(1)
      ifound=0
      do l=1,num_input_lines
        call parse_string(input_lines(l),substrings,nstrings,MX_SUBSTRINGS)
        if(nstrings .ge. NUMBER_FIELDS) then
          head_sub1=substrings(1)
          if(head_sub1 .eq. head_name) then
            read(substrings(2),*) ndx
            if(ndx .ge.1 .and. ndx. le. isize) then
              ifound=ifound+1
              read(substrings(3),*) integer(ndx)
              iout(1)=ndx
              iout(2)=integer(ndx)
              call int_write_line(name,iout,2,string)
              call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
            endif
            if(ifound .eq. isize) return
          endif
        endif
      enddo
      call input_not_found1(head_name,ifound,isize)
      return
      end
c-----------------------------------------------------------------------

      subroutine real_read1(real,isize,name,iunit_in,iunit_out)
      include 'general.h'
      include 'io.h'
      pointer (input_buffer_pointer,input_lines)
      character*132 input_lines(*)
      dimension real(isize)
      character*132 substrings(MX_SUBSTRINGS),string,head_name,head_sub1
      character name*(*)
      parameter(NUMBER_FIELDS=3,NUMBER_AGREE=2)
      dimension iout(NUMBER_AGREE)
      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif
      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'real_read1: nstrings. eq. 0'
        stop 'real_read1: nstrings. eq. 0'
      endif
      head_name=substrings(1)
      ifound=0
      do l=1,num_input_lines
        call parse_string(input_lines(l),substrings,nstrings,MX_SUBSTRINGS)
        if(nstrings .ge. NUMBER_FIELDS) then
          head_sub1=substrings(1)
          if(head_sub1 .eq. head_name) then
            read(substrings(2),*) ndx
            if(ndx .ge.1 .and. ndx. le. isize) then
              ifound=ifound+1
              read(substrings(3),*) real(ndx)
              iout(1)=ndx
              call real_write_line(name,iout,real(ndx),1,string)
              call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
            endif
            if(ifound .eq. isize) return
          endif
        endif
      enddo
      call input_not_found1(head_name,ifound,isize)
      return
      end
c-----------------------------------------------------------------------

      subroutine int_read2(integer,iphys1,ifrom,ito,name,iunit_in,iunit_out)
      include 'general.h'
      include 'io.h'
      pointer (input_buffer_pointer,input_lines)
      character*132 input_lines(*)
      dimension integer(iphys1,*)
      character*132 substrings(MX_SUBSTRINGS),string,head_name,head_sub1
      character name*(*)
      parameter(NUMBER_FIELDS=4,NUMBER_AGREE=3)
      dimension iout(NUMBER_AGREE),ifrom(2),ito(2)
      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif
      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'int_read1: nstrings. eq. 0'
        stop 'int_read1: nstrings. eq. 0'
      endif
      head_name=substrings(1)
      ifound=0
      number=(ito(1)-ifrom(1)+1)*(ito(2)-ifrom(2)+1)
      do l=1,num_input_lines
        call parse_string(input_lines(l),substrings,nstrings,MX_SUBSTRINGS)
        if(nstrings .ge. NUMBER_FIELDS) then
          head_sub1=substrings(1)
          if(head_sub1 .eq. head_name) then
            read(substrings(2),*) ndx1
            read(substrings(3),*) ndx2
            if(ndx1 .ge. ifrom(1) .and. ndx1 .le. ito(1) .and.
     &        ndx2 .ge. ifrom(2) .and. ndx2 .le. ito(2)) then
              ifound=ifound+1
              read(substrings(4),*) integer(ndx1,ndx2)
              iout(1)=ndx1
              iout(2)=ndx2
              iout(3)=integer(ndx1,ndx2)
              call int_write_line(name,iout,3,string)
              call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
            endif
            if(ifound .eq. number) return
          endif
        endif
      enddo
      call input_not_found1(head_name,ifound,number)
      return
      end
c-----------------------------------------------------------------------

      subroutine real_read2(real,iphys1,ifrom,ito,name,iunit_in,iunit_out)
      include 'general.h'
      include 'io.h'
      pointer (input_buffer_pointer,input_lines)
      character*132 input_lines(*)
      dimension real(iphys1,*)
      character*132 substrings(MX_SUBSTRINGS),string,head_name,head_sub1
      character name*(*)
      parameter(NUMBER_FIELDS=4,NUMBER_AGREE=3)
      dimension iout(NUMBER_AGREE),ifrom(2),ito(2)
      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif
      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'real_read1: nstrings. eq. 0'
        stop 'real_read1: nstrings. eq. 0'
      endif
      head_name=substrings(1)
      ifound=0
      number=(ito(1)-ifrom(1)+1)*(ito(2)-ifrom(2)+1)
      do l=1,num_input_lines
        call parse_string(input_lines(l),substrings,nstrings,MX_SUBSTRINGS)
        if(nstrings .ge. NUMBER_FIELDS) then
          head_sub1=substrings(1)
          if(head_sub1 .eq. head_name) then
            read(substrings(2),*) ndx1
            read(substrings(3),*) ndx2
            if(ndx1 .ge. ifrom(1) .and. ndx1 .le. ito(1) .and.
     &        ndx2 .ge. ifrom(2) .and. ndx2 .le. ito(2)) then
               ifound=ifound+1
              read(substrings(4),*) real(ndx1,ndx2)
              iout(1)=ndx1
              iout(2)=ndx2
              call real_write_line(name,iout,real(ndx1,ndx2),2,string)
              call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
            endif
            if(ifound .eq. number) return
          endif
        endif
      enddo
      call input_not_found1(head_name,ifound,number)
      return
      end
c-----------------------------------------------------------------------

      block data io_setup_block
      include 'io.h'
      data io_setup/.false./
      end
c-----------------------------------------------------------------------

      subroutine int_write_line(string_in,iout,number_out,string_out)
      include 'general.h'
      include 'io.h'
      character*(*) string_in,string_out
      character*132 string_tmp
      character*4 ifield_length_char
      dimension iout(number_out)
      string_out=string_in
      do i=1,number_out
        ainteger=max(abs(iout(i)),1)
        ifield_length=log10(ainteger)+2
        if(ifield_length .ge. 10000) then
          print *,'int_read0: ifield_length .ge. 10000'
          stop 'int_read0: ifield_length .ge. 10000'
        endif
        write(ifield_length_char,'(i4)')ifield_length
        write(string_tmp,'(i'//ifield_length_char//','';'')') iout(i)
        call appendto(string_out,string_tmp)
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine int_write0(integer,name,iunit_in,iunit_out)
      include 'general.h'
      include 'io.h'
c      pointer (input_buffer_pointer,input_lines)
c      character*132 input_lines(*)
      character*132 substrings(MX_SUBSTRINGS),string
      character name*(*)
      dimension iout(1)
      parameter(NUMBER_FIELDS=2,NUMBER_AGREE=1)

      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif

      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'int_read0: nstrings. eq. 0'
        stop 'int_read0: nstrings. eq. 0'
      endif
      iout(1)=integer
      call int_write_line(name,iout,1,string)
      call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
      return
      end
c-----------------------------------------------------------------------

      subroutine real_write0(real,name,iunit_in,iunit_out)
      include 'general.h'
      include 'io.h'
c      pointer (input_buffer_pointer,input_lines)
c      character*132 input_lines(*)
      character*132 substrings(MX_SUBSTRINGS),string
      character name*(*)
      parameter(NUMBER_FIELDS=2,NUMBER_AGREE=1)

      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif

      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'real_read0: nstrings. eq. 0'
        stop 'real_read0: nstrings. eq. 0'
      endif
      call real_write_line(name,0,real,0,string)
      call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
      return
      end
c-----------------------------------------------------------------------

      subroutine int_write1(integers,ifrom,ito,name,iunit_in,iunit_out)
c write/overwrite integer array to output buffer
c This routine was not written for speed !!!!!!!
c name = identifier
c integers(1,...) = integer array to be written
c ifrom = first element to be written integers(ifrom)
c ito = last element to be written integers(ito)
c iunit_in,iunit_out = not used
      include 'general.h'
      include 'io.h'
      dimension integers(*)
      character name*(*)
      character*132 substrings(MX_SUBSTRINGS),string
      dimension iout(2)
      parameter(NUMBER_FIELDS=3,NUMBER_AGREE=2)

      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif
      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'int_read1: nstrings. eq. 0'
        stop 'int_read1: nstrings. eq. 0'
      endif
      do i=ifrom,ito
        iout(1)=i
        iout(2)=integers(i)
        call int_write_line(name,iout,2,string)
        call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine int_write2(integers,iphys1,ifrom,ito,name,iunit_in,iunit_out)
c write/overwrite integer array to output buffer
c This routine was not written for speed !!!!!!!
c integers(1:iphys1,*) = integer array to be written
c iphys1 = physical dimension of first index of integers
c ifrom(.) = first element to be written integers(ifrom(1),ifrom(2))
c ito = last element to be written integers(ito(1),ito(2))
c name = identifier
c iunit_in,iunit_out = not used
      include 'general.h'
      include 'io.h'
      dimension integers(iphys1,*)
      character name*(*)
      character*132 substrings(MX_SUBSTRINGS),string
      dimension ifrom(2),ito(2)
      dimension iout(3)
      parameter(NUMBER_FIELDS=4,NUMBER_AGREE=3)

      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif
      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'int_read1: nstrings. eq. 0'
        stop 'int_read1: nstrings. eq. 0'
      endif
      do i1=ifrom(1),ito(1)
        do i2=ifrom(2),ito(2)
          iout(1)=i1
          iout(2)=i2
          iout(3)=integers(i1,i2)
          call int_write_line(name,iout,3,string)
c         print *,'int2_write:',string
          call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine real_write1(reals,ifrom,ito,name,iunit_in,iunit_out)
c write/overwrite integer array to output buffer
c This routine was not written for speed !!!!!!!
c name = identifier
c reals(1,...) = real array to be written
c ifrom = first element to be written integers(ifrom)
c ito = last element to be written integers(ifrom)
c iunit_in,iunit_out = not used
      include 'general.h'
      include 'io.h'
      dimension reals(ito)
c     pointer (ioutput_buffer_pointer,output_lines)
      character name*(*)
c     character*132 output_lines(*)
      character*132 substrings(MX_SUBSTRINGS),string,head_name
      dimension iout(1)
      parameter(NUMBER_FIELDS=3,NUMBER_AGREE=2)

      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif

      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'real_read1: nstrings. eq. 0'
        stop 'real_read1: nstrings. eq. 0'
      endif
      head_name=substrings(1)
      do i=ifrom,ito
        iout(1)=i
        out=reals(i)
        call real_write_line(name,iout,out,1,string)
        call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine real_write2(reals,iphys1,ifrom,ito,name,iunit_in,iunit_out)
c write/overwrite integer array to output buffer
c This routine was not written for speed !!!!!!!
c reals(1:iphys1,*) = real array to be written
c iphys1 = physical dimension of first index of integers
c ifrom(.) = first element to be written integers(ifrom(1),ifrom(2))
c ito = last element to be written integers(ito(1),ito(2))
c name = identifier
c iunit_in,iunit_out = not used
      include 'general.h'
      include 'io.h'
      dimension reals(iphys1,*)
c pointer (ioutput_buffer_pointer,output_lines)
      character name*(*)
c character*132 output_lines(*)
      character*132 substrings(MX_SUBSTRINGS),string,head_name
      dimension ifrom(2),ito(2)
      dimension iout(2)
      parameter(NUMBER_FIELDS=4,NUMBER_AGREE=3)

      if(.not. io_setup) then
        call initialize_io(iunit_in,iunit_out,.false.)
      endif

      call parse_string(name,substrings,nstrings,MX_SUBSTRINGS)
      if(nstrings. eq. 0) then
        print *,'real_read1: nstrings. eq. 0'
        stop 'real_read1: nstrings. eq. 0'
      endif
      head_name=substrings(1)
      do i1=ifrom(1),ito(1)
        do i2=ifrom(2),ito(2)
          iout(1)=i1
          iout(2)=i2
          out=reals(i1,i2)
          call real_write_line(name,iout,out,2,string)
          call replace(string,NUMBER_FIELDS,NUMBER_AGREE)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine wrong_input(string1,string2)
      character*(*) string1,string2
      print *,'wrong_input: found '
      print '((a))',string1
      print *,'expected'
      print '((a))',string2
      stop 'wrong_input'
      end
c-----------------------------------------------------------------------

      subroutine input_not_found(string)
      character*(*) string
      print '((a))','input_not_found:',string
      stop 'input_not_found: see output'
      end
c-----------------------------------------------------------------------

      subroutine initialize_io(iunit_in,iunit_out,i2o)
c initialization
c iunit_in = read input form this unit
c iunit_out = write to this unit
c i2o = if(i2o) copy input to output
      include 'general.h'
      include 'io.h'
      pointer (input_buffer_pointer,input_lines),(ioutput_buffer_pointer,output_lines)
      character*132 input_lines(*),output_lines(*)
      logical i2o
      num_input_lines=0
 10   continue
      num_input_lines=num_input_lines+1
      read(iunit_in,'((a))',end=20)
      goto 10

 20   continue
      num_input_lines=num_input_lines-1
      mx_input_lines=num_input_lines+100
      input_buffer_pointer=malloc(mx_input_lines*LENGTH_STRINGS)
      if(input_buffer_pointer .le. 0) then
         print '((a))','initialize_io: malloc failed for input'
         stop 'initialize_io: malloc failed for input'
      endif
      mx_output_lines=mx_input_lines
      ioutput_buffer_pointer=malloc(mx_output_lines*LENGTH_STRINGS)
      if(input_buffer_pointer .le. 0) then
         print '((a))','initialize_io: malloc failed for output'
         stop 'initialize_io: malloc failed for output'
      endif
      num_input_lines=0
      rewind(iunit_in)

 30   continue
      num_input_lines=num_input_lines+1
      read(iunit_in,'((a))',end=40) input_lines(num_input_lines)
      if(i2o) output_lines(num_input_lines)=input_lines(num_input_lines)
      goto 30
 40   continue
      num_input_lines=num_input_lines-1
      if(i2o) then
        num_output_lines=num_input_lines
      else
        num_output_lines=0
      endif
      io_setup=.true.
      return
      end
c-----------------------------------------------------------------------

      subroutine input_not_found1(head_name,ifound,isize)
      character*(*) head_name
      print '(2(a))','input_not_found1: ',head_name(1:index(head_name,' '))
      print '(''found '',(i),'' of '',(i))',ifound,isize
      return
      end
c-----------------------------------------------------------------------

      subroutine expand_buffer(ibuffer_pointer)
      include 'general.h'
      include 'io.h'
      pointer (ibuffer_pointer,lines),(new_buffer_pointer,new_lines)
      character*132 lines(*),new_lines(*)
      if(ibuffer_pointer .eq. input_buffer_pointer) then
        mx_input_lines=mx_input_lines+IEXPAND_BUFFER
        num_lines=num_input_lines
        mx_lines=mx_input_lines
      else
        mx_output_lines=mx_output_lines+IEXPAND_BUFFER
        num_lines=num_output_lines
        mx_lines=mx_output_lines
      endif
      new_buffer_pointer=malloc(mx_lines*LENGTH_STRINGS)
      if(new_buffer_pointer .le. 0) then
        print '((a))','expand_buffer: malloc failed'
        stop 'expand_buffer: malloc failed'
      endif
      do i=1,num_lines
        new_lines(i)=lines(i)
      enddo
      call free(ibuffer_pointer)
      if(ibuffer_pointer .eq. input_buffer_pointer) then
        input_buffer_pointer=new_buffer_pointer
      elseif(ibuffer_pointer .eq. ioutput_buffer_pointer) then
        ioutput_buffer_pointer=new_buffer_pointer
      else
        print *,'expand_buffer: wrong ibuffer_pointer'
        stop 'expand_buffer: wrong ibuffer_pointer'
      endif
      return
      end
c-----------------------------------------------------------------------

      subroutine real_write_line(string_in,iout,out,idim,string_out)
      include 'general.h'
      include 'io.h'
      character*(*) string_in,string_out
      character*132 string_tmp
      character*4 ifield_length_char
      dimension iout(idim)
      string_out=string_in
      do i=1,idim
        ainteger=max(abs(iout(i)),1)
        ifield_length=log10(ainteger)+2
        if(ifield_length .ge. 10000) then
         print *,'int_read0: ifield_length .ge. 10000'
          stop 'int_read0: ifield_length .ge. 10000'
        endif
        write(ifield_length_char,'(i4)')ifield_length
        write(string_tmp,'(i'//ifield_length_char//','';'')') iout(i)
        call appendto(string_out,string_tmp)
      enddo
      write(string_tmp,'(g23.15,(a))')out,TERMINATOR
      call appendto(string_out,string_tmp)
      return
      end
c-----------------------------------------------------------------------

      subroutine appendto(string1,string2)
      include 'general.h'
      include 'io.h'
      character*(*) string1,string2
      do i=1,len(string1)
        iend1=i-1
        if(string1(i:i) .eq. TERMINATOR) goto 10
      enddo
      print *,'appendto: end of string1 (',TERMINATOR,') not found in',string1
      call abort('appendto: end of string1 not found')
 10   continue
      do i=1,len(string2)
        iend2=i-1
        if(string2(i:i) .eq. TERMINATOR) goto 20
      enddo
      print *,'appendto: end of string2 (',TERMINATOR,') not found',string2
      stop 'appendto: end of string2 not found'
 20   continue
      if(iend1+iend2+2 .gt. LENGTH_STRINGS) then
        print *,'appendto: iend1+iend2 .gt. LENGTH_STRINGS',iend1,iend2,LENGTH_STRINGS
        call abort('appendto: iend1+iend2 .gt. LENGTH_STRINGS')
      endif
      string1=string1(1:iend1)//SEPARATOR//string2(1:iend2)//TERMINATOR
      return
      end
c-----------------------------------------------------------------------

      subroutine print_buffer(ibuffer_pointer,iunit_out)
      include 'general.h'
      include 'io.h'
      external precedes
      pointer (ibuffer_pointer,lines)
      character*132 lines(*)
      if(ibuffer_pointer .eq. input_buffer_pointer) then
        num_lines=num_input_lines
      elseif(ibuffer_pointer .eq. ioutput_buffer_pointer) then
        num_lines=num_output_lines
      else
        stop 'print_buffer: invalid pointer value'
      endif
      call sort_strings(num_lines,lines,1,precedes)
      do i=1,num_lines
        length=index(lines(i),TERMINATOR)
        if(length .gt. LENGTH_STRINGS) then
          print *,'print_buffer: length .gt. LENGTH_STRINGS',length,LENGTH_STRINGS
          print *,'print_buffer: lines(i)',lines(i)
          call abort('print_buffer: length .gt. LENGTH_STRINGS')
c         stop 'print_buffer: length .gt. LENGTH_STRINGS'
        endif
        write(iunit_out,'((a))') lines(i)(1:index(lines(i),TERMINATOR))
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine replace(string,number_fields,number_agree)
c replace string in output buffer or add if not present
c string = string to be replaced
c number_fields = total number of fields
c number_agree = number of fields of agreement required for replacement
c note: this routine only replaces the string that occurs first
      include 'general.h'
      include 'io.h'
      character*(*) string
      character*132 substrings1(MX_SUBSTRINGS),substrings2(MX_SUBSTRINGS)
      logical agree
      character*132 output_lines(*)
      pointer (ioutput_buffer_pointer,output_lines)
      call parse_string(string,substrings1,n1,MX_SUBSTRINGS)
      if(n1 .ge. number_agree) then
        do i=1,num_output_lines
          call parse_string(output_lines(i),substrings2,n2,MX_SUBSTRINGS)
          if(n2 .ge. number_agree) then
            do j=1,number_agree
              if(substrings1(j) .ne. substrings2(j)) goto 10
            enddo
            line=i
            goto 20
          endif
 10       continue
        enddo
        if(num_output_lines .eq. mx_output_lines) call expand_buffer(ioutput_buffer_pointer)
        num_output_lines=num_output_lines+1
        line=num_output_lines
      endif
 20   continue
      length=index(string,TERMINATOR)
      if(length .gt. LENGTH_STRINGS) then
        print *,'replace: length .gt. LENGTH_STRINGS',length,LENGTH_STRINGS
        stop 'replace: length .gt. LENGTH_STRINGS'
        print *,'string =',string
      endif
      if(length .lt. 0) then
        print *,'replace: length .lt. 0'
        stop 'replace: length .lt. 0'
        print *,'string =',string
      endif
      if(len(string(1:length)) .gt. len(output_lines(line))) then
        print *,'len(string(1:length)) .gt. len(output_lines(line))',
     &    output_lines(line),' ',string(1:length)
      endif
      if(line .gt. mx_output_lines) then
        print *,'line .gt. mx_output_lines',line,mx_output_lines
      endif

      output_lines(line)=string(1:length)
      return
      end
c-----------------------------------------------------------------------

      logical function precedes(string1,string2)
c define ordering on strings
      include 'io.h'
c      include 'general.h'
c      include 'parameters.h'
      character*132 substrings1(MX_SUBSTRINGS),substrings2(MX_SUBSTRINGS)
      character*(*) string1,string2
      call parse_string(string1,substrings1,nstrings1,MX_SUBSTRINGS)
      call parse_string(string2,substrings2,nstrings2,MX_SUBSTRINGS)
      if(nstrings1 .lt. nstrings2) then
        precedes=.true.
        return
      elseif(nstrings1 .gt. nstrings2) then
        precedes=.false.
        return
      else !same number of substrings
        if(llt(substrings1(1),substrings2(1))) then
          precedes=.true.
          return
        elseif(lgt(substrings1(1),substrings2(1))) then
          precedes=.false.
          return
        else !same first substring
          do i=2,nstrings1-1
            read(substrings1(i),*)index1
            read(substrings2(i),*)index2
            if(index1 .lt. index2) then
              precedes=.true.
              return
            elseif(index1 .gt. index2) then
              precedes=.false.
              return
            endif
          enddo
        endif
      endif
      print '((a))','precedes: list contains same item more than once'
      print '((a))',string1
      print '((a))',string2
      stop 'precedes: list contains same item more than once'
      end
c-----------------------------------------------------------------------

      subroutine parse_string(string,substrings,n,max_substrings)
c chop string into substrings
c string = string (in)
c substrings = substrings (out)
c n = their number
c max_substrings = maximal number of substrings allowed
      include 'general.h'
      include 'io.h'
      character*(*) string,substrings(*)
      length=index(string,TERMINATOR)-1
      n=0
      if(length .lt. 0) length=len(string)
      ibegin=1
 10   continue
      do while(ibegin .le. length .and. string(ibegin:ibegin) .eq. SEPARATOR)
        ibegin=ibegin+1
      enddo
      if(ibegin .gt. length) return
      iend=ibegin+1
      do while(iend .le. length .and. string(iend:iend) .ne. SEPARATOR)
        iend=iend+1
      enddo
      iend=iend-1
      n=n+1
      if(n .gt. max_substrings) then
        print '((a))','parse_string: n .gt. max_substrings'
        print *,'string: ',string(1:min(length,LENGTH_STRINGS))
        call abort('parse_string: n .gt. max_substrings')
      endif
      substrings(n)=string(ibegin:iend)//TERMINATOR
      ibegin=iend+1
      goto 10
      end


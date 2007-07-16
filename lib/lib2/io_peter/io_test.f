      include 'general.h'
      include 'io.h'
      dimension iarray(10),a(10),b(10),c(10)
      dimension jarray(10,2)
      dimension array(10,2)
      dimension ifrom(2),ito(2)
      character*20 string, substrings(20)
      data iarray/10*0/, a/10*0/, b/10*0/,c/10*0/,jarray/20*0/,array/20*0/
      iunit_in=2
      iunit_out=1
      call int_read0(iseed,'iseed;',iunit_in,iunit_out)
      call int_read0(jseed,'jseed;',iunit_in,iunit_out)
      call int_write0(jseed+1,'jseed;',iunit_in,iunit_out)
      call real_read0(seed,'seed;',iunit_in,iunit_out)
      call real_write0(seed+1,'seed;',iunit_in,iunit_out)
      call real_write0(seed+2,'seed;',iunit_in,iunit_out)
      number=4
      call int_read1(iarray,number,'iarray;',iunit_in,iunit_out)
      print *,'iarray:',(iarray(j),j=1,number)
      call real_read0(a,'a;',iunit_in,iunit_out)
      call real_read0(c,'c;',iunit_in,iunit_out)
      call real_read1(b,number,'b;',iunit_in,iunit_out)
      print *,'b',(b(i),i=1,number)
      do i=1,number
        b(i)=-b(i)**2
        iarray(i)=-iarray(i)**2
      enddo
      call real_write1(b,1,number,'b;',iunit_in,iunit_out)
      call int_write1(iarray,1,number,'iarray;',iunit_in,iunit_out)
      do i=1,number
        b(i)=-b(i)
        iarray(i)=-iarray(i)
      enddo
      call real_write1(b,1,number,'b;',iunit_in,iunit_out)
      call int_write1(iarray,1,number,'iarray;',iunit_in,iunit_out)
      ifrom(1)=1
      ifrom(2)=1
      ito(1)=number
      ito(2)=2
      call int_read2(jarray,10,ifrom,ito,'jarray;',iunit_in,iunit_out)
      do i=ifrom(1),ito(1)
        do j=ifrom(2),ito(2)
          jarray(i,j)=jarray(i,j)+10
        enddo
      enddo
      call int_write2(jarray,10,ifrom,ito,'jarray;',iunit_in,iunit_out)
      call real_read2(array,10,ifrom,ito,'array;',iunit_in,iunit_out)
      do i=ifrom(1),ito(1)
        do j=ifrom(2),ito(2)
          array(i,j)=array(i,j)+10
        enddo
      enddo
      call real_write2(array,10,ifrom,ito,'array;',iunit_in,iunit_out)
      call print_buffer(ioutput_buffer_pointer,iunit_out)
      end

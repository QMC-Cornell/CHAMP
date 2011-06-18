      subroutine den2dwrt(passes)
c Written by A.D.Guclu, modified by Cyrus Umrigar for MPI
c routine to print out 2d-density related quantities
c called by finwrt from vmc,dmc,dmc_elec

      use mpi_mod
      use dets_mod
      use contr3_mod
      use pairden_mod
      use fourier_mod
      implicit real*8(a-h,o-z)
      character*20 file1,file2,file3,file4,file5,file6

      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys

c verify the normalization later...
      delx=1/delxi
      if(icoosys.eq.1) then
        del1=delx
        del2=delx
        nax1=NAX
        nax2=NAX
      else
        del1=1/delradi
        del2=1/delti
        nax1=nmeshr
        nax2=nmesht
      endif
      term=1/(passes*del1*del2)

      if(ifixe.gt.0) then          ! fixed electron pair-density
        if(ifixe.le.nup) then
          file1='pairden_ut'
          file2='pairden_ud'
          file3='pairden_uu'
         else
          file1='pairden_dt'
          file2='pairden_dd'
          file3='pairden_du'
        endif

        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
         else
          open(41,file='/dev/null',status='unknown')
          open(42,file='/dev/null',status='unknown')
          open(43,file='/dev/null',status='unknown')
        endif

c verify the normalization later...
        do in1=-nax1,nax1
          do in2=-nax2,nax2
            write(41,'(2g19.8,g19.8)') in1*del1,in2*del2,den2d_t(in1,in2)*term
            write(42,'(2g19.8,g19.8)') in1*del1,in2*del2,den2d_d(in1,in2)*term
            write(43,'(2g19.8,g19.8)') in1*del1,in2*del2,den2d_u(in1,in2)*term
          enddo
c following spaces are for gnuplot convention:
          write(41,*)
          write(42,*)
          write(43,*)
        enddo
        close(41)
        close(42)
        close(43)
      endif

      if(ifixe.eq.-1 .or. ifixe.eq.-3) then          ! 2d density
        if(index(mode,'vmc').ne.0) then
          file1='den2d_t_vmc'
          file2='den2d_d_vmc'
          file3='den2d_u_vmc'
          file4='pot_t_vmc'
          file5='pot_d_vmc'
          file6='pot_u_vmc'
         else
          file1='den2d_t_dmc'
          file2='den2d_d_dmc'
          file3='den2d_u_dmc'
          file4='pot_t_dmc'
          file5='pot_d_dmc'
          file6='pot_u_dmc'
        endif
        
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
          open(44,file=file4,status='unknown')
          open(45,file=file5,status='unknown')
          open(46,file=file6,status='unknown')
         else
          open(41,file='/dev/null',status='unknown')
          open(42,file='/dev/null',status='unknown')
          open(43,file='/dev/null',status='unknown')
          open(44,file='/dev/null',status='unknown')
          open(45,file='/dev/null',status='unknown')
          open(46,file='/dev/null',status='unknown')
        endif

c verify the normalization later...
        do in1=-nax1,nax1
          do in2=-nax2,nax2
            write(41,'(2g19.8,g19.8)') in1*del1,in2*del2,den2d_t(in1,in2)*term
            write(42,'(2g19.8,g19.8)') in1*del1,in2*del2,den2d_d(in1,in2)*term
            write(43,'(2g19.8,g19.8)') in1*del1,in2*del2,den2d_u(in1,in2)*term
            write(44,'(2g19.8,g19.8)') in1*del1,in2*del2,pot_ee2d_t(in1,in2)/passes
            write(45,'(2g19.8,g19.8)') in1*del1,in2*del2,pot_ee2d_d(in1,in2)/passes
            write(46,'(2g19.8,g19.8)') in1*del1,in2*del2,pot_ee2d_u(in1,in2)/passes
          enddo
          write(41,*)
          write(42,*)
          write(43,*)
          write(44,*)
          write(45,*)
          write(46,*)
        enddo
        close(41)
        close(42)
        close(43)
        close(44)
        close(45)
        close(46)
      endif

      if(ifixe.le.-2) then          ! full pair-density

c up electron:
        if(nup.gt.0) then
          file1='pairden_ut'
          file2='pairden_ud'
          file3='pairden_uu'

          if(idtask.eq.0) then
            open(41,file=file1,status='unknown')
            open(42,file=file2,status='unknown')
            open(43,file=file3,status='unknown')
           else
            open(41,file='/dev/null',status='unknown')
            open(42,file='/dev/null',status='unknown')
            open(43,file='/dev/null',status='unknown')
          endif
          do in0=0,NAX
            r0=in0*delx
            if(in0.ge.nint(xfix(1)*delxi) .and. in0.le.nint(xfix(2)*delxi)) then
              write(41,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              write(42,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              write(43,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              do in1=-NAX,NAX
                do in2=-NAX,NAX
                  write(41,'(2g19.8,g19.8)') in1*delx,in2*delx,xx0probut(in0,in1,in2)*term
                  write(42,'(2g19.8,g19.8)') in1*delx,in2*delx,xx0probud(in0,in1,in2)*term
                  write(43,'(2g19.8,g19.8)') in1*delx,in2*delx,xx0probuu(in0,in1,in2)*term
                enddo
                write(41,*)
                write(42,*)
                write(43,*)
              enddo
            endif
          enddo
          close(41)
          close(42)
          close(43)
        endif

c down electron:
        if(ndn.gt.0) then
          file1='pairden_dt'
          file2='pairden_dd'
          file3='pairden_du'

          if(idtask.eq.0) then
            open(41,file=file1,status='unknown')
            open(42,file=file2,status='unknown')
            open(43,file=file3,status='unknown')
           else
            open(41,file='/dev/null',status='unknown')
            open(42,file='/dev/null',status='unknown')
            open(43,file='/dev/null',status='unknown')
          endif

          do in0=0,NAX
            r0=in0*delx
            if(in0.ge.nint(xfix(1)*delxi) .and. in0.le.nint(xfix(2)*delxi)) then
              write(41,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              write(42,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              write(43,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              do in1=-NAX,NAX
                do in2=-NAX,NAX
                  write(41,'(2g19.8,g19.8)') in1*delx,in2*delx,xx0probdt(in0,in1,in2)*term
                  write(42,'(2g19.8,g19.8)') in1*delx,in2*delx,xx0probdd(in0,in1,in2)*term
                  write(43,'(2g19.8,g19.8)') in1*delx,in2*delx,xx0probdu(in0,in1,in2)*term
                enddo
                write(41,*)
                write(42,*)
                write(43,*)
              enddo
            endif
          enddo
          close(41)
          close(42)
          close(43)
        endif

      endif

      if(ifourier.eq.1 .or. ifourier.eq.3) then          ! 2d fourier t. of the density
        file1='fourierrk_t'
        file2='fourierrk_d'
        file3='fourierrk_u'
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
         else
          open(41,file='/dev/null',status='unknown')
          open(42,file='/dev/null',status='unknown')
          open(43,file='/dev/null',status='unknown')
        endif

c verify the normalization later...
        term=1/(passes*delx)
        do in1=1,NAX
          ri=(in1*1.d0-0.5d0)*delx
          ri2=ri*ri
          do in2=0,NAK1
            fk=delk1*in2
            write(41,'(2g19.8,g19.8)') ri,fk,fourierrk_t(in1,in2)*term/ri2
            write(42,'(2g19.8,g19.8)') ri,fk,fourierrk_d(in1,in2)*term/ri2
            write(43,'(2g19.8,g19.8)') ri,fk,fourierrk_u(in1,in2)*term/ri2
          enddo
          write(41,*)
          write(42,*)
          write(43,*)
        enddo
        close(41)
        close(42)
        close(43)
      endif

      if(ifourier.eq.2 .or. ifourier.eq.3) then
        file1='fourierkk_t'
        file2='fourierkk_d'
        file3='fourierkk_u'
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
         else
          open(41,file='/dev/null',status='unknown')
          open(42,file='/dev/null',status='unknown')
          open(43,file='/dev/null',status='unknown')
        endif

c verify the normalization later...
        do in1=-NAK2,NAK2
          do in2=-NAK2,NAK2
            write(41,'(2g19.8,g19.8)') in1*delk2,in2*delk2,fourierkk_t(in1,in2)/passes
            write(42,'(2g19.8,g19.8)') in1*delk2,in2*delk2,fourierkk_d(in1,in2)/passes
            write(43,'(2g19.8,g19.8)') in1*delk2,in2*delk2,fourierkk_u(in1,in2)/passes
          enddo
          write(41,*)
          write(42,*)
          write(43,*)
        enddo
        close(41)
        close(42)
        close(43)
      endif

      return
      end

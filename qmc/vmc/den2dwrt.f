      subroutine den2dwrt(passes)
c Written by A.D.Guclu, modified by Cyrus Umrigar for MPI
c routine to print out 2d-density related quantities
c called by finwrt from vmc,dmc,dmc_elec

      use mpi_mod
      use dets_mod
      use contr3_mod
      use pairden_mod
      use fourier_mod
      use contrl_per_mod
      use zigzag_mod
      use const_mod
      implicit real*8(a-h,o-z)
      character*20 file1,file2,file3,file4,file5,file6

      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
c verify the normalization later...
c      delx=1/delxi    ! doesn't work now that delxi is an array
      if(icoosys.eq.1) then
        del1=1/delxi(1)
        del2=1/delxi(2)
        dely=del2 !used for pair density - vary transverse coord ("y")
        delxt=del1  !used for zz pair density
        nax1=NAX
        nax2=NAX
      else
        del1=1/delradi
        del2=1/delti
        dely=del1 !used for pair density - vary transverse coord ("r")
        delxt=del2
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
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
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
cc        open(41,file='/dev/null',status='unknown')
cc        open(42,file='/dev/null',status='unknown')
cc        open(43,file='/dev/null',status='unknown')
c         open(41,file='junk41',status='unknown')
c         open(42,file='junk42',status='unknown')
c         open(43,file='junk43',status='unknown')
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
          open(44,status='scratch')
          open(45,status='scratch')
          open(46,status='scratch')
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
          if(index(mode,'vmc').ne.0) then
            file1='pairden_ut_vmc'
            file2='pairden_ud_vmc'
            file3='pairden_uu_vmc'
          else
            file1='pairden_ut_dmc'
            file2='pairden_ud_dmc'
            file3='pairden_uu_dmc'
          endif

          if(idtask.eq.0) then
            open(41,file=file1,status='unknown')
            open(42,file=file2,status='unknown')
            open(43,file=file3,status='unknown')
           else
            open(41,status='scratch')
            open(42,status='scratch')
            open(43,status='scratch')
          endif
          do in0=0,NAX
            if(icoosys.eq.1) then
              r0=in0*dely
              if(in0.lt.nint(xfix(1)/dely) .or. in0.gt.nint(xfix(2)/dely)) cycle
              write(41,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              write(42,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              write(43,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
            else
              r0=in0*dely + xfix(1)
              if(r0.gt.xfix(2)) exit
              write(41,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0-rmean
              write(42,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0-rmean
              write(43,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0-rmean
            endif
            do in1=-NAX,NAX
              do in2=-NAX,NAX
                write(41,'(2g19.8,g19.8)') in1*del1,in2*del2,xx0probut(in0,in1,in2)*term
                write(42,'(2g19.8,g19.8)') in1*del1,in2*del2,xx0probud(in0,in1,in2)*term
                write(43,'(2g19.8,g19.8)') in1*del1,in2*del2,xx0probuu(in0,in1,in2)*term
              enddo
              write(41,*)
              write(42,*)
              write(43,*)
            enddo
          enddo
          close(41)
          close(42)
          close(43)
        endif

c down electron:
        if(ndn.gt.0) then
          if(index(mode,'vmc').ne.0) then
            file1='pairden_dt_vmc'
            file2='pairden_dd_vmc'
            file3='pairden_du_vmc'
          else
            file1='pairden_dt_dmc'
            file2='pairden_dd_dmc'
            file3='pairden_du_dmc'
          endif

          if(idtask.eq.0) then
            open(41,file=file1,status='unknown')
            open(42,file=file2,status='unknown')
            open(43,file=file3,status='unknown')
           else
            open(41,status='scratch')
            open(42,status='scratch')
            open(43,status='scratch')
          endif

          do in0=0,NAX
            if(icoosys.eq.1) then
              r0=in0*dely
              if(in0.lt.nint(xfix(1)/dely) .or. in0.gt.nint(xfix(2)/dely)) cycle
              write(41,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              write(42,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
              write(43,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0
            else
              r0=in0*dely + xfix(1)
              if(r0.gt.xfix(2)) exit
              write(41,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0-rmean
              write(42,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0-rmean
              write(43,'(''# Grid point:'',i4,''  r0 ='',g19.8)') in0,r0-rmean
            endif
            do in1=-NAX,NAX
              do in2=-NAX,NAX
                write(41,'(2g19.8,g19.8)') in1*del1,in2*del2,xx0probdt(in0,in1,in2)*term
                write(42,'(2g19.8,g19.8)') in1*del1,in2*del2,xx0probdd(in0,in1,in2)*term
                write(43,'(2g19.8,g19.8)') in1*del1,in2*del2,xx0probdu(in0,in1,in2)*term
              enddo
              write(41,*)
              write(42,*)
              write(43,*)
            enddo
          enddo
          close(41)
          close(42)
          close(43)
        endif

      endif

      if(ifourier.eq.1 .or. ifourier.eq.3) then          ! 2d fourier t. of the density
        if(index(mode,'vmc').ne.0) then
          file1='fourierrk_t_vmc'
          file2='fourierrk_d_vmc'
          file3='fourierrk_u_vmc'
        else
          file1='fourierrk_t_dmc'
          file2='fourierrk_d_dmc'
          file3='fourierrk_u_dmc'
        endif
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
         else
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
        endif

c verify the normalization later...
        if(iperiodic.eq.1) then
          term = 1/(passes*dely)
          naxmin = -NAX
          naxmax = NAX
        elseif(rring.eq.0.d0) then
          term=1/(passes*del1)
          naxmin = 1
          naxmax = NAX
        else ! for rings, use same r_min and r_max as density
          term=delradi/passes
          naxmin = -nmeshr
          naxmax = nmeshr
        endif
        do in1=naxmin,naxmax
          if(iperiodic.eq.1) then
            ri = in1*dely
            ri2 = 1.0  ! no need to normalize by 1/r^2 for wires
          elseif(rring.eq.0.d0) then
            ri=(in1*1.d0-0.5d0)*delx
            ri2 = ri*ri
          else
            ri = in1/delradi + rmean
            ri2 = ri*ri
          endif
          do in2=0,nmeshk1
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
        if(index(mode,'vmc').ne.0) then
          file1='fourierkk_t_vmc'
          file2='fourierkk_d_vmc'
          file3='fourierkk_u_vmc'
        else
          file1='fourierkk_t_dmc'
          file2='fourierkk_d_dmc'
          file3='fourierkk_u_dmc'
        endif
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(43,file=file3,status='unknown')
         else
          open(41,status='scratch')
          open(42,status='scratch')
          open(43,status='scratch')
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
      
      if(izigzag.gt.0) then
        if(index(mode,'vmc').ne.0) then
          file1='zzcorr_vmc'
          file2='zzcorrij_vmc'
          file5='znncorr_vmc'
          file6='zn2ncorr_vmc'
        else
          file1='zzcorr_dmc'
          file2='zzcorrij_dmc'
          file5='znncorr_dmc'
          file6='zn2ncorr_dmc'
        endif
        if(idtask.eq.0) then
          open(41,file=file1,status='unknown')
          open(42,file=file2,status='unknown')
          open(45,file=file5,status='unknown')
          open(46,file=file6,status='unknown')
         else
          open(41,status='scratch')
          open(42,status='scratch')
          open(45,status='scratch')
          open(46,status='scratch')
        endif
        do in2=0,nax2
          write(41,'(g19.8,g19.8)') in2*delxt,zzcorr(in2)/passes
          write(45,'(g19.8,g19.8)') in2*delxt*10./dble(nelec),znncorr(in2)/passes
          write(46,'(g19.8,g19.8)') in2*delxt*10./dble(nelec),zn2ncorr(in2)/passes
        enddo
        do ine = 0,nelec-1
          write(42,'(i8,g19.8)') ine,zzcorrij(ine)/passes
        enddo
        
        if(izigzag.eq.2) then
          if(index(mode,'vmc').ne.0) then
            file3='zzpairden_t_vmc'
            file4='zzpairdenij_t_vmc'
          else
            file3='zzpairden_t_dmc'
            file4='zzpairdenij_t_dmc'
          endif
          if(idtask.eq.0) then
            open(43,file=file3,status='unknown')
            open(44,file=file4,status='unknown')
           else
            open(43,status='scratch')
            open(44,status='scratch')
          endif
          
          do in1=-nax1,nax1
            do in2 = -nax2,nax2
              write(43,'(2g19.8,g19.8)') in1*zzdelyr,in2*delxt,zzpairden_t(in1,in2)*term
            enddo
            do ine = 0,nelec-1
              write(44,'(g19.8,i8,g19.8)') in1*zzdelyr,ine,zzpairdenij_t(in1,ine)/(passes*del1)
            enddo
            write(43,*)
            write(44,*)
          enddo
          close(43)
          close(44)
        endif
        close(41)
        close(42)
        close(45)
        close(46)
      endif

      return
      end

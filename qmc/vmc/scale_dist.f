      subroutine set_scale_dist(ipr,iw)
c Written by Cyrus Umrigar

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      parameter (half=0.5d0,third=1.d0/3.d0)

      common /dim/ ndim
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /bparm/ nspin2b,nocuspb

c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 have zero 2nd and 3rd derivatives at 0
c isc = 6,7 are exponential and power law scalings resp. that have
c       zero, value and 1st 2 derivatives at cut_jas.
c isc = 12,14, are  are similar 2,4 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c isc = 16,17, iflag=1,2 are the same as 6,7 resp.
c       This is used to scale r_en in J_en and r_ee in J_ee for solids.
c isc = 16,17, iflag=3 are similar 6,7 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_en in J_een for solids.
c isc = 16,17, iflag=4 have infinite range and are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_ee in J_een for solids.
c isc = 8,10 are the same as 2,4 resp. but multiplied by scalek so that at very
c       small scalek the normalizations of the basis functions while optimizing
c       is not many orders of magnitude different.  No longer needed because I
c       now use the diagonal of the overlap to normalize the basis functions.

c Evaluate constants that need to be reset if scalek is being varied.
c Warning: At present we are assuming that the same scalek is used
c for primary and secondary wavefns.  Otherwise c1_jas6i, c1_jas6, c2_jas6
c should be dimensioned to MWF
c Note val_cutjas is set isc=6,7 when isc=16,17 because isc=6,7 scalings are
c used for J_en and J_ee when isc=16,17.
      if(isc.eq.6 .or. isc.eq.16) then

        scalek2=scalek(iw)*scalek(iw)
        term=third*scalek(iw)*cutjas_en
        val_cutjas_en=exp(-term)
        cutjasi_en=1/cutjas_en
        c1_jas6i_en=1-val_cutjas_en
        c1_jas6_en(iw)=1/c1_jas6i_en
        c2_jas6_en(iw)=val_cutjas_en*c1_jas6_en(iw)
        asymp_r_en(iw)=c1_jas6i_en/scalek(iw)
        dasymp_r_en(iw)=((term+1)*val_cutjas_en-1)/scalek2
        d2asymp_r_en(iw)=(2-(term*term+2*term+2)*val_cutjas_en)/(scalek2*scalek(iw))
        if(ipr.eq.1) write(6,'(''cutjas_en,c1_jas6_en,c2_jas6_en='',f10.6,9d12.5)')
     &       cutjas_en,c1_jas6_en(iw),c2_jas6_en(iw)

        term=third*scalek(1)*cutjas_ee
        val_cutjas_ee=exp(-term)
        cutjasi_ee=1/cutjas_ee
        c1_jas6i_ee=1-val_cutjas_ee
        c1_jas6_ee(iw)=1/c1_jas6i_ee
        c2_jas6_ee(iw)=val_cutjas_ee*c1_jas6_ee(iw)
        asymp_r_ee(iw)=c1_jas6i_ee/scalek(iw)
        dasymp_r_ee(iw)=((term+1)*val_cutjas_ee-1)/scalek2
        d2asymp_r_ee(iw)=(2-(term*term+2*term+2)*val_cutjas_ee)/(scalek2*scalek(iw))
        if(ipr.eq.1) write(6,'(''cutjas_ee,c1_jas6_ee,c2_jas6_ee='',f10.6,9d12.5)')
     &       cutjas_ee,c1_jas6_ee(iw),c2_jas6_ee(iw)

       elseif(isc.eq.7 .or. isc.eq.17) then

        val_cutjas_en=1/(1+third*scalek(iw)*cutjas_en)
        cutjasi_en=1/cutjas_en
        c1_jas6i_en=1-val_cutjas_en
        c1_jas6_en(iw)=1/c1_jas6i_en
        c2_jas6_en(iw)=val_cutjas_en*c1_jas6_en(iw)
        asymp_r_en(iw)=c1_jas6i_en/scalek(iw)
        dasymp_r_en(iw)=-asymp_r_en(iw)*asymp_r_en(iw)
        d2asymp_r_en(iw)=-2*asymp_r_en(iw)*dasymp_r_en(iw)
        if(ipr.eq.1) write(6,'(''cutjas_en,c1_jas6_en,c2_jas6_en='',f10.6,9d12.5)')
     &       cutjas_en,c1_jas6_en(iw),c2_jas6_en(iw)

        val_cutjas_ee=1/(1+third*scalek(iw)*cutjas_ee)
        cutjasi_ee=1/cutjas_ee
        c1_jas6i_ee=1-val_cutjas_ee
        c1_jas6_ee(iw)=1/c1_jas6i_ee
        c2_jas6_ee(iw)=val_cutjas_ee*c1_jas6_ee(iw)
        asymp_r_ee(iw)=c1_jas6i_ee/scalek(iw)
        dasymp_r_ee(iw)=-asymp_r_ee(iw)*asymp_r_ee(iw)
        d2asymp_r_ee(iw)=-2*asymp_r_ee(iw)*dasymp_r_ee(iw)
        if(ipr.eq.1) write(6,'(''cutjas_ee,c1_jas6_ee,c2_jas6_ee='',f10.6,9d12.5)')
     &       cutjas_ee,c1_jas6_ee(iw),c2_jas6_ee(iw)

       elseif(isc.eq.8 .or. isc.eq.10) then

        cutjas_en=1.d99
        val_cutjas_en=0
        cutjasi_en=0
        c1_jas6i_en=1
        c1_jas6_en(iw)=1
        c2_jas6_en(iw)=0
        asymp_r_en(iw)=c1_jas6i_en
        dasymp_r_en(iw)=0
        d2asymp_r_en(iw)=0

        cutjas_ee=1.d99
        val_cutjas_ee=0
        cutjasi_ee=0
        c1_jas6i_ee=1
        c1_jas6_ee(iw)=1
        c2_jas6_ee(iw)=0
        asymp_r_ee(iw)=asymp_r_en(iw)
        dasymp_r_ee(iw)=dasymp_r_en(iw)
        d2asymp_r_ee(iw)=d2asymp_r_en(iw)

       else

        cutjas_en=1.d99
        val_cutjas_en=0
        cutjasi_en=0
        c1_jas6i_en=1
        c1_jas6_en(iw)=1
        c2_jas6_en(iw)=0
        asymp_r_en(iw)=c1_jas6i_en/scalek(iw)
        dasymp_r_en(iw)=-1/(scalek(iw)*scalek(iw))
        d2asymp_r_en(iw)=2/(scalek(iw)*scalek(iw)*scalek(iw))

        cutjas_ee=1.d99
        val_cutjas_ee=0
        cutjasi_ee=0
        c1_jas6i_ee=1
        c1_jas6_ee(iw)=1
        c2_jas6_ee(iw)=0
        asymp_r_ee(iw)=asymp_r_en(iw)
        dasymp_r_ee(iw)=dasymp_r_en(iw)
        d2asymp_r_ee(iw)=d2asymp_r_en(iw)

      endif

c Calculate asymptotic value of A and B terms
      do 10 it=1,nctype
        asymp_jasa(it,iw)=a4(1,it,iw)*asymp_r_en(iw)/(1+a4(2,it,iw)*asymp_r_en(iw))
        dasymp_jasa(it,iw)=a4(1,it,iw)/(1+a4(2,it,iw)*asymp_r_en(iw))**2
        d2asymp_jasa(it,iw)=-2*dasymp_jasa(it,iw)*a4(2,it,iw)/(1+a4(2,it,iw)*asymp_r_en(iw))
        if(isc.eq.8 .or. isc.eq.10) then
          asymp_jasa(it,iw)=asymp_jasa(it,iw)/scalek(iw)
          dasymp_jasa(it,iw)=dasymp_jasa(it,iw)/scalek(iw)
          d2asymp_jasa(it,iw)=d2asymp_jasa(it,iw)/scalek(iw)
        endif
        do 10 iord=2,norda
          asymp_jasa(it,iw)=asymp_jasa(it,iw)+a4(iord+1,it,iw)*asymp_r_en(iw)**iord
          dasymp_jasa(it,iw)=dasymp_jasa(it,iw)+iord*a4(iord+1,it,iw)*asymp_r_en(iw)**(iord-1)
   10     d2asymp_jasa(it,iw)=d2asymp_jasa(it,iw)+iord*(iord-1)*a4(iord+1,it,iw)*asymp_r_en(iw)**(iord-2)

      if(ijas.eq.4) then
        do 20 i=1,2
          if(i.eq.1) then
            sspinn=1
            isp=1
           else
            if(nspin2b.eq.1.and.nocuspb.eq.0) then
              if(ndim.eq.3) then
                sspinn=half
               elseif(ndim.eq.2) then
                sspinn=third
              endif
             else
              sspinn=1
            endif
            isp=nspin2b
          endif
          asymp_jasb(i,iw)=sspinn*b(1,isp,iw)*asymp_r_ee(iw)/(1+b(2,isp,iw)*asymp_r_ee(iw))
          dasymp_jasb(i,iw)=sspinn*b(1,isp,iw)/(1+b(2,isp,iw)*asymp_r_ee(iw))**2
          d2asymp_jasb(i,iw)=-2*dasymp_jasb(i,iw)*b(2,isp,iw)/(1+b(2,isp,iw)*asymp_r_ee(iw))
          if(isc.eq.8 .or. isc.eq.10) then
            asymp_jasb(i,iw)=asymp_jasb(i,iw)/scalek(iw)
            dasymp_jasb(i,iw)=dasymp_jasb(i,iw)/scalek(iw)
            d2asymp_jasb(i,iw)=d2asymp_jasb(i,iw)/scalek(iw)
          endif
          do 20 iord=2,nordb
            asymp_jasb(i,iw)=asymp_jasb(i,iw)+b(iord+1,isp,iw)*asymp_r_ee(iw)**iord
            dasymp_jasb(i,iw)=dasymp_jasb(i,iw)+iord*b(iord+1,isp,iw)*asymp_r_ee(iw)**(iord-1)
   20       d2asymp_jasb(i,iw)=d2asymp_jasb(i,iw)+iord*(iord-1)*b(iord+1,isp,iw)*asymp_r_ee(iw)**(iord-2)
       elseif(ijas.eq.5) then
        do 35 i=1,2
          if(i.eq.1) then
            sspinn=1
            isp=1
           else
            if(nspin2b.eq.1.and.nocuspb.eq.0) then
              if(ndim.eq.3) then
                sspinn=half
               elseif(ndim.eq.2) then
                sspinn=third
              endif
             else
              sspinn=1
            endif
            isp=nspin2b
          endif
          asymp_jasb(i,iw)=b(1,isp,iw)*asymp_r_ee(iw)/(1+b(2,isp,iw)*asymp_r_ee(iw))
          do 30 iord=2,nordb
   30       asymp_jasb(i,iw)=asymp_jasb(i,iw)+b(iord+1,isp,iw)*asymp_r_ee(iw)**iord
   35     asymp_jasb(i,iw)=sspinn*asymp_jasb(i,iw)
      endif
      if((ijas.eq.4.or.ijas.eq.5).and.ipr.ge.-1) then
        write(6,'(''ifr,asymp_r_en,asymp_r_ee='',i2,2f10.6)') iw,asymp_r_en(iw),asymp_r_ee(iw)
        write(6,'(''ifr,asympa='',i2,10f10.6)') iw,(asymp_jasa(it,iw),it=1,nctype)
        write(6,'(''ifr,asympb='',i2,10f10.6)') iw,(asymp_jasb(i,iw),i=1,2)
        write(6,'(''ifr,dasymp_r_en,dasymp_r_ee='',i2,2f14.6)') iw,dasymp_r_en(iw),dasymp_r_ee(iw)
        write(6,'(''ifr,dasympa='',i2,10f10.6)') iw,(dasymp_jasa(it,iw),it=1,nctype)
        write(6,'(''ifr,dasympb='',i2,10f10.6)') iw,(dasymp_jasb(i,iw),i=1,2)
        write(6,'(''ifr,d2asymp_r_en,d2asymp_r_ee='',i2,2f14.6)') iw,d2asymp_r_en(iw),d2asymp_r_ee(iw)
        write(6,'(''ifr,d2asympa='',i2,10f10.6)') iw,(d2asymp_jasa(it,iw),it=1,nctype)
        write(6,'(''ifr,d2asympb='',i2,10f10.6)') iw,(d2asymp_jasb(i,iw),i=1,2)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_dist(r,rr,iflag)
c Written by Cyrus Umrigar
c Scale interparticle distances.

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      parameter (zero=0.d0,one=1.d0,half=0.5d0,third=1.d0/3.d0,d4b3=4.d0/3.d0)

      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 have zero 2nd and 3rd derivatives at 0
c isc = 6,7 are exponential and power law scalings resp. that have
c       zero, value and 1st 2 derivatives at cut_jas.
c isc = 8,10 are like 2 and 4 but multiplied by scalek.
c isc = 12,14, are  are similar 2,4 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c isc = 16,17, iflag=1,2 are the same as 6,7 resp.
c       This is used to scale r_en in J_en and r_ee in J_ee for solids.
c isc = 16,17, iflag=3 are similar 6,7 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_en in J_een for solids.
c isc = 16,17, iflag=4 have infinite range and are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_ee in J_een for solids.
c isc = 8,10 are the same as 2,4 resp. but multiplied by scalek so that at very
c       small scalek the normalizations of the basis functions while optimizing
c       is not many orders of magnitude different.  No longer needed because I
c       now use the diagonal of the overlap to normalize the basis functions.

      if(iflag.eq.1.or.iflag.eq.3) then
        cutjas=cutjas_en
        cutjasi=cutjasi_en
        c1_jas6=c1_jas6_en(iwf)
        c2_jas6=c2_jas6_en(iwf)
        asymp_r=asymp_r_en(iwf)
       elseif(iflag.eq.2.or.iflag.eq.4) then
        cutjas=cutjas_ee
        cutjasi=cutjasi_ee
        c1_jas6=c1_jas6_ee(iwf)
        c2_jas6=c2_jas6_ee(iwf)
        asymp_r=asymp_r_ee(iwf)
      endif

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2 .or. (isc.eq.12.and.(iflag.eq.1.or.iflag.eq.2))) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
	 elseif(isc.eq.3) then
	  rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
	  exprij=dexp(-rsc-half*rsc2-third*rsc3)
	  rr=(one-exprij)/scalek(iwf)
	 elseif(isc.eq.4 .or. (isc.eq.14.and.(iflag.eq.1.or.iflag.eq.2))) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=r*deni
 	 elseif(isc.eq.5) then
	  deni=one/(one+(scalek(iwf)*r)**3)**third
	  rr=r*deni
         elseif(isc.eq.6 .or. (isc.eq.16.and.(iflag.eq.1.or.iflag.eq.2))) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r*term))/scalek(iwf)
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r*term)-c2_jas6
            endif
          endif
         elseif(isc.eq.7 .or. (isc.eq.17.and.(iflag.eq.1.or.iflag.eq.2))) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            deni=1/(1+scalek(iwf)*r*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r*term*deni
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
	    endif
	  endif
         elseif(isc.eq.8) then
           rr=(one-dexp(-scalek(iwf)*r))
         elseif(isc.eq.10) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=scalek(iwf)*r*deni
         elseif(isc.eq.16 .and. iflag.ge.3) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r**2*term))/scalek(iwf)
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r**2*term)-c2_jas6
            endif
          endif
         elseif(isc.eq.12 .and. iflag.ge.3) then
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=(1-dexp(-scalek(iwf)*r**2))/scalek(iwf)
           elseif(ijas.eq.6) then
            rr=c1_jas6*dexp(-scalek(iwf)*r**2)-c2_jas6
          endif
         elseif(isc.eq.17 .and. iflag.ge.3) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            deni=1/(1+scalek(iwf)*r**2*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r**2*term*deni
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
	    endif
	  endif
         elseif(isc.eq.14 .and. iflag.ge.3) then
          deni=1/(1+scalek(iwf)*r**2)
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=r**2*deni
           elseif(ijas.eq.6) then
            rr=c1_jas6*deni-c2_jas6
	  endif
	endif
       else
        rr=r
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_dist1(r,rr,dd1,iflag)
c Written by Cyrus Umrigar
c Scale interparticle distances and calculate the 1st derivative
c of the scaled distances wrt the unscaled ones for calculating the
c gradient and laplacian.

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      parameter (zero=0.d0,one=1.d0,half=0.5d0,third=1.d0/3.d0,d4b3=4.d0/3.d0)

      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 have zero 2nd and 3rd derivatives at 0
c isc = 6,7 are exponential and power law scalings resp. that have
c       zero, value and 1st 2 derivatives at cut_jas.
c isc = 12,14, are  are similar 2,4 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c isc = 16,17, iflag=1,2 are the same as 6,7 resp.
c       This is used to scale r_en in J_en and r_ee in J_ee for solids.
c isc = 16,17, iflag=3 are similar 6,7 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_en in J_een for solids.
c isc = 16,17, iflag=4 have infinite range and are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_ee in J_een for solids.

      if(iflag.eq.1.or.iflag.eq.3) then
        cutjas=cutjas_en
        cutjasi=cutjasi_en
        c1_jas6=c1_jas6_en(iwf)
        c2_jas6=c2_jas6_en(iwf)
        asymp_r=asymp_r_en(iwf)
       elseif(iflag.eq.2.or.iflag.eq.4) then
        cutjas=cutjas_ee
        cutjasi=cutjasi_ee
        c1_jas6=c1_jas6_ee(iwf)
        c2_jas6=c2_jas6_ee(iwf)
        asymp_r=asymp_r_ee(iwf)
      endif

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2 .or. (isc.eq.12.and.(iflag.eq.1.or.iflag.eq.2))) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
          dd1=one-scalek(iwf)*rr
	 elseif(isc.eq.3) then
	  rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
	  exprij=dexp(-rsc-half*rsc2-third*rsc3)
	  rr=(one-exprij)/scalek(iwf)
	  dd1=(one+rsc+rsc2)*exprij
	 elseif(isc.eq.4 .or. (isc.eq.14.and.(iflag.eq.1.or.iflag.eq.2))) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=r*deni
	  dd1=deni*deni
 	 elseif(isc.eq.5) then
	  deni=one/(one+(scalek(iwf)*r)**3)**third
	  rr=r*deni
	  dd1=deni**4
         elseif(isc.eq.6 .or. (isc.eq.16.and.(iflag.eq.1.or.iflag.eq.2))) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r*term))/scalek(iwf)
c             rr_plus_c2=rr+c2_jas6
              dd1=(one-scalek(iwf)*rr)*term2
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r*term)-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*term2
            endif
          endif
         elseif(isc.eq.7 .or. (isc.eq.17.and.(iflag.eq.1.or.iflag.eq.2))) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            deni=1/(1+scalek(iwf)*r*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r*term*deni
              dd1=term2*deni*deni
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*deni*term2
	    endif
	  endif
         elseif(isc.eq.8) then
           rr=(one-dexp(-scalek(iwf)*r))
           dd1=(one-rr)*scalek(iwf)
         elseif(isc.eq.10) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=scalek(iwf)*r*deni
          dd1=deni*deni*scalek(iwf)
         elseif(isc.eq.16 .and. iflag.ge.3) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r**2*term))/scalek(iwf)
              dd1=2*r*(one-scalek(iwf)*rr)*term2
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r**2*term)-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*term2
            endif
          endif
         elseif(isc.eq.12 .and. iflag.ge.3) then
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=(1-dexp(-scalek(iwf)*r**2))/scalek(iwf)
            dd1=2*r*(one-scalek(iwf)*rr)
           elseif(ijas.eq.6) then
            rr=c1_jas6*dexp(-scalek(iwf)*r**2)-c2_jas6
c           rr_plus_c2=rr+c2_jas6
c           dd1=-scalek(iwf)*rr_plus_c2*term2
          endif
         elseif(isc.eq.17 .and. iflag.ge.3) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            deni=1/(1+scalek(iwf)*r**2*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r**2*term*deni
              dd1=2*r*term2*deni*deni
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*deni*term2
	    endif
	  endif
         elseif(isc.eq.14 .and. iflag.ge.3) then
          deni=1/(1+scalek(iwf)*r**2)
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=r**2*deni
            dd1=2*r*deni*deni
           elseif(ijas.eq.6) then
            rr=c1_jas6*deni-c2_jas6
c           rr_plus_c2=rr+c2_jas6
c           dd1=-scalek(iwf)*rr_plus_c2*deni*term2
	  endif
	endif
       else
        rr=r
        dd1=one
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_dist2(r,rr,dd1,dd2,iflag)
c Written by Cyrus Umrigar
c Scale interparticle distances and calculate the 1st and 2nd derivs
c of the scaled distances wrt the unscaled ones for calculating the
c gradient and laplacian.

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      parameter (zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,third=1.d0/3.d0,d4b3=4.d0/3.d0)

      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 have zero 2nd and 3rd derivatives at 0
c isc = 6,7 are exponential and power law scalings resp. that have
c       zero, value and 1st 2 derivatives at cut_jas.
c isc = 12,14, are  are similar 2,4 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c isc = 16,17, iflag=1,2 are the same as 6,7 resp.
c       This is used to scale r_en in J_en and r_ee in J_ee for solids.
c isc = 16,17, iflag=3 are similar 6,7 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_en in J_een for solids.
c isc = 16,17, iflag=4 have infinite range and are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_ee in J_een for solids.

      if(iflag.eq.1.or.iflag.eq.3) then
        cutjas=cutjas_en
        cutjasi=cutjasi_en
        c1_jas6=c1_jas6_en(iwf)
        c2_jas6=c2_jas6_en(iwf)
        asymp_r=asymp_r_en(iwf)
       elseif(iflag.eq.2.or.iflag.eq.4) then
        cutjas=cutjas_ee
        cutjasi=cutjasi_ee
        c1_jas6=c1_jas6_ee(iwf)
        c2_jas6=c2_jas6_ee(iwf)
        asymp_r=asymp_r_ee(iwf)
      endif

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2 .or. (isc.eq.12.and.(iflag.eq.1.or.iflag.eq.2))) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
          dd1=one-scalek(iwf)*rr
          dd2=-scalek(iwf)*dd1
	 elseif(isc.eq.3) then
	  rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
	  exprij=dexp(-rsc-half*rsc2-third*rsc3)
	  rr=(one-exprij)/scalek(iwf)
	  dd1=(one+rsc+rsc2)*exprij
	  dd2=-scalek(iwf)*rsc2*(3+2*rsc+rsc2)*exprij
	 elseif(isc.eq.4 .or. (isc.eq.14.and.(iflag.eq.1.or.iflag.eq.2))) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=r*deni
	  dd1=deni*deni
	  dd2=-two*scalek(iwf)*deni*dd1
 	 elseif(isc.eq.5) then
	  deni=one/(one+(scalek(iwf)*r)**3)**third
	  rr=r*deni
	  dd1=deni**4
	  dd2=-4*(scalek(iwf)*r)**2*scalek(iwf)*dd1*dd1/deni
         elseif(isc.eq.6 .or. (isc.eq.16.and.(iflag.eq.1.or.iflag.eq.2))) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
            dd2=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r*term))/scalek(iwf)
c             rr_plus_c2=rr+c2_jas6
              dd1=(one-scalek(iwf)*rr)*term2
              dd2=-scalek(iwf)*dd1*term2+2*(one-scalek(iwf)*rr)*cutjasi*(r_by_cut-1)
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r*term)-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*term2
c             dd2=-scalek(iwf)*(dd1*term2+rr_plus_c2*2*cutjasi*(r_by_cut-1))
            endif
          endif
         elseif(isc.eq.7 .or. (isc.eq.17.and.(iflag.eq.1.or.iflag.eq.2))) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
            dd2=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            deni=1/(1+scalek(iwf)*r*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r*term*deni
              dd1=term2*deni*deni
              dd2=-2*deni**2*(scalek(iwf)*deni*term2**2-cutjasi*(r_by_cut-1))
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*deni*term2
c             dd2=-2*scalek(iwf)*deni*(dd1*term2+rr_plus_c2*cutjasi*(r_by_cut-1))
	    endif
	  endif
         elseif(isc.eq.8) then
           rr=(one-dexp(-scalek(iwf)*r))
           dd1=(one-rr)*scalek(iwf)
           dd2=-scalek(iwf)*dd1
         elseif(isc.eq.10) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=scalek(iwf)*r*deni
          dd1=deni*deni*scalek(iwf)
          dd2=-two*scalek(iwf)*dd1*deni
         elseif(isc.eq.16 .and. iflag.ge.3) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
            dd2=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r**2*term))/scalek(iwf)
              dd1=2*r*(one-scalek(iwf)*rr)*term2
c             dd2=-scalek(iwf)*dd1*term2+2*(one-scalek(iwf)*rr)*cutjasi*(r_by_cut-1)
              dd2=-2*(scalek(iwf)*r*dd1*term2
     &        -(one-scalek(iwf)*rr)*(1-4*r_by_cut+3*r_by_cut2))
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r**2*term)-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*term2
c             dd2=-scalek(iwf)*(dd1*term2+rr_plus_c2*2*cutjasi*(r_by_cut-1))
            endif
          endif
         elseif(isc.eq.12 .and. iflag.ge.3) then
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=(1-dexp(-scalek(iwf)*r**2))/scalek(iwf)
            dd1=2*r*(one-scalek(iwf)*rr)
            dd2=-2*(scalek(iwf)*r*dd1-(one-scalek(iwf)*rr))
           elseif(ijas.eq.6) then
            rr=c1_jas6*dexp(-scalek(iwf)*r**2)-c2_jas6
c           rr_plus_c2=rr+c2_jas6
c           dd1=-scalek(iwf)*rr_plus_c2*term2
c           dd2=-scalek(iwf)*(dd1*term2+rr_plus_c2*2*cutjasi*(r_by_cut-1))
          endif
         elseif(isc.eq.17 .and. iflag.ge.3) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
            dd2=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            deni=1/(1+scalek(iwf)*r**2*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r**2*term*deni
              dd1=2*r*term2*deni*deni
c             dd2=-2*deni**2*(scalek(iwf)*deni*term2**2-cutjasi*(r_by_cut-1))
              dd2=-2*deni**2*(4*scalek(iwf)*r**2*deni*term2**2
     &        -(1-4*r_by_cut+3*r_by_cut2))
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*deni*term2
c             dd2=-2*scalek(iwf)*deni*(dd1*term2+rr_plus_c2*cutjasi*(r_by_cut-1))
	    endif
	  endif
         elseif(isc.eq.14 .and. iflag.ge.3) then
          deni=1/(1+scalek(iwf)*r**2)
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=r**2*deni
            dd1=2*r*deni*deni
            dd2=-2*deni**2*(4*scalek(iwf)*r**2*deni-1)
           elseif(ijas.eq.6) then
            rr=c1_jas6*deni-c2_jas6
c           rr_plus_c2=rr+c2_jas6
c           dd1=-scalek(iwf)*rr_plus_c2*deni*term2
c           dd2=-2*scalek(iwf)*deni*(dd1*term2+rr_plus_c2*cutjasi*(r_by_cut-1))
	  endif
	endif
       else
        rr=r
        dd1=one
        dd2=0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine switch_scale(rr,iflag)
c Written by Cyrus Umrigar
c Switch scaling for ijas=4,5 from that appropriate for A,B terms to
c that appropriate for C terms, for dist.

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      if(iflag.eq.1.or.iflag.eq.3) then
        c1_jas6=c1_jas6_en(iwf)
       elseif(iflag.eq.2.or.iflag.eq.4) then
        c1_jas6=c1_jas6_ee(iwf)
      endif

      if(isc.eq.8 .or. isc.eq.10) then
        rr=1-c1_jas6*rr
      else
        rr=1-c1_jas6*scalek(iwf)*rr
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine switch_scale1(rr,dd1,iflag)
c Written by Cyrus Umrigar
c Switch scaling for ijas=4,5 from that appropriate for A,B terms to
c that appropriate for C terms, for dist and 1st deriv.

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      if(iflag.eq.1.or.iflag.eq.3) then
        c1_jas6=c1_jas6_en(iwf)
       elseif(iflag.eq.2.or.iflag.eq.4) then
        c1_jas6=c1_jas6_ee(iwf)
      endif

      if(isc.eq.8 .or. isc.eq.10) then
        rr=1-c1_jas6*rr
        dd1=-c1_jas6*dd1
      else
        rr=1-c1_jas6*scalek(iwf)*rr
        dd1=-c1_jas6*scalek(iwf)*dd1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine switch_scale2(rr,dd1,dd2,iflag)
c Written by Cyrus Umrigar
c Switch scaling for ijas=4,5 from that appropriate for A,B terms to
c that appropriate for C terms, for dist and 1st two derivs.

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      if(iflag.eq.1.or.iflag.eq.3) then
        c1_jas6=c1_jas6_en(iwf)
       elseif(iflag.eq.2.or.iflag.eq.4) then
        c1_jas6=c1_jas6_ee(iwf)
      endif

      if(isc.eq.8 .or. isc.eq.10) then
        rr=1-c1_jas6*rr
        dd1=-c1_jas6*dd1
        dd2=-c1_jas6*dd2
      else
        rr=1-c1_jas6*scalek(iwf)*rr
        dd1=-c1_jas6*scalek(iwf)*dd1
        dd2=-c1_jas6*scalek(iwf)*dd2
      endif

      return
      end
c-----------------------------------------------------------------------

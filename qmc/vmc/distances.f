      subroutine distances(x,pe,pei)
c Written by Cyrus Umrigar
c calculate interparticle distances
      use control_mod
      use deriv_orb_mod
      use eloc_mod
      use atom_mod
      use const_mod
      use dim_mod
      use pseudo_mod
      use contrl_per_mod
      use distance_mod
      use periodic_1d_mod
      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /dotcenter/ dot_bump_height, dot_bump_radius, dot_bump_radius_inv2
      common /wire/ wire_w,wire_length,wire_length2,wire_radius2, wire_potential_cutoff,wire_prefactor,wire_root1
      common /angularpert/ ang_perturb,amp_perturb,shrp_perturb,iperturb
c     common /compferm/ emagv,nv,idot
      common /jel_sph1/ dn_background,rs_jel,radius_b !RM

c Warning: temporary

      dimension x(3,*)

c  pe from nucleus-nucleus repulsion
      pe=pecent
      if(iperiodic.eq.0) then

c Calculate e-N inter-particle distances
        pe_en=0.d0 !JT
        do 26 ic=1,ncent
          do 26 i=1,nelec
            r_en(i,ic)=0
            do 25 k=1,ndim
              rvec_en(k,i,ic)=x(k,i)-cent(k,ic)
   25         r_en(i,ic)=r_en(i,ic)+rvec_en(k,i,ic)**2
            r_en(i,ic)=dsqrt(r_en(i,ic))
            if(nloc.eq.0) pe_en=pe_en-znuc(iwctype(ic))/r_en(i,ic)
c           if(nloc.eq.-1) pe_en=pe_en+0.5d0*(w0*r_en(i,ic))**2
            if((nloc.eq.-1 .or. nloc.eq.-5) .and. ic.eq.1) then  
c             pe_en=pe_en+0.5d0*(w0*r_en(i,ic))**2
c             emag=emag+0.125d0*(bext*r_en(i,ic))**2
              pe_en=pe_en+0.5d0*(we*(r_en(i,ic)-rring))**2
              if(nloc.eq.-5 .and. r_en(i,ic).lt.dot_bump_radius) then
                 pe_en = pe_en + dot_bump_height*dexp(1.0d0 - 1.0d0/(1.0d0 - dot_bump_radius_inv2*(r_en(i,ic)**2)))
              endif
              if(iperturb.eq.1) then
c note that the perturbation potential is defined only for -pi<theta<pi  (it's not periodic)
                theta=datan2(x(2,i),x(1,i))
                pe_en=pe_en+amp_perturb*(tanh(shrp_perturb*(theta+ang_perturb))
     &                            -tanh(shrp_perturb*(theta-ang_perturb)))
               elseif(iperturb.eq.2) then               ! this is for testing
                theta=datan2(x(2,i),x(1,i))
                pe_en=pe_en+2.d0*amp_perturb*dexp(-0.5d0*(theta/(ang_perturb))**2)
              endif
            endif
            if(nloc.eq.-2) pe_en=pe_en+p1*rvec_en(1,i,ic)**4+p2*rvec_en(2,i,ic)**4
     &      -2*p3*(rvec_en(1,i,ic)*rvec_en(2,i,ic))**2
     &      +p4*(rvec_en(1,i,ic)-rvec_en(2,i,ic))*rvec_en(1,i,ic)*rvec_en(2,i,ic)*r_en(i,ic)
c Contribution from Jellium to the potential energy. A temporary patch which should be more smart? ! RM
            if(nloc.eq.-3) then ! jellium RM
              if(r_en(i,ic).ge.radius_b)then
                p_bg=-dn_background/r_en(i,ic)
               else
                p_bg=-0.5d0*(dn_background/radius_b)*(3.d0-(r_en(i,ic)/radius_b)**2)
              endif
              pe=pe+p_bg
! MS The potential energy between Z and e
              pe_en=pe_en-znuc(iwctype(ic))/r_en(i,ic)
            endif
            if(nloc.eq.-4) then  !  quantum wire
              pe_y=0.5d0*(wire_w*x(2,i))**2  !  y-direction
              pe_en=pe_en+pe_y
              if(iperturb.eq.1) then
c note that in wires, ang_perturb is the semi-width of the perturbation in the x-direction
                pe_en=pe_en+amp_perturb*(tanh(shrp_perturb*(x(1,i)+ang_perturb))
     &                            -tanh(shrp_perturb*(x(1,i)-ang_perturb)))
              endif
c      These are used to calculate confining, x-direction potential:
              xshift=x(1,i)+wire_length*0.5d0
              wire_root2 = dsqrt(wire_radius2 + xshift**2)
              wire_root3 = dsqrt(wire_radius2 + (wire_length - xshift)**2)
c      The x-direction (confining) potential:
              pe_x=wire_prefactor * 
     &((1.d0 - 2*wire_potential_cutoff**2) * wire_length2 
     & + wire_length * (2*wire_potential_cutoff*wire_root1-wire_root3-2*xshift)
     &+ xshift * (wire_root3 + 2 * xshift - wire_root2) 
     &+ wire_radius2 * dlog( ((wire_potential_cutoff * wire_length 
     &+ wire_root1)**2) 
     &/ ((wire_length + wire_root3 - xshift)*(xshift + wire_root2))  ) )
              pe_en=pe_en+pe_x
c              write(6,*) 'wire_root1, wire_root2, wire_root3 = ', wire_root1, wire_root2, wire_root3
c              write(6,*) 'wire_radius2, wire_length, wire_prefactor ', wire_radius2, wire_length, wire_prefactor
c              write(6,*) 'i,x(1,i),x(2,i),pe_x,pe_y=',i,x(1,i),x(2,i),pe_x,pe_y
            endif
   26   continue

c Calculate e-e inter-particle distances
        pe_ee=0.d0
        ij=0
        do 29 i=2,nelec
          do 29 j=1,i-1
            ij=ij+1
            r_ee(ij)=0
            do 28 k=1,ndim
              rvec_ee(k,ij)=x(k,i)-x(k,j)
   28         r_ee(ij)=r_ee(ij)+rvec_ee(k,ij)**2
            r_ee(ij)=dsqrt(r_ee(ij))
            pe_ee=pe_ee+1/r_ee(ij)
   29   continue

        pe=pe+pe_en+pe_ee !JT
        pei=pe_ee

       elseif(iperiodic.eq.1)then ! periodic in 1d
         pe_en=0.d0
         do ic=1,ncent
          do i=1,nelec
            do k=1,ndim
              rvec_en(k,i,ic)=x(k,i)-cent(k,ic)
            enddo
            call find_image_1d(rvec_en(:,i,ic), r_en(i,ic)) ! modulo math
            if(nloc.eq.-4) then  !  quantum wire
              pe_y=0.5d0*(wire_w*x(2,i))**2  !  y-direction
              pe_en=pe_en+pe_y
              if(iperturb.eq.1) then
c note that in wires, ang_perturb is the semi-width of the perturbation in the x-direction
                pe_en=pe_en+amp_perturb*(tanh(shrp_perturb*(x(1,i)+ang_perturb))
     &                            -tanh(shrp_perturb*(x(1,i)-ang_perturb)))
              endif
            endif
          enddo
         enddo
         call pot_ee_ewald_1d(x,pe_ee)
         pe = pe+pe_en+pe_ee
         pei = pe_ee

       else

        call pot_en_ewald(x,pe_en)
        call pot_ee_ewald(x,pe_ee)
        pe=pe+pe_en+pe_ee
        pei=pe_ee

c     write(6,*) 'in distances'
c     write(6,'(''r_en(i,j)'',9f9.5)') ((r_en(i,j),i=1,nelec),j=1,2)
c     write(6,'(''r_ee(ij)'',9f9.5)') (r_ee(ij),ij=1,nelec*(nelec-1)/2)
c     write(6,'(''rvec_ee(k,ij)'',9f12.4)') ((rvec_ee(k,ij),k=1,ndim),ij=1,nelec*(nelec-1)/2)
      if(ipr.ge.3) write(6,'(''pe,pe_en(loc),pe_ee'',11f9.5)') pe,pe_en,pe_ee

      endif

      eloc_pot_loc = pe
      call object_modified_by_index (eloc_pot_loc_index)  !JT
      call object_modified_by_index (pe_ee_index)  !JT
      call object_modified_by_index (pe_en_index)  !JT
      call object_modified_by_index (r_en_index) !JT
      call object_modified_by_index (rvec_en_index) !JT

      return
      end
c-----------------------------------------------------------------------
      subroutine distancese(iel,x)
c Written by Cyrus Umrigar
c calculate distances of electron iel to all other particles
      use control_mod
      use atom_mod
      use const_mod
      use dim_mod
      use contrl_per_mod
      use distance_mod
      use distance_sav_mod
      use periodic_mod
      use periodic_1d_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*)

c Calculate e-N inter-particle distances
      do 30 ic=1,ncent
        r_en_sav(ic)=r_en(iel,ic)
        do 10 k=1,ndim
          rshift_sav(k,ic)=rshift(k,iel,ic)
          rvec_en_sav(k,ic)=rvec_en(k,iel,ic)
   10     rvec_en(k,iel,ic)=x(k,iel)-cent(k,ic)
        if(iperiodic.eq.0) then
          r_en(iel,ic)=0
          do 20 k=1,ndim
   20       r_en(iel,ic)=r_en(iel,ic)+rvec_en(k,iel,ic)**2
          r_en(iel,ic)=dsqrt(r_en(iel,ic))
        elseif(iperiodic.eq.1)then
          call find_image_1d(rvec_en(:,iel,ic),r_en(iel,ic)) ! modulo math
        else
          call find_image4(rshift(1,iel,ic),rvec_en(1,iel,ic),r_en(iel,ic),rlatt,rlatt_inv)
        endif
   30 continue

c Calculate e-e inter-particle distances
      do 60 jj=1,nelec

        if(jj.eq.iel) goto 60
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif
        ij=((i-1)*(i-2))/2+j

        r_ee_sav(jj)=r_ee(ij)
        do 40 k=1,ndim
          rvec_ee_sav(k,jj)=rvec_ee(k,ij)
   40     rvec_ee(k,ij)=x(k,i)-x(k,j)
        if(iperiodic.eq.0) then
          r_ee(ij)=0
          do 50 k=1,ndim
   50       r_ee(ij)=r_ee(ij)+rvec_ee(k,ij)**2
          r_ee(ij)=dsqrt(r_ee(ij))
         elseif(iperiodic.eq.1) then
            call find_image_1d(rvec_ee(:,ij), r_ee(ij))
         else
          call find_image3(rvec_ee(1,ij),r_ee(ij),rlatt_sim,rlatt_sim_inv)
        endif
   60 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine distancese_restore(iel)
c Written by Cyrus Umrigar
c restore interparticle distances (called if move rejected)
      use atom_mod
      use const_mod
      use dim_mod
      use distance_mod
      use distance_sav_mod
      implicit real*8(a-h,o-z)


c Restore e-N inter-particle distances
      do 25 ic=1,ncent
        r_en(iel,ic)=r_en_sav(ic)
        do 25 k=1,ndim
          rshift(k,iel,ic)=rshift_sav(k,ic)
   25     rvec_en(k,iel,ic)=rvec_en_sav(k,ic)

c Restore e-e inter-particle distances
      do 29 jj=1,nelec

        if(jj.eq.iel) goto 29
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif
        ij=((i-1)*(i-2))/2+j

        r_ee(ij)=r_ee_sav(jj)
        do 28 k=1,ndim
   28     rvec_ee(k,ij)=rvec_ee_sav(k,jj)
   29 continue

      return
      end
c-----------------------------------------------------------------------

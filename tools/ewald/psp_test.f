c write out a model pseudopotential that is cubic or quartic at the
c origin and goes as -1/r at infinity

      implicit real*8(a-h,o-z)
  
      write(6,'(''Input: z,a (v(inf)=-z/r, v(0)=-z*a)'')')
      read(5,*) Z,a
      write(6,'(''Input: n,rmax,arg_ps'')')
      read(5,*) n,rmax,arg_ps
      write(6,'(''Input: n_deriv, (the 1st non-zero deriv)'')')
      read(5,*) n_deriv
      r0_ps=rmax/(arg_ps**(n-1)-1)

      if(n_deriv.eq.1) write(3,'(f4.1,'' 3.87'',f7.3,'' Z, rcut, rmax, (linear at r=0, v(r)=z*(dexp(-a*r)-1)/r'')') Z,rmax
      if(n_deriv.eq.2) write(3,'(f4.1,'' 3.87'',f7.3,'' Z, rcut, rmax,'',
     & '' (quadratic at r=0, v(r)=z*(dexp(-a*r*(1+a*r/2))-1)/r'')') Z,rmax
      if(n_deriv.eq.3) write(3,'(f4.1,'' 3.87'',f7.3,'' Z, rcut, rmax,'',
     & '' (cubic at r=0, v(r)=z*(dexp(-a*r*(1+(a*r/2)*(1+a*r*2/3)))-1)/r'')') Z,rmax
      if(n_deriv.eq.4) write(3,'(f4.1,'' 3.87'',f7.3, '' Z, rcut, rmax,'',
     &'' (quartic at r=0, v(r)=z*(dexp(-a*r*(1+(a*r/2)*(1+(2*(a*r)/3)*(1+3*a*r/4))))-1)/r'')') Z,rmax

      write(3,'(i5,f8.4,'' n,arg_ps'')') n,arg_ps
      do 10 i=1,n
        r=max(r0_ps*(arg_ps**(i-1)-1),1.d-9)
        if(r.gt.1.d-3) then
          if(n_deriv.eq.1) f=z*(dexp(-a*r)-1)/r
          if(n_deriv.eq.2) f=z*(dexp(-a*r*(1+a*r/2))-1)/r
          if(n_deriv.eq.3) f=z*(dexp(-a*r*(1+(a*r/2)*(1+a*r*2/3)))-1)/r
          if(n_deriv.eq.4) f=z*(dexp(-a*r*(1+(a*r/2)*(1+(2*(a*r)/3)*(1+3*a*r/4))))-1)/r
         else
          if(n_deriv.eq.1) f=-z*a*(1-(a*r)**1/2+(a*r)**2/6-(a*r)**3/24)
          if(n_deriv.eq.2) f=-z*a*(1-(a*r)**2/3+(a*r)**3/12+(a*r)**4/20)
          if(n_deriv.eq.3) f=-z*a*(1-(a*r)**3/4+(a*r)**4/20+(a*r)**5/30)
          if(n_deriv.eq.4) f=-z*a*(1-(a*r)**4/5+(a*r)**5/30+(a*r)**6/42)
        endif
   10   write(3,'(1p9d20.12)') r, f, z/r

      stop
      end

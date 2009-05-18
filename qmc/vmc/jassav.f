      subroutine jassav(iel)
c Written by Claudia Filippi

      use const_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'

!JT      include 'force.h'

!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /jasn/ fsn(MELEC,MELEC),fijn(3,MELEC,MELEC)
     &,d2ijn(MELEC,MELEC),d2n,fsumn,fjn(3,MELEC)
      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)

      fsumo=fsumn
      do 10 i=1,nelec
        fjo(1,i)=fjn(1,i)
        fjo(2,i)=fjn(2,i)
  10    fjo(3,i)=fjn(3,i)

      do 20 j=1,iel
  20    fso(iel,j)=fsn(iel,j)

      do 30 j=iel+1,nelec
  30    fso(j,iel)=fsn(j,iel)

      do 40 j=1,nelec
        fijo(1,iel,j)=fijn(1,iel,j)
        fijo(2,iel,j)=fijn(2,iel,j)
        fijo(3,iel,j)=fijn(3,iel,j)
        fijo(1,j,iel)=fijn(1,j,iel)
        fijo(2,j,iel)=fijn(2,j,iel)
  40    fijo(3,j,iel)=fijn(3,j,iel)

      return
      end

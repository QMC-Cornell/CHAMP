      subroutine walksav_jas(iw)
! Written by Claudia Filippi

      use const_mod
      use jaso_mod
      use jasow_mod
      implicit real*8(a-h,o-z)

      fsumow(iw)=fsumo

      do 10 i=1,nelec
        lapjow(i,iw)=lapjo(i)
        fjow(1,i,iw)=fjo(1,i)
        fjow(2,i,iw)=fjo(2,i)
  10    fjow(3,i,iw)=fjo(3,i)

      do 20 i=2,nelec
        do 20 j=1,i-1
        lapjijow(i,j,iw)=lapjijo(i,j)
        lapjijow(j,i,iw)=lapjijo(j,i)
        fsow(i,j,iw)=fso(i,j)
        fijow(1,i,j,iw)=fijo(1,i,j)
        fijow(2,i,j,iw)=fijo(2,i,j)
        fijow(3,i,j,iw)=fijo(3,i,j)
        fijow(1,j,i,iw)=fijo(1,j,i)
        fijow(2,j,i,iw)=fijo(2,j,i)
  20    fijow(3,j,i,iw)=fijo(3,j,i)

      do 25 i=1,nelec
        lapjijow(i,i,iw)=lapjijo(i,i)
        fsow(i,i,iw)=fso(i,i)
        fijow(1,i,i,iw)=fijo(1,i,i)
        fijow(2,i,i,iw)=fijo(2,i,i)
  25    fijow(3,i,i,iw)=fijo(3,i,i)

      return

      entry walkstrjas(iw)

      fsumo=fsumow(iw)

      do 30 i=1,nelec
        lapjo(i)=lapjow(i,iw)
        fjo(1,i)=fjow(1,i,iw)
        fjo(2,i)=fjow(2,i,iw)
  30    fjo(3,i)=fjow(3,i,iw)

      do 40 i=2,nelec
        do 40 j=1,i-1
        lapjijo(i,j)=lapjijow(i,j,iw)
        lapjijo(j,i)=lapjijow(j,i,iw)
        fso(i,j)=fsow(i,j,iw)
        fijo(1,i,j)=fijow(1,i,j,iw)
        fijo(2,i,j)=fijow(2,i,j,iw)
        fijo(3,i,j)=fijow(3,i,j,iw)
        fijo(1,j,i)=fijow(1,j,i,iw)
        fijo(2,j,i)=fijow(2,j,i,iw)
  40    fijo(3,j,i)=fijow(3,j,i,iw)

      do 45 i=1,nelec
        lapjijo(i,i)=lapjijow(i,i,iw)
        lapjijo(i,i)=lapjijow(i,i,iw)
        fso(i,i)=fsow(i,i,iw)
        fijo(1,i,i)=fijow(1,i,i,iw)
        fijo(2,i,i)=fijow(2,i,i,iw)
  45    fijo(3,i,i)=fijow(3,i,i,iw)

      return

      entry splitjjas(iw,iw2)

      fsumow(iw2)=fsumow(iw)

      do 50 i=1,nelec
        lapjow(i,iw2)=lapjow(i,iw)
        fjow(1,i,iw2)=fjow(1,i,iw)
        fjow(2,i,iw2)=fjow(2,i,iw)
  50    fjow(3,i,iw2)=fjow(3,i,iw)

      do 60 i=2,nelec
        do 60 j=1,i-1
        lapjijow(i,j,iw2)=lapjijow(i,j,iw)
        lapjijow(j,i,iw2)=lapjijow(j,i,iw)
        fsow(i,j,iw2)=fsow(i,j,iw)
        fijow(1,i,j,iw2)=fijow(1,i,j,iw)
        fijow(2,i,j,iw2)=fijow(2,i,j,iw)
        fijow(3,i,j,iw2)=fijow(3,i,j,iw)
        fijow(1,j,i,iw2)=fijow(1,j,i,iw)
        fijow(2,j,i,iw2)=fijow(2,j,i,iw)
  60    fijow(3,j,i,iw2)=fijow(3,j,i,iw)

      do 65 i=1,nelec
        lapjijow(i,i,iw2)=lapjijow(i,i,iw)
        lapjijow(i,i,iw2)=lapjijow(i,i,iw)
        fsow(i,i,iw2)=fsow(i,i,iw)
        fijow(1,i,i,iw2)=fijow(1,i,i,iw)
        fijow(2,i,i,iw2)=fijow(2,i,i,iw)
  65    fijow(3,i,i,iw2)=fijow(3,i,i,iw)

      return
      end

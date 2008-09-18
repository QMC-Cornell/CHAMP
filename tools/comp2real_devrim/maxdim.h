c MAXNELEC=maximum number of electrons
c MAXNDET=maximum number of determinants after expansion, but before reduction
c          ~< 2^MAXNELEC
c MAXQN=number of quantum states = 4 as usually (n,m,l,s)
c MAXS,MAXN,MAXL,MAXM = maximum for absolute value of quantum numbers.
      parameter        (MAXNELEC=15)
      parameter        (MAXNDET=2**MAXNELEC)
      parameter        (MAXQN=4)
      parameter        (MAXS=1,MAXN=11,MAXL=10,MAXM=10)

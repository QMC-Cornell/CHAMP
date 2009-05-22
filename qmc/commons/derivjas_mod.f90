module derivjas_mod

 use constants_mod
 implicit none
 save

 double precision gvalue(MPARMJ),g(3,MELEC,MPARMJ),d2g(MPARMJ)
 double precision go(MELEC,MELEC,MPARMJ)

end module derivjas_mod

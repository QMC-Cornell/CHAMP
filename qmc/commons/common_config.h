      double precision xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)               &
      ,vnew(3,MELEC),psi2o(MFORCE),psi2n(MFORCE),eold(MFORCE),enew(MFORCE)     &
      ,peo,pen,peio,pein,tjfn,tjfo,psido,psijo                                           &
      ,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)               &
      ,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)           &
      ,delttn(MELEC)
      integer nearesto(MELEC),nearestn(MELEC)

      common /config/ xold,xnew,vold,vnew,psi2o,psi2n,eold,enew    &
      ,peo,pen,peio,pein,tjfn,tjfo,psido,psijo,rmino,rminn,rvmino,rvminn     &
      ,rminon,rminno,rvminon,rvminno                               &
      ,nearesto,nearestn,delttn

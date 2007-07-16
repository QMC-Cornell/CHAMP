      double precision slmui(MMAT_DIM,MDETUD),slmdi(MMAT_DIM,MDETUD) 
      double precision fpu(3,MMAT_DIM,MDETUD),fpd(3,MMAT_DIM,MDETUD)
      double precision fppu(MMAT_DIM,MDETUD),fppd(MMAT_DIM,MDETUD)
      double precision detu(MDETUD),detd(MDETUD)
      double precision ddeti_deti(3,MELEC,MDETUD),d2edeti_deti(MELEC,MDETUD)
      double precision deti_det(MPARMD),ddeti_det(3,MELEC,MPARMD)
      double precision d2deti_det(MPARMD),d2det_det,detij_det(MPARMD,MPARMD)

      common /slater/ slmui,slmdi,fpu,fpd,fppu,fppd,detu,detd,ddeti_deti,d2edeti_deti,deti_det,ddeti_det,d2deti_det,d2det_det,detij_det

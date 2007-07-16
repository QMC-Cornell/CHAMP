      double precision wtgen(0:MFPRD1),ff(0:MFPRD1),eoldw(MWALK,MFORCE)
      double precision pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE)
      double precision wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod
      integer nwalk

      common /branch/ wtgen,ff,eoldw,pwt,wthist,wt,eigv,eest,wdsumo,wgdsumo,fprod,nwalk


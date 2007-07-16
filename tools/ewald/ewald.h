c NCOEFX     max number of coefficients in short-range polynomial for optimal Ewald 
c NPX        is the number of factors of (r/cutr-1) we multiply polynomial by (not used anymore)
c MKPTS      maximum number of k-pts.  This can be at most the ratio of the simulation to primitive cell,
c            volumes, V_sim/V.  However, since some k-pts have two indep. states (most if V_sim/V>>1)
c            the number of k-pts is between V_sim/(2V) and V_sim/V.  I have not tried to distinguish
c            between arrays that could be dimensioned to nkpts and those that need V_sim/V, but
c            have used the larger number for all.
c IVOL_RATIO is the ratio of the simulation to primitive cell volumes
c            However, since (cutg_sim/cutg) and (cutg_sim_big/cutg_big) can be chosen to
c            be (V_sim/V)^(-1/3), this ratio can just be set to 1 and I could remove this parameter.
c IBIG_RATIO is the ratio of the number of vectors or norms before the optimal Ewald separation
c            to the number after the separation. This ratio is (cutg_big/cutg)^3
c            and (cutg_sim_big/cutg_sim)^3 for the primitive and simulation cells resp.
c            This ratio is accurate for vectors and an overestimate for norms.
c NSYM       is the ratio of the number of vectors to the number of norms
c            and depends on the symmetry of the lattice.  Since many stars have less
c            vectors than the number of symmetry operations, this number can be set
c            somewhat smaller than the number of sym. ops.  e.g., for cubic one can
c            try 32 instead of 48.
c     parameter (NCOEFX=10, NPX=4, IVOL_RATIO=36, IBIG_RATIO=20, NSYM=48
c     parameter (NCOEFX=10, NPX=4, IVOL_RATIO=10, IBIG_RATIO=15, NSYM=8
c     parameter (NCOEFX=10, NPX=4, IVOL_RATIO=28, IBIG_RATIO=15, NSYM=8
c     parameter (NCOEFX=20, NPX=4, MKPTS=32, IVOL_RATIO=1, IBIG_RATIO=15, NSYM=32
c     parameter (NCOEFX=20, NPX=4, MKPTS=64, IVOL_RATIO=1, IBIG_RATIO=10, NSYM=32

      parameter (NCOEFX=20, NPX=4, MKPTS=32, IVOL_RATIO=1, IBIG_RATIO=15, NSYM=32
     &,NGNORMX=190, NGVECX=NGNORMX*NSYM, NGVEC2X=2*NGVECX, NG1DX=60
     &,NGNORM_SIMX=NGNORMX*IVOL_RATIO, NGVEC_SIMX=NGVECX*IVOL_RATIO
     &,NGNORM_BIGX=IBIG_RATIO*NGNORMX, NGVEC_BIGX=IBIG_RATIO*NGVECX
     &,NGNORM_SIM_BIGX=IBIG_RATIO*NGNORM_SIMX, NGVEC_SIM_BIGX=IBIG_RATIO*NGVEC_SIMX)

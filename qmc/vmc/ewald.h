! NCOEFX     max number of coefficients in short-range polynomial for optimal Ewald 
! NPX        is the number of factors of (r/cutr-1) we multiply polynomial by (not used anymore)
! MKPTS      maximum number of k-pts.  This can be at most the ratio of the simulation to primitive cell,
!            volumes, V_sim/V.  However, since some k-pts have two indep. states (most if V_sim/V>>1)
!            the number of k-pts is between V_sim/(2V) and V_sim/V.  I have not tried to distinguish
!            between arrays that could be dimensioned to nkpts and those that need V_sim/V, but
!            have used the larger number for all.
! IVOL_RATIO is the ratio of the simulation to primitive cell volumes
!            However, since (cutg_sim/cutg) and (cutg_sim_big/cutg_big) can be chosen to
!            be (V_sim/V)^(-1/3), this ratio can just be set to 1 and I could remove this parameter.
! IBIG_RATIO is the ratio of the number of vectors or norms before the optimal Ewald separation
!            to the number after the separation. This ratio is (cutg_big/cutg)^3
!            and (cutg_sim_big/cutg_sim)^3 for the primitive and simulation cells resp.
!            This ratio is accurate for vectors and an overestimate for norms.
! NSYM       is the ratio of the number of vectors to the number of norms
!            and depends on the symmetry of the lattice.  Since many stars have less
!            vectors than the number of symmetry operations, this number can be set
!            somewhat smaller than the number of sym. ops.  e.g., for cubic one can
!            try 32 instead of 48.
!     parameter (NCOEFX=10, NPX=4, IVOL_RATIO=36, IBIG_RATIO=20, NSYM=48
!     parameter (NCOEFX=10, NPX=4, IVOL_RATIO=10, IBIG_RATIO=15, NSYM=8
!     parameter (NCOEFX=10, NPX=4, IVOL_RATIO=28, IBIG_RATIO=15, NSYM=8
!     parameter (NCOEFX=20, NPX=4, MKPTS=32, IVOL_RATIO=1, IBIG_RATIO=15, NSYM=32
!     parameter (NCOEFX=20, NPX=4, MKPTS=64, IVOL_RATIO=1, IBIG_RATIO=10, NSYM=32
      parameter(NCOEFX=20, NPX=4, MKPTS=32, IVOL_RATIO=1, IBIG_RATIO=15, NSYM=32)
      parameter(NGNORMX=500, NGVECX=NGNORMX*NSYM, NGVEC2X=2*NGVECX, NG1DX=60)
      parameter(NGNORM_SIMX=NGNORMX*IVOL_RATIO, NGVEC_SIMX=NGVECX*IVOL_RATIO)
      parameter(NGNORM_BIGX=IBIG_RATIO*NGNORMX, NGVEC_BIGX=IBIG_RATIO*NGVECX)
      parameter(NGNORM_SIM_BIGX=IBIG_RATIO*NGNORM_SIMX, NGVEC_SIM_BIGX=IBIG_RATIO*NGVEC_SIMX)

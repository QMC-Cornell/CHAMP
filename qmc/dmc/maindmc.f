      subroutine maindmc
      use dmc_mod
      implicit real*8(a-h,o-z)

      call dmc_init
      call opt_wf_dmc

      end

      double precision wsum,w_acc_sum,wfsum,wgsum(MFORCE),wg_acc_sum,wdsum   
      double precision wgdsum, wsum1(MFORCE),w_acc_sum1,wfsum1,wgsum1(MFORCE),wg_acc_sum1
      double precision wdsum1, esum,efsum,egsum(MFORCE),esum1(MFORCE),efsum1,egsum1(MFORCE)
      double precision ei1sum,ei2sum,ei3sum,pesum(MFORCE),peisum(MFORCE),tpbsum(MFORCE),tjfsum(MFORCE),r2sum
      double precision risum,tausum(MFORCE)

      common /estsum_dmc/ wsum,w_acc_sum,wfsum,wgsum,wg_acc_sum,wdsum, &
        wgdsum, wsum1,w_acc_sum1,wfsum1,wgsum1,wg_acc_sum1,        &
        wdsum1, esum,efsum,egsum,esum1,efsum1,egsum1,              &
        ei1sum,ei2sum,ei3sum,pesum,peisum,tpbsum,tjfsum,r2sum,           &
        risum,tausum

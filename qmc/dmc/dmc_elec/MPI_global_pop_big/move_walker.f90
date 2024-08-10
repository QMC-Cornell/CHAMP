      subroutine send_walker(irecv,walker)
! Written by Claudia Filippi

# if defined (MPI)
      use all_tools_mod
      use const_mod
      use forcepar_mod
      use force_dmc_mod
      use jacobsave_mod
      use forcest_dmc_mod
      use config_dmc_mod
      use branch_mod
      use stats_mod
      use age_mod
      use velratio_mod
      use walker_type_mod
      implicit real*8(a-h,o-z)
      type(walker_t), target, intent(inout) :: walker

      dimension istatus(MPI_STATUS_SIZE)

      if (ipr.LE.-8) then
          write(6,'(''SEND walker '',1i4,'' to process '',1i3,'' with psid,psij,wt= '',1e16.8,2f16.8)') &
          nwalk,irecv,walker%psi%det,walker%psi%jas,walker%weight
          write(6,'(''x= '',100f16.8)') walker%x
      endif

      call mpi_isend(walker%weight,1,mpi_double_precision,irecv,1,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%x,size(walker%x),mpi_double_precision,irecv,2,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%age,1,mpi_integer,irecv,3,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)

      call mpi_isend(walker%psi%det,1,mpi_double_precision,irecv,4,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%jas,1,mpi_double_precision,irecv,5,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%grad,size(walker%psi%grad),mpi_double_precision,irecv,6,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%lapl,1,mpi_double_precision,irecv,7,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%hess,size(walker%psi%hess),mpi_double_precision,irecv,8,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%eloc,1,mpi_double_precision,irecv,9,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%ekinpb,1,mpi_double_precision,irecv,10,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%ekinjf,1,mpi_double_precision,irecv,11,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%epot,1,mpi_double_precision,irecv,12,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%epot_ee,1,mpi_double_precision,irecv,13,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      
      call mpi_isend(walker%psi%fsum,1,mpi_double_precision,irecv,14,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%d2,1,mpi_double_precision,irecv,15,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%fj,size(walker%psi%fj),mpi_double_precision,irecv,16,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%fs,size(walker%psi%fs),mpi_double_precision,irecv,17,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%fij,size(walker%psi%fij),mpi_double_precision,irecv,18,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%d2ij,size(walker%psi%d2ij),mpi_double_precision,irecv,19,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%lapj,size(walker%psi%lapj),mpi_double_precision,irecv,20,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%lapjij,size(walker%psi%lapjij),mpi_double_precision,irecv,21,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)

      call mpi_isend(walker%psi%orb,size(walker%psi%orb),mpi_double_precision,irecv,22,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%dorb,size(walker%psi%dorb),mpi_double_precision,irecv,23,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%ddorb,size(walker%psi%ddorb),mpi_double_precision,irecv,24,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%aiup,size(walker%psi%aiup),mpi_double_precision,irecv,25,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%aidn,size(walker%psi%aidn),mpi_double_precision,irecv,26,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%tup,size(walker%psi%tup),mpi_double_precision,irecv,27,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%tdn,size(walker%psi%tdn),mpi_double_precision,irecv,28,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%detup,size(walker%psi%detup),mpi_double_precision,irecv,29,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%detdn,size(walker%psi%detdn),mpi_double_precision,irecv,30,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%invup,size(walker%psi%invup),mpi_double_precision,irecv,31,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%invdn,size(walker%psi%invdn),mpi_double_precision,irecv,32,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%yup,size(walker%psi%yup),mpi_double_precision,irecv,33,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%ydn,size(walker%psi%ydn),mpi_double_precision,irecv,34,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%deta_up,1,mpi_double_precision,irecv,35,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%deta_dn,1,mpi_double_precision,irecv,36,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%chi,1,mpi_double_precision,irecv,37,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      
      call mpi_isend(walker%psi%rvec_en,size(walker%psi%rvec_en),mpi_double_precision,irecv,38,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%rvec_ee,size(walker%psi%rvec_ee),mpi_double_precision,irecv,39,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%quadr,size(walker%psi%quadr),mpi_double_precision,irecv,40,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%quadx,size(walker%psi%quadx),mpi_double_precision,irecv,41,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%iwfragelec,size(walker%psi%iwfragelec),mpi_double_precision,irecv,42,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%enefrag,size(walker%psi%enefrag),mpi_double_precision,irecv,43,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%r_en,size(walker%psi%r_en),mpi_double_precision,irecv,44,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)
      call mpi_isend(walker%psi%r_ee,size(walker%psi%r_ee),mpi_double_precision,irecv,45,MPI_COMM_WORLD,irequest,ierr)
      call MPI_Request_free(irequest,ierr)

!      call mpi_isend(wt(nwalk),1,mpi_double_precision,irecv,1 &
!     &,MPI_COMM_WORLD,irequest,ierr)
!      call MPI_Request_free(irequest,ierr)
!      call mpi_isend(iage(nwalk),1,mpi_integer,irecv,2 &
!     &,MPI_COMM_WORLD,irequest,ierr)
!      call MPI_Request_free(irequest,ierr)
!
!! MPI_Request_free is faster than MPI_wait and is OK if there is a barrier
!! before variables are reused.
!      itag=2
!      do 15 ifr=1,nforce
!        call mpi_isend(ajacold(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+1,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(eoldw(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+2,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(psidow(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+3,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(psijow(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+4,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(peow(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+5,MPI_COMM_WORLD,irequest,ierr)
!        call mpi_isend(peiow(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+5,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(d2ow(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+6,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(pwt(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+7,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(fratio(nwalk,ifr),1,mpi_double_precision,irecv &
!     &  ,itag+8,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(voldw(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
!     &  ,irecv,itag+9,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!        call mpi_isend(xoldw(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
!     &  ,irecv,itag+10,MPI_COMM_WORLD,irequest,ierr)
!        itag=itag+10
!        call MPI_Request_free(irequest,ierr)
!        do 15 ip=0,nwprod-1
!        itag=itag+1
!        call mpi_isend(wthist(nwalk,ip,ifr),1,mpi_double_precision,irecv &
!     &  ,itag,MPI_COMM_WORLD,irequest,ierr)
!        call MPI_Request_free(irequest,ierr)
!   15 continue
!
!!     call send_det(itag,irecv)
!!     call send_jas(itag,irecv)
!
!!     nwalk=nwalk-1

      return

      entry recv_walker(isend,walker)

      call mpi_recv(walker%weight,1,mpi_double_precision,isend,1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%x,size(walker%x),mpi_double_precision,isend,2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%age,1,mpi_integer,isend,3,MPI_COMM_WORLD,istatus,ierr)

      call mpi_recv(walker%psi%det,1,mpi_double_precision,isend,4,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%jas,1,mpi_double_precision,isend,5,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%grad,size(walker%psi%grad),mpi_double_precision,isend,6,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%lapl,1,mpi_double_precision,isend,7,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%hess,size(walker%psi%hess),mpi_double_precision,isend,8,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%eloc,1,mpi_double_precision,isend,9,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%ekinpb,1,mpi_double_precision,isend,10,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%ekinjf,1,mpi_double_precision,isend,11,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%epot,1,mpi_double_precision,isend,12,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%epot_ee,1,mpi_double_precision,isend,13,MPI_COMM_WORLD,istatus,ierr)
      
      call mpi_recv(walker%psi%fsum,1,mpi_double_precision,isend,14,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%d2,1,mpi_double_precision,isend,15,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%fj,size(walker%psi%fj),mpi_double_precision,isend,16,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%fs,size(walker%psi%fs),mpi_double_precision,isend,17,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%fij,size(walker%psi%fij),mpi_double_precision,isend,18,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%d2ij,size(walker%psi%d2ij),mpi_double_precision,isend,19,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%lapj,size(walker%psi%lapj),mpi_double_precision,isend,20,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%lapjij,size(walker%psi%lapjij),mpi_double_precision,isend,21,MPI_COMM_WORLD,istatus,ierr)

      call mpi_recv(walker%psi%orb,size(walker%psi%orb),mpi_double_precision,isend,22,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%dorb,size(walker%psi%dorb),mpi_double_precision,isend,23,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%ddorb,size(walker%psi%ddorb),mpi_double_precision,isend,24,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%aiup,size(walker%psi%aiup),mpi_double_precision,isend,25,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%aidn,size(walker%psi%aidn),mpi_double_precision,isend,26,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%tup,size(walker%psi%tup),mpi_double_precision,isend,27,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%tdn,size(walker%psi%tdn),mpi_double_precision,isend,28,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%detup,size(walker%psi%detup),mpi_double_precision,isend,29,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%detdn,size(walker%psi%detdn),mpi_double_precision,isend,30,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%invup,size(walker%psi%invup),mpi_double_precision,isend,31,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%invdn,size(walker%psi%invdn),mpi_double_precision,isend,32,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%yup,size(walker%psi%yup),mpi_double_precision,isend,33,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%ydn,size(walker%psi%ydn),mpi_double_precision,isend,34,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%deta_up,1,mpi_double_precision,isend,35,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%deta_dn,1,mpi_double_precision,isend,36,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%chi,1,mpi_double_precision,isend,37,MPI_COMM_WORLD,istatus,ierr)
      
      call mpi_recv(walker%psi%rvec_en,size(walker%psi%rvec_en),mpi_double_precision,isend,38,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%rvec_ee,size(walker%psi%rvec_ee),mpi_double_precision,isend,39,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%quadr,size(walker%psi%quadr),mpi_double_precision,isend,40,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%quadx,size(walker%psi%quadx),mpi_double_precision,isend,41,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%iwfragelec,size(walker%psi%iwfragelec),mpi_double_precision,isend,42,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%enefrag,size(walker%psi%enefrag),mpi_double_precision,isend,43,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%r_en,size(walker%psi%r_en),mpi_double_precision,isend,44,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(walker%psi%r_ee,size(walker%psi%r_ee),mpi_double_precision,isend,45,MPI_COMM_WORLD,istatus,ierr)

!     nwalk=nwalk+1

!      call mpi_recv(wt(nwalk),1,mpi_double_precision,isend,1 &
!     &,MPI_COMM_WORLD,istatus,ierr)
!      call mpi_recv(iage(nwalk),1,mpi_integer,isend,2 &
!     &,MPI_COMM_WORLD,istatus,ierr)
!
!      itag=2
!      do 25 ifr=1,nforce
!        call mpi_recv(ajacold(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+1,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(eoldw(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+2,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(psidow(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+3,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(psijow(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+4,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(peow(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+5,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(peiow(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+5,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(d2ow(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+6,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(pwt(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+7,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(fratio(nwalk,ifr),1,mpi_double_precision,isend &
!     &  ,itag+8,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(voldw(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
!     &  ,isend,itag+9,MPI_COMM_WORLD,istatus,ierr)
!        call mpi_recv(xoldw(1,1,nwalk,ifr),3*nelec,mpi_double_precision &
!     &  ,isend,itag+10,MPI_COMM_WORLD,istatus,ierr)
!        itag=itag+10
!        do 25 ip=0,nwprod-1
!        itag=itag+1
!        call mpi_recv(wthist(nwalk,ip,ifr),1,mpi_double_precision,isend &
!     &  ,itag,MPI_COMM_WORLD,istatus,ierr)
!   25 continue

!     call recv_det(itag,isend)
!     call recv_jas(itag,isend)
      if (ipr.LE.-8) then
          write(6,'(''RECV walker '',1i4,'' from process '',1i3,'' with psid,psij,wt= '',1e16.8,2f16.8)') &
          nwalk,isend,walker%psi%det,walker%psi%jas,walker%weight
          write(6,'(''x= '',100f16.8)') walker%x
      endif

# endif
      return
      end

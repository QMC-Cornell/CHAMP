      subroutine detsav(iel)
! Written by Claudia Filippi, modified by Cyrus Umrigar

      use dete_mod
      use dorb_mod
      use slatn_mod
      use orbe_mod
      use dets_mod
      use slater_mod
      implicit real*8(a-h,o-z)

      orb(iel,:) = orbe(:) !TA
      dorb(1,iel,:) = dorbe(1,:)
      dorb(2,iel,:) = dorbe(2,:)
      dorb(3,iel,:) = dorbe(3,:)
      ddorb(iel,:) = ddorbe(:)

      return
      end

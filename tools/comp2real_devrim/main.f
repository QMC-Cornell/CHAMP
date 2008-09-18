***********************************************************************
      program psi_comp2real
c Written by Devrim Guclu for Cyrus Umrigar
*----------------------------------------------------------------------
* title: Psi complexe to real
* file: main.f filem.f order.f comp2real.f
* author: A.D.Guclu
* date: nov 2003
* modification:
* description: finds the real and imag. parts of a linear comb. of
*              Slater determinants, by converting the complex orbitals
*              into real orbitals. The results can be used for QMC
*              calculations
*----------------------------------------------------------------------
* LIST OF VARIABLES:
*
* job_name = job_name used to define input/output files
* ndim     = number of spatial dimensions
* nelec    = number of electrons
* ndet     = number of determinants (which increases during expansion)
* det(idet,in,qnx) = quantum number "qnx" of in^{th} electron in the
*                    idet^{th} determinant
* cdet(idet)    = coefficient of the idet^{th} determinant.
*----------------------------------------------------------------------

      include         'maxdim.h'
      include         'qnumbers.h'

      character*35    job_name

      integer         ndim,nelec,ndet
      integer         det(MAXNDET,MAXNELEC,MAXQN)
      complex*16      cdet(MAXNDET)

      common/psi_n/   cdet,ndim,nelec,ndet,det

      call get_job_name(job_name)

      call read_input(job_name)

c update the sign of coefficents for 3D system
      call sign_3D

c define an ordered list of single electron states (orbitals)
      call def_ses

c the big expansion:
      call expand_det

c order the orbitals within each determinant, and get rid of
c null determinants due to Pauli pr.
      call order_det

c get rid of doubled determinants:
      call reduce_det

      call write_output(job_name)

      stop
      end

      include 'constants.f'
      real(8) :: R,rho,pi0,pi1(8),pis(group)
      common/para_1/R,rho,pi0,pi1,pis
      logical :: H0
      common/h0_vs_h1/H0
      real(8) :: m(0:2,group),n(0:1,group)
      common/para_2/m,n
      character(len=3) :: model
      common/para_3/model
      integer :: mmt(8),nnt(8),mtot,ntot,tot
      common/mi_ni/mmt,nnt,mtot,ntot,tot
      integer :: m_vec(group),n_vec(group)
      common/mn_vec/m_vec,n_vec
      integer :: seed
      common/seed/seed

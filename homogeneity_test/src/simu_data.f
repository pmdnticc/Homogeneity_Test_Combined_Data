c$$$      program simu
c$$$      implicit none
c$$$!      include 'constants.f'
c$$$      include 'parameters.f'
c$$$      integer :: i
c$$$      call read_input
c$$$      call simu_data
c$$$      do i=1,group
c$$$         write(*,*) 'group ', i
c$$$         write(*,*) 'm is: ', m(:,i)
c$$$         write(*,*) 'n is: ', n(:,i)
c$$$      end do
c$$$      end program simu
      

      subroutine simu_data
      implicit none
!      include 'constants.f'
      include 'parameters.f'
      real(8) :: prob(0:2,group),rr,omr!,aa(mp0),bb(np0)
      integer :: i,j,k
      real(8), allocatable :: aa(:),bb(:)
c$$$      integer :: sn
c$$$      integer,allocatable :: seed(:)

      if(H0) then
         pis=pi0
      else
         pis(1:group)=pi1(1:group)
      end if
      
c$$$      call random_seed(size=sn)
c$$$      allocate(seed(sn))
c$$$      call random_seed(get=seed)
c$$$      write(*,*) 'show seed: ', seed
      select case(model)
      case('Ros')
         rr=R
c$$$         prob(0,:)=R*pis(:)**2-2.d0*pis(:)+1.d0
c$$$         prob(1,:)=2.d0*pis(:)*(1.d0-R*pis(:))
c$$$         prob(2,:)=R*pis(:)**2
      case('Don')
         rr=rho
c$$$         omr=1.d0-rho
c$$$         prob(0,:)=rho*(1.d0-pis(:))+omr*(1.d0-pis(:))**2
c$$$         prob(1,:)=2.d0*omr*pis(:)*(1.d0-pis(:))
c$$$         prob(2,:)=rho*pis(:)+omr*pis(:)**2
      case default
         write(*,*) "models haven't been implemented ..."//
     .        'routine simu_data'
         return
      end select
      call prob_gen(prob,rr)
c$$$      print*, 'generated prob: '
c$$$      print*, prob(0,:)
c$$$      print*, prob(1,:)
c$$$      print*, prob(2,:)
C$$$      Seed=111111
c$$$      call random_seed(put=seed)
      m=0
      n=0
      do i=1,group
 100     continue
         if(allocated(aa)) deallocate(aa)
         if(allocated(bb)) deallocate(bb)
         allocate(aa(m_vec(i)))
         allocate(bb(n_vec(i)))
         call random_number(aa) ! gen random #'s for each group
         call random_number(bb)
         do j=1,m_vec(i)!mp0
            if(aa(j)<=prob(0,i)) then
               m(0,i)=m(0,i)+1
            else if(aa(j)>prob(0,i)+prob(1,i)) then
               m(2,i)=m(2,i)+1
            end if
         end do
!         m(1,i)=mp(i)-m(0,i)-m(2,i)
         m(1,i)=m_vec(i)-m(0,i)-m(2,i)
         if(m(1,i)<0) then
            m(:,i)=0
            go to 100
         end if
         do j=0,2
            if(m(j,i)==0) m(j,i)=0.5d0
         end do
         do k=1,n_vec(i)!np0
            if(bb(k)<=1.d0-pis(i)) n(0,i)=n(0,i)+1
         end do
!         n(1,i)=np(i)-n(0,i)
         n(1,i)=n_vec(i)-n(0,i)
         do k=0,1
            if(n(k,i)==0) n(k,i)=0.5d0
         end do
      end do
      deallocate(aa)
      deallocate(bb)
      return
      end subroutine simu_data




!     generate probabilities of bilateral data based on models
      subroutine prob_gen(prob,rr)
      implicit none
      include 'parameters.f'
      real(8) :: prob(0:2,group),rr,omr
      integer :: i,j
      select case(model)
      case('Ros')
         prob(0,:)=rr*pis(:)**2-2.d0*pis(:)+1.d0
         prob(1,:)=2.d0*pis(:)*(1.d0-rr*pis(:))
         prob(2,:)=rr*pis(:)**2
      case('Don')
         omr=1.d0-rr
         prob(0,:)=rr*(1.d0-pis(:))+omr*(1.d0-pis(:))**2
         prob(1,:)=2.d0*omr*pis(:)*(1.d0-pis(:))
         prob(2,:)=rr*pis(:)+omr*pis(:)**2
      case default
         write(*,*) "models haven't been implemented ..."//
     .        'routine prob_gen'
         return
      end select
      do i=1,group
         do j=0,2
            if(prob(j,i)==0.d0)
     .           prob(j,i)=1.d0-prob(mod(j+1,3),i)-prob(mod(j+2,3),i)
         end do
      end do
      end subroutine prob_gen

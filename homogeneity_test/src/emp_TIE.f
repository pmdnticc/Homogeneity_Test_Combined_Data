!     compute empirical type I err rate
      program emp_TIE
      implicit none
      include 'parameters.f'
      real(8) :: rej(3),test(3),pv(3),df,tie(3)
      integer :: i,j,ist
c$$$      model='Don'
c$$$      rho=0.4d0
c$$$      R=(1.d0-pi0)*rho/pi0+1.d0
      call read_input
      pis=pi0
      df=dble(group-1)
      rej=0.d0
      do i=1,iter
 200     call simu_data
         call TS(test,ist)
         if(ist/=0) go to 200
         do j=1,3
            call get_pv(pv(j),test(j),df)
            if(pv(j)<alf) rej(j)=rej(j)+1.d0
         end do
      end do
      tie=rej/iter
      write(*,*) 'The empirical TIE rate is (LR, Wald, Score): '
      write(*,*) tie
      return
      end program emp_TIE

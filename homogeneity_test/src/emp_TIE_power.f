!     compute empirical type I err rate under H0
!     & power under H1
      program emp_TIE_power
      implicit none
      include 'parameters.f'
      include 'head_tail.f'
      real(8) :: rej(3),test(3),pv(3),df,tie(3)
      integer :: i,j,ist(3),seed2,seed3,count(3)
      character(len=2) :: hyp
      character(len=3) :: iseed
c$$$      model='Don'
c$$$      rho=0.4d0
c$$$      R=(1.d0-pi0)*rho/pi0+1.d0
      call read_input
      hyp='H1'
      if(H0) hyp='H0'
      write(iseed,'(I3)') seed
      seed2=seed+1
      seed3=seed+2
      if(H0) then
         pis=pi0
      else
         pis(1:group)=pi1(1:group)
      end if
      
!      pis=pi0
      df=dble(group-1)
      rej=0.d0
      header=.true.
      tailer=.false.
c$$$      call system('mkdir -pv ../data/sas/run_'//iseed)
c$$$      call system('mkdir -pv ../data/fort/run_'//iseed)
c$$$      open(unit=seed,file='../data/sas/run_'//iseed
c$$$     .     //'/run_'//iseed//'.sas',status='unknown')
c$$$      open(unit=seed2,file='../data/fort/run_'//iseed
c$$$     .     //'/run_'//iseed//'.1',status='unknown')
c$$$      open(unit=seed3,file='../data/fort/run_'//iseed
c$$$     .     //'/run_'//iseed//'.2',status='unknown')
      call system('mkdir -pv ../data/sas')
      call system('mkdir -pv ../data/fort')
      open(unit=seed,file='../data/sas/run_'//iseed//'_'//model
     .     //'_'//hyp//'.sas',status='unknown')
      open(unit=seed2,file='../data/fort/run_'//iseed//'_'//model
     .     //'_'//hyp//'.1',status='unknown')
      if(H0) 
     .     open(unit=seed3,file='../data/fort/run_'//iseed//'_'//model
     .     //'_'//hyp//'.2',status='unknown')
      call data_to_sas(seed)
      if(H0) then
         write(seed2,*) "#Empirical Type I Error Rate (alf=0.05)"
      else
         write(seed2,*) "#Power (alf=0.05)"
      end if
      if(H0) write(seed3,*) "#*********** P-values **********"
      write(seed2,*) "#  LR        Wald      Score"
      if(H0) write(seed3,*) "#  LR        Wald      Score"
      header=.false.
      count=iter
      do i=1,iter
         call simu_data
!         print*, "m === ", m
!         print*, "n === ", n
         call TS(test,ist)
!         print*, "test scores: ", test
         call data_to_sas(seed)
         do j=1,3
            if(ist(j)/=0) then
               count(j)=count(j)-1
               cycle
            end if
            call get_pv(pv(j),test(j),df)
            if(pv(j)<alf) then
               rej(j)=rej(j)+1.d0
            end if
         end do
         if(H0) write(seed3,'(3(F10.4))') test
      end do
      tailer=.true.
      call data_to_sas(seed)
      call data_stack_sas(seed)
      call invoke_genmod(seed)
      tie=rej/count
      print*, "count: ", count
      write(seed2,'(3(F10.4))') tie
      if(H0) then
         write(*,*) 'The empirical TIE rate is (LR, Wald, Score): '
      else
         write(*,*) 'The power is (LR, Wald, Score): '
      end if
      write(*,*) tie
      
      close(seed)
      close(seed2)
      if(H0) close(seed3)
      return
      end program emp_TIE_power

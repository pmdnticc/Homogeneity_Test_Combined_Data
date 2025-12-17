!     real data analysis
!     TWO examples with group=2
!     Ex.1: OME by Mandel et al. (1982)
!     Ex.2: 2nd example from Liang et al. (2023)
!     [study conducted at the First Affiliated Hospital of Xiamen University]
!     ONE example with group=4
!     Ex. ocular exam by Rosner (1982)
      program real_example
      implicit none
      include 'parameters.f'
      real(8) :: res(4)
      select case(group)
      case(2)
!     data structure for Ex.1
         m(0,1)=9
         m(0,2)=7
         m(1,1)=7
         m(1,2)=5
         m(2,1)=23
         m(2,2)=13
         n(0,1)=20
         n(0,2)=19
         n(1,1)=34
         n(1,2)=36
         model='Ros'
         call get_statistic
         call goodness_of_fit(res,'com')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
         model='Don'
         call get_statistic
         call goodness_of_fit(res,'com')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
!     data structure for Ex.2
         m(0,1)=20
         m(0,2)=7
         m(1,1)=10
         m(1,2)=13
         m(2,1)=2
         m(2,2)=2
         n(0,1)=3
         n(0,2)=3
         n(1,1)=0
         n(1,2)=0
         model='Ros'
         call get_statistic
         call goodness_of_fit(res,'com')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
         model='Don'
         call get_statistic
         call goodness_of_fit(res,'com')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
!     Ex.3&4 for goodness-of-fit
!     data source: Goodness-of-fit tests for correlated paired binary data by Tang et al. (2010)
!     Table 6 (leg)
         write(*,*) ">>>>>>>>>>>>>>>> X-check AIC <<<<<<<<<<<<<<<<<<"
         m(0,1)=6
         m(1,1)=1
         m(2,1)=32
         m(0,2)=9
         m(1,2)=2
         m(2,2)=15
         n=0
         model='Ros'
         call get_statistic
         call goodness_of_fit(res,'bin')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
         model='Don'
         call get_statistic
         call goodness_of_fit(res,'bin')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
!     Table 7 (eye)
         m(0,1)=8
         m(1,1)=1
         m(2,1)=1
         m(0,2)=0
         m(1,2)=1
         m(2,2)=2
         n=0
         model='Ros'
         call get_statistic
         call goodness_of_fit(res,'bin')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
         model='Don'
         call get_statistic
         call goodness_of_fit(res,'bin')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'         
      case(4)
!     data structure for Ex w/ group=4
         m(0,1)=15
         m(0,2)=7
         m(0,3)=3
         m(0,4)=67
         m(1,1)=6
         m(1,2)=5
         m(1,3)=2
         m(1,4)=24
         m(2,1)=7
         m(2,2)=9
         m(2,3)=14
         m(2,4)=57
         n=0
         model='Ros'
         call get_statistic
         call goodness_of_fit(res,'bin')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
         model='Don'
         call get_statistic
         call goodness_of_fit(res,'bin')
         write(*,*) '>>>>>>>>>>>> Goodness Of Fit <<<<<<<<<<<<'
         write(*,*) 'G2, X2, X2_adj, AIC: ', res
         write(*,*) '>>>>>>>>>>>> End <<<<<<<<<<<<'
      end select
      return
      end program real_example





      subroutine get_statistic
      implicit none
      include 'parameters.f'
      real(8) :: test(3),pv(3),df,rr0,rr,pi00
      integer :: i,j,ist
      df=dble(group-1)
      write(*,*) ">>>>>>>>>>>>>>> START <<<<<<<<<<<<<<<"
      write(*,*) "Under Model: ", model
      call mle0(pi00,rr0)
      if(model=='Ros') then
         write(*,*) "pi0, R0 = ", pi00,rr0
         write(*,*) "Rho0 = ", pi00*(rr0-1.d0)/(1.d0-pi00)
      else if(model=='Don') then
         write(*,*) "pi0, Rho0 = ", pi00,rr0
      end if
      call mle(pis,rr)
      write(*,*) "pis = ", (pis(i),i=1,group)
      if(model=='Ros') then
         write(*,*) "R = ", rr
         write(*,*) "Rho = ", (pis(i)*(rr-1.d0)/(1.d0-pis(i)),i=1,group)
      else if(model=='Don') then
         write(*,*) "Rho = ", rr
      end if
      call TS(test,ist)
      do j=1,3
         call get_pv(pv(j),test(j),df)
      end do
      write(*,*) "Test Statistics: LR, Wald, Score"
      write(*,*) test
      write(*,*) "P-values: "
      write(*,*) pv
      write(*,*) ">>>>>>>>>>>>>>> END <<<<<<<<<<<<<<<"
      return
      end subroutine get_statistic



!     goodness-of-fit: return p-values for G2,X2,X2_adj, &AIC
      subroutine goodness_of_fit(res,tag)
      implicit none
      include 'parameters.f'
      character(len=3) :: tag !tag='com' for combined data; otherwise if bilateral only
      real(8) :: res(4),df,rr0,rr,pi00,prob(0:2,group)
      real(8) :: Em(0:2,group),En(0:1,group),mi(group),ni(group)
      real(8) :: Gsq,Xsq,Xsq_adj,AIC
      integer :: i,j
      real(8) :: log_fac
      external log_fac
      df=dble(group-1)
      call mle0(pi00,rr0)
      call mle(pis,rr)
      call prob_gen(prob,rr)
      mi(:)=m(0,:)+m(1,:)+m(2,:)
      ni(:)=n(0,:)+n(1,:)
      do i=1,group
         Em(0:2,i)=mi(i)*prob(0:2,i)
         En(1,i)=ni(i)*pis(i)
         En(0,i)=ni(i)*(1.d0-pis(i))
      end do
      Gsq=0.d0
      Xsq=0.d0
      Xsq_adj=0.d0
      AIC=0.d0
!     LR test stat: Gsq
      do i=1,group
         do j=0,2
            if(m(j,i)>0.d0) Gsq=Gsq+2.d0*m(j,i)*log(m(j,i)/Em(j,i))
            Xsq=Xsq+(m(j,i)-Em(j,i))**2/Em(j,i)
            Xsq_adj=Xsq_adj+(abs(m(j,i)-Em(j,i))-0.5d0)**2/Em(j,i)
            AIC=AIC+m(j,i)*log(prob(j,i))
         end do
!         AIC=AIC+log(gamma(mi(i)+1)/gamma(m(0,i)+1)/gamma(m(1,i)+1)
!     .        /gamma(m(2,i)+1))
!         AIC=AIC+log_fac(mi(i))-log_fac(m(0,i))-log_fac(m(1,i))
!     .        -log_fac(m(2,i))
      end do
      if(tag=='com') then
         do i=1,group
            do j=0,1
               if(n(j,i)>0.d0) Gsq=Gsq+2.d0*n(j,i)*log(n(j,i)/En(j,i))
               Xsq=Xsq+(n(j,i)-En(j,i))**2/En(j,i)
               Xsq_adj=Xsq_adj+(abs(n(j,i)-En(j,i))-0.5d0)**2/En(j,i)
            end do
            AIC=AIC+n(0,i)*log(1.d0-pis(i))+n(1,i)*log(pis(i))
!     .           +log(gamma(ni(i)+1)/gamma(n(0,i)+1)/gamma(n(1,i)+1))
!     .           +log_fac(ni(i))-log_fac(n(0,i))-log_fac(n(1,i))
         end do
      end if
      AIC=2.d0*(group+1)-2.d0*AIC
      call get_pv(res(1),Gsq,df)
      call get_pv(res(2),Xsq,df)
      call get_pv(res(3),Xsq_adj,df)
      res(4)=AIC
      return
      end subroutine goodness_of_fit



!     function calculating log(x!) using sum(log(i)) (i=1,...,x)
!     avoid using factorial directly in case of x! = infty numerically
      real(8) function log_fac(x)
      implicit none
      real(8) :: x,sum
      integer :: i
      sum=0.d0
      do i=1,int(x)
         sum=sum+log(dble(i))
      end do
      log_fac=sum
      end function log_fac

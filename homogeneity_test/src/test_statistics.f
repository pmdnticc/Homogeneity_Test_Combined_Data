!     test statistics: LRT, Wald, Score
      subroutine TS(test,ist)
      implicit none
      include 'parameters.f'
      include 'info.f'
      real(8) :: test(3),ll,rr,ll0,rr0,pi00
      real(8) :: fi(group+1,group+1)
      integer :: ist(3),i
      ist=0
      test=0.d0
      call mle0(pi00,rr0)
!      print*, "pi00 is ", pi00
!     one can rely on this step of call mle0 that returns rr<0
!     w/ a flag ist=-1, since within routine mle, it also calls
!     mle0 at the 1st step. Iniitial rr<0 if same configurations of m & n.
      if(rr0<0.d0.and.model=='Ros'.or.
     .     abs(rr0)>1.d0.and.model=='Don') then
         ist(1)=-1
         ist(3)=-1
      end if
!      if(model=='Don'.and.sign(1.d0,rho)/=sign(1.d0,rr0)) then
!         ist(1)=-1
!         ist(3)=-1
!      end if
      pis=pi00
      call loglikelihood(ll0,rr0)
!      print*, 'll0, rr0 ===> ', ll0,rr0
      call mle(pis,rr)
      if(rr<0.d0.and.model=='Ros'.or.
     .     abs(rr)>1.d0.and.model=='Don') then
         ist(1)=-1
         ist(2)=-1
      end if
!      print*, "pis is: ", pis
      call loglikelihood(ll,rr)
!      print*, 'll, rr ===> ', ll,rr
!     LRT
      test(1)=2.d0*(ll-ll0)
      call fisher_info(fi,pis,rr)
      if(info/=0) ist(2)=-1
!      print*, 'fisher: ', fi
!     Wald
      call wald_test(test(2),pis,rr,fi)
      if(info/=0) ist(2)=-1
!     Score
      pis=pi00
      call fisher_info(fi,pis,rr0)
      if(info/=0) ist(3)=-1
!      print*, 'pis: ', pis
!      print*, 'rr0: ', rr0
      call score_test(test(3),pis,rr0,fi)
      if(info/=0) ist(3)=-1
!      print*, 'Score: ', test(3)
!      if(test(1)<0.d0.or.test(2)<0.d0.or.test(3)<0.d0) go to 101
      do i=1,3
         if(test(i)<0.d0.or.test(i)/=test(i)) ist(i)=-1
      end do
      return
 101  ist=-1
      end subroutine TS



!     log-likelihood
      subroutine loglikelihood(ll,rr)
      implicit none
      include 'parameters.f'
      real(8) :: ll,rr,prob(0:2,group)
      integer :: i
      call prob_gen(prob,rr)
      ll=0.d0
      do i=1,group
!         print*, "p_ri: ",prob(0:2,i)
         ll=ll+m(0,i)*log(prob(0,i))+m(1,i)*log(prob(1,i))
     .        +m(2,i)*log(prob(2,i))+n(0,i)*log(1.d0-pis(i))
     .        +n(1,i)*log(pis(i))
      end do
      return
      end subroutine loglikelihood


!     use routine cdfchi from cdflib to get p-value based on chisq-dist
      subroutine get_pv(pv,x,df)
      implicit none
      real(8) :: pv,x,df,p,q,bound
      integer :: index,status
      parameter(index=1)
      call cdfchi(index,p,q,x,df,status,bound)
      if(status/=0) then
         write(*,*) 'something is not right ...'
         write(*,*) 'x is ', x
         write(*,*) 'cdfchi returns STATUS: ', status
         stop
      end if
      pv=q
      return
      end subroutine get_pv




!     return Wald ts
      subroutine wald_test(tw,ph,rh,fi)
      implicit none
      include 'parameters.f'
      include 'info.f'
      real(8) :: tw,ph(group),rh
      real(8) :: fi(group+1,group+1),ifi(group+1,group+1)
      real(8) :: beta(group+1),Cmatrix(group-1,group+1)
      real(8) :: Kmatrix(group-1,group-1),invK(group-1,group-1),kij
      integer :: i,j,i1,i2
      interface
         function inv(A) result(Ainv)
         real(8), dimension(:,:), intent(in) :: A
         real(8), dimension(size(A,1),size(A,2)) :: Ainv
         end function inv
      end interface

      beta(1:group)=ph(1:group)
      beta(group+1)=rh
      Cmatrix=0.d0
      do i=1,group-1
         Cmatrix(i,i)=1.d0
         Cmatrix(i,i+1)=-1.d0
      end do
!      call fisher_info(fi,ph,rh)
      ifi=inv(fi)
      do i=1,group-1
         do j=1,group-1
            kij=0.d0
            do i1=1,group+1
               do i2=1,group+1
                  kij=kij+Cmatrix(i,i1)*ifi(i1,i2)*Cmatrix(j,i2)
               end do
            end do
            Kmatrix(i,j)=kij
         end do
      end do
      invK=inv(Kmatrix)
      tw=0.d0
      do i=1,group+1
         do i1=1,group-1
            do i2=1,group-1
               do j=1,group+1
                  tw=tw+beta(i)*Cmatrix(i1,i)*
     .                 invK(i1,i2)*Cmatrix(i2,j)*beta(j)
               end do
            end do
         end do
      end do
      return
      end subroutine wald_test




!     return Score ts
      subroutine score_test(ts,ph,rh,fi)
      implicit none
      include 'parameters.f'
      include 'info.f'
      real(8) :: ts,ph(group),rh
      real(8) :: U_vec(group+1)
      real(8) :: fi(group+1,group+1),ifi(group+1,group+1)
      integer :: i,j
      interface
         function inv(A) result(Ainv)
         real(8), dimension(:,:), intent(in) :: A
         real(8), dimension(size(A,1),size(A,2)) :: Ainv
         end function inv
      end interface
      U_vec=0.d0
      do i=1,group
         if(model=='Ros') then
            U_vec(i)=2.d0*m(0,i)*(rh*ph(i)-1.d0)
     .           /(rh*ph(i)**2-2.d0*ph(i)+1.d0)
     .           +m(1,i)*(2.d0*rh*ph(i)-1.d0)/ph(i)/(rh*ph(i)-1)
     .           +2.d0*m(2,i)/ph(i)
     .           -n(0,i)/(1.d0-ph(i))
     .           +n(1,i)/ph(i)
         else if(model=='Don') then
            U_vec(i)=m(0,i)*(2.d0+2.d0*rh*ph(i)-2.d0*ph(i)-rh)
     .           /(ph(i)-1.d0)/(rh*ph(i)-ph(i)+1.d0)
     .           +m(1,i)*(2.d0*ph(i)-1.d0)/ph(i)/(ph(i)-1.d0)
     .           +m(2,i)*(rh+2.d0*ph(i)-2.d0*rh*ph(i))
     .           /ph(i)/(rh+ph(i)-rh*ph(i))
     .           -n(0,i)/(1.d0-ph(i))
     .           +n(1,i)/ph(i)
         else
            stop "model hasn't been implemented ... score_test"
         end if
      end do
!      call fisher_info(fi,ph,rh)
      ifi=inv(fi)
      ts=0.d0
      do i=1,group+1
         do j=1,group+1
            ts=ts+U_vec(i)*ifi(i,j)*U_vec(j)
         end do
      end do
      return
      end subroutine score_test



c$$$!     test routine wald_test & score_test
c$$$      program test
c$$$      implicit none
c$$$      include 'parameters.f'
c$$$      real(8) :: rr,fi(group+1,group+1),tw,ts
c$$$      model='Ros'
c$$$      rho=0.4d0
c$$$      R=(1.d0-pi0)*rho/pi0+1.d0
c$$$      pis=pi0
c$$$      call simu_data
c$$$      call mle(pis,rr)
c$$$      call fisher_info(fi,pis,rr)
c$$$      call wald_test(tw,pis,rr,fi)
c$$$      print*, 'Tw2 = ', tw
c$$$      call score_test(ts,pis,rr,fi)
c$$$      print*, 'Ts2 = ', ts
c$$$      end program test
      

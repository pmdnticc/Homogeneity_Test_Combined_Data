c$$$!     test routine mle0 & mle
c$$$      program test_mle
c$$$      implicit none
c$$$      include 'parameters.f'
c$$$      real(8) :: ph,rh
c$$$      model='Don'
c$$$      rho=0.2d0
c$$$      R=1.4d0
c$$$      pis=pi0
c$$$      call simu_data
c$$$      call mle0(ph,rh)
c$$$      write(*,*) 'returned pi & r by mle0: ', ph, rh
c$$$      call mle(pis,rh)
c$$$      write(*,*) 'returned pi & r by mle : ', pis, rh
c$$$      return
c$$$      end program test_mle

      


      
!     calculate constrained (routine mle0)
!     and unconstraind MLE (routine mle)
!     expressions are of same structure for Rosner & Donner's model
      subroutine mle0(pihat,rhat)
      implicit none
      include 'parameters.f'
      integer :: m0,m1,m2,n0,n1,mt,nt,nn,i
      real(8) :: pihat,rhat,a,c,theta
      m0=0
      m1=0
      m2=0
      n0=0
      n1=0
      do i=1,group
         m0=m0+m(0,i)
         m1=m1+m(1,i)
         m2=m2+m(2,i)
         n0=n0+n(0,i)
         n1=n1+n(1,i)
      end do
      mt=m0+m1+m2
      nt=n0+n1
      nn=mt+nt
      a=n0+5.d0*n1+2.d0*m0+3.d0*m1+4.d0*m2
      c=m0*(3.d0*n1+m1+2.d0*m2)+m1*(m1+3.d0*m2+n0)+m2*(2.d0*m2+n0)
     .     +n1*(2.d0*n0+4.d0*n1+5.d0*m1+6.d0*m2)
      theta=18.d0*a*c-2.d0*a**3-108.d0*nn*n1*(n1+m1+m2)
      theta=theta/2.d0/dsqrt((a**2-6.d0*c)**3)
      theta=dacos(theta)/3.d0
      pihat=a
     .     +dsqrt(a**2-6.d0*c)*(dcos(theta)-dsqrt(3.d0)*dsin(theta))
      pihat=pihat/6.d0/nn
c$$$      if(pihat<1.d-3) then
c$$$         print*, 'pi = 0 occurs'
c$$$         print*, 'm0: ', m(0,:)
c$$$         print*, 'm1: ', m(1,:)
c$$$         print*, 'm2: ', m(2,:)
c$$$         print*, 'n0: ', n(0,:)
c$$$         print*, 'n1: ', n(1,:)
c$$$      end if
      if(model=='Ros') then
         rhat=2.d0*nn*pihat**2-pihat*(2.d0*mt+n0+3.d0*n1+m1)+n1+m1
         rhat=rhat/pihat/(n1-pihat*(2.d0*mt+n0+3.d0*n1-2.d0*nn*pihat))
      else if(model=='Don') then
         rhat=(m2-m0)*(mt+nt)*pihat**2
     .        - (2.d0*m2*(mt+n1)-m0*n1+m2*n0)*pihat
     .        + m2*(mt+n1)
         rhat=rhat/((m2-m0)*(mt+nt)*pihat**2
     .        + ((m0-m2)*(mt+n1)+m0*nt)*pihat
     .        - m0*n1)
      else
         write(*,*) 'The model '//model//' has not been implemented!'
         stop
      end if
!      if(rhat<0.d0) print*, 'negative rhat: ', rhat, pihat
      return
      end subroutine mle0






      


      subroutine mle(pihat,rhat)
      implicit none
      include 'parameters.f'
      real(8) :: pihat(group),rhat
      real(8) :: pihat0,rhat0,dlri,dlri2
      real(8), parameter :: tiny=1.d-6
      integer :: i,count
      if(model=='Ros') rhat=1.d2
      if(model=='Don') rhat=1.d0
      count=0
      call mle0(pihat0,rhat0)
      if(rhat0<0.d0.and.model=='Ros') rhat0=1.d0 !goto 102
      if(abs(rhat0)>1.d0.and.model=='Don') rhat0=0.d0 !goto 104
!      if(rhat0<0.d0) rhat=0.d0
!      print*, 'pihat0, rhat0, initial: ', pihat0, rhat0
      if(abs(rhat0-rhat)<=tiny) rhat0=0.5d0
      do while (abs(rhat0-rhat)>tiny)
         rhat=rhat0
         call find_pi(pihat,rhat)
!         print*, 'pi & rho in iterations: ', pihat,rhat
         dlri=0.d0
         dlri2=0.d0
         do i=1,group
            if(model=='Ros') then
               dlri=dlri+m(2,i)/rhat+pihat(i)**2*m(0,i)
     .              /(rhat*pihat(i)**2-2.d0*pihat(i)+1.d0)
     .              +pihat(i)*m(1,i)/(rhat*pihat(i)-1.d0)
               dlri2=dlri2-m(2,i)/rhat**2-pihat(i)**4*m(0,i)
     .              /(rhat*pihat(i)**2-2.d0*pihat(i)+1.d0)**2
     .              -pihat(i)**2*m(1,i)/(rhat*pihat(i)-1.d0)**2
            else if(model=='Don') then
               dlri=dlri+pihat(i)*m(0,i)/(rhat*pihat(i)-pihat(i)+1.d0)
     .              +m(1,i)/(rhat-1.d0)
     .              +(1.d0-pihat(i))*m(2,i)
     .              /(rhat+pihat(i)-rhat*pihat(i))
               dlri2=dlri2
     .              -pihat(i)**2*m(0,i)/(rhat*pihat(i)-pihat(i)+1.d0)**2
     .              -m(1,i)/(rhat-1.d0)**2
     .              -(pihat(i)-1.d0)**2*m(2,i)
     .              /(rhat+pihat(i)-rhat*pihat(i))**2
            else
               write(*,*) 'The model '//model
     .              //' has not been implemented! ... routine mle'
               stop               
            end if
         end do
         rhat0=rhat-1.d0/dlri2*dlri
!         if(rhat0<0.d0.and.model=='Ros') go to 103
!         if(abs(rhat0)>1.d0.and.model=='Don') rhat0=sign(.99d0,rhat0)
         count=count+1
         if(count>2000) then
            print*, 'too many iterations: ', count
!            exit
            go to 100
         end if
      end do
      return
 100  print*, 'del_R is: ', abs(rhat0-rhat), rhat0, rhat
      print*, 'm0: ', m(0,:)
      print*, 'm1: ', m(1,:)
      print*, 'm2: ', m(2,:)
      print*, 'n0: ', n(0,:)
      print*, 'n1: ', n(1,:)
      print*, 'pihat: ', pihat(:)
      return
 102  print*, 'initial R0 is negative: ', rhat0
      return
 103  print*, 'iterative R0 is negative: ', rhat0
 104  print*, 'correlation rho0 is out of range: ', rhat0
      end subroutine mle



      

!     find pi under unconstrained MLEs
      subroutine find_pi(pix,rr)
      implicit none
      include 'parameters.f'
      real(8) :: rr,pix(group),err
      integer :: i
      real(8), parameter :: tiny=1.d-12
      real(8) :: quartic_func
      external quartic_func
      real(8) :: xaa,xbb,xcc,xdd,xee
      common/coe_quartic/xaa,xbb,xcc,xdd,xee
      complex(8) :: pii(0:2,group)
      real(8) :: pix2(2,group)
      select case(model)
      case('Ros')
         do i=1,group
            xaa=rr**2*(2.d0*(m(0,i)+m(1,i)+m(2,i))+n(0,i)+n(1,i))
            xbb=-(
     .           + (4.d0*rr+2.d0*rr**2)*m(0,i)
     .           + (5.d0*rr+2.d0*rr**2)*m(1,i)
     .           + (6.d0*rr+2.d0*rr**2)*m(2,i)
     .           + (3.d0*rr)*n(0,i)
     .           + (3.d0*rr+rr**2)*n(1,i)
     .           )
            xcc=rr*(4.d0*m(0,i)+7.d0*m(1,i)+8.d0*m(2,i)+n(0,i)
     .           +4.d0*n(1,i))
     .           +2.d0*(m(0,i)+m(1,i)+m(2,i)+n(0,i)+n(1,i))
     .           +2.d0*m(2,i)
            xdd=-(
     .           +2.d0*m(0,i)+(3.d0+2.d0*rr)*m(1,i)
     .           +(6.d0+2.d0*rr)*m(2,i)
     .           +n(0,i)+(3.d0+rr)*n(1,i)
     .           )
            xee=m(1,i)+2.d0*m(2,i)+n(1,i)
c$$$            print*, 'coefficients '
c$$$            print*, 'aa: ', xaa
c$$$            print*, 'bb: ', xbb
c$$$            print*, 'cc: ', xcc
c$$$            print*, 'dd: ', xdd
c$$$            print*, 'ee: ', xee
!            call dzero(tiny,min(1.d0/rr-tiny,1.d0-tiny),pix(i),
!     .           err,1.d-16,100000,quartic_func)
            call rt4(pix(i),xbb/xaa,xcc/xaa,xdd/xaa,xee/xaa)
c$$$            print*, 'root by rt4: ', pix(i)
         end do
      case('Don')
c$$$         if(rr>0.d0) then
c$$$         xaa=0.d0
c$$$         do i=1,group
c$$$            xbb=(rr-1.d0)**2*(2.d0*(m(0,i)+m(1,i)+m(2,i))+n(0,i)+n(1,i))
c$$$            xcc=-(rr-1.d0)*(
c$$$     .           + m(0,i)*(3.d0*rr-2.d0)+3.d0*m(1,i)*(rr-1.d0)
c$$$     .           + m(2,i)*(3.d0*rr-4.d0)+(n(0,i)+2.d0*n(1,i))*(rr-1.d0)
c$$$     .           )
c$$$            xdd=m(0,i)*rr*(rr-2.d0)+m(1,i)*(rr**2-4.d0*rr+1.d0)
c$$$     .           + m(2,i)*(rr**2-4.d0*rr+2.d0)-n(0,i)*rr
c$$$     .           + n(1,i)*(rr**2-3.d0*rr+1.d0)
c$$$            xee=rr*(m(1,i)+m(2,i)+n(1,i))
c$$$            call dzero(0.d0,1.d0,pix(i),err,1.d-16,100000
c$$$     .           ,quartic_func)
c$$$         end do
c$$$         else 
            call find_pi_rt3(pii,rr,pix2)
            pix(:)=dreal(pii(2,:))
c$$$         end if
!         print*, 'find_pi_rt3: '
!         print*, 'pii: ', pii
!         print*, 'pix2: ', pix2(1,:)
      case default
         write(*,*) "models haven't been implemented ... "//
     .        'routine find_pi'
         return
      end select
      return
      end subroutine find_pi
      



      


c$$$!     cross check the outcomes from routines pi_exp & find_pi_rt3
c$$$!     & root solved by dzero
c$$$      program xcheck_pi
c$$$      implicit none
c$$$      include 'parameters.f'
c$$$      complex(8) :: pii(group),xpii(0:2,group)
c$$$      real(8) :: pix(2,group),pix2(2,group)
c$$$      integer :: i
c$$$      rho=0.2d0
c$$$      model='Don'
c$$$      pis=pi0
c$$$      call simu_data
c$$$!      call pi_exp(pii,rho)
c$$$      call find_pi_rt3(xpii,rho,pix)
c$$$      call find_pi(pix2,rho)
c$$$      do i=1,group
c$$$!         write(*,*) 'pi_exp: ', i, pii(i)
c$$$         write(*,*) 'pi_rt3: ', i, xpii(:,i)
c$$$         write(*,*) '1 rt by dzero: ', pix(:,i)
c$$$         write(*,*) 'use quartic func in dzero: ', pix2(:,i)
c$$$      end do
c$$$      model='Ros'
c$$$      R=1.4d0
c$$$      call simu_data
c$$$      call find_pi(pix2,R)
c$$$      do i=1,group
c$$$         write(*,*) 'Rosner model (', i, ')'
c$$$         write(*,*) 'use quartic func in dzero: ', pix2(:,i)
c$$$      end do
c$$$      end program xcheck_pi




!     define a quartic func & use routine dzero to solve the roots
      real(8) function quartic_func(x)
      implicit none
      real(8) :: x
      real(8) :: xaa,xbb,xcc,xdd,xee
      common/coe_quartic/xaa,xbb,xcc,xdd,xee
      quartic_func=xaa*x**4+xbb*x**3+xcc*x**2+xdd*x+xee
      return
      end function quartic_func
      


      

!     define cubic function to be used by routine dzero
      real(8) function cubic_func(x)
      implicit none
      real(8) :: x
      real(8) :: xa,xb,xc,xd
      common/coe_cubic/xa,xb,xc,xd
      cubic_func=xa*x**3+xb*x**2+xc*x+xd
      return
      end function cubic_func



      
      
!     use direct expression for pi_i
      subroutine pi_exp(pii,rr)
      implicit none
      include 'parameters.f'
      real(8) :: rr,b(group),c(group),mt,nt,m1,m2,n1
      complex(8) :: pii(group),theta(group)
      integer :: i
      mt=0.d0
      nt=0.d0
      m1=0.d0
      m2=0.d0
      n1=0.d0
      do i=1,group
         mt=mt+m(0,i)+m(1,i)+m(2,i)
         nt=nt+n(0,i)+n(1,i)
         m1=m1+m(1,i)
         m2=m2+m(2,i)
         n1=n1+n(1,i)
      end do
      do i=1,group
         b(i)=rr*(3.d0*mt+nt+n1)
     .        -(2.d0*m(0,i)+3.d0*m(1,i)+4.d0*m(2,i)+n(0,i)+2.d0*n(1,i))
         c(i)=rr**2*(mt+n1)
     .        -rr*(2.d0*m(0,i)+4.d0*m(1,i)+4.d0*m(2,i)+n(0,i)
     .        +3.d0*n(1,i))+m1+2.d0*m2+n1
         theta(i)=2.d0*b(i)**3-9.d0*(2.d0*mt+nt)*b(i)*c(i)
     .        +27.d0*(2.d0*mt+nt)**2*(m1+m2+n1)*(1.d0-rr)*rr
         theta(i)=theta(i)/2.d0
     .        /cmplx(b(i)**2-3.d0*c(i)*(2.d0*mt+nt))**1.5d0
         theta(i)=acos(theta(i))/3.d0
         pii(i)=b(i)
     .        -sqrt(cmplx(b(i)**2-3.d0*c(i)*(2.d0*mt+nt)))*
     .        (cos(theta(i))-dsqrt(3.d0)*sin(theta(i)))
         pii(i)=pii(i)/(2.d0*mt+nt)/(rr-1.d0)/3.d0
c$$$         write(*,*) '******** debugging ********'
c$$$         write(*,*) 'index i: ', i
c$$$         write(*,*) 'b = ', b(i)
c$$$         write(*,*) 'c = ', c(i)
c$$$         write(*,*) 'theta=', theta(i)
c$$$         write(*,*) 'pii = ', pii(i)
      end do
      return
      end subroutine pi_exp
      


      

!     use routine rt3 to find pi_i
      subroutine find_pi_rt3(pii,rr,pix)
      implicit none
      include 'parameters.f'
      real(8) :: rr,aa,bb,cc,dd,pix(2,group)
      complex(8) :: pii(0:2,group),xx(0:2)
      integer :: i
      real(8) :: cubic_func
      external cubic_func
      real(8) :: xa,xb,xc,xd
      common/coe_cubic/xa,xb,xc,xd
      do i=1,group
         aa=(rr-1.d0)**2*(2.d0*(m(0,i)+m(1,i)+m(2,i))+n(0,i)+n(1,i))
         bb=-(rr-1.d0)*(
     .        + m(0,i)*(3.d0*rr-2.d0)+3.d0*m(1,i)*(rr-1.d0)
     .        + m(2,i)*(3.d0*rr-4.d0)+(n(0,i)+2.d0*n(1,i))*(rr-1.d0)
     .        )
         cc=m(0,i)*rr*(rr-2.d0)+m(1,i)*(rr**2-4.d0*rr+1.d0)
     .        + m(2,i)*(rr**2-4.d0*rr+2.d0)-n(0,i)*rr
     .        + n(1,i)*(rr**2-3.d0*rr+1.d0)
         dd=rr*(m(1,i)+m(2,i)+n(1,i))
         call rt3(xx,aa,bb,cc,dd)
         pii(:,i)=xx(:)
         xa=aa
         xb=bb
         xc=cc
         xd=dd
!         call dzero(0.d0,1.d0,pix(1,i),pix(2,i),1.d-16,100000
!     .        ,cubic_func)
      end do
      return
      end subroutine find_pi_rt3


      

!     roots of cubic eqation: aa*x^3+bb*x^2+cc*x+dd=0
!     there are 3 roots which could be complex
!     we later on will focus on one of the (real) roots in [0,1]
      subroutine rt3(x,aa,bb,cc,dd)
      implicit none
      real(8) :: aa,bb,cc,dd,d0,d1
      complex(8) :: x(0:2),r,ir,ksi
      integer :: i
      d0=bb**2-3.d0*aa*cc
      d1=2.d0*bb**3-9.0*aa*bb*cc+27.d0*aa**2*dd
      r=0.5d0*(d1+sqrt(cmplx(d1**2-4.d0*d0**3,0.d0,8)))
      r=r**(1.d0/3.d0)
      ir=0.5d0*(d1-sqrt(cmplx(d1**2-4.d0*d0**3,0.d0,8)))
      ir=ir**(1.d0/3.d0)
      ksi=cmplx(-1.d0,sqrt(3.d0),8)/2.d0
      do i=0,2
         x(i)=bb+ksi**i*r+ir/ksi**i
      end do
      x=-x/3.d0/aa
      return
      end subroutine rt3

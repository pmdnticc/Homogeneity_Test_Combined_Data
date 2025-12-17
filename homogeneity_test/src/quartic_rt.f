c$$$!     compare the root by analytic and numeric solution, given rt in [0,1]
c$$$      program rt_cmp
c$$$      implicit none
c$$$      real(8) :: xl,xh,x,err,quartic_func
c$$$      integer :: ncall
c$$$      external quartic_func
c$$$      real(8) :: b,c,d,e
c$$$      common/coe/b,c,d,e
c$$$      b=-2.d0
c$$$      c=1.d0
c$$$      d=-1.d0
c$$$      e=1.d0
c$$$      xl=0.d0
c$$$      xh=1.d0
c$$$      ncall=100000
c$$$      call dzero(xl,xh,x,err,1.d-8,ncall,quartic_func)
c$$$      write(*,*) '****** use dzero to find the root ******'
c$$$      write(*,*) 'root is: ', x
c$$$      write(*,*) 'error is: ', err
c$$$      call rt4(x,b,c,d,e)
c$$$      write(*,*) '****** use rt4 to find the root ******'
c$$$      write(*,*) 'root is: ', x
c$$$      return
c$$$      end program rt_cmp





      
!     analytic solution of one real root of quartic equation
!     x^4 + b x^3 + c x^2 + d x + e = 0
!     such that x in [0,1]
      subroutine rt4(x,b,c,d,e)
      implicit none
      real(8) :: x,b,c,d,e,d0,d1,p,q,r,s,theta,s1
      d0=c**2-3.d0*b*d+12.d0*e
      d1=2.d0*c**3-9.d0*b*c*d+27.d0*b**2*e+27.d0*d**2-72.d0*c*e
      p=(c-3.d0*b**2/8.d0)
      q=(b**3-4.d0*b*c+8.d0*d)/8.d0
      r=((d1+dsqrt(d1**2-4*d0**3))/2.d0)**(1.d0/3.d0)
      s=dsqrt(-2.d0*p/3.d0+(r+d0/r)/3.d0)/2.d0
      if(d1**2-4.d0*d0**3<0.d0) then
         theta=dacos(d1/dsqrt(4.d0*d0**3))/3.d0
         s1=dsqrt(-2.d0*p/3.d0+2.d0*dsqrt(d0)*dcos(theta)/3.d0)/2.d0
         x=-b/4.d0-s1-dsqrt(-4.d0*s1**2-2.d0*p+q/s1)/2.d0
      elseif(q>0.d0) then
         x=-b/4.d0-s-dsqrt(-4.d0*s**2-2.d0*p+q/s)/2.d0
      else
         x=-b/4.d0+s-dsqrt(-4.d0*s**2-2.d0*p-q/s)/2.d0
      end if
      return
      end subroutine rt4



c$$$!     define a quartic func & use routine dzero to solve the roots
c$$$      real(8) function quartic_func(x)
c$$$      implicit none
c$$$      real(8) :: x
c$$$      real(8) :: b,c,d,e
c$$$      common/coe/b,c,d,e
c$$$      quartic_func=x**4+b*x**3+c*x**2+d*x+e
c$$$      return
c$$$      end function quartic_func

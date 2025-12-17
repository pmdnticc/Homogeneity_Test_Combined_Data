!     returned fisher's information matrix
      subroutine fisher_info(fi,ph,rh)
      implicit none
      include 'parameters.f'
      include 'info.f'
      real(8) :: fi(group+1,group+1),ph(group),rh,mi,ni
      real(8) :: den1,den2
      real(8), parameter :: tiny=1.d-18,big=1.d16
      integer :: i,j
      logical :: debug=.false.!.true.
      fi=0.d0
      select case(model)
      case('Ros')
         do i=1,group
            mi=m(0,i)+m(1,i)+m(2,i)
            ni=n(0,i)+n(1,i)
            den1=1.d0-rh*ph(i)
            den2=rh*ph(i)**2-2.d0*ph(i)+1.d0
            if(den1==0.d0) den1=tiny
            if(den2==0.d0) den2=tiny
            fi(i,i)=2.d0*mi*(2.d0*rh**2*ph(i)**2-rh*ph(i)**2
     .           -2.d0*rh*ph(i)+1.d0)/ph(i)/den1/den2
     .           +ni/ph(i)/(1.d0-ph(i))
            if(fi(i,i)<0.d0.or.fi(i,i).ne.fi(i,i).or.
     .           fi(i,i)>big) info=-1
            fi(i,group+1)=-2.d0*(1.d0-rh)*ph(i)**2*mi
     .           /den1/den2
            if(fi(i,group+1).ne.fi(i,group+1).or.
     .           fi(i,group+1)>big)
     .           info=-1
            fi(group+1,i)=fi(i,group+1)
            fi(group+1,group+1)=fi(group+1,group+1)
     .           +ph(i)**2*mi*(rh*ph(i)-2.d0*ph(i)+1.d0)/rh
     .           /den1/den2
            if(fi(group+1,group+1)<0.d0.or.
     .           fi(group+1,group+1).ne.fi(group+1,group+1).or.
     .           fi(group+1,group+1)>big)
     .           info=-1
         end do
      case('Don')
         do i=1,group
            mi=m(0,i)+m(1,i)+m(2,i)
            ni=n(0,i)+n(1,i)
c$$$            if(debug) then ! fi(i,group+1), fi(group+1,group+1) differ from the correct ones
c$$$               fi(i,i)=
c$$$     .              +mi*(2.d0*ph(i)*(1.d0-ph(i))*(rh-1.d0)*(2.d0*rh
c$$$     .              -1.d0)+rh*(2.d0-rh))/ph(i)/(1.d0-ph(i))
c$$$     .              /(rh+ph(i)-rh*ph(i))/(rh*ph(i)-ph(i)+1.d0)
c$$$     .              +ni/ph(i)/(1.d0-ph(i))
c$$$               fi(i,group+1)=
c$$$     .              -mi*rh*(2.d0*ph(i)-1.d0)/(rh+ph(i)-rh*ph(i))
c$$$     .              /(rh*ph(i)-ph(i)+1.d0)
c$$$               fi(group+1,i)=fi(i,group+1)
c$$$               fi(group+1,group+1)=fi(group+1,group+1)
c$$$     .              -mi*ph(i)*(rh+1.d0)*(1.d0-ph(i))
c$$$     .              /(rh+ph(i)-rh*ph(i))/(rh*ph(i)-ph(i)+1.d0)
c$$$               go to 111
c$$$            end if
c$$$            fi(i,i)=mi*(2.d0*ph(i)*(1.d0-ph(i))*(rh-1.d0)*(2.d0*rh-1.d0)
c$$$     .           +rh*(2.d0-rh))/ph(i)/(1.d0-ph(i))/(rh+ph(i)-rh*ph(i))
c$$$     .           /(rh*ph(i)-ph(i)+1.d0)
c$$$     .           +ni/ph(i)/(1.d0-ph(i))
c$$$            fi(i,group+1)=mi*rh*(2.d0*ph(i)-1.d0)/(rh+ph(i)-rh*ph(i))
c$$$     .           /(rh*ph(i)-ph(i)+1.d0)
c$$$            fi(group+1,i)=fi(i,group+1)
c$$$            fi(group+1,group+1)=fi(group+1,group+1)
c$$$     .           +mi*ph(i)*(rh+1.d0)*(1.d0-ph(i))/(1.d0-rh)
c$$$     .           /(rh+ph(i)-rh*ph(i))/(rh*ph(i)-ph(i)+1.d0)
            den1=rh+ph(i)-rh*ph(i)
            den2=rh*ph(i)-ph(i)+1.d0
            if(den1==0.d0) den1=tiny
            if(den2==0.d0) den2=tiny
            fi(i,i)=mi*(2.d0*ph(i)*(1.d0-ph(i))*(rh-1.d0)*(2.d0*rh-1.d0)
     .           +rh*(2.d0-rh))/ph(i)/(1.d0-ph(i))/den1
     .           /den2
     .           +ni/ph(i)/(1.d0-ph(i))
            if(fi(i,i)<0.d0.or.fi(i,i).ne.fi(i,i).or.
     .           fi(i,i)>big) info=-1
            fi(i,group+1)=mi*rh*(2.d0*ph(i)-1.d0)/den1
     .           /den2
            if(fi(i,group+1).ne.fi(i,group+1).or.
     .           fi(i,group+1)>big)
     .           info=-1
            fi(group+1,i)=fi(i,group+1)
            fi(group+1,group+1)=fi(group+1,group+1)
     .           +mi*ph(i)*(rh+1.d0)*(1.d0-ph(i))/(1.d0-rh)
     .           /den1/den2
            if(fi(group+1,group+1)<0.d0.or.
     .           fi(group+1,group+1).ne.fi(group+1,group+1).or.
     .           fi(group+1,group+1)>big)
     .           info=-1
c$$$ 111        continue
         end do
      case default
         stop "models haven't bee implemented ... fisher_info"
      end select
c$$$      if(info/=0) then
c$$$         print*, "Fisher Info issue"
c$$$         print*, "ph: ", ph
c$$$         print*, "rh: ", rh
c$$$      end if
      return
      end subroutine fisher_info




c$$$!     test of fisher_info
c$$$      program test_fisher
c$$$      implicit none
c$$$      include 'parameters.f'
c$$$      real(8) :: rr,fi(group+1,group+1),inv_fi(group+1,group+1)
c$$$c$$$      real(8) :: inv
c$$$c$$$      external inv
c$$$      Interface
c$$$         function inv(A) result(Ainv)
c$$$         real(8), dimension(:,:), intent(in) :: A
c$$$         real(8), dimension(size(A,1),size(A,2)) :: Ainv
c$$$         end function inv
c$$$      end interface
c$$$      model='Ros'
c$$$      R=1.4d0
c$$$      rho=0.4d0
c$$$      call simu_data
c$$$      call mle(pis,rr)
c$$$      call fisher_info(fi,pis,rr)
c$$$      write(*,*) "Fisher's information is: "
c$$$      write(*,*) fi
c$$$      inv_fi=inv(fi)
c$$$      write(*,*) "inverse of FI: "
c$$$      write(*,*) inv_fi
c$$$      return
c$$$      end program test_fisher

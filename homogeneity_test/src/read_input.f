      subroutine read_input
      implicit none
      include 'parameters.f'
      character(len=100) :: buffer,label
      integer :: posi,posf
      integer, parameter :: fh=15
      integer :: ios=0
      integer :: line=0
!      real(8) :: mmt(8),nnt(8)
      integer :: i
  
      open(fh,file="./input.DAT")

      do while(ios==0)
         read(fh,"(A)",iostat=ios) buffer
         if(ios==0) then
            line=line+1

            posi=scan(buffer,"[")
            posf=scan(buffer,"]")
            label=buffer(posi+1:posf-1)
            buffer=buffer(1:posi-1)

            select case(label)
               case("model")
                  read(buffer,*,iostat=ios) model
                  write(6,*) "The model is: ", model
               case("H0")
                  read(buffer,*,iostat=ios) H0
                  write(6,*) "under H0? ", H0
               case("pi0 under H0")
                  read(buffer,*,iostat=ios) pi0
                  if(H0) write(6,*) "pi0 is: ", pi0
               case("pi11")
                  read(buffer,*,iostat=ios) pi1(1)
                  if(.not.H0) write(6,*) "pi11 is: ", pi1(1)
               case("pi12")
                  read(buffer,*,iostat=ios) pi1(2)
                  if(.not.H0) write(6,*) "pi12 is: ", pi1(2)
               case("pi13")
                  read(buffer,*,iostat=ios) pi1(3)
                  if((.not.H0).and.(group>=3))
     .                 write(6,*) "pi13 is: ", pi1(3)
               case("pi14")
                  read(buffer,*,iostat=ios) pi1(4)
                  if((.not.H0).and.(group>=4))
     .                 write(6,*) "pi14 is: ", pi1(4)
               case("pi15")
                  read(buffer,*,iostat=ios) pi1(5)
                  if((.not.H0).and.(group>=5))
     .                 write(6,*) "pi15 is: ", pi1(5)
               case("pi16")
                  read(buffer,*,iostat=ios) pi1(6)
                  if((.not.H0).and.(group>=6))
     .                 write(6,*) "pi16 is: ", pi1(6)
               case("pi17")
                  read(buffer,*,iostat=ios) pi1(7)
                  if((.not.H0).and.(group>=7))
     .                 write(6,*) "pi17 is: ", pi1(7)
               case("pi18")
                  read(buffer,*,iostat=ios) pi1(8)
                  if((.not.H0).and.(group>=8))
     .                 write(6,*) "pi18 is: ", pi1(8)
               case("R0")
                  read(buffer,*,iostat=ios) R
                  write(6,*) "R0: ", R
               case("rho0")
                  read(buffer,*,iostat=ios) rho
                  write(6,*) "rho0: ", rho
               case("m1")
                  read(buffer,*,iostat=ios) mmt(1)
                  write(6,*) "m1: ", mmt(1)
               case("m2")
                  read(buffer,*,iostat=ios) mmt(2)
                  write(6,*) "m2: ", mmt(2)
               case("m3")
                  read(buffer,*,iostat=ios) mmt(3)
                  if(group>=3) write(6,*) "m3: ", mmt(3)
               case("m4")
                  read(buffer,*,iostat=ios) mmt(4)
                  if(group>=4) write(6,*) "m4: ", mmt(4)
               case("m5")
                  read(buffer,*,iostat=ios) mmt(5)
                  if(group>=5) write(6,*) "m5: ", mmt(5)
               case("m6")
                  read(buffer,*,iostat=ios) mmt(6)
                  if(group>=6) write(6,*) "m6: ", mmt(6)
               case("m7")
                  read(buffer,*,iostat=ios) mmt(7)
                  if(group>=7) write(6,*) "m7: ", mmt(7)
               case("m8")
                  read(buffer,*,iostat=ios) mmt(8)
                  if(group>=8) write(6,*) "m8: ", mmt(8)
               case("n1")
                  read(buffer,*,iostat=ios) nnt(1)
                  write(6,*) "n1: ", nnt(1)
               case("n2")
                  read(buffer,*,iostat=ios) nnt(2)
                  write(6,*) "n2: ", nnt(2)
               case("n3")
                  read(buffer,*,iostat=ios) nnt(3)
                  if(group>=3) write(6,*) "n3: ", nnt(3)
               case("n4")
                  read(buffer,*,iostat=ios) nnt(4)
                  if(group>=4) write(6,*) "n4: ", nnt(4)
               case("n5")
                  read(buffer,*,iostat=ios) nnt(5)
                  if(group>=5) write(6,*) "n5: ", nnt(5)
               case("n6")
                  read(buffer,*,iostat=ios) nnt(6)
                  if(group>=6) write(6,*) "n6: ", nnt(6)
               case("n7")
                  read(buffer,*,iostat=ios) nnt(7)
                  if(group>=7) write(6,*) "n7: ", nnt(7)
               case("n8")
                  read(buffer,*,iostat=ios) nnt(8)
                  if(group>=8) write(6,*) "n8: ", nnt(8)
               case("seed")
                  read(buffer,*,iostat=ios) seed
                  write(6,*) "file seed: ", seed
!               case default
!                  write(6,*) "skipping other parameters reading"
            end select

         end if

      end do
      
      model=trim(model)
!     # of obs for bilateral (mtot) & unilateral (ntot)
!     tot # of observations: tot=mtot+ntot      
      mtot=0
      ntot=0
      do i=1,group
         m_vec(i)=mmt(i)
         n_vec(i)=nnt(i)
         mtot=mtot+m_vec(i)
         ntot=ntot+n_vec(i)
      end do
      tot=mtot+ntot
      
      close(fh)
  
      end subroutine read_input



c$$$      program test_read
c$$$      implicit none
c$$$      include 'parameters.f'
c$$$      call read_input
c$$$      end program test_read

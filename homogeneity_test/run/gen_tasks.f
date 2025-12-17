      program gen_tasks
      implicit none
      character(len=3) :: model,fig,run,cmp,iseed
      character(len=2) :: hyp
      character(len=1) :: tab=char(9),sng
      real(8) :: rho(3),RR(3)
      real(8) :: ppi0(2),ppi1_g2(2,2),ppi1_g3(2,3),ppi1_g4(2,4),
     .     ppi1_g5(2,5),ppi1_g6(2,6),ppi1_g7(2,7),ppi1_g8(2,8)
      integer :: i,i0,j,k,l,gen_seed,ng
      integer :: m_g2(3,2),m_g3(3,3),m_g4(3,4),m_g5(3,5),m_g6(3,6),
     .     m_g7(3,7),m_g8(3,8),n_g2(3,2),n_g3(3,3),n_g4(3,4),
     .     n_g5(3,5),n_g6(3,6),n_g7(3,7),n_g8(3,8)
!     m_gx(3,x)/n_gx(3,x) -- configuration for m/n when group=x;
!     first entry has 3 dimensions: 1&2 for balanced
!     and 3 for unbalanced sample sizes

      write(*,*) "generate configurations?"
      read(*,*) fig
      if(fig(1:1)=='y'.or.fig(1:1)=='Y') then
         write(*,*) "What is the model?"
 102     write(*,*) "type 'Don' for Donner's model; "
     .        //"'Ros' for Rosner's model"
         read(*,*) model
         if(model/='Don'.and.model/='Ros') go to 102
         write(*,*) "empirical type I error rates or powers?"
 103     write(*,*) "type 'H0' for TIE rates; 'H1' for powers"
         read(*,*) hyp
         if(hyp/='H0'.and.hyp/='H1') go to 103
      else
         write(*,*) "...no generating configurations..."
         write(*,*) "run excutives?"
         read(*,*) run
         if(run(1:1)=='y'.or.run(1:1)=='Y') then
            go to 101
         else
            write(*,*) "nothing is done!!!"
            go to 100
         end if
      end if

!     pi_0 under H0
      ppi0=[0.2d0,0.5d0]![0.3d0,0.5d0]
!     pi_1 under H1A & H1B
      ppi1_g2(1,:)=[0.25d0,0.4d0]
      ppi1_g2(2,:)=[0.2d0,0.4d0]
      ppi1_g3(1,:)=[0.25d0,0.3d0,0.4d0]
      ppi1_g3(2,:)=[0.2d0,0.3d0,0.4d0]
      ppi1_g4(1,:)=[0.25d0,0.3d0,0.35d0,0.4d0]
      ppi1_g4(2,:)=[0.2d0,0.2d0,0.4d0,0.4d0]
      ppi1_g5(1,:)=[0.25d0,0.3d0,0.3d0,0.35d0,0.4d0]
      ppi1_g5(2,:)=[0.2d0,0.2d0,0.3d0,0.4d0,0.4d0]
      ppi1_g6(1,:)=[0.25d0,0.3d0,0.3d0,0.35d0,0.35d0,0.4d0]
      ppi1_g6(2,:)=[0.2d0,0.2d0,0.3d0,0.3d0,0.4d0,0.4d0]
      ppi1_g7(1,:)=[0.25d0,0.25d0,0.3d0,0.3d0,0.35d0,0.35d0,0.4d0]
      ppi1_g7(2,:)=[0.2d0,0.2d0,0.3d0,0.3d0,0.3d0,0.4d0,0.4d0]
      ppi1_g8(1,:)=[0.25d0,0.3d0,0.35d0,0.4d0,0.25d0,0.3d0,0.35d0,0.4d0]
      ppi1_g8(2,:)=[0.2d0,0.2d0,0.4d0,0.4d0,0.2d0,0.2d0,0.4d0,0.4d0]
!      if(model=='Ros') then
         RR=[1.2d0,1.5d0,1.8d0]!rho=[1.2d0,1.5d0,2.d0]
!      else if(model=='Don') then
         rho=[0.1d0,0.5d0,0.7d0]![0.4d0,0.5d0,0.7d0]
!      else
!         write(*,*) "model not recognized!!!"
!         go to 102
!      end if
!     configurations of m and n
      m_g2(1,:)=20
      m_g2(2,:)=40
      n_g2(1,:)=20
      n_g2(2,:)=40
      
      m_g3(1,:)=20
      m_g3(2,:)=40
      n_g3(1,:)=20
      n_g3(2,:)=40
      
      m_g4(1,:)=20
      m_g4(2,:)=40
      n_g4(1,:)=20
      n_g4(2,:)=40
      
      m_g5(1,:)=20
      m_g5(2,:)=40
      n_g5(1,:)=20
      n_g5(2,:)=40
      
      m_g6(1,:)=20
      m_g6(2,:)=40
      n_g6(1,:)=20
      n_g6(2,:)=40
      
      m_g7(1,:)=20
      m_g7(2,:)=40
      n_g7(1,:)=20
      n_g7(2,:)=40
      
      m_g8(1,:)=20
      m_g8(2,:)=40
      n_g8(1,:)=20
      n_g8(2,:)=40
      
      m_g2(3,1)=20
      m_g2(3,2)=40
      n_g2(3,1)=20
      n_g2(3,2)=40
      
      m_g3(3,1)=20
      m_g3(3,2)=30
      m_g3(3,3)=40
      n_g3(3,1)=20
      n_g3(3,2)=30
      n_g3(3,3)=40

      m_g4(3,1)=20
      m_g4(3,2)=20
      m_g4(3,3)=40
      m_g4(3,4)=40
      n_g4(3,1)=20
      n_g4(3,2)=20
      n_g4(3,3)=40
      n_g4(3,4)=40

      m_g5(3,1)=20
      m_g5(3,2)=20
      m_g5(3,3)=30
      m_g5(3,4)=40
      m_g5(3,5)=40
      n_g5(3,1)=20
      n_g5(3,2)=20
      n_g5(3,3)=30
      n_g5(3,4)=40
      n_g5(3,5)=40
      
      m_g6(3,1)=20
      m_g6(3,2)=20
      m_g6(3,3)=30
      m_g6(3,4)=30
      m_g6(3,5)=40
      m_g6(3,6)=40
      n_g6(3,1)=20
      n_g6(3,2)=20
      n_g6(3,3)=30
      n_g6(3,4)=30
      n_g6(3,5)=40
      n_g6(3,6)=40

      m_g7(3,1)=20
      m_g7(3,2)=20
      m_g7(3,3)=30
      m_g7(3,4)=30
      m_g7(3,5)=30
      m_g7(3,6)=40
      m_g7(3,7)=40
      n_g7(3,1)=20
      n_g7(3,2)=20
      n_g7(3,3)=30
      n_g7(3,4)=30
      n_g7(3,5)=30
      n_g7(3,6)=40
      n_g7(3,7)=40

      m_g8(3,1)=20
      m_g8(3,2)=20
      m_g8(3,3)=30
      m_g8(3,4)=30
      m_g8(3,5)=40
      m_g8(3,6)=40
      m_g8(3,7)=50
      m_g8(3,8)=50
      n_g8(3,1)=20
      n_g8(3,2)=20
      n_g8(3,3)=30
      n_g8(3,4)=30
      n_g8(3,5)=40
      n_g8(3,6)=40
      n_g8(3,7)=50
      n_g8(3,8)=50

!     
      
      gen_seed=100

      write(*,*) "Compile excutables w/ corresponding # of groups???"
      read(*,*) cmp
      if(cmp(1:1)=='y'.or.cmp(1:1)=='Y') then
      call system('mkdir -pv ./Bin') ! assume current folder is "run"
      call system('for grp in {2..8}; do sed s/group=.*/group=$grp/g '
     .     //'< ../Inc/const')
      call system('cd ..; '//
     .     'for i in {2..8}; do sed s/group=.*/group=$i/g '//
     .     '< Inc/constants.f-save > Inc/constants.f; make clean; '//
     .     'make; cp Bin/main.exe run/Bin/main_g$i.exe; done')
      else
         write(*,*) "No excutables are (re)compiled this time!!!"
      end if
      write(*,*) "*************************************************"
      write(*,*) "Find excutables in folder 'Bin' under current DIR"
      write(*,*) "Symbolic links are generated when creating config"
      write(*,*) "*************************************************"
      write(*,*) "Now writing configurations with group = 2, ..., 8"
!      do ng=2,8 ! # of groups
      do i0=1,3
         ng=2**i0
         write(sng,'(I1)') ng
         do i=1,size(ppi0)
            do j=1,size(rho)
!     convert rho0 notation to R0 in Rosner's model
               if(model=='Ros'.and.hyp=='H0')
     .              RR(j)=(1-ppi0(i))*rho(j)/ppi0(i)+1.d0
               do k=1,2 ! for H1A & H1B
                  do l=1,3 ! even & uneven sample size
                     write(iseed,'(I3)') gen_seed
                     call system('mkdir -pv ./run_'//iseed//'_'//model
     .                    //'_'//hyp)
!     create symbolic links of the excutives "main_g$i.exe" in each run_seed folder
                     call system('ln -s ../Bin/main_g'//sng//'.exe '
     .                    //'./run_'//iseed//'_'//model
     .                    //'_'//hyp//'/main.exe')
                     open(unit=90,file='./run_'//iseed//'_'//model
     .                    //'_'//hyp//'/input.DAT',status='unknown')
                     write(90,'(4A)') "'"//model//"'",tab,tab,"[model]"
                     if(hyp=='H0') then
                        write(90,'(4A)') ".true.",tab,tab,"[H0]"
                     else
                        write(90,'(4A)') ".false.",tab,tab,"[H0]"
                     end if
                     write(90,'(F4.2,3A)') ppi0(i),tab,tab,
     .                    "[pi0 under H0]"
                     write(90,'(A)') "[pi1 under H1; "//
     .                    "values may differ by groups (upto 8)]"
                     select case(ng)
                       case(2)
                          write(90,'(F4.2,3A)') ppi1_g2(k,1),tab,tab,
     .                         "[pi11]"
                          write(90,'(F4.2,3A)') ppi1_g2(k,2),tab,tab,
     .                         "[pi12]"
                       case(3)
                          write(90,'(F4.2,3A)') ppi1_g3(k,1),tab,tab,
     .                         "[pi11]"
                          write(90,'(F4.2,3A)') ppi1_g3(k,2),tab,tab,
     .                         "[pi12]"
                          write(90,'(F4.2,3A)') ppi1_g3(k,3),tab,tab,
     .                         "[pi13]"
                       case(4)
                          write(90,'(F4.2,3A)') ppi1_g4(k,1),tab,tab,
     .                         "[pi11]"
                          write(90,'(F4.2,3A)') ppi1_g4(k,2),tab,tab,
     .                         "[pi12]"
                          write(90,'(F4.2,3A)') ppi1_g4(k,3),tab,tab,
     .                         "[pi13]"
                          write(90,'(F4.2,3A)') ppi1_g4(k,4),tab,tab,
     .                         "[pi14]"
                       case(5)
                          write(90,'(F4.2,3A)') ppi1_g5(k,1),tab,tab,
     .                         "[pi11]"
                          write(90,'(F4.2,3A)') ppi1_g5(k,2),tab,tab,
     .                         "[pi12]"
                          write(90,'(F4.2,3A)') ppi1_g5(k,3),tab,tab,
     .                         "[pi13]"
                          write(90,'(F4.2,3A)') ppi1_g5(k,4),tab,tab,
     .                         "[pi14]"
                          write(90,'(F4.2,3A)') ppi1_g5(k,5),tab,tab,
     .                         "[pi15]"
                       case(6)
                          write(90,'(F4.2,3A)') ppi1_g6(k,1),tab,tab,
     .                         "[pi11]"
                          write(90,'(F4.2,3A)') ppi1_g6(k,2),tab,tab,
     .                         "[pi12]"
                          write(90,'(F4.2,3A)') ppi1_g6(k,3),tab,tab,
     .                         "[pi13]"
                          write(90,'(F4.2,3A)') ppi1_g6(k,4),tab,tab,
     .                         "[pi14]"
                          write(90,'(F4.2,3A)') ppi1_g6(k,5),tab,tab,
     .                         "[pi15]"
                          write(90,'(F4.2,3A)') ppi1_g6(k,6),tab,tab,
     .                         "[pi16]"
                       case(7)
                          write(90,'(F4.2,3A)') ppi1_g7(k,1),tab,tab,
     .                         "[pi11]"
                          write(90,'(F4.2,3A)') ppi1_g7(k,2),tab,tab,
     .                         "[pi12]"
                          write(90,'(F4.2,3A)') ppi1_g7(k,3),tab,tab,
     .                         "[pi13]"
                          write(90,'(F4.2,3A)') ppi1_g7(k,4),tab,tab,
     .                         "[pi14]"
                          write(90,'(F4.2,3A)') ppi1_g7(k,5),tab,tab,
     .                         "[pi15]"
                          write(90,'(F4.2,3A)') ppi1_g7(k,6),tab,tab,
     .                         "[pi16]"
                          write(90,'(F4.2,3A)') ppi1_g7(k,7),tab,tab,
     .                         "[pi17]"
                       case(8)
                          write(90,'(F4.2,3A)') ppi1_g8(k,1),tab,tab,
     .                         "[pi11]"
                          write(90,'(F4.2,3A)') ppi1_g8(k,2),tab,tab,
     .                         "[pi12]"
                          write(90,'(F4.2,3A)') ppi1_g8(k,3),tab,tab,
     .                         "[pi13]"
                          write(90,'(F4.2,3A)') ppi1_g8(k,4),tab,tab,
     .                         "[pi14]"
                          write(90,'(F4.2,3A)') ppi1_g8(k,5),tab,tab,
     .                         "[pi15]"
                          write(90,'(F4.2,3A)') ppi1_g8(k,6),tab,tab,
     .                         "[pi16]"
                          write(90,'(F4.2,3A)') ppi1_g8(k,7),tab,tab,
     .                         "[pi17]"
                          write(90,'(F4.2,3A)') ppi1_g8(k,8),tab,tab,
     .                         "[pi18]"
                     end select
                     write(90,'(F4.2,3A)') RR(j),tab,tab,"[R0]"
                     write(90,'(F4.2,3A)') rho(j),tab,tab,"[rho0]"
                     write(90,*)
                     write(90,*)
                     write(90,'(A)') "[input parameters containing "
     .                    //"even and uneven cases]"
                     write(90,'(A)') "[individual cell in mi & ni; "//
     .                    "# of groups is upto 8]"
                     select case(ng)
                       case(2)
                          write(90,'(I2,3A)') m_g2(l,1),tab,tab,"[m1]"
                          write(90,'(I2,3A)') m_g2(l,2),tab,tab,"[m2]"
                          write(90,'(I2,3A)') n_g2(l,1),tab,tab,"[n1]"
                          write(90,'(I2,3A)') n_g2(l,2),tab,tab,"[n2]"
                       case(3)
                          write(90,'(I2,3A)') m_g3(l,1),tab,tab,"[m1]"
                          write(90,'(I2,3A)') m_g3(l,2),tab,tab,"[m2]"
                          write(90,'(I2,3A)') m_g3(l,3),tab,tab,"[m3]"
                          write(90,'(I2,3A)') n_g3(l,1),tab,tab,"[n1]"
                          write(90,'(I2,3A)') n_g3(l,2),tab,tab,"[n2]"
                          write(90,'(I2,3A)') n_g3(l,3),tab,tab,"[n3]"
                       case(4)
                          write(90,'(I2,3A)') m_g4(l,1),tab,tab,"[m1]"
                          write(90,'(I2,3A)') m_g4(l,2),tab,tab,"[m2]"
                          write(90,'(I2,3A)') m_g4(l,3),tab,tab,"[m3]"
                          write(90,'(I2,3A)') m_g4(l,4),tab,tab,"[m4]"
                          write(90,'(I2,3A)') n_g4(l,1),tab,tab,"[n1]"
                          write(90,'(I2,3A)') n_g4(l,2),tab,tab,"[n2]"
                          write(90,'(I2,3A)') n_g4(l,3),tab,tab,"[n3]"
                          write(90,'(I2,3A)') n_g4(l,4),tab,tab,"[n4]"
                       case(5)
                          write(90,'(I2,3A)') m_g5(l,1),tab,tab,"[m1]"
                          write(90,'(I2,3A)') m_g5(l,2),tab,tab,"[m2]"
                          write(90,'(I2,3A)') m_g5(l,3),tab,tab,"[m3]"
                          write(90,'(I2,3A)') m_g5(l,4),tab,tab,"[m4]"
                          write(90,'(I2,3A)') m_g5(l,5),tab,tab,"[m5]"
                          write(90,'(I2,3A)') n_g5(l,1),tab,tab,"[n1]"
                          write(90,'(I2,3A)') n_g5(l,2),tab,tab,"[n2]"
                          write(90,'(I2,3A)') n_g5(l,3),tab,tab,"[n3]"
                          write(90,'(I2,3A)') n_g5(l,4),tab,tab,"[n4]"
                          write(90,'(I2,3A)') n_g5(l,5),tab,tab,"[n5]"
                       case(6)
                          write(90,'(I2,3A)') m_g6(l,1),tab,tab,"[m1]"
                          write(90,'(I2,3A)') m_g6(l,2),tab,tab,"[m2]"
                          write(90,'(I2,3A)') m_g6(l,3),tab,tab,"[m3]"
                          write(90,'(I2,3A)') m_g6(l,4),tab,tab,"[m4]"
                          write(90,'(I2,3A)') m_g6(l,5),tab,tab,"[m5]"
                          write(90,'(I2,3A)') m_g6(l,6),tab,tab,"[m6]"
                          write(90,'(I2,3A)') n_g6(l,1),tab,tab,"[n1]"
                          write(90,'(I2,3A)') n_g6(l,2),tab,tab,"[n2]"
                          write(90,'(I2,3A)') n_g6(l,3),tab,tab,"[n3]"
                          write(90,'(I2,3A)') n_g6(l,4),tab,tab,"[n4]"
                          write(90,'(I2,3A)') n_g6(l,5),tab,tab,"[n5]"
                          write(90,'(I2,3A)') n_g6(l,6),tab,tab,"[n6]"
                       case(7)
                          write(90,'(I2,3A)') m_g7(l,1),tab,tab,"[m1]"
                          write(90,'(I2,3A)') m_g7(l,2),tab,tab,"[m2]"
                          write(90,'(I2,3A)') m_g7(l,3),tab,tab,"[m3]"
                          write(90,'(I2,3A)') m_g7(l,4),tab,tab,"[m4]"
                          write(90,'(I2,3A)') m_g7(l,5),tab,tab,"[m5]"
                          write(90,'(I2,3A)') m_g7(l,6),tab,tab,"[m6]"
                          write(90,'(I2,3A)') m_g7(l,7),tab,tab,"[m7]"
                          write(90,'(I2,3A)') n_g7(l,1),tab,tab,"[n1]"
                          write(90,'(I2,3A)') n_g7(l,2),tab,tab,"[n2]"
                          write(90,'(I2,3A)') n_g7(l,3),tab,tab,"[n3]"
                          write(90,'(I2,3A)') n_g7(l,4),tab,tab,"[n4]"
                          write(90,'(I2,3A)') n_g7(l,5),tab,tab,"[n5]"
                          write(90,'(I2,3A)') n_g7(l,6),tab,tab,"[n6]"
                          write(90,'(I2,3A)') n_g7(l,7),tab,tab,"[n7]"
                       case(8)
                          write(90,'(I2,3A)') m_g8(l,1),tab,tab,"[m1]"
                          write(90,'(I2,3A)') m_g8(l,2),tab,tab,"[m2]"
                          write(90,'(I2,3A)') m_g8(l,3),tab,tab,"[m3]"
                          write(90,'(I2,3A)') m_g8(l,4),tab,tab,"[m4]"
                          write(90,'(I2,3A)') m_g8(l,5),tab,tab,"[m5]"
                          write(90,'(I2,3A)') m_g8(l,6),tab,tab,"[m6]"
                          write(90,'(I2,3A)') m_g8(l,7),tab,tab,"[m7]"
                          write(90,'(I2,3A)') m_g8(l,8),tab,tab,"[m8]"
                          write(90,'(I2,3A)') n_g8(l,1),tab,tab,"[n1]"
                          write(90,'(I2,3A)') n_g8(l,2),tab,tab,"[n2]"
                          write(90,'(I2,3A)') n_g8(l,3),tab,tab,"[n3]"
                          write(90,'(I2,3A)') n_g8(l,4),tab,tab,"[n4]"
                          write(90,'(I2,3A)') n_g8(l,5),tab,tab,"[n5]"
                          write(90,'(I2,3A)') n_g8(l,6),tab,tab,"[n6]"
                          write(90,'(I2,3A)') n_g8(l,7),tab,tab,"[n7]"
                          write(90,'(I2,3A)') n_g8(l,8),tab,tab,"[n8]"
                     end select
                     write(90,*)
                     write(90,'(A)') "[file extension: "
     .                    //"distinguished by # of seed]"
                     write(90,'(I3,3A)') gen_seed,tab,tab,"[seed]"
                     close(90)
                     gen_seed=gen_seed+3
                  end do
               end do
            end do
         end do
      end do

 100  continue
      return
 101  continue
      write(*,*) "Now excute the excutables and generate data ..."
      call system('for i in run_*; do echo run in $i; cd $i; '//
     .     './main.exe; cd ..; done')
      write(*,*) "DONE!!!"
      return
      end program gen_tasks

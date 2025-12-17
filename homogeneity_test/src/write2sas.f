!     write simulated data to sas datastep
      subroutine data_to_sas(un)
      implicit none
      include 'parameters.f'
      include 'head_tail.f'
      integer :: i,un
      character(len=2) kk
      character(len=3) :: unn
      write(kk,'(I2)') 5*group
      write(unn,'(I3)') un
      if(header) then
         write(un,*) "data data_ori_"//unn//";"
         write(un,*) "input n1 - n"//adjustl(kk)//";"
         write(un,*) "cards;"
      else if (tailer) then
         write(un,*) ";"
         write(un,*) "run;"
         write(un,*)
         write(un,*)
         write(un,*)
      else
         write(un,'(3(I3,2X))') (int(m(0:2,i)),i=1,group)
         write(un,'(2(I3,2X))') (int(n(0:1,i)),i=1,group)
      end if
      return
      end subroutine data_to_sas


c$$$      program test_sas
c$$$      implicit none
c$$$      include 'parameters.f'
c$$$      integer :: io=100
c$$$      open(unit=io,file='../data/sas/test_1.sas',status='unknown',
c$$$     .     access='append')
c$$$      call data_stack_sas(io)
c$$$      call invoke_genmod(io)
c$$$      return
c$$$      end program test_sas

      

!     write stacked data into sas
      subroutine data_stack_sas(un)
      implicit none
      include 'parameters.f'
      integer :: i,j,k,l,un
      character(len=1) :: ii
      character(len=2) :: jj,kk
      character(len=3) :: unn
      character(len=1), dimension(2) :: res1
      character(len=1), dimension(6) :: res2
      res1=[character(len=1) :: "0","1"]
      res2=[character(len=1) :: "0","0","1","0","1","1"]
      write(unn,'(I3)') un
      write(un,*) "data data_stacked_"//unn//";"
      write(un,*) "set data_ori_"//unn//";"
!     bilateral stacked data
      do i=1,group
         write(ii,'(I1)') i
         l=1
         do j=3*i-2,3*i
            write(jj,'(I2)') j
            write(un,*) "id="//adjustl(jj)//"; n=n"//adjustl(jj)
     .       //"; res="//res2(l)//
     .           "; group="//ii//"; simcase=_n_; output;"
            l=l+1
            write(un,*) "id="//adjustl(jj)//"; n=n"//adjustl(jj)
     .           //"; res="//res2(l)//
     .           "; group="//ii//"; simcase=_n_; output;"
            l=l+1
         end do
      end do
!     unilateral stacked data
      k=3*group+1
      do i=1,group
         write(ii,'(I1)') i
         do l=1,2
            write(kk,'(I2)') k
            write(un,*) "id="//adjustl(kk)//"; n=n"//adjustl(kk)
     .           //"; res="//res1(l)//
     .           "; group="//ii//"; simcase=_n_; output;"
            k=k+1
         end do
      end do
      write(un,*) "drop n1 - n"//adjustl(kk)//";"
      write(un,*) "run;"
      write(un,*)
      write(un,*)
      write(un,*)

      return
      end subroutine data_stack_sas



!     run proc genmod in SAS with the stacked data
      subroutine invoke_genmod(un)
      implicit none
      include 'parameters.f'
      integer :: un
      character(len=3) :: unn
      character(len=180), dimension(7) :: g_contr
      g_contr=[character(len=180) :: 
     .     'group 1 -1;',       ! g=2
     .     'group 1 -1 0, group 0 1 -1;', ! g=3
     .     'group 1 -1 0 0, group 0 1 -1 0, group 0 0 1 -1;', ! g=4
     .     'group 1 -1 0 0 0, group 0 1 -1 0 0, group 0 0 1 -1 0, '
     .     //'group 0 0 0 1 -1;', ! g=5
     .     'group 1 -1 0 0 0 0, group 0 1 -1 0 0 0, 
     .     '//'group 0 0 1 -1 0 0, group 0 0 0 1 -1 0, '
     .     //'group 0 0 0 0 1 -1;', ! g=6
     .     'group 1 -1 0 0 0 0 0, group 0 1 -1 0 0 0 0, '
     .     //'group 0 0 1 -1 0 0 0, group 0 0 0 1 -1 0 0, '
     .     //'group 0 0 0 0 1 -1 0, group 0 0 0 0 0 1 -1;', ! g=7
     .     'group 1 -1 0 0 0 0 0 0, group 0 1 -1 0 0 0 0 0, ' 
     .     //'group 0 0 1 -1 0 0 0 0, group 0 0 0 1 -1 0 0 0, '
     .     //'group 0 0 0 0 1 -1 0 0, group 0 0 0 0 0 1 -1 0, '
     .     //'group 0 0 0 0 0 0 1 -1;' ! g=8
     .     ]
      write(unn,'(I3)') un
      write(un,*) "ods graphics off;"
      write(un,*) "ods exclude all;"
      write(un,*) "ods noresults;"
      write(un,*) "options nonotes;"
      write(un,*) "ods output Contrasts=contrasts_"//unn//";"
      write(un,*) "proc genmod data=data_stacked_"//unn//" descending;"
      write(un,*) "freq n;"
      write(un,*) "class group id;"
      write(un,*) "model res=group / link=logit dist=bin;"
      write(un,*) "repeated subject=id / type=un corrw;"
      write(un,*) "by simcase;"
      write(un,*) "contrast 'group' "//trim(g_contr(group-1))
      write(un,*) "run; quit;"
      write(un,*) "ods graphics on;"
      write(un,*) "ods exclude none;"
      write(un,*) "ods results;"
      write(un,*)
      write(un,*)
      write(un,*)
!     write the result from Contrasts to local file
      write(un,*) "filename doc 'C:\Users\cxma\Documents\jiazhou\';"
      write(un,*) "data _null_;"
      write(un,*) "set contrasts_"//unn//";"
      if(H0) then
         write(un,*) "file doc('p-value-H0_sas_"//unn//".out');"
      else
         write(un,*) "file doc('p-value-H1_sas_"//unn//".out');"
      end if
      write(un,*) "put ProbChisq;"
      write(un,*) "run;"
      write(un,*)
      write(un,*)
      write(un,*)
      return
      end subroutine invoke_genmod

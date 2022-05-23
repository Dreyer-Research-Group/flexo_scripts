program make_Makb
  implicit none
  integer :: nlines_1,nlines_2,ilin1,ilin2
  real :: makb1,makb2
  real, allocatable :: p1(:,:),p2(:,:)
  
  ! TODO:
  ! Here we assume Mz. Neet to generalize
  
  open(1,file="./P1/derivatives_piezo.dat")
  open(2,file="./P2/derivatives_piezo.dat")

  open(3,file="Makb.dat")
  write(3,*) "M_\alpha,\kappa\beta"
  write(3,*) "kappa  alpha  beta    M_akb"


  ! skip first line, check number of lines
  read(1,*); read(1,*) nlines_1; read(1,*)
  read(2,*); read(2,*) nlines_2; read(2,*)

  if (nlines_1 /=nlines_2) then
     write(*,*) "Different number of derivatives!!"
     stop
  end if

  ! d P_alpha,kappa beta / d q_gamma, in units of e*Bohr
  ! p1(1,:)=kappa, p1(2,:)=alpha, p1(3,:)=beta, p1(4,:)=gamma
  allocate(p1(5,nlines_1))
  allocate(p2(5,nlines_1))

  ! read in piezo
  do ilin1=1,nlines_1
     read(1,*) p1(:,ilin1)
     read(2,*) p2(:,ilin1)
  end do

  close(1)
  close(2)

  ! Construnct antisymmetric part
  do ilin1=1,nlines_1

     makb1=0.0; makb2=0.0
     
     ! Code golf version of levi civita 
     makb1=0.5*(3-1)*(1-p1(2,ilin1))*(p1(2,ilin1)-3)*p1(5,ilin1)

     if (abs(makb1) < 1.0d-20) cycle

     do ilin2=1,nlines_1
        
        ! Same \beta and \kappa
        if (p1(1,ilin1)==p2(1,ilin2).and.p1(3,ilin1)==p2(3,ilin2)) then

           makb2=0.5*(3-2)*(2-p2(2,ilin2))*(p2(2,ilin2)-3)*p2(5,ilin2)

           ! Only one other term, no need to keep looping
           if (abs(makb2) > 1.0d-20) exit

        end if
     
     end do

     ! TEST
     !write(*,*) 'kappa',p1(1,ilin1)
     !write(*,'(2i5,a14,3i5,e12.4e2,a14,3i5,e12.4e2)') int(p1(1,ilin1)),int(p1(3,ilin1)), 'first term:',3,1,int(p1(2,ilin1)),makb1,&
     !     & 'second term:',3,2,int(p2(2,ilin2)),makb2
     !write(*,*) ''

     
     !if (makb1+makb2  > 1.0d-20) 
     write(3,'(3i5,e18.8e2)') int(p1(1,ilin1)),3,int(p1(3,ilin1)),(makb1+makb2)*0.5

  end do


  close(3)

end program make_Makb

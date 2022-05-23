program sum_dist_col
  implicit none
  integer :: nlines,ilines,ilines2,linecount
  integer, allocatable :: at(:),dir(:),linetst(:)
  real, allocatable :: qvec(:,:),pind(:,:),outlines(:,:)
  

  open(1,file="totg_col.dat")
  read(1,*) nlines
  
  ! Read in quatities
  allocate(at(nlines))
  allocate(dir(nlines))
  allocate(qvec(3,nlines))
  allocate(pind(6,nlines))
  
  do ilines=1,nlines
     read(1,*) at(ilines),dir(ilines),qvec(:,ilines),pind(:,ilines)
  end do

  close(1)

  ! Sum up lines with the same q and direction
  allocate(linetst(nlines))
  allocate(outlines(10,nlines))
  linetst(:)=1
  outlines(:,:)=0.

  linecount=0
  do ilines=1,nlines

     ! Test if line has been summed already
     if (linetst(ilines)==0) cycle
     linetst(ilines)=0

     linecount=linecount+1
     
     outlines(1,linecount)=real(dir(ilines))
     outlines(2:4,linecount)=qvec(:,ilines)
     outlines(5:10,linecount)=pind(:,ilines)

     ! Sum lines
     do ilines2=ilines+1,nlines
        if (dir(ilines2)==dir(ilines).and.all(qvec(:,ilines2)==qvec(:,ilines))) then
           linetst(ilines2)=0
           outlines(5:10,linecount)=outlines(5:10,linecount)+pind(:,ilines2)
        end if
     end do
  end do

  ! Write output
  open(1,file="sum_totg_col.dat")
  do ilines=1,linecount
     write(1,"(i5,9e18.8e2)") int(outlines(1,ilines)),outlines(2:10,ilines)
  end do
  close(1)


end program sum_dist_col

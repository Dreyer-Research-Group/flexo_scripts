program flexo_shear_deriv
  implicit none
  integer :: nlines,ilines,iat,alpha,beta,gamma,lambda,ii
  integer, allocatable :: at(:),dir(:),linetst(:)
  real,parameter :: piezo_conv=57.21493204, flexo_conv=3.027682609,pz_prt_tol=1.0d-20
  real :: vol,sum_flexo
  real, allocatable :: qvec(:,:),pind(:,:),flexo(:)
  integer :: pz_convert
  
  ! TODO:
  ! Right now just long and trans (i.e., 1 q direction)
  ! Right now one beta at a time
  ! Right now Just 3 q pts per atom


  ! Read in the components:
  write(*,*) "For \mu_{\alpha\beta,\gamma\lambda}, type indicies"
  read(*,*) alpha, beta, gamma, lambda

  open(1,file="totg_col.dat")
  read(1,*) nlines,vol
  
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
  allocate(flexo(at(nlines)))
  flexo = 0.
 
  ! Assumes we have 4 qpts for each atom/direction. Should be:
  ! (-1,-1), (-1,1), (1,-1), (1,1)
  do ilines=1,nlines,4
        
     if (dir(ilines) /= beta) then
        write(*,*) "Skipping line", ilines, "since",dir(ilines), "/= beta",beta
        cycle
     end if

     flexo(at(ilines))=(1./2.)*(flexo_conv/vol)* &
          & (pind(2*alpha,ilines)+pind(2*alpha,ilines+3)-pind(2*alpha,ilines+1)-pind(2*alpha,ilines+2))/ &
          & (4.0*qvec(gamma,ilines+3)*qvec(lambda,ilines+3))
                 
     ! TEST
     !write(*,*) ilines,at(ilines),ialpha,dir(ilines),iq,jq,flexo(at(ilines),ialpha,dir(ilines),iq,jq)
     !stop
              
  end do ! ilines

  ! Write output
  open(1,file="shear_deriv_flexo.dat")
  write(1,*) "d2 P_alpha,kappa beta / d q_gamma d q_lambda" 
  write(1,*) "      kappa     alpha      beta    gamma     lambda     mu"
  sum_flexo=0.
  do iat=1,at(nlines)
     write(1,'(5i10,e18.8e2)') iat,alpha,beta,gamma,lambda,flexo(iat)
     sum_flexo=sum_flexo+flexo(iat)
  end do ! iat
  write(*,*) "Sum:",sum_flexo
  write(1,*) "Sum:",sum_flexo

  close(1)



end program flexo_shear_deriv

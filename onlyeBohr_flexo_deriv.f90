! SAME AS flexo_derivatives.f90 BUT NO OPTION FOR C/m^2
program flexo_derivatives
  implicit none
  integer :: nlines,ilines,iq,jq,ialpha,iat,idir,z_ct,ii,natom,gamma(3)
  integer, allocatable :: at(:),dir(:),linetst(:)
  real(16),parameter :: piezo_conv=57.21493204, flexo_conv=3.027682609,pz_prt_tol=1.0d-30
  real(16) :: vol
  real(16), allocatable :: qvec(:,:),pind(:,:),piezo(:,:,:,:),flexo(:,:,:,:,:),z_akb(:,:)
  real(16), allocatable :: zion(:)
  integer :: pz_convert
  
  ! TODO:
  ! Right now just long and trans (i.e., 1 q direction)
  ! Right now one beta at a time
  ! Right now Just 3 q pts per atom


  ! Should be convert to C/m^2 ?
!  write(*,*) "For the piezo derivative, (1) for C/m^2, (2) for e Bohr"
!  read(*,*) pz_convert

!  if (pz_convert==1) then
!     write(*,*) "Piezo in units of C/m^2"
!  elseif (pz_convert==2) then
!     write(*,*) "Piezo in units of e*Bohr"
!  else
!     write(*,*) "Error: must enter 1 or 2!"
!     stop
!  end if

  ! To identify qhich q's are included
  gamma(:)=0
     
  pz_convert=2


  open(1,file="totg_col.dat")
  read(1,*) nlines,vol,natom
  
  allocate(zion(natom))
  rewind(1)
  read(1,*) nlines,vol,natom,zion(:)

  
  ! Read in quatities
  allocate(at(nlines))
  allocate(dir(nlines))
  allocate(qvec(3,nlines))
  allocate(pind(6,nlines))
  allocate(z_akb(4,nlines))

  do ilines=1,nlines
     read(1,*) at(ilines),dir(ilines),qvec(:,ilines),pind(:,ilines)
  end do

  close(1)

  ! Sum up lines with the same q and direction
  allocate(piezo(at(nlines),3,3,3))
  allocate(flexo(at(nlines),3,3,3,3))
  piezo=0. ; flexo = 0.
  z_ct=1

  ! Assumes we have 3 qpts for each atom/direction
  do ilines=1,nlines,3
     
     ! Extract (real) BEC tensor
     do ii=1,3
        z_akb(1,z_ct)=at(ilines)
        z_akb(2,z_ct)=dir(ilines)
        z_akb(3,z_ct)=ii
        z_akb(4,z_ct)=pind(2*ii,ilines)

        z_ct=z_ct+1
     end do
        


     do ialpha=1,3 ! P_\alpha
        do iq=1,3 

           ! Check that calculation was done at given q direction
           ! Right now only for one q component
           if (abs(qvec(iq,ilines+1)) > 1.0d-6 ) then

              ! Identify that we have qpts in this dir
              gamma(iq)=1

              if (pz_convert==1) then
                 piezo(at(ilines),ialpha,dir(ilines),iq)=(piezo_conv/vol)* &
                      & ((-3./2.)*pind(2*(ialpha-1)+1,ilines)+2.*pind(2*(ialpha-1)+1,ilines+1) &
                      & +(-1./2.)*pind(2*(ialpha-1)+1,ilines+2))/qvec(iq,ilines+1)
              
              else
                 piezo(at(ilines),ialpha,dir(ilines),iq)=& !TEST, factor of vol??
                      & ((-3./2.)*pind(2*(ialpha-1)+1,ilines)+2.*pind(2*(ialpha-1)+1,ilines+1) &
                      & +(-1./2.)*pind(2*(ialpha-1)+1,ilines+2))/qvec(iq,ilines+1)

                 ! TEST
                 !if (at(ilines)==1.and.ialpha==1.and.dir(ilines)==3) then
                 !   write(*,*) "PIEZO",pind(2*(ialpha-1)+1,ilines),pind(2*(ialpha-1)+1,ilines+1),pind(2*(ialpha-1)+1,ilines+2)
                 !end if



              end if

              ! TEST
              !if (ialpha==3.and.iq==2) then
                 !write(*,*) ilines,ialpha,iq,at(ilines),dir(ilines),piezo(at(ilines),ialpha,dir(ilines),iq)
              !   write(*,*) pind(2*(ialpha-1)+1,ilines),pind(2*(ialpha-1)+1,ilines+1), &
              !        & pind(2*(ialpha-1)+1,ilines+2),qvec(iq,ilines+1)
              !   stop
              !end if
              
           end if

           do jq=1,3
              ! Check that calculation was done at given q direction
              if (abs(qvec(jq,ilines+1)) > 1.0d-6 .and. abs(qvec(iq,ilines+1)) > 1.0d-6) then
                 flexo(at(ilines),ialpha,dir(ilines),iq,jq)=(1./2.)*(flexo_conv/vol)* &
                      & (pind(2*ialpha,ilines)-2.*pind(2*ialpha,ilines+1)+pind(2*ialpha,ilines+2))/ &
                      & (qvec(iq,ilines+1)*qvec(jq,ilines+1))
                 
                 ! TEST
                 !write(*,*) pind(2*ialpha,ilines),2.*pind(2*ialpha,ilines+1),pind(2*ialpha,ilines+2)
                 !write(*,*) flexo(at(ilines),ialpha,dir(ilines),iq,jq) 
                 !ilines,at(ilines),ialpha,dir(ilines),iq,jq,flexo(at(ilines),ialpha,dir(ilines),iq,jq)
                 !stop
              end if
              
           end do ! jq
        end do !iq
     end do ! ialpha
  end do ! ilines

!write(*,*) flexo(1,1,1,1,1)
!  stop

  ! Write output
  open(1,file="derivatives_piezo.dat")
  if (pz_convert==1) then
     write(1,*) "d P_alpha,kappa beta / d q_gamma, in units of C/m^2"
  else
     write(1,*) "d P_alpha,kappa beta / d q_gamma, in units of e*Bohr"
  end if
  write(1,*) nlines !CHECK THIS
  write(1,*) "      kappa     alpha      beta    gamma     e"
  do iat=1,at(nlines)
     do ialpha=1,3
        do idir=1,3
           do iq=1,3
              if (gamma(iq)==1) then ! Make sure we have qpts in this direction
                 !if (abs(piezo(iat,ialpha,idir,iq)) > pz_prt_tol) then
                 write(1,'(4i10,e18.8e2)') iat,ialpha,idir,iq,piezo(iat,ialpha,idir,iq)
                !end if
              end if
              ! TEST
              !if (iat==1.and.ialpha==1.and.idir==3) then
              !   write(*,'(4i10,e18.8e2)') iat,ialpha,idir,iq,piezo(iat,ialpha,idir,iq)
              !end if

           end do !iq
        end do ! idir
     end do ! ialpha
  end do ! iat
  close(1)


  ! Write output
  open(1,file="derivatives_flexo.dat")
  write(1,*) "d2 P_alpha,kappa beta / d q_gamma d q_lambda" 
  write(1,*) "      kappa     alpha      beta    gamma     lambda     mu"
  do iat=1,at(nlines)
     do ialpha=1,3
        do idir=1,3
           do iq=1,3
              do jq=1,3
              if (abs(flexo(iat,ialpha,idir,iq,jq)) > 1.0d-10) then
                 write(1,'(5i10,e18.8e2)') iat,ialpha,idir,iq,jq,flexo(iat,ialpha,idir,iq,jq)
              end if
              end do !jq
           end do !iq
        end do ! idir
     end do ! ialpha
  end do ! iat
  close(1)

  ! Write output
  open(1,file="Z_akb.dat")
  write(1,*) "P(0)_\alpha\kappa\beta" 
  write(1,*) "      kappa     alpha      beta     Z*"
  do ii=1,nlines
     
     if (int(z_akb(2,ii))==int(z_akb(3,ii))) then
        write(1,'(3i10,e18.8e2)') int(z_akb(1,ii)),int(z_akb(2,ii)),int(z_akb(3,ii)),z_akb(4,ii)+zion(int(z_akb(1,ii)))
     else
        write(1,'(3i10,e18.8e2)') int(z_akb(1,ii)),int(z_akb(2,ii)),int(z_akb(3,ii)),z_akb(4,ii)
     end if

  end do
  close(1)



end program flexo_derivatives

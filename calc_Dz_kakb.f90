! Implement the levi civita symbol. Code golf lol
integer function eps(ii,jj,kk)
  implicit none
  integer,intent(in) :: ii,jj,kk

  eps=int(0.5*(ii-jj)*(jj-kk)*(kk-ii))

end function eps


program calc_Dz_kakb
  implicit none
  integer :: nlines_dm,nlines_dz,nlines_p0,ilin,found,natom
  integer :: kappa,alpha,kappa1,beta,gamma,ll,jj
  integer :: eps
  real :: d_gam_m,d_gam_p0,del_d_gam
  real, allocatable :: dm(:,:),dz(:,:),p0(:,:),xcart(:,:)
  logical :: round

  ! TODO:
  gamma=3 ! Here we assume Mz. Neet to generalize

  ! Run in Dz folder

  open(1,file="dMdt.dat")
  open(2,file="dZdt.dat")
  open(3,file="P0.dat")

  ! skip first line, check number of lines
  read(1,*); read(1,*) nlines_dm; read(1,*)
  read(2,*); read(2,*) nlines_dz; read(2,*)
  read(3,*); read(3,*) nlines_p0; read(3,*)

  ! allocate arrays
  allocate(dm(6,nlines_dm))
  allocate(dz(6,nlines_dz))
  allocate(p0(7,nlines_p0))

  ! read in info from files
  do ilin=1,max(nlines_dm,nlines_dz,nlines_p0)

     if (ilin<=nlines_dm) read(1,*) dm(:,ilin)
     if (ilin<=nlines_dz) read(2,*) dz(:,ilin)
     if (ilin<=nlines_p0) read(3,*) p0(:,ilin)
  end do

  close(1);close(2);close(3)

  ! Output file
  open(1,file="D_g_kakb.dat")
  write(1,*) "Dynamical matrix from Geometric magneization, D^\gamma_\kappa\alpha,\kappa'\beta"
  write(1,*) "  gamma kappa alpha kappa' beta        D_M             D_P0             delta D           total D"


  ! Get cartesian coodinates verus kappa
  allocate(xcart(natom,3))
  do ilin=1,nlines_p0
     kappa=int(p0(1,ilin))
     xcart(kappa,:)=p0(5:7,ilin)
  end do

  ! Calculate D
  natom=maxval(dm(1,:))
  do kappa=1,natom
     do alpha=1,3
        do kappa1=1,natom
           do beta=1,3

              !TEST
              write(*,*) 'kappa',kappa,'alpha',alpha,'kappa1',kappa1,'beta',beta
              
              ! Different components
              d_gam_m=0.0; d_gam_p0=0.0; del_d_gam=0.0

              ! First two terms: dM^z_\kappa'beta/dt_\kappa\alpha and dM^z_\kappa'beta/dt_\kappa\alpha
              found=0
              do ilin=1,nlines_dm
                 if (dm(1,ilin)==kappa.and.dm(2,ilin)==alpha.and.dm(3,ilin)==kappa1.and.dm(5,ilin)==beta) then

                    d_gam_m=d_gam_m+dm(6,ilin)
                    found=found+1

                    ! TEST
                    !write(*,*) "M1 is",int(dm(1:5,ilin)),dm(6,ilin)
                 end if
                    
                 if (dm(1,ilin)==kappa1.and.dm(2,ilin)==beta.and.dm(3,ilin)==kappa.and.dm(5,ilin)==alpha) then

                    d_gam_m=d_gam_m-dm(6,ilin)
                    found=found+1

                    ! TEST
                    !write(*,*) "M2 is",int(dm(1:5,ilin)),dm(6,ilin)

                 end if

                 if (found==2) exit

              end do ! ilin

              ! Make sure we got both terms
              if (found<2) then
                 write(*,*) "Could not find M terms for:",kappa,alpha,kappa1,beta
                 stop
              end if


              ! Second term: 0.5*\delta_kk' (eps^\gamma\alpha\l*P^(0)_l\kappa\beta-eps^\gamma\beta\l*P^(0)_l\kappa\alpha)
              !found=0                               
              if (kappa==kappa1) then

                 do ilin=1,nlines_p0

                    if (p0(1,ilin)==kappa.and.p0(3,ilin)==beta) then
                       ll=int(p0(2,ilin))
                       d_gam_p0=d_gam_p0 + 0.5*eps(gamma,alpha,ll)*p0(4,ilin)
                    end if

                    if (p0(1,ilin)==kappa.and.p0(3,ilin)==alpha) then
                       ll=int(p0(2,ilin))
                       d_gam_p0=d_gam_p0 - 0.5*eps(gamma,beta,ll)*p0(4,ilin)

                    end if

                 end do
              end if

              ! Electrical anharmonicity: \Delta D^z=0.5*eps(\gamma,j,l)*(t_\kappa'j*dP^(0)_l,\kappa'\beta)/dt_\kappa\alpha 
              !                                                          -t_\kappaj*dP^(0)_l,\kappa\alpha)/dt_\kappa'\beta)                
              do ilin=1, nlines_dz

                 if (dz(1,ilin)==kappa.and.dz(2,ilin)==alpha.and.dz(3,ilin)==kappa1.and.dz(5,ilin)==beta) then

                    do jj=1,3
                       ll=int(dz(4,ilin))
                       del_d_gam = del_d_gam + 0.5*eps(gamma,jj,ll)*xcart(kappa1,jj)*dz(6,ilin)
                    end do

                 end if

                 if (dz(1,ilin)==kappa1.and.dz(2,ilin)==beta.and.dz(3,ilin)==kappa.and.dz(5,ilin)==alpha) then

                    do jj=1,3
                       ll=int(dz(4,ilin))
                       del_d_gam = del_d_gam - 0.5*eps(gamma,jj,ll)*xcart(kappa,jj)*dz(6,ilin)
                    end do

                 end if

              end do


              ! Write to output file
              round=.true.
              if (round) then
                 write(1,'(5i6,4f18.8)') 3, kappa, alpha,kappa1,beta,d_gam_m,d_gam_p0, del_d_gam,d_gam_m+d_gam_p0+del_d_gam
              else
                 write(1,'(5i6,4e18.8e2)') 3, kappa, alpha,kappa1,beta,d_gam_m,d_gam_p0, del_d_gam,d_gam_m+d_gam_p0+del_d_gam
              end if

           end do ! beta
        end do ! kappa1
     end do ! alpha
  end do ! kappa

  close(1)



end program calc_Dz_kakb

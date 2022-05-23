! Calculate the LM part of the flexoelectric coefficients by taking numerical derivatives of the dynamical matrix
! Right now only works for Type I/II longitudinal and Type I transverse (equivalent to type II shear)
!
! Cyrus Dreyer, Stony Brook University and Flatiron CCQ
!
! To compile: gfortran -lblas -llapack -o twoddyn_pi.x twoddyn_psdoinv.f90
program twoddyn
  implicit none
  real, parameter :: pi = 4.*atan(1.),tol=1.e-10,flexo_conv=3.027682609
  integer, parameter :: alpha(6) = (/ 1, 2, 3, 1, 1, 2 /)
  integer, parameter :: beta(6)  = (/ 1, 2, 3, 2, 3, 3 /)
  integer :: ii,jj,ll,mm,nel,nqpt,mel,idir1,iat1,idir2,iat2,nat,nn,rho
  integer :: index(4),gam,qorder(4),kinvext
  integer, allocatable :: derivdir(:,:)
  real,allocatable :: cqpt(:,:),tau(:,:),amu(:),qmag(:),delq(:),kmatrix(:,:)
  real :: dynre,dynimg,rprim(3,3),acell(3),gprim(3,3),detinv,qpt(3),phasefac,dum,amutot,rebec,imgbec
  real :: ucvol
  complex :: dynel,phase,twoddm,oneddm
  complex,allocatable :: dynmat(:,:,:,:,:),sqbra(:,:,:,:,:),that(:,:,:,:,:),ttot(:,:,:,:),rndbra(:,:,:,:,:)
  complex,allocatable :: flxintstr(:,:,:,:,:),mul(:,:,:,:),bec(:,:,:),lambda(:,:,:,:),flexo_t(:,:,:,:,:)
  complex,allocatable :: piezointrel(:,:,:,:),phione(:,:,:,:,:),phitwo(:,:,:,:,:,:)

  !*************************************
  ! Read in data from DFPT calculation *
  !*************************************

  ! Dynamical matrix and other info
  open(1,file="dyn.dat")
  read(1,*) nqpt,nat

  !atomic weights
  allocate(amu(nat))
  read(1,*)
  read (1,*) amu(:)
  amutot=sum(amu)
  amu(:)=amu(:)/amutot
  
  ! Lattice info:
  read(1,*)
  read(1,*) acell(:)
  read(1,*) ucvol
  do ii=1,3
     read(1,*) rprim(ii,:)
     rprim(ii,:)=rprim(ii,:)*acell(ii)
  end do

  
  ! get gprim
  detinv = 1/(rprim(1,1)*rprim(2,2)*rprim(3,3) - rprim(1,1)*rprim(2,3)*rprim(3,2)&
       - rprim(1,2)*rprim(2,1)*rprim(3,3) + rprim(1,2)*rprim(2,3)*rprim(3,1)&
       + rprim(1,3)*rprim(2,1)*rprim(3,2) - rprim(1,3)*rprim(2,2)*rprim(3,1))
  gprim(1,1) = +detinv * (rprim(2,2)*rprim(3,3) - rprim(2,3)*rprim(3,2))
  gprim(2,1) = -detinv * (rprim(2,1)*rprim(3,3) - rprim(2,3)*rprim(3,1))
  gprim(3,1) = +detinv * (rprim(2,1)*rprim(3,2) - rprim(2,2)*rprim(3,1))
  gprim(1,2) = -detinv * (rprim(1,2)*rprim(3,3) - rprim(1,3)*rprim(3,2))
  gprim(2,2) = +detinv * (rprim(1,1)*rprim(3,3) - rprim(1,3)*rprim(3,1))
  gprim(3,2) = -detinv * (rprim(1,1)*rprim(3,2) - rprim(1,2)*rprim(3,1))
  gprim(1,3) = +detinv * (rprim(1,2)*rprim(2,3) - rprim(1,3)*rprim(2,2))
  gprim(2,3) = -detinv * (rprim(1,1)*rprim(2,3) - rprim(1,3)*rprim(2,1))
  gprim(3,3) = +detinv * (rprim(1,1)*rprim(2,2) - rprim(1,2)*rprim(2,1))

  gprim=gprim*2*pi

  
  ! get atomic positions
  allocate(tau(nat,3))
  read(1,*)
  do ii=1,nat
     read(1,*) tau(ii,:)
  end do

  !Dynamical matrix
  allocate(dynmat(nqpt,3,nat,3,nat))
  allocate(cqpt(nqpt,3))
  allocate(qmag(nqpt))

  allocate(derivdir(nqpt,6))
  allocate(delq(6))

  derivdir=0; dynmat=0.; qmag=0.; delq=0.

  !Read in dynmat
  read(1,*)
  do ii=1,nqpt
     read(1,*) qpt(:)

     ! Identify which derivative this qpt will be used for
     ! count backward so that, e.g. (0 1 1) not used for z and y derivs
     if (abs(qpt(1))<tol.and.abs(qpt(2))<tol.and.abs(qpt(3))<tol) then 
        gam = ii
     else
        do jj=6,1,-1         
           if (abs(qpt(alpha(jj)))>tol.and.abs(qpt(beta(jj)))>tol.and.derivdir(ii,jj)==0) then
              derivdir(ii,jj)=1 
              
              !TEST: 
              !write(*,*) ii,jj, alpha(jj),beta(jj), qpt(:)
              !write(*,*) ii,jj, derivdir(ii,jj)
             
              exit
           end if
        end do
     end if

     ! qpt in cartesian coordinates
     cqpt(ii,:)=matmul(gprim(:,:),qpt(:))
     qmag(ii)=dot_product(cqpt(ii,:),cqpt(ii,:))

     do jj=1,3*nat*3*nat
        read(1,*) idir1,iat1,idir2,iat2,dynre,dynimg

        ! compute phase
        phasefac=dot_product(cqpt(ii,:),(tau(iat2,:)-tau(iat1,:)))
        phase=cmplx(cos(phasefac),sin(phasefac))

        dynel=cmplx(dynre,dynimg)*phase
        dynmat(ii,idir1,iat1,idir2,iat2)=dynel

     end do   

  end do

  close(1)

  ! read in BEC. These should be put in dyn.dat at some point
  allocate(bec(nat,3,3))

  open(1,file="bec.dat")
  read(1,*)

  do ii=1,nat
     do jj=1,3
        do mm=1,3
          read(1,*) rebec,imgbec
          bec(ii,mm,jj)=cmplx(rebec,imgbec)
       end do
    end do
 end do

 close(1)


 ! q's should be in the correct order! See below
 ! qdel taken from the second element (after gamma)
 ! IS THIS ARRAY THE WRONG ORDER OF DIMENSIONS???
 do jj=1,6
    do ii=1,nqpt  

       ! TEST
       !write(*,*) ii,jj,derivdir(ii,jj),delq(jj)

       if (derivdir(ii,jj)==1.and.delq(jj)==0) then
          delq(jj)=abs(cqpt(ii,alpha(jj)))*abs(cqpt(ii,beta(jj)))
          
          !TEST
          !write(*,*) ii,jj,alpha(jj),delq(jj)
          !stop
          
          exit
       end if
    end do
 end do

 !TEST
write(*,*) "delq",delq(:)

 ! ********************************************
 ! Get pseudoinverse of force constant matrix *
 ! ********************************************
 
 allocate(kmatrix(3*nat,3*nat))

 ! External file provided?
 kinvext=0
 if (kinvext==1) then
    open(1,file="kinv.dat")
    do ii=1,nat*3
!       do jj=1,nat*3
          read(1,*) kmatrix(ii,:)
!       end do
    end do
    close(1)
 else

    call psdoinv(real(real(dynmat(gam,:,:,:,:))),kmatrix,nat)

 end if

 !TEST: print pseudoinverse
! do ii=1,3*nat
!    write(*,'(24e12.4e2)') kmatrix(ii,:)
! end do
 

  
  ! Allocations. Several don't need to be allocatable, lol
  allocate(sqbra(3,3,3,3,nat))
  allocate(lambda(3,3,3,nat))
  allocate(that(3,3,3,3,nat))
  allocate(ttot(3,3,3,3))
  allocate(flxintstr(3,3,3,3,nat))
  allocate(mul(3,3,3,3))
  allocate(flexo_t(3,3,3,3,nat))
  allocate(rndbra(3,3,3,3,nat))
  allocate(piezointrel(3,3,3,nat))
  allocate(phione(3,nat,3,nat,3))
  allocate(phitwo(3,nat,3,nat,3,3))
  sqbra=0.;lambda=0.;ttot=0.;that=0.;flxintstr=0.;mul=0.;flexo_t=0.
  rndbra=0.;piezointrel=0.;phione=0.;phitwo=0.

  !*********************************************
  !             Take derivatives               *
  !*********************************************
  do idir1=1,6
     if (maxval(derivdir(:,idir1))==0) cycle
     write(*,*) "taking the ",alpha(idir1),beta(idir1),"derivative"

     ! Find elements
     qorder=0
     do nn=1,nqpt
        if (derivdir(nn,idir1)==1.and.qorder(1)==0) then 
           qorder(1)=nn
        else if (derivdir(nn,idir1)==1.and.qorder(2)==0) then
           qorder(2)=nn
        else if (derivdir(nn,idir1)==1.and.qorder(3)==0) then
           qorder(3)=nn
        else if (derivdir(nn,idir1)==1.and.qorder(4)==0) then
           qorder(4)=nn
           exit
        end if
     end do

     do ii=1,3
        do jj=1,nat
           do ll=1,3
              do mm=1,nat
                 
                 if (idir1<=3) then ! Longitudinal/transverse (Type I)
    
                    ! SECOND DERIVATIVES ***************************************************************************

                    ! First order
                    !twoddm=-0.5*(dynmat(gam,ii,jj,ll,mm)-2.*dynmat(qorder(1),ii,jj,ll,mm) &
                    !     & +dynmat(qorder(2),ii,jj,ll,mm))/delq(idir1)**2     

                    ! Second order 
                    !twoddm=-0.5*(2.*dynmat(gam,ii,jj,ll,mm)-5.*dynmat(qorder(1),ii,jj,ll,mm) &
                    !     & +4.*dynmat(qorder(2),ii,jj,ll,mm)-dynmat(qorder(3),ii,jj,ll,mm))/delq(idir1)**2

                    ! Second order central, order: Gam, +1, -1, 2, -2
                    !twoddm=-0.5*(-(5./2.)*dynmat(gam,ii,jj,ll,mm)+(4./3.)*dynmat(qorder(1),ii,jj,ll,mm) &
                    !     & +(4./3.)*dynmat(qorder(2),ii,jj,ll,mm)-(1./12.)*dynmat(qorder(3),ii,jj,ll,mm) & 
                    !     & -(1./12.)*dynmat(qorder(4),ii,jj,ll,mm))/delq(idir1)
                    
                    ! First order, central
                    twoddm=-0.5*(dynmat(qorder(1),ii,jj,ll,mm)-2.*dynmat(gam,ii,jj,ll,mm) &
                         & +dynmat(qorder(2),ii,jj,ll,mm))/delq(idir1)!**2     

                    ! **********************************************************************************************


                    ! FIRST DERIVATIVES ****************************************************************************

                    ! First order
                    !oneddm=(-dynmat(gam,ii,jj,ll,mm)+dynmat(qorder(1),ii,jj,ll,mm))/delq(idir1)

                    ! Second order
                    !oneddm=(-1.5*dynmat(gam,ii,jj,ll,mm)+2.*dynmat(qorder(1),ii,jj,ll,mm) &
                    !     & -0.5*dynmat(qorder(2),ii,jj,ll,mm))/delq(idir1)
                    
                    ! Third order
                    !oneddm=((-11./6.)*dynmat(gam,ii,jj,ll,mm)+3.*dynmat(qorder(1),ii,jj,ll,mm) &
                    !     & -1.5*dynmat(qorder(2),ii,jj,ll,mm)+(1./3.)*dynmat(qorder(3),ii,jj,ll,mm))/delq(idir1)

                    ! First order, central
                    ! This is best for getting accurate first derivatives
                    oneddm=(-0.5*dynmat(qorder(1),ii,jj,ll,mm)+0.5*dynmat(qorder(2),ii,jj,ll,mm))/sqrt(delq(idir1))

                    oneddm=cmplx(real(aimag(oneddm)),real(real(oneddm)))
                    ! **********************************************************************************************
                    
                    lambda(ii,ll,alpha(idir1),jj)=lambda(ii,ll,alpha(idir1),jj)+oneddm
                    phione(ii,jj,ll,mm,idir1)=oneddm
                    !write(*,'(5i5,2e12.4e2)') ii,jj,ll,mm,idir1,phione(ii,jj,ll,mm,idir1)

                 else ! Shear (Type I)
                    ! The order needs to be: -1-1,11,-11,1-1
                    twoddm=-0.5*0.25*(dynmat(qorder(1),ii,jj,ll,mm)+dynmat(qorder(2),ii,jj,ll,mm) &
                         & -dynmat(qorder(3),ii,jj,ll,mm)-dynmat(qorder(4),ii,jj,ll,mm))/(delq(idir1))!**2)
                    
                 end if                 

                 ! Store second derivative of dynamical matrix
                 phitwo(ii,jj,ll,mm,alpha(idir1),beta(idir1))=twoddm

                 ! Sum for square bracket
                 sqbra(ii,ll,alpha(idir1),beta(idir1),jj)=sqbra(ii,ll,alpha(idir1),beta(idir1),jj)-twoddm
                
                 
              end do
           end do
        end do
     end do
  end do


  ! *******************************************
  ! * Assemble internal stress/strain tensors *
  ! *******************************************

  ! A lot of this should be combined, too many multi dimensional loops and arrays.
  ! In any case, the code is so fast anyway...

  ! Calculate the piezo internal strain tensor, \Gamma
!  do idir1=1,3 ! gamma,lambda
!     iat1=0
!     do jj=1,nat ! kappa
!        do ii=1,3 ! alpha
!           iat1=iat1+1
!           do ll=1,3 ! beta
!              iat2=0
!              do mm=1,nat ! kappa^prime
!                 do rho=1,3 ! rho
!                    iat2=iat2+1               

!                    piezointrel(ii,ll,alpha(idir1),jj)=piezointrel(ii,ll,alpha(idir1),jj)+ &
!                         & lambda(rho,ll,alpha(idir1),mm)*kmatrix(iat1,iat2)


!                    if (jj==1.and.ii==2.and.ll==2) write(*,'(10e20.12e2)') kmatrix(iat1,iat2)

!                 end do
!              end do
!           end do
!        end do
!     end do
!  end do
  
  ! Calculate the piezo internal strain tensor, \Gamma, Try two
  do idir1=1,3 ! gamma,lambda
     iat1=0
     do jj=1,nat ! kappa
        do ii=1,3 ! alpha
           iat1=iat1+1
           do ll=1,3 ! beta
              iat2=0
              do mm=1,nat ! kappa^prime
                 do rho=1,3 ! lambda
                    iat2=iat2+1               

                    piezointrel(ii,ll,alpha(idir1),jj)=piezointrel(ii,ll,alpha(idir1),jj)+ &
                         & lambda(rho,ll,alpha(idir1),mm)*kmatrix(iat1,iat2)

                    !TEST
                    !if (jj==1.and.ii==2.and.ll==2.and.rho==2.and.idir1==2) write(*,'(10e20.12e2)') kmatrix(iat1,iat2)

                 end do
              end do
           end do
        end do
     end do
  end do



  !TEST, (1111)^1
  !write(*,'(10e12.4e2)') real(piezointrel(2,1,1,:))
  !write(*,'(10e12.4e2)') real(phione(1,1,2,:,1))

  !TEST, (2222)^1
  !write(*,'(10e12.4e2)') real(piezointrel(2,2,2,:))
  !write(*,'(10e12.4e2)') real(phione(2,1,2,:,2))
  
  !TEST, (1122)^1
  !write(*,'(10e20.12e2)') real(piezointrel(2,2,2,:))
  !write(*,'(10e20.12e2)') real(lambda(2,2,2,:))
  

  do ii=1,3
     do ll=1,3
        do idir1=1,3
           do jj=1,nat
!              write(*,'(4i5,2e12.4e2)') ii,ll,idir1,jj,piezointrel(ii,ll,idir1,jj)
           end do
        end do
     end do
  end do


  ! Calculate round brackets
  do ii=1,nat ! kappa
     do jj=1,3 ! alpha
        do idir1=1,3 ! lambda
           do ll=1,3 ! beta
              do idir2=1,3 ! gamma
                 do mm=1,nat ! kappa prime
                    do rho=1,3 ! rho

                       rndbra(jj,idir1,ll,idir2,ii)=rndbra(jj,idir1,ll,idir2,ii)&
                            & + phione(jj,ii,rho,mm,idir1)*piezointrel(rho,ll,idir2,mm)
                       
                       !if (real(rndbra(jj,idir1,ll,idir2,ii))>1.0) then 
                       !   write(*,'(5i5,4e12.4e2)') jj,idir1,ll,idir2,ii,rndbra(jj,idir1,ll,idir2,ii),phione(jj,ii,rho,mm,idir1)
                       !end if

                       
                       
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do



  ! Get flexo force-response tensor, T, and Ttot
  do jj=1,nat ! kappa
     do ii=1,3 ! alpha
        do ll=1,3 ! beta
           do mm=1,3 ! gamma
              do rho=1,3 ! lambda
                 flexo_t(ii,ll,mm,rho,jj)=sqbra(ii,ll,mm,rho,jj) &
                      & +0.5*(rndbra(ii,mm,ll,rho,jj)+rndbra(ii,rho,ll,mm,jj))
                 ttot(ii,ll,mm,rho)=ttot(ii,ll,mm,rho)+flexo_t(ii,ll,mm,rho,jj)
              end do
           end do
        end do
     end do
  end do
  

  ! Get \hat{T}
  do idir1=1,6
     do ii=1,3
        do ll=1,3
           do jj=1,nat
              that(ii,ll,alpha(idir1),beta(idir1),jj)=sqbra(ii,ll,alpha(idir1),beta(idir1),jj) &
                   & -amu(jj)*ttot(ii,ll,alpha(idir1),beta(idir1))    
           end do
        end do
     end do
  end do

  ! Calculate the flexo internal strain tensor
  !iat1=0;iat2=0
  do idir1=1,6 ! gamma,lambda
     iat1=0
     do jj=1,nat ! kappa
        do ii=1,3 ! alpha
           iat1=iat1+1
           do ll=1,3 ! beta
              iat2=0
              do mm=1,nat ! kappa^prime
                 do rho=1,3 ! rho
                    iat2=iat2+1               
                    !kmat(ii,jj,rho,mm)=kmatrix(iat1,iat2)
                    flxintstr(ii,ll,alpha(idir1),beta(idir1),jj)=flxintstr(ii,ll,alpha(idir1),beta(idir1),jj) &
                         & +that(rho,ll,alpha(idir1),beta(idir1),mm)*kmatrix(iat1,iat2)                    

                    !TEST
                    !if (ii==1.and.ll==3.and.jj==1.and.idir1==3) then
                    !   write(*,'(2i5,2e12.4e2)') rho,mm,real(that(rho,ll,alpha(idir1),beta(idir1),mm)),real(kmatrix(iat1,iat2))
                    !end if

                 end do
              end do
           end do
        end do
     end do
  end do

!stop

! Calculate lattice contribution to the flexo coefficients
! Probably can be combined with the above loops, TODO
do idir1=1,6
   do ii=1,3 ! alpha
      do jj=1,3 ! beta

         do iat1=1,nat ! kappa
            do rho=1,3   
               mul(ii,jj,alpha(idir1),beta(idir1))=mul(ii,jj,alpha(idir1),beta(idir1)) &
                    & +flexo_conv*flxintstr(rho,jj,alpha(idir1),beta(idir1),iat1)*bec(iat1,ii,rho)/ucvol
               
               !TEST
               !if (idir1==3.and.ii==1.and.jj==1) then
               !   write(*,*) ii,jj,iat1,rho,mul(ii,jj,alpha(idir1),beta(idir1)),bec(iat1,ii,rho)
               !end if

            end do
         end do

      end do
   end do
end do



! ****************
! * Write output *
! ****************


open(2,file="flexoforce.dat")

! [ab,gd]^k
write(2,*) "[ab,gd]^k"
do idir1=1,6
   if (maxval(derivdir(:,idir1))==0) cycle
   do jj=1,nat
      write(2,"(a10,i5)") "kappa: ",jj
      do ii=1,3
         do ll=ii,3
!         do ll=1,3
           write(2,"(4i5,2e12.4e2)") ii,ll,alpha(idir1),beta(idir1),sqbra(ii,ll,alpha(idir1),beta(idir1),jj)
         end do
      end do
   end do
end do

!(ab,gd)^k
write(2,*) "(ab,gd)^k"
do idir1=1,6
   if (maxval(derivdir(:,idir1))==0) cycle
   do jj=1,nat
      write(2,"(a10,i5)") "kappa: ",jj
      do ii=1,3
!         do ll=ii,3
         do ll=1,3
            write(2,"(4i5,2e12.4e2)") ii,ll,alpha(idir1),beta(idir1),rndbra(ii,ll,alpha(idir1),beta(idir1),jj)
         end do
      end do
   end do
end do

! N^k_ab,gd
write(2,*)
write(2,*) "N^k_ab,gd"
do idir1=1,6
   if (maxval(derivdir(:,idir1))==0) cycle
   do jj=1,nat
      write(2,"(a10,i5)") "kappa: ",jj
      do ii=1,3
         do ll=ii,3
            write(2,"(4i5,2e12.4e2)") ii,ll,alpha(idir1),beta(idir1),flxintstr(ii,ll,alpha(idir1),beta(idir1),jj)
         end do
      end do
   end do
end do

! mul^I_ab,gl
write(2,*)
write(2,*) "mul^I_ab,gl"
do idir1=1,6
   if (maxval(derivdir(:,idir1))==0) cycle
   do ii=1,3
      do ll=ii,3
         write(2,"(4i5,2e12.4e2)") ii,ll,alpha(idir1),beta(idir1),mul(ii,ll,alpha(idir1),beta(idir1))
      end do
   end do
end do



close(2)




end program twoddyn




subroutine psdoinv(blkval,kmatrix,natom)
  
  implicit none

!Arguments----------------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real,intent(in) :: blkval(3,natom,3,natom)
 real,intent(out) :: kmatrix(3*natom,3*natom)

!Local variables------------------------------------
!scalars
 integer :: idirA,idirB,ier,ii1,ii2,ipertA,ipertB,ivarA,ivarB,ii,jj,asr
 character(len=500) :: direction,message
!arrays
 real, parameter :: tol6 = 1.e-6, tol8 = 1.e-8, tol5=1e-5
 real :: Amatr(3*natom-3,3*natom-3),Apmatr(3*natom,3*natom)
 real :: Cmatr(3*natom-3,3*natom-3)
 real :: Cpmatr(3*natom,3*natom),Nmatr(3*natom,3*natom),deviation(3,6)
 real :: d2cart(2,3*natom,3*natom)
 real :: d2cart_tmp(3,natom,3,natom),d2asr(3,natom,3,natom)
!for ZHPEV
 double precision :: zhpev2p(3*3*natom-2),eigvalp(3*natom)
 double precision ::  zhpev2(3*3*natom-5),eigval(3*natom-3)
 complex*16 :: Bpmatr((3*natom*(3*natom+1))/2),eigvecp(3*natom,3*natom),zhpev1p(2*3*natom-1)
 complex*16 :: Bmatr(((3*natom-3)*(3*natom-2))/2),eigvec(3*natom-3,3*natom-3),zhpev1(2*3*natom-4)
!***************************************************************

 !***************************
 ! Simple enforcement of ASR*
 !***************************
 asr=1
 d2asr=0
 d2cart_tmp=blkval
 if (asr==1) then
    ! Get matrix to subtract
    do ii=1,natom
       do idirA=1,3
          do idirB=1,3
             do jj=1,natom
                d2asr(idirA,ii,idirB,ii)=&
                     &         d2asr(idirA,ii,idirB,ii)+&
                     &         blkval(idirA,ii,idirB,jj) 
             end do
          end do
       end do
    end do
    ! Remove d2asr
    do jj=1,natom
       do idirB=1,3
          do ii=1,natom
             do idirA=1,3
                d2cart_tmp(idirA,ii,idirB,jj)= d2cart_tmp(idirA,ii,idirB,jj) - d2asr(idirA,ii,idirB,jj)
             end do
          end do
       end do
    end do
 end if


!try to get the displacement response internal strain tensor
!first need the inverse of force constant matrix
 d2cart = 0.
 do ipertA=1,natom
   do ii1=1,3
     ivarA=ii1+3*(ipertA-1)
     do ipertB=1,natom
       do ii2=1,3
         ivarB=ii2+3*(ipertB-1)
         d2cart(1,ivarA,ivarB)=d2cart_tmp(ii1,ipertA,ii2,ipertB) !blkval(ii1,ipertA,ii2,ipertB)
       end do
     end do
   end do
 end do

 kmatrix = d2cart(1,:,:)
 Apmatr(:,:)=kmatrix(:,:)

 !TEST
! write(*,*) "DynMat"
!  do ivarA=1,3*natom
!     write(*,'(24e20.12e2)') kmatrix(ivarA,:)
!  end do
! stop
     ! write(*,'(15e12.4e2)') blkval(:,1,1,1)

 
 Nmatr(:,:)=0.
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     if (mod(ivarA,3)==0 .and. mod(ivarB,3)==0)then
       Nmatr(ivarA,ivarB)=1.0
     end if
     if (mod(ivarA,3)==1 .and. mod(ivarB,3)==1)then
       Nmatr(ivarA,ivarB)=1.0
     end if
     if (mod(ivarA,3)==2 .and. mod(ivarB,3)==2)then
       Nmatr(ivarA,ivarB)=1.0
     end if
   end do
 end do

!TEST
! do ivarB=1,3*natom
!      write(*,'(24f10.1)') Nmatr(:,ivarB)
!end do
!stop

!starting the pseudoinervering processes
!then get the eigenvectors of the big matrix,give values to matrixBp
 Bpmatr=0.
 ii1=1
 do ivarA=1,3*natom
   do ivarB=1,ivarA
     Bpmatr(ii1)=cmplx(Nmatr(ivarB,ivarA))
     ii1=ii1+1
   end do
 end do


!TEST
! do ivarB=1,300
!    write(*,'(24e12.4e2)') real(Bpmatr(ivarB))
! end do
! stop


!Bpmatr(2,:) is the imaginary part of the force matrix
!then call the subroutines CHPEV and ZHPEV to get the eigenvectors
 call ZHPEV ('V','U',3*natom,Bpmatr,eigvalp,eigvecp,3*natom,zhpev1p,zhpev2p,ier)

 !TEST
! write (*,*) "eigenval"
!write(*,'(24e12.4e2)') eigvalp
! write (*,*) "eigenvec"
! do ivarB=1,3*natom
!      write(*,'(24e12.4e2)') real(eigvecp(:,ivarB))
!end do
!stop


!DEBUG
!the eigenval and eigenvec
!write(std_out,'(/,a,/)')'the eigenvalues and eigenvectors'
!do ivarA=1,3*natom
!write(std_out,'(/)')
!write(std_out,'(es16.6)')eigvalp(ivarA)
!end do
!do ivarA=1,3*natom
!write(std_out,'(/)')
!do ivarB=1,3*natom
!write(std_out,'(es16.6)')eigvecp(1,ivarB,ivarA)
!end do
!end do
!ENDDEBUG

!Then do the multiplication to get the reduced matrix,in two steps
!After this the force constant matrix is decouple in two bloks,
!acoustic and optical ones
 Cpmatr(:,:)=0.0
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     do ii1=1,3*natom
       Cpmatr(ivarA,ivarB)=Cpmatr(ivarA,ivarB)+real(eigvecp(ii1,ivarA))*Apmatr(ii1,ivarB)
     end do
   end do
 end do

 Apmatr(:,:)=0.0
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     do ii1=1,3*natom
       Apmatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)+Cpmatr(ivarA,ii1)*real(eigvecp(ii1,ivarB))
     end do
   end do
 end do

!TEST
! do ivarB=1,3*natom 
!      write(*,'(15e12.4e2)') real(Apmatr(:,ivarB))
!   end do
!stop

!DEBUG
!the blok diago
!write(std_out,'(/,a,/)')'matrixAp'
!do ivarA=1,3*natom
!write(std_out,'(/)')
!do ivarB=1,3*natom
!write(std_out,'(es16.6)')Apmatr(ivarA,ivarB)
!end do
!end do
!ENDDEBUG

!Check the last three eigenvalues whether too large or not
 ivarB=0
 do ivarA=3*natom-2,3*natom
   if (ABS(Apmatr(ivarA,ivarA))>tol5)then
     ivarB=1
     write(*,*) "Acoustic sum rule violation met"
     write(*,*) Apmatr(ivarA,ivarA)
!     stop
   end if
 end do

!Give the value of reduced matrix form Apmatr to Amatr
 do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
     Amatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)
   end do
 end do

!Now the reduced matrix is in the matrixA, the convert it
!first give the give the value of matixB from matrixA
 ii1=1
 do ivarA=1,3*natom-3
   do ivarB=1,ivarA
     Bmatr(ii1)=Amatr(ivarB,ivarA)
     ii1=ii1+1
   end do
 end do

!Call the subroutines CHPEV and ZHPEV to get the eigenvectors and the eigenvalues
 call ZHPEV ('V','U',3*natom-3,Bmatr,eigval,eigvec,3*natom-3,zhpev1,zhpev2,ier)

!Check the unstable phonon modes, if the first is negative then print
!warning message
 if(eigval(1)<-1.0*tol8)then
    write(*,*) 'Error 1:Unstable eigenvalue detected in force constant matrix at Gamma point',eigval(1)
    !stop
end if

!Do the matrix muplication to get pseudoinverse inverse matrix
 Cmatr(:,:)=0.0
 Amatr(:,:)=0.0
 do ivarA=1,3*natom-3
   Cmatr(ivarA,ivarA)=1.0/eigval(ivarA)
 end do

 do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
     do ii1=1,3*natom-3
       Amatr(ivarA,ivarB)=Amatr(ivarA,ivarB)+real(eigvec(ivarA,ii1))*Cmatr(ii1,ivarB)
     end do
   end do
 end do

!TEST 
! do ivarB=1,3*natom-3 
!    write(*,'(12e12.4e2)') Amatr(:,ivarB) 
! end do
!stop


!The second mulplication
 Cmatr(:,:)=0.0
 do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
     do ii1=1,3*natom-3
       Cmatr(ivarA,ivarB)=Cmatr(ivarA,ivarB)+ Amatr(ivarA,ii1)*real(eigvec(ivarB,ii1))
     end do
   end do
 end do

!DEBUG
!write(std_out,'(/,a,/)')'the pseudo inverse of the force matrix'
!do ivarA=1,3*natom
!write(std_out,'(/)')
!do ivarB=1,3*natom
!write(std_out,'(es16.6)')Cmatr(ivarA,ivarB)
!end do
!end do
!ENDDEBUG

!So now the inverse of the reduced matrix is in the matrixC
!now do another mulplication to get the pseudoinverse of the original
 Cpmatr(:,:)=0.0
 Apmatr(:,:)=0.0
 do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
     Cpmatr(ivarA,ivarB)=Cmatr(ivarA,ivarB)
   end do
 end do

!Now times the eigvecp
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     do ii1=1,3*natom
       Apmatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)+real(eigvecp(ivarA,ii1))*&
&       Cpmatr(ii1,ivarB)
     end do
   end do
 end do
 Cpmatr(:,:)=0.0
 do ivarA=1,3*natom
   do ivarB=1,3*natom
     do ii1=1,3*natom
       Cpmatr(ivarA,ivarB)=Cpmatr(ivarA,ivarB)+ Apmatr(ivarA,ii1)*real(eigvecp(ivarB,ii1))
     end do
   end do
 end do


!TEST
!do ivarA=1,3*natom
!   write(*,'(15e12.4e2)') kmatrix(:,ivarA)
!end do

!Now the inverse in in Cpmatr
 kmatrix(:,:)=Cpmatr(:,:)
!transfer the inverse of k-matrix back to the k matrix
!so now the inverse of k matrix is in the kmatrix
!ending the part for pseudoinversing the K matrix

!TEST
!write(*,*) "PI of DynMat"
!do ivarA=1,3*natom
!   write(*,'(24e20.12e2)') kmatrix(:,ivarA)
!end do
!stop


end subroutine psdoinv
!!***




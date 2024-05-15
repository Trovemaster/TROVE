!
!  This unit defines all specific routines for a fouratomic molecule of ch3oh type
!
module mol_zxy3
  use accuracy
  use moltype
  use lapack

  implicit none

  public ML_b0_zxy3,ML_coordinate_transform_zxy3,ML_rotsymmetry_ZXY3,ML_symmetry_transformation_ZXY3
  !
  private
 
  integer(ik), parameter :: verbose  = 5                       ! Verbosity level
  !
  contains
  !
  ! Procedures to define  ch3oh 
  ! 
  function ML_coordinate_transform_zxy3(src,ndst,direct) result (dst)
    !
    implicit none
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: dsrc(size(src))
    real(ark)                 :: tau1,tau2,tau3,x,y,z,phi1,phi2,phi3,alpha12,alpha13,alpha14,alpha23,alpha24,alpha34
    real(ark)                 :: cosbeta,beta312,beta412,beta413,cosa34
    character(len=cl)         :: txt
    !
    if (verbose>=7) write(out,"('ML_coordinate_transform_zxy3/start')") 
    !
    if (direct) then 
       !
       dsrc(:) = src(:) - molec%local_eq(:)
       !
    else 
       !
       dsrc(:) = src(:)
       !
    endif
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('MLcoordinate_transform_func: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'MLcoordinate_transform_func - bad coord. type'
       !
    case('R-BETA-TAU')
       !
       if (direct) then
         !
         tau2 = src( 8) 
         tau3 = src( 9)
         phi1 = tau2-tau3
         phi2 = 2.0_ark*pi-tau2
         phi3 = tau3
         !
         phi1 = tau2-tau3
         phi2 = 2.0_ark*pi-tau2
         phi3 = tau3
         !
         dst(1:7) = dsrc(1:7)
         dst(8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
         dst(9) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
         !
       else
         !
         x = 2.0_ark*pi/sqrt(3.0_ark)
         y = src(8)
         z = src(9)
         !
         phi1 =(sqrt(2.0_ark)*x+2.0_ark*y                )/sqrt(6.0_ark)
         phi2 =(sqrt(2.0_ark)*x-        y+sqrt(3.0_ark)*z)/sqrt(6.0_ark)
         phi3 =(sqrt(2.0_ark)*x-        y-sqrt(3.0_ark)*z)/sqrt(6.0_ark)
         !
         tau3  = phi3
         tau2 = 2.0_ark*pi-phi2
         !
         dst(8) = tau2
         dst(9) = tau3
         !
      endif
       !
    case('R-BETA-SYM')
       !
       if (direct) then
          !
          alpha12 = src(5)
          alpha13 = src(6)
          alpha23 = src(7)
          alpha14 = src(8)
          alpha24 = src(9)
          !
          if (size(src)==10) then 
            !
            alpha34 = src(10)
            !
            cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
            beta312 = aacos(cosbeta,txt)
            !
            cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
            beta412 = aacos(cosbeta,txt)
            !
            cosbeta = (cos(alpha34)-cos(alpha13)*cos(alpha14) )/(sin(alpha13)*sin(alpha14))
            beta413 = aacos(cosbeta,txt)
            !
          else
            !
            !alpha34 = calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
            !
            cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
            beta312 = aacos(cosbeta,txt)
            !
            cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
            beta412 = aacos(cosbeta,txt)
            !
            cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
            alpha34 = aacos(cosa34,txt)
            !
          endif 
          !
          dst(1)=dsrc(1)
          dst(2)=dsrc(2)
          dst(3)=dsrc(3)
          dst(4)=dsrc(4)
          dst(5)=dsrc(5)
          dst(6)=dsrc(6)
          dst(7)=dsrc(8)
          !
          phi1 = beta413 !2.0_ark*pi-(beta312+beta412)
          phi3 = beta312
          phi2 = beta412
          !
          dst(8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
          dst(9) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
          !
       else
          !
          !  write(out,*) dsrc
          !
          dst(1:4)=dsrc(1:4) + molec%local_eq(1:4)
          !
          call from_sym2alpha(src(5:9),dst(5:9),alpha34)
          !
          alpha12 = dst(5)
          alpha13 = dst(6)
          alpha23 = dst(7)
          alpha14 = dst(8)
          alpha24 = dst(9)
          !
          dst(5) = alpha12
          dst(6) = alpha13
          dst(7) = alpha23
          dst(8) = alpha14
          dst(9) = alpha24
          !
          if (size(dst)==10) then 
            !
            dst(10) = alpha34
            !
          endif 
          !
          !
       endif
       !
    end select 
    !
    if (verbose>=7) write(out,"('ML_coordinate_transform_zxy3/end')") 
    !
    !
  end function ML_coordinate_transform_zxy3





  ! Here we define structural parameters a0 for ABCD molecule,
  !
  subroutine ML_b0_zxy3(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
     !
     implicit none
     !
     integer(ik), intent(in) :: Npoints,Natoms
     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in) ,optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional  :: rho_ref
     real(ark),   intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
     !
     real(ark)    :: CM_shift,reYX
     real(ark)    :: mZ,mY,mX1,mX2,mX3,Mtotal,rho
     integer(ik) ::  i,i0,ix

      if (verbose>=4) write(out,"('ML_b0_zxy3/start')") 
      !
      if (size(molec%req)/=4) then
        write(out,"('ML_b0_zxy3: Nbonds must be 4 in this routine, not  ',i9)") size(molec%req)
        stop 'ML_b0_zxy3: wrong Nbonds '
      endif 
      !
      if (molec%Natoms/=5) then
        write(out,"('ML_b0_zxy3: Natoms must be 5 in this routine, not  ',i9)") molec%Natoms
        stop 'ML_b0_zxy3: wrong Natoms '
      endif 
      !
      mZ  = molec%AtomMasses(1)
      my  = molec%AtomMasses(2)
      !
      mX1 = molec%AtomMasses(3) 
      mX2 = molec%AtomMasses(4) 
      mX3 = molec%AtomMasses(5) 
      !
      Mtotal = mX1+mX2+mX3
      !
      reYX = molec%req(2)
      rho = molec%alphaeq(1)
      !
      !rho = pi-asin(2.0_ark/sqrt(3.0_ark)*sin(alpha/2.0_ark))
      !
      b0(1,1,0) = 0
      b0(1,2,0) = 0
      b0(1,3,0) = 0
      b0(2,1,0) = 0
      b0(2,2,0) = 0
      b0(2,3,0) = molec%req(1)
      b0(3,1,0) = reYX*sin(rho)
      b0(3,2,0) = 0
      b0(3,3,0) = reYX*cos(rho)
      b0(4,1,0) = -reYX*sin(rho)/2.0_ark
      b0(4,2,0) = sqrt(3.0_ark)*reYX*sin(rho)/2.0_ark
      b0(4,3,0) = reYX*cos(rho)
      b0(5,1,0) = -reYX*sin(rho)/2.0_ark
      b0(5,2,0) = -sqrt(3.0_ark)*reYX*sin(rho)/2.0_ark
      b0(5,3,0) = reYX*cos(rho)
      !
      do i = 1,3
        CM_shift = sum(b0(:,i,0)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
        b0(:,i,0) = b0(:,i,0) - CM_shift
      enddo 
      !
      ! We can also simulate the reference structure at each point of rho, if npoints presented
      !
      if (Npoints/=0) then
         !
         rho_ref = rho
         ! 
         if (.not.present(rho_borders).or..not.present(rho_ref)) then  
            write(out,"('ML_b0_zxy3: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_zxy3: rho_borders or rho_ref not specified '
         endif
         !
         do i = 0,npoints
            !
            rho = rho_i(i)
            !
            b0(1,1,i) = 0
            b0(1,2,i) = 0
            b0(1,3,i) = molec%req(1)
            b0(2,1,i) = 0
            b0(2,2,i) = 0
            b0(2,3,i) = 0
            b0(3,1,i) = reYX*sin(rho)
            b0(3,2,i) = 0
            b0(3,3,i) = reYX*cos(rho)
            b0(4,1,i) = -reYX*sin(rho)/2.0_ark
            b0(4,2,i) = sqrt(3.0_ark)*reYX*sin(rho)/2.0_ark
            b0(4,3,i) = reYX*cos(rho)
            b0(5,1,i) = -reYX*sin(rho)/2.0_ark
            b0(5,2,i) = -sqrt(3.0_ark)*reYX*sin(rho)/2.0_ark
            b0(5,3,i) = reYX*cos(rho)
            !
            do ix = 1,3
              CM_shift = sum(b0(:,ix,0)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
              b0(:,ix,i) = b0(:,ix,i) - CM_shift
            enddo 
            !
         enddo
         !
         !
      endif
      !
      !
      if (verbose>=3) then 
        !
        do i = 0,npoints
           !
           write(out,"(i6)") molec%natoms
           !
           write(out,"(/a,3x,3f14.8)") trim(molec%zmatrix(1)%name),b0(1,:,0) 
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(2)%name),b0(2,:,0)
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(3)%name),b0(3,:,0)
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(4)%name),b0(4,:,0)
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(5)%name),b0(5,:,0)
           !
        enddo
        !
      endif
      !
      if (verbose>=4) write(out,"('ML_b0_zxy3/end')") 
      !
  end subroutine ML_b0_zxy3




  ! Here we define the symmetry transformation of the Nmodes coordinates according the symmetry operations
  !
  subroutine ML_symmetry_transformation_ZXY3(ioper,natoms,src,dst)
    !
    integer(ik),intent(in)    :: ioper  ! group operation  
    integer(ik),intent(in)    :: natoms
    real(ark),intent(in)      :: src(natoms)
    real(ark),intent(out)     :: dst(natoms)
    real(ark)                 :: repres(12,2,2),a,b,e,o
    !
    integer(ik) :: nsrc
    !
    if (verbose>=6) write(out,"('ML_symmetry_transformation_ZXY3/start')") 
    !
    nsrc = size(src)
    !
    select case(trim(molec%coords_transform))
      !
    case default
       write (out,"('ML_coordinate_transform_ZXY2: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'ML_coordinate_transform_ZXY2 - bad coord. type'
       !
    case('R-BETA-SYM')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_ZXY3: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_ZXY3 - bad symm. type'
          !
       case('C3V','C3V(M)')
           !
         select case(ioper)
           !
         case (1) ! identity 
           !
           dst = src
           !
         case (3) ! (123)
           !
           dst(2) = src(3)
           dst(3) = src(4)
           dst(4) = src(2)
           !
           dst(5) = src(6)
           dst(6) = src(7)
           dst(7) = src(5)
           !
         case (2) ! (321)
           !
           dst(2) = src(4)
           dst(3) = src(2)
           dst(4) = src(3)
           !
           dst(5) = src(7)
           dst(6) = src(5)
           dst(7) = src(6)
           !
         case (6) ! (12)
           !
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(4)
           !
           dst(5) = src(6)
           dst(6) = src(5)
           dst(7) = src(7)
           !
         case (5) ! (13)
           !
           dst(2) = src(4)
           dst(3) = src(3)
           dst(4) = src(2)
           !
           dst(5) = src(7)
           dst(6) = src(6)
           dst(7) = src(5)
           !
         case (4) ! (23)
           !
           dst(2) = src(2)
           dst(3) = src(4)
           dst(4) = src(3)
           !
           dst(5) = src(5)
           dst(6) = src(7)
           dst(7) = src(6)
           !
         case default
           !
           write (out,"('symmetry_transformation_local: operation ',i8,' unknown')") ioper
           stop 'symmetry_transformation_local - bad operation. type'
           !
         end select 
         !
         a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
         !
         repres ( 1,:,:)= reshape((/ e, o,  & 
                                     o, e/),(/2,2/))
         !
         repres ( 3,:,:)= reshape((/-a,-b,  &
                                     b,-a/),(/2,2/))
         !
         repres ( 2,:,:)= reshape((/-a, b,  &
                                    -b,-a/),(/2,2/))
         !
         repres ( 4,:,:)= reshape((/ e, o,  &
                                     o,-e/),(/2,2/))
         !
         repres ( 6,:,:)= reshape((/-a, b,  &
                                     b, a/),(/2,2/))
         !
         repres ( 5,:,:)= reshape((/-a,-b,  &
                                    -b, a/),(/2,2/))
         !
         dst(8) = repres(ioper,1,1)*src(8)+repres(ioper,1,2)*src(9)
         dst(9) = repres(ioper,2,1)*src(8)+repres(ioper,2,2)*src(9)
         !
       end select
       !
    end select
    !
    if (verbose>=6) write(out,"('ML_symmetry_transformation_ZXY3/end')") 
    !
  end subroutine ML_symmetry_transformation_ZXY3


  ! Here we find the rotational symmetry from J,K
  !
  subroutine ML_rotsymmetry_ZXY3(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_ZXY3/start')") 
    !
    !
    select case(trim(molec%symmetry))
    case default
       write (out,"('ML_rotsymmetry_ZXY3: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_rotsymmetry_ZXY3 - bad symm. type'
       !
    case('C3V','C3V(M)')
       !
       gamma = 0 
       ideg = 1 
       !
       if (mod(K+3,3)==0.and.tau==0) gamma = 1 !; return
       if (mod(K+3,3)==0.and.tau==1) gamma = 2 !; return
       !
       if (mod(K+3,3)/=0.and.tau==0) then 
          gamma = 3 ; ideg = 1 
       endif 
       if (mod(K+3,3)/=0.and.tau==1) then 
          gamma = 3 ; ideg = 2
       endif 
       !
    end select
    !
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_ZXY3/end')") 
    !
    !
  end subroutine ML_rotsymmetry_ZXY3


  subroutine from_sym2alpha(s,local,alpha34)
     !
     ! obtaining 5 angles from 5 internal sym. coordinates
     !
     real(ark),intent(in)  :: s(5)
     real(ark),intent(out) :: local(5),alpha34
     !
     real(ark) :: s_t(5),eps(5)
     real(ark) :: rjacob(5,5),s_l(5),s_r(5)
     real(ark) :: am(5,5),bm(5),cm(5),alpha34_
     !
     real(rk) :: a(5,5),b(5,1)

     real(ark) :: stadev_old,stadev,ssq,stadev_best,h
     !
     integer(ik) :: iter,itmax,k,irow,icolumn,ierror
       !
       local(:) = molec%local_eq(5:)
       !
       iter = 0
       stadev_old = 1.e10
       stadev    =  1.e10
       !
       stadev_best = sqrt(small_)*0.001_ark
       itmax = 500
       !
       ! Initial value for alpha10
       !
       outer_loop: & 
       do while( iter<itmax .and. stadev>stadev_best )   
         !
         iter = iter + 1
         ssq=0
         !
         ! Caclulate the function 
         !
         call calc_sym_from_alpha(local,s_r,alpha34)
         !
         eps(:) = s(:)-s_r(:)
         !
         do k = 1,5
           !
           h = 1.e-3*abs(local(k)) ; if (h<1e-12) h = 1e-7
           !
           local(k) = local(k) + h 
           !
           call calc_sym_from_alpha(local,s_r,alpha34_)
           !
           local(k) = local(k) - h - h 
           !
           call calc_sym_from_alpha(local,s_l,alpha34_)
           !
           rjacob(:,k)  = ( s_r(:)-s_l(:))/h*0.5_ark
           !
           local(k) = local(k) + h
           !
         enddo 
         !
         ssq=sqrt(sum(eps(:)**2))
         !
         ! We constract a set of linear equations A x = B
         !
         ! form A matrix 
         !
         do irow=1,5       !==== row-...... ====!
           do icolumn=1,irow    !==== column-....====!
             am(irow,icolumn)=sum(rjacob(:,icolumn)*rjacob(:,irow))
             am(icolumn,irow)=am(irow,icolumn)
           enddo
         enddo
         !
         ! form B matrix 
         !
         do irow=1,5       !==== row-...... ====!
           bm(irow)=sum(eps(:)*rjacob(:,irow))
         enddo   
         !
         ! Solve the set of linear equations 
         !
         call MLlinurark(5,am,bm(:),cm,ierror)
         !
         if (ierror>0) then
           !
           a = am ; b(:,1) = bm(:) 
           !
           call lapack_gelss(a(:,:),b(:,:))
           !
           cm(:) = b(:,1)
           !
         endif
         !
         local(:)=local(:)+cm(:)
         !
         stadev=ssq/sqrt(5.0_ark)
         !
         stadev_old=stadev
         !
       enddo  outer_loop ! --- iter
       !
       if (iter==itmax) then
          write(out,"('from_sym2alpha: could not find solution after ',i8,' iterations')") iter
          stop 'from_sym2alphaII: could not find solution'
       endif 

       !
       contains 


    subroutine calc_sym_from_alpha(src,dst,alpha34)

      real(ark),intent(in)  :: src(5)
      real(ark),intent(out) :: dst(5),alpha34
      real(ark)             :: alpha12,alpha13,alpha23,alpha14,alpha24,cosbeta,beta312,beta412,beta413,cosa34,phi1,phi2,phi3
      character(len=cl)     :: txt
       !
       txt = 'calc_sym_from_alpha'
       !
       alpha12 = src(1)
       alpha13 = src(2)
       alpha23 = src(3)
       alpha14 = src(4)
       alpha24 = src(5)
       !
       !alpha34 = calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
       !
       cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
       beta312 = aacos(cosbeta,txt)
       !
       cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
       beta412 = aacos(cosbeta,txt)
       !
       cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
       alpha34 = aacos(cosa34,txt)
       !
       cosbeta = (cos(alpha34)-cos(alpha13)*cos(alpha14) )/(sin(alpha13)*sin(alpha14))
       beta413 = aacos(cosbeta,txt)
       !
       dst(1)=alpha12-molec%local_eq(5)
       dst(2)=alpha13-molec%local_eq(6)
       dst(3)=alpha14-molec%local_eq(8)
       !
       phi1 = beta413 ! 2.0_ark*pi-(beta312+beta412)
       phi3 = beta312
       phi2 = beta412
       !
       dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
       dst(5) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
       !
      end subroutine calc_sym_from_alpha
      !
  end subroutine from_sym2alpha


end module mol_zxy3

!
!  This unit defines all specific routines for a fouratomic molecule of XY3 type
!
module mol_xy3
  use accuracy
  use moltype
  use lapack

  implicit none

  public ML_b0_XY3,ML_coordinate_transform_XY3,ML_symmetry_transformation_XY3,ML_rotsymmetry_XY3,ML_MEP_NH3

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level

  contains


!
! Here we define structural parameters for rigid XY3 molecule
!
  subroutine ML_b0_XY3(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
     !
     !
     integer(ik), intent(in) :: Npoints,Natoms
     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in) ,optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional :: rho_ref
     real(ark),   intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
     !
     real(ark)                :: c0(2,2),a0(molec%Natoms,3),cosr,sinr,re14,mH1,mH2,mH3,mX,rho,Mtotal,cosphi,phi,Inert(3),req(6)
     real(ark)                :: transform1(3,3),alpha,rho0,CM_shift,tmat(3,3),transform(3,3),transform0(3,3),a_t
     real(ark)                :: sinrho,hx,hy,unit(3)
     integer(ik)             :: i,n,ix,jx,ioper
     character(cl)           :: method
     !
     real(ark)    :: rho_ark,Inert0(3),resid(3,3)
     !
     integer(ik)             :: ijk(3,2),i0,j0,ipar,istep
     !
      if (verbose>=4) write(out,"('ML_b0_XY3/start')") 


      if (size(molec%req)/=3) then
        write(out,"('ML_b0_XY3: Nbonds must be 3 in this routine, not  ',i8)") size(molec%req)
        stop 'ML_b0_XY3: wrong Nbonds '
      endif 

      if (molec%Natoms/=4) then
        write(out,"('ML_b0_XY3: Natoms must be 4 in this routine, not  ',i8)") molec%Natoms
        stop 'ML_b0_XY3: wrong Natoms '
      endif 

      if (molec%req(1)/=molec%req(2).or.molec%req(1)/=molec%req(3).or.molec%req(2)/=molec%req(3)) then
        write(out,"('ML_b0_XY3: req-s must be equal: ',3f14.6)") molec%req(1:3)
        stop 'ML_b0_XY3: req-s must be equal: '
      endif 
      !
      alpha = molec%alphaeq(1)
      if (any(molec%alphaeq(:)/=alpha)) then
        write(out,"('ML_b0_XY3: alphaeq-s must be equal: ',3f14.6)") molec%alphaeq(:)
        stop 'ML_b0_XY3: alphaeq-s must be equal: '
      endif 
      !
      mH1 = molec%AtomMasses(2) 
      mH2 = molec%AtomMasses(3) 
      mH3 = molec%AtomMasses(4) 
      mX  = molec%AtomMasses(1)
      !
      Mtotal = mH1+mH2+mH3
      !
      if (any(molec%AtomMasses(2:4)>molec%AtomMasses(1))) then
        write(out,"('ML_b0_XY3: masses-s are given in wrong order, must be M m m m: ',4f14.6)") molec%AtomMasses(:)
        !stop 'ML_b0_XY3: ,masses are in wrong order'
      endif 
      !
      re14 = molec%req(1)
      alpha = molec%alphaeq(1)
      rho = pi-asin(2.0_ark/sqrt(3.0_ark)*sin(alpha/2.0_ark))

      select case(trim(molec%meptype))
           !
      case('MEP_H3OP')
           !
           rho_ark = rho
           !
           re14 = ML_mep_oh3p(rho_ark)
           alpha = 2.0_ark/3.0_ark*pi
           !
           !rho = 0
           !
      case('MEP_LAF3')
           !
           rho_ark = rho
           !
           re14 = ML_mep_oh3p(rho_ark)
           alpha = 2.0_ark/3.0_ark*pi
           !
      end select 
      !
      call generate_structure(rho,b0(1:4,1:3,0))
      !
      !cosr = cos(rho)
      !sinr = sin(rho)
      !
      !b0(2,1,0) = re14*sinr
      !b0(2,2,0) = 0
      !b0(2,3,0) = mX*re14*cosr/(Mtotal+mX)
      !b0(3,1,0) = -re14*sinr/2.0_ark
      !b0(3,2,0) = sqrt(3.0_ark)*re14*sinr/2.0_ark
      !b0(3,3,0) = mX*re14*cosr/(Mtotal+mX)
      !b0(4,1,0) = -re14*sinr/2.0_ark
      !b0(4,2,0) = -sqrt(3.0_ark)*re14*sinr/2.0_ark
      !b0(4,3,0) = mX*re14*cosr/(Mtotal+mX)
      !b0(1,1,0) = 0
      !b0(1,2,0) = 0
      !b0(1,3,0) = -Mtotal*re14*cosr/(Mtotal+mX)
      !
      ! Find center of mass
      !
      if (any(molec%AtomMasses(2:4)/=mH1)) then 
        !
        do n = 1,3
          CM_shift = sum(b0(:,n,0)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
          b0(:,n,0) = b0(:,n,0) - CM_shift
        enddo 
        !
        call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,0),transform)
        !
        Inert(1) = sum(molec%AtomMasses(:)*( b0(:,2,0)**2+ b0(:,3,0)**2) )
        Inert(2) = sum(molec%AtomMasses(:)*( b0(:,1,0)**2+ b0(:,3,0)**2) )
        Inert(3) = sum(molec%AtomMasses(:)*( b0(:,1,0)**2+ b0(:,2,0)**2) )
        !
        ! Second Eckart equation
        ! 
        do ix = 1,3 
           do jx = 1,3 
              if (ix/=jx) then  
                 !
                 a_t =  sum(molec%AtomMasses(:)*b0(:,ix,0)*b0(:,jx,0) )
                 !
                 if (abs(a_t)>100.0_rk*small_) then 
                     write(out,"('ML_b0_XY3: b0 is not a solution of Eckart 2 for ix,jx =',2i4,d18.8)") ix,jx,a_t
                     stop 'ML_b0_XY3: b0 is not solution of Eckart2'
                 endif
              endif
           enddo
           !
        enddo
        !
      endif
      !
      ! We can also simulate the reference structure at each point of rho, if npoints presented
      !
      if (Npoints/=0) then
         ! 
         if (.not.present(rho_borders).or..not.present(rho_ref)) then  
            write(out,"('ML_b0_XY3: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_XY3: rho_borders or rho_ref not specified '
         endif
         !
         rho_ref = 0 
         rho0 = pi*0.5_ark
         !
         if (trim(molec%coords_transform)=='R-2D-RHO') then 
            !
            rho0 = 0
            rho_ref = pi*0.5_ark
            !
         endif 
         !
         ! define the rho-type coordinate 
         !
         transform0 = 0 
         !
         forall (ix=1:3) transform0(ix,ix) = 1.0_ark
         !
         call generate_structure(rho0,a0(1:4,1:3))
         !
         ! symmetric XY3 case 
         !
         if (all(molec%AtomMasses(2:4)==mH1)) then 
           !
           do i = 0,npoints
              !
              rho = rho0+rho_i(i)
              !
              re14 = molec%req(1)
              !
              select case(trim(molec%coords_transform))
              !
              case('R-S-DELTA-MEP','R-PHI-DELTA-MEP','R-SYMPHI-DELTA-MEP')
                 !
                 select case(trim(molec%meptype))
                 case default
                      write (out,"('ML_b0_XY3: MEP type ',a,' unknown')") trim(molec%meptype)
                      stop 'ML_b0_XY3 - bad MEP type'
                      !
                 case('MEP_H3OP')
                      !
                      re14 = ML_mep_oh3p(rho)
                      !
                 end select 
                 !
              end select
              !
              call generate_structure(rho,b0(1:4,1:3,i))
              !
           enddo
           !
         else
           !
           ! non symmetric XY3 case (isotopologues)
           !
           i0 = npoints/2
           !
           re14 = molec%req(1)
           !
           do ipar = 0,1
             !
             istep = (-1)**(ipar+2)
             !
             i = npoints/2
             !
             rho = rho0+rho_i(i)
             !
             call generate_structure(rho,b0(1:4,1:3,i))
             !
             a0(:,:) = b0(:,:,i)
             !
             call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,transform0)
             !
             Inert0(1) = sum(molec%AtomMasses(:)*( a0(:,2)**2+ a0(:,3)**2) )
             Inert0(2) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,3)**2) )
             Inert0(3) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,2)**2) )
             !
             !do ix = 1,molec%Natoms
             !   a0(ix,:) = matmul(transpose(transform),b0(ix,:,i))
             !enddo
             !
             !b0(:,:,i) = a0(:,:)
             !
             do j0 = 1,npoints/2
                !
                i = i + istep
                !
                rho = rho0+rho_i(i)
                !
                call generate_structure(rho,b0(1:4,1:3,i))
                !
                do ix = 1,molec%Natoms
                   a0(ix,:) = matmul(transpose(transform0),b0(ix,:,i))
                enddo
                !
                call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,transform1)
                !
                transform = matmul(transform1,transform0)
                !
                do ix = 1,molec%Natoms
                   a0(ix,:) = matmul(transpose(transform),b0(ix,:,i))
                enddo
                !
                b0(:,:,i) = a0(:,:)
                !
                do ix = 1,3 
                   do jx = 1,3 
                      resid(ix,jx) =  sum(molec%AtomMasses(:)*b0(:,ix,i)*b0(:,jx,i) )
                   enddo
                   !
                enddo
                !
                Inert(1) = sum(molec%AtomMasses(:)*( b0(:,2,i)**2+ b0(:,3,i)**2) )
                Inert(2) = sum(molec%AtomMasses(:)*( b0(:,1,i)**2+ b0(:,3,i)**2) )
                Inert(3) = sum(molec%AtomMasses(:)*( b0(:,1,i)**2+ b0(:,2,i)**2) )
                !
                select case (trim(molec%coords_transform))
                  !
                case ('R-A2-A3-TAU')
                  !
                  phi = rho_i(i)
                  !
                  hx = -sqrt(1.0_ark-2.0_ark*cos(phi) )/(-1.0_ark+cos(phi))
                  hy = -cos(phi)/(-1.0_ark+cos(phi))
                  !
                  alpha = atan2(hx,hy)
                  !
                  sinrho = 2.0_ark*sin(alpha*0.5_ark)/sqrt(3.0_ark)
                  !
                  if ( sinrho>=1.0_ark ) then 
                     rho = 0.5_ark*pi
                  else 
                     rho = asin(sinrho)
                  endif 
                  !
                  if (phi>pi) rho = pi - rho
                  !
                end select
                !
                ijk(1,1) = 2 ; ijk(1,2) = 3
                ijk(2,1) = 3 ; ijk(2,2) = 1
                ijk(3,1) = 1 ; ijk(3,2) = 2
                !
                method = 'ULEN'
                !
                a0 = b0(:,:,i)
                !
                !call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i),transform1,method=method)
                !
                loop_xyz : do ix =1,3
                  !
                  do ioper = 1,4
                    !
                    phi = real(ioper-1,ark)*pi*0.5_ark
                    !
                    tmat(:,:) = 0 
                    !
                    tmat(ix,ix) = 1.0_ark
                    !
                    tmat(ijk(ix,1),ijk(ix,1)) = cos(phi)
                    tmat(ijk(ix,1),ijk(ix,2)) = sin(phi)
                    tmat(ijk(ix,2),ijk(ix,1)) =-sin(phi)
                    tmat(ijk(ix,2),ijk(ix,2)) = cos(phi)
                    !
                    transform1 = matmul(tmat,transform)
                    !
                    forall(ix=1:3) unit(ix) = dot_product(transform0(ix,:),transform1(ix,:))
                    !
                    unit = unit - 1.0_ark
                    !
                    if ( all( abs( unit(:) )<1.0e-1  ) ) exit loop_xyz
                    !
                  enddo
                enddo loop_xyz
                !
                transform0 = matmul(tmat,transform)
                !
                forall(ix=1:4) b0(ix,:,i) = matmul(transpose(tmat),b0(ix,:,i))
                !
                Inert(1) = sum(molec%AtomMasses(:)*( b0(:,2,i)**2+ b0(:,3,i)**2) )
                Inert(2) = sum(molec%AtomMasses(:)*( b0(:,1,i)**2+ b0(:,3,i)**2) )
                Inert(3) = sum(molec%AtomMasses(:)*( b0(:,1,i)**2+ b0(:,2,i)**2) )
                !
                write(out,"(f15.8,1x,3f15.8)") rho,Inert(1:3)
                !
             enddo
             !
           enddo  
           !
         endif
         !
      endif 
      !
      do i = 0,npoints
         !
         ! Second Eckart equation
         ! 
         do ix = 1,3 
            do jx = 1,3 
               if (ix/=jx) then  
                  !
                  a_t =  sum(molec%AtomMasses(:)*b0(:,ix,i)*b0(:,jx,i) )
                  !
                  if (abs(a_t)>100.0_ark*small_) then 
                      write(out,"('ML_b0_XY3: b0 is not a solution of Eckart 2 for ix,jx =',2i4,d18.8)") ix,jx,a_t
                      stop 'ML_b0_XY3: b0 is not solution of Eckart2'
                  endif
               endif
            enddo
            !
         enddo
         !
         if (verbose>=3) then 
           !
           write(out,"(i6)") molec%natoms
           !
           write(out,"(/a,3x,3f14.8)") trim(molec%zmatrix(1)%name),b0(1,:,i) 
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(2)%name),b0(2,:,i)
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(3)%name),b0(3,:,i)
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(4)%name),b0(4,:,i)
           !
         endif
         !
         if (verbose>=5) then 
            write(out,"('b0',i4,12f12.8)") i,b0(:,:,i)
         endif
         ! 
      enddo
      !
      if (verbose>=4) write(out,"('ML_b0_XY3/end')") 
      !
      contains 
      !
      subroutine generate_structure(rho,b0)
        !
        real(ark),intent(in)  :: rho
        real(ark),intent(out) :: b0(4,3)
        real(ark) :: sinr,cosr,CM_shift,theta
         !
         select case(trim(molec%coords_transform))
           !
           case('R-THETA-TAU')
              !
              theta = molec%alphaeq(2)
              !
              b0(1,1) = 0.0_ark
              b0(1,2) = 0.0_ark
              b0(1,3) = 0.0_ark
              !
              b0(2,1) = 0.0_ark
              b0(2,2) = 0.0_ark
              b0(2,3) = molec%req(1)
              !
              b0(3,1) = molec%req(2)*sin(theta)*cos(rho*0.5_ark)
              b0(3,2) =-molec%req(2)*sin(theta)*sin(rho*0.5_ark)
              b0(3,3) = molec%req(2)*cos(theta)
              !
              b0(4,1) = molec%req(2)*sin(theta)*cos(rho*0.5_ark)
              b0(4,2) = molec%req(2)*sin(theta)*sin(rho*0.5_ark)
              b0(4,3) = molec%req(2)*cos(theta)            
              !
           case default 
              !
              cosr = cos(rho)
              sinr = sin(rho)
              !
              if (all(molec%AtomMasses(2:4)==mH1)) then 
                !
                b0(2,1) = re14*sinr
                b0(2,2) = 0.0_ark
                b0(2,3) = mX*re14*cosr/(Mtotal+mX)
                b0(3,1) = -re14*sinr/2.0_ark
                b0(3,2) = sqrt(3.0_ark)*re14*sinr/2.0_ark
                b0(3,3) = mX*re14*cosr/(Mtotal+mX)
                b0(4,1) = -re14*sinr/2.0_ark
                b0(4,2) = -sqrt(3.0_ark)*re14*sinr/2.0_ark
                b0(4,3) = mX*re14*cosr/(Mtotal+mX)
                b0(1,1) = 0.0_ark
                b0(1,2) = 0.0_ark
                b0(1,3) = -Mtotal*re14*cosr/(Mtotal+mX)
                !
              else
                !
                select case(trim(molec%frame))
                  !
                case default  !  R1-Z-R2-X-Y 
                  !
                  b0(2,3) = re14*sinr
                  b0(2,1) = 0.0_ark
                  b0(2,2) = mX*re14*cosr/(Mtotal+mX)
                  b0(3,3) = -re14*sinr/2.0_ark
                  b0(3,1) = sqrt(3.0_ark)*re14*sinr/2.0_ark
                  b0(3,2) = mX*re14*cosr/(Mtotal+mX)
                  b0(4,3) = -re14*sinr/2.0_ark
                  b0(4,1) = -sqrt(3.0_ark)*re14*sinr/2.0_ark
                  b0(4,2) = mX*re14*cosr/(Mtotal+mX)
                  b0(1,3) = 0.0_ark
                  b0(1,1) = 0.0_ark
                  b0(1,2) = -Mtotal*re14*cosr/(Mtotal+mX)
                  !
                case ('R1-X-R2-Y-Z')
                  !
                  b0(2,1) = re14*sinr
                  b0(2,2) = 0.0_ark
                  b0(2,3) = mX*re14*cosr/(Mtotal+mX)
                  b0(3,1) = -re14*sinr/2.0_ark
                  b0(3,2) = sqrt(3.0_ark)*re14*sinr/2.0_ark
                  b0(3,3) = mX*re14*cosr/(Mtotal+mX)
                  b0(4,1) = -re14*sinr/2.0_ark
                  b0(4,2) = -sqrt(3.0_ark)*re14*sinr/2.0_ark
                  b0(4,3) = mX*re14*cosr/(Mtotal+mX)
                  b0(1,1) = 0.0_ark
                  b0(1,2) = 0.0_ark
                  b0(1,3) = -Mtotal*re14*cosr/(Mtotal+mX)
                  !
                case ('R1-X-R2-Z-Y')
                  !
                  b0(2,1) = re14*sinr
                  b0(2,3) = 0.0_ark
                  b0(2,2) =-mX*re14*cosr/(Mtotal+mX)
                  b0(3,1) = -re14*sinr/2.0_ark
                  b0(3,3) = sqrt(3.0_ark)*re14*sinr/2.0_ark
                  b0(3,2) =-mX*re14*cosr/(Mtotal+mX)
                  b0(4,1) = -re14*sinr/2.0_ark
                  b0(4,3) = -sqrt(3.0_ark)*re14*sinr/2.0_ark
                  b0(4,2) =-mX*re14*cosr/(Mtotal+mX)
                  b0(1,1) = 0.0_ark
                  b0(1,3) = 0.0_ark
                  b0(1,2) = Mtotal*re14*cosr/(Mtotal+mX)
                  !
                end select
                !
              endif
              !
          end select 
          !
          do n = 1,3 
            CM_shift = sum(b0(:,n)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
            b0(:,n) = b0(:,n) - CM_shift
          enddo 
          !      
      end subroutine

  end subroutine ML_b0_XY3

  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  ! Important! The internal angles numbering defers from the conventional 
  ! numbering when each alpha_i is opposite to r_i.
  ! The internal angles are adjusted to the bond of the same numer. 
  ! This is defined by the way we set up zmatrix. 
  ! I.e. internal alpha1,alpha2,alpha3 correspond to alpha3,alpha2,alpha1 (conventional)
  ! alpha3-> 4
  ! alpha2-> 5
  ! alpha1-> 6
  !
  function ML_coordinate_transform_XY3(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct
    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: dsrc(size(src))
    real(ark)                 :: cosalpha,c_t(molec%Nmodes,molec%Nmodes)
    real(ark)                 :: alpha,sinrho,s6,alpha1,alpha2,alpha3,tau_2,r1,r2,r3,norm_2,sindelta
    real(ark)                 :: fact_sign,cosbeta,delta,xt2,xt3,qa,qb,phi3,phi2,phi1,re14
    real(ark)                 :: beta,sinphi,phi12,phi13,phi23,sinbeta1,sinbeta2,sinbeta3,tau,costau
    !
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_XY3/start')") 
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
    nsrc = size(src)
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('MLcoordinate_transform_func: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'MLcoordinate_transform_func - bad coord. type'
       !
    case('C3V-SYMMETRY')
       !
       c_t =  reshape( &
     (/1.0_ark/sqrt(3.0_ark), 1.0_ark/sqrt(3.0_ark), 1.0_ark/sqrt(3.0_ark), 0.0_ark, 0.0_ark, 0.0_ark   , &     ! 1
       0.0_ark      , 0.0_ark  , 0.0_ark , 1.0_ark/sqrt(3.0_ark), 1.0_ark/sqrt(3.0_ark), 1.0_ark/sqrt(3.0_ark)    , &     ! 2
       2.0_ark/sqrt(6.0_ark),-1.0_ark/sqrt(6.0_ark),-1.0_ark/sqrt(6.0_ark), 0.0_ark    , 0.0_ark  , 0.0_ark     , &     ! 3
       0.0_ark   , 1.0_ark/sqrt(2.0_ark) ,-1.0_ark/sqrt(2.0_ark) , 0.0_ark  , 0.0_ark, 0.0_ark , &     ! 4
       0.0_ark  , 0.0_ark   , 0.0_ark , 2.0_ark/sqrt(6.0_ark)   ,-1.0_ark/sqrt(6.0_ark),-1.0_ark/sqrt(6.0_ark), &     ! 5
       0.0_ark  , 0.0_ark   , 0.0_ark , 0.0_ark        , 1.0_ark/sqrt(2.0_ark),-1.0_ark/sqrt(2.0_ark)    /), &  ! 6
       (/6,6/))
       !
       if (direct) then 
           dst = matmul(c_t,dsrc)
       else
           dst = matmul(transpose(c_t),dsrc)
           dst(:) = dst(:) +  molec%local_eq(:)
       endif
       !
    case('NORMAL')
       !
       if (direct) then 
           dst = src
       else
           dst = src
           !dst(:) = dst(:) +  molec%local_eq(:)
       endif
       !
    case('R-S-RHO') 
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(6)-dsrc(5)-dsrc(4) )
          dst(5) = 1.0_ark/sqrt(2.0_ark)*(                 dsrc(5)-dsrc(4) )
          !
          alpha=(src(4)+src(5)+src(6))/3.0_ark
          !
          sinrho = 2.0_ark*sin(alpha*0.5_ark)/sqrt(3.0_ark)
          !
          !rho=pi-dasin(2.0_rk*sin(alpha*0.5_rk)/sqrt(3.0_rk))
          if ( sinrho>1.0_ark+sqrt(small_) ) then 
             !dst(6) = 0.5_ark*pi
             !
             write (out,"('MLcoordinate_transform_func: sinrho>1: ',f18.8)") sinrho
             write (out,"('Consider change difftype ')")
             stop 'MLcoordinate_transform_func - bad sinrho'
             !
          elseif ( sinrho>=1.0_ark) then 
             !
             dst(6) = 0.5_ark*pi
             !
          else 
             if (size(src)==7) then
                if (src(7)<0.0_ark) then 
                   dst(6) = asin(sinrho)
                else
                   dst(6) = pi - asin(sinrho)
                endif
             else
                dst(6) = pi - asin(sinrho)
             endif
          endif 
          !
      else ! not direct
          !
          dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
          !
          sinrho = sin(abs(src(6)))
          !
          alpha = 2.0_ark*asin(sqrt(3.0_ark)*0.5_ark*sinrho)
          !
          s6 = alpha*sqrt(3.0_ark)
          !
          dst(6) =(sqrt(2.0_ark)*s6+2.0_ark*src(4)                     )/sqrt(6.0_ark)
          dst(5) =(sqrt(2.0_ark)*s6-        src(4)+sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
          dst(4) =(sqrt(2.0_ark)*s6-        src(4)-sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
          !
          if (size(dst)==7) then 
             !
             alpha1 = dst(6)
             alpha2 = dst(5)
             alpha3 = dst(4)
             !
             tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
                    +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)
             !
             if ( tau_2<-sqrt(small_) ) then 
                !
                write (out,"('MLcoordinate_transform_func: tau**2<0: ',f18.8)") tau_2
                stop 'MLcoordinate_transform_func - tau**2<0'
                !
             elseif ( tau_2<0.0_ark) then 
                dst(7) = 0
             else 
                dst(7) = sign(1.0_ark,src(6))*sqrt(tau_2)
             endif
             !
          endif 
          !
      endif
       !
       !
    case('R-2D-RHO')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          Qa = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(6)-dsrc(5)-dsrc(4) )
          Qb = 1.0_ark/sqrt(2.0_ark)*(                 dsrc(5)-dsrc(4) )
          !
          dst(4) = (Qa+Qb)/sqrt(2.0_ark)
          dst(5) = (Qa-Qb)/sqrt(2.0_ark)   ! atan2(Qa,Qb)
          !
          alpha=(src(4)+src(5)+src(6))/3.0_ark
          !
          sinrho = 2.0_ark*sin(alpha*0.5_ark)/sqrt(3.0_ark)
          !
          !rho=pi-dasin(2.0_rk*sin(alpha*0.5_rk)/sqrt(3.0_rk))
          if ( sinrho>1.0_ark+sqrt(small_) ) then 
             !dst(6) = 0.5_ark*pi
             !
             write (out,"('MLcoordinate_transform_func: sinrho>1: ',f18.8)") sinrho
             write (out,"('Consider change difftype ')")
             stop 'MLcoordinate_transform_func - bad sinrho'
             !
          elseif ( sinrho>=1.0_ark) then 
             !
             dst(6) = 0.5_ark*pi
             !
          else 
             if (size(src)==7) then
                if (src(7)<0.0_ark) then 
                   dst(6) = asin(sinrho)
                else
                   dst(6) = pi - asin(sinrho)
                endif
             else
                dst(6) = pi - asin(sinrho)
             endif
          endif 
          !
      else ! not direct
          !
          dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
          !
          sinrho = sin(src(6))
          !
          alpha = 2.0_ark*asin(sqrt(3.0_ark)*0.5_ark*sinrho)
          !
          s6 = alpha*sqrt(3.0_ark)
          !
          !Qa = dsrc(4)*cos(dsrc(5))/(cos(dsrc(5))+sin(dsrc(5)))
          !Qb = dsrc(4)*sin(dsrc(5))/(cos(dsrc(5))+sin(dsrc(5)))
          !
          Qa = (dsrc(4)+dsrc(5))/sqrt(2.0_ark)
          Qb = (dsrc(4)-dsrc(5))/sqrt(2.0_ark)
          !
          dst(6) =(sqrt(2.0_ark)*s6+2.0_ark*Qa                 )/sqrt(6.0_ark)
          dst(5) =(sqrt(2.0_ark)*s6-        Qa+sqrt(3.0_ark)*Qb)/sqrt(6.0_ark)
          dst(4) =(sqrt(2.0_ark)*s6-        Qa-sqrt(3.0_ark)*Qb)/sqrt(6.0_ark)
          !
      endif
      !
      case('R-ALPHA')
       !
       if (size(src)/=6) then 
          write(out,"('MLcoordinate_transform_func: r-alpha  works only with 6 coords')")
          stop 'MLcoordinate_transform_func: r-alpha  works only with 6 coords'
       endif 
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          dst(6) = dsrc(4)
          dst(5) = dsrc(5)
          dst(4) = dsrc(6)
          !
      else ! not direct
          !
          dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
          dst(6) = dsrc(4)+molec%local_eq(4)
          dst(5) = dsrc(5)+molec%local_eq(5)
          dst(4) = dsrc(6)+molec%local_eq(6)
          !
      endif
       !
    case('R-S-DELTA')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          if (size(src)==7) then
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
              dst(5) = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )
              !
              dst(6) = src(7)
              !
          elseif (size(src)==6) then
              !
              alpha3 = src(4)
              alpha2 = src(5)
              !
              tau_2 = src(6)**2
              !
              fact_sign = -1.0_ark
              !
              if (alpha2<pi*0.5_rk.and.alpha3<pi*0.5_rk) fact_sign = 1.0_ark
              !     
              cosalpha=cos(alpha2)*cos(alpha3) + fact_sign* &
                       sqrt(cos(alpha2)**2*cos(alpha3)**2+1.0_ark-cos(alpha3)**2-cos(alpha2)**2-tau_2)
              !
              if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLcoordinate_transform_func: cosalpha>1: ',f18.8)") cosalpha
                 stop 'MLcoordinate_transform_func - bad cosalpha'
                 !
              elseif ( cosalpha>=1.0_ark) then 
                 alpha1 = 0.0_ark
              else 
                 alpha1 = acos(cosalpha)
              endif
              !
              dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
              dst(5) = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )
              !
              dst(6) = src(6)
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          !
          if (size(dst)==7) then
             !
             dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             !
             delta = src(6)
             !
             !call find_alpha_from_sindelta_zbrac(src(4),src(5),sin(delta),alpha1,alpha2,alpha3)
             !
             call find_alpha_from_sindelta(src(4),src(5),sin(delta),alpha1,alpha2,alpha3)
             !
             !dst(4:5) = dsrc(4:5)
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
             dst(7) = src(6)
             !
          elseif (size(src)==6) then 
             !
             write(out,"('MLcoordinate_transform_func: inverse for r-s-delta with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for r-s-delta with 6 coords not implemenetd'
             !
          endif 
          !
      endif
      !
    case('R-S-DELTA-MEP')
       !
       if (direct) then
          !
          if (size(src)==7) then
              !
              select case(trim(molec%meptype))
              !
              case('MEP_H3OP') 
                   !
                   tau = src(7) + pi*0.5_ark
                   !
                   re14 = ML_mep_oh3p(tau)
                   !
              end select 
              ! 
              dst(1:3) = src(1:3)-re14
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
              dst(5) = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )
              !
              dst(6) = src(7)
              !
          elseif (size(src)==6) then
              !
              alpha3 = src(4)
              alpha2 = src(5)
              !
              tau_2 = src(6)**2
              !
              fact_sign = -1.0_ark
              !
              if (alpha2<pi*0.5_rk.and.alpha3<pi*0.5_rk) fact_sign = 1.0_ark
              !     
              cosalpha=cos(alpha2)*cos(alpha3) + fact_sign* &
                       sqrt(cos(alpha2)**2*cos(alpha3)**2+1.0_ark-cos(alpha3)**2-cos(alpha2)**2-tau_2)
              !
              if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLcoordinate_transform_func: cosalpha>1: ',f18.8)") cosalpha
                 stop 'MLcoordinate_transform_func - bad cosalpha'
                 !
              elseif ( cosalpha>=1.0_ark) then 
                 alpha1 = 0.0_ark
              else 
                 alpha1 = acos(cosalpha)
              endif
              !
              dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
              dst(5) = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )
              !
              dst(6) = src(6)
              !
              select case(trim(molec%meptype))
              !
              case('MEP_H3OP') 
                   !
                   tau = src(6) + pi*0.5_ark
                   !
                   re14 = ML_mep_oh3p(tau)
                   !
              end select 
              ! 
              dst(1:3) = src(1:3)-re14
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          if (size(dst)==7) then
             !
             delta = src(6)
             !
             !call find_alpha_from_sindelta_zbrac(src(4),src(5),sin(delta),alpha1,alpha2,alpha3)
             !
             call find_alpha_from_sindelta(src(4),src(5),sin(delta),alpha1,alpha2,alpha3)
             !
             !dst(4:5) = dsrc(4:5)
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
             dst(7) = src(6)
             !
             select case(trim(molec%meptype))
             !
             case('MEP_H3OP') 
                  !
                  tau = delta + pi*0.5_ark
                  !
                  re14 = ML_mep_oh3p(tau)
                  !
             end select 
             ! 
             dst(1:3) = src(1:3)+re14
             !
          elseif (size(src)==6) then 
             !
             write(out,"('MLcoordinate_transform_func: inverse for r-s-delta with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for r-s-delta with 6 coords not implemenetd'
             !
          endif 
          !
      endif
       !
    case('R-PHI-DELTA')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          if (size(src)==7) then
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              delta = src(7)
              !
              dst(6) = src(7)
              !
              sinphi = sin(alpha1*0.5_ark)/cos(delta)
              phi1 = asin(sinphi)*2.0_ark
              !
              sinphi = sin(alpha2*0.5_ark)/cos(delta)
              phi2 = asin(sinphi)*2.0_ark
              !
              sinphi = sin(alpha3*0.5_ark)/cos(delta)
              phi3 = asin(sinphi)*2.0_ark
              !
              dst(4) = phi2-2.0_ark/3.0_ark*pi
              dst(5) = phi3-2.0_ark/3.0_ark*pi
              !
              dst(6) = src(7)
              !
              beta = pi*0.5_ark+dst(6)
              !
              ! check if beta is correct
              !
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi1 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha1)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi1, alpha1 dont agree ',3f12.6)") beta,phi1,alpha1
                stop 'xy3-coords, beta, phi1, alpha1 dont agree'
                !
              endif 
              !
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi2 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha2)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi2, alpha2 dont agree ',3f12.6)") beta,phi2,alpha2
                stop 'xy3-coords, beta, phi2, alpha2 dont agree'
                !
              endif 
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi3 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha3)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi3, alpha3 dont agree ',3f12.6)") beta,phi3,alpha3
                stop 'xy3-coords, beta, phi3, alpha3 dont agree'
                !
              endif 
              !
              continue
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          if (size(dst)==7) then
             !
             dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             !
             delta = src(6)
             !
             phi2 = src(4)+2.0_ark/3.0_ark*pi
             !
             phi3 = src(5)+2.0_ark/3.0_ark*pi
             !
             phi1 = 2.0_ark*pi-phi2-phi3
             !
             alpha1 = asin(sin(phi1*0.5_ark)*cos(delta))*2.0_ark
             alpha2 = asin(sin(phi2*0.5_ark)*cos(delta))*2.0_ark
             alpha3 = asin(sin(phi3*0.5_ark)*cos(delta))*2.0_ark
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
             dst(7) = src(6)
             !
          else
             !
             write(out,"('MLcoordinate_transform_func: inverse for R-SYMPHI-DELTA with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for R-SYMPHI-DELTA with 6 coords not implemenetd'
             !
          endif 
          !
      endif
       !
    case('R-PHI-DELTA-MEP')
       !
       if (direct) then
          !
          if (size(src)==7) then
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              delta = src(7)
              !
              select case(trim(molec%meptype))
              !
              case('MEP_H3OP') 
                   !
                   tau = src(7) + pi*0.5_ark
                   !
                   re14 = ML_mep_oh3p(tau)
                   !
              end select 
              ! 
              dst(1:3) = src(1:3)-re14
              dst(6) = src(7)
              !
              sinphi = sin(alpha1*0.5_ark)/cos(delta)
              phi1 = asin(sinphi)*2.0_ark
              !
              sinphi = sin(alpha2*0.5_ark)/cos(delta)
              phi2 = asin(sinphi)*2.0_ark
              !
              sinphi = sin(alpha3*0.5_ark)/cos(delta)
              phi3 = asin(sinphi)*2.0_ark
              !
              !phi1 = acos( (cos(alpha1)-sin(delta)**2)/cos(delta)**2 )
              !phi2 = acos( (cos(alpha2)-sin(delta)**2)/cos(delta)**2 )
              !phi3 = acos( (cos(alpha3)-sin(delta)**2)/cos(delta)**2 )
              !
              dst(4) = phi2-2.0_ark/3.0_ark*pi
              dst(5) = phi3-2.0_ark/3.0_ark*pi
              !
              beta = pi*0.5_ark+dst(6)
              !
              ! check if beta is correct
              !
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi1 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha1)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi1, alpha1 dont agree ',3f12.6)") beta,phi1,alpha1
                stop 'xy3-coords, beta, phi1, alpha1 dont agree'
                !
              endif 
              !
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi2 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha2)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi2, alpha2 dont agree ',3f12.6)") beta,phi2,alpha2
                stop 'xy3-coords, beta, phi2, alpha2 dont agree'
                !
              endif 
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi3 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha3)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi3, alpha3 dont agree ',3f12.6)") beta,phi3,alpha3
                stop 'xy3-coords, beta, phi3, alpha3 dont agree'
                !
              endif 
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          if (size(dst)==7) then
             !
             delta = src(6)
             !
             phi2 = src(4)+2.0_ark/3.0_ark*pi
             !
             phi3 = src(5)+2.0_ark/3.0_ark*pi
             !
             phi1 = 2.0_ark*pi-phi2-phi3
             !
             alpha1 = asin(sin(phi1*0.5_ark)*cos(delta))*2.0_ark
             alpha2 = asin(sin(phi2*0.5_ark)*cos(delta))*2.0_ark
             alpha3 = asin(sin(phi3*0.5_ark)*cos(delta))*2.0_ark
             !
             !alpha1 = acos( cos(phi1)*cos(delta)**2+sin(delta)**2 )
             !alpha2 = acos( cos(phi2)*cos(delta)**2+sin(delta)**2 )
             !alpha3 = acos( cos(phi3)*cos(delta)**2+sin(delta)**2 )
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
             dst(7) = src(6)
             !
             select case(trim(molec%meptype))
             !
             case('MEP_H3OP') 
                  !
                  tau = delta + pi*0.5_ark
                  !
                  re14 = ML_mep_oh3p(tau)
                  !
             end select 
             ! 
             dst(1:3) = src(1:3)+re14
             !
          else
             !
             write(out,"('MLcoordinate_transform_func: inverse for R-SYMPHI-DELTA with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for R-SYMPHI-DELTA with 6 coords not implemenetd'
             !
          endif 
          !
      endif
       !
    case('R-SYMPHI-DELTA')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          if (size(src)==7) then
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              dst(6) = src(7)
              !
              delta = src(7)
              !
              sinphi = sin(alpha1*0.5_ark)/cos(delta)
              phi1 = asin(sinphi)*2.0_ark
              !
              sinphi = sin(alpha2*0.5_ark)/cos(delta)
              phi2 = asin(sinphi)*2.0_ark
              !
              sinphi = sin(alpha3*0.5_ark)/cos(delta)
              phi3 = asin(sinphi)*2.0_ark
              !
              dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
              dst(5) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
              !
              dst(6) = delta
              !
              beta = pi*0.5_ark+delta
              !
              ! check if beta is correct
              !
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi1 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha1)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi1, alpha1 dont agree ',3f12.6)") beta,phi1,alpha1
                stop 'xy3-coords, beta, phi1, alpha1 dont agree'
                !
              endif 
              !
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi2 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha2)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi2, alpha2 dont agree ',3f12.6)") beta,phi2,alpha2
                stop 'xy3-coords, beta, phi2, alpha2 dont agree'
                !
              endif 
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi3 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha3)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi3, alpha3 dont agree ',3f12.6)") beta,phi3,alpha3
                stop 'xy3-coords, beta, phi3, alpha3 dont agree'
                !
              endif 
              !
              continue
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          if (size(dst)==7) then
             !
             dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             !
             delta = src(6)
             dst(7) = delta
             !
             phi1 = (sqrt(6.0_ark)*src(4)+2.0_ark*pi)/3.0_ark
             !
             phi2 = (sqrt(2.0_ark)*src(5)+2.0_ark*pi-phi1)/2.0_ark
             !
             phi3 = 2.0_ark*pi-phi1-phi2
             !
             alpha1 = asin(sin(phi1*0.5_ark)*cos(delta))*2.0_ark
             alpha2 = asin(sin(phi2*0.5_ark)*cos(delta))*2.0_ark
             alpha3 = asin(sin(phi3*0.5_ark)*cos(delta))*2.0_ark
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
          else
             !
             write(out,"('MLcoordinate_transform_func: inverse for R-SYMPHI-DELTA with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for R-SYMPHI-DELTA with 6 coords not implemenetd'
             !
          endif 
          !
      endif
       !
    case('R-SYMPHI-DELTA-MEP')
       !
       if (direct) then
          !
          if (size(src)==7) then
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              dst(6) = src(7)
              !
              delta = src(7)
              !
              select case(trim(molec%meptype))
              !
              case('MEP_H3OP') 
                   !
                   tau = src(7) + pi*0.5_ark
                   !
                   re14 = ML_mep_oh3p(tau)
                   !
              end select
              !
              dst(1:3) = src(1:3)-re14
              !
              sinphi = sin(alpha1*0.5_ark)/cos(delta)
              phi1 = asin(sinphi)*2.0_ark
              !
              sinphi = sin(alpha2*0.5_ark)/cos(delta)
              phi2 = asin(sinphi)*2.0_ark
              !
              sinphi = sin(alpha3*0.5_ark)/cos(delta)
              phi3 = asin(sinphi)*2.0_ark
              !
              dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
              dst(5) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
              !
              dst(6) = delta
              !
              beta = pi*0.5_ark+delta
              !
              ! check if beta is correct
              !
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi1 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha1)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi1, alpha1 dont agree ',3f12.6)") beta,phi1,alpha1
                stop 'xy3-coords, beta, phi1, alpha1 dont agree'
                !
              endif 
              !
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi2 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha2)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi2, alpha2 dont agree ',3f12.6)") beta,phi2,alpha2
                stop 'xy3-coords, beta, phi2, alpha2 dont agree'
                !
              endif 
              cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi3 )
              alpha = acos(cosalpha)
              if (abs(alpha-alpha3)>sqrt(small_)) then
                !
                write(out,"('xy3-coords, beta, phi3, alpha3 dont agree ',3f12.6)") beta,phi3,alpha3
                stop 'xy3-coords, beta, phi3, alpha3 dont agree'
                !
              endif 
              !
              continue
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          if (size(dst)==7) then
             !
             delta = src(6)
             dst(7) = delta
             !
             select case(trim(molec%meptype))
             !
             case('MEP_H3OP') 
                  !
                  tau = delta + pi*0.5_ark
                  !
                  re14 = ML_mep_oh3p(tau)
                  !
             end select 
             !
             dst(1:3) = src(1:3)+re14
             !
             phi1 = (sqrt(6.0_ark)*src(4)+2.0_ark*pi)/3.0_ark
             !
             phi2 = (sqrt(2.0_ark)*src(5)+2.0_ark*pi-phi1)/2.0_ark
             !
             phi3 = 2.0_ark*pi-phi1-phi2
             !
             alpha1 = asin(sin(phi1*0.5_ark)*cos(delta))*2.0_ark
             alpha2 = asin(sin(phi2*0.5_ark)*cos(delta))*2.0_ark
             alpha3 = asin(sin(phi3*0.5_ark)*cos(delta))*2.0_ark
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
          else
             !
             write(out,"('MLcoordinate_transform_func: inverse for R-SYMPHI-DELTA with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for R-SYMPHI-DELTA with 6 coords not implemenetd'
             !
          endif 
          !
      endif
       !
    case('R-A2A3-DELTA')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          if (size(src)==7) then
              !
              delta = src(7)
              !
              dst(5) = dsrc(4)
              dst(4) = dsrc(5)
              dst(6) = src(7)
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          if (size(dst)==7) then
             !
             dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             !
             delta = src(6)
             !
             alpha2 = src(4)+molec%local_eq(5)
             alpha3 = src(5)+molec%local_eq(4)
             !
             sinphi = sin(alpha2*0.5_ark)/cos(delta)
             phi2 = asin(sinphi)*2.0_ark
             !
             sinphi = sin(alpha3*0.5_ark)/cos(delta)
             phi3 = asin(sinphi)*2.0_ark
             !
             phi1 = 2.0_ark*pi-phi2-phi3
             !
             alpha1 = asin(sin(phi1*0.5_ark)*cos(delta))*2.0_ark
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
             dst(7) = src(6)
             !
          else
             !
             write(out,"('MLcoordinate_transform_func: inverse for R-A2A3-DELTA with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for R-A2A3-DELTA with 6 coords not implemenetd'
             !
          endif 
          !
      endif

       !
    case('R-SYMPHI-TAU')
       !
       if (molec%ncoords/=6) then
         !
         write(out,"('MLcoordinate_transform_func: inverse for R-SYMPHI-TAU with ',i6,' coords not implemenetd')") molec%ncoords
         stop 'MLcoordinate_transform_func: inverse for R-SYMPHI-TAU with 6 coords not implemenetd'
         !
       endif 
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          !
          alpha3 = src(4)
          alpha2 = src(5)
          tau    = src(6)
          !
          cosalpha = cos(alpha2)*cos(alpha3) + sin(alpha2)*sin(alpha3)*cos(tau)
          !
          if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('MLcoordinate_transform_func: cosalpha>1: ',f18.8)") cosalpha
             stop 'MLcoordinate_transform_func - bad cosalpha'
             !
          endif 
          !
          alpha1 = acos(cosalpha)
          !
          costau =  (cos(alpha1)-cos(alpha2)*cos(alpha3))/(sin(alpha2)*sin(alpha3) )
          !
          tau = acos(costau)
          !
          !dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
          !dst(5) = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )

          ! tau
          !dst(6) = src(7)
          !
          tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
                 +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)
          !
          r1 = src(1)
          r2 = src(2)
          r3 = src(3)
          !
          norm_2 = sin(alpha3)**2+sin(alpha2)**2+sin(alpha1)**2+&
                   2.0_ark*cos(alpha3)*cos(alpha1)-2.0_ark*cos(alpha2)+ &
                   2.0_ark*cos(alpha2)*cos(alpha3)-2.0_ark*cos(alpha1)+ &
                   2.0_ark*cos(alpha2)*cos(alpha1)-2.0_ark*cos(alpha3)
                   !
          sindelta = sqrt(tau_2)/sqrt(norm_2)
          !
          if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('MLcoordinate_transform_func: sindelta>1: ',f18.8)") sindelta
             write (out,"('Consider change difftype ')")
             stop 'MLcoordinate_transform_func - bad sindelta'
             !
          elseif ( sindelta>=1.0_ark) then 
             !
             delta = 0.0_ark
             !
          else 
             !
             delta = asin(sindelta)
             !
             if (src(6)>pi) delta = -delta
             !
          endif 
          !
          sinphi = sin(alpha1*0.5_ark)/cos(delta)
          phi1 = asin(sinphi)*2.0_ark
          !
          sinphi = sin(alpha2*0.5_ark)/cos(delta)
          phi2 = asin(sinphi)*2.0_ark
          !
          sinphi = sin(alpha3*0.5_ark)/cos(delta)
          phi3 = asin(sinphi)*2.0_ark
          !
          dst(4) = phi2  ! 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
          dst(5) = phi3  ! 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
          !
          dst(6) = delta
          !
          beta = pi*0.5_ark+delta
          !
          ! check if beta is correct
          !
          cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi1 )
          alpha = acos(cosalpha)
          if (abs(alpha-alpha1)>sqrt(small_)) then
            !
            write(out,"('xy3-coords, beta, phi1, alpha1 dont agree ',3f12.6)") beta,phi1,alpha1
            stop 'xy3-coords, beta, phi1, alpha1 dont agree'
            !
          endif 
          !
          cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi2 )
          alpha = acos(cosalpha)
          if (abs(alpha-alpha2)>sqrt(small_)) then
            !
            write(out,"('xy3-coords, beta, phi2, alpha2 dont agree ',3f12.6)") beta,phi2,alpha2
            stop 'xy3-coords, beta, phi2, alpha2 dont agree'
            !
          endif 
          cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi3 )
          alpha = acos(cosalpha)
          if (abs(alpha-alpha3)>sqrt(small_)) then
            !
            write(out,"('xy3-coords, beta, phi3, alpha3 dont agree ',3f12.6)") beta,phi3,alpha3
            stop 'xy3-coords, beta, phi3, alpha3 dont agree'
            !
          endif 
          !
          continue
          !
          !
      else    ! not direct
          !
          !
          dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
          !
          delta = src(6)
          !
          phi2 = src(4)  ! (sqrt(6.0_ark)*src(4)+2.0_ark*pi)/3.0_ark
          !
          phi3 = src(5)  ! (sqrt(2.0_ark)*src(5)+2.0_ark*pi-phi1)/2.0_ark
          !
          phi1 = 2.0_ark*pi-phi2-phi3
          !
          beta = pi*0.5_ark-delta
          !
          cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi1 )
          alpha1 = acos(cosalpha)
          !
          cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi2 )
          alpha2 = acos(cosalpha)
          !
          cosalpha = cos( beta )**2 + sin( beta )**2*cos( phi3 )
          alpha3 = acos(cosalpha)
          !
          !alpha1 = asin(sin(phi1*0.5_ark)*cos(delta))*2.0_ark
          !alpha2 = asin(sin(phi2*0.5_ark)*cos(delta))*2.0_ark
          !alpha3 = asin(sin(phi3*0.5_ark)*cos(delta))*2.0_ark
          !
          dst(4) = alpha3
          dst(5) = alpha2
          !
          costau =  (cos(alpha1)-cos(alpha2)*cos(alpha3))/(sin(alpha2)*sin(alpha3) )
          !
          tau = acos(costau)
          !
          if (delta<0) tau = 2.0_ark*pi-tau
          !
          dst(6) = tau
          !
      endif
       !
    case('SYM-DELTA')
       !
       if (direct) then
          ! 
          !dst(1:3) = dsrc(1:3)

          dst(1) = 1.0_ark/sqrt(3.0_ark)*(         dsrc(1)+dsrc(2)+dsrc(3) )
          dst(2) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(1)-dsrc(2)-dsrc(3) )
          dst(3) = 1.0_ark/sqrt(2.0_ark)*(                 dsrc(2)-dsrc(3) )


          if (size(src)==7) then
              !
              !dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(4)-dsrc(5)-dsrc(6) )
              !dst(5) = 1.0_ark/sqrt(2.0_ark)*(                 dsrc(5)-dsrc(6) )
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
              dst(5) = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )

              ! tau
              dst(6) = src(7)
              !
              !tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
              !       +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)
              !r1 = src(1)
              !r2 = src(2)
              !r3 = src(3)
              !
              !norm_2 = sin(alpha3)**2+sin(alpha2)**2+sin(alpha1)**2+&
              !         2.0_ark*cos(alpha3)*cos(alpha1)-2.0_ark*cos(alpha2)+ &
              !         2.0_ark*cos(alpha2)*cos(alpha3)-2.0_ark*cos(alpha1)+ &
              !         2.0_ark*cos(alpha2)*cos(alpha1)-2.0_ark*cos(alpha3)
              !sindelta = src(7)/sqrt(norm_2)
              !
              !if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
              !   !
              !   write (out,"('MLcoordinate_transform_func: sindelta>1: ',f18.8)") sindelta
              !   write (out,"('Consider change difftype ')")
              !   stop 'MLcoordinate_transform_func - bad sindelta'
              !   !
              !elseif ( sindelta>=1.0_ark) then 
              !   !
              !   dst(6) = 0.0_ark
              !   !
              !else 
              !   dst(6) = asin(sindelta)
              !   !
              !endif 
              !
              !if (src(7)<0) then 
              !  dst(2)=dsrc(3)
              !  dst(3)=dsrc(2)
              !  dst(5)=-dst(5)
              !endif
              !
          elseif (size(src)==6) then
              !
              write(out,"('MLcoordinate_transform_func:  for r-s-delta with 6 coords not implemenetd')")
              stop 'MLcoordinate_transform_func: for r-s-delta with 6 coords not implemenetd'
              !
              alpha3 = src(4)
              alpha2 = src(5)
              !
              tau_2 = src(6)**2
              !
              fact_sign = -1.0_ark
              !
              if (alpha2<pi*0.5_rk.and.alpha3<pi*0.5_rk) fact_sign = 1.0_ark
              !     
              cosalpha=cos(alpha2)*cos(alpha3) + fact_sign* &
                       sqrt(cos(alpha2)**2*cos(alpha3)**2+1.0_ark-cos(alpha3)**2-cos(alpha2)**2-tau_2)
              !
              if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLcoordinate_transform_func: cosalpha>1: ',f18.8)") cosalpha
                 stop 'MLcoordinate_transform_func - bad cosalpha'
                 !
              elseif ( cosalpha>=1.0_ark) then 
                 alpha1 = 0.0_ark
              else 
                 alpha1 = acos(cosalpha)
              endif
              !
              dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
              dst(5) = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )
              !
              norm_2 = sin(alpha3)**2+sin(alpha2)**2+sin(alpha1)**2+&
                       2._ark*cos(alpha3)*cos(alpha1)-2._ark*cos(alpha2)+ &
                       2._ark*cos(alpha2)*cos(alpha3)-2._ark*cos(alpha1)+ &
                       2._ark*cos(alpha2)*cos(alpha1)-2._ark*cos(alpha3)
              sindelta = src(6)/sqrt(norm_2)
              !
              if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLcoordinate_transform_func: sindelta>1: ',f18.8)") sindelta
                 write (out,"('Consider change difftype ')")
                 stop 'MLcoordinate_transform_func - bad sindelta'
                 !
              elseif ( sindelta>=1.0_ark) then 
                 !
                 dst(6) = 0.0_ark
                 !
              else 
                 dst(6) = asin(sindelta)
                 !
              endif 
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          dst(1) =(sqrt(2.0_ark)*src(1)+2.0_ark*src(2)                     )/sqrt(6.0_ark)
          dst(2) =(sqrt(2.0_ark)*src(1)-        src(2)+sqrt(3.0_ark)*src(3))/sqrt(6.0_ark)
          dst(3) =(sqrt(2.0_ark)*src(1)-        src(2)-sqrt(3.0_ark)*src(3))/sqrt(6.0_ark)
          !
          !dst(1:3) = dst(1:3)+molec%local_eq(1:3)
          !
          if (size(dst)==7) then
             !
             delta = src(6)
             !
             call find_alpha_from_sindelta(src(4),src(5),sin(delta),alpha1,alpha2,alpha3)
             !
             !dst(4:5) = dsrc(4:5)
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
             dst(7) = src(6)
             !
             !dst(4) =(sqrt(2.0_ark)*s6+2.0_ark*src(4)                     )/sqrt(6.0_ark)
             !dst(6) =(sqrt(2.0_ark)*s6-        src(4)+sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
             !dst(5) =(sqrt(2.0_ark)*s6-        src(4)-sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
             !
             !tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
             !         +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)
             !
             !if ( tau_2<-sqrt(small_) ) then 
             !   !
             !   write (out,"('MLcoordinate_transform_func: tau**2<0: ',f18.8)") tau_2
             !   stop 'MLcoordinate_transform_func - tau**2<0'
             !   !
             !elseif ( tau_2<0.0_ark) then 
             !   dst(7) = 0
             !else 
             !   dst(7) = sign(1.0_ark,delta)*sqrt(tau_2)
             !endif
             !
             !
             !s6 = alpha1*sqrt(3.0_ark)-sqrt(2.0_ark)*src(4)
             !
             !
          elseif (size(src)==6) then 
             !
             !dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             write(out,"('MLcoordinate_transform_func: inverse for sym-delta with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for sym-delta with 6 coords not implemenetd'
             !dst(6) = dsrc(6)+molec%local_eq(6)
             !dst(4:5) = dsrc(4:5)
             !
             !call find_alpha_from_tau(src(4),src(5),sin(src(6)),alpha1)
             !
             !s6 = alpha1*sqrt(3.0_ark)-sqrt(2.0_ark)*src(4)
             !
             !dst(4) =(sqrt(2.0_ark)*s6+2.0_ark*src(4)                     )/sqrt(6.0_ark)
             !dst(5) =(sqrt(2.0_ark)*s6-        src(4)+sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
             !dst(4) =(sqrt(2.0_ark)*s6-        src(4)-sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
             !
          endif 
          !
      endif
      !       !
    case('R-2D-DELTA')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)

          !dst(1) = 1.0_ark/sqrt(3.0_ark)*(         dsrc(1)+dsrc(2)+dsrc(3) )
          !dst(2) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(1)-dsrc(2)-dsrc(3) )
          !dst(3) = 1.0_ark/sqrt(2.0_ark)*(                 dsrc(2)-dsrc(3) )

          if (size(src)==7) then
              !
              !dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(4)-dsrc(5)-dsrc(6) )
              !dst(5) = 1.0_ark/sqrt(2.0_ark)*(                 dsrc(5)-dsrc(6) )
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              Qa = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
              Qb = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )
              !
              dst(4) = (Qa+Qb)/sqrt(2.0_ark)
              dst(5) = (Qa-Qb)/sqrt(2.0_ark)   ! atan2(Qa,Qb)
              !
              ! tau
              !dst(6) = src(7)
              !
              tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
                     +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)
              r1 = src(1)
              r2 = src(2)
              r3 = src(3)
              !
              norm_2 = sin(alpha3)**2+sin(alpha2)**2+sin(alpha1)**2+&
                       2.0_ark*cos(alpha3)*cos(alpha1)-2.0_ark*cos(alpha2)+ &
                       2.0_ark*cos(alpha2)*cos(alpha3)-2.0_ark*cos(alpha1)+ &
                       2.0_ark*cos(alpha2)*cos(alpha1)-2.0_ark*cos(alpha3)
              sindelta = src(7)/sqrt(norm_2)
              !
              if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLcoordinate_transform_func: sindelta>1: ',f18.8)") sindelta
                 write (out,"('Consider change difftype ')")
                 stop 'MLcoordinate_transform_func - bad sindelta'
                 !
              elseif ( sindelta>=1.0_ark) then 
                 !
                 dst(6) = 0.0_ark
                 !
              else 
                 dst(6) = asin(sindelta)
                 !
              endif 
              !
              !if (src(7)<0) then 
              !  dst(2)=dsrc(3)
              !  dst(3)=dsrc(2)
              !  dst(5)=-dst(5)
              !endif
              !
          elseif (size(src)==6) then
              !
              write(out,"('MLcoordinate_transform_func: inverse for r-s-delta with 6 coords not implemenetd')")
              stop 'MLcoordinate_transform_func: inverse for r-s-delta with 6 coords not implemenetd'
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          !
          if (size(dst)==7) then
             !
             dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             !
             delta = src(6)
             !
             Qa = (dsrc(4)+dsrc(5))/sqrt(2.0_ark)
             Qb = (dsrc(4)-dsrc(5))/sqrt(2.0_ark)
             !
             call find_alpha_from_sindelta(Qa,Qb,sin(delta),alpha1,alpha2,alpha3)
             !
             !dst(4:5) = dsrc(4:5)
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
             dst(7) = src(6)
             !
             !dst(4) =(sqrt(2.0_ark)*s6+2.0_ark*src(4)                     )/sqrt(6.0_ark)
             !dst(6) =(sqrt(2.0_ark)*s6-        src(4)+sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
             !dst(5) =(sqrt(2.0_ark)*s6-        src(4)-sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
             !
             !tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
             !         +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)
             !
             !if ( tau_2<-sqrt(small_) ) then 
             !   !
             !   write (out,"('MLcoordinate_transform_func: tau**2<0: ',f18.8)") tau_2
             !   stop 'MLcoordinate_transform_func - tau**2<0'
             !   !
             !elseif ( tau_2<0.0_ark) then 
             !   dst(7) = 0
             !else 
             !   dst(7) = sign(1.0_ark,delta)*sqrt(tau_2)
             !endif
             !
             !
             !s6 = alpha1*sqrt(3.0_ark)-sqrt(2.0_ark)*src(4)
             !
             !
          elseif (size(src)==6) then 
             !
             !dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             write(out,"('MLcoordinate_transform_func: inverse for r-2d-delta with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for r-2d-delta with 6 coords not implemenetd'
             !dst(6) = dsrc(6)+molec%local_eq(6)
             !dst(4:5) = dsrc(4:5)
             !
             !call find_alpha_from_tau(src(4),src(5),sin(src(6)),alpha1)
             !
             !s6 = alpha1*sqrt(3.0_ark)-sqrt(2.0_ark)*src(4)
             !
             !dst(4) =(sqrt(2.0_ark)*s6+2.0_ark*src(4)                     )/sqrt(6.0_ark)
             !dst(5) =(sqrt(2.0_ark)*s6-        src(4)+sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
             !dst(4) =(sqrt(2.0_ark)*s6-        src(4)-sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
             !
          endif 
          !
      endif
       !
    case('R-S4-RHO-PHI-DELTA')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          if (size(src)==7) then
              !
              alpha3 = src(4)
              alpha2 = src(5)
              alpha1 = src(6)
              !
              Qa = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha1-alpha2-alpha3 )
              Qb = 1.0_ark/sqrt(2.0_ark)*(                alpha2-alpha3 )
              !
              dst(4) = (Qa+Qb)*0.5_ark
              dst(5) = atan2(Qa,Qb)
              !
              tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
                     +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)
              r1 = src(1)
              r2 = src(2)
              r3 = src(3)
              !
              norm_2 = sin(alpha3)**2+sin(alpha2)**2+sin(alpha1)**2+&
                       2.0_ark*cos(alpha3)*cos(alpha1)-2.0_ark*cos(alpha2)+ &
                       2.0_ark*cos(alpha2)*cos(alpha3)-2.0_ark*cos(alpha1)+ &
                       2.0_ark*cos(alpha2)*cos(alpha1)-2.0_ark*cos(alpha3)
              sindelta = src(7)/sqrt(norm_2)
              !
              if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLcoordinate_transform_func: sindelta>1: ',f18.8)") sindelta
                 write (out,"('Consider change difftype ')")
                 stop 'MLcoordinate_transform_func - bad sindelta'
                 !
              elseif ( sindelta>=1.0_ark) then 
                 !
                 dst(6) = 0.0_ark
                 !
              else 
                 dst(6) = asin(sindelta)
                 !
              endif 
              !
              !
          elseif (size(src)==6) then
              !
              write(out,"('MLcoordinate_transform_func: inverse for R-S4-RHO-PHI-DELTA with 6 coords not implemenetd')")
              stop 'MLcoordinate_transform_func: inverse for R-S4-RHO-PHI-DELTA with 6 coords not implemenetd'
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          !
          if (size(dst)==7) then
             !
             dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             !
             delta = src(6)
             !
             Qa = 2.0_ark*src(4)/(1.0_ark+tan(src(5)) )
             Qb = 2.0_ark*tan(src(5))*src(4)/(1.0_ark+tan(src(5)) )
             !
             call find_alpha_from_sindelta(Qa,Qb,sin(delta),alpha1,alpha2,alpha3)
             !
             !dst(4:5) = dsrc(4:5)
             !
             dst(6) = alpha1
             dst(5) = alpha2
             dst(4) = alpha3
             !
             dst(7) = src(6)
             !
          elseif (size(src)==6) then 
             !
             write(out,"('MLcoordinate_transform_func: inverse for r-2d-delta with 6 coords not implemenetd')")
             stop 'MLcoordinate_transform_func: inverse for r-2d-delta with 6 coords not implemenetd'
             !
          endif 
          !
      endif
          !
    case('R-A2-A3-TAU')
       !
       if (direct) then
          ! 
          dst(1:5) = dsrc(1:5)
          !
          dst(6) = src(6)
          !
      else ! not direct
          !
          dst(1:5) = dsrc(1:5)+molec%local_eq(1:5)
          !
          dst(6) = src(6)
          !
      endif
          !
       !
       !
    case('R-B1-B2-BETA')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          dst(4:6) = dsrc(4:6)
          !
          !dst(1) = 1.0_ark/sqrt(3.0_ark)*(         dsrc(1)+dsrc(2)+dsrc(3) )
          !dst(2) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(1)-dsrc(2)-dsrc(3) )
          !dst(3) = 1.0_ark/sqrt(2.0_ark)*(                 dsrc(2)-dsrc(3) )
          !
          alpha3 = src(4)
          alpha2 = src(5)
          alpha1 = src(6)
          !
          xt2 = (sin(alpha1*0.5_ark)**2+sin(alpha3*0.5_ark)**2-sin(alpha2*0.5_ark)**2)/ &
                ( 2.0_ark*sin(alpha3*0.5_ark)*sin(alpha1*0.5_ark) )
          xt3 = (sin(alpha2*0.5_ark)**2+sin(alpha1*0.5_ark)**2-sin(alpha3*0.5_ark)**2)/ &
                ( 2.0_ark*sin(alpha1*0.5_ark)*sin(alpha2*0.5_ark) )
          !
          dst(4) = 2.0_ark*acos(xt2)-2.0_ark*pi/3.0_ark
          dst(5) = 2.0_ark*acos(xt3)-2.0_ark*pi/3.0_ark
          !
          !cosbeta = src(7)/sin(alpha1)
          !if ( abs(cosbeta)>1.0_ark+sqrt(small_) ) then 
          !   !
          !   ierror = 1
          !   return
          !   !
          !elseif ( cosbeta>=1.0_ark) then 
          !   dst(6) = 0.0_ark
          !else 
          !   dst(6) = acos(cosbeta)
          !endif 
          !
          cosbeta = src(7)/&
                sqrt(-cos(alpha3)**2+2.0_ark*cos(alpha3)*cos(alpha2)-2.0_ark*cos(alpha3)+ &
                2.0_ark*cos(alpha3)*cos(alpha1)-2.0_ark*cos(alpha1)+2*cos(alpha1)*cos(alpha2)-cos(alpha2)**2-&
                cos(alpha1)**2-2.0_ark*cos(alpha2)+3.0_ark)
          !
          if ( abs(cosbeta)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('MLcoordinate_transform_func: cosbeta>1: ',f18.8)") cosbeta
             write (out,"('Consider change difftype ')")
             stop 'MLcoordinate_transform_func - bad cosbeta'
             !
          elseif ( cosbeta>=1.0_ark) then 
             !
             dst(6) = 0.0_ark
             !
          else 
             dst(6) = acos(cosbeta)
          endif 
          !
      else ! not direct
          !
          write (out,"('MLcoordinate_transform_func: r-b1-b2-beta inverse coord. transformation is not supported.')") 
          stop 'MLcoordinate_transform_func - bad inverse transformation'
          !
      endif
          !
    case('R-S-LINRHO')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(4)-dsrc(5)-dsrc(6) )
          dst(5) = 1.0_ark/sqrt(2.0_ark)*(               dsrc(5)-dsrc(6) )
          !
          alpha=(src(4)+src(5)+src(6))/3.0_ark
          !
          sinrho = 2.0_ark*sin(alpha*0.5_ark)/sqrt(3.0_ark)
          !
          !rho=pi-dasin(2.0_rk*sin(alpha*0.5_rk)/sqrt(3.0_rk))
          if ( sinrho>=1.0_ark ) then 
             dst(6) = 0.5_ark*pi-molec%local_eq(6)
          else 
             dst(6) = pi - asin(sinrho)-molec%local_eq(6)
          endif 
          !
      else ! not direct
          !
          dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
          !
          sinrho = sin(src(6)+molec%local_eq(6))
          !
          alpha = 2.0_ark*asin(sqrt(3.0_ark)*0.5_ark*sinrho)
          !
          s6 = alpha*sqrt(3.0_ark)
          !
          dst(4) =(sqrt(2.0_ark)*s6+2.0_ark*src(4)                    )/sqrt(6.0_ark)
          dst(5) =(sqrt(2.0_ark)*s6-        src(4)+sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
          dst(6) =(sqrt(2.0_ark)*s6-        src(4)-sqrt(3.0_ark)*src(5))/sqrt(6.0_ark)
          !
      endif
      !
    end select
    !
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_XY3/end')") 
    !
    !
  end function ML_coordinate_transform_XY3
  !
  subroutine find_alpha_from_sindelta(S4,S5,sindelta,alpha1,alpha2,alpha3)

    real(ark),intent(in)  :: S4,S5,sindelta
    real(ark),intent(out) :: alpha1,alpha2,alpha3

    real(ark) :: eps,s6_old,f
    real(ark) :: rjacob,dx

    real(ark) :: stadev_old,stability,stadev,ssq,stadev_best,fac_sign,s6,h,dx0
    !
    integer(ik) :: iter,itmax,i
    !
   
    h = 1e-3
    rjacob = 0 
    iter = 0
    stadev_old = 1.e10
    stability =  1.e10
    stadev    =  1.e10

    stadev_best = sqrt(small_)*10.0_ark
    itmax = 1000
    !
    ! Initial value for alpha1
    !
    alpha1 = 2.0_ark*pi/3.0_ark
    !
    s6 = alpha1*sqrt(3.0_ark)
    !
    fac_sign = sign(1._ark,sindelta)
    !
    outer_loop: & 
    do while( iter<itmax .and. stadev>stadev_best )   
       !
       iter = iter + 1
       ssq=0
       !
       ! Caclulate the function 
       !
       f = calc_s2sindelta2(s4,s5,s6)
       !
       eps = f - sindelta**2
       !
       ssq=abs(eps)
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       rjacob  = ( calc_s2sindelta2(s4,s5,s6+h)-calc_s2sindelta2(s4,s5,s6-h) )/h*0.5_ark
       !
       !
       if (itmax>=0) then
         !
         dx = eps/rjacob
         !
         stadev=sqrt(ssq)
         !
         s6_old=s6
         !   
         ! Update the pot. parameters to the new values 
         !
         !s6 = s6 - dx 
         !
         do i =1,10
           !
           f = calc_s2sindelta2(s4,s5,s6-dx*real(i,ark)/10._ark)
           !
           if (abs(f - sindelta**2)>ssq) exit
           !
           ssq = abs(f - sindelta**2)
           !
           dx0 = dx*real(i,ark)/10._ark
           !
         enddo
         !
         s6 = s6 - dx0
         !      
         stability=abs( (stadev-stadev_old)/stadev )
         stadev_old=stadev
             
      else   ! Itmax = 0, so there is no fit
             ! only straightforward calculations 
             !
         stadev_old=stadev
         stadev=sqrt(ssq)
         !
      endif 
      !
    enddo  outer_loop ! --- iter
    !
    if (iter==itmax) then
       write(out,"('find_alpha_from_tau: could not find solution after ',i8,' iterations')") iter
       stop 'find_alpha_from_tau: could not find solution'
    endif 
    !
    alpha1 =(sqrt(2.0_ark)*s6+2.0_ark*s4                 )/sqrt(6.0_ark)
    alpha2 =(sqrt(2.0_ark)*s6-        s4+sqrt(3.0_ark)*s5)/sqrt(6.0_ark)
    alpha3 =(sqrt(2.0_ark)*s6-        s4-sqrt(3.0_ark)*s5)/sqrt(6.0_ark)
    !
  end subroutine find_alpha_from_sindelta
  !
  ! This function computes tau**2 from s4,s5,s6
  ! 
  function calc_s2sindelta2(s4,s5,s6) result (sindelta_2)

    real(ark),intent(in)  :: S4,S5,S6
    real(ark)             :: alpha1,alpha2,alpha3,sindelta_2,tau_2,norm_2
     !
     alpha1 =(sqrt(2.0_ark)*s6+2.0_ark*s4                 )/sqrt(6.0_ark)
     alpha2 =(sqrt(2.0_ark)*s6-        s4+sqrt(3.0_ark)*s5)/sqrt(6.0_ark)
     alpha3 =(sqrt(2.0_ark)*s6-        s4-sqrt(3.0_ark)*s5)/sqrt(6.0_ark)
     !
     tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
            +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)

     norm_2 = sin(alpha3)**2+sin(alpha2)**2+sin(alpha1)**2+&
              2._ark*cos(alpha3)*cos(alpha1)-2._ark*cos(alpha2)+ &
              2._ark*cos(alpha2)*cos(alpha3)-2._ark*cos(alpha1)+ &
              2._ark*cos(alpha2)*cos(alpha1)-2._ark*cos(alpha3)
     !
     !if ( tau_2/norm_2<-sqrt(small_) ) then 
     !   !
     !   write (out,"('calc_s2sindelta: tau**2<0: ',f18.8)") tau_2
     !   stop 'calc_s2sindelta - tau**2<0'
     !   !
     !elseif ( tau_2/norm_2<0.0_ark) then 
     !   sindelta_2 = 0
     !else 
     sindelta_2 = tau_2/norm_2
     !endif
     !
     !
  end function calc_s2sindelta2




  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  subroutine ML_symmetry_transformation_XY3(ioper,nmodes,src,dst)
    !
    integer(ik),intent(in)    :: ioper  ! group operation  
    integer(ik),intent(in)    :: nmodes
    real(ark),intent(in)      :: src(nmodes)
    real(ark),intent(out)     :: dst(nmodes)
    !
    real(ark)         :: repres(12,2,2),a,b,e,o,r_t(2,2),phi1,phi2,phi3,s1,s2,s3,p2,p3
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY3/start')") 
    !
    nsrc = size(src)
    !
    if (size(src)/=6) then
       !
       write (out,"('ML_symmetry_transformation_XY3: bad number of coords ',i8)") size(src)
       !
       stop 'ML_symmetry_transformation_XY3 bad number of coords'
       !
    endif 
    !
    select case(trim(molec%symmetry))
    case default
       write (out,"('ML_symmetry_transformation_XY3: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_symmetry_transformation_XY3 - bad symm. type'
       !
    case('C3V','C3V(M)')
       !
       select case(trim(molec%coords_transform))
           !
           case default
           write (out,"('ML_symmetry_transformation_XY3: coord_transf ',a,' unknown')") trim(molec%coords_transform)
           stop 'ML_symmetry_transformation_XY3 - bad coord. type'
           !
       case('R-ALPHA')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst = src

           case (3) ! (132)

             dst(1) = src(2)
             dst(2) = src(3)
             dst(3) = src(1)
             dst(4) = src(5)
             dst(5) = src(6)
             dst(6) = src(4)

           case (2) ! (123)

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
             dst(4) = src(6)
             dst(5) = src(4)
             dst(6) = src(5)
 
           case (6) ! (12)

             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(3)
             dst(4) = src(5)
             dst(5) = src(4)
             dst(6) = src(6)

           case (5) ! (13)

             dst(1) = src(3)
             dst(2) = src(2)
             dst(3) = src(1)
             dst(4) = src(6)
             dst(5) = src(5)
             dst(6) = src(4)


           case (4) ! (23)

             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)
             dst(4) = src(4)
             dst(5) = src(6)
             dst(6) = src(5)

           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
           end select 
           !
       case('NORMAL','C3V-SYMMETRY')
           !
           dst = src
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
           dst(3) = repres(ioper,1,1)*src(3)+repres(ioper,1,2)*src(4)
           dst(4) = repres(ioper,2,1)*src(3)+repres(ioper,2,2)*src(4)
           !
           dst(5) = repres(ioper,1,1)*src(5)+repres(ioper,1,2)*src(6)
           dst(6) = repres(ioper,2,1)*src(5)+repres(ioper,2,2)*src(6)

           !
       case('R-S-RHO','R-S-DELTA')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst = src

           case (3) ! (123)

             dst(1) = src(2)
             dst(2) = src(3)
             dst(3) = src(1)

           case (2) ! (321)

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
 
           case (6) ! (12)

             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(3)

           case (5) ! (13)

             dst(1) = src(3)
             dst(2) = src(2)
             dst(3) = src(1)

           case (4) ! (23)

             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)

           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
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
           dst(4) = repres(ioper,1,1)*src(4)+repres(ioper,1,2)*src(5)
           dst(5) = repres(ioper,2,1)*src(4)+repres(ioper,2,2)*src(5)
           !
           dst(6) = src(6) 
           !
       case('R-2D-RHO')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst = src

           case (3) ! (123)

             dst(1) = src(2)
             dst(2) = src(3)
             dst(3) = src(1)

           case (2) ! (321)

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
 
           case (6) ! (12)

             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(3)

           case (5) ! (13)

             dst(1) = src(3)
             dst(2) = src(2)
             dst(3) = src(1)

           case (4) ! (23)

             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)

           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
           end select 
           !
           a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
           !
           repres ( 1,:,:)= reshape((/ e, o,  & 
                                       o, e/),(/2,2/))
           !
           repres ( 2,:,:)= reshape((/-a,-b,  &
                                       b,-a/),(/2,2/))
           !
           repres ( 3,:,:)= reshape((/-a, b,  &
                                      -b,-a/),(/2,2/))
           !
           repres ( 4,:,:)= reshape((/ e, o,  &
                                       o,-e/),(/2,2/))
           !
           repres ( 5,:,:)= reshape((/-a, b,  &
                                       b, a/),(/2,2/))
           !
           repres ( 6,:,:)= reshape((/-a,-b,  &
                                      -b, a/),(/2,2/))
           !
           r_t = reshape((/ sqrt(0.5_ark),-sqrt(0.5_ark),  & 
                            sqrt(0.5_ark), sqrt(0.5_ark)/),(/2,2/))
           !
           repres(ioper,:,:) = matmul(matmul(transpose(r_t),repres(ioper,:,:)),r_t)
           !
           dst(4) = repres(ioper,1,1)*src(4)+repres(ioper,1,2)*src(5)
           dst(5) = repres(ioper,2,1)*src(4)+repres(ioper,2,2)*src(5)
           !
           dst(6) = src(6) 
           !
           !
       case('R-SYMPHI-DELTA','R-SYMPHI-TAU','R-S-DELTA-MEP')
           !
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst = src

           case (3) ! (123)

             dst(1) = src(2)
             dst(2) = src(3)
             dst(3) = src(1)

           case (2) ! (321)

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
 
           case (6) ! (12)

             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(3)

           case (5) ! (13)

             dst(1) = src(3)
             dst(2) = src(2)
             dst(3) = src(1)

           case (4) ! (23)

             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)

           end select 
           !
           a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
           !
           repres ( 1,:,:)= reshape((/ e, o,  & 
                                       o, e/),(/2,2/))
           !
           repres ( 2,:,:)= reshape((/-a,-b,  &
                                       b,-a/),(/2,2/))
           !
           repres ( 3,:,:)= reshape((/-a, b,  &
                                      -b,-a/),(/2,2/))
           !
           repres ( 4,:,:)= reshape((/ e, o,  &
                                       o,-e/),(/2,2/))
           !
           repres ( 5,:,:)= reshape((/-a, b,  &
                                       b, a/),(/2,2/))
           !
           repres ( 6,:,:)= reshape((/-a,-b,  &
                                      -b, a/),(/2,2/))
           !
           dst(4) = repres(ioper,1,1)*src(4)+repres(ioper,1,2)*src(5)
           dst(5) = repres(ioper,2,1)*src(4)+repres(ioper,2,2)*src(5)
           !
           dst(6) = src(6) 
           ! 
       end select  
       !
    case('D3H','D3H(M)')
       !
       select case(trim(molec%coords_transform))
           !
           case default
           write (out,"('ML_symmetry_transformation_XY3: coord_transf ',a,' unknown')") trim(molec%coords_transform)
           stop 'ML_symmetry_transformation_XY3 - bad coord. type'
           !
       case('R-S-DELTA','R-SYMPHI-DELTA','R-SYMPHI-TAU','R-S-DELTA-MEP')
           !
           a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
           !
           !
           repres ( 1,:,:)= reshape((/ e, o,  & 
                                       o, e/),(/2,2/))
           !
           repres ( 2,:,:)= reshape((/-a,-b,  &
                                       b,-a/),(/2,2/))
           !
           repres ( 3,:,:)= reshape((/-a, b,  &
                                      -b,-a/),(/2,2/))
           !
           repres ( 4,:,:)= reshape((/ e, o,  &
                                       o,-e/),(/2,2/))
           !
           repres ( 5,:,:)= reshape((/-a, b,  &
                                       b, a/),(/2,2/))
           !
           repres ( 6,:,:)= reshape((/-a,-b,  &
                                      -b, a/),(/2,2/))
           !
           repres ( 7,:,:) = repres (1,:,:)
           repres ( 8,:,:) = repres (2,:,:)
           repres ( 9,:,:) = repres (3,:,:)
           repres (10,:,:) = repres (4,:,:)
           repres (11,:,:) = repres (5,:,:)
           repres (12,:,:) = repres (6,:,:)

           !
           dst(4) = repres(ioper,1,1)*src(4)+repres(ioper,1,2)*src(5)
           dst(5) = repres(ioper,2,1)*src(4)+repres(ioper,2,2)*src(5)
           !
           dst(6) = src(6) ; if (ioper>=4.and.ioper<=9) dst(6) = -dst(6)
           !
           select case(ioper)
           !
           case (1, 7) ! identity,(23)*

             dst(1:3) = src(1:3)
             !
           case (3, 9) ! (123),(123)*

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
             !
           case (2, 8) ! (321),(321)*

             dst(1) = src(2)
             dst(2) = src(3)
             dst(3) = src(1)
 
           case (4,10) ! (23),(23)*

             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)

           case (6,12) ! (13),(13)*

             dst(1) = src(3)
             dst(2) = src(2)
             dst(3) = src(1)

           case (5,11) ! (12),(12)*

             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(3)

           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
           end select 
           !
       case('SYM-DELTA')
           !
           a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
           !
           repres ( 1,:,:)= reshape((/ e, o,  & 
                                       o, e/),(/2,2/))
           !
           repres ( 2,:,:)= reshape((/-a, b,  &
                                      -b,-a/),(/2,2/))
           !
           repres ( 3,:,:)= reshape((/-a,-b,  &
                                       b,-a/),(/2,2/))
           !
           repres ( 4,:,:)= reshape((/ e, o,  &
                                       o,-e/),(/2,2/))
           !
           repres ( 5,:,:)= reshape((/-a,-b,  &
                                      -b, a/),(/2,2/))
           !
           repres ( 6,:,:)= reshape((/-a, b,  &
                                       b, a/),(/2,2/))
           !
           repres ( 7,:,:) = repres (1,:,:)
           repres ( 8,:,:) = repres (2,:,:)
           repres ( 9,:,:) = repres (3,:,:)
           repres (10,:,:) = repres (4,:,:)
           repres (11,:,:) = repres (5,:,:)
           repres (12,:,:) = repres (6,:,:)
           !
           dst(1) = src(1)
           !
           dst(2) = repres(ioper,1,1)*src(2)+repres(ioper,1,2)*src(3)
           dst(3) = repres(ioper,2,1)*src(2)+repres(ioper,2,2)*src(3)
           !
           repres ( 1,:,:)= reshape((/ e, o,  & 
                                       o, e/),(/2,2/))
           !
           repres ( 2,:,:)= reshape((/-a,-b,  &
                                       b,-a/),(/2,2/))
           !
           repres ( 3,:,:)= reshape((/-a, b,  &
                                      -b,-a/),(/2,2/))
           !
           repres ( 4,:,:)= reshape((/ e, o,  &
                                       o,-e/),(/2,2/))
           !
           repres ( 5,:,:)= reshape((/-a, b,  &
                                       b, a/),(/2,2/))
           !
           repres ( 6,:,:)= reshape((/-a,-b,  &
                                      -b, a/),(/2,2/))
           !
           repres ( 7,:,:) = repres (1,:,:)
           repres ( 8,:,:) = repres (2,:,:)
           repres ( 9,:,:) = repres (3,:,:)
           repres (10,:,:) = repres (4,:,:)
           repres (11,:,:) = repres (5,:,:)
           repres (12,:,:) = repres (6,:,:)

           !
           dst(4) = repres(ioper,1,1)*src(4)+repres(ioper,1,2)*src(5)
           dst(5) = repres(ioper,2,1)*src(4)+repres(ioper,2,2)*src(5)


           !
           !
           ! new !!!! 
           !

           a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
           !
           !
           repres ( 1,:,:)= reshape((/ e, o,  & 
                                       o, e/),(/2,2/))
           !
           repres ( 2,:,:)= reshape((/-a,-b,  &
                                       b,-a/),(/2,2/))
           !
           repres ( 3,:,:)= reshape((/-a, b,  &
                                      -b,-a/),(/2,2/))
           !
           repres ( 4,:,:)= reshape((/ e, o,  &
                                       o,-e/),(/2,2/))
           !
           repres ( 5,:,:)= reshape((/-a, b,  &
                                       b, a/),(/2,2/))
           !
           repres ( 6,:,:)= reshape((/-a,-b,  &
                                      -b, a/),(/2,2/))
           !
           repres ( 7,:,:) = repres (1,:,:)
           repres ( 8,:,:) = repres (2,:,:)
           repres ( 9,:,:) = repres (3,:,:)
           repres (10,:,:) = repres (4,:,:)
           repres (11,:,:) = repres (5,:,:)
           repres (12,:,:) = repres (6,:,:)
           !
           dst(2) = repres(ioper,1,1)*src(2)+repres(ioper,1,2)*src(3)
           dst(3) = repres(ioper,2,1)*src(2)+repres(ioper,2,2)*src(3)
           !
           dst(4) = repres(ioper,1,1)*src(4)+repres(ioper,1,2)*src(5)
           dst(5) = repres(ioper,2,1)*src(4)+repres(ioper,2,2)*src(5)
           !
           dst(6) = src(6) ; if (ioper>=4.and.ioper<=9) dst(6) = -dst(6)

           !
           !dst(6) = src(6) ; if (ioper>=4.and.ioper<=9) dst(6) = -dst(6)
           !           
       case('R-2D-DELTA')
           !
           a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
           !
           repres ( 1,:,:)= reshape((/ e, o,  & 
                                       o, e/),(/2,2/))
           !
           repres ( 2,:,:)= reshape((/-a,-b,  &
                                       b,-a/),(/2,2/))
           !
           repres ( 3,:,:)= reshape((/-a, b,  &
                                      -b,-a/),(/2,2/))
           !
           repres ( 4,:,:)= reshape((/ e, o,  &
                                       o,-e/),(/2,2/))
           !
           repres ( 5,:,:)= reshape((/-a, b,  &
                                       b, a/),(/2,2/))
           !
           repres ( 6,:,:)= reshape((/-a,-b,  &
                                      -b, a/),(/2,2/))
           !
           repres (10,:,:) = repres (1,:,:)
           repres (11,:,:) = repres (2,:,:)
           repres (12,:,:) = repres (3,:,:)
           repres ( 7,:,:) = repres (4,:,:)
           repres ( 8,:,:) = repres (5,:,:)
           repres ( 9,:,:) = repres (6,:,:)
           !
           r_t = reshape((/ sqrt(0.5_ark),-sqrt(0.5_ark),  & 
                            sqrt(0.5_ark), sqrt(0.5_ark)/),(/2,2/))
           !
           repres(ioper,:,:) = matmul(matmul(transpose(r_t),repres(ioper,:,:)),r_t)
           !
           dst(4) = repres(ioper,1,1)*src(4)+repres(ioper,1,2)*src(5)
           dst(5) = repres(ioper,2,1)*src(4)+repres(ioper,2,2)*src(5)
           !
           dst(6) = src(6) ; if (ioper>=4.and.ioper<=9) dst(6) = -dst(6)
           !
           select case(ioper)
           !
           case (1,7) ! identity 

             dst(1:3) = src(1:3)
             !
           case (2,8) ! (321)

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
             !
           case (3,9) ! (123)

             dst(1) = src(2)
             dst(2) = src(3)
             dst(3) = src(1)
 
           case (4,10) ! (23)

             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)

           case (5,11) ! (13)

             dst(1) = src(3)
             dst(2) = src(2)
             dst(3) = src(1)

           case (6,12) ! (12)

             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(3)

           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
           end select 
           !           
       case('R-PHI-DELTA')
           !
           a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
           !
           !
           repres ( 1,:,:)= reshape((/ e, o,  & 
                                       o, e/),(/2,2/))
           !
           repres ( 2,:,:)= reshape((/-a,-b,  &
                                       b,-a/),(/2,2/))
           !
           repres ( 3,:,:)= reshape((/-a, b,  &
                                      -b,-a/),(/2,2/))
           !
           repres ( 4,:,:)= reshape((/ e, o,  &
                                       o,-e/),(/2,2/))
           !
           repres ( 5,:,:)= reshape((/-a, b,  &
                                       b, a/),(/2,2/))
           !
           repres ( 6,:,:)= reshape((/-a,-b,  &
                                      -b, a/),(/2,2/))
           !
           repres ( 7,:,:) = repres (1,:,:)
           repres ( 8,:,:) = repres (2,:,:)
           repres ( 9,:,:) = repres (3,:,:)
           repres (10,:,:) = repres (4,:,:)
           repres (11,:,:) = repres (5,:,:)
           repres (12,:,:) = repres (6,:,:)
           !
           phi2 = src(4) + 2.0_ark/3.0_ark*pi
           phi3 = src(5) + 2.0_ark/3.0_ark*pi
           !
           phi1 = 2.0_ark*pi-phi2-phi3
           !
           s1 = 2.0_ark*pi/sqrt(3.0_ark)
           s2 = 1.0_ark/sqrt(6.0_ark)*(2.0_ark*phi1-phi2-phi3)
           s3 = 1.0_ark/sqrt(2.0_ark)*(phi2-phi3)
           !
           p2 = repres(ioper,1,1)*s2+repres(ioper,1,2)*s3
           p3 = repres(ioper,2,1)*s2+repres(ioper,2,2)*s3
           !
           phi1   =(sqrt(2.0_ark)*s1+2.0_ark*p2                 )/sqrt(6.0_ark)
           phi2   =(sqrt(2.0_ark)*s1-        p2+sqrt(3.0_ark)*p3)/sqrt(6.0_ark)
           phi3   =(sqrt(2.0_ark)*s1-        p2-sqrt(3.0_ark)*p3)/sqrt(6.0_ark)
           !
           dst(4) = phi2 - 2.0_ark/3.0_ark*pi
           dst(5) = phi3 - 2.0_ark/3.0_ark*pi
           !
           dst(6) = src(6) ; if (ioper>=4.and.ioper<=9) dst(6) = -dst(6)
           !
           select case(ioper)
           !
           case (1,7) ! identity 

             dst(1:3) = src(1:3)
             !
           case (2,8) ! (321)

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
             !
           case (3,9) ! (123)

             dst(1) = src(2)
             dst(2) = src(3)
             dst(3) = src(1)
 
           case (4,10) ! (23)

             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)

           case (5,11) ! (13)

             dst(1) = src(3)
             dst(2) = src(2)
             dst(3) = src(1)

           case (6,12) ! (12)

             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(3)

           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
           end select 
           !
       end select
       !
    case('CS','CS(M)')
       !
       select case(trim(molec%coords_transform))
           !
           case default
           write (out,"('ML_symmetry_transformation_XY3: coord_transf ',a,' unknown')") trim(molec%coords_transform)
           stop 'ML_symmetry_transformation_XY3 - bad coord. type'
           !
       case('R-ALPHA')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst = src

           case (2) ! (sigma)

             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)
             dst(4) = src(4)
             dst(5) = src(6)
             dst(6) = src(5)

           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
           end select 
           !
       case('R-S-DELTA','R-PHI-DELTA','R-A2-A3-TAU','R-S-DELTA-MEP','R-PHI-DELTA-MEP','R-A2A3-DELTA')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst = src

           case (2) ! (sigma)

             dst(1) = src(1)
             dst(2) = src(2)
             dst(3) = src(3)
             dst(4) = src(4)
             dst(5) = src(5)
             dst(6) =-src(6)

           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
           end select            
           !
       case('XXXXXX')
           !
           select case(ioper)
           !
           case (1) ! identity 
             !
             dst = src
             !
           case (2) ! (sigma)
             !
             dst(1:5) = src(1:5)
             dst(6) = -src(6)
             !
           case default

             write (out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
           end select 
           !
       case('NORMAL')
           !
           dst = src
           !
           if(ioper==2) then
             !
             dst(4) =-src(4)
             dst(6) =-src(6)
             !
           endif
           !
       end select
       !
    case('C2V','C2V(M)')
       !
       select case(trim(molec%coords_transform))
           !
           case default
           write (out,"('ML_symmetry_transformation_XY3: coord_transf ',a,' unknown')") trim(molec%coords_transform)
           stop 'ML_symmetry_transformation_XY3 - bad coord. type'
           !
       case('R-A2-A3-TAU','R-PHI-DELTA','R-S-DELTA-MEP','R-PHI-DELTA-MEP','R-A2A3-DELTA')
          !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) =-src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) =-src(6)

         case (4) ! (12)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         end select 
         !
       case('R-SYMPHI-DELTA-MEP','R-SYMPHI-DELTA','R-S-DELTA')
          !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(4)
           dst(5) =-src(5)
           dst(6) =-src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) =-src(6)

         case (4) ! (12)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(4)
           dst(5) =-src(5)
           dst(6) = src(6)

         end select 
         !
       end select
       !
    end select
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY3/end')") 
    !
  end subroutine ML_symmetry_transformation_XY3


  ! Here we find the rotational symmetry from J,K
  !
  subroutine ML_rotsymmetry_XY3(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_XY3/start')") 
    !
    !
    select case(trim(molec%symmetry))
    case default
       write (out,"('ML_rotsymmetry_XY3: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_rotsymmetry_XY3 - bad symm. type'
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
    case('D3H','D3H(M)')
       !
       gamma = 0 
       ideg = 1 
       !
       if     (mod(K+3,3)==0.and.tau==0.and.mod(k+2,2)==0) then 
          gamma = 1 
       elseif (mod(K+3,3)==0.and.tau==1.and.mod(k+2,2)==0) then 
          gamma = 2
       elseif (mod(K+3,3)==0.and.tau==0.and.mod(k+2,2)/=0) then 
          gamma = 5
       elseif (mod(K+3,3)==0.and.tau==1.and.mod(k+2,2)/=0) then 
          gamma = 4
       elseif (mod(K+3,3)/=0.and.tau==0.and.mod(k+2,2)==0) then 
          gamma = 3 ; ideg = 1 
       elseif (mod(K+3,3)/=0.and.tau==1.and.mod(k+2,2)==0) then 
          gamma = 3 ; ideg = 2
       elseif (mod(K+3,3)/=0.and.tau==0.and.mod(k+2,2)/=0) then 
          gamma = 6 ; ideg = 1 
       elseif (mod(K+3,3)/=0.and.tau==1.and.mod(k+2,2)/=0) then 
          gamma = 6 ; ideg = 2
       else
          !
          write(out,"('ML_rotsymmetry_XY3-D3h: illegal j,k,tau - ',3i8)") j,k,tau
          !
       endif 
       !
    case('CS','CS(M)')
       !
       gamma = 0 
       ideg = 1 
       !
       if (mod(tau+2,2)==0) gamma = 1 !; return
       if (mod(tau+2,2)/=0) gamma = 2 !; return
       !
       !ideg = 1 
       !
       !if (mod(K+3,3)==0.and.tau==0) gamma = 1 !; return
       !if (mod(K+3,3)==0.and.tau==1) gamma = 2 !; return
       !
       !if (mod(K+3,3)/=0.and.tau==0) then 
       !   gamma = 1 ; ideg = 1 
       !endif 
       !if (mod(K+3,3)/=0.and.tau==1) then 
       !   gamma = 2 ; ideg = 1
       !endif 
       !
    case('C2V','C2V(M)')
       !
       !
       select case(trim(molec%frame))
          !
       case default
          !
          gamma = 0 
          ideg = 1
          !
          ! R1-X-R2-Z-Y
          ! 
          ! for the planer configuration, x is along the symmetry axis, z is in the plane and y is orhogonal to the plane
          !
          if (mod(K+2,2)==0.and.tau==0) gamma = 1 !; 1;1
          if (mod(K+2,2)==0.and.tau==1) gamma = 3 !; 2;2
          if (mod(K+2,2)/=0.and.tau==0) gamma = 4 !; 3;4
          if (mod(K+2,2)/=0.and.tau==1) gamma = 2 !; 4;3
          !
       case('R1-X-R2-Y-Z')
          ! 
          ! for the planer configuration, x is along the symmetry axis, y is in the plane and y is orhogonal to the plane
          !
          gamma = 0 
          ideg = 1
          if (mod(K+2,2)==0.and.tau==0) gamma = 1 !;  1;1;1;1;1
          if (mod(K+2,2)==0.and.tau==1) gamma = 4 !;  2;2;3;3;4
          if (mod(K+2,2)/=0.and.tau==0) gamma = 3 !;  3;4;4;2;2
          if (mod(K+2,2)/=0.and.tau==1) gamma = 2 !;  4;3;2;4;3
          !
       case('R1-Z-R2-X-Y')
          ! 
          ! for the planer configuration, z is along the symmetry axis, x is in the plane and y is orhogonal to the plane
          !
          gamma = 0 
          ideg = 1
          if (mod(K+2,2)==0.and.tau==0) gamma = 1 !; return
          if (mod(K+2,2)==0.and.tau==1) gamma = 2 !; return
          if (mod(K+2,2)/=0.and.tau==0) gamma = 4 !; return
          if (mod(K+2,2)/=0.and.tau==1) gamma = 3 !; return
          !
       end select 
       !
    case('C','C(M)')
       !
       gamma = 1
       ideg = 1 
       !
    end select
    !
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_XY3/end')") 
    !
    !
  end subroutine ML_rotsymmetry_XY3

  !
  ! Defining MEP function for OH3+ molecule
  !
  function ML_MEP_OH3P(x)  result(dst)

   real(ark),intent(in)   ::  x
   real(ark)              ::  dst
   real(ark)              ::  rf(0:6),y
   integer(ik)          ::  k(0:6) = (/0,1,2,3,4,5,6/)
     !
     rf(0:6)  = molec%mep_params(1:7) 
     !
     y=sin(x)-1.0_ark
     !
     dst =  sum(rf(:)*y**k(:))
     ! 
  end function ML_MEP_OH3P


  !
  ! Defining MEP function for NH3 molecule
  !
  function ML_MEP_NH3(x)  result(dst)

   real(ark),intent(in)   ::  x
   real(ark)              ::  dst(1:molec%ncoords)
   real(ark)              ::  rf(0:6),y
   integer(ik)          ::  k(0:6) = (/0,1,2,3,4,5,6/)
     !
     rf(0:6)  = molec%mep_params(1:7) 
     !
     y=cos(x)-1.0_ark
     !
     dst(1:3) =  sum(rf(:)*y**k(:))
     dst(4:) = 2.0_ark*asin( 0.5_ark*sqrt(3.0_ark)*cos(x) )

     dst(molec%ncoords) = x
     ! 
  end function ML_MEP_NH3



end module mol_xy3
    
    
    
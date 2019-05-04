!
!  This unit defines all specific routines for a triatomic molecule of XY2 type
!
module mol_xy2
  use accuracy
  use moltype
  use symmetry,only : sym

  implicit none

  public ML_b0_XY2,ML_coordinate_transform_XY2,ML_symmetry_transformation_XY2,ML_rotsymmetry_XY2
  !
  private
  ! 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  !
  contains
  !
  !
  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  function ML_coordinate_transform_XY2(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: c_t(3,3),dsrc(size(src)),r_e1(2),r_e2(2),radau_e1(2),radau_e2(2),radau_e(2)
    real(ark)                 :: m1,m2,m3,a0,r_t1(2),r_t2(2),radau_t1(2),radau_t2(2),cosalpha,rho,alpha,sinalpha
    real(ark)                 :: r1,r2,r3,R,r12,r12e,r13e,alpha1,alpha2,theta,sintheta,q0,q,q1,q2,alpha0
    !
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_XY2/start')") 
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
       write (out,"('ML_coordinate_transform_XY2: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'ML_coordinate_transform_XY2 - bad coord. type'
       !
    case('C2V-SYMMETRY')
       !
       c_t =  reshape( &
     (/1.0_ark/sqrt2, 1.0_ark/sqrt2      , 0.0_ark    , &     ! 1
       1.0_ark/sqrt2,-1.0_ark/sqrt2      , 0.0_ark    , &     ! 3
       0.0_ark      , 0.0_ark            , 1.0_ark/)  , &     ! 2
       (/3,3/))
       !
       if (direct) then 
           dst = matmul(c_t,dsrc)
       else
           dst = matmul(transpose(c_t),dsrc)
           dst(:) = dst(:) +  molec%local_eq(:)
       endif
       !
    case('R-RHO','R-RHO-Z')
       !
       if (direct) then 
          dst(1:2) = dsrc(1:2)
          dst(3) =  pi-src(3)
       else
          dst(1:2) = src(1:2)+molec%local_eq(1:2)
          dst(3) = pi-src(3)
       endif
       !
    case('R-RHO-HALF')
       !
       if (direct) then 
          dst(1:2) = dsrc(1:2)
          dst(3) =  (pi-src(3))*0.5_ark
       else
          dst(1:2) = src(1:2)+molec%local_eq(1:2)
          dst(3) = pi-src(3)*2.0_ark
       endif
       !
    case('R-PHI-RHO','R-PHI-RHO-Z')
       !
       if (direct) then 
          dst(1:2) = dsrc(1:2)
          dst(3) =  src(4)
          dst(4) =  src(3)
       else
          dst(1:2) = src(1:2)+molec%local_eq(1:2)
          dst(3) = src(4)
          dst(4) = src(3)
       endif
       !
    case('R1-R2-Y+X')
       !
       if (direct) then 
          !
          dst(1:2) = dsrc(1:2)
          !
          dst(3) =-src(4)-molec%local_eq(3)
          dst(4) = src(3)-molec%local_eq(4)
          !
       else
          !
          dst(1:2) = src(1:2)+molec%local_eq(1:2)
          !
          dst(3) = src(4)+molec%local_eq(4)
          dst(4) =-src(3)-molec%local_eq(3)
          !
       endif       
       !
    case('R-PHI1-PHI2','R-PHI1-PHI2-Z')
       !
       if (direct) then 
          dst(1:4) = dsrc(1:4)
       else
          dst(1:4) = src(1:4)+molec%local_eq(1:4)
       endif
       !
    case('R-ALPHA-THETA-Z')
       !
       alpha0 = pi-asin( sqrt( sin(molec%taueq(1))**2+sin(molec%taueq(2))**2 ))
       q0 = sqrt(molec%taueq(1)**2+molec%taueq(2)**2 )
       !
       if (direct) then 
          dst(1:2) = dsrc(1:2)
          q1 = src(3)
          q2 = src(4)
          dst(3) = sqrt(q1**2+q2**2) - q0
          dst(4) = atan2(q2,q1)
          !
       else
          q = src(3)+q0
          theta = src(4)
          dst(1:2) = src(1:2)+molec%local_eq(1:2)
          dst(3) = q*cos(theta)
          dst(4) = q*sin(theta)
       endif
       !
    case('R-S1-S2','R-S1-S2-Z')
       !
       if (direct) then 
          dst(1:2) = dsrc(1:2)
          !
          dst(3) = (dsrc(3)-dsrc(4))/sqrt(2.0_ark)
          dst(4) = (dsrc(3)+dsrc(4))/sqrt(2.0_ark)
          !
       else
          dst(1:2) = src(1:2)+molec%local_eq(1:2)
          !
          dst(3) = (src(3)+src(4))/sqrt(2.0_ark)+molec%local_eq(3)
          dst(4) = (src(4)-src(3))/sqrt(2.0_ark)+molec%local_eq(4)
          !
       endif
       !
    case('R-PHI1','R-PHI1-Z')
       !
       if (direct) then 
          dst(1:3) = dsrc(1:3)
       else
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
       endif
       !
    case('R12-R')
       !
       if (direct) then 
          !
          r1 = src(1) ; r2 = src(2) ;  alpha = src(3)
          !
          r12 = sqrt(r1**2+r2**2-2.0_ark*r1*r2*cos(alpha))
          !
          !sinalpha = r1/r12*sin(alpha)
          !alpha1 = 0
          !if (sinalpha>small_) alpha1 = asin(sinalpha)
          !!
          !sinalpha = r2/r12*sin(alpha)
          !alpha2 = 0
          !if (sinalpha>small_) alpha2 = asin(sinalpha)
          !
          R = 0.5_ark*sqrt(r1**2+r2**2+2.0_ark*r1*r2*cos(alpha))
          !
          !sintheta = r1/R*sin(alpha2)
          !theta = 0
          !
          !if ( abs(sintheta)>1.0_ark+sqrt(small_) ) then 
          !   !
          !   write (out,"('MLpoten_xy2_R_R12_theta: sintheta>1: ',f18.8)") sintheta
          !   stop 'MLpoten_xy2_R_R12_theta - bad sintheta'
          !   !
          !elseif ( sintheta>=1.0_ark) then 
          !   theta = 0.5_ark*pi
          !else 
          !   theta = asin(sintheta)
          !endif
          !
          r12e = ML_MEP_xy2_R12_R(R)
		  !
		  !r12e = sqrt(2.0_ark*molec%local_eq(1)**2-2.0_ark*molec%local_eq(1)**2*cos(molec%local_eq(3)))
          !
          dst(1) = r12 -r12e
          dst(2) = (r1-r2)*0.5_ark
          dst(3) = R
          !
       else
          !
          R = src(3)
          !
          r12e = ML_MEP_xy2_R12_R(R)
		  !
		  !r12e = sqrt(2.0_ark*molec%local_eq(1)**2-2.0_ark*molec%local_eq(1)**2*cos(molec%local_eq(3)))
          !
          !theta = src(2) + pi*0.5_ark
          !
          r12 = src(1) + r12e
          !
          r1 = sqrt(R**2-src(2)**2+0.25_ark*r12**2)+src(2)
          r2 = sqrt(R**2-src(2)**2+0.25_ark*r12**2)-src(2)
          !
          !r1 = sqrt(R**2+0.25_ark*r12**2-R*r12*cos(theta))
          !r2 = sqrt(R**2+0.25_ark*r12**2+R*r12*cos(theta))
          !
          cosalpha = (r1**2+r2**2-r12**2)/(2.0_ark*r1*r2)
          !
          if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('ML_coordinate_transform_XY2: cosalpha>1: ',f18.8)") cosalpha
             stop 'ML_coordinate_transform_XY2 - bad cosalpha'
             !
          elseif ( cosalpha>=1.0_ark) then 
             alpha = 0
          else 
             alpha = acos(cosalpha)
          endif
          !
          dst(1) = r1
          dst(2) = r2
          dst(3) = alpha
          !
       endif
       !
    case('R-R-R')
       !
       if (direct) then 
          !
          r1 = src(1) ; r2 = src(2) ;  alpha = src(3)
          !
          r3 = sqrt(r1**2+r2**2-2.0_ark*r1*r2*cos(alpha))
          !
          dst(1) = r1-molec%local_eq(1)
          dst(2) = r2-molec%local_eq(1)
          dst(3) = r3-molec%local_eq(1)
          !
       else
          !
          r1 = src(1)+molec%local_eq(1)
          r2 = src(2)+molec%local_eq(1)
          r3 = src(3)+molec%local_eq(1)
          !R = src(3)+molec%local_eq(1)
          !
          cosalpha = (r1**2+r2**2-r3**2)/(2.0_ark*r1*r2)
          !
          dst(1) = r1
          dst(2) = r2
          dst(3) = acos(cosalpha)
          !
       endif
       !
    case('R12-RHO')
       !
       if (direct) then 
          !
          r1 = src(1) ; r2 = src(2) ;  alpha = src(3)
          !
          r12e = ML_MEP_xy2_R12_ALPHA(alpha)
          !
          dst(1) = r1-r12e
          dst(2) = r2-r12e
          dst(3) = pi-alpha
          !
       else
          !
          alpha = pi-src(3)
          !
          r12e = ML_MEP_xy2_R12_ALPHA(alpha)
          !
          r1 = src(1) + r12e
          r2 = src(2) + r12e
          !
          dst(1) = r1
          dst(2) = r2
          dst(3) = alpha
          !
       endif
       !
    case('R13-RHO')
       !
       if (direct) then 
          !
          r1 = src(1) ; r2 = src(2) ;  alpha = src(3)
          !
          r13e = ML_MEP_xy2_R13_ALPHA(alpha)
          r12e = 0.5_ark*r13e/sin(alpha*0.5_ark)
          !
          dst(1) = r1-r12e
          dst(2) = r2-r12e
          dst(3) = pi-alpha
          !
       else
          !
          alpha = pi-src(3)
          !
          r13e = ML_MEP_xy2_R13_ALPHA(alpha)
          r12e = 0.5_ark*r13e/sin(alpha*0.5_ark)
          !
          r1 = src(1) + r12e
          r2 = src(2) + r12e
          !
          dst(1) = r1
          dst(2) = r2
          dst(3) = alpha
          !
       endif
       !
    case('R12-R-LIN')
       !
       if (direct) then 
          !
          r1 = src(1) ; r2 = src(2) ;  alpha = src(3)
          !
          r12 = sqrt(r1**2+r2**2-2.0_ark*r1*r2*cos(alpha))
          !
          R = 0.5_ark*sqrt(r1**2+r2**2+2.0_ark*r1*r2*cos(alpha))
		  !
		  !r1x = 0.5_ark*r1*r2*sin(alpha)/R
		  !r1y = 0.5_ark*r1*(r1+r2*cos(alpha))/R
		  !r2x =-0.5_ark*r1*r2*sin(alpha)/R
		  !r2y = 0.5_ark*r2*(r2+r1*cos(alpha))/R
		  !
          r12e = ML_MEP_xy2_R12_R(R)
		  !
		  !r12e = sqrt(2.0_ark*molec%local_eq(1)**2-2.0_ark*molec%local_eq(1)**2*cos(molec%local_eq(3)))
          !
          !dst(1) = r1x-0.5_ark*r12e-(r2x+0.5_ark*r12e)
          dst(2) = (r1-r2)*0.5_ark
          dst(3) = R
          !
       else
          !
          R = src(3)
          !
          r12e = ML_MEP_xy2_R12_R(R)
		  !
		  !r12e = sqrt(2.0_ark*molec%local_eq(1)**2-2.0_ark*molec%local_eq(1)**2*cos(molec%local_eq(3)))
          !
          !theta = src(2) + pi*0.5_ark
          !
          r12 = src(1) + r12e
          !
          r1 = sqrt(R**2-src(2)**2+0.25_ark*r12**2)+src(2)
          r2 = sqrt(R**2-src(2)**2+0.25_ark*r12**2)-src(2)
          !
          !r1 = sqrt(R**2+0.25_ark*r12**2-R*r12*cos(theta))
          !r2 = sqrt(R**2+0.25_ark*r12**2+R*r12*cos(theta))
          !
          cosalpha = (r1**2+r2**2-r12**2)/(2.0_ark*r1*r2)
          !
          if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('ML_coordinate_transform_XY2: cosalpha>1: ',f18.8)") cosalpha
             stop 'ML_coordinate_transform_XY2 - bad cosalpha'
             !
          elseif ( cosalpha>=1.0_ark) then 
             alpha = 0
          else 
             alpha = acos(cosalpha)
          endif
          !
          dst(1) = r1
          dst(2) = r2
          dst(3) = alpha
          !
       endif
        !
    case('R-EXPRHO')
       !
       if (direct) then 
	      !
          dst(1:2) = dsrc(1:2)
		  !
		  rho = max(pi-src(3),small_)
		  !
          dst(3) = log(rho)
		  !
		  !dst(3) = (1.0_ark+rho**2)/rho
		  !
       else
	      !
          dst(1:2) = src(1:2)+molec%local_eq(1:2)
		  !
          rho = exp(src(3))
		  !
		  !rho = (src(3))**(0.1)
          !
          !rho = 0.5_ark*(src(3)-sqrt( src(3)**2-4.0_ark ))
		  !
          dst(3) = pi-rho
		  !
       endif
       !
    case('RADAU')
       !
       m1 = molec%AtomMasses(2) ; m2 = molec%AtomMasses(3) ; m3 = molec%AtomMasses(1)
       !
       a0=sqrt(m3/(m1+m2+m3))
       !
       r_e1(1) = molec%local_eq(1)
       r_e1(2) = 0
       !
       r_e2(1) = molec%local_eq(2)*cos(molec%local_eq(3))
       r_e2(2) = molec%local_eq(2)*sin(molec%local_eq(3))
	   !
       radau_e1(:) = ( 1.0_ark+(a0-1.0_ark)*m1/(m1+m2) )*r_e1(:)+(a0-1.0_ark)*m2/(m1+m2)*r_e2(:)
       radau_e2(:) = ( 1.0_ark+(a0-1.0_ark)*m2/(m1+m2) )*r_e2(:)+(a0-1.0_ark)*m1/(m1+m2)*r_e1(:)
	   !
       radau_e(1) = sqrt( sum( radau_e1(:)**2 )  )
       radau_e(2) = sqrt( sum( radau_e2(:)**2 )  )
       !
       if (direct) then 
          !
		  alpha1 = src(3)
          !
          r_t1(1) = src(1)
          r_t1(2) = 0
          !
          r_t2(1) = src(2)*cos(alpha1)
          r_t2(2) = src(2)*sin(alpha1)
          !
          radau_t1(:) = ( 1.0_ark+(a0-1.0_ark)*m1/(m1+m2) )*r_t1(:)+(a0-1.0_ark)*m2/(m1+m2)*r_t2(:)
          radau_t2(:) = ( 1.0_ark+(a0-1.0_ark)*m2/(m1+m2) )*r_t2(:)+(a0-1.0_ark)*m1/(m1+m2)*r_t1(:)
          !
          dst(1) = sqrt( sum( radau_t1(:)**2 )  )
          dst(2) = sqrt( sum( radau_t2(:)**2 )  )

          cosalpha = sum(radau_t1(:)*radau_t2(:) )/( dst(1)*dst(2) )
          !
          if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('ML_coordinate_transform_XY2: cosalpha>1: ',f18.8)") cosalpha
             stop 'ML_coordinate_transform_XY2 - bad cosalpha'
             !
          elseif ( cosalpha>=1.0_ark) then 
             dst(3) = pi
          else 
             dst(3) = pi-acos(cosalpha)
          endif
          !
          dst(1) = dst(1)-radau_e(1)
          dst(2) = dst(2)-radau_e(2)
          !
          !dst(3) = 1.0_ark+cosalpha
          !
       else
          !
          !alpha1 = acos(src(3)-1.0_ark)
		  alpha1 = pi-src(3)
		  !
		  dsrc(1:2) = src(1:2) + radau_e(1:2)
          !
          radau_t1(1) = dsrc(1)
          radau_t1(2) = 0
          !
          radau_t2(1) = dsrc(2)*cos(alpha1)
          radau_t2(1) = dsrc(2)*sin(alpha1)
          !
          r_t1(:) = ( 1.0_ark+(1.0_ark/a0-1.0_ark)*m1/(m1+m2) )*radau_t1(:)+(1.0_ark/a0-1.0_ark)*m2/(m1+m2)*radau_t2(:)
          r_t2(:) = ( 1.0_ark+(1.0_ark/a0-1.0_ark)*m2/(m1+m2) )*radau_t2(:)+(1.0_ark/a0-1.0_ark)*m1/(m1+m2)*radau_t1(:)
          !
          dst(1) = sqrt( sum( r_t1(:)**2 )  )
          dst(2) = sqrt( sum( r_t2(:)**2 )  )

          cosalpha = sum(r_t1(:)*r_t2(:) )/( dst(1)*dst(2) )
          !
          if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('ML_coordinate_transform_XY2: cosalpha>1: ',f18.8)") cosalpha
             stop 'ML_coordinate_transform_XY2 - bad cosalpha'
             !
          elseif ( cosalpha>=1.0_ark) then 
             dst(3) = 0.0_ark
          else 
             dst(3) = acos(cosalpha)
          endif
          !
       endif
       !
    case('R-LINRHO')
       !
       if (direct) then 
          dst(1:2) = dsrc(1:2)
          dst(3) =  pi-src(3)-molec%local_eq(3)
       else
          dst(1:2) = src(1:2)+molec%local_eq(1:2)
          dst(3) = pi-(src(3)+molec%local_eq(3))
       endif
       !
    end select
    !
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_XY2/end')") 
    !
    !
  end function ML_coordinate_transform_XY2


  ! Here we define structural parameters for an XY2 molecule,
  ! a0/b0
  ! which determine the reference geometry
  !
  subroutine ML_b0_XY2(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
     !
     integer(ik),intent(in)  :: Npoints,Natoms
     !
     real(ark),   intent(out) :: b0(molec%Natoms,3,0:Npoints)
     real(ark),   intent(in),optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional :: rho_ref
     real(ark),   intent(in),optional  :: rho_borders(2)  ! rhomim, rhomax - borders
     !
     real(ark)                :: re13,m1,m2,m3,rho,m,rho0,re23,u2,u3,u23,e,r12
     real(ark)                :: alpha,alphae_h,a02,cosa,a,b,g1,g2,rot(3,3),alphaeq
     real(ark)               :: rho_ark,alpha_ark,transform(3,3),CM_shift,a_t,a0(molec%Natoms,3),phi
     real(ark)               :: tmat(3,3),transform0(3,3),transform1(3,3),sinrho,hx,hy,unit(3),tanphi,Inert(3),tau
     integer(ik)             :: i,Nbonds,ix,jx
     integer(ik)             :: ijk(3,2),ioper,Nangles,iatom
     character(cl)           :: method
      !
      !
      if (verbose>=4) write(out,"('ML_b0_XY2/start')")
      !
      if (natoms/=molec%Natoms) stop 'Illegal number of atoms'
      !
      Nbonds = molec%Nbonds
      Nangles = molec%Nangles
      !
      if (Nbonds/=2) then
        write(out,"('ML_b0_XY2: Nbonds must be 2 in this routine, not  ',i9)") Nbonds
        stop 'ML_a0_XY2: wrong Nbonds '
      endif
      ! 
      if (Natoms/=3) then
        write(out,"('ML_b0_XY2: Natoms must be 3 in this routine, not  ',i9)") Natoms
        stop 'ML_a0_XY2: wrong Natoms '
      endif 

      !if (molec%req(1)/=molec%req(2)) then
      !  write(out,"('ML_b0_XY2: req-s must be equal: ',3f14.6)") molec%req(1:2)
      !  stop 'ML_a0_XY2: req-s must be equal: '
      !endif 
      !
      re13 = molec%req(1)
      re23 = molec%req(2)
      !
      if (Nangles>0) then
        alphaeq = molec%alphaeq(1)
      elseif (molec%Ndihedrals>1) then 
        !alphaeq = pi-(molec%taueq(1)**2+molec%taueq(2)**2) !asin(sin(molec%taueq(1))*sqrt(2.0_ark))
        !
        alphaeq = pi-asin(sqrt(molec%taueq(1)**2+molec%taueq(2)**2))
        !
      else
        alphaeq = pi-molec%taueq(1)
      endif 
      alphae_h = alphaeq*0.5_ark
      !
      m1 = molec%AtomMasses(1) ; m2 = molec%AtomMasses(2) ; m3 = molec%AtomMasses(3)
      !
      m = m1+m2+m3
      !
      u3 =m3*(m2+m1)*re13**2
      u2 =m2*(m1+m3)*re23**2
      u23=m3*m2*re13*re23
      !
      b0(1,1,0) = 0.0_ark
      b0(1,2,0) = 0.0_ark
      b0(1,3,0) = 2.0_ark*m3/m*re13*cos(alphae_h)
      !
      b0(2,1,0) =-re13*sin(alphae_h)
      b0(2,2,0) = 0.0_ark
      b0(2,3,0) =-m1/m*re13*cos(alphae_h)
      !
      b0(3,1,0) = re13*sin(alphae_h)
      b0(3,2,0) = 0.0_ark
      b0(3,3,0) =-m1/m*re13*cos(alphae_h)
      !
      rho = pi - alphaeq
      !
      ! Sayvetz : 
      !
      !e = rho*0.5_ark+(u3-u2)/sqrt((u3+u2)**2-4.0_ark*u23**2)*atan(sqrt((u3+u2-2.0_ark*u23)/(u3+u2+2*u23))*tan(rho*0.5_ark))
      !
      !b0(3,2,0) = 0
      !b0(1,2,0) = 0
      !b0(2,2,0) = 0
      !
      !b0(3,3,0) = (m2+m1)/m*re13*sin(rho-e)-m2/m*re23*sin(e)
      !b0(1,3,0) = -m3/m*re13*sin(rho-e)-m2/m*re23*sin(e)
      !b0(2,3,0) = (m1+m3)/m*re23*sin(e)-m3/m*re13*sin(rho-e)
      !
      !b0(3,1,0) =-(m2+m1)/m*re13*cos(rho-e)-m2/m*re23*cos(e)
      !b0(1,1,0) = m3/m*re13*cos(rho-e)-m2/m*re23*cos(e)
      !b0(2,1,0) = (m1+m3)/m*re23*cos(e)+m3/m*re13*cos(rho-e)


      select case(trim(molec%coords_transform))
      case default
         !
      case('R-RHO-Z','R-PHI-RHO-Z')
         !
         if (Nangles>0) then
           alphaeq = molec%alphaeq(1)
         elseif (molec%Ndihedrals>1) then 
           alphaeq = pi+(-molec%taueq(1)+molec%taueq(2)) !asin(sin(molec%taueq(1))*sqrt(2.0_ark))
         else
           alphaeq = pi-molec%taueq(1)
         endif
         alphae_h = alphaeq*0.5_ark
         !
         b0(1,3,0) = 0.0_ark
         b0(1,2,0) = 0.0_ark
         b0(1,1,0) = 2.0_ark*m3/m*re13*cos(alphae_h)
         !
         b0(2,3,0) =-re13*sin(alphae_h)
         b0(2,2,0) = 0.0_ark
         b0(2,1,0) =-m1/m*re13*cos(alphae_h)
         !
         b0(3,3,0) = re13*sin(alphae_h)
         b0(3,2,0) = 0.0_ark
         b0(3,1,0) =-m1/m*re13*cos(alphae_h)
         !
      case('R1-R2-Y+X')
         !
         if (Nangles>0) then
           alphaeq = molec%alphaeq(1)
         elseif (molec%Ndihedrals>1) then 
           alphaeq = pi-asin(sqrt(molec%taueq(1)**2+molec%taueq(2)**2))
         else
           alphaeq = pi-molec%taueq(1)
         endif
         !
         tau = pi*0.25_ark
         alphae_h = alphaeq*0.5_ark
         !
         b0(1,3,0) = 0.0_ark
         b0(1,2,0) = 2.0_ark*m3/m*re13*cos(alphae_h)*sin(tau)
         b0(1,1,0) = 2.0_ark*m3/m*re13*cos(alphae_h)*cos(tau)
         !
         b0(2,3,0) =-re13*sin(alphae_h)
         b0(2,2,0) =-m1/m*re13*cos(alphae_h)*sin(tau)
         b0(2,1,0) =-m1/m*re13*cos(alphae_h)*cos(tau)
         !
         b0(3,3,0) = re13*sin(alphae_h)
         b0(3,2,0) =-m1/m*re13*cos(alphae_h)*sin(tau)
         b0(3,1,0) =-m1/m*re13*cos(alphae_h)*cos(tau)
         !
      case('R-PHI1-PHI2-Z','R-PHI1-Z','R-S1-S2-Z','R-ALPHA-THETA-Z')
         !
         if (Nangles>0) then
           alphaeq = molec%alphaeq(1)
         elseif (molec%Ndihedrals>1) then
           !alphaeq = pi-(molec%taueq(1)+molec%taueq(2)) !asin(sin(molec%taueq(1))*sqrt(2.0_ark))
           alphaeq = pi-asin( sqrt( sin(molec%taueq(1))**2+sin(molec%taueq(2))**2 ))
         else
           alphaeq = pi-molec%taueq(1)
         endif 
         alphae_h = alphaeq*0.5_ark
         !
         b0(1,3,0) = 0.0_ark
         b0(1,2,0) = 0.0_ark
         b0(1,1,0) = 2.0_ark*m3/m*re13*cos(alphae_h)
         !
         b0(2,3,0) =-re13*sin(alphae_h)
         b0(2,2,0) = 0.0_ark
         b0(2,1,0) =-m1/m*re13*cos(alphae_h)
         !
         b0(3,3,0) = re13*sin(alphae_h)
         b0(3,2,0) = 0.0_ark
         b0(3,1,0) =-m1/m*re13*cos(alphae_h)
         !
         !a0(:,:) = b0(:,:,0)
         !
         !b0(:,1,0) = sqrt(0.5_ark)*( a0(:,1)+a0(:,2) )
         !b0(:,2,0) = sqrt(0.5_ark)*( a0(:,1)-a0(:,2) )
         !
      !case('R-R-R')
      !   !
      !   alphaeq = acos((re13**2+re13**2-re13**2)/(2.0_ark*re13*re13))
      !   !
      !   b0(1,3,0) = 0.0_ark
      !   b0(1,2,0) = 0.0_ark
      !   b0(1,1,0) = 2.0_ark*m3/m*re13*cos(alphaeq)
      !   !
      !   b0(2,3,0) =-re13*sin(alphaeq)
      !   b0(2,2,0) = 0.0_ark
      !   b0(2,1,0) =-m1/m*re13*cos(alphaeq)
      !   !
      !   b0(3,3,0) = re13*sin(alphaeq)
      !   b0(3,2,0) = 0.0_ark
      !   b0(3,1,0) =-m1/m*re13*cos(alphaeq)
      !   !
      end select
      !
      if (any(molec%AtomMasses(2:3)/=m2).or.re13/=re23) then
        !
        b0(1,1,0) = 0
        b0(1,2,0) = 0
        b0(1,3,0) = 0
        !
        b0(2,1,0) = 0
        b0(2,2,0) = 0
        b0(2,3,0) = re13
        !
        b0(3,1,0) = re23*sin(alphaeq)
        b0(3,2,0) = 0.0_ark
        b0(3,3,0) = re23*cos(alphaeq)
        !
        do i = 1,3
          CM_shift = sum(b0(:,i,0)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
          b0(:,i,0) = b0(:,i,0) - CM_shift
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

        !call MLorienting_a0(Natoms,molec%AtomMasses,b0(:,:,0))
        
        !write(out,"('ML_b0_XY2: masses-s are given in wrong order, must be M m m: ',3f14.6)") molec%AtomMasses(:)
        !stop 'ML_b0_XY2: ,masses are in wrong order'
        !
      endif 
      !
      !
      ! define the rho-type coordinate 
      !
      ! We can also simulate the reference structure at each point of rho, if npoints presented
      !
      if (Npoints/=0) then
         ! 
         !if (.not.present(rho_borders).or..not.present(rho_ref)) then  
         !   write(out,"('ML_b0_XY2: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
         !   stop 'ML_b0_XY2: rho_borders or rho_ref not specified '
         !endif
         !
         select case(trim(molec%coords_transform))
         case default
            write (out,"('ML_b0_XY2: coord. type ',a,' unknown')") trim(molec%coords_transform)
            stop 'ML_b0_XY2 - bad coord. type'
            !
         case('LINEAR')
            !
            rho_ref = alphaeq
            rho0 = rho_ref
            a02 = (m1/m)
            !
         case('DX-RHO')
            !
            rho_ref = alphaeq
            rho0 = 0.0
            a02 = (m1/m)
            !
         case('R-RHO','R-EXPRHO','R-RHO-Z','R12-R','R12-RHO','R13-RHO','R-PHI-RHO','R-PHI-RHO-Z','R-PHI1-PHI2-Z','R-PHI1-Z',&
              'R-S1-S2-Z','R-ALPHA-THETA-Z')
            !
            rho_ref = 0.0_ark
            rho0 = 0.0_ark
            a02 = (m1/m)
            !
         case('R-RHO-HALF')
            !
            rho_ref = 0.0_ark
            rho0 = 0.0_ark
            a02 = (m1/m)
            !
         case('RADAU')
            !
            rho_ref = 0.0_ark
            rho0 = 0.0_ark
            a02 = (m1/m)
            !
            a = sqrt(m3 / m)
            b = m2 / (m1+m2)
            g1 = 1.0_ark - a / (a+b-a*b)
            g2 = 1.0_ark - a / (1.0_ark-b+a*b)
            !
         end select
         !
         if ( rho0+rho_borders(1)<0.0.or.rho0+rho_borders(2)>pi) then  
            write(out,"('ML_b0_XY2: rho_borders exceed the [0,pi] interval: ',2f12.4)") (rho0+rho_borders(1:2))*rad
            !stop 'ML_b0_XY2: rho_borders exceed 0..pi '
         endif

         do i = 0,npoints
            !
            rho = rho0+rho_i(i)
            !
            alpha = rho
            !
            select case(trim(molec%coords_transform))
              case('R-RHO','R12-RHO','R13-RHO','R-RHO-Z','R-PHI-RHO','R-PHI-RHO-Z')
               alpha = pi-rho
              case('R-RHO-HALF')
               alpha = pi-rho*2.0_ark
            end select 
            !
            if(trim(molec%coords_transform)=='R-EXPRHO') alpha = pi-exp(rho)
            !
            !if(trim(molec%coords_transform)=='R-EXPRHO') alpha = pi-log(rho)
			!
			!if(trim(molec%coords_transform)=='R-EXPRHO') alpha = pi-(1.0_ark+rho**2)/(max(rho,small_))
			!
            if(trim(molec%coords_transform)=='RADAU') then 
              !
              !
              !f1= 1.0_ark/g1
              !f2= 1.0_ark/g2
              !f12= 1.0_ark - f1*f2
              !p1= r1*(1.0_ark-f1)/(g2*f12)
              !p2= r2*(1.0_ark-f2)/(g1*f12)
              !s1= r1-p1
              !s2= r2-p2
              !q1= sqrt(p1*p1 + s2*s2 + 2.0_ark*p1*s2*xcos)/(1.0_ark-g1)
              !q2= sqrt(p2*p2 + s1*s1 + 2.0_ark*p2*s1*xcos)/(1.0_ark-g2)
              !q3= sqrt(p1*p1 + p2*p2 - 2.0_ark*p1*p2*xcos)
              !cost = (q1*q1 + q2*q2 - q3*q3)/(2.0_ark*q1*q2)
              !theta = acos(cost)
              !
              if ( abs(rho-1.0_rk)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('ML_b0_XY2: cosrho = z-1>1: ',f18.8)") cosa
                 !stop 'ML_b0_XY2 - bad cosrho=z-1>1'
                 !
              elseif ( rho-1.0_ark>=1.0_ark) then 
                 alpha = 0.0_ark
              else 
                 alpha = acos(rho-1.0_ark)
              endif
              !
			  alpha = pi-rho
			  !
              !alpha = acos(rho-1.0_ark)
              !
              cosa = -(a02*cos(alpha)+cos(alpha)+1.0_ark-a02)/(-a02+a02*cos(alpha)-1.0_ark-cos(alpha))
              !
              if ( abs(cosa)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('ML_b0_XY2: cosalpha>1: ',f18.8)") cosa
                 stop 'ML_b0_XY2 - bad cosalpha'
                 !
              elseif ( cosa>=1.0_ark) then 
                 alpha = 0.0_ark
              else 
                 alpha = acos(cosa)
              endif
              !
              !
            endif 
            !
            if(trim(molec%coords_transform)=='R12-R') then
              !
              rho_ark = rho
              !
              r12 = ML_MEP_xy2_R12_R(rho_ark)
			  !
			  !r12 = sqrt(2.0_ark*molec%local_eq(1)**2-2.0_ark*molec%local_eq(1)**2*cos(molec%local_eq(3)))
              !
              alpha = atan2(r12*0.5_ark,rho)*2.0_ark
              !
              if (alpha>small_) re13 = rho/cos(alpha*0.5_ark)
              !
            endif
            !
            if(trim(molec%coords_transform)=='R12-RHO') then
              !
              alpha_ark = alpha
              !
              re13 = ML_MEP_xy2_R12_ALPHA(alpha_ark)
              !
            endif
            !
            if(trim(molec%coords_transform)=='R13-RHO') then
              !
              alpha_ark = alpha
              !
              r12 = ML_MEP_xy2_R13_ALPHA(alpha_ark)
              !
              re13 = 0.5_ark*r12/sin(alpha*0.5_ark)
              !
            endif
            !
            alphae_h = alpha*0.5_ark
            !
            b0(1,1,i) = 0.0_ark
            b0(1,2,i) = 0.0_ark
            b0(1,3,i) = 2.0_ark*m3/m*re13*cos(alphae_h)
            !
            b0(2,1,i) =-re13*sin(alphae_h)
            b0(2,2,i) = 0.0_ark
            b0(2,3,i) =-m1/m*re13*cos(alphae_h)
            !
            b0(3,1,i) = re13*sin(alphae_h)
            b0(3,2,i) = 0.0_ark
            b0(3,3,i) =-m1/m*re13*cos(alphae_h)
            !
            select case(trim(molec%coords_transform))
               !
            case('R-RHO-Z','R-PHI-RHO-Z')
               !
               b0(1,3,i) = 0.0_ark
               b0(1,2,i) = 0.0_ark
               b0(1,1,i) = 2.0_ark*m3/m*re13*cos(alphae_h)
               !
               b0(2,3,i) =-re13*sin(alphae_h)
               b0(2,2,i) = 0.0_ark
               b0(2,1,i) =-m1/m*re13*cos(alphae_h)
               !
               b0(3,3,i) = re13*sin(alphae_h)
               b0(3,2,i) = 0.0_ark
               b0(3,1,i) =-m1/m*re13*cos(alphae_h)
               !
               !rho_ark = rho
               !
               !rot = ML_euler_rotait(rho_ark,0.0_ark,0.0_ark)
               !do ix = 1,Natoms
               !  b0(ix,:,i) = matmul(transpose(rot),b0(ix,:,i))
               !enddo
               !
            case('R-PHI1-PHI2-Z','R-PHI1-Z','R-S1-S2-Z')
               !
               if (Nangles>0) then
                 alphaeq = molec%alphaeq(1)
               elseif (molec%Ndihedrals>1) then
                 !alphaeq = pi-(molec%taueq(1)+molec%taueq(2)) !asin(sin(molec%taueq(1))*sqrt(2.0_ark))
                 alphaeq = pi-asin( sqrt( sin(rho)**2+sin(rho)**2 ))
               else
                 alphaeq = pi-rho
               endif 
               alphae_h = alphaeq*0.5_ark
               !
               b0(1,3,i) = 0.0_ark
               b0(1,2,i) = 0.0_ark
               b0(1,1,i) = 2.0_ark*m3/m*re13*cos(alphae_h)
               !
               b0(2,3,i) =-re13*sin(alphae_h)
               b0(2,2,i) = 0.0_ark
               b0(2,1,i) =-m1/m*re13*cos(alphae_h)
               !
               b0(3,3,i) = re13*sin(alphae_h)
               b0(3,2,i) = 0.0_ark
               b0(3,1,i) =-m1/m*re13*cos(alphae_h)
               !
            case('R-ALPHA-THETA-Z')
               !
               if (Nangles>0) then
                 alphaeq = molec%alphaeq(1)
               elseif (molec%Ndihedrals>1) then
                 !alphaeq = pi-(molec%taueq(1)+molec%taueq(2)) !asin(sin(molec%taueq(1))*sqrt(2.0_ark))
                 !
                 tau = rho
                 !
                 alphaeq = pi-asin( sqrt( sin(tau)**2+sin(tau)**2 ))
                 !
               else
                 alphaeq = pi-rho
               endif 
               alphae_h = alphaeq*0.5_ark
               !
               b0(1,3,i) = 0.0_ark
               b0(1,2,i) = 0.0_ark
               b0(1,1,i) = 2.0_ark*m3/m*re13*cos(alphae_h)
               !
               b0(2,3,i) =-re13*sin(alphae_h)
               b0(2,2,i) = 0.0_ark
               b0(2,1,i) =-m1/m*re13*cos(alphae_h)
               !
               b0(3,3,i) = re13*sin(alphae_h)
               b0(3,2,i) = 0.0_ark
               b0(3,1,i) =-m1/m*re13*cos(alphae_h)
               !
            end select
            !
            !if (any(molec%AtomMasses(2:3)/=m2)) then
            !  !
            !  b0(1,1,i) = 0.0_ark
            !  b0(1,2,i) = 0.0_ark
            !  b0(1,3,i) = 0.0_ark
            !  !
            !  b0(2,1,i) =-re13*cos(alphae_h)
            !  b0(2,2,i) = 0.0_ark
            !  b0(2,3,i) =-re13*sin(alphae_h)
            !  !
            !  b0(3,1,i) =-re13*cos(alphae_h)
            !  b0(3,2,i) = 0.0_ark
            !  b0(3,3,i) = re13*sin(alphae_h)
            !  !
            !  do ix = 1,3
            !    CM_shift = sum(b0(:,ix,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
            !    b0(:,ix,i) = b0(:,ix,i) - CM_shift
            !  enddo 
            !  !
            !  ijk(1,1) = 2 ; ijk(1,2) = 3
            !  ijk(2,1) = 3 ; ijk(2,2) = 1
            !  ijk(3,1) = 1 ; ijk(3,2) = 2
            !  !
            !  method = 'ULEN'
            !  !
            !  a0 = b0(:,:,i)
            !  !
            !  call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i),transform,method=method)
            !  !
            !  loop_xyz : do ix =1,3
            !    !
            !    do ioper = 1,4
            !      !
            !      phi = real(ioper-1,ark)*pi*0.5_ark
            !      !
            !      tmat = 0 
            !      !
            !      tmat(ix,ix) = 1.0_ark
            !      !
            !      tmat(ijk(ix,1),ijk(ix,1)) = cos(phi)
            !      tmat(ijk(ix,1),ijk(ix,2)) = sin(phi)
            !      tmat(ijk(ix,2),ijk(ix,1)) =-sin(phi)
            !      tmat(ijk(ix,2),ijk(ix,2)) = cos(phi)
            !      !
            !      transform1 = matmul(tmat,transform)
            !      !
            !      forall(ix=1:3) unit(ix) = dot_product(transform0(ix,:),transform1(ix,:))
            !      !
            !      unit = unit - 1.0_ark
            !      !
            !      if ( all( abs( unit(:) )<1.0e-1  ) ) exit loop_xyz
            !      !
            !    enddo
            !  enddo loop_xyz
            !  !
            !  transform0 = matmul(tmat,transform)
            !  !
            !  forall(ix=1:3) b0(ix,:,i) = matmul(transpose(tmat),b0(ix,:,i))
            !  !
            !endif
            !
            if (any(molec%AtomMasses(2:3)/=m2)) then
              !
              b0(1,1,0) = 0
              b0(1,2,0) = 0
              b0(1,3,0) = 0
              !
              b0(2,1,0) = 0
              b0(2,2,0) = 0
              b0(2,3,0) = re13
              !
              b0(3,1,0) = re23*sin(alpha)
              b0(3,2,0) = 0.0_ark
              b0(3,3,0) = re23*cos(alpha)
              !
              !b0(1,1,i) = re13*cos(alphae_h)*(m2+m3)/(m1+m2+m3)
              !b0(1,2,i) = 0.0_ark
              !b0(1,3,i) = re13*sin(alphae_h)*(m2-m3)/(m1+m2+m3)
              !
              !b0(2,1,i) = -re13*cos(alphae_h)*m1/(m1+m2+m3)
              !b0(2,2,i) = 0.0_ark
              !b0(2,3,i) = -re13*sin(alphae_h)*(m1+2.0_ark*m3)/(m1+m2+m3)
              !
              !b0(3,1,i) = -re13*cos(alphae_h)*m1/(m1+m2+m3)
              !b0(3,2,i) = 0.0_ark
              !b0(3,3,i) = re13*sin(alphae_h)*(m1+2.0_ark*m2)/(m1+m2+m3)

              do ix = 1,3
                CM_shift = sum(b0(:,ix,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
                b0(:,ix,i) = b0(:,ix,i) - CM_shift
              enddo 
              !
              !tanphi = -0.5_ark*(m1*sin(alpha)*(m2-m3))/(m2*m1*cos(alpha)-2*m2*m3+2*m2*m3*cos(alpha)+m1*m3*cos(alpha))
              !
              !phi = atan(tanphi)
              !
              !print *,i,phi
              !
              !phi =-0.25_ark*atan(m1*sin(alpha)*(m2-m3)/(2.0_ark*m2*m3*cos(alpha)-2.0_ark*m2*m3+m1*m3*cos(alpha)+m1*m2*cos(alpha)))
              !
              !tmat = 0 ; tmat(2,2) = 1.0_ark
              !
              !tmat(1,1) = cos(phi)
              !tmat(1,3) =-sin(phi)
              !tmat(3,1) = sin(phi)
              !tmat(3,3) = cos(phi)
              !
              !forall(ix=1:3) b0(ix,:,i) = matmul((tmat),b0(ix,:,i))
              !
              !transform1 = matmul(tmat,transform)
              !
              !do ix = 1,3
              !  CM_shift = sum(b0(:,ix,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
              !  b0(:,ix,i) = b0(:,ix,i) - CM_shift
              !enddo 
              !
              call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i),transform)
              ! 
              do ix = 1,3 
                 do jx = 1,3 
                    if (ix/=jx) then  
                       !
                       a_t =  sum(molec%AtomMasses(:)*b0(:,ix,i)*b0(:,jx,i) )
                       !
                       if (abs(a_t)>100.0_rk*(small_)) then 
                           write(out,"('ML_b0_XY3: b0 is not a solution of Eckart 2 for ix,jx =',3i4,d18.8)") ix,jx,i,a_t
                           stop 'ML_b0_XY3: b0 is not solution of Eckart2'
                       endif
                    endif
                 enddo
                 !
              enddo
              !
            endif
            !
            !if (verbose>=4) then 
            !  write(out,"(i6)") molec%natoms
            !  !
            !  write(out,"(/'O',3x,3f14.8)") b0(1,:,i)*bohr
            !  write(out,"( 'D',3x,3f14.8)") b0(2,:,i)*bohr
            !  write(out,"( 'H',3x,3f14.8)") b0(3,:,i)*bohr
            !  !
            !endif
            !
         enddo
         !
      else
         !
         !if (verbose>=4) then 
         !  write(out,"(i6)") molec%natoms
         !  !
         !  write(out,"(/'O',3x,3f14.8)") b0(1,:,0)
         !  write(out,"( 'H',3x,3f14.8)") b0(2,:,0)
         !  write(out,"( 'H',3x,3f14.8)") b0(3,:,0)
         !  !
         !endif
         !
      endif
      !
      ! For special angles of a linear molecule
      !
      if (Nangles==0) then
         !
         do i = 0,npoints
            !
            !phi = pi*0.25_ark
            !
            !tmat = 0 
            !
            !tmat(1,1) = 1.0_ark
            !tmat(2,2) = 1.0_ark
            !tmat(3,3) = 1.0_ark
            !
            !tmat(2,2) = cos(phi)
            !tmat(2,3) = sin(phi)
            !tmat(3,2) =-sin(phi)
            !tmat(3,3) = cos(phi)
            !
            !do iatom = 1,molec%natoms
            !  b0(iatom,:,i) = matmul(tmat,b0(iatom,:,i))
            !enddo 
            !
            if (verbose>=4) then 
              write(out,"(i6)") molec%natoms
              !
              write(out,"(/'O',3x,3f14.8)") b0(1,:,i) !*bohr
              write(out,"( 'H',3x,3f14.8)") b0(2,:,i) !*bohr
              write(out,"( 'H',3x,3f14.8)") b0(3,:,i) !*bohr
              !
            endif
            !
         enddo
         !
      endif
      !
      if (verbose>=4) write(out,"('ML_b0_XY2/end')") 

  end subroutine ML_b0_XY2
  !
  !
  ! Here we define the symmetry transformation of the Nmodes coordinates according the symmetry operations
  !
  subroutine ML_symmetry_transformation_XY2(ioper,natoms,src,dst)
    !
    integer(ik),intent(in)    :: ioper  ! group operation  
    integer(ik),intent(in)    :: natoms
    real(ark),intent(in)      :: src(natoms)
    real(ark),intent(out)     :: dst(natoms)
    !
    integer(ik) :: nsrc,Nrot,irot,N_Cn,ioper_,NC2
    real(ark)   :: phi,q1,q2,r,theta,phi_n,a,b,e,o,repres(sym%Noper,2,2),qx,qy
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY2/start')") 
    !
    nsrc = size(src)
    !
    select case(trim(molec%symmetry))
      !
    case('C','C(M)')
      !
      dst = src
      !
      return
      !
    end select
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('ML_symmetry_transformation_XY2. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'ML_coordinate_transform_XY2 - bad coord. type'
       !
    case('R-RHO','R-EXPRHO','RADAU','R-RHO-Z','R12-RHO','R13-RHO','R-PHI1','R-PHI1-Z','R-RHO-HALF')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_XY2 - bad symm. type'
          !
       case('C2V','C2V(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)

          case (3) ! (E*)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)

          case (4) ! (12)*

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       case('C2VN')
          !
          irot = mod(ioper-1,4)+1
          !
          select case(irot)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)

          case (3) ! (E*)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)

          case (4) ! (12)*

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)

          case default
            write (out,"('ML_symmetry_transformation_XY2: irot ',i8,' unknown')") irot
            stop 'ML_symmetry_transformation_XY2 - bad irot. type'
 
          end select 
          !
       case('CS','CS(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       case('DNH','DNH(M)')
          !
          ! Number of eq. rotations 
          Nrot = sym%N
          !
          ! Number of Cn classes 
          N_Cn = sym%N/2
          !
          ! for odd Dnh groups
          !
          if (mod(sym%N,2)==1) then
              stop 'ML_symmetry_transformation_XY2: DNH odd for R1-R2-Y+X is not implemented' 
          endif
          !
          !qx= 1
          !qy= 1
          !q2x= 2
          !q2y= -2
          !
          NC2 = sym%N/2
          N_Cn = sym%N/2-1
          !
          if (ioper==1) then ! E 
            !
            dst = src
            !
          elseif (ioper<=1+2*N_Cn+1) then ! Cn x 2 x(n/2-1), C2 only once
            !
            ioper_ =ioper-1 
            irot = (ioper_+1)/2
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            !
          elseif (ioper<=2+2*N_Cn+NC2) then !  C'2
            !
            irot =ioper-(2+2*N_Cn)-1
            !
            phi_n = phi*irot*2.0_ark
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            !
          elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
            !
            irot =ioper-(2+2*N_Cn+NC2)-1
            !
            !phi_n =(irot*2.0_ark-1.0_ark)*phi*0.5_ark
            !
            phi_n = phi*(2*irot+1)
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            !
          elseif (ioper==3+2*N_Cn+2*NC2) then ! i
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            !
          elseif (ioper<=3+4*N_Cn+2*NC2) then !  2xnxSn
            !
            ioper_ =ioper-(3+2*N_Cn+2*NC2)
            irot = (ioper_+1)/2
            !
            phi_n = phi*irot
            !
            ! Second oper in a class is with negative phi
            if (mod(ioper_,2)==0)  phi_n = -phi_n
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            !
          elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigmah
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            !
          elseif (ioper<=4+4*N_Cn+3*NC2) then !  sigmav
            !
            irot = ioper-(4+4*N_Cn+2*NC2)-1
            !
            phi_n = phi*irot*2.0_ark
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            !
          elseif (ioper<=4+4*N_Cn+4*NC2) then !  sigmad
            !
            irot = ioper-(4+4*N_Cn+3*NC2)-1
            !
            phi_n = (2*irot+1)*phi
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            !
          else
            !
            write (out,"('ML_symmetry_transformation_abcd  in Dinfty: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_abcd Dinfty - bad operation. type'
            !         
          endif 
          !
       end select
       !
    case('R12-R')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_XY2 - bad symm. type'
          !
       case('C2V','C2V(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(1)
            dst(2) =-src(2)
            dst(3) = src(3)

          case (3) ! (E*)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)

          case (4) ! (12)*

            dst(1) = src(1)
            dst(2) =-src(2)
            dst(3) = src(3)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       end select
       !
    case('R-R-R')
       !
       select case(trim(molec%symmetry))
       !
       case default
          write (out,"('ML_symmetry_transformation_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_XY2 - bad symm. type'
          !
       case('CS','CS(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       case('D3H','D3H(M)')
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

             write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
             stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
           end select 
          !
       end select
       !
    case('R1-R2-Y+X')
       !
       select case(trim(molec%symmetry))
         !
       case default
         write (out,"('ML_symmetry_transformation_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
         stop 'ML_symmetry_transformation_XY2 - bad symm. type'
         !
       case('DNH','DNH(M)')
          !
          ! Number of eq. rotations 
          Nrot = sym%N
          !
          ! Number of Cn classes 
          N_Cn = sym%N/2
          !
          phi = 2.0_ark*pi/real(Nrot,ark)
          !
          qx= src(3)
          qy= src(4)
          !
          ! for odd Dnh groups
          !
          if (mod(sym%N,2)==1) then
              stop 'ML_symmetry_transformation_XY2: DNH odd for R1-R2-Y+X is not implemented' 
          endif
          !
          !qx= 1
          !qy= 1
          !q2x= 2
          !q2y= -2
          !
          NC2 = sym%N/2
          N_Cn = sym%N/2-1
          !
          if (ioper==1) then ! E 
            !
            dst = src
            !
          elseif (ioper<=1+2*N_Cn+1) then ! Cn x 2 x(n/2-1), C2 only once
            !
            ioper_ =ioper-1 
            irot = (ioper_+1)/2
            !
            phi_n = phi*irot
            !
            ! Second oper in a class is with negative phi
            if ( mod(ioper_,2)==0 ) phi_n = -phi_n
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = qx*cos(phi_n)-qy*sin(phi_n)
            dst(4) = qx*sin(phi_n)+qy*cos(phi_n)
            !
          elseif (ioper<=2+2*N_Cn+NC2) then !  C'2
            !
            irot =ioper-(2+2*N_Cn)-1
            !
            phi_n = phi*irot*2.0_ark
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = qx*cos(phi_n)+qy*sin(phi_n)
            dst(4) = qx*sin(phi_n)-qy*cos(phi_n)
            !
          elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
            !
            irot =ioper-(2+2*N_Cn+NC2)-1
            !
            !phi_n =(irot*2.0_ark-1.0_ark)*phi*0.5_ark
            !
            phi_n = phi*(2*irot+1)
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = qx*cos(phi_n)+qy*sin(phi_n)
            dst(4) = qx*sin(phi_n)-qy*cos(phi_n)
            !
          elseif (ioper==3+2*N_Cn+2*NC2) then ! i
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-qx
            dst(4) =-qy
            !
          elseif (ioper<=3+4*N_Cn+2*NC2) then !  2xnxSn
            !
            ioper_ =ioper-(3+2*N_Cn+2*NC2)
            irot = (ioper_+1)/2
            !
            phi_n = phi*irot
            !
            ! Second oper in a class is with negative phi
            if (mod(ioper_,2)==0)  phi_n = -phi_n
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =  (qx*cos(phi_n)-qy*sin(phi_n))
            dst(4) =  (qx*sin(phi_n)+qy*cos(phi_n))
            !
          elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigmah
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = qx
            dst(4) = qy
            !
          elseif (ioper<=4+4*N_Cn+3*NC2) then !  sigmav
            !
            irot = ioper-(4+4*N_Cn+2*NC2)-1
            !
            phi_n = phi*irot*2.0_ark
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) =  (qx*cos(phi_n)+qy*sin(phi_n))
            dst(4) =  (qx*sin(phi_n)-qy*cos(phi_n))
            !
          elseif (ioper<=4+4*N_Cn+4*NC2) then !  sigmad
            !
            irot = ioper-(4+4*N_Cn+3*NC2)-1
            !
            phi_n = (2*irot+1)*phi
            !
            !irot = ioper-(4+4*N_Cn+3*NC2)
            !phi_n = (-phi*0.5+phi*irot)*2.0_ark
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) =  (qx*cos(phi_n)+qy*sin(phi_n))
            dst(4) =  (qx*sin(phi_n)-qy*cos(phi_n))
            !
          else
            !
            write (out,"('ML_symmetry_transformation_abcd  in Dinfty: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_abcd Dinfty - bad operation. type'
            !         
          endif
          !
       end select 
       !
    case('R-PHI-RHO-Z','R-PHI-RHO')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_XY2 - bad symm. type'
          !
       case('C2V','C2V(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            dst(4) = src(4)

          case (3) ! (E*)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            dst(4) = src(4)

          case (4) ! (12)*

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            dst(4) = src(4)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       case('CS','CS(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            dst(4) = src(4)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       end select
          !
       !
    case('R-PHI1-PHI2','R-PHI1-PHI2-Z')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_XY2 - bad symm. type'
          !
       case('C2V','C2V(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)
            dst(4) =-src(4)

          case (3) ! (E*)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) =-src(4)
            dst(4) =-src(3)

          case (4) ! (12)*

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(4)
            dst(4) = src(3)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       case('C2H','C2H(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (C2c)

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)
            dst(4) = src(4)

          case (3) ! (sigma-ab)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            dst(4) =-src(4)

          case (4) ! (i)

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)
            dst(4) =-src(4)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select           
          !
       case('D2H','D2H(M)')
          !
          select case(ioper)
            !
          case (1) ! E 
            !
            dst = src
            !
          case (2) ! (C2a) -> Ra
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) =-src(3)
            dst(4) =-src(4)
            !
          case (3) ! (C2b)
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            dst(4) =-src(4)
            !
          case (4) ! (C2c)
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)
            dst(4) = src(4)
            !
          case (5) ! E* -> sigma_ab
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            dst(4) =-src(4)
            !
          case (6) ! sigma_ac
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) =-src(3)
            dst(4) = src(4)
            !
          case (7) ! sigma_bc
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            dst(4) = src(4)
            !
          case (8) ! i
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)
            dst(4) =-src(4)
            
          case default
            !
            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
            !
          end select
          !
       case('D3H','D3H(M)')
          !
          a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
          !
          repres ( 1,:,:)= sym%irr(6, 1)%repres
          repres ( 2,:,:)= sym%irr(6, 2)%repres
          repres ( 3,:,:)= sym%irr(6, 3)%repres
          repres ( 4,:,:)= sym%irr(6, 4)%repres
          repres ( 5,:,:)= sym%irr(6, 5)%repres
          repres ( 6,:,:)= sym%irr(6, 6)%repres
          repres ( 7,:,:)= sym%irr(6, 1)%repres
          repres ( 8,:,:)= sym%irr(6, 2)%repres
          repres ( 9,:,:)= sym%irr(6, 3)%repres
          repres (10,:,:)= sym%irr(6, 4)%repres
          repres (11,:,:)= sym%irr(6, 5)%repres
          repres (12,:,:)= sym%irr(6, 6)%repres
          !
          dst(3) = repres(ioper,1,1)*src(3)+repres(ioper,1,2)*src(4)
          dst(4) = repres(ioper,2,1)*src(3)+repres(ioper,2,2)*src(4)
          !
          select case(ioper)
          !
          case (1:6) ! identity,(23)*

            dst(1:2) = src(1:2)
            !
          case (7:12) ! (123),(123)*
            !
            dst(1) = src(2)
            dst(2) = src(1)
            !
          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       case('CS','CS(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            dst(4) = src(4)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       case('DINFTYH','DINFTYH(M)')
          !
          Nrot = sym%Nelements(3)
          !
          phi = 2.0_ark*pi/real(Nrot,ark)
          !
          q1 = src(3)
          q2 = src(4)
          !
          r = sqrt(q1**2+q2**2)
          theta = atan2(q2,q1)
          !
          if (ioper==1) then ! E 
            !
            dst = src
            !
          elseif (ioper==2) then !  Cinf
            !
            !irot = (ioper-1)
            phi_n = phi
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = q1*cos(phi_n)-q2*sin(phi_n)
            dst(4) = q1*sin(phi_n)+q2*cos(phi_n)
            !
          elseif (ioper==3) then !  2*Cinf
            !
            !irot = ioper-(1+Nrot)
            !
            phi_n = phi
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = q1*cos(-phi_n)-q2*sin(-phi_n)
            dst(4) = q1*sin(-phi_n)+q2*cos(-phi_n)
            !
          elseif (ioper<=3+Nrot) then !  sigmav
            !
            irot = ioper-3
            !
            phi_n = phi*irot
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = q1*cos(phi_n)+q2*sin(phi_n)
            dst(4) =-q1*sin(phi_n)+q2*cos(phi_n)
            !
          elseif (ioper==3+Nrot+1) then ! i
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)
            dst(4) =-src(4)
            !
          elseif (ioper==3+Nrot+2) then !  Sinf
            !
            !irot = ioper-(2+4*Nrot)
            !
            phi_n = phi
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = q1*cos(phi_n)-q2*sin(phi_n)
            dst(4) = q1*sin(phi_n)+q2*cos(phi_n)
            !
          elseif (ioper==3+Nrot+3) then !  2 Sinf
            !
            !irot = ioper-(2+5*Nrot)
            phi_n = phi
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = q1*cos(-phi_n)-q2*sin(-phi_n)
            dst(4) = q1*sin(-phi_n)+q2*cos(-phi_n)
            !
          elseif (ioper<=3+Nrot+3+Nrot) then !  C'2
            !
            irot = ioper-(3+Nrot+3)
            !
            phi_n = phi*irot
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = q1*cos(phi_n)+q2*sin(phi_n)
            dst(4) =-q1*sin(phi_n)+q2*cos(phi_n)
            !
          else
            !
            write (out,"('ML_symmetry_transformation_XY2  in Dinfty: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 Dinfty - bad operation. type'
 
          endif 
          !
       case('DNH','DNH(M)')
          !
          ! Number of eq. rotations 
          Nrot = sym%N
          !
          ! Number of Cn classes 
          N_Cn = sym%N/2
          !
          phi = 2.0_ark*pi/real(Nrot,ark)
          !
          q1 = src(3)
          q2 = src(4)
          !
          r = sqrt(q1**2+q2**2)
          theta = atan2(q2,q1)
          !
          if (mod(sym%N,2)==1) then 
             !
             if (ioper==1) then ! E 
               !
               dst = src
               !
             elseif (ioper<=1+2*N_Cn) then !  Cinf
               !
               ioper_ =ioper-1 
               irot = (ioper_+1)/2
               !
               phi_n = phi*irot
               !
               ! Second oper in a class is with negative phi
               if ( mod(ioper_,2)==0 ) phi_n = -phi_n
               !
               dst(1) = src(1)
               dst(2) = src(2)
               dst(3) = q1*cos(phi_n)-q2*sin(phi_n)
               dst(4) = q1*sin(phi_n)+q2*cos(phi_n)
               !
             elseif (ioper<=1+2*N_Cn+Nrot) then !  C'2
               !
               irot = ioper-(1+2*N_Cn)
               !
               phi_n = phi*irot
               !
               dst(1) = src(2)
               dst(2) = src(1)
               dst(3) = q1*cos(phi_n)+q2*sin(phi_n)
               dst(4) = q1*sin(phi_n)-q2*cos(phi_n)
               !
             elseif (ioper==1+2*N_Cn+Nrot+1) then ! sigmah
               !
               dst(1) = src(2)
               dst(2) = src(1)
               dst(3) = src(3)
               dst(4) = src(4)
               !
             elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn) then !  Sinf
               !
               ioper_ = ioper-(1+2*N_Cn+Nrot+1)
               !
               irot = (ioper_+1)/2
               !
               phi_n = phi*irot
               !
               ! Second oper in a class is with negative phi
               if (mod(ioper_,2)==0)  phi_n = -phi_n
               !
               dst(1) = src(2)
               dst(2) = src(1)
               dst(3) = (q1*cos(phi_n)-q2*sin(phi_n))
               dst(4) = (q1*sin(phi_n)+q2*cos(phi_n))
               !
             elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then !  C'2
               !
               irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
               !
               phi_n = phi*irot
               !
               dst(1) = src(1)
               dst(2) = src(2)
               dst(3) = ( q1*cos(phi_n)+q2*sin(phi_n))
               dst(4) = ( q1*sin(phi_n)-q2*cos(phi_n))
               !
             else
               !
               write (out,"('ML_symmetry_transformation_XY2  in Dinfty: operation ',i8,' unknown')") ioper
               stop 'ML_symmetry_transformation_XY2 Dinfty - bad operation. type'
         
             endif 
             !
             !do ioper_=1,sym%Noper/2
             !  repres ( ioper_,            :,:)= sym%irr(5, ioper_)%repres
             !  repres ( ioper_+sym%Noper/2,:,:)= sym%irr(5, ioper_)%repres
             !enddo
             !
             do ioper_=1,sym%Noper
               repres ( ioper_,            :,:)= sym%irr(5, ioper_)%repres
               !repres ( ioper_+sym%Noper/2,:,:)= sym%irr(5, ioper_)%repres
             enddo
             !
             !dst(3) = repres(ioper,1,1)*src(3)+repres(ioper,1,2)*src(4)
             !dst(4) = repres(ioper,2,1)*src(3)+repres(ioper,2,2)*src(4)
             !
          else ! even n 
             !
             NC2 = sym%N/2
             N_Cn = sym%N/2-1
             !
             qx= src(3)
             qy= src(4)
             !
             if (ioper==1) then ! E 
               !
               dst = src
               !
             elseif (ioper<=1+2*N_Cn+1) then ! Cn x 2 x(n/2-1), C2 only once
               !
               ioper_ =ioper-1 
               irot = (ioper_+1)/2
               !
               phi_n = phi*irot
               !
               ! Second oper in a class is with negative phi
               if ( mod(ioper_,2)==0 ) phi_n = -phi_n
               !
               dst(1) = src(1)
               dst(2) = src(2)
               dst(3) = qx*cos(phi_n)-qy*sin(phi_n)
               dst(4) = qx*sin(phi_n)+qy*cos(phi_n)
               !
             elseif (ioper<=2+2*N_Cn+NC2) then !  C'2
               !
               irot =ioper-(2+2*N_Cn)-1
               !
               phi_n = phi*irot*2.0_ark
               !
               dst(1) = src(2)
               dst(2) = src(1)
               dst(3) = qx*cos(phi_n)+qy*sin(phi_n)
               dst(4) = qx*sin(phi_n)-qy*cos(phi_n)
               !
             elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
               !
               irot =ioper-(2+2*N_Cn+NC2)-1
               !
               !phi_n =(irot*2.0_ark-1.0_ark)*phi*0.5_ark
               !
               phi_n = phi*(2*irot+1)
               !
               dst(1) = src(2)
               dst(2) = src(1)
               dst(3) = qx*cos(phi_n)+qy*sin(phi_n)
               dst(4) = qx*sin(phi_n)-qy*cos(phi_n)
               !
             elseif (ioper==3+2*N_Cn+2*NC2) then ! i
               !
               dst(1) = src(2)
               dst(2) = src(1)
               dst(3) =-qx
               dst(4) =-qy
               !
             elseif (ioper<=3+4*N_Cn+2*NC2) then !  2xnxSn
               !
               ioper_ =ioper-(3+2*N_Cn+2*NC2)
               irot = (ioper_+1)/2
               !
               phi_n = phi*irot
               !
               ! Second oper in a class is with negative phi
               if (mod(ioper_,2)==0)  phi_n = -phi_n
               !
               dst(1) = src(2)
               dst(2) = src(1)
               dst(3) =  (qx*cos(phi_n)-qy*sin(phi_n))
               dst(4) =  (qx*sin(phi_n)+qy*cos(phi_n))
               !
             elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigmah
               !
               dst(1) = src(2)
               dst(2) = src(1)
               dst(3) = qx
               dst(4) = qy
               !
             elseif (ioper<=4+4*N_Cn+3*NC2) then !  sigmav
               !
               irot = ioper-(4+4*N_Cn+2*NC2)-1
               !
               phi_n = phi*irot*2.0_ark
               !
               dst(1) = src(1)
               dst(2) = src(2)
               dst(3) =  (qx*cos(phi_n)+qy*sin(phi_n))
               dst(4) =  (qx*sin(phi_n)-qy*cos(phi_n))
               !
             elseif (ioper<=4+4*N_Cn+4*NC2) then !  sigmad
               !
               irot = ioper-(4+4*N_Cn+3*NC2)-1
               !
               phi_n = (2*irot+1)*phi
               !
               !irot = ioper-(4+4*N_Cn+3*NC2)
               !phi_n = (-phi*0.5+phi*irot)*2.0_ark
               !
               dst(1) = src(1)
               dst(2) = src(2)
               dst(3) =  (qx*cos(phi_n)+qy*sin(phi_n))
               dst(4) =  (qx*sin(phi_n)-qy*cos(phi_n))
               !
             else
               !
               write (out,"('ML_symmetry_transformation_abcd  in Dinfty: operation ',i8,' unknown')") ioper
               stop 'ML_symmetry_transformation_abcd Dinfty - bad operation. type'
         
             endif 
             !
          endif
          !
       end select 
       !
    case('R-ALPHA-THETA-Z')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_XY2 - bad symm. type'
          !
       case('DNH','DNH(M)')
          !
          ! Number of eq. rotations 
          Nrot = sym%N
          !
          ! Number of Cn classes 
          N_Cn = sym%N/2
          !
          phi = 2.0_ark*pi/real(Nrot,ark)
          !
          r = src(4)
          theta = src(3)
          !
          if (ioper==1) then ! E 
            !
            dst = src
            !
          elseif (ioper<=1+2*N_Cn) then !  Cinf
            !
            ioper_ =ioper-1 
            irot = (ioper_+1)/2
            !
            phi_n = phi*irot
            !
            ! Second oper in a class is with negative phi
            if ( mod(ioper_,2)==0 ) phi_n = -phi_n
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)+phi_n
            dst(4) = src(4)
            !
          elseif (ioper<=1+2*N_Cn+Nrot) then !  C'2
            !
            irot = ioper-(1+2*N_Cn)
            !
            phi_n = phi*irot
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)+phi_n
            dst(4) = src(4)
            !
          elseif (ioper==1+2*N_Cn+Nrot+1) then ! sigmah
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            dst(4) = src(4)
            !
          elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn) then !  Sinf
            !
            ioper_ = ioper-(1+2*N_Cn+Nrot+1)
            !
            irot = (ioper_+1)/2
            !
            phi_n = phi*irot
            !
            ! Second oper in a class is with negative phi
            if (mod(ioper_,2)==0)  phi_n = -phi_n
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)+phi_n
            dst(4) = src(4)
            !
          elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then !  C'2
            !
            irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
            !
            phi_n = phi*irot
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) =-src(3)+phi_n
            dst(4) = src(4)
            !
          else
            !
            write (out,"('ML_symmetry_transformation_XY2  in Dinfty: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 Dinfty - bad operation. type'
 
          endif 
          !
          !do ioper_=1,sym%Noper/2
          !  repres ( ioper_,            :,:)= sym%irr(5, ioper_)%repres
          !  repres ( ioper_+sym%Noper/2,:,:)= sym%irr(5, ioper_)%repres
          !enddo
          !
          do ioper_=1,sym%Noper
            repres ( ioper_,            :,:)= sym%irr(5, ioper_)%repres
            !repres ( ioper_+sym%Noper/2,:,:)= sym%irr(5, ioper_)%repres
          enddo
          !
          dst(3) = repres(ioper,1,1)*src(3)+repres(ioper,1,2)*src(4)
          dst(4) = repres(ioper,2,1)*src(3)+repres(ioper,2,2)*src(4)
          !
       end select        
       !
    case('R-S1-S2','R-S1-S2-Z')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_XY2 - bad symm. type'
          !
       case('C2V','C2V(M)')
          !
          select case(ioper)

          case (1) ! E 

            dst = src

          case (2) ! (12)

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)
            dst(4) =-src(4)

          case (3) ! (E*)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            dst(4) =-src(4)

          case (4) ! (12)*

            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) =-src(3)
            dst(4) = src(4)

          case default

            write (out,"('ML_symmetry_transformation_XY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 - bad operation. type'
 
          end select 
          !
       case('DNH','DNH(M)')
          !
          ! Number of eq. rotations 
          Nrot = sym%N
          !
          ! Number of Cn classes 
          N_Cn = sym%N/2
          !
          phi = 2.0_ark*pi/real(Nrot,ark)
          !
          !q1 = src(3)
          !q2 = src(4)
          !
          q1 = (src(3)+src(4))/sqrt(2.0_ark)
          q2 = (src(4)-src(3))/sqrt(2.0_ark)
          !
          r = sqrt(q1**2+q2**2)
          theta = atan2(q2,q1)
          !
          if (ioper==1) then ! E 
            !
            dst = src
            !
          elseif (ioper<=1+2*N_Cn) then !  Cinf
            !
            ioper_ =ioper-1 
            irot = (ioper_+1)/2
            !
            phi_n = phi*irot
            !
            ! Second oper in a class is with negative phi
            if ( mod(ioper_,2)==0 ) phi_n = -phi_n
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = q1*cos(phi_n)-q2*sin(phi_n)
            dst(4) = q1*sin(phi_n)+q2*cos(phi_n)
            !
          elseif (ioper<=1+2*N_Cn+Nrot) then !  C'2
            !
            irot = ioper-(1+2*N_Cn)
            !
            phi_n = phi*irot
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = q1*cos(phi_n)+q2*sin(phi_n)
            dst(4) = q1*sin(phi_n)-q2*cos(phi_n)
            !
          elseif (ioper==1+2*N_Cn+Nrot+1) then ! sigmah
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = src(3)
            dst(4) = src(4)
            !
          elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn) then !  Sinf
            !
            ioper_ = ioper-(1+2*N_Cn+Nrot+1)
            !
            irot = (ioper_+1)/2
            !
            phi_n = phi*irot
            !
            ! Second oper in a class is with negative phi
            if (mod(ioper_,2)==0)  phi_n = -phi_n
            !
            dst(1) = src(2)
            dst(2) = src(1)
            dst(3) = (q1*cos(phi_n)-q2*sin(phi_n))
            dst(4) = (q1*sin(phi_n)+q2*cos(phi_n))
            !
          elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then !  C'2
            !
            irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
            !
            phi_n = phi*irot
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = ( q1*cos(phi_n)+q2*sin(phi_n))
            dst(4) = ( q1*sin(phi_n)-q2*cos(phi_n))
            !
          else
            !
            write (out,"('ML_symmetry_transformation_XY2  in Dinfty: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_XY2 Dinfty - bad operation. type'
 
          endif 
          !
          !do ioper_=1,sym%Noper/2
          !  repres ( ioper_,            :,:)= sym%irr(5, ioper_)%repres
          !  repres ( ioper_+sym%Noper/2,:,:)= sym%irr(5, ioper_)%repres
          !enddo
          !
          do ioper_=1,sym%Noper
            repres ( ioper_,            :,:)= sym%irr(5, ioper_)%repres
            !repres ( ioper_+sym%Noper/2,:,:)= sym%irr(5, ioper_)%repres
          enddo
          !
          !dst(3) = repres(ioper,1,1)*src(3)+repres(ioper,1,2)*src(4)
          !dst(4) = repres(ioper,2,1)*src(3)+repres(ioper,2,2)*src(4)
          !
       end select       
       !
    end select
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY2/end')") 
    !
  end subroutine ML_symmetry_transformation_XY2


  ! Here we find the rotational symmetry from J,K
  !
  subroutine ML_rotsymmetry_XY2(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
    integer(ik),intent(out) :: gamma,ideg
    integer(ik)             :: N,k_,l,N_Cn
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_XY2/start')") 
    !
    ! the trivial case 
    select case(trim(molec%symmetry))
    case('C','C(M)')
         gamma = 1
         ideg = 1
         return
    end select 
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      select case(trim(molec%symmetry))
      case default
         !
         write (out,"('ML_rotsymmetry_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
         stop 'ML_rotsymmetry_XY2 - bad symm. type'
         !
      case('C2V','C2V(M)')
         !
         gamma = 0 
         ideg = 1
         if (mod(K+2,2)==0.and.tau==0) gamma = 1 !; return
         if (mod(K+2,2)==0.and.tau==1) gamma = 2 !; return
         if (mod(K+2,2)/=0.and.tau==0) gamma = 4 !; return
         if (mod(K+2,2)/=0.and.tau==1) gamma = 3 !; return
         !
      case('C2H','C2H(M)')
         !
         gamma = 0 
         ideg = 1
         if (mod(K+2,2)==0.and.tau==0) gamma = 1 !1 !1 !1 !1 !1 ! ; return
         if (mod(K+2,2)==0.and.tau==1) gamma = 3 !2 !4 !3 !4 !2 ! ; return
         if (mod(K+2,2)/=0.and.tau==0) gamma = 2 !3 !3 !4 !2 !4 ! ; return
         if (mod(K+2,2)/=0.and.tau==1) gamma = 4 !4 !2 !2 !3 !3 ! ; return
         !
      case('CS','CS(M)')
         !
         gamma = 0 
         ideg = 1
         !
         if (molec%AtomMasses(2)/=molec%AtomMasses(3).or.molec%req(1)/=molec%req(2)) then
           !
           if (mod(tau+2,2)==0) gamma = 1 !; return
           if (mod(tau+2,2)/=0) gamma = 2 !; return
           !
         else
           !
           if (mod(K+tau+2,2)==0) gamma = 1 !; return
           if (mod(K+tau+2,2)/=0) gamma = 2 !; return
           !
         endif 
         !
      case('D2H(M)')
         !
         gamma = 0 
         ideg = 1
         if (mod(K+2,2)==0.and.tau==0) gamma = 1 !1 !1 !; return
         if (mod(K+2,2)==0.and.tau==1) gamma = 7 !3 !3 !; return
         if (mod(K+2,2)/=0.and.tau==0) gamma = 3 !7 !7 !; return
         if (mod(K+2,2)/=0.and.tau==1) gamma = 5 !5 !5 !; return
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
            write(out,"('ML_rotsymmetry_XY2-D3h: illegal j,k,tau - ',3i8)") j,k,tau
            !
         endif 
         !
      case('DNH','DNH(M)')
         !
         gamma = 0 
         ideg = 1
         !
         N = sym%N
         N_Cn = sym%N/2
         k_ = mod(K+N_Cn,N_Cn)
         l = k_ ; if (k_>N_Cn) l = sym%N-k_
         !
         if (mod(sym%N,2)==1) then
            !
            if (mod(K+N_Cn,N_Cn)==0) then
               !
               if     (tau==0.and.mod(k+2,2)==0) then 
                  gamma = 1 
               elseif (tau==1.and.mod(k+2,2)==0) then 
                  gamma = 2
               elseif (tau==0.and.mod(k+2,2)/=0) then 
                  gamma = 4
               elseif (tau==1.and.mod(k+2,2)/=0) then 
                  gamma = 3
               else
                  stop 'ML_rotsymmetry_abcd-Dnh: illegal k,tau (K mod N  = 0)'
               endif
               !
            elseif (tau<=1.and.k<=j) then
               !
               ideg = 1 ! tau +1
               if (mod(k+tau,2)/=0) ideg = 2
               !
               if     (mod(k+2,2)==0) then 
                   gamma = 4+2*l-1
               else
                  gamma = 4+2*l
               endif
               !
            else
                 stop 'ML_rotsymmetry_abcd-Dnh: illegal k,tau (K mod N  /= 0)'
            endif
            !
         else ! even Dnh
            !
            if (mod(K+N_Cn,N_Cn)==0) then
               !
               if     (tau==0.and.mod(k+2,2)==0) then 
                  gamma = 1 
               elseif (tau==1.and.mod(k+2,2)==0) then 
                  gamma = 2
               elseif (tau==0.and.mod(k+2,2)/=0.and.mod(N_Cn,2)/=0) then 
                  gamma = 4
               elseif (tau==1.and.mod(k+2,2)/=0.and.mod(N_Cn,2)/=0) then 
                  gamma = 3
               else
                  stop 'ML_rotsymmetry_abcd-Dnh: illegal k,tau (K mod N  = 0)'
               endif
               !
            elseif (tau<=1.and.k<=j) then
               !
               !ideg = tau +1
               !
               ideg = 1
               !
               if (mod(k+tau,2)/=0) ideg = 2
               !
               gamma = 8+2*l-1
               !
               !if     (mod(k+2,2)==0) then 
               !    gamma = 8+2*l
               !else
               !    gamma = 8+2*l-1
               !endif
               !
            else
                 stop 'ML_rotsymmetry_abcd-Dnh: illegal k,tau (K mod N  /= 0)'
            endif
            !
         endif
         !
      end select
      !
    case('R-RHO-Z')
      !
      select case(trim(molec%symmetry))
      case default
         !
         write (out,"('ML_rotsymmetry_XY2: symmetry ',a,' unknown')") trim(molec%symmetry)
         stop 'ML_rotsymmetry_XY2 - bad symm. type'
         !
      case('CS','CS(M)')
         !
         gamma = 0 
         ideg = 1
         !
         if (molec%AtomMasses(2)/=molec%AtomMasses(3).or.molec%req(1)/=molec%req(2)) then
           !
           if (mod(tau+2,2)==0) gamma = 1 !; return
           if (mod(tau+2,2)/=0) gamma = 2 !; return
           !
         else
           !
           if (mod(K+tau+2,2)==0) gamma = 1 !; return
           if (mod(K+tau+2,2)/=0) gamma = 2 !; return
           !
         endif 
         !
      case('C2V','C2V(M)')
         !
         gamma = 0 
         ideg = 1
         if (mod(K+2,2)==0.and.tau==0) gamma = 1 !; return
         if (mod(K+2,2)==0.and.tau==1) gamma = 3 !; return
         if (mod(K+2,2)/=0.and.tau==0) gamma = 4 !; return
         if (mod(K+2,2)/=0.and.tau==1) gamma = 2 !; return
         !
      case('C2VN')
         !
         gamma = 0 
         ideg = 1
         if (mod(K+2,2)==0.and.tau==0) gamma = 1+4*K
         if (mod(K+2,2)==0.and.tau==1) gamma = 4+4*K
         if (mod(K+2,2)/=0.and.tau==0) gamma = 2+4*K
         if (mod(K+2,2)/=0.and.tau==1) gamma = 3+4*K
         !
      case('DNH','DNH(M)')
         !
         gamma = 0 
         ideg = 1
         !
         N = sym%N
         N_Cn = sym%N/2
         k_ = mod(K+N_Cn,N_Cn)
         l = k_ ; if (k_>N_Cn) l = sym%N-k_
         !
         if (mod(sym%N,2)==1) then
            !
            if (mod(K+N_Cn,N_Cn)==0) then
               !
               if     (tau==0.and.mod(k+2,2)==0) then 
                  gamma = 1 
               elseif (tau==1.and.mod(k+2,2)==0) then 
                  gamma = 2
               elseif (tau==0.and.mod(k+2,2)/=0) then 
                  gamma = 4
               elseif (tau==1.and.mod(k+2,2)/=0) then 
                  gamma = 3
               else
                  stop 'ML_rotsymmetry_abcd-Dnh: illegal k,tau (K mod N  = 0)'
               endif
               !
            elseif (tau<=1.and.k<=j) then
               !
               ideg = 1 ! tau +1
               if (mod(k+tau,2)/=0) ideg = 2
               !
               if     (mod(k+2,2)==0) then 
                   gamma = 4+2*l-1
               else
                  gamma = 4+2*l
               endif
               !
            else
                 stop 'ML_rotsymmetry_abcd-Dnh: illegal k,tau (K mod N  /= 0)'
            endif
            !
         else ! even Dnh
            !
            if (mod(K+N_Cn,N_Cn)==0) then
               !
               if     (tau==0.and.mod(k+2,2)==0) then 
                  gamma = 1 
               elseif (tau==1.and.mod(k+2,2)==0) then 
                  gamma = 2
               elseif (tau==0.and.mod(k+2,2)/=0.and.mod(N_Cn,2)/=0) then 
                  gamma = 4
               elseif (tau==1.and.mod(k+2,2)/=0.and.mod(N_Cn,2)/=0) then 
                  gamma = 3
               else
                  stop 'ML_rotsymmetry_abcd-Dnh: illegal k,tau (K mod N  = 0)'
               endif
               !
            elseif (tau<=1.and.k<=j) then
               !
               !ideg = tau +1
               !
               ideg = 1
               !
               if (mod(k+tau,2)/=0) ideg = 2
               !
               gamma = 8+2*l-1
               !
               !if     (mod(k+2,2)==0) then 
               !    gamma = 8+2*l
               !else
               !    gamma = 8+2*l-1
               !endif
               !
            else
                 stop 'ML_rotsymmetry_abcd-Dnh: illegal k,tau (K mod N  /= 0)'
            endif
            !
         endif
         !
      end select
      !
    end select 
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_XY2/end')") 
    !
    !
  end subroutine ML_rotsymmetry_XY2


  !
  ! Defining MEP function for CH2 molecule
  !
  function ML_MEP_xy2_R12_R(x)  result(f)

   real(ark),intent(in) ::  x
   real(ark)            ::  f,r12_inf,r12_0,a

     r12_inf = molec%mep_params(1)
     r12_0   = molec%mep_params(2)
     a       = molec%mep_params(3)

     f = r12_inf+(r12_0-r12_inf)*exp(-a*sqrt(x)**9)
 
  end function ML_MEP_xy2_R12_R


  !
  ! Defining MEP function for CH2 molecule
  !
  function ML_MEP_xy2_R12_alpha(x)  result(f)

   real(ark),intent(in) ::  x
   real(ark)            ::  f,a0,a1,a2,a3,a4,a5,a6,y

        a0      = molec%mep_params(1)
        a1      = molec%mep_params(2)
        a2      = molec%mep_params(3)
        a3      = molec%mep_params(4)
        a4      = molec%mep_params(5)
        a5      = molec%mep_params(6)
        a6      = molec%mep_params(7)
        !
        y = 1.0d0/(x+small_)
        !
        f = a0+a1*y+a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6
 
  end function ML_MEP_xy2_R12_alpha


  !
  ! Defining MEP function for CH2 molecule
  !
  function ML_MEP_xy2_R13_alpha(x)  result(f)

   real(ark),intent(in) ::  x
   real(ark)            ::  f,a0,a1,a2,a3,a4,a5,a6,y,x0

        a0      = molec%mep_params(1)
        a1      = molec%mep_params(2)
        a2      = molec%mep_params(3)
        a3      = molec%mep_params(4)
        a4      = molec%mep_params(5)
        a5      = molec%mep_params(6)
        a6      = molec%mep_params(7)
        !
        x0  =  molec%alphaeq(1)
        !
        y = cos(x)-cos(x0)
        !
        f = a0+a1*y+a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6
 
  end function ML_MEP_xy2_R13_alpha


  
  end module mol_xy2

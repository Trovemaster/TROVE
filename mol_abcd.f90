!
!  This unit defines all specific routines for a fouratomic molecule of ABCD type
!
module mol_abcd
  use accuracy
  use moltype
  use lapack
  use pot_abcd
  use symmetry,only : sym

  implicit none

  public ML_b0_ABCD,ML_coordinate_transform_abcd,ML_rotsymmetry_abcd,ML_symmetry_transformation_abcd
  !
  private
  !
  integer(ik), parameter :: verbose = 4  ! Verbosity level
  !
  contains
  !
  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter 
  ! as conjugate momenta coordinates
  !
  function ML_coordinate_transform_abcd(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: dsrc(size(src))
    real(ark)                 :: rho,r_eq(6)
    !
    integer(ik) :: nsrc
    !
    if (verbose>=6) write(out,"('ML_coordinate_transform_abcd/start')") 
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
    case('NORMAL')
       !
       if (direct) then 
           dst = src
       else
           dst = src
           !dst(:) = dst(:) +  molec%local_eq(:)
       endif
       !
    case('R-ALPHA-TAU')
       !
       if (direct) then
          ! 
          dst(1:5) = dsrc(1:5)
          dst(6) = src(6)
          !
          !dst(6) = mod(dst(6)+2.0_ark*pi,2.0_ark*pi)
          !
      else ! not direct
          !
          dst(1:5) = dsrc(1:5)+molec%local_eq(1:5)
          dst(6) = src(6)
          !
      endif
       !
    case('R-ALPHA-R-TAU')
       !
       if (direct) then
          ! 
          dst(3) = dsrc(1)
          dst(1:2) = dsrc(2:3)
          dst(4:5) = dsrc(4:5)
          dst(6) = src(6)
          !
          !dst(6) = mod(dst(6)+2.0_ark*pi,2.0_ark*pi)
          !
      else ! not direct
          !
          dst(1) = dsrc(3)+molec%local_eq(1)
          dst(2:3) = dsrc(1:2)+molec%local_eq(2:3)
          dst(4:5) = dsrc(4:5)+molec%local_eq(4:5)
          dst(6) = src(6)
          !
      endif
       !
    case('R-ALPHA-TAU-REF')
       !
       rho = src(6)
       !
       r_eq(1:6) = ML_MEP_ABCD_tau_ref(rho)
       !
       !rref(1,0:4)  = molec%force(1:5) 
       !rref(2,0:4)  = molec%force(6:10) 
       !rref(3,0:4)  = molec%force(11:15) 
       !rref(4,0:4)  = molec%force(16:20)/rad 
       !rref(5,0:4)  = molec%force(21:25)/rad  

       !k(1)=1; k(2)=2; k(3)=3; k(4)=4
       !
       !r_eq(1) = rref(1,0)+sum(rref(1,1:4)*(1.0_rk+cos(rho))**k(1:4))
       !r_eq(2) = rref(2,0)+sum(rref(2,1:4)*(1.0_rk+cos(rho))**k(1:4))
       !r_eq(3) = rref(3,0)+sum(rref(3,1:4)*(1.0_rk+cos(rho))**k(1:4))
       !r_eq(4) = rref(4,0)+sum(rref(4,1:4)*(1.0_rk+cos(rho))**k(1:4))
       !r_eq(5) = rref(5,0)+sum(rref(5,1:4)*(1.0_rk+cos(rho))**k(1:4))
       !r_eq(6) = rho
       !
       if (direct) then
          !
          dst(1:5) = src(1:5) - r_eq(1:5)
          !
          dst(6) = src(6)
          !
          !dst(6) = mod(dst(6)+2.0_ark*pi,2.0_ark*pi)
          !
      else ! not direct
          !
          dst(1:5) = src(1:5) + r_eq(1:5)
          dst(6) = src(6)
          !
      endif
       !
    case('R-S-TAU-REF')
       !
       rho = src(6)
       !
       r_eq(1:6) = ML_MEP_ABCD_tau_ref(rho)
       !
       if (direct) then
          !
          !dst(1:3) = src(1:3)-r_eq(1:3)
          !
          dst(1) = src(1)-r_eq(1)
          !
          dst(2) = 0.5_ark*(src(2)+src(3)-r_eq(2)-r_eq(3))
          dst(3) = 0.5_ark*(src(2)-src(3))
          !
          dst(4) = 0.5_ark*(src(4)+src(5)-r_eq(4)-r_eq(5))
          dst(5) = 0.5_ark*(src(4)-src(5))
          !
          dst(6) = src(6)
          !
      else ! not direct
          !
          dst(1) = src(1)+r_eq(1)
          !
          dst(2) = src(2)+src(3)+0.5_ark*(r_eq(2)+r_eq(3))
          dst(3) = src(2)-src(3)+0.5_ark*(r_eq(2)+r_eq(3))
          !
          dst(4) = src(4)+src(5)+0.5_ark*(r_eq(4)+r_eq(5))
          dst(5) = src(4)-src(5)+0.5_ark*(r_eq(4)+r_eq(5))
          !
          dst(6) = src(6)
          !
      endif
       !
    case('R-S-TAU')
       !
       rho = src(6)
       !
       r_eq(1:5) = molec%local_eq(1:5)
       !
       if (direct) then
          !

          dst(1:3) = src(1:3)-r_eq(1:3)
          !
          dst(4) = 0.5_ark*(src(4)+src(5)-r_eq(4)-r_eq(5))
          dst(5) = 0.5_ark*(src(4)-src(5))
          !
          dst(6) = src(6)
          !
      else ! not direct
          !
          dst(1:3) = src(1:3)+r_eq(1:3)
          !
          dst(4) = src(4)+src(5)+0.5_ark*(r_eq(4)+r_eq(5))
          dst(5) = src(4)-src(5)+0.5_ark*(r_eq(4)+r_eq(5))
          !
          dst(6) = src(6)
          !
      endif
       !
    case('R-PHI-TAU')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          dst(4) =  dsrc(4)
          dst(5) = -dsrc(5)
          !
          dst(6) = src(6)
          !
          !dst(6) = mod(dst(6)+2.0_ark*pi,2.0_ark*pi)
          !
      else ! not direct
          !
          dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
          !
          dst(4) = dsrc(4)+molec%local_eq(4)
          dst(5) =-dsrc(5)+molec%local_eq(4)
          !
          dst(6) = src(6)
          !
      endif
       !
    case('R-PHI-TAU-REF')
       !
       rho = src(6)
       !
       r_eq(1:6) = ML_MEP_ABCD_tau_ref(rho)
       !
       if (direct) then
          !
          dst(1:3) = src(1:3) - r_eq(1:3)
          !
          dst(4) =  (src(4) - r_eq(4))
          dst(5) = -(src(5) - r_eq(5))
          !
          dst(6) = src(6)
          !
          !dst(6) = mod(dst(6)+2.0_ark*pi,2.0_ark*pi)
          !
      else ! not direct
          !
          dst(1:3) = src(1:3) + r_eq(1:3)
          !
          dst(4) = src(4) + r_eq(4)
          dst(5) =-src(5) + r_eq(5)
          !
          dst(6) = src(6)
          !
      endif
        !
    case('R-ALPHA-TAU-1D')
       !
       if (direct) then
          ! 
          dst(1) = dsrc(6)
          !
          !dst(6) = mod(dst(6)+2.0_ark*pi,2.0_ark*pi)
          !
      else ! not direct
          !
          dst(1:5) = molec%local_eq(1:5)
          !
          dst(6) = dsrc(1)+molec%local_eq(6)
          !
      endif
        !
    case('R-R1-R2-T1-T2-T3-T4')
       !
       if (direct) then 
          !
          dst(1:) = dsrc(1:)
          !
       else
          dst(1:) = src(1:)+molec%local_eq(1:)
          !
       endif
       !
    case('R-R1-R2-TX-TY-TX-TY')
       !
       if (direct) then 
          !
          dst(1:3) = dsrc(1:3)
          !
          dst(4) =-dsrc(5)
          dst(5) = dsrc(4)
          dst(6) =-dsrc(7)
          dst(7) = dsrc(6)
          !
       else
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          dst(4) = src(5)+molec%local_eq(5)
          dst(5) =-src(4)-molec%local_eq(4)
          dst(6) = src(7)+molec%local_eq(7)
          dst(7) =-src(6)-molec%local_eq(6)
          !
       endif
       !
    case('R-R1-R2+TX-TY+TX-TY')
       !
       if (direct) then 
          !
          dst(1:3) = dsrc(1:3)
          !
          dst(4) = dsrc(5)
          dst(5) =-dsrc(4)
          dst(6) = dsrc(7)
          dst(7) =-dsrc(6)
          !
       else
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          dst(4) =-src(5)-molec%local_eq(5)
          dst(5) = src(4)+molec%local_eq(4)
          dst(6) =-src(7)-molec%local_eq(7)
          dst(7) = src(6)+molec%local_eq(6)
          !
       endif
       !
    case('R-R1-R2-TX+TY+TX-TY')
       !
       if (direct) then 
          !
          dst(1:3) = dsrc(1:3)
          !
          dst(4) =-dsrc(5)
          dst(5) = dsrc(4)
          dst(6) = dsrc(7)
          dst(7) =-dsrc(6)
          !
       else
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          dst(4) = src(5)+molec%local_eq(5)
          dst(5) =-src(4)-molec%local_eq(4)
          dst(6) =-src(7)-molec%local_eq(7)
          dst(7) = src(6)+molec%local_eq(6)
          !
       endif
       !
    case('R-R1-R2-TX+TY-TX+TY')
       !
       if (direct) then 
          !
          dst(1:3) = dsrc(1:3)
          !
          dst(4) =-dsrc(5)
          dst(5) = dsrc(4)
          dst(6) =-dsrc(7)
          dst(7) = dsrc(6)
          !
       else
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          dst(4) = src(5)+molec%local_eq(5)
          dst(5) =-src(4)-molec%local_eq(4)
          dst(6) = src(7)+molec%local_eq(7)
          dst(7) =-src(6)-molec%local_eq(6)
          !
       endif       
       !
    case('R-R1-R2+Y-X-Y+X')
       !
       if (direct) then 
          !
          dst(1:3) = dsrc(1:3)
          !
          dst(4) = dsrc(5)
          dst(5) =-dsrc(4)
          dst(6) =-dsrc(7)
          dst(7) = dsrc(6)
          !
       else
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          dst(4) =-src(5)-molec%local_eq(5)
          dst(5) = src(4)+molec%local_eq(4)
          dst(6) = src(7)+molec%local_eq(7)
          dst(7) =-src(6)-molec%local_eq(6)
          !
       endif
       !
    case('R-R1-R2-X-Y+X+Y')
       !
       if (direct) then 
          !
          dst(1:3) = dsrc(1:3)
          !
          dst(4) =-dsrc(4)
          dst(5) =-dsrc(5)
          dst(6) = dsrc(6)
          dst(7) = dsrc(7)
          !
       else
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          dst(4) =-src(4)-molec%local_eq(4)
          dst(5) =-src(5)-molec%local_eq(5)
          dst(6) = src(6)+molec%local_eq(6)
          dst(7) = src(7)+molec%local_eq(7)
          !
       endif
       !
    case('R-R1-R2+X+Y-X-Y')
       !
       if (direct) then 
          !
          dst(1:3) = dsrc(1:3)
          !
          dst(4) = dsrc(4)
          dst(5) = dsrc(5)
          dst(6) =-dsrc(6)
          dst(7) =-dsrc(7)
          !
       else
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          dst(4) = src(4)+molec%local_eq(4)
          dst(5) = src(5)+molec%local_eq(5)
          dst(6) =-src(6)-molec%local_eq(6)
          dst(7) =-src(7)-molec%local_eq(7)
          !
       endif       
       !
    case('R1-R2-Y+X+Y-X-R')
       !
       if (direct) then 
          !
          dst(7) = src(1)-molec%local_eq(1)
          dst(1:2) = src(2:3)-molec%local_eq(2:3)
          !
          dst(3) =-src(5)
          dst(4) = src(4)
          dst(5) = src(7)
          dst(6) =-src(6)
          !
       else
          !
          dst(1) = src(7)+molec%local_eq(7)
          dst(2:3) = src(1:2)+molec%local_eq(1:2)
          !
          dst(4) = src(4)+molec%local_eq(4)
          dst(5) =-src(3)-molec%local_eq(3)
          dst(6) =-src(6)-molec%local_eq(6)
          dst(7) = src(5)+molec%local_eq(5)
          !
       endif       
          !
    case('R-R1-R2-Y+X+Y-X')
       !
       if (direct) then 
          !
          dst(1:3) = dsrc(1:3)
          !
          dst(4) =-dsrc(5)
          dst(5) = dsrc(4)
          dst(6) = dsrc(7)
          dst(7) =-dsrc(6)
          !
       else
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          dst(4) = src(5)+molec%local_eq(5)
          dst(5) =-src(4)-molec%local_eq(4)
          dst(6) =-src(7)-molec%local_eq(7)
          dst(7) = src(6)+molec%local_eq(6)
          !
       endif       
       !
    case('R-R1-R2-X-Y-X-Y')
       !
       if (direct) then 
          !
          dst(1) = dsrc(1)
          dst(2) = dsrc(2)
          dst(3) = dsrc(3)
          !
          dst(4) = dsrc(4)
          dst(5) = dsrc(5)
          dst(6) = dsrc(6)
          dst(7) = dsrc(7)
          !
       else
          !
          dst(1) = src(1)+molec%local_eq(1)
          dst(2) = src(2)+molec%local_eq(2)
          dst(3) = src(3)+molec%local_eq(3)
          dst(4) = src(4)+molec%local_eq(4)
          dst(5) = src(5)+molec%local_eq(5)
          dst(6) = src(6)+molec%local_eq(6)
          dst(7) = src(7)+molec%local_eq(7)
          !
       endif
       !
    case('R-R1-R2+X+Y-X+Y')
       !
       if (direct) then 
          !
          dst(1) = dsrc(1)
          dst(2) = dsrc(2)
          dst(3) = dsrc(3)
          !
          dst(4) = dsrc(4)
          dst(5) = dsrc(5)
          dst(6) =-dsrc(6)
          dst(7) = dsrc(7)
          !
       else
          !
          dst(1) = src(1)+molec%local_eq(1)
          dst(2) = src(2)+molec%local_eq(2)
          dst(3) = src(3)+molec%local_eq(3)
          dst(4) = src(4)+molec%local_eq(4)
          dst(5) = src(5)+molec%local_eq(5)
          dst(6) =-src(6)+molec%local_eq(6)
          dst(7) = src(7)+molec%local_eq(7)
          !
       endif
       !
    case('R-Z1-Z2-X-Y-X-Y')
       !
       if (direct) then 
          !
          dst(1) = dsrc(1)
          dst(2) =-dsrc(2)
          dst(3) = dsrc(3)
          !
          dst(4) = dsrc(4)
          dst(5) = dsrc(5)
          dst(6) = dsrc(6)
          dst(7) = dsrc(7)
          !
       else
          !
          dst(1) = src(1)+molec%local_eq(1)
          dst(2) =-src(2)+molec%local_eq(2)
          dst(3) = src(3)+molec%local_eq(3)
          dst(4) = src(4)+molec%local_eq(4)
          dst(5) = src(5)+molec%local_eq(5)
          dst(6) = src(6)+molec%local_eq(6)
          dst(7) = src(7)+molec%local_eq(7)
          !
       endif
       !
    case('R-Z1+Z2-Y-X-Y-X')
       !
       if (direct) then 
          !
          dst(1) = dsrc(1)
          dst(2) =-dsrc(2)
          dst(3) = dsrc(3)
          !
          dst(4) = dsrc(5)
          dst(5) = dsrc(4)
          dst(6) = dsrc(7)
          dst(7) = dsrc(6)
          !
       else
          !
          dst(1) = src(1)+molec%local_eq(1)
          dst(2) =-src(2)+molec%local_eq(2)
          dst(3) = src(3)+molec%local_eq(3)
          dst(4) = src(5)+molec%local_eq(5)
          dst(5) = src(4)+molec%local_eq(4)
          dst(6) = src(7)+molec%local_eq(7)
          dst(7) = src(6)+molec%local_eq(6)
          !
       endif       
       !
    case('R-Z1+Z2+X+Y-X-Y')
       !
       if (direct) then 
          !
          dst(1) = dsrc(1)
          dst(2) =-dsrc(2)
          dst(3) = dsrc(3)
          !
          dst(4) = dsrc(4)
          dst(5) = dsrc(5)
          dst(6) =-dsrc(6)
          dst(7) =-dsrc(7)
          !
       else
          !
          dst(1) = src(1)+molec%local_eq(1)
          dst(2) =-src(2)+molec%local_eq(2)
          dst(3) = src(3)+molec%local_eq(3)
          dst(4) = src(4)+molec%local_eq(4)
          dst(5) = src(5)+molec%local_eq(5)
          dst(6) =-src(6)+molec%local_eq(6)
          dst(7) =-src(7)+molec%local_eq(7)
          !
       endif
       !
    case('R-Z1+Z2-X-Y+X+Y')
       !
       if (direct) then 
          !
          dst(1) = dsrc(1)
          dst(2) =-dsrc(2)
          dst(3) = dsrc(3)
          !
          dst(4) =-dsrc(4)
          dst(5) =-dsrc(5)
          dst(6) = dsrc(6)
          dst(7) = dsrc(7)
          !
       else
          !
          dst(1) = src(1)+molec%local_eq(1)
          dst(2) =-src(2)+molec%local_eq(2)
          dst(3) = src(3)+molec%local_eq(3)
          !
          dst(4) =-src(4)-molec%local_eq(4)
          dst(5) =-src(5)-molec%local_eq(5)
          dst(6) = src(6)+molec%local_eq(6)
          dst(7) = src(7)+molec%local_eq(7)
          !
       endif       
       !
    end select
    !
    !
    if (verbose>=6) write(out,"('ML_coordinate_transform_abcd/end')") 
    !
    !
  end function ML_coordinate_transform_abcd


  ! Here we define structural parameters a0 for ABCD molecule,
  !
  subroutine ML_b0_ABCD(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
     !
     !
     integer(ik),intent(in) :: Npoints,Natoms
     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in) ,optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional :: rho_ref
     real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
     !
     real(ark)    ::  rho,transform(3,3),c(3,3),a0(molec%Natoms,3),CM_shift,re1,re2,re3,ae1,ae2,r_eq(6),phi,ayz,ayy,azz
     integer(ik)  ::  n,i,ix,jx,i0,in,i1,ipar,istep,j0,Nangles,Nbonds
     !
     real(ark)    :: Inert0(3),Inert(3),Inert1(3),a(3,3),b(3,1),x(5),a0_tt,a0_t(molec%Natoms,3)
     real(ark)                        :: rho_ark
      !
      if (verbose>=4) write(out,"('ML_b0_ABCD/start')") 
      !
      if (size(molec%req)/=3) then
        write(out,"('ML_b0_ABCD: Nbonds must be 3 in this routine, not  ',i9)") size(molec%req)
        stop 'ML_b0_ABCD: wrong Nbonds '
      endif 

      if (molec%Natoms/=4) then
        write(out,"('ML_b0_ABCD: Natoms must be 4 in this routine, not  ',i9)") molec%Natoms
        stop 'ML_b0_ABCD: wrong Natoms '
      endif 
      !
      Nbonds = molec%Nbonds
      Nangles = molec%Nangles
      !
      re1   = molec%req(1)
      re2   = abs(molec%req(2))
      re3   = abs(molec%req(3))
      !
      rho = 0 
      !
      if (Nangles>1) then
        !
        ae1   = molec%alphaeq(1)
        ae2   = molec%alphaeq(2)
        !
        rho = molec%taueq(1)
        !
      elseif (Nangles>0.and.molec%Ndihedrals>1) then
        !
        ae1   = molec%alphaeq(1)
        ae2   = pi+(-molec%taueq(1)+molec%taueq(2)) 
        !
      elseif (Nangles==0.and.molec%Ndihedrals>3) then 
        !
        ae1   = pi+(-molec%taueq(1)+molec%taueq(2)) 
        ae2   = pi+(-molec%taueq(3)+molec%taueq(4))
        !
      else
        write(out,"('ML_b0_ABCD: case undefined')")
        stop 'ML_b0_ABCD not defined'
      endif 
      !
      if (trim(molec%coords_transform)=='R-ALPHA-TAU-REF'.or.trim(molec%coords_transform)=='R-S-TAU-REF'.or.&
          trim(molec%coords_transform)=='R-PHI-TAU-REF') then 
         !
         rho_ark = rho
         !
         r_eq(1:6) = ML_MEP_ABCD_tau_ref(rho_ark)
         !
         re1   = r_eq(1)
         re2   = r_eq(2)
         re3   = r_eq(3)
         ae1   = r_eq(4)
         ae2   = r_eq(5)
         !
         if (verbose>=5) then 
           write(out,"(i8,f18.8,' r1 r2 r3 a1 a2 ',5f16.8)") i,rho*180./pi,re1,re2,re3,ae1,ae2
         endif
         !
      endif 
      !
      if (molec%AtomMasses(1)==molec%AtomMasses(2).and.molec%AtomMasses(3)==molec%AtomMasses(4).and..false.) then
        !
        phi = rho*0.5_ark

        a0(1,1) = 0.0_ark
        a0(1,2) = 0.0_ark
        a0(1,3) = re1*0.5_ark
        !
        a0(2,1) = 0.0_ark
        a0(2,2) = 0.0_ark
        a0(2,3) = -re1*0.5_ark
        !
        a0(3,1) = re2*sin(ae1)*cos(phi)
        a0(3,2) =-re2*sin(ae1)*sin(phi)
        a0(3,3) =-re2*cos(ae1)+re1*0.5_ark
        !
        a0(4,1) = re2*sin(ae1)*cos(phi)
        a0(4,2) = re2*sin(ae1)*sin(phi)
        a0(4,3) = re2*cos(ae1)-re1*0.5_ark
        !
        do n = 1,3 
          CM_shift = sum(a0(:,n)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
          a0(:,n) = a0(:,n) - CM_shift
        enddo
        !
        ayz =-sum(molec%AtomMasses(:)*a0(:,2)*a0(:,3) )
        ayy = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,3)**2) )
        azz = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,2)**2) )
        !
        phi = atan2(2.0_ark*ayz,ayy-azz)
        !
        !phi = atan(ayz/(ayy-azz))
        !
        transform = 0 
        !
        transform(1,1) = 1.0_ark
        transform(2,2) = cos(phi)
        transform(3,3) = cos(phi)
        transform(2,3) = sin(phi)
        transform(3,2) =-sin(phi)
        !
        do n = 1,molec%Natoms
           b0(n,:,0) = matmul(transpose(transform),a0(n,:))
        enddo
        !
        ayz =-sum(molec%AtomMasses(:)*b0(:,2,0)*b0(:,3,0) )
        ayy = sum(molec%AtomMasses(:)*( b0(:,1,0)**2+ b0(:,3,0)**2) )
        azz = sum(molec%AtomMasses(:)*( b0(:,1,0)**2+ b0(:,2,0)**2) )
        !
        call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,0),transform)
        !
      else 
        !
        a0(3,2) = re2*sin(ae1)
        a0(3,1) = 0.0_ark
        a0(3,3) = re2*cos(ae1)
        !
        a0(1,2) = 0.0_ark
        a0(1,1) = 0.0_ark
        a0(1,3) = 0.0_ark
        !
        a0(2,2) = 0.0_ark
        a0(2,1) = 0.0_ark
        a0(2,3) = re1
        !
        a0(4,2) = re3*sin(ae2)*cos(rho)
        a0(4,1) = re3*sin(ae2)*sin(rho)
        a0(4,3) = re1-re3*cos(ae2)
        !
        call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,transform)
        !
        a0(3,2) = re2*sin(ae1)*sin(rho*0.5_ark)
        a0(3,1) = re2*sin(ae2)*cos(rho*0.5_ark)
        a0(3,3) = re2*cos(ae1)
        !
        a0(1,2) = 0.0_ark
        a0(1,1) = 0.0_ark
        a0(1,3) = 0.0_ark
        !
        a0(2,2) = 0.0_ark
        a0(2,1) = 0.0_ark
        a0(2,3) = re1
        !
        a0(4,2) =-re3*sin(ae2)*sin(rho*0.5_ark)
        a0(4,1) = re3*sin(ae2)*cos(rho*0.5_ark)
        a0(4,3) = re1-re3*cos(ae2)        
        !
        do n = 1,3 
          CM_shift = sum(a0(:,n)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
          a0(:,n) = a0(:,n) - CM_shift
        enddo 
        !        
      endif 
      !
      b0(:,:,0) = a0(:,:)
      !
      ! We can also simulate the reference structure at each point of rho, if npoints presented
      !
      if (Npoints/=0) then
         ! 
         if (.not.present(rho_borders).or..not.present(rho_ref)) then  
            write(out,"('ML_b0_ABCD: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_ABCD: rho_borders or rho_ref not specified '
         endif
         !
         rho_ref = molec%taueq(1)
         !
         do i = 0,npoints
            !
            rho = rho_i(i)
            !
            re1   = molec%req(1)
            re2   = molec%req(2)
            re3   = molec%req(3)
            ae1   = molec%alphaeq(1)
            ae2   = molec%alphaeq(2)
            !
            if (trim(molec%coords_transform)=='R-ALPHA-TAU-REF'.or.trim(molec%coords_transform)=='R-S-TAU-REF'.or.&
                trim(molec%coords_transform)=='R-PHI-TAU-REF') then 
               !
               rho_ark = rho
               !
               r_eq(1:6) = ML_MEP_ABCD_tau_ref(rho_ark)
               !
               re1   = r_eq(1)
               re2   = r_eq(2)
               re3   = r_eq(3)
               ae1   = r_eq(4)
               ae2   = r_eq(5)
               !
               if (verbose>=5) then 
                 write(out,"(i8,f18.8,' r1 r2 r3 a1 a2 ',5f16.8)") i,rho*180./pi,re1,re2,re3,ae1,ae2
               endif
               !
            endif 
            !
            b0(3,2,i) = re2*sin(ae1)
            b0(3,1,i) = 0.0_ark
            b0(3,3,i) = re2*cos(ae1)
            !
            b0(1,2,i) = 0.0_ark
            b0(1,1,i) = 0.0_ark
            b0(1,3,i) = 0.0_ark
            !
            b0(2,2,i) = 0.0_ark
            b0(2,1,i) = 0.0_ark
            b0(2,3,i) = re1
            !
            b0(4,2,i) = re3*sin(ae2)*cos(rho)
            b0(4,1,i) = re3*sin(ae2)*sin(rho)
            b0(4,3,i) = re1-re3*cos(ae2)
            !
            ! Find center of mass
            !
            do n = 1,3 
              CM_shift = sum(b0(:,n,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
              b0(:,n,i) = b0(:,n,i) - CM_shift
            enddo 
            !
            do n = 1,molec%Natoms
               b0(n,:,i) = matmul(transpose(transform),b0(n,:,i))
            enddo
            !
            if (verbose>=8) then 
              write(out,"('b0',i8,12f12.8)") i,b0(:,:,i)
            endif
            !
         enddo
         !
         if (molec%AtomMasses(1)==molec%AtomMasses(2).and.molec%AtomMasses(3)==molec%AtomMasses(4).and..false.) then
            !
            do i = 0,npoints
               !
               rho = rho_i(i)
               !
               re1   = molec%req(1)
               re2   = molec%req(2)
               re3   = molec%req(3)
               ae1   = molec%alphaeq(1)
               ae2   = molec%alphaeq(2)
               !
               phi = rho*0.5_ark

               b0(1,1,i) = 0.0_ark
               b0(1,2,i) = 0.0_ark
               b0(1,3,i) = re1*0.5_ark

               b0(2,1,i) = 0.0_ark
               b0(2,2,i) = 0.0_ark
               b0(2,3,i) = -re1*0.5_ark

               b0(3,1,i) = re2*sin(ae1)*cos(phi)
               b0(3,2,i) =-re2*sin(ae1)*sin(phi)
               b0(3,3,i) =-re2*cos(ae1)+re1*0.5_ark

               b0(4,1,i) = re2*sin(ae1)*cos(phi)
               b0(4,2,i) = re2*sin(ae1)*sin(phi)
               b0(4,3,i) = re2*cos(ae1)-re1*0.5_ark
               !
               do n = 1,3 
                 CM_shift = sum(b0(:,n,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
                 b0(:,n,i) = b0(:,n,i) - CM_shift
               enddo 
               !
            enddo
            !
            !do i = 0,npoints
            !   !
            !   rho = rho_i(i)
            !   !
            !   phi = rho*0.5_ark+pi*0.5_ark
            !   !
            !   transform = 0 
            !   !
            !   transform(3,3) = 1.0_ark
            !   !
            !   transform(1,1) = cos(phi)
            !   transform(1,2) = sin(phi)
            !   transform(2,1) =-sin(phi)
            !   transform(2,2) = cos(phi)
            !   !
            !   forall(ix=1:4) b0(ix,:,i) = matmul(transpose(transform),b0(ix,:,i))
            !   !
            !   continue
            !   !
            !enddo
            !
         else 
           !
           i0 = npoints/2
           !
           do ipar = 0,1
             !
             istep = (-1)**(ipar+2)
             !
             i = npoints/2

             a0(:,:) = b0(:,:,i)
             !
             call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,transform)
             !
             Inert0(1) = sum(molec%AtomMasses(:)*( a0(:,2)**2+ a0(:,3)**2) )
             Inert0(2) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,3)**2) )
             Inert0(3) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,2)**2) )
             !
             do ix = 1,molec%Natoms
                a0(ix,:) = matmul(transpose(transform),b0(ix,:,i))
             enddo
             b0(:,:,i) = a0(:,:)
             !
             do j0 = 1,npoints/2
                !
                !write(out,'("j0 = ",i8)') j0
                !
                i = i + istep
                !
                do ix = 1,molec%Natoms
                   a0(ix,:) = matmul(transpose(transform),b0(ix,:,i))
                enddo
                !
                Inert(1) = sum(molec%AtomMasses(:)*( b0(:,2,i)**2+ b0(:,3,i)**2) )
                Inert(2) = sum(molec%AtomMasses(:)*( b0(:,1,i)**2+ b0(:,3,i)**2) )
                Inert(3) = sum(molec%AtomMasses(:)*( b0(:,1,i)**2+ b0(:,2,i)**2) )
                !
                ! Second Eckart equation
                ! 
                do ix = 1,3 
                   do jx = 1,3 
                      a(ix,jx) =  sum(molec%AtomMasses(:)*a0(:,ix)*a0(:,jx) )
                   enddo
                   !
                enddo
                !
                call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,a)
                !
                ! Found coordinate transformation "c"
                !
                transform = matmul(transform,real(a,kind=rk))
                !
                ! Transformation of a0 
                !
                do ix = 1,molec%Natoms
                   a0(ix,:) = matmul(transpose(transform),b0(ix,:,i))
                enddo
                !
                do ix = 1,3 
                   do jx = 1,3 
                      a(ix,jx) =  sum(molec%AtomMasses(:)*a0(:,ix)*a0(:,jx) )
                   enddo
                   !
                enddo
                !
                Inert(1) = sum(molec%AtomMasses(:)*( a0(:,2)**2+ a0(:,3)**2) )
                Inert(2) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,3)**2) )
                Inert(3) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,2)**2) )
                !
                Inert1 = 1
                !
                if (abs(i-i0)>2) then 
                  !
                  do in = 1,molec%Natoms
                     do ix = 1,3
                       !
                       N = min(4,abs(i-i0))
                       !
                       x(1:N+1) = rho_i(i-istep*N:i:istep)
                       !
                       call extrapolate(N,x,b0(in,ix,i-istep*N:i-istep:istep),a0_tt)
                       a0_t(in,ix) = a0_tt
                     enddo 
                  enddo
                  !
                  !
                  Inert1(1) = sum(molec%AtomMasses(:)*( a0_t(:,2)**2+ a0_t(:,3)**2) )
                  Inert1(2) = sum(molec%AtomMasses(:)*( a0_t(:,1)**2+ a0_t(:,3)**2) )
                  Inert1(3) = sum(molec%AtomMasses(:)*( a0_t(:,1)**2+ a0_t(:,2)**2) )
                  !
                endif 
                !
                if (all(abs(Inert(:)-Inert0(:))<5.0).or.all(abs(Inert(:)-Inert1(:))<5.0)) then 
                  Inert0 = Inert
                  b0(:,:,i) = a0(:,:)
                  cycle 
                endif 
                !
                if (verbose>=4) then 
                  write(out,"('b0',i4,12f12.8)") i,b0(:,:,i)
                endif
                !
             enddo
             !    
           enddo
           !
         endif
         !
      endif
      !
      !
      !b0 = cshift(b0(:,:,:),1,dim=2)
      !
      !b0(:,3,:) = -b0(:,3,:)

      !
      if (verbose>=3) then
        !
        do i = 0,npoints
           !
           write(out,"(i6)") molec%natoms
           !
           write(out,"(/'C',3x,3f14.8)") b0(1,:,i)
           write(out,"( 'C',3x,3f14.8)") b0(2,:,i)
           write(out,"( 'H',3x,3f14.8)") b0(3,:,i)
           write(out,"( 'H',3x,3f14.8)") b0(4,:,i)
           !
        enddo 
        !
        write(out,"('Coordinates:')")
        !
        do i = 0,npoints
          !
          write(out,"('b0',i4,12f12.8)") i,b0(:,:,i)
          !
        enddo 
      endif
      !
      if (verbose>=4) write(out,"('ML_b0_ABCD/end')") 
      !
    contains 
     !
     ! extrapolation at borders
     !
     subroutine extrapolate(N,x,src,dst)
     !
     integer(ik),intent(in)   :: N
     real(ark),intent(inout)  :: src(1:N),x(1:N+1)
     real(ark),intent(out)    :: dst

     integer(ik)        :: i1,i2
     real(rk)           :: a(N,N),b(N,1)
        !
        !
        do i1 = 1,N
           !
           a(i1,1) = 1.0_ark
           !
           b(i1,1) = src(i1)
           !
           do i2 = 2,N
             !
             a(i1,i2) = x(i1)**(i2-1)
             !
           enddo
        enddo
        !
        !  lapack_gelss 
        ! 
        call lapack_gelss(a(:,:),b(:,:))
        !
        dst = b(1,1)
        !
        do i1 = 2,N
          !
          dst = dst + real(b(i1,1),ark)*x(N+1)**(i1-1)
          !
        enddo
        !
   end subroutine extrapolate


  end subroutine ML_b0_ABCD


  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  subroutine ML_symmetry_transformation_abcd(ioper,natoms,src,dst)
    !
    integer(ik),intent(in)    :: ioper  ! group operation  
    integer(ik),intent(in)    :: natoms
    real(ark),intent(in)      :: src(natoms)
    real(ark),intent(out)     :: dst(natoms)
    logical                   :: extended = .false.
    !
    integer(ik) :: nsrc,Nrot,N_Cn,ioper_,irot,NC2
    real(ark)   :: q1x,q2x,q1y,q2y,phi_n,phi,repres(sym%Noper,2,2)
    !
    if (verbose>=7) write(out,"('ML_symmetry_transformation_abcd/start')") 
    !
    nsrc = size(src)
    !
    if(trim(molec%symmetry)=='C'.or.trim(molec%symmetry)=='C(M)') then 
       dst = src
       return 
    endif
    !
    if (molec%rho_border(2)>2.0_ark*pi) extended = .true.
    !
    select case(trim(molec%coords_transform))
        !
        case default
        write (out,"('ML_symmetry_transformation_abcd: coord_transf ',a,' unknown')") trim(molec%coords_transform)
        stop 'ML_symmetry_transformation_abcd - bad coord. type'
       !
    case('LINEAR','X-XE')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_abcd - bad symm. type'
          !
       case('CS','CS(M)')
         !
         select case(ioper)
         !
         case (1) ! identity 
           !
           dst = src
           !
         case (2) ! (E*)
           !
           dst(1:5) = src(1:5)
           dst(6) =-src(6)
           !
         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
         !
       case('D2H(M)')
           !
         select case(ioper)
           !
         case (1) ! E 
           !
           dst = src
           !
         case (2) ! (C2z) -> Rz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = src(6)
           !
         case (3) ! (C2x)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)
           !
         case (4) ! (C2y)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)
           !
         case (5) ! E* -> sigma_xz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) =-src(6)

         case (6) ! sigma_yz

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) =-dst(6)

         case (7) ! sigma_xy

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) =-dst(6)

         case (8) ! i

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) =-dst(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 

       case('C2H','C2H(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) =-src(6)
           !
         case (4) ! (12)(34)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) =-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       case('C2V','C2V(M)')
         !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) =-src(6)
           !
         case (4) ! (12)(34)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) =-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       end select
        !
    case('R-ALPHA-TAU-REF','R-ALPHA-TAU')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_abcd - bad symm. type'
          !
       case('C','C(M)')
         !
         dst = src
         !
       case('CS','CS(M)')
         !
         select case(ioper)
         !
         case (1) ! identity 

           dst = src

         case (2) ! (E*)

           dst = src
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
          !
       case('G4','G4(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6)
           !
           if (src(6)>2.0_ark*pi) dst(6) = 4.0_ark*pi-src(6)+4.0_ark*pi

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6)
           !
           if (src(6)>2.0_ark*pi) dst(6) = 4.0_ark*pi-src(6)+4.0_ark*pi

         case (4) ! (12)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       case('CS(EM)')
           !
         select case(ioper)
           !
         case (1) ! E 
           !
           dst = src
           !
         case (2) ! (C2z) -> Rz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi+src(6) 
           if (extended) then
              if(dst(6)>4.0_ark*pi) dst(6) = dst(6) - 4.0_ark*pi
           else
              dst(6) = src(6)
           endif
           !
         case (3) ! E* -> sigma_xz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 4.0_ark*pi-src(6)
           if (extended) then 
              dst(6) = dst(6)
           else
              dst(6) = 2.0_ark*pi-src(6)
           endif
           !
         case (4) ! sigma_yz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6) !; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           if (extended) then 
              if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           else
              dst(6) = dst(6)
           endif
           !
         case default
           !
           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
           !
         end select 
           !
       case('C2V','C2V(M)')
           !
         select case(ioper)
           !
         case (1) ! E 
           !
           dst = src
           !
         case (2) ! (12)(34)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6)

         case (4) ! (12)(34)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
          !
       case('C2H','C2H(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)
           !dst(6) = 4.0_ark*pi-src(6)
           !if (src(6)>2.0_ark*pi) then 
           ! dst(6) = src(6)+2.0_ark*pi
           !endif

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 4.0_ark*pi-src(6) !; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           !dst(6) = 2.0_ark*pi-src(6) ; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           !
           !dst(6) = 2.0_ark*pi-src(6)
           !if (src(6)>2.0_ark*pi) dst(6) = 6.0_ark*pi-src(6)
           !
           !if (src(6)>2.0_ark*pi) then 
            !dst(1) = src(1)
            !dst(2) = src(3)
            !dst(3) = src(2)
            !dst(4) = src(5)
            !dst(5) = src(4)
           ! dst(6) = 6.0_ark*pi-src(6)
           !endif
           !
         case (4) ! (12)(34)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           !dst(6) = 2.0_ark*pi-src(6) ; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           dst(6) = 4.0_ark*pi-src(6) ! ; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           ! 
           !dst(6) = src(6)
           !dst(6) = 2.0_ark*pi-src(6)
           !if (src(6)>2.0_ark*pi) dst(6) = 6.0_ark*pi-src(6)

           !
           !if (src(6)>2.0_ark*pi) then 
            !dst(1) = src(1)
            !dst(2) = src(2)
            !dst(3) = src(3)
            !dst(4) = src(4)
            !dst(5) = src(5)
            !dst(6) = 6.0_ark*pi-src(6)
            !dst(6) = src(6)-2.0_ark*pi
           !endif




         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 

          !
       case('G4(EM)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34) <-> a

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case (3) ! (E*) <->  b

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 4.0_ark*pi-src(6)

         case (4) ! (12)(34)* <->  ab 

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
          dst(6) = 4.0_ark*pi-src(6)

         case (5) ! E'

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi+src(6) ; if (dst(6)>0) dst(6) = dst(6) - 4.0_ark*pi

         case (6) ! E'a

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi+src(6) ; if (dst(6)>0) dst(6) = dst(6) - 4.0_ark*pi

         case (7) ! E'b

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6) ; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi

         case (8) ! E'ab 

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6) ; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       case('D2H(M)')
           !
         select case(ioper)
           !
         case (1) ! E 
           !
           dst = src
           !
         case (2) ! (C2z) -> Rz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi+src(6) 
           if (extended) then
              if(dst(6)>4.0_ark*pi) dst(6) = dst(6) - 4.0_ark*pi
           else
              dst(6) = src(6)
           endif
           !
         case (3) ! (C2x)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)
           !
         case (4) ! (C2y)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi+src(6) ! ; if (dst(6)>0) dst(6) = dst(6) - 4.0_ark*pi
           if (extended) then 
              if(dst(6)>4.0_ark*pi) dst(6) = dst(6) - 4.0_ark*pi
           else
              dst(6) = src(6)
           endif
           !
         case (5) ! E* -> sigma_xz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 4.0_ark*pi-src(6)
           if (extended) then 
              dst(6) = dst(6)
           else
              dst(6) = 2.0_ark*pi-src(6)
           endif

         case (6) ! sigma_yz

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6) !; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           if (extended) then 
              if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           else
              dst(6) = dst(6)
           endif

         case (7) ! sigma_xy

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 4.0_ark*pi-src(6)
           if (extended) then 
              dst(6) = dst(6)
           else
              dst(6) = 2.0_ark*pi-src(6)
           endif

         case (8) ! i

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6) !; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           if (extended) then 
              if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           else
              dst(6) = dst(6)
           endif

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       end select
        !
        !
    case('R-ALPHA-R-TAU')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_abcd - bad symm. type'
          !
       case('C','C(M)')
         !
         dst = src
         !
       case('CS','CS(M)')
         !
         select case(ioper)
         !
         case (1) ! identity 

           dst = src

         case (2) ! (E*)

           dst = src
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
          !
       case('G4','G4(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6)
           !
           if (src(6)>2.0_ark*pi) dst(6) = 4.0_ark*pi-src(6)+4.0_ark*pi

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6)
           !
           if (src(6)>2.0_ark*pi) dst(6) = 4.0_ark*pi-src(6)+4.0_ark*pi

         case (4) ! (12)*

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
          !
       case('C2V','C2V(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6)

         case (4) ! (12)(34)*

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
          !
       case('C2H','C2H(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6)

         case (4) ! (12)(34)*

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 

          !
       case('G4(EM)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34) <-> a

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case (3) ! (E*) <->  b

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6)

         case (4) ! (12)(34)* <->  ab 

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6)

         case (5) ! E'

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 4.0_ark*pi-src(6)

         case (6) ! E'a

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 4.0_ark*pi-src(6)

         case (7) ! E'b

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi+src(6)

         case (8) ! E'ab 

           dst(1) = src(2)
           dst(2) = src(1)
           dst(3) = src(3)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi+src(6)


         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       end select
       !
    case('R-S-TAU-REF','R-S-TAU')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_abcd - bad symm. type'
          !
       case('CS','CS(M)')
         !
         select case(ioper)
         !
         case (1) ! identity 

           dst = src

         case (2) ! (E*)

           dst = src
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       case('C2H','C2H(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) =-src(3)
           dst(4) = src(4)
           dst(5) =-src(5)
           dst(6) = src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6)

         case (4) ! (12)(34)*

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) =-src(3)
           dst(4) = src(4)
           dst(5) =-src(5)
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       end select
       !
    case('R-R1-R2-T1-T2-T3-T4')
       ! 
       select case(trim(molec%symmetry))
       case default
         write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
         stop 'ML_symmetry_transformation_abcd - bad symm. type'
         !
       case('CS','CS(M)')
          !
          dst = src
          !
          N_Cn = sym%N/2
          Nrot = sym%N
          !
          if (ioper>=1+2*N_Cn+Nrot+1) then
            !
            dst(1) = src(1)
            dst(2) = src(3)
            dst(3) = src(2)
            dst(6) = src(4)
            dst(7) = src(5)
            dst(4) = src(6)
            dst(5) = src(7)
            !
          endif 
          !
       case('C2V','C2V(M)')
          !
          select case(ioper)
           !
          case (1) ! E 
           !
           dst = src
           !
          case (2) ! (12)(34)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(6)
           dst(5) =-src(7)
           dst(6) = src(4)
           dst(7) =-src(5)
           !
          case (3) ! (E*)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(6)
           dst(5) = src(7)
           dst(6) = src(4)
           dst(7) = src(5)
           !
          case (4) ! (12)(34)*
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) =-src(5)
           dst(6) = src(6)
           dst(7) =-src(7)
           !
          case default
           !
           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
            !
          end select 
           !
       case('D2H(M)')
           !
         select case(ioper)
           !
         case (1) ! E 
           !
           dst = src
           !
         case (2) ! (C2z) -> Rz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) =-src(4)
           dst(5) =-src(5)
           dst(6) =-src(6)
           dst(7) = src(7)
           !
         case (3) ! (C2x)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) =-src(6)
           dst(5) = src(7)
           dst(6) =-src(4)
           dst(7) = src(5)
           !
         case (4) ! (C2y)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(6)
           dst(5) =-src(7)
           dst(6) = src(4)
           dst(7) =-src(5)
           !
         case (5) ! E* -> sigma_xz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) =-src(4)
           dst(5) = src(5)
           dst(6) =-src(6)
           dst(7) = src(7)

         case (6) ! sigma_yz

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) =-src(5)
           dst(6) = src(6)
           dst(7) =-src(7)

         case (7) ! sigma_xy

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(6)
           dst(5) = src(7)
           dst(6) = src(4)
           dst(7) = src(5)

         case (8) ! i

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) =-src(6)
           dst(5) =-src(7)
           dst(6) =-src(4)
           dst(7) =-src(5)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
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
          q1x= src(4)
          q1y= src(5)
          q2x= src(6)
          q2y= src(7)
          !
          ! for odd Dnh groups
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
                dst(3) = src(3)
                dst(4) = q1x*cos(phi_n)-q1y*sin(phi_n)
                dst(5) = q1x*sin(phi_n)+q1y*cos(phi_n)
                dst(6) = q2x*cos(phi_n)-q2y*sin(phi_n)
                dst(7) = q2x*sin(phi_n)+q2y*cos(phi_n)
                !
              elseif (ioper<=1+2*N_Cn+Nrot) then !  C'2
                !
                irot = ioper-(1+2*N_Cn)
                !
                phi_n = phi*irot
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                !
              elseif (ioper==1+2*N_Cn+Nrot+1) then ! sigmah
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = src(4)
                dst(7) = src(5)
                dst(4) = src(6)
                dst(5) = src(7)
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
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = (q1x*cos(phi_n)-q1y*sin(phi_n))
                dst(7) = (q1x*sin(phi_n)+q1y*cos(phi_n))
                dst(4) = (q2x*cos(phi_n)-q2y*sin(phi_n))
                dst(5) = (q2x*sin(phi_n)+q2y*cos(phi_n))
                !
              elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then !  sigmav
                !
                irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
                !
                phi_n = phi*irot
                !
                dst(1) = src(1)
                dst(2) = src(2)
                dst(3) = src(3)
                dst(4) = ( q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) = ( q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) = ( q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) = ( q2x*sin(phi_n)-q2y*cos(phi_n))
                !
              else
                !
                write (out,"('ML_symmetry_transformation_abcd  in Dinfty: operation ',i8,' unknown')") ioper
                stop 'ML_symmetry_transformation_abcd Dinfty - bad operation. type'
         
              endif 
              !
              !do ioper_=1,sym%Noper/2
              !  repres ( ioper_,            :,:)= sym%irr(5, ioper_)%repres
              !  repres ( ioper_+sym%Noper/2,:,:)= sym%irr(6, ioper_)%repres
              !enddo
              !
              !do ioper_=1,sym%Noper
              !  repres ( ioper_,            :,:)= sym%irr(5, ioper_)%repres
              !  !repres ( ioper_+sym%Noper/2,:,:)= sym%irr(5, ioper_)%repres
              !enddo
              !
              !dst(4) = repres(ioper,1,1)*src(4)+repres(ioper,1,2)*src(5)
              !dst(5) = repres(ioper,2,1)*src(4)+repres(ioper,2,2)*src(5)
              !dst(6) = repres(ioper,1,1)*src(6)+repres(ioper,1,2)*src(7)
              !dst(7) = repres(ioper,2,1)*src(6)+repres(ioper,2,2)*src(7)
              !
          else ! for even Dnh groups
              !
              !q1x= 1
              !q1y= 1
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
                dst(3) = src(3)
                dst(4) = q1x*cos(phi_n)-q1y*sin(phi_n)
                dst(5) = q1x*sin(phi_n)+q1y*cos(phi_n)
                dst(6) = q2x*cos(phi_n)-q2y*sin(phi_n)
                dst(7) = q2x*sin(phi_n)+q2y*cos(phi_n)
                !
              elseif (ioper<=2+2*N_Cn+NC2) then !  C'2
                !
                irot =ioper-(2+2*N_Cn)-1
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                !
              elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
                !
                irot =ioper-(2+2*N_Cn+NC2)-1
                !
                !phi_n =(irot*2.0_ark-1.0_ark)*phi*0.5_ark
                !
                phi_n = phi*(2*irot+1)
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                !
              elseif (ioper==3+2*N_Cn+2*NC2) then ! i
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) =-q2x
                dst(5) =-q2y
                dst(6) =-q1x
                dst(7) =-q1y
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
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) =  (q2x*cos(phi_n)-q2y*sin(phi_n))
                dst(5) =  (q2x*sin(phi_n)+q2y*cos(phi_n))
                dst(6) =  (q1x*cos(phi_n)-q1y*sin(phi_n))
                dst(7) =  (q1x*sin(phi_n)+q1y*cos(phi_n))
                !
              elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigmah
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x
                dst(5) = q2y
                dst(6) = q1x
                dst(7) = q1y
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
                dst(4) =  (q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) =  (q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) =  (q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) =  (q2x*sin(phi_n)-q2y*cos(phi_n))
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
                dst(3) = src(3)
                dst(4) =  (q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) =  (q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) =  (q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) =  (q2x*sin(phi_n)-q2y*cos(phi_n))
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
    case('R-R1-R2-TX-TY-TX-TY','R-R1-R2+TX-TY+TX-TY','R-R1-R2+Y-X-Y+X','R-R1-R2-Y+X+Y-X','R-R1-R2-X-Y+X+Y',&
         'R-Z1+Z2-X-Y+X+Y','R-Z1+Z2+X+Y-X-Y','R-Z1+Z2-Y-X-Y-X','R-R1-R2+X+Y-X-Y','R-R1-R2-TX+TY-TX+TY',&
         'R-R1-R2-TX+TY+TX-TY')
       ! 
       select case(trim(molec%symmetry))
       case default
         write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
         stop 'ML_symmetry_transformation_abcd - bad symm. type'
         !
       case('CS','CS(M)')
          !
          dst = src
          N_Cn = sym%N/2
          Nrot = sym%N
          !
          if (ioper>=1+2*N_Cn+Nrot+1) then
            !
            dst(1) = src(1)
            dst(2) = src(3)
            dst(3) = src(2)
            dst(6) = src(4)
            dst(7) = src(5)
            dst(4) = src(6)
            dst(5) = src(7)
            !
          endif 
          !
       case('C2V','C2V(M)')
          !
          select case(ioper)
           !
          case (1) ! E 
           !
           dst = src
           !
          case (2) ! (12)(34)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(6)
           dst(5) =-src(7)
           dst(6) = src(4)
           dst(7) =-src(5)
           !
          case (3) ! (E*)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(6)
           dst(5) = src(7)
           dst(6) = src(4)
           dst(7) = src(5)
           !
          case (4) ! (12)(34)*
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) =-src(5)
           dst(6) = src(6)
           dst(7) =-src(7)
           !
          case default
           !
           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
            !
          end select 
           !
       case('D2H(M)')
           !
         select case(ioper)
           !
         case (1) ! E 
           !
           dst = src
           !
         case (2) ! (C2z) -> Rz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) =-src(4)
           dst(5) =-src(5)
           dst(6) =-src(6)
           dst(7) = src(7)
           !
         case (3) ! (C2x)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) =-src(6)
           dst(5) = src(7)
           dst(6) =-src(4)
           dst(7) = src(5)
           !
         case (4) ! (C2y)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(6)
           dst(5) =-src(7)
           dst(6) = src(4)
           dst(7) =-src(5)
           !
         case (5) ! E* -> sigma_xz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) =-src(4)
           dst(5) = src(5)
           dst(6) =-src(6)
           dst(7) = src(7)

         case (6) ! sigma_yz

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) =-src(5)
           dst(6) = src(6)
           dst(7) =-src(7)

         case (7) ! sigma_xy

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(6)
           dst(5) = src(7)
           dst(6) = src(4)
           dst(7) = src(5)

         case (8) ! i

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) =-src(6)
           dst(5) =-src(7)
           dst(6) =-src(4)
           dst(7) =-src(5)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
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
          q1x= src(4)
          q1y= src(5)
          q2x= src(6)
          q2y= src(7)
          !
          ! for odd Dnh groups
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
                dst(3) = src(3)
                dst(4) = q1x*cos(phi_n)-q1y*sin(phi_n)
                dst(5) = q1x*sin(phi_n)+q1y*cos(phi_n)
                dst(6) = q2x*cos(phi_n)-q2y*sin(phi_n)
                dst(7) = q2x*sin(phi_n)+q2y*cos(phi_n)
                !
              elseif (ioper<=1+2*N_Cn+Nrot) then !  C'2
                !
                irot = ioper-(1+2*N_Cn)
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                !
              elseif (ioper==1+2*N_Cn+Nrot+1) then ! sigmah
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = src(4)
                dst(7) = src(5)
                dst(4) = src(6)
                dst(5) = src(7)
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
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = (q1x*cos(phi_n)-q1y*sin(phi_n))
                dst(7) = (q1x*sin(phi_n)+q1y*cos(phi_n))
                dst(4) = (q2x*cos(phi_n)-q2y*sin(phi_n))
                dst(5) = (q2x*sin(phi_n)+q2y*cos(phi_n))
                !
              elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then !  sigmav
                !
                irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(1) = src(1)
                dst(2) = src(2)
                dst(3) = src(3)
                dst(4) = ( q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) = ( q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) = ( q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) = ( q2x*sin(phi_n)-q2y*cos(phi_n))
                !
              else
                !
                write (out,"('ML_symmetry_transformation_abcd  in Dinfty: operation ',i8,' unknown')") ioper
                stop 'ML_symmetry_transformation_abcd Dinfty - bad operation. type'
         
              endif 
              !
          else ! for even Dnh groups
              !
              !q1x= 1
              !q1y= 1
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
                dst(3) = src(3)
                dst(4) = q1x*cos(phi_n)-q1y*sin(phi_n)
                dst(5) = q1x*sin(phi_n)+q1y*cos(phi_n)
                dst(6) = q2x*cos(phi_n)-q2y*sin(phi_n)
                dst(7) = q2x*sin(phi_n)+q2y*cos(phi_n)
                !
              elseif (ioper<=2+2*N_Cn+NC2) then !  C'2
                !
                irot =ioper-(2+2*N_Cn)-1
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                !
              elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
                !
                irot =ioper-(2+2*N_Cn+NC2)-1
                !
                !phi_n =(irot*2.0_ark-1.0_ark)*phi*0.5_ark
                !
                phi_n = phi*(2*irot+1)
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                !
              elseif (ioper==3+2*N_Cn+2*NC2) then ! i
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) =-q2x
                dst(5) =-q2y
                dst(6) =-q1x
                dst(7) =-q1y
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
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) =  (q2x*cos(phi_n)-q2y*sin(phi_n))
                dst(5) =  (q2x*sin(phi_n)+q2y*cos(phi_n))
                dst(6) =  (q1x*cos(phi_n)-q1y*sin(phi_n))
                dst(7) =  (q1x*sin(phi_n)+q1y*cos(phi_n))
                !
              elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigmah
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x
                dst(5) = q2y
                dst(6) = q1x
                dst(7) = q1y
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
                dst(4) =  (q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) =  (q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) =  (q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) =  (q2x*sin(phi_n)-q2y*cos(phi_n))
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
                dst(3) = src(3)
                dst(4) =  (q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) =  (q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) =  (q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) =  (q2x*sin(phi_n)-q2y*cos(phi_n))
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
    case('R1-R2-Y+X+Y-X-R')
       ! 
       select case(trim(molec%symmetry))
       case default
         write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
         stop 'ML_symmetry_transformation_abcd - bad symm. type'
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
          q1x= src(3)
          q1y= src(4)
          q2x= src(5)
          q2y= src(6)
          !
          ! for odd Dnh groups
          !
          if (mod(sym%N,2)==1) then 
              !
              write(out,"('R1-R2-Y+X+Y-X-R: Dnh for odd n is not implemented yet')")
              stop 'R1-R2-Y+X+Y-X-R: Dnh for odd n is not implemented yet'
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
                dst(3) = src(3)
                dst(4) = q1x*cos(phi_n)-q1y*sin(phi_n)
                dst(5) = q1x*sin(phi_n)+q1y*cos(phi_n)
                dst(6) = q2x*cos(phi_n)-q2y*sin(phi_n)
                dst(7) = q2x*sin(phi_n)+q2y*cos(phi_n)
                !
              elseif (ioper<=1+2*N_Cn+Nrot) then !  C'2
                !
                irot = ioper-(1+2*N_Cn)
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                !
              elseif (ioper==1+2*N_Cn+Nrot+1) then ! sigmah
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = src(4)
                dst(7) = src(5)
                dst(4) = src(6)
                dst(5) = src(7)
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
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(6) = (q1x*cos(phi_n)-q1y*sin(phi_n))
                dst(7) = (q1x*sin(phi_n)+q1y*cos(phi_n))
                dst(4) = (q2x*cos(phi_n)-q2y*sin(phi_n))
                dst(5) = (q2x*sin(phi_n)+q2y*cos(phi_n))
                !
              elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then !  sigmav
                !
                irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(1) = src(1)
                dst(2) = src(2)
                dst(3) = src(3)
                dst(4) = ( q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) = ( q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) = ( q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) = ( q2x*sin(phi_n)-q2y*cos(phi_n))
                !
              else
                !
                write (out,"('ML_symmetry_transformation_abcd  in Dinfty: operation ',i8,' unknown')") ioper
                stop 'ML_symmetry_transformation_abcd Dinfty - bad operation. type'
         
              endif 
              !
          else ! for even Dnh groups
              !
              !q1x= 1
              !q1y= 1
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
                dst(7) = src(7)
                dst(1) = src(1)
                dst(2) = src(2)
                dst(3) = q1x*cos(phi_n)-q1y*sin(phi_n)
                dst(4) = q1x*sin(phi_n)+q1y*cos(phi_n)
                dst(5) = q2x*cos(phi_n)-q2y*sin(phi_n)
                dst(6) = q2x*sin(phi_n)+q2y*cos(phi_n)
                !
              elseif (ioper<=2+2*N_Cn+NC2) then !  C'2
                !
                irot =ioper-(2+2*N_Cn)-1
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(7) = src(7)
                dst(1) = src(2)
                dst(2) = src(1)
                dst(3) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(4) = q2x*sin(phi_n)-q2y*cos(phi_n)
                dst(5) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(6) = q1x*sin(phi_n)-q1y*cos(phi_n)
                !
              elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
                !
                irot =ioper-(2+2*N_Cn+NC2)-1
                !
                !phi_n =(irot*2.0_ark-1.0_ark)*phi*0.5_ark
                !
                phi_n = phi*(2*irot+1)
                !
                dst(7) = src(7)
                dst(1) = src(2)
                dst(2) = src(1)
                dst(3) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(4) = q2x*sin(phi_n)-q2y*cos(phi_n)
                dst(5) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(6) = q1x*sin(phi_n)-q1y*cos(phi_n)
                !
              elseif (ioper==3+2*N_Cn+2*NC2) then ! i
                !
                dst(7) = src(7)
                dst(1) = src(2)
                dst(2) = src(1)
                dst(3) =-q2x
                dst(4) =-q2y
                dst(5) =-q1x
                dst(6) =-q1y
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
                dst(7) = src(7)
                dst(1) = src(2)
                dst(2) = src(1)
                dst(3) =  (q2x*cos(phi_n)-q2y*sin(phi_n))
                dst(4) =  (q2x*sin(phi_n)+q2y*cos(phi_n))
                dst(5) =  (q1x*cos(phi_n)-q1y*sin(phi_n))
                dst(6) =  (q1x*sin(phi_n)+q1y*cos(phi_n))
                !
              elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigmah
                !
                dst(7) = src(7)
                dst(1) = src(2)
                dst(2) = src(1)
                dst(3) = q2x
                dst(4) = q2y
                dst(5) = q1x
                dst(6) = q1y
                !
              elseif (ioper<=4+4*N_Cn+3*NC2) then !  sigmav
                !
                irot = ioper-(4+4*N_Cn+2*NC2)-1
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(7) = src(7)
                dst(1) = src(1)
                dst(2) = src(2)
                dst(3) =  (q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(4) =  (q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(5) =  (q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(6) =  (q2x*sin(phi_n)-q2y*cos(phi_n))
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
                dst(7) = src(7)
                dst(1) = src(1)
                dst(2) = src(2)
                dst(3) =  (q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(4) =  (q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(5) =  (q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(6) =  (q2x*sin(phi_n)-q2y*cos(phi_n))
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
    case('R-R1-R2-X-Y-X-Y','R-Z1-Z2-X-Y-X-Y','R-R1-R2+X+Y-X+Y')
       ! 
       select case(trim(molec%symmetry))
       case default
         write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
         stop 'ML_symmetry_transformation_abcd - bad symm. type'
         !
       case('CS','CS(M)')
          !
          dst = src
          !
          N_Cn = sym%N/2
          Nrot = sym%N
          !
          if (ioper>=1+2*N_Cn+Nrot+1) then
            !
            dst(1) = src(1)
            dst(2) = src(3)
            dst(3) = src(2)
            dst(6) = src(4)
            dst(7) = src(5)
            dst(4) = src(6)
            dst(5) = src(7)
            !
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
          q1x= src(4)
          q1y= src(5)
          q2x= src(6)
          q2y= src(7)
          !
          ! for odd Dnh groups
          !
          if (mod(sym%N,2)==1) then 
              !
              write(out,"('symmetries odd-DNH for R-R1-R2-X-Y-X-Y not implemented')")
              !
              stop 'odd DNH for R-R1-R2-X-Y-X-Y not implemented'
              !
          else ! for even Dnh groups
              !
              !q1x= 1
              !q1y= 1
              !q2x= 2
              !q2y= 2
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
                dst(3) = src(3)
                dst(4) = q1x*cos(phi_n)-q1y*sin(phi_n)
                dst(5) = q1x*sin(phi_n)+q1y*cos(phi_n)
                dst(6) = q2x*cos(phi_n)-q2y*sin(phi_n)
                dst(7) = q2x*sin(phi_n)+q2y*cos(phi_n)
                !
              elseif (ioper<=2+2*N_Cn+NC2) then !  C'2
                !
                irot =ioper-(2+2*N_Cn)-1
                !
                phi_n = phi*irot*2.0_ark
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                !
              elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
                !
                irot =ioper-(2+2*N_Cn+NC2)-1
                !
                !phi_n =(irot*2.0_ark-1.0_ark)*phi*0.5_ark
                !
                phi_n = phi*(2*irot+1)
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x*cos(phi_n)+q2y*sin(phi_n)
                dst(5) = q2x*sin(phi_n)-q2y*cos(phi_n)
                dst(6) = q1x*cos(phi_n)+q1y*sin(phi_n)
                dst(7) = q1x*sin(phi_n)-q1y*cos(phi_n)
                !
              elseif (ioper==3+2*N_Cn+2*NC2) then ! i
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) =-q2x
                dst(5) =-q2y
                dst(6) =-q1x
                dst(7) =-q1y
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
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) =  (q2x*cos(phi_n)-q2y*sin(phi_n))
                dst(5) =  (q2x*sin(phi_n)+q2y*cos(phi_n))
                dst(6) =  (q1x*cos(phi_n)-q1y*sin(phi_n))
                dst(7) =  (q1x*sin(phi_n)+q1y*cos(phi_n))
                !
              elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigmah
                !
                dst(1) = src(1)
                dst(2) = src(3)
                dst(3) = src(2)
                dst(4) = q2x
                dst(5) = q2y
                dst(6) = q1x
                dst(7) = q1y
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
                dst(4) =  (q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) =  (q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) =  (q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) =  (q2x*sin(phi_n)-q2y*cos(phi_n))
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
                dst(3) = src(3)
                dst(4) =  (q1x*cos(phi_n)+q1y*sin(phi_n))
                dst(5) =  (q1x*sin(phi_n)-q1y*cos(phi_n))
                dst(6) =  (q2x*cos(phi_n)+q2y*sin(phi_n))
                dst(7) =  (q2x*sin(phi_n)-q2y*cos(phi_n))
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
    case('R-PHI-TAU-REF','R-PHI-TAU')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_abcd - bad symm. type'
          !
       case('CS','CS(M)')
         !
         select case(ioper)
         !
         case (1) ! identity 

           dst = src

         case (2) ! (E*)

           dst = src
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       case('C2H','C2H(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) =-src(5)
           dst(5) =-src(4)
           dst(6) = src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6)

         case (4) ! (12)(34)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) =-src(5)
           dst(5) =-src(4)
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_abcd: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_abcd - bad operation. type'
 
         end select 
         !
       end select
       !
    end select
    !
    if (verbose>=7) write(out,"('ML_symmetry_transformation_abcd/end')") 
    !
    !
  end subroutine ML_symmetry_transformation_abcd


  ! Here we find the rotational symmetry from J,K
  !
  subroutine ML_rotsymmetry_abcd(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
    integer(ik),intent(out) :: gamma,ideg
    !
    integer(ik) :: N,N_Cn,K_,L
    !
    if (verbose>=7) write(out,"('ML_rotsymmetry_abcd/start')") 
    !
    !
    select case(trim(molec%symmetry))
    case default
       write (out,"('ML_rotsymmetry_abcd: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_rotsymmetry_abcd - bad symm. type'
       !
    case('C','C(M)')
       !
       gamma = 1
       ideg = 1 
       !
    case('CS','CS(M)')
       !
       gamma = 0 
       ideg = 1 
       !
       if (mod(K+tau+2,2)==0) gamma = 1 !; return
       if (mod(K+tau+2,2)/=0) gamma = 2 !; return
       !
    case('G4','G4(M)')
       !
       gamma = 0
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma = 1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma = 3 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma = 4 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma = 2 !; return
       !
    case('C2H(M)','C2H')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma =1 !1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma =3 !3 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma =3 !4 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma =1 !2 !; return
       !
    case('CS(EM)')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma =1 !1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma =2 !3 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma =4 !4 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma =3 !2 !; return
       !
    case('C2V','C2V(M)')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma = 1 !1 !1 !1 !1 !1 !1 ; return
       if (mod(K+2,2)==0.and.tau==1) gamma = 3 !2 !4 !4 !3 !3 !2 ; return
       if (mod(K+2,2)/=0.and.tau==0) gamma = 4 !4 !3 !2 !4 !2 !3 ; return
       if (mod(K+2,2)/=0.and.tau==1) gamma = 2 !3 !2 !3 !2 !4 !4 ; return
       !
    case('G4(EM)')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma = 1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma = 2 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma = 3 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma = 4 !; return
       !
       !if (mod(K+2,2)==0.and.tau==0) gamma = 1 !; return
       !if (mod(K+2,2)==0.and.tau==1) gamma = 2 !; return
       !if (mod(K+2,2)/=0.and.tau==0) gamma = 3 !; return
       !if (mod(K+2,2)/=0.and.tau==1) gamma = 4 !; return
       !
    case('D2H(M)')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma = 1 !1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma = 3 !3 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma = 5 !7 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma = 7 !5 !; return
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
    !
    if (verbose>=7) write(out,"('ML_rotsymmetry_abcd/end')") 
    !
    !
  end subroutine ML_rotsymmetry_abcd



end module mol_abcd

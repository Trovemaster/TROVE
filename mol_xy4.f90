!
!  This unit defines all specific routines for a fiveatomic molecule of XY4 type
!
module mol_xy4
  use accuracy
  use moltype
  use lapack
  use symmetry,     only : sym
  use pot_xy4, only : ML_XY4_calc_alpha34

  implicit none

  public ML_b0_XY4,ML_coordinate_transform_XY4,ML_symmetry_transformation_XY4,ML_rotsymmetry_XY4
  !public MLpoten_xy4_Bowman2000,MLpoten_xy4_ZZZ

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains

  ! Roman
  !
  ! Here we define structural parameters for rigid XY4 molecule,
  ! a0 and Amat,
  ! which determine the normal coordinates 
  !
  subroutine ML_b0_XY4(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)

     integer(ik),intent(in)  :: Npoints,Natoms

     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in),optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional :: rho_ref
     real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders

     integer(ik)             :: i
     !
     real(ark)                :: re14,mH,mX
     real(ark)                :: alphae,alpha2,cosae
     
      if (verbose>=4) write(out,"(/'ML_b0_XY4/start')") 


      if (size(molec%req)/=4) then
        write(out,"('ML_b0_XY4: Nbonds must be 4 in this routine, not  ',i)") size(molec%req)
        stop 'ML_b0_XY4: wrong Nbonds '
      endif 

      if (molec%Natoms/=5) then
        write(out,"('ML_b0_XY4: Natoms must be 5 in this routine, not  ',i)") molec%Natoms
        stop 'ML_b0_XY4: wrong Natoms '
      endif 

      if (molec%req(1)/=molec%req(2).or.molec%req(1)/=molec%req(3).or.molec%req(1)/=molec%req(4)) then
        write(out,"('ML_b0_XY4: req-s must be equal: ',4f14.6)") molec%req(1:4)
        stop 'ML_b0_XY4: req-s must be equal: '
      endif 
      !
      re14 = molec%req(1)
      alphae = molec%alphaeq(1)
      alpha2 = alphae*0.5_ark
      if (any(molec%alphaeq(:)/=alphae)) then
        write(out,"('ML_b0_XY4: alphaeq-s must be equal: ',4f14.6)") molec%alphaeq(:)
        stop 'ML_b0_XY4: alphaeq-s must be equal: '
      endif 
      !
      mH = molec%AtomMasses(2) ; mX = molec%AtomMasses(1)
      !
      if (any(molec%AtomMasses(2:5)/=mH)) then
        write(out,"('ML_b0_XY4: masses-s are given in wrong order, must be M m m m m: ',5f14.6)") molec%AtomMasses(:)
        stop 'ML_b0_XY4: ,masses are in wrong order'
      endif 
      !
      cosae = -1.0_ark/3.0_ark
      molec%local_eq(5:10) = acos(cosae)


    select case(trim(molec%coords_transform))
    case default
       !
       !
       ! 
       !b0(1,1,0) = 0.0_ark
       !b0(1,2,0) = 0.0_ark
       !b0(1,3,0) = 0.0_ark
       !b0(2,1,0) = re14/sqrt3
       !b0(2,2,0) = re14/sqrt3
       !b0(2,3,0) = re14/sqrt3
       !b0(3,1,0) =-re14/sqrt3
       !b0(3,2,0) = re14/sqrt3
       !b0(3,3,0) =-re14/sqrt3
       !b0(4,1,0) = re14/sqrt3
       !b0(4,2,0) =-re14/sqrt3
       !b0(4,3,0) =-re14/sqrt3
       !b0(5,1,0) =-re14/sqrt3
       !b0(5,2,0) =-re14/sqrt3
       !b0(5,3,0) = re14/sqrt3
       !
       b0(1,1,0) = 0.0_ark
       b0(1,2,0) = 0.0_ark
       b0(1,3,0) = 0.0_ark
       b0(2,1,0) =-re14/sqrt3
       b0(2,2,0) = re14/sqrt3
       b0(2,3,0) = re14/sqrt3
       b0(3,1,0) =-re14/sqrt3
       b0(3,2,0) =-re14/sqrt3
       b0(3,3,0) =-re14/sqrt3
       b0(4,1,0) = re14/sqrt3
       b0(4,2,0) = re14/sqrt3
       b0(4,3,0) =-re14/sqrt3
       b0(5,1,0) = re14/sqrt3
       b0(5,2,0) =-re14/sqrt3
       b0(5,3,0) = re14/sqrt3
       !
    case('R-ALPHA')
       !
       !
       b0(1,1,0) = 0.0_ark
       b0(1,2,0) = 0.0_ark
       b0(1,3,0) = 0.0_ark
       b0(2,1,0) = re14*sin(alpha2)
       b0(2,2,0) = 0
       b0(2,3,0) = re14*cos(alpha2)
       b0(3,1,0) =-re14*sin(alpha2)
       b0(3,2,0) = 0
       b0(3,3,0) = re14*cos(alpha2)
       b0(4,1,0) = 0
       b0(4,2,0) = re14*sin(alpha2)
       b0(4,3,0) =-re14*cos(alpha2)
       b0(5,1,0) = 0
       b0(5,2,0) =-re14*sin(alpha2)
       b0(5,3,0) =-re14*cos(alpha2)
       !
    end select



      !b0(1,1,0) = 0.0_ark
      !b0(1,2,0) = 0.0_ark
      !b0(1,3,0) = 0.0_ark
      !b0(4,1,0) = re14/sqrt3
      !b0(4,2,0) =-re14/sqrt3
      !b0(4,3,0) =-re14/sqrt3
      !b0(5,1,0) =-re14/sqrt3
      !b0(5,2,0) = re14/sqrt3
      !b0(5,3,0) =-re14/sqrt3
      !b0(2,1,0) = re14/sqrt3
      !b0(2,2,0) = re14/sqrt3
      !b0(2,3,0) = re14/sqrt3
      !b0(3,1,0) =-re14/sqrt3
      !b0(3,2,0) =-re14/sqrt3
      !b0(3,3,0) = re14/sqrt3

! Roman CH4
!      b0(1,1,0) = 0.0_rk
!      b0(1,2,0) = 0.0_rk
!      b0(1,3,0) = 0.0_rk
!      b0(2,1,0) = re14
!      b0(2,2,0) = 0.0_rk
!      b0(2,3,0) = 0.0_rk
!      b0(3,1,0) = re14*cos(alphae)
!      b0(3,2,0) = re14*sin(alphae)
!      b0(3,3,0) = 0.0_rk
!      b0(4,1,0) = b0(3,1,0)
!      b0(4,2,0) = re14*cos(alphae)*(1.0_rk-cos(alphae))/sin(alphae)
!      b0(4,3,0) = re14*sqrt(sin(alphae)**2-(cos(alphae)*(1.0_rk-cos(alphae))/sin(alphae))**2)
!      b0(5,1,0) = b0(4,1,0)
!      b0(5,2,0) = b0(4,2,0)
!      b0(5,3,0) = -b0(4,3,0)


      !
! Roman end
      ! We can also simulate the reference structure at each point of rho, if npoints presented
      !
      if (Npoints/=0) then
         !
         !rho_i = 0
         !js947: comment this because it messes up the interface and is not
         !important if the program is about to stop
         ! 
         write(out,"('ML_b0_XY4: non-rigid  is not implemented yet')") 
         stop 'ML_b0_XY4: non-rigid  is not implemented yet '
         !
         if (.not.present(rho_borders).or..not.present(rho_ref)) then  
            write(out,"('ML_b0_XY4: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_XY4: rho_borders or rho_ref not specified '
         endif
         !
         rho_ref = 0 
         !
      endif
      !
      if (verbose>=3) then 
         !
         write(out,"(i6)") molec%natoms
         !
         write(out,"(/'C',3x,3f14.8)") b0(1,:,0)
         write(out,"( 'H',3x,3f14.8)") b0(2,:,0)
         write(out,"( 'H',3x,3f14.8)") b0(3,:,0)
         write(out,"( 'H',3x,3f14.8)") b0(4,:,0)
         write(out,"( 'H',3x,3f14.8)") b0(5,:,0)
         !
      endif
      !
      if (verbose>=4) write(out,"('ML_b0_XY4/end')") 

  end subroutine ML_b0_XY4


  !
  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  function ML_coordinate_transform_XY4(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: dsrc(size(src))
    real(ark)                 :: alpha,alpha12,alpha13,alpha14,alpha23,alpha24,alpha34,s2a,s2b,beta312,beta412
    real(ark)                 :: cosa23,cosa24,cosa34,cosbeta,s(5)
    character(len=cl)         :: txt
    !
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"(/'ML_coordinate_transform_XY4/start')") 
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
    txt = 'ML_coord_XY4'
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
       endif
       !
    case('R-ALPHA')
       !
       if (direct) then
          !
          dst(1)=dsrc(1)
          dst(2)=dsrc(2)
          dst(3)=dsrc(3)
          dst(4)=dsrc(4)
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
            !alpha = molec%local_eq(10)
            !
          else
            !
            alpha34 = ML_XY4_calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
            !
            !alpha   = calc_alpha34(molec%local_eq(5),molec%local_eq(6),molec%local_eq(8),molec%local_eq(7),molec%local_eq(9))
            !
          endif
          !
          dst(5) = alpha12 - molec%local_eq(5)
          dst(6) = alpha14 - molec%local_eq(8)
          dst(7) = alpha23 - molec%local_eq(7)
          dst(8) = alpha13 - molec%local_eq(6)
          dst(9) = alpha24 - molec%local_eq(9)
          !
       else
          !
          dst(1:4)=dsrc(1:4) + molec%local_eq(1:4)
          !
          !if (size(dst)==10) then 
          !  !
          !  alpha   = calc_alpha34(molec%local_eq(5),molec%local_eq(6),molec%local_eq(8),molec%local_eq(7),molec%local_eq(9))
          !  !
          !endif
          !
          alpha12 = src(5)+molec%local_eq(5)
          alpha14 = src(6)+molec%local_eq(8)
          alpha23 = src(7)+molec%local_eq(7)
          alpha13 = src(8)+molec%local_eq(6)
          alpha24 = src(9)+molec%local_eq(9)
          !
          !alpha   = calc_alpha34(molec%local_eq(5),molec%local_eq(6),molec%local_eq(8),molec%local_eq(7),molec%local_eq(9))
          !
          !alpha23 = calc_alpha34(alpha14,alpha24,alpha34,alpha12,alpha13) ! alpha23

          alpha34 = ML_XY4_calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
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
       endif
       !
    case('SYM-S')
       !
       if (direct) then
          !
          alpha=dsrc(7)
          dsrc(7)=dsrc(8)
          dsrc(8)=alpha

          alpha=find_alpha34(src(5),src(6),src(7),src(8),src(9))

          alpha=alpha-molec%local_eq(5)
          dst(1)=(dsrc(1)+dsrc(2)+dsrc(3)+dsrc(4))*0.5_ark
          dst(2)=(2.0_ark*dsrc(5)-dsrc(6)-dsrc(7)-dsrc(8)-dsrc(9)+2.0_ark*alpha)/sqrt(12.0_ark)
          dst(3)=(dsrc(6)-dsrc(7)-dsrc(8)+dsrc(9))*0.5_ark
          dst(4)=(dsrc(1)-dsrc(2)+dsrc(3)-dsrc(4))*0.5_ark
          dst(5)=(dsrc(1)-dsrc(2)-dsrc(3)+dsrc(4))*0.5_ark
          dst(6)=(dsrc(1)+dsrc(2)-dsrc(3)-dsrc(4))*0.5_ark
          dst(7)=(dsrc(9)-dsrc(6))/sqrt(2.0_ark)
          dst(8)=(dsrc(8)-dsrc(7))/sqrt(2.0_ark)
          dst(9)=(alpha-dsrc(5))/sqrt(2.0_ark)
          !
       else
          !
          !  write(out,*) dsrc
          dst(1)=(dsrc(1)+dsrc(4)+dsrc(5)+dsrc(6))*0.5_ark+molec%local_eq(1)
          dst(2)=(dsrc(1)-dsrc(4)-dsrc(5)+dsrc(6))*0.5_ark+molec%local_eq(2)
          dst(3)=(dsrc(1)+dsrc(4)-dsrc(5)-dsrc(6))*0.5_ark+molec%local_eq(3)
          dst(4)=(dsrc(1)-dsrc(4)+dsrc(5)-dsrc(6))*0.5_ark+molec%local_eq(4)
          call find_alpha_for_XY4(dst(:),dsrc(:),dst(1),dst(2),dst(3),dst(4))
          !
          alpha=dst(7)
          dst(7)=dst(8)
          dst(8)=alpha

          !
       endif
       !
    case('SYM-S2')
       !
       if (direct) then
          !
          alpha12 = src(5)
          alpha13 = src(6)
          alpha23 = src(7)
          alpha14 = src(8)
          alpha24 = src(9)


          alpha34 = ML_XY4_calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)

          dst(1)=(dsrc(1)+dsrc(2)+dsrc(3)+dsrc(4))*0.5_ark
          dst(4)=(dsrc(1)-dsrc(2)+dsrc(3)-dsrc(4))*0.5_ark
          dst(5)=(dsrc(1)-dsrc(2)-dsrc(3)+dsrc(4))*0.5_ark
          dst(6)=(dsrc(1)+dsrc(2)-dsrc(3)-dsrc(4))*0.5_ark
          !
          dst(2)=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
          dst(3)=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
          !
          dst(7)=(alpha24-alpha13)/sqrt(2.0_ark)
          dst(8)=(alpha23-alpha14)/sqrt(2.0_ark)
          dst(9)=(alpha34-alpha12)/sqrt(2.0_ark)
          !
       else
          !
          !  write(out,*) dsrc
          !
          dst(1)=(dsrc(1)+dsrc(4)+dsrc(5)+dsrc(6))*0.5_ark+molec%local_eq(1)
          dst(2)=(dsrc(1)-dsrc(4)-dsrc(5)+dsrc(6))*0.5_ark+molec%local_eq(2)
          dst(3)=(dsrc(1)+dsrc(4)-dsrc(5)-dsrc(6))*0.5_ark+molec%local_eq(3)
          dst(4)=(dsrc(1)-dsrc(4)+dsrc(5)-dsrc(6))*0.5_ark+molec%local_eq(4)
          !
          dsrc(5) = src(2)
          dsrc(6) = src(3)
          dsrc(7) = src(7)
          dsrc(8) = src(8)
          dsrc(9) = src(9)
          !
          call from_sym2alphaII(dsrc(5:9),dst(5:9),alpha34)
          !
          !dst(5) = 2.0_ark/sqrt(12.0_ark)*src(2)-1.0_ark/sqrt(2.0_ark)*src(9)
          !dst(6) =-1.0_ark/sqrt(12.0_ark)*src(2)+0.5_ark*src(3)-1.0_ark/sqrt(2.0_ark)*src(7)
          !dst(7) =-1.0_ark/sqrt(12.0_ark)*src(2)-0.5_ark*src(3)+1.0_ark/sqrt(2.0_ark)*src(8)
          !dst(8) =-1.0_ark/sqrt(12.0_ark)*src(2)-0.5_ark*src(3)-1.0_ark/sqrt(2.0_ark)*src(8)
          !dst(9) =-1.0_ark/sqrt(12.0_ark)*src(2)+0.5_ark*src(3)+1.0_ark/sqrt(2.0_ark)*src(7)
          !
          !dst(5:9) = dst(5:9) + molec%local_eq(5:9)
          !
          !
       endif
       !
    case('R-SYM')
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
          dst(5)=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
          dst(6)=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
          dst(7)=(alpha24-alpha13)/sqrt(2.0_ark)
          dst(8)=(alpha23-alpha14)/sqrt(2.0_ark)
          dst(9)=(alpha34-alpha12)/sqrt(2.0_ark)
          !
       else
          !
          !  write(out,*) dsrc
          !
          dst(1:4)=dsrc(1:4) + molec%local_eq(1:4)
          !
          call from_sym2alphaII(dsrc(5:9),dst(5:9),alpha34)
          !
          alpha12 = dst(5)
          alpha13 = dst(6)
          alpha14 = dst(7)
          alpha23 = dst(8)
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
            !cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
            !beta312 = aacos(cosbeta,txt)
            !
            !cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
            !beta412 = aacos(cosbeta,txt)
            !
            !cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
            !alpha34 = aacos(cosa34,txt)
            !
            dst(10) = alpha34
            !
            !dst(10) = calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
            !
          endif 
          !
          !
       endif


    case('R-SYM-II')
       !

       !
       if (direct) then
          !
          alpha12 = src(5)
          alpha13 = src(6)
          alpha23 = src(7)
          alpha14 = src(8)
          alpha24 = src(9)

           
       !
       !minor   = reshape((/1.0_ark,cos(alpha12),cos(alpha13),&
       !                   cos(alpha12),1.0_ark,cos(alpha23),
       !                   cos(alpha14),cos(alpha24),cos(alpha34) /), 3,3)


          alpha34 = ML_XY4_calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)

          dst(1)=dsrc(1)
          dst(2)=dsrc(2)
          dst(3)=dsrc(3)
          dst(4)=dsrc(4)
          s2a=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
          s2b=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
          dst(7)=(alpha24-alpha13)/sqrt(2.0_ark)
          dst(8)=(alpha23-alpha14)/sqrt(2.0_ark)
          dst(9)=(alpha34-alpha12)/sqrt(2.0_ark)
          !
          dst(5) = sqrt(0.5_ark)*(s2a+s2b)
          dst(6) = sqrt(0.5_ark)*(s2a-s2b)
          !
       else
          !
          !  write(out,*) dsrc
          !
          dst(1:4)=dsrc(1:4) + molec%local_eq(1:4)
          !
          s2a = sqrt(0.5_ark)*(dsrc(5)+dsrc(6))
          s2b = sqrt(0.5_ark)*(dsrc(5)-dsrc(6))
          !
          dsrc(5)=s2a
          dsrc(6)=s2b
          call from_sym2alphaII(dsrc(5:9),dst(5:9),alpha34)
          !
          !
       endif

       !
    case('R-2-SYM')
       !
       if (direct) then
          !
          alpha12 = src(5)
          alpha13 = src(6)
          alpha14 = src(7)
          beta312 = src(8)
          beta412 = src(9)
          !
          cosa23 = cos(alpha12)*cos(alpha13)+cos(beta312)*sin(alpha12)*sin(alpha13)
          cosa24 = cos(alpha12)*cos(alpha14)+cos(beta412)*sin(alpha12)*sin(alpha14)
          cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
          !
          alpha23 = aacos(cosa23,txt)
          alpha24 = aacos(cosa24,txt)
          alpha34 = aacos(cosa34,txt)
          !
          alpha=find_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
          !
          alpha=ML_XY4_calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
          !
          if (abs(alpha-alpha34)>sqrt(small_)) then 
           !
           continue
           !
          endif
          !
          dst(1)=dsrc(1)
          dst(2)=dsrc(2)
          dst(3)=dsrc(3)
          dst(4)=dsrc(4)
          dst(5)=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
          dst(6)=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
          dst(7)=(alpha24-alpha13)/sqrt(2.0_ark)
          dst(8)=(alpha23-alpha14)/sqrt(2.0_ark)
          dst(9)=(alpha34-alpha12)/sqrt(2.0_ark)
          !
       else
          !
          !  write(out,*) dsrc
          !
          dst(1:4)=dsrc(1:4) + molec%local_eq(1:4)
          !
          !dsrc(2)=dsrc(5)
          !dsrc(3)=dsrc(6)
          !
          call from_sym2alphaII(dsrc(5:9),dst(5:9),alpha34)
          !
          alpha12 = dst(5)
          alpha13 = dst(6)
          alpha14 = dst(7)
          alpha23 = dst(8)
          alpha24 = dst(9)
          !
          cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
          beta312 = aacos(cosbeta,txt)
          !
          cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
          beta412 = aacos(cosbeta,txt)
          !
          dst(8) = beta312
          dst(9) = beta412
          !
       endif
       !
    case('R-S')
       !
       if (direct) then
          !
          alpha=find_alpha34(src(5),src(6),src(7),src(8),src(9))
          alpha=alpha-molec%local_eq(5)
          dst(1)=dsrc(1)
          dst(2)=dsrc(2)
          dst(3)=dsrc(3)
          dst(4)=dsrc(4)
          dst(5)=(2.0_ark*dsrc(5)-dsrc(6)-dsrc(7)-dsrc(8)-dsrc(9)+2.0_ark*alpha)/sqrt(12.0_ark)
          dst(6)=(dsrc(6)-dsrc(7)-dsrc(8)+dsrc(9))*0.5_ark
          dst(7)=(dsrc(9)-dsrc(6))/sqrt(2.0_ark)
          dst(8)=(dsrc(8)-dsrc(7))/sqrt(2.0_ark)
          dst(9)=(alpha-dsrc(5))/sqrt(2.0_ark)
          !
       else
          !
          !  write(out,*) dsrc
          !
          dsrc(1:4)=dsrc(1:4) + molec%local_eq(1:4)
          !
          !dsrc(2)=dsrc(5)
          !dsrc(3)=dsrc(6)
          call find_alpha_for_XY4(dst(:),dsrc(:),dsrc(1),dsrc(2),dsrc(3),dsrc(4))
          !
          !
       endif


       !
    case('R-SYM-F-E')
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
          dst(5)=(alpha24-alpha13)/sqrt(2.0_ark)
          dst(6)=(alpha23-alpha14)/sqrt(2.0_ark)
          dst(7)=(alpha34-alpha12)/sqrt(2.0_ark)
          dst(8)=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
          dst(9)=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
          !
       else
          !
          !  write(out,*) dsrc
          !
          dst(1:4)=dsrc(1:4) + molec%local_eq(1:4)
          !
          s(1:2) = dsrc(8:9)
          s(3:5) = dsrc(5:7)
          !
          call from_sym2alphaII(s(1:5),dst(5:9),alpha34)
          !
          alpha12 = dst(5)
          alpha13 = dst(6)
          alpha14 = dst(7)
          alpha23 = dst(8)
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
            !cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
            !beta312 = aacos(cosbeta,txt)
            !
            !cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
            !beta412 = aacos(cosbeta,txt)
            !
            !cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
            !alpha34 = aacos(cosa34,txt)
            !
            dst(10) = alpha34
            !
            !dst(10) = calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
            !
          endif 
          !
          !
       endif

       !
    end select
    !
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_XY4/end')") 
    !
    contains 

! Roman
    !


    function calc_delta_XY4(r1,r2,r3,r4,s5,s6,s7,s8,s9,s10,dsrt) result (delta_2)

      real(ark),intent(in)  :: s10
      real(ark) delta_2
      real(ark),intent(out) :: s5,s6,s7,s8,s9
      real(ark),intent(in)  :: dsrt(:),r1,r2,r3,r4
      
    



        s5=s10-dsrt(9)
        s6=(2.0_ark*dsrc(3)-sqrt(12.0_ark)*dsrc(2)-2.0_ark*sqrt(2.0_ark)*dsrc(7)+2.0_ark*(s5+s10))*0.25_ark
        s7=(-2.0_ark*dsrc(3)-sqrt(12.0_ark)*dsrc(2)-2.0_ark*sqrt(2.0_ark)*dsrc(8)+2.0_ark*(s5+s10))*0.25_ark
        s8=(-2.0_ark*dsrc(3)-sqrt(12.0_ark)*dsrc(2)+2.0_ark*sqrt(2.0_ark)*dsrc(8)+2.0_ark*(s5+s10))*0.25_ark
        s9=(2.0_ark*dsrc(3)-sqrt(12.0_ark)*dsrc(2)+2.0_ark*sqrt(2.0_ark)*dsrc(7)+2.0_ark*(s5+s10))*0.25_ark
        
        delta_2=find_alpha34(s5,s6,s7,s8,s9)-s10
       !
       !
    end function calc_delta_XY4
    !
    !
    !
    subroutine find_alpha_for_XY4(dst,dsrc,r1,r2,r3,r4)
    ! obtaining 6 angles from 5 internal coordinates+4 lenghts
                !5 alpha12
                !6 alpha13
                !7 alpha14
                !8 alpha23
                !9 alpha24
                !10 alpha34
                !dst(4)=(dsrt(1)-dsrt(4)+dsrt(5)-dsrt(6))*0.5_ark
        real(ark),intent(in)  :: dsrc(:),r1,r2,r3,r4
        real(ark),intent(out) :: dst(:)
        real(ark) :: eps,s10_old,f
        real(ark) :: rjacob,dx,s10

        real(ark) :: stadev_old,stability,stadev,ssq,stadev_best,fac_sign
    !
        integer(ik) :: iter,itmax
             
    !
    rjacob = 0 
    iter = 0
    stadev_old = 1.e10
    stability =  1.e10
    stadev    =  1.e10
    !
    dst(5:9) = 0
    !
    !
    dst(1)=r1
    dst(2)=r2
    dst(3)=r3
    dst(4)=r4
    !
    stadev_best = sqrt(small_)*10.0_ark
    itmax = 500
    !
    ! Initial value for alpha10
    !
! Roman CH4
    s10 = 1.9106332362490185563277142050315_ark
! Roman end
    !
    !s6 = alpha1*sqrt(3.0_ark)
    !
    !fac_sign = sign(1._ark,sindelta)
    !
!    write(out,*) dsrc
    outer_loop: & 
    do while( iter<itmax .and. stadev>stadev_best )   
       !
       iter = iter + 1
       ssq=0
       !
       ! Caclulate the function 
       !
       f = calc_delta_XY4(dst(1),dst(2),dst(3),dst(4),dst(5),dst(6),dst(7),dst(8),dst(9),s10,dsrc)
       !
       eps = f
       !
       ssq=abs(eps)
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       rjacob  = ( calc_delta_XY4(dst(1),dst(2),dst(3),dst(4),dst(5),dst(6),dst(7),dst(8),dst(9),s10+1e-7,dsrc)    &
                  -calc_delta_XY4(dst(1),dst(2),dst(3),dst(4),dst(5),dst(6),dst(7),dst(8),dst(9),s10-1e-7,dsrc) )/1e-7*0.5_ark
       !
       !
       if (itmax>=0) then
         !
         dx = eps/rjacob
         !
         stadev=sqrt(ssq)
         !
         s10_old=s10
!         write(out,*) s10
         !   
         ! Update the pot. parameters to the new values 
         !
         s10 = s10 - dx 
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

    end subroutine find_alpha_for_XY4
    !


    !
    !
  end function ML_coordinate_transform_XY4
  !

  subroutine from_sym2alpha(s,local)
        !
        ! obtaining 5 angles from 5 internal sym. coordinates
        !
        real(ark),intent(in)  :: s(5)
        real(ark),intent(out) :: local(5)
        !
        real(ark) :: eps,alpha34_old
        real(ark) :: rjacob,dx,stability

        real(ark) :: stadev_old,stadev,ssq,stadev_best,&
                     alpha12,alpha13,alpha14,alpha23,alpha24,alpha34,h
        !
        integer(ik) :: iter,itmax
             
 
        alpha34 = molec%local_eq(5)
        !
        !
        iter = 0
        stadev_old = 1.e10
        stadev    =  1.e10
        !
        !
        stadev_best = sqrt(small_)*1.0_ark
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
           eps = calc_alpha34fromsym(s,alpha34)
           !
           h = 1.e-3*abs(alpha34) ; if (h<1e-12) h = 1e-7
           !
           rjacob  = ( calc_alpha34fromsym(s,alpha34+h)-calc_alpha34fromsym(s,alpha34-h) )/h*0.5_ark
         !
         !
         dx = eps/rjacob
         !
         stadev=sqrt(eps**2)
         !
         alpha34_old=alpha34
         !   
         ! Update the pot. parameters to the new values 
         !
         alpha34 = alpha34 - dx 
         !      
         stability=abs( (stadev-stadev_old)/stadev )
         stadev_old=stadev
           !
           !
           !
           !ssq = sum( ( st(:)-s(:) )**2 )
           !
           !stadev=sqrt(ssq)
           !
           !stadev_old=stadev
         !
       enddo  outer_loop ! --- iter
       !
       if (iter==itmax) then
          write(out,"('from_sym2alpha: could not find solution after ',i8,' iterations')") iter
          stop 'from_sym2alpha: could not find solution'
       endif 
       !
       eps = calc_alpha34fromsym(s,alpha34)
       !
       !alpha34_old = calc_alpha34(alpha12,alpha13,alpha23,alpha14,alpha24)
       !
       local(1) = alpha12
       local(2) = alpha13
       local(3) = alpha23
       local(4) = alpha14
       local(5) = alpha24
       ! local(6) = alpha34
       !
       contains 


    function calc_alpha34fromsym(s,alpha34) result (ssq)

      real(ark),intent(in)  :: s(5),alpha34
      real(ark)             :: s2a,s2b,s4x,s4y,s4z,alpha34_t,ssq
        !
        S2a = s(1)
        S2b = s(2)
        S4x = s(3)
        S4y = s(4)
        S4z = s(5)
        !
        alpha12 = alpha34-sqrt(2.0_ark)*S4z
        alpha13 = -0.5_ark*sqrt(3.0_ark)*S2a+alpha34-0.5_ark*sqrt(2.0_ark)*S4z+0.5_ark*S2b-0.5_ark*sqrt(2.0_ark)*S4x
        alpha23 = -0.5_ark*sqrt(3.0_ark)*S2a+0.5_ark*sqrt(2.0_ark)*S4y+alpha34-0.5_ark*sqrt(2.0_ark)*S4z-0.5_ark*S2b
        alpha14 = -0.5_ark*sqrt(3.0_ark)*S2a-0.5_ark*sqrt(2.0_ark)*S4y+alpha34-0.5_ark*sqrt(2.0_ark)*S4z-0.5_ark*S2b
        alpha24 = -0.5_ark*sqrt(3.0_ark)*S2a+alpha34-0.5_ark*sqrt(2.0_ark)*S4z+0.5_ark*S2b+0.5_ark*sqrt(2.0_ark)*S4x
        !
        !alpha12 =  0.5_ark*(sqrt(2.0_ark)*alpha34-2.0_ark*S4z)*sqrt(2.0_ark)
        !alpha13 = -0.5_ark*sqrt(3.0_ark)*S2a+alpha34-0.5_ark*sqrt(2.0_ark)*S4z+0.5_ark*S2b-0.5_ark*sqrt(2.0_ark)*S4x 
        !alpha23 = -0.5_ark*(3.0_ark)*S2a+0.5_ark*sqrt(2.0_ark)*S4y+alpha34-0.5_ark*sqrt(2.0_ark)*S4z-0.5_ark*S2b 
        !alpha14 = -0.5_ark*(3.0_ark)*S2a-0.5_ark*sqrt(2.0_ark)*S4y+alpha34-0.5_ark*sqrt(2.0_ark)*S4z-0.5_ark*S2b 
        !alpha24 = -0.5_ark*(3.0_ark)*S2a+alpha34-0.5_ark*sqrt(2.0_ark)*S4z+0.5_ark*S2b+0.5_ark*sqrt(2.0_ark)*S4x 
        !
           ssq  = 1.0_ark-cos(alpha34)**2-      &
                    cos(alpha23)**2+            &
                    2.0_ark*cos(alpha23)*cos(alpha24)*cos(alpha34)-   &
                    cos(alpha24)**2-   &
                    cos(alpha12)**2+   &
                    cos(alpha12)**2*cos(alpha34)**2+  &
                    2.0_ark*cos(alpha12)*cos(alpha23)*cos(alpha13)-  &
                    2.0_ark*cos(alpha12)*cos(alpha23)*cos(alpha14)*cos(alpha34)-  &
                    2.0_ark*cos(alpha12)*cos(alpha24)*cos(alpha13)*cos(alpha34)+  &
                    2.0_ark*cos(alpha12)*cos(alpha24)*cos(alpha14)-  &
                    cos(alpha13)**2+  &
                    2.0_ark*cos(alpha13)*cos(alpha14)*cos(alpha34)+  &
                    cos(alpha13)**2*cos(alpha24)**2-  &
                    2.0_ark*cos(alpha13)*cos(alpha24)*cos(alpha14)*cos(alpha23)-  &
                    cos(alpha14)**2+cos(alpha14)**2*cos(alpha23)**2  

        !
        !alpha34_t = calc_alpha34(alpha12,alpha13,alpha23,alpha14,alpha24)
        !
        !
    end function calc_alpha34fromsym

       !

  end subroutine from_sym2alpha
  !
  !
  subroutine from_sym2alphaII(s,local,alpha34)
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
       local(:) = molec%local_eq(5)
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
          write(out,"('from_sym2alphaII: could not find solution after ',i8,' iterations')") iter
          stop 'from_sym2alphaII: could not find solution'
       endif 

       !
       contains 


    subroutine calc_sym_from_alpha(src,dst,alpha34)

      real(ark),intent(in)  :: src(5)
      real(ark),intent(out) :: dst(5),alpha34
      real(ark)             :: alpha12,alpha13,alpha23,alpha14,alpha24,cosbeta,beta312,beta412,cosa34
      character(len=cl)     :: txt
       !
       txt = 'calc_sym_from_alpha'
       !
       alpha12 = src(1)
       alpha13 = src(2)
       alpha14 = src(3)
       alpha23 = src(4)
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
       dst(1)=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
       dst(2)=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
       dst(3)=(alpha24-alpha13)/sqrt(2.0_ark)
       dst(4)=(alpha23-alpha14)/sqrt(2.0_ark)
       dst(5)=(alpha34-alpha12)/sqrt(2.0_ark)
       !
      end subroutine calc_sym_from_alpha
      !
  end subroutine from_sym2alphaII
  !
  !
  function find_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24) result (s10_f)
      !
      real(ark),intent(in)  ::alpha12,alpha13,alpha14,alpha23,alpha24
      real(ark) s10_f
      real(ark) xtmp(4),ytmp(4),ztmp(4)

        xtmp(1)=1.0_ark
        ytmp(1)=0
        ztmp(1)=0
        xtmp(2)=cos(alpha12)
        ytmp(2)=sin(alpha12)
        ztmp(2)=0


        xtmp(3)=cos(alpha13)
        ytmp(3)=(cos(alpha23)-xtmp(2)*xtmp(3))/ytmp(2)
        ztmp(3)=sqrt(1.0_ark-xtmp(3)**2-ytmp(3)**2)
        

        xtmp(4)=cos(alpha14)
        ytmp(4)=(cos(alpha24)-xtmp(2)*xtmp(4))/ytmp(2)
        ztmp(4)=-sqrt(1.0_ark-xtmp(4)**2-ytmp(4)**2)

        s10_f=acos((xtmp(3)*xtmp(4)+ytmp(3)*ytmp(4)+ztmp(3)*ztmp(4)))

  end function find_alpha34



  recursive subroutine ML_symmetry_transformation_XY4(ioper,natoms,src,dst)
    !
    integer(ik),intent(in)    :: ioper  ! group operation  
    integer(ik),intent(in)    :: natoms
    real(ark),intent(in)      :: src(natoms)
    real(ark),intent(out)     :: dst(natoms)
    !
    real(ark),dimension(size(src)) :: tmp
    !
    integer(ik)  :: tn(24,2)
    integer(ik) :: nsrc

    tn   = reshape(   &
                         (/0, 0, 2, 5, 3, 3,  6,21, 8, 2, 2, 2, 2,21,21,21,21,21,21,21,0,21,21,21,&
                           0, 0, 2, 5, 8, 4,  6,13, 8, 7, 8, 4,21, 3,11,10, 9, 2, 5, 7,0,12, 4, 6/),(/24,2/))
    !
    if (verbose>=5) write(out,"(/'ML_symmetry_transformation_XY4/start')") 
    !
    nsrc = size(src)
    !
    select case(trim(molec%symmetry))
    case default
       write (out,"('ML_symmetry_transformation_XY4: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_symmetry_transformation_XY4 - bad symm. type'
       !
    case('TD','TD(M)')
       !
       select case(trim(molec%coords_transform))
           !
           case default
           write (out,"('ML_symmetry_transformation_XY4: coord_transf ',a,' unknown')") trim(molec%coords_transform)
           stop 'ML_symmetry_transformation_XY4 - bad coord. type'
           !
       case('R-SYM','NORMAL')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst(:) = src(:)
             !
             return
             !
           case (2) ! (123)

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
             dst(4) = src(4)
             !
             !dst(5) = sum(sym%irr(3,2)%repres(1,1:2)*src(5:6))
             !dst(6) = sum(sym%irr(3,2)%repres(2,1:2)*src(5:6))
             !
             dst(5:6) = matmul((sym%irr(3,2)%repres(1:2,1:2)),src(5:6))
             !
             !dst(5) =-0.5_ark*src(5)              +0.5_ark*sqrt(3.0_ark)*src(6)
             !dst(6) =-0.5_ark*sqrt(3.0_ark)*src(5)-0.5_ark*src(6)
             !
             dst(7:9) = matmul((sym%irr(5,2)%repres(1:3,1:3)),src(7:9))
             !
             !dst(7) =-src(8)
             !dst(8) =-src(9)
             !dst(9) = src(7)
             !
             return
             !
           case (21) ! (14)*
             !
             dst(1) = src(4)
             dst(2) = src(2)
             dst(3) = src(3)
             dst(4) = src(1)
             !
             !dst(5) = sum(sym%irr(3,21)%repres(1,1:2)*src(5:6))
             !dst(6) = sum(sym%irr(3,21)%repres(2,1:2)*src(5:6))
             !
             dst(5:6) = matmul((sym%irr(3,21)%repres(1:2,1:2)),src(5:6))
             !
             !dst(5) =-0.5_ark*src(5)              +0.5_ark*sqrt(3.0_ark)*src(6)
             !dst(6) = 0.5_ark*sqrt(3.0_ark)*src(5)+0.5_ark*src(6)
             !
             !dst(7) =-src(9)
             !dst(8) = src(8)
             !dst(9) =-src(7)
             !
             dst(7:9) = matmul((sym%irr(5,21)%repres(1:3,1:3)),src(7:9))
             !
             return
             !
           end select 
           !
        case('SYM-S2')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst(:) = src(:)
             !
             return
             !
           case (2) ! (123)

             dst(1) = src(1)
             dst(2) =-0.5_ark*src(2)               +0.5_ark*sqrt(3.0_ark)*src(3)
             dst(3) =-0.5_ark*sqrt(3.0_ark)*src(2)  -0.5_ark*src(3)
             dst(4) =-src(6)
             dst(5) =-src(4)
             dst(6) =-src(5)
             dst(7) =-src(9)
             dst(8) =-src(7)
             dst(9) = src(8)
             !
             return
             !
           case (21) ! (14)*
             !
             dst(1) = src(1)
             dst(2) =-0.5_ark*src(2)               -0.5_ark*sqrt(3.0_ark)*src(3)
             dst(3) =-0.5_ark*sqrt(3.0_ark)*src(2)  +0.5_ark*src(3)
             dst(4) = src(5)
             dst(5) = src(4)
             dst(6) = src(6)
             dst(7) = src(8)
             dst(8) = src(7)
             dst(9) = src(9)
             !
             return
             !
           end select
           !
       case('R-SYM-F-E')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst(:) = src(:)
             !
             return
             !
           case (2) ! (123)

             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
             dst(4) = src(4)
             !
             !dst(5) = sum(sym%irr(3,2)%repres(1,1:2)*src(5:6))
             !dst(6) = sum(sym%irr(3,2)%repres(2,1:2)*src(5:6))
             !
             dst(8:9) = matmul((sym%irr(3,2)%repres(1:2,1:2)),src(8:9))
             !
             !dst(5) =-0.5_ark*src(5)              +0.5_ark*sqrt(3.0_ark)*src(6)
             !dst(6) =-0.5_ark*sqrt(3.0_ark)*src(5)-0.5_ark*src(6)
             !
             dst(5:7) = matmul((sym%irr(5,2)%repres(1:3,1:3)),src(5:7))
             !
             !dst(7) =-src(8)
             !dst(8) =-src(9)
             !dst(9) = src(7)
             !
             return
             !
           case (21) ! (14)*
             !
             dst(1) = src(4)
             dst(2) = src(2)
             dst(3) = src(3)
             dst(4) = src(1)
             !
             !dst(5) = sum(sym%irr(3,21)%repres(1,1:2)*src(5:6))
             !dst(6) = sum(sym%irr(3,21)%repres(2,1:2)*src(5:6))
             !
             dst(8:9) = matmul((sym%irr(3,21)%repres(1:2,1:2)),src(8:9))
             !
             !dst(5) =-0.5_ark*src(5)              +0.5_ark*sqrt(3.0_ark)*src(6)
             !dst(6) = 0.5_ark*sqrt(3.0_ark)*src(5)+0.5_ark*src(6)
             !
             !dst(7) =-src(9)
             !dst(8) = src(8)
             !dst(9) =-src(7)
             !
             dst(5:7) = matmul((sym%irr(5,21)%repres(1:3,1:3)),src(5:7))
             !
             return
             !
           end select 
           !
       end select 
       !
       if (all(tn(ioper,:)/=0)) then 
          call ML_symmetry_transformation_XY4(tn(ioper,1),natoms,src,tmp)
          call ML_symmetry_transformation_XY4(tn(ioper,2),natoms,tmp,dst)
       endif 
       !
    case('C2V','C2V(M)')
       !
       select case(trim(molec%coords_transform))
           !
           case default
           write (out,"('ML_symmetry_transformation_XY4: coord_transf ',a,' unknown')") trim(molec%coords_transform)
           stop 'ML_symmetry_transformation_XY4 - bad coord. type'
           !
       case('R-SYM')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst(:) = src(:)
             !
             return
             !
           case (2) ! (14)(23)

             dst(1) = src(4)
             dst(2) = src(3)
             dst(3) = src(2)
             dst(4) = src(1)

             dst(5) = src(5)
             dst(6) = src(6)
             !
             dst(7) =-src(7)
             dst(8) =-src(8)
             dst(9) = src(9)
             !
             return
             !
           case (4) ! (14)*
             !
             dst(1) = src(4)
             dst(2) = src(2)
             dst(3) = src(3)
             dst(4) = src(1)
             dst(5) = -0.5_ark*src(5)                +0.5_ark*sqrt(3.0_ark)*src(6)
             dst(6) =  0.5_ark*sqrt(3.0_ark)*src(5)  +0.5_ark*src(6)
             !
             dst(7) = -src(9)
             dst(8) = -src(8)
             dst(9) = -src(7)
             !
             return
             !
           case (3) ! (123)
             !
             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
             dst(4) = src(4)
             !
             dst(5) = -0.5_ark*src(5)                +0.5_ark*sqrt(3.0_ark)*src(6)
             dst(6) = -0.5_ark*sqrt(3.0_ark)*src(5)  -0.5_ark*src(6)
             !
             dst(7) = src(7)
             dst(8) = src(9)
             dst(9) = src(8)             
             !
           case (5) ! (23)*
             !
             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)
             dst(4) = src(4)
             dst(5) = src(5)
             dst(6) =-src(6)
             !
             dst(7) = src(7)
             dst(8) = src(9)
             dst(9) = src(8)             
             !
             return
             !
           end select
           !
       case('R-SYM-F-E')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst(:) = src(:)
             !
             return
             !
           case (2) ! (14)(23)

             dst(1) = src(4)
             dst(2) = src(3)
             dst(3) = src(2)
             dst(4) = src(1)

             dst(8) = src(8)
             dst(9) = src(9)
             !
             dst(5) =-src(5)
             dst(6) =-src(6)
             dst(7) = src(7)
             !
             return
             !
           case (4) ! (14)*
             !
             dst(1) = src(4)
             dst(2) = src(2)
             dst(3) = src(3)
             dst(4) = src(1)
             dst(8) = -0.5_ark*src(8)                +0.5_ark*sqrt(3.0_ark)*src(9)
             dst(9) =  0.5_ark*sqrt(3.0_ark)*src(8)  +0.5_ark*src(9)
             !
             dst(5) = -src(7)
             dst(6) = -src(6)
             dst(7) = -src(5)
             !
             return
             !
           case (3) ! (123)
             !
             dst(1) = src(3)
             dst(2) = src(1)
             dst(3) = src(2)
             dst(4) = src(4)
             !
             dst(8) = -0.5_ark*src(8)                +0.5_ark*sqrt(3.0_ark)*src(9)
             dst(9) = -0.5_ark*sqrt(3.0_ark)*src(8)  -0.5_ark*src(9)
             !
             dst(5) = src(5)
             dst(6) = src(7)
             dst(7) = src(6)             
             !
           case (5) ! (23)*
             !
             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)
             dst(4) = src(4)
             dst(7) = src(7)
             dst(8) =-src(8)
             !
             dst(5) = src(5)
             dst(6) = src(7)
             dst(7) = src(6)             
             !
             return
             !
           end select           
           !
        case('R-ALPHA')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst(:) = src(:)
             !
             return
             !
           case (2) ! (12)(34)

             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(4)
             dst(4) = src(3)
             !
             dst(5) = src(5)
             !
             dst(6) = src(7)
             dst(7) = src(6)
             dst(8) = src(9)
             dst(9) = src(8)
             !
             return
             !
           case (3) ! (12)*
             !
             dst(1) = src(2)
             dst(2) = src(1)
             dst(3) = src(3)
             dst(4) = src(4)
             !
             dst(5) = src(5)
             !
             dst(6) = src(9)
             dst(7) = src(8)
             dst(8) = src(7)
             dst(9) = src(6)
             !
           case (4) ! (34)*

             dst(1) = src(1)
             dst(2) = src(2)
             dst(3) = src(4)
             dst(4) = src(3)
             !
             dst(5) = src(5)
             !
             dst(6) = src(8)
             dst(7) = src(9)
             dst(8) = src(6)
             dst(9) = src(7)
             !
           end select
           !
        case('R-ALPHA-??')
           !
           select case(ioper)
           !
           case (1) ! identity 

             dst(:) = src(:)
             !
             return
             !
           case (2) ! (14)(23)

             dst(1) = src(4)
             dst(2) = src(3)
             dst(3) = src(2)
             dst(4) = src(1)
             !
             dst(5) = src(5)
             !
             dst(6) = src(7)
             dst(7) = src(6)
             dst(8) = src(9)
             dst(9) = src(8)
             !
             return
             !
           case (3) ! (14)*
             !
             dst(1) = src(4)
             dst(2) = src(2)
             dst(3) = src(3)
             dst(4) = src(1)
             !
             dst(5) = src(5)
             !
             dst(6) = src(9)
             dst(7) = src(8)
             dst(8) = src(7)
             dst(9) = src(6)
             !
           case (4) ! (23)*
             !
             dst(1) = src(1)
             dst(2) = src(3)
             dst(3) = src(2)
             dst(4) = src(4)
             !
             dst(5) = src(5)
             !
             dst(6) = src(8)
             dst(7) = src(9)
             dst(8) = src(6)
             dst(9) = src(7)
             !
           end select
           !
       end select  
       !
    end select
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY4/end')") 
    !
  end subroutine ML_symmetry_transformation_XY4


  ! Here we find the rotational symmetry from J,K
  !
  subroutine ML_rotsymmetry_XY4(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_XY4/start')") 
    !
    !
    select case(trim(molec%symmetry))
    case default
       write (out,"('ML_rotsymmetry_XY4: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_rotsymmetry_XY4 - bad symm. type'
       !
    case('TD','TD(M)')
       !
       write(out,"('The Td symmetry cannot be defined as k-based. Change rotsym to euler!')")
       stop 'Change rotsym to euler!'
       !
    case('C2V','C2V(M)')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma = 1 !1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma = 2 !3 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma = 3 !4 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma = 4 !2 !; return
       !
    case('CS','CS(M)')
       !
       gamma = 0 
       ideg = 1
       !
       if (mod(K+tau+2,2)==0) gamma = 1 !; return
       if (mod(K+tau+2,2)/=0) gamma = 2 !; return
       !
    end select
    !
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_XY4/end')") 
    !
    !
  end subroutine ML_rotsymmetry_XY4

end module mol_xy4

!
!  This unit defines all specific routines for a fouratomic molecule of ZXY2 type
!
module mol_zxy2
  use accuracy
  use moltype

  implicit none

  public ML_b0_ZXY2,ML_coordinate_transform_ZXY2,ML_b0_SOHF
  public ML_coordinate_transform_SOHF,ML_symmetry_transformation_ZXY2
  public ML_rotsymmetry_ZXY2
  public ML_MEP_zxy2_R_rho,ML_MEP_zxy2_rho_coeff
  !public MLpoten_zxy2_mlt,MLpoten_zxy2_mep_r_alpha_rho_powers,MLpotential_zxy2_mep_r_alpha_rho_powers
  !public MLdms2xyz_zxy2_symadap_powers,MLpoten_zxy2_andrey_coeff,MLpoten_h2cs_tz_damp1
  !public MLpoten_h2cs_damp,MLpoten_sohf,MLpoten_zxy2_andrey_01
  private
 
  integer(ik), parameter :: verbose     = 3                          ! Verbosity level


  contains


  !
  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  function ML_coordinate_transform_ZXY2(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: dsrc(size(src))
    real(ark)                 :: rho1,rho2,rho3,alpha1,alpha2,alpha3,cosdelta,cosalpha,delta,req(6)
    real(ark)                 :: theta1,theta2,theta,phi,theta3,costheta3,tx1,ty1,tx2,ty2,tx0,t1,x1,x2,cosphi
    real(ark)                 :: norm_2,tau,sindelta
    !
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_ZXY2/start')") 
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
    case('SYMM-C2V')
       !
       !
       if (direct) then
          ! 
          rho1 = (src(1) - molec%local_eq(1))/(src(1))
          rho2 = (src(2) - molec%local_eq(2))/(src(2))
          rho3 = (src(3) - molec%local_eq(3))/(src(3))
          !
          dst(1) = (rho2+rho3)/sqrt(2.0_ark)
          dst(2) =  rho1
          dst(3) = (dsrc(4)+dsrc(5))/sqrt(2.0_ark)
          dst(4) =  dsrc(6)
          dst(5) = (rho2-rho3)/sqrt(2.0_ark)
          dst(6) = (dsrc(4)-dsrc(5))/sqrt(2.0_ark)
          !
       else    ! not direct
          !
          dst(6) = dsrc(4)+molec%local_eq(6)
          !
          rho2 = (src(1)+src(5))/sqrt(2.0_ark)
          rho3 = (src(1)-src(5))/sqrt(2.0_ark)
          rho1 = src(2)
          !
          dst(1) =  molec%local_eq(1)/(1.0_ark-rho1)
          dst(2) =  molec%local_eq(2)/(1.0_ark-rho2)
          dst(3) =  molec%local_eq(3)/(1.0_ark-rho3)
          !
          dst(4) = (src(3)+src(6))/sqrt(2.0_ark)+molec%local_eq(4)
          dst(5) = (src(3)-src(6))/sqrt(2.0_ark)+molec%local_eq(5)
          !
      endif
          !
          !
    case('R-S-TAU')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          dst(4) = (dsrc(4)+dsrc(5))/sqrt(2.0_ark)
          dst(5) = (dsrc(4)-dsrc(5))/sqrt(2.0_ark)
          !
          dst(6) =  src(6)
          !
          !
       else    ! not direct
          !
          dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
          dst(6) = src(6)
          !
          dst(4) = (src(4)+src(5))/sqrt(2.0_ark)+molec%local_eq(4)
          dst(5) = (src(4)-src(5))/sqrt(2.0_ark)+molec%local_eq(5)
          !
      endif
       !
    case('R-THETA-TAU')
       !
       if (direct) then
         !
         if (size(src)==7) then
           !
           dst(1:3) = dsrc(1:3)
           dst(4:5) = dsrc(4:5)
           dst(6)   = dsrc(7)
           !
         else
           ! 
           dst(:) = dsrc(:)
           !
         endif 
          !
       else    ! not direct
          !
         if (size(dst)==7) then
           !
           dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
           dst(4:5) = dsrc(4:5)+molec%local_eq(4:5)
           dst(7)   = dsrc(6)+molec%local_eq(7)
           !
         else
           ! 
           dst(:) = dsrc(:)+molec%local_eq(:)
           !
         endif 
          !
      endif
       !
    case('R-THETA-TAU-MEP')
       !
       delta = src(6)
       !
       req(:) = ML_MEP_zxy2_R_rho(delta)
       req(6) = pi
       !
       !cosdelta = cos(delta) + 1.0_ark
       !
       !req(1) = molec%local_eq(1)+sum(molec%force(2:5)*cosdelta**molec%pot_ind(1,2:5))
       !req(2) = molec%local_eq(2)+sum(molec%force(7:10)*cosdelta**molec%pot_ind(2,7:10))
       !req(3) = req(2)
       !req(4) = molec%local_eq(4)+sum(molec%force(12:15)*cosdelta**molec%pot_ind(4,12:15))
       !req(5) = req(4)
       !req(6) = pi 
       !
       if (direct) then
          !
          dst(:) = src(:)-req(:)
          !
      else ! not direct
          !
          dst(:) = src(:)+req(:)
          !
      endif
      !
    case('R-THETA-DELTA1')
       !
       ! delta1 is the angle between r1 and [r2,r3]/(r2*r3*sin(alpha32))
       !
       if (direct) then
          ! 
          dst(1:5) = dsrc(1:5)
          !
          if (size(src)==7) then
              !
              alpha3 = src(6)
              !
              cosdelta = src(7)/sin(alpha3)
              !
              if ( abs(cosdelta)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLcoordinate_transform_func: cosdelta>1: ',f18.8)") cosdelta
                 write (out,"('Consider change difftype ')")
                 stop 'MLcoordinate_transform_func - bad cosdelta'
                 !
              elseif ( cosdelta>=1.0_ark) then 
                 !
                 dst(6) = 0.0_ark
                 !
              else 
                 dst(6) = acos(cosdelta)
                 !
              endif 
              !
          else
              !
              write (out,"('MLcoordinate_transform_func: wrong number of coordinates for R-THETA-DELTA: ',i6)") size(src)
              stop 'MLcoordinate_transform_func - wrong number of coords'
              !
          endif
          !
      else    ! not direct
          !
          if (size(dst)==7) then
             !
             dst(1:5) = dsrc(1:5)+molec%local_eq(1:5)
             !
             alpha1 = dsrc(4)+molec%local_eq(4)
             alpha2 = dsrc(5)+molec%local_eq(5)
             !
             delta = src(6)
             !
             cosalpha=0.5_ark*(-2._ark*cos(alpha2)*cos(alpha1)+2._ark*sqrt(cos(alpha2)**2*cos(alpha1)**2 &
             -2._ark*cos(delta)**2+cos(delta)**2*cos(alpha1)**2+cos(delta)**2*cos(alpha2)**2 &
             +cos(delta)**4+1._ark-cos(alpha1)**2-cos(alpha2)**2))/(cos(delta)**2-1._ark)                  
             !
             if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
                !
                write (out,"('MLcoordinate_transform_func: cosalpha>1: ',f18.8)") cosalpha
                write (out,"('Consider change difftype ')")
                stop 'MLcoordinate_transform_func - bad cosalpha'
                !
             elseif ( cosalpha>=1.0_ark) then 
                !
                alpha3 = 0.0_ark
                !
             else 
                alpha3 = acos(cosalpha)
                !
             endif 
             !
             dst(6) = alpha3
             !
             dst(7) = cos(delta)*sin(alpha3)
             !
         endif 
      endif 
      !
    case('R-THETA-DELTA')
       !
       if (direct) then
          ! 
          dst(1:3) = dsrc(1:3)
          !
          if (size(src)==7) then
              !
              !dst(4) = dsrc(4)
              !dst(5) = dsrc(5)
              !dst(6) = src(7)
              !
              !theta1 = src(4)
              !theta2 = src(5)
              !theta3 = src(6)
              !
              !cosphi = (-cos(theta2)*cos(theta1)+cos(theta3))/(sin(theta2)*sin(theta1))
              !
              !if ( abs(cosphi)>1.0_ark+sqrt(small_) ) then 
              !   !
              !   write (out,"('MLcoordinate_transform_func: cosphi>1: ',f18.8)") cosphi
              !   stop 'MLcoordinate_transform_func - bad cosphi'
              !   !
              !elseif ( cosphi>=1.0_ark) then 
              !   !
              !   phi = 0.0_ark
              !   !
              !elseif ( cosphi<=-1.0_ark) then 
              !   !
              !   phi = pi
              !   !
              !else 
              !   !
              !   phi = acos(cosphi)
              !   !
              !endif
              !
              !if (src(7)>0) phi = 2.0_ark*pi-phi
              !
              dst(4) = dsrc(4)
              dst(5) = dsrc(5)
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
          !
          if (size(dst)==7) then
             !
             dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
             !
             theta1 = dsrc(4)+molec%local_eq(4)
             theta2 = dsrc(5)+molec%local_eq(5)
             !
             delta = src(6)
             !
             !
             t1 = ( cos(theta2)-1.0_ark )*( cos(theta2)+2.0_ark*cos(delta)**2-1.0_ark )*( cos(theta1)-1.0_ark )*( cos(theta1)-1.0_ark+2.0_ark*cos(delta)**2 )
             !
             if ( t1<-sqrt(small_) ) then
                write (out,"('MLcoordinate_transform_func: sqrt(-1) for t1= ',f18.8)") t1
                stop 'MLcoordinate_transform_func - bad sqrt(-1)'
             elseif ( t1<sqrt(small_) ) then 
                t1 = 0
             endif
             !
             x1 = -(-1.0_ark+cos(delta)**2+cos(theta2)-cos(theta2)*cos(delta)**2+cos(theta1)-cos(theta1)*cos(delta)**2-cos(theta2)*cos(theta1)+cos(theta1)*cos(theta2)*cos(delta)**2+sqrt( t1 ) )/(cos(delta)**2*sin(theta1)*sin(theta2))
             !
             t1 = ( cos(theta1)-1.0_ark )*( cos(theta1)-1.0_ark+2.0_ark*cos(delta)**2 )*( 2.0_ark*cos(theta2)*cos(delta)**2-2.0_ark*cos(delta)**2+1.0_ark-2.0_ark*cos(theta2)+cos(theta2)**2 ) 
             !
             if ( t1<-sqrt(small_) ) then
                write (out,"('MLcoordinate_transform_func: sqrt(-1) for t1= ',f18.8)") t1
                stop 'MLcoordinate_transform_func - bad sqrt(-1)'
             elseif ( t1<sqrt(small_) ) then 
                t1 = 0
             endif
             !
             x2 = -(-1.0_ark+cos(delta)**2+cos(theta2)-cos(theta2)*cos(delta)**2+cos(theta1)-cos(theta1)*cos(delta)**2-cos(theta2)*cos(theta1)+cos(theta1)*cos(theta2)*cos(delta)**2-sqrt( t1 ) )/(cos(delta)**2*sin(theta1)*sin(theta2))
             !
             cosphi = x1
             !
             costheta3 =sin(theta2)*cosphi*sin(theta1)+cos(theta2)*cos(theta1)
 
             !
             !if ( abs(cosphi)>1.0_ark+sqrt(small_) ) then 
             !   !
             !   write (out,"('MLcoordinate_transform_func: costheta3>1: ',f18.8)") costheta3
             !   stop 'MLcoordinate_transform_func - bad costheta3'
             !   !
             !elseif ( cosphi>=1.0_ark) then 
             !   !
             !   phi = 0.0_ark
             !   !
             !else 
             !   phi = acos(cosphi)
             !   !
             !endif
             !
             if ( abs(costheta3)>1.0_ark+sqrt(small_) ) then 
                !
                write (out,"('MLcoordinate_transform_func: costheta3>1: ',f18.8)") costheta3
                stop 'MLcoordinate_transform_func - bad costheta3'
                !
             elseif ( costheta3>=1.0_ark) then 
                !
                theta3 = 0.0_ark
                !
             else 
                theta3 = acos(costheta3)
                !
             endif
             !
             !tau = sin(theta1)*sin(theta2)*sin(phi)
             !
             !norm_2 = 3.0_ark-cos(theta3)**2-cos(theta2)**2-cos(theta1)**2+2.0_ark*cos(theta3)*cos(theta2)-2.0_ark*cos(theta2)+&
             !          2.0_ark*cos(theta2)*cos(theta3)-2.0_ark*cos(theta2)+2.0_ark*cos(theta2)*cos(theta2)-2.0_ark*cos(theta3)
             !!
             !if (norm_2<small_.and.abs(tau)<small_) then 
             !  !
             !  sindelta = 1.0_ark
             !  !
             !elseif(norm_2<small_) then 
             !  !
             !  write (out,"('MLcoordinate_transform_func: norm2 = ',f18.8', delta = infty!')") norm_2
             !  stop 'MLcoordinate_transform_func - bad norm2'
             !  !
             !else
             !  !
             !  sindelta = tau/sqrt(norm_2)
             !  !
             !endif
             !!
             !if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
             !   !
             !   write (out,"('MLcoordinate_transform_func: sindelta>1: ',f18.8)") sindelta
             !   stop 'MLcoordinate_transform_func - bad sindelta'
             !   !
             !elseif ( sindelta>=1.0_ark) then 
             !   !
             !   delta = 0.0_ark
             !   !
             !else 
             !   !
             !   delta = asin(sindelta)
             !   !
             !endif 
             !
             dst(4) = theta1
             dst(5) = theta2
             dst(6) = theta3
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
    end select
    !
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_ZXY2/end')") 
    !
    !
  end function ML_coordinate_transform_ZXY2





  !
  !
  ! Here we define structural parameters a0 for rigid XY twoatomic molecule.
  !
  subroutine ML_b0_ZXY2(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)

     integer(ik),intent(in)  :: Npoints,Natoms

     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in),optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional :: rho_ref
     real(ark),   intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders

     real(ark)               :: a0(molec%Natoms,3),CM_shift,tau,cosalpha_2,costau,delta,cosrho,req(1:3),alphaeq(1:2),tau_
     real(ark)               :: xieq(6),rho,theta,sint_2,theta12,beta
     integer(ik)             :: Nbonds,i,n


      if (verbose>=4) write(out,"('ML_b0_ZXY2/start')") 

      Nbonds = size(molec%req)

      if (Nbonds/=3) then
        write(out,"('ML_b0_ZXY2: Nbonds must be 3 in this routine, not  ',i9)") Nbonds
        stop 'ML_b0_ZXY2: wrong Nbonds '
      endif 

      if (molec%Natoms/=4) then
        write(out,"('ML_b0_ZXY2: Natoms must be 4 in this routine, not  ',i9)") molec%Natoms
        stop 'ML_b0_ZXY2: wrong Natoms '
      endif 

      select case(molec%dihedtype(1))
      case default
         write (out,"('ML_b0_ZXY2: dihedral type ',i4,' is not working here')") molec%dihedtype(1)
         stop 'ML_b0_ZXY2 - bad dihedral type'
         !
      case( 1) 
         !
         tau = pi 
         !
      case(-2) 
         !
         tau = molec%taueq(1)
         !
      case(2) 
         !
         tau = molec%taueq(1)
         !
      end select 
      !
      if (trim(molec%coords_transform)=='R-THETA-TAU-MEP') then 
         !
         !if (trim(molec%potentype)/='POTEN_ZXY2_MEP_R_ALPHA_RHO_POWERS') then
         ! write (out,"('ML_b0_ZXY2: potential function',a,' cannot be used with coord type ',a)") trim(molec%potentype),trim(molec%coords_transform)
         ! stop 'ML_b0_ZXY2 - bad PES type'
         !endif 
         !
         !
         rho = tau
         !
         xieq(:) = ML_MEP_zxy2_R_rho(rho)
         !
         req(1:3)     = xieq(1:3)
         alphaeq(1:2) = xieq(1:2)
         !
      endif 
      !
      a0(1,1) = 0.0_ark
      a0(1,2) = 0.0_ark
      a0(1,3) = 0.0_ark
      !
      a0(2,1) = 0.0_ark
      a0(2,2) = 0.0_ark
      a0(2,3) = molec%req(1)
      !
      a0(3,1) = molec%req(2)*sin(molec%alphaeq(1))*cos(tau*0.5_ark)
      a0(3,2) =-molec%req(2)*sin(molec%alphaeq(1))*sin(tau*0.5_ark)
      a0(3,3) = molec%req(2)*cos(molec%alphaeq(1))
      !
      a0(4,1) = molec%req(2)*sin(molec%alphaeq(1))*cos(tau*0.5_ark)
      a0(4,2) = molec%req(2)*sin(molec%alphaeq(1))*sin(tau*0.5_ark)
      a0(4,3) = molec%req(2)*cos(molec%alphaeq(1))
      !
      do n = 1,3 
        CM_shift = sum(a0(:,n)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
        a0(:,n) = a0(:,n) - CM_shift
      enddo 
      !
      b0(:,:,0) = a0(:,:)
      !
      !
      if (Npoints/=0) then
         ! 
         if (.not.present(rho_borders).or..not.present(rho_ref)) then  
            write(out,"('ML_b0_ZXY2: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_ZXY2: rho_borders or rho_ref not specified '
         endif
         !
         rho_ref = 0 ! pi
         !
         do i = 0,npoints
            !
            theta = molec%alphaeq(1)
            !
            tau = rho_i(i)
            !
            select case(molec%dihedtype(1))
               !
            case( -1,1) 
               !
               beta = rho_i(i)+pi*0.5_ark
               !
               sint_2 = sin(molec%alphaeq(1)*0.5_ark)*sin(beta)
               !
               if ( abs(sint_2)>1.0_ark+sqrt(small_) ) then 
                  !
                  write (out,"('ML_b0_ZXY2: sint_2>1: ',f18.8)") sint_2
                  stop 'ML_b0_ZXY2 - bad sint_2'
                  !
               elseif ( sint_2>=1.0_ark) then 
                  theta = pi
               else 
                  theta = 2.0_ark*asin(sint_2)
               endif
               !
               theta12 = 2.0_ark*(pi-molec%alphaeq(1))
               !
               sint_2 = sin(molec%alphaeq(1))*sin(beta)
               !
               if ( abs(sint_2)>1.0_ark+sqrt(small_) ) then 
                  !
                  write (out,"('ML_b0_ZXY2: sint_2>1: ',f18.8)") sint_2
                  stop 'ML_b0_ZXY2 - bad sint_2'
                  !
               elseif ( sint_2>=1.0_ark) then 
                  theta12 = pi
               else 
                  theta12 = 2.0_ark*asin(sint_2)
               endif
               !
               costau = (-cos(theta)*cos(theta)+cos(theta12))/(sin(theta)*sin(theta))
               !
               if ( abs(costau)>1.0_ark+sqrt(small_) ) then 
                  !
                  write (out,"('ML_b0_ZXY2: costau>1: ',f18.8)") costau
                  stop 'ML_b0_ZXY2 - bad costau'
                  !
               elseif ( costau>=1.0_ark) then 
                  tau = 0.0_ark
               elseif ( costau<=-1.0_ark) then 
                  tau = pi
               else 
                  tau = acos(costau)
               endif
               !
               if (beta>pi*0.5_ark) tau = 2.0_ark*pi-tau
               !
            case(-2,2) 
               !
               tau = rho_i(i)+pi
               !
            end select
            !
            if (trim(molec%coords_transform)=='R-THETA-TAU-MEP') then 
               !
               rho = tau
               !
               xieq(:) = ML_MEP_zxy2_R_rho(rho)
               !
               req(1:3)     = xieq(1:3)
               alphaeq(1:2) = xieq(1:2)
               !
               !cosrho = cos(tau) + 1.0_ark
               !
               !req(1)     = molec%req(1)+sum(molec%force(2:5 )*cosrho**molec%pot_ind(1,2:5 ))
               !req(2)     = molec%req(2)+sum(molec%force(7:10)*cosrho**molec%pot_ind(2,7:10))
               !req(3)     = req(2)
               !alphaeq(1) = molec%alphaeq(1)+sum(molec%force(12:15)*cosrho**molec%pot_ind(4,12:15))
               !alphaeq(2) = alphaeq(1)
               !
               !
            endif 
            !
            b0(1,1,i) = 0.0_ark
            b0(1,2,i) = 0.0_ark
            b0(1,3,i) = 0.0_ark
            !
            b0(2,1,i) = 0.0_ark
            b0(2,2,i) = 0.0_ark
            b0(2,3,i) = molec%req(1)
            !
            b0(3,1,i) = molec%req(2)*sin(theta)*cos(tau*0.5_ark)
            b0(3,2,i) =-molec%req(2)*sin(theta)*sin(tau*0.5_ark)
            b0(3,3,i) = molec%req(2)*cos(theta)
            !
            b0(4,1,i) = molec%req(2)*sin(theta)*cos(tau*0.5_ark)
            b0(4,2,i) = molec%req(2)*sin(theta)*sin(tau*0.5_ark)
            b0(4,3,i) = molec%req(2)*cos(theta)
            !
            call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i))
            !
            ! Find center of mass
            !
            do n = 1,3 
              CM_shift = sum(b0(:,n,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
              b0(:,n,i) = b0(:,n,i) - CM_shift
            enddo 
            !
            !do n = 1,molec%Natoms
            !   b0(n,:,i) = matmul(transpose(transform),b0(n,:,i))
            !enddo
            !
            if (verbose>=3) then 
              write(out,"(i6)") molec%natoms
              !
              write(out,"(/'C',3x,3f14.8)") b0(1,:,i)
              write(out,"( 'O',3x,3f14.8)") b0(2,:,i)
              write(out,"( 'H',3x,3f14.8)") b0(3,:,i)
              write(out,"( 'H',3x,3f14.8)") b0(4,:,i)
              !
            endif
            !
            if (verbose>=4) then 
              write(out,"('b0',i4,12f12.8)") i,b0(:,:,i)
            endif
            !
         enddo
         !
      endif
      !
      if (verbose>=4) write(out,"('ML_b0_ZXY2/end')") 

  end subroutine ML_b0_ZXY2







  !
  ! Defining MEP function for H2CO molecule
  !
  function ML_MEP_zxy2_R_rho(x)  result(dst)

   real(ark),intent(in)   ::  x
   real(ark)              ::  dst(6)
   real(ark)              ::  cosrho
     !
     cosrho = cos(x) + 1.0_ark
     !
     dst(1)     = sum(molec%mep_params(1:5 )*cosrho**molec%mep_ind(1,1:5 ))
     dst(2)     = sum(molec%mep_params(6:10)*cosrho**molec%mep_ind(2,6:10))
     dst(3)     = dst(2)
     dst(4)     = sum(molec%mep_params(11:15)*cosrho**molec%mep_ind(4,11:15))
     dst(5)     = dst(4)
     dst(6)     = x
     !
  end function ML_MEP_zxy2_R_rho




  !
  ! Defining MEP function for H2CS molecule
  !
  function ML_MEP_zxy2_rho_coeff(x)  result(dst)

   real(ark),intent(in)   ::  x
   real(ark)              ::  dst(6)
   real(ark)              ::  cs,rad
     !
     rad=pi/180.0_ark
     ! 
     cs=1.0_ark+cos(x)
     !
     dst(1) = molec%mep_params(1)+molec%mep_params(2)*cs+molec%mep_params(3)*cs**2+molec%mep_params(4)*cs**3+molec%mep_params( 5)*cs**4
     dst(2) = molec%mep_params(6)+molec%mep_params(7)*cs+molec%mep_params(8)*cs**2+molec%mep_params(9)*cs**3+molec%mep_params(10)*cs**4
     dst(3) = dst(2)
     dst(4) = molec%mep_params(11)*rad+molec%mep_params(12)*cs+molec%mep_params(13)*cs**2+molec%mep_params(14)*cs**3+molec%mep_params(15)*cs**4
     dst(5) = dst(4)
     dst(6) = x
     !
  end function ML_MEP_zxy2_rho_coeff




  ! Here we define structural parameters a0 for ABCD molecule,
  !
  subroutine ML_b0_SOHF(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
     !
     integer(ik),intent(in) :: Npoints,Natoms
     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in) ,optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional :: rho_ref
     real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
     !
     real(ark)    :: transform(3,3),a0(molec%Natoms,3),CM_shift
     real(ark)    :: r_at(1:molec%ncoords),a0_ark(molec%Natoms,3)
     real(ark)    :: tau,alpha_SOF,alpha_HOS,alpha_HOF,phi,cosa,tau_2,cosphi
     !
     integer(ik)  :: i
      !
      if (verbose>=4) write(out,"('ML_b0_SOHF/start')") 
      !
      if (size(molec%req)/=3) then
        write(out,"('ML_b0_SOHF: Nbonds must be 3 in this routine, not  ',i9)") size(molec%req)
        stop 'ML_b0_SOHF: wrong Nbonds '
      endif 

      if (molec%Natoms/=4) then
        write(out,"('ML_b0_SOHF: Natoms must be 4 in this routine, not  ',i9)") molec%Natoms
        stop 'ML_b0_SOHF: wrong Natoms '
      endif 
      !
      r_at(:) = molec%local_eq(:)
      !
      if (molec%ncoords==7) then
        !
        r_at(1:7) = molec%local_eq(1:7)
        !
      else
        !
        write(out,"('ML_b0_SOHF: Ncoords must be 7 in this routine, not  ',i9)") molec%ncoords
        !stop 'ML_b0_SOHF: wrong Ncoords '
        !
      endif


      !
      alpha_SOF = r_at(4)
      alpha_HOF = r_at(5)
      alpha_HOS = r_at(6)
      !
      tau_2 = 1.0_ark-cos(alpha_SOF)**2-cos(alpha_HOF)**2-cos(alpha_HOS)**2 & 
               +2.0_ark*cos(alpha_SOF)*cos(alpha_HOF)*cos(alpha_HOS)
      !
      phi = sqrt(tau_2)
      !
      cosphi =  ( cos(alpha_HOS)-cos(alpha_HOF)*cos(alpha_SOF) )/( sin(alpha_HOF)*sin(alpha_SOF) )
      !
      if ( abs(cosphi)>1.0_ark+sqrt(small_) ) then 
         !
         write (out,"('ML_coordinate_transform_SOHF: cosphi>1: ',f18.8)") cosphi
         write (out,"('Consider change difftype ')")
         stop 'ML_coordinate_transform_SOHF - bad cosphi'
         !
      elseif ( abs(cosphi)>=1.0_ark) then 
         !
         phi = 0.0_ark
         !
      else 
         phi = acos(cosphi)
         !
      endif 
      !
      r_at(7) =  phi
      !
      call MLfromlocal2cartesian(1,r_at,a0_ark)      
      !
      a0 = a0_ark
      !
      call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,transform)
      !
      b0(:,:,0) = a0(:,:)
      !
      if (verbose>=3) then 
        write(out,"(i6)") molec%natoms
        !
        write(out,"(/'O',3x,3f14.8)") b0(1,:,0)
        write(out,"( 'F',3x,3f14.8)") b0(2,:,0)
        write(out,"( 'S',3x,3f14.8)") b0(3,:,0)
        write(out,"( 'H',3x,3f14.8)") b0(4,:,0)
      endif

      if (Npoints/=0) then
         ! 
         if (.not.present(rho_borders).or..not.present(rho_ref)) then  
            write(out,"('ML_b0_SOHF: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_SOHF: rho_borders or rho_ref not specified '
         endif
         !
         do i = 0,npoints
           !
           rho_ref = pi*0.5_ark
           !
           r_at(1:6) = molec%local_eq(1:6)
           !
           phi =rho_i(i)
           !
           alpha_SOF = r_at(4)
           alpha_HOF = r_at(5)
           alpha_HOS = r_at(6)
           ! 
           cosa = cos(alpha_HOF)*cos(alpha_SOF) + sin(alpha_HOF)*sin(alpha_SOF)*cos(phi)
           !
           alpha_HOF = acos(cosa)
           !
           r_at(6) = alpha_HOF
           !
           tau_2 = 1.0_ark-cos(alpha_SOF)**2-cos(alpha_HOF)**2-cos(alpha_HOS)**2 & 
                    +2.0_ark*cos(alpha_SOF)*cos(alpha_HOF)*cos(alpha_HOS)
           !
           !phi = sign(1.0_ark,phi-pi)*sqrt(tau_2)
           !
           r_at(7) = phi
           !
           call MLfromlocal2cartesian(1,r_at,a0_ark)      
           !
           a0 = a0_ark
           !
           call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,transform)
           !
           b0(:,:,i) = a0(:,:)
           !
           if (verbose>=4) then 
             write(out,"(i6)") molec%natoms
             !
             write(out,"(/'O',3x,3f14.8)") b0(1,:,i)
             write(out,"( 'F',3x,3f14.8)") b0(2,:,i)
             write(out,"( 'S',3x,3f14.8)") b0(3,:,i)
             write(out,"( 'H',3x,3f14.8)") b0(4,:,i)
             !
           endif
           !
         enddo
         !
      endif 
      !
      if (verbose>=4) write(out,"('ML_b0_SOHF/start')") 

  end subroutine ML_b0_SOHF




  function ML_coordinate_transform_SOHF(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct
    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: dsrc(size(src))
    real(ark)                 :: cosalpha,c_t(molec%Nmodes,molec%Nmodes)
    real(ark)                 :: alpha,sinrho,s6,alpha1,alpha2,alpha3,tau_2,r1,r2,r3,norm_2,sindelta
    real(ark)                 :: alpha_SOF,alpha_HOS,alpha_HOF,cosphi,phi,cosa
    !
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_SOHF/start')") 
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
       write (out,"('ML_coordinate_transform_SOHF: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'ML_coordinate_transform_SOHF - bad coord. type'
       !
    case('NORMAL')
       !
       if (direct) then 
           dst = src
       else
           dst = src
       endif
       !
    case('R-S-PHI')
       !
       if (direct) then
          !
          if (size(src)/=7) then
            !
            write(out,"('ML_coordinate_transform_SOHF, R-S-DELTA: Ncoords must be 7 in this routine, not  ',i9)") molec%ncoords
            stop 'ML_coordinate_transform_SOHF, R-S-DELTA: wrong Ncoords '
            !
          endif 
          ! 
          dst(1:3) = dsrc(1:3)
          !
          !dst(4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*dsrc(4)-dsrc(5)-dsrc(6) )
          !dst(5) = 1.0_ark/sqrt(2.0_ark)*(                 dsrc(5)-dsrc(6) )
          !
          alpha_SOF = src(4)
          alpha_HOF = src(5)
          alpha_HOS = src(6)
          !
          dst(4) = dsrc(4) !alpha_SOF
          dst(5) = dsrc(5) !alpha_HOF
          !
          cosphi =  ( cos(alpha_HOS)-cos(alpha_HOF)*cos(alpha_SOF) )/( sin(alpha_HOF)*sin(alpha_SOF) )
          !
          if ( abs(cosphi)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('ML_coordinate_transform_SOHF: cosphi>1: ',f18.8)") cosphi
             write (out,"('Consider change difftype ')")
             stop 'ML_coordinate_transform_SOHF - bad cosphi'
             !
          elseif ( abs(cosphi)>=1.0_ark) then 
             !
             phi = 0.0_ark
             !
          else 
             phi = acos(cosphi)
             !
             if (src(7)>0) phi = phi + pi
             !
          endif 
          !
          dst(6) =  phi
          !
      else    ! not direct
          !
          !
          if (size(dst)/=7) then
            !
            write(out,"('ML_coordinate_transform_SOHF, R-S-DELTA: Ncoords must be 7 in this routine, not  ',i9)") molec%ncoords
            stop 'ML_coordinate_transform_SOHF, R-S-DELTA : wrong Ncoords '
            !
          endif
          !
          dst(1:3) = src(1:3)+molec%local_eq(1:3)
          !
          alpha_SOF = src(4)+molec%local_eq(4)
          alpha_HOF = src(5)+molec%local_eq(5)
          !
          phi = src(6)
          ! 
          cosa = cos(alpha_HOF)*cos(alpha_SOF) + sin(alpha_HOF)*sin(alpha_SOF)*cos(phi)
          !
          if ( abs(cosa)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('ML_coordinate_transform_SOHF: cosa>1: ',f18.8)") cosa
             write (out,"('Consider change difftype ')")
             stop 'ML_coordinate_transform_SOHF - bad cosa'
             !
          elseif ( abs(cosa)>=1.0_ark) then 
             !
             alpha_HOS = 0.0_ark
             !
          else 
             alpha_HOS = acos(cosa)
             !
          endif 
          !
          dst(4) = alpha_SOF
          dst(5) = alpha_HOF
          dst(6) = alpha_HOS

          tau_2 = 1.0_ark-cos(alpha_SOF)**2-cos(alpha_HOF)**2-cos(alpha_HOS)**2 & 
                   +2.0_ark*cos(alpha_SOF)*cos(alpha_HOF)*cos(alpha_HOS)
          !
          if ( tau_2<-sqrt(small_) ) then 
             !
             write (out,"('MLcoordinate_transform_func: tau**2<0: ',f18.8)") tau_2
             stop 'MLcoordinate_transform_func - tau**2<0'
             !
          elseif ( tau_2<0.0_ark) then 
             dst(7) = 0
          else 
             dst(7) = sign(1.0_ark,phi-pi)*sqrt(tau_2)
          endif
          !
       endif 
       !
    end select
    !
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_SOHF/end')") 
    !
    !
  end function ML_coordinate_transform_SOHF


  ! Here we define the symmetry transformation of the Nmodes coordinates according the symmetry operations
  !
  subroutine ML_symmetry_transformation_ZXY2(ioper,natoms,src,dst)
    !
    integer(ik),intent(in)    :: ioper  ! group operation  
    integer(ik),intent(in)    :: natoms
    real(ark),intent(in)      :: src(natoms)
    real(ark),intent(out)     :: dst(natoms)
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_ZXY2/start')") 
    !
    nsrc = size(src)
    !
    select case(trim(molec%coords_transform))
    !
    case default
       write (out,"('ML_coordinate_transform_ZXY2: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'ML_coordinate_transform_ZXY2 - bad coord. type'
       !
    case('SYMM-C2V')
       !
       select case(trim(molec%symmetry))
          !
       case default
          write (out,"('ML_symmetry_transformation_ZXY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_ZXY2 - bad symm. type'
          !
       case('CS','CS(M)')
          !
          select case(ioper)
            !
          case (1) ! E 
            !
            dst = src
            !
          case (2) ! (12)
            !
            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            dst(4) = src(4)
            dst(5) =-src(5)
            dst(6) =-src(6)

          case default
            !
          end select
          !
       end select
       !
    case('R-THETA-TAU','R-THETA-TAU-MEP','R-THETA-DELTA')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_ZXY2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_ZXY2 - bad symm. type'
          !
       case('C2V','C2V(M)')
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
            dst(6) = -src(6)

          case (3) ! (E*)

            dst(1) = src(1)
            dst(2) = src(2)
            dst(3) = src(3)
            dst(4) = src(4)
            dst(5) = src(5)
            dst(6) = -src(6)

          case (4) ! (12)*

            dst(1) = src(1)
            dst(2) = src(3)
            dst(3) = src(2)
            dst(4) = src(5)
            dst(5) = src(4)
            dst(6) = src(6)

          case default

            write (out,"('ML_symmetry_transformation_ZXY2: operation ',i8,' unknown')") ioper
            stop 'ML_symmetry_transformation_ZXY2 - bad operation. type'
            !
          end select 
          !
       end select
        !
    end select
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY2/end')") 
    !
  end subroutine ML_symmetry_transformation_ZXY2


  ! Here we find the rotational symmetry from J,K
  !
  subroutine ML_rotsymmetry_ZXY2(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_ZXY2/start')") 
    !
    !
    select case(trim(molec%symmetry))
    case default
       write (out,"('ML_rotsymmetry_ZXY2: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_rotsymmetry_ZXY2 - bad symm. type'
       !
    case('C2V','C2V(M)')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma = 1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma = 2 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma = 3 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma = 4 !; return
      !
    end select

    !
    !
    if (verbose>=5) write(out,"('ML_rotsymmetry_ZXY2/end')") 
    !
    !
  end subroutine ML_rotsymmetry_ZXY2




end module mol_zxy2

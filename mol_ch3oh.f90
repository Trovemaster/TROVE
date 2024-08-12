!
!  This unit defines all specific routines for a fouratomic molecule of ch3oh type
!
module mol_ch3oh
  use accuracy
  use moltype
  use lapack
  use pot_ch3oh, only : ML_CH3OH_MEP_geometry

  implicit none

  public ML_b0_ch3oh,ML_coordinate_transform_ch3oh,ML_symmetry_transformation_ch3oh,ML_rotsymmetry_CH3OH
  !
  private
 
  integer(ik), parameter :: verbose     = 5                       ! Verbosity level
  !
  contains
  !
  ! Procedures to define  ch3oh 
  ! 
  function ML_coordinate_transform_ch3oh(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: dsrc(size(src))
    real(ark)                 :: tau,rref(12,0:4),r_eq(12),r_t(12),xi_eq(12)
    integer(ik)               :: i0,nsrc, k(0:4),n(0:4),icoord
    real(ark)                 :: beta(3),alpha(3),cosa(3),phi1,phi2,phi3,chi1,chi2,chi3
    real(ark)                 :: alpha0(3),cosa0(3),phi(3),chi(3),coschi(3)
    real(ark)                 :: t1,t2,t3,theta12,theta23,theta13,a1,a2,tbar
    !
    if (verbose>=7) write(out,"('ML_coordinate_transform_ch3oh/start')") 
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
    case('R-ALPHA-THETA-TAU')
      !
      if (direct) then 
        !
        !for stretches and 'alpha' bends just subtract equilibrium coordinates
        dst(1:12) = src(1:12)-molec%local_eq(1:12)
        !
        t1 = src(10)
        t2 = src(11)
        t3 = src(12)
        !
        ! subtract equilbrium theta values to make a1/a2 zero at equilibrium
        ! and ensure consistent transfroms
        !
        if (t2-t1<small_) t2 = t2 + 2.0_ark*pi
        if (t3-t2<small_) t3 = t3 + 2.0_ark*pi
        !
        theta12 = mod(t2-t1+2.0_ark*pi,2.0_ark*pi)
        theta23 = mod(t3-t2+2.0_ark*pi,2.0_ark*pi)
        theta13 = mod(t1-t3+2.0_ark*pi,2.0_ark*pi)
        !
        a1  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
        a2  = (                   theta13 - theta12 )/sqrt(2.0_ark)
        !
        tbar = (t1 + t2 + t3-2.0_ark*pi)/3.0_ark
        !
        dst(10) = a1
        dst(11) = a2
        dst(12) = tbar 
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        dst(1:12) = src(1:12)+molec%local_eq(1:12)
        !
        A1 = src(10) 
        A2 = src(11) 
        tbar = src(12)
        !
        t1 =                    tbar+sqrt(2.0_ark)/3.0_ark*A2 
        t2 = 2.0_ark/3.0_ark*Pi+tbar-sqrt(6.0_ark)/6.0_ark*A1-sqrt(2.0_ark)/6.0_ark*A2
        t3 = 4.0_ark/3.0_ark*Pi+tbar+sqrt(6.0_ark)/6.0_ark*A1-sqrt(2.0_ark)/6.0_ark*A2
        !
        dst(10) =  mod(t1+4.0_ark*pi,4.0_ark*pi)
        dst(11) =  mod(t2+4.0_ark*pi,4.0_ark*pi)
        dst(12) =  mod(t3+4.0_ark*pi,4.0_ark*pi)
        !
        continue
        !
      endif
      !
    case('R-TAU-CHI')
       !
       !
       call ML_CH3OH_MEP_geometry(tau,r_eq(1),r_eq(2),r_eq(6),r_eq(3),r_eq(4),r_eq(5),beta(1),beta(2),beta(3),&
                                          alpha(1),alpha(2),alpha(3),phi(1),phi(2),phi(3),chi(1),chi(2),chi(3))
       !
       phi1   = src(10)
       phi2   = src(11)
       phi3   = src(12)
       !
       if (phi2>2.0_ark/3.0_ark*pi)  phi2 = phi2-2.0_ark*pi
       !
       tau = 1.0_ark/3.0_ark*( phi1+phi2+phi3 )
       !
       xi_eq(1:2) = r_eq(1:2)

       xi_eq( 3) = 1.0_ark/sqrt(3.0_ark)*(         r_eq(3)+r_eq(4)+r_eq(5) )
       xi_eq( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*r_eq(3)-r_eq(4)-r_eq(5) )
       xi_eq( 5) = 1.0_ark/sqrt(2.0_ark)*(                 r_eq(4)-r_eq(5) )
       !
       xi_eq(6) = r_eq(6)
       !
       xi_eq( 7) = 1.0_ark/sqrt(3.0_ark)*(         beta(1)+beta(2)+beta(3) )
       xi_eq( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*beta(1)-beta(2)-beta(3) )
       xi_eq( 9) = 1.0_ark/sqrt(2.0_ark)*(                 beta(2)-beta(3) )
       !
       xi_eq(10) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*chi(1)-chi(2)-chi(3) ) 
       xi_eq(11) = 1.0_ark/sqrt(2.0_ark)*(                chi(2)-chi(3) ) 
       xi_eq(12) = tau
       !
       if (direct) then
          !
          dst(1:2) = src(1:2)
          !
          dst( 3) = 1.0_ark/sqrt(3.0_ark)*(         src(3)+src(4)+src(5) )
          dst( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(3)-src(4)-src(5) )
          dst( 5) = 1.0_ark/sqrt(2.0_ark)*(                src(4)-src(5) )
          !
          dst( 6) = src(6)
          !
          dst( 7) = 1.0_ark/sqrt(3.0_ark)*(         src(7)+src(8)+src(9) )
          dst( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(7)-src(8)-src(9) )
          dst( 9) = 1.0_ark/sqrt(2.0_ark)*(                src(8)-src(9) )
          !
          beta(1) = src(7)
          beta(2) = src(8)
          beta(3) = src(9)
          !
          cosa(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi2 - phi3 )
          cosa(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi3 - phi1 )
          cosa(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi1 - phi2 )
          !
          if ( any(abs(cosa(:))>1.0_ark+sqrt(small_))) then 
              !
              write (out,"('ML_coordinate_transform_ch3oh: cosalpha>1: ',3f18.8)") cosa(1:3)
              write (out,"('Consider change coordinate type ')")
              stop 'MLcoordinate_transform_func - bad cos(alpha)'
          endif 
          !
          cosa(:) = max(-1.0_ark,cosa(:))
          cosa(:) = min( 1.0_ark,cosa(:))
          !
          alpha(:) = acos(cosa(:))
          !
          coschi(1) =  ( cos(alpha(1))-cos(beta(2))*cos(beta(3)) )/( sin(beta(2))*sin(beta(3)) )
          coschi(2) =  ( cos(alpha(2))-cos(beta(3))*cos(beta(1)) )/( sin(beta(3))*sin(beta(1)) )
          coschi(3) =  ( cos(alpha(3))-cos(beta(2))*cos(beta(1)) )/( sin(beta(2))*sin(beta(1)) )
          !
          chi(:) = acos(coschi(:))
          !
          dst(10) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*chi(1)-chi(2)-chi(3) )
          dst(11) = 1.0_ark/sqrt(2.0_ark)*(                chi(2)-chi(3) )
          dst(12) = tau
          !
          dst(1:11) = dst(1:11)-xi_eq(1:11)
          !
      else ! not direct
          !
          r_t(12) = tau
          !
          !xi_eq(7) = xi_eq(7)/sqrt(2.0_ark)-(alpha0(1)+alpha0(2)+alpha0(3))/sqrt(6.0_ark)
          !
          r_t(1:11) = src(1:11)+xi_eq(1:11)
          !
          dst(1:2) = r_t(1:2)
          !
          dst(3) =(sqrt(2.0_ark)*r_t(3)+2.0_ark*r_t(4)                     )/sqrt(6.0_ark)
          dst(4) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)+sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          dst(5) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)-sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          !
          dst(6) = r_t(6)
          !
          dst( 7) =(sqrt(2.0_ark)*r_t(7)+2.0_ark*r_t(8)                     )/sqrt(6.0_ark)
          dst( 8) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)+sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          dst( 9) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)-sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          !
          beta(1) = dst(7)
          beta(2) = dst(8)
          beta(3) = dst(9)
          !
          call find_phi_from_beta(beta(1:3),r_t(10:12),phi(1:3))
          !
          dst(10:12) = phi(1:3)
          !
      endif
       !
    case('R-ALPHA-S-TAU-REF')
       !
       rref( 1,0:4)  = molec%force( 1: 5) 
       rref( 2,0:4)  = molec%force( 6:10) 
       rref( 3,0:4)  = molec%force(11:15)  
       rref( 4,0:4)  = molec%force(16:20)
       rref( 5,0:4)  = molec%force(21:25) 
       !
       rref( 6,0:4)  = molec%force(26:30)/rad
       rref( 7,0:4)  = molec%force(31:35)/rad 
       rref( 8,0:4)  = molec%force(36:40)/rad  
       rref( 9,0:4)  = molec%force(41:45)/rad 
       !
       rref(10,0:4)  = molec%force(46:50)/rad 
       rref(11,0:4)  = molec%force(51:55)/rad  
       rref(12,0:4)  = molec%force(56:60)/rad 
       !
       !
       if (direct) then
          !
          dsrc(10) = src(10)
          if (dsrc(10)>2.0_ark*pi-sqrt(small_)) dsrc(10) = dsrc(10)-2.0_ark*pi
          if (dsrc(10)<-sqrt(small_)) dsrc(10) = dsrc(10)+2.0_ark*pi
          !
          dsrc(11) = src(11) 
          if (dsrc(11)-rref(11,1)>2.0_ark*pi-sqrt(small_)) dsrc(11) = dsrc(11)-2.0_ark*pi
          if (dsrc(11)-rref(11,1)<-sqrt(small_)) dsrc(11) = dsrc(11)+2.0_ark*pi
          dsrc(12) = src(12) 
          if (dsrc(12)-rref(12,1)>2.0_ark*pi-sqrt(small_)) dsrc(12) = dsrc(12)-2.0_ark*pi
          if (dsrc(12)-rref(12,1)<-sqrt(small_)) dsrc(12) = dsrc(12)+2.0_ark*pi
          !
          tau = 1.0_ark/3.0_ark*( dsrc(10)+dsrc(11)+dsrc(12) )
          !
          !tau = tau - molec%local_eq(10)
          !
          !tau = src(10)
          !
       else  ! not direct
          !
          tau = src(12)
          !
       endif
       !
       i0 = 0 
       !
       do icoord = 1,9
         !
         n(0:4) = molec%pot_ind(11,i0+1:i0+5)
         k(0:4) = molec%pot_ind(12,i0+1:i0+5)
         !
         if (all(n(:)==0)) then 
           !
           xi_eq(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*cos(real(k(1:4),ark)*tau))
           !
         else
           !
           xi_eq(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*sin(real(n(1:4),ark)*tau))
           !
         endif 
         !
         i0 = i0+5
         !
       enddo
       !
       xi_eq(10:11) = rref(10:11,0)
       xi_eq(12) = rref(12,0) + tau
       !
       !xi_eq(1:6) = r_t(1:6)
       !
       !xi_eq(8) =(sqrt(2.0_ark)*r_t(7)+2.0_ark*r_t(8)                     )/sqrt(6.0_ark)
       !xi_eq(9) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)+sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
       !xi_eq(7) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)-sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
       !
       !xi_eq(10:11) = r_t(10:11) !tau+molec%local_eq(10)-molec%local_eq(11)
       !xi_eq(12) = 0
       !
       !xi_eq(11) = tau+molec%local_eq(10)-molec%local_eq(12)
       !
       if (direct) then
          !
          dst(1:2) = src(1:2)
          !
          dst( 3) = 1.0_ark/sqrt(3.0_ark)*(         src(3)+src(4)+src(5) )
          dst( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(3)-src(4)-src(5) )
          dst( 5) = 1.0_ark/sqrt(2.0_ark)*(                src(4)-src(5) )
          !
          dst( 6) = src(6)
          !
          dst( 7) = 1.0_ark/sqrt(3.0_ark)*(         src(7)+src(8)+src(9) )
          dst( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(7)-src(8)-src(9) )
          dst( 9) = 1.0_ark/sqrt(2.0_ark)*(                src(8)-src(9) )
          !
          !dsrc(10) = src(10)
          !
          !dsrc(11) = src(11) 
          !if (dsrc(11)-molec%local_eq(11)>2.0_ark*pi-sqrt(small_)) dsrc(11) = dsrc(11)-2.0_ark*pi
          !if (dsrc(11)-molec%local_eq(11)<-sqrt(small_)) dsrc(11) = dsrc(11)+2.0_ark*pi
          !dsrc(12) = src(12) 
          !if (dsrc(12)-molec%local_eq(12)>2.0_ark*pi-sqrt(small_)) dsrc(12) = dsrc(12)-2.0_ark*pi
          !if (dsrc(12)-molec%local_eq(12)<-sqrt(small_)) dsrc(12) = dsrc(12)+2.0_ark*pi
          !
          !dsrc(11) = mod(src(11),rref(11,0))+rref(11,0)
          !dsrc(12) = mod(src(12),rref(12,0))+rref(12,0)
          !
          dst(10) = 1.0_ark/6.0_ark*( 2.0_ark*dsrc(10)-dsrc(11)-dsrc(12) )
          if (dst(10)>pi) dst(10) = 2.0_ark*pi-dst(10)
          dst(11) = 1.0_ark/2.0_ark*(                  dsrc(11)-dsrc(12) )
          if (dst(11)>pi) dst(11) =-(2.0_ark*pi-dst(11))
          dst(12) = tau
          if (dst(12)>pi) dst(12) = 2.0_ark*pi-dst(12)
          !
          dst(1:11) = dst(1:11)-xi_eq(1:11)
          !

          !
      else ! not direct
          !
          r_t(12) = tau
          !
          r_t(1:11) = src(1:11)+xi_eq(1:11)
          !
          dst(1:2) = r_t(1:2)
          !
          dst(3) =(sqrt(2.0_ark)*r_t(3)+2.0_ark*r_t(4)                     )/sqrt(6.0_ark)
          dst(4) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)+sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          dst(5) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)-sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          !
          dst(6) = r_t(6)
          !
          dst(7) =(sqrt(2.0_ark)*r_t(7)+2.0_ark*r_t(8)                     )/sqrt(6.0_ark)
          dst(8) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)+sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          dst(9) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)-sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          !
          dst(10) =(r_t(12)+2.0_ark*r_t(10)        )
          dst(11) =(r_t(12)-        r_t(10)+r_t(11))
          dst(12) =(r_t(12)-        r_t(10)-r_t(11))
          !
          if (dst(10)>2.0_ark*pi-sqrt(small_)) dst(10) = dst(10)-2.0_ark*pi
          if (dst(10)<-sqrt(small_)) dst(10) = dst(10)+2.0_ark*pi
          !if (dst(10)>pi) dst(10) = 2.0_ark*pi-dst(10)
          !
          if (dst(11)-molec%local_eq(11)>2.0_ark*pi-sqrt(small_)) dst(11) = dst(11)-2.0_ark*pi
          if (dst(11)-molec%local_eq(11)<-sqrt(small_)) dst(11) = dst(11)+2.0_ark*pi
          !if (dst(11)-molec%local_eq(11)>pi) dst(11) = 2.0_ark*pi-dst(11)
          !
          if (dst(12)-molec%local_eq(12)>2.0_ark*pi-sqrt(small_)) dst(12) = dst(12)-2.0_ark*pi
          if (dst(12)-molec%local_eq(12)<-sqrt(small_)) dst(12) = dst(12)+2.0_ark*pi
          !if (dst(12)-molec%local_eq(12)>pi) dst(12) = -(2.0_ark*pi-dst(12))
          !
      endif
       !
    case('R-TAU-BOWMAN-REF')
       !
       if (direct) then
         !
         alpha(3) = src( 9)
         alpha(2) = src(11)
         !
         beta(1) = src(7)
         beta(2) = src(8)
         beta(3) = src(10)
         !
         cosa(2) =  ( cos(alpha(2))-cos(beta(3))*cos(beta(1)) )/( sin(beta(3))*sin(beta(1)) )
         cosa(3) =  ( cos(alpha(3))-cos(beta(2))*cos(beta(1)) )/( sin(beta(2))*sin(beta(1)) )
         !
         chi2 =-acos(cosa(2))
         chi3 = acos(cosa(3))
         !
         chi1 = 2.0_ark*pi-(-chi2+chi3)
         !
         tau = src(12)+(chi2+chi3)/3.0_ark
         !
       else
         !
         tau = src(12)
         !
       endif
       !
       rref( 1,0:4)  = molec%force( 1: 5) 
       rref( 2,0:4)  = molec%force( 6:10) 
       rref( 3,0:4)  = molec%force(11:15)  
       rref( 4,0:4)  = molec%force(16:20)
       rref( 5,0:4)  = molec%force(21:25) 
       !
       rref( 6,0:4)  = molec%force(26:30)/rad
       rref( 7,0:4)  = molec%force(31:35)/rad 
       rref( 8,0:4)  = molec%force(36:40)/rad  
       rref( 9,0:4)  = molec%force(41:45)/rad 
       !
       rref(10,0:4)  = molec%force(46:50)/rad 
       rref(11,0:4)  = molec%force(51:55)/rad  
       rref(12,0:4)  = molec%force(56:60)/rad 
       !
       if (direct) then
         !
         alpha(3) = src( 9)
         alpha(2) = src(11)
         !
         beta(1) = src(7)
         beta(2) = src(8)
         beta(3) = src(10)
         !
         cosa(2) =  ( cos(alpha(2))-cos(beta(3))*cos(beta(1)) )/( sin(beta(3))*sin(beta(1)) )
         cosa(3) =  ( cos(alpha(3))-cos(beta(2))*cos(beta(1)) )/( sin(beta(2))*sin(beta(1)) )
         !
         chi2 =-acos(cosa(2))
         chi3 = acos(cosa(3))
         !
         chi1 = 2.0_ark*pi-(-chi2+chi3)
         !
         tau = src(12)+(chi2+chi3)/3.0_ark
         !
       else
         !
         tau = src(12)
         !
       endif
       !
       i0 = 0 
       !
       do icoord = 1,9
         !
         n(0:4) = molec%pot_ind(11,i0+1:i0+5)
         k(0:4) = molec%pot_ind(12,i0+1:i0+5)
         !
         if (all(n(:)==0)) then 
           !
           xi_eq(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*cos(real(k(1:4),ark)*tau))
           !
         else
           !
           xi_eq(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*sin(real(n(1:4),ark)*tau))
           !
         endif 
         !
         i0 = i0+5
         !
       enddo
       !
       xi_eq(10:11) = rref(10:11,0)
       xi_eq(12) = rref(12,0) + tau
       !
       beta(1) = (sqrt(2.0_ark)*xi_eq(7)+2.0_ark*xi_eq(8)                       )/sqrt(6.0_ark)
       beta(2) = (sqrt(2.0_ark)*xi_eq(7)-        xi_eq(8)+sqrt(3.0_ark)*xi_eq(9))/sqrt(6.0_ark)
       beta(3) = (sqrt(2.0_ark)*xi_eq(7)-        xi_eq(8)-sqrt(3.0_ark)*xi_eq(9))/sqrt(6.0_ark)
       !
       phi1   = (xi_eq(12)+2.0_ark*xi_eq(10)          )
       phi2   = (xi_eq(12)-        xi_eq(10)+xi_eq(11))
       phi3   = (xi_eq(12)-        xi_eq(10)-xi_eq(11))
       !
       cosa0(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi2 - phi3 )
       cosa0(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi3 - phi1 )
       cosa0(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi1 - phi2 )
       !
       if ( any(abs(cosa0(:))>1.0_ark+sqrt(small_))) then 
           !
           write (out,"('ML_coordinate_transform_ch3oh: cosalpha0>1: ',3f18.8)") cosa0(1:3)
           write (out,"('Consider change coordinate type ')")
           stop 'MLcoordinate_transform_func - bad cos(alpha)'
       endif 
       !
       cosa0(:) = max(-1.0_ark,cosa0(:))
       cosa0(:) = min( 1.0_ark,cosa0(:))
       !
       alpha0(:) = acos(cosa0(:))
       !
       !xi_eq(10) = (alpha0(2)+alpha0(3))/sqrt(2.0_ark)
       !xi_eq(11) = (alpha0(2)-alpha0(3))/sqrt(2.0_ark)
       !
       xi_eq(10) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha0(1)-alpha0(2)-alpha0(3) ) 
       xi_eq(11) = 1.0_ark/sqrt(2.0_ark)*(                   alpha0(2)-alpha0(3) ) 
       !
       if (direct) then
          !
          dst(1:2) = src(1:2)
          !
          dst( 3) = 1.0_ark/sqrt(3.0_ark)*(         src(3)+src(4)+src(5) )
          dst( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(3)-src(4)-src(5) )
          dst( 5) = 1.0_ark/sqrt(2.0_ark)*(                src(4)-src(5) )
          !
          dst( 6) = src(6)
          !
          dst( 7) = 1.0_ark/sqrt(3.0_ark)*(         src(7)+src(8)+src(10) )
          dst( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(7)-src(8)-src(10) )
          dst( 9) = 1.0_ark/sqrt(2.0_ark)*(                src(8)-src(10) )
          !
          !dsrc(10) = src(10)
          !
          !dsrc(11) = src(11) 
          !if (dsrc(11)-molec%local_eq(11)>2.0_ark*pi-sqrt(small_)) dsrc(11) = dsrc(11)-2.0_ark*pi
          !if (dsrc(11)-molec%local_eq(11)<-sqrt(small_)) dsrc(11) = dsrc(11)+2.0_ark*pi
          !dsrc(12) = src(12) 
          !if (dsrc(12)-molec%local_eq(12)>2.0_ark*pi-sqrt(small_)) dsrc(12) = dsrc(12)-2.0_ark*pi
          !if (dsrc(12)-molec%local_eq(12)<-sqrt(small_)) dsrc(12) = dsrc(12)+2.0_ark*pi
          !
          !dsrc(11) = mod(src(11),rref(11,0))+rref(11,0)
          !dsrc(12) = mod(src(12),rref(12,0))+rref(12,0)
          !
          !dst(7)   = dst(7)  /sqrt(2.0_ark)-(alpha(1) +alpha(2) +alpha(3) )/sqrt(6.0_ark)
          !xi_eq(7) = xi_eq(7)/sqrt(2.0_ark)-(alpha0(1)+alpha0(2)+alpha0(3))/sqrt(6.0_ark)
          !
          !dst(10) = (src( 9)+src(11))/sqrt(2.0_ark)
          !dst(11) = (src( 9)-src(11))/sqrt(2.0_ark)
          !
          beta(1) = src(7)
          beta(2) = src(8)
          beta(3) = src(10)
          !
          cosa(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( chi1 )
          !
          alpha(1) = acos(cosa(1))
          !
          dst(10) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha(1)-alpha(2)-alpha(3) )
          dst(11) = 1.0_ark/sqrt(2.0_ark)*(                  alpha(2)-alpha(3) )
          dst(12) = tau
          !
          dst(1:11) = dst(1:11)-xi_eq(1:11)
          !
      else ! not direct
          !
          r_t(12) = tau
          !
          !xi_eq(7) = xi_eq(7)/sqrt(2.0_ark)-(alpha0(1)+alpha0(2)+alpha0(3))/sqrt(6.0_ark)
          !
          r_t(1:11) = src(1:11)+xi_eq(1:11)
          !
          dst(1:2) = r_t(1:2)
          !
          dst(3) =(sqrt(2.0_ark)*r_t(3)+2.0_ark*r_t(4)                     )/sqrt(6.0_ark)
          dst(4) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)+sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          dst(5) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)-sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          !
          dst(6) = r_t(6)
          !
          dst( 7) =(sqrt(2.0_ark)*r_t(7)+2.0_ark*r_t(8)                     )/sqrt(6.0_ark)
          dst( 8) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)+sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          dst(10) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)-sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          !
          beta(1) = dst(7)
          beta(2) = dst(8)
          beta(3) = dst(10)
          !
          call find_alpha_from_beta(beta(1:3),r_t(10:12),dst(11),dst(9),dst(12))
          !
      endif

       !
    case('R-1TAU-2BETA-REF')
       !
       rref( 1,0:4)  = molec%force( 1: 5) 
       rref( 2,0:4)  = molec%force( 6:10) 
       rref( 3,0:4)  = molec%force(11:15)  
       rref( 4,0:4)  = molec%force(16:20)
       rref( 5,0:4)  = molec%force(21:25) 
       !
       rref( 6,0:4)  = molec%force(26:30)/rad
       rref( 7,0:4)  = molec%force(31:35)/rad 
       rref( 8,0:4)  = molec%force(36:40)/rad  
       rref( 9,0:4)  = molec%force(41:45)/rad 
       !
       rref(10,0:4)  = molec%force(46:50)/rad 
       rref(11,0:4)  = molec%force(51:55)/rad  
       rref(12,0:4)  = molec%force(56:60)/rad 
       !
       !
       tau = src(12)
       !
       i0 = 0 
       !
       do icoord = 1,9
         !
         n(0:4) = molec%pot_ind(11,i0+1:i0+5)
         k(0:4) = molec%pot_ind(12,i0+1:i0+5)
         !
         if (all(n(:)==0)) then 
           !
           xi_eq(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*cos(real(k(1:4),ark)*tau))
           !
         else
           !
           xi_eq(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*sin(real(n(1:4),ark)*tau))
           !
         endif 
         !
         i0 = i0+5
         !
       enddo
       !
       xi_eq(10:11) = rref(10:11,0)
       xi_eq(12) = rref(12,0) + tau
       !
       beta(1) = (sqrt(2.0_ark)*xi_eq(7)+2.0_ark*xi_eq(8)                       )/sqrt(6.0_ark)
       beta(2) = (sqrt(2.0_ark)*xi_eq(7)-        xi_eq(8)+sqrt(3.0_ark)*xi_eq(9))/sqrt(6.0_ark)
       beta(3) = (sqrt(2.0_ark)*xi_eq(7)-        xi_eq(8)-sqrt(3.0_ark)*xi_eq(9))/sqrt(6.0_ark)
       !
       phi1   = (xi_eq(12)+2.0_ark*xi_eq(10)          )
       phi2   = (xi_eq(12)-        xi_eq(10)+xi_eq(11))
       phi3   = (xi_eq(12)-        xi_eq(10)-xi_eq(11))
       !
       cosa0(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi2 - phi3 )
       cosa0(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi3 - phi1 )
       cosa0(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi1 - phi2 )
       !
       if ( any(abs(cosa0(:))>1.0_ark+sqrt(small_))) then 
           !
           write (out,"('ML_coordinate_transform_ch3oh: cosalpha0>1: ',3f18.8)") cosa0(1:3)
           write (out,"('Consider change coordinate type ')")
           stop 'MLcoordinate_transform_func - bad cos(alpha)'
       endif 
       !
       cosa0(:) = max(-1.0_ark,cosa0(:))
       cosa0(:) = min( 1.0_ark,cosa0(:))
       !
       alpha0(:) = acos(cosa0(:))
       !
       !xi_eq(10) = (alpha0(2)+alpha0(3))/sqrt(2.0_ark)
       !xi_eq(11) = (alpha0(2)-alpha0(3))/sqrt(2.0_ark)
       !
       xi_eq(10) = alpha0(3)
       xi_eq(11) = alpha0(2)
       !
       if (direct) then
          !
          dst(1:2) = src(1:2)
          !
          dst( 3) = 1.0_ark/sqrt(3.0_ark)*(         src(3)+src(4)+src(5) )
          dst( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(3)-src(4)-src(5) )
          dst( 5) = 1.0_ark/sqrt(2.0_ark)*(                src(4)-src(5) )
          !
          dst( 6) = src(6)
          !
          dst( 7) = 1.0_ark/sqrt(3.0_ark)*(         src(7)+src(8)+src(10) )
          dst( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(7)-src(8)-src(10) )
          dst( 9) = 1.0_ark/sqrt(2.0_ark)*(                src(8)-src(10) )
          !
          !dsrc(10) = src(10)
          !
          !dsrc(11) = src(11) 
          !if (dsrc(11)-molec%local_eq(11)>2.0_ark*pi-sqrt(small_)) dsrc(11) = dsrc(11)-2.0_ark*pi
          !if (dsrc(11)-molec%local_eq(11)<-sqrt(small_)) dsrc(11) = dsrc(11)+2.0_ark*pi
          !dsrc(12) = src(12) 
          !if (dsrc(12)-molec%local_eq(12)>2.0_ark*pi-sqrt(small_)) dsrc(12) = dsrc(12)-2.0_ark*pi
          !if (dsrc(12)-molec%local_eq(12)<-sqrt(small_)) dsrc(12) = dsrc(12)+2.0_ark*pi
          !
          !dsrc(11) = mod(src(11),rref(11,0))+rref(11,0)
          !dsrc(12) = mod(src(12),rref(12,0))+rref(12,0)
          !
          !dst(7)   = dst(7)  /sqrt(2.0_ark)-(alpha(1) +alpha(2) +alpha(3) )/sqrt(6.0_ark)
          !xi_eq(7) = xi_eq(7)/sqrt(2.0_ark)-(alpha0(1)+alpha0(2)+alpha0(3))/sqrt(6.0_ark)
          !
          !dst(10) = (src( 9)+src(11))/sqrt(2.0_ark)
          !dst(11) = (src( 9)-src(11))/sqrt(2.0_ark)
          !
          dst(10) = src( 9)
          dst(11) = src(11)
          dst(12) = tau
          !
          dst(1:11) = dst(1:11)-xi_eq(1:11)
          !
      else ! not direct
          !
          r_t(12) = tau
          !
          !xi_eq(7) = xi_eq(7)/sqrt(2.0_ark)-(alpha0(1)+alpha0(2)+alpha0(3))/sqrt(6.0_ark)
          !
          r_t(1:11) = src(1:11)+xi_eq(1:11)
          !
          dst(1:2) = r_t(1:2)
          !
          dst(3) =(sqrt(2.0_ark)*r_t(3)+2.0_ark*r_t(4)                     )/sqrt(6.0_ark)
          dst(4) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)+sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          dst(5) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)-sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          !
          dst(6) = r_t(6)
          !
          dst( 7) =(sqrt(2.0_ark)*r_t(7)+2.0_ark*r_t(8)                     )/sqrt(6.0_ark)
          dst( 8) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)+sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          dst(10) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)-sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)

          !dst( 9) = (r_t(10)+r_t(11))/sqrt(2.0_ark)
          !dst(11) = (r_t(10)-r_t(11))/sqrt(2.0_ark)
          !
          dst( 9) = r_t(10)
          dst(11) = r_t(11)
          dst(12) = tau
          !
          !call find_alpha_for_ch3oh(dst(7:9),r_t(10:12),dst(10:12))
          !call find_alpha_for_ch3oh_II(r_t(7:12),dst(7:12))
          !
          !dst(10) =(r_t(12)+2.0_ark*r_t(10)        )
          !dst(11) =(r_t(12)-        r_t(10)+r_t(11))
          !dst(12) =(r_t(12)-        r_t(10)-r_t(11))
          !
          !dst(10:12) = mod(dst(10:12)+2.0_ark*pi,2.0_ark*pi)+molec%local_eq(10:12)
          !
      endif
      !
    case('R-ALPHA-S-TAU-BETA-REF')
       !
       rref( 1,0:4)  = molec%force( 1: 5) 
       rref( 2,0:4)  = molec%force( 6:10) 
       rref( 3,0:4)  = molec%force(11:15)  
       rref( 4,0:4)  = molec%force(16:20)
       rref( 5,0:4)  = molec%force(21:25) 
       !
       rref( 6,0:4)  = molec%force(26:30)/rad
       rref( 7,0:4)  = molec%force(31:35)/rad 
       rref( 8,0:4)  = molec%force(36:40)/rad  
       rref( 9,0:4)  = molec%force(41:45)/rad 
       !
       rref(10,0:4)  = molec%force(46:50)/rad 
       rref(11,0:4)  = molec%force(51:55)/rad  
       rref(12,0:4)  = molec%force(56:60)/rad 
       !
       !
       if (direct) then
          !
          !dsrc(10) = src(10)
          !if (dsrc(10)>2.0_ark*pi-sqrt(small_)) dsrc(10) = dsrc(10)-2.0_ark*pi
          !if (dsrc(10)<-sqrt(small_)) dsrc(10) = dsrc(10)+2.0_ark*pi
          !
          !dsrc(11) = src(11) 
          !if (dsrc(11)-rref(11,1)>2.0_ark*pi-sqrt(small_)) dsrc(11) = dsrc(11)-2.0_ark*pi
          !if (dsrc(11)-rref(11,1)<-sqrt(small_)) dsrc(11) = dsrc(11)+2.0_ark*pi
          !dsrc(12) = src(12) 
          !if (dsrc(12)-rref(12,1)>2.0_ark*pi-sqrt(small_)) dsrc(12) = dsrc(12)-2.0_ark*pi
          !if (dsrc(12)-rref(12,1)<-sqrt(small_)) dsrc(12) = dsrc(12)+2.0_ark*pi
          !
          tau = src(10) ! 1.0_ark/3.0_ark*( src(10)+src(11)+src(12) )
          !
          !tau = 1.0_ark/3.0_ark*( src(10)+src(11)+src(12) )
          !
          !tau = mod(tau+2.0_ark*pi,2.0_ark*pi/3.0_ark)
          !
          beta(1) = src(7)
          beta(2) = src(8)
          beta(3) = src(9)
          !
          phi1   = src(10)
          phi2   = src(11)
          phi3   = src(12)
          !
          cosa(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi2 - phi3 )
          cosa(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi3 - phi1 )
          cosa(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi1 - phi2 )

          !
          if ( any(abs(cosa(:))>1.0_ark+sqrt(small_))) then 
              !
              write (out,"('ML_coordinate_transform_ch3oh: cosalpha>1: ',3f18.8)") cosa(1:3)
              write (out,"('Consider change coordinate type ')")
              stop 'MLcoordinate_transform_func - bad cos(alpha)'
          endif 
          !
          cosa(:) = max(-1.0_ark,cosa(:))
          cosa(:) = min( 1.0_ark,cosa(:))
          !
          alpha(:) = acos(cosa(:))
          !
       else  ! not direct
          !
          tau = src(12)
          !
       endif
       !
       i0 = 0 
       !
       do icoord = 1,9
         !
         n(0:4) = molec%pot_ind(11,i0+1:i0+5)
         k(0:4) = molec%pot_ind(12,i0+1:i0+5)
         !
         if (all(n(:)==0)) then 
           !
           xi_eq(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*cos(real(k(1:4),ark)*tau))
           !
         else
           !
           xi_eq(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*sin(real(n(1:4),ark)*tau))
           !
         endif 
         !
         i0 = i0+5
         !
       enddo
       !
       xi_eq(10:11) = rref(10:11,0)
       xi_eq(12) = rref(12,0) + tau
       !
       beta(1) = (sqrt(2.0_ark)*xi_eq(7)+2.0_ark*xi_eq(8)                       )/sqrt(6.0_ark)
       beta(2) = (sqrt(2.0_ark)*xi_eq(7)-        xi_eq(8)+sqrt(3.0_ark)*xi_eq(9))/sqrt(6.0_ark)
       beta(3) = (sqrt(2.0_ark)*xi_eq(7)-        xi_eq(8)-sqrt(3.0_ark)*xi_eq(9))/sqrt(6.0_ark)
       !
       phi1   = (xi_eq(12)+2.0_ark*xi_eq(10)          )
       phi2   = (xi_eq(12)-        xi_eq(10)+xi_eq(11))
       phi3   = (xi_eq(12)-        xi_eq(10)-xi_eq(11))
       !
       cosa0(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi2 - phi3 )
       cosa0(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi3 - phi1 )
       cosa0(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi1 - phi2 )
       !
       if ( any(abs(cosa0(:))>1.0_ark+sqrt(small_))) then 
           !
           write (out,"('ML_coordinate_transform_ch3oh: cosalpha0>1: ',3f18.8)") cosa0(1:3)
           write (out,"('Consider change coordinate type ')")
           stop 'MLcoordinate_transform_func - bad cos(alpha)'
       endif 
       !
       cosa0(:) = max(-1.0_ark,cosa0(:))
       cosa0(:) = min( 1.0_ark,cosa0(:))
       !
       alpha0(:) = acos(cosa0(:))
       !
       xi_eq(10) = (2.0_ark*alpha0(1)-alpha0(2)-alpha0(3))/sqrt(6.0_ark)
       xi_eq(11) = (                  alpha0(2)-alpha0(3))/sqrt(2.0_ark)
       !
       if (direct) then
          !
          dst(1:2) = src(1:2)
          !
          dst( 3) = 1.0_ark/sqrt(3.0_ark)*(         src(3)+src(4)+src(5) )
          dst( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(3)-src(4)-src(5) )
          dst( 5) = 1.0_ark/sqrt(2.0_ark)*(                src(4)-src(5) )
          !
          dst( 6) = src(6)
          !
          dst( 7) = 1.0_ark/sqrt(3.0_ark)*(         src(7)+src(8)+src(9) )
          dst( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*src(7)-src(8)-src(9) )
          dst( 9) = 1.0_ark/sqrt(2.0_ark)*(                src(8)-src(9) )
          !
          !dsrc(10) = src(10)
          !
          !dsrc(11) = src(11) 
          !if (dsrc(11)-molec%local_eq(11)>2.0_ark*pi-sqrt(small_)) dsrc(11) = dsrc(11)-2.0_ark*pi
          !if (dsrc(11)-molec%local_eq(11)<-sqrt(small_)) dsrc(11) = dsrc(11)+2.0_ark*pi
          !dsrc(12) = src(12) 
          !if (dsrc(12)-molec%local_eq(12)>2.0_ark*pi-sqrt(small_)) dsrc(12) = dsrc(12)-2.0_ark*pi
          !if (dsrc(12)-molec%local_eq(12)<-sqrt(small_)) dsrc(12) = dsrc(12)+2.0_ark*pi
          !
          !dsrc(11) = mod(src(11),rref(11,0))+rref(11,0)
          !dsrc(12) = mod(src(12),rref(12,0))+rref(12,0)
          !
          !dst(7)   = dst(7)  /sqrt(2.0_ark)-(alpha(1) +alpha(2) +alpha(3) )/sqrt(6.0_ark)
          !xi_eq(7) = xi_eq(7)/sqrt(2.0_ark)-(alpha0(1)+alpha0(2)+alpha0(3))/sqrt(6.0_ark)
          !
          dst(10) = (2.0_ark*alpha(1)-alpha(2)-alpha(3))/sqrt(6.0_ark)
          dst(11) = (                alpha(2)-alpha(3))/sqrt(2.0_ark)
          dst(12) = tau
          !
          dst(1:11) = dst(1:11)-xi_eq(1:11)
          !
      else ! not direct
          !
          r_t(12) = tau
          !
          !xi_eq(7) = xi_eq(7)/sqrt(2.0_ark)-(alpha0(1)+alpha0(2)+alpha0(3))/sqrt(6.0_ark)
          !
          r_t(1:11) = src(1:11)+xi_eq(1:11)
          !
          dst(1:2) = r_t(1:2)
          !
          dst(3) =(sqrt(2.0_ark)*r_t(3)+2.0_ark*r_t(4)                     )/sqrt(6.0_ark)
          dst(4) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)+sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          dst(5) =(sqrt(2.0_ark)*r_t(3)-        r_t(4)-sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
          !
          dst(6) = r_t(6)
          !
          dst(7) =(sqrt(2.0_ark)*r_t(7)+2.0_ark*r_t(8)                     )/sqrt(6.0_ark)
          dst(8) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)+sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          dst(9) =(sqrt(2.0_ark)*r_t(7)-        r_t(8)-sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
          !
          call find_alpha_for_ch3oh(dst(7:9),r_t(10:12),dst(10:12))
          !call find_alpha_for_ch3oh_II(r_t(7:12),dst(7:12))
          !
          !dst(10) =(r_t(12)+2.0_ark*r_t(10)        )
          !dst(11) =(r_t(12)-        r_t(10)+r_t(11))
          !dst(12) =(r_t(12)-        r_t(10)-r_t(11))
          !
          !dst(10:12) = mod(dst(10:12)+2.0_ark*pi,2.0_ark*pi)+molec%local_eq(10:12)
          !
      endif
      !
    end select
    !
    !
    if (verbose>=7) write(out,"('ML_coordinate_transform_ch3oh/end')") 
    !
    !
  end function ML_coordinate_transform_ch3oh




  !
  subroutine find_alpha_from_beta(beta,xi,alpha2,alpha3,tau)
    ! obtaining 3 angles phi from 3 beta_i, tau, s_a, and s_b 
    real(ark),intent(in)  :: beta(3),xi(3)
    real(ark),intent(out) :: alpha2,alpha3,tau

    real(ark) :: eps(3)
    real(ark) :: rjacob(3,3),dx(3),phi_l(3),phi_r(3),xi_t(3),phi(3)

    real(ark) :: hstep
    real(rk)  :: stadev_old,stability,stadev,stadev_best,ssq,a(3,3),b(3,1)
    !
    integer(ik) :: iter,itmax,i1,irow,icolumn
         
    !
    rjacob = 0 
    iter = 0
    stadev_old = 1.e10
    stability =  1.e10
    stadev    =  1.e10
    !
    hstep = 1e-4
    !
    phi(1) = molec%local_eq( 9)
    phi(2) = molec%local_eq(11)
    phi(3) = molec%local_eq(12)
    xi_t = xi
    !xi_t(3) = mod(xi_t(3)+2.0_ark*pi,2.0_ark*pi)
    !
    stadev_best = sqrt(small_)*0.01_ark
    itmax = 500
    !
    outer_loop: & 
    do while( iter<itmax .and. stadev>stadev_best )   
       !
       iter = iter + 1
       ssq=0
       !
       ! Caclulate the function 
       !
       eps(:) = calc_xi(beta,phi)-xi_t(:)
       !
       ssq=sqrt(sum(eps(:)**2))
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       do i1 = 1,3
         !
         phi_r = phi ; phi_l = phi
         phi_r(i1) = phi(i1) + hstep
         phi_l(i1) = phi(i1) - hstep
         !
         rjacob(:,i1)  = ( calc_xi(beta,phi_r)-calc_xi(beta,phi_l))/hstep*0.5_ark
         !
       enddo
       !
       if (itmax>=0) then
         !
         ! form A matrix 
         do irow=1,3         !==== row-...... ====!
           do icolumn=1,irow      !==== column-....====!
               a(irow,icolumn)=sum( rjacob(:,irow)*rjacob(:,icolumn))
               a(icolumn,irow) = a(irow,icolumn)
           enddo
         enddo
         ! form B matrix 
         do irow=1,3       !==== row-...... ====!
           b(irow,1)=sum(eps(:)*rjacob(:,irow))
         enddo   
         !
         call lapack_gelss(a(:,:),b(:,:))
         !
         dx(:) = b(:,1)
         !
         stadev=sqrt(ssq)
         !
         phi(1:3) = phi(1:3) - dx(1:3)
         !
         !phi(:) = mod(phi(:)+2.0_ark*pi,2.0_ark*pi)
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
     write(out,"('find_alpha_from_beta: could not find solution after ',i8,' iterations')") iter
     stop 'find_alpha_from_beta: could not find solution'
    endif 
    !
    alpha2 = phi(1)
    alpha3 = phi(2)
    tau    = phi(3)
    !


  contains


  function calc_xi(beta,x) result (xi)

    real(ark),intent(in)  :: beta(3),x(3)
    real(ark)             :: xi(3)
    real(ark)             :: cosa(3),phi(3),tau,alpha(3),chi1,chi2,chi3
       


     alpha(2) = x(1)
     alpha(3) = x(2)
     tau = x(3)
     !
     cosa(2) =  ( cos(alpha(2))-cos(beta(3))*cos(beta(1)) )/( sin(beta(3))*sin(beta(1)) )
     cosa(3) =  ( cos(alpha(3))-cos(beta(2))*cos(beta(1)) )/( sin(beta(2))*sin(beta(1)) )
     !
     chi2 =-acos(cosa(2))
     chi3 = acos(cosa(3))
     !
     chi1 = 2.0_ark*pi-(-chi2+chi3)
     !
     cosa(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( chi1 )
     !
     alpha(1) = acos(cosa(1))
     !
     xi(1) = (2.0_ark*alpha(1)-alpha(2)-alpha(3))/sqrt(6.0_ark)
     xi(2) = (                 alpha(2)-alpha(3))/sqrt(2.0_ark)
     xi(3) = tau+(chi2+chi3)/3.0_ark
     !
     !
  end function calc_xi


  end subroutine find_alpha_from_beta
  !




  !
  subroutine find_phi_from_beta(beta,xi,phi)
    ! obtaining 3 angles phi from 3 beta_i, tau, s_a, and s_b 
    real(ark),intent(in)  :: beta(3),xi(3)
    real(ark),intent(out) :: phi(3)

    real(ark) :: eps(3)
    real(ark) :: rjacob(3,3),dx(3),phi_l(3),phi_r(3),xi_t(3)

    real(ark) :: hstep
    real(rk)  :: stadev_old,stability,stadev,stadev_best,ssq,a(3,3),b(3,1)
    !
    integer(ik) :: iter,itmax,i1,irow,icolumn
         
    !
    rjacob = 0 
    iter = 0
    stadev_old = 1.e10
    stability =  1.e10
    stadev    =  1.e10
    !
    hstep = 1e-4
    !
    phi(1) = molec%local_eq(10)
    phi(2) = molec%local_eq(11)-2.0_ark*pi
    phi(3) = molec%local_eq(12)
    xi_t = xi
    !xi_t(3) = mod(xi_t(3)+2.0_ark*pi,2.0_ark*pi)
    !
    stadev_best = sqrt(small_)*0.01_ark
    itmax = 500
    !
    outer_loop: & 
    do while( iter<itmax .and. stadev>stadev_best )   
       !
       iter = iter + 1
       ssq=0
       !
       ! Caclulate the function 
       !
       eps(:) = calc_xi(beta,phi)-xi_t(:)
       !
       ssq=(sum(eps(:)**2))
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       do i1 = 1,3
         !
         phi_r = phi ; phi_l = phi
         phi_r(i1) = phi(i1) + hstep
         phi_l(i1) = phi(i1) - hstep
         !
         rjacob(:,i1)  = ( calc_xi(beta,phi_r)-calc_xi(beta,phi_l))/hstep*0.5_ark
         !
       enddo
       !
       if (itmax>=0) then
         !
         ! form A matrix 
         do irow=1,3         !==== row-...... ====!
           do icolumn=1,irow      !==== column-....====!
               a(irow,icolumn)=sum( rjacob(:,irow)*rjacob(:,icolumn))
               a(icolumn,irow) = a(irow,icolumn)
           enddo
         enddo
         ! form B matrix 
         do irow=1,3       !==== row-...... ====!
           b(irow,1)=sum(eps(:)*rjacob(:,irow))
         enddo   
         !
         call lapack_gelss(a(:,:),b(:,:))
         !
         dx(:) = b(:,1)
         !
         stadev=sqrt(ssq/3.0_ark)
         !
         phi(1:3) = phi(1:3) - dx(1:3)
         !
         !phi(:) = mod(phi(:)+2.0_ark*pi,2.0_ark*pi)
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
     write(out,"('find_alpha_from_beta: could not find solution after ',i8,' iterations')") iter
     stop 'find_alpha_from_beta: could not find solution'
    endif 
    !
  contains


  function calc_xi(beta,x) result (xi)

    real(ark),intent(in)  :: beta(3),x(3)
    real(ark)             :: xi(3)
    real(ark)             :: cosa(3),phi(3),tau,alpha(3),chi(3),coschi(3)
       


     phi(:) = x(:)

     cosa(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi(2) - phi(3) )
     cosa(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi(3) - phi(1) )
     cosa(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi(1) - phi(2) )
     !
     if ( any(abs(cosa(:))>1.0_ark+sqrt(small_))) then 
         !
         write (out,"('ML_coordinate_transform_ch3oh: cosalpha>1: ',3f18.8)") cosa(1:3)
         write (out,"('Consider change coordinate type ')")
         stop 'MLcoordinate_transform_func - bad cos(alpha)'
     endif 
     !
     cosa(:) = max(-1.0_ark,cosa(:))
     cosa(:) = min( 1.0_ark,cosa(:))
     !
     alpha(:) = acos(cosa(:))
     !
     coschi(1) =  ( cos(alpha(1))-cos(beta(2))*cos(beta(3)) )/( sin(beta(2))*sin(beta(3)) )
     coschi(2) =  ( cos(alpha(2))-cos(beta(3))*cos(beta(1)) )/( sin(beta(3))*sin(beta(1)) )
     coschi(3) =  ( cos(alpha(3))-cos(beta(2))*cos(beta(1)) )/( sin(beta(2))*sin(beta(1)) )
     !
     chi(:) = acos(coschi(:))
     !
     xi(1) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*chi(1)-chi(2)-chi(3) )
     xi(2) = 1.0_ark/sqrt(2.0_ark)*(                chi(2)-chi(3) )
     xi(3) = 1.0_ark/3.0_ark*( phi(1)+phi(2)+phi(3) )

     !
  end function calc_xi


  end subroutine find_phi_from_beta



  !
  subroutine find_alpha_for_ch3oh(beta,xi,phi)
    ! obtaining 3 angles phi from 3 beta_i, tau, s_a, and s_b 
    real(ark),intent(in)  :: beta(3),xi(3)
    real(ark),intent(out) :: phi(3)

    real(ark) :: eps(3)
    real(ark) :: rjacob(3,3),dx(3),phi_l(3),phi_r(3),xi_t(3)

    real(ark) :: hstep
    real(rk)  :: stadev_old,stability,stadev,stadev_best,ssq,a(3,3),b(3,1)
    !
    integer(ik) :: iter,itmax,i1,irow,icolumn
         
    !
    rjacob = 0 
    iter = 0
    stadev_old = 1.e10
    stability =  1.e10
    stadev    =  1.e10
    !
    hstep = 1e-4
    !
    phi(1:3) = molec%local_eq(10:12)
    xi_t = xi
    !xi_t(3) = mod(xi_t(3)+2.0_ark*pi,2.0_ark*pi)
    !
    stadev_best = sqrt(small_)*10.0_ark
    itmax = 500
    !
    outer_loop: & 
    do while( iter<itmax .and. stadev>stadev_best )   
       !
       iter = iter + 1
       ssq=0
       !
       ! Caclulate the function 
       !
       eps(:) = calc_xi(beta,phi)-xi_t(:)
       !
       ssq=(sum(eps(:)**2))
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       do i1 = 1,3
         !
         phi_r = phi ; phi_l = phi
         phi_r(i1) = phi(i1) + hstep
         phi_l(i1) = phi(i1) - hstep
         !
         rjacob(:,i1)  = ( calc_xi(beta,phi_r)-calc_xi(beta,phi_l))/hstep*0.5_ark
         !
       enddo
       !
       if (itmax>=0) then
         !
         ! form A matrix 
         do irow=1,3         !==== row-...... ====!
           do icolumn=1,irow      !==== column-....====!
               a(irow,icolumn)=sum( rjacob(:,irow)*rjacob(:,icolumn))
               a(icolumn,irow) = a(irow,icolumn)
           enddo
         enddo
         ! form B matrix 
         do irow=1,3       !==== row-...... ====!
           b(irow,1)=sum(eps(:)*rjacob(:,irow))
         enddo   
         !
         call lapack_gelss(a(:,:),b(:,:))
         !
         dx(:) = b(:,1)
         !
         stadev=sqrt(ssq)
         !
         phi(1:3) = phi(1:3) - dx(1:3)
         !
         !phi(:) = mod(phi(:)+2.0_ark*pi,2.0_ark*pi)
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
     write(out,"('find_alpha_for_ch3oh: could not find solution after ',i8,' iterations')") iter
     stop 'find_alpha_for_ch3oh: could not find solution'
    endif 


  contains


  function calc_xi(beta,phi) result (xi)

    real(ark),intent(in)  :: beta(3),phi(3)
    real(ark)             :: xi(3)
    real(ark)             :: cosa(3),alpha(3)
       
    
      cosa(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi(2) - phi(3) )
      cosa(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi(3) - phi(1) )
      cosa(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi(1) - phi(2) )

      if ( any(abs(cosa(:))>1.0_ark+sqrt(small_))) then 
          !
          write (out,"('ML_coordinate_transform_ch3oh: cosalpha>1: ',3f18.8)") cosa(1:3)
          write (out,"('Consider change coordinate type ')")
          stop 'MLcoordinate_transform_func - bad cos(alpha)'
      endif 
      !
      cosa(:) = max(-1.0_ark,cosa(:))
      cosa(:) = min( 1.0_ark,cosa(:))
      !
      alpha(:) = acos(cosa(:))
      !
      xi(1) = (2.0_ark*alpha(1)-alpha(2)-alpha(3))/sqrt(6.0_ark)
      xi(2) = (                alpha(2)-alpha(3))/sqrt(2.0_ark)
      xi(3) = phi(1) ! 1.0_ark/3.0_ark*( phi(1)+phi(2)+phi(3) )
      !
      !xi(3) = mod(xi(3)+2.0_ark*pi,2.0_ark*pi/3.0_ark)
      !if (xi(3)>2.0_ark*pi/3.0_ark-small_) xi(3) = 2.0_ark*pi-xi(3)
      !
      !
  end function calc_xi


  end subroutine find_alpha_for_ch3oh
  !




  !
  subroutine find_alpha_for_ch3oh_II(xi,phi)
    ! obtaining 3 angles phi from 3 beta_i, tau, s_a, and s_b 
    real(ark),intent(in)  :: xi(6)
    real(ark),intent(out) :: phi(6)

    real(ark) :: eps(6)
    real(ark) :: rjacob(6,6),dx(6),phi_l(6),phi_r(6),xi_t(6)

    real(ark) :: hstep
    real(rk)  :: stadev_old,stability,stadev,stadev_best,ssq,a(6,6),b(6,1)
    !
    integer(ik) :: iter,itmax,i1,irow,icolumn
         
    !
    rjacob = 0 
    iter = 0
    stadev_old = 1.e10
    stability =  1.e10
    stadev    =  1.e10
    !
    hstep = 1e-4
    !
    phi(1:6) = molec%local_eq(7:12)
    xi_t(1:6) = xi(1:6)
    !xi_t(4)   = xi(1)
    !xi_t(3) = mod(xi_t(3)+2.0_ark*pi,2.0_ark*pi)
    !
    stadev_best = 1e-2 ! sqrt(small_)*10.0_ark
    itmax = 500
    !
    outer_loop: & 
    do while( iter<itmax .and. stadev>stadev_best )   
       !
       iter = iter + 1
       ssq=0
       !
       ! Caclulate the function 
       !
       eps(:) = calc_xi(phi)-xi_t(:)
       !
       ssq=sqrt(sum(eps(:)**2))
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       do i1 = 1,3
         !
         phi_r = phi ; phi_l = phi
         phi_r(i1) = phi(i1) + hstep
         phi_l(i1) = phi(i1) - hstep
         !
         rjacob(:,i1)  = ( calc_xi(phi_r)-calc_xi(phi_l))/hstep*0.5_ark
         !
       enddo
       !
       if (itmax>=0) then
         !
         ! form A matrix 
         do irow=1,3         !==== row-...... ====!
           do icolumn=1,irow      !==== column-....====!
               a(irow,icolumn)=sum( rjacob(:,irow)*rjacob(:,icolumn))
               a(icolumn,irow) = a(irow,icolumn)
           enddo
         enddo
         ! form B matrix 
         do irow=1,3       !==== row-...... ====!
           b(irow,1)=sum(eps(:)*rjacob(:,irow))
         enddo   
         !
         call lapack_gelss(a(:,:),b(:,:))
         !
         dx(:) = b(:,1)
         !
         stadev=sqrt(ssq)
         !
         phi(1:6) = phi(1:6) - dx(1:6)
         !
         !phi(:) = mod(phi(:)+2.0_ark*pi,2.0_ark*pi)
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
     write(out,"('find_alpha_for_ch3oh_II: could not find solution after ',i8,' iterations')") iter
     stop 'find_alpha_for_ch3oh_II: could not find solution'
    endif 



  contains


  function calc_xi(phi_) result (xi)

    real(ark),intent(in)  :: phi_(6)
    real(ark)             :: xi(6),beta(3),phi(3)
    real(ark)             :: cosa(3),alpha(3)
       
      beta(1:3) = phi_(1:3)
      phi(1:3)   = phi_(4:6)
      !
      cosa(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi(2) - phi(3) )
      cosa(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi(3) - phi(1) )
      cosa(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi(1) - phi(2) )

      if ( any(abs(cosa(:))>1.0_ark+sqrt(small_))) then 
          !
          write (out,"('ML_coordinate_transform_ch3oh: cosalpha>1: ',3f18.8)") cosa(1:3)
          write (out,"('Consider change coordinate type ')")
          stop 'MLcoordinate_transform_func - bad cos(alpha)'
      endif 
      !
      cosa(:) = max(-1.0_ark,cosa(:))
      cosa(:) = min( 1.0_ark,cosa(:))
      !
      alpha(:) = acos(cosa(:))
      !
      xi(1) = (beta(1)+beta(2)+beta(3)-(alpha(1)+alpha(2)+alpha(3)) )/sqrt(6.0_ark)
      xi(2) = (2.0_ark*beta(1)-beta(2)-beta(3) )/sqrt(6.0_ark)
      xi(3) = (                 beta(2)-beta(3) )/sqrt(2.0_ark)

      xi(4) = (2.0_ark*alpha(1)-alpha(2)-alpha(3))/sqrt(6.0_ark)
      xi(5) = (                alpha(2)-alpha(3))/sqrt(2.0_ark)
      xi(6) = phi(1) ! 1.0_ark/3.0_ark*( phi(1)+phi(2)+phi(3) )
      !
      !xi(3) = mod(xi(3)+2.0_ark*pi,2.0_ark*pi/3.0_ark)
      !if (xi(3)>2.0_ark*pi/3.0_ark-small_) xi(3) = 2.0_ark*pi-xi(3)
      !
      !
  end function calc_xi



  end subroutine find_alpha_for_ch3oh_II



  ! Here we define structural parameters a0 for ABCD molecule,
  !
  subroutine ML_b0_ch3oh(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
     !
     !
     integer(ik), intent(in) :: Npoints,Natoms
     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in) ,optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional  :: rho_ref
     real(ark),   intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
     !
     real(ark)    ::  transform(3,3),c(3,3),a0(molec%Natoms,3),CM_shift,rref(12,0:4),re1,re2,re3,ae1,ae2
     real(ark)    ::  r_t(1:12),Inert(3)
     real(ark)    :: r_at(1:12),a0_ark(molec%Natoms,3),beta1,beta2,beta3,phi1,phi2,phi3,alpha0(3),cosa0(3)
     real(ark)    :: tau,alpha1,alpha2,alpha3,tau1,tau2,tau3,chi1,chi2,chi3,r_eq(1:12)

     integer(ik) ::  k(0:4),n(0:4),i,i0,icoord

      if (verbose>=4) write(out,"('ML_b0_ch3oh/start')") 
      !
      if (size(molec%req)/=5) then
        write(out,"('ML_b0_ch3oh: Nbonds must be 5 in this routine, not  ',i9)") size(molec%req)
        stop 'ML_b0_ch3oh: wrong Nbonds '
      endif 
      !
      if (molec%Natoms/=6) then
        write(out,"('ML_b0_ch3oh: Natoms must be 6 in this routine, not  ',i9)") molec%Natoms
        stop 'ML_b0_ch3oh: wrong Natoms '
      endif 
      !
      !
      tau = molec%taueq(1)
      !
      r_at = molec%local_eq
      !
      !
      !call ML_CH3OH_MEP_geometry(tau,r_eq(1),r_eq(2),r_eq(6),r_eq(3),r_eq(4),r_eq(5),&
      !                                         beta1,beta2,beta3,alpha1,alpha2,alpha3,tau1,tau2,tau3,chi1,chi2,chi3)
      !
      select case(trim(molec%frame))
      case default
         write (out,"('ML_b0_ch3oh: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'ML_b0_ch3oh - bad coord. type'
         !
      case('R-TAU-BOWMAN-REF')
         !
         r_at(1:6) = r_eq(1:6)
         !
         r_at(7)  = beta1
         r_at(8)  = beta2
         r_at(9)  = alpha2
         r_at(10) = beta3
         r_at(11) = alpha3
         r_at(12) = tau
         !
      case('R-TAU-CHI')
         !
         r_at(1:6) = r_eq(1:6)
         !
         r_at(7)  = beta1
         r_at(8)  = beta2
         r_at(9)  = beta3
         r_at(10) = tau1
         r_at(11) = tau2
         r_at(12) = tau3
         !
      case('R-ALPHA-THETA-TAU','R-ALPHA-THETA--TAU')
         !
         r_at = molec%local_eq
         !
      end select 
      !
      call MLfromlocal2cartesian(-1,r_at,a0_ark)
      !
      a0 = a0_ark
      !
      call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,transform)
      !
      b0(:,:,0) = a0(:,:)
      !
      ! We can also simulate the reference structure at each point of rho, if npoints presented
      !
      if (Npoints/=0) then
         ! 
         if (.not.present(rho_borders).or..not.present(rho_ref)) then  
            write(out,"('ML_b0_ch3oh: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_ch3oh: rho_borders or rho_ref not specified '
         endif
         !
         rho_ref = pi
         !
         r_eq = molec%local_eq
         !
         do i = 0,npoints
            !
            tau = rho_i(i)
            !
            !
            !call ML_CH3OH_MEP_geometry(tau,r_eq(1),r_eq(2),r_eq(6),r_eq(3),r_eq(4),r_eq(5),&
            !                                       beta1,beta2,beta3,alpha1,alpha2,alpha3,tau1,tau2,tau3,chi1,chi2,chi3)
            !
            select case(trim(molec%frame))
            case default
               write (out,"('ML_b0_ch3oh: coord. type ',a,' unknown')") trim(molec%coords_transform)
               stop 'ML_b0_ch3oh - bad coord. type'
               !
            case('R-TAU-BOWMAN-REF')

            !if (trim(molec%potentype)=='POTEN_CH3OH_REF'.or.trim(molec%potentype)=='R-ALPHA-S-TAU-BETA-REF'&
            !.or.trim(molec%potentype)=='R-TAU-BOWMAN-REF'.or.trim(molec%potentype)=='R-1TAU-2BETA-REF') then 

               !
               r_at(1:6) = r_eq(1:6)
               !
               r_at(7)  = beta1
               r_at(8)  = beta2
               r_at(9)  = alpha2
               r_at(10) = beta3
               r_at(11) = alpha3
               r_at(12) = tau
               !
            case('R-TAU-CHI')
               !
               r_at(1:6) = r_eq(1:6)
               !
               r_at(7)  = beta1
               r_at(8)  = beta2
               r_at(9)  = beta3
               r_at(10) = tau1
               r_at(11) = tau2
               r_at(12) = tau3
               !
            case('R-ALPHA-THETA-TAU')
               !
               r_eq(10) = tau
               r_eq(11) = tau+molec%local_eq(10)+molec%local_eq(11)
               r_eq(12) = tau+molec%local_eq(10)+molec%local_eq(12)
               !
               r_at = r_eq
               !
            case('R-ALPHA-THETA--TAU')
               !
               r_eq(10) = -tau+molec%local_eq(10)
               r_eq(11) = -tau+molec%local_eq(10)+molec%local_eq(11)
               r_eq(12) = -tau+molec%local_eq(10)+molec%local_eq(12)
               !
               r_at = r_eq
               !
            end select
            !
            if (verbose>=6) then 
              write(out,"('r0',i8,12f12.8)") i,r_at(1:12)
            endif 
            !
            call MLfromlocal2cartesian(-1,r_at,a0_ark)
            !
            b0(:,:,i) = a0_ark(:,:)
            !
            call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i),transform)
            !
            Inert(1) = sum(molec%AtomMasses(:)*( a0_ark(:,2)**2+ a0_ark(:,3)**2) )
            Inert(2) = sum(molec%AtomMasses(:)*( a0_ark(:,1)**2+ a0_ark(:,3)**2) )
            Inert(3) = sum(molec%AtomMasses(:)*( a0_ark(:,1)**2+ a0_ark(:,2)**2) )
            !
            if (verbose>=6) then 
              write(out,"('i0',3f12.8)") Inert(1:3)
            endif
            !
            !if (verbose>=5) then 
            !  write(out,"('b0',i8,18f12.8)") i,b0(:,:,i)
            !endif

         enddo
         !
         if (verbose>=4) then 
           do i = 0,npoints
              write(out,"(i6)") molec%natoms
              !
              write(out,"(/a,3x,3f14.8)") trim(molec%zmatrix(1)%name),b0(1,:,i) 
              write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(2)%name),b0(2,:,i)
              write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(3)%name),b0(3,:,i)
              write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(3)%name),b0(4,:,i)
              write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(3)%name),b0(5,:,i)
              write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(3)%name),b0(6,:,i)
              !
           enddo
         endif
         !
         do i = 0,npoints
            if (verbose>=6) then 
              write(out,"('b0',i8,18f12.8)") i,( (b0(i0,icoord,i),icoord=1,3),i0=1,6 )
            endif
         enddo 
         !
      endif
      !
      !
      if (verbose>=4) write(out,"('ML_b0_ch3oh/end')") 
      !
  end subroutine ML_b0_ch3oh


  subroutine ML_symmetry_transformation_CH3OH(ioper, nmodes, src, dst)
    !
    ! Symmetry transformation rules of coordinates
    !
    integer(ik), intent(in)  ::  ioper
    integer(ik), intent(in)  ::  nmodes
    real(ark), intent(in)    ::  src(1:nmodes)
    real(ark), intent(out)   ::  dst(1:nmodes)
    real(ark) :: a,b,p
    !
    a = 0.5_ark
    b = 0.5_ark*sqrt(3.0_ark)
    p = 2.0_ark*pi/3.0_ark
    if (verbose>=6) write(out, '(/a)') 'ML_symmetry_transformation_CH3OH/start'
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_symmetry_transformation_CH3OH error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_symmetry_transformation_CH3OH error: bad coordinate type'
      !
    case('R-ALPHA-THETA-TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_symmetry_transformation_CH3OH error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop
        !
      case('C3V(M)-working')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_CH3OH error: symmetry operation ', ioper, 'is unknown'
          stop
          !
        case (1) ! E
          !
          dst(1:12) = src(1:12)
          !
        case (2) ! (123)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(4)
          dst(4) = src(5)
          dst(5) = src(3)
          !
          dst(6) = src(6)
          !
          dst(7) = src(8)
          dst(8) = src(9)
          dst(9) = src(7)
          !
          dst(10) = -a*src(10) + b*src(11)
          dst(11) = -b*src(10) - a*src(11)
          !
          dst(12) = mod(src(12) + p,2.0_ark*pi)
          !
        case (3) !(132)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(3)
          dst(5) = src(4)
          !
          dst(6) = src(6)
          !
          dst(7) = src(9)
          dst(8) = src(7)
          dst(9) = src(8)
          !
          dst(10) = -a*src(10) - b*src(11)
          dst(11) = +b*src(10) - a*src(11)
          !
          dst(12) = mod(src(12) + 2.0_ark*p,2.0_ark*pi)
          !
        case (4) ! (32)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(3)
          dst(4) = src(5)
          dst(5) = src(4)
          !
          dst(6) = src(6)
          !
          dst(7) = src(7)
          dst(8) = src(9)
          dst(9) = src(8)
          !
          dst(10) = src(10)
          dst(11) = -src(11)
          dst(12) = -src(12)+2.0_ark*pi
          !
        case (5) ! (12)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(5)
          !
          dst(6) = src(6)
          !
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(9)
          !
          dst(10) = -a*src(10) + b*src(11)
          dst(11) = +b*src(10) + a*src(11)
          dst(12) = 2.0_ark*pi-mod(src(12)+p,2.0_ark*pi)
          !
        case (6) ! (13)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(4)
          dst(5) = src(3)
          !
          dst(6) = src(6)
          !
          dst(7) = src(9)
          dst(8) = src(8)
          dst(9) = src(7)
          !
          dst(10) = -a*src(10) - b*src(11)
          dst(11) = -b*src(10) + a*src(11)
          !
          dst(12) = 2.0_ark*pi-mod(src(12)+2.0_ark*p,2.0_ark*pi)
          !
        end select
        !
      case('C3V(M)')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_CH3OH error: symmetry operation ', ioper, 'is unknown'
          stop
          !
        case (1) ! E
          !
          dst(1:12) = src(1:12)
          !
        case (3) ! (123)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(4)
          dst(4) = src(5)
          dst(5) = src(3)
          !
          dst(6) = src(6)
          !
          dst(7) = src(8)
          dst(8) = src(9)
          dst(9) = src(7)
          !
          dst(10) = -a*src(10) + b*src(11)
          dst(11) = -b*src(10) - a*src(11)
          !
          dst(12) = mod(src(12) + p,2.0_ark*pi)
          !
        case (2) !(132)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(3)
          dst(5) = src(4)
          !
          dst(6) = src(6)
          !
          dst(7) = src(9)
          dst(8) = src(7)
          dst(9) = src(8)
          !
          dst(10) = -a*src(10) - b*src(11)
          dst(11) = +b*src(10) - a*src(11)
          !
          dst(12) = mod(src(12) + 2.0_ark*p,2.0_ark*pi)
          !
        case (4) ! (32)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(3)
          dst(4) = src(5)
          dst(5) = src(4)
          !
          dst(6) = src(6)
          !
          dst(7) = src(7)
          dst(8) = src(9)
          dst(9) = src(8)
          !
          dst(10) = src(10)
          dst(11) = -src(11)
          dst(12) = -src(12)+2.0_ark*pi
          !
        case (6) ! (12)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(5)
          !
          dst(6) = src(6)
          !
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(9)
          !
          dst(10) = -a*src(10) + b*src(11)
          dst(11) = +b*src(10) + a*src(11)
          dst(12) = 2.0_ark*pi-mod(src(12)+2.0_ark*p,2.0_ark*pi)
          !
        case (5) ! (13)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(4)
          dst(5) = src(3)
          !
          dst(6) = src(6)
          !
          dst(7) = src(9)
          dst(8) = src(8)
          dst(9) = src(7)
          !
          dst(10) = -a*src(10) - b*src(11)
          dst(11) = -b*src(10) + a*src(11)
          !
          dst(12) = 2.0_ark*pi-mod(src(12)+p,2.0_ark*pi)
          !
        end select
        !
      case('C3V(M)-2')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_CH3OH error: symmetry operation ', ioper, 'is unknown'
          stop
          !
        case (1) ! E
          !
          dst(1:12) = src(1:12)
          !
        case (3) ! (123)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(4)
          dst(4) = src(5)
          dst(5) = src(3)
          !
          dst(6) = src(6)
          !
          dst(7) = src(8)
          dst(8) = src(9)
          dst(9) = src(7)
          !
          dst(10) = -a*src(10) + b*src(11)
          dst(11) = -b*src(10) - a*src(11)
          !
          dst(12) = mod(src(12) + p,2.0_ark*pi)
          !
        case (2) !(132)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(3)
          dst(5) = src(4)
          !
          dst(6) = src(6)
          !
          dst(7) = src(9)
          dst(8) = src(7)
          dst(9) = src(8)
          !
          dst(10) = -a*src(10) - b*src(11)
          dst(11) = +b*src(10) - a*src(11)
          !
          dst(12) = mod(src(12) + 2.0_ark*p,2.0_ark*pi)
          !
        case (4) ! (32)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(3)
          dst(4) = src(5)
          dst(5) = src(4)
          !
          dst(6) = src(6)
          !
          dst(7) = src(7)
          dst(8) = src(9)
          dst(9) = src(8)
          !
          dst(10) = src(10)
          dst(11) = -src(11)
          dst(12) = -src(12)+2.0_ark*pi
          !
        case (6) ! (12)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(5)
          !
          dst(6) = src(6)
          !
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(9)
          !
          dst(10) = -a*src(10) + b*src(11)
          dst(11) = +b*src(10) + a*src(11)
          dst(12) = 2.0_ark*pi-mod(src(12)+p,2.0_ark*pi)
          !
        case (5) ! (13)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(4)
          dst(5) = src(3)
          !
          dst(6) = src(6)
          !
          dst(7) = src(9)
          dst(8) = src(8)
          dst(9) = src(7)
          !
          dst(10) = -a*src(10) - b*src(11)
          dst(11) = -b*src(10) + a*src(11)
          !
          dst(12) = 2.0_ark*pi-mod(src(12)-p,2.0_ark*pi)
          !
        end select
        !
      end select
      !
    end select
    !
    if (verbose>=6) write(out, '(/a)') 'ML_symmetry_transformation_CH3OH/end'
    !
  end subroutine ML_symmetry_transformation_CH3OH


  subroutine ML_rotsymmetry_CH3OH(J,K,tau,gamma,ideg)
    !
    ! Symmetry transformation rules of TROVE rotational functions
    !
    integer(ik),intent(in) :: J,K,tau
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_CH3OH/start'
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_rotsymmetry_CH3OH error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_rotsymmetry_CH3OH error: bad coordinate type'
      !
    case('R-ALPHA-THETA-TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_rotsymmetry_CH3OH: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_CH3OH error: bad symmetry type'
        !
      case('C','C(M)')
        !
        gamma = 1
        ideg = 1
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
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_CH3OH/end'
    !
  end subroutine ML_rotsymmetry_CH3OH


end module mol_ch3oh

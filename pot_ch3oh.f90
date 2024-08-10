!
!  This unit defines all specific routines for a fouratomic molecule of ch3oh type
!
module pot_ch3oh
  use accuracy
  use moltype
  use lapack
  !
  implicit none
  !
  public MLpoten_ch3oh_ref,ML_CH3OH_MEP_geometry,MLpoten_ch3oh_sym
  !
  private
  !
  integer(ik), parameter :: verbose     = 5                       ! Verbosity level
  !
  contains

  subroutine ML_CH3OH_MEP_geometry(tau,rco,roh,betacoh,r1,r2,r3,beta1,beta2,beta3,alpha1,alpha2,alpha3,& 
                                                  tau1,tau2,tau3,chi1,chi2,chi3)

     real(ark),intent(in) :: tau
     real(ark),intent(out) :: rco,roh,betacoh,r1,r2,r3,beta1,beta2,beta3,alpha1,alpha2,alpha3,tau1,tau2,tau3,chi1,chi2,chi3
     integer(ik)  :: i0,n(0:5),k(0:5),icoord
     real(ark)    :: rref(12,0:4),r_t(12),r_eq(12)
     real(ark)    :: cosa0(3),coschi(3)
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
     i0 = 0 
     !
     do icoord = 1,9
       !
       n(0:4) = molec%pot_ind(11,i0+1:i0+5)
       k(0:4) = molec%pot_ind(12,i0+1:i0+5)
       !
       if (all(n(:)==0)) then 
         !
         r_t(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*cos(real(k(1:4),ark)*tau))
         !
       else
         !
         r_t(icoord) = rref(icoord,0)+sum(rref(icoord,1:4)*sin(real(n(1:4),ark)*tau))
         !
       endif 
       !
       i0 = i0+5
       !
     enddo
     !
     r_t(10:11) = rref(10:11,0)
     !
     r_t(12) = tau !+molec%local_eq(12)
     !
     rco = r_t(1)
     roh = r_t(2)
     !
     r1 =(sqrt(2.0_ark)*r_t(3)+2.0_ark*r_t(4)                     )/sqrt(6.0_ark)
     r2 =(sqrt(2.0_ark)*r_t(3)-        r_t(4)+sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
     r3 =(sqrt(2.0_ark)*r_t(3)-        r_t(4)-sqrt(3.0_ark)*r_t(5))/sqrt(6.0_ark)
     !
     betacoh = r_t(6)
     !
     beta1 =(sqrt(2.0_ark)*r_t(7)+2.0_ark*r_t(8)                     )/sqrt(6.0_ark)
     beta2 =(sqrt(2.0_ark)*r_t(7)-        r_t(8)+sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
     beta3 =(sqrt(2.0_ark)*r_t(7)-        r_t(8)-sqrt(3.0_ark)*r_t(9))/sqrt(6.0_ark)
     !
     tau1   = (r_t(12)+2.0_ark*r_t(10)        )
     tau2   = (r_t(12)-        r_t(10)+r_t(11))
     tau3   = (r_t(12)-        r_t(10)-r_t(11))
     !
     cosa0(1) = cos( beta2 )*cos( beta3 ) + sin( beta2 )*sin( beta3 )*cos( tau2 - tau3 )
     cosa0(2) = cos( beta3 )*cos( beta1 ) + sin( beta3 )*sin( beta1 )*cos( tau3 - tau1 )
     cosa0(3) = cos( beta1 )*cos( beta2 ) + sin( beta1 )*sin( beta2 )*cos( tau1 - tau2 )
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
     alpha1 = acos(cosa0(1))
     alpha2 = acos(cosa0(2))
     alpha3 = acos(cosa0(3))
     !
     coschi(1) =  ( cos(alpha1)-cos(beta2)*cos(beta3) )/( sin(beta2)*sin(beta3) )
     coschi(2) =  ( cos(alpha2)-cos(beta3)*cos(beta1) )/( sin(beta3)*sin(beta1) )
     coschi(3) =  ( cos(alpha3)-cos(beta2)*cos(beta1) )/( sin(beta2)*sin(beta1) )
     !
     chi1 = acos(coschi(1))
     chi2 = acos(coschi(2))
     chi3 = acos(coschi(3))
     !
  end subroutine ML_CH3OH_MEP_geometry
  !
  !
  ! Procedures to define  ch3oh 
  !
  function mlpoten_ch3oh_ref(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          ::  i,k(12),i0,n(0:4),m(0:4),icoord
   real(ark)    :: xi(12),y(12),a0(5),tau,rref(12,0:4),r_t(12),xi_t(12),r_eq(12)
   real(ark)     :: beta(3),alpha0(3),cosa0(3),phi1,phi2,phi3,alpha(3),chi1,chi2,chi3,cosa(3)
   real(ark)     :: beta0(3),phi0(3),chi0(3),coschi(3),chi(3)
   !
   if (verbose>=6) write(out,"('MLpoten_ch3oh_ref/start')") 
   !
   select case(trim(molec%coords_transform))
   case default
      write (out,"('MLcoordinate_transform_func: coord. type ',a,' unknown')") trim(molec%coords_transform)
      stop 'mlpoten_ch3oh_ref - bad coord. type'
      !
   case('R-TAU-BOWMAN-REF','R-ALPHA-THETA-TAU')
        !
        alpha(3) = local( 9)
        alpha(2) = local(11)
        !
        beta(1) = local(7)
        beta(2) = local(8)
        beta(3) = local(10)
        !
        cosa(2) =  ( cos(alpha(2))-cos(beta(3))*cos(beta(1)) )/( sin(beta(3))*sin(beta(1)) )
        cosa(3) =  ( cos(alpha(3))-cos(beta(2))*cos(beta(1)) )/( sin(beta(2))*sin(beta(1)) )
        !
        chi2 =-acos(cosa(2))
        chi3 = acos(cosa(3))
        !
        chi1 = 2.0_ark*pi-(-chi2+chi3)
        !
        tau = local(12)+(chi2+chi3)/3.0_ark
        !
   case('R-TAU-CHI')
        !
        phi1   = local(10)
        phi2   = local(11)
        phi3   = local(12)
        !
        if (phi2>2.0_ark/3.0_ark*pi)  phi2 = phi2-2.0_ark*pi
        !
        tau = 1.0_ark/3.0_ark*( phi1+phi2+phi3 )
        !
   end select 
   !
   call ML_CH3OH_MEP_geometry(tau,r_eq(1),r_eq(2),r_eq(6),r_eq(3),r_eq(4),r_eq(5),beta0(1),beta0(2),beta0(3),&
                              alpha0(1),alpha0(2),alpha0(3),phi0(1),phi0(2),phi0(3),chi0(1),chi0(2),chi0(3))


   !
   select case(trim(molec%coords_transform))
   case default
      write (out,"('mlpoten_ch3oh_ref: coord. type ',a,' unknown')") trim(molec%coords_transform)
      stop 'mlpoten_ch3oh_ref - bad coord. type'
      !
   case('R-TAU-CHI')
      !
      r_t(1:2) = r_eq(1:2)
      !
      r_t( 3) = 1.0_ark/sqrt(3.0_ark)*(         r_eq(3)+r_eq(4)+r_eq(5) )
      r_t( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*r_eq(3)-r_eq(4)-r_eq(5) )
      r_t( 5) = 1.0_ark/sqrt(2.0_ark)*(                 r_eq(4)-r_eq(5) )
      !
      r_t(6) = r_eq(6)
      !
      r_t( 7) = 1.0_ark/sqrt(3.0_ark)*(         beta0(1)+beta0(2)+beta0(3) )
      r_t( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*beta0(1)-beta0(2)-beta0(3) )
      r_t( 9) = 1.0_ark/sqrt(2.0_ark)*(                  beta0(2)-beta0(3) )
      !
      r_t(10) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*chi0(1)-chi0(2)-chi0(3) ) 
      r_t(11) = 1.0_ark/sqrt(2.0_ark)*(                 chi0(2)-chi0(3) ) 
      r_t(12) = tau
      !
      y(1:2) = local(1:2)
      !
      y( 3) = 1.0_ark/sqrt(3.0_ark)*(         local(3)+local(4)+local(5) )
      y( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*local(3)-local(4)-local(5) )
      y( 5) = 1.0_ark/sqrt(2.0_ark)*(                  local(4)-local(5) )
      !
      y( 6) = local(6)
      !
      y( 7) = 1.0_ark/sqrt(3.0_ark)*(         local(7)+local(8)+local(9) )
      y( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*local(7)-local(8)-local(9) )
      y( 9) = 1.0_ark/sqrt(2.0_ark)*(                  local(8)-local(9) )
      !
      beta(1) = local(7)
      beta(2) = local(8)
      beta(3) = local(9)
      !
      cosa(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( phi2 - phi3 )
      cosa(2) = cos( beta(3) )*cos( beta(1) ) + sin( beta(3) )*sin( beta(1) )*cos( phi3 - phi1 )
      cosa(3) = cos( beta(1) )*cos( beta(2) ) + sin( beta(1) )*sin( beta(2) )*cos( phi1 - phi2 )
      !
      if ( any(abs(cosa(:))>1.0_ark+sqrt(small_))) then 
          !
          write (out,"('ML_coordinate_transform_ch3oh: cosalpha>1: ',3f18.8)") cosa(1:3)
          write (out,"('Consider change coordinate type ')")
          stop 'mlpoten_ch3oh_ref - bad cos(alpha)'
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
      y(10) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*chi(1)-chi(2)-chi(3) )
      y(11) = 1.0_ark/sqrt(2.0_ark)*(                chi(2)-chi(3) )
      y(12) = tau
      !
      y(1:11) = y(1:11)-r_t(1:11)
      !
   case('R-TAU-BOWMAN-REF','R-ALPHA-THETA-TAU')
      !
      r_t(1:2) = r_eq(1:2)

      r_t( 3) = 1.0_ark/sqrt(3.0_ark)*(         r_eq(3)+r_eq(4)+r_eq(5) )
      r_t( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*r_eq(3)-r_eq(4)-r_eq(5) )
      r_t( 5) = 1.0_ark/sqrt(2.0_ark)*(                 r_eq(4)-r_eq(5) )
      !
      r_t(6) = r_eq(6)
      !
      r_t( 7) = 1.0_ark/sqrt(3.0_ark)*(         beta0(1)+beta0(2)+beta0(3) )
      r_t( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*beta0(1)-beta0(2)-beta0(3) )
      r_t( 9) = 1.0_ark/sqrt(2.0_ark)*(                  beta0(2)-beta0(3) )
      !
      r_t(10) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha0(1)-alpha0(2)-alpha0(3) )
      r_t(11) = 1.0_ark/sqrt(2.0_ark)*(                   alpha0(2)-alpha0(3) )
      r_t(12) = tau
      !
      y(1:2) = local(1:2)
      !
      y( 3) = 1.0_ark/sqrt(3.0_ark)*(         local(3)+local(4)+local(5) )
      y( 4) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*local(3)-local(4)-local(5) )
      y( 5) = 1.0_ark/sqrt(2.0_ark)*(                  local(4)-local(5) )
      !
      y( 6) = local(6)
      !
      y( 7) = 1.0_ark/sqrt(3.0_ark)*(         local(7)+local(8)+local(9) )
      y( 8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*local(7)-local(8)-local(9) )
      y( 9) = 1.0_ark/sqrt(2.0_ark)*(                  local(8)-local(9) )
      !
      beta(1) = local(7)
      beta(2) = local(8)
      beta(3) = local(9)
      !
      chi1 = 2.0_ark*pi+chi2-chi3
      !
      cosa0(1) = cos( beta(2) )*cos( beta(3) ) + sin( beta(2) )*sin( beta(3) )*cos( chi1 )
      !
      alpha(1) = acos(cosa0(1))
      !
      y(10) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*alpha(1)-alpha(2)-alpha(3) )
      y(11) = 1.0_ark/sqrt(2.0_ark)*(                  alpha(2)-alpha(3) )
      y(12) = tau
      !
      y(1:11) = y(1:11)-r_t(1:11)
      !
   end select 
   !
   a0(1:5) = molec%specparam(1:5)
   !
   y(12) = tau
   !
   f = 0
   i0 = 61
   do i = i0,molec%parmax
      !
      k(:) = molec%pot_ind(:,i)
      !
      xi(1:11) = y(1:11)**k(1:11)
      xi(12) = cos(real(k(12),ark)*tau)
      !
      !do i1 = 1,12
      !   xi(i1) =1.0_ark ; if (k(i1)/=0) xi(i1) = y(i1)**k(i1)
      !enddo
      !
      f = f + force(i)*product(xi(:))
      !
   enddo
   !
   f = f*1.0e-11/planck/vellgt
   !
   if (verbose>=6) write(out,"('MLpoten_ch3oh_ref/end')") 
 
 end function MLpoten_ch3oh_ref
 
 
function MLpoten_ch3oh_sym(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(12),y(12)
  real(ark) :: chi(12,6),term,deg
  integer(ik) :: ioper,ipower(12),i,nparam_eq,iparam,parmax,itau
  integer,parameter :: Nsym = 6, Ndeg = 12
  logical :: do_cos
  !
  real(ark) :: rCO,ROH,RCOe,ROHe,alpha_COHe,alpha_OCHe,deltae,RH1,RH2,RH3,RH1E
  real(ark) :: b1,b3,b4,b5,theta1,theta2,theta3,theta12,theta23,theta13,a1,a2,tau
  !
  nparam_eq = 8
  deg = pi/180.0_ark
  !
  parmax = size(force)
  !
  f = 0
  !
  select case(trim(molec%coords_transform))
    !
  case default
    !
    write(out, '(/a,1x,a,1x,a)') &
    'MLpoten_ch3oh_sym error', trim(molec%coords_transform), 'is unknown'
    stop 'MLpoten_ch3oh_sym error error: bad coordinate type'
      !
  case('R-ALPHA-THETA-TAU')
    !
    rCOe      = force(1)
    rOHe      = force(2)
    rH1e      = force(3)
    !
    alpha_COHe = force(4)*deg
    alpha_OCHe = force(5)*deg
    deltae = 2.0_ark*pi/3.0_ark
    !
    a1       = force(6)
    a2       = force(7)
    b1       = force(8)
    !
    rCO      = local(1)
    rOH      = local(2)
    rH1      = local(3)
    rH2      = local(4)
    rH3      = local(5)
    !
    ! C-O
    xi(1)=1.0_ark-exp(-a1*(rCO-rCOe))
    xi(2)=1.0_ark-exp(-a2*(rOH-rOHe))
    !CH3
    xi(3)=1.0_ark-exp(-b1*(rH1-rH1e))
    xi(4)=1.0_ark-exp(-b1*(rH2-rH1e))
    xi(5)=1.0_ark-exp(-b1*(rH3-rH1e))
    !
    !OCH
    xi(6) = local(6)- alpha_COHe
    !
    ! alphas
    !H2-C-H3,H1-C-H2,H1-C-H3
    xi(7) = local(7)- alpha_OCHe
    xi(8) = local(8)- alpha_OCHe
    xi(9) = local(9)- alpha_OCHe
    !
    theta1 = local(10)
    theta1 = mod(theta1+2.0_ark*pi,2.0_ark*pi)
    !
    theta2 = local(11)
    theta3 = local(12)
    !
    if (theta2-theta1<small_) theta2 = theta2 + 2.0_ark*pi
    if (theta3-theta2<small_) theta3 = theta3 + 2.0_ark*pi
    !
    theta12 = mod(theta2-theta1+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(theta3-theta2+2.0_ark*pi,2.0_ark*pi)
    !theta13 = mod(theta1-theta3+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    xi(10)   = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
    xi(11)   = (                   theta13 - theta12 )/sqrt(2.0_ark)
    !
    if (theta1<-sqrt(small_)) then
      theta1 = theta1 + 2.0_ark*pi
    endif
    if (theta1>2.0_ark*pi+sqrt(small_)) then
      theta1 = theta1 - 2.0_ark*pi
    endif
    !
    tau = (theta1 + theta2 + theta3-2.0_ark*pi)/3.0_ark
    !
    xi(12) = tau
    !
    call ML_symmetry_transformation_CH3OH_II(nsym,xi,chi,ndeg)
    !
    do iparam = nparam_eq+1,parmax
      !
      do_cos = .true.
      !      
      itau = molec%pot_ind(12,iparam)
      !    
      if (itau<0) then
        do_cos = .false.
        itau = abs(itau)
      endif
      !
      term = 0
      !
      do ioper =1,Nsym
        !
        y(1:11) = chi(1:11,ioper)**molec%pot_ind(1:11,iparam)
        if (do_cos) then
          !print*,itau,chi(21,ioper)
          y(12) = cos(real(itau,8)*chi(12,ioper))
        else
          y(12) = sin(real(itau,8)*chi(12,ioper))
        endif
        !
        term = term + product(y(:))
        !
      enddo
      !
      f = f + term*force(iparam)/6.0_ark
      !
    enddo   
    !
  end select 
  !  
end function  MLpoten_ch3oh_sym
 
 
 
 subroutine ML_symmetry_transformation_CH3OH(nsym,src,dst,ndeg)
    implicit none
    !
    integer(ik),intent(in)    :: nsym,ndeg  ! group operation  
    real(ark),intent(in)      :: src(1:ndeg)
    real(ark),intent(out)     :: dst(1:ndeg,6)
    !
    real(ark)         :: repres(nsym,ndeg,ndeg),a,b,e,o
    integer(rk)       :: ioper
    !
    if (verbose>=6) write(out,"('ML_symmetry_transformation_CH3OH/start')")
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    !
    repres(:,12,12) = 1.0_ark
    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !d9
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    !
    !C3+/(132)
    repres(2,3,5) = 1.0_ark
    repres(2,4,3) = 1.0_ark
    repres(2,5,4) = 1.0_ark
    !
    repres(2,7,9) = 1.0_ark
    repres(2,8,7) = 1.0_ark
    repres(2,9,8) = 1.0_ark
    !
    repres(2,10,10) = -a
    repres(2,10,11) = -b
    repres(2,11,10) =  b
    repres(2,11,11) = -a
    !
    !C3-/(93)
    !
    repres(3,3,4) = 1.0_ark
    repres(3,4,5) = 1.0_ark
    repres(3,5,3) = 1.0_ark
    !
    repres(3,7,8) = 1.0_ark
    repres(3,8,9) = 1.0_ark
    repres(3,9,7) = 1.0_ark
    !
    repres(3,10,10) = -a
    repres(3,10,11) =  b
    repres(3,11,10) = -b
    repres(3,11,11) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,7,7) = 1.0_ark
    repres(4,8,9) = 1.0_ark
    repres(4,9,8) = 1.0_ark
    !
    repres(4,10,10) =  1.0_ark
    repres(4,11,11) = -1.0_ark
    !
    !C2'/(9)->(34)
    repres(5,3,4) = 1.0_ark
    repres(5,4,3) = 1.0_ark
    repres(5,5,5) = 1.0_ark
    !
    repres(5,7,8)  = 1.0_ark
    repres(5,8,7)  = 1.0_ark
    repres(5,9,9)  = 1.0_ark
    !
    repres(5,10,10) = -a
    repres(5,10,11) =  b
    repres(5,11,10) =  b
    repres(5,11,11) =  a
    !
    !(13)->(35)
    repres(6,3,5) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,3) = 1.0_ark
    !
    repres(6,7,9) = 1.0_ark
    repres(6,8,8) = 1.0_ark
    repres(6,9,7) = 1.0_ark
    !
    repres(6,10,10) = -a
    repres(6,10,11) = -b
    repres(6,11,10) = -b
    repres(6,11,11) =  a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo
    !
    dst(12,1) = src(12)
    dst(12,2) = src(12)+2.0_ark*pi/3.0_ark
    dst(12,3) = src(12)-2.0_ark*pi/3.0_ark
    dst(12,4) =-src(12)
    dst(12,5) =-src(12)-2.0_ark*pi/3.0_ark
    dst(12,6) =-src(12)+2.0_ark*pi/3.0_ark
    !
  end subroutine ML_symmetry_transformation_CH3OH


 subroutine ML_symmetry_transformation_CH3OH_II(nsym,src,dst,ndeg)
    implicit none
    !
    integer(ik),intent(in)    :: nsym,ndeg  ! group operation  
    real(ark),intent(in)      :: src(1:ndeg)
    real(ark),intent(out)     :: dst(1:ndeg,6)
    !
    real(ark)         :: repres(nsym,ndeg,ndeg),a,b,e,o
    integer(rk)       :: ioper
    !
    if (verbose>=6) write(out,"('ML_symmetry_transformation_CH3OH_II/start')")
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    !
    repres(:,12,12) = 1.0_ark
    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !d9
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    !
    !C3+/(132)
    repres(2,3,5) = 1.0_ark
    repres(2,4,3) = 1.0_ark
    repres(2,5,4) = 1.0_ark
    !
    repres(2,7,9) = 1.0_ark
    repres(2,8,7) = 1.0_ark
    repres(2,9,8) = 1.0_ark
    !
    repres(2,10,10) = -a
    repres(2,10,11) = -b
    repres(2,11,10) =  b
    repres(2,11,11) = -a
    !
    !C3-/(93)
    !
    repres(3,3,4) = 1.0_ark
    repres(3,4,5) = 1.0_ark
    repres(3,5,3) = 1.0_ark
    !
    repres(3,7,8) = 1.0_ark
    repres(3,8,9) = 1.0_ark
    repres(3,9,7) = 1.0_ark
    !
    repres(3,10,10) = -a
    repres(3,10,11) =  b
    repres(3,11,10) = -b
    repres(3,11,11) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,7,7) = 1.0_ark
    repres(4,8,9) = 1.0_ark
    repres(4,9,8) = 1.0_ark
    !
    repres(4,10,10) =  1.0_ark
    repres(4,11,11) = -1.0_ark
    !
    !C2'/(9)->(34)
    repres(5,3,4) = 1.0_ark
    repres(5,4,3) = 1.0_ark
    repres(5,5,5) = 1.0_ark
    !
    repres(5,7,8)  = 1.0_ark
    repres(5,8,7)  = 1.0_ark
    repres(5,9,9)  = 1.0_ark
    !
    repres(5,10,10) = -a
    repres(5,10,11) =  b
    repres(5,11,10) =  b
    repres(5,11,11) =  a
    !
    !(13)->(35)
    repres(6,3,5) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,3) = 1.0_ark
    !
    repres(6,7,9) = 1.0_ark
    repres(6,8,8) = 1.0_ark
    repres(6,9,7) = 1.0_ark
    !
    repres(6,10,10) = -a
    repres(6,10,11) = -b
    repres(6,11,10) = -b
    repres(6,11,11) =  a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo
    !
    dst(12,1) = src(12)
    dst(12,2) = src(12)-2.0_ark*pi/3.0_ark
    dst(12,3) = src(12)+2.0_ark*pi/3.0_ark
    dst(12,4) =-src(12)
    dst(12,5) =-src(12)+2.0_ark*pi/3.0_ark
    dst(12,6) =-src(12)-2.0_ark*pi/3.0_ark
    !
    if (verbose>=6) write(out,"('ML_symmetry_transformation_CH3OH_II/end')")
    !
  end subroutine ML_symmetry_transformation_CH3OH_II


 !
end module pot_ch3oh

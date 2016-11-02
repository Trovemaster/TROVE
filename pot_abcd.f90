!
!  This unit defines all specific routines for a fouratomic molecule of ABCD type
!
module pot_abcd
  use accuracy
  use moltype
  use lapack

  implicit none

  public MLloc2pqr_abcd,MLdms2pqr_abcd,ML_MEP_ABCD_tau_ref,MLpoten_h2o2_koput,MLpoten_hsoh
  public mlpoten_hsoh_ref,MLdms2xyz_abcd,MLpoten_h2o2_koput_morse,MLdms_hooh_MB,MLpoten_h2o2_koput_unique,MLpoten_v_c2h2_katy,MLpoten_v_c2h2_mlt
  public MLpoten_c2h2_morse,MLpoten_c2h2_7,MLpoten_c2h2_7_xyz,MLpoten_c2h2_streymills,MLdms_HCCH_MB
  
  private
 
  integer(ik), parameter :: verbose     = 4                        ! Verbosity level


  contains
  
  !
  ! Defining potential energy function 
  !
  ! This type is for ABCD-molecules, Morse+cos+Taylor type of expansion 
  ! with respect to the symmetrized GD-coordinates
  ! V =  sum_{ind} f_ind prod(xi(:)**ind(:)) 
  ! with xi(1:3) = 1-exp(-a(1:3)*(r(1:3)-re(1:3)))
  !      xi(4:6) = cos-type or taylor type 
  !
  function MLpoten_hsoh(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          ::  i,k(6),i1
   real(ark)    :: xi(6),y(6),beta(3),rho,req(6)

   if (verbose>=6) write(out,"('MLpoten_hsoh/start')") 
   !
   beta(1:3) = molec%specparam(1:3)
   !
   rho = local(6)
   !
   select case(trim(molec%meptype))
   case default
        !
        req(1:3) = molec%req(1:3)
        req(4:5) = molec%alphaeq(1:2)
        !
   case('MEP_ABCD_TAU_REF') 
        !
        req(1:6) =  ML_MEP_ABCD_tau_ref(rho)
        !
   end select
   !
   y(1:3)=1.0_ark-exp(-beta(1:3)*(local(1:3)-req(1:3)))
   y(4)= local(4)-req(4)
   y(5)= local(5)-req(5)
   !
   !rho = mod(rho+2.0_ark*pi,2.0_ark*pi)
   !
   !if (rho>pi) then 
   !   rho = 2.0_ark*pi-rho
   !endif 
   !
   y(6) = cos(rho) - cos(molec%taueq(1))
   !
   f = 0
   do i = 1,molec%parmax
      k(:) = molec%pot_ind(:,i)
      do i1 = 1,6
         !xi(i1) =1.0_ark ; if (k(i1)/=0) 
         !
         xi(i1) = y(i1)**k(i1)
         !
      enddo
      !
      f = f + force(i)*product(xi(:))
      !
   enddo

   if (verbose>=6) write(out,"('MLpoten_hsoh/end')") 
 
 end function MLpoten_hsoh


  !
  ! Defining potential energy function 
  !
  ! This type is for ABCD-molecules, Morse+cos+Taylor type of expansion 
  ! with respect to the symmetrized GD-coordinates
  ! V =  sum_{ind} f_ind prod(xi(:)**ind(:)) 
  ! with xi(1:3) = 1-exp(-a(1:3)*(r(1:3)-re(1:3)))
  !      xi(4:6) = cos-type or taylor type 
  ! This is for the expansion around a non-rigid reference configuration 
  !
  function mlpoten_hsoh_ref(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          ::  i,k(6),i1,i0
   real(ark)    :: xi(6),y(6),beta(3),rho,re(6),rref(5,0:4)

   if (verbose>=6) write(out,"('MLpoten_hsoh_ref/start')") 
   !
   rho = local(6)
   !
   if (rho>pi) then 
      rho = 2.0_ark*pi-rho
   endif 
   !
   !re(1:6) = ML_MEP_ABCD_tau_ref(rho)
   !
   rref(1,0:4)  = force(1:5) 
   rref(2,0:4)  = force(6:10) 
   rref(3,0:4)  = force(11:15) 
   rref(4,0:4)  = force(16:20)/rad 
   rref(5,0:4)  = force(21:25)/rad
   
   k(1)=1; k(2)=2; k(3)=3; k(4)=4
   
   !re(1)    = rref(1,0)+sum(rref(1,1:4)*cos(rho)**k(1:4))
   !re(2)    = rref(2,0)+sum(rref(2,1:4)*cos(rho)**k(1:4))
   !re(3)    = rref(3,0)+sum(rref(3,1:4)*cos(rho)**k(1:4))
   !re(4)    = rref(4,0)+sum(rref(4,1:4)*cos(rho)**k(1:4))
   !re(5)    = rref(5,0)+sum(rref(5,1:4)*cos(rho)**k(1:4))
   
   re(1)    = rref(1,0)+sum(rref(1,1:4)*(1.0_ark+cos(rho))**k(1:4))
   re(2)    = rref(2,0)+sum(rref(2,1:4)*(1.0_ark+cos(rho))**k(1:4))
   re(3)    = rref(3,0)+sum(rref(3,1:4)*(1.0_ark+cos(rho))**k(1:4))
   re(4)    = rref(4,0)+sum(rref(4,1:4)*(1.0_ark+cos(rho))**k(1:4))
   re(5)    = rref(5,0)+sum(rref(5,1:4)*(1.0_ark+cos(rho))**k(1:4))
   !
   beta(1:3) = molec%specparam(1:3)
   !
   y(1)=1.0_ark-exp(-beta(1)*(local(1)-re(1)))
   y(2)=1.0_ark-exp(-beta(2)*(local(2)-re(2)))
   y(3)=1.0_ark-exp(-beta(3)*(local(3)-re(3)))
   y(4)= local(4)-(4)
   y(5)= local(5)-(5)
   !
   y(6) = cos(rho) - cos(molec%taueq(1))

   f = 0
   i0 = 26
   do i = i0,molec%parmax
      k(:) = molec%pot_ind(:,i)
      do i1 = 1,6
         xi(i1) =1.0_ark ; if (k(i1)/=0) xi(i1) = y(i1)**k(i1)
      enddo
      !
      f = f + force(i)*product(xi(:))
      !
   enddo

   if (verbose>=6) write(out,"('MLpoten_hsoh_ref/end')") 
 
 end function MLpoten_hsoh_ref
 !


!     
!     Jacek Koput
!     updated 12/01/2010
!     the potential energy surface of small-amplitude
!     vibrations and internal rotation of HOOH:
!     q1 - OO, q2 - OHa, q3 - OHb, q4 - OOHa, q5 - OOHb, q6 -HOOH
!     bond stretches in the SPF coordinates
!     energy in atomic units !

function MLpoten_h2o2_koput(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
      !
      integer(ik)          ::  i,i1,i2,i3,i4,i5,i6
      real(ark)    :: x1,x2,x3,x4,x5,x6,e1,e2,e4,e6,vpot,cphi,q1,q2,q3,q4,q5,q6,pd

!     
      real(ark),parameter :: tocm = 219474.63067_ark
      !
      x1 = local(1)
      x2 = local(2)
      x3 = local(3)
      x4 = local(4)
      x5 = local(5)
      x6 = local(6)
      !
      pd=pi/180.0_ark
      e1=force(1)
      e2=force(2)
      e4=force(3)*pd
      e6=force(4)*pd
      !
      q1=(x1-e1)/x1
      q2=(x2-e2)/x2
      q3=(x3-e2)/x3
      q4=(x4-e4)
      q5=(x5-e4)
      q6=x6
      !
      vpot=0.0_ark
      !
      do i=5,molec%parmax
        !
        i1=molec%pot_ind(1,i)
        i2=molec%pot_ind(2,i)
        i3=molec%pot_ind(3,i)
        i4=molec%pot_ind(4,i)
        i5=molec%pot_ind(5,i)
        i6=molec%pot_ind(6,i)
        !
        cphi=1.0_ark
        if (i6.ne.0) cphi=cos(real(i6,ark)*q6)
        !
        vpot=vpot+force(i)*q1**i1*q2**i2*q3**i3*q4**i4*q5**i5*cphi
        !
      enddo
      !
      f = tocm*vpot
      !
  end function MLpoten_h2o2_koput


!     
!     Jacek Koput
!     updated 12/01/2010
!     the potential energy surface of small-amplitude
!     vibrations and internal rotation of HOOH:
!     q1 - OO, q2 - OHa, q3 - OHb, q4 - OOHa, q5 - OOHb, q6 -HOOH
!     bond stretches in the SPF coordinates
!     energy in atomic units !

function MLpoten_h2o2_koput_unique(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
      !
      integer(ik)          ::  i,i1,i2,i3,i4,i5,i6
      real(ark)    :: x1,x2,x3,x4,x5,x6,e1,e2,e4,e6,vpot,cphi,q1,q2,q3,q4,q5,q6,pd

!     
      real(ark),parameter :: tocm = 219474.63067_ark
      !
      x1 = local(1)
      x2 = local(2)
      x3 = local(3)
      x4 = local(4)
      x5 = local(5)
      x6 = local(6)
      !
      pd=pi/180.0_ark
      e1=force(1)
      e2=force(2)
      e4=force(3)*pd
      e6=force(4)*pd
      !
      q1=(x1-e1)/x1
      q2=(x2-e2)/x2
      q3=(x3-e2)/x3
      q4=(x4-e4)
      q5=(x5-e4)
      q6=x6
      !
      vpot=0.0_ark
      !
      do i=5,molec%parmax
        !
        i1=molec%pot_ind(1,i)
        i2=molec%pot_ind(2,i)
        i3=molec%pot_ind(3,i)
        i4=molec%pot_ind(4,i)
        i5=molec%pot_ind(5,i)
        i6=molec%pot_ind(6,i)
        !
        cphi=1.0_ark
        !
        if (i6.ne.0) cphi=cos(real(i6,ark)*q6)
        !
        vpot=vpot+force(i)*q1**i1*q2**i2*q3**i3*q4**i4*q5**i5*cphi
        !
        if (i2/=i3.or.i4/=i5) then
          !
          i2 = molec%pot_ind(3,i)
          i3 = molec%pot_ind(2,i)
          i4 = molec%pot_ind(5,i)
          i5 = molec%pot_ind(4,i)
          !
          !vpot=vpot+force(i)*product(y(1:6))
          !
          vpot=vpot+force(i)*q1**i1*q2**i2*q3**i3*q4**i4*q5**i5*cphi
          !
        endif
        !
      enddo
      !
      f = tocm*vpot
      !
  end function MLpoten_h2o2_koput_unique




function MLpoten_h2o2_koput_morse(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
      !
      integer(ik)          ::  i,i1,i2,i3,i4,i5,i6,k_ind(6)
      real(ark)    :: x1,x2,x3,x4,x5,x6,e1,e2,e4,e6,vpot,cphi,q(6),y(6),a1,a2,pd
      !     
      real(ark),parameter :: tocm = 219474.63067_ark
      !
      x1 = local(1)
      x2 = local(2)
      x3 = local(3)
      x4 = local(4)
      x5 = local(5)
      x6 = local(6)
      !
      pd=pi/180.0_ark
      e1=force(1)
      e2=force(2)
      e4=force(3)*pd
      e6=force(4)*pd
      !
      a1 = force(5)
      a2 = force(6)
      !
      q(1)=1.0d0-exp(-a1*(x1-e1))
      q(2)=1.0d0-exp(-a2*(x2-e2))
      q(3)=1.0d0-exp(-a2*(x3-e2))
      q(4)=(x4-e4)
      q(5)=(x5-e4)
      q(6)=x6
      !
      vpot=0.0_ark
      !
      do i=7,molec%parmax
        !
        y(1:5) = q(1:5)**molec%pot_ind(1:5,i)
        y(6)=cos(real(molec%pot_ind(6,i),4)*q(6))
        !
        vpot=vpot+force(i)*product(y(1:6))
        !
        if (molec%pot_ind(2,i)/=molec%pot_ind(3,i).or.molec%pot_ind(4,i)/=molec%pot_ind(5,i)) then
          !
          k_ind(2) = molec%pot_ind(3,i)
          k_ind(3) = molec%pot_ind(2,i)
          k_ind(4) = molec%pot_ind(5,i)
          k_ind(5) = molec%pot_ind(4,i)
          !
          y(2:5) = q(2:5)**k_ind(2:5)
          !
          vpot=vpot+force(i)*product(y(1:6))
          !
        endif
        !
      enddo
      !
      f = vpot
      !
  end function MLpoten_h2o2_koput_morse


function MLpoten_c2h2_morse(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
      !
      integer(ik)          ::  i,i1,i2,i3,i4,i5,i6,k_ind(6)
      real(ark)    :: x1,x2,x3,x4,x5,x6,e1,e2,e4,e6,vpot,cphi,q(6),y(6),a1,a2,pd,rc1c2,rc1h1,rc2h2,delta1x,delta1y,delta2x,delta2y,tau
      real(ark)    :: alpha1,alpha2,sinalpha2,sinalpha1,tau1,tau2,v1(3),v2(3),v3(3),v12(3),v31(3),r21,r31,n3(3),sintau,costau,cosalpha1,cosalpha2
      character(len=cl)  :: txt = 'MLpoten_c2h2_morse'
      !
      integer(ik)  :: Nangles 
      !
      Nangles = molec%Nangles
      !
      pd=pi/180.0_ark
      e1=force(1)
      e2=force(2)
      e4=force(3)*pd
      e6=force(4)*pd
      !
      a1 = force(5)
      a2 = force(6)
      !
      rc1c2    = local(1)
      rc1h1    = local(2)
      rc2h2    = local(3)
      !
      if (molec%zmatrix(3)%connect(4)==101) then 
         !
         delta1x = local(4)
         delta1y = local(5)
         !
         delta2x = local(6)
         delta2y = local(7)
         !
         tau1 = pi-atan2(sin(delta1y),-sin(delta1x))
         tau2 = pi-atan2(sin(delta2y),-sin(delta2x))
         !
         !tau1 = pi-atan2((delta1y),-(delta1x))
         !tau2 = pi-atan2((delta2y),-(delta2x))
         !
         tau = tau2-tau1
         !
         sinalpha1 = sqrt(sin(delta1x)**2+sin(delta1y)**2)
         sinalpha2 = sqrt(sin(delta2x)**2+sin(delta2y)**2)
         !
         !sinalpha1 = sqrt((delta1x)**2+(delta1y)**2)
         !sinalpha2 = sqrt((delta2x)**2+(delta2y)**2)
         !
         alpha1 = pi-asin(sinalpha1)
         alpha2 = pi-asin(sinalpha2)
         !
         x1 = rc1c2
         x3 = rc1h1
         x2 = rc2h2
         x5 = alpha1
         x4 = alpha2
         !
         x6 = cos(tau)*sin(alpha1)*sin(alpha2)
         !
      elseif (molec%zmatrix(3)%connect(4)==102) then 
         !
         delta1x = local(4)
         delta1y = local(5)
         !
         delta2x = local(6)
         delta2y = local(7)
         !
         !tau1 = pi-atan2(sin(delta1y),-sin(delta1x))
         !tau2 = pi-atan2(sin(delta2y),-sin(delta2x))
         !
         tau1 = pi-atan2((delta1y),-(delta1x))
         tau2 = pi-atan2((delta2y),-(delta2x))
         !
         tau = tau2-tau1
         !
         !sinalpha1 = sqrt(sin(delta1x)**2+sin(delta1y)**2)
         !sinalpha2 = sqrt(sin(delta2x)**2+sin(delta2y)**2)
         !
         sinalpha1 = sqrt((delta1x)**2+(delta1y)**2)
         sinalpha2 = sqrt((delta2x)**2+(delta2y)**2)
         !
         alpha1 = pi-asin(sinalpha1)
         alpha2 = pi-asin(sinalpha2)
         !
         x1 = rc1c2
         x3 = rc1h1
         x2 = rc2h2
         x5 = alpha1
         x4 = alpha2
         x6 = tau
         !
         x6 = cos(tau)*sin(alpha1)*sin(alpha2)
         ! 
      elseif (molec%zmatrix(3)%connect(4)==103) then
         !
         v1(:) = xyz(2,:)-xyz(1,:)
         v2(:) = xyz(3,:)-xyz(1,:)
         v3(:) = xyz(4,:)-xyz(2,:)
         !
         x1 = sqrt(sum(v1(:)**2))
         x2 = sqrt(sum(v2(:)**2))
         x3 = sqrt(sum(v3(:)**2))
         !
         cosalpha1 = sum( v1(:)*v2(:))/(x1*x2)
         cosalpha2 = sum(-v1(:)*v3(:))/(x1*x3)
         !
         x4 = aacos(cosalpha1,txt)
         x5 = aacos(cosalpha2,txt)
         !
         v12(:) = MLvector_product(v2(:),v1(:))
         v31(:) = MLvector_product(v3(:),v1(:))
         !
         v12(:) = v12(:)/(x1*x2)
         v31(:) = v31(:)/(x1*x3)
         !
         x6 = sum(v12*v31)
         !
         r21 = sqrt(sum(v12(:)**2))
         r31 = sqrt(sum(v31(:)**2))
         !
         tau = 0
         !
         x6 = 0
         !
         if (r21>sqrt(small_).and.r31>sqrt(small_)) then 
           !
           v12(:) = v12(:)/r21
           v31(:) = v31(:)/r31
           !
           n3(:) = MLvector_product(v12(:),v31(:))
           !
           sintau =  sqrt( sum( n3(:)**2 ) )
           costau =  sum( v12(:)*v31(:) )
           !
           tau = atan2(sintau,costau)
           tau = aacos(costau,txt)
           !
           !q(6)= x6
           !
           x6 = tau
           !
           !if (abs(x6-q(6))>sqrt(small_)) then
           !  !
           !  continue 
           !  !
           !endif 
           
           !
         endif
         !
         !cosalpha1 =  v2(3)/rc1h1
         !cosalpha2 = -v3(3)/rc2h2
         !
         !x4 = aacos(cosalpha1,txt)
         !x5 = aacos(cosalpha2,txt)
         !
         !x6 = (v2(1)*v3(1)+v2(2)*v3(2))/(rc1h1*rc2h2)
         !
         !e4 = cos(e4)
         !
         !
         !n3(:) = MLvector_product(v12(:),v31(:))
         !
         !sintau =  sqrt( sum( n3(:)**2 ) )
         !costau = -sum( v12(:)*v31(:) )
         !
         !tau = atan2(sintau,costau)
         !
      else 
        !
        write(out,"('MLpoten_v_c2h2_katy: only designed for 7 coordinates, not ',i6)") ncoords
        stop 'only designed for 7 coordinates ' 
        !
      endif
      !
      q(1)=1.0_ark-exp(-a1*(x1-e1))
      q(2)=1.0_ark-exp(-a2*(x2-e2))
      q(3)=1.0_ark-exp(-a2*(x3-e2))
      !
      q(4)=(x4-e4)
      q(5)=(x5-e4)
      q(6)=x6
      !
      f = 0 
      !
      vpot=0.0_ark
      !
      do i=7,molec%parmax
        !
        if ((molec%pot_ind(4,i)==0.or.molec%pot_ind(5,i)==0).and.molec%pot_ind(6,i)/=0) return
        !
        y(1:5) = q(1:5)**molec%pot_ind(1:5,i)
        !
        y(6) = 1.0_ark
        !
        if (molec%pot_ind(4,i)/=0.and.molec%pot_ind(5,i)/=0) y(6)=cos(real(molec%pot_ind(6,i),ark)*q(6))
        !
        !y(6)=cos(real(molec%pot_ind(6,i),4)*q(6))
        !
        !if (molec%pot_ind(6,i)/=0) y(6)=y(6)*sin(x4)*sin(x5)
        !
        !y(6) =q(6)**molec%pot_ind(6,i)
        !
        !y(6)=(cos(x6)*sin(x4)*sin(x5))**molec%pot_ind(6,i)
        !
        vpot=vpot+force(i)*product(y(1:6))
        !
        if (molec%pot_ind(2,i)/=molec%pot_ind(3,i).or.molec%pot_ind(4,i)/=molec%pot_ind(5,i)) then
          !
          k_ind(2) = molec%pot_ind(3,i)
          k_ind(3) = molec%pot_ind(2,i)
          k_ind(4) = molec%pot_ind(5,i)
          k_ind(5) = molec%pot_ind(4,i)
          !
          y(2:5) = q(2:5)**k_ind(2:5)
          !
          vpot=vpot+force(i)*product(y(1:6))
          !
        endif
        !
      enddo
      !
      f = vpot
      !
  end function MLpoten_c2h2_morse



function MLpoten_c2h2_morse_kappa(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
      !
      integer(ik)          ::  i,i1,i2,i3,i4,i5,i6,k_ind(6)
      real(ark)    :: x1,x2,x3,x4,x5,x6,e1,e2,e4,e6,vpot,cphi,q(6),y(6),a1,a2,pd,rc1c2,rc1h1,rc2h2,delta1x,delta1y,delta2x,delta2y,tau
      real(ark)    :: alpha1,alpha2,sinalpha2,sinalpha1,tau1,tau2,v1(3),v2(3),v3(3),v12(3),v31(3),r21,r31,n3(3),sintau,costau,cosalpha1,cosalpha2
      character(len=cl)  :: txt = 'MLpoten_c2h2_morse'
      !
      integer(ik)  :: Nangles 
      !
      Nangles = molec%Nangles
      !
      pd=pi/180.0_ark
      e1=force(1)
      e2=force(2)
      e4=force(3)*pd
      e6=force(4)*pd
      !
      a1 = force(5)
      a2 = force(6)
      !
      rc1c2    = local(1)
      rc1h1    = local(2)
      rc2h2    = local(3)
      !
      if (molec%zmatrix(3)%connect(4)==101) then 
         !
         delta1x = local(4)
         delta1y = local(5)
         !
         delta2x = local(6)
         delta2y = local(7)
         !
         tau1 = pi-atan2(sin(delta1y),-sin(delta1x))
         tau2 = pi-atan2(sin(delta2y),-sin(delta2x))
         !
         !tau1 = pi-atan2((delta1y),-(delta1x))
         !tau2 = pi-atan2((delta2y),-(delta2x))
         !
         tau = tau2-tau1
         !
         sinalpha1 = sqrt(sin(delta1x)**2+sin(delta1y)**2)
         sinalpha2 = sqrt(sin(delta2x)**2+sin(delta2y)**2)
         !
         !sinalpha1 = sqrt((delta1x)**2+(delta1y)**2)
         !sinalpha2 = sqrt((delta2x)**2+(delta2y)**2)
         !
         alpha1 = pi-asin(sinalpha1)
         alpha2 = pi-asin(sinalpha2)
         !
         x1 = rc1c2
         x3 = rc1h1
         x2 = rc2h2
         x5 = alpha1
         x4 = alpha2
         !
         x6 = cos(tau)*sin(alpha1)*sin(alpha2)
         !
      elseif (molec%zmatrix(3)%connect(4)==102) then 
         !
         delta1x = local(4)
         delta1y = local(5)
         !
         delta2x = local(6)
         delta2y = local(7)
         !
         !tau1 = pi-atan2(sin(delta1y),-sin(delta1x))
         !tau2 = pi-atan2(sin(delta2y),-sin(delta2x))
         !
         tau1 = pi-atan2((delta1y),-(delta1x))
         tau2 = pi-atan2((delta2y),-(delta2x))
         !
         tau = tau2-tau1
         !
         !sinalpha1 = sqrt(sin(delta1x)**2+sin(delta1y)**2)
         !sinalpha2 = sqrt(sin(delta2x)**2+sin(delta2y)**2)
         !
         sinalpha1 = sqrt((delta1x)**2+(delta1y)**2)
         sinalpha2 = sqrt((delta2x)**2+(delta2y)**2)
         !
         alpha1 = pi-asin(sinalpha1)
         alpha2 = pi-asin(sinalpha2)
         !
         x1 = rc1c2
         x3 = rc1h1
         x2 = rc2h2
         x5 = alpha1
         x4 = alpha2
         x6 = tau
         !
         x6 = cos(tau)*sin(alpha1)*sin(alpha2)
         ! 
      elseif (molec%zmatrix(3)%connect(4)==103) then
         !
         v1(:) = xyz(2,:)-xyz(1,:)
         v2(:) = xyz(3,:)-xyz(1,:)
         v3(:) = xyz(4,:)-xyz(2,:)
         !
         x1 = sqrt(sum(v1(:)**2))
         x2 = sqrt(sum(v2(:)**2))
         x3 = sqrt(sum(v3(:)**2))
         !
         cosalpha1 = sum( v1(:)*v2(:))/(x1*x2)
         cosalpha2 = sum(-v1(:)*v3(:))/(x1*x3)
         !
         x4 = aacos(cosalpha1,txt)
         x5 = aacos(cosalpha2,txt)
         !
         v12(:) = MLvector_product(v2(:),v1(:))
         v31(:) = MLvector_product(v3(:),v1(:))
         !
         v12(:) = v12(:)/(x1*x2)
         v31(:) = v31(:)/(x1*x3)
         !
         x6 = sum(v12*v31)
         !
         r21 = sqrt(sum(v12(:)**2))
         r31 = sqrt(sum(v31(:)**2))
         !
         if (r21>sqrt(small_).and.r31>sqrt(small_)) then 
           !
           v12(:) = v12(:)/r21
           v31(:) = v31(:)/r31
           !
           n3(:) = MLvector_product(v12(:),v31(:))
           !
           sintau =  sqrt( sum( n3(:)**2 ) )
           costau = -sum( v12(:)*v31(:) )
           !
           tau = atan2(sintau,costau)
           !
           !q(6)= x6
           !
           x6 = cos(tau)*sin(x4)*sin(x5)
           !
           !if (abs(x6-q(6))>sqrt(small_)) then
           !  !
           !  continue 
           !  !
           !endif 
           
           !
         endif
         !
         !cosalpha1 =  v2(3)/rc1h1
         !cosalpha2 = -v3(3)/rc2h2
         !
         !x4 = aacos(cosalpha1,txt)
         !x5 = aacos(cosalpha2,txt)
         !
         !x6 = (v2(1)*v3(1)+v2(2)*v3(2))/(rc1h1*rc2h2)
         !
         !e4 = cos(e4)
         !
         !
         !n3(:) = MLvector_product(v12(:),v31(:))
         !
         !sintau =  sqrt( sum( n3(:)**2 ) )
         !costau = -sum( v12(:)*v31(:) )
         !
         !tau = atan2(sintau,costau)
         !
      else 
        !
        write(out,"('MLpoten_v_c2h2_katy: only designed for 7 coordinates, not ',i6)") ncoords
        stop 'only designed for 7 coordinates ' 
        !
      endif
      !
      q(1)=1.0_ark-exp(-a1*(x1-e1))
      q(2)=1.0_ark-exp(-a2*(x2-e2))
      q(3)=1.0_ark-exp(-a2*(x3-e2))
      !
      q(4)=(x4-e4)
      q(5)=(x5-e4)
      q(6)=x6
      !
      vpot=0.0_ark
      !
      do i=7,molec%parmax
        !
        y(1:6) = q(1:6)**molec%pot_ind(1:6,i)
        !
        !y(6)=cos(real(molec%pot_ind(6,i),4)*q(6))
        !
        !if (molec%pot_ind(6,i)/=0) y(6)=y(6)*sin(x4)*sin(x5)
        !
        y(6) =q(6)**molec%pot_ind(6,i)
        !
        !y(6)=(cos(x6)*sin(x4)*sin(x5))**molec%pot_ind(6,i)
        !
        vpot=vpot+force(i)*product(y(1:6))
        !
        if (molec%pot_ind(2,i)/=molec%pot_ind(3,i).or.molec%pot_ind(4,i)/=molec%pot_ind(5,i)) then
          !
          k_ind(2) = molec%pot_ind(3,i)
          k_ind(3) = molec%pot_ind(2,i)
          k_ind(4) = molec%pot_ind(5,i)
          k_ind(5) = molec%pot_ind(4,i)
          !
          y(2:5) = q(2:5)**k_ind(2:5)
          !
          vpot=vpot+force(i)*product(y(1:6))
          !
        endif
        !
      enddo
      !
      f = vpot
      !
  end function MLpoten_c2h2_morse_kappa



 !===============================================================================
 !                   electric dipole moment section
 !===============================================================================




 function MLloc2pqr_abcd(r) result(f)
    !
    !define cartesian coordinates of atoms in space-fixed system
    !
    ! s  = (0,0,0)
    ! so = (0,0,+z)
    ! sh = (0,+y,z)
    !
    !       +z
    !       |
    !       |
    !       |
    !       | H
    !       |/
    !       O
    !       |
    !       S-------- +x
    !      / 
    !     H
    !    /
    !   /
    !  /
    ! +y
    !
    real(ark), intent(in) :: r(molec%ncoords)
    real(ark)             :: f(molec%natoms, 3)
    !
    integer(ik)           :: icart
    real(ark)             :: rho, a0(molec%natoms, 3), cm
    !
    rho = r(6)
    !
    a0(1, 1) = 0.0_ark
    a0(1, 2) = 0.0_ark
    a0(1, 3) = 0.0_ark
    !
    a0(2, 1) = 0.0_ark
    a0(2, 2) = 0.0_ark
    a0(2, 3) = r(1)
    !
    a0(3, 1) = 0.0_ark
    a0(3, 2) = r(2) * sin(r(4))
    a0(3, 3) = r(2) * cos(r(4))
    !
    a0(4, 1) = r(3) * sin(r(5)) * sin(rho)
    a0(4, 2) = r(3) * sin(r(5)) * cos(rho)
    a0(4, 3) = r(1) - r(3) * cos(r(5))
    !
    do icart = 1, 3
       cm = sum(a0(1:4, icart) * molec%atommasses(1:4)) / sum(molec%atommasses(1:4))
       a0(1:4, icart) = a0(1:4, icart) - cm
    end do
    !
    f(1:molec%natoms, 1:3) = a0(1:molec%natoms, 1:3)
    !
 end function MLloc2pqr_abcd

 function MLdms2pqr_abcd(r) result(f)
    !
    !define cartesian components of the dipole moment in space-fixed system
    !
    ! mu_x has transformation properties of sin(tau)
    ! mu_y has transformation properties of cos(tau)
    !
    real(ark), intent(in) :: r(molec%ncoords)
    real(ark)             :: f(3)
    !
    integer(ik)           :: imu, iterm
    real(ark)             :: r1, r2, r3, alpha1, alpha2, tau
    real(ark)             :: re1(1:3),re2(1:3),re3(1:3),alphae1(1:3),alphae2(1:3),taue(1:3), &
                             beta1(1:3),beta2(1:3),beta3(1:3), y(molec%ncoords,1:3), mu(3), xi(molec%ncoords), tau_
    !
    r1     = r(1)
    r2     = r(2)
    r3     = r(3)
    alpha1 = r(4)
    alpha2 = r(5)
    tau    = r(6)
    !
    if (tau > pi) then
       tau_ = 2.0*pi - tau
    else
       tau_ = tau
    end if
    !
    re1(1:3)     = extF%coef(1,1:3)
    re2(1:3)     = extF%coef(2,1:3)
    re3(1:3)     = extF%coef(3,1:3)
    alphae1(1:3) = extF%coef(4,1:3)/rad 
    alphae2(1:3) = extF%coef(5,1:3)/rad 
    taue(1:3)    = extF%coef(6,1:3)/rad 
    !
    beta1(1:3)   = extF%coef(7,1:3)
    beta2(1:3)   = extF%coef(8,1:3)
    beta3(1:3)   = extF%coef(9,1:3)
    !
    y(1,:) = (r1 - re1(:)) * exp(-beta1(:) * (r1 - re1(:)) ** 2)
    y(2,:) = (r2 - re2(:)) * exp(-beta2(:) * (r2 - re2(:)) ** 2)
    y(3,:) = (r3 - re3(:)) * exp(-beta3(:) * (r3 - re3(:)) ** 2)
    y(4,:) = (alpha1 - alphae1(:))
    y(5,:) = (alpha2 - alphae2(:))
    y(6,:) = cos(tau_) - cos(taue(:))
    !
    mu = 0
    !
    do imu = 1, 3
       !
       do iterm = 10, extF%nterms(imu)
          !
          xi(1:molec%ncoords) = y(1:molec%ncoords,imu) ** extF%term(1:molec%ncoords, iterm, imu)
          !
          mu(imu) = mu(imu) + extF%coef(iterm, imu) * product(xi(1:molec%ncoords))
          !
       end do
       !
    end do
    !
    mu(1) = mu(1)*sin(tau)
    !
    !if (tau > pi) then
    !   mu(1) = -mu(1)
    !end if
    !
    f(1:3) = mu(1:3)
    !
 end function MLdms2pqr_abcd



 !
 !define cartesian components of the dipole moment in space-fixed system
 !
 ! mu_x has transformation properties of sin(tau)
 ! mu_y has transformation properties of cos(tau)
 ! function MLdms2xyz_abcd(xyz,r) result(f)
 !
 recursive subroutine MLdms2xyz_abcd(rank,ncoords,natoms,r,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  r(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: imu, iterm
    real(ark)             :: r1, r2, r3, alpha1, alpha2, tau, e1(3),e2(3),e3(3),tmat(3,3)
    real(ark)             :: re1(1:3),re2(1:3),re3(1:3),alphae1(1:3),alphae2(1:3),taue(1:3), &
                             beta1(1:3),beta2(1:3),beta3(1:3), y(molec%ncoords,1:3), mu(3), xi(molec%ncoords), tau_

    !  
    e3(:) =  xyz(2,:)-xyz(1,:)
    !
    e3(:) = e3(:)/sqrt(sum(e3(:)**2))
    !
    e1 = MLvector_product(xyz(3,:)-xyz(1,:),e3)
    !
    e1(:) = e1(:)/sqrt(sum(e1(:)**2))
    !
    e2(:) =  MLvector_product(e3,e1)
    !
    tmat(:,1) = e1
    tmat(:,2) = e2
    tmat(:,3) = e3
    !
    r1     = r(1)
    r2     = r(2)
    r3     = r(3)
    alpha1 = r(4)
    alpha2 = r(5)
    tau    = r(6)
    !
    if (tau > pi) then
       tau_ = 2.0*pi - tau
    else
       tau_ = tau
    end if
    !
    re1(1:3)     = extF%coef(1,1:3)
    re2(1:3)     = extF%coef(2,1:3)
    re3(1:3)     = extF%coef(3,1:3)
    alphae1(1:3) = extF%coef(4,1:3)/rad 
    alphae2(1:3) = extF%coef(5,1:3)/rad 
    taue(1:3)    = extF%coef(6,1:3)/rad 
    !
    beta1(1:3)   = extF%coef(7,1:3)
    beta2(1:3)   = extF%coef(8,1:3)
    beta3(1:3)   = extF%coef(9,1:3)
    !
    y(1,:) = (r1 - re1(:)) * exp(-beta1(:) * (r1 - re1(:)) ** 2)
    y(2,:) = (r2 - re2(:)) * exp(-beta2(:) * (r2 - re2(:)) ** 2)
    y(3,:) = (r3 - re3(:)) * exp(-beta3(:) * (r3 - re3(:)) ** 2)
    y(4,:) = (alpha1 - alphae1(:))
    y(5,:) = (alpha2 - alphae2(:))
    y(6,:) = cos(tau_) - cos(taue(:))
    !
    mu = 0
    !
    do imu = 1, 3
       !
       do iterm = 10, extF%nterms(imu)
          !
          xi(1:molec%ncoords) = y(1:molec%ncoords,imu) ** extF%term(1:molec%ncoords, iterm, imu)
          !
          mu(imu) = mu(imu) + extF%coef(iterm, imu) * product(xi(1:molec%ncoords))
          !
       end do
       !
    end do
    !
    mu(1) = mu(1)*sin(tau)
    !
    !if (tau > pi) then
    !   mu(1) = -mu(1)
    !end if
	!
	!write(out,"('mu = ',3f18.8)") mu(1:3)
    !
    f(1:3) = matmul((tmat),mu)
    !
    !f(1:3) = mu(1:3)
    !
 end subroutine MLdms2xyz_abcd

  !
  ! Defining MEP function for CH2 molecule
  !
  function ML_MEP_ABCD_tau_ref(x)  result(dst)

   real(ark),intent(in)   ::  x
   real(ark)              ::  dst(6)
   real(ark)              ::  rho,r(5,0:4),xi
   integer(ik)            ::  k(4)
     !
     rho = x
     !
     r(1,0:4)  = molec%mep_params(1:5) 
     r(2,0:4)  = molec%mep_params(6:10) 
     r(3,0:4)  = molec%mep_params(11:15) 
     r(4,0:4)  = molec%mep_params(16:20)
     r(5,0:4)  = molec%mep_params(21:25)
     !
     k(1)=1; k(2)=2; k(3)=3; k(4)=4
     !
     xi = cos(rho) - cos(molec%taueq(1))
     !
     dst(1) = r(1,0)+sum(r(1,1:4)*xi**k(1:4))
     dst(2) = r(2,0)+sum(r(2,1:4)*xi**k(1:4))
     dst(3) = r(3,0)+sum(r(3,1:4)*xi**k(1:4))
     dst(4) = r(4,0)+sum(r(4,1:4)*xi**k(1:4))
     dst(5) = r(5,0)+sum(r(5,1:4)*xi**k(1:4))
     dst(6) = rho
 
  end function ML_MEP_ABCD_tau_ref



  !
 !define cartesian components of the dipole moment in space-fixed system
 !
 ! mu_x has transformation properties of cos(tau) and r2+r3, a1+a2
 ! mu_y has transformation properties of cos(tau) and r2-r3  a1-a2
 ! mu_z has transformation properties of sin(tau) and r2-r3  a1-a2
 !
 recursive subroutine MLdms_hooh_MB(rank,ncoords,natoms,r,xyz,f)
    !
    implicit none
    !
    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  r(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: imu, iterm, ind(1:molec%ncoords)
    real(ark)             :: mu_t,f_t,r21,r31,r1, r2, r3, alpha1, alpha2, tau, e1(3),e2(3),e3(3),tmat(3,3),n1(3),n2(3),n3(3),v12(3),v31(3)
    real(ark)             :: re1(1:3),re2(1:3),alphae(1:3),e0(3),costau, &
                             beta1(1:3),beta2(1:3),y(molec%ncoords,1:3), mu(3), xi(molec%ncoords), tau_,sintau,r0,tau1,tau2,x1,x2,y1,y2,tau_sign
    !
    integer(ik),parameter :: lspace = 150
    integer(ik) :: ierror,rank0
    real(rk)    :: dip_rk(3, 1), tmat_rk(3, 3), tsing(3), wspace(lspace),tol = -1.0d-12
    character(len=cl)  :: txt

    !
    !write(out,"(i6)") molec%natoms
    !
    !write(out,"(/'O',3x,3f14.8)") xyz(1,:)
    !write(out,"( 'O',3x,3f14.8)") xyz(2,:)
    !write(out,"( 'H',3x,3f14.8)") xyz(3,:)
    !write(out,"( 'H',3x,3f14.8)") xyz(4,:)
    !
    e1(:) = xyz(1,:)-xyz(2,:)
    e2(:) = xyz(3,:)-xyz(1,:)
    e3(:) = xyz(4,:)-xyz(2,:)
    !
    r1 = sqrt(sum(e1(:)**2))
    r2 = sqrt(sum(e2(:)**2))
    r3 = sqrt(sum(e3(:)**2))
    !
    alpha1 = acos(sum(-e1(:)*e2(:))/(r1*r2))
    alpha2 = acos(sum( e1(:)*e3(:))/(r1*r3))
    !
    v12(:) = MLvector_product(e1(:),e2(:))
    v31(:) = MLvector_product(e3(:),e1(:))
    r21 = sqrt(sum(v12(:)**2))
    r31 = sqrt(sum(v31(:)**2))
    v12(:) = v12(:)/r21
    v31(:) = v31(:)/r31
    !
    n3(:) = MLvector_product(v12(:),v31(:))
    !
    sintau =  sqrt( sum( n3(:)**2 ) )
    costau = -sum( v12(:)*v31(:) )
    !
    tau = atan2(sintau,costau)
    !
    e0 = MLvector_product(v12(:),v31(:))
    !
    tau_sign = -sum( e1(:)*e0(:) )
    !
    if (tau_sign<-small_a) then 
       !
       tau = 2.0_ark*pi-tau
       !
    endif
    !
    if ( tau<-small_) then 
      tau  = mod(tau+2.0_ark*pi,2.0_ark*pi)
    endif
    !
    !txt='MLdms_hooh_MB: illegal sin'
    !tau_ = aasin(sintau,txt)
    !
    n3(:) = e1(:)/sqrt(sum(e1(:)**2))
    !
    e2(:) =  MLvector_product(v12(:),n3(:))
    e2(:) =  e2(:)/sqrt(sum(e2(:)**2))
    e3(:) = -MLvector_product(v31(:),n3(:))
    e3(:) =  e3(:)/sqrt(sum(e3(:)**2))

    !
    n1(:) = (v12(:) + v31(:))
    e1(:) = e2(:) + e3(:)
    !
    if (e2(2)>0.or.e3(2)>0) then
      e1 = -e1
    endif
    !
    if (v12(2)>0.or.v31(2)>0) then
      n1 = -n1
    endif
    !
    if (sum(n1(:)**2)>1e-1.and.sum(e1(:)**2)>1e-1) then
      !
      n1(:) = n1(:)/sqrt(sum(n1(:)**2))
      e1(:) = e1(:)/sqrt(sum(e1(:)**2))
      !
      if (any(abs(n1-e1)>1e-12)) then
        !
        write(out,"('MLdms_hooh_MB: n1<>e1 : n1 = ',3g12.5,' e1 = ',3g12.5)") n1,e1
        stop 'MLdms_hooh_MB: n1<>e1 '
        !
      endif 
      !
    elseif (sum(n1(:)**2)>1e-1) then 
      !
      n1(:) = n1(:)/sqrt(sum(n1(:)**2))
      !
    elseif (sum(e1(:)**2)>1e-1) then
      !
      n1(:) = e1(:)
      n1(:) = n1(:)/sqrt(sum(n1(:)**2))
      !
    else 
      !
      write(out,"('MLdms_hooh_MB: both n1 and e1 are  0 : n1 = ',3g12.5,' e1 = ',3g12.5)") n1,e1
      stop 'MLdms_hooh_MB: n1 = e1 = 0 '
      !
    endif
    !
    !n1(:) =  MLvector_product(v12(:),e1(:))
    !n1(:) = n1(:)/sqrt(sum(n1(:)**2))
    !
    n2(:) = MLvector_product(n3(:),n1(:))
    n2(:) = n2(:)/sqrt(sum(n2(:)**2))
    !
    tmat(1,:) = n1(:)
    tmat(2,:) = n2(:)
    tmat(3,:) = n3(:)
    !
    r1 = r(1)
    r2 = r(2)
    r3 = r(3)
    alpha1 = r(4)
    alpha2 = r(5)
    tau_ = r(6)
    !
    tau_ = tau
    !
    if (v12(2)>small_) tau_ = 2.0_ark*pi + tau
    !
    re1(1:3)     = extF%coef(1,1:3)
    re2(1:3)     = extF%coef(2,1:3)
    alphae(1:3) = extF%coef(3,1:3)/rad 
    !
    beta1(1:3)   = extF%coef(4,1:3)
    beta2(1:3)   = extF%coef(5,1:3)
    !
    y(1,:) = (r1 - re1(:)) * exp(-beta1(:) * (r1 - re1(:)) ** 2)
    y(2,:) = (r2 - re2(:)) * exp(-beta2(:) * (r2 - re2(:)) ** 2)
    y(3,:) = (r3 - re2(:)) * exp(-beta2(:) * (r3 - re2(:)) ** 2)
    y(4,:) = (alpha1 - alphae(:))
    y(5,:) = (alpha2 - alphae(:))
    !
    y(6,:) = cos(tau_)
    !
    !y(6,:) = cos(2.0_ark*tau_)
    !
    mu = 0
    !
    do imu = 1, 3
       !
       do iterm =  7, extF%nterms(imu)
          !
          ind(1:6) = extF%term(1:6, iterm, imu)
          xi(1:6) = y(1:6,imu) ** ind(1:6)
          !
          mu_t = product(xi(1:molec%ncoords))
          !
          if (ind(2)/=ind(3).or.ind(4)/=ind(5)) then 
            !
            ind(2) = extF%term(3, iterm, imu)
            ind(3) = extF%term(2, iterm, imu)
            ind(4) = extF%term(5, iterm, imu)
            ind(5) = extF%term(4, iterm, imu)
            !
            xi(2:5) = y(2:5,imu) ** ind(2:5)
            !
            f_t = 1.0_ark
            if (imu/=1)  f_t = -1.0_ark
            !
            mu_t = mu_t + f_t*product(xi(1:molec%ncoords))
            !
          endif
          !
          mu(imu) = mu(imu) + extF%coef(iterm, imu)*mu_t
          !
       end do
       !
    end do
    !
    mu(1) = mu(1)*cos(tau_*0.5_ark)
    mu(2) = mu(2)*sin(tau_*0.5_ark)
    !
    tmat = transpose(tmat)
    !
    call MLlinurark(3,tmat,mu,f,ierror)
    !
    if (ierror>0) then
      !
      tmat_rk = tmat
      dip_rk(:,1) = mu(:)
      !
      call dgelss(3,3,1,tmat_rk,3,dip_rk,3,tsing,tol,rank0,wspace,lspace,ierror)
      !
      f(:) = real(dip_rk(:,1),ark)
      !
      if (ierror>0) then
        !
        print *,ierror,tmat,mu
        write(out,"('MLdms_hooh_MB: dgelss error = ',i)") ierror
        stop 'MLdms_hooh_MB: dgelss error'
        !
      endif
      !
    endif
    !
    f(1:3) = matmul(tmat,mu)
    !
    !f = mu
    !
  end subroutine MLdms_hooh_MB

! HCCH Dipole using R,r1,r2,alpha1,alpha2,tau coordinates
!define cartesian components of the dipole moment in space-fixed system
 !
 ! mu_x has transformation properties of cos(tau) and r2+r3, a1+a2
 ! mu_y has transformation properties of cos(tau) and r2-r3  a1-a2
 ! mu_z has transformation properties of sin(tau) and r2-r3  a1-a2
 !
 recursive subroutine MLdms_HCCH_MB(rank,ncoords,natoms,r,xyz,f)
    !
    implicit none
    !
    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  r(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: imu, iterm, ind(1:molec%ncoords),dimen
    real(ark)             :: mu_t,f_t,r21,r31,r1, r2, r3, alpha1, alpha2, tau, e1(3),e2(3),e3(3),tmat(3,3),tmatout(3,3),n1(3),n2(3),n3(3),xyz_(natoms,3)
    real(ark)             :: re1(1:3),re2(1:3),alphae(1:3),e0(3),costau, &
                             beta1(1:3),beta2(1:3),y(molec%ncoords,1:3), mu(3), xi(molec%ncoords), tau_,sintau,r0,tau1,tau2,tau_sign
    !
    integer(ik),parameter :: lspace = 150
    integer(ik) :: ierror,rank0,info
    real(rk)    :: dip_rk(3, 1), tmat_rk(3, 3), tsing(3), wspace(lspace),tol = -1.0d-12
    real(ark)    :: cosalpha1,cosalpha2,x1(3),x2(3),x3(3),n30,v2(3),v3(3),v20,v30,n0(3),e20,e30
    character(len=cl)  :: txt
    !
    !
    txt = 'Error: MLdms_HCCH_MB'
    !
    !xyz_(1,1)=-0.017734242
    !xyz_(1,2)=-0.0183397
    !xyz_(1,3)=-1.134903609
    !xyz_(2,1)=0.041027921
    !xyz_(2,2)=0.006163171
    !xyz_(2,3)=1.13187384
    !xyz_(3,1)=-0.013711712
    !xyz_(3,2)=0.125514884
    !xyz_(3,3)=-3.132837078
    !xyz_(4,1)=-0.263864708
    !xyz_(4,2)=0.019585304
    !xyz_(4,3)=3.168940968
    !-0.06208661    0.03255714    0.00256929     -77.209450006   334.7231997638  0.9999997247    14576   1.2     1.06    1.09    175     170     110
    !
    x1(:) = xyz(1,:)-xyz(2,:)
    x2(:) = xyz(3,:)-xyz(1,:)
    x3(:) = xyz(4,:)-xyz(2,:)
    !
    r1 = sqrt(sum(x1(:)**2))
    r2 = sqrt(sum(x2(:)**2))
    r3 = sqrt(sum(x3(:)**2))
    !
    cosalpha1 = sum(-x1(:)*x2(:))/(r1*r2)
    cosalpha2 = sum( x1(:)*x3(:))/(r1*r3)
    !
    alpha1 = aacos(cosalpha1,txt)
    alpha2 = aacos(cosalpha2,txt)
    !
    n3(:) = -x1(:)/r1
    !
    n30 = sum(n3*n3)
    !
    n3(:) = n3/sqrt(n30)
    !
    v2(:) = MLvector_product(x2(:),n3(:))
    !
    v3(:) = MLvector_product(n3(:),x3(:))
    !
    v20 = sum(v2*v2)
    v30 = sum(v3*v3)

    if (abs(r(4))>1e-4   ) then
      continue
    endif


    if (abs(r(5))>1e-4   ) then
      continue
    endif

    if (abs(r(6))>1e-4  ) then
      continue
    endif


    if (abs(r(4))>1e-4 .and. abs(r(6))>1e-4  ) then
      continue
    endif


    if (abs(r(5))>1e-4 .and. abs(r(7))>1e-4  ) then
      continue
    endif

    !
    !setting tau
    if (v20>sqrt(small_).and.v30>sqrt(small_)) then !both v2 and v3 defined
       !
       v2 = v2/sqrt(v20)
       v3 = v3/sqrt(v30)
       !
       costau =-sum(v2*v3)
       !
       tau = aacos(costau,txt)
       !
    elseif (v20<sqrt(small_).and.v30>sqrt(small_)) then !v2 undefined
       !
       v3 = v3/sqrt(v30)
       tau=0
       !
    elseif (v20>sqrt(small_).and.v30<sqrt(small_)) then !v3 undefined
       v2 = v2/sqrt(v20)
       tau=0
       !
    else !both undefined
       !
       tau = 0
       !
    endif
    !End of setting tau
    !
    e2(:) = -MLvector_product(v2(:),n3(:))
    !
    e3(:) = MLvector_product(v3(:),n3(:))
    !
    e20 = sum(e2*e2)
    e30 = sum(e3*e3)
    !
    n1 = 0
    n2 = 0
    !
    !setting n1 and n2
    if (abs(tau-pi)<sqrt(small_).or.abs(tau+pi)<sqrt(small_)) then !tau=pi planar
      !
      n2(:)=-e2(:)
      !
    elseif(e30>small_.and.e20>small_) then
      !
      n1(:) =  (e2(:) + e3(:))
      n2(:) =  (e2(:) - e3(:))
      !
    elseif (e20<small_.and.e30>small_) then !v3 defined
      !
      n1(:) = e3(:)
      !
    elseif (e30<small_.and.e20>small_) then !v2 defined
      !
      n1(:) = e2(:)
      !
    elseif (e30<small_.and.e20<small_) then  !both v2 and v3 undefined
      !
      !define arbitrary direction for n1
      !
      n1 = (/ 1.0_ark,0.0_ark,0.0_ark /)
      !
    else 
      !
      stop 'MLdms_HCCH_MB: you should not be here'
      !
    endif
    !
    !if (e2(2)>0.or.e3(2)>0) then
    !  n1 = -n1
    !endif
    !
    !if (v2(2)>0.or.v3(2)>0) then
    !  n2 = -n2
    !endif
    !
    if (sum(n1(:)*n1(:))>small_) then !n1 defined
      !
      n1(:) = n1(:)/SQRT(sum(n1(:)*n1(:)))
      n2(:) = MLvector_product(n3(:),n1(:))
      !
    elseif (sum(n2(:)**2)>small_) then !n2 defined
      !
      n2(:) = n2(:)/SQRT(sum(n2(:)*n2(:)))
      n1(:) = MLvector_product(n2(:),n3(:))
      !
    else !both n1 and n2 are still undefined
      !
      write(6,"(' both n1 and n2 are  0 : n1 = ',3g12.5,' n2 = ',3g12.5)") n1,n2
      stop ' n1 = n2 = 0 '
      !
    endif
    !
    tmat(1,:) = n1(:)
    tmat(2,:) = n2(:)
    tmat(3,:) = n3(:)
    !
    tau_ = tau
    !
    re1(1:3)     = extF%coef(1,1:3)
    re2(1:3)     = extF%coef(2,1:3)
    alphae(1:3) = extF%coef(3,1:3)/rad
    !
    beta1(1:3)   = extF%coef(4,1:3)
    beta2(1:3)   = extF%coef(5,1:3)
    !
    y(1,:) = (r1 - re1(:)) * exp(-beta1(:) * (r1 - re1(:)) ** 2)
    y(2,:) = (r2 - re2(:)) * exp(-beta2(:) * (r2 - re2(:)) ** 2)
    y(3,:) = (r3 - re2(:)) * exp(-beta2(:) * (r3 - re2(:)) ** 2)
    !
    y(4,:) = (alphae(:) - alpha1)
    y(5,:) = (alphae(:) - alpha2)
    !
    !y(4,:) = sin(alphae(:)) - sin(alpha1)
    !y(5,:) = sin(alphae(:)) - sin(alpha2)
    !
    y(6,:) = cos(tau_)
    !
    mu=0
    mu_t=0
    !
    do imu = 1, 3
       !
       do iterm =  7, extF%nterms(imu)
          !
          ind(1:6) = extF%term(1:6, iterm, imu)
          xi(1:6) = y(1:6,imu) ** ind(1:6)
          !
          mu_t = product(xi(1:6))
          !
          if (ind(2)/=ind(3).or.ind(4)/=ind(5)) then
            !
            ind(2) = extF%term(3, iterm, imu)
            ind(3) = extF%term(2, iterm, imu)
            ind(4) = extF%term(5, iterm, imu)
            ind(5) = extF%term(4, iterm, imu)
            !
            xi(2:5) = y(2:5,imu) ** ind(2:5)
            !
            f_t = 1.0_ark
            if (imu/=1)  f_t = -1.0_ark
            !
            mu_t = mu_t + f_t*product(xi(1:6))
            !
            !
          endif
          !
          !if (imu==1) mu(imu) = mu(1)*cos(tau_*0.5_ark)
          !if (imu==2) mu(imu) = mu(2)*sin(tau_*0.5_ark)
          !
          mu(imu) = mu(imu) + extF%coef(iterm, imu)*mu_t
          ! 
          !if (imu==1) mu(imu) = mu(1)*cos(tau_*0.5_ark)
          !if (imu==2) mu(imu) = mu(2)*sin(tau_*0.5_ark)
          !
       end do
       !
    end do
    !
    mu(1) = mu(1)*cos(tau_*0.5_ark)
    mu(2) = mu(2)*sin(tau_*0.5_ark)
    !
    !write(out,"('tmat start')")
    !
    !if (verbose>=6) write(out,"('tmat')")
    !if (verbose>=6) write(out,*) tmat

    tmat = transpose(tmat)

    !if (verbose>=6) write(out,"('tmat trans')")
    !if (verbose>=6) write(out,*) tmat
    !if (verbose>=6) write(out,"('trans tmat')") tmat
    !tmat = inv(tmat)
    !if (verbose>=6) write(out,"('inv tmat')") tmat
    !
    dimen=3
    !
    !call MLinvmatark(tmat(1:dimen,1:dimen),tmatout(1:dimen,1:dimen),dimen,info)
    !
    f(1:3) = matmul(tmat,mu)
    !
    continue 
    !
    if (abs(r(4))>1e-4   ) then
      continue
    endif


    if (abs(r(5))>1e-4   ) then
      continue
    endif

    if (abs(r(6))>1e-4  ) then
      continue
    endif


    if (abs(r(4))>1e-4 .and. abs(r(6))>1e-4  ) then
      continue
    endif


    if (abs(r(5))>1e-4 .and. abs(r(7))>1e-4  ) then
      continue
    endif
    !
    !if (verbose>=6) write(out,"('f is tmat times mu')") f
    !if (verbose>=6) write(out,"('mu')") mu
    !
    !f(1:2) = f(2:1:-1)
    !
    end subroutine MLdms_HCCH_MB





function MLpoten_v_c2h2_katy(ncoords,natoms,local,xyz,force) result(f)
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
 !
  real(ark)            :: s1,s2,s3,s4,s5,s6,v_c2h2
 real(ark)            :: rech,recc,rehh,a1,a2,a3,a11,a22,a33,a44,a55,a66
 real(ark)            :: a12,a13,a23,a166,a266,a344,a355,a366,a6666,vi
 real(ark)            :: gamma1,gamma2,gamma3
 real(ark)            :: rc1h1,rc2h2,rc1c2,rc2h1,rc1h2,rh1h2
 !ch2
 real(ark)            :: s1_ch2,s2_ch2,s3_ch2,v_ch2
 real(ark)            :: rech_ch2,rehh_ch2,a1_ch2,a2_ch2,a3_ch2,a11_ch2,a22_ch2,a33_ch2
 real(ark)            :: a12_ch2,a13_ch2,a23_ch2,vi_ch2,gamma1_ch2,gamma2_ch2,gamma3_ch2
 real(ark)            :: rch1_ch2,rch2_ch2,rhh_ch2
 !c2h
 real(ark)            :: s1_c2h,s2_c2h,s3_c2h,v_c2h
 real(ark)            :: rehc_c2h,recc_c2h,a1_c2h,a2_c2h,a3_c2h,a11_c2h,a22_c2h,a33_c2h
 real(ark)            :: a12_c2h,a13_c2h,a23_c2h,vi_c2h,gamma1_c2h
 real(ark)            :: rhc_c2h,rcc_c2h,rhcd_c2h
 !h2
 real(ark)            :: s1_h2,v_h2
 real(ark)            :: re_h2,de_h2,a1_h2,a2_h2,a3_h2
 real(ark)            :: rloc_h2,rho_h2
 !c2
 real(ark)            :: s1_c2,v_c2
 real(ark)            :: re_c2,de_c2,a1_c2,a2_c2,a3_c2
 real(ark)            :: rloc_c2,rho_c2
 !ch
 real(ark)            :: s1_ch1,v_ch1,rloc_ch1,rloc_ch2
 real(ark)            :: re_ch,de_ch,a1_ch,a2_ch,a3_ch
 real(ark)            :: rloc_ch,rho_ch,s_ch_1,s_ch_2,v_ch_1,v_ch_2
 
 real(ark)    :: x(4,3),delta1x,delta1y,delta2x,delta2y,v0
 integer(ik)  :: Nangles 
 !
 if (verbose>=6) write(out,"('MLpoten_v_c2h2_katy/start')")
 !
 Nangles = molec%Nangles
 !
 rc1c2    = local(1)
 rc1h1    = local(2)
 rc2h2    = local(3)
 !
 if (Nangles==0.and.molec%Ndihedrals>3) then 
    !
    delta1x = local(4)
    delta1y = local(5)
    !
    delta2x = local(6)
    delta2y = local(7)
    !
    x = 0 
    !
    x(2,1) = 0
    x(2,2) = 0
    x(2,3) = rc1c2
    !
    x(3,1) = -rc1h1*sin(delta1y)
    x(3,2) =  rc1h1*sin(delta1x)
    x(3,3) = -rc1h1*sqrt(cos(delta1x)**2-sin(delta1y)**2)
    !
    x(4,1) = -rc2h2*sin(delta2y)
    x(4,2) =  rc2h2*sin(delta2x)
    x(4,3) =  rc1c2+rc2h2*sqrt(cos(delta2x)**2-sin(delta2y)**2)
    !
    rc1h2 = sqrt( sum( (x(1,:)-x(4,:))**2 ) )
    rc2h1 = sqrt( sum( (x(2,:)-x(3,:))**2 ) )
    rh1h2 = sqrt( sum( (x(3,:)-x(4,:))**2 ) )
   !
 else 
   !
   write(out,"('MLpoten_v_c2h2_katy: only designed for 7 coordinates, not ',i6)") ncoords
   stop 'only designed for 7 coordinates ' 
   !
 endif
 !
 rech      = force(1)
 recc      = force(2)
 rehh      = force(3)
 a1        = force(4)
 a2        = force(5)
 a3        = force(6)
 a11       = force(7)
 a22       = force(8) 
 a33       = force(9)
 a44       = force(10)
 a55       = force(11)
 a66       = force(12)
 a12       = force(13)
 a13       = force(14)
 a23       = force(15)
 a166      = force(16)
 a266      = force(17)
 a344      = force(18)
 a355      = force(19)
 a366      = force(20)
 a6666     = force(21)
 vi        = force(22)
 gamma1    = force(23)
 gamma2    = force(24)
 gamma3    = force(25)
!ch2
 rech_ch2   = force(26)
 rehh_ch2   = force(27)
 a1_ch2     = force(28)
 a2_ch2     = force(29)
 a3_ch2     = force(30)
 a11_ch2    = force(31)
 a22_ch2    = force(32)
 a33_ch2    = force(33)
 a12_ch2    = force(34)
 a13_ch2    = force(35)
 a23_ch2    = force(36)
 vi_ch2     = force(37)
 gamma1_ch2 = force(38)
 gamma2_ch2 = force(39)
 gamma3_ch2 = force(40)
!c2h
 rehc_c2h   = force(41)
 recc_c2h   = force(42)
 a1_c2h     = force(43)
 a2_c2h     = force(44)
 a3_c2h     = force(45)
 a11_c2h    = force(46)
 a22_c2h    = force(47)
 a33_c2h    = force(48)
 a12_c2h    = force(49)
 a13_c2h    = force(50)
 a23_c2h    = force(51)
 vi_c2h     = force(52)
 gamma1_c2h = force(53)
!h2
 a1_h2   = force(54)
 a2_h2   = force(55)
 a3_h2   = force(56)
 de_h2   = force(57)
 re_h2   = force(58)
!c2
 a1_c2   = force(59)
 a2_c2   = force(60)
 a3_c2   = force(61) 
 de_c2   = force(62)
 re_c2   = force(63) 
!ch
 a1_ch   = force(64)
 a2_ch   = force(65)
 a3_ch   = force(66) 
 de_ch   = force(67)
 re_ch   = force(68) 

!ch2
 !rch1_ch2 = local(7)
 !rch2_ch2 = local(8)
 !rhh_ch2  = local(9)
!c2h
 !rhc_c2h  = local(10)
 !rcc_c2h  = local(11)
 !rhcd_c2h = local(12)
!h2
 rloc_h2  = rh1h2
!c2
 rloc_c2  = rc1c2 
!ch
 rloc_ch1  = rc1h1
 rloc_ch2  = rc2h2

!c2h2
         s1 = rc1h1+rc2h2+rc2h1+rc1h2-(4.0_ark*rech)
         s2 = rc1c2-recc
         s3 = rh1h2-rehh
         s4 = rc1h1+rc2h2-rc2h1-rc1h2
         s5 = rc1h1-rc2h2+rc2h1-rc1h2
         s6 = rc1h1-rc2h2-rc2h1+rc1h2
!ch2
     !    s1_ch2 = rch1_ch2-rech_ch2
     !    s2_ch2 = rch2_ch2-rech_ch2
     !    s3_ch2 = rhh_ch2-rehh_ch2
!c2h
     !    s1_c2h = rhc_c2h-rehc_c2h
     !    s2_c2h = rcc_c2h-recc_c2h
     !    s3_c2h = rhcd_c2h-rehc_c2h
!h2
         s1_h2 = rloc_h2-re_h2
!c2
         s1_c2 = rloc_c2-re_c2
!ch
         s_ch_1 = rloc_ch1-re_ch
         s_ch_2 = rloc_ch2-re_ch

!c2h2
          v_c2h2 = 0
    !     v_c2h2 = (1.0_ark-tanh(gamma1*s1*0.25_ark) )*(1.0_ark-tanh(gamma2*s2*0.5_ark))*(1.0_ark-tanh(gamma3*s3*0.5_ark))* &
    !      vi*(1+(a1*s1)+(a2*s2)+(a3*s3)+(a11*s1*s1)+(a22*s2*s2)+(a33*s3*s3)+(a44*s4*s4)+(a55*s5*s5)+(a66*s6*s6)+ &
    !      (a12*s1*s2)+(a13*s1*s3)+(a23*s2*s3)+(a166*s1*s6*s6)+(a266*s2*s6*s6)+(a344*s3*s4*s4)+(a355*s3*s5*s5)+(a366*s3*s6*s6)+(a6666*s6**6))

!ch2
          v_ch2 = 0
    !     v_ch2 = (1.0_ark-tanh(gamma1_ch2*s1_ch2*0.5_ark))*(1.0_ark-tanh(gamma2_ch2*s2_ch2*0.5_ark))*(1.0_ark-tanh(gamma3_ch2*s3_ch2*0.5_ark))* &
    !      vi_ch2*(1+(a1_ch2*s1_ch2)+(a2_ch2*s2_ch2)+(a3_ch2*s3_ch2)+(a11_ch2*s1_ch2*s1_ch2)+(a22_ch2*s2_ch2*s2_ch2)+(a33_ch2*s3_ch2*s3_ch2)+(a12_ch2*s1_ch2*s2_ch2)+(a13_ch2*s1_ch2*s3_ch2)+(a23_ch2*a2_ch2*a3_ch2))

!c2h
          v_c2h= 0
    !     v_c2h = vi_c2h*(1.0_ark+(a1_c2h*s1_c2h)+(a2_c2h*s2_c2h)+(a3_c2h*s3_c2h)+(a11_c2h*s1_c2h*s1_c2h)+(a22_c2h*s2_c2h*s2_c2h)+(a33_c2h*s3_c2h*s3_c2h)+(a12_c2h*s1_c2h*s2_c2h)+(a13_c2h*s1_c2h*s3_c2h)+(a23_c2h*s2_c2h*s3_c2h))*(1.0_ark-tanh(gamma1_c2h*(s1_c2h+s2_c2h+s3_c2h)*0.5_ark/sqrt3))

!h2
         v_h2 = (-de_h2*(1.0_ark+(a1_h2*s1_h2)+(a2_h2*s1_h2**2)+(a3_h2*s1_h2**3)*exp(-a1_h2*s1_h2)))

!c2
         v_c2 = (-de_c2*(1.0_ark+(a1_c2*s1_c2)+(a2_c2*s1_c2**2)+(a3_c2*a1_c2**3)*exp(-a1_c2*s1_c2)))

!ch
         v_ch_1 = (-de_ch*(1.0_ark+(a1_ch*s_ch_1)+(a2_ch*s_ch_1**2)+(a3_ch*s_ch_1**3)*exp(-a1_ch*s_ch_1)))
         v_ch_2 = (-de_ch*(1.0_ark+(a1_ch*s_ch_2)+(a2_ch*s_ch_2**2)+(a3_ch*s_ch_2**3)*exp(-a1_ch*s_ch_2)))

!full pot function

         f = v_h2+v_c2+4.0_ark*(v_ch_1+v_ch_2)+(2.0_ark*v_c2h)+(2.0_ark*v_ch2)+v_c2h2
         !
         v0 = force(1)
         !
         f = f - v0
         !
         f = f*1.0e-11/planck/vellgt

 if (verbose>=6) write(out,"('MLpoten_v_c2h2_katy/end')")

 end function MLpoten_v_c2h2_katy


function MLpoten_v_c2h2_mlt(ncoords,natoms,local,xyz,force) result(f)
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
 !ch
 real(ark)            :: s1,s2,s3,s4,s5,s6,s7
 real(ark)            :: f11,f22,f33,f44,f55,f66,f77,f12
 real(ark)            :: rc1h1,rc2h2,rc1c2,v0,rech,recc

 real(ark)    :: x(4,3),theta1x,theta1y,theta2x,theta2y
 integer(ik)  :: Nangles
 !
 if (verbose>=6) write(out,"('MLpoten_v_c2h2_mlt/start')")
 !
 Nangles = molec%Nangles
 !
 rc1c2    = local(1)
 rc1h1    = local(2)
 rc2h2    = local(3)
 !
 if (Nangles==0.and.molec%Ndihedrals>3) then
    !
    theta1x = local(4)
    theta1y = local(5)
    !
    theta2x = local(6)
    theta2y = local(7)
    !
    x = 0
    !
    x(2,1) = 0
    x(2,2) = 0
    x(2,3) = rc1c2
    !
    x(3,1) = -rc1h1*sin(theta1x)
    x(3,2) =  rc1h1*sin(theta1y)
    x(3,3) = -(sqrt((rc1h1**2)-((rc1h1*sin(theta1x))**2)-((rc1h1*sin(theta1y))**2)))
    !
    x(4,1) = -rc2h2*sin(theta2x)
    x(4,2) =  rc2h2*sin(theta2y)
    x(4,3) =  rc1c2+(sqrt((rc2h2**2)-((rc2h2*sin(theta2x))**2)-((rc2h2*sin(theta2y))**2)))
   !
 else
   !
   write(out,"('MLpoten_v_c2h2_mlt: only designed for 7 coordinates, not ',i6)") ncoords
   stop 'only designed for 7 coordinates '
   !
 endif
 f11   = force(2)
 f22   = force(3)
 f33   = force(4)
 f44   = force(5)
 f55   = force(6)
 f66   = force(7)
 f77   = force(8)
 f12   = force(9)
 rech  = force(10)
 recc  = force(11)

 s1=(rc1h1+rc2h2-2.0_ark*rech)/(sqrt(2.0_ark))
 s2=rc1c2-recc
 s3=(rc1h1-rc2h2)/(sqrt(2.0_ark))
 s4=(theta1x-theta2x)/(sqrt(2.0_ark))
 s5=(theta1y-theta2y)/(sqrt(2.0_ark))
 s6=(theta1x+theta2x)/(sqrt(2.0_ark))
 s7=(theta1y+theta2y)/(sqrt(2.0_ark))

!full pot function

         f = 0.5_ark*((f11*s1*s1)+(f22*s2*s2)+(f33*s3*s3)+(f44*s4*s4)+(f55*s5*s5)+(f66*s6*s6)+(f77*s7*s7)+(f12*((s1*s2)+(s1*s2))))
         !f = 0.5_ark*((f11*s1*s1)+(f33*s3*s3))
         !
         v0 = force(1)
         !
         f = f + v0
         !
        f = f*1.0e-11/planck/vellgt
        !
        !f = theta1x
        !

 if (verbose>=6) write(out,"('MLpoten_v_c2h2_mlt/end')")

 end function MLpoten_v_c2h2_mlt


subroutine potC2H2_diff_V(n,y1,y2,y3,y4,y5,y6,y7,dF)
    !
    implicit none
    !
    integer(ik),intent(in)  :: n
    real(ark),intent(in)  :: y1,y2,y3,y4,y5,y6,y7
    real(ark),intent(out) :: dF(n)

      !
      dF(1:4) = 0
      !
      !dF(1) = y4
      !dF(2) = y5
      !dF(3) = y6
      !dF(4) = y7
      !
      dF(5) = 1.0_ark
      dF(6) = y3+y2
      dF(7) = y1
      dF(8) = y6**2+y5**2+y4**2+y7**2  
      dF(9) = y5*y7+y4*y6
      dF(10) = y2**2+y3**2
      dF(11) = y2*y3
      dF(12) = (y3+y2)*y1
      dF(13) = y1**2
      dF(14) = (y7**2+y6**2)*y2+(y4**2+y5**2)*y3
      dF(15) = (y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3
      dF(16) = (y4**2+y5**2)*y2+(y7**2+y6**2)*y3
      dF(17) = y2*y3**2+y2**2*y3
      dF(18) = y2**3+y3**3
      dF(19) = (y5*y7+y4*y6)*y1
      dF(20) = (y7**2+y5**2+y6**2+y4**2)*y1
      dF(21) = y1*y2*y3
      dF(22) = (y3**2+y2**2)*y1
      dF(23) = (y3+y2)*y1**2
      dF(24) = y1**3
      dF(25) = y7**4+y5**4+y6**4+y4**4
      dF(26) = y4*y5*y6*y7
      dF(27) = y4**2*y7**2+y5**2*y6**2
      dF(28) = y4**2*y6**2+y5**2*y7**2
      dF(29) = y4**2*y5*y7+(y5**2*y6+y6*y7**2)*y4+y5*y6**2*y7
      dF(30) = y4**2*y5**2+y6**2*y7**2
      dF(31) = y4**3*y6+y5*y7**3+y4*y6**3+y5**3*y7
      dF(32) = (y4**2+y5**2)*y2**2+(y6**2+y7**2)*y3**2
      dF(33) = (y4*y6+y5*y7)*y2**2+(y4*y6+y5*y7)*y3**2
      dF(34) = (y6**2+y7**2)*y2**2+(y4**2+y5**2)*y3**2
      dF(35) = y3**4+y2**4
      dF(36) = (y4*y6+y5*y7)*y3*y2
      dF(37) = (y4**2+y5**2+y6**2+y7**2)*y3*y2
      dF(38) = y2**2*y3**2
      dF(39) = y2**3*y3+y2*y3**3
      dF(40) = ((y6**2+y7**2)*y2+(y4**2+y5**2)*y3)*y1
      dF(41) = ((y4*y6+y5*y7)*y2+(y4*y6+y5*y7)*y3)*y1
      dF(42) = ((y4**2+y5**2)*y2+(y6**2+y7**2)*y3)*y1
      dF(43) = (y2**2*y3+y2*y3**2)*y1
      dF(44) = (y2**3+y3**3)*y1
      dF(45) = (y4*y6+y5*y7)*y1**2
      dF(46) = (y4**2+y5**2+y6**2+y7**2)*y1**2
      dF(47) = y1**2*y2*y3
      dF(48) = (y2**2+y3**2)*y1**2
      dF(49) = (y2+y3)*y1**3
      dF(50) = y1**4
      dF(51) = y2*y4**2*y5**2+y3*y6**2*y7**2
      dF(52) = (y5**4+y4**4)*y2+(y7**4+y6**4)*y3
      dF(53) = (y4**3*y6+y5**3*y7)*y2+(y4*y6**3+y5*y7**3)*y3
      dF(54) = (y4**2*y5*y7+y4*y5**2*y6)*y2+(y4*y6*y7**2+y5*y6**2*y7)*y3
      dF(55) = (y4**2*y6**2+y5**2*y7**2)*y2+(y4**2*y6**2+y5**2*y7**2)*y3
      dF(56) = (y5**2*y6**2+y4**2*y7**2)*y2+(y5**2*y6**2+y4**2*y7**2)*y3
      dF(57) = (y4**2+y5**2)*y2**3+(y6**2+y7**2)*y3**3
      dF(58) = (y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3
      dF(59) = (y6**2+y7**2)*y2**3+(y4**2+y5**2)*y3**3
      dF(60) = y3**5+y2**5
      dF(61) = y2*y6**2*y7**2+y3*y4**2*y5**2
      dF(62) = (y7**4+y6**4)*y2+(y5**4+y4**4)*y3
      dF(63) = (y4*y6*y7**2+y5*y6**2*y7)*y2+(y4**2*y5*y7+y4*y5**2*y6)*y3
      dF(64) = (y4*y6**3+y5*y7**3)*y2+(y4**3*y6+y5**3*y7)*y3
      dF(65) = y2*y4*y5*y6*y7+y3*y4*y5*y6*y7
      dF(66) = (y4**2+y5**2)*y3*y2**2+(y6**2+y7**2)*y3**2*y2
      dF(67) = (y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2
      dF(68) = (y6**2+y7**2)*y3*y2**2+(y4**2+y5**2)*y3**2*y2
      dF(69) = y2*y3**4+y2**4*y3
      dF(70) = y2**3*y3**2+y2**2*y3**3
      dF(71) = y1*y4*y5*y6*y7
      dF(72) = (y5**2*y6**2+y4**2*y7**2)*y1
      dF(73) = (y4**2*y6**2+y5**2*y7**2)*y1
      dF(74) = (y4**2*y5*y7+(y6*y7**2+y5**2*y6)*y4+y5*y6**2*y7)*y1
      dF(75) = (y6**2*y7**2+y4**2*y5**2)*y1
      dF(76) = (y5*y7**3+y4*y6**3+y5**3*y7+y4**3*y6)*y1
      dF(77) = (y5**4+y7**4+y6**4+y4**4)*y1
      dF(78) = (y5*y7+y4*y6)*y3*y2*y1
      dF(79) = (y5**2+y7**2+y6**2+y4**2)*y3*y2*y1
      dF(80) = ((y6**2+y7**2)*y2**2+(y4**2+y5**2)*y3**2)*y1
      dF(81) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1
      dF(82) = ((y4**2+y5**2)*y2**2+(y6**2+y7**2)*y3**2)*y1
      dF(83) = y1*y2**2*y3**2
      dF(84) = (y2**3*y3+y2*y3**3)*y1
      dF(85) = (y2**4+y3**4)*y1
      dF(86) = ((y6**2+y7**2)*y2+(y4**2+y5**2)*y3)*y1**2
      dF(87) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**2
      dF(88) = ((y4**2+y5**2)*y2+(y6**2+y7**2)*y3)*y1**2
      dF(89) = (y2**2*y3+y2*y3**2)*y1**2
      dF(90) = (y2**3+y3**3)*y1**2
      dF(91) = (y5*y7+y4*y6)*y1**3
      dF(92) = (y5**2+y7**2+y6**2+y4**2)*y1**3
      dF(93) = y1**3*y2*y3
      dF(94) = (y3**2+y2**2)*y1**3
      dF(95) = (y2+y3)*y1**4
      dF(96) = y1**5
      dF(97) = y7**6+y4**6+y6**6+y5**6
      dF(98) = y6**2*y7**4+y4**2*y5**4+y4**4*y5**2+y6**4*y7**2
      dF(99) = y5*y6**2*y7**3+y4**2*y5**3*y7+y4**3*y5**2*y6+y4*y6**3*y7**2
      dF(100) = y5**2*y6**4+y4**4*y7**2+y5**4*y6**2+y4**2*y7**4
      dF(101) = y5**3*y7**3+y4**3*y6**3
      dF(102) = y5**3*y6**2*y7+y4**2*y5*y7**3+y4**3*y6*y7**2+y4*y5**2*y6**3
      dF(103) = y5**4*y7**2+y4**2*y6**4+y4**4*y6**2+y5**2*y7**4
      dF(104) = y5**5*y7+y4**5*y6+y4*y6**5+y5*y7**5
      dF(105) = y4*y5**2*y6*y7**2+y4**2*y5*y6**2*y7
      dF(106) = y4**4*y5*y7+(y5**4*y6+y6*y7**4)*y4+y5*y6**4*y7
      dF(107) = ((y6**2+y7**2)*y5**2+y6**2*y7**2)*y4**2+y5**2*y6**2*y7**2
      dF(108) = y4**3*y5*y6*y7+(y5**3*y6*y7+(y6*y7**3+y6**3*y7)*y5)*y4
      dF(109) = (y4**2+y5**2)*y2**4+(y6**2+y7**2)*y3**4
      dF(110) = (y6**2+y7**2)*y2**4+(y4**2+y5**2)*y3**4
      dF(111) = (y5*y7+y4*y6)*y2**4+(y5*y7+y4*y6)*y3**4
      dF(112) = y3**6+y2**6
      dF(113) = (y7**4+y6**4+y4**4+y5**4)*y3*y2
      dF(114) = (y6**2*y7**2+y4**2*y5**2)*y3*y2
      dF(115) = (y5*y7**3+y4*y6**3+y4**3*y6+y5**3*y7)*y3*y2
      dF(116) = (y4**2*y5*y7+(y6*y7**2+y5**2*y6)*y4+y5*y6**2*y7)*y3*y2
      dF(117) = (y5**2*y7**2+y4**2*y6**2)*y3*y2
      dF(118) = (y5**2*y6**2+y4**2*y7**2)*y3*y2
      dF(119) = y2*y3*y4*y5*y6*y7
      dF(120) = (y4**2+y5**2)*y3*y2**3+(y6**2+y7**2)*y3**3*y2
      dF(121) = (y5*y7+y4*y6)*y3*y2**3+(y5*y7+y4*y6)*y3**3*y2
      dF(122) = (y6**2+y7**2)*y3*y2**3+(y4**2+y5**2)*y3**3*y2
      dF(123) = y2*y3**5+y2**5*y3
      dF(124) = y2**2*y6**2*y7**2+y3**2*y4**2*y5**2
      dF(125) = (y7**4+y6**4)*y2**2+(y4**4+y5**4)*y3**2
      dF(126) = (y5*y6**2*y7+y4*y6*y7**2)*y2**2+(y4**2*y5*y7+y4*y5**2*y6)*y3**2
      dF(127) = (y4*y6**3+y5*y7**3)*y2**2+(y5**3*y7+y4**3*y6)*y3**2
      dF(128) = y2**2*y4*y5*y6*y7+y3**2*y4*y5*y6*y7
      dF(129) = (y5**2*y6**2+y4**2*y7**2)*y2**2+(y5**2*y6**2+y4**2*y7**2)*y3**2
      dF(130) = (y5**2*y7**2+y4**2*y6**2)*y2**2+(y5**2*y7**2+y4**2*y6**2)*y3**2
      dF(131) = (y4**2*y5*y7+y4*y5**2*y6)*y2**2+(y5*y6**2*y7+y4*y6*y7**2)*y3**2
      dF(132) = y2**2*y4**2*y5**2+y3**2*y6**2*y7**2
      dF(133) = (y5**3*y7+y4**3*y6)*y2**2+(y4*y6**3+y5*y7**3)*y3**2
      dF(134) = (y4**4+y5**4)*y2**2+(y7**4+y6**4)*y3**2
      dF(135) = (y5*y7+y4*y6)*y3**2*y2**2
      dF(136) = (y4**2+y7**2+y6**2+y5**2)*y3**2*y2**2
      dF(137) = y2**2*y3**4+y2**4*y3**2
      dF(138) = y2**3*y3**3
      dF(139) = (y2*y6**2*y7**2+y3*y4**2*y5**2)*y1
      dF(140) = ((y7**4+y6**4)*y2+(y4**4+y5**4)*y3)*y1
      dF(141) = ((y5*y6**2*y7+y4*y6*y7**2)*y2+(y4**2*y5*y7+y4*y5**2*y6)*y3)*y1
      dF(142) = ((y4*y6**3+y5*y7**3)*y2+(y5**3*y7+y4**3*y6)*y3)*y1
      dF(143) = (y2*y4*y5*y6*y7+y3*y4*y5*y6*y7)*y1
      dF(144) = ((y5**2*y6**2+y4**2*y7**2)*y2+(y5**2*y6**2+y4**2*y7**2)*y3)*y1
      dF(145) = ((y5**2*y7**2+y4**2*y6**2)*y2+(y5**2*y7**2+y4**2*y6**2)*y3)*y1
      dF(146) = ((y4**2*y5*y7+y4*y5**2*y6)*y2+(y5*y6**2*y7+y4*y6*y7**2)*y3)*y1
      dF(147) = (y2*y4**2*y5**2+y3*y6**2*y7**2)*y1
      dF(148) = ((y5**3*y7+y4**3*y6)*y2+(y4*y6**3+y5*y7**3)*y3)*y1
      dF(149) = ((y4**4+y5**4)*y2+(y7**4+y6**4)*y3)*y1
      dF(150) = ((y6**2+y7**2)*y3*y2**2+(y4**2+y5**2)*y3**2*y2)*y1
      dF(151) = ((y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2)*y1
      dF(152) = ((y4**2+y5**2)*y3*y2**2+(y6**2+y7**2)*y3**2*y2)*y1
      dF(153) = ((y6**2+y7**2)*y2**3+(y4**2+y5**2)*y3**3)*y1
      dF(154) = ((y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3)*y1
      dF(155) = ((y4**2+y5**2)*y2**3+(y6**2+y7**2)*y3**3)*y1
      dF(156) = (y2**3*y3**2+y2**2*y3**3)*y1
      dF(157) = (y2**4*y3+y2*y3**4)*y1
      dF(158) = (y2**5+y3**5)*y1
      dF(159) = y1**2*y4*y5*y6*y7
      dF(160) = (y5**2*y6**2+y4**2*y7**2)*y1**2
      dF(161) = (y5**2*y7**2+y4**2*y6**2)*y1**2
      dF(162) = (y4**2*y5*y7+(y6*y7**2+y5**2*y6)*y4+y5*y6**2*y7)*y1**2
      dF(163) = (y6**2*y7**2+y4**2*y5**2)*y1**2
      dF(164) = (y5*y7**3+y4*y6**3+y4**3*y6+y5**3*y7)*y1**2
      dF(165) = (y7**4+y6**4+y4**4+y5**4)*y1**2
      dF(166) = (y5*y7+y4*y6)*y3*y2*y1**2
      dF(167) = (y4**2+y7**2+y6**2+y5**2)*y3*y2*y1**2
      dF(168) = ((y6**2+y7**2)*y2**2+(y4**2+y5**2)*y3**2)*y1**2
      dF(169) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1**2
      dF(170) = ((y4**2+y5**2)*y2**2+(y6**2+y7**2)*y3**2)*y1**2
      dF(171) = y1**2*y2**2*y3**2
      dF(172) = (y2*y3**3+y2**3*y3)*y1**2
      dF(173) = (y2**4+y3**4)*y1**2
      dF(174) = (y3**3+y2**3)*y1**3
      dF(175) = ((y6**2+y7**2)*y2+(y4**2+y5**2)*y3)*y1**3
      dF(176) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**3
      dF(177) = ((y4**2+y5**2)*y2+(y6**2+y7**2)*y3)*y1**3
      dF(178) = (y2**2*y3+y2*y3**2)*y1**3
      dF(179) = (y5*y7+y4*y6)*y1**4
      dF(180) = (y4**2+y7**2+y6**2+y5**2)*y1**4
      dF(181) = y1**4*y2*y3
      dF(182) = (y2**2+y3**2)*y1**4
      dF(183) = (y2+y3)*y1**5
      dF(184) = y1**6

      dF(185) = y3**3*y6**2*y7**2+y2**3*y4**2*y5**2
      dF(186) = (y4*y5**2*y6+y4**2*y5*y7)*y2**3+(y4*y6*y7**2+y5*y6**2*y7)*y3**3
      dF(187) = (y5**3*y7+y4**3*y6)*y2**3+(y4*y6**3+y5*y7**3)*y3**3
      dF(188) = y3**3*y4*y5*y6*y7+y2**3*y4*y5*y6*y7
      dF(189) = (y4*y6*y7**2+y5*y6**2*y7)*y2**3+(y4*y5**2*y6+y4**2*y5*y7)*y3**3
      dF(190) = (y5**2*y6**2+y4**2*y7**2)*y2**3+(y5**2*y6**2+y4**2*y7**2)*y3**3
      dF(191) = (y5**2*y7**2+y4**2*y6**2)*y2**3+(y5**2*y7**2+y4**2*y6**2)*y3**3
      dF(192) = y3**3*y4**2*y5**2+y2**3*y6**2*y7**2
      dF(193) = (y4*y6**3+y5*y7**3)*y2**3+(y5**3*y7+y4**3*y6)*y3**3
      dF(194) = (y6**4+y7**4)*y2**3+(y4**4+y5**4)*y3**3
      dF(195) = (y6**2+y7**2)*y2**5+(y5**2+y4**2)*y3**5
      dF(196) = y3**7+y2**7
      dF(197) = (y6**2*y7**4+y6**4*y7**2)*y2+(y4**2*y5**4+y4**4*y5**2)*y3
      dF(198) = (y6**6+y7**6)*y2+(y5**6+y4**6)*y3
      dF(199) = (y4*y6*y7**4+y5*y6**4*y7)*y2+(y4*y5**4*y6+y4**4*y5*y7)*y3
      dF(200) = (y4*y6**3*y7**2+y5*y6**2*y7**3)*y2+(y4**2*y5**3*y7+y4**3*y5**2*y6)*y3
      dF(201) = (y4*y6**5+y5*y7**5)*y2+(y5**5*y7+y4**5*y6)*y3
      dF(202) = (y6*y7**3+y6**3*y7)*y5*y4*y2+(y4*y5**3*y6*y7+y4**3*y5*y6*y7)*y3
      dF(203) = (y5**2*y6**4+y4**2*y7**4)*y2+(y5**4*y6**2+y4**4*y7**2)*y3
      dF(204) = (y4**2*y6**2*y7**2+y5**2*y6**2*y7**2)*y2+(y6**2+y7**2)*y5**2*y4**2*y3
      dF(205) = (y5**2*y7**4+y4**2*y6**4)*y2+(y5**4*y7**2+y4**4*y6**2)*y3
      dF(206) = (y4**2*y5*y7**3+y4*y5**2*y6**3)*y2+(y4**3*y6*y7**2+y5**3*y6**2*y7)*y3
      dF(207) = (y4**2*y5*y6**2*y7+y4*y5**2*y6*y7**2)*y2+(y4**2*y5*y6**2*y7+y4*y5**2*y6*y7**2)*y3
      dF(208) = (y6**2+y7**2)*y5**2*y4**2*y2+(y4**2*y6**2*y7**2+y5**2*y6**2*y7**2)*y3
      dF(209) = (y4**3*y6*y7**2+y5**3*y6**2*y7)*y2+(y4**2*y5*y7**3+y4*y5**2*y6**3)*y3
      dF(210) = (y4**3*y6**3+y5**3*y7**3)*y2+(y4**3*y6**3+y5**3*y7**3)*y3
      dF(211) = (y4*y5**3*y6*y7+y4**3*y5*y6*y7)*y2+(y6*y7**3+y6**3*y7)*y5*y4*y3
      dF(212) = (y4**2*y5**3*y7+y4**3*y5**2*y6)*y2+(y4*y6**3*y7**2+y5*y6**2*y7**3)*y3
      dF(213) = (y5**4*y6**2+y4**4*y7**2)*y2+(y5**2*y6**4+y4**2*y7**4)*y3
      dF(214) = (y5**4*y7**2+y4**4*y6**2)*y2+(y5**2*y7**4+y4**2*y6**4)*y3
      dF(215) = (y4*y5**4*y6+y4**4*y5*y7)*y2+(y4*y6*y7**4+y5*y6**4*y7)*y3
      dF(216) = (y4**2*y5**4+y4**4*y5**2)*y2+(y6**2*y7**4+y6**4*y7**2)*y3
      dF(217) = (y5**5*y7+y4**5*y6)*y2+(y4*y6**5+y5*y7**5)*y3
      dF(218) = (y5**6+y4**6)*y2+(y6**6+y7**6)*y3
      dF(219) = (y5**2+y4**2)*y3*y2**4+(y6**2+y7**2)*y3**4*y2
      dF(220) = (y4*y6+y5*y7)*y3*y2**4+(y4*y6+y5*y7)*y3**4*y2
      dF(221) = (y6**2+y7**2)*y3*y2**4+(y5**2+y4**2)*y3**4*y2
      dF(222) = y2*y3**6+y2**6*y3
      dF(223) = (y6**4+y7**4)*y3*y2**2+(y4**4+y5**4)*y3**2*y2
      dF(224) = y2**2*y3*y6**2*y7**2+y2*y3**2*y4**2*y5**2
      dF(225) = (y4*y6*y7**2+y5*y6**2*y7)*y3*y2**2+(y4*y5**2*y6+y4**2*y5*y7)*y3**2*y2
      dF(226) = (y4*y6**3+y5*y7**3)*y3*y2**2+(y5**3*y7+y4**3*y6)*y3**2*y2
      dF(227) = y2**2*y3*y4*y5*y6*y7+y2*y3**2*y4*y5*y6*y7
      dF(228) = (y4*y5**2*y6+y4**2*y5*y7)*y3*y2**2+(y4*y6*y7**2+y5*y6**2*y7)*y3**2*y2
      dF(229) = (y5**2*y6**2+y4**2*y7**2)*y3*y2**2+(y5**2*y6**2+y4**2*y7**2)*y3**2*y2
      dF(230) = (y5**2*y7**2+y4**2*y6**2)*y3*y2**2+(y5**2*y7**2+y4**2*y6**2)*y3**2*y2
      dF(231) = y2**2*y3*y4**2*y5**2+y2*y3**2*y6**2*y7**2
      dF(232) = (y5**3*y7+y4**3*y6)*y3*y2**2+(y4*y6**3+y5*y7**3)*y3**2*y2
      dF(233) = (y4**4+y5**4)*y3*y2**2+(y6**4+y7**4)*y3**2*y2
      dF(234) = (y5**2+y4**2)*y3**2*y2**3+(y6**2+y7**2)*y3**3*y2**2
      dF(235) = (y4*y6+y5*y7)*y3**2*y2**3+(y4*y6+y5*y7)*y3**3*y2**2
      dF(236) = (y6**2+y7**2)*y3**2*y2**3+(y5**2+y4**2)*y3**3*y2**2
      dF(237) = (y4**4+y5**4)*y2**3+(y6**4+y7**4)*y3**3
      dF(238) = y2**3*y3**4+y2**4*y3**3
      dF(239) = (y5**2+y4**2)*y2**5+(y6**2+y7**2)*y3**5
      dF(240) = (y4*y6+y5*y7)*y2**5+(y4*y6+y5*y7)*y3**5
      dF(241) = y2**5*y3**2+y2**2*y3**5
      dF(242) = (y4**4*y5*y7+(y6*y7**4+y5**4*y6)*y4+y5*y6**4*y7)*y1
      dF(243) = (y4*y6**5+y5**5*y7+y5*y7**5+y4**5*y6)*y1
      dF(244) = (y4**3*y5*y6*y7+(y5**3*y6*y7+(y6*y7**3+y6**3*y7)*y5)*y4)*y1
      dF(245) = (y4**2*y5*y6**2*y7+y4*y5**2*y6*y7**2)*y1
      dF(246) = (((y6**2+y7**2)*y5**2+y6**2*y7**2)*y4**2+y5**2*y6**2*y7**2)*y1
      dF(247) = (y4*y6**3*y7**2+y4**2*y5**3*y7+y5*y6**2*y7**3+y4**3*y5**2*y6)*y1
      dF(248) = (y6**2*y7**4+y6**4*y7**2+y4**2*y5**4+y4**4*y5**2)*y1
      dF(249) = (y4**3*y6*y7**2+y4*y5**2*y6**3+y5**3*y6**2*y7+y4**2*y5*y7**3)*y1
      dF(250) = (y4**3*y6**3+y5**3*y7**3)*y1
      dF(251) = (y5**4*y6**2+y4**2*y7**4+y4**4*y7**2+y5**2*y6**4)*y1
      dF(252) = (y5**4*y7**2+y4**4*y6**2+y5**2*y7**4+y4**2*y6**4)*y1
      dF(253) = (y4**6+y5**6+y6**6+y7**6)*y1
      dF(254) = ((y4*y5**2*y6+y4**2*y5*y7)*y2**2+(y4*y6*y7**2+y5*y6**2*y7)*y3**2)*y1
      dF(255) = (y6**2*y7**2+y4**2*y5**2)*y3*y2*y1
      dF(256) = (y5*y7**3+y4*y6**3+y4**3*y6+y5**3*y7)*y3*y2*y1
      dF(257) = (y4**2*y5*y7+(y5**2*y6+y6*y7**2)*y4+y5*y6**2*y7)*y3*y2*y1
      dF(258) = y1*y2*y3*y4*y5*y6*y7
      dF(259) = (y5**2*y6**2+y4**2*y7**2)*y3*y2*y1
      dF(260) = (y5**2*y7**2+y4**2*y6**2)*y3*y2*y1
      dF(261) = (y4**4+y6**4+y7**4+y5**4)*y3*y2*y1
      dF(262) = (y2**2*y6**2*y7**2+y3**2*y4**2*y5**2)*y1
      dF(263) = ((y6**4+y7**4)*y2**2+(y4**4+y5**4)*y3**2)*y1
      dF(264) = ((y5**2*y6**2+y4**2*y7**2)*y2**2+(y5**2*y6**2+y4**2*y7**2)*y3**2)*y1
      dF(265) = ((y5**3*y7+y4**3*y6)*y2**2+(y4*y6**3+y5*y7**3)*y3**2)*y1
      dF(266) = ((y4**4+y5**4)*y2**2+(y6**4+y7**4)*y3**2)*y1
      dF(267) = ((y4*y6*y7**2+y5*y6**2*y7)*y2**2+(y4*y5**2*y6+y4**2*y5*y7)*y3**2)*y1
      dF(268) = ((y4*y6**3+y5*y7**3)*y2**2+(y5**3*y7+y4**3*y6)*y3**2)*y1
      dF(269) = (y2**2*y4*y5*y6*y7+y3**2*y4*y5*y6*y7)*y1
      dF(270) = ((y5**2*y7**2+y4**2*y6**2)*y2**2+(y5**2*y7**2+y4**2*y6**2)*y3**2)*y1
      dF(271) = (y2**2*y4**2*y5**2+y3**2*y6**2*y7**2)*y1
      dF(272) = (y6**2+y7**2+y5**2+y4**2)*y3**2*y2**2*y1
      dF(273) = (y4*y6+y5*y7)*y3**2*y2**2*y1
      dF(274) = ((y6**2+y7**2)*y3*y2**3+(y5**2+y4**2)*y3**3*y2)*y1
      dF(275) = ((y4*y6+y5*y7)*y3*y2**3+(y4*y6+y5*y7)*y3**3*y2)*y1
      dF(276) = ((y5**2+y4**2)*y3*y2**3+(y6**2+y7**2)*y3**3*y2)*y1
      dF(277) = y1*y2**3*y3**3
      dF(278) = ((y6**2+y7**2)*y2**4+(y5**2+y4**2)*y3**4)*y1
      dF(279) = ((y4*y6+y5*y7)*y2**4+(y4*y6+y5*y7)*y3**4)*y1
      dF(280) = ((y5**2+y4**2)*y2**4+(y6**2+y7**2)*y3**4)*y1
      dF(281) = (y2**2*y3**4+y2**4*y3**2)*y1
      dF(282) = (y2*y3**5+y2**5*y3)*y1
      dF(283) = (y2**6+y3**6)*y1
      dF(284) = (y3*y6**2*y7**2+y2*y4**2*y5**2)*y1**2
      dF(285) = ((y4*y6*y7**2+y5*y6**2*y7)*y2+(y4*y5**2*y6+y4**2*y5*y7)*y3)*y1**2
      dF(286) = ((y6**4+y7**4)*y2+(y4**4+y5**4)*y3)*y1**2
      dF(287) = (y2*y6**2*y7**2+y3*y4**2*y5**2)*y1**2
      dF(288) = ((y4*y6**3+y5*y7**3)*y2+(y5**3*y7+y4**3*y6)*y3)*y1**2
      dF(289) = ((y5**2*y7**2+y4**2*y6**2)*y2+(y5**2*y7**2+y4**2*y6**2)*y3)*y1**2
      dF(290) = (y2*y4*y5*y6*y7+y3*y4*y5*y6*y7)*y1**2
      dF(291) = ((y4*y5**2*y6+y4**2*y5*y7)*y2+(y4*y6*y7**2+y5*y6**2*y7)*y3)*y1**2
      dF(292) = ((y5**2*y6**2+y4**2*y7**2)*y2+(y5**2*y6**2+y4**2*y7**2)*y3)*y1**2
      dF(293) = ((y5**3*y7+y4**3*y6)*y2+(y4*y6**3+y5*y7**3)*y3)*y1**2
      dF(294) = ((y4**4+y5**4)*y2+(y6**4+y7**4)*y3)*y1**2
      dF(295) = ((y6**2+y7**2)*y3*y2**2+(y5**2+y4**2)*y3**2*y2)*y1**2
      dF(296) = ((y4*y6+y5*y7)*y3*y2**2+(y4*y6+y5*y7)*y3**2*y2)*y1**2
      dF(297) = ((y5**2+y4**2)*y3*y2**2+(y6**2+y7**2)*y3**2*y2)*y1**2
      dF(298) = ((y6**2+y7**2)*y2**3+(y5**2+y4**2)*y3**3)*y1**2
      dF(299) = ((y4*y6+y5*y7)*y2**3+(y4*y6+y5*y7)*y3**3)*y1**2
      dF(300) = ((y5**2+y4**2)*y2**3+(y6**2+y7**2)*y3**3)*y1**2
      dF(301) = (y2**2*y3**3+y2**3*y3**2)*y1**2
      dF(302) = (y2*y3**4+y2**4*y3)*y1**2
      dF(303) = (y3**5+y2**5)*y1**2
      dF(304) = (y4**4+y6**4+y7**4+y5**4)*y1**3
      dF(305) = (y6**2*y7**2+y4**2*y5**2)*y1**3
      dF(306) = (y5**2*y6**2+y4**2*y7**2)*y1**3
      dF(307) = (y4**2*y5*y7+(y5**2*y6+y6*y7**2)*y4+y5*y6**2*y7)*y1**3
      dF(308) = y1**3*y4*y5*y6*y7
      dF(309) = (y5**2*y7**2+y4**2*y6**2)*y1**3
      dF(310) = (y5*y7**3+y4*y6**3+y4**3*y6+y5**3*y7)*y1**3
      dF(311) = ((y4*y6+y5*y7)*y2**2+(y4*y6+y5*y7)*y3**2)*y1**3
      dF(312) = (y6**2+y7**2+y5**2+y4**2)*y3*y2*y1**3
      dF(313) = (y4*y6+y5*y7)*y3*y2*y1**3
      dF(314) = (y2**3*y3+y2*y3**3)*y1**3
      dF(315) = ((y6**2+y7**2)*y2**2+(y5**2+y4**2)*y3**2)*y1**3
      dF(316) = ((y5**2+y4**2)*y2**2+(y6**2+y7**2)*y3**2)*y1**3
      dF(317) = y1**3*y2**2*y3**2
      dF(318) = (y2**4+y3**4)*y1**3
      dF(319) = ((y4*y6+y5*y7)*y2+(y4*y6+y5*y7)*y3)*y1**4
      dF(320) = ((y6**2+y7**2)*y2+(y5**2+y4**2)*y3)*y1**4
      dF(321) = ((y5**2+y4**2)*y2+(y6**2+y7**2)*y3)*y1**4
      dF(322) = (y2**2*y3+y2*y3**2)*y1**4
      dF(323) = (y2**3+y3**3)*y1**4
      dF(324) = (y6**2+y7**2+y5**2+y4**2)*y1**5
      dF(325) = (y4*y6+y5*y7)*y1**5
      dF(326) = (y2**2+y3**2)*y1**5
      dF(327) = y1**5*y2*y3
      dF(328) = (y2+y3)*y1**6
      dF(329) = y1**7
      !
      dF(330) = y4**8+y7**8+y5**8+y6**8
      dF(331) = y6**4*y7**4+y4**4*y5**4
      dF(332) = (y6**3*y7**2+y5**2*y6**3)*y4**3+y4**2*y5**3*y7**3+y5**3*y6**2*y7**3
      dF(333) = y4*y5**2*y6**5+y4**2*y5*y7**5+y4**5*y6*y7**2+y5**5*y6**2*y7
      dF(334) = y4**6*y5*y7+(y5**6*y6+y6*y7**6)*y4+y5*y6**6*y7
      dF(335) = y4*y5*y6**3*y7**3+y4**3*y5**3*y6*y7
      dF(336) = y4**3*y5**2*y6*y7**2+(y5*y6**2*y7**3+y5**3*y6**2*y7)*y4**2+y4*y5**2*y6**3*y7**2
      dF(337) = y4*y5**3*y6*y7**3+y4**3*y5*y6**3*y7
      dF(338) = y4*y5**3*y6**3*y7+y4**3*y5*y6*y7**3
      dF(339) = y4**4*y5*y6**2*y7+y4**2*y5*y6**4*y7+(y5**2*y6*y7**4+y5**4*y6*y7**2)*y4
      dF(340) = y4**2*y5**2*y6**2*y7**2
      dF(341) = y4**4*y7**4+y5**4*y6**4
      dF(342) = y4**4*y6**2*y7**2+(y6**4+y7**4)*y5**2*y4**2+y5**4*y6**2*y7**2
      dF(343) = y4**4*y6**4+y5**4*y7**4
      dF(344) = y4**3*y6*y7**4+y4**4*y5*y7**3+y4*y5**4*y6**3+y5**3*y6**4*y7
      dF(345) = y4**4*y5**2*y7**2+(y6**2*y7**4+y5**4*y6**2)*y4**2+y5**2*y6**4*y7**2
      dF(346) = y4**4*y5**2*y6**2+(y6**4*y7**2+y5**4*y7**2)*y4**2+y5**2*y6**2*y7**4
      dF(347) = y5*y6**4*y7**3+y4**3*y5**4*y6+y4*y6**3*y7**4+y4**4*y5**3*y7
      dF(348) = y4**5*y6**3+y5**3*y7**5+y5**5*y7**3+y4**3*y6**5
      dF(349) = y4**5*y5*y6*y7+(y5**5*y6*y7+(y6*y7**5+y6**5*y7)*y5)*y4
      dF(350) = y4*y6**5*y7**2+y4**5*y5**2*y6+y4**2*y5**5*y7+y5*y6**2*y7**5
      dF(351) = y5**6*y6**2+y4**6*y7**2+y4**2*y7**6+y5**2*y6**6
      dF(352) = y4**6*y6**2+y4**2*y6**6+y5**6*y7**2+y5**2*y7**6
      dF(353) = y4**6*y5**2+y6**6*y7**2+y6**2*y7**6+y4**2*y5**6
      dF(354) = y4**7*y6+y4*y6**7+y5**7*y7+y5*y7**7
      dF(355) = (y4**3*y5**2*y6+y4**2*y5**3*y7)*y2**2+(y5*y6**2*y7**3+y4*y6**3*y7**2)*y3**2
      dF(356) = y3**4*y6**2*y7**2+y2**4*y4**2*y5**2
      dF(357) = (y7**6+y6**6+y5**6+y4**6)*y3*y2
      dF(358) = (y4**4*y5**2+y4**2*y5**4+y6**4*y7**2+y6**2*y7**4)*y3*y2
      dF(359) = (y5**5*y7+y4*y6**5+y4**5*y6+y5*y7**5)*y3*y2
      dF(360) = (y4*y6**3*y7**2+y4**2*y5**3*y7+y4**3*y5**2*y6+y5*y6**2*y7**3)*y3*y2
      dF(361) = (y4**4*y5*y7+(y6*y7**4+y5**4*y6)*y4+y5*y6**4*y7)*y3*y2
      dF(362) = (y5**4*y7**2+y4**2*y6**4+y4**4*y6**2+y5**2*y7**4)*y3*y2
      dF(363) = (((y7**2+y6**2)*y5**2+y6**2*y7**2)*y4**2+y5**2*y6**2*y7**2)*y3*y2
      dF(364) = (y5**4*y6**2+y4**4*y7**2+y4**2*y7**4+y5**2*y6**4)*y3*y2
      dF(365) = (y4**3*y6**3+y5**3*y7**3)*y3*y2
      dF(366) = (y4*y5**2*y6**3+y5**3*y6**2*y7+y4**2*y5*y7**3+y4**3*y6*y7**2)*y3*y2
      dF(367) = (y4**3*y5*y6*y7+(y5**3*y6*y7+(y6*y7**3+y6**3*y7)*y5)*y4)*y3*y2
      dF(368) = (y4*y5**2*y6*y7**2+y4**2*y5*y6**2*y7)*y3*y2
      dF(369) = y2*y3**3*y6**2*y7**2+y2**3*y3*y4**2*y5**2
      dF(370) = (y5**2*y7**2+y4**2*y6**2)*y3*y2**3+(y5**2*y7**2+y4**2*y6**2)*y3**3*y2
      dF(371) = (y5**2+y4**2)*y3*y2**5+(y7**2+y6**2)*y3**5*y2
      dF(372) = (y5*y7+y4*y6)*y3*y2**5+(y5*y7+y4*y6)*y3**5*y2
      dF(373) = (y7**2+y6**2)*y3*y2**5+(y5**2+y4**2)*y3**5*y2
      dF(374) = y2*y3**7+y2**7*y3
      dF(375) = (y7**6+y6**6)*y2**2+(y4**6+y5**6)*y3**2
      dF(376) = (y6**4*y7**2+y6**2*y7**4)*y2**2+(y4**2*y5**4+y4**4*y5**2)*y3**2
      dF(377) = (y5*y7**5+y4*y6**5)*y2**2+(y4**5*y6+y5**5*y7)*y3**2
      dF(378) = (y5*y6**2*y7**3+y4*y6**3*y7**2)*y2**2+(y4**3*y5**2*y6+y4**2*y5**3*y7)*y3**2
      dF(379) = (y5**2*y7**4+y4**2*y6**4)*y2**2+(y4**4*y6**2+y5**4*y7**2)*y3**2
      dF(380) = (y5**2*y6**2*y7**2+y4**2*y6**2*y7**2)*y2**2+(y7**2+y6**2)*y5**2*y4**2*y3**2
      dF(381) = (y5**2*y6**4+y4**2*y7**4)*y2**2+(y4**4*y7**2+y5**4*y6**2)*y3**2
      dF(382) = (y5**3*y6**2*y7+y4**3*y6*y7**2)*y2**2+(y4**2*y5*y7**3+y4*y5**2*y6**3)*y3**2
      dF(383) = (y4**5*y6+y5**5*y7)*y2**2+(y5*y7**5+y4*y6**5)*y3**2
      dF(384) = (y4*y6*y7**4+y5*y6**4*y7)*y2**2+(y4*y5**4*y6+y4**4*y5*y7)*y3**2
      dF(385) = (y6*y7**3+y6**3*y7)*y5*y4*y2**2+(y4**3*y5*y6*y7+y4*y5**3*y6*y7)*y3**2
      dF(386) = (y4**2*y5*y7**3+y4*y5**2*y6**3)*y2**2+(y5**3*y6**2*y7+y4**3*y6*y7**2)*y3**2
      dF(387) = (y4*y5**2*y6*y7**2+y4**2*y5*y6**2*y7)*y2**2+(y4*y5**2*y6*y7**2+y4**2*y5*y6**2*y7)*y3**2
      dF(388) = (y7**2+y6**2)*y5**2*y4**2*y2**2+(y5**2*y6**2*y7**2+y4**2*y6**2*y7**2)*y3**2
      dF(389) = (y4**2*y5**4+y4**4*y5**2)*y2**2+(y6**4*y7**2+y6**2*y7**4)*y3**2
      dF(390) = (y4**3*y6**3+y5**3*y7**3)*y2**2+(y4**3*y6**3+y5**3*y7**3)*y3**2
      dF(391) = (y4**3*y5*y6*y7+y4*y5**3*y6*y7)*y2**2+(y6*y7**3+y6**3*y7)*y5*y4*y3**2
      dF(392) = (y4**4*y7**2+y5**4*y6**2)*y2**2+(y5**2*y6**4+y4**2*y7**4)*y3**2
      dF(393) = (y4**4*y6**2+y5**4*y7**2)*y2**2+(y5**2*y7**4+y4**2*y6**4)*y3**2
      dF(394) = (y4*y5**4*y6+y4**4*y5*y7)*y2**2+(y4*y6*y7**4+y5*y6**4*y7)*y3**2
      dF(395) = (y4**6+y5**6)*y2**2+(y7**6+y6**6)*y3**2
      dF(396) = (y6**4+y4**4+y7**4+y5**4)*y3**2*y2**2
      dF(397) = (y4**2*y5*y7+(y6*y7**2+y5**2*y6)*y4+y5*y6**2*y7)*y3**2*y2**2
      dF(398) = (y5**2*y7**2+y4**2*y6**2)*y3**2*y2**2
      dF(399) = (y4*y6**3+y5*y7**3+y5**3*y7+y4**3*y6)*y3**2*y2**2
      dF(400) = y2**2*y3**2*y4*y5*y6*y7
      dF(401) = (y4**2*y7**2+y5**2*y6**2)*y3**2*y2**2
      dF(402) = (y4**2*y5**2+y6**2*y7**2)*y3**2*y2**2
      dF(403) = (y5**2+y4**2)*y3**2*y2**4+(y7**2+y6**2)*y3**4*y2**2
      dF(404) = (y5*y7+y4*y6)*y3**2*y2**4+(y5*y7+y4*y6)*y3**4*y2**2
      dF(405) = (y7**2+y6**2)*y3**2*y2**4+(y5**2+y4**2)*y3**4*y2**2
      dF(406) = y2**2*y3**6+y2**6*y3**2
      dF(407) = (y6**4+y7**4)*y3*y2**3+(y4**4+y5**4)*y3**3*y2
      dF(408) = y2*y3**3*y4**2*y5**2+y2**3*y3*y6**2*y7**2
      dF(409) = (y4*y6*y7**2+y5*y6**2*y7)*y3*y2**3+(y4**2*y5*y7+y4*y5**2*y6)*y3**3*y2
      dF(410) = (y4**2*y7**2+y5**2*y6**2)*y3*y2**3+(y4**2*y7**2+y5**2*y6**2)*y3**3*y2
      dF(411) = (y4**3*y6+y5**3*y7)*y3*y2**3+(y4*y6**3+y5*y7**3)*y3**3*y2
      dF(412) = (y4*y6**3+y5*y7**3)*y3*y2**3+(y4**3*y6+y5**3*y7)*y3**3*y2
      dF(413) = y2*y3**3*y4*y5*y6*y7+y2**3*y3*y4*y5*y6*y7
      dF(414) = (y4**2*y5*y7+y4*y5**2*y6)*y3*y2**3+(y4*y6*y7**2+y5*y6**2*y7)*y3**3*y2
      dF(415) = (y4**4+y5**4)*y3*y2**3+(y6**4+y7**4)*y3**3*y2
      dF(416) = (y5**2+y4**2+y7**2+y6**2)*y3**3*y2**3
      dF(417) = (y5*y7+y4*y6)*y3**3*y2**3
      dF(418) = y3**4*y4**2*y5**2+y2**4*y6**2*y7**2
      dF(419) = (y6**4+y7**4)*y2**4+(y4**4+y5**4)*y3**4
      dF(420) = (y5**2*y7**2+y4**2*y6**2)*y2**4+(y5**2*y7**2+y4**2*y6**2)*y3**4
      dF(421) = (y4*y6*y7**2+y5*y6**2*y7)*y2**4+(y4**2*y5*y7+y4*y5**2*y6)*y3**4
      dF(422) = (y4*y6**3+y5*y7**3)*y2**4+(y4**3*y6+y5**3*y7)*y3**4
      dF(423) = y2**4*y4*y5*y6*y7+y3**4*y4*y5*y6*y7
      dF(424) = (y4**2*y5*y7+y4*y5**2*y6)*y2**4+(y4*y6*y7**2+y5*y6**2*y7)*y3**4
      dF(425) = (y4**2*y7**2+y5**2*y6**2)*y2**4+(y4**2*y7**2+y5**2*y6**2)*y3**4
      dF(426) = (y4**3*y6+y5**3*y7)*y2**4+(y4*y6**3+y5*y7**3)*y3**4
      dF(427) = (y4**4+y5**4)*y2**4+(y6**4+y7**4)*y3**4
      dF(428) = y2**4*y3**4
      dF(429) = y2**5*y3**3+y2**3*y3**5
      dF(430) = (y7**2+y6**2)*y2**6+(y5**2+y4**2)*y3**6
      dF(431) = (y5*y7+y4*y6)*y2**6+(y5*y7+y4*y6)*y3**6
      dF(432) = (y5**2+y4**2)*y2**6+(y7**2+y6**2)*y3**6
      dF(433) = y2**8+y3**8
      dF(434) = ((y7**2+y6**2)*y5**2*y4**2*y2+(y5**2*y6**2*y7**2+y4**2*y6**2*y7**2)*y3)*y1
      dF(435) = ((y4**4*y6**2+y5**4*y7**2)*y2+(y5**2*y7**4+y4**2*y6**4)*y3)*y1
      dF(436) = ((y5**3*y6**2*y7+y4**3*y6*y7**2)*y2+(y4**2*y5*y7**3+y4*y5**2*y6**3)*y3)*y1
      dF(437) = ((y4*y5**2*y6*y7**2+y4**2*y5*y6**2*y7)*y2+(y4*y5**2*y6*y7**2+y4**2*y5*y6**2*y7)*y3)*y1
      dF(438) = ((y5*y6**2*y7**3+y4*y6**3*y7**2)*y2+(y4**3*y5**2*y6+y4**2*y5**3*y7)*y3)*y1
      dF(439) = ((y5**2*y7**4+y4**2*y6**4)*y2+(y4**4*y6**2+y5**4*y7**2)*y3)*y1
      dF(440) = ((y6**4*y7**2+y6**2*y7**4)*y2+(y4**2*y5**4+y4**4*y5**2)*y3)*y1
      dF(441) = ((y7**6+y6**6)*y2+(y4**6+y5**6)*y3)*y1
      dF(442) = ((y6**4+y7**4)*y2**3+(y4**4+y5**4)*y3**3)*y1
      dF(443) = ((y5*y7**5+y4*y6**5)*y2+(y4**5*y6+y5**5*y7)*y3)*y1
      dF(444) = ((y4*y6*y7**4+y5*y6**4*y7)*y2+(y4*y5**4*y6+y4**4*y5*y7)*y3)*y1
      dF(445) = ((y5**2*y6**2*y7**2+y4**2*y6**2*y7**2)*y2+(y7**2+y6**2)*y5**2*y4**2*y3)*y1
      dF(446) = ((y4**3*y6**3+y5**3*y7**3)*y2+(y4**3*y6**3+y5**3*y7**3)*y3)*y1
      dF(447) = ((y6*y7**3+y6**3*y7)*y5*y4*y2+(y4**3*y5*y6*y7+y4*y5**3*y6*y7)*y3)*y1
      dF(448) = ((y5**2*y6**4+y4**2*y7**4)*y2+(y4**4*y7**2+y5**4*y6**2)*y3)*y1
      dF(449) = ((y4**2*y5*y7**3+y4*y5**2*y6**3)*y2+(y5**3*y6**2*y7+y4**3*y6*y7**2)*y3)*y1
      dF(450) = ((y4**3*y5*y6*y7+y4*y5**3*y6*y7)*y2+(y6*y7**3+y6**3*y7)*y5*y4*y3)*y1
      dF(451) = ((y4**3*y5**2*y6+y4**2*y5**3*y7)*y2+(y5*y6**2*y7**3+y4*y6**3*y7**2)*y3)*y1
      dF(452) = ((y4**4*y7**2+y5**4*y6**2)*y2+(y5**2*y6**4+y4**2*y7**4)*y3)*y1
      dF(453) = ((y4*y5**4*y6+y4**4*y5*y7)*y2+(y4*y6*y7**4+y5*y6**4*y7)*y3)*y1
      dF(454) = ((y4**2*y5**4+y4**4*y5**2)*y2+(y6**4*y7**2+y6**2*y7**4)*y3)*y1
      dF(455) = ((y4**5*y6+y5**5*y7)*y2+(y5*y7**5+y4*y6**5)*y3)*y1
      dF(456) = ((y4**6+y5**6)*y2+(y7**6+y6**6)*y3)*y1
      dF(457) = ((y5**2+y4**2)*y3*y2**4+(y7**2+y6**2)*y3**4*y2)*y1
      dF(458) = ((y5*y7+y4*y6)*y3*y2**4+(y5*y7+y4*y6)*y3**4*y2)*y1
      dF(459) = ((y7**2+y6**2)*y3*y2**4+(y5**2+y4**2)*y3**4*y2)*y1
      dF(460) = (y2*y3**6+y2**6*y3)*y1
      dF(461) = ((y6**4+y7**4)*y3*y2**2+(y4**4+y5**4)*y3**2*y2)*y1
      dF(462) = (y2**2*y3*y6**2*y7**2+y2*y3**2*y4**2*y5**2)*y1
      dF(463) = ((y4*y6**3+y5*y7**3)*y3*y2**2+(y4**3*y6+y5**3*y7)*y3**2*y2)*y1
      dF(464) = ((y5**2*y7**2+y4**2*y6**2)*y3*y2**2+(y5**2*y7**2+y4**2*y6**2)*y3**2*y2)*y1
      dF(465) = ((y4**2*y7**2+y5**2*y6**2)*y3*y2**2+(y4**2*y7**2+y5**2*y6**2)*y3**2*y2)*y1
      dF(466) = ((y4**3*y6+y5**3*y7)*y3*y2**2+(y4*y6**3+y5*y7**3)*y3**2*y2)*y1
      dF(467) = ((y4**4+y5**4)*y3*y2**2+(y6**4+y7**4)*y3**2*y2)*y1
      dF(468) = ((y4*y6*y7**2+y5*y6**2*y7)*y3*y2**2+(y4**2*y5*y7+y4*y5**2*y6)*y3**2*y2)*y1
      dF(469) = (y2**2*y3*y4*y5*y6*y7+y2*y3**2*y4*y5*y6*y7)*y1
      dF(470) = ((y4**2*y5*y7+y4*y5**2*y6)*y3*y2**2+(y4*y6*y7**2+y5*y6**2*y7)*y3**2*y2)*y1
      dF(471) = (y2**2*y3*y4**2*y5**2+y2*y3**2*y6**2*y7**2)*y1
      dF(472) = (y2**3*y6**2*y7**2+y3**3*y4**2*y5**2)*y1
      dF(473) = ((y4*y6*y7**2+y5*y6**2*y7)*y2**3+(y4**2*y5*y7+y4*y5**2*y6)*y3**3)*y1
      dF(474) = ((y4*y6**3+y5*y7**3)*y2**3+(y4**3*y6+y5**3*y7)*y3**3)*y1
      dF(475) = (y2**3*y4*y5*y6*y7+y3**3*y4*y5*y6*y7)*y1
      dF(476) = ((y4**2*y7**2+y5**2*y6**2)*y2**3+(y4**2*y7**2+y5**2*y6**2)*y3**3)*y1
      dF(477) = ((y5**2*y7**2+y4**2*y6**2)*y2**3+(y5**2*y7**2+y4**2*y6**2)*y3**3)*y1
      dF(478) = ((y4**2*y5*y7+y4*y5**2*y6)*y2**3+(y4*y6*y7**2+y5*y6**2*y7)*y3**3)*y1
      dF(479) = (y2**3*y4**2*y5**2+y3**3*y6**2*y7**2)*y1
      dF(480) = ((y4**3*y6+y5**3*y7)*y2**3+(y4*y6**3+y5*y7**3)*y3**3)*y1
      dF(481) = ((y4**4+y5**4)*y2**3+(y6**4+y7**4)*y3**3)*y1
      dF(482) = ((y7**2+y6**2)*y3**2*y2**3+(y5**2+y4**2)*y3**3*y2**2)*y1
      dF(483) = ((y5*y7+y4*y6)*y3**2*y2**3+(y5*y7+y4*y6)*y3**3*y2**2)*y1
      dF(484) = ((y5**2+y4**2)*y3**2*y2**3+(y7**2+y6**2)*y3**3*y2**2)*y1
      dF(485) = (y2**4*y3**3+y2**3*y3**4)*y1
      dF(486) = ((y7**2+y6**2)*y2**5+(y5**2+y4**2)*y3**5)*y1
      dF(487) = ((y5*y7+y4*y6)*y2**5+(y5*y7+y4*y6)*y3**5)*y1
      dF(488) = ((y5**2+y4**2)*y2**5+(y7**2+y6**2)*y3**5)*y1
      dF(489) = (y2**5*y3**2+y2**2*y3**5)*y1
      dF(490) = (y2**7+y3**7)*y1
      dF(491) = (y4*y5**2*y6*y7**2+y4**2*y5*y6**2*y7)*y1**2
      dF(492) = (((y7**2+y6**2)*y5**2+y6**2*y7**2)*y4**2+y5**2*y6**2*y7**2)*y1**2
      dF(493) = (y4*y5**2*y6**3+y5**3*y6**2*y7+y4**2*y5*y7**3+y4**3*y6*y7**2)*y1**2
      dF(494) = (y4**3*y6**3+y5**3*y7**3)*y1**2
      dF(495) = (y4**3*y5*y6*y7+(y5**3*y6*y7+(y6*y7**3+y6**3*y7)*y5)*y4)*y1**2
      dF(496) = (y4*y6**3*y7**2+y4**2*y5**3*y7+y4**3*y5**2*y6+y5*y6**2*y7**3)*y1**2
      dF(497) = (y5**4*y6**2+y4**4*y7**2+y4**2*y7**4+y5**2*y6**4)*y1**2
      dF(498) = (y5**4*y7**2+y4**2*y6**4+y4**4*y6**2+y5**2*y7**4)*y1**2
      dF(499) = (y4**4*y5*y7+(y6*y7**4+y5**4*y6)*y4+y5*y6**4*y7)*y1**2
      dF(500) = (y4**4*y5**2+y4**2*y5**4+y6**4*y7**2+y6**2*y7**4)*y1**2
      dF(501) = (y5**5*y7+y4*y6**5+y4**5*y6+y5*y7**5)*y1**2
      dF(502) = (y7**6+y6**6+y5**6+y4**6)*y1**2
      dF(503) = (y3**2*y6**2*y7**2+y2**2*y4**2*y5**2)*y1**2
      dF(504) = ((y4**2*y7**2+y5**2*y6**2)*y2**2+(y4**2*y7**2+y5**2*y6**2)*y3**2)*y1**2
      dF(505) = ((y4*y6**3+y5*y7**3)*y2**2+(y4**3*y6+y5**3*y7)*y3**2)*y1**2
      dF(506) = ((y6**4+y7**4)*y2**2+(y4**4+y5**4)*y3**2)*y1**2
      dF(507) = ((y5**2+y4**2)*y2**4+(y7**2+y6**2)*y3**4)*y1**2
      dF(508) = ((y5*y7+y4*y6)*y2**4+(y5*y7+y4*y6)*y3**4)*y1**2
      dF(509) = ((y7**2+y6**2)*y2**4+(y5**2+y4**2)*y3**4)*y1**2
      dF(510) = (y3**6+y2**6)*y1**2
      dF(511) = (y4**2*y5**2+y6**2*y7**2)*y3*y2*y1**2
      dF(512) = (y6**4+y4**4+y7**4+y5**4)*y3*y2*y1**2
      dF(513) = (y4**2*y5*y7+(y6*y7**2+y5**2*y6)*y4+y5*y6**2*y7)*y3*y2*y1**2
      dF(514) = (y5**2*y7**2+y4**2*y6**2)*y3*y2*y1**2
      dF(515) = (y4**2*y7**2+y5**2*y6**2)*y3*y2*y1**2
      dF(516) = (y4*y6**3+y5*y7**3+y5**3*y7+y4**3*y6)*y3*y2*y1**2
      dF(517) = y1**2*y2*y3*y4*y5*y6*y7
      dF(518) = ((y5**2+y4**2)*y3*y2**3+(y7**2+y6**2)*y3**3*y2)*y1**2
      dF(519) = ((y5*y7+y4*y6)*y3*y2**3+(y5*y7+y4*y6)*y3**3*y2)*y1**2
      dF(520) = ((y7**2+y6**2)*y3*y2**3+(y5**2+y4**2)*y3**3*y2)*y1**2
      dF(521) = (y2**5*y3+y2*y3**5)*y1**2
      dF(522) = (y3**2*y4**2*y5**2+y2**2*y6**2*y7**2)*y1**2
      dF(523) = ((y4*y6*y7**2+y5*y6**2*y7)*y2**2+(y4**2*y5*y7+y4*y5**2*y6)*y3**2)*y1**2
      dF(524) = (y3**2*y4*y5*y6*y7+y2**2*y4*y5*y6*y7)*y1**2
      dF(525) = ((y5**2*y7**2+y4**2*y6**2)*y2**2+(y5**2*y7**2+y4**2*y6**2)*y3**2)*y1**2
      dF(526) = ((y4**2*y5*y7+y4*y5**2*y6)*y2**2+(y4*y6*y7**2+y5*y6**2*y7)*y3**2)*y1**2
      dF(527) = ((y4**3*y6+y5**3*y7)*y2**2+(y4*y6**3+y5*y7**3)*y3**2)*y1**2
      dF(528) = ((y4**4+y5**4)*y2**2+(y6**4+y7**4)*y3**2)*y1**2
      dF(529) = (y5*y7+y4*y6)*y3**2*y2**2*y1**2
      dF(530) = (y5**2+y4**2+y7**2+y6**2)*y3**2*y2**2*y1**2
      dF(531) = y1**2*y2**3*y3**3
      dF(532) = (y2**2*y3**4+y2**4*y3**2)*y1**2
      dF(533) = ((y4**4+y5**4)*y2+(y6**4+y7**4)*y3)*y1**3
      dF(534) = (y2*y4**2*y5**2+y3*y6**2*y7**2)*y1**3
      dF(535) = ((y4**2*y5*y7+y4*y5**2*y6)*y2+(y4*y6*y7**2+y5*y6**2*y7)*y3)*y1**3
      dF(536) = ((y5**2*y7**2+y4**2*y6**2)*y2+(y5**2*y7**2+y4**2*y6**2)*y3)*y1**3
      dF(537) = ((y6**4+y7**4)*y2+(y4**4+y5**4)*y3)*y1**3
      dF(538) = ((y4**2*y7**2+y5**2*y6**2)*y2+(y4**2*y7**2+y5**2*y6**2)*y3)*y1**3
      dF(539) = (y2**5+y3**5)*y1**3
      dF(540) = (y3*y4**2*y5**2+y2*y6**2*y7**2)*y1**3
      dF(541) = ((y4*y6*y7**2+y5*y6**2*y7)*y2+(y4**2*y5*y7+y4*y5**2*y6)*y3)*y1**3
      dF(542) = ((y4*y6**3+y5*y7**3)*y2+(y4**3*y6+y5**3*y7)*y3)*y1**3
      dF(543) = (y3*y4*y5*y6*y7+y2*y4*y5*y6*y7)*y1**3
      dF(544) = ((y4**3*y6+y5**3*y7)*y2+(y4*y6**3+y5*y7**3)*y3)*y1**3
      dF(545) = (y2**4*y3+y2*y3**4)*y1**3
      dF(546) = ((y7**2+y6**2)*y3*y2**2+(y5**2+y4**2)*y3**2*y2)*y1**3
      dF(547) = ((y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2)*y1**3
      dF(548) = ((y5**2+y4**2)*y3*y2**2+(y7**2+y6**2)*y3**2*y2)*y1**3
      dF(549) = (y2**3*y3**2+y2**2*y3**3)*y1**3
      dF(550) = ((y7**2+y6**2)*y2**3+(y5**2+y4**2)*y3**3)*y1**3
      dF(551) = ((y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3)*y1**3
      dF(552) = ((y5**2+y4**2)*y2**3+(y7**2+y6**2)*y3**3)*y1**3
      dF(553) = (y4*y6**3+y5*y7**3+y5**3*y7+y4**3*y6)*y1**4
      dF(554) = y1**4*y4*y5*y6*y7
      dF(555) = (y4**2*y7**2+y5**2*y6**2)*y1**4
      dF(556) = (y5**2*y7**2+y4**2*y6**2)*y1**4
      dF(557) = (y4**2*y5*y7+(y6*y7**2+y5**2*y6)*y4+y5*y6**2*y7)*y1**4
      dF(558) = (y4**2*y5**2+y6**2*y7**2)*y1**4
      dF(559) = (y6**4+y4**4+y7**4+y5**4)*y1**4
      dF(560) = ((y5**2+y4**2)*y2**2+(y7**2+y6**2)*y3**2)*y1**4
      dF(561) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1**4
      dF(562) = ((y7**2+y6**2)*y2**2+(y5**2+y4**2)*y3**2)*y1**4
      dF(563) = (y2**4+y3**4)*y1**4
      dF(564) = (y5**2+y4**2+y7**2+y6**2)*y3*y2*y1**4
      dF(565) = (y5*y7+y4*y6)*y3*y2*y1**4
      dF(566) = (y2**3*y3+y2*y3**3)*y1**4
      dF(567) = y1**4*y2**2*y3**2
      dF(568) = (y2**3+y3**3)*y1**5
      dF(569) = ((y7**2+y6**2)*y2+(y5**2+y4**2)*y3)*y1**5
      dF(570) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**5
      dF(571) = ((y5**2+y4**2)*y2+(y7**2+y6**2)*y3)*y1**5
      dF(572) = (y2*y3**2+y2**2*y3)*y1**5
      dF(573) = (y5*y7+y4*y6)*y1**6
      dF(574) = (y5**2+y4**2+y7**2+y6**2)*y1**6
      dF(575) = y1**6*y2*y3
      dF(576) = (y2**2+y3**2)*y1**6
      dF(577) = (y2+y3)*y1**7
      dF(578) = y1**8
      !
  end subroutine  potC2H2_diff_V


subroutine potC2H2_D8h_diff_V_tmp(n,y1,y2,y3,y4,y5,y6,y7,dF)
    !
    implicit none
    !
    integer(ik),intent(in)  :: n
    real(ark),intent(in)  :: y1,y2,y3,y4,y5,y6,y7
    real(ark),intent(out) :: dF(n)

      !
      dF(:) = 0
      dF(5) = 1.0_ark
      dF(6) = y3+y2
      dF(7) = y1

      dF(8) = y6**2+y5**2+y7**2+y4**2
      dF(9) = y5*y7+y4*y6
      dF(10) = y2*y3
      dF(11) = y3**2+y2**2
      dF(12) = (y3+y2)*y1
      dF(13) = y1**2

      dF(14) = (y6**2+y7**2)*y2+(y5**2+y4**2)*y3
      dF(15) = (y4*y6+y5*y7)*y2+(y4*y6+y5*y7)*y3
      dF(16) = (y5**2+y4**2)*y2+(y6**2+y7**2)*y3
      dF(17) = y2**2*y3+y2*y3**2
      dF(18) = y2**3+y3**3
      dF(19) = (y4*y6+y5*y7)*y1
      dF(20) = (y4**2+y5**2+y6**2+y7**2)*y1
      dF(21) = y1*y2*y3
      dF(22) = (y2**2+y3**2)*y1
      dF(23) = (y2+y3)*y1**2
      dF(24) = y1**3

      dF(25) = -2._ark*y4*y5*y6*y7+y4**2*y7**2+y5**2*y6**2
      dF(26) = 2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2
      dF(27) = y4**3*y6+y4**2*y5*y7+(y5**2*y6+y6*y7**2+y6**3)*y4+y5**3*y7+(y6**2*y7+y7**3)*y5
      dF(28) = y5**4+y7**4+y6**4+y4**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2
      dF(29) = (y4*y6+y5*y7)*y3*y2
      dF(30) = (y4**2+y6**2+y7**2+y5**2)*y3*y2
      dF(31) = (y6**2+y7**2)*y2**2+(y4**2+y5**2)*y3**2
      dF(32) = (y4*y6+y5*y7)*y2**2+(y4*y6+y5*y7)*y3**2
      dF(33) = (y4**2+y5**2)*y2**2+(y6**2+y7**2)*y3**2
      dF(34) = y2**2*y3**2
      dF(35) = y2**3*y3+y2*y3**3
      dF(36) = y2**4+y3**4
      dF(37) = ((y6**2+y7**2)*y2+(y4**2+y5**2)*y3)*y1
      dF(38) = ((y4*y6+y5*y7)*y2+(y4*y6+y5*y7)*y3)*y1
      dF(39) = ((y4**2+y5**2)*y2+(y6**2+y7**2)*y3)*y1
      dF(40) = (y2**2*y3+y2*y3**2)*y1
      dF(41) = (y3**3+y2**3)*y1
      dF(42) = (y4*y6+y5*y7)*y1**2
      dF(43) = (y4**2+y6**2+y7**2+y5**2)*y1**2
      dF(44) = y1**2*y2*y3
      dF(45) = (y2**2+y3**2)*y1**2
      dF(46) = (y2+y3)*y1**3
      dF(47) = y1**4

      dF(48) = (y4**4+2._ark*y4**2*y5**2+y5**4)*y2+(2._ark*y6**2*y7**2+y6**4+y7**4)*y3
      dF(49) = (2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y2+(2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3
      dF(50) = (2._ark*y6**2*y7**2+y6**4+y7**4)*y2+(y4**4+2._ark*y4**2*y5**2+y5**4)*y3
      dF(51) = (y4*y5**2*y6+y4**3*y6+y5**3*y7+y4**2*y5*y7)*y2+((y6*y7**2+y6**3)*y4+(y6**2*y7+y7**3)*y5)*y3
      dF(52) = ((y6*y7**2+y6**3)*y4+(y6**2*y7+y7**3)*y5)*y2+(y4*y5**2*y6+y4**3*y6+y5**3*y7+y4**2*y5*y7)*y3
      dF(53) = (-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y2+(-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y3
      dF(54) = (y6**2+y7**2)*y3*y2**2+(y4**2+y5**2)*y3**2*y2
      dF(55) = (y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2
      dF(56) = (y4**2+y5**2)*y3*y2**2+(y6**2+y7**2)*y3**2*y2
      dF(57) = (y6**2+y7**2)*y2**3+(y4**2+y5**2)*y3**3
      dF(58) = (y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3
      dF(59) = (y4**2+y5**2)*y2**3+(y6**2+y7**2)*y3**3
      dF(60) = y2**3*y3**2+y2**2*y3**3
      dF(61) = y2**4*y3+y2*y3**4
      dF(62) = y2**5+y3**5
      dF(63) = (-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y1
      dF(64) = (2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y1
      dF(65) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y6**2*y7+y7**3)*y5)*y1
      dF(66) = (y4**4+y7**4+2._ark*y4**2*y5**2+y6**4+y5**4+2._ark*y6**2*y7**2)*y1
      dF(67) = (y5*y7+y4*y6)*y3*y2*y1
      dF(68) = (y4**2+y5**2+y7**2+y6**2)*y3*y2*y1
      dF(69) = ((y6**2+y7**2)*y2**2+(y4**2+y5**2)*y3**2)*y1
      dF(70) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1
      dF(71) = ((y4**2+y5**2)*y2**2+(y6**2+y7**2)*y3**2)*y1
      dF(72) = y1*y2**2*y3**2
      dF(73) = (y2**3*y3+y2*y3**3)*y1
      dF(74) = (y2**4+y3**4)*y1
      dF(75) = ((y6**2+y7**2)*y2+(y4**2+y5**2)*y3)*y1**2
      dF(76) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**2
      dF(77) = ((y4**2+y5**2)*y2+(y6**2+y7**2)*y3)*y1**2
      dF(78) = (y2**2*y3+y2*y3**2)*y1**2
      dF(79) = (y2**3+y3**3)*y1**2
      dF(80) = (y5*y7+y4*y6)*y1**3
      dF(81) = (y4**2+y5**2+y7**2+y6**2)*y1**3
      dF(82) = y1**3*y2*y3
      dF(83) = (y2**2+y3**2)*y1**3
      dF(84) = (y2+y3)*y1**4
      dF(85) = y1**5

      dF(86) = y4**3*y6**3+3._ark*y4**2*y5*y6**2*y7+3._ark*y4*y5**2*y6*y7**2+y5**3*y7**3
      dF(87) = y4**4*y6**2+2._ark*y4**3*y5*y6*y7+((y6**2+y7**2)*y5**2+y6**4+y6**2*y7**2)*y4**2+(2._ark*y5**3*y6*y7+(2._ark*y6**3*y7+2._ark*y6*y7**3)*y5)*y4+y5**4*y7**2+(y6**2*y7**2+y7**4)*y5**2
      dF(88) = y4**3*y6*y7**2+(y7**3-2._ark*y6**2*y7)*y5*y4**2+(y6**3-2._ark*y6*y7**2)*y5**2*y4+y5**3*y6**2*y7
      dF(89) = (-y6**2+y7**2)*y4**4-4._ark*y4**3*y5*y6*y7+(-y6**4+y7**4)*y4**2+(-4._ark*y5**3*y6*y7+(-4._ark*y6**3*y7-4._ark*y6*y7**3)*y5)*y4+(-y7**2+y6**2)*y5**4+(y6**4-y7**4)*y5**2
      dF(90) = y4**5*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7+(y6**5+2._ark*y6**3*y7**2+y5**4*y6+y6*y7**4)*y4+y5**5*y7+(y7**5+y6**4*y7+2._ark*y6**2*y7**3)*y5
      dF(91) = y5**6+y6**6+y7**6+3._ark*y4**2*y5**4+y4**6+3._ark*y4**4*y5**2+3._ark*y6**4*y7**2+3._ark*y6**2*y7**4
      dF(92) = (y4**2+y5**2)*y2**4+(y6**2+y7**2)*y3**4
      dF(93) = (y5*y7+y4*y6)*y2**4+(y5*y7+y4*y6)*y3**4
      dF(94) = (y6**2+y7**2)*y2**4+(y4**2+y5**2)*y3**4
      dF(95) = (y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3*y2
      dF(96) = (y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3*y2
      dF(97) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y3*y2
      dF(98) = (y4**4+y7**4+y6**4+y5**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y3*y2
      dF(99) = (y4**2+y5**2)*y3*y2**3+(y6**2+y7**2)*y3**3*y2
      dF(100) = (y5*y7+y4*y6)*y3*y2**3+(y5*y7+y4*y6)*y3**3*y2
      dF(101) = (y6**2+y7**2)*y3*y2**3+(y4**2+y5**2)*y3**3*y2
      dF(102) = y2*y3**5+y2**5*y3
      dF(103) = (y7**4+y6**4+2._ark*y6**2*y7**2)*y2**2+(2._ark*y4**2*y5**2+y5**4+y4**4)*y3**2
      dF(104) = ((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y2**2+(y4**2*y5*y7+y4*y5**2*y6+y4**3*y6+y5**3*y7)*y3**2
      dF(105) = (y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y2**2+(y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3**2
      dF(106) = (y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y2**2+(y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3**2
      dF(107) = (y4**2*y5*y7+y4*y5**2*y6+y4**3*y6+y5**3*y7)*y2**2+((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3**2
      dF(108) = (2._ark*y4**2*y5**2+y5**4+y4**4)*y2**2+(y7**4+y6**4+2._ark*y6**2*y7**2)*y3**2
      dF(109) = (y5*y7+y4*y6)*y3**2*y2**2
      dF(110) = (y4**2+y5**2+y7**2+y6**2)*y3**2*y2**2
      dF(111) = y2**3*y3**3
      dF(112) = y2**4*y3**2+y2**2*y3**4
      dF(113) = y2**6+y3**6
      dF(114) = ((y4**2*y5*y7+y4*y5**2*y6+y4**3*y6+y5**3*y7)*y2+((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3)*y1
      dF(115) = ((y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3)*y1
      dF(116) = ((y7**4+y6**4+2._ark*y6**2*y7**2)*y2+(2._ark*y4**2*y5**2+y5**4+y4**4)*y3)*y1
      dF(117) = (((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y2+(y4**2*y5*y7+y4*y5**2*y6+y4**3*y6+y5**3*y7)*y3)*y1
      dF(118) = ((y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y2+(y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3)*y1
      dF(119) = ((y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y2+(y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3)*y1
      dF(120) = ((2._ark*y4**2*y5**2+y5**4+y4**4)*y2+(y7**4+y6**4+2._ark*y6**2*y7**2)*y3)*y1
      dF(121) = ((y6**2+y7**2)*y3*y2**2+(y4**2+y5**2)*y3**2*y2)*y1
      dF(122) = ((y4**2+y5**2)*y3*y2**2+(y6**2+y7**2)*y3**2*y2)*y1
      dF(123) = ((y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2)*y1
      dF(124) = ((y6**2+y7**2)*y2**3+(y4**2+y5**2)*y3**3)*y1
      dF(125) = ((y4**2+y5**2)*y2**3+(y6**2+y7**2)*y3**3)*y1
      dF(126) = (y2**3*y3**2+y2**2*y3**3)*y1
      dF(127) = (y2**4*y3+y2*y3**4)*y1
      dF(128) = (y2**5+y3**5)*y1
      dF(129) = (y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y1**2
      dF(130) = (y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y1**2
      dF(131) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y1**2
      dF(132) = (y4**4+y7**4+y6**4+y5**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y1**2
      dF(133) = ((y4**2+y5**2)*y2**2+(y6**2+y7**2)*y3**2)*y1**2
      dF(134) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1**2
      dF(135) = ((y6**2+y7**2)*y2**2+(y4**2+y5**2)*y3**2)*y1**2
      dF(136) = (y3**4+y2**4)*y1**2
      dF(137) = (y5*y7+y4*y6)*y3*y2*y1**2
      dF(138) = (y4**2+y5**2+y7**2+y6**2)*y3*y2*y1**2
      dF(139) = y1**2*y2**2*y3**2
      dF(140) = (y2**3*y3+y2*y3**3)*y1**2
      dF(141) = ((y4**2+y5**2)*y2+(y6**2+y7**2)*y3)*y1**3
      dF(142) = ((y6**2+y7**2)*y2+(y4**2+y5**2)*y3)*y1**3
      dF(143) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**3
      dF(144) = (y3**3+y2**3)*y1**3
      dF(145) = (y2*y3**2+y2**2*y3)*y1**3
      dF(146) = (y5*y7+y4*y6)*y1**4
      dF(147) = (y4**2+y5**2+y7**2+y6**2)*y1**4
      dF(148) = (y3**2+y2**2)*y1**4
      dF(149) = y1**4*y2*y3
      dF(150) = (y2+y3)*y1**5
      dF(151) = y1**6

      dF(152) = (y5**2+y4**2)*y2**5+(y7**2+y6**2)*y3**5
      dF(153) = (y5*y7+y4*y6)*y2**5+(y5*y7+y4*y6)*y3**5
      dF(154) = (y7**2+y6**2)*y2**5+(y5**2+y4**2)*y3**5
      dF(155) = y3**7+y2**7
      dF(156) = (y7**6+3._ark*y6**4*y7**2+3._ark*y6**2*y7**4+y6**6)*y2+(3._ark*y4**2*y5**4+y5**6+y4**6+3._ark*y4**4*y5**2)*y3
      dF(157) = ((y6**5+y6*y7**4+2._ark*y6**3*y7**2)*y4+(y7**5+y6**4*y7+2._ark*y6**2*y7**3)*y5)*y2+(y4**5*y6+y5**5*y7+y4**4*y5*y7+2._ark*y4**2*y5**3*y7+y4*y5**4*y6+2._ark*y4**3*y5**2*y6)*y3
      dF(158) = (y4**3*y6**3+y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+3._ark*y4*y5**2*y6*y7**2)*y2+(y4**3*y6**3+y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+3._ark*y4*y5**2*y6*y7**2)*y3
      dF(159) = ((y7**4-y6**4)*y4**2+(-4._ark*y6*y7**3-4._ark*y6**3*y7)*y5*y4+(y6**4-y7**4)*y5**2)*y2+((y7**2-y6**2)*y4**4-4._ark*y4**3*y5*y6*y7-4._ark*y4*y5**3*y6*y7+(-y7**2+y6**2)*y5**4)*y3
      dF(160) = ((y6**2*y7**2+y6**4)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y6**2*y7**2+y7**4)*y5**2)*y2+(y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y3
      dF(161) = (y4**3*y6*y7**2+(y7**3-2._ark*y6**2*y7)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y2+(y4**3*y6*y7**2+(y7**3-2._ark*y6**2*y7)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y3
      dF(162) = (y4**4*y7**2-2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2-2._ark*y4*y5**3*y6*y7+y5**4*y6**2)*y2+((y6**2*y7**2+y7**4)*y4**2+(-2._ark*y6**3*y7-2._ark*y6*y7**3)*y5*y4+(y6**2*y7**2+y6**4)*y5**2)*y3
      dF(163) = ((-y7**2+y6**2)*y4**4+4._ark*y4**3*y5*y6*y7+4._ark*y4*y5**3*y6*y7+(y7**2-y6**2)*y5**4)*y2+((y6**4-y7**4)*y4**2+(4._ark*y6**3*y7+4._ark*y6*y7**3)*y5*y4+(y7**4-y6**4)*y5**2)*y3
      dF(164) = (y4**5*y6+y5**5*y7+y4**4*y5*y7+2._ark*y4**2*y5**3*y7+y4*y5**4*y6+2._ark*y4**3*y5**2*y6)*y2+((y6**5+y6*y7**4+2._ark*y6**3*y7**2)*y4+(y7**5+y6**4*y7+2._ark*y6**2*y7**3)*y5)*y3
      dF(165) = (3._ark*y4**2*y5**4+y5**6+y4**6+3._ark*y4**4*y5**2)*y2+(y7**6+3._ark*y6**4*y7**2+3._ark*y6**2*y7**4+y6**6)*y3
      dF(166) = (y5**2+y4**2)*y3*y2**4+(y7**2+y6**2)*y3**4*y2
      dF(167) = (y5*y7+y4*y6)*y3*y2**4+(y5*y7+y4*y6)*y3**4*y2
      dF(168) = (y7**2+y6**2)*y3*y2**4+(y5**2+y4**2)*y3**4*y2
      dF(169) = y2*y3**6+y2**6*y3
      dF(170) = (y7**4+y6**4+2._ark*y6**2*y7**2)*y3*y2**2+(y5**4+y4**4+2._ark*y4**2*y5**2)*y3**2*y2
      dF(171) = ((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3*y2**2+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3**2*y2
      dF(172) = (y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3*y2**2+(y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3**2*y2
      dF(173) = (y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3*y2**2+(y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3**2*y2
      dF(174) = (y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3*y2**2+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3**2*y2
      dF(175) = (y5**4+y4**4+2._ark*y4**2*y5**2)*y3*y2**2+(y7**4+y6**4+2._ark*y6**2*y7**2)*y3**2*y2
      dF(176) = (y5*y7+y4*y6)*y3**2*y2**3+(y5*y7+y4*y6)*y3**3*y2**2
      dF(177) = (y7**2+y6**2)*y3**2*y2**3+(y5**2+y4**2)*y3**3*y2**2
      dF(178) = y2**2*y3**5+y2**5*y3**2
      dF(179) = (y7**4+y6**4+2._ark*y6**2*y7**2)*y2**3+(y5**4+y4**4+2._ark*y4**2*y5**2)*y3**3
      dF(180) = (y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y2**3+(y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3**3
      dF(181) = (y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y2**3+(y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3**3
      dF(182) = (y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y2**3+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3**3
      dF(183) = (y5**4+y4**4+2._ark*y4**2*y5**2)*y2**3+(y7**4+y6**4+2._ark*y6**2*y7**2)*y3**3
      dF(184) = ((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y2**3+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3**3
      dF(185) = (y5**2+y4**2)*y3**2*y2**3+(y7**2+y6**2)*y3**3*y2**2
      dF(186) = y2**4*y3**3+y2**3*y3**4
      dF(187) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+((y7**2+y6**2)*y5**2+y6**2*y7**2+y6**4)*y4**2+(2._ark*y5**3*y6*y7+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5)*y4+y5**4*y7**2+(y6**2*y7**2+y7**4)*y5**2)*y1
      dF(188) = (y4**3*y6*y7**2+(y7**3-2._ark*y6**2*y7)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y1
      dF(189) = (y4**3*y6**3+y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+3._ark*y4*y5**2*y6*y7**2)*y1
      dF(190) = (y4**4*y7**2-2._ark*y4**3*y5*y6*y7+((y7**2+y6**2)*y5**2+y6**2*y7**2+y7**4)*y4**2+(-2._ark*y5**3*y6*y7+(-2._ark*y6**3*y7-2._ark*y6*y7**3)*y5)*y4+y5**4*y6**2+(y6**2*y7**2+y6**4)*y5**2)*y1
      dF(191) = (y4**5*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7+(y6**5+y5**4*y6+y6*y7**4+2._ark*y6**3*y7**2)*y4+y5**5*y7+(y7**5+y6**4*y7+2._ark*y6**2*y7**3)*y5)*y1
      dF(192) = (y4**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2+3._ark*y4**2*y5**4+3._ark*y4**4*y5**2+y6**6+y7**6+y5**6)*y1
      dF(193) = ((y5**2+y4**2)*y2**4+(y7**2+y6**2)*y3**4)*y1
      dF(194) = ((y5*y7+y4*y6)*y2**4+(y5*y7+y4*y6)*y3**4)*y1
      dF(195) = ((y7**2+y6**2)*y2**4+(y5**2+y4**2)*y3**4)*y1
      dF(196) = (y3**6+y2**6)*y1
      dF(197) = (y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3*y2*y1
      dF(198) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y5**2*y6+y6**3)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y3*y2*y1
      dF(199) = (y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3*y2*y1
      dF(200) = (y4**4+y7**4+y5**4+y6**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y3*y2*y1
      dF(201) = ((y7**2+y6**2)*y3*y2**3+(y5**2+y4**2)*y3**3*y2)*y1
      dF(202) = ((y5*y7+y4*y6)*y3*y2**3+(y5*y7+y4*y6)*y3**3*y2)*y1
      dF(203) = (y2**5*y3+y2*y3**5)*y1
      dF(204) = ((y7**4+y6**4+2._ark*y6**2*y7**2)*y2**2+(y5**4+y4**4+2._ark*y4**2*y5**2)*y3**2)*y1
      dF(205) = (((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y2**2+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3**2)*y1
      dF(206) = ((y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y2**2+(y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3**2)*y1
      dF(207) = ((y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y2**2+(y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3**2)*y1
      dF(208) = ((y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y2**2+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3**2)*y1
      dF(209) = ((y5**4+y4**4+2._ark*y4**2*y5**2)*y2**2+(y7**4+y6**4+2._ark*y6**2*y7**2)*y3**2)*y1
      dF(210) = (y6**2+y7**2+y4**2+y5**2)*y3**2*y2**2*y1
      dF(211) = (y5*y7+y4*y6)*y3**2*y2**2*y1
      dF(212) = (y2**4*y3**2+y2**2*y3**4)*y1
      dF(213) = ((y5**2+y4**2)*y3*y2**3+(y7**2+y6**2)*y3**3*y2)*y1
      dF(214) = y1*y2**3*y3**3
      dF(215) = ((y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y2+(y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y3)*y1**2
      dF(216) = (((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y2+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3)*y1**2
      dF(217) = ((y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y2+(y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y3)*y1**2
      dF(218) = ((y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3)*y1**2
      dF(219) = ((y7**2+y6**2)*y2**3+(y5**2+y4**2)*y3**3)*y1**2
      dF(220) = ((y7**4+y6**4+2._ark*y6**2*y7**2)*y2+(y5**4+y4**4+2._ark*y4**2*y5**2)*y3)*y1**2
      dF(221) = ((y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y2+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3)*y1**2
      dF(222) = ((y5**4+y4**4+2._ark*y4**2*y5**2)*y2+(y7**4+y6**4+2._ark*y6**2*y7**2)*y3)*y1**2
      dF(223) = ((y7**2+y6**2)*y3*y2**2+(y5**2+y4**2)*y3**2*y2)*y1**2
      dF(224) = ((y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2)*y1**2
      dF(225) = ((y5**2+y4**2)*y3*y2**2+(y7**2+y6**2)*y3**2*y2)*y1**2
      dF(226) = ((y5**2+y4**2)*y2**3+(y7**2+y6**2)*y3**3)*y1**2
      dF(227) = (y2**2*y3**3+y2**3*y3**2)*y1**2
      dF(228) = (y2**4*y3+y2*y3**4)*y1**2
      dF(229) = (y2**5+y3**5)*y1**2
      dF(230) = (y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7)*y1**3
      dF(231) = (y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7)*y1**3
      dF(232) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y5**2*y6+y6**3)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y1**3
      dF(233) = (y4**4+y7**4+y5**4+y6**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y1**3
      dF(234) = (y5*y7+y4*y6)*y3*y2*y1**3
      dF(235) = (y6**2+y7**2+y4**2+y5**2)*y3*y2*y1**3
      dF(236) = ((y7**2+y6**2)*y2**2+(y5**2+y4**2)*y3**2)*y1**3
      dF(237) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1**3
      dF(238) = ((y5**2+y4**2)*y2**2+(y7**2+y6**2)*y3**2)*y1**3
      dF(239) = y1**3*y2**2*y3**2
      dF(240) = (y2**3*y3+y2*y3**3)*y1**3
      dF(241) = (y2**4+y3**4)*y1**3
      dF(242) = (y3**3+y2**3)*y1**4
      dF(243) = ((y7**2+y6**2)*y2+(y5**2+y4**2)*y3)*y1**4
      dF(244) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**4
      dF(245) = ((y5**2+y4**2)*y2+(y7**2+y6**2)*y3)*y1**4
      dF(246) = (y2*y3**2+y2**2*y3)*y1**4
      dF(247) = (y6**2+y7**2+y4**2+y5**2)*y1**5
      dF(248) = (y5*y7+y4*y6)*y1**5
      dF(249) = y1**5*y2*y3
      dF(250) = (y2**2+y3**2)*y1**5
      dF(251) = (y3+y2)*y1**6
      dF(252) = y1**7

      dF(253) = y5**4*y6**4+y4**4*y7**4-4._ark*y4*y5**3*y6**3*y7+6._ark*y4**2*y5**2*y6**2*y7**2-4._ark*y4**3*y5*y6*y7**3
      dF(254) = y5**8+y7**8+14._ark*y6**4*y7**4+y6**8+14._ark*y4**4*y5**4+y4**8
      dF(255) = y4**6*y5*y7/3._ark+y4**5*y5**2*y6-4._ark/3._ark*y4**4*y5**3*y7-4._ark/3._ark*y4**3*y5**4*y6+y4**2*y5**5*y7+(y6**5*y7**2-4._ark/3._ark*y6**3*y7**4+y5**6*y6/3._ark+y6*y7**6/3._ark)*y4+(-4._ark/3._ark*y6**4*y7**3+y6**2*y7**5+y6**6*y7/3._ark)*y5
      dF(256) = y4**7*y6+7._ark*y4**4*y5**3*y7+7._ark*y4**3*y5**4*y6+(y6**7+7._ark*y6**3*y7**4)*y4+y5**7*y7+(y7**7+7._ark*y6**4*y7**3)*y5
      dF(257) = (y6**3*y7-y6*y7**3)*y5*y4**3+(y6*y7**3-y6**3*y7)*y5**3*y4
      dF(258) = 3._ark*y4**5*y6*y7**2+y4**4*y5*y7**3+((5._ark*y6**3-9._ark*y6*y7**2)*y5**2+y6*y7**4+5._ark*y6**3*y7**2)*y4**3+((5._ark*y7**3-9._ark*y6**2*y7)*y5**3+(3._ark*y7**5-9._ark*y6**2*y7**3)*y5)*y4**2+(y5**4*y6**3+(-9._ark*y6**3*y7**2+3._ark*y6**5)*y5**2)*y4+3._ark*y5**5*y6**2*y7+(y6**4*y7+5._ark*y6**2*y7**3)*y5**3
      dF(259) = y4**6*y7**2+3._ark*y4**4*y5**2*y6**2-8._ark*y4**3*y5**3*y6*y7+(y7**6+3._ark*y6**4*y7**2+3._ark*y5**4*y7**2)*y4**2-8._ark*y4*y5*y6**3*y7**3+y5**6*y6**2+(y6**6+3._ark*y6**2*y7**4)*y5**2
      dF(260) = y4**6*y6**2+3._ark*y4**5*y5*y6*y7+3._ark*y4**4*y5**2*y6**2+2._ark*y4**3*y5**3*y6*y7+(3._ark*y6**4*y7**2+y6**6+3._ark*y5**4*y7**2)*y4**2+(3._ark*y5**5*y6*y7+(3._ark*y6*y7**5+2._ark*y6**3*y7**3+3._ark*y6**5*y7)*y5)*y4+y5**6*y7**2+(3._ark*y6**2*y7**4+y7**6)*y5**2
      dF(261) = y4**5*y6*y7**2+y4**4*y5*y6**2*y7+((2._ark*y6**3-4._ark*y6*y7**2)*y5**2+2._ark*y6**3*y7**2)*y4**3+((2._ark*y7**3-4._ark*y6**2*y7)*y5**3+(y6**4*y7-4._ark*y6**2*y7**3+y7**5)*y5)*y4**2+(y5**4*y6*y7**2+(y6*y7**4+y6**5-4._ark*y6**3*y7**2)*y5**2)*y4+2._ark*y5**3*y6**2*y7**3+y5**5*y6**2*y7
      dF(262) = y4**4*y6**2*y7**2+(-4._ark*y6**2*y7**2+y7**4+y6**4)*y5**2*y4**2+y5**4*y6**2*y7**2
      dF(263) = y4**2*y5**6-2._ark*y4**4*y5**4+y6**6*y7**2-2._ark*y6**4*y7**4+y6**2*y7**6+y4**6*y5**2
      dF(264) = (-3._ark*y6*y7**2+y6**3)*y4**5+((-5._ark*y6**3+15._ark*y6*y7**2)*y5**2+y6**5-5._ark*y6**3*y7**2)*y4**3+((-5._ark*y7**3+15._ark*y6**2*y7)*y5**3+(-3._ark*y7**5+15._ark*y6**2*y7**3)*y5)*y4**2+(15._ark*y6**3*y7**2-3._ark*y6**5)*y5**2*y4+(-3._ark*y6**2*y7+y7**3)*y5**5+(y7**5-5._ark*y6**2*y7**3)*y5**3
      dF(265) = y4**4*y6**4+4._ark*y4*y5**3*y6**3*y7+6._ark*y4**2*y5**2*y6**2*y7**2+4._ark*y4**3*y5*y6*y7**3+y5**4*y7**4
      dF(266) = -y4**5*y5*y6*y7+(-y6**2+y7**2)*y5**2*y4**4+2._ark*y4**3*y5**3*y6*y7+((-y7**2+y6**2)*y5**4+y6**2*y7**4-y6**4*y7**2)*y4**2+(-y5**5*y6*y7+(-y6*y7**5+2._ark*y6**3*y7**3-y6**5*y7)*y5)*y4+(-y6**2*y7**4+y6**4*y7**2)*y5**2
      dF(267) = (3._ark*y4**4*y5**2+3._ark*y4**2*y5**4+y5**6+y4**6)*y2**2+(y6**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2+y7**6)*y3**2
      dF(268) = ((-y6**4+y7**4)*y4**2+(-4._ark*y6**3*y7-4._ark*y6*y7**3)*y5*y4+(-y7**4+y6**4)*y5**2)*y2**2+((-y6**2+y7**2)*y4**4-4._ark*y4**3*y5*y6*y7-4._ark*y4*y5**3*y6*y7+(-y7**2+y6**2)*y5**4)*y3**2
      dF(269) = ((y6**5+2._ark*y6**3*y7**2+y6*y7**4)*y4+(y6**4*y7+y7**5+2._ark*y6**2*y7**3)*y5)*y2**2+(y5**5*y7+y4**5*y6+y4*y5**4*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7)*y3**2
      dF(270) = (y6**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2+y7**6)*y2**2+(3._ark*y4**4*y5**2+3._ark*y4**2*y5**4+y5**6+y4**6)*y3**2
      dF(271) = (y5**2+y4**2)*y2**6+(y7**2+y6**2)*y3**6
      dF(272) = (y7**2+y6**2)*y2**6+(y5**2+y4**2)*y3**6
      dF(273) = y3**8+y2**8
      dF(274) = (y7**6+y6**6+y5**6+y4**6+3._ark*y4**4*y5**2+3._ark*y4**2*y5**4+3._ark*y6**4*y7**2+3._ark*y6**2*y7**4)*y3*y2
      dF(275) = (y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2)*y3*y2
      dF(276) = ((-y6**2+y7**2)*y4**4-4._ark*y4**3*y5*y6*y7+(-y6**4+y7**4)*y4**2+(-4._ark*y5**3*y6*y7+(-4._ark*y6**3*y7-4._ark*y6*y7**3)*y5)*y4+(-y7**2+y6**2)*y5**4+(-y7**4+y6**4)*y5**2)*y3*y2
      dF(277) = (y4**5*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7+(y6**5+y6*y7**4+y5**4*y6+2._ark*y6**3*y7**2)*y4+y5**5*y7+(y6**4*y7+y7**5+2._ark*y6**2*y7**3)*y5)*y3*y2
      dF(278) = (y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y3*y2
      dF(279) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+((y7**2+y6**2)*y5**2+y6**4+y6**2*y7**2)*y4**2+(2._ark*y5**3*y6*y7+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5)*y4+y5**4*y7**2+(y7**4+y6**2*y7**2)*y5**2)*y3*y2
      dF(280) = (y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3*y2**3+((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3**3*y2
      dF(281) = (y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3*y2**3+(y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3**3*y2
      dF(282) = (y5*y7+y4*y6)*y3*y2**5+(y5*y7+y4*y6)*y3**5*y2
      dF(283) = (y7**2+y6**2)*y3*y2**5+(y5**2+y4**2)*y3**5*y2
      dF(284) = ((y6**4+y6**2*y7**2)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y2**2+(y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y3**2
      dF(285) = (y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2)*y2**2+(y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2)*y3**2
      dF(286) = (y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y2**2+(y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y3**2
      dF(287) = (y4**4*y7**2-2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2-2._ark*y4*y5**3*y6*y7+y5**4*y6**2)*y2**2+((y7**4+y6**2*y7**2)*y4**2+(-2._ark*y6**3*y7-2._ark*y6*y7**3)*y5*y4+(y6**4+y6**2*y7**2)*y5**2)*y3**2
      dF(288) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y2**2+((y6**4+y6**2*y7**2)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y3**2
      dF(289) = (y5**5*y7+y4**5*y6+y4*y5**4*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7)*y2**2+((y6**5+2._ark*y6**3*y7**2+y6*y7**4)*y4+(y6**4*y7+y7**5+2._ark*y6**2*y7**3)*y5)*y3**2
      dF(290) = (y6**4+y5**4+y7**4+y4**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y3**2*y2**2
      dF(291) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y6**2*y7+y7**3)*y5)*y3**2*y2**2
      dF(292) = (y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3**2*y2**2
      dF(293) = (y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3**2*y2**2
      dF(294) = (y5**2+y4**2)*y3**2*y2**4+(y7**2+y6**2)*y3**4*y2**2
      dF(295) = (y6**4+y7**4+2._ark*y6**2*y7**2)*y3*y2**3+(2._ark*y4**2*y5**2+y5**4+y4**4)*y3**3*y2
      dF(296) = ((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3*y2**3+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3**3*y2
      dF(297) = (y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3*y2**3+(y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3**3*y2
      dF(298) = (2._ark*y4**2*y5**2+y5**4+y4**4)*y3*y2**3+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**3*y2
      dF(299) = (y5*y7+y4*y6)*y3**3*y2**3
      dF(300) = (y4**2+y5**2+y7**2+y6**2)*y3**3*y2**3
      dF(301) = (y6**4+y7**4+2._ark*y6**2*y7**2)*y2**4+(2._ark*y4**2*y5**2+y5**4+y4**4)*y3**4
      dF(302) = (2._ark*y4**2*y5**2+y5**4+y4**4)*y2**4+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**4
      dF(303) = ((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y2**4+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3**4
      dF(304) = (y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y2**4+(y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3**4
      dF(305) = (y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y2**4+(y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3**4
      dF(306) = (y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y2**4+((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3**4
      dF(307) = (y7**2+y6**2)*y3**2*y2**4+(y5**2+y4**2)*y3**4*y2**2
      dF(308) = (y5*y7+y4*y6)*y3**2*y2**4+(y5*y7+y4*y6)*y3**4*y2**2
      dF(309) = y2**4*y3**4
      dF(310) = (y5**2+y4**2)*y3*y2**5+(y7**2+y6**2)*y3**5*y2
      dF(311) = y2**5*y3**3+y2**3*y3**5
      dF(312) = (y5*y7+y4*y6)*y2**6+(y5*y7+y4*y6)*y3**6
      dF(313) = y2**6*y3**2+y2**2*y3**6
      dF(314) = y2**7*y3+y2*y3**7
      dF(315) = ((3._ark*y4**4*y5**2+3._ark*y4**2*y5**4+y5**6+y4**6)*y2+(y6**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2+y7**6)*y3)*y1
      dF(316) = ((y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2)*y2+(y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2)*y3)*y1
      dF(317) = ((y6**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2+y7**6)*y2+(3._ark*y4**4*y5**2+3._ark*y4**2*y5**4+y5**6+y4**6)*y3)*y1
      dF(318) = (((y6**4+y6**2*y7**2)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y2+(y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y3)*y1
      dF(319) = (((y6**5+2._ark*y6**3*y7**2+y6*y7**4)*y4+(y6**4*y7+y7**5+2._ark*y6**2*y7**3)*y5)*y2+(y5**5*y7+y4**5*y6+y4*y5**4*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7)*y3)*y1
      dF(320) = ((y6**4+y7**4+2._ark*y6**2*y7**2)*y2**3+(2._ark*y4**2*y5**2+y5**4+y4**4)*y3**3)*y1
      dF(321) = (((y7**4+y6**2*y7**2)*y4**2+(-2._ark*y6**3*y7-2._ark*y6*y7**3)*y5*y4+(y6**4+y6**2*y7**2)*y5**2)*y2+(y4**4*y7**2-2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2-2._ark*y4*y5**3*y6*y7+y5**4*y6**2)*y3)*y1
      dF(322) = ((y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y2+(y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y3)*y1
      dF(323) = (((-y6**2+y7**2)*y4**4-4._ark*y4**3*y5*y6*y7-4._ark*y4*y5**3*y6*y7+(-y7**2+y6**2)*y5**4)*y2+((-y6**4+y7**4)*y4**2+(-4._ark*y6**3*y7-4._ark*y6*y7**3)*y5*y4+(-y7**4+y6**4)*y5**2)*y3)*y1
      dF(324) = ((y5**5*y7+y4**5*y6+y4*y5**4*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7)*y2+((y6**5+2._ark*y6**3*y7**2+y6*y7**4)*y4+(y6**4*y7+y7**5+2._ark*y6**2*y7**3)*y5)*y3)*y1
      dF(325) = ((y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y2+((y6**4+y6**2*y7**2)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y3)*y1
      dF(326) = ((y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3*y2**2+(y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3**2*y2)*y1
      dF(327) = ((y6**4+y7**4+2._ark*y6**2*y7**2)*y3*y2**2+(2._ark*y4**2*y5**2+y5**4+y4**4)*y3**2*y2)*y1
      dF(328) = (((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3*y2**2+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3**2*y2)*y1
      dF(329) = ((y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3*y2**2+(y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3**2*y2)*y1
      dF(330) = ((y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3*y2**2+((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3**2*y2)*y1
      dF(331) = ((2._ark*y4**2*y5**2+y5**4+y4**4)*y3*y2**2+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**2*y2)*y1
      dF(332) = ((y5**2+y4**2)*y3**2*y2**3+(y7**2+y6**2)*y3**3*y2**2)*y1
      dF(333) = ((y5*y7+y4*y6)*y3**2*y2**3+(y5*y7+y4*y6)*y3**3*y2**2)*y1
      dF(334) = ((y7**2+y6**2)*y3**2*y2**3+(y5**2+y4**2)*y3**3*y2**2)*y1
      dF(335) = ((y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y2**3+(y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3**3)*y1
      dF(336) = ((y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y2**3+(y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3**3)*y1
      dF(337) = ((y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y2**3+((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3**3)*y1
      dF(338) = (((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y2**3+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3**3)*y1
      dF(339) = ((2._ark*y4**2*y5**2+y5**4+y4**4)*y2**3+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**3)*y1
      dF(340) = (y2**3*y3**4+y2**4*y3**3)*y1
      dF(341) = ((y7**2+y6**2)*y3*y2**4+(y5**2+y4**2)*y3**4*y2)*y1
      dF(342) = ((y5*y7+y4*y6)*y3*y2**4+(y5*y7+y4*y6)*y3**4*y2)*y1
      dF(343) = ((y5**2+y4**2)*y3*y2**4+(y7**2+y6**2)*y3**4*y2)*y1
      dF(344) = ((y7**2+y6**2)*y2**5+(y5**2+y4**2)*y3**5)*y1
      dF(345) = ((y5*y7+y4*y6)*y2**5+(y5*y7+y4*y6)*y3**5)*y1
      dF(346) = ((y5**2+y4**2)*y2**5+(y7**2+y6**2)*y3**5)*y1
      dF(347) = (y2**5*y3**2+y2**2*y3**5)*y1
      dF(348) = (y2**6*y3+y2*y3**6)*y1
      dF(349) = (y2**7+y3**7)*y1
      dF(350) = (y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y1**2
      dF(351) = (y5**3*y7**3+3._ark*y4**2*y5*y6**2*y7+y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2)*y1**2
      dF(352) = (y4**4*y7**2-2._ark*y4**3*y5*y6*y7+((y7**2+y6**2)*y5**2+y7**4+y6**2*y7**2)*y4**2+(-2._ark*y5**3*y6*y7+(-2._ark*y6**3*y7-2._ark*y6*y7**3)*y5)*y4+y5**4*y6**2+(y6**4+y6**2*y7**2)*y5**2)*y1**2
      dF(353) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+((y7**2+y6**2)*y5**2+y6**4+y6**2*y7**2)*y4**2+(2._ark*y5**3*y6*y7+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5)*y4+y5**4*y7**2+(y7**4+y6**2*y7**2)*y5**2)*y1**2
      dF(354) = (y4**5*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7+(y6**5+y6*y7**4+y5**4*y6+2._ark*y6**3*y7**2)*y4+y5**5*y7+(y6**4*y7+y7**5+2._ark*y6**2*y7**3)*y5)*y1**2
      dF(355) = (y7**6+y6**6+y5**6+y4**6+3._ark*y4**4*y5**2+3._ark*y4**2*y5**4+3._ark*y6**4*y7**2+3._ark*y6**2*y7**4)*y1**2
      dF(356) = ((y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y2**2+(y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3**2)*y1**2
      dF(357) = ((y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y2**2+(y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3**2)*y1**2
      dF(358) = ((y5**2+y4**2)*y2**4+(y7**2+y6**2)*y3**4)*y1**2
      dF(359) = (y6**4+y5**4+y7**4+y4**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y3*y2*y1**2
      dF(360) = (y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3*y2*y1**2
      dF(361) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y6**2*y7+y7**3)*y5)*y3*y2*y1**2
      dF(362) = (y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3*y2*y1**2
      dF(363) = ((y6**4+y7**4+2._ark*y6**2*y7**2)*y2**2+(2._ark*y4**2*y5**2+y5**4+y4**4)*y3**2)*y1**2
      dF(364) = (((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y2**2+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3**2)*y1**2
      dF(365) = ((y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y2**2+((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3**2)*y1**2
      dF(366) = ((2._ark*y4**2*y5**2+y5**4+y4**4)*y2**2+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**2)*y1**2
      dF(367) = (y5*y7+y4*y6)*y3**2*y2**2*y1**2
      dF(368) = (y4**2+y5**2+y7**2+y6**2)*y3**2*y2**2*y1**2
      dF(369) = ((y7**2+y6**2)*y3*y2**3+(y5**2+y4**2)*y3**3*y2)*y1**2
      dF(370) = ((y5*y7+y4*y6)*y3*y2**3+(y5*y7+y4*y6)*y3**3*y2)*y1**2
      dF(371) = ((y5**2+y4**2)*y3*y2**3+(y7**2+y6**2)*y3**3*y2)*y1**2
      dF(372) = y1**2*y2**3*y3**3
      dF(373) = ((y7**2+y6**2)*y2**4+(y5**2+y4**2)*y3**4)*y1**2
      dF(374) = ((y5*y7+y4*y6)*y2**4+(y5*y7+y4*y6)*y3**4)*y1**2
      dF(375) = (y2**2*y3**4+y2**4*y3**2)*y1**2
      dF(376) = (y2**5*y3+y2*y3**5)*y1**2
      dF(377) = (y2**6+y3**6)*y1**2
      dF(378) = ((2._ark*y4**2*y5**2+y5**4+y4**4)*y2+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3)*y1**3
      dF(379) = ((y5**2+y4**2)*y2**3+(y7**2+y6**2)*y3**3)*y1**3
      dF(380) = ((y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3)*y1**3
      dF(381) = ((y7**2+y6**2)*y2**3+(y5**2+y4**2)*y3**3)*y1**3
      dF(382) = (y3**5+y2**5)*y1**3
      dF(383) = ((y6**4+y7**4+2._ark*y6**2*y7**2)*y2+(2._ark*y4**2*y5**2+y5**4+y4**4)*y3)*y1**3
      dF(384) = (((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y2+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3)*y1**3
      dF(385) = ((y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y2+(y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3)*y1**3
      dF(386) = ((y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y2+(y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y3)*y1**3
      dF(387) = ((y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y2+((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3)*y1**3
      dF(388) = ((y5**2+y4**2)*y3*y2**2+(y7**2+y6**2)*y3**2*y2)*y1**3
      dF(389) = ((y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2)*y1**3
      dF(390) = ((y7**2+y6**2)*y3*y2**2+(y5**2+y4**2)*y3**2*y2)*y1**3
      dF(391) = (y2*y3**4+y2**4*y3)*y1**3
      dF(392) = (y2**2*y3**3+y2**3*y3**2)*y1**3
      dF(393) = (y5**2*y7**2+y4**2*y6**2+2._ark*y4*y5*y6*y7)*y1**4
      dF(394) = (y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y1**4
      dF(395) = (y6**4+y5**4+y7**4+y4**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y1**4
      dF(396) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y6**2*y7+y7**3)*y5)*y1**4
      dF(397) = (y5*y7+y4*y6)*y3*y2*y1**4
      dF(398) = (y4**2+y5**2+y7**2+y6**2)*y3*y2*y1**4
      dF(399) = ((y7**2+y6**2)*y2**2+(y5**2+y4**2)*y3**2)*y1**4
      dF(400) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1**4
      dF(401) = ((y5**2+y4**2)*y2**2+(y7**2+y6**2)*y3**2)*y1**4
      dF(402) = y1**4*y2**2*y3**2
      dF(403) = (y2**3*y3+y2*y3**3)*y1**4
      dF(404) = (y2**4+y3**4)*y1**4
      dF(405) = (y3**3+y2**3)*y1**5
      dF(406) = ((y7**2+y6**2)*y2+(y5**2+y4**2)*y3)*y1**5
      dF(407) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**5
      dF(408) = ((y5**2+y4**2)*y2+(y7**2+y6**2)*y3)*y1**5
      dF(409) = (y2*y3**2+y2**2*y3)*y1**5
      dF(410) = (y4**2+y5**2+y7**2+y6**2)*y1**6
      dF(411) = (y5*y7+y4*y6)*y1**6
      dF(412) = y1**6*y2*y3
      dF(413) = (y2**2+y3**2)*y1**6
      dF(414) = (y2+y3)*y1**7
      dF(415) = y1**8
      !
  end subroutine  potC2H2_D8h_diff_V_tmp


subroutine potC2H2_D8h_diff_V(n,y1,y2,y3,y4,y5,y6,y7,dF)
    !
    implicit none
    !
    integer(ik),intent(in)  :: n
    real(ark),intent(in)  :: y1,y2,y3,y4,y5,y6,y7
    real(ark),intent(out) :: dF(n)
      !
      dF(1) = 0._ark      
      dF(2) = 0._ark
      dF(3) = 0._ark
      dF(4) = 0._ark
      dF(5) = 1._ark
      dF(6) = y3+y2
      dF(7) = y1
      dF(8) = y4*y6+y5*y7
      dF(9) = y6**2+y4**2+y5**2+y7**2
      dF(10) = y2*y3
      dF(11) = y2**2+y3**2
      dF(12) = (y2+y3)*y1
      dF(13) = y1**2
      dF(14) = (y6**2+y7**2)*y2+(y4**2+y5**2)*y3
      dF(15) = (y4*y6+y5*y7)*y2+(y4*y6+y5*y7)*y3
      dF(16) = (y4**2+y5**2)*y2+(y6**2+y7**2)*y3
      dF(17) = y2**2*y3+y2*y3**2
      dF(18) = y2**3+y3**3
      dF(19) = (y4*y6+y5*y7)*y1
      dF(20) = (y4**2+y6**2+y5**2+y7**2)*y1
      dF(21) = y1*y2*y3
      dF(22) = (y2**2+y3**2)*y1
      dF(23) = (y2+y3)*y1**2
      dF(24) = y1**3
      dF(25) = y4**2*y7**2+y5**2*y6**2-2._ark*y4*y5*y6*y7
      dF(26) = y4**2*y6**2+y5**2*y7**2+2._ark*y4*y5*y6*y7
      dF(27) = y4**3*y6+y4**2*y5*y7+(y5**2*y6+y6**3+y6*y7**2)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5
      dF(28) = y5**4+y4**4+y7**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2+y6**4
      dF(29) = (y4*y6+y5*y7)*y3*y2
      dF(30) = (y4**2+y7**2+y6**2+y5**2)*y3*y2
      dF(31) = (y7**2+y6**2)*y2**2+(y5**2+y4**2)*y3**2
      dF(32) = (y4*y6+y5*y7)*y2**2+(y4*y6+y5*y7)*y3**2
      dF(33) = (y5**2+y4**2)*y2**2+(y7**2+y6**2)*y3**2
      dF(34) = y2**2*y3**2
      dF(35) = y2**3*y3+y2*y3**3
      dF(36) = y2**4+y3**4
      dF(37) = ((y7**2+y6**2)*y2+(y5**2+y4**2)*y3)*y1
      dF(38) = ((y4*y6+y5*y7)*y2+(y4*y6+y5*y7)*y3)*y1
      dF(39) = ((y5**2+y4**2)*y2+(y7**2+y6**2)*y3)*y1
      dF(40) = (y2**2*y3+y2*y3**2)*y1
      dF(41) = (y2**3+y3**3)*y1
      dF(42) = (y4**2+y7**2+y6**2+y5**2)*y1**2
      dF(43) = (y4*y6+y5*y7)*y1**2
      dF(44) = y1**2*y2*y3
      dF(45) = (y2**2+y3**2)*y1**2
      dF(46) = (y2+y3)*y1**3
      dF(47) = y1**4
      !
      dF(48) = (y5**2+y4**2)*y2**3+(y7**2+y6**2)*y3**3
      dF(49) = (y4*y6+y5*y7)*y2**3+(y4*y6+y5*y7)*y3**3
      dF(50) = (y7**2+y6**2)*y2**3+(y5**2+y4**2)*y3**3
      dF(51) = y2**5+y3**5
      dF(52) = (y6**4+y7**4+2._ark*y6**2*y7**2)*y2+(y4**4+2._ark*y4**2*y5**2+y5**4)*y3
      dF(53) = ((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y2+(y5**3*y7+y4*y5**2*y6+y4**2*y5*y7+y4**3*y6)*y3
      dF(54) = (y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y2+(y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y3
      dF(55) = (y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y2+(y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y3
      dF(56) = (y5**3*y7+y4*y5**2*y6+y4**2*y5*y7+y4**3*y6)*y2+((y6**3+y6*y7**2)*y4+(y6**2*y7+y7**3)*y5)*y3
      dF(57) = (y4**4+2._ark*y4**2*y5**2+y5**4)*y2+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3
      dF(58) = (y5**2+y4**2)*y3*y2**2+(y7**2+y6**2)*y3**2*y2
      dF(59) = (y4*y6+y5*y7)*y3*y2**2+(y4*y6+y5*y7)*y3**2*y2
      dF(60) = (y7**2+y6**2)*y3*y2**2+(y5**2+y4**2)*y3**2*y2
      dF(61) = y2*y3**4+y2**4*y3
      dF(62) = y2**3*y3**2+y2**2*y3**3
      dF(63) = (y6**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2+y7**4+y4**4+y5**4)*y1
      dF(64) = (y4**3*y6+y4**2*y5*y7+(y5**2*y6+y6**3+y6*y7**2)*y4+y5**3*y7+(y6**2*y7+y7**3)*y5)*y1
      dF(65) = (y4**2*y7**2-2._ark*y4*y5*y6*y7+y5**2*y6**2)*y1
      dF(66) = (y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y1
      dF(67) = (y4*y6+y5*y7)*y3*y2*y1
      dF(68) = (y7**2+y5**2+y6**2+y4**2)*y3*y2*y1
      dF(69) = (y2**3*y3+y2*y3**3)*y1
      dF(70) = ((y7**2+y6**2)*y2**2+(y5**2+y4**2)*y3**2)*y1
      dF(71) = ((y4*y6+y5*y7)*y2**2+(y4*y6+y5*y7)*y3**2)*y1
      dF(72) = ((y5**2+y4**2)*y2**2+(y7**2+y6**2)*y3**2)*y1
      dF(73) = y1*y2**2*y3**2
      dF(74) = (y2**4+y3**4)*y1
      dF(75) = ((y5**2+y4**2)*y2+(y7**2+y6**2)*y3)*y1**2
      dF(76) = ((y4*y6+y5*y7)*y2+(y4*y6+y5*y7)*y3)*y1**2
      dF(77) = ((y7**2+y6**2)*y2+(y5**2+y4**2)*y3)*y1**2
      dF(78) = (y2*y3**2+y2**2*y3)*y1**2
      dF(79) = (y3**3+y2**3)*y1**2
      dF(80) = (y7**2+y5**2+y6**2+y4**2)*y1**3
      dF(81) = (y4*y6+y5*y7)*y1**3
      dF(82) = y1**3*y2*y3
      dF(83) = (y3**2+y2**2)*y1**3
      dF(84) = (y3+y2)*y1**4
      dF(85) = y1**5
      !
      dF(86) = y5**3*y7**3+y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2+3._ark*y4**2*y5*y6**2*y7
      dF(87) = y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(y6**3-2._ark*y6*y7**2)*y5**2*y4+y5**3*y6**2*y7
      dF(88) = y4**4*y6**2+2._ark*y4**3*y5*y6*y7+((y6**2+y7**2)*y5**2+y6**2*y7**2+y6**4)*y4**2+(2._ark*y5**3*y6*y7+(2._ark*y6**3*y7+2._ark*y6*y7**3)*y5)*y4+y5**4*y7**2+(y6**2*y7**2+y7**4)*y5**2
      dF(89) = y4**4*y7**2-2._ark*y4**3*y5*y6*y7+((y6**2+y7**2)*y5**2+y6**2*y7**2+y7**4)*y4**2+(-2._ark*y5**3*y6*y7+(-2._ark*y6*y7**3-2._ark*y6**3*y7)*y5)*y4+y5**4*y6**2+(y6**2*y7**2+y6**4)*y5**2
      dF(90) = y4**5*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7+(y5**4*y6+2._ark*y6**3*y7**2+y6*y7**4+y6**5)*y4+y5**5*y7+(y6**4*y7+2._ark*y6**2*y7**3+y7**5)*y5
      dF(91) = y4**6+y5**6+y7**6+y6**6+3._ark*y6**4*y7**2+3._ark*y4**4*y5**2+3._ark*y6**2*y7**4+3._ark*y4**2*y5**4
      dF(92) = (y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y2**2+(y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y3**2
      dF(93) = ((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y2**2+(y5**3*y7+y4*y5**2*y6+y4**3*y6+y4**2*y5*y7)*y3**2
      dF(94) = (y7**4+2._ark*y6**2*y7**2+y6**4)*y2**2+(y5**4+y4**4+2._ark*y4**2*y5**2)*y3**2
      dF(95) = (y5**3*y7+y4*y5**2*y6+y4**3*y6+y4**2*y5*y7)*y2**2+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3**2
      dF(96) = (y5**2*y6**2-2._ark*y4*y5*y6*y7+y4**2*y7**2)*y2**2+(y5**2*y6**2-2._ark*y4*y5*y6*y7+y4**2*y7**2)*y3**2
      dF(97) = (y4**2+y5**2)*y2**4+(y6**2+y7**2)*y3**4
      dF(98) = (y5*y7+y4*y6)*y2**4+(y5*y7+y4*y6)*y3**4
      dF(99) = (y6**2+y7**2)*y2**4+(y4**2+y5**2)*y3**4
      dF(100) = y3**6+y2**6
      dF(101) = (y5**2*y6**2-2._ark*y4*y5*y6*y7+y4**2*y7**2)*y3*y2
      dF(102) = (y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y3*y2
      dF(103) = (y4**3*y6+y4**2*y5*y7+(y5**2*y6+y6*y7**2+y6**3)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y3*y2
      dF(104) = (y4**4+y7**4+y6**4+y5**4+2._ark*y6**2*y7**2+2._ark*y4**2*y5**2)*y3*y2
      dF(105) = (y5*y7+y4*y6)*y3*y2**3+(y5*y7+y4*y6)*y3**3*y2
      dF(106) = (y5**4+y4**4+2._ark*y4**2*y5**2)*y2**2+(y7**4+2._ark*y6**2*y7**2+y6**4)*y3**2
      dF(107) = (y5*y7+y4*y6)*y3**2*y2**2
      dF(108) = (y4**2+y7**2+y5**2+y6**2)*y3**2*y2**2
      dF(109) = y2**2*y3**4+y2**4*y3**2
      dF(110) = (y6**2+y7**2)*y3*y2**3+(y4**2+y5**2)*y3**3*y2
      dF(111) = (y4**2+y5**2)*y3*y2**3+(y6**2+y7**2)*y3**3*y2
      dF(112) = y2**3*y3**3
      dF(113) = y2**5*y3+y2*y3**5
      dF(114) = ((y5**4+y4**4+2._ark*y4**2*y5**2)*y2+(y7**4+2._ark*y6**2*y7**2+y6**4)*y3)*y1
      dF(115) = ((y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y2+(y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y3)*y1
      dF(116) = (((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y2+(y5**3*y7+y4*y5**2*y6+y4**3*y6+y4**2*y5*y7)*y3)*y1
      dF(117) = ((y7**4+2._ark*y6**2*y7**2+y6**4)*y2+(y5**4+y4**4+2._ark*y4**2*y5**2)*y3)*y1
      dF(118) = ((y5**2*y6**2-2._ark*y4*y5*y6*y7+y4**2*y7**2)*y2+(y5**2*y6**2-2._ark*y4*y5*y6*y7+y4**2*y7**2)*y3)*y1
      dF(119) = ((y5**3*y7+y4*y5**2*y6+y4**3*y6+y4**2*y5*y7)*y2+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3)*y1
      dF(120) = ((y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2)*y1
      dF(121) = ((y6**2+y7**2)*y3*y2**2+(y4**2+y5**2)*y3**2*y2)*y1
      dF(122) = ((y4**2+y5**2)*y3*y2**2+(y6**2+y7**2)*y3**2*y2)*y1
      dF(123) = (y2**2*y3**3+y2**3*y3**2)*y1
      dF(124) = ((y6**2+y7**2)*y2**3+(y4**2+y5**2)*y3**3)*y1
      dF(125) = ((y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3)*y1
      dF(126) = ((y4**2+y5**2)*y2**3+(y6**2+y7**2)*y3**3)*y1
      dF(127) = (y2**4*y3+y2*y3**4)*y1
      dF(128) = (y2**5+y3**5)*y1
      dF(129) = (y4**3*y6+y4**2*y5*y7+(y5**2*y6+y6*y7**2+y6**3)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y1**2
      dF(130) = (y5**2*y6**2-2._ark*y4*y5*y6*y7+y4**2*y7**2)*y1**2
      dF(131) = (y4**2*y6**2+2._ark*y4*y5*y6*y7+y5**2*y7**2)*y1**2
      dF(132) = (y4**4+y7**4+y6**4+y5**4+2._ark*y6**2*y7**2+2._ark*y4**2*y5**2)*y1**2
      dF(133) = (y4**2+y7**2+y5**2+y6**2)*y3*y2*y1**2
      dF(134) = (y5*y7+y4*y6)*y3*y2*y1**2
      dF(135) = ((y6**2+y7**2)*y2**2+(y4**2+y5**2)*y3**2)*y1**2
      dF(136) = ((y4**2+y5**2)*y2**2+(y6**2+y7**2)*y3**2)*y1**2
      dF(137) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1**2
      dF(138) = y1**2*y2**2*y3**2
      dF(139) = (y2**3*y3+y2*y3**3)*y1**2
      dF(140) = (y2**4+y3**4)*y1**2
      dF(141) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**3
      dF(142) = ((y6**2+y7**2)*y2+(y4**2+y5**2)*y3)*y1**3
      dF(143) = ((y4**2+y5**2)*y2+(y6**2+y7**2)*y3)*y1**3
      dF(144) = (y2**2*y3+y2*y3**2)*y1**3
      dF(145) = (y2**3+y3**3)*y1**3
      dF(146) = (y4**2+y7**2+y5**2+y6**2)*y1**4
      dF(147) = (y5*y7+y4*y6)*y1**4
      dF(148) = (y3**2+y2**2)*y1**4
      dF(149) = y1**4*y2*y3
      dF(150) = (y3+y2)*y1**5
      dF(151) = y1**6

      dF(152) = (y4*y5**4*y6+y4**4*y5*y7+y5**5*y7+y4**5*y6+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7)*y2+((y6*y7**4+2._ark*y6**3*y7**2+y6**5)*y4+(y7**5+2._ark*y6**2*y7**3+y6**4*y7)*y5)*y3
      dF(153) = ((y7**2-y6**2)*y4**4-4._ark*y4**3*y5*y6*y7-4._ark*y4*y5**3*y6*y7+(y6**2-y7**2)*y5**4)*y2+((y7**4-y6**4)*y4**2+(-4._ark*y6*y7**3-4._ark*y6**3*y7)*y5*y4+(y6**4-y7**4)*y5**2)*y3
      dF(154) = ((y6**2*y7**2+y6**4)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y2+(y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y6**2+y7**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y3
      dF(155) = (y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y2**3+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3**3
      dF(156) = (-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y2**3+(-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y3**3
      dF(157) = (2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y2**3+(2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y3**3
      dF(158) = (y5**2+y4**2)*y2**5+(y6**2+y7**2)*y3**5
      dF(159) = (y4*y6+y5*y7)*y2**5+(y4*y6+y5*y7)*y3**5
      dF(160) = (y6**2+y7**2)*y2**5+(y5**2+y4**2)*y3**5
      dF(161) = y3**7+y2**7
      dF(162) = (y6**6+y7**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2)*y2+(3._ark*y4**4*y5**2+3._ark*y4**2*y5**4+y5**6+y4**6)*y3
      dF(163) = ((y6*y7**4+2._ark*y6**3*y7**2+y6**5)*y4+(y7**5+2._ark*y6**2*y7**3+y6**4*y7)*y5)*y2+(y4*y5**4*y6+y4**4*y5*y7+y5**5*y7+y4**5*y6+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7)*y3
      dF(164) = ((y7**4+y6**2*y7**2)*y4**2+(-2._ark*y6**3*y7-2._ark*y6*y7**3)*y5*y4+(y6**2*y7**2+y6**4)*y5**2)*y2+(y4**4*y7**2-2._ark*y4**3*y5*y6*y7+(y6**2+y7**2)*y5**2*y4**2-2._ark*y4*y5**3*y6*y7+y5**4*y6**2)*y3
      dF(165) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y6**2+y7**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y2+((y6**2*y7**2+y6**4)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y3
      dF(166) = (y4**3*y6*y7**2+(y7**3-2._ark*y6**2*y7)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y2+(y4**3*y6*y7**2+(y7**3-2._ark*y6**2*y7)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y3
      dF(167) = (y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2+3._ark*y4**2*y5*y6**2*y7+y5**3*y7**3)*y2+(y4**3*y6**3+3._ark*y4*y5**2*y6*y7**2+3._ark*y4**2*y5*y6**2*y7+y5**3*y7**3)*y3
      dF(168) = (3._ark*y4**4*y5**2+3._ark*y4**2*y5**4+y5**6+y4**6)*y2+(y6**6+y7**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2)*y3
      dF(169) = (y4**4+2._ark*y4**2*y5**2+y5**4)*y3*y2**2+(y7**4+2._ark*y6**2*y7**2+y6**4)*y3**2*y2
      dF(170) = (-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y3*y2**2+(-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y3**2*y2
      dF(171) = (y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3*y2**2+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3**2*y2
      dF(172) = ((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3*y2**2+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3**2*y2
      dF(173) = (y5**2+y4**2)*y3*y2**4+(y6**2+y7**2)*y3**4*y2
      dF(174) = (y4*y6+y5*y7)*y3*y2**4+(y4*y6+y5*y7)*y3**4*y2
      dF(175) = (y6**2+y7**2)*y3*y2**4+(y5**2+y4**2)*y3**4*y2
      dF(176) = (y7**4+2._ark*y6**2*y7**2+y6**4)*y3*y2**2+(y4**4+2._ark*y4**2*y5**2+y5**4)*y3**2*y2
      dF(177) = (2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y3*y2**2+(2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y3**2*y2
      dF(178) = (y5**2+y4**2)*y3**2*y2**3+(y6**2+y7**2)*y3**3*y2**2
      dF(179) = y2**5*y3**2+y2**2*y3**5
      dF(180) = (y7**4+2._ark*y6**2*y7**2+y6**4)*y2**3+(y4**4+2._ark*y4**2*y5**2+y5**4)*y3**3
      dF(181) = (y4**4+2._ark*y4**2*y5**2+y5**4)*y2**3+(y7**4+2._ark*y6**2*y7**2+y6**4)*y3**3
      dF(182) = ((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y2**3+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3**3
      dF(183) = (y6**2+y7**2)*y3**2*y2**3+(y5**2+y4**2)*y3**3*y2**2
      dF(184) = (y4*y6+y5*y7)*y3**2*y2**3+(y4*y6+y5*y7)*y3**3*y2**2
      dF(185) = y2**3*y3**4+y2**4*y3**3
      dF(186) = y2**6*y3+y2*y3**6
      dF(187) = (y4**5*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7+(y6**5+y5**4*y6+y6*y7**4+2._ark*y6**3*y7**2)*y4+y5**5*y7+(y7**5+2._ark*y6**2*y7**3+y6**4*y7)*y5)*y1
      dF(188) = (y4*y5**2*y6*y7**2+y4**2*y5*y6**2*y7+y5**3*y7**3/3._ark+y4**3*y6**3/3._ark)*y1
      dF(189) = ((y7**2-y6**2)*y4**4-4._ark*y4**3*y5*y6*y7+(y7**4-y6**4)*y4**2+(-4._ark*y5**3*y6*y7+(-4._ark*y6*y7**3-4._ark*y6**3*y7)*y5)*y4+(y6**2-y7**2)*y5**4+(y6**4-y7**4)*y5**2)*y1
      dF(190) = ((y6*y7**2+2._ark/3._ark*y6**3)*y4**3+y4**2*y5*y7**3+y4*y5**2*y6**3+(2._ark/3._ark*y7**3+y6**2*y7)*y5**3)*y1
      dF(191) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+((y6**2+y7**2)*y5**2+y6**2*y7**2+y6**4)*y4**2+(2._ark*y5**3*y6*y7+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5)*y4+y5**4*y7**2+(y7**4+y6**2*y7**2)*y5**2)*y1
      dF(192) = (y4**6+y7**6+y6**6+y5**6+3._ark*y4**4*y5**2+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2+3._ark*y4**2*y5**4)*y1
      dF(193) = ((-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y2**2+(-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y3**2)*y1
      dF(194) = ((y5**2+y4**2)*y2**4+(y6**2+y7**2)*y3**4)*y1
      dF(195) = ((y4*y6+y5*y7)*y2**4+(y4*y6+y5*y7)*y3**4)*y1
      dF(196) = ((y6**2+y7**2)*y2**4+(y5**2+y4**2)*y3**4)*y1
      dF(197) = (y7**4+y5**4+y4**4+y6**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y3*y2*y1
      dF(198) = (y4**3*y6+y4**2*y5*y7+(y5**2*y6+y6**3+y6*y7**2)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y3*y2*y1
      dF(199) = (2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y3*y2*y1
      dF(200) = (-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y3*y2*y1
      dF(201) = ((y5**2+y4**2)*y3*y2**3+(y6**2+y7**2)*y3**3*y2)*y1
      dF(202) = ((y4*y6+y5*y7)*y3*y2**3+(y4*y6+y5*y7)*y3**3*y2)*y1
      dF(203) = ((y6**2+y7**2)*y3*y2**3+(y5**2+y4**2)*y3**3*y2)*y1
      dF(204) = (y2*y3**5+y2**5*y3)*y1
      dF(205) = ((y7**4+2._ark*y6**2*y7**2+y6**4)*y2**2+(y4**4+2._ark*y4**2*y5**2+y5**4)*y3**2)*y1
      dF(206) = (((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y2**2+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3**2)*y1
      dF(207) = ((2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y2**2+(2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y3**2)*y1
      dF(208) = ((y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y2**2+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3**2)*y1
      dF(209) = ((y4**2*y5**2+y5**4/2._ark+y4**4/2._ark)*y2**2+(y6**2*y7**2+y7**4/2._ark+y6**4/2._ark)*y3**2)*y1
      dF(210) = (y4*y6+y5*y7)*y3**2*y2**2*y1
      dF(211) = (y7**2+y5**2+y6**2+y4**2)*y3**2*y2**2*y1
      dF(212) = (y2**2*y3**4+y2**4*y3**2)*y1
      dF(213) = y1*y2**3*y3**3
      dF(214) = (y3**6+y2**6)*y1
      dF(215) = ((-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y2+(-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y3)*y1**2
      dF(216) = (((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y2+(y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y3)*y1**2
      dF(217) = ((y7**4+2._ark*y6**2*y7**2+y6**4)*y2+(y4**4+2._ark*y4**2*y5**2+y5**4)*y3)*y1**2
      dF(218) = ((y6**2+y7**2)*y2**3+(y5**2+y4**2)*y3**3)*y1**2
      dF(219) = ((y4*y5**2*y6+y4**2*y5*y7+y4**3*y6+y5**3*y7)*y2+((y6*y7**2+y6**3)*y4+(y7**3+y6**2*y7)*y5)*y3)*y1**2
      dF(220) = ((2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y2+(2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y3)*y1**2
      dF(221) = ((y4**4+2._ark*y4**2*y5**2+y5**4)*y2+(y7**4+2._ark*y6**2*y7**2+y6**4)*y3)*y1**2
      dF(222) = ((y4*y6+y5*y7)*y3*y2**2+(y4*y6+y5*y7)*y3**2*y2)*y1**2
      dF(223) = ((y6**2+y7**2)*y3*y2**2+(y5**2+y4**2)*y3**2*y2)*y1**2
      dF(224) = (y2*y3**4+y2**4*y3)*y1**2
      dF(225) = ((y5**2+y4**2)*y3*y2**2+(y6**2+y7**2)*y3**2*y2)*y1**2
      dF(226) = ((y4*y6+y5*y7)*y2**3+(y4*y6+y5*y7)*y3**3)*y1**2
      dF(227) = ((y5**2+y4**2)*y2**3+(y6**2+y7**2)*y3**3)*y1**2
      dF(228) = (y2**2*y3**3+y2**3*y3**2)*y1**2
      dF(229) = (y2**5+y3**5)*y1**2
      dF(230) = (-2._ark*y4*y5*y6*y7+y5**2*y6**2+y4**2*y7**2)*y1**3
      dF(231) = (y4**3*y6+y4**2*y5*y7+(y5**2*y6+y6**3+y6*y7**2)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y1**3
      dF(232) = (y7**4+y5**4+y4**4+y6**4+2._ark*y4**2*y5**2+2._ark*y6**2*y7**2)*y1**3
      dF(233) = (2._ark*y4*y5*y6*y7+y4**2*y6**2+y5**2*y7**2)*y1**3
      dF(234) = ((y5**2+y4**2)*y2**2+(y6**2+y7**2)*y3**2)*y1**3
      dF(235) = ((y6**2+y7**2)*y2**2+(y5**2+y4**2)*y3**2)*y1**3
      dF(236) = ((y4*y6+y5*y7)*y2**2+(y4*y6+y5*y7)*y3**2)*y1**3
      dF(237) = (y7**2+y5**2+y6**2+y4**2)*y3*y2*y1**3
      dF(238) = (y4*y6+y5*y7)*y3*y2*y1**3
      dF(239) = y1**3*y2**2*y3**2
      dF(240) = (y2**3*y3+y2*y3**3)*y1**3
      dF(241) = (y2**4+y3**4)*y1**3
      dF(242) = ((y5**2+y4**2)*y2+(y6**2+y7**2)*y3)*y1**4
      dF(243) = ((y4*y6+y5*y7)*y2+(y4*y6+y5*y7)*y3)*y1**4
      dF(244) = ((y6**2+y7**2)*y2+(y5**2+y4**2)*y3)*y1**4
      dF(245) = (y3**3+y2**3)*y1**4
      dF(246) = (y2*y3**2+y2**2*y3)*y1**4
      dF(247) = (y7**2+y5**2+y6**2+y4**2)*y1**5
      dF(248) = (y4*y6+y5*y7)*y1**5
      dF(249) = (y3**2+y2**2)*y1**5
      dF(250) = y1**5*y2*y3
      dF(251) = (y2+y3)*y1**6
      dF(252) = y1**7
      !
      dF(253) = y7**8+y5**8+y6**8+4._ark*y4**2*y5**6+4._ark*y4**6*y5**2+4._ark*y6**6*y7**2+6._ark*y6**4*y7**4+6._ark*y4**4*y5**4+y4**8+4._ark*y6**2*y7**6
      dF(254) = y4**7*y6+y4**6*y5*y7+3._ark*y4**5*y5**2*y6+3._ark*y4**4*y5**3*y7+3._ark*y4**3*y5**4*y6+3._ark*y4**2*y5**5*y7+(y6*y7**6+y5**6*y6+y6**7+3._ark*y6**5*y7**2+3._ark*y6**3*y7**4)*y4+y5**7*y7+(3._ark*y6**4*y7**3+y7**7+y6**6*y7+3._ark*y6**2*y7**5)*y5
      dF(255) = (-y6**2/3._ark+2._ark/3._ark*y7**2)*y4**6-2._ark*y4**5*y5*y6*y7+y4**4*y5**2*y7**2-4._ark*y4**3*y5**3*y6*y7+(2._ark/3._ark*y7**6+y5**4*y6**2-y6**6/3._ark+y6**2*y7**4)*y4**2+(-2._ark*y5**5*y6*y7+(-2._ark*y6*y7**5-2._ark*y6**5*y7-4._ark*y6**3*y7**3)*y5)*y4+(-y7**2/3._ark+2._ark/3._ark*y6**2)*y5**6+(-y7**6/3._ark+y6**4*y7**2+2._ark/3._ark*y6**6)*y5**2
      dF(256) = 4._ark*y4**3*y5*y6**3*y7+6._ark*y4**2*y5**2*y6**2*y7**2+y4**4*y6**4+4._ark*y4*y5**3*y6*y7**3+y5**4*y7**4
      dF(257) = (2._ark*y6**2*y7**2+y7**4)*y4**4-4._ark*y4**3*y5*y6**3*y7+(2._ark*y6**4-2._ark*y6**2*y7**2+2._ark*y7**4)*y5**2*y4**2-4._ark*y4*y5**3*y6*y7**3+(2._ark*y6**2*y7**2+y6**4)*y5**4
      dF(258) = (-y7**2/3._ark+2._ark/3._ark*y6**2)*y4**6+2._ark*y4**5*y5*y6*y7+y4**4*y5**2*y6**2+4._ark*y4**3*y5**3*y6*y7+(y6**4*y7**2+y5**4*y7**2+2._ark/3._ark*y6**6-y7**6/3._ark)*y4**2+(2._ark*y5**5*y6*y7+(4._ark*y6**3*y7**3+2._ark*y6*y7**5+2._ark*y6**5*y7)*y5)*y4+(-y6**2/3._ark+2._ark/3._ark*y7**2)*y5**6+(y6**2*y7**4-y6**6/3._ark+2._ark/3._ark*y7**6)*y5**2
      dF(259) = y4**4*y6**2*y7**2/2._ark+(-y6**3*y7+y6*y7**3)*y5*y4**3+(y6**4/2._ark-2._ark*y6**2*y7**2+y7**4/2._ark)*y5**2*y4**2+(-y6*y7**3+y6**3*y7)*y5**3*y4+y5**4*y6**2*y7**2/2._ark
      dF(260) = (2._ark/3._ark*y6**3+y6*y7**2)*y4**5+y4**4*y5*y7**3+((5._ark/3._ark*y6**3+y6*y7**2)*y5**2+5._ark/3._ark*y6**3*y7**2+y6*y7**4+2._ark/3._ark*y6**5)*y4**3+((y6**2*y7+5._ark/3._ark*y7**3)*y5**3+(y7**5+y6**2*y7**3)*y5)*y4**2+(y5**4*y6**3+(y6**5+y6**3*y7**2)*y5**2)*y4+(2._ark/3._ark*y7**3+y6**2*y7)*y5**5+(5._ark/3._ark*y6**2*y7**3+y6**4*y7+2._ark/3._ark*y7**5)*y5**3
      dF(261) = (-y6**3/3._ark-y6*y7**2)*y4**5+(y6**2*y7-y7**3)*y5*y4**4+(-4._ark/3._ark*y6**3*y7**2-y6**5/3._ark-y6*y7**4-4._ark/3._ark*y5**2*y6**3)*y4**3+(-4._ark/3._ark*y5**3*y7**3+(y6**4*y7-y7**5)*y5)*y4**2+((-y6**3+y6*y7**2)*y5**4+(-y6**5+y6*y7**4)*y5**2)*y4+(-y6**2*y7-y7**3/3._ark)*y5**5+(-y6**4*y7-4._ark/3._ark*y6**2*y7**3-y7**5/3._ark)*y5**3
      dF(262) = (2._ark*y4**2*y5**3*y7+2._ark*y4**3*y5**2*y6+y4**4*y5*y7+y5**5*y7+y4*y5**4*y6+y4**5*y6)*y2**2+((y6**5+y6*y7**4+2._ark*y6**3*y7**2)*y4+(2._ark*y6**2*y7**3+y6**4*y7+y7**5)*y5)*y3**2
      dF(263) = (y4**3*y6**3+y5**3*y7**3+3._ark*y4*y5**2*y6*y7**2+3._ark*y4**2*y5*y6**2*y7)*y2**2+(y4**3*y6**3+y5**3*y7**3+3._ark*y4*y5**2*y6*y7**2+3._ark*y4**2*y5*y6**2*y7)*y3**2
      dF(264) = ((y6**2*y7**2+y6**4)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y2**2+(y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y3**2
      dF(265) = ((y6**5+y6*y7**4+2._ark*y6**3*y7**2)*y4+(2._ark*y6**2*y7**3+y6**4*y7+y7**5)*y5)*y2**2+(2._ark*y4**2*y5**3*y7+2._ark*y4**3*y5**2*y6+y4**4*y5*y7+y5**5*y7+y4*y5**4*y6+y4**5*y6)*y3**2
      dF(266) = (y4**4+y5**4+2._ark*y4**2*y5**2)*y2**4+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**4
      dF(267) = (2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y2**4+(2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3**4
      dF(268) = (y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y2**4+(y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y3**4
      dF(269) = ((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y2**4+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3**4
      dF(270) = (y4**5*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7+(y6*y7**4+2._ark*y6**3*y7**2+y6**5+y5**4*y6)*y4+y5**5*y7+(2._ark*y6**2*y7**3+y6**4*y7+y7**5)*y5)*y3*y2
      dF(271) = (-y4**3*y6*y7**2/2._ark+(-y7**3/2._ark+y6**2*y7)*y5*y4**2+(y6*y7**2-y6**3/2._ark)*y5**2*y4-y5**3*y6**2*y7/2._ark)*y3*y2
      dF(272) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+((y7**2+y6**2)*y5**2+y6**2*y7**2+y6**4)*y4**2+(2._ark*y5**3*y6*y7+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5)*y4+y5**4*y7**2+(y7**4+y6**2*y7**2)*y5**2)*y3*y2
      dF(273) = ((3._ark/2._ark*y6*y7**2+y6**3)*y4**3+3._ark/2._ark*y4**2*y5*y7**3+3._ark/2._ark*y4*y5**2*y6**3+(y7**3+3._ark/2._ark*y6**2*y7)*y5**3)*y3*y2
      dF(274) = ((y7**2-y6**2)*y4**4-4._ark*y4**3*y5*y6*y7+(y7**4-y6**4)*y4**2+(-4._ark*y5**3*y6*y7+(-4._ark*y6**3*y7-4._ark*y6*y7**3)*y5)*y4+(-y7**2+y6**2)*y5**4+(y6**4-y7**4)*y5**2)*y3*y2
      dF(275) = (y7**6+y6**6+y5**6+y4**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2+3._ark*y4**4*y5**2+3._ark*y4**2*y5**4)*y3*y2
      dF(276) = (y4**4+y5**4+2._ark*y4**2*y5**2)*y3*y2**3+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**3*y2
      dF(277) = (y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3*y2**3+((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3**3*y2
      dF(278) = (2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3*y2**3+(2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3**3*y2
      dF(279) = (y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y3*y2**3+(y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y3**3*y2
      dF(280) = ((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3*y2**3+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3**3*y2
      dF(281) = (y6**4+y7**4+2._ark*y6**2*y7**2)*y3*y2**3+(y4**4+y5**4+2._ark*y4**2*y5**2)*y3**3*y2
      dF(282) = (y5**2+y4**2)*y3*y2**5+(y7**2+y6**2)*y3**5*y2
      dF(283) = (y5*y7+y4*y6)*y3*y2**5+(y5*y7+y4*y6)*y3**5*y2
      dF(284) = (y7**2+y6**2)*y3*y2**5+(y5**2+y4**2)*y3**5*y2
      dF(285) = y2*y3**7+y2**7*y3
      dF(286) = (3._ark*y6**4*y7**2+y6**6+3._ark*y6**2*y7**4+y7**6)*y2**2+(3._ark*y4**2*y5**4+y5**6+3._ark*y4**4*y5**2+y4**6)*y3**2
      dF(287) = (3._ark*y4**2*y5**4+y5**6+3._ark*y4**4*y5**2+y4**6)*y2**2+(3._ark*y6**4*y7**2+y6**6+3._ark*y6**2*y7**4+y7**6)*y3**2
      dF(288) = ((y7**4-y6**4)*y4**2+(-4._ark*y6**3*y7-4._ark*y6*y7**3)*y5*y4+(y6**4-y7**4)*y5**2)*y2**2+((y7**2-y6**2)*y4**4-4._ark*y4**3*y5*y6*y7-4._ark*y4*y5**3*y6*y7+(-y7**2+y6**2)*y5**4)*y3**2
      dF(289) = (y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y2**2+(y4**3*y6*y7**2+(-2._ark*y6**2*y7+y7**3)*y5*y4**2+(-2._ark*y6*y7**2+y6**3)*y5**2*y4+y5**3*y6**2*y7)*y3**2
      dF(290) = (y4**4*y7**2-2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2-2._ark*y4*y5**3*y6*y7+y5**4*y6**2)*y2**2+((y7**4+y6**2*y7**2)*y4**2+(-2._ark*y6*y7**3-2._ark*y6**3*y7)*y5*y4+(y6**2*y7**2+y6**4)*y5**2)*y3**2
      dF(291) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y2**2+((y6**2*y7**2+y6**4)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y3**2
      dF(292) = (y4*y5*y6*y7-y5**2*y6**2/2._ark-y4**2*y7**2/2._ark)*y3**2*y2**2
      dF(293) = ((y7**2+y6**2)*y4**2+(y7**2+y6**2)*y5**2)*y3**2*y2**2
      dF(294) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y3**2*y2**2
      dF(295) = (2._ark*y4**2*y5**2+2._ark*y6**2*y7**2+y4**4+y5**4+y6**4+y7**4)*y3**2*y2**2
      dF(296) = (y5**2+y4**2)*y3**2*y2**4+(y7**2+y6**2)*y3**4*y2**2
      dF(297) = (y7**2+y6**2)*y3**2*y2**4+(y5**2+y4**2)*y3**4*y2**2
      dF(298) = y2**2*y3**6+y2**6*y3**2
      dF(299) = (y6**2+y7**2+y5**2+y4**2)*y3**3*y2**3
      dF(300) = (y5*y7+y4*y6)*y3**3*y2**3
      dF(301) = (y6**4+y7**4+2._ark*y6**2*y7**2)*y2**4+(y4**4+y5**4+2._ark*y4**2*y5**2)*y3**4
      dF(302) = (y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y2**4+((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3**4
      dF(303) = (y5*y7+y4*y6)*y3**2*y2**4+(y5*y7+y4*y6)*y3**4*y2**2
      dF(304) = y2**4*y3**4
      dF(305) = y2**5*y3**3+y2**3*y3**5
      dF(306) = (y7**2+y6**2)*y2**6+(y5**2+y4**2)*y3**6
      dF(307) = (y5*y7+y4*y6)*y2**6+(y5*y7+y4*y6)*y3**6
      dF(308) = (y5**2+y4**2)*y2**6+(y7**2+y6**2)*y3**6
      dF(309) = y2**8+y3**8
      dF(310) = (((y6**4-y7**4)*y4**2+(4._ark*y6**3*y7+4._ark*y6*y7**3)*y5*y4+(y7**4-y6**4)*y5**2)*y2+((-y7**2+y6**2)*y4**4+4._ark*y4**3*y5*y6*y7+4._ark*y4*y5**3*y6*y7+(y7**2-y6**2)*y5**4)*y3)*y1
      dF(311) = ((3._ark*y6**4*y7**2+y6**6+3._ark*y6**2*y7**4+y7**6)*y2+(3._ark*y4**2*y5**4+y5**6+3._ark*y4**4*y5**2+y4**6)*y3)*y1
      dF(312) = (((2._ark/3._ark*y6**3+y6*y7**2)*y4**3+y4**2*y5*y7**3+y4*y5**2*y6**3+(2._ark/3._ark*y7**3+y6**2*y7)*y5**3)*y2+((2._ark/3._ark*y6**3+y6*y7**2)*y4**3+y4**2*y5*y7**3+y4*y5**2*y6**3+(2._ark/3._ark*y7**3+y6**2*y7)*y5**3)*y3)*y1
      dF(313) = ((y4**4*y6**2+2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2+2._ark*y4*y5**3*y6*y7+y5**4*y7**2)*y2+((y6**2*y7**2+y6**4)*y4**2+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5*y4+(y7**4+y6**2*y7**2)*y5**2)*y3)*y1
      dF(314) = ((y4**4*y7**2-2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2-2._ark*y4*y5**3*y6*y7+y5**4*y6**2)*y2+((y7**4+y6**2*y7**2)*y4**2+(-2._ark*y6*y7**3-2._ark*y6**3*y7)*y5*y4+(y6**2*y7**2+y6**4)*y5**2)*y3)*y1
      dF(315) = ((3._ark*y4**2*y5**4+y5**6+3._ark*y4**4*y5**2+y4**6)*y2+(3._ark*y6**4*y7**2+y6**6+3._ark*y6**2*y7**4+y7**6)*y3)*y1
      dF(316) = (((y6**3*y7**2+y6**5/2._ark+y6*y7**4/2._ark)*y4+(y6**2*y7**3+y7**5/2._ark+y6**4*y7/2._ark)*y5)*y2+(y4*y5**4*y6/2._ark+y4**5*y6/2._ark+y4**4*y5*y7/2._ark+y5**5*y7/2._ark+y4**3*y5**2*y6+y4**2*y5**3*y7)*y3)*y1
      dF(317) = (((y7**4+y6**2*y7**2)*y4**2+(-2._ark*y6*y7**3-2._ark*y6**3*y7)*y5*y4+(y6**2*y7**2+y6**4)*y5**2)*y2+(y4**4*y7**2-2._ark*y4**3*y5*y6*y7+(y7**2+y6**2)*y5**2*y4**2-2._ark*y4*y5**3*y6*y7+y5**4*y6**2)*y3)*y1
      dF(318) = ((y4**2*y5*y6**2*y7+y4**3*y6**3/3._ark+y4*y5**2*y6*y7**2+y5**3*y7**3/3._ark)*y2+(y4**2*y5*y6**2*y7+y4**3*y6**3/3._ark+y4*y5**2*y6*y7**2+y5**3*y7**3/3._ark)*y3)*y1
      dF(319) = ((2._ark*y4**2*y5**3*y7+2._ark*y4**3*y5**2*y6+y4**4*y5*y7+y5**5*y7+y4*y5**4*y6+y4**5*y6)*y2+((y6**5+y6*y7**4+2._ark*y6**3*y7**2)*y4+(2._ark*y6**2*y7**3+y6**4*y7+y7**5)*y5)*y3)*y1
      dF(320) = ((y5*y7+y4*y6)*y3*y2**4+(y5*y7+y4*y6)*y3**4*y2)*y1
      dF(321) = (y2*y3**6+y2**6*y3)*y1
      dF(322) = ((y6**4+y7**4+2._ark*y6**2*y7**2)*y3*y2**2+(y4**4+y5**4+2._ark*y4**2*y5**2)*y3**2*y2)*y1
      dF(323) = (((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3*y2**2+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3**2*y2)*y1
      dF(324) = ((y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y3*y2**2+(y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y3**2*y2)*y1
      dF(325) = ((2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3*y2**2+(2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3**2*y2)*y1
      dF(326) = ((y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3*y2**2+((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3**2*y2)*y1
      dF(327) = ((y4**4+y5**4+2._ark*y4**2*y5**2)*y3*y2**2+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**2*y2)*y1
      dF(328) = ((y5**2+y4**2)*y3**2*y2**3+(y7**2+y6**2)*y3**3*y2**2)*y1
      dF(329) = ((y5*y7+y4*y6)*y3**2*y2**3+(y5*y7+y4*y6)*y3**3*y2**2)*y1
      dF(330) = ((y6**4+y7**4+2._ark*y6**2*y7**2)*y2**3+(y4**4+y5**4+2._ark*y4**2*y5**2)*y3**3)*y1
      dF(331) = (((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y2**3+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3**3)*y1
      dF(332) = ((y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y2**3+(y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y3**3)*y1
      dF(333) = ((2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y2**3+(2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3**3)*y1
      dF(334) = ((y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y2**3+((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3**3)*y1
      dF(335) = ((y4**4+y5**4+2._ark*y4**2*y5**2)*y2**3+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**3)*y1
      dF(336) = ((y7**2+y6**2)*y3**2*y2**3+(y5**2+y4**2)*y3**3*y2**2)*y1
      dF(337) = ((y7**2+y6**2)*y3*y2**4+(y5**2+y4**2)*y3**4*y2)*y1
      dF(338) = ((y5**2+y4**2)*y3*y2**4+(y7**2+y6**2)*y3**4*y2)*y1
      dF(339) = (y2**4*y3**3+y2**3*y3**4)*y1
      dF(340) = ((y7**2+y6**2)*y2**5+(y5**2+y4**2)*y3**5)*y1
      dF(341) = ((y5*y7+y4*y6)*y2**5+(y5*y7+y4*y6)*y3**5)*y1
      dF(342) = ((y5**2+y4**2)*y2**5+(y7**2+y6**2)*y3**5)*y1
      dF(343) = (y2**5*y3**2+y2**2*y3**5)*y1
      dF(344) = (y2**7+y3**7)*y1
      dF(345) = (y7**6+y6**6+y5**6+y4**6+3._ark*y6**2*y7**4+3._ark*y6**4*y7**2+3._ark*y4**4*y5**2+3._ark*y4**2*y5**4)*y1**2
      dF(346) = (y4**5*y6+y4**4*y5*y7+2._ark*y4**3*y5**2*y6+2._ark*y4**2*y5**3*y7+(y6*y7**4+2._ark*y6**3*y7**2+y6**5+y5**4*y6)*y4+y5**5*y7+(2._ark*y6**2*y7**3+y6**4*y7+y7**5)*y5)*y1**2
      dF(347) = (y4**4*y6**2+2._ark*y4**3*y5*y6*y7+((y7**2+y6**2)*y5**2+y6**2*y7**2+y6**4)*y4**2+(2._ark*y5**3*y6*y7+(2._ark*y6*y7**3+2._ark*y6**3*y7)*y5)*y4+y5**4*y7**2+(y7**4+y6**2*y7**2)*y5**2)*y1**2
      dF(348) = (-y4**3*y6*y7**2/2._ark+(-y7**3/2._ark+y6**2*y7)*y5*y4**2+(y6*y7**2-y6**3/2._ark)*y5**2*y4-y5**3*y6**2*y7/2._ark)*y1**2
      dF(349) = ((3._ark/2._ark*y6*y7**2+y6**3)*y4**3+3._ark/2._ark*y4**2*y5*y7**3+3._ark/2._ark*y4*y5**2*y6**3+(y7**3+3._ark/2._ark*y6**2*y7)*y5**3)*y1**2
      dF(350) = ((y7**2-y6**2)*y4**4-4._ark*y4**3*y5*y6*y7+(y7**4-y6**4)*y4**2+(-4._ark*y5**3*y6*y7+(-4._ark*y6**3*y7-4._ark*y6*y7**3)*y5)*y4+(-y7**2+y6**2)*y5**4+(y6**4-y7**4)*y5**2)*y1**2
      dF(351) = ((y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y2**2+((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3**2)*y1**2
      dF(352) = (y3**6+y2**6)*y1**2
      dF(353) = (2._ark*y4**2*y5**2+2._ark*y6**2*y7**2+y4**4+y5**4+y6**4+y7**4)*y3*y2*y1**2
      dF(354) = (y4*y5*y6*y7+y5**2*y7**2/2._ark+y4**2*y6**2/2._ark)*y3*y2*y1**2
      dF(355) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y3*y2*y1**2
      dF(356) = ((y7**2+y6**2)*y4**2+(y7**2+y6**2)*y5**2)*y3*y2*y1**2
      dF(357) = (y2*y3**5+y2**5*y3)*y1**2
      dF(358) = ((y6**4+y7**4+2._ark*y6**2*y7**2)*y2**2+(y4**4+y5**4+2._ark*y4**2*y5**2)*y3**2)*y1**2
      dF(359) = (((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y2**2+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3**2)*y1**2
      dF(360) = ((y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y2**2+(y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y3**2)*y1**2
      dF(361) = ((2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y2**2+(2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3**2)*y1**2
      dF(362) = ((y4**4+y5**4+2._ark*y4**2*y5**2)*y2**2+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3**2)*y1**2
      dF(363) = (y5*y7+y4*y6)*y3**2*y2**2*y1**2
      dF(364) = (y6**2+y7**2+y5**2+y4**2)*y3**2*y2**2*y1**2
      dF(365) = (y2**2*y3**4+y2**4*y3**2)*y1**2
      dF(366) = ((y7**2+y6**2)*y3*y2**3+(y5**2+y4**2)*y3**3*y2)*y1**2
      dF(367) = ((y5*y7+y4*y6)*y3*y2**3+(y5*y7+y4*y6)*y3**3*y2)*y1**2
      dF(368) = ((y5**2+y4**2)*y3*y2**3+(y7**2+y6**2)*y3**3*y2)*y1**2
      dF(369) = y1**2*y2**3*y3**3
      dF(370) = ((y7**2+y6**2)*y2**4+(y5**2+y4**2)*y3**4)*y1**2
      dF(371) = ((y5*y7+y4*y6)*y2**4+(y5*y7+y4*y6)*y3**4)*y1**2
      dF(372) = ((y5**2+y4**2)*y2**4+(y7**2+y6**2)*y3**4)*y1**2
      dF(373) = ((y4**4+y5**4+2._ark*y4**2*y5**2)*y2+(y6**4+y7**4+2._ark*y6**2*y7**2)*y3)*y1**3
      dF(374) = ((y6**4+y7**4+2._ark*y6**2*y7**2)*y2+(y4**4+y5**4+2._ark*y4**2*y5**2)*y3)*y1**3
      dF(375) = (((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y2+(y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y3)*y1**3
      dF(376) = ((y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y2+(y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y3)*y1**3
      dF(377) = ((2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y2+(2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y3)*y1**3
      dF(378) = ((y4**3*y6+y5**3*y7+y4**2*y5*y7+y4*y5**2*y6)*y2+((y6**3+y6*y7**2)*y4+(y7**3+y6**2*y7)*y5)*y3)*y1**3
      dF(379) = ((y7**2+y6**2)*y3*y2**2+(y5**2+y4**2)*y3**2*y2)*y1**3
      dF(380) = ((y5*y7+y4*y6)*y3*y2**2+(y5*y7+y4*y6)*y3**2*y2)*y1**3
      dF(381) = ((y5**2+y4**2)*y3*y2**2+(y7**2+y6**2)*y3**2*y2)*y1**3
      dF(382) = ((y7**2+y6**2)*y2**3+(y5**2+y4**2)*y3**3)*y1**3
      dF(383) = ((y5*y7+y4*y6)*y2**3+(y5*y7+y4*y6)*y3**3)*y1**3
      dF(384) = ((y5**2+y4**2)*y2**3+(y7**2+y6**2)*y3**3)*y1**3
      dF(385) = (y2**3*y3**2+y2**2*y3**3)*y1**3
      dF(386) = (y2**4*y3+y2*y3**4)*y1**3
      dF(387) = (y2**5+y3**5)*y1**3
      dF(388) = (y5**2*y6**2+y4**2*y7**2-2._ark*y4*y5*y6*y7)*y1**4
      dF(389) = (2._ark*y4*y5*y6*y7+y5**2*y7**2+y4**2*y6**2)*y1**4
      dF(390) = (y4**3*y6+y4**2*y5*y7+(y6*y7**2+y6**3+y5**2*y6)*y4+y5**3*y7+(y7**3+y6**2*y7)*y5)*y1**4
      dF(391) = (2._ark*y4**2*y5**2+2._ark*y6**2*y7**2+y4**4+y5**4+y6**4+y7**4)*y1**4
      dF(392) = ((y7**2+y6**2)*y2**2+(y5**2+y4**2)*y3**2)*y1**4
      dF(393) = (y5*y7+y4*y6)*y3*y2*y1**4
      dF(394) = (y6**2+y7**2+y5**2+y4**2)*y3*y2*y1**4
      dF(395) = ((y5*y7+y4*y6)*y2**2+(y5*y7+y4*y6)*y3**2)*y1**4
      dF(396) = ((y5**2+y4**2)*y2**2+(y7**2+y6**2)*y3**2)*y1**4
      dF(397) = y1**4*y2**2*y3**2
      dF(398) = (y2**3*y3+y2*y3**3)*y1**4
      dF(399) = (y2**4+y3**4)*y1**4
      dF(400) = ((y7**2+y6**2)*y2+(y5**2+y4**2)*y3)*y1**5
      dF(401) = ((y5*y7+y4*y6)*y2+(y5*y7+y4*y6)*y3)*y1**5
      dF(402) = ((y5**2+y4**2)*y2+(y7**2+y6**2)*y3)*y1**5
      dF(403) = (y2*y3**2+y2**2*y3)*y1**5
      dF(404) = (y3**3+y2**3)*y1**5
      dF(405) = (y5*y7+y4*y6)*y1**6
      dF(406) = (y6**2+y7**2+y5**2+y4**2)*y1**6
      dF(407) = y1**6*y2*y3
      dF(408) = (y2**2+y3**2)*y1**6
      dF(409) = (y2+y3)*y1**7
      dF(410) = y1**8
      !
  end subroutine  potC2H2_D8h_diff_V


 function MLpoten_c2h2_7_xyz(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
    !
    integer(ik),parameter :: n = 578
    integer(ik) :: i,k,nmax
    real(ark) :: dF(n)
      !
      integer(ik)  ::  i1,i2,i3,i4,i5,i6,k_ind(6)
      real(ark)    :: x1,x2,x3,x4,x5,x6,e1,e2,e4,e6,vpot,cphi,q(6),y(6),a1,a2,pd,rc1c2,rc1h1,rc2h2,delta1x,delta1y,delta2x,delta2y,tau
      real(ark)    :: alpha1,alpha2,sinalpha2,sinalpha1,tau1,tau2,b1(3),b0(3),b2(3),t1,t0,t2,w1(3),w2(3),cosalpha2,sindelta1x,sindelta1y,sindelta2x,sindelta2y,y1,y2,y3,y4,y5,y6,y7,c1(3),c0(3),c2(3)
      real(ark)    :: r_na(4,3)
      integer(ik)  :: Nangles 
      !
      character(len=cl)  :: txt = 'MLpoten_c2h2_7'
      !
      Nangles = molec%Nangles
      !
      pd=pi/180.0_ark
      e1=molec%force(1)
      e2=molec%force(2)
      e4=pi
      e6=pi
      !
      a1 = molec%force(3)
      a2 = molec%force(4)
      !
      x1    = local(1)
      x2    = local(2)
      x3    = local(3)
      !
      if (molec%zmatrix(3)%connect(4)==101) then 
         !
         stop 'MLpoten_c2h2_7_xyz 101 is not implemented'
         !
         y4 = local(4)
         y5 = local(5)
         !
         y6 = local(6)
         y7 = local(7)
         !
      elseif (molec%zmatrix(3)%connect(4)>=102.and.trim(molec%coords_transform)=='R-R1-R2-TX-TY-TX-TY') then
         !
         call MLfromlocal2cartesian(1_ik,local,r_na)
         !
         c1(:) = r_na(3,:)-r_na(1,:)
         c0(:) = r_na(2,:)-r_na(1,:)
         c2(:) = r_na(4,:)-r_na(2,:)
         !
         x2 =  sqrt(sum(c1(:)**2))
         x1 =  sqrt(sum(c0(:)**2))
         x3 =  sqrt(sum(c2(:)**2))
         !
         b1 =  c1(:)/x2
         b0 =  c0(:)/x1
         b2 =  c2(:)/x3
         !
         !y2 = -c1(3)-e2
         !y3 =  c2(3)-e2
         !
         x2 = -c1(3)
         x3 =  c2(3)
         !
         y4 = c1(1)
         y5 = c1(2)
         y6 = c2(1)
         y7 = c2(2)
         !
         w1(:) = MLvector_product(b1,b0)
         w2(:) = MLvector_product(b2,b0)
         !
         cosalpha2 = sum(b0(:)*b2(:))
         !
         alpha2 = aacos(-cosalpha2,txt)
         !
         !if (alpha1>pi*0.5d0) w1 = -w1
         !if (alpha2>pi*0.5d0) w2 = -w2
         !
         sindelta1x = -w1(1)
         sindelta1y = -w1(2)
         sindelta2x = -w2(1)
         sindelta2y = -w2(2)
         !
         !sindelta1x = -w1(2)
         !sindelta1y = -w1(1)
         !sindelta2x = -w2(2)
         !sindelta2y =  w2(1)
         !
         !y4 = aasin(sindelta1x,txt)
         !y5 = aasin(sindelta1y,txt)
         !y6 = aasin(sindelta2x,txt)
         !y7 = aasin(sindelta2y,txt)
         !
         y4 =-sindelta1y
         y5 = sindelta1x
         y6 = sindelta2y
         y7 =-sindelta2x
         !
      elseif (molec%zmatrix(3)%connect(4)==103.and.trim(molec%coords_transform)=='R-Z1-Z2-X-Y-X-Y') then 
         !
         y4 = local(4)
         y5 = local(5)
         !
         y6 = local(6)
         y7 = local(7)
         !
         call MLfromlocal2cartesian(1_ik,local,r_na)
         !
         c1(:) = r_na(3,:)-r_na(1,:)
         c0(:) = r_na(2,:)-r_na(1,:)
         c2(:) = r_na(4,:)-r_na(2,:)
         !
         x2 = -c1(3)
         x3 =  c2(3)
         !
         y4 = c1(1)
         y5 = c1(2)
         y6 = c2(1)
         y7 = c2(2)
         !
      elseif (molec%zmatrix(3)%connect(4)==103.and.trim(molec%coords_transform)=='R-R1-R2-X-Y-X-Y') then 
         !
         y4 = local(4)
         y5 = local(5)
         !
         y6 = local(6)
         y7 = local(7)
         !
         call MLfromlocal2cartesian(1_ik,local,r_na)
         !
         c1(:) = r_na(3,:)-r_na(1,:)
         c0(:) = r_na(2,:)-r_na(1,:)
         c2(:) = r_na(4,:)-r_na(2,:)
         !
         x2 = -c1(3)
         x3 =  c2(3)
         !
         y4 = c1(1)
         y5 = c1(2)
         y6 = c2(1)
         y7 = c2(2)
         !
      else 
        !
        write(out,"('MLpoten_c2h2_7_xyz: only designed for zmatrix-connect( =103 ',i)") molec%zmatrix(3)%connect(4)
        stop 'only designed for zmat=103' 
        !
      endif
      !
      !y4    = local(4)
      !y5    = local(5)
      !y6    = local(6)
      !y7    = local(7)
      !
      y1=1.0_ark-exp(-a1*(x1-e1))
      !
      !y2=1.0_ark-exp(-a2*(x2-e2))
      !y3=1.0_ark-exp(-a2*(x3-e2))
      !
      y2=( x2-e2)
      y3=( x3-e2)
      !
      call potC2H2_D8h_diff_V(n,y1,y2,y3,y4,y5,y6,y7,dF) 
      !
      f = 0
      !
      nmax = min(size(force),molec%parmax)
      !
      do i = 5,nmax
        !
        !k = molec%pot_ind(1,i)
        !
        f = f + force(i)*dF(i)
       !
      enddo
      !
 end function MLpoten_c2h2_7_xyz



 function MLpoten_c2h2_7(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
    !
    integer(ik),parameter :: n = 415
    integer(ik) :: i,k,nmax
    real(ark) :: dF(n)
      !
      integer(ik)  ::  i1,i2,i3,i4,i5,i6,k_ind(6)
      real(ark)    :: x1,x2,x3,x4,x5,x6,e1,e2,e4,e6,vpot,cphi,q(6),y(6),a1,a2,pd,rc1c2,rc1h1,rc2h2,delta1x,delta1y,delta2x,delta2y,tau
      real(ark)    :: alpha1,alpha2,sinalpha2,sinalpha1,tau1,tau2,b1(3),b0(3),b2(3),t1,t0,t2,w1(3),w2(3),cosalpha2,sindelta1x,sindelta1y,sindelta2x,sindelta2y,y1,y2,y3,y4,y5,y6,y7,c1(3),c0(3),c2(3)
      real(ark)    :: r_na(4,3)
      integer(ik)  :: Nangles 
      !
      character(len=cl)  :: txt = 'MLpoten_c2h2_7'
      !
      Nangles = molec%Nangles
      !
      pd=pi/180.0_ark
      e1=molec%force(1)
      e2=molec%force(2)
      e4=pi
      e6=pi
      !
      a1 = molec%force(3)
      a2 = molec%force(4)
      !
      x1    = local(1)
      x2    = local(2)
      x3    = local(3)
      !
      if (molec%zmatrix(3)%connect(4)==101) then 
         !
         y4 = local(4)
         y5 = local(5)
         !
         y6 = local(6)
         y7 = local(7)
         !
      elseif (molec%zmatrix(3)%connect(4)==102) then
         !
         call MLfromlocal2cartesian(1_ik,local,r_na)
         !
         b1(:) = r_na(3,:)-r_na(1,:)
         b0(:) = r_na(2,:)-r_na(1,:)
         b2(:) = r_na(4,:)-r_na(2,:)
         !
         x2 =  sqrt(sum(b1(:)**2))
         x1 =  sqrt(sum(b0(:)**2))
         x3 =  sqrt(sum(b2(:)**2))
         !
         b1 =  b1(:)/x2
         b0 =  b0(:)/x1
         b2 =  b2(:)/x3
         !
         w1(:) = MLvector_product(b1,b0)
         w2(:) = MLvector_product(b2,b0)
         !
         !sindelta1x = w1(1)
         !sindelta1y = w1(2)
         !sindelta2x = w2(1)
         !sindelta2y = w2(2)
         !
         !sindelta1x = -w1(2)
         !sindelta1y = -w1(1)
         !sindelta2x = -w2(2)
         !sindelta2y =  w2(1)
         !
         !y4 = aasin(sindelta1x,txt)
         !y5 = aasin(sindelta1y,txt)
         !y6 = aasin(sindelta2x,txt)
         !y7 = aasin(sindelta2y,txt)
         !
         y4 =-w1(2)
         y5 = w1(1)
         y6 =-w2(2)
         y7 = w2(1)
         !
      else 
         !
         call MLfromlocal2cartesian(1_ik,local,r_na)
         !
         !
         b1(:) = r_na(3,:)-r_na(1,:)
         b0(:) = r_na(2,:)-r_na(1,:)
         b2(:) = r_na(4,:)-r_na(2,:)
         !
         x2 =  sqrt(sum(b1(:)**2))
         x1 =  sqrt(sum(b0(:)**2))
         x3 =  sqrt(sum(b2(:)**2))
         !
         b1 =  b1(:)/x2
         b0 =  b0(:)/x1
         b2 =  b2(:)/x3
         !
         w1(:) = MLvector_product(b1,b0)
         w2(:) = MLvector_product(b2,b0)
         !
         !sindelta1x = -w1(2)
         !sindelta1y = -w1(1)
         !sindelta2x = -w2(2)
         !sindelta2y =  w2(1)
         !
         !y4 = aasin(sindelta1x,txt)
         !y5 = aasin(sindelta1y,txt)
         !y6 = aasin(sindelta2x,txt)
         !y7 = aasin(sindelta2y,txt)
         !
         y4 =-w1(2)
         y5 = w1(1)
         y6 =-w2(2)
         y7 = w2(1)
         !
      endif
      !
      !y4    = local(4)
      !y5    = local(5)
      !y6    = local(6)
      !y7    = local(7)
      !
      y1=1.0_ark-exp(-a1*(x1-e1))
      y2=1.0_ark-exp(-a2*(x2-e2))
      y3=1.0_ark-exp(-a2*(x3-e2))
      !
      call potC2H2_D8h_diff_V(n,y1,y2,y3,y4,y5,y6,y7,dF)
      !
      f = 0
      !
      nmax = min(size(force),molec%parmax)
      !
      do i = 5,nmax
        !
        !k = molec%pot_ind(1,i)
        !
        f = f + force(i)*dF(i)
       !
      enddo
      !
 end function MLpoten_c2h2_7



function MLpoten_c2h2_streymills(ncoords,natoms,local,xyz,force) result(f)
   !
   implicit none
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          ::  i,i1,i2,i3,i4,i5,i6,k_ind(6)
   real(ark)    :: x1,x2,x3,x4,x5,x6,e1,e2,e3,vpot,cphi,q(6),y(6),pd,rc1c2,rc1h1,rc2h2,delta1x,delta1y,delta2x,delta2y,tau,y1,y2,y3,y4,y5,y6
   real(ark)    :: alpha1,alpha2,sinalpha2,sinalpha1,tau1,tau2,v1(3),v2(3),v3(3),v12(3),v31(3),r21,r31,n3(3),sintau,costau,cosalpha1,cosalpha2,e4
   character(len=cl)  :: txt = 'MLpoten_c2h2_8'
   !
   integer(ik)  :: Nangles
   integer(ik),parameter :: n = 44
   integer(ik) :: k,nmax
   real(ark) :: dF(n)
      !
      Nangles = molec%Nangles
      !
      e1=force(1)
      e2=force(2)
      e4=pi
      !
      rc1c2    = local(1)
      rc1h1    = local(2)
      rc2h2    = local(3)
      !
      if (molec%zmatrix(3)%connect(4)==101) then 
         !
         delta1x = local(4)
         delta1y = local(5)
         !
         delta2x = local(6)
         delta2y = local(7)
         !
         tau1 = pi-atan2(sin(delta1y),-sin(delta1x))
         tau2 = pi-atan2(sin(delta2y),-sin(delta2x))
         !
         !tau1 = pi-atan2((delta1y),-(delta1x))
         !tau2 = pi-atan2((delta2y),-(delta2x))
         !
         tau = tau2-tau1
         !
         sinalpha1 = sqrt(sin(delta1x)**2+sin(delta1y)**2)
         sinalpha2 = sqrt(sin(delta2x)**2+sin(delta2y)**2)
         !
         !sinalpha1 = sqrt((delta1x)**2+(delta1y)**2)
         !sinalpha2 = sqrt((delta2x)**2+(delta2y)**2)
         !
         alpha1 = pi-asin(sinalpha1)
         alpha2 = pi-asin(sinalpha2)
         !
         x1 = rc1c2
         x3 = rc1h1
         x2 = rc2h2
         x5 = alpha1
         x4 = alpha2
         !
         costau = cos(tau)
         !
      elseif (molec%zmatrix(3)%connect(4)==103) then
         !
         v1(:) = xyz(2,:)-xyz(1,:)
         v2(:) = xyz(3,:)-xyz(1,:)
         v3(:) = xyz(4,:)-xyz(2,:)
         !
         x1 = sqrt(sum(v1(:)**2))
         x2 = sqrt(sum(v2(:)**2))
         x3 = sqrt(sum(v3(:)**2))
         !
         cosalpha1 = sum( v1(:)*v2(:))/(x1*x2)
         cosalpha2 = sum(-v1(:)*v3(:))/(x1*x3)
         !
         x4 = aacos(cosalpha1,txt)
         x5 = aacos(cosalpha2,txt)
         !
         v12(:) = MLvector_product(v2(:),v1(:))
         v31(:) = MLvector_product(v3(:),v1(:))
         !
         v12(:) = v12(:)/(x1*x2)
         v31(:) = v31(:)/(x1*x3)
         !
         x6 = sum(v12*v31)
         !
         r21 = sqrt(sum(v12(:)**2))
         r31 = sqrt(sum(v31(:)**2))
         !
         tau = 0
         !
         costau = 1.0_ark
         !
         x6 = 0
         !
         if (r21>sqrt(small_).and.r31>sqrt(small_)) then
           !
           v12(:) = v12(:)/r21
           v31(:) = v31(:)/r31
           !
           n3(:) = MLvector_product(v12(:),v31(:))
           !
           sintau =  sqrt( sum( n3(:)**2 ) )
           costau =  sum( v12(:)*v31(:) )
           !
           tau = atan2(sintau,costau)
           tau = aacos(costau,txt)
           !
           x6 = tau
           !
           !
         endif
         !
      else
        !
        write(out,"('MLpoten_c2h2_8: only designed for zmatrix-connect( =103 ',i)") molec%zmatrix(3)%connect(4)
        stop 'only designed for zmat=103'
        !
      endif
      !
      y1=(x1-e1)
      y2=(x2-e2)
      y3=(x3-e2)
      !
      y4=x4-e4
      y5=x5-e4
      y6=costau
      sintau = sin(tau)
      !
      dF(1:4) = 0
      !
      dF(5) = 0.5_ark*((y2*y2)+(y3*y3))
      dF(6) = 0.5_ark*(y1*y1)
      dF(7) = y2*y3
      dF(8) = (y2+y3)*y1
      dF(9) = (1.0_ark/6.0_ark)*((y2*y2*y2)+(y3*y3*y3))
      dF(10) = 0.5_ark*((y2*y2*y3)+(y3*y3*y2))
      dF(11) = 0.5_ark*(((y2*y2)+(y3*y3))*y1)
      dF(12) = y2*y3*y1
      dF(13) = 0.5_ark*((y2+y3)*y1*y1)
      dF(14) = (1.0_ark/6.0_ark)*(y1*y1*y1)
      dF(15) = (1.0_ark/24.0_ark)*((y2*y2*y2*y2)+(y3*y3*y3*y3))
      dF(16) = (1.0_ark/6.0_ark)*((y2*y2*y2*y3)+(y2*y3*y3*y3))
      dF(17) = (1.0_ark/4.0_ark)*(y2*y2*y3*y3)
      dF(18) = (1.0_ark/6.0_ark)*((y2*y2*y2)+(y3*y3*y3))*y1
      dF(19) = 0.5_ark*((y2*y2*y3)+(y3*y3*y2))*y1
      dF(20) = 0.25_ark*(((y2*y2)+(y3*y3))*y1*y1)
      dF(21) = 0.5_ark*y2*y3*y1*y1
      dF(22) = (1.0_ark/6.0_ark)*((y2+y3)*y1*y1*y1)
      dF(23) = (1.0_ark/24.0_ark)*(y1*y1*y1*y1)
      dF(24) = 0.25_ark*((y2*y2*y4*y4)+(y3*y3*y5*y5))
      dF(25) = 0.25_ark*((y2*y2*y5*y5)+(y3*y3*y4*y4))
      dF(26) = 0.5_ark*((y2*y2)+(y3*y3))*y4*y5*y6
      dF(27) = 0.5_ark*((y4*y4)+(y5*y5))
      dF(28) = y4*y5*y6
      dF(29) = 0.5_ark*((y2*y4*y4)+(y3*y5*y5))
      dF(30) = (y2+y3)*y4*y5*y6
      dF(31) = 0.5_ark*((y2*y5*y5)+(y3*y4*y4))
      dF(32) = 0.5_ark*(y1*((y4*y4)+(y5*y5)))
      dF(33) = y1*y4*y5*y6
      dF(34) = 0.5_ark*y2*y3*((y4*y4)+(y5*y5))
      dF(35) = y2*y3*y4*y5*y6
      dF(36) = 0.5_ark*(((y2*y4*y4)+(y3*y5*y5))*y1)
      dF(37) = 0.5_ark*(((y2*y5*y5)+(y3*y4*y4))*y1)
      dF(38) = (y2+y3)*y4*y5*y6*y1
      dF(39) = 0.25_ark*(y1*y1*((y4*y4)+(y5*y5)))
      dF(40) = y1*y1*y4*y5*y6
      dF(41) = (1.0_ark/24.0_ark)*((y4*y4*y4*y4)+(y5*y5*y5*y5))
      dF(42) = (1.0_ark/6.0_ark)*((y4*y4)+(y5*y5))*y4*y5*y6
      dF(43) = 0.25_ark*(y4*y4*y5*y5)
      dF(44) = -0.25_ark*(y4*y4*y5*y5)*sintau
      !
      nmax = min(size(force),molec%parmax,size(dF))
      !
      f = sum(force(5:nmax)*dF(5:nmax))
      !
      f = f*1.0e-11/planck/vellgt
      !
end function MLpoten_c2h2_streymills



end module pot_abcd

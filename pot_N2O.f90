!
!  This unit is for a user defined potential: N2O PES 
!
module pot_user
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xyz

  implicit none

  public MLdipole,MLpoten,ML_MEP

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  !
  integer(ik) :: npol(5),order(5,5)
  real(ark)   :: gam1(5,5),x0(5,5),r0(5,5),beta(5,5)
  integer(ik) :: np(5)
  !
  contains
  !
  function ML_MEP(dim,rho)  result(f)

   integer(ik),intent(in) ::  dim
   real(ark),intent(in)   ::  rho
   real(ark)              ::  f(dim)
   !
   if (dim/=3) stop 'Illegal size of the function - must be 3'
   !
   f(:) = molec%local_eq(:)
   f(molec%Ncoords) = rho

  end function ML_MEP



 !
 ! Defining potential energy function (built for SO2)

 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f,f1,f2
   integer(ik)            ::  nparam
   !
   f1 = MLpoten_n2opotlongrange(ncoords,natoms,local,xyz,force)
   !
   nparam = int(force(1))
   !
   f2 = MLpoten_xyz_N2O_Zobov(ncoords,natoms,local,xyz,force(nparam+1:))
   !
   f = f1+f2
   !
   !
 end function MLpoten
 !


 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
    !
    real(ark) :: xyz_(natoms,3), r1, r2,  n1(3), n2(3), n3(3), tmat(3,3), CN_CM(3),&
                 x(natoms,3),xyz0(natoms,3),cos_theta,alpha,mu(3),u1(3),u2(3),u3(3),bigr,smallr

    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      xyz0 = MLloc2pqr_xyz(local)
      !
      x(1,:) = xyz0(2,:) - xyz0(1,:)
      x(2,:) = xyz0(3,:) - xyz0(1,:)
      !
    else
      !
      x(1,:) = xyz(2,:) - xyz(1,:)
      x(2,:) = xyz(3,:) - xyz(1,:)
      !
    endif
    !
    r1 = sqrt(sum(x(1,:)**2))
    r2 = sqrt(sum(x(2,:)**2))
    !
    !
    !
    stop 'MLdipole N2O is undefined'
    !
  end subroutine MLdipole




  function MLpoten_xyz_N2O_Zobov(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f
   !
   real(ark)            :: r1,r2,alpha,xcos,v,v1,v2,v3,v4,v5,v6,v7,v8

   real(ark)            :: aa1,aa2,re1,re2,alphae,xst,xs1,xs2,v0
   integer(ik)          :: N,Nparams
   real(ark)            :: rho,xp,vp1
   real(ark)            :: ReNN_ang,ReNO_ang
   !
   if (verbose>=6) write(out,"('MLpoten_xyz_tyuterev/start')")
   !
   r1  = local(1) ;  r2    = local(2) ;  alpha = local(3)
   !
   if (molec%Nangles==0) then
     stop 'MLpoten_xyz_N2O_Zobov: ilegal number of bond anlges 0'
   endif 
   !
   ReNN_ang    = force( 1)
   ReNO_ang    = force( 2)
   alphae      = force( 3)*pi/180.0_ark
   aa1         = force( 4)
   aa2         = force( 5)
   !
   ! calculate potential energy function values
   !    
   !
   Nparams = size(force)
   !
   rho=pi-alpha
   !   
   xst= 1.0_ark - cos(rho)
   xs1= 1.0_ark - exp(-aa1*(r1-ReNN_ang))
   xs2= 1.0_ark - exp(-aa2*(r2-ReNO_ang))
   !
   N = 5
   !
  v0= force(N+1)  *xs1**0*xs2**0*xst**0
 vp1= force(N+2)  *xs1**1*xs2**0*xst**0& 
     +force(N+3)  *xs1**0*xs2**1*xst**0& 
     +force(N+4)  *xs1**0*xs2**0*xst**1& 
     +force(N+5)  *xs1**2*xs2**0*xst**0& 
     +force(N+6)  *xs1**1*xs2**1*xst**0& 
     +force(N+7)  *xs1**1*xs2**0*xst**1& 
     +force(N+8)  *xs1**0*xs2**2*xst**0& 
     +force(N+9)  *xs1**0*xs2**1*xst**1& 
     +force(N+10) *xs1**0*xs2**0*xst**2     ! end of 2
 
 if (Nparams>N+10) then 
   vp1=vp1+&     
     +force(N+11) *xs1**3*xs2**0*xst**0& 
     +force(N+12) *xs1**2*xs2**1*xst**0& 
     +force(N+13) *xs1**2*xs2**0*xst**1& 
     +force(N+14) *xs1**1*xs2**2*xst**0& 
     +force(N+15) *xs1**1*xs2**1*xst**1& 
     +force(N+16) *xs1**1*xs2**0*xst**2& 
     +force(N+17) *xs1**0*xs2**3*xst**0& 
     +force(N+18) *xs1**0*xs2**2*xst**1& 
     +force(N+19) *xs1**0*xs2**1*xst**2& 
     +force(N+20) *xs1**0*xs2**0*xst**3&     !  end of 3
     +force(N+21) *xs1**4*xs2**0*xst**0& 
     +force(N+22) *xs1**3*xs2**1*xst**0& 
     +force(N+23) *xs1**3*xs2**0*xst**1& 
     +force(N+24) *xs1**2*xs2**2*xst**0& 
     +force(N+25) *xs1**2*xs2**1*xst**1& 
     +force(N+26) *xs1**2*xs2**0*xst**2& 
     +force(N+27) *xs1**1*xs2**3*xst**0& 
     +force(N+28) *xs1**1*xs2**2*xst**1&
     +force(N+29) *xs1**1*xs2**1*xst**2& 
     +force(N+30) *xs1**1*xs2**0*xst**3& 
     +force(N+31) *xs1**0*xs2**4*xst**0& 
     +force(N+32) *xs1**0*xs2**3*xst**1& 
     +force(N+33) *xs1**0*xs2**2*xst**2& 
     +force(N+34) *xs1**0*xs2**1*xst**3& 
     +force(N+35) *xs1**0*xs2**0*xst**4   !  end of 4
 endif
 !
 if (Nparams>N+35) then 
   vp1=vp1+&          
     +force(N+36) *xs1**5*xs2**0*xst**0& 
     +force(N+37) *xs1**4*xs2**1*xst**0& 
     +force(N+38) *xs1**4*xs2**0*xst**1& 
     +force(N+39) *xs1**3*xs2**2*xst**0& 
     +force(N+40) *xs1**3*xs2**1*xst**1& 
     +force(N+41) *xs1**3*xs2**0*xst**2& 
     +force(N+42) *xs1**2*xs2**3*xst**0& 
     +force(N+43) *xs1**2*xs2**2*xst**1& 
     +force(N+44) *xs1**2*xs2**1*xst**2& 
     +force(N+45) *xs1**2*xs2**0*xst**3& 
     +force(N+46) *xs1**1*xs2**4*xst**0& 
     +force(N+47) *xs1**1*xs2**3*xst**1& 
     +force(N+48) *xs1**1*xs2**2*xst**2& 
     +force(N+49) *xs1**1*xs2**1*xst**3& 
     +force(N+50) *xs1**1*xs2**0*xst**4&  !
     +force(N+51) *xs1**0*xs2**5*xst**0& 
     +force(N+52) *xs1**0*xs2**4*xst**1& 
     +force(N+53) *xs1**0*xs2**3*xst**2& 
     +force(N+54) *xs1**0*xs2**2*xst**3& 
     +force(N+55) *xs1**0*xs2**1*xst**4& 
     +force(N+56) *xs1**0*xs2**0*xst**5   !  end of 5
 endif
 !
 if (Nparams>N+56) then 
   vp1=vp1+&          
     +force(N+57) *xs1**6*xs2**0*xst**0& 
     +force(N+58) *xs1**5*xs2**1*xst**0& 
     +force(N+59) *xs1**5*xs2**0*xst**1& 
     +force(N+60) *xs1**4*xs2**2*xst**0& 
     +force(N+61) *xs1**4*xs2**1*xst**1& 
     +force(N+62) *xs1**4*xs2**0*xst**2& 
     +force(N+63) *xs1**3*xs2**3*xst**0& 
     +force(N+64) *xs1**3*xs2**2*xst**1& 
     +force(N+65) *xs1**3*xs2**1*xst**2& 
     +force(N+66) *xs1**3*xs2**0*xst**3& 
     +force(N+67) *xs1**2*xs2**4*xst**0& 
     +force(N+68) *xs1**2*xs2**3*xst**1& 
     +force(N+69) *xs1**2*xs2**2*xst**2& 
     +force(N+70) *xs1**2*xs2**1*xst**3& 
     +force(N+71) *xs1**2*xs2**0*xst**4& 
     +force(N+72) *xs1**1*xs2**5*xst**0& 
     +force(N+73) *xs1**1*xs2**4*xst**1& 
     +force(N+74) *xs1**1*xs2**3*xst**2& 
     +force(N+75) *xs1**1*xs2**2*xst**3& 
     +force(N+76) *xs1**1*xs2**1*xst**4& 
     +force(N+77) *xs1**1*xs2**0*xst**5& 
     +force(N+78) *xs1**0*xs2**6*xst**0& 
     +force(N+79) *xs1**0*xs2**5*xst**1& 
     +force(N+80) *xs1**0*xs2**4*xst**2& 
     +force(N+81) *xs1**0*xs2**3*xst**3& 
     +force(N+82) *xs1**0*xs2**2*xst**4& 
     +force(N+83) *xs1**0*xs2**1*xst**5& 
     +force(N+84) *xs1**0*xs2**0*xst**6   !  enf of 6
 endif
 f=v0+vp1

 end function MLpoten_xyz_N2O_Zobov


  function MLpoten_n2opotlongrange(ncoords,natoms,local,xyz,force) result(f)

     integer(ik),intent(in) ::  ncoords,natoms
     real(ark),intent(in)   ::  local(ncoords)
     real(ark),intent(in)   ::  xyz(natoms,3)
     real(ark),intent(in)   ::  force( :)
     real(ark)              ::  f
     integer(ik) :: iopttmp(3)

     real(ark) :: coe,emin
     real(ark) :: icoe

     real(ark) :: De1,De_1,De2,De_2,De3,De_3
     real(ark) :: edp1,edp2,edp3,edp4,edp5,edp6
     real(ark) :: rnnref,rnoref,alpha1,alpha1b,alpha2,alpha2b
     !
     integer(ik) :: i,k,Ncoe
     real(ark)   :: r12ref,Ae1,Ae2,edamp2,edamp4,edamp5,edamp6,alphae
     !
     real(ark)   :: rmin,rminbohr,alpha,rref,a2b
     !
     real(ark)   :: str1,str2,enetmp1,atp1,enetmp2,edamp11,edamp12,edamp1,V
     real(ark)   :: v0,rrco1,rrco2,r1,r2,a3,xx1,xx2,rco1,rco2,ang,dstr1,dstr2,sumstr2,sumstr4,angref,angx,bdamp2,bdamp4
     !
     real(ark)   :: v1(3),v2(3),x1,x2,cosalpha1,rnn,rno,ang2,sumx0,rno2,enetmp2B
     real(ark)   :: etmp2,anno,ax,enetmp3,etmp1
        !
        !
        rnn  = local(1) ;  rno    = local(2) ;  alpha = local(3)
        !
        alphae  = molec%alphaeq(1)
        !
        if (molec%Nangles==0) then
          stop 'MLpoten_xyz_N2O_Zobov: ilegal number of bond anlges 0'
        endif 
	    !
        Ncoe    = int(force(1))
        !
        rnnref  = force( 2)
        rnoref  = force( 3)
        alphae  = force( 4)
        alpha1  = force( 5)
        alpha2  = force( 6)
        De1     = force( 7) 
        De_1    = force( 8) 
        De2     = force( 9) 
        De_2    = force(10) 
        De3     = force(11) 
        De_3    = force(12) 
        edamp2  = force(13)
        edamp4  = force(14)
        edamp5  = force(15)
        edamp6  = force(16)
        edp1    = force(17)
        edp2    = force(18)
        edp3    = force(19)
        edp4    = force(20)
        edp5    = force(21)
        edp6    = force(22)
        !
        Emin    = force(23)
        !
        r1=rnn-rnnref
        r2=rno-rnoref
        !
        a3=sin(alpha)
        !
	    a3=a3*a3
        !
        v0=0
        do i=24,Ncoe
          v0=v0+force(i)*r1**molec%pot_ind(1,i)*r2**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
        end do
        !
        sumstr2=(rnn-rnnref)**2+(rno-rnoref)**2
        sumstr4=(rnn-rnnref)**4+(rno-rnoref)**4
        ang2=pi-alpha
        angref=150.0_ark*pi/180.0_ark-pi
        !angx=min(-abs(alpha*180.0_ark/pi)-angref,0.0_ark)
        !bdamp2=angx**2;bdamp4=angx**4        !
        bdamp2=ang2**2; bdamp4=ang2**4
        etmp2=exp(edp1*sumstr2+edp2*sumstr4+edp3*bdamp2+edp4*bdamp4)

	    sumx0=V0
	    !
        V0=V0*etmp2
        !
! calc etmp1
        sumstr2=(rnn-rnnref)**2+(rno-rnoref)**2
        sumstr4=(rnn-rnnref)**4+(rno-rnoref)**4
        enetmp1=De1*(1.0_ark-exp(-alpha1*(rnn-rnnref)))**2 + De_1*(1.0_ark-exp(-alpha1*(rnn-rnnref)))**4
        enetmp2=De2*(1.0_ark-exp(-alpha2*(rno-rnoref)))**2 + De_2*(1.0_ark-exp(-alpha2*(rno-rnoref)))**4
        !
        rno2=rnn*rnn+rno*rno-2.0_ark*rnn*rno*cos(alpha)
        rno2=sqrt(rno2)
        enetmp2B=De2*(1-exp(-alpha2*(rno2-rnoref)))**2 + De_2*(1.0_ark-exp(-alpha2*(rno2-rnoref)))**4
        anno=alpha; ax=(pi-anno)/2.0_ark
        enetmp3=De3*sin(ax)**2 + De_3*sin(ax)**4  !bending simulation
        edamp1=exp(edp5*sumstr2+0.0_ark*edp6*sumstr4)
        enetmp3=edamp1*enetmp3
        !
        etmp1=enetmp1+enetmp2+enetmp3 !+ enetmp2B
        !
        V0=V0+etmp1
	    !
        f=V0-emin*219474.63067d0
        !
      end function MLpoten_n2opotlongrange




 !returns electric dipole moment cartesian coordinates in the user-defined frame for locals specified
 !
 recursive subroutine MLdms_xyz_Kolya(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: k
    real(ark)             :: y1,y2,y3, mu(3),ux(3),uy(3),uz(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,re1,re2,ae,a0,b0
    real(ark)             :: p(3:extF%nterms(1)),q(5:extF%nterms(2)),v0,v1,v2,v3,v4,v5,v6,v7,v8,xyz0(natoms,3)
    !
    ! xyz are undefined for the local case
    !
    !write(out,"('MLdms2pqr_xyz_coeff is temporally deactivated as a bisector frame, use DIPOLE_PQR_XYZ_Z-FRAME instead')")
    !stop 'MLdms2pqr_xyz_coeff is temporally deactivated as a bisector frame, use DIPOLE_PQR_XYZ_Z-FRAME instead'
    !
    if (all(abs(xyz)<small_)) then 
      !
      xyz0 = MLloc2pqr_xyz(local)
      !
      x(1,:) = xyz0(2,:) - xyz0(1,:)
      x(2,:) = xyz0(3,:) - xyz0(1,:)
      !
    else
      !
      x(1,:) = xyz(2,:) - xyz(1,:)
      x(2,:) = xyz(3,:) - xyz(1,:)
      !
    endif
    !
    r1 = sqrt(sum(x(1,:)**2))
    r2 = sqrt(sum(x(2,:)**2))
    !
    n1 = x(1,:)/r1
    n2 = x(2,:)/r2
    !
    alpha = acos(sum(n1*n2))
    !
    ux = n1 + n2
    if (sum(ux(:)**2)<sqrt(small_)) ux = (/-1.0_ark,0.0_ark,0.0_ark/)
    !
    ux = ux / sqrt(sum(ux(:)**2))
    !
    uy = MLvector_product(n1,ux)
    !
    uy = uy / sqrt(sum(uy(:)**2))
    !
    uz = MLvector_product(ux,uy)
    !
    tmat(1, :) = ux
    tmat(2, :) = uy
    tmat(3, :) = uz
    !
    re1 = extF%coef(1,1)
    re2 = extF%coef(2,1)
    ae = extF%coef(3,1)*pi/180.0_ark
    a0 = extF%coef(4,1)
    b0 = extF%coef(5,1)
    !
    y1 = (r1 - re1)*exp(-a0*(r1-re1)**2)
    y2 = (r2 - re2)*exp(-b0*(r1-re2)**2)
    y3 = cos(alpha) - cos(ae)
    !
    k = extF%nterms(1)
    !
    p(7:k) = extF%coef(7:k,1)
    !
 v0 = p(  7)*y1**0*y2**0*y3**0
 v1 = p(  8)*y1**0*y2**0*y3**1& 
    + p(  9)*y1**1*y2**0*y3**0& 
    + p( 10)*y1**0*y2**1*y3**0
 v2 = p( 11)*y1**0*y2**0*y3**2& 
    + p( 12)*y1**1*y2**0*y3**1& 
    + p( 13)*y1**0*y2**1*y3**1& 
    + p( 14)*y1**1*y2**1*y3**0& 
    + p( 15)*y1**2*y2**0*y3**0& 
    + p( 16)*y1**0*y2**2*y3**0
 v3 = p( 17)*y1**0*y2**0*y3**3& 
    + p( 18)*y1**1*y2**0*y3**2& 
    + p( 19)*y1**0*y2**1*y3**2& 
    + p( 20)*y1**1*y2**1*y3**1& 
    + p( 21)*y1**2*y2**0*y3**1& 
    + p( 22)*y1**0*y2**2*y3**1& 
    + p( 23)*y1**2*y2**1*y3**0& 
    + p( 24)*y1**1*y2**2*y3**0& 
    + p( 25)*y1**3*y2**0*y3**0& 
    + p( 26)*y1**0*y2**3*y3**0
  v4 =p( 27)*y1**0*y2**0*y3**4& 
    + p( 28)*y1**1*y2**0*y3**3& 
    + p( 29)*y1**0*y2**1*y3**3& 
    + p( 30)*y1**1*y2**1*y3**2& 
    + p( 31)*y1**2*y2**0*y3**2& 
    + p( 32)*y1**0*y2**2*y3**2& 
    + p( 33)*y1**2*y2**1*y3**1& 
    + p( 34)*y1**1*y2**2*y3**1& 
    + p( 35)*y1**2*y2**2*y3**0& 
    + p( 36)*y1**3*y2**0*y3**1& 
    + p( 37)*y1**0*y2**3*y3**1& 
    + p( 38)*y1**3*y2**1*y3**0& 
    + p( 39)*y1**1*y2**3*y3**0& 
    + p( 40)*y1**4*y2**0*y3**0& 
    + p( 41)*y1**0*y2**4*y3**0
  v5 =p( 42)*y1**0*y2**0*y3**5& 
    + p( 43)*y1**1*y2**0*y3**4& 
    + p( 44)*y1**0*y2**1*y3**4& 
    + p( 45)*y1**1*y2**1*y3**3& 
    + p( 46)*y1**2*y2**0*y3**3& 
    + p( 47)*y1**0*y2**2*y3**3& 
    + p( 48)*y1**2*y2**1*y3**2& 
    + p( 49)*y1**1*y2**2*y3**2& 
    + p( 50)*y1**2*y2**2*y3**1& 
    + p( 51)*y1**3*y2**0*y3**2& 
    + p( 52)*y1**0*y2**3*y3**2& 
    + p( 53)*y1**3*y2**1*y3**1& 
    + p( 54)*y1**1*y2**3*y3**1& 
    + p( 55)*y1**3*y2**2*y3**0& 
    + p( 56)*y1**2*y2**3*y3**0& 
    + p( 57)*y1**4*y2**0*y3**1& 
    + p( 58)*y1**0*y2**4*y3**1& 
    + p( 59)*y1**4*y2**1*y3**0& 
    + p( 60)*y1**1*y2**4*y3**0& 
    + p( 61)*y1**5*y2**0*y3**0& 
    + p( 62)*y1**0*y2**5*y3**0
  v6 =p( 63)*y1**0*y2**0*y3**6& 
    + p( 64)*y1**1*y2**0*y3**5& 
    + p( 65)*y1**0*y2**1*y3**5& 
    + p( 66)*y1**1*y2**1*y3**4& 
    + p( 67)*y1**2*y2**0*y3**4& 
    + p( 68)*y1**0*y2**2*y3**4& 
    + p( 69)*y1**2*y2**1*y3**3& 
    + p( 70)*y1**1*y2**2*y3**3& 
    + p( 71)*y1**2*y2**2*y3**2& 
    + p( 72)*y1**3*y2**0*y3**3& 
    + p( 73)*y1**0*y2**3*y3**3& 
    + p( 74)*y1**3*y2**1*y3**2& 
    + p( 75)*y1**1*y2**3*y3**2& 
    + p( 76)*y1**3*y2**2*y3**1& 
    + p( 77)*y1**2*y2**3*y3**1& 
    + p( 78)*y1**3*y2**3*y3**0& 
    + p( 79)*y1**4*y2**0*y3**2& 
    + p( 80)*y1**0*y2**4*y3**2& 
    + p( 81)*y1**4*y2**1*y3**1& 
    + p( 82)*y1**1*y2**4*y3**1& 
    + p( 83)*y1**4*y2**2*y3**0& 
    + p( 84)*y1**2*y2**4*y3**0& 
    + p( 85)*y1**5*y2**0*y3**1& 
    + p( 86)*y1**0*y2**5*y3**1& 
    + p( 87)*y1**5*y2**1*y3**0& 
    + p( 88)*y1**1*y2**5*y3**0& 
    + p( 89)*y1**6*y2**0*y3**0& 
    + p( 90)*y1**0*y2**6*y3**0
 v7 = p( 91)*y1**0*y2**0*y3**7& 
    + p( 92)*y1**1*y2**0*y3**6& 
    + p( 93)*y1**0*y2**1*y3**6& 
    + p( 94)*y1**1*y2**1*y3**5& 
    + p( 95)*y1**2*y2**0*y3**5& 
    + p( 96)*y1**0*y2**2*y3**5& 
    + p( 97)*y1**2*y2**1*y3**4& 
    + p( 98)*y1**1*y2**2*y3**4& 
    + p( 99)*y1**2*y2**2*y3**3& 
    + p(100)*y1**3*y2**0*y3**4& 
    + p(101)*y1**0*y2**3*y3**4& 
    + p(102)*y1**3*y2**1*y3**3& 
    + p(103)*y1**1*y2**3*y3**3& 
    + p(104)*y1**3*y2**2*y3**2& 
    + p(105)*y1**2*y2**3*y3**2& 
    + p(106)*y1**3*y2**3*y3**1& 
    + p(107)*y1**4*y2**0*y3**3& 
    + p(108)*y1**0*y2**4*y3**3& 
    + p(109)*y1**4*y2**1*y3**2& 
    + p(110)*y1**1*y2**4*y3**2& 
    + p(111)*y1**4*y2**2*y3**1& 
    + p(112)*y1**2*y2**4*y3**1& 
    + p(113)*y1**4*y2**3*y3**0& 
    + p(114)*y1**3*y2**4*y3**0& 
    + p(115)*y1**5*y2**0*y3**2& 
    + p(116)*y1**0*y2**5*y3**2& 
    + p(117)*y1**5*y2**1*y3**1& 
    + p(118)*y1**1*y2**5*y3**1& 
    + p(119)*y1**5*y2**2*y3**0& 
    + p(120)*y1**2*y2**5*y3**0& 
    + p(121)*y1**6*y2**0*y3**1& 
    + p(122)*y1**0*y2**6*y3**1& 
    + p(123)*y1**6*y2**1*y3**0& 
    + p(124)*y1**1*y2**6*y3**0& 
    + p(125)*y1**7*y2**0*y3**0& 
    + p(126)*y1**0*y2**7*y3**0
 v8 = p(127)*y1**0*y2**0*y3**8& 
    + p(128)*y1**1*y2**0*y3**7& 
    + p(129)*y1**0*y2**1*y3**7& 
    + p(130)*y1**1*y2**1*y3**6& 
    + p(131)*y1**2*y2**0*y3**6& 
    + p(132)*y1**0*y2**2*y3**6& 
    + p(133)*y1**2*y2**1*y3**5& 
    + p(134)*y1**1*y2**2*y3**5& 
    + p(135)*y1**2*y2**2*y3**4& 
    + p(136)*y1**3*y2**0*y3**5& 
    + p(137)*y1**0*y2**3*y3**5& 
    + p(138)*y1**3*y2**1*y3**4& 
    + p(139)*y1**1*y2**3*y3**4& 
    + p(140)*y1**3*y2**2*y3**3& 
    + p(141)*y1**2*y2**3*y3**3& 
    + p(142)*y1**3*y2**3*y3**2& 
    + p(143)*y1**4*y2**0*y3**4& 
    + p(144)*y1**0*y2**4*y3**4& 
    + p(145)*y1**4*y2**1*y3**3& 
    + p(146)*y1**1*y2**4*y3**3& 
    + p(147)*y1**4*y2**2*y3**2& 
    + p(148)*y1**2*y2**4*y3**2& 
    + p(149)*y1**4*y2**3*y3**1& 
    + p(150)*y1**3*y2**4*y3**1& 
    + p(151)*y1**4*y2**4*y3**0& 
    + p(152)*y1**5*y2**0*y3**3& 
    + p(153)*y1**0*y2**5*y3**3& 
    + p(154)*y1**5*y2**1*y3**2& 
    + p(155)*y1**1*y2**5*y3**2& 
    + p(156)*y1**5*y2**2*y3**1& 
    + p(157)*y1**2*y2**5*y3**1& 
    + p(158)*y1**5*y2**3*y3**0& 
    + p(159)*y1**3*y2**5*y3**0& 
    + p(160)*y1**6*y2**0*y3**2& 
    + p(161)*y1**0*y2**6*y3**2& 
    + p(162)*y1**6*y2**1*y3**1& 
    + p(163)*y1**1*y2**6*y3**1& 
    + p(164)*y1**6*y2**2*y3**0& 
    + p(165)*y1**2*y2**6*y3**0& 
    + p(166)*y1**7*y2**0*y3**1& 
    + p(167)*y1**0*y2**7*y3**1& 
    + p(168)*y1**7*y2**1*y3**0& 
    + p(169)*y1**1*y2**7*y3**0& 
    + p(170)*y1**8*y2**0*y3**0& 
    + p(171)*y1**0*y2**8*y3**0
    !
    mu(3) = v0+v1+v2+v3+v4+v5+v6+v7+v8
    !
    re1 = extF%coef(1,1)
    re2 = extF%coef(2,1)
    ae = extF%coef(3,1)*pi/180.0_ark
    a0 = extF%coef(4,1)
    b0 = extF%coef(5,1)
    !
    y1 = (r1 - re1)*exp(-a0*(r1-re1)**2)
    y2 = (r2 - re2)*exp(-b0*(r1-re2)**2)
    !
    k = extF%nterms(2) 
    !
    q(5:k) = extF%coef(5:k,2)
    !
    v0 = q(  7)*y1**0*y2**0*y3**0
    v1 = q(  8)*y1**0*y2**0*y3**1& 
       + q(  9)*y1**1*y2**0*y3**0& 
       + q( 10)*y1**0*y2**1*y3**0
    v2 = q( 11)*y1**0*y2**0*y3**2& 
       + q( 12)*y1**1*y2**0*y3**1& 
       + q( 13)*y1**0*y2**1*y3**1& 
       + q( 14)*y1**1*y2**1*y3**0& 
       + q( 15)*y1**2*y2**0*y3**0& 
       + q( 16)*y1**0*y2**2*y3**0
    v3 = q( 17)*y1**0*y2**0*y3**3& 
       + q( 18)*y1**1*y2**0*y3**2& 
       + q( 19)*y1**0*y2**1*y3**2& 
       + q( 20)*y1**1*y2**1*y3**1& 
       + q( 21)*y1**2*y2**0*y3**1& 
       + q( 22)*y1**0*y2**2*y3**1& 
       + q( 23)*y1**2*y2**1*y3**0& 
       + q( 24)*y1**1*y2**2*y3**0& 
       + q( 25)*y1**3*y2**0*y3**0& 
       + q( 26)*y1**0*y2**3*y3**0
     v4 =q( 27)*y1**0*y2**0*y3**4& 
       + q( 28)*y1**1*y2**0*y3**3& 
       + q( 29)*y1**0*y2**1*y3**3& 
       + q( 30)*y1**1*y2**1*y3**2& 
       + q( 31)*y1**2*y2**0*y3**2& 
       + q( 32)*y1**0*y2**2*y3**2& 
       + q( 33)*y1**2*y2**1*y3**1& 
       + q( 34)*y1**1*y2**2*y3**1& 
       + q( 35)*y1**2*y2**2*y3**0& 
       + q( 36)*y1**3*y2**0*y3**1& 
       + q( 37)*y1**0*y2**3*y3**1& 
       + q( 38)*y1**3*y2**1*y3**0& 
       + q( 39)*y1**1*y2**3*y3**0& 
       + q( 40)*y1**4*y2**0*y3**0& 
       + q( 41)*y1**0*y2**4*y3**0
     v5 =q( 42)*y1**0*y2**0*y3**5& 
       + q( 43)*y1**1*y2**0*y3**4& 
       + q( 44)*y1**0*y2**1*y3**4& 
       + q( 45)*y1**1*y2**1*y3**3& 
       + q( 46)*y1**2*y2**0*y3**3& 
       + q( 47)*y1**0*y2**2*y3**3& 
       + q( 48)*y1**2*y2**1*y3**2& 
       + q( 49)*y1**1*y2**2*y3**2& 
       + q( 50)*y1**2*y2**2*y3**1& 
       + q( 51)*y1**3*y2**0*y3**2& 
       + q( 52)*y1**0*y2**3*y3**2& 
       + q( 53)*y1**3*y2**1*y3**1& 
       + q( 54)*y1**1*y2**3*y3**1& 
       + q( 55)*y1**3*y2**2*y3**0& 
       + q( 56)*y1**2*y2**3*y3**0& 
       + q( 57)*y1**4*y2**0*y3**1& 
       + q( 58)*y1**0*y2**4*y3**1& 
       + q( 59)*y1**4*y2**1*y3**0& 
       + q( 60)*y1**1*y2**4*y3**0& 
       + q( 61)*y1**5*y2**0*y3**0& 
       + q( 62)*y1**0*y2**5*y3**0
     v6 =q( 63)*y1**0*y2**0*y3**6& 
       + q( 64)*y1**1*y2**0*y3**5& 
       + q( 65)*y1**0*y2**1*y3**5& 
       + q( 66)*y1**1*y2**1*y3**4& 
       + q( 67)*y1**2*y2**0*y3**4& 
       + q( 68)*y1**0*y2**2*y3**4& 
       + q( 69)*y1**2*y2**1*y3**3& 
       + q( 70)*y1**1*y2**2*y3**3& 
       + q( 71)*y1**2*y2**2*y3**2& 
       + q( 72)*y1**3*y2**0*y3**3& 
       + q( 73)*y1**0*y2**3*y3**3& 
       + q( 74)*y1**3*y2**1*y3**2& 
       + q( 75)*y1**1*y2**3*y3**2& 
       + q( 76)*y1**3*y2**2*y3**1& 
       + q( 77)*y1**2*y2**3*y3**1& 
       + q( 78)*y1**3*y2**3*y3**0& 
       + q( 79)*y1**4*y2**0*y3**2& 
       + q( 80)*y1**0*y2**4*y3**2& 
       + q( 81)*y1**4*y2**1*y3**1& 
       + q( 82)*y1**1*y2**4*y3**1& 
       + q( 83)*y1**4*y2**2*y3**0& 
       + q( 84)*y1**2*y2**4*y3**0& 
       + q( 85)*y1**5*y2**0*y3**1& 
       + q( 86)*y1**0*y2**5*y3**1& 
       + q( 87)*y1**5*y2**1*y3**0& 
       + q( 88)*y1**1*y2**5*y3**0& 
       + q( 89)*y1**6*y2**0*y3**0& 
       + q( 90)*y1**0*y2**6*y3**0
    v7 = q( 91)*y1**0*y2**0*y3**7& 
       + q( 92)*y1**1*y2**0*y3**6& 
       + q( 93)*y1**0*y2**1*y3**6& 
       + q( 94)*y1**1*y2**1*y3**5& 
       + q( 95)*y1**2*y2**0*y3**5& 
       + q( 96)*y1**0*y2**2*y3**5& 
       + q( 97)*y1**2*y2**1*y3**4& 
       + q( 98)*y1**1*y2**2*y3**4& 
       + q( 99)*y1**2*y2**2*y3**3& 
       + q(100)*y1**3*y2**0*y3**4& 
       + q(101)*y1**0*y2**3*y3**4& 
       + q(102)*y1**3*y2**1*y3**3& 
       + q(103)*y1**1*y2**3*y3**3& 
       + q(104)*y1**3*y2**2*y3**2& 
       + q(105)*y1**2*y2**3*y3**2& 
       + q(106)*y1**3*y2**3*y3**1& 
       + q(107)*y1**4*y2**0*y3**3& 
       + q(108)*y1**0*y2**4*y3**3& 
       + q(109)*y1**4*y2**1*y3**2& 
       + q(110)*y1**1*y2**4*y3**2& 
       + q(111)*y1**4*y2**2*y3**1& 
       + q(112)*y1**2*y2**4*y3**1& 
       + q(113)*y1**4*y2**3*y3**0& 
       + q(114)*y1**3*y2**4*y3**0& 
       + q(115)*y1**5*y2**0*y3**2& 
       + q(116)*y1**0*y2**5*y3**2& 
       + q(117)*y1**5*y2**1*y3**1& 
       + q(118)*y1**1*y2**5*y3**1& 
       + q(119)*y1**5*y2**2*y3**0& 
       + q(120)*y1**2*y2**5*y3**0& 
       + q(121)*y1**6*y2**0*y3**1& 
       + q(122)*y1**0*y2**6*y3**1& 
       + q(123)*y1**6*y2**1*y3**0& 
       + q(124)*y1**1*y2**6*y3**0& 
       + q(125)*y1**7*y2**0*y3**0& 
       + q(126)*y1**0*y2**7*y3**0
    v8 = q(127)*y1**0*y2**0*y3**8& 
       + q(128)*y1**1*y2**0*y3**7& 
       + q(129)*y1**0*y2**1*y3**7& 
       + q(130)*y1**1*y2**1*y3**6& 
       + q(131)*y1**2*y2**0*y3**6& 
       + q(132)*y1**0*y2**2*y3**6& 
       + q(133)*y1**2*y2**1*y3**5& 
       + q(134)*y1**1*y2**2*y3**5& 
       + q(135)*y1**2*y2**2*y3**4& 
       + q(136)*y1**3*y2**0*y3**5& 
       + q(137)*y1**0*y2**3*y3**5& 
       + q(138)*y1**3*y2**1*y3**4& 
       + q(139)*y1**1*y2**3*y3**4& 
       + q(140)*y1**3*y2**2*y3**3& 
       + q(141)*y1**2*y2**3*y3**3& 
       + q(142)*y1**3*y2**3*y3**2& 
       + q(143)*y1**4*y2**0*y3**4& 
       + q(144)*y1**0*y2**4*y3**4& 
       + q(145)*y1**4*y2**1*y3**3& 
       + q(146)*y1**1*y2**4*y3**3& 
       + q(147)*y1**4*y2**2*y3**2& 
       + q(148)*y1**2*y2**4*y3**2& 
       + q(149)*y1**4*y2**3*y3**1& 
       + q(150)*y1**3*y2**4*y3**1& 
       + q(151)*y1**4*y2**4*y3**0& 
       + q(152)*y1**5*y2**0*y3**3& 
       + q(153)*y1**0*y2**5*y3**3& 
       + q(154)*y1**5*y2**1*y3**2& 
       + q(155)*y1**1*y2**5*y3**2& 
       + q(156)*y1**5*y2**2*y3**1& 
       + q(157)*y1**2*y2**5*y3**1& 
       + q(158)*y1**5*y2**3*y3**0& 
       + q(159)*y1**3*y2**5*y3**0& 
       + q(160)*y1**6*y2**0*y3**2& 
       + q(161)*y1**0*y2**6*y3**2& 
       + q(162)*y1**6*y2**1*y3**1& 
       + q(163)*y1**1*y2**6*y3**1& 
       + q(164)*y1**6*y2**2*y3**0& 
       + q(165)*y1**2*y2**6*y3**0& 
       + q(166)*y1**7*y2**0*y3**1& 
       + q(167)*y1**0*y2**7*y3**1& 
       + q(168)*y1**7*y2**1*y3**0& 
       + q(169)*y1**1*y2**7*y3**0& 
       + q(170)*y1**8*y2**0*y3**0& 
       + q(171)*y1**0*y2**8*y3**0
    !
    mu(1) =(v0+v1+v2+v3+v4+v5+v6+v7+v8)*sin(pi-alpha)
    !
    mu(2) = 0
    !
    f(1:3) = matmul(tmat,mu)
    !
 end subroutine MLdms_xyz_Kolya




end module pot_user

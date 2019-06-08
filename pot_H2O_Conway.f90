!
!  This unit is for a user defined potential 
!
module pot_user
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xy2

  implicit none

  public MLdipole,MLpoten,ML_MEP

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  !
 contains
 !
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
 !
 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)

    integer(ik)           :: k,imu,iterm
    real(ark)             :: y(3), mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,re,ae
    real(ark)             :: xyz0(natoms,3),xi(3),mu_t,cos_theta

    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      select case(trim(molec%coords_transform))
      case default
         write (out,"('MLdipole_h2o_lpt2011: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLdipole_h2o_lpt2011 - bad coord. type'
      case('R-RHO-Z')
         !
         xyz0 = MLloc2pqr_xy2(local)
         !
         x(1,:) = xyz0(2,:) - xyz0(1,:)
         x(2,:) = xyz0(3,:) - xyz0(1,:)
         !
      end select
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
    n1 = x(1,:) / r1
    n2 = x(2,:) / r2
    !
    alpha = acos(sum(n1*n2))
    !
    cos_theta = sum(n1*n2)
    !
    r1 = r1/bohr ; r2 = r2/bohr
    !
    call DIPSA(mu(2), mu(1), r1, r2, cos_theta)
    !
    mu(3) = 0
    !
    u1 = n1 + n2
    u2 = n2 - n1
    !
    u1 = u1 / sqrt(sum(u1(:)**2))
    u2 = u2 / sqrt(sum(u2(:)**2))
    !
    u3 = MLvector_product(u1,u2)
    !
    tmat(1, :) = u1
    tmat(2, :) = u2
    tmat(3, :) = u3
    !
    f(1:3) = matmul(mu,tmat)
    !
    f = f/0.393430_ark ! to debye
   !
 end subroutine MLdipole
 !
 !
 ! Defining potential energy function (built for SO2)

 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   f = MLpoten_xy2_dmbe(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !
 !

  !
  ! Defining potential energy function 
  !
  ! This is a Varandas type (DMBE) PES for XY2 molecules. It simply uses the external routine,
  ! defined in J. Chem. Phys. 118, 2637 (2003).
  !
  function MLpoten_xy2_dmbe(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)               :: r12,r32,alpha,xcos,v
   real(ark),parameter     :: tocm = 219474.63067_ark
   !
   if (verbose>=6) write(out,"('MLpoten_xy2_tennys/start')") 


     !write (out,"('MLpoten_xy2_dmbe: is turned off')")
     !stop 'MLpoten_xy2_dmbe: PES is turned off'


     r12 = local(1)/bohr ; r32 = local(2)/bohr ;  alpha = local(3)
     !
     if (molec%AtomMasses(1)>15.8_rk.and.molec%AtomMasses(1)<18.0_rk) then
       !
       ! whater
       ! 
       xcos = cos(alpha)
       !
       call potv(v,r12,r32,xcos)
       !v = 0
       v = v*tocm
       !
       !v = v + MLpoten_xy2_bubukina(ncoords,natoms,local,xyz,force)
       !
     elseif(molec%AtomMasses(1)>31.9_rk.and.molec%AtomMasses(1)<36.0_rk) then 
       !
       ! h2s
       !
       !call potv_h2s(v,r12,r32,alpha)
       v = 0 
       !
     else
        !
        write (out,"('MLpoten_xy2_dmbe: for these atoms ',3f12.6,' PES is not provided.')") molec%AtomMasses(1:3)
        stop 'MLpoten_xy2_dmbe: PES is not provided'
        !
     endif
     ! 
     f = v
     !
     if (verbose>=6) write(out,"('MLpoten_xy2_dmbe/end')") 
 
 end function MLpoten_xy2_dmbe



 recursive subroutine MLdms2pqr_dvr3d(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: i,ik(1:3)
    real(ark)             :: b(molec%ncoords), mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),alpha
    real(rk)              :: r1,r2,xcos,dipp,dipq
    !
    x(1,:) = xyz(2,:) - xyz(1,:)
    x(2,:) = xyz(3,:) - xyz(1,:)
    !
    r1 = sqrt(sum(x(1,:)**2))
    r2 = sqrt(sum(x(2,:)**2))
    !
    n1 = x(1,:) / r1
    n2 = x(2,:) / r2
    !
    alpha = acos(sum(n1*n2))
    !
    u1 = n1 + n2
    u2 = n2 - n1
    !
    u1 = u1 / sqrt(sum(u1(:)**2))
    u2 = u2 / sqrt(sum(u2(:)**2))
    !
    u3 = MLvector_product(u1,u2)
    !
    tmat(1, :) = u1
    tmat(2, :) = u2
    tmat(3, :) = u3
    !
    r1 = local(1)/bohr ; r2 = local(2)/bohr ;  alpha = local(3)
    ! 
    xcos = cos(alpha)
    !
    !call dipd(dipp,r1,r2,xcos,1)
    !
    !call dipd(dipq,r1,r2,xcos,0)
    !
    mu(1) = dipq
    mu(2) = dipp
    mu(3) = 0 
    !
    f(1:3) = matmul(mu,tmat)
    !
 end subroutine MLdms2pqr_dvr3d
 !

      SUBROUTINE potv(V,R1,R2,xcos)

      IMPLICIT real(ark) (A-H,O-Y),logical(z)
      
      !COMMON /MASS/ XMASS(3),G1,G2,xmassr(3)
      real(ark) :: RZ=0.9576257_ark,RHO=75.48992362_ark
      real(ark) :: TOANG= 0.5291772_ark
      real(ark) ::  X1= 1.0_ark,X0=0.0_ark,TINY=9.0D-15,X2=2.0_ark

      !pi=acos(-1._ark)
      G1 = 0

      IF (G1 .EQ. X0) THEN
         Q1 = R1
         Q2 = R2
         THETA = ACOS(XCOS)
      ELSE IF (G2 .EQ. X0) THEN
         XX = R1 * G1
         YY = R1 * (X1 - G1)
         IF (R2 .EQ. X0 .OR. XCOS .GE. (X1 - TINY)) THEN
            Q1 = ABS(XX - R2)
            Q2 = (YY + R2)
            COST = -X1
         ELSE IF (XCOS .LE. (TINY - X1)) THEN
            Q1 = (XX + R2)
            Q2 = ABS(YY + R2)
            COST = X1
         ELSE
            Q1 = SQRT(XX*XX + R2*R2 - X2*XX*R2*XCOS)
            Q2 = SQRT(YY*YY + R2*R2 + X2*YY*R2*XCOS)
            COST = (Q1**2 + Q2**2 - R1**2) / (X2 * Q1 * Q2)
         ENDIF
         THETA = ACOS(COST)
      ELSE
         F1= X1/G1
         F2= X1/G2
         F12= X1 - F1*F2
         P1= R1*(X1-F1)/(G2*F12)
         P2= R2*(X1-F2)/(G1*F12)
         S1= R1-P1
         S2= R2-P2
         Q1= SQRT(P1*P1 + S2*S2 + X2*P1*S2*XCOS)/(X1-G1)
         Q2= SQRT(P2*P2 + S1*S1 + X2*P2*S1*XCOS)/(X1-G2)
         Q3= SQRT(P1*P1 + P2*P2 - X2*P1*P2*XCOS)
         COST = (Q1*Q1 + Q2*Q2 - Q3*Q3)/(X2*Q1*Q2)
         THETA = ACOS(COST)
      ENDIF

      THETAeq = (180._ark - RHO)*pi/180._ark
      Req= RZ/TOANG
      Y1 = R1 - Req
      Y2 = Cos(THETA) - Cos(THETAeq)
      Y3 = R2 - Req

      S1 = (Y1 + Y3)/2.0_ark
      S2 =  Y2
      S3 = (Y1 - Y3)/2.0_ark 

      if (Q1.gt.4.0 .and. Q2.gt.4.0 .and. THETA .gt.2.25 ) then
      THETANEW  = 2.25_ark
      else
      THETANEW = THETA 
      end if

       CALL PESleq6(VQ,Q1,Q2,THETA)
       CALL BREITB3lin(Vbr,Q1,Q2,THETA)
       CALL pesd2x(Vd,Q1,Q2,THETA)
       CALL bodc16(V16,Q1,Q2,THETA)
       CALL bodc18(V18,Q1,Q2,THETA)

      xmaso16=15.990526_ark
      xmaso17=16.9947425_ark
      xmaso18=17.9947714_ark
      xmasd=2.013553214_ark
      xmash=1.00727647_ark
      v1=(xmaso18*xmaso16*(v16-v18)/(xmaso18-xmaso16))/xmaso16+&
         (2*xmash*(xmaso18*v18-xmaso16*v16)/(2*(xmaso18-xmaso16)))/xmash

       CALL rel(Vr,Q1,Q2,THETA)	
       CALL cvps(Vcvps,Q1,Q2,THETA)

      att1=Q1*0.5291772_ark
      att2=Q2*0.5291772_ark

      call poten(vp,att1,att2,THETA)
      v=vp/219474.624_ark  +Vcvps +vr+vbr+vd+VQ+V1
      !print*,att1,att2,THETA,vp,Vcvps
! This is where the ab initio potential gets defined from the bits and pieces

       !v=Vcvps+vr+vbr+vd+VQ+V1
!
! vp   is the extrapolated ICMRCI basic potential
! vr   is the relativistic (MVD1?) incremental potential
! vbr  is the Breit increment 
! vd   is the MVD2 increment
! v1   is the BODC increment
! vq   is the QED increment

      end SUBROUTINE potv
      
      SUBROUTINE BODC16(V,R1,R2,THETA)
      IMPLICIT real(ark) (A-H,O-Z)
      DIMENSION FT(69)

      real(ark),parameter :: C0 = 0.0_ark,SCALE = 1.0D-6
      integer(ik),parameter :: nv = 69
      real(ark),parameter :: ZERO = 0.0_ark
      real(ark),parameter :: RZ = 0.9576257_ark,RHO=75.48992362_ark
      real(ark),parameter :: TOANG = 0.5291772
      real(ark),parameter :: cv(69)  = (/2785.4109260806,-90.0463360612,48.5226385100,&
      201.6994168988,-19.4625405248,217.8639531832,&
      -48.9272565544,-76.0803927294,25.9675092231,-65.2206704082,&
      60.7284555026,-386.5173198637,-2.0352745279,3.5405842322,&
      -32.1830475213,126.8885218953,60.2606780157,134.6300248947,&
      -14.0991224631,560.8792812623,27.8464660146,-67.1464124186,&
      -181.5695164733,43.7287769150,-150.7468776877,-103.2652779984,&
      -159.6432245839,-24.8956181482,185.4825335893,-231.4775497546,&
      -51.9208759775,3.7066140140,-212.4422129941,207.1884490703,&
      383.1659239100,50.8641728660,35.1754939408,127.2280510484,&
      -154.4544312699,55.5787967758,282.6216945851,116.5606405651,&
      5.4433254144,-107.1461094167,-173.8895717556,-26.5859674990,&
      -560.8756697840,-237.7109212157,143.9462552048,-592.3478209334,&
      0.0,-198.6468835005,-19.9674473372,-14.1731270087,193.4510720304,&
      4.6347021028,32.9502486772,-221.1685318724,26.4090449111,&
      -268.432837385,-147.1422366151,133.5465868568,363.9554096142,&
      673.9006484856,214.9454642643,40.7735822438,65.2246188257,&
      173.0708970426,1.9795259929/)
      integer(ik) :: i5,i

      THETAeq = (180._ark - RHO)*3.141592654/180._ark
      Req= RZ/toang
      S1 = R1 - Req
      S2 = cos(THETA) - cos(THETAeq)
      S3 = R2 - Req

      Y1 = (S1 + S3)/2.0_ark
      Y2 = S2
      Y3 = (S1 - S3)/2.0_ark

      DO 88 I5=1,NV
 88   FT(I5)=0.0_ark

      FT(1)=1.0_ark

      FT(2)=Y1
      FT(3)=Y2

      FT(4)=Y1**2
      FT(5)=Y2**2
      FT(6)=Y3**2
      FT(7)=Y1*Y2

      FT(8)=Y1**3
      FT(9)=Y2**3
      FT(10)=Y1**2*Y2
      FT(11)=Y2**2*Y1
      FT(12)=Y3**2*Y1
      FT(13)=Y3**2*Y2

      FT(14)=Y1**4
      FT(15)=Y2**4
      FT(16)=Y3**4
      FT(17)=Y1**3*Y2
      FT(18)=Y2**3*Y1
      FT(19)=Y1**2*Y2**2
      FT(20)=Y1**2*Y3**2
      FT(21)=Y2**2*Y3**2
      FT(22)=Y3**2*Y1*Y2

      FT(23)=Y1**5
      FT(24)=Y2**5
      FT(25)=Y1**4*Y2
      FT(26)=Y2**4*Y1
      FT(27)=Y3**4*Y1
      FT(28)=Y3**4*Y2
      FT(29)=Y1**3*Y2**2
      FT(30)=Y1**3*Y3**2
      FT(31)=Y2**3*Y1**2
      FT(32)=Y2**3*Y3**2
      FT(33)=Y1**2*Y2*Y3**2
      FT(34)=Y1*Y2**2*Y3**2

      FT(35)=Y1**6
      FT(36)=Y2**6
      FT(37)=Y3**6
      FT(38)=Y1**5*Y2
      FT(39)=Y2**5*Y1
      FT(40)=Y1**4*Y2**2
      FT(41)=Y2**4*Y1**2
      FT(42)=Y2**4*Y3**2
      FT(43)=Y3**4*Y2**2
      FT(44)=Y1**4*Y3**2
      FT(45)=Y3**4*Y1**2
      FT(46)=Y3**4*Y1*Y2
      FT(47)=Y1**3*Y2**3
      FT(48)=Y2**3*Y1**2*Y3**2
      FT(49)=Y1**3*Y3**2*Y2
      FT(50)=Y2**3*Y3**2*Y1
      FT(51)=Y1**2*Y2**2*Y1**2


      FT(52)=Y1**7
      FT(53)=Y2**7
      FT(54)=Y1**6*Y2
      FT(55)=Y2**6*Y1
      FT(56)=Y3**6*Y1
      FT(57)=Y3**6*Y2
      FT(58)=Y1**5*Y2**2
      FT(59)=Y1**5*Y3**2
      FT(60)=Y2**5*Y1**2
      FT(61)=Y2**5*Y3**2
      FT(62)=Y1**4*Y2*Y3**2
      FT(63)=Y1**4*Y2**3
      FT(64)=Y2**4*Y1*Y3**2
      FT(65)=Y2**4*Y1**3
      FT(66)=Y3**4*Y1*Y2**2
      FT(67)=Y3**4*Y2*Y1**2
      FT(68)=Y3**4*Y1**3
      FT(69)=Y3**4*Y2**3

      V=ZERO
      DO 40 I=1,NV
   40 V=V+CV(I)*FT(I)
      V=C0+SCALE*V
      RETURN
      END SUBROUTINE BODC16

      SUBROUTINE BODC18(V,R1,R2,THETA)
      IMPLICIT real(ark) (A-H,O-Z)
      integer(ik),parameter :: nv = 117
      DIMENSION FT(nv)

      real(ark),parameter :: C0=0.0_ark,SCALE=1.0D-6
      real(ark),parameter ::  ZERO = 0.0_ark
      real(ark),parameter ::  RZ = 0.9576257,RHO = 75.48992362
      real(ark),parameter ::  TOANG =0.5291772
      real(ark),parameter ::  cv(nv) = (/2520.8349350538,-91.1209324309,49.9198880161,&
       201.5986877584,-27.0761946057,213.0632778654,-42.3461775982,&
       -89.5081758955,3.1590606555,-82.9610483797,91.2462747271,&
       -345.6845385214,2.9125633339,-9.2743637228,19.1249589270,&
       126.8439173732,102.4748093797,29.1206937761,-0.4903235370,&
       548.5647508851,57.4425043137,-65.3115364766,-32.7831814238,&
       111.9709105628,-278.4048997815,-231.4930911427,-323.4608284981,&
       -98.0602324965,181.1676886771,-282.2900662875,167.5274873935,&
       47.8266275957,-213.7402930366,20.8555467237,169.6474418937,&
       -62.1263131255,71.9040085206,72.1511749317,123.9432784777,&
       173.6777915640,58.5107154936,-27.8257872910,99.6002434065,&
       -135.2045098416,145.2475743775,-23.3633794048,-599.5375919902,&
       -540.4989859645,275.5466454040,-502.8009780288,167.2307803194,&
       -109.3097362364,-77.1423552840,492.0803282409,313.6574722992,&
       -22.4885664677,141.5175472800,-1134.9499411146,88.6368503668,&
       -587.9451428738,-247.5329792460,-4.9333140627,322.2004639448,&
       916.1768044740,479.0204491692,93.9437933859,107.5055433254,&
       0.0,63.1640462956,191.6712842623,1.6032078252,78.8696064123,&
       1.0174252672,-377.6585243760,-192.1212462644,653.6712620983,&
       0.0,261.5131272695,146.4973633837,0.0,-137.1867814770,&
       -98.7015480745,0.0,0.0,0.0,0.0,-97.1438825386,0.0,0.0,0.0,0.0,&
       0.0,22.7027546667,0.0,0.0,0.0,-341.9933632553,0.0,0.0,&
       -111.8520186501,0.0,0.0,0.0,0.0,-48.2716070600,0.0,0.0,0.0,0.0,&
       -5.0447226741,0.0,0.0,0.0,0.0,0.0,0.0,144.1528855579/)

      integer(ik) :: i5,i
       
      THETAeq = (180._ark - RHO)*3.141592654/180._ark
      Req= RZ/toang
      S1 = R1 - Req
      S2 = cos(THETA) - cos(THETAeq)
      S3 = R2 - Req

      Y1 = (S1 + S3)/2.0_ark
      Y2 = S2
      Y3 = (S1 - S3)/2.0_ark

!   Potential in Y1,Y2,Y3. Full set of
!   derivatives for each order (1-4) and upper terms for others.

      DO 88 I5=1,NV
 88   FT(I5)=0.0_ark

      FT(1)=1.0_ark

      FT(2)=Y1
      FT(3)=Y2

      FT(4)=Y1**2
      FT(5)=Y2**2
      FT(6)=Y3**2
      FT(7)=Y1*Y2

      FT(8)=Y1**3
      FT(9)=Y2**3
      FT(10)=Y1**2*Y2
      FT(11)=Y2**2*Y1
      FT(12)=Y3**2*Y1
      FT(13)=Y3**2*Y2

      FT(14)=Y1**4
      FT(15)=Y2**4
      FT(16)=Y3**4
      FT(17)=Y1**3*Y2
      FT(18)=Y2**3*Y1
      FT(19)=Y1**2*Y2**2
      FT(20)=Y1**2*Y3**2
      FT(21)=Y2**2*Y3**2
      FT(22)=Y3**2*Y1*Y2

      FT(23)=Y1**5
      FT(24)=Y2**5
      FT(25)=Y1**4*Y2
      FT(26)=Y2**4*Y1
      FT(27)=Y3**4*Y1
      FT(28)=Y3**4*Y2
      FT(29)=Y1**3*Y2**2
      FT(30)=Y1**3*Y3**2
      FT(31)=Y2**3*Y1**2
      FT(32)=Y2**3*Y3**2
      FT(33)=Y1**2*Y2*Y3**2
      FT(34)=Y1*Y2**2*Y3**2

      FT(35)=Y1**6
      FT(36)=Y2**6
      FT(37)=Y3**6
      FT(38)=Y1**5*Y2
      FT(39)=Y2**5*Y1
      FT(40)=Y1**4*Y2**2
      FT(41)=Y2**4*Y1**2
      FT(42)=Y2**4*Y3**2
      FT(43)=Y3**4*Y2**2
      FT(44)=Y1**4*Y3**2
      FT(45)=Y3**4*Y1**2
      FT(46)=Y3**4*Y1*Y2
      FT(47)=Y1**3*Y2**3
      FT(48)=Y2**3*Y1**2*Y3**2
      FT(49)=Y1**3*Y3**2*Y2
      FT(50)=Y2**3*Y3**2*Y1
      FT(51)=Y1**2*Y2**2*Y1**2


      FT(52)=Y1**7
      FT(53)=Y2**7
      FT(54)=Y1**6*Y2
      FT(55)=Y2**6*Y1
      FT(56)=Y3**6*Y1
      FT(57)=Y3**6*Y2
      FT(58)=Y1**5*Y2**2
      FT(59)=Y1**5*Y3**2
      FT(60)=Y2**5*Y1**2
      FT(61)=Y2**5*Y3**2
      FT(62)=Y1**4*Y2*Y3**2
      FT(63)=Y1**4*Y2**3
      FT(64)=Y2**4*Y1*Y3**2
      FT(65)=Y2**4*Y1**3
      FT(66)=Y3**4*Y1*Y2**2
      FT(67)=Y3**4*Y2*Y1**2
      FT(68)=Y3**4*Y1**3
      FT(69)=Y3**4*Y2**3
      FT(70)=Y1**3*Y3**2*Y2**2

      FT(71)=Y1**8
      FT(72)=Y2**8
      FT(73)=Y3**8
      FT(74)=Y1**7*Y2
      FT(75)=Y2**7*Y1
      FT(76)=Y1**6*Y2**2
      FT(77)=Y1**6*Y3**2
      FT(78)=Y2**6*Y1**2
      FT(79)=Y2**6*Y3**2
      FT(80)=Y3**6*Y1**2
      FT(81)=Y3**6*Y2**2
      FT(82)=Y3**6*Y2*Y1
      FT(83)=Y1**5*Y3**2*Y2
      FT(84)=Y1**5*Y2**3
      FT(85)=Y2**5*Y1**3
      FT(86)=Y2**5*Y3**2*Y1
      FT(87)=Y1**4*Y2**4
      FT(88)=Y1**4*Y3**4
      FT(89)=Y1**4*Y2**2*Y3**2
      FT(90)=Y2**4*Y3**4
      FT(91)=Y2**4*Y3**2*Y1**2
      FT(92)=Y3**4*Y1**2*Y2**2
      FT(93)=Y3**4*Y1*Y2**3
      FT(94)=Y3**4*Y1**3*Y2
      FT(95)=Y1**3*Y2**3*Y3**2
      
      FT(96)=Y2**9
      FT(97)=Y2**9*Y1*Y3**4
      FT(98)=Y2**9*Y1**5

      FT(99)=Y2**10
      FT(100)=Y2**10*Y1**2*Y3**2
      FT(101)=Y2**10*Y3**4
      FT(102)=Y2**10*Y1**4
      FT(103)=Y2**11
      FT(104)=Y2**11*Y1*Y3**2
      FT(105)=Y2**11*Y1**3

      FT(106)=Y2**12
      FT(107)=Y2**12*Y1**2
      FT(108)=Y2**12*Y3**2
      
      FT(109)=Y2**13
      FT(110)=Y2**13*Y1
      
      FT(111)=Y2**14

      FT(112)=Y2**8*Y1**6
      FT(113)=Y2**8*Y3**6
      FT(114)=Y2**8*Y1**4*Y3**2
      FT(115)=Y2**8*Y3**4*Y1**2
      
      FT(116)=Y2**7*Y1**7
      FT(117)=Y2**7*Y1**5*Y3**2

      V=ZERO
      DO 40 I=1,NV
   40 V=V+CV(I)*FT(I)
!     SCALE AND SHIFT THE ZERO
      V=C0+SCALE*V
      RETURN
      END SUBROUTINE BODC18

      SUBROUTINE rel(V,R1,R2,THETA)
      IMPLICIT real(ark) (A-H,O-Z)
      integer(ik),parameter :: nv = 82
      DIMENSION FT(nv)
      !
      real(ark),parameter :: ZERO = 0.0_ark
      real(ark),parameter ::  RZ = 0.9576257, RHO = 75.48992362
      real(ark),parameter ::  TOANG = 0.5291772
      real(ark),parameter ::  cv(nv) =(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.0519929240,&
      -0.0000574005,-0.0001860844,-0.0006152872,0.0000879624,&
      -0.0004932809,-0.0000565213,0.0006886853,-0.0001027470,&
      0.0003081795,0.0001078370,0.0019925277,0.0001286772,&
      -0.0004511478,0.0000072157,-0.0004439080,-0.0001416903,&
      0.0000984644,-0.0003092456,-0.0029118771,-0.0000562847,&
      0.0000377062,0.0002610717,0.0000367337,0.0,-0.0002565517,&
      0.0008262448,0.0,-0.0001036317,0.0025090209,0.0002105614,&
      -0.0000204982,0.0,-0.0000414402,-0.0000676532,0.0000918240,&
      -0.0000400594,0.0,-0.0006468368,0.0001397619,0.0005356017,&
      0.0001585601,-0.0001899817,-0.0015914516,-0.0002822918,&
      0.0,0.0004567782,0.0,0.0001064879,0.0,0.0,0.0,-0.0001250236,&
      0.0000077559,0.0007064063,0.0,0.0,0.0,0.0006102666,-0.0004987995,&
      0.0,-0.0001943536,-0.0002855510,0.0,-0.0002176976,0.0002410759,&
      -0.0001644075,-0.0001182853,0.0,0.0,-0.0000272857,0.0/)
      integer(ik) :: i5,i


      THETAeq = (180._ark - RHO)*3.141592654/180._ark
      Req= RZ/toang
      S1 = R1 - Req
      S2 = cos(THETA) - cos(THETAeq)
      S3 = R2 - Req

      Y1 = (S1 + S3)/2.0_ark
      Y2 = S2
      Y3 = (S1 - S3)/2.0_ark

      DO 88 I5=1,NV
      FT(I5)=0.0_ark
 88   continue
      FT(11)=1.0_ark
      FT(12)=Y1
      FT(13)=Y2
      FT(14)=Y1**2
      FT(15)=Y2**2
      FT(16)=Y3**2
      FT(17)=Y1*Y2
      FT(18)=Y1**3
      FT(19)=Y2**3
      FT(20)=Y1**2*Y2
      FT(21)=Y2**2*Y1
      FT(22)=Y3**2*Y1
      FT(23)=Y3**2*Y2
      FT(24)=Y1**4
      FT(25)=Y2**4
      FT(26)=Y3**4
      FT(27)=Y1**3*Y2
      FT(28)=Y2**3*Y1
      FT(29)=Y1**2*Y2**2
      FT(30)=Y1**2*Y3**2
      FT(31)=Y2**2*Y3**2
      FT(32)=Y3**2*Y1*Y2
      FT(33)=Y1**5
      FT(34)=Y2**5
      FT(35)=Y1**4*Y2
      FT(36)=Y2**4*Y1
      FT(37)=Y3**4*Y1
      FT(38)=Y3**4*Y2
      FT(39)=Y1**3*Y2**2
      FT(40)=Y1**3*Y3**2
      FT(41)=Y2**3*Y1**2
      FT(42)=Y2**3*Y3**2
      FT(43)=Y1**2*Y2*Y3**2
      FT(44)=Y1*Y2**2*Y3**2
      FT(45)=Y1**6
      FT(46)=Y2**6
      FT(47)=Y3**6
      FT(48)=Y1**5*Y2
      FT(49)=Y2**5*Y1
      FT(50)=Y1**4*Y2**2
      FT(51)=Y2**4*Y1**2
      FT(52)=Y2**4*Y3**2
      FT(53)=Y3**4*Y2**2
      FT(54)=Y1**4*Y3**2
      FT(55)=Y3**4*Y1**2
      FT(56)=Y3**4*Y1*Y2
      FT(57)=Y1**3*Y2**3
      FT(58)=Y2**3*Y1**3
      FT(59)=Y1**3*Y3**2*Y2
      FT(60)=Y2**3*Y3**2*Y1
      FT(61)=Y1**2*Y2**2*Y1**2
      FT(62)=Y1**7
      FT(63)=Y2**7
      FT(64)=Y1**6*Y2
      FT(65)=Y2**6*Y1
      FT(66)=Y3**6*Y1
      FT(67)=Y3**6*Y2
      FT(68)=Y1**5*Y2**2
      FT(69)=Y1**5*Y3**2
      FT(70)=Y2**5*Y1**2
      FT(71)=Y2**5*Y3**2
      FT(72)=Y1**4*Y2*Y3**2
      FT(73)=Y1**4*Y2**3
      FT(74)=Y2**4*Y1*Y3**2
      FT(75)=Y2**4*Y1**3
      FT(76)=Y3**4*Y1*Y2**2
      FT(77)=Y3**4*Y2*Y1**2
      FT(78)=Y3**4*Y1**3
      FT(79)=Y3**4*Y2**3
      FT(80)=Y1**3*Y3**2*Y2**2
      FT(81)=Y1**8
      FT(82)=Y2**8

      V=ZERO
      DO 40 I=12,NV
      V=V+CV(I)*FT(I)
 40   continue
      RETURN
      END SUBROUTINE rel
      
      SUBROUTINE BREITB3lin(VR,X,Y,Z)
      IMPLICIT real(ark) (A-H,O-Z)

      real(ark),parameter :: SCALE=1.0D-6
! This subroutine contains corrections to the water NBO PES due to the BREIT 
! term.
! [see for istance HM Quiney et al Chem. Phys. Lett. 290 (1998) 473  , 
! Bethe and Salpheter, "Quantum mechanics of one and two-electron atoms"]
! The corrections have been computed on a grid based on the 325 point grid 
! from P&S (see J. chem. phys., 106 (1997) 4618).
! Moreover, a few extra points have been added , as well as a cut in the radial
! coordinates (lines 1-2),  in order to account for high bending modes.
! The final grid used contains 293 points.
! Then the points have been fitted with a polynomial in X,Y and Z, using a 
! Mathematica script.
! The basis set used for the electronic calculations is the set called 'B' 
!  provided by  H.M. Quiney (see Chem. Phys. Lett. ).
!
! Those corrections have been computed by Paolo 
! email:  paolo@theory.phys.ucl.ac.uk  .

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  The inputs are in u.a., and X,Y are the distances of the H atoms from
!  the oxygen, and Z is the angle HOH in radiants. The final result is in 
!  Hartree. 

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!   Potential  symmetrical in X,Y,Z. Full set of
!   derivatives for each order (1-4) and upper terms for others

      Seq=1.8240445_ark
      e=1.809645_ark
      h=0.5_ark

      if (x.gt.3.5d0) x=3.5_ark
      if (y.gt.3.5d0) y=3.5_ark

      xm=x-y
      xp=x+y
      xp2=-e+h*xp
      zm=-Seq+z

        v1=7577.371974812879+54.16425693884002*xm**2+&
      6.774810680339735*xm**4-0.7125792752605449*xm**6-&
      212.3502737780579*xp2-&
      119.8400228527018*xm**2*xp2-29.18501086262345*xm**4*xp2+&
      221.1271759652365*xp2**2+72.97844564936518*xm**2*xp2**2+&
      89.0829000401668*xm**4*xp2**2-180.6736318109525*xp2**3+&
      321.1679099895901*xm**2*xp2**3-86.9772033725647*xm**4*xp2**3+&
      81.467133012329*xp2**4-890.357807854604*xm**2*xp2**4+&
      57.52112099193129*xp2**5+676.0925261740228*xm**2*xp2**5-&
      58.228375665496*xp2**6+1.323061103513093*zm-&
      0.5270121336782832*xm**2*zm+0.01985434592556156*xp2*&
      zm+7.715533877300036*xm**2*xp2*zm-6.229144779356836*xm**4*xp2*&
      zm-0.1710351457311522*xp2**2*zm+9.92240371274571*xm**4*xp2**2*&
      zm+10.47917821735698*xp2**3*zm-27.67291046310705*xm**2*xp2**4*&
      zm+15.42094531691969*xp2**6*zm+7.749223096520598*zm**2-&
      0.2147443261198608*xm**2*zm**2-1.164069372530965*xm**4*zm**2+&
      10.24148015693046*xp2*zm**2-2.364801830233726*xm**2*xp2*zm**2
       v2=3.405435470980123*xm**4*xp2*zm**2+11.53659470986242*xp2**2*&
      zm**2+40.71096970108562*xp2**3*zm**2-65.49114275587444*xp2**4*&
      zm**2+0.5246601257334035*zm**3+1.025008298074623*xm**2*zm**3+&
      13.57824254391274*xp2*zm**3-7.469419914001589*xp2**2*&
      zm**3-33.70757112970705*xp2**3*zm**3+30.20514216972833*xp2**4*&
      zm**3-10.53913543447923*zm**4-0.6159136295163627*xm**2*zm**4-&
      19.56431274355461*xp2*zm**4-20.81965238209867*xp2**2*zm**4+&
      5.998958874987758*xp2**3*zm**4-9.44711265431818*zm**5-&
      22.55622148750276*xp2*zm**5+16.30440168684215*xp2**2*&
      zm**5+19.20675957512514*zm**6+19.78080962673524*xp2*&
      zm**6+8.08923849773384*zm**7-10.68490632273025*zm**8

         VR=V1+V2
 
!      SCALE AND SHIFT THE ZERO
      VR=VR*scale

      RETURN
      END SUBROUTINE BREITB3lin

        SUBROUTINE PESd2x(VR,x,y,z)
        IMPLICIT real(ark) (A-H,O-Z)
   
! This subroutine contains corrections to the water NBO PES due to the Darwin 
! 2 electrons term. Those corrections have been computed by Gyorgy 
! tarczay@para.chem.elte.hu
! Those corrections are computed on the P&S grid of 325 points.
!  (see J. chem. phys., 106 (1997) 4618), to which a few x points have been 
! added in irder to account to high bending modes. The final grid containd 341
!  points.

!  The input are in u.a., and X,Y are the distances of the H atoms from
!  the oxygen, and Z is the angle HOH in radiants. The final result is in 
!  Hartree.  

      Seq=1.8240445_ark
      pi = acos(-1.0_ark)
      e=1.809645_ark
      h=0.5_ark
      x1=x
      y1=y
      if (x.gt.3.5_ark) x1=3.5_ark
      if (y.gt.3.5_ark) y1=3.5_ark

      xm=x1-y1
      xp=x1+y1
      xp2=-e+h*xp
      zm=-Seq+z

!   Potential  symmetrical in S1,S2,S3. Full set of
!   derivatives for each order (1-4) and upper terms for others
!
         v1=-3263.067522028298-9.17958228916955*xm**2-&
      0.5165611032611932*xm**4-0.5157212949525876*xm**6+&
      18.31161578215203*xp2+&
      30.14193751791963*xm**2*xp2-13.62543868575853*xm**4*xp2-&
      43.37232019119388*xp2**2-50.50353364079294*xm**2*xp2**2+&
      70.36441193443143*xm**4*xp2**2+27.43935454999898*xp2**3+&
      123.751990625258*xm**2*xp2**3-76.80240321256033*xm**4*xp2**3-&
      9.50017804016001*xp2**4-363.4487347625543*xm**2*xp2**4+&
      113.1940248029751*xp2**5+376.6560011408163*xm**2*xp2**5-&
      164.6523756673548*xp2**6+9.16256842998227*zm&
       -1.22230095639504*xm**2*zm+1.33032571356463*xp2*&
      zm-0.94822119654751*xm**2*xp2*zm+0.7645470802285307*xm**4*xp2*&
      zm-11.77270680473595*xp2**2*zm-0.4065994514809928*xm**4*xp2**2*&
      zm-2.113651214829342*xp2**3*zm-3.653921741665064*xm**2*xp2**4*&
      zm+26.53983199106825*xp2**6*zm+3.099164302936567*zm**2
      v2=-0.4668245990549825*xm**2*zm**2+0.05845413180128318*xm**4*zm**2&
      +2.708722250876111*xp2*zm**2-2.578482367020144*xm**2*xp2*zm**2-&
      0.1605392233404811*xm**4*xp2*zm**2-10.57780429022803*xp2**2*&
      zm**2-3.496293826717189*xp2**3*zm**2+23.46280699747645*xp2**4*&
      zm**2+1.8547816858377*zm**3-0.4003662844685243*xm**2*zm**3+&
      3.040229985315839*xp2*zm**3-4.955739113923876*xp2**2*&
      zm**3+14.05364889791468*xp2**3*zm**3-21.6926320924828*xp2**4*&
      zm**3-1.321464834042384*zm**4+2.298844571392118*xm**2*zm**4-&
      2.633405421645483*xp2*zm**4+20.97178840867901*xp2**2*&
      zm**4-32.18658937476802*xp2**3*zm**4-0.5992225949734171*zm**5+&
      2.059827452250273*xp2*zm**5+0.6453850286056735*xp2**2*&
      zm**5-0.4620689505336259*zm**6+0.7465042626807512*xp2*&
      zm**6-0.1254018119377959*zm**7+0.01947721364782498*zm**8

      VR=V1+V2

!     SCALE AND SHIFT THE ZERO
      VR=VR*1.0d-6
        RETURN
        END SUBROUTINE PESd2x

        SUBROUTINE PESleq6(V1,x,y,z)
        IMPLICIT real(ark) (A-H,O-Z)

        V1=-3514.850376005703_ark+1.189315215135138_ark*(x-y)**2+&
       0.824157459531989_ark*(x-y)**4+0.03853108456851828_ark*(x-y)**6+&
       12.83265590340491_ark*(-1.809659_ark+0.5_ark*(x+y))-&
       9.51736455454466_ark*(x-y)**2*(-1.809659_ark+0.5_ark*(x+y))-&
       3.027576695974858_ark*(x-y)**4*(-1.809659_ark+0.5_ark*(x+y))+&
       10.94033338777717_ark*(-1.809659_ark+0.5_ark*(x+y))**2+&
       15.53332877554612_ark*(x-y)**2*(-1.809659_ark+0.5_ark*(x+y))**2+&
       6.063907309056958_ark*(x-y)**4*(-1.809659_ark+0.5_ark*(x+y))**2-&
       13.79644533708051_ark*(-1.809659_ark+0.5_ark*(x+y))**3-&
       26.67549601926293_ark*(x-y)**2*(-1.809659_ark+0.5_ark*(x+y))**3-&
       6.200836894255189_ark*(x-y)**4*(-1.809659_ark+0.5_ark*(x+y))**3+&
       5.688103460541242_ark*(-1.809659_ark+0.5_ark*(x+y))**4+&
       52.9835500898771_ark*(x-y)**2*(-1.809659_ark+0.5_ark*(x+y))**4-&
       11.88910471926647_ark*(-1.809659_ark+0.5_ark*(x+y))**5-&
       43.99657824332825_ark*(x-y)**2*(-1.809659_ark+0.5_ark*(x+y))**5+&
       15.874512160015_ark*(-1.809659_ark+0.5*(x+y))**6-&
       8.60706112101134_ark*(-1.824045_ark+z)+&
       1.264485336667462_ark*(x-y)**2*(-1.824045_ark+z)-&
       0.915127202947929_ark*(-1.809659_ark+0.5_ark*(x+y))*(-1.824045_ark+z)
       v1=v1+  0.6556566908758441_ark*(x-y)**2*(-1.809659_ark+0.5_ark*(x+y))*&
        (-1.824045_ark+z)-&
       0.813078328219753_ark*(x-y)**4*(-1.809659_ark+0.5_ark*(x+y))*&
        (-1.824045_ark+z)+&
       12.42234678727481_ark*(-1.809659_ark+0.5_ark*(x+y))**2*&
        (-1.824045_ark+z)+&
       2.805488560712774_ark*(x-y)**4*(-1.809659_ark+0.5_ark*(x+y))**2*&
        (-1.824045_ark+z)-&
       4.937250627623143_ark*(-1.809659_ark+0.5_ark*(x+y))**3*&
        (-1.824045_ark+z)-&
       3.095201035295474_ark*(x-y)**2*(-1.809659_ark+0.5_ark*(x+y))**4*&
        (-1.824045_ark+z)+&
       18.85309150691318_ark*(-1.809659_ark+0.5_ark*(x+y))**6*&
        (-1.824045_ark+z)-&
       3.209703208057476_ark*(-1.824045_ark+z)**2+&
       0.5360421552708203_ark*(x-y)**2*(-1.824045_ark+z)**2-&
       0.263467844585989_ark*(x-y)**4*(-1.824045_ark+z)**2-&
       1.13298516075929_ark*(-1.809659_ark+0.5_ark*(x+y))*(-1.824045_ark+z)**2
       v1=v1-  0.06909229322445753_ark*(x-y)**2*(-1.809659_ark+0.5_ark*(x+y))*&
        (-1.824045_ark+z)**2+&
       1.649063526503709_ark*(x-y)**4*(-1.809659_ark+0.5_ark*(x+y))*&
        (-1.824045_ark+z)**2+&
       3.603611347474725_ark*(-1.809659_ark+0.5_ark*(x+y))**2*&
        (-1.824045_ark+z)**2+&
       3.757418764813337_ark*(-1.809659_ark+0.5_ark*(x+y))**3*&
        (-1.824045_ark+z)**2+&
       4.607672502246032_ark*(-1.809659_ark+0.5_ark*(x+y))**4*&
        (-1.824045_ark+z)**2-&
       0.7490414640610651_ark*(-1.824045_ark+z)**3-&
       0.0888181500794012_ark*(x-y)**2*(-1.824045_ark+z)**3-&
       5.334303151299991_ark*(-1.809659_ark+0.5_ark*(x+y))*&
        (-1.824045_ark+z)**3+&
       1.37948603262339_ark*(-1.809659_ark+0.5_ark*(x+y))**2*&
        (-1.824045_ark+z)**3+&
       11.24395154910416_ark*(-1.809659_ark+0.5_ark*(x+y))**3*&
        (-1.824045_ark+z)**3 -&
       17.85690001161674_ark*(-1.809659_ark+0.5_ark*(x+y))**4*(-1.824045_ark+z)**3
       v1=v1+  0.7694433624551493_ark*(-1.824045_ark+z)**4-&
       0.939662303404418_ark*(x-y)**2*(-1.824045_ark+z)**4-&
       2.296000209594694_ark*(-1.809659_ark+0.5_ark*(x+y))*&
        (-1.824045_ark+z)**4-&
       4.514249057965571_ark*(-1.809659_ark+0.5_ark*(x+y))**2*&
        (-1.824045_ark+z)**4-&
       2.324765391545952_ark*(-1.809659_ark+0.5_ark*(x+y))**3*&
        (-1.824045_ark+z)**4+&
       0.223711667169141_ark*(-1.824045_ark+z)**5+&
       1.164515013150094_ark*(-1.809659_ark+0.5_ark*(x+y))*&
        (-1.824045_ark+z)**5-&
       2.825913168656484_ark*(-1.809659_ark+0.5_ark*(x+y))**2*&
        (-1.824045_ark+z)**5+&
       0.4811142779617512_ark*(-1.824045_ark+z)**6+&
       1.292817090808966_ark*(-1.809659_ark+0.5_ark*(x+y))*&
        (-1.824045_ark+z)**6+&
       0.1657130839026308_ark*(-1.824045_ark+z)**7-&
       0.02192338698614548_ark*(-1.824045_ark+z)**8

        v1=v1/1000000.0_ark
        RETURN
        END SUBROUTINE PESleq6


      SUBROUTINE cvps(V,rij1,rij2,rij3)
      implicit real(16) (a-h,o-z)
      integer(ik) :: i,j
!
!     pes for h2o,
!     Harry Partridge and David W. Schwenke, J. Chem. Phys.,
!     submitted Nov. 8, 1996.
!     rij(i,1)& rij(i,2) are oh distances in au
!     rij(i,3) is hoh angle in rad
!     v(i) is pes in au
!     n is number of geometries
!     mass dependent factors are included. the nuclear masses
!     should be passed to this program using the array xm in
!     common potmcm. xm(1) is the
!     mass of the hydrogen associated with rij(i,1), and xm(2)
!     is the mass of the hydrogen associated with rij(i,2).
!     all masses are in au.
!
      real(ark) :: c5z(245),cbasis(245),ccore(245),&
               crest(245),idx(245,3),fmat(15,3),cmass(9),idxm(9,3),rad
!      common/potrot/fact1,fact2,c1,s1,icoord,xm(2),xmx,iperm   
! $$$      common/potmcm/xm(2)
!
!     expansion indicies
!
      integer(ik),save :: ifirst = 0 

       data (idx(i,1),i=1,245)/                                   &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, &
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
      2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
      3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, &
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
      4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, &
      4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, &
      6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
      6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, &
      5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, &
      7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, &
      6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, &
      9, 9, 9, 9, 9/
       data (idx(i,2),i=1,245)/                                  &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,&
      1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
      2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,&
      2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3,&
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
      2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,&
      3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,&
      1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3,&
      2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4,&
      4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,&
      2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,&
      4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1,&
      1, 1, 1, 1, 1/                                              
      data (idx(i,3),i=1,245)/                                   &
      1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 1, 2, 3, 4, 5,&
      6, 7, 8, 9,10,11,12,13,14, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,&
     12,13, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1, 2, 3, 4, 5,&
      6, 7, 8, 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 1,&
      2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,&
     11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8,&
      9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 7, 8,&
      9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9,&
      1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2,&
      3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6,&
      7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3,&
      4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2,&
      3, 4, 5, 6, 7/
!
!     expansion coefficients for 5z ab initio data
!
       data (c5z(i),i=1,45)/                                         &
      4.2278462684916D+04, 4.5859382909906D-02, 9.4804986183058D+03, &
      7.5485566680955D+02, 1.9865052511496D+03, 4.3768071560862D+02, &
      1.4466054104131D+03, 1.3591924557890D+02,-1.4299027252645D+03, &
      6.6966329416373D+02, 3.8065088734195D+03,-5.0582552618154D+02, &
     -3.2067534385604D+03, 6.9673382568135D+02, 1.6789085874578D+03, &
     -3.5387509130093D+03,-1.2902326455736D+04,-6.4271125232353D+03, &
     -6.9346876863641D+03,-4.9765266152649D+02,-3.4380943579627D+03, &
      3.9925274973255D+03,-1.2703668547457D+04,-1.5831591056092D+04, &
      2.9431777405339D+04, 2.5071411925779D+04,-4.8518811956397D+04, &
     -1.4430705306580D+04, 2.5844109323395D+04,-2.3371683301770D+03, &
      1.2333872678202D+04, 6.6525207018832D+03,-2.0884209672231D+03, &
     -6.3008463062877D+03, 4.2548148298119D+04, 2.1561445953347D+04, &
     -1.5517277060400D+05, 2.9277086555691D+04, 2.6154026873478D+05, &
     -1.3093666159230D+05,-1.6260425387088D+05, 1.2311652217133D+05, &
     -5.1764697159603D+04, 2.5287599662992D+03, 3.0114701659513D+04/ 
                                                                      
      data (c5z(i),i=46,90)/                                         &
     -2.0580084492150D+03, 3.3617940269402D+04, 1.3503379582016D+04, &
     -1.0401149481887D+05,-6.3248258344140D+04, 2.4576697811922D+05, &
      8.9685253338525D+04,-2.3910076031416D+05,-6.5265145723160D+04, &
      8.9184290973880D+04,-8.0850272976101D+03,-3.1054961140464D+04, &
     -1.3684354599285D+04, 9.3754012976495D+03,-7.4676475789329D+04, &
     -1.8122270942076D+05, 2.6987309391410D+05, 4.0582251904706D+05, &
     -4.7103517814752D+05,-3.6115503974010D+05, 3.2284775325099D+05, &
      1.3264691929787D+04, 1.8025253924335D+05,-1.2235925565102D+04, &
     -9.1363898120735D+03,-4.1294242946858D+04,-3.4995730900098D+04, &
      3.1769893347165D+05, 2.8395605362570D+05,-1.0784536354219D+06, &
     -5.9451106980882D+05, 1.5215430060937D+06, 4.5943167339298D+05, &
     -7.9957883936866D+05,-9.2432840622294D+04, 5.5825423140341D+03, &
      3.0673594098716D+03, 8.7439532014842D+04, 1.9113438435651D+05, &
     -3.4306742659939D+05,-3.0711488132651D+05, 6.2118702580693D+05, &
     -1.5805976377422D+04,-4.2038045404190D+05, 3.4847108834282D+05/ 
                                                                      
      data (c5z(i),i=91,135)/                                        &
     -1.3486811106770D+04, 3.1256632170871D+04, 5.3344700235019D+03, &
      2.6384242145376D+04, 1.2917121516510D+05,-1.3160848301195D+05, &
     -4.5853998051192D+05, 3.5760105069089D+05, 6.4570143281747D+05, &
     -3.6980075904167D+05,-3.2941029518332D+05,-3.5042507366553D+05, &
      2.1513919629391D+03, 6.3403845616538D+04, 6.2152822008047D+04, &
     -4.8805335375295D+05,-6.3261951398766D+05, 1.8433340786742D+06, &
      1.4650263449690D+06,-2.9204939728308D+06,-1.1011338105757D+06, &
      1.7270664922758D+06, 3.4925947462024D+05,-1.9526251371308D+04, &
     -3.2271030511683D+04,-3.7601575719875D+05, 1.8295007005531D+05, &
      1.5005699079799D+06,-1.2350076538617D+06,-1.8221938812193D+06, &
      1.5438780841786D+06,-3.2729150692367D+03, 1.0546285883943D+04, &
     -4.7118461673723D+04,-1.1458551385925D+05, 2.7704588008958D+05, &
      7.4145816862032D+05,-6.6864945408289D+05,-1.6992324545166D+06, &
      6.7487333473248D+05, 1.4361670430046D+06,-2.0837555267331D+05, &
      4.7678355561019D+05,-1.5194821786066D+04,-1.1987249931134D+05/ 
                                                                     
                                                                      
      data (c5z(i),i=136,180)/                                       &
      1.3007675671713D+05, 9.6641544907323D+05,-5.3379849922258D+05, &
     -2.4303858824867D+06, 1.5261649025605D+06, 2.0186755858342D+06, &
     -1.6429544469130D+06,-1.7921520714752D+04, 1.4125624734639D+04, &
     -2.5345006031695D+04, 1.7853375909076D+05,-5.4318156343922D+04, &
     -3.6889685715963D+05, 4.2449670705837D+05, 3.5020329799394D+05, &
      9.3825886484788D+03,-8.0012127425648D+05, 9.8554789856472D+04, &
      4.9210554266522D+05,-6.4038493953446D+05,-2.8398085766046D+06, &
      2.1390360019254D+06, 6.3452935017176D+06,-2.3677386290925D+06, &
     -3.9697874352050D+06,-1.9490691547041D+04, 4.4213579019433D+04, &
      1.6113884156437D+05,-7.1247665213713D+05,-1.1808376404616D+06, &
      3.0815171952564D+06, 1.3519809705593D+06,-3.4457898745450D+06, &
      2.0705775494050D+05,-4.3778169926622D+05, 8.7041260169714D+03, &
      1.8982512628535D+05,-2.9708215504578D+05,-8.8213012222074D+05, &
      8.6031109049755D+05, 1.0968800857081D+06,-1.0114716732602D+06, &
      1.9367263614108D+05, 2.8678295007137D+05,-9.4347729862989D+04/ 
                                                                      
      data (c5z(i),i=181,225)/                                       &
      4.4154039394108D+04, 5.3686756196439D+05, 1.7254041770855D+05, &
     -2.5310674462399D+06,-2.0381171865455D+06, 3.3780796258176D+06, &
      7.8836220768478D+05,-1.5307728782887D+05,-3.7573362053757D+05, &
      1.0124501604626D+06, 2.0929686545723D+06,-5.7305706586465D+06, &
     -2.6200352535413D+06, 7.1543745536691D+06,-1.9733601879064D+04, &
      8.5273008477607D+04, 6.1062454495045D+04,-2.2642508675984D+05, &
      2.4581653864150D+05,-9.0376851105383D+05,-4.4367930945690D+05, &
      1.5740351463593D+06, 2.4563041445249D+05,-3.4697646046367D+03, &
     -2.1391370322552D+05, 4.2358948404842D+05, 5.6270081955003D+05, &
     -8.5007851251980D+05,-6.1182429537130D+05, 5.6690751824341D+05, &
     -3.5617502919487D+05,-8.1875263381402D+02,-2.4506258140060D+05, &
      2.5830513731509D+05, 6.0646114465433D+05,-6.9676584616955D+05, &
      5.1937406389690D+05, 1.7261913546007D+05,-1.7405787307472D+04, &
     -3.8301842660567D+05, 5.4227693205154D+05, 2.5442083515211D+06, &
     -1.1837755702370D+06,-1.9381959088092D+06,-4.0642141553575D+05/ 
                                                                     
                                                                      
      data (c5z(i),i=226,245)/                                       &
      1.1840693827934D+04,-1.5334500255967D+05, 4.9098619510989D+05, &
      6.1688992640977D+05, 2.2351144690009D+05,-1.8550462739570D+06, &
      9.6815110649918D+03,-8.1526584681055D+04,-8.0810433155289D+04, &
      3.4520506615177D+05, 2.5509863381419D+05,-1.3331224992157D+05, &
     -4.3119301071653D+05,-5.9818343115856D+04, 1.7863692414573D+03, &
      8.9440694919836D+04,-2.5558967650731D+05,-2.2130423988459D+04, &
      4.4973674518316D+05,-2.2094939343618D+05/                      
!                                                                     
!    expansion coefficients for basis correction                     
!                                                                     
      data (cbasis(i),i=1,45)/                                       &
      6.9770019624764D-04,-2.4209870001642D+01, 1.8113927151562D+01, &
      3.5107416275981D+01,-5.4600021126735D+00,-4.8731149608386D+01, &
      3.6007189184766D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -7.7178474355102D+01,-3.8460795013977D+01,-4.6622480912340D+01, &
      5.5684951167513D+01, 1.2274939911242D+02,-1.4325154752086D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-6.0800589055949D+00, &
      8.6171499453475D+01,-8.4066835441327D+01,-5.8228085624620D+01, &
      2.0237393793875D+02, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      3.3525582670313D+02, 7.0056962392208D+01,-4.5312502936708D+01/ 
                                                                     
                                                                      
      data (cbasis(i),i=46,90)/                                      &
     -3.0441141194247D+02, 2.8111438108965D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.2983583774779D+02, 3.9781671212935D+01, &
     -6.6793945229609D+01,-1.9259805675433D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-8.2855757669957D+02,-5.7003072730941D+01, &
     -3.5604806670066D+01, 9.6277766002709D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 8.8645622149112D+02,-7.6908409772041D+01, &
      6.8111763314154D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                     
                                                                      
      data (cbasis(i),i=91,135)/                                     &
      2.5090493428062D+02,-2.3622141780572D+02, 5.8155647658455D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 2.8919570295095D+03, &
     -1.7871014635921D+02,-1.3515667622500D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-3.6965613754734D+03, 2.1148158286617D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-1.4795670139431D+03, &
      3.6210798138768D+02, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -5.3552886800881D+03, 3.1006384016202D+02, 0.0000000000000D+00/ 
                                                                     
                                                                     
                                                                     
                                                                      
      data (cbasis(i),i=136,180)/                                    &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 1.6241824368764D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 4.3764909606382D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 1.0940849243716D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 3.0743267832931D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (cbasis(i),i=181,225)/                                    &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (cbasis(i),i=226,245)/                                    &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00/                      
!                                                                     
!    expansion coefficients for core correction                      
!                                                                     
      data (ccore(i),i=1,45)/                                        &
      2.4332191647159D-02,-2.9749090113656D+01, 1.8638980892831D+01, &
     -6.1272361746520D+00, 2.1567487597605D+00,-1.5552044084945D+01, &
      8.9752150543954D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -3.5693557878741D+02,-3.0398393196894D+00,-6.5936553294576D+00, &
      1.6056619388911D+01, 7.8061422868204D+01,-8.6270891686359D+01, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-3.1688002530217D+01, &
      3.7586725583944D+01,-3.2725765966657D+01,-5.6458213299259D+00, &
      2.1502613314595D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      5.2789943583277D+02,-4.2461079404962D+00,-2.4937638543122D+01/ 
                                                                     
                                                                      
      data (ccore(i),i=46,90)/                                       &
     -1.1963809321312D+02, 2.0240663228078D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-6.2574211352272D+02,-6.9617539465382D+00, &
     -5.9440243471241D+01, 1.4944220180218D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.2851139918332D+03,-6.5043516710835D+00, &
      4.0410829440249D+01,-6.7162452402027D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 1.0031942127832D+03, 7.6137226541944D+01, &
     -2.7279242226902D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (ccore(i),i=91,135)/                                      &
     -3.3059000871075D+01, 2.4384498749480D+01,-1.4597931874215D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 1.6559579606045D+03, &
      1.5038996611400D+02,-7.3865347730818D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.9738401290808D+03,-1.4149993809415D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-1.2756627454888D+02, &
      4.1487702227579D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -1.7406770966429D+03,-9.3812204399266D+01, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (ccore(i),i=136,180)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.1890301282216D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 2.3723447727360D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.0279968223292D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 5.7153838472603D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                      
      data (ccore(i),i=181,225)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                      
      data (ccore(i),i=226,245)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00/                      
!                                                                     
!    expansion coefficients for v rest                               
!                                                                     
      data (crest(i),i=1,45)/                                        &
      0.0000000000000D+00,-4.7430930170000D+00,-1.4422132560000D+01, &
     -1.8061146510000D+01, 7.5186735000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -2.7962099800000D+02, 1.7616414260000D+01,-9.9741392630000D+01, &
      7.1402447000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-7.8571336480000D+01, &
      5.2434353250000D+01, 7.7696745000000D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      1.7799123760000D+02, 1.4564532380000D+02, 2.2347226000000D+02/ 
                                                                      
      data (crest(i),i=46,90)/                                       &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-4.3823284100000D+02,-7.2846553000000D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-2.6752313750000D+02, 3.6170310000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                      
      data (crest(i),i=91,135)/                                      &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (crest(i),i=136,180)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (crest(i),i=181,225)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (crest(i),i=226,245)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00/
!
!     expansion indicies for mass correction
!
       data idxm/1,2,1,1,3,2,1,2,1,&
                2,1,1,3,1,2,2,1,1, &
                1,1,2,1,1,1,2,2,3/
!
!     expansion coefficients for mass correction
!
       data cmass/ -8.3554183D+00,3.7036552D+01,-5.2722136D+00,&
           1.6843857D+01,-7.0929741D+01,5.5380337D+00,-2.9962997D+01,&
           1.3637682D+02,-3.0530195d+00/
!
!     two body parameters
!
       data reoh,thetae,b1,roh,alphaoh,deoh,phh1,phh2/0.958649_ark,&
      104.3475_ark,2.0_ark,0.9519607159623009_ark,2.587949757553683_ark,&
      42290.92019288289_ark,16.94879431193463_ark,12.66426998162947_ark/
!
!     scaling factors for contributions to emperical potential
!
        data f5z,fbasis,fcore,frest/0.0_ark,0.0_ark,-1.0_ark,0.0_ark/
      if(ifirst.eq.0)then
       ifirst=1
    1  format(/1x,'pes for h2o',&
             /1x,'by Harry Partridge and David W. Schwenke',&
             /1x,'submitted to J. Chem. Phys. Nov. 8, 1996')
   56  format(/1x,'parameters before adjustment')
   55  format(/1x,'two body potential parameters:',&
             /1x,'hh: phh1 = ',f10.1,' phh2 = ',f5.2,&
             /1x,'oh: deoh = ',f10.1,' alpha = ',f7.4,&
             ' re = ',f7.4)
    4  format(/1x,'three body parameters:',&
             /1x,'reoh = ',f10.4,' thetae = ',f10.4,&
             /1x,'betaoh = ',f10.4,&
             /1x,'    i    j    k',7x,'c5z',9x,'cbasis',10x,'ccore',&
             10x,'crest')
       do 2 i=1,245
    2  continue
!
!     remove mass correction from vrest
!
       xmh=1836.152697_ark
       xmhi=1.0_ark/xmh
       xmd=3670.483031_ark
       fact=1._ark/((1._ark/xmd)-(1._ark/xmh))
   65  format(/1x,'parameters for delta v hdo ',&
            /1x,'    i    j    k')
       do 60 i=1,9
        cmass(i)=cmass(i)*fact
        corr=cmass(i)*xmhi
        if(idxm(i,1).eq.idxm(i,2))corr=corr*0.5_ark
        do 61 j=1,245
         if(idx(j,1).eq.idxm(i,1).and.idx(j,2).eq.idxm(i,2).and.&
           idx(j,3).eq.idxm(i,3))then
          crest(j)=crest(j)-corr
          go to 62
         end if
   61   continue
   62   continue
        do 63 j=1,245
         if(idx(j,2).eq.idxm(i,1).and.idx(j,1).eq.idxm(i,2).and.&
           idx(j,3).eq.idxm(i,3))then
          crest(j)=crest(j)-corr
          go to 64
         end if
   63   continue
   64   continue
   60  continue
       xm1=1._ark/1836.152697_ark
       xm2=1._ark/1836.152697_ark
!
!     adjust parameters using scale factors
!
   57  format(/1x,'adjusting parameters using scale factors ',&
             /1x,'f5z =    ',f11.8,                           &
             /1x,'fbasis = ',f11.8,                           &
             /1x,'fcore =  ',f11.8,                           &
             /1x,'frest =  ',f11.8)
       phh1=phh1*f5z
       deoh=deoh*f5z
       do 59 i=1,245
        c5z(i)=f5z*c5z(i)+fbasis*cbasis(i)+fcore*ccore(i)&
            +frest*crest(i)
   59  continue
   58  format(/1x,'three body parameters:',         &
             /1x,'reoh = ',f10.4,' thetae = ',f10.4,&
             /1x,'betaoh = ',f10.4,                 &
             /1x,'    i    j    k   cijk',          &
             /(1x,3i5,1pe15.7))
       do 66 i=1,9
        cmass(i)=cmass(i)*frest
   66  continue
   76  format(/1x,'mass correction factors ',&
             /1x,'    i    j    k   cijk',   &
             /(1x,3i5,1pe15.7))
!
!     convert parameters from 1/cm, angstrom to a.u.
!
       reoh=reoh/0.529177249_ark
       b1=b1*0.529177249_ark*0.529177249_ark
       do 3 i=1,245
        c5z(i)=c5z(i)*4.556335e-6_ark
    3  continue
       do 67 i=1,9
        cmass(i)=cmass(i)*4.556335e-6_ark
   67  continue
       rad=acos(-1.0_ark)/180.0_ark
       ce=cos(thetae*rad)
       phh1=phh1*exp(phh2)
       phh1=phh1*4.556335e-6_ark
       phh2=phh2*0.529177249_ark
       deoh=deoh*4.556335e-6_ark
       roh=roh/0.529177249_ark
       alphaoh=alphaoh*0.529177249_ark
       c5z(1)=c5z(1)*2.0_ark
      end if
       x1=(rij1-reoh)/reoh
       x2=(rij2-reoh)/reoh
       x3=cos(rij3)-ce
       rhh=sqrt(rij1**2+rij2**2 &
           -2._ark*rij1*rij2*cos(rij3))
       vhh=phh1*exp(-phh2*rhh)
       ex=exp(-alphaoh*(rij1-roh))
       voh1=deoh*ex*(ex-2._ark)
       ex=exp(-alphaoh*(rij2-roh))
       voh2=deoh*ex*(ex-2._ark)
       fmat(1,1)=1._ark
       fmat(1,2)=1._ark
       fmat(1,3)=1._ark
       do 10 j=2,15
        fmat(j,1)=fmat(j-1,1)*x1
        fmat(j,2)=fmat(j-1,2)*x2
        fmat(j,3)=fmat(j-1,3)*x3
   10  continue
       v=0
       do 12 j=2,245
        term=c5z(j)*(fmat(idx(j,1),1)*fmat(idx(j,2),2)&
                         +fmat(idx(j,2),1)*fmat(idx(j,1),2))&
                         *fmat(idx(j,3),3)
        v=v+term
   12  continue
       v1=0
       v2=0
       do 13 j=1,9
        v1=v1+cmass(j)*fmat(idxm(j,1),1)*fmat(idxm(j,2),2)&
            *fmat(idxm(j,3),3)
        v2=v2+cmass(j)*fmat(idxm(j,2),1)*fmat(idxm(j,1),2)&
            *fmat(idxm(j,3),3)
   13  continue
       v=v+xm1*v1+xm2*v2

       v=v*exp(-b1*((rij1-reoh)**2+(rij2-reoh)**2))&
            +c5z(1)&
            +voh1+voh2+vhh

      return
      end SUBROUTINE cvps



      SUBROUTINE cvps_v1(V,rij1,rij2,rij3)
      implicit real(ark) (a-h,o-z)

!
!     pes for h2o,
!     Harry Partridge and David W. Schwenke, J. Chem. Phys.,
!     submitted Nov. 8, 1996.
!     rij(i,1)& rij(i,2) are oh distances in au
!     rij(i,3) is hoh angle in rad
!     v(i) is pes in au
!     n is number of geometries
!     mass dependent factors are included. the nuclear masses
!     should be passed to this program using the array xm in
!     common potmcm. xm(1) is the
!     mass of the hydrogen associated with rij(i,1), and xm(2)
!     is the mass of the hydrogen associated with rij(i,2).
!     all masses are in au.
!
      real(ark) c5z(245),cbasis(245),ccore(245),&
               crest(245),fmat(15,3),cmass(9),rad
      integer(ik) :: idx(245,3),idxm(9,3),I,j
!      common/potrot/fact1,fact2,c1,s1,icoord,xm(2),xmx,iperm   
! $$$      common/potmcm/xm(2)
!
!     expansion indicies
!

       real(ark) :: reoh = 0.958649_ark,&
      thetae=104.3475_ark,b1=2.0_ark,roh=0.9519607159623009_ark,alphaoh=2.587949757553683_ark,&
      deoh=42290.92019288289_ark,phh1=16.94879431193463_ark,phh2=12.66426998162947_ark


      real(ark) :: f5z=0.0_ark,fbasis=0.0_ark,fcore=-1.0_ark,frest=0.0_ark


      integer(ik) :: ifirst = 0

      data (idx(i,1),i=1,245)/                              &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, &
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
      2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
      3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, &
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
      4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, &
      4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, &
      6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
      6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, &
      5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, &
      7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, &
      6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, &
      9, 9, 9, 9, 9/
       data (idx(i,2),i=1,245)/                                  &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,&
      1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
      2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,&
      2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3,&
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
      2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,&
      3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,&
      1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3,&
      2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4,&
      4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,&
      2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,&
      4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1,&
      1, 1, 1, 1, 1/                                              
      data (idx(i,3),i=1,245)/                                   &
      1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 1, 2, 3, 4, 5,&
      6, 7, 8, 9,10,11,12,13,14, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,&
     12,13, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1, 2, 3, 4, 5,&
      6, 7, 8, 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 1,&
      2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,&
     11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8,&
      9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 7, 8,&
      9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9,&
      1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2,&
      3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6,&
      7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3,&
      4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2,&
      3, 4, 5, 6, 7/
!
!     expansion coefficients for 5z ab initio data
!
       data (c5z(i),i=1,45)/                                         &
      4.2278462684916D+04, 4.5859382909906D-02, 9.4804986183058D+03, &
      7.5485566680955D+02, 1.9865052511496D+03, 4.3768071560862D+02, &
      1.4466054104131D+03, 1.3591924557890D+02,-1.4299027252645D+03, &
      6.6966329416373D+02, 3.8065088734195D+03,-5.0582552618154D+02, &
     -3.2067534385604D+03, 6.9673382568135D+02, 1.6789085874578D+03, &
     -3.5387509130093D+03,-1.2902326455736D+04,-6.4271125232353D+03, &
     -6.9346876863641D+03,-4.9765266152649D+02,-3.4380943579627D+03, &
      3.9925274973255D+03,-1.2703668547457D+04,-1.5831591056092D+04, &
      2.9431777405339D+04, 2.5071411925779D+04,-4.8518811956397D+04, &
     -1.4430705306580D+04, 2.5844109323395D+04,-2.3371683301770D+03, &
      1.2333872678202D+04, 6.6525207018832D+03,-2.0884209672231D+03, &
     -6.3008463062877D+03, 4.2548148298119D+04, 2.1561445953347D+04, &
     -1.5517277060400D+05, 2.9277086555691D+04, 2.6154026873478D+05, &
     -1.3093666159230D+05,-1.6260425387088D+05, 1.2311652217133D+05, &
     -5.1764697159603D+04, 2.5287599662992D+03, 3.0114701659513D+04/ 
                                                                      
      data (c5z(i),i=46,90)/                                         &
     -2.0580084492150D+03, 3.3617940269402D+04, 1.3503379582016D+04, &
     -1.0401149481887D+05,-6.3248258344140D+04, 2.4576697811922D+05, &
      8.9685253338525D+04,-2.3910076031416D+05,-6.5265145723160D+04, &
      8.9184290973880D+04,-8.0850272976101D+03,-3.1054961140464D+04, &
     -1.3684354599285D+04, 9.3754012976495D+03,-7.4676475789329D+04, &
     -1.8122270942076D+05, 2.6987309391410D+05, 4.0582251904706D+05, &
     -4.7103517814752D+05,-3.6115503974010D+05, 3.2284775325099D+05, &
      1.3264691929787D+04, 1.8025253924335D+05,-1.2235925565102D+04, &
     -9.1363898120735D+03,-4.1294242946858D+04,-3.4995730900098D+04, &
      3.1769893347165D+05, 2.8395605362570D+05,-1.0784536354219D+06, &
     -5.9451106980882D+05, 1.5215430060937D+06, 4.5943167339298D+05, &
     -7.9957883936866D+05,-9.2432840622294D+04, 5.5825423140341D+03, &
      3.0673594098716D+03, 8.7439532014842D+04, 1.9113438435651D+05, &
     -3.4306742659939D+05,-3.0711488132651D+05, 6.2118702580693D+05, &
     -1.5805976377422D+04,-4.2038045404190D+05, 3.4847108834282D+05/ 
                                                                      
      data (c5z(i),i=91,135)/                                        &
     -1.3486811106770D+04, 3.1256632170871D+04, 5.3344700235019D+03, &
      2.6384242145376D+04, 1.2917121516510D+05,-1.3160848301195D+05, &
     -4.5853998051192D+05, 3.5760105069089D+05, 6.4570143281747D+05, &
     -3.6980075904167D+05,-3.2941029518332D+05,-3.5042507366553D+05, &
      2.1513919629391D+03, 6.3403845616538D+04, 6.2152822008047D+04, &
     -4.8805335375295D+05,-6.3261951398766D+05, 1.8433340786742D+06, &
      1.4650263449690D+06,-2.9204939728308D+06,-1.1011338105757D+06, &
      1.7270664922758D+06, 3.4925947462024D+05,-1.9526251371308D+04, &
     -3.2271030511683D+04,-3.7601575719875D+05, 1.8295007005531D+05, &
      1.5005699079799D+06,-1.2350076538617D+06,-1.8221938812193D+06, &
      1.5438780841786D+06,-3.2729150692367D+03, 1.0546285883943D+04, &
     -4.7118461673723D+04,-1.1458551385925D+05, 2.7704588008958D+05, &
      7.4145816862032D+05,-6.6864945408289D+05,-1.6992324545166D+06, &
      6.7487333473248D+05, 1.4361670430046D+06,-2.0837555267331D+05, &
      4.7678355561019D+05,-1.5194821786066D+04,-1.1987249931134D+05/ 
                                                                     
                                                                      
      data (c5z(i),i=136,180)/                                       &
      1.3007675671713D+05, 9.6641544907323D+05,-5.3379849922258D+05, &
     -2.4303858824867D+06, 1.5261649025605D+06, 2.0186755858342D+06, &
     -1.6429544469130D+06,-1.7921520714752D+04, 1.4125624734639D+04, &
     -2.5345006031695D+04, 1.7853375909076D+05,-5.4318156343922D+04, &
     -3.6889685715963D+05, 4.2449670705837D+05, 3.5020329799394D+05, &
      9.3825886484788D+03,-8.0012127425648D+05, 9.8554789856472D+04, &
      4.9210554266522D+05,-6.4038493953446D+05,-2.8398085766046D+06, &
      2.1390360019254D+06, 6.3452935017176D+06,-2.3677386290925D+06, &
     -3.9697874352050D+06,-1.9490691547041D+04, 4.4213579019433D+04, &
      1.6113884156437D+05,-7.1247665213713D+05,-1.1808376404616D+06, &
      3.0815171952564D+06, 1.3519809705593D+06,-3.4457898745450D+06, &
      2.0705775494050D+05,-4.3778169926622D+05, 8.7041260169714D+03, &
      1.8982512628535D+05,-2.9708215504578D+05,-8.8213012222074D+05, &
      8.6031109049755D+05, 1.0968800857081D+06,-1.0114716732602D+06, &
      1.9367263614108D+05, 2.8678295007137D+05,-9.4347729862989D+04/ 
                                                                      
      data (c5z(i),i=181,225)/                                       &
      4.4154039394108D+04, 5.3686756196439D+05, 1.7254041770855D+05, &
     -2.5310674462399D+06,-2.0381171865455D+06, 3.3780796258176D+06, &
      7.8836220768478D+05,-1.5307728782887D+05,-3.7573362053757D+05, &
      1.0124501604626D+06, 2.0929686545723D+06,-5.7305706586465D+06, &
     -2.6200352535413D+06, 7.1543745536691D+06,-1.9733601879064D+04, &
      8.5273008477607D+04, 6.1062454495045D+04,-2.2642508675984D+05, &
      2.4581653864150D+05,-9.0376851105383D+05,-4.4367930945690D+05, &
      1.5740351463593D+06, 2.4563041445249D+05,-3.4697646046367D+03, &
     -2.1391370322552D+05, 4.2358948404842D+05, 5.6270081955003D+05, &
     -8.5007851251980D+05,-6.1182429537130D+05, 5.6690751824341D+05, &
     -3.5617502919487D+05,-8.1875263381402D+02,-2.4506258140060D+05, &
      2.5830513731509D+05, 6.0646114465433D+05,-6.9676584616955D+05, &
      5.1937406389690D+05, 1.7261913546007D+05,-1.7405787307472D+04, &
     -3.8301842660567D+05, 5.4227693205154D+05, 2.5442083515211D+06, &
     -1.1837755702370D+06,-1.9381959088092D+06,-4.0642141553575D+05/ 
                                                                     
                                                                      
      data (c5z(i),i=226,245)/                                       &
      1.1840693827934D+04,-1.5334500255967D+05, 4.9098619510989D+05, &
      6.1688992640977D+05, 2.2351144690009D+05,-1.8550462739570D+06, &
      9.6815110649918D+03,-8.1526584681055D+04,-8.0810433155289D+04, &
      3.4520506615177D+05, 2.5509863381419D+05,-1.3331224992157D+05, &
     -4.3119301071653D+05,-5.9818343115856D+04, 1.7863692414573D+03, &
      8.9440694919836D+04,-2.5558967650731D+05,-2.2130423988459D+04, &
      4.4973674518316D+05,-2.2094939343618D+05/                      
!                                                                     
!    expansion coefficients for basis correction                     
!                                                                     
      data (cbasis(i),i=1,45)/                                       &
      6.9770019624764D-04,-2.4209870001642D+01, 1.8113927151562D+01, &
      3.5107416275981D+01,-5.4600021126735D+00,-4.8731149608386D+01, &
      3.6007189184766D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -7.7178474355102D+01,-3.8460795013977D+01,-4.6622480912340D+01, &
      5.5684951167513D+01, 1.2274939911242D+02,-1.4325154752086D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-6.0800589055949D+00, &
      8.6171499453475D+01,-8.4066835441327D+01,-5.8228085624620D+01, &
      2.0237393793875D+02, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      3.3525582670313D+02, 7.0056962392208D+01,-4.5312502936708D+01/ 
                                                                     
                                                                      
      data (cbasis(i),i=46,90)/                                      &
     -3.0441141194247D+02, 2.8111438108965D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.2983583774779D+02, 3.9781671212935D+01, &
     -6.6793945229609D+01,-1.9259805675433D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-8.2855757669957D+02,-5.7003072730941D+01, &
     -3.5604806670066D+01, 9.6277766002709D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 8.8645622149112D+02,-7.6908409772041D+01, &
      6.8111763314154D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                     
                                                                      
      data (cbasis(i),i=91,135)/                                     &
      2.5090493428062D+02,-2.3622141780572D+02, 5.8155647658455D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 2.8919570295095D+03, &
     -1.7871014635921D+02,-1.3515667622500D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-3.6965613754734D+03, 2.1148158286617D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-1.4795670139431D+03, &
      3.6210798138768D+02, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -5.3552886800881D+03, 3.1006384016202D+02, 0.0000000000000D+00/ 
                                                                     
                                                                     
                                                                     
                                                                      
      data (cbasis(i),i=136,180)/                                    &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 1.6241824368764D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 4.3764909606382D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 1.0940849243716D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 3.0743267832931D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (cbasis(i),i=181,225)/                                    &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (cbasis(i),i=226,245)/                                    &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00/                      
!                                                                     
!    expansion coefficients for core correction                      
!                                                                     
      data (ccore(i),i=1,45)/                                        &
      2.4332191647159D-02,-2.9749090113656D+01, 1.8638980892831D+01, &
     -6.1272361746520D+00, 2.1567487597605D+00,-1.5552044084945D+01, &
      8.9752150543954D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -3.5693557878741D+02,-3.0398393196894D+00,-6.5936553294576D+00, &
      1.6056619388911D+01, 7.8061422868204D+01,-8.6270891686359D+01, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-3.1688002530217D+01, &
      3.7586725583944D+01,-3.2725765966657D+01,-5.6458213299259D+00, &
      2.1502613314595D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      5.2789943583277D+02,-4.2461079404962D+00,-2.4937638543122D+01/ 
                                                                     
                                                                      
      data (ccore(i),i=46,90)/                                       &
     -1.1963809321312D+02, 2.0240663228078D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-6.2574211352272D+02,-6.9617539465382D+00, &
     -5.9440243471241D+01, 1.4944220180218D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.2851139918332D+03,-6.5043516710835D+00, &
      4.0410829440249D+01,-6.7162452402027D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 1.0031942127832D+03, 7.6137226541944D+01, &
     -2.7279242226902D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (ccore(i),i=91,135)/                                      &
     -3.3059000871075D+01, 2.4384498749480D+01,-1.4597931874215D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 1.6559579606045D+03, &
      1.5038996611400D+02,-7.3865347730818D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.9738401290808D+03,-1.4149993809415D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-1.2756627454888D+02, &
      4.1487702227579D+01, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -1.7406770966429D+03,-9.3812204399266D+01, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (ccore(i),i=136,180)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.1890301282216D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 2.3723447727360D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-1.0279968223292D+03, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 5.7153838472603D+02, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                      
      data (ccore(i),i=181,225)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                      
      data (ccore(i),i=226,245)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00/                      
!                                                                     
!    expansion coefficients for v rest                               
!                                                                     
      data (crest(i),i=1,45)/                                        &
      0.0000000000000D+00,-4.7430930170000D+00,-1.4422132560000D+01, &
     -1.8061146510000D+01, 7.5186735000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
     -2.7962099800000D+02, 1.7616414260000D+01,-9.9741392630000D+01, &
      7.1402447000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00,-7.8571336480000D+01, &
      5.2434353250000D+01, 7.7696745000000D+01, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      1.7799123760000D+02, 1.4564532380000D+02, 2.2347226000000D+02/ 
                                                                      
      data (crest(i),i=46,90)/                                       &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-4.3823284100000D+02,-7.2846553000000D+02, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00,-2.6752313750000D+02, 3.6170310000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                      
      data (crest(i),i=91,135)/                                      &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (crest(i),i=136,180)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (crest(i),i=181,225)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/ 
                                                                     
                                                                      
      data (crest(i),i=226,245)/                                     &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00, &
      0.0000000000000D+00, 0.0000000000000D+00/
!
!     expansion indicies for mass correction
!
       data idxm/1,2,1,1,3,2,1,2,1,&
                2,1,1,3,1,2,2,1,1, &
                1,1,2,1,1,1,2,2,3/
!
!     expansion coefficients for mass correction
!
       data cmass/ -8.3554183D+00,3.7036552D+01,-5.2722136D+00,&
           1.6843857D+01,-7.0929741D+01,5.5380337D+00,-2.9962997D+01,&
           1.3637682D+02,-3.0530195d+00/
!
!     two body parameters
!
!
!     scaling factors for contributions to emperical potential
!

      if(ifirst.eq.0)then
 !      ifirst=1
 !   1  format(/1x,'pes for h2o',&
 !            /1x,'by Harry Partridge and David W. Schwenke',&
 !            /1x,'submitted to J. Chem. Phys. Nov. 8, 1996')
 !  56  format(/1x,'parameters before adjustment')
 !  55  format(/1x,'two body potential parameters:',&
 !            /1x,'hh: phh1 = ',f10.1,' phh2 = ',f5.2,&
 !            /1x,'oh: deoh = ',f10.1,' alpha = ',f7.4,&
 !            ' re = ',f7.4)
 !   4  format(/1x,'three body parameters:',&
 !            /1x,'reoh = ',f10.4,' thetae = ',f10.4,&
 !            /1x,'betaoh = ',f10.4,&
 !            /1x,'    i    j    k',7x,'c5z',9x,'cbasis',10x,'ccore',&
 !            10x,'crest')
 !      do 2 i=1,245
 !   2  continue
!
!     remove mass correction from vrest
!
       xmh=1836.152697_ark
       xmhi=1.0_ark/xmh
       xmd=3670.483031_ark
       fact=1._ark/((1._ark/xmd)-(1._ark/xmh))
   65  format(/1x,'parameters for delta v hdo ',&
            /1x,'    i    j    k')
       do 60 i=1,9
        cmass(i)=cmass(i)*fact
        corr=cmass(i)*xmhi
        if(idxm(i,1).eq.idxm(i,2))corr=corr*0.5_ark
        do 61 j=1,245
         if(idx(j,1).eq.idxm(i,1).and.idx(j,2).eq.idxm(i,2).and.&
           idx(j,3).eq.idxm(i,3))then
          crest(j)=crest(j)-corr
          go to 62
         end if
   61   continue
   62   continue
        do 63 j=1,245
         if(idx(j,2).eq.idxm(i,1).and.idx(j,1).eq.idxm(i,2).and.&
           idx(j,3).eq.idxm(i,3))then
          crest(j)=crest(j)-corr
          go to 64
         end if
   63   continue
   64   continue
   60  continue
       xm1=1._ark/1836.152697_ark
       xm2=1._ark/1836.152697_ark
!
!     adjust parameters using scale factors
!
   57  format(/1x,'adjusting parameters using scale factors ',&
             /1x,'f5z =    ',f11.8,                           &
             /1x,'fbasis = ',f11.8,                           &
             /1x,'fcore =  ',f11.8,                           &
             /1x,'frest =  ',f11.8)
       phh1=phh1*f5z
       deoh=deoh*f5z
       do 59 i=1,245
        c5z(i)=f5z*c5z(i)+fbasis*cbasis(i)+fcore*ccore(i)&
            +frest*crest(i)
   59  continue
   58  format(/1x,'three body parameters:',         &
             /1x,'reoh = ',f10.4,' thetae = ',f10.4,&
             /1x,'betaoh = ',f10.4,                 &
             /1x,'    i    j    k   cijk',          &
             /(1x,3i5,1pe15.7))
       do 66 i=1,9
        cmass(i)=cmass(i)*frest
   66  continue
   76  format(/1x,'mass correction factors ',&
             /1x,'    i    j    k   cijk',   &
             /(1x,3i5,1pe15.7))
!
!     convert parameters from 1/cm, angstrom to a.u.
!
       reoh=reoh/0.529177249_ark
       b1=b1*0.529177249_ark*0.529177249_ark
       do 3 i=1,245
        c5z(i)=c5z(i)*4.556335e-6_ark
    3  continue
       do 67 i=1,9
        cmass(i)=cmass(i)*4.556335e-6_ark
   67  continue
       rad=acos(-1.0_ark)/180.0_ark
       ce=cos(thetae*rad)
       phh1=phh1*exp(phh2)
       phh1=phh1*4.556335e-6_ark
       phh2=phh2*0.529177249_ark
       deoh=deoh*4.556335e-6_ark
       roh=roh/0.529177249_ark
       alphaoh=alphaoh*0.529177249_ark
       c5z(1)=c5z(1)*2.0_ark
      end if
       x1=(rij1-reoh)/reoh
       x2=(rij2-reoh)/reoh
       x3=cos(rij3)-ce
       rhh=sqrt(rij1**2+rij2**2 &
           -2._ark*rij1*rij2*cos(rij3))
       vhh=phh1*exp(-phh2*rhh)
       ex=exp(-alphaoh*(rij1-roh))
       voh1=deoh*ex*(ex-2._ark)
       ex=exp(-alphaoh*(rij2-roh))
       voh2=deoh*ex*(ex-2._ark)
       fmat(1,1)=1._ark
       fmat(1,2)=1._ark
       fmat(1,3)=1._ark
       do 10 j=2,15
        fmat(j,1)=fmat(j-1,1)*x1
        fmat(j,2)=fmat(j-1,2)*x2
        fmat(j,3)=fmat(j-1,3)*x3
   10  continue
       v=0
       do 12 j=2,245
        term=c5z(j)*(fmat(idx(j,1),1)*fmat(idx(j,2),2)&
                         +fmat(idx(j,2),1)*fmat(idx(j,1),2))&
                         *fmat(idx(j,3),3)
        v=v+term
   12  continue
       v1=0
       v2=0
       do 13 j=1,9
        v1=v1+cmass(j)*fmat(idxm(j,1),1)*fmat(idxm(j,2),2)&
            *fmat(idxm(j,3),3)
        v2=v2+cmass(j)*fmat(idxm(j,2),1)*fmat(idxm(j,1),2)&
            *fmat(idxm(j,3),3)
   13  continue
       v=v+xm1*v1+xm2*v2

       v=v*exp(-b1*((rij1-reoh)**2+(rij2-reoh)**2))&
            +c5z(1)&
            +voh1+voh2+vhh

      return
      end SUBROUTINE cvps_v1

      subroutine poten(v,r1,r2,th)

      implicit none
      real(ark),intent(in) ::  r1,r2,th
      real(ark),intent(out)::  v
      real(ark)            :: reoh,thetae,b1,roh,alphaoh
      real(ark)            :: phh2,ut,x0,ut2,x02,xs1,xs2,xst,&
      rr1,rr2,xep1,xep2,xep3
      real(ark)            :: rhh,vhh,vpb1,vpb2,v0,vp1,vp2,vp3,&
      vp,vps1,vps2,y1,y2,y12,y22,voh1,voh2

      reoh  = 0.958649_ark
      thetae= 104.3475_ark
      b1=2.0_ark
      roh=0.951961_ark
      alphaoh=2.587949757553683_ark
      phh2=6.70164303995_ark
      ut=20._ark
      x0=2.5_ark
      ut2=20._ark
      x02=2.5_ark
        
      thetae=thetae*3.14159265358979312_ark*.00555555555555555555_ark

      xs1=(r1+r2)*0.5_ark-reoh
      xs2=(r1-r2)*0.5_ark
      xst=cos(th)-cos(thetae)

      rr1=r1-roh
      rr2=r2-roh

      xep1=exp(-2._ark*alphaoh*rr1)-2._ark*exp(-alphaoh*rr1)+1._ark
      xep2=exp(-2._ark*alphaoh*rr2)-2._ark*exp(-alphaoh*rr2)+1._ark
      xep3=exp(-b1*(rr1**2+rr2**2))
      rhh=sqrt(r1**2+r2**2-2._ark*r1*r2*cos(th))
      vhh=0.900642240911285975e6_ark*exp(-phh2*rhh)
      vpb1=0.518622556959170834e5_ark*xep1
      vpb2=0.518622556959170834e5_ark*xep2

      v0 = -0.42928021494267e+02_ark*xs1**0*xs2**0*xst**0
     vp1 = -0.89507721126119e+02_ark*xs1**0*xs2**0*xst**1              &
     -0.87092714754886E+04_ark*xs1**1*xs2**0*xst**0                     &
     +0.18555919554430E+05_ark*xs1**0*xs2**0*xst**2                     &
     -0.22925394921245E+06_ark*xs1**0*xs2**2*xst**0                     &
     -0.25435205159307E+05_ark*xs1**1*xs2**0*xst**1                     &
     -0.24218598651361E+06_ark*xs1**2*xs2**0*xst**0                     &
     +0.93530733075594E+03_ark*xs1**0*xs2**0*xst**3                     &
     -0.20899025564266E+05_ark*xs1**0*xs2**2*xst**1                     &
     -0.96989161752558E+04_ark*xs1**1*xs2**0*xst**2                     &
     +0.21885768944145E+07_ark*xs1**1*xs2**2*xst**0                     &
     +0.23437895856756E+05_ark*xs1**2*xs2**0*xst**1                     &
     +0.70093387436425E+06_ark*xs1**3*xs2**0*xst**0                     &
     +0.35885579708482E+04_ark*xs1**0*xs2**0*xst**4                     &
     +0.56490957166597E+05_ark*xs1**0*xs2**2*xst**2                     &
     -0.20017103092868E+07_ark*xs1**0*xs2**4*xst**0                     &
     -0.10570227592471E+05_ark*xs1**1*xs2**0*xst**3                     &
     -0.16691503150751E+05_ark*xs1**1*xs2**2*xst**1                     &
     +0.69141739428769E+05_ark*xs1**2*xs2**0*xst**2                     &
     -0.83293916823201E+07_ark*xs1**2*xs2**2*xst**0                     &
     -0.72636651757492E+05_ark*xs1**3*xs2**0*xst**1                     &
     -0.20426551065586E+07_ark*xs1**4*xs2**0*xst**0                     &
     +0.28062216488162E+03_ark*xs1**0*xs2**0*xst**5                     &
     +0.21549280957266E+04_ark*xs1**0*xs2**2*xst**3                     &
     -0.25683270446184E+05_ark*xs1**0*xs2**4*xst**1                     &
     +0.28441886109915E+04_ark*xs1**1*xs2**0*xst**4                     &
     -0.43182839988294E+05_ark*xs1**1*xs2**2*xst**2                     &
     +0.13134477476588E+08_ark*xs1**1*xs2**4*xst**0                     &
     -0.60836641332073E+04_ark*xs1**2*xs2**0*xst**3                     &
     -0.87481928936921E+05_ark*xs1**2*xs2**2*xst**1                     &
     -0.30990616415060E+05_ark*xs1**3*xs2**0*xst**2                     &
     +0.20290673964095E+08_ark*xs1**3*xs2**2*xst**0                     &
     +0.17703228499377E+05_ark*xs1**4*xs2**0*xst**1                     &
     +0.37866131292554E+07_ark*xs1**5*xs2**0*xst**0                     &
     +0.10242815610045E+03_ark*xs1**0*xs2**0*xst**6                     &
     +0.22514283131421E+05_ark*xs1**0*xs2**2*xst**4                     &
     +0.12031585456000E+06_ark*xs1**0*xs2**4*xst**2                     &
     -0.63040312191604E+07_ark*xs1**0*xs2**6*xst**0                     &
     +0.42288805070414E+03_ark*xs1**1*xs2**0*xst**5                     &
     -0.73804912826702E+05_ark*xs1**1*xs2**2*xst**3                     &
     +0.58368057796765E+05_ark*xs1**1*xs2**4*xst**1                     &
     -0.38837920975341E+03_ark*xs1**2*xs2**0*xst**4                     &
     +0.34503289504885E+06_ark*xs1**2*xs2**2*xst**2                     &
     -0.37831345076438E+08_ark*xs1**2*xs2**4*xst**0                     &
     +0.37293966483751E+04_ark*xs1**3*xs2**0*xst**3                     &
     -0.10619720711173E+06_ark*xs1**3*xs2**2*xst**1                     &
     +0.25704546621260E+05_ark*xs1**4*xs2**0*xst**2                     &
     -0.37682650264356E+08_ark*xs1**4*xs2**2*xst**0                     &
     +0.74993318397304E+05_ark*xs1**5*xs2**0*xst**1                     &
     -0.65551970043981E+07_ark*xs1**6*xs2**0*xst**0                     &
     -0.58711372389116E+03_ark*xs1**0*xs2**0*xst**7                     &
     -0.52599265336948E+05_ark*xs1**0*xs2**2*xst**5                     &
     -0.18401270013519E+06_ark*xs1**0*xs2**4*xst**3                     &
     -0.32489685212066E+06_ark*xs1**0*xs2**6*xst**1                     &
     +0.15087691116288E+05_ark*xs1**1*xs2**0*xst**6                     &
     +0.11023845357016E+06_ark*xs1**1*xs2**2*xst**4                     &
     -0.42021818220258E+06_ark*xs1**1*xs2**4*xst**2                     &
     +0.31696716507181E+08_ark*xs1**1*xs2**6*xst**0                     &
     -0.74924087286956E+05_ark*xs1**2*xs2**0*xst**5                     &
     -0.23207574225873E+04_ark*xs1**2*xs2**2*xst**3                     &
     -0.10860704762756E+06_ark*xs1**2*xs2**4*xst**1                     &
     +0.16494892377660E+06_ark*xs1**3*xs2**0*xst**4                     &
     -0.99307653490068E+06_ark*xs1**3*xs2**2*xst**2                     &
     +0.68372807732485E+08_ark*xs1**3*xs2**4*xst**0                     &
     -0.72192743968295E+05_ark*xs1**4*xs2**0*xst**3                     &
     +0.18242143618460E+07_ark*xs1**4*xs2**2*xst**1                     &
     +0.12198797535213E+06_ark*xs1**5*xs2**0*xst**2                     &
     +0.50352162447823E+08_ark*xs1**5*xs2**2*xst**0                     &
     -0.20073330708130E+06_ark*xs1**6*xs2**0*xst**1                     &
     +0.90942929153787E+07_ark*xs1**7*xs2**0*xst**0                     &
     +0.21582347445035E+04_ark*xs1**0*xs2**0*xst**8                     &
     +0.36648586520692E+05_ark*xs1**0*xs2**2*xst**6                     &
     +0.15779759583354E+05_ark*xs1**0*xs2**4*xst**4                     &
     -0.20334050535353E+06_ark*xs1**0*xs2**6*xst**2                     &
     -0.11420250004310E+08_ark*xs1**0*xs2**8*xst**0                     &
     -0.21794937410225E+05_ark*xs1**1*xs2**0*xst**7                     &
     +0.20811252208119E+05_ark*xs1**1*xs2**2*xst**5                     &
     +0.11630454025798E+07_ark*xs1**1*xs2**4*xst**3                     &
     -0.26762775135417E+07_ark*xs1**1*xs2**6*xst**1                     &
     +0.62062769406155E+05_ark*xs1**2*xs2**0*xst**6                     &
     -0.29607240633187E+06_ark*xs1**2*xs2**2*xst**4                     &
     -0.88550613198782E+03_ark*xs1**2*xs2**4*xst**2                     &
     -0.53011250897316E+08_ark*xs1**2*xs2**6*xst**0                     &
     -0.58016761944746E+05_ark*xs1**3*xs2**0*xst**5                     &
     +0.10778020346865E+07_ark*xs1**3*xs2**2*xst**3                     &
     -0.47549560490496E+05_ark*xs1**3*xs2**4*xst**1                     &
     +0.21144023010349E+05_ark*xs1**4*xs2**0*xst**4                     &
     +0.15001022817007E+07_ark*xs1**4*xs2**2*xst**2                     &
     -0.64501924207419E+08_ark*xs1**4*xs2**4*xst**0                     &
     -0.14626394795061E+06_ark*xs1**5*xs2**0*xst**3                     &
     -0.66285625567564E+05_ark*xs1**5*xs2**2*xst**1                     &
     +0.31990537126595E+06_ark*xs1**6*xs2**0*xst**2                      
     vp2 = -0.46302498058330E+08_ark*xs1**6*xs2**2*xst**0               &
     +0.13170608713587E+06_ark*xs1**7*xs2**0*xst**1                     &
     -0.11158698289618E+08_ark*xs1**8*xs2**0*xst**0                     &
     -0.10890241724383E+04_ark*xs1**0*xs2**0*xst**9                     &
     +0.11247008313740E+06_ark*xs1**0*xs2**2*xst**7                     &
     +0.40543010150418E+06_ark*xs1**0*xs2**4*xst**5                     &
     -0.24242500577626E+05_ark*xs1**0*xs2**6*xst**3                     &
     +0.38176989518620E+07_ark*xs1**0*xs2**8*xst**1                     &
     -0.16250026069957E+05_ark*xs1**1*xs2**0*xst**8                     &
     -0.35068589136732E+06_ark*xs1**1*xs2**2*xst**6                     &
     +0.70301568960516E+06_ark*xs1**1*xs2**4*xst**4                     &
     +0.16191272634392E+05_ark*xs1**1*xs2**6*xst**2                     &
     +0.75512446037555E+04_ark*xs1**1*xs2**8*xst**0                     &
     +0.14520823268926E+06_ark*xs1**2*xs2**0*xst**7                     &
     +0.85941012656887E+06_ark*xs1**2*xs2**2*xst**5                     &
     -0.53580379579664E+07_ark*xs1**2*xs2**4*xst**3                     &
     +0.12582241754270E+08_ark*xs1**2*xs2**6*xst**1                     &
     -0.37131351301740E+06_ark*xs1**3*xs2**0*xst**6                     &
     -0.47064889925520E+06_ark*xs1**3*xs2**2*xst**4                     &
     +0.10422212766918E+08_ark*xs1**3*xs2**4*xst**2                     &
     -0.73996717187917E+08_ark*xs1**3*xs2**6*xst**0                     &
     +0.25925093192681E+06_ark*xs1**4*xs2**0*xst**5                     &
     -0.34056704087033E+07_ark*xs1**4*xs2**2*xst**3                     &
     -0.17870477975560E+05_ark*xs1**4*xs2**4*xst**1                     &
     -0.29019095537882E+06_ark*xs1**5*xs2**0*xst**4                     &
     +0.10456976435106E+08_ark*xs1**5*xs2**2*xst**2                     &
     -0.71012740172848E+08_ark*xs1**5*xs2**4*xst**0                     &
     -0.64758758551951E+06_ark*xs1**6*xs2**0*xst**3                     &
     -0.34346175712972E+08_ark*xs1**6*xs2**2*xst**1                     &
     +0.12569402016238E+06_ark*xs1**7*xs2**0*xst**2                     &
     +0.18473754829902E+08_ark*xs1**7*xs2**2*xst**0                     &
     -0.35467027561537E+06_ark*xs1**8*xs2**0*xst**1                     &
     +0.97227821072430E+07_ark*xs1**9*xs2**0*xst**0                     &
     -0.21596820455183E+04_ark*xs1**0*xs2**0*xst**10                    &
     -0.62753040823736E+05_ark*xs1**0*xs2**2*xst**8                     &
     -0.15122707809662E+06_ark*xs1**0*xs2**4*xst**6                     &
     +0.28605050468368E+05_ark*xs1**0*xs2**6*xst**4                     &
     +0.64423453064932E+07_ark*xs1**0*xs2**8*xst**2                     &
     +0.13276959048016E+05_ark*xs1**0*xs2**10*xst**0                    &
     +0.24070706739293E+05_ark*xs1**1*xs2**0*xst**9                     &
     +0.81596833267886E+05_ark*xs1**1*xs2**2*xst**7                     &
     -0.28559429016717E+07_ark*xs1**1*xs2**4*xst**5                     &
     +0.17758228492946E+07_ark*xs1**1*xs2**6*xst**3                     &
     +0.47618499828633E+08_ark*xs1**1*xs2**8*xst**1                     &
     -0.14338513851675E+06_ark*xs1**2*xs2**0*xst**8                     &
     +0.88046523013615E+06_ark*xs1**2*xs2**2*xst**6                     &
     +0.43632961649509E+05_ark*xs1**2*xs2**4*xst**4                     &
     +0.40855428706867E+08_ark*xs1**2*xs2**6*xst**2                     &
     +0.15590052161742E+09_ark*xs1**2*xs2**8*xst**0                     &
     +0.23839537263228E+06_ark*xs1**3*xs2**0*xst**7                     &
     -0.30986527433149E+07_ark*xs1**3*xs2**2*xst**5                     &
     +0.81912344699208E+07_ark*xs1**3*xs2**4*xst**3                     &
     -0.14524552903320E+09_ark*xs1**3*xs2**6*xst**1                     &
     -0.27611309901574E+06_ark*xs1**4*xs2**0*xst**6                     &
     +0.27086691982660E+07_ark*xs1**4*xs2**2*xst**4                     &
     -0.48766159858853E+08_ark*xs1**4*xs2**4*xst**2                     &
     +0.64814047369542E+09_ark*xs1**4*xs2**6*xst**0                     &
     +0.10270479409000E+07_ark*xs1**5*xs2**0*xst**5                     &
     -0.61301540395406E+07_ark*xs1**5*xs2**2*xst**3                     &
     +0.48355391607449E+07_ark*xs1**5*xs2**4*xst**1                     &
     -0.16950330704990E+07_ark*xs1**6*xs2**0*xst**4                     &
     -0.36875318775114E+08_ark*xs1**6*xs2**2*xst**2                     &
     +0.32375436315231E+09_ark*xs1**6*xs2**4*xst**0                     &
     +0.20890194939845E+07_ark*xs1**7*xs2**0*xst**3                     &
     +0.85634762713006E+08_ark*xs1**7*xs2**2*xst**1                     &
     -0.28453916840284E+07_ark*xs1**8*xs2**0*xst**2                     &
     -0.70157397242334E+05_ark*xs1**8*xs2**2*xst**0                     &
     +0.16590097494102E+06_ark*xs1**9*xs2**0*xst**1                     &
     -0.46100413467387E+07_ark*xs1**10*xs2**0*xst**0                    &
     +0.12990573820879E+04_ark*xs1**0*xs2**0*xst**11                    &
     -0.81875873232894E+05_ark*xs1**0*xs2**2*xst**9                     &
     -0.58457711341525E+06_ark*xs1**0*xs2**4*xst**7                     &
     +0.11599991081938E+07_ark*xs1**0*xs2**6*xst**5                     &
     -0.81818416797062E+07_ark*xs1**0*xs2**8*xst**3                     &
     -0.38544363585196E+08_ark*xs1**0*xs2**10*xst**1                    &
     +0.14380329527025E+05_ark*xs1**1*xs2**0*xst**10                    &
     +0.45495216378088E+06_ark*xs1**1*xs2**2*xst**8                     &
     +0.80001586002240E+06_ark*xs1**1*xs2**4*xst**6                     &
     +0.35942614195816E+05_ark*xs1**1*xs2**6*xst**4                     &
     -0.92442593621709E+08_ark*xs1**1*xs2**8*xst**2                     &
     +0.43271316229640E+08_ark*xs1**1*xs2**10*xst**0                    &
     -0.17425299721619E+06_ark*xs1**2*xs2**0*xst**9                      
     vp3 = -0.14264905533936E+07_ark*xs1**2*xs2**2*xst**7               &
     +0.40545757182569E+07_ark*xs1**2*xs2**4*xst**5                     &
     +0.17589297120178E+05_ark*xs1**2*xs2**6*xst**3                     &
     -0.10276064043441E+09_ark*xs1**2*xs2**8*xst**1                     &
     +0.56784453470325E+06_ark*xs1**3*xs2**0*xst**8                     &
     +0.18439783794232E+07_ark*xs1**3*xs2**2*xst**6                     &
     -0.12864927415382E+08_ark*xs1**3*xs2**4*xst**4                     &
     -0.10451573963676E+09_ark*xs1**3*xs2**6*xst**2                     &
     -0.43651679281752E+09_ark*xs1**3*xs2**8*xst**0                     &
     -0.73109417390594E+06_ark*xs1**4*xs2**0*xst**7                     &
     -0.17007594294248E+05_ark*xs1**4*xs2**2*xst**5                     &
     +0.59610187456335E+08_ark*xs1**4*xs2**4*xst**3                     &
     +0.49938979992303E+09_ark*xs1**4*xs2**6*xst**1                     &
     +0.68154313082103E+06_ark*xs1**5*xs2**0*xst**6                     &
     -0.59994121510676E+07_ark*xs1**5*xs2**2*xst**4                     &
     +0.58530065307103E+04_ark*xs1**5*xs2**4*xst**2                     &
     -0.14136428128750E+10_ark*xs1**5*xs2**6*xst**0                     &
     +0.94980927547676E+05_ark*xs1**6*xs2**0*xst**5                     &
     +0.50803781774926E+08_ark*xs1**6*xs2**2*xst**3                     &
     -0.17599150824124E+05_ark*xs1**6*xs2**4*xst**1                     &
     +0.24656468761383E+07_ark*xs1**7*xs2**0*xst**4                     &
     +0.27067872985736E+08_ark*xs1**7*xs2**2*xst**2                     &
     -0.30519784547178E+09_ark*xs1**7*xs2**4*xst**0                     &
     +0.99849050749970E+05_ark*xs1**8*xs2**0*xst**3                     &
     -0.58311940782272E+08_ark*xs1**8*xs2**2*xst**1                     &
     +0.31058062408740E+07_ark*xs1**9*xs2**0*xst**2                     &
     +0.33552612558884E+05_ark*xs1**9*xs2**2*xst**0                     &
     +0.59946019331995E+03_ark*xs1**0*xs2**0*xst**12                    &
     +0.66867268726893E+05_ark*xs1**0*xs2**2*xst**10                    &
     +0.30700337437714E+06_ark*xs1**0*xs2**4*xst**8                     &
     -0.38988317039688E+06_ark*xs1**0*xs2**6*xst**6                     &
     +0.65081830594207E+04_ark*xs1**0*xs2**8*xst**4                     &
     +0.33025369843350E+08_ark*xs1**0*xs2**10*xst**2                    &
     -0.17164908243459E+08_ark*xs1**0*xs2**12*xst**0                    &
     -0.12894279764727E+05_ark*xs1**1*xs2**0*xst**11                    &
     -0.33457489081673E+06_ark*xs1**1*xs2**2*xst**9                     &
     +0.10267124768913E+07_ark*xs1**1*xs2**4*xst**7                     &
     -0.39836458160981E+07_ark*xs1**1*xs2**6*xst**5                     &
     +0.34049617803398E+08_ark*xs1**1*xs2**8*xst**3                     &
     -0.21563575639396E+04_ark*xs1**1*xs2**10*xst**1                    &
     +0.14601592475540E+06_ark*xs1**2*xs2**0*xst**10                    &
     +0.59457412855332E+06_ark*xs1**2*xs2**2*xst**8                     &
     -0.23885890796883E+07_ark*xs1**2*xs2**4*xst**6                     &
     +0.50371609889570E+07_ark*xs1**2*xs2**6*xst**4                     &
     +0.96234635627910E+08_ark*xs1**2*xs2**8*xst**2                     &
     -0.88806691563152E+08_ark*xs1**2*xs2**10*xst**0                    &
     -0.43381503884347E+06_ark*xs1**3*xs2**0*xst**9                     &
     +0.33693022101345E+05_ark*xs1**3*xs2**2*xst**7                     &
     +0.88181499653639E+05_ark*xs1**3*xs2**4*xst**5                     &
     -0.44813806477506E+08_ark*xs1**3*xs2**6*xst**3                     &
     +0.11963770272834E+09_ark*xs1**3*xs2**8*xst**1                     &
     +0.59083684026351E+06_ark*xs1**4*xs2**0*xst**8                     &
     +0.12302081530284E+06_ark*xs1**4*xs2**2*xst**6                     &
     +0.18693218016832E+08_ark*xs1**4*xs2**4*xst**4                     &
     +0.13670141780925E+09_ark*xs1**4*xs2**6*xst**2                     &
     +0.41735524044315E+09_ark*xs1**4*xs2**8*xst**0                     &
     -0.11246548344715E+07_ark*xs1**5*xs2**0*xst**7                     &
     +0.17460150703783E+07_ark*xs1**5*xs2**2*xst**5                     &
     -0.11089446386451E+09_ark*xs1**5*xs2**4*xst**3                     &
     -0.58155496994476E+09_ark*xs1**5*xs2**6*xst**1                     &
     +0.14311855842965E+07_ark*xs1**6*xs2**0*xst**6                     &
     +0.26916434283893E+06_ark*xs1**6*xs2**2*xst**4                     &
     +0.10382617106543E+09_ark*xs1**6*xs2**4*xst**2                     &
     +0.10618975562842E+10_ark*xs1**6*xs2**6*xst**0                     &
     -0.34896791481610E+07_ark*xs1**7*xs2**0*xst**5                     &
     -0.50144884976014E+08_ark*xs1**7*xs2**2*xst**3

       vp=vp1+vp2+vp3

       vps1=42395.535333_ark*xep1
       vps2=42395.535333_ark*xep2

       y1=1._ark/(1._ark+exp(ut*(x0-r1)))
       y2=1._ark/(1._ark+exp(ut*(x0-r2)))
       y12=1._ark/(1._ark+exp(ut2*(x02-r1)))
       y22=1._ark/(1._ark+exp(ut2*(x02-r2)))

       vp=vp*xep3*(1.0_ark-y12)*(1.0_ark-y22)
       voh1=vpb1*(1.0_ark-y1)+y1*vps1
       voh2=vpb2*(1.0_ark-y2)+y2*vps2

       v=v0+vp+voh1+voh2+vhh

      end subroutine poten


subroutine DIPSA(muY, muX, r1, r2, cos_theta)

implicit none
real(ark), parameter :: reoh=1.8141
real(ark), intent(in) :: r1, r2, cos_theta
real(ark), intent(out) :: muX, muY
real(ark) :: xs1, xs2, xst, theta,thetae
real(ark) :: muX_0,muX_1,muX_2,muX_3,muY_1,muY_2,muY_3
!**********************************************************************************************


!***************                                                     
! FITTING VARIABLES                                                  
theta = acos(cos_theta)                                             
!***************
thetae=104.52_ark
thetae=thetae*PI/180.0_ark
xs1=(r1 + r2)/2.0_ark - reoh 
xs2=(r2 - r1)
xst=(theta)/(thetae)
!***************


   muX_0= 0.139000343502E+00_ark*xs1** 0*xs2**0*xst** 0
   muX_1= 0.138574906994E+01_ark*xs1** 0*xs2**0*xst** 1 + &
         -0.300882192853E+01_ark*xs1** 0*xs2**0*xst** 2 + &  
          0.397879622480E+01_ark*xs1** 0*xs2**0*xst** 3 + &  
         -0.291587470288E+01_ark*xs1** 0*xs2**0*xst** 4 + &  
          0.118249243990E+01_ark*xs1** 0*xs2**0*xst** 5 + & 
         -0.208246984853E+00_ark*xs1** 0*xs2**0*xst** 6 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst** 7 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst** 8 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst** 9 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**10 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**0*xst**18 + &  
          0.290040185310E+00_ark*xs1** 0*xs2**2*xst** 0 + &  
         -0.123168219804E+01_ark*xs1** 0*xs2**2*xst** 1 + &  
          0.149029869605E+01_ark*xs1** 0*xs2**2*xst** 2 + &  
          0.189765834651E+00_ark*xs1** 0*xs2**2*xst** 3 + &  
         -0.163365957223E+01_ark*xs1** 0*xs2**2*xst** 4 + &  
          0.116165044956E+01_ark*xs1** 0*xs2**2*xst** 5 + &  
         -0.258202260189E+00_ark*xs1** 0*xs2**2*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**10 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**11 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**12 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**13 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**16 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**2*xst**18 + &  
         -0.517088486321E-02_ark*xs1** 0*xs2**4*xst** 0 + &  
          0.585510437756E-01_ark*xs1** 0*xs2**4*xst** 1 + &  
         -0.937057483869E-01_ark*xs1** 0*xs2**4*xst** 2 + &  
          0.392174690408E-01_ark*xs1** 0*xs2**4*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**10 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**11 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**12 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**13 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**14 + & 
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**15 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**16 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**17 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**4*xst**18 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**10 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**11 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**12 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**13 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**14 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**15 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**16 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**17 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**6*xst**18 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**10 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**11 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**12 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**13 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**14 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**15 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**16 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**17 + &  
          0.000000000000E+00_ark*xs1** 0*xs2**8*xst**18 + &  
         -0.305455565784E+00_ark*xs1** 1*xs2**0*xst** 0 + &  
          0.143225159740E+01_ark*xs1** 1*xs2**0*xst** 1 + &  
         -0.179873351133E+01_ark*xs1** 1*xs2**0*xst** 2 + &  
          0.361285014190E-01_ark*xs1** 1*xs2**0*xst** 3 + &  
          0.121381033278E+01_ark*xs1** 1*xs2**0*xst** 4 + &  
          0.350219823675E+00_ark*xs1** 1*xs2**0*xst** 5 + &  
         -0.821617366273E+00_ark*xs1** 1*xs2**0*xst** 6 + &  
         -0.780012193064E+00_ark*xs1** 1*xs2**0*xst** 7 + &  
          0.284406102097E+00_ark*xs1** 1*xs2**0*xst** 8 + &  
          0.904068611959E+00_ark*xs1** 1*xs2**0*xst** 9 + &  
          0.526099119854E-01_ark*xs1** 1*xs2**0*xst**10 + &  
         -0.100629139951E+01_ark*xs1** 1*xs2**0*xst**11 + &  
          0.621691151627E+00_ark*xs1** 1*xs2**0*xst**12 + &  
         -0.116400665247E+00_ark*xs1** 1*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**0*xst**18 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst** 0 + & 
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst** 1 + & 
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst** 2 + & 
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst** 3 + & 
         -0.304372058999E+00_ark*xs1** 1*xs2**2*xst** 4 + & 
          0.299504907067E+00_ark*xs1** 1*xs2**2*xst** 5 + &  
          0.168893508692E+00_ark*xs1** 1*xs2**2*xst** 6 + &  
         -0.151642390849E+00_ark*xs1** 1*xs2**2*xst** 7 + &  
         -0.200741006292E+00_ark*xs1** 1*xs2**2*xst** 8 + &  
          0.190040561011E+00_ark*xs1** 1*xs2**2*xst** 9 + &  
         -0.398120120846E-01_ark*xs1** 1*xs2**2*xst**10 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst**11 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst**12 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst**13 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst**16 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**2*xst**18 + &  
         -0.159015064368E+00_ark*xs1** 1*xs2**4*xst** 0 + &  
          0.278888849567E+00_ark*xs1** 1*xs2**4*xst** 1 + &  
         -0.122466328779E+00_ark*xs1** 1*xs2**4*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**10 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**11 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**12 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**13 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**14 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**15 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**16 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**17 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**4*xst**18 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**10 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**11 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**12 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**13 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**14 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**15 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**16 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**17 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**6*xst**18 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst** 0 + &  
          0.704498014410E-02_ark*xs1** 1*xs2**8*xst** 1 + &  
         -0.125075640552E-01_ark*xs1** 1*xs2**8*xst** 2 + &  
          0.551468799951E-02_ark*xs1** 1*xs2**8*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**10 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**11 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**12 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**13 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**14 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**15 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**16 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**17 + &  
          0.000000000000E+00_ark*xs1** 1*xs2**8*xst**18 + &  
          0.622363248261E-02_ark*xs1** 2*xs2**0*xst** 0 + &  
         -0.909962620269E-01_ark*xs1** 2*xs2**0*xst** 1 + &  
         -0.208347403599E+00_ark*xs1** 2*xs2**0*xst** 2 + &  
          0.280120089326E+00_ark*xs1** 2*xs2**0*xst** 3 + &  
          0.256008073086E+00_ark*xs1** 2*xs2**0*xst** 4 + &  
         -0.139742857685E+00_ark*xs1** 2*xs2**0*xst** 5 + &  
         -0.344164584577E+00_ark*xs1** 2*xs2**0*xst** 6 + &  
         -0.135670673774E+00_ark*xs1** 2*xs2**0*xst** 7 + &  
          0.222346353864E+00_ark*xs1** 2*xs2**0*xst** 8 + &  
          0.285171468697E+00_ark*xs1** 2*xs2**0*xst** 9 + &  
         -0.788565777214E-01_ark*xs1** 2*xs2**0*xst**10 + &  
         -0.338552671613E+00_ark*xs1** 2*xs2**0*xst**11 + &  
          0.251924640052E+00_ark*xs1** 2*xs2**0*xst**12 + &  
         -0.518003936041E-01_ark*xs1** 2*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**0*xst**14 + &   
          0.000000000000E+00_ark*xs1** 2*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**0*xst**18 + &  
         -0.165949619855E-01_ark*xs1** 2*xs2**2*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst** 4 + &  
         -0.886382641194E-02_ark*xs1** 2*xs2**2*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**10 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**11 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**12 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**13 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**16 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**2*xst**18 + &  
         -0.998370308325E-01_ark*xs1** 2*xs2**4*xst** 0 + &  
          0.385931713053E+00_ark*xs1** 2*xs2**4*xst** 1 + &  
          0.542057636567E-02_ark*xs1** 2*xs2**4*xst** 2 + &  
         -0.215158536386E+00_ark*xs1** 2*xs2**4*xst** 3 + &  
         -0.148084597240E+00_ark*xs1** 2*xs2**4*xst** 4 + &  
          0.277951230475E-01_ark*xs1** 2*xs2**4*xst** 5 + &  
          0.107041579749E+00_ark*xs1** 2*xs2**4*xst** 6 + &  
          0.230703510631E-01_ark*xs1** 2*xs2**4*xst** 7 + &  
         -0.347110605651E-01_ark*xs1** 2*xs2**4*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**10 + &   
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**11 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**12 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**13 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**14 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**15 + &   
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**16 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**17 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**4*xst**18 + &  
          0.186602502199E+00_ark*xs1** 2*xs2**6*xst** 0 + &  
         -0.473031863902E+00_ark*xs1** 2*xs2**6*xst** 1 + &  
          0.156712006318E+00_ark*xs1** 2*xs2**6*xst** 2 + &  
          0.291681391439E+00_ark*xs1** 2*xs2**6*xst** 3 + &  
         -0.188889482586E+00_ark*xs1** 2*xs2**6*xst** 4 + &  
          0.158592455531E-01_ark*xs1** 2*xs2**6*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**10 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**11 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**12 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**13 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**14 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**15 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**16 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**17 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**6*xst**18 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 7 + &   
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**10 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**11 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**12 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**13 + &   
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**14 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**15 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**16 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**17 + &  
          0.000000000000E+00_ark*xs1** 2*xs2**8*xst**18 + &  
          0.204287827229E+00_ark*xs1** 3*xs2**0*xst** 0 + &   
         -0.533159873307E+00_ark*xs1** 3*xs2**0*xst** 1 + &   
          0.137916441651E+00_ark*xs1** 3*xs2**0*xst** 2 + &  
          0.407363661916E+00_ark*xs1** 3*xs2**0*xst** 3 + &  
          0.677423729406E-01_ark*xs1** 3*xs2**0*xst** 4 + &  
         -0.320550695638E+00_ark*xs1** 3*xs2**0*xst** 5 + &  
         -0.314366143945E+00_ark*xs1** 3*xs2**0*xst** 6 + &  
          0.555624903519E-01_ark*xs1** 3*xs2**0*xst** 7 + &  
          0.361543494572E+00_ark*xs1** 3*xs2**0*xst** 8 + &  
          0.158021897815E+00_ark*xs1** 3*xs2**0*xst** 9 + &  
         -0.411077313819E+00_ark*xs1** 3*xs2**0*xst**10 + &  
          0.134478905438E+00_ark*xs1** 3*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**0*xst**18 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 8   
   muX_2= 0.000000000000E+00_ark*xs1** 3*xs2**2*xst** 9 + &
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**10 + & 
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**11 + & 
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**12 + & 
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**13 + & 
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**16 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**2*xst**18 + &  
          0.799596397128E-01_ark*xs1** 3*xs2**4*xst** 0 + &  
         -0.306381016877E-01_ark*xs1** 3*xs2**4*xst** 1 + &  
         -0.320171775471E+00_ark*xs1** 3*xs2**4*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst** 5 + &  
          0.462720883432E+00_ark*xs1** 3*xs2**4*xst** 6 + &  
          0.184501395661E-01_ark*xs1** 3*xs2**4*xst** 7 + &  
         -0.465593360089E+00_ark*xs1** 3*xs2**4*xst** 8 + &  
          0.184574969420E+00_ark*xs1** 3*xs2**4*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**10 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**11 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**12 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**13 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**14 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**15 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**16 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**17 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**4*xst**18 + &  
         -0.689069174614E-01_ark*xs1** 3*xs2**6*xst** 0 + &  
         -0.316877166195E-01_ark*xs1** 3*xs2**6*xst** 1 + &  
          0.318293598146E+00_ark*xs1** 3*xs2**6*xst** 2 + &  
          0.594445500995E-01_ark*xs1** 3*xs2**6*xst** 3 + &  
         -0.482594663477E+00_ark*xs1** 3*xs2**6*xst** 4 + &  
          0.222047829350E+00_ark*xs1** 3*xs2**6*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**10 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**11 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**12 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**13 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**14 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**15 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**16 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**17 + &  
          0.000000000000E+00_ark*xs1** 3*xs2**6*xst**18 + &  
          0.710583732361E-02_ark*xs1** 4*xs2**0*xst** 0 + &  
          0.731512912010E-02_ark*xs1** 4*xs2**0*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**10 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**0*xst**18 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst** 6 + &  
          0.132888641167E+00_ark*xs1** 4*xs2**2*xst** 7 + &  
         -0.672463756135E-01_ark*xs1** 4*xs2**2*xst** 8 + &  
         -0.108328088594E-01_ark*xs1** 4*xs2**2*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**10 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**11 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**12 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**13 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**16 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**2*xst**18 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 4 + &   
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**10 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**11 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**12 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**13 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**14 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**15 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**16 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**17 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**4*xst**18 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**10 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**11 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**12 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**13 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**14 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**15 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**16 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**17 + &  
          0.000000000000E+00_ark*xs1** 4*xs2**6*xst**18 + &  
          0.310491121396E-01_ark*xs1** 5*xs2**0*xst** 0 + &   
          0.106935184487E+00_ark*xs1** 5*xs2**0*xst** 1 + &  
         -0.169245445911E-02_ark*xs1** 5*xs2**0*xst** 2 + &  
         -0.165613457554E+00_ark*xs1** 5*xs2**0*xst** 3 + &  
         -0.132673705194E+00_ark*xs1** 5*xs2**0*xst** 4 + &  
          0.172085738190E+00_ark*xs1** 5*xs2**0*xst** 5 + &   
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**10 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**0*xst**18 + &  
         -0.300660763288E+00_ark*xs1** 5*xs2**2*xst** 0 + &  
          0.622406405000E+00_ark*xs1** 5*xs2**2*xst** 1 + &  
          0.350823666478E-02_ark*xs1** 5*xs2**2*xst** 2 + &  
         -0.321849807044E+00_ark*xs1** 5*xs2**2*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**10 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**11 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**12 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**13 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**16 + &   
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**2*xst**18 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 3 + &   
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst** 9 + &   
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**10 + &   
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**11 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**12 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**13 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**14 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**15 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**16 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**17 + &  
          0.000000000000E+00_ark*xs1** 5*xs2**4*xst**18 + &  
         -0.174485954068E+00_ark*xs1** 6*xs2**0*xst** 0 + &  
          0.538495570600E-01_ark*xs1** 6*xs2**0*xst** 1 + &  
          0.423006244628E-02_ark*xs1** 6*xs2**0*xst** 2 + &  
          0.186470338707E-01_ark*xs1** 6*xs2**0*xst** 3 + &  
          0.186308233494E+00_ark*xs1** 6*xs2**0*xst** 4 + &  
          0.264248866768E+00_ark*xs1** 6*xs2**0*xst** 5 + &  
         -0.311734472051E+00_ark*xs1** 6*xs2**0*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst**10 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**0*xst**17   
   muX_3= 0.000000000000E+00_ark*xs1** 6*xs2**0*xst**18 + &
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 0 + & 
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 1 + & 
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 2 + & 
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 3 + & 
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 4 + & 
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**10 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**11 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**12 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**13 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**16 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**2*xst**18 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**10 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**11 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**12 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**13 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**14 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**15 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**16 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**17 + &  
          0.000000000000E+00_ark*xs1** 6*xs2**4*xst**18 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**10 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**0*xst**18 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**10 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**11 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**12 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**13 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**16 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 7*xs2**2*xst**18 + &  
          0.169009710399E+00_ark*xs1** 8*xs2**0*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst** 2 + &  
         -0.862073253838E+00_ark*xs1** 8*xs2**0*xst** 3 + &  
          0.632703184040E+00_ark*xs1** 8*xs2**0*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**10 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**14 + &   
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**0*xst**18 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**10 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**11 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**12 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**13 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**14 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**15 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**16 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**17 + &  
          0.000000000000E+00_ark*xs1** 8*xs2**2*xst**18 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 0 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 1 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 2 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 3 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 4 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 5 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 6 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 7 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 8 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst** 9 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**10 + &   
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**15 + &   
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1** 9*xs2**0*xst**18 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 0 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 1 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 2 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 3 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 4 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 5 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 6 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 7 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 8 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst** 9 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**10 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**11 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**12 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**13 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**14 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**15 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**16 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**17 + &  
          0.000000000000E+00_ark*xs1**10*xs2**0*xst**18  

                                                                        
!********************************************************************** **********************
muY_1=-0.285025557823E+00_ark*xs1** 0*xs2** 1*xst**0 +  &
       0.173926920509E+01_ark*xs1** 0*xs2** 1*xst** 1 +  &   
      -0.292195508049E+01_ark*xs1** 0*xs2** 1*xst** 2 +  &   
       0.128729522627E+01_ark*xs1** 0*xs2** 1*xst** 3 +  &   
       0.537564488620E+00_ark*xs1** 0*xs2** 1*xst** 4 +  &   
      -0.593456967270E+00_ark*xs1** 0*xs2** 1*xst** 5 +  &   
       0.289151059413E-01_ark*xs1** 0*xs2** 1*xst** 6 +  &   
       0.521985276981E-01_ark*xs1** 0*xs2** 1*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 1*xst**18 +  &   
       0.128500362604E+00_ark*xs1** 0*xs2** 3*xst** 0 +  &   
      -0.440994978773E+00_ark*xs1** 0*xs2** 3*xst** 1 +  &   
       0.456629056581E+00_ark*xs1** 0*xs2** 3*xst** 2 +  &   
      -0.143364107140E+00_ark*xs1** 0*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 3*xst**18 +  &   
       0.111734117326E-01_ark*xs1** 0*xs2** 5*xst** 0 +  &   
      -0.227062348625E-01_ark*xs1** 0*xs2** 5*xst** 1 +  &   
       0.127802587954E-01_ark*xs1** 0*xs2** 5*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 5*xst**18 +  &   
       0.929529605962E-04_ark*xs1** 0*xs2** 7*xst** 0 +  &   
      -0.380947184546E-04_ark*xs1** 0*xs2** 7*xst** 1 +  &   
      -0.730036686921E-03_ark*xs1** 0*xs2** 7*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 7*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 0*xs2** 9*xst**18 +  &   
      -0.123929332186E+00_ark*xs1** 1*xs2** 1*xst** 0 +  &   
       0.114800455836E+01_ark*xs1** 1*xs2** 1*xst** 1 +  &   
      -0.282171335901E+01_ark*xs1** 1*xs2** 1*xst** 2 +  &   
       0.269205302958E+01_ark*xs1** 1*xs2** 1*xst** 3 +  &   
       0.137466845362E+01_ark*xs1** 1*xs2** 1*xst** 4 +  &   
      -0.335131147774E+01_ark*xs1** 1*xs2** 1*xst** 5 +  &   
      -0.111612795996E+01_ark*xs1** 1*xs2** 1*xst** 6 +  &   
       0.480892471688E+01_ark*xs1** 1*xs2** 1*xst** 7 +  &   
      -0.315602373516E+01_ark*xs1** 1*xs2** 1*xst** 8 +  &   
       0.664408659597E+00_ark*xs1** 1*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 1*xst**18 +  &   
       0.353903344003E+00_ark*xs1** 1*xs2** 3*xst** 0 +  &   
      -0.103869582792E+01_ark*xs1** 1*xs2** 3*xst** 1 +  &   
       0.100357404164E+01_ark*xs1** 1*xs2** 3*xst** 2 +  &   
      -0.305465827499E+00_ark*xs1** 1*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 3*xst**18 +  &   
      -0.346916593681E-01_ark*xs1** 1*xs2** 5*xst** 0 +  &   
      -0.810393741779E-01_ark*xs1** 1*xs2** 5*xst** 1 +  &   
       0.378281880477E+00_ark*xs1** 1*xs2** 5*xst** 2 +  &   
      -0.373298806455E+00_ark*xs1** 1*xs2** 5*xst** 3 +  &   
       0.111383875645E+00_ark*xs1** 1*xs2** 5*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 5*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 7*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 1*xs2** 9*xst**18 +  &   
      -0.107602847167E+01_ark*xs1** 2*xs2** 1*xst** 0 +  &   
       0.410853015153E+01_ark*xs1** 2*xs2** 1*xst** 1 +  &   
      -0.401029772799E+01_ark*xs1** 2*xs2** 1*xst** 2 +  &   
      -0.121517945186E+01_ark*xs1** 2*xs2** 1*xst** 3 +  &   
       0.258379864700E+01_ark*xs1** 2*xs2** 1*xst** 4 +  &   
       0.160277350147E+01_ark*xs1** 2*xs2** 1*xst** 5 +  &   
      -0.249276779140E+01_ark*xs1** 2*xs2** 1*xst** 6 +  &   
       0.673604774339E+00_ark*xs1** 2*xs2** 1*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 1*xst**18 +  &   
       0.483396368056E+00_ark*xs1** 2*xs2** 3*xst** 0 +  &   
      -0.239992884318E+01_ark*xs1** 2*xs2** 3*xst** 1 +  &   
       0.259594110620E+01_ark*xs1** 2*xs2** 3*xst** 2 +  &   
       0.193353870186E+01_ark*xs1** 2*xs2** 3*xst** 3 +  &   
      -0.283042612989E+01_ark*xs1** 2*xs2** 3*xst** 4 +  &   
      -0.296551288827E+01_ark*xs1** 2*xs2** 3*xst** 5 +  &   
       0.457399095125E+01_ark*xs1** 2*xs2** 3*xst** 6 +  &   
      -0.141766701938E+01_ark*xs1** 2*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 3*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 5*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 7*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 2*xs2** 9*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2** 1 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2** 2 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2** 3 +  &   
      -0.457936836965E-01_ark*xs1** 3*xs2** 1*xs2** 4 +  &   
       0.206729030729E-01_ark*xs1** 3*xs2** 1*xs2** 5 +  &   
       0.783058621632E-01_ark*xs1** 3*xs2** 1*xs2** 6 +  &   
      -0.412587437113E-01_ark*xs1** 3*xs2** 1*xs2** 7 +  &   
      -0.138555397864E-01_ark*xs1** 3*xs2** 1*xs2** 8 +  &   
       0.131859699808E-01_ark*xs1** 3*xs2** 1*xs2** 9 +  &   
      -0.228974684356E-02_ark*xs1** 3*xs2** 1*xs2**10 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2**11 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2**12 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2**13 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2**14 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2**15 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2**16 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2**17 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 1*xs2**18      
muY_2=  0.000000000000E+00_ark*xs1** 3*xs2** 3*xs2**0 +  &
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 3*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 5*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 7*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 3*xs2** 9*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst** 2 +  &   
      -0.704432161254E+00_ark*xs1** 4*xs2** 1*xst** 3 +  &   
       0.686877393192E+00_ark*xs1** 4*xs2** 1*xst** 4 +  &   
       0.514499990474E-01_ark*xs1** 4*xs2** 1*xst** 5 +  &   
      -0.935490678324E-01_ark*xs1** 4*xs2** 1*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 1*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 3*xst**10 +  &   
       0.927233073942E+00_ark*xs1** 4*xs2** 3*xst**11 +  &   
      -0.391738675843E+00_ark*xs1** 4*xs2** 3*xst**12 +  &   
      -0.978373295350E+00_ark*xs1** 4*xs2** 3*xst**13 +  &   
      -0.471757787232E+00_ark*xs1** 4*xs2** 3*xst**14 +  &   
       0.782559125569E+00_ark*xs1** 4*xs2** 3*xst**15 +  &   
       0.110315175361E+01_ark*xs1** 4*xs2** 3*xst**16 +  &   
      -0.125985374824E+01_ark*xs1** 4*xs2** 3*xst**17 +  &   
       0.319095000991E+00_ark*xs1** 4*xs2** 3*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst** 4 +  &   
      -0.130028319922E+01_ark*xs1** 4*xs2** 5*xst** 5 +  &   
       0.257711587708E+01_ark*xs1** 4*xs2** 5*xst** 6 +  &   
      -0.172559831008E+01_ark*xs1** 4*xs2** 5*xst** 7 +  &   
       0.388782058836E+00_ark*xs1** 4*xs2** 5*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 5*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 7*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 4*xs2** 9*xst**18 +  &   
      -0.553803841029E-01_ark*xs1** 5*xs2** 1*xst** 0 +  &   
      -0.645855979808E+00_ark*xs1** 5*xs2** 1*xst** 1 +  &   
       0.741561375126E+00_ark*xs1** 5*xs2** 1*xst** 2 +  &   
       0.146438064036E+01_ark*xs1** 5*xs2** 1*xst** 3 +  &   
       0.447305334912E+00_ark*xs1** 5*xs2** 1*xst** 4 +  &   
      -0.111039860408E+01_ark*xs1** 5*xs2** 1*xst** 5 +  &   
      -0.143936660007E+01_ark*xs1** 5*xs2** 1*xst** 6 +  &   
       0.782018974967E-01_ark*xs1** 5*xs2** 1*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst** 8 +  &   
       0.152752666154E+01_ark*xs1** 5*xs2** 1*xst** 9 +  &   
      -0.717187838369E+00_ark*xs1** 5*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 1*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst** 0 +  &   
       0.569482779272E-01_ark*xs1** 5*xs2** 3*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst** 8 +  &   
      -0.358315655015E-01_ark*xs1** 5*xs2** 3*xst** 9 +  &   
       0.480953858956E-01_ark*xs1** 5*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 3*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 5*xst**18 +  &   
      -0.283004954425E-01_ark*xs1** 5*xs2** 7*xst** 0 +  &   
       0.642230604849E-01_ark*xs1** 5*xs2** 7*xst** 1 +  &   
      -0.339047768664E-01_ark*xs1** 5*xs2** 7*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 5*xs2** 7*xst**18 +  &   
       0.870024767838E+00_ark*xs1** 6*xs2** 1*xst** 0 +  &   
      -0.151259000995E+01_ark*xs1** 6*xs2** 1*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 1*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 3*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 5*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 6*xs2** 7*xst**18      
muY_3= 0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**0 +  &
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst** 4 +  &   
       0.478706056806E+00_ark*xs1** 7*xs2** 1*xst** 5 +  &   
       0.166166641828E-01_ark*xs1** 7*xs2** 1*xst** 6 +  &   
       0.241450145421E+00_ark*xs1** 7*xs2** 1*xst** 7 +  &   
      -0.378089571567E+00_ark*xs1** 7*xs2** 1*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 1*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xs2**13 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xs2**14 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xs2**15 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xs2**16 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xs2**17 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 3*xs2**18 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 0 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 1 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 2 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 3 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 4 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 5 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 6 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 7 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 8 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2** 9 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2**10 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2**11 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xs2**12 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 7*xs2** 5*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 1*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 3*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 8*xs2** 5*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 1*xst**18 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 0 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 1 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 2 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1** 9*xs2** 3*xst**18 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 0 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 1 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 2 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 3 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 4 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 5 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 6 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 7 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 8 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst** 9 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**10 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**11 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**12 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**13 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**14 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**15 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**16 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**17 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 1*xst**18 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 0 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 1 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 2 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 3 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 4 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 5 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 6 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 7 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 8 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst** 9 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**10 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**11 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**12 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**13 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**14 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**15 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**16 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**17 +  &   
       0.000000000000E+00_ark*xs1**10*xs2** 3*xst**18      
                                                                     
!************************************
muX = muX_0 + muX_1 + muX_2 + muX_3
!************************************

!************************************
muY = muY_1 + muY_2 + muY_3
!************************************

muX =(PI-theta)*muX
muY = muY

!************************************

end subroutine DIPSA



end module pot_user

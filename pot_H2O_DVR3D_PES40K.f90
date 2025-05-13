!
!  This unit is for a user defined potential 
!
module pot_user
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xy2

  implicit none

  public MLdipole,MLpoten,ML_MEP,MLpoten_name

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
 ! Check the potential name 
 subroutine MLpoten_name(name)
   !
   character(len=cl),intent(in) ::  name
   character(len=cl),parameter ::  poten_name = 'DVR3D_PES40'
   ! 
   if (poten_name/=trim(name)) then
     write(out,"('a,a,a,a')") 'Wrong Potential ',trim(name),'; should be ',trim(poten_name)
   endif
   !
   write(out,"(a,x,a)") '  Using USER-type PES ',trim(poten_name)
   !
 end subroutine MLpoten_name
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
   real(ark)               :: r12,r32,alpha,xcos,v
   real(ark),parameter     :: tocm = 219474.63067_ark

     r12 = local(1)/bohr ; r32 = local(2)/bohr ;  alpha = local(3)
     !
     !r12 = 0.9586490/bohr
     !r32 = 0.963217334478864291542644713544619/bohr
     !alpha = 3.14159265358979323846264338327950
     ! 
     xcos = cos(alpha)
     !
     call potv(v,r12,r32,xcos,force)
     !
     v = v*tocm
     ! 
     f = v
     !
 end function MLpoten
 !
 !

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

 SUBROUTINE potv(V,R1,R2,xcos,xp)
       
    IMPLICIT real(ark) (A-H,O-Y),logical(z)
    real(ark),intent(in) ::  xp(:)       
    !
    real(ark) :: RZ=0.9576257_ark,RHO=75.48992362_ark
    real(ark) :: TOANG= 0.5291772_ark
    real(ark) ::  X1= 1.0_ark,X0=0.0_ark,TINY=9.0D-15,X2=2.0_ark
    real(ark) :: xpBu(1:300),a,b
    integer(ik),parameter :: npropin = 246
    !
    !character(cl) :: type_calc = 'RADAU'
    !
    character(cl) :: type_calc = 'G1.EQ.X0'
      !
      G1 = 0
      !
      !select case (trim(type_calc))
      !  !
      !case ('G1.EQ.X0') 
      !IF (G1.EQ.X0) THEN
      !   !stop 'G1=0 options is not implemented'
      !   g1 = 0
      !   Q1 = R1
      !   Q2 = R2
      !   Q3 = 0 
      !   THETA = ACOS(XCOS)
      !case ('G2.EQ.X0') 
      !!ELSE IF (G2 .EQ. X0) THEN
      !   stop 'G2=0 options is not implemented'
      !   g2 = 0
      !   !
      !   XX = R1 * G1
      !   YY = R1 * (X1 - G1)
      !   IF (R2 .EQ. X0 .OR. XCOS .GE. (X1 - TINY)) THEN
      !      Q1 = ABS(XX - R2)
      !      Q2 = (YY + R2)
      !      COST = -X1
      !   ELSE IF (XCOS .LE. (TINY - X1)) THEN
      !      Q1 = (XX + R2)
      !      Q2 = ABS(YY + R2)
      !      COST = X1
      !   ELSE
      !      Q1 = SQRT(XX*XX + R2*R2 - X2*XX*R2*XCOS)
      !      Q2 = SQRT(YY*YY + R2*R2 + X2*YY*R2*XCOS)
      !      COST = (Q1**2 + Q2**2 - R1**2) / (X2 * Q1 * Q2)
      !   ENDIF
      !   THETA = ACOS(COST)
      !   !
      !case('RADAU')
      !   !
      !   a = sqrt(molec%AtomMasses(1) / (molec%AtomMasses(1)+molec%AtomMasses(2)+molec%AtomMasses(3)))
      !   b = molec%AtomMasses(2) / (molec%AtomMasses(2)+molec%AtomMasses(3))
      !   g1 = x1 - a / (a+b-a*b)
      !   g2 = x1 - a / (x1-b+a*b)
      !   !         
      !   F1= X1/G1
      !   F2= X1/G2
      !   F12= X1 - F1*F2
      !   P1= R1*(X1-F1)/(G2*F12)
      !   P2= R2*(X1-F2)/(G1*F12)
      !   S1= R1-P1
      !   S2= R2-P2
      !   Q1= SQRT(P1*P1 + S2*S2 + X2*P1*S2*XCOS)/(X1-G1)
      !   Q2= SQRT(P2*P2 + S1*S1 + X2*P2*S1*XCOS)/(X1-G2)
      !   Q3= SQRT(P1*P1 + P2*P2 - X2*P1*P2*XCOS)
      !   COST = (Q1*Q1 + Q2*Q2 - Q3*Q3)/(X2*Q1*Q2)
      !   !
      !   if (COST<-1.0_ark+small_) then 
      !     THETA = pi
      !   else
      !     THETA = ACOS(COST)
      !   endif
      !   !
      !end select 
      !
      Q1 = R1
      Q2 = R2
      Q3 = 0 
      THETA = ACOS(XCOS)

      !open(unit=32,status='old',form='formatted',file='pot.fit39')
      !read(32,*) npropin
      !do i = 1,npropin
      !   read(32,*) i_t, xp(i)
      !end do
      !
      !close(unit=32)
      !
      !if (npropin>size(xp)) then 
      !  write(6,"('wifin: Too many parameters in pot.fit39: ',&
      !                    i8,' max = ',i8)") npropin,size(xp)
      !  stop 'wifin: Too many parameters in pot.fit39'
      !endif
      !
      att1=Q1*0.5291772d0
      att2=Q2*0.5291772d0
      !
      call potvBuPokaz(Vlow,R1,R2,xcos,q1,q2,q3,theta)
      call poten(xp(1:npropin),Vup,att1,att2,THETA)
      Vlow=Vlow*219474.634d0-779.9304d0

      E0 = 35000.d0
      de = Vup - E0
      
      ga0 = 1.d0 /130.d0
      ga1 = 1.d0 / (130.d0**3)

      ga = ga0 + ga1 * (de**2)
      
      ffu = 0.5_ark * (1.0_ark + tanh(ga * de))
      ffl = 0.5_ark * (1.0_ark + tanh(-ga * de))
      Vv  = ffl * Vlow + ffu * Vup
      !
      ! 219474.6313708
      v = Vv / 219474.634d0

  end SUBROUTINE potv
      


  SUBROUTINE potvBuPokaz(V,R1,R2,xcos,q1,q2,q3,theta)

      !
      !COMMON /MASS/ XMASS(3),G1,G2,xmassr(3)
      !
      IMPLICIT real(ark) (A-H,O-Y),logical(z)
      !
      real(ark),intent(in) :: r1,r2,xcos,q1,q2,q3,theta
      real(ark),intent(out):: v
      real(ark) :: vn
      dimension pv(32),u(32)
      dimension ht(50)
      DATA RZ/0.9576257   /,RHO/75.48992362  /
      DATA TOANG/0.5291772/, CMTOAU/219474.624/
      DATA A/2.22600/
      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/
      integer(ik) :: i

      DATA pv/0.999963011541667,-0.000831585900216,0.000440024019720,&
             -0.000790390661341,0.000682751133154,-0.000787433811192,&
             0.000137727815420,-0.000051565403109,-0.000122703443403,&
             0.003387614188353,-0.000271677251731,-0.004052838735161,&
             0.001707147193552,-0.001805888500169,-0.006788069703931,&
             0.001208782633980,0.001108778944947,0.009632702226484,&
            -0.009093431064531,-0.000836179586411,0.001775570609778,&
             0.000051292593874,0.006623947296453,-0.001827850224422,&
            -0.012626322010510,-0.003304654368357,0.010437389017074,&
             0.000227788649284,0.008030675091008,&
            -0.000723724819511,0.002063490878186,0.001356187797029/

!        Equilibrium


      THETAeq = (180.d0 - RHO)*3.141592654/180.d0
      Req= RZ/TOANG
      Y1 = Q1 - Req
      Y2 = Cos(THETA) - Cos(THETAeq)
      Y3 = Q2 - Req

      S1 = (Y1 + Y3)/2.0d0
      S2 =  Y2
      S3 = (Y1 - Y3)/2.0d0 
!   Potential  symmetrical in S1,S2,S3. Full set of
!   derivatives for each order (1-4) and upper terms for others.

      HT(1)=1.0d0
      ht(2)=s1*s2
      ht(3)=s1**2
      ht(4)=s1**3 
      ht(5)=s3**2
      ht(6)=s1**4
      ht(7)=s2**2
      ht(8)=s2**3
      ht(9)=s2**4
      ht(10)=s1**2*s2
      ht(11)=s2**2*s1
      ht(12)=s3**2*s1
      ht(13)=s3**2*s2
      ht(14)=s3**4
      ht(15)=s1**3*S2
      ht(16)=S2**3*s1
      ht(17)=s1**2*s2**2
      ht(18)=s1**2*s3**2
      ht(19)=s1*s2*s3**2
      ht(20)=s2**2*s3**2
      ht(21)=s1**5
      ht(22)=s2**5
      ht(23)=s1**4*s2
      ht(24)=s1**3*s2**2
      ht(25)=s1**3*s3**2
      ht(26)=s1**2*s2**3
      ht(27)=s1**2*s2*s3**2
      ht(28)=s1*s2**4      
      ht(29)=s1*s3**4      
      ht(30)=s2*s3**4
      ht(31)=s1*s2**2*s3**2
      ht(32)=s2**3*s3**2


      vt = 0.0d0
      do i = 1,32
         u(i) = (ht(i) * pv(i))
         vt   = vt + u(i)
      enddo

      call pantsPokaz(vn,r1,r2,xcos,q1,q2,q3,theta,vbn)

! This is where the ab initio and empirical potentials are defined
! vn: ab initio
! vbn: empirically determined barrier to linearity
! vt: empirical morphing of the ab initio PES

      v = vn*vt
      end SUBROUTINE potvBuPokaz

      SUBROUTINE pantsPokaz(V,R1,R2,xcos,q1_,q2_,q3_,theta,vbn)

      implicit real(ark) (a-h,o-y),logical(z)
      !
      real(ark),intent(in) :: r1,r2,xcos,q1_,q2_,q3_,theta
      real(ark),intent(out):: v,vbn
      !
      dimension hp(84),hc(84),pp(84),uz(84)
      integer(ik) :: jp(84)
      dimension rij(3)

      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/

  
      q1 = q1_
      q2 = q2_
      !
      vbn=1.0
      if(q1.gt.4.0) q1=4.0
      if(q2.gt.4.0) q2=4.0
      !
      CALL BREITB3lin(Vbr,Q1,Q2,THETA)
      CALL PESleq6(VQ,Q1,Q2,THETA)
      CALL pesd2x(Vd,Q1,Q2,THETA)
      CALL bodc16(V1,Q1,Q2,THETA)
      !
!     rij -- bohrs, THETA -- rad
      rij(1) = Q1
      rij(2) = Q2
      rij(3) = THETA
      !
      Call rel(Vr,rij(1),rij(2),rij(3))	
      !
	  call pots(vp,rij(1),rij(2),rij(3))
      !
!     rij -- angstrems
      rij(1) = Q1*0.5291772
      rij(2) = Q2*0.5291772
      rij(3) = THETA
      !
      Call wifinBuPokaz(Vai2,rij(1),rij(2),rij(3))
      !
!     v -- hartries, vai2 -- 1/cm  | should be 219474.6313708 
      v=vai2/219474.624+vp+vr+vbr+vd+v1+vq
      !
      RETURN
      END SUBROUTINE pantsPokaz

      SUBROUTINE wifinBuPokaz(v,r1,r2,th)

      implicit real(ark) (a-h,o-z)

!
! In here almost all the points have been included, 342 in total
!
!
      

           reoh=0.958649d0
           thetae=104.3475d0
           b1=2.0d0
           roh=0.951961d0
           alphaoh=2.587949757553683d0
           phh2=6.70164303995d0
           t0=0.01d0
           ut=20.d0
           x0=2.5d0
           ut2=20.d0
           x02=2.5d0
           
         thetae=thetae*.314159265358979312d01*.00555555555555555555d0

         xs1=(r1+r2)*0.5d0-reoh
         xs2=(r1-r2)*0.5d0
         xst=cos(th)-cos(thetae)

         rs=sqrt(0.5d0)*(r1+r2)
         rm=sqrt(0.5d0)*(r1-r2)

         rr1=r1-roh
         rr2=r2-roh

         xep1=exp(-2.d0*alphaoh*rr1)-2.d0*exp(-alphaoh*rr1)+1.d0
         xep2=exp(-2.d0*alphaoh*rr2)-2.d0*exp(-alphaoh*rr2)+1.d0
         xep3=exp(-b1*(rr1**2+rr2**2))
         rhh=sqrt(r1**2+r2**2-2.d0*r1*r2*cos(th))
         vhh=0.900642240911285975d6*exp(-phh2*rhh)
         vpb1=0.518622556959170834d5*xep1
         vpb2=0.518622556959170834d5*xep2

      v0=-0.69833729793860E+02*xs1**0*xs2**0*xst**0
      vp1=-0.88931405097186E+02*xs1**0*xs2**0*xst**1&
     -0.87065167087729E+04*xs1**1*xs2**0*xst**0     &
     +0.18561833211448E+05*xs1**0*xs2**0*xst**2     &
     -0.22925714045631E+06*xs1**0*xs2**2*xst**0     &
     -0.25473078391522E+05*xs1**1*xs2**0*xst**1     &
     -0.24221902469321E+06*xs1**2*xs2**0*xst**0     &
     +0.91104690191314E+03*xs1**0*xs2**0*xst**3     &
     -0.20846258713362E+05*xs1**0*xs2**2*xst**1     &
     -0.97259666765996E+04*xs1**1*xs2**0*xst**2     &
     +0.21895229646374E+07*xs1**1*xs2**2*xst**0     &
     +0.23321398548491E+05*xs1**2*xs2**0*xst**1     &
     +0.70135119326976E+06*xs1**3*xs2**0*xst**0     &
     +0.35115546540196E+04*xs1**0*xs2**0*xst**4     &
     +0.56216397977173E+05*xs1**0*xs2**2*xst**2     &
     -0.20038134529604E+07*xs1**0*xs2**4*xst**0     &
     -0.10382409164897E+05*xs1**1*xs2**0*xst**3     &
     -0.15791163967577E+05*xs1**1*xs2**2*xst**1     &
     +0.69131950382604E+05*xs1**2*xs2**0*xst**2     &
     -0.83355495109342E+07*xs1**2*xs2**2*xst**0     &
     -0.72835013274781E+05*xs1**3*xs2**0*xst**1     &
     -0.20429129369404E+07*xs1**4*xs2**0*xst**0     &
     +0.55792697169003E+03*xs1**0*xs2**0*xst**5     &
     +0.24482126719682E+04*xs1**0*xs2**2*xst**3     &
     -0.27399255063928E+05*xs1**0*xs2**4*xst**1     &
     +0.30011447505506E+04*xs1**1*xs2**0*xst**4     &
     -0.40810705479272E+05*xs1**1*xs2**2*xst**2     &
     +0.13148302968514E+08*xs1**1*xs2**4*xst**0     &
     -0.59362784536737E+04*xs1**2*xs2**0*xst**3     &
     -0.86061422326869E+05*xs1**2*xs2**2*xst**1     &
     -0.32156522858267E+05*xs1**3*xs2**0*xst**2     &
     +0.20293315578616E+08*xs1**3*xs2**2*xst**0     &
     +0.19710424146155E+05*xs1**4*xs2**0*xst**1     &
     +0.37828666711928E+07*xs1**5*xs2**0*xst**0     &
     +0.53407786571380E+03*xs1**0*xs2**0*xst**6     &
     +0.24159882636341E+05*xs1**0*xs2**2*xst**4     &
     +0.12125087320692E+06*xs1**0*xs2**4*xst**2     &
     -0.63007628874696E+07*xs1**0*xs2**6*xst**0     &
     -0.30411762363758E+03*xs1**1*xs2**0*xst**5     &
     -0.78383097659001E+05*xs1**1*xs2**2*xst**3     &
     +0.58781557929045E+05*xs1**1*xs2**4*xst**1     &
     +0.13375434109158E+03*xs1**2*xs2**0*xst**4     &
     +0.33893910748383E+06*xs1**2*xs2**2*xst**2     &
     -0.37826889997246E+08*xs1**2*xs2**4*xst**0     &
     +0.62361486440787E+04*xs1**3*xs2**0*xst**3     &
     -0.10596388642131E+06*xs1**3*xs2**2*xst**1     &
     +0.28702790231058E+05*xs1**4*xs2**0*xst**2     &
     -0.37680124674448E+08*xs1**4*xs2**2*xst**0     &
     +0.76453502516977E+05*xs1**5*xs2**0*xst**1     &
     -0.65507833178799E+07*xs1**6*xs2**0*xst**0     &
     -0.19940650191011E+04*xs1**0*xs2**0*xst**7     &
     -0.53454319818075E+05*xs1**0*xs2**2*xst**5     &
     -0.18482082317071E+06*xs1**0*xs2**4*xst**3     &
     -0.32482357937925E+06*xs1**0*xs2**6*xst**1     &
     +0.14750414212517E+05*xs1**1*xs2**0*xst**6     &
     +0.10992166354556E+06*xs1**1*xs2**2*xst**4     &
     -0.42127450223572E+06*xs1**1*xs2**4*xst**2     &
     +0.31697813513995E+08*xs1**1*xs2**6*xst**0     &
     -0.76319535999048E+05*xs1**2*xs2**0*xst**5     &
     -0.39996355835188E+04*xs1**2*xs2**2*xst**3     &
     -0.10848426889348E+06*xs1**2*xs2**4*xst**1     &
     +0.16031984625558E+06*xs1**3*xs2**0*xst**4     &
     -0.99440969626463E+06*xs1**3*xs2**2*xst**2     &
     +0.68374298978030E+08*xs1**3*xs2**4*xst**0     &
     -0.71454371570541E+05*xs1**4*xs2**0*xst**3     &
     +0.18244007909156E+07*xs1**4*xs2**2*xst**1     &
     +0.12228399205468E+06*xs1**5*xs2**0*xst**2     &
     +0.50353404665938E+08*xs1**5*xs2**2*xst**0     &
     -0.20016358839030E+06*xs1**6*xs2**0*xst**1     &
     +0.90962316465138E+07*xs1**7*xs2**0*xst**0     &
     +0.10149096018617E+04*xs1**0*xs2**0*xst**8     &
     +0.32566540304056E+05*xs1**0*xs2**2*xst**6     &
     +0.17026834226150E+05*xs1**0*xs2**4*xst**4     &
     -0.20338737696627E+06*xs1**0*xs2**6*xst**2     &
     -0.11419910468431E+08*xs1**0*xs2**8*xst**0     &
     -0.20146404554912E+05*xs1**1*xs2**0*xst**7     &
     +0.20490886901277E+05*xs1**1*xs2**2*xst**5     &
     +0.11627387113119E+07*xs1**1*xs2**4*xst**3     &
     -0.26762513146301E+07*xs1**1*xs2**6*xst**1     &
     +0.63330899794849E+05*xs1**2*xs2**0*xst**6     &
     -0.29823655377502E+06*xs1**2*xs2**2*xst**4     &
     -0.11572269239294E+04*xs1**2*xs2**4*xst**2     &
     -0.53010926363911E+08*xs1**2*xs2**6*xst**0     &
     -0.58084107449611E+05*xs1**3*xs2**0*xst**5     &
     +0.10772705558304E+07*xs1**3*xs2**2*xst**3     &
     -0.47494756934226E+05*xs1**3*xs2**4*xst**1     &
     +0.21811000394523E+05*xs1**4*xs2**0*xst**4     &
     +0.14998583534579E+07*xs1**4*xs2**2*xst**2     &
     -0.64501462536687E+08*xs1**4*xs2**4*xst**0     &
     -0.14570114276761E+06*xs1**5*xs2**0*xst**3     &
     -0.66198424020631E+05*xs1**5*xs2**2*xst**1     &
     +0.32017806514222E+06*xs1**6*xs2**0*xst**2      
     vp2=-0.46302044457094E+08*xs1**6*xs2**2*xst**0 &
     +0.13193687870007E+06*xs1**7*xs2**0*xst**1     &
     -0.11157788556702E+08*xs1**8*xs2**0*xst**0     &
     +0.20805390750196E+04*xs1**0*xs2**0*xst**9     &
     +0.10977717065665E+06*xs1**0*xs2**2*xst**7     &
     +0.40544491766050E+06*xs1**0*xs2**4*xst**5     &
     -0.24296852330663E+05*xs1**0*xs2**6*xst**3     &
     +0.38177105134606E+07*xs1**0*xs2**8*xst**1     &
     -0.16390212024836E+05*xs1**1*xs2**0*xst**8     &
     -0.35305243614455E+06*xs1**1*xs2**2*xst**6     &
     +0.70273286813606E+06*xs1**1*xs2**4*xst**4     &
     +0.16135959231803E+05*xs1**1*xs2**6*xst**2     &
     +0.76303520865603E+04*xs1**1*xs2**8*xst**0     &
     +0.14575708768624E+06*xs1**2*xs2**0*xst**7     &
     +0.85864519904820E+06*xs1**2*xs2**2*xst**5     &
     -0.53581284988223E+07*xs1**2*xs2**4*xst**3     &
     +0.12582250565590E+08*xs1**2*xs2**6*xst**1     &
     -0.37156076835264E+06*xs1**3*xs2**0*xst**6     &
     -0.47107091715944E+06*xs1**3*xs2**2*xst**4     &
     +0.10422174650334E+08*xs1**3*xs2**4*xst**2     &
     -0.73996618932697E+08*xs1**3*xs2**6*xst**0     &
     +0.25959907471094E+06*xs1**4*xs2**0*xst**5     &
     -0.34057704907398E+07*xs1**4*xs2**2*xst**3     &
     -0.17848038023478E+05*xs1**4*xs2**4*xst**1     &
     -0.29021508562549E+06*xs1**5*xs2**0*xst**4     &
     +0.10456945081877E+08*xs1**5*xs2**2*xst**2     &
     -0.71012605744803E+08*xs1**5*xs2**4*xst**0     &
     -0.64745699046240E+06*xs1**6*xs2**0*xst**3     &
     -0.34346142250957E+08*xs1**6*xs2**2*xst**1     &
     +0.12578165297978E+06*xs1**7*xs2**0*xst**2     &
     +0.18473904269184E+08*xs1**7*xs2**2*xst**0     &
     -0.35458625565554E+06*xs1**8*xs2**0*xst**1     &
     +0.97231092585162E+07*xs1**9*xs2**0*xst**0     &
     -0.97725063032968E+03*xs1**0*xs2**0*xst**10    &
     -0.58810552328534E+05*xs1**0*xs2**2*xst**8     &
     -0.15097015389582E+06*xs1**0*xs2**4*xst**6     &
     +0.28619491235330E+05*xs1**0*xs2**6*xst**4     &
     +0.64423406036478E+07*xs1**0*xs2**8*xst**2     &
     +0.13304147310031E+05*xs1**0*xs2**10*xst**0    &
     +0.24818228591588E+05*xs1**1*xs2**0*xst**9     &
     +0.82145427857098E+05*xs1**1*xs2**2*xst**7     &
     -0.28560408228532E+07*xs1**1*xs2**4*xst**5     &
     +0.17758056545809E+07*xs1**1*xs2**6*xst**3     &
     +0.47618501940831E+08*xs1**1*xs2**8*xst**1     &
     -0.14362715517470E+06*xs1**2*xs2**0*xst**8     &
     +0.87939534136632E+06*xs1**2*xs2**2*xst**6     &
     +0.43542816466615E+05*xs1**2*xs2**4*xst**4     &
     +0.40855419184907E+08*xs1**2*xs2**6*xst**2     &
     +0.15590054335119E+09*xs1**2*xs2**8*xst**0     &
     +0.23847828237922E+06*xs1**3*xs2**0*xst**7     &
     -0.30988931738056E+07*xs1**3*xs2**2*xst**5     &
     +0.81912153955218E+07*xs1**3*xs2**4*xst**3     &
     -0.14524552518400E+09*xs1**3*xs2**6*xst**1     &
     -0.27542471524966E+06*xs1**4*xs2**0*xst**6     &
     +0.27085838951528E+07*xs1**4*xs2**2*xst**4     &
     -0.48766162840619E+08*xs1**4*xs2**4*xst**2     &
     +0.64814050166329E+09*xs1**4*xs2**6*xst**0     &
     +0.10272921778939E+07*xs1**5*xs2**0*xst**5     &
     -0.61301714047235E+07*xs1**5*xs2**2*xst**3     &
     +0.48355468765097E+07*xs1**5*xs2**4*xst**1     &
     -0.16950012608773E+07*xs1**6*xs2**0*xst**4     &
     -0.36875319697041E+08*xs1**6*xs2**2*xst**2     &
     +0.32375440102365E+09*xs1**6*xs2**4*xst**0     &
     +0.20890611536445E+07*xs1**7*xs2**0*xst**3     &
     +0.85634774262523E+08*xs1**7*xs2**2*xst**1     &
     -0.28453560365440E+07*xs1**8*xs2**0*xst**2     &
     -0.70110715851666E+05*xs1**8*xs2**2*xst**0     &
     +0.16593108812323E+06*xs1**9*xs2**0*xst**1     &
     -0.46099243127663E+07*xs1**10*xs2**0*xst**0    &
     -0.12157030586644E+04*xs1**0*xs2**0*xst**11    &
     -0.82857116167377E+05*xs1**0*xs2**2*xst**9     &
     -0.58451292296010E+06*xs1**0*xs2**4*xst**7     &
     +0.11599886074607E+07*xs1**0*xs2**6*xst**5     &
     -0.81818449411873E+07*xs1**0*xs2**8*xst**3     &
     -0.38544362571641E+08*xs1**0*xs2**10*xst**1    &
     +0.13772894688675E+05*xs1**1*xs2**0*xst**10    &
     +0.45365215594235E+06*xs1**1*xs2**2*xst**8     &
     +0.79988903784077E+06*xs1**1*xs2**4*xst**6     &
     +0.35926515007257E+05*xs1**1*xs2**6*xst**4     &
     -0.92442595978415E+08*xs1**1*xs2**8*xst**2     &
     +0.43271322021484E+08*xs1**1*xs2**10*xst**0    &
     -0.17400692445498E+06*xs1**2*xs2**0*xst**9      
     vp3=-0.14268059220845E+07*xs1**2*xs2**2*xst**7 &
     +0.40545378149964E+07*xs1**2*xs2**4*xst**5     &
     +0.17585280543815E+05*xs1**2*xs2**6*xst**3     &
     -0.10276063972665E+09*xs1**2*xs2**8*xst**1     &
     +0.56807964581015E+06*xs1**3*xs2**0*xst**8     &
     +0.18437826095465E+07*xs1**3*xs2**2*xst**6     &
     -0.12864944253983E+08*xs1**3*xs2**4*xst**4     &
     -0.10451574046722E+09*xs1**3*xs2**6*xst**2     &
     -0.43651678669326E+09*xs1**3*xs2**8*xst**0     &
     -0.73083658333429E+06*xs1**4*xs2**0*xst**7     &
     -0.17056814638312E+05*xs1**4*xs2**2*xst**5     &
     +0.59610184323397E+08*xs1**4*xs2**4*xst**3     &
     +0.49938980132233E+09*xs1**4*xs2**6*xst**1     &
     +0.68167414002953E+06*xs1**5*xs2**0*xst**6     &
     -0.59994275353066E+07*xs1**5*xs2**2*xst**4     &
     +0.58536371084931E+04*xs1**5*xs2**4*xst**2     &
     -0.14136428052005E+10*xs1**5*xs2**6*xst**0     &
     +0.95044891506254E+05*xs1**6*xs2**0*xst**5     &
     +0.50803779370631E+08*xs1**6*xs2**2*xst**3     &
     -0.17596700209434E+05*xs1**6*xs2**4*xst**1     &
     +0.24656515735182E+07*xs1**7*xs2**0*xst**4     &
     +0.27067874384387E+08*xs1**7*xs2**2*xst**2     &
     -0.30519783492194E+09*xs1**7*xs2**4*xst**0     &
     +0.99859651431732E+05*xs1**8*xs2**0*xst**3     &
     -0.58311937002508E+08*xs1**8*xs2**2*xst**1     &
     +0.31058183048947E+07*xs1**9*xs2**0*xst**2     &
     +0.33566859753505E+05*xs1**9*xs2**2*xst**0     &
     +0.28534012133708E+03*xs1**0*xs2**0*xst**12    &
     +0.73181259963006E+05*xs1**0*xs2**2*xst**10    &
     +0.30712308016695E+06*xs1**0*xs2**4*xst**8     &
     -0.38988311746447E+06*xs1**0*xs2**6*xst**6     &
     +0.65071010834709E+04*xs1**0*xs2**8*xst**4     &
     +0.33025369570808E+08*xs1**0*xs2**10*xst**2    &
     -0.17164906136195E+08*xs1**0*xs2**12*xst**0    &
     -0.13913272952356E+05*xs1**1*xs2**0*xst**11    &
     -0.33399040011605E+06*xs1**1*xs2**2*xst**9     &
     +0.10266852525562E+07*xs1**1*xs2**4*xst**7     &
     -0.39836520885698E+07*xs1**1*xs2**6*xst**5     &
     +0.34049616900886E+08*xs1**1*xs2**8*xst**3     &
     -0.21561729380190E+04*xs1**1*xs2**10*xst**1    &
     +0.14552323333248E+06*xs1**2*xs2**0*xst**10    &
     +0.59410384696083E+06*xs1**2*xs2**2*xst**8     &
     -0.23886244165805E+07*xs1**2*xs2**4*xst**6     &
     +0.50371572986824E+07*xs1**2*xs2**6*xst**4     &
     +0.96234635335528E+08*xs1**2*xs2**8*xst**2     &
     -0.88806690050628E+08*xs1**2*xs2**10*xst**0    &
     -0.43369177090809E+06*xs1**3*xs2**0*xst**9     &
     +0.33593315021506E+05*xs1**3*xs2**2*xst**7     &
     +0.88172686594762E+05*xs1**3*xs2**4*xst**5     &
     -0.44813807185195E+08*xs1**3*xs2**6*xst**3     &
     +0.11963770299610E+09*xs1**3*xs2**8*xst**1     &
     +0.59120204513323E+06*xs1**4*xs2**0*xst**8     &
     +0.12298199958256E+06*xs1**4*xs2**2*xst**6     &
     +0.18693214990721E+08*xs1**4*xs2**4*xst**4     &
     +0.13670141790264E+09*xs1**4*xs2**6*xst**2     &
     +0.41735524210803E+09*xs1**4*xs2**8*xst**0     &
     -0.11245306922000E+07*xs1**5*xs2**0*xst**7     &
     +0.17460051518372E+07*xs1**5*xs2**2*xst**5     &
     -0.11089446429012E+09*xs1**5*xs2**4*xst**3     &
     -0.58155496949524E+09*xs1**5*xs2**6*xst**1     &
     +0.14312315939245E+07*xs1**6*xs2**0*xst**6     &
     +0.26916166982085E+06*xs1**6*xs2**2*xst**4     &
     +0.10382617150293E+09*xs1**6*xs2**4*xst**2     &
     +0.10618975583625E+10*xs1**6*xs2**6*xst**0     &
     -0.34896597462504E+07*xs1**7*xs2**0*xst**5     &
     -0.50144885137780E+08*xs1**7*xs2**2*xst**3
        vp=vp1+vp2+vp3

         vps1=42395.535333d0*xep1
         vps2=42395.535333d0*xep2

         y1=1.d0/(1.d0+exp(ut*(x0-r1)))
         y2=1.d0/(1.d0+exp(ut*(x0-r2)))
         y12=1.d0/(1.d0+exp(ut2*(x02-r1)))
         y22=1.d0/(1.d0+exp(ut2*(x02-r2)))

         vp=vp*xep3*(1-y12)*(1-y22)
         voh1=vpb1*(1-y1)+y1*vps1
         voh2=vpb2*(1-y2)+y2*vps2

        v=v0+vp+voh1+voh2+vhh

        return
        end SUBROUTINE wifinBuPokaz


      
      
      SUBROUTINE BODC16(V,R1,R2,THETA)
      IMPLICIT real(ark) (A-H,O-Z)
      DIMENSION FT(69)

      real(ark),parameter :: C0 = 0.0_ark,SCALE = 1.0e-6_ark
      integer(ik),parameter :: nv = 69
      real(ark),parameter :: ZERO = 0.0_ark
      real(ark),parameter :: RZ = 0.9576257_ark,RHO=75.48992362_ark
      real(ark),parameter :: TOANG = 0.5291772_ark
      real(ark),parameter :: cv(69)  = (/2785.4109260806_ark,-90.0463360612_ark,48.5226385100_ark,&
      201.6994168988_ark,-19.4625405248_ark,217.8639531832_ark,&
      -48.9272565544_ark,-76.0803927294_ark,25.9675092231_ark,-65.2206704082_ark,&
      60.7284555026_ark,-386.5173198637_ark,-2.0352745279_ark,3.5405842322_ark,&
      -32.1830475213_ark,126.8885218953_ark,60.2606780157_ark,134.6300248947_ark,&
      -14.0991224631_ark,560.8792812623_ark,27.8464660146_ark,-67.1464124186_ark,&
      -181.5695164733_ark,43.7287769150_ark,-150.7468776877_ark,-103.2652779984_ark,&
      -159.6432245839_ark,-24.8956181482_ark,185.4825335893_ark,-231.4775497546_ark,&
      -51.9208759775_ark,3.7066140140_ark,-212.4422129941_ark,207.1884490703_ark,&
      383.1659239100_ark,50.8641728660_ark,35.1754939408_ark,127.2280510484_ark,&
      -154.4544312699_ark,55.5787967758_ark,282.6216945851_ark,116.5606405651_ark,&
      5.4433254144_ark,-107.1461094167_ark,-173.8895717556_ark,-26.5859674990_ark,&
      -560.8756697840_ark,-237.7109212157_ark,143.9462552048_ark,-592.3478209334_ark,&
      0.0_ark,-198.6468835005_ark,-19.9674473372_ark,-14.1731270087_ark,193.4510720304_ark,&
      4.6347021028_ark,32.9502486772_ark,-221.1685318724_ark,26.4090449111_ark,&
      -268.432837385_ark,-147.1422366151_ark,133.5465868568_ark,363.9554096142_ark,&
      673.9006484856_ark,214.9454642643_ark,40.7735822438_ark,65.2246188257_ark,&
      173.0708970426_ark,1.9795259929_ark/)
      integer(ik) :: i5,i

      THETAeq = (180._ark - RHO)*3.141592654_ark/180._ark
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

      real(ark),parameter :: C0=0.0_ark,SCALE=1.0e-6_ark
      real(ark),parameter ::  ZERO = 0.0_ark
      real(ark),parameter ::  RZ = 0.9576257_ark,RHO = 75.48992362_ark
      real(ark),parameter ::  TOANG =0.5291772_ark
      real(ark),parameter ::  cv(nv) = (/2520.8349350538_ark,-91.1209324309_ark,49.9198880161_ark,&
       201.5986877584_ark,-27.0761946057_ark,213.0632778654_ark,-42.3461775982_ark,&
       -89.5081758955_ark,3.1590606555_ark,-82.9610483797_ark,91.2462747271_ark,&
       -345.6845385214_ark,2.9125633339_ark,-9.2743637228_ark,19.1249589270_ark,&
       126.8439173732_ark,102.4748093797_ark,29.1206937761_ark,-0.4903235370_ark,&
       548.5647508851_ark,57.4425043137_ark,-65.3115364766_ark,-32.7831814238_ark,&
       111.9709105628_ark,-278.4048997815_ark,-231.4930911427_ark,-323.4608284981_ark,&
       -98.0602324965_ark,181.1676886771_ark,-282.2900662875_ark,167.5274873935_ark,&
       47.8266275957_ark,-213.7402930366_ark,20.8555467237_ark,169.6474418937_ark,&
       -62.1263131255_ark,71.9040085206_ark,72.1511749317_ark,123.9432784777_ark,&
       173.6777915640_ark,58.5107154936_ark,-27.8257872910_ark,99.6002434065_ark,&
       -135.2045098416_ark,145.2475743775_ark,-23.3633794048_ark,-599.5375919902_ark,&
       -540.4989859645_ark,275.5466454040_ark,-502.8009780288_ark,167.2307803194_ark,&
       -109.3097362364_ark,-77.1423552840_ark,492.0803282409_ark,313.6574722992_ark,&
       -22.4885664677_ark,141.5175472800_ark,-1134.9499411146_ark,88.6368503668_ark,&
       -587.9451428738_ark,-247.5329792460_ark,-4.9333140627_ark,322.2004639448_ark,&
       916.1768044740_ark,479.0204491692_ark,93.9437933859_ark,107.5055433254_ark,&
       0.0_ark,63.1640462956_ark,191.6712842623_ark,1.6032078252_ark,78.8696064123_ark,&
       1.0174252672_ark,-377.6585243760_ark,-192.1212462644_ark,653.6712620983_ark,&
       0.0_ark,261.5131272695_ark,146.4973633837_ark,0.0_ark,-137.1867814770_ark,&
       -98.7015480745_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,-97.1438825386_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,&
       0.0_ark,22.7027546667_ark,0.0_ark,0.0_ark,0.0_ark,-341.9933632553_ark,0.0_ark,0.0_ark,&
       -111.8520186501_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,-48.2716070600_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,&
       -5.0447226741_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,144.1528855579_ark/)

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
      real(ark),parameter ::  RZ = 0.9576257_ark, RHO = 75.48992362_ark
      real(ark),parameter ::  TOANG = 0.5291772_ark
      real(ark),parameter ::  cv(nv) =(/0.0_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,&
      0.0_ark,0.0_ark,0.0_ark,0.0_ark,0.0_ark,-0.0519929240_ark,&
      -0.0000574005_ark,-0.0001860844_ark,-0.0006152872_ark,0.0000879624_ark,&
      -0.0004932809_ark,-0.0000565213_ark,0.0006886853_ark,-0.0001027470_ark,&
      0.0003081795_ark,0.0001078370_ark,0.0019925277_ark,0.0001286772_ark,&
      -0.0004511478_ark,0.0000072157_ark,-0.0004439080_ark,-0.0001416903_ark,&
      0.0000984644_ark,-0.0003092456_ark,-0.0029118771_ark,-0.0000562847_ark,&
      0.0000377062_ark,0.0002610717_ark,0.0000367337_ark,0.0_ark,-0.0002565517_ark,&
      0.0008262448_ark,0.0_ark,-0.0001036317_ark,0.0025090209_ark,0.0002105614_ark,&
      -0.0000204982_ark,0.0_ark,-0.0000414402_ark,-0.0000676532_ark,0.0000918240_ark,&
      -0.0000400594_ark,0.0_ark,-0.0006468368_ark,0.0001397619_ark,0.0005356017_ark,&
      0.0001585601_ark,-0.0001899817_ark,-0.0015914516_ark,-0.0002822918_ark,&
      0.0_ark,0.0004567782_ark,0.0_ark,0.0001064879_ark,0.0_ark,0.0_ark,0.0_ark,-0.0001250236_ark,&
      0.0000077559_ark,0.0007064063_ark,0.0_ark,0.0_ark,0.0_ark,0.0006102666_ark,-0.0004987995_ark,&
      0.0_ark,-0.0001943536_ark,-0.0002855510_ark,0.0_ark,-0.0002176976_ark,0.0002410759_ark,&
      -0.0001644075_ark,-0.0001182853_ark,0.0_ark,0.0_ark,-0.0000272857_ark,0.0_ark/)
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

      real(ark),parameter :: SCALE=1e-6_ark
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
         v1=-3263.067522028298_ark-9.17958228916955_ark*xm**2-&
      0.5165611032611932_ark*xm**4-0.5157212949525876_ark*xm**6+&
      18.31161578215203_ark*xp2+&
      30.14193751791963_ark*xm**2*xp2-13.62543868575853_ark*xm**4*xp2-&
      43.37232019119388_ark*xp2**2-50.50353364079294_ark*xm**2*xp2**2+&
      70.36441193443143_ark*xm**4*xp2**2+27.43935454999898_ark*xp2**3+&
      123.751990625258_ark*xm**2*xp2**3-76.80240321256033_ark*xm**4*xp2**3-&
      9.50017804016001_ark*xp2**4-363.4487347625543_ark*xm**2*xp2**4+&
      113.1940248029751_ark*xp2**5+376.6560011408163_ark*xm**2*xp2**5-&
      164.6523756673548_ark*xp2**6+9.16256842998227_ark*zm&
       -1.22230095639504_ark*xm**2*zm+1.33032571356463_ark*xp2*&
      zm-0.94822119654751_ark*xm**2*xp2*zm+0.7645470802285307_ark*xm**4*xp2*&
      zm-11.77270680473595_ark*xp2**2*zm-0.4065994514809928_ark*xm**4*xp2**2*&
      zm-2.113651214829342_ark*xp2**3*zm-3.653921741665064_ark*xm**2*xp2**4*&
      zm+26.53983199106825_ark*xp2**6*zm+3.099164302936567_ark*zm**2
      v2=-0.4668245990549825_ark*xm**2*zm**2+0.05845413180128318_ark*xm**4*zm**2&
      +2.708722250876111_ark*xp2*zm**2-2.578482367020144_ark*xm**2*xp2*zm**2-&
      0.1605392233404811_ark*xm**4*xp2*zm**2-10.57780429022803_ark*xp2**2*&
      zm**2-3.496293826717189_ark*xp2**3*zm**2+23.46280699747645_ark*xp2**4*&
      zm**2+1.8547816858377_ark*zm**3-0.4003662844685243_ark*xm**2*zm**3+&
      3.040229985315839_ark*xp2*zm**3-4.955739113923876_ark*xp2**2*&
      zm**3+14.05364889791468_ark*xp2**3*zm**3-21.6926320924828_ark*xp2**4*&
      zm**3-1.321464834042384_ark*zm**4+2.298844571392118_ark*xm**2*zm**4-&
      2.633405421645483_ark*xp2*zm**4+20.97178840867901_ark*xp2**2*&
      zm**4-32.18658937476802_ark*xp2**3*zm**4-0.5992225949734171_ark*zm**5+&
      2.059827452250273_ark*xp2*zm**5+0.6453850286056735_ark*xp2**2*&
      zm**5-0.4620689505336259_ark*zm**6+0.7465042626807512_ark*xp2*&
      zm**6-0.1254018119377959_ark*zm**7+0.01947721364782498_ark*zm**8

      VR=V1+V2

!     SCALE AND SHIFT THE ZERO
      VR=VR*1.0e-6_ark
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




      subroutine poten(xp,v,r1,r2,th)

      implicit none
      integer,parameter           ::  npropin = 246

      real(ark),intent(in) ::  r1,r2,th
      real(ark),intent(out)::  v
      real(ark)            :: reoh,thetae,b1,roh,alphaoh,alpha
      real(ark)            :: phh2,ut,x0,ut2,x02,xs1,xs2,xst,&
      rr1,rr2,xep1,xep2,xep3
      real(ark)            :: rhh,vhh,vpb1,vpb2,v0,vp1,vp2,vp3,&
      vp,vps1,vps2,y1,y2,y12,y22,voh1,voh2


      real(ark) ::  xp(1:npropin)
      integer(ik) :: i 
      !
           reoh=0.9586d0
           thetae=104.48d0
           b1=2.15d0
           alphaoh=0.0d0
           phh2=6.20164303995d0

      thetae=thetae*.314159265358979312d01*.00555555555555555555d0

      xs1=(r1+r2)*0.5d0-reoh
      xs2=(r1-r2)*0.5d0
      xst=cos(th)-cos(thetae)

      rr1=r1-reoh
      rr2=r2-reoh

         alpha=2.2668d0
         xep1=(exp(-2.d0*alpha*rr1)-2.d0*exp(-alpha*rr1))* &
     (0.4389830771267D+05                                    &
     -0.8606164763375D+03*rr1                                &           
     +0.1746549380317D+05*rr1**2                             &
     -0.3326053761996D+05*rr1**3                             &
     +0.4988743744461D+04*rr1**4                             &
     +0.8751979626215D+04*rr1**5                             &
     -0.3770596827415D+04*rr1**6                             &
     +0.4270854122634D+03*rr1**7) + 0.4389830771267D+05      
                                                              
        xep2=(exp(-2.d0*alpha*rr2)-2.d0*exp(-alpha*rr2))*  &
     (0.4389830771267D+05                                    &
     -0.8606164763375D+03*rr2                                &
     +0.1746549380317D+05*rr2**2                             &
     -0.3326053761996D+05*rr2**3                             &
     +0.4988743744461D+04*rr2**4                             &
     +0.8751979626215D+04*rr2**5                             &
     -0.3770596827415D+04*rr2**6                             &
     +0.4270854122634D+03*rr2**7) + 0.4389830771267D+05      
                                                              
        xep3=exp(-b1*((r1-reoh)**2+(r2-reoh)**2))            
        rhh=sqrt(r1**2+r2**2-2.d0*r1*r2*cos(th))            
        vhh=0.820894739261131734D+06*exp(-phh2*rhh)         
                                                              
     v0 = xp(  1)  *xs1**0 *xs2**0 *xst**0                    
     vp1=+xp(  2)  *xs1**0 *xs2**0 *xst**1                   &
     +xp(  3)  *xs1**1 *xs2**0 *xst**0                       &
     +xp(  4)  *xs1**0 *xs2**0 *xst**2                       &
     +xp(  5)  *xs1**0 *xs2**2 *xst**0                       &
     +xp(  6)  *xs1**1 *xs2**0 *xst**1                       &
     +xp(  7)  *xs1**2 *xs2**0 *xst**0                       &
     +xp(  8)  *xs1**0 *xs2**0 *xst**3                       &
     +xp(  9)  *xs1**0 *xs2**2 *xst**1                       &
     +xp( 10)  *xs1**1 *xs2**0 *xst**2                       &
     +xp( 11)  *xs1**1 *xs2**2 *xst**0                       &
     +xp( 12)  *xs1**2 *xs2**0 *xst**1                       &
     +xp( 13)  *xs1**3 *xs2**0 *xst**0                       &
     +xp( 14)  *xs1**0 *xs2**0 *xst**4                       &
     +xp( 15)  *xs1**0 *xs2**2 *xst**2                       &
     +xp( 16)  *xs1**0 *xs2**4 *xst**0                       &
     +xp( 17)  *xs1**1 *xs2**0 *xst**3                       &
     +xp( 18)  *xs1**1 *xs2**2 *xst**1                       &
     +xp( 19)  *xs1**2 *xs2**0 *xst**2                       &
     +xp( 20)  *xs1**2 *xs2**2 *xst**0                       &
     +xp( 21)  *xs1**3 *xs2**0 *xst**1                       &
     +xp( 22)  *xs1**4 *xs2**0 *xst**0                       &
     +xp( 23)  *xs1**0 *xs2**0 *xst**5                       &
     +xp( 24)  *xs1**0 *xs2**2 *xst**3                       &
     +xp( 25)  *xs1**0 *xs2**4 *xst**1                       &
     +xp( 26)  *xs1**1 *xs2**0 *xst**4                       &
     +xp( 27)  *xs1**1 *xs2**2 *xst**2                       &
     +xp( 28)  *xs1**1 *xs2**4 *xst**0                       &
     +xp( 29)  *xs1**2 *xs2**0 *xst**3                       &
     +xp( 30)  *xs1**2 *xs2**2 *xst**1                       &
     +xp( 31)  *xs1**3 *xs2**0 *xst**2                       &
     +xp( 32)  *xs1**3 *xs2**2 *xst**0                       &
     +xp( 33)  *xs1**4 *xs2**0 *xst**1                       &
     +xp( 34)  *xs1**5 *xs2**0 *xst**0                       &
     +xp( 35)  *xs1**0 *xs2**0 *xst**6                       &
     +xp( 36)  *xs1**0 *xs2**2 *xst**4                       &
     +xp( 37)  *xs1**0 *xs2**4 *xst**2                       &
     +xp( 38)  *xs1**0 *xs2**6 *xst**0                       &
     +xp( 39)  *xs1**1 *xs2**0 *xst**5                       &
     +xp( 40)  *xs1**1 *xs2**2 *xst**3                       &
     +xp( 41)  *xs1**1 *xs2**4 *xst**1                       &
     +xp( 42)  *xs1**2 *xs2**0 *xst**4                       &
     +xp( 43)  *xs1**2 *xs2**2 *xst**2                       &
     +xp( 44)  *xs1**2 *xs2**4 *xst**0                       &
     +xp( 45)  *xs1**3 *xs2**0 *xst**3                       &
     +xp( 46)  *xs1**3 *xs2**2 *xst**1                       &
     +xp( 47)  *xs1**4 *xs2**0 *xst**2                       &
     +xp( 48)  *xs1**4 *xs2**2 *xst**0                       &
     +xp( 49)  *xs1**5 *xs2**0 *xst**1                       &
     +xp( 50)  *xs1**6 *xs2**0 *xst**0                       &
     +xp( 51)  *xs1**0 *xs2**0 *xst**7                       &
     +xp( 52)  *xs1**0 *xs2**2 *xst**5                       &
     +xp( 53)  *xs1**0 *xs2**4 *xst**3                       &
     +xp( 54)  *xs1**0 *xs2**6 *xst**1                       &
     +xp( 55)  *xs1**1 *xs2**0 *xst**6                       &
     +xp( 56)  *xs1**1 *xs2**2 *xst**4                       &
     +xp( 57)  *xs1**1 *xs2**4 *xst**2                       &
     +xp( 58)  *xs1**1 *xs2**6 *xst**0                       &
     +xp( 59)  *xs1**2 *xs2**0 *xst**5                       &
     +xp( 60)  *xs1**2 *xs2**2 *xst**3                       &
     +xp( 61)  *xs1**2 *xs2**4 *xst**1                       &
     +xp( 62)  *xs1**3 *xs2**0 *xst**4                       &
     +xp( 63)  *xs1**3 *xs2**2 *xst**2                       &
     +xp( 64)  *xs1**3 *xs2**4 *xst**0                       &
     +xp( 65)  *xs1**4 *xs2**0 *xst**3                       &
     +xp( 66)  *xs1**4 *xs2**2 *xst**1                       &
     +xp( 67)  *xs1**5 *xs2**0 *xst**2                       &
     +xp( 68)  *xs1**5 *xs2**2 *xst**0                       &
     +xp( 69)  *xs1**6 *xs2**0 *xst**1                       &
     +xp( 70)  *xs1**7 *xs2**0 *xst**0                       &
     +xp( 71)  *xs1**0 *xs2**0 *xst**8                       &
     +xp( 72)  *xs1**0 *xs2**2 *xst**6                       &
     +xp( 73)  *xs1**0 *xs2**4 *xst**4                       &
     +xp( 74)  *xs1**0 *xs2**6 *xst**2                       &
     +xp( 75)  *xs1**0 *xs2**8 *xst**0                       &
     +xp( 76)  *xs1**1 *xs2**0 *xst**7                       &
     +xp( 77)  *xs1**1 *xs2**2 *xst**5                       &
     +xp( 78)  *xs1**1 *xs2**4 *xst**3                       &
     +xp( 79)  *xs1**1 *xs2**6 *xst**1                       &
     +xp( 80)  *xs1**2 *xs2**0 *xst**6                       &
     +xp( 81)  *xs1**2 *xs2**2 *xst**4                       &
     +xp( 82)  *xs1**2 *xs2**4 *xst**2                       &
     +xp( 83)  *xs1**2 *xs2**6 *xst**0                       &
     +xp( 84)  *xs1**3 *xs2**0 *xst**5                       &
     +xp( 85)  *xs1**3 *xs2**2 *xst**3                       &
     +xp( 86)  *xs1**3 *xs2**4 *xst**1                       &
     +xp( 87)  *xs1**4 *xs2**0 *xst**4                       &
     +xp( 88)  *xs1**4 *xs2**2 *xst**2                       &
     +xp( 89)  *xs1**4 *xs2**4 *xst**0                       &
     +xp( 90)  *xs1**5 *xs2**0 *xst**3                       &
     +xp( 91)  *xs1**5 *xs2**2 *xst**1                       &
     +xp( 92)  *xs1**6 *xs2**0 *xst**2                        
     vp2=+xp( 93)  *xs1**6 *xs2**2 *xst**0                   &
     +xp( 94)  *xs1**7 *xs2**0 *xst**1                       &
     +xp( 95)  *xs1**8 *xs2**0 *xst**0                       &
     +xp( 96)  *xs1**0 *xs2**0 *xst**9                       &
     +xp( 97)  *xs1**0 *xs2**2 *xst**7                       &
     +xp( 98)  *xs1**0 *xs2**4 *xst**5                       &
     +xp( 99)  *xs1**0 *xs2**6 *xst**3                       &
     +xp(100)  *xs1**0 *xs2**8 *xst**1                       &
     +xp(101)  *xs1**1 *xs2**0 *xst**8                       &
     +xp(102)  *xs1**1 *xs2**2 *xst**6                       &
     +xp(103)  *xs1**1 *xs2**4 *xst**4                       &
     +xp(104)  *xs1**1 *xs2**6 *xst**2                       &
     +xp(105)  *xs1**1 *xs2**8 *xst**0                       &
     +xp(106)  *xs1**2 *xs2**0 *xst**7                       &
     +xp(107)  *xs1**2 *xs2**2 *xst**5                       &
     +xp(108)  *xs1**2 *xs2**4 *xst**3                       &
     +xp(109)  *xs1**2 *xs2**6 *xst**1                       &
     +xp(110)  *xs1**3 *xs2**0 *xst**6                       &
     +xp(111)  *xs1**3 *xs2**2 *xst**4                       &
     +xp(112)  *xs1**3 *xs2**4 *xst**2                       &
     +xp(113)  *xs1**3 *xs2**6 *xst**0                       &
     +xp(114)  *xs1**4 *xs2**0 *xst**5                       &
     +xp(115)  *xs1**4 *xs2**2 *xst**3                       &
     +xp(116)  *xs1**4 *xs2**4 *xst**1                       &
     +xp(117)  *xs1**5 *xs2**0 *xst**4                       &
     +xp(118)  *xs1**5 *xs2**2 *xst**2                       &
     +xp(119)  *xs1**5 *xs2**4 *xst**0                       &
     +xp(120)  *xs1**6 *xs2**0 *xst**3                       &
     +xp(121)  *xs1**6 *xs2**2 *xst**1                       &
     +xp(122)  *xs1**7 *xs2**0 *xst**2                       &
     +xp(123)  *xs1**7 *xs2**2 *xst**0                       &
     +xp(124)  *xs1**8 *xs2**0 *xst**1                       &
     +xp(125)  *xs1**9 *xs2**0 *xst**0                       &
     +xp(126)  *xs1**0 *xs2**2 *xst**8                       &
     +xp(127)  *xs1**0 *xs2**4 *xst**6                       &
     +xp(128)  *xs1**0 *xs2**6 *xst**4                       &
     +xp(129)  *xs1**0 *xs2**8 *xst**2                       &
     +xp(130)  *xs1**1 *xs2**0 *xst**9                       &
     +xp(131)  *xs1**1 *xs2**2 *xst**7                       &
     +xp(132)  *xs1**1 *xs2**4 *xst**5                       &
     +xp(133)  *xs1**1 *xs2**6 *xst**3                       &
     +xp(134)  *xs1**1 *xs2**8 *xst**1                       &
     +xp(135)  *xs1**2 *xs2**0 *xst**8                       &
     +xp(136)  *xs1**2 *xs2**2 *xst**6                       &
     +xp(137)  *xs1**2 *xs2**4 *xst**4                       &
     +xp(138)  *xs1**2 *xs2**6 *xst**2                       &
     +xp(139)  *xs1**3 *xs2**0 *xst**7                       &
     +xp(140)  *xs1**3 *xs2**2 *xst**5                       &
     +xp(141)  *xs1**3 *xs2**4 *xst**3                       &
     +xp(142)  *xs1**3 *xs2**6 *xst**1                       &
     +xp(143)  *xs1**4 *xs2**0 *xst**6                       &
     +xp(144)  *xs1**4 *xs2**2 *xst**4                       &
     +xp(145)  *xs1**4 *xs2**4 *xst**2                       &
     +xp(146)  *xs1**5 *xs2**0 *xst**5                       &
     +xp(147)  *xs1**5 *xs2**2 *xst**3                       &
     +xp(148)  *xs1**5 *xs2**4 *xst**1                       &
     +xp(149)  *xs1**6 *xs2**0 *xst**4                       &
     +xp(150)  *xs1**6 *xs2**2 *xst**2                       &
     +xp(151)  *xs1**6 *xs2**4 *xst**0                       &
     +xp(152)  *xs1**7 *xs2**0 *xst**3                       &
     +xp(153)  *xs1**7 *xs2**2 *xst**1                       &
     +xp(154)  *xs1**8 *xs2**0 *xst**2                       &
     +xp(155)  *xs1**8 *xs2**2 *xst**0                       &
     +xp(156)  *xs1**9 *xs2**0 *xst**1                       &
     +xp(157)  *xs1**10*xs2**0 *xst**0                       &
     +xp(158)  *xs1**0 *xs2**2 *xst**9                       &
     +xp(159)  *xs1**0 *xs2**4 *xst**7                       &
     +xp(160)  *xs1**0 *xs2**6 *xst**5                       &
     +xp(161)  *xs1**0 *xs2**8 *xst**3                       &
     +xp(162)  *xs1**0 *xs2**10*xst**1                       &
     +xp(163)  *xs1**1 *xs2**0 *xst**10                      &
     +xp(164)  *xs1**1 *xs2**2 *xst**8                       &
     +xp(165)  *xs1**1 *xs2**4 *xst**6                       &
     +xp(166)  *xs1**1 *xs2**6 *xst**4                       &
     +xp(167)  *xs1**1 *xs2**8 *xst**2                       &
     +xp(168)  *xs1**1 *xs2**10*xst**0                       &
     +xp(169)  *xs1**2 *xs2**0 *xst**9                       &
     +xp(170)  *xs1**2 *xs2**2 *xst**7                       &
     +xp(171)  *xs1**2 *xs2**4 *xst**5                       &
     +xp(172)  *xs1**2 *xs2**6 *xst**3                       &
     +xp(173)  *xs1**2 *xs2**8 *xst**1                       &
     +xp(174)  *xs1**3 *xs2**0 *xst**8                        
     vp3=+xp(175)  *xs1**3 *xs2**2 *xst**6                   &
     +xp(176)  *xs1**3 *xs2**4 *xst**4                       &
     +xp(177)  *xs1**3 *xs2**6 *xst**2                       &
     +xp(178)  *xs1**3 *xs2**8 *xst**0                       &
     +xp(179)  *xs1**4 *xs2**0 *xst**7                       &
     +xp(180)  *xs1**4 *xs2**2 *xst**5                       &
     +xp(181)  *xs1**4 *xs2**4 *xst**3                       &
     +xp(182)  *xs1**4 *xs2**6 *xst**1                       &
     +xp(183)  *xs1**5 *xs2**0 *xst**6                       &
     +xp(184)  *xs1**5 *xs2**2 *xst**4                       &
     +xp(185)  *xs1**5 *xs2**4 *xst**2                       &
     +xp(186)  *xs1**5 *xs2**6 *xst**0                       &
     +xp(187)  *xs1**6 *xs2**0 *xst**5                       &
     +xp(188)  *xs1**6 *xs2**2 *xst**3                       &
     +xp(189)  *xs1**6 *xs2**4 *xst**1                       &
     +xp(190)  *xs1**7 *xs2**0 *xst**4                       &
     +xp(191)  *xs1**7 *xs2**2 *xst**2                       &
     +xp(192)  *xs1**8 *xs2**0 *xst**3                       &
     +xp(193)  *xs1**8 *xs2**2 *xst**1                       &
     +xp(194)  *xs1**9 *xs2**0 *xst**2                       &
     +xp(195)  *xs1**9 *xs2**2 *xst**0                       &
     +xp(196)  *xs1**10*xs2**0 *xst**1                       &
     +xp(197)  *xs1**11*xs2**0 *xst**0                       &
     +xp(198)  *xs1**0 *xs2**0 *xst**2                       &
     +xp(199)  *xs1**0 *xs2**2 *xst**0                       &
     +xp(200)  *xs1**0 *xs2**4 *xst**8                       &
     +xp(201)  *xs1**0 *xs2**6 *xst**6                       &
     +xp(202)  *xs1**0 *xs2**8 *xst**4                       &
     +xp(203)  *xs1**0 *xs2**10*xst**2                       &
     +xp(204)  *xs1**0 *xs2**12*xst**0                       &
     +xp(205)  *xs1**1 *xs2**0 *xst**11                      &
     +xp(206)  *xs1**1 *xs2**2 *xst**9                       &
     +xp(207)  *xs1**1 *xs2**4 *xst**7                       &
     +xp(208)  *xs1**1 *xs2**6 *xst**5                       &
     +xp(209)  *xs1**1 *xs2**8 *xst**3                       &
     +xp(210)  *xs1**1 *xs2**10*xst**1                       &
     +xp(211)  *xs1**2 *xs2**0 *xst**10                      &
     +xp(212)  *xs1**2 *xs2**2 *xst**8                       &
     +xp(213)  *xs1**2 *xs2**4 *xst**6                       &
     +xp(214)  *xs1**2 *xs2**6 *xst**4                       &
     +xp(215)  *xs1**2 *xs2**8 *xst**2                       &
     +xp(216)  *xs1**2 *xs2**10*xst**0                       &
     +xp(217)  *xs1**3 *xs2**0 *xst**9                       &
     +xp(218)  *xs1**3 *xs2**2 *xst**7                       &
     +xp(219)  *xs1**3 *xs2**4 *xst**5                       &
     +xp(220)  *xs1**3 *xs2**6 *xst**3                       &
     +xp(221)  *xs1**3 *xs2**8 *xst**1                       &
     +xp(222)  *xs1**4 *xs2**0 *xst**8                       &
     +xp(223)  *xs1**4 *xs2**2 *xst**6                       &
     +xp(224)  *xs1**4 *xs2**4 *xst**4                       &
     +xp(225)  *xs1**4 *xs2**6 *xst**2                       &
     +xp(226)  *xs1**4 *xs2**8 *xst**0                       &
     +xp(227)  *xs1**5 *xs2**0 *xst**7                       &
     +xp(228)  *xs1**5 *xs2**2 *xst**5                       &
     +xp(229)  *xs1**5 *xs2**4 *xst**3                       &
     +xp(230)  *xs1**5 *xs2**6 *xst**1                       &
     +xp(231)  *xs1**6 *xs2**0 *xst**6                       &
     +xp(232)  *xs1**6 *xs2**2 *xst**4                       &
     +xp(233)  *xs1**6 *xs2**4 *xst**2                       &
     +xp(234)  *xs1**6 *xs2**6 *xst**0                       &
     +xp(235)  *xs1**7 *xs2**0 *xst**5                       &
     +xp(236)  *xs1**7 *xs2**2 *xst**3                       &
     +xp(237)  *xs1**7 *xs2**4 *xst**1                       &
     +xp(238)  *xs1**8 *xs2**0 *xst**4                       &
     +xp(239)  *xs1**8 *xs2**2 *xst**2                       &
     +xp(240)  *xs1**8 *xs2**4 *xst**0                       &
     +xp(241)  *xs1**9 *xs2**0 *xst**3                       &
     +xp(242)  *xs1**9 *xs2**2 *xst**1                       &
     +xp(243)  *xs1**10*xs2**0 *xst**2                       &
     +xp(244)  *xs1**10*xs2**2 *xst**0                       &
     +xp(245)  *xs1**11*xs2**0 *xst**1                       &
     +xp(246)  *xs1**12*xs2**0 *xst**0          
                                       
       vp=vp1+vp2+vp3

        v=v0+vp*xep3+vhh+xep1+xep2

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



      SUBROUTINE POTS(V,rij1,rij2,rij3)
      implicit real(ark) (a-h,o-y),logical(z)
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
      real(rk)    :: rad
      integer(ik) :: i,idx(245,3),idxm(9,3),ifirst,j
      dimension c5z(245),cbasis(245),ccore(245),&
               crest(245),fmat(15,3),cmass(9)

      dimension c5z_(245),&
               crest_(245),cmass_(9)

!
!     expansion indicies
!
       data (idx(i,1),i=1,245)/                                         &
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,        &
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,        &
     2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,        &
     3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,        &
     3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,        &
     4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,        &
     4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,        &
     6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5,        &
     6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5,        &
     5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7,        &
     7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6,        &
     6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9,        &
     9, 9, 9, 9, 9/                                                      
     data (idx(i,2),i=1,245)/                                           &
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,        &
     1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,        &
     2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,        &
     2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3,        &
     3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,        &
     2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,        &
     3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,        &
     1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3,        &
     2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4,        &
     4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,        &
     2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,        &
     4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1,        &
     1, 1, 1, 1, 1/                                                      
     data (idx(i,3),i=1,245)/                                           &
     1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 1, 2, 3, 4, 5,        &
     6, 7, 8, 9,10,11,12,13,14, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,        &
    12,13, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1, 2, 3, 4, 5,        &
     6, 7, 8, 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 1,        &
     2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,        &
    11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8,        &
     9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 7, 8,        &
     9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9,        &
     1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2,        &
     3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6,        &
     7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3,        &
     4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2,        &
     3, 4, 5, 6, 7/                                                    
!                                                                        
!   expansion coefficients for 5z ab initio data                       
!                                                                        
     data (c5z_(i),i=1,45)/                                              &
     4.2278462684916D+04, 4.5859382909906D-02, 9.4804986183058D+03,     &
     7.5485566680955D+02, 1.9865052511496D+03, 4.3768071560862D+02,     &
     1.4466054104131D+03, 1.3591924557890D+02,-1.4299027252645D+03,     &
     6.6966329416373D+02, 3.8065088734195D+03,-5.0582552618154D+02,     &
    -3.2067534385604D+03, 6.9673382568135D+02, 1.6789085874578D+03,     &
    -3.5387509130093D+03,-1.2902326455736D+04,-6.4271125232353D+03,     &
    -6.9346876863641D+03,-4.9765266152649D+02,-3.4380943579627D+03,     &
     3.9925274973255D+03,-1.2703668547457D+04,-1.5831591056092D+04,     &
     2.9431777405339D+04, 2.5071411925779D+04,-4.8518811956397D+04,     &
    -1.4430705306580D+04, 2.5844109323395D+04,-2.3371683301770D+03,     &
     1.2333872678202D+04, 6.6525207018832D+03,-2.0884209672231D+03,     &
    -6.3008463062877D+03, 4.2548148298119D+04, 2.1561445953347D+04,     &
    -1.5517277060400D+05, 2.9277086555691D+04, 2.6154026873478D+05,     &
    -1.3093666159230D+05,-1.6260425387088D+05, 1.2311652217133D+05,     &
    -5.1764697159603D+04, 2.5287599662992D+03, 3.0114701659513D+04/    
                                                                         
     data (c5z_(i),i=46,90)/                                             &
    -2.0580084492150D+03, 3.3617940269402D+04, 1.3503379582016D+04,     &
    -1.0401149481887D+05,-6.3248258344140D+04, 2.4576697811922D+05,     &
     8.9685253338525D+04,-2.3910076031416D+05,-6.5265145723160D+04,     &
     8.9184290973880D+04,-8.0850272976101D+03,-3.1054961140464D+04,     &
    -1.3684354599285D+04, 9.3754012976495D+03,-7.4676475789329D+04,     &
    -1.8122270942076D+05, 2.6987309391410D+05, 4.0582251904706D+05,     &
    -4.7103517814752D+05,-3.6115503974010D+05, 3.2284775325099D+05,     &
     1.3264691929787D+04, 1.8025253924335D+05,-1.2235925565102D+04,     &
    -9.1363898120735D+03,-4.1294242946858D+04,-3.4995730900098D+04,     &
     3.1769893347165D+05, 2.8395605362570D+05,-1.0784536354219D+06,     &
    -5.9451106980882D+05, 1.5215430060937D+06, 4.5943167339298D+05,     &
    -7.9957883936866D+05,-9.2432840622294D+04, 5.5825423140341D+03,     &
     3.0673594098716D+03, 8.7439532014842D+04, 1.9113438435651D+05,     &
    -3.4306742659939D+05,-3.0711488132651D+05, 6.2118702580693D+05,     &
    -1.5805976377422D+04,-4.2038045404190D+05, 3.4847108834282D+05/    
                                                                         
     data (c5z_(i),i=91,135)/                                            &
    -1.3486811106770D+04, 3.1256632170871D+04, 5.3344700235019D+03,     &
     2.6384242145376D+04, 1.2917121516510D+05,-1.3160848301195D+05,     &
    -4.5853998051192D+05, 3.5760105069089D+05, 6.4570143281747D+05,     &
    -3.6980075904167D+05,-3.2941029518332D+05,-3.5042507366553D+05,     &
     2.1513919629391D+03, 6.3403845616538D+04, 6.2152822008047D+04,     &
    -4.8805335375295D+05,-6.3261951398766D+05, 1.8433340786742D+06,     &
     1.4650263449690D+06,-2.9204939728308D+06,-1.1011338105757D+06,     &
     1.7270664922758D+06, 3.4925947462024D+05,-1.9526251371308D+04,     &
    -3.2271030511683D+04,-3.7601575719875D+05, 1.8295007005531D+05,     &
     1.5005699079799D+06,-1.2350076538617D+06,-1.8221938812193D+06,     &
     1.5438780841786D+06,-3.2729150692367D+03, 1.0546285883943D+04,     &
    -4.7118461673723D+04,-1.1458551385925D+05, 2.7704588008958D+05,     &
     7.4145816862032D+05,-6.6864945408289D+05,-1.6992324545166D+06,     &
     6.7487333473248D+05, 1.4361670430046D+06,-2.0837555267331D+05,     &
     4.7678355561019D+05,-1.5194821786066D+04,-1.1987249931134D+05/    
                                                                       
                                                                         
     data (c5z_(i),i=136,180)/                                           &
     1.3007675671713D+05, 9.6641544907323D+05,-5.3379849922258D+05,     &
    -2.4303858824867D+06, 1.5261649025605D+06, 2.0186755858342D+06,     &
    -1.6429544469130D+06,-1.7921520714752D+04, 1.4125624734639D+04,     &
    -2.5345006031695D+04, 1.7853375909076D+05,-5.4318156343922D+04,     &
    -3.6889685715963D+05, 4.2449670705837D+05, 3.5020329799394D+05,     &
     9.3825886484788D+03,-8.0012127425648D+05, 9.8554789856472D+04,     &
     4.9210554266522D+05,-6.4038493953446D+05,-2.8398085766046D+06,     &
     2.1390360019254D+06, 6.3452935017176D+06,-2.3677386290925D+06,     &
    -3.9697874352050D+06,-1.9490691547041D+04, 4.4213579019433D+04,     &
     1.6113884156437D+05,-7.1247665213713D+05,-1.1808376404616D+06,     &
     3.0815171952564D+06, 1.3519809705593D+06,-3.4457898745450D+06,     &
     2.0705775494050D+05,-4.3778169926622D+05, 8.7041260169714D+03,     &
     1.8982512628535D+05,-2.9708215504578D+05,-8.8213012222074D+05,     &
     8.6031109049755D+05, 1.0968800857081D+06,-1.0114716732602D+06,     &
     1.9367263614108D+05, 2.8678295007137D+05,-9.4347729862989D+04/    
                                                                         
     data (c5z_(i),i=181,225)/                                           &
     4.4154039394108D+04, 5.3686756196439D+05, 1.7254041770855D+05,     &
    -2.5310674462399D+06,-2.0381171865455D+06, 3.3780796258176D+06,     &
     7.8836220768478D+05,-1.5307728782887D+05,-3.7573362053757D+05,     &
     1.0124501604626D+06, 2.0929686545723D+06,-5.7305706586465D+06,     &
    -2.6200352535413D+06, 7.1543745536691D+06,-1.9733601879064D+04,     &
     8.5273008477607D+04, 6.1062454495045D+04,-2.2642508675984D+05,     &
     2.4581653864150D+05,-9.0376851105383D+05,-4.4367930945690D+05,     &
     1.5740351463593D+06, 2.4563041445249D+05,-3.4697646046367D+03,     &
    -2.1391370322552D+05, 4.2358948404842D+05, 5.6270081955003D+05,     &
    -8.5007851251980D+05,-6.1182429537130D+05, 5.6690751824341D+05,     &
    -3.5617502919487D+05,-8.1875263381402D+02,-2.4506258140060D+05,     &
     2.5830513731509D+05, 6.0646114465433D+05,-6.9676584616955D+05,     &
     5.1937406389690D+05, 1.7261913546007D+05,-1.7405787307472D+04,     &
    -3.8301842660567D+05, 5.4227693205154D+05, 2.5442083515211D+06,     &
    -1.1837755702370D+06,-1.9381959088092D+06,-4.0642141553575D+05/    
                                                                       
                                                                         
     data (c5z_(i),i=226,245)/                                           &
     1.1840693827934D+04,-1.5334500255967D+05, 4.9098619510989D+05,     &
     6.1688992640977D+05, 2.2351144690009D+05,-1.8550462739570D+06,     &
     9.6815110649918D+03,-8.1526584681055D+04,-8.0810433155289D+04,     &
     3.4520506615177D+05, 2.5509863381419D+05,-1.3331224992157D+05,     &
    -4.3119301071653D+05,-5.9818343115856D+04, 1.7863692414573D+03,     &
     8.9440694919836D+04,-2.5558967650731D+05,-2.2130423988459D+04,     &
     4.4973674518316D+05,-2.2094939343618D+05/                         
!                                                                        
!   expansion coefficients for basis correction                        
!                                                                        
     data (cbasis(i),i=1,45)/                                           &
     6.9770019624764D-04,-2.4209870001642D+01, 1.8113927151562D+01,     &
     3.5107416275981D+01,-5.4600021126735D+00,-4.8731149608386D+01,     &
     3.6007189184766D+01, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
    -7.7178474355102D+01,-3.8460795013977D+01,-4.6622480912340D+01,     &
     5.5684951167513D+01, 1.2274939911242D+02,-1.4325154752086D+02,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00,-6.0800589055949D+00,     &
     8.6171499453475D+01,-8.4066835441327D+01,-5.8228085624620D+01,     &
     2.0237393793875D+02, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     3.3525582670313D+02, 7.0056962392208D+01,-4.5312502936708D+01/    
                                                                       
                                                                         
     data (cbasis(i),i=46,90)/                                          &
    -3.0441141194247D+02, 2.8111438108965D+02, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-1.2983583774779D+02, 3.9781671212935D+01,     &
    -6.6793945229609D+01,-1.9259805675433D+02, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-8.2855757669957D+02,-5.7003072730941D+01,     &
    -3.5604806670066D+01, 9.6277766002709D+01, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 8.8645622149112D+02,-7.6908409772041D+01,     &
     6.8111763314154D+01, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                       
                                                                       
                                                                         
     data (cbasis(i),i=91,135)/                                         &
     2.5090493428062D+02,-2.3622141780572D+02, 5.8155647658455D+02,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 2.8919570295095D+03,     &
    -1.7871014635921D+02,-1.3515667622500D+02, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-3.6965613754734D+03, 2.1148158286617D+02,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00,-1.4795670139431D+03,     &
     3.6210798138768D+02, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
    -5.3552886800881D+03, 3.1006384016202D+02, 0.0000000000000D+00/    
                                                                       
                                                                       
                                                                       
                                                                         
     data (cbasis(i),i=136,180)/                                        &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 1.6241824368764D+03, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 4.3764909606382D+03, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 1.0940849243716D+03, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 3.0743267832931D+03, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                       
                                                                         
     data (cbasis(i),i=181,225)/                                        &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                       
                                                                         
     data (cbasis(i),i=226,245)/                                        &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00/                         
!                                                                        
!   expansion coefficients for core correction                         
!                                                                        
     data (ccore(i),i=1,45)/                                            &
     2.4332191647159D-02,-2.9749090113656D+01, 1.8638980892831D+01,     &
    -6.1272361746520D+00, 2.1567487597605D+00,-1.5552044084945D+01,     &
     8.9752150543954D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
    -3.5693557878741D+02,-3.0398393196894D+00,-6.5936553294576D+00,     &
     1.6056619388911D+01, 7.8061422868204D+01,-8.6270891686359D+01,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00,-3.1688002530217D+01,     &
     3.7586725583944D+01,-3.2725765966657D+01,-5.6458213299259D+00,     &
     2.1502613314595D+01, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     5.2789943583277D+02,-4.2461079404962D+00,-2.4937638543122D+01/    
                                                                       
                                                                         
     data (ccore(i),i=46,90)/                                           &
    -1.1963809321312D+02, 2.0240663228078D+02, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-6.2574211352272D+02,-6.9617539465382D+00,     &
    -5.9440243471241D+01, 1.4944220180218D+01, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-1.2851139918332D+03,-6.5043516710835D+00,     &
     4.0410829440249D+01,-6.7162452402027D+01, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 1.0031942127832D+03, 7.6137226541944D+01,     &
    -2.7279242226902D+01, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                       
                                                                         
     data (ccore(i),i=91,135)/                                          &
    -3.3059000871075D+01, 2.4384498749480D+01,-1.4597931874215D+02,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 1.6559579606045D+03,     &
     1.5038996611400D+02,-7.3865347730818D+01, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-1.9738401290808D+03,-1.4149993809415D+02,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00,-1.2756627454888D+02,     &
     4.1487702227579D+01, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
    -1.7406770966429D+03,-9.3812204399266D+01, 0.0000000000000D+00/    
                                                                       
                                                                         
     data (ccore(i),i=136,180)/                                         &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-1.1890301282216D+03, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 2.3723447727360D+03, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-1.0279968223292D+03, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 5.7153838472603D+02, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                         
     data (ccore(i),i=181,225)/                                         &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                         
     data (ccore(i),i=226,245)/                                         &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00/                         
!                                                                        
!   expansion coefficients for v rest                                  
!                                                                        
     data (crest_(i),i=1,45)/                                            &
     0.0000000000000D+00,-4.7430930170000D+00,-1.4422132560000D+01,     &
    -1.8061146510000D+01, 7.5186735000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
    -2.7962099800000D+02, 1.7616414260000D+01,-9.9741392630000D+01,     &
     7.1402447000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00,-7.8571336480000D+01,     &
     5.2434353250000D+01, 7.7696745000000D+01, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     1.7799123760000D+02, 1.4564532380000D+02, 2.2347226000000D+02/    
                                                                         
     data (crest_(i),i=46,90)/                                           &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-4.3823284100000D+02,-7.2846553000000D+02,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00,-2.6752313750000D+02, 3.6170310000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                         
     data (crest_(i),i=91,135)/                                          &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                       
                                                                         
     data (crest_(i),i=136,180)/                                         &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                       
                                                                         
     data (crest_(i),i=181,225)/                                         &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00/    
                                                                       
                                                                         
     data (crest_(i),i=226,245)/                                         &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00, 0.0000000000000D+00,     &
     0.0000000000000D+00, 0.0000000000000D+00/                         
!                                                                        
!   expansion indicies for mass correction                             
!                                                                        
     data idxm/1,2,1,1,3,2,1,2,1,                                       &
               2,1,1,3,1,2,2,1,1,                                       &
               1,1,2,1,1,1,2,2,3/                                      
!                                                                        
!   expansion coefficients for mass correction                         
!                                                                        
     data cmass_/ -8.3554183D+00,3.7036552D+01,-5.2722136D+00,           &
          1.6843857D+01,-7.0929741D+01,5.5380337D+00,-2.9962997D+01,    &
          1.3637682D+02,-3.0530195d+00/                                
!                                                                        
!   two body parameters                                                
!                                                                        
     data reoh_,thetae_,b1_,roh_,alphaoh_,deoh_,phh1_,phh2_/0.958649d0,         &
          104.3475d0,2.0d0,0.9519607159623009d0,2.587949757553683d0,    &
          42290.92019288289d0,16.94879431193463d0,12.66426998162947d0/
!
!     scaling factors for contributions to emperical potential
!
        data f5z,fbasis,fcore,frest/0.0d0,0.0d0,-1.0d0,0.0d0/
      !save
      data ifirst/0/
      !
      !if(ifirst.eq.0)then
       !ifirst=1
    !1  format(/1x,'pes for h2o',&
    !         /1x,'by Harry Partridge and David W. Schwenke',&
    !         /1x,'submitted to J. Chem. Phys. Nov. 8, 1996')
    !56  format(/1x,'parameters before adjustment')
    !55  format(/1x,'two body potential parameters:',    &
    !         /1x,'hh: phh1 = ',f10.1,' phh2 = ',f5.2,  &
    !         /1x,'oh: deoh = ',f10.1,' alpha = ',f7.4, &
    !         ' re = ',f7.4)
    !4  format(/1x,'three body parameters:',                                &
    !        /1x,'reoh = ',f10.4,' thetae = ',f10.4,                        &
    !        /1x,'betaoh = ',f10.4,                                         &
    !        /1x,'    i    j    k',7x,'c5z',9x,'cbasis',10x,'ccore',        &
    !        10x,'crest')
    !   do 2 i=1,245
    !2  continue
!
!     remove mass correction from vrest
!
       !
       c5z = c5z_
       crest = crest_
       cmass = cmass_
       reoh = reoh_
       thetae = thetae_
       b1 = b1_
       roh = roh_
       alphaoh = alphaoh_
       deoh = deoh_
       phh1 = phh1_
       phh2 = phh2_
       !
       xmh=1836.152697d0
       xmhi=1d0/xmh
       xmd=3670.483031d0
       fact=1.d0/((1.d0/xmd)-(1.d0/xmh))
   !65  format(/1x,'parameters for delta v hdo ',&
   !         /1x,'    i    j    k')
       do 60 i=1,9
        cmass(i)=cmass(i)*fact
        corr=cmass(i)*xmhi
        if(idxm(i,1).eq.idxm(i,2))corr=corr*0.5d0
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

       xm1=1.d0/1836.152697d0
       xm2=1.d0/1836.152697d0
!
!     adjust parameters using scale factors
!
   !57  format(/1x,'adjusting parameters using scale factors ',&
   !          /1x,'f5z =    ',f11.8,                           &
   !          /1x,'fbasis = ',f11.8,                           &
   !          /1x,'fcore =  ',f11.8,                           &
   !          /1x,'frest =  ',f11.8)
       phh1=phh1*f5z
       deoh=deoh*f5z
       do 59 i=1,245
        c5z(i)=f5z*c5z(i)+fbasis*cbasis(i)+fcore*ccore(i)&
            +frest*crest(i)
   59  continue

   !58  format(/1x,'three body parameters:',          &
   !          /1x,'reoh = ',f10.4,' thetae = ',f10.4, &     
   !          /1x,'betaoh = ',f10.4,                  &
   !          /1x,'    i    j    k   cijk',           &
   !          /(1x,3i5,1pe15.7))
       do 66 i=1,9
        cmass(i)=cmass(i)*frest
   66  continue

   !76  format(/1x,'mass correction factors ', &
   !          /1x,'    i    j    k   cijk', &
   !          /(1x,3i5,1pe15.7))
!
!     convert parameters from 1/cm, angstrom to a.u.
!
       reoh=reoh/0.529177249d0
       b1=b1*0.529177249d0*0.529177249d0
       do i=1,245
        c5z(i)=c5z(i)*4.556335d-6 
       enddo
       do i=1,9
        cmass(i)=cmass(i)*4.556335d-6
       enddo
       rad=acos(-1.0_ark)/180.0_ark
       ce=cos(thetae*rad)
       phh1=phh1*exp(phh2)
       phh1=phh1*4.556335d-6
       phh2=phh2*0.529177249d0
       deoh=deoh*4.556335d-6
       roh=roh/0.529177249d0
       alphaoh=alphaoh*0.529177249d0
       c5z(1)=c5z(1)*2d0
       !
      !end if
       !
       x1=(rij1-reoh)/reoh
       x2=(rij2-reoh)/reoh
       x3=cos(rij3)-ce
       rhh=sqrt(rij1**2+rij2**2 &
           -2d0*rij1*rij2*cos(rij3))
       vhh=phh1*exp(-phh2*rhh)
       ex=exp(-alphaoh*(rij1-roh))
       voh1=deoh*ex*(ex-2d0)
       ex=exp(-alphaoh*(rij2-roh))
       voh2=deoh*ex*(ex-2d0)
       fmat(1,1)=1d0
       fmat(1,2)=1d0
       fmat(1,3)=1d0
       !
       do j=2,15
        fmat(j,1)=fmat(j-1,1)*x1
        fmat(j,2)=fmat(j-1,2)*x2
        fmat(j,3)=fmat(j-1,3)*x3
       enddo
       !
       v=0d0
       do j=2,245
        term=c5z(j)*(fmat(idx(j,1),1)*fmat(idx(j,2),2)          &
                         +fmat(idx(j,2),1)*fmat(idx(j,1),2))    &
                         *fmat(idx(j,3),3)
        v=v+term
       enddo
       !
       v1=0d0
       v2=0d0
       do j=1,9
        v1=v1+cmass(j)*fmat(idxm(j,1),1)*fmat(idxm(j,2),2) &
            *fmat(idxm(j,3),3)                              
       v2=v2+cmass(j)*fmat(idxm(j,2),1)*fmat(idxm(j,1),2)  &
            *fmat(idxm(j,3),3)
       enddo
       v=v+xm1*v1+xm2*v2

       v=v*exp(-b1*((rij1-reoh)**2+(rij2-reoh)**2)) &
          +c5z(1)                                   &
          +voh1+voh2+vhh

      return
      end SUBROUTINE POTS


end module pot_user

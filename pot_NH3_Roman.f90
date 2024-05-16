!
!  This unit is for a user defined potential 
!
module pot_user
  use accuracy
  use moltype

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

 ! Check the potential name 
 subroutine MLpoten_name(name)
   !
   character(len=cl),intent(in) ::  name
   character(len=cl),parameter ::  poten_name = 'GENERAL'
   ! 
   if (poten_name/=trim(name)) then
     write(out,"('a,a,a,a')") 'Wrong Potential ',trim(name),'; should be ',trim(poten_name)
   endif
   !
   write(out,"('a')") '  Using USER-tpye PES ',trim(poten_name)
   !
 end subroutine MLpoten_name
 !

 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
       !
       f = 0
       !
  end subroutine MLdipole

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
   f = MLpoten_xy3_morbid_Roman_10(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !
 !!!!start change
  function ML_MEP_OH3P(x)  result(dst)

   real(ark),intent(in)   ::  x
   real(ark)              ::  dst
   real(ark)              ::  rf(0:6),y
   integer(ik)          ::  k(0:6) = (/0,1,2,3,4,5,6/)
     !
     rf(0:6)  = molec%mep_params(1:7) 
     !
     y=sin(x)-1.0_ark
     !
     dst =  sum(rf(:)*y**k(:))
     ! 
  end function ML_MEP_OH3P
 !!!!end change
 !
  ! Defining potential energy function 
  !
  ! This type is for XY3-molecules, MORBID type of expansion #10
  ! see, e.g., JCP 117, 11265 (2002) for the defienition 
  !
  function MLpoten_xy3_morbid_Roman_10(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  force2(304)
   real(ark)              ::  f,f1,f2
   !
   integer(ik)          :: N

   real(ark)            ::  r14,r24,r34,alpha1,alpha2,alpha3,alpha0
   real(ark)            ::  y1,y2,y3,y4,y5
   real(ark)            ::  v0,rhoe,v1,v2,v3,v4,v5,v6,rho,re14,aa1,tau
   real(ark)            ::  sinrho,alpha,coro,tau_2,cosalpha,sinphi,phi1,phi2,phi3,delta
      !
      if (verbose>=6) write(out,"('MLpoten_xy3_morbid_10/start')") 
      !
      re14    = molec%req(1)
      !
      alpha0  = molec%alphaeq(1)
      rhoe    = pi-asin(2.0_ark*sin(alpha0*0.5_ark)/sqrt(3.0_ark))
      !
      aa1     = molec%specparam(1)
      !
      r14    = local(1) ;  r24    = local(2) ;  r34    = local(3)
      !
      if (size(local)==6.and.molec%Ndihedrals==1) then 
         !
         select case(trim(molec%coords_transform))
         case default
            !
            write (out,"('MLpoten_xy3_morbid_10: coord. type ',a,' unknown')") trim(molec%coords_transform)
            stop 'MLpoten_xy3_morbid_10 - bad coord. type'
            !

         case('R-SYMPHI-TAU')
            !
            !
            alpha3 = local(4)
            alpha2 = local(5)
            !
            tau    = local(6)
            !
            cosalpha = cos(alpha2)*cos(alpha3) + sin(alpha2)*sin(alpha3)*cos(tau)
            !
            if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
               !
               write (out,"('MLpoten_xy3_morbid_45760: cosalpha>1: ',f18.8)") cosalpha
               stop 'MLpoten_xy3_morbid_45760 - bad cosalpha'
               !
            endif 
            !
            alpha1 = acos(cosalpha)
            !
         case('R-S-DELTA','R-2D-DELTA','R-S-RHO','SYM-DELTA','R-SYMPHI-DELTA','R-PHI-DELTA','R-S-DELTA-MEP','R-PHI-DELTA-MEP','R-SYMPHI-DELTA-MEP','R-A2A3-DELTA')
            !
            alpha3 = local(4)
            alpha2 = local(5)
            !
            delta = local(6)
            !
            sinphi = sin(alpha2*0.5_ark)/cos(delta)
            phi2 = asin(sinphi)*2.0_ark
            phi2 = mod(phi2+2.0_ark*pi,2.0_ark*pi)
            !
            sinphi = sin(alpha3*0.5_ark)/cos(delta)
            phi3 = asin(sinphi)*2.0_ark
            phi3 = mod(phi3+2.0_ark*pi,2.0_ark*pi)
            phi1 = 2.0_ark*pi-phi2-phi3
            !     
            cosalpha = cos( delta )**2+sin(delta)**2*cos(phi1)
            !
            if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
               !
               write (out,"('MLpoten_xy3_morbid_45760: cosalpha>1: ',f18.8)") cosalpha
               stop 'MLpoten_xy3_morbid_45760 - bad cosalpha'
               !
            elseif ( cosalpha>=1.0_ark) then 
               alpha1 = 0.0_ark
            else 
               alpha1 = acos(cosalpha)
            endif
            !
         case('R-A2-A3-TAU','R-THETA-TAU')
            !
            !
            alpha2 = local(4)
            alpha3 = local(5)
            !
            tau    = local(6)
            !
            cosalpha = cos(alpha2)*cos(alpha3) + sin(alpha2)*sin(alpha3)*cos(tau)
            !
            if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
               !
               write (out,"('MLpoten_xy3_morbid_45760: cosalpha>1: ',f18.8)") cosalpha
               stop 'MLpoten_xy3_morbid_45760 - bad cosalpha'
               !
            endif 
            !
            alpha1 = acos(cosalpha)

            !
         end select


      elseif (size(local)==7.and.molec%Ndihedrals==1) then
         !
         alpha3 = local(4)
         alpha2 = local(5)
         alpha1 = local(6)
         tau    = local(7)
         ! 
      elseif (size(local)==6.and.molec%Ndihedrals==0) then 
         !
         alpha3 = local(4)
         alpha2 = local(5)
         alpha1 = local(6)
         !
         tau = sqrt(1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
                   +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3) )
         !
      else
        write (out,"('MLpoten_xy3_morbid_delta: bad locals; local = ',30f18.8)") local(:)
        stop 'MLpoten_xy3_morbid_delta - bad local coordinates'
      endif
      !
      alpha=(alpha1+alpha2+alpha3)/3.0_ark
      !
      if ( 2.0_ark*sin(alpha*0.5_ark)/sqrt(3.0_ark).ge.1.0_ark ) then 
        sinrho=1.0_ark ; rho = 0.5_ark*pi
      else 
        sinrho = 2.0_ark*sin(alpha*0.5_ark)/sqrt(3.0_ark)
        !
        rho=pi-asin(sinrho)
        !
      endif
      !
      coro=(sin(rhoe)-sinrho)
      !
      if (trim(molec%potentype)=='POTEN_XY3_MORBID_10_MEP') re14 = ML_mep_oh3p(rho)
      !
      y1=1.0_ark-exp(-aa1*(r14-re14))
      y2=1.0_ark-exp(-aa1*(r24-re14))
      y3=1.0_ark-exp(-aa1*(r34-re14))
      !
      y4=(2.0_ark*alpha1-alpha2-alpha3)/sqrt(6.0_ark)
      y5=(alpha2-alpha3)/sqrt(2.0_ark)
      !
      f1 = poten_xy3_morbid_10(y1,y2,y3,y4,y5,coro,force)

      force2(:)=0.0_ark

force2(1)=-4501.81203917_ark
force2(2)=-19881.26329389_ark
force2(3)=-23974.57727407_ark
force2(4)=-9334.08769568_ark
force2(6)=-322.99732624_ark
force2(12)=5732.64717612_ark
force2(13)=13470.10273810_ark
force2(14)=10831.73624903_ark
force2(15)=2696.17194994_ark
force2(19)=-390.35384268_ark
force2(26)=-199529.16354149_ark
force2(27)=-744261.41099013_ark
force2(28)=-1006934.10723694_ark
force2(29)=-576639.05043630_ark
force2(30)=-114345.27125559_ark
force2(31)=-28897.81099237_ark
force2(32)=-123077.58667199_ark
force2(33)=-196547.81666239_ark
force2(34)=-139124.46991904_ark
force2(35)=-36788.63210179_ark
force2(36)=1959.66818301_ark
force2(37)=1306.96225957_ark
force2(38)=-3189.04443438_ark
force2(39)=-2357.13560627_ark
force2(40)=-1180.58134135_ark
force2(41)=-1405.50631046_ark
force2(42)=624.03401645_ark
force2(43)=853.31266981_ark
force2(44)=42446.96444091_ark
force2(45)=139604.39540308_ark
force2(46)=152760.79415013_ark
force2(47)=55622.89083767_ark
force2(48)=1928.88517275_ark
force2(49)=3306.00669222_ark
force2(50)=1337.21817462_ark
force2(52)=14899.53219932_ark
force2(53)=49070.51150982_ark
force2(54)=53871.52094367_ark
force2(55)=19720.93341366_ark
force2(56)=8002.31666938_ark
force2(57)=26001.07795769_ark
force2(58)=28166.98161352_ark
force2(59)=10165.90720340_ark
force2(60)=-7373.48068441_ark
force2(61)=-23883.57111100_ark
force2(62)=-25628.40548148_ark
force2(63)=-9138.11334295_ark
force2(64)=-7446.34665931_ark
force2(65)=-24161.39395965_ark
force2(66)=-26037.86030725_ark
force2(67)=-9324.29441909_ark
force2(68)=-2170.24940951_ark
force2(69)=-4331.72359507_ark
force2(70)=-2131.17374574_ark
force2(71)=939.09024088_ark
force2(72)=1915.74690059_ark
force2(73)=985.54062584_ark
force2(74)=-4.02427755_ark

     coro=-sinrho
      y1=r14!-re14
      y2=r24!-re14
      y3=r34!-re14
!
      f2 = poten_xy3_morbid_10(y1,y2,y3,y4,y5,coro,force2)
!      !
      f=f1+f2*force(305)
!      !
      if (verbose>=6) write(out,"('MLpoten_xy3_morbid_10/end')") 
 
 end function MLpoten_xy3_morbid_Roman_10



 function poten_xy3_morbid_10(y1,y2,y3,y4,y5,coro,force) result (f)
   !
   real(ark),intent(in) ::  y1,y2,y3,y4,y5,coro
   real(ark),intent(in) ::  force(:)
   real(ark)            ::  f
   !
   integer(ik)          :: N

   real(ark)            ::  v0,v1,v2,v3,v4,v5,v6


   real(ark)            ::  &
      fea1  ,                                      &
      fea11  ,fea12  ,fea14  ,fea44  ,                     &
      fea111 ,fea112 ,fea114 ,fea123 ,                     &
      fea124 ,fea144 ,fea155 ,fea455 ,                     &
      fea1111,fea1112,fea1114,fea1122,                     &
      fea1123,fea1124,fea1125,fea1144,                     &
      fea1155,fea1244,fea1255,fea1444,                     &
      fea1455,fea4444
   real(ark)            ::  &
      fea44444 ,fea33455 ,fea33445 ,fea33345 ,fea33344 ,&
      fea33334 ,fea33333 ,fea25555 ,fea24455 ,fea24445 ,fea23333 ,&
      fea13455 ,fea13445 ,fea13345 ,fea12355 ,fea11334 ,fea11333 ,&
      fea11255 ,fea11245 ,fea11234 ,fea11233 ,fea11135 ,fea11134 ,&
      fea11123 ,fea555555,fea444444,fea335555,fea334455,fea334445,&
      fea333555,fea333333,fea244555,fea244455,fea233445,fea233444,&
      fea233345,fea233344,fea233335,fea223355,fea222335,fea222334,&
      fea222333,fea222255,fea222245,fea222233,fea222224,fea145555,&
      fea134444,fea133444,fea133345,fea133334,fea133333,fea124555,&
      fea124455,fea123455,fea123345,fea113555,fea113345,fea112355,&
      fea112335,fea112233,fea111444,fea111234,fea111233,fea111123
      !
   real(ark)            ::  s1,s2,s3,s4,s5,tau
      !
   real(ark)            ::  &      
      ve  ,            &
      f0a1,f1a,f2a,f3a,f4a,f5a,f6a,f7a,f8a, &
      f1a1,f2a1,f3a1,f4a1,f5a1,f6a1,  &
      f0a11,f1a11,f2a11,f3a11,f4a11, &
      f0a12,f1a12,f2a12,f3a12,f4a12, &
      f0a14,f1a14,f2a14,f3a14,f4a14, &
      f0a44,f1a44,f2a44,f3a44,f4a44, &
      f0a111,f1a111,f2a111,f3a111  , &
      f0a112,f1a112,f2a112,f3a112  , &
      f0a114,f1a114,f2a114,f3a114  , &
      f0a123,f1a123,f2a123,f3a123  , &
      f0a124,f1a124,f2a124,f3a124  , &
      f0a144,f1a144,f2a144,f3a144  , &
      f0a155,f1a155,f2a155,f3a155  , &
      f0a455,f1a455,f2a455,f3a455  , &
      f0a1111,f1a1111,f2a1111      , &
      f0a1112,f1a1112,f2a1112      , &
      f0a1114,f1a1114,f2a1114      , &
      f0a1122,f1a1122,f2a1122      , &
      f0a1123,f1a1123,f2a1123      , &
      f0a1124,f1a1124,f2a1124      , &
      f0a1125,f1a1125,f2a1125      , &
      f0a1144,f1a1144,f2a1144      , &
      f0a1155,f1a1155,f2a1155      , &
      f0a1244,f1a1244,f2a1244      , &
      f0a1255,f1a1255,f2a1255      , &
      f0a1444,f1a1444,f2a1444      , &
      f0a1455,f1a1455,f2a1455      , &
      f0a4444,f1a4444,f2a4444      , &
      f0a44444 ,f1a44444 

   real(ark)            ::  &      
      f2a44444 ,f0a33455 ,f1a33455 ,f2a33455 ,f0a33445 ,f1a33445 ,&
      f2a33445 ,f0a33345 ,f1a33345 ,f2a33345 ,f0a33344 ,f1a33344 ,&
      f2a33344 ,f0a33334 ,f1a33334 ,f2a33334 ,f0a33333 ,f1a33333 ,&
      f2a33333 ,f0a25555 ,f1a25555 ,f2a25555 ,f0a24455 ,f1a24455 ,&
      f2a24455 ,f0a24445 ,f1a24445 ,f2a24445 ,f0a23333 ,f1a23333 ,&
      f2a23333 ,f0a13455 ,f1a13455 ,f2a13455 ,f0a13445 ,f1a13445 ,&
      f2a13445 ,f0a13345 ,f1a13345 ,f2a13345 ,f0a12355 ,f1a12355 ,&
      f2a12355 ,f0a11334 ,f1a11334 ,f2a11334 ,f0a11333 ,f1a11333 ,&
      f2a11333 ,f0a11255 ,f1a11255 ,f2a11255 ,f0a11245 ,f1a11245 ,&
      f2a11245 ,f0a11234 ,f1a11234 ,f2a11234 ,f0a11233 ,f1a11233 ,&
      f2a11233 ,f0a11135 ,f1a11135 ,f2a11135 ,f0a11134 ,f1a11134 ,&
      f2a11134 ,f0a11123 ,f1a11123 ,f2a11123 ,f0a555555,f1a555555


   real(ark)            ::  &      
      f2a555555,f0a444444,f1a444444,f2a444444,f0a335555,f1a335555,&
      f2a335555,f0a334455,f1a334455,f2a334455,f0a334445,f1a334445,&
      f2a334445,f0a333555,f1a333555,f2a333555,f0a333333,f1a333333,&
      f2a333333,f0a244555,f1a244555,f2a244555,f0a244455,f1a244455,&
      f2a244455,f0a233445,f1a233445,f2a233445,f0a233444,f1a233444,&
      f2a233444,f0a233345,f1a233345,f2a233345,f0a233344,f1a233344,&
      f2a233344,f0a233335,f1a233335,f2a233335,f0a223355,f1a223355,&
      f2a223355,f0a222335,f1a222335,f2a222335,f0a222334,f1a222334,&
      f2a222334,f0a222333,f1a222333,f2a222333,f0a222255,f1a222255,&
      f2a222255,f0a222245,f1a222245,f2a222245,f0a222233,f1a222233,&
      f2a222233,f0a222224,f1a222224,f2a222224,f0a145555,f1a145555,&
      f2a145555,f0a134444,f1a134444,f2a134444,f0a133444,f1a133444,&
      f2a133444,f0a133345,f1a133345,f2a133345,f0a133334,f1a133334,&
      f2a133334,f0a133333,f1a133333,f2a133333,f0a124555,f1a124555,&
      f2a124555,f0a124455,f1a124455,f2a124455,f0a123455,f1a123455,&
      f2a123455,f0a123345,f1a123345,f2a123345,f0a113555,f1a113555,&
      f2a113555,f0a113345,f1a113345,f2a113345,f0a112355,f1a112355,&
      f2a112355,f0a112335,f1a112335,f2a112335,f0a112233,f1a112233,&
      f2a112233,f0a111444,f1a111444,f2a111444,f0a111234,f1a111234,&
      f2a111234,f0a111233,f1a111233,f2a111233,f0a111123,f1a111123,&
      f2a111123


      N = size(force)
      !
      ve         = force(  1)
      f1a        = force(  2)
      f2a        = force(  3)
      f3a        = force(  4)
      f4a        = force(  5)
      f5a        = force(  6)
      f6a        = force(  7)
      f7a        = force(  8)
      f0a1       = force(  9)
      f1a1       = force( 10)
      f2a1       = force( 11)
      f3a1       = force( 12)
      f4a1       = force( 13)
      f5a1       = force( 14)
      f6a1       = force( 15)
      f0a11      = force( 16)
      f1a11      = force( 17)
      f2a11      = force( 18)
      f3a11      = force( 19)
      f4a11      = force( 20)
      f0a12      = force( 21)
      f1a12      = force( 22)
      f2a12      = force( 23)
      f3a12      = force( 24)
      f4a12      = force( 25)
      f0a14      = force( 26)
      f1a14      = force( 27)
      f2a14      = force( 28)
      f3a14      = force( 29)
      f4a14      = force( 30)
      f0a44      = force( 31)
      f1a44      = force( 32)
      f2a44      = force( 33)
      f3a44      = force( 34)
      f4a44      = force( 35)

      v0=ve+f1a*coro+f2a*coro**2+f3a*coro**3+f4a*coro**4+f5a*coro**5 &
            +f6a*coro**6+f7a*coro**7  !  +f8a*coro**8

      fea1= f0a1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4+f5a1*coro**5+f6a1*coro**6
      !
      fea11=   f0a11+f1a11*coro+f2a11*coro**2+f3a11*coro**3+f4a11*coro**4
      fea12=   f0a12+f1a12*coro+f2a12*coro**2+f3a12*coro**3+f4a12*coro**4
      fea14=   f0a14+f1a14*coro+f2a14*coro**2+f3a14*coro**3+f4a14*coro**4
      fea44=   f0a44+f1a44*coro+f2a44*coro**2+f3a44*coro**3+f4a44*coro**4
      !
      !
      v1 = (y3+y2+y1)*fea1
      !
      v2 = (y2*y3+y1*y3+y1*y2)*fea12                                                                 &
       +(y2**2+y3**2+y1**2)*fea11                                                                    &
       +(-sqrt(3.0_ark)*y3*y5/2.0_ark-y3*y4/2.0_ark+y1*y4+sqrt(3.0_ark)*y2*y5/2.0_ark-y2*y4/2.0_ark)*fea14 &
       +(y5**2+y4**2)*fea44
      !
      v3 = 0 ; v4 = 0 ; v5 = 0 ; v6 = 0 
      !
      if (N>35) then  
        f0a111     = force( 36)
        f1a111     = force( 37)
        f2a111     = force( 38)
        f3a111     = force( 39)
        f0a112     = force( 40)
        f1a112     = force( 41)
        f2a112     = force( 42)
        f3a112     = force( 43)
        f0a114     = force( 44)
        f1a114     = force( 45)
        f2a114     = force( 46)
        f3a114     = force( 47)
        f0a123     = force( 48)
        f1a123     = force( 49)
        f2a123     = force( 50)
        f3a123     = force( 51)
        f0a124     = force( 52)
        f1a124     = force( 53)
        f2a124     = force( 54)
        f3a124     = force( 55)
        f0a144     = force( 56)
        f1a144     = force( 57)
        f2a144     = force( 58)
        f3a144     = force( 59)
        f0a155     = force( 60)
        f1a155     = force( 61)
        f2a155     = force( 62)
        f3a155     = force( 63)
        f0a455     = force( 64)
        f1a455     = force( 65)
        f2a455     = force( 66)
        f3a455     = force( 67)
        !
        fea111= f0a111+f1a111*coro+f2a111*coro**2+f3a111*coro**3
        fea112= f0a112+f1a112*coro+f2a112*coro**2+f3a112*coro**3
        fea114= f0a114+f1a114*coro+f2a114*coro**2+f3a114*coro**3
        fea123= f0a123+f1a123*coro+f2a123*coro**2+f3a123*coro**3
        fea124= f0a124+f1a124*coro+f2a124*coro**2+f3a124*coro**3
        fea144= f0a144+f1a144*coro+f2a144*coro**2+f3a144*coro**3
        fea155= f0a155+f1a155*coro+f2a155*coro**2+f3a155*coro**3
        fea455= f0a455+f1a455*coro+f2a455*coro**2+f3a455*coro**3
        !
        v3 = (y1*y3*y4+y1*y2*y4-2.0_ark*y2*y3*y4+sqrt(3.0_ark)*y1*y2*y5-sqrt(3.0_ark)*y1*y3*y5)*fea124   &
         +(3.0_ark/4.0_ark*y3*y4**2-sqrt(3.0_ark)*y3*y4*y5/2.0_ark+y1*y5**2+y2*y5**2/4.0_ark               & 
         +3.0_ark/4.0_ark*y2*y4**2+sqrt(3.0_ark)*y2*y4*y5/2.0_ark+y3*y5**2/4.0_ark)*fea155                 &
         +(y2*y3**2+y1*y3**2+y1**2*y3+y1*y2**2+y2**2*y3+y1**2*y2)*fea112+                             &
         (-y4**3/3.0_ark+y4*y5**2)*fea455+fea123*y1*y2*y3                                              &
         +(y1*y4**2+3.0_ark/4.0_ark*y3*y5**2+3.0_ark/4.0_ark*y2*y5**2+y2*y4**2/4.0_ark                     &
         -sqrt(3.0_ark)*y2*y4*y5/2.0_ark+sqrt(3.0_ark)*y3*y4*y5/2.0_ark+y3*y4**2/4.0_ark)*fea144           &
         +(y3**3+y2**3+y1**3)*fea111+(-y2**2*y4/2.0_ark-y3**2*y4/2.0_ark+sqrt(3.0_ark)*y2**2*y5/2.0_ark   & 
         +y1**2*y4-sqrt(3.0_ark)*y3**2*y5/2.0_ark)*fea114
         !
      endif
      !
      if (N>67) then  
        !
        f0a1111    = force( 68)
        f1a1111    = force( 69)
        f2a1111    = force( 70)
        f0a1112    = force( 71)
        f1a1112    = force( 72)
        f2a1112    = force( 73)
        f0a1114    = force( 74)
        f1a1114    = force( 75)
        f2a1114    = force( 76)
        f0a1122    = force( 77)
        f1a1122    = force( 78)
        f2a1122    = force( 79)
        f0a1123    = force( 80)
        f1a1123    = force( 81)
        f2a1123    = force( 82)
        f0a1124    = force( 83)
        f1a1124    = force( 84)
        f2a1124    = force( 85)
        f0a1125    = force( 86)
        f1a1125    = force( 87)
        f2a1125    = force( 88)
        f0a1144    = force( 89)
        f1a1144    = force( 90)
        f2a1144    = force( 91)
        f0a1155    = force( 92)
        f1a1155    = force( 93)
        f2a1155    = force( 94)
        f0a1244    = force( 95)
        f1a1244    = force( 96)
        f2a1244    = force( 97)
        f0a1255    = force( 98)
        f1a1255    = force( 99)
        f2a1255    = force(100)
        f0a1444    = force(101)
        f1a1444    = force(102)
        f2a1444    = force(103)
        f0a1455    = force(104)
        f1a1455    = force(105)
        f2a1455    = force(106)
        f0a4444    = force(107)
        f1a4444    = force(108)
        f2a4444    = force(109)
        !
        fea1111= f0a1111+f1a1111*coro+f2a1111*coro**2
        fea1112= f0a1112+f1a1112*coro+f2a1112*coro**2
        fea1114= f0a1114+f1a1114*coro+f2a1114*coro**2
        fea1122= f0a1122+f1a1122*coro+f2a1122*coro**2
        fea1123= f0a1123+f1a1123*coro+f2a1123*coro**2
        fea1124= f0a1124+f1a1124*coro+f2a1124*coro**2
        fea1125= f0a1125+f1a1125*coro+f2a1125*coro**2
        fea1144= f0a1144+f1a1144*coro+f2a1144*coro**2
        fea1155= f0a1155+f1a1155*coro+f2a1155*coro**2
        fea1244= f0a1244+f1a1244*coro+f2a1244*coro**2
        fea1255= f0a1255+f1a1255*coro+f2a1255*coro**2
        fea1444= f0a1444+f1a1444*coro+f2a1444*coro**2
        fea1455= f0a1455+f1a1455*coro+f2a1455*coro**2
        fea4444= f0a4444+f1a4444*coro+f2a4444*coro**2
        !
        s2 = (y4**4+y5**4+2.0_ark*y4**2*y5**2)*fea4444+(3.0_ark/8.0_ark*sqrt(3.0_ark)*&
         y2*y5**3-3.0_ark/8.0_ark*sqrt(3.0_ark)*y3*y4**2*y5-3.0_ark/8.0_ark*sqrt(3.0_ark)*y3*&
         y5**3-9.0_ark/8.0_ark*y2*y4*y5**2-y3*y4**3/8.0_ark-y2*y4**3/8.0_ark-9.0_ark/8.0_ark*&
         y3*y4*y5**2+y1*y4**3+3.0_ark/8.0_ark*sqrt(3.0_ark)*y2*y4**2*y5)*fea1444 &
         +(3.0_ark/4.0_ark*y2**2*y4**2+3.0_ark/4.0_ark*y3**2*y4**2+y1**2*y5**2+y3**2*y5**2/4.0_ark &
         -sqrt(3.0_ark)*y3**2*y4*y5/2.0_ark+sqrt(3.0_ark)*y2**2*y4*y5/2.0_ark+y2**2&
         *y5**2/4.0_ark)*fea1155 
         s1 = s2+(y3**2*y4**2/4.0_ark+3.0_ark/4.0_ark*y3**2*y5**2+y1**2*y4**2+y2**2*&
         y4**2/4.0_ark+sqrt(3.0_ark)*y3**2*y4*y5/2.0_ark-sqrt(3.0_ark)*y2**2*y4*y5/2.0_ark&
         +3.0_ark/4.0_ark*y2**2*y5**2)*fea1144+(y1**3*y4+sqrt(3.0_ark)*y2**3*y5/2.0_ark&
         -sqrt(3.0_ark)*y3**3*y5/2.0_ark-y2**3*y4/2.0_ark-y3**3*y4/2.0_ark)*fea1114&
         +(y2**4+y1**4+y3**4)*fea1111+(sqrt(3.0_ark)*y1*y3*y4*y5+3.0_ark/2.0_ark*y2*y3*y5**2&
         -y2*y3*y4**2/2.0_ark+y1*y2*y4**2-sqrt(3.0_ark)*y1*y2*y4*y5+y1*y3*y4**2)*fea1244 
         !
        s2 = s1+(y1*y3*y5**2+y1*y2*y5**2-sqrt(3.0_ark)*y1*y3*y4*y5-y2*y3*y5**&
         2/2.0_ark+3.0_ark/2.0_ark*y2*y3*y4**2+sqrt(3.0_ark)*y1*y2*y4*y5)*fea1255+&
         (-y1*y3**2*y4/2.0_ark+y1**2*y3*y4-sqrt(3.0_ark)*y1*y3**2*y5/2.0_ark-sqrt(3.0_ark)*y2&
         *y3**2*y5/2.0_ark+y1**2*y2*y4+sqrt(3.0_ark)*y2**2*y3*y5/2.0_ark-y2**2*y3*y4&
         /2.0_ark+sqrt(3.0_ark)*y1*y2**2*y5/2.0_ark-y2*y3**2*y4/2.0_ark-y1*y2**2*y4/2.0_ark&
         )*fea1124+(y1**2*y2*y5+sqrt(3.0_ark)*y1*y3**2*y4/2.0_ark+sqrt(3.0_ark)*y1*&
         y2**2*y4/2.0_ark-sqrt(3.0_ark)*y2*y3**2*y4/2.0_ark-sqrt(3.0_ark)*y2**2*y3*y4/2.0_ark&
         -y2**2*y3*y5/2.0_ark+y2*y3**2*y5/2.0_ark-y1*y3**2*y5/2.0_ark+y1*y2**2*y5&
         /2.0_ark-y1**2*y3*y5)*fea1125 
         !
        v4 = s2+(y2*y3**3+y1**3*y3+y1**3*y2+y1*y2**3+y1*y3**3+y2**3*y3)*fea1112+&
         (y2**2*y3**2+y1**2*y3**2+y1**2*y2**2)*fea1122+(y1*y2**2*y3&
         +y1**2*y2*y3+y1*y2*y3**2)*fea1123+(5.0_ark/8.0_ark*y2*y4*y5**2+sqrt(3.0_ark)*&
         y2*y5**3/8.0_ark-sqrt(3.0_ark)*y3*y4**2*y5/8.0_ark+sqrt(3.0_ark)*y2*y4**2*y5/8.0_ark&
         -3.0_ark/8.0_ark*y2*y4**3+y1*y4*y5**2-sqrt(3.0_ark)*y3*y5**3/8.0_ark&
         +5.0_ark/8.0_ark*y3*y4*y5**2-3.0_ark/8.0_ark*y3*y4**3)*fea1455
        !
      endif
      !
      if (N>109) then  
        !
        f0a44444   = force(110)
        f1a44444   = force(111)
        f2a44444   = force(112)
        f0a33455   = force(113)
        f1a33455   = force(114)
        f2a33455   = force(115)
        f0a33445   = force(116)
        f1a33445   = force(117)
        f2a33445   = force(118)
        f0a33345   = force(119)
        f1a33345   = force(120)
        f2a33345   = force(121)
        f0a33344   = force(122)
        f1a33344   = force(123)
        f2a33344   = force(124)
        f0a33334   = force(125)
        f1a33334   = force(126)
        f2a33334   = force(127)
        f0a33333   = force(128)
        f1a33333   = force(129)
        f2a33333   = force(130)
        f0a25555   = force(131)
        f1a25555   = force(132)
        f2a25555   = force(133)
        f0a24455   = force(134)
        f1a24455   = force(135)
        f2a24455   = force(136)
        f0a24445   = force(137)
        f1a24445   = force(138)
        f2a24445   = force(139)
        f0a23333   = force(140)
        f1a23333   = force(141)
        f2a23333   = force(142)
        f0a13455   = force(143)
        f1a13455   = force(144)
        f2a13455   = force(145)
        f0a13445   = force(146)
        f1a13445   = force(147)
        f2a13445   = force(148)
        f0a13345   = force(149)
        f1a13345   = force(150)
        f2a13345   = force(151)
        f0a12355   = force(152)
        f1a12355   = force(153)
        f2a12355   = force(154)
        f0a11334   = force(155)
        f1a11334   = force(156)
        f2a11334   = force(157)
        f0a11333   = force(158)
        f1a11333   = force(159)
        f2a11333   = force(160)
        f0a11255   = force(161)
        f1a11255   = force(162)
        f2a11255   = force(163)
        f0a11245   = force(164)
        f1a11245   = force(165)
        f2a11245   = force(166)
        f0a11234   = force(167)
        f1a11234   = force(168)
        f2a11234   = force(169)
        f0a11233   = force(170)
        f1a11233   = force(171)
        f2a11233   = force(172)
        f0a11135   = force(173)
        f1a11135   = force(174)
        f2a11135   = force(175)
        f0a11134   = force(176)
        f1a11134   = force(177)
        f2a11134   = force(178)
        f0a11123   = force(179)
        f1a11123   = force(180)
        f2a11123   = force(181)

        fea44444 = f0a44444  + f1a44444 *coro+ f2a44444 *coro**2
        fea33455 = f0a33455  + f1a33455 *coro+ f2a33455 *coro**2
        fea33445 = f0a33445  + f1a33445 *coro+ f2a33445 *coro**2
        fea33345 = f0a33345  + f1a33345 *coro+ f2a33345 *coro**2
        fea33344 = f0a33344  + f1a33344 *coro+ f2a33344 *coro**2
        fea33334 = f0a33334  + f1a33334 *coro+ f2a33334 *coro**2
        fea33333 = f0a33333  + f1a33333 *coro+ f2a33333 *coro**2
        fea25555 = f0a25555  + f1a25555 *coro+ f2a25555 *coro**2
        fea24455 = f0a24455  + f1a24455 *coro+ f2a24455 *coro**2
        fea24445 = f0a24445  + f1a24445 *coro+ f2a24445 *coro**2
        fea23333 = f0a23333  + f1a23333 *coro+ f2a23333 *coro**2
        fea13455 = f0a13455  + f1a13455 *coro+ f2a13455 *coro**2
        fea13445 = f0a13445  + f1a13445 *coro+ f2a13445 *coro**2
        fea13345 = f0a13345  + f1a13345 *coro+ f2a13345 *coro**2
        fea12355 = f0a12355  + f1a12355 *coro+ f2a12355 *coro**2
        fea11334 = f0a11334  + f1a11334 *coro+ f2a11334 *coro**2
        fea11333 = f0a11333  + f1a11333 *coro+ f2a11333 *coro**2
        fea11255 = f0a11255  + f1a11255 *coro+ f2a11255 *coro**2
        fea11245 = f0a11245  + f1a11245 *coro+ f2a11245 *coro**2
        fea11234 = f0a11234  + f1a11234 *coro+ f2a11234 *coro**2
        fea11233 = f0a11233  + f1a11233 *coro+ f2a11233 *coro**2
        fea11135 = f0a11135  + f1a11135 *coro+ f2a11135 *coro**2
        fea11134 = f0a11134  + f1a11134 *coro+ f2a11134 *coro**2
        fea11123 = f0a11123  + f1a11123 *coro+ f2a11123 *coro**2
        !    
        s3 = (y4**5-2.0_ark*y4**3*y5**2-3.0_ark*y4*y5**4)*fea44444+(-4.0_ark*y3*y4*&
        y5**3*sqrt(3.0_ark)+9.0_ark*y1*y4**2*y5**2-3.0_ark/2.0_ark*y1*y4**4+4.0_ark*y2*y4&
        *y5**3*sqrt(3.0_ark)+3.0_ark*y2*y4**4+5.0_ark/2.0_ark*y1*y5**4+3.0_ark*y3*y4**4+&
        y2*y5**4+y3*y5**4)*fea25555+(-y2*y4**4+y3*y4**2*y5**2-2.0_ark*y2*y4*y5&
        **3*sqrt(3.0_ark)-y3*y4**4-7.0_ark/2.0_ark*y1*y4**2*y5**2-3.0_ark/4.0_ark*y1*y5**4&
        +2.0_ark*y3*y4*y5**3*sqrt(3.0_ark)+y2*y4**2*y5**2+5.0_ark/4.0_ark*y1*y4**4)*fea24455 
        !
        s2 = s3+(y2*y4**3*y5-3.0_ark*y3*y4*y5**3+2.0_ark/3.0_ark*y3*y4**4*sqrt(3.0_ark&
        )+3.0_ark/4.0_ark*y1*y5**4*sqrt(3.0_ark)+3.0_ark*y2*y4*y5**3-&
        7.0_ark/12.0_ark*y1*y4**4*sqrt(3.0_ark)+3.0_ark/2.0_ark*y1*y4**2*y5**2*sqrt(3.0_ark)-y3*y4**3*y5&
        +2.0_ark/3.0_ark*y2*y4**4*sqrt(3.0_ark))*fea24445+(-y2**2*y5**3+y3**2*y4**2*y5+ &
        y3**2*y5**3+4.0_ark/9.0_ark*y2**2*y4**3*sqrt(3.0_ark)-5.0_ark/9.0_ark*y1**2*y4**3*&
        sqrt(3.0_ark)+4.0_ark/9.0_ark*y3**2*y4**3*sqrt(3.0_ark)-y2**2*y4**2*y5&
        -y1**2*y4*y5**2*sqrt(3.0_ark))*fea33445+(y3**2*y4*y5**2-y1**2*y4**3/3.0_ark&
        -y3**2*y4**3/3.0_ark+y1**2*y4*y5**2+y2**2*y4*y5**2-y2**2*y4**3/3.0_ark)*fea33455 
        !
        s1 = s2+(-y2**3*y4*y5+y3**3*y4*y5+y2**3*y5**2*sqrt(3.0_ark)/3.0_ark+y1**&
        3*y4**2*sqrt(3.0_ark)/2.0_ark+y3**3*y5**2*sqrt(3.0_ark)/3.0_ark- & 
        y1**3*y5**2*sqrt(3.0_ark)/6.0_ark)*fea33345+(y3**3*y4**2+y3**3*y5**2+y2**3*y4**2+y2**3&
        *y5**2+y1**3*y5**2+y1**3*y4**2)*fea33344+(y3**4*y4+sqrt(3.0_ark)*y3**&
        4*y5+y2**4*y4-2.0_ark*y1**4*y4-sqrt(3.0_ark)*y2**4*y5)*fea33334+(y2**5+ &
        y3**5+y1**5)*fea33333+(-4.0_ark/9.0_ark*y1*y2*y4**3*sqrt(3.0_ark)-y1*y2*y5**3+ &
        y1*y3*y4**2*y5+y2*y3*y4*y5**2*sqrt(3.0_ark)-y1*y2*y4**2*y5+5.0_ark/9.0_ark&
        *y2*y3*y4**3*sqrt(3.0_ark)-4.0_ark/9.0_ark*y1*y3*y4**3*sqrt(3.0_ark)+y1*y3*y5&
        **3)*fea13445+(y2*y3*y4*y5**2+y1*y2*y4*y5**2-y2*y3*y4**3/3.0_ark- & 
        y1*y2*y4**3/3.0_ark-y1*y3*y4**3/3.0_ark+y1*y3*y4*y5**2)*fea13455 
        
        s3 = s1+(y1**2*y3*y5**2+y2**2*y3*y4**2+y2**2*y3*y5**2+y1*y2**2*y5**2+&
        y1**2*y2*y5**2+y1*y2**2*y4**2+y2*y3**2*y4**2+y1*y3**2*y4**2+&
        y1**2*y3*y4**2+y1**2*y2*y4**2+y1*y3**2*y5**2+y2*y3**2*y5**2)*fea11255&
        +(2.0_ark/3.0_ark*y1**2*y3*y4**2*sqrt(3.0_ark)+y1*y3**2*y5**2*sqrt(3.0_ark)/2.0_ark+&
        y1*y2**2*y5**2*sqrt(3.0_ark)/2.0_ark+y2**2*y3*y5**2*sqrt(3.0_ark)/2.0_ark-&
        y1*y2**2*y4*y5+y2*y3**2*y4*y5+y1*y3**2*y4*y5-y2**2*y3*y4*y5+y2*y3**&
        2*y4**2*sqrt(3.0_ark)/6.0_ark+y1*y3**2*y4**2*sqrt(3.0_ark)/6.0_ark+y1*y2**2*y4&
        **2*sqrt(3.0_ark)/6.0_ark+2.0_ark/3.0_ark*y1**2*y2*y4**2*sqrt(3.0_ark)+&
        y2*y3**2*y5**2*sqrt(3.0_ark)/2.0_ark+y2**2*y3*y4**2*sqrt(3.0_ark)/6.0_ark)*fea13345 
        s4 = s3+(y1**2*y2*y4*y5+y1**2*y3*y4**2*sqrt(3.0_ark)/3.0_ark+y1**2*y2*y4&
        **2*sqrt(3.0_ark)/3.0_ark-y1*y2**2*y4**2*sqrt(3.0_ark)/6.0_ark+y2*y3**2*y4*y5-&
        y2**2*y3*y4*y5-y1**2*y3*y4*y5+y2*y3**2*y4**2*sqrt(3.0_ark)/3.0_ark+y1*y2&
        **2*y5**2*sqrt(3.0_ark)/2.0_ark-y1*y3**2*y4**2*sqrt(3.0_ark)/6.0_ark+y2**2*y3*&
        y4**2*sqrt(3.0_ark)/3.0_ark+y1*y3**2*y5**2*sqrt(3.0_ark)/2.0_ark)*fea11245 
        s2 = s4+(-y1**3*y2*y5+y1**3*y3*y5+y2**3*y3*y5/2.0_ark-y1*y2**3*y4*sqrt(3.0_ark)/2.0_ark-&
        y1*y2**3*y5/2.0_ark-y2*y3**3*y5/2.0_ark+y1*y3**3*y5/2.0_ark+y2&
        **3*y3*y4*sqrt(3.0_ark)/2.0_ark+y2*y3**3*y4*sqrt(3.0_ark)/2.0_ark-y1*y3**3*y4*&
        sqrt(3.0_ark)/2.0_ark)*fea11135+(y1**3*y3*y4-y2**3*y3*y4/2.0_ark+y1**3*y2*y4-&
        y2*y3**3*y4/2.0_ark-y1*y3**3*y4/2.0_ark+y1*y2**3*y5*sqrt(3.0_ark)/2.0_ark+y2&
        **3*y3*y5*sqrt(3.0_ark)/2.0_ark-y2*y3**3*y5*sqrt(3.0_ark)/2.0_ark-y1*y2**3*y4/&
        2.0_ark-y1*y3**3*y5*sqrt(3.0_ark)/2.0_ark)*fea11134 
        
        v5 = s2+(y1*y2**4+y1**4*y3+y1**4*y2+y2**4*y3+y2*y3**4+y1*y3**4)*fea23333+&
        (-2.0_ark*y2**2*y3**2*y4+y1**2*y2**2*y4-sqrt(3.0_ark)*y1**2*y3**2&
        *y5+sqrt(3.0_ark)*y1**2*y2**2*y5+y1**2*y3**2*y4)*fea11334+(y1**2*y3**&
        3+y1**3*y3**2+y2**2*y3**3+y1**2*y2**3+y1**3*y2**2+y2**3*y3**2)*fea11333+&
        (y1*y2*y3*y4**2+y1*y2*y3*y5**2)*fea12355+(-y1*y2*y3**2*y4/2.0_ark-&
        y1*y2**2*y3*y4/2.0_ark-sqrt(3.0_ark)*y1*y2*y3**2*y5/2.0_ark+y1**2*y2*y3*y4+&
        sqrt(3.0_ark)*y1*y2**2*y3*y5/2.0_ark)*fea11234+(y1*y2**3*y3+y1*y2*y3**3+&
        y1**3*y2*y3)*fea11123+(y1**2*y2**2*y3+y1*y2**2*y3**2+y1**2*y2*y3**2)*fea11233
        !
      endif
      !
      if (N>181) then  
        !
        f0a555555  = force(182)
        f1a555555  = force(183)
        f2a555555  = force(184)
        f0a444444  = force(185)
        f1a444444  = force(186)
        f2a444444  = force(187)
        f0a335555  = force(188)
        f1a335555  = force(189)
        f2a335555  = force(190)
        f0a334455  = force(191)
        f1a334455  = force(192)
        f2a334455  = force(193)
        f0a334445  = force(194)
        f1a334445  = force(195)
        f2a334445  = force(196)
        f0a333555  = force(197)
        f1a333555  = force(198)
        f2a333555  = force(199)
        f0a333333  = force(200)
        f1a333333  = force(201)
        f2a333333  = force(202)
        f0a244555  = force(203)
        f1a244555  = force(204)
        f2a244555  = force(205)
        f0a244455  = force(206)
        f1a244455  = force(207)
        f2a244455  = force(208)
        f0a233445  = force(209)
        f1a233445  = force(210)
        f2a233445  = force(211)
        f0a233444  = force(212)
        f1a233444  = force(213)
        f2a233444  = force(214)
        f0a233345  = force(215)
        f1a233345  = force(216)
        f2a233345  = force(217)
        f0a233344  = force(218)
        f1a233344  = force(219)
        f2a233344  = force(220)
        f0a233335  = force(221)
        f1a233335  = force(222)
        f2a233335  = force(223)
        f0a223355  = force(224)
        f1a223355  = force(225)
        f2a223355  = force(226)
        f0a222335  = force(227)
        f1a222335  = force(228)
        f2a222335  = force(229)
        f0a222334  = force(230)
        f1a222334  = force(231)
        f2a222334  = force(232)
        f0a222333  = force(233)
        f1a222333  = force(234)
        f2a222333  = force(235)
        f0a222255  = force(236)
        f1a222255  = force(237)
        f2a222255  = force(238)
        f0a222245  = force(239)
        f1a222245  = force(240)
        f2a222245  = force(241)
        f0a222233  = force(242)
        f1a222233  = force(243)
        f2a222233  = force(244)
        f0a222224  = force(245)
        f1a222224  = force(246)
        f2a222224  = force(247)
        f0a145555  = force(248)
        f1a145555  = force(249)
        f2a145555  = force(250)
        f0a134444  = force(251)
        f1a134444  = force(252)
        f2a134444  = force(253)
        f0a133444  = force(254)
        f1a133444  = force(255)
        f2a133444  = force(256)
        f0a133345  = force(257)
        f1a133345  = force(258)
        f2a133345  = force(259)
        f0a133334  = force(260)
        f1a133334  = force(261)
        f2a133334  = force(262)
        f0a133333  = force(263)
        f1a133333  = force(264)
        f2a133333  = force(265)
        f0a124555  = force(266)
        f1a124555  = force(267)
        f2a124555  = force(268)
        f0a124455  = force(269)
        f1a124455  = force(270)
        f2a124455  = force(271)
        f0a123455  = force(272)
        f1a123455  = force(273)
        f2a123455  = force(274)
        f0a123345  = force(275)
        f1a123345  = force(276)
        f2a123345  = force(277)
        f0a113555  = force(278)
        f1a113555  = force(279)
        f2a113555  = force(280)
        f0a113345  = force(281)
        f1a113345  = force(282)
        f2a113345  = force(283)
        f0a112355  = force(284)
        f1a112355  = force(285)
        f2a112355  = force(286)
        f0a112335  = force(287)
        f1a112335  = force(288)
        f2a112335  = force(289)
        f0a112233  = force(290)
        f1a112233  = force(291)
        f2a112233  = force(292)
        f0a111444  = force(293)
        f1a111444  = force(294)
        f2a111444  = force(295)
        f0a111234  = force(296)
        f1a111234  = force(297)
        f2a111234  = force(298)
        f0a111233  = force(299)
        f1a111233  = force(300)
        f2a111233  = force(301)
        f0a111123  = force(302)
        f1a111123  = force(303)
        f2a111123  = force(304)

        fea555555= f0a555555 + f1a555555*coro+ f2a555555*coro**2
        fea444444= f0a444444 + f1a444444*coro+ f2a444444*coro**2
        fea335555= f0a335555 + f1a335555*coro+ f2a335555*coro**2
        fea334455= f0a334455 + f1a334455*coro+ f2a334455*coro**2
        fea334445= f0a334445 + f1a334445*coro+ f2a334445*coro**2
        fea333555= f0a333555 + f1a333555*coro+ f2a333555*coro**2
        fea333333= f0a333333 + f1a333333*coro+ f2a333333*coro**2
        fea244555= f0a244555 + f1a244555*coro+ f2a244555*coro**2
        fea244455= f0a244455 + f1a244455*coro+ f2a244455*coro**2
        fea233445= f0a233445 + f1a233445*coro+ f2a233445*coro**2
        fea233444= f0a233444 + f1a233444*coro+ f2a233444*coro**2
        fea233345= f0a233345 + f1a233345*coro+ f2a233345*coro**2
        fea233344= f0a233344 + f1a233344*coro+ f2a233344*coro**2
        fea233335= f0a233335 + f1a233335*coro+ f2a233335*coro**2
        fea223355= f0a223355 + f1a223355*coro+ f2a223355*coro**2
        fea222335= f0a222335 + f1a222335*coro+ f2a222335*coro**2
        fea222334= f0a222334 + f1a222334*coro+ f2a222334*coro**2
        fea222333= f0a222333 + f1a222333*coro+ f2a222333*coro**2
        fea222255= f0a222255 + f1a222255*coro+ f2a222255*coro**2
        fea222245= f0a222245 + f1a222245*coro+ f2a222245*coro**2
        fea222233= f0a222233 + f1a222233*coro+ f2a222233*coro**2
        fea222224= f0a222224 + f1a222224*coro+ f2a222224*coro**2
        fea145555= f0a145555 + f1a145555*coro+ f2a145555*coro**2
        fea134444= f0a134444 + f1a134444*coro+ f2a134444*coro**2
        fea133444= f0a133444 + f1a133444*coro+ f2a133444*coro**2
        fea133345= f0a133345 + f1a133345*coro+ f2a133345*coro**2
        fea133334= f0a133334 + f1a133334*coro+ f2a133334*coro**2
        fea133333= f0a133333 + f1a133333*coro+ f2a133333*coro**2
        fea124555= f0a124555 + f1a124555*coro+ f2a124555*coro**2
        fea124455= f0a124455 + f1a124455*coro+ f2a124455*coro**2
        fea123455= f0a123455 + f1a123455*coro+ f2a123455*coro**2
        fea123345= f0a123345 + f1a123345*coro+ f2a123345*coro**2
        fea113555= f0a113555 + f1a113555*coro+ f2a113555*coro**2
        fea113345= f0a113345 + f1a113345*coro+ f2a113345*coro**2
        fea112355= f0a112355 + f1a112355*coro+ f2a112355*coro**2
        fea112335= f0a112335 + f1a112335*coro+ f2a112335*coro**2
        fea112233= f0a112233 + f1a112233*coro+ f2a112233*coro**2
        fea111444= f0a111444 + f1a111444*coro+ f2a111444*coro**2
        fea111234= f0a111234 + f1a111234*coro+ f2a111234*coro**2
        fea111233= f0a111233 + f1a111233*coro+ f2a111233*coro**2
        fea111123= f0a111123 + f1a111123*coro+ f2a111123*coro**2
        !
        s3 = (y2**3*y4**3*sqrt(3.0_ark)-y2**3*y4**2*y5+y3**3*y4**2*y5-&
        5.0_ark/3.0_ark*y2**3*y4*y5**2*sqrt(3.0_ark)+y3**3*y4**3*sqrt(3.0_ark)-5.0_ark/3.0_ark*y3**&
        3*y4*y5**2*sqrt(3.0_ark)-y2**3*y5**3+y3**3*y5**3-8.0_ark/3.0_ark*y1**3*y4*y5**2*sqrt(3.0_ark))*fea333555+&
        (y1**4*y5**2*sqrt(3.0_ark)/2.0_ark+y2**4*y4*y5+y2**4*y4**2*sqrt(3.0_ark)/3.0_ark+&
        y3**4*y4**2*sqrt(3.0_ark)/3.0_ark-y3**4*y4&
        *y5-y1**4*y4**2*sqrt(3.0_ark)/6.0_ark)*fea222245+(y1*y3**5+y1*y2**5+y2**&
        5*y3+y1**5*y3+y1**5*y2+y2*y3**5)*fea133333+(y1**4*y3*y4-2.0_ark*y2**4&
        *y3*y4+y1**4*y2*y4+y1*y2**4*y5*sqrt(3.0_ark)+y1*y3**4*y4-2.0_ark*y2*y3**&
        4*y4+y1**4*y2*y5*sqrt(3.0_ark)-y1*y3**4*y5*sqrt(3.0_ark)-y1**4*y3*y5*sqrt(3.0_ark)+&
        y1*y2**4*y4)*fea133334+(-y1*y2*y3*y4**3/3.0_ark+y1*y2*y3*y4*y5**2)*fea123455 
        
        s4 = s3+(2.0_ark/3.0_ark*sqrt(3.0_ark)*y1*y2**2*y3**2*y4-y1**2*y2**2*y3*y5-&
        sqrt(3.0_ark)*y1**2*y2**2*y3*y4/3.0_ark+y1**2*y2*y3**2*y5-&
        sqrt(3.0_ark)*y1**2*y2*y3**2*y4/3.0_ark)*fea112335+(y1*y2**2*y3*y5**2+y1*y2*y3**2*y5**2+&
        y1*y2*y3**2*y4**2+y1*y2**2*y3*y4**2+y1**2*y2*y3*y4**2+y1**2*y2*y3*y5**2)*fea112355 
        
        s2 = s4+(y2**3*y3**2*y5-y1**3*y2**2*y5/2.0_ark-y1**2*y3**3*y5/2.0_ark-y2&
        **2*y3**3*y5+y1**3*y2**2*y4*sqrt(3.0_ark)/2.0_ark-y1**2*y2**3*y4*sqrt(3.0_ark)/2.0_ark+&
        y1**3*y3**2*y5/2.0_ark+y1**2*y2**3*y5/2.0_ark+y1**3*y3**2*y4*sqrt(3.0_ark)/2.0_ark-&
        y1**2*y3**3*y4*sqrt(3.0_ark)/2.0_ark)*fea222335+(-y1**2*y2&
        **2*y5**2*sqrt(3.0_ark)/2.0_ark-y1**2*y3**2*y5**2*sqrt(3.0_ark)/2.0_ark-y1**2*&
        y2**2*y4**2*sqrt(3.0_ark)/6.0_ark-y1**2*y2**2*y4*y5-2.0_ark/3.0_ark*y2**2*y3**&
        2*y4**2*sqrt(3.0_ark)+y1**2*y3**2*y4*y5-y1**2*y3**2*y4**2*sqrt(3.0_ark)/&
        6.0_ark)*fea113345+(y2**2*y3**2*y5**2+y2**2*y3**2*y4**2+y1**2*y2**2*y5**2+&
        y1**2*y3**2*y4**2+y1**2*y3**2*y5**2+y1**2*y2**2*y4**2)*fea223355 
        
        s3 = s2+(y1*y2*y3**2*y4**2*sqrt(3.0_ark)/6.0_ark+y1*y2*y3**2*y4*y5+y1*y2&
        *y3**2*y5**2*sqrt(3.0_ark)/2.0_ark+2.0_ark/3.0_ark*y1**2*y2*y3*y4**2*sqrt(3.0_ark&
        )-y1*y2**2*y3*y4*y5+y1*y2**2*y3*y4**2*sqrt(3.0_ark)/6.0_ark+y1*y2**2*y3*&
        y5**2*sqrt(3.0_ark)/2.0_ark)*fea123345+(-y1**3*y2**2*y5*sqrt(3.0_ark)/2.0_ark-&
        y1**3*y2**2*y4/2.0_ark-y1**3*y3**2*y4/2.0_ark-y1**2*y2**3*y4/2.0_ark+y1**3*&
        y3**2*y5*sqrt(3.0_ark)/2.0_ark-y1**2*y3**3*y4/2.0_ark+y2**3*y3**2*y4-y1**2*&
        y2**3*y5*sqrt(3.0_ark)/2.0_ark+y2**2*y3**3*y4+y1**2*y3**3*y5*sqrt(3.0_ark)/&
        2.0_ark)*fea222334+(3.0_ark*y3**2*y4**4+5.0_ark/2.0_ark*y1**2*y5**4+y2**2*y5**&
        4+3.0_ark*y2**2*y4**4-4.0_ark*y3**2*y4*y5**3*sqrt(3.0_ark)+y3**2*y5**4+9.0_ark&
        *y1**2*y4**2*y5**2-3.0_ark/2.0_ark*y1**2*y4**4+4.0_ark*y2**2*y4*y5**3*sqrt(&
        3.0_ark))*fea335555+(y1**3*y2**3+y1**3*y3**3+y2**3*y3**3)*fea222333 
        
        s4 = s3+(y3*y4**5/5.0_ark-y2*y4**4*y5*sqrt(3.0_ark)/2.0_ark-2.0_ark/5.0_ark*y1*y4&
        **5-2.0_ark*y1*y4**3*y5**2-3.0_ark/10.0_ark*y2*y5**5*sqrt(3.0_ark)+y3*y4**3*y5&
        **2+y3*y4**4*y5*sqrt(3.0_ark)/2.0_ark+y2*y4**3*y5**2+3.0_ark/10.0_ark*y3*y5**5&
        *sqrt(3.0_ark)+y2*y4**5/5.0_ark)*fea244455+(y2**5*y4-2.0_ark*y1**5*y4-sqrt(&
        3.0_ark)*y2**5*y5+y3**5*y4+sqrt(3.0_ark)*y3**5*y5)*fea222224 
        
        s5 = s4+(-y3*y5**5*sqrt(3.0_ark)/5.0_ark+y2*y5**5*sqrt(3.0_ark)/5.0_ark+y1*y4*&
        y5**4-7.0_ark/15.0_ark*y2*y4**5+y2*y4**4*y5*sqrt(3.0_ark)/3.0_ark-y3*y4**4*y5*&
        sqrt(3.0_ark)/3.0_ark+y3*y4*y5**4+y2*y4*y5**4+2.0_ark*y1*y4**3*y5**2-7.0_ark/15.0_ark*y3*y4**5-&
        y1*y4**5/15.0_ark)*fea145555 
        
        s1 = s5+(-sqrt(3.0_ark)*y1*y2*y3**3*y5/2.0_ark+y1**3*y2*y3*y4+sqrt(3.0_ark)&
        *y1*y2**3*y3*y5/2.0_ark-y1*y2**3*y3*y4/2.0_ark-y1*y2*y3**3*y4/2.0_ark)*fea111234+&
        (y3*y4**4*y5/3.0_ark+y3*y4**5*sqrt(3.0_ark)/18.0_ark-y2*y4**4*y5/3.0_ark&
        -y2*y4*y5**4*sqrt(3.0_ark)/2.0_ark-y3*y4**2*y5**3+2.0_ark/9.0_ark*y1*y4**5*sqrt(3.0_ark)+&
        y2*y4**5*sqrt(3.0_ark)/18.0_ark+y2*y4**2*y5**3-2.0_ark/3.0_ark*y1*y4**&
        3*y5**2*sqrt(3.0_ark)-y3*y4*y5**4*sqrt(3.0_ark)/2.0_ark)*fea244555+(y1*y2*y4**2*y5**2-&
        3.0_ark/4.0_ark*y2*y3*y4**4-y1*y2*y5**4-y1*y3*y5**4+5.0_ark/4.0_ark*y2*y3*y5**4+&
        y1*y3*y4**2*y5**2-7.0_ark/2.0_ark*y2*y3*y4**2*y5**2-2.0_ark*y1&
        *y2*y4**3*y5*sqrt(3.0_ark)+2.0_ark*y1*y3*y4**3*y5*sqrt(3.0_ark))*fea124455 
        
        s3 = s1+(y2**6+y1**6+y3**6)*fea333333+(y1*y2**4*y3+y1**4*y2*y3+y1*&
        y2*y3**4)*fea111123+fea112233*y1**2*y2**2*y3**2+(y1**4*y4**2+y2**4&
        *y4**2+y2**4*y5**2+y3**4*y4**2+y1**4*y5**2+y3**4*y5**2)*fea222255 
        s4 = s3+(3.0_ark*y1*y3*y5**4+y1*y3*y4**4+9.0_ark*y2*y3*y4**2*y5**2-3.0_ark/&
        2.0_ark*y2*y3*y5**4-4.0_ark*y1*y3*y4**3*y5*sqrt(3.0_ark)+y1*y2*y4**4+&
        4.0_ark*y1*y2*y4**3*y5*sqrt(3.0_ark)+3.0_ark*y1*y2*y5**4+5.0_ark/2.0_ark*y2*y3*y4**4)*fea134444+&
        (-y1*y3**2*y5**3*sqrt(3.0_ark)/3.0_ark-7.0_ark/3.0_ark*y1**2*y3*y4*y5&
        **2+5.0_ark/3.0_ark*y1*y2**2*y4**2*y5*sqrt(3.0_ark)-13.0_ark/3.0_ark*y2**2*y3*y4*&
        y5**2-4.0_ark/3.0_ark*y2*y3**2*y5**3*sqrt(3.0_ark)-7.0_ark/3.0_ark*y1**2*y2*y4*y5&
        **2-16.0_ark/3.0_ark*y1*y3**2*y4*y5**2+4.0_ark/3.0_ark*y1**2*y3*y4**2*y5*sqrt(&
        3.0_ark)+4.0_ark/3.0_ark*y2**2*y3*y5**3*sqrt(3.0_ark)+3.0_ark*y1**2*y2*y4**3+&
        y2*y3**2*y4**3+y1*y2**2*y5**3*sqrt(3.0_ark)/3.0_ark+y2**2*y3*y4**3-13.0_ark/3.0_ark&
        *y2*y3**2*y4*y5**2-5.0_ark/3.0_ark*y1*y3**2*y4**2*y5*sqrt(3.0_ark)-&
        4.0_ark/3.0_ark*y1**2*y2*y4**2*y5*sqrt(3.0_ark)+3.0_ark*y1**2*y3*y4**3-16.0_ark/3.0_ark*y1*&
        y2**2*y4*y5**2)*fea233444 
        
        s5 = s4+(2.0_ark*y1*y3**2*y5**3+4.0_ark*y2*y3**2*y5**3+4.0_ark*y2**2*y3*y4*&
        y5**2*sqrt(3.0_ark)-2.0_ark*y1*y2**2*y5**3+y1**2*y3*y4*y5**2*sqrt(3.0_ark)+&
        6.0_ark*y1*y3**2*y4**2*y5-6.0_ark*y1*y2**2*y4**2*y5-3.0_ark*y1**2*y3*y4**2*&
        y5+y1**2*y2*y4*y5**2*sqrt(3.0_ark)+4.0_ark*y1*y3**2*y4*y5**2*sqrt(3.0_ark)-&
        3.0_ark*y1**2*y2*y4**3*sqrt(3.0_ark)-4.0_ark*y2**2*y3*y5**3+3.0_ark*y1**2*y2*y4**2*y5-&
        y1**2*y2*y5**3+y1**2*y3*y5**3-3.0_ark*y1**2*y3*y4**3*sqrt(3.0_ark&
        )+4.0_ark*y2*y3**2*y4*y5**2*sqrt(3.0_ark)+4.0_ark*y1*y2**2*y4*y5**2*sqrt(3.0_ark))*fea113555 
        
        s2 = s5+(-2.0_ark/3.0_ark*y3**2*y4**4*sqrt(3.0_ark)-3.0_ark/2.0_ark*y1**2*y4**2*y5**2*sqrt(3.0_ark)-&
        3.0_ark/4.0_ark*y1**2*y5**4*sqrt(3.0_ark)-y2**2*y4**3*y5+&
        7.0_ark/12.0_ark*y1**2*y4**4*sqrt(3.0_ark)+y3**2*y4**3*y5+3.0_ark*y3**2*y4*y5**3&
        -2.0_ark/3.0_ark*y2**2*y4**4*sqrt(3.0_ark)-3.0_ark*y2**2*y4*y5**3)*fea334445+(&
        -3.0_ark*y1*y3*y4**3*y5+2.0_ark/3.0_ark*y1*y2*y5**4*sqrt(3.0_ark)-y1*y3*y4*y5**3+&
        2.0_ark/3.0_ark*y1*y3*y5**4*sqrt(3.0_ark)+3.0_ark*y1*y2*y4**3*y5-7.0_ark/12.0_ark&
        *y2*y3*y5**4*sqrt(3.0_ark)+3.0_ark/2.0_ark*y2*y3*y4**2*y5**2*sqrt(3.0_ark)+y1*&
        y2*y4*y5**3+3.0_ark/4.0_ark*y2*y3*y4**4*sqrt(3.0_ark))*fea124555+(2.0_ark*y3**&
        2*y4*y5**3*sqrt(3.0_ark)-7.0_ark/2.0_ark*y1**2*y4**2*y5**2+y2**2*y4**2*y5**&
        2-y2**2*y4**4-y3**2*y4**4-2.0_ark*y2**2*y4*y5**3*sqrt(3.0_ark)-3.0_ark/4.0_ark&
        *y1**2*y5**4+5.0_ark/4.0_ark*y1**2*y4**4+y3**2*y4**2*y5**2)*fea334455 
        s3 = s2+(-6.0_ark*y4**2*y5**4+9.0_ark*y4**4*y5**2+y5**6)*fea555555+(y2*y3**3*y4**2+&
        y2*y3**3*y5**2+y1*y3**3*y4**2+y1*y2**3*y4**2+y1**3*y2*y4**2+&
        y1*y2**3*y5**2+y1**3*y3*y5**2+y1**3*y3*y4**2+y1**3*y2*y5**2+y2**3*y3*y4**2+&
        y1*y3**3*y5**2+y2**3*y3*y5**2)*fea233344+(y1*y2**3*y5**2*sqrt(3.0_ark)/6.0_ark-&
        y2**3*y3*y5**2*sqrt(3.0_ark)/3.0_ark-y2*y3**3*y5**2&
        *sqrt(3.0_ark)/3.0_ark+y1**3*y2*y4*y5-y1**3*y2*y5**2*sqrt(3.0_ark)/3.0_ark-&
        y1**3*y3*y4*y5-y1**3*y3*y5**2*sqrt(3.0_ark)/3.0_ark-y1*y3**3*y4**2*sqrt(3.0_ark&
        )/2.0_ark+y1*y3**3*y5**2*sqrt(3.0_ark)/6.0_ark-y2**3*y3*y4*y5+y2*y3**3*y4*&
        y5-y1*y2**3*y4**2*sqrt(3.0_ark)/2.0_ark)*fea233345+(-3.0_ark*y2**3*y4*y5**2&
        +y3**3*y4**3-3.0_ark*y3**3*y4*y5**2-3.0_ark*y1**3*y4*y5**2+y2**3*y4**3+&
        y1**3*y4**3)*fea111444+(y1*y2**3*y3**2+y1**3*y2**2*y3+y1**2*y2**3*y3+&
        y1*y2**2*y3**3+y1**2*y2*y3**3+y1**3*y2*y3**2)*fea111233 
        
        s4 = s3+(9.0_ark*y4**2*y5**4-6.0_ark*y4**4*y5**2+y4**6)*fea444444+(-5.0_ark&
        /3.0_ark*y1*y2**2*y4**2*y5*sqrt(3.0_ark)+y1*y2**2*y4**3-4.0_ark/3.0_ark*y1**2*&
        y3*y4**2*y5*sqrt(3.0_ark)-2.0_ark*y1**2*y2*y4**3-y1*y2**2*y5**3*sqrt(3.0_ark&
        )/3.0_ark+4.0_ark/3.0_ark*y2**2*y3*y4*y5**2-4.0_ark/3.0_ark*y2**2*y3*y5**3*sqrt(&
        3.0_ark)-2.0_ark*y1**2*y3*y4**3+7.0_ark/3.0_ark*y1*y2**2*y4*y5**2-2.0_ark/3.0_ark*y1&
        **2*y3*y4*y5**2+y1*y3**2*y4**3+4.0_ark/3.0_ark*y2*y3**2*y5**3*sqrt(3.0_ark)&
        +y1*y3**2*y5**3*sqrt(3.0_ark)/3.0_ark+4.0_ark/3.0_ark*y1**2*y2*y4**2*y5*sqrt(3.0_ark)&
        +4.0_ark/3.0_ark*y2*y3**2*y4*y5**2+5.0_ark/3.0_ark*y1*y3**2*y4**2*y5*sqrt(&
        3.0_ark)-2.0_ark/3.0_ark*y1**2*y2*y4*y5**2+7.0_ark/3.0_ark*y1*y3**2*y4*y5**2)*fea133444 
        
        s5 = s4+(-y1**3*y2*y4*y5+2.0_ark/3.0_ark*y2**3*y3*y5**2*sqrt(3.0_ark)+y1*y3&
        **3*y4**2*sqrt(3.0_ark)/2.0_ark+y1**3*y3*y4**2*sqrt(3.0_ark)/2.0_ark+y1**3*y3*&
        y5**2*sqrt(3.0_ark)/6.0_ark+y1**3*y2*y5**2*sqrt(3.0_ark)/6.0_ark+y1**3*y3*y4*y5+&
        y1*y2**3*y5**2*sqrt(3.0_ark)/6.0_ark+y1**3*y2*y4**2*sqrt(3.0_ark)/2.0_ark+&
        2.0_ark/3.0_ark*y2*y3**3*y5**2*sqrt(3.0_ark)-y1*y2**3*y4*y5+y1*y2**3*y4**2*sqrt(3.0_ark)/2.0_ark+&
        y1*y3**3*y5**2*sqrt(3.0_ark)/6.0_ark+y1*y3**3*y4*y5)*fea133345 

        v6 = s5+(-y2**2*y3*y4**2*y5+y1**2*y3*y4*y5**2*sqrt(3.0_ark)/3.0_ark+y2*y3**2*y4**2*y5+&
        y2*y3**2*y5**3-y1*y2**2*y5**3+4.0_ark/3.0_ark*y2**2*y3*y4*&
        y5**2*sqrt(3.0_ark)+4.0_ark/3.0_ark*y2*y3**2*y4*y5**2*sqrt(3.0_ark)-y1*y2**2*y4**2*y5+&
        4.0_ark/3.0_ark*y1*y3**2*y4*y5**2*sqrt(3.0_ark)-y2**2*y3*y5**3+y1*y3**2*y5**3+&
        y1**2*y2*y4*y5**2*sqrt(3.0_ark)/3.0_ark-y1**2*y2*y4**3*sqrt(3.0_ark)&
        +y1*y3**2*y4**2*y5-y1**2*y3*y4**3*sqrt(3.0_ark)+4.0_ark/3.0_ark*y1*y2**&
        2*y4*y5**2*sqrt(3.0_ark))*fea233445+(y2*y3**4*y4*sqrt(3.0_ark)-y1**4*y2*&
        y5+y2**4*y3*y4*sqrt(3.0_ark)-y1**4*y3*y4*sqrt(3.0_ark)+y2*y3**4*y5-2.0_ark*&
        y1*y2**4*y5+2.0_ark*y1*y3**4*y5-y1**4*y2*y4*sqrt(3.0_ark)+y1**4*y3*y5-y2&
        **4*y3*y5)*fea233335+(y2**2*y3**4+y1**4*y3**2+y1**2*y2**4+y2**4*y3&
        **2+y1**2*y3**4+y1**4*y2**2)*fea222233
        !
      endif
      !
      f =  v0+v1+v2+v3+v4+v5+v6


  end function poten_xy3_morbid_10

end module pot_user

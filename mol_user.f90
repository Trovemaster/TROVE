!
!  This unit is for a user defined potential 
!
module mol_user
  use accuracy
  use moltype
  use lapack
  use symmetry,     only : sym

  implicit none

  public MLpotential,MLdipole,MLpoten,ML_MEP,MLpoten_xy3_morbid_10

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  !interface MLpotential
  !  module procedure MLpotential_xy3_morbid_10
  !end interface 


  interface MLpoten
    module procedure MLpoten_xy3_morbid_10
  end interface 


  interface MLdipole
    module procedure MLdms2xyz_xy3_mb4
  end interface ! MLdipole

  !procedure (template) :: MLpoten_xy3_morbid_10


  contains




  !
  ! Defining potential energy function 
  !
  ! This type is for XY3-molecules, MORBID type of expansion #10
  ! see, e.g., JCP 117, 11265 (2002) for the defienition 
  !
  function MLpoten_xy3_morbid_10(local,xyz)  result(f)

   real(ark),intent(in) ::  local(:)
   real(ark),intent(in) ::  xyz(:,:)
   real(ark)            ::  f

   f = MLpotential_xy3_morbid_10(local,molec%force)

 end function MLpoten_xy3_morbid_10


  function ML_MEP(dim,rho)  result(f)

   integer(ik),intent(in) ::  dim
   real(ark),intent(in)   ::  rho
   real(ark)              ::  f(dim)

   f(1:3) = molec%req(1:3)
   f(4:5) = molec%alphaeq(1:2)
   f(6) = rho

 end function ML_MEP


  !
  ! Defining potential energy function 
  !
  ! This type is for XY3-molecules, MORBID type of expansion #10
  ! see, e.g., JCP 117, 11265 (2002) for the defienition 
  !
  function MLpotential(local,xyz,force)  result(f)

   real(ark),intent(in) ::  local(:)
   real(ark),intent(in) ::  xyz(:,:)
   real(ark),intent(in) ::  force(:)
   real(ark)            ::  f

   f = MLpotential_xy3_morbid_10(local,force)


 end function MLpotential




  !
  ! Defining potential energy function 
  !
  ! This type is for XY3-molecules, MORBID type of expansion #10
  ! see, e.g., JCP 117, 11265 (2002) for the defienition 
  !
  function MLpotential_xy3_morbid_10(local,force)  result(f)

   real(ark),intent(in) ::  local(:)
   real(ark),intent(in) ::  force(:)
   real(ark)            ::  f
   !
   integer(ik)          :: N

   real(ark)            ::  r14,r24,r34,alpha1,alpha2,alpha3,alpha0
   real(ark)            ::  y1,y2,y3,y4,y5
   real(ark)            ::  v0,rhoe,v1,v2,v3,v4,v5,v6,rho,re14,aa1,tau
   real(ark)            ::  sinrho,alpha,coro,tau_2,cosalpha,sinphi,phi1,phi2,phi3,delta
      !
      if (verbose>=6) write(out,"('MLpotential_xy3_morbid_10/start')") 
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
            write (out,"('MLpotential_xy3_morbid_10: coord. type ',a,' unknown')") trim(molec%coords_transform)
            stop 'MLpotential_xy3_morbid_10 - bad coord. type'
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
         case('R-S-DELTA','R-2D-DELTA','R-S-RHO','SYM-DELTA','R-SYMPHI-DELTA','R-PHI-DELTA','R-S-DELTA-MEP','R-PHI-DELTA-MEP',&
              'R-SYMPHI-DELTA-MEP','R-A2A3-DELTA')
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
            phi1 = 2.0_ark-phi2-phi3
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
         case('R-A2-A3-TAU')
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
      !if (trim(molec%potentype)=='POTEN_XY3_MORBID_10_MEP') re14 = ML_mep_oh3p(rho)
      !
      y1=1.0_ark-exp(-aa1*(r14-re14))
      y2=1.0_ark-exp(-aa1*(r24-re14))
      y3=1.0_ark-exp(-aa1*(r34-re14))
      !
      y4=(2.0_ark*alpha1-alpha2-alpha3)/sqrt(6.0_ark)
      y5=(alpha2-alpha3)/sqrt(2.0_ark)
      !
      f = poten_xy3_morbid_10(y1,y2,y3,y4,y5,coro,force)
      !
      if (verbose>=6) write(out,"('MLpotential_xy3_morbid_10/end')") 
 
 end function MLpotential_xy3_morbid_10


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




 function MLdms2xyz_xy3_mb4(xyz) result(f)

! return electric dipole moment vector component for xy3 nonplanar molecule
! in the Eckart frame
!
! see:
! S. N. Yurchenko, W. Thiel, M. Carvajal, H. Lin, and P. Jensen, 
! Adv. Quant. Chem. 48 (2005) 209

    !
    real(ark), intent(in) :: xyz(molec%natoms, 3)
    real(ark)             :: f(3)
    !
    integer(ik) parmax, ierror, i, rank0
    real(ark) charge, param(extF%nterms(1) - 2), x(molec%natoms,3), y(molec%natoms,3), &
              r1, r2, r3, a1, a2, a3, tau, local(7), dip(4), d_t(4), ec, vcprm, beta
    real(rk) edip(4, 1), tmat(4, 3), tsing(3), wspace(32 * 3),smat(3,3)
    !
    real(ark) f_t, bmat(3), levich(3, 3, 3), rcond, cm(3),v12(3),v23(3),v31(3),n1(3),n2(3),n3(3)
    integer(ik) ix, iy, iz, parmaxB,lspace
    !

    !
    ! elementary charge in coulombs
    !
    ec    = 1.60217733e-19_ark
    vcprm = 8.854187817e-12_ark
    !
    ! elementary charge in electrostatic units (e.s.u.)
    !
    ec = ec / sqrt(4.0e-09_ark * pi * vcprm)
    !
    ! elementary charge in debye/angstrom = 1.0e8 * debye/cm = 1.0e-10 * e.s.u.
    !
    ec = 1.0e+10_ark * ec
    !
    charge          = extF%coef(1, 1)
    !
    do i = 1, 3
       y(i,:) = xyz(i + 1,:)
    end do
    y(4,:) = xyz(1,:)
    !
    forall(i=1:3) cm(i) = ( sum(xyz(1:4,i)*molec%atommasses(1:4)) )/sum(molec%atommasses(1:4))
    !
    forall(i=1:4) x(i,1:3) = y(i,1:3) - y(4,1:3)
    !
    r1 = sqrt(sum(x(1,:)**2))
    r2 = sqrt(sum(x(2,:)**2))
    r3 = sqrt(sum(x(3,:)**2))
    !
    a1 = acos(sum(x(2,:)*x(3,:))/(r2*r3))
    a2 = acos(sum(x(1,:)*x(3,:))/(r1*r3))
    a3 = acos(sum(x(2,:)*x(1,:))/(r2*r1))
    !
    local(1) = r1
    local(2) = r2
    local(3) = r3
    !
    local(4) = a1
    local(5) = a2
    local(6) = a3
    !
    parmax          = extF%coef(2, 1) 
    param(1:parmax) = extF%coef(3:parmax + 2, 1)
    !
    dip(1) = MLdms2loc_xy3(1, parmax, param, local(1:6))
    dip(2) = MLdms2loc_xy3(2, parmax, param, local(1:6))
    dip(3) = MLdms2loc_xy3(3, parmax, param, local(1:6))
    !
    dip(1:3) = dip(1:3) + ec * cm(1:3) * charge
    !
    edip(1:3, 1) = real(dip(1:3), kind = rk)
    !
    tmat(1, :) = real(x(1,:) / r1, kind = rk)
    tmat(2, :) = real(x(2,:) / r2, kind = rk)
    tmat(3, :) = real(x(3,:) / r3, kind = rk)
    !
    lspace = size(wspace)
    !
    call dgelss(3,3,1,tmat(1:3,:),3,edip(1:3,1),3,tsing,-1.0d-12,rank0,wspace,lspace,ierror)
    !
    !if (abs(a1+a2+a3-2.0_ark*pi)<=0.1) then  
      !
      tmat(1, :) = real(x(1,:) / r1, kind = rk)
      tmat(2, :) = real(x(2,:) / r2, kind = rk)
      tmat(3, :) = real(x(3,:) / r3, kind = rk)
      !
      edip(1:3, 1) = real(dip(1:3), kind = rk)
      !
      v12 = MLvector_product(x(1,:)/r1,x(2,:)/r2)
      v23 = MLvector_product(x(2,:)/r2,x(3,:)/r3)
      v31 = MLvector_product(x(3,:)/r3,x(1,:)/r1)
      !
      n3 = v12 + v23 + v31
      !
      n3 = n3/sqrt(sum(n3(:)**2))
      !
      beta = acos(sum(x(1,:)*n3(:))/r1)
      !
      local(7) = beta
      !
      parmaxB          = extF%coef(parmax+3, 1) 
      param(1:parmaxB) = extF%coef(parmax + 4:extF%nterms(1),1)
      !
      dip(4) = MLdms2loc_A_xy3(parmaxB, param(1:parmaxB), local(1:7))
      !
      edip(4, 1) = real(dip(4), kind = rk)
      !
      tmat(4, :) = real(n3(:), kind = rk)
      !
      call dgelss(4, 3, 1, tmat, 4, edip, 4, tsing,-1.0d-12, rank0, wspace, lspace , ierror)
      !
    !endif 
    !
    f(1:3) = real(edip(1:3, 1), kind = ark)
    !
    f(1:3) = f(1:3) + ec * cm(1:3) * charge
    !
 end function MLdms2xyz_xy3_mb4
 !


 function MLdms2loc_A_xy3(parmax, param, local) result(f)

! return electric dipole moment vector component for xy3 molecule in the space-fixed frame
!

    integer, intent(in)   ::  parmax
    real(ark), intent(in) ::  param(parmax)
    real(ark), intent(in) ::  local(7)
    real(ark)             ::  f

    real(ark) r14, r24, r34, alpha1, alpha2, alpha3, rhobar
    real(ark) rhoe
    real(ark) re14, beta,de,y1,y2,y3,y4,y5,alpha,rho,sinrho,cosrho,drho
    real(ark) v0,v1,v2,v3,v4,v5,v6




    real(ark)                                        &
      fea    ,fea1  ,                                    &
      fea11  ,fea12  ,fea14  ,fea44  ,                   &
      fea111 ,fea112 ,fea114 ,fea123 ,                   &
      fea124 ,fea144 ,fea155 ,fea455 ,                   &
      fea1111,fea1112,fea1114,fea1122,                   &
      fea1123,fea1124,fea1125,fea1144,                   &
      fea1155,fea1244,fea1255,fea1444,                   &
      fea1455,fea4444
    real(ark) ::   &
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

    real(ark) :: s1,s2,s3,s4,s5,s6
                    
    real(ark) :: &      
      f1a,f2a,f3a,f4a,f5a,f6a,f7a,f0a1, &
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
      f0a44444 ,f1a44444 ,           &
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
      f2a11134 ,f0a11123 ,f1a11123 ,f2a11123 ,f0a555555,f1a555555,&
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



   rhoe   = param(  1)* pi / 180.0_ark 

   re14       = param( 2)
   beta       = param( 3)**2
   de         = param( 4)

   f1a        = param(  5)
   f2a        = param(  6)
   f3a        = param(  7)
   f4a        = param(  8)
   f5a        = param(  9)
   f6a        = param( 10)
   f7a        = param( 11)
   f0a1       = param( 12)
   f1a1       = param( 13)
   f2a1       = param( 14)
   f3a1       = param( 15)
   f4a1       = param( 16)
   f5a1       = param( 17)
   f6a1       = param( 18)
   f0a11      = param( 19)
   f1a11      = param( 20)
   f2a11      = param( 21)
   f3a11      = param( 22)
   f4a11      = param( 23)
   f0a12      = param( 24)
   f1a12      = param( 25)
   f2a12      = param( 26)
   f3a12      = param( 27)
   f4a12      = param( 28)
   f0a14      = param( 29)
   f1a14      = param( 30)
   f2a14      = param( 31)
   f3a14      = param( 32)
   f4a14      = param( 33)
   f0a44      = param( 34)
   f1a44      = param( 35)
   f2a44      = param( 36)
   f3a44      = param( 37)
   f4a44      = param( 38)
   f0a111     = param( 39)
   f1a111     = param( 40)
   f2a111     = param( 41)
   f3a111     = param( 42)
   f0a112     = param( 43)
   f1a112     = param( 44)
   f2a112     = param( 45)
   f3a112     = param( 46)
   f0a114     = param( 47)
   f1a114     = param( 48)
   f2a114     = param( 49)
   f3a114     = param( 50)
   f0a123     = param( 51)
   f1a123     = param( 52)
   f2a123     = param( 53)
   f3a123     = param( 54)
   f0a124     = param( 55)
   f1a124     = param( 56)
   f2a124     = param( 57)
   f3a124     = param( 58)
   f0a144     = param( 59)
   f1a144     = param( 60)
   f2a144     = param( 61)
   f3a144     = param( 62)
   f0a155     = param( 63)
   f1a155     = param( 64)
   f2a155     = param( 65)
   f3a155     = param( 66)
   f0a455     = param( 67)
   f1a455     = param( 68)
   f2a455     = param( 69)
   f3a455     = param( 70)
   f0a1111    = param( 71)
   f1a1111    = param( 72)
   f2a1111    = param( 73)
   f0a1112    = param( 74)
   f1a1112    = param( 75)
   f2a1112    = param( 76)
   f0a1114    = param( 77)
   f1a1114    = param( 78)
   f2a1114    = param( 79)
   f0a1122    = param( 80)
   f1a1122    = param( 81)
   f2a1122    = param( 82)
   f0a1123    = param( 83)
   f1a1123    = param( 84)
   f2a1123    = param( 85)
   f0a1124    = param( 86)
   f1a1124    = param( 87)
   f2a1124    = param( 88)
   f0a1125    = param( 89)
   f1a1125    = param( 90)
   f2a1125    = param( 91)
   f0a1144    = param( 92)
   f1a1144    = param( 93)
   f2a1144    = param( 94)
   f0a1155    = param( 95)
   f1a1155    = param( 96)
   f2a1155    = param( 97)
   f0a1244    = param( 98)
   f1a1244    = param( 99)
   f2a1244    = param(100)
   f0a1255    = param(101)
   f1a1255    = param(102)
   f2a1255    = param(103)
   f0a1444    = param(104)
   f1a1444    = param(105)
   f2a1444    = param(106)
   f0a1455    = param(107)
   f1a1455    = param(108)
   f2a1455    = param(109)
   f0a4444    = param(110)
   f1a4444    = param(111)
   f2a4444    = param(112)
   !
   !
   if (parmax>112) then 
     !
     f0a44444   = param(113)
     f1a44444   = param(114)
     f2a44444   = param(115)
     f0a33455   = param(116)
     f1a33455   = param(117)
     f2a33455   = param(118)
     f0a33445   = param(119)
     f1a33445   = param(120)
     f2a33445   = param(121)
     f0a33345   = param(122)
     f1a33345   = param(123)
     f2a33345   = param(124)
     f0a33344   = param(125)
     f1a33344   = param(126)
     f2a33344   = param(127)
     f0a33334   = param(128)
     f1a33334   = param(129)
     f2a33334   = param(130)
     f0a33333   = param(131)
     f1a33333   = param(132)
     f2a33333   = param(133)
     f0a25555   = param(134)
     f1a25555   = param(135)
     f2a25555   = param(136)
     f0a24455   = param(137)
     f1a24455   = param(138)
     f2a24455   = param(139)
     f0a24445   = param(140)
     f1a24445   = param(141)
     f2a24445   = param(142)
     f0a23333   = param(143)
     f1a23333   = param(144)
     f2a23333   = param(145)
     f0a13455   = param(146)
     f1a13455   = param(147)
     f2a13455   = param(148)
     f0a13445   = param(149)
     f1a13445   = param(150)
     f2a13445   = param(151)
     f0a13345   = param(152)
     f1a13345   = param(153)
     f2a13345   = param(154)
     f0a12355   = param(155)
     f1a12355   = param(156)
     f2a12355   = param(157)
     f0a11334   = param(158)
     f1a11334   = param(159)
     f2a11334   = param(160)
     f0a11333   = param(161)
     f1a11333   = param(162)
     f2a11333   = param(163)
     f0a11255   = param(164)
     f1a11255   = param(165)
     f2a11255   = param(166)
     f0a11245   = param(167)
     f1a11245   = param(168)
     f2a11245   = param(169)
     f0a11234   = param(170)
     f1a11234   = param(171)
     f2a11234   = param(172)
     f0a11233   = param(173)
     f1a11233   = param(174)
     f2a11233   = param(175)
     f0a11135   = param(176)
     f1a11135   = param(177)
     f2a11135   = param(178)
     f0a11134   = param(179)
     f1a11134   = param(180)
     f2a11134   = param(181)
     f0a11123   = param(182)
     f1a11123   = param(183)
     f2a11123   = param(184)
     f0a555555  = param(185)
     f1a555555  = param(186)
     f2a555555  = param(187)
     f0a444444  = param(188)
     f1a444444  = param(189)
     f2a444444  = param(190)
     f0a335555  = param(191)
     f1a335555  = param(192)
     f2a335555  = param(193)
     f0a334455  = param(194)
     f1a334455  = param(195)
     f2a334455  = param(196)
     f0a334445  = param(197)
     f1a334445  = param(198)
     f2a334445  = param(199)
     f0a333555  = param(200)
     f1a333555  = param(201)
     f2a333555  = param(202)
     f0a333333  = param(203)
     f1a333333  = param(204)
     f2a333333  = param(205)
     f0a244555  = param(206)
     f1a244555  = param(207)
     f2a244555  = param(208)
     f0a244455  = param(209)
     f1a244455  = param(210)
     f2a244455  = param(211)
     f0a233445  = param(212)
     f1a233445  = param(213)
     f2a233445  = param(214)
     f0a233444  = param(215)
     f1a233444  = param(216)
     f2a233444  = param(217)
     f0a233345  = param(218)
     f1a233345  = param(219)
     f2a233345  = param(220)
     f0a233344  = param(221)
     f1a233344  = param(222)
     f2a233344  = param(223)
     f0a233335  = param(224)
     f1a233335  = param(225)
     f2a233335  = param(226)
     f0a223355  = param(227)
     f1a223355  = param(228)
     f2a223355  = param(229)
     f0a222335  = param(230)
     f1a222335  = param(231)
     f2a222335  = param(232)
     f0a222334  = param(233)
     f1a222334  = param(234)
     f2a222334  = param(235)
     f0a222333  = param(236)
     f1a222333  = param(237)
     f2a222333  = param(238)
     f0a222255  = param(239)
     f1a222255  = param(240)
     f2a222255  = param(241)
     f0a222245  = param(242)
     f1a222245  = param(243)
     f2a222245  = param(244)
     f0a222233  = param(245)
     f1a222233  = param(246)
     f2a222233  = param(247)
     f0a222224  = param(248)
     f1a222224  = param(249)
     f2a222224  = param(250)
     f0a145555  = param(251)
     f1a145555  = param(252)
     f2a145555  = param(253)
     f0a134444  = param(254)
     f1a134444  = param(255)
     f2a134444  = param(256)
     f0a133444  = param(257)
     f1a133444  = param(258)
     f2a133444  = param(259)
     f0a133345  = param(260)
     f1a133345  = param(261)
     f2a133345  = param(262)
     f0a133334  = param(263)
     f1a133334  = param(264)
     f2a133334  = param(265)
     f0a133333  = param(266)
     f1a133333  = param(267)
     f2a133333  = param(268)
     f0a124555  = param(269)
     f1a124555  = param(270)
     f2a124555  = param(271)
     f0a124455  = param(272)
     f1a124455  = param(273)
     f2a124455  = param(274)
     f0a123455  = param(275)
     f1a123455  = param(276)
     f2a123455  = param(277)
     f0a123345  = param(278)
     f1a123345  = param(279)
     f2a123345  = param(280)
     f0a113555  = param(281)
     f1a113555  = param(282)
     f2a113555  = param(283)
     f0a113345  = param(284)
     f1a113345  = param(285)
     f2a113345  = param(286)
     f0a112355  = param(287)
     f1a112355  = param(288)
     f2a112355  = param(289)
     f0a112335  = param(290)
     f1a112335  = param(291)
     f2a112335  = param(292)
     f0a112233  = param(293)
     f1a112233  = param(294)
     f2a112233  = param(295)
     f0a111444  = param(296)
     f1a111444  = param(297)
     f2a111444  = param(298)
     f0a111234  = param(299)
     f1a111234  = param(300)
     f2a111234  = param(301)
     f0a111233  = param(302)
     f1a111233  = param(303)
     f2a111233  = param(304)
     f0a111123  = param(305)
     f1a111123  = param(306)
     f2a111123  = param(307)
     !
   endif 

   r14    = local(1) ;  r24    = local(2) ;  r34    = local(3)
   alpha1 = local(4) ;  alpha2 = local(5) ;  alpha3 = local(6)


  y4=(2.0_ark*alpha1-alpha2-alpha3)/sqrt(6.0_ark)
  y5=(alpha2-alpha3)/sqrt(2.0_ark)

  alpha=(alpha1+alpha2+alpha3)/3.0_ark
  !
  !rho=pi-asin(2.0_ark*sin(alpha*0.50_ark)/sqrt(3.0_ark))

  if ( 2.0_ark*sin(alpha*0.50_ark)/sqrt(3.0_ark).ge.1.00_ark ) then 
    sinrho=1.0_ark 
  else 
    sinrho = 2.0_ark*sin(alpha*0.50_ark)/sqrt(3.0_ark)
  endif 
  !
  cosrho = sqrt(1.0_ark-sinrho**2)*sign(1.0_ark,cos(local(7)))
  drho=(sin(rhoe)-sinrho)
  !
  rhobar = local(7)
  !
  !cosrho = cos(rhobar)
  !drho=(sin(rhoe)-sin(rhobar))
  !
  !drho=-sqrt(1.0_ark-sinrho**2)
  !
  y1=1.00_ark*(r14-re14) *exp(-beta*(r14-re14)**2)
  y2=1.00_ark*(r24-re14) *exp(-beta*(r24-re14)**2)
  y3=1.00_ark*(r34-re14) *exp(-beta*(r34-re14)**2)
  !
  v0=de+f1a*drho+f2a*drho**2+f3a*drho**3+f4a*drho**4+f5a*drho**5 &
        +f6a*drho**6+f7a*drho**7

  fea1= f0a1+f1a1*drho+f2a1*drho**2+f3a1*drho**3+f4a1*drho**4+f5a1*drho**5+f6a1*drho**6

  fea11=   f0a11+f1a11*drho+f2a11*drho**2+f3a11*drho**3+f4a11*drho**4
  fea12=   f0a12+f1a12*drho+f2a12*drho**2+f3a12*drho**3+f4a12*drho**4
  fea14=   f0a14+f1a14*drho+f2a14*drho**2+f3a14*drho**3+f4a14*drho**4
  fea44=   f0a44+f1a44*drho+f2a44*drho**2+f3a44*drho**3+f4a44*drho**4
 
  fea111= f0a111+f1a111*drho+f2a111*drho**2+f3a111*drho**3
  fea112= f0a112+f1a112*drho+f2a112*drho**2+f3a112*drho**3
  fea114= f0a114+f1a114*drho+f2a114*drho**2+f3a114*drho**3
  fea123= f0a123+f1a123*drho+f2a123*drho**2+f3a123*drho**3
  fea124= f0a124+f1a124*drho+f2a124*drho**2+f3a124*drho**3
  fea144= f0a144+f1a144*drho+f2a144*drho**2+f3a144*drho**3
  fea155= f0a155+f1a155*drho+f2a155*drho**2+f3a155*drho**3
  fea455= f0a455+f1a455*drho+f2a455*drho**2+f3a455*drho**3

  fea1111= f0a1111+f1a1111*drho+f2a1111*drho**2
  fea1112= f0a1112+f1a1112*drho+f2a1112*drho**2
  fea1114= f0a1114+f1a1114*drho+f2a1114*drho**2
  fea1122= f0a1122+f1a1122*drho+f2a1122*drho**2
  fea1123= f0a1123+f1a1123*drho+f2a1123*drho**2
  fea1124= f0a1124+f1a1124*drho+f2a1124*drho**2
  fea1125= f0a1125+f1a1125*drho+f2a1125*drho**2
  fea1144= f0a1144+f1a1144*drho+f2a1144*drho**2
  fea1155= f0a1155+f1a1155*drho+f2a1155*drho**2
  fea1244= f0a1244+f1a1244*drho+f2a1244*drho**2
  fea1255= f0a1255+f1a1255*drho+f2a1255*drho**2
  fea1444= f0a1444+f1a1444*drho+f2a1444*drho**2
  fea1455= f0a1455+f1a1455*drho+f2a1455*drho**2
  fea4444= f0a4444+f1a4444*drho+f2a4444*drho**2


  if (parmax>112) then 
     !
     fea44444 = f0a44444  + f1a44444 *drho+ f2a44444 *drho**2
     fea33455 = f0a33455  + f1a33455 *drho+ f2a33455 *drho**2
     fea33445 = f0a33445  + f1a33445 *drho+ f2a33445 *drho**2
     fea33345 = f0a33345  + f1a33345 *drho+ f2a33345 *drho**2
     fea33344 = f0a33344  + f1a33344 *drho+ f2a33344 *drho**2
     fea33334 = f0a33334  + f1a33334 *drho+ f2a33334 *drho**2
     fea33333 = f0a33333  + f1a33333 *drho+ f2a33333 *drho**2
     fea25555 = f0a25555  + f1a25555 *drho+ f2a25555 *drho**2
     fea24455 = f0a24455  + f1a24455 *drho+ f2a24455 *drho**2
     fea24445 = f0a24445  + f1a24445 *drho+ f2a24445 *drho**2
     fea23333 = f0a23333  + f1a23333 *drho+ f2a23333 *drho**2
     fea13455 = f0a13455  + f1a13455 *drho+ f2a13455 *drho**2
     fea13445 = f0a13445  + f1a13445 *drho+ f2a13445 *drho**2
     fea13345 = f0a13345  + f1a13345 *drho+ f2a13345 *drho**2
     fea12355 = f0a12355  + f1a12355 *drho+ f2a12355 *drho**2
     fea11334 = f0a11334  + f1a11334 *drho+ f2a11334 *drho**2
     fea11333 = f0a11333  + f1a11333 *drho+ f2a11333 *drho**2
     fea11255 = f0a11255  + f1a11255 *drho+ f2a11255 *drho**2
     fea11245 = f0a11245  + f1a11245 *drho+ f2a11245 *drho**2
     fea11234 = f0a11234  + f1a11234 *drho+ f2a11234 *drho**2
     fea11233 = f0a11233  + f1a11233 *drho+ f2a11233 *drho**2
     fea11135 = f0a11135  + f1a11135 *drho+ f2a11135 *drho**2
     fea11134 = f0a11134  + f1a11134 *drho+ f2a11134 *drho**2
     fea11123 = f0a11123  + f1a11123 *drho+ f2a11123 *drho**2
     fea555555= f0a555555 + f1a555555*drho+ f2a555555*drho**2
     fea444444= f0a444444 + f1a444444*drho+ f2a444444*drho**2
     fea335555= f0a335555 + f1a335555*drho+ f2a335555*drho**2
     fea334455= f0a334455 + f1a334455*drho+ f2a334455*drho**2
     fea334445= f0a334445 + f1a334445*drho+ f2a334445*drho**2
     fea333555= f0a333555 + f1a333555*drho+ f2a333555*drho**2
     fea333333= f0a333333 + f1a333333*drho+ f2a333333*drho**2
     fea244555= f0a244555 + f1a244555*drho+ f2a244555*drho**2
     fea244455= f0a244455 + f1a244455*drho+ f2a244455*drho**2
     fea233445= f0a233445 + f1a233445*drho+ f2a233445*drho**2
     fea233444= f0a233444 + f1a233444*drho+ f2a233444*drho**2
     fea233345= f0a233345 + f1a233345*drho+ f2a233345*drho**2
     fea233344= f0a233344 + f1a233344*drho+ f2a233344*drho**2
     fea233335= f0a233335 + f1a233335*drho+ f2a233335*drho**2
     fea223355= f0a223355 + f1a223355*drho+ f2a223355*drho**2
     fea222335= f0a222335 + f1a222335*drho+ f2a222335*drho**2
     fea222334= f0a222334 + f1a222334*drho+ f2a222334*drho**2
     fea222333= f0a222333 + f1a222333*drho+ f2a222333*drho**2
     fea222255= f0a222255 + f1a222255*drho+ f2a222255*drho**2
     fea222245= f0a222245 + f1a222245*drho+ f2a222245*drho**2
     fea222233= f0a222233 + f1a222233*drho+ f2a222233*drho**2
     fea222224= f0a222224 + f1a222224*drho+ f2a222224*drho**2
     fea145555= f0a145555 + f1a145555*drho+ f2a145555*drho**2
     fea134444= f0a134444 + f1a134444*drho+ f2a134444*drho**2
     fea133444= f0a133444 + f1a133444*drho+ f2a133444*drho**2
     fea133345= f0a133345 + f1a133345*drho+ f2a133345*drho**2
     fea133334= f0a133334 + f1a133334*drho+ f2a133334*drho**2
     fea133333= f0a133333 + f1a133333*drho+ f2a133333*drho**2
     fea124555= f0a124555 + f1a124555*drho+ f2a124555*drho**2
     fea124455= f0a124455 + f1a124455*drho+ f2a124455*drho**2
     fea123455= f0a123455 + f1a123455*drho+ f2a123455*drho**2
     fea123345= f0a123345 + f1a123345*drho+ f2a123345*drho**2
     fea113555= f0a113555 + f1a113555*drho+ f2a113555*drho**2
     fea113345= f0a113345 + f1a113345*drho+ f2a113345*drho**2
     fea112355= f0a112355 + f1a112355*drho+ f2a112355*drho**2
     fea112335= f0a112335 + f1a112335*drho+ f2a112335*drho**2
     fea112233= f0a112233 + f1a112233*drho+ f2a112233*drho**2
     fea111444= f0a111444 + f1a111444*drho+ f2a111444*drho**2
     fea111234= f0a111234 + f1a111234*drho+ f2a111234*drho**2
     fea111233= f0a111233 + f1a111233*drho+ f2a111233*drho**2
     fea111123= f0a111123 + f1a111123*drho+ f2a111123*drho**2
    !
  endif 
    


      v1 = (y3+y2+y1)*fea1

      v2 = (y2*y3+y1*y3+y1*y2)*fea12                                                                 &
       +(y2**2+y3**2+y1**2)*fea11                                                                    &
       +(-sqrt(3.0_ark)*y3*y5/2.0_ark-y3*y4/2.0_ark+y1*y4+sqrt(3.0_ark)*y2*y5/2.0_ark-y2*y4/2.0_ark)*fea14 &
       +(y5**2+y4**2)*fea44

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

      v5 = 0 ;  v6 = 0 
      !
      if (parmax>112) then 
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
       f = ( v0+(v1+v2+v3+v4+v5+v6) )*cosrho
       !
 end function MLdms2loc_A_xy3


 !
 !--------------------------------------------------------------------------
 !orig from DMS_light.f90
 !
 ! return electric dipole moment vector component for xy3 nonplanar molecule
 ! in the molecular-bond frame
 !
 ! see:
 ! S. N. Yurchenko, W. Thiel, M. Carvajal, H. Lin, and P. Jensen, 
 ! Adv. Quant. Chem. 48 (2005) 209
 !
 !--------------------------------------------------------------------------
 !
 function MLdms2loc_xy3(ix, parmax, param, local) result(f)

! return electric dipole moment vector component for xy3 molecule in the space-fixed frame
!
! see:
! S. N. Yurchenko, W. Thiel, M. Carvajal, H. Lin, and P. Jensen, 
! Adv. Quant. Chem. 48 (2005) 209

    integer, intent(in)   ::  ix, parmax
    real(ark), intent(in) ::  param(parmax)
    real(ark), intent(in) ::  local(6)
    real(ark)             ::  f

    real(ark) r14, r24, r34, alpha1, alpha2, alpha3
    real(ark) xi1, xi2, xi3, xi4, xi5, xi6
    real(ark) alphae
    real(ark) alphaedg
    real(ark) re14, beta, t1, t2, t3, t4, t5, t6
    real(ark) F5,    F4,    F3,    F1,    F56,   F55,   F46,   F44,   &
              F36,   F35,   F34,   F33,   F23,   F16,   F14,   F13,   &
              F11,   F556,  F555,  F466,  F456,  F445,  F444,  F344,  &
              F334,  F333,  F266,  F256,  F255,  F246,  F245,  F235,  &
              F234,  F226,  F225,  F223,  F156,  F155,  F146,  F144,  &
              F136,  F135,  F133,  F124,  F123,  F115,  F114,  F112,  &
              F111,  F6666, F5666, F5566, F4556, F4555, F4466, F4456, &
              F4445, F4444, F3666, F3566, F3556, F3555, F3445, F3335, &
              F2466, F2456, F2455, F2445, F2444, F2366, F2356, F2345, &
              F2344, F2336, F2335, F2334, F2333, F2266, F2256, F2255, &
              F2246, F2245, F2244, F2233, F2225, F2224, F2222, F1666, &
              F1566, F1466, F1456, F1445, F1444, F1366, F1355, F1335, &
              F1334, F1256, F1246, F1245, F1244, F1235, F1234, F1233, &
              F1225, F1222, F1156, F1155, F1146, F1144, F1136, F1126, &
              F1124, F1123, F1122, F1115, F1114, F1112, F1111
    real(ark) F56666, &
              F55666,F55555,F46666,F45666,F45566,F44666,F44566,                &
              F44466,F44456,F44445,F44444,F36666,F35666,F35566,                &
              F35556,F35555,F34666,F34566,F34556,F34555,F34466,                &
              F34456,F34455,F34446,F34445,F34444,F33666,F33566,                &
              F33556,F33555,F33466,F33455,F33446,F33444,F33366,                &
              F33356,F33355,F33346,F33345,F33344,F33336,F33335,                &
              F33333,F23556,F23555,F23466,F23456,F23446,F23444,                &
              F23366,F23356,F23355,F23346,F23344,F23336,F23335,                &
              F23334,F23333,F22456,F22446,F22346,F22336,F22334,                &
              F22333,F22224,F15566,F15556,F15555,F14666,F14556,                &
              F14466,F14456,F14445,F14444,F13666,F13566,F13556,                &
              F13456,F13446,F13445,F13356,F13346,F13345,F13344,                &
              F13335,F13334,F13333,F12666,F12466,F12455,F12444,                &
              F12356,F12355,F12345,F12344,F12336,F12334,F12333,                &
              F12266,F12255,F12236,F12233,F12225,F11556,F11555,                &
              F11456,F11455,F11446,F11444,F11356,F11345,F11336,                &
              F11335,F11333,F11266,F11255,F11245,F11244,F11236,                &
              F11234,F11224,F11223,F11166,F11156,F11145,F11144,                &
              F11136,F11135,F11133,F11124,F11123,F11116,F11114,                &
              F11112,F11111
    real(ark) F666666,F555666,F555566,F555556,F466666,F456666,                  &
              F455666,F445666,F445566,F445555,F444666,F444566,F444466,F444456,  &
              F444446,F444444,F266666,F256666,F255666,F255566,F255556,F255555,  &
              F246666,F245666,F245566,F245556,F245555,F244666,F244566,F244556,  &
              F244555,F244466,F244456,F244455,F244446,F244445,F244444,F235566,  &
              F235556,F235555,F234666,F234566,F234466,F234456,F234446,F234444,  &
              F233666,F233466,F233455,F233336,F226666,F225666,F225566,F225556,  &
              F225555,F224666,F224566,F224556,F224555,F224466,F224456,F224455,  &
              F224446,F224445,F224444,F223666,F223566,F223556,F223456,F223446,  &
              F223445,F223444,F223366,F223356,F223346,F223344,F223333,F222666,  &
              F222566,F222556,F222555,F222466,F222456,F222455,F222446,F222445,  &
              F222444,F222366,F222356,F222355,F222346,F222345,F222344,F222336,  &
              F222335,F222334,F222333,F222266,F222256,F222255,F222246,F222245,  &
              F222244,F222236,F222234,F222226,F222225,F222224,F222223,F222222,  &
              F166666,F155666,F155556,F146666,F145666,F145566,F144666,F144566,  &
              F144466,F144456,F144446,F144444,F134666,F134566,F134556,F134555,  &
              F134456,F134455,F133566,F133556,F126666,F125666,F125566,F125556,  &
              F125555,F124455,F124446,F124445,F124444,F123566,F123555,F123466,  &
              F123456,F123446,F123444,F123366,F123356,F123355,F122666,F122555,  &
              F122466,F122456,F122455,F122446,F122445,F122444,F122346,F122345,  &
              F122344,F122336,F122334,F122333,F122266,F122256,F122255,F122246,  &
              F122245,F122244,F122236,F122235,F122234,F122226,F122225,F122224,  &
              F122223,F122222,F115666,F115566,F115555,F114666,F114556,F114466,  &
              F114456,F114446,F114444,F113336,F113335,F113334,F113333,F112666,  &
              F112566,F112556,F112555,F112466,F112456,F112455,F112446,F112445,  &
              F112444,F112366,F112356,F112345,F112344,F112333,F112266,F112256,  &
              F112255,F112246,F112245,F112244,F112236,F112235,F112234,F112233,  &
              F111666,F111556,F111466,F111456,F111446,F111444,F111356,F111346,  &
              F111344,F111336,F111335,F111334,F111266,F111255,F111246,F111236,  &
              F111234,F111223,F111222,F111166,F111156,F111145,F111144,F111136,  &
              F111133,F111126,F111124,F111123,F111116,F111114,F111112,F111111 
    real(ark) s1, s2, s3, s4, s5, gamma, delta, mu0
    real(ark) ta1, ta2, ta3


    alphaedg = param( 1)
    alphae   = alphaedg * pi / 180.0_ark

    re14    =  param(  2)
    beta    =  param(  3)**2
    gamma   =  param(  4)**2
    delta   =  param(  5)
    mu0     =  param(  6)
    F1      =  param(    7)
    F3      =  param(    8)
    F4      =  param(    9)
    F5      =  param(   10)

   if (parmax>=11) then
      F11     =  param(   11)
      F13     =  param(   12)
      F14     =  param(   13)
      F16     =  param(   14)
      F23     =  param(   15)
      F33     =  param(   16)
      F34     =  param(   17)
      F35     =  param(   18)
      F36     =  param(   19)
      F44     =  param(   20)
      F46     =  param(   21)
      F55     =  param(   22)
      F56     =  param(   23)
   endif
   if (parmax>=24) then
      F111    =  param(   24)
      F112    =  param(   25)
      F114    =  param(   26)
      F115    =  param(   27)
      F123    =  param(   28)
      F124    =  param(   29)
      F133    =  param(   30)
      F135    =  param(   31)
      F136    =  param(   32)
      F144    =  param(   33)
      F146    =  param(   34)
      F155    =  param(   35)
      F156    =  param(   36)
      F223    =  param(   37)
      F225    =  param(   38)
      F226    =  param(   39)
      F234    =  param(   40)
      F235    =  param(   41)
      F245    =  param(   42)
      F246    =  param(   43)
      F255    =  param(   44)
      F256    =  param(   45)
      F266    =  param(   46)
      F333    =  param(   47)
      F334    =  param(   48)
      F344    =  param(   49)
      F444    =  param(   50)
      F445    =  param(   51)
      F456    =  param(   52)
      F466    =  param(   53)
      F555    =  param(   54)
      F556    =  param(   55)
   endif
   if (parmax>=56) then
      F1111   =  param(   56)
      F1112   =  param(   57)
      F1114   =  param(   58)
      F1115   =  param(   59)
      F1122   =  param(   60)
      F1123   =  param(   61)
      F1124   =  param(   62)
      F1126   =  param(   63)
      F1136   =  param(   64)
      F1144   =  param(   65)
      F1146   =  param(   66)
      F1155   =  param(   67)
      F1156   =  param(   68)
      F1222   =  param(   69)
      F1225   =  param(   70)
      F1233   =  param(   71)
      F1234   =  param(   72)
      F1235   =  param(   73)
      F1244   =  param(   74)
      F1245   =  param(   75)
      F1246   =  param(   76)
      F1256   =  param(   77)
      F1334   =  param(   78)
      F1335   =  param(   79)
      F1355   =  param(   80)
      F1366   =  param(   81)
      F1444   =  param(   82)
      F1445   =  param(   83)
      F1456   =  param(   84)
      F1466   =  param(   85)
      F1566   =  param(   86)
      F1666   =  param(   87)
      F2222   =  param(   88)
      F2224   =  param(   89)
      F2225   =  param(   90)
      F2233   =  param(   91)
      F2244   =  param(   92)
      F2245   =  param(   93)
      F2246   =  param(   94)
      F2255   =  param(   95)
      F2256   =  param(   96)
      F2266   =  param(   97)
      F2333   =  param(   98)
      F2334   =  param(   99)
      F2335   =  param(  100)
      F2336   =  param(  101)
      F2344   =  param(  102)
      F2345   =  param(  103)
      F2356   =  param(  104)
      F2366   =  param(  105)
      F2444   =  param(  106)
      F2445   =  param(  107)
      F2455   =  param(  108)
      F2456   =  param(  109)
      F2466   =  param(  110)
      F3335   =  param(  111)
      F3445   =  param(  112)
      F3555   =  param(  113)
      F3556   =  param(  114)
      F3566   =  param(  115)
      F3666   =  param(  116)
      F4444   =  param(  117)
      F4445   =  param(  118)
      F4456   =  param(  119)
      F4466   =  param(  120)
      F4555   =  param(  121)
      F4556   =  param(  122)
      F5566   =  param(  123)
      F5666   =  param(  124)
      F6666   =  param(  125)
   endif
   if (parmax>=126) then
      F11111  =  param(  126)
      F11112  =  param(  127)
      F11114  =  param(  128)
      F11116  =  param(  129)
      F11123  =  param(  130)
      F11124  =  param(  131)
      F11133  =  param(  132)
      F11135  =  param(  133)
      F11136  =  param(  134)
      F11144  =  param(  135)
      F11145  =  param(  136)
      F11156  =  param(  137)
      F11166  =  param(  138)
      F11223  =  param(  139)
      F11224  =  param(  140)
      F11234  =  param(  141)
      F11236  =  param(  142)
      F11244  =  param(  143)
      F11245  =  param(  144)
      F11255  =  param(  145)
      F11266  =  param(  146)
      F11333  =  param(  147)
      F11335  =  param(  148)
      F11336  =  param(  149)
      F11345  =  param(  150)
      F11356  =  param(  151)
      F11444  =  param(  152)
      F11446  =  param(  153)
      F11455  =  param(  154)
      F11456  =  param(  155)
      F11555  =  param(  156)
      F11556  =  param(  157)
      F12225  =  param(  158)
      F12233  =  param(  159)
      F12236  =  param(  160)
      F12255  =  param(  161)
      F12266  =  param(  162)
      F12333  =  param(  163)
      F12334  =  param(  164)
      F12336  =  param(  165)
      F12344  =  param(  166)
      F12345  =  param(  167)
      F12355  =  param(  168)
      F12356  =  param(  169)
      F12444  =  param(  170)
      F12455  =  param(  171)
      F12466  =  param(  172)
      F12666  =  param(  173)
      F13333  =  param(  174)
      F13334  =  param(  175)
      F13335  =  param(  176)
      F13344  =  param(  177)
      F13345  =  param(  178)
      F13346  =  param(  179)
      F13356  =  param(  180)
      F13445  =  param(  181)
      F13446  =  param(  182)
      F13456  =  param(  183)
      F13556  =  param(  184)
      F13566  =  param(  185)
      F13666  =  param(  186)
      F14444  =  param(  187)
      F14445  =  param(  188)
      F14456  =  param(  189)
      F14466  =  param(  190)
      F14556  =  param(  191)
      F14666  =  param(  192)
      F15555  =  param(  193)
      F15556  =  param(  194)
      F15566  =  param(  195)
      F22224  =  param(  196)
      F22333  =  param(  197)
      F22334  =  param(  198)
      F22336  =  param(  199)
      F22346  =  param(  200)
      F22446  =  param(  201)
      F22456  =  param(  202)
      F23333  =  param(  203)
      F23334  =  param(  204)
      F23335  =  param(  205)
      F23336  =  param(  206)
      F23344  =  param(  207)
      F23346  =  param(  208)
      F23355  =  param(  209)
      F23356  =  param(  210)
      F23366  =  param(  211)
      F23444  =  param(  212)
      F23446  =  param(  213)
      F23456  =  param(  214)
      F23466  =  param(  215)
      F23555  =  param(  216)
      F23556  =  param(  217)
      F33333  =  param(  218)
      F33335  =  param(  219)
      F33336  =  param(  220)
      F33344  =  param(  221)
      F33345  =  param(  222)
      F33346  =  param(  223)
      F33355  =  param(  224)
      F33356  =  param(  225)
      F33366  =  param(  226)
      F33444  =  param(  227)
      F33446  =  param(  228)
      F33455  =  param(  229)
      F33466  =  param(  230)
      F33555  =  param(  231)
      F33556  =  param(  232)
      F33566  =  param(  233)
      F33666  =  param(  234)
      F34444  =  param(  235)
      F34445  =  param(  236)
      F34446  =  param(  237)
      F34455  =  param(  238)
      F34456  =  param(  239)
      F34466  =  param(  240)
      F34555  =  param(  241)
      F34556  =  param(  242)
      F34566  =  param(  243)
      F34666  =  param(  244)
      F35555  =  param(  245)
      F35556  =  param(  246)
      F35566  =  param(  247)
      F35666  =  param(  248)
      F36666  =  param(  249)
      F44444  =  param(  250)
      F44445  =  param(  251)
      F44456  =  param(  252)
      F44466  =  param(  253)
      F44566  =  param(  254)
      F44666  =  param(  255)
      F45566  =  param(  256)
      F45666  =  param(  257)
      F46666  =  param(  258)
      F55555  =  param(  259)
      F55666  =  param(  260)
      F56666  =  param(  261)
   endif
   if (parmax>=262) then
      F111111 =  param(  262)
      F111112 =  param(  263)
      F111114 =  param(  264)
      F111116 =  param(  265)
      F111123 =  param(  266)
      F111124 =  param(  267)
      F111126 =  param(  268)
      F111133 =  param(  269)
      F111136 =  param(  270)
      F111144 =  param(  271)
      F111145 =  param(  272)
      F111156 =  param(  273)
      F111166 =  param(  274)
      F111222 =  param(  275)
      F111223 =  param(  276)
      F111234 =  param(  277)
      F111236 =  param(  278)
      F111246 =  param(  279)
      F111255 =  param(  280)
      F111266 =  param(  281)
      F111334 =  param(  282)
      F111335 =  param(  283)
      F111336 =  param(  284)
      F111344 =  param(  285)
      F111346 =  param(  286)
      F111356 =  param(  287)
      F111444 =  param(  288)
      F111446 =  param(  289)
      F111456 =  param(  290)
      F111466 =  param(  291)
      F111556 =  param(  292)
      F111666 =  param(  293)
      F112233 =  param(  294)
      F112234 =  param(  295)
      F112235 =  param(  296)
      F112236 =  param(  297)
      F112244 =  param(  298)
      F112245 =  param(  299)
      F112246 =  param(  300)
      F112255 =  param(  301)
      F112256 =  param(  302)
      F112266 =  param(  303)
      F112333 =  param(  304)
      F112344 =  param(  305)
      F112345 =  param(  306)
      F112356 =  param(  307)
      F112366 =  param(  308)
      F112444 =  param(  309)
      F112445 =  param(  310)
      F112446 =  param(  311)
      F112455 =  param(  312)
      F112456 =  param(  313)
      F112466 =  param(  314)
      F112555 =  param(  315)
      F112556 =  param(  316)
      F112566 =  param(  317)
      F112666 =  param(  318)
      F113333 =  param(  319)
      F113334 =  param(  320)
      F113335 =  param(  321)
      F113336 =  param(  322)
      F114444 =  param(  323)
      F114446 =  param(  324)
      F114456 =  param(  325)
      F114466 =  param(  326)
      F114556 =  param(  327)
      F114666 =  param(  328)
      F115555 =  param(  329)
      F115566 =  param(  330)
      F115666 =  param(  331)
      F122222 =  param(  332)
      F122223 =  param(  333)
      F122224 =  param(  334)
      F122225 =  param(  335)
      F122226 =  param(  336)
      F122234 =  param(  337)
      F122235 =  param(  338)
      F122236 =  param(  339)
      F122244 =  param(  340)
      F122245 =  param(  341)
      F122246 =  param(  342)
      F122255 =  param(  343)
      F122256 =  param(  344)
      F122266 =  param(  345)
      F122333 =  param(  346)
      F122334 =  param(  347)
      F122336 =  param(  348)
      F122344 =  param(  349)
      F122345 =  param(  350)
      F122346 =  param(  351)
      F122444 =  param(  352)
      F122445 =  param(  353)
      F122446 =  param(  354)
      F122455 =  param(  355)
      F122456 =  param(  356)
      F122466 =  param(  357)
      F122555 =  param(  358)
      F122666 =  param(  359)
      F123355 =  param(  360)
      F123356 =  param(  361)
      F123366 =  param(  362)
      F123444 =  param(  363)
      F123446 =  param(  364)
      F123456 =  param(  365)
      F123466 =  param(  366)
      F123555 =  param(  367)
      F123566 =  param(  368)
      F124444 =  param(  369)
      F124445 =  param(  370)
      F124446 =  param(  371)
      F124455 =  param(  372)
      F125555 =  param(  373)
      F125556 =  param(  374)
      F125566 =  param(  375)
      F125666 =  param(  376)
      F126666 =  param(  377)
      F133556 =  param(  378)
      F133566 =  param(  379)
      F134455 =  param(  380)
      F134456 =  param(  381)
      F134555 =  param(  382)
      F134556 =  param(  383)
      F134566 =  param(  384)
      F134666 =  param(  385)
      F144444 =  param(  386)
      F144446 =  param(  387)
      F144456 =  param(  388)
      F144466 =  param(  389)
      F144566 =  param(  390)
      F144666 =  param(  391)
      F145566 =  param(  392)
      F145666 =  param(  393)
      F146666 =  param(  394)
      F155556 =  param(  395)
      F155666 =  param(  396)
      F166666 =  param(  397)
      F222222 =  param(  398)
      F222223 =  param(  399)
      F222224 =  param(  400)
      F222225 =  param(  401)
      F222226 =  param(  402)
      F222234 =  param(  403)
      F222236 =  param(  404)
      F222244 =  param(  405)
      F222245 =  param(  406)
      F222246 =  param(  407)
      F222255 =  param(  408)
      F222256 =  param(  409)
      F222266 =  param(  410)
      F222333 =  param(  411)
      F222334 =  param(  412)
      F222335 =  param(  413)
      F222336 =  param(  414)
      F222344 =  param(  415)
      F222345 =  param(  416)
      F222346 =  param(  417)
      F222355 =  param(  418)
      F222356 =  param(  419)
      F222366 =  param(  420)
      F222444 =  param(  421)
      F222445 =  param(  422)
      F222446 =  param(  423)
      F222455 =  param(  424)
      F222456 =  param(  425)
      F222466 =  param(  426)
      F222555 =  param(  427)
      F222556 =  param(  428)
      F222566 =  param(  429)
      F222666 =  param(  430)
      F223333 =  param(  431)
      F223344 =  param(  432)
      F223346 =  param(  433)
      F223356 =  param(  434)
      F223366 =  param(  435)
      F223444 =  param(  436)
      F223445 =  param(  437)
      F223446 =  param(  438)
      F223456 =  param(  439)
      F223556 =  param(  440)
      F223566 =  param(  441)
      F223666 =  param(  442)
      F224444 =  param(  443)
      F224445 =  param(  444)
      F224446 =  param(  445)
      F224455 =  param(  446)
      F224456 =  param(  447)
      F224466 =  param(  448)
      F224555 =  param(  449)
      F224556 =  param(  450)
      F224566 =  param(  451)
      F224666 =  param(  452)
      F225555 =  param(  453)
      F225556 =  param(  454)
      F225566 =  param(  455)
      F225666 =  param(  456)
      F226666 =  param(  457)
      F233336 =  param(  458)
      F233455 =  param(  459)
      F233466 =  param(  460)
      F233666 =  param(  461)
      F234444 =  param(  462)
      F234446 =  param(  463)
      F234456 =  param(  464)
      F234466 =  param(  465)
      F234566 =  param(  466)
      F234666 =  param(  467)
      F235555 =  param(  468)
      F235556 =  param(  469)
      F235566 =  param(  470)
      F244444 =  param(  471)
      F244445 =  param(  472)
      F244446 =  param(  473)
      F244455 =  param(  474)
      F244456 =  param(  475)
      F244466 =  param(  476)
      F244555 =  param(  477)
      F244556 =  param(  478)
      F244566 =  param(  479)
      F244666 =  param(  480)
      F245555 =  param(  481)
      F245556 =  param(  482)
      F245566 =  param(  483)
      F245666 =  param(  484)
      F246666 =  param(  485)
      F255555 =  param(  486)
      F255556 =  param(  487)
      F255566 =  param(  488)
      F255666 =  param(  489)
      F256666 =  param(  490)
      F266666 =  param(  491)
      F444444 =  param(  492)
      F444446 =  param(  493)
      F444456 =  param(  494)
      F444466 =  param(  495)
      F444566 =  param(  496)
      F444666 =  param(  497)
      F445555 =  param(  498)
      F445566 =  param(  499)
      F445666 =  param(  500)
      F455666 =  param(  501)
      F456666 =  param(  502)
      F466666 =  param(  503)
      F555556 =  param(  504)
      F555566 =  param(  505)
      F555666 =  param(  506)
      F666666 =  param(  507)
   endif



!  define local coordinates

   r14 = local(1)
   r24 = local(2)
   r34 = local(3)
   !
   alpha1 = local(4)
   alpha2 = local(5)
   alpha3 = local(6)

!  define coordinates xi

   ta1 = cos(alpha1)-cos(alphae)
   ta2 = cos(alpha2)-cos(alphae)
   ta3 = cos(alpha3)-cos(alphae)

     select case ( ix )
     case default
       write (6,"(' dip. order component',i8)") ix
       stop 'dip. order component'
     case (1)
       xi1=(r14-re14) *exp(-beta*(r14-re14)**2)
       xi2=(r24-re14) *exp(-beta*(r24-re14)**2)
       xi3=(r34-re14) *exp(-beta*(r34-re14)**2)

       xi4=ta1
       xi5=ta2
       xi6=ta3
       !
     case (2)
       xi1=(r24-re14) *exp(-beta*(r24-re14)**2)
       xi2=(r34-re14) *exp(-beta*(r34-re14)**2)
       xi3=(r14-re14) *exp(-beta*(r14-re14)**2)

       xi4=ta2
       xi5=ta3
       xi6=ta1
       ! 
     case (3)
       xi1=(r34-re14) *exp(-beta*(r34-re14)**2)
       xi2=(r14-re14) *exp(-beta*(r14-re14)**2)
       xi3=(r24-re14) *exp(-beta*(r24-re14)**2)

       xi4=ta3
       xi5=ta1
       xi6=ta2
       ! 
     end select


     t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0
     if (parmax>=1) then 
      t1 = F1*xi1+(xi2+xi3)*F3+(xi6+xi5)*F5+F4*xi4
     endif 

     if (parmax>=11) then 
      t2 = F11*xi1**2+F14*xi1*xi4+F23*xi2*xi3+(xi1*xi3+xi1*xi2)*F13+(xi1*xi5+xi1*xi6)*F16+&
        (xi2*xi4+xi3*xi4)*F34+(xi3**2+xi2**2)*F33+(xi3*xi6+xi2*xi5)*F36+(xi2*xi6+xi3*xi5)*F35+&
        F44*xi4**2+(xi4*xi6+xi4*xi5)*F46+(xi6**2+xi5**2)*F55+F56*xi5*xi6
     endif                                                              
                                                                        
     if (parmax>=24) then                                                 
        s1 = (xi1**2*xi3+xi1**2*xi2)*F112+(xi3*xi4*xi6+xi2*xi4*xi5)*F245+&
        (xi2*xi5*xi6+xi3*xi6*xi5)*F256+(xi1**2*xi5+xi1**2*xi6)*F115+&
        (xi3*xi5**2+xi2*xi6**2)*F266+(xi2**2*xi6+xi3**2*xi5)*F226+&
        (xi1*xi4*xi6+xi1*xi4*xi5)*F146+(xi1*xi6**2+xi1*xi5**2)*F155+&
        (xi2*xi5**2+xi3*xi6**2)*F255+(xi1*xi2*xi5+xi1*xi3*xi6)*F136+&
        (xi1*xi3**2+xi1*xi2**2)*F133+F144*xi1*xi4**2+F114*xi1**2*xi4+&
        F234*xi2*xi3*xi4+F123*xi1*xi2*xi3+F456*xi4*xi5*xi6
        t3 = s1+(xi1*xi3*xi4+xi1*xi2*xi4)*F124+(xi6**3+xi5**3)*F555+&
        (xi1*xi3*xi5+xi1*xi2*xi6)*F135+(xi3*xi4*xi5+xi2*xi4*xi6)*F246+&
        (xi2**2*xi4+xi3**2*xi4)*F334+(xi3**2*xi6+xi2**2*xi5)*F225+&
        (xi3*xi4**2+xi2*xi4**2)*F344+F156*xi1*xi5*xi6+(xi3*xi2*xi6+xi2*xi3*xi5)*F235+&
        (xi3**3+xi2**3)*F333+F444*xi4**3+(xi4*xi5**2+xi4*xi6**2)*F466+&
        (xi3**2*xi2+xi2**2*xi3)*F223+F111*xi1**3+(xi5**2*xi6+xi6**2*xi5)*F556+&
        (xi4**2*xi6+xi4**2*xi5)*F445
     endif                                                              
                                                                        
     if (parmax>=56) then                                                       
        s2 = (xi1*xi2**3+xi1*xi3**3)*F1222+ (xi3*xi2**3+xi2*xi3**3)*F2333+ &
        (xi3*xi5**3+xi2*xi6**3)*F3555+                                     &
        (xi1*xi2*xi6**2+xi1*xi3*xi5**2)*F1355+                             &
        (xi3*xi6**3+xi2*xi5**3)*F3666+                                     &
        (xi2*xi4*xi5*xi6+xi3*xi4*xi6*xi5)*F2456+                           &
        (xi2**3*xi6+xi3**3*xi5)*F3335+                                     &
        (xi3*xi4**2*xi6+xi2*xi4**2*xi5)*F2445+                             &
        (xi1*xi6*xi5**2+xi1*xi5*xi6**2)*F1566+                             &
        (xi4**2*xi5**2+xi4**2*xi6**2)*F4466+                               &
        (xi3*xi4*xi6**2+xi2*xi4*xi5**2)*F2455+ (xi5**4+xi6**4)*F6666+      &
        (xi1*xi3*xi6*xi5+xi1*xi2*xi5*xi6)*F1256+                           &
        (xi1*xi3*xi2*xi6+xi1*xi2*xi3*xi5)*F1235+                           &
        (xi1**2*xi3*xi6+xi1**2*xi2*xi5)*F1136+                             &
        (xi3**3*xi4+xi2**3*xi4)*F2224+ (xi4*xi5**3+xi4*xi6**3)*F4555       
        
        s1 = s2+ (xi2*xi3*xi4*xi5+xi3*xi2*xi4*xi6)*F2345+                  &
        (xi2*xi3**2*xi4+xi3*xi2**2*xi4)*F2334+                             &
        (xi1*xi3*xi2**2+xi1*xi2*xi3**2)*F1233+                             &
        (xi1*xi3*xi4*xi5+xi1*xi2*xi4*xi6)*F1246+                           &
        (xi1*xi2*xi4*xi5+xi1*xi3*xi4*xi6)*F1245+                           &
        (xi2**2*xi6**2+xi3**2*xi5**2)*F2266+                               &
        (xi1*xi3**2*xi6+xi1*xi2**2*xi5)*F1225+                             &
        (xi2**2*xi4*xi6+xi3**2*xi4*xi5)*F2246 +F1444*xi1*xi4**3            &
        +F5566*xi5**2*xi6**2 +F2233*xi2**2*xi3**2 +F1114*xi1**3*xi4+       &
        (xi3*xi4**2*xi5+xi2*xi4**2*xi6)*F3445+                             &
        (xi1**2*xi3*xi4+xi1**2*xi2*xi4)*F1124+                             &
        (xi1**3*xi6+xi1**3*xi5)*F1115+ (xi3**4+xi2**4)*F2222+              &
        (xi2**2*xi4**2+xi3**2*xi4**2)*F2244 +F4444*xi4**4 

        s2 = s1          &
        +F1144*xi1**2*xi4**2+ (xi2*xi6*xi5**2+xi3*xi5*xi6**2)*F3566+       &
        (xi1*xi3**2*xi5+xi1*xi2**2*xi6)*F1335+                             &
        (xi6*xi5**3+xi5*xi6**3)*F5666+                                     &
        (xi1*xi4*xi6**2+xi1*xi4*xi5**2)*F1466 +F2356*xi2*xi3*xi5*xi6       &
        +F1234*xi1*xi2*xi3*xi4+ (xi1*xi2*xi5**2+xi1*xi3*xi6**2)*F1366+     &
        (xi1**2*xi4*xi6+xi1**2*xi4*xi5)*F1146+                             &
        (xi1*xi4**2*xi5+xi1*xi4**2*xi6)*F1445+                             &
        (xi3*xi2*xi5**2+xi2*xi3*xi6**2)*F2366+                             &
        (xi3*xi5**2*xi6+xi2*xi6**2*xi5)*F3556 +F1456*xi1*xi4*xi5*xi6+      &
        (xi2*xi3**2*xi6+xi3*xi2**2*xi5)*F2336 +F1111*xi1**4+               &
        (xi4*xi5**2*xi6+xi4*xi6**2*xi5)*F4556+                             &
        (xi3*xi2**2*xi6+xi2*xi3**2*xi5)*F2335 

        t4 = s2+         &
        (xi2*xi4**3+xi3*xi4**3)*F2444+ (xi1*xi5**3+xi1*xi6**3)*F1666+      &
        (xi1*xi2**2*xi4+xi1*xi3**2*xi4)*F1334+                             &
        (xi2**2*xi4*xi5+xi3**2*xi4*xi6)*F2245+                             &
        (xi1**2*xi3**2+xi1**2*xi2**2)*F1122 +F4456*xi4**2*xi5*xi6+         &
        (xi1**2*xi5**2+xi1**2*xi6**2)*F1155+                               &
        (xi2*xi4*xi6**2+xi3*xi4*xi5**2)*F2466+                             &
        (xi1*xi2*xi4**2+xi1*xi3*xi4**2)*F1244+                             &
        (xi2**3*xi5+xi3**3*xi6)*F2225+ (xi4**3*xi5+xi4**3*xi6)*F4445       &
        +F1156*xi1**2*xi5*xi6+ (xi1**2*xi2*xi6+xi1**2*xi3*xi5)*F1126+      &
        (xi1**3*xi3+xi1**3*xi2)*F1112+                                     &
        (xi2**2*xi5*xi6+xi3**2*xi6*xi5)*F2256 +F2344*xi2*xi3*xi4**2+       &
        (xi2**2*xi5**2+xi3**2*xi6**2)*F2255 +F1123*xi1**2*xi2*xi3           
     endif                                                                  


   if (parmax>=126) then
     s3 = (xi3*xi4**2*xi5*xi6+ xi2*xi4**2*xi6*xi5)*F34456+                 &
     (xi3**2*xi5*xi6**2+ xi2**2*xi6*xi5**2)*F33566+ (xi2*xi4**3*xi5+       &
     xi3*xi4**3*xi6)*F34446+ (xi1**3*xi3*xi4+ xi1**3*xi2*xi4)*F11124+      &
     (xi2**3*xi6*xi5+ xi3**3*xi5*xi6)*F33356+ (xi3**2*xi6**3+              &
     xi2**2*xi5**3)*F33666+ (xi3**2*xi2**2*xi5+                            &
     xi2**2*xi3**2*xi6)*F22336+ (xi4**2*xi6*xi5**2+                        &
     xi4**2*xi5*xi6**2)*F44566+ (xi3**4*xi5+ xi2**4*xi6)*F33335+           &
     (xi3*xi2**2*xi5**2+ xi2*xi3**2*xi6**2)*F23366+                        &
     (xi1*xi2*xi6**2*xi5+ xi1*xi3*xi5**2*xi6)*F13556+                      &
     (xi3**2*xi5**2*xi6+ xi2**2*xi6**2*xi5)*F33556+                        &
     (xi3*xi2*xi4**2*xi5+ xi2*xi3*xi4**2*xi6)*F23446+                      &
     (xi3*xi4*xi5*xi6**2+ xi2*xi4*xi6*xi5**2)*F34566+                      &
     F11123*xi1**3*xi2*xi3+ F15566*xi1*xi5**2*xi6**2+                      &
     F12233*xi1*xi2**2*xi3**2                                              
     s2 = s3+ F22334*xi2**2*xi3**2*xi4+ (xi3**2*xi5**3+                    &
     xi2**2*xi6**3)*F33555+ (xi2**2*xi4*xi5**2+                            &
     xi3**2*xi4*xi6**2)*F33466+ F44456*xi4**3*xi5*xi6+                     &
     F11444*xi1**2*xi4**3+ (xi3*xi2*xi4*xi5**2+                            &
     xi2*xi3*xi4*xi6**2)*F23466+ F23444*xi2*xi3*xi4**3+                    &
     F45566*xi4*xi5**2*xi6**2+ F11156*xi1**3*xi5*xi6+ (xi6*xi5**4+         &
     xi5*xi6**4)*F56666+ (xi3*xi5**4+ xi2*xi6**4)*F35555+                  &
     (xi3**3*xi6**2+ xi2**3*xi5**2)*F33366+ F11144*xi1**3*xi4**2+          &
     (xi2*xi4*xi6**3+ xi3*xi4*xi5**3)*F34555+ F11114*xi1**4*xi4+           &
     F14444*xi1*xi4**4+ (xi2*xi5**4+ xi3*xi6**4)*F36666
     s3 = (xi3*xi4**2*xi6**2+ xi2*xi4**2*xi5**2)*F34466+                   &
     (xi2**2*xi4*xi6**2+ xi3**2*xi4*xi5**2)*F33455+                        &
     (xi3*xi2**2*xi4**2+ xi2*xi3**2*xi4**2)*F23344+ (xi2*xi4*xi5**3+       &
     xi3*xi4*xi6**3)*F34666+ (xi1*xi2**2*xi5**2+                           &
     xi1*xi3**2*xi6**2)*F12255+ (xi3*xi2**3*xi5+                           &
     xi2*xi3**3*xi6)*F23336+ (xi2**5+ xi3**5)*F33333+                      &
     (xi1*xi3*xi5*xi6**2+ xi1*xi2*xi6*xi5**2)*F13566+                      &
     (xi2**3*xi4*xi6+ xi3**3*xi4*xi5)*F33345+ (xi1*xi2**3*xi6+             &
     xi1*xi3**3*xi5)*F13335+ (xi2*xi6**3*xi5+ xi3*xi5**3*xi6)*F35556+      &
     (xi3*xi2**2*xi6*xi5+ xi2*xi3**2*xi5*xi6)*F23356+                      &
     (xi3*xi2*xi6**2*xi5+ xi2*xi3*xi5**2*xi6)*F23556+                      &
     (xi1*xi3**2*xi5**2+ xi1*xi2**2*xi6**2)*F12266+ (xi2*xi4**4+           &
     xi3*xi4**4)*F34444+ (xi6**5+ xi5**5)*F55555+ s2
     s4 = s3+ (xi3**2*xi4**2*xi6+ xi2**2*xi4**2*xi5)*F33446+               &
     (xi2**2*xi4*xi5*xi6+ xi3**2*xi4*xi6*xi5)*F22456+                      &
     (xi2*xi6*xi5**3+ xi3*xi5*xi6**3)*F35666+ (xi3*xi5**2*xi6**2+          &
     xi2*xi6**2*xi5**2)*F35566+ (xi3*xi2**3*xi6+                           &
     xi2*xi3**3*xi5)*F23335+ (xi1**4*xi3+ xi1**4*xi2)*F11112+              &
     (xi2**3*xi4*xi5+ xi3**3*xi4*xi6)*F33346+ (xi3*xi4*xi5**2*xi6+         &
     xi2*xi4*xi6**2*xi5)*F34556
     s1 = s4+ (xi2*xi3**2*xi4*xi6+ xi3*xi2**2*xi4*xi5)*F23346+             &
     (xi1**3*xi2*xi5+ xi1**3*xi3*xi6)*F11136+ (xi1**3*xi6**2+              &
     xi1**3*xi5**2)*F11166+ (xi1*xi2**2*xi4*xi6+                           &
     xi1*xi3**2*xi4*xi5)*F13345+ (xi1*xi3**4+ xi1*xi2**4)*F13333+          &
     (xi3*xi4**2*xi5**2+ xi2*xi4**2*xi6**2)*F34455+                        &
     (xi1*xi2*xi3**2*xi6+ xi1*xi3*xi2**2*xi5)*F12336+                      &
     (xi1**2*xi2**2*xi3+ xi1**2*xi3**2*xi2)*F11223+                        &
     (xi1**2*xi4*xi5**2+ xi1**2*xi4*xi6**2)*F11455+                        &
     (xi1*xi3*xi4*xi5*xi6+ xi1*xi2*xi4*xi6*xi5)*F13456
      s3 = s1+ (xi1*xi3**2*xi4**2+ xi1*xi2**2*xi4**2)*F13344+              &
     (xi1**2*xi2*xi4*xi5+ xi1**2*xi3*xi4*xi6)*F11245+ F11111*xi1**5+       &
     (xi1**2*xi2**2*xi6+ xi1**2*xi3**2*xi5)*F11335+                        &
     (xi1*xi3*xi2*xi4*xi6+ xi1*xi2*xi3*xi4*xi5)*F12345+                    &
     (xi1**2*xi3**2*xi6+ xi1**2*xi2**2*xi5)*F11336+ (xi1*xi2*xi4**3+       &
     xi1*xi3*xi4**3)*F12444+ (xi1*xi4*xi6**2*xi5+                          &
     xi1*xi4*xi5**2*xi6)*F14556+ F44444*xi4**5+ (xi1*xi2*xi4**2*xi5+       &
     xi1*xi3*xi4**2*xi6)*F13446+ (xi1*xi3*xi2**3+                          &
     xi1*xi2*xi3**3)*F12333+ (xi4**3*xi5**2+ xi4**3*xi6**2)*F44466+        &
     (xi3**2*xi4**2*xi5+ xi2**2*xi4**2*xi6)*F22446+                        &
     F12356*xi1*xi2*xi3*xi5*xi6+ (xi4**4*xi5+ xi4**4*xi6)*F44445+          &
     (xi1**4*xi6+ xi1**4*xi5)*F11116
     s4 = s3+ (xi1**2*xi3*xi2*xi5+ xi1**2*xi2*xi3*xi6)*F11236+             &
     (xi2**4*xi4+ xi3**4*xi4)*F22224+ (xi1*xi3*xi2*xi6**2+                 &
     xi1*xi2*xi3*xi5**2)*F12355+ (xi1**2*xi2*xi6**2+                       &
     xi1**2*xi3*xi5**2)*F11266+ (xi1**2*xi4**2*xi6+                        &
     xi1**2*xi4**2*xi5)*F11446+ F23456*xi2*xi3*xi4*xi5*xi6+                &
     (xi1**2*xi2*xi5**2+ xi1**2*xi3*xi6**2)*F11255+                        &
     (xi1**2*xi3*xi4*xi5+ xi1**2*xi2*xi4*xi6)*F11345
     s2 = s4+ (xi1*xi3*xi4*xi6**2+ xi1*xi2*xi4*xi5**2)*F12455+             &
     F12344*xi1*xi2*xi3*xi4**2+ (xi4*xi5*xi6**3+                           &
     xi4*xi6*xi5**3)*F45666+ (xi1*xi4*xi6**3+ xi1*xi4*xi5**3)*F14666+      &
     (xi1*xi5**4+ xi1*xi6**4)*F15555+ (xi1*xi2**2*xi3*xi6+                 &
     xi1*xi3**2*xi2*xi5)*F12236+ (xi1**2*xi6**2*xi5+                       &
     xi1**2*xi5**2*xi6)*F11556+ (xi2*xi3**2*xi5**2+                        &
     xi3*xi2**2*xi6**2)*F23355+ (xi1*xi2*xi5**3+                           &
     xi1*xi3*xi6**3)*F13666
     s4 = s2+ (xi1*xi3*xi4**2*xi5+ xi1*xi2*xi4**2*xi6)*F13445+             &
     (xi1*xi4**2*xi6**2+ xi1*xi4**2*xi5**2)*F14466+                        &
     F14456*xi1*xi4**2*xi5*xi6+ F11456*xi1**2*xi4*xi5*xi6+                 &
     (xi1**2*xi3*xi4**2+ xi1**2*xi2*xi4**2)*F11244+                        &
     F11234*xi1**2*xi2*xi3*xi4+ (xi1*xi2**2*xi6*xi5+                       &
     xi1*xi3**2*xi5*xi6)*F13356+ (xi1**2*xi5**3+                           &
     xi1**2*xi6**3)*F11555
     s3 = s4+ (xi1*xi2**3*xi5+ xi1*xi3**3*xi6)*F12225+                     &
     (xi1*xi3**3*xi4+ xi1*xi2**3*xi4)*F13334+ (xi1*xi5**3*xi6+             &
     xi1*xi6**3*xi5)*F15556+ (xi1*xi2*xi4*xi6**2+                          &
     xi1*xi3*xi4*xi5**2)*F12466+ (xi1**2*xi2**2*xi4+                       &
     xi1**2*xi3**2*xi4)*F11224+ (xi1**3*xi2*xi6+                           &
     xi1**3*xi3*xi5)*F11135+ (xi1**2*xi3*xi5*xi6+                          &
     xi1**2*xi2*xi6*xi5)*F11356+ (xi1*xi4**3*xi6+                          &
     xi1*xi4**3*xi5)*F14445+ (xi1**2*xi3**3+ xi1**2*xi2**3)*F11333
     s4 = s3+ (xi1*xi2*xi6**3+ xi1*xi3*xi5**3)*F12666+                     &
     (xi1*xi2**2*xi4*xi5+ xi1*xi3**2*xi4*xi6)*F13346+ (xi3**3*xi5**2+      &
     xi2**3*xi6**2)*F33355+ (xi1*xi3*xi2**2*xi4+                           &
     xi1*xi2*xi3**2*xi4)*F12334+ (xi1**3*xi3**2+                           &
     xi1**3*xi2**2)*F11133+ (xi2**2*xi4**3+ xi3**2*xi4**3)*F33444+         &
     (xi3**2*xi2*xi4*xi5+ xi2**2*xi3*xi4*xi6)*F22346+ (xi5**2*xi6**3+      &
     xi6**2*xi5**3)*F55666
     t5 = s4+ (xi2**3*xi4**2+ xi3**3*xi4**2)*F33344+                       &
     (xi4**2*xi6**3+ xi4**2*xi5**3)*F44666+ (xi1**3*xi4*xi6+               &
     xi1**3*xi4*xi5)*F11145+ (xi3*xi2**3*xi4+ xi2*xi3**3*xi4)*F23334+      &
     (xi3*xi2**4+ xi2*xi3**4)*F23333+ (xi4*xi6**4+                         &
     xi4*xi5**4)*F46666+ (xi3**2*xi2**3+ xi2**2*xi3**3)*F22333+            &
     (xi3**4*xi6+ xi2**4*xi5)*F33336+ (xi3*xi4**3*xi5+                     &
     xi2*xi4**3*xi6)*F34445+ (xi2*xi3*xi5**3+ xi3*xi2*xi6**3)*F23555
   endif
                                                                           &
   if(parmax>=262) then                                                        
     s4 = (xi1*xi2**2*xi4**3+ xi1*xi3**2*xi4**3)*F122444+                  &
     (xi2**2*xi3*xi4**2*xi5+ xi3**2*xi2*xi4**2*xi6)*F223445+               &
     (xi3**2*xi2**2*xi4*xi5+ xi2**2*xi3**2*xi4*xi6)*F223346+               &
     (xi1**3*xi3**2*xi6+ xi1**3*xi2**2*xi5)*F111336+                       &
     (xi1*xi3**3*xi4**2+ xi1*xi2**3*xi4**2)*F122244+                       &
     (xi3**2*xi4*xi6**2*xi5+ xi2**2*xi4*xi5**2*xi6)*F224556+               &
     (xi3**3*xi5**3+ xi2**3*xi6**3)*F222666+ (xi1*xi3*xi2*xi6*xi5**2+      &
     xi1*xi2*xi3*xi5*xi6**2)*F123566+ (xi1*xi3**2*xi2**2*xi5+              &
     xi1*xi2**2*xi3**2*xi6)*F122336+ F234444*xi2*xi3*xi4**4+               &
     (xi2**3*xi3*xi4**2+ xi3**3*xi2*xi4**2)*F222344+                       &
     (xi1*xi4*xi5*xi6**3+ xi1*xi4*xi6*xi5**3)*F145666+                     &
     (xi1*xi4**3*xi6**2+ xi1*xi4**3*xi5**2)*F144466+                       &
     (xi1**2*xi3**3*xi5+ xi1**2*xi2**3*xi6)*F113335+                       &
     (xi2*xi3*xi4**2*xi6**2+ xi3*xi2*xi4**2*xi5**2)*F234466
     s3 = s4+ (xi1*xi3**2*xi5**2*xi6+                                      &
     xi1*xi2**2*xi6**2*xi5)*F133556+ F115566*xi1**2*xi5**2*xi6**2+         &
     F223344*xi2**2*xi3**2*xi4**2+ (xi1*xi3*xi4**2*xi5**2+                 &
     xi1*xi2*xi4**2*xi6**2)*F134455+ (xi1*xi4**2*xi5**3+                   &
     xi1*xi4**2*xi6**3)*F144666+ F444456*xi4**4*xi5*xi6+                   &
     (xi2*xi4**2*xi5**3+ xi3*xi4**2*xi6**3)*F244555+                       &
     F111123*xi1**4*xi2*xi3+ F111156*xi1**4*xi5*xi6+                       &
     F112233*xi1**2*xi2**2*xi3**2+ (xi1**2*xi2**3*xi5+                     &
     xi1**2*xi3**3*xi6)*F113336+ (xi1**2*xi2*xi3*xi6**2+                   &
     xi1**2*xi3*xi2*xi5**2)*F112366+ (xi3*xi4**4*xi5+                      &
     xi2*xi4**4*xi6)*F244446+ (xi3**2*xi4*xi5**3+                          &
     xi2**2*xi4*xi6**3)*F224666+ (xi1*xi4**4*xi5+                          &
     xi1*xi4**4*xi6)*F144446 
     s4 = F445566*xi4**2*xi5**2*xi6**2+ (xi1*xi6**4*xi5+                   &
     xi1*xi5**4*xi6)*F155556+ (xi1*xi2**2*xi5**3+                          &
     xi1*xi3**2*xi6**3)*F122555+ F144444*xi1*xi4**5+                       &
     F222333*xi2**3*xi3**3+ (xi1**4*xi2*xi5+ xi1**4*xi3*xi6)*F111136+      &
     (xi3**3*xi2*xi4*xi5+ xi2**3*xi3*xi4*xi6)*F222346+                     &
     F111144*xi1**4*xi4**2+ F114444*xi1**2*xi4**4+                         &
     F555666*xi5**3*xi6**3+ F111114*xi1**5*xi4+                            &
     F111444*xi1**3*xi4**3+ (xi3**3*xi4**3+ xi2**3*xi4**3)*F222444+        &
     (xi1*xi2**2*xi6*xi5**2+ xi1*xi3**2*xi5*xi6**2)*F133566+               &
     (xi2**4*xi4**2+ xi3**4*xi4**2)*F222244+ (xi2**4*xi5*xi6+              &
     xi3**4*xi6*xi5)*F222256
     s2 = s4+ (xi3**3*xi4*xi5**2+ xi2**3*xi4*xi6**2)*F222466+              &
     (xi3*xi4*xi6**2*xi5**2+ xi2*xi4*xi5**2*xi6**2)*F245566+               &
     (xi3**2*xi2*xi4*xi6*xi5+ xi2**2*xi3*xi4*xi5*xi6)*F223456+             &
     (xi3**3*xi4*xi6**2+ xi2**3*xi4*xi5**2)*F222455+ (xi3**4*xi2*xi4+      &
     xi2**4*xi3*xi4)*F222234+ (xi3**3*xi2*xi4*xi6+                         &
     xi2**3*xi3*xi4*xi5)*F222345+ s3+ (xi3*xi2**2*xi4*xi5**2+              &
     xi2*xi3**2*xi4*xi6**2)*F233466+ (xi2**4*xi6**2+                       &
     xi3**4*xi5**2)*F222266+ (xi3*xi2*xi4*xi5**3+                          &
     xi2*xi3*xi4*xi6**3)*F234666+ F111111*xi1**6+                          &
     (xi1*xi3**2*xi4**2*xi5+ xi1*xi2**2*xi4**2*xi6)*F122446+               &
     (xi1*xi2**2*xi3*xi4*xi6+ xi1*xi3**2*xi2*xi4*xi5)*F122346+             &
     (xi4**4*xi5**2+ xi4**4*xi6**2)*F444466+ (xi1**4*xi3**2+               &
     xi1**4*xi2**2)*F111133+ (xi3**2*xi2**4+ xi2**2*xi3**4)*F223333
     s4 = s2+ (xi2**4*xi5**2+ xi3**4*xi6**2)*F222255+                      &
     (xi1**2*xi6*xi5**3+ xi1**2*xi5*xi6**3)*F115666+ (xi1**4*xi4*xi6+      &
     xi1**4*xi4*xi5)*F111145+ (xi3*xi4*xi6**4+                             &
     xi2*xi4*xi5**4)*F245555+ (xi1*xi3**2*xi5**3+                          &
     xi1*xi2**2*xi6**3)*F122666+ (xi1*xi2*xi3*xi4**2*xi6+                  &
     xi1*xi3*xi2*xi4**2*xi5)*F123446+ (xi2**2*xi4**3*xi5+                  &
     xi3**2*xi4**3*xi6)*F224445+ (xi4*xi5**2*xi6**3+                       &
     xi4*xi6**2*xi5**3)*F455666+ (xi1**2*xi2*xi5**2*xi6+                   &
     xi1**2*xi3*xi6**2*xi5)*F112556+ (xi1**2*xi3*xi4*xi5**2+               &
     xi1**2*xi2*xi4*xi6**2)*F112466+ (xi1**2*xi2*xi4**2*xi5+               &
     xi1**2*xi3*xi4**2*xi6)*F112445+ (xi1*xi2**3*xi3*xi4+                  &
     xi1*xi3**3*xi2*xi4)*F122234+ (xi1*xi3*xi2**2*xi6**2+                  &
     xi1*xi2*xi3**2*xi5**2)*F123355+ (xi1*xi2**2*xi4*xi5**2+               &
     xi1*xi3**2*xi4*xi6**2)*F122455
     s5 = s4+ (xi1*xi2*xi4**4+ xi1*xi3*xi4**4)*F124444+                    &
     (xi1*xi2**4*xi5+ xi1*xi3**4*xi6)*F122225+ (xi2**3*xi5*xi6**2+         &
     xi3**3*xi6*xi5**2)*F222566+ (xi1*xi3**4*xi5+                          &
     xi1*xi2**4*xi6)*F122226+ (xi1**2*xi3*xi4*xi6**2+                      &
     xi1**2*xi2*xi4*xi5**2)*F112455+ (xi2*xi3*xi4*xi5*xi6**2+              &
     xi3*xi2*xi4*xi6*xi5**2)*F234566+ (xi2**4*xi3*xi6+                     &
     xi3**4*xi2*xi5)*F222236
     s3 = s5+ (xi1**4*xi3*xi5+ xi1**4*xi2*xi6)*F111126+                    &
     (xi3*xi4*xi5**4+ xi2*xi4*xi6**4)*F246666+ (xi2*xi4**3*xi5**2+         &
     xi3*xi4**3*xi6**2)*F244455+ (xi1*xi3*xi2*xi6**3+                      &
     xi1*xi2*xi3*xi5**3)*F123555+ (xi1*xi3*xi2**2*xi6*xi5+                 &
     xi1*xi2*xi3**2*xi5*xi6)*F123356+ (xi2**2*xi4**4+                      &
     xi3**2*xi4**4)*F224444+ (xi1**2*xi4*xi6**2*xi5+                       &
     xi1**2*xi4*xi5**2*xi6)*F114556+ (xi1**2*xi3*xi6**3+                   &
     xi1**2*xi2*xi5**3)*F112555+ (xi1**2*xi3*xi2**3+                       &
     xi1**2*xi2*xi3**3)*F112333
     s5 = s3+ (xi1**2*xi4*xi6**3+ xi1**2*xi4*xi5**3)*F114666+              &
     (xi1*xi3*xi4**2*xi5*xi6+ xi1*xi2*xi4**2*xi6*xi5)*F134456+             &
     (xi1*xi2**2*xi3*xi4*xi5+ xi1*xi3**2*xi2*xi4*xi6)*F122345+             &
     (xi1**3*xi2*xi6*xi5+ xi1**3*xi3*xi5*xi6)*F111356+                     &
     (xi2*xi4*xi5**3*xi6+ xi3*xi4*xi6**3*xi5)*F245556+                     &
     (xi1**2*xi3*xi4**2*xi5+ xi1**2*xi2*xi4**2*xi6)*F112446+               &
     (xi1**5*xi5+ xi1**5*xi6)*F111116
     s4 = s5+ (xi3**3*xi2**2*xi4+ xi2**3*xi3**2*xi4)*F222334+              &
     (xi1**3*xi3*xi2*xi5+ xi1**3*xi2*xi3*xi6)*F111236+                     &
     (xi1**2*xi2**2*xi4*xi5+ xi1**2*xi3**2*xi4*xi6)*F112245+               &
     (xi3*xi6**2*xi5**3+ xi2*xi5**2*xi6**3)*F255666+                       &
     (xi1**3*xi3*xi4**2+ xi1**3*xi2*xi4**2)*F111344+                       &
     (xi1**2*xi3**2*xi2*xi5+ xi1**2*xi2**2*xi3*xi6)*F112236+               &
     (xi1**3*xi3*xi4*xi5+ xi1**3*xi2*xi4*xi6)*F111246+                     &
     (xi1*xi3*xi4*xi5*xi6**2+ xi1*xi2*xi4*xi6*xi5**2)*F134566
     s5 = s4+ (xi2**2*xi3*xi4**2*xi6+                                      &
     xi3**2*xi2*xi4**2*xi5)*F223446+ (xi1**2*xi3**2*xi4**2+                &
     xi1**2*xi2**2*xi4**2)*F112244+ (xi1**2*xi5**4+                        &
     xi1**2*xi6**4)*F115555+ (xi1**2*xi2*xi3*xi4*xi5+                      &
     xi1**2*xi3*xi2*xi4*xi6)*F112345+ (xi1*xi3*xi4**3*xi5+                 &
     xi1*xi2*xi4**3*xi6)*F124446+ F234456*xi2*xi3*xi4**2*xi5*xi6+          &
     (xi1**2*xi2*xi5*xi6**2+ xi1**2*xi3*xi6*xi5**2)*F112566+               &
     (xi2*xi6**5+ xi3*xi5**5)*F266666
     s1 = s5+ (xi3*xi6**5+ xi2*xi5**5)*F255555+ (xi3**3*xi6**3+            &
     xi2**3*xi5**3)*F222555+ (xi1*xi3*xi5**4+                              &
     xi1*xi2*xi6**4)*F126666+ (xi3*xi4*xi6*xi5**3+                         &
     xi2*xi4*xi5*xi6**3)*F245666+ (xi1*xi3**2*xi2**3+                      &
     xi1*xi2**2*xi3**3)*F122333+ (xi2**3*xi4*xi5*xi6+                      &
     xi3**3*xi4*xi6*xi5)*F222456+ (xi2*xi4**2*xi5*xi6**2+                  &
     xi3*xi4**2*xi6*xi5**2)*F244566+ (xi1*xi3*xi4**2*xi6**2+               &
     xi1*xi2*xi4**2*xi5**2)*F124455+ (xi1*xi3*xi6**3*xi5+                  &
     xi1*xi2*xi5**3*xi6)*F125556
     s4 = s1+ (xi1**3*xi4**2*xi6+ xi1**3*xi4**2*xi5)*F111446+              &
     (xi2**2*xi3*xi6**3+ xi3**2*xi2*xi5**3)*F223666+                       &
     (xi1**2*xi3**2*xi6**2+ xi1**2*xi2**2*xi5**2)*F112255+                 &
     (xi1**3*xi3**3+ xi1**3*xi2**3)*F111222+ (xi1**3*xi5**2*xi6+           &
     xi1**3*xi6**2*xi5)*F111556+ (xi1*xi3**2*xi2*xi4**2+                   &
     xi1*xi2**2*xi3*xi4**2)*F122344+ (xi1**2*xi2**4+                       &
     xi1**2*xi3**4)*F113333+ (xi3*xi2*xi4**3*xi5+                          &
     xi2*xi3*xi4**3*xi6)*F234446+ (xi1**2*xi3*xi4**3+                      &
     xi1**2*xi2*xi4**3)*F112444+ F223356*xi2**2*xi3**2*xi5*xi6+            &
     F112356*xi1**2*xi2*xi3*xi5*xi6+ (xi1**2*xi3**2*xi2*xi4+               &
     xi1**2*xi2**2*xi3*xi4)*F112234+ (xi1**4*xi6**2+                       &
     xi1**4*xi5**2)*F111166+ (xi1**5*xi3+ xi1**5*xi2)*F111112
     s3 = s4+ F123456*xi1*xi2*xi3*xi4*xi5*xi6+                             &
     (xi1*xi2**3*xi4*xi5+ xi1*xi3**3*xi4*xi6)*F122245+                     &
     (xi1**3*xi6**3+ xi1**3*xi5**3)*F111666+ (xi1**2*xi2*xi6**3+           &
     xi1**2*xi3*xi5**3)*F112666+ (xi1**2*xi2**3*xi4+                       &
     xi1**2*xi3**3*xi4)*F113334+ (xi3**2*xi4**2*xi6**2+                    &
     xi2**2*xi4**2*xi5**2)*F224455+ (xi2**2*xi3**2*xi6**2+                 &
     xi3**2*xi2**2*xi5**2)*F223366+ (xi3**2*xi6**4+                        &
     xi2**2*xi5**4)*F225555+ F112344*xi1**2*xi2*xi3*xi4**2+                &
     (xi1**2*xi3**2*xi4*xi5+ xi1**2*xi2**2*xi4*xi6)*F112246+               &
     F111234*xi1**3*xi2*xi3*xi4+ F111456*xi1**3*xi4*xi5*xi6+               &
     (xi1**3*xi3*xi6**2+ xi1**3*xi2*xi5**2)*F111255+ (xi6**4*xi5**2+       &
     xi5**4*xi6**2)*F555566+ (xi1**3*xi2**2*xi6+                           &
     xi1**3*xi3**2*xi5)*F111335+ F235566*xi2*xi3*xi5**2*xi6**2
     s4 = s3+ (xi1**3*xi3*xi5**2+ xi1**3*xi2*xi6**2)*F111266+              &
     (xi1**4*xi3*xi4+ xi1**4*xi2*xi4)*F111124+ (xi3*xi2**2*xi5**3+         &
     xi2*xi3**2*xi6**3)*F233666+ (xi2*xi4**3*xi5*xi6+                      &
     xi3*xi4**3*xi6*xi5)*F244456+ (xi4**2*xi6*xi5**3+                      &
     xi4**2*xi5*xi6**3)*F445666+ (xi1*xi2*xi5**2*xi6**2+                   &
     xi1*xi3*xi6**2*xi5**2)*F125566+ (xi2**3*xi4**2*xi6+                   &
     xi3**3*xi4**2*xi5)*F222446+ F144456*xi1*xi4**3*xi5*xi6+               &
     (xi4**3*xi5*xi6**2+ xi4**3*xi6*xi5**2)*F444566+                       &
     (xi1**3*xi3*xi4*xi6+ xi1**3*xi2*xi4*xi5)*F111346+                     &
     F145566*xi1*xi4*xi5**2*xi6**2+ (xi2*xi3*xi5**4+                       &
     xi3*xi2*xi6**4)*F235555+ F123444*xi1*xi2*xi3*xi4**3+                  &
     (xi1*xi3**4*xi2+ xi1*xi2**4*xi3)*F122223+ (xi4**3*xi5**3+             &
     xi4**3*xi6**3)*F444666
     s2 = s4+ (xi3**2*xi4**2*xi5**2+                                       &
     xi2**2*xi4**2*xi6**2)*F224466+ F444444*xi4**6+                        &
     F122334*xi1*xi2**2*xi3**2*xi4+ (xi1**3*xi2**2*xi3+                    &
     xi1**3*xi3**2*xi2)*F111223+ (xi1*xi2**5+ xi1*xi3**5)*F122222+         &
     (xi3*xi4**3*xi5**2+ xi2*xi4**3*xi6**2)*F244466+ (xi2**4*xi4*xi6+      &
     xi3**4*xi4*xi5)*F222246+ F114456*xi1**2*xi4**2*xi5*xi6+               &
     (xi4*xi6**5+ xi4*xi5**5)*F466666+ (xi2**2*xi5**2*xi6**2+              &
     xi3**2*xi6**2*xi5**2)*F225566+ (xi2*xi3**4*xi6+                       &
     xi3*xi2**4*xi5)*F233336+ (xi6**5*xi5+ xi5**5*xi6)*F555556+            &
     (xi1*xi3*xi2*xi4*xi5**2+ xi1*xi2*xi3*xi4*xi6**2)*F123466+             &
     (xi6**6+ xi5**6)*F666666+ (xi2*xi4**2*xi5**2*xi6+                     &
     xi3*xi4**2*xi6**2*xi5)*F244556+ (xi1*xi2**2*xi4*xi6**2+               &
     xi1*xi3**2*xi4*xi5**2)*F122466
     s4 = s2+ (xi1*xi6**2*xi5**3+ xi1*xi5**2*xi6**3)*F155666+              &
     (xi2*xi4**2*xi6**3+ xi3*xi4**2*xi5**3)*F244666+                       &
     (xi1*xi3**2*xi4**2*xi6+ xi1*xi2**2*xi4**2*xi5)*F122445+               &
     (xi1*xi4*xi5**4+ xi1*xi4*xi6**4)*F146666+ (xi2**3*xi5**2*xi6+         &
     xi3**3*xi6**2*xi5)*F222556+ (xi3**2*xi2*xi4**3+                       &
     xi2**2*xi3*xi4**3)*F223444+ (xi2*xi5**4*xi6+                          &
     xi3*xi6**4*xi5)*F255556+ (xi1*xi3*xi4*xi5**2*xi6+                     &
     xi1*xi2*xi4*xi6**2*xi5)*F134556+ (xi1**3*xi3**2*xi4+                  &
     xi1**3*xi2**2*xi4)*F111334+ (xi3**2*xi6*xi5**3+                       &
     xi2**2*xi5*xi6**3)*F225666+ (xi4**5*xi6+ xi4**5*xi5)*F444446+         &
     (xi3*xi2**2*xi4*xi6**2+ xi2*xi3**2*xi4*xi5**2)*F233455+               &
     (xi3*xi6*xi5**4+ xi2*xi5*xi6**4)*F256666+ (xi1**3*xi4*xi5**2+         &
     xi1**3*xi4*xi6**2)*F111466
     s5 = s4+ (xi2**5*xi5+ xi3**5*xi6)*F222225+                            &
     (xi1**2*xi3*xi4*xi6*xi5+ xi1**2*xi2*xi4*xi5*xi6)*F112456+             &
     (xi1*xi6**5+ xi1*xi5**5)*F166666+ (xi3**2*xi2*xi6*xi5**2+             &
     xi2**2*xi3*xi5*xi6**2)*F223566+ (xi1*xi4**2*xi6*xi5**2+               &
     xi1*xi4**2*xi5*xi6**2)*F144566+ (xi2**2*xi4**2*xi5*xi6+               &
     xi3**2*xi4**2*xi6*xi5)*F224456+ (xi1*xi3**3*xi5**2+                   &
     xi1*xi2**3*xi6**2)*F122266
     s3 = s5+ (xi3**2*xi2*xi6**2*xi5+                                      &
     xi2**2*xi3*xi5**2*xi6)*F223556+ (xi1*xi3**3*xi6*xi5+                  &
     xi1*xi2**3*xi5*xi6)*F122256+ (xi1*xi3*xi4*xi6**3+                     &
     xi1*xi2*xi4*xi5**3)*F134666+ (xi1*xi2**3*xi5**2+                      &
     xi1*xi3**3*xi6**2)*F122255+ (xi2*xi5**3*xi6**2+                       &
     xi3*xi6**3*xi5**2)*F255566+ (xi4**2*xi5**4+                           &
     xi4**2*xi6**4)*F445555+ (xi3**5*xi5+ xi2**5*xi6)*F222226+             &
     (xi1*xi3*xi4**3*xi6+ xi1*xi2*xi4**3*xi5)*F124445+                     &
     (xi2**2*xi4*xi5*xi6**2+ xi3**2*xi4*xi6*xi5**2)*F224566

     s4 = s3+ (xi2*xi4**4*xi5+ xi3*xi4**4*xi6)*F244445+                    &
     (xi1**2*xi4**3*xi5+ xi1**2*xi4**3*xi6)*F114446+                       &
     (xi1**2*xi3**2*xi5**2+ xi1**2*xi2**2*xi6**2)*F112266+                 &
     (xi3*xi2*xi6**3*xi5+ xi2*xi3*xi5**3*xi6)*F235556+                     &
     (xi3**3*xi2*xi6**2+ xi2**3*xi3*xi5**2)*F222355+                       &
     (xi1**2*xi4**2*xi6**2+ xi1**2*xi4**2*xi5**2)*F114466+                 &
     (xi3**3*xi2*xi5**2+ xi2**3*xi3*xi6**2)*F222366+ (xi2**5*xi4+          &
     xi3**5*xi4)*F222224+ (xi2**3*xi4**2*xi5+                              &
     xi3**3*xi4**2*xi6)*F222445
     s4 = s4 + (xi2**4*xi4*xi5+                                            &
     xi3**4*xi4*xi6)*F222245+ (xi1*xi3**4*xi4+                             &
     xi1*xi2**4*xi4)*F122224+ (xi3**3*xi2**2*xi6+                          &
     xi2**3*xi3**2*xi5)*F222335+ (xi1*xi3*xi6*xi5**3+                      &
     xi1*xi2*xi5*xi6**3)*F125666+ (xi1*xi2*xi4*xi6**3+                     &
     xi1*xi3*xi4*xi5**3)*F134555+ (xi2**3*xi3*xi5*xi6+                     &
     xi3**3*xi2*xi6*xi5)*F222356

     s5 = s4+ (xi2**6+ xi3**6)*F222222+ (xi1*xi2**2*xi4*xi5*xi6+           &
     xi1*xi3**2*xi4*xi6*xi5)*F122456+ (xi4*xi6*xi5**4+                     &
     xi4*xi5*xi6**4)*F456666+ (xi1*xi3*xi6**4+                             &
     xi1*xi2*xi5**4)*F125555+ (xi1*xi3**3*xi2*xi5+                         &
     xi1*xi2**3*xi3*xi6)*F122236+ (xi1*xi2**3*xi3*xi5+                     &
     xi1*xi3**3*xi2*xi6)*F122235+ (xi2**2*xi4*xi5**3+                      &
     xi3**2*xi4*xi6**3)*F224555+ (xi1*xi3**3*xi4*xi5+                      &
     xi1*xi2**3*xi4*xi6)*F122246
     t6 = s5+ (xi3**2*xi4**3*xi5+ xi2**2*xi4**3*xi6)*F224446+              &
     (xi2**2*xi6**4+ xi3**2*xi5**4)*F226666+ (xi2**3*xi3**2*xi6+           &
     xi3**3*xi2**2*xi5)*F222336+ (xi1**2*xi3**2*xi2*xi6+                   &
     xi1**2*xi2**2*xi3*xi5)*F112235+ (xi2**2*xi5**3*xi6+                   &
     xi3**2*xi6**3*xi5)*F225556+ (xi1**2*xi2**2*xi5*xi6+                   &
     xi1**2*xi3**2*xi6*xi5)*F112256+ (xi1*xi3*xi2**2*xi5**2+               &
     xi1*xi2*xi3**2*xi6**2)*F123366+ (xi2**5*xi3+                          &
     xi3**5*xi2)*F222223+ (xi2*xi4**5+ xi3*xi4**5)*F244444
   endif


   f = (mu0+t1+t2+t3+t4+t5+t6)

 end function MLdms2loc_xy3



end module mol_user

!
!  This unit defines all specific routines for a triatomic molecule of XY2 type
!
module pot_xy2
  use accuracy
  use moltype
  !use pot_h216o

  implicit none

  public MLpoten_xy2_morbid
  public MLpoten_xy2_halonen_I,MLpoten_xy2_dmbe,MLdms2pqr_xy2,MLloc2pqr_xy2,MLpoten_xy2_tyuterev
  public MLpoten_h2o_tennyson,MLpoten_xy2_schwenke,MLpoten_c3_mladenovic
  public MLpoten_SO2_pes_8d,MLpoten_so2_damp,MLpoten_co2_ames1,MLpoten_so2_ames1,MLpoten_c3_R_theta
  public MLpoten_xy2_tyuterev_damp,MLdms2pqr_xy2_coeff,MLpoten_xy2_mlt_co2,MLpoten_h2s_dvr3d,MLdipole_so2_ames1,MLdipole_ames1,&
         MLdipole_xy2_lorenzo,MLdms2pqr_xy2_linear,MLpoten_xy2_bubukina,MLpoten_xyz_tyuterev,MLdms2pqr_xyz_coeff
  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains



  !
  ! Defining potential energy function 
  !
  ! This is a MORBID type PES for XY2 molecules, Taylor expansion... 
  ! with respect to the GD-coordinates
  ! V = ..........
  !
  function MLpoten_xy2_morbid(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)            :: r12,r32,alpha,rhoe,re13

   real(ark)  ::  ve    ,re12  ,aa1   ,f1    ,f11   ,f13,&
         f111  ,f113  ,f1111 ,f1113 ,f1133 ,fa1   ,     &
         fa2   ,fa3   ,fa4   ,fa5   ,fa6   ,fa7   ,     &
         fa8   ,f1a1  ,f2a1  ,f3a1  ,f4a1  ,f1a11 ,     &
         f2a11 ,f3a11 ,f1a13 ,f2a13 ,f3a13 ,f1a111,     &
         f2a111,f1a113,f2a113,f1a1111,f1a1113,f1a1133

   real(ark) :: fe1,fe3,fe11,fe33,fe13,fe111,fe333,          &
         fe113,fe133,fe1111,fe3333,fe1113,fe1333,fe1133 ,   &
         f1a3  ,f2a3  ,f3a3  ,f4a3  ,f1a33 ,f2a33 ,f3a33 ,  &
         f1a333,f2a333,f1a133,f2a133,f1a3333,f1a1333 ,      & 
         f3,f33  ,f333 ,f133 ,f3333,f1333      

   real(ark)  :: y1,y3,v,coro,v0,a1(3),a0(3),a2(3),t1,t2,w(3),cosalpha,xyz0(3,3)
   integer(ik) :: vtype
   character(len=cl)     :: txt
   
   if (verbose>=6) write(out,"('MLpoten_xy2_morbid/start')") 
     !
     vtype = 2
     !
     r12 = local(1) ; r32 = local(2) ;  alpha = local(3)
     !
     if (size(force)/=36) then  
       !
       write(out,"('MLpoten_xy2_morbid: The new form of this potential must include re,rhoe,a and 36 parameters')")
       stop 'MLpoten_xy2_morbid: Illegal number of potential parameteres'
       !
     endif 
     !
     re12 = force( 1)
     rhoe = force( 2)*pi/180.0_ark
     aa1  = force( 3)
     !
     !if  (molec%Ndihedrals>1) then
     !  !
     !  alpha = pi-asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
     !  !
     !endif 
     !
     !if  (molec%Ndihedrals>1) then
     !  !
     !  !alpha = pi-asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
     !  !
     !  !
     !  !
     !  alpha = pi-asin( sqrt( (local(3))**2+(local(4))**2 ))
     !  !
     !endif 
     !
     select case(trim(molec%coords_transform))
     case default
        !
        alpha = local(3)
        !
     case('R-PHI1-PHI2','R-PHI1-PHI2-Z')
        !
        a1(:) = xyz(2,:) - xyz(1,:)
        a2(:) = xyz(3,:) - xyz(1,:)
        !
        t1 =  sqrt(sum(a1(:)**2))
        t2 =  sqrt(sum(a2(:)**2))
        !
        a1 =  a1(:)/t1
        a2 =  a2(:)/t2
        !
        w(:) = MLvector_product(a1,a2)
        !
        cosalpha = sum(a1(:)*a2(:))
        !
        txt = "pot_xy2"
        !
        alpha = aacos(cosalpha,txt)
        !
        if (alpha<sqrt(small_)) alpha = pi-asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
        !
        !alpha = pi-asin( sqrt( (local(3))**2+(local(4))**2 ))
        !
     case('R1-R2-Y+X')
        !
        call MLfromlocal2cartesian(1_ik,local,xyz0)
        !
        a1(:) = xyz0(3,:) - xyz0(1,:)
        a2(:) = xyz0(2,:) - xyz0(1,:)
        !
        t1 =  sqrt(sum(a1(:)**2))
        t2 =  sqrt(sum(a2(:)**2))
        !
        a1 =  a1(:)/t1
        a2 =  a2(:)/t2
        !
        w(:) = MLvector_product(a1,a2)
        !
        cosalpha = sum(a1(:)*a2(:))
        !
        txt = "pot_xy2"
        !
        alpha = aacos(cosalpha,txt)
        !
        !if (alpha<sqrt(small_)) alpha = pi-asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
        !
     end select 
     !
     !if (molec%Nangles==0) then
     !  if (molec%Ndihedrals>1) then
     !     alpha = pi-asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
     !  else
     !     alpha = pi-asin( sqrt( sin(local(3))**2))
     !  endif 
     !endif 
     !
     if (rhoe==0.0_ark) vtype = 1
     !
     ve        = force( 4)
     fa1       = force( 5)
     fa2       = force( 6)
     fa3       = force( 7)
     fa4       = force( 8)
     fa5       = force( 9)
     fa6       = force(10)
     fa7       = force(11)
     fa8       = force(12)
     f1a1      = force(13)
     f2a1      = force(14)
     f3a1      = force(15)
     f4a1      = force(16)
     f11       = force(17)
     f1a11     = force(18)
     f2a11     = force(19)
     f3a11     = force(20)
     f13       = force(21)
     f1a13     = force(22)
     f2a13     = force(23)
     f3a13     = force(24)
     f111      = force(25)
     f1a111    = force(26)
     f2a111    = force(27)
     f113      = force(28)
     f1a113    = force(29)
     f2a113    = force(30)
     f1111     = force(31)
     f1a1111   = force(32)
     f1113     = force(33)
     f1a1113   = force(34)
     f1133     = force(35)
     f1a1133   = force(36)

     f1 = 0.0_ark
     f3  = f1
     f33  = f11
     f333  = f111
     f133  = f113
     f3333  = f1111
     f1333  = f1113


     f1a3    = f1a1
     f2a3    = f2a1
     f3a3    = f3a1
     f4a3    = f4a1
     f1a33   = f1a11 
     f2a33   = f2a11
     f3a33   = f3a11
     f1a333  = f1a111
     f2a333  = f2a111
     f1a133  = f1a113
     f2a133  = f2a113
     f1a3333 = f1a1111
     f1a1333 = f1a1113


!
! calculate potential energy function values
!
     if (vtype.eq.1) then
       coro=1.0_ark+cos(alpha)
     else if (vtype==2) then 
       coro=cos(rhoe)+cos(alpha)
     else
       stop 'vtype = 1 or 2 only!'
     endif 
     !
     !coro=cos(rhoe)+cos(alpha)
     !
     !coro= pi-alpha-rhoe
     !
     !
     !if(trim(molec%coords_transform)=='R12-RHO') then
     !  !
     !  re12 = ML_MEP_xy2_R12_ALPHA(alpha)
     !  !
     !elseif(trim(molec%coords_transform)=='R13-RHO') then
     !  !
     !  re13 = ML_MEP_xy2_R13_ALPHA(alpha)
     !  re12 = 0.5_ark*re13/sin(alpha*0.5_ark)
     !  !
     !endif
     !
     y1=1.0_ark-exp(-aa1*(r12-re12))
     y3=1.0_ark-exp(-aa1*(r32-re12))

!
! calculate potential energy function values
!

     v0= ve+fa1*coro+fa2*coro**2+fa3*coro**3+fa4*coro**4+fa5*coro**5+fa6*coro**6+fa7*coro**7+fa8*coro**8
     fe1= f1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4
     fe3= f3+f1a3*coro+f2a3*coro**2+f3a3*coro**3+f4a3*coro**4
     fe11= f11+f1a11*coro+f2a11*coro**2+f3a11*coro**3
     fe33= f33+f1a33*coro+f2a33*coro**2+f3a33*coro**3
     fe13= f13+f1a13*coro+f2a13*coro**2+f3a13*coro**3
     fe111= f111+f1a111*coro+f2a111*coro**2
     fe333= f333+f1a333*coro+f2a333*coro**2
     fe113= f113+f1a113*coro+f2a113*coro**2
     fe133= f133+f1a133*coro+f2a133*coro**2
     fe1111= f1111+f1a1111*coro
     fe3333= f3333+f1a3333*coro
     fe1113= f1113+f1a1113*coro
     fe1333= f1333+f1a1333*coro
     fe1133= f1133+f1a1133*coro
     v     =  v0+fe1*y1+fe3*y3                            &
              +fe11*y1**2+fe33*y3**2+fe13*y1*y3           &
              +fe111*y1**3+fe333*y3**3+fe113*y1**2*y3     &
              +fe133*y1*y3**2                             &
              +fe1111*y1**4+fe3333*y3**4+fe1113*y1**3*y3  &
              +fe1333*y1*y3**3+fe1133*y1**2*y3**2

     f = v

     !f = alpha

     if (verbose>=6) write(out,"('MLpoten_xy2_morbid/end')") 
 
 end function MLpoten_xy2_morbid



 function MLpoten_xy2_halonen_I(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)    :: r1,r2,da,T2,T1,y1,y2,alpha
   real(ark)    :: r0,alpha0

   real(ark)    :: v1,v0,v5,v4,v2,v3,ar,de, frr1,ftt,fttt,ftttt,frrr1,frt,frrt,frtt,frr1t,frrtt

      if (verbose>=6) write(out,"('MLpoten_xy3_halonen_I/start')") 


      r0      =  molec%req(1)
      alpha0  =  molec%alphaeq(1)
      ar      = molec%specparam(1)

      de      =  force(  1)
      T1      =  force(  2)
      T2      =  force(  3)

      ftt     =  force(  4)
      fttt    =  force(  5)
      ftttt   =  force(  6)


      frr1    =  force(  7)
      frrr1   =  force(  8)
      frt     =  force(  9)
      frtt    =  force( 10)
      frrt    =  force( 11)
      frr1t   =  force( 12)
      frrtt   =  force( 13)
      !
      r1 = local(1) ;  r2  = local(2) ;  alpha = local(3)
      !
      r1 = r1-r0
      r2 = r2-r0
      !
      da = alpha-alpha0
      !
      y1=1.0_ark-exp(-ar*(r1))
      y2=1.0_ark-exp(-ar*(r2))

      !
      !------------------------------------------------------
      !
      v0 = De*(y1**2+y2**2)+frr1/ar**2*(y1*y2) + &
           T1*(y1**3+y2**3)+&
           T2*(y1**4+y2**4)+&
           frt*(y1+y2)*da/ar+&
           0.5_ark*(frrr1/ar**3+frr1/ar**2)*(y1**2*y2+y2**2*y1)

      v1 = 0.5_ark*ftt*da**2+&
           fttt*da**3/6.0_ark+&
           ftttt*da**4/24.0_ark

      v2 = frt/ar*(y1+y2)*da

      v3 = 0.5_ark*frtt/ar*(y1+y2)*da**2

      v4 = 0.25_ark*(frrtt/ar**2+frtt/ar)*(y1**2+y2**2)*da**2


      v5 = 0.5_ark*(frrt/ar**2+frt/ar)*(y1**2+y2**2)*da+frr1t/ar**2*y1*y2*da

      !
      f =  (v0+v1+v2+v3+v4+v5)*1.0e-11/planck/vellgt !  50340.359783704819122_ark


      if (verbose>=6) write(out,"('MLpoten_xy3_halonen_I/end')") 
 
  end function MLpoten_xy2_halonen_I



 function MLpoten_xy2_mlt_co2(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)    :: r1,r2,da,y1,y2,rho

   real(ark)    :: v0,r0,rho0,frr,frr1,faa,frrr,frrr1,faar,frrrr,frrrr1,frrr1r1,faarr,faarr1,faaaa  

      if (verbose>=6) write(out,"('MLpoten_xy2_mlt_co2/start')") 


      !
      r1 = local(1) ;  r2  = local(2) ;  rho = pi-local(3)
      !
      if  (molec%Ndihedrals>1) then
        !
        rho = asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
        !
      endif 
      !
      !r0      =  molec%req(1)
      !alpha0  =  molec%alphaeq(1)
      !
      v0      =  force(  1)
      r0      =  force(  2)
      rho0    =  force(  3)
      frr     =  force(  4)
      frr1    =  force(  5)
      faa     =  force(  6)
      frrr    =  force(  7)
      frrr1   =  force(  8)
      faar    =  force(  9)
      frrrr   =  force( 10)
      frrrr1  =  force( 11)
      frrr1r1 =  force( 12)
      faarr   =  force( 13)
      faarr1  =  force( 14)
      faaaa   =  force( 15)
      !
      y1 = r1-r0
      y2 = r2-r0
      !
      da = rho-rho0
      !
      f  = v0+&
           frr*(y1**2+y2**2)/2.0_ark+&
           frr1*(y1*y2) + &
           faa*da**2/2.0_ark+&
           frrr*(y1**3+y2**3)/6.0_ark+&
           frrr1*(y1**2*y2+y2**2*y1)/2.0_ark+&
           faar*(y1+y2)*da**2/2.0_ark+&
           frrrr*(y1**4+y2**4)/24.0_ark+&
           frrrr1*(y1**3*y2+y2**3*y1)/6.0_ark+&
           frrr1r1*(y1**2*y2**2)/4.0_ark+&
           faarr*(y1**2+y2**2)*da**2/4.0_ark+&
           faarr1*(y1*y2)*da**2/2.0_ark+&
           faaaa*da**4/24.0_ark

      !
      f =  f*1.0e-11/planck/vellgt !  50340.359783704819122_ark


      if (verbose>=6) write(out,"('MLpoten_xy2_mlt_co2/end')") 
 
  end function MLpoten_xy2_mlt_co2



  !
  ! Defining potential energy function
  !
  function MLpoten_xy2_tyuterev(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)            :: r12,r32,alpha,xcos,v,v1,v2,v3,v4,v5,v6,v7,v8

   real(ark)            :: aa1,re12,alphae,xst,y1,y2,y3,xs1,xs2,v0,vp1,vp2,vp3
   real(ark)            :: g1,g2,b1,b2,rhh,vhh
   integer(ik)          :: N
   real(ark)            :: ycos,v_t,q1,q2
   real(ark)            :: a1(3),a2(3),t1,t2,w(3),cosalpha
   character(len=cl)    :: txt
   !
   if (verbose>=6) write(out,"('MLpoten_xy2_tyuterev/start')")
   !
   r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
   !
   re12    = molec%req(1)
   aa1     = molec%specparam(1)
   !
   if (molec%Nangles>0) then
     alphae  = molec%alphaeq(1)
   else
     alphae = pi+(-molec%taueq(1)+molec%taueq(2))
     alpha = pi - local(3) 
   endif 
   !
   select case(trim(molec%coords_transform))
   case default
      !
      alpha = local(3)
      !
   case('R-PHI1-PHI2','R-PHI1-PHI2-Z')
      !
      a1(:) = xyz(2,:) - xyz(1,:)
      a2(:) = xyz(3,:) - xyz(1,:)
      !
      t1 =  sqrt(sum(a1(:)**2))
      t2 =  sqrt(sum(a2(:)**2))
      !
      a1 =  a1(:)/t1
      a2 =  a2(:)/t2
      !
      w(:) = MLvector_product(a1,a2)
      !
      cosalpha = sum(a1(:)*a2(:))
      !
      txt = "pot_xy2"
      !
      alpha = aacos(cosalpha,txt)
      !
      if (alpha<sqrt(small_)) alpha = pi-asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
      !
      !alpha = pi-asin( sqrt( (local(3))**2+(local(4))**2 ))
      !
   end select 
   !
   b1   = force(1)
   b2   = force(2)
   g1   = force(3)
   g2   = force(4)
   !
   rhh=sqrt(r12**2+r32**2-2._ark*r12*r32*cos(alpha))
   vhh=b1*exp(-g1*rhh)+b2*exp(-g2*rhh**2)
   !
   ! calculate potential energy function values
   !
   y1=1.0_ark-exp(-aa1*(r12-re12))
   y2=1.0_ark-exp(-aa1*(r32-re12))
   !
   y3=cos(alpha)-cos(alphae)
   !
   N = size(force)
   !
   v4 = 0 ; v5 = 0 ; v6 = 0 ; v7 = 0 ; v8 = 0
   !
 v0 = force(5)*y1**0*y2**0*y3**0
 v1 = force(6)*y1**0*y2**0*y3**1&
    + force(7)*y1**1*y2**0*y3**0&
    + force(7)*y1**0*y2**1*y3**0
 v2 = force(8)*y1**0*y2**0*y3**2&
    + force(9)*y1**1*y2**0*y3**1&
    + force(9)*y1**0*y2**1*y3**1&
    + force(10)*y1**1*y2**1*y3**0&
    + force(11)*y1**2*y2**0*y3**0&
    + force(11)*y1**0*y2**2*y3**0
 v3 = force(12)*y1**0*y2**0*y3**3&
    + force(13)*y1**1*y2**0*y3**2&
    + force(13)*y1**0*y2**1*y3**2&
    + force(14)*y1**1*y2**1*y3**1&
    + force(15)*y1**2*y2**0*y3**1&
    + force(15)*y1**0*y2**2*y3**1&
    + force(16)*y1**2*y2**1*y3**0&
    + force(16)*y1**1*y2**2*y3**0&
    + force(17)*y1**3*y2**0*y3**0&
    + force(17)*y1**0*y2**3*y3**0

 if (N>18) then
  v4 = force(18)*y1**0*y2**0*y3**4&
    + force(19)*y1**1*y2**0*y3**3&
    + force(19)*y1**0*y2**1*y3**3&
    + force(20)*y1**1*y2**1*y3**2&
    + force(21)*y1**2*y2**0*y3**2&
    + force(21)*y1**0*y2**2*y3**2&
    + force(22)*y1**2*y2**1*y3**1&
    + force(22)*y1**1*y2**2*y3**1&
    + force(23)*y1**2*y2**2*y3**0&
    + force(24)*y1**3*y2**0*y3**1&
    + force(24)*y1**0*y2**3*y3**1&
    + force(25)*y1**3*y2**1*y3**0&
    + force(25)*y1**1*y2**3*y3**0&
    + force(26)*y1**4*y2**0*y3**0&
    + force(26)*y1**0*y2**4*y3**0
endif

 if (N>26) then
  v5 = force(27)*y1**0*y2**0*y3**5&
    + force(28)*y1**1*y2**0*y3**4&
    + force(28)*y1**0*y2**1*y3**4&
    + force(29)*y1**1*y2**1*y3**3&
    + force(30)*y1**2*y2**0*y3**3&
    + force(30)*y1**0*y2**2*y3**3&
    + force(31)*y1**2*y2**1*y3**2&
    + force(31)*y1**1*y2**2*y3**2&
    + force(32)*y1**2*y2**2*y3**1&
    + force(33)*y1**3*y2**0*y3**2&
    + force(33)*y1**0*y2**3*y3**2&
    + force(34)*y1**3*y2**1*y3**1&
    + force(34)*y1**1*y2**3*y3**1&
    + force(35)*y1**3*y2**2*y3**0&
    + force(35)*y1**2*y2**3*y3**0&
    + force(36)*y1**4*y2**0*y3**1&
    + force(36)*y1**0*y2**4*y3**1&
    + force(37)*y1**4*y2**1*y3**0&
    + force(37)*y1**1*y2**4*y3**0&
    + force(38)*y1**5*y2**0*y3**0&
    + force(38)*y1**0*y2**5*y3**0
endif

 if (N>38) then
  v6 = force(39)*y1**0*y2**0*y3**6&
    + force(40)*y1**1*y2**0*y3**5&
    + force(40)*y1**0*y2**1*y3**5&
    + force(41)*y1**1*y2**1*y3**4&
    + force(42)*y1**2*y2**0*y3**4&
    + force(42)*y1**0*y2**2*y3**4&
    + force(43)*y1**2*y2**1*y3**3&
    + force(43)*y1**1*y2**2*y3**3&
    + force(44)*y1**2*y2**2*y3**2&
    + force(45)*y1**3*y2**0*y3**3&
    + force(45)*y1**0*y2**3*y3**3&
    + force(46)*y1**3*y2**1*y3**2&
    + force(46)*y1**1*y2**3*y3**2&
    + force(47)*y1**3*y2**2*y3**1&
    + force(47)*y1**2*y2**3*y3**1&
    + force(48)*y1**3*y2**3*y3**0&
    + force(49)*y1**4*y2**0*y3**2&
    + force(49)*y1**0*y2**4*y3**2&
    + force(50)*y1**4*y2**1*y3**1&
    + force(50)*y1**1*y2**4*y3**1&
    + force(51)*y1**4*y2**2*y3**0&
    + force(51)*y1**2*y2**4*y3**0&
    + force(52)*y1**5*y2**0*y3**1&
    + force(52)*y1**0*y2**5*y3**1&
    + force(53)*y1**5*y2**1*y3**0&
    + force(53)*y1**1*y2**5*y3**0&
    + force(54)*y1**6*y2**0*y3**0&
    + force(54)*y1**0*y2**6*y3**0
 endif

 if (N>54) then
 v7 = force(55)*y1**0*y2**0*y3**7&
    + force(56)*y1**1*y2**0*y3**6&
    + force(56)*y1**0*y2**1*y3**6&
    + force(57)*y1**1*y2**1*y3**5&
    + force(58)*y1**2*y2**0*y3**5&
    + force(58)*y1**0*y2**2*y3**5&
    + force(59)*y1**2*y2**1*y3**4&
    + force(59)*y1**1*y2**2*y3**4&
    + force(60)*y1**2*y2**2*y3**3&
    + force(61)*y1**3*y2**0*y3**4&
    + force(61)*y1**0*y2**3*y3**4&
    + force(62)*y1**3*y2**1*y3**3&
    + force(62)*y1**1*y2**3*y3**3&
    + force(63)*y1**3*y2**2*y3**2&
    + force(63)*y1**2*y2**3*y3**2&
    + force(64)*y1**3*y2**3*y3**1&
    + force(65)*y1**4*y2**0*y3**3&
    + force(65)*y1**0*y2**4*y3**3&
    + force(66)*y1**4*y2**1*y3**2&
    + force(66)*y1**1*y2**4*y3**2&
    + force(67)*y1**4*y2**2*y3**1&
    + force(67)*y1**2*y2**4*y3**1&
    + force(68)*y1**4*y2**3*y3**0&
    + force(68)*y1**3*y2**4*y3**0&
    + force(69)*y1**5*y2**0*y3**2&
    + force(69)*y1**0*y2**5*y3**2&
    + force(70)*y1**5*y2**1*y3**1&
    + force(70)*y1**1*y2**5*y3**1&
    + force(71)*y1**5*y2**2*y3**0&
    + force(71)*y1**2*y2**5*y3**0&
    + force(72)*y1**6*y2**0*y3**1&
    + force(72)*y1**0*y2**6*y3**1&
    + force(73)*y1**6*y2**1*y3**0&
    + force(73)*y1**1*y2**6*y3**0&
    + force(74)*y1**7*y2**0*y3**0&
    + force(74)*y1**0*y2**7*y3**0
 endif

 if (N>74) then
 v8 = force(75)*y1**0*y2**0*y3**8&
    + force(76)*y1**1*y2**0*y3**7&
    + force(76)*y1**0*y2**1*y3**7&
    + force(77)*y1**1*y2**1*y3**6&
    + force(78)*y1**2*y2**0*y3**6&
    + force(78)*y1**0*y2**2*y3**6&
    + force(79)*y1**2*y2**1*y3**5&
    + force(79)*y1**1*y2**2*y3**5&
    + force(80)*y1**2*y2**2*y3**4&
    + force(81)*y1**3*y2**0*y3**5&
    + force(81)*y1**0*y2**3*y3**5&
    + force(82)*y1**3*y2**1*y3**4&
    + force(82)*y1**1*y2**3*y3**4&
    + force(83)*y1**3*y2**2*y3**3&
    + force(83)*y1**2*y2**3*y3**3&
    + force(84)*y1**3*y2**3*y3**2&
    + force(85)*y1**4*y2**0*y3**4&
    + force(85)*y1**0*y2**4*y3**4&
    + force(86)*y1**4*y2**1*y3**3&
    + force(86)*y1**1*y2**4*y3**3&
    + force(87)*y1**4*y2**2*y3**2&
    + force(87)*y1**2*y2**4*y3**2&
    + force(88)*y1**4*y2**3*y3**1&
    + force(88)*y1**3*y2**4*y3**1&
    + force(89)*y1**4*y2**4*y3**0&
    + force(90)*y1**5*y2**0*y3**3&
    + force(90)*y1**0*y2**5*y3**3&
    + force(91)*y1**5*y2**1*y3**2&
    + force(91)*y1**1*y2**5*y3**2&
    + force(92)*y1**5*y2**2*y3**1&
    + force(92)*y1**2*y2**5*y3**1&
    + force(93)*y1**5*y2**3*y3**0&
    + force(93)*y1**3*y2**5*y3**0&
    + force(94)*y1**6*y2**0*y3**2&
    + force(94)*y1**0*y2**6*y3**2&
    + force(95)*y1**6*y2**1*y3**1&
    + force(95)*y1**1*y2**6*y3**1&
    + force(96)*y1**6*y2**2*y3**0&
    + force(96)*y1**2*y2**6*y3**0&
    + force(97)*y1**7*y2**0*y3**1&
    + force(97)*y1**0*y2**7*y3**1&
    + force(98)*y1**7*y2**1*y3**0&
    + force(98)*y1**1*y2**7*y3**0&
    + force(99)*y1**8*y2**0*y3**0&
    + force(99)*y1**0*y2**8*y3**0
endif
    !
    !th1 = 0.5d0*( 1.0d0-tanh( 0.0001_ark*( v0+v1+v2-50000_ark ) ) )
    !
    !f=(v0+v1+v2)+(v3+v4+v5+v6+v7+v8)*th1+vhh
    !
    !vhh = (v0+v1+v2)*(1.0_ark/r12**12+1.0_ark/r32**12+1.0_ark/rhh**12)*b1
    !
    f=v0+v1+v2+v3+v4+v5+v6+v7+v8+vhh
    !
    q1 = local(1)/bohr ; q2 = local(2)/bohr ;  ycos = cos(alpha)
    !
    v_t = 0
    !
    !call potv(v_t,q1,q2,ycos)
    !
    f = f + v_t*219474.630670_ark
    !
    if (verbose>=6) write(out,"('MLpoten_xy2_tyuterev/end')")

 end function MLpoten_xy2_tyuterev




function  MLpoten_h2o_tennyson(ncoords,natoms,local,xyz,force) result(v) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  v
   !
      integer(ik),parameter           ::  npropin = 240
      real(ark)            ::  r1,r2,th
      real(ark)            ::  xp(1:npropin)
      real(ark)            :: reoh,thetae,b1,roh,alphaoh
      real(ark)            :: phh2,t0,ut,x0,ut2,x02,xs1,xs2,xst,rs,rm,rr1,rr2,xep1,xep2,xep3
      real(ark)            :: rhh,vhh,vpb1,vpb2,v0,vp1,vp2,vp3,vp,vps1,vps2,y1,y2,y12,y22,voh1,voh2
      integer(ik) :: i
      !
      !
      xp(1:npropin) = force(1:npropin)
      !
      r1 = local(1) !/bohr
      r2 = local(2) !/bohr
      th = local(3)
      !
      reoh  = 0.958649_ark
      thetae= 104.3475_ark
      b1=2.0_ark
      roh=0.951961_ark
      alphaoh=2.587949757553683_ark
      phh2=6.70164303995_ark
      t0=0.01_ark
      ut=20._ark
      x0=2.5_ark
      ut2=20._ark
      x02=2.5_ark
        
      thetae=thetae*3.14159265358979312*.00555555555555555555_ark

      xs1=(r1+r2)*0.5_ark-reoh
      xs2=(r1-r2)*0.5_ark
      xst=cos(th)-cos(thetae)

      rs=sqrt(0.5_ark)*(r1+r2)
      rm=sqrt(0.5_ark)*(r1-r2)

      rr1=r1-roh
      rr2=r2-roh

      xep1=exp(-2._ark*alphaoh*rr1)-2._ark*exp(-alphaoh*rr1)+1._ark
      xep2=exp(-2._ark*alphaoh*rr2)-2._ark*exp(-alphaoh*rr2)+1._ark
      xep3=exp(-b1*(rr1**2+rr2**2))
      rhh=sqrt(r1**2+r2**2-2._ark*r1*r2*cos(th))
      vhh=0.900642240911285975d6*exp(-phh2*rhh)
      vpb1=0.518622556959170834d5*xep1
      vpb2=0.518622556959170834d5*xep2


  v0= xp(1)  *xs1**0*xs2**0*xst**0
 vp1= xp(2)  *xs1**0*xs2**0*xst**1&
     +xp(3)  *xs1**1*xs2**0*xst**0&
     +xp(4)  *xs1**0*xs2**0*xst**2&
     +xp(5)  *xs1**0*xs2**2*xst**0&
     +xp(6)  *xs1**1*xs2**0*xst**1&
     +xp(7)  *xs1**2*xs2**0*xst**0&
     +xp(8)  *xs1**0*xs2**0*xst**3&
     +xp(9)  *xs1**0*xs2**2*xst**1&
     +xp(10) *xs1**1*xs2**0*xst**2&
     +xp(11) *xs1**1*xs2**2*xst**0&
     +xp(12) *xs1**2*xs2**0*xst**1&
     +xp(13) *xs1**3*xs2**0*xst**0&
     +xp(14) *xs1**0*xs2**0*xst**4&
     +xp(15) *xs1**0*xs2**2*xst**2&
     +xp(16) *xs1**0*xs2**4*xst**0&
     +xp(17) *xs1**1*xs2**0*xst**3&
     +xp(18) *xs1**1*xs2**2*xst**1&
     +xp(19) *xs1**2*xs2**0*xst**2&
     +xp(20) *xs1**2*xs2**2*xst**0&
     +xp(21) *xs1**3*xs2**0*xst**1&
     +xp(22) *xs1**4*xs2**0*xst**0&
     +xp(23) *xs1**0*xs2**0*xst**5&
     +xp(24) *xs1**0*xs2**2*xst**3&
     +xp(25) *xs1**0*xs2**4*xst**1&
     +xp(26) *xs1**1*xs2**0*xst**4&
     +xp(27) *xs1**1*xs2**2*xst**2&
     +xp(28) *xs1**1*xs2**4*xst**0&
     +xp(29) *xs1**2*xs2**0*xst**3&
     +xp(30) *xs1**2*xs2**2*xst**1&
     +xp(31) *xs1**3*xs2**0*xst**2&
     +xp(32) *xs1**3*xs2**2*xst**0&
     +xp(33) *xs1**4*xs2**0*xst**1&
     +xp(34) *xs1**5*xs2**0*xst**0&
     +xp(35) *xs1**0*xs2**0*xst**6&
     +xp(36) *xs1**0*xs2**2*xst**4&
     +xp(37) *xs1**0*xs2**4*xst**2&
     +xp(38) *xs1**0*xs2**6*xst**0&
     +xp(39) *xs1**1*xs2**0*xst**5&
     +xp(40) *xs1**1*xs2**2*xst**3&
     +xp(41) *xs1**1*xs2**4*xst**1&
     +xp(42) *xs1**2*xs2**0*xst**4&
     +xp(43) *xs1**2*xs2**2*xst**2&
     +xp(44) *xs1**2*xs2**4*xst**0&
     +xp(45) *xs1**3*xs2**0*xst**3&
     +xp(46) *xs1**3*xs2**2*xst**1&
     +xp(47) *xs1**4*xs2**0*xst**2&
     +xp(48) *xs1**4*xs2**2*xst**0&
     +xp(49) *xs1**5*xs2**0*xst**1&
     +xp(50) *xs1**6*xs2**0*xst**0&
     +xp(51) *xs1**0*xs2**0*xst**7&
     +xp(52) *xs1**0*xs2**2*xst**5&
     +xp(53) *xs1**0*xs2**4*xst**3&
     +xp(54) *xs1**0*xs2**6*xst**1&
     +xp(55) *xs1**1*xs2**0*xst**6&
     +xp(56) *xs1**1*xs2**2*xst**4&
     +xp(57) *xs1**1*xs2**4*xst**2&
     +xp(58) *xs1**1*xs2**6*xst**0&
     +xp(59) *xs1**2*xs2**0*xst**5&
     +xp(60) *xs1**2*xs2**2*xst**3&
     +xp(61) *xs1**2*xs2**4*xst**1&
     +xp(62) *xs1**3*xs2**0*xst**4&
     +xp(63) *xs1**3*xs2**2*xst**2&
     +xp(64) *xs1**3*xs2**4*xst**0&
     +xp(65) *xs1**4*xs2**0*xst**3&
     +xp(66) *xs1**4*xs2**2*xst**1&
     +xp(67) *xs1**5*xs2**0*xst**2&
     +xp(68) *xs1**5*xs2**2*xst**0&
     +xp(69) *xs1**6*xs2**0*xst**1&
     +xp(70) *xs1**7*xs2**0*xst**0&
     +xp(71) *xs1**0*xs2**0*xst**8&
     +xp(72) *xs1**0*xs2**2*xst**6&
     +xp(73) *xs1**0*xs2**4*xst**4&
     +xp(74) *xs1**0*xs2**6*xst**2&
     +xp(75) *xs1**0*xs2**8*xst**0&
     +xp(76) *xs1**1*xs2**0*xst**7&
     +xp(77) *xs1**1*xs2**2*xst**5&
     +xp(78) *xs1**1*xs2**4*xst**3&
     +xp(79) *xs1**1*xs2**6*xst**1&
     +xp(80) *xs1**2*xs2**0*xst**6&
     +xp(81) *xs1**2*xs2**2*xst**4&
     +xp(82) *xs1**2*xs2**4*xst**2&
     +xp(83) *xs1**2*xs2**6*xst**0&
     +xp(84) *xs1**3*xs2**0*xst**5&
     +xp(85) *xs1**3*xs2**2*xst**3&
     +xp(86) *xs1**3*xs2**4*xst**1&
     +xp(87) *xs1**4*xs2**0*xst**4&
     +xp(88) *xs1**4*xs2**2*xst**2&
     +xp(89) *xs1**4*xs2**4*xst**0&
     +xp(90) *xs1**5*xs2**0*xst**3&
     +xp(91) *xs1**5*xs2**2*xst**1&
     +xp(92) *xs1**6*xs2**0*xst**2
 vp2= xp(93) *xs1**6*xs2**2*xst**0&
     +xp(94) *xs1**7*xs2**0*xst**1&
     +xp(95) *xs1**8*xs2**0*xst**0&
     +xp(96) *xs1**0*xs2**0*xst**9&
     +xp(97) *xs1**0*xs2**2*xst**7&
     +xp(98) *xs1**0*xs2**4*xst**5&
     +xp(99) *xs1**0*xs2**6*xst**3&
     +xp(100)*xs1**0*xs2**8*xst**1&
     +xp(101)*xs1**1*xs2**0*xst**8&
     +xp(102)*xs1**1*xs2**2*xst**6&
     +xp(103)*xs1**1*xs2**4*xst**4&
     +xp(104)*xs1**1*xs2**6*xst**2&
     +xp(105)*xs1**1*xs2**8*xst**0&
     +xp(106)*xs1**2*xs2**0*xst**7&
     +xp(107)*xs1**2*xs2**2*xst**5&
     +xp(108)*xs1**2*xs2**4*xst**3&
     +xp(109)*xs1**2*xs2**6*xst**1&
     +xp(110)*xs1**3*xs2**0*xst**6&
     +xp(111)*xs1**3*xs2**2*xst**4&
     +xp(112)*xs1**3*xs2**4*xst**2&
     +xp(113)*xs1**3*xs2**6*xst**0&
     +xp(114)*xs1**4*xs2**0*xst**5&
     +xp(115)*xs1**4*xs2**2*xst**3&
     +xp(116)*xs1**4*xs2**4*xst**1&
     +xp(117)*xs1**5*xs2**0*xst**4&
     +xp(118)*xs1**5*xs2**2*xst**2&
     +xp(119)*xs1**5*xs2**4*xst**0&
     +xp(120)*xs1**6*xs2**0*xst**3&
     +xp(121)*xs1**6*xs2**2*xst**1&
     +xp(122)*xs1**7*xs2**0*xst**2&
     +xp(123)*xs1**7*xs2**2*xst**0&
     +xp(124)*xs1**8*xs2**0*xst**1&
     +xp(125)*xs1**9*xs2**0*xst**0&
     +xp(126)*xs1**0*xs2**0*xst**10&
     +xp(127)*xs1**0*xs2**2*xst**8&
     +xp(128)*xs1**0*xs2**4*xst**6&
     +xp(129)*xs1**0*xs2**6*xst**4&
     +xp(130)*xs1**0*xs2**8*xst**2&
     +xp(131)*xs1**0*xs2**10*xst**0&
     +xp(132)*xs1**1*xs2**0*xst**9&
     +xp(133)*xs1**1*xs2**2*xst**7&
     +xp(134)*xs1**1*xs2**4*xst**5&
     +xp(135)*xs1**1*xs2**6*xst**3&
     +xp(136)*xs1**1*xs2**8*xst**1&
     +xp(137)*xs1**2*xs2**0*xst**8&
     +xp(138)*xs1**2*xs2**2*xst**6&
     +xp(139)*xs1**2*xs2**4*xst**4&
     +xp(140)*xs1**2*xs2**6*xst**2&
     +xp(141)*xs1**2*xs2**8*xst**0&
     +xp(142)*xs1**3*xs2**0*xst**7&
     +xp(143)*xs1**3*xs2**2*xst**5&
     +xp(144)*xs1**3*xs2**4*xst**3&
     +xp(145)*xs1**3*xs2**6*xst**1&
     +xp(146)*xs1**4*xs2**0*xst**6&
     +xp(147)*xs1**4*xs2**2*xst**4&
     +xp(148)*xs1**4*xs2**4*xst**2&
     +xp(149)*xs1**4*xs2**6*xst**0&
     +xp(150)*xs1**5*xs2**0*xst**5&
     +xp(151)*xs1**5*xs2**2*xst**3&
     +xp(152)*xs1**5*xs2**4*xst**1&
     +xp(153)*xs1**6*xs2**0*xst**4&
     +xp(154)*xs1**6*xs2**2*xst**2&
     +xp(155)*xs1**6*xs2**4*xst**0&
     +xp(156)*xs1**7*xs2**0*xst**3&
     +xp(157)*xs1**7*xs2**2*xst**1&
     +xp(158)*xs1**8*xs2**0*xst**2&
     +xp(159)*xs1**8*xs2**2*xst**0&
     +xp(160)*xs1**9*xs2**0*xst**1&
     +xp(161)*xs1**10*xs2**0*xst**0&
     +xp(162)*xs1**0*xs2**0*xst**11&
     +xp(163)*xs1**0*xs2**2*xst**9&
     +xp(164)*xs1**0*xs2**4*xst**7&
     +xp(165)*xs1**0*xs2**6*xst**5&
     +xp(166)*xs1**0*xs2**8*xst**3&
     +xp(167)*xs1**0*xs2**10*xst**1&
     +xp(168)*xs1**1*xs2**0*xst**10&
     +xp(169)*xs1**1*xs2**2*xst**8&
     +xp(170)*xs1**1*xs2**4*xst**6&
     +xp(171)*xs1**1*xs2**6*xst**4&
     +xp(172)*xs1**1*xs2**8*xst**2&
     +xp(173)*xs1**1*xs2**10*xst**0&
     +xp(174)*xs1**2*xs2**0*xst**9
 vp3= xp(175)*xs1**2*xs2**2*xst**7&
     +xp(176)*xs1**2*xs2**4*xst**5&
     +xp(177)*xs1**2*xs2**6*xst**3&
     +xp(178)*xs1**2*xs2**8*xst**1&
     +xp(179)*xs1**3*xs2**0*xst**8&
     +xp(180)*xs1**3*xs2**2*xst**6&
     +xp(181)*xs1**3*xs2**4*xst**4&
     +xp(182)*xs1**3*xs2**6*xst**2&
     +xp(183)*xs1**3*xs2**8*xst**0&
     +xp(184)*xs1**4*xs2**0*xst**7&
     +xp(185)*xs1**4*xs2**2*xst**5&
     +xp(186)*xs1**4*xs2**4*xst**3&
     +xp(187)*xs1**4*xs2**6*xst**1&
     +xp(188)*xs1**5*xs2**0*xst**6&
     +xp(189)*xs1**5*xs2**2*xst**4&
     +xp(190)*xs1**5*xs2**4*xst**2&
     +xp(191)*xs1**5*xs2**6*xst**0&
     +xp(192)*xs1**6*xs2**0*xst**5&
     +xp(193)*xs1**6*xs2**2*xst**3&
     +xp(194)*xs1**6*xs2**4*xst**1&
     +xp(195)*xs1**7*xs2**0*xst**4&
     +xp(196)*xs1**7*xs2**2*xst**2&
     +xp(197)*xs1**7*xs2**4*xst**0&
     +xp(198)*xs1**8*xs2**0*xst**3&
     +xp(199)*xs1**8*xs2**2*xst**1&
     +xp(200)*xs1**9*xs2**0*xst**2&
     +xp(201)*xs1**9*xs2**2*xst**0&
     +xp(202)*xs1**0*xs2**0*xst**12&
     +xp(203)*xs1**0*xs2**2*xst**10&
     +xp(204)*xs1**0*xs2**4*xst**8&
     +xp(205)*xs1**0*xs2**6*xst**6&
     +xp(206)*xs1**0*xs2**8*xst**4&
     +xp(207)*xs1**0*xs2**10*xst**2&
     +xp(208)*xs1**0*xs2**12*xst**0&
     +xp(209)*xs1**1*xs2**0*xst**11&    
     +xp(210)*xs1**1*xs2**2*xst**9&
     +xp(211)*xs1**1*xs2**4*xst**7&
     +xp(212)*xs1**1*xs2**6*xst**5&
     +xp(213)*xs1**1*xs2**8*xst**3&
     +xp(214)*xs1**1*xs2**10*xst**1&
     +xp(215)*xs1**2*xs2**0*xst**10&
     +xp(216)*xs1**2*xs2**2*xst**8&
     +xp(217)*xs1**2*xs2**4*xst**6&
     +xp(218)*xs1**2*xs2**6*xst**4&
     +xp(219)*xs1**2*xs2**8*xst**2&
     +xp(220)*xs1**2*xs2**10*xst**0&
     +xp(221)*xs1**3*xs2**0*xst**9&
     +xp(222)*xs1**3*xs2**2*xst**7&
     +xp(223)*xs1**3*xs2**4*xst**5&
     +xp(224)*xs1**3*xs2**6*xst**3&
     +xp(225)*xs1**3*xs2**8*xst**1&
     +xp(226)*xs1**4*xs2**0*xst**8&
     +xp(227)*xs1**4*xs2**2*xst**6&
     +xp(228)*xs1**4*xs2**4*xst**4&
     +xp(229)*xs1**4*xs2**6*xst**2&
     +xp(230)*xs1**4*xs2**8*xst**0&
     +xp(231)*xs1**5*xs2**0*xst**7&
     +xp(232)*xs1**5*xs2**2*xst**5&
     +xp(233)*xs1**5*xs2**4*xst**3&
     +xp(234)*xs1**5*xs2**6*xst**1&
     +xp(235)*xs1**6*xs2**0*xst**6&
     +xp(236)*xs1**6*xs2**2*xst**4&
     +xp(237)*xs1**6*xs2**4*xst**2&
     +xp(238)*xs1**6*xs2**6*xst**0&
     +xp(239)*xs1**7*xs2**0*xst**5&
     +xp(240)*xs1**7*xs2**2*xst**3


     vp=vp1+vp2+vp3
     !
     vps1=42395.535333_ark*xep1
     vps2=42395.535333_ark*xep2

     y1=1._ark/(1._ark+exp(ut*(x0-r1)))
     y2=1._ark/(1._ark+exp(ut*(x0-r2)))
     y12=1._ark/(1._ark+exp(ut2*(x02-r1)))
     y22=1._ark/(1._ark+exp(ut2*(x02-r2)))

     vp=vp*xep3*(1-y12)*(1-y22)
     voh1=vpb1*(1-y1)+y1*vps1
     voh2=vpb2*(1-y2)+y2*vps2
     !
     v=v0+vp+voh1+voh2+vhh
     !
  end function  MLpoten_h2o_tennyson



  !
  ! Defining potential energy function 
  !
  function MLpoten_xy2_tyuterev_damp(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)            :: r12,r32,alpha,xcos,v,v1,v2,v3,v4,v5,v6,v7,v8

   real(ark)            :: aa1,re12,alphae,xst,y1,y2,y3,xs1,xs2,v0,vp1,vp2,vp3, y(3)
   real(ark)            :: g1,g2,b1,b2,de,beta,damp2,damp4,vdamp,vlong,dampa
   integer(ik)          :: N
   real(ark)            :: ycos,v_t,q1,q2
   !
   if (verbose>=6) write(out,"('MLpoten_xy2_tyuterev_damp/start')") 
   !
   r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
   !
   re12    = molec%req(1)
   alphae  = molec%alphaeq(1)
   aa1     = molec%specparam(1)
   !
   De     = force(1)
   beta   = force(2)
   damp2  = force(3)
   damp4  = force(4)
   dampa  = damp2
   !
   y(1) = r12-re12
   y(2) = r32-re12
   y(3) = alpha-alphae
   !
   vlong = De*(1.0_ark-exp(-beta*y(1)))**2 + De*(1.0_ark-exp(-beta*y(2)))**2
   vdamp = exp( -damp2*( y(1)**2+y(2)**2)-0.0*damp4*( y(1)**4+y(2)**4 )-0.0d0*dampa*y(3)**2 )
   !
   ! calculate potential energy function values
   !
   y1=r12-re12 !1.0d+00-exp(-aa1*(r12-re12))
   y2=r32-re12 !1.0d+00-exp(-aa1*(r32-re12))
   !
   y3=cos(alpha)-cos(alphae)
   !
   N = size(force)
   !
   v4 = 0 ; v5 = 0 ; v6 = 0 ; v7 = 0 ; v8 = 0
   !
 v0 = force(5)*y1**0*y2**0*y3**0
 v1 = force(6)*y1**0*y2**0*y3**1& 
    + force(7)*y1**1*y2**0*y3**0& 
    + force(7)*y1**0*y2**1*y3**0
 v2 = force(8)*y1**0*y2**0*y3**2& 
    + force(9)*y1**1*y2**0*y3**1& 
    + force(9)*y1**0*y2**1*y3**1& 
    + force(10)*y1**1*y2**1*y3**0& 
    + force(11)*y1**2*y2**0*y3**0& 
    + force(11)*y1**0*y2**2*y3**0
 v3 = force(12)*y1**0*y2**0*y3**3& 
    + force(13)*y1**1*y2**0*y3**2& 
    + force(13)*y1**0*y2**1*y3**2& 
    + force(14)*y1**1*y2**1*y3**1& 
    + force(15)*y1**2*y2**0*y3**1& 
    + force(15)*y1**0*y2**2*y3**1& 
    + force(16)*y1**2*y2**1*y3**0& 
    + force(16)*y1**1*y2**2*y3**0& 
    + force(17)*y1**3*y2**0*y3**0& 
    + force(17)*y1**0*y2**3*y3**0

 if (N>18) then 
  v4 = force(18)*y1**0*y2**0*y3**4&
    + force(19)*y1**1*y2**0*y3**3& 
    + force(19)*y1**0*y2**1*y3**3& 
    + force(20)*y1**1*y2**1*y3**2& 
    + force(21)*y1**2*y2**0*y3**2& 
    + force(21)*y1**0*y2**2*y3**2& 
    + force(22)*y1**2*y2**1*y3**1& 
    + force(22)*y1**1*y2**2*y3**1& 
    + force(23)*y1**2*y2**2*y3**0&
    + force(24)*y1**3*y2**0*y3**1& 
    + force(24)*y1**0*y2**3*y3**1& 
    + force(25)*y1**3*y2**1*y3**0& 
    + force(25)*y1**1*y2**3*y3**0& 
    + force(26)*y1**4*y2**0*y3**0& 
    + force(26)*y1**0*y2**4*y3**0
endif

 if (N>26) then 
  v5 = force(27)*y1**0*y2**0*y3**5& 
    + force(28)*y1**1*y2**0*y3**4& 
    + force(28)*y1**0*y2**1*y3**4& 
    + force(29)*y1**1*y2**1*y3**3& 
    + force(30)*y1**2*y2**0*y3**3& 
    + force(30)*y1**0*y2**2*y3**3& 
    + force(31)*y1**2*y2**1*y3**2& 
    + force(31)*y1**1*y2**2*y3**2& 
    + force(32)*y1**2*y2**2*y3**1& 
    + force(33)*y1**3*y2**0*y3**2& 
    + force(33)*y1**0*y2**3*y3**2& 
    + force(34)*y1**3*y2**1*y3**1& 
    + force(34)*y1**1*y2**3*y3**1& 
    + force(35)*y1**3*y2**2*y3**0& 
    + force(35)*y1**2*y2**3*y3**0& 
    + force(36)*y1**4*y2**0*y3**1& 
    + force(36)*y1**0*y2**4*y3**1& 
    + force(37)*y1**4*y2**1*y3**0& 
    + force(37)*y1**1*y2**4*y3**0& 
    + force(38)*y1**5*y2**0*y3**0& 
    + force(38)*y1**0*y2**5*y3**0
endif

 if (N>38) then 
  v6 = force(39)*y1**0*y2**0*y3**6& 
    + force(40)*y1**1*y2**0*y3**5&
    + force(40)*y1**0*y2**1*y3**5& 
    + force(41)*y1**1*y2**1*y3**4& 
    + force(42)*y1**2*y2**0*y3**4& 
    + force(42)*y1**0*y2**2*y3**4& 
    + force(43)*y1**2*y2**1*y3**3& 
    + force(43)*y1**1*y2**2*y3**3& 
    + force(44)*y1**2*y2**2*y3**2& 
    + force(45)*y1**3*y2**0*y3**3&
    + force(45)*y1**0*y2**3*y3**3& 
    + force(46)*y1**3*y2**1*y3**2& 
    + force(46)*y1**1*y2**3*y3**2& 
    + force(47)*y1**3*y2**2*y3**1& 
    + force(47)*y1**2*y2**3*y3**1& 
    + force(48)*y1**3*y2**3*y3**0& 
    + force(49)*y1**4*y2**0*y3**2& 
    + force(49)*y1**0*y2**4*y3**2& 
    + force(50)*y1**4*y2**1*y3**1& 
    + force(50)*y1**1*y2**4*y3**1& 
    + force(51)*y1**4*y2**2*y3**0& 
    + force(51)*y1**2*y2**4*y3**0& 
    + force(52)*y1**5*y2**0*y3**1& 
    + force(52)*y1**0*y2**5*y3**1& 
    + force(53)*y1**5*y2**1*y3**0& 
    + force(53)*y1**1*y2**5*y3**0& 
    + force(54)*y1**6*y2**0*y3**0& 
    + force(54)*y1**0*y2**6*y3**0
 endif

 if (N>54) then 
 v7 = force(55)*y1**0*y2**0*y3**7& 
    + force(56)*y1**1*y2**0*y3**6& 
    + force(56)*y1**0*y2**1*y3**6& 
    + force(57)*y1**1*y2**1*y3**5& 
    + force(58)*y1**2*y2**0*y3**5& 
    + force(58)*y1**0*y2**2*y3**5& 
    + force(59)*y1**2*y2**1*y3**4& 
    + force(59)*y1**1*y2**2*y3**4& 
    + force(60)*y1**2*y2**2*y3**3& 
    + force(61)*y1**3*y2**0*y3**4& 
    + force(61)*y1**0*y2**3*y3**4& 
    + force(62)*y1**3*y2**1*y3**3& 
    + force(62)*y1**1*y2**3*y3**3& 
    + force(63)*y1**3*y2**2*y3**2&
    + force(63)*y1**2*y2**3*y3**2& 
    + force(64)*y1**3*y2**3*y3**1& 
    + force(65)*y1**4*y2**0*y3**3& 
    + force(65)*y1**0*y2**4*y3**3& 
    + force(66)*y1**4*y2**1*y3**2& 
    + force(66)*y1**1*y2**4*y3**2& 
    + force(67)*y1**4*y2**2*y3**1& 
    + force(67)*y1**2*y2**4*y3**1&
    + force(68)*y1**4*y2**3*y3**0& 
    + force(68)*y1**3*y2**4*y3**0& 
    + force(69)*y1**5*y2**0*y3**2& 
    + force(69)*y1**0*y2**5*y3**2& 
    + force(70)*y1**5*y2**1*y3**1& 
    + force(70)*y1**1*y2**5*y3**1& 
    + force(71)*y1**5*y2**2*y3**0& 
    + force(71)*y1**2*y2**5*y3**0& 
    + force(72)*y1**6*y2**0*y3**1& 
    + force(72)*y1**0*y2**6*y3**1& 
    + force(73)*y1**6*y2**1*y3**0& 
    + force(73)*y1**1*y2**6*y3**0& 
    + force(74)*y1**7*y2**0*y3**0& 
    + force(74)*y1**0*y2**7*y3**0
 endif

 if (N>74) then 
 v8 = force(75)*y1**0*y2**0*y3**8& 
    + force(76)*y1**1*y2**0*y3**7& 
    + force(76)*y1**0*y2**1*y3**7& 
    + force(77)*y1**1*y2**1*y3**6& 
    + force(78)*y1**2*y2**0*y3**6& 
    + force(78)*y1**0*y2**2*y3**6& 
    + force(79)*y1**2*y2**1*y3**5& 
    + force(79)*y1**1*y2**2*y3**5& 
    + force(80)*y1**2*y2**2*y3**4& 
    + force(81)*y1**3*y2**0*y3**5& 
    + force(81)*y1**0*y2**3*y3**5& 
    + force(82)*y1**3*y2**1*y3**4& 
    + force(82)*y1**1*y2**3*y3**4& 
    + force(83)*y1**3*y2**2*y3**3& 
    + force(83)*y1**2*y2**3*y3**3& 
    + force(84)*y1**3*y2**3*y3**2& 
    + force(85)*y1**4*y2**0*y3**4& 
    + force(85)*y1**0*y2**4*y3**4&
    + force(86)*y1**4*y2**1*y3**3& 
    + force(86)*y1**1*y2**4*y3**3& 
    + force(87)*y1**4*y2**2*y3**2& 
    + force(87)*y1**2*y2**4*y3**2& 
    + force(88)*y1**4*y2**3*y3**1& 
    + force(88)*y1**3*y2**4*y3**1& 
    + force(89)*y1**4*y2**4*y3**0& 
    + force(90)*y1**5*y2**0*y3**3&
    + force(90)*y1**0*y2**5*y3**3& 
    + force(91)*y1**5*y2**1*y3**2& 
    + force(91)*y1**1*y2**5*y3**2& 
    + force(92)*y1**5*y2**2*y3**1& 
    + force(92)*y1**2*y2**5*y3**1& 
    + force(93)*y1**5*y2**3*y3**0& 
    + force(93)*y1**3*y2**5*y3**0& 
    + force(94)*y1**6*y2**0*y3**2& 
    + force(94)*y1**0*y2**6*y3**2& 
    + force(95)*y1**6*y2**1*y3**1& 
    + force(95)*y1**1*y2**6*y3**1& 
    + force(96)*y1**6*y2**2*y3**0& 
    + force(96)*y1**2*y2**6*y3**0& 
    + force(97)*y1**7*y2**0*y3**1& 
    + force(97)*y1**0*y2**7*y3**1& 
    + force(98)*y1**7*y2**1*y3**0& 
    + force(98)*y1**1*y2**7*y3**0& 
    + force(99)*y1**8*y2**0*y3**0& 
    + force(99)*y1**0*y2**8*y3**0
endif


    f=(v0+v1+v2+v3+v4+v5+v6+v7+v8)*vdamp+vlong
    !
    !f = vlong + vshort*vdamp
    !
    !q1 = local(1)/bohr ; q2 = local(2)/bohr ;  ycos = cos(alpha)
    !
    !v_t = 0 
    !
    !call potv(v_t,q1,q2,ycos)
    !
    !f = f + v_t*219474.630670_ark
    !
    if (verbose>=6) write(out,"('MLpoten_xy2_tyuterev_damp/end')") 
 
 end function MLpoten_xy2_tyuterev_damp



 function MLpoten_SO2_pes_8d(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   !real(ark),intent(in) ::  force(:)

   integer(ik) ifst,nparams
   real(rk) y1, y2, y3, V0, V1, V2, V3, V_damp, V_long, re, ae, De, beta, damp2, damp4
   real(rk) F000,  F100,  F200,  F300,  F400,  F500,  F600,  F700,  F800,  F010,  F020,  F030,  F040, &
       F050,  F060,  F070,  F080,  F110,  F210,  F310,  F410,  F510,  F610,  F710,  F011,  F120, &
       F220,  F320,  F420,  F520,  F620,  F021,  F022,  F130,  F230,  F330,  F430,  F530,  F031, &
       F032,  F033,  F140,  F240,  F340,  F440,  F041,  F042,  F043,  F044,  F150,  F250,  F350, &
       F051,  F052,  F053,  F160,  F260,  F061,  F062,  F170,  F071,  F111,  F211,  F311,  F411, &
       F511,  F611,  F121,  F221,  F321,  F421,  F521,  F122,  F222,  F322,  F422,  F131,  F231, &
       F331,  F431,  F132,  F232,  F332,  F133,  F233,  F141,  F241,  F341,  F142,  F242,  F143, &
       F151,  F251,  F152,  F161

   real(ark), parameter :: rad = 57.29577951308232088_ark

    nparams = size(force)

    re     = force(1)
    ae     = force(2)/rad
    De     = force(3)
    beta   = force(4)
    damp2  = force(5)
    damp4  = force(6)
    ifst   = 6
    
    y1 = local(1)-re
    y2 = local(2)-re
    y3 = ( local(3)-ae )
    
    V_long = De*(1.0d0-exp(-beta*y1))**2 + De*(1.0d0-exp(-beta*y2))**2
    V_damp = exp(-damp2*y1**2 -damp4*y1**4 -damp2*y2**2 -damp4*y2**4)
    
    if (nparams>=(1+ifst)) then
     F000 = force(1+ifst)
     V0=F000
    endif
    if (nparams>=(17+ifst)) then
     F100 = force(2+ifst)
     F200 = force(3+ifst)
     F300 = force(4+ifst)
     F400 = force(5+ifst)
     F500 = force(6+ifst)
     F600 = force(7+ifst)
     F700 = force(8+ifst)
     F800 = force(9+ifst)
     F010 = force(10+ifst)
     F020 = force(11+ifst)
     F030 = force(12+ifst)
     F040 = force(13+ifst)
     F050 = force(14+ifst)
     F060 = force(15+ifst)
     F070 = force(16+ifst)
     F080 = force(17+ifst)
     V1=F100*y3 + F200*y3**2 + F300*y3**3 + F400*y3**4 + F500*y3**5 + F600*y3**6 + &
    F700*y3**7 + F800*y3**8 + F010*(y1 + y2) + F020*(y1**2 + y2**2) + F030*(y1**3 + &
    y2**3) + F040*(y1**4 + y2**4) + F050*(y1**5 + y2**5) + F060*(y1**6 + y2**6) + &
    F070*(y1**7 + y2**7) + F080*(y1**8 + y2**8)
    endif
    if (nparams>=(61+ifst)) then
    F110 = force(18+ifst)
    F210 = force(19+ifst)
    F310 = force(20+ifst)
    F410 = force(21+ifst)
    F510 = force(22+ifst)
    F610 = force(23+ifst)
    F710 = force(24+ifst)
    F011 = force(25+ifst)
    F120 = force(26+ifst)
    F220 = force(27+ifst)
    F320 = force(28+ifst)
    F420 = force(29+ifst)
    F520 = force(30+ifst)
    F620 = force(31+ifst)
    F021 = force(32+ifst)
    F022 = force(33+ifst)
    F130 = force(34+ifst)
    F230 = force(35+ifst)
    F330 = force(36+ifst)
    F430 = force(37+ifst)
    F530 = force(38+ifst)
    F031 = force(39+ifst)
    F032 = force(40+ifst)
    F033 = force(41+ifst)
    F140 = force(42+ifst)
    F240 = force(43+ifst)
    F340 = force(44+ifst)
    F440 = force(45+ifst)
    F041 = force(46+ifst)
    F042 = force(47+ifst)
    F043 = force(48+ifst)
    F044 = force(49+ifst)
    F150 = force(50+ifst)
    F250 = force(51+ifst)
    F350 = force(52+ifst)
    F051 = force(53+ifst)
    F052 = force(54+ifst)
    F053 = force(55+ifst)
    F160 = force(56+ifst)
    F260 = force(57+ifst)
    F061 = force(58+ifst)
    F062 = force(59+ifst)
    F170 = force(60+ifst)
    F071 = force(61+ifst)
    V2=F011*y1*y2 + F022*y1**2*y2**2 + F033*y1**3*y2**3 + F044*y1**4*y2**4 + &
   F110*y3*(y1 + y2) + F210*y3**2*(y1 + y2) + F310*y3**3*(y1 + y2) + F410*y3**4*(y1 &
   + y2) + F510*y3**5*(y1 + y2) + F610*y3**6*(y1 + y2) + F710*y3**7*(y1 + y2) + &
   F120*y3*(y1**2 + y2**2) + F220*y3**2*(y1**2 + y2**2) + F320*y3**3*(y1**2 + &
   y2**2) + F420*y3**4*(y1**2 + y2**2) + F520*y3**5*(y1**2 + y2**2) + &
   F620*y3**6*(y1**2 + y2**2) + F021*(y1**2*y2 + y1*y2**2) + F130*y3*(y1**3 + &
   y2**3) + F230*y3**2*(y1**3 + y2**3) + F330*y3**3*(y1**3 + y2**3) + &
   F430*y3**4*(y1**3 + y2**3) + F530*y3**5*(y1**3 + y2**3) + F031*(y1**3*y2 + &
   y1*y2**3) + F032*(y1**3*y2**2 + y1**2*y2**3) + F140*y3*(y1**4 + y2**4) + &
   F240*y3**2*(y1**4 + y2**4) + F340*y3**3*(y1**4 + y2**4) + F440*y3**4*(y1**4 + &
   y2**4) + F041*(y1**4*y2 + y1*y2**4) + F042*(y1**4*y2**2 + y1**2*y2**4) + &
   F043*(y1**4*y2**3 + y1**3*y2**4) + F150*y3*(y1**5 + y2**5) + F250*y3**2*(y1**5 + &
   y2**5) + F350*y3**3*(y1**5 + y2**5) + F051*(y1**5*y2 + y1*y2**5) + &
   F052*(y1**5*y2**2 + y1**2*y2**5) + F053*(y1**5*y2**3 + y1**3*y2**5) + &
   F160*y3*(y1**6 + y2**6) + F260*y3**2*(y1**6 + y2**6) + F061*(y1**6*y2 + &
   y1*y2**6) + F062*(y1**6*y2**2 + y1**2*y2**6) + F170*y3*(y1**7 + y2**7) + &
   F071*(y1**7*y2 + y1*y2**7)
   endif
   if (nparams>=(95+ifst)) then
    F111 = force(62+ifst)
    F211 = force(63+ifst)
    F311 = force(64+ifst)
    F411 = force(65+ifst)
    F511 = force(66+ifst)
    F611 = force(67+ifst)
    F121 = force(68+ifst)
    F221 = force(69+ifst)
    F321 = force(70+ifst)
    F421 = force(71+ifst)
    F521 = force(72+ifst)
    F122 = force(73+ifst)
    F222 = force(74+ifst)
    F322 = force(75+ifst)
    F422 = force(76+ifst)
    F131 = force(77+ifst)
    F231 = force(78+ifst)
    F331 = force(79+ifst)
    F431 = force(80+ifst)
    F132 = force(81+ifst)
    F232 = force(82+ifst)
    F332 = force(83+ifst)
    F133 = force(84+ifst)
    F233 = force(85+ifst)
    F141 = force(86+ifst)
    F241 = force(87+ifst)
    F341 = force(88+ifst)
    F142 = force(89+ifst)
    F242 = force(90+ifst)
    F143 = force(91+ifst)
    F151 = force(92+ifst)
    F251 = force(93+ifst)
    F152 = force(94+ifst)
    F161 = force(95+ifst)
    V3=F111*y3*y1*y2 + F211*y3**2*y1*y2 + F311*y3**3*y1*y2 + F411*y3**4*y1*y2 + &
   F511*y3**5*y1*y2 + F611*y3**6*y1*y2 + F122*y3*y1**2*y2**2 + &
   F222*y3**2*y1**2*y2**2 + F322*y3**3*y1**2*y2**2 + F422*y3**4*y1**2*y2**2 + &
   F133*y3*y1**3*y2**3 + F233*y3**2*y1**3*y2**3 + F121*y3*(y1**2*y2 + y1*y2**2) + &
   F221*y3**2*(y1**2*y2 + y1*y2**2) + F321*y3**3*(y1**2*y2 + y1*y2**2) + &
   F421*y3**4*(y1**2*y2 + y1*y2**2) + F521*y3**5*(y1**2*y2 + y1*y2**2) + &
   F131*y3*(y1**3*y2 + y1*y2**3) + F231*y3**2*(y1**3*y2 + y1*y2**3) + &
   F331*y3**3*(y1**3*y2 + y1*y2**3) + F431*y3**4*(y1**3*y2 + y1*y2**3) + &
   F132*y3*(y1**3*y2**2 + y1**2*y2**3) + F232*y3**2*(y1**3*y2**2 + y1**2*y2**3) + &
   F332*y3**3*(y1**3*y2**2 + y1**2*y2**3) + F141*y3*(y1**4*y2 + y1*y2**4) + &
   F241*y3**2*(y1**4*y2 + y1*y2**4) + F341*y3**3*(y1**4*y2 + y1*y2**4) + &
   F142*y3*(y1**4*y2**2 + y1**2*y2**4) + F242*y3**2*(y1**4*y2**2 + y1**2*y2**4) + &
   F143*y3*(y1**4*y2**3 + y1**3*y2**4) + F151*y3*(y1**5*y2 + y1*y2**5) + &
   F251*y3**2*(y1**5*y2 + y1*y2**5) + F152*y3*(y1**5*y2**2 + y1**2*y2**5) + &
   F161*y3*(y1**6*y2 + y1*y2**6)
   endif
   
   f = (V0+V1+V2+V3)*V_damp + V_long
   
 end function MLpoten_SO2_pes_8d
   
     
  !
  ! Defining potential energy function 
  !
  function MLpoten_xy2_schwenke(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)            :: r12,r32,alpha,xcos,v,v1,v2,v3,v4,v5,v6,v7,v8

   real(ark)            :: aa1,re12,alphae,xst,y1,y2,y3,xs1,xs2,v0,r0,rr1,rr2
   real(ark)            :: rhh,vhh,D,a,AA,b,va1,va2,beta,dump
   integer(ik)          :: N
   !
   if (verbose>=6) write(out,"('MLpoten_xy2_schwenke/start')") 
   !
   r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
   !
   re12    = molec%req(1)
   alphae  = molec%alphaeq(1)
   a       = molec%specparam(1)
   !
   D   = force(1)
   r0  = force(2)
   AA  = force(3)
   b   = force(4)
   !
   rhh=sqrt(r12**2+r32**2-2.d0*r12*r32*cos(alpha))
   vhh=AA*exp(-b*rhh)
   !
   rr1 = r12-r0
   rr2 = r32-r0
   !
   Va1=D*(exp(-2.d0*a*rr1)-2.d0*exp(-a*rr1)+1.d0)
   Va2=D*(exp(-2.d0*a*rr2)-2.d0*exp(-a*rr2)+1.d0)
   !
   beta = 2.0d0
   !
   ! calculate potential energy function values
   !
   dump=exp(-beta*(r12-re12)**2-beta*(r32-re12)**2)
   !
   y1 = (r12-re12)/re12
   y2 = (r32-re12)/re12
   !
   y3=cos(alpha)-cos(alphae)
   !
   N = size(force)
   !
   v4 = 0 ; v5 = 0 ; v6 = 0 ; v7 = 0 ; v8 = 0
   !
 v0 = force(5)*y1**0*y2**0*y3**0
 v1 = force(6)*y1**0*y2**0*y3**1& 
    + force(7)*y1**1*y2**0*y3**0& 
    + force(7)*y1**0*y2**1*y3**0
 v2 = force(8)*y1**0*y2**0*y3**2& 
    + force(9)*y1**1*y2**0*y3**1& 
    + force(9)*y1**0*y2**1*y3**1& 
    + force(10)*y1**1*y2**1*y3**0& 
    + force(11)*y1**2*y2**0*y3**0& 
    + force(11)*y1**0*y2**2*y3**0
 v3 = force(12)*y1**0*y2**0*y3**3& 
    + force(13)*y1**1*y2**0*y3**2& 
    + force(13)*y1**0*y2**1*y3**2& 
    + force(14)*y1**1*y2**1*y3**1& 
    + force(15)*y1**2*y2**0*y3**1& 
    + force(15)*y1**0*y2**2*y3**1& 
    + force(16)*y1**2*y2**1*y3**0& 
    + force(16)*y1**1*y2**2*y3**0& 
    + force(17)*y1**3*y2**0*y3**0& 
    + force(17)*y1**0*y2**3*y3**0

 if (N>18) then 
  v4 = force(18)*y1**0*y2**0*y3**4& 
    + force(19)*y1**1*y2**0*y3**3& 
    + force(19)*y1**0*y2**1*y3**3& 
    + force(20)*y1**1*y2**1*y3**2& 
    + force(21)*y1**2*y2**0*y3**2& 
    + force(21)*y1**0*y2**2*y3**2& 
    + force(22)*y1**2*y2**1*y3**1& 
    + force(22)*y1**1*y2**2*y3**1& 
    + force(23)*y1**2*y2**2*y3**0& 
    + force(24)*y1**3*y2**0*y3**1& 
    + force(24)*y1**0*y2**3*y3**1& 
    + force(25)*y1**3*y2**1*y3**0& 
    + force(25)*y1**1*y2**3*y3**0& 
    + force(26)*y1**4*y2**0*y3**0& 
    + force(26)*y1**0*y2**4*y3**0
endif

 if (N>26) then 
  v5 = force(27)*y1**0*y2**0*y3**5& 
    + force(28)*y1**1*y2**0*y3**4& 
    + force(28)*y1**0*y2**1*y3**4& 
    + force(29)*y1**1*y2**1*y3**3& 
    + force(30)*y1**2*y2**0*y3**3& 
    + force(30)*y1**0*y2**2*y3**3& 
    + force(31)*y1**2*y2**1*y3**2& 
    + force(31)*y1**1*y2**2*y3**2& 
    + force(32)*y1**2*y2**2*y3**1& 
    + force(33)*y1**3*y2**0*y3**2& 
    + force(33)*y1**0*y2**3*y3**2& 
    + force(34)*y1**3*y2**1*y3**1& 
    + force(34)*y1**1*y2**3*y3**1& 
    + force(35)*y1**3*y2**2*y3**0& 
    + force(35)*y1**2*y2**3*y3**0& 
    + force(36)*y1**4*y2**0*y3**1& 
    + force(36)*y1**0*y2**4*y3**1& 
    + force(37)*y1**4*y2**1*y3**0& 
    + force(37)*y1**1*y2**4*y3**0& 
    + force(38)*y1**5*y2**0*y3**0& 
    + force(38)*y1**0*y2**5*y3**0
endif

 if (N>38) then 
  v6 = force(39)*y1**0*y2**0*y3**6& 
    + force(40)*y1**1*y2**0*y3**5& 
    + force(40)*y1**0*y2**1*y3**5& 
    + force(41)*y1**1*y2**1*y3**4& 
    + force(42)*y1**2*y2**0*y3**4& 
    + force(42)*y1**0*y2**2*y3**4& 
    + force(43)*y1**2*y2**1*y3**3& 
    + force(43)*y1**1*y2**2*y3**3& 
    + force(44)*y1**2*y2**2*y3**2& 
    + force(45)*y1**3*y2**0*y3**3& 
    + force(45)*y1**0*y2**3*y3**3& 
    + force(46)*y1**3*y2**1*y3**2& 
    + force(46)*y1**1*y2**3*y3**2& 
    + force(47)*y1**3*y2**2*y3**1& 
    + force(47)*y1**2*y2**3*y3**1& 
    + force(48)*y1**3*y2**3*y3**0& 
    + force(49)*y1**4*y2**0*y3**2& 
    + force(49)*y1**0*y2**4*y3**2& 
    + force(50)*y1**4*y2**1*y3**1& 
    + force(50)*y1**1*y2**4*y3**1& 
    + force(51)*y1**4*y2**2*y3**0& 
    + force(51)*y1**2*y2**4*y3**0& 
    + force(52)*y1**5*y2**0*y3**1& 
    + force(52)*y1**0*y2**5*y3**1& 
    + force(53)*y1**5*y2**1*y3**0& 
    + force(53)*y1**1*y2**5*y3**0& 
    + force(54)*y1**6*y2**0*y3**0& 
    + force(54)*y1**0*y2**6*y3**0
 endif

 if (N>54) then 
 v7 = force(55)*y1**0*y2**0*y3**7& 
    + force(56)*y1**1*y2**0*y3**6& 
    + force(56)*y1**0*y2**1*y3**6& 
    + force(57)*y1**1*y2**1*y3**5& 
    + force(58)*y1**2*y2**0*y3**5& 
    + force(58)*y1**0*y2**2*y3**5& 
    + force(59)*y1**2*y2**1*y3**4& 
    + force(59)*y1**1*y2**2*y3**4& 
    + force(60)*y1**2*y2**2*y3**3& 
    + force(61)*y1**3*y2**0*y3**4& 
    + force(61)*y1**0*y2**3*y3**4& 
    + force(62)*y1**3*y2**1*y3**3& 
    + force(62)*y1**1*y2**3*y3**3& 
    + force(63)*y1**3*y2**2*y3**2& 
    + force(63)*y1**2*y2**3*y3**2& 
    + force(64)*y1**3*y2**3*y3**1& 
    + force(65)*y1**4*y2**0*y3**3& 
    + force(65)*y1**0*y2**4*y3**3& 
    + force(66)*y1**4*y2**1*y3**2& 
    + force(66)*y1**1*y2**4*y3**2& 
    + force(67)*y1**4*y2**2*y3**1& 
    + force(67)*y1**2*y2**4*y3**1& 
    + force(68)*y1**4*y2**3*y3**0& 
    + force(68)*y1**3*y2**4*y3**0& 
    + force(69)*y1**5*y2**0*y3**2& 
    + force(69)*y1**0*y2**5*y3**2& 
    + force(70)*y1**5*y2**1*y3**1& 
    + force(70)*y1**1*y2**5*y3**1& 
    + force(71)*y1**5*y2**2*y3**0& 
    + force(71)*y1**2*y2**5*y3**0& 
    + force(72)*y1**6*y2**0*y3**1& 
    + force(72)*y1**0*y2**6*y3**1& 
    + force(73)*y1**6*y2**1*y3**0& 
    + force(73)*y1**1*y2**6*y3**0& 
    + force(74)*y1**7*y2**0*y3**0& 
    + force(74)*y1**0*y2**7*y3**0
 endif

 if (N>74) then 
 v8 = force(75)*y1**0*y2**0*y3**8& 
    + force(76)*y1**1*y2**0*y3**7& 
    + force(76)*y1**0*y2**1*y3**7& 
    + force(77)*y1**1*y2**1*y3**6& 
    + force(78)*y1**2*y2**0*y3**6& 
    + force(78)*y1**0*y2**2*y3**6& 
    + force(79)*y1**2*y2**1*y3**5& 
    + force(79)*y1**1*y2**2*y3**5& 
    + force(80)*y1**2*y2**2*y3**4& 
    + force(81)*y1**3*y2**0*y3**5& 
    + force(81)*y1**0*y2**3*y3**5& 
    + force(82)*y1**3*y2**1*y3**4& 
    + force(82)*y1**1*y2**3*y3**4& 
    + force(83)*y1**3*y2**2*y3**3& 
    + force(83)*y1**2*y2**3*y3**3& 
    + force(84)*y1**3*y2**3*y3**2& 
    + force(85)*y1**4*y2**0*y3**4& 
    + force(85)*y1**0*y2**4*y3**4& 
    + force(86)*y1**4*y2**1*y3**3& 
    + force(86)*y1**1*y2**4*y3**3& 
    + force(87)*y1**4*y2**2*y3**2& 
    + force(87)*y1**2*y2**4*y3**2& 
    + force(88)*y1**4*y2**3*y3**1& 
    + force(88)*y1**3*y2**4*y3**1& 
    + force(89)*y1**4*y2**4*y3**0& 
    + force(90)*y1**5*y2**0*y3**3& 
    + force(90)*y1**0*y2**5*y3**3& 
    + force(91)*y1**5*y2**1*y3**2& 
    + force(91)*y1**1*y2**5*y3**2& 
    + force(92)*y1**5*y2**2*y3**1& 
    + force(92)*y1**2*y2**5*y3**1& 
    + force(93)*y1**5*y2**3*y3**0& 
    + force(93)*y1**3*y2**5*y3**0& 
    + force(94)*y1**6*y2**0*y3**2& 
    + force(94)*y1**0*y2**6*y3**2& 
    + force(95)*y1**6*y2**1*y3**1& 
    + force(95)*y1**1*y2**6*y3**1& 
    + force(96)*y1**6*y2**2*y3**0& 
    + force(96)*y1**2*y2**6*y3**0& 
    + force(97)*y1**7*y2**0*y3**1& 
    + force(97)*y1**0*y2**7*y3**1& 
    + force(98)*y1**7*y2**1*y3**0& 
    + force(98)*y1**1*y2**7*y3**0& 
    + force(99)*y1**8*y2**0*y3**0& 
    + force(99)*y1**0*y2**8*y3**0
endif


    f=v0+(v1+v2+v3+v4+v5+v6+v7+v8)*dump+vhh+va1+va2

    !
    if (verbose>=6) write(out,"('MLpoten_xy2_schwenke/end')") 
 
 end function MLpoten_xy2_schwenke




  !
  ! Defining potential energy function (built for SO2)


  function MLpoten_so2_damp(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
      !
      integer(ik)          :: iterm,k_ind(3),ifst
      real(ark)            :: y(3),xi(3),vshort,vlong,vdamp
      real(ark)            :: De,ae,re,beta,damp2,damp4,th1,vshort2,temp
      !
      re     = molec%req(1)
      ae     = molec%alphaeq(1)
      !
      De     = force(1)
      beta   = force(2)
      damp2  = force(3)
      damp4  = force(4)
      !
      ifst   = 5
      !
      y(1) = local(1)-re
      y(2) = local(2)-re
      y(3) = local(3)-ae
      !
      vlong = De*(1.0_ark-exp(-beta*y(1)))**2 + De*(1.0_ark-exp(-beta*y(2)))**2
      vdamp = exp( -damp2*( y(1)**2+y(2)**2)-damp4*( y(1)**4+y(2)**4 ) )
      !
      vshort2 = 0.0_ark
      !
      vshort = 0.0_ark
      !
      do iterm=ifst,size(force)
       !
       xi(1:3) = y(1:3)**molec%pot_ind(1:3,iterm)
       !
       temp = force(iterm)*product(xi(1:3))
       !
       if (molec%pot_ind(1,iterm)/=molec%pot_ind(2,iterm)) then
         !
         k_ind(1) = molec%pot_ind(2,iterm)
         k_ind(2) = molec%pot_ind(1,iterm)
         k_ind(3) = molec%pot_ind(3,iterm)
         !
         xi(1:3) = y(1:3)**k_ind(1:3)
         !
         temp = temp + force(iterm)*product(xi(1:3))
         !
       endif
       !
       if (sum(molec%pot_ind(1:3,iterm))>2) then 
         !
         vshort = vshort + temp
         !
       else 
         !
         vshort2 = vshort2 + temp
         !
       endif
       !
      enddo
      !
      temp = vshort2*vdamp+vlong
      !
      th1 = 0.5d0*( 1.0d0-tanh( 0.0001_ark*( temp-40000.0_ark ) ) )
      !
      f = temp + vshort*th1*vdamp
      !
      !f = vlong + vshort*vdamp
      !
  end function MLpoten_so2_damp


!
!   Defining potential energy function (built for CO2, AMES1, )
!   Ames-1 PES for CO2, see details in
!    Xinchuan Huang, David W. Schwenke, Sergey A. Tashkun, and Timothy J. Lee
!    J. Chem. Phys. 2012 (submitted)
!    
! We refined Ames-0 PES [CCSD(T)/cc-pCVQZ + 0.5*ACPF correction] 
! with selected 12C16O2 rovibrational energy levels
!
! Accuracy statistics of 6873 energy levels computed on Ames-1 PES, 
! compared to purely experimentally determined energies. 
! Overall rms = 0.0156 cm-1. Unit for Energy and rms is cm-1


  function MLpoten_co2_ames1(ncoords,natoms,local,xyz,force) result(f)
   !
   implicit none
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          :: i,k,Ncoe
   real(ark)            :: r12ref,alpha2,De1,De2,Ae1,Ae2,edamp2,edamp4,edamp5,edamp6
   !
   real(ark)            :: emin,rmin,rminbohr,alpha,rref,a2b
   !
   real(ark)            :: str1,str2,enetmp1,atp1,enetmp2,edamp11,edamp12,edamp1,V
   real(ark)            :: v0,rrco1,rrco2,r1,r2,a3,xx1,xx2,rco1,rco2,ang,dstr1,dstr2,sumstr2,sumstr4,angref,angx,bdamp2,bdamp4
   !
   real(ark)            ::   v1(3),v2(3),x1,x2,cosalpha1
   character(len=cl)    :: txt = 'MLpoten_c3_R_theta'
        !
        r12ref  = force(1)
        alpha2  = force(2)
        De1     = force(3)
        De2     = force(4)
        Ae1     = force(5)
        Ae2     = force(6)
        edamp2  = force(7)
        edamp4  = force(8)
        edamp5  = force(9)
        edamp6  = force(10)
        Emin    = force(11)
        Rmin    = force(12)
        rminbohr= force(13)
        alpha   = force(14)
        rref    = force(15)
        !
        rco1 = local(1)
        rco2 = local(2)
        !
        rrco1=rco1/0.529177249_ark
        rrco2=rco2/0.529177249_ark
        !
        ang = pi-local(3)
        if  (molec%Ndihedrals>1) then
          !
          ang = asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
          !
        endif 
        !
        if  (molec%Ndihedrals>1) then
          !
          ang = asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
          !
          v1(:) = xyz(2,:)-xyz(1,:)
          v2(:) = xyz(3,:)-xyz(1,:)
          !
          x1 = sqrt(sum(v1(:)**2))
          x2 = sqrt(sum(v2(:)**2))
          !
          cosalpha1 = sum( v1(:)*v2(:))/(x1*x2)
          !
          !alpha = aacos(cosalpha1,txt)
          ang  = aacos(cosalpha1,txt)-pi
          !
        endif
        !
        r1=1._ark-exp(-alpha*(rco1-r12ref))
        r2=1._ark-exp(-alpha*(rco2-r12ref))
        a3=cos(ang)
        !
        Ncoe = molec%parmax
        !
        v0=0
        do i=16,Ncoe
           !
           v0=v0+force(i)*r1**molec%pot_ind(1,i)*r2**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
           if(molec%pot_ind(1,i).ne.molec%pot_ind(2,i))then
             v0=v0+force(i)*r2**molec%pot_ind(1,i)*r1**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
           end if
           !
        end do
        !
        xx1=rco1; xx2=rco2
        dstr1=(xx1-r12ref)
        dstr2=(xx2-r12ref)
        sumstr2=dstr1**2+dstr2**2
        sumstr4=dstr1**4+dstr2**4
        !
        angref=150.0_ark*pi/180.0_ark-pi
        !
        angx=min(-abs(ang)-angref,0.0_ark)
        bdamp2=angx**2;bdamp4=angx**4
        bdamp2=ang**2; bdamp4=ang**4
        !
        !angx=min(-abs(ang)+angref,0.d0)
        !bdamp2=angx**2;bdamp4=angx**4
        !bdamp2=ang**2; bdamp4=ang**4
        !
        V0=V0*exp(edamp2*sumstr2+edamp4*sumstr4+edamp5*bdamp2+edamp6*bdamp4)
        !
        a2b=0.529177249_ark
        str1=1._ark-exp(-alpha2*(xx1-r12ref))
        str2=1._ark-exp(-alpha2*(xx2-r12ref))
        enetmp1=De1*(str1**2+str2**2)+De2*(str1**4+str2**4)
        enetmp1=enetmp1/219474.63067_ark
        !
        atp1=acos(-1._ark)-abs(ang)
        enetmp2=Ae1*(1._ark+cos(atp1))**2+Ae2*(1._ark+cos(atp1))**4
        enetmp2=enetmp2/219474.63067_ark
        !
        sumstr2=(xx1-r12ref)**2+(xx2-r12ref)**2
        sumstr4=(xx1-r12ref)**4+(xx2-r12ref)**4
        edamp11=-0.5_ark; edamp12=-0.5_ark
        edamp1=exp(0.2_ark*edamp11*sumstr2+0*edamp12*sumstr4)
        enetmp2=edamp1*enetmp2
        !
        V0=V0+enetmp1+enetmp2
        !
        V=V0-emin
        !
        f = V*219474.63067_ark
        !
  end function MLpoten_co2_ames1





  function MLpoten_so2_ames1(ncoords,natoms,local,xyz,force) result(f)
   !
   implicit none
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          :: i,k,Ncoe
   real(ark)            :: r12ref,alpha2,De1,De2,Ae1,Ae2,edamp2,edamp4,edamp5,edamp6
   !
   real(ark)            :: emin,rmin,rminbohr,alpha,rref,a2b
   !
   real(ark)            :: str1,str2,enetmp1,atp1,enetmp2,edamp11,edamp12,edamp1,V
   real(ark)            :: v0,rrso1,rrso2,r1,r2,a3,xx1,xx2,rso1,rso2,ang,dstr1,dstr2,sumstr2,sumstr4,angref,angx,bdamp2,bdamp4,aref
        !
        r12ref  = force(1)
        aref    = force(2)
        alpha2  = force(3)
        De1     = force(4)
        De2     = force(5)
        Ae1     = force(6)
        Ae2     = force(7)
        edamp2  = force(8)
        edamp4  = force(9)
        edamp5  = force(10)
        edamp6  = force(11)
        Emin    = force(12)
        Rmin    = force(13)
        rminbohr= force(14)
        alpha   = force(15)
        rref    = force(16)
        !
        rref = r12ref
        !
        a2b=0.529177249_ark
        !
        rrso1 = local(1)
        rrso2 = local(2)
        !
        !rso1 = rrso1/a2b
        !rso2 = rrso2/a2b
        !
		!rrso1=rso1*a2b
		!rrso2=rso2*a2b
        !
        ang = local(3)
        !
        if  (molec%Ndihedrals>1) then
          !
          ang = pi - asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
          !
        endif 
        !
        r1=rrso1-rref
	    r2=rrso2-rref
        a3=cos(ang)-cos(aref)
        !
        Ncoe = molec%parmax
        !
        v0=0
        do i=17,Ncoe
           !
           v0=v0+force(i)*r1**molec%pot_ind(1,i)*r2**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
           if(molec%pot_ind(1,i).ne.molec%pot_ind(2,i))then
             v0=v0+force(i)*r2**molec%pot_ind(1,i)*r1**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
           end if
           !
        end do
        !
        r12ref=rref
        xx1=rrso1; xx2=rrso2
        dstr1=(xx1-r12ref)
        dstr2=(xx2-r12ref)
        sumstr2=dstr1**2+dstr2**2
        sumstr4=dstr1**4+dstr2**4
        !
        angref=160.d0*acos(-1.d0)/180.d0-aref
        !
        angx=min(-abs(ang-aref)+angref,0.d0)
        !
        bdamp2=angx**2;bdamp4=angx**4
        !
        V0=V0*exp(edamp2*sumstr2+edamp4*sumstr4+edamp5*bdamp2+edamp6*bdamp4)
        a2b=0.529177249_ark 
        str1=1.0_ark-exp(-alpha*(rrso1-rref)/a2b)
        str2=1.0_ark-exp(-alpha*(rrso2-rref)/a2b)
        enetmp1=De1*(str1**2+str2**2)+De2*(str1**4+str2**4)
        enetmp1=enetmp1/219474.63067d0

        atp1=ang
        enetmp2=Ae1*(cos(atp1)-cos(aref))**2+Ae2*(cos(atp1)-cos(aref))**4
        enetmp2=enetmp2/219474.63067_ark

        sumstr2=(xx1-r12ref)**2+(xx2-r12ref)**2
        sumstr4=(xx1-r12ref)**4+(xx2-r12ref)**4
        edamp11=-0.5d0; edamp12=-0.5d0
        edamp1=exp(0.4d0*edamp11*sumstr2+0.d0*edamp12*sumstr4)
        enetmp2=edamp1*enetmp2

        V0=V0+enetmp1+enetmp2

        V=V0-emin

        !
        f = V*219474.63067d0
        !
  end function MLpoten_so2_ames1



 !returns electric dipole moment cartesian coordinates in the user-defined frame for locals specified
 !
 recursive subroutine MLdipole_so2_ames1(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)            :: i,Ndcoe
    real(ark)              :: tmat(3,3),n1(3),n2(3),n3(3),mu(3)
    !
    real(ark)   :: dcoe(1000),idcoe(1000,3)
    real(ark)   :: x(3,3),r1,r2,r3,str1,str2,bend3
    real(ark)   :: afunc(1000),afunc1(1000),afunc2(1000)
    !
    do i=1,3
      x(i,:)=xyz(i,:)-xyz(1,:)
    end do
    !    
    r1=sqrt(sum(x(2,:)**2)) !r12 distance
    r2=sqrt(sum(x(3,:)**2)) !r13 distance
    r3=sqrt(sum((x(2,:)-x(3,:))**2)) !r23 distance
    ! 
    !str1=r1-1.15958d0 !Emil: shifting to 0 when in equilibrium geometry?!
    !str2=r2-1.15958d0
    !
    str1=r1-1.43108_ark
    str2=r2-1.43108_ark
    !
    !bend3=1.d0+(r1**2+r2**2-r3**2)/(2*r1*r2) !Emil@: [0;2] <- cos(alpha)+1
    !
    bend3=(r1**2+r2**2-r3**2)/(2.0_ark*r1*r2)
    bend3=max(-1.0_ark,min(1.0_ark,bend3))-cos(2.0825435274_ark)
    !
    Ndcoe = extF%nterms(1)
    !
    afunc1=0; afunc2=0; afunc=0
    do i=1,Ndcoe !Emil: loop over all geometries
      afunc1(i)=str1**extF%term(1,i,1)*str2**extF%term(2,i,1)*bend3**extF%term(3,i,1)
      afunc2(i)=str2**extF%term(1,i,1)*str1**extF%term(2,i,1)*bend3**extF%term(3,i,1)
    end do
    !
    do i=1,3
      afunc(1:Ndcoe)=afunc1(1:Ndcoe)*x(2,i)+afunc2(1:Ndcoe)*x(3,i)
      mu(i)=dot_product(extF%coef(1:Ndcoe,1),afunc(1:Ndcoe))
    end do
    !
    f = mu*2.541765_ark
    !
 end subroutine MLdipole_so2_ames1



 !returns electric dipole moment cartesian coordinates in the user-defined frame for locals specified
 !
 recursive subroutine MLdipole_ames1(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)            :: i,Ndcoe
    real(ark)              :: tmat(3,3),n1(3),n2(3),n3(3),mu(3),xyz0(3,3)
    !
    real(ark)   :: dcoe(1000),idcoe(1000,3)
    real(ark)   :: x(3,3),r1,r2,r3,str1,str2,bend3,re,ae
    real(ark)   :: afunc(1000),afunc1(1000),afunc2(1000)
    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      select case(trim(molec%coords_transform))
      case default
         write (out,"('MLdipole_ames1: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLdipole_ames1 - bad coord. type'
      case('R-RHO-Z','R-RHO-Z-ECKART')
         !
         xyz0 = MLloc2pqr_xy2(local)
         !
      end select
      !
    else
      !
      xyz0 = xyz
      !
    endif
    !
    do i=1,3
      x(i,:)=xyz0(i,:)-xyz0(1,:)
    end do
    !    
    r1=sqrt(sum(x(2,:)**2)) !r12 distance
    r2=sqrt(sum(x(3,:)**2)) !r13 distance
    r3=sqrt(sum((x(2,:)-x(3,:))**2)) !r23 distance
    ! 
    !str1=r1-1.15958d0 !Emil: shifting to 0 when in equilibrium geometry?!
    !str2=r2-1.15958d0
    !
    re = extF%coef(1,1)
    ae = extF%coef(2,1)*pi/180.0_ark
    !
    str1=r1-re
    str2=r2-re
    !
    !bend3=1.d0+(r1**2+r2**2-r3**2)/(2*r1*r2) !Emil@: [0;2] <- cos(alpha)+1
    !
    bend3=(r1**2+r2**2-r3**2)/(2.0_ark*r1*r2)
    bend3=max(-1.0_ark,min(1.0_ark,bend3))-cos(ae)
    !
    Ndcoe = extF%nterms(1)-2
    !
    afunc1=0; afunc2=0; afunc=0
    do i=1,Ndcoe !Emil: loop over all geometries
      afunc1(i)=str1**extF%term(1,i+2,1)*str2**extF%term(2,i+2,1)*bend3**extF%term(3,i+2,1)
      afunc2(i)=str2**extF%term(1,i+2,1)*str1**extF%term(2,i+2,1)*bend3**extF%term(3,i+2,1)
    end do
    !
    do i=1,3
      afunc(1:Ndcoe)=afunc1(1:Ndcoe)*x(2,i)+afunc2(1:Ndcoe)*x(3,i)
      mu(i)=dot_product(extF%coef(2+1:2+Ndcoe,1),afunc(1:Ndcoe))
    end do
    !
    f = mu/0.393430307_ark
    !
 end subroutine MLdipole_ames1



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
       !call potv(v,r12,r32,xcos)
       v = 0
       !v = v*tocm
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



  function MLpoten_xy2_bubukina(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f

   integer,parameter           ::  npropin = 240

   real(ark) ::  r1,r2,th
   real(ark)            :: reoh,thetae,b1,roh,alphaoh
   real(ark)            :: phh2,t0,ut,x0,ut2,x02,xs1,xs2,xst,rs,rm,rr1,rr2,xep1,xep2,xep3
   real(ark)            :: rhh,vhh,vpb1,vpb2,v0,vp1,vp2,vp3,vp,vps1,vps2,y1,y2,y12,y22,voh1,voh2
   integer(4) :: i 
      !
      r1 = local(1) ; r2 = local(2) ;  th = local(3)
      !
      reoh  = 0.958649_ark
      thetae= 104.3475_ark
      b1=2.0_ark
      roh=0.951961_ark
      alphaoh=2.587949757553683_ark
      phh2=6.70164303995_ark
      t0=0.01_ark
      ut=20._ark
      x0=2.5_ark
      ut2=20._ark
      x02=2.5_ark
        
      thetae=thetae*3.14159265358979312_ark*.00555555555555555555_ark

      xs1=(r1+r2)*0.5_ark-reoh
      xs2=(r1-r2)*0.5_ark
      xst=cos(th)-cos(thetae)

      rs=sqrt(0.5_ark)*(r1+r2)
      rm=sqrt(0.5_ark)*(r1-r2)

      rr1=r1-roh
      rr2=r2-roh

      xep1=exp(-2._ark*alphaoh*rr1)-2._ark*exp(-alphaoh*rr1)+1._ark
      xep2=exp(-2._ark*alphaoh*rr2)-2._ark*exp(-alphaoh*rr2)+1._ark
      xep3=exp(-b1*(rr1**2+rr2**2))
      rhh=sqrt(r1**2+r2**2-2._ark*r1*r2*cos(th))
      vhh=0.900642240911285975d6*exp(-phh2*rhh)
      vpb1=0.518622556959170834d5*xep1
      vpb2=0.518622556959170834d5*xep2


  v0= force(1)  *xs1**0*xs2**0*xst**0
 vp1= force(2)  *xs1**0*xs2**0*xst**1&
     +force(3)  *xs1**1*xs2**0*xst**0&
     +force(4)  *xs1**0*xs2**0*xst**2&
     +force(5)  *xs1**0*xs2**2*xst**0&
     +force(6)  *xs1**1*xs2**0*xst**1&
     +force(7)  *xs1**2*xs2**0*xst**0&
     +force(8)  *xs1**0*xs2**0*xst**3&
     +force(9)  *xs1**0*xs2**2*xst**1&
     +force(10) *xs1**1*xs2**0*xst**2&
     +force(11) *xs1**1*xs2**2*xst**0&
     +force(12) *xs1**2*xs2**0*xst**1&
     +force(13) *xs1**3*xs2**0*xst**0&
     +force(14) *xs1**0*xs2**0*xst**4&
     +force(15) *xs1**0*xs2**2*xst**2&
     +force(16) *xs1**0*xs2**4*xst**0&
     +force(17) *xs1**1*xs2**0*xst**3&
     +force(18) *xs1**1*xs2**2*xst**1&
     +force(19) *xs1**2*xs2**0*xst**2&
     +force(20) *xs1**2*xs2**2*xst**0&
     +force(21) *xs1**3*xs2**0*xst**1&
     +force(22) *xs1**4*xs2**0*xst**0&
     +force(23) *xs1**0*xs2**0*xst**5&
     +force(24) *xs1**0*xs2**2*xst**3&
     +force(25) *xs1**0*xs2**4*xst**1&
     +force(26) *xs1**1*xs2**0*xst**4&
     +force(27) *xs1**1*xs2**2*xst**2&
     +force(28) *xs1**1*xs2**4*xst**0&
     +force(29) *xs1**2*xs2**0*xst**3&
     +force(30) *xs1**2*xs2**2*xst**1&
     +force(31) *xs1**3*xs2**0*xst**2&
     +force(32) *xs1**3*xs2**2*xst**0&
     +force(33) *xs1**4*xs2**0*xst**1&
     +force(34) *xs1**5*xs2**0*xst**0&
     +force(35) *xs1**0*xs2**0*xst**6&
     +force(36) *xs1**0*xs2**2*xst**4&
     +force(37) *xs1**0*xs2**4*xst**2&
     +force(38) *xs1**0*xs2**6*xst**0&
     +force(39) *xs1**1*xs2**0*xst**5&
     +force(40) *xs1**1*xs2**2*xst**3&
     +force(41) *xs1**1*xs2**4*xst**1&
     +force(42) *xs1**2*xs2**0*xst**4&
     +force(43) *xs1**2*xs2**2*xst**2&
     +force(44) *xs1**2*xs2**4*xst**0&
     +force(45) *xs1**3*xs2**0*xst**3&
     +force(46) *xs1**3*xs2**2*xst**1&
     +force(47) *xs1**4*xs2**0*xst**2&
     +force(48) *xs1**4*xs2**2*xst**0&
     +force(49) *xs1**5*xs2**0*xst**1&
     +force(50) *xs1**6*xs2**0*xst**0&
     +force(51) *xs1**0*xs2**0*xst**7&
     +force(52) *xs1**0*xs2**2*xst**5&
     +force(53) *xs1**0*xs2**4*xst**3&
     +force(54) *xs1**0*xs2**6*xst**1&
     +force(55) *xs1**1*xs2**0*xst**6&
     +force(56) *xs1**1*xs2**2*xst**4&
     +force(57) *xs1**1*xs2**4*xst**2&
     +force(58) *xs1**1*xs2**6*xst**0&
     +force(59) *xs1**2*xs2**0*xst**5&
     +force(60) *xs1**2*xs2**2*xst**3&
     +force(61) *xs1**2*xs2**4*xst**1&
     +force(62) *xs1**3*xs2**0*xst**4&
     +force(63) *xs1**3*xs2**2*xst**2&
     +force(64) *xs1**3*xs2**4*xst**0&
     +force(65) *xs1**4*xs2**0*xst**3&
     +force(66) *xs1**4*xs2**2*xst**1&
     +force(67) *xs1**5*xs2**0*xst**2&
     +force(68) *xs1**5*xs2**2*xst**0&
     +force(69) *xs1**6*xs2**0*xst**1&
     +force(70) *xs1**7*xs2**0*xst**0&
     +force(71) *xs1**0*xs2**0*xst**8&
     +force(72) *xs1**0*xs2**2*xst**6&
     +force(73) *xs1**0*xs2**4*xst**4&
     +force(74) *xs1**0*xs2**6*xst**2&
     +force(75) *xs1**0*xs2**8*xst**0&
     +force(76) *xs1**1*xs2**0*xst**7&
     +force(77) *xs1**1*xs2**2*xst**5&
     +force(78) *xs1**1*xs2**4*xst**3&
     +force(79) *xs1**1*xs2**6*xst**1&
     +force(80) *xs1**2*xs2**0*xst**6&
     +force(81) *xs1**2*xs2**2*xst**4&
     +force(82) *xs1**2*xs2**4*xst**2&
     +force(83) *xs1**2*xs2**6*xst**0&
     +force(84) *xs1**3*xs2**0*xst**5&
     +force(85) *xs1**3*xs2**2*xst**3&
     +force(86) *xs1**3*xs2**4*xst**1&
     +force(87) *xs1**4*xs2**0*xst**4&
     +force(88) *xs1**4*xs2**2*xst**2&
     +force(89) *xs1**4*xs2**4*xst**0&
     +force(90) *xs1**5*xs2**0*xst**3&
     +force(91) *xs1**5*xs2**2*xst**1&
     +force(92) *xs1**6*xs2**0*xst**2
 vp2= force(93) *xs1**6*xs2**2*xst**0&
     +force(94) *xs1**7*xs2**0*xst**1&
     +force(95) *xs1**8*xs2**0*xst**0&
     +force(96) *xs1**0*xs2**0*xst**9&
     +force(97) *xs1**0*xs2**2*xst**7&
     +force(98) *xs1**0*xs2**4*xst**5&
     +force(99) *xs1**0*xs2**6*xst**3&
     +force(100)*xs1**0*xs2**8*xst**1&
     +force(101)*xs1**1*xs2**0*xst**8&
     +force(102)*xs1**1*xs2**2*xst**6&
     +force(103)*xs1**1*xs2**4*xst**4&
     +force(104)*xs1**1*xs2**6*xst**2&
     +force(105)*xs1**1*xs2**8*xst**0&
     +force(106)*xs1**2*xs2**0*xst**7&
     +force(107)*xs1**2*xs2**2*xst**5&
     +force(108)*xs1**2*xs2**4*xst**3&
     +force(109)*xs1**2*xs2**6*xst**1&
     +force(110)*xs1**3*xs2**0*xst**6&
     +force(111)*xs1**3*xs2**2*xst**4&
     +force(112)*xs1**3*xs2**4*xst**2&
     +force(113)*xs1**3*xs2**6*xst**0&
     +force(114)*xs1**4*xs2**0*xst**5&
     +force(115)*xs1**4*xs2**2*xst**3&
     +force(116)*xs1**4*xs2**4*xst**1&
     +force(117)*xs1**5*xs2**0*xst**4&
     +force(118)*xs1**5*xs2**2*xst**2&
     +force(119)*xs1**5*xs2**4*xst**0&
     +force(120)*xs1**6*xs2**0*xst**3&
     +force(121)*xs1**6*xs2**2*xst**1&
     +force(122)*xs1**7*xs2**0*xst**2&
     +force(123)*xs1**7*xs2**2*xst**0&
     +force(124)*xs1**8*xs2**0*xst**1&
     +force(125)*xs1**9*xs2**0*xst**0&
     +force(126)*xs1**0*xs2**0*xst**10&
     +force(127)*xs1**0*xs2**2*xst**8&
     +force(128)*xs1**0*xs2**4*xst**6&
     +force(129)*xs1**0*xs2**6*xst**4&
     +force(130)*xs1**0*xs2**8*xst**2&
     +force(131)*xs1**0*xs2**10*xst**0&
     +force(132)*xs1**1*xs2**0*xst**9&
     +force(133)*xs1**1*xs2**2*xst**7&
     +force(134)*xs1**1*xs2**4*xst**5&
     +force(135)*xs1**1*xs2**6*xst**3&
     +force(136)*xs1**1*xs2**8*xst**1&
     +force(137)*xs1**2*xs2**0*xst**8&
     +force(138)*xs1**2*xs2**2*xst**6&
     +force(139)*xs1**2*xs2**4*xst**4&
     +force(140)*xs1**2*xs2**6*xst**2&
     +force(141)*xs1**2*xs2**8*xst**0&
     +force(142)*xs1**3*xs2**0*xst**7&
     +force(143)*xs1**3*xs2**2*xst**5&
     +force(144)*xs1**3*xs2**4*xst**3&
     +force(145)*xs1**3*xs2**6*xst**1&
     +force(146)*xs1**4*xs2**0*xst**6&
     +force(147)*xs1**4*xs2**2*xst**4&
     +force(148)*xs1**4*xs2**4*xst**2&
     +force(149)*xs1**4*xs2**6*xst**0&
     +force(150)*xs1**5*xs2**0*xst**5&
     +force(151)*xs1**5*xs2**2*xst**3&
     +force(152)*xs1**5*xs2**4*xst**1&
     +force(153)*xs1**6*xs2**0*xst**4&
     +force(154)*xs1**6*xs2**2*xst**2&
     +force(155)*xs1**6*xs2**4*xst**0&
     +force(156)*xs1**7*xs2**0*xst**3&
     +force(157)*xs1**7*xs2**2*xst**1&
     +force(158)*xs1**8*xs2**0*xst**2&
     +force(159)*xs1**8*xs2**2*xst**0&
     +force(160)*xs1**9*xs2**0*xst**1&
     +force(161)*xs1**10*xs2**0*xst**0&
     +force(162)*xs1**0*xs2**0*xst**11&
     +force(163)*xs1**0*xs2**2*xst**9&
     +force(164)*xs1**0*xs2**4*xst**7&
     +force(165)*xs1**0*xs2**6*xst**5&
     +force(166)*xs1**0*xs2**8*xst**3&
     +force(167)*xs1**0*xs2**10*xst**1&
     +force(168)*xs1**1*xs2**0*xst**10&
     +force(169)*xs1**1*xs2**2*xst**8&
     +force(170)*xs1**1*xs2**4*xst**6&
     +force(171)*xs1**1*xs2**6*xst**4&
     +force(172)*xs1**1*xs2**8*xst**2&
     +force(173)*xs1**1*xs2**10*xst**0&
     +force(174)*xs1**2*xs2**0*xst**9
 vp3= force(175)*xs1**2*xs2**2*xst**7&
     +force(176)*xs1**2*xs2**4*xst**5&
     +force(177)*xs1**2*xs2**6*xst**3&
     +force(178)*xs1**2*xs2**8*xst**1&
     +force(179)*xs1**3*xs2**0*xst**8&
     +force(180)*xs1**3*xs2**2*xst**6&
     +force(181)*xs1**3*xs2**4*xst**4&
     +force(182)*xs1**3*xs2**6*xst**2&
     +force(183)*xs1**3*xs2**8*xst**0&
     +force(184)*xs1**4*xs2**0*xst**7&
     +force(185)*xs1**4*xs2**2*xst**5&
     +force(186)*xs1**4*xs2**4*xst**3&
     +force(187)*xs1**4*xs2**6*xst**1&
     +force(188)*xs1**5*xs2**0*xst**6&
     +force(189)*xs1**5*xs2**2*xst**4&
     +force(190)*xs1**5*xs2**4*xst**2&
     +force(191)*xs1**5*xs2**6*xst**0&
     +force(192)*xs1**6*xs2**0*xst**5&
     +force(193)*xs1**6*xs2**2*xst**3&
     +force(194)*xs1**6*xs2**4*xst**1&
     +force(195)*xs1**7*xs2**0*xst**4&
     +force(196)*xs1**7*xs2**2*xst**2&
     +force(197)*xs1**7*xs2**4*xst**0&
     +force(198)*xs1**8*xs2**0*xst**3&
     +force(199)*xs1**8*xs2**2*xst**1&
     +force(200)*xs1**9*xs2**0*xst**2&
     +force(201)*xs1**9*xs2**2*xst**0&
     +force(202)*xs1**0*xs2**0*xst**12&
     +force(203)*xs1**0*xs2**2*xst**10&
     +force(204)*xs1**0*xs2**4*xst**8&
     +force(205)*xs1**0*xs2**6*xst**6&
     +force(206)*xs1**0*xs2**8*xst**4&
     +force(207)*xs1**0*xs2**10*xst**2&
     +force(208)*xs1**0*xs2**12*xst**0&
     +force(209)*xs1**1*xs2**0*xst**11&    
     +force(210)*xs1**1*xs2**2*xst**9&
     +force(211)*xs1**1*xs2**4*xst**7&
     +force(212)*xs1**1*xs2**6*xst**5&
     +force(213)*xs1**1*xs2**8*xst**3&
     +force(214)*xs1**1*xs2**10*xst**1&
     +force(215)*xs1**2*xs2**0*xst**10&
     +force(216)*xs1**2*xs2**2*xst**8&
     +force(217)*xs1**2*xs2**4*xst**6&
     +force(218)*xs1**2*xs2**6*xst**4&
     +force(219)*xs1**2*xs2**8*xst**2&
     +force(220)*xs1**2*xs2**10*xst**0&
     +force(221)*xs1**3*xs2**0*xst**9&
     +force(222)*xs1**3*xs2**2*xst**7&
     +force(223)*xs1**3*xs2**4*xst**5&
     +force(224)*xs1**3*xs2**6*xst**3&
     +force(225)*xs1**3*xs2**8*xst**1&
     +force(226)*xs1**4*xs2**0*xst**8&
     +force(227)*xs1**4*xs2**2*xst**6&
     +force(228)*xs1**4*xs2**4*xst**4&
     +force(229)*xs1**4*xs2**6*xst**2&
     +force(230)*xs1**4*xs2**8*xst**0&
     +force(231)*xs1**5*xs2**0*xst**7&
     +force(232)*xs1**5*xs2**2*xst**5&
     +force(233)*xs1**5*xs2**4*xst**3&
     +force(234)*xs1**5*xs2**6*xst**1&
     +force(235)*xs1**6*xs2**0*xst**6&
     +force(236)*xs1**6*xs2**2*xst**4&
     +force(237)*xs1**6*xs2**4*xst**2&
     +force(238)*xs1**6*xs2**6*xst**0&
     +force(239)*xs1**7*xs2**0*xst**5&
     +force(240)*xs1**7*xs2**2*xst**3

       vp=vp1+vp2+vp3
       !
       vps1=42395.535333_ark*xep1
       vps2=42395.535333_ark*xep2

       y1=1._ark/(1._ark+exp(ut*(x0-r1)))
       y2=1._ark/(1._ark+exp(ut*(x0-r2)))
       y12=1._ark/(1._ark+exp(ut2*(x02-r1)))
       y22=1._ark/(1._ark+exp(ut2*(x02-r2)))

       vp=vp*xep3*(1-y12)*(1-y22)
       voh1=vpb1*(1-y1)+y1*vps1
       voh2=vpb2*(1-y2)+y2*vps2


       f=v0+vp+voh1+voh2+vhh

      end function MLpoten_xy2_bubukina







 !===============================================================================
 !                   electric dipole moment section
 !===============================================================================

 !returns electric dipole moment cartesian coordinates in the user-defined frame for locals specified
 recursive subroutine MLdms2pqr_xy2(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: i,ik(1:3)
    real(ark)             :: b(molec%ncoords), mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,dipp,dipq
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
    !
    b(1) = r1 - molec%req(1)
    b(2) = r2 - molec%req(2)
    b(3) = cos(alpha) - cos(molec%alphaeq(1))
    !
    dipp = 0
    !
    do i=1,extF%nterms(1)
       !
       ik(:) = extF%term(:,i,1)
       !
       dipp = dipp + b(1)**ik(1)*b(2)**ik(2)*b(3)**ik(3)*extF%coef(i,1)
       !
       if (ik(1)/=ik(2)) then
         !
         dipp = dipp - b(1)**ik(2)*b(2)**ik(1)*b(3)**ik(3)*extF%coef(i,1)
         !
       endif
       !
    enddo
    !
    mu(2) = dipp
    !
    dipq = 0 
    !
    do i=1,extF%nterms(2)
       !
       ik(:) = extF%term(:,i,2)
       !
       dipq = dipq + b(1)**ik(1)*b(2)**ik(2)*b(3)**ik(3)*extF%coef(i,2)
       !
       if (ik(1)/=ik(2)) then
         !
         dipq = dipq + b(1)**ik(2)*b(2)**ik(1)*b(3)**ik(3)*extF%coef(i,2)
         !
       endif
       !
    enddo
    !
    mu(1) = dipq * sin(alpha)
    !
    mu(3) = 0
    !
    f(1:3) = matmul(mu,tmat)
    !
 end subroutine MLdms2pqr_xy2
 !

 !returns electric dipole moment cartesian coordinates in the user-defined frame for locals specified
 !
 recursive subroutine MLdms2pqr_xy2_coeff(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: k
    real(ark)             :: y1,y2,y3, mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,re,ae,b0
    real(ark)             :: p(3:extF%nterms(1)),q(5:extF%nterms(2)),v0,v1,v2,v3,v4,v5,v6,v7,v8,xyz0(natoms,3)
    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      select case(trim(molec%coords_transform))
      case default
         write (out,"('MLdms2pqr_xy2_coeff: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLdms2pqr_xy2_coeff - bad coord. type'
      case('R-RHO-Z','R-RHO-Z-ECKART')
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
    re = extF%coef(1,1)
    ae = extF%coef(2,1)*pi/180.0_ark
    b0 = extF%coef(3,1)
    !
    y1 = (r1 - re)*exp(-b0*(r1-re)**2)
    y2 = (r2 - re)*exp(-b0*(r1-re)**2)
    y3 = cos(alpha) - cos(ae)
    !
    k = extF%nterms(1)
    !
    p(3:k) = extF%coef(3:k,1)
    !
    v1 = p( 3)*y1**1*y2**0*y3**0& 
       - p( 3)*y1**0*y2**1*y3**0
    v2 = p( 4)*y1**1*y2**0*y3**1& 
       - p( 4)*y1**0*y2**1*y3**1&  
       + p( 5)*y1**2*y2**0*y3**0& 
       - p( 5)*y1**0*y2**2*y3**0
       !
    v3 = p( 6)*y1**1*y2**0*y3**2& 
       - p( 6)*y1**0*y2**1*y3**2& 
       + p( 7)*y1**2*y2**0*y3**1& 
       - p( 7)*y1**0*y2**2*y3**1& 
       + p( 8)*y1**2*y2**1*y3**0& 
       - p( 8)*y1**1*y2**2*y3**0& 
       + p( 9)*y1**3*y2**0*y3**0& 
       - p( 9)*y1**0*y2**3*y3**0
       !
     v4 =p(10)*y1**1*y2**0*y3**3& 
       - p(10)*y1**0*y2**1*y3**3& 
       + p(11)*y1**2*y2**0*y3**2& 
       - p(11)*y1**0*y2**2*y3**2& 
       + p(12)*y1**2*y2**1*y3**1& 
       - p(12)*y1**1*y2**2*y3**1& 
       + p(13)*y1**3*y2**0*y3**1& 
       - p(13)*y1**0*y2**3*y3**1& 
       + p(14)*y1**3*y2**1*y3**0& 
       - p(14)*y1**1*y2**3*y3**0& 
       + p(15)*y1**4*y2**0*y3**0& 
       - p(15)*y1**0*y2**4*y3**0
       !
     v5 =p(16)*y1**1*y2**0*y3**4& 
       - p(16)*y1**0*y2**1*y3**4& 
       + p(17)*y1**2*y2**0*y3**3& 
       - p(17)*y1**0*y2**2*y3**3& 
       + p(18)*y1**2*y2**1*y3**2& 
       - p(18)*y1**1*y2**2*y3**2& 
       + p(19)*y1**3*y2**0*y3**2& 
       - p(19)*y1**0*y2**3*y3**2& 
       + p(20)*y1**3*y2**1*y3**1& 
       - p(20)*y1**1*y2**3*y3**1& 
       + p(21)*y1**3*y2**2*y3**0& 
       - p(21)*y1**2*y2**3*y3**0& 
       + p(22)*y1**4*y2**0*y3**1& 
       - p(22)*y1**0*y2**4*y3**1& 
       + p(23)*y1**4*y2**1*y3**0& 
       - p(23)*y1**1*y2**4*y3**0& 
       + p(24)*y1**5*y2**0*y3**0& 
       - p(24)*y1**0*y2**5*y3**0
       !
     v6 =p(25)*y1**1*y2**0*y3**5& 
       - p(25)*y1**0*y2**1*y3**5& 
       + p(26)*y1**2*y2**0*y3**4& 
       - p(26)*y1**0*y2**2*y3**4& 
       + p(27)*y1**2*y2**1*y3**3& 
       - p(27)*y1**1*y2**2*y3**3& 
       + p(28)*y1**3*y2**0*y3**3& 
       - p(28)*y1**0*y2**3*y3**3& 
       + p(29)*y1**3*y2**1*y3**2& 
       - p(29)*y1**1*y2**3*y3**2& 
       + p(30)*y1**3*y2**2*y3**1& 
       - p(30)*y1**2*y2**3*y3**1& 
       + p(31)*y1**4*y2**0*y3**2& 
       - p(31)*y1**0*y2**4*y3**2& 
       + p(32)*y1**4*y2**1*y3**1& 
       - p(32)*y1**1*y2**4*y3**1& 
       + p(33)*y1**4*y2**2*y3**0& 
       - p(33)*y1**2*y2**4*y3**0& 
       + p(34)*y1**5*y2**0*y3**1& 
       - p(34)*y1**0*y2**5*y3**1& 
       + p(35)*y1**5*y2**1*y3**0& 
       - p(35)*y1**1*y2**5*y3**0& 
       + p(36)*y1**6*y2**0*y3**0& 
       - p(36)*y1**0*y2**6*y3**0
       !
    v7 = p(37)*y1**1*y2**0*y3**6& 
       - p(37)*y1**0*y2**1*y3**6& 
       + p(38)*y1**2*y2**0*y3**5& 
       - p(38)*y1**0*y2**2*y3**5& 
       + p(39)*y1**2*y2**1*y3**4& 
       - p(39)*y1**1*y2**2*y3**4& 
       + p(40)*y1**3*y2**0*y3**4& 
       - p(40)*y1**0*y2**3*y3**4& 
       + p(41)*y1**3*y2**1*y3**3& 
       - p(41)*y1**1*y2**3*y3**3& 
       + p(42)*y1**3*y2**2*y3**2& 
       - p(42)*y1**2*y2**3*y3**2& 
       + p(43)*y1**4*y2**0*y3**3& 
       - p(43)*y1**0*y2**4*y3**3& 
       + p(44)*y1**4*y2**1*y3**2& 
       - p(44)*y1**1*y2**4*y3**2& 
       + p(45)*y1**4*y2**2*y3**1& 
       - p(45)*y1**2*y2**4*y3**1& 
       + p(46)*y1**4*y2**3*y3**0& 
       - p(46)*y1**3*y2**4*y3**0& 
       + p(47)*y1**5*y2**0*y3**2& 
       - p(47)*y1**0*y2**5*y3**2& 
       + p(48)*y1**5*y2**1*y3**1& 
       - p(48)*y1**1*y2**5*y3**1& 
       + p(49)*y1**5*y2**2*y3**0& 
       - p(49)*y1**2*y2**5*y3**0& 
       + p(50)*y1**6*y2**0*y3**1& 
       - p(50)*y1**0*y2**6*y3**1& 
       + p(51)*y1**6*y2**1*y3**0& 
       - p(51)*y1**1*y2**6*y3**0& 
       + p(52)*y1**7*y2**0*y3**0& 
       - p(52)*y1**0*y2**7*y3**0
       !
    v8 = p(53)*y1**1*y2**0*y3**7& 
       - p(53)*y1**0*y2**1*y3**7& 
       + p(54)*y1**2*y2**0*y3**6& 
       - p(54)*y1**0*y2**2*y3**6& 
       + p(55)*y1**2*y2**1*y3**5& 
       - p(55)*y1**1*y2**2*y3**5& 
       + p(56)*y1**3*y2**0*y3**5& 
       - p(56)*y1**0*y2**3*y3**5& 
       + p(57)*y1**3*y2**1*y3**4& 
       - p(57)*y1**1*y2**3*y3**4& 
       + p(58)*y1**3*y2**2*y3**3& 
       - p(58)*y1**2*y2**3*y3**3& 
       + p(59)*y1**4*y2**0*y3**4& 
       - p(59)*y1**0*y2**4*y3**4& 
       + p(60)*y1**4*y2**1*y3**3& 
       - p(60)*y1**1*y2**4*y3**3& 
       + p(61)*y1**4*y2**2*y3**2& 
       - p(61)*y1**2*y2**4*y3**2& 
       + p(62)*y1**4*y2**3*y3**1& 
       - p(62)*y1**3*y2**4*y3**1& 
       + p(63)*y1**5*y2**0*y3**3& 
       - p(63)*y1**0*y2**5*y3**3& 
       + p(64)*y1**5*y2**1*y3**2& 
       - p(64)*y1**1*y2**5*y3**2& 
       + p(65)*y1**5*y2**2*y3**1& 
       - p(65)*y1**2*y2**5*y3**1& 
       + p(66)*y1**5*y2**3*y3**0& 
       - p(66)*y1**3*y2**5*y3**0& 
       + p(67)*y1**6*y2**0*y3**2& 
       - p(67)*y1**0*y2**6*y3**2& 
       + p(68)*y1**6*y2**1*y3**1& 
       - p(68)*y1**1*y2**6*y3**1& 
       + p(69)*y1**6*y2**2*y3**0& 
       - p(69)*y1**2*y2**6*y3**0& 
       + p(70)*y1**7*y2**0*y3**1& 
       - p(70)*y1**0*y2**7*y3**1& 
       + p(71)*y1**7*y2**1*y3**0& 
       - p(71)*y1**1*y2**7*y3**0& 
       + p(72)*y1**8*y2**0*y3**0& 
       - p(72)*y1**0*y2**8*y3**0
    !
    mu(2) = v1+v2+v3+v4+v5+v6+v7+v8
    !
    re = extF%coef(1,2)
    ae = extF%coef(2,2)*pi/180.0_ark
    b0 = extF%coef(3,2)
    !
    y1 = (r1 - re)*exp(-b0*(r1-re)**2)
    y2 = (r2 - re)*exp(-b0*(r1-re)**2)
    y3 = cos(alpha) - cos(ae)
    !
    k = extF%nterms(2) 
    !
    q(5:k) = extF%coef(5:k,2)
    !
    v0 = q(5)*y1**0*y2**0*y3**0
    v1 = q(6)*y1**0*y2**0*y3**1& 
       + q(7)*y1**1*y2**0*y3**0& 
       + q(7)*y1**0*y2**1*y3**0
    v2 = q(8)*y1**0*y2**0*y3**2& 
       + q(9)*y1**1*y2**0*y3**1& 
       + q(9)*y1**0*y2**1*y3**1& 
       + q(10)*y1**1*y2**1*y3**0& 
       + q(11)*y1**2*y2**0*y3**0& 
       + q(11)*y1**0*y2**2*y3**0
    v3 = q(12)*y1**0*y2**0*y3**3& 
       + q(13)*y1**1*y2**0*y3**2& 
       + q(13)*y1**0*y2**1*y3**2& 
       + q(14)*y1**1*y2**1*y3**1& 
       + q(15)*y1**2*y2**0*y3**1& 
       + q(15)*y1**0*y2**2*y3**1& 
       + q(16)*y1**2*y2**1*y3**0& 
       + q(16)*y1**1*y2**2*y3**0& 
       + q(17)*y1**3*y2**0*y3**0& 
       + q(17)*y1**0*y2**3*y3**0
       !
     v4 = q(18)*y1**0*y2**0*y3**4& 
       + q(19)*y1**1*y2**0*y3**3& 
       + q(19)*y1**0*y2**1*y3**3& 
       + q(20)*y1**1*y2**1*y3**2& 
       + q(21)*y1**2*y2**0*y3**2& 
       + q(21)*y1**0*y2**2*y3**2& 
       + q(22)*y1**2*y2**1*y3**1& 
       + q(22)*y1**1*y2**2*y3**1& 
       + q(23)*y1**2*y2**2*y3**0& 
       + q(24)*y1**3*y2**0*y3**1& 
       + q(24)*y1**0*y2**3*y3**1& 
       + q(25)*y1**3*y2**1*y3**0& 
       + q(25)*y1**1*y2**3*y3**0& 
       + q(26)*y1**4*y2**0*y3**0& 
       + q(26)*y1**0*y2**4*y3**0
       !
     v5 = q(27)*y1**0*y2**0*y3**5& 
       + q(28)*y1**1*y2**0*y3**4& 
       + q(28)*y1**0*y2**1*y3**4& 
       + q(29)*y1**1*y2**1*y3**3& 
       + q(30)*y1**2*y2**0*y3**3& 
       + q(30)*y1**0*y2**2*y3**3& 
       + q(31)*y1**2*y2**1*y3**2& 
       + q(31)*y1**1*y2**2*y3**2& 
       + q(32)*y1**2*y2**2*y3**1& 
       + q(33)*y1**3*y2**0*y3**2& 
       + q(33)*y1**0*y2**3*y3**2& 
       + q(34)*y1**3*y2**1*y3**1& 
       + q(34)*y1**1*y2**3*y3**1& 
       + q(35)*y1**3*y2**2*y3**0& 
       + q(35)*y1**2*y2**3*y3**0& 
       + q(36)*y1**4*y2**0*y3**1& 
       + q(36)*y1**0*y2**4*y3**1& 
       + q(37)*y1**4*y2**1*y3**0& 
       + q(37)*y1**1*y2**4*y3**0& 
       + q(38)*y1**5*y2**0*y3**0& 
       + q(38)*y1**0*y2**5*y3**0
       !
    v6 = q(39)*y1**0*y2**0*y3**6& 
       + q(40)*y1**1*y2**0*y3**5& 
       + q(40)*y1**0*y2**1*y3**5& 
       + q(41)*y1**1*y2**1*y3**4& 
       + q(42)*y1**2*y2**0*y3**4& 
       + q(42)*y1**0*y2**2*y3**4& 
       + q(43)*y1**2*y2**1*y3**3& 
       + q(43)*y1**1*y2**2*y3**3& 
       + q(44)*y1**2*y2**2*y3**2& 
       + q(45)*y1**3*y2**0*y3**3& 
       + q(45)*y1**0*y2**3*y3**3& 
       + q(46)*y1**3*y2**1*y3**2& 
       + q(46)*y1**1*y2**3*y3**2& 
       + q(47)*y1**3*y2**2*y3**1& 
       + q(47)*y1**2*y2**3*y3**1& 
       + q(48)*y1**3*y2**3*y3**0& 
       + q(49)*y1**4*y2**0*y3**2& 
       + q(49)*y1**0*y2**4*y3**2& 
       + q(50)*y1**4*y2**1*y3**1& 
       + q(50)*y1**1*y2**4*y3**1& 
       + q(51)*y1**4*y2**2*y3**0& 
       + q(51)*y1**2*y2**4*y3**0& 
       + q(52)*y1**5*y2**0*y3**1& 
       + q(52)*y1**0*y2**5*y3**1& 
       + q(53)*y1**5*y2**1*y3**0& 
       + q(53)*y1**1*y2**5*y3**0& 
       + q(54)*y1**6*y2**0*y3**0& 
       + q(54)*y1**0*y2**6*y3**0
       !
    v7 = q(55)*y1**0*y2**0*y3**7& 
       + q(56)*y1**1*y2**0*y3**6& 
       + q(56)*y1**0*y2**1*y3**6& 
       + q(57)*y1**1*y2**1*y3**5& 
       + q(58)*y1**2*y2**0*y3**5& 
       + q(58)*y1**0*y2**2*y3**5& 
       + q(59)*y1**2*y2**1*y3**4& 
       + q(59)*y1**1*y2**2*y3**4& 
       + q(60)*y1**2*y2**2*y3**3& 
       + q(61)*y1**3*y2**0*y3**4& 
       + q(61)*y1**0*y2**3*y3**4& 
       + q(62)*y1**3*y2**1*y3**3& 
       + q(62)*y1**1*y2**3*y3**3& 
       + q(63)*y1**3*y2**2*y3**2& 
       + q(63)*y1**2*y2**3*y3**2& 
       + q(64)*y1**3*y2**3*y3**1& 
       + q(65)*y1**4*y2**0*y3**3& 
       + q(65)*y1**0*y2**4*y3**3& 
       + q(66)*y1**4*y2**1*y3**2& 
       + q(66)*y1**1*y2**4*y3**2& 
       + q(67)*y1**4*y2**2*y3**1& 
       + q(67)*y1**2*y2**4*y3**1& 
       + q(68)*y1**4*y2**3*y3**0& 
       + q(68)*y1**3*y2**4*y3**0& 
       + q(69)*y1**5*y2**0*y3**2& 
       + q(69)*y1**0*y2**5*y3**2& 
       + q(70)*y1**5*y2**1*y3**1& 
       + q(70)*y1**1*y2**5*y3**1& 
       + q(71)*y1**5*y2**2*y3**0& 
       + q(71)*y1**2*y2**5*y3**0& 
       + q(72)*y1**6*y2**0*y3**1& 
       + q(72)*y1**0*y2**6*y3**1& 
       + q(73)*y1**6*y2**1*y3**0& 
       + q(73)*y1**1*y2**6*y3**0& 
       + q(74)*y1**7*y2**0*y3**0& 
       + q(74)*y1**0*y2**7*y3**0
       !
    v8 = q(75)*y1**0*y2**0*y3**8& 
       + q(76)*y1**1*y2**0*y3**7& 
       + q(76)*y1**0*y2**1*y3**7& 
       + q(77)*y1**1*y2**1*y3**6& 
       + q(78)*y1**2*y2**0*y3**6& 
       + q(78)*y1**0*y2**2*y3**6& 
       + q(79)*y1**2*y2**1*y3**5& 
       + q(79)*y1**1*y2**2*y3**5& 
       + q(80)*y1**2*y2**2*y3**4& 
       + q(81)*y1**3*y2**0*y3**5& 
       + q(81)*y1**0*y2**3*y3**5& 
       + q(82)*y1**3*y2**1*y3**4& 
       + q(82)*y1**1*y2**3*y3**4& 
       + q(83)*y1**3*y2**2*y3**3& 
       + q(83)*y1**2*y2**3*y3**3& 
       + q(84)*y1**3*y2**3*y3**2& 
       + q(85)*y1**4*y2**0*y3**4& 
       + q(85)*y1**0*y2**4*y3**4& 
       + q(86)*y1**4*y2**1*y3**3& 
       + q(86)*y1**1*y2**4*y3**3& 
       + q(87)*y1**4*y2**2*y3**2& 
       + q(87)*y1**2*y2**4*y3**2& 
       + q(88)*y1**4*y2**3*y3**1& 
       + q(88)*y1**3*y2**4*y3**1& 
       + q(89)*y1**4*y2**4*y3**0& 
       + q(90)*y1**5*y2**0*y3**3& 
       + q(90)*y1**0*y2**5*y3**3& 
       + q(91)*y1**5*y2**1*y3**2& 
       + q(91)*y1**1*y2**5*y3**2& 
       + q(92)*y1**5*y2**2*y3**1& 
       + q(92)*y1**2*y2**5*y3**1& 
       + q(93)*y1**5*y2**3*y3**0& 
       + q(93)*y1**3*y2**5*y3**0& 
       + q(94)*y1**6*y2**0*y3**2& 
       + q(94)*y1**0*y2**6*y3**2& 
       + q(95)*y1**6*y2**1*y3**1& 
       + q(95)*y1**1*y2**6*y3**1& 
       + q(96)*y1**6*y2**2*y3**0& 
       + q(96)*y1**2*y2**6*y3**0& 
       + q(97)*y1**7*y2**0*y3**1& 
       + q(97)*y1**0*y2**7*y3**1& 
       + q(98)*y1**7*y2**1*y3**0& 
       + q(98)*y1**1*y2**7*y3**0& 
       + q(99)*y1**8*y2**0*y3**0& 
       + q(99)*y1**0*y2**8*y3**0
       !
    mu(1) =(v0+v1+v2+v3+v4+v5+v6+v7+v8)*sin(alpha)
    !
    mu(3) = 0
    !
    f(1:3) = matmul(mu,tmat)
    !
 end subroutine MLdms2pqr_xy2_coeff




 !returns electric dipole moment cartesian coordinates in the user-defined frame for locals specified
 !
 recursive subroutine MLdms2pqr_xy2_linear(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: k
    real(ark)             :: y1,y2,y3, mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,re,ae,b0
    real(ark)             :: p(3:extF%nterms(1)),q(5:extF%nterms(2)),v0,v1,v2,v3,v4,v5,v6,v7,v8,xyz0(natoms,3)
    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      select case(trim(molec%coords_transform))
      case default
         write (out,"('MLdms2pqr_xy2_linear: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLdms2pqr_xy2_linear - bad coord. type'
      case('R-RHO-Z','R-RHO-Z-ECKART')
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
    re = extF%coef(1,1)
    ae = extF%coef(2,1)*pi/180.0_ark
    b0 = extF%coef(3,1)
    !
    y1 = (r1 - re)*exp(-b0*(r1-re)**2)
    y2 = (r2 - re)*exp(-b0*(r1-re)**2)
    !
    y3 = alpha-ae
    !
    k = extF%nterms(1)
    !
    p(3:k) = extF%coef(3:k,1)
    !
    v1 = p( 3)*y1**1*y2**0*y3**0& 
       - p( 3)*y1**0*y2**1*y3**0
    v2 = p( 4)*y1**1*y2**0*y3**1& 
       - p( 4)*y1**0*y2**1*y3**1&  
       + p( 5)*y1**2*y2**0*y3**0& 
       - p( 5)*y1**0*y2**2*y3**0
       !
    v3 = p( 6)*y1**1*y2**0*y3**2& 
       - p( 6)*y1**0*y2**1*y3**2& 
       + p( 7)*y1**2*y2**0*y3**1& 
       - p( 7)*y1**0*y2**2*y3**1& 
       + p( 8)*y1**2*y2**1*y3**0& 
       - p( 8)*y1**1*y2**2*y3**0& 
       + p( 9)*y1**3*y2**0*y3**0& 
       - p( 9)*y1**0*y2**3*y3**0
       !
     v4 =p(10)*y1**1*y2**0*y3**3& 
       - p(10)*y1**0*y2**1*y3**3& 
       + p(11)*y1**2*y2**0*y3**2& 
       - p(11)*y1**0*y2**2*y3**2& 
       + p(12)*y1**2*y2**1*y3**1& 
       - p(12)*y1**1*y2**2*y3**1& 
       + p(13)*y1**3*y2**0*y3**1& 
       - p(13)*y1**0*y2**3*y3**1& 
       + p(14)*y1**3*y2**1*y3**0& 
       - p(14)*y1**1*y2**3*y3**0& 
       + p(15)*y1**4*y2**0*y3**0& 
       - p(15)*y1**0*y2**4*y3**0
       !
     v5 =p(16)*y1**1*y2**0*y3**4& 
       - p(16)*y1**0*y2**1*y3**4& 
       + p(17)*y1**2*y2**0*y3**3& 
       - p(17)*y1**0*y2**2*y3**3& 
       + p(18)*y1**2*y2**1*y3**2& 
       - p(18)*y1**1*y2**2*y3**2& 
       + p(19)*y1**3*y2**0*y3**2& 
       - p(19)*y1**0*y2**3*y3**2& 
       + p(20)*y1**3*y2**1*y3**1& 
       - p(20)*y1**1*y2**3*y3**1& 
       + p(21)*y1**3*y2**2*y3**0& 
       - p(21)*y1**2*y2**3*y3**0& 
       + p(22)*y1**4*y2**0*y3**1& 
       - p(22)*y1**0*y2**4*y3**1& 
       + p(23)*y1**4*y2**1*y3**0& 
       - p(23)*y1**1*y2**4*y3**0& 
       + p(24)*y1**5*y2**0*y3**0& 
       - p(24)*y1**0*y2**5*y3**0
       !
     v6 =p(25)*y1**1*y2**0*y3**5& 
       - p(25)*y1**0*y2**1*y3**5& 
       + p(26)*y1**2*y2**0*y3**4& 
       - p(26)*y1**0*y2**2*y3**4& 
       + p(27)*y1**2*y2**1*y3**3& 
       - p(27)*y1**1*y2**2*y3**3& 
       + p(28)*y1**3*y2**0*y3**3& 
       - p(28)*y1**0*y2**3*y3**3& 
       + p(29)*y1**3*y2**1*y3**2& 
       - p(29)*y1**1*y2**3*y3**2& 
       + p(30)*y1**3*y2**2*y3**1& 
       - p(30)*y1**2*y2**3*y3**1& 
       + p(31)*y1**4*y2**0*y3**2& 
       - p(31)*y1**0*y2**4*y3**2& 
       + p(32)*y1**4*y2**1*y3**1& 
       - p(32)*y1**1*y2**4*y3**1& 
       + p(33)*y1**4*y2**2*y3**0& 
       - p(33)*y1**2*y2**4*y3**0& 
       + p(34)*y1**5*y2**0*y3**1& 
       - p(34)*y1**0*y2**5*y3**1& 
       + p(35)*y1**5*y2**1*y3**0& 
       - p(35)*y1**1*y2**5*y3**0& 
       + p(36)*y1**6*y2**0*y3**0& 
       - p(36)*y1**0*y2**6*y3**0
       !
    v7 = p(37)*y1**1*y2**0*y3**6& 
       - p(37)*y1**0*y2**1*y3**6& 
       + p(38)*y1**2*y2**0*y3**5& 
       - p(38)*y1**0*y2**2*y3**5& 
       + p(39)*y1**2*y2**1*y3**4& 
       - p(39)*y1**1*y2**2*y3**4& 
       + p(40)*y1**3*y2**0*y3**4& 
       - p(40)*y1**0*y2**3*y3**4& 
       + p(41)*y1**3*y2**1*y3**3& 
       - p(41)*y1**1*y2**3*y3**3& 
       + p(42)*y1**3*y2**2*y3**2& 
       - p(42)*y1**2*y2**3*y3**2& 
       + p(43)*y1**4*y2**0*y3**3& 
       - p(43)*y1**0*y2**4*y3**3& 
       + p(44)*y1**4*y2**1*y3**2& 
       - p(44)*y1**1*y2**4*y3**2& 
       + p(45)*y1**4*y2**2*y3**1& 
       - p(45)*y1**2*y2**4*y3**1& 
       + p(46)*y1**4*y2**3*y3**0& 
       - p(46)*y1**3*y2**4*y3**0& 
       + p(47)*y1**5*y2**0*y3**2& 
       - p(47)*y1**0*y2**5*y3**2& 
       + p(48)*y1**5*y2**1*y3**1& 
       - p(48)*y1**1*y2**5*y3**1& 
       + p(49)*y1**5*y2**2*y3**0& 
       - p(49)*y1**2*y2**5*y3**0& 
       + p(50)*y1**6*y2**0*y3**1& 
       - p(50)*y1**0*y2**6*y3**1& 
       + p(51)*y1**6*y2**1*y3**0& 
       - p(51)*y1**1*y2**6*y3**0& 
       + p(52)*y1**7*y2**0*y3**0& 
       - p(52)*y1**0*y2**7*y3**0
       !
    v8 = p(53)*y1**1*y2**0*y3**7& 
       - p(53)*y1**0*y2**1*y3**7& 
       + p(54)*y1**2*y2**0*y3**6& 
       - p(54)*y1**0*y2**2*y3**6& 
       + p(55)*y1**2*y2**1*y3**5& 
       - p(55)*y1**1*y2**2*y3**5& 
       + p(56)*y1**3*y2**0*y3**5& 
       - p(56)*y1**0*y2**3*y3**5& 
       + p(57)*y1**3*y2**1*y3**4& 
       - p(57)*y1**1*y2**3*y3**4& 
       + p(58)*y1**3*y2**2*y3**3& 
       - p(58)*y1**2*y2**3*y3**3& 
       + p(59)*y1**4*y2**0*y3**4& 
       - p(59)*y1**0*y2**4*y3**4& 
       + p(60)*y1**4*y2**1*y3**3& 
       - p(60)*y1**1*y2**4*y3**3& 
       + p(61)*y1**4*y2**2*y3**2& 
       - p(61)*y1**2*y2**4*y3**2& 
       + p(62)*y1**4*y2**3*y3**1& 
       - p(62)*y1**3*y2**4*y3**1& 
       + p(63)*y1**5*y2**0*y3**3& 
       - p(63)*y1**0*y2**5*y3**3& 
       + p(64)*y1**5*y2**1*y3**2& 
       - p(64)*y1**1*y2**5*y3**2& 
       + p(65)*y1**5*y2**2*y3**1& 
       - p(65)*y1**2*y2**5*y3**1& 
       + p(66)*y1**5*y2**3*y3**0& 
       - p(66)*y1**3*y2**5*y3**0& 
       + p(67)*y1**6*y2**0*y3**2& 
       - p(67)*y1**0*y2**6*y3**2& 
       + p(68)*y1**6*y2**1*y3**1& 
       - p(68)*y1**1*y2**6*y3**1& 
       + p(69)*y1**6*y2**2*y3**0& 
       - p(69)*y1**2*y2**6*y3**0& 
       + p(70)*y1**7*y2**0*y3**1& 
       - p(70)*y1**0*y2**7*y3**1& 
       + p(71)*y1**7*y2**1*y3**0& 
       - p(71)*y1**1*y2**7*y3**0& 
       + p(72)*y1**8*y2**0*y3**0& 
       - p(72)*y1**0*y2**8*y3**0
    !
    mu(2) = v1+v2+v3+v4+v5+v6+v7+v8
    !
    re = extF%coef(1,2)
    ae = extF%coef(2,2)*pi/180.0_ark
    b0 = extF%coef(3,1)
    !
    y1 = (r1 - re)*exp(-b0*(r1-re)**2)
    y2 = (r2 - re)*exp(-b0*(r1-re)**2)
    !
    y3 = alpha-ae
    !
    k = extF%nterms(2) 
    !
    q(5:k) = extF%coef(5:k,2)
    !
    v0 = q(5)*y1**0*y2**0*y3**0
    v1 = q(6)*y1**0*y2**0*y3**1& 
       + q(7)*y1**1*y2**0*y3**0& 
       + q(7)*y1**0*y2**1*y3**0
    v2 = q(8)*y1**0*y2**0*y3**2& 
       + q(9)*y1**1*y2**0*y3**1& 
       + q(9)*y1**0*y2**1*y3**1& 
       + q(10)*y1**1*y2**1*y3**0& 
       + q(11)*y1**2*y2**0*y3**0& 
       + q(11)*y1**0*y2**2*y3**0
    v3 = q(12)*y1**0*y2**0*y3**3& 
       + q(13)*y1**1*y2**0*y3**2& 
       + q(13)*y1**0*y2**1*y3**2& 
       + q(14)*y1**1*y2**1*y3**1& 
       + q(15)*y1**2*y2**0*y3**1& 
       + q(15)*y1**0*y2**2*y3**1& 
       + q(16)*y1**2*y2**1*y3**0& 
       + q(16)*y1**1*y2**2*y3**0& 
       + q(17)*y1**3*y2**0*y3**0& 
       + q(17)*y1**0*y2**3*y3**0
       !
     v4 = q(18)*y1**0*y2**0*y3**4& 
       + q(19)*y1**1*y2**0*y3**3& 
       + q(19)*y1**0*y2**1*y3**3& 
       + q(20)*y1**1*y2**1*y3**2& 
       + q(21)*y1**2*y2**0*y3**2& 
       + q(21)*y1**0*y2**2*y3**2& 
       + q(22)*y1**2*y2**1*y3**1& 
       + q(22)*y1**1*y2**2*y3**1& 
       + q(23)*y1**2*y2**2*y3**0& 
       + q(24)*y1**3*y2**0*y3**1& 
       + q(24)*y1**0*y2**3*y3**1& 
       + q(25)*y1**3*y2**1*y3**0& 
       + q(25)*y1**1*y2**3*y3**0& 
       + q(26)*y1**4*y2**0*y3**0& 
       + q(26)*y1**0*y2**4*y3**0
       !
     v5 = q(27)*y1**0*y2**0*y3**5& 
       + q(28)*y1**1*y2**0*y3**4& 
       + q(28)*y1**0*y2**1*y3**4& 
       + q(29)*y1**1*y2**1*y3**3& 
       + q(30)*y1**2*y2**0*y3**3& 
       + q(30)*y1**0*y2**2*y3**3& 
       + q(31)*y1**2*y2**1*y3**2& 
       + q(31)*y1**1*y2**2*y3**2& 
       + q(32)*y1**2*y2**2*y3**1& 
       + q(33)*y1**3*y2**0*y3**2& 
       + q(33)*y1**0*y2**3*y3**2& 
       + q(34)*y1**3*y2**1*y3**1& 
       + q(34)*y1**1*y2**3*y3**1& 
       + q(35)*y1**3*y2**2*y3**0& 
       + q(35)*y1**2*y2**3*y3**0& 
       + q(36)*y1**4*y2**0*y3**1& 
       + q(36)*y1**0*y2**4*y3**1& 
       + q(37)*y1**4*y2**1*y3**0& 
       + q(37)*y1**1*y2**4*y3**0& 
       + q(38)*y1**5*y2**0*y3**0& 
       + q(38)*y1**0*y2**5*y3**0
       !
    v6 = q(39)*y1**0*y2**0*y3**6& 
       + q(40)*y1**1*y2**0*y3**5& 
       + q(40)*y1**0*y2**1*y3**5& 
       + q(41)*y1**1*y2**1*y3**4& 
       + q(42)*y1**2*y2**0*y3**4& 
       + q(42)*y1**0*y2**2*y3**4& 
       + q(43)*y1**2*y2**1*y3**3& 
       + q(43)*y1**1*y2**2*y3**3& 
       + q(44)*y1**2*y2**2*y3**2& 
       + q(45)*y1**3*y2**0*y3**3& 
       + q(45)*y1**0*y2**3*y3**3& 
       + q(46)*y1**3*y2**1*y3**2& 
       + q(46)*y1**1*y2**3*y3**2& 
       + q(47)*y1**3*y2**2*y3**1& 
       + q(47)*y1**2*y2**3*y3**1& 
       + q(48)*y1**3*y2**3*y3**0& 
       + q(49)*y1**4*y2**0*y3**2& 
       + q(49)*y1**0*y2**4*y3**2& 
       + q(50)*y1**4*y2**1*y3**1& 
       + q(50)*y1**1*y2**4*y3**1& 
       + q(51)*y1**4*y2**2*y3**0& 
       + q(51)*y1**2*y2**4*y3**0& 
       + q(52)*y1**5*y2**0*y3**1& 
       + q(52)*y1**0*y2**5*y3**1& 
       + q(53)*y1**5*y2**1*y3**0& 
       + q(53)*y1**1*y2**5*y3**0& 
       + q(54)*y1**6*y2**0*y3**0& 
       + q(54)*y1**0*y2**6*y3**0
       !
    v7 = q(55)*y1**0*y2**0*y3**7& 
       + q(56)*y1**1*y2**0*y3**6& 
       + q(56)*y1**0*y2**1*y3**6& 
       + q(57)*y1**1*y2**1*y3**5& 
       + q(58)*y1**2*y2**0*y3**5& 
       + q(58)*y1**0*y2**2*y3**5& 
       + q(59)*y1**2*y2**1*y3**4& 
       + q(59)*y1**1*y2**2*y3**4& 
       + q(60)*y1**2*y2**2*y3**3& 
       + q(61)*y1**3*y2**0*y3**4& 
       + q(61)*y1**0*y2**3*y3**4& 
       + q(62)*y1**3*y2**1*y3**3& 
       + q(62)*y1**1*y2**3*y3**3& 
       + q(63)*y1**3*y2**2*y3**2& 
       + q(63)*y1**2*y2**3*y3**2& 
       + q(64)*y1**3*y2**3*y3**1& 
       + q(65)*y1**4*y2**0*y3**3& 
       + q(65)*y1**0*y2**4*y3**3& 
       + q(66)*y1**4*y2**1*y3**2& 
       + q(66)*y1**1*y2**4*y3**2& 
       + q(67)*y1**4*y2**2*y3**1& 
       + q(67)*y1**2*y2**4*y3**1& 
       + q(68)*y1**4*y2**3*y3**0& 
       + q(68)*y1**3*y2**4*y3**0& 
       + q(69)*y1**5*y2**0*y3**2& 
       + q(69)*y1**0*y2**5*y3**2& 
       + q(70)*y1**5*y2**1*y3**1& 
       + q(70)*y1**1*y2**5*y3**1& 
       + q(71)*y1**5*y2**2*y3**0& 
       + q(71)*y1**2*y2**5*y3**0& 
       + q(72)*y1**6*y2**0*y3**1& 
       + q(72)*y1**0*y2**6*y3**1& 
       + q(73)*y1**6*y2**1*y3**0& 
       + q(73)*y1**1*y2**6*y3**0& 
       + q(74)*y1**7*y2**0*y3**0& 
       + q(74)*y1**0*y2**7*y3**0
       !
    v8 = q(75)*y1**0*y2**0*y3**8& 
       + q(76)*y1**1*y2**0*y3**7& 
       + q(76)*y1**0*y2**1*y3**7& 
       + q(77)*y1**1*y2**1*y3**6& 
       + q(78)*y1**2*y2**0*y3**6& 
       + q(78)*y1**0*y2**2*y3**6& 
       + q(79)*y1**2*y2**1*y3**5& 
       + q(79)*y1**1*y2**2*y3**5& 
       + q(80)*y1**2*y2**2*y3**4& 
       + q(81)*y1**3*y2**0*y3**5& 
       + q(81)*y1**0*y2**3*y3**5& 
       + q(82)*y1**3*y2**1*y3**4& 
       + q(82)*y1**1*y2**3*y3**4& 
       + q(83)*y1**3*y2**2*y3**3& 
       + q(83)*y1**2*y2**3*y3**3& 
       + q(84)*y1**3*y2**3*y3**2& 
       + q(85)*y1**4*y2**0*y3**4& 
       + q(85)*y1**0*y2**4*y3**4& 
       + q(86)*y1**4*y2**1*y3**3& 
       + q(86)*y1**1*y2**4*y3**3& 
       + q(87)*y1**4*y2**2*y3**2& 
       + q(87)*y1**2*y2**4*y3**2& 
       + q(88)*y1**4*y2**3*y3**1& 
       + q(88)*y1**3*y2**4*y3**1& 
       + q(89)*y1**4*y2**4*y3**0& 
       + q(90)*y1**5*y2**0*y3**3& 
       + q(90)*y1**0*y2**5*y3**3& 
       + q(91)*y1**5*y2**1*y3**2& 
       + q(91)*y1**1*y2**5*y3**2& 
       + q(92)*y1**5*y2**2*y3**1& 
       + q(92)*y1**2*y2**5*y3**1& 
       + q(93)*y1**5*y2**3*y3**0& 
       + q(93)*y1**3*y2**5*y3**0& 
       + q(94)*y1**6*y2**0*y3**2& 
       + q(94)*y1**0*y2**6*y3**2& 
       + q(95)*y1**6*y2**1*y3**1& 
       + q(95)*y1**1*y2**6*y3**1& 
       + q(96)*y1**6*y2**2*y3**0& 
       + q(96)*y1**2*y2**6*y3**0& 
       + q(97)*y1**7*y2**0*y3**1& 
       + q(97)*y1**0*y2**7*y3**1& 
       + q(98)*y1**7*y2**1*y3**0& 
       + q(98)*y1**1*y2**7*y3**0& 
       + q(99)*y1**8*y2**0*y3**0& 
       + q(99)*y1**0*y2**8*y3**0
       !
    mu(1) =(v0+v1+v2+v3+v4+v5+v6+v7+v8)
    !
    mu(3) = 0
    !
    f(1:3) = matmul(mu,tmat)
    !
 end subroutine MLdms2pqr_xy2_linear



 function MLloc2pqr_xy2(r) result(f)

 !return cartesian coordinates of atoms in the user-defined frame for locals specified

    real(ark), intent(in) :: r(molec%ncoords)
    real(ark)             :: f(molec%natoms, 3)

    integer(ik)           :: icart
    real(ark)             :: a0(molec%natoms, 3), cm,mX,mY,R1,R2,alpha,alpha_2,v

    a0 = 0
    !
    select case(trim(molec%coords_transform))
       !
    case('R-RHO-Z')
       !
       a0(2, 1) = -r(1) * cos(r(3)*0.5_ark)
       a0(2, 3) =  r(1) * sin(r(3)*0.5_ark)
       !
       a0(3, 1) = -r(2) * cos(r(3)*0.5_ark)
       a0(3, 3) = -r(2) * sin(r(3)*0.5_ark)
       !
    case('R-RHO-Z-ECKART')
       !
       R1 = r(1)
       R2 = r(2)
       alpha = r(3)
       alpha_2 = alpha*0.5_ark
       v = -atan(cos(alpha_2)*(R1-R2)/(sin(alpha_2)*(R1+R2)))
       !
       mX  = molec%atommasses(1)
       mY  = molec%atommasses(2)
       !
       a0(1,1) =  mY*(cos(v)*cos(alpha_2)*R1+cos(v)*cos(alpha_2)*R2+sin(v)*sin(alpha_2)*R1-sin(v)*sin(alpha_2)*R2)/(mX+2*mY)
       a0(1,3) =  mY*(sin(v)*cos(alpha_2)*R1+sin(v)*cos(alpha_2)*R2-cos(v)*sin(alpha_2)*R1+cos(v)*sin(alpha_2)*R2)/(mX+2*mY)
       a0(2,1) =  -(cos(v)*cos(alpha_2)*R1*mX+cos(v)*cos(alpha_2)*R1*mY-cos(v)*cos(alpha_2)*mY*R2+sin(v)*sin(alpha_2)*R1*mX+sin(v)*sin(alpha_2)*R1*mY+sin(v)*sin(alpha_2)*mY*R2)/(mX+2*mY)
       a0(2,3) =  -(sin(v)*cos(alpha_2)*R1*mX+sin(v)*cos(alpha_2)*R1*mY-sin(v)*cos(alpha_2)*mY*R2-cos(v)*sin(alpha_2)*R1*mX-cos(v)*sin(alpha_2)*R1*mY-cos(v)*sin(alpha_2)*mY*R2)/(mX+2*mY)
       a0(3,1) =  -(cos(v)*cos(alpha_2)*R2*mX+cos(v)*cos(alpha_2)*mY*R2-cos(v)*cos(alpha_2)*R1*mY-sin(v)*sin(alpha_2)*R2*mX-sin(v)*sin(alpha_2)*mY*R2-sin(v)*sin(alpha_2)*R1*mY)/(mX+2*mY)
       a0(3,3) =  -(sin(v)*cos(alpha_2)*R2*mX+sin(v)*cos(alpha_2)*mY*R2-sin(v)*cos(alpha_2)*R1*mY+cos(v)*sin(alpha_2)*R2*mX+cos(v)*sin(alpha_2)*mY*R2+cos(v)*sin(alpha_2)*R1*mY)/(mX+2*mY)
       !
    case default 
       stop 'MLloc2pqr_xy2: illegla coordinate type'
    end select
    
    do icart = 1, 3
       cm = sum(a0(1:molec%natoms, icart) * molec%atommasses(1:molec%natoms)) / sum(molec%atommasses(1:molec%natoms))
       a0(1:molec%natoms, icart) = a0(1:molec%natoms, icart) - cm
    end do

    f(1:molec%natoms, 1:3) = a0(1:molec%natoms, 1:3)

 end function MLloc2pqr_xy2


 function MLloc2pqr_xyz(r) result(f)

 !return cartesian coordinates of atoms in the user-defined frame for locals specified

    real(ark), intent(in) :: r(molec%ncoords)
    real(ark)             :: f(molec%natoms, 3)

    integer(ik)           :: icart
    real(ark)             :: a0(molec%natoms, 3), cm
    !
    a0 = 0
    !
    select case(trim(molec%coords_transform))
       !
    case('R1-Z-R2-ALPHA','R1-Z-R2-RHO')
       !
       a0(2, 1) =  0
       a0(2, 3) =  r(1)
       !
       a0(3, 1) =  r(2) * sin(r(3))
       a0(3, 3) =  r(2) * cos(r(3))
       !
    case default 
       stop 'MLloc2pqr_xyz: illegla coordinate type'
    end select
    
    do icart = 1, 3
       cm = sum(a0(1:molec%natoms, icart) * molec%atommasses(1:molec%natoms)) / sum(molec%atommasses(1:molec%natoms))
       a0(1:molec%natoms, icart) = a0(1:molec%natoms, icart) - cm
    end do

    f(1:molec%natoms, 1:3) = a0(1:molec%natoms, 1:3)

 end function MLloc2pqr_xyz



  !
  ! Defining potential energy function 
  !    LARGE SCALE AB INITIO CALCULATIONS FOR C3
  !    M. MLADENOVIC ET AL.
  !    J.CHEM. PHYS. VOL.101 N0.7 1994
  !
  function MLpoten_c3_mladenovic(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik),parameter  :: N = 52
   !
   real(ark)    ::   P1, P2, th, prod,re, sum, q1, q2,theta
   integer(ik)  ::   k(1:52), l(1:52), m(1:52), j
     !
     k = (/2,3,4,5,6,0,0,0,1,2,1,2,3,4,0,0,0,0,0,0,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2/)
     l = (/0,0,0,0,0,2,4,6,2,2,4,4,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,4,4,4,4,2,2,2,2,2,2,2,2/)
     m = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,4,6,8,1,1,2,2,2,2,4,4,4,4,6,6,6,6,8,8,8,8,2,4,6,8,2,4,6,8,2,4,6,8,2,4,6,8/)
     !
     re = force( 1)/bohr
     !
     q1=local(1)/bohr
     q2=local(2)/bohr
     theta = local(3)
     !
     if  (molec%Ndihedrals>1) then
       !
       theta = asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
       !
     endif 
     !
     !re=2.453016666d0
     !
     P1=(q1+q2-2.0_ark*re)*sqrt(0.5_ark)
     P2=(q1-q2)*sqrt(0.5_ark)
     !
     sum=0
     do j=1, n
        prod=force(j+1)*(P1**k(j))*(P2**l(j))*(theta**m(j))
        sum=sum+prod
     end do
     !
     f=sum*1.0e-11/planck/vellgt
 
   end function MLpoten_c3_mladenovic


  !
  ! B. Schröder and P. Sebald, JCP 2016 http://dx.doi.org/10.1063/1.4940780]
  !
  function MLpoten_c3_R_theta(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   real(ark)    ::   re,alphae,vpot, q(3),theta,y(3),p1,p2,v1(3),v2(3),alpha,x1,x2,cosalpha1
   integer(ik)  ::   i,k_ind(2)
   real(ark),parameter :: tocm = 219474.63067_ark
   character(len=cl)  :: txt = 'MLpoten_c3_R_theta'
     !
     re = force( 1)
     alphae = force(2)
     !
     q(1)=local(1)-re
     q(2)=local(2)-re
     q(3) = (pi-local(3))-alphae*pi/180.0_ark
     !
     !p1=local(1)-re
     !p2=local(2)-re
     !
     !q(1) = (p1+p2)*sqrt(0.5_ark)/bohr
     !q(2) = (p1-p2)*sqrt(0.5_ark)/bohr
     !
     q(1:2) = q(1:2)/bohr
     !
     if  (molec%Ndihedrals>1) then
       !
       theta = asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
       q(3) = theta
       !
       v1(:) = xyz(2,:)-xyz(1,:)
       v2(:) = xyz(3,:)-xyz(1,:)
       !
       x1 = sqrt(sum(v1(:)**2))
       x2 = sqrt(sum(v2(:)**2))
       !
       cosalpha1 = sum( v1(:)*v2(:))/(x1*x2)
       !
       alpha = aacos(cosalpha1,txt)
       theta  = pi- alpha
       q(3) = theta
       !
     endif
     !
     vpot = 0 
     !
     do i=3,molec%parmax
       !
       y(1:3) = q(1:3)**molec%pot_ind(1:3,i)
       !
       vpot=vpot+force(i)*product(y(1:3))
       !
       if (molec%pot_ind(1,i)/=molec%pot_ind(2,i)) then
         !
         k_ind(1) = molec%pot_ind(2,i)
         k_ind(2) = molec%pot_ind(1,i)
         !
         y(1:2) = q(1:2)**k_ind(1:2)
         !
         vpot=vpot+force(i)*product(y(1:3))
         !
       endif
       !
     enddo
     !
     !f=vpot*1.0e-11/planck/vellgt
     !
     f=vpot*tocm
     
 
   end function MLpoten_c3_R_theta




   function MLpoten_h2s_dvr3d(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   !
!     caculation d2 at level ccsd(t)/aug-cc-pvqz
!     and d3 , d4 at level mp2/aug-cc-pvqz .
!     geometry from hoy , mills , strey.
!     units: aJ and Angstrom

!     rz = oh equilibrium value
!     rho = equilibrium value of pi - bond  angle(theta)

!     mh = 1.007825 aem      mo = 15.99491 aem      ms = 31.97210 aem

      real(ark),parameter :: pre=50340.2_ark,ridber=109737.309_ark,hartree=2*ridber/pre
      real(ark),parameter :: tocm = 219474.63067_ark
      !
      real(ark) :: x1,a,rz,rho1,d11 ,d12 ,d13 ,d22 ,d23 ,d33 ,d111,d112,&
                  d113,d122,d123,d133,d222,d223,d233,d333,d1111,d1112,d1113,d1122,d1123,&
                  d1133,d1222,d1223,d1233,d1333,d2222,d2223,d2233,d2333,d3333,dr,ds,y1,y2,y3,theta
      real(ark) :: rho,v2,v3,v4,v

      !
      !toang = 0.5291772
      !
      ! cmtoau = 219474.624
      x1 = 1.0_ark
      !
      a    = force( 1)
      rz   = force( 2)
      rho1 = force( 3)   !  180.0-92.376
      !
      d11  = force( 4)
      d12  = force( 5)
      d13  = force( 6)
      d22  = force( 7)
      d23  = force( 8)
      d33  = force( 9)
      d111 = force(10)
      d112 = force(11)
      d113 = force(12)
      d122 = force(13)
      d123 = force(14)
      d133 = force(15)
      d222 = force(16)
      d223 = force(17)
      d233 = force(18)
      d333 = force(19)
      d1111= force(20)
      d1112= force(21)
      d1113= force(22)
      d1122= force(23)
      d1123= force(24)
      d1133= force(25)
      d1222= force(26)
      d1223= force(27)
      d1233= force(28)
      d1333= force(29)
      d2222= force(30)
      d2223= force(31)
      d2233= force(32)
      d2333= force(33)
      d3333= force(34)
                    
!     find value for dr and ds
      !dr = toang*q1 - rz
      !ds = toang*q2 - rz
      dr = local(1) - rz
      ds = local(2) - rz
      theta = local(3)

!     transform to morse coordinates
      y1 = x1 - exp(-a * dr)
      y2 = x1 - exp(-a * ds)
      rho=rho1*pi/180.0_ark
      y3 = -(cos(theta) + cos(rho))

!     now for the potential
      v2 =  d11*y1*y1         &
       + 2*d12*y1*y2          &
       + 2*d13*y1*y3          &
       +   d22*y2*y2          &
       + 2*d23*y2*y3          &
       +   d33*y3*y3           
                                
     v3 =  d111*y1*y1*y1      &
       + 3*d112*y1*y1*y2      &
       + 3*d113*y1*y1*y3      &
       + 3*d122*y1*y2*y2      &
       + 6*d123*y1*y2*y3      &
       + 3*d133*y1*y3*y3      &
       +   d222*y2*y2*y2      &
       + 3*d223*y2*y2*y3      &
       + 3*d233*y2*y3*y3      &
       +   d333*y3*y3*y3       
                                
     v4 =  d1111*y1*y1*y1*y1  &
       + 4*d1112*y1*y1*y1*y2  &
       + 4*d1113*y1*y1*y1*y3  &
       + 6*d1122*y1*y1*y2*y2  &
       +12*d1123*y1*y1*y2*y3  &
       + 6*d1133*y1*y1*y3*y3  &
       + 4*d1222*y1*y2*y2*y2  &
       +12*d1223*y1*y2*y2*y3  &
       +12*d1233*y1*y2*y3*y3  &
       + 4*d1333*y1*y3*y3*y3  &
       +   d2222*y2*y2*y2*y2  &
       + 4*d2223*y2*y2*y2*y3  &
       + 6*d2233*y2*y2*y3*y3  &
       + 4*d2333*y2*y3*y3*y3  &
       +   d3333*y3*y3*y3*y3
       !
       v = v2/2.0_ark +v3/6.0_ark + v4/24.0_ark
       !
       f = v*1.0e-11/planck/vellgt ! /hartree*tocm
       !
 
 end function MLpoten_h2s_dvr3d



 recursive subroutine MLdipole_xy2_lorenzo(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: k,imu,iterm
    real(ark)             :: y(3), mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,re,ae
    real(ark)             :: xyz0(natoms,3),xi(3),mu_t
    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      select case(trim(molec%coords_transform))
      case default
         write (out,"('MLdipole_xy2_lorenzo: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLdipole_xy2_lorenzo - bad coord. type'
      case('R-RHO-Z','R-RHO-Z-ECKART')
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
    re = extF%coef(1,1)
    ae = extF%coef(2,1)*pi/180.0_ark
    !
    y(1) = ( r1 + r2 )*0.5_ark - re
    y(2) = ( r1 - r2 )
    y(3) = ae - alpha
    !
    mu = 0
    !
    do imu = 1, 2
       !
       do iterm =  3, extF%nterms(imu)
          !
          xi(1:3) = y(1:3) ** extF%term(1:3, iterm, imu)
          !
          mu_t = product(xi(1:3))
          !
          mu(imu) = mu(imu) + extF%coef(iterm, imu)*mu_t
          !
       end do
       !
    end do
    !
    mu(3) = 0
    !
    f(1:3) = matmul(mu,tmat)
    !
 end subroutine MLdipole_xy2_lorenzo



  !
  ! Defining potential energy function
  !
  function MLpoten_xyz_tyuterev(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)            :: r1,r2,alpha,xcos,v,v1,v2,v3,v4,v5,v6,v7,v8

   real(ark)            :: aa1,aa2,re1,re2,alphae,xst,y1,y2,y3,xs1,xs2,v0,vp1,vp2,vp3
   real(ark)            :: g1,g2,b1,b2,rhh,vhh
   integer(ik)          :: N
   real(ark)            :: ycos,v_t,q1,q2
   real(ark)            :: a1(3),a2(3),t1,t2,w(3),cosalpha
   character(len=cl)    :: txt
   !
   if (verbose>=6) write(out,"('MLpoten_xyz_tyuterev/start')")
   !
   r1  = local(1) ;  r2    = local(2) ;  alpha = local(3)
   !
   if (molec%Nangles>0) then
     alphae  = molec%alphaeq(1)
   else
     alphae = pi+(-molec%taueq(1)+molec%taueq(2))
     alpha = pi - local(3) 
   endif 
   !
   re1    = force(1)
   re2    = force(2)
   alphae = force(3)
   aa1  = force(4)
   aa2  = force(5)
   !
   b1   = force(6)
   b2   = force(7)
   g1   = force(8)
   g2   = force(9)
   !
   rhh=sqrt(r1**2+r2**2-2._ark*r1*r2*cos(alpha))
   vhh=b1*exp(-g1*rhh)+b2*exp(-g2*rhh**2)
   !
   ! calculate potential energy function values
   !
   y1=1.0_ark-exp(-aa1*(r1-re1))
   y2=1.0_ark-exp(-aa1*(r2-re2))
   !
   y3=cos(alpha)-cos(alphae)
   !
   N = size(force)
   !
   v4 = 0 ; v5 = 0 ; v6 = 0 ; v7 = 0 ; v8 = 0
   !
   v0 = force( 10)*y1**0*y2**0*y3**0
   v1 = force( 11)*y1**0*y2**0*y3**1& 
      + force( 12)*y1**1*y2**0*y3**0& 
      + force( 13)*y1**0*y2**1*y3**0
   v2 = force( 14)*y1**0*y2**0*y3**2& 
      + force( 15)*y1**1*y2**0*y3**1& 
      + force( 16)*y1**0*y2**1*y3**1& 
      + force( 17)*y1**1*y2**1*y3**0& 
      + force( 18)*y1**2*y2**0*y3**0& 
      + force( 19)*y1**0*y2**2*y3**0
   v3 = force( 20)*y1**0*y2**0*y3**3& 
      + force( 21)*y1**1*y2**0*y3**2& 
      + force( 22)*y1**0*y2**1*y3**2& 
      + force( 23)*y1**1*y2**1*y3**1& 
      + force( 24)*y1**2*y2**0*y3**1& 
      + force( 25)*y1**0*y2**2*y3**1& 
      + force( 26)*y1**2*y2**1*y3**0& 
      + force( 27)*y1**1*y2**2*y3**0& 
      + force( 28)*y1**3*y2**0*y3**0& 
      + force( 29)*y1**0*y2**3*y3**0
   v4 = force( 30)*y1**0*y2**0*y3**4& 
      + force( 31)*y1**1*y2**0*y3**3& 
      + force( 32)*y1**0*y2**1*y3**3& 
      + force( 33)*y1**1*y2**1*y3**2& 
      + force( 34)*y1**2*y2**0*y3**2& 
      + force( 35)*y1**0*y2**2*y3**2& 
      + force( 36)*y1**2*y2**1*y3**1& 
      + force( 37)*y1**1*y2**2*y3**1& 
      + force( 38)*y1**2*y2**2*y3**0& 
      + force( 39)*y1**3*y2**0*y3**1& 
      + force( 40)*y1**0*y2**3*y3**1& 
      + force( 41)*y1**3*y2**1*y3**0& 
      + force( 42)*y1**1*y2**3*y3**0& 
      + force( 43)*y1**4*y2**0*y3**0& 
      + force( 44)*y1**0*y2**4*y3**0
   v5 = force( 45)*y1**0*y2**0*y3**5& 
      + force( 46)*y1**1*y2**0*y3**4& 
      + force( 47)*y1**0*y2**1*y3**4& 
      + force( 48)*y1**1*y2**1*y3**3& 
      + force( 49)*y1**2*y2**0*y3**3& 
      + force( 50)*y1**0*y2**2*y3**3& 
      + force( 51)*y1**2*y2**1*y3**2& 
      + force( 52)*y1**1*y2**2*y3**2& 
      + force( 53)*y1**2*y2**2*y3**1& 
      + force( 54)*y1**3*y2**0*y3**2& 
      + force( 55)*y1**0*y2**3*y3**2& 
      + force( 56)*y1**3*y2**1*y3**1& 
      + force( 57)*y1**1*y2**3*y3**1& 
      + force( 58)*y1**3*y2**2*y3**0& 
      + force( 59)*y1**2*y2**3*y3**0& 
      + force( 60)*y1**4*y2**0*y3**1& 
      + force( 61)*y1**0*y2**4*y3**1& 
      + force( 62)*y1**4*y2**1*y3**0& 
      + force( 63)*y1**1*y2**4*y3**0& 
      + force( 64)*y1**5*y2**0*y3**0& 
      + force( 65)*y1**0*y2**5*y3**0
   v6 = force( 66)*y1**0*y2**0*y3**6& 
      + force( 67)*y1**1*y2**0*y3**5& 
      + force( 68)*y1**0*y2**1*y3**5& 
      + force( 69)*y1**1*y2**1*y3**4& 
      + force( 70)*y1**2*y2**0*y3**4& 
      + force( 71)*y1**0*y2**2*y3**4& 
      + force( 72)*y1**2*y2**1*y3**3& 
      + force( 73)*y1**1*y2**2*y3**3& 
      + force( 74)*y1**2*y2**2*y3**2& 
      + force( 75)*y1**3*y2**0*y3**3& 
      + force( 76)*y1**0*y2**3*y3**3& 
      + force( 77)*y1**3*y2**1*y3**2& 
      + force( 78)*y1**1*y2**3*y3**2& 
      + force( 79)*y1**3*y2**2*y3**1& 
      + force( 80)*y1**2*y2**3*y3**1& 
      + force( 81)*y1**3*y2**3*y3**0& 
      + force( 82)*y1**4*y2**0*y3**2& 
      + force( 83)*y1**0*y2**4*y3**2& 
      + force( 84)*y1**4*y2**1*y3**1& 
      + force( 85)*y1**1*y2**4*y3**1& 
      + force( 86)*y1**4*y2**2*y3**0& 
      + force( 87)*y1**2*y2**4*y3**0& 
      + force( 88)*y1**5*y2**0*y3**1& 
      + force( 89)*y1**0*y2**5*y3**1& 
      + force( 90)*y1**5*y2**1*y3**0& 
      + force( 91)*y1**1*y2**5*y3**0& 
      + force( 92)*y1**6*y2**0*y3**0& 
      + force( 93)*y1**0*y2**6*y3**0
   v7 = force( 94)*y1**0*y2**0*y3**7&
      + force( 95)*y1**1*y2**0*y3**6&
      + force( 96)*y1**0*y2**1*y3**6&
      + force( 97)*y1**1*y2**1*y3**5&
      + force( 98)*y1**2*y2**0*y3**5&
      + force( 99)*y1**0*y2**2*y3**5&
      + force(100)*y1**2*y2**1*y3**4&
      + force(101)*y1**1*y2**2*y3**4&
      + force(102)*y1**2*y2**2*y3**3&
      + force(103)*y1**3*y2**0*y3**4&
      + force(104)*y1**0*y2**3*y3**4&
      + force(105)*y1**3*y2**1*y3**3&
      + force(106)*y1**1*y2**3*y3**3&
      + force(107)*y1**3*y2**2*y3**2&
      + force(108)*y1**2*y2**3*y3**2&
      + force(109)*y1**3*y2**3*y3**1&
      + force(110)*y1**4*y2**0*y3**3&
      + force(111)*y1**0*y2**4*y3**3&
      + force(112)*y1**4*y2**1*y3**2&
      + force(113)*y1**1*y2**4*y3**2&
      + force(114)*y1**4*y2**2*y3**1&
      + force(115)*y1**2*y2**4*y3**1&
      + force(116)*y1**4*y2**3*y3**0&
      + force(117)*y1**3*y2**4*y3**0&
      + force(118)*y1**5*y2**0*y3**2&
      + force(119)*y1**0*y2**5*y3**2&
      + force(120)*y1**5*y2**1*y3**1&
      + force(121)*y1**1*y2**5*y3**1&
      + force(122)*y1**5*y2**2*y3**0&
      + force(123)*y1**2*y2**5*y3**0&
      + force(124)*y1**6*y2**0*y3**1&
      + force(125)*y1**0*y2**6*y3**1&
      + force(126)*y1**6*y2**1*y3**0&
      + force(127)*y1**1*y2**6*y3**0&
      + force(128)*y1**7*y2**0*y3**0&
      + force(129)*y1**0*y2**7*y3**0
   v8 = force(130)*y1**0*y2**0*y3**8&
      + force(131)*y1**1*y2**0*y3**7&
      + force(132)*y1**0*y2**1*y3**7&
      + force(133)*y1**1*y2**1*y3**6&
      + force(134)*y1**2*y2**0*y3**6&
      + force(135)*y1**0*y2**2*y3**6&
      + force(136)*y1**2*y2**1*y3**5&
      + force(137)*y1**1*y2**2*y3**5&
      + force(138)*y1**2*y2**2*y3**4&
      + force(139)*y1**3*y2**0*y3**5&
      + force(140)*y1**0*y2**3*y3**5&
      + force(141)*y1**3*y2**1*y3**4&
      + force(142)*y1**1*y2**3*y3**4&
      + force(143)*y1**3*y2**2*y3**3&
      + force(144)*y1**2*y2**3*y3**3&
      + force(145)*y1**3*y2**3*y3**2&
      + force(146)*y1**4*y2**0*y3**4&
      + force(147)*y1**0*y2**4*y3**4&
      + force(148)*y1**4*y2**1*y3**3&
      + force(149)*y1**1*y2**4*y3**3&
      + force(150)*y1**4*y2**2*y3**2&
      + force(151)*y1**2*y2**4*y3**2&
      + force(152)*y1**4*y2**3*y3**1&
      + force(153)*y1**3*y2**4*y3**1&
      + force(154)*y1**4*y2**4*y3**0&
      + force(155)*y1**5*y2**0*y3**3&
      + force(156)*y1**0*y2**5*y3**3&
      + force(157)*y1**5*y2**1*y3**2&
      + force(158)*y1**1*y2**5*y3**2&
      + force(159)*y1**5*y2**2*y3**1&
      + force(160)*y1**2*y2**5*y3**1&
      + force(161)*y1**5*y2**3*y3**0&
      + force(162)*y1**3*y2**5*y3**0&
      + force(163)*y1**6*y2**0*y3**2&
      + force(164)*y1**0*y2**6*y3**2&
      + force(165)*y1**6*y2**1*y3**1&
      + force(166)*y1**1*y2**6*y3**1&
      + force(167)*y1**6*y2**2*y3**0&
      + force(168)*y1**2*y2**6*y3**0&
      + force(169)*y1**7*y2**0*y3**1&
      + force(170)*y1**0*y2**7*y3**1&
      + force(171)*y1**7*y2**1*y3**0&
      + force(172)*y1**1*y2**7*y3**0&
      + force(173)*y1**8*y2**0*y3**0&
      + force(174)*y1**0*y2**8*y3**0
    !
    f=v0+v1+v2+v3+v4+v5+v6+v7+v8+vhh
    !
    if (verbose>=6) write(out,"('MLpoten_xyz_tyuterev/end')")

 end function MLpoten_xyz_tyuterev


 !returns electric dipole moment cartesian coordinates in the user-defined frame for locals specified
 !
 recursive subroutine MLdms2pqr_xyz_coeff(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: k
    real(ark)             :: y1,y2,y3, mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,re1,re2,ae,a0,b0
    real(ark)             :: p(3:extF%nterms(1)),q(5:extF%nterms(2)),v0,v1,v2,v3,v4,v5,v6,v7,v8,xyz0(natoms,3)
    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      select case(trim(molec%coords_transform))
      case default
         write (out,"('MLdms2pqr_xyz_coeff: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLdms2pqr_xyz_coeff - bad coord. type'
      case('R1-Z-R2-ALPHA','R1-Z-R2-RHO')
         !
         xyz0 = MLloc2pqr_xyz(local)
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
    mu(2) = v0+v1+v2+v3+v4+v5+v6+v7+v8
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
    mu(3) = 0
    !
    f(1:3) = matmul(mu,tmat)
    !
 end subroutine MLdms2pqr_xyz_coeff


end module pot_xy2

!
!  This unit defines all specific routines for a fouratomic molecule of ZXY2 type
!
module pot_zxy2
  use accuracy
  use moltype

  implicit none

  public MLpoten_sohf,MLpoten_zxy2_andrey_01,MLpoten_zxy2_mep_r_alpha_rho_powers
  public MLdms2xyz_zxy2_symadap_powers,ML_MEP_zxy2_R_rho,MLpoten_zxy2_andrey_coeff,ML_MEP_zxy2_rho_coeff,MLpoten_h2cs_tz_damp1
  public MLpoten_h2cs_damp,MLpoten_zxy2_mlt,MLpoten_h2cs_damp_scaling,MLpoten_zxy2_morse_cos,&
         MLpoten_zxy2_mep_r_alpha_rho_powers_iso
  private
 
  integer(ik), parameter :: verbose     = 3                          ! Verbosity level


  contains

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
     dst(1) = molec%mep_params(1)+molec%mep_params(2)*cs+molec%mep_params(3)*cs**2+molec%mep_params(4)*cs**3 &
              +molec%mep_params( 5)*cs**4
     dst(2) = molec%mep_params(6)+molec%mep_params(7)*cs+molec%mep_params(8)*cs**2+molec%mep_params(9)*cs**3 &
              +molec%mep_params(10)*cs**4
     dst(3) = dst(2)
     dst(4) = molec%mep_params(11)*rad+molec%mep_params(12)*cs+molec%mep_params(13)*cs**2+molec%mep_params(14)*cs**3+&
              molec%mep_params(15)*cs**4
     dst(5) = dst(4)
     dst(6) = x
     !
  end function ML_MEP_zxy2_rho_coeff





function MLpoten_h2cs_damp(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
 !
 integer(ik)          :: iterm,k_ind(6)
 real(ark)            :: y(6),xi(6),vshort,vlong,vdamp
 real(ark)            :: re1,re2,ae,De1,De2,De3,beta1,beta2,beta3,damp1,damp2
 !
 !
 re1   = force(1)
 re2   = force(2)
 ae    = force(3)*pi/180.0_ark
 De1   = force(4)
 beta1 = force(5)
 De2   = force(6)
 beta2 = force(7)
 damp1 = force(8)
 damp2 = force(9)
 !
 ! transformation to the standard form 
 !
 xi = from_local_to_r1r2r3a1a2tau(local,6)
 !
 y(1) = xi(1)-re1
 y(2) = xi(2)-re2
 y(3) = xi(3)-re2
 y(4) = xi(4)-ae
 y(5) = xi(5)-ae
 y(6) = 1.0_ark+cos(xi(6))
 !
 vlong = De1*(1.0_ark-exp(-beta1*y(1)))**2 + &
        De2*(1.0_ark-exp(-beta2*y(2)))**2 + &
        De2*(1.0_ark-exp(-beta2*y(3)))**2
 !
 !
 vdamp = exp( -damp1*y(1)**4 -damp2*y(2)**4 -damp2*y(3)**4 )
 !
 !
 vshort = 0.0_ark
 !
 do iterm = 10, molec%parmax
  !
  xi(1:6) = y(1:6)**molec%pot_ind(1:6,iterm)
  !
  vshort = vshort + force(iterm)*product(xi(1:6))
  !
  if (molec%pot_ind(2,iterm)/=molec%pot_ind(3,iterm).or.molec%pot_ind(4,iterm)/=molec%pot_ind(5,iterm)) then
    !
    k_ind(1) = molec%pot_ind(1,iterm)
    k_ind(2) = molec%pot_ind(3,iterm)
    k_ind(3) = molec%pot_ind(2,iterm)
    k_ind(4) = molec%pot_ind(5,iterm)
    k_ind(5) = molec%pot_ind(4,iterm)
    k_ind(6) = molec%pot_ind(6,iterm)
    !
    xi(1:6) = y(1:6)**k_ind(1:6)
    !
    vshort = vshort + force(iterm)*product(xi(1:6))
    !
  endif
  !
 enddo
 !
 f = vlong + vshort*vdamp
 !
end function MLpoten_h2cs_damp





function MLpoten_h2cs_damp_scaling(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
 !
 integer(ik)          :: iterm,k_ind(6)
 real(ark)            :: y(6),xi(6),vshort,vlong,vdamp
 real(ark)            :: re1,re2,ae,De1,De2,De3,beta1,beta2,beta3,damp1,damp2,s(6),sc1,sc2,sc3
 !
 !
 re1   = force(1)
 re2   = force(2)
 ae    = force(3)*pi/180.0_ark
 De1   = force(4)
 beta1 = force(5)
 De2   = force(6)
 beta2 = force(7)
 damp1 = force(8)
 damp2 = force(9)
 !
 ! transformation to the standard form 
 !
 xi = from_local_to_r1r2r3a1a2tau(local,6)
 !
 s(1) = re1/molec%local_eq(1)
 s(2) = re2/molec%local_eq(2)
 s(3) = re2/molec%local_eq(3)
 s(4) = ae/molec%local_eq(4)
 s(5) = ae/molec%local_eq(5)
 !
 sc1 = force(10)
 sc2 = force(11)
 sc3 = force(12)
 !
 y(1) = xi(1)*sc1-re1
 y(2) = xi(2)*sc2-re2
 y(3) = xi(3)*sc2-re2
 y(4) = xi(4)*sc3-ae
 y(5) = xi(5)*sc3-ae
 y(6) = 1.0_ark+cos(xi(6))
 !
 vlong = De1*(1.0_ark-exp(-beta1*y(1)))**2 + &
        De2*(1.0_ark-exp(-beta2*y(2)))**2 + &
        De2*(1.0_ark-exp(-beta2*y(3)))**2
 !
 !
 vdamp = exp( -damp1*y(1)**4 -damp2*y(2)**4 -damp2*y(3)**4 )
 !
 !
 vshort = 0.0_ark
 !
 do iterm = 13, molec%parmax
  !
  xi(1:6) = y(1:6)**molec%pot_ind(1:6,iterm)
  !
  vshort = vshort + force(iterm)*product(xi(1:6))
  !
  if (molec%pot_ind(2,iterm)/=molec%pot_ind(3,iterm).or.molec%pot_ind(4,iterm)/=molec%pot_ind(5,iterm)) then
    !
    k_ind(1) = molec%pot_ind(1,iterm)
    k_ind(2) = molec%pot_ind(3,iterm)
    k_ind(3) = molec%pot_ind(2,iterm)
    k_ind(4) = molec%pot_ind(5,iterm)
    k_ind(5) = molec%pot_ind(4,iterm)
    k_ind(6) = molec%pot_ind(6,iterm)
    !
    xi(1:6) = y(1:6)**k_ind(1:6)
    !
    vshort = vshort + force(iterm)*product(xi(1:6))
    !
  endif
  !
 enddo
 !
 f = vlong + vshort*vdamp
 !
end function MLpoten_h2cs_damp_scaling



function MLpoten_zxy2_morse_cos(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
 !
 integer(ik)          :: iterm,k_ind(6)
 real(ark)            :: y(6),xi(6),v
 real(ark)            :: re1,re2,ae,De1,De2,De3,a1,a2
 !
 !
 re1   = force(1)
 re2   = force(2)
 ae    = force(3)*pi/180.0_ark
 a1 = force(4)
 a2 = force(5)
 !
 ! transformation to the standard form 
 !
 y(1) = 1.0_ark-exp(-a1*(local(1)-re1)) 
 y(2) = 1.0_ark-exp(-a2*(local(2)-re2))
 y(3) = 1.0_ark-exp(-a2*(local(3)-re2))
 y(4) = local(4)-ae
 y(5) = local(5)-ae
 y(6) = 1.0_ark+cos(local(6))
 !
 v = 0.0_ark
 !
 do iterm = 6, molec%parmax
  !
  if (abs(force(iterm))<small_) cycle
  !
  xi(1:6) = y(1:6)**molec%pot_ind(1:6,iterm)
  !
  v = v + force(iterm)*product(xi(1:6))
  !
  if (molec%pot_ind(2,iterm)/=molec%pot_ind(3,iterm).or.molec%pot_ind(4,iterm)/=molec%pot_ind(5,iterm)) then
    !
    k_ind(1) = molec%pot_ind(1,iterm)
    k_ind(2) = molec%pot_ind(3,iterm)
    k_ind(3) = molec%pot_ind(2,iterm)
    k_ind(4) = molec%pot_ind(5,iterm)
    k_ind(5) = molec%pot_ind(4,iterm)
    k_ind(6) = molec%pot_ind(6,iterm)
    !
    xi(1:6) = y(1:6)**k_ind(1:6)
    !
    v = v + force(iterm)*product(xi(1:6))
    !
  endif
  !
 enddo
 !
 f = v
 !
end function MLpoten_zxy2_morse_cos



  function from_local_to_r1r2r3a1a2tau(src,ndst) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst

    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: phi,cosphi,theta1,theta2,theta3
    !
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('from_local_to_r1r2r3a1a2tau/start')") 
    !
    nsrc = size(src)
    !
    if (nsrc==ndst) dst = src
    !
    select case(trim(molec%coords_transform))
       !
    case('R-THETA-DELTA')
      !
      dst(1:3) = src(1:3)
      !
      theta1 = src(4)
      theta2 = src(5)
      theta3 = src(6)
      !
      cosphi = (-cos(theta2)*cos(theta1)+cos(theta3))/(sin(theta2)*sin(theta1))
      !
      if ( abs(cosphi)>1.0_ark+sqrt(small_) ) then 
         !
         write (out,"('from_local_to_r1r2r3a1a2tau: cosphi>1: ',f18.8)") cosphi
         stop 'from_local_to_r1r2r3a1a2tau - bad cosphi'
         !
      elseif ( cosphi>=1.0_ark) then 
         !
         phi = 0.0_ark
         !
      elseif ( cosphi<=-1.0_ark) then 
         !
         phi = pi
         !
      else 
         !
         phi = acos(cosphi)
         !
      endif 
      !
      if (src(7)>0) phi = 2.0_ark*pi-phi
      !
      dst(4) = theta1
      dst(5) = theta2
      dst(6) = phi
      !
    end select
    !
    !
    if (verbose>=5) write(out,"('from_local_to_r1r2r3a1a2tau/end')") 
    !
    !
  end function from_local_to_r1r2r3a1a2tau
  !
  !
  !
  function MLpoten_h2cs_tz_damp1(ncoords,natoms,x,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  x(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
      !
      integer(ik)          :: iterm,k_ind(6)
      real(ark)            :: y(6),xi(6),vshort,vlong,vdamp,local(6)
      real(ark)            :: De1,De2,De3,beta1,beta2,beta3,damp1,damp2,damp3,damp4
      !
      local = from_local_to_r1r2r3a1a2tau(x,6)
      !
      y(1:5) = local(1:5)-molec%local_eq(1:5)
      y(6)   = 1.0_ark+cos(local(6))
      !
      vshort = 0.0_ark
      !
      do iterm=1,size(force)
       !
       xi(1:6) = y(1:6)**molec%pot_ind(1:6,iterm)
       !
       vshort = vshort + force(iterm)*product(xi(1:6))
       !
       if (molec%pot_ind(2,iterm)/=molec%pot_ind(3,iterm).or.molec%pot_ind(4,iterm)/=molec%pot_ind(5,iterm)) then
         !
         k_ind(1) = molec%pot_ind(1,iterm)
         k_ind(2) = molec%pot_ind(3,iterm)
         k_ind(3) = molec%pot_ind(2,iterm)
         k_ind(4) = molec%pot_ind(5,iterm)
         k_ind(5) = molec%pot_ind(4,iterm)
         k_ind(6) = molec%pot_ind(6,iterm)
         !
         xi(1:6) = y(1:6)**k_ind(1:6)
         !
         vshort = vshort + force(iterm)*product(xi(1:6))
         !
       endif
       !
      enddo
      !
      De1   = 665036.229040111299_ark
      De2   = 917124.424633309245_ark
      De3   =     96.160229704135_ark
      beta1 =      1.0_ark
      beta2 =      1.0_ark
      beta3 =      2.0_ark
      !
      vlong = De1*(1.0_ark-exp(-beta1*y(1)))**2 + &
             De2*(1.0_ark-exp(-beta2*y(2)))**2 + &
             De2*(1.0_ark-exp(-beta2*y(3)))**2 + &
             De3*(1.0_ark-exp(-beta3*y(4)))**2 + &
             De3*(1.0_ark-exp(-beta3*y(5)))**2
      !
      damp1 = 0.5_ark
      damp2 = 0.5_ark
      damp3 = 0.5_ark
      damp4 = 0.5_ark
      !
      vdamp = exp( -damp1*y(1)**2-damp2*y(1)**4 &
                  -damp3*y(2)**2-damp4*y(2)**4 &
                  -damp3*y(3)**2-damp4*y(3)**4 )
      !
      f = vlong + vshort*vdamp
      !
  end function MLpoten_h2cs_tz_damp1



  !
  ! Defining potential energy function (built for H2CS)

  function MLpoten_zxy2_andrey_coeff(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   f = h2co_pes_6m_morse_delta_cos_mep(local,force)
   !
 end function MLpoten_zxy2_andrey_coeff
 !
 !
 !
 !
 function h2co_pes_6m_morse_delta_cos_mep(x,param) result(f)
  !
  real(ark),intent(in) ::  x(:)
  real(ark),intent(in) ::  param(:)
  !
  real(ark)              :: f,coord(6)
  !
  integer(ik) nparams,ifst
  real(ark)   am1,am2,cs,y1,y2,y3,y4,y5,y6,V0,V1,V2,V3,V4,V5,V6,xieq(6),rho
  !
  real(ark) F000000
  real(ark) F100000,F200000,F300000,F400000,F500000,F600000,F010000,F020000,F030000,F040000,F050000,F060000, &
            F000100,F000200,F000300,F000400,F000500,F000600,F000001,F000002,F000003,F000004,F000005,F000006
  real(ark) F110000,F210000,F310000,F410000,F510000,F011000,F120000,F220000,F320000,F420000,F021000,F022000, &
            F130000,F230000,F330000,F031000,F032000,F033000,F140000,F240000,F041000,F042000,F150000,F051000, &
            F100100,F200100,F300100,F400100,F500100,F010100,F001100,F020100,F002100,F030100,F003100,F040100, &
            F004100,F050100,F005100,F000110,F100200,F200200,F300200,F400200,F010200,F001200,F020200,F002200, &
            F030200,F003200,F040200,F004200,F000210,F000220,F100300,F200300,F300300,F010300,F001300,F020300, &
            F002300,F030300,F003300,F000310,F000320,F000330,F100400,F200400,F010400,F001400,F020400,F002400, &
            F000410,F000420,F100500,F010500,F001500,F000510,F100001,F200001,F300001,F400001,F500001,F010001, &
            F020001,F030001,F040001,F050001,F000101,F000201,F000301,F000401,F000501,F100002,F200002,F300002, &
            F400002,F010002,F020002,F030002,F040002,F000102,F000202,F000302,F000402,F100003,F200003,F300003, &
            F010003,F020003,F030003,F000103,F000203,F000303,F100004,F200004,F010004,F020004,F000104,F000204, &
            F100005,F010005,F000105
  real(ark) F111000,F211000,F311000,F411000,F121000,F221000,F321000,F122000,F222000,F131000,F231000,F132000, &
            F141000,F110100,F101100,F210100,F201100,F310100,F301100,F410100,F401100,F011100,F120100,F102100, &
            F220100,F202100,F320100,F302100,F021100,F012100,F022100,F130100,F103100,F230100,F203100,F031100, &
            F013100,F032100,F023100,F140100,F104100,F041100,F014100,F100110,F200110,F300110,F400110,F010110, &
            F020110,F030110,F040110,F110200,F101200,F210200,F201200,F310200,F301200,F011200,F120200,F102200, &
            F220200,F202200,F021200,F012200,F022200,F130200,F103200,F031200,F013200,F100210,F200210,F300210, &
            F010210,F001210,F020210,F002210,F030210,F003210,F100220,F200220,F010220,F020220,F110300,F101300, &
            F210300,F201300,F011300,F120300,F102300,F021300,F012300,F100310,F200310,F010310,F001310,F020310, &
            F002310,F100320,F010320,F001320,F110400,F101400,F011400,F100410,F010410,F001410,F110001,F210001, &
            F310001,F410001,F011001,F120001,F220001,F320001,F021001,F022001,F130001,F230001,F031001,F032001, &
            F140001,F041001,F100101,F200101,F300101,F400101,F010101,F001101,F020101,F002101,F030101,F003101, &
            F040101,F004101,F000111,F100201,F200201,F300201,F010201,F001201,F020201,F002201,F030201,F003201, &
            F000211,F000221,F100301,F200301,F010301,F001301,F020301,F002301,F000311,F000321,F100401,F010401, &
            F001401,F000411,F110002,F210002,F310002,F011002,F120002,F220002,F021002,F022002,F130002,F031002, &
            F100102,F200102,F300102,F010102,F001102,F020102,F002102,F030102,F003102,F000112,F100202,F200202, &
            F010202,F001202,F020202,F002202,F000212,F000222,F100302,F010302,F001302,F000312,F110003,F210003, &
            F011003,F120003,F021003,F100103,F200103,F010103,F001103,F020103,F002103,F000113,F100203,F010203, &
            F001203,F000213,F110004,F011004,F100104,F010104,F001104,F000114
  real(ark) F111100,F211100,F311100,F121100,F112100,F221100,F212100,F122100,F131100,F113100,F110110,F210110, &
            F310110,F011110,F120110,F220110,F021110,F022110,F130110,F031110,F111200,F211200,F121200,F112200, &
            F110210,F101210,F210210,F201210,F011210,F120210,F102210,F021210,F012210,F110220,F011220,F111300, &
            F110310,F101310,F011310,F111001,F211001,F311001,F121001,F221001,F122001,F131001,F110101,F101101, &
            F210101,F201101,F310101,F301101,F011101,F120101,F102101,F220101,F202101,F021101,F012101,F022101, &
            F130101,F103101,F031101,F013101,F100111,F200111,F300111,F010111,F020111,F030111,F110201,F101201, &
            F210201,F201201,F011201,F120201,F102201,F021201,F012201,F100211,F200211,F010211,F001211,F020211, &
            F002211,F100221,F010221,F110301,F101301,F011301,F100311,F010311,F001311,F111002,F211002,F121002, &
            F110102,F101102,F210102,F201102,F011102,F120102,F102102,F021102,F012102,F100112,F200112,F010112, &
            F020112,F110202,F101202,F011202,F100212,F010212,F001212,F111003,F110103,F101103,F011103,F100113, &
            F010113
  real(ark) F111110,F211110,F121110,F111210,F111101,F211101,F121101,F112101,F110111,F210111,F011111,F120111, &
            F021111,F111201,F110211,F101211,F011211,F111102,F110112,F011112
  real(ark) F111111
   !
   coord = from_local_to_r1r2r3a1a2tau(x,6)
   !
   nparams=molec%parmax
   !
   rho  = coord(6)
   !
   cs=1.0_ark+cos(rho)
   !
   xieq(1:6) = ML_MEP_zxy2_rho_coeff(rho)
   !
   am1 = molec%specparam(1)
   am2 = molec%specparam(2)
   ! 
   y1=1.0_ark-exp(-am1*(coord(1)-xieq(1)))
   y2=1.0_ark-exp(-am2*(coord(2)-xieq(2)))
   y3=1.0_ark-exp(-am2*(coord(3)-xieq(3)))
   !
   y4=coord(4)-xieq(4)
   y5=coord(5)-xieq(5)
   y6=cs
   !
   V0=0.0_ark
   V1=0.0_ark
   V2=0.0_ark
   V3=0.0_ark
   V4=0.0_ark
   V5=0.0_ark
   V6=0.0_ark
   !
   ifst=0
   !
   if (nparams>=(1+ifst)) then
    !
    F000000 = param(1+ifst)
    V0=F000000
    !
   endif
   !
   if (nparams>=(25+ifst)) then
     !
     F100000 = param(2+ifst)
     F200000 = param(3+ifst)
     F300000 = param(4+ifst)
     F400000 = param(5+ifst)
     F500000 = param(6+ifst)
     F600000 = param(7+ifst)
     F010000 = param(8+ifst)
     F020000 = param(9+ifst)
     F030000 = param(10+ifst)
     F040000 = param(11+ifst)
     F050000 = param(12+ifst)
     F060000 = param(13+ifst)
     F000100 = param(14+ifst)
     F000200 = param(15+ifst)
     F000300 = param(16+ifst)
     F000400 = param(17+ifst)
     F000500 = param(18+ifst)
     F000600 = param(19+ifst)
     F000001 = param(20+ifst)
     F000002 = param(21+ifst)
     F000003 = param(22+ifst)
     F000004 = param(23+ifst)
     F000005 = param(24+ifst)
     F000006 = param(25+ifst)
     !
     V1=F100000*y1 + F200000*y1**2 + F300000*y1**3 + F400000*y1**4 + F500000*y1**5 + F600000*y1**6 + &
     F010000*(y2 + y3) + F020000*(y2**2 + y3**2) + F030000*(y2**3 + y3**3) + F040000*(y2**4 + y3**4) + &
     F050000*(y2**5 + y3**5) + F060000*(y2**6 + y3**6) + F000100*(y4 + y5) + F000200*(y4**2 + y5**2) + &
     F000300*(y4**3 + y5**3) + F000400*(y4**4 + y5**4) + F000500*(y4**5 + y5**5) + F000600*(y4**6 + y5**6) + &
     F000001*y6 + F000002*y6**2 + F000003*y6**3 + F000004*y6**4 + F000005*y6**5 + F000006*y6**6
     !
   endif
   !
   if (nparams>=(148+ifst)) then
    !
     F110000 = param(26+ifst)
     F210000 = param(27+ifst)
     F310000 = param(28+ifst)
     F410000 = param(29+ifst)
     F510000 = param(30+ifst)
     F011000 = param(31+ifst)
     F120000 = param(32+ifst)
     F220000 = param(33+ifst)
     F320000 = param(34+ifst)
     F420000 = param(35+ifst)
     F021000 = param(36+ifst)
     F022000 = param(37+ifst)
     F130000 = param(38+ifst)
     F230000 = param(39+ifst)
     F330000 = param(40+ifst)
     F031000 = param(41+ifst)
     F032000 = param(42+ifst)
     F033000 = param(43+ifst)
     F140000 = param(44+ifst)
     F240000 = param(45+ifst)
     F041000 = param(46+ifst)
     F042000 = param(47+ifst)
     F150000 = param(48+ifst)
     F051000 = param(49+ifst)
     F100100 = param(50+ifst)
     F200100 = param(51+ifst)
     F300100 = param(52+ifst)
     F400100 = param(53+ifst)
     F500100 = param(54+ifst)
     F010100 = param(55+ifst)
     F001100 = param(56+ifst)
     F020100 = param(57+ifst)
     F002100 = param(58+ifst)
     F030100 = param(59+ifst)
     F003100 = param(60+ifst)
     F040100 = param(61+ifst)
     F004100 = param(62+ifst)
     F050100 = param(63+ifst)
     F005100 = param(64+ifst)
     F000110 = param(65+ifst)
     F100200 = param(66+ifst)
     F200200 = param(67+ifst)
     F300200 = param(68+ifst)
     F400200 = param(69+ifst)
     F010200 = param(70+ifst)
     F001200 = param(71+ifst)
     F020200 = param(72+ifst)
     F002200 = param(73+ifst)
     F030200 = param(74+ifst)
     F003200 = param(75+ifst)
     F040200 = param(76+ifst)
     F004200 = param(77+ifst)
     F000210 = param(78+ifst)
     F000220 = param(79+ifst)
     F100300 = param(80+ifst)
     F200300 = param(81+ifst)
     F300300 = param(82+ifst)
     F010300 = param(83+ifst)
     F001300 = param(84+ifst)
     F020300 = param(85+ifst)
     F002300 = param(86+ifst)
     F030300 = param(87+ifst)
     F003300 = param(88+ifst)
     F000310 = param(89+ifst)
     F000320 = param(90+ifst)
     F000330 = param(91+ifst)
     F100400 = param(92+ifst)
     F200400 = param(93+ifst)
     F010400 = param(94+ifst)
     F001400 = param(95+ifst)
     F020400 = param(96+ifst)
     F002400 = param(97+ifst)
     F000410 = param(98+ifst)
     F000420 = param(99+ifst)
     F100500 = param(100+ifst)
     F010500 = param(101+ifst)
     F001500 = param(102+ifst)
     F000510 = param(103+ifst)
     F100001 = param(104+ifst)
     F200001 = param(105+ifst)
     F300001 = param(106+ifst)
     F400001 = param(107+ifst)
     F500001 = param(108+ifst)
     F010001 = param(109+ifst)
     F020001 = param(110+ifst)
     F030001 = param(111+ifst)
     F040001 = param(112+ifst)
     F050001 = param(113+ifst)
     F000101 = param(114+ifst)
     F000201 = param(115+ifst)
     F000301 = param(116+ifst)
     F000401 = param(117+ifst)
     F000501 = param(118+ifst)
     F100002 = param(119+ifst)
     F200002 = param(120+ifst)
     F300002 = param(121+ifst)
     F400002 = param(122+ifst)
     F010002 = param(123+ifst)
     F020002 = param(124+ifst)
     F030002 = param(125+ifst)
     F040002 = param(126+ifst)
     F000102 = param(127+ifst)
     F000202 = param(128+ifst)
     F000302 = param(129+ifst)
     F000402 = param(130+ifst)
     F100003 = param(131+ifst)
     F200003 = param(132+ifst)
     F300003 = param(133+ifst)
     F010003 = param(134+ifst)
     F020003 = param(135+ifst)
     F030003 = param(136+ifst)
     F000103 = param(137+ifst)
     F000203 = param(138+ifst)
     F000303 = param(139+ifst)
     F100004 = param(140+ifst)
     F200004 = param(141+ifst)
     F010004 = param(142+ifst)
     F020004 = param(143+ifst)
     F000104 = param(144+ifst)
     F000204 = param(145+ifst)
     F100005 = param(146+ifst)
     F010005 = param(147+ifst)
     F000105 = param(148+ifst)
     !
     V2=F011000*y2*y3 + F022000*y2**2*y3**2 + F033000*y2**3*y3**3 + F110000*y1*(y2 + y3) + &
     F210000*y1**2*(y2 + y3) + F310000*y1**3*(y2 + y3) + F410000*y1**4*(y2 + y3) + F510000*y1**5*(y2 + y3) + &
     F120000*y1*(y2**2 + y3**2) + F220000*y1**2*(y2**2 + y3**2) + F320000*y1**3*(y2**2 + y3**2) + &
     F420000*y1**4*(y2**2 + y3**2) + F021000*(y2**2*y3 + y2*y3**2) + F130000*y1*(y2**3 + y3**3) + &
     F230000*y1**2*(y2**3 + y3**3) + F330000*y1**3*(y2**3 + y3**3) + F031000*(y2**3*y3 + y2*y3**3) + &
     F032000*(y2**3*y3**2 + y2**2*y3**3) + F140000*y1*(y2**4 + y3**4) + F240000*y1**2*(y2**4 + y3**4) + &
     F041000*(y2**4*y3 + y2*y3**4) + F042000*(y2**4*y3**2 + y2**2*y3**4) + F150000*y1*(y2**5 + y3**5) + &
     F051000*(y2**5*y3 + y2*y3**5) + F000110*y4*y5 + F000220*y4**2*y5**2 + F000330*y4**3*y5**3 + &
     F100100*y1*(y4 + y5) + F200100*y1**2*(y4 + y5) + F300100*y1**3*(y4 + y5) + F400100*y1**4*(y4 + y5) + &
     F500100*y1**5*(y4 + y5) + F001100*(y3*y4 + y2*y5) + F002100*(y3**2*y4 + y2**2*y5) + &
     F003100*(y3**3*y4 + y2**3*y5) + F004100*(y3**4*y4 + y2**4*y5) + F005100*(y3**5*y4 + y2**5*y5) + &
     F010100*(y2*y4 + y3*y5) + F020100*(y2**2*y4 + y3**2*y5) + F030100*(y2**3*y4 + y3**3*y5) + &
     F040100*(y2**4*y4 + y3**4*y5) + F050100*(y2**5*y4 + y3**5*y5) + F100200*y1*(y4**2 + y5**2) + &
     F200200*y1**2*(y4**2 + y5**2) + F300200*y1**3*(y4**2 + y5**2) + F400200*y1**4*(y4**2 + y5**2) + &
     F001200*(y3*y4**2 + y2*y5**2) + F002200*(y3**2*y4**2 + y2**2*y5**2) + F003200*(y3**3*y4**2 + y2**3*y5**2) + &
     F004200*(y3**4*y4**2 + y2**4*y5**2) + F010200*(y2*y4**2 + y3*y5**2) + F020200*(y2**2*y4**2 + y3**2*y5**2) + &
     F030200*(y2**3*y4**2 + y3**3*y5**2) + F040200*(y2**4*y4**2 + y3**4*y5**2) + F000210*(y4**2*y5 + y4*y5**2) + &
     F100300*y1*(y4**3 + y5**3) + F200300*y1**2*(y4**3 + y5**3) + F300300*y1**3*(y4**3 + y5**3) + &
     F001300*(y3*y4**3 + y2*y5**3) + F002300*(y3**2*y4**3 + y2**2*y5**3) + F003300*(y3**3*y4**3 + y2**3*y5**3) + &
     F010300*(y2*y4**3 + y3*y5**3) + F020300*(y2**2*y4**3 + y3**2*y5**3) + F030300*(y2**3*y4**3 + y3**3*y5**3) + &
     F000310*(y4**3*y5 + y4*y5**3) + F000320*(y4**3*y5**2 + y4**2*y5**3) + F100400*y1*(y4**4 + y5**4) + &
     F200400*y1**2*(y4**4 + y5**4) + F001400*(y3*y4**4 + y2*y5**4) + F002400*(y3**2*y4**4 + y2**2*y5**4) + &
     F010400*(y2*y4**4 + y3*y5**4) + F020400*(y2**2*y4**4 + y3**2*y5**4) + F000410*(y4**4*y5 + y4*y5**4) + &
     F000420*(y4**4*y5**2 + y4**2*y5**4) + F100500*y1*(y4**5 + y5**5) + F001500*(y3*y4**5 + y2*y5**5) + &
     F010500*(y2*y4**5 + y3*y5**5) + F000510*(y4**5*y5 + y4*y5**5) + F100001*y1*y6 + F200001*y1**2*y6 + &
     F300001*y1**3*y6 + F400001*y1**4*y6 + F500001*y1**5*y6 + F010001*(y2 + y3)*y6 + F020001*(y2**2 + y3**2)*y6 + &
     F030001*(y2**3 + y3**3)*y6 + F040001*(y2**4 + y3**4)*y6 + F050001*(y2**5 + y3**5)*y6 + F000101*(y4 + y5)*y6 + &
     F000201*(y4**2 + y5**2)*y6 + F000301*(y4**3 + y5**3)*y6 + F000401*(y4**4 + y5**4)*y6 + F000501*(y4**5 + y5**5)*y6 + &
     F100002*y1*y6**2 + F200002*y1**2*y6**2 + F300002*y1**3*y6**2 + F400002*y1**4*y6**2 + F010002*(y2 + y3)*y6**2 + &
     F020002*(y2**2 + y3**2)*y6**2 + F030002*(y2**3 + y3**3)*y6**2 + F040002*(y2**4 + y3**4)*y6**2 + &
     F000102*(y4 + y5)*y6**2 + F000202*(y4**2 + y5**2)*y6**2 + F000302*(y4**3 + y5**3)*y6**2 + &
     F000402*(y4**4 + y5**4)*y6**2 + F100003*y1*y6**3 + F200003*y1**2*y6**3 + F300003*y1**3*y6**3 + &
     F010003*(y2 + y3)*y6**3 + F020003*(y2**2 + y3**2)*y6**3 + F030003*(y2**3 + y3**3)*y6**3 + F000103*(y4 + y5)*y6**3 + &
     F000203*(y4**2 + y5**2)*y6**3 + F000303*(y4**3 + y5**3)*y6**3 + F100004*y1*y6**4 + F200004*y1**2*y6**4 + &
     F010004*(y2 + y3)*y6**4 + F020004*(y2**2 + y3**2)*y6**4 + F000104*(y4 + y5)*y6**4 + F000204*(y4**2 + y5**2)*y6**4 + &
     F100005*y1*y6**5 + F010005*(y2 + y3)*y6**5 + F000105*(y4 + y5)*y6**5
     !
   endif
   if (nparams>=(360+ifst)) then
     !
     F111000 = param(149+ifst)
     F211000 = param(150+ifst)
     F311000 = param(151+ifst)
     F411000 = param(152+ifst)
     F121000 = param(153+ifst)
     F221000 = param(154+ifst)
     F321000 = param(155+ifst)
     F122000 = param(156+ifst)
     F222000 = param(157+ifst)
     F131000 = param(158+ifst)
     F231000 = param(159+ifst)
     F132000 = param(160+ifst)
     F141000 = param(161+ifst)
     F110100 = param(162+ifst)
     F101100 = param(163+ifst)
     F210100 = param(164+ifst)
     F201100 = param(165+ifst)
     F310100 = param(166+ifst)
     F301100 = param(167+ifst)
     F410100 = param(168+ifst)
     F401100 = param(169+ifst)
     F011100 = param(170+ifst)
     F120100 = param(171+ifst)
     F102100 = param(172+ifst)
     F220100 = param(173+ifst)
     F202100 = param(174+ifst)
     F320100 = param(175+ifst)
     F302100 = param(176+ifst)
     F021100 = param(177+ifst)
     F012100 = param(178+ifst)
     F022100 = param(179+ifst)
     F130100 = param(180+ifst)
     F103100 = param(181+ifst)
     F230100 = param(182+ifst)
     F203100 = param(183+ifst)
     F031100 = param(184+ifst)
     F013100 = param(185+ifst)
     F032100 = param(186+ifst)
     F023100 = param(187+ifst)
     F140100 = param(188+ifst)
     F104100 = param(189+ifst)
     F041100 = param(190+ifst)
     F014100 = param(191+ifst)
     F100110 = param(192+ifst)
     F200110 = param(193+ifst)
     F300110 = param(194+ifst)
     F400110 = param(195+ifst)
     F010110 = param(196+ifst)
     F020110 = param(197+ifst)
     F030110 = param(198+ifst)
     F040110 = param(199+ifst)
     F110200 = param(200+ifst)
     F101200 = param(201+ifst)
     F210200 = param(202+ifst)
     F201200 = param(203+ifst)
     F310200 = param(204+ifst)
     F301200 = param(205+ifst)
     F011200 = param(206+ifst)
     F120200 = param(207+ifst)
     F102200 = param(208+ifst)
     F220200 = param(209+ifst)
     F202200 = param(210+ifst)
     F021200 = param(211+ifst)
     F012200 = param(212+ifst)
     F022200 = param(213+ifst)
     F130200 = param(214+ifst)
     F103200 = param(215+ifst)
     F031200 = param(216+ifst)
     F013200 = param(217+ifst)
     F100210 = param(218+ifst)
     F200210 = param(219+ifst)
     F300210 = param(220+ifst)
     F010210 = param(221+ifst)
     F001210 = param(222+ifst)
     F020210 = param(223+ifst)
     F002210 = param(224+ifst)
     F030210 = param(225+ifst)
     F003210 = param(226+ifst)
     F100220 = param(227+ifst)
     F200220 = param(228+ifst)
     F010220 = param(229+ifst)
     F020220 = param(230+ifst)
     F110300 = param(231+ifst)
     F101300 = param(232+ifst)
     F210300 = param(233+ifst)
     F201300 = param(234+ifst)
     F011300 = param(235+ifst)
     F120300 = param(236+ifst)
     F102300 = param(237+ifst)
     F021300 = param(238+ifst)
     F012300 = param(239+ifst)
     F100310 = param(240+ifst)
     F200310 = param(241+ifst)
     F010310 = param(242+ifst)
     F001310 = param(243+ifst)
     F020310 = param(244+ifst)
     F002310 = param(245+ifst)
     F100320 = param(246+ifst)
     F010320 = param(247+ifst)
     F001320 = param(248+ifst)
     F110400 = param(249+ifst)
     F101400 = param(250+ifst)
     F011400 = param(251+ifst)
     F100410 = param(252+ifst)
     F010410 = param(253+ifst)
     F001410 = param(254+ifst)
     F110001 = param(255+ifst)
     F210001 = param(256+ifst)
     F310001 = param(257+ifst)
     F410001 = param(258+ifst)
     F011001 = param(259+ifst)
     F120001 = param(260+ifst)
     F220001 = param(261+ifst)
     F320001 = param(262+ifst)
     F021001 = param(263+ifst)
     F022001 = param(264+ifst)
     F130001 = param(265+ifst)
     F230001 = param(266+ifst)
     F031001 = param(267+ifst)
     F032001 = param(268+ifst)
     F140001 = param(269+ifst)
     F041001 = param(270+ifst)
     F100101 = param(271+ifst)
     F200101 = param(272+ifst)
     F300101 = param(273+ifst)
     F400101 = param(274+ifst)
     F010101 = param(275+ifst)
     F001101 = param(276+ifst)
     F020101 = param(277+ifst)
     F002101 = param(278+ifst)
     F030101 = param(279+ifst)
     F003101 = param(280+ifst)
     F040101 = param(281+ifst)
     F004101 = param(282+ifst)
     F000111 = param(283+ifst)
     F100201 = param(284+ifst)
     F200201 = param(285+ifst)
     F300201 = param(286+ifst)
     F010201 = param(287+ifst)
     F001201 = param(288+ifst)
     F020201 = param(289+ifst)
     F002201 = param(290+ifst)
     F030201 = param(291+ifst)
     F003201 = param(292+ifst)
     F000211 = param(293+ifst)
     F000221 = param(294+ifst)
     F100301 = param(295+ifst)
     F200301 = param(296+ifst)
     F010301 = param(297+ifst)
     F001301 = param(298+ifst)
     F020301 = param(299+ifst)
     F002301 = param(300+ifst)
     F000311 = param(301+ifst)
     F000321 = param(302+ifst)
     F100401 = param(303+ifst)
     F010401 = param(304+ifst)
     F001401 = param(305+ifst)
     F000411 = param(306+ifst)
     F110002 = param(307+ifst)
     F210002 = param(308+ifst)
     F310002 = param(309+ifst)
     F011002 = param(310+ifst)
     F120002 = param(311+ifst)
     F220002 = param(312+ifst)
     F021002 = param(313+ifst)
     F022002 = param(314+ifst)
     F130002 = param(315+ifst)
     F031002 = param(316+ifst)
     F100102 = param(317+ifst)
     F200102 = param(318+ifst)
     F300102 = param(319+ifst)
     F010102 = param(320+ifst)
     F001102 = param(321+ifst)
     F020102 = param(322+ifst)
     F002102 = param(323+ifst)
     F030102 = param(324+ifst)
     F003102 = param(325+ifst)
     F000112 = param(326+ifst)
     F100202 = param(327+ifst)
     F200202 = param(328+ifst)
     F010202 = param(329+ifst)
     F001202 = param(330+ifst)
     F020202 = param(331+ifst)
     F002202 = param(332+ifst)
     F000212 = param(333+ifst)
     F000222 = param(334+ifst)
     F100302 = param(335+ifst)
     F010302 = param(336+ifst)
     F001302 = param(337+ifst)
     F000312 = param(338+ifst)
     F110003 = param(339+ifst)
     F210003 = param(340+ifst)
     F011003 = param(341+ifst)
     F120003 = param(342+ifst)
     F021003 = param(343+ifst)
     F100103 = param(344+ifst)
     F200103 = param(345+ifst)
     F010103 = param(346+ifst)
     F001103 = param(347+ifst)
     F020103 = param(348+ifst)
     F002103 = param(349+ifst)
     F000113 = param(350+ifst)
     F100203 = param(351+ifst)
     F010203 = param(352+ifst)
     F001203 = param(353+ifst)
     F000213 = param(354+ifst)
     F110004 = param(355+ifst)
     F011004 = param(356+ifst)
     F100104 = param(357+ifst)
     F010104 = param(358+ifst)
     F001104 = param(359+ifst)
     F000114 = param(360+ifst)
     !
     V3=F111000*y1*y2*y3 + F211000*y1**2*y2*y3 + F311000*y1**3*y2*y3 + F411000*y1**4*y2*y3 + F122000*y1*y2**2*y3**2 + &
     F222000*y1**2*y2**2*y3**2 + F121000*y1*(y2**2*y3 + y2*y3**2) + F221000*y1**2*(y2**2*y3 + y2*y3**2) + &
     F321000*y1**3*(y2**2*y3 + y2*y3**2) + F131000*y1*(y2**3*y3 + y2*y3**3) + F231000*y1**2*(y2**3*y3 + y2*y3**3) + &
     F132000*y1*(y2**3*y3**2 + y2**2*y3**3) + F141000*y1*(y2**4*y3 + y2*y3**4) + F100110*y1*y4*y5 + &
     F200110*y1**2*y4*y5 + F300110*y1**3*y4*y5 + F400110*y1**4*y4*y5 + F010110*(y2 + y3)*y4*y5 + &
     F020110*(y2**2 + y3**2)*y4*y5 + F030110*(y2**3 + y3**3)*y4*y5 + F040110*(y2**4 + y3**4)*y4*y5 + &
     F100220*y1*y4**2*y5**2 + F200220*y1**2*y4**2*y5**2 + F010220*(y2 + y3)*y4**2*y5**2 + &
     F020220*(y2**2 + y3**2)*y4**2*y5**2 + F011100*y2*y3*(y4 + y5) + F022100*y2**2*y3**2*(y4 + y5) + &
     F101100*y1*(y3*y4 + y2*y5) + F201100*y1**2*(y3*y4 + y2*y5) + F301100*y1**3*(y3*y4 + y2*y5) + &
     F401100*y1**4*(y3*y4 + y2*y5) + F102100*y1*(y3**2*y4 + y2**2*y5) + F202100*y1**2*(y3**2*y4 + y2**2*y5) + &
     F302100*y1**3*(y3**2*y4 + y2**2*y5) + F103100*y1*(y3**3*y4 + y2**3*y5) + F203100*y1**2*(y3**3*y4 + y2**3*y5) + &
     F104100*y1*(y3**4*y4 + y2**4*y5) + F110100*y1*(y2*y4 + y3*y5) + F210100*y1**2*(y2*y4 + y3*y5) + &
     F310100*y1**3*(y2*y4 + y3*y5) + F410100*y1**4*(y2*y4 + y3*y5) + F012100*(y2*y3**2*y4 + y2**2*y3*y5) + &
     F013100*(y2*y3**3*y4 + y2**3*y3*y5) + F014100*(y2*y3**4*y4 + y2**4*y3*y5) + F120100*y1*(y2**2*y4 + y3**2*y5) + &
     F220100*y1**2*(y2**2*y4 + y3**2*y5) + F320100*y1**3*(y2**2*y4 + y3**2*y5) + F021100*(y2**2*y3*y4 + y2*y3**2*y5) + &
     F023100*(y2**2*y3**3*y4 + y2**3*y3**2*y5) + F130100*y1*(y2**3*y4 + y3**3*y5) + F230100*y1**2*(y2**3*y4 + y3**3*y5) + &
     F031100*(y2**3*y3*y4 + y2*y3**3*y5) + F032100*(y2**3*y3**2*y4 + y2**2*y3**3*y5) + F140100*y1*(y2**4*y4 + y3**4*y5) + &
     F041100*(y2**4*y3*y4 + y2*y3**4*y5) + F011200*y2*y3*(y4**2 + y5**2) + F022200*y2**2*y3**2*(y4**2 + y5**2) + &
     F101200*y1*(y3*y4**2 + y2*y5**2) + F201200*y1**2*(y3*y4**2 + y2*y5**2) + F301200*y1**3*(y3*y4**2 + y2*y5**2) + &
     F102200*y1*(y3**2*y4**2 + y2**2*y5**2) + F202200*y1**2*(y3**2*y4**2 + y2**2*y5**2) + &
     F103200*y1*(y3**3*y4**2 + y2**3*y5**2) + F110200*y1*(y2*y4**2 + y3*y5**2) + F210200*y1**2*(y2*y4**2 + y3*y5**2) + &
     F310200*y1**3*(y2*y4**2 + y3*y5**2) + F012200*(y2*y3**2*y4**2 + y2**2*y3*y5**2) + &
     F013200*(y2*y3**3*y4**2 + y2**3*y3*y5**2) + F120200*y1*(y2**2*y4**2 + y3**2*y5**2) + &
     F220200*y1**2*(y2**2*y4**2 + y3**2*y5**2) + F021200*(y2**2*y3*y4**2 + y2*y3**2*y5**2) + &
     F130200*y1*(y2**3*y4**2 + y3**3*y5**2) + F031200*(y2**3*y3*y4**2 + y2*y3**3*y5**2) + &
     F100210*y1*(y4**2*y5 + y4*y5**2) + F200210*y1**2*(y4**2*y5 + y4*y5**2) + F300210*y1**3*(y4**2*y5 + y4*y5**2) + &
     F001210*(y3*y4**2*y5 + y2*y4*y5**2) + F002210*(y3**2*y4**2*y5 + y2**2*y4*y5**2) + &
     F003210*(y3**3*y4**2*y5 + y2**3*y4*y5**2) + F010210*(y2*y4**2*y5 + y3*y4*y5**2) + &
     F020210*(y2**2*y4**2*y5 + y3**2*y4*y5**2) + F030210*(y2**3*y4**2*y5 + y3**3*y4*y5**2) + &
     F011300*y2*y3*(y4**3 + y5**3) + F101300*y1*(y3*y4**3 + y2*y5**3) + F201300*y1**2*(y3*y4**3 + y2*y5**3) + &
     F102300*y1*(y3**2*y4**3 + y2**2*y5**3) + F110300*y1*(y2*y4**3 + y3*y5**3) + F210300*y1**2*(y2*y4**3 + y3*y5**3) + &
     F012300*(y2*y3**2*y4**3 + y2**2*y3*y5**3) + F120300*y1*(y2**2*y4**3 + y3**2*y5**3) + &
     F021300*(y2**2*y3*y4**3 + y2*y3**2*y5**3) + F100310*y1*(y4**3*y5 + y4*y5**3) + F200310*y1**2*(y4**3*y5 + y4*y5**3) + &
     F001310*(y3*y4**3*y5 + y2*y4*y5**3) + F002310*(y3**2*y4**3*y5 + y2**2*y4*y5**3) + F010310*(y2*y4**3*y5 + y3*y4*y5**3) + &
     F020310*(y2**2*y4**3*y5 + y3**2*y4*y5**3) + F100320*y1*(y4**3*y5**2 + y4**2*y5**3) + &
     F001320*(y3*y4**3*y5**2 + y2*y4**2*y5**3) + F010320*(y2*y4**3*y5**2 + y3*y4**2*y5**3) + &
     F011400*y2*y3*(y4**4 + y5**4) + F101400*y1*(y3*y4**4 + y2*y5**4) + F110400*y1*(y2*y4**4 + y3*y5**4) + &
     F100410*y1*(y4**4*y5 + y4*y5**4) + F001410*(y3*y4**4*y5 + y2*y4*y5**4) + F010410*(y2*y4**4*y5 + y3*y4*y5**4) + &
     F011001*y2*y3*y6 + F022001*y2**2*y3**2*y6 + F110001*y1*(y2 + y3)*y6 + F210001*y1**2*(y2 + y3)*y6 + &
     F310001*y1**3*(y2 + y3)*y6 + F410001*y1**4*(y2 + y3)*y6 + F120001*y1*(y2**2 + y3**2)*y6 + F220001*y1**2*(y2**2 + y3**2)*y6 + &
     F320001*y1**3*(y2**2 + y3**2)*y6 + F021001*(y2**2*y3 + y2*y3**2)*y6 + F130001*y1*(y2**3 + y3**3)*y6 + &
     F230001*y1**2*(y2**3 + y3**3)*y6 + F031001*(y2**3*y3 + y2*y3**3)*y6 + F032001*(y2**3*y3**2 + y2**2*y3**3)*y6 + &
     F140001*y1*(y2**4 + y3**4)*y6 + F041001*(y2**4*y3 + y2*y3**4)*y6 + F000111*y4*y5*y6 + F000221*y4**2*y5**2*y6 + &
     F100101*y1*(y4 + y5)*y6 + F200101*y1**2*(y4 + y5)*y6 + F300101*y1**3*(y4 + y5)*y6 + F400101*y1**4*(y4 + y5)*y6 + &
     F001101*(y3*y4 + y2*y5)*y6 + F002101*(y3**2*y4 + y2**2*y5)*y6 + F003101*(y3**3*y4 + y2**3*y5)*y6 + &
     F004101*(y3**4*y4 + y2**4*y5)*y6 + F010101*(y2*y4 + y3*y5)*y6 + F020101*(y2**2*y4 + y3**2*y5)*y6 + &
     F030101*(y2**3*y4 + y3**3*y5)*y6 + F040101*(y2**4*y4 + y3**4*y5)*y6 + F100201*y1*(y4**2 + y5**2)*y6 + &
     F200201*y1**2*(y4**2 + y5**2)*y6 + F300201*y1**3*(y4**2 + y5**2)*y6 + F001201*(y3*y4**2 + y2*y5**2)*y6 + &
     F002201*(y3**2*y4**2 + y2**2*y5**2)*y6 + F003201*(y3**3*y4**2 + y2**3*y5**2)*y6 + F010201*(y2*y4**2 + y3*y5**2)*y6 + &
     F020201*(y2**2*y4**2 + y3**2*y5**2)*y6 + F030201*(y2**3*y4**2 + y3**3*y5**2)*y6 + F000211*(y4**2*y5 + y4*y5**2)*y6 + &
     F100301*y1*(y4**3 + y5**3)*y6 + F200301*y1**2*(y4**3 + y5**3)*y6 + F001301*(y3*y4**3 + y2*y5**3)*y6 + &
     F002301*(y3**2*y4**3 + y2**2*y5**3)*y6 + F010301*(y2*y4**3 + y3*y5**3)*y6 + F020301*(y2**2*y4**3 + y3**2*y5**3)*y6 + &
     F000311*(y4**3*y5 + y4*y5**3)*y6 + F000321*(y4**3*y5**2 + y4**2*y5**3)*y6 + F100401*y1*(y4**4 + y5**4)*y6 + &
     F001401*(y3*y4**4 + y2*y5**4)*y6 + F010401*(y2*y4**4 + y3*y5**4)*y6 + F000411*(y4**4*y5 + y4*y5**4)*y6 + &
     F011002*y2*y3*y6**2 + F022002*y2**2*y3**2*y6**2 + F110002*y1*(y2 + y3)*y6**2 + F210002*y1**2*(y2 + y3)*y6**2 + &
     F310002*y1**3*(y2 + y3)*y6**2 + F120002*y1*(y2**2 + y3**2)*y6**2 + F220002*y1**2*(y2**2 + y3**2)*y6**2 + &
     F021002*(y2**2*y3 + y2*y3**2)*y6**2 + F130002*y1*(y2**3 + y3**3)*y6**2 + F031002*(y2**3*y3 + y2*y3**3)*y6**2 + &
     F000112*y4*y5*y6**2 + F000222*y4**2*y5**2*y6**2 + F100102*y1*(y4 + y5)*y6**2 + F200102*y1**2*(y4 + y5)*y6**2 + &
     F300102*y1**3*(y4 + y5)*y6**2 + F001102*(y3*y4 + y2*y5)*y6**2 + F002102*(y3**2*y4 + y2**2*y5)*y6**2 + &
     F003102*(y3**3*y4 + y2**3*y5)*y6**2 + F010102*(y2*y4 + y3*y5)*y6**2 + F020102*(y2**2*y4 + y3**2*y5)*y6**2 + &
     F030102*(y2**3*y4 + y3**3*y5)*y6**2 + F100202*y1*(y4**2 + y5**2)*y6**2 + F200202*y1**2*(y4**2 + y5**2)*y6**2 + &
     F001202*(y3*y4**2 + y2*y5**2)*y6**2 + F002202*(y3**2*y4**2 + y2**2*y5**2)*y6**2 + F010202*(y2*y4**2 + y3*y5**2)*y6**2 + &
     F020202*(y2**2*y4**2 + y3**2*y5**2)*y6**2 + F000212*(y4**2*y5 + y4*y5**2)*y6**2 + F100302*y1*(y4**3 + y5**3)*y6**2 + &
     F001302*(y3*y4**3 + y2*y5**3)*y6**2 + F010302*(y2*y4**3 + y3*y5**3)*y6**2 + F000312*(y4**3*y5 + y4*y5**3)*y6**2 + &
     F011003*y2*y3*y6**3 + F110003*y1*(y2 + y3)*y6**3 + F210003*y1**2*(y2 + y3)*y6**3 + F120003*y1*(y2**2 + y3**2)*y6**3 + &
     F021003*(y2**2*y3 + y2*y3**2)*y6**3 + F000113*y4*y5*y6**3 + F100103*y1*(y4 + y5)*y6**3 + F200103*y1**2*(y4 + y5)*y6**3 + &
     F001103*(y3*y4 + y2*y5)*y6**3 + F002103*(y3**2*y4 + y2**2*y5)*y6**3 + F010103*(y2*y4 + y3*y5)*y6**3 + &
     F020103*(y2**2*y4 + y3**2*y5)*y6**3 + F100203*y1*(y4**2 + y5**2)*y6**3 + F001203*(y3*y4**2 + y2*y5**2)*y6**3 + &
     F010203*(y2*y4**2 + y3*y5**2)*y6**3 + F000213*(y4**2*y5 + y4*y5**2)*y6**3 + F011004*y2*y3*y6**4 + &
     F110004*y1*(y2 + y3)*y6**4 + F000114*y4*y5*y6**4 + F100104*y1*(y4 + y5)*y6**4 + F001104*(y3*y4 + y2*y5)*y6**4 + &
     F010104*(y2*y4 + y3*y5)*y6**4
     !
   endif
   if (nparams>=(481+ifst)) then
     !
     F111100 = param(361+ifst)
     F211100 = param(362+ifst)
     F311100 = param(363+ifst)
     F121100 = param(364+ifst)
     F112100 = param(365+ifst)
     F221100 = param(366+ifst)
     F212100 = param(367+ifst)
     F122100 = param(368+ifst)
     F131100 = param(369+ifst)
     F113100 = param(370+ifst)
     F110110 = param(371+ifst)
     F210110 = param(372+ifst)
     F310110 = param(373+ifst)
     F011110 = param(374+ifst)
     F120110 = param(375+ifst)
     F220110 = param(376+ifst)
     F021110 = param(377+ifst)
     F022110 = param(378+ifst)
     F130110 = param(379+ifst)
     F031110 = param(380+ifst)
     F111200 = param(381+ifst)
     F211200 = param(382+ifst)
     F121200 = param(383+ifst)
     F112200 = param(384+ifst)
     F110210 = param(385+ifst)
     F101210 = param(386+ifst)
     F210210 = param(387+ifst)
     F201210 = param(388+ifst)
     F011210 = param(389+ifst)
     F120210 = param(390+ifst)
     F102210 = param(391+ifst)
     F021210 = param(392+ifst)
     F012210 = param(393+ifst)
     F110220 = param(394+ifst)
     F011220 = param(395+ifst)
     F111300 = param(396+ifst)
     F110310 = param(397+ifst)
     F101310 = param(398+ifst)
     F011310 = param(399+ifst)
     F111001 = param(400+ifst)
     F211001 = param(401+ifst)
     F311001 = param(402+ifst)
     F121001 = param(403+ifst)
     F221001 = param(404+ifst)
     F122001 = param(405+ifst)
     F131001 = param(406+ifst)
     F110101 = param(407+ifst)
     F101101 = param(408+ifst)
     F210101 = param(409+ifst)
     F201101 = param(410+ifst)
     F310101 = param(411+ifst)
     F301101 = param(412+ifst)
     F011101 = param(413+ifst)
     F120101 = param(414+ifst)
     F102101 = param(415+ifst)
     F220101 = param(416+ifst)
     F202101 = param(417+ifst)
     F021101 = param(418+ifst)
     F012101 = param(419+ifst)
     F022101 = param(420+ifst)
     F130101 = param(421+ifst)
     F103101 = param(422+ifst)
     F031101 = param(423+ifst)
     F013101 = param(424+ifst)
     F100111 = param(425+ifst)
     F200111 = param(426+ifst)
     F300111 = param(427+ifst)
     F010111 = param(428+ifst)
     F020111 = param(429+ifst)
     F030111 = param(430+ifst)
     F110201 = param(431+ifst)
     F101201 = param(432+ifst)
     F210201 = param(433+ifst)
     F201201 = param(434+ifst)
     F011201 = param(435+ifst)
     F120201 = param(436+ifst)
     F102201 = param(437+ifst)
     F021201 = param(438+ifst)
     F012201 = param(439+ifst)
     F100211 = param(440+ifst)
     F200211 = param(441+ifst)
     F010211 = param(442+ifst)
     F001211 = param(443+ifst)
     F020211 = param(444+ifst)
     F002211 = param(445+ifst)
     F100221 = param(446+ifst)
     F010221 = param(447+ifst)
     F110301 = param(448+ifst)
     F101301 = param(449+ifst)
     F011301 = param(450+ifst)
     F100311 = param(451+ifst)
     F010311 = param(452+ifst)
     F001311 = param(453+ifst)
     F111002 = param(454+ifst)
     F211002 = param(455+ifst)
     F121002 = param(456+ifst)
     F110102 = param(457+ifst)
     F101102 = param(458+ifst)
     F210102 = param(459+ifst)
     F201102 = param(460+ifst)
     F011102 = param(461+ifst)
     F120102 = param(462+ifst)
     F102102 = param(463+ifst)
     F021102 = param(464+ifst)
     F012102 = param(465+ifst)
     F100112 = param(466+ifst)
     F200112 = param(467+ifst)
     F010112 = param(468+ifst)
     F020112 = param(469+ifst)
     F110202 = param(470+ifst)
     F101202 = param(471+ifst)
     F011202 = param(472+ifst)
     F100212 = param(473+ifst)
     F010212 = param(474+ifst)
     F001212 = param(475+ifst)
     F111003 = param(476+ifst)
     F110103 = param(477+ifst)
     F101103 = param(478+ifst)
     F011103 = param(479+ifst)
     F100113 = param(480+ifst)
     F010113 = param(481+ifst)
     !
     V4=F011110*y2*y3*y4*y5 + F022110*y2**2*y3**2*y4*y5 + F110110*y1*(y2 + y3)*y4*y5 + F210110*y1**2*(y2 + y3)*y4*y5 + &
     F310110*y1**3*(y2 + y3)*y4*y5 + F120110*y1*(y2**2 + y3**2)*y4*y5 + F220110*y1**2*(y2**2 + y3**2)*y4*y5 + &
     F021110*(y2**2*y3 + y2*y3**2)*y4*y5 + F130110*y1*(y2**3 + y3**3)*y4*y5 + F031110*(y2**3*y3 + y2*y3**3)*y4*y5 + &
     F011220*y2*y3*y4**2*y5**2 + F110220*y1*(y2 + y3)*y4**2*y5**2 + F111100*y1*y2*y3*(y4 + y5) + &
     F211100*y1**2*y2*y3*(y4 + y5) + F311100*y1**3*y2*y3*(y4 + y5) + F122100*y1*y2**2*y3**2*(y4 + y5) + &
     F112100*y1*(y2*y3**2*y4 + y2**2*y3*y5) + F212100*y1**2*(y2*y3**2*y4 + y2**2*y3*y5) + &
     F113100*y1*(y2*y3**3*y4 + y2**3*y3*y5) + F121100*y1*(y2**2*y3*y4 + y2*y3**2*y5) + &
     F221100*y1**2*(y2**2*y3*y4 + y2*y3**2*y5) + F131100*y1*(y2**3*y3*y4 + y2*y3**3*y5) + &
     F111200*y1*y2*y3*(y4**2 + y5**2) + F211200*y1**2*y2*y3*(y4**2 + y5**2) + F112200*y1*(y2*y3**2*y4**2 + y2**2*y3*y5**2) + &
     F121200*y1*(y2**2*y3*y4**2 + y2*y3**2*y5**2) + F011210*y2*y3*(y4**2*y5 + y4*y5**2) + &
     F101210*y1*(y3*y4**2*y5 + y2*y4*y5**2) + F201210*y1**2*(y3*y4**2*y5 + y2*y4*y5**2) + &
     F102210*y1*(y3**2*y4**2*y5 + y2**2*y4*y5**2) + F110210*y1*(y2*y4**2*y5 + y3*y4*y5**2) + &
     F210210*y1**2*(y2*y4**2*y5 + y3*y4*y5**2) + F012210*(y2*y3**2*y4**2*y5 + y2**2*y3*y4*y5**2) + &
     F120210*y1*(y2**2*y4**2*y5 + y3**2*y4*y5**2) + F021210*(y2**2*y3*y4**2*y5 + y2*y3**2*y4*y5**2) + &
     F111300*y1*y2*y3*(y4**3 + y5**3) + F011310*y2*y3*(y4**3*y5 + y4*y5**3) + F101310*y1*(y3*y4**3*y5 + y2*y4*y5**3) + &
     F110310*y1*(y2*y4**3*y5 + y3*y4*y5**3) + F111001*y1*y2*y3*y6 + F211001*y1**2*y2*y3*y6 + F311001*y1**3*y2*y3*y6 + &
     F122001*y1*y2**2*y3**2*y6 + F121001*y1*(y2**2*y3 + y2*y3**2)*y6 + F221001*y1**2*(y2**2*y3 + y2*y3**2)*y6 + &
     F131001*y1*(y2**3*y3 + y2*y3**3)*y6 + F100111*y1*y4*y5*y6 + F200111*y1**2*y4*y5*y6 + F300111*y1**3*y4*y5*y6 + &
     F010111*(y2 + y3)*y4*y5*y6 + F020111*(y2**2 + y3**2)*y4*y5*y6 + F030111*(y2**3 + y3**3)*y4*y5*y6 + &
     F100221*y1*y4**2*y5**2*y6 + F010221*(y2 + y3)*y4**2*y5**2*y6 + F011101*y2*y3*(y4 + y5)*y6 + F022101*y2**2*y3**2*(y4+y5)*y6+&
     F101101*y1*(y3*y4 + y2*y5)*y6 + F201101*y1**2*(y3*y4 + y2*y5)*y6 + F301101*y1**3*(y3*y4 + y2*y5)*y6 + &
     F102101*y1*(y3**2*y4 + y2**2*y5)*y6 + F202101*y1**2*(y3**2*y4 + y2**2*y5)*y6 + F103101*y1*(y3**3*y4 + y2**3*y5)*y6 + &
     F110101*y1*(y2*y4 + y3*y5)*y6 + F210101*y1**2*(y2*y4 + y3*y5)*y6 + F310101*y1**3*(y2*y4 + y3*y5)*y6 + &
     F012101*(y2*y3**2*y4 + y2**2*y3*y5)*y6 + F013101*(y2*y3**3*y4 + y2**3*y3*y5)*y6 + F120101*y1*(y2**2*y4 + y3**2*y5)*y6 + &
     F220101*y1**2*(y2**2*y4 + y3**2*y5)*y6 + F021101*(y2**2*y3*y4 + y2*y3**2*y5)*y6 + F130101*y1*(y2**3*y4 + y3**3*y5)*y6 + &
     F031101*(y2**3*y3*y4 + y2*y3**3*y5)*y6 + F011201*y2*y3*(y4**2 + y5**2)*y6 + F101201*y1*(y3*y4**2 + y2*y5**2)*y6 + &
     F201201*y1**2*(y3*y4**2 + y2*y5**2)*y6 + F102201*y1*(y3**2*y4**2 + y2**2*y5**2)*y6 + F110201*y1*(y2*y4**2 + y3*y5**2)*y6 + &
     F210201*y1**2*(y2*y4**2 + y3*y5**2)*y6 + F012201*(y2*y3**2*y4**2 + y2**2*y3*y5**2)*y6 + &
     F120201*y1*(y2**2*y4**2 + y3**2*y5**2)*y6 + F021201*(y2**2*y3*y4**2 + y2*y3**2*y5**2)*y6 + &
     F100211*y1*(y4**2*y5 + y4*y5**2)*y6 + F200211*y1**2*(y4**2*y5 + y4*y5**2)*y6 + F001211*(y3*y4**2*y5 + y2*y4*y5**2)*y6 + &
     F002211*(y3**2*y4**2*y5 + y2**2*y4*y5**2)*y6 + F010211*(y2*y4**2*y5 + y3*y4*y5**2)*y6 + &
     F020211*(y2**2*y4**2*y5 + y3**2*y4*y5**2)*y6 + F011301*y2*y3*(y4**3 + y5**3)*y6 + F101301*y1*(y3*y4**3 + y2*y5**3)*y6 + &
     F110301*y1*(y2*y4**3 + y3*y5**3)*y6 + F100311*y1*(y4**3*y5 + y4*y5**3)*y6 + F001311*(y3*y4**3*y5 + y2*y4*y5**3)*y6 + &
     F010311*(y2*y4**3*y5 + y3*y4*y5**3)*y6 + F111002*y1*y2*y3*y6**2 + F211002*y1**2*y2*y3*y6**2 + &
     F121002*y1*(y2**2*y3 + y2*y3**2)*y6**2 + F100112*y1*y4*y5*y6**2 + F200112*y1**2*y4*y5*y6**2 + &
     F010112*(y2 + y3)*y4*y5*y6**2 + F020112*(y2**2 + y3**2)*y4*y5*y6**2 + F011102*y2*y3*(y4 + y5)*y6**2 + &
     F101102*y1*(y3*y4 + y2*y5)*y6**2 + F201102*y1**2*(y3*y4 + y2*y5)*y6**2 + F102102*y1*(y3**2*y4 + y2**2*y5)*y6**2 + &
     F110102*y1*(y2*y4 + y3*y5)*y6**2 + F210102*y1**2*(y2*y4 + y3*y5)*y6**2 + F012102*(y2*y3**2*y4 + y2**2*y3*y5)*y6**2 + &
     F120102*y1*(y2**2*y4 + y3**2*y5)*y6**2 + F021102*(y2**2*y3*y4 + y2*y3**2*y5)*y6**2 + F011202*y2*y3*(y4**2 + y5**2)*y6**2 + &
     F101202*y1*(y3*y4**2 + y2*y5**2)*y6**2 + F110202*y1*(y2*y4**2 + y3*y5**2)*y6**2 + F100212*y1*(y4**2*y5 + y4*y5**2)*y6**2 + &
     F001212*(y3*y4**2*y5 + y2*y4*y5**2)*y6**2 + F010212*(y2*y4**2*y5 + y3*y4*y5**2)*y6**2 + F111003*y1*y2*y3*y6**3 + &
     F100113*y1*y4*y5*y6**3 + F010113*(y2 + y3)*y4*y5*y6**3 + F011103*y2*y3*(y4 + y5)*y6**3 + F101103*y1*(y3*y4 + y2*y5)*y6**3 + &
     F110103*y1*(y2*y4 + y3*y5)*y6**3
     !
   endif
   if (nparams>=(501+ifst)) then
    !
     F111110 = param(482+ifst)
     F211110 = param(483+ifst)
     F121110 = param(484+ifst)
     F111210 = param(485+ifst)
     F111101 = param(486+ifst)
     F211101 = param(487+ifst)
     F121101 = param(488+ifst)
     F112101 = param(489+ifst)
     F110111 = param(490+ifst)
     F210111 = param(491+ifst)
     F011111 = param(492+ifst)
     F120111 = param(493+ifst)
     F021111 = param(494+ifst)
     F111201 = param(495+ifst)
     F110211 = param(496+ifst)
     F101211 = param(497+ifst)
     F011211 = param(498+ifst)
     F111102 = param(499+ifst)
     F110112 = param(500+ifst)
     F011112 = param(501+ifst)
     V5=F111110*y1*y2*y3*y4*y5 + F211110*y1**2*y2*y3*y4*y5 + F121110*y1*(y2**2*y3 + y2*y3**2)*y4*y5 + &
     F111210*y1*y2*y3*(y4**2*y5 + y4*y5**2) + F011111*y2*y3*y4*y5*y6 + F110111*y1*(y2 + y3)*y4*y5*y6 + &
     F210111*y1**2*(y2 + y3)*y4*y5*y6 + F120111*y1*(y2**2 + y3**2)*y4*y5*y6 + F021111*(y2**2*y3 + y2*y3**2)*y4*y5*y6 + &
     F111101*y1*y2*y3*(y4 + y5)*y6 + F211101*y1**2*y2*y3*(y4 + y5)*y6 + F112101*y1*(y2*y3**2*y4 + y2**2*y3*y5)*y6 + &
     F121101*y1*(y2**2*y3*y4 + y2*y3**2*y5)*y6 + F111201*y1*y2*y3*(y4**2 + y5**2)*y6 + F011211*y2*y3*(y4**2*y5 + y4*y5**2)*y6 + &
     F101211*y1*(y3*y4**2*y5 + y2*y4*y5**2)*y6 + F110211*y1*(y2*y4**2*y5 + y3*y4*y5**2)*y6 + F011112*y2*y3*y4*y5*y6**2 + &
     F110112*y1*(y2 + y3)*y4*y5*y6**2 + F111102*y1*y2*y3*(y4 + y5)*y6**2
     !
   endif
   if (nparams>=(502+ifst)) then
     F111111 = param(502+ifst)
     V6=F111111*y1*y2*y3*y4*y5*y6
   endif
   
   f=V0+V1+V2+V3+V4+V5+V6
 
 end function h2co_pes_6m_morse_delta_cos_mep
 
 
  
 function MLpoten_zxy2_mep_r_alpha_rho_powers(ncoords,natoms,x,xyz,force) result(f) 
  !
  integer(ik),intent(in) ::  ncoords,natoms
  real(ark),intent(in)   ::  x(ncoords)
  real(ark),intent(in)   ::  xyz(natoms,3)
  real(ark),intent(in)   ::  force(:)
  real(ark)              ::  f
  !
  integer(ik)          :: iterm,k_ind(6)
  real(ark)            :: y(6),xi(6),req(3),alphaeq(2),cosrho,rho,xieq(6),local(6)
  !
  local = from_local_to_r1r2r3a1a2tau(x,6)
  !
  rho = local(6)
  !
  cosrho = cos(rho) + 1.0_ark
  !
  ! reference
  !
  xieq(1:6) = ML_MEP_zxy2_R_rho(rho)
  !
  !xieq(1)     = sum(force(1:5 )*cosrho**molec%pot_ind(1,1:5 ))
  !xieq(2)     = sum(force(6:10)*cosrho**molec%pot_ind(2,6:10))
  !xieq(3)     = xieq(2)
  !xieq(4)     = sum(force(11:15)*cosrho**molec%pot_ind(4,11:15))
  !xieq(5)     = xieq(4)
  !
 ! expansion functions
 !
 y(1:3) = 1.0_ark-exp(-(local(1:3)-xieq(1:3)))
 y(4:5) = local(4:5)-xieq(4:5)
 y(6)   = cosrho
 !
 ! def potential energy
 !
 f = 0.0_ark
 !
 do iterm = 1,size(force)
   !
   xi(1:6) = y(1:6)**molec%pot_ind(1:6,iterm)
   !
   f = f + force(iterm)*product(xi(1:6))
   !
   if (molec%pot_ind(2,iterm)/=molec%pot_ind(3,iterm).or.molec%pot_ind(4,iterm)/=molec%pot_ind(5,iterm)) then 
     !
     k_ind(1) = molec%pot_ind(1,iterm)
     k_ind(2) = molec%pot_ind(3,iterm)
     k_ind(3) = molec%pot_ind(2,iterm)
     k_ind(4) = molec%pot_ind(5,iterm)
     k_ind(5) = molec%pot_ind(4,iterm)
     k_ind(6) = molec%pot_ind(6,iterm)
     !
     xi(1:6) = y(1:6)**k_ind(1:6)
     !
     f = f + force(iterm)*product(xi(1:6))
     !
   endif
   !
 enddo
 !
end function MLpoten_zxy2_mep_r_alpha_rho_powers


 function MLpoten_zxy2_mep_r_alpha_rho_powers_iso(ncoords,natoms,x,xyz,force) result(f) 
  !
  integer(ik),intent(in) ::  ncoords,natoms
  real(ark),intent(in)   ::  x(ncoords)
  real(ark),intent(in)   ::  xyz(natoms,3)
  real(ark),intent(in)   ::  force(:)
  real(ark)              ::  f
  !
  integer(ik)          :: N,N0,N_iso
  real(ark)            :: M0,M_main,V_iso,V,M_iso
  !
  N0 = force(1)
  V = MLpoten_zxy2_mep_r_alpha_rho_powers(ncoords,natoms,x,xyz,force(1+1:1+N0))
  !
  N = N0+1
  N_iso  = force(N+1)
  M_main = force(N+2)
  M_iso  = force(N+3)
  V_iso = MLpoten_zxy2_mep_r_alpha_rho_powers(ncoords,natoms,x,xyz,force(N+4:N+3+N_iso))
  !
  f = V + V_iso*(M_main-M_iso)/M_main
  !
end function MLpoten_zxy2_mep_r_alpha_rho_powers_iso


  !
  ! Defining potential energy function 
  !
  function MLpoten_zxy2_andrey_01(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f

    integer(ik)        :: N
    real(ark) thetae,reco,rech,betao,betah,f0,rh1,rh2,rco,alpha1,alpha2,tau
    real(ark) xi1, xi2, xi3, xi4, xi5, xi6
    real(ark) alphae
    real(ark) alphaedg
    real(ark) t1, t2, t3, t4, t5, t6, t0
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
    real(ark) s1, s2, s3, s4, s5
    real(ark) ta1, ta2, ta3
    !
    N = size(force)
    !
    thetae=force( 1)*pi/180.0_ark
    !
    alphae = pi-thetae*0.5_ark
    !
    reco      =  force(  2)
    rech      =  force(  3)
    betao     =  force(  4)
    betah     =  force(  5)
    F0      =  force(    6)
    F1      =  force(    7)
    F3      =  force(    8)
    F4      =  force(    9)
    F5      =  force(   10)

    if (N>=11) then
      F11     =  force(   11)
      F13     =  force(   12)
      F14     =  force(   13)
      F16     =  force(   14)
      F23     =  force(   15)
      F33     =  force(   16)
      F34     =  force(   17)
      F35     =  force(   18)
      F36     =  force(   19)
      F44     =  force(   20)
      F46     =  force(   21)
      F55     =  force(   22)
      F56     =  force(   23)
    endif
    if (N>=24) then
      F111    =  force(   24)
      F112    =  force(   25)
      F114    =  force(   26)
      F115    =  force(   27)
      F123    =  force(   28)
      F124    =  force(   29)
      F133    =  force(   30)
      F135    =  force(   31)
      F136    =  force(   32)
      F144    =  force(   33)
      F146    =  force(   34)
      F155    =  force(   35)
      F156    =  force(   36)
      F223    =  force(   37)
      F225    =  force(   38)
      F226    =  force(   39)
      F234    =  force(   40)
      F235    =  force(   41)
      F245    =  force(   42)
      F246    =  force(   43)
      F255    =  force(   44)
      F256    =  force(   45)
      F266    =  force(   46)
      F333    =  force(   47)
      F334    =  force(   48)
      F344    =  force(   49)
      F444    =  force(   50)
      F445    =  force(   51)
      F456    =  force(   52)
      F466    =  force(   53)
      F555    =  force(   54)
      F556    =  force(   55)
    endif
    if (N>=56) then
      F1111   =  force(   56)
      F1112   =  force(   57)
      F1114   =  force(   58)
      F1115   =  force(   59)
      F1122   =  force(   60)
      F1123   =  force(   61)
      F1124   =  force(   62)
      F1126   =  force(   63)
      F1136   =  force(   64)
      F1144   =  force(   65)
      F1146   =  force(   66)
      F1155   =  force(   67)
      F1156   =  force(   68)
      F1222   =  force(   69)
      F1225   =  force(   70)
      F1233   =  force(   71)
      F1234   =  force(   72)
      F1235   =  force(   73)
      F1244   =  force(   74)
      F1245   =  force(   75)
      F1246   =  force(   76)
      F1256   =  force(   77)
      F1334   =  force(   78)
      F1335   =  force(   79)
      F1355   =  force(   80)
      F1366   =  force(   81)
      F1444   =  force(   82)
      F1445   =  force(   83)
      F1456   =  force(   84)
      F1466   =  force(   85)
      F1566   =  force(   86)
      F1666   =  force(   87)
      F2222   =  force(   88)
      F2224   =  force(   89)
      F2225   =  force(   90)
      F2233   =  force(   91)
      F2244   =  force(   92)
      F2245   =  force(   93)
      F2246   =  force(   94)
      F2255   =  force(   95)
      F2256   =  force(   96)
      F2266   =  force(   97)
      F2333   =  force(   98)
      F2334   =  force(   99)
      F2335   =  force(  100)
      F2336   =  force(  101)
      F2344   =  force(  102)
      F2345   =  force(  103)
      F2356   =  force(  104)
      F2366   =  force(  105)
      F2444   =  force(  106)
      F2445   =  force(  107)
      F2455   =  force(  108)
      F2456   =  force(  109)
      F2466   =  force(  110)
      F3335   =  force(  111)
      F3445   =  force(  112)
      F3555   =  force(  113)
      F3556   =  force(  114)
      F3566   =  force(  115)
      F3666   =  force(  116)
      F4444   =  force(  117)
      F4445   =  force(  118)
      F4456   =  force(  119)
      F4466   =  force(  120)
      F4555   =  force(  121)
      F4556   =  force(  122)
      F5566   =  force(  123)
      F5666   =  force(  124)
      F6666   =  force(  125)
    endif
    if (N>=126) then
      F11111  =  force(  126)
      F11112  =  force(  127)
      F11114  =  force(  128)
      F11116  =  force(  129)
      F11123  =  force(  130)
      F11124  =  force(  131)
      F11133  =  force(  132)
      F11135  =  force(  133)
      F11136  =  force(  134)
      F11144  =  force(  135)
      F11145  =  force(  136)
      F11156  =  force(  137)
      F11166  =  force(  138)
      F11223  =  force(  139)
      F11224  =  force(  140)
      F11234  =  force(  141)
      F11236  =  force(  142)
      F11244  =  force(  143)
      F11245  =  force(  144)
      F11255  =  force(  145)
      F11266  =  force(  146)
      F11333  =  force(  147)
      F11335  =  force(  148)
      F11336  =  force(  149)
      F11345  =  force(  150)
      F11356  =  force(  151)
      F11444  =  force(  152)
      F11446  =  force(  153)
      F11455  =  force(  154)
      F11456  =  force(  155)
      F11555  =  force(  156)
      F11556  =  force(  157)
      F12225  =  force(  158)
      F12233  =  force(  159)
      F12236  =  force(  160)
      F12255  =  force(  161)
      F12266  =  force(  162)
      F12333  =  force(  163)
      F12334  =  force(  164)
      F12336  =  force(  165)
      F12344  =  force(  166)
      F12345  =  force(  167)
      F12355  =  force(  168)
      F12356  =  force(  169)
      F12444  =  force(  170)
      F12455  =  force(  171)
      F12466  =  force(  172)
      F12666  =  force(  173)
      F13333  =  force(  174)
      F13334  =  force(  175)
      F13335  =  force(  176)
      F13344  =  force(  177)
      F13345  =  force(  178)
      F13346  =  force(  179)
      F13356  =  force(  180)
      F13445  =  force(  181)
      F13446  =  force(  182)
      F13456  =  force(  183)
      F13556  =  force(  184)
      F13566  =  force(  185)
      F13666  =  force(  186)
      F14444  =  force(  187)
      F14445  =  force(  188)
      F14456  =  force(  189)
      F14466  =  force(  190)
      F14556  =  force(  191)
      F14666  =  force(  192)
      F15555  =  force(  193)
      F15556  =  force(  194)
      F15566  =  force(  195)
      F22224  =  force(  196)
      F22333  =  force(  197)
      F22334  =  force(  198)
      F22336  =  force(  199)
      F22346  =  force(  200)
      F22446  =  force(  201)
      F22456  =  force(  202)
      F23333  =  force(  203)
      F23334  =  force(  204)
      F23335  =  force(  205)
      F23336  =  force(  206)
      F23344  =  force(  207)
      F23346  =  force(  208)
      F23355  =  force(  209)
      F23356  =  force(  210)
      F23366  =  force(  211)
      F23444  =  force(  212)
      F23446  =  force(  213)
      F23456  =  force(  214)
      F23466  =  force(  215)
      F23555  =  force(  216)
      F23556  =  force(  217)
      F33333  =  force(  218)
      F33335  =  force(  219)
      F33336  =  force(  220)
      F33344  =  force(  221)
      F33345  =  force(  222)
      F33346  =  force(  223)
      F33355  =  force(  224)
      F33356  =  force(  225)
      F33366  =  force(  226)
      F33444  =  force(  227)
      F33446  =  force(  228)
      F33455  =  force(  229)
      F33466  =  force(  230)
      F33555  =  force(  231)
      F33556  =  force(  232)
      F33566  =  force(  233)
      F33666  =  force(  234)
      F34444  =  force(  235)
      F34445  =  force(  236)
      F34446  =  force(  237)
      F34455  =  force(  238)
      F34456  =  force(  239)
      F34466  =  force(  240)
      F34555  =  force(  241)
      F34556  =  force(  242)
      F34566  =  force(  243)
      F34666  =  force(  244)
      F35555  =  force(  245)
      F35556  =  force(  246)
      F35566  =  force(  247)
      F35666  =  force(  248)
      F36666  =  force(  249)
      F44444  =  force(  250)
      F44445  =  force(  251)
      F44456  =  force(  252)
      F44466  =  force(  253)
      F44566  =  force(  254)
      F44666  =  force(  255)
      F45566  =  force(  256)
      F45666  =  force(  257)
      F46666  =  force(  258)
      F55555  =  force(  259)
      F55666  =  force(  260)
      F56666  =  force(  261)
    endif
    if (N>=262) then
      F111111 =  force(  262)
      F111112 =  force(  263)
      F111114 =  force(  264)
      F111116 =  force(  265)
      F111123 =  force(  266)
      F111124 =  force(  267)
      F111126 =  force(  268)
      F111133 =  force(  269)
      F111136 =  force(  270)
      F111144 =  force(  271)
      F111145 =  force(  272)
      F111156 =  force(  273)
      F111166 =  force(  274)
      F111222 =  force(  275)
      F111223 =  force(  276)
      F111234 =  force(  277)
      F111236 =  force(  278)
      F111246 =  force(  279)
      F111255 =  force(  280)
      F111266 =  force(  281)
      F111334 =  force(  282)
      F111335 =  force(  283)
      F111336 =  force(  284)
      F111344 =  force(  285)
      F111346 =  force(  286)
      F111356 =  force(  287)
      F111444 =  force(  288)
      F111446 =  force(  289)
      F111456 =  force(  290)
      F111466 =  force(  291)
      F111556 =  force(  292)
      F111666 =  force(  293)
      F112233 =  force(  294)
      F112234 =  force(  295)
      F112235 =  force(  296)
      F112236 =  force(  297)
      F112244 =  force(  298)
      F112245 =  force(  299)
      F112246 =  force(  300)
      F112255 =  force(  301)
      F112256 =  force(  302)
      F112266 =  force(  303)
      F112333 =  force(  304)
      F112344 =  force(  305)
      F112345 =  force(  306)
      F112356 =  force(  307)
      F112366 =  force(  308)
      F112444 =  force(  309)
      F112445 =  force(  310)
      F112446 =  force(  311)
      F112455 =  force(  312)
      F112456 =  force(  313)
      F112466 =  force(  314)
      F112555 =  force(  315)
      F112556 =  force(  316)
      F112566 =  force(  317)
      F112666 =  force(  318)
      F113333 =  force(  319)
      F113334 =  force(  320)
      F113335 =  force(  321)
      F113336 =  force(  322)
      F114444 =  force(  323)
      F114446 =  force(  324)
      F114456 =  force(  325)
      F114466 =  force(  326)
      F114556 =  force(  327)
      F114666 =  force(  328)
      F115555 =  force(  329)
      F115566 =  force(  330)
      F115666 =  force(  331)
      F122222 =  force(  332)
      F122223 =  force(  333)
      F122224 =  force(  334)
      F122225 =  force(  335)
      F122226 =  force(  336)
      F122234 =  force(  337)
      F122235 =  force(  338)
      F122236 =  force(  339)
      F122244 =  force(  340)
      F122245 =  force(  341)
      F122246 =  force(  342)
      F122255 =  force(  343)
      F122256 =  force(  344)
      F122266 =  force(  345)
      F122333 =  force(  346)
      F122334 =  force(  347)
      F122336 =  force(  348)
      F122344 =  force(  349)
      F122345 =  force(  350)
      F122346 =  force(  351)
      F122444 =  force(  352)
      F122445 =  force(  353)
      F122446 =  force(  354)
      F122455 =  force(  355)
      F122456 =  force(  356)
      F122466 =  force(  357)
      F122555 =  force(  358)
      F122666 =  force(  359)
      F123355 =  force(  360)
      F123356 =  force(  361)
      F123366 =  force(  362)
      F123444 =  force(  363)
      F123446 =  force(  364)
      F123456 =  force(  365)
      F123466 =  force(  366)
      F123555 =  force(  367)
      F123566 =  force(  368)
      F124444 =  force(  369)
      F124445 =  force(  370)
      F124446 =  force(  371)
      F124455 =  force(  372)
      F125555 =  force(  373)
      F125556 =  force(  374)
      F125566 =  force(  375)
      F125666 =  force(  376)
      F126666 =  force(  377)
      F133556 =  force(  378)
      F133566 =  force(  379)
      F134455 =  force(  380)
      F134456 =  force(  381)
      F134555 =  force(  382)
      F134556 =  force(  383)
      F134566 =  force(  384)
      F134666 =  force(  385)
      F144444 =  force(  386)
      F144446 =  force(  387)
      F144456 =  force(  388)
      F144466 =  force(  389)
      F144566 =  force(  390)
      F144666 =  force(  391)
      F145566 =  force(  392)
      F145666 =  force(  393)
      F146666 =  force(  394)
      F155556 =  force(  395)
      F155666 =  force(  396)
      F166666 =  force(  397)
      F222222 =  force(  398)
      F222223 =  force(  399)
      F222224 =  force(  400)
      F222225 =  force(  401)
      F222226 =  force(  402)
      F222234 =  force(  403)
      F222236 =  force(  404)
      F222244 =  force(  405)
      F222245 =  force(  406)
      F222246 =  force(  407)
      F222255 =  force(  408)
      F222256 =  force(  409)
      F222266 =  force(  410)
      F222333 =  force(  411)
      F222334 =  force(  412)
      F222335 =  force(  413)
      F222336 =  force(  414)
      F222344 =  force(  415)
      F222345 =  force(  416)
      F222346 =  force(  417)
      F222355 =  force(  418)
      F222356 =  force(  419)
      F222366 =  force(  420)
      F222444 =  force(  421)
      F222445 =  force(  422)
      F222446 =  force(  423)
      F222455 =  force(  424)
      F222456 =  force(  425)
      F222466 =  force(  426)
      F222555 =  force(  427)
      F222556 =  force(  428)
      F222566 =  force(  429)
      F222666 =  force(  430)
      F223333 =  force(  431)
      F223344 =  force(  432)
      F223346 =  force(  433)
      F223356 =  force(  434)
      F223366 =  force(  435)
      F223444 =  force(  436)
      F223445 =  force(  437)
      F223446 =  force(  438)
      F223456 =  force(  439)
      F223556 =  force(  440)
      F223566 =  force(  441)
      F223666 =  force(  442)
      F224444 =  force(  443)
      F224445 =  force(  444)
      F224446 =  force(  445)
      F224455 =  force(  446)
      F224456 =  force(  447)
      F224466 =  force(  448)
      F224555 =  force(  449)
      F224556 =  force(  450)
      F224566 =  force(  451)
      F224666 =  force(  452)
      F225555 =  force(  453)
      F225556 =  force(  454)
      F225566 =  force(  455)
      F225666 =  force(  456)
      F226666 =  force(  457)
      F233336 =  force(  458)
      F233455 =  force(  459)
      F233466 =  force(  460)
      F233666 =  force(  461)
      F234444 =  force(  462)
      F234446 =  force(  463)
      F234456 =  force(  464)
      F234466 =  force(  465)
      F234566 =  force(  466)
      F234666 =  force(  467)
      F235555 =  force(  468)
      F235556 =  force(  469)
      F235566 =  force(  470)
      F244444 =  force(  471)
      F244445 =  force(  472)
      F244446 =  force(  473)
      F244455 =  force(  474)
      F244456 =  force(  475)
      F244466 =  force(  476)
      F244555 =  force(  477)
      F244556 =  force(  478)
      F244566 =  force(  479)
      F244666 =  force(  480)
      F245555 =  force(  481)
      F245556 =  force(  482)
      F245566 =  force(  483)
      F245666 =  force(  484)
      F246666 =  force(  485)
      F255555 =  force(  486)
      F255556 =  force(  487)
      F255566 =  force(  488)
      F255666 =  force(  489)
      F256666 =  force(  490)
      F266666 =  force(  491)
      F444444 =  force(  492)
      F444446 =  force(  493)
      F444456 =  force(  494)
      F444466 =  force(  495)
      F444566 =  force(  496)
      F444666 =  force(  497)
      F445555 =  force(  498)
      F445566 =  force(  499)
      F445666 =  force(  500)
      F455666 =  force(  501)
      F456666 =  force(  502)
      F466666 =  force(  503)
      F555556 =  force(  504)
      F555566 =  force(  505)
      F555666 =  force(  506)
      F666666 =  force(  507)
    endif

    rco    = local(1) ; rh1    = local(2)  ; rh2   = local(3)
    alpha1 = local(4) ; alpha2 = local(5)  ; tau   = local(6)

    xi1=1.0_ark-exp(-betao*(rco-reco))
    xi2=1.0_ark-exp(-betah*(rh1-rech))
    xi3=1.0_ark-exp(-betah*(rh2-rech))
    !
    xi4 = cos(tau)+1.0_ark
    !xi5 = cos(alpha1)-cos(alphae)
    !xi6 = cos(alpha2)-cos(alphae)
    !
    !xi4 = (tau-pi)**2
    xi5 = alpha1-alphae
    xi6 = alpha2-alphae

    t0 = F0
    !
    t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0
    if (N>=1) then 
     t1 = F1*xi1+(xi2+xi3)*F3+(xi6+xi5)*F5+F4*xi4
    endif 

    if (N>=11) then 
     t2 = F11*xi1**2+F14*xi1*xi4+F23*xi2*xi3+(xi1*xi3+xi1*xi2)*F13+(xi1*xi5+xi1*xi6)*F16+&
       (xi2*xi4+xi3*xi4)*F34+(xi3**2+xi2**2)*F33+(xi3*xi6+xi2*xi5)*F36+(xi2*xi6+xi3*xi5)*F35+&
       F44*xi4**2+(xi4*xi6+xi4*xi5)*F46+(xi6**2+xi5**2)*F55+F56*xi5*xi6
    endif                                                              
                                                                       
    if (N>=24) then                                                 
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
                                                                       
    if (N>=56) then                                                       
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


   if (N>=126) then
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

   if(N>=262) then                                                        
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


   f = (t0+t1+t2+t3+t4+t5+t6)

  end function MLpoten_zxy2_andrey_01

  !
  ! Defining potential energy function 
  !
  ! This type is for XY3-molecules, Taylor type II of expansion 
  ! with respect to the symmetrized GD-coordinates
  ! V =  sum_{i<=j} f_{i,j} S_i S_j + sum_{i<=j<=k} f_{i,j,k} S_i S_j S_k ...
  ! with S2(A1) = (a1+a2+a3) /sqrt(3)
  ! The fortran code has been generated by Maple routine "RecalcPoten.mws"
  ! For example, NH3 PES is given by this expression
  !
  function MLpoten_zxy2_mlt(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f

   real(ark)    :: r1,r2,r3,r4,r5,r6,rho1,rho2,rho3
   real(ark)    :: rCH1,rCH2,rCO,theta1,theta2,tau,theta3,costau

   real(ark)    :: s1,v2,v3,v4, &
      reCH,reCO,thetae,taue,   &
      w11  ,w12  ,w13  ,w22  ,w23  ,w33  ,w44  ,w55  ,w56  , &
      w66  ,w111 ,w112 ,w113 ,w122 ,w123 ,w133 ,w144 ,w155 , &
      w156 ,w166 ,w222 ,w223 ,w233 ,w244 ,w255 ,w256 ,w266 , &
      w333 ,w344 ,w355 ,w356 ,w366 ,w1111,w1112,w1113,w1122, &
      w1123,w1133,w1144,w1155,w1156,w1166,w1222,w1223,w1233, &
      w1244,w1255,w1256,w1266,w1333,w1344,w1355,w1358,w1366, &
      w2222,w2223,w2233,w2244,w2255,w2256,w2266,w2333,w2344, &
      w2355,w2356,w2366,w3333,w3344,w3355,w3350,w3366,w4444, &
      w4455,w4456,w4466,w5555,w5556,w5568,w5666,w6666


      if (verbose>=6) write(out,"('MLpoten_zxy2_mlt/start')") 

      reCO    =  molec%req(1)
      reCH    =  molec%req(2)
      thetae  =  molec%alphaeq(1)
      taue    =  molec%taueq(1)

      w11   =  force(  1)
      w12   =  force(  2)
      w13   =  force(  3)
      w22   =  force(  4)
      w23   =  force(  5)
      w33   =  force(  6)
      w44   =  force(  7)
      w55   =  force(  8)
      w56   =  force(  9)
      w66   =  force( 10)
      w111  =  force( 11)
      w112  =  force( 12)
      w113  =  force( 13)
      w122  =  force( 14)
      w123  =  force( 15)
      w133  =  force( 16)
      w144  =  force( 17)
      w155  =  force( 18)
      w156  =  force( 19)
      w166  =  force( 20)
      w222  =  force( 21)
      w223  =  force( 22)
      w233  =  force( 23)
      w244  =  force( 24)
      w255  =  force( 25)
      w256  =  force( 26)
      w266  =  force( 27)
      w333  =  force( 28)
      w344  =  force( 29)
      w355  =  force( 30)
      w356  =  force( 31)
      w366  =  force( 32)
      w1111 =  force( 33)
      w1112 =  force( 34)
      w1113 =  force( 35)
      w1122 =  force( 36)
      w1123 =  force( 37)
      w1133 =  force( 38)
      w1144 =  force( 39)
      w1155 =  force( 40)
      w1156 =  force( 41)
      w1166 =  force( 42)
      w1222 =  force( 43)
      w1223 =  force( 44)
      w1233 =  force( 45)
      w1244 =  force( 46)
      w1255 =  force( 47)
      w1256 =  force( 48)
      w1266 =  force( 49)
      w1333 =  force( 50)
      w1344 =  force( 51)
      w1355 =  force( 52)
      w1358 =  force( 53)
      w1366 =  force( 54)
      w2222 =  force( 55)
      w2223 =  force( 56)
      w2233 =  force( 57)
      w2244 =  force( 58)
      w2255 =  force( 59)
      w2256 =  force( 60)
      w2266 =  force( 61)
      w2333 =  force( 62)
      w2344 =  force( 63)
      w2355 =  force( 64)
      w2356 =  force( 65)
      w2366 =  force( 66)
      w3333 =  force( 67)
      w3344 =  force( 68)
      w3355 =  force( 69)
      w3350 =  force( 70)
      w3366 =  force( 71)
      w4444 =  force( 72)
      w4455 =  force( 73)
      w4456 =  force( 74)
      w4466 =  force( 75)
      w5555 =  force( 76)
      w5556 =  force( 77)
      w5568 =  force( 78)
      w5666 =  force( 79)
      w6666 =  force( 80)
      !
      theta1 = local(4)
      theta2 = local(5)
      !
      select case(molec%dihedtype(1))
      case default
         write (out,"('ML_b0_ZXY2: dihedral type ',i4,' is not working here')") molec%dihedtype(1)
         stop 'ML_b0_ZXY2 - bad dihedral type'
         !
      case( 1) 
         !
         theta3 = local(6)
         !
         costau = -(cos(theta1)*cos(theta2)-cos(theta3))/(sin(theta1)*sin(theta2))
         !
         if ( abs(costau)>1.0_ark+sqrt(small_) ) then 
            !
            write (out,"('MLpoten_zxy2_mlt: costau>1: ',f18.8)") costau
            stop 'MLpoten_zxy2_mlt - bad costau'
            !
         elseif ( costau>=1.0_ark) then 
            !
            tau = 0.0_ark
            !
         elseif ( costau<=-1.0_ark) then 
            !
            tau = pi
            !
         else
            ! 
            tau = acos(costau)
            !
         endif
         !
         taue = pi
         !
         !costau_2 = (1._ark-cos(theta1)**2-cos(theta2)**2+cos(theta2)**2*cos(theta1)**2-local(7)**2)/&
         !           (1._ark-cos(theta1)**2-cos(theta2)**2+cos(theta2)**2*cos(theta1)**2)

         !tau = acos(-sqrt(costau_2))
         !
         !continue 
         !
      case(-2) 
         !
         tau    = local(6)
         taue = molec%taueq(1)
         !
      end select 
      !
      rCO    = local(1)
      rCH1   = local(2)
      rCH2   = local(3)
     
      !
      rho1 = (rCH1 - reCH)/(rCH1)
      rho2 = (rCH2 - reCH)/(rCH2)
      rho3 = (rCO  - reCO)/(rCO )
      !
      r1 = (rho1+rho2)/sqrt(2.0_ark)
      r2 =  rho3
      r3 = (theta1+theta2-2.0_ark*thetae)/sqrt(2.0_ark)
      r4 =  tau-taue
      r5 = (rho1-rho2)/sqrt(2.0_ark)
      r6 = (theta1-theta2)/sqrt(2.0_ark)
      !
      !------------------------------------------------------
      !
      v2 = w11*r1**2+w12*r2*r1+w22*r2**2+w13*r3*r1+w23*r3*r2+w33*r3**2+w44*r4**2+ & 
           w55*r5**2+w56*r6*r5+w66*r6**2

      !
      v3 = w111*r1**3+w112*r2*r1**2+w122*r2**2*r1+w222*r2**3+w113*r3*r1**2+&
           w123*r3*r2*r1+w223*r3*r2**2+w133*r3**2*r1+w233*r3**2*r2+w333*r3**3+&
           w144*r4**2*r1+w244*r4**2*r2+w344*r4**2*r3+w155*r5**2*r1+w255*r5**2*r2+&
           w355*r5**2*r3+w156*r6*r5*r1+w256*r6*r5*r2+w356*r6*r5*r3+w166*r6**2*r1+&
           w266*r6**2*r2+w366*r6**2*r3
      !
      !
      s1 = w1111*r1**4+w2222*r2**4+w3333*r3**4+w4444*r4**4+w5555*r5**4+w5556*r6*r5**3+          &
      w1112*r2*r1**3+w1123*r3*r2*r1**2+w1223*r3*r2**2*r1+w1233*r3**2*r2*r1+w1244*r4**2*r2*r1+   &
      w1344*r4**2*r3*r1+w2344*r4**2*r3*r2+w1255*r5**2*r2*r1+w1355*r5**2*r3*r1+w2355*r5**2*r3*r2+&
      w1222*r2**3*r1+w1113*r3*r1**3+w2244*r4**2*r2**2+w2255*r5**2*r2**2+w3344*r4**2*r3**2+      &
      w6666*r6**4

      v4 = s1+w1156*r6*r5*r1**2+w2256*r6*r5*r2**2+w4456*r6*r5*r4**2+w1266*r6**2*r2*r1+          &
      w1366*r6**2*r3*r1+w2366*r6**2*r3*r2+w1256*r6*r5*r2*r1+w2356*r6*r5*r3*r2+                  &
      w2233*r3**2*r2**2+w1166*r6**2*r1**2+w2333*r3**3*r2+w2223*r3*r2**3+w1155*r5**2*r1**2+      &
      w1144*r4**2*r1**2+w3355*r5**2*r3**2+w3366*r6**2*r3**2+w1333*r3**3*r1+w1133*r3**2*r1**2+   &
      w2266*r6**2*r2**2+w4466*r6**2*r4**2+w1122*r2**2*r1**2+w5666*r6**3*r5+w4455*r5**2*r4**2 
      
      f =  (v2+v3+v4)*real(1.0e-11/planck/vellgt,ark) !  50340.359783704819122_rk

      if (verbose>=6) write(out,"('MLpoten_zxy2_mlt/end')") 
 
 end function MLpoten_zxy2_mlt



 !function MLdms2xyz_zxy2_symadap_powers(xyz) result(f)
 recursive subroutine MLdms2xyz_zxy2_symadap_powers(rank,ncoords,natoms,local,xyz,f)
    !
    implicit none
    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)      :: i,imu,iterm,lwork,info,nsv
    real(ark)        :: x(molec%natoms,3),r_CO,r_CH1,r_CH2,a_H1CO,a_H2CO,tau,re_CO(3),re_CH(3),ae_HCO(3),&
                        y(6,3),xi(6),dip(3),N1(3),N2(3),NA1(3),NB1(3),NB2(3),rho,costheta,costau,tmat_ark(3,3),&
                        dip0(3),xyz0(natoms,3)
    double precision :: tmat(3,3),dipd(3,1),work(64*3),sv(3),svtol
    !
    real(ark)        :: x0(molec%natoms,3)
    !
    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      xyz0 = MLloc2pqr_zxy2(local)
      !
      forall(i=1:4) x(i,1:3)=xyz0(i,1:3)-xyz0(1,1:3) ! 1:C  2:O  3:H1  4:H2
      !
    else
      !
      ! C atom is in the origin
      !
      forall(i=1:4) x(i,1:3)=xyz(i,1:3)-xyz(1,1:3) ! 1:C  2:O  3:H1  4:H2
      !
    endif
    !
    ! internal coordinates
    !
    r_CO   = sqrt( sum(x(2,1:3)**2) )
    r_CH1  = sqrt( sum(x(3,1:3)**2) )
    r_CH2  = sqrt( sum(x(4,1:3)**2) )
    a_H1CO = acos( sum(x(3,1:3)*x(2,1:3))/(r_CH1*r_CO) )
    a_H2CO = acos( sum(x(4,1:3)*x(2,1:3))/(r_CH2*r_CO) )
    !
    N1(1:3)= MLvector_product( x(3,1:3)/r_CH1, x(2,1:3)/r_CO ) ! vec perp to H1CO plane
    N1(1:3)= N1(1:3)/sqrt( sum(N1(1:3)**2) )
    !
    N2(1:3)= MLvector_product( x(2,1:3)/r_CO, x(4,1:3)/r_CH2 ) ! vec perp to H2CO plane
    N2(1:3)= N2(1:3)/sqrt( sum(N2(1:3)**2) )
    !
    ! symmetry-adapted unit vectors
    !
    NA1(1:3) = x(2,1:3)/r_CO
    NA1(1:3) = NA1(1:3)/sqrt( sum(NA1(1:3)**2) )
    !
    NB1(1:3) = N1(1:3)+N2(1:3)
    NB1(1:3) = NB1(1:3)/sqrt( sum(NB1(1:3)**2) )
    !
    NB2(1:3) = MLvector_product( NA1(1:3),NB1(1:3) )
    NB2(1:3) = NB2(1:3)/sqrt( sum(NB2(1:3)**2) )
    !
    ! dihedral book angle
    !
    costau = sum(N1(1:3)*N2(1:3))
    !
    if ( abs(costau)>1.0_ark+sqrt(small_) ) then 
       !
       write (out,"('ML_coordinate_transform_ZXY2: costau>1: ',f18.8)") costau
       stop 'ML_coordinate_transform_ZXY2 - bad costau'
       !
    elseif ( costau>=1.0_ark) then 
       tau = 0
    else 
       tau = acos(costau)
    endif
    !
    costheta = sum(N1(1:3)*NB2(1:3))
    if ( sign(1.0_ark,costheta )<0.0) then
     rho = pi - tau
    else
     rho = pi + tau
    endif
    !
    ! expansion functions
    !
    re_CO(1:3)  = extF%coef(1,1:3)
    re_CH(1:3)  = extF%coef(2,1:3)
    ae_HCO(1:3) = extF%coef(3,1:3)/real(rad,kind=ark)
    !
    y(1,1:3) = r_CO-re_CO(1:3)
    y(2,1:3) = r_CH1-re_CH(1:3)
    y(3,1:3) = r_CH2-re_CH(1:3)
    y(4,1:3) = a_H1CO-ae_HCO(1:3)
    y(5,1:3) = a_H2CO-ae_HCO(1:3)
    y(6,1:3) = cos(rho)+1.0_ark
    !
    ! symmetry-adapted dipole moment (in the order a1 b1 b2)
    !
    dip(1:3)=0.0_ark
    !
    do imu=1,3
     do iterm=4,extF%nterms(imu) !(iterm=1:3 are equilibrium constants)
       !
       xi(1:6)=y(1:6,imu)**extF%term(1:6,iterm,imu)
       !
       dip(imu)=dip(imu)+extF%coef(iterm,imu)*product(xi(1:6))
       !
     enddo
    enddo
    !
    !dip(3) = -dip(3)
    !
    dip(2)=dip(2)*sin(rho)
    !
    ! transform dipole moment from symmetrized to xyz
    !
    tmat_ark(1,1:3) = NA1(1:3)
    tmat_ark(2,1:3) = NB1(1:3)
    tmat_ark(3,1:3) = NB2(1:3)
    !
    dip0 = dip
    !
    call MLlinurark(3,tmat_ark,dip(:),f,info)
    !
    if (info/=0) then 
      !
      dipd(1:3,1)= dip0(1:3)
      tmat = tmat_ark
      !
      lwork=size(work)
      svtol=-1.0d-12
      call dgelss(3,3,1,tmat,3,dipd,3,sv,svtol,nsv,work,lwork,info)
      if (info/=0) then
       write(out,'(/a,1x,i3)') 'MLdms2xyz_zxy2_symadap_powers error: SV decomposition failed, info=',info
       stop 'MLdms2xyz_zxy2_symadap_powers error: SV decomposition failed'
      endif
      !
      f(1:3)=real(dipd(1:3,1),kind=ark)
      !
    endif
    !
  end subroutine MLdms2xyz_zxy2_symadap_powers



recursive subroutine MLdms2xyz_zxy2_symadap_powers_tmp(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)      :: i,imu,iterm,lwork,info,nsv
    real(ark)        :: x(molec%natoms,3),r_CO,r_CH1,r_CH2,a_H1CO,a_H2CO,tau,re_CO(3),re_CH(3),ae_HCO(3),&
                        y(6,3),xi(6),dip(3),N1(3),N2(3),NA1(3),NB1(3),NB2(3),rho,costheta,costau
    double precision :: tmat(3,3),dipd(3,1),work(64*3),sv(3),svtol
    !
    integer(ik) :: IOunit
    logical :: dms_debug = .false., ifopened
    character(cl) :: IOname, filename
    !
    ! C atom is in the origin
    !
    forall(i=1:4) x(i,1:3)=xyz(i,1:3)-xyz(1,1:3) ! 1:C  2:O  3:H1  4:H2
    !
    ! internal coordinates
    !
    r_CO   = sqrt( sum(x(2,1:3)**2) )
    r_CH1  = sqrt( sum(x(3,1:3)**2) )
    r_CH2  = sqrt( sum(x(4,1:3)**2) )
    a_H1CO = acos( sum(x(3,1:3)*x(2,1:3))/(r_CH1*r_CO) )
    a_H2CO = acos( sum(x(4,1:3)*x(2,1:3))/(r_CH2*r_CO) )
    !
    N1(1:3)= MLvector_product( x(3,1:3)/r_CH1, x(2,1:3)/r_CO ) ! vec perp to H1CO plane
    N1(1:3)= N1(1:3)/sqrt( sum(N1(1:3)**2) )
    !
    N2(1:3)= MLvector_product( x(2,1:3)/r_CO, x(4,1:3)/r_CH2 ) ! vec perp to H2CO plane
    N2(1:3)= N2(1:3)/sqrt( sum(N2(1:3)**2) )
    !
    ! symmetry-adapted unit vectors
    !
    NA1(1:3) = x(2,1:3)/r_CO
    NA1(1:3) = NA1(1:3)/sqrt( sum(NA1(1:3)**2) )
    !
    NB1(1:3) = N1(1:3)+N2(1:3)
    NB1(1:3) = NB1(1:3)/sqrt( sum(NB1(1:3)**2) )
    !
    NB2(1:3) = MLvector_product( NA1(1:3),NB1(1:3) )
    NB2(1:3) = NB2(1:3)/sqrt( sum(NB2(1:3)**2) )
    !
    ! dihedral book angle
    !
    costau = sum(N1(1:3)*N2(1:3))
    !
    if ( abs(costau)>1.0_ark+sqrt(small_) ) then 
       !
       write (out,"('ML_coordinate_transform_ZXY2: costau>1: ',f18.8)") costau
       stop 'ML_coordinate_transform_ZXY2 - bad costau'
       !
    elseif ( costau>=1.0_ark) then 
       tau = 0
    else 
       tau = acos(costau)
    endif
    !
    costheta = sum(N1(1:3)*NB2(1:3))
    if (sign(1.0_ark,costheta)<0.0_ark) then
     rho = pi - tau
    else
     rho = pi + tau
    endif
    !
    ! expansion functions
    !
    re_CO(1:3)  = extF%coef(1,1:3)
    re_CH(1:3)  = extF%coef(2,1:3)
    ae_HCO(1:3) = extF%coef(3,1:3)/real(rad,kind=ark)
    !
    y(1,1:3) = r_CO-re_CO(1:3)
    y(2,1:3) = r_CH1-re_CH(1:3)
    y(3,1:3) = r_CH2-re_CH(1:3)
    y(4,1:3) = a_H1CO-ae_HCO(1:3)
    y(5,1:3) = a_H2CO-ae_HCO(1:3)
    y(6,1:3) = cos(rho)+1.0_ark
    !
    ! symmetry-adapted dipole moment (in the order a1 b1 b2)
    !
    dip(1:3)=0.0_ark
    !
    do imu=1,3
     do iterm=4,extF%nterms(imu) !(iterm=1:3 are equilibrium constants)
       !
       xi(1:6)=y(1:6,imu)**extF%term(1:6,iterm,imu)
       !
       dip(imu)=dip(imu)+extF%coef(iterm,imu)*product(xi(1:6))
       !
     enddo
    enddo
    !
    !dip(3) = -dip(3)
    !
    dip(2)=dip(2)*sin(rho)
    !
    ! transform dipole moment from symmetrized to xyz
    !
    tmat(1,1:3) = real(NA1(1:3),kind=rk)
    tmat(2,1:3) = real(NB1(1:3),kind=rk)
    tmat(3,1:3) = real(NB2(1:3),kind=rk)
    !
    dipd(1:3,1)=real(dip(1:3),kind=rk)
    !
    lwork=size(work)
    svtol=-1.0d-12
    call dgelss(3,3,1,tmat,3,dipd,3,sv,svtol,nsv,work,lwork,info)
    if (info/=0) then
     write(out,'(/a,1x,i3)') 'MLdms2xyz_zxy2_symadap_powers error: SV decomposition failed, info=',info
     stop 'MLdms2xyz_zxy2_symadap_powers error: SV decomposition failed'
    endif
    !
    f(1:3)=real(dipd(1:3,1),kind=ark)
    !
  end subroutine MLdms2xyz_zxy2_symadap_powers_tmp

  !
  ! Defining potential energy function 
  !
  !
  function MLpoten_SOHF(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f

   real(ark)    :: r1,r2,r3,v2,v0,reF,reS,reH,ae_SOF,ae_HOF,ae_HOS
   real(ark)    :: w1,w2,w3,w4,w5,f6,y6,phi,alpha1,alpha2,tau,alpha_SOF,alpha_HOS,alpha_HOF,cosphi


    if (verbose>=6) write(out,"('MLpoten_SOHF/start')") 
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('MLpoten_SOHF: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'MLpoten_SOHF - bad coord. type'
       !
    case('R-S-PHI')
      !
      reF    =  molec%req(1)
      reS    =  molec%req(2)
      reH    =  molec%req(2)
      ae_SOF  = molec%alphaeq(1)
      ae_HOF  = molec%alphaeq(2)
      ae_HOS  = molec%alphaeq(3)
      !
      w1   =  force(  1)
      w2   =  force(  2)
      w3   =  force(  3)
      w4   =  force(  4)
      w5   =  force(  5)
      f6   =  force(  6)
      !
      r1    = local(1)-reF
      r2    = local(2)-reS
      r3    = local(3)-reH
      !
      alpha1 = local(4)-ae_SOF
      alpha2 = local(5)-ae_HOF

      alpha_SOF = local(4)
      alpha_HOF = local(5)
      alpha_HOS = local(6)
      tau = local(7)
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
      y6 = sin(molec%taueq(1))-sin(phi)
      !
      v2 = w1*r1**2+w2*r2**2+w3*r3**2+w4*alpha1**2+w5*alpha2**2
      !
      v0 = f6*y6**2
      !
      f =  v0+v2
      !
    end select 
    !
    if (verbose>=6) write(out,"('MLpoten_SOHF/end')") 
 
 end function MLpoten_SOHF
 !


 function MLloc2pqr_zxy2(r) result(f)

 !return cartesian coordinates of atoms in the user-defined frame for locals specified

    real(ark), intent(in) :: r(molec%ncoords)
    real(ark)             :: f(molec%natoms, 3)

    integer(ik)           :: icart
    real(ark)             :: a0(molec%natoms, 3), cm,alpha1,alpha2,tau_2

    a0 = 0    !
    alpha1 = r(4)
    alpha2 = r(5)
    tau_2 = r(6)*0.5_ark
    !
    select case(trim(molec%coords_transform))
       !
    case('R-THETA-TAU')
       ! 
       a0(2,3) = r(1)
       !
       a0(3,1) = r(2)*sin(alpha1)*cos(tau_2)
       a0(3,2) =-r(2)*sin(alpha1)*sin(tau_2)
       a0(3,3) = r(2)*cos(alpha1)
       !
       a0(4,1) = r(3)*sin(alpha2)*cos(tau_2)
       a0(4,2) = r(3)*sin(alpha2)*sin(tau_2)
       a0(4,3) = r(3)*cos(alpha2)
       !
    case default 
       write(out,"('MLloc2pqr_xy2: illegal coordinate type',a)") trim(molec%coords_transform)
       stop 'MLloc2pqr_xy2: illegal coordinate type'
    end select
    
    do icart = 1, 3
       cm = sum(a0(1:molec%natoms, icart) * molec%atommasses(1:molec%natoms)) / sum(molec%atommasses(1:molec%natoms))
       a0(1:molec%natoms, icart) = a0(1:molec%natoms, icart) - cm
    end do
    !
    f(1:molec%natoms, 1:3) = a0(1:molec%natoms, 1:3)
    !
 end function MLloc2pqr_zxy2



  !
end module pot_zxy2

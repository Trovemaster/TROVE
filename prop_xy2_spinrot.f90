! Contains external functions for spin-rotation tensor of XY2-type molecule.

module prop_xy2_spinrot
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xy2
  use timer

  implicit none

  public prop_xy2_spin_rotation_bisector, prop_xy2_spin_rotation_bisector_nonlin, &
         TEST_prop_xy2_spin_rotation_bisector_nonlin, prop_xy2_gtensor_bisector, &
         prop_xy2_gtens_nuclear_bisector, prop_xy2_gtens_electronic_bisector,prop_xy2_gtensor_bisector_rank3


contains


! Nuclear rotational g-tensor for XY2-type molecule:
!      g_nuc = 1/(2*c) * Ie/Im,
!      where Im is inertia tensor and Ie is Im with atomic masses
!      replaced by atomic charges.
!      The returned values of g_nuc are in units of Debye/hbar.
!
! The corresponding nuclear contribution to magnetic moment:
!      mu_nuc = g_nuc * J

subroutine prop_xy2_gtens_nuclear_bisector(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha, rho, m0, m1, e0, e1, g(3,3), n1(3), n2(3), &
               x(natoms,3), rho_over_sinrhohalf
  real(ark), parameter  :: rho_threshold = 0.01_rk

  if (rank/=5) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_gtens_nuclear_bisector: rank of the input tensor =', rank, ', expected 5'
    stop
  endif

  ! xyz are undefined for the local case

  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_gtens_nuclear_bisector: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_gtens_nuclear_bisector - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  ! internal coordinates

  xyz0 = x(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = x(iatom,:) - xyz0(:)
  enddo

  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))

  n1 = xyz_(2,:)/r1
  n2 = xyz_(3,:)/r2

  alpha = aacos(sum(n1*n2))

  rho = pi-alpha

  ! g_nuc = 1/2 * Ie/Im

  m0 = molec%atomMasses(1) ! mass of atom X in XY2 molecule
  m1 = molec%atomMasses(2) ! mass of atom Y

  e0 = extF%coef(2,1) ! charge of atom X (in AU)
  e1 = extF%coef(3,1) ! charge of atom Y (in AU)
! NOTE: extF%coef(1,:) defines the power of rho-singularity for each tensor element
!       for this function extF%coef(1,1:5) must be equal to (/0,-1,0,0,0/)

  if (rho>rho_threshold) then
    rho_over_sinrhohalf = rho/sin(rho*0.5_ark)
  else
    rho_over_sinrhohalf = 2.0_ark + rho**2/12.0_ark + (7.0_ark*rho**4)/2880.0_ark &
                        + (31.0_ark*rho**6)/483840.0_ark + (127.0_ark*rho**8)/77414400.0_ark &
                        + (73.0_ark*rho**10)/1751777280.0_ark
  endif

  g = 0
  g(1,1) = (e1*(e1*m0*(r1 + r2)**2 + e0*(-(m1*(r1 - r2)**2) + 2.0_ark*m0*r1*r2))) &
           / (4.0_ark*(e0 + 2.0_ark*e1)*m0*m1*r1*r2)
  g(1,3) = (e1*(e1*m0 - e0*m1)*(r1 - r2)*(r1 + r2)*cos(rho*0.5_ark)*rho_over_sinrhohalf) &
           / (4.0_ark*(e0 + 2.0_ark*e1)*m0*m1*r1*r2)
  g(2,2) = (e1*(m0 + 2.0_ark*m1)*((e0 + e1)*(r1**2 + r2**2)+ 2.0_ark*e1*r1*r2*cos(rho))) &
           / (2.0_ark*(e0 + 2.0_ark*e1)*m1*((m0 + m1)*(r1**2 + r2**2) + 2.0_ark*m1*r1*r2*cos(rho)))
  g(3,1) = -(e1*(e1*m0 - e0*m1)*(r1 - r2)*(r1 + r2)*tan(rho*0.5_ark)) &
           / (4.0_ark*(e0 + 2.0_ark*e1)*m0*m1*r1*r2)
  g(3,3) = (e1*(-(e1*m0*(r1 - r2)**2) + e0*(2.0_ark*m0*r1*r2 + m1*(r1 + r2)**2))) &
           / (4.0_ark*(e0 + 2.0_ark*e1)*m0*m1*r1*r2)

  g = g * 1.0175071201112751e-05 ! to have units of g_nuc in Debye/hbar

  f = (/g(1,1), g(1,3), g(2,2), g(3,1), g(3,3)/) ! the second component g(1,3) must be divided by rho

end subroutine prop_xy2_gtens_nuclear_bisector



! Electronic rotational g-tensor for XY2-type molecule:
!      g_el = 1/hbar * mu_N(CGS) * 0.5 * (I*g + (I*g)^T) * G_rot,
!      where mu_N is nuclear magneton in CGS units, g is rotational g-tensor
!      from the electronic structure calculations, I is inertia tensor,
!      and G_rot is rotational kinetic energy matrix (same units as I^{-1})
!      The returned values of g_el are in units of Debye/hbar.
!
! The corresponding electronic contribution to magnetic moment:
!      mu_el = g_el * J

subroutine prop_xy2_gtens_electronic_bisector(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha, rho, m0, m1, e0, e1, g(3,3,-1:0), &
               n1(3), n2(3), x(natoms,3), Grot(3,3), gxx, gxz, gyy, gzx, gzz, mx, my, mu, &
               RhoOverSinRho, Rho2OverSin2RhoHalf, muxx, muxz, muzz, muyy, p1, p2
  real(ark), parameter :: muN = 5.050783699e-6 ! nuclear magneton in units of Debye
  real(ark), parameter  :: rho_threshold = 0.01_rk

  if (rank/=5) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_gtens_electronic_bisector: rank of the input tensor =', rank, ', expected 5'
    stop
  endif

  ! xyz are undefined for the local case

  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_gtens_electronic_bisector: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_gtens_electronic_bisector - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  ! internal coordinates

  xyz0 = x(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = x(iatom,:) - xyz0(:)
  enddo

  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))

  n1 = xyz_(2,:)/r1
  n2 = xyz_(3,:)/r2

  alpha = aacos(sum(n1*n2))

  rho = pi - alpha

  ! compute fitted elements of I*g matrix

  ! The first element in extF%coef(1,:) defines the power of rho-singularity (0 or -1)
  ! of the corresponding element in the resulting (output) g-tensor and not fitted I*g tensor
  ! (which is nonsingular), the expansion coefficients for I*g tensor start from the index no. 2.

  ! "_rhopow_min_one" in the name of the function means that the power of the rho coordinate
  ! in the expansion was reduced by one (since for these functions all expansion coefficients
  ! at rho**0 are equal to zero), meaning that the returned values must be multiplied by rho

  ! For bisector frame I*g elements (1,3) and (3,1) are equal to each other

  gxx =                fit_xy2_sr_A1(extF%nterms(1)-1, extF%coef(2:extF%nterms(1),1), (/r1, r2, alpha/))
  gxz = fit_xy2_sr_B2_rhopow_min_one(extF%nterms(2)-1, extF%coef(2:extF%nterms(2),2), (/r1, r2, alpha/))
  gyy =                fit_xy2_sr_A1(extF%nterms(3)-1, extF%coef(2:extF%nterms(3),3), (/r1, r2, alpha/))
  gzz = fit_xy2_sr_A1_rhopow_min_one(extF%nterms(5)-1, extF%coef(2:extF%nterms(5),5), (/r1, r2, alpha/))

  ! compute electronic g-tensor = 0.5 * (I*g + (I*g)^T) * G_rot

  if (rho>rho_threshold) then
    RhoOverSinRho = rho/sin(rho)
    Rho2OverSin2RhoHalf = rho**2/sin(rho*0.5_ark)**2
  else
    RhoOverSinRho = 1.0_ark + rho**2/6.0_ark + (7.0_ark*rho**4)/360.0_ark + (31.0_ark*rho**6)/15120.0_ark &
                  + (127.0_ark*rho**8)/604800.0_ark + (73.0_ark*rho**10)/3.42144e6_ark
    Rho2OverSin2RhoHalf = 4.0_ark + rho**2/3.0_ark + rho**4/60.0_ark + rho**6/1512.0_ark &
                        + rho**8/43200.0_ark + rho**10/1.33056e6_ark
  endif

  mx = molec%atomMasses(1)
  my = molec%atomMasses(2)
  !
  ! mu -> 1/mu
  mu = 1.0_ark/(1.0_ark/mx + 1.0_ark/my)

  g = 0
  !
  g(1,1,0)  = (gxz*mx*(0.5_ark*r1**2 - 0.5_ark*r2**2)*RhoOverSinRho + &
               gxx*(0.25_ark*mx*r1**2 - 0.5_ark*mu*r1*r2 + 0.25_ark*mx*r2**2)/cos(rho*0.5_ark)**2)/(mu*mx*r1**2*r2**2)
  g(1,3,-1) = (gxz*(0.25_ark*mx*r1**2 + 0.5_ark*mu*r1*r2 + 0.25_ark*mx*r2**2)*Rho2OverSin2RhoHalf + &
               gxx*mx*(0.5_ark*r1**2 - 0.5_ark*r2**2)*RhoOverSinRho)/(mu*mx*r1**2*r2**2)
  g(2,2,0)  = (gyy*(0.25_ark*mx*r1**2 + 0.25_ark*mx*r2**2 - 0.5_ark*mu*r1*r2*cos(rho)))/(mu*mx*r1**2*r2**2)
  g(3,1,0)  = (gzz*mx*(0.5_ark*r1**2 - 0.5_ark*r2**2)*RhoOverSinRho + gxz*(0.25_ark*mx*r1**2 - 0.5_ark*mu*r1*r2 + 0.25_ark*mx*r2**2)/cos(rho*0.5_ark)**2)/(mu*mx*r1**2*r2**2)
  ! 
  g(3,3,-1) = (gzz*(0.25_ark*mx*r1**2 + 0.5_ark*mu*r1*r2 + 0.25_ark*mx*r2**2)*Rho2OverSin2RhoHalf + gxz*mx*(0.5_ark*r1**2 - 0.5_ark*r2**2)*rho*RhoOverSinRho)/(mu*mx*r1**2*r2**2)

  !
  !p1 = 1.0_ark/r1
  !p2 = 1.0_ark/r2
  !
  !muxx = ( 0.25_ark*(mX+mY)/(mX*mY)*(p1**2+p2**2)- 0.5_ark/mX/(r1*r2) )/cos(rho*0.5_ark)**2
  !muyy = 0.25_ark*(mX+mY)/(mX*mY)*(p1**2+p2**2)- 0.5_ark/mX*cos(rho*0.5_ark)*p1*p2
  !muxz  = 0.5_ark*(p2**2 - p1**2)*(mX+mY)/mX/mY
  !muzz = 0.25_ark*(mX+mY)/mX/mY*(p1**2+p2**2)+.5_ark/mX/(r1*r2)
  
  !g(1,1,0)  = gxx*muxx+gxz*muxz*RhoOverSinRho
  !g(1,3,-1) = gxx*muxz*RhoOverSinRho + gxz*muzz*Rho2OverSin2RhoHalf
  !g(2,2,0)  = gyy*muyy
  !g(3,1,0)  = gxz*muxx+gzz*muxz*RhoOverSinRho
  !g(3,3,-1) = gzz*muzz*Rho2OverSin2RhoHalf+gxz*muxz*RhoOverSinRho*rho

  f = (/g(1,1,0), g(1,3,-1), g(2,2,0), g(3,1,0), g(3,3,-1)/)* muN

end subroutine prop_xy2_gtens_electronic_bisector


! Rotational g-tensor for quasi-linear molecule, like H2O, in the bisector
! frame, where some of the elements become singular at linear geometry.

subroutine prop_xy2_gtensor_bisector_tmp(rank, ncoords, natoms, local, xyz, f)
  !
  implicit none 
  !
  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)
  integer(ik) :: iatom
  !
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha, rho, a, b, c, d
  real(ark) :: g(3,3), mat(3,3), g_out(5,-2:0), e1(3), e2(3), e3(3), x(natoms,3), mY, mX
  real(ark) :: rho_over_sinrho, rho2_over_sinrho, rho2_over_sin2rhohalf
  real(ark),parameter  :: rho_threshold = 0.01_rk
  integer(ik) :: icentre
  !
  if (rank/=5) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_gtensor_bisector: rank of the input tensor =', rank, ', expected 5'
    stop
  endif
  !
  ! xyz are undefined for the local case
  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_gtensor_bisector: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_gtensor_bisector - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  xyz0 = x(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = x(iatom,:) - xyz0(:)
  enddo
  !
  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))
  !
  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2
  !
  alpha = aacos(sum(e1*e2))
  !
  rho = pi-alpha
  !
  ! fitted tensor elements
  !
  g = 0
  !
  g(1,1) = fit_xy2_sr_A1(extF%nterms(1)-2, extF%coef(3:extF%nterms(1),1), (/r1, r2, alpha/))
  g(2,2) = fit_xy2_sr_A1(extF%nterms(2)-2, extF%coef(3:extF%nterms(2),2), (/r1, r2, alpha/))
  g(3,3) = fit_xy2_sr_A1_rhopow_min_one(extF%nterms(3)-2, extF%coef(3:extF%nterms(3),3), (/r1, r2, alpha/))
  g(1,3) = fit_xy2_sr_B2_rhopow_min_one(extF%nterms(4)-2, extF%coef(3:extF%nterms(4),4), (/r1, r2, alpha/))
  g(3,1) = fit_xy2_sr_B2_rhopow_min_one(extF%nterms(5)-2, extF%coef(3:extF%nterms(5),5), (/r1, r2, alpha/))
  !
  ! transform with the inverse inertia tensor
  !
  if (rho>rho_threshold) then
    !
    rho_over_sinrho = rho/sin(rho)
    rho2_over_sinrho = rho**2/sin(rho)
    rho2_over_sin2rhohalf = rho**2/sin(rho*0.5_ark)**2
    !
  else
    !
    rho_over_sinrho = 1.0_ark + rho**2/6.0_ark + (7.0_ark*rho**4)/360.0_ark + (31.0_ark*rho**6)/15120.0_ark &
                    + (127.0_ark*rho**8)/604800.0_ark + (73.0_ark*rho**10)/3.42144e6_ark
    rho2_over_sinrho = rho + rho**3/6.0_ark + (7.0_ark*rho**5)/360.0_ark + (31.0_ark*rho**7)/15120.0_ark &
                     + (127._ark*rho**9)/604800.0_ark
    rho2_over_sin2rhohalf = 4.0_ark + rho**2/3.0_ark + rho**4/60.0_ark + rho**6/1512.0_ark &
                          + rho**8/43200.0_ark + rho**10/1.33056e6_ark
    !
  endif
  !
  ! inverse tensor of inertia
  !
  mX = molec%atomMasses(1)
  mY = molec%atomMasses(2)
  ! I^-1 = [[a,0,b/sin(rho)],[0,c,0],[b/sin(rho),0,d/sin(rho/2)^2]]
  a = ((mY*(r1 - r2)**2 + mX*(r1**2 + r2**2))*(1.0_ark/cos(rho/2.0_ark))**2)/(4.0_ark*mY*mX*r1**2*r2**2)
  b = ((mY + mX)*(r1**2 - r2**2))/(2.0_ark*mY*mX*r1**2*r2**2)
  c = (2.0_ark*mY + mX)/(mY*((mY + mX)*(r1**2 + r2**2) + 2.0_ark*mY*r1*r2*cos(rho)))
  d = ((mY*(r1 + r2)**2 + mX*(r1**2 + r2**2)))/(4.0_ark*mY*mX*r1**2*r2**2)
  !
  ! 1,1
  g_out(1,0) = g(1,1)*a + g(1,3)*b*rho_over_sinrho
  g_out(1,-1) = 0
  g_out(1,-2) = 0
  !
  ! 1,3
  g_out(2,0) = 0
  g_out(2,-1) = g(1,3)*d*rho2_over_sin2rhohalf + g(1,1)*b*rho_over_sinrho
  g_out(2,-2) = 0
  !
  ! 2,2
  g_out(3,0) = g(2,2)*c
  g_out(3,-1) = 0
  g_out(3,-2) = 0
  !
  ! 3,1
  g_out(4,0) = g(3,1)*a*rho + g(3,3)*b*rho_over_sinrho
  g_out(4,-1) = 0
  g_out(4,-2) = 0
  !
  ! 3,3
  g_out(5,0) = 0
  g_out(5,-1) = g(3,3)*d*rho2_over_sin2rhohalf + g(3,1)*b*rho2_over_sinrho
  g_out(5,-2) = 0
  !
  f = (/g_out(1,0), -g_out(2,-1), g_out(3,0), -g_out(4,0), g_out(5,-1)/)
  ! we multiply g(1,3) and g(3,1) here with -1 because of the opposite direction
  ! of the bisector axis (x-axis) in TROVE and in the actual ab initio data 
  !
end subroutine prop_xy2_gtensor_bisector_tmp





! Rotational g-tensor for quasi-linear molecule, like H2O, in the bisector
! frame, where some of the elements become singular at linear geometry.

subroutine prop_xy2_gtensor_bisector_tmp2(rank, ncoords, natoms, local, xyz, f)
  !
  implicit none 
  !
  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)
  integer(ik) :: iatom
  !
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha, rho, a, b, c, d
  real(ark) :: g(3,3), mat(3,3), g_out(5,-2:0), e1(3), e2(3), e3(3), x(natoms,3), mY, mX
  real(ark) :: rho_over_sinrho, rho2_over_sinrho, rho2_over_sin2rhohalf
  real(ark),parameter  :: rho_threshold = 0.01_rk
  integer(ik) :: icentre
  !
  if (rank/=5) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_gtensor_bisector: rank of the input tensor =', rank, ', expected 5'
    stop
  endif
  !
  ! xyz are undefined for the local case
  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_gtensor_bisector: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_gtensor_bisector - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  xyz0 = x(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = x(iatom,:) - xyz0(:)
  enddo
  !
  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))
  !
  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2
  !
  alpha = aacos(sum(e1*e2))
  !
  rho = pi-alpha
  !
  ! fitted tensor elements
  !
  g = 0
  !
  g(1,1) = fit_xy2_sr_A1(extF%nterms(1)-2, extF%coef(3:extF%nterms(1),1), (/r1, r2, alpha/))
  g(2,2) = fit_xy2_sr_A1(extF%nterms(2)-2, extF%coef(3:extF%nterms(2),2), (/r1, r2, alpha/))
  g(3,3) = fit_xy2_sr_A1_rhopow_min_one(extF%nterms(3)-2, extF%coef(3:extF%nterms(3),3), (/r1, r2, alpha/))
  !
  g(1,3) = fit_xy2_sr_B2(extF%nterms(4)-2, extF%coef(3:extF%nterms(4),4), (/r1, r2, alpha/))
  !
  g(3,1) = fit_xy2_sr_B2_rhopow_min_one(extF%nterms(5)-2, extF%coef(3:extF%nterms(5),5), (/r1, r2, alpha/))
  !
  ! transform with the inverse inertia tensor
  !
  if (rho>rho_threshold) then
    !
    rho_over_sinrho = rho/sin(rho)
    rho2_over_sinrho = rho**2/sin(rho)
    rho2_over_sin2rhohalf = rho**2/sin(rho*0.5_ark)**2
    !
  else
    !
    rho_over_sinrho = 1.0_ark + rho**2/6.0_ark + (7.0_ark*rho**4)/360.0_ark + (31.0_ark*rho**6)/15120.0_ark &
                    + (127.0_ark*rho**8)/604800.0_ark + (73.0_ark*rho**10)/3.42144e6_ark
    rho2_over_sinrho = rho + rho**3/6.0_ark + (7.0_ark*rho**5)/360.0_ark + (31.0_ark*rho**7)/15120.0_ark &
                     + (127._ark*rho**9)/604800.0_ark
    rho2_over_sin2rhohalf = 4.0_ark + rho**2/3.0_ark + rho**4/60.0_ark + rho**6/1512.0_ark &
                          + rho**8/43200.0_ark + rho**10/1.33056e6_ark
    !
  endif
  !
  ! inverse tensor of inertia
  !
  mX = molec%atomMasses(1)
  mY = molec%atomMasses(2)
  ! I^-1 = [[a,0,b/sin(rho)],[0,c,0],[b/sin(rho),0,d/sin(rho/2)^2]]
  a = ((mY*(r1 - r2)**2 + mX*(r1**2 + r2**2))*(1.0_ark/cos(rho/2.0_ark))**2)/(4.0_ark*mY*mX*r1**2*r2**2)
  b = ((mY + mX)*(r1**2 - r2**2))/(2.0_ark*mY*mX*r1**2*r2**2)
  c = (2.0_ark*mY + mX)/(mY*((mY + mX)*(r1**2 + r2**2) + 2.0_ark*mY*r1*r2*cos(rho)))
  d = ((mY*(r1 + r2)**2 + mX*(r1**2 + r2**2)))/(4.0_ark*mY*mX*r1**2*r2**2)
  !
  ! 1,1
  g_out(1,0) =  0 !g(1,1)*a + g(1,3)*b*rho_over_sinrho
  g_out(1,-1) = 0
  g_out(1,-2) = 0
  !
  ! 1,3
  g_out(2,0) = 0
  g_out(2,-1) = g(1,3) ! g(1,3)*d*rho2_over_sin2rhohalf + g(1,1)*b*rho_over_sinrho
  g_out(2,-2) = 0
  !
  ! 2,2
  g_out(3,0) =  0 ! g(2,2)*c
  g_out(3,-1) = 0
  g_out(3,-2) = 0
  !
  ! 3,1
  g_out(4,0) =  0 ! g(3,1)*a*rho + g(3,3)*b*rho_over_sinrho
  g_out(4,-1) = 0
  g_out(4,-2) = 0
  !
  ! 3,3
  g_out(5,0) = 0
  g_out(5,-1) = 0 ! g(3,3)*d*rho2_over_sin2rhohalf + g(3,1)*b*rho2_over_sinrho
  g_out(5,-2) = 0
  !
  f = (/g_out(1,0), g_out(2,-1), g_out(3,0), -g_out(4,0), g_out(5,-1)/)
  ! we multiply g(1,3) and g(3,1) here with -1 because of the opposite direction
  ! of the bisector axis (x-axis) in TROVE and in the actual ab initio data 
  !
end subroutine prop_xy2_gtensor_bisector_tmp2

! Rotational g-tensor for quasi-linear molecule, like H2O, in the bisector
! frame, where some of the elements become singular at linear geometry.

subroutine prop_xy2_gtensor_bisector(rank, ncoords, natoms, local, xyz, f)
  !
  implicit none 
  !
  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)
  integer(ik) :: iatom
  !
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha, rho, a, b, c, d
  real(ark) :: g(3,3), mat(3,3), g_out(5,-2:0), e1(3), e2(3), e3(3), x(natoms,3), mY, mX, mP
  real(ark) :: rho_over_sinrho, rho2_over_sinrho, rho2_over_sin2rhohalf,rho_over_sinrho_2,gtxz
  real(ark),parameter  :: rho_threshold = 0.01_rk
  integer(ik) :: icentre

  integer(ik)           :: k
  real(ark)             :: y1,y2,y3,re,ae,b0,muxz,muzz,Izz,g_zz,gxx,gyy,g1y,g2y,g3y
  real(ark)             :: p(3:extF%nterms(4)-2),v0,v1,v2,v3,v4,v5,v6,v7,v8,q(5:extF%nterms(3))

  !
  if (rank/=5) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_gtensor_bisector: rank of the input tensor =', rank, ', expected 5'
    stop
  endif
  !
  ! xyz are undefined for the local case
  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_gtensor_bisector: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_gtensor_bisector - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  xyz0 = x(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = x(iatom,:) - xyz0(:)
  enddo
  !
  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))
  !
  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2
  !
  alpha = aacos(sum(e1*e2))
  !
  rho = pi-alpha
  !
  ! fitted tensor elements
  !
  g = 0
  !
  g(1,1) = fit_xy2_sr_A1(extF%nterms(1)-2, extF%coef(3:extF%nterms(1),1), (/r1, r2, alpha/))
  g(2,2) = fit_xy2_sr_A1(extF%nterms(2)-2, extF%coef(3:extF%nterms(2),2), (/r1, r2, alpha/))
  g(3,3) = fit_xy2_sr_A1_rhopow_min_one(extF%nterms(3)-2, extF%coef(3:extF%nterms(3),3), (/r1, r2, alpha/))
  !g(1,3) = fit_xy2_sr_B2_rhopow_min_one(extF%nterms(4)-2, extF%coef(3:extF%nterms(4),4), (/r1, r2, alpha/))
  g(3,1) = fit_xy2_sr_B2_rhopow_min_one(extF%nterms(5)-2, extF%coef(3:extF%nterms(5),5), (/r1, r2, alpha/))
  !
  ! transform with the inverse inertia tensor
  !
  if (rho>rho_threshold) then
    !
    rho_over_sinrho = rho/sin(rho)
    rho_over_sinrho_2 = rho/sin(rho*0.5_ark)
    rho2_over_sinrho = rho**2/sin(rho)
    rho2_over_sin2rhohalf = rho**2/sin(rho*0.5_ark)**2
    !
  else
    !
    rho_over_sinrho = 1.0_ark + rho**2/6.0_ark + (7.0_ark*rho**4)/360.0_ark + (31.0_ark*rho**6)/15120.0_ark &
                    + (127.0_ark*rho**8)/604800.0_ark + (73.0_ark*rho**10)/3.42144e6_ark
    rho_over_sinrho_2 = 0.5_ark*rho-1._ark/48._ark*rho**3+1._ark/3840._ark*rho**5-1._ark/645120._ark*rho**7
    rho2_over_sinrho = rho + rho**3/6.0_ark + (7.0_ark*rho**5)/360.0_ark + (31.0_ark*rho**7)/15120.0_ark &
                     + (127._ark*rho**9)/604800.0_ark
    rho2_over_sin2rhohalf = 4.0_ark + rho**2/3.0_ark + rho**4/60.0_ark + rho**6/1512.0_ark &
                          + rho**8/43200.0_ark + rho**10/1.33056e6_ark
    !
  endif
  !
  mP = 1.007276_ark
  !
  mX = molec%atomMasses(1)
  mY = molec%atomMasses(2)
  !
  ! nuclear part of g[x,z] from Cybilski 1994 
  gtxz = -mP*4.0_ark*cos(rho/2.0_ark)*rho_over_sinrho_2*(&
            -11.0_ark*mY**2*r1**2+11.0_ark*mY**2*r2**2-4.0_ark*r1**2*mX**2&
            -8.0_ark*r1**2*mX*mY+4*r2**2*mX**2+8.0_ark*r2**2*mX*mY)/&
          ( (mX+2.0_ark*mY)*mY*(mX*r1**2+mX*r2**2+mY*r1**2-2.0_ark*mY*r1*r2+mY*r2**2) )

  !
  Izz = (mY*(mX*r2**2+mX*r1**2-2.0_ark*mY*r1*r2+mY*r1**2+mY*r2**2))/(mX+2.0_ark*mY)/cos(rho*0.5_ark)**2/4.0_ark
  !
  g_zz =   cos(rho*0.5_ark)**2/mY*(r1**2+r2**2)
  !
  re = extF%coef(3,4) 
  ae = extF%coef(4,4)*pi/180.0_ark
  !b0 = extF%coef(3,1)
  !
  y1 = (r1 - re)
  y2 = (r2 - re)
  y3 = cos(alpha) - cos(ae)
  !
  k = extF%nterms(4)
  !
  p(3:k-2) = extF%coef(5:k,4)
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
  muxz = v1+v2+v3+v4+v5+v6+v7+v8

  re = extF%coef(3,3) 
  ae = extF%coef(4,3)*pi/180.0_ark
  !b0 = extF%coef(3,1)
  !
  y1 = (r1 - re)
  y2 = (r2 - re)
  y3 = cos(alpha) - cos(ae)


    k = extF%nterms(3) 
    !
    q(5:k) = extF%coef(5:k,3)
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
    muzz =(v0+v1+v2+v3+v4+v5+v6+v7+v8)


  !
  gxx = .25_ark*(mX+mY)/cos(rho*0.5_ark)**2/mX/mY*(1.0_ark/r1**2+1.0_ark/r2**2)-.5_ark/cos(rho*0.5_ark)**2/mX/r1*r2
  !
  gyy = .25_ark*(mX+mY)/mX/mY*(1.0_ark/r1**2+1.0_ark/r2**2)-.5_ark*(2.0_ark*cos(rho*0.5_ark)**2-1.0_ark)/mX/r1*r2
  !
  g1y =  -.5_ark*sin(rho)/mX/r2
  g2y =   .5_ark*sin(rho)/mX/r1
  g3y =   .5_ark*sin(rho)*(mX+mY)/mX/mY*(1.0_ark/r1**2-1.0_ark/r2**2)
  !
  ! inverse tensor of inertia
  !
  mX = molec%atomMasses(1)
  mY = molec%atomMasses(2)
  ! I^-1 = [[a,0,b/sin(rho)],[0,c,0],[b/sin(rho),0,d/sin(rho/2)^2]]
  a = ((mY*(r1 - r2)**2 + mX*(r1**2 + r2**2))*(1.0_ark/cos(rho/2.0_ark))**2)/(4.0_ark*mY*mX*r1**2*r2**2)
  b = ((mY + mX)*(r1**2 - r2**2))/(2.0_ark*mY*mX*r1**2*r2**2)
  c = (2.0_ark*mY + mX)/(mY*((mY + mX)*(r1**2 + r2**2) + 2.0_ark*mY*r1*r2*cos(rho)))
  d = ((mY*(r1 + r2)**2 + mX*(r1**2 + r2**2)))/(4.0_ark*mY*mX*r1**2*r2**2)
  !
  ! 1,1
  g_out(1,0) =  g(2,2)/gyy*g1y !g(1,1)*a + g(1,3)*b*rho_over_sinrho
  g_out(1,-1) = 0
  g_out(1,-2) = 0
  !
  ! 1,3
  g_out(2,0)  = g(2,2)/gyy*g2y    !  muzz !
  g_out(2,-1) = 0  ! muxz !*Izz*g_zz*2.0_ark  ! gtxz ! g(1,3)*d*rho2_over_sin2rhohalf + g(1,1)*b*rho_over_sinrho
  g_out(2,-2) = 0
  !
  ! 2,2
  g_out(3,0)  = muzz/gyy*g3y   !g(2,2)*c
  g_out(3,-1) = 0
  g_out(3,-2) = 0
  !
  ! 3,1
  g_out(4,0) =  0 !g(3,1)*a*rho + g(3,3)*b*rho_over_sinrho
  g_out(4,-1) = 0
  g_out(4,-2) = 0
  !
  ! 3,3
  g_out(5,0) =  0 ! muzz
  g_out(5,-1) = 0 ! g(3,3)*d*rho2_over_sin2rhohalf + g(3,1)*b*rho2_over_sinrho
  g_out(5,-2) = 0
  !
  f = (/g_out(1,0), g_out(2,-1), g_out(3,0), -g_out(4,0), g_out(5,-1)/)
  ! we multiply g(1,3) and g(3,1) here with -1 because of the opposite direction
  ! of the bisector axis (x-axis) in TROVE and in the actual ab initio data 
  !
end subroutine prop_xy2_gtensor_bisector




! Rotational g-tensor for quasi-linear molecule, like H2O, in the bisector
! frame, where some of the elements become singular at linear geometry.

subroutine prop_xy2_gtensor_bisector_rank3(rank, ncoords, natoms, local, xyz, f)
  !
  implicit none 
  !
  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)
  integer(ik) :: iatom
  !
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha, rho, a, b, c, d
  real(ark) :: g(3,3), mat(3,3), g_out(5,-2:0), e1(3), e2(3), e3(3), x(natoms,3), mY, mX, mP
  real(ark) :: rho_over_sinrho, rho2_over_sinrho, rho2_over_sin2rhohalf,rho_over_sinrho_2,gtxz
  real(ark),parameter  :: rho_threshold = 0.01_rk
  integer(ik) :: icentre

  integer(ik)           :: k
  real(ark)             :: y1,y2,y3,re,ae,b0,muxz,muzz,Izz,g_zz,gxx,gyy,g1y,g2y,g3y
  real(ark)             :: v0,v1,v2,v3,v4,v5,v6,v7,v8,q(3:extF%nterms(3))
  real(ark), parameter :: muN = 5.050783699e-6 ! nuclear magneton in units of Debye
  !
  !
  if (rank/=3) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_gtensor_bisector: rank of the input tensor =', rank, ', expected 5'
    stop
  endif
  !
  ! xyz are undefined for the local case
  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_gtensor_bisector: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_gtensor_bisector - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  xyz0 = x(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = x(iatom,:) - xyz0(:)
  enddo
  !
  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))
  !
  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2
  !
  alpha = aacos(sum(e1*e2))
  !
  rho = pi-alpha
  !
  ! fitted tensor elements
  !
  g = 0

  !
  ! transform with the inverse inertia tensor
  !
  if (rho>rho_threshold) then
    !
    rho_over_sinrho = rho/sin(rho)
    rho_over_sinrho_2 = rho/sin(rho*0.5_ark)
    rho2_over_sinrho = rho**2/sin(rho)
    rho2_over_sin2rhohalf = rho**2/sin(rho*0.5_ark)**2
    !
  else
    !
    rho_over_sinrho = 1.0_ark + rho**2/6.0_ark + (7.0_ark*rho**4)/360.0_ark + (31.0_ark*rho**6)/15120.0_ark &
                    + (127.0_ark*rho**8)/604800.0_ark + (73.0_ark*rho**10)/3.42144e6_ark
    rho_over_sinrho_2 = 0.5_ark*rho-1._ark/48._ark*rho**3+1._ark/3840._ark*rho**5-1._ark/645120._ark*rho**7
    rho2_over_sinrho = rho + rho**3/6.0_ark + (7.0_ark*rho**5)/360.0_ark + (31.0_ark*rho**7)/15120.0_ark &
                     + (127._ark*rho**9)/604800.0_ark
    rho2_over_sin2rhohalf = 4.0_ark + rho**2/3.0_ark + rho**4/60.0_ark + rho**6/1512.0_ark &
                          + rho**8/43200.0_ark + rho**10/1.33056e6_ark
    !
  endif
  !
  mP = 1.007276_ark
  !
  mX = molec%atomMasses(1)
  mY = molec%atomMasses(2)
  !
  ! nuclear part of g[x,z] from Cybilski 1994 
  gtxz = -mP*4.0_ark*cos(rho/2.0_ark)*rho_over_sinrho_2*(&
            -11.0_ark*mY**2*r1**2+11.0_ark*mY**2*r2**2-4.0_ark*r1**2*mX**2&
            -8.0_ark*r1**2*mX*mY+4*r2**2*mX**2+8.0_ark*r2**2*mX*mY)/&
          ( (mX+2.0_ark*mY)*mY*(mX*r1**2+mX*r2**2+mY*r1**2-2.0_ark*mY*r1*r2+mY*r2**2) )

  !
  Izz = (mY*(mX*r2**2+mX*r1**2-2.0_ark*mY*r1*r2+mY*r1**2+mY*r2**2))/(mX+2.0_ark*mY)/cos(rho*0.5_ark)**2/4.0_ark
  !
  g_zz =   cos(rho*0.5_ark)**2/mY*(r1**2+r2**2)

  re = extF%coef(3,3) 
  ae = extF%coef(4,3)*pi/180.0_ark
  !b0 = extF%coef(3,1)
  !
  y1 = (r1 - re)
  y2 = (r2 - re)
  y3 = cos(alpha) - cos(ae)


    k = extF%nterms(3) 
    !
    q(5:k) = extF%coef(5:k,3)
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
    muzz =(v0+v1+v2+v3+v4+v5+v6+v7+v8)


  !
  gxx = .25_ark*(mX+mY)/cos(rho*0.5_ark)**2/mX/mY*(1.0_ark/r1**2+1.0_ark/r2**2)-.5_ark/cos(rho*0.5_ark)**2/mX/r1*r2
  !
  gyy = .25_ark*(mX+mY)/mX/mY*(1.0_ark/r1**2+1.0_ark/r2**2)-.5_ark*(2.0_ark*cos(rho*0.5_ark)**2-1.0_ark)/mX/r1*r2
  !
  g1y =  -.5_ark*sin(rho)/mX/r2
  g2y =   .5_ark*sin(rho)/mX/r1
  g3y =   .5_ark*sin(rho)*(mX+mY)/mX/mY*(1.0_ark/r1**2-1.0_ark/r2**2)
  !
  ! inverse tensor of inertia
  !
  mX = molec%atomMasses(1)
  mY = molec%atomMasses(2)
  ! I^-1 = [[a,0,b/sin(rho)],[0,c,0],[b/sin(rho),0,d/sin(rho/2)^2]]
  a = ((mY*(r1 - r2)**2 + mX*(r1**2 + r2**2))*(1.0_ark/cos(rho/2.0_ark))**2)/(4.0_ark*mY*mX*r1**2*r2**2)
  b = ((mY + mX)*(r1**2 - r2**2))/(2.0_ark*mY*mX*r1**2*r2**2)
  c = (2.0_ark*mY + mX)/(mY*((mY + mX)*(r1**2 + r2**2) + 2.0_ark*mY*r1*r2*cos(rho)))
  d = ((mY*(r1 + r2)**2 + mX*(r1**2 + r2**2)))/(4.0_ark*mY*mX*r1**2*r2**2)
  !
  ! 1,1
  g_out(1,0) =  0 !g(1,1)*a + g(1,3)*b*rho_over_sinrho
  g_out(1,-1) = 0
  g_out(1,-2) = 0
  !
  ! 1,3
  g_out(2,0)  = muzz/gyy*g3y    !  muzz !
  g_out(2,-1) = 0  ! muxz !*Izz*g_zz*2.0_ark  ! gtxz ! g(1,3)*d*rho2_over_sin2rhohalf + g(1,1)*b*rho_over_sinrho
  g_out(2,-2) = 0
  !
  ! 2,2
  g_out(3,0)  = 0   !g(2,2)*c
  g_out(3,-1) = 0
  g_out(3,-2) = 0
  !
  f = (/g_out(1,0), g_out(2,0), g_out(3,0)/)*muN*0.5_ark
  ! we multiply g(1,3) and g(3,1) here with -1 because of the opposite direction
  ! of the bisector axis (x-axis) in TROVE and in the actual ab initio data 
  !
end subroutine prop_xy2_gtensor_bisector_rank3



!###################################################################################################


! Spin-rotation tensor for quasi-linear molecule, like H2O, in the bisector
! frame, where some of the elements become singular at linear geometry.

subroutine prop_xy2_spin_rotation_bisector(rank, ncoords, natoms, local, xyz, f)
  !
  implicit none 
  !
  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)
  integer(ik) :: iatom
  !
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha, rho, a, b, c, d
  real(ark) :: csr(3,3), mat(3,3), c_out(5,-2:0), e1(3), e2(3), e3(3), x(natoms,3), mH, mO
  real(ark) :: rho_over_sinrho, rho2_over_sinrho, rho2_over_sin2rhohalf
  real(ark),parameter  :: rho_threshold = 0.01_rk
  integer(ik) :: icentre
  !
  if (rank/=5) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_spin_rotation_bisector: rank of the spin-rotation tensor =', rank, ', expected 5'
    stop
  endif
  !
  ! xyz are undefined for the local case
  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_spin_rotation_bisector: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_spin_rotation_bisector - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  xyz0 = x(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = x(iatom,:) - xyz0(:)
  enddo
  !
  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))
  !
  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2
  !
  alpha = aacos(sum(e1*e2))
  !
  rho = pi-alpha
  !
  icentre = extF%coef(2,1) ! = 1 or 2, determines on which Y-atom spin-rotation tensor is centered
  !
  ! fitted tensor elements
  !
  csr = 0
  !
  if (icentre==1) then
      csr(1,1) = fit_xy2_sr(extF%nterms(1)-2, extF%coef(3:extF%nterms(1),1), (/r1, r2, alpha/))
      csr(2,2) = fit_xy2_sr(extF%nterms(2)-2, extF%coef(3:extF%nterms(2),2), (/r1, r2, alpha/))
      csr(3,3) = fit_xy2_sr_rhopow_min_one(extF%nterms(3)-2, extF%coef(3:extF%nterms(3),3), (/r1, r2, alpha/))
      csr(1,3) = fit_xy2_sr_rhopow_min_one(extF%nterms(4)-2, extF%coef(3:extF%nterms(4),4), (/r1, r2, alpha/))
      csr(3,1) = fit_xy2_sr_rhopow_min_one(extF%nterms(5)-2, extF%coef(3:extF%nterms(5),5), (/r1, r2, alpha/))
  elseif (icentre==2) then
      csr(1,1) = fit_xy2_sr(extF%nterms(1)-2, extF%coef(3:extF%nterms(1),1), (/r2, r1, alpha/))
      csr(2,2) = fit_xy2_sr(extF%nterms(2)-2, extF%coef(3:extF%nterms(2),2), (/r2, r1, alpha/))
      csr(3,3) = fit_xy2_sr_rhopow_min_one(extF%nterms(3)-2, extF%coef(3:extF%nterms(3),3), (/r2, r1, alpha/))
      csr(1,3) = -fit_xy2_sr_rhopow_min_one(extF%nterms(4)-2, extF%coef(3:extF%nterms(4),4), (/r2, r1, alpha/))
      csr(3,1) = -fit_xy2_sr_rhopow_min_one(extF%nterms(5)-2, extF%coef(3:extF%nterms(5),5), (/r2, r1, alpha/))
  else
      write(out, '(a,1x,i3)') 'prop_xy2_spin_rotation_bisector error: icentre /= (1 or 2)'
      stop 'prop_xy2_spin_rotation_bisector - bad icentre parameter' 
  endif
  !
  ! transform with the inverse inertia tensor
  !
  if (rho>rho_threshold) then
    !
    rho_over_sinrho = rho/sin(rho)
    rho2_over_sinrho = rho**2/sin(rho)
    rho2_over_sin2rhohalf = rho**2/sin(rho*0.5_ark)**2
    !
  else
    !
    rho_over_sinrho = 1.0_ark + rho**2/6.0_ark + (7.0_ark*rho**4)/360.0_ark + (31.0_ark*rho**6)/15120.0_ark &
                    + (127.0_ark*rho**8)/604800.0_ark + (73.0_ark*rho**10)/3.42144e6_ark
    rho2_over_sinrho = rho + rho**3/6.0_ark + (7.0_ark*rho**5)/360.0_ark + (31.0_ark*rho**7)/15120.0_ark &
                     + (127.0_ark*rho**9)/604800.0_ark
    rho2_over_sin2rhohalf = 4.0_ark + rho**2/3.0_ark + rho**4/60.0_ark + rho**6/1512.0_ark &
                          + rho**8/43200.0_ark + rho**10/1.33056e6_ark
    !
  endif
  !
  ! inverse tensor of inertia
  !
  mO = molec%atomMasses(1) ! masses used in the fitting for H2O: mO = 15.994914630, mH = 1.007825035
  mH = molec%atomMasses(2)
  ! I^-1 = [[a,0,b/sin(rho)],[0,c,0],[b/sin(rho),0,d/sin(rho/2)^2]]
  a = ((mH*(r1 - r2)**2 + mO*(r1**2 + r2**2))*(1.0_ark/cos(rho/2.0_ark))**2)/(4.0_ark*mH*mO*r1**2*r2**2)
  b = ((mH + mO)*(r1**2 - r2**2))/(2.0_ark*mH*mO*r1**2*r2**2)
  c = (2*mH + mO)/(mH*((mH + mO)*(r1**2 + r2**2) + 2*mH*r1*r2*cos(rho)))
  d = ((mH*(r1 + r2)**2 + mO*(r1**2 + r2**2)))/(4.0_ark*mH*mO*r1**2*r2**2)
  !
  ! 1,1
  c_out(1,0) = csr(1,1)*a + csr(1,3)*b*rho_over_sinrho
  c_out(1,-1) = 0
  c_out(1,-2) = 0
  !
  ! 1,3
  c_out(2,0) = 0
  c_out(2,-1) = csr(1,3)*d*rho2_over_sin2rhohalf + csr(1,1)*b*rho_over_sinrho
  c_out(2,-2) = 0
  !
  ! 2,2
  c_out(3,0) = csr(2,2)*c
  c_out(3,-1) = 0
  c_out(3,-2) = 0
  !
  ! 3,1
  c_out(4,0) = csr(3,1)*a*rho + csr(3,3)*b*rho_over_sinrho
  c_out(4,-1) = 0
  c_out(4,-2) = 0
  !
  ! 3,3
  c_out(5,0) = 0
  c_out(5,-1) = csr(3,3)*d*rho2_over_sin2rhohalf + csr(3,1)*b*rho2_over_sinrho
  c_out(5,-2) = 0
  !
  f = (/c_out(1,0), c_out(2,-1), c_out(3,0), c_out(4,0), c_out(5,-1)/)
  !
end subroutine prop_xy2_spin_rotation_bisector


!###################################################################################################


function fit_xy2_sr(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams 
  real(ark), intent(in) :: params(nparams), coords(3)
  real(ark) :: f

  real(ark) :: y1, y2, y3, alpha, rho, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad ! obsolete
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha  = coords(3)
  rho    = pi - alpha
  y3     = rho !1.0_ark - cos(rho)

  f0 = params(4)
  f1 = xy2_func_n1_d6( (/y1,y2,y3/), params(5:22)  )  ! nparams = 18
  f2 = xy2_func_n2_d6( (/y1,y2,y3/), params(23:67) )  ! nparams = 45
  f3 = xy2_func_n3_d6( (/y1,y2,y3/), params(68:87) )  ! nparams = 20

  f = f0 + f1 + f2 + f3

end function fit_xy2_sr


!###################################################################################################


function fit_xy2_sr_rhopow_min_one(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams 
  real(ark), intent(in) :: params(nparams), coords(3)
  real(ark) :: f

  real(ark) :: y1, y2, y3, alpha, rho, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad ! obsolete
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha  = coords(3)
  rho    = pi - alpha
  y3     = rho

  f0 = 0
  f1 = xy2_func_n1_d6_rhopow_min_one( (/y1,y2,y3/), params(5:22)  )  ! nparams = 18
  f2 = xy2_func_n2_d6_rhopow_min_one( (/y1,y2,y3/), params(23:67) )  ! nparams = 45
  f3 = xy2_func_n3_d6_rhopow_min_one( (/y1,y2,y3/), params(68:87) )  ! nparams = 20

  f = f0 + f1 + f2 + f3

end function fit_xy2_sr_rhopow_min_one


!###################################################################################################


function fit_xy2_sr_A1(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(3)
  real(ark) :: f

  real(ark) :: y1, y2, y3, alpha, rho, rad, f0, f1, f2, f3, req, alphaeq, beta, rhoe

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad ! obsolete
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha  = coords(3)
  rho    = pi - alpha
  rhoe   = pi - alphaeq
  y3     = rho

  f0 = params(4)
  f1 = xy2_func_a1_n1_d6( (/y1,y2,y3/), params(5:16)  )  ! nparams = 12
  f2 = xy2_func_a1_n2_d6( (/y1,y2,y3/), params(17:40) )  ! nparams = 24
  f3 = xy2_func_a1_n3_d6( (/y1,y2,y3/), params(41:53) )  ! nparams = 13

  f = f0 + f1 + f2 + f3

end function fit_xy2_sr_A1


!###################################################################################################


function fit_xy2_sr_A1_rhopow_min_one(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(3)
  real(ark) :: f

  real(ark) :: y1, y2, y3, alpha, rho, rad, f0, f1, f2, f3, req, alphaeq, beta, rhoe

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad ! obsolete
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha  = coords(3)
  rho    = pi - alpha
  rhoe   = pi - alphaeq
  y3     = rho

  f0 = params(4)
  f1 = xy2_func_a1_n1_d6_rhopow_min_one( (/y1,y2,y3/), params(5:16)  )  ! nparams = 12
  f2 = xy2_func_a1_n2_d6_rhopow_min_one( (/y1,y2,y3/), params(17:40) )  ! nparams = 24
  f3 = xy2_func_a1_n3_d6_rhopow_min_one( (/y1,y2,y3/), params(41:53) )  ! nparams = 13

  f = f0 + f1 + f2 + f3

end function fit_xy2_sr_A1_rhopow_min_one


!###################################################################################################



function fit_xy2_sr_B2(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(3)
  real(ark) :: f

  real(ark) :: y1, y2, y3, alpha, rho, rad, f0, f1, f2, f3, req, alphaeq, beta, rhoe

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad ! obsolete
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha  = coords(3)
  rho    = pi - alpha
  rhoe   = pi - alphaeq
  y3     = rho

  f0 = params(4)
  f1 = xy2_func_b2_n1_d6((/y1,y2,y3/), params(5:10))    ! nparams = 6
  f2 = xy2_func_b2_n2_d6((/y1,y2,y3/), params(11:32))   ! nparams = 22
  f3 = xy2_func_b2_n3_d6((/y1,y2,y3/), params(33:44))   ! nparams = 12

  f = f0 + f1 + f2 + f3

end function fit_xy2_sr_B2





function fit_xy2_sr_B2_rhopow_min_one(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(3)
  real(ark) :: f

  real(ark) :: y1, y2, y3, alpha, rho, rad, f0, f1, f2, f3, req, alphaeq, beta, rhoe

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad ! obsolete
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha  = coords(3)
  rho    = pi - alpha
  rhoe   = pi - alphaeq
  y3     = rho

  f0 = params(4)
  f1 = xy2_func_b2_n1_d6_rhopow_min_one((/y1,y2,y3/), params(5:10))    ! nparams = 6
  f2 = xy2_func_b2_n2_d6_rhopow_min_one((/y1,y2,y3/), params(11:32))   ! nparams = 22
  f3 = xy2_func_b2_n3_d6_rhopow_min_one((/y1,y2,y3/), params(33:44))   ! nparams = 12

  f = f0 + f1 + f2 + f3

end function fit_xy2_sr_B2_rhopow_min_one


!###################################################################################################


! Spin-rotation tensor for non-linear molecule, like H2S, in the bisector frame.
! Care is not taken about the behaviour at near-linear geometries.

recursive subroutine prop_xy2_spin_rotation_bisector_nonlin(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom, icentre
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha1, e1(3), e2(3), coords(3), sr(3,3)

  icentre = nint(extf%coef(1,1)) ! icentre=1 means centre on Y1 and icentre=2 - on Y2

  if (rank/=9) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_spin_rotation_bisector_nonlin error: rank of the external field vector =', rank, &
      ', expected 9'
    stop
  endif

  xyz0 = xyz(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = xyz(iatom,:) - xyz0(:)
  enddo

  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))

  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2

  alpha1 = aacos(sum(e1*e2))

  if (abs(alpha1-pi)<0.0001) &
      stop 'prop_xy2_spin_rotation_bisector_nonlin error: valence bond angle is 180 degrees &
           (does not work for linear molecule)'

  if (icentre==1) then
    coords = (/r1,r2,alpha1/)
  elseif (icentre==2) then
    coords = (/r2,r1,alpha1/)
  else
    stop 'prop_xy2_spin_rotation_bisector_nonlin error: illegal value for "icentre" (can take values 1 or 2)'
  endif

  sr = 0
  sr(1,1) = fit_xy2_nosym(extF%nterms(1)-1, extF%coef(2:extF%nterms(1),1), coords) ! first parameter defines centre on Y1 or Y2, see above
  sr(2,2) = fit_xy2_nosym(extF%nterms(2), extF%coef(1:extF%nterms(2),2), coords)
  sr(3,3) = fit_xy2_nosym(extF%nterms(3), extF%coef(1:extF%nterms(3),3), coords)
  sr(1,3) = fit_xy2_nosym(extF%nterms(4), extF%coef(1:extF%nterms(4),4), coords)
  sr(3,1) = fit_xy2_nosym(extF%nterms(5), extF%coef(1:extF%nterms(5),5), coords)

  if (icentre==2) then
    sr(1,3) = -sr(1,3)
    sr(3,1) = -sr(3,1)
  endif

  ! the order of elements:
  ! 11 ==> zz
  ! 22 ==> yy
  ! 33 ==> xx
  ! 13 ==> zx
  ! 31 ==> xz

  !               xx     xy      xz        yx        yy     yz       zx       zy       zz
  f(1:9) = (/ sr(3,3), 0.0_ark, sr(3,1), 0.0_ark, sr(2,2), 0.0_ark, sr(1,3), 0.0_ark, sr(1,1) /)

end subroutine prop_xy2_spin_rotation_bisector_nonlin


!###################################################################################################


! Subroutine to test prop_xy2_spin_rotation_bisector_nonline, if it is able to reproduce
! the original ab initio Cartesian components

subroutine TEST_prop_xy2_spin_rotation_bisector_nonlin(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iounit, iatom, npoints, ipoint, info, ithread, icentre, i, ii
  real(ark) :: xyz_(natoms,3), sr_inp(9), sr_calc(9), rms(9), enr, sr(3,3)
  character(cl) :: fname
  integer(ik), external :: omp_get_thread_num

  ithread = omp_get_thread_num()
  print*, ithread
  icentre = nint(extf%coef(1,1))

  if (ithread==0) then

  write(out, '(/a)') 'TEST_prop_xy2_spin_rotation_bisector_nonlin/start: test spin-rotation tensor transformation'

  if (icentre==1) then
    fname = 'SR1.dat'
  elseif (icentre==2) then
    fname = 'SR2.dat'
  endif

  call IOstart(fname, iounit)

  open(iounit,form='formatted',file=fname,iostat=info)
  if (info/=0) then
    write(out, '(/a,a,a)') 'TEST_prop_xy2_spin_rotation_bisector_nonlin error: file "', trim(fname), '" not found'
    stop
  endif

  ipoint = 0
  rms = 0
  do
    read(iounit,*,iostat=info) (xyz_(iatom,1:3), iatom=1, natoms), enr, (sr(i,1:3), i=1,3)
    !                   xx      xy      xz        yx        yy     yz       zx       zy       zz
    sr_inp(1:9) = (/ sr(3,3), sr(3,2), sr(3,1), sr(2,3), sr(2,2), sr(2,1), sr(1,3), sr(1,2), sr(1,1) /)
    if (info/=0) exit
    if (enr>10000) cycle
    ipoint = ipoint + 1
    call prop_xy2_spin_rotation_bisector_nonlin(rank, ncoords, natoms, local, xyz_, sr_calc(1:9))
    write(out, '(1x,i6,1x,f10.4,1x,9(5x,f12.4,1x,f12.4,1x,f12.4))') &
      ipoint, enr, (sr_inp(i), sr_calc(i), sr_inp(i)-sr_calc(i), i=1, 9)
    rms = rms + (sr_inp-sr_calc)**2
  enddo
  npoints = ipoint
  close(iounit)

  rms = sqrt(rms/real(npoints,ark))

  write(out, '(/1x,a,9(1x,f12.4))') 'rms =', rms

  write(out, '(/a)') 'TEST_prop_xy2_spin_rotation_bisector_nonlin/done'

  endif
  f = 0

  !$omp barrier
  stop

end subroutine TEST_prop_xy2_spin_rotation_bisector_nonlin


!###################################################################################################


function fit_xy2_nosym(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(:)
  real(ark) :: f

  real(ark) :: r1, r2, alpha1, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = real(pi,ark)/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad
  beta    = params(3)

  r1     = (coords(1)-req) * exp(-beta*(coords(1)-req)**2)
  r2     = (coords(2)-req) * exp(-beta*(coords(2)-req)**2)
  alpha1 = coords(3)-alphaeq

  f0 = params(4)
  f1 = xy2_func_n1_d6( (/r1,r2,alpha1/), params(5:22)  )  ! nparams = 18
  f2 = xy2_func_n2_d6( (/r1,r2,alpha1/), params(23:67) )  ! nparams = 45
  f3 = xy2_func_n3_d6( (/r1,r2,alpha1/), params(68:87) )  ! nparams = 20

  f = f0 + f1 + f2 + f3

end function fit_xy2_nosym


!###################################################################################################


function xy2_func_n1_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(18)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(18)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1
f(2) = r1**2
f(3) = r1**3
f(4) = r1**4
f(5) = r1**5
f(6) = r1**6
f(7) = r2
f(8) = r2**2
f(9) = r2**3
f(10) = r2**4
f(11) = r2**5
f(12) = r2**6
f(13) = a1
f(14) = a1**2
f(15) = a1**3
f(16) = a1**4
f(17) = a1**5
f(18) = a1**6
v = sum(f*params)
end function xy2_func_n1_d6


!###################################################################################################


function xy2_func_n2_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(45)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(45)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2
f(2) = r1*r2**2
f(3) = r1*r2**3
f(4) = r1*r2**4
f(5) = r1*r2**5
f(6) = r1**2*r2
f(7) = r1**2*r2**2
f(8) = r1**2*r2**3
f(9) = r1**2*r2**4
f(10) = r1**3*r2
f(11) = r1**3*r2**2
f(12) = r1**3*r2**3
f(13) = r1**4*r2
f(14) = r1**4*r2**2
f(15) = r1**5*r2
f(16) = a1*r1
f(17) = a1**2*r1
f(18) = a1**3*r1
f(19) = a1**4*r1
f(20) = a1**5*r1
f(21) = a1*r1**2
f(22) = a1**2*r1**2
f(23) = a1**3*r1**2
f(24) = a1**4*r1**2
f(25) = a1*r1**3
f(26) = a1**2*r1**3
f(27) = a1**3*r1**3
f(28) = a1*r1**4
f(29) = a1**2*r1**4
f(30) = a1*r1**5
f(31) = a1*r2
f(32) = a1**2*r2
f(33) = a1**3*r2
f(34) = a1**4*r2
f(35) = a1**5*r2
f(36) = a1*r2**2
f(37) = a1**2*r2**2
f(38) = a1**3*r2**2
f(39) = a1**4*r2**2
f(40) = a1*r2**3
f(41) = a1**2*r2**3
f(42) = a1**3*r2**3
f(43) = a1*r2**4
f(44) = a1**2*r2**4
f(45) = a1*r2**5
v = sum(f*params)
end function xy2_func_n2_d6


!###################################################################################################


function xy2_func_n3_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(20)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(20)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = a1*r1*r2
f(2) = a1**2*r1*r2
f(3) = a1**3*r1*r2
f(4) = a1**4*r1*r2
f(5) = a1*r1*r2**2
f(6) = a1**2*r1*r2**2
f(7) = a1**3*r1*r2**2
f(8) = a1*r1*r2**3
f(9) = a1**2*r1*r2**3
f(10) = a1*r1*r2**4
f(11) = a1*r1**2*r2
f(12) = a1**2*r1**2*r2
f(13) = a1**3*r1**2*r2
f(14) = a1*r1**2*r2**2
f(15) = a1**2*r1**2*r2**2
f(16) = a1*r1**2*r2**3
f(17) = a1*r1**3*r2
f(18) = a1**2*r1**3*r2
f(19) = a1*r1**3*r2**2
f(20) = a1*r1**4*r2
v = sum(f*params)
end function xy2_func_n3_d6


!###################################################################################################


function xy2_func_n1_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(18)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(18)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = 0
f(2) = 0
f(3) = 0
f(4) = 0
f(5) = 0
f(6) = 0
f(7) = 0
f(8) = 0
f(9) = 0
f(10) = 0
f(11) = 0
f(12) = 0
f(13) = 1.0_ark
f(14) = a1**1
f(15) = a1**2
f(16) = a1**3
f(17) = a1**4
f(18) = a1**5
v = sum(f*params)
end function xy2_func_n1_d6_rhopow_min_one


!###################################################################################################


function xy2_func_n2_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(45)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(45)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = 0
f(2) = 0
f(3) = 0
f(4) = 0
f(5) = 0
f(6) = 0
f(7) = 0
f(8) = 0
f(9) = 0
f(10) = 0
f(11) = 0
f(12) = 0
f(13) = 0
f(14) = 0
f(15) = 0
f(16) = r1
f(17) = a1**1*r1
f(18) = a1**2*r1
f(19) = a1**3*r1
f(20) = a1**4*r1
f(21) = r1**2
f(22) = a1**1*r1**2
f(23) = a1**2*r1**2
f(24) = a1**3*r1**2
f(25) = r1**3
f(26) = a1**1*r1**3
f(27) = a1**2*r1**3
f(28) = r1**4
f(29) = a1**1*r1**4
f(30) = r1**5
f(31) = r2
f(32) = a1**1*r2
f(33) = a1**2*r2
f(34) = a1**3*r2
f(35) = a1**4*r2
f(36) = r2**2
f(37) = a1**1*r2**2
f(38) = a1**2*r2**2
f(39) = a1**3*r2**2
f(40) = r2**3
f(41) = a1**1*r2**3
f(42) = a1**2*r2**3
f(43) = r2**4
f(44) = a1**1*r2**4
f(45) = r2**5
v = sum(f*params)
end function xy2_func_n2_d6_rhopow_min_one


!###################################################################################################


function xy2_func_n3_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(20)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(20)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2
f(2) = a1**1*r1*r2
f(3) = a1**2*r1*r2
f(4) = a1**3*r1*r2
f(5) = r1*r2**2
f(6) = a1**1*r1*r2**2
f(7) = a1**2*r1*r2**2
f(8) = r1*r2**3
f(9) = a1**1*r1*r2**3
f(10) = r1*r2**4
f(11) = r1**2*r2
f(12) = a1**1*r1**2*r2
f(13) = a1**2*r1**2*r2
f(14) = r1**2*r2**2
f(15) = a1**1*r1**2*r2**2
f(16) = r1**2*r2**3
f(17) = r1**3*r2
f(18) = a1**1*r1**3*r2
f(19) = r1**3*r2**2
f(20) = r1**4*r2
v = sum(f*params)
end function xy2_func_n3_d6_rhopow_min_one


!###################################################################################################


function xy2_func_a1_n1_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(12)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(12)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1 + r2
f(2) = r1**2 + r2**2
f(3) = r1**3 + r2**3
f(4) = r1**4 + r2**4
f(5) = r1**5 + r2**5
f(6) = r1**6 + r2**6
f(7) = a1
f(8) = a1**2
f(9) = a1**3
f(10) = a1**4
f(11) = a1**5
f(12) = a1**6
v = sum(f*params)
end function xy2_func_a1_n1_d6


!###################################################################################################


function xy2_func_a1_n2_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(24)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(24)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2
f(2) = r1*r2*(r1 + r2)
f(3) = r1*r2*(r1**2 + r2**2)
f(4) = r1*r2*(r1**3 + r2**3)
f(5) = r1*r2*(r1**4 + r2**4)
f(6) = r1**2*r2**2
f(7) = r1**2*r2**2*(r1 + r2)
f(8) = r1**2*r2**2*(r1**2 + r2**2)
f(9) = r1**3*r2**3
f(10) = a1*(r1 + r2)
f(11) = a1**2*(r1 + r2)
f(12) = a1**3*(r1 + r2)
f(13) = a1**4*(r1 + r2)
f(14) = a1**5*(r1 + r2)
f(15) = a1*(r1**2 + r2**2)
f(16) = a1**2*(r1**2 + r2**2)
f(17) = a1**3*(r1**2 + r2**2)
f(18) = a1**4*(r1**2 + r2**2)
f(19) = a1*(r1**3 + r2**3)
f(20) = a1**2*(r1**3 + r2**3)
f(21) = a1**3*(r1**3 + r2**3)
f(22) = a1*(r1**4 + r2**4)
f(23) = a1**2*(r1**4 + r2**4)
f(24) = a1*(r1**5 + r2**5)
v = sum(f*params)
end function xy2_func_a1_n2_d6


!###################################################################################################


function xy2_func_a1_n3_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(13)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(13)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = a1*r1*r2
f(2) = a1**2*r1*r2
f(3) = a1**3*r1*r2
f(4) = a1**4*r1*r2
f(5) = a1*r1*r2*(r1 + r2)
f(6) = a1**2*r1*r2*(r1 + r2)
f(7) = a1**3*r1*r2*(r1 + r2)
f(8) = a1*r1*r2*(r1**2 + r2**2)
f(9) = a1**2*r1*r2*(r1**2 + r2**2)
f(10) = a1*r1*r2*(r1**3 + r2**3)
f(11) = a1*r1**2*r2**2
f(12) = a1**2*r1**2*r2**2
f(13) = a1*r1**2*r2**2*(r1 + r2)
v = sum(f*params)
end function xy2_func_a1_n3_d6


!###################################################################################################


function xy2_func_b2_n1_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(6)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(6)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = -r1 + r2
f(2) = r1**2 - r2**2
f(3) = r1**3 - r2**3
f(4) = -r1**4 + r2**4
f(5) = r1**5 - r2**5
f(6) = -r1**6 + r2**6
v = sum(f*params)
end function xy2_func_b2_n1_d6


!###################################################################################################


function xy2_func_b2_n2_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(22)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(22)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2*(r1 - r2)
f(2) = r1*r2*(-r1**2 + r2**2)
f(3) = r1*r2*(r1**3 - r2**3)
f(4) = r1*r2*(r1**4 - r2**4)
f(5) = r1**2*r2**2*(r1 - r2)
f(6) = r1**2*r2**2*(r1**2 - r2**2)
f(7) = r1**2*r2**2*(-r1**2 + r2**2)
f(8) = a1*(-r1 + r2)
f(9) = a1**2*(-r1 + r2)
f(10) = a1**3*(r1 - r2)
f(11) = a1**4*(r1 - r2)
f(12) = a1**5*(r1 - r2)
f(13) = a1*(-r1**2 + r2**2)
f(14) = a1**2*(-r1**2 + r2**2)
f(15) = a1**3*(-r1**2 + r2**2)
f(16) = a1**4*(r1**2 - r2**2)
f(17) = a1*(-r1**3 + r2**3)
f(18) = a1**2*(r1**3 - r2**3)
f(19) = a1**3*(-r1**3 + r2**3)
f(20) = a1*(r1**4 - r2**4)
f(21) = a1**2*(-r1**4 + r2**4)
f(22) = a1*(r1**5 - r2**5)
v = sum(f*params)
end function xy2_func_b2_n2_d6


!###################################################################################################


function xy2_func_b2_n3_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(12)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(12)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = a1*r1*r2*(r1 - r2)
f(2) = a1**2*r1*r2*(r1 - r2)
f(3) = a1**3*r1*r2*(r1 - r2)
f(4) = a1*r1*r2*(r1**2 - r2**2)
f(5) = a1**2*r1*r2*(r1**2 - r2**2)
f(6) = a1*r1*r2*(-r1**3 + r2**3)
f(7) = a1*r1*r2*(-r1 + r2)
f(8) = a1**2*r1*r2*(-r1 + r2)
f(9) = a1**3*r1*r2*(-r1 + r2)
f(10) = a1*r1**2*r2**2*(-r1 + r2)
f(11) = a1*r1*r2*(-r1**2 + r2**2)
f(12) = a1**2*r1*r2*(-r1**2 + r2**2)
v = sum(f*params)
end function xy2_func_b2_n3_d6


!###################################################################################################


function xy2_func_a1_n1_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(12)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(12)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1 + r2
f(2) = r1**2 + r2**2
f(3) = r1**3 + r2**3
f(4) = r1**4 + r2**4
f(5) = r1**5 + r2**5
f(6) = r1**6 + r2**6
f(7) = 1.0_ark
f(8) = a1**1
f(9) = a1**2
f(10) = a1**3
f(11) = a1**4
f(12) = a1**5
v = sum(f*params)
end function xy2_func_a1_n1_d6_rhopow_min_one


!###################################################################################################


function xy2_func_a1_n2_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(24)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(24)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2
f(2) = r1*r2*(r1 + r2)
f(3) = r1*r2*(r1**2 + r2**2)
f(4) = r1*r2*(r1**3 + r2**3)
f(5) = r1*r2*(r1**4 + r2**4)
f(6) = r1**2*r2**2
f(7) = r1**2*r2**2*(r1 + r2)
f(8) = r1**2*r2**2*(r1**2 + r2**2)
f(9) = r1**3*r2**3
f(10) = (r1 + r2)
f(11) = a1**1*(r1 + r2)
f(12) = a1**2*(r1 + r2)
f(13) = a1**3*(r1 + r2)
f(14) = a1**4*(r1 + r2)
f(15) = (r1**2 + r2**2)
f(16) = a1**1*(r1**2 + r2**2)
f(17) = a1**2*(r1**2 + r2**2)
f(18) = a1**3*(r1**2 + r2**2)
f(19) = (r1**3 + r2**3)
f(20) = a1**1*(r1**3 + r2**3)
f(21) = a1**2*(r1**3 + r2**3)
f(22) = (r1**4 + r2**4)
f(23) = a1**1*(r1**4 + r2**4)
f(24) = (r1**5 + r2**5)
v = sum(f*params)
end function xy2_func_a1_n2_d6_rhopow_min_one


!###################################################################################################


function xy2_func_a1_n3_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(13)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(13)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2
f(2) = a1**1*r1*r2
f(3) = a1**2*r1*r2
f(4) = a1**3*r1*r2
f(5) = r1*r2*(r1 + r2)
f(6) = a1**1*r1*r2*(r1 + r2)
f(7) = a1**2*r1*r2*(r1 + r2)
f(8) = r1*r2*(r1**2 + r2**2)
f(9) = a1**1*r1*r2*(r1**2 + r2**2)
f(10) = r1*r2*(r1**3 + r2**3)
f(11) = r1**2*r2**2
f(12) = a1**1*r1**2*r2**2
f(13) = r1**2*r2**2*(r1 + r2)
v = sum(f*params)
end function xy2_func_a1_n3_d6_rhopow_min_one


!###################################################################################################


function xy2_func_b2_n1_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(6)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(6)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = -r1 + r2
f(2) = r1**2 - r2**2
f(3) = r1**3 - r2**3
f(4) = -r1**4 + r2**4
f(5) = r1**5 - r2**5
f(6) = -r1**6 + r2**6
v = sum(f*params)
end function xy2_func_b2_n1_d6_rhopow_min_one


!###################################################################################################


function xy2_func_b2_n2_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(22)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(22)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2*(r1 - r2)
f(2) = r1*r2*(-r1**2 + r2**2)
f(3) = r1*r2*(r1**3 - r2**3)
f(4) = r1*r2*(r1**4 - r2**4)
f(5) = r1**2*r2**2*(r1 - r2)
f(6) = r1**2*r2**2*(r1**2 - r2**2)
f(7) = r1**2*r2**2*(-r1**2 + r2**2)
f(8) = (-r1 + r2)
f(9) = a1**1*(-r1 + r2)
f(10) = a1**2*(r1 - r2)
f(11) = a1**3*(r1 - r2)
f(12) = a1**4*(r1 - r2)
f(13) = (-r1**2 + r2**2)
f(14) = a1**1*(-r1**2 + r2**2)
f(15) = a1**2*(-r1**2 + r2**2)
f(16) = a1**3*(r1**2 - r2**2)
f(17) = (-r1**3 + r2**3)
f(18) = a1**1*(r1**3 - r2**3)
f(19) = a1**2*(-r1**3 + r2**3)
f(20) = (r1**4 - r2**4)
f(21) = a1**1*(-r1**4 + r2**4)
f(22) = (r1**5 - r2**5)
v = sum(f*params)
end function xy2_func_b2_n2_d6_rhopow_min_one


!###################################################################################################


function xy2_func_b2_n3_d6_rhopow_min_one(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(12)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(12)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2*(r1 - r2)
f(2) = a1**1*r1*r2*(r1 - r2)
f(3) = a1**2*r1*r2*(r1 - r2)
f(4) = r1*r2*(r1**2 - r2**2)
f(5) = a1**1*r1*r2*(r1**2 - r2**2)
f(6) = r1*r2*(-r1**3 + r2**3)
f(7) = r1*r2*(-r1 + r2)
f(8) = a1**1*r1*r2*(-r1 + r2)
f(9) = a1**2*r1*r2*(-r1 + r2)
f(10) = r1**2*r2**2*(-r1 + r2)
f(11) = r1*r2*(-r1**2 + r2**2)
f(12) = a1**1*r1*r2*(-r1**2 + r2**2)
v = sum(f*params)
end function xy2_func_b2_n3_d6_rhopow_min_one


!###################################################################################################


subroutine matrix_pseudoinverse_ark(m, n, mat, invmat, info)

  integer(ik), intent(in) :: m, n
  real(ark), intent(in) :: mat(m,n)
  real(ark), intent(out) :: invmat(n,m)
  integer(ik), intent(out), optional :: info

  integer(ik) lwork, info_, i, j
  double precision work1(1), matd(m,n), matu(m,m), matvt(n,n), invmatd(n,m), mat_d(n,m), sv(n), &
                   tmat(n,m), tol
  double precision, allocatable :: work(:)

  tol = 1.0d-08 !epsilon(1.0d0)

  matd = dble(mat)

  lwork = -1
  call dgesvd('A', 'A', m, n, matd, m, sv, matu, m, matvt, n, work1, lwork, info_)
  lwork = int(work1(1), kind=ik)

  allocate(work(lwork), stat=info_)
  if (info_/=0) then
    write(out, '(/a,1x,i6)') &
      'matrix_pseudoinverse_ark error: failed to allocate workspace array required for SVD, size =', &
      lwork
    stop
  endif

  call dgesvd('A', 'A', m, n, matd, m, sv, matu, m, matvt, n, work, lwork, info_)

  if (info_/=0) then
    write(out, '(/a,1x,i6)') 'matrix_pseudoinverse_ark error: SVD failed, info =', info_
    stop
  endif

  if (present(info)) info = 0

  mat_d = 0.0
  do i=1, n
    if (sv(i)>=tol) then
      mat_d(i,i) = 1.0d0/sv(i)
    else
      if (present(info)) then
        info = i
        mat_d(i,i) = 0.0
      else
        write(out, '(/a)') 'matrix_pseudoinverse_ark error: matrix is singular'
        write(out, '(a)') 'matrix:'
        do j=1, m
          write(out, '(<n>(1x,f10.6))') mat(j,1:n)
        enddo
        write(out, '(a)') 'singular elements:'
        do j=1, n
          write(out, '(1x,f)') sv(j)
        enddo
        stop 'STOP'
      endif
    endif
  enddo

  call dgemm('N', 'T', n, m, m, 1.0d0, mat_d(1:n,1:m), n, matu(1:m,1:m), m, 0.0d0, tmat(1:n,1:m), n)
  call dgemm('T', 'N', n, m, n, 1.0d0, matvt(1:n,1:n), n, tmat(1:n,1:m), n, 0.0d0, invmatd(1:n,1:m), n)

  invmat = real(invmatd,kind=ark)

  deallocate(work)

end subroutine matrix_pseudoinverse_ark

end module prop_xy2_spinrot

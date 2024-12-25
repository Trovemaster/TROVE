!
!  This unit defines all kinetic routines for a triatomic molecule of the XY2 type
!
module kin_xy2
  use accuracy
  use moltype

  implicit none

  public MLkinetic_xy2_bisect_EKE,MLkinetic_xyz_bisect_EKE,MLkinetic_xy2_bisect_EKE_sinrho,&
         MLkinetic_xy2_Radau_bisect_EKE,MLkinetic_xyz_EKE_sinrho,MLkinetic_xyz_bond_EKE,MLkinetic_xyz_bond_EKE_r2,&
         MLkinetic_xyz_Radau_EKE,MLkinetic_compact_xy2_bisect_EKE_rigid,MLkinetic_compact_xyz_alpha_bond2_EKE_rigid
  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains



  !
  ! Defining kinetic energy function 
  ! This is EKE generated using Maple for a bisecting frame with bond-length-angle and 
  ! removed singularity by combining U rho with dG/drho and multiplying muzz by rho^2
  ! and muxz and muyz by rho. The vecinity of zero (singularity) is expanded wrt rho, up to the 7th order 
  !
  subroutine MLkinetic_xy2_bisect_EKE(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY,rho_2
   real(ark),parameter  :: rho_threshold = 0.01_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_xy2_bisect_EKE-error: can be used with non-rigid case only')")
       stop 'MLkinetic_xy2_bisect_EKE can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     !
     rho_2 = rho*0.5_ark
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  (mX+mY)/mX/mY
     g_vib(1,2,1) =  -cos(rho)/mX
     g_vib(1,3,2) =  sin(rho)/mX
     g_vib(2,1,1) =  -cos(rho)/mX
     g_vib(2,2,1) =  (mX+mY)/mX/mY
     g_vib(2,3,3) =  sin(rho)/mX
     g_vib(3,1,2) =  sin(rho)/mX
     g_vib(3,2,3) =  sin(rho)/mX
     g_vib(3,3,4) =  (mX+mY)/mX/mY
     g_vib(3,3,5) =  2.0_ark*cos(rho)/mX
     g_vib(3,3,6) =  (mX+mY)/mX/mY
     !
     g_rot(1,1,4) =  .25_ark*(mX+mY)/cos(rho_2)**2/mX/mY
     g_rot(1,1,5) =  -.5_ark/cos(rho_2)**2/mX
     g_rot(1,1,6) =  .25_ark*(mX+mY)/cos(rho_2)**2/mX/mY
     g_rot(2,2,4) =  .25_ark*(mX+mY)/mX/mY
     g_rot(2,2,5) =  -.5_ark*(2.0_ark*cos(rho_2)**2-1.0_ark)/mX
     g_rot(2,2,6) =  .25_ark*(mX+mY)/mX/mY
     !
     g_cor(1,2,2) =  -.5_ark*sin(rho)/mX
     g_cor(2,2,3) =  .5_ark*sin(rho)/mX
     g_cor(3,2,4) =  -.5_ark*(mX+mY)/mX/mY
     g_cor(3,2,6) =  .5_ark*(mX+mY)/mX/mY
     !
     if (rho>rho_threshold) then
        !
        g_rot(1,3,4) =  .25_ark*rho*(mX+mY)/cos(rho_2)/sin(rho_2)/mX/mY 
        g_rot(1,3,6) =  -.25_ark*rho*(mX+mY)/cos(rho_2)/sin(rho_2)/mX/mY 
        g_rot(3,1,4) =  .25_ark*rho*(mX+mY)/cos(rho_2)/sin(rho_2)/mX/mY
        g_rot(3,1,6) =  -.25_ark*rho*(mX+mY)/cos(rho_2)/sin(rho_2)/mX/mY
        g_rot(3,3,4) =  .25_ark*rho**2*(mX+mY)/sin(rho_2)**2/mX/mY
        g_rot(3,3,5) =  .5_ark*rho**2/mX/sin(rho_2)**2
        g_rot(3,3,6) =  .25_ark*rho**2*(mX+mY)/sin(rho_2)**2/mX/mY
        !
        pseudo(4) =  .125_ark*(-mX*cos(rho)**2-mY*cos(rho)**2+rho**2*mY*cos(rho)**2+rho**2*mX*cos(rho)**2-2.0_ark *rho**2*mX+mX-&
                     2.0_ark*rho**2*mY+mY)/sin(rho)**2/mX/rho/mY 
        pseudo(5) =  -.250_ark *(cos(rho)**3+rho**2*cos(rho)**3+2._ark *rho*sin(rho)*cos(rho)**2-cos(rho)-&
                     2._ark *rho*sin(rho))/sin(rho)**2/mX/rho 
        pseudo(6) =  .125_ark*(-mX*cos(rho)**2-mY*cos(rho)**2+rho**2*mY*cos(rho)**2+rho**2*mX*cos(rho)**2-2.0_ark*rho**2*mX+mX-&
                     2.0_ark*rho**2*mY+mY)/sin(rho)**2/mX/rho/mY 
        !
     else
        !
        ! expansion around rho=0
        !
        g_rot(1,3,4) = (2._ark*mX+2._ark*mY)/mX/mY/4._ark+(mX/3._ark+mY/3._ark)/mX/mY*rho**2/4._ark+(7._ark/180._ark*mX+&
                        7._ark/180._ark*mY)/mX/mY*rho**4/4._ark+(31._ark/7560._ark*mX+31._ark/7560._ark*mY)/mX/mY*rho**6/4._ark
        g_rot(1,3,6) = -(2._ark*mX+2._ark*mY)/mX/mY/4._ark-(mX/3._ark+mY/3._ark)/mX/mY*rho**2/4._ark-(7._ark/180._ark*mX+&
                        7._ark/180._ark*mY)/mX/mY*rho**4/4._ark-(31._ark/7560._ark*mX+31._ark/7560._ark*mY)/mX/mY*rho**6/4._ark
        g_rot(3,1,4) = (2._ark*mX+2._ark*mY)/mX/mY/4._ark+(mX/3._ark+mY/3._ark)/mX/mY*rho**2/4._ark+(7._ark/180._ark*mX+&
                        7._ark/180._ark*mY)/mX/mY*rho**4/4._ark+(31._ark/7560._ark*mX+31._ark/7560._ark*mY)/mX/mY*rho**6/4._ark
        g_rot(3,1,6) = -(2._ark*mX+2._ark*mY)/mX/mY/4._ark-(mX/3._ark+mY/3._ark)/mX/mY*rho**2/4._ark-(7._ark/180._ark*mX+&
                        7._ark/180._ark*mY)/mX/mY*rho**4/4._ark-(31._ark/7560._ark*mX+31._ark/7560._ark*mY)/mX/mY*rho**6/4._ark
        g_rot(3,3,4) = -(-4._ark*mX-4._ark*mY)/mX/mY/4._ark-(-mX/3._ark-mY/3._ark)/mX/mY*rho**2/4._ark-&
                        (-mX/60._ark-mY/60._ark)/mX/mY*rho**4/4._ark
        g_rot(3,3,5) = 2._ark/mX+1._ark/mX*rho**2/6._ark+1._ark/mX*rho**4/120._ark
        g_rot(3,3,6) = -(-4._ark*mX-4._ark*mY)/mX/mY/4._ark-(-mX/3._ark-mY/3._ark)/mX/mY*rho**2/4._ark-&
                       (-mX/60._ark-mY/60._ark)/mX/mY*rho**4/4._ark
        !
        pseudo(4) = -(4._ark/3._ark*mY+4._ark/3._ark*mX)/mX/mY*rho/8._ark-(mY/15._ark+mX/15._ark)/mX/mY*rho**3/8._ark
        pseudo(5) =  2._ark/3._ark/mX*rho-11._ark/60._ark/mX*rho**3
        pseudo(6) = -(4._ark/3._ark*mY+4._ark/3._ark*mX)/mX/mY*rho/8._ark-(mY/15._ark+mX/15._ark)/mX/mY*rho**3/8._ark
        !
     endif

     !
   end subroutine  MLkinetic_xy2_bisect_EKE





  !
  ! Defining kinetic energy function 
  ! This is EKE generated using Maple for a bisecting frame with bond-length-angle and 
  ! removed singularity by combining U sin(rho) with dG/drho and multiplying muzz by sin(rho)^2
  ! and muxz and muyz by rho. The vecinity of zero (singularity) is expanded wrt rho, up to the 7th order 
  !
  subroutine MLkinetic_xy2_bisect_EKE_sinrho(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY,rho_2
   real(ark),parameter  :: rho_threshold = 0.01_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_xy2_bisect_EKE_sinrho-error: can be used with non-rigid case only')")
       stop 'MLkinetic_xy2_bisect_EKE_sinrho can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     !
     rho_2 = rho*0.5_ark
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  (mX+mY)/mX/mY
     g_vib(1,2,1) =  -cos(rho)/mX
     g_vib(1,3,2) =  sin(rho)/mX
     g_vib(2,1,1) =  -cos(rho)/mX
     g_vib(2,2,1) =  (mX+mY)/mX/mY
     g_vib(2,3,3) =  sin(rho)/mX
     g_vib(3,1,2) =  sin(rho)/mX
     g_vib(3,2,3) =  sin(rho)/mX
     g_vib(3,3,4) =  (mX+mY)/mX/mY
     g_vib(3,3,5) =  2.0_ark*cos(rho)/mX
     g_vib(3,3,6) =  (mX+mY)/mX/mY
     !
     g_rot(1,1,4) =  .25_ark*(mX+mY)/cos(rho_2)**2/mX/mY
     g_rot(1,1,5) =  -.5_ark/cos(rho_2)**2/mX
     g_rot(1,1,6) =  .25_ark*(mX+mY)/cos(rho_2)**2/mX/mY
     g_rot(2,2,4) =  .25_ark*(mX+mY)/mX/mY
     g_rot(2,2,5) =  -.5_ark*(2.0_ark*cos(rho_2)**2-1.0_ark)/mX
     g_rot(2,2,6) =  .25_ark*(mX+mY)/mX/mY
     !
     g_cor(1,2,2) =  -.5_ark*sin(rho)/mX
     g_cor(2,2,3) =  .5_ark*sin(rho)/mX
     g_cor(3,2,4) =  -.5_ark*(mX+mY)/mX/mY
     g_cor(3,2,6) =  .5_ark*(mX+mY)/mX/mY
     !
     pseudo(5) =  sin(rho)*cos(rho)/mX
     !
     g_rot(1,3,4) = -0.5_ark*(mX+mY)/(mX*mY)
     g_rot(1,3,6) = -0.5_ark*(mX+mY)/(mX*mY)
     g_rot(3,1,4) = -0.5_ark*(mX+mY)/(mX*mY)
     g_rot(3,1,6) = -0.5_ark*(mX+mY)/(mX*mY)
     !
     g_rot(3,3,4) = cos(rho_2)**2*(mX+mY)/(mX*mY)
     g_rot(3,3,6) = cos(rho_2)**2*(mX+mY)/(mX*mY)
     g_rot(3,3,5) = 2.0_ark*cos(rho_2)**2/mX
     !
   end subroutine  MLkinetic_xy2_bisect_EKE_sinrho



  !
  ! Defining kinetic energy function: Radau coordinates 
  ! This is EKE generated using Maple for a bisecting frame with bond-length-angle and 
  ! removed singularity by combining U rho with dG/drho and multiplying muzz by rho^2
  ! and muxz and muyz by rho. The vecinity of zero (singularity) is expanded wrt rho, up to the 7th order 
  !
  subroutine MLkinetic_xy2_Radau_bisect_EKE(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY,rho_2,rho_over_sinrho,rho2_over_sinrho2,pseudo_fac
   real(ark),parameter  :: rho_threshold = 0.01_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_xy2_Radau_bisect_EKE-error: can be used with non-rigid case only')")
       stop 'MLkinetic_xy2_Radau_bisect_EKE can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     !
     rho_2 = rho*0.5_ark
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  1.0_ark/mY
     g_vib(2,2,1) =  1.0_ark/mY
     g_vib(3,3,4) =  1.0_ark/mY
     g_vib(3,3,6) =  1.0_ark/mY
     !
     g_rot(1,1,4) =  0.25_ark/cos(rho_2)**2/mY
     g_rot(1,1,6) =  0.25_ark/cos(rho_2)**2/mY
     !
     g_rot(2,2,4) =  0.25_ark/cos(rho_2)**2/mY
     g_rot(2,2,6) =  0.25_ark/cos(rho_2)**2/mY
     !     
     g_cor(3,2,4) =   .5_ark/mY
     g_cor(3,2,6) =  -.5_ark/mY
     !
     if (rho>rho_threshold) then
        !
        rho_over_sinrho = rho/sin(rho)
        rho2_over_sinrho2 = rho**2/sin(rho)**2
        pseudo_fac  =  1.0_ark/rho - rho - rho/sin(rho)**2
        !
     else
        !
        ! expansion around rho=0
        !
        rho_over_sinrho = 1.0_ark+1.0_ark/6.0_ark*rho**2+7.0_ark/360.0_ark*rho**4+31.0_ark/15120.0_ark*rho**6+&
                         127.0_ark/604800.0_ark*rho**8
        rho2_over_sinrho2 = 1.0_ark+1.0_ark/3.0_ark*rho**2+1.0_ark/15.0_ark*rho**4+2.0_ark/189.0_ark*rho**6+1.0_ark/675.0_ark*rho**8
        pseudo_fac = -4.0_ark/3.0_ark*rho-1.0_ark/15.0_ark*rho**3-2.0_ark/189.0_ark*rho**5-1.0_ark/675.0_ark*rho**7-&
                      2.0_ark/10395.0_ark*rho**9
        !
     endif
     !
     g_rot(1,3,4) =  -0.5_ark*rho_over_sinrho/mY 
     g_rot(1,3,6) =   0.5_ark*rho_over_sinrho/mY 
     !
     g_rot(3,1,4) =  -0.5_ark*rho_over_sinrho/mY 
     g_rot(3,1,6) =   0.5_ark*rho_over_sinrho/mY 
     !
     g_rot(3,3,4) =   cos(rho_2)**2/mY*rho2_over_sinrho2
     g_rot(3,3,6) =   cos(rho_2)**2/mY*rho2_over_sinrho2
     !
     pseudo(4) =  .125_ark/mY*pseudo_fac
     pseudo(6) =  .125_ark/mY*pseudo_fac
     !
   end subroutine  MLkinetic_xy2_Radau_bisect_EKE



!
  ! Defining kinetic energy function 
  ! This is EKE generated using Maple for a bond frame with bond-length-angle and 
  ! removed singularity by combining U rho with dG/drho and multiplying muzz by rho^2
  ! and muxz and muyz by rho. The vecinity of zero (singularity) is expanded wrt rho, up to the 7th order 
  !
  subroutine MLkinetic_xyz_bond_EKE(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY,mZ
   real(ark),parameter  :: rho_threshold = 0.02_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_xyz_bond_EKE-error: can be used with non-rigid case only')")
       stop 'MLkinetic_xyz_bond_EKE can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     mZ = molec%AtomMasses(3)
     !
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  (mX+mY)/mX/mY 
     g_vib(1,2,1) =  -cos(rho)/mX 
     g_vib(1,3,2) =  sin(rho)/mX 
     g_vib(2,1,1) =  -cos(rho)/mX 
     g_vib(2,2,1) =  (mX+mZ)/mX/mZ 
     g_vib(2,3,3) =  sin(rho)/mX 
     g_vib(3,1,2) =  sin(rho)/mX 
     g_vib(3,2,3) =  sin(rho)/mX 
     g_vib(3,3,4) =  (mX+mZ)/mX/mZ 
     g_vib(3,3,5) =  2.0_ark*cos(rho)/mX 
     g_vib(3,3,6) =  (mX+mY)/mX/mY 
     g_rot(1,1,6) =  (mX+mY)/mX/mY 
     g_rot(2,2,6) =  (mX+mY)/mX/mY 
     g_cor(2,2,3) =  -sin(rho)/mX 
     g_cor(3,2,5) =  -cos(rho)/mX 
     g_cor(3,2,6) =  -(mX+mY)/mX/mY 
     !
     if (rho>rho_threshold) then
        !
        g_rot(1,3,5) =  1.0_ark/sin(rho)/mX*rho 
        g_rot(1,3,6) =  rho*cos(rho)*(mX+mY)/sin(rho)/mX/mY 
        g_rot(3,1,5) =  1.0_ark /sin(rho)/mX*rho 
        g_rot(3,1,6) =  rho*cos(rho)*(mX+mY)/sin(rho)/mX/mY 
        g_rot(3,3,4) =  rho**2*(mX+mZ)/sin(rho)**2/mX/mZ 
        g_rot(3,3,5) =  2.0_ark*cos(rho)*rho**2/mX/sin(rho)**2 
        g_rot(3,3,6) =  rho**2*cos(rho)**2*(mX+mY)/mX/mY/sin(rho)**2 
        !
        pseudo(4) =  .125_ark *(mZ+mX-2.0_ark*rho**2*mZ-2.0_ark*rho**2*mX+rho**2*mZ*cos(rho)**2-mZ*cos(rho)**2-&
                     mX*cos(rho)**2+rho**2*mX*cos(rho)**2)/mZ/mX/rho/sin(rho)**2 
        pseudo(5) =  -.25_ark*(-2.0_ark*sin(rho)*rho+2.0_ark*rho*sin(rho)*cos(rho)**2-cos(rho)+cos(rho)**3+&
                     rho**2*cos(rho)**3)/mX/rho/sin(rho)**2 
        pseudo(6) =  .125_ark*(-2.0_ark*rho**2*mY+mY+mX+rho**2*mY*cos(rho)**2-2.0_ark*rho**2*mX+rho**2*mX*cos(rho)**2-&
                     mY*cos(rho)**2-mX*cos(rho)**2)/mY/mX/rho/sin(rho)**2 
        !
     else
        !
        ! expansion around rho=0
        !
        g_rot(1,3,5) = 1._ark/mX+1._ark/mX*rho**2/6._ark+7._ark/360._ark/mX*rho**4+31._ark/15120._ark/mX*rho**6
        g_rot(1,3,6) = (mX+mY)/mX/mY-(mX+mY)/mX/mY*rho**2/3._ark-(mX+mY)/mX/mY*rho**4/45._ark-2._ark/945._ark*(mX+mY)/mX/mY*rho**6
        g_rot(3,1,5) = 1._ark/mX+1._ark/mX*rho**2/6._ark+7._ark/360._ark/mX*rho**4+31._ark/15120._ark/mX*rho**6
        g_rot(3,1,6) = (mX+mY)/mX/mY-(mX+mY)/mX/mY*rho**2/3._ark-(mX+mY)/mX/mY*rho**4/45._ark-2._ark/945._ark*(mX+mY)/mX/mY*rho**6
        g_rot(3,3,4) = (mX+mZ)/mX/mZ+(mX+mZ)/mX/mZ*rho**2/3._ark+(mX+mZ)/mX/mZ*rho**4/15._ark+2._ark/189._ark*(mX+mZ)/mX/mZ*rho**6
        g_rot(3,3,5) = 2._ark/mX-1._ark/mX*rho**2/3._ark-7._ark/60._ark/mX*rho**4-31._ark/1512._ark/mX*rho**6
        g_rot(3,3,6) = (mX+mY)/mX/mY-2._ark/3._ark*(mX+mY)/mX/mY*rho**2+(mX+mY)/mX/mY*rho**4/15._ark+&
                       2._ark/189._ark*(mX+mY)/mX/mY*rho**6
        pseudo(4) = -(mX+mZ)/mX/mZ*rho/6._ark-(mX+mZ)/mX/mZ*rho**3/120._ark-(mX+mZ)/mX/mZ*rho**5/756._ark
        pseudo(5) = 2._ark/3._ark/mX*rho-11._ark/60._ark/mX*rho**3+127._ark/7560._ark/mX*rho**5
        pseudo(6) = -(mX+mY)/mX/mY*rho/6._ark-(mX+mY)/mX/mY*rho**3/120._ark-(mX+mY)/mX/mY*rho**5/756._ark
        !
     endif

     !
   end subroutine  MLkinetic_xyz_bond_EKE





  ! Defining kinetic energy function 
  ! The same as MLkinetic_xyz_bond_EKE but with r2 as the z-bond with the negative direction and r1 in the postive (top-top) x-z corner 
  ! 
  ! This is EKE generated using Maple for a bond frame with bond-length-angle and 
  ! removed singularity by combining U rho with dG/drho and multiplying muzz by rho^2
  ! and muxz and muyz by rho. The vecinity of zero (singularity) is expanded wrt rho, up to the 7th order 
  !
  subroutine MLkinetic_xyz_bond_EKE_r2(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY,mZ
   real(ark),parameter  :: rho_threshold = 0.02_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_xyz_bond_EKE-error: can be used with non-rigid case only')")
       stop 'MLkinetic_xyz_bond_EKE can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     mZ = molec%AtomMasses(3)
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  (mX+mY)/mX/mY 
     g_vib(1,2,1) =  -cos(rho)/mX 
     g_vib(1,3,2) =  sin(rho)/mX 
     g_vib(2,1,1) =  -cos(rho)/mX 
     g_vib(2,2,1) =  (mX+mZ)/mX/mZ 
     g_vib(2,3,3) =  sin(rho)/mX 
     g_vib(3,1,2) =  sin(rho)/mX 
     g_vib(3,2,3) =  sin(rho)/mX 
     g_vib(3,3,4) =  (mX+mZ)/mX/mZ 
     g_vib(3,3,5) =  2.0_ark*cos(rho)/mX 
     g_vib(3,3,6) =  (mX+mY)/mX/mY 
     !
     g_rot(1,1,4) =  (mX+mZ)/mX/mZ 
     g_rot(2,2,4) =  (mX+mZ)/mX/mZ 
     !
     g_cor(1,2,2) =  -sin(rho)/mX 
     g_cor(3,2,4) =  -(mX+mZ)/mX/mZ 
     g_cor(3,2,5) =  -cos(rho)/mX 


     !
     if (rho>rho_threshold) then
        !
        g_rot(1,3,4) =  rho*cos(rho)*(mX+mZ)/sin(rho)/mX/mZ 
        g_rot(1,3,5) =  1.0_ark/sin(rho)/mX*rho 
        g_rot(3,1,4) =  rho*cos(rho)*(mX+mZ)/sin(rho)/mX/mZ 
        g_rot(3,1,5) =  1.0_ark/sin(rho)/mX*rho 
        g_rot(3,3,4) =  rho**2*cos(rho)**2*(mX+mZ)/sin(rho)**2/mX/mZ 
        g_rot(3,3,5) =  2.0_ark*cos(rho)*rho**2/mX/sin(rho)**2 
        g_rot(3,3,6) =  rho**2*(mX+mY)/sin(rho)**2/mX/mY 
        !
        pseudo(4) =  .125_ark*(-2.0_ark*rho**2*mX-mZ*cos(rho)**2+mZ+rho**2*mZ*cos(rho)**2-&
                               2.0_ark*rho**2*mZ-mX*cos(rho)**2+mX+rho**2*mX*cos(rho)**2)/mZ/mX/rho/sin(rho)**2 
        pseudo(5) =  -.25_ark*(2.0_ark*rho*sin(rho)*cos(rho)**2-2.0_ark*sin(rho)*rho-&
                               cos(rho)+cos(rho)**3+rho**2*cos(rho)**3)/mX/rho/sin(rho)**2 
        pseudo(6) =  .125_ark*(-2.0_ark*rho**2*mX-2.0_ark*rho**2*mY-mY*cos(rho)**2+mY+rho**2*mX*cos(rho)**2+&
                               rho**2*mY*cos(rho)**2-mX*cos(rho)**2+mX)/mY/mX/rho/sin(rho)**2 
        !
     else
        !
        ! expansion around rho=0
        !
        g_rot(1,3,4) =  (mX+mZ)/mX/mZ-(mX+mZ)/mX/mZ*rho**2/3._ark-(mX+mZ)/mX/mZ*rho**4/45._ark-2._ark/945._ark*(mX+mZ)/mX/mZ*rho**6
        g_rot(1,3,5) =  1._ark/mX+1._ark/mX*rho**2/6._ark+7._ark/360._ark/mX*rho**4+31._ark/15120._ark/mX*rho**6
        g_rot(3,1,4) =  (mX+mZ)/mX/mZ-(mX+mZ)/mX/mZ*rho**2/3._ark-(mX+mZ)/mX/mZ*rho**4/45._ark-2._ark/945._ark*(mX+mZ)/mX/mZ*rho**6
        g_rot(3,1,5) =  1._ark/mX+1._ark/mX*rho**2/6._ark+7._ark/360._ark/mX*rho**4+31._ark/15120._ark/mX*rho**6
        g_rot(3,3,4) =  (mX+mZ)/mX/mZ-2._ark/3._ark*(mX+mZ)/mX/mZ*rho**2+(mX+mZ)/mX/mZ*rho**4/15._ark+&
                        2._ark/189._ark*(mX+mZ)/mX/mZ*rho**6
        g_rot(3,3,5) =  2._ark/mX-1._ark/mX*rho**2/3._ark-7._ark/60._ark/mX*rho**4-31._ark/1512._ark/mX*rho**6
        g_rot(3,3,6) =  (mX+mY)/mX/mY+(mX+mY)/mX/mY*rho**2/3._ark+(mX+mY)/mX/mY*rho**4/15._ark+2._ark/189._ark*(mX+mY)/mX/mY*rho**6
        !       
        pseudo(4) =  -(mX+mZ)/mX/mZ*rho/6._ark-(mX+mZ)/mX/mZ*rho**3/120._ark-(mX+mZ)/mX/mZ*rho**5/756._ark
        pseudo(5) =   2._ark/3._ark/mX*rho-11._ark/60._ark/mX*rho**3+127._ark/7560._ark/mX*rho**5
        pseudo(6) =  -(mX+mY)/mX/mY*rho/6._ark-(mX+mY)/mX/mY*rho**3/120._ark-(mX+mY)/mX/mY*rho**5/756._ark
        !
     endif
     !
   end subroutine  MLkinetic_xyz_bond_EKE_r2


  ! Defining kinetic energy function: Radau XYZ
  ! This is EKE generated using Maple for a bond frame with bond-length-angle and 
  ! removed singularity by combining U rho with dG/drho and multiplying muzz by rho^2
  ! and muxz and muyz by rho. The vecinity of zero (singularity) is expanded wrt rho, up to the 7th order 
  !
  subroutine MLkinetic_xyz_Radau_EKE(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY,mZ,rho_2
   real(ark),parameter  :: rho_threshold = 0.02_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_xyz_bond_EKE-error: can be used with non-rigid case only')")
       stop 'MLkinetic_xyz_bond_EKE can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     mZ = molec%AtomMasses(3)
     !
     rho_2 = 0.5_ark*rho
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  1.0_ark/mY 
     g_vib(2,2,1) =  1.0_ark/mZ 
     g_vib(3,3,4) =  1.0_ark/mZ 
     g_vib(3,3,6) =  1.0_ark/mY 
     g_rot(1,1,4) =  .250_ark/mZ/cos(rho_2)**2 
     g_rot(1,1,6) =  .250_ark/mY/cos(rho_2)**2 
     g_rot(2,2,4) =  .250_ark/mZ 
     g_rot(2,2,6) =  .250_ark/mY 

     g_cor(3,2,4) =  -.500_ark/mZ 
     g_cor(3,2,6) =  .500_ark/mY 


     if (rho>rho_threshold) then
        !
        g_rot(1,3,4) =  .250_ark*rho/mZ/sin(rho_2)/cos(rho_2) 
        g_rot(1,3,6) =  -.250_ark*rho/mY/sin(rho_2)/cos(rho_2) 
        g_rot(3,1,4) =  .250_ark*rho/mZ/sin(rho_2)/cos(rho_2) 
        g_rot(3,1,6) =  -.250_ark*rho/mY/sin(rho_2)/cos(rho_2) 
        g_rot(3,3,4) =  .125_ark*rho**2/mZ/(1.-cos(rho)) 
        g_rot(3,3,6) =  .125_ark*rho**2/mY/(1.-cos(rho)) 
        !
        pseudo(4) =  .125_ark*(-2._ark*rho**2+1._ark+rho**2*cos(rho)**2-cos(rho)**2)/rho/mZ/sin(rho)**2 
        pseudo(6) =  .125_ark*(-2._ark*rho**2+1._ark+rho**2*cos(rho)**2-cos(rho)**2)/rho/mY/sin(rho)**2 
        !
     else
        !
        ! expansion around rho=0
        !
       g_rot(1,3,4) = 1.0_ark/mZ/2.0_ark+1.0_ark/mZ*rho**2/12.0_ark+7.0_ark/720.0_ark/mZ*rho**4+31.0_ark/30240.0_ark/mZ*rho**6
       g_rot(1,3,6) = -1.0_ark/mY/2.0_ark-1.0_ark/mY*rho**2/12.0_ark-7.0_ark/720.0_ark/mY*rho**4-31.0_ark/30240.0_ark/mY*rho**6
       g_rot(3,1,4) = 1.0_ark/mZ/2.0_ark+1.0_ark/mZ*rho**2/12.0_ark+7.0_ark/720.0_ark/mZ*rho**4+31.0_ark/30240.0_ark/mZ*rho**6
       g_rot(3,1,6) = -1.0_ark/mY/2.0_ark-1.0_ark/mY*rho**2/12.0_ark-7.0_ark/720.0_ark/mY*rho**4-31.0_ark/30240.0_ark/mY*rho**6
       g_rot(3,3,4) = 1.0_ark/mZ+1.0_ark/mZ*rho**2/12.0_ark+1.0_ark/mZ*rho**4/240.0_ark+1.0_ark/mZ*rho**6/6048.0_ark
       g_rot(3,3,6) = 1.0_ark/mY+1.0_ark/mY*rho**2/12.0_ark+1.0_ark/mY*rho**4/240.0_ark+1.0_ark/mY*rho**6/6048.0_ark
       pseudo(4) = -1.0_ark/mZ*rho/6.0_ark-1.0_ark/mZ*rho**3/120.0_ark-1.0_ark/mZ*rho**5/756.0_ark
       pseudo(6) = -1.0_ark/mY*rho/6.0_ark-1.0_ark/mY*rho**3/120.0_ark-1.0_ark/mY*rho**5/756.0_ark
        !
     endif

     !
   end subroutine  MLkinetic_xyz_Radau_EKE




  ! Defining kinetic energy function 
  ! This is EKE generated using Maple for a bisecting frame with bond-length-angle and 
  ! removed singularity by combining U rho with dG/drho and multiplying muzz by sinrho^2
  ! and muxz and muyz by sinrho. 
  !
  subroutine MLkinetic_xyz_EKE_sinrho(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY,mZ
   real(ark),parameter  :: rho_threshold = 0.02_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_xyz_EKE_sinrho-error: can be used with non-rigid case only')")
       stop 'MLkinetic_xyz_EKE_sinrho can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     mZ = molec%AtomMasses(3)
     !
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  (mX+mY)/mX/mY 
     g_vib(1,2,1) =  -cos(rho)/mX 
     g_vib(1,3,2) =  sin(rho)/mX 
     g_vib(2,1,1) =  -cos(rho)/mX 
     g_vib(2,2,1) =  (mX+mZ)/mX/mZ 
     g_vib(2,3,3) =  sin(rho)/mX 
     g_vib(3,1,2) =  sin(rho)/mX 
     g_vib(3,2,3) =  sin(rho)/mX 
     g_vib(3,3,4) =  (mX+mZ)/mX/mZ 
     g_vib(3,3,5) =  2.0_ark*cos(rho)/mX 
     g_vib(3,3,6) =  (mX+mY)/mX/mY 
     !
     g_rot(1,1,6) =  (mX+mY)/mX/mY 
     g_rot(2,2,6) =  (mX+mY)/mX/mY 
     g_cor(2,2,3) =  -sin(rho)/mX 
     !
     g_cor(3,2,5) =  -cos(rho)/mX 
     g_cor(3,2,6) =  -(mX+mY)/mX/mY 
     !
     g_rot(1,3,5) =  1.0_ark/mX
     g_rot(1,3,6) =  cos(rho)*(mX+mY)/mX/mY 
     !
     g_rot(3,1,5) =  1.0_ark/mX
     g_rot(3,1,6) =  cos(rho)*(mX+mY)/mX/mY 
     !
     g_rot(3,3,4) =  (mX+mZ)/mX/mZ 
     g_rot(3,3,5) =  2.0_ark*cos(rho)/mX
     g_rot(3,3,6) =  cos(rho)**2*(mX+mY)/mX/mY
     !
     pseudo(4) =  sin(rho)*cos(rho)/(mX)
     !
   end subroutine  MLkinetic_xyz_EKE_sinrho



  ! Defining kinetic energy function 
  ! This is EKE generated using Maple for a bisect non-symmetric frame with bond-length-angle and 
  ! removed singularity by combining U rho with dG/drho and multiplying muzz by rho^2
  ! and muxz and muyz by rho. The vecinity of zero (singularity) is expanded wrt rho, up to the 7th order 
  ! The geometry was assumed to have atom 1 atom at the centre, atom 2 is in the bottom-left and atom 3 is in the bottom right, 
  ! with negative z2 and z3: 
  ! x1= 0
  ! y1= 0
  ! z1= 0
  ! x2= -r1 cos(alpha/2)
  ! y2=  0
  ! z2= -r1 sin(alpha/2)
  ! x3= -r2 cos(alpha/2)
  ! y3=  0
  ! z3=  r2 sin(alpha/2)
  !
  subroutine MLkinetic_xyz_bisect_EKE(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY1,mY2,rho_2
   real(ark),parameter  :: rho_threshold = 0.02_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_xyz_bisect_EKE-error: can be used with non-rigid case only')")
       stop 'MLkinetic_xyz_bisect_EKE can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY1 = molec%AtomMasses(2)
     mY2 = molec%AtomMasses(3)
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     rho_2 = rho*0.5_ark
     !
     g_vib(1,1,1) =  (mX+mY1)/mX/mY1 
     g_vib(1,2,1) =  -cos(rho)/mX 
     g_vib(1,3,2) =  sin(rho)/mX 
     g_vib(2,1,1) =  -cos(rho)/mX 
     g_vib(2,2,1) =  (mX+mY2)/mY2/mX 
     g_vib(2,3,3) =  sin(rho)/mX 
     g_vib(3,1,2) =  sin(rho)/mX 
     g_vib(3,2,3) =  sin(rho)/mX 
     g_vib(3,3,4) =  (mX+mY2)/mY2/mX 
     g_vib(3,3,5) =  2.0_ark*cos(rho)/mX 
     g_vib(3,3,6) =  (mX+mY1)/mX/mY1 
     !     
     g_rot(1,1,4) =  0.25_ark*(mX+mY2)/cos(rho_2)**2/mX/mY2      
     g_rot(1,1,5) =  -0.5_ark/cos(rho_2)**2/mX 
     g_rot(1,1,6) =  0.25_ark*(mX+mY1)/cos(rho_2)**2/mX/mY1 
     g_rot(2,2,4) =  0.25_ark*(mX+mY2)/mY2/mX 
     g_rot(2,2,5) =  -0.5_ark*(2.0_ark*cos(rho_2)**2-1._ark)/mX 
     g_rot(2,2,6) =  0.25_ark*(mX+mY1)/mX/mY1 
     !
     g_cor(1,2,2) =  -0.5_ark*sin(rho)/mX 
     g_cor(2,2,3) =  0.5_ark*sin(rho)/mX 
     g_cor(3,2,4) =  -0.5_ark*(mX+mY2)/mY2/mX 
     g_cor(3,2,6) =  0.5_ark*(mX+mY1)/mX/mY1 
     !     
     if (rho>rho_threshold) then
       !
       g_rot(1,3,4) =  0.250_ark*rho*(mX+mY2)/cos(rho_2)/sin(rho_2)/mX/mY2 
       g_rot(1,3,6) =  -0.250_ark*rho*(mX+mY1)/cos(rho_2)/sin(rho_2)/mX/mY1 
       g_rot(3,1,4) =  0.250_ark*rho*(mX+mY2)/cos(rho_2)/sin(rho_2)/mX/mY2 
       g_rot(3,1,6) =  -0.250_ark*rho*(mX+mY1)/cos(rho_2)/sin(rho_2)/mX/mY1 
       g_rot(3,3,4) =  0.250_ark*rho**2*(mX+mY2)/sin(rho_2)**2/mX/mY2 
       g_rot(3,3,5) =  0.5_ark*rho**2/mX/sin(rho_2)**2 
       g_rot(3,3,6) =  0.250_ark*rho**2*(mX+mY1)/sin(rho_2)**2/mX/mY1 
       !
       pseudo(4) =  0.125_ark*(rho**2*mY2*cos(rho)**2+rho**2*mX*cos(rho)**2-mX*cos(rho)**2-mY2*cos(rho)**2-&
                    2.0_ark*rho**2*mX-2.0_ark*rho**2*mY2+mY2+mX)/sin(rho)**2/mY2/mX/rho 
       pseudo(5) =  -0.25_ark*(cos(rho)**3+rho**2*cos(rho)**3+2.0_ark*rho*sin(rho)*cos(rho)**2-cos(rho)-&
                    2.0_ark*rho*sin(rho))/sin(rho)**2/mX/rho 
       pseudo(6) =  0.125_ark*(rho**2*mX*cos(rho)**2+rho**2*mY1*cos(rho)**2-mY1*cos(rho)**2-mX*cos(rho)**2+mY1-&
                        2.0_ark*rho**2*mX+mX-2.0_ark*rho**2*mY1)/sin(rho)**2/mX/mY1/rho 
       !
     else
       !
       g_rot(1,3,4) =  (2.0_ark*mX+2.0_ark*mY2)/mX/mY2/4.0_ark+(mX/3.0_ark+mY2/3.0_ark)/mX/mY2*rho**2/4.0_ark+&
                       (7.0_ark/180.0_ark*mX+7.0_ark/180.0_ark*mY2)/mX/mY2*rho**4/4.0_ark+&
                       (31.0_ark/7560.0_ark*mX+31.0_ark/7560.0_ark*mY2)/mX/mY2*rho**6/4.0_ark
       g_rot(1,3,6) =  -(2.0_ark*mX+2.0_ark*mY1)/mX/mY1/4.0_ark-(mX/3.0_ark+mY1/3.0_ark)/mX/mY1*rho**2/4.0_ark-&
                       (7.0_ark/180.0_ark*mX+7.0_ark/180.0_ark*mY1)/mX/mY1*rho**4/4.0_ark-(31.0_ark/7560.0_ark*mX+&
                       31.0_ark/7560.0_ark*mY1)/mX/mY1*rho**6/4.0_ark
       g_rot(3,1,4) =  (2.0_ark*mX+2.0_ark*mY2)/mX/mY2/4.0_ark+(mX/3.0_ark+mY2/3.0_ark)/mX/mY2*rho**2/4.0_ark+&
                       (7.0_ark/180.0_ark*mX+7.0_ark/180.0_ark*mY2)/mX/mY2*rho**4/4.0_ark+(31.0_ark/7560.0_ark*mX+&
                       31.0_ark/7560.0_ark*mY2)/mX/mY2*rho**6/4.0_ark
       g_rot(3,1,6) =  -(2.0_ark*mX+2.0_ark*mY1)/mX/mY1/4.0_ark-(mX/3.0_ark+mY1/3.0_ark)/mX/mY1*rho**2/4.0_ark-&
                       (7.0_ark/180.0_ark*mX+7.0_ark/180.0_ark*mY1)/mX/mY1*rho**4/4.0_ark-(31.0_ark/7560.0_ark*mX+&
                       31.0_ark/7560.0_ark*mY1)/mX/mY1*rho**6/4.0_ark
       g_rot(3,3,4) =  -(-4.0_ark*mX-4.0_ark*mY2)/mX/mY2/4.0_ark-(-mX/3.0_ark-mY2/3.0_ark)/mX/mY2*rho**2/4.0_ark-&
                       (-mX/60.0_ark-mY2/60.0_ark)/mX/mY2*rho**4/4.0_ark
       g_rot(3,3,5) =  2.0_ark/mX+1.0_ark/mX*rho**2/6.0_ark+1.0_ark/mX*rho**4/120.0_ark
       g_rot(3,3,6) =  -(-4.0_ark*mX-4.0_ark*mY1)/mX/mY1/4.0_ark-(-mX/3.0_ark-mY1/3.0_ark)/mX/mY1*rho**2/4.0_ark-&
                       (-mX/60.0_ark-mY1/60.0_ark)/mX/mY1*rho**4/4.0_ark
       !
       pseudo(4) =  -(4.0_ark/3.0_ark*mY2+4.0_ark/3.0_ark*mX)/mY2/mX*rho/8.0_ark-(mY2/15.0_ark+mX/15.0_ark)/mY2/mX*rho**3/8.0_ark
       pseudo(5) =    2.0_ark/3.0_ark/mX*rho-11.0_ark/60.0_ark/mX*rho**3
       pseudo(6) =  -(4.0_ark/3.0_ark*mX+4.0_ark/3.0_ark*mY1)/mX/mY1*rho/8.0_ark-(mX/15.0_ark+mY1/15.0_ark)/mX/mY1*rho**3/8.0_ark
       !     
     endif
     !
   end subroutine  MLkinetic_xyz_bisect_EKE
   
   

  !
  !
  ! Defining kinetic energy function: sparse representation, rigid congiguration for XY2 (away from singularity)
  ! This is EKE generated using Maple for a bisecting frame with bond-length-angle. 
  !
  subroutine MLkinetic_compact_xy2_bisect_EKE_rigid(nmodes,rho,ntermmax,ng_vib,ng_rot,ng_cor,npseudo,&
                                                    g_vib,g_rot,g_cor,pseudo,ig_vib,ig_rot,ig_cor,ipseudo)
   !
   use accuracy
   !
   integer(ik),intent(in) ::  nmodes
   real(ark),intent(in)   ::  rho
   integer(ik),intent(in) ::  ntermmax
   integer(ik),intent(inout) ::  ng_vib(nmodes,nmodes),ng_rot(3,3),ng_cor(nmodes,3),npseudo
   real(ark),intent(out)     ::  g_vib(nmodes,nmodes,ntermmax),g_rot(3,3,ntermmax),g_cor(nmodes,3,ntermmax),pseudo(ntermmax)
   integer(ik),intent(out)   ::  ig_vib(nmodes,nmodes,ntermmax,nmodes),ig_rot(3,3,ntermmax,nmodes),&
                                 ig_cor(nmodes,3,ntermmax,nmodes),ipseudo(ntermmax,nmodes)
   !
   !type(FLpolynomT),pointer ::  g_vib(:,:) 
   !
   real(ark)            :: mX,mY
   integer(ik) :: info,Nterms,NMax
     !
     if (manifold==1) then
       write(out,"('MLkinetic_compact_xy2_bisect_EKE_rigid-error: can be used with rigid case only')")
       stop 'MLkinetic_compact_xy2_bisect_EKE_rigid can be used only with npoints=0'
     endif
     !
     NMax = 13
     !
     if (size(pseudo,dim=1)<NMax) then
       write(out,"('MLkinetic_compact_xy2_bisect_EKE_rigid-error: The NKinOrder is too small, increas to',i4)") NMax
       stop 'MLkinetic_compact_xy2_bisect_EKE_rigid The NKinOrder is too small'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     !
     Ng_vib = 0
     Ng_rot = 0
     Ng_cor = 0
     !
     g_vib = 0
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     Ng_vib(1,1) = 1     
     Ng_vib(1,2) = 2
     Ng_vib(1,3) = 1
     Ng_vib(2,1) = 2
     Ng_vib(2,2) = 1
     Ng_vib(2,3) = 1
     Ng_vib(3,1) = 1
     Ng_vib(3,2) = 1
     Ng_vib(3,3) = 4
     Ng_rot(1,1) = 3
     Ng_rot(1,3) = 2
     Ng_rot(2,2) = 4
     Ng_rot(3,1) = 2
     Ng_rot(3,3) = 3
     Ng_cor(1,2) = 1
     Ng_cor(2,2) = 1
     Ng_cor(3,2) = 2
     Npseudo = 13
     !
     ig_vib(1,1,1,:) = (/0,0,0/)
     ig_vib(1,2,1,:) = (/0,0,0/)
     ig_vib(1,2,2,:) = (/0,0,1/)
     ig_vib(1,3,1,:) = (/0,2,4/)
     ig_vib(2,1,1,:) = (/0,0,0/)
     ig_vib(2,1,2,:) = (/0,0,1/)
     ig_vib(2,2,1,:) = (/0,0,0/)
     ig_vib(2,3,1,:) = (/2,0,4/)
     ig_vib(3,1,1,:) = (/0,2,4/)
     ig_vib(3,2,1,:) = (/2,0,4/)
     ig_vib(3,3,1,:) = (/0,1,0/)
     ig_vib(3,3,2,:) = (/1,0,0/)
     ig_vib(3,3,3,:) = (/2,2,0/)
     ig_vib(3,3,4,:) = (/2,2,1/)
     ig_rot(1,1,1,:) = (/0,1,3/)
     ig_rot(1,1,2,:) = (/1,0,3/)
     ig_rot(1,1,3,:) = (/2,2,3/)
     ig_rot(1,3,1,:) = (/0,1,5/)
     ig_rot(1,3,2,:) = (/1,0,5/)
     ig_rot(2,2,1,:) = (/0,1,0/)
     ig_rot(2,2,2,:) = (/1,0,0/)
     ig_rot(2,2,3,:) = (/2,2,0/)
     ig_rot(2,2,4,:) = (/2,2,1/)
     ig_rot(3,1,1,:) = (/0,1,5/)
     ig_rot(3,1,2,:) = (/1,0,5/)
     ig_rot(3,3,1,:) = (/0,1,2/)
     ig_rot(3,3,2,:) = (/1,0,2/)
     ig_rot(3,3,3,:) = (/2,2,2/)
     ig_cor(1,2,1,:) = (/0,2,4/)
     ig_cor(2,2,1,:) = (/2,0,4/)
     ig_cor(3,2,1,:) = (/0,1,0/)
     ig_cor(3,2,2,:) = (/1,0,0/)
     ipseudo(1,:) = (/0,1,0/)
     ipseudo(2,:) = (/1,0,0/)
     ipseudo(3,:) = (/2,2,0/)
     ipseudo(4,:) = (/2,2,1/)
     ipseudo(5,:) = (/0,1,2/)
     ipseudo(6,:) = (/1,0,2/)
     ipseudo(7,:) = (/0,1,3/)
     ipseudo(8,:) = (/1,0,3/)
     ipseudo(9,:) = (/2,2,2/)
     ipseudo(10,:) = (/0,1,6/)
     ipseudo(11,:) = (/1,0,6/)
     ipseudo(12,:) = (/2,2,3/)
     ipseudo(13,:) = (/2,2,6/)
     !
     g_vib(1,1,1) =  (mX+mY)/mX/mY
     g_vib(1,2,1) =  -1.0_ark/mX
     g_vib(1,2,2) =  2.0_ark/mX
     g_vib(1,3,1) =  -1.0_ark/mX
     g_vib(2,1,1) =  -1.0_ark/mX
     g_vib(2,1,2) =  2.0_ark/mX
     g_vib(2,2,1) =  (mX+mY)/mX/mY
     g_vib(2,3,1) =  -1.0_ark/mX
     g_vib(3,1,1) =  -1.0_ark/mX
     g_vib(3,2,1) =  -1.0_ark/mX
     g_vib(3,3,1) =  (mX+mY)/mX/mY
     g_vib(3,3,2) =  (mX+mY)/mX/mY
     g_vib(3,3,3) =  2.0_ark/mX
     g_vib(3,3,4) =  -4.0_ark/mX
     g_rot(1,1,1) =  .250_ark*(mX+mY)/mX/mY
     g_rot(1,1,2) =  .250_ark*(mX+mY)/mX/mY
     g_rot(1,1,3) =  -.5_ark/mX
     g_rot(1,3,1) =  .5_ark*(mX+mY)/mX/mY
     g_rot(1,3,2) =  -.5_ark*(mX+mY)/mX/mY
     g_rot(2,2,1) =  .25_ark*(mX+mY)/mX/mY
     g_rot(2,2,2) =  .25_ark*(mX+mY)/mX/mY
     g_rot(2,2,3) =  -.5_ark/mX
     g_rot(2,2,4) =  1.0_ark/mX
     g_rot(3,1,1) =  .5_ark*(mX+mY)/mX/mY
     g_rot(3,1,2) =  -.5_ark*(mX+mY)/mX/mY
     g_rot(3,3,1) =  .25_ark*(mX+mY)/mX/mY
     g_rot(3,3,2) =  .25_ark*(mX+mY)/mX/mY
     g_rot(3,3,3) =  .5_ark/mX
     g_cor(1,2,1) =  -.5_ark/mX
     g_cor(2,2,1) =  .5_ark/mX
     g_cor(3,2,1) =  .5_ark*(mX+mY)/mX/mY
     g_cor(3,2,2) =  -.5_ark*(mX+mY)/mX/mY
     pseudo(1) =  -.03125_ark*(6.0_ark*mY+6.0_ark*mX)/mX/mY
     pseudo(2) =  -.03125_ark*(6.0_ark*mY+6.0_ark*mX)/mX/mY
     pseudo(3) =  .375_ark/mX
     pseudo(4) =  -.5_ark/mX
     pseudo(5) =  -.03125_ark*(mX+mY)/mX/mY
     pseudo(6) =  -.03125_ark*(mX+mY)/mX/mY
     pseudo(7) =  .03125_ark*(mX+mY)/mX/mY
     pseudo(8) =  .03125_ark*(mX+mY)/mX/mY
     pseudo(9) =  -.0625_ark/mX
     pseudo(10) =  -.0625_ark*(mX+mY)/mX/mY
     pseudo(11) =  -.0625_ark*(mX+mY)/mX/mY
     pseudo(12) =  -.0625_ark/mX
     pseudo(13) =  .125_ark/mX
     !
   end subroutine  MLkinetic_compact_xy2_bisect_EKE_rigid



  !
  !
  ! Defining kinetic energy function: sparse representation, rigid congiguration for XYZ (away from singularity)
  ! with z as the second bond and alpha as the angle cordinate 
  ! The KEO was generated using Maple using the 2nd bond-length-angle frame, see XYZ_KEO_analytic_v1.mws
  !
  subroutine MLkinetic_compact_xyz_alpha_bond2_EKE_rigid(nmodes,rho,ntermmax,ng_vib,ng_rot,ng_cor,npseudo,&
                                                    g_vib,g_rot,g_cor,pseudo,ig_vib,ig_rot,ig_cor,ipseudo)
   !
   use accuracy
   !
   integer(ik),intent(in) ::  nmodes
   real(ark),intent(in)   ::  rho
   integer(ik),intent(in) ::  ntermmax
   integer(ik),intent(inout) ::  ng_vib(nmodes,nmodes),ng_rot(3,3),ng_cor(nmodes,3),npseudo
   real(ark),intent(out)     ::  g_vib(nmodes,nmodes,ntermmax),g_rot(3,3,ntermmax),g_cor(nmodes,3,ntermmax),pseudo(ntermmax)
   integer(ik),intent(out)   ::  ig_vib(nmodes,nmodes,ntermmax,nmodes),ig_rot(3,3,ntermmax,nmodes),&
                                 ig_cor(nmodes,3,ntermmax,nmodes),ipseudo(ntermmax,nmodes)
   !
   !type(FLpolynomT),pointer ::  g_vib(:,:) 
   !
   real(ark)            :: mX,mY,mZ
   integer(ik) :: info,Nterms,NMax

   integer(ik),parameter :: nlines = 18
   character(len=wl) :: constructor(nlines) = (/&
      'Mode 1 2',&
      '1 1 -1 I 1 1',&
      '2 1 -2 I 1 1',&
      'Mode 2 2',&
      '1 1 -1 I 1 1',&
      '2 1 -2 I 1 1',&
      'Mode 3 10',&
      ' 1 1 1 cos 1.0 1',&
      ' 2 1 1 sin 1.0 1',&
      ' 3 1 2 cos 1.0 1',&
      ' 4 1 2 sin 1.0 1',&
      ' 5 1 1 csc 1.0 1',&
      ' 6 1 2 csc 1.0 1',&
      ' 7 2 1 cos 1.0 1 2 csc 1.0 1',& 
      ' 8 2 1 cos 1.0 1 1 csc 1.0 1',& 
      ' 9 2 2 cos 1.0 1 2 csc 1.0 1',& 
      '10 2 3 cos 1.0 1 2 csc 1.0 1',& 
      'end'/)

     
     !
     if (manifold==1) then
       write(out,"('MLkinetic_compact_xyz_alpha_bond2_EKE_rigid-error: can be used with rigid case only')")
       stop 'MLkinetic_compact_xyz_alpha_bond2_EKE_rigid can be used only with npoints=0'
     endif
     !
     NMax = 8
     !
     if (size(pseudo,dim=1)<NMax) then
       write(out,"('MLkinetic_compact_xyz_alpha_bond2_EKE_rigid-error: The NKinOrder is too small, increas to',i4)") NMax
       stop 'MLkinetic_compact_xyz_alpha_bond2_EKE_rigid The NKinOrder is too small'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(2)
     mz = molec%AtomMasses(3)
     !
     Ng_vib = 0
     Ng_rot = 0
     Ng_cor = 0
     !
     g_vib = 0
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     Ng_vib(1,1) = 2
     Ng_vib(1,2) = 1
     Ng_vib(1,3) = 1
     Ng_vib(2,1) = 1
     Ng_vib(2,2) = 1
     Ng_vib(2,3) = 1
     Ng_vib(3,1) = 1
     Ng_vib(3,2) = 1
     Ng_vib(3,3) = 4
     Ng_rot(1,1) = 1
     Ng_rot(1,3) = 2
     Ng_rot(2,2) = 1
     Ng_rot(3,1) = 2
     Ng_rot(3,3) = 3
     Ng_cor(1,2) = 1
     Ng_cor(3,2) = 2
     Npseudo = 8

     !
     ig_vib(1,1,1,:) = (/0,0,3/)
     ig_vib(1,1,2,:) = (/0,0,4/)
     ig_vib(1,2,1,:) = (/0,0,1/)
     ig_vib(1,3,1,:) = (/0,2,2/)
     ig_vib(2,1,1,:) = (/0,0,1/)
     ig_vib(2,2,1,:) = (/0,0,0/)
     ig_vib(2,3,1,:) = (/2,0,2/)
     ig_vib(3,1,1,:) = (/0,2,2/)
     ig_vib(3,2,1,:) = (/2,0,2/)
     ig_vib(3,3,1,:) = (/0,1,0/)
     ig_vib(3,3,2,:) = (/1,0,3/)
     ig_vib(3,3,3,:) = (/1,0,4/)
     ig_vib(3,3,4,:) = (/2,2,1/)
     ig_rot(1,1,1,:) = (/0,1,0/)
     ig_rot(1,3,1,:) = (/0,1,8/)
     ig_rot(1,3,2,:) = (/2,2,5/)
     ig_rot(2,2,1,:) = (/0,1,0/)
     ig_rot(3,1,1,:) = (/0,1,8/)
     ig_rot(3,1,2,:) = (/2,2,5/)
     ig_rot(3,3,1,:) = (/1,0,6/)
     ig_rot(3,3,2,:) = (/0,1,9/)
     ig_rot(3,3,3,:) = (/2,2,7/)
     ig_cor(1,2,1,:) = (/0,2,2/)
     ig_cor(3,2,1,:) = (/0,1,0/)
     ig_cor(3,2,2,:) = (/2,2,1/)
     !
     ipseudo(1,:) = (/0,1,0/)
     ipseudo(2,:) = (/1,0,0/)
     ipseudo(3,:) = (/1,0,3/)
     ipseudo(4,:) = (/1,0,4/)
     ipseudo(5,:) = (/1,0,6/)
     ipseudo(6,:) = (/0,1,9/)
     ipseudo(7,:) = (/1,0,9/)
     ipseudo(8,:) = (/2,2,10/)
     !
     g_vib(1,1,1) =  (mX+mY)/mY/mX
     g_vib(1,1,2) =  (mX+mY)/mY/mX
     g_vib(1,2,1) =  1._ark/mX
     g_vib(1,3,1) =  -1._ark/mX
     g_vib(2,1,1) =  1._ark/mX
     g_vib(2,2,1) =  (mX+mZ)/mX/mZ
     g_vib(2,3,1) =  -1._ark/mX
     g_vib(3,1,1) =  -1._ark/mX
     g_vib(3,2,1) =  -1._ark/mX
     g_vib(3,3,1) =  (mX+mZ)/mX/mZ
     g_vib(3,3,2) =  (mX+mY)/mX/mY
     g_vib(3,3,3) =  (mX+mY)/mX/mY
     g_vib(3,3,4) =  -2._ark/mX
     g_rot(1,1,1) =  (mX+mZ)/mX/mZ
     g_rot(1,3,1) =  -1.*(mX+mZ)/mX/mZ
     g_rot(1,3,2) =  1._ark/mX
     g_rot(2,2,1) =  (mX+mZ)/mX/mZ
     g_rot(3,1,1) =  -1._ark*(mX+mZ)/mX/mZ
     g_rot(3,1,2) =  1._ark/mX
     g_rot(3,3,1) =  (mX+mY)/mY/mX
     g_rot(3,3,2) =  (mX+mZ)/mZ/mX
     g_rot(3,3,3) =  -2._ark/mX
     g_cor(1,2,1) =  -1._ark/mX
     g_cor(3,2,1) =  (mX+mZ)/mX/mZ
     g_cor(3,2,2) =  -1._ark/mX
     pseudo(1) =  -.125_ark*(2.*mZ+2.*mX)/mX/mZ
     pseudo(2) =  .125_ark*(-2.*mX-2.*mY)/mY/mX
     pseudo(3) =  -.125_ark*(mX+mY)/mY/mX
     pseudo(4) =  -.125_ark*(mX+mY)/mY/mX
     pseudo(5) =  .125_ark*(mX+mY)/mY/mX
     pseudo(6) =  -.125_ark*(mX+mZ)/mX/mZ
     pseudo(7) =  -.250_ark*(mX+mY)/mY/mX
     pseudo(8) =  .250_ark/mX

     !
     call read_basic_function_constructor(nlines,constructor)
     !
   end subroutine  MLkinetic_compact_xyz_alpha_bond2_EKE_rigid
   
   

end module kin_xy2
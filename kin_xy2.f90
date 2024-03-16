!
!  This unit defines all kinetic routines for a triatomic molecule of the XY2 type
!
module kin_xy2
  use accuracy
  use moltype

  implicit none

  public MLkinetic_xy2_bisect_EKE,MLkinetic_xyz_bisect_EKE,MLkinetic_xy2_bisect_EKE_sinrho,&
         MLkinetic_xy2_Radau_bisect_EKE,MLkinetic_xyz_EKE_sinrho,MLkinetic_xyz_bond_EKE,MLkinetic_xyz_bond_EKE_r2
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
        g_rot(3,3,5) = 2._ark/mX-1._ark/mX*rho**2/3._ark-7._ark/60._ark/mX*rho**4-31._ark/1512.D0/mX*rho**6
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
  ! The same as MLkinetic_xyz_bond_EKE but with r1 and r2 swapped in order to allow for r2 to be the z-bond
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
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  (mX+mZ)/mX/mZ 
     g_vib(1,2,1) =  -cos(rho)/mX 
     g_vib(1,3,2) =  sin(rho)/mX 
     g_vib(2,1,1) =  -cos(rho)/mX 
     g_vib(2,2,1) =  (mX+mY)/mX/mY 
     g_vib(2,3,3) =  sin(rho)/mX 
     g_vib(3,1,2) =  sin(rho)/mX 
     g_vib(3,2,3) =  sin(rho)/mX 
     g_vib(3,3,4) =  (mX+mY)/mX/mY 
     g_vib(3,3,5) =  2.0_ark*cos(rho)/mX 
     g_vib(3,3,6) =  (mX+mZ)/mX/mZ 
     !
     g_rot(1,1,4) =  (mX+mY)/mX/mY 
     g_rot(2,2,4) =  (mX+mY)/mX/mY 
     !
     g_cor(1,2,2) =  -sin(rho)/mX 
     g_cor(3,2,4) =  -(mX+mY)/mX/mY 
     g_cor(3,2,5) =  -cos(rho)/mX 

     !
     if (rho>rho_threshold) then
        !
        g_rot(1,3,4) =  rho*cos(rho)*(mX+mY)/sin(rho)/mX/mY 
        g_rot(1,3,5) =  1.0_ark/sin(rho)/mX*rho 
        g_rot(3,1,4) =  rho*cos(rho)*(mX+mY)/sin(rho)/mX/mY 
        g_rot(3,1,5) =  1.0_ark/sin(rho)/mX*rho 
        g_rot(3,3,4) =  rho**2*cos(rho)**2*(mX+mY)/mX/mY/sin(rho)**2
        g_rot(3,3,5) =  2.0_ark*cos(rho)*rho**2/mX/sin(rho)**2 
        g_rot(3,3,6) =  rho**2*(mX+mZ)/sin(rho)**2/mX/mZ 
        !
        pseudo(4) =  .125_ark*(-mY*cos(rho)**2+mY+mX+rho**2*mX*cos(rho)**2-2.0_ark*rho**2*mX+rho**2*mY*cos(rho)**2-&
                               2.0_ark*rho**2*mY-cos(rho)**2*mX)/mY/mX/rho/sin(rho)**2 
        pseudo(5) =  -.25_ark*(rho**2*cos(rho)**3-cos(rho)+cos(rho)**3+2.0_ark*rho*sin(rho)*cos(rho)**2-&
                               2.0_ark*sin(rho)*rho)/mX/rho/sin(rho)**2 
        pseudo(6) =  .125_ark*(-2.0_ark*rho**2*mX-2.0_ark*rho**2*mZ+mX+mZ+rho**2*mX*cos(rho)**2+rho**2*mZ*cos(rho)**2-&
                               mZ*cos(rho)**2-cos(rho)**2*mX)/mZ/mX/rho/sin(rho)**2 
        !
     else
        !
        ! expansion around rho=0
        !
        g_rot(1,3,4) =  (mX+mY)/mX/mY-(mX+mY)/mX/mY*rho**2/3.0_ark-(mX+mY)/mX/mY*rho**4/45.0_ark-2.0_ark/945.0_ark*(mX+mY)/mX/mY*rho**6
        g_rot(1,3,5) =  1.0_ark/mX+1.0_ark/mX*rho**2/6.0_ark+7.0_ark/360.0_ark/mX*rho**4+31.0_ark/15120.0_ark/mX*rho**6
        g_rot(3,1,4) =  (mX+mY)/mX/mY-(mX+mY)/mX/mY*rho**2/3.0_ark-(mX+mY)/mX/mY*rho**4/45.0_ark-2.0_ark/945.0_ark*(mX+mY)/mX/mY*rho**6
        g_rot(3,1,5) =  1.0_ark/mX+1.0_ark/mX*rho**2/6.0_ark+7.0_ark/360.0_ark/mX*rho**4+31.0_ark/15120.0_ark/mX*rho**6
        g_rot(3,3,4) =  (mX+mY)/mX/mY-2.0_ark/3.0_ark*(mX+mY)/mX/mY*rho**2+(mX+mY)/mX/mY*rho**4/15.0_ark+2.0_ark/189.0_ark*(mX+mY)/mX/mY*rho**6
        g_rot(3,3,5) =  2.0_ark/mX-1.0_ark/mX*rho**2/3.0_ark-7.0_ark/60.0_ark/mX*rho**4-31.0_ark/1512.0_ark/mX*rho**6
        g_rot(3,3,6) =  (mX+mZ)/mX/mZ+(mX+mZ)/mX/mZ*rho**2/3.0_ark+(mX+mZ)/mX/mZ*rho**4/15.0_ark+2.0_ark/189.0_ark*(mX+mZ)/mX/mZ*rho**6
      
        pseudo(4) = -(mX+mY)/mY/mX*rho/6.0_ark-(mX+mY)/mY/mX*rho**3/120.0_ark-(mX+mY)/mY/mX*rho**5/756.0_ark
        pseudo(5) =  2.0_ark/3.0_ark/mX*rho-11.0_ark/60.0_ark/mX*rho**3+127.0_ark/7560.0_ark/mX*rho**5
        pseudo(6) =  -(mX+mZ)/mX/mZ*rho/6.0_ark-(mX+mZ)/mX/mZ*rho**3/120.0_ark-(mX+mZ)/mX/mZ*rho**5/756.0_ark
       !
     endif

     !
   end subroutine  MLkinetic_xyz_bond_EKE_r2





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

end module kin_xy2
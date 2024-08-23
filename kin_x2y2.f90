!
!  This unit defines all kinetic routines for a triatomic molecule of the X2Y2 type
!
module kin_x2y2
  use accuracy
  use moltype
  use timer

  implicit none

  public MLkinetic_x2y2_bisect_EKE_sinrho,MLkinetic_compact_x2y2_bisect_EKE_sinrho_rigid
  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains


  !
  ! Defining kinetic energy function 
  ! This is EKE generated using Maple for a bisecting frame with bond-length-angle and 
  ! removed singularity by combining U sin(theta1)sin(theta2) with dG/d theta_i and multiplying muzz and G66 by sin(theta1)^2 sin(theta2)^2 
  ! and muxz and muyz by U sin(theta1)sin(theta2. 
  !
  subroutine MLkinetic_x2y2_bisect_EKE_sinrho(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   use accuracy
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
   !
   real(ark)            :: mX,mY,rho_2
   real(ark),parameter  :: rho_threshold = 0.01_rk
     !
     if (manifold/=1) then
       write(out,"('MLkinetic_x2y2_bisect_EKE_sinrho-error: can be used with non-rigid case only')")
       stop 'MLkinetic_x2y2_bisect_EKE_sinrho can be used only with npoints>0'
     endif
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(3)
     !
     rho_2 = rho*0.5_ark
     !
     g_vib = 0 
     g_rot = 0
     g_cor = 0
     pseudo = 0
     !
     g_vib(1,1,1) =  2._ark/mX
     g_vib(1,2,9) =  -1./mX
     g_vib(1,3,7) =  -1._ark/mX
     g_vib(1,4,39) =  1._ark/mX
     g_vib(1,5,29) =  1._ark/mX
     g_vib(2,1,9) =  -1._ark/mX
     g_vib(2,2,1) =  (mY+mX)/mY/mX
     g_vib(2,4,53) =  1._ark/mX
     g_vib(2,5,53) =  -1._ark*cos(rho)/mX
     g_vib(2,6,112) =  2._ark*cos(rho_2)*sin(rho_2)/mX
     g_vib(2,6,670) =  -2._ark*cos(rho_2)*sin(rho_2)/mX
     g_vib(3,1,7) =  -1._ark/mX
     g_vib(3,3,1) =  (10._ark*mX**2*mY+7._ark*mY**2*mX+mY**3+4._ark*mX**3)/(2._ark*mX+mY)**2/mX/mY
     g_vib(3,4,52) =  -1._ark*cos(rho)/mX
     g_vib(3,5,52) =  1/mX
     g_vib(3,6,114) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(3,6,669) =  -2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(4,1,39) =  1._ark/mX
     g_vib(4,2,53) =  1._ark/mX
     g_vib(4,3,52) =  -1._ark*cos(rho)/mX
     g_vib(4,4,5) =  (mY+mX)/mX/mY
     g_vib(4,4,6) =  2./mX
     g_vib(4,4,420) =  2./mX
     g_vib(4,5,6) =  -2.*cos(rho)/mX
     g_vib(4,5,402) =  -1.*cos(rho)/mX
     g_vib(4,5,420) =  -1.*cos(rho)/mX
     g_vib(4,6,93) =  4.*cos(rho_2)*sin(rho_2)/mX
     g_vib(4,6,226) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(4,6,2594) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(5,1,29) =  1/mX
     g_vib(5,2,53) =  -1.*cos(rho)/mX
     g_vib(5,3,52) =  1/mX
     g_vib(5,4,6) =  -2.*cos(rho)/mX
     g_vib(5,4,402) =  -1.*cos(rho)/mX
     g_vib(5,4,420) =  -1.*cos(rho)/mX
     g_vib(5,5,4) =  (4.*mX**3+7.*mY**2*mX+10.*mX**2*mY+mY**3)/(2.*mX+mY)**2/mX/mY
     g_vib(5,5,6) =  2./mX
     g_vib(5,5,402) =  2./mX
     g_vib(5,6,94) =  4.*cos(rho_2)*sin(rho_2)/mX
     g_vib(5,6,234) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(5,6,2528) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,2,112) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,2,670) =  -2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,3,114) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,3,669) =  -2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,4,93) =  4.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,4,226) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,4,2594) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,5,94) =  4.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,5,234) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,5,2528) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(6,6,4) =  (10.*mX**2*mY+7.*mY**2*mX+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_vib(6,6,5) =  (mX+mY)/mY/mX
     g_vib(6,6,65) =  -1.*(10.*mX**2*mY+7.*mY**2*mX+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_vib(6,6,72) =  -1.*(mX+mY)/mX/mY
     g_vib(6,6,92) =  2./mX
     g_vib(6,6,95) =  2./mX
     g_vib(6,6,402) =  2./mX
     g_vib(6,6,420) =  2./mX
     g_vib(6,6,586) =  -4./mX
     g_vib(6,6,1512) =  4.*cos(rho)/mX
     g_vib(6,6,1691) =  -2./mX
     g_vib(6,6,1692) =  2.*cos(rho)/mX
     g_vib(6,6,1738) =  2.*cos(rho)/mX
     g_vib(6,6,1739) =  -2./mX

     g_rot(1,1,6) =  2./mX
     g_rot(1,3,93) =  -1.*cos(rho_2)/mX
     g_rot(1,3,94) =  cos(rho_2)/mX
     g_rot(1,3,226) =  -.5000000000*cos(rho_2)/mX
     g_rot(1,3,234) =  .5000000000*cos(rho_2)/mX
     g_rot(2,2,6) =  2./mX
     g_rot(2,3,93) =  sin(rho_2)/mX
     g_rot(2,3,94) =  sin(rho_2)/mX
     g_rot(2,3,226) =  .5000000000*sin(rho_2)/mX
     g_rot(2,3,234) =  .5000000000*sin(rho_2)/mX
     g_rot(3,1,93) =  -1.*cos(rho_2)/mX
     g_rot(3,1,94) =  cos(rho_2)/mX
     g_rot(3,1,226) =  -.5000000000*cos(rho_2)/mX
     g_rot(3,1,234) =  .5000000000*cos(rho_2)/mX
     g_rot(3,2,93) =  sin(rho_2)/mX
     g_rot(3,2,94) =  sin(rho_2)/mX
     g_rot(3,2,226) =  .5000000000*sin(rho_2)/mX
     g_rot(3,2,234) =  .5000000000*sin(rho_2)/mX
     g_rot(3,3,4) =  .2500000000*(7.*mY**2*mX+10.*mX**2*mY+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_rot(3,3,5) =  .2500000000*(mX+mY)/mY/mX
     g_rot(3,3,65) =  -.2500000000*(7.*mY**2*mX+10.*mX**2*mY+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_rot(3,3,72) =  -.2500000000*(mX+mY)/mY/mX
     g_rot(3,3,92) =  .5000000000/mX
     g_rot(3,3,95) =  .5000000000/mX
     g_rot(3,3,402) =  .5000000000/mX
     g_rot(3,3,420) =  .5000000000/mX
     g_rot(3,3,586) =  -1./mX
     g_rot(3,3,1512) =  -1.*cos(rho)/mX
     g_rot(3,3,1691) =  -.5000000000/mX
     g_rot(3,3,1692) =  -.5000000000*cos(rho)/mX
     g_rot(3,3,1738) =  -.5000000000*cos(rho)/mX
     g_rot(3,3,1739) =  -.5000000000/mX

     g_cor(2,1,53) =  sin(rho_2)/mX
     g_cor(2,2,53) =  -1._ark*cos(rho_2)/mX
     g_cor(2,3,112) =  -1._ark*cos(rho_2)*sin(rho_2)/mX
     g_cor(2,3,670) =  cos(rho_2)*sin(rho_2)/mX
     g_cor(3,1,52) =  sin(rho_2)/mX
     g_cor(3,2,52) =  cos(rho_2)/mX
     g_cor(3,3,114) =  cos(rho_2)*sin(rho_2)/mX
     g_cor(3,3,669) =  -1._ark*cos(rho_2)*sin(rho_2)/mX
     g_cor(4,1,6) =  2._ark*sin(rho_2)/mX
     g_cor(4,1,420) =  sin(rho_2)/mX
     g_cor(4,2,6) =  -2._ark*cos(rho_2)/mX
     g_cor(4,2,420) =  -1._ark*cos(rho_2)/mX
     g_cor(4,3,93) =  -2._ark*cos(rho_2)*sin(rho_2)/mX
     g_cor(4,3,226) =  -1._ark*cos(rho_2)*sin(rho_2)/mX
     g_cor(4,3,2594) =  -1._ark*cos(rho_2)*sin(rho_2)/mX
     g_cor(5,1,6) =  2._ark*sin(rho_2)/mX
     g_cor(5,1,402) =  sin(rho_2)/mX
     g_cor(5,2,6) =  2._ark*cos(rho_2)/mX
     g_cor(5,2,402) =  cos(rho_2)/mX
     g_cor(5,3,94) =  2._ark*cos(rho_2)*sin(rho_2)/mX
     g_cor(5,3,234) =  cos(rho_2)*sin(rho_2)/mX
     g_cor(5,3,2528) =  cos(rho_2)*sin(rho_2)/mX
     g_cor(6,1,93) =  2._ark*cos(rho_2)/mX
     g_cor(6,1,94) =  2._ark*cos(rho_2)/mX
     g_cor(6,1,226) =  cos(rho_2)/mX
     g_cor(6,1,234) =  cos(rho_2)/mX
     g_cor(6,2,93) =  -2._ark*sin(rho_2)/mX
     g_cor(6,2,94) =  2._ark*sin(rho_2)/mX
     g_cor(6,2,226) =  -1._ark*sin(rho_2)/mX
     g_cor(6,2,234) =  sin(rho_2)/mX
     g_cor(6,3,4) =  -.5_ark*(10.*mX**2*mY+7.*mY**2*mX+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_cor(6,3,5) =  .5_ark*(mX+mY)/mX/mY
     g_cor(6,3,65) =  .5_ark*(10.*mX**2*mY+7.*mY**2*mX+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_cor(6,3,72) =  -.5_ark*(mX+mY)/mX/mY
     g_cor(6,3,92) =  -1._ark/mX
     g_cor(6,3,95) =  1._ark/mX
     g_cor(6,3,402) =  -1._ark/mX
     g_cor(6,3,420) =  1._ark/mX
     g_cor(6,3,1691) =  1._ark/mX
     g_cor(6,3,1739) =  -1._ark/mX


     pseudo(185) =  2.5_ark*cos(rho)/mX
     pseudo(1689) =  1.0_ark/mX
     pseudo(1690) =  1.25_ark*cos(rho)/mX
     pseudo(1740) =  1.25_ark*cos(rho)/mX
     pseudo(1741) =  1.0_ark/mX

     !
   end subroutine  MLkinetic_x2y2_bisect_EKE_sinrho



  !
  ! Defining kinetic energy function: sparse representation, rigid congiguration 
  ! This is EKE generated using Maple for a bisecting frame with bond-length-angle and 
  ! removed singularity by combining U sin(theta1)sin(theta2) with dG/d theta_i and multiplying muzz and G66 by sin(theta1)^2 sin(theta2)^2 
  ! and muxz and muyz by U sin(theta1)sin(theta2. 
  !
  subroutine MLkinetic_compact_x2y2_bisect_EKE_sinrho_rigid(rho,nmodes,ntermmax,ng_vib,ng_rot,ng_cor,npseudo,&
                                                            g_vib,g_rot,g_cor,pseudo,ig_vib,ig_rot,ig_cor,ipseudo)
   !
   use accuracy
   !
   real(ark),intent(in)   ::  rho
   integer(ik),intent(in) ::  nmodes
   integer(ik),intent(in) ::  ntermmax
   integer(ik),intent(out)::  ng_vib(nmodes,nmodes),ng_rot(3,3),ng_cor(nmodes,3),npseudo
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,ntermmax),g_rot(3,3,ntermmax),g_cor(nmodes,3,ntermmax),pseudo(ntermmax)
   integer(ik),intent(out)::  ig_vib(nmodes,nmodes,ntermmax,nmodes),ig_rot(3,3,ntermmax,nmodes),&
                              ig_cor(nmodes,3,ntermmax,nmodes),ipseudo(ntermmax,nmodes)
   !
   !type(FLpolynomT),pointer ::  g_vib(:,:) 
   !
   real(ark)            :: mX,mY
   integer(ik) :: info,Nterms
   logical :: check_sizes  = .false.
     !
     if (manifold==1) then
       write(out,"('MLkinetic_compact_x2y2_bisect_EKE_sinrho_rigid-error: can be used with rigid case only')")
       stop 'MLkinetic_compact_x2y2_bisect_EKE_sinrho_rigid can be used only with npoints=0'
     endif
     !
     check_sizes  = .false.
     !
     if (npseudo<0) check_sizes  = .true.
     !
     mX = molec%AtomMasses(1)
     mY = molec%AtomMasses(3)
     !
     Ng_vib = 0
     Ng_rot = 0
     Ng_cor = 0
     !
     Ng_vib(1,1) = 1
     Ng_vib(1,2) = 1
     Ng_vib(1,3) = 1
     Ng_vib(1,4) = 1
     Ng_vib(1,5) = 1
     Ng_vib(2,1) = 1
     Ng_vib(2,2) = 1
     Ng_vib(2,4) = 1
     Ng_vib(2,5) = 2
     Ng_vib(2,6) = 2
     Ng_vib(3,1) = 1
     Ng_vib(3,3) = 1
     Ng_vib(3,4) = 2
     Ng_vib(3,5) = 1
     Ng_vib(3,6) = 2
     Ng_vib(4,1) = 1
     Ng_vib(4,2) = 1
     Ng_vib(4,3) = 2
     Ng_vib(4,4) = 3
     Ng_vib(4,5) = 6
     Ng_vib(4,6) = 3
     Ng_vib(5,1) = 1
     Ng_vib(5,2) = 2
     Ng_vib(5,3) = 1
     Ng_vib(5,4) = 6
     Ng_vib(5,5) = 3
     Ng_vib(5,6) = 3
     Ng_vib(6,2) = 2
     Ng_vib(6,3) = 2
     Ng_vib(6,4) = 3
     Ng_vib(6,5) = 3
     Ng_vib(6,6) = 17
     Ng_rot(1,1) = 1
     Ng_rot(1,3) = 4
     Ng_rot(2,2) = 1
     Ng_rot(2,3) = 4
     Ng_rot(3,1) = 4
     Ng_rot(3,2) = 4
     Ng_rot(3,3) = 17
     Ng_cor(2,1) = 1
     Ng_cor(2,2) = 1
     Ng_cor(2,3) = 2
     Ng_cor(3,1) = 1
     Ng_cor(3,2) = 1
     Ng_cor(3,3) = 2
     Ng_cor(4,1) = 2
     Ng_cor(4,2) = 2
     Ng_cor(4,3) = 3
     Ng_cor(5,1) = 2
     Ng_cor(5,2) = 2
     Ng_cor(5,3) = 3
     Ng_cor(6,1) = 4
     Ng_cor(6,2) = 4
     Ng_cor(6,3) = 10
     Npseudo     = 8
     !
     Npseudo = 8
     !
     if (check_sizes) return 
     !
     ig_vib(1,1,1,:) = (/0,0,0,0,0,0/)
     ig_vib(1,1,1,:) = (/0,0,0,0,0,0/)
     ig_vib(1,2,1,:) = (/0,0,0,2,0,0/)
     ig_vib(1,3,1,:) = (/0,0,0,0,2,0/)
     ig_vib(1,4,1,:) = (/0,2,0,1,0,0/)
     ig_vib(1,5,1,:) = (/0,0,2,0,1,0/)
     ig_vib(2,1,1,:) = (/0,0,0,2,0,0/)
     ig_vib(2,2,1,:) = (/0,0,0,0,0,0/)
     ig_vib(2,4,1,:) = (/2,0,0,1,0,0/)
     ig_vib(2,5,1,:) = (/2,0,0,1,0,0/)
     ig_vib(2,5,2,:) = (/2,0,0,1,0,3/)
     ig_vib(2,6,1,:) = (/2,0,0,0,2,4/)
     ig_vib(2,6,2,:) = (/2,0,0,3,2,4/)
     ig_vib(3,1,1,:) = (/0,0,0,0,2,0/)
     ig_vib(3,3,1,:) = (/0,0,0,0,0,0/)
     ig_vib(3,4,1,:) = (/2,0,0,0,1,0/)
     ig_vib(3,4,2,:) = (/2,0,0,0,1,3/)
     ig_vib(3,5,1,:) = (/2,0,0,0,1,0/)
     ig_vib(3,6,1,:) = (/2,0,0,2,0,4/)
     ig_vib(3,6,2,:) = (/2,0,0,2,3,4/)
     ig_vib(4,1,1,:) = (/0,2,0,1,0,0/)
     ig_vib(4,2,1,:) = (/2,0,0,1,0,0/)
     ig_vib(4,3,1,:) = (/2,0,0,0,1,0/)
     ig_vib(4,3,2,:) = (/2,0,0,0,1,3/)
     ig_vib(4,4,1,:) = (/0,1,0,0,0,0/)
     ig_vib(4,4,2,:) = (/1,0,0,0,0,0/)
     ig_vib(4,4,3,:) = (/2,2,0,2,0,0/)
     ig_vib(4,5,1,:) = (/1,0,0,0,0,0/)
     ig_vib(4,5,2,:) = (/1,0,0,0,0,3/)
     ig_vib(4,5,3,:) = (/2,0,2,0,2,0/)
     ig_vib(4,5,4,:) = (/2,2,0,2,0,0/)
     ig_vib(4,5,5,:) = (/2,0,2,0,2,3/)
     ig_vib(4,5,6,:) = (/2,2,0,2,0,3/)
     ig_vib(4,6,1,:) = (/1,0,0,1,2,4/)
     ig_vib(4,6,2,:) = (/2,0,2,1,0,4/)
     ig_vib(4,6,3,:) = (/2,2,0,4,2,4/)
     ig_vib(5,1,1,:) = (/0,0,2,0,1,0/)
     ig_vib(5,2,1,:) = (/2,0,0,1,0,0/)
     ig_vib(5,2,2,:) = (/2,0,0,1,0,3/)
     ig_vib(5,3,1,:) = (/2,0,0,0,1,0/)
     ig_vib(5,4,1,:) = (/1,0,0,0,0,0/)
     ig_vib(5,4,2,:) = (/1,0,0,0,0,3/)
     ig_vib(5,4,3,:) = (/2,0,2,0,2,0/)
     ig_vib(5,4,4,:) = (/2,2,0,2,0,0/)
     ig_vib(5,4,5,:) = (/2,0,2,0,2,3/)
     ig_vib(5,4,6,:) = (/2,2,0,2,0,3/)
     ig_vib(5,5,1,:) = (/0,0,1,0,0,0/)
     ig_vib(5,5,2,:) = (/1,0,0,0,0,0/)
     ig_vib(5,5,3,:) = (/2,0,2,0,2,0/)
     ig_vib(5,6,1,:) = (/1,0,0,2,1,4/)
     ig_vib(5,6,2,:) = (/2,2,0,0,1,4/)
     ig_vib(5,6,3,:) = (/2,0,2,2,4,4/)
     ig_vib(6,2,1,:) = (/2,0,0,0,2,4/)
     ig_vib(6,2,2,:) = (/2,0,0,3,2,4/)
     ig_vib(6,3,1,:) = (/2,0,0,2,0,4/)
     ig_vib(6,3,2,:) = (/2,0,0,2,3,4/)
     ig_vib(6,4,1,:) = (/1,0,0,1,2,4/)
     ig_vib(6,4,2,:) = (/2,0,2,1,0,4/)
     ig_vib(6,4,3,:) = (/2,2,0,4,2,4/)
     ig_vib(6,5,1,:) = (/1,0,0,2,1,4/)
     ig_vib(6,5,2,:) = (/2,2,0,0,1,4/)
     ig_vib(6,5,3,:) = (/2,0,2,2,4,4/)
     ig_vib(6,6,1,:) = (/0,0,1,0,0,0/)
     ig_vib(6,6,2,:) = (/0,1,0,0,0,0/)
     ig_vib(6,6,3,:) = (/0,0,1,3,0,0/)
     ig_vib(6,6,4,:) = (/0,1,0,0,3,0/)
     ig_vib(6,6,5,:) = (/1,0,0,0,3,0/)
     ig_vib(6,6,6,:) = (/1,0,0,3,0,0/)
     ig_vib(6,6,7,:) = (/2,0,2,0,2,0/)
     ig_vib(6,6,8,:) = (/2,2,0,2,0,0/)
     ig_vib(6,6,9,:) = (/1,0,0,3,3,0/)
     ig_vib(6,6,10,:) = (/1,0,0,4,4,0/)
     ig_vib(6,6,11,:) = (/2,0,2,3,2,0/)
     ig_vib(6,6,12,:) = (/2,0,2,4,1,0/)
     ig_vib(6,6,13,:) = (/2,2,0,1,4,0/)
     ig_vib(6,6,14,:) = (/2,2,0,2,3,0/)
     ig_vib(6,6,15,:) = (/1,0,0,4,4,3/)
     ig_vib(6,6,16,:) = (/2,0,2,4,1,3/)
     ig_vib(6,6,17,:) = (/2,2,0,1,4,3/)
     ig_rot(1,1,1,:) = (/1,0,0,0,0,0/)
     ig_rot(1,3,1,:) = (/1,0,0,1,2,2/)
     ig_rot(1,3,2,:) = (/1,0,0,2,1,2/)
     ig_rot(1,3,3,:) = (/2,0,2,1,0,2/)
     ig_rot(1,3,4,:) = (/2,2,0,0,1,2/)
     ig_rot(2,2,1,:) = (/1,0,0,0,0,0/)
     ig_rot(2,3,1,:) = (/1,0,0,1,2,1/)
     ig_rot(2,3,2,:) = (/1,0,0,2,1,1/)
     ig_rot(2,3,3,:) = (/2,0,2,1,0,1/)
     ig_rot(2,3,4,:) = (/2,2,0,0,1,1/)
     ig_rot(3,1,1,:) = (/1,0,0,1,2,2/)
     ig_rot(3,1,2,:) = (/1,0,0,2,1,2/)
     ig_rot(3,1,3,:) = (/2,0,2,1,0,2/)
     ig_rot(3,1,4,:) = (/2,2,0,0,1,2/)
     ig_rot(3,2,1,:) = (/1,0,0,1,2,1/)
     ig_rot(3,2,2,:) = (/1,0,0,2,1,1/)
     ig_rot(3,2,3,:) = (/2,0,2,1,0,1/)
     ig_rot(3,2,4,:) = (/2,2,0,0,1,1/)
     ig_rot(3,3,1,:) = (/0,0,1,0,0,0/)
     ig_rot(3,3,2,:) = (/0,1,0,0,0,0/)
     ig_rot(3,3,3,:) = (/0,0,1,3,0,0/)
     ig_rot(3,3,4,:) = (/0,1,0,0,3,0/)
     ig_rot(3,3,5,:) = (/1,0,0,0,3,0/)
     ig_rot(3,3,6,:) = (/1,0,0,3,0,0/)
     ig_rot(3,3,7,:) = (/2,0,2,0,2,0/)
     ig_rot(3,3,8,:) = (/2,2,0,2,0,0/)
     ig_rot(3,3,9,:) = (/1,0,0,3,3,0/)
     ig_rot(3,3,10,:) = (/1,0,0,4,4,0/)
     ig_rot(3,3,11,:) = (/2,0,2,3,2,0/)
     ig_rot(3,3,12,:) = (/2,0,2,4,1,0/)
     ig_rot(3,3,13,:) = (/2,2,0,1,4,0/)
     ig_rot(3,3,14,:) = (/2,2,0,2,3,0/)
     ig_rot(3,3,15,:) = (/1,0,0,4,4,3/)
     ig_rot(3,3,16,:) = (/2,0,2,4,1,3/)
     ig_rot(3,3,17,:) = (/2,2,0,1,4,3/)
     ig_cor(2,1,1,:) = (/2,0,0,1,0,1/)
     ig_cor(2,2,1,:) = (/2,0,0,1,0,2/)
     ig_cor(2,3,1,:) = (/2,0,0,0,2,4/)
     ig_cor(2,3,2,:) = (/2,0,0,3,2,4/)
     ig_cor(3,1,1,:) = (/2,0,0,0,1,1/)
     ig_cor(3,2,1,:) = (/2,0,0,0,1,2/)
     ig_cor(3,3,1,:) = (/2,0,0,2,0,4/)
     ig_cor(3,3,2,:) = (/2,0,0,2,3,4/)
     ig_cor(4,1,1,:) = (/1,0,0,0,0,1/)
     ig_cor(4,1,2,:) = (/2,2,0,2,0,1/)
     ig_cor(4,2,1,:) = (/1,0,0,0,0,2/)
     ig_cor(4,2,2,:) = (/2,2,0,2,0,2/)
     ig_cor(4,3,1,:) = (/1,0,0,1,2,4/)
     ig_cor(4,3,2,:) = (/2,0,2,1,0,4/)
     ig_cor(4,3,3,:) = (/2,2,0,4,2,4/)
     ig_cor(5,1,1,:) = (/1,0,0,0,0,1/)
     ig_cor(5,1,2,:) = (/2,0,2,0,2,1/)
     ig_cor(5,2,1,:) = (/1,0,0,0,0,2/)
     ig_cor(5,2,2,:) = (/2,0,2,0,2,2/)
     ig_cor(5,3,1,:) = (/1,0,0,2,1,4/)
     ig_cor(5,3,2,:) = (/2,2,0,0,1,4/)
     ig_cor(5,3,3,:) = (/2,0,2,2,4,4/)
     ig_cor(6,1,1,:) = (/1,0,0,1,2,2/)
     ig_cor(6,1,2,:) = (/1,0,0,2,1,2/)
     ig_cor(6,1,3,:) = (/2,0,2,1,0,2/)
     ig_cor(6,1,4,:) = (/2,2,0,0,1,2/)
     ig_cor(6,2,1,:) = (/1,0,0,1,2,1/)
     ig_cor(6,2,2,:) = (/1,0,0,2,1,1/)
     ig_cor(6,2,3,:) = (/2,0,2,1,0,1/)
     ig_cor(6,2,4,:) = (/2,2,0,0,1,1/)
     ig_cor(6,3,1,:) = (/0,0,1,0,0,0/)
     ig_cor(6,3,2,:) = (/0,1,0,0,0,0/)
     ig_cor(6,3,3,:) = (/0,0,1,3,0,0/)
     ig_cor(6,3,4,:) = (/0,1,0,0,3,0/)
     ig_cor(6,3,5,:) = (/1,0,0,0,3,0/)
     ig_cor(6,3,6,:) = (/1,0,0,3,0,0/)
     ig_cor(6,3,7,:) = (/2,0,2,0,2,0/)
     ig_cor(6,3,8,:) = (/2,2,0,2,0,0/)
     ig_cor(6,3,9,:) = (/2,0,2,3,2,0/)
     ig_cor(6,3,10,:) = (/2,2,0,2,3,0/)
     ipseudo(1,:) = (/1,0,0,2,2,0/)
     ipseudo(2,:) = (/1,0,0,2,2,3/)
     ipseudo(3,:) = (/2,0,2,1,4,0/)
     ipseudo(4,:) = (/2,0,2,2,3,0/)
     ipseudo(5,:) = (/2,2,0,3,2,0/)
     ipseudo(6,:) = (/2,2,0,4,1,0/)
     ipseudo(7,:) = (/2,0,2,2,3,3/)
     ipseudo(8,:) = (/2,2,0,3,2,3/)
     !
     g_vib(1,1,1) =  2./mX
     g_vib(1,2,1) =  -1./mX
     g_vib(1,3,1) =  -1./mX
     g_vib(1,4,1) =  1/mX
     g_vib(1,5,1) =  1/mX
     g_vib(2,1,1) =  -1./mX
     g_vib(2,2,1) =  (mX+mY)/mY/mX
     g_vib(2,4,1) =  1/mX
     g_vib(2,5,1) =  1/mX
     g_vib(2,5,2) =  -2./mX
     g_vib(2,6,1) =  2./mX
     g_vib(2,6,2) =  -2./mX
     g_vib(3,1,1) =  -1./mX
     g_vib(3,3,1) =  (mY**3+10.*mX**2*mY+7.*mY**2*mX+4.*mX**3)/(2.*mX+mY)**2/mX/mY
     g_vib(3,4,1) =  1/mX
     g_vib(3,4,2) =  -2./mX
     g_vib(3,5,1) =  1/mX
     g_vib(3,6,1) =  2./mX
     g_vib(3,6,2) =  -2./mX
     g_vib(4,1,1) =  1/mX
     g_vib(4,2,1) =  1/mX
     g_vib(4,3,1) =  1/mX
     g_vib(4,3,2) =  -2./mX
     g_vib(4,4,1) =  (mY+mX)/mX/mY
     g_vib(4,4,2) =  2./mX
     g_vib(4,4,3) =  2./mX
     g_vib(4,5,1) =  2./mX
     g_vib(4,5,2) =  -4./mX
     g_vib(4,5,3) =  1/mX
     g_vib(4,5,4) =  1/mX
     g_vib(4,5,5) =  -2./mX
     g_vib(4,5,6) =  -2./mX
     g_vib(4,6,1) =  4./mX
     g_vib(4,6,2) =  2./mX
     g_vib(4,6,3) =  2./mX
     g_vib(5,1,1) =  1/mX
     g_vib(5,2,1) =  1/mX
     g_vib(5,2,2) =  -2./mX
     g_vib(5,3,1) =  1/mX
     g_vib(5,4,1) =  2./mX
     g_vib(5,4,2) =  -4./mX
     g_vib(5,4,3) =  1/mX
     g_vib(5,4,4) =  1/mX
     g_vib(5,4,5) =  -2./mX
     g_vib(5,4,6) =  -2./mX
     g_vib(5,5,1) =  (7.*mY**2*mX+10.*mX**2*mY+mY**3+4.*mX**3)/(2.*mX+mY)**2/mX/mY
     g_vib(5,5,2) =  2./mX
     g_vib(5,5,3) =  2./mX
     g_vib(5,6,1) =  4./mX
     g_vib(5,6,2) =  2./mX
     g_vib(5,6,3) =  2./mX
     g_vib(6,2,1) =  2./mX
     g_vib(6,2,2) =  -2./mX
     g_vib(6,3,1) =  2./mX
     g_vib(6,3,2) =  -2./mX
     g_vib(6,4,1) =  4./mX
     g_vib(6,4,2) =  2./mX
     g_vib(6,4,3) =  2./mX
     g_vib(6,5,1) =  4./mX
     g_vib(6,5,2) =  2./mX
     g_vib(6,5,3) =  2./mX
     g_vib(6,6,1) =  -1.*(-7.*mY**2*mX-10.*mX**2*mY-4.*mX**3-1.*mY**3)/(2.*mX+mY)**2/mY/mX
     g_vib(6,6,2) =  -1.*(-1.*mX-1.*mY)/mY/mX
     g_vib(6,6,3) =  -1.*(mY**3+4.*mX**3+7.*mY**2*mX+10.*mX**2*mY)/(2.*mX+mY)**2/mY/mX
     g_vib(6,6,4) =  -1.*(mX+mY)/mX/mY
     g_vib(6,6,5) =  2./mX
     g_vib(6,6,6) =  2./mX
     g_vib(6,6,7) =  2./mX
     g_vib(6,6,8) =  2./mX
     g_vib(6,6,9) =  -4./mX
     g_vib(6,6,10) =  -4./mX
     g_vib(6,6,11) =  -2./mX
     g_vib(6,6,12) =  -2./mX
     g_vib(6,6,13) =  -2./mX
     g_vib(6,6,14) =  -2./mX
     g_vib(6,6,15) =  8./mX
     g_vib(6,6,16) =  4./mX
     g_vib(6,6,17) =  4./mX
     g_rot(1,1,1) =  2./mX
     g_rot(1,3,1) =  -1./mX
     g_rot(1,3,2) =  1/mX
     g_rot(1,3,3) =  -.500/mX
     g_rot(1,3,4) =  .500/mX
     g_rot(2,2,1) =  2./mX
     g_rot(2,3,1) =  1/mX
     g_rot(2,3,2) =  1/mX
     g_rot(2,3,3) =  .500/mX
     g_rot(2,3,4) =  .500/mX
     g_rot(3,1,1) =  -1./mX
     g_rot(3,1,2) =  1/mX
     g_rot(3,1,3) =  -.500/mX
     g_rot(3,1,4) =  .500/mX
     g_rot(3,2,1) =  1/mX
     g_rot(3,2,2) =  1/mX
     g_rot(3,2,3) =  .500/mX
     g_rot(3,2,4) =  .500/mX
     g_rot(3,3,1) =  -.250*(-10.*mX**2*mY-7.*mY**2*mX-4.*mX**3-1.*mY**3)/(2.*mX+mY)**2/mX/mY
     g_rot(3,3,2) =  -.250*(-1.*mX-1.*mY)/mX/mY
     g_rot(3,3,3) =  -.250*(4.*mX**3+10.*mX**2*mY+7.*mY**2*mX+mY**3)/(2.*mX+mY)**2/mX/mY
     g_rot(3,3,4) =  -.250*(mX+mY)/mX/mY
     g_rot(3,3,5) =  .500/mX
     g_rot(3,3,6) =  .500/mX
     g_rot(3,3,7) =  .500/mX
     g_rot(3,3,8) =  .500/mX
     g_rot(3,3,9) =  -1./mX
     g_rot(3,3,10) =  1/mX
     g_rot(3,3,11) =  -.500/mX
     g_rot(3,3,12) =  .500/mX
     g_rot(3,3,13) =  .500/mX
     g_rot(3,3,14) =  -.500/mX
     g_rot(3,3,15) =  -2./mX
     g_rot(3,3,16) =  -1./mX
     g_rot(3,3,17) =  -1./mX
     g_cor(2,1,1) =  1/mX
     g_cor(2,2,1) =  -1./mX
     g_cor(2,3,1) =  -1./mX
     g_cor(2,3,2) =  1/mX
     g_cor(3,1,1) =  1/mX
     g_cor(3,2,1) =  1/mX
     g_cor(3,3,1) =  1/mX
     g_cor(3,3,2) =  -1./mX
     g_cor(4,1,1) =  2./mX
     g_cor(4,1,2) =  1/mX
     g_cor(4,2,1) =  -2./mX
     g_cor(4,2,2) =  -1./mX
     g_cor(4,3,1) =  -2./mX
     g_cor(4,3,2) =  -1./mX
     g_cor(4,3,3) =  -1./mX
     g_cor(5,1,1) =  2./mX
     g_cor(5,1,2) =  1/mX
     g_cor(5,2,1) =  2./mX
     g_cor(5,2,2) =  1/mX
     g_cor(5,3,1) =  2./mX
     g_cor(5,3,2) =  1/mX
     g_cor(5,3,3) =  1/mX
     g_cor(6,1,1) =  2./mX
     g_cor(6,1,2) =  2./mX
     g_cor(6,1,3) =  1/mX
     g_cor(6,1,4) =  1/mX
     g_cor(6,2,1) =  -2./mX
     g_cor(6,2,2) =  2./mX
     g_cor(6,2,3) =  -1./mX
     g_cor(6,2,4) =  1/mX
     g_cor(6,3,1) =  .500*(-7.*mY**2*mX-10.*mX**2*mY-1.*mY**3-4.*mX**3)/(2.*mX+mY)**2/mX/mY
     g_cor(6,3,2) =  -.500*(-1.*mX-1.*mY)/mX/mY
     g_cor(6,3,3) =  .500*(mY**3+4.*mX**3+10.*mX**2*mY+7.*mY**2*mX)/(2.*mX+mY)**2/mX/mY
     g_cor(6,3,4) =  -.500*(mX+mY)/mY/mX
     g_cor(6,3,5) =  -1./mX
     g_cor(6,3,6) =  1/mX
     g_cor(6,3,7) =  -1./mX
     g_cor(6,3,8) =  1/mX
     g_cor(6,3,9) =  1/mX
     g_cor(6,3,10) =  -1./mX
     pseudo(1) =  -2.50/mX
     pseudo(2) =  5./mX
     pseudo(3) =  1/mX
     pseudo(4) =  -1.25/mX
     pseudo(5) =  -1.25/mX
     pseudo(6) =  1/mX
     pseudo(7) =  2.50/mX
     pseudo(8) =  2.50/mX
     !
   end subroutine  MLkinetic_compact_x2y2_bisect_EKE_sinrho_rigid





end module kin_x2y2
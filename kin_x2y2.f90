!
!  This unit defines all kinetic routines for a triatomic molecule of the X2Y2 type
!
module kin_x2y2
  use accuracy
  use moltype

  implicit none

  public MLkinetic_x2y2_bisect_EKE_sinrho
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
     g_vib(1,1,1) =  2./mX
     g_vib(1,2,9) =  -1./mX
     g_vib(1,3,7) =  -1./mX
     g_vib(1,4,39) =  1/mX
     g_vib(1,5,29) =  1/mX
     g_vib(2,1,9) =  -1./mX
     g_vib(2,2,1) =  (mY+mX)/mY/mX
     g_vib(2,4,53) =  1/mX
     g_vib(2,5,53) =  -1.*cos(rho)/mX
     g_vib(2,6,112) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(2,6,670) =  -2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(3,1,7) =  -1./mX
     g_vib(3,3,1) =  (10.*mX**2*mY+7.*mY**2*mX+mY**3+4.*mX**3)/(2.*mX+mY)**2/mX/mY
     g_vib(3,4,52) =  -1.*cos(rho)/mX
     g_vib(3,5,52) =  1/mX
     g_vib(3,6,114) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(3,6,669) =  -2.*cos(rho_2)*sin(rho_2)/mX
     g_vib(4,1,39) =  1/mX
     g_vib(4,2,53) =  1/mX
     g_vib(4,3,52) =  -1.*cos(rho)/mX
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
     g_cor(2,2,53) =  -1.*cos(rho_2)/mX
     g_cor(2,3,112) =  -1.*cos(rho_2)*sin(rho_2)/mX
     g_cor(2,3,670) =  cos(rho_2)*sin(rho_2)/mX
     g_cor(3,1,52) =  sin(rho_2)/mX
     g_cor(3,2,52) =  cos(rho_2)/mX
     g_cor(3,3,114) =  cos(rho_2)*sin(rho_2)/mX
     g_cor(3,3,669) =  -1.*cos(rho_2)*sin(rho_2)/mX
     g_cor(4,1,6) =  2.*sin(rho_2)/mX
     g_cor(4,1,420) =  sin(rho_2)/mX
     g_cor(4,2,6) =  -2.*cos(rho_2)/mX
     g_cor(4,2,420) =  -1.*cos(rho_2)/mX
     g_cor(4,3,93) =  -2.*cos(rho_2)*sin(rho_2)/mX
     g_cor(4,3,226) =  -1.*cos(rho_2)*sin(rho_2)/mX
     g_cor(4,3,2594) =  -1.*cos(rho_2)*sin(rho_2)/mX
     g_cor(5,1,6) =  2.*sin(rho_2)/mX
     g_cor(5,1,402) =  sin(rho_2)/mX
     g_cor(5,2,6) =  2.*cos(rho_2)/mX
     g_cor(5,2,402) =  cos(rho_2)/mX
     g_cor(5,3,94) =  2.*cos(rho_2)*sin(rho_2)/mX
     g_cor(5,3,234) =  cos(rho_2)*sin(rho_2)/mX
     g_cor(5,3,2528) =  cos(rho_2)*sin(rho_2)/mX
     g_cor(6,1,93) =  2.*cos(rho_2)/mX
     g_cor(6,1,94) =  2.*cos(rho_2)/mX
     g_cor(6,1,226) =  cos(rho_2)/mX
     g_cor(6,1,234) =  cos(rho_2)/mX
     g_cor(6,2,93) =  -2.*sin(rho_2)/mX
     g_cor(6,2,94) =  2.*sin(rho_2)/mX
     g_cor(6,2,226) =  -1.*sin(rho_2)/mX
     g_cor(6,2,234) =  sin(rho_2)/mX
     g_cor(6,3,4) =  -.5000000000*(10.*mX**2*mY+7.*mY**2*mX+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_cor(6,3,5) =  .5000000000*(mX+mY)/mX/mY
     g_cor(6,3,65) =  .5000000000*(10.*mX**2*mY+7.*mY**2*mX+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_cor(6,3,72) =  -.5000000000*(mX+mY)/mX/mY
     g_cor(6,3,92) =  -1./mX
     g_cor(6,3,95) =  1/mX
     g_cor(6,3,402) =  -1./mX
     g_cor(6,3,420) =  1/mX
     g_cor(6,3,1691) =  1/mX
     g_cor(6,3,1739) =  -1./mX


     pseudo(185) =  2.500000000*cos(rho)/mX
     pseudo(1689) =  1/mX
     pseudo(1690) =  1.250000000*cos(rho)/mX
     pseudo(1740) =  1.250000000*cos(rho)/mX
     pseudo(1741) =  1/mX

     !
   end subroutine  MLkinetic_x2y2_bisect_EKE_sinrho






end module kin_x2y2
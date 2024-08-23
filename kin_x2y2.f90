!
!  This unit defines all kinetic routines for a triatomic molecule of the X2Y2 type
!
module kin_x2y2
  use accuracy
  use moltype
  use timer

  implicit none

  public MLkinetic_x2y2_bisect_EKE_sinrho,MLkinetic_compact_x2y2_bisect_EKE_sinrho
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
  ! Defining kinetic energy function: sparse representation. 
  ! This is EKE generated using Maple for a bisecting frame with bond-length-angle and 
  ! removed singularity by combining U sin(theta1)sin(theta2) with dG/d theta_i and multiplying muzz and G66 by sin(theta1)^2 sin(theta2)^2 
  ! and muxz and muyz by U sin(theta1)sin(theta2. 
  !
  subroutine MLkinetic_compact_x2y2_bisect_EKE_sinrho(nmodes,rho,g_vib,g_rot,g_cor,pseudo,ig_vib,ig_rot,ig_cor,ipseudo,&
                                                     ng_vib,ng_rot,ng_cor,npseudo)
   !
   integer(ik),intent(in) ::  nmodes
   real(ark),intent(in)   ::  rho
   real(ark),allocatable,intent(out)    ::  g_vib(:,:,:),g_rot(:,:,:),g_cor(:,:,:),pseudo(:)
   integer(ik),allocatable,intent(out)  ::  ig_vib(:,:,:),ig_rot(:,:,:),ig_cor(:,:,:),ipseudo(:)
   integer(ik),intent(out)  ::  ng_vib(nmodes,nmodes),ng_rot(3,3),ng_cor(nmodes,3),npseudo
   !
   !type(FLpolynomT),pointer ::  g_vib(:,:) 
   !
   real(ark)            :: mX,mY,rho_2
   real(ark),parameter  :: rho_threshold = 0.01_rk
   integer(ik) :: N_g_rot(3,3),N_g_cor(6,3),N_g_vib(6,6),N_pseudo,info,i1,i2,Nterms
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
     Ng_vib(1,1) = 1
     Ng_vib(1,2) = 1
     Ng_vib(1,3) = 1
     Ng_vib(1,4) = 1
     Ng_vib(1,5) = 1
     Ng_vib(1,6) = 0
     Ng_vib(2,1) = 1
     Ng_vib(2,2) = 1
     Ng_vib(2,3) = 0
     Ng_vib(2,4) = 1
     Ng_vib(2,5) = 1
     Ng_vib(2,6) = 2
     Ng_vib(3,1) = 1
     Ng_vib(3,2) = 0
     Ng_vib(3,3) = 1
     Ng_vib(3,4) = 1
     Ng_vib(3,5) = 1
     Ng_vib(3,6) = 2
     Ng_vib(4,1) = 1
     Ng_vib(4,2) = 1
     Ng_vib(4,3) = 1
     Ng_vib(4,4) = 3
     Ng_vib(4,5) = 3
     Ng_vib(4,6) = 3
     Ng_vib(5,1) = 1
     Ng_vib(5,2) = 1
     Ng_vib(5,3) = 1
     Ng_vib(5,4) = 3
     Ng_vib(5,5) = 3
     Ng_vib(5,6) = 3
     Ng_vib(6,1) = 0
     Ng_vib(6,2) = 2
     Ng_vib(6,3) = 2
     Ng_vib(6,4) = 3
     Ng_vib(6,5) = 3
     Ng_vib(6,6) = 14
     !
     Ng_rot(1,1) = 1
     Ng_rot(1,2) = 0
     Ng_rot(1,3) = 4
     Ng_rot(2,1) = 0
     Ng_rot(2,2) = 1
     Ng_rot(2,3) = 4
     Ng_rot(3,1) = 4
     Ng_rot(3,2) = 4
     Ng_rot(3,3) = 14
     !
     Ng_cor(1,1) = 0
     Ng_cor(1,2) = 0
     Ng_cor(1,3) = 0
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
     !
     Npseudo = 5
     !
     Nterms = maxval(N_g_vib)
     !
     if (.not.allocated(g_vib)) then 
        !
        if (Nterms>0) then
            allocate(g_vib(nmodes,nmodes,Nterms),ig_vib(nmodes,nmodes,Nterms),stat=info)
          else
            allocate(g_vib(nmodes,nmodes,0:0),ig_vib(nmodes,nmodes,0:0),stat=info)
        endif
        call ArrayStart('kinetic_on_grid-fields',info,size(g_vib),kind(g_vib))
        call ArrayStart('kinetic_on_grid-fields',info,size(ig_vib),kind(ig_vib))
        !
        Nterms = maxval(N_g_cor)
        !
        if (Nterms>0) then
            allocate(g_cor(nmodes,3,Nterms),ig_cor(nmodes,3,Nterms),stat=info)
          else
            allocate(g_cor(nmodes,3,0:0),ig_cor(nmodes,3,0:0),stat=info)
        endif
        call ArrayStart('kinetic_on_grid-fields',info,size(g_cor),kind(g_cor))
        call ArrayStart('kinetic_on_grid-fields',info,size(ig_cor),kind(ig_cor))
        !
        Nterms = maxval(N_g_rot)
        !
        if (Nterms>0) then
            allocate(g_rot(3,3,Nterms),ig_rot(3,3,Nterms),stat=info)
          else
            allocate(g_rot(3,3,0:0),ig_rot(3,3,0:0),stat=info)
        endif
        call ArrayStart('kinetic_on_grid-fields',info,size(g_rot),kind(g_rot))
        call ArrayStart('kinetic_on_grid-fields',info,size(ig_rot),kind(ig_rot))
        !
        Nterms  = N_pseudo
        !
        if (Nterms>0) then
            allocate(pseudo(Nterms),ipseudo(Nterms),stat=info)
          else
            allocate(pseudo(0:0),ipseudo(0:0),stat=info)
        endif
        call ArrayStart('kinetic_on_grid-fields',info,size(pseudo),kind(pseudo))
        call ArrayStart('kinetic_on_grid-fields',info,size(ipseudo),kind(ipseudo))
        !
     endif
     !
     g_vib(1,1,1) =  2./mX
     g_vib(1,2,1) =  -1./mX
     g_vib(1,3,1) =  -1./mX
     g_vib(1,4,1) =  1/mX
     g_vib(1,5,1) =  1/mX
     g_vib(2,1,1) =  -1./mX
     g_vib(2,2,1) =  (mX+mY)/mY/mX
     g_vib(2,4,1) =  1/mX
     g_vib(2,5,1) =  -1.*cos(rho)/mX
     g_vib(2,6,1) =  1/mX*sin(rho)
     g_vib(2,6,2) =  -1./mX*sin(rho)
     g_vib(3,1,1) =  -1./mX
     g_vib(3,3,1) =  (10.*mX**2*mY+7.*mY**2*mX+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_vib(3,4,1) =  -1.*cos(rho)/mX
     g_vib(3,5,1) =  1/mX
     g_vib(3,6,1) =  1/mX*sin(rho)
     g_vib(3,6,2) =  -1./mX*sin(rho)
     g_vib(4,1,1) =  1/mX
     g_vib(4,2,1) =  1/mX
     g_vib(4,3,1) =  -1.*cos(rho)/mX
     g_vib(4,4,1) =  (mY+mX)/mX/mY
     g_vib(4,4,2) =  2./mX
     g_vib(4,4,3) =  2./mX
     g_vib(4,5,1) =  -2.*cos(rho)/mX
     g_vib(4,5,2) =  -1.*cos(rho)/mX
     g_vib(4,5,3) =  -1.*cos(rho)/mX
     g_vib(4,6,1) =  2./mX*sin(rho)
     g_vib(4,6,2) =  1/mX*sin(rho)
     g_vib(4,6,3) =  1/mX*sin(rho)
     g_vib(5,1,1) =  1/mX
     g_vib(5,2,1) =  -1.*cos(rho)/mX
     g_vib(5,3,1) =  1/mX
     g_vib(5,4,1) =  -2.*cos(rho)/mX
     g_vib(5,4,2) =  -1.*cos(rho)/mX
     g_vib(5,4,3) =  -1.*cos(rho)/mX
     g_vib(5,5,1) =  (mY**3+4.*mX**3+7.*mY**2*mX+10.*mX**2*mY)/(2.*mX+mY)**2/mX/mY
     g_vib(5,5,2) =  2./mX
     g_vib(5,5,3) =  2./mX
     g_vib(5,6,1) =  2./mX*sin(rho)
     g_vib(5,6,2) =  1/mX*sin(rho)
     g_vib(5,6,3) =  1/mX*sin(rho)
     g_vib(6,2,1) =  1/mX*sin(rho)
     g_vib(6,2,2) =  -1./mX*sin(rho)
     g_vib(6,3,1) =  1/mX*sin(rho)
     g_vib(6,3,2) =  -1./mX*sin(rho)
     g_vib(6,4,1) =  2./mX*sin(rho)
     g_vib(6,4,2) =  1/mX*sin(rho)
     g_vib(6,4,3) =  1/mX*sin(rho)
     g_vib(6,5,1) =  2./mX*sin(rho)
     g_vib(6,5,2) =  1/mX*sin(rho)
     g_vib(6,5,3) =  1/mX*sin(rho)
     g_vib(6,6,1) =  (mY**3+4.*mX**3+7.*mY**2*mX+10.*mX**2*mY)/(2.*mX+mY)**2/mX/mY
     g_vib(6,6,2) =  (mX+mY)/mY/mX
     g_vib(6,6,3) =  -1.*(mY**3+4.*mX**3+7.*mY**2*mX+10.*mX**2*mY)/(2.*mX+mY)**2/mX/mY
     g_vib(6,6,4) =  -1.*(mX+mY)/mY/mX
     g_vib(6,6,5) =  2./mX
     g_vib(6,6,6) =  2./mX
     g_vib(6,6,7) =  2./mX
     g_vib(6,6,8) =  2./mX
     g_vib(6,6,9) =  -4./mX
     g_vib(6,6,10) =  4.*cos(rho)/mX
     g_vib(6,6,11) =  -2./mX
     g_vib(6,6,12) =  2.*cos(rho)/mX
     g_vib(6,6,13) =  2.*cos(rho)/mX
     g_vib(6,6,14) =  -2./mX
     !
     g_rot(1,1,1) =  2./mX
     g_rot(1,3,1) =  -1.*cos(rho_2)/mX
     g_rot(1,3,2) =  cos(rho_2)/mX
     g_rot(1,3,3) =  -.500*cos(rho_2)/mX
     g_rot(1,3,4) =  .500*cos(rho_2)/mX
     g_rot(2,2,1) =  2./mX
     g_rot(2,3,1) =  sin(rho_2)/mX
     g_rot(2,3,2) =  sin(rho_2)/mX
     g_rot(2,3,3) =  .500*sin(rho_2)/mX
     g_rot(2,3,4) =  .500*sin(rho_2)/mX
     g_rot(3,1,1) =  -1.*cos(rho_2)/mX
     g_rot(3,1,2) =  cos(rho_2)/mX
     g_rot(3,1,3) =  -.500*cos(rho_2)/mX
     g_rot(3,1,4) =  .500*cos(rho_2)/mX
     g_rot(3,2,1) =  sin(rho_2)/mX
     g_rot(3,2,2) =  sin(rho_2)/mX
     g_rot(3,2,3) =  .500*sin(rho_2)/mX
     g_rot(3,2,4) =  .500*sin(rho_2)/mX
     g_rot(3,3,1) =  .250*(7.*mY**2*mX+10.*mX**2*mY+4.*mX**3+mY**3)/(2.*mX+mY)**2/mY/mX
     g_rot(3,3,2) =  .250*(mX+mY)/mY/mX
     g_rot(3,3,3) =  -.250*(7.*mY**2*mX+10.*mX**2*mY+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_rot(3,3,4) =  -.250*(mX+mY)/mY/mX
     g_rot(3,3,5) =  .500/mX
     g_rot(3,3,6) =  .500/mX
     g_rot(3,3,7) =  .500/mX
     g_rot(3,3,8) =  .500/mX
     g_rot(3,3,9) =  -1./mX
     g_rot(3,3,10) =  -1.*cos(rho)/mX
     g_rot(3,3,11) =  -.500/mX
     g_rot(3,3,12) =  -.500*cos(rho)/mX
     g_rot(3,3,13) =  -.500*cos(rho)/mX
     g_rot(3,3,14) =  -.500/mX
     !
     g_cor(2,1,1) =  sin(rho_2)/mX
     g_cor(2,2,1) =  -1.*cos(rho_2)/mX
     g_cor(2,3,1) =  -.500/mX*sin(rho)
     g_cor(2,3,2) =  .500/mX*sin(rho)
     g_cor(3,1,1) =  sin(rho_2)/mX
     g_cor(3,2,1) =  cos(rho_2)/mX
     g_cor(3,3,1) =  .500/mX*sin(rho)
     g_cor(3,3,2) =  -.500/mX*sin(rho)
     g_cor(4,1,1) =  2.*sin(rho_2)/mX
     g_cor(4,1,2) =  sin(rho_2)/mX
     g_cor(4,2,1) =  -2.*cos(rho_2)/mX
     g_cor(4,2,2) =  -1.*cos(rho_2)/mX
     g_cor(4,3,1) =  -1./mX*sin(rho)
     g_cor(4,3,2) =  -.500/mX*sin(rho)
     g_cor(4,3,3) =  -.500/mX*sin(rho)
     g_cor(5,1,1) =  2.*sin(rho_2)/mX
     g_cor(5,1,2) =  sin(rho_2)/mX
     g_cor(5,2,1) =  2.*cos(rho_2)/mX
     g_cor(5,2,2) =  cos(rho_2)/mX
     g_cor(5,3,1) =  1/mX*sin(rho)
     g_cor(5,3,2) =  .500/mX*sin(rho)
     g_cor(5,3,3) =  .500/mX*sin(rho)
     g_cor(6,1,1) =  2.*cos(rho_2)/mX
     g_cor(6,1,2) =  2.*cos(rho_2)/mX
     g_cor(6,1,3) =  cos(rho_2)/mX
     g_cor(6,1,4) =  cos(rho_2)/mX
     g_cor(6,2,1) =  -2.*sin(rho_2)/mX
     g_cor(6,2,2) =  2.*sin(rho_2)/mX
     g_cor(6,2,3) =  -1.*sin(rho_2)/mX
     g_cor(6,2,4) =  sin(rho_2)/mX
     g_cor(6,3,1) =  -.500*(10.*mX**2*mY+7.*mY**2*mX+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_cor(6,3,2) =  .500*(mX+mY)/mX/mY
     g_cor(6,3,3) =  .500*(7.*mY**2*mX+10.*mX**2*mY+4.*mX**3+mY**3)/(2.*mX+mY)**2/mX/mY
     g_cor(6,3,4) =  -.500*(mX+mY)/mY/mX
     g_cor(6,3,5) =  -1./mX
     g_cor(6,3,6) =  1/mX
     g_cor(6,3,7) =  -1./mX
     g_cor(6,3,8) =  1/mX
     g_cor(6,3,9) =  1/mX
     g_cor(6,3,10) =  -1./mX
     pseudo(1) =  2.50*cos(rho)/mX
     pseudo(2) =  1/mX
     pseudo(3) =  1.25*cos(rho)/mX
     pseudo(4) =  1.25*cos(rho)/mX
     pseudo(5) =  1/mX
     !
   end subroutine  MLkinetic_compact_x2y2_bisect_EKE_sinrho





end module kin_x2y2
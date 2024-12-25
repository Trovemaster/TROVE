!
!  This unit defines all kinetic routines for a triatomic molecule of the X2Y2 type
!
module kin_abcd
  use accuracy
  use moltype
  use timer

  implicit none

  public MLkinetic_abcd_EKE_z_alpha1_singular
  
  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains


  ! Defining kinetic energy function 
  ! This is EKE generated using Maple for Y1X1X2Y2 with the x axis in the X2Y2-X1 plane, z is along X1-X2 and 
  ! alpha =Y1-X1-Y2 as the singular coordinate (can become 180) 
  ! singularity is removed by combining U sin(alpha1) with dG/d alpha_1 and multiplying G66 and G44 by sin(alpha1)^2.
  !
  subroutine MLkinetic_abcd_EKE_z_alpha1_singular(nmodes,rho,ntermmax,ng_vib,ng_rot,ng_cor,npseudo,&
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
   real(ark)            :: mX1,mX2,mY1,mY2
   integer(ik) :: info,Nterms
   integer(ik),parameter :: nlines = 33
   character(len=wl) :: constructor(nlines) = (/&
      'Mode 1 2',&
      '1 1 -1 I 1 1',&
      '2 1 -2 I 1 1',&
      'Mode 2 2',&
      '1 1 -1 I 1 1',&
      '2 1 -2 I 1 1',&
      'Mode 3 2',&
      '1 1 -1 I 1 1',&
      '2 1 -2 I 1 1',&
      'Mode 4 8',&
      ' 1 1 1 sin 1 1',&
      ' 2 1 1 cos 1 1',&
      ' 3 1 1 csc 1.0 1',&
      ' 4 1 2 csc 1.0 1',&
      ' 5 1 1 cot 1.0 1',&
      ' 6 1 2 cot 1.0 1',&
      ' 7 2 1 cos 1.0 1  2 csc 1.0 1',& 
      ' 8 1 1 cos 2.0 1',&
      'Mode 5 8',&
      ' 1 1 1 sin 1 1',&
      ' 2 1 1 cos 1 1',&
      ' 3 1 1 csc 1.0 1',&
      ' 4 1 2 csc 1.0 1',&
      ' 5 1 1 cot 1.0 1',&
      ' 6 1 2 cot 1.0 1',&
      ' 7 2 1 cos 1.0 1  2 csc 1.0 1',& 
      ' 8 1 1 cos 2.0 1',&
      'Mode 6 4',&
      ' 1 1 1 sin 1 1',&
      ' 2 1 1 cos 1 1',&
      ' 3 1 2 cos 1 1',&
      ' 4 1 1 sin 2 1',&
      'end'/)

     !
     mX1 = molec%AtomMasses(1)
     mX2 = molec%AtomMasses(2)
     mY1 = molec%AtomMasses(3)
     mY2 = molec%AtomMasses(4)
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
     Ng_vib(1,2) = 1
     Ng_vib(1,3) = 1
     Ng_vib(1,4) = 1
     Ng_vib(1,5) = 1
     Ng_vib(2,1) = 1
     Ng_vib(2,2) = 1
     Ng_vib(2,4) = 1
     Ng_vib(2,5) = 1
     Ng_vib(2,6) = 1
     Ng_vib(3,1) = 1
     Ng_vib(3,3) = 1
     Ng_vib(3,4) = 1
     Ng_vib(3,5) = 1
     Ng_vib(3,6) = 1
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
     Ng_vib(6,2) = 1
     Ng_vib(6,3) = 1
     Ng_vib(6,4) = 3
     Ng_vib(6,5) = 3
     Ng_vib(6,6) = 9
     Ng_rot(1,1) = 1
     Ng_rot(1,3) = 2
     Ng_rot(2,2) = 1
     Ng_rot(3,1) = 2
     Ng_rot(3,3) = 3
     Ng_cor(2,1) = 1
     Ng_cor(2,2) = 1
     Ng_cor(2,3) = 1
     Ng_cor(3,2) = 1
     Ng_cor(4,1) = 2
     Ng_cor(4,2) = 2
     Ng_cor(4,3) = 3
     Ng_cor(5,2) = 2
     Ng_cor(6,1) = 4
     Ng_cor(6,2) = 2
     Ng_cor(6,3) = 6
     Npseudo = 12
     !
     !Npseudo = 8
     !
     ig_vib(1,1,1,:) = (/0,0,0,0,0,0/)
     ig_vib(1,2,1,:) = (/0,0,0,2,0,0/)
     ig_vib(1,3,1,:) = (/0,0,0,0,2,0/)
     ig_vib(1,4,1,:) = (/0,1,0,1,0,0/)
     ig_vib(1,5,1,:) = (/0,0,1,0,1,0/)
     ig_vib(2,1,1,:) = (/0,0,0,2,0,0/)
     ig_vib(2,2,1,:) = (/0,0,0,0,0,0/)
     ig_vib(2,4,1,:) = (/1,0,0,1,0,0/)
     ig_vib(2,5,1,:) = (/1,0,0,1,0,2/)
     ig_vib(2,6,1,:) = (/1,0,0,1,5,1/)
     ig_vib(3,1,1,:) = (/0,0,0,0,2,0/)
     ig_vib(3,3,1,:) = (/0,0,0,0,0,0/)
     ig_vib(3,4,1,:) = (/1,0,0,0,1,2/)
     ig_vib(3,5,1,:) = (/1,0,0,0,1,0/)
     ig_vib(3,6,1,:) = (/1,0,0,5,1,1/)
     ig_vib(4,1,1,:) = (/0,1,0,1,0,0/)
     ig_vib(4,2,1,:) = (/1,0,0,1,0,0/)
     ig_vib(4,3,1,:) = (/1,0,0,0,1,2/)
     ig_vib(4,4,1,:) = (/0,2,0,0,0,0/)
     ig_vib(4,4,2,:) = (/2,0,0,0,0,0/)
     ig_vib(4,4,3,:) = (/1,1,0,2,0,0/)
     ig_vib(4,5,1,:) = (/2,0,0,0,0,2/)
     ig_vib(4,5,2,:) = (/1,0,1,0,2,2/)
     ig_vib(4,5,3,:) = (/1,1,0,2,0,2/)
     ig_vib(4,6,1,:) = (/1,0,1,0,3,1/)
     ig_vib(4,6,2,:) = (/2,0,0,0,5,1/)
     ig_vib(4,6,3,:) = (/1,1,0,2,5,1/)
     ig_vib(5,1,1,:) = (/0,0,1,0,1,0/)
     ig_vib(5,2,1,:) = (/1,0,0,1,0,2/)
     ig_vib(5,3,1,:) = (/1,0,0,0,1,0/)
     ig_vib(5,4,1,:) = (/2,0,0,0,0,2/)
     ig_vib(5,4,2,:) = (/1,0,1,0,2,2/)
     ig_vib(5,4,3,:) = (/1,1,0,2,0,2/)
     ig_vib(5,5,1,:) = (/0,0,2,0,0,0/)
     ig_vib(5,5,2,:) = (/2,0,0,0,0,0/)
     ig_vib(5,5,3,:) = (/1,0,1,0,2,0/)
     ig_vib(5,6,1,:) = (/1,1,0,3,0,1/)
     ig_vib(5,6,2,:) = (/2,0,0,5,0,1/)
     ig_vib(5,6,3,:) = (/1,0,1,5,2,1/)
     ig_vib(6,2,1,:) = (/1,0,0,1,5,1/)
     ig_vib(6,3,1,:) = (/1,0,0,5,1,1/)
     ig_vib(6,4,1,:) = (/1,0,1,0,3,1/)
     ig_vib(6,4,2,:) = (/2,0,0,0,5,1/)
     ig_vib(6,4,3,:) = (/1,1,0,2,5,1/)
     ig_vib(6,5,1,:) = (/1,1,0,3,0,1/)
     ig_vib(6,5,2,:) = (/2,0,0,5,0,1/)
     ig_vib(6,5,3,:) = (/1,0,1,5,2,1/)
     ig_vib(6,6,1,:) = (/0,0,2,0,4,0/)
     ig_vib(6,6,2,:) = (/0,2,0,4,0,0/)
     ig_vib(6,6,3,:) = (/2,0,0,0,6,0/)
     ig_vib(6,6,4,:) = (/2,0,0,6,0,0/)
     ig_vib(6,6,5,:) = (/1,0,1,0,7,0/)
     ig_vib(6,6,6,:) = (/1,1,0,7,0,0/)
     ig_vib(6,6,7,:) = (/1,0,1,5,3,2/)
     ig_vib(6,6,8,:) = (/1,1,0,3,5,2/)
     ig_vib(6,6,9,:) = (/2,0,0,5,5,2/)
     ig_rot(1,1,1,:) = (/2,0,0,0,0,0/)
     ig_rot(1,3,1,:) = (/1,0,1,0,3,0/)
     ig_rot(1,3,2,:) = (/2,0,0,0,5,0/)
     ig_rot(2,2,1,:) = (/2,0,0,0,0,0/)
     ig_rot(3,1,1,:) = (/1,0,1,0,3,0/)
     ig_rot(3,1,2,:) = (/2,0,0,0,5,0/)
     ig_rot(3,3,1,:) = (/0,0,2,0,4,0/)
     ig_rot(3,3,2,:) = (/2,0,0,0,6,0/)
     ig_rot(3,3,3,:) = (/1,0,1,0,7,0/)
     ig_cor(2,1,1,:) = (/1,0,0,1,0,1/)
     ig_cor(2,2,1,:) = (/1,0,0,1,0,2/)
     ig_cor(2,3,1,:) = (/1,0,0,1,5,1/)
     ig_cor(3,2,1,:) = (/1,0,0,0,1,0/)
     ig_cor(4,1,1,:) = (/2,0,0,0,0,1/)
     ig_cor(4,1,2,:) = (/1,1,0,2,0,1/)
     ig_cor(4,2,1,:) = (/2,0,0,0,0,2/)
     ig_cor(4,2,2,:) = (/1,1,0,2,0,2/)
     ig_cor(4,3,1,:) = (/1,0,1,0,3,1/)
     ig_cor(4,3,2,:) = (/2,0,0,0,5,1/)
     ig_cor(4,3,3,:) = (/1,1,0,2,5,1/)
     ig_cor(5,2,1,:) = (/2,0,0,0,0,0/)
     ig_cor(5,2,2,:) = (/1,0,1,0,2,0/)
     ig_cor(6,1,1,:) = (/1,0,1,0,3,0/)
     ig_cor(6,1,2,:) = (/1,1,0,3,0,2/)
     ig_cor(6,1,3,:) = (/2,0,0,0,5,0/)
     ig_cor(6,1,4,:) = (/2,0,0,5,0,2/)
     ig_cor(6,2,1,:) = (/1,1,0,3,0,1/)
     ig_cor(6,2,2,:) = (/2,0,0,5,0,1/)
     ig_cor(6,3,1,:) = (/0,0,2,0,4,0/)
     ig_cor(6,3,2,:) = (/2,0,0,0,6,0/)
     ig_cor(6,3,3,:) = (/1,0,1,0,7,0/)
     ig_cor(6,3,4,:) = (/1,0,1,5,3,2/)
     ig_cor(6,3,5,:) = (/1,1,0,3,5,2/)
     ig_cor(6,3,6,:) = (/2,0,0,5,5,2/)
     !
     ipseudo(1,:) = (/0,0,2,1,0,0/)
     ipseudo(2,:) = (/2,0,0,1,0,0/)
     ipseudo(3,:) = (/1,0,1,1,2,0/)
     ipseudo(4,:) = (/0,0,2,1,4,0/)
     ipseudo(5,:) = (/1,0,1,2,1,2/)
     ipseudo(6,:) = (/2,0,0,1,4,0/)
     ipseudo(7,:) = (/1,0,1,2,3,2/)
     ipseudo(8,:) = (/1,1,0,0,5,2/)
     ipseudo(9,:) = (/1,0,1,1,7,0/)
     ipseudo(10,:) = (/1,1,0,8,0,0/)
     ipseudo(11,:) = (/2,0,0,2,5,2/)
     ipseudo(12,:) = (/1,1,0,8,5,2/)
     !
     g_vib(1,1,1) =  (mX1+mX2)/mX1/mX2
     g_vib(1,2,1) =  1/mX1
     g_vib(1,3,1) =  1/mX2
     g_vib(1,4,1) =  -1./mX1
     g_vib(1,5,1) =  -1./mX2
     g_vib(2,1,1) =  1/mX1
     g_vib(2,2,1) =  (mX1+mY1)/mX1/mY1
     g_vib(2,4,1) =  -1./mX1
     g_vib(2,5,1) =  1/mX1
     g_vib(2,6,1) =  -1./mX1
     g_vib(3,1,1) =  1/mX2
     g_vib(3,3,1) =  (mX2+mY2)/mX2/mY2
     g_vib(3,4,1) =  1/mX2
     g_vib(3,5,1) =  -1./mX2
     g_vib(3,6,1) =  -1./mX2
     g_vib(4,1,1) =  -1./mX1
     g_vib(4,2,1) =  -1./mX1
     g_vib(4,3,1) =  1/mX2
     g_vib(4,4,1) =  (mX1+mY1)/mX1/mY1
     g_vib(4,4,2) =  (mX1+mX2)/mX1/mX2
     g_vib(4,4,3) =  -2./mX1
     g_vib(4,5,1) =  -1.*(mX1+mX2)/mX1/mX2
     g_vib(4,5,2) =  1/mX2
     g_vib(4,5,3) =  1/mX1
     g_vib(4,6,1) =  -1./mX2
     g_vib(4,6,2) =  (mX1+mX2)/mX1/mX2
     g_vib(4,6,3) =  -1./mX1
     g_vib(5,1,1) =  -1./mX2
     g_vib(5,2,1) =  1/mX1
     g_vib(5,3,1) =  -1./mX2
     g_vib(5,4,1) =  -1.*(mX1+mX2)/mX1/mX2
     g_vib(5,4,2) =  1/mX2
     g_vib(5,4,3) =  1/mX1
     g_vib(5,5,1) =  (mX2+mY2)/mX2/mY2
     g_vib(5,5,2) =  (mX1+mX2)/mX2/mX1
     g_vib(5,5,3) =  -2./mX2
     g_vib(5,6,1) =  -1./mX1
     g_vib(5,6,2) =  (mX1+mX2)/mX1/mX2
     g_vib(5,6,3) =  -1./mX2
     g_vib(6,2,1) =  -1./mX1
     g_vib(6,3,1) =  -1./mX2
     g_vib(6,4,1) =  -1./mX2
     g_vib(6,4,2) =  (mX1+mX2)/mX1/mX2
     g_vib(6,4,3) =  -1./mX1
     g_vib(6,5,1) =  -1./mX1
     g_vib(6,5,2) =  (mX1+mX2)/mX1/mX2
     g_vib(6,5,3) =  -1./mX2
     g_vib(6,6,1) =  (mX2+mY2)/mX2/mY2
     g_vib(6,6,2) =  (mX1+mY1)/mX1/mY1
     g_vib(6,6,3) =  (mX1+mX2)/mX1/mX2
     g_vib(6,6,4) =  (mX1+mX2)/mX1/mX2
     g_vib(6,6,5) =  -2./mX2
     g_vib(6,6,6) =  -2./mX1
     g_vib(6,6,7) =  -2./mX2
     g_vib(6,6,8) =  -2./mX1
     g_vib(6,6,9) =  2.*(mX1+mX2)/mX2/mX1
     g_rot(1,1,1) =  (mX1+mX2)/mX1/mX2
     g_rot(1,3,1) =  1/mX2
     g_rot(1,3,2) =  -1.*(mX1+mX2)/mX2/mX1
     g_rot(2,2,1) =  (mX1+mX2)/mX1/mX2
     g_rot(3,1,1) =  1/mX2
     g_rot(3,1,2) =  -1.*(mX1+mX2)/mX2/mX1
     g_rot(3,3,1) =  (mX2+mY2)/mX2/mY2
     g_rot(3,3,2) =  (mX1+mX2)/mX1/mX2
     g_rot(3,3,3) =  -2./mX2
     g_cor(2,1,1) =  -1./mX1
     g_cor(2,2,1) =  1/mX1
     g_cor(2,3,1) =  1/mX1
     g_cor(3,2,1) =  -1./mX2
     g_cor(4,1,1) =  (mX1+mX2)/mX1/mX2
     g_cor(4,1,2) =  -1./mX1
     g_cor(4,2,1) =  -1.*(mX1+mX2)/mX1/mX2
     g_cor(4,2,2) =  1/mX1
     g_cor(4,3,1) =  1/mX2
     g_cor(4,3,2) =  -1.*(mX1+mX2)/mX2/mX1
     g_cor(4,3,3) =  1/mX1
     g_cor(5,2,1) =  (mX1+mX2)/mX1/mX2
     g_cor(5,2,2) =  -1./mX2
     g_cor(6,1,1) =  -1./mX2
     g_cor(6,1,2) =  -1./mX1
     g_cor(6,1,3) =  (mX1+mX2)/mX1/mX2
     g_cor(6,1,4) =  (mX1+mX2)/mX1/mX2
     g_cor(6,2,1) =  -1./mX1
     g_cor(6,2,2) =  (mX1+mX2)/mX1/mX2
     g_cor(6,3,1) =  -1.*(mX2+mY2)/mX2/mY2
     g_cor(6,3,2) =  -1.*(mX1+mX2)/mX1/mX2
     g_cor(6,3,3) =  2./mX2
     g_cor(6,3,4) =  1/mX2
     g_cor(6,3,5) =  1/mX1
     g_cor(6,3,6) =  -1.*(mX1+mX2)/mX1/mX2
     pseudo(1) =  -.1250*(mX2+mY2)/mX2/mY2
     pseudo(2) =  -.1250*(mX1+mX2)/mX2/mX1
     pseudo(3) =  -.2500/mX2
     pseudo(4) =  -.1250*(mX2+mY2)/mX2/mY2
     pseudo(5) =  -.2500/mX2
     pseudo(6) =  -.1250*(mX1+mX2)/mX1/mX2
     pseudo(7) =  -.2500/mX2
     pseudo(8) =  -.5000/mX1
     pseudo(9) =  .2500/mX2
     pseudo(10) =  -.5000/mX1
     pseudo(11) =  .2500*(mX1+mX2)/mX1/mX2
     pseudo(12) =  .2500/mX1
     !
     call read_basic_function_constructor(nlines,constructor)
     !
   end subroutine  MLkinetic_abcd_EKE_z_alpha1_singular


end module kin_abcd
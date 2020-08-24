!
!  This unit makes a bridge between different molecular types and the main program
!
module molecules
  use accuracy
  use lapack
  use moltype
  use mol_xy
  use mol_xy2
  use mol_xy3
  use mol_xy4, only  : ML_coordinate_transform_XY4, ML_symmetry_transformation_XY4, ML_b0_XY4, ML_rotsymmetry_XY4
  use mol_zxy2, only : ML_coordinate_transform_ZXY2, ML_symmetry_transformation_ZXY2, ML_b0_ZXY2, ML_rotsymmetry_ZXY2,&
                       ML_coordinate_transform_SOHF,ML_b0_SOHF
  use mol_zxy3, only : ML_coordinate_transform_ZXY3, ML_symmetry_transformation_ZXY3, ML_b0_ZXY3, ML_rotsymmetry_ZXY3
  use mol_ch3oh
  use mol_abcd
  use mol_c2h4
  use mol_c2h6
  use mol_c3h6
  
  use pot_xy2
  use pot_xy3
  use pot_abcd
  use pot_zxy2
  use pot_zxy3
  use pot_xy4
  use pot_ch3oh
  use pot_c2h4, only : ML_dipole_c2h4_4m_dummy,MLpoten_c2h4_88, MLpoten_c2h4_lee,MLpoten_c2h4_886666
  use pot_c2h6, only : MLpoten_c2h6_88,MLpoten_c2h6_88_cos3tau,MLpoten_c2h6_88_cos3tau_142536,&
                       MLpoten_c2h6_88_cos3tau_sym,MLpoten_c2h6_Duncan,&
                       MLpoten_c2h6_88_cos3tau_G36,ML_alpha_C2H6_zero_order
  use pot_c3h6, only : MLpoten_c3h6_harmtest,MLpoten_c3h6_sym_II
  !
  use prop_xy2,      only : prop_xy2_qmom_sym,MLdipole_h2o_lpt2011
  use prop_xy2_quad, only : prop_xy2_qmom_bisect_frame,TEST_xy2_qmom_bisect_frame
  use prop_xy2_spinrot, only : prop_xy2_spin_rotation_bisector, prop_xy2_spin_rotation_bisector_nonlin, &
                               TEST_prop_xy2_spin_rotation_bisector_nonlin, prop_xy2_gtensor_bisector,&
                               prop_xy2_gtens_nuclear_bisector, prop_xy2_grot_electronic_bisector, &
                               prop_xy2_gcor_electronic_bisector
  use prop_xy2_spinspin, only : prop_xy2_spinspin_dipoleYY
  !
  use kin_xy2, only  : MLkinetic_xy2_bisect_EKE,MLkinetic_xyz_bisect_EKE,MLkinetic_xy2_bisect_EKE_sinrho,&
                       MLkinetic_xy2_Radau_bisect_EKE
  !
  use pot_user, only : MLdipole,MLpoten,ML_MEP
  !
  use symmetry,     only : sym
  !
  implicit none
  public MLequilibrium_xyz, &
         MLcoord_direct,MLcoord_invert,MLcoord_firstderiv, & 
         MLcoordinate_transform_func,ML_check_steps4coordinvert, &
         MLequilibrium_xyz_1d,diff_2d_4points,ML_diffs,MLsymmetry_transform_func
  public MLpotenfunc,MLrotsymmetry_func,ratint,polint,MLinvmat,diff_2d_4points_ark,diff_3d_6points,polintark,&
         MLextF_func,MLpotentialfunc,ML_MEPfunc,MLinvmatark,MLratintark,MLrotsymmetry_generate_CII,&
         MLrotsymmetry_generate,MLkineticfunc
  public MLdefine_potenfunc,MLcoordinate_transform_func_define,MLextF_func_define,MLdefine_kinetic_subroutine
         !
  public MOrepres_arkT,ddlmn_conj,dlmn,Phi_rot,calc_phirot
   !
  private
   !
   integer(ik), parameter :: verbose     = 4                          ! Verbosity level
   !
   type MOrepres_arkT
      real(ark),pointer     :: repres(:,:,:)  
   end type MOrepres_arkT
   !
   !procedure (MLtemplate_poten),pointer :: MLpotenfunc => null()
   !
   procedure (MLtemplate_potential),pointer :: MLpotentialfunc => null()
   procedure (MLtemplate_coord_transform),pointer :: MLcoordinate_transform_func => null()
   procedure (MLtemplate_b0),pointer :: MLequilibrium_xyz  => null()
   procedure (MLtemplate_extF),pointer :: MLextF_func => null()
   procedure (MLtemplate_symmetry_transformation),pointer :: MLsymmetry_transform_func => null()
   procedure (MLtemplate_rotsymmetry),pointer :: MLrotsymmetry_func => null()
   procedure (MLtemplate_kinetic),pointer :: MLkineticfunc => null()
   !
  contains

  !
  ! Defining potential energy function 
  !

  !
  ! Defining potential energy function 
  !
  function MLpotenfunc(local,xyz)  result(f)

   real(ark),intent(in) ::  local(:)
   real(ark),intent(in) ::  xyz(:,:)
   real(ark)            ::  f
     !
     f=MLpotentialfunc(molec%ncoords,molec%natoms,local,xyz,molec%force)
     !
  end function MLpotenfunc
  !
  !
  subroutine MLdefine_potenfunc
   !
   if (verbose>=6) write(out,"(/'MLdefine_potenfunc/start')") 

    select case(trim(molec%potentype))
         !
    case default
         !
         write (out,"('MLdefine_potenfunc: potential type ',a,' unknown')") trim(molec%potentype)
         stop 'MLdefine_potenfunc - bad potential'
         !
    case('POTEN_XY3_MLT') 
         !
         MLpotentialfunc => MLpoten_xy3_mlt
         !
    case('POTEN_XY3_MORBID_10','POTEN_XY3_MORBID_10_MEP') 
         !
         MLpotentialfunc => MLpoten_xy3_morbid_10
         !
    case('POTEN_XY3_MORBID_11') 
         !
         MLpotentialfunc => MLpoten_xy3_morbid_11
         !
    case('POTEN_XY3_MORBID_MORPHING') 
         !
         MLpotentialfunc => MLpoten_xy3_morbid_morphing
         !
    case('POTEN_XY2_SCHWENKE') 
         !
         MLpotentialfunc => MLpoten_xy2_schwenke
         !
    case('POTEN_XY3_MORBID_45760') 
         !
         MLpotentialfunc => MLpoten_xy3_morbid_45760 
         !
    case('POTEN_XY3_MORBID_DELTA') 
         !
         MLpotentialfunc => MLpoten_xy3_morbid_delta
         !
    case('POTEN_XY3_HSL')
         !
         MLpotentialfunc => MLpoten_xy3_HSL
         !
    case('POTEN_XY3_HANDY')
         !
         MLpotentialfunc => MLpoten_xy3_handy
         !
    case('POTEN_XY3_D3H')
         !
         MLpotentialfunc => MLpoten_xy3_d3h
         !
    case('POTEN_XY3_MLT_II')
         !
         MLpotentialfunc => MLpoten_xy3_mlt_II
         !
    case('POTEN_XY3_SEARS')
         !
         MLpotentialfunc => MLpoten_xy3_sears
         !
    case('POTEN_ZXY2_MLT')
         !
         MLpotentialfunc => MLpoten_zxy2_mlt
         !
    case('POTEN_ZXY2_ANDREY_01')
         !
         MLpotentialfunc => MLpoten_zxy2_andrey_01
         !
    case('POTEN_ZXY2_MEP_R_ALPHA_RHO_POWERS')
         !
         MLpotentialfunc => MLpoten_zxy2_mep_r_alpha_rho_powers
         !
    case('POTEN_ZXY2_MEP_R_ALPHA_RHO_COEFF')
         !
         MLpotentialfunc => MLpoten_zxy2_andrey_coeff
         !
    case('POTEN_H2CS_TZ1')
         !
         MLpotentialfunc => MLpoten_h2cs_tz_damp1
         !
    case('POTEN_H2CS_DAMP')
         !
         MLpotentialfunc => MLpoten_h2cs_damp
         !
    case('POTEN_H2CS_DAMP_SCALING')
         !
         MLpotentialfunc => MLpoten_h2cs_damp_scaling
         !
    case('POTEN_ABCD') 
         !
         MLpotentialfunc => MLpoten_hsoh
         !
    case('POTEN_ABCD_REF') 
         !
         MLpotentialfunc => MLpoten_hsoh_ref
         !
     case('POTEN_H2O2_KOPUT') 
         !
         MLpotentialfunc => MLpoten_h2o2_koput
         !
     case('POTEN_H2O2_KOPUT_MORSE') 
         !
         MLpotentialfunc => MLpoten_h2o2_koput_morse
         !
     case('POTEN_C2H2_MORSE') 
         !
         MLpotentialfunc => MLpoten_c2h2_morse
         !
     case('POTEN_C2H2_7') 
         !
         MLpotentialfunc => MLpoten_c2h2_7
         !
     case('POTEN_C2H2_7_Q1Q2Q3Q4') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_q1q2q3q4
         !
     case('POTEN_C2H2_7_Q2Q1Q4Q3_LINEARIZED') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_q2q1q4q3_linearized
         !
     case('POTEN_C2H2_7_Q2Q1Q4Q3_LINEAR_MORPHED') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_q2q1q4q3_linearized_morphing
         !
     case('POTEN_C2H2_7_Q2Q1Q4Q3') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_q2q1q4q3         
         !
     case('POTEN_C2H2_7_415') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_415
         !
     case('POTEN_C2H2_7_XYZ') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_xyz
         !
     case('POTEN_C2H2_7_R_RR_NNNN') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_r_rr_nnnn
         !
     case('POTEN_C2H2_7_R_RR_XY') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_r_rr_xy
         !
     case('POTEN_C2H2_7_R_ZZ_NNNN') 
         !
         MLpotentialfunc => MLpoten_c2h2_7_r_zz_nnnn
         !
     case('POTEN_C2H2_STREYMILLS') 
         !
         MLpotentialfunc =>  MLpoten_c2h2_streymills
         !
     case('POTEN_H2O2_KOPUT_UNIQUE') 
         !
         MLpotentialfunc => MLpoten_h2o2_koput_unique
         !
     case('POTEN_C2H2_KATY') 
         !
         MLpotentialfunc => MLpoten_v_c2h2_katy
         !
     case('POTEN_C2H2_MLT') 
         !
         MLpotentialfunc => MLpoten_v_c2h2_mlt
         !
     case('POTEN_P2H2_MORSE_COS') 
         !
         MLpotentialfunc => MLpoten_p2h2_morse_cos
         !
    case('POTEN_XY2_MORBID') 
         !
         MLpotentialfunc => MLpoten_xy2_morbid
         !
    case('POTEN_XY2_DMBE') 
         !
         MLpotentialfunc => MLpoten_xy2_dmbe
         !
    case('POTEN_H2O_BUBUKINA') 
         !
         MLpotentialfunc => MLpoten_xy2_bubukina
         !
    case('POTEN_XY2_TYUTEREV') 
         !
         MLpotentialfunc => MLpoten_xy2_tyuterev
         !
    case('POTEN_XY2_TYUTEREV_ALPHA') 
         !
         MLpotentialfunc => MLpoten_xy2_tyuterev_alpha
         !
    case('POTEN_XY2_MORSE_COS') 
         !
         MLpotentialfunc => MLpoten_xy2_morse_cos
         !
    case('POTEN_XYZ_KOPUT') 
         !
         MLpotentialfunc => MLpoten_xyz_Koput
         !
    case('POTEN_XYZ_TYUTEREV') 
         !
         MLpotentialfunc => MLpoten_xyz_tyuterev
         !
    case('POTEN_SO2_PES_8D') 
         !
         MLpotentialfunc => MLpoten_SO2_pes_8d
         !
    case('POTEN_XY2_TYUTEREV_DAMP') 
         !
         MLpotentialfunc => MLpoten_xy2_tyuterev_damp
         !
    case('POTEN_SO2_DAMP') 
         !
         MLpotentialfunc => MLpoten_so2_damp
         !
    case('POTEN_CO2_AMES1') 
         !
         MLpotentialfunc => MLpoten_co2_ames1
         !
    case('POTEN_SO2_AMES1') 
         !
         MLpotentialfunc => MLpoten_so2_ames1
         !
    case('POTEN_H2O_TENNYSON') 
         !
         MLpotentialfunc => MLpoten_h2o_tennyson
         !
    case('POTEN_C3_MLADENOVIC') 
         !
         MLpotentialfunc => MLpoten_c3_mladenovic
         !
    case('POTEN_C3_R_THETA') 
         !
         MLpotentialfunc => MLpoten_c3_R_theta
         !
    case('POTEN_H2S_DVR3D') 
         !
         MLpotentialfunc => MLpoten_h2s_dvr3d
         !
    case('POTEN_XY2_HALONEN_I') 
         !
         MLpotentialfunc => MLpoten_xy2_halonen_I
         !
    case('POTEN_XY2_MLT_CO2') 
         !
         MLpotentialfunc => MLpoten_xy2_mlt_co2
         !
    case('POTEN_XY_MORSE','POTEN_XY_DUNHAM','POTEN_XY_SPF','POTEN_XY_SPF_H2') 
         !
         MLpotentialfunc => MLpoten_xy_gener
         !
    case('POTEN_XY4_ZZZ') 
         !
         MLpotentialfunc => MLpoten_xy4_ZZZ
         !
    case('POTEN_USER','USER','GENERAL','POTEN') 
         !
         MLpotentialfunc => MLpoten
         !
    case('POTEN_XY4_NIKITIN') 
         !
         MLpotentialfunc => MLpoten_xy4_Nikitin
         !
    case('POTEN_XY4_BOWMAN2000') 
         !
         MLpotentialfunc => MLpoten_xy4_Bowman2000
         !
    case('POTEN_CH3OH_REF') 
         !
         MLpotentialfunc => MLpoten_ch3oh_ref
         !
    case('POTEN_SOHF') 
         !
         MLpotentialfunc => MLpoten_sohf
         !
    case('POTEN_OH3P_MEP') 
         !
         MLpotentialfunc => MLpoten_oh3p_mep
         !
    case('POTEN_C2H4_88') 
         !
         MLpotentialfunc => MLpoten_c2h4_88
         !
    case('POTEN_C2H4_886666') 
         !
         MLpotentialfunc => MLpoten_c2h4_886666
         !
    case('POTEN_C2H6_88') 
         !
         MLpotentialfunc => MLpoten_c2h6_88
         !
    case('POTEN_C2H6_88_COS3TAU') 
         !
         MLpotentialfunc => MLpoten_c2h6_88_cos3tau
         !
    case('POTEN_C2H6_88_COS3TAU_142536') 
         !
         MLpotentialfunc => MLpoten_c2h6_88_cos3tau_142536
         !
    case('POTEN_C2H6_88_COS3TAU_SYM') 
         !
         MLpotentialfunc => MLpoten_c2h6_88_cos3tau_sym
         !
    case('POTEN_C2H6_88_COS3TAU_G36') 
         !
         MLpotentialfunc => MLpoten_c2h6_88_cos3tau_G36
         !
    case('POTEN_C2H6_DUNCAN') 
         !
         MLpotentialfunc => MLpoten_c2h6_Duncan
         !
    case('POTEN_ZXY3_SYM') 
         !
         MLpotentialfunc => MLpoten_zxy3_sym
         !
    case('POTEN_ZXY3_NIKITIN') 
         !
         MLpotentialfunc => MLpoten_zxy3_Nikitin
         !
    case('POTEN_C3H6_HARMTEST') 
         !
         MLpotentialfunc => MLpoten_c3h6_harmtest
         !
    case('POTEN_C3H6_SYM') 
         !
         MLpotentialfunc => MLpoten_c3h6_sym_II
         !
    end select
    !
   if (verbose>=6) write(out,"('MLdefine_potenfunc/end')") 
 
end subroutine MLdefine_potenfunc


  !
  subroutine MLdefine_kinetic_subroutine
   !
   if (verbose>=6) write(out,"(/'MLdefine_kinetic_subroutine/start')") 
    !
    select case(trim(molec%kinetic_type))
    case default
         !
         write (out,"('MLdefine_kinetic_subroutine: kinetic type ',a,' unknown')") trim(molec%kinetic_type)
         stop 'MLdefine_kinetic_subroutine - bad kinetic'
         !
    case('KINETIC_XY2_EKE_BISECT') 
         !
         MLkineticfunc => MLkinetic_xy2_bisect_EKE
         !
    case('KINETIC_XY2_EKE_RADAU_BISECT') 
         !
         MLkineticfunc => MLkinetic_xy2_Radau_bisect_EKE
         !
    case('KINETIC_XY2_EKE_BISECT_SINRHO') 
         !
         MLkineticfunc => MLkinetic_xy2_bisect_EKE_sinrho
         !
    case('KINETIC_XYZ_EKE_BISECT') 
         !
         MLkineticfunc => MLkinetic_xyz_bisect_EKE
         !
    case('GENERAL') 
         !
         MLkineticfunc => MLkinetic_dummy
         !
    end select
    !
   if (verbose>=6) write(out,"('MLdefine_kinetic_subroutine/end')") 
 
  end subroutine MLdefine_kinetic_subroutine


  ! A dummy kinetic energy function 
  !
  subroutine MLkinetic_dummy(nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
   !
   integer(ik),intent(in) ::  nmodes,Nterms
   real(ark),intent(in)   ::  rho
   real(ark),intent(out)  ::  g_vib(nmodes,nmodes,Nterms),g_rot(3,3,Nterms),g_cor(nmodes,3,Nterms),pseudo(Nterms)
     !
     write(out,"('MLkinetic_MLkinetic_dummy: If you are here you use LOCAL but KINETIC is undefined')")
     stop 'MLkinetic_MLkinetic_dummy: KINETIC is undefined but it should be to be used with LOCAL'
     !
     g_vib=0
     g_rot=0
     g_cor=0
     pseudo=0
     !
   end subroutine  MLkinetic_dummy

  !
  ! Defining MEP function 
  !
  !
  function ML_MEPfunc(nsize,src)  result(dst)

   integer(ik),intent(in) :: nsize
   real(ark),intent(in) ::  src
   real(ark)            ::  dst(1:nsize)

   if (verbose>=6) write(out,"(/'ML_MEPfunc/start')") 

    select case(trim(molec%meptype))
    case default
         write (out,"('ML_MEPfunc: MEP type ',a,' unknown')") trim(molec%meptype)
         stop 'ML_MEPfunc - bad MEP type'

    !case('MEP_XY2_R12_R') 
    !     !
    !     dst(1) =  ML_MEP_xy2_R12_R(src)
    !     !
    !case('MEP_XY2_R12_ALPHA') 
    !     !
    !     dst(1) =  ML_MEP_xy2_R12_alpha(src)
    !     !
    !case('MEP_ABCD_TAU-REF') 
    !     !
    !     dst(1:6) =  ML_MEP_ABCD_tau_ref(src)
    !     !
    !case('MEP_NH3')
    !     !
    !     dst(1:nsize) =  ML_mep_nh3(src)
    !     !
    !case('MEP_ZXY2_R_RHO')
    !     !
    !     dst(1:nsize) =  ML_MEP_ZXY2_R_RHO(src)
    !     !
    !case('MEP_ZXY2_R_COEFF')
    !     !
    !     dst(1:nsize) =  ML_MEP_zxy2_rho_coeff(src)
    !     !
    case('MEP_USER','GENERIC')
         !
         dst(1:nsize) =  ML_MEP(nsize,src)
         !
    end select 

   if (verbose>=6) write(out,"('ML_MEPfunc/end')") 
 
end function ML_MEPfunc

  !
  ! Defining external field function 
  !
  !
  recursive subroutine MLextF_func_define
   !
   if (verbose>=6) write(out,"(/'MLextF_func_define/start')") 
    !
    select case (trim(extF%ftype))
    !
    case default
       !
       write(out, '(/2(1x, a))') 'MLextF_func error: unknown type of extF:', trim(extF%ftype)
       stop 'MLextF_func error: unknown type of extF'
       !
    ! pq space-fixed frame of P.Jensen for xy2-type molecule (J.Mol.Spectr.132 (1988) 429)
    !
    case('XY2_PQ')
        !
        MLextF_func => MLdms2pqr_xy2
        !
    case('XY2_PQ_COEFF')
        !
        MLextF_func => MLdms2pqr_xy2_coeff
        !
    case('XY2_PQ_LINEAR')
        !
        MLextF_func => MLdms2pqr_xy2_linear
        !
    case('DIPOLE_SO2_AMES1')
        !
        MLextF_func => MLdipole_so2_ames1
        !
    case('DIPOLE_AMES1')
        !
        MLextF_func => MLdipole_ames1
        !
    case('DIPOLE_XY2_LORENZO')
        !
        MLextF_func => MLdipole_xy2_lorenzo
        !
    case('DIPOLE_H2O_LPT2011')
        !
        MLextF_func => MLdipole_h2o_lpt2011
        !
    case('DIPOLE_PQR_XYZ')
        !
        MLextF_func => MLdms2pqr_xyz_coeff
        !
    case('DIPOLE_BISECT_S1S2T_XYZ')
        !
        MLextF_func => MLdipole_bisect_s1s2theta_xy2
        !
    case('XY2_QMOM_SYM')
        !
        MLextF_func => prop_xy2_qmom_sym
        !
    case('XY2_QMOM_BISECT_FRAME')
        !
        MLextF_func => prop_xy2_qmom_bisect_frame
        !
    case('TEST_XY2_QMOM_BISECT_FRAME')
        !
        MLextF_func => TEST_xy2_qmom_bisect_frame
        !
    case('XY2_SR-BISECT-NONLIN')
        !
        MLextF_func => prop_xy2_spin_rotation_bisector_nonlin
        !
    case('TEST_XY2_SR-BISECT-NONLIN')
        !
        MLextF_func => TEST_prop_xy2_spin_rotation_bisector_nonlin
        !
    case('XY2_SR-BISECT')
        !
        MLextF_func =>  prop_xy2_spin_rotation_bisector
        !
    case('XY2_SS_DIPOLE_YY')
        !
        MLextF_func =>  prop_xy2_spinspin_dipoleYY
        !
    case('XY2_G-BISECT')
        !
        MLextF_func =>  prop_xy2_gtensor_bisector
        !
    case('XY2_G-ROT-ELEC')
        !
        MLextF_func =>  prop_xy2_grot_electronic_bisector
        !
    case('XY2_G-COR-ELEC')
        !
        MLextF_func =>  prop_xy2_gcor_electronic_bisector
        !
    case('XY2_G-TENS-NUC')
        !
        MLextF_func =>  prop_xy2_gtens_nuclear_bisector
        !
    case('XY3_MB')
        !
        MLextF_func => MLdms2xyz_xy3_mb
        !
    case('XY3_MB4')
        !
        MLextF_func => MLdms2xyz_xy3_mb4
        !
    case('XY3_SYMMB')
       !
       MLextF_func => MLdms2xyz_xy3_symmb
       !
       ! space-fixed frame for abcd-type (hsoh) molecule: ba = +x, bc = +z
       !
    case('ABCD')
       !
       MLextF_func => MLdms2xyz_abcd
       !
       ! Molecular Bond representaion of the DMS of HOOH
       !
    case('HOOH_MB')
       !
       MLextF_func => MLdms_hooh_MB
       !
       ! Molecular Bond representaion of alpha of HOOH
       !
    case('HOOH_ALPHA_MB')
       !
       MLextF_func => MLalpha_hooh_MB
       !
    case('HPPH_MB')
       !
       MLextF_func => MLdms_hpph_MB
       !
       ! Molecular Bond representaion of the DMS of HCCH
       !
    case('HCCH_MB')
       !
       MLextF_func => MLdms_hcch_MB
       !
    case('HCCH_DMS_7D')
       !
       MLextF_func => MLdms_HCCH_7D
       !
    case('HCCH_DMS_7D_7ORDER')
       !
       MLextF_func => MLdms_HCCH_7D_7ORDER
       !
    case('HCCH_DMS_7D_7ORDER_LINEAR')
       !
       MLextF_func => MLdms_HCCH_7D_7ORDER_linear
       !
    case('HCCH_DMS_7D_LOCAL')
       !
       MLextF_func => MLdms_HCCH_7D_LOCAL
       !
    case('COORDINATES')
       !
       MLextF_func => MLextF_coordinates
       !
    case('XY3_NSS_MB')
       !
       MLextF_func => MLspinspin_xy3_mb
       !
    case('ZXY2_SYMADAP')
       !
       MLextF_func => MLdms2xyz_zxy2_symadap_powers
       !
    case('ZXY3_SYM')
       !
       MLextF_func => MLdms2xyz_zxy3_sym
       !
    case('DIPOLE_C2H4_4M') 
       !
       MLextF_func => ML_dipole_c2h4_4m_dummy  ! dummy dipole does not work
       !
    case('ALPHA_C2H6_ZERO') 
       !
       MLextF_func => ML_alpha_C2H6_zero_order  ! alpha polarizablity of a zero order type 
       !
    case('DIPOLE','USER','GENERAL','DIPOLE_USER')
       !
       MLextF_func => MLdipole
       !
    case('POTENTIAL','POTEN')
        !
        MLextF_func => MLextF_potential
        !
    end select
    ! 
   if (verbose>=6) write(out,"('MLextF_func_define/end')") 
 
  end subroutine MLextF_func_define
  !

  recursive subroutine MLextF_coordinates(rank,ncoords,natoms,local,xyz,mu)

   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  mu(rank)
     !
     if (rank/=molec%ncoords) stop 'rank /= ncoords for COORDS in MLextF_func'
     !
     mu  = local
     !
  end subroutine MLextF_coordinates



  recursive subroutine MLextF_potential(rank,ncoords,natoms,local,xyz,mu)

   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  mu(rank)
   !
   integer(ik)            :: i,imu
   real(ark)              :: force(molec%parmax) ! force(rank)
     !
     do i = 1,rank
       !
       if (extF%ifit(1,i)==0) then 
         force(i) = extF%coef(1,i)
       else
         force(i) = 0
       endif
       !
     enddo
     !
     force(rank:molec%parmax) = 0
     !
     mu = 0 
     !
     do imu = 1,rank
       !
       !force = 0 
       !
       if (extF%ifit(1,imu)==0) cycle
       !
       force(imu) = 1.0_ark
       !
       mu(imu) = MLpotentialfunc(ncoords,natoms,local,xyz,force)
       !
       force(imu) = 0.0_ark
       !
     enddo
    !
  end subroutine MLextF_potential



  !
  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates: linear or normal case
  !
  function ML_coordinate_transform_linear(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_linear/start')") 
    !
    if (direct) then 
        dst = src(:) - molec%local_eq(:)
    else
        dst(:) = src(:) +  molec%local_eq(:)
    endif
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_linear/end')") 
    !
    !
  end function ML_coordinate_transform_linear
  !
  function ML_coordinate_transform_dx_rho(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    !
    integer(ik) :: nsrc
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_dx_rho/start')") 
    !
    if (direct) then 
        dst = src(:) - molec%local_eq(:)
        dst(molec%Nmodes) = src(molec%Nmodes)
    else
        dst(:) = src(:) +  molec%local_eq(:)
        dst(molec%Nmodes) = src(molec%Nmodes)
    endif
    !
    if (verbose>=5) write(out,"('ML_coordinate_transform_dx_rho/end')") 
    !
    !
  end function ML_coordinate_transform_dx_rho
  !
  !
  !
  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  subroutine MLcoordinate_transform_func_define
    !
    !
    if (verbose>=5) write(out,"(/'MLcoordinate_transform_func_define/start')") 
    !
    ! Identical coordinate transformation is a general case - works identically for all molecules
    !
    select case(trim(molec%moltype))
    case default
         write (out,"('MLcoordinate_transform_func_define: molecule type ',a,' unknown')") trim(molec%moltype)
         stop 'MLcoordinate_transform_func_define - bad molecule'

    case('ABCD') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_ABCD
         MLequilibrium_xyz => ML_b0_ABCD
         MLsymmetry_transform_func => ML_symmetry_transformation_ABCD
         MLrotsymmetry_func => ML_rotsymmetry_ABCD
         !
    case('C2H4') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_C2H4
         MLequilibrium_xyz => ML_b0_C2H4
         MLsymmetry_transform_func => ML_symmetry_transformation_C2H4
         MLrotsymmetry_func => ML_rotsymmetry_C2H4
         !
    case('C2H6') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_C2H6
         MLequilibrium_xyz => ML_b0_C2H6
         MLsymmetry_transform_func => ML_symmetry_transformation_C2H6
         MLrotsymmetry_func => ML_rotsymmetry_C2H6
         !
    case('C3H6') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_C3H6
         MLequilibrium_xyz => ML_b0_C3H6
         MLsymmetry_transform_func => ML_symmetry_transformation_C3H6
         MLrotsymmetry_func => ML_rotsymmetry_C3H6
         !
    case('XY3') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_XY3
         MLequilibrium_xyz => ML_b0_XY3
         MLsymmetry_transform_func => ML_symmetry_transformation_XY3
         MLrotsymmetry_func => ML_rotsymmetry_XY3
         !
    case('XY4') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_XY4
         MLequilibrium_xyz => ML_b0_XY4
         MLsymmetry_transform_func => ML_symmetry_transformation_XY4
         MLrotsymmetry_func => ML_rotsymmetry_XY4
         !
    case('ZXY2') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_ZXY2
         MLequilibrium_xyz => ML_b0_ZXY2
         MLsymmetry_transform_func => ML_symmetry_transformation_ZXY2
         MLrotsymmetry_func => ML_rotsymmetry_ZXY2
         !
    case('ZXY3') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_ZXY3
         MLequilibrium_xyz => ML_b0_ZXY3
         MLsymmetry_transform_func => ML_symmetry_transformation_ZXY3
         MLrotsymmetry_func => ML_rotsymmetry_ZXY3
         !
    case('XY2') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_XY2
         MLequilibrium_xyz => ML_b0_XY2
         MLsymmetry_transform_func => ML_symmetry_transformation_XY2
         MLrotsymmetry_func => ML_rotsymmetry_XY2
         !
    case('CH3OH') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_ch3oh
         MLequilibrium_xyz => ML_b0_ch3oh
         MLsymmetry_transform_func => MLsymmetry_transform_C
         MLrotsymmetry_func => ML_rotsymmetry_C
         !
    case('SOHF') 
         !
         MLcoordinate_transform_func =>  ML_coordinate_transform_SOHF
         MLequilibrium_xyz => ML_b0_SOHF
         MLsymmetry_transform_func => MLsymmetry_transform_C
         MLrotsymmetry_func => ML_rotsymmetry_C
         !
    case('XY')
         !
         MLequilibrium_xyz => ML_b0_XY
         MLsymmetry_transform_func => MLsymmetry_transform_C
         MLrotsymmetry_func => ML_rotsymmetry_C
         !
    end select
    !
    ! for some special cases we can redefine the coord_transform to a simpler and faster function
    !
    select case(trim(trim(molec%coords_transform)))
       !
    case('LINEAR','HARMONIC','X-XE')
       !
       MLcoordinate_transform_func => ML_coordinate_transform_linear
       !
    case('DX-RHO')
       !
       MLcoordinate_transform_func => ML_coordinate_transform_dx_rho
       !
    case default
       !
       if (trim(molec%moltype)=='XY') then  
          write (out,"('MLcoordinate_transform_func_define: coord. transformation for XY ',a,' unknown')") &
                trim(molec%coords_transform)
          stop 'MLcoordinate_transform_func_define - bad coord. type'
       endif
       !
    end select
    !
    if (verbose>=5) write(out,"('MLcoordinate_transform_func_define/end')") 
    !
    !
  end subroutine MLcoordinate_transform_func_define
  !
  !
  !
  subroutine MLsymmetry_transform_C(ioper,natoms,src,dst)
    !
    integer(ik),intent(in)    :: ioper  ! group operation  
    integer(ik),intent(in)    :: natoms
    real(ark),intent(in)      :: src(natoms)
    real(ark),intent(out)     :: dst(natoms)
    !
    ! trivial case of no symmetry
    !
    select case(trim(molec%symmetry))
    case default
       !
       write (out,"('MLsymmetry_transform_C: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'MLsymmetry_transform_C - bad symm. type'
       !
    case('C','C(M)')
       !
       dst = src(:)
       !
    end select
    !
  end subroutine MLsymmetry_transform_C
  !
  !
  subroutine ML_rotsymmetry_C(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
    integer(ik),intent(out) :: gamma,ideg

    ! trivial case of no symmetry
    !
    select case(trim(molec%symmetry))
    case default
       !
       write (out,"('ML_rotsymmetry_C: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_rotsymmetry_C - bad symm. type'
       !
    case('C','C(M)')
       !
       gamma = 1
       ideg = 1
       !
    end select
    !
  end subroutine ML_rotsymmetry_C
  !
  !
  ! define the symmetry of the rotational basis function 
  !
  subroutine MLrotsymmetry_func_(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau  ! group operation  

    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out,"(/'MLrotsymmetry_func/start')") 
    !
    ! trivial case of no symmetry
    !
    if (trim(molec%symmetry)=='C'.or.trim(molec%symmetry)=='C(M)') then
       !
       gamma = 1
       ideg = 1
       return
       !
    endif
    !
    select case(trim(molec%moltype))
    case default
         write (out,"('MLrotsymmetry_func: molecule type ',a,' unknown')") trim(molec%moltype)
         stop 'MLrotsymmetry_func - bad molecule'
         !
    case('XY2') 
         !
         call  ML_rotsymmetry_XY2(J,K,tau,gamma,ideg)
         !
    case('XY3') 
         !
         call  ML_rotsymmetry_XY3(J,K,tau,gamma,ideg)
         !
    case('XY4') 
         !
         call  ML_rotsymmetry_XY4(J,K,tau,gamma,ideg)
         !
    case('ZXY2') 
         !
         call  ML_rotsymmetry_ZXY2(J,K,tau,gamma,ideg)
         !
    case('ABCD') 
         !
         call  ML_rotsymmetry_abcd(J,K,tau,gamma,ideg)
         !
    end select 
    !
    if (verbose>=5) write(out,"('MLrotsymmetry_func/end')") 
    !
    !
  end subroutine MLrotsymmetry_func_





  !
  ! define the symmetry of the rotational basis function 
  !
  subroutine MLrotsymmetry_generate(J,iverbose,count_index,eigenvects,Ncount,Ntot)
    !
    integer(ik),intent(in)   :: J,iverbose
    integer(ik),intent(out)  :: count_index(2*j+1,2*j+1),Ncount,Ntot(:)
    real(rk),   intent(out)  :: eigenvects(2*j+1,2*j+1)
    real(ark)                :: theta,phi,chi,coeff,vect_t(2*j+1)
    !
    integer(ik)              :: i1,k1,n1,i2,k2,n2,tau1,tau2,i,ioper,sigma,alloc,iroot,icount,ideg,nroots
    integer(ik)              :: izeta,irep,jroot,irep_t
    complex(ark)             :: temp,d
    complex(ark),allocatable :: R(:,:),T(:,:),Q(:,:)
    real(ark),allocatable    :: tmat(:,:,:)

    type(MOrepres_arkT),allocatable   :: irr(:)

    !
    integer(ik)              :: count_degen(2*j+1),Ntotal(1:sym%Nrepresen)
    integer(ik)              :: isym
    !
    if (iverbose>=5) write(out,"(/'MLrotsymmetry_generate/start')") 
    !
    ! trivial case of no symmetry
    !
    if (trim(molec%symmetry)=='C'.or.trim(molec%symmetry)=='C(M)') then
       !
       !gamma = 1
       !ideg = 1
       return
       !
    endif
    !
    ! dimension of the problem
    !
    nroots  = 2*j+1
    !
    ! Obtain transformation-matrix 
    !
    allocate(tmat(sym%Noper,nroots,nroots),R(nroots,nroots),T(nroots,nroots),Q(nroots,nroots),stat=alloc)
    if (alloc/=0) stop 'MLrotsymmetry_generate - out of memory'

    !
    ! R is a transformation to the Wang functions
    !
    R = 0 
    !
    if (mod(j,2)==0) then 
      !
      r(1,1) = cmplx(1.0_ark,0.0_ark)
      !
    else
      !
      r(1,1) = cmplx(0.0_ark,-1.0_ark)
      !
    endif
    !
    i1 = 2
    do k1 = 1,j
       !
       sigma = mod(k1,3)
       !
       R(i1  ,i1  ) = cmplx(1.0_ark/sqrt(2.0_ark)           ,0.0_ark                               )
       R(i1  ,i1+1) = cmplx(0.0_ark                         ,-(-1.0_ark)**sigma/sqrt(2.0_ark)      )
       R(i1+1,i1  ) = cmplx((-1.0_ark)**(j+k1)/sqrt(2.0_ark),0.0_ark                               )
       R(i1+1,i1+1) = cmplx(0.0_ark                         ,(-1.0_ark)**(j+k1+sigma)/sqrt(2.0_ark))
       !
       i1 = i1 + 2
       !
    enddo
    !
    ! transformation to the wang rot. basis 
    !
    !
    do ioper = 1,sym%Noper
      !
      i1 = 0
      do k1 = 0,j
         do tau1 = 0,min(1,k1)
           !if (k1==0.and.mod(j+tau1,2)/=0) cycle
           i1 = i1 + 1
           !
           n1=k1*(-1)**tau1
           !
           i2 = 0
           do k2 = 0,j
              do tau2 = 0,min(1,k2)
                !if (k2==0.and.mod(j+tau2,2)/=0) cycle
                i2 = i2 + 1
                !
                n2=k2*(-1)**tau2
                !
                phi   = sym%euler(ioper,1)
                theta = sym%euler(ioper,2)
                chi   = sym%euler(ioper,3)
                !
                T(i2,i1) =  ddlmn_conj(j,n2,n1,phi,theta,chi)
                !
              enddo
           enddo
           !
         enddo
      enddo
      !
      Q = matmul( transpose( conjg(R) ),matmul( T,R ) )
      !
      tmat(ioper,:,:) = real(Q,ark)
      !
    enddo
    !
    allocate (irr(sym%Nrepresen),stat=alloc)
    !
    call MLsymmetrization_rot(nroots,tmat,iverbose,Ntotal)
    !
    do isym = 1,sym%Nrepresen
      !
      allocate (irr(isym)%repres(max(Ntotal(isym),1),sym%degen(isym),nroots),stat=alloc)
      call ArrayStart('irr(isym)%coeffs',alloc,size(irr(isym)%repres),kind(irr(isym)%repres))
      !
      irr(isym)%repres = 0
      !
    enddo
    !
    call MLsymmetrization_rot(nroots,tmat,iverbose,Ntotal,irr)
    !
    iroot = 0
    icount = 0
    !
    do isym = 1,sym%Nrepresen
      !
      do i1 = 1,Ntotal(isym)
        !
        icount = icount + 1
        !
        do ideg = 1,sym%degen(isym)
          !
          iroot = iroot + 1
          ! 
          count_index(icount,ideg) = iroot
          count_degen(icount) = ideg
          !
          eigenvects(:,iroot) = irr(isym)%repres(i1,ideg,:)
          !
        enddo
        !
      enddo
      !
    enddo
    !
    do ioper = 1,sym%Noper
      !
      Q(:,:) = matmul( ( transpose(eigenvects) ),matmul( tmat(ioper,:,:),(eigenvects) ) )
      !
      tmat(ioper,:,:) = real(Q)
      !
    enddo
    !
    Ncount = icount
    Ntot = Ntotal
    !
    deallocate(tmat,R,T,Q,irr)
    !
    if (iverbose>=5) write(out,"('MLrotsymmetry_generate/end')") 
    !
    !
  end subroutine MLrotsymmetry_generate

   !
   !                               Symmetrization.  
   !
   !
   ! Construction of irreducible representation - to project the rotational 
   ! basis set into the symmetrized representaion
   !

  subroutine MLsymmetrization_rot(dimen,transform,iverbose,Ntotal,irr)

   integer(ik),intent(in) :: dimen
   real(ark),intent(in)   :: transform(:,:,:)
   integer(ik),intent(in) :: iverbose
   integer(ik),intent(out):: Ntotal(1:sym%Nrepresen)
   type(MOrepres_arkT),optional,intent(out)  :: irr(1:sym%Nrepresen)
   !
   integer(ik)             :: alloc,isym,ioper,ideg,iclasses,Nirr(sym%Nrepresen)
   integer(ik)             :: i0,ielem_t
   real(ark)               :: chi(sym%Nclasses),Nirr_rk(sym%Nrepresen),f_t,g_t
   real(ark),allocatable   :: projector(:,:,:,:)
   real(ark),allocatable   :: vect(:,:),t_vect(:,:)
   !
   integer(ik)             :: Nelem,Ndeg,ielem,jelem,felem,try_elem,jdeg,kdeg,jsym
   !
   allocate (projector(sym%Nrepresen,sym%Noper,sym%maxdegen,sym%maxdegen),stat=alloc)
   call ArrayStart('projector-rot',alloc,size(projector),kind(projector))
   !
   do isym =1,sym%Nrepresen
     !
     do ioper =1,sym%Noper
       !
       Ndeg = sym%degen(isym)
       !
       do ideg = 1,Ndeg
         !
         projector(isym,ioper,ideg,1:Ndeg) = sym%irr(isym,ioper)%repres(ideg,1:Ndeg)*real(sym%degen(isym),ark)/real(sym%Noper,ark)
         !
       enddo
     enddo 
   enddo 
   !
   ! We symmetrized the contracted basis set 
   ! by calculating irreducible represntations of the product of Nclasses the symmetrized components 
   ! Here the field 'represent' defines how the product of the contracted functions 
   ! transform with the group operations, i.e. the representation of basis set in the 
   ! contracted basis. 
   ! The field 'transform' will define the transformation from the contracted basis set 
   ! to the symmetrized basis set. 
   !
   ! first time only to estimate the total numbers (Ntotal) of irred. representations 
   !
   if (iverbose>=5.and..not.present(irr)) write(out,"('Counting the total number of irr. representations...')") 
   !
   Ntotal = 0 
   !
   if (iverbose>=5) then 
      write(out,"(//'Rotational irreducible representations')")
      write(out,"('  Find numbers of irr. representations...')")
   endif 
   !
   ! calculate the characters 
   !
   chi = 0
   ioper = 1
   !
   do iclasses =1,sym%Nclasses
     !
     do ideg = 1,dimen
       !
       chi(iclasses) = chi(iclasses) + transform(ioper,ideg,ideg)
       !
     enddo 
     !
     if (iverbose>=6) then 
       write(out,"('ioper = ',i5,'; character: ',f12.8)") ioper,chi(iclasses)
     endif 
     !
     ioper = ioper+sym%Nelements(iclasses)
     !
   enddo 
   !
   ! estimate the number of the irreducible representasions 
   !
   do isym =1,sym%Nrepresen
      !
      Nirr_rk(isym) = sum(real(sym%Nelements(:)*sym%characters(isym,:),ark)*chi(:))/real(sym%Noper,ark)
      Nirr(isym) = nint(Nirr_rk(isym),ik)
      !
   enddo 
   !
   if (all(Nirr(:)==0)) then 
     write(out,"('No symmetry defined, Nirr = 0')") 
     stop 'No symmetry defined, Nirr = 0'
   endif 
   !
   if (iverbose>=5) then 
     write(out,"('Number of irr. rep-s = ',30f12.4)") Nirr_rk(1:min(30,size(Nirr_rk)))
   endif 
   !
   Ntotal(:) = Nirr(:)
   !
   if (.not.present(irr)) return 
   !
   ! Now when we now the sizes we can allocate the matrices with all information 
   ! about the symmetry 
   !
   if (iverbose>=3) write(out,"('Generating the irr. representations...')") 
   !
   if (iverbose>=3) then 
      write(out,"(/'Total number of irr. representations: ',40i10)") Ntotal
   endif 
   !
   if (iverbose>=6) then 
      write(out,"(//'Construct IRREPS for each rotational level')")
   endif 
   !
   Nelem = dimen
   !
   allocate(vect(Nelem,Nelem),t_vect(Nelem,Nelem),stat=alloc)
   call ArrayStart('vect',alloc,size(vect),kind(vect))
   call ArrayStart('vect',alloc,size(t_vect),kind(t_vect))
   !
   do isym = 1,sym%Nrepresen
     !
     if (Nirr(isym)==0) cycle
     !
     Ndeg = sym%degen(isym)
     !
     ! Now we orthogonalize these vectors among themself, removing repeated ones
     !
     felem = 0
     try_elem = 0 
     !
     elem_loop2: do while (try_elem<Nelem.and.felem<Nirr(isym))
       !
       try_elem = try_elem + 1
       !
       do ideg = 1,Ndeg
         do jelem = 1,Nelem
           !
           t_vect(ideg,jelem) =  sum(projector(isym,1:sym%Noper,ideg,ideg)*transform(1:sym%Noper,try_elem,jelem))
           !
         enddo
       enddo
       !
       ielem_t = 0 
       !
       elem_loop: do while (ielem_t<Ndeg.and.felem<Nirr(isym))
         !
         ielem_t = ielem_t + 1
         !
         vect(1,1:Nelem) = t_vect(ielem_t,1:Nelem)
         !
         ! initial check for non-vanish 
         !
         f_t = sum(vect(1,:)**2)
         ! Continue for non-zero vectors
         !
         if (f_t>small_) then 
           !
           vect(1,:) = vect(1,:)/sqrt(f_t)
           !
           else
           !
           cycle elem_loop
           !
         endif
         !
         ! Schmidt orthogonalization
         !
         do ielem = 1,felem
           do  ideg= 1,sym%degen(isym)
             !
             f_t = sum(irr(isym)%repres(ielem,ideg,1:Nelem)*vect(1,1:Nelem))
             !
             vect(1,1:Nelem) = vect(1,1:Nelem)- irr(isym)%repres(ielem,ideg,1:Nelem)*f_t
             !
           enddo
         enddo
         !
         ! normalization 
         !
         f_t = sum(vect(1,:)**2)
         !
         ! Continue for non-zero vectors
         !
         if (f_t>small_) then 
           !
           vect(1,:) = vect(1,:)/sqrt(f_t)
           !
           else
           !
           cycle elem_loop
           !
         endif
         !
         ! Reconstructing other degenerate components
         !
         do jdeg= 2,Ndeg
           !
           vect(jdeg,1:Nelem) = 0
           !
           do ioper = 1,sym%Noper
            !
            do jelem = 1,Nelem
             !
             vect(jdeg,1:Nelem) =  vect(jdeg,1:Nelem) + &
                                   projector(isym,ioper,1,jdeg)*&
                                   transform(ioper,jelem,1:Nelem)*vect(1,jelem)
                                   !sym%irr(isym,ioper)%repres(ielem_t,jdeg)*&
                                   !real(sym%degen(isym),ark)/real(sym%Noper,ark)

            enddo

           enddo
           !
           !vect(jdeg,1:Nelem) = t_vect(jdeg,1:Nelem)
           !
         enddo
         !
         ! normalize 
         !
         do jdeg =2,Ndeg
           !
           f_t = sum(vect(jdeg,:)**2)
           !
           ! Continue for non-zero vectors
           !
           if (f_t>small_) then 
             !
             vect(jdeg,:) = vect(jdeg,:)/sqrt(f_t)
             !
           else
             !
             cycle elem_loop
             !
           endif
           !
           ! Orthogonalization with each other 
           !
           do kdeg = 1,jdeg-1
             !
             f_t = sum(vect(kdeg,:)*vect(jdeg,:))
             !
             vect(jdeg,:) = vect(jdeg,:)-vect(kdeg,:)*f_t
             !
           enddo
           !
         enddo
         !
         ! check if the found vector does transorm correctly with the group
         !         
         do ioper = 1,sym%Noper
           !
           do jdeg = 1,Ndeg
             !
             do jelem = 1,Nelem
               !
               f_t = sum( vect(jdeg,1:Nelem)*transform(ioper,1:Nelem,jelem) )
               !
               g_t = sum( sym%irr(isym,ioper)%repres(jdeg,1:Ndeg)*vect(1:Ndeg,jelem) )
               !
               ! Continue if this two quanttaties are the same 
               !
               if (abs(f_t-g_t)>100.0*sqrt(small_)) then 
                 !
                 continue
                 !
                 cycle elem_loop
                 !
               endif 
               !
             enddo
             !
           enddo
           !
         enddo
         !
         ! check if the found vector is really unique, by comparing also with other symmetries 
         !
         do jsym = 1,isym
           !
           do ielem = 1,min(Nirr(jsym),felem)
             !
             do kdeg =1,Ndeg
               !
               do jdeg =1,sym%degen(jsym)
                 !
                 if (all(abs(vect(kdeg,1:Nelem) - &
                             irr(jsym)%repres(ielem,jdeg,1:Nelem))<sqrt(small_))) cycle elem_loop
                 !
               enddo 
               !
             enddo 
             !
           enddo
           !
         enddo
         !
         ! if we are through all obstacles, we are finally there
         !
         felem = felem + 1
         !
         irr(isym)%repres(felem,1:Ndeg,1:Nelem) = vect(1:Ndeg,1:Nelem)
         !
        enddo elem_loop
        !
     enddo elem_loop2
     !
     ! Check if all representations have been found 
     !
     if (felem/=Nirr(isym)) then 
        write(out,"('degenerate_projectors: Not all irr. representations have been found for')")
        write(out,"(' isym = ',i4,' : only ',i7,' out of ',i7)") isym,felem,Nirr(isym)
        stop 'degenerate_projectors: Not all irr!. representations have been found'
     endif 
     !
  enddo 
  !
  deallocate (t_vect,vect,projector)
  call ArrayStop('vect')
  call ArrayStop('projector-rot')
  !
  if (iverbose>=6) then 
    !
    do isym =1,sym%Nrepresen
      !
      do i0 = 1,Nirr(isym)
          !
          do ideg = 1,sym%degen(isym)
             write(out,"(' isym =',i4,' ideg =',i4,' irr. repres.: ',30f18.8)") &
                           isym,ideg,irr(isym)%repres(i0,ideg,1:dimen)
          enddo 
        !
      enddo
      !
    enddo
    !
  endif
  !
  end subroutine MLsymmetrization_rot
  !
  !
  !
  ! define the symmetry of the rotational basis function 
  !
  subroutine MLrotsymmetry_generate_CII(J,count_index,count_degen,eigenval,Ncount)
    !
    integer(ik),intent(in)   :: J
    integer(ik),intent(out)  :: count_index(2*j+1,2*j+1),count_degen(2*j+1),eigenval(2*j+1),Ncount
    real(ark)                :: transform(2*j+1,sym%maxdegen,2*j+1)
    !
    real(ark)                :: theta,phi,chi,coeff,eigenvects(2*j+1,2*j+1),vect_t(2*j+1),e_t
    !
    integer(ik)              :: i1,k1,n1,i2,k2,n2,tau1,tau2,i,ioper,sigma,alloc,iroot,icount,ideg,jdeg,nroots
    integer(ik)              :: izeta,irep,jroot,count_sym(2*j+1),irep_t
    complex(ark)             :: temp,d
    complex(ark),allocatable :: R(:,:),T(:,:),Q(:,:)
    real(8),allocatable      :: da(:,:),db(:,:)
    real(ark),allocatable    :: a(:,:),b(:)
    !
    integer(ik)              :: Nirr(sym%Nrepresen),Nirr_elem
    integer(ik)              :: try_elem,felem,gamma,ielem,jelem,ielem_t,Ndeg,info,Nelem,kdeg
    real(ark)                :: t_vect(2*j+1,2*j+1),vect(2*j+1,2*j+1),fnormal,f_t,g_t
    real(ark),allocatable    :: tmat(:,:,:),rmat(:,:,:)
    !
    if (verbose>=5) write(out,"(/'MLrotsymmetry_generate_CII/start')") 
    !
    ! trivial case of no symmetry
    !
    if (trim(molec%symmetry)=='C'.or.trim(molec%symmetry)=='C(M)') then
       !
       !gamma = 1
       !ideg = 1
       return
       !
    endif
    !
    ! dimension of the problem
    !
    nroots  = 2*j+1
    !
    ! Obtain transformation-matrix 
    !
    allocate(R(nroots,nroots),T(nroots,nroots),Q(nroots,nroots),a(nroots,nroots),b(nroots),da(nroots,nroots),&
             db(nroots,1),stat=alloc)
    if (alloc/=0) stop 'MLrotsymmetry_generate_CII - out of memory'
    !
    allocate(tmat(sym%Noper,nroots,nroots),rmat(sym%Noper,nroots,nroots),stat=alloc)
    if (alloc/=0) stop 'MLrotsymmetry_generate - out of memory'
    !

    !
    ! R is a transformation to the Wang functions
    !
    R = 0 
    !
    if (mod(j,2)==0) then 
      !
      r(1,1) = cmplx(1.0_ark,0.0_ark)
      !
    else
      !
      r(1,1) = cmplx(0.0_ark,-1.0_ark)
      !
    endif
    !
    i1 = 2
    do k1 = 1,j
       !
       sigma = mod(k1,3)
       !
       R(i1  ,i1  ) = cmplx(1.0_ark/sqrt(2.0_ark)           ,0.0_ark                               )
       R(i1  ,i1+1) = cmplx(0.0_ark                         ,-(-1.0_ark)**sigma/sqrt(2.0_ark)      )
       R(i1+1,i1  ) = cmplx((-1.0_ark)**(j+k1)/sqrt(2.0_ark),0.0_ark                               )
       R(i1+1,i1+1) = cmplx(0.0_ark                         ,(-1.0_ark)**(j+k1+sigma)/sqrt(2.0_ark))
       !
       i1 = i1 + 2
       !
    enddo
    !
    ! transformation to the wang rot. basis 
    !
    !
    do ioper = 1,sym%Noper
      !
      i1 = 0
      do k1 = 0,j
         do tau1 = 0,1
           if (k1==0.and.mod(j+tau1,2)/=0) cycle
           i1 = i1 + 1
           !
           n1=k1*(-1)**tau1
           !
           i2 = 0
           do k2 = 0,j
              do tau2 = 0,1
                if (k2==0.and.mod(j+tau2,2)/=0) cycle
                i2 = i2 + 1
                !
                n2=k2*(-1)**tau2
                !
                phi   = sym%euler(ioper,1)
                theta = sym%euler(ioper,2)
                chi   = sym%euler(ioper,3)
                !
                T(i1,i2) =  ddlmn_conj(j,n1,n2,phi,theta,chi)
                !
              enddo
           enddo
           !
         enddo
      enddo
      !
      Q = matmul( transpose( conjg(R) ),matmul( T,R ) )
      !
      tmat(ioper,:,:) = real(Q,ark)
      !
    enddo
    !
    ! construct the CII matrix
    !
    do i1 = 1,nroots
       do i2 = 1,nroots
         !
         a(i1,i2) = sum( sym%CII%coeff( 1:sym%CII%Noper )*tmat( sym%CII%ioper( 1:sym%CII%Noper ),i1,i2 ) )
         !
       enddo
    enddo
    !
    da = real(a,8)
    !
    call lapack_syev(da,db(:,1))
    !
    eigenvects = real(da,ark)
    !
    eigenval(:) = nint(db(:,1),ik)
    !
    ! Estimation of the number of the unique levels. 
    !
    !
    count_index = 0
    count_index(1,1) = 1
    count_degen(1) = 1
    !
    icount = 1
    ideg = 1
    !
    do iroot = 1,nroots
      !
      izeta = 0
      !
      do irep = 1,sym%Nrepresen
        !
        do ideg = 1,sym%degen(irep)
          !
          izeta = izeta + 1
          !
          if (sym%CII%izeta(izeta)==eigenval(iroot) )  then
            !
            count_sym(iroot) = irep
            !
            exit
            !
          endif 
          !
        enddo
        !
      enddo
      !
    enddo
    !
    ! sort eigenvectros according with symmetry
    !
    !
    do iroot = 1,nroots
      !
      irep = count_sym(iroot)
      !
      Nirr(irep) = Nirr(irep) + 1
      !
      do jroot = iroot+1,nroots
        !
        if ( count_sym(jroot)<irep ) then
         !
         vect_t(:) = eigenvects(:,iroot)
         eigenvects(:,iroot) = eigenvects(:,jroot)
         eigenvects(:,jroot) = vect_t(:)
         !
         irep_t = count_sym(iroot)
         count_sym(iroot) = count_sym(jroot)
         count_sym(jroot) = irep_t
         !
         e_t = eigenval(iroot)
         eigenval(iroot) = eigenval(jroot)
         eigenval(jroot) = e_t
         !
        endif
        !
      enddo
      !
    enddo    
    !
    ! count degeneracies and unique levels 
    !
    count_index = 0
    count_index(1,1) = 1
    count_degen(1) = 1
    !
    icount = 1
    Nirr = 0
    ideg = 1
    !
    Nirr(count_sym(1)) = 1
    !
    do iroot = 2,nroots
      !
      irep = count_sym(iroot)
      !
      Nirr(irep) = Nirr(irep) + 1
      !
      if ( count_sym(iroot-1)==irep ) then
        !
        ideg = ideg + 1
        !
      else
        !
        ideg = 1
        icount = icount + 1
        !
      endif
      !
      count_index(icount,ideg) = iroot
      count_degen(icount) = ideg
      !
    enddo
    !
    Ncount = icount
    !
    ! transformtaion matrix of the eigenvectors
    !
    do ioper = 1,sym%Noper
      !
      !rmat(ioper,:,:) = matmul(transpose(eigenvects),tmat(ioper,:,:))
      !
      do iroot = 1,nroots
        !
        do i2 = 1,nroots
         !
         rmat(ioper,iroot,i2) = sum(eigenvects(:,iroot)*tmat(ioper,:,i2))
         !
        enddo
        !
      enddo
      !
    enddo
    !
    !
    ! Construct the irreps for the vectors with multiplicity
    !
    Nirr_elem = 0 
    !
    gamma_loop : do icount = 1,Ncount
      !
      gamma = count_sym(count_index(icount,1))
      !
      !if (Nirr(gamma)==0) cycle gamma_loop
      !
      ! Degeneracy of the term
      !
      Ndeg = sym%degen(gamma)
      Nelem = Nroots
      !
      do ideg=1,Ndeg
        !
        iroot = count_index(icount,ideg)
        !
        vect(ideg,1:Nelem) = eigenvects(1:Nelem,iroot)
        !
      enddo
      !
      try_elem = 0 
      felem = 0 
      t_vect = vect
      !
      elem_loop2: do while (try_elem<Nelem.and.felem<Nirr(gamma))
        !
        try_elem = try_elem + 1
        !
        ielem_t=0 
        !
        elem_loop: do while (ielem_t<Ndeg.and.felem<Nirr(gamma))
          !
          info = 0 
          !
          ! transformation of vector to irr. repres. 
          !
          vect = 0 
          !
          ielem_t = ielem_t + 1
          !
          do ioper = 1,sym%Noper
             !
             vect(1:Nelem,1:Nelem) =  vect(1:Nelem,1:Nelem) + &
                                      sym%irr(gamma,ioper)%repres(ielem_t,1)*&
                                      real(sym%degen(gamma),ark)/real(sym%Noper,ark)*&
                                      rmat(ioper,1:Nelem,1:Nelem)
            !
            !
          enddo
          !
          t_vect(1,:) = vect(try_elem,:)
          !
          !
          if (verbose>=6) then 
            do kdeg = 1,Nelem
               write(out,"('1,kdeg,t_vect: ',i8,f18.10)") kdeg,t_vect(1,kdeg)
            enddo
          endif 
          !
          ! Now we normalize these vectors 
          !
          fnormal = dot_product(t_vect(1,:),t_vect(1,:))

          if (fnormal>small_) then 
            !
            t_vect(1,:) = t_vect(1,:)/sqrt(fnormal)
            !
            else
             !
             if (verbose>=6) write(out,"('try_elem,Nelem,felem,Nirr: ',4i8)") try_elem,Nelem,felem,Nirr(gamma)
             !
             cycle elem_loop
             !
          endif
          !
          ! Schmidt orthogonalization
          !
          do ielem = 1,Nirr_elem
            do  ideg= 1,sym%maxdegen
              !
              f_t = sum(transform(ielem,ideg,1:Nelem)*t_vect(1,1:Nelem))
              !
              t_vect(1,1:Nelem) = t_vect(1,1:Nelem)- transform(ielem,ideg,1:Nelem)*f_t

              if (abs(f_t)>0.9) then 
                 if (verbose>=5) &
                   write(out,"('non-orthgon: i,gamma,ielem,ideg: ',4i6,f20.8)") Nirr_elem+1,gamma,ielem,ideg,f_t
                 info = 1
                 cycle elem_loop 
              endif
              !
              fnormal = dot_product(t_vect(1,1:Nelem),t_vect(1,1:Nelem))
              t_vect(1,1:Nelem) = t_vect(1,1:Nelem)/sqrt(fnormal)
              f_t = sum(transform(ielem,ideg,1:Nelem)*t_vect(1,1:Nelem))
              !
              if (abs(f_t)>sqrt(small_)) then 
                 if (verbose>=5) &
                   write(out,"('non-orthgon: i,gamma,ielem,ideg: ',4i6,f20.8)") Nirr_elem+1,gamma,ielem,ideg,f_t
                 info = 1
                 cycle elem_loop 
              endif
              !
            enddo
            !
          enddo
          !
          ! Now we normalize these vectors again 
          !
          fnormal = dot_product(t_vect(1,:),t_vect(1,:))
          !
          if (fnormal>small_) then 
            !
            t_vect(1,:) = t_vect(1,:)/sqrt(fnormal)
            !
            else
             !
             cycle elem_loop
             !
          endif
          !
          ! Re-constructing other degenerate components from the first one
          !
          do  ideg= 2,Ndeg
            !
            t_vect(ideg,1:Nelem) = 0
            !
            do ioper = 1,sym%Noper
               !
               t_vect(ideg,1:Nelem) =   t_vect(ideg,1:Nelem) + &
                                        sym%irr(gamma,ioper)%repres(ielem_t,ideg)*&
                                        real(sym%degen(gamma),ark)/real(sym%Noper,ark)*&
                                        rmat(ioper,try_elem,1:Nelem)
               !
             enddo
             !
          enddo

          !jdeg = 1
          !do  ideg= 1,Ndeg
          !  !
          !  if (ideg==try_elem) cycle
          !  !
          !  jdeg = jdeg + 1
          !  !
          !  t_vect(jdeg,:) = vect(ideg,:)
          !  !
          !enddo



          !
          ! Now we normalize these vectors 
          !
          do jdeg= 2,Ndeg
            !
            fnormal = dot_product(t_vect(jdeg,:),t_vect(jdeg,:))
            !
            if (fnormal>small_) then 
              !
              t_vect(jdeg,:) = t_vect(jdeg,:)/sqrt(fnormal)
              !
            endif
            !
          enddo
          !
          ! Check orthogonality and apply the Schmidt orthogonalization
          !
    !      do  jdeg= 1,Ndeg
    !        do ielem = 1,Nirr_elem
    !          do  ideg= 1,sym%maxdegen
    !            !
    !            f_t = sum(transform(ielem,ideg,1:Nelem)*t_vect(jdeg,1:Nelem))
    !            !
    !            if (abs(f_t)>0.9) then 
    !               if (verbose>=5) &
    !                 write(out,"('non-orthgon: i,gamma,jdeg,ielem,kdeg: ',5i6,f20.8)") Nirr_elem+1,gamma,jdeg,ielem,ideg,f_t
    !               info = 1
    !               cycle elem_loop 
    !            endif
    !            !
    !            t_vect(jdeg,1:Nelem) = t_vect(jdeg,1:Nelem)- transform(ielem,ideg,1:Nelem)*f_t
    !            fnormal = dot_product(t_vect(jdeg,1:Nelem),t_vect(jdeg,1:Nelem))
    !            t_vect(jdeg,1:Nelem) = t_vect(jdeg,1:Nelem)/sqrt(fnormal)
    !            f_t = sum(transform(ielem,ideg,1:Nelem)*t_vect(jdeg,1:Nelem))
    !            !
    !            if (abs(f_t)>sqrt(small_)) then 
    !               if (verbose>=5) &
    !               write(out,"('non-orthgon: i,gamma,jdeg,ielem,ideg: ',5i6,f20.8)") Nirr_elem+1,gamma,jdeg,ielem,ideg,f_t
    !                info = 1
    !                cycle elem_loop 
    !            endif
    !            
    !          enddo
    !        enddo
    !      enddo
          !
          ! check if the found vector does transorm correctly with the group
          !
          do jdeg = 1,Ndeg
            do kdeg = 1,Nelem
              !
              if (verbose>=6) then 
                !
                write(out,"('jdeg,kdeg,t_vect: ',2i8,f18.10)") jdeg,kdeg,t_vect(jdeg,kdeg)
                !
              endif 
              !         
              do ioper = 1,sym%Noper
                !
                f_t = sum( t_vect(jdeg,1:Nelem)*rmat(ioper,kdeg,1:Nelem))
                !
                g_t = sum( sym%irr(gamma,ioper)%repres(jdeg,1:Ndeg)*t_vect(1:Ndeg,kdeg) )
                !
                ! Continue if these two quantaties are the same 
                !
                if (abs(f_t-g_t)>0.01) then 
                   info = 1
                   if (verbose>=5) then 
                     write(out,"('jdeg,kdeg,ioper,f,g',3i6,2f20.8)") jdeg,kdeg,ioper,f_t,g_t
                   else
                     cycle elem_loop 
                   endif 
                endif
                !
              enddo
              !
            enddo
            !
          enddo
          !
          if (info==1) cycle elem_loop 
          !
          ! if we are through all obstacles, we are finally there
          !
          felem = felem + 1
          !
          if (info==0) then  
            !
            Nirr_elem = Nirr_elem + 1
            transform(Nirr_elem,1:Ndeg,1:Nelem) = t_vect(1:Ndeg,1:Nelem) 
            !try_elem = 0 
            t_vect = 0 
            !
          else 
            !
            Nirr(gamma) = -info 
            !
          endif 
          !
        enddo elem_loop
        !
      enddo elem_loop2
      !
    enddo  gamma_loop
    !
    !
    deallocate(R,T,Q,da,db,a,b,tmat,rmat)
    !
    select case(trim(molec%moltype))
    case default
         write (out,"('MLrotsymmetry_generate_CII: molecule type ',a,' unknown')") trim(molec%moltype)
         stop 'MLrotsymmetry_generate_CII - bad molecule'
         !
         !
    case('XY4') 
         !
         !
    end select 
    !
    if (verbose>=5) write(out,"('MLrotsymmetry_generate_CII/end')") 
    !
    !
  end subroutine MLrotsymmetry_generate_CII



   function dlmn(j,n,m,t,error) result (d)
     integer(ik),intent(in) :: j,n,m
     integer(ik),intent(out) :: error
     real(ark),intent(in)   :: t
     real(ark)              :: ss,prefactor,d,f,g,g1,g2,x,sinx,cosx,gh,h,hs,cos_sin,hlog
     integer(ik)            :: s
     
         error = 0
         s  =0
         ss =0
         !
         x = 0.5_ark*t ; cosx = cos(x) ; sinx = sin(x)
         !
         g2 = 0.5_ark*(faclog(j-m)+faclog(j+m)+faclog(j-n)+faclog(j+n))
         !
         do s = 0,min(j-m,j-n)
            !
            if (j-m-s>=0.and.j-n-s>=0.and.m+n+s>=0) then
              !
              !f=1/((j-m-s)!*(j-n-s)!*s!*(m+n+s)!);
              !
              g1 = faclog(j-m-s)+faclog(j-n-s)+faclog(s)+faclog(m+n+s)
              !
              g = g2-g1
              !
              cos_sin = cosx**(2*s+m+n)*sinx**(2*j-2*s-m-n)*(-1.0_ark)**s
              hs = sign(1.0_ark,cos_sin)
              h = abs(cos_sin)
              !
              if (h>0) then
                !
                hlog = log(h)
                !
                gh = g + hlog
                !
                if (gh>max_exp) then 
                  !write(out,"('dlmn error:  the value is too large for exp:',g18.9)")  g
                  error = 1
                  d = 0
                  return 
                  !stop 'dlmn error:  the value is too large for exp'
                endif
                !
                f = exp(gh)*hs
                !
                ss = ss + f
                !
              endif
              !
              !f = exp(-g)
              !
            endif
            !
         enddo
         !
         !g = faclog(j-m)+faclog(j+m)+faclog(j-n)+faclog(j+n)
         !
         !prefactor=exp(g)
         !
         d=(-1.0_ark)**(j-n)*ss  !*sqrt(prefactor)
         !
   end function dlmn



   !
   ! D function conjugated
   !
   function ddlmn_conj(j,n,m,phi,theta,chi,info) result (d)
     integer(ik),intent(in) :: j,m,n
     integer(ik),optional,intent(out) :: info
     real(ark),intent(in)   :: phi,chi,theta
     complex(ark)           :: f,d
     real(ark)              :: prefactor,x,y
     integer(ik)            :: error
     
         !
         prefactor=(-1.0_ark)**(n-m)*dlmn(j,-n,-m,theta,error)
         !
         x = real(n,ark)*phi ; y = real(m,ark)*chi
         !
         f=cmplx(cos(x+y),sin(x+y))
         !
         d = f*prefactor
         !
         if (present(info)) info = error
         if (.not.present(info).and.error/=0) then
            info = error
            write(out,"('dlmn error:  the value is too large for exp ')")
            stop 'dlmn error:  the value is too large for exp'
         endif
         !
   end function ddlmn_conj

     function Phi_rot(j,m,theta,phi) result (d)
     integer(ik),intent(in) :: j,m
     real(ark),intent(in)   :: phi,theta
     complex(ark)           :: f,d
     real(ark)              :: prefactor,x,y,g
     
         x = cos(theta)
         !
         g=plgndr_s(j,abs(m),x)
         !
         y = real(m,ark)*phi
         !
         f=cmplx(cos(y),sin(y))
         !
         prefactor = 0.5_ark*(faclog(j-abs(m))-faclog(j+abs(m)))
         !
         if (m<0) prefactor = prefactor + faclog(j-abs(m))-faclog(j+abs(m))
         !
         prefactor = exp(prefactor)*sqrt(real(2*j+1,ark)/(4.0_ark*pi))/sqrt(2.0_ark*pi)
         !
         d = g*f*prefactor
         !
         if (m<0) d=d*(-1.0_ark)**(-m)
         !
   end function Phi_rot



   function plgndr_s(l,m,x)
     implicit none
     integer(ik), intent(in) :: l,m
     real(ark), intent(in)   :: x
     real(ark) :: plgndr_s
     !Computes the associated Legendre polynomial Pm
     !l(x). Here m and l are integers satisfying
     ! 0 <= m <= l and -1 <= x <= 1
     integer(ik) :: ll
     real(ark) :: pll,pmm,pmmp1,somx2
       !
       if (m <0 .or. m > l  .or.  abs(x) > 1.0) stop 'plgndr_s args'
       !
       pmm=1.0_ark !Compute Pmm .
       !
       if (m > 0) then
           somx2=sqrt((1.0_ark-x)*(1.0_ark+x))
           pmm=product(arth_ark(1.0_ark,2.0_ark,m))*somx2**m
           if (mod(m,2) == 1) pmm=-pmm
       end if
       if (l == m) then
         plgndr_s=pmm
       else
         pmmp1=x*(2*m+1)*pmm ! Compute P^m_m+1.
         if (l == m+1) then
           plgndr_s=pmmp1
         else !Compute Pm_l , l>m + 1.
           do ll=m+2,l
              pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
              pmm=pmmp1
              pmmp1=pll
           end do
           plgndr_s=pll
         end if
       end if
   end function plgndr_s   
   
	function arth_ark(first,increment,n)
	real(ark), intent(in) :: first,increment
	integer(ik), intent(in) :: n
	real(ark), dimension(n) :: arth_ark
	integer(ik) :: k,k2
	integer(ik), parameter :: NPAR_ARTH=16,NPAR2_ARTH=8
	real(ark) :: temp
	if (n > 0) arth_ark(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_ark(k)=arth_ark(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_ark(k)=arth_ark(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_ark(k+1:min(k2,n))=temp+arth_ark(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
	END function arth_ark
   
   
    function calc_phirot(j,m,k,theta,phi) result (phirot)
     !
     integer(ik),intent(in) :: j,k,m
     real(ark),intent(in) :: theta,phi
     complex(rk)         :: phirot,f
     !
     integer(ik) :: s0,IJkm,ifact
     real(ark)   :: XJkm,fn_t,ct,st,minus,fact_t,y,x,cosx,sinx
        !
        x = 0.5_ark*theta ; cosx = cos(x) ; sinx = sin(x)
        !
        y = real(k,ark)*phi
        !
        f=cmplx(cos(y),sin(y))
        !
        Xjkm = faclog(J+m)+faclog(J-m)+faclog(J+k)+faclog(J-k)
        Xjkm = 0.5_ark*Xjkm
        !
        fn_t = 0
        do s0 = max(0,k-m),min(j-m,j+k)
           !
           !if (mod(s0,2)/=0) minus = -1.0_ark
           !
           !ct = 1.0_ark
           !if (2*J+k-m-2*s0/=0) 
           ct= cosx**(2*J+k-m-2*s0)
           !
           !st0 = 1.0_ark
           !if (m-k+2*s0/=0)     
           st=(-sinx)**(m-k+2*s0)
           !
           fact_t = faclog(s0)+faclog(J-m-s0)+faclog(m-k+s0)+faclog(J+k-s0)
           fact_t = Xjkm-fact_t
           !
           fact_t = exp(fact_t)
           !
           fn_t = fn_t+ (-1.0_ark)**s0*ct*st*fact_t
           !
        enddo 
        !
        phirot=fn_t*f*sqrt(real(2*J+1,ark)/(8.0_ark*pi**2))

     end function calc_phirot
   


!
! Here we define b0 equilibrium molecular cartesian coordinates in the xyz system 
! in the case of manifold rank = 1
!
  subroutine MLequilibrium_xyz_1d(Npoints,arho_borders,rhostep,periodic,rho_ref,b0,db0,rho_i)

     integer(ik),intent(in)  :: Npoints
     real(ark)   ,intent(in)  :: arho_borders(2)  ! rhomim, rhomax - borders
     real(ark)   ,intent(in)  :: rhostep
     logical    ,intent(in)   :: periodic
     real(ark)   ,intent(out)  :: rho_ref
     !
     real(ark),  intent(out) ::  b0(molec%Natoms,3,0:Npoints)
     real(ark),  intent(out) ::  db0(molec%Natoms,3,0:Npoints,3)
     real(ark),   intent(out) ::  rho_i(0:Npoints)
     !
     real(ark)                ::  rho_borders(2)  ! rhomim, rhomax - borders 
     !
     integer(ik)             :: i,ix

    if (verbose>=4) write(out,"(/'MLequilibrium_xyz_1d/start: a0 equilibrium cartesian coords')") 
    
    !case('GAUSS-LEGENDRE') 
    !
    !call gaulegf(rho_borders(1),rho_borders(2),rho_i,weight_i,Npoints+1)
    !
    !case('SIMSON-RULE')
    !
    !weight_i = 1.0_rk
    !
    rho_borders = arho_borders
    !
    do i = 0,npoints
       !
       rho_i(i) = real(i,kind=rk)*rhostep+rho_borders(1)
       !
    enddo
    !
    call MLequilibrium_xyz(Npoints,molec%Natoms,b0,rho_i,rho_ref,rho_borders)
    !
    ! computing derivatives only if the rank of the non-rigid manifold is 1
    !
    if (npoints/=0) then 
       !
       do i = 1,molec%Natoms
         do ix = 1,3
           !
           call diff_2d_4points_ark(Npoints,arho_borders,b0(i,ix,0:npoints),periodic,0_ik,db0(i,ix,0:npoints,1),&
                                    db0(i,ix,0:npoints,2))
           !
           call diff_3d_6points(Npoints,arho_borders,b0(i,ix,0:npoints),periodic,db0(i,ix,0:npoints,3))
           !
         enddo
       enddo
       !
       if (verbose>=5) then 
         write(out,"(/'Derivatives of b:')")
         do i = 0,npoints
           write(out,"(i8,18f12.8)") i,db0(:,:,i,1)
         enddo
       endif
       !
    else
        !
        db0 = 0
        !
    endif
    !
    if (verbose>=4) write(out,"('MLequilibrium_xyz_1d/end')") 
    !
  end subroutine MLequilibrium_xyz_1d
  !
  subroutine diff_2d_4points(Npoints,rho_b,f,periodic,d1f,d2f)

   integer(ik),intent(in)  :: Npoints
   real(rk)   ,intent(in)  :: f(0:Npoints)
   real(rk)   ,intent(in)  :: rho_b(2)
   logical    ,intent(in)  :: periodic
   !
   real(rk),   intent(out) ::  d1f(0:Npoints)
   real(rk),   optional    ::  d2f(0:Npoints)
   !
   real(ark)               :: fark(0:Npoints),d1fark(0:Npoints),d2fark(0:Npoints)

   real(ark)               :: rhostep,d1t,d2t,x,dy
   integer(ik)             :: ipoint,i,kl,kr
   integer(ik),parameter   :: Nextrap = 4
   integer(ik),parameter   :: Mextrap = 4
   real(ark)                :: v_t(-Nextrap:Nextrap),g_t(1:Nextrap)
     !
     fark = f
     !
     rhostep = real((rho_b(2)-rho_b(1)),ark)/real(npoints,kind=ark)
     !
     do ipoint = 0,npoints
        !
        kl = max(0,ipoint-4) ;   kr = min(Npoints,ipoint+4)
        !
        call ML_diffs_ark(kl,kr,ipoint,fark(kl:kr),rhostep,d1t,d2t) 
        !
        d1fark(ipoint) = d1t
        !
        if (present(d2f)) then
           !
           d2fark(ipoint) = d2t
           !
        endif
        !
     enddo
     !
     do i = 1,Nextrap
       !
       ipoint = Mextrap-1+i
       v_t(i) = rho_b(1)+real(ipoint,ark)*rhostep
       ipoint = Npoints-Mextrap-Nextrap+i
       g_t(i) = rho_b(1)+real(ipoint,ark)*rhostep
       !
     enddo

     do i = 0,Mextrap-1
       !
       ipoint = i
       x = rho_b(1)+real(ipoint,ark)*rhostep
       call polintark(v_t(1:Nextrap), d1fark(Mextrap:Mextrap+Nextrap-1), x, d1t, dy)
       d1fark(ipoint) = d1t
       !
       if (present(d2f)) & 
          call polintark(v_t(1:Nextrap), d2fark(Mextrap:Mextrap+Nextrap-1), x, d2fark(ipoint), dy)
       !
       ipoint = Npoints-Mextrap+1+i
       x = rho_b(1)+real(ipoint,ark)*rhostep
       !
       call polintark(g_t(1:Nextrap), d1fark(Npoints-Mextrap-Nextrap+1:Npoints-Mextrap), x, d1t, dy)
       !
       d1fark(ipoint) = d1t
       !
       if (present(d2f)) & 
          call polintark(g_t(1:Nextrap), d2fark(Npoints-Mextrap-Nextrap+1:Npoints-Mextrap), x, d2fark(ipoint), dy)
       !
     enddo
     !
     if (periodic) then
       !
       do ipoint = 0,3
         !
         v_t(0:4) = f(ipoint:ipoint+4)
         !
         do i = -4,-1
           !
           kl = mod(max(ipoint+i,npoints+ipoint+i),npoints)
           kl = mod(npoints+ipoint+i,npoints)
           !kl = mod(ipoint+i,npoints)
           !
           v_t(i) = f(kl)
           !
         enddo
         !
         call ML_diffs_ark(-4,4,0,v_t(-4:4),rhostep,d1t,d2t) 
         !
         d1fark(ipoint) = d1t
         !
         if (present(d2f)) then
            !
            d2fark(ipoint) = d2t
            !
         endif
         !
       enddo
       !
       do ipoint = npoints-3,npoints
         !
         v_t(-4:0) = f(ipoint-4:ipoint)
         !
         do i = 1,4
           !
           kl = mod(ipoint+i,npoints)
           !
           v_t(i) = f(kl)
           !
         enddo
         !
         call ML_diffs_ark(-4,4,0,v_t(-4:4),rhostep,d1t,d2t) 
         !
         d1fark(ipoint) = d1t
         !
         if (present(d2f)) then
            !
            d2fark(ipoint) = d2t
            !
         endif
         !
       enddo
       !
     endif 
     !
     d1f = d1fark 
     !
     if (present(d2f)) d2f = d2fark
     !
  end subroutine diff_2d_4points



  subroutine diff_2d_4points_ark(Npoints,rho_b,f,periodic,reflect,d1f,d2f)

   integer(ik) ,intent(in)  :: Npoints
   real(ark)   ,intent(in)  :: f(0:Npoints)
   real(ark)   ,intent(in)  :: rho_b(2)
   logical     ,intent(in)  :: periodic
   integer(ik),intent(in)   :: reflect
   !
   real(ark),   intent(out) ::  d1f(0:Npoints)
   real(ark),   optional    ::  d2f(0:Npoints)
   !
   real(ark)               :: rhostep,d1t,d2t,x,dy
   integer(ik)             :: ipoint,i,kl,kr
   integer(ik),parameter   :: Nextrap = 4
   integer(ik),parameter   :: Mextrap = 4
   real(ark)                :: v_t(-Nextrap:Nextrap),g_t(1:Nextrap)
     !
     !
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
     !
     do ipoint = 0,npoints
        !
        kl = max(0,ipoint-4) ;   kr = min(Npoints,ipoint+4)
        !
        call ML_diffs_ark(kl,kr,ipoint,f(kl:kr),rhostep,d1t,d2t) 
        !
        d1f(ipoint) = d1t
        !
        if (present(d2f)) then
           !
           d2f(ipoint) = d2t
           !
        endif
        !
     enddo
     !
     do i = 1,Nextrap
       !
       ipoint = Mextrap-1+i
       v_t(i) = rho_b(1)+real(ipoint,ark)*rhostep
       ipoint = Npoints-Mextrap-Nextrap+i
       g_t(i) = rho_b(1)+real(ipoint,ark)*rhostep
       !
     enddo

     do i = 0,Mextrap-1
       !
       ipoint = i
       x = rho_b(1)+real(ipoint,ark)*rhostep
       call polintark(v_t(1:Nextrap), d1f(Mextrap:Mextrap+Nextrap-1), x, d1t, dy)
       d1f(ipoint) = d1t
       !
       if (present(d2f)) & 
          call polintark(v_t(1:Nextrap), d2f(Mextrap:Mextrap+Nextrap-1), x, d2f(ipoint), dy)
       !
       ipoint = Npoints-Mextrap+1+i
       x = rho_b(1)+real(ipoint,ark)*rhostep
       !
       call polintark(g_t(1:Nextrap), d1f(Npoints-Mextrap-Nextrap+1:Npoints-Mextrap), x, d1t, dy)
       !
       d1f(ipoint) = d1t
       !
       if (present(d2f)) & 
          call polintark(g_t(1:Nextrap), d2f(Npoints-Mextrap-Nextrap+1:Npoints-Mextrap), x, d2f(ipoint), dy)
       !
     enddo
     !
     if (periodic) then
       !
       do ipoint = 0,3
         !
         v_t(0:4) = f(ipoint:ipoint+4)
         !
         do i = -4,-1
           !
           kl = mod(npoints+ipoint+i,npoints)
           !
           v_t(i) = f(kl)
           !
           if (reflect/=0.and.kl>4) v_t(i) = v_t(i)*real(reflect,ark)
           !
         enddo
         !
         call ML_diffs_ark(-4,4,0,v_t(-4:4),rhostep,d1t,d2t) 
         !
         d1f(ipoint) = d1t
         !
         if (present(d2f)) then
            !
            d2f(ipoint) = d2t
            !
         endif
         !
       enddo
       !
       do ipoint = npoints-3,npoints
         !
         v_t(-4:0) = f(ipoint-4:ipoint)
         !
         do i = 1,4
           !
           kl = mod(ipoint+i,npoints)
           !
           v_t(i) = f(kl)
           !
           if (reflect/=0.and.kl<npoints-4) v_t(i) = v_t(i)*real(reflect,ark)
           !
         enddo
         !
         call ML_diffs_ark(-4,4,0,v_t(-4:4),rhostep,d1t,d2t) 
         !
         d1f(ipoint) = d1t
         !
         if (present(d2f)) then
            !
            d2f(ipoint) = d2t
            !
         endif
         !
       enddo
       !
     elseif (reflect/=0) then
       !
       do ipoint = npoints-3,npoints
         !
         v_t(-4:0) = f(ipoint-4:ipoint)
         !
         do i = 1,4
           !
           if (ipoint+i<=npoints) then
             kl = ipoint+i
             v_t(i) = f(kl)
           else
             !
             kl = mod(ipoint+i,npoints)
             v_t(i) = f(npoints-kl)*real(reflect,ark)
             !
           endif 
           !
         enddo
         !
         call ML_diffs_ark(-4,4,0,v_t(-4:4),rhostep,d1t,d2t) 
         !
         d1f(ipoint) = d1t
         !
         if (present(d2f)) then
            !
            d2f(ipoint) = d2t
            !
         endif
         !
       enddo
       !
     endif 
     !
  end subroutine diff_2d_4points_ark




  subroutine diff_3d_6points(Npoints,rho_b,f,periodic,d3f)

   integer(ik) ,intent(in)  :: Npoints
   real(ark)   ,intent(in)  :: f(0:Npoints)
   real(ark)   ,intent(in)  :: rho_b(2)
   logical    ,intent(in)   :: periodic
   !
   real(ark)   ,intent(out) :: d3f(0:Npoints)
   !
   real(ark)               :: rhostep,d3t,x,dy
   integer(ik)             :: ipoint,i,kl,kr
   integer(ik),parameter   :: Nextrap = 6
   integer(ik),parameter   :: Mextrap = 6
   real(ark)               :: v_t(-Nextrap:Nextrap),g_t(1:Nextrap)
     !
     rhostep = real((rho_b(2)-rho_b(1)),ark)/real(npoints,kind=ark)
     !
     do ipoint = 0,npoints
        !
        kl = max(0,ipoint-6) ;   kr = min(Npoints,ipoint+6)
        !
        call ML_diffs_3d_ark(kl,kr,ipoint,f(kl:kr),rhostep,d3t) 
        !
        d3f(ipoint) = d3t
        !
     enddo
     !
     do i = 1,Nextrap
       !
       ipoint = Mextrap-1+i
       v_t(i) = rho_b(1)+real(ipoint,ark)*rhostep
       ipoint = Npoints-Mextrap-Nextrap+i
       g_t(i) = rho_b(1)+real(ipoint,ark)*rhostep
       !
     enddo

     do i = 0,Mextrap-1
       !
       ipoint = i
       x = rho_b(1)+real(ipoint,ark)*rhostep
       call polintark(v_t(1:Nextrap), d3f(Mextrap:Mextrap+Nextrap-1), x, d3t, dy)
       d3f(ipoint) = d3t
       !
       ipoint = Npoints-Mextrap+1+i
       x = rho_b(1)+real(ipoint,ark)*rhostep
       !
       call polintark(g_t(1:Nextrap), d3f(Npoints-Mextrap-Nextrap+1:Npoints-Mextrap), x, d3t, dy)
       !
       d3f(ipoint) = d3t
       !
     enddo
     !
     if (periodic) then
       !
       do ipoint = 0,5
         !
         v_t(0:6) = f(ipoint:ipoint+6)
         !
         do i = -6,-1
           !
           kl = mod(max(ipoint+i,npoints+ipoint+i),npoints)
           kl = mod(npoints+ipoint+i,npoints)
           !kl = mod(ipoint+i,npoints)
           !
           v_t(i) = f(kl)
           !
         enddo
         !
         call ML_diffs_3d_ark(-6,6,0,v_t(-6:6),rhostep,d3t) 
         !
         d3f(ipoint) = d3t
         !
       enddo
       !
       do ipoint = npoints-5,npoints
         !
         v_t(-6:0) = f(ipoint-6:ipoint)
         !
         do i = 1,6
           !
           kl = mod(ipoint+i,npoints)
           !
           v_t(i) = f(kl)
           !
         enddo
         !
         call ML_diffs_3d_ark(-6,6,0,v_t(-6:6),rhostep,d3t) 
         !
         d3f(ipoint) = d3t
         !
       enddo
       !
     endif 
     !

  end subroutine diff_3d_6points



subroutine ratint(xa, ya, x, y, dy)

  real(rk), dimension(:), intent(in) :: xa, ya
  real(rk), intent(in) :: x
  real(rk), intent(out) :: y, dy

  ! given arrays xa and ya of length n, and given a value of x, this routine
  ! returns a value of y and an accuracy estimate dy. the value returned is 
  ! that of the diagonal rational function, evaluated at x, that passes
  ! through the n points (xa_i, ya_i), i = 1... n.

  integer(ik) :: m, n, ns
  real(rk), dimension(size(xa)) :: c, d, dd, h, t
  real(rk), parameter :: tiny=1.0e-25_rk

  n = size(xa)
  if (size(ya)/=n) then 
      stop  'ratint: wrong sizes'
  endif 
  h = xa - x
  ns = minloc(abs(h),dim=1)
  y = ya(ns)
  if (x == xa(ns)) then
     dy = 0.0
     return
  end if
  c = ya
  d = ya + tiny               ! the tiny is needed to prevent 0/0
  ns = ns - 1
  do m=1, n-1
     t(1:n-m) = (xa(1:n-m)-x) * d(1:n-m)/h(1+m:n)   ! h will never be 0
     dd(1:n-m) = t(1:n-m) - c(2:n-m+1)
     !
     if (any(dd(1:n-m) == 0.0)) then                   ! interpolating function
          write(out,"('failure in ratint, has a pole here')")
          stop 'failure in ratint'
          !call nrerror('failure in ratint')         ! has a pole here
     endif 
     !
     dd(1:n-m) = (c(2:n-m+1) - d(1:n-m))/dd(1:n-m)
     d(1:n-m) = c(2:n-m+1) * dd(1:n-m)
     c(1:n-m) = t(1:n-m) * dd(1:n-m)
     if(2*ns < n-m) then
        dy = c(ns+1)
     else
        dy = d(ns)
        ns = ns-1
     end if
     y=y+dy
  end do

end subroutine ratint



subroutine MLratintark(xa, ya, x, y, dy)

  real(ark), dimension(:), intent(in) :: xa, ya
  real(ark), intent(in) :: x
  real(ark), intent(out) :: y, dy

  ! given arrays xa and ya of length n, and given a value of x, this routine
  ! returns a value of y and an accuracy estimate dy. the value returned is 
  ! that of the diagonal rational function, evaluated at x, that passes
  ! through the n points (xa_i, ya_i), i = 1... n.

  integer(ik) :: m, n, ns
  real(ark), dimension(size(xa)) :: c, d, dd, h, t
  !real(ark), parameter :: tiny=small_a

  n = size(xa)
  if (size(ya)/=n) then 
      stop  'MLratintark: wrong sizes'
  endif 
  h = xa - x
  ns = minloc(abs(h),dim=1)
  y = ya(ns)
  if (x == xa(ns)) then
     dy = 0
     return
  end if
  c = ya
  d = ya + small_a             ! is needed to prevent 0/0
  ns = ns - 1
  do m=1, n-1
     t(1:n-m) = (xa(1:n-m)-x) * d(1:n-m)/h(1+m:n)   ! h will never be 0
     dd(1:n-m) = t(1:n-m) - c(2:n-m+1)
     !
     if (any(dd(1:n-m) == 0.0)) then                   ! interpolating function
          write(out,"('failure in ratint, has a pole here')")
          stop 'failure in MLratintark'
     endif 
     !
     dd(1:n-m) = (c(2:n-m+1) - d(1:n-m))/dd(1:n-m)
     d(1:n-m) = c(2:n-m+1) * dd(1:n-m)
     c(1:n-m) = t(1:n-m) * dd(1:n-m)
     if(2*ns < n-m) then
        dy = c(ns+1)
     else
        dy = d(ns)
        ns = ns-1
     end if
     y=y+dy
  end do

end subroutine MLratintark



subroutine polint(xa, ya, x, y, dy)

  !use nrtype; use nrutil, only : assert_eq, iminloc, nrerror
  implicit none
  double precision, dimension(:), intent(in) :: xa, ya
  double precision, intent(in) :: x
  double precision, intent(out) :: y, dy

  ! given arrays xa and ya of length n, and given a value x, this routine
  ! returns a value y, and an error estimate dy. if p(x) is the polynomial
  ! of degree n - 1 such that p(x a_i) = y a_i, i = 1, ..., n, then the
  ! returned value y = p(x).

  integer(ik) :: m, n, ns
  real(rk), dimension(size(xa)) :: c, d, den, ho

  n = size(xa)
  if (size(ya)/=n) then 
      stop  'polint: wrong sizes'
  endif 

  !n = assert_eq(size(xa), size(ya), 'polint')
  c = ya                          ! initialize the tableau of c's and d's
  d = ya
  ho = xa - x
  ns = minloc(abs(x - xa),dim=1)  ! find index ns of closest table entry
  y = ya(ns)                      ! initial approximation to y.
  ns = ns - 1
  do m = 1, n - 1                        ! for each column of the tableau
     den(1:n-m) = ho(1:n-m) - ho(1+m:n)  ! we loop over c's and d's and
     !if (any(den(1:n-m) == 0.0)) &       ! update them
     !     call nrerror('polint: calculation failure')

     if (any(den(1:n-m) == 0.0)) then                   ! interpolating function
          write(out,"('failure in polint, has a pole here')")
          stop 'failure in polint'
          !call nrerror('failure in polint')         ! has a pole here
     endif 
     ! this error can occur only if two input xa's are (to within roundoff)
     ! identical.

     den(1:n - m) = (c(2:n-m+1) - d(1:n-m))/den(1:n-m)
     d(1:n-m) = ho(1+m:n) * den(1:n-m)   ! here c's and d's get updated
     c(1:n-m) = ho(1:n-m) * den(1:n-m)
     if (2 * ns < n-m) then       ! after each column in the tableau is
        dy=c(ns+1)                ! completed decide, which correction
     else                         ! c or d we add to y. we take the
        dy=d(ns)                  ! straightest line through the tableau
        ns=ns-1                   ! to its apex. the partial approximations
     end if                       ! are thus centred on x. the last dy
     y = y+dy                     ! is the measure of error.
  end do

end subroutine polint



recursive subroutine polintark(xa, ya, x, y, dy)

  !use nrtype; use nrutil, only : assert_eq, iminloc, nrerror
  implicit none
  real(ark), dimension(:), intent(in) :: xa, ya
  real(ark), intent(in) :: x
  real(ark), intent(out) :: y, dy

  !
  ! given arrays xa and ya of length n, and given a value x, this routine
  ! returns a value y, and an error estimate dy. if p(x) is the polynomial
  ! of degree n - 1 such that p(x a_i) = y a_i, i = 1, ..., n, then the
  ! returned value y = p(x).

  integer(ik) :: m, n, ns
  real(ark), dimension(size(xa)) :: c, d, den, ho

  if (verbose>=6) write(out,"('polintark...')")

  n = size(xa)
  if (size(ya)/=n) then 
      stop  'polint: wrong sizes'
  endif 

  !n = assert_eq(size(xa), size(ya), 'polint')
  c = ya                          ! initialize the tableau of c's and d's
  d = ya
  ho = xa - x
  ns = minloc(abs(x - xa),dim=1)  ! find index ns of closest table entry
  y = ya(ns)                      ! initial approximation to y.
  ns = ns - 1
  do m = 1, n - 1                        ! for each column of the tableau
     den(1:n-m) = ho(1:n-m) - ho(1+m:n)  ! we loop over c's and d's and
     !if (any(den(1:n-m) == 0.0)) &       ! update them
     !     call nrerror('polint: calculation failure')

     if (any(den(1:n-m) == 0.0)) then                   ! interpolating function
          write(out,"('failure in polint, has a pole here')")
          stop 'failure in polint'
          !call nrerror('failure in polint')         ! has a pole here
     endif 
     ! this error can occur only if two input xa's are (to within roundoff)
     ! identical.

     den(1:n - m) = (c(2:n-m+1) - d(1:n-m))/den(1:n-m)
     d(1:n-m) = ho(1+m:n) * den(1:n-m)   ! here c's and d's get updated
     c(1:n-m) = ho(1:n-m) * den(1:n-m)
     if (2 * ns < n-m) then       ! after each column in the tableau is
        dy=c(ns+1)                ! completed decide, which correction
     else                         ! c or d we add to y. we take the
        dy=d(ns)                  ! straightest line through the tableau
        ns=ns-1                   ! to its apex. the partial approximations
     end if                       ! are thus centred on x. the last dy
     y = y+dy                     ! is the measure of error.
  end do
  !
  if (verbose>=6) write(out,"('polintark.')")
  !

end subroutine polintark




  subroutine ML_diffs(n1,n2,n0,v,rhostep,d1,d2)
    !
    integer(ik),intent(in) :: n1,n2,n0
    real(rk),intent(in) :: v(n1:n2),rhostep
    real(rk),intent(out) :: d1,d2
    !
    real(ark)             :: v_t(-4:4)
    !
    !n1 = lbound(v,dim=1) ; n2 = ubound(v,dim=1)
    !
    if (n1>n0.or.n0>n2) then 
       !
       write (out,"('ML_diffs: must be n1<=n0<=n2, you give: ',3i8)") n1,n0,n2
       stop 'ML_diffs: wrong indexes n1,n0,n2'
       !
    endif 
    !
    if ( n0-n1==4 .and. n2-n0==4 ) then 
       !
       v_t(-4:4) = v(n1:n2)
       !
       d1 = (-v_t( 2)/12.0_ark+2.0_ark/3.0_ark*v_t( 1) & 
             +v_t(-2)/12.0_ark-2.0_ark/3.0_ark*v_t(-1) )/rhostep
       !
       d2 = ( v_t( 4)/144.0_ark-v_t( 3)/9.0_ark+v_t( 2)*4.0_ark/9.0_ark+v_t( 1)/9.0_ark  & 
           -v_t(0)*65.0_ark/72.0_ark  &
           +v_t(-4)/144.0_ark-v_t(-3)/9.0_ark+v_t(-2)*4.0_ark/9.0_ark+v_t(-1)/9.0_ark) &
           /rhostep**2
       !
    elseif ( n0-n1>=2 .and. n2-n0>=2 ) then
       !
       v_t(-2:2) = v(n0-2:n0+2)
       !
       d1 = (-v_t( 2)/12.0_ark+2.0_ark/3.0_ark*v_t( 1) & 
             +v_t(-2)/12.0_ark-2.0_ark/3.0_ark*v_t(-1) )/rhostep
       !
       d2 = (v_t( 2) + v_t(-2) - 2.0_ark*v_t(0) )/(rhostep**2*4.0_ark)
       !
    elseif ( n0-n1>=1 .and. n2-n0>=1 ) then
       !
       v_t(-1:1) = v(n0-1:n0+1)
       !
       d1 = ( v_t( 1) - v_t(-1) )/(rhostep*2.0_ark)
       !
       d2 = (v_t( 1) + v_t(-1) - 2.0_ark*v_t(0) )/(rhostep**2)
       !
    elseif ( n0-n1==0 .and. n2-n0>=2 ) then
       !
       v_t(0:2) = v(n0:n0+2)
       !
       d1 = ( v_t( 1) - v_t(0) )/rhostep
       !
       d2 = (v_t( 0) - 2.0_ark*v_t( 1) + v_t( 2) )/(rhostep**2)
       !
    elseif ( n0-n1>=2 .and. n2-n0==0 ) then
       !
       v_t(-2:0) = v(n0-2:n0)
       !
       d1 = ( v_t(0) - v_t(-1) )/rhostep
       !
       d2 = ( v_t( 0) - 2.0_ark*v_t(-1) + v_t(-2) )/(rhostep**2)
       !
    else
       !
       write (out,"('ML_diffs: wrong n1,n0,n2: ',3i8)") n1,n0,n2
       stop 'ML_diffs: wrong indexes n1,n0,n2'
       !
    endif
    !
  end subroutine ML_diffs



  subroutine ML_diffs_ark(n1,n2,n0,v,rhostep,d1,d2)
    !
    integer(ik),intent(in) :: n1,n2,n0
    real(ark),intent(in) :: v(n1:n2),rhostep
    real(ark),intent(out) :: d1,d2
    !
    real(ark)             :: v_t(-4:4)
    !
    !n1 = lbound(v,dim=1) ; n2 = ubound(v,dim=1)
    !
    if (n1>n0.or.n0>n2) then 
       !
       write (out,"('ML_diffs_ark: must be n1<=n0<=n2, you give: ',3i8)") n1,n0,n2
       stop 'ML_diffs_ark: wrong indexes n1,n0,n2'
       !
    endif 
    !
    if ( n0-n1==4 .and. n2-n0==4 ) then 
       !
       v_t(-4:4) = v(n1:n2)
       !
       d1 = (-v_t( 2)/12.0_ark+2.0_ark/3.0_ark*v_t( 1) & 
             +v_t(-2)/12.0_ark-2.0_ark/3.0_ark*v_t(-1) )/rhostep
       !
       d2 = ( v_t( 4)/144.0_ark-v_t( 3)/9.0_ark+v_t( 2)*4.0_ark/9.0_ark+v_t( 1)/9.0_ark  & 
           -v_t(0)*65.0_ark/72.0_ark  &
           +v_t(-4)/144.0_ark-v_t(-3)/9.0_ark+v_t(-2)*4.0_ark/9.0_ark+v_t(-1)/9.0_ark) &
           /rhostep**2
       !
    elseif ( n0-n1>=2 .and. n2-n0>=2 ) then
       !
       v_t(-2:2) = v(n0-2:n0+2)
       !
       d1 = (-v_t( 2)/12.0_ark+2.0_ark/3.0_ark*v_t( 1) & 
             +v_t(-2)/12.0_ark-2.0_ark/3.0_ark*v_t(-1) )/rhostep
       !
       d2 = (v_t( 2) + v_t(-2) - 2.0_ark*v_t(0) )/(rhostep**2*4.0_ark)
       !
    elseif ( n0-n1>=1 .and. n2-n0>=1 ) then
       !
       v_t(-1:1) = v(n0-1:n0+1)
       !
       d1 = ( v_t( 1) - v_t(-1) )/(rhostep*2.0_ark)
       !
       d2 = (v_t( 1) + v_t(-1) - 2.0_ark*v_t(0) )/(rhostep**2)
       !
    elseif ( n0-n1==0 .and. n2-n0>=2 ) then
       !
       v_t(0:2) = v(n0:n0+2)
       !
       d1 = ( v_t( 1) - v_t(0) )/rhostep
       !
       d2 = (v_t( 0) - 2.0_ark*v_t( 1) + v_t( 2) )/(rhostep**2)
       !
    elseif ( n0-n1>=2 .and. n2-n0==0 ) then
       !
       v_t(-2:0) = v(n0-2:n0)
       !
       d1 = ( v_t(0) - v_t(-1) )/rhostep
       !
       d2 = ( v_t( 0) - 2.0_ark*v_t(-1) + v_t(-2) )/(rhostep**2)
       !
    else
       !
       write (out,"('ML_diffs_ark: wrong n1,n0,n2: ',3i8)") n1,n0,n2
       stop 'ML_diffs_ark: wrong indexes n1,n0,n2'
       !
    endif
    !
  end subroutine ML_diffs_ark




  subroutine ML_diffs_3d_ark(n1,n2,n0,v,rhostep,d3)
    !
    integer(ik),intent(in) :: n1,n2,n0
    real(ark),intent(in) :: v(n1:n2),rhostep
    real(ark),intent(out) :: d3
    !
    real(ark)             :: v_t(-6:6)
    !
    !n1 = lbound(v,dim=1) ; n2 = ubound(v,dim=1)
    !
    if (n1>n0.or.n0>n2) then 
       !
       write (out,"('ML_diffs_3d_ark: must be n1<=n0<=n2, you give: ',3i8)") n1,n0,n2
       stop 'ML_diffs_3d_ark: wrong indexes n1,n0,n2'
       !
    endif 
    !
    if ( n0-n1==6 .and. n2-n0==6 ) then 
       !
       v_t(-6:6) = v(n1:n2)

       d3 = -1.0_ark/1728.0_ark*(v_t(6)-387.0_ark*v_t(2)-192.0_ark*v_t(-4)+192.0_ark*v_t(4)+&
            387.0_ark*v_t(-2)-24.0_ark*v_t(5)+&
            24.0_ark*v_t(-5)-v_t(-6)+488.0_ark*v_t(-3)-1584.0_ark*v_t(-1)-&
            488.0_ark*v_t(3)+1584.0_ark*v_t(1))/rhostep**3
       !
    elseif ( n0-n1>=3 .and. n2-n0>=3 ) then
       !
       v_t(-3:3) = v(n0-3:n0+3)
       !
       d3 = -0.125*(-v_t(3)+3.0_ark*v_t(1)-3.0_ark*v_t(-1)+v_t(-3))/rhostep**3
       !
    elseif ( n2-n0>=3 ) then
       !
       v_t(0:3) = v(n0:n0+3)
       !
       d3  = (v_t(3)-3.0_ark*v_t(2)+3.0_ark*v_t(1)-v_t(0))/rhostep**3
       !
    elseif ( n0-n1>=3  ) then
       !
       v_t(-3:0) = v(n0-3:n0)
       !
       d3  = (v_t(0)-3.0_ark*v_t(-1)+3.0_ark*v_t(-2)-v_t(-3))/rhostep**3
       !
    else
       !
       write (out,"('ML_diffs_3d_ark: wrong n1,n0,n2: ',3i8)") n1,n0,n2
       stop 'ML_diffs_3d_ark: wrong indexes n1,n0,n2'
       !
    endif
    !
  end subroutine ML_diffs_3d_ark



!
! We use the following procedure, which is assumed 
! to define the molecule:
! equilibrium geometry, potential function, masses, and something else     
!
  !
  ! Defining the rho-coordinate 
  !
  function MLcoord_direct(x,itype,imode,iorder)  result(v)

   real(ark),intent(in)   ::  x
   integer(ik),intent(in) :: itype
   integer(ik),intent(in) :: imode
   integer(ik),optional   :: iorder
   real(ark)              ::  rhoe,v,amorse
     !
     if (verbose>=6) write(out,"(/'MLcoord_direct/start')") 
     !
     select case(trim(molec%coordinates(itype,imode)))
     case default
        write (out,"('MLcoord_direct: coordtype ',a,' unknown')") trim(molec%coordinates(itype,imode))
        stop 'MLcoord_direct - bad coordtype'
     case('COSTAU') 
        !
        v = cos(x)
        !
     case('COSTAU2') 
        !
        v = cos(x)**2
        !
     case('COSRHO') 
        !
        !rhoe = molec%specparam(imode)
        !
        rhoe = molec%chi_eq(imode)
        !
        v = cos(rhoe)-cos(x)
        !
     case('COSX') 
        !
        v = cos(x)
        !
     case('1-COSX') 
        !
        v = 1.0_ark-cos(x)
        !
     case('SINX') 
        !
        v = sin(x)
        !
     case('SINRHO') 
        !
        !rhoe = molec%specparam(imode)
        !
        rhoe = molec%chi_eq(imode)
        !
        v = sin(rhoe)-sin(x)
        !
     case('LINCOSRHO') 
        !
        !rhoe = molec%specparam(imode)
        !
        rhoe = molec%chi_eq(imode)
        !
        v = cos(rhoe)-cos(x+rhoe)
        !
     case('LINSINRHO') 
        !
        !rhoe = molec%specparam(imode)
        !
        rhoe = molec%chi_eq(imode)
        !
        v = sin(rhoe)-sin(x+rhoe)
        !
     case('LINEAR','HARMONIC','NORMAL','SAME') 
        !
        v = x
        !
     case('X-XE') 
        !
        v = x-molec%chi_eq(imode)
        !
     case('MORSE') 
        !
        amorse = molec%specparam(imode)
        !
        v = 1.0_ark-exp( -amorse*( x ) )
        !
     case('RATIONAL') 
        !
        v = x
        !
     case('BOND-LENGTH', 'ANGLE', 'DIHEDRAL')
        !
        v = x
        !
     end select
     !
     if (present(iorder)) then 
       !
       select case(trim(molec%coordinates(itype,imode)))
          !
       case default
          !
          v = v**iorder 
          !
       case('RATIONAL') 
          !
          if (iorder<0) stop 'MLcoord_direct error: negative iorder'
          !
          select case(iorder) 
            !
          case (0)
            !
            v = 1.0_ark
            !
          case (1)
            !
            v = 1.0_ark/(molec%local_eq(imode)+x)
            !
          case (2)
            !
            v = 1.0_ark/(molec%local_eq(imode)+x)**2
            !
          case default
            !
            v = x**(iorder-2)
            !
          end select 
          !
        case('BOND-LENGTH')
          !
          if (iorder < 0) then
            print*,'MLcoord_direct error: negative iorder'
            stop 'MLcoord_direct error: negative iorder'
          endif
          !
          select case(iorder)
            !
            case(0)
              !
              v = 1.0_ark
              !
            case(1)
              !
              v = 1.0_ark/(molec%local_eq(imode) + x)
              !
            case(2)
              !
              v = 1.0_ark/(molec%local_eq(imode) + x)**2
              !
            case default
              !
              v = 1.0_ark
              !
          end select 
          !
        case('ANGLE')
           !
           if(iorder < 0) then 
              print*, 'MLcoord_direct error: negative iorder'
              stop 'MLcoord_direct error: negative iorder'
           endif
           !
          select case(iorder)
            !
            case(0)          
              !
              v = 1.0_ark
              !
            case(1)
              !
              v = Cos(molec%local_eq(imode) + x) 
              !
            case(2)
              !
              v = 1.0_ark/Tan(molec%local_eq(imode) + x)
              !
            case(3)
              !
              v = 1.0_ark/Sin(molec%local_eq(imode) + x)
              !
            case(4)
              !
              v = Sin(molec%local_eq(imode) + x)
              !
            case(5)
              !
              v = 1.0_ark/Tan(molec%local_eq(imode) + x)**2
              !
            case(6)
              !
              v = 1.0_ark/(Sin(molec%local_eq(imode) +x)*Tan(molec%local_eq(imode) + x))
              !
            case(7)
              !
              v = 1.0_ark/(Sin(molec%local_eq(imode) + x))**2
              !
            case default
              !
              v = 1.0_ark
              !
           end select
           !
        case('DIHEDRAL')
           !
           if(iorder < 0) then 
              print*, 'MLcoord_direct error: negative iorder'
              stop 'MLcoord_direct error: negative iorder'
           endif
            !
           select case(iorder) 
            !
           case(0)
              !
              v = 1.0_ark
              !
           case(1)
              !
              v = Cos((molec%local_eq(imode) +x)/2.0_ark)
              !
           case(2)
              !
              v = Sin((molec%local_eq(imode) +x)/2.0_ark)
              !
            case(3)
              !
              v = Cos((molec%local_eq(imode) +x)/2.0_ark)**2
              !
            case(4)
              !
              v = Cos((molec%local_eq(imode) + x)/2.0_ark)*Sin((molec%local_eq(imode)+x)/2.0_ark)
              !
            case(5)
              !
              v = Sin((molec%local_eq(imode) + x)/2.0_ark)**2
              !
            case default
              !
              v = 1.0_ark
              !
          end select
          !
       end select
       !
     endif 
     !
     if (verbose>=6) write(out,"('MLcoord_direct/end')") 
     !
 end function MLcoord_direct

  !
  ! Defining the rho-coordinate 
  !
  function MLcoord_firstderiv(xi,imode)  result(dxi_dchi)

   real(ark),intent(in)    ::  xi(1:molec%Nmodes)
   integer(ik),intent(in)  :: imode

   integer(ik)            :: itype
   real(ark)            ::  rhoe,dxi_dchi,amorse,chi

     if (verbose>=6) write(out,"(/'MLcoord_firstderiv/start')") 
     !
     ! We use this routine only for the kinetic energy expansion
     !
     itype = 1
     !
     chi = MLcoord_invert(xi,itype,imode)
     !
     select case(trim(molec%coordinates(itype,imode)))
     case default
        write (out,"('MLcoord_firstderiv: coordtype ',a,' unknown')") trim(molec%coordinates(itype,imode))
        stop 'MLcoord_firstderiv - bad coordtype'
     case('COSTAU') 
        !
        dxi_dchi = -sin(chi)
        !
     case('COSTAU2') 
        !
        dxi_dchi = -2.0_ark*cos(chi)*sin(chi)
        !
     case('COSRHO') 
        !
        dxi_dchi = -sin(chi)
        !
     case('SINRHO') 
        !
        dxi_dchi = cos(chi)
        !
     case('LINCOSRHO') 
        !
        rhoe = molec%chi_eq(imode)
        !
        dxi_dchi = -sin(rhoe+chi)
        !
     case('LINSINRHO') 
        !
        rhoe = molec%chi_eq(imode)
        !
        dxi_dchi = cos(rhoe+chi)
        !
     case('LINEAR','HARMONIC','NORMAL','X-XE','SAME') 
        !
        dxi_dchi = 1.0_rk
        !
     case('MORSE') 
        !
        amorse = molec%specparam(imode)
        !
        dxi_dchi = amorse*exp( -amorse*( chi ) )
        !
     end select 
     !
     if (verbose>=6) write(out,"('MLcoord_firstderiv/end')") 
 
 end function MLcoord_firstderiv


 !
 ! Coordinate transformation from linearized coordinates to lin-function coordinates, 
 ! where function = 'MORSE', 'csrho', ...
 ! It is a one-dimensional transformation so far
 !

 function MLcoord_invert(xi,itype,imode) result (chi)

   real(ark),intent(in) :: xi(1:molec%Nmodes)
   integer(ik),intent(in) :: itype
   integer(ik),intent(in) :: imode
   
   real(ark) :: chi,amorse,rhoe

   if (verbose>=6) write(out,"(/'MLcoord_invert/start')")
    
   !write (out,"('imode=',i)") imode

   select case(trim(molec%coordinates(itype,imode)))
   case default
        write (out,"('MLcoord_invert: coordinate type ',a,' unknown')") trim(molec%coordinates(itype,imode))
        stop 'MLcoord_invert - bad coordinate-type'
   case('MORSE') 
        !
        ! xi(i) = 1/( i*amorse )
        !
        amorse = molec%specparam(imode)
        !
        chi = -log(1.0_ark-xi(imode))/amorse
        !
    case('HARMONIC','NORMAL','LINEAR','SAME') 
        !
        chi = xi(imode)
        !
    case('X-XE') 
        !
        chi = xi(imode)+ molec%chi_eq(imode)
        !
    case('COSRHO') 
        !
        rhoe = molec%chi_eq(imode)
        !
        if (abs(cos(rhoe)-xi(imode))>1.0_ark) then 
            !
            write (out,"('MLcoord_invert: |cos x| = cos rhoe -xi  > 1: ',f18.8)") cos(rhoe)-xi(imode)
            write (out,"('Consider change difftype ')")
            stop 'MLcoord_invert - bad cos x'
            !
        endif 
        !
        chi =acos( cos(rhoe)-xi(imode))
        !
     case('COSTAU') 
        !
        if (abs(xi(imode))>1.0_ark) then 
            !
            write (out,"('MLcoord_invert: |cos tau| =  xi  > 1: ',f18.8)") xi(imode)
            write (out,"('Consider change difftype ')")
            stop 'MLcoord_invert - bad cos x'
            !
        endif 
        !
        chi =acos( xi(imode) )
        !
     case('COSX') 
        !
        if (abs(xi(imode))>1.0_ark) then 
            !
            write (out,"('MLcoord_invert: |cos x| =  xi  > 1: ',f18.8)") xi(imode)
            write (out,"('Consider change difftype ')")
            stop 'MLcoord_invert - bad cos x'
            !
        endif 
        !
        chi =acos( xi(imode) )
        !
     case('1-COSX') 
        !
        if (xi(imode)<0.0_ark) then 
            !
            write (out,"('MLcoord_invert: 1-cos x <0: ',f18.8)") xi(imode)
            write (out,"('Consider change difftype ')")
            stop 'MLcoord_invert - bad 1- cos x'
            !
        endif 
        !
        chi =acos( 1.0_ark-xi(imode) )
        !
     case('COSTAU2') 
        !
        if ((xi(imode))<0.0_ark) then 
            !
            write (out,"('MLcoord_invert: costau**2 < 0 : ',f18.8)") xi(imode)
            stop 'MLcoord_invert - costau**2 negative'
            !
        endif 
        !
        if (sqrt(xi(imode))>1.0_ark) then 
            !
            write (out,"('MLcoord_invert: costau**2 > 1 : ',f18.8)") sqrt(xi(imode))
            stop 'MLcoord_invert - costau**2 >1'
            !
        endif 
        !
        chi =acos( sqrt(xi(imode)) )
        !
     case('SINRHO')
        !
        rhoe = molec%chi_eq(imode)
        !
        if (abs(sin(rhoe)-xi(imode))>1.0_ark) then 
            !
            write (out,"('MLcoord_invert: |sin x| = sin rhoe -xi  > 1: ',f18.8)") sin(rhoe)-xi(imode)
            write (out,"('Consider change difftype ')")
            stop 'MLcoord_invert - bad sin x'
            !
        endif 
        !
        chi = asin( sin(rhoe)-xi(imode))
        !
        if (rhoe>0.5_ark*pi.and.chi<0.5_ark*pi) chi = pi - chi
        !
     case('LINCOSRHO') 
        !
        rhoe = molec%chi_eq(imode)
        !
        if (abs(cos(rhoe)-xi(imode))>1.0_ark) then 
            !
            write (out,"('MLcoord_invert: |cos x| = cos rhoe -xi_lin  > 1: ',f18.8)") cos(rhoe)-xi(imode)
            write (out,"('Consider change difftype ')")
            stop 'MLcoord_invert - bad cos x'
            !
        endif 
        !
     case('LINSINRHO') 
        !
        rhoe = molec%chi_eq(imode)
        !
        if (abs(sin(rhoe)-xi(imode))>1.0_ark) then 
            !
            write (out,"('MLcoord_invert: |sin x| = sin rhoe -xi_lin  > 1: ',f18.8)") sin(rhoe)-xi(imode)
            write (out,"('Consider change difftype ')")
            stop 'MLcoord_invert - bad sin x'
            !
        endif 
        !
    end select 
   
   if (verbose>=6) write(out,"('MLcoord_invert/end')") 
   
    
 end function MLcoord_invert


 subroutine ML_check_steps4coordinvert(xi,itype,imode,fstep)

   real(ark),intent(in) :: xi(1:molec%Nmodes)
   integer(ik),intent(in) :: itype,imode
   real(ark),intent(out) :: fstep(2)
   real(ark) :: rhoe

   if (verbose>=6) write(out,"(/'ML_check_steps4coordinvert/start')") 
   !
   fstep(1:2) = 1.0_ark
   !
   select case(trim(molec%coordinates(itype,imode)))
   case default
        write (out,"('ML_check_steps4coordinvert: coordinate type ',a,' unknown')") trim(molec%coordinates(itype,imode))
        stop 'ML_check_steps4coordinvert - bad coordinate-type'
   case('MORSE') 
        !
        fstep(1:2) = 1.0_ark
        !
    case('HARMONIC','NORMAL','LINEAR','SAME','X-XE') 
        !
        fstep(1:2) = 1.0_ark
        !
     case('COSRHO') 
        !
        rhoe = molec%chi_eq(imode)
        !
        if (abs( cos(rhoe)-xi(imode)-1.0_ark )<small_) then 
            !
            fstep(2) = 0.0_ark
            !
        endif 
        if (abs( cos(rhoe)-xi(imode)+1.0_ark )<small_) then 
            !
            fstep(1) = 0.0_ark
            !
        endif 
        !
     case('COSTAU') 
        !
        if (abs( xi(imode)-1.0_ark )>small_) then 
            !
            fstep(2) = 0.0_ark
            !
        endif 
        if (abs( xi(imode)+1.0_ark )<small_) then 
            !
            fstep(1) = 0.0_ark
            !
        endif 
        !
     case('COSTAU2') 
        !
        if (abs( xi(imode) )<small_) then 
            !
            fstep(1) = 0.0_ark
            !
        endif 
        if (abs( xi(imode)**2-1.0_ark )<small_) then 
            !
            fstep(2) = 0.0_ark
            !
        endif 
        !
     case('SINRHO')
        !
        rhoe = molec%chi_eq(imode)
        !
        if (abs( sin(rhoe)-xi(imode)-1.0_ark )>small_) then 
            !
            fstep(2) = 0.0_ark
            !
        endif 
        if (abs( sin(rhoe)-xi(imode)+1.0_ark )<small_) then 
            !
            fstep(1) = 0.0_ark
            !
        endif 
        !
     case('LINCOSRHO') 
        !
        rhoe = molec%chi_eq(imode)
        !
        if (abs( cos(rhoe)-xi(imode)-1.0_ark )>small_) then 
            !
            fstep(2) = 0.0_ark
            !
        endif 
        if (abs( cos(rhoe)-xi(imode)+1.0_ark )<small_) then 
            !
            fstep(1) = 0.0_ark
            !
        endif 
        !
     case('LINSINRHO') 
        !
        rhoe = molec%chi_eq(imode)
        !
        if (abs( sin(rhoe)-xi(imode)-1.0_ark )>small_) then 
            !
            fstep(2) = 0.0_ark
            !
        endif 
        if (abs( sin(rhoe)-xi(imode)+1.0_ark )<small_) then 
            !
            fstep(1) = 0.0_ark
            !
        endif
        !
     case('RATIONAL')
        !
        fstep(1:2) = 1.0_ark
        !
     case('BOND-LENGTH', 'ANGLE', 'DIHEDRAL')
        !
        write(out, "('ML_check_steps4coordinvert (BLAD) ','a', 'not applicable')")  trim(molec%coordinates(itype,imode))
        stop 'ML_check_steps4coordinvert - bad coordinate-type'
        !
    end select 
    !
    if (fstep(1)+fstep(2)==0.0_ark) then 
       write (out,"('ML_check_steps4coordinvert: no numerical derivatives allowed around point ',f18.8)") xi(imode)
       write (out,"('imode -  ',i8)") imode
       stop 'ML_check_steps4coordinvert - bad point for derivatives'
    endif
    !
    if (verbose>=6) write(out,"('ML_check_steps4coordinvert/end')") 
    !    
 end subroutine ML_check_steps4coordinvert


! gauleg.f90     P145 Numerical Recipes in Fortran
! compute x(i) and w(i)  i=1,n  Legendre ordinates and weights
! on interval -1.0 to 1.0 (length is 2.0)
! use ordinates and weights for Gauss Legendre integration
!
subroutine gaulegf(x1, x2, x, w, n)
  integer(ik), intent(in) :: n
  real(rk), intent(in) :: x1, x2
  real(rk), dimension(n), intent(out) :: x, w
  integer(ik) :: i, j, m
  real(rk) :: p1, p2, p3, pp, xl, xm, z, z1
  real(rk), parameter :: eps=3.d-14
      
  m = (n+1)/2
  xm = 0.5_rk*(x2+x1)
  xl = 0.5_rk*(x2-x1)
  do i=1,m
    z = cos(3.141592654_rk*(i-0.25_rk)/(n+0.5_rk))
    z1 = 0.0
    do while(abs(z-z1) .gt. eps)
      p1 = 1.0_rk
      p2 = 0.0_rk
      do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0_rk*j-1.0_rk)*z*p2-(j-1.0_rk)*p3)/j
      end do
      pp = n*(z*p1-p2)/(z*z-1.0_rk)
      z1 = z
      z = z1 - p1/pp
    end do
    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = (2.0_rk*xl)/((1.0_rk-z*z)*pp*pp)
    w(n+1-i) = w(i)
  end do
end subroutine gaulegf

  !
  !
  !
  subroutine MLinvmat(al,ai,dimen,ierr)
  integer,intent(in)   :: dimen
  real(rk),intent(in)  :: al(dimen,dimen)
  real(rk),intent(out) :: ai(dimen,dimen)
  integer(ik),intent(out) :: ierr
  real(rk)             :: h(dimen),p,q
  integer(ik)          :: i1,i2,k,i,j,k8,k9
      

    ierr = 0
    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        !
        if (abs(p)<small_) then 
          !
          ierr = i
          !
          return
          !
        endif 
        !
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine MLinvmat



  subroutine MLinvmatark(al,ai,dimen,ierr)
  integer,intent(in)   :: dimen
  real(ark),intent(in)  :: al(:,:)
  real(ark),intent(out) :: ai(:,:)
  integer(ik),intent(out) :: ierr
  real(ark)             :: h(dimen),p,q
  integer(ik)           :: i1,i2,k,i,j,k8,k9
      

    ierr = 0
    ai=al
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        !
        if (abs(p)<small_a) then 
          !
          ierr = i
          !
          return
          !
        endif 
        !
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0_ark/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine MLinvmatark



end module molecules

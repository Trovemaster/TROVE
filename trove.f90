!
!  TROVE program: theoretical rovibrational energies
!  version 08.10.2006
!
  module tp_module

    use accuracy
    use fields
    use perturbation
    use symmetry
    use timer
    use moltype, only: intensity 
    use dipole, only: dm_tranint,dm_analysis_density
    use refinement, only : refinement_by_fit,external_expectation_values
    use tran, only : TRconvert_matel_j0_eigen,TRconvert_repres_J0_to_contr

    implicit none

!
! Defining the calculations
!
    !
    integer, parameter    :: verbose  = 2 ! Verbosity level    
    !
    contains
    !
    ! Here we do the TROVE calculations 
    ! 
    subroutine ptmain
      !
      integer(ik) :: NPTorder,Natoms,Nmodes,Npolyads
      integer(ik) :: j
      !type(FLbasissetT),pointer  :: basisset(:)     ! Basis set specifications: range and type

      !
      !type(FLbasissetT)          :: rotbasis            ! Rotational basis set specifications
      !
      ! Begin with constants intitialization
      !
      call TimerStart('TROVE')
      !
      call accuracyInitialize
      !
      if (job%verbose>=4) then 
        write(out,"('spacing around 1.0 is ',d18.8)") spacing(1.0d0)
      endif 
      !
      !
      ! Here we define the molecular structure and potential function  
      !
      ! make sure that size(poten)>= size(pseudo)
      !
      call FLReadInput(NPTorder,Npolyads,Natoms,Nmodes,j)
      !
      ! Rotation will be treated only for a given J
      !
      !j = basisset(0)%range(1)
      !
      call FLsetMolecule
      !
      ! Intensity calculations
      !
      if (intensity%do) then 
         !
         !call FLbsetInit
         !call FLinitilize_Kinetic
         !call FLinitilize_Potential
         !call FLinit_External_field
         call PTinit(NPTorder,Nmodes,Npolyads)
         !
         if (job%rotsym_do) call PT_conctracted_rotational_bset(j)
         !
         call dm_tranint
         !
         return
         !
      endif
      !
      if (action%fitting) then 
         !
         call PTinit(NPTorder,Nmodes,Npolyads)
         !
         !call FLinitilize_Kinetic
         !call FLinitilize_Potential
         !call FLinit_External_field
         !call FLbsetInit
         !
         call refinement_by_fit
         !
         if (job%verbose>=2) then 
           write(out,"(/'The End of TROVE-Fit')") 
         endif 
         !
         return
         !
      endif
      !       
      ! Here we define the kinetic operator functions g_vib,g_rot,g_cor,pseudo and also potential function 
      !
      call FLinitilize_Kinetic
      !
      ! Here we initialize the potential energy field
      !
      if (job%Potential_Simple) then
         call FLinitilize_Potential_II
      else 
         call FLinitilize_Potential
      end if 
      !
      ! Here we initialize the dipole moment field (if needed)
      !
      call FLinit_External_field
      !
      ! Store the expansion coefficients of Hamiltonian
      !
      if (trim(trove%IO_hamiltonian)=='SAVE') call FLcheck_point_Hamiltonian('HAMILTONIAN_SAVE')
      !
      !
      ! Classical analysis of the Hamiltonian
      !
      if (analysis%classical) then 
        !
        call FL_rotation_energy_surface 
        return
        !
      endif
      !
      ! Initialization of the basis set 
      !
      call FLbsetInit
      !       
      ! Here we initialize the PT elements, such as Nclasses, Nspecies, etc
      !
      call PTinit(NPTorder,Nmodes,Npolyads)
      !
      ! Analysis of the density 
      !
      if (analysis%extF) then  
        call external_expectation_values
        return
      endif 
      !
      if (analysis%density) then 
        call dm_analysis_density
        !call PTanalysis_density(j)
        return
      endif
      !
      ! Copy matrix elements from the FIELD to PT modules  
      !
      call PTget_primitive_matelements(j) 
      !
      ! Restoring  the contracted basis set vectors from the check point: 
      !
      if (trim(job%IOcontr_action)=='READ') then
        !
        call PTcheck_point_contracted_space('READ') 
        !
      elseif (action%convert_vibme) then
        !
        call TRconvert_repres_J0_to_contr(j)
        !
      else
        !
        call PTcontracted_prediagonalization(j)
        !
      endif 
      !
      ! The rotational part of the contracted basis set to finish its constraction:
      !
      call PT_conctracted_rotational_bset(j)
      !
      call PTsymmetrization(j)
      !
      !if (trove%DVR) call PTDVR_contracted_bases(j)
      !
      ! This type of diagonalization is obsolete: if (job%global_contract) call PTDiagonalize_hamiltonian_symm(j,nroots)
      ! 
      ! Convert the J=0 basis set and mat.elements to the contracted represent. 
      !
      if (action%convert_vibme) then 
         call TRconvert_matel_j0_eigen(j)
         return 
      endif
      !
      if (job%contrci_me_fast) then 
        !
        call PTstore_contr_matelem(j)
        !
        call PTcontracted_matelem_class_fast(j) 
        !
      else
        !
        call PTcontracted_matelem_class(j)
        !
      endif
      !
      call PThamiltonian_contract(j)
      !
      ! convert to j=0 representation as part of the first step j=0 calculation
      !
      if (job%convert_model_j0) then
        !
        call TRconvert_repres_J0_to_contr(j)
        call TRconvert_matel_j0_eigen(j)
        !
      endif
      !
      call TimerStop('TROVE')
      !
      call MemoryReport
      !
      call TimerReport
      !       
      if (job%verbose>=2) then 
        write(out,"(/'End of TROVE')") 
      endif 
      !
    end subroutine ptmain 


  end module tp_module
  !
  program driver
    !
    use tp_module
    !call ompaffinity()

    call ptmain

  end program driver

!
!  S.N. Yurchenko, s.yurchenko@ucl.ac.uk
!

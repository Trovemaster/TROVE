!
!  Various  fields: kinetic, potential, pseudopotential 
!
module fields
   use accuracy
   use molecules
   use lapack
   use me_str
   use me_bnd, only : ME_box,ME_Fourier,ME_Legendre,ME_Associate_Legendre,ME_sinrho_polynomial,ME_sinrho_polynomial_k,&
                      ME_sinrho_polynomial_k_switch,ME_sinrho_polynomial_muzz,ME_legendre_polynomial_k,&
                      ME_laguerre_k,ME_laguerre_simple_k
   use me_numer
   use me_rot
   use timer
   use moltype
   use symmetry , only : SymmetryInitialize,sym

   ! use perturbation

   implicit none
   private
   public FLHarmonicEnergy,FLQindex,FLsetMolecule,FLinitilize_Kinetic,FLinitilize_Potential_II,FLpolynomT
   public FLbsetInit,FLbasissetT,FLmatrixelements
   public FLpt_orders_distribution,FLReadInput,FLrotation,FLenergy_zero,FLfingerprint
   public FLmatrixelements_single_iterm,FLread_fields_dimensions,FLread_rot_matelem
   public FLIndexQ,FLfree_primitive_objects,FL_exclude_specific_modes
   public FLread_coeff_matelem,FLinitilize_Potential_original
   public FLcalc_poten_kinet_dvr,job,FLcalcsT,FLenercutT,FLeigenfile,FLinitilize_Potential,FLinit_External_field_andrey
   public FLextF_coeffs,FL_rotation_energy_surface,FLextF_matelem,FLread_iorder_send
   public jobt, trove, bset, analysis, action, FLL2_coeffs, FLread_fields_dimension_field,FLread_IndexQ_field
   !
   public BaisSetT,Basis1DT,FL_fdf,FLNmodes,FLanalysisT,FLresT,FLpartfunc,FLactionT,FLfinitediffs,FLpoten_linearized,FLread_ZPE
   !
   public j0fit,fitting,FLfittingT,FLobsT,FLread_extF_rank,FLcoeffs2dT,FLpoten4xi,FLfinitediffs_2d
   public FLcheck_point_Hamiltonian,FLinitilize_Potential_Andrey,FLinit_External_field,FLpoten_linearized_dchi,&
          FLDVR_gmat_dvr,FLfromcartesian2local,FLcoeffprunT,coeffprun
   
   public choose, sum_choose, powers_from_index
   !
   !FLfrom_local2chi_by_fit
   !
!
!   Type to define Z-matrix
!


   !type  FLZmatrixT
   !   character(len=cl)    :: name         ! Identifying name of atom (no effect on anything)
   !   integer(ik)          :: connect(4)   ! z-matrix connections
   !end type FLZmatrixT



!
!   Type that will define all fields: kinetic, potential etc.
!
   type FLpolynomT
      character(len=cl)    :: name         ! Identifying name of the polynom
      integer(ik),pointer  :: iorder(:)    ! distribution of poten-coeffs over different PT-orders, 
                                           ! size(iorder)=size(field)
      integer(ik)          :: Orders       ! Max. expansion order 
      integer(ik)          :: Ncoeff       ! Number of expansion coeffs.
      integer(ik)          :: Npoints      ! Number of expansion centers.
      real(ark),pointer    :: field(:,:)   ! Expansion parameters
      real(ark),pointer    :: me(:,:,:)    ! 1D-numerof type matrix elements 
      !
      !integer(ik)          :: SNterms      ! Number of expansion coeffs in the sparse representation
      integer(ik),pointer  :: ifromsparse(:) ! a accounting-index from isparse to icoeff 
      !integer(ik),pointer  :: itosparse(:) ! a accounting-index from isparse to icoeff 
      integer(ik),pointer  :: IndexQ(:,:)    ! This is to store FLIndexQ for each object individually 
      logical :: sparse = .false.            ! indicates if the field in the sparse representation 
      !
   end type FLpolynomT

!

   type  FLenercutT 
     !
     real(rk)  :: general        ! Energy cutoff for everything
     real(rk)  :: primt          ! Energy cutoff for primitives
     real(rk)  :: contr         ! Energy cutoff for contratced basis functions 
     real(rk)  :: matelem         !  Energy cutoff for contratced matrix elements in the matelem section
     real(rk)  ::  DeltaE=0      ! Energgy cutoff based on teh difference between the row-Energy and column-Energy
     integer(ik) :: polyad     ! Polyad cutoff for primit. functs. in contratced basis functions 
     !
   end type  FLenercutT  
   !
   type FLcoeffs2dT                      
      integer(ik),pointer :: icoeff(:,:)
      real(rk),pointer    :: fcoeff(:)
   end type FLcoeffs2dT
   !
   ! files with the eigenvectors 
   !
   type  FLeigenfile
     !
     character(len=cl)  :: filebase
     character(len=cl)  :: dscr       ! file with fingeprints and descriptions of each levels + energy values
     character(len=cl)  :: primitives ! file with the primitive quantum numbres   
     character(len=cl)  :: vectors    ! eigenvectors stored here 
     character(len=cl)  :: dvr        ! eigenvectors in dvr representation stored here 
     !
   end type  FLeigenfile 


   type  FLbasissetT 
     character(len=cl)          :: type         ! Identifying type    of the basis functions 
     character(len=cl)          :: coord_kinet  ! Identifying type    of the basis functions 
     character(len=cl)          :: coord_poten  ! Identifying type    of the basis functions 
     integer                    :: model        ! Applied for contraction of basis sets, reduce or mormal model will be used
     character(len=2 )          :: dim          ! Identifying dimensionality of the basis 
     integer(ik)                :: species      ! Identifying the spicies 
     integer(ik)                :: class        ! Identifying the class  
     integer(ik)                :: range(2)     ! Index range for the given mode 
     real(ark)                  :: res_coeffs   ! contribution of each quanta into a polyd 
     integer(ik)                :: npoints      ! optional integer, e.g. npoints for numerov
     real(ark)                  :: borders(2)   ! optional real, integration borders for numerov 
     logical                    :: periodic     ! change to the periodic boundary condition
     integer(ik)                :: iperiod      ! the period for the periodic tretament 
     character(len=cl)          :: dvr          ! Identifying type    of the dvr representation
     integer(ik)                :: dvrpoints    ! number of dvr-integration points for each mode 
     logical                    :: postprocess  ! Post-diagonalization of the contracted basis set
     logical                    :: Lvib         ! Using the angualr vibrational momentum for symmetrization and building contracted classes
     logical                    :: check_sym    ! check that the corresponding 1D Hamiltonians from a class are identical in 
   end type  FLbasissetT 

   type FLcoeffprunT
    !
    real(rk)           :: contribution_threshold = 1e-4_rk
    !
   end type FLcoeffprunT




!
! We collect all relevant data of the molecule under MoleculeT - type structure 
!
   type JobT
      character(len=cl)   :: Moltype        ! Identifying type of the Molecule (e.g. XY3)
      character(len=cl),pointer :: Coordinates(:,:) ! Identifying the coordinate system, e.g. 'normal'; kinetic and potential part can be different
      character(len=cl)         :: internal_coords  ! type of the internal coordinates: linear,harmonic,local
      character(len=cl)         :: coords_transform ! type of the coordinate transformation: linear,morse-xi-s-rho,...
      character(len=cl)         :: symmetry         ! molecular symmetry
      integer(ik)::  Natoms               ! Number of atoms
      integer(ik)::  Nmodes               ! Number of modes = 3*Natoms-6
      integer(ik)::  Nmodes_e             ! Effective number of modes: when rank = 1 Nmodes_e = Nmodes-1
      integer(ik)::  Nmodes_n             ! non-rigid motion mode, given numerically
      integer(ik),pointer ::  manifold_rank(:) ! The size of the expansion center manifold: 0 for 0d, "1" for 1d
      integer(ik)::  Nbonds               ! Number of bonds
      integer(ik)::  NAngles              ! Number of angles
      integer(ik)::  NDihedrals           ! Number of dihedral angles type 1 and 2
      integer(ik)::  Ncoords              ! Total number of gdc (Nbonds+Nangles+Ndihedrals)
      integer(ik)::  Npoints              ! Number of points of the 1d numerical representation in case of rank=1
      integer(ik),pointer ::  bonds(:,:)  ! Bond connections
      integer(ik),pointer ::  angles(:,:) ! Angles connections
      integer(ik),pointer ::  dihedrals(:,:)  ! Dihedral Angles connections type 1
      integer(ik),pointer ::  dihedtype(:)  ! Dihedral Angles connections type 1
      type(MLZmatrixT),pointer  :: zmatrix(:) ! 
      integer(ik)::  jmax                 ! Angular momentum quantum number - maximal value 
      !
      real(ark),pointer   ::  a0(:,:)      ! Cartesian coordinates at the equilibrium  x_Na = a0(N,a),
                                          ! a = x,y,z (1,2,3)
      real(ark),pointer   ::  b0(:,:,:)    ! General Cartesian coordinates at the equilibrium  x_Na = b0(N,a,i),
                                          ! a = x,y,z (1,2,3), i = 0..Npoints -> in case of manifoldrank=1
      real(ark),pointer   ::  db0(:,:,:,:) ! derivative of b0 wrt rho 
      !
      real(ark),pointer   ::  rho_i(:)     ! Integration points for the Gaussian/Simson integration rules
      real(ark)           ::  rho_ref
      integer(ik)        ::  iPotmin      ! Point of the minimum in the Numerov-integration representation
      !real(rk),pointer   ::  weight_i(:) ! Integration weightd
      !
      real(ark),pointer ::  Amatrho(:,:,:,:) ! A-matrix: x_Na = a0_Na + \sum_l A_Nal xi_l,  for the rank=1 case
      real(ark),pointer ::  dAmatrho(:,:,:,:,:)! dA-matrix: derivatives of Amat wrt rho
      real(ark),pointer ::  Bmatrho(:,:,:,:) ! B-matrix: xi_l  = a0_Na + \sum_{Na} B_lNa (x_Na-a_Na)
      real(ark),pointer ::  dBmatrho(:,:,:,:,:) ! dB-matrix
      real(ark),pointer ::  mass(:)        ! Mass of a nuclear
      real(ark),pointer ::  req(:)         ! Equilibrium bond lengths
      real(ark),pointer ::  alphaeq(:)     ! Equilibrium angles
      real(ark),pointer ::  taueq(:)       ! Equilibrium dihedral angles
      real(ark),pointer ::  chi_ref(:,:)   ! Reference values for the internal coordinates chi 
      real(ark),pointer ::  chi_eq(:)      ! Equilibrium values for the internal coordinates chi 
      real(ark),pointer ::  chi0_ref(:)    ! reference-equilibrium values for the internal coordinates chi 
      real(ark),pointer ::  local_eq(:)    ! Equilibrium values for the internal coordinates
      real(ark),pointer ::  local_ref(:)   ! reference values for the internal coordinates
      real(ark),pointer ::  specparam(:)   ! Special parameters for special cases. For example, for amorse
      real(ark),pointer ::  fdstep(:)      ! finite difference element for the numerical differentiation for every mode 
      !real(rk),pointer ::  coordtransform(:,:) ! Linear transformation of the linearized local coordinates 
      integer(ik):: NKinOrder     ! Max order in the kinetic   energy expansion
      integer(ik):: NPotOrder     ! Max order in the potential energy expansion
      integer(ik):: NExtOrder     ! Max order in the external function expansion

      real(rk),pointer:: PotPolyad(:)  ! polyad coefficients for the the potential energy expansion

      logical :: DVR = .false.    ! Process the matrix elements by the DVR (Gaussian Quadtrature) grid-integration 
      logical :: FBR = .true.     ! Process the matrix elements by the DVR (Gaussian Quadtrature) grid-integration 

      logical :: smolyak = .false.  ! Smolyak Gaussian Quadtratures
      character(len=cl) :: smolyak_rule = 'LG'

      integer(ik) ::  Ncoeff             ! Number of parameters in fields-arrays (determined from MaxOrder)
      integer(ik) ::  MaxOrder           ! Max(NKinOrder,NPotOrder)
      integer(ik),pointer ::  RangeOrder(:)      ! polynomials ranges at different orders: fields(RangeOrder(N-1)..RangeOrder(N))
      type(FLpolynomT),pointer ::  poten         ! Potential parameters
      type(FLpolynomT),pointer ::  g_vib(:,:)    ! Vibrational part of Kinetic factor G
      type(FLpolynomT),pointer ::  g_rot(:,:)    ! Rotational  part of Kinetic factor G
      type(FLpolynomT),pointer ::  g_cor(:,:)    ! Coriolis part of Kinetic factor G
      type(FLpolynomT),pointer ::  pseudo        ! pseudo-potential part
      type(FLpolynomT),pointer ::  L2_vib(:,:)   ! Vibrational angular momentum L2
      real(ark),pointer        ::  imat_s(:)        ! singular compnent of moments of inertia
      !
      type(FLpolynomT),pointer ::  extF(:)       ! External field
      !
      character(len=cl)        ::  potentype    ='GENERAL'  ! Type of the potential energy function representaion
      character(len=cl)        ::  kinetic_type ='GENERAL'  ! Type of the kinetic energy representaion
      real(ark),pointer        ::  qwforce(:,:,:)  ! qwadratic force constants f(l,m) in case of rank=1
      real(ark),pointer        ::  omega(:)      ! Harmonic frequencies 
      real(ark),pointer        ::  coord_f(:)      ! conversion factor to the standard coordinate (normal, morse ...)
      !
      integer(ik)              ::  lincoord = 0  ! degenerated cartesian coord. x,y, or z (1,2,3) or none (0)
                                                 !  in case of linear molecules 
      real(ark)                :: rho_border(2)  ! rhomim, rhomax - borders
      real(ark)                :: rhostep        ! step size
      integer(ik)              :: numerpoints = -1   ! step size
      integer(ik)              :: iothreads = 1          ! Number of threads used for IO
      logical                  :: periodic       ! periodic boundary condition
      !
      logical                  :: sing_at_rho_0  !  we use this parameter to check whether 
                                                   ! if some functions are infinity at rho=0  and the basis function has to have a singularity at rho=0
      character(len=cl)   :: IO_hamiltonian = 'NONE'  ! we can either SAVE or READ the Hamiltonian objects
      character(len=cl)   :: IO_potential   = 'NONE'  ! we can either SAVE or READ the potential objects
      character(len=cl)   :: IO_kinetic     = 'NONE'  ! we can either SAVE or READ the kinetic objects
      character(len=cl)   :: IO_basisset    = 'NONE'  ! we can either SAVE or READ the primitive basis set
      character(len=cl)   :: IO_ext_coeff   = 'NONE'  ! we can either SAVE or READ the potential objects
      character(len=cl)   :: IO_contrCI   = 'NONE'  ! we can either SAVE or READ the CI contracted matrix elemenents by classes
      character(len=cl)   :: IO_primitive_hamiltonian = 'NONE' ! we can either SAVE or READ the primive matrix elements
      character(len=cl)   :: chk_fname = 'prim_bset.chk'   ! file name to store the primitive basis set data 
      character(len=cl)   :: chk_numerov_fname = 'numerov_bset.chk'   ! file name to store the numerov basis set data 
      character(len=cl)   :: chk_hamil_fname   = 'hamiltonian.chk'   ! file name to store the primitive basis set data 
      character(len=cl)   :: chk_poten_fname   = 'potential.chk'     ! file name to store the expansion parameters 
      character(len=cl)   :: chk_kinet_fname   = 'kinetic.chk'       ! file name to store the expansion parameters 
      character(len=cl)   :: chk_external_fname   = 'external.chk'       ! file name to store the expansion parameters 
      character(len=cl)   :: chk_contrib_fname = 'contribution.chk'  ! file name to store the highest contribution to eigenfunctions of contracted basis set 
      !
      logical             :: separate_store = .false.             ! if want to store the Hamiltonian chk also into separate files
      logical             :: separate_convert  = .false.          ! convert hamiltonian.chk to potential.chk and kinetic.chk
      logical             :: checkpoint_iorder  = .false.
      logical             :: sparse  = .false.                    ! A sparse representation of fields 
      logical             :: triatom_sing_resolve = .false.
      integer(ik)         :: krot = 0  ! The value of the krot quantum number (reference or maximal) to generate non-rigid basis sets
      integer(ik)         :: kmax = 0  ! The value of the kmax quantum number (maximal) to generate non-rigid basis sets
      !
   end type JobT
   !
   type FLpartfunc
     !
     real(rk)   :: value      ! The value of the partition function
     real(rk)   :: temperature = 300._rk ! The temperature for the partition function in the matexp mode
     real(rk)   :: zpe         = 0       ! The ZPE for the partition function in the matexp mode
     real(rk),pointer   :: gns(:)        ! The stat. weights for the partition function in the matexp mode
     !
   end type FLpartfunc
   !
   type FLcalcsT
      !
      real(rk)            :: PTthreshold    ! with PT threshold we control if PT is valid, i.e. F/Delta_ener << 1
      real(rk)            :: enercut        ! energy cut for the basis set 
      real(rk)            :: potencut=1e6   ! potential energy cut for the dvr-grid points
      real(ark)           :: ZPE = -epsilon(1.0_rk)  ! Zero point energy 
      logical,pointer     :: isym_do(:)     ! process or not the symmetry in question
      !
      logical             :: sym_C  = .false.   ! if symmetry = C 
      !
      real(rk)            :: erange(2) =(/-0.1,huge(1.0)/) !  energy range 
      type(FLenercutT)    :: enercutoff     ! energy cut for the basis set 
      integer(ik)         :: pot_pt_shift    
      integer(ik)         :: Npolyads_prim  ! maximal polyad number for the contraction
      integer(ik)         :: stored_size = 0  ! dimension of the vector to be stored to RAM
      integer(ik)         :: swap_size = 1e4  ! dimension of the vector to be swaped to disk 
      real(rk)            :: compress = 1.0_rk ! the commression factor for the contracted basis functions in the mat. elem. calc.
      integer(ik)         :: Npolyads_contr ! maximal polyad number for the contraction
      real(rk)            :: coeff_thresh   ! primitve bs-function threshold to exclude quantum with small coeffs
      real(rk)            :: exp_coeff_thresh = 1e-24 ! threshold to remove vanishing expansion coefficients of the Hamiltonian
      real(rk)            :: zeroerror      ! allowed error in the zero order solution 
      integer(ik)         :: iwork          ! maximal allowed size for the matrix to be diagonalaized
      real(rk)            :: factor = 1     ! factor for nroots defining the vector space in seudv diagonalization
      real(rk)            :: pt_ener_thresh = -1  ! energy threshold for the PT corrections, in 1/(Ei-Ej), where |Ei-Ej|<pt_ener_thresh
      real(rk)            :: ener_thresh = -1  ! energy threshold for the basis set reduction using contributions to the primitive energies estimated from the matrix elements by PT
      integer(ik)         :: maxiter = 1000 ! maximal number of iterations in arpack 
      real(rk)            :: tolerance = 0  ! tolerance for arpack diagonalization, 0 means the machine accuracy
      logical             :: restart=.false. ! restart the eigensolutions using the stored eigenvectors of a lower dimension
      integer(ik)         :: PTDeltaQuanta  ! Every new Perturb. order brings increment of the the MaxPolyad(iorder)
      character(len=cl)   :: PTtype         ! Defines type of the perturbationa theory - standard or diagonal
      character(len=cl)   :: diagonalizer
      character(len=cl)   :: mat_readwrite
      character(len=cl)   :: orthogonalizer  = 'GRAM-SCHMIDT'
      real(rk)            :: cluster = -1     ! the factor to porepare clustering of the basis set
      integer(ik)         :: swap_after       ! The Hamiltonian will be stoted, read, and transformed vector by vector
      real(rk)            :: upper_ener       ! upper energy limit for the eigenvalues to be found in variational diagonalization with syevr
      real(rk)            :: thresh           ! thresh of general use
      real(rk)            :: max_swap_size              ! maximal allowed disk space 
      character(len=cl)   :: IOeigen_action = 'NONE'   ! we can either SAVE or READ the eigenfunctions from an external file
      character(len=cl)   :: IOcontr_action = 'NONE'   ! we can either SAVE or READ the contracted eigenfunctions from an external file
      character(len=cl)   :: IOcontr_ = 'NONE'         ! we can either SAVE or READ the incoplete contracted eigenfunctions in DVR representaion from an external file
      character(len=cl)   :: IOkinet_action = 'NONE'   ! we can either SAVE or READ the vib. matrix elem. of the kinetik rot. part from an external file
      character(len=cl)   :: IOextF_action = 'NONE'    ! we can either SAVE or READ the vib. matrix elem. of the kinetik rot. part from an external file
      character(len=cl)   :: IOfitpot_action = 'NONE'  ! we can either SAVE or READ the vib. matrix elem. of the kinetik rot. part from an external file
      character(len=cl)   :: IOj0matel_action = 'NONE' ! we can either SAVE or READ the j=0 eigen.vib. matrix elem. of the kinetik rot. part from an external file
      character(len=cl)   :: IOj0ext_action = 'NONE'   ! we can either SAVE or READ the j=0 eigen.vib. matrix elem. of the kinetik rot. part from an external file
      character(len=cl)   :: IOj0contr_action = 'NONE' ! convert the contracted basis set into the j0-representaion
      character(len=cl)   :: IOswap_matelem = 'NONE'   ! we can swap (using divide) or join the contracted mat. elements from an external file
      character(len=cl)   :: IOdvr_prim = 'NONE'       ! we store and read the DVR primitive objects 
      logical             :: IOvector_symm = .true.    ! store eigen-vectors in the symmetrized iired. representation 
      logical             :: IOeigen_compress = .false.   ! compress the computed eigenvectors using the threshold factor coeff_thresh
      logical             :: IOmatelem_divide = .false.   ! divide the matelem checkpoint into pieces 
      logical             :: IOmatelem_stitch = .false.   ! stitch the matelem checkpoints from split pieces 
      logical             :: IOmatelem_split  = .false.   ! split the matelem checkpoint into pieces (can be different from divide)
      logical             :: IOExtF_divide = .false.      ! divide the ExtF checkpoint into pieces 
      logical             :: IOextF_stitch = .false.      ! stitch the ExF part during the J=0 conversion
      logical             :: IOfitpot_divide = .false.    ! divide the ExtF checkpoint into pieces 
      character(len=cl)   :: compress_file  = 'compress'  ! the file name for storing the compressed eigenvectors will start with this character
      logical             :: matelem_append               ! append the matrix elements after the record (basis function) 
      logical             :: IOmatelem_dump               ! dump into a temperal file an for append later 
      integer(ik)         :: iappend                      ! Record in the matelem matrix to append after
      logical             :: extmatelem_append               ! append the matrix elements after the record (basis function) 
      logical             :: IOextmatelem_dump               ! dump into a temperal file an for append later 
      integer(ik)         :: iextappend                      ! Record in the matelem matrix to append after
      !
      type(FLeigenfile)   :: eigenfile
      type(FLeigenfile)   :: contrfile
      character(len=cl)   :: kinetmat_file  = 'contr_matelem.chk'  ! file where we store the contracted mat. elemens of the Hamiltonian
      character(len=cl)   :: kinetmat_format = 'NEW'               ! in case an old format of checkpointing is used
      character(len=cl)   :: dvr_chkfile  = 'dvr_prim.chk'  ! file where we store the dvr-objects
      character(len=cl)   :: extFmat_file   = 'contr_extfield.chk'  ! file where we store the external field contr. mat. elemens
      character(len=cl)   :: fitpot_file    = 'fitpot_mat'  ! file where we store the external field contr. mat. elemens
      character(len=cl)   :: kineteigen_file = 'j0_matelem.chk'  ! file where we store the eigenvec. mat. elemens of the kinet. part
      character(len=cl)   :: exteigen_file   = 'j0_extfield.chk'  ! file where we store the eigenvec. mat. elemens of the external function
      character(len=cl)   :: matrix_file     = 'matrix'  ! file where we store the contracted mat. elemens of the Hamiltonian
      character(len=cl)   :: extmat_suffix   = 'extmatelem'  ! filename suffix  for storing extF matrix elements of the Hamiltonian
      character(len=cl)   :: j0extmat_suffix  = 'j0_extmatelem'  ! filename suffix  for storing j=0 extF matrix elements of the Hamiltonian
      character(len=cl)   :: matelem_suffix   = 'matelem'  ! filename suffix  for storing matrix elements of the Hamiltonian
      character(len=cl)   :: j0matelem_suffix  = 'j0_matelem'  ! filename suffix  for storing j=0 matrix elements of the Hamiltonian
      character(len=cl)   :: tdm_file          = 'j0_tdm'      ! filename name for vibrational j=0 transition dipole moments for replacement 

       
      real(rk)            :: TMcutoff = epsilon(1.0_rk )   ! threshold to select basis set based on the TM or vibrational intensities
      real(rk)            :: TMenermin  = 0     ! Minimal energy to apply the TMcutoff for 
      logical             :: TMpruning = .false. ! TM-prune if TMcutoff > 0
      character(len=cl)   :: solution_file     = 'solution_'  ! file where we store the contracted mat. elemens of the Hamiltonian
      real(rk)            :: degen_threshold         ! a threshold to find degenerate values
      logical             :: test_diag        ! switch on the test (memory demanding) diagonalization of the hamiltonian 
      character(len=cl)   :: sym_group        ! symmetry 
      integer(ik)         :: verbose          ! soft verbose level (from input, not precompiled)
      logical             :: vib_contract     ! whether the vibrational contraction is used 
                                              ! Vibrational contraction  - we use the computed J=0 basis functions for rovibr. problem)
      logical             :: eigen_contract = .false.   ! whether we contract based on the contracted basis coefficients of the eigenfunctions
      logical             :: fast = .true.    ! fast and exensive calculation of the contracted matrix elements 
      logical             :: vib_rot_contr = .false.    ! the contracted basis is computed using vibration-rotation scheme where 
                                              !  the vibrational indeces run first and K runs last in contrast to the default rot-vib scheme 
                                              !  where K is the first index and the matrix is build as K-blocks. 
      logical             :: sparse = .false. ! to switch on sparse matrix processing
      !
      type(FLbasissetT),pointer  :: bset(:)   ! Basis set specifications: range and type
      real(rk),pointer    :: symm_toler(:) ! tolerance that decides whether the symmetry transformation matrix 
                                              ! has been properly recostracted at the sample point, i.e. the transformed function
                                              ! coincides with the sampe function at the transformed sample point. 
      integer(ik)         :: msample_points = 40  ! number of sample points for determination of the symmetry transformational properties of the contr. solution
      integer(ik)         :: msample_attempts = 100 ! maximal number of attempts of samplings in the symmetry determinations
      type(FLpartfunc)    :: partfunc         ! The parameters for the partition function calcs.
      integer(ik)         :: iswap(2)  = (/1,1000/)    ! compute and checkpoint to the disk the matrix elements for iswap(1) ... iswap(2)
      logical             :: rotsym_do = .false. ! for 'K-BASED' or 'NONE' the rotational symmetry is defined based on the K and tau quanta (default) 
                                                 ! for 'EULER-BASED' the rotational symmetry is defined based on the euler angles transformations
      logical             :: contrci_me_fast = .false.
      integer(ik)         :: MaxVibMomentum_contr    ! maximal L (vibang) for the contraction
      logical		      :: ignore_vectors = .false.
      logical             :: convert_model_j0   = .false. ! convert to J=0 representation as part of the 1st step J=0 calculation
      logical             :: exomol_format  = .false.     ! exomol format of intensity output 
      logical :: Potential_Simple = .false. ! This is simple finite differences type of the potential expansion
                                            ! the default is to exand to N+1 and set the N+1 terms to zero, which is more accurate
      logical             :: triatom_sing_resolve = .false.
      logical,pointer     :: select_gamma(:)! the diagonalization will be done only for selected gamma's
      integer(ik),pointer :: nroots(:) ! number of the roots to be found in variational diagonalization with syevr
      integer(ik)         :: lincoord=0 ! a singularity axis 1,2,3 if present, otherwise 0 
      !
   end type FLcalcsT



   type FLobsT
     !
     integer(ik) :: Jrot
     integer(ik) :: symmetry
     integer(ik) :: N
     real(rk)    :: energy
     real(rk)    :: weight
     integer(ik),pointer :: quanta(:)
     !
   end type FLobsT


   type FLfittingT
     !
     logical              :: run
     integer(ik)          :: j_list(1:100) = -1
     integer(ik)          :: iparam(1:2) = (/1,1000000/)
     integer(ik)          :: itermax = 500
     integer(ik)          :: Nenergies = 1
     real(rk)             :: factor = 1.0_rk
     real(rk)             :: target_rms = 1e-8
     real(rk)             :: robust = 0
     real(rk)             :: watson = 1.0d0
     character(len=cl)    :: method = 'FAST'
     character(len=cl)    :: geom_file = 'pot.fit'
     character(len=cl)    :: output_file = 'fitting'
     character(len=cl)    :: fit_type = 'DGELSS'      ! to switch between fitting methods. 
     real(rk)             :: threshold_coeff = -1e-18
     real(rk)             :: threshold_lock  = -1e-18
     real(rk)             :: threshold_obs_calc  = -1e-16
     real(rk)             :: fit_scale  = 0.4
     type(FLobsT),pointer :: obs(:)
     !
   end type FLfittingT



   type FLresT
     !
     integer(ik) :: Ntheta = 0
     integer(ik) :: Nphi =0 
     integer(ik) :: Ntau = 0
     real(rk)    :: theta = 0
     real(rk)    :: phi = 0
     real(rk)    :: tau = 0
     real(rk)    :: theta1 = 0
     integer(ik) :: Itermax = 200
     character(len=cl) :: type = 'GLOBAL_SEARCH'      ! to switch between fitting methods. 
     !
   end type FLresT

   type FLanalysisT
      !
      logical             :: density = .false.
      logical             :: classical = .false.
      logical             :: rotation_energy_surface = .false.
      logical             :: rotation_density = .false.
      type(FLresT)        :: res
      logical             :: reduced_density = .false.
      logical             :: print_vector = .false.
      logical             :: rotation_matrix = .false.
      logical             :: extF = .false.
      logical             :: check_Hamiltonian = .false.
      integer(ik)         :: dens_list(1:100) = -1      ! List of eigenvalues for the reduced density analysis 
      integer(ik)         :: j_list(1:100) = -1
      integer(ik)         :: sym_list(1:100) = -1
      real(ark)           :: threshold = 1e-8     ! threshold to print out eige-coefficients 
      logical             :: reducible_eigen_contribution = .false. 
      !
   end type FLanalysisT
   !
   !
   type FLactionT
     !
     logical :: fitting       = .false.
     logical :: band_fitting  = .false.
     logical :: convert_vibme = .false.
     logical :: intensity     = .false.
     !
   end type FLactionT



!
!--------------------------------------------------------------
! Bsis set section 
!--------------------------------------------------------------
!


!
!  1D basis set Type
!
   type Basis1DT
      integer(ik)         :: imodes     ! How many modes under this type
      integer(ik),pointer :: mode(:)    ! Which modes under the current type (1..trove%Nmodes)
      integer(ik)         :: Size       ! Size of the 1D basis set 
      integer(ik)         :: Order      ! X**order : max order of magnitude 
      character(len=cl)   :: name       ! Identifying name of the basis functions 
      character(len=cl)   :: type       ! Identifying type of the basis functions 
      real(ark)           :: params(3)  ! few useful parameters, e.g. coef_norm (harmonic), or De/a (morse)
      real(ark), pointer  :: ener0(:)   ! Zero-order energy 
      real(ark), pointer  :: matelements   (:,:,:,:)  ! Matrix elemens <a|x**k*p**l|b>; or (l,k,a,b)
                                                      ! l - 0:2 
                                                      ! k - 0:Order
                                                      ! a,b - 0:Size
   end type Basis1DT





!
!  general basis type definition
!
   type BaisSetT
      integer(ik)                :: n_bset1D_max   ! Number of 1D basis types requested at initialization 
      integer(ik)                :: n_bset_max     ! Total number of basis types:  n_bset1D_max
      integer(ik)                :: Nbset1D        ! Number of 1D bset
      type(FLbasissetT), pointer :: dscr(:)        ! Initial description of the basis set for every mode
      type(Basis1DT), pointer    :: bs1D(:)        ! Simple 1D bset 
      type(Basis1DT), pointer    :: rot            ! Rotational bset 
      !
   end type BaisSetT
!
!  The Basis set itself
!
   type(BaisSetT) , save    :: bset
   type(JobT), save         :: trove

   integer(ik),allocatable,save :: FLIndexQ(:,:)    ! Forward  Relations between 1D arraz and Modes-Dimension array 
   !integer(ik),allocatable,save :: FLIndexQ_legatee(:,:)    ! Addresses to the previous FLIndexQ(:,:) 
   !                                                         !and helps to calculate matrix elements faster

   !real(rk),allocatable,save    :: mat_legatee(:,:)    ! Stored values of already computed mat. elements 


   integer, parameter       :: verbose     = 2    ! Verbosity level
   integer, parameter       :: difftype    = 2    ! differential type: two points or four points finite differences
   !
   integer(ik), save        :: FLNelements        ! number of all real elements we use to count used memory 
   !
   ! This parameter is for switching on and off the numerical representation for the last chi(xi)=rho variable:
   ! if it is "1", this special variable will be given as a 1D table, as well all functions will be represented 
   ! as expansions around each point rho. Thus the manifold has the rank "1"
   ! If it is "0" - all functions are given by the expansions around the 0d manifold, i.e. one point
   !
   character(len=cl),parameter :: axis_system = 'Eckart'
   !
   !
   ! If we solve only vibrational problem, there is no point to keep around the rotational and Coriolis 
   ! kinetic parts. We turn them off by the following logical parameter
   !
   logical :: FLrotation = .false.
   !
   ! External field matrix elements calculations: default values
   !
   logical :: FLextF_matelem = .false.
   logical :: FLextF_coeffs = .false.
   logical :: FLL2_coeffs = .false.
   !
   ! check and fix the discontiunity of non-rigid fields
   logical :: FL_iron_field_out = .false.
   !
   ! current  rho-point 
   integer(ik),save              :: FLirho
   !
   ! potential function is speciafied and stored here 
   !
   real(ark),allocatable :: force(:)
   character(len=16),allocatable :: forcename(:)
   integer(ik),allocatable :: pot_ind(:,:)
   integer(ik),allocatable :: ifit(:)
   logical :: bset_initialized = .false.    !  defaul value is false - to make sure that we do not use the basis set before it has been initialized
   !
   type(FLcalcsT),save           :: job
   type(FLanalysisT),save        :: analysis
   type(FLfittingT),save         :: fitting     ! objects defining the fitting to the observed energies
   type(FLfittingT),save         :: j0fit       ! objects defining the refinement of the band centers
   type(FLcoeffprunT)            :: coeffprun   ! object for pruning basis funcitons 
   !
   ! Andrey's temporale measure: it is not adviceble to use an enviroment variable for that. 
   ! it can bbe very bad for the parallelization. 
   ! 
   integer(ik)              :: FLcurr_imu
   integer(ik)              :: FLNmodes    ! Number of modes
 
   type(FLactionT),save               :: action   ! defines dfifferent actions to perform

   real(ark) :: fd_step_Bmat=1e-4  ! finite difference parameter used for Bmat differentiation numerically 
   real(rk)  :: symm_toler_defaul = 1e-8  ! defaul value for symm_toler 

  contains

!
! Routine to read an input file
!
  subroutine FLReadInput(NPTorder,Npolyads,Natoms,Nmodes,Jrot)
   !
   use input
   !
   integer(ik),intent(out) :: NPTorder,Npolyads,Nmodes,Jrot
   !
   ! parameters used with the "syevr" diagonalizer  
   !
   ! Here we define default values of input parameters 
   !
   integer(ik), parameter :: NPTorder_ =  0     ! Max Perturbation order 
   integer(ik), parameter :: NKinOrder_ = 2     ! Max order in the kinetic   energy expansion
   integer(ik), parameter :: NPotOrder_ = 4     ! Max order in the potential energy expansion
   integer(ik), parameter :: NExtOrder_ = 4     ! Max order in the external function expansion
   !
   ! define the molecule 
   !
   integer(ik), parameter :: Natoms_= 0       ! Number of atoms
   integer(ik), parameter :: Nmodes_= 0       ! Number of modes = 3*Natoms-6
   !
    character(len=cl),parameter  :: Moltype_ ='ABCDEFG'   ! Identifying type of the Molecule (e.g. XY3)
   !
   real(rk), parameter :: fdstep_ = 0.005   ! finite difference element for the numerical differentiation
   !
   ! Perturbation theory parameters 
   !
   integer(ik), parameter :: PTDeltaQuanta_ = 3 ! Every new Perturb. order brings increment of the the MaxPolyad(iorder)
                                                ! The increasment for harmonic bs-case is 3 quanta vs PTorder. 
                                                ! Taking into account res_coeffs - different quanta have  different 
                                                ! contributions into the polyds number "P" - 
                                                ! we obtain: PTDeltaQuanta = 3*max(res_coeffs(:))
                                                ! standard value for harmonic basis set and PH3-type of polyads (nu_s = nu_b):
                                                ! PTDeltaQuanta = 6 
                                                !
   real(rk), parameter    :: PTthreshold_  = 0.9 ! with PT threshold we control if PT is valid, i.e. F/Delta_ener << 1 
   real(rk), parameter    :: PTzeroerror_  = 1.e-2 ! allowed error in the zero order solution 
   integer(ik), parameter :: pot_pt_shift_ = 2  ! Zero order for the potential function starts from the harmonic approximation (=2) 
   !
   ! Basis set parameters 
   !
   !
   integer(ik),parameter  :: Npolyads_  =   80   ! maximal polyad number, i.e. how many polyads we calculate  
   integer(ik), parameter :: iwork_     = 1      ! maximal size of the matrix for the variational diagonalization
                                                 ! used as parameter NCV required by ARPACK. default value 1 means ncv = 21/10*nroots
   !
   integer(ik),parameter  :: manifold_ = 1       !
   real(rk),parameter :: res_coeffs_ = 1.0 ! This defines the polyads or resonanses: P = sum( res%coeffs(:)*nu(:) )
   character(len=cl),parameter :: coords_ = 'LINEAR' ! default internal coordinates 
   !
   type(FLbasissetT),parameter :: vibbasisset_= FLbasissetT('HARMONIC','HARMONIC','HARMONIC',1000,'1D',100,1,(/0,20/),&
                                                             1.0,0,(/-1.0,1.0/),.false.,0,'NUMEROV-POL',0,.false.,.false.,.true.)
   type(FLbasissetT),parameter :: rotbasisset_= FLbasissetT('JKTAU', 'xxxxxx','xxxxxx',1000,'1D',0,0,(/0,0 /),0.0,0,(/0.0,0.0/),&
                                                            .false.,0,'xxxxxx',0,.false.,.false.,.true.)

   integer(ik),parameter :: nroots_=1e6
   real(rk),parameter    :: uv_syevr_=1e9
   integer(ik),parameter :: swap_after_=20000 ! The Hamiltonian will be stoted, read, and transformed vector by vector
                                              ! if the size exceeds this limit, default value 
   real(rk),parameter    :: max_swap_size_=1e6 ! Maximal limit for the swap file (primitive matrix)
                                              ! if the size exceeds this limit, default value 
   real(rk),parameter    :: enermax_=1e6      ! default value for the energy cut of the variational matrix
   real(rk),parameter    :: coeff_thresh_= -tiny(1.0_rk) ! primitve bs-function threshold to exclude quantum witn small coeffs
   !
   logical :: eof,zmat_defined,basis_defined,equil_defined,pot_defined,symmetry_defined,extF_defined,refer_defined,chk_defined
   logical :: kinetic_defined
   character(len=cl) :: Molecule,pot_coeff_type,exfF_coeff_type,chk_type
   character(len=wl) :: w
   real(rk)    :: lfact,f_t
   integer(ik) :: i,iatom,imode,natoms,alloc,Nparam,iparam,i_t,i_tt
   integer(ik) :: Nbonds,Nangles,Ndihedrals,j,ispecies,imu,iterm,Ncoords,icoords
   character(len=4) :: char_j
   integer :: arg_status, arg_length, arg_unit
   character(:), allocatable :: arg
   character(len=cl) :: my_fmt !format for I/O specification
   !
   !
   ! default values: 
   !
   Natoms = Natoms_
   Nmodes = Nmodes_
   trove%Ncoords = 0 
   Ncoords = 0 
   !
   trove%moltype = Moltype_
   manifold = manifold_
   job%IOeigen_action = 'NONE'
   job%IOcontr_action = 'NONE'
   job%IOkinet_action = 'NONE'
   !
   job%iwork      = iwork_
   !
   NPTorder = NPTorder_
   trove%NKinOrder = NKinOrder_
   trove%NPotOrder = NPotOrder_
   trove%NExtOrder = NExtOrder_
   !
   job%PTDeltaQuanta = PTDeltaQuanta_
   job%PTthreshold   = PTthreshold_
   job%zeroerror   = PTzeroerror_
   job%pot_pt_shift  = pot_pt_shift_
   Npolyads      = Npolyads_

   job%enercut =  enermax_
   job%enercutoff%general = enermax_
   job%enercutoff%contr   = enermax_
   job%enercutoff%primt   = enermax_
   job%enercutoff%matelem = enermax_


   job%coeff_thresh = coeff_thresh_
   job%Npolyads_contr = Npolyads_
   job%Npolyads_prim  = Npolyads_
   job%vib_contract = .false.
   !
   trove%internal_coords  = 'LINEARIZED'
   trove%coords_transform = 'LINEAR'
   job%PTtype = 'POWERS'
   if (manifold==1) job%PTtype = 'DIAGONAL'
   !
   zmat_defined  = .false.
   basis_defined = .false.
   equil_defined = .false.
   refer_defined = .false.
   chk_defined = .false.
   pot_defined   = .false.
   kinetic_defined = .false.
   extF_defined  = .false.
   symmetry_defined = .false.
   Nparam = 1 
   !
   job%diagonalizer = 'SYEV'
   !
   job%swap_after = swap_after_
   job%max_swap_size = max_swap_size_
   job%upper_ener = uv_syevr_
   job%thresh = small_
   trove%symmetry = 'C(M)'
   job%sym_group = trove%symmetry
   job%degen_threshold = 10.0**(-4-(ark/8))
   job%test_diag = .false.
   job%verbose = 2
   !
   job%eigenfile%filebase   = 'eigen'
   job%eigenfile%dscr       = 'eigen_descr'
   job%eigenfile%primitives = 'eigen_quanta'
   job%eigenfile%vectors    = 'eigen_vectors'
   !
   job%contrfile%dscr       = 'contr_descr.chk'
   job%contrfile%primitives = 'contr_quanta.chk'
   job%contrfile%vectors    = 'contr_vectors.chk'
   job%contrfile%dvr        = 'contr_dvr.chk'
   !
   job%dvr_chkfile          = 'dvr_prim.chk'
   !
   j0fit%method = 'SLOW'
   !
   ! Intensity block of initial values 
   !
   intensity%action = 'NONE'
   intensity%temperature = 300.0_rk
   intensity%part_func = 1
   intensity%ZPE = -small_
   intensity%freq_window = (/-0.1,1e6/)
   intensity%erange_low = (/-0.1,1e6/)
   intensity%erange_upp = (/-0.1,1e6/)
   !
   !
   ! MEP
   !
   molec%meptype = ''
   !
   ! default value of the rank of the external field is 3, which stands for the dipole moment function.  
   !
   extF%rank = 3

   !
   ! Rotational quantum -  default values
   !
   Jrot = 0
   !
   arg_unit = 5
   call get_command_argument(1, status=arg_status, length=arg_length)
   if (arg_status == 0) then
     allocate( character(arg_length) :: arg )
     call get_command_argument(1, status=arg_status, value=arg)
     if (arg_status == 0) then
       open(newunit=arg_unit, file=arg, status='old')
     end if
   end if
   !
   call input_options(echo_lines=.true.,error_flag=1, default_unit=arg_unit)
   !
   ! read the general input 
   !
   do
       call read_line(eof) ; if (eof) exit
       call readu(w)
       select case(w)
       case("STOP","FINISH","END")
         exit
       case("")
         !
         !print "(1x)"    !  Echo blank lines
         !
       case ("PTORDER")
         !
         call readi(NPTorder)
         !
       case ("KINORDER")
         !
         call readi(trove%NKinOrder)
         !
       case ("POTORDER")
         !
         call readi(trove%NPotOrder)
         !
       case ("POTPOLYAD","POTPOLYADS")
         !
         if (Nmodes == 0) call report("POTPOLYAD cannot appear before NMODES",.true.)
         !
         call readf(trove%PotPolyad(1))
         !
         trove%PotPolyad(2:Nmodes) = trove%PotPolyad(1)
         !
         do i=2,min(Nitems-1,Nmodes)
            call readf(trove%PotPolyad(i))
         enddo
         !
       case ("POTENTIAL_SIMPLE")
         ! This is simple finite differences type of the potential expansion
         ! the default is to exand to N+1 and set the N+1 terms to zero, which is more accurate
         !
         job%Potential_Simple = .true.
         !
       case ("DVR")
         !
         trove%dvr = .true.
         trove%fbr = .false.
         !
         if (Nitems>1) then 
            call readu(w)
            !
            if (trim(w)=="SMOLYAK") trove%smolyak = .true.
            !
            if (Nitems>2) then 
              call readu(trove%smolyak_rule)
            endif
            !
         endif
         !
       case ("FBR")
         !
         trove%dvr = .false.
         trove%fbr = .true.
         !
       case ("NATOMS")
         !
         call readi(trove%Natoms) 
         !
         Natoms = trove%Natoms
         !
       case ("MEM","MEMORY")
         !
         call readf(memory_limit)
         !
         call readu(w)
         !
         select case(w)
             !
           case default 
             !
             call report("Unexpected argument in MEMORY",.true.)
             !
           case("TB","T")
             !
             memory_limit = memory_limit*1024_rk
             !
           case("GB","G")
             !
             memory_limit = memory_limit
             !
           case("MB","M")
             !
             memory_limit = memory_limit/1024_rk
             !
           case("KB","K")
             !
             memory_limit = memory_limit/1024_rk**2
             !
           case("B")
             !
             memory_limit = memory_limit/1024_rk**3
             !
         end select
         !
       case ("NMODES")
         !
         call readi(trove%Nmodes) 
         Nmodes = trove%Nmodes
         ! 
         ! Enviroment variable to be seen outside the module
         !
         FLNmodes = trove%Nmodes
         !
         ! Allocation of bset
         !
         allocate (job%bset(0:Nmodes),trove%fdstep(Nmodes),stat=alloc)
         if (alloc/=0) then
             write (out,"(' Error ',i9,' trying to allocate bs and fdstep')") alloc
             stop 'FLinput, bs and fdstep - out of memory'
         end if
         !
         allocate (extF%intcoords(1:Nmodes),extF%fdstep(Nmodes),stat=alloc)
         if (alloc/=0) then
            write (out,"(' Error ',i9,' allocating matix extF%intcoords ')") alloc
            stop 'FLinput - cannot allocate extF%intcoords'
         end if
         !
         allocate (job%symm_toler(Nmodes),stat=alloc)
         if (alloc/=0) then
             write (out,"(' Error ',i9,' trying to allocate symm_toler')") alloc
             stop 'FLinput, symm_toler - out of memory'
         end if
         !
         ! defafault values for job%symm_toler
         !
         job%symm_toler = symm_toler_defaul
         !
         ! default values for bset and fdstep 
         !
         job%bset(1:trove%Nmodes) = vibbasisset_
         job%bset(0)              = rotbasisset_
         job%bset(0)%range(1) = 0
         job%bset(0)%range(2) = 0
         !
         trove%fdstep = fdstep_
         extF%fdstep = fdstep_
         !
         ! Default values 
         !
         do i=0,Nmodes
            job%bset(i)%species = i
         enddo
         !
         job%iswap(1) = 1 ; job%iswap(2) =  (3+Nmodes)*3
         !
         allocate (trove%PotPolyad(1:Nmodes),stat=alloc)
         !
         trove%PotPolyad = 1.0_rk
         !
       case ("ENERCUT")
         !
         call readf(job%enercut)
         job%enercutoff%general = job%enercut
         job%enercutoff%contr   = job%enercut
         job%enercutoff%primt   = job%enercut
         job%enercutoff%matelem = job%enercut
         !
       case("PRIMITIVES")
         !
         call readu(w)
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
           !
           case ("NUMERPOINTS")
             !
             call readi(trove%numerpoints)
             !
           case ("NPOLYADS")
             !
             call readi(Npolyads)
             !
             job%Npolyads_prim  = Npolyads
             job%Npolyads_contr = Npolyads
             job%MaxVibMomentum_contr = Npolyads
             !
           case ("ENERCUT")
             !
             call readf(job%enercut)
             job%enercutoff%general = job%enercut
             job%enercutoff%contr   = job%enercut
             job%enercutoff%primt   = job%enercut
             job%enercutoff%matelem = job%enercut
             !
           case ("POTENCUT")
             !
             call readf(job%potencut)
             !
           case default
             !
             call report ("Unrecognized unit name "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (trim(w)/="".and.trim(w)/="END") then 
            !
            write (out,"('FLinput: wrong last line in PRIMITIVES =',a)") trim(w)
            stop 'FLinput - illegal last line in PRIMITIVES'
            !
         endif 
         !
       case("CONTRACTION")
         !
         call readu(w)
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
             !
           case ("ENERCUT")
             !
             call readf(job%enercutoff%contr)
             !
             job%erange(2) = job%enercutoff%contr
             !
           case ("ENERCUT_PRIMITIVE","ENERCUT_PRIMIT","ENERCUT_PRIM")
             !
             call readf(job%enercutoff%primt)
             !
           case ("ENERCUT_MATELEM","ENERMAX_MATELEM")
             !
             call readf(job%enercutoff%matelem)
             !
           case ("ENERCUT_DELTA","ENERMAX_DELTA")
             !
             call readf(job%enercutoff%DeltaE)
             !
           case ("SYMM_TOLER")
             !
             !call readf(job%symm_toler)
             !
             if (nitems-1==1) then 
                call readf(f_t)
                job%symm_toler(:) = f_t
                call read_line(eof) ; if (eof) exit
                call readu(w)
                cycle 
             endif 
             !
             if (nitems-1/=Nmodes) then 
                write (out,"('FLinput: wrong number elements in job%symm_toler : ',i8)") nitems-1
                stop 'FLinput - illigal number of job%symm_toler'
             endif 
             !
             do i =1,Nmodes
                !
                call readf(job%symm_toler(i))
                !
             end do
             !
           case ("SAMPLE_POINTS")
             !
             call readi(job%msample_points)
             !
           case ("SAMPLE_ATTEMPTS")
             !
             call readi(job%msample_attempts)
             !
           case ("COEFF_THRESH","COEFF-THRESH")
             !
             call readf(job%coeff_thresh)
             !
           case ("EXP_COEFF_THRESH")
             !
             call readf(job%exp_coeff_thresh)
             !
           case ("ENERGY-THRESHOLD","ENERGY_THRESH","ENERGY-THRESH")
             !
             call readf(job%ener_thresh)
             !
           case ("DEGEN_THRESH","DEGENERACY")
             !
             call readf(job%degen_threshold)
             !
           case ("NPOLYADS")
             !
             call readi(job%Npolyads_contr)
             !
           case ("N_VIBMOMENT","NVIBMOMENT")
             !
             call readi(job%MaxVibMomentum_contr)
             !
           case ("NPOLYADS_PRIM","NPOLYADS_PRIMITIVE")
             !
             call readi(job%Npolyads_prim)
             !
           case ("VIBRATIONAL")
             !
             job%vib_contract = .true.
             !
           case ("EIGEN_COEFF","EIGEN_PRUNING")
             ! 
             job%eigen_contract = .true.
             !
           case("EIGEN_PRUNING_THRES")
             !
             call readf(coeffprun%contribution_threshold)
             !
           case ("FAST_CI","FAST","FAST-CI")
             !
             job%contrci_me_fast = .true.
             trove%IO_contrCI = "NONE"
             !
           case ("VIBINTENSITY","TM","TM_CUTOFF")
             !
             call readf(job%TMcutoff)
             if (job%TMcutoff>0.0_ark) job%TMpruning = .true.
             !
           case ("TM_ENERMIN")
             !
             call readf(job%TMenermin)
             !
           case ("TM_PRINUNG","TMPRINUNG")
             !
             job%TMpruning = .true.
             !
           case ("CLUSTER")
             !
             call readf(job%cluster)
             !
           case ("VIB-ROT")
             !
             job%vib_rot_contr = .true.
             !
             job%IOmatelem_split = .true.
             job%IOextF_divide  = .true.
             !
             job%iswap(1) = 0
             job%iswap(2) = 1e6
             !
           case('SWAP_SIZE')
             !
             call readi(job%swap_size)
             !
           case('STORED_SIZE','STORED')
             !
             call readi(job%stored_size)
             !
           case('COMPRESSION','COMPRESS')
             !
             call readf(job%compress)
             !
           case ("MODEL")
             !
             call readu(w)
             !
             select case(w)
               !
             case ("NONE")
               !
               job%vib_contract = .false.
               !
             case ("VIBRATIONAL","J=0")
               !
               job%vib_contract = .true.
               job%fast = .false.
               !
             case default
               !
               call report ("Unrecognized unit name "//trim(w),.true.)
               !
             end select
             !
           case('ERANGE')
             !
             call readf(job%erange(2))
             !
           case ("ROTSYM")
             !
             call readu(w)
             !
             select case(w)
               !
             case ("NONE","K-BASED")
               !
               job%rotsym_do = .false.
               !
             case ("EULER","EULER-BASED")
               !
               job%rotsym_do = .true.
               !
             case default
               !
               call report ("Unrecognized unit name "//trim(w),.true.)
               !
             end select
             !
           case default
             !
             call report ("Unrecognized unit name "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (trim(w)/="".and.trim(w)/="END") then 
            !
            write (out,"('FLinput: wrong last line in CONTRACTION =',a)") trim(w)
            stop 'FLinput - illegal last line in CONTRACTION'
            !
         endif 
         !
       case ("VERBOSE")
         !
         call readi(job%verbose)
         !
       case ("PTDELTANU")
         !
         call readi(job%PTDeltaQuanta)
         !
       case ("PTTHRESH")
         !
         call readf(job%PTthreshold)
         !
       case ("PTZERO")
         !
         call readf(job%zeroerror)
         !
       case ("MATRIXTYPE")
         !
         ! obsolete! 
         !
         call readu(w)
         !
       case ("PTPOTSHIFT")
         !
         call readi(job%pot_pt_shift)
         !
       case ("NO-SPARSE")
         !
         trove%sparse = .false.
         !
       case ("SPARSE")
         !
         trove%sparse = .true.
         !
       case ("IRON-OUT","IRONOUT","IRON-FIELD-OUT")
         !
         FL_iron_field_out = .true.
         !
       case ("DIAGONALIZER")
         !
         call readu(w)
         !
         select case (trim(w))
         !
         case('SYEV','SYEVR','SYEVX')
           !
           job%diagonalizer = trim(w)
           !
         case default
           !
           if (trim(w)/="") &
              call report ("Unrecognized unit name "//trim(w),.true.)
           !
         end select
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
           !
           case('SYEV','SYEVR','SYEVD','SYEVX','SEUPD','PSEUPD','SYEV4','SYEVR4','SYEVR-4TO8','SYEV-4TO8','SYEV-KROT',&
                'SYEVR-KROT','SYEV-E-PT','SYEVR-E-PT','SYEV-BS-PT','SYEVR-BS-PT','SYEV-KROT-BS-PT','SYEVR-KROT-BS-PT',&
                'SYEV-KROT-EN-PT','SYEV-KROT-EN','SYEVR-KROT-EN-PT','ULEN','DSTEVX','DSTEVR','DSTEVD','DSTEV','DSTEVX-P',&
                'DSTEVR-P','DSTEVD-P','DSTEV-P','DSYEV-ILP','DSYEVR-ILP','DSYEVD-ILP','DSYEVX-ILP','READ-EIGEN',&
                'SYEV0','PROPACK','READ-EIGEN-GRID','PLASMA_DSYTRDX','ENERGY=DIAGONAL','NO-DIAGONALIZATION','READ-EIGEN-GRID4',&
                'READ-ENERGIES')
             !
             job%diagonalizer = trim(w)
             !
             if (trim(w)=='SEUPD'.or.trim(w)=='PSEUPD'.or.trim(w)=='PROPACK') then 
                job%sparse = .true.
             endif
             !
             if ( w(1:10)=='READ-EIGEN'.and.nitems==2) then 
                !
                call readu(w)
                job%solution_file = trim(w)
                !
             endif
             !
             if (trim(w)=='READ-ENERGIES') then 
                !
                job%ignore_vectors = .true.
                !
             endif
             !
           case('MATEXP','MATEXP_SPARSE','MATEXP-SPARSE')
             !
             job%diagonalizer = w
             !
             if (trim(w)=='MATEXP_SPARSE'.or.trim(w)=='MATEXP-SPARSE') then 
                job%diagonalizer ='MATEXP-SPARSE' 
                job%sparse = .true.
             endif
             !
             if (.not.symmetry_defined) then 
                !
                write (out,"('FLinput: MATEXP cannot appear before symmetry is defined')") 
                stop 'FLinput - MATEXP defined before symmetry'
                !
             endif 
             !
             allocate(job%partfunc%gns(sym%Nrepresen),stat=alloc)
             if (alloc/=0) stop 'FLinput, job*partfunc%gns - out of memory'
             !
             job%partfunc%gns = 1.0_rk
             !
           case('READ','SAVE','STORE','SAVE-LOWER','STORE-LOWER','STORE_CHEAP','READ-LOWER')
             !
             job%mat_readwrite = w
             !
             if ( nitems==2) then 
                call readu(w)
                job%matrix_file = trim(w)
             endif
             !
           case ("GAMMA")
             !
             if (.not.symmetry_defined) then 
                !
                write (out,"('FLinput: keyword GAMMA in CONTRACT cannot appear before symmetry is defined')") 
                stop 'FLinput - symmetry should be defined before GAMMA '
                !
             endif
             !
             if (Nitems>sym%Nrepresen+1.or.Nitems==1) then  
               write (out,"('FLinput: illegal number of irreps in gamma in DIAGONALIZER: ',i7)") Nitems-1
               stop 'FLinput - illegal number of gammas in DIAGONALIZER'
             endif
             !
             i = 0
             job%select_gamma = .false.
             !
             do while (item<Nitems.and.i<size(job%select_gamma))
                !
                i = i + 1
                !
                call readu(w)
                !
                if (trim(w)/="-") then
                  !
                  read(w,*) i_t
                  job%select_gamma(i_t) = .true.
                  !
                else
                  !
                  call readi(i_tt)
                  !
                  do while (i_t<i_tt.and.i<size(job%select_gamma))
                    !
                    i_t = i_t + 1
                    job%select_gamma(i_t) = .true.
                    i = i + 1
                    !
                  enddo
                  i = i - 1
                  !
                endif
                !
             enddo
             !
             !call readi(i_t)
             !
             !do while (i_t/=0.and.i<=size(job%select_gamma))
             !   !
             !   i = i+1
             !   !
             !   job%select_gamma(i_t) = .true.
             !   !
             !   call readi(i_t)
             !   !
             !enddo 
             !
           case ("GRAM-SCHMIDT","SVD","SCHMIDT")
             !
             job%orthogonalizer =trim(w)
             !
           case ("SLOW")
             !
             job%fast = .false.
             !
           case ("FAST")
             !
             job%fast = .true.
             !
           case ("SPARSE")
             !
             job%sparse = .true.
             !
           case ("NCV","IWORK","IFACTOR")
             !
             call readf(job%factor)
             !
           case ("PT-THRESHOLD","PT-THRESH")
             !
             call readf(job%pt_ener_thresh)
             !
           case ("RESTART")
             !
             call readu(w)
             !
             if (trim(w)=='YES') then
               job%restart = .true.
             elseif(trim(w)/='NO') then
               write (out,"('FLinput: wrong RESTART key, can be only YES or NO, NOT ',a)") trim(w)
               stop 'FLinput - wrong RESTART key'
             endif
             !
           case("SWAP_AFTER")
             !
             call readi(job%swap_after)
             !
           case("MAX_SWAP_SIZE")
             !
             call readf(job%max_swap_size)
             !
           case("NROOTS")
             !
             !call readi(job%nroots)
             !
             if (.not.symmetry_defined) then 
                !
                write (out,"('FLinput: keyword nroots cannot appear before symmetry is defined')") 
                stop 'FLinput - symmetry should be defined before NROOTS '
                !
             endif
             !
             if (nitems-1==1) then 
                !
                call readi(i)
                job%nroots(:) = i
                !
             else
               !
               if (nitems-1>sym%Nrepresen) then 
                  !
                  write (out,"('FLinput: too many entries in roots (>sym%Nrepresen): ',i8)") nitems-1
                  stop 'FLinput - illigal number of entries in nroots'
                  !
               endif 
               !
               do i =1,nitems-1
                  !
                  call readi(job%nroots(i))
                  !
               end do
               !
             endif
             !
           case("MAXITER")
             !
             call readi(job%maxiter)
             !
           case("TOLERANCE","TOL")
             !
             call readf(job%tolerance)
             !
           case("UPLIMIT","ENERMAX","ENERCUT")
             !
             call readf(job%upper_ener)
             !
           case("THRESHOLD","THRESH")
             !
             call readf(job%thresh)
             !
           case("TEST")
             !
             job%test_diag = .true.
             !
           case("TEMPERATURE")
             !
             call readf(job%partfunc%temperature)
             !
           case("ZPE")
             !
             call readf(job%partfunc%zpe)
             job%zpe = job%partfunc%zpe
             !
           case("GNS")
             !
             i = 0
             !
             job%isym_do = .false.
             !
             do while (item<Nitems.and.i<sym%Nrepresen)
               !
               i = i + 1
               call readf(job%partfunc%gns(i))
               if (job%partfunc%gns(i)>small_) job%isym_do(i) = .true.
               !
             enddo
             !
           case default
             !
             call report ("Unrecognized unit name "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (trim(w)/="".and.trim(w)/="END") then 
            !
            write (out,"('FLinput: wrong last line in DIAGONALIZER =',a)") trim(w)
            stop 'FLinput - illigal last line in DIAGONALIZER'
            !
         endif 
         !
         if (job%upper_ener/=uv_syevr_.and.any(job%nroots/=nroots_)) then 
             !
             write (out,"('FLinput: contradicted definition in SYEVR:')")
             write (out,"('         both nroots and uplimit defined. ')")
             write (out,"('         the program will disregard uplimit and be based on nroots')")
             job%upper_ener = uv_syevr_
             !
         endif
         !
       case ("DSTEP")
         !
         if (Nmodes==0) then 
            !
            write (out,"('FLinput: DSTEP cannot appear before NMODES')") 
            stop 'FLinput - DSTEP defined befor NMODES'
            !
         endif
         !
         ! If all dsteps are the same - only one number in the input can be given. 
         !
         if (nitems-1==1) then 
            !
            call readf(f_t)
            trove%fdstep(:) = f_t
            cycle 
            !
         endif 
         !
         if (nitems-1/=Nmodes) then 
            !
            write (out,"('FLinput: wrong number elements in dstep : ',i8)") nitems-1
            stop 'FLinput - illigal number of dsteps'
            !
         endif 
         !
         do i =1,Nmodes
            !
            call readf(trove%fdstep(i))
            !
         end do
         !
         extF%fdstep = trove%fdstep
         !
       case("COORDS")
         !
         call readu(w)
         !
         trove%internal_coords = trim(w)
         !
       case("TRANSFORM","XI")
         !
         call readu(w)
         !
         trove%coords_transform =    trim(w)
         !
       case("MOLTYPE")
         !
         call readu(w)
         !
         trove%moltype = trim(w)
         !
       case("IOTHREADS")
         !
         call readi(trove%iothreads)
         !
       case("MOLECULE") ! optional, makes no effect 
         !
         call readu(w)
         !
         Molecule = trim(w)
         !
       case("PTTYPE") 
         !
         call readu(w)
         !
         job%PTtype = trim(w)
         !
       case("SINGULAR-AXIS") 
         !
         call readi(trove%lincoord)
         job%lincoord = trove%lincoord
         !
       case("SYMGROUP","SYMMETRY","SYMM","SYM","SYM_GROUP") 
         !
         call readu(w)
         !
         if (nitems>2) call readi(sym%N)
         !
         trove%symmetry = trim(w)
         job%sym_group = trove%symmetry
         !
         if (trim(job%sym_group)=="C".or.trim(job%sym_group)=="C(M)") job%sym_C = .true.
         !
         ! Initialize the group symmetry 
         !
         call SymmetryInitialize(job%sym_group)
         !
         symmetry_defined = .true.
         !
         allocate(job%isym_do(sym%Nrepresen),stat=alloc)
         if (alloc/=0)  stop 'FLinput, isym_do - out of memory'

         allocate(job%select_gamma(sym%Nrepresen),stat=alloc)
         if (alloc/=0)  stop 'FLinput, select_gamma - out of memory'
         !
         job%select_gamma = .true.
         !
         allocate(job%nroots(sym%Nrepresen),stat=alloc)
         call ArrayStart('nroots:sym',alloc,size(job%nroots),kind(job%nroots))
         !
         job%nroots = nroots_
         !
         job%isym_do = .true.
         !
       case("REFER-CONF","REFER-CONFIGURATION")
         !
         call readu(w)
         !
         select case(w)
         !
         case("NON-RIGID")
           !
           manifold = 1
           !
         case("RIGID")
           !
           manifold = 0
           !
         case default
           !
           call report ("Unrecognized unit name "//trim(w),.true.)
           !
         end select
         !
       case("ZMAT")
         !
         iatom = 0
         !
         if (Natoms==0) then 
            !
            write (out,"('FLinput: ZMAT cannot appear before NATOMS')") 
            stop 'FLinput - ZMAT defined befor NATOMS'
            !
         endif 
         !
         allocate (trove%zmatrix(Natoms),trove%mass(Natoms),stat=alloc)
         if (alloc/=0) then
             write (out,"(' Error ',i9,' trying to allocate zmat')") alloc
             stop 'FLinput, zmat - out of memory'
         end if
         !
         call read_line(eof) ; if (eof) exit
         !
         call readu(w)
         !
         do while (trim(w)/="".and.iatom<Natoms.and.trim(w)/="END")
           !
           iatom = iatom+1
           !
           trove%zmatrix(iatom)%connect(1:4) = 0 
           !
           trove%zmatrix(iatom)%name = trim(w)
           !
           if (nitems-1>5) then 
              !
              write (out,"('FLinput: too many  columns in Z-matrix for atom',i8,': ',i8)") iatom,nitems-1
              stop 'FLinput - too many columns in Z-mat'
              !
           endif 
           !
           i=0
           !
           do while (item < nitems-1)
              !
              i=i+1
              !
              call readi( trove%zmatrix(iatom)%connect(i) )
              !
           end do
           !
           call readf(trove%mass(iatom))
           !
           call read_line(eof) ; if (eof) exit
           !
           call readu(w)
           !
         enddo 
         !
         if (iatom/=Natoms.or.(trim(w)/="".and.trim(w)/="END")) then 
            !
            write (out,"('FLinput: wrong number of rows in Zmat for Natoms =',i8,': ',i8)") Natoms,iatom
            stop 'FLinput - illigal number of rows in Zmat'
            !
         endif
         !
         ! Compute number of coordinates 
         !
         Nbonds  = Natoms-1
         !
         nangles = 0 
         Ndihedrals = 0
         !
         do iatom = 3,Natoms
            !
            nangles = nangles + 1
            !
            if (iatom>=4.or.trove%zmatrix(iatom)%connect(3)/=0) then
               !
               J = trove%zmatrix(iatom)%connect(4)
               !
               select case (J) 
                  !
               case(-1,0)
                  !
                  NAngles = NAngles + 1 
                  !
               case(1)
                  !
                  NAngles = NAngles + 1 
                  Ndihedrals = Ndihedrals + 1
                  !
               case(2,202,302,402) 
                  !
                  Ndihedrals = Ndihedrals + 1
                  !
               case(-2,-202,-302,-402) 
                  !
                  Ndihedrals = Ndihedrals + 1
                  !
               case(3:100)
                  !
                  NAngles = NAngles + 2
                  !
               case(101,103,105) 
                  !
                  ! special case of angles for a linear molecule
                  !
                  NAngles = NAngles - 1
                  Ndihedrals = Ndihedrals + 2
                  !
               case(102,104,106,107,108) 
                  !
                  ! special case of an angle for a linear molecule
                  !
                  NAngles = NAngles - 1
                  Ndihedrals = Ndihedrals + 2
                  !
               end select 
            endif
            !
         enddo
         !
         trove%Ncoords = Nbonds+Nangles+Ndihedrals
         Ncoords = trove%Ncoords
         !
         allocate (trove%local_eq(trove%Ncoords),trove%specparam(trove%Ncoords),extF%geom_ref(trove%Ncoords),&
                   trove%local_ref(trove%Ncoords),stat=alloc)
         if (alloc/=0) then
             write (out,"(' Error ',i9,' trying to allocate local_eq')") alloc
             stop 'FLinput, local_eq - out of memory'
         end if
         !
         extF%geom_ref = 0
         trove%local_eq = 0
         trove%specparam = 0
         !
         zmat_defined = .true.
         !
         ! Basis set section 
         !
       case("BASIS")
         !
         if (Nmodes==0) then 
            !
            write (out,"('FLinput: BASIS cannot appear before NMODES')") 
            stop 'FLinput - BASIS defined befor NMODES'
            !
         endif 
         !
         imode = 0
         !
         call read_line(eof) ; if (eof) exit
         !call readu(w)
         !
         do while (trim(w)/="".and.imode<Nmodes.and.trim(w)/="END")
           !
           call readi(i)
           !
           ! this is for the rotational basis set
           !
           if (i==0) then 
              !
              call readu(w)
              !
              select case(w)
              !
              case("JKTAU")
                !
                job%bset(0)%type = trim(w)
                !
                ! objects from the same classes are pre-diagonalized together
                !
                job%bset(0)%class = 0 
                !
                FLrotation = .true.
                !
              case default
                !
                call report ("Unrecognized unit name "//trim(w),.true.)
                !
              end select
              !
           else
              !
              imode = imode+1
              !
              job%bset(imode)%class = i
              !
              if (job%bset(imode-1)%class-i>1) then 
                write(out,"('FLinput: illegal identification of the mode  :',i8)") i
                write(out,"('         the increament > 1')")
                stop 'llegal identification of the mode'
              endif 
              !
              call readu(w)
              !
              job%bset(imode)%type = trim(w)
              !
              !if (trim(job%bset(imode)%dvr)=='xxxxxx') job%bset(imode)%dvr) = 'NUMEROV-POL'
              !
              call readu(w)
              job%bset(imode)%coord_kinet = trim(w)
              !
              ! default value for the extF-expansion
              !
              extF%intcoords(imode)= job%bset(imode)%coord_kinet
              !
              call readu(w)
              job%bset(imode)%coord_poten = trim(w)
              !
              select case (trim(job%bset(imode)%type)) 
                 !
              case ('NUMEROV','BOX','LAGUERRE','FOURIER','LEGENDRE','SINRHO','LAGUERRE-K') 
                 !
              case default 
                 !
                 job%bset(imode)%coord_kinet = job%bset(imode)%type
                 job%bset(imode)%coord_poten = job%bset(imode)%type
                 !
              end select 
              !
              if (trim(job%bset(imode)%type)=='HARMONIC'.and.trove%dvr) then 
                 !
                 if (trim(job%bset(imode)%dvr)=='NUMEROV-POL') job%bset(imode)%dvr = 'HERMITE'
                 !
              endif
              !
              if (trim(job%bset(imode)%type)=='LEGENDRE'.or.trim(job%bset(imode)%type)=='SINRHO'.or.&
                  trim(job%bset(imode)%type)=='LAGUERRE-K') then 
                 trove%triatom_sing_resolve = .true.
                 job%triatom_sing_resolve = .true.
              endif
              !
           endif
           !
           do while (trim(w)/="".and.item<Nitems)
              !
              call readu(w)
              !
              select case(w)
              !
              case("RANGE")
                !
                call readi(job%bset(imode)%range(1))
                call readi(job%bset(imode)%range(2))
                !
                ! in case the range(2) is given for imode=0 we use Jrot as range(2) in order 
                !
                if (imode==0.and.Jrot/=0) then 
                  !
                  job%bset(imode)%range(1) = Jrot
                  job%bset(imode)%range(2) = 2*jrot
                  !
                endif 
                !
                if (job%bset(imode)%dvrpoints==0) job%bset(imode)%dvrpoints = job%bset(imode)%range(2)+1
                !
              case("JROT")
                call readi(Jrot)
                !
                ! we use range(1) to store the Jrot value 
                !
                job%bset(imode)%range(1) = Jrot
                !job%bset(imode)%range(2) = 2*Jrot
                !
              case("KROT","K")
                !
                call readi(i_t)
                !
                ! we use range(1) to store the Jrot value 
                !
                job%bset(imode)%range(2) = i_t
                trove%krot = i_t
                if (trove%kmax==0) trove%kmax = i_t
                !
              case("KMAX")
                !
                call readi(i_t)
                !
                ! we use range(1) to store the Jrot value 
                !
                trove%kmax = i_t
                !
              case("OVRLP","DVRPOINTS","DPOINTS")
                !
                call readi(job%bset(imode)%dvrpoints)
                !
                if (job%bset(imode)%dvrpoints==0) then 
                  write(out,"('illegal number of dvrpoins:',i7)") job%bset(imode)%dvrpoints
                  stop 'input: illegal number of dvrpoins'
                endif 
                !
              case("REDUCED","RED","R","J")
                !
                call readi(job%bset(imode)%model)
                !
              case("DVR")
                !
                call readu(job%bset(imode)%dvr)
                !
              case("RESC","RESCOEF","COEFF")
                !
                call readf(job%bset(imode)%res_coeffs)
                !
              case("POINTS")
                !
                call readi(job%bset(imode)%npoints) 
                !
              case("PERIODIC","PERIOD","P")
                !
                call readi(job%bset(imode)%iperiod)
                !
                if (job%bset(imode)%iperiod>0) job%bset(imode)%periodic =.true.
                !
              case("POST","POSTPROCESS")
                !
                job%bset(imode)%postprocess =.true.
                !
              case("LVIB","VIB_MOMENTUM")
                !
                job%bset(imode)%lvib =.true.
                !
                FLl2_coeffs = .true.
                !
              case("NOCHECK")
                !
                job%bset(imode)%check_sym =.false.
                !
              case("BORDERS")
                !
                call readf(job%bset(imode)%borders(1)) 
                call readf(job%bset(imode)%borders(2))
                !
              case("ANGSTROM")
                !
                lfact= 1.0_rk
                job%bset(imode)%borders(:) = job%bset(imode)%borders(:)*lfact
                !
              case("BOHR")
                !
                lfact=bohr
                job%bset(imode)%borders(:) = job%bset(imode)%borders(:)*lfact
                !
              case("DEG","DEGREE","DEGREES")
                !
                lfact=1.0_rk/rad
                job%bset(imode)%borders(:) = job%bset(imode)%borders(:)*lfact
                !
              case default
                !
                call report ("Unrecognized unit name "//trim(w),.true.)
                !
              end select
              !
           enddo 
           !
           call read_line(eof) ; if (eof) exit
           !
         end do
         !
         ! special case of Assoc Legendre or SINRHO-polynomials 
         !
         if (trim(job%bset(Nmodes)%type)=='LEGENDRE'.or.trim(job%bset(Nmodes)%type)=='SINRHO'.or.&
            trim(job%bset(Nmodes)%type)=='LAGUERRE-K') then 
            !
            if (.not.trove%triatom_sing_resolve ) then
             write(out,"(a)") '   LEGEDRE or SINRHO or LAGUERRE-K types assume a singular triatomic molecule with 3 modes.'
            endif
            !
            job%bset(Nmodes)%range(2) = (job%bset(Nmodes)%range(2)+1)*(job%bset(0)%range(2)+1)-1
            job%bset(Nmodes)%res_coeffs = job%bset(imode)%res_coeffs/real((job%bset(0)%range(2)+1),ark)
         endif 
         !
         call readu(w)
         !
         if (imode/=Nmodes.or.(trim(w)/="".and.trim(w)/="END")) then 
            !
            write (out,"('FLinput: wrong number of rows in Basis for Nmodes =',i8,': ',i8)") Nmodes,imode
            stop 'FLinput - illigal number of rows in BASIS'
            !
         endif 
         !
         ! Introducing equivalent species that will share the same basis set 
         ! 
         ispecies = 0
         job%bset(0)%species = 0
         ! 
         do i=1,Nmodes
            !
            if (job%bset(i)%type       ==job%bset(i-1)%type       .and.job%bset(i)%dim        ==job%bset(i-1)%dim.and.&
                job%bset(i)%coord_kinet==job%bset(i-1)%coord_kinet.and.job%bset(i)%coord_poten==job%bset(i-1)%coord_poten.and.&
                job%bset(i)%class      ==job%bset(i-1)%class      .and.job%bset(i)%dvrpoints  ==job%bset(i-1)%dvrpoints.and.&
                job%bset(i)%range(1)   ==job%bset(i-1)%range(1)   .and.job%bset(i)%range(2)   ==job%bset(i-1)%range(2).and.&
                job%bset(i)%borders(1) ==job%bset(i-1)%borders(1) .and.job%bset(i)%borders(2) ==job%bset(i-1)%borders(2).and.&
                job%bset(i)%res_coeffs ==job%bset(i-1)%res_coeffs .and.job%bset(i)%npoints    ==job%bset(i-1)%npoints .and.&
                (job%bset(i)%periodic.eqv.job%bset(i-1)%periodic)   .and.job%bset(i)%iperiod    ==job%bset(i-1)%iperiod ) then
                !
                job%bset(i)%species = ispecies
                !
            else 
                ispecies = ispecies + 1
                job%bset(i)%species = ispecies
            endif 
            !
         enddo
         !
         basis_defined = .true.
         !
       case("EQUIL","EQUILIBRIUM")
         !
         if (Nmodes==0) then 
            !
            write (out,"('FLinput: EQUIL cannot appear before NMODES')") 
            stop 'FLinput - EQUIL defined befor NMODES'
            !
         endif 
         !
         imode = 0
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.imode<trove%Ncoords.and.trim(w)/="END")
           !
           imode = imode+1
           !
           call readi( i_t  )  ! reading an integer weight, not used so far 
           !
           call readf( trove%local_eq(imode)  )
           !
           if (nitems==4) then
             !
             call readu(w)
             !
             lfact = 1.0_rk
             !
             select case(w)
             !
             case("ANGSTROM")
               !
               lfact= 1.0_rk
               !
             case("BOHR")
               !
               lfact=bohr
               !
             case("DEG","DEGREE","DEGREES")
               !
               lfact=1.0_rk/rad
               !
             case default
               !
               call report ("Unrecognized unit name "//trim(w),.true.)
               !
             end select
             !
             trove%local_eq(imode) = trove%local_eq(imode)*lfact
             !
           endif 
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (imode/=trove%Ncoords.or.(trim(w)/="".and.trim(w)/="END")) then 
            !
            write (out,"('FLinput: wrong number of rows in EQUL for trove%Ncoords =',i8,': ',i8)") trove%Ncoords,imode
            stop 'FLinput - illigal number of rows in EQUIL'
            !
         endif 
         !
         equil_defined = .true.
         !
       case("REFER","REFERENCE")
         !
         if (Nmodes==0) then 
            !
            write (out,"('FLinput: REFER cannot appear before NMODES')") 
            stop 'FLinput - REFER defined befor NMODES'
            !
         endif 
         !
         imode = 0
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.imode<trove%Ncoords.and.trim(w)/="END")
           !
           imode = imode+1
           !
           call readi( i_t  )  ! reading an integer weight, not used so far 
           !
           call readf( trove%local_ref(imode)  )
           !
           if (nitems==4) then
             !
             call readu(w)
             !
             lfact = 1.0_rk
             !
             select case(w)
             !
             case("ANGSTROM")
               !
               lfact= 1.0_rk
               !
             case("BOHR")
               !
               lfact=bohr
               !
             case("DEG","DEGREE","DEGREES")
               !
               lfact=1.0_rk/rad
               !
             case default
               !
               call report ("Unrecognized unit name "//trim(w),.true.)
               !
             end select
             !
             trove%local_ref(imode) = trove%local_ref(imode)*lfact
             !
           endif 
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (imode/=trove%Ncoords.or.(trim(w)/="".and.trim(w)/="END")) then 
            !
            write (out,"('FLinput: wrong number of rows in EQUL for trove%Ncoords =',i8,': ',i8)") trove%Ncoords,imode
            stop 'FLinput - illigal number of rows in EQUIL'
            !
         endif 
         !
         refer_defined = .true.
         !
       case("ZPE")
         !
         call readf(job%zpe)
         !
       case("CHECK_POINT","CHECKPOINT")
         !
         call readu(w)
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         chk_defined = .true.
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
           !
           case('HAMILTONIAN')
             !
             call readu(trove%IO_hamiltonian)
             !
             call check_read_save_none(trove%IO_hamiltonian,w)
             if (trim(trove%IO_hamiltonian)=='SEPARATE'.or.trim(trove%IO_hamiltonian)=='ASCII') then
                trove%separate_store = .true.
             endif 
             !
             if (nitems>=3) then
               !
               call readu(w)
               !
               if (trim(w)=='SEPARATE'.or.trim(w)=='ASCII') then
                  trove%separate_store = .true.
               elseif (trim(w)=='CONVERT') then
                  trove%separate_convert = .true.
                  trove%separate_store = .false.
                  if (trove%separate_store)  trove%IO_hamiltonian = 'SAVE'
               elseif (trim(w)=='CONVERT-BACK') then
                  trove%separate_convert = .true.
                  trove%separate_store = .true.
                  trove%IO_hamiltonian = 'SAVE'
               else
                   call report (" Illegal record after hamilton",.true.)
               endif
               !
             endif 
             !
           case('PRIMITIVE_HAMILTONIAN','PRIM_HAMILTONIAN','PRIM_MATELEM')
             !
             call readu(trove%IO_primitive_hamiltonian)
             !
             call check_read_save_none(trove%IO_primitive_hamiltonian,w)
             !
           case('POTENTIAL','POTEN')
             !
             call readu(trove%IO_potential)
             !
             call check_read_save_none(trove%IO_potential,w)
             !
             if (trim(w)=='SEPARATE') then
                trove%separate_store = .true.
             endif 
             !
             if (trim(trove%IO_potential)=='SAVE') trove%IO_hamiltonian = 'SAVE'
             !
           case('KINET','KINETIC')
             !
             call readu(trove%IO_kinetic)
             !
             call check_read_save_none(trove%IO_potential,w)
             !
             if (trim(w)=='SEPARATE') then
                trove%separate_store = .true.
             endif 
             !
             if (trim(trove%IO_potential)=='SAVE') trove%IO_hamiltonian = 'SAVE'
             !
             if (trim(trove%IO_kinetic)=='SAVE'.and.trove%separate_store) then
                trove%IO_hamiltonian = 'SAVE'
             endif 
             !
           case('BASIS_SET')
             !
             call readu(trove%IO_basisset)
             !
           case('CONTR_CI','CONTR-CI','CONTRCI','CI')
             !
             call readu(trove%IO_contrCI)
             !
             call check_read_save_none(trove%IO_basisset,w)
             !
             !if (trim(trove%IO_contrCI) == "NONE") trove%IO_contrCI = "SAVE"
             !
           case('EIGENFUNC','VECTORS')
             !
             call readu(w)
             !
             job%IOeigen_action = trim(w)
             !
             if (trim(w)/='READ'.and.trim(w)/='SAVE'.and.trim(w)/='APPEND'.and.trim(w)/='NONE') then 
               !
               write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT'
               !
             endif
             !
             ! additional argument 
             !
             if (nitems>=3) then
               !
               call readu(w)
               !
               select case (trim(w))
               !
               case('COMPRESS')
                 !
                 select case (trim(job%IOeigen_action))
                 case ('APPEND', 'SAVE')
                   continue
                 case default
                  call report("Unexpected 3d argument in eigenfuc",.true.)
                 end select
                 !
                 job%IOeigen_compress = .true.
                 !
                 if (nitems>=4) call readf(job%compress)
                 !
               case ('CONVERT')
                 !
                 if (all(trim(job%IOeigen_action)/=(/'SAVE','READ'/))) then
                  call report("Unexpected 3d argument in eigenfuc",.true.)
                 endif
                 !
                 job%convert_model_j0 = .true.
                 !
                 if (all(trim(job%IOeigen_action)==(/'SAVE'/))) then
                   job%IOj0contr_action = 'SAVE'
                 endif
                 !
                 if (job%vib_contract) then
                   !
                   write(out,"('EIGENFUNC SAVE CONVERT cannot be used with MODEL J=0')")
                   call report("illegal usage of  eigenfuc",.true.)
                   !
                 endif
                 !
               case default
                 !
                 write(out,"('illegal keyword after EIGENFUNC XXXX ')")
                 call report("illegal keyword in eigenfuc xxxx",.true.)
                 !
               end select
               ! 
             endif
             !
           case('CONTRACT','CONTR','CONTRACTED')
             !
             call readu(w)
             !
             if (all(trim(w)/=(/'READ','SAVE','NONE'/))) then 
               !
               write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT'
               !
             endif 
             !
             job%IOcontr_action = trim(w)
             !
             if (trim(w)=='SAVE'.and.job%vib_contract) then
               job%IOj0contr_action = 'SAVE'
             endif
             !
             if (nitems>=3) then
               !
               call readu(w)
               !
               if (trim(w)=='CONVERT') then
                  trove%separate_convert = .true.
               else
                   call report (" Illegal record after CONTRACT",.true.)
               endif
               !
             endif 
             !
           case('MATELEM')
             !
             call readu(w)
             !
             select case (trim(w))
             case ('READ','SAVE','NONE','VIB','CONVERT','APPEND')
               continue
             case default
               write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT'
             end select
             !
             if (job%contrci_me_fast.and.trim(w)=='NONE') then 
               !
               !w = 'SAVE'
               !
             endif
             !
             job%IOkinet_action = trim(w)
             !
             if (trim(w)=='CONVERT') then 
                !
                job%vib_contract = .true.
                action%convert_vibme = .true.
                !
                !job%IOj0matel_action = trim('SAVE')
                !
             endif 
             !
             if (trim(w)=='APPEND') then 
                !
                job%matelem_append = .true.
                !
                if (Nitems>3) then
                   call readi(job%iappend)
                else
                   call report (" The reccord to append after is missing",.true.)
                endif
                !
                job%IOkinet_action = trim('SAVE')
                !
             endif 
             !
             if (trim(w)=='VIB') then
               !
               call readu(chk_type)
               !
               job%IOkinet_action = trim(w)//'_'//trim(chk_type) 
               !
             endif 
             !
             if (job%vib_contract) then 
               !
               job%IOj0matel_action = trim(w)
               !
               if (trim(w)=='CONVERT') job%IOj0matel_action = trim('SAVE')
               !
               !if ( any( trim(w)==(/'SAVE','APPEND'/) ) ) action%convert_vibme = .true.
               !
             endif 
             !
             ! an addional key to specify SAVE/READ 
             !
             if (Nitems>2.or.(job%matelem_append.and.Nitems>3)) then
               call readu(w)
               !
               select case (trim(w))
               case ('DIVIDE','SPLIT','STITCH','COLLECT','NON-SPLIT')
                 continue
               case default
                 call report ("Unrecognized unit name (CAN BE SPLIT OR STITCH) "//trim(w),.true.)
               end select
               !
               if (trim(w)=='NON-SPLIT') cycle
               !
               job%iswap(1) = 0
               job%iswap(2) = 12
               !
               if (trim(w)=='DIVIDE') then 
                 job%IOmatelem_divide = .true.
                 job%iswap(1) = 1
                 job%iswap(2) = (trove%Nmodes+3)*3+trove%Nmodes**2+1
               endif
               !
               job%IOmatelem_split  = .true.
               !
               if (job%contrci_me_fast) then 
                 job%iswap(1) = 0
                 job%iswap(2) = 1e6
               endif
               !
               if (Nitems>3) then
                  call readi(job%iswap(1))
                  call readi(job%iswap(2))
               endif
               !
               if (trim(w)=='STITCH'.or.trim(w)=='COLLECT') job%iswap(:)=0
               !
               if (trim(w)=='DIVIDE') job%iswap(:)=0
               !
             endif
             !
             if (job%contrci_me_fast.and..not.job%IOmatelem_split) then
               write(out,"('Read-input: for fast-ci use MATELEM XXXX split (XXXX=read,save,none)')")
               call report ("For fast-ci MATELEM must be split",.true.)
             endif 
             !
             if (item<Nitems) then 
               call readu(w)
               if (trim(w)=='DUMP') job%IOmatelem_dump = .true.
               !
             endif
             !
             if (item<Nitems) then
                call readi(job%iappend)
             endif
             !
           case('DVR')
             !
             call readu(w)
             !
             if (all(trim(w)/=(/'READ','SAVE','NONE'/))) then
               !
               write (out,"('FLinput: illegal key in CHECK_POINT_DVR :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT_DVR'
               !
             endif 
             !
             job%IOdvr_prim = trim(w)
             !
           case('J0_MATELEM')
             !
             call readu(w)
             !
             select case (trim(w))
             case ('READ','SAVE','APPEND','NONE')
               continue
             case default
               write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT'
             end select
             !
             job%IOj0matel_action = trim(w)
             !
             select case (trim(w))
             case ('SAVE', 'APPEND')
               action%convert_vibme = .true.
             end select
             !
             !if (trim(job%IOj0matel_action)=='SAVE') action%convert_vibme = .true.
             !
             ! an addional key to specify SAVE/READ 
             !
             if (Nitems>2) then
               call readu(w)
               !
               select case (trim(w))
               case ('DIVIDE', 'SPLIT')
                 continue
               case default
                 call report ("Unrecognized unit name (<>DIVIDE) "//trim(w),.true.)
               end select
               !
               job%IOmatelem_split = .true.
               !
               job%iswap(1) = 0
               job%iswap(2) = 1e6
               !
               if (Nitems>3) then
                  call readi(job%iswap(1))
                  call readi(job%iswap(2))
               endif
             endif
             !
           case('VECTORS_SYMM','VEC_SYMM','VECTOR_SYMM')
             !
             job%IOvector_symm = .true.
             !
             !job%IOvector_symm = trim(w)
             !
             !if (trim(w)/='READ'.and.trim(w)/='SAVE'.and.trim(w)/='NONE') then 
             !  !
             !  write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
             !  stop 'FLinput -illegal key in CHECK_POINT'
             !  !
             !endif 
             !
           case('J0_EXTERNAL','J0_EXTMATELEM')
             !
             call readu(w)
             !
             select case (trim(w))
             case ('READ','SAVE','APPEND','NONE','DIVIDE','SPLIT')
               continue
             case default
               write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT'
             end select
             !
             job%IOj0ext_action = trim(w)
             !
             select case (trim(job%IOj0ext_action))
             case ('SAVE','APPEND','DIVIDE')
                action%convert_vibme = .true.
                FLextF_matelem = .true.
             end select
             !
             ! an addional key to specify SAVE/READ 
             !
             if (Nitems>2) then
               call readu(w)
               !
               select case (trim(w))
               case ('DIVIDE', 'SPLIT')
                 continue
               case default
                 call report ("Unrecognized unit name (<>DIVIDE) "//trim(w),.true.)
               end select
               !
               job%IOextF_divide = .true.
               !
               if (Nitems>3) then
                  call readi(fitting%iparam(1))
                  call readi(fitting%iparam(2))
               endif
               !
             endif
             !
           case('J0_CONTR','J0_CONTRACT','J0_CONTRACTED')
             !
             call readu(w)
             !
             select case (trim(w))
             case ('READ','SAVE','APPEND','NONE')
               continue
             case default
               write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT'
             end select
             !
             job%IOj0contr_action = trim(w)
             !
             if (any(trim(job%IOj0contr_action)==(/'SAVE'/))) action%convert_vibme = .true.
             !
           case('EXTERNAL')
             !
             call readu(w)
             !
             if (trim(w)/='READ'.and.trim(w)/='SAVE'.and.trim(w)/='NONE') then 
               !
               write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT-EXTERNAL'
               !
             endif 
             trove%IO_ext_coeff = trim(w)
             !
             if (trim(trove%IO_hamiltonian)=='READ'.and.trim(w)=='SAVE') then 
               trove%IO_hamiltonian = 'NONE'
               trove%IO_potential   = 'READ'
               trove%IO_kinetic     = 'READ'
               trove%IO_basisset    = 'SAVE'
             endif 
             !
             if (trim(trove%IO_hamiltonian)=='SAVE'.and.trim(w)/='NONE') then 
               trove%IO_ext_coeff   = 'SAVE'
             endif
             !
             if (trim(trove%IO_ext_coeff)=='SAVE') trove%IO_hamiltonian = 'SAVE'
             !
             if (any(trim(trove%IO_ext_coeff)==(/'SAVE','READ'/))) FLextF_coeffs = .true.
             !
             if (trove%separate_convert.and.trove%separate_store ) then
                  trove%IO_ext_coeff   = 'READ'
             endif 
             !
           case('EXTMATELEM')
             !
             call readu(w)
             !
             select case (trim(w))
             case ('READ','SAVE','MATELEM','NONE','CONVERT','APPEND','STITCH')
               continue
             case default
               write (out,"('FLinput: illegal key in EXTMATELEM :',a)") trim(w)
               stop 'FLinput -illegal key in EXTMATELEM'
             end select
             !
             job%IOextF_action = trim(w)
             !
             if (trim(w)=='CONVERT') then 
                !
                job%vib_contract = .true.
                action%convert_vibme = .true.
                !
                !job%IOj0matel_action = trim('SAVE')
                !
             endif 
             !
             if (job%vib_contract) then 
               !
               job%IOj0ext_action = trim(w)
               !
               if (trim(w)=='CONVERT') job%IOj0ext_action = 'SAVE'
               !
               if (job%IOj0ext_action=='SAVE') FLextF_matelem = .true.
               !
             endif 
             !
             if (trim(w)=='APPEND') then 
                !
                job%extmatelem_append = .true.
                !
                if (Nitems>3) then
                   call readi(job%iextappend)
                else
                   call report (" The reccord to append after is missing",.true.)
                endif
                !
                job%IOextF_action = trim('SAVE')
                !
             endif 
             !
             if (trim(job%IOextF_action)=='SAVE') FLextF_matelem = .true.
             !
             ! an addional key to specify SAVE/READ 
             !
             if (Nitems>2.or.(job%extmatelem_append.and.Nitems>3)) then
               !
               call readu(w)
               !
               select case (trim(w))
                 !
               case ('DIVIDE', 'SPLIT')
                 !
                 job%IOextF_divide = .true.
                 !
               case ('STITCH')
                 !
                 job%IOextF_divide = .true.
                 job%IOextF_stitch = .true.
                 !
               case ( 'DUMP' )
                 !
                 if (trim(w)=='DUMP') job%IOextmatelem_dump = .true.
                 !
               case default
                 !
                 call report ("Unrecognized unit name (<>DIVIDE) "//trim(w),.true.)
                 !
               end select
               !
               if (job%IOextF_divide.and.Nitems>3) then
                  call readi(fitting%iparam(1))
                  call readi(fitting%iparam(2))
               endif
               !
             endif
             !
             if (item<Nitems) then 
               call readu(w)
               if (trim(w)=='DUMP') job%IOextmatelem_dump = .true.
               !
             endif
             !
           case('FORMAT')
             !
             call readu(job%kinetmat_format)
             !
           case('FIT_POTEN','FIT_POT','FITPOTEN')
             !
             call readu(w)
             !
             select case (trim(w))
             case ('READ','SAVE','NONE','VIB','DIVIDE','REWRITE','JOIN','SPLIT')
               continue
             case default
               !
               write (out,"('FLinput: illegal key in CHECK_POINT :',a)") trim(w)
               stop 'FLinput -illegal key in CHECK_POINT'
               !
             end select
             !
             if (trim(w)=='DIVIDE') w = 'SPLIT'
             !
             job%IOfitpot_action = trim(w)
             !
             if (trim(job%IOfitpot_action)/='NONE') then 
                job%IOextF_action = 'READ'
                !
                if (job%vib_contract) then 
                  job%IOj0ext_action = 'READ'
                endif
                !
             endif
             !
             ! an addional key to specify SAVE/READ 
             !
             if (Nitems>2) then
               call readu(w)
               select case (trim(w))
               case ('DIVIDE', 'SPLIT')
                 continue
               case default
                 call report ("Unrecognized unit name (<>DIVIDE) "//trim(w),.true.)
               end select
               !
               job%IOextF_divide = .true.
               job%IOfitpot_divide = .true.
               !
               if (Nitems>3) then
                  call readi(fitting%iparam(1))
                  call readi(fitting%iparam(2))
               endif
               !
             endif
             !
           case('CHK_FILE')
             !
             call readu(w)
             call locase(w)
             !
             trove%chk_fname = trim(w)
             !
           case('NUMEROV_FILE')
             !
             call readu(w)
             call locase(w)
             !
             trove%chk_numerov_fname = trim(w)
             !
           case('HAMILTONIAN_FILE')
             !
             call readu(w)
             call locase(w)
             !
             trove%chk_hamil_fname = trim(w)
             !
           case('SWAP','SWAPMATELEM','MATELEMSWAP')
             !
             call readu(w)
             !
             select case (trim(w))
             case ('NONE','DIVIDE','JOIN','SAVE','READ','SPLIT')
               continue
             case default
               !
               write (out,"('FLinput: illegal key in CONTRSWAP :',a)") trim(w)
               stop 'FLinput -illegal key in CONTRSWAP'
               !
             end select
             !
             job%IOswap_matelem = trim(w)
             !
           case('ISWAP')
             !
             call readi(job%iswap(1))
             call readi(job%iswap(2))
             !
           case('CONTR_FILE')
             !
             call readu(w)
             call locase(w)
             !
             job%contrfile%dscr       = trim(w)//'_descr.chk'
             job%contrfile%primitives = trim(w)//'_quanta.chk'
             job%contrfile%vectors    = trim(w)//'_vectors.chk'
             job%contrfile%dvr        = trim(w)//'_dvr.chk'
             !
           case('EIGEN_FILE')
             !
             call readu(w)
             call locase(w)
             !
             job%eigenfile%filebase   = trim(w)
             !
             job%eigenfile%dscr       = trim(w)//'_descr'
             job%eigenfile%primitives = trim(w)//'_quanta'
             job%eigenfile%vectors    = trim(w)//'_vectors'
             !
           case default
             !
             call report ("Unrecognized unit name "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (trim(w)/="".and.trim(w)/="END") then 
            !
            write (out,"('FLinput: wrong last line in CHECK_POINTS =',a)") trim(w)
            stop 'FLinput - illegal last line in CHECK_POINTS'
            !
         endif 
         !
       case("ANALYSIS")
         !
         call readu(w)
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
             !
           case('NONE')
             !
             analysis%density = .false.
             analysis%classical = .false.
             !
           case('DENSITY')
             !
             analysis%reduced_density = .true.
             analysis%density = .true.
             !
           case("PRINT_CONTRIBUTION")
             analysis%reducible_eigen_contribution = .true.   
             i = 0
             do while(trim(w)/="".and.item<Nitems.and.i<100)
                i = i + 1
                call readi(analysis%j_list(i)) 
             enddo 
             analysis%density = .true.
             !
           case('PRINT_VECTOR')
             !
             analysis%print_vector = .true.
             !
             if (nitems>1) call readf(analysis%threshold)
             !
           case('ROT_MATRIX')
             !
             analysis%rotation_matrix = .true.
             analysis%density = .true.
             !
           case('EXTERNAL')
             !
             analysis%extF = .true.
             analysis%density = .true.
             !
           case('TEST_HAMILTONIAN')
             !
             analysis%rotation_matrix = .true.
             analysis%density = .true.
             !
           case('RES_TYPE')
             !
             call readu(analysis%res%type)
             analysis%classical = .true.
             !
           case('RES')
             !
             analysis%rotation_energy_surface = .true.
             analysis%classical = .true.
             !
             do while (trim(w)/="".and.item<Nitems)
                !
                call readu(w)
                !
                select case(w)
                    !
                  case('NTHETA')
                    !
                    call readi(analysis%res%ntheta)
                    !
                  case('NPHI')
                    !
                    call readi(analysis%res%nphi)
                    !
                  case('NTAU')
                    !
                    call readi(analysis%res%ntau)
                    !
                  case('THETA')
                    !
                    call readf(analysis%res%theta)
                    !
                  case('PHI')
                    !
                    call readf(analysis%res%phi)
                    !
                  case('TAU')
                    !
                    call readf(analysis%res%tau)
                    !
                  case('THETA_RESTART')
                    !
                    call readf(analysis%res%theta1)
                    !
                end select
                !
             enddo
             !
           case('ROT-DENSITY','ROT_DENSITY')
             !
             analysis%rotation_density = .true.
             !
             do while (trim(w)/="".and.item<Nitems)
                !
                call readu(w)
                !
                select case(w)
                    !
                  case('THETA')
                    !
                    call readf(analysis%res%theta) ; analysis%res%theta = analysis%res%theta/rad
                    !
                  case('PHI')
                    !
                    call readf(analysis%res%phi)  ; analysis%res%phi = analysis%res%phi/rad
                    !
                  case('NTHETA')
                    !
                    call readi(analysis%res%ntheta)
                    !
                  case('NPHI')
                    !
                    call readi(analysis%res%nphi)
                    !
                end select
                !
             enddo
             !
           case('DENSITY_LIST','LIST')
             !
             do while (trim(w)/="".and.item<Nitems)
                !
                call readu(w)
                !
                i = 0
                do while (trim(w)/="MODES".and.item<Nitems.and.i<100)
                  !
                  i = i + 1 
                  !
                  read(w,*) analysis%j_list(i)
                  call readi(analysis%sym_list(i))
                  call readi(analysis%dens_list(i))
                  !
                  call readu(w)
                  !
                enddo
                !
                if (trim(w)/='MODES') then 
                  call report ("Illegal use of DENSITY_LIST"//trim(w),.true.)
                endif
                !
                i_t = i
                !
                do i = i_t+1,i_t+3
                  !
                  if (item<Nitems) then 
                    !
                    call readi(analysis%dens_list(i))
                    !
                  else
                    !
                    analysis%dens_list(i) = analysis%dens_list(i-1)
                    !
                  endif 
                  !
                enddo
                !
             enddo
             !
           case default
             !
             call report ("Unrecognized unit name "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (trim(w)/="".and.trim(w)/="END") then 
            !
            write (out,"('FLinput: wrong last line in CHECK_POINTS =',a)") trim(w)
            stop 'FLinput - illegal last line in CHECK_POINTS'
            !
         endif 
         !
       case("KIN","KINET","KINETIC")
         !
         if (kinetic_defined) then 
            write (out,"(' Error: trying to read KINETIC  second time!')") 
            stop 'FLinput - reading KINETIC  second time'
         endif 
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
              !
           case("KIN_TYPE","KINETIC_TYPE")
              !
              call readu(w)
              !
              trove%kinetic_type = trim(w)
              !
              case default
               !
               call report ("Unrecognized unit name "//trim(w),.true.)
               !
           end select 
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         kinetic_defined = .true.
         !
       case("POT","POTEN","POTENTIAL","PES")
         !
         if (pot_defined) then 
            write (out,"(' Error: trying to read POTEN  second time!')") 
            stop 'FLinput - reading POTEN  second time'
         endif 
         !
         if (Nmodes==0) then 
            !
            write (out,"('FLinput: POTEN cannot appear before NMODES')") 
            stop 'FLinput - POTEN defined befor NMODES'
            !
         endif 
         !
         pot_coeff_type = 'LIST'
         !
         ! read Nparam and Type of PES
         !
         do i = 1,3
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
            select case(w)
            !
            case("POT_TYPE")
              !
              call readu(w)
              !
              trove%potentype = trim(w)
              !
            case("NPARAM")
              !
              call readi(Nparam)
              !
            case("COEFF")
              !
              call readu(w)
              !
              pot_coeff_type = trim(w)
              !
            case default
              !
              call report ("Unrecognized unit name "//trim(w),.true.)
              !
            end select
            !
         enddo
         !
         ! Allocation of the pot. parameters 
         !
         allocate (force(Nparam),forcename(Nparam),pot_ind(1:trove%Ncoords,Nparam),ifit(Nparam),stat=alloc)
         if (alloc/=0) then
            write (out,"(' Error ',i9,' allocating matix force,forcename ')") alloc
            stop 'FLinput - cannot allocate force,forcename'
         end if
         !
         iparam = 0 
         pot_ind = 0
         ifit = 0  
         force = 0 
         !
         forcename = 'dummy'
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.iparam<Nparam.and.trim(w)/="END")
            !
            iparam = iparam+1 
            !
            select case(trim(pot_coeff_type))
            !
            case("LIST")
              !
              forcename(iparam)=trim(w)
              call readi(ifit(iparam))
              call readf(force(iparam))
              !
            case("POWERS")
              !
              if (nitems<trove%Ncoords+3) then 
                 !
                 write (out,"('FLinput: wrong number of records in POTEN on row=',i6,'()')") iparam
                 stop 'FLinput - illigal number of records in POTEN'
                 !
              endif
              !
              forcename(iparam)=trim(w)
              !
              do i=1,trove%Ncoords
                 call readi(pot_ind(i,iparam))
              enddo
              !
              call readi(ifit(iparam))
              call readf(force(iparam))
              !
              ! trick to prevent writing error
              !
              !if (all(pot_ind(1:Ncoords,iparam)<10)) then
              !  write(forcename(iparam),"('f',<Ncoords>i1)") pot_ind(1:Ncoords,iparam)
              !else
              !  write(forcename(iparam),"('f',<Ncoords>i1)") 0,pot_ind(2:Ncoords,iparam)
              !endif
              !
              if (any(pot_ind(:,iparam)<0)) then 
                  write(out,"('FLinput: negative POTEN powers on row',i8)") iparam
                    stop 'FLinput: wrong indexes '
              endif 
              !
            case default
              !
              call report ("Unrecognized unit name "//trim(w),.true.)
              !
            end select
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
         enddo 
         !
         if ((iparam/=Nparam.and.trim(forcename(1))=='dummy').or.(trim(w)/="".and.trim(w)/="END")) then 
            !
            write (out,"('FLinput: wrong number of rows in POTEN ',i6,'()')") iparam,Nparam
            stop 'FLinput - illigal number of rows in POTEN'
            !
         endif 
         !
         pot_defined = .true.
         !
       case("MEP","RPH")
         !
         ! read Nparam and Type of PES
         !
         pot_coeff_type = 'LIST'
         !
         do i = 1,3
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
            select case(w)
            !
            case("MEP_TYPE")
              !
              call readu(w)
              !
              molec%meptype = trim(w)
              !
            case("NPARAM")
              !
              call readi(Nparam)
              !
              molec%N_meppars = Nparam 
              !
            case("COEFF")
              !
              call readu(w)
              !
              pot_coeff_type = trim(w)
              !
            case default
              !
              call report ("Unrecognized unit name "//trim(w),.true.)
              !
            end select
            !
         enddo
         !
         ! Allocation of the pot. parameters 
         !
         allocate (molec%mep_params(Nparam),stat=alloc)
         if (alloc/=0) then
            write (out,"(' Error ',i9,' allocating matix mep_params ')") alloc
            stop 'FLinput - cannot allocate mep_params'
         end if
         !
         if (trim(pot_coeff_type)=="POWERS") then 
           allocate (molec%mep_ind(1:trove%Ncoords,Nparam))
         else
           allocate (molec%mep_ind(1,Nparam))
         endif
         !
         iparam = 0 
         molec%mep_ind = 0
         molec%mep_params = 0 
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.iparam<Nparam.and.trim(w)/="END")
            !
            iparam = iparam+1 
            !
            !call readf(molec%mep_params(iparam))
            !
            select case(trim(pot_coeff_type))
            !
            case("LIST")
              !
              call readi(molec%mep_ind(1,iparam))
              call readf(molec%mep_params(iparam))
              !
            case("POWERS")
              !
              if (nitems<trove%Ncoords+2) then 
                 !
                 write (out,"('FLinput: wrong number of records in POTEN on row=',i6,'()')") iparam
                 stop 'FLinput - illigal number of records in POTEN'
                 !
              endif 
              !
              do i=1,trove%Ncoords
                 !
                 call readi(molec%mep_ind(i,iparam))
                 !
              enddo
              !
              call readf(molec%mep_params(iparam))
              !
              if (any(molec%mep_ind(:,iparam)<0)) then 
                  write(out,"('FLinput: negative MEP powers on row',i8)") iparam
                    stop 'FLinput: wrong MEP-indexes '
              endif 
              !
            case default
              !
              call report ("Unrecognized unit name "//trim(w),.true.)
              !
            end select
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
         enddo 
         !
         if ((iparam/=Nparam.and.trim(w)/="".and.trim(w)/="END")) then 
            !
            write (out,"('FLinput: wrong number of rows in MEP ',i6,'()')") iparam,Nparam
            stop 'FLinput - illigal number of rows in MEP'
            !
         endif 
         !
       case("SPECPARAM","SPECPARAMETERS","SPEC","PARAMETERS")
         !
         if (trove%Ncoords==0) then 
            !
            write (out,"('FLinput: SPECPARAM cannot appear before Zmat')") 
            stop 'FLinput - SPECPARAM defined befor Zmat'
            !
         endif 
         !
         iparam = 0 
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         trove%specparam = 0
         !         !
         do while (trim(w)/="".and.iparam<trove%Ncoords.and.trim(w)/="END")
            !
            iparam = iparam+1
            !
            call readi( i_t  )  ! reading an integer weight, not used so far 
            !
            call readf( trove%specparam(iparam)  )
            !
            if (nitems==4) then
              !
              call readu(w)
              !
              lfact = 1.0_rk
              !
              select case(w)
              !
              case("ANGSTROM")
                !
                lfact= 1.0_rk
                !
              case("BOHR")
                !
                lfact=bohr
                !
              case("DEG","DEGREE","DEGREES")
                !
                lfact=1.0_rk/rad
                !
              case default
                !
                call report ("Unrecognized unit name "//trim(w),.true.)
                !
              end select
              !
              trove%specparam(iparam) = trove%specparam(iparam)*lfact
              !
            endif
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
         enddo 
         !
         if (iparam>trove%Ncoords.or.(trim(w)/="".and.trim(w)/="END")) then 
            !
            write (out,"('FLinput: wrong number of rows in SPECPARAM ',2i6,'()')") iparam,trove%Ncoords
            stop 'FLinput - illigal number of rows in SPECPARAM'
            !
         endif 
         !
         ! Intensity section 
         !
       case("INTENSITY")
         !
         if (Nmodes==0) then 
            !
            write (out,"('FLinput: INTENSITY cannot appear before NMODES')") 
            stop 'FLinput - INTENSITY defined befor NMODES'
            !
         endif 
         !
         if (.not.symmetry_defined) then 
            !
            write (out,"('FLinput: INTENSITY cannot appear before symmetry is defined')") 
            stop 'FLinput - INTENSITY defined before symmetry'
            !
         endif 
         !
         allocate(intensity%gns(sym%Nrepresen),intensity%isym_pairs(sym%Nrepresen),intensity%v_low(trove%nmodes,2),&
                  intensity%v_upp(trove%nmodes,2),stat=alloc)
         if (alloc/=0) then
             write (out,"(' Error ',i9,' trying to allocate qwrange')") alloc
             stop 'FLinput, qwrange - out of memory'
         end if
         !
         ! defauls values
         !
         intensity%gns = 1
         forall(i=1:sym%Nrepresen) intensity%isym_pairs(i) = 1
         intensity%v_low(:,1) = 0 ; intensity%v_low(:,2) = job%bset(1:)%range(2)
         intensity%v_upp(:,1) = 0 ; intensity%v_upp(:,2) = job%bset(1:)%range(2)
         !
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
           !
           case('NONE','ABSORPTION','EMISSION','TM','DIPOLE-TM','RAMAN','POLARIZABILITY-TM','PARTFUNC')
             !
             intensity%action = trim(w)
             !
             if (trim(intensity%action)=='DIPOLE-TM') intensity%action = 'TM'
             !
             if (trim(intensity%action)=='TM'      .or.trim(intensity%action)=='ABSORPTION'.or.&
                 trim(intensity%action)=='EMISSION'.or.trim(intensity%action)=='PARTFUNC') then 
               !
               action%intensity = .true.
               intensity%do = .true.
                 !
               job%IOextF_action = 'READ'
               !
               if (job%vib_contract) then 
                 !
                 job%IOj0ext_action = 'READ'
                 !
               endif 
               !
             endif
             !
           case('THRESH_INTES','THRESH_TM','THRESH_INTENS')
             !
             call readf(intensity%threshold%intensity)
             !
           case('THRESH_LINE','THRESH_LINESTRENGHT','THRESH_LINESTRENGTH')
             !
             call readf(intensity%threshold%linestrength)
             !
           case('THRESH_COEFF','THRESH_COEFFICIENTS')
             !
             call readf(intensity%threshold%coeff)
             !
           case('TEMPERATURE')
             !
             call readf(intensity%temperature)
             !
           case('EXOMOL')
             !
             job%exomol_format = .true.
             !
           case('QSTAT','PARTITION','PART_FUNC')
             !
             call readf(intensity%part_func)
             !
           case('PRUNING')
             !              
             intensity%pruning = .true.
             !
           case('TDM_REPLACE','DIPOLE_REPLACE','DIPOLE_SCALE')
             !              
             intensity%tdm_replace = .true.
             !
           case('OUTPUT')
             !
             call readu(w)
             !              
             if (trim(w)=='SHORT') then 
               intensity%output_short = .true.
             elseif(trim(w)=='LONG') then
               intensity%output_short = .false.
             else
               call report ("Illegal OUTPUT value (expected SHORT or LONG) "//trim(w),.true.)
             endif
             !
           case('GNS')
             !
             i = 0
             !
             intensity%gns = 0
             if( trim(intensity%action)=='TM') intensity%gns = 1 
             job%select_gamma = .false.
             job%isym_do = .false.
             !
             do while (item<Nitems.and.i<sym%Nrepresen)
               !
               i = i + 1
               call readf(intensity%gns(i))
               if (intensity%gns(i)>small_) then 
                 !
                 job%isym_do(i) = .true.
                 job%select_gamma(i) = .true.
                 !
               endif
               !
             enddo
             !
             !if (i/=sym%Nrepresen) then 
             !  !
             !  write (out,"('FLinput: illegal number entries in gns: ',i8,' /= ',i8)") i,sym%Nrepresen
             !  stop 'FLinput - illegal number entries in gns'
             !  !
             !endif 
             !
           case('SELECTION','SELECTION_RULES','SELECT','PAIRS')
             !
             ! default values are by pairs: 1 1 2 2 3 3 4 4 ...
             do i=1,sym%Nrepresen,2
               intensity%isym_pairs(i  ) = (i+1)/2
               if (i+1<=sym%Nrepresen) intensity%isym_pairs(i+1) = (i+1)/2
             enddo
             !
             if( trim(intensity%action)=='TM') intensity%isym_pairs = 1 
             !
             i = 0
             !
             do while (item<Nitems.and.i<sym%Nrepresen)
               !
               i = i + 1
               !
               call readi(intensity%isym_pairs(i))
               !
             enddo
             !
             if (sym%Nrepresen<=8.and.i/=sym%Nrepresen) then 
               write (out,"('FLinput: illegal number entries in SELECTION for Nrepresen<=8',i8,' /= ',i8)") i,sym%Nrepresen
               stop 'FLinput - illegal number entries in SELECTION, Nentries<>Nrepresen'
             endif 
             !
           case('ZPE')
             !
             call readf(intensity%zpe)
             job%zpe = intensity%zpe
             !
           case('SWAP_SIZE')
             !
             call readi(intensity%swap_size)
             !
           case('SWAP')
             !
             call readu(intensity%swap)
             !
           case('SYMMETRY')
             !
             call readu(w)
             !
             if (trim(w)=='REDUCED') then
               intensity%reduced = .true.
             else
               call report ("Unrecognized SYMMETRY value "//trim(w),.true.)
             endif
             !
           case('FACTOR')
             !
             call readf(intensity%factor)
             !
           case('INCREMENT')
             !
             call readi(intensity%int_increm)
             !
           case('WALLCLOCK','WALL')
             !
             call readf(intensity%wallclock)
             !
           case('J')
             !
             call readi(intensity%j(1))
             call readi(intensity%j(2))
             !
           case('FREQ-WINDOW','FREQ','NU','FREQUENCY')
             !
             call readf(intensity%freq_window(1))
             call readf(intensity%freq_window(2))
             !
           case('ENERGY')
             !
             call readu(w)
             !
             do while (trim(w)/="")
                !
                select case(w)
                !
                case("LOWER","LOW","L")
                  !
                  call readf(intensity%erange_low(1))
                  call readf(intensity%erange_low(2))
                  !
                case("UPPER","UPP","UP","U")
                  !
                  call readf(intensity%erange_upp(1))
                  call readf(intensity%erange_upp(2))
                  !
                end select 
                !
                call readu(w)
                !
             enddo 
             !
           case('QUANTA','V')
             !
             call readi(imode)
             !
             call readu(w)
             !
             do while (trim(w)/="")
                !
                select case(w)
                !
                case("LOWER","LOW","L")
                  !
                  call readi(intensity%v_low(imode,1))
                  call readi(intensity%v_low(imode,2))
                  !
                case("UPPER","UPP","UP","U")
                  !
                  call readi(intensity%v_upp(imode,1))
                  call readi(intensity%v_upp(imode,2))
                  !
                end select 
                !
                call readu(w)
                !
             enddo 
             !
             case default
               !
               call report ("Unrecognized unit name "//trim(w),.true.)
             !
           end select 
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         if (trim(intensity%action) == 'ABSORPTION'.or.trim(intensity%action) == 'EMISSION') then 
           !
           ! define the selection pairs by gns if not yet defined
           !
           i_t = 0
           !
           do i = 1,sym%Nrepresen
             !
             if (intensity%isym_pairs(i)/=0) cycle
             !
             do j = 1,sym%Nrepresen
               !
               if (i/=j.and.intensity%isym_pairs(j)==0.and.intensity%gns(i)==intensity%gns(j)) then 
                 !
                 i_t = i_t + 1
                 !
                 intensity%isym_pairs(i) = i_t
                 intensity%isym_pairs(j) = i_t
                 !
               endif 
               !
             enddo
             !
           enddo
           !
         endif 
         !
         job%erange(1) = min(intensity%erange_low(1),intensity%erange_upp(1))
         job%erange(2) = max(intensity%erange_low(2),intensity%erange_upp(2))
         !
       case("FITTING")
         !
         call readu(w)
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         action%fitting = .true.
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
             !
           case('NONE')
             !
             action%fitting = .false.
             !
           case('J_LIST','JROTLIST','J-LIST','JLIST')
             !
             i = 0
             do while (item<Nitems.and.i<100)
                !
                i = i + 1 
                !
                call readi(fitting%j_list(i))
                !
             enddo
             !
           case('ITMAX','ITERMAX','ITER')
             !
             call readi(fitting%itermax)
             !
           case('ROBUST')
             !
             call readf(fitting%robust)
             !
           case('WATSON')
             !
             call readf(fitting%watson)
             !
           case('TARGET_RMS')
             !
             call readf(fitting%target_rms)
             !
           case('METHOD')
             !
             call readu(fitting%method)
             !
           case('FIT_TYPE')
             !
             call readu(fitting%fit_type)
             !
           case('THRESH_ASSIGN','THRESH_REASSIGN','LOCK')
             !
             call readf(fitting%threshold_lock)
             !
           case('THRESH_OBS-CALC') 
             ! switch off weights for residuals larger than THRESH_OBS-CALC
             !
             call readf(fitting%threshold_obs_calc)
             !
           case('FIT_SCALE') 
             ! parameter to scale the correction dx for fitting x = x + scale*dx, default = 0.4
             !
             call readf(fitting%fit_scale)
             !
           case("IPARAM")
             !
             call readi(fitting%iparam(1))
             call readi(fitting%iparam(2))
             !
           case("SYM","SYMM","SYMMETRY","SYMMETRIES","GNS")
             !
             job%isym_do = .false.
             job%select_gamma = .false.
             !
             do while (item<Nitems.and.item-1<sym%Nrepresen)
               !
               call readi(i)
               if (i>0) job%isym_do(i) = .true.
               if (i>0) job%select_gamma(i) = .true.
               !
             enddo
             !
           case('ENERCUT')
             !
             call readf(job%erange(2))
             !
           case('GEOMETRIES')
             !
             call readl(fitting%geom_file)
             !
           case('OUTPUT')
             !
             call readl(fitting%output_file)
             !
           case('FIT_FACTOR')
             !
             call readf(fitting%factor)
             !
           case('THRESH_COEFF','THRESH_COEFFICIENTS')
             !
             call readf(fitting%threshold_coeff)
             !
           case('OBS','OBS_ENERGIES')
             !
             call readi(fitting%Nenergies)
             !
             allocate (fitting%obs(1:fitting%Nenergies),stat=alloc)
             if (alloc/=0) then
               write (out,"(' Error ',i8,' initializing obs. energy related arrays')") alloc
               stop 'obs. energy arrays - alloc'
             end if
             !
             do i = 1,fitting%Nenergies
               allocate(fitting%obs(i)%quanta(0:trove%nmodes),stat=alloc)
               if (alloc/=0) then
                 write (out,"(' Error ',i8,' initializing obs%quanta')") alloc
                 stop 'initializing obs%quanta - alloc'
                end if
             enddo
             !
             i = 0
             !
             call read_line(eof) ; if (eof) exit
             call readu(w) 
             !
             do while (trim(w)/="END".and.i<fitting%Nenergies)
                !
                i = i + 1
                !
                if (nitems<trove%Nmodes+5) then 
                   !
                   write (out,"('FLinput: wrong number of records in obs_fitting_energies on row=',i6,'()')") i
                   stop 'FLinput - illigal number of records in obs_fitting_energies'
                   !
                elseif(nitems==trove%Nmodes+5) then
                   !
                   write (out,"('FLinput: old format in number of records in obs_fitting_energies on row=',i6)") i
                   write (out,"('Please include K-quantum number as the column after energies')")
                   stop 'FLinput - illigal number of records in obs_fitting_energies, missing K?'
                   !
                endif 
                !
                read(w,*) fitting%obs(i)%Jrot
                !
                call readi(fitting%obs(i)%symmetry) 
                call readi(fitting%obs(i)%N)
                call readf(fitting%obs(i)%energy)
                !
                do j=0,trove%nmodes
                  call readi( fitting%obs(i)%quanta(j) )
                enddo
                !
                call readf(fitting%obs(i)%weight)
                !
                call read_line(eof) ; if (eof) exit
                call readu(w)
                !
             enddo
             !
           case('J0FIT','J0FITTING')
             !
             action%band_fitting = .true.
             !
             call readi(j0fit%Nenergies)
             !
             allocate (j0fit%obs(1:j0fit%Nenergies),stat=alloc)
             if (alloc/=0) then
               write (out,"(' Error ',i8,' initializing obs. energy related arrays')") alloc
               stop 'obs. energy arrays - alloc'
             end if
             !
             do i = 1,j0fit%Nenergies
               allocate(j0fit%obs(i)%quanta(0:trove%nmodes),stat=alloc)
               if (alloc/=0) then
                 write (out,"(' Error ',i8,' initializing obs%quanta')") alloc
                 stop 'initializing obs%quanta - alloc'
                end if
             enddo
             !
             i = 0
             !
             call read_line(eof) ; if (eof) exit
             call readu(w) 
             !
             do while (trim(w)/="END".and.i<j0fit%Nenergies)
                !
                i = i + 1
                !
                if (nitems<trove%Nmodes+4) then 
                   !
                   write (out,"('FLinput: wrong number of records in j0fit_energies on row=',i6,'()')") i
                   stop 'FLinput - illigal number of records in j0fit_energies'
                   !
                endif 
                !
                !read(w,"(i)") j0fit%obs(i)%Jrot
                !
                read(w,*) j0fit%obs(i)%symmetry
                !
                call readi(j0fit%obs(i)%N)
                call readf(j0fit%obs(i)%energy)
                !
                do j=1,trove%nmodes
                  call readi( j0fit%obs(i)%quanta(j) )
                enddo
                !
                call readf(j0fit%obs(i)%weight)
                !
                call read_line(eof) ; if (eof) exit
                call readu(w)
                !
             enddo
             !
           case default
             !
             call report ("Unrecognized unit name "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (trim(w)/="".and.trim(w)/="END") then 
            !
            write (out,"('FLinput: wrong last line in CHECK_POINTS =',a)") trim(w)
            stop 'FLinput - illegal last line in CHECK_POINTS'
            !
         endif 
         !
       case("DMS","DIPOLE","EXTERNAL")
         !
         extF_defined = .true.
         !
         if (Nmodes==0.or.Ncoords==0) then 
            !
            write (out,"('FLinput: EXTERNAL cannot appear before NMODES and Ncoords defined')") 
            stop 'FLinput - EXTERNAL defined befor NMODES'
            !
         endif 
         !
         ! read Nparam and Type of PES
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
         do while(w(1:5)/='PARAM')
            !
            select case(w)
              !
            case("RANK")
              !
              call readi(extF%rank)
              !
              if (all(fitting%iparam(:)==(/1,1000000/))) fitting%iparam = (/1,extF%rank/)
              !
            case("DMS_TYPE","TYPE")
              !
              call readu(w)
              !
              extF%ftype = trim(w)
              !
           case("IPARAM")
             !
             if (fitting%iparam(2)==1000000) then 
               write(out,"('it is illegal to define iparam before rank in external field')")
               stop 'iparam is being defined before rank in external field'
             endif
             !
             if (fitting%iparam(1)==1)   call readi(fitting%iparam(1))
             if (fitting%iparam(1)==1e6) call readi(fitting%iparam(2))
              !
            case("NPARAM")
              !
              if (nitems/=extF%rank+1.and.nitems/=2) then
                 write (out,"('wrong number of records in  EXTF-NPARAM, neither  ',i5,', nor  1')") extF%rank
                 stop 'FLinput - wrong number of records in  EXTF-NPARAM'
              end if
              !
              allocate(extF%nterms(1:extF%rank),extF%maxord(1:extF%rank),stat=alloc)
              if (alloc/=0) stop 'FLinput - cannot allocate extF%nterms'
              !
              do i=1,min(extF%rank,nitems-1)
                call readi(extF%nterms(i))
              enddo
              !
              extF%nterms(nitems:extF%rank) = extF%nterms(1)
              !
              Nparam = maxval(extF%nterms(:),dim=1)
              !
              if (Nparam<1) then
                 write (out,"('wrong value of EXTF-NPARAM = ',i5,', must be >0')") Nparam
                 stop 'FLinput - wrong value of EXTF-NPARAM'
              end if
              !
            case("COEFF")
              !
              call readu(w)
              !
              exfF_coeff_type = trim(w)
              !
            case("THRESHOLD")
              !
              call readf(extF%matelem_threshold)
              !
            case("COORDS")
              !
              if (nitems-1==1) then 
                 !
                 call readu(w)
                 extF%intcoords(:) = trim(w)
                 !
                 call read_line(eof) ; if (eof) exit
                 call readu(w)
                 cycle 
                 !
              endif 
              !
              if (nitems/=trove%Nmodes+1) then
                 write (out,"('wrong number of records in  EXTF-COORDS for trove%nmodes = ',i5)") trove%nmodes
                 stop 'FLinput - wrong number of records in  EXTF-COORDS'
              end if
              !
              do imode=1,trove%nmodes
                 call readu(extF%intcoords(imode))
              enddo
              !
            case("REF_GEOM","GEOM_REF")
              !
              if (nitems<Ncoords+1) then
                 write (out,"('wrong number of records in  GEOM_REF for trove%Ncoords = ',i5)") trove%Ncoords
                 stop 'FLinput - wrong number of records in  GEOM_REF'
              end if
              !
              imode = 0
              !
              !do while (imode<Ncoords.and.i<nitems)
              !
              do i=2,nitems
                !
                call readu(w)
                !
                lfact = 1.0_rk
                !
                select case(w)
                !
                case("ANGSTROM")
                  !
                  lfact= 1.0_rk
                  !
                case("BOHR")
                  !
                  lfact=bohr
                  !
                case("DEG","DEGREE","DEGREES")
                  !
                  lfact=1.0_rk/rad
                  !
                case default
                  !
                  imode = imode + 1
                  !
                  if (imode>Ncoords) then
                     write (out,"('too many records in  GEOM_REF for trove%Ncoords = ',i5)") trove%Ncoords
                     stop 'FLinput - wrong number of records in  GEOM_REF'
                  end if
                  !
                  read(w,*) extF%geom_ref(imode)
                  !
                end select
                !
                extF%geom_ref(imode) = extF%geom_ref(imode)*lfact
                !
              enddo
              !
            case ("DSTEP_BMAT")
              !
              call readf(fd_step_Bmat)
              !
            case ("DSTEP")
              !
              ! If all dsteps are the same - only one number in the input can be given. 
              !
              if (nitems-1==1) then 
                 !
                 call readf(f_t)
                 extF%fdstep(:) = f_t
                 call read_line(eof) ; if (eof) exit
                 call readu(w)
                 cycle 
                 !
              endif 
              !
              if (nitems-1/=Nmodes) then 
                 !
                 write (out,"('FLinput: wrong number elements in extF%dstep : ',i8)") nitems-1
                 stop 'FLinput - illigal number of extF%dstep'
                 !
              endif 
              !
              do i =1,Nmodes
                 !
                 call readf(f_t)
                 extF%fdstep(i) = f_t
                 !
              end do
              !
            case("DIPORDER","ORDER","ORDERS")
              !
              if (nitems/=extF%rank+1.and.nitems/=2) then
                 write (out,"('wrong number of records in  EXTF-ORDER for rank = ',i5)") extF%rank
                 stop 'FLinput - wrong number of records in  EXTF-ORDER'
              end if
              !
              do i=1,min(extF%rank,nitems-1)
                call readi(extF%maxord(i))
              enddo
              !
              extF%maxord(nitems:extF%rank) = extF%maxord(1)
              !
              trove%NExtOrder = maxval(extF%maxord(:),dim=1)
              !
            case default
              !
              call report ("Unrecognized unit name "//trim(w),.true.)
              !
            end select
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
         enddo
         !
         ! Allocation of the pot. parameters 
         !
         allocate (extF%coef(Nparam,extF%rank),extF%name(Nparam,extF%rank),& 
                   extF%term(1:trove%Ncoords,Nparam,extF%rank),&
                   extF%ifit(Nparam,extF%rank),stat=alloc)
         if (alloc/=0) then
            write (out,"(' Error ',i9,' allocating matix extF%coef ')") alloc
            stop 'FLinput - cannot allocate extF%coef'
         end if
         !
         iparam = 0
         extF%coef = 0
         extF%name = 'xxxx'
         extF%term = -1
         extF%ifit = 1
         !
         Nparam = sum(extF%nterms(:))
         !
         do imu = 1,extF%rank
           !
           do iterm = 1, extF%nterms(imu)
             !
             iparam = iparam+1 
             !
             call read_line(eof) ; if (eof) exit
             call readu(w)
             !
             if (trim(w)=="".or.trim(w)=="END") then
                 write (out,"('FLinput: wrong number of rows in DMS, imu = ',i6,' i = ',i8)") imu,iterm
                 stop 'FLinput - illigal number of rows in DMS'
             endif
             !
             select case(trim(exfF_coeff_type))
             !
             case("LIST")
               !
               extF%name(iterm,imu)=trim(w)
               !
               call readi(i_t) ; extF%ifit(iterm,imu) = i_t
               call readf(f_t) ; extF%coef(iterm,imu) = f_t
               !
             case("POWERS")
               !
               if (nitems<trove%Ncoords+3) then 
                  !
                  write (out,"('FLinput: wrong number of records in extF on row=',i6,'()')") iparam
                  stop 'FLinput - illigal number of records in extF'
                  !
               endif 
               !
               do i=1,trove%Ncoords
                  !
                  call readi(extF%term(i,iterm,imu))
                  !
               enddo
               !
               call readi(i_t); extF%ifit(iterm,imu) = i_t
               call readf(f_t); extF%coef(iterm,imu) = f_t
               !
               write(my_fmt,'(a,i0,a)') "(a,",Ncoords,"i1)"
               !
               write(extF%name(iterm,imu),my_fmt) 'f',(mod(extF%term(i,iterm,imu),10),i=1,trove%Ncoords)
               !
               if (any(extF%term(:,iterm,imu)<0)) then 
                   write(out,"('FLinput: negative extF powers on row',i8)") iparam
                   stop 'FLinput: wrong extF indexes '
               endif 
               !
             case default
               !
               call report ("Unrecognized unit name "//trim(w),.true.)
               !
             end select
             !
           enddo
           !
         enddo
         !
         call read_line(eof) ; if (eof) exit
         call readu(w)
         !
       case default
         call report ("Principal keyword "//trim(w)//" not recognized",.true.)
       end select
       !
   end do
   !
   select case (trim(job%IOextF_action))
     !
   case ('SAVE', 'DIVIDE', 'SPLIT')
     if (trim(trove%IO_ext_coeff) == 'NONE') then
      FLextF_coeffs = .true.
      !
      !write(out,"('FLinput - EXTF-coeffs are not defined but the extmatelem are to be computed')")
      !stop 'FLinput - EXTF-coeffs are not defined but the extmatelem are to be computed'
      !
    end if
   end select
   !
   if (trim(job%IOextF_action)=='SAVE') then 
       if (trim(job%IOcontr_action)  =='NONE') job%IOcontr_action   = 'READ'
       if (trim(job%IOj0matel_action)=='NONE'.and.job%vib_contract) job%IOj0matel_action = 'READ'
   end if 
   !
   if (.not.symmetry_defined) then 
      !
      call SymmetryInitialize(job%sym_group)
      !
      symmetry_defined = .true.
      !
   endif
   !
   if (job%verbose>=6) then
     !
     call print_symmetries
     !    
   endif
   !
   if (trim(trove%symmetry)=='C2VN'.and.sym%N<job%bset(0)%range(2)) then
      write (out,"('FLinput: The C2VN number',i5,' must be defined and equal to (or <) krot',i5)") sym%N,job%bset(0)%range(2)
      stop 'FLinput - The C2VN number is undefined or too small'
   endif
    !
   if (.not.refer_defined) then 
      !
      trove%local_ref = trove%local_eq
      !
   endif
   !
   if (.not.basis_defined) then 
      !
      write (out,"('FLinput: The basis set is not defined')") 
      stop 'FLinput - basis set is not defined'
      !
   endif 
   !
   if (.not.chk_defined) then 
      !
      write (out,"('FLinput: The check_point section has not been defined')") 
      stop 'FLinput - The check_point section has not been defined'
      !
   endif 
   !
   if (.not.zmat_defined) then 
      !
      write (out,"('FLinput: ZMAT is not defined')") 
      stop 'FLinput - ZMAT is not defined'
      !
   endif 
   !
   if (.not.pot_defined) then 
      !
      write (out,"('FLinput: POTEN is not defined')") 
      stop 'FLinput - POTEN is not defined'
      !
   endif 
   !
   if (FLextF_matelem.and..not.extF_defined) then 
      !
      write (out,"('FLinput: External function has to be defined for external mat-elem calcs. ')") 
      stop 'FLinput - EXTF is not defined'
      !
   endif 
   !
   if (action%fitting.and..not.extF_defined.and..not.action%band_fitting) then 
      !
      write (out,"('FLinput: External function has to be defined for performing fitting ')") 
      stop 'FLinput - EXTF is not defined'
      !
   endif 
   !
   !
   if (.not.equil_defined) then 
      !
      write (out,"('FLinput: EQUIL is not defined')") 
      stop 'FLinput - EQUIL is not defined'
      !
   endif 
   !
   if (trim(job%PTtype)/='DIAGONAL'.and.manifold==1.and.NPTorder>0) then 
      !
      write (out,"('FLinput: PTtype  not compatible with non-rigid bender:',a)") trim(job%PTtype)
      stop 'FLinput - wrong PTtype'
      !
   endif
   !
   if (FLextF_coeffs.and..not.extF_defined) then 
      write (out,"('FLinput: External/Dipole field has to be specified when EXTERNAL checkpoint is used',a)")
      stop 'FLinput - External/Dipole field has to be specified for the EXTERNAL checkpoint'
   endif 
   !
   ! For the transformation to the j=0 basis set representation we require 
   ! to work only with all modes as one class in the contr. vibrational representaion, i.e.
   ! the vibr. Hamiltonian is assumed to be diagonal in this representaion:
   !
   if ( any( (/character(len=wl) :: trim(job%IOj0ext_action),trim(job%IOj0matel_action)/) /='NONE' ) ) then 
      !
      job%vib_contract = .true.
      !
      do i=1,Nmodes
         job%bset(i)%class = 1
      enddo
      !
      if (any(trim(job%IOj0ext_action) == (/'READ'/) ) ) then 
        job%extFmat_file  = job%exteigen_file
        job%extmat_suffix = job%j0extmat_suffix
      endif 
      !
      if (trim(job%IOj0matel_action) == 'READ') then 
        job%kinetmat_file = job%kineteigen_file
        job%matelem_suffix = job%j0matelem_suffix
      endif 
      !
      if (trim(job%IOeigen_action)=='SAVE'.and.action%convert_vibme) then 
        job%IOeigen_action = 'NONE'
        write(out,"('FLReadInput: It is illegal to save eigenvectors during the J=0 convertion; EIGENFUNC changed to NONE')")
      endif
      !
      if (trim(job%IOcontr_action)=='READ'.and.action%convert_vibme) then 
        job%IOcontr_action = 'NONE'
        write(out,"('FLReadInput: It is illegal to <CONTRACT READ> during the J=0 convertion; CONTRACT changed to NONE')")
      endif
      !
      if (trim(job%IOeigen_action)=='SAVE'.or.trim(job%IOeigen_action)=='APPEND'.or.trim(job%IOj0matel_action)=='READ'.or.&
               trim(job%IOj0ext_action)=='READ') then 
        !
        !if (trim(job%eigenfile%filebase)/='eigen') then 
          !
          job%eigenfile%filebase   = 'j0eigen'
          job%eigenfile%dscr       = 'j0eigen_descr'
          job%eigenfile%primitives = 'j0eigen_quanta'
          job%eigenfile%vectors    = 'j0eigen_vectors'
          !if (job%IOvector_symm) job%eigenfile%vectors    = 'j0eigen_vectors'
          !
        !endif 
        !
        job%contrfile%dscr       = 'j0'//trim(job%contrfile%dscr)
        job%contrfile%primitives = 'j0'//trim(job%contrfile%primitives)
        job%contrfile%vectors    = 'j0'//trim(job%contrfile%vectors)
        job%contrfile%dvr        = 'j0'//trim(job%contrfile%dvr)
        !
      endif 
      !
   endif
   !
   if ( job%convert_model_j0.and.job%bset(0)%range(1)>0 ) then
     write (out,"('EIGENFUNC SAVE CONVERT cannot be used with jrot>0')") jrot
     stop 'EIGENFUNC SAVE CONVERT cannot be used with jrot>0'
   endif
   !
   if ( trove%triatom_sing_resolve .and. (trove%Natoms/=3 .or. trove%Nmodes/=3 ) ) then
     write(out,"('Input error: LEGENDRE or SINRHO are currently only working with Nmodes=Natoms=3')") 
     stop 'Illegal usage of LEGENDRE or SINRHO'
   endif
   !
   trove%jmax = jrot 
   !
   write(char_j,"(i4)") jrot
   !
   job%eigenfile%dscr       = trim(job%eigenfile%dscr)//trim(adjustl(char_j))       !//'.chk'
   job%eigenfile%primitives = trim(job%eigenfile%primitives)//trim(adjustl(char_j)) !//'.chk'
   job%eigenfile%vectors    = trim(job%eigenfile%vectors)//trim(adjustl(char_j))    !//'.chk'
   !
   ! Check if everything defined 
   !
   if (job%verbose>=4) write(out,"('FLReadInput/end')")  
   !
   contains
   !
   subroutine print_symmetries
     integer(ik) :: igamma,iclass,ioper,ielem
     !
     write(out,"(/'Symmetry:',a)") trim(sym%group)
     !    
     write(out,"(/'Characters')")
     !
     do igamma = 1,sym%Nrepresen
       do iclass = 1,sym%Nclasses
          write(out,"(i4,1x,i4,1x,f16.8)") igamma,iclass,sym%characters(igamma,iclass)
      enddo 
     enddo
     !
     write(out,"(/'Irreps:')")
     !
     do igamma = 1,sym%Nrepresen
       !
       do ioper = 1,sym%Noper
          do ielem = 1,sym%degen(igamma)
            write(out,"(i4,1x,i4,1x,i4,1x,10f16.8)") igamma,ioper,ielem,sym%irr(igamma,ioper)%repres(ielem,:)
          enddo
       enddo 
       !
     enddo 
     !
   end subroutine print_symmetries
   !
end subroutine FLReadInput

subroutine check_read_save_none(w,place)
  !
  character(len=cl),intent(in) :: w,place
  !
  select case(trim(w))
  !
  case('NONE','SAVE','READ','SEPARATE','ASCII')
    !
    continue
    !
  case default 
    !
    write (out,"('FLinput: illegal key ',a,' in section ',a)") trim(w),trim(place)
    stop 'FLinput - illegal key '
    !
  end select 

end subroutine check_read_save_none


!
! Initilizing the molecule 
!
  subroutine FLsetMolecule

    !
    integer(ik) :: NPotOrder,NKinOrder,PotOrderShift,Natoms,Nmodes
    integer(ik) :: bonds(trove%Natoms-1,2)
    integer(ik) :: angles((trove%Natoms-3)*2+2,3)
    integer(ik) :: dihedrals(0:max(trove%Natoms,0),4) ! Dihedral Angles connections type 
    integer(ik) :: dihedtype(0:max(trove%Natoms,0))
    integer(ik) :: Ndihedrals ! number of dihedral angles of type 1 and type 2
    !
    integer(ik) :: alloc,io,ibond,n_t
    integer(ik) :: Nbonds,Nangles,i1
    !
    integer(ik) :: Kindex(trove%Nmodes),Nmodes_e,k1,k2,imode,iterm,jterm,dm2,irho,x1
    real(ark)   :: amorse,masses(trove%Natoms)
    real(ark)   :: ar_t(trove%Ncoords),a0_ark(trove%Natoms,3),chi(trove%Nmodes)
    real(ark)   :: b0_(trove%Natoms,3,0:0)
    real(ark)   :: rho_(0:0)
    real(ark)   :: rho_ref_
    real(ark)   :: rho_b_(2)
    !
    real(ark)   :: f(trove%Nmodes,trove%Nmodes)
    logical     :: dir
    real(ark)   :: step(2,trove%Nmodes),factor,df,rho_eq,Inertm(3)
    character(len=cl) :: my_fmt !format for I/O specification
    !
    if (job%verbose>=4) write(out,"(/'FLsetMolecule/start')")   
    !
    NKinOrder    = trove%NKinOrder
    NPotOrder    = trove%NPotOrder
    PotOrderShift= job%pot_pt_shift
    Nmodes       = trove%Nmodes
    Natoms       = trove%Natoms
    !
    ! convert z-matrix into the rbond,balpha,dalpha - connections
    !
    call zmat_to_bonds(bonds,angles,dihedrals,dihedtype,Nbonds,Nangles,Ndihedrals)
    !
    trove%Nbonds  = Nbonds
    trove%NAngles = NAngles
    trove%NDihedrals = Ndihedrals   
    !
    trove%Ncoords = Nbonds+Nangles+Ndihedrals
    !
    if (job%verbose>=2 .and. max(trove%Natoms-1,0)/=Ndihedrals ) then 
       write(out,"('Warning: number of dihedrals is not (Natoms-3): ',2i8)") Ndihedrals,max(trove%Natoms-1,0)
       !stop 'FLsetMolecule: Wrong number of dihedrals'
    endif
    ! 
    if (trove%Nmodes/=3*trove%Natoms-6.and.trove%Nmodes/=3*trove%Natoms-5) then 
       write(out,"('Warning: Number of modes, neither 3n-5, nor 3n-6, but ',i8)") trove%Nmodes
       stop 'FLsetMolecule: Wrong number of modes'
    endif 
    !
    ! define maximal value of NPotOrder and NKinOrder
    trove%MaxOrder = max(NPotOrder,NKinOrder,trove%NExtOrder)
    !trove%NPotOrder = trove%MaxOrder
    ! 
    ! Allocation of the molecular structure parameters matrixes:
    !
    allocate (trove%coordinates(3,Nmodes),trove%manifold_rank(Nmodes), &
              trove%bonds(trove%Nbonds,2),trove%angles(trove%NAngles,3), &
              trove%dihedtype(0:Ndihedrals),trove%dihedrals(0:Ndihedrals,4),stat=alloc)

    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate bonds and angles connections')") alloc
        stop 'FLsetMolecule, bonds and angles connections - out of memory'
    end if
    !
    !if (size(trove%bonds)/=size(bonds).or.size(trove%angles)/=size(angles) ) then 
    !    write(out,"('FLsetMolecule: ',i4,' inconsitent sizes of matrix atoms and angles in input: ',4i4)") &
    !                 ibond,trove%angles(ibond,:)
    !    stop 'FLsetMolecule - bad sizes of bonds and angles matrixes'
    !endif
    ! 
    !
    trove%bonds(1:Nbonds,:) = bonds(1:Nbonds,:)
    trove%angles(1:Nangles,:) = angles(1:Nangles,:)
    trove%dihedrals(:,:) = dihedrals(:,:)
    trove%dihedtype(:) = dihedtype(:)
    !
    ! We define the coordinates 
    !
    trove%Coordinates(1,:) = job%bset(1:Nmodes)%coord_kinet
    trove%Coordinates(2,:) = job%bset(1:Nmodes)%coord_poten
    trove%Coordinates(3,:) = extF%intcoords(1:Nmodes)
    !
    trove%sing_at_rho_0 = .false.
    !
    ! For the non-rigid configuration we reserve the last mode (imode=Nmodes) 
    ! for the large-amplitude motion, which will be treated explicitely 
    !
    select case(manifold)
    !
    case (0) ! rigid case, expansion around the rank=0 manifold
      !
      !     Nmodes_e - rigid motion modes, given as an expansion 
      trove%Nmodes_e  = Nmodes
      !
      !     Nmodes_n- non-rigid motion mode, given numerically
      trove%Nmodes_n  = 0
      trove%manifold_rank(:) = 0
      trove%npoints=0
      trove%rhostep = 0
      trove%rho_border(1:2) = 0
      trove%periodic = job%bset(Nmodes)%periodic
      !
    case (1) ! non-rigid case, expansion around a rank=1 manifold
      !
      trove%Nmodes_e  = Nmodes-1
      trove%Nmodes_n  = Nmodes
      trove%manifold_rank(:) = 0
      trove%manifold_rank(Nmodes) = 1
      trove%npoints   = job%bset(Nmodes)%npoints
      trove%rho_border(1:2) = job%bset(Nmodes)%borders(1:2)
      trove%rhostep  = (trove%rho_border(2)-trove%rho_border(1))/real(trove%npoints,kind=rk)
      trove%periodic = job%bset(Nmodes)%periodic
      !
    end select
    !
    if (any( trove%Coordinates(:,:)=='NORMAL').and..not. all(trove%Coordinates(:,:)=='NORMAL') ) then 
       write(out,"('Wrong defenition of coordinates:')")
       write(out,"('they must be either all normal or none')")
       write(out,"('The coordinates (first 30) are: ',30(a,2x))") ( trim(trove%Coordinates(1,i1)),i1=1,min(30,Nmodes) )
       write(out,"('The coordinates (first 30) are: ',30(a,2x))") ( trim(trove%Coordinates(2,i1)),i1=1,min(30,Nmodes) )
       stop 'FLsetMolecule: Wrong defenition of coordinates'
    endif 
    !
    if ( trove%Coordinates(1,1)=='NORMAL'.and.manifold/=0 ) then 
       write(out,"('Wrong defenition of coordinates:')")
       write(out,"('normal can be used only for the rigid bender case (mamifold=0)')")
       stop 'FLsetMolecule: normal is not a proper choice for the non-rigid bender'
    endif 
    !
    if ( Nmodes>1.and.(job%bset(Nmodes-1)%species==job%bset(Nmodes)%species) ) then 
        write(out,"('Wrong definition of the bst%species: ')")
        write(out,"('the last mode has to be separated from all others')")
        !stop 'FLsetMolecule: wrong definition of the last bst%species'
    endif
    !
    ! Check if all res_coeffs are not zero
    ! 
    if ( any(job%bset(1:Nmodes)%res_coeffs<small_) ) then 
        write(out,"('Wrong defenition of the bst%res_coeffs:')")
        write(out,"('They cannot be zero: ',30f12.2)") ( job%bset(i1)%res_coeffs,i1=1,min(30,Nmodes) )
        stop 'FLsetMolecule: wrong definition of the last bst%res_coeffs'
    endif
    !
    ! Check if the definition of bonds and angles is physical
    !
    do ibond = 1,trove%NAngles
       if (Nangles/=0) then 
          n_t = trove%angles(ibond,1)
          if (any(trove%angles(ibond,2:3)==n_t).or.trove%angles(ibond,2)==trove%angles(ibond,3)) then 
             write(out,"('FLsetMolecule: ',i4,' angle connections belong to the same atom: ',3i4)") & 
                          ibond,trove%angles(ibond,:)
             stop 'FLsetMolecule - bad angle connections'
          endif
       endif 
    enddo 

    do ibond = 1,trove%Nbonds
       if (trove%bonds(ibond,1)==trove%bonds(ibond,2)) then 
          write(out,"('FLsetMolecule: ',i4,' bond connections belong to the same atom: ',2i4)") & 
                       ibond,trove%bonds(ibond,:)
          stop 'FLsetMolecule - bad angle connections'
       endif
    enddo 
    !
    do ibond = 1,Ndihedrals
        n_t = trove%dihedrals(ibond,1)
        do io = 2,4 
          if (trove%dihedrals(ibond,io)==n_t) then
             if (any(trove%zmatrix(n_t)%connect(4)==(/101,102,103,104,105,106,107,108/))) cycle
             write(out,"('FLsetMolecule: ',i4,' dihedral angle connections belong to the same atom: ',3i4)") & 
                          ibond,trove%dihedrals(ibond,:)
             stop 'FLsetMolecule - bad dihedrals angle connections'
          endif
        enddo
    enddo
    ! 
    ! Allocation of the molecular structure parameters matrixes:
    !
    allocate (trove%a0(Natoms,3), & 
              trove%Amatrho(Natoms,3,Nmodes,0:trove%Npoints),  &
              trove%dAmatrho(Natoms,3,Nmodes,0:trove%Npoints,3),  &
              trove%Bmatrho(Nmodes,Natoms,3,0:trove%Npoints),  &
              trove%dBmatrho(Nmodes,Natoms,3,0:trove%Npoints,2),  &
              trove%req(Nbonds),trove%chi_ref(Nmodes,0:trove%Npoints),trove%chi_eq(Nmodes), &
              trove%chi0_ref(Nmodes),trove%alphaeq(NAngles),trove%taueq(Ndihedrals),stat=alloc)


    call ArrayStart('trove%a0'      ,alloc,size(trove%a0)      ,kind(trove%a0))
    call ArrayStart('trove%Amatrho' ,alloc,size(trove%Amatrho) ,kind(trove%Amatrho))
    call ArrayStart('trove%dAmatrho',alloc,size(trove%dAmatrho),kind(trove%dAmatrho))
    call ArrayStart('trove%Bmatrho' ,alloc,size(trove%Bmatrho) ,kind(trove%Bmatrho))
    call ArrayStart('trove%dBmatrho',alloc,size(trove%dBmatrho),kind(trove%dBmatrho))
    !
    trove%req(1:Nbonds) = trove%local_eq(1:Nbonds)
    trove%alphaeq(1:Nangles) = trove%local_eq(Nbonds+1:Nbonds+Nangles)
    trove%taueq(1:Ndihedrals) = trove%local_eq(Nbonds+Nangles+1:trove%Ncoords)
    !
    ! Allocation of the molecular structure parameters matrix in the case of manifold rank=1
    !
    allocate (trove%b0(Natoms,3,0:trove%npoints), &
              trove%db0(Natoms,3,0:trove%npoints,3),&
              trove%rho_i(0:trove%Npoints),stat=alloc)

    call ArrayStart('trove%b0',alloc,size(trove%b0),kind(trove%b0))
    call ArrayStart('trove%db0',alloc,size(trove%db0),kind(trove%db0))
    call ArrayStart('trove%rho_i',alloc,size(trove%rho_i),kind(trove%rho_i))
    !
    ! be verbose
    ! 
    if (job%verbose>=3) then
      write(out,"('/Molecular parameters:')")
      write(out,"('Number of atoms (Natoms):',i5)")   Natoms
      write(out,"('Number of bonds (Nbonds):',i5)")   Nbonds
      write(out,"('Number of angles (Nangles):',i5)") Nangles
      write(out,"('Number of dihedral angles (Ndihedrals):',i5)") Ndihedrals
      !
      write(out,"('The molecular type (Moltype):',a)") trim(trove%Moltype)
      !
      write(my_fmt,'(a,i0,a)') "(a,",Nmodes,"(1x,a))"
      !
      write(out,my_fmt) 'The kinetic   coordinates (Molecule):',( trim(trove%Coordinates(1,i1)),i1=1,Nmodes)
      write(out,my_fmt) 'The potential coordinates (Molecule):',( trim(trove%Coordinates(2,i1)),i1=1,Nmodes)
      write(out,my_fmt) 'The externalF coordinates (Molecule):',( trim(trove%Coordinates(3,i1)),i1=1,Nmodes)
      !
      write(out,"('fdstep:',40f14.6)") trove%fdstep(1:min(40,trove%Nmodes))
      !
      write(out,"('The coordinates type (internal_coords):',a)") trim(trove%internal_coords)
      !
      write(out,"('The coordinate-transformation type (coords_transform):',a)") trim(trove%coords_transform)
      !
      write(out,"('Bonds  connections are ')")  
      do ibond = 1,trove%Nbonds
       write(out,"(i4,' - ',2i4)") ibond,trove%bonds(ibond,:)
      enddo 
      !
      if (Nangles/=0) then 
         write(out,"('Angles  connections are ')")  
         do ibond = 1,trove%Nangles
            write(out,"(i4,' - ',3i4)") ibond,trove%angles(ibond,:)
         enddo 
      endif
      !
      !
      if (Nangles/=0) then 
         write(out,"('Angles  connections are ')")  
         do ibond = 1,trove%Nangles
            write(out,"(i4,' - ',3i4)") ibond,trove%angles(ibond,:)
         enddo 
      endif
      !
      if (trove%Npoints/=0) then 
         write(out,"('The coordinate ',i3,' is given explicitly by 1d numerical table of size',i8)") Nmodes,trove%Npoints
      endif
      !
      write(out,"(/'Calculations control parameters:')")
      write(out,"('Maximal order of the kinetic   energy expansion (NKinOrder):',i5)") NKinOrder
      write(out,"('Maximal order of the potential energy expansion (NPotOrder):',i5)") NPotOrder
      !
      write(out,"('A small number:',d20.10)") small_
      !
      !
    endif 
    !
    Nmodes_e = trove%Nmodes_e
    !
    ! define the ranges of the index within the array field betwen the different orders
    !
    allocate (trove%RangeOrder(0:trove%MaxOrder+2),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate RangeOrder')") alloc
        stop 'FLsetMolecule, RangeOrder - out of memory'
    end if
    !
    ! Here we estimate the maximal number of coefficients in all kind of fields
    ! FLQindex will compute just it, when the second optional argument is omitted
    !
    call TimerStart('FLQindex-1')
    !
    do io = 0,trove%MaxOrder+2
      !
      Kindex = 0 
      Kindex(1) = io 
      trove%RangeOrder(io) = FLQindex(Nmodes_e,Kindex)
      !
    enddo
    call TimerStop('FLQindex-1')
    !
    trove%Ncoeff = trove%RangeOrder(trove%MaxOrder)
    io = trove%RangeOrder(trove%MaxOrder+2)
    !
    ! Definition of the indexing: relations between the Nmode-D and 1D arrays 
    ! are stored in FLIndexQ (forward) and deifned by Qindex routine (backward)
    !
    allocate (FLIndexQ(trove%Nmodes,io),stat=alloc)
    call ArrayStart('FLIndexQ',alloc,size(FLIndexQ),kind(FLIndexQ))
    !
    !allocate (FLIndexQ_legatee(trove%Nmodes,io),stat=alloc)
    !call ArrayStart('FLIndexQ_legatee',alloc,size(FLIndexQ_legatee),kind(FLIndexQ_legatee))

    !allocate (mat_legatee(trove%Nmodes,io),stat=alloc)
    !call ArrayStart('mat_legatee',alloc,size(mat_legatee),kind(mat_legatee))
    !
    ! Final call: defining the matrix FLIndexQ
    !
    Kindex = 0 ; Kindex(1) = trove%MaxOrder+2
    !
    call TimerStart('FLQindex-2')
    !
    io = FLQindex(Nmodes_e,Kindex,FLIndexQ)
    !
    call TimerStop('FLQindex-2')
    !
    ! Defining FLindexQ_legatee
    !
    !FLindexQ_legatee(:,:) = -1 ! no inheritance
    !
    !dm2 = size(FLIndexQ,dim=2)
    !
    ! Searching inheritance with maximum initial coincidence 
    !
    !do imode=trove%Nmodes-1,1,-1
    !  !
    !  do iterm=1,dm2-1
    !    !
    !    do jterm=iterm+1,dm2
    !      !
    !      if(FLindexQ_legatee(1,jterm)==-1.and.all(FLIndexQ(1:imode,iterm)==FLIndexQ(1:imode,jterm))) then
    !          !
    !          FLindexQ_legatee(1:imode,jterm)=0 ! We do not need to calculate corresponding part in mat_element
    !          FLindexQ_legatee(imode+1,jterm)=iterm ! The point of starting calculation from another stream
    !          FLindexQ_legatee(imode+2:trove%Nmodes,jterm)=jterm ! Calculate in an ordinar stream
    !          !
    !      endif
    !    enddo ! jterm
    !    !
    !  enddo ! iterm
    !  !
    !enddo ! imode
    !
    !do jterm=1,dm2
    !  !
    !  if (FLindexQ_legatee(1,jterm)==-1) then 
    !     !
    !     FLindexQ_legatee(2:trove%Nmodes,jterm)=jterm
    !     FLindexQ_legatee(1,jterm) = 1
    !     !
    !  elseif (FLindexQ_legatee(1,jterm)==0) then 
    !     !
    !     FLindexQ_legatee(1,jterm) = 1
    !     !
    !     do imode=2,trove%Nmodes
    !       !
    !       if (FLindexQ_legatee(imode,jterm)/=0) then
    !          FLindexQ_legatee(1,jterm) = imode
    !          exit
    !       endif 
    !       !
    !     enddo
    !     !
    !  endif 
    !  !
    !enddo
    !
    !do jterm=1,dm2
    !  !
    !  imode = 1
    !  !
    !  do while(imode<trove%Nmodes.and.FLindexQ_legatee(imode,jterm)==0)
    !    !
    !    imode = imode + 1
    !    !
    !  enddo
    !  !
    !  FLindexQ_legatee(1,jterm)=imode
    !  !
    !enddo
    !
    !if (job%verbose>=7) then
    !  !
    !  write(out,"(/'Inheritance matrix')")
    !  do iterm = 1,dm2 
    !    write(out,"(2i8,30i8)") io,iterm,(FLindexQ_legatee(io,iterm),io=1,min(30,trove%Nmodes))
    !    !
    !  enddo
    !  !
    !endif 
    !
    !
    ! be verbose
    ! 
    !if (job%verbose>=1) then
    !  !
    !  !write(out,"(/'RangeOrder vs order ',30i8 )") (io               ,io=0,min(30,trove%MaxOrder))
    !  !write(out,"( 'RangeOrder:         ',30i8/)") (trove%RangeOrder(io),io=0,min(30,trove%MaxOrder))
    !  !
    !endif
    !
    if (job%verbose>=7) then
      !
      write(out,"(/'FLIndexQ:')")
      do k1=1,size(FLIndexQ,dim=2)
         write(out,"(30i8)") (FLIndexQ(io,k1),io=1,min(30,trove%Nmodes))
      enddo
      !
    endif
    !
    ! Allocation of the quadratic force constants matrix qwforce 
    !
    allocate (trove%omega(Nmodes),trove%coord_f(Nmodes),stat=alloc)
    call ArrayStart('trove%omega',alloc,size(trove%omega),kind(trove%omega))
    call ArrayStart('trove%omega',alloc,size(trove%coord_f),kind(trove%coord_f))
    !
    trove%coord_f = 1
    !
    allocate (trove%qwforce(trove%Nmodes,trove%Nmodes,0:trove%npoints),stat=alloc)
    call ArrayStart('trove%qwforce',alloc,size(trove%qwforce),kind(trove%qwforce))
    !
    ! define all molecular paramerers here 
    !
    call MLinitialize_molec(  trove%Moltype,trove%Coordinates,trove%coords_transform,&
                              trove%Nbonds,trove%Nangles,trove%Ndihedrals,&
                              trove%dihedtype,&
                              trove%mass,trove%local_eq,&
                              force,forcename,ifit,pot_ind,trove%specparam,trove%potentype,trove%kinetic_type,&
                              trove%symmetry,trove%rho_border,trove%zmatrix)
    !
    ! define the potential function method
    !
    call MLdefine_potenfunc
    !
    ! define the dipole function method
    !
    call MLextF_func_define
    !
    ! define the kinetic energy method
    !
    call MLdefine_kinetic_subroutine
    !
    ! define the coordinate-transformation method
    !
    call MLcoordinate_transform_func_define
    !
    deallocate(force,forcename,pot_ind,ifit)
    !
    ! Morse coordinates expansion is a special case: amorse has to be defined, for example. 
    ! if it is not - we stop the calculations 
    !
    do k1=1,NModes
       !
       if ( trim(trove%Coordinates(1,k1))=='MORSE'.or.trim(trove%Coordinates(2,k1))=='MORSE' ) then 
          !
          !amorse = trove%specparam(k1)
          !
          if (all(trove%specparam(:)<small_)) then 
            !
            amorse=1.0_rk   ! alternative option 
            !
            write (out,"('FLsetMoleculet: amorse is undefined')")
            stop 'FLsetMolecule - bad amorse'
            !
          endif 
          !
        endif 
       !
    enddo 
    !
    ! masses for the internal use 
    !
    masses = trove%mass

    !
    ! define the equilibrium chi parameters 
    !
    ar_t = trove%local_eq 
    !
    dir = .true.
    !
    trove%chi_eq(:) = MLcoordinate_transform_func(ar_t,Nmodes,dir)
    !
    call MLequilibrium_chi(trove%chi_eq(:))
    !
    trove%chi_ref(:,0) = trove%chi_eq(:)
    !
    ! define the reference chi parameters 
    !
    ar_t = trove%local_ref
    !
    dir = .true.
    !
    trove%chi0_ref(:) = MLcoordinate_transform_func(ar_t,Nmodes,dir)
    !
    ! reference geometry for expansion of the external function (if not given in input)
    !
    if (all(extF%geom_ref(:)<sqrt(small_))) extF%geom_ref = trove%local_eq
    !
    ! define the equilibrium Cartesian geometry a0
    !
    b0_(:,:,0) = trove%a0(:,:)
    call MLequilibrium_xyz(0_ik,molec%Natoms,b0_)
    trove%a0 = b0_(:,:,0)
    !
    ! check the definition of the local generalized coordinates
    !
    a0_ark = trove%a0
    !
    call FLfromcartesian2local(a0_ark,ar_t)
    !
    do imode = 1,trove%Ncoords
       if (abs(ar_t(imode)-trove%local_eq(imode))>1e-6) then
          !
          if (trove%local_eq(imode)==2.0_ark*pi.and.abs(ar_t(imode))<sqrt(small_))  cycle
          !
          write(out,"('FLsetMolecule-ERROR:')")
          write(out,"('Input eq. value ',f18.12,' for ',i4,'-th coordinate is different from the calculated ',f18.12)") & 
                       trove%local_eq(imode),imode,ar_t(imode)
          write(out,"('You can copy and paste this number into the input.')")
          stop 'FLsetMolecule - equilibrium or a0 are wrong'
       endif
    enddo 
    !
    ! Quadratic force constants that will be needed to construct the Amat constants 
    ! will be calculated in terms of the internal coordinates chi
    !
    ! Here we proceed with the finite differences for the qudratic potential fucntion at the equilibrium
    !
    f = 0
    !
    step(1,:) = trove%fdstep(:)
    step(2,:) = trove%fdstep(:)
    !
    do k1=1,trove%Nmodes
       do k2=1,k1
          !
          kindex = 0 ; kindex(k1) = kindex(k1)+1 ; kindex(k2) = kindex(k2)+1
          ! 
          ! Factorial factors to convert from 1/cm to aJ
          !
          factor = real(1.0e11_rk*planck*vellgt,ark)
          !
          ! Here we calculate derivatives by the finite differences 
          ! of the function "poten_local" (or poten_xi)
          ! with respect to local coordinates r1^k1 r2^k2 r3^k3 ...
          ! at the equilibrium given by q_eq,
          ! while k1,k2,k3... are stored in "kindex".
          ! fdstep defines the finite differences spacings 
          !
          df = FLfinitediffs(kindex,poten_chi,trove%chi_ref(:,0),step) 
          !
          f(k1,k2) = df*factor
          f(k2,k1) = df*factor
          !
       enddo 
    enddo 
    !
    !f = matmul(matmul((trove%coordtransform),f),transpose(trove%coordtransform))
    !
    f = f*0.5_ark
    !
    trove%qwforce(:,:,0) = f
    !
    if (job%verbose >= 5) then 
        write(out,"('Normal quadratic pot. parameteres   qwforce:')")
        !
        do k1 = 1,Nmodes
           do k2 = 1,k1
              write(out,"(20x,2i5,d18.8)") k1,k2,trove%qwforce(k1,k2,0)
           enddo
        enddo
        !
    endif 
    !
    ! First we calculate the equilibrium values for the internal coordinates chi
    !
    !
    ! Here we do similar cacluations of the quadratic potential function but for 
    !
    ! Here we compute the structure parameter matrix b0 in case of manifold rank=1
    !
    select case (manifold)
      !
    case (0)
      !
      ! The trick is to use for the 0D case the same routines developed for rhe 1D case
      !
      trove%b0(:,:,0) = trove%a0(:,:)
      trove%rho_border = trove%chi_ref(trove%Nmodes,0)
      trove%rho_i(0) = trove%chi_ref(trove%Nmodes,0)
      !
      trove%ipotmin = 0 
      !
      extF%irho_ref= 0
      !
      ! Determine is this is a linear-type molecule and the singular axis lincoord = x,y,z
      !
      !trove%lincoord = 0
      !
      !if (trove%lincoord==0) then
      !  if (all(trove%a0(:,1)==0.0_rk).and.all(trove%a0(:,2)==0.0_rk)) trove%lincoord = 3
      !  if (all(trove%a0(:,2)==0.0_rk).and.all(trove%a0(:,3)==0.0_rk)) trove%lincoord = 1
      !  if (all(trove%a0(:,1)==0.0_rk).and.all(trove%a0(:,3)==0.0_rk)) trove%lincoord = 2
      !endif
      !
      Inertm(1) = sum( trove%mass(:)*( trove%a0(:,2)**2+ trove%a0(:,3)**2 ) )
      Inertm(2) = sum( trove%mass(:)*( trove%a0(:,1)**2+ trove%a0(:,3)**2 ) )
      Inertm(3) = sum( trove%mass(:)*( trove%a0(:,1)**2+ trove%a0(:,2)**2 ) )
      !
      ! Check if all Inertia moments are non zero 
      !
      if (trove%lincoord==0) then
        do x1 = 1,3
          if (Inertm(x1)<sqrt(small_)) then 
             trove%lincoord=x1
             job%lincoord = trove%lincoord
          endif
        enddo
      endif
      !
      ! This can be a linear-type molecule also when the number of modes = 3N-5
      !
      if (Nmodes==(3*Natoms-5).and.trove%lincoord==0) then
          trove%lincoord = minloc(Inertm,dim=1)
          job%lincoord = trove%lincoord
      endif
      !
    case (1)
      !
      ar_t = extF%geom_ref
      dir = .true.
      chi(:) = MLcoordinate_transform_func(ar_t,trove%Nmodes,dir)
      extF%irho_ref = mod(nint( ( chi(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
      !
      rho_eq = trove%chi_ref(Nmodes,0)
      trove%ipotmin = mod(nint(( rho_eq-trove%rho_border(1) )  &
                         /trove%rhostep,kind=ik)+trove%npoints,trove%npoints)
      dir = .true.
      !
      call MLequilibrium_xyz_1d(trove%Npoints,trove%rho_border,trove%rhostep,trove%periodic,trove%rho_ref,trove%b0,&
                                trove%db0,trove%rho_i)
      !
      if (job%verbose>=6)  write(out,"(2x,'rho values:')") 
      !
      do irho=0,trove%Npoints
        !
        !write(out,"('irho = ',i8)") irho
        !
        a0_ark(:,:) = trove%b0(:,:,irho)
        !
        call FLfromcartesian2local(a0_ark,ar_t)
        !
        trove%chi_ref(:,irho) = MLcoordinate_transform_func(ar_t,Nmodes,dir)
        !
        if (job%verbose>=6)  write(out,"(i5,1x,f15.7)") irho,trove%chi_ref(trove%Nmodes,irho)*180.0_ark/pi
        !
        ! Check for the linearity 
        !
        Inertm(1) = sum( trove%mass(:)*( trove%b0(:,2,irho)**2+ trove%b0(:,3,irho)**2 ) )
        Inertm(2) = sum( trove%mass(:)*( trove%b0(:,1,irho)**2+ trove%b0(:,3,irho)**2 ) )
        Inertm(3) = sum( trove%mass(:)*( trove%b0(:,1,irho)**2+ trove%b0(:,2,irho)**2 ) )
        !
        ! Check if all Inertia moments are non zero
        !
        !write(out,"('irho = ',i8,3f16.8)") irho, Inertm
        !
        if (trove%lincoord==0) then
          do x1 = 1,3
            if (Inertm(x1)<sqrt(small_)) trove%lincoord=x1
          enddo
        endif
        !
      enddo 
      !
    end select
    !
    !if (trim(trove%IO_hamiltonian)=='SAVE') call FLcheck_point_Hamiltonian('AMAT_SAVE')
    !
    ! the rank = 1 manifold case
    !
    if (manifold/=0.and.trove%Coordinates(1,1)=='NORMAL') then
       !
       write (out,"('FLsetMolecule: Normal  is not for the 1d non-rigid case')")
       stop 'FLsetMolecule: normal is not used for non-rigid'
       !
    endif 
    !
    !
    if (job%verbose>=4) write(out,"('FLsetMolecule/end')")   
    !
  end subroutine FLsetMolecule
!
!
! convert z-matrix into the rbond,balpha,dalpha - connections
!
  subroutine zmat_to_bonds(bonds,angles,dihedrals,dihedtype,Nbonds,Nangles,Ndihedrals)
    !
    integer(ik),intent(out) :: bonds(trove%Natoms-1,2)
    integer(ik),intent(out) :: angles((trove%Natoms-3)*2+2,3)
    integer(ik),intent(out) :: dihedrals(0:max(trove%Natoms,0),4) ! Dihedral Angles connections type 
    integer(ik),intent(out) :: dihedtype(0:max(trove%Natoms,0))
    integer(ik),intent(out) :: Ndihedrals     ! number of dihedral angles
    integer(ik),intent(out) :: Nbonds,Nangles
    integer(ik) ::  Natoms,iatom,J,kappa,zeta
    !
    Natoms = trove%Natoms
    !
    do iatom = 2,Natoms
       bonds(iatom-1,1) = iatom 
       bonds(iatom-1,2) = trove%zmatrix(iatom)%connect(1)
    enddo
    !
    Nbonds  = Natoms-1
    angles  = 0 
    nangles = 0 
    dihedtype = 0
    !
    Ndihedrals = 0
    dihedrals = 0 
    !
    do iatom = 3,Natoms
       !
       nangles = nangles + 1
       !
       angles(nangles,1) = iatom 
       angles(nangles,2) = trove%zmatrix(iatom)%connect(1)
       angles(nangles,3) = trove%zmatrix(iatom)%connect(2)
       !
       if (iatom>=4.or.trove%zmatrix(iatom)%connect(3)/=0) then
          !
          J = trove%zmatrix(iatom)%connect(4)
          !
          select case (J) 
            !
          case default
             !
             write (out,"('zmat_to_bonds: illegal J value in Z-matrix for atom',i8,': ',i8)") iatom,J
             stop 'zmat_to_bonds - illegal J value in Z-matrix'
             !
          case(-1,0)
             !
             ! J=0 -> beta = alpha(p0,p1,p3)
             !
             NAngles = NAngles + 1 
             !
             angles(NAngles,1) = iatom 
             angles(NAngles,2) = trove%zmatrix(iatom)%connect(1)
             angles(NAngles,3) = trove%zmatrix(iatom)%connect(3)
             !
          case(1) 
             !
             ! J=1 -> beta = dihedral-beta(p0,p1,p2,p3) - type 1
             !
             NAngles = NAngles + 1 
             !
             angles(NAngles,1) = iatom 
             angles(NAngles,2) = trove%zmatrix(iatom)%connect(1)
             angles(NAngles,3) = trove%zmatrix(iatom)%connect(3)
             Ndihedrals = Ndihedrals + 1 
             !
             dihedtype(Ndihedrals) = J
             !
             dihedrals(Ndihedrals,1) = iatom 
             dihedrals(Ndihedrals,2) = trove%zmatrix(iatom)%connect(1)
             dihedrals(Ndihedrals,3) = trove%zmatrix(iatom)%connect(2)
             dihedrals(Ndihedrals,4) = trove%zmatrix(iatom)%connect(3)
             !
          case(-2,2,-202,202,-302,302,-402,402) 
             !
             ! J=2 -> beta = dihedral-beta(p0,p1,p2,p3) - type 2
             !
             Ndihedrals = Ndihedrals + 1 
             !
             dihedtype(Ndihedrals) = J
             !
             dihedrals(Ndihedrals,1) = iatom 
             dihedrals(Ndihedrals,2) = trove%zmatrix(iatom)%connect(1)
             dihedrals(Ndihedrals,3) = trove%zmatrix(iatom)%connect(2)
             dihedrals(Ndihedrals,4) = trove%zmatrix(iatom)%connect(3)
             !
          case(3:100) 
             !
             ! J>2 -> beta1 = alpha(p0,p1,p3)
             !
             NAngles = NAngles + 1 
             !
             angles(NAngles,1) = iatom 
             angles(NAngles,2) = trove%zmatrix(iatom)%connect(1)
             angles(NAngles,3) = trove%zmatrix(iatom)%connect(3)
             !
             ! J>2 -> beta2 = alpha(p0,p1,J)
             !
             NAngles = NAngles + 1 
             !
             angles(NAngles,1) = iatom 
             angles(NAngles,2) = trove%zmatrix(iatom)%connect(1)
             angles(NAngles,3) = trove%zmatrix(iatom)%connect(4)
             !
             !
          case(101,103,104,105,106,107,108) 
             !
             ! special case of angles for a linear molecule
             !
             NAngles = NAngles - 1
             !
             ! trove%zmatrix(iatom)%connect(3) is a Cartesian component the special agnles are build around, 
             ! e.g. for ..ct(3)=z, the  x and y ref-vectors will be used to define the special angles as 
             ! (i dot [e1xe2])
             !
             zeta = trove%zmatrix(iatom)%connect(3)
             !
             if (all(zeta/=(/1,2,3/))) then
               !
               write (out,"('zmat_to_bonds: illegal zeta = ',i4,' of 3d dihedral for the linear angle of the atom ',i4,'  ')") &
                      kappa,iatom
               stop 'zmat_to_bonds - illegal zeta'
               !
             endif
             !
             do kappa = 1,3
               !
               if (kappa==zeta) cycle
               !
               Ndihedrals = Ndihedrals + 1
               !
               dihedtype(Ndihedrals) = J
               !
               dihedrals(Ndihedrals,1) = iatom 
               dihedrals(Ndihedrals,2) = trove%zmatrix(iatom)%connect(1)
               dihedrals(Ndihedrals,3) = trove%zmatrix(iatom)%connect(2)
               dihedrals(Ndihedrals,4) = kappa
               !
             enddo
             !
          case(102) 
             !
             ! special case of an angle for a linear molecule
             !
             NAngles = NAngles - 1
             !
             ! trove%zmatrix(iatom)%connect(3) is a Cartesian component the special agnles are build around, 
             ! e.g. for ..ct(3)=z, the  x and y ref-vectors will be used to define the special angles as 
             ! (i dot [e1xe2])
             !
             zeta = trove%zmatrix(iatom)%connect(3)
             !
             if (all(zeta/=(/1,2,3/))) then
               !
               write (out,"('zmat_to_bonds: illegal zeta = ',i4,' of 3d dihed for the linear angle of atom ',i4,'  ')") kappa,iatom
               stop 'zmat_to_bonds - illegal zeta'
               !
             endif
             !
             do kappa = 1,3
               !
               if (kappa==zeta) cycle
               !
               Ndihedrals = Ndihedrals + 1
               !
               dihedtype(Ndihedrals) = J
               !
               dihedrals(Ndihedrals,1) = iatom 
               dihedrals(Ndihedrals,2) = trove%zmatrix(iatom)%connect(1)
               dihedrals(Ndihedrals,3) = trove%zmatrix(iatom)%connect(2)
               dihedrals(Ndihedrals,4) = kappa
               !
             enddo
             !
             !kappa=zeta
             !
             !Ndihedrals = Ndihedrals + 1
             !
             !dihedtype(Ndihedrals) = J
             !
             !dihedrals(Ndihedrals,1) = iatom 
             !dihedrals(Ndihedrals,2) = trove%zmatrix(iatom)%connect(1)
             !dihedrals(Ndihedrals,3) = trove%zmatrix(iatom)%connect(2)
             !dihedrals(Ndihedrals,4) = kappa
             !
          end select 
          !
       endif
       !
    enddo

  end subroutine zmat_to_bonds

  !
  ! Here we initialize the kinetic operator fields: g_vib, g_rot, and g_cor
  !
  subroutine FLinitilize_Kinetic
    !
    integer(ik) :: alloc,Tcoeff,Tcoeff1,k1,k2,iatom,imode,i,ix
    integer(ik) :: Nmodes,Natoms,Nbonds,Nangles,Ndihedrals,io
    integer(ik) :: Kindex(trove%Nmodes),irho,npoints,Nmodes_e,NKinOrder_
    type(FLpolynomT),pointer     :: fl
    type(FLpolynomT),pointer     :: s_vib(:,:,:),s_rot(:,:,:)
    character(len=2)             :: txt1,txt2

    real(ark)                    :: hstep,diferror

    !real(ark)                     :: masses(trove%Natoms)
    !real(ark)                     :: df,factor,amorse
    !real(ark)                     :: ar_t(trove%Ncoords),a0_ark(trove%Natoms,3)
    !real(ark)                     :: step(2,trove%Nmodes),f(trove%Nmodes,trove%Nmodes)
    character(len=cl)             :: dir,fname
    !
    if (job%verbose>=4) write(out,"(/'FLinitilize_Kinetic/start')")   
    !
    ! If the kinetic operator fields have been stored we can just read them from the hard disk and leave...
    !
    if (trim(trove%IO_hamiltonian)=='READ'.or.&
        !trim(trove%IO_potential)=='READ'.or.&
        trim(trove%IO_kinetic)=='READ') then 
        !
        call FLcheck_point_Hamiltonian('KINETIC_READ')
        !
        call print_kinetic
        !
        if (trove%sparse) call compact_sparse_kinetic
        !
        return 
        !
    endif
    !
    call TimerStart('Kinetic')
    ! 
    ! Parameters for the internal use 
    !
    Nmodes  = trove%Nmodes
    Natoms  = trove%Natoms
    Nbonds  = trove%Nbonds
    NAngles = trove%Nangles
    Npoints = trove%Npoints
    Ndihedrals = trove%Ndihedrals
    !
    Nmodes_e = trove%Nmodes_e
    !
    ! define the ranges of the index within the array field betwen the different orders
    !
    ! Definition of 
    ! s_vib vector = (d xi_l / d r_Na )           !
    ! s_rot vector = (d rot_angles / d r_Na )     !
    ! r_na   vector = r_na - cartesian coordinates !
    !
    !
    ! Calculate Ncoeff: number of elements in the kinetic fields for NkinOrder
    !
    Tcoeff  = trove%RangeOrder(trove%NKinOrder+2)
    Tcoeff1 = trove%RangeOrder(trove%NKinOrder+1)
    !do k1 = 1,size(s_rot)
    !   fl => s_rot(k1)
    !   call polynom_initialization(fl,trove%NKinOrder+2,Tcoeff,Npoints,'s_rot')
    !enddo
    !
    if (trove%internal_coords/='LOCAL') then
      !
      ! Allocation 
      !
      allocate (s_vib(Nmodes,Natoms,3),s_rot(3,Natoms,3),stat=alloc)
      if (alloc/=0) then
          write (out,"(' Error ',i9,' trying to allocate s_vib, s_rot-fields')") alloc
          stop 'FLinitilize_Kinetic, s_vib, s_rot-fields -  out of memory'
      end if
      !
      ! Vibrational vector s_vib : 
      !
      do imode = 1,Nmodes
        do iatom = 1,Natoms
          do ix = 1,3
             fl => s_vib(imode,iatom,ix)
             call polynom_initialization(fl,trove%NKinOrder+2,Tcoeff,Npoints,'s_vib')
          enddo
        enddo
      enddo
      !
      ! Rotational vector s_rot : 
      !
      do imode = 1,3
        do iatom = 1,Natoms
          do ix = 1,3
             fl => s_rot(imode,iatom,ix)
             call polynom_initialization(fl,trove%NKinOrder+1,Tcoeff1,Npoints,'s_rot')
          enddo
        enddo
      enddo
      !
      ! Vibrational angular momentum
      !
      if (FLl2_coeffs) then
         !
         allocate (trove%L2_vib(Nmodes,Nmodes),stat=alloc)
         if (alloc/=0) then
             write (out,"(' Error ',i9,' trying to allocate L2_vib')") alloc
             stop 'FLinitilize_Kinetic, L2_vib - out of memory'
         end if
         !
         NKinOrder_ = 2
         !
         Tcoeff  = trove%RangeOrder(2)
         !
         do k1 = 1,Nmodes
            do k2 = 1,Nmodes
               fl => trove%L2_vib(k1,k2)
               call polynom_initialization(fl,NKinOrder_,Tcoeff,Npoints,'L2_vib')
            enddo
         enddo
         !
      endif      
      !
      ! Generate the Amat/Lmat matrix with first derivatives of cartesian coordinates 
      ! with respect to the local linearized coordinates or thier combinations 
      ! Lmat is generated for the Normal-coordinates, when the diagonal potential part 
      ! is required. In all other cases we work with Lmat = Amat with non-diagonal 
      !  potential representaion
      !
      call Lmat_generation1d
      !
      ! Testing the new routine
      !
      ! call fromlocal2cartesian((/trove%req,trove%alphaeq,trove%taueq/),trove%a0)
      !
      ! We calculate the s_vib vector = (d xi_l / d r_Na )
      !
      !call s_vib_s_rot_polynom1d(s_vib,s_rot)
      !
      call s_vib_s_rot_Sorensen(s_vib,s_rot)
      !
    endif
    !-------------------------------------
    ! We come to the g-s fields definition 
    !-------------------------------------
    !
    if (job%verbose>=5) then
      write (out,"(' Kinetic operator fields need ',f9.3,' Mbytes of memory (plus a bit)')") &
             real(rk*Tcoeff,kind=rk)*real(Nmodes**2+3*Nmodes+9+1,kind=rk)/(1024.0_rk**2)
    end if
    !
    allocate (trove%g_vib(Nmodes,Nmodes),trove%g_rot(3,3),trove%g_cor(Nmodes,3),trove%pseudo,stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate g-fields')") alloc
        stop 'FLinitilize_Kinetic, g-fields - out of memory'
    end if
    !
    ! Vibrational part g_vib : Nmodes x Nmodes matrix 
    !
    ! Calculate Ncoeff: number of elements in the kinetic fields 
    !
    Tcoeff  = trove%RangeOrder(trove%NKinOrder)
    !
    do k1 = 1,Nmodes
       do k2 = 1,Nmodes
          fl => trove%g_vib(k1,k2)
          call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_vib')
       enddo
    enddo
    !
    ! Rotational part g_rot : 3 x 3 matrix 
    !
    do k1 = 1,3
       do k2 = 1,3
          fl => trove%g_rot(k1,k2)
          call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_rot')
       enddo
    enddo
    !
    ! Coriolis part g_cor : Nmodes x 3 matrix 
    !
    do k1 = 1,Nmodes
       do k2 = 1,3
          fl => trove%g_cor(k1,k2)
          call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_cor')
       enddo
    enddo
    !
    ! Pseudo-potential function field initialization 
    !
    fl => trove%pseudo
    call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'pseudo')
    !
    if (trove%internal_coords=='LOCAL') then
       !
       call compute_kinetic_on_rho_grid(Tcoeff)
       !
    else
       !
       ! With s_vib and s_rot we can calculate the kinetic energy operator expansions, namely:
       ! g_vib, g_rot, and the first part of the pseudopotential pseudo1
       ! 
       call gmat_polynom(s_vib,s_rot,trove%g_vib,trove%g_rot,trove%g_cor,trove%pseudo)
       !
       ! We can destroy some fields
       !
       call ArrayStop('s_vib')
       call ArrayStop('s_rot')
       deallocate(s_vib,s_rot)
       !
    endif
    !
    ! be verbose
    ! 
    if (job%verbose>=1) then
        write(out,"('Molecular parameters:')")
        if (Natoms<40) then
          write(out,"('Masses:',40f20.10)") trove%mass(:)
        endif
        !
        if (Nbonds<40) then
           write(out,"('Equilibrium bond lengths    ',40f18.8)") trove%req(:)
        endif
        !
        if (Nangles<40) then
           write(out,"('Equilibrium interbond angles',40f18.8)") trove%alphaeq(:)*180.0_rk/pi
        endif
        !
        if (Nangles<40) then
           write(out,"('Equilibrium dihedral angles',40f18.8)") trove%taueq(:)*180.0_rk/pi
        endif
        !
        write(out,"(/'a0 matrix (equilibrium cartesian coordinates in the xyz-system):')") 
        do iatom = 1,trove%Natoms
           write(out,"(3f18.8)") trove%a0(iatom,:)
        enddo 
        !
        if (trove%lincoord/=0) write(out,"('Molecular is linear lying along ',i4,' axis')") trove%lincoord

    endif 
    !
    if (job%verbose>=2) then
        !
        if (trove%omega(1)>small_) then
           !
           write(out,"(/'Quadratic force constants parameters  in the geometrically defined coordinates :')")
           !
           write(out,"('     k1      k2          const')")
           !
           write(out,"('Harmonic frequencies :')")
           !
           do k1 = 1,trove%Nmodes
              write(out,"(i8,f18.6)") k1,trove%omega(k1)
           enddo
           !
        endif
        !
        write(out,"(1x,i4,'-points central finite difference formula is used.')") difftype
        !
        hstep = epsilon(1.0_ark)**(1.0_rk/(trove%NPotOrder+difftype))
        diferror = hstep**(difftype-1)
        write(out,"(' Optimal spacing / actual spacing : ',2d18.8)") sum(trove%fdstep)/trove%Nmodes,hstep
        write(out,"(' Estimation for the finite differences computational error: ',d18.8)") diferror
        !
    endif
    !
    ! check the smoothness of the follwing field and fix if necessary 
    !
    !
    do k1 = 1,3
       do k2 = 1,3
          fl => trove%g_rot(k1,k2)
          write(txt1,"(i2)") k1
          write(txt2,"(i2)") k2
          call check_field_smoothness(fl,'CHECK',npoints,'g_rot'//trim(adjustl(txt1))//'_'//trim(adjustl(txt2)))
          if (FL_iron_field_out) &
               call check_field_smoothness(fl,'FIX',npoints,'g_rot'//trim(adjustl(txt1))//'_'//trim(adjustl(txt2)))
          !
       enddo
    enddo
    !
    do k1 = 1,Nmodes
       do k2 = 1,3
          fl => trove%g_cor(k1,k2)
          write(txt1,"(i2)") k1
          write(txt2,"(i2)") k2
          call check_field_smoothness(fl,'CHECK',npoints,'g_cor'//trim(adjustl(txt1))//'_'//trim(adjustl(txt2)))
          if (FL_iron_field_out) & 
              call check_field_smoothness(fl,'FIX',npoints,'g_cor'//trim(adjustl(txt1))//'_'//trim(adjustl(txt2)))
       enddo
    enddo
    !
    do k1 = 1,Nmodes
       do k2 = 1,Nmodes
          fl => trove%g_vib(k1,k2)
          write(txt1,"(i2)") k1
          write(txt2,"(i2)") k2
          call check_field_smoothness(fl,'CHECK',npoints,'g_vib'//trim(adjustl(txt1))//'_'//trim(adjustl(txt2)))
          if (FL_iron_field_out) &
             call check_field_smoothness(fl,'FIX',npoints,'g_vib'//trim(adjustl(txt1))//'_'//trim(adjustl(txt2)))
       enddo
    enddo
    !
    fl => trove%pseudo
    call check_field_smoothness(fl,'CHECK',npoints,'pseudo')
    if (FL_iron_field_out) call check_field_smoothness(fl,'FIX',npoints,'pseudo')
    !
    if (FLl2_coeffs) then
      !
      do k1 = 1,Nmodes
         do k2 = 1,Nmodes
            fl => trove%L2_vib(k1,k2)
            write(txt1,"(i2)") k1
            write(txt2,"(i2)") k2
            call check_field_smoothness(fl,'CHECK',npoints,'L2_vib'//trim(adjustl(txt1))//'_'//trim(adjustl(txt2)))
            if (FL_iron_field_out) & 
               call check_field_smoothness(fl,'FIX',npoints,'L2_vib'//trim(adjustl(txt1))//'_'//trim(adjustl(txt2)))
         enddo
      enddo
      !
    endif
    !
    call print_kinetic
    !
    if (trove%sparse) call compact_sparse_kinetic

    if (trim(trove%IO_kinetic)=='SAVE'.and.&
        trove%separate_store.and..not.trove%separate_convert) then
          !
          ! save kinetic.chk as ASCII
          !
          call FLcheck_point_Hamiltonian('KINETIC_SAVE_SPARSE')
          !      
    endif 
    !
    if (job%verbose>=4) call MemoryReport
    !
    call TimerStop('Kinetic')
    !
    call TimerReport
    !
    if (job%verbose>=4) write(out,"('FLinitilize_Kinetic/end')")   
    !
  end subroutine FLinitilize_Kinetic


    !
    subroutine print_kinetic
       !
       implicit none
       !
       integer(ik) :: k1,k2,i,irho,iatom,imode,Nmodes,Npoints
       !
       Nmodes = trove%Nmodes
       Npoints = trove%Npoints
       !
       if (job%verbose>=6.or.(job%verbose>=6.and.manifold==0)) then
           write(out,"(/'Kinetic parameteres    irho  k1    k2   i    g_vib             g_cor            g_rot:')")
           !
           do k1 = 1,Nmodes
              do k2 = 1,Nmodes
                 do i = 1,trove%g_vib(k1,k2)%Ncoeff
                    do irho = 0,Npoints
                      if (k2<=3.and.k1<=3) then 
                        write(out,"(20x,3i5,1x,i7,3e18.8)") irho,k1,k2,i,trove%g_vib(k1,k2)%field(i,irho), &
                                                      trove%g_cor(k1,k2)%field(i,irho),trove%g_rot(k1,k2)%field(i,irho)
                      elseif (k2<=3) then
                        write(out,"(20x,3i5,1x,i7,2e18.8)") irho,k1,k2,i,trove%g_vib(k1,k2)%field(i,irho),&
                                   trove%g_cor(k1,k2)%field(i,irho)
                      else 
                        write(out,"(20x,3i5,1x,i7,e18.8)") irho,k1,k2,i,trove%g_vib(k1,k2)%field(i,irho)
                      endif
                    enddo
                 enddo
              enddo
           enddo
           !
           if (job%verbose>6) then 
             !
             write(out,"(/'b0 and db0 matrices:')") 
             do iatom = 1,trove%Natoms
                   !
                   if (manifold==0) then
                     write(out,"(i5,3f18.8)") iatom,trove%b0(iatom,:,0)
                   else
                     !
                     do irho = 0,Npoints
                       !
                       write(out,"(2i5,12f18.8)") iatom,irho,trove%b0(iatom,:,irho),trove%db0(iatom,:,irho,1),&
                             trove%db0(iatom,:,irho,2),trove%db0(iatom,:,irho,3)
                       !
                     enddo
                   endif 
             enddo
             !
             write(out,"(/'Amatrho matrices:')") 
             do iatom = 1,trove%Natoms
                do imode = 1,Nmodes
                   do irho = 0,Npoints
                     !
                     write(out,"(3i5,3f18.8)") iatom,imode,irho,trove%Amatrho(iatom,:,imode,irho)
                     !
                   enddo
               enddo
             enddo
             !
             write(out,"(/'dAmatrho matrices:')") 
             do iatom = 1,trove%Natoms
                do imode = 1,Nmodes
                   do irho = 0,Npoints
                     !
                     write(out,"(i5,9f18.8)") irho,trove%dAmatrho(iatom,:,imode,irho,1),trove%dAmatrho(iatom,:,imode,irho,2),&
                           trove%dAmatrho(iatom,:,imode,irho,3)
                     !
                   enddo
               enddo
             enddo
             !
             write(out,"(/'Bmatrho matrices:')") 
             do imode = 1,Nmodes
               do iatom = 1,trove%Natoms
                   do irho = 0,Npoints
                     !
                     write(out,"(3i5,3f18.8)") imode,iatom,irho,trove%Bmatrho(imode,iatom,:,irho)
                     !
                   enddo
               enddo
             enddo
             !
             if (manifold/=0) then
               write(out,"(/'dBmatrho matrices:')") 
               do iatom = 1,trove%Natoms
                  do imode = 1,Nmodes
                     do irho = 0,Npoints
                       !
                       write(out,"(i5,6f18.8)") irho,trove%dBmatrho(imode,iatom,:,irho,1),trove%dBmatrho(imode,iatom,:,irho,2)
                       !
                     enddo
                 enddo
               enddo
               !
             endif
             !
           endif
           !
           if (trove%Coordinates(1,1)=='NORMAL') then
              !
              write(out,"('Normal quadratic pot. parameteres   qwforce:')")
              !
              do k1 = 1,Nmodes
                 do k2 = 1,k1
                    do irho = 0,Npoints
                       write(out,"(20x,3i5,d18.8)") irho,k1,k2,trove%qwforce(k1,k2,irho)
                    enddo
                 enddo
              enddo
              !
           endif
           !
       endif ! job%verbose 
       !
    end subroutine print_kinetic


    subroutine compact_sparse_kinetic
       !
       implicit none
       !
       integer(ik) :: k1,k2,i,irho,iatom,imode,Nmodes,Npoints
       type(FLpolynomT),pointer     :: fl,gl
       !
       if (job%verbose>=4) write(out,"('Compacting the kinetic energy matrices into a sparse representation ...')")
       !
       Nmodes = trove%Nmodes
       Npoints = trove%Npoints
       !
       do k1 = 1,Nmodes
          do k2 = 1,Nmodes
             !
             if (k1==Nmodes.and.k2==Nmodes) cycle
             !
             fl => trove%g_vib(k1,k2)
             !
             call FLCompact_a_field_sparse(fl,"g_vib")
             !
          enddo
       enddo
       !
       do k1 = 1,Nmodes
          do k2 = 1,3
             !
             fl => trove%g_cor(k1,k2)
             call FLCompact_a_field_sparse(fl,"g_cor")
             !
          enddo
       enddo
       !
       do k1 = 1,3
          do k2 = 1,3
             !
             fl => trove%g_rot(k1,k2)
             !
             if (trove%triatom_sing_resolve.and.(k1==3.and.k2==3)) cycle
             !
             call FLCompact_a_field_sparse(fl,"g_rot")
             !
          enddo
       enddo
       !
       fl => trove%g_vib(Nmodes,Nmodes)
       !
       if (.not.trove%triatom_sing_resolve) then 
         call FLCompact_and_combine_fields_sparse(fl,"g_vib",trove%pseudo,"pseudo")
         !
       else
         !
         gl => trove%g_rot(3,3)
         !
         call FLCompact_and_combine_three_fields_sparse(fl,"g_vib",gl,"g_rot",trove%pseudo,"pseudo")
         !
       endif
       !
       if (FLl2_coeffs) then
         !
         do k1 = 1,Nmodes
            do k2 = 1,Nmodes
             fl => trove%L2_vib(k1,k2)
             call FLCompact_a_field_sparse(fl,"L2_vib")
            enddo
         enddo
         !
       endif
       !
       if (job%verbose>=4) write(out,"('... done!')")
       !
    end subroutine compact_sparse_kinetic





  !
  ! This procedure is to check the smoothness of a field wrt to rho and 
  ! fix it if necessary 
  !
  subroutine check_field_smoothness(object,action,npoints,msg)
    !
    !real(ark),intent(in) :: field(:)
    type(FLpolynomT),pointer :: object
    !
    integer(ik),intent(in) :: npoints
    character(len=*),intent(in) :: action,msg
    integer(ik)     :: irho,N_
    integer(ik),parameter :: N_max = 10,Nattempt_max = 1000
    integer(ik)                  :: i,i1,i2,k,ioutlier
    logical                      :: outliers
    real(ark)                    :: rho_(1:N_max),func_(1:N_max),rho,fval,foutlier,fcorr,df
    integer(ik)                  :: iattempt,Nattempt,ioutlier_max,Ncoeff,iterm,info
    integer(ik),allocatable      :: outlier(:)

      if (job%verbose>=6) write(out,"('  check field smoothness for ',a)") msg
      !
      if (manifold==0) return
      !
      select case (action)
      !
      case ('CHECK')
       !
       Nattempt = 1
       !
      case ('FIX')
       !
       Nattempt = Nattempt_max
       !
      end select 
      !
      N_ = min(N_max,npoints)
      allocate (outlier(0:npoints),stat=info)
      call ArrayStart('check_field_smoothness',info,size(outlier),kind(outlier))
      outlier = 0
      !
      Ncoeff = object%Ncoeff
      !
      !$omp parallel do private(iterm,iattempt,outliers,foutlier,ioutlier_max,ioutlier,irho,i1,i2,k,i,rho_,func_,rho,fval,df,fcorr) shared(outlier) schedule(guided)
      do iterm = 1, Ncoeff
        !
        if (job%verbose>=6) write(out,"('  iterm = ',i5)") iterm
        !
        iattempt = 0 ; outliers = .true.
        !
        do while (outliers.and.iattempt<Nattempt)
          !
          iattempt = iattempt + 1
          !
          outliers = .false.
          foutlier = 0
          ioutlier_max = -1
          ioutlier = 0
          !
          if (job%verbose>=7) write(out,"('  iattempt = ',i3)") iattempt
          !
          do irho = 0, npoints
            !
            if (trove%periodic) then 
              !
              i1 = irho-N_/2
              i2 = irho+N_/2
              !
              k = 0
              do i = i1,i2
                 if (i==irho.or.k>N_) cycle
                 k = k + 1
                 rho_(k) = trove%rho_border(1)+i*trove%rhostep
                 !
                 if (i<0) then 
                   func_(k) = object%field(iterm,npoints+i)
                 elseif(i>npoints) then
                   func_(k) = object%field(iterm,i-npoints)
                 else
                   func_(k) = object%field(iterm,i)
                 endif
              enddo
              !
            else
              !
              i1 = max(irho-N_/2,0)
              i2 = min(irho+N_/2,npoints)
              !
              k = 0
              do i = i1,i2
                 if (i==irho.or.k>N_) cycle
                 !if (outlier(i)==1) cycle
                 k = k + 1
                 rho_(k) = trove%rho_border(1)+i*trove%rhostep
                 func_(k) = object%field(iterm,i)
              enddo
              !
            endif
            !
            !if ( k<1.or.k>N_max ) then 
            !   write(out,"('check_field_smoothness is out of range, k = ',i5,' <> [ 0,',i2,']')") k,N_max
            !   stop 'check_field_smoothness is out of range'
            !endif
            !
            !k = min(max(1,k),N_max)
            !
            rho = trove%rho_border(1)+irho*trove%rhostep
            !
            call polintark(rho_(1:k),func_(1:k),rho,fval,df)
            !
            if ( abs(fval-object%field(iterm,irho))>max(abs(object%field(iterm,irho))*1e4*sqrt(small_),1e3*sqrt(small_)) ) then
               !
               outliers = .true.
               outlier(irho) = 1
               !
               ioutlier = ioutlier + 1
               !
               if ( abs(fval-object%field(iterm,irho))>foutlier ) then 
                  ioutlier_max = irho
                  foutlier = abs(fval-object%field(iterm,irho))
                  fcorr = fval
               endif
               !
            endif
            !
          end do
          !
          if (ioutlier_max>=0) then 
            !
            if (job%verbose>=5) write(out,"('check_field_smoothness: ',a)",advance='NO') msg
            if (job%verbose>=5) write(out,"('; an outlier found for iterm =  ',i6,' at i = ',i5,': ',e18.11,' vs ',e18.11)") &
                                      iterm,ioutlier_max,object%field(iterm,ioutlier_max),fcorr
            !
            if (Nattempt>1) then 
               if (job%verbose>=5) write(out,"(e18.11,' will be replaced by the extrapolated value = ',e18.11)") &
                                   object%field(iterm,ioutlier_max),fcorr
               object%field(iterm,ioutlier_max) = fcorr
            endif
            !
          endif 
          !
        enddo
        !
        if ( iattempt==Nattempt.and.outliers.and.Nattempt>1) then
          write(out,"('check_field_smoothness: ',a)") msg
          write(out,"('      Too many outliers, it was impossible to fix after ',i9,' attempts for iterm = ',i5 )") Nattempt,iterm
          !stop 'check_field_smoothness: too many outliers'
        endif
        !
      enddo
      !$omp end parallel do
      !
      deallocate(outlier)
      call ArrayStop('check_field_smoothness')
      !
 end subroutine check_field_smoothness



  !
  ! This procedure is to compute kinetic fields on the rho-grid
  !
  subroutine compute_kinetic_on_rho_grid(Nterms)
    !
    !
    integer(ik),intent(in) :: Nterms
    real(ark),allocatable :: g_vib(:,:,:),g_rot(:,:,:),g_cor(:,:,:),pseudo(:)
    integer(ik)  :: k1,k2,irho,npoints,info,Nmodes
    real(ark)    :: rho,factor
      !
      Nmodes = trove%Nmodes
      !
      ! Conversion factor to the cm-1 units 
      !
      factor = real(planck,ark)*real(avogno,ark)*real(1.0d+16,kind=ark)/(4.0_ark*pi*pi*real(vellgt,ark))
      !
      allocate (g_vib(Nmodes,Nmodes,Nterms),stat=info)
      call ArrayStart('kinetic_on_grid-fields',info,size(g_vib),kind(g_vib))
      allocate (g_cor(Nmodes,Nmodes,Nterms),stat=info)
      call ArrayStart('kinetic_on_grid-fields',info,size(g_cor),kind(g_cor))
      allocate (g_rot(Nmodes,Nmodes,Nterms),stat=info)
      call ArrayStart('kinetic_on_grid-fields',info,size(g_rot),kind(g_rot))
      allocate (pseudo(Nterms),stat=info)
      call ArrayStart('kinetic_on_grid-fields',info,size(pseudo),kind(pseudo))
      !
      npoints = trove%npoints 
      !
      do irho = 0, npoints
         !
         rho = trove%rho_i(irho)
         !
         call MLkineticfunc(trove%nmodes,Nterms,rho,g_vib,g_rot,g_cor,pseudo)
         !
         do k1 = 1,Nmodes
            do k2 = 1,Nmodes
               trove%g_vib(k1,k2)%field(:,irho) = g_vib(k1,k2,:)*factor
            enddo
         enddo
         !
         do k1 = 1,Nmodes
            do k2 = 1,3
               trove%g_cor(k1,k2)%field(:,irho) = g_cor(k1,k2,:)*factor
            enddo
         enddo
         !
         do k1 = 1,3
            do k2 = 1,3
               trove%g_rot(k1,k2)%field(:,irho) = g_rot(k1,k2,:)*factor
            enddo
         enddo
         !
         trove%pseudo%field(:,irho) = pseudo(:)*factor
        !
      enddo
      !
      !
      deallocate(g_vib,g_rot,g_cor,pseudo)
      call ArrayStop('kinetic_on_grid-fields')
      !
 end subroutine compute_kinetic_on_rho_grid



! The procedire where we calculate Amat akin l_nal paramiters
! using general recurcive procedure:
! we solve the Eckart equations numerically
! There are M x M equations for M variables Amat (A_Nal parameters), 
! where M=3*Natoms*Nmodes  (N=1..Natoms; a = x,y,z, l = 1..Nmodes)
! (3 Eckart equation for the mass center )         x Nmodes
! (3  Eckart equations for the intertia mass axes) x Nmodes
! Nmodes*(Nmodesd+1)/2 equations for the orthonormality  conditions 
! Nmodes*(Nmodes-1)/2 equations for the orthogonality of the quadratic potential part
! The rest equations we get from expanding the local coordinates r_l 
! Equations are
! 1. sum_{N,a} m_N Amat(N,a,l) = 0 
! 2. sum_{N,a,b} m_N eps_{a,b,g} a_{N,b} Amat(N,g,l) = 0 
! 3. sum_{N,a} lmat_{l,N,a} Amat{N,a,m} = delta_{l,m}
! 4. f(i,j) = 0 for i/=j
! where matrix B is inverse transformtaion to the matrix A
!
! The routine also takes into account the explicitly defined 
! rank=1 variable in the numerical representation
!
  subroutine Lmat_generation1d

    real(ark)   :: a0(trove%Natoms,3)

    real(ark)   :: a_t,b_t(3)
    !
    integer(ik) ::  Nmodes,Natoms,Nbonds,Nangles
    integer(ik) ::  iatom,imode
    integer(ik) ::  jatom,ix,jx,kx,jmode,kmode,lfolds(trove%Nmodes),foldness(trove%Nmodes,trove%Nmodes)
    integer(ik) ::  ieq,ivar
    integer(ik) ::  alloc,dimen  , lmode, ifolds(trove%Nmodes),ifold
    integer(ik) ::  Nequat,lincoord,Npoints,ierror
    !
    real(ark),allocatable  ::  Bmat(:,:,:),Bmat_t(:,:,:),Amat(:,:,:)
    real(ark),allocatable  ::  bm(:) 
    real(ark),allocatable  ::  Tmat(:,:) 

    real(ark),allocatable  ::  f(:,:),c(:,:),b(:,:),a(:,:)
    !
    real(rk),allocatable   ::  db(:,:),da(:,:)
    !
    real(ark)                    :: astep(2)
    !
    real(ark)   :: factor,aJ_over_hc,rhostep
    real(ark)   :: q_eq(trove%Ncoords)
    real(ark)   :: chi_eq(trove%Nmodes)
    real(ark)   :: coordtransform(trove%Nmodes,trove%Ncoords)
    real(ark)   :: step(2,trove%Ncoords)
    !
    real(ark)   :: xna_step(trove%Natoms*3,2),xna(trove%Natoms*3)
    integer(ik) :: nindex(trove%Natoms*3)
    !
    character(len=cl)  :: job_is
    logical ::  warning_b0,warning_amat,dir
    !
    integer(ik) :: kindex(trove%Ncoords),irho
    !
    if (job%verbose>=2) write(out,"(/'Lmat_generation1d/start  ')") 

    ! for easy reference 
    !
    Nmodes = trove%Nmodes
    Natoms = trove%Natoms
    Nbonds = trove%Nbonds
    Nangles = trove%Nangles
    Npoints = trove%Npoints
    !
    lincoord = 0
    if ( Nmodes==3*Nmodes-5 ) lincoord = trove%lincoord
    !
    !lincoord = 0
    !
    ! these will check if there are some minor problems with definition of b0 and amat
    ! e.g. not all Eckarts are satisfied 
    !
    warning_b0   = .false.
    warning_amat = .false.
    !
    factor = planck*avogno*real(1.0d+16,kind=rk)/(4.0_ark*pi*pi*vellgt)
    aJ_over_hc = 1.0e-11/planck/vellgt

    !
    ! define number of equation that will be solved 
    !
    allocate (f(Nmodes,Nmodes),stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i9,' trying to allocate array f')") alloc
       stop 'Lmat_generation1d, f  - out of memory'
    end if
    !
    Nequat = 6+Nmodes-min(1,lincoord)
    !
    allocate (Amat(Nmodes,3,Nmodes), &
              Bmat(Nmodes,Natoms,3),Bmat_t(trove%Ncoords,Natoms,3),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate array for Bmat')") alloc
        stop 'Lmat_generation1d, Bmat - out of memory'
    end if

    ! 
    ! we will also neet a Tmat matryix and bm vector for T x =b linear equation 
    ! Tmat is a M by M matrix, where M is a number of equations and variables. 
    ! Since all equations can  be solved for every mode independently, 
    ! M = 3 + 3 + Nmodes-max(1,licoord)
    ! 
    allocate (Tmat(Nequat,Nequat),Bm(Nequat),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate array for Tmat and Bm')") alloc
        stop 'Lmat_generation1d, Bm ant Tmat - out of memory'
    end if

    !
    ! The linear equation will be solve with Lapack at double precision 
    ! therefor we will need a double precision matrix A and vector b 
    ! as analogs of Tmat and Bm, respectively 
    !
    allocate (a(Nequat,Nequat),b(Nequat,1),c(Nmodes,Nmodes),da(Nequat,Nequat),db(Nequat,1),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate array for a and b')") alloc
        stop 'Lmat_generation1d, a and b  - out of memory'
    end if
    !
    ! It will be needed to transform numerically coordinates:
    !
    q_eq(:) = trove%local_eq(:)
    !
    rhostep = trove%rhostep
    !
    ! it is gonna be a big loop over all Npoints+1
    !
    rho_loop : do irho = 0,Npoints
      !
      if (job%verbose>=6) write(out,"('irho = ',i7)") irho
      !
      Amat = 0
      !
      a0(1:Natoms,1:3) = trove%b0(1:Natoms,1:3,irho)
      !
      ! Before we start we can check a0(N,a) if it was introduced correctly by the user 
      ! and it sutisfies the Eckart conditions 
      !
      ! First Eckart equation
      ! 
      do ix = 1,3 
          !
          a_t = sum(trove%mass(:)*a0(:,ix))
          if (abs(a_t)>10.0_ark**(-rk)) then 
               write(out,"('Lmat_generation1d: b0 is not a solution of Eckart 1 for  ix =',i4,d18.8)") ix,a_t
                  stop 'Lmat_generation1d: b0 is not solution of Eckart1'
          endif
          !
      enddo
      !
      ! Second Eckart equation
      ! 
      do ix = 1,3 
         do jx = 1,3 
            if (ix/=jx) then  
               a_t =  sum(trove%mass(:)*a0(:,ix)*a0(:,jx) )

               if (abs(a_t)>100.0_rk*small_) then 
                   if (job%verbose>=5) write(out,"('b0 is not good at Eckart 2 for ix,jx =',2i4,d18.8)") ix,jx,a_t
                   !
                   warning_b0 = .true.
                   !
                   !if (trove%internal_coords/='LOCAL') then
                   !   stop 'Lmat_generation1d: b0 is not solution of Eckart2'
                   !endif
               endif
            endif
         enddo
         !
      enddo
      !
      ! Defining Bmat matrix
      !
      call Bmat_generation(a0,Bmat_t)
      !
      ! Find linear coordinate transformation from r to chi
      !
      chi_eq(:) = trove%chi_ref(:,irho)
      !
      !chi_eq(Nmodes) = trove%rho_i(irho)
      !
      dir = .false.
      !
      q_eq = MLcoordinate_transform_func(chi_eq,size(q_eq),dir)
      !  
      step = trove%fdstep(1)
      !
      !do imode=1,trove%Nmodes
      !   !
      !   call ML_check_steps4coordinvert(q_eq,2,imode,astep)
      !   !
      !   step(:,imode) = trove%fdstep(imode)*astep(:)
      !   !
      !enddo
      !
      !job_is = 'cartesian2local'
      !
      !ieq = 0 
      !
      !do iatom = 1,Natoms
      !   do ix = 1,3
      !      !
      !      ieq = ieq + 1
      !      !
      !      xna(ieq) = trove%b0(iatom,ix,irho)
      !      !
      !   enddo
      !enddo
      !
      !xna_step = trove%fdstep(1)
      !
      !ieq = 0
      !
      !do iatom = 1,Natoms
      !   do ix = 1,3
      !      !
      !      ieq = ieq + 1
      !      !
      !      nindex = 0 ; nindex(ieq) = 1
      !      !
      !      Bmat_t(:,iatom,ix) = FLvect_finitediffs(job_is,trove%Ncoords,nindex,xna,xna_step,irho)
      !      !
      !   enddo
      !enddo
      !
      ! for the derivatives d r_na / d xi we will need the vector function 'MLcoordinate_transform_func'
      ! which will be called by envoiking the following option 
      !
      job_is = 'local2chi'
      !
      do imode =1,trove%Ncoords
         !
         kindex = 0 ; kindex(imode) = 1
         !
         coordtransform(:,imode) = FLvect_finitediffs(job_is,trove%Ncoords,kindex,q_eq,step,irho)
         !
         !chi2 = MLcoordinate_transform_func(q2,size(chi2),dir)
         !
      enddo
      !
      do iatom = 1,Natoms
         do ix = 1,3
            Bmat(:,iatom,ix) = matmul(coordtransform,Bmat_t(:,iatom,ix))
         enddo
      enddo
      !
      do imode = 1,Nmodes
         !
         if (imode==trove%Nmodes_n) then 
            !
            trove%Bmatrho(imode,:,:,irho) = trove%db0(:,:,irho,1)
            !
         else
            !
            trove%Bmatrho(imode,:,:,irho) = Bmat(imode,:,:)
            !
         endif
         !
      enddo
      !
      ! Quadratic force constants in local coordinates are stored in "f"
      !
      !f(1:Nmodes,1:Nmodes) = trove%qwforce(1:Nmodes,1:Nmodes,irho)
      !
      ! Constructing a system of linear equations T x = b, where x = Amat
      ! The equations can be solved for every mode indipendently,  
      ! i.e. Tmat doesn't depend on the internal coordinates 
      !
      Amat = 0 
      Tmat = 0 
      ieq  = 0 
      !
      ! First Eckart equation
      !
      do ix = 1,3
         !
         !if (ix==lincoord) cycle
         !
         ieq = ieq+1
         !
         ivar = 0 
         do jatom = 1,Natoms  
            do jx = 1,3
               !
               !if (jx==lincoord) cycle
               ! 
               ivar = ivar+1
               if (jx==ix) then 
                  Tmat(ieq,ivar) = sqrt(trove%mass(jatom))
               endif
             enddo
         enddo
         !
      enddo
      !
      ! Second Eckart equation
      !
      do ix = 1,3
         ! 
         if (ix/=lincoord) then 
            !
            ieq = ieq+1
            !
            ivar = 0 
            do jatom = 1,Natoms  
               do jx = 1,3
                  !
                  !if (jx==lincoord) cycle 
                  !
                  ivar = ivar+1
                  Tmat(ieq,ivar) = sqrt(trove%mass(jatom))*sum( epsil(ix,:,jx)*a0(jatom,:) )
               enddo
            enddo
            !
         endif
         !
      enddo
      !
      ! Equation for the coordinates  B*A = delta
      !
      do jmode = 1,Nmodes
         !
         ieq = ieq+1
         !
         ivar = 0 
         do jatom = 1,Natoms  
            !
            do jx = 1,3
               !
               !if (jx==lincoord) cycle
               !
               ivar = ivar+1
               !
               !   or Sayvetz equation
               if (jmode==trove%Nmodes_n) then 
                  !
                  Tmat(ieq,ivar) = trove%db0(jatom,jx,irho,1)*sqrt(trove%mass(jatom))
                  !
               else
                  !
                  Tmat(ieq,ivar) = Bmat(jmode,jatom,jx)/sqrt(trove%mass(jatom))
                  !
               endif 
               !
            enddo
            !
         enddo
         !
      enddo
      !
      !
      ! Right side of the equation is different for different modes 
      !
      dimen = ieq
      !
      if (dimen/=3*Natoms.or.dimen/=6+Nmodes-min(lincoord,1).or.dimen/=size(a,dim=1)) then 
          write(out,"('Lmat_generation1d: the size of Tmat contradicts the number of eqs or vars:',3i7)") 6+Nmodes,ieq,ivar
          !stop 'Lmat_generation1d: ieq is inconsistent'
      endif 
      !
      Amat = 0
      !
      do imode = 1,Nmodes 
         !
         ! The only place we the right side of the equation appears is in the orthogonalization equations:
         ! bm(i) = delta(i,6+imode)
         !
         bm = 0.0_ark  ;  bm(Nequat-Nmodes+imode) = 1.0_ark
         !
         a = Tmat
         b(:,1) = bm(:)
         !
         call MLlinurark(Nequat,a,b(:,1),bm,ierror)
         !
         if (ierror/=0) then
           !
           da = a
           db = b
           !
           call lapack_gelss(da(:,:),db(:,:))
           !
           bm(:) = db(:,1)
           !
         end if
         !
         ! The solution written in bm -> Amat
         !
         ivar= 0
         !
         do jatom = 1,Natoms  
            do jx = 1,3
               !if (jx==lincoord) cycle 
               ivar = ivar+1
               Amat(jatom,jx,imode) = bm(ivar)
             enddo
         enddo
         !
      enddo 
      !
      ! Internal test for Amat: 
      ! check if the Eckart equations are satisfied
      !
      ! First Eckart equation
      ! 
      do imode = 1,Nmodes
          do ix = 1,3
             !if (ix==lincoord) cycle 
             a_t = sum(sqrt(trove%mass(1:Natoms))*Amat(1:Natoms,ix,imode))
             if (abs(a_t)>1000.0_rk*sqrt(small_)) then 
                 write(out,"('Lmat_generation1d: Eckart 1 is not =0 for irho = ',i6,',imode = ',i4,', ix =',i4,d18.8)") & 
                                                        irho,imode,ix,a_t
                 stop 'Lmat_generation1d: Eckart 1 is not solved'
             endif
          enddo
          !
      enddo
      !
      ! Second Eckart equation
      ! 
      do imode = 1,Nmodes
          do ix = 1,3
             !
             if (trove%lincoord==ix) cycle
             ! 
             a_t = 0.0_ark
             do jx = 1,3 
                do kx = 1,3 
                   a_t = a_t + sum(sqrt(trove%mass(1:Natoms))*a0(1:Natoms,jx)*epsil(ix,jx,kx)*Amat(1:Natoms,kx,imode) )
                enddo
             enddo
             if (abs(a_t)>1000.0_rk*sqrt(small_)) then 
                 write(out,"('Lmat_generation1d:')") 
                 write(out,"('Eckart 2 is not =0 for irho = ',i4,', imode = ',i4,', ix =',i4,d18.8)") irho,imode,ix,a_t
                 warning_amat = .true.
                 if (trove%internal_coords/='LOCAL') then
                    stop 'Lmat_generation1d: Eckart 2 is not solved'
                 endif
             endif
          enddo
          !
      enddo
      !
      ! Third equation: orthogonality to Bmat or Sayvetz
      ! 
      do imode = 1,Nmodes
         do jmode = 1,Nmodes
            !
            a_t = 0.0_ark
            !
            if (jmode==imode) a_t = -1.0_ark
            !
            do ix = 1,3
               !
               !if (trove%lincoord==ix) cycle
               !
               !   or Sayvetz equation
               !
               if (jmode==trove%Nmodes_n) then 
                  !
                  a_t = a_t + sum(trove%db0(:,ix,irho,1)*sqrt(trove%mass(1:Natoms))*Amat(1:Natoms,ix,imode) )
                  !
               else
                  !
                  a_t = a_t + sum(Bmat(jmode,1:Natoms,ix)/sqrt(trove%mass(1:Natoms))*Amat(1:Natoms,ix,imode) )
                  !
               endif 
               !
            enddo
            !
            if (abs(a_t)>100.0_rk*sqrt(small_)) then 
               !
               ! we do not worry much if this happnens when 
               ! all derivatives Bmat(jmode,:,ix)=0  
               !
               if (jmode==imode.and.abs(a_t+1.0_rk)<sqrt(small_).and.sum(Bmat(jmode,:,:)**2)<sqrt(small_)) then 
                 !
                 write(out,"('Lmat_generation1d: we are aware that all Bmat for irho=',i4,',jmode = ',i4,' are zero')") irho,jmode
                 write(out,"('               and orthogon. eq. 3 is not = 0 at imode = ',i4,', jmode =',i4,d18.8)") imode,jmode,a_t
                 !
                 !trove%sing_at_rho_0 = .true.
                 !
               else
                 !
                 write(out,"('Lmat_generation1d: Orthogon. eq. 3 is not = 0, irho = ',i4,',imode = ',i4,', jmode =',i4,d18.8)") &
                              irho,imode,jmode,a_t
                 if (trove%internal_coords/='LOCAL') then
                     stop 'Lmat_generation1d: Orthon.eq. 3 is not solved'
                 endif 
                 !
               endif 
            endif
            !
         enddo
         !
      enddo
      !
      if (trove%Coordinates(1,1)=='NORMAL') then 
         !
         ! Store the result 
         !
         trove%Amatrho(1:Natoms,1:3,1:Nmodes,irho) = Amat(1:Natoms,1:3,1:Nmodes)
         !
         ! Quadratic potential part - must be also diagonal for the normal-coordinates-case
         !
         call Lmat_qwforce
         !
      else
         !
         ! In general, we just get rid of 1/sqrt(mass) in the definition of Lmat:
         ! by multiplying the Amat by 1/sqrt(omega)
         !
         ! factor = planck*avogno*real(1.0d+16,kind=rk)/(4.0_rk*pi*pi*vellgt)
         !
         do imode =1,Nmodes
            do iatom =1,Natoms
               Amat(iatom,:,imode) = Amat(iatom,:,imode)/sqrt(trove%mass(iatom))
            enddo
         enddo
         !
      endif 
      !
      ! Store the result 
      !
      trove%Amatrho(1:Natoms,1:3,1:Nmodes,irho) = Amat(1:Natoms,1:3,1:Nmodes)
      !
      !if (job%verbose>=6) then 
      !  write(out,"(i8,18(<Nmodes>f16.8))") irho,(Amat(1:Natoms,ix,1:Nmodes),ix=1,3)
      !endif
      !
    enddo rho_loop
    !
    ! do some reporting 
    !
    if (job%verbose>=6) then
      write(out,"(a)") 'irho, iatom, imode, Amat:'
      do iatom =1,Natoms
        do imode =1,Nmodes
           do irho = 0,Npoints
              write(out,"(i8,1x,i3,1x,i3,1x,3(f16.8))") irho,iatom,imode,(trove%Amatrho(iatom,ix,imode,irho),ix=1,3)
           end do
         end do
      end do
    endif
    !
    ! Generate the matrix with derivatives of Amatrho wrt to rho 
    ! only when it is an 1D case 
    !
    trove%dAmatrho = 0
    !
    if (manifold/=0) then 
       !
       do iatom = 1,Natoms
         !
         do ix = 1,3
           !
           do imode = 1,Nmodes
              !
              call diff_2d_4points_ark( Npoints,trove%rho_border,trove%Amatrho(iatom,ix,imode,0:Npoints  ),&
                                   job%bset(imode)%periodic,0_ik,trove%dAmatrho(iatom,ix,imode,0:Npoints,1),&
                                                            trove%dAmatrho(iatom,ix,imode,0:Npoints,2))
              call diff_3d_6points( Npoints,trove%rho_border,trove%Amatrho(iatom,ix,imode,0:Npoints  ),&
                                   job%bset(imode)%periodic,trove%dAmatrho(iatom,ix,imode,0:Npoints,3))
                                       !
              call diff_2d_4points_ark(Npoints,trove%rho_border,trove%Bmatrho(imode,iatom,ix,0:Npoints  ),&
                                  job%bset(imode)%periodic,0_ik,trove%dBmatrho(imode,iatom,ix,0:Npoints,1),&
                                                           trove%dBmatrho(imode,iatom,ix,0:Npoints,2))
              !
           enddo
           !
         enddo
         !
       enddo
       !
    endif 
    !
    if ( warning_amat ) then
       write(out,"('Warning-Lmat_generation1d: some Eckarts were not fullfilled for Amat')")
       write(out,"('                           turn on verbose=5 to see the details')")
    endif 
    !
    if ( warning_b0 ) then
       write(out,"('Warning-Lmat_generation1d: some Eckarts were not fullfilled for b0')")
       write(out,"('                           turn on verbose=5 to see the details')")
    endif 
    !
    ! cleaning up 
    !
    deallocate(a,b,c,da,db)
    !
    deallocate(Tmat,bm,Bmat,Bmat_t,f)

    if (job%verbose>=2) write(out,"('Lmat_generation1d/end  ')") 
    !
    !
    ! Lmat has to provide the quadratic part potential field to be diagonal
    ! in case of the normal coordinate representaion 
    ! We diagonalize it by rotating Lmat 
    !
    !
   contains 
    !
    subroutine Lmat_qwforce
      !
      real(rk),allocatable :: f_t(:,:),a(:,:),b(:,:),c(:,:)
      real(ark)   :: omega(trove%Nmodes),step(2,trove%Nmodes),xi_eq(trove%Nmodes),df
      integer(ik) :: imode,k1,k2
      !
      step(1,:) = trove%fdstep(:)
      step(2,:) = trove%fdstep(:)

      !
      do imode=1,trove%Nmodes_e
         !
         xi_eq(imode) = MLcoord_direct(trove%chi_ref(imode,0),2,imode)
         !
      enddo
      !
      do k1=1,trove%Nmodes
         do k2=1,k1
            !
            kindex = 0 ; kindex(k1) = kindex(k1)+1 ; kindex(k2) = kindex(k2)+1
            !
            df = FLfinitediffs(kindex,poten_chi,xi_eq,step)


            !df = FLfinitediffs(kindex,poten_chi,trove%chi_ref(:,0),step) 
            !
            f(k1,k2) = df/aJ_over_hc
            f(k2,k1) = df/aJ_over_hc
            !
         enddo 
      enddo 
      !
      !f = matmul(matmul((trove%coordtransform),f),transpose(trove%coordtransform))
      !
      f = f*0.5_ark
      !
      ! Orthogonalization by diagonalization 
      !
      ! The diagonalization will be done with Lapack at double precision 
      ! therefore we will need a double precision matrix A and vector b 
      !
      allocate (a(Nmodes,Nmodes),b(Nmodes,1),c(Nmodes,Nmodes),f_t(Nmodes,Nmodes),stat=alloc)
      if (alloc/=0) then
          write (out,"(' Error ',i9,' trying to allocate array for a and b')") alloc
          stop 'Lmat_generation, a and b  - out of memory'
      end if
      !
      do imode = 1,Nmodes
         do jmode = 1,Nmodes
            a_t = sum(sum(Amat(:,:,jmode)*Amat(:,:,imode),dim=1))       !
            a(imode,jmode) =  real(a_t,kind=rk)
         enddo
         !
      enddo
      !
      call lapack_syev(a,b(:,1))
      !
      ! Building up the orthogonal transformtaion "c"
      !
      do imode = 1,Nmodes
         do jmode = 1,Nmodes
            c(imode,jmode) = real(a(imode,jmode)/sqrt(b(jmode,1)),kind=rk)
         enddo
      enddo
      !
      ! Transformation to an orthogonal Amat
      !
      do iatom = 1,Natoms
         do ix = 1,3 
            Amat(iatom,ix,:) = matmul(Amat(iatom,ix,:),c)
         enddo
      enddo
      !
      !
      ! Check for the orthogonality sum_{Na} Amat_{Na l} Amat_{Na m} = 0 for m/=l
      !
      do imode = 1,Nmodes
         !
         do jmode = 1,Nmodes
            !
            a_t = sum(sum(Amat(:,:,jmode)*Amat(:,:,imode),dim=1))       !
            !
            if (abs(a_t-1.0_rk)>100.0_rk*small_.and.imode==jmode) then 
               write(out,"('Lmat_generation: Negative or zero diagonal element for Amat x Amat for imode = ',i4,d18.8)") imode,a_t
               stop 'mat_generation: Amat is not Orthogonal'
            endif
            !
            if (abs(a_t)>100.0_rk*small_.and.imode/=jmode) then 
                 write(out,"('Lmat_generation: Orthogon. eq. 4 is not = 0 for imode = ',i4,', jmode =',i4,d18.8)") imode,jmode,a_t
                 stop 'mat_generation: Amat is not Orthogonal'
            endif
         enddo
         !
      enddo
      !
      ! Quadratic constants in the  lineraized coordinates, defined by matrix Amat
      !
      f_t = matmul(matmul(transpose(c),f),c)
      !
      ! Diagonalization of "a" 
      !
      a = real(f_t,kind=rk)
      !
      call lapack_syev(a,b(:,1))
      !
      ! Found coordinate transformation "c"
      !
      c = real(a,kind=rk)
      !
      !a = matmul(transpose(c),c)
      ! 
      ! Transformation to an orthogonal Amat
      !
      do iatom = 1,Natoms
         do ix = 1,3 
            Amat(iatom,ix,:) = matmul(Amat(iatom,ix,:),c)
         enddo
      enddo
      !
      !trove%Amat = Amat
      !
      ! Check for the orthogonality sum_{Na} Amat_{Na l} Amat_{Na m} = 0 for m/=l
      !
      !
      !
      ! Check for the orthogonality sum_{Na} Amat_{Na l} Amat_{Na m} = 0 for m/=l
      !
      do imode = 1,Nmodes
         !
         do jmode = 1,Nmodes
            !
            a_t = sum(sum(Amat(:,:,jmode)*Amat(:,:,imode),dim=1))
            !
            if (abs(a_t-1.0_rk)>100.0_rk*small_.and.imode==jmode) then 
               write(out,"('Lmat_generation: Negative or zero diagonal element for Amat x Amat for imode = ',i4,d18.8)") imode,a_t
               stop 'mat_generation: Amat is not Orthogonal'
            endif
            !
            if (abs(a_t)>100.0_rk*small_.and.imode/=jmode) then 
                 write(out,"('Lmat_generation: Orthogon. eq. 4 is not = 0 for imode = ',i4,', jmode =',i4,d18.8)") imode,jmode,a_t
                 stop 'mat_generation: Amat is not Orthogonal'
            endif
         enddo
         !
      enddo
      !
      ! Checking quadratic constants diagonality 
      !
      f_t = matmul(matmul(transpose(c),f_t),c)
      !
      do imode =1,Nmodes
        omega(imode) = sqrt(2.0_ark*f_t(imode,imode)*factor*aJ_over_hc)
      enddo
      !
      do imode = 1,Nmodes
         !
         if (f_t(imode,imode)<small_) then 
              write(out,"('Lmat_gener: diagonal quadtrat. force constants are negative for imode = ',i4,d18.8)") imode,&
                    f_t(imode,imode)
              stop 'mat_generation: negative quadratic force constants'
         endif
         !
         do jmode = 1,Nmodes
            !
            if (abs(f_t(imode,jmode))/sqrt(f_t(imode,imode)*f_t(jmode,jmode))>1000.0_rk*small_.and.imode/=jmode) then 
                 write(out,"('Lmat_gener: quadtrat. force constants  not orthogonal for imode =',i4,', jmode =',i4,d18.8)") & 
                              imode,jmode,f_t(imode,jmode)
                 stop 'mat_generation: quadratic force constants are not orthogonal'
            endif
         enddo
         !
      enddo
      !
      ! Sorting up coordinate numbering to the standard ordering 
      ! i.e. 1-fold,2-fold,3-fold ... Within the same foldness  omega-s are given in the increasing order 
      !
      ifolds = 0
      lmode = 0
      lfolds = 0 
      foldness = 0 
      do imode = 1,Nmodes  
            kmode = 1  
            lmode = lmode +1  ; 
            ifolds(imode) = lmode
            do jmode = 1,Nmodes  
               if (abs(f_t(imode,imode)-f_t(jmode,jmode))/sqrt(f_t(imode,imode)*f_t(jmode,jmode))< &
               100.0_rk*sqrt(small_).and. &
                  imode/=jmode) then 
                  kmode = kmode + 1
                  ifolds(jmode) = lmode
               endif
            enddo 
         lfolds(kmode) = lfolds(kmode) +1
         foldness(kmode,lfolds(kmode)) = imode
      enddo 
      !
      !
      c = 0 
      imode = 0
      do ifold=1,Nmodes
         do kmode = lfolds(ifold),1,-1
           !
           imode = imode +1
           jmode = foldness(ifold,kmode)
           !
           c(jmode,imode) = 1.0_ark
           !
         enddo
      enddo 
      !
      f_t = matmul(matmul(transpose(c),f_t),c)
      !
      do imode =1,Nmodes
        omega(imode) = sqrt(2.0_ark*f_t(imode,imode)*factor*aJ_over_hc)
      enddo
      ! 
      ! Transformation to an orthogonal Amat
      !
      do iatom = 1,Natoms
         do ix = 1,3 
            Amat(iatom,ix,:) = matmul(Amat(iatom,ix,:),c)
         enddo
      enddo
      !
      ! Check if the Amat we get is the same as the usual lmat-matrix 
      ! in terms of the potential function expansion
      ! i.e. the quadratic PES is orthogonal 
      !
      ! We are building up the coordinate transformation from scratch, i.e. from Bmat through Amat
      !
      do imode = 1,Nmodes
         !
         do jmode = 1,Nmodes
            do ix = 1,3
               b_t(ix)=sum(Bmat(imode,:,ix)*Amat(:,ix,jmode)/sqrt(trove%mass(:)))
            enddo
            !
            ! the transformations
            !
            c(imode,jmode) = sum(b_t(:))
         enddo
         !
      enddo
      !
      ! New quadratic force constants in local coordinates are stored in "f"
      !
      f_t = matmul(matmul(transpose(c),f),c)
      !
      ! We double check once again that "f" is diagonal:
      !
      do imode = 1,Nmodes
         !
         do jmode = 1,Nmodes
            !
            if (abs(f_t(imode,jmode))/sqrt(f_t(imode,imode)*f_t(jmode,jmode))>1000.0_rk*small_.and.imode/=jmode) then 
                 write(out,"('Lmat_gener: quadtratic force constants  not orthogonal for imode = ',i4,', jmode =',i4,d18.8)") & 
                              imode,jmode,f_t(imode,jmode)
                 stop 'mat_generation: quadratic forcre constants are not orthogonal'
            endif
         enddo
         !
      enddo
      !
      ! here we calculate  harmonic frequencies omega 
      ! We need to adjust the dimensions to get 1/cm 
      ! 1) devide by h c
      ! 2) transform from Q to q by Q = q sqrt(hbar/2pi c omega)
      ! 3) devide by Angstrom^2 = 10^-8 and by avogno [gram]
      !
      do imode =1,Nmodes
        omega(imode) = sqrt(2.0_ark*f_t(imode,imode)*factor*aJ_over_hc)
      enddo
      !
      ! Now we can transform the quadratic potential parameters to the normal coordinates form 
      !
      do imode =1,Nmodes
         do jmode =1,Nmodes
            f_t(imode,jmode) = 2.0_ark*f_t(imode,jmode)*factor/sqrt(omega(imode)*omega(jmode))*aJ_over_hc
         enddo
      enddo
      !
      ! Check if the obtained normal force constants agree with those given at the input
      !
      !do imode =1,Nmodes
      !  do jmode =1,Nmodes
      !      powers = 0  
      !      powers(imode) = 1 ; powers(jmode) = powers(jmode) + 1
      !      ivar = FLQindex(powers)
      !      if (abs(f_t(imode,jmode)-trove%poten%field(ivar))> &
      !              100.0*small_*sqrt(abs(f_t(imode,imode)*f_t(jmode,jmode)))) then
      !         write(out,"('Lmat_generation: Normal qwforce /= poten field for modes = ',2i4,2d18.8)") imode,jmode,f_t(imode,jmode),trove%poten%field(ivar)
      !         stop 'Lmat_generation: poten /= qwforce'
      !      endif
      !   enddo
      !enddo
      !
      ! Now we complete the coordinate transformation to the dimensionless normal coordinates 
      ! by multiplying the Amat by 1/sqrt(omega) anb 1/sqrt(m)
      !
      factor = planck*avogno*real(1.0d+16,kind=rk)/(4.0_ark*pi*pi*vellgt)
      !factor = planck/(4.0_rk*pi*pi*vellgt)
      !
      do imode =1,Nmodes
         do iatom =1,Natoms
            Amat(iatom,:,imode) = Amat(iatom,:,imode)*sqrt(factor/omega(imode))/sqrt(trove%mass(iatom))
         enddo
      enddo
      !   
      !
      ! Store the result 
      !
      trove%omega = omega
      !
      ! Last check for the orthogonality sum_{Na} Amat_{Na l} Amat_{Na m} = 0 for m/=l
      !
      do imode = 1,Nmodes
         !
         do jmode = 1,Nmodes
            !
            a_t = 0.0_ark
            do iatom = 1,Natoms
               a_t = a_t + sum(Amat(iatom,:,jmode)*Amat(iatom,:,imode)*trove%mass(iatom))
            enddo
            !
            !if (abs(a_t-1.0_rk)>100.0_rk*small_.and.imode==jmode) then 
            !   write(out,"('Lmat_generation: Negative or zero diagonal element for Amat x Amat for imode = ',i4,d18.8)") imode,a_t
            !     stop 'mat_generation: Amat is not Orthogonal'
            !endif
            !
            if (abs(a_t)>1000.0_rk*small_.and.imode/=jmode) then 
                 write(out,"('Lmat_generation: Orthogon. eq. 4 is not = 0 for imode = ',i4,', jmode =',i4,d18.8)") imode,jmode,a_t
                 stop 'Lmat_generation: Amat is not orthogonal'
            endif
         enddo
         !
      enddo
      !

      !
      ! Check if the Amat we get is the same as the usual lmat-matrix 
      ! in terms of the potential function expansion
      ! i.e. the quadratic PES is orthogonal 
      !
      ! We building up the coordinate transformation from scratch, i.e. from Bmat through Amat
      !
      do imode = 1,Nmodes
         !
         do jmode = 1,Nmodes
            do ix = 1,3
               b_t(ix)=sum(Bmat(imode,:,ix)*Amat(:,ix,jmode))
            enddo
            !
            ! the transformations
            !
            c(imode,jmode) = sum(b_t(:))
            !
         enddo
         !
      enddo
      !
      ! New quadratic force constants in local coordinates are stored in "f"
      !
      f_t = matmul(matmul(transpose(c),f),c)*aJ_over_hc*2.0_ark
      !
      deallocate(a,b,c,f_t)
      !
    end subroutine Lmat_qwforce

    !
  end subroutine Lmat_generation1d



  !
  ! Here we initialize the potential energy field
  !
  subroutine FLinitilize_Potential_original
    !
    integer(ik) :: alloc,iM,imode,i,iterm
    integer(ik) :: Nmodes,Natoms,Nbonds,Nangles,Npoints,irho,Kindex(trove%Nmodes)
    type(FLpolynomT),pointer     :: fl
    real(ark)                    :: q_eq(trove%Nmodes),xi_eq(trove%Nmodes)
    real(ark)                    :: df,factor
    real(ark)                    :: step(2,trove%Nmodes),rhostep
    real(ark)                    :: astep(2)

    if (job%verbose>=2) write(out,"(/'FLinitilize_Potential_original/start')")   
    !
    call TimerStart('Potential')
    ! 
    ! Parameters for the internal use 
    !
    Nmodes  = trove%Nmodes
    Natoms  = trove%Natoms
    Nbonds  = trove%Nbonds
    NAngles = trove%Nangles
    Npoints = trove%Npoints
    !
    ! The potential energy function within the calculations must be given 
    ! in the normal coordinates. However at the input it is suposed to be 
    ! in the geometrically defined coordinates. 
    ! Therefore 
    ! To initialize potential energy function ("poten"-object) we perform 
    ! 1. Coordinate transformation from GDC to normal coordinates 
    ! 2. Calculate derivatives (normal force constants) using finite difference method
    !
    !
    ! Allocation of the poten-field
    !
    allocate (trove%poten,stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate poten-field')") alloc
        stop 'FLinitilize_Potential_original, poten-field - out of memory'
    end if
    !
    fl => trove%poten
    call polynom_initialization(fl,trove%NPotOrder,trove%RangeOrder(trove%NPotOrder),Npoints,'poten')
    !
    ! The potential function is expanded at its minimum, which is zero for the normal coordinates
    !
    q_eq(:) = trove%chi_ref(:,0)
    !
    ! step size 
    !
    rhostep = trove%rhostep
    !
    astep = 1.0_ark
    !
    do imode=1,trove%Nmodes_e
       !
       !call ML_check_steps4coordinvert(q_eq,2,imode,astep)
       !
       step(:,imode) = trove%fdstep(imode)*astep(:)
       !
       xi_eq(imode) = MLcoord_direct(q_eq(imode),2,imode)
       !
    enddo
    !
    if (job%verbose>=2) then
      !
      write(out,"('generated step-left :',40f14.6)") step(1,1:min(40,trove%Nmodes_e))
      write(out,"('generated step-right:',40f14.6)") step(2,1:min(40,trove%Nmodes_e))
      !
    endif
    !
    do irho=0,Npoints
      !
      FLirho = irho
      !
      if (job%verbose>=2.and.mod(irho,max(Npoints/20,1))==0) write(out,"('irho= ',i5)") irho
      !
      if (manifold/=0) then 
         !
         q_eq(Nmodes) = trove%rho_i(irho)
         xi_eq(Nmodes) = MLcoord_direct(q_eq(Nmodes),2,Nmodes)
      endif 
      !
      !$omp parallel do private(iterm,kindex,factor,imode,iM,df) schedule(guided)
      do iterm = 1,fl%Ncoeff
         !
         kindex(:) = FLIndexQ(:,iterm)
         !
         if (job%verbose>=5) write(out,"(10i8)") iterm,kindex(:)
         !
         ! Factorial factors to convert derivatives to the force constants 
         !
         factor = 1.0_ark
         do imode = 1,Nmodes
           do iM = 1,kindex(imode)
              factor = factor*real(iM,kind=rk)
           enddo
         enddo 
         !
         if (trove%internal_coords=='LOCAL') then 
            !
            ! Here we calculate derivatives by the finite differences 
            ! of the function "poten_local"
            ! with respect to local coordinates r1^k1 r2^k2 r3^k3 ...
            ! at the equilibrium given by q_eq,
            ! while k1,k2,k3... are stored in "kindex".
            ! fdstep defines the finite differences spacings 
            !
            df = FLfinitediffs(kindex,poten_xi,xi_eq,step(2:1:-1,:)) 
            !
         else
            !
            ! Here we calculate derivatives by the finite differences 
            ! of the function "poten_normal"
            ! with respect to normal coordinates q^k1 q^k2 q^k3 ...
            ! at the equilibrium given by q_eq, which is zero, 
            ! while k1,k2,k3... are stored in "kindex".
            ! fdstep defines finite differences spacings 
            !
            df = FLfinitediffs(kindex,poten_normal,xi_eq,step) 
            !
         endif 
         !
         trove%poten%field(iterm,irho) = df/factor
         !
      enddo
      !$omp end parallel do
      !
    enddo
    !
    if (job%verbose>=4.or.(job%verbose>=2.and.manifold==0)) then
       !
       write(out,"('pseudo-potential and potential parameteres:')")
       !
       do i = 1,max(trove%pseudo%Ncoeff,trove%poten%Ncoeff)
          !
          if (i<=min(trove%pseudo%Ncoeff,trove%poten%Ncoeff)) then 
             !
             do irho=0,Npoints,1
                write(out,"(20x,2i5,2g24.8,30i4)") irho,i,trove%pseudo%field(i,irho),trove%poten%field(i,irho),&
                                                 (FLIndexQ(imode,i),imode=1,min(30,Nmodes))
             enddo
             !
          elseif (i<trove%pseudo%Ncoeff) then 
             !
             do irho=0,Npoints,1
                !
                write(out,"(20x,2i5,g24.8,18x,30i4)") irho,i,trove%pseudo%field(i,irho),(FLIndexQ(imode,i),imode=1,min(30,Nmodes))
                !
             enddo
             !
          else
             !
             do irho=0,Npoints,1
                write(out,"(20x,2i5,18x,g24.8,30i4)") irho,i,trove%poten%field(i,irho),(FLIndexQ(imode,i),imode=1,min(30,Nmodes))
             enddo
             !
          endif
          !
       enddo
       !
    endif
    !
    if (trove%sparse) then
      !
      call FLCompact_a_field_sparse(trove%poten,"poten")
      !
      if (job%verbose>=5.or.(job%verbose>=2.and.manifold==0)) then
        !
        write(out,"('After compacting:')")
        call print_poten
      endif
      !
    endif
    !
    if (trim(trove%IO_potential)=='SAVE'.and.&
        trove%separate_store.and..not.trove%separate_convert) then
          !
          ! save potential.chk as ASCII
          !
          call FLcheck_point_Hamiltonian('POTENTIAL_SAVE_SPARSE')
          !      
    endif 
    !
    call TimerStop('Potential')
    !
    call TimerReport
    !
    if (job%verbose>=4) call MemoryReport
    !
    if (job%verbose>=2) write(out,"('FLinitilize_Potential_original/end')")   
    !
  end subroutine FLinitilize_Potential_original


  !
  ! Here we initialize the potential energy field
  !
  subroutine FLinitilize_Potential_II
    !
    integer(ik) :: alloc,iM,imode,i,iterm,N_potpoint
    integer(ik) :: Nmodes,Natoms,Nbonds,Nangles,Npoints,irho,Kindex(trove%Nmodes)
    type(FLpolynomT),pointer     :: fl
    real(ark)                    :: q_eq(trove%Nmodes),xi_eq(trove%Nmodes),xi_et(trove%Nmodes)
    real(ark)                    :: df,factor
    real(ark)                    :: step(2,trove%Nmodes),rhostep
    real(ark)                    :: astep(2),xi_t(trove%Nmodes),aq_et(trove%Nmodes)
    integer(ik)                  :: istep(2,trove%Nmodes)
    integer(ik),allocatable      :: ipoint_address(:,:)
    real(ark),allocatable        :: pot_points(:),poten_t(:)
    !
    if (trim(trove%IO_hamiltonian)=='READ'.or.&
        trim(trove%IO_potential)=='READ') then 
        !
        if (trim(trove%IO_kinetic)/='READ'.and.trim(trove%IO_hamiltonian)/='READ'.and..not.trove%separate_store) &
            call FLcheck_point_Hamiltonian('KINETIC_SKIP') 
        !
        call FLcheck_point_Hamiltonian('POTENTIAL_READ')
        !
        if (job%verbose>=6.or.(job%verbose>=2.and.manifold==0).or.(job%verbose>=5.and.trove%sparse)) then
           call print_poten
        endif 
        !
        if (trove%sparse) then
           !
          call FLCompact_a_field_sparse(trove%poten,"poten")
          !
          if (job%verbose>=5.or.(job%verbose>=2.and.manifold==0)) then
            !
            write(out,"('After compacting:')")
            call print_poten
          endif
          !
        endif
        !
        return 
        !
    endif
    !
    if (job%verbose>=2) write(out,"(/'FLinitilize_Potential_II/start')")   

    if (job%verbose>=4) write(out,"('  This is a simple finite differences type of the potential expansion. ')") 
    !
    call TimerStart('Potential')
    ! 
    ! Parameters for the internal use 
    !
    Nmodes  = trove%Nmodes
    Natoms  = trove%Natoms
    Nbonds  = trove%Nbonds
    NAngles = trove%Nangles
    Npoints = trove%Npoints
    !
    ! The potential energy function within the calculations must be given 
    ! in the normal coordinates. However at the input it is suposed to be 
    ! in the geometrically defined coordinates. 
    ! Therefore 
    ! To initialize potential energy function ("poten"-object) we perform 
    ! 1. Coordinate transformation from GDC to normal coordinates 
    ! 2. Calculate derivatives (normal force constants) using finite difference method
    !
    !
    ! Allocation of the poten-field
    !
    allocate (trove%poten,stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate poten-field')") alloc
        stop 'FLinitilize_Potential_II, poten-field - out of memory'
    end if
    !
    fl => trove%poten
    call polynom_initialization(fl,trove%NPotOrder,trove%RangeOrder(trove%NPotOrder),Npoints,'poten')
    !
    ! The potential function is expanded at its minimum, which is zero for the normal coordinates
    !
    q_eq(:)= trove%chi_ref(:,0)
    !
    ! step size 
    !
    rhostep = trove%rhostep
    !
    istep(:,:) = 1
    !
    astep = 1.0_ark
    !
    do imode=1,trove%Nmodes_e
       !
       !call ML_check_steps4coordinvert(q_eq,2,imode,astep)
       !
       step(:,imode) = trove%fdstep(imode)*astep(:)
       !
       if (astep(1)<small_) istep(1,imode) = 0 
       if (astep(2)<small_) istep(2,imode) = 0 
       !
       xi_eq(imode) = MLcoord_direct(q_eq(imode),2,imode)
       !
    enddo
    !
    xi_eq(Nmodes) = MLcoord_direct(q_eq(Nmodes),2,Nmodes)
    !
    ! Obtain adresses of all points needed for the finite derivatives
    ! 1. Estimate total number of all points required:
    !
    do imode = 1,trove%Nmodes
      !
      kindex(imode) = min(trove%NPotOrder,maxval(FLIndexQ(imode,:),dim=1))
      !
    enddo 
    !
    call FLdiffsadresses(kindex,istep,N_potpoint)
    !
    allocate (ipoint_address(N_potpoint,Nmodes),poten_t(trove%RangeOrder(trove%NPotOrder)),stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i9,' trying to allocate array ipoint_address')") alloc
       stop 'FLinitilize_Potential_II, ipoint_address  - out of memory'
    end if
    !
    call FLdiffsadresses(kindex,istep,N_potpoint,ipoint_address)
    !
    if (job%verbose>=2) then
      !
      write(out,"('generated step-left :',40f14.6)") step(1,1:min(40,trove%Nmodes_e))
      write(out,"('generated step-right:',40f14.6)") step(2,1:min(40,trove%Nmodes_e))
      !
    endif
    !
    !omp parallel private(pot_points,alloc)
    allocate (pot_points(N_potpoint),stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i9,' trying to allocate array pot_points')") alloc
       stop 'FLinitilize_Potential_II, pot_points  - out of memory'
    end if
    !
    xi_et = xi_eq
    !
    !omp do private(iterm,irho,kindex,factor,imode,iM,df,i) schedule(dynamic)
    do irho=0,Npoints
      !
      if (job%verbose>=2.and.mod(irho,max(Npoints/20,1))==0) write(out,"('irho= ',i5)") irho
      !
      if (manifold/=0) then 
         !
         aq_et(:) = trove%chi_ref(:,irho)
         !
         do imode = 1,Nmodes
           !
           xi_et(imode) = MLcoord_direct(aq_et(imode),1,imode)
           !
         enddo
         !
         !q_eq(Nmodes)  = trove%rho_i(irho)
         !xi_et(Nmodes) = MLcoord_direct(q_eq(Nmodes),2,Nmodes)
         !
      endif 
      !
      !$omp parallel do private(iterm,xi_t) shared(pot_points) schedule(dynamic)
      do iterm = 1,N_potpoint
         !
         xi_t(:) = xi_et(:) + trove%fdstep(:)*real(ipoint_address(iterm,:),ark)
         !
         if (trove%internal_coords=='LOCAL') then 
            !
            !stop 'check local potential here'
            !
            pot_points(iterm) = poten_xi(xi_t)
            !
         else
            !
            pot_points(iterm) = FLpoten_linearized(xi_t,irho)
            !
         endif 
         !
      enddo 
      !$omp end parallel do
      !
      do iterm = trove%poten%Ncoeff,1,-1
         !
         kindex(:) = FLIndexQ(:,iterm)
         !
         if (job%verbose>=6) write(out,"(10i8)") iterm,kindex(:)
         !
         ! Factorial factors to convert derivatives to the force constants 
         !
         factor = 1.0_ark
         do imode = 1,Nmodes
           do iM = 1,kindex(imode)
              factor = factor*real(iM,kind=ark)
           enddo
         enddo 
         !
         ! Here we calculate derivatives by the finite differences 
         ! of the function "poten_normal"
         ! with respect to normal coordinates q^k1 q^k2 q^k3 ...
         ! at the equilibrium given by q_eq, which is zero, 
         ! while k1,k2,k3... are stored in "kindex".
         ! fdstep defines finite differences spacings 
         !
         df = FLfinitediffs_precomp(kindex,ipoint_address,pot_points,trove%fdstep,istep) 
         !
         poten_t(iterm) = df/factor
         !
         !$omp parallel do private(i,xi_t,factor) shared(pot_points) schedule(dynamic)
         do i=1,N_potpoint
            !
            xi_t(:) = trove%fdstep(:)*real(ipoint_address(i,:),ark)
            !
            factor = product(xi_t(:)**kindex(:))
            !
            pot_points(i) = pot_points(i) - poten_t(iterm)*factor
            !
         enddo
         !$omp end parallel do
         !
      enddo
      !
      trove%poten%field(:,irho) = poten_t(:)
      !
    enddo
    !omp end do
    !
    deallocate(pot_points)
    !omp end parallel
    !
    deallocate(ipoint_address,poten_t)
    !
    if (job%verbose>=4.or.(job%verbose>=2.and.manifold==0)) then
       !
       write(out,"('pseudo-potential and potential parameteres:')")
       !
       do i = 1,max(trove%pseudo%Ncoeff,trove%poten%Ncoeff)
          !
          if (i<=min(trove%pseudo%Ncoeff,trove%poten%Ncoeff)) then 
             !
             do irho=0,Npoints,1
                write(out,"(20x,2i5,2f24.8,30i4)") irho,i,trove%pseudo%field(i,irho),trove%poten%field(i,irho),&
                                                 (FLIndexQ(imode,i),imode=1,min(30,Nmodes))
             enddo
             !
          elseif (i<trove%pseudo%Ncoeff) then 
             !
             do irho=0,Npoints,1
                !
                write(out,"(20x,2i5,f24.8,18x,30i4)") irho,i,trove%pseudo%field(i,irho),(FLIndexQ(imode,i),imode=1,min(30,Nmodes))
                !
             enddo
             !
          else
             !
             do irho=0,Npoints,1
                write(out,"(20x,2i5,24x,f24.8,30i4)") irho,i,trove%poten%field(i,irho),(FLIndexQ(imode,i),imode=1,min(30,Nmodes))
             enddo
             !
          endif
          !
       enddo
       !
    endif
    !
    fl => trove%poten
    call check_field_smoothness(fl,'CHECK',npoints,'poten')
    if (FL_iron_field_out) call check_field_smoothness(fl,'FIX',npoints,'poten')
    !
    if (job%verbose>=6.or.(job%verbose>=2.and.manifold==0).or.(job%verbose>=5.and.trove%sparse)) then
       call print_poten
    endif
    !
    if (trove%sparse) then
      !
      call FLCompact_a_field_sparse(trove%poten,"poten")
      !
      if (job%verbose>=5.or.(job%verbose>=2.and.manifold==0)) then
        !
        write(out,"('After compacting:')")
        call print_poten
      endif
      !
    endif
    !
    if (trim(trove%IO_potential)=='SAVE'.and.&
        trove%separate_store.and..not.trove%separate_convert) then
          !
          ! save potential.chk as ASCII
          !
          call FLcheck_point_Hamiltonian('POTENTIAL_SAVE_SPARSE')
          !      
    endif 
    !
    call TimerStop('Potential')
    !
    call TimerReport
    !
    if (job%verbose>=4) call MemoryReport
    !
    if (job%verbose>=2) write(out,"('FLinitilize_Potential_II/end')")   
    !
  end subroutine FLinitilize_Potential_II



   subroutine FLinitilize_Potential  
    !
    integer(ik) :: alloc,alloc_p,iM,imode,i,iterm,N_potpoint
    integer(ik) :: Nmodes,Natoms,Nbonds,Nangles,Npoints,irho,Kindex(trove%Nmodes)
    type(FLpolynomT),pointer     :: fl
    real(ark)                    :: q_eq(trove%Nmodes),xi_eq(trove%Nmodes),xi_et(trove%Nmodes)
    real(ark)                    :: df,factor
    real(ark)                    :: step(2,trove%Nmodes),rhostep
    real(ark)                    :: astep(2),xi_t(trove%Nmodes),aq_et(trove%Nmodes)
    integer(ik)                  :: istep(2,trove%Nmodes),jpar
    integer(ik),allocatable      :: ipoint_address(:,:)
    real(ark),allocatable        :: pot_points(:),poten_t(:)
    integer(ik)                  :: par(0:trove%Nmodes+1) ! Parity of each mode in derivatives calculation
    logical                      :: par_i
    character(len=cl)            :: my_fmt  !format for I/O specification

    if (job%verbose>=2) write(out,"(/'FLinitilize_Potential/start')")   
    if (job%verbose>=4) write(out,"('  Default expansion with finite diff, extended to N+1 with (N+1)s terms is set to zero')")
    !
    ! If the potentil function has been stored we can just read it from the hard disk and leave...
    !
    if (trim(trove%IO_hamiltonian)=='READ'.or.&
        trim(trove%IO_potential)=='READ') then 
        !
        if (trim(trove%IO_kinetic)/='READ'.and.trim(trove%IO_hamiltonian)/='READ'.and..not.trove%separate_store) &
            call FLcheck_point_Hamiltonian('KINETIC_SKIP') 
        !
        call FLcheck_point_Hamiltonian('POTENTIAL_READ')
        !
        if (job%verbose>=6.or.(job%verbose>=2.and.manifold==0).or.(job%verbose>=5.and.trove%sparse)) then
           call print_poten
        endif 
        !
        if (trove%sparse) then
           !
          call FLCompact_a_field_sparse(trove%poten,"poten")
          !
          if (job%verbose>=5.or.(job%verbose>=2.and.manifold==0)) then
            !
            write(out,"('After compacting:')")
            call print_poten
          endif
          !
        endif
        !
        return 
        !
    endif
    !
    call TimerStart('Potential')
    ! 
    ! Parameters for the internal use 
    !
    Nmodes  = trove%Nmodes
    Natoms  = trove%Natoms
    Nbonds  = trove%Nbonds
    NAngles = trove%Nangles
    Npoints = trove%Npoints
    !
    ! The potential energy function within the calculations must be given 
    ! in the normal coordinates. However at the input it is suposed to be 
    ! in the geometrically defined coordinates. 
    ! Therefore 
    ! To initialize potential energy function ("poten"-object) we perform 
    ! 1. Coordinate transformation from GDC to normal coordinates 
    ! 2. Calculate derivatives (normal force constants) using finite difference method
    !
    !
    ! Allocation of the poten-field
    !
    allocate (trove%poten,stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate poten-field')") alloc
        stop 'FLinitilize_Potential, poten-field - out of memory'
    end if
    !
    fl => trove%poten
    call polynom_initialization(fl,trove%NPotOrder,trove%RangeOrder(trove%NPotOrder),Npoints,'poten')
    !
    ! The potential function is expanded at its minimum, which is zero for the normal coordinates
    !
    q_eq(:)= trove%chi_ref(:,0)
    !
    ! step size 
    !
    rhostep = trove%rhostep
    !
    istep(:,:) = 1
    !
    do imode=1,trove%Nmodes_e
       !
       call ML_check_steps4coordinvert(q_eq,2,imode,astep)
       !
       step(:,imode) = trove%fdstep(imode)*astep(:)
       !
       if (astep(1)<small_) istep(1,imode) = 0 
       if (astep(2)<small_) istep(2,imode) = 0 
       !
       xi_eq(imode) = MLcoord_direct(q_eq(imode),2,imode)
       !
    enddo
    !
    !
    if (job%verbose>=2) then
      !
      write(out,"('generated step-left :',40f14.6)") step(1,1:min(40,trove%Nmodes_e))
      write(out,"('generated step-right:',40f14.6)") step(2,1:min(40,trove%Nmodes_e))
      !
    endif
    !
    xi_eq(Nmodes) = MLcoord_direct(q_eq(Nmodes),2,Nmodes)
    !
    ! Initialize parity
    !
    do imode=1,trove%Nmodes
      !
      if (istep(1,imode)==0.or.istep(2,imode)==0) then
        !
        par(imode)=-1 ! For this mode only one parity is possible
        !
      else
        !
        par(imode)=0  ! For this mode both parities are possible
        !
      endif
      !
    enddo
    !
    par(trove%Nmodes+1)=0
    !
    ! Start a parity cycle. All posible parities must be processed
    !
    if (job%verbose>=2) write(out,"(/'Number of different types of finite diff. points :',i8)") 2**trove%nmodes
    !
    ! Count different types of parity 
    !
    jpar = 0 
    write(my_fmt,'(a,i0,a)') "(i8,a,",Nmodes,"i4)"
    !
    do while(par(trove%Nmodes+1)==0)
      !
      jpar = jpar + 1
      !
      if (job%verbose>=3) write(out,my_fmt) jpar,' -> ',par(1:trove%Nmodes)
      !
      do imode = 1,trove%Nmodes
        !
        kindex(imode) = min(trove%NPotOrder,maxval(FLIndexQ(imode,:),dim=1))
        !
      enddo 
      !
      ! Obtain adresses of all points needed for the finite derivatives
      !
      ! 1. Estimate total number of all points required for a current parity
      !
      call FLdiffsadresses_roman(kindex,istep,N_potpoint,par)
      !
      if (N_potpoint/=0) then
        !
        allocate (ipoint_address(N_potpoint,Nmodes),stat=alloc)
        if (alloc/=0) then
           write (out,"(' Error ',i9,' trying to allocate array ipoint_address')") alloc
           stop 'FLinitilize_Potential, ipoint_address  - out of memory'
        end if
        !
        call FLdiffsadresses_roman(kindex,istep,N_potpoint,par,ipoint_address)
        !
        !$omp parallel private(pot_points,poten_t,alloc_p)
        allocate (pot_points(N_potpoint),poten_t(trove%RangeOrder(trove%NPotOrder)),stat=alloc_p)
        if (alloc_p/=0) then
           write (out,"(' Error ',i9,' trying to allocate array pot_points')") alloc_p
           stop 'FLinitilize_Potential, pot_points  - out of memory'
        end if
        !
        !$omp do private(irho,aq_et,imode,xi_et,iterm,xi_t,kindex,factor,par_i,iM,df,i) schedule(dynamic)
        do irho=0,Npoints
          !
          if (job%verbose>=6.and.mod(irho,max(Npoints/100,1))==0.and.Npoints/=0) write(out,"('irho= ',i5)") irho
          !
          xi_et = xi_eq
          !
          if (manifold/=0) then 
             !
             aq_et(:) = trove%chi_ref(:,irho)
             !
             do imode = 1,Nmodes
               !
               xi_et(imode) = MLcoord_direct(aq_et(imode),2,imode)
               !
             enddo
             !
          endif 
          !
          do iterm = 1,N_potpoint
             !
             xi_t(:) = xi_et(:) + trove%fdstep(:)*real(ipoint_address(iterm,:),ark)
             !
             if (trove%internal_coords=='LOCAL') then 
                !
                !stop 'check local potential here'
                !
                pot_points(iterm) = poten_xi(xi_t)
                !
             else
                !
                pot_points(iterm) = FLpoten_linearized(xi_t,irho)
                !
             endif 
             !
          enddo 
          !
          do iterm = trove%poten%Ncoeff,1,-1
           !
           kindex(:) = FLIndexQ(:,iterm)
           !
           par_i = par_check(par,kindex)
           !
           ! If parity of kindex (polynom coefficients and derivatives) is the 
           ! same as parity of points
           !
           if (par_i) then 
               !
               if (job%verbose>=7) write(out,"(10i8)") iterm,kindex(:)
               !
               ! Factorial factors to convert derivatives to the force constants 
               !
               factor = 1.0_ark
               do imode = 1,Nmodes
                 do iM = 1,kindex(imode)
                    factor = factor*real(iM,kind=ark)
                 enddo
               enddo 
               !
               ! Here we calculate derivatives by the finite differences 
               ! of the function "poten_normal"
               ! with respect to normal coordinates q^k1 q^k2 q^k3 ...
               ! at the equilibrium given by q_eq, which is zero, 
               ! while k1,k2,k3... are stored in "kindex".
               ! fdstep defines finite differences spacings 
               !
               df = FLfinitediffs_precomp(kindex,ipoint_address,pot_points,trove%fdstep,istep) 
               !
               poten_t(iterm) = df/factor
               !
               !if (job%verbose>=6) then
               !  !
               !  write(out,"('iterm,poten_t :',i8,f14.6)") iterm,poten_t(iterm)
               !  !
               !endif
               !
               !omp parallel do private(i,xi_t,factor) shared(pot_points) schedule(dynamic)
               do i=1,N_potpoint
                  !
                  xi_t(:) = trove%fdstep(:)*real(ipoint_address(i,:),ark)
                  !
                  factor = product(xi_t(:)**kindex(:))
                  !
                  pot_points(i) = pot_points(i) - poten_t(iterm)*factor
                  !
                  !if (job%verbose>=6) then
                  !  !
                  !  write(out,"('i,pot_points :',i8,f14.6)") i,pot_points(i)
                  !  !
                  !endif
                  !
               enddo
               !omp end parallel do
               !
               trove%poten%field(iterm,irho) = poten_t(iterm)
               !
           endif ! parity check
           !
          enddo
          !
          !trove%poten%field(:,irho) = poten_t(:)
          !
        enddo
        !$omp end do
        !
        deallocate(pot_points,poten_t)
        !$omp end parallel
        !
        deallocate(ipoint_address)
        !
      endif ! N_potpoint/=0
      !
      ! Switching the parity
      !
      call par_switch(par)
      !
    enddo ! parithy cycle
    !
    fl => trove%poten
    call check_field_smoothness(fl,'CHECK',npoints,'poten')
    if (FL_iron_field_out) call check_field_smoothness(fl,'FIX',npoints,'poten')
    !
    if (job%verbose>=6.or.(job%verbose>=2.and.manifold==0).or.(job%verbose>=5.and.trove%sparse)) then
       call print_poten
    endif
    !
    if (trove%sparse) then
      !
      call FLCompact_a_field_sparse(trove%poten,"poten")
      !
      if (job%verbose>=5.or.(job%verbose>=2.and.manifold==0)) then
        !
        write(out,"('After compacting:')")
        call print_poten
      endif
      !
    endif
    !
    if (trim(trove%IO_potential)=='SAVE'.and.&
        trove%separate_store.and..not.trove%separate_convert) then
          !
          ! save potential.chk as ASCII
          !
          call FLcheck_point_Hamiltonian('POTENTIAL_SAVE_SPARSE')
          !      
    endif 
    !
    call TimerStop('Potential')
    !
    call TimerReport
    !
    if (job%verbose>=4) call MemoryReport
    !
    if (job%verbose>=2) write(out,"('FLinitilize_Potential/end')")  
    
  contains
    
    subroutine par_switch(par)
      !
      integer(ik),intent(inout)   :: par(0:trove%Nmodes+1)
      integer(ik)              :: imode
      !
      imode = 0
      !
      do while(par(imode)<=0.or.imode==0)
        !
        imode = imode + 1
        !
        ! if we possible we switch the parity of this mode
        !
        if (par(imode)>=0) par(imode) = 1 - par(imode)
        !
      enddo 
      !
    end subroutine par_switch



    logical function par_check(par,kindex)
      !
      integer(ik),intent(in)   :: par(0:trove%Nmodes+1),kindex(trove%Nmodes)
      integer(ik)              :: imode
      logical :: ch0
      !
      ch0 = .true.
      !
      do imode = 1,trove%Nmodes
        !
        ch0 = ch0 .and. (par(imode)==-1.or.mod(par(imode)+kindex(imode),2)==0)
        !
      enddo
      !
      par_check = ch0
      !
    end function par_check

    
  end subroutine FLinitilize_Potential


    subroutine print_poten
     !
     integer(ik) :: i,irho,imode
        !
        write(out,"('pseudo-potential parameteres:')")
        !
        do i = 1,trove%pseudo%Ncoeff
           !
           do irho=0,trove%Npoints,1
              !
              write(out,"(20x,i5,1x,i7,f24.8,18x,30i4)") irho,i,trove%pseudo%field(i,irho),(trove%pseudo%IndexQ(imode,i),&
                    imode=1,min(30,trove%Nmodes))
              !
           enddo
           !
        enddo
        !
        write(out,"('potential parameteres:')")
        !
        do i = 1,trove%poten%Ncoeff
           !
           do irho=0,trove%Npoints,1
              write(out,"(20x,i5,1x,i7,f24.8,18x,30i4)")  irho,i,trove%poten%field(i,irho),&
                                                         (trove%poten%IndexQ(imode,i),imode=1,min(30,trove%Nmodes))
           enddo
           !
        enddo
        !     
    end subroutine print_poten 


  subroutine FLinit_External_field  
    !
    integer(ik) :: alloc,alloc_p,iM,imode,i,iterm,N_extFpoint
    integer(ik) :: Nmodes,Natoms,Nbonds,Nangles,Npoints,irho,Kindex(trove%Nmodes)
    type(FLpolynomT),pointer     :: fl
    real(ark)                    :: q_eq(trove%Nmodes),xi_eq(trove%Nmodes),xi_et(trove%Nmodes)
    real(ark)                    :: df,factor
    real(ark)                    :: step(2,trove%Nmodes),rhostep
    real(ark)                    :: astep(2),xi_t(trove%Nmodes),aq_et(trove%Nmodes)
    integer(ik)                  :: istep(2,trove%Nmodes),jpar,rank,nterms,nterms_max,imu,Ncoeff,i1,i2,N_,k,ioutlier
    integer(ik),allocatable      :: ipoint_address(:,:)
    real(ark),allocatable        :: extF_points(:,:),extF_t(:)
    integer(ik)                  :: par(0:trove%Nmodes+1) ! Parity of each mode in derivatives calculation
    logical                      :: par_i,outliers
    integer(ik),parameter        :: Ni = 10, Nattempt = 100
    real(ark)                    :: rho_(1:Ni),func_(1:Ni),rho,fval,foutlier,fcorr
    integer(ik)                  :: iattempt,ioutlier_max
    character(len=4)             :: txt
    character(len=cl)            :: my_fmt !format for I/O specification
    !
    ! If the potentil function has been stored we can just read it from the hard disk and leave...
    !
    if (trim(trove%IO_ext_coeff)=='READ') then 
      !
      call FLcheck_point_Hamiltonian('EXTERNAL_READ') 
      !
      call print_external
      !
      if (trove%sparse) call compact_sparse_external
      !
      return
      !
    endif
    !
    npoints   = trove%Npoints
    nmodes    = trove%Nmodes
    rank      = extF%rank
    !
    if (.not.FLextF_coeffs)  return
    !
    if (job%verbose>=3) write(out,"(/'Expansion of the external field ...')")   
    !
    if (job%verbose>=1) then
      write (out,"(' External fields need ',f9.3,' Mbytes of memory (plus a bit)')") &
             real(rk*product(extF%maxord(:)),kind=rk)/(1024.0_rk**2)
    end if
    !
    allocate (trove%extF(extF%rank),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate extF-field')") alloc
        stop 'FLinit_External_field, extF-field - out of memory'
    end if
    !
    if (job%verbose>=5)  write(out, '(/1x, a)') 'generate expansion terms'
    !
    nterms_max = 1
    !
    do imu = 1,extF%rank
       !
       ! Number of the expansion coefficients
       !
       nterms  = trove%RangeOrder(extF%maxord(imu))
       !
       fl => trove%extF(imu)
       call polynom_initialization(fl,extF%maxord(imu),nterms,Npoints,'extF')
       !
       nterms_max = max(nterms_max,nterms)
       !
       write(my_fmt,'(a,i0,a)') "(/2(1x, a),",Nmodes,"(1x, i3))"
       !
       if (job%verbose>=7) then
          write(out, '(/1x, a, 25x, i3/1x, a, 4x, i3/1x, a, 1x, i5)') &
          'imu', imu, 'maximum expansion degree',trove%NExtOrder,'number of expansion terms', nterms
          write(out,my_fmt) 'iterm', 'imode->', (imode, imode = 1, nmodes)
       end if
       !
       !print expansion terms
       !
       write(my_fmt,'(a,i0,a)') "(3x, i6, 8x,",Nmodes,"(1x, i3))"
       !
       if (job%verbose>=5) then
          do iterm = 1, trove%extF(imu)%Ncoeff
             write(out,my_fmt) iterm, (FLIndexQ(imode,iterm), imode = 1, nmodes)
          end do
       end if
       !
    enddo
    !
    call TimerStart('External')
    !
    ! The potential function is expanded at its minimum, which is zero for the normal coordinates
    !
    q_eq(:)= trove%chi_ref(:,0)
    !
    ! step size 
    !
    rhostep = trove%rhostep
    !
    istep(:,:) = 1
    !
    do imode=1,trove%Nmodes_e
       !
       call ML_check_steps4coordinvert(q_eq,3,imode,astep)
       !
       step(:,imode) = extF%fdstep(imode)*astep(:)
       !
       if (astep(1)<small_) istep(1,imode) = 0 
       if (astep(2)<small_) istep(2,imode) = 0 
       !
       xi_eq(imode) = MLcoord_direct(q_eq(imode),3,imode)
       !
    enddo
    !
    !
    if (job%verbose>=2) then
      !
      write(out,"('generated step-left :',40f14.6)") step(1,1:min(40,trove%Nmodes_e))
      write(out,"('generated step-right:',40f14.6)") step(2,1:min(40,trove%Nmodes_e))
      !
    endif
    !
    xi_eq(Nmodes) = MLcoord_direct(q_eq(Nmodes),3,Nmodes)
    !
    ! Initialize parity
    !
    do imode=1,trove%Nmodes
      !
      if (istep(1,imode)==0.or.istep(2,imode)==0) then
        !
        par(imode)=-1 ! For this mode only one parity is possible
        !
      else
        !
        par(imode)=0  ! For this mode both parities are possible
        !
      endif
      !
    enddo
    !
    par(trove%Nmodes+1)=0
    !
    ! Start a parity cycle. All posible parities must be processed
    !
    if (job%verbose>=2) write(out,"(/'Number of different types of finite diff. points :',i8)") 2**trove%nmodes
    !
    ! Count different types of parity 
    !
    jpar = 0 
    !
    write(my_fmt,'(a,i0,a)') "(i8,a,",Nmodes,"i4)"
    !
    do while(par(trove%Nmodes+1)==0)
      !
      jpar = jpar + 1
      !
      if (job%verbose>=2) write(out,my_fmt) jpar,' -> ',par(1:trove%Nmodes)
      !
      do imode = 1,trove%Nmodes
        !
        kindex(imode) = min(trove%NExtOrder,maxval(FLIndexQ(imode,:),dim=1))
        !
      enddo 
      !
      ! Obtain adresses of all points needed for the finite derivatives
      !
      ! 1. Estimate total number of all points required for a current parity
      !
      call FLdiffsadresses_roman(kindex,istep,N_extFpoint,par)
      !
      if (N_extFpoint/=0) then
        !
        allocate (ipoint_address(N_extFpoint,Nmodes),stat=alloc)
        if (alloc/=0) then
           write (out,"(' Error ',i9,' trying to allocate array ipoint_address')") alloc
           stop 'FLinit_External_field, ipoint_address  - out of memory'
        end if
        !
        call FLdiffsadresses_roman(kindex,istep,N_extFpoint,par,ipoint_address)
        !
        !$omp parallel private(extF_points,extF_t,alloc_p)
        allocate (extF_points(N_extFpoint,rank),extF_t(nterms_max),stat=alloc_p)
        if (alloc_p/=0) then
           write (out,"(' Error ',i9,' trying to allocate array extF_points')") alloc_p
           stop 'FLinit_External_field, extF_points  - out of memory'
        end if
        !
        !$omp do private(irho,aq_et,imode,xi_et,iterm,xi_t,kindex,factor,par_i,iM,df,i,imu) schedule(dynamic)
        do irho=0,Npoints
          !
          if (job%verbose>=5.and.mod(irho,max(Npoints/100,1))==0.and.Npoints/=0) write(out,"('irho= ',i5)") irho
          !
          xi_et = xi_eq
          !
          if (manifold/=0) then 
             !
             aq_et(:) = trove%chi_ref(:,irho)
             !
             do imode = 1,Nmodes
               !
               xi_et(imode) = MLcoord_direct(aq_et(imode),3,imode)
               !
             enddo
             !
          endif 
          !
          do iterm = 1,N_extFpoint
             !
             xi_t(:) = xi_et(:) + extF%fdstep(:)*real(ipoint_address(iterm,:),ark)
             !
             call dms4xi(rank,nmodes,irho,xi_t,extF_points(iterm,:))
             !
          enddo 
          !
          do imu = 1,rank
            !
            do iterm = trove%extF(imu)%Ncoeff,1,-1
             !
             kindex(:) = FLIndexQ(:,iterm)
             !
             par_i = par_check(par,kindex)
             !
             ! If parity of kindex (polynom coefficients and derivatives) is the 
             ! same as parity of points
             !
             if (par_i) then 
                 !
                 if (job%verbose>=6) write(out,"(10i8)") iterm,kindex(:)
                 !
                 ! Factorial factors to convert derivatives to the force constants 
                 !
                 factor = 1.0_ark
                 do imode = 1,Nmodes
                   do iM = 1,kindex(imode)
                      factor = factor*real(iM,kind=ark)
                   enddo
                 enddo 
                 !
                 ! Here we calculate derivatives by the finite differences 
                 ! of the function "poten_normal"
                 ! with respect to normal coordinates q^k1 q^k2 q^k3 ...
                 ! at the equilibrium given by q_eq, which is zero, 
                 ! while k1,k2,k3... are stored in "kindex".
                 ! fdstep defines finite differences spacings 
                 !
                 df = FLfinitediffs_precomp(kindex,ipoint_address,extF_points(:,imu),extF%fdstep,istep)
                 !
                 extF_t(iterm) = df/factor
                 !
                 do i=1,N_extFpoint
                    !
                    xi_t(:) = extF%fdstep(:)*real(ipoint_address(i,:),ark)
                    !
                    factor = product(xi_t(:)**kindex(:))
                    !
                    extF_points(i,imu) = extF_points(i,imu) - extF_t(iterm)*factor
                    !
                 enddo
                 !
                 trove%extF(imu)%field(iterm,irho) = extF_t(iterm)
                 !
             endif ! parity check
             !
            enddo
            !
          enddo
          !
        enddo
        !$omp end do
        !
        deallocate(extF_points,extF_t)
        !$omp end parallel
        !
        deallocate(ipoint_address)
        !
      endif ! N_extFpoint/=0
      !
      ! Switching the parity
      !
      call par_switch(par)
      !
    enddo ! parithy cycle
    !
    ! report coefficients 
    !
    call print_external  
    !
    if (trove%sparse) call compact_sparse_external
    !
    ! check the smoothness and fix if necessary 
    !
    do imu = 1, extF%rank
      !
      write(txt,"(i4)") imu
      !
      fl => trove%extF(imu)
      !
      call check_field_smoothness(fl,'CHECK',npoints,'FLinit_External_field'//trim(txt))
      !
      !call check_field_smoothness(fl,'FIX',npoints,'FLinit_External_field'//trim(txt))
      !
    enddo
    !
    if (trim(trove%IO_ext_coeff)=='SAVE'.and.&
        trove%separate_store.and..not.trove%separate_convert) then
          !
          ! save external.chk as ASCII
          !
          call FLcheck_point_Hamiltonian('EXTERNAL_SAVE_SPARSE')
          !      
    endif 
    !
    call TimerStop('External')
    !
    call TimerReport
    !
    if (job%verbose>=4) call MemoryReport
    !
    if (job%verbose>=3) then
       write(out, '(/a)') '...done!'
    endif
    
  contains
  
  
  
   subroutine print_external  
  
      integer(ik) :: imu,irho,Ncoeff,npoints,nterms_max
      character(len=cl) :: my_fmt !format for I/O specification

      !print expansion coefficients and expansion points
      !
      if (job%verbose<5) return 
      !
      nterms_max = trove%extF(1)%Ncoeff
      Ncoeff = nterms_max
      npoints   = trove%Npoints
      !
      write(my_fmt,'(a,i0,a)') "(/1x, a, 1x, a, 1x, a,",Ncoeff,"(8x, i5))"
      !
      write(out,my_fmt) 'ipoint', 'imu', 'iterm->', (iterm, iterm = 1,nterms_max)
      !
      write(my_fmt,'(a,i0,a)') "(1x, i6, 2x, i4, 8x,",Ncoeff,"(1x, es12.4))"
      !
      do imu = 1, extF%rank
        do irho = 0, npoints
            !
            Ncoeff = trove%extF(imu)%Ncoeff
            !
            write(out,my_fmt) irho, imu, (trove%extF(imu)%field(iterm,irho), iterm = 1,Ncoeff)
         end do
      end do

     end subroutine print_external  


    subroutine compact_sparse_external
       !
       implicit none
       !
       integer(ik) :: imu
       type(FLpolynomT),pointer     :: fl
       !
       do imu = 1, extF%rank
          !
          fl => trove%extF(imu)
          if (trove%sparse) call FLCompact_a_field_sparse(fl,"extF")
          !
       enddo
       !
    end subroutine compact_sparse_external
  
    
    subroutine par_switch(par)
      !
      integer(ik),intent(inout)   :: par(0:trove%Nmodes+1)
      integer(ik)              :: imode
      !
      imode = 0
      !
      do while(par(imode)<=0.or.imode==0)
        !
        imode = imode + 1
        !
        ! if we possible we switch the parity of this mode
        !
        if (par(imode)>=0) par(imode) = 1 - par(imode)
        !
      enddo 
      !
    end subroutine par_switch


    logical function par_check(par,kindex)
      !
      integer(ik),intent(in)   :: par(0:trove%Nmodes+1),kindex(trove%Nmodes)
      integer(ik)              :: imode
      logical :: ch0
      !
      ch0 = .true.
      !
      do imode = 1,trove%Nmodes
        !
        ch0 = ch0 .and. (par(imode)==-1.or.mod(par(imode)+kindex(imode),2)==0)
        !
      enddo
      !
      par_check = ch0
      !
    end function par_check

    
  end subroutine FLinit_External_field


   !
   ! Finite differences scheme, which uses the precumputed potential energy points 
   !

   function FLfinitediffs_precomp(itarget,ipoint_address,pot_points,step,istep) result(f)

      integer(ik),intent(in) :: itarget(trove%Nmodes)
      integer(ik),intent(in) :: ipoint_address(:,:)
      real(ark),intent(in)   :: pot_points(:)
      real(ark),intent(in)    :: step(trove%Nmodes)
      integer(ik),intent(in) :: istep(2,trove%Nmodes)
      !
      integer(ik)            :: i_potpoints
      real(ark)              :: f,astep(trove%Nmodes)
      integer(ik) :: imode,Nmodes
      integer(ik) :: isearch(trove%Nmodes)
      integer     :: jindex(trove%Nmodes)
      character(len=cl) :: my_fmt !format for I/O specification
      !
      Nmodes = trove%Nmodes
      !
      write(my_fmt,'(a,i0,a)') "(a,",Nmodes,"i4)"
      !
      if (verbose>=6) write(out,my_fmt) 'FLfinitediffs/start: finite derivatives for k=',itarget(1:Nmodes)
      !
      i_potpoints = size(ipoint_address,dim=1)
      !
      Isearch =  Itarget
      jindex = 0
      !
      astep = step
      !
      ! here we go recursively untill we hit the bottom, i.e. untill isearch = 0 
      !
      f = fdiff(isearch,ipoint_address,pot_points,i_potpoints,astep,istep,jindex) 

      if (verbose>=6) write(out,"('FLfinitediffs/end')") 
      !
    contains 

     recursive function fdiff(isearch,ipoint_address,pot_points,i_potpoints,astep,istep,jindex) result (f)

       integer(ik),intent(inout) :: isearch(trove%Nmodes)
       integer(ik),intent(in) :: ipoint_address(:,:)
       real(ark),intent(in)   :: pot_points(:)
       integer(ik),intent(in) :: i_potpoints
       real(ark),intent(in)   :: astep(trove%Nmodes)
       integer(ik),intent(in) :: istep(2,trove%Nmodes)
       integer,intent(in)     :: jindex(trove%Nmodes)
       integer                :: kindex(trove%Nmodes),sig
       !
       real(ark)              :: f,h,f1,f2,f3,f4
       integer(ik)            :: k01(trove%Nmodes),k02(trove%Nmodes),k03(trove%Nmodes),k04(trove%Nmodes)
       integer(ik)            :: imode,jmode,index1,index2,iterm
       !
       ! find first active coordinate 
       !
       imode = 1
       !
       do while(Isearch(imode)==0.and.imode<trove%Nmodes)
         imode = imode+1   
       enddo
       !
       if (imode == trove%Nmodes.and.isearch(imode)==0) then
          !
          index1=1
          index2=i_potpoints+1
          sig=2
          ! searching
          do while(sig/=0)
           iterm=(index1+index2)/2
           kindex(:)=ipoint_address(iterm,:)
           sig=0
           do jmode=1,trove%Nmodes
            if(sig==0.and.jindex(jmode)>kindex(jmode)) sig=1
            if(sig==0.and.jindex(jmode)<kindex(jmode)) sig=-1
           enddo
           if (sig== 1) index1=iterm
           if (sig==-1) index2=iterm
          enddo
          !
          f = pot_points(iterm)
          !
       else 
          !
          ! if not all derivatives have been evaluated, we continue to the next lower recursive level
          !
          ! here we go one level lower for the coordinate imode 
          !
          Isearch(imode) = Isearch(imode)-1
          !
          ! we use simple central formula for the finite difference derivatives  
          !
          select case(difftype)
          case default
             !
             write (out,"('FLfinitediffs: difftype ',i8,' unknown')") difftype
             stop 'FLfinitediffs - bad difftype'
             !
          case(2)
             !
             k02 = jindex ; k01 = jindex 
             k01(imode) = jindex(imode) - istep(1,imode)
             k02(imode) = jindex(imode) + istep(2,imode)
             !
             f1 = fdiff(isearch,ipoint_address,pot_points,i_potpoints,astep,istep,k01)
             f2 = fdiff(isearch,ipoint_address,pot_points,i_potpoints,astep,istep,k02)
             !
             h = astep(imode)*real(istep(1,imode)+istep(2,imode),ark)
             !
             f = ( f2-f1 )/(h)
             !
             !f_t = 0.5_ark*(f2-f1)
             !f=log(f_t+sqrt(f_t**2+1.0_ark))/astep(1,imode)
             !
!          case(4)
             !
!             x1 = x ; x2 = x ;  x3 = x ; x4 = x 
             !
!             x1(imode) = x(imode) - astep(1,imode)
!             x2(imode) = x(imode) + astep(2,imode)
!             x3(imode) = x(imode) - astep(1,imode)*2.0_ark
!             x4(imode) = x(imode) + astep(2,imode)*2.0_ark
             !
!             h = astep(1,imode)+astep(2,imode)
             !
!             f1 = fdiff(isearch,x1,func)
!             f2 = fdiff(isearch,x2,func)
!             f3 = fdiff(isearch,x3,func)
!             f4 = fdiff(isearch,x4,func)
             !
             !write (out,"('FLfinitediffs-4: You should check the 4-point equation first!')") 
             !stop 'FLfinitediffs-4: Check the 4-point equation!'
             !
!             f = (-f4+8.0_ark*f2+f3-8.0_ark*f1 )/(6.0_ark*h)

             !f = (-f4/12.0_ark+2.0_ark/3.0_ark*f2 & 
             !     +f3/12.0_ark-2.0_ark/3.0_ark*f1 )/h

             !
          end select
          !
          ! here we are back to the current level and the "imode"-coordinate is also back 
          !
          Isearch(imode) = Isearch(imode)+1
          !
       endif 
       !
     end function fdiff
 

    !
   end function FLfinitediffs_precomp


  !
  !
  !
  subroutine Bmat_generation(a0,Bmat)
    !
    real(ark),intent(in)   ::  a0(trove%Natoms,3)
    real(ark),intent(out)  ::  Bmat(trove%Ncoords,trove%Natoms,3) 

    real(ark)   :: rcon(trove%Natoms,trove%Natoms)
    real(ark)   :: acon(trove%Natoms,trove%Natoms,trove%Natoms)
    real(ark)   :: r1,r2,r3,cosa,sina,ddelta_t(trove%Natoms,3)
    real(ark)   :: dvec1(3),dvec2(3),Bmat_t(trove%Natoms,3)
    real(ark)   :: da_dxna(trove%Natoms,trove%Natoms,trove%Natoms,trove%Natoms,3) 
    real(ark)   :: vec1,vec2,B,d1_dy_i_x(3),d2_dy_i_x(3),a1,a2,a3,b_t(3)

    integer(ik) ::  Nmodes,Natoms,Nbonds,Nangles,Ndihedrals
    integer(ik) ::  ibond,iangle,ida(trove%Natoms,trove%Natoms,trove%Natoms),iatom
    integer(ik) ::  n1,n2,n3,n4,n0,ix,iy,iz,iq,J,m,n,o,p,a,mnop(4),i,kappa,icoord,k1,k2
    !
    real(ark)   :: a_t1(3),a_t2(3),a_t3(3),delta,deltaf(trove%Natoms,trove%Natoms),zetaf(trove%Natoms,&
                   trove%Natoms,trove%Natoms),dB(3,3),tau_sign,dnorm_da1,dnorm_da2,dnorm_da3,norm_2
    real(ark)   :: cosa1,cosa2,cosa3,sina1,sina2,sina3,cosu,cosv,sinu,sinv,u(3),v(3),w(3),cart_vec(3,3),sindelta,phi,e1(3),e2(3),f_t
    !
    if (job%verbose>=7) write(out,"(/'Bmat_generation/start  ')") 


    !
    ! for (102) and (103,104) cases fo numerical derivatives 
    ! instead of more acurate analytical 
    !
    if ( any(trove%dihedtype(:)==101).or.any(trove%dihedtype(:)==102).or.any(trove%dihedtype(:)==103).or.&
         any(trove%dihedtype(:)==104).or.any(trove%dihedtype(:)==105).or.any(trove%dihedtype(:)==106).or.&
         any(trove%dihedtype(:)==107).or.any(trove%dihedtype(:)==108) ) then
       !
       do icoord = 1,trove%Ncoords
        !
        call diff_local2cartesian(icoord,a0,Bmat_t)
        !
        Bmat(icoord,:,:) = Bmat_t(:,:)
        !
      enddo
      !
      return 
      !
    endif
    !
    ! for easy reference 
    !
    Nmodes = trove%Nmodes
    Natoms = trove%Natoms
    Nbonds = trove%Nbonds
    Nangles = trove%Nangles
    Ndihedrals = trove%Ndihedrals
    !
    ! defining the delta function
    !
    deltaf = 0
    zetaf = 0
    !
    do ix=1,trove%Natoms
      !
      deltaf(ix,ix) = 1.0_ark
      !
    enddo
    !
    do ix=1,trove%Natoms
      do iy=1,trove%Natoms
        do iz=1,trove%Natoms
          zetaf(ix,iy,iz) = deltaf(ix,iy) - deltaf(ix,iz)
        enddo
      enddo
    enddo
    !
    ! Three cartesian axes
    !
    cart_vec(1,:)  = (/1.0_ark,0.0_ark,0.0_ark/)
    cart_vec(2,:)  = (/0.0_ark,1.0_ark,0.0_ark/)
    cart_vec(3,:)  = (/0.0_ark,0.0_ark,1.0_ark/)
    !
    ida = 0 ! check if the corresponding derivative d alpha(1,2,3)/d Xna has been defined
    !
    ! Defining connections for bonds and angles
    !
    rcon = 0
    !
    do ibond = 1,Nbonds
       n1 = trove%bonds(ibond,1)
       n2 = trove%bonds(ibond,2)
       rcon(n1,n2) = sqrt( sum( ( a0(n1,:)-a0(n2,:) )**2 ) )
       rcon(n2,n1) = rcon(n1,n2)
    enddo 
    !
    acon = 0
    do iangle = 1,Nangles
       n1 = trove%angles(iangle,1)
       n0 = trove%angles(iangle,2)
       n2 = trove%angles(iangle,3)
       !
       do ix =1,3
          a_t1(ix) = a0(n1,ix) - a0(n0,ix)
          a_t2(ix) = a0(n2,ix) - a0(n0,ix)
       enddo
       !
       cosa = sum(a_t1(:)*a_t2(:) )/( rcon(n1,n0)*rcon(n2,n0) )
       !
       acon(n1,n0,n2) = acos(cosa) 
       acon(n2,n0,n1) = acon(n1,n0,n2)
    enddo
    !
    ! Defining Bmat matrix
    !
    Bmat = 0
    !
    ! 1. stretching coordinates 
    !
    do ibond = 1,Nbonds
       n1 = trove%bonds(ibond,1)
       n2 = trove%bonds(ibond,2)
       Bmat(ibond,n1,:) = (a0(n1,:) - a0(n2,:))/rcon(n1,n2)
       Bmat(ibond,n2,:) = -Bmat(ibond,n1,:)
    enddo
    !
    ! 2. bending coordinates 
    !
    do iangle = 1,Nangles
         !
         n1 = trove%angles(iangle,1)
         n0 = trove%angles(iangle,2)
         n2 = trove%angles(iangle,3)
         r1 =  rcon(n1,n0)
         r2 =  rcon(n2,n0)
         cosa = cos( acon(n1,n0,n2) )
         sina = sin( acon(n1,n0,n2) )
        !
        ! First we define d cos(alpha) / d x_na
        !
        Bmat(Nbonds+iangle,n1,:) = ( a0(n2,:)-a0(n0,:) )/( r1*r2 ) -( a0(n1,:)-a0(n0,:) )*cosa/r1**2
        Bmat(Nbonds+iangle,n2,:) = ( a0(n1,:)-a0(n0,:) )/( r1*r2 ) -( a0(n2,:)-a0(n0,:) )*cosa/r2**2


        Bmat(Nbonds+iangle,n0,:) =-( a0(n1,:)-a0(n0,:) )/( r1*r2 )*( 1.0_ark-cosa*r2/r1 )  &
                                  -( a0(n2,:)-a0(n0,:) )/( r1*r2 )*( 1.0_ark-cosa*r1/r2 )

        !
        ! Second, we transform to the d alpha / d x_na = -1/sin(alpha) d cos(alpha) / d x_na
        !
        Bmat(Nbonds+iangle,n1,:) = -Bmat(Nbonds+iangle,n1,:)/sina
        Bmat(Nbonds+iangle,n0,:) = -Bmat(Nbonds+iangle,n0,:)/sina
        Bmat(Nbonds+iangle,n2,:) = -Bmat(Nbonds+iangle,n2,:)/sina
        !
        da_dxna(n1,n0,n2,:,:) = Bmat(Nbonds+iangle,:,:)
        da_dxna(n2,n0,n1,:,:) = Bmat(Nbonds+iangle,:,:)
        !
        ida(n1,n0,n2) = 1 ;  ida(n2,n0,n1) = 1
        !
    enddo
    !
    ! 3. Dihedral coordinates
    !
    do iangle = 1,Ndihedrals
       !
       J = trove%dihedtype(iangle)
       !
       !select case (abs(J)) 
       !
       select case (J) 
       !
       case(1) ! type 1 
          !
          n1 = trove%dihedrals(iangle,1)
          n4 = trove%dihedrals(iangle,2)
          n3 = trove%dihedrals(iangle,3)
          n2 = trove%dihedrals(iangle,4)
          !
          a_t1(:) = a0(n1,:) - a0(n4,:)
          a_t2(:) = a0(n2,:) - a0(n4,:)
          a_t3(:) = a0(n3,:) - a0(n4,:)
          !
          delta = 0 
          !
          do ix =1,3
            do iy =1,3
              do iz =1,3
                 delta = delta + epsil(ix,iy,iz)*a_t1(ix)*a_t2(iy)*a_t3(iz)
              enddo
            enddo
          enddo
          !
          delta = delta/( rcon(n1,n4)*rcon(n2,n4)*rcon(n3,n4) )
          !
          r1 =  rcon(n1,n4)
          r2 =  rcon(n2,n4)
          r3 =  rcon(n3,n4)
          !
          ddelta_t = 0 
          !
          do iy =1,3
             do iz =1,3
                ddelta_t(n1,:) = ddelta_t(n1,:) + epsil(:,iy,iz)*a_t2(iy)*a_t3(iz) 
                ddelta_t(n2,:) = ddelta_t(n2,:) + epsil(:,iy,iz)*a_t3(iy)*a_t1(iz) 
                ddelta_t(n3,:) = ddelta_t(n3,:) + epsil(:,iy,iz)*a_t1(iy)*a_t2(iz) 
             enddo
          enddo
          !
          Bmat(Nbonds+Nangles+iangle,n1,:) = ddelta_t(n1,:)/(r1*r2*r3)-a_t1(:)/r1**2*delta
          Bmat(Nbonds+Nangles+iangle,n2,:) = ddelta_t(n2,:)/(r1*r2*r3)-a_t2(:)/r2**2*delta
          Bmat(Nbonds+Nangles+iangle,n3,:) = ddelta_t(n3,:)/(r1*r2*r3)-a_t3(:)/r3**2*delta
          !
          Bmat(Nbonds+Nangles+iangle,n4,:) =-Bmat(Nbonds+Nangles+iangle,n1,:)   & 
                                            -Bmat(Nbonds+Nangles+iangle,n2,:)   & 
                                            -Bmat(Nbonds+Nangles+iangle,n3,:)

          !
          cosa1 = sum(a_t2(:)*a_t3(:) )/( r2*r3 )
          cosa2 = sum(a_t1(:)*a_t2(:) )/( r1*r2 )
          cosa3 = sum(a_t1(:)*a_t3(:) )/( r1*r3 )
          !
          sina1 = sqrt(1.0_ark-cosa1**2)
          sina2 = sqrt(1.0_ark-cosa2**2)
          sina3 = sqrt(1.0_ark-cosa3**2)
          !
          norm_2 = 3.0_ark-cosa3**2-cosa2**2-cosa1**2+2.0_ark*cosa3*cosa1-2.0_ark*cosa2+&
                   2.0_ark*cosa2*cosa3-2.0_ark*cosa1+2.0_ark*cosa2*cosa1-2.0_ark*cosa3
          !
          dnorm_da1 = 2.0_ark**sina1*( 1.0_ark+cosa1-cosa2-cosa3 )
          dnorm_da2 = 2.0_ark**sina2*( 1.0_ark+cosa2-cosa3-cosa1 )
          dnorm_da3 = 2.0_ark**sina3*( 1.0_ark+cosa3-cosa1-cosa2 )
          !
          if(norm_2<small_) then 
            !
            write (out,"('Bmat_generation: norm2 = ',f18.8', delta = infinity!')") norm_2
            stop 'Bmat_generation - bad norm2'
            !
          endif
          !
          if (any((/ida(n2,n4,n3),ida(n1,n4,n2),ida(n3,n4,n1)/)==0)) then 
            write (out,"('Bmat_gener: dihedral: not all deriv of alpha wrt Xna was defined: n0,n1,n2,n3 =  ',4i3,' ida =  ',3i4)") &
                   n0,n1,n2,n3,ida(n2,n0,n3),ida(n1,n0,n2),ida(n3,n0,n1)
            stop 'Bmat_generation - undefined angle derivatives'
          endif 
          !
          Bmat(Nbonds+Nangles+iangle,:,:) = Bmat(Nbonds+Nangles+iangle,:,:)/sqrt(norm_2)-&
              0.5_ark/sqrt(norm_2)**3*( &
              dnorm_da1*da_dxna(n2,n4,n3,:,:)+dnorm_da2*da_dxna(n1,n4,n2,:,:)+dnorm_da3*da_dxna(n3,n4,n1,:,:) )
          !
          !
       case(-2,2) ! type 2   B = (a*b)/(|a|*|b|), a = [y1 times y2]; b = [y2 times y3]
          !
          n1 = trove%dihedrals(iangle,1)
          n2 = trove%dihedrals(iangle,2)
          n3 = trove%dihedrals(iangle,3)
          n4 = trove%dihedrals(iangle,4)
          !
          if (J>0) then 
             !
             u(:) = a0(n4,:) - a0(n3,:)
             w(:) = a0(n2,:) - a0(n3,:)
             v(:) = a0(n1,:) - a0(n2,:)
             !
             m = n4 ; o = n3 ; p = n2 ; n = n1
             !
          else
             !
             u(:) = a0(n4,:) - a0(n3,:)
             w(:) = a0(n3,:) - a0(n2,:)
             v(:) = a0(n1,:) - a0(n2,:)
             !
             m = n4 ; o = n2 ; p = n3 ; n = n1
             !
          endif
          !
          mnop(1) = m; mnop(2) = o; mnop(3) = p; mnop(4) = n
          !
          r1 =  sqrt(sum(u(:)**2))
          r2 =  sqrt(sum(w(:)**2))
          r3 =  sqrt(sum(v(:)**2))
          !
          u = u/r1
          w = w/r2
          v = v/r3
          !
          cosu = sum( u(:)*w(:) )
          cosv =-sum( v(:)*w(:) )
          !
          sinu = sqrt(1.0_ark-cosu**2)
          sinv = sqrt(1.0_ark-cosv**2)
          !
          dvec1(:) = 0 
          dvec2(:) = 0 
          !
          do iy =1,3
            do iz =1,3
               dvec1(:) = dvec1(:) + epsil(:,iy,iz)*u(iy)*w(iz)
               dvec2(:) = dvec2(:) + epsil(:,iy,iz)*v(iy)*w(iz)
            enddo
          enddo
          !
          vec1 = sqrt( sum(dvec1(:)**2) )
          vec2 = sqrt( sum(dvec2(:)**2) )
          !
          B = sum( dvec1(:)*dvec2(:) )/( sinu*sinv )
          !
          tau_sign = 0
          !
          do iy =1,3
            do iz =1,3
               tau_sign = tau_sign + sum(epsil(:,iy,iz)*a_t2(:)*dvec1(iy)*dvec2(iz))
            enddo
          enddo
          !
          if (abs(tau_sign)<small_a) tau_sign = small_a
          tau_sign = sign(1.0_ark,tau_sign)
          !
          do i = 1,4
            !
            a = mnop(i)
            !
            b_t  = zetaf(a,m,o)*dvec1(:)/sinu**2/r1 + zetaf(a,p,n)*dvec2(:)/sinv**2/r3 + & 
                   zetaf(a,o,p)*( dvec1(:)*cosu/sinu**2/r2 + dvec2(:)*cosv/sinv**2/r2 )

            Bmat(Nbonds+Nangles+iangle,a,:) = b_t*tau_sign
            !
          enddo
          !
       case(202,-202)
          !
          call diff_local2cartesian(Nbonds+Nangles+iangle,a0,Bmat_t,fmod=2.0_ark*pi)
          !
          Bmat(Nbonds+Nangles+iangle,:,:) = Bmat_t(:,:)
          !
       case(402,-402)
          !
          call diff_local2cartesian(Nbonds+Nangles+iangle,a0,Bmat_t,fmod=2.0_ark*pi)
          !
          Bmat(Nbonds+Nangles+iangle,:,:) = Bmat_t(:,:)
          !
       case(-302,302) ! type 2   B = (a*b)/(|a|*|b|), a = [y1 times y2]; b = [y2 times y3]
          !
          n1 = trove%dihedrals(iangle,1)
          n2 = trove%dihedrals(iangle,2)
          n3 = trove%dihedrals(iangle,3)
          n4 = trove%dihedrals(iangle,4)
          !
          if (J>0) then 
             !
             u(:) = a0(n4,:) - a0(n3,:)
             w(:) = a0(n2,:) - a0(n3,:)
             v(:) = a0(n1,:) - a0(n2,:)
             !
             m = n4 ; o = n3 ; p = n2 ; n = n1
             !
          else
             !
             u(:) = a0(n4,:) - a0(n3,:)
             w(:) = a0(n3,:) - a0(n2,:)
             v(:) = a0(n1,:) - a0(n2,:)
             !
             m = n4 ; o = n2 ; p = n3 ; n = n1
             !
          endif
          !
          mnop(1) = m; mnop(2) = o; mnop(3) = p; mnop(4) = n
          !
          r1 =  sqrt(sum(u(:)**2))
          r2 =  sqrt(sum(w(:)**2))
          r3 =  sqrt(sum(v(:)**2))
          !
          u = u/r1
          w = w/r2
          v = v/r3
          !
          cosu = sum( u(:)*w(:) )
          cosv =-sum( v(:)*w(:) )
          !
          sinu = sqrt(1.0_ark-cosu**2)
          sinv = sqrt(1.0_ark-cosv**2)
          !
          dvec1(:) = 0 
          dvec2(:) = 0 
          !
          do iy =1,3
            do iz =1,3
               dvec1(:) = dvec1(:) + epsil(:,iy,iz)*u(iy)*w(iz)
               dvec2(:) = dvec2(:) + epsil(:,iy,iz)*v(iy)*w(iz)
            enddo
          enddo
          !
          vec1 = sqrt( sum(dvec1(:)**2) )
          vec2 = sqrt( sum(dvec2(:)**2) )
          !
          B = sum( dvec1(:)*dvec2(:) )/( sinu*sinv )
          !
          tau_sign = 0
          !
          do iy =1,3
            do iz =1,3
               tau_sign = tau_sign + sum(epsil(:,iy,iz)*a_t2(:)*dvec1(iy)*dvec2(iz))
            enddo
          enddo
          !
          if (abs(tau_sign)<100.0_ark*small_a) tau_sign = B
          tau_sign = sign(1.0_ark,tau_sign)
          !
          do i = 1,4
            !
            a = mnop(i)
            !
            b_t  = zetaf(a,m,o)*dvec1(:)/sinu**2/r1 + zetaf(a,p,n)*dvec2(:)/sinv**2/r3 + & 
                   zetaf(a,o,p)*( dvec1(:)*cosu/sinu**2/r2 + dvec2(:)*cosv/sinv**2/r2 )

            Bmat(Nbonds+Nangles+iangle,a,:) = b_t*tau_sign
            !
          enddo
          !
       case(101) ! The special bond-angles for the linear molecule case 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = a0(n1,:) - a0(n0,:)
          a_t2(:) = a0(n2,:) - a0(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          !
          u =  a_t1(:)/r1
          v =  a_t2(:)/r2
          !
          w(:) = MLvector_product(u,v)
          !
          ! special angle is -arcsin( kappa . [uxv] ), 
          ! i.e. the scalar product kappa.w is a kappa component of w
          !
          sindelta = w(kappa)
          !
          if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('FLfromcartesian2local: sindelta>1: ',f18.8)") sindelta
             write (out,"('Consider change difftype ')")
             stop 'FLfromcartesian2local - bad sindelta'
             !
          elseif ( abs(sindelta)>=1.0_ark) then 
             !
             phi = 0.0_ark
             !
          else
             ! 
             phi =  -asin(sindelta)
             !
          endif
          !
          icoord = Nbonds+Nangles+iangle
          !
          !Bmat(icoord,n1,:) =-( e2(:) + sin(phi)*e1(:))/(r1*cos(phi))
          !Bmat(icoord,n2,:) = ( e1(:) - sin(phi)*e2(:))/(r2*cos(phi))
          !Bmat(icoord,n0,:) =-( Bmat(icoord,n1,:) + Bmat(icoord,n2,:) )
          !
          Bmat(icoord,n1,:) = ( e2(:) - sin(phi)*u(:))/(r1*cos(phi))
          Bmat(icoord,n2,:) =-( e1(:) + sin(phi)*v(:))/(r2*cos(phi))
          Bmat(icoord,n0,:) =-( Bmat(icoord,n1,:) + Bmat(icoord,n2,:) )
          !
       case(102) ! The special bond-angle for the linear molecule case confined to the plane
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = a0(n1,:) - a0(n0,:)
          a_t2(:) = a0(n2,:) - a0(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          !
          u =  a_t1(:)/r1
          v =  a_t2(:)/r2
          !
          w(:) = MLvector_product(u,v)
          !
          ! special angle is -arcsin( kappa . [uxv] ), 
          ! i.e. the scalar product kappa.w is a kappa component of w
          !
          sindelta = w(kappa)
          !
          !if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
          !   !
          !   write (out,"('FLfromcartesian2local: sindelta>1: ',f18.8)") sindelta
          !   write (out,"('Consider change difftype ')")
          !   stop 'FLfromcartesian2local - bad sindelta'
          !   !
          !elseif ( abs(sindelta)>=1.0_ark) then 
          !   !
          !   phi = 0.0_ark
          !   !
          !else 
          !   phi = -asin(sindelta)
          !   !
          !endif
          !
          e1(:) = -MLvector_product(cart_vec(kappa,:),u)
          e2(:) = -MLvector_product(cart_vec(kappa,:),v)
          !
          icoord = Nbonds+Nangles+iangle
          !
          !Bmat(icoord,n1,:) = ( e2(:) - sin(phi)*u(:))/(r1*cos(phi))
          !Bmat(icoord,n2,:) =-( e1(:) + sin(phi)*v(:))/(r2*cos(phi))
          !Bmat(icoord,n0,:) =-( Bmat(icoord,n1,:) + Bmat(icoord,n2,:) )
          !
          Bmat(icoord,n1,:) = e2(:)/r1 - sindelta*u(:)/r1
          Bmat(icoord,n2,:) =-e1(:)/r2 + sindelta*v(:)/r2
          !
          Bmat(icoord,n0,:) =-( Bmat(icoord,n1,:) + Bmat(icoord,n2,:) )
          !
          call diff_local2cartesian(icoord,a0,Bmat_t)
          !
          Bmat(icoord,:,:) = Bmat_t(:,:)
          !
       case(103) ! The special bond-angles for the linear molecule case 
          !
          ! This derivative is now done numerically, the rest can be skipped 
          !
          !cycle 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = a0(n1,:) - a0(n0,:)
          a_t2(:) = a0(n2,:) - a0(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          !
          u =  a_t1(:)/r1
          !
          phi = u(kappa)
          !
          !e1(:) = MLvector_product(cart_vec(kappa,:),u)
          !e2(:) = MLvector_product(cart_vec(kappa,:),v)
          !
          icoord = Nbonds+Nangles+iangle
          !
          !Bmat(icoord,n1,:) = deltaf(kappa,:)/r1 - phi/r1**3*a_t1(:)
          !Bmat(icoord,n2,:) =-Bmat(icoord,n1,:)
          !
          !Bmat(icoord,n1,:) = ( deltaf(kappa,:) - u(kappa)*u(:) )/r1
          !Bmat(icoord,n0,:) = -Bmat(icoord,n1,:)
          !
          call diff_local2cartesian(icoord,a0,Bmat_t)
          !
          Bmat(icoord,:,:) = Bmat_t(:,:)
          !
          do ibond = 1,Nbonds
             !
             k1 = trove%bonds(ibond,1)
             k2 = trove%bonds(ibond,2)
             !
             if (k1/=n1.and.k2/=n1) cycle
             !
             call   diff_local2cartesian(ibond,a0,Bmat_t)
             Bmat(ibond,:,:) = Bmat_t(:,:)
             !
          enddo
          !
       case(-5) ! type -2   B = (a*b)/(|a|*|b|), a = [y1 times y2]; b = [y2 times y3]
          !
          n1 = trove%dihedrals(iangle,1)
          n2 = trove%dihedrals(iangle,2)
          n3 = trove%dihedrals(iangle,3)
          n4 = trove%dihedrals(iangle,4)
          !
          if (J<0) then 
             !
             a_t1(:) = a0(n4,:) - a0(n3,:)
             a_t2(:) = a0(n2,:) - a0(n3,:)
             a_t3(:) = a0(n2,:) - a0(n1,:)
             !
          else
             !
             a_t1(:) = a0(n4,:) - a0(n3,:)
             a_t2(:) = a0(n2,:) - a0(n3,:)
             a_t3(:) = a0(n2,:) - a0(n1,:)
             !
          endif 
          !
          a1 = sqrt(sum(a_t1(:)**2))
          a2 = sqrt(sum(a_t2(:)**2))
          a3 = sqrt(sum(a_t3(:)**2))
          !
          dvec1(:) = 0 
          dvec2(:) = 0 
          !
          do iy =1,3
            do iz =1,3
               dvec1(:) = dvec1(:) + epsil(:,iy,iz)*a_t1(iy)*a_t2(iz)
               dvec2(:) = dvec2(:) + epsil(:,iy,iz)*a_t2(iy)*a_t3(iz)
            enddo
          enddo
          !
          !if (J<0) dvec2 = -dvec2
          !
          vec1 = sqrt( sum(dvec1(:)**2) )
          vec2 = sqrt( sum(dvec2(:)**2) )
          !
          B = sum( dvec1(:)*dvec2(:) )/( vec1*vec2 )
          !
          tau_sign = 0
          !
          do iy =1,3
            do iz =1,3
               tau_sign = tau_sign - sum(epsil(:,iy,iz)*a_t2(:)*dvec1(iy)*dvec2(iz))
            enddo
          enddo
          !
          if (abs(tau_sign)<small_a) tau_sign = small_a
          tau_sign = sign(1.0_ark,tau_sign)
          !
          Bmat_t = 0 
          !
          Bmat_t(3,:) =  a2/vec1**2*dvec1(:)
          Bmat_t(4,:) = -a2/vec2**2*dvec2(:)
          Bmat_t(1,:) = ( sum( a_t3(:)*a_t2(:) )/a2**2-1.0_ark )*Bmat_t(4,:) &
                         -sum( a_t1(:)*a_t2(:) )/a2**2*Bmat_t(3,:)
          Bmat_t(2,:) = ( sum( a_t1(:)*a_t2(:) )/a2**2-1.0_ark )*Bmat_t(3,:) &
                         -sum( a_t3(:)*a_t2(:) )/a2**2*Bmat_t(4,:)
          !
          Bmat(Nbonds+Nangles+iangle,n1,:) = Bmat_t(4,:)*tau_sign
          Bmat(Nbonds+Nangles+iangle,n2,:) = Bmat_t(1,:)*tau_sign
          Bmat(Nbonds+Nangles+iangle,n3,:) = Bmat_t(2,:)*tau_sign
          Bmat(Nbonds+Nangles+iangle,n4,:) = Bmat_t(3,:)*tau_sign
          !
       end select 
       !
    enddo
    !
    if (job%verbose>=7) write(out,"('Bmat_generation/stop')") 
    !
    contains 

    subroutine diff_local2cartesian(icoord,a0,Bmat,fmod)
 
      integer(ik),intent(in):: icoord
      real(ark),intent(in)  ::  a0(trove%Natoms,3)
      real(ark),intent(out) ::  Bmat(trove%Natoms,3) 
      real(ark),intent(in),optional  ::  fmod ! a fmod-value which be subtructed from a large angle value to prevent discontinuity  
                                              ! of angles in finite differences 

      real(ark)             :: xi_p(trove%Ncoords),xi_m(trove%Ncoords),deltax,xna(trove%Natoms,3)
      real(ark)             :: xi_pp(trove%Ncoords),xi_mm(trove%Ncoords),dx,ddx
      integer(ik)           :: iatom,ix

          !
          xna = a0
          !
          deltax = fd_step_Bmat
          !
          do iatom = 1,trove%Natoms
             do ix = 1,3
                !
                xna(iatom,ix)  = a0(iatom,ix) + deltax
                call FLfromcartesian2local(xna,xi_p)
                !
                xna(iatom,ix)  = a0(iatom,ix) - deltax
                call FLfromcartesian2local(xna,xi_m)
                !
                xna(iatom,ix)  = a0(iatom,ix) + deltax*2.0_ark
                call FLfromcartesian2local(xna,xi_pp)
                !
                xna(iatom,ix)  = a0(iatom,ix) - deltax*2.0_ark
                call FLfromcartesian2local(xna,xi_mm)
                !
                dx  = xi_p(icoord) -xi_m(icoord)
                ddx = xi_pp(icoord)-xi_mm(icoord)
                !
                if (present(fmod)) then 
                  if (dx> fmod*0.5_ark) dx  = dx-fmod
                  if (dx<-fmod*0.5_ark) dx  = dx+fmod
                  if (ddx> fmod*0.5_ark) ddx  = ddx-fmod
                  if (ddx<-fmod*0.5_ark) ddx  = ddx+fmod
                endif 
                !
                !f = (-f4/12.0_ark+2.0_ark/3.0_ark*f2 & 
                !     +f3/12.0_ark-2.0_ark/3.0_ark*f1 )/h
                !
                !Bmat(iatom,ix) = 0.5_ark*(xi_p(icoord)-xi_m(icoord))/deltax
                !
                Bmat(iatom,ix) = (  ( -ddx )/12.0_ark &
                                      +2.0_ark/3.0_ark*( dx ) )/deltax
                !
                xna(iatom,ix)  = a0(iatom,ix)
                !
             enddo
          enddo    

    end subroutine diff_local2cartesian
    !
 end subroutine Bmat_generation

! Here we calculate the vibrational part of the jacobian: 
! s_vib vector = (d xi_l / d r_Na ) 
! defined at every point rho
! We use the iterative-recursive procedure from MP2005
! This version is for transformation to the normal coordinates 
!
  subroutine s_vib_s_rot_polynom1d(s_vib,s_rot)

    type(FLpolynomT),intent(out) :: s_vib(trove%Nmodes,trove%Natoms,3)
    type(FLpolynomT),intent(out) :: s_rot(3,trove%Natoms,3)

    real(ark):: chi_eq(trove%Nmodes)
    real(ark),allocatable     :: dr_na_dq(:,:,:,:),r_na_xi(:,:,:)

    type(FLpolynomT),pointer :: fl

    integer :: Nordersmax,iNcoeff2,int,alloc,alloc_p,n1,iNcoeff1

    integer(ik) ::  q1,q2,x1,k0
    integer(ik) ::  dimen,ieq
    integer(ik) ::  k(trove%Nmodes)
    integer(ik) ::  n,nmodes,Natoms
    integer(ik) ::  imode,jmode,lincoord,Nequat,Npoints,irho
    !
    real(ark),allocatable  ::  s_mat_t(:,:,:,:)

    real(ark)  :: xistep(trove%Nmodes)
    real(ark) :: dr_t(trove%Natoms*3),dr_tt(trove%Natoms*3)   !   d r/ d rho and  d2 r/ d rho^2
    real(ark)  :: xi_eq(1:trove%Nmodes) 
    real(ark) :: chi_t1,astep(2),axi_eq(1:trove%Nmodes),dxi,factor
    real(ark)                     :: step(2,trove%Nmodes),rhostep
    character(len=cl)            :: job_is,fname


    if (job%verbose>=4) write(out,"(/'s_vib_s_rot_polynom1d/start  ')") 
    !
    call TimerStart('s_vib_s_rot')
    !
    ! Sumstitution for easier reference within the procedure 
    !
    Nmodes = trove%Nmodes
    Natoms = trove%Natoms
    Npoints = trove%Npoints
    lincoord = trove%lincoord
    !
    ! define number of equation that will be solved 
    !
    Nequat = 6+Nmodes-min(1,lincoord)
    !
    iNcoeff1= trove%RangeOrder(trove%NKinOrder+1) 
    iNcoeff2= trove%RangeOrder(trove%NKinOrder+2) 

    !k = 0 ; k(1) = trove%NKinOrder+2+1
    !iNcoeff3 = FLQindex(trove%Nmodes_e,k)
    !
    Nordersmax = trove%NKinOrder+2

    if ( Nordersmax/=s_vib(1,1,1)%orders ) then 
       write(out,"(' s_vib_s_rot_polynom1d: sizes of s_vib and r_na are different: ',2i9)") s_vib(1,1,1)%orders,Nordersmax
       stop 's_vib_s_rot_polynom1d: sizes of s_vib and r_na are different'
    endif
    !
    !
    ! Here we finish the introduction and start the actual calculations 
    !
    rhostep = trove%rhostep
    !
    !
    if (trove%internal_coords=='LOCAL') then 
      !
      write(out,"(' s_vib_s_rot_polynom1d: LOCAL-option has been deactivated, use LINEAR')") 
      stop 's_vib_s_rot_polynom1d: LOCAL option has been deactivated'
      !
    endif
    ! 
    chi_eq(:) = trove%chi_ref(:,0)
    xi_eq = 0 
    axi_eq= 0
    do q1 = 1,trove%Nmodes
      !
      xi_eq(q1) = MLcoord_direct(chi_eq(q1),1,q1)
      chi_t1 = chi_eq(q1)+trove%fdstep(q1)
      xistep(q1) = MLcoord_direct(chi_t1,1,q1)-xi_eq(q1)
      !
    enddo
    !
    do imode=1,trove%Nmodes
       call ML_check_steps4coordinvert(axi_eq,1,imode,astep)
       step(:,imode) = xistep(imode)*astep(:)
    enddo
    !
    if (job%verbose>=2) then
      write(out,"('generated step-left :',40f14.6)") step(1,1:min(40,trove%Nmodes))
      write(out,"('generated step-right:',40f14.6)") step(2,1:min(40,trove%Nmodes))
    endif
    !
    ! We will need Tmat matryix and bm vector for the T x =b linear equation 
    ! Tmat is a M by M matrix, where M is a number of equations and variables. 
    ! Since all equations can  be solved for every mode independently, 
    ! the size of the problem dimen x dimen, where  dimen = 3*Natoms
    !
    dimen = Nequat
    !
    ! The linear equation will be solve with Lapack at double precision 
    ! therefor we will need a double precision matrix A and vector b 
    ! as analogs of Tmat and Bm, respectively 
    !
    if (job%verbose>=2) write(out,"(/'Calculating the s-matrix...')") 
    !
    ! it is gonna be a big loop over all Npoints+1
    !
    !$omp parallel private(s_mat_t,dr_na_dq,r_na_xi,alloc_p)
    allocate (s_mat_t(max(Nmodes,3),Natoms,3,iNcoeff2),&
              dr_na_dq(Natoms,3,Nmodes,iNcoeff2),r_na_xi(Natoms,3,iNcoeff2),stat=alloc_p)
    if (alloc_p/=0) then
        write (out,"(' Error ',i9,' trying to allocate arrays s_mat_t,r_na_xi,dr_na_dq')") alloc_p
        stop 's_vib_s_rot_polynom1d, s_mat_t - out of memory'
    end if
    !
    r_na_xi = 0 
    dr_na_dq = 0 
    !
    !$omp do private(irho,k0,imode,jmode,n1,x1,ieq,k,factor,job_is,chi_eq,xi_eq,dr_t,dr_tt) schedule(guided)
    do irho = 0,Npoints
      !
      if (job%verbose>=2.and.mod(irho,max(Npoints/20,1))==0) write(out,"('irho= ',i5)") irho
      ! 
      ! Defining r_na expansion with respect to the linearized (normal) coordinates 
      !
      ! The zero order part 
      !
      r_na_xi(:,:,1) = trove%b0(:,:,irho)
      !
      !FLirho = irho
      !
      do k0 = 2,iNcoeff2
         !
         k(:) = FLIndexQ(:,k0)
         !
         ! Factorial factors to convert derivatives to the force constants 
         !
         factor = 1.0_ark
         do imode = 1,Nmodes
           do jmode = 1,k(imode)
              factor = factor*real(jmode,kind=ark)
           enddo
         enddo 
         !
         ! rho point as an expansion center
         !
         chi_eq(:) = trove%chi_ref(:,irho)
         !
         do jmode = 1,Nmodes
           !
           xi_eq(jmode) = MLcoord_direct(chi_eq(jmode),1,jmode)
           !
         enddo
         !
         ! Here we first calculate the expansion of r_na in terms of xi 
         ! and then the same for d r_na/ d rho
         !
         job_is = 'xi2cartesian'
         !
         call FLfinitediffs_vect(job_is,Nmodes,trove%Natoms*3,k,xi_eq,step,irho,dr_t)
         !
         if (manifold/=0) then 
            !
            job_is = 'xi2dcartesian_drho'
            !
            call FLfinitediffs_vect(job_is,Nmodes,trove%Natoms*3,k,xi_eq,step,irho,dr_tt)
            !
         endif 
         !
         ieq = 0
         do n1 = 1,Natoms
            do x1 = 1,3
              !
              ieq = ieq +1 
              r_na_xi(n1,x1,k0) = dr_t(ieq)/factor
              !
              if (verbose>=6)  write(out,"('r_na_xi (',3i5,')',g20.12)") n1,x1,k0,r_na_xi(n1,x1,k0)
              !
              if (manifold/=0) then 
                !
                dr_na_dq(n1,x1,trove%Nmodes_n,k0) = dr_tt(ieq)/factor
                !
                if (verbose>=6)  write(out,"('dr_na_dq (',3i5,')',g20.12)") n1,x1,k0,dr_na_dq(n1,x1,trove%Nmodes_n,k0)
                !
              endif 
              !
            enddo
         enddo
         !
      enddo
      ! 
      ! Defining d r_na wrt q expansion in terms of the linearized (normal) coordinates 
      !
      do jmode = 1,Nmodes
         !
         if (jmode==trove%Nmodes_n) then 
            !
            dr_na_dq(:,:,jmode,1) = trove%db0(:,:,irho,1)
            !
         else
            !
            dr_na_dq(:,:,jmode,1) = trove%Amatrho(:,:,jmode,irho)
            !
         endif
         !
      enddo
      !
      ! 
      ! Since the rotational and vibrational S-matrix are been obtained using very similar procedures,
      ! we combine into the same procedure: 
      ! icase = 1 ->  vibrational, icase =2 -> rotational 
      ! the major difference is in the summation limits: either q1 = 1..Nmodes or x,y,z
      !
      call s_vib_s_rot_solve(1,Nmodes,trove%NKinOrder+2,iNcoeff2,trove%b0(:,:,irho),r_na_xi(:,:,:),&
                             dr_na_dq(:,:,:,:),s_mat_t(1:Nmodes,:,:,:))
      !
      do jmode = 1,Nmodes
         do n1 = 1,Natoms
            do x1 = 1,3
              s_vib(jmode,n1,x1)%field(1:iNcoeff2,irho) = s_mat_t(jmode,n1,x1,1:iNcoeff2)
            enddo 
         enddo
      enddo
      !
      call s_vib_s_rot_solve(2,     3,trove%NKinOrder+1,iNcoeff1,trove%b0(:,:,irho), r_na_xi(:,:,:),&
                            dr_na_dq(:,:,:,:),s_mat_t(1:3,:,:,:))
      !
      do jmode = 1,3
         do n1 = 1,Natoms
            do x1 = 1,3
              s_rot(jmode,n1,x1)%field(1:iNcoeff1,irho) = s_mat_t(jmode,n1,x1,1:iNcoeff1)
            enddo 
         enddo
      enddo
      !
    enddo
    !$omp end do
    !
    !write(out,"(/'vib')") 
    !do irho = 0,Npoints
    !  if (job%verbose>=4) then 
    !    write(out,"(i8,18(<Nmodes>f16.8))") irho,((s_vib(Nmodes,n1,x1)%field(1,irho),n1=1,Natoms),x1=1,3)
    !  endif
    !enddo
    !
    deallocate(s_mat_t,stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to deallocate arrays s_mat_t')") alloc
        stop 's_vib_s_rot_polynom1d, s_mat_t - out of memory'
    end if
    deallocate(r_na_xi,stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to deallocate arrays r_na_xi')") alloc
        stop 's_vib_s_rot_polynom1d, s_mat_t - out of memory'
    end if
    deallocate(dr_na_dq,stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to deallocate array dr_na_dq')") alloc
        stop 's_vib_s_rot_polynom1d, s_mat_t - out of memory'
    end if
    !$omp end parallel
    !
    if (job%verbose>=2) write(out,"('...done!')") 
    !
    !call ArrayStop('dr_na_dq-r_na_xi')
    !
    ! clean up! 
    !
    call TimerStop('s_vib_s_rot')
    !
    if (job%verbose>=4) write(out,"('s_vib_s_rot_polynom1d/end  ')") 
    !
    !

  end subroutine s_vib_s_rot_polynom1d



! Here we calculate the vibrational part of the jacobian: 
! s_vib vector = (d xi_l / d r_Na ) 
! defined at every point rho
! We use the iterative-recursive procedure from MP2005
! This version is for transformation to the normal coordinates 
!
  subroutine s_vib_s_rot_Sorensen(s_vib,s_rot)

    type(FLpolynomT),intent(inout) :: s_vib(trove%Nmodes,trove%Natoms,3)
    type(FLpolynomT),intent(inout) :: s_rot(3,trove%Natoms,3)
    !
    integer :: alloc,alloc_p,n1,iNcoeff1,iNcoeff2
    !
    integer(ik) ::  x1,x2,k0,k1,Ng
    integer(ik) ::  n,nmodes,Natoms,iatom,ix,jx,kx
    integer(ik) ::  imode,jmode,kmode,mmode,Npoints,irho,lincoord,Tcoeff,iswitch
    !
    real(ark)    :: Inertm(3),r_t(0:trove%Npoints),f_t,df_t,factor
    !
    real(ark)   :: z_t,rho_switch
    character(len=cl)            :: job_is
    !
    real(ark),allocatable     ::  s_mat_t(:,:,:,:),eta(:,:,:),tmat(:,:,:,:),Jmat(:,:,:),c(:,:,:),dzeta(:,:,:),dzeta_(:,:,:)
    integer(ik),allocatable  ::  powers(:)
    !
    if (job%verbose>=4) write(out,"(/'s_vib_s_rot_Sorensen/start  ')") 
    !
    call TimerStart('s_vib_s_rot')
    !
    ! Sumstitution for easier reference within the procedure 
    !
    Nmodes = trove%Nmodes
    Natoms = trove%Natoms
    Npoints = trove%Npoints
    !
    iNcoeff1= trove%RangeOrder(trove%NKinOrder+1) 
    iNcoeff2= trove%RangeOrder(trove%NKinOrder+2)
    !
    ! to convert L2 to cm-1
    !
    factor = real(planck,ark)*real(avogno,ark)*real(1.0d+16,kind=ark)/(4.0_ark*pi*pi*real(vellgt,ark))
    ! 
    !factor = real(planck,ark)/(4.0_ark*pi*pi*real(vellgt,ark))
    !
    factor = 1.0_ark
    !
    ! moments of inertia
    !
    allocate(trove%imat_s(0:Npoints))
    !
    !chi_eq(:) = trove%chi_ref(:,0)
    !xi_eq = 0 
    !axi_eq= 0
    !do q1 = 1,trove%Nmodes
    !  !
    !  xi_eq(q1) = MLcoord_direct(chi_eq(q1),1,q1)
    !  chi_t1 = chi_eq(q1)+trove%fdstep(q1)
    !  xistep(q1) = MLcoord_direct(chi_t1,1,q1)-xi_eq(q1)
    !  !
    !enddo
    !
    !do imode=1,trove%Nmodes
    !   call ML_check_steps4coordinvert(axi_eq,1,imode,astep)
    !   step(:,imode) = xistep(imode)*astep(:)
    !enddo
    !
    if ( trove%NKinOrder+2/=s_vib(1,1,1)%orders ) then 
       write(out,"(' s_vib_s_rot_Sorensen: sizes of s_vib and r_na are different: ',2i9)") s_vib(1,1,1)%orders,trove%NKinOrder+2
       stop 's_vib_s_rot_Sorensen: sizes of s_vib and r_na are different'
    endif 
    !
    if (trove%internal_coords=='LOCAL') then 
      !
      write(out,"(' s_vib_s_rot_Sorensen: LOCAL-option has been deactivated, sorry ')") 
      stop 's_vib_s_rot_Sorensen: LOCAL option has been deactivated'
      !
    endif
    ! 
    if (job%verbose>=2) write(out,"(/'Calculating the s-matrix...')") 
    !
    ! it is gonna be a big loop over all Npoints+1
    !
    do n1 = 1,Natoms
       do x1 = 1,3
         do imode = 1,trove%Nmodes
           s_vib(imode,n1,x1)%field(:,:) = 0 
         enddo
         do imode = 1,3
           s_rot(imode,n1,x1)%field(:,:) = 0 
         enddo
      enddo
    enddo
    !
    !lincoord=0
    lincoord = trove%lincoord
    !
    !$omp parallel private(s_mat_t,eta,tmat,Jmat,c,dzeta,dzeta_,powers,alloc_p)
    allocate (s_mat_t(4,Natoms,3,iNcoeff2),eta(4,4,iNcoeff2),tmat(trove%Natoms,3,4,0:trove%Nmodes),&
              Jmat(4,4,0:trove%Nmodes),c(4,trove%Natoms,3),dzeta(4,trove%Nmodes,trove%Nmodes),&
              dzeta_(3,trove%Nmodes,trove%Nmodes),powers(trove%Nmodes),stat=alloc_p)
    if (alloc_p/=0) then
        write (out,"(' Error ',i9,' trying to allocate arrays s_mat_t')") alloc_p
        stop 's_vib_s_rot_Sorensen, s_mat_t - out of memory'
    end if
    !
    rho_switch = trove%specparam(Nmodes)
    !
    iswitch = 0 
    if (Npoints>0) then 
      iswitch = mod(nint( ( rho_switch-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
    endif
    !
    !$omp do private(irho,lincoord,Inertm,x1,x2,Ng,ix,n1,imode,z_t,k0,jmode,jx,kx,k1,kmode,mmode) schedule(guided)
    do irho = 0,Npoints
      !
      lincoord = trove%lincoord
      if (Nmodes == 3*natoms-6.and.irho>iswitch) lincoord = 0
      !
      if (job%verbose>=2.and.mod(irho,max(Npoints/10,1))==0) write(out,"('irho= ',i5)") irho
      !
      ! Inertia moments, zero order (equilibrium)
      !
      Inertm(1) = sum( trove%mass(:)*( trove%b0(:,2,irho)**2+ trove%b0(:,3,irho)**2 ) )
      Inertm(2) = sum( trove%mass(:)*( trove%b0(:,1,irho)**2+ trove%b0(:,3,irho)**2 ) )
      Inertm(3) = sum( trove%mass(:)*( trove%b0(:,1,irho)**2+ trove%b0(:,2,irho)**2 ) )
      !
      trove%imat_s(irho) = Inertm(1)
      !
      !do x1 = 1,3
      !    if (Inertm(x1)<5e-2) lincoord=x1
      !enddo
      !
      !if (irho<15) lincoord = trove%lincoord
      !
      tmat = 0
      !
      ! Check if all Inertia moments are non zero 
      !
      do x1 = 1,3
        if (Inertm(x1)<sqrt(small_).and.lincoord/=x1) then 
           trove%sing_at_rho_0 = .true.
        endif
        !
        do x2 = x1+1,x1-1
          !
          if (sum( trove%mass(:)*trove%b0(:,x1,irho)*trove%b0(:,x2,irho) )>sqrt(small_)) then 
             write(out,"('s_rot_solve: illegal use of Sorensen, non diagonal Inertia moment non zero:',2i4,f20.10)") &
                   x1,x2,trove%mass(:)*trove%b0(:,x1,irho)*trove%b0(:,x2,irho)
             stop 'non zero a non diagonal Inertia moment'
          endif
          !
        enddo
      enddo
      !
      ! Maximal number of the linear terms 
      !
      Ng = 0
      !
      do ix = 1,3
         !
         if (ix/=lincoord) then 
           !
           Ng = Ng + 1
           !
           do n1 = 1,trove%Natoms
              do x1 = 1,3
                 !
                 tmat(n1,x1,Ng,0) = sum(epsil(x1,ix,:)*trove%b0(n1,:,irho))
                 c(Ng,n1,x1) = trove%mass(n1)*tmat(n1,x1,Ng,0)
                 !
                 do imode = 1,trove%Nmodes_e
                   !
                   tmat(n1,x1,Ng,imode) = sum(epsil(x1,ix,:)*trove%Amatrho(n1,:,imode,irho))
                   !
                 enddo
                 !
              enddo
           enddo
           !
         endif
         !
      enddo
      !
      if (manifold/=0) then 
         !
         Ng = Ng + 1
         !
         do n1 = 1,trove%Natoms
            do x1 = 1,3
               !
               tmat(n1,x1,Ng,0) = trove%db0(n1,x1,irho,1)
               c(Ng,n1,x1) = trove%mass(n1)*tmat(n1,x1,Ng,0)
               !
               do imode = 1,trove%Nmodes_e
                 !
                 tmat(n1,x1,Ng,imode) = trove%dAmatrho(n1,x1,imode,irho,1)
                 !
               enddo
               !
            enddo
         enddo
         !
      endif
      !
      ! We must control the relations involving the translational t-vectors
      !
      do x1 = 1,Ng
         !
         z_t = sum(c(x1,:,1)+c(x1,:,2)+c(x1,:,3))
         !
         if (abs(z_t)>sqrt(small_)) then 
             write(out,"('s_rot_solve: the translational relation for c is not zero for g = ',i4,': ',f20.10)") x1,z_t
             stop 'the translational relation for c is not zero'
         endif
         !
      enddo
      !
      ! First we define the  J-matrix (see Sorensen)
      !
      Jmat = 0
      !
      do x1 = 1,Ng
         !
         do x2 = 1,Ng
            !
            Jmat(x1,x2,0) = sum( trove%mass(:)*( tmat(:,1,x1,0)*tmat(:,1,x2,0) &
                                                +tmat(:,2,x1,0)*tmat(:,2,x2,0) &
                                                +tmat(:,3,x1,0)*tmat(:,3,x2,0) ) )
            !
            do imode = 1,trove%Nmodes_e
               !
               Jmat(x1,x2,imode) =  sum( trove%mass(:)*( tmat(:,1,x1,0)*tmat(:,1,x2,imode) + &
                                                         tmat(:,2,x1,0)*tmat(:,2,x2,imode) + &
                                                         tmat(:,3,x1,0)*tmat(:,3,x2,imode) ) )
               !
            enddo
            !
         enddo
         !
      enddo
      !
      !if (lincoord/=0) then 
      !  !
      !  Jmat(lincoord,lincoord,0) = Jmat(lincoord,lincoord,1)
      !  !
      !  Jmat(lincoord,lincoord,1) = 0
      !  !
      !endif 
      !
      call s_rot_solve(Ng,trove%NKinOrder+2,iNcoeff2,jmat(1:Ng,1:Ng,0:trove%Nmodes),eta(1:Ng,1:Ng,1:iNcoeff2))
      !
      do n1 = 1,Natoms
         do x2 = 1,3
            do k0 = 1,iNcoeff2
               do x1 = 1,Ng
                  s_mat_t(x1,n1,x2,k0) = sum( eta(x1,1:Ng,k0)*c(1:Ng,n1,x2) )
               enddo
               !
            enddo
         enddo 
         !
      enddo
      !
      do n1 = 1,Natoms
         !
         x1 = 0
         do ix = 1,3
           if (ix/=lincoord) then 
              !
              x1 = x1 + 1
              do x2 = 1,3
                s_rot(ix,n1,x2)%field(1:iNcoeff1,irho) = s_mat_t(x1,n1,x2,1:iNcoeff1)
              enddo
              !
           endif
         enddo 
         !
      enddo
      !
      if (manifold/=0) then 
         do n1 = 1,Natoms
            do x2 = 1,3
               !
               s_vib(trove%Nmodes_n,n1,x2)%field(:,irho) =  s_mat_t(Ng,n1,x2,:)
               !
            enddo 
            !
         enddo
      endif
      !
      x1 = 0
      !
      dzeta = 0
      ! 
      do x2 = 1,3
         if (ix/=lincoord) then 
            !
            x1 = x1 + 1
            do imode = 1,trove%Nmodes_e
               do jmode = 1,trove%Nmodes_e
                 !
                 z_t = 0
                 do jx = 1,3 
                    do kx = 1,3 
                       z_t = z_t +epsil(x2,jx,kx)* sum( trove%Amatrho(:,jx,imode,irho)*trove%Bmatrho(jmode,:,kx,irho) )
                    enddo
                 enddo
                 !
                 dzeta(x1,imode,jmode) = z_t
                 !
               enddo 
            enddo
         endif
         !
      enddo
      !
      if (manifold/=0) then 
        do imode = 1,trove%Nmodes_e
           do jmode = 1,trove%Nmodes_e
              !
              dzeta(Ng,imode,jmode) = sum( & ! trove%mass(:)*( 
                                           trove%dAmatrho(:,1,imode,irho,1)*trove%Bmatrho(jmode,:,1,irho)+&
                                           trove%dAmatrho(:,2,imode,irho,1)*trove%Bmatrho(jmode,:,2,irho)+&
                                           trove%dAmatrho(:,3,imode,irho,1)*trove%Bmatrho(jmode,:,3,irho)   )
              !
           enddo 
           !
        enddo
      endif
      !
      do jmode = 1,trove%Nmodes_e
         !
         do n1 = 1,Natoms
            do x1 = 1,3
              s_vib(jmode,n1,x1)%field(1,irho) = trove%Bmatrho(jmode,n1,x1,irho)
            enddo 
         enddo 
         !
      enddo
      !
      do imode = 1,trove%Nmodes_e
        !
        do n1 = 1,Natoms
           !
           do x1 = 1,3
             !
             do k0 = 1,iNcoeff1
                !
                powers(:) = FLIndexQ(:,k0)
                !
                do jmode = 1,trove%Nmodes_e
                   !
                   powers(jmode) = powers(jmode) + 1
                   !
                   k1 = FLQindex(trove%Nmodes_e,powers)
                   !
                   s_vib(imode,n1,x1)%field(k1,irho) = s_vib(imode,n1,x1)%field(k1,irho) & 
                                                      -sum( dzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,k0) )
                   !
                   powers(jmode) = powers(jmode) - 1
                   !
                 enddo
                 !
              enddo
              !
            enddo
            !
         enddo
      enddo
      !
      ! Total Vibrational angular momentum L^2
      ! This version of Corriolis dzeta_ does not skip the singular axis
      !
      if (FLl2_coeffs) then
         !
         dzeta_ = 0
         ! 
         do x2 = 1,3
            !
            do imode = 1,trove%Nmodes_e
               do jmode = 1,trove%Nmodes_e
                 !
                 z_t = 0
                 do jx = 1,3 
                    do kx = 1,3 
                       z_t = z_t +epsil(x2,jx,kx)* sum( trove%Amatrho(:,jx,imode,irho)*trove%Bmatrho(jmode,:,kx,irho) )
                    enddo
                 enddo
                 !
                 dzeta_(x2,imode,jmode) = z_t
                 !
               enddo 
            enddo
            !
         enddo
         !
         do imode = 1,trove%Nmodes_e
           !
           do jmode = 1,trove%Nmodes_e
              !
              powers = 0
              !
              do kmode = 1,trove%Nmodes_e
                !
                powers(kmode) = powers(kmode) + 1
                !
                do mmode = 1,trove%Nmodes_e
                   !
                   powers(mmode) = powers(mmode) + 1             
                   !
                   k1 = FLQindex(trove%Nmodes_e,powers)
                   !
                   trove%L2_vib(imode,jmode)%field(k1,irho) = trove%L2_vib(imode,jmode)%field(k1,irho) - & 
                                                              factor*sum( dzeta_(3:3,imode,kmode)*dzeta_(3:3,jmode,mmode) )
                   !
                   powers(mmode) = powers(mmode) - 1
                   !
                 enddo
                 !
                 powers(kmode) = powers(kmode) - 1
                 !
               enddo
               !
            enddo
            !
         enddo      
         !
      endif
      !
    enddo
    !$omp end do
    !
    !write(out,"(/'vib')") 
    !do irho = 0,Npoints
    !  if (job%verbose>=6) then 
    !    write(out,"(i8,18(<Nmodes>f16.8))") irho,((s_vib(Nmodes,n1,x1)%field(1,irho),n1=1,Natoms),x1=1,3)
    !  endif
    !enddo
    !
    !
    deallocate(s_mat_t,eta,tmat,Jmat,c,dzeta,dzeta_,powers)
    !$omp end parallel
    !
    if (trove%lincoord/=0.and..false.) then 
      !
      do imode = 1,3
        !
        do n1 = 1,Natoms
           !
           do x1 = 1,3
             !
             do k1 = 1,iNcoeff1
                !
                do irho = 0,Npoints
                  !
                  r_t(irho) = trove%rho_border(1)+trove%rhostep*irho
                  !
                enddo
                !
                do irho = 1,30
                   !
                   call MLratintark(r_t(31:100),s_rot(imode,n1,x1)%field(k1,31:100),r_t(irho),f_t,df_t)
                   !
                   s_rot(imode,n1,x1)%field(k1,irho) = f_t
                   !
                enddo
                !
                !continue
                !
             enddo
             !
           enddo
           !
        enddo
      enddo
      !
    endif
    !
    if (job%verbose>=2) write(out,"('...done!')") 
    !
    call TimerStop('s_vib_s_rot')
    !
    if (job%verbose>=4) write(out,"('s_vib_s_rot_Sorensen/end  ')") 
    !
  end subroutine s_vib_s_rot_Sorensen


  !
  ! This is a part of s_rot_polynom calcualtions, which  does the recursive Sorensen procedure
  !
  recursive subroutine s_rot_solve(dimen,maxorder,iNcoeff,jmat,eta)
    ! 
    integer(ik),intent(in)      :: dimen,maxorder,iNcoeff
    !
    real(ark),intent(in)  :: Jmat(dimen,dimen,0:trove%Nmodes)
    real(ark),intent(out) :: eta(dimen,dimen,iNcoeff)
    !
    integer(ik)            :: Nmodes,x1,x2,jNcoeff1,jNcoeff2,k0,k1,imode
    integer(ik)            :: n0,n,int
    integer(ik)            :: rank,iwork,info,alloc

    real(ark)              :: b(dimen),a(dimen,dimen)

    real(rk),allocatable    :: db(:,:),da(:,:)
    real(ark),allocatable    :: b_t(:,:,:)
    real(rk),allocatable    :: s_t(:),work(:)
    integer(ik),allocatable :: powers(:)
       !
       Nmodes = trove%Nmodes
       iwork = 50*dimen
       !
       allocate (db(dimen,1),da(dimen,dimen),s_t(dimen),work(iwork),powers(trove%Nmodes),b_t(dimen,dimen,iNcoeff),stat=alloc)
       if (alloc/=0) stop 's_rot_solve, powers,s_t,b_t - out of memory'
       !
       eta  = 0
       !
       eta(:,:,1) = jmat(:,:,0)
       !
       call MLinvmatark(jmat(1:dimen,1:dimen,0),eta(1:dimen,1:dimen,1),dimen,info)
       !
       if (info/=0) then
         !
         da = eta(:,:,1)
         !
         call lapack_ginverse(da)
         !
         eta(:,:,1) = da
         !
       endif
       !
       ! Higher order expansion:
       !
       jNcoeff1 = 0
       jNcoeff2 = 0
       !
       b_t = 0
       !
       do n  = 0,maxorder-1
          !
          jNcoeff1= jNcoeff2+1
          jNcoeff2= trove%RangeOrder(n)
          !
          do k0 = jNcoeff1,jNcoeff2
             !
             do x1 = 1,dimen
                !
                do x2 = 1,dimen
                   !
                   powers(:) = FLIndexQ(:,k0)
                   !
                   do imode = 1,trove%Nmodes_e
                      !
                      powers(imode) = powers(imode) + 1
                      !
                      k1 = FLQindex(trove%Nmodes_e,powers)
                      !
                      b_t(x1,x2,k1) = b_t(x1,x2,k1) + sum(eta(x1,:,k0)*Jmat(:,x2,imode))
                      !
                      powers(imode) = powers(imode) - 1
                      !
                   enddo
                   !
                enddo
                !
             enddo
             !
          enddo
          !
          do int = jNcoeff2+1,trove%RangeOrder(n+1)
             !
             do x1 = 1,dimen
                !
                a(:,:) = transpose(jmat(:,:,0))
                b(:) = -b_t(x1,:,int)
                !
                call MLlinurark(dimen,a,b,eta(x1,:,int),info)
                !
                if (info/=0) then
                  !
                  da = a
                  db(:,1) = b(:)
                  !
                  call lapack_gelss(da(:,:),db(:,:))
                  !
                  ! gelss - solves a linera equation by least squares method 
                  ! 
                  call dgelss(dimen,dimen,1,da(:,:),dimen,db(:,1),dimen,s_t,-1.0d-12, rank, work, iwork, info)
                  !
                  if (info/=0) then
                    write (out,"('s_rot_solve: dgelss returned ',i8)") info
                    stop 's_rot_solve - dgelss failed'
                  end if
                  !
                  eta(x1,:,int) = db(:,1)
                  !
                end if
                !
              enddo
              !
          enddo
          !
       enddo
       !
       deallocate (db,da,s_t,work,powers)
       !
  end subroutine s_rot_solve



!
! Here we calculate the jacobian matrices s_vib and s_rot for a given geometry. 
!
  subroutine s_vib_s_rot_dvr(dchi,irho,s_vib,s_rot)


    real(ark),intent(in)    :: dchi(trove%Nmodes)
    integer(ik),intent(in) :: irho
    real(ark),intent(out)   :: s_vib(trove%Nmodes,trove%Natoms,3,0:trove%Nmodes,0:trove%Nmodes)
    real(ark),intent(out)   :: s_rot(3,trove%Natoms,3,0:trove%Nmodes,0:trove%Nmodes)
    !
    integer(ik) :: n1,Ne,info
    !
    integer(ik) ::  x1,x2,k0,k1,Ng
    integer(ik) ::  n,nmodes,Natoms,iatom,ix,jx,kx
    integer(ik) ::  imode,jmode,ilin,q1,q2
    !
    real(ark)    :: Inertm(3)
    !
    real(ark)   :: z_t
    !
    real(ark)    :: s_mat_t(4,trove%Natoms,3,0:trove%Nmodes,0:trove%Nmodes)
    real(ark)    :: tmat0(trove%Nmodes),tmat1(trove%Nmodes,trove%Nmodes),tmat2(trove%Nmodes)
    real(ark)    :: r(trove%Natoms,3),dr(trove%Natoms,3),ddr(trove%Natoms,3),dddr(trove%Natoms,3)
    real(ark)    :: Jmat(4,4,0:trove%Nmodes,0:trove%Nmodes)
    real(ark)    :: eta(4,4,0:trove%Nmodes,0:trove%Nmodes)
    real(ark)    :: c(4,trove%Natoms,3),dc(4,trove%Natoms,3),ddc(4,trove%Natoms,3)
    real(ark)    :: dzeta(4,trove%Nmodes,trove%Nmodes),ddzeta(4,trove%Nmodes,trove%Nmodes),dddzeta(4,trove%Nmodes,trove%Nmodes)
    real(rk)     :: deta(4,4)
    !
    !call TimerStart('s_vib_s_rot_dvr')
    !
    ! Sumstitution for easier reference within the procedure 
    !
    Nmodes = trove%Nmodes
    Natoms = trove%Natoms
    Ne     = trove%Nmodes_e
    !
    s_vib = 0
    s_rot = 0
    !
    do n1 = 1,trove%Natoms
       do ix = 1,3
          !
          r(n1,ix) = trove%b0(n1,ix,irho)    + sum(  trove%Amatrho (n1,ix,1:Ne,irho  )*dchi(1:Ne))
          !
          if (manifold/=0) then 
              !
              dr(n1,ix) = trove%db0(n1,ix,irho,1) + sum(  trove%dAmatrho(n1,ix,1:Ne,irho,1)*dchi(1:Ne))
             ddr(n1,ix) = trove%db0(n1,ix,irho,2) + sum(  trove%dAmatrho(n1,ix,1:Ne,irho,2)*dchi(1:Ne))
            dddr(n1,ix) = trove%db0(n1,ix,irho,3) + sum(  trove%dAmatrho(n1,ix,1:Ne,irho,3)*dchi(1:Ne))
          endif 
          !
       enddo
    enddo
    !
    !call FLfromcartesian2local(xyz,r)
    !
    !
    !
    ! Inertia moments, zero order (equilibrium)
    !
    !Inertm(1) = sum( trove%mass(:)*( trove%b0(:,2,irho)**2+ trove%b0(:,3,irho)**2 ) )
    !Inertm(2) = sum( trove%mass(:)*( trove%b0(:,1,irho)**2+ trove%b0(:,3,irho)**2 ) )
    !Inertm(3) = sum( trove%mass(:)*( trove%b0(:,1,irho)**2+ trove%b0(:,2,irho)**2 ) )
    !
    ilin = trove%lincoord
    !
    ! Check if all Inertia moments are non zero 
    !
    !do x1 = 1,3
    !  !
    !  do x2 = x1+1,x1-1
    !    !
    !    if (sum( trove%mass(:)*trove%b0(:,x1,irho)*trove%b0(:,x2,irho) )>sqrt(small_)) then 
    !       write(out,"('s_mat-dvr: illegal use of Sorensen, non diagonal Inertia moment non zero:',2i4,f20.10)") x1,x2,trove%mass(:)*trove%b0(:,x1,irho)*trove%b0(:,x2,irho)
    !       stop 's-mat-dvr: non zero a non diagonal Inertia moment'
    !    endif
    !    !
    !  enddo
    !enddo
    !
    ! Maximal number of the linear terms 
    !

    Jmat = 0
    !
    do x1 = 1,3
       !
       do x2 = 1,3
          !
          ! Jmat Sorensen matrix 
          !
          Jmat(x1,x2,0,0) = jmat_gg(x1,x2,trove%b0(:,:,irho),r(:,:))
          !
          !if (Jmat(x1,x2,0,0)<0) then 
          !   write(out,"('s-mat-dvr: the Jmat is negative for ',2i4,': ',f20.10)") x1,x2,Jmat(x1,x2,0,0)
          !   stop 's-mat-dvr: the Jmat is negative for'
          !endif 
          !
          !
          ! deriv. of Jmat wrt chi_i
          !
          do imode = 1,Ne
            !
            Jmat(x1,x2,imode,0) = jmat_gg(x1,x2,trove%b0(:,:,irho),trove%Amatrho(:,:,imode,irho))
            !
          enddo
          !
          if (manifold/=0) then 
            !
            ! deriv. of Jmat wrt rho
            !
            Jmat(x1,x2,Nmodes,0) = jmat_gg( x1,x2,trove%db0(:,:,irho,1), r(:,:) ) + &
                                   jmat_gg( x1,x2,trove%b0 (:,:,irho  ),dr(:,:) )
            !
            ! second deriv. of Jmat wrt rho^2
            !
            Jmat(x1,x2,Nmodes,Nmodes) =   jmat_gg( x1,x2,trove%db0(:,:,irho,2),  r(:,:) ) + &
                                  2.0_ark*jmat_gg( x1,x2,trove%db0(:,:,irho,1), dr(:,:) ) + &
                                          jmat_gg( x1,x2,trove%b0 (:,:,irho  ),ddr(:,:) )
            !
            ! second deriv. of Jmat wrt rho and chi_i
            !
            do imode = 1,Ne
              !
              Jmat(x1,x2,Nmodes,imode) = jmat_gg( x1,x2,trove%db0(:,:,irho,1),trove%Amatrho (:,:,imode,irho  ) ) + &
                                         jmat_gg( x1,x2,trove%b0 (:,:,irho  ),trove%dAmatrho(:,:,imode,irho,1) )
              !
              Jmat(x1,x2,imode,Nmodes) = Jmat(x1,x2,Nmodes,imode)
              !
            enddo
            !
          endif
          !
       enddo
       !
    enddo


    if (manifold/=0) then 
      !
      do x1 = 1,3
          !
          ! Jmat Sorensen matrix 
          !
          Jmat(4,x1,0,0) = jmat_grho(x1,trove%db0(:,:,irho,1),r(:,:))
          !
          !
          ! deriv. of Jmat wrt chi_i
          !
          do imode = 1,Ne
            !
            Jmat(4,x1,imode,0) = jmat_grho(x1,trove%db0(:,:,irho,1),trove%Amatrho(:,:,imode,irho))
            !
          enddo
          !
          ! deriv. of Jmat wrt rho
          !
          Jmat(4,x1,Nmodes,0) = jmat_grho( x1,trove%db0(:,:,irho,2), r(:,:) ) + &
                                jmat_grho( x1,trove%db0(:,:,irho,1),dr(:,:) )
          !
          ! second deriv. of Jmat wrt rho^2
          !
          Jmat(4,x1,Nmodes,Nmodes) =   jmat_grho( x1,trove%db0(:,:,irho,3),  r(:,:) ) + &
                               2.0_ark*jmat_grho( x1,trove%db0(:,:,irho,2), dr(:,:) ) + &
                                       jmat_grho( x1,trove%db0(:,:,irho,1),ddr(:,:) )
          !
          ! second deriv. of Jmat wrt rho and chi_i
          !
          do imode = 1,Ne
            !
            Jmat(4,x1,Nmodes,imode) = jmat_grho( x1,trove%db0(:,:,irho,2),trove%Amatrho (:,:,imode,irho  ) ) + &
                                      jmat_grho( x1,trove%db0(:,:,irho,1),trove%dAmatrho(:,:,imode,irho,1) )
            !
            Jmat(4,x1,imode,Nmodes) = Jmat(4,x1,Nmodes,imode)
            !
          enddo
          !
          ! interchange x1 and 4:
          !
          ! Jmat Sorensen matrix 
          !
          Jmat(x1,4,0,0) = jmat_grho(x1,dr(:,:),trove%b0(:,:,irho))
          !
          !
          ! deriv. of Jmat wrt chi_i
          !
          do imode = 1,Ne
            !
            Jmat(x1,4,imode,0) = jmat_grho(x1,trove%dAmatrho(:,:,imode,irho,1),trove%b0(:,:,irho))
            !
          enddo
          !
          ! deriv. of Jmat wrt rho
          !
          Jmat(x1,4,Nmodes,0) = jmat_grho( x1, dr(:,:),trove%db0(:,:,irho,1)) + &
                                jmat_grho( x1,ddr(:,:),trove%b0 (:,:,irho)  )
          !
          ! second deriv. of Jmat wrt rho^2
          !
          Jmat(x1,4,Nmodes,Nmodes) =   jmat_grho( x1,  dr(:,:),trove%db0(:,:,irho,2) ) + &
                               2.0_ark*jmat_grho( x1, ddr(:,:),trove%db0(:,:,irho,1) ) + &
                                       jmat_grho( x1,dddr(:,:),trove%b0 (:,:,irho  ) )
          !
          ! second deriv. of Jmat wrt rho and chi_i
          !
          do imode = 1,Ne
            !
            Jmat(x1,4,Nmodes,imode) = jmat_grho( x1,trove%dAmatrho(:,:,imode,irho,1),trove%db0(:,:,irho,1) ) + &
                                      jmat_grho( x1,trove%dAmatrho(:,:,imode,irho,2),trove%b0 (:,:,irho  ) )
            !
            Jmat(x1,4,imode,Nmodes) = Jmat(x1,4,Nmodes,imode)
            !
          enddo
          !
          !Jmat(x1,4,:,:) = Jmat(4,x1,:,:)
          !
      enddo
      !
      ! the last rho,rho term 
      !
      ! Jmat Sorensen matrix 
      !
      Jmat(4,4,0,0) = jmat_rhorho(trove%db0(:,:,irho,1),dr(:,:))
      !
      !
      ! deriv. of Jmat wrt chi_i
      !
      do imode = 1,Ne
        !
        Jmat(4,4,imode,0) = jmat_rhorho(trove%db0(:,:,irho,1),trove%dAmatrho(:,:,imode,irho,1))
        !
      enddo
      !
      ! deriv. of Jmat wrt rho
      !
      Jmat(4,4,Nmodes,0) = jmat_rhorho( trove%db0(:,:,irho,2), dr(:,:) ) + &
                           jmat_rhorho( trove%db0(:,:,irho,1),ddr(:,:) )
      !
      ! second deriv. of Jmat wrt rho^2
      !
      Jmat(4,4,Nmodes,Nmodes) =   jmat_rhorho( trove%db0(:,:,irho,3),  dr(:,:) ) + &
                          2.0_ark*jmat_rhorho( trove%db0(:,:,irho,2), ddr(:,:) ) + &
                                  jmat_rhorho( trove%db0(:,:,irho,1),dddr(:,:) )
      !
      ! second deriv. of Jmat wrt rho and chi_i
      !
      do imode = 1,Ne
        !
        Jmat(4,4,Nmodes,imode) = jmat_rhorho( trove%db0(:,:,irho,2),trove%dAmatrho(:,:,imode,irho,1) ) + &
                                 jmat_rhorho( trove%db0(:,:,irho,1),trove%dAmatrho(:,:,imode,irho,2) )
        !
        Jmat(4,4,imode,Nmodes) = Jmat(4,4,Nmodes,imode)
        !
      enddo
      !
    endif 

    if (ilin/=0) then 
      !
      do imode = 0,Nmodes
         do jmode = 0,Nmodes
           !
           Jmat(ilin:,ilin:,imode,jmode) = eoshift(Jmat(ilin:,ilin:,imode,jmode),1,dim=1)
           Jmat(ilin:,ilin:,imode,jmode) = eoshift(Jmat(ilin:,ilin:,imode,jmode),1,dim=2)
           !
         enddo
         !
      enddo
      !
    endif 
    !
    Ng = 0
    !
    do ix = 1,3
       !
       if (ix/=ilin) then 
         !
         Ng = Ng + 1
         !
         do n1 = 1,trove%Natoms
            do x1 = 1,3
               !
               c(Ng,n1,x1) = trove%mass(n1)*sum(epsil(x1,ix,:)*trove% b0(n1,:,irho  ))
               !
               if (manifold/=0) then 
                  !
                  dc(Ng,n1,x1) = trove%mass(n1)*sum(epsil(x1,ix,:)*trove%db0(n1,:,irho,1))
                 ddc(Ng,n1,x1) = trove%mass(n1)*sum(epsil(x1,ix,:)*trove%db0(n1,:,irho,2))
               endif 
               !
            enddo
         enddo
         !
       endif 
       !
    enddo
    !
    if (manifold/=0) then
       !
       Ng = Ng + 1
       !
       do n1 = 1,trove%Natoms
          do x1 = 1,3
             !
             c(Ng,n1,x1) = trove%mass(n1)*trove%db0(n1,x1,irho,1)
             !
             if (manifold/=0) then 
                !
                dc(Ng,n1,x1) = trove%mass(n1)*trove%db0(n1,x1,irho,2)
               ddc(Ng,n1,x1) = trove%mass(n1)*trove%db0(n1,x1,irho,3)
             endif 
             !
          enddo
       enddo
    endif
    !
    ! We must control the relations involving the translational t-vectors
    !
    do x1 = 1,Ng
       !
       z_t = sum(c(x1,:,1)+c(x1,:,2)+c(x1,:,3))
       !
       if (abs(z_t)>sqrt(small_)) then 
           write(out,"('s-mat-dvr: the translational relation for c is not zero for g = ',i4,': ',f20.10)") x1,z_t
           stop 's-mat-dvr: the translational relation for c is not zero'
       endif
       !
    enddo
    !
    ! Check if the Jmat is symmetric
    !
    do ix = 1,Ng
       do jx = 1,Ng
         !
         if ( any( abs( Jmat(ix,jx,:,:)-Jmat(jx,ix,:,:) ) >10000.0_ark*sqrt(small_)  ) ) then 
           write(out,"('s-mat-dvr: jmat is not symmetric for ',i4,',',i4)") ix,jx
           write(out,"(100f20.10)") Jmat(ix,jx,:,:)
           write(out,"(100f20.10)") Jmat(jx,ix,:,:)
           !
           stop 's-mat-dvr: jmat is not symmetric'
           !
         endif 
         !
       enddo
       !
    enddo
    !
    !eta = 0
    !
    !eta(1:Ng,1:Ng,0,0) = jmat(1:Ng,1:Ng,0,0)
    !
    call MLinvmatark(jmat(1:Ng,1:Ng,0,0),eta(1:Ng,1:Ng,0,0),Ng,info)
    !
    if (info/=0) then 
     !
     deta(1:Ng,1:Ng) = eta(1:Ng,1:Ng,0,0)
     call lapack_ginverse(deta(1:Ng,1:Ng))
     !
     eta(1:Ng,1:Ng,0,0) = deta(1:Ng,1:Ng)
     !
    endif
    !
    do imode = 1,Nmodes
       !
       eta(1:Ng,1:Ng,imode,0) = -matmul( matmul( eta(1:Ng,1:Ng,0,0),jmat(1:Ng,1:Ng,imode,0) ),eta(1:Ng,1:Ng,0,0) )
       !
       do jmode = 1,imode
         !
         eta(1:Ng,1:Ng,imode,jmode) = - ( matmul ( matmul( eta(1:Ng,1:Ng,jmode,0),jmat(1:Ng,1:Ng,imode,0    ) ),&
                                         eta(1:Ng,1:Ng,0    ,0) )  & 
                                         +matmul ( matmul( eta(1:Ng,1:Ng,0    ,0),jmat(1:Ng,1:Ng,imode,0    ) ),&
                                         eta(1:Ng,1:Ng,jmode,0) )  &
                                         +matmul ( matmul( eta(1:Ng,1:Ng,0,    0),jmat(1:Ng,1:Ng,imode,jmode) ),&
                                         eta(1:Ng,1:Ng,0    ,0) ) )
         ! 
         eta(1:Ng,1:Ng,jmode,imode) = eta(1:Ng,1:Ng,imode,jmode)
         !
       enddo
       !
    enddo
    !
    s_mat_t = 0
    !
    do n1 = 1,Natoms
       do x2 = 1,3
          do x1 = 1,Ng
             !
             s_mat_t(x1,n1,x2,0,0) = sum( eta(x1,1:Ng,0,0)*c(1:Ng,n1,x2) )
             !
             ! deriv. of smat wrt chi_i
             !
             do imode = 1,Ne
               !
               s_mat_t(x1,n1,x2,imode,0) = sum( eta(x1,1:Ng,imode,0)*c(1:Ng,n1,x2) )
               !
               ! second deriv. of smat wrt chi_i and chi_j
               !
               do jmode = 1,imode
                 !
                 s_mat_t(x1,n1,x2,imode,jmode) = sum( eta(x1,1:Ng,imode,jmode)*c(1:Ng,n1,x2) )
                 s_mat_t(x1,n1,x2,jmode,imode) = s_mat_t(x1,n1,x2,imode,jmode)
                 !
               enddo
               !
             enddo
             !
             if (manifold/=0) then 
               !
               ! deriv. of smat wrt rho
               !
               s_mat_t(x1,n1,x2,Nmodes,0) = sum( eta(x1,1:Ng,Nmodes,0)*c(1:Ng,n1,x2) ) + & 
                                            sum( eta(x1,1:Ng,0,0)*dc(1:Ng,n1,x2) )
               !
               ! second deriv. of smat wrt rho and chi_i
               !
               do imode = 1,Ne
                 !
                 s_mat_t(x1,n1,x2,Nmodes,imode) = sum( eta(x1,1:Ng,imode,Nmodes)*c(1:Ng,n1,x2) ) + &
                                                  sum( eta(x1,1:Ng,imode,0)*dc(1:Ng,n1,x2) )
                 !
                 s_mat_t(x1,n1,x2,imode,Nmodes) = s_mat_t(x1,n1,x2,Nmodes,imode)
                 !
               enddo
               !
               ! second deriv. of smat wrt rho^2
               !
               s_mat_t(x1,n1,x2,Nmodes,Nmodes) =        sum( eta(x1,1:Ng,Nmodes,Nmodes)*c(1:Ng,n1,x2) ) + & 
                                                2.0_ark*sum( eta(x1,1:Ng,Nmodes,0)*dc(1:Ng,n1,x2) )    + &
                                                        sum( eta(x1,1:Ng,0,0)*ddc(1:Ng,n1,x2) ) 
               !
             endif
          enddo
       enddo 
    enddo
    !
    do n1 = 1,Natoms
       !
       x1 = 0
       do ix = 1,3
         if (ix/=ilin) then 
            !
            x1 = x1 + 1
            do x2 = 1,3
              !
              s_rot(ix,n1,x2,:,:)= s_mat_t(x1,n1,x2,:,:)
              !
            enddo
            !
         endif
       enddo 
       !
    enddo
    !
    if (manifold/=0) then 
       do n1 = 1,Natoms
          do x2 = 1,3
             !
             s_vib(trove%Nmodes_n,n1,x2,:,:) =  s_mat_t(Ng,n1,x2,:,:)
             !
          enddo 
          !
       enddo
    endif
    !
    x1 = 0
    !
    dzeta = 0
    ! 
    do ix = 1,3
       if (ix/=ilin) then 
          !
          x1 = x1 + 1
          do imode = 1,Ne
             do jmode = 1,Ne
               !
               z_t = 0
               do jx = 1,3 
                  do kx = 1,3 
                     z_t = z_t +epsil(ix,jx,kx)* sum( trove%Amatrho(:,jx,imode,irho)*trove%Bmatrho(jmode,:,kx,irho) )
                  enddo
               enddo
               !
               dzeta(x1,imode,jmode) = z_t
               !
               if (manifold/=0) then
                 !
                 z_t = 0
                 do jx = 1,3 
                    do kx = 1,3 
                       z_t = z_t +epsil(ix,jx,kx)*( sum( trove%dAmatrho(:,jx,imode,irho,1)*trove% Bmatrho(jmode,:,kx,irho ) ) + &
                                                    sum( trove% Amatrho(:,jx,imode,irho  )*trove%dBmatrho(jmode,:,kx,irho,1) ) )
                    enddo
                 enddo
                 !
                 ddzeta(x1,imode,jmode) = z_t
                 !
                 z_t = 0
                 do jx = 1,3 
                    do kx = 1,3 
                       z_t = z_t +epsil(ix,jx,kx)*(   sum( trove%dAmatrho(:,jx,imode,irho,2)*trove% Bmatrho(jmode,:,kx,irho  ) ) + &
                                              2.0_ark*sum( trove%dAmatrho(:,jx,imode,irho,1)*trove%dBmatrho(jmode,:,kx,irho,1) ) + &
                                                      sum( trove% Amatrho(:,jx,imode,irho  )*trove%dBmatrho(jmode,:,kx,irho,2) ) )
                    enddo
                 enddo
                 !
                 dddzeta(x1,imode,jmode) = z_t
                 !
               endif
               !
             enddo 
          enddo
       endif
       !
       if (manifold/=0) then 
         do imode = 1,Ne
            do jmode = 1,Ne
                 !
                 dzeta(Ng,imode,jmode) = sum( &
                                              trove%dAmatrho(:,1,imode,irho,1)*trove%Bmatrho(jmode,:,1,irho)+&
                                              trove%dAmatrho(:,2,imode,irho,1)*trove%Bmatrho(jmode,:,2,irho)+&
                                              trove%dAmatrho(:,3,imode,irho,1)*trove%Bmatrho(jmode,:,3,irho)   )
                                              !
               ddzeta(Ng,imode,jmode) = sum( &
                                             trove%dAmatrho(:,1,imode,irho,2)*trove% Bmatrho(jmode,:,1,irho  )+&
                                             trove%dAmatrho(:,2,imode,irho,2)*trove% Bmatrho(jmode,:,2,irho  )+&
                                             trove%dAmatrho(:,3,imode,irho,2)*trove% Bmatrho(jmode,:,3,irho  )+&
                                             trove%dAmatrho(:,1,imode,irho,1)*trove%dBmatrho(jmode,:,1,irho,1)+&
                                             trove%dAmatrho(:,2,imode,irho,1)*trove%dBmatrho(jmode,:,2,irho,1)+&
                                             trove%dAmatrho(:,3,imode,irho,1)*trove%dBmatrho(jmode,:,3,irho,1)   )
                                             !
              dddzeta(Ng,imode,jmode) = sum( &
                                               trove%dAmatrho(:,1,imode,irho,3)*trove% Bmatrho(jmode,:,1,irho  )+&
                                               trove%dAmatrho(:,2,imode,irho,3)*trove% Bmatrho(jmode,:,2,irho  )+&
                                               trove%dAmatrho(:,3,imode,irho,3)*trove% Bmatrho(jmode,:,3,irho  )+&
                                       2.0_ark*trove%dAmatrho(:,1,imode,irho,2)*trove%dBmatrho(jmode,:,1,irho,1)+&
                                       2.0_ark*trove%dAmatrho(:,2,imode,irho,2)*trove%dBmatrho(jmode,:,2,irho,1)+&
                                       2.0_ark*trove%dAmatrho(:,3,imode,irho,2)*trove%dBmatrho(jmode,:,3,irho,1)+&  
                                               trove%dAmatrho(:,1,imode,irho,1)*trove%dBmatrho(jmode,:,1,irho,2)+&
                                               trove%dAmatrho(:,2,imode,irho,1)*trove%dBmatrho(jmode,:,2,irho,2)+&
                                               trove%dAmatrho(:,3,imode,irho,1)*trove%dBmatrho(jmode,:,3,irho,2) )
            enddo 
            !
         enddo
       endif
       !
    enddo
    !
    do imode = 1,Ne
      !
      do n1 = 1,Natoms
         !
         do x1 = 1,3
           !
           do jmode = 1,Ne
              !
              tmat0(jmode) = -sum( dzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,0,0) )
              !
           enddo
           !
           s_vib(imode,n1,x1,0,0) = trove%Bmatrho(imode,n1,x1,irho)+sum( tmat0(1:Ne)*dchi(1:Ne) )
           !
           do q1 = 1,Ne
              !
              do jmode = 1,Ne
                 !
                 tmat1(jmode,q1) = -sum( dzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,q1,0) )
                 !
              enddo
              !
              s_vib(imode,n1,x1,q1,0) = sum( tmat1(1:Ne,q1)*dchi(1:Ne) )+tmat0(q1)
              !
              do q2 = 1,q1
                !
                do jmode = 1,Ne
                   !
                   tmat2(jmode) = -sum( dzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,q1,q2) )
                   !
                enddo
                !
                s_vib(imode,n1,x1,q1,q2) = sum( tmat2(1:Ne)*dchi(1:Ne) )+tmat1(q1,q2)+tmat1(q2,q1)
                s_vib(imode,n1,x1,q2,q1) = s_vib(imode,n1,x1,q1,q2)
                !
              enddo
           enddo
           !
           if (manifold/=0) then
             !
             ! d s / d rho 
             !
             do jmode = 1,Ne
                !
                tmat0(jmode) = -sum( ddzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,0,0)+dzeta(1:Ng,jmode,imode)*&
                                                              s_mat_t(1:Ng,n1,x1,Nmodes,0) )
                !
             enddo
             !
             s_vib(imode,n1,x1,Nmodes,0) = trove%dBmatrho(imode,n1,x1,irho,1)+sum( tmat0(1:Ne)*dchi(1:Ne) )
             !
             ! d^2 s / d rho^2 
             !
             do jmode = 1,Ne
                !
                tmat1(jmode,Nmodes) = -sum(  dddzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,0,0) + & 
                                      2.0_ark*ddzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,Nmodes,0) + &
                                               dzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,Nmodes,Nmodes) )
                !
             enddo
             !
             s_vib(imode,n1,x1,Nmodes,Nmodes) = trove%dBmatrho(imode,n1,x1,irho,2)+sum( tmat1(1:Ne,Nmodes)*dchi(1:Ne) )
             !
             ! d^2 s / d rho d xi_i 
             !
             do q1 = 1,Ne
                !
                do jmode = 1,Ne
                   !
                   tmat1(jmode,q1) = -sum(ddzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,q1,0) + & 
                                           dzeta(1:Ng,jmode,imode)*s_mat_t(1:Ng,n1,x1,q1,Nmodes) )
                   !
                enddo
                !
                s_vib(imode,n1,x1,q1,Nmodes) = sum( tmat1(1:Ne,q1)*dchi(1:Ne) )+tmat0(q1)
                s_vib(imode,n1,x1,Nmodes,q1) = s_vib(imode,n1,x1,q1,Nmodes)
                !
             enddo
             !
           endif
           !
          enddo
          !
       enddo
    enddo
    !
    !call TimerStop('s_vib_s_rot_dvr')
    !
    continue
    !
    contains 
    !
    function vec_scal_vec(a,b,c,d) result (v)
      !
      real(ark),intent(in) :: a(3),b(3),c(3),d(3)
      real(ark) :: v
      !
      v = sum(a*c)*sum(b*d)-sum(a*d)*sum(b*c)
      !
    end function vec_scal_vec
    !
    function jmat_rhorho(a,b) result (v)
      !
      real(ark),intent(in)   :: a(:,:),b(:,:)
      real(ark)              :: v
       !
       v = sum( trove%mass(:)*( a(:,1)*b(:,1)+a(:,2)*b(:,2)+a(:,3)*b(:,3) ) )
       !
    end function jmat_rhorho
    !
    function jmat_gg(g1,g2,a,b) result (v)
      !
      integer(ik),intent(in) :: g1,g2
      real(ark),intent(in)   :: a(:,:),b(:,:)
      real(ark)              :: v
       !
       v = -sum( trove%mass(:)*b(:,g1)*a(:,g2) )
       !
       if (g1==g2) v = v + sum( trove%mass(:)*( a(:,1)*b(:,1)+a(:,2)*b(:,2)+a(:,3)*b(:,3) ) )
       !
    end function jmat_gg
    !
    function jmat_grho(g,a,b) result (v)
      !
      integer(ik),intent(in) :: g
      real(ark),intent(in)   :: a(:,:),b(:,:)
      real(ark)              :: v
      integer(ik)            :: n,ix
       !
       v = 0 
       !
       do n = 1,trove%Natoms
         !
         do ix  = 1,3
           !
           v = v + sum(trove%mass(n)*a(n,:)*epsil(:,g,ix)*b(n,ix) )
           !
         enddo
         !
       enddo
       !
    end function jmat_grho
    !
    function vector_product(v1,v2) result (v)
      !
      real(ark),intent(in) :: v1(3),v2(3)
      real(ark) :: v(3)
      !
      v(1) = v1(2)*v2(3)-v1(3)*v2(2)
      !
      v(2) = v1(3)*v2(1)-v1(1)*v2(3)
      !
      v(3) = v1(1)*v2(2)-v1(2)*v2(1)
      !
    end function vector_product
    !
  end subroutine s_vib_s_rot_dvr




!
! Here we calculate the jacobian matrices s_vib and s_rot for a given geometry. 
!
  subroutine s_vib_s_rot_local(dchi,irho,s_vib,s_rot)


    real(ark),intent(in)    :: dchi(trove%Nmodes)
    integer(ik),intent(in) :: irho
    real(ark),intent(out)   :: s_vib(trove%Nmodes,trove%Natoms,3)
    real(ark),intent(out)   :: s_rot(3,trove%Natoms,3)
    !
    integer(ik) :: n1,Ne,info
    !
    integer(ik) ::  x1,x2,k0,k1,Ng
    integer(ik) ::  n,nmodes,Natoms,iatom,ix,jx,kx
    integer(ik) ::  imode,jmode,ilin,q1,q2
    !
    real(ark)    :: Inertm(3)
    !
    real(ark)   :: z_t
    !
    real(ark)    :: s_mat_t(4,trove%Natoms,3)
    real(ark)    :: r(trove%Natoms,3),dr(trove%Natoms,3)
    real(ark)    :: Jmat(4,4)
    real(ark)    :: eta(4,4)
    real(ark)    :: c(4,trove%Natoms,3)
    real(rk)     :: deta(4,4)
    !
    !call TimerStart('s_vib_s_rot_local')
    !
    ! Sumstitution for easier reference within the procedure 
    !
    Nmodes = trove%Nmodes
    Natoms = trove%Natoms
    Ne     = trove%Nmodes_e
    !
    s_vib = 0
    s_rot = 0
    !
    do n1 = 1,trove%Natoms
       do ix = 1,3
          !
          r(n1,ix) = trove%b0(n1,ix,irho)    + sum(  trove%Amatrho (n1,ix,1:Ne,irho  )*dchi(1:Ne))
          !
          !if (manifold/=0) then 
          !    !
          !    dr(n1,ix) = trove%db0(n1,ix,irho,1) + sum(  trove%dAmatrho(n1,ix,1:Ne,irho,1)*dchi(1:Ne))
          !   ddr(n1,ix) = trove%db0(n1,ix,irho,2) + sum(  trove%dAmatrho(n1,ix,1:Ne,irho,2)*dchi(1:Ne))
          !  dddr(n1,ix) = trove%db0(n1,ix,irho,3) + sum(  trove%dAmatrho(n1,ix,1:Ne,irho,3)*dchi(1:Ne))
          !endif 
          !
       enddo
    enddo
    !
    ilin = trove%lincoord

    Jmat = 0
    !
    do x1 = 1,3
       !
       do x2 = 1,3
          !
          ! Jmat Sorensen matrix 
          !
          Jmat(x1,x2) = jmat_gg(x1,x2,trove%b0(:,:,irho),r(:,:))
          !
          !
       enddo
       !
    enddo


    if (manifold/=0) then 
      !
      do x1 = 1,3
          !
          ! Jmat Sorensen matrix 
          !
          Jmat(4,x1) = jmat_grho(x1,trove%db0(:,:,irho,1),r(:,:))
          !
          !
          ! interchange x1 and 4:
          !
          ! Jmat Sorensen matrix 
          !
          Jmat(x1,4) = jmat_grho(x1,dr(:,:),trove%b0(:,:,irho))
          !
          !
      enddo
      !
      ! the last rho,rho term 
      !
      ! Jmat Sorensen matrix 
      !
      Jmat(4,4) = jmat_rhorho(trove%db0(:,:,irho,1),dr(:,:))
      !
    endif 

    if (ilin/=0) then 
      !
      Jmat(ilin:,ilin:) = eoshift(Jmat(ilin:,ilin:),1,dim=1)
      Jmat(ilin:,ilin:) = eoshift(Jmat(ilin:,ilin:),1,dim=2)
      !
    endif 
    !
    Ng = 0
    !
    do ix = 1,3
       !
       if (ix/=ilin) then 
         !
         Ng = Ng + 1
         !
         do n1 = 1,trove%Natoms
            do x1 = 1,3
               !
               c(Ng,n1,x1) = trove%mass(n1)*sum(epsil(x1,ix,:)*trove% b0(n1,:,irho  ))
               !
            enddo
         enddo
         !
       endif 
       !
    enddo
    !
    if (manifold/=0) then
       !
       Ng = Ng + 1
       !
       do n1 = 1,trove%Natoms
          do x1 = 1,3
             !
             c(Ng,n1,x1) = trove%mass(n1)*trove%db0(n1,x1,irho,1)
             !
          enddo
       enddo
    endif
    !
    ! We must control the relations involving the translational t-vectors
    !
    do x1 = 1,Ng
       !
       z_t = sum(c(x1,:,1)+c(x1,:,2)+c(x1,:,3))
       !
       if (abs(z_t)>sqrt(small_)) then 
           write(out,"('s-mat-dvr: the translational relation for c is not zero for g = ',i4,': ',f20.10)") x1,z_t
           stop 's-mat-dvr: the translational relation for c is not zero'
       endif
       !
    enddo
    !
    ! Check if the Jmat is symmetric
    !
    do ix = 1,Ng
       do jx = 1,Ng
         !
         if ( abs( Jmat(ix,jx)-Jmat(jx,ix) ) >10000.0_ark*sqrt(small_) ) then 
           write(out,"('s-mat-dvr: jmat is not symmetric for ',i4,',',i4)") ix,jx
           write(out,"(f20.10)") Jmat(ix,jx)
           write(out,"(f20.10)") Jmat(jx,ix)
           !
           stop 's-mat-dvr: jmat is not symmetric'
           !
         endif 
         !
       enddo
       !
    enddo
    !
    call MLinvmatark(jmat(1:Ng,1:Ng),eta(1:Ng,1:Ng),Ng,info)
    !
    if (info/=0) then 
     !
     deta(1:Ng,1:Ng) = eta(1:Ng,1:Ng)
     call lapack_ginverse(deta(1:Ng,1:Ng))
     !
     eta(1:Ng,1:Ng) = deta(1:Ng,1:Ng)
     !
    endif
    !
    s_mat_t = 0
    !
    do n1 = 1,Natoms
       do x2 = 1,3
          do x1 = 1,Ng
             !
             s_mat_t(x1,n1,x2) = sum( eta(x1,1:Ng)*c(1:Ng,n1,x2) )
             !
          enddo
       enddo 
    enddo
    !
    do n1 = 1,Natoms
       !
       x1 = 0
       do ix = 1,3
         if (ix/=ilin) then 
            !
            x1 = x1 + 1
            do x2 = 1,3
              !
              s_rot(ix,n1,x2)= s_mat_t(x1,n1,x2)
              !
            enddo
            !
         endif
       enddo 
       !
    enddo
    !
    if (manifold/=0) then 
       do n1 = 1,Natoms
          do x2 = 1,3
             !
             s_vib(trove%Nmodes_n,n1,x2) =  s_mat_t(Ng,n1,x2)
             !
          enddo 
          !
       enddo
    endif
    !
    do imode = 1,Ne
      !
      do n1 = 1,Natoms
         !
         do x1 = 1,3
           !
           !s_vib(imode,n1,x1) = trove%Bmatrho(imode,n1,x1,irho)+sum( tmat0(1:Ne)*dchi(1:Ne) )
           !
         enddo
         !
       enddo
    enddo
    !
    contains 
    !
    function jmat_rhorho(a,b) result (v)
      !
      real(ark),intent(in)   :: a(:,:),b(:,:)
      real(ark)              :: v
       !
       v = sum( trove%mass(:)*( a(:,1)*b(:,1)+a(:,2)*b(:,2)+a(:,3)*b(:,3) ) )
       !
    end function jmat_rhorho
    !
    function jmat_gg(g1,g2,a,b) result (v)
      !
      integer(ik),intent(in) :: g1,g2
      real(ark),intent(in)   :: a(:,:),b(:,:)
      real(ark)              :: v
       !
       v = -sum( trove%mass(:)*b(:,g1)*a(:,g2) )
       !
       if (g1==g2) v = v + sum( trove%mass(:)*( a(:,1)*b(:,1)+a(:,2)*b(:,2)+a(:,3)*b(:,3) ) )
       !
    end function jmat_gg
    !
    function jmat_grho(g,a,b) result (v)
      !
      integer(ik),intent(in) :: g
      real(ark),intent(in)   :: a(:,:),b(:,:)
      real(ark)              :: v
      integer(ik)            :: n,ix
       !
       v = 0 
       !
       do n = 1,trove%Natoms
         !
         do ix  = 1,3
           !
           v = v + sum(trove%mass(n)*a(n,:)*epsil(:,g,ix)*b(n,ix) )
           !
         enddo
         !
       enddo
       !
    end function jmat_grho
    !
  end subroutine s_vib_s_rot_local



  !
  ! Here we initialize the kinetic operator fields: g_vib, g_rot, and g_cor
  subroutine FLDVR_gmat_dvr(dchi,irho,g_vib,g_rot,g_cor,pseudo)

    real(ark),intent(in)    :: dchi(trove%Nmodes)
    integer(ik),intent(in) :: irho

    real(ark),intent(out) ::  g_vib(trove%Nmodes,trove%Nmodes)    ! Vibrational part of Kinetic factor G
    real(ark),intent(out) ::  g_rot(3,3)              ! Rotational  part of Kinetic factor G
    real(ark),intent(out) ::  g_cor(trove%Nmodes,3)        ! Coriolis part of Kinetic factor G
    real(ark),intent(out) ::  pseudo                  ! Coriolis part of Kinetic factor G

    real(ark)            :: s_vib(trove%Nmodes,trove%Natoms,3,0:trove%Nmodes,0:trove%Nmodes)
    real(ark)            :: s_rot(3,trove%Natoms,3,0:trove%Nmodes,0:trove%Nmodes)
    real(ark)            :: d_srot_dsvib(trove%Nmodes*trove%Natoms*3+3*trove%Natoms*3)

    real(ark)            :: d1s_vib(trove%Nmodes,trove%Natoms,3,trove%Nmodes)
    real(ark)            :: d2s_vib(trove%Nmodes,trove%Natoms,3,trove%Nmodes,trove%Nmodes)
    real(ark)            :: d1s_rot(3,trove%Natoms,3,trove%Nmodes)


    real(ark) :: pseudo_t(4),astep(trove%Nmodes,2)

    real(ark)    :: masses_i(1:trove%Natoms)
    real(ark)    :: factor
    real(ark)   :: chi_eq(trove%Nmodes)
    integer(ik) :: q1,q2,n1,x1,x2,x0,powers(trove%Nmodes)
    integer(ik) :: Natoms,Nmodes,i,isize,imode,jmode
    character(len=cl)  :: job_is
    !
    Natoms = trove%Natoms
    Nmodes = trove%Nmodes
    !
    ! Conversion factor to the cm-1 units 
    !
    factor = real(planck,ark)*real(avogno,ark)*real(1.0d+16,ark)/(4.0_ark*pi*pi*real(vellgt,ark))
    !
    ! Masses for easy reference 
    masses_i = 1.0_ark/trove%mass
    !
    !dchi = chi - trove%chi_ref(:,irho)
    !
    call s_vib_s_rot_dvr(dchi,irho,s_vib,s_rot)
    !
    g_vib = 0
    g_rot = 0 
    g_cor = 0 
    pseudo_t = 0
    !
    ! All kinetic parts are obtained as scalar products of two polynoms: 
    ! G_vib(q1,q2)  = -1/2*factor*sum_{n1} s_vib * s_vib /mass(n1)
    ! G_rot(q1,q2)  = -1/2*factor*sum_{n1} s_rot * s_rot /mass(n1)
    ! C_cor(q1,q2)  = -1/2*factor*sum_{n1} s_vib * s_rot /mass(n1)
    ! Only pseudopotential function (part 1) is obtained as a vector product  of two s_rot polynoms: 
    ! V_pseudi1(q1,q2) = 1/4*factor*sum_{n1} [s_rot x s_rot] /mass(n1)
    !
    !
    ! Vibrational part 
    !
    !omp parallel do private(q1,q2,n1,x1) shared(g_vib) schedule(dynamic)
    do q1 = 1,Nmodes
       do q2 = 1,Nmodes
          do n1 = 1,Natoms
             do x1 = 1,3
                !
                g_vib(q1,q2) = g_vib(q1,q2)+masses_i(n1)*s_vib(q1,n1,x1,0,0)*s_vib(q2,n1,x1,0,0)
                !
             enddo
          enddo 
       enddo
    enddo
    !omp end parallel do
    !
    g_vib = factor*g_vib
    !
    if (FLrotation) then
       ! 
       !write (out,"('Please make sure that you are ready to run the rotational part')") alloc
       !stop 'Is the rotational kinetic part ready?'
       !
       ! Rotational part 
       !
       !omp parallel do private(q1,q2,n1,x1) shared(g_rot) schedule(dynamic)
       do q1 = 1,3
          do q2 = 1,3
             do n1 = 1,Natoms
                do x1 = 1,3
                   !
                   g_rot(q1,q2) =g_rot(q1,q2)+masses_i(n1)*s_rot(q1,n1,x1,0,0)*s_rot(q2,n1,x1,0,0)
                   !
                enddo
             enddo 
          enddo
       enddo
       !omp end parallel do 
       !
       g_rot = factor*g_rot
       !
       ! Coriolis part 
       !
       !omp parallel do private(q1,q2,n1,x1) shared(g_cor) schedule(dynamic)
       do q1 = 1,Nmodes
          do q2 = 1,3
             do n1 = 1,Natoms
                do x1 = 1,3
                   !
                   g_cor(q1,q2) =g_cor(q1,q2)+masses_i(n1)*s_vib(q1,n1,x1,0,0)*s_rot(q2,n1,x1,0,0)
                   !
                enddo
             enddo 
          enddo
       enddo
       !omp end parallel do 
       !
       g_cor = factor*g_cor
       !
    endif
    !
    ! Pseudopotential function: part 1  
    !
    !omp parallel do private(n1,x0,x1,x2,q1,q2) shared(pseudo_t) schedule(dynamic)
    do n1 = 1,Natoms
       do x0 = 1,3
          do x1 = 1,3
             do x2 = 1,3
                do q1 = 1,3
                   do q2 = 1,3
                     !
                     pseudo_t(1) =pseudo_t(1)+0.125_ark*masses_i(n1)*s_rot(q1,n1,x1,0,0)*s_rot(q2,n1,x2,0,0)* &
                                              epsil(x0,q1,x1)*epsil(x0,q2,x2)
                   enddo
                enddo 
             enddo
          enddo 
       enddo
    enddo
    !omp end parallel do 
    !
    !
    !U2:=simplify(
    ! 1/4*sum(sum(add(add(add(
    ! epsilon[x0,g0,y0]*S[k0][N0,y0]*Ds[g0][N0,x0][k0]/mm[N0]
    !    ,x0=XYZ),y0=XYZ),g0=XYZ),k0=1..6),N0=1..4));
    !
    ! Pseudopotential fucntion: part 2
    !
    !omp parallel do private(n1,x0,q1,q2,x1) shared(pseudo_t) schedule(dynamic)
    do n1 = 1,Natoms
       do x0 = 1,3
          do q1 = 1,3
             do q2 = 1,Nmodes
                !
                do x1 = 1,3
                   ! +
                   pseudo_t(2) =pseudo_t(2)+ 0.25_ark*s_rot(q1,n1,x0,q2,0)*s_vib(q2,n1,x1,0,0)*epsil(x0,q1,x1)*masses_i(n1)
                   !
                enddo 
            enddo
         enddo 
      enddo
    enddo
    !omp end parallel do 
    !
    !U3:=simplify(sum(sum(sum(add(
    !> 1/4*1/mm[N0]*( S[q01][N0,x0]*DDS[q02][N0,x0][q02,q01]+1/2*DS[q01][N0,x0][q01]*DS[q02][N0,x0][q02] )
    !>   ,x0=XYZ),N0=1..4),q01=1..6),q02=1..6 )):
    !
    ! Pseudopotential fucntion: part 3a
    !
    !omp parallel do private(n1,x0,q1,q2) shared(pseudo_t) schedule(dynamic)
    do n1 = 1,Natoms
       do x0 = 1,3
          do q1 = 1,Nmodes
             do q2 = 1,Nmodes
                !
                pseudo_t(3) = pseudo_t(3)-0.25_ark*s_vib(q1,n1,x0,0,0)*s_vib(q2,n1,x0,q1,q2)*masses_i(n1)
                !
             enddo
          enddo 
       enddo
    enddo
    !omp end parallel do 
    !
    ! Pseudopotential fucntion: part 3b
    !
    !omp parallel do private(n1,x0,q1,q2) shared(pseudo_t) schedule(dynamic)
    do n1 = 1,Natoms
       do x0 = 1,3
          do q1 = 1,Nmodes
             do q2 = 1,Nmodes
                !
                pseudo_t(4) =pseudo_t(4)-0.125_ark*s_vib(q1,n1,x0,q1,0)*s_vib(q2,n1,x0,q2,0)*masses_i(n1)
                !
             enddo
          enddo 
       enddo
    enddo
    !omp end parallel do 
    !
    pseudo = ( pseudo_t(1)  + pseudo_t(2) + pseudo_t(3) +pseudo_t(4) )*factor
    !
    !if (job%verbose>=7) write(out,"(i,4g18.8)") irho,pseudo_t(1:4)
    !
    !
  end subroutine FLDVR_gmat_dvr




  function s_vib_s_rot_dvr_vec(chi,irho,Nsize) result (f)
    !
    real(ark),intent(in)    :: chi(trove%Nmodes)
    integer(ik),intent(in) :: irho
    integer(ik),intent(in)    :: Nsize

    real(rk),dimension(Nsize) :: f

    real(rk)               :: s_vib(trove%Nmodes,trove%Natoms,3)
    real(rk)               :: s_rot(3,trove%Natoms,3)
    !
    integer(ik)            :: i,q1,n1,x1
      !
      !
      !call s_vib_s_rot_dvr(chi,irho,s_vib,s_rot)
      !
      i = 0
      !
      do q1 = 1,3
         do n1 = 1,trove%Natoms
             do x1 = 1,3
                !
                i = i + 1
                !
                f(i) = s_rot(q1,n1,x1)
                !
             enddo
         enddo
      enddo
      !
      do q1 = 1,trove%Nmodes
         do n1 = 1,trove%Natoms
             do x1 = 1,3
                !
                i = i + 1
                !
                f(i) = s_vib(q1,n1,x1)
                !
             enddo
         enddo
      enddo
    !
  end function s_vib_s_rot_dvr_vec


  !
  ! This is a part of s_vib_s_rot_polynom1d calcualtions, which 
  ! does the actual recursive solving of the linear equations. 
  !
  recursive subroutine s_vib_s_rot_solve(icase,qmax,maxorder,iNcoeff,b0,r_na_xi,dr_na_dq,s_mat)
    ! 
    integer(ik),intent(in)      :: icase,qmax,maxorder,iNcoeff
    !
    real(ark),intent(in) :: b0(:,:),r_na_xi(:,:,:)
    real(ark),intent(in) :: dr_na_dq(:,:,:,:)
    real(ark),intent(out) :: s_mat(:,:,:,:)
    !
    integer(ik)            :: Nmodes,linmode,imode,Nequat,dimen,alloc,numeq,x1,i1,i2
    integer(ik)            :: numvar,n0,x0,y0,jmode,n,n1,qt,b_t_size,int,iNcoeff_t
    real(ark),allocatable   :: bm(:),Tmat(:,:),b(:,:),a(:,:)
    real(ark),allocatable  ::  db(:,:),da(:,:)
    !
    real(ark),allocatable   :: b_t(:,:)
    real(ark),allocatable   :: x_1t(:,:),x_2t(:,:),x_3t(:,:)
    real(ark),allocatable   :: s_t(:),work(:)
    integer(ik)            :: rank,iwork,info
       !
       Nmodes = trove%Nmodes
       Nequat = 6+Nmodes-min(1,trove%lincoord)
       dimen = Nequat
       iwork = 50*dimen
       b_t_size = 3+3+trove%Nmodes-min(trove%lincoord,1)
       !
       allocate (Tmat(dimen,dimen),Bm(dimen),b(dimen,1),a(dimen,dimen),& 
                 s_t(dimen),work(iwork),db(dimen,1),da(dimen,dimen),stat=alloc)
       if (alloc/=0) then
            write (out,"(' Error ',i9,' trying to allocate array for Tmat,Bm,a,b')") alloc
            stop 's_vib_s_rot_polynom1d, Bm ant Tmat - out of memory'
       end if
       !
       allocate(b_t(b_t_size,iNcoeff),stat=alloc)
       !
       allocate (x_1t(iNcoeff,0:0),stat=alloc)
       allocate (x_2t(iNcoeff,0:0),stat=alloc)
       allocate (x_3t(iNcoeff,0:0),stat=alloc)
       !
       do x1 = 1,3
          if (x1/=trove%lincoord) then 
             !
             numeq = numeq +1 
             do n1 = 1,trove%Natoms
                do x0 = 1,3
                   do y0 = 1,3
                      !
                      x_2t(1:iNcoeff,0) = r_na_xi(n1,y0,1:iNcoeff)
                      !
                   enddo
                enddo
             enddo
             !
          endif
       enddo
       !
       x_1t = 0 ; x_2t = 0 ; x_3t = 0
       !
       !allocate(b_t(b_t_size),s_t,s_1t,x_1t,x_2t,x_3t,stat=alloc)
       !
       s_mat = 0 
       !
       select case (icase) 
       case(1) 
         linmode=0 
         if (qmax /= Nmodes) then 
           write(out,"('Illegal qmax: ',2i8)") qmax,Nmodes
           stop 'Illegal qmax'
         endif 
       case(2) 
         linmode=trove%lincoord 
         if (qmax /= 3) then 
           write(out,"('Illegal qmax: ',2i8)") qmax,3
           stop 'Illegal qmax'
         endif 
       end select 
       !
       ! Run the loop over different modes
       !
       do imode = 1,qmax
          !
          if (imode/=linmode) then 
             !
             !if (verbose>=2) write(out,"('imode = ',i5)") imode
             !
             tmat = 0
             bm = 0

             ! a) translational part (3 members: x,y,z)
             !
             numeq = 0
             do x1 = 1,3
                numeq = numeq +1 
                numvar= 0
                do n0 = 1,trove%Natoms
                   do x0 = 1,3
                      numvar= numvar+1
                      if (x0==x1) then 
                         tmat(numeq,numvar) = 1.0_ark
                      endif 
                   enddo
                enddo
             enddo

             ! b) rotational part (3 members: x,y,z)       
             !
             do x1 = 1,3
                if (x1/=trove%lincoord) then 
                   numeq = numeq +1 
                   if (icase==2.and.imode==x1) bm(numeq) = 1.0_ark
                   !
                   numvar= 0
                   do n0 = 1,trove%Natoms
                      do x0 = 1,3
                         numvar= numvar+1
                         tmat(numeq,numvar) = sum(epsil(x1,:,x0)*b0(n0,:))
                      enddo
                   enddo
                endif 
             enddo

             ! c) Vibrational part (Nmodes members)
             !
             do jmode = 1,Nmodes
               !
               numeq = numeq +1 
               if (icase==1.and.imode==jmode) bm(numeq) = 1.0_ark
               !
               numvar= 0
               do n0 = 1,trove%Natoms
                  do x0 = 1,3
                     numvar= numvar+1
                     tmat(numeq,numvar) = dr_na_dq(n0,x0,jmode,1)  ! amat(n0,x0,jmode)
                  enddo
               enddo
               !
             enddo
             !
             ! -------------------------------------------
             ! Solve linear equation with Lapack procedure
             ! -------------------------------------------
             !
             if (dimen/=numeq.or.dimen/=numvar) then 
                write(out,"('s_vib_s_rot_polynom1d: the size of Tmat contradicts the number of eqs or vars:',3i7)") & 
                             dimen,numeq,numvar
                write(out,"('                     icase,imode = :',2i7)") icase,imode
                stop 's_vib_s_rot_polynom1d: ieq is inconsistent'
             endif 
             !
             da = tmat
             db(:,1) = bm(:)
             !
             ! lapack_gelss - solves a linear equation by least squares method 
             ! 
             call dgelss(dimen,dimen,1,da(:,:),dimen,db(:,1),dimen, &
                         s_t,-1.0d-12, rank, work, iwork, info)
             !
             if (info/=0) then
                write (out,"('s_vib_s_rot_polynom1d: dgelss returned ',i8)") info
                stop 's_vib_s_rot_polynom1d - dgelss failed'
             end if
             !
             ! call lapack_gelss(a(:,:),b(:,:))
             !
             bm(:) = db(:,1)
             !
             ! check if the solution is not trivial 
             !
             if (sum(bm(:)**2)<small_ ) then 
                write(out,"(/'s_vib_s_rot_polynom1d: the trivial solition found for imode: ',i7)") imode
                write(out,"('                       this point will make a naught in basis functions.'/)") 
                !
                ! 
                ! if this happens we do not give up, instead 
                ! we make ensure a naught for the basis functions at this point
                !
                trove%sing_at_rho_0 = .true.
                !
                !stop 's_vib_s_rot_polynom1d: ieq is inconsistent'
                !
             endif 
             !
             ! The solution written in bm -> s_vib_0
             !
             numvar= 0
             do n0 = 1,trove%Natoms
                do x0 = 1,3
                   !
                   numvar= numvar+1
                   s_mat(imode,n0,x0,1) = bm(numvar)
                   !
                enddo
             enddo
             !
             ! Higher order solutions s_vib_n obtained with the same Tmat and bm_n
             ! bm_n_x1 = sum_{x0,y0} sum_{0}^{n-1} epsil(x0,x1,y0)*s_vib(n1,x0,l)*r_na(n1,y0,n-l)
             !
             ! bm_n is not zero only at the rotational part b) unless it is a local coordinates transformation 
             !
             do n  = 1,maxorder
                !
                b_t = 0
                bm  = 0
                numeq = 3

                ! b) rotational part (3 members: x,y,z)       
                !
                do x1 = 1,3
                   if (x1/=trove%lincoord) then 
                      !
                      numeq = numeq +1 
                      do n1 = 1,trove%Natoms
                         do x0 = 1,3
                            do y0 = 1,3
                               ! 
                               x_1t(1:iNcoeff,0) = s_mat(imode,n1,x0,1:iNcoeff)
                               !
                               x_2t(1:iNcoeff,0) = r_na_xi(n1,y0,1:iNcoeff)
                               !
                               call product_of_polynoms_simple(iNcoeff,x_1t,iNcoeff,x_2t,iNcoeff_t,x_3t,n)
                               !
                               b_t(numeq,:) = b_t(numeq,:) + x_3t(:,0)*epsil(x1,y0,x0)
                               ! sign +/- 
                               !
                            enddo
                         enddo
                      enddo
                      !
                   endif
                enddo
                !
                ! c) Vibrational part (Nmodes members)
                !
                !if (trove%internal_coords=='LOCAL') then 
                !
                do jmode = 1,Nmodes
                  numeq = numeq +1 
                  !
                  do n1 = 1,trove%Natoms
                     do x0 = 1,3
                        !
                        x_1t(1:iNcoeff,0) = s_mat(imode,n1,x0,1:iNcoeff)
                        !
                        x_2t(1:iNcoeff,0) = dr_na_dq(n1,x0,jmode,1:iNcoeff)
                        !
                        call product_of_polynoms_simple(iNcoeff,x_1t,iNcoeff,x_2t,iNcoeff_t,x_3t,n)
                        !
                        b_t(numeq,:) = b_t(numeq,:) + x_3t(:,0)
                        !
                     enddo
                  enddo
                enddo
                !
                ! Solving the same equation Tmat s = b with new b_n at every term of xi1^k1 xi2^k2 ..., k1+k2+k3+..k_modes = n 
                !
                do int = trove%RangeOrder(n-1)+1,trove%RangeOrder(n)
                   do qt = 1,b_t_size
                      bm(qt) =-b_t(qt,int)
                      !
                   enddo 
                   !
                   !  Solve it
                   !
                   a = tmat
                   b(:,1) = bm(:)
                   !
                   ! gelss - solves a linera equation by least squares method 
                   ! 
                   call dgelss(dimen,dimen,1,a(:,:),dimen,b(:,1),dimen, &
                        s_t,-1.0d-12, rank, work, iwork, info)

                   !
                   if (info/=0) then
                     write (out,"('s_vib_s_rot_polynom1d: dgelss returned ',i8)") info
                     stop 's_vib_s_rot_polynom1d - dgelss failed'
                   end if
                   !
                   !call lapack_gelss(a(:,:),b(:,:))
                   !
                   bm(:) = b(:,1)
                   !
                   numvar= 0
                   do n0 = 1,trove%Natoms
                      do x0 = 1,3
                         !
                         numvar= numvar+1
                         s_mat(imode,n0,x0,int) = bm(numvar)
                         !
                      enddo
                   enddo
                enddo
             enddo
             !
          endif 
          !
       enddo
       !
       deallocate(Tmat,bm,a,b,da,db)
       deallocate(s_t,work,x_1t,x_2t,x_3t,b_t)
       !
       !call ArrayStop('dr_na_dq')
       !
     end subroutine s_vib_s_rot_solve

  !
  ! Here we initialize the kinetic operator fields: g_vib, g_rot, and g_cor
  subroutine gmat_polynom(s_vib,s_rot,g_vib,g_rot,g_cor,pseudo)


    type(FLpolynomT),intent(in) :: s_vib(trove%Nmodes,trove%Natoms,3)
    type(FLpolynomT),intent(in) :: s_rot(3,trove%Natoms,3)

    type(FLpolynomT),intent(inout) ::  g_vib(trove%Nmodes,trove%Nmodes)    ! Vibrational part of Kinetic factor G
    type(FLpolynomT),intent(inout) ::  g_rot(3,3)              ! Rotational  part of Kinetic factor G
    type(FLpolynomT),intent(inout) ::  g_cor(trove%Nmodes,3)        ! Coriolis part of Kinetic factor G
    type(FLpolynomT),intent(inout) ::  pseudo                  ! Coriolis part of Kinetic factor G

    real(ark)    :: masses(1:trove%Natoms),rhostep,f_t,df_t
    real(ark),allocatable :: s_1t(:,:),s_2t(:,:),s_3t(:,:),s_4t(:,:),s_5t(:,:)
    type(FLpolynomT),pointer :: fl
    real(ark)    ::  factor
    real(ark),allocatable   ::  r_t(:),pseudo_t(:,:,:)
    integer(ik) :: q1,q2,n1,x1,x2,x0,iterm,gN,gN1,gN2,gN_t,gO,kindex(trove%Nmodes)
    integer(ik) :: Natoms,Nmodes,Npoints,i0,alloc
    !
    if (job%verbose>=2) write(out,"(/'gmat_polynom/start')") 
    !

    Natoms = trove%Natoms
    Nmodes = trove%Nmodes
    Npoints = trove%Npoints
    !
    ! Here we go! 
    !
    !
    ! Conversion factor to the cm-1 units 
    !
    factor = real(planck,ark)*real(avogno,ark)*real(1.0d+16,kind=ark)/(4.0_ark*pi*pi*real(vellgt,ark))
    !
    !factor = 1.0_rk
    !
    ! Masses for easy reference 
    masses = trove%mass
    !
    ! How many coeffs in g_xxx fields
    !
    gO = g_vib(1,1)%Orders
    gN = g_vib(1,1)%Ncoeff
    gN1 = trove%RangeOrder(max(gO+1,0))
    gN2 = trove%RangeOrder(max(gO+2,0))
    !
    rhostep = trove%rhostep
    allocate(pseudo_t(4,gN,0:Npoints),stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i9,' trying to allocate pseudo1 arrays')") alloc
       stop 'gmat_polynom: s_t - out of memory'
    endif
    !
    pseudo_t = 0
    !
    ! All kinetic parts are obtained as scalar products of two polynoms: 
    ! G_vib(q1,q2)  = -1/2*factor*sum_{n1} s_vib * s_vib /mass(n1)
    ! G_rot(q1,q2)  = -1/2*factor*sum_{n1} s_rot * s_rot /mass(n1)
    ! C_cor(q1,q2)  = -1/2*factor*sum_{n1} s_vib * s_rot /mass(n1)
    ! Only pseudopotential function (part 1) is obtained as a vector product  of two s_rot polynoms: 
    ! V_pseudi1(q1,q2) = 1/4*factor*sum_{n1} [s_rot x s_rot] /mass(n1)
    !
    !$omp sections private(s_1t,s_2t,s_3t,s_4t,s_5t,r_t,kindex)
    !
    ! Allocating two temporaly arrays s_1t, s2t and s_3t
    !
    allocate(s_1t(gN,0:Npoints),s_2t(gN,0:Npoints),s_3t(gN,0:Npoints),&
             s_4t(gN1,0:Npoints),s_5t(gN2,0:Npoints),r_t(0:Npoints),stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i9,' trying to allocate s_t arrays')") alloc
       stop 'gmat_polynom: s_t - out of memory'
    endif 
    !
    !$omp section
    !
    ! Vibrational part 
    !
    do q1 = 1,Nmodes
       do q2 = 1,Nmodes
          do n1 = 1,Natoms
             do x1 = 1,3
                !
                s_1t(1:gN,:) = s_vib(q1,n1,x1)%field(1:gN,:)
                s_2t(1:gN,:) = s_vib(q2,n1,x1)%field(1:gN,:)
                !
                call product_of_polynoms_simple(gN,s_1t,gN,s_2t,gN_t,s_3t)
                !
                g_vib(q1,q2)%field(1:gN,:) = g_vib(q1,q2)%field(1:gN,:)+factor/masses(n1)*s_3t(1:gN,:)
                !
             enddo
          enddo 
       enddo
    enddo
    !
    !$omp section
    !
    if (job%verbose>=4) write(out,"('g_vib... Done!')") 
    !
    if (FLrotation) then
       ! 
       !write (out,"('Please make sure that you are ready to run the rotational part')") alloc
       !stop 'Is the rotational kinetic part ready?'
       !
       ! Rotational part 
       !
       do q1 = 1,3
          do q2 = 1,3
             do n1 = 1,Natoms
                do x1 = 1,3
                   !
                   s_1t(1:gN,:) = s_rot(q1,n1,x1)%field(1:gN,:)
                   s_2t(1:gN,:) = s_rot(q2,n1,x1)%field(1:gN,:)
                   !
                   call product_of_polynoms_simple(gN,s_1t,gN,s_2t,gN_t,s_3t)
                   !
                   g_rot(q1,q2)%field(1:gN,:) =g_rot(q1,q2)%field(1:gN,:)+factor/masses(n1)*s_3t(1:gN,:)
                   !
                enddo
             enddo 
          enddo
       enddo
       !
       if (job%verbose>=4) write(out,"('g_rot... Done!')") 
       !
       ! Coriolis part 
       !
       do q1 = 1,Nmodes
          do q2 = 1,3
             do n1 = 1,Natoms
                do x1 = 1,3
                   !
                   s_1t(1:gN,:) = s_vib(q1,n1,x1)%field(1:gN,:)
                   s_2t(1:gN,:) = s_rot(q2,n1,x1)%field(1:gN,:)
                   !
                   call product_of_polynoms_simple(gN,s_1t,gN,s_2t,gN_t,s_3t)
                   !
                   g_cor(q1,q2)%field(1:gN,:) =g_cor(q1,q2)%field(1:gN,:)+factor/masses(n1)*s_3t(1:gN,:)
                   !
                enddo
             enddo 
          enddo
       enddo
       !
       if (job%verbose>=4) write(out,"('g_cor... Done!')") 
       !
    endif
    !
    !$omp section
    !
    ! Pseudopotential function: part 1  
    !
    do n1 = 1,Natoms
       do x0 = 1,3
          do x1 = 1,3
             do x2 = 1,3
                do q1 = 1,3
                   do q2 = 1,3
                     !
                     s_1t(1:gN,:) = s_rot(q1,n1,x1)%field(1:gN,:)
                     s_2t(1:gN,:) = s_rot(q2,n1,x2)%field(1:gN,:)
                     !
                     call product_of_polynoms_simple(gN,s_1t,gN,s_2t,gN_t,s_3t)
                     !
                     pseudo_t(1,1:gN,:) =pseudo_t(1,1:gN,:)+0.125_ark*factor*s_3t(1:gN,:)* &
                                                                epsil(x0,q1,x1)*epsil(x0,q2,x2)/masses(n1)
                   enddo
                enddo 
             enddo
          enddo 
       enddo
    enddo
    !
    !fl => trove%pseudo
    !
    if (job%verbose>=4) write(out,"('pseudo1... Done!')") 
    !
    !$omp section
    !
    !U2:=simplify(
    ! 1/4*sum(sum(add(add(add(
    ! epsilon[x0,g0,y0]*S[k0][N0,y0]*Ds[g0][N0,x0][k0]/mm[N0]
    !    ,x0=XYZ),y0=XYZ),g0=XYZ),k0=1..6),N0=1..4));
    !
    ! Pseudopotential fucntion: part 2
    !
    do n1 = 1,Natoms
       do x0 = 1,3
          do q1 = 1,3
             do q2 = 1,Nmodes
                !
                if (q2==trove%Nmodes_n) then 
                   !
                   do i0 = 1,gN
                      !
                      call diff_2d_4points_ark(Npoints,trove%rho_border,&
                                           s_rot(q1,n1,x0)%field(i0,0:npoints),trove%periodic,0_ik,s_1t(i0,0:npoints))
                      !
                   enddo
                   !
                else
                   !
                   s_4t(1:gN1,:) = s_rot(q1,n1,x0)%field(1:gN1,:)
                   !
                   kindex = 0 ; kindex(q2) = 1
                   !
                   call deriv_of_polynoms_simple(gO+1,gN1,s_4t,kindex,gN,s_1t)
                   !
                endif
                !
                do x1 = 1,3
                   !
                   s_2t(1:gN,:) = s_vib(q2,n1,x1)%field(1:gN,:)
                   !
                   call product_of_polynoms_simple(gN,s_2t,gN,s_1t,gN_t,s_3t)
                   ! +
                   pseudo_t(2,1:gN,:) =pseudo_t(2,1:gN,:)+ & 
                                                0.25_ark*factor*s_3t(1:gN,:)*epsil(x0,q1,x1)/masses(n1)
                enddo 
            enddo
         enddo 
      enddo
    enddo
    !
    if (job%verbose>=4) write(out,"('pseudo2... Done!')") 
    !
    !$omp section
    !
    !U3:=simplify(sum(sum(sum(add(
    !> 1/4*1/mm[N0]*( S[q01][N0,x0]*DDS[q02][N0,x0][q02,q01]+1/2*DS[q01][N0,x0][q01]*DS[q02][N0,x0][q02] )
    !>   ,x0=XYZ),N0=1..4),q01=1..6),q02=1..6 )):
    !
    ! Pseudopotential fucntion: part 3a
    !
    do n1 = 1,Natoms
       do x0 = 1,3
          do q1 = 1,Nmodes
             do q2 = 1,Nmodes
                !
                if (q1==trove%Nmodes_n.and.q2==trove%Nmodes_n) then 
                   !
                   do i0 = 1,gN
                      !
                      call diff_2d_4points_ark(Npoints,trove%rho_border,&
                                s_vib(q2,n1,x0)%field(i0,0:npoints),trove%periodic,0_ik,r_t(0:npoints),s_2t(i0,0:npoints))
                      !
                   enddo
                   !
                elseif (q1==trove%Nmodes_n.or.q2==trove%Nmodes_n) then
                   !
                   s_4t(1:gN1,:) = s_vib(q2,n1,x0)%field(1:gN1,:)
                   !
                   kindex = 0 ; kindex( min(q1,q2) ) = 1 
                   !
                   call deriv_of_polynoms_simple(gO+1,gN1,s_4t,kindex,gN,s_1t)
                   !
                   do i0 = 1,gN
                      !
                      call diff_2d_4points_ark(Npoints,trove%rho_border,s_1t(i0,0:npoints),trove%periodic,0_ik,s_2t(i0,0:npoints))
                      !
                   enddo
                   !
                else
                   !
                   s_5t(1:gN2,:) = s_vib(q2,n1,x0)%field(1:gN2,:)
                   !
                   kindex = 0 ; kindex(q1) = 1 ; kindex(q2) = kindex(q2) + 1
                   !
                   call deriv_of_polynoms_simple(gO+2,gN2,s_5t,kindex,gN,s_2t)
                   !
                endif
                !
                s_1t(1:gN,:) = s_vib(q1,n1,x0)%field(1:gN,:)
                !
                call product_of_polynoms_simple(gN,s_2t,gN,s_1t,gN_t,s_3t)
                !
                pseudo_t(3,1:gN,:) = pseudo_t(3,1:gN,:)-0.25_ark*factor*s_3t(1:gN,:)/masses(n1)
                !
             enddo
          enddo 
       enddo
    enddo
    !
    !$omp section
    !
    if (job%verbose>=4) write(out,"('pseudo3a... Done!')") 
    !
    ! Pseudopotential fucntion: part 3b
    !
    do n1 = 1,Natoms
       do x0 = 1,3
          do q1 = 1,Nmodes
             !
             if (q1==trove%Nmodes_n) then 
                !
                do i0 = 1,gN
                   !
                   call diff_2d_4points_ark(Npoints,trove%rho_border,&
                                        s_vib(q1,n1,x0)%field(i0,0:npoints),trove%periodic,0_ik,s_1t(i0,0:npoints))
                   !
                enddo
                !
             else
                !
                s_4t(1:gN1,:) = s_vib(q1,n1,x0)%field(1:gN1,:)
                !
                kindex = 0 ; kindex(q1) = 1
                !
                call deriv_of_polynoms_simple(gO+1,gN1,s_4t,kindex,gN,s_1t)
                !
             endif
             !
             do q2 = 1,Nmodes
                !
                if (q2==trove%Nmodes_n) then 
                   !
                   do i0 = 1,gN
                      !
                       call diff_2d_4points_ark(Npoints,trove%rho_border,&
                                            s_vib(q2,n1,x0)%field(i0,0:npoints),trove%periodic,0_ik,s_2t(i0,0:npoints))
                      !
                   enddo
                   !
                else
                   !
                   s_4t(1:gN1,:) = s_vib(q2,n1,x0)%field(1:gN1,:)
                   !
                   kindex = 0 ; kindex(q2) = 1 
                   !
                   call deriv_of_polynoms_simple(gO+1,gN1,s_4t,kindex,gN,s_2t)
                   !
                endif
                !
                call product_of_polynoms_simple(gN,s_1t,gN,s_2t,gN_t,s_3t)
                !
                pseudo_t(4,1:gN,:) =pseudo_t(4,1:gN,:)-0.125_ark*factor*s_3t(1:gN,:)/masses(n1)
                !
             enddo
          enddo 
       enddo
    enddo
    !
    if (job%verbose>=4) write(out,"('pseudo4a... Done!')") 
    !
    pseudo%field(:,:) = pseudo_t(1,:,:)  + pseudo_t(2,:,:) + pseudo_t(3,:,:) +pseudo_t(4,:,:)
    !
    write(out,"(/'pseudo')") 
    do i0 = 0,Npoints
      !
      !if (i0<34) then
      !  pseudo%field(1:gN,i0) = -0.125_ark*g_cor(1,1)%field(1:gN,i0)
      !endif 
      !
      if (job%verbose>=4.and.trim(molec%moltype)/='XY') then 
        !
        write(out,"(i8,8g16.8)") i0,pseudo_t(1,1,i0),pseudo_t(2,1,i0),pseudo_t(3,1,i0),pseudo_t(4,1,i0),&
              g_rot(1,1)%field(1,i0),g_rot(2,2)%field(1,i0),g_rot(3,3)%field(1,i0),g_vib(3,3)%field(1,i0)
        !
      endif
    enddo
    !
    !pseudo%field(:,:) =  pseudo_t(1,:,:) +pseudo_t(2,:,:)
    !
    if (trove%lincoord/=0) then 
      !
      !
      !pseudo%field = 0
      !
      !pseudo%field(1:gN,:) = -0.125_ark*g_cor(2,2)%field(1:gN,:)
      !
      !do x0=1,3 
      !  !
      !  if (x0/=trove%lincoord) then 
      !    !
      !    pseudo%field(1:gN,:) = pseudo%field(1:gN,:)-0.5_ark*0.125_ark*g_cor(x0,x0)%field(1:gN,:)
      !    !
      !  endif 
      !  !
      !enddo
      !
    endif 
    !
    !
    deallocate(r_t)
    !
    deallocate(s_1t,s_2t,s_3t,s_4t,s_5t)
    !$omp end sections 
    !
    !
    if (trim(molec%coords_transform)=='R-RHO'.and..false.) then 
      !
      allocate(s_1t(gN,0:Npoints),s_2t(gN,0:Npoints),s_3t(gN,0:Npoints),&
               s_4t(gN1,0:Npoints),s_5t(gN2,0:Npoints),r_t(0:Npoints),stat=alloc)
      if (alloc/=0) then
         write (out,"(' Error ',i9,' trying to allocate s_t arrays')") alloc
         stop 'gmat_polynom: s_t - out of memory'
      endif 
      !
      !trove%sing_at_rho_0(1) = 1
      !
      ! Pseudopotential fucntion: part 1-2sing
      !
      pseudo_t = 0 
      !
      do n1 = 1,Natoms
         do x0 = 1,3
            do q1 = 1,Nmodes
               !
               if (q1==trove%Nmodes_n) then 
                  !
                  do i0 = 1,gN
                     !
                     call diff_2d_4points_ark(Npoints,trove%rho_border,&
                                          s_vib(q1,n1,x0)%field(i0,0:npoints),trove%periodic,0_ik,s_1t(i0,0:npoints))
                     !
                     call diff_2d_4points_ark(Npoints,trove%rho_border,&
                                          s_vib(trove%Nmodes_n,n1,x0)%field(i0,0:npoints),trove%periodic,0_ik,s_2t(i0,0:npoints))
                     !
                  enddo
                  !
               else
                  !
                  s_5t(1:gN1,:) = s_vib(q1,n1,x0)%field(1:gN1,:)
                  !
                  kindex = 0 ; kindex(q1) = 1
                  !
                  call deriv_of_polynoms_simple(gO+1,gN1,s_5t,kindex,gN,s_1t)
                  !
                  s_5t(1:gN1,:) = s_vib(trove%Nmodes_n,n1,x0)%field(1:gN1,:)
                  !
                  call deriv_of_polynoms_simple(gO+1,gN1,s_5t,kindex,gN,s_2t)
                  !
               endif
               !
               s_3t(1:gN,:) = s_vib(q1            ,n1,x0)%field(1:gN,:)
               s_4t(1:gN,:) = s_vib(trove%Nmodes_n,n1,x0)%field(1:gN,:)
               !
               call product_of_polynoms_simple(gN,s_1t,gN,s_4t,gN_t,s_5t)
               !
               call product_of_polynoms_simple(gN,s_2t,gN,s_3t,gN_t,s_1t)
               !
               pseudo_t(1,1:gN,:) = pseudo_t(1,1:gN,:)+factor*(s_1t(1:gN,:)+s_5t(1:gN,:))/masses(n1)
               !
            enddo 
         enddo
      enddo
      !
      s_1t(1,:) = sqrt(trove%imat_s(:)/factor)
      !
      call diff_2d_4points_ark(Npoints,trove%rho_border,s_1t(1,:),trove%periodic,0_ik,s_2t(1,:))
      call diff_2d_4points_ark(Npoints,trove%rho_border,s_1t(1,:),trove%periodic,0_ik,r_t(:),s_3t(1,:))
      !
      do i0 = 0,Npoints
        !
        r_t(i0) = (trove%rho_border(1)+rhostep*i0)
        !
      enddo
      !
      forall(i0 = 1:gN) pseudo_t(1,i0,0:Npoints) = pseudo_t(1,i0,0:Npoints)*s_2t(1,0:Npoints)*s_1t(1,0:Npoints)
      !
      forall(i0 = 1:gN) pseudo%field(i0,0:Npoints) = pseudo%field(i0,0:Npoints)*s_1t(1,0:Npoints)*s_1t(1,0:Npoints)
      !
      forall(i0 = 1:gN) pseudo_t(2,i0,0:Npoints) = pseudo_t(2,i0,0:Npoints)+&
                        g_vib(Nmodes,Nmodes)%field(i0,0:Npoints)*s_3t(1,0:Npoints)*s_1t(1,0:Npoints)
      !
      pseudo%field(:,:) =  pseudo%field(:,:) -0.5_ark*(pseudo_t(1,:,:) + pseudo_t(2,:,:))
      !
      do q1 = 1,3
         do q2 = 1,3
             !
             forall(i0 = 1:gN) g_rot(q1,q2)%field(i0,0:Npoints) = g_rot(q1,q2)%field(i0,0:Npoints)*s_1t(1,0:Npoints)*&
                                                                  s_1t(1,0:Npoints)
             !
         enddo
      enddo
      !
      do q1 = 1,Nmodes
         do q2 = 1,3
             !
             forall(i0 = 1:gN) g_cor(q1,q2)%field(i0,0:Npoints) = g_cor(q1,q2)%field(i0,0:Npoints)*s_1t(1,0:Npoints)*&
                                                                  s_1t(1,0:Npoints)
             !
         enddo
      enddo
      !
      do q1 = 1,Nmodes
         do q2 = 1,Nmodes
             !
             forall(i0 = 1:gN) g_vib(q1,q2)%field(i0,0:Npoints) = g_vib(q1,q2)%field(i0,0:Npoints)*s_1t(1,0:Npoints)*&
                                                                  s_1t(1,0:Npoints)
             !
         enddo
      enddo
      !
      !write(out,"(/'pseudo*rho: pseudo%field(1,i0),pseudo_t(1,1,i0),pseudo_t(2,1,i0),g_rot(1,1)%field(1,i0),g_rot(2,2)%field(1,i0),g_rot(3,3)%field(1,i0),g_vib(3,3)%field(1,i0)')") 
      do i0 = 0,Npoints
        !
        if (job%verbose>=4) then 
          !
          write(out,"(i8,8f22.8)") i0,r_t(i0),pseudo%field(1,i0),pseudo_t(1,1,i0),pseudo_t(2,1,i0),g_rot(1,1)%field(1,i0),&
                                   g_rot(2,2)%field(1,i0),g_rot(3,3)%field(1,i0),g_vib(3,3)%field(1,i0)
          !
        endif
      enddo
      !
    endif 
    !
    if (trove%sing_at_rho_0) then 
      !
      return
      !
      do i0 = 0,10
        !
        r_t(i0) = trove%rho_border(1)+rhostep*i0
        !
      enddo
      !
      do iterm = 1,1 ! iNcoeff 
        !
        do i0 = 0,4
          !
          call MLratintark(r_t(5:8),pseudo%field(iterm,5:8),r_t(i0),f_t,df_t)
          !
          pseudo%field(iterm,i0) = f_t
          !
        enddo
        !
      enddo
      !
      do q1 = 1,3
         do q2 = 1,3
           !
           do iterm = 1,1 ! iNcoeff 
             !
             do i0 = 0,4
               !
               call MLratintark(r_t(5:9),g_rot(q1,q2)%field(iterm,5:9),r_t(i0),f_t,df_t)
               !
               g_rot(q1,q2)%field(iterm,i0) = f_t
               !
             enddo
             !
           enddo 
           !
         enddo
      enddo
      !
      deallocate(s_1t,s_2t,s_3t,s_4t,s_5t)
      !
    endif 

    !pseudo%field = 0
    !
    deallocate(pseudo_t)
    !
    !call ArrayStop('s_1t')
    !call ArrayStop('s_2t')
    !call ArrayStop('s_3t')
    !call ArrayStop('s_4t')
    !call ArrayStop('s_5t')
    !
    if (job%verbose>=2) write(out,"('gmat_polynom/done!')") 
    !
  end subroutine gmat_polynom
  !
  ! Scalar product of two polinomial forms 
  ! A polynom is defifed by the its order (Norder) and a number of modes (Nmodes), and, of course, 
  ! by its coeffitients fields(:), given as a 1D matrix (i=1..Ncoeff)
  ! The correlation between indexes k1,k2,k3,...k_Nmodes -> i are determined from recursive function FLQindex
  ! The inverse transformation must be prestored in  FLIndexQ(l,i), where l=1,2,3,...,NModes
  ! 
   subroutine product_of_polynoms_fields(src1,src2,dst,level)

    type(FLpolynomT),intent(in)        :: src1,src2
    type(FLpolynomT),intent(out)       :: dst
    !integer(ik), intent(in),optional   :: Nrho
    integer(ik), intent(in), optional  :: level    ! We need only result at "level"=0..Norder
    !
    !integer                      :: irho1,irho2,irhodst
    !
    integer(ik) :: k1(trove%Nmodes),k2(trove%Nmodes),kdst(trove%Nmodes),n1,iopt,i1,i2,ndst,idst,Norder
    integer(ik) :: kindex(trove%Nmodes)
    integer(ik) :: index1,index2,sig,jmode
    !real(rk)    :: f_tmp
    !

      if (verbose>=7) write (out,"('Product of two polynoms ',a,' and ',a,' goes to ',a)") &
                                     trim(src1%name),trim(src2%name),trim(dst%name)
 
      !
      ! Define the order of the dst-polynom as the order of the smallest polynom 
      !
      !
      Norder = min(src1%Orders,src2%Orders)
      !
      ! debuging purposes
      !
      if (verbose>=7.and.size(dst%field)<trove%RangeOrder(Norder)) then 
        write (out,"('Error in product_of_polynoms_fields: size of dst-polynom is too small: ',i8)") size(dst%field)
        stop  'Error in product_of_polynoms_fields: order of dst-polynom is too small'
      endif
      !
      if (verbose>=7.and..not.allocated(FLIndexQ)) then 
        write (out,"('product_of_polynoms_fields: FLIndexQ has not been predefined ')")
        stop  'product_of_polynoms_fields: FLIndexQ must have been already defined!'
      endif


      dst%Orders = Norder
      dst%Ncoeff = min(src1%Ncoeff,src2%Ncoeff)


      dst%field = 0

      if (present(level)) then 

        if (verbose>=7.and.Norder<level-1) then 
            write(out,"('product_of_polynoms_fields-> level-1 is too big, bigger than Norder:',2i8)") level-1,Norder
            stop 'product_of_polynoms_fields-> level-1 is too big'
        endif 
        !
        iopt = 1
        !
      else 
        !
        iopt = 0
        !
      endif 
      !
      ! indexes of dst polynom
      !
      do idst = 1,dst%Ncoeff
         kdst(:) = FLIndexQ(:,idst)
         ndst = sum(kdst(:))
         !
         ! if level is presented we only need the fields for iorder = level
         ! besides the summation must be taken only up to level-1
         !
         if (iopt==0.or.ndst==level) then
            !
            ! indexes for summation: i2 = 1..Ncoeff, k2 = (k1,k2,k3,k4...), n2 = k1+k2+k3..
            !
            i1 = 0 
            k1 = 0
            n1 = 0 
            !f_tmp = 0.0_rk
            !
            do while(n1<=ndst-iopt.and.i1<src1%Ncoeff)
               i1 = i1 + 1
               k1(:) = FLIndexQ(:,i1)
               n1 = sum(k1(:))

               k2(:) = kdst(:)-k1(:)
               if (all(k2(:)>=0)) then 
                 ! i2 = FLQindex(trove%Nmodes_e,k2)

                 index1=1
                 index2=src2%Ncoeff
                 sig=2
                 !
                 !i2 = 0
                 !
                 do while(sig/=0)
                  !
                  !i2 = i2 + 1
                  !
                  i2=(index1+index2+sig)/2
                  !iterm = max(1,iterm) ; iterm = min(icoeff2,iterm)
                  kindex(:)=FLIndexQ(:,i2)
                  sig=0
                  !
                  if (sum(kindex)>sum(k2)) sig =-1
                  if (sum(kindex)<sum(k2)) sig = 1
                  !
                  if (sig==0) then 
                    !
                    !if (sig==0.and.any((kindex(1,trove%Nmodes_e)<k2(1,trove%Nmodes_e))) sig = 1
                    !if (sig==0.and.any((kindex(1,trove%Nmodes_e)>k2(1,trove%Nmodes_e))) sig =-1
                    !
                    do jmode=1,trove%Nmodes_e
                     !
                     if(kindex(jmode)<k2(jmode)) then 
                       sig= 1 
                     elseif(kindex(jmode)>k2(jmode)) then 
                       sig=-1
                     endif 
                     !
                     if (sig/=0) exit
                     !
                    enddo
                    !
                  endif 
                  !
                  if (sig== 1) index1=i2
                  if (sig==-1) index2=i2
                  !
                 enddo

                 if (i2<=src2%Ncoeff) then
                   ! 
                   !dst%field(idst,irhodst)=dst%field(idst,irhodst) + src1%field(i1,irho1)*src2%field(i2,irho2)
                   !
                   dst%field(idst,:)=dst%field(idst,:) + src1%field(i1,:)*src2%field(i2,:)
                   !
                 endif 
               endif 
            enddo
            !
         endif
      enddo

   end subroutine product_of_polynoms_fields



   subroutine product_of_polynoms_simple(icoeff1,src1,icoeff2,src2,icoeff,dst,level)

    integer(ik),intent(in)     :: icoeff1,icoeff2
    real(ark),intent(in)        :: src1(:,:),src2(:,:)
    integer(ik),intent(out)    :: icoeff
    real(ark),intent(out)       :: dst(:,:)
    !integer(ik), intent(in),optional   :: Nrho
    integer(ik), intent(in), optional  :: level    ! We need only result at "level"=0..Norder
    !
    integer(ik) :: k1(trove%Nmodes),k2(trove%Nmodes),kdst(trove%Nmodes),n1,iopt,i1,i2,ndst,idst,Norder
    integer(ik) :: kindex(trove%Nmodes)
    integer(ik) :: index1,index2,sig,iterm,jmode
    !

      !
      ! Define the order of the dst-polynom as the order of the smallest polynom 
      !
      icoeff = min(icoeff1,icoeff2)

      dst = 0

      if (present(level)) then 
        iopt = 1
      else 
        iopt = 0
      endif 
      !
      ! indexes of dst polynom
      !
      do idst = 1,icoeff
         kdst(:) = FLIndexQ(:,idst)
         ndst = sum(kdst(:))
         !
         ! if level is presented we only need the fields for iorder = level
         ! besides the summation must be taken only up to level-1
         !
         if (iopt==0.or.ndst==level) then
            !
            ! indexes for summation: i2 = 1..Ncoeff, k2 = (k1,k2,k3,k4...), n2 = k1+k2+k3..
            !
            i1 = 0 
            k1 = 0
            n1 = 0 
            !
            do while(n1<=ndst-iopt.and.i1<icoeff1)
               i1 = i1 + 1
               k1(:) = FLIndexQ(:,i1)
               n1 = sum(k1(:))

               k2(:) = kdst(:)-k1(:)
               if (all(k2(:)>=0)) then 
                 !
                 index1=1
                 index2=icoeff2
                 sig=2
                 !
                 do while(sig/=0)
                  !
                  iterm=min((index1+index2+sig)/2,icoeff2)
                  kindex(:)=FLIndexQ(:,iterm)
                  sig=0
                  !
                  if (sum(kindex)>sum(k2)) sig =-1
                  if (sum(kindex)<sum(k2)) sig = 1
                  !
                  if (sig==0) then 
                    !
                    do jmode=1,trove%Nmodes_e
                     !
                     if(kindex(jmode)<k2(jmode)) then 
                       sig= 1 
                     elseif(kindex(jmode)>k2(jmode)) then 
                       sig=-1
                     endif 
                     !
                     if (sig/=0) exit
                     !
                    enddo
                    !
                  endif 
                  !
                  if (sig== 1) index1=iterm
                  if (sig==-1) index2=iterm
                  !
                 enddo
                 !
                 i2 = iterm
                 !
                 if (i2<=icoeff2) then
                   !
                   dst(idst,:)=dst(idst,:) + src1(i1,:)*src2(i2,:)
                   !
                 endif 
               endif 
            enddo
            !
         endif
      enddo

   end subroutine product_of_polynoms_simple


  !
  ! Partial derivative of the polinomial form, given by indexes Kindex
  ! A polynom is defifed by the its order (Norder) and a number of modes (Nmodes), and, of course, 
  ! by its coeffitients fields(:), given as a 1D matrix (i=1..Ncoeff)
  ! The correlation between indexes k1,k2,k3,...k_Nmodes -> i are determined from recursive function FLQindex
  ! The inverse transformation must be prestored in  FLIndexQ(l,i), where l=1,2,3,...,NModes
  ! 
   subroutine deriv_of_polynoms_simple(osrc,isrc,src,kindex,idst,dst)

    integer(ik),intent(in)      :: osrc,isrc,idst
    real(ark),intent(in)  :: src(:,:)
    integer(ik)     ,intent(in)  :: kindex(:)
    real(ark),intent(out) :: dst(:,:)

    integer(ik) :: Norder

    integer(ik) :: kdst(trove%Nmodes),ksrc(trove%Nmodes),imode,i1,jdst,jsrc,ndst,DerivOrder
    real(ark)   :: factor
    integer(ik) :: index1,index2,sig,iterm,jmode
    integer(ik) :: k_t(trove%Nmodes)
    !

      !if (verbose>=7) write (out,"('Derivative of the polynom ',a,' goes to ',a)") trim(src%name),trim(dst%name)

      !
      ! Define the order of the dst-polynom as the order of the smallest polynom 
      !
      !
      !if (.not.allocated(FLIndexQ)) then 
      !  write (out,"('gmat_polynom: FLIndexQ has not been predefined ')")
      !  stop  'deriv_of_polynoms: FLIndexQ must have been already defined!'
      !endif
      !
      dst = 0
      !
      if (size(src,dim=1)<isrc.or.size(dst,dim=1)<idst) then
        write (out,"('deriv_of_polynoms: illegal sizes of dst or src : ',4i8)") size(src),isrc,size(dst),idst
        stop  'Error in deriv_of_polynoms: order of dst-polynom is too small'
      endif
      !
      ! Depending on the derivative order the polznom order is reduced 
      !
      DerivOrder = sum(kindex(:))
      !
      Norder = max(0,osrc - DerivOrder)
      !
      if (idst<trove%RangeOrder(Norder)) then
        write (out,"('Error in deriv_of_polynoms: size of dst-polynom is too small (isrc,idst,Norder): ',i8)") isrc,idst,Norder
        stop  'Error in deriv_of_polynoms: order of dst-polynom is too small'
      endif
      !
      do jdst = 1,idst
         !
         kdst(:) = FLIndexQ(:,jdst)
         ndst = sum(kdst(:))
         !
         ksrc(:) = kdst(:)+kindex(:)
         !
         !jsrc = FLQindex(trove%Nmodes_e,ksrc)
         !
         !jsrc = FLQindex(trove%Nmodes_e,ksrc)
         !
         index1=1
         index2=isrc
         sig=2
         !
         do while(sig/=0)
          !
          iterm=(index1+index2+sig)/2
          k_t(:)=FLIndexQ(:,iterm)
          sig=0
          !
          if (sum(k_t)>sum(ksrc)) sig =-1
          if (sum(k_t)<sum(ksrc)) sig = 1
          !
          if (sig==0) then 
            !
            do jmode=1,trove%Nmodes_e
             !
             if(k_t(jmode)<ksrc(jmode)) then 
               sig= 1 
             elseif(k_t(jmode)>ksrc(jmode)) then 
               sig=-1
             endif 
             !
             if (sig/=0) exit
             !
            enddo
            !
          endif 
          !
          if (sig== 1) index1=iterm
          if (sig==-1) index2=iterm
          !
         enddo
         !
         jsrc = iterm
         !
         factor = 1.0_ark
         do imode = 1,trove%Nmodes
           i1 = 0
           do while ( i1<kindex(imode) )
             i1 = i1+1
             factor = factor*(kdst(imode)+i1)
           enddo
         enddo
         !
         if (jsrc<=isrc) then 
            dst(jdst,:) = src(jsrc,:)*factor
         endif
      enddo

   end subroutine deriv_of_polynoms_simple


  !
  ! Partial derivative of the polinomial form, given by indexes Kindex
  ! A polynom is defifed by the its order (Norder) and a number of modes (Nmodes), and, of course, 
  ! by its coeffitients fields(:), given as a 1D matrix (i=1..Ncoeff)
  ! The correlation between indexes k1,k2,k3,...k_Nmodes -> i are determined from recursive function FLQindex
  ! The inverse transformation must be prestored in  FLIndexQ(l,i), where l=1,2,3,...,NModes
  ! 
   subroutine deriv_of_polynoms(src,kindex,dst)

    type(FLpolynomT),intent(in)  :: src
    integer(ik)     ,intent(in)  :: kindex(:)
    type(FLpolynomT),intent(out) :: dst


    integer(ik) :: kdst(trove%Nmodes),ksrc(trove%Nmodes),imode,i1,jdst,jsrc,ndst,DerivOrder,Norder
    real(rk)    :: factor
    !

      if (verbose>=7) write (out,"('Derivative of the polynom ',a,' goes to ',a)") trim(src%name),trim(dst%name)

      !
      ! Define the order of the dst-polynom as the order of the smallest polynom 
      !
      !
      if (.not.allocated(FLIndexQ)) then 
        write (out,"('gmat_polynom: FLIndexQ has not been predefined ')")
        stop  'deriv_of_polynoms: FLIndexQ must have been already defined!'
      endif
      !
      dst%field = 0
      !
      ! Depending on the derivative order the polznom order is reduced 
      !
      DerivOrder = sum(kindex(:))
      !
      !
      Norder = max(0,src%Orders - DerivOrder)
      !
      if (size(dst%field)<trove%RangeOrder(Norder)) then 
        write (out,"('Error in deriv_of_polynoms: size of dst-polynom is too small (src,dst): ',i8)") size(dst%field)
        stop  'Error in deriv_of_polynoms: order of dst-polynom is too small'
      endif
      !
      dst%Orders = Norder
      dst%Ncoeff = trove%RangeOrder(Norder)

      ! indexes of dst polynom
      do jdst = 1,dst%Ncoeff
         !
         kdst(:) = FLIndexQ(:,jdst)
         ndst = sum(kdst(:))
         !
         ksrc(:) = kdst(:)+kindex(:)
         !
         jsrc = FLQindex(trove%Nmodes_e,ksrc)
         !
         !
         factor = 1.0_ark
         do imode = 1,trove%Nmodes
           i1 = 0
           do while ( i1<kindex(imode) )
             i1 = i1+1
             factor = factor*(kdst(imode)+i1)
           enddo
         enddo
         !
         if (jsrc<=src%Ncoeff) then 
            dst%field(jdst,:) = src%field(jsrc,:)*factor
         endif
      enddo

   end subroutine deriv_of_polynoms



   ! 
   ! Relations between 1D array index j and Modes-Dimension array indexes i1,i2,i3... for
   ! coeffs(i1,i2,i3..) x1^i1 x2^i2 x3^i3 ...  -> coeffs(j) x1^i1 x2^i2 x3^i3 ...
   ! Input: Itarget is a vector of i1,i2,i3,.... indexes 
   ! Output: FLQindex = j is the number that corrsponds to vector of i1,i2,i3,...
   ! FLIndexQ(k,j) [optional] gives the backward correspondence j -> i_k for k-th mode
   !
   recursive function FLQindex(Nmodes,Itarget,FLIndexQ_) result (isum)
     ! integer(ik),intent(in) :: k(trove%Nmodes)
     !
     integer(ik),intent(in) :: Nmodes
     integer(ik),intent(in) :: Itarget(:)
     integer(ik),intent(out), optional  :: FLIndexQ_(:,:)
     integer(ik)            :: Isearch(size(Itarget))
     integer(ik) :: n0,n,imodes,isum,dm1,dm2,nsize
     logical     :: flag_go
      !
      imodes = 1
      !
      ! For the checks only 
      !
      dm1 = 0 
      dm2 = 0 
      !
      nsize = size(Itarget,dim=1)
      !
      if (verbose>=4.and.present(FLIndexQ_)) then 
        dm1 = size(FLIndexQ_,dim=1)
        dm2 = size(FLIndexQ_,dim=2)
        if (trove%Nmodes/=dm1) then 
            write(out,"('FLQindex: number of modes in FLIndexQ /= Nmodes:',2I8)") dm1,trove%Nmodes
            stop 'number of modes in FLIndexQ /= Nmodes'
        endif 
      endif 
      !
      isum = 0 
      N = sum(Itarget(1:Nmodes))
      !
      Isearch = 0 
      flag_go = .true.
      !
      n0 = 0 
      if (Nmodes==1) then 
        isum = N+1
        if (present(FLIndexQ_)) then 
          do n0=0,N
             FLIndexQ_(1,n0+1) = n0
          enddo
        endif 
      elseif(Nmodes==0) then 
        isum = 1
        if (present(FLIndexQ_))  FLIndexQ_ = 0
      else
        do while (n0<=N.and.flag_go)
           !
           if (present(FLIndexQ_)) then 
             !
             call fsum(Nmodes,nsize,n0,n0,imodes,flag_go,Itarget,Isearch,isum,FLIndexQ_)
             !
           else
             !
             call fsum(Nmodes,nsize,n0,n0,imodes,flag_go,Itarget,Isearch,isum)
             !
           endif
           !
           n0 = n0 +1
           !
        enddo
      endif
      !
      contains
      !
      recursive subroutine fsum(Nmodes,nsize,N,Norder,imodes,flag_go,Itarget,Isearch,isum,FLIndexQ_)
        integer(ik),intent(in) :: Nmodes,nsize,N,Norder,imodes
        logical                :: flag_go
        integer(ik),intent(in) :: Itarget(nsize)
        integer(ik),intent(inout) :: isum
        integer(ik),intent(out), optional  :: FLIndexQ_(:,:)
        integer(ik)            :: Isearch(nsize),dm2
        integer(ik) :: i0 
        !
        i0 = 0 
        do while (i0<=N.and.flag_go)
           !
           Isearch(imodes) = i0
           if (imodes == Nmodes-1) then
              Isearch(imodes+1) = Norder-sum(Isearch(1:Nmodes-1))
              !
              if (sum(nint(Isearch(1:Nmodes)*trove%PotPolyad(1:Nmodes)))<=Norder) then
                 !
                 isum = isum +1
                 !
                 if (present(FLIndexQ_)) then 
                    !
                    if (verbose>=4) then 
                       dm2 = size(FLIndexQ_,dim=2)
                       if (isum>dm2) then 
                         write(out,"('FLQindex: isum > size of FLIndexQ:',2i8)") isum,dm2
                         stop 'FLQindex: isum > size of FLIndexQ'
                       endif 
                    endif 
                    !
                    FLIndexQ_(:,isum) = Isearch(:)
                    !
                 endif 
                 !
              endif
              !
              if (all(Itarget(:)==Isearch(:))) flag_go = .false.
              !
           else
              !
              if (present(FLIndexQ_)) then 
                !
                call fsum(Nmodes,nsize,N-i0,Norder,imodes+1_ik,flag_go,Itarget,Isearch,isum,FLIndexQ_)
                !
              else
                !
                call fsum(Nmodes,nsize,N-i0,Norder,imodes+1_ik,flag_go,Itarget,Isearch,isum)
                !
              endif

              !
           endif 
           !
           i0 = i0 +1 
           !
        enddo
        !
        Isearch(imodes:Nmodes)=0
        !
      end subroutine fsum
      !
   end function FLQindex

!
! Initilizing the basis set 
!
  subroutine FLbsetInit
    !
    integer(ik)               :: numax(0:trove%Nmodes) 
    integer(ik)               :: ibset_total  
    integer(ik)               :: i, alloc, nmodes,imode,itype,maxnumax
    integer(ik)               :: isubsp1D
    type(Basis1DT),pointer    ::  bs1

    if (job%verbose>=5) write(out,"(/'FLbsetInit/start')")
    !
    if (job%verbose>=3) write(out,"(/'Initilizing the basis set ...')")
    ! 
    !
    nmodes = trove%Nmodes
    !
    numax(:) = job%bset(:)%range(2)
    ! 
    !
    ! define number of different basis set types 
    !
    isubsp1D = 0 
    do imode = 1,nmodes
       if ( job%bset(imode)%dim=='1D'.and.isubsp1D/=job%bset(imode)%species ) then 
            isubsp1D = isubsp1D +1 

            if ( isubsp1D/=job%bset(imode)%species ) then 
                 write(out,"('FLbsetInit: the species number of the basis set ',i6,', mode ',i4,' not valid')") & 
                              job%bset(imode)%species,imode
                 stop 'FLbsetInit: bad basis set species number'
            endif
            !
       endif
       if ( job%bset(imode)%dim/='1D' ) then 
          write(out,"('FLbsetInit: the dimensionality of the basis set ',a,' not valid')") job%bset(imode)%dim
          stop 'FLbsetInit: bad basis set dimensionality'
       endif
    enddo 
    !
    do imode = 1,trove%Nmodes
       if (job%bset(imode)%type/='NUMEROV') then 
         if ( job%bset(imode)%coord_kinet/=job%bset(imode)%coord_poten  ) then 
            write(out,"('FLbsetInit: Wrong defenition of coordinates:')")
            write(out,"('the kinetic and potential parts must be the same for the non/numerov-basis sets, while')")
            write(out,"('The kinetic   coordinate is :',a)") trim(job%bset(imode)%coord_kinet)
            write(out,"('The potential coordinate is :',a)") trim(job%bset(imode)%coord_poten)            
            stop 'FLbsetInit: Wrong defenition of coordinates'
         endif 
       endif
    enddo
    !
    bset%n_bset1D_max = isubsp1D
    !
    ibset_total = isubsp1D
    bset%n_bset_max = ibset_total

    if (ibset_total>trove%Nmodes) then 
      write(out,"('FLbsetInit: the number of basis set types is more that Nmodes:',2i10)") ibset_total,trove%Nmodes
      stop 'FLbsetInit: bad number of basis set types'
    endif 

    !
    ! Check if the bset spicies number gradualy increasing 
    !
    if (ibset_total>trove%Nmodes) then 
      write(out,"('FLbsetInit: the number of basis set types is more that Nmodes:',2i10)") ibset_total,trove%Nmodes
      stop 'FLbsetInit: bad number of basis set types'
    endif 
    !
    if (job%verbose>=5) write(out,"('Allocation of bset')") 
    !
    allocate (bset%bs1D(isubsp1D),bset%rot,stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i9,' initializing the basi set, ibset_1D = ',i9)") &
              alloc, isubsp1D
       stop 'BsetInit - alloc'
    end if
    !
    if (job%verbose>=5) write(out,"('Allocation of bset%dscr')") 
    !
    allocate (bset%dscr(0:trove%Nmodes),stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i9,' initializing the main basi set type, dscr',i8)") alloc
       stop 'BsetInit - alloc'
    end if
    !
    ! Store the initial basis set information with the bset object as descr
    !
    bset%dscr(0:Nmodes) = job%bset(0:Nmodes)
    !
    ! Inititialization of the 1d type 
    !
    if (job%verbose>=5) write(out,"('Inititialization of the 1d type')") 
    !
    do itype = 1,isubsp1D
       allocate(bset%bs1D(itype)%mode(trove%Nmodes),stat=alloc)
       if (alloc/=0) then
          write (out,"(' Error ',i9,' initializing bs%mode for 1d type ',i8)")  alloc, itype
          stop 'BsetInit - alloc'
       end if
       bset%bs1D(itype)%mode = 0
    enddo
    !
    ! define the basis set specifics: dimensions, which modes belong to every BS type 
    !
    bset%bs1D(:)%imodes = 0
    !
    !if (job%verbose>=6) write(out,"('give a name to the basis set 1d')")
    !
    do imode = 1,nmodes
       !
       do itype = 1,isubsp1D
           bs1 => bset%bs1D(itype)
           if ( itype==job%bset(imode)%species.and.job%bset(imode)%dim=='1D') then 
              bs1%imodes = bs1%imodes +1
              bs1%mode(bs1%imodes) = imode
           endif 
       enddo
       !
    enddo
    !
    ! give a name to the basis set (take the type of the 1st mode in bset%bs1D)
    !
    !if (job%verbose>=6) write(out,"('give a name to the basis set 1d')")
    !
    bset%rot%type = job%bset(0)%type
    !
    do itype = 1,isubsp1D
       bs1 => bset%bs1D(itype)
       imode = bs1%mode(1)
       write(bs1%name,"(a2,a,' #',i4)") trim(job%bset(imode)%dim),trim(job%bset(imode)%type),job%bset(imode)%species
       bs1%type = job%bset(imode)%type
    enddo 
    !
    ! Verbose: 
    !
    if (job%verbose>=4) then 
       do itype = 1,isubsp1D
          !bs1 => bset%bs1D(itype)
          if (bset%bs1D(itype)%imodes<=30) then 
             write(out,"('1d basis set, species ',i4,' for the modes ',30i5)") &
                  itype,(bset%bs1D(itype)%mode(imode),imode=1,bset%bs1D(itype)%imodes)
          endif
       enddo
       !
    endif

    !
    ! Inititialization of the rotation type 
    !
    if (job%verbose>=6) write(out,"('Inititialization of the rotational type')")
    !
    allocate(bset%rot%mode(1:1),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' initializing bs%mode for rot type ')")  alloc
        stop 'BsetInit-rot-mode-alloc'
    end if
    !
    bset%rot%type = job%bset(0)%type
    bset%rot%name = 'ROTATIONAL'
    bset%rot%imodes = 1  ! Obsolete
    bset%rot%mode   = 0  ! Obsolete
    !
    ! Complete the initilization of the 1d basis sets 
    ! generating the primitive matrix elements 
    !
    ! If the basis set has been stored we can just read it from the hard disk and leave
    !
    if (trim(trove%IO_hamiltonian)=='READ'.or.&
        trim(trove%IO_basisset)=='READ') then 
        !
        call FLcheck_point_Hamiltonian('BASIS_READ') 
        !
    else
      !
      do itype = 1,isubsp1D
         maxnumax = 0
         do i = 1,bset%bs1D(itype)%imodes
            imode = bset%bs1D(itype)%mode(i)
            !call PTgetsizeandorder(imode,nu_t)
            maxnumax = max(numax(imode),maxnumax)
         enddo 
         !
         ! Here we generate the basis sets and all matrix elements 
         !
         call FLbset1DNew(itype,maxnumax)
         !
      enddo
      !
      if (FLrotation) call FLRotbset_new(job%bset(0)%range(1))
      !
    endif
    !
    bset_initialized = .true.
    !
    if (trim(trove%IO_hamiltonian)=='SAVE'.or.trim(trove%IO_basisset)=='SAVE') then
      !
      call FLcheck_point_Hamiltonian('BASIS_SAVE') 
      !
    endif
    !
    if (job%verbose>=4) call MemoryReport
    !
    if (job%verbose>=4) write(out,"('...done!')") 
    !
   end subroutine FLbsetInit

   !
   ! read-write the primitive basis set information: 
   ! all primitive matrix elements;
   ! all objects defining the Hamiltonan expansion:
   ! Amat, Amat_rho, g_mat, poten, pseudo, and so on. 
   !
   subroutine FLcheck_point_Hamiltonian(action)

    character(len=*), intent(in) :: action ! 'SAVE' or 'READ'
    integer(ik)                  :: ntypes_stored,itype_stored(0:trove%Nmodes)
    !
    call TimerStart('FLcheck_point_Hamiltonian')
    select case (action)
      case default
        write (out,"(' FLcheck_point_Hamiltonian - action ',a,' is not valid')") trim(action)
        stop 'FLcheck_point_Hamiltonian - bogus command'
      case ('BASIS_SAVE')
        call basissetSave
        call numerovSave
      case ('AMAT_SAVE')
        call AmatBmatSave
      case ('HAMILTONIAN_SAVE')
        !
        call HamiltonianSave
        !
        !call HamiltonianSave_2007
        !
        if (trove%separate_store.and..not.trove%separate_convert) then
          !
          ! use kinetic.chk and potential.chk as ASCII
          !
          !call KineticSave_ASCII
          !call PotentialSave_ASCII
          !call ExternalSave_ASCII
          !
        endif
        !
      case ('BASIS_READ')
        call basisRestore
        call numerovRestore
      case ('KINETIC_READ')
        !
        call checkpointRestore_kinetic
        !
        if (trove%separate_convert.and..not.trove%separate_store) call KineticSave_ASCII
        !
        if (trove%separate_store) call checkpointRestore_kinetic_ascii
        !
      case ('KINETIC_SAVE_SPARSE')
        call AmatBmatSave
        call KineticSave_ASCII
      case ('POTENTIAL_SAVE_SPARSE')
        call PotentialSave_ASCII
      case ('EXTERNAL_SAVE_SPARSE')
        call ExternalSave_ASCII
      case ('KINETIC_SKIP')
        call checkpointSkip_kinetic
      case ('POTENTIAL_READ')
        !
        if (trove%separate_store) then 
          call checkpointRestore_potential_ascii
        else
          call checkpointRestore_potential
        endif
        !
        if (trove%separate_convert.and..not.trove%separate_store) call PotentialSave_ASCII
        !
      case ('EXTERNAL_READ')
        !
        if (trove%separate_store) then
          call checkpointRestore_external_ascii
        else
          call checkpointRestore_external
        endif
        !
        if (trove%separate_convert.and..not.trove%separate_store) call ExternalSave_ASCII
        !
    end select
    call TimerStop('FLcheck_point_Hamiltonian')

    contains

      subroutine basissetSave

        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,itype,imode,i
        type(Basis1DT),pointer    ::  bs
        type(FLpolynomT),pointer    :: fl
        integer(ik)        :: Nmodes,k1,k2,imu
        !
        if (job%verbose>=3) write(out,"(/'Store primitive matrix elements...')")
        !
        !unitfname ='Check point of the Hamiltonian'
        unitfname ='Check point of the basis set'
        call IOStart(trim(unitfname),chkptIO)
        !
        open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=trove%chk_fname)
        !
        if (trove%sparse) then 
           write(chkptIO) 'Start SparseBas'
        else
           write(chkptIO) 'Start Basis set'
        endif
        !
        write(chkptIO) bset%n_bset1D_max
        !
        loop_bset : do itype = 0,bset%n_bset1D_max
           !
           if (itype==0) then 
             bs => bset%rot
             if (.not.FLrotation) cycle 
           else
             bs => bset%bs1D(itype)
           endif
           !
           do i = 1,bs%imodes
             !
             imode = bs%mode(i)
             !
             write(chkptIO) itype
             write(chkptIO) bs%name
             write(chkptIO) bs%type
             write(chkptIO) bs%imodes
             write(chkptIO) bs%mode(i)
             write(chkptIO) bs%size
             write(chkptIO) bs%order
             write(chkptIO) bs%ener0
             write(chkptIO) bs%params
             write(chkptIO) bs%matelements
             !
           enddo 
           !
        enddo loop_bset
        !
        ! Techically the data below belong to the Hamiltonian fields poten, pseudo, g_vib, ..
        ! However this part of these fields contains the matrix elements related 
        ! to the last mode 'Nmodes'.  It is needed for the final Hamiltonian 
        ! constraction, which is why we store this information as well. 
        !
        itype = bset%n_bset1D_max
        bs => bset%bs1D(itype)
        Nmodes = trove%Nmodes
        !
        if (all(bs%mode(:)/=Nmodes)) then 
          write (out,"(' check_point_Hamiltonian: the last species does not contain the last mode',30i7)") bs%mode(:)
          stop 'check_point_Hamiltonian - bset type mismatch'
        endif 
        !
        if (trove%sparse) then 
          write(chkptIO) 'sparse'
          write(chkptIO) job%exp_coeff_thresh
        endif
        !
        write(chkptIO) 'g_vib'
        if (.not.trove%sparse) write(chkptIO) trove%g_vib(1,1)%Ncoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,Nmodes
              !
              fl => trove%g_vib(k1,k2)
              !
              if (trove%sparse) write(chkptIO) fl%Ncoeff
              !
              if (fl%Ncoeff>0) write(chkptIO) fl%me 
              !
           enddo
        enddo
        !
        write(chkptIO) 'g_rot'
        if (.not.trove%sparse) write(chkptIO) trove%g_rot(1,1)%Ncoeff
        !
        do k1 = 1,3
           do k2 = 1,3
              !
              fl => trove%g_rot(k1,k2)
              if (trove%sparse) write(chkptIO) fl%Ncoeff
              if (fl%Ncoeff>0) write(chkptIO) fl%me 
              !
           enddo
        enddo
        !
        write(chkptIO) 'g_cor'
        if (.not.trove%sparse) write(chkptIO) trove%g_cor(1,1)%Ncoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,3
              !
              fl => trove%g_cor(k1,k2)
              if (trove%sparse) write(chkptIO) fl%Ncoeff
              if (fl%Ncoeff>0) write(chkptIO) fl%me  
              !
           enddo
        enddo
        !
        write(chkptIO) 'poten'
        write(chkptIO) trove%poten%Ncoeff
        !
        fl => trove%poten
        !
        write(chkptIO) trove%poten%me 
        !
        if (FLL2_coeffs) then
          !
          write(chkptIO) 'L2_vib'
          if (.not.trove%sparse) write(chkptIO) trove%L2_vib(1,1)%Ncoeff
          !
          do k1 = 1,Nmodes
             do k2 = 1,Nmodes
                !
                fl => trove%L2_vib(k1,k2)
                if (trove%sparse) write(chkptIO) fl%Ncoeff
                !
                if (fl%Ncoeff>0) write(chkptIO) fl%me 
                !
             enddo
          enddo
          !
        endif 
        !
        if (FLextF_coeffs) then
          !
          write(chkptIO) 'ext_f'
          !
          do imu = 1,extF%rank
             !
             fl => trove%extF(imu)
             !
             write(chkptIO) fl%Ncoeff
             if (fl%Ncoeff>0) write(chkptIO) fl%me 
             !
          enddo
          !
        endif
        !
        write(chkptIO) 'End Basis set'
        close(chkptIO,status='keep')
        !
      end subroutine basissetSave
      !
      subroutine basisRestore

        character(len=15) :: buf
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,itype,alloc,imode,jmax
        type(Basis1DT),pointer    ::  bs
        type(FLpolynomT),pointer    :: fl
        integer(ik)          :: Nmodes,Npoints,k1,k2,Tcoeff,isize,jmax_t,i,imu
        integer(hik)         :: totalsize 
        real(ark),allocatable :: matelements_t(:,:,:,:),ener0_t(:)
        real(ark)             :: b_param_t(1:3)
        logical :: create_new_rot_basis = .false.
        real(rk) :: exp_coeff_thresh

        unitfname ='Check point of the basis set'
        !unitfname ='Check point of the Hamiltonian'
        call IOStart(trim(unitfname),chkptIO)
        open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=trove%chk_fname)
        !
        read(chkptIO) buf
        if (trove%sparse.and.buf=='Start Basis set') then
          write (out,"(' Checkpoint non-sparse file',a,' is used  for sparse job: ',a)") trove%chk_fname, buf
          stop 'check_point_Hamiltonian - non-sparse used for sparse job (1)'
        end if
        !
        if (.not.trove%sparse.and.buf=='Start SparseBas') then
          write (out,"(' Checkpoint sparse file ',a,' ised for non-sparse job: ',a)") trove%chk_fname, buf
          stop 'check_point_Hamiltonian - sparse used for non-sparse job (2)'
        end if

        if (trove%sparse.and.buf/='Start SparseBas'.or.(.not.trove%sparse.and.buf/='Start Basis set')) then
          write (out,"(' Checkpoint file ',a,' has bogus header: ',a)") trove%chk_fname, buf
          stop 'check_point_Hamiltonian - bogus file format (3)'
        end if

        !
        read(chkptIO) ntypes_stored
        !
        do itype = 0,bset%n_bset1D_max
           !
           if (itype==0) then 
             bs => bset%rot
             if (.not.FLrotation) cycle 
             imode = 0
           else
             bs => bset%bs1D(itype)
             imode = bset%bs1D(itype)%mode(1)
           endif 
           !
           do i = 1,bs%imodes
             !
             imode = bs%mode(i)
             !
             read(chkptIO) itype_stored(imode)
             !
             read(chkptIO) bs%name
             !
             read(chkptIO) bs%type
             !
             if (trim(bs%type)/=trim(job%bset(imode)%type)) then
               write (out,"(' check_point_Hamiltonian:  bset type mismatch: ',a,2x,a)") trim(bs%type),trim(job%bset(imode)%type)
               stop 'check_point_Hamiltonian - bset type mismatch'
             end if
             !
             read(chkptIO) bs%imodes
             !
             if (bs%imodes<1.or.bs%imodes>trove%Nmodes) then
               write (out,"(' check_point_Hamiltonian: parameter mismatch, imodes : ',i8)") bs%imodes 
               stop 'check_point_Hamiltonian - parameter mismatch, imodes'
             end if
             !
             read(chkptIO) bs%mode(i)
             read(chkptIO) bs%size
             read(chkptIO) bs%order
             !
             if (any(bs%mode(1:bs%imodes)<0).or.any(bs%mode(1:bs%imodes)>trove%Nmodes)) then
               write (out,"(' check_point_Hamiltonian:  parameter mismatch, modes ',30i7)") bs%mode(1:bs%imodes)
               stop 'check_point_Hamiltonian - parameter mismatch, imodes'
             end if
             !
             if (itype>0.and.bs%size<job%bset(imode)%range(2)) then
               write (out,"(' check_point_Hamiltonian:  parameter mismatch, size ',2i7)") bs%size,job%bset(imode)%range(2)
               stop 'check_point_Hamiltonian - parameter mismatch, size'
             end if
             !
             if (itype>0.and.bs%order/=trove%MaxOrder) then
               write (out,"(' check_point_Hamiltonian:  parameter mismatch, order ',2i7)") bs%order,trove%MaxOrder
               write (out,"(' Some of orders (KinOrder, PotOrder, or external) is inconsistent with the checkpointed value')")
               stop 'check_point_Hamiltonian - parameter mismatch, order'
             end if
             !
             if (itype==0) then 
               !
               jmax = job%bset(0)%range(1)
               jmax_t = bs%order
               !
               if (bs%size/=(jmax+1)*(jmax+2)/2-1.or.jmax_t/=jmax) then
                 !
                 create_new_rot_basis = .true.
                 !
                 !write (out,"(' check_point_Hamiltonian:  parameter mismatch, size ',2i7)") bs%size,job%bset(imode)%range(2)
                 !
                 ! Inititialization of the rotation type 
                 !
                 if (job%verbose>=5) write(out,"('Inititialization of the rotational type')")
                 !
                 allocate(bset%rot%mode(1:1),stat=alloc)
                 if (alloc/=0) then
                     write (out,"(' Error ',i9,' initializing bs%mode for rot type ')")  alloc
                     stop 'BsetInit-rot-mode-alloc'
                 end if
                 !
                 bset%rot%type = job%bset(0)%type
                 bset%rot%name = 'ROTATIONAL'
                 bset%rot%imodes = 1  ! Obsolete
                 bset%rot%mode   = 0  ! Obsolete
                 !
                 isize = (jmax_t+1)*(jmax_t+2)/2-1
                 !
                 if (FLrotation) call FLRotbset_new(job%bset(0)%range(1))
                 !
               else
                 !
                 isize = bs%size
                 allocate (bs%matelements(7,(jmax+1)*(jmax+2)/2,-2:2,0:1),bs%ener0(0:isize),stat=alloc)
                 call ArrayStart('bs%matelements-rot',alloc,size(bs%matelements),kind(bs%matelements))
                 call ArrayStart('bs%ener0-rot',alloc,size(bs%ener0),kind(bs%ener0))
                 !
               endif
               !
               allocate (matelements_t(7,isize+1,-2:2,0:1),ener0_t(0:isize),stat=alloc)
               !
             else
               !
               isize = job%bset(imode)%range(2)
               !
               allocate (matelements_t(-1:3,0:trove%MaxOrder,0:bs%size,0:bs%size),ener0_t(0:bs%size),stat=alloc)
               !
             endif 
             !
             read(chkptIO) ener0_t(0:)
             read(chkptIO) b_param_t(1:3)
             read(chkptIO) matelements_t
             !
             if (itype==0) then 
               !
               if (.not.create_new_rot_basis) then 
                 bs%matelements = matelements_t
                 bs%ener0(0:isize) = ener0_t(0:isize)
               endif 
               !
             elseif(i==1) then 
               !
               ! Only the first imode within the givem type is unique  
               !
               allocate (bs%matelements(-1:3,0:trove%MaxOrder,0:isize,0:isize),bs%ener0(0:isize),stat=alloc)
               call ArrayStart('bs%matelements',alloc,1_ik,kind(bs%matelements),size(bs%matelements,kind=hik))
               call ArrayStart('bs%ener0',alloc,size(bs%ener0),kind(bs%ener0))
               !
               bs%matelements(-1:3,0:trove%MaxOrder,0:isize,0:isize) = matelements_t(-1:3,0:trove%MaxOrder,0:isize,0:isize)
               bs%params(1:3) = b_param_t(1:3)
               bs%ener0(0:isize) = ener0_t(0:isize)
               !
             endif
             !
             deallocate(matelements_t,ener0_t)
             !
           enddo
           !
        enddo
        !
        ! Now we read the matrix elements related to the last mode
        !
        itype = bset%n_bset1D_max
        bs => bset%bs1D(itype)
        Nmodes = trove%Nmodes
        Npoints = trove%Npoints
        !
        if (all(bs%mode(:)/=Nmodes)) then 
          write (out,"(' check_point_Hamiltonian: the last species does not contain the last mode',30i7)") bs%mode(:)
          stop 'check_point_Hamiltonian - bset type mismatch'
        endif 
        !
        if (.not.associated(trove%g_vib).or..not.associated(trove%g_rot).or. &
            .not.associated(trove%g_cor)) then 
           !
           allocate (trove%poten,trove%g_vib(Nmodes,Nmodes),trove%g_rot(3,3),trove%g_cor(Nmodes,3),trove%pseudo,stat=alloc)
           if (alloc/=0) then
               write (out,"('chk_Restore_kin-Error ',i9,' trying to allocate g-fields')") alloc
               stop 'chk_Restore_kin, g-fields - out of memory'
           end if
           !
        endif
        !
        ! sparse check 
        !
        if (trove%sparse) then
          !
          read(chkptIO) buf(1:6)
          !
          if (buf(1:6)/='sparse') then
            write (out,"(' Checkpoint file ',a,' in g_vib for sparse has bogus label ',a,'sparse expected')") &
                   trove%chk_fname, buf(1:6)
            stop 'check_point_Hamiltonian - bogus file format g_vib sparse'
          end if
          !
          read(chkptIO) exp_coeff_thresh
          !
          if ( abs(exp_coeff_thresh-job%exp_coeff_thresh)>small_ ) then
            !
            write(out,"(a,e18.10,a,e18.10,a)") & 
                        'basisRestore-gvib: exp_coeff_thresh used ',job%exp_coeff_thresh>exp_coeff_thresh,&
                        ' is incompatible with stored ',', change threshold or redo BASIS_SET SAVE'
            stop 'basisRestore-gvib: exp_coeff_thresh is incompatible with BASIS, change threshold or BASIS_SET SAVE'
            !
          endif
          !
        endif
        !
        read(chkptIO) buf(1:5)
        !
        !if (.not.trove%sparse.and.trim(trove%IO_basisset)=="READ") then
        !  write(out,"('basisRestore: A sparse option was used to save kinetic.chk, please switch SPARSE on by adding SPARSE to the inpiut')")
        !  stop 'basisRestore error SPARSE should be switched on'
        !endif
        !
        if (buf(1:5)/='g_vib') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_vib ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format g_vib'
        end if
        !
        if (.not.trove%sparse) read(chkptIO) Tcoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,Nmodes
              !
              fl => trove%g_vib(k1,k2)
              !
              if (trove%sparse) then 
                read(chkptIO) Tcoeff
                !
                if (fl%Ncoeff/= Tcoeff) then 
                  write (out,"(' Checkpoint file ',a,':  Ncoeff (basis) in g_vib disagree with ncoeff of field',2i4,1x,2I8)") &
                              k1,k2,fl%Ncoeff,Tcoeff
                  write (out,"('Consider switching BASIS_SET SAVE')")
                  stop 'check_point_Hamiltonian - Ncoeff (basis) in g_vib disagree with ncoeff of field'
                end if 
                !
              endif
              !
              fl%Ncoeff = Tcoeff
              !
              if (fl%Ncoeff>0) then 
                 !
                 allocate (fl%me(fl%Ncoeff,0:bs%Size,0:bs%Size),stat=alloc)
                 call ArrayStart('trove%g_vib%me',alloc,1_ik,kind(fl%me),size(fl%me,kind=hik))
                 read(chkptIO) fl%me     !(fl%Ncoeff,0:bs%Size,0:bs%Size)
                 !
              endif
              !
              !read(chkptIO) fl%iorder !(fl%Ncoeff)
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='g_rot') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_rot ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format g_rot'
        end if
        !
        if (.not.trove%sparse) read(chkptIO) Tcoeff
        !
        do k1 = 1,3
           do k2 = 1,3
              !
              fl => trove%g_rot(k1,k2)
              !
              if (trove%sparse) then 
                read(chkptIO) Tcoeff
                !
                if (fl%Ncoeff/= Tcoeff) then 
                  write (out,"(' Checkpoint file ',a,':  Ncoeff (basis) in g_rot disagree with ncoeff of field',2i4,1x,2I8)") &
                                trim(trove%chk_fname),k1,k2,fl%Ncoeff,Tcoeff
                  write (out,"(' Check exp_coeff_thresh and consider switching BASIS_SET SAVE')")
                  stop 'check_point_Hamiltonian - Ncoeff (basis) in g_rot disagree with ncoeff of field'
                end if 
              endif
              !
              fl%Ncoeff = Tcoeff
              !
              if (fl%Ncoeff>0) then 
                 allocate (fl%me(fl%Ncoeff,0:bs%Size,0:bs%Size),stat=alloc)
                 call ArrayStart('trove%g_rot%me',alloc,1_ik,kind(fl%me),size(fl%me,kind=hik))
                 read(chkptIO) fl%me     !(fl%Ncoeff,0:bs%Size,0:bs%Size)
              endif
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='g_cor') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_cor ',a)") trove%chk_fname, buf(1:6)
          stop 'check_point_Hamiltonian - bogus file format g_cor'
        end if
        !
        if (.not.trove%sparse) read(chkptIO) Tcoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,3
              !
              fl => trove%g_cor(k1,k2)
              !
              if (trove%sparse) then 
                read(chkptIO) Tcoeff
                !
                if (fl%Ncoeff/= Tcoeff) then 
                  write (out,"(' Checkpoint file ',a,':  Ncoeff (basis) in g_cor disagree with ncoeff of field',2i4,1x,2I8)") &
                          trove%chk_fname,k1,k2,fl%Ncoeff,Tcoeff
                  write (out,"(' Check exp_coeff_thresh and consider switching BASIS_SET SAVE')")
                  stop 'check_point_Hamiltonian - Ncoeff (basis) in g_cor disagree with ncoeff of field'
                end if 
              endif
              !
              fl%Ncoeff = Tcoeff
              !
              if (fl%Ncoeff>0) then 
                allocate (fl%me(fl%Ncoeff,0:bs%Size,0:bs%Size),stat=alloc)
                call ArrayStart('trove%g_cor%me',alloc,1_ik,kind(fl%me),size(fl%me,kind=hik))
                read(chkptIO) fl%me
              endif
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='poten') then
          write (out,"(' Checkpoint file ',a,' has bogus label poten ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format poten'
        end if
        !
        read(chkptIO) Tcoeff
        !
        if (trove%poten%Ncoeff/=Tcoeff) then 
          write (out,"(' Checkpoint file ',a,':  Ncoeff (basis) in poten  disagrees with ncoeff of field',2i4,1x,2I8)") &
                        trim(trove%chk_fname),k1,k2,fl%Ncoeff,Tcoeff
          write (out,"(' Check exp_coeff_thresh and consider switching BASIS_SET SAVE')")
          stop 'check_point_Hamiltonian - Ncoeff (basis) in poten disagree with ncoeff of field'
        end if 
        !
        if (Tcoeff==0) then 
          write (out,"(' Checkpoint Ncoeff (basis) in poten==0')") 
          stop ' Checkpoint Ncoeff (basis) in poten==0'
        end if 
        !
        trove%poten%Ncoeff = Tcoeff
        !
        allocate (trove%poten%me(trove%poten%Ncoeff,0:bs%Size,0:bs%Size),stat=alloc)
        call ArrayStart('trove%poten%me',alloc,1_ik,kind(fl%me),size(fl%me,kind=hik))
        !
        fl => trove%poten
        read(chkptIO) fl%me  
        !
        !
        ! L2 field function 
        !
        if (FLL2_coeffs) then
          !
          if (.not.associated(trove%L2_vib)) then 
            allocate (trove%L2_vib(Nmodes,Nmodes),stat=alloc)
            if (alloc/=0)  stop 'chk_Restore_L2, L2_vib-fields - out of memory'
          endif 
          !
          read(chkptIO) buf(1:6)
          if (buf(1:6)/='L2_vib') then
            write (out,"(' Checkpoint file ',a,' has bogus label L2_vib ',a)") trove%chk_fname, buf(1:6)
            stop 'check_point_Hamiltonian - bogus file format L2_vib'
          end if
          !
          if (.not.trove%sparse) read(chkptIO) Tcoeff
          !
          do k1 = 1,Nmodes
             do k2 = 1,Nmodes
                !
                fl => trove%L2_vib(k1,k2)
                if (trove%sparse) then 
                  read(chkptIO) Tcoeff
                  !
                  if (fl%Ncoeff/= Tcoeff) then 
                    write (out,"(' Checkpoint file ',a,':  Ncoeff (basis) in L2vib disagree with ncoeff of field',2i4,1x,2I8)") &
                                   trim(trove%chk_fname),k1,k2,fl%Ncoeff,Tcoeff
                    write (out,"(' Check exp_coeff_thresh and consider switching BASIS_SET SAVE')")
                    stop 'check_point_Hamiltonian - Ncoeff (basis) in L2vib disagree with ncoeff of field'
                  end if 
                endif
                !
                fl%Ncoeff = Tcoeff
                !
                if (fl%Ncoeff>0) then 
                  allocate (fl%me(fl%Ncoeff,0:bs%Size,0:bs%Size),stat=alloc)
                  call ArrayStart('trove%L2_vib%me',alloc,1_ik,kind(fl%me),size(fl%me,kind=hik))
                  read(chkptIO) fl%me
                endif
                !
             enddo
          enddo
          !
        endif 
        !
        ! External field function 
        !
        if (FLextF_coeffs) then
          !
          if (.not.associated(trove%extF)) then 
            allocate (trove%extF(extF%rank),stat=alloc)
            if (alloc/=0)  stop 'chk_Restore_extF, extF-fields - out of memory'
          endif 
          !
          read(chkptIO) buf(1:5)
          if (buf(1:5)/='ext_f') then
            write (out,"(' Checkpoint file ',a,' has bogus label ext_f ',a)") trove%chk_fname, buf(1:6)
            stop 'check_point_Hamiltonian - bogus file format ext_f'
          end if
          !
          do imu = 1,extF%rank
             !
             fl => trove%extF(imu)
             !
             !if (trove%sparse) then 
             read(chkptIO) Tcoeff
             !
             if (fl%Ncoeff/= Tcoeff) then 
               write (out,"(' Checkpoint file ',a,':  Ncoeff (basis) in extF disagree with ncoeff of field',i4,1x,2I8)") &
                              trim(trove%chk_fname),imu,fl%Ncoeff,Tcoeff
               write (out,"(' Check exp_coeff_thresh Consider switching BASIS_SET SAVE')")
               stop 'check_point_Hamiltonian - Ncoeff (basis) in extF disagree with ncoeff of field'
             end if 
             !
             fl%Ncoeff = Tcoeff
             !
             if (fl%Ncoeff>0) then 
                allocate (fl%me(fl%Ncoeff,0:bs%Size,0:bs%Size),stat=alloc)
                totalsize = int(fl%Ncoeff,hik)*int(bs%Size+1_hik,hik)*int(bs%Size+1_hik,hik)
                if (job%verbose>=6) write(out,"('Allocating trove%exfF%me:',i12,' x ',i9,' x ',i9,' = ',i14)") &
                                    fl%Ncoeff,bs%Size+1,bs%Size+1,totalsize
                !
                call ArrayStart('trove%exfF%me',alloc,1_ik,kind(fl%me),totalsize)
                read(chkptIO) fl%me
             endif
             !
          enddo
          !
        endif 
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='End B'.and.buf(1:5)/='ext_f'.and.buf(1:5)/='L2_vib') then
          write (out,"(' Checkpoint file ',a,' has bogus final header: ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format (2)'
        end if
        close(chkptIO,status='keep')
        !
        call MemoryReport
        !
      end subroutine basisRestore
      !
      ! we use this routine to pack all numerov basis functions into a single check point file together.
      ! we read all records from the scratch files and copy the into one place. 
      !
      subroutine numerovSave
        !
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,io_slot,itype,imode
        type(Basis1DT),pointer ::  bs
        integer(ik)        :: rec_len,vrec,v,maxpoints,npoints,alloc,bs_size,i
        real(ark) ,allocatable     :: bs_funct(:),dbs_funct(:)
        !
        maxpoints = 0
        !
        do itype = 1,bset%n_bset1D_max
           !
           bs => bset%bs1D(itype)
           !
           imode = bs%mode(1)
           !
           maxpoints = max(maxpoints,job%bset(imode)%npoints)
           !
        enddo 
        !
        !
        allocate (bs_funct(0:maxpoints),dbs_funct(0:maxpoints),stat=alloc)
        if (alloc/=0) then
           write (out,"('numerovSave: Error ',i9,' trying to allocate bs_funct-field')") alloc
           stop 'numerovSave, bs_funct -  out of memory'
        end if
        !
        inquire(iolength=rec_len) bs_funct(:),dbs_funct(:)
        !
        unitfname ='Check point of the numerov wave functions'
        call IOStart(trim(unitfname),chkptIO)
        open(chkptIO,action='write',status='replace',file=trove%chk_numerov_fname,access='direct',recl=rec_len)
        !
        do itype = 1,bset%n_bset1D_max
          !
          bs => bset%bs1D(itype)
          !
          imode = bs%mode(1)
          !
          if( bset%dscr(imode)%type/='NUMEROV'.and.bset%dscr(imode)%type/='FOURIER'.and.&
              bset%dscr(imode)%type/='LEGENDRE'.and.bset%dscr(imode)%type/='SINRHO'.and.&
              bset%dscr(imode)%type/='LAGUERRE-K') cycle
          !
          write(unitfname,"('Numerov basis set # ',i6)") imode
          ! get the i/o unit with stored numerov basis functions
          call IOStart(trim(unitfname),io_slot)
          !open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
          !
          npoints = job%bset(imode)%npoints
          !
          bs_size  = job%bset(imode)%range(2) 
          !
          do v = 0,bs_size+1
            !
            if(v>bs_size.and.bset%dscr(imode)%type/='NUMEROV') cycle
            !
            ! copy:
            read (io_slot,rec=v+1) (bs_funct(i),i=0,npoints),(dbs_funct(i),i=0,npoints)
            !
            ! and paste:
            vrec = itype+v*bset%n_bset1D_max
            !
            write (chkptIO,rec=vrec) (bs_funct(i),i=0,npoints),(dbs_funct(i),i=0,npoints)
            !
          enddo
           !
        enddo 
        !
        deallocate (bs_funct,dbs_funct)
        !
        close(chkptIO,status='keep')
        !
      end subroutine numerovSave


      !
      ! here we unpack the numerov wave function and copy into designated records
      !
      subroutine numerovRestore
        !
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,io_slot,itype,imode
        type(Basis1DT),pointer      ::  bs
        integer(ik)        :: rec_len,vrec,v,maxpoints,npoints,alloc,bs_size,i,itype_t
        real(ark) ,allocatable         :: bs_funct(:),dbs_funct(:)
        !
        maxpoints = 0
        !
        do itype = 1,bset%n_bset1D_max
           !
           bs => bset%bs1D(itype)
           imode = bs%mode(1)
           maxpoints = max(maxpoints,job%bset(imode)%npoints)
           !
        enddo 
        !
        allocate (bs_funct(0:maxpoints),dbs_funct(0:maxpoints),stat=alloc)
        if (alloc/=0) then
           write (out,"('numerovRestore: Error ',i9,' trying to allocate bs_funct-field')") alloc
           stop 'numerovRestore, bs_funct -  out of memory'
        end if
        !
        inquire(iolength=rec_len) bs_funct(:),dbs_funct(:)
        !
        unitfname ='Check point of the numerov wave functions'
        call IOStart(trim(unitfname),chkptIO)
        open(chkptIO,action='read',status='old',file=trove%chk_numerov_fname,&
                     access='direct',recl=rec_len)
        !
        do itype = 1,bset%n_bset1D_max
          !
          bs => bset%bs1D(itype)
          !
          imode = bs%mode(1)
          !
          itype_t = itype_stored(imode)
          !
          if( bset%dscr(imode)%type/='NUMEROV'.and.bset%dscr(imode)%type/='FOURIER'.and.&
              bset%dscr(imode)%type/='LEGENDRE'.and.bset%dscr(imode)%type/='SINRHO'.and.&
              bset%dscr(imode)%type/='LAGUERRE-K') cycle
          !
          npoints = job%bset(imode)%npoints
          !
          bs_size  = job%bset(imode)%range(2) 
          !
          inquire(iolength=rec_len) bs_funct(0:npoints),dbs_funct(0:npoints)
          !
          write(unitfname,"('Numerov basis set # ',i6)") imode
          ! get the i/o unit with stored numerov basis functions
          call IOStart(trim(unitfname),io_slot)
          open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
          !
          do v = 0,bs_size+1
            !
            if(v>bs_size.and.bset%dscr(imode)%type/='NUMEROV') cycle
            !
            ! copy:
            vrec = itype_t+v*ntypes_stored !bset%n_bset1D_max
            !
            read (chkptIO,rec=vrec) (bs_funct(i),i=0,npoints),(dbs_funct(i),i=0,npoints)
            !
            ! and paste:
            write (io_slot,rec=v+1) (bs_funct(i),i=0,npoints),(dbs_funct(i),i=0,npoints)
            !
          enddo
           !
        enddo 
        !
        deallocate (bs_funct,dbs_funct)
        !
        close(chkptIO,status='keep')
        !
      end subroutine numerovRestore
      !

      subroutine AmatBmatSave
        !
        character(len=cl)  :: unitfname
        integer(ik)        :: Nmodes,k1,k2,chkptIO,chkptIO_pot,chkptIO_kin,chkptIO_ext,i,iterm,npoints
        type(FLpolynomT),pointer    :: fl
        logical     :: i_opened
        !
        if (job%verbose>=3) write(out,"(/'Store Amat and Bmat objects (sparse) ...')")
        !
        if (trove%separate_store.and.trim(trove%internal_coords)=='LOCAL') return
        if (.not.trove%separate_store.or.trove%separate_convert) return
        !
        unitfname ='Check point of the Hamiltonian'
        call IOStart(trim(unitfname),chkptIO)
        !
        inquire (chkptIO,opened=i_opened)
        !
        if (i_opened) close(chkptIO)
        !
        open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=trove%chk_hamil_fname)
        !
        write(chkptIO) 'Start Hamiltonian objects'
        !
        Nmodes = trove%Nmodes
        Npoints = trove%npoints
        !
        write(chkptIO) 'Amatrho'
        write(chkptIO) trove%Amatrho
        !
        write(chkptIO) 'dAmatrho'
        write(chkptIO) trove%dAmatrho
        !
        write(chkptIO) 'Bmatrho'
        write(chkptIO) trove%Bmatrho
        !
        write(chkptIO) 'dBmatrho'
        write(chkptIO) trove%dBmatrho
        !
        close(chkptIO,status='keep')
        !
      end subroutine AmatBmatSave


      subroutine HamiltonianSave

        character(len=cl)  :: unitfname
        integer(ik)        :: Nmodes,k1,k2,chkptIO,chkptIO_pot,chkptIO_kin,chkptIO_ext,i,iterm,npoints
        type(FLpolynomT),pointer    :: fl
        logical     :: i_opened
        !
        if (job%verbose>=3) write(out,"(/'Store all objects defining the Hamiltonian...')")
        !
        if (trove%separate_store.and.trim(trove%internal_coords)=='LOCAL') return
        !
        if (trove%separate_store.and..not.trove%separate_convert) return
        !
        unitfname ='Check point of the Hamiltonian'
        call IOStart(trim(unitfname),chkptIO)
        !
        inquire (chkptIO,opened=i_opened)
        !
        if (i_opened) close(chkptIO)
        !
        open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=trove%chk_hamil_fname)
        !
        write(chkptIO) 'Start Hamiltonian objects'
        !
        Nmodes = trove%Nmodes
        Npoints = trove%npoints
        !
        write(chkptIO) 'Amatrho'
        write(chkptIO) trove%Amatrho
        !
        write(chkptIO) 'dAmatrho'
        write(chkptIO) trove%dAmatrho
        !
        write(chkptIO) 'Bmatrho'
        write(chkptIO) trove%Bmatrho
        !
        write(chkptIO) 'dBmatrho'
        write(chkptIO) trove%dBmatrho
        !
        write(chkptIO) 'g_vib'
        write(chkptIO) trove%g_vib(1,1)%Ncoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,Nmodes
              !
              write(chkptIO) trove%g_vib(k1,k2)%field
              write(chkptIO) trove%g_vib(k1,k2)%iorder
              !
           enddo
        enddo
        !
        write(chkptIO) 'g_rot'
        write(chkptIO) trove%g_rot(1,1)%Ncoeff
        !
        do k1 = 1,3
           do k2 = 1,3
              !
              write(chkptIO) trove%g_rot(k1,k2)%field
              write(chkptIO) trove%g_rot(k1,k2)%iorder
              !
           enddo
        enddo
        !
        write(chkptIO) 'g_cor'
        write(chkptIO) trove%g_cor(1,1)%Ncoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,3
              !
              write(chkptIO) trove%g_cor(k1,k2)%field
              write(chkptIO) trove%g_cor(k1,k2)%iorder
              !
           enddo
        enddo
        !
        fl => trove%pseudo
        !
        write(chkptIO) 'pseudo'
        write(chkptIO) trove%pseudo%Ncoeff
        write(chkptIO) trove%pseudo%field
        write(chkptIO) trove%pseudo%iorder
        !
        if (FLl2_coeffs) then
          !
          write(chkptIO) 'L2_vib'
          write(chkptIO) trove%L2_vib(1,1)%Ncoeff
          !
          do k1 = 1,Nmodes
             do k2 = 1,Nmodes
                !
                write(chkptIO) trove%L2_vib(k1,k2)%field
                write(chkptIO) trove%L2_vib(k1,k2)%iorder
                !
             enddo
          enddo
        endif
        !
        !
        ! Potential energy part
        !
        write(chkptIO) 'poten'
        write(chkptIO) trove%poten%Ncoeff
        write(chkptIO) trove%poten%field
        write(chkptIO) trove%poten%iorder
        !
        write(chkptIO) 'End Hamiltonian objects'
        !
        if (trim(trove%IO_ext_coeff)=='SAVE') then
          !
          write(chkptIO) 'extF'
          write(chkptIO) trove%extF(1)%Ncoeff
          !
          do k1 = 1,extF%rank
             !
             write(chkptIO) trove%extF(k1)%field
             write(chkptIO) trove%extF(k1)%iorder
             !
          enddo
          !
          write(chkptIO) 'End External object'
          !
        endif
        !
        close(chkptIO,status='keep')
        !
      end subroutine HamiltonianSave

      !
      ! this is the original version of HamiltonianSave uspporting a differnt, old format of hamiltonian.chk
      !
      subroutine HamiltonianSave_2007

        character(len=cl)  :: unitfname
        integer(ik)        :: Nmodes,k1,k2,chkptIO
        type(FLpolynomT),pointer    :: fl
        logical     :: i_opened
        !
        if (job%verbose>=3) write(out,"(/'Store all objects defining the Hamiltonian...')")
        !
        unitfname ='Check point of the Hamiltonian'
        call IOStart(trim(unitfname),chkptIO)
        !
        inquire (chkptIO,opened=i_opened)
        !
        if (i_opened) close(chkptIO)
        !
        open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=trove%chk_hamil_fname)
        !
        write(chkptIO) 'Start Hamiltonian objects'
        !
        Nmodes = trove%Nmodes
        !
        !
        write(chkptIO) 'Amatrho'
        write(chkptIO) trove%Amatrho
        !
        write(chkptIO) 'dAmatrho'
        write(chkptIO) trove%dAmatrho
        !
        write(chkptIO) 'Bmatrho'
        write(chkptIO) trove%Bmatrho
        !
        write(chkptIO) 'dBmatrho'
        write(chkptIO) trove%dBmatrho
        !
        write(chkptIO) 'g_vib'
        write(chkptIO) trove%g_vib(1,1)%Ncoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,Nmodes
              !
              write(chkptIO) trove%g_vib(k1,k2)%field
              write(chkptIO) trove%g_vib(k1,k2)%iorder
              !
           enddo
        enddo
        !
        write(chkptIO) 'g_rot'
        write(chkptIO) trove%g_rot(1,1)%Ncoeff
        !
        do k1 = 1,3
           do k2 = 1,3
              !
              write(chkptIO) trove%g_rot(k1,k2)%field
              write(chkptIO) trove%g_rot(k1,k2)%iorder
              !
           enddo
        enddo
        !
        write(chkptIO) 'g_cor'
        write(chkptIO) trove%g_cor(1,1)%Ncoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,3
              !
              write(chkptIO) trove%g_cor(k1,k2)%field
              write(chkptIO) trove%g_cor(k1,k2)%iorder
              !
           enddo
        enddo
        !
        fl => trove%pseudo
        !
        write(chkptIO) 'pseudo'
        write(chkptIO) trove%pseudo%Ncoeff
        write(chkptIO) trove%pseudo%field
        write(chkptIO) trove%pseudo%iorder
        !
        if (FLl2_coeffs) then
          !
          write(chkptIO) 'L2_vib'
          write(chkptIO) trove%L2_vib(1,1)%Ncoeff
          !
          do k1 = 1,Nmodes
             do k2 = 1,Nmodes
                !
                write(chkptIO) trove%L2_vib(k1,k2)%field
                write(chkptIO) trove%L2_vib(k1,k2)%iorder
                !
             enddo
          enddo
        endif
        !
        write(chkptIO) 'poten'
        write(chkptIO) trove%poten%Ncoeff
        write(chkptIO) trove%poten%field
        write(chkptIO) trove%poten%iorder
        !
        write(chkptIO) 'End Hamiltonian objects'
        !
        if (trim(trove%IO_ext_coeff)/='NONE') then
          !
          write(chkptIO) 'extF'
          write(chkptIO) trove%extF(1)%Ncoeff
          !
          do k1 = 1,extF%rank
             !
             write(chkptIO) trove%extF(k1)%field
             write(chkptIO) trove%extF(k1)%iorder
             !
          enddo
          !
          write(chkptIO) 'End External object'          
          !
        endif
        !
        close(chkptIO,status='keep')
        !
      end subroutine HamiltonianSave_2007
      !

      subroutine KineticSave_ASCII

        character(len=cl)  :: unitfname
        integer(ik)        :: Nmodes,k1,k2,chkptIO_kin,i,iterm,npoints,Ncoeff
        type(FLpolynomT),pointer    :: fl
        logical     :: i_opened
        !
        if (job%verbose>=3) write(out,"(/'Store all objects as ASCII ...')")
        !
        if (trim(trove%IO_kinetic)/='SAVE') return
        !
        unitfname ='Check point of the kinetic'
        call IOStart(trim(unitfname),chkptIO_kin)
        inquire (chkptIO_kin,opened=i_opened)
        if (i_opened) close(chkptIO_kin)
        open(chkptIO_kin,action='write',position='rewind',status='replace',file=trove%chk_kinet_fname)
        !
        Ncoeff = trove%RangeOrder(trove%NKinOrder)
        !
        write(chkptIO_kin,"(3i9,10x,a)") trove%g_vib(1,1)%Npoints,trove%g_vib(1,1)%orders,Ncoeff,"<- g_vib Npoints,Norder,Ncoeff"
        Nmodes = trove%Nmodes
        !
        do k1 = 1,Nmodes
          do k2 = 1,Nmodes
            !
            fl => trove%g_vib(k1,k2) 
            !
            call write_ascii(k1,k2,fl%Ncoeff,fl%Npoints,chkptIO_kin,fl%ifromsparse,fl%field)
            !
          enddo
        enddo
        !
        write(chkptIO_kin,"(i11,1x,i5,1x,i8,1x,i8,1x,e15.8,' <- End')") 987654321,0,0,0,0.0_ark
        !
        write(chkptIO_kin,"(3i9,10x,a)") trove%g_rot(1,1)%Npoints,trove%g_rot(1,1)%orders,Ncoeff,"<- g_rot Npoints,Norder,Ncoeff"
        !
        do k1 = 1,3
          do k2 = 1,3
            !
            fl => trove%g_rot(k1,k2) 
            !
            call write_ascii(k1,k2,fl%Ncoeff,fl%Npoints,chkptIO_kin,fl%ifromsparse,fl%field)
            !
          enddo
        enddo
        !
        write(chkptIO_kin,"(i11,1x,i5,1x,i8,1x,i8,1x,e18.11,' <- End')") 987654321,0,0,0,0.0_ark
        !
        write(chkptIO_kin,"(3i9,10x,a)") trove%g_cor(1,1)%Npoints,trove%g_cor(1,1)%orders,Ncoeff,"<- g_cor Npoints,Norder,Ncoeff"
        !
        do k1 = 1,Nmodes
          do k2 = 1,3
            !
            fl => trove%g_cor(k1,k2) 
            !
            call write_ascii(k1,k2,fl%Ncoeff,fl%Npoints,chkptIO_kin,fl%ifromsparse,fl%field)
            !
          enddo
        enddo
        !
        write(chkptIO_kin,"(i11,1x,i5,1x,i8,1x,i8,1x,e18.11,' <- End')") 987654321,0,0,0,0.0_ark
        !
        fl => trove%pseudo
        !
        write(chkptIO_kin,"(3i9,10x,a)") trove%pseudo%Npoints,trove%pseudo%Orders,Ncoeff,"<- pseudo Npoints,Norder,Ncoeff"
        !
        call write_ascii(0_ik,0_ik,fl%Ncoeff,fl%Npoints,chkptIO_kin,fl%ifromsparse,fl%field)
        !
        write(chkptIO_kin,"(i11,1x,i5,1x,i8,1x,i8,1x,e18.11,' <- End')") 987654321,0,0,0,0.0_ark
        !
        if (FLl2_coeffs) then
          !
          Ncoeff = trove%RangeOrder(2)
          !
          write(chkptIO_kin,"(3i9,10x,a)") trove%L2_vib(1,1)%Npoints,trove%L2_vib(1,1)%Orders,Ncoeff,&
                                           "<- L2_vib Npoints,Norder,Ncoeff"
          !
          do k1 = 1,Nmodes
            do k2 = 1,Nmodes
              !
              fl => trove%L2_vib(k1,k2) 
              !
              call write_ascii(k1,k2,fl%Ncoeff,fl%Npoints,chkptIO_kin,fl%ifromsparse,fl%field)
              !
            enddo
          enddo
          !
          write(chkptIO_kin,"(i11,1x,i5,1x,i8,1x,i8,1x,e18.11,' <- End')") 987654321,0,0,0,0.0_ark
          !
        endif
        ! 
        ! record indicating and/or specifying the sparse threshold used for expansion 
        ! coefficients 
        !
        if (trove%sparse) then 
          write(chkptIO_kin,"(e15.8,5x,' <- sparse threshold used')") job%exp_coeff_thresh
        else
          write(chkptIO_kin,"(e15.8,5x,' <- no sparse threshold was used')") 0.0d0
        endif
        !
        write(chkptIO_kin,"(a)") 'End of kinetic'
        !
        close(chkptIO_kin,status='keep')
        !
      end subroutine KineticSave_ASCII
      !
      !
      subroutine PotentialSave_ASCII

        character(len=cl)  :: unitfname
        integer(ik)        :: Nmodes,chkptIO_pot,i,iterm,npoints,Ncoeff
        type(FLpolynomT),pointer    :: fl
        logical     :: i_opened
        !
        if (job%verbose>=3) write(out,"(/'Store all objects as ASCII ...')")
        !
        Ncoeff = trove%RangeOrder(trove%NPotOrder)
        !
        unitfname ='Check point of the potential'
        call IOStart(trim(unitfname),chkptIO_pot)
        inquire (chkptIO_pot,opened=i_opened)
        if (i_opened) close(chkptIO_pot)
        open(chkptIO_pot,action='write',position='rewind',status='replace',file=trove%chk_poten_fname)
        !
        write(chkptIO_pot,"(3i9,10x,a)") trove%poten%Npoints,trove%poten%Orders,Ncoeff,"<- Npoints,Norder,Ncoeff"
        !
        fl => trove%poten
        !
        do iterm = 1,fl%Ncoeff
          !
          do i = 0,fl%Npoints
            !
            !if (abs(fl%field(iterm,i))>job%exp_coeff_thresh) then
              !
              write(chkptIO_pot,"(i8,1x,i8,1x,e23.16)") fl%ifromsparse(iterm),i,fl%field(iterm,i)
              !
            !endif
            !
          enddo
          !
        enddo
        !
        write(chkptIO_pot,"(i11,1x,i5,e18.11,' <- End')") 987654321,0,0.0_ark
        !
        if (trove%sparse) then 
          write(chkptIO_pot,"(e15.8,5x,' <- sparse threshold used')") job%exp_coeff_thresh
        else
          write(chkptIO_pot,"(e15.8,5x,' <- no sparse threshold was used')") 0.0d0
        endif
        !
        write(chkptIO_pot,"(a)") 'End of potential'
        !
        close(chkptIO_pot,status='keep')
        !
      end subroutine PotentialSave_ASCII
      !
      !
      subroutine ExternalSave_ASCII

        character(len=cl)  :: unitfname
        integer(ik)        :: Nmodes,k1,k2,chkptIO_ext,i,iterm,npoints,Ncoeff
        type(FLpolynomT),pointer    :: fl
        logical     :: i_opened
        !
        if (job%verbose>=3) write(out,"(/'Store all objects as ASCII ...')")
        !
        if (.not.trove%separate_convert.and.trim(trove%IO_ext_coeff)/='SAVE') return
        !
        Ncoeff = trove%RangeOrder(trove%NExtOrder)
        !
        unitfname ='Check point of the external'
        call IOStart(trim(unitfname),chkptIO_ext)
        inquire (chkptIO_ext,opened=i_opened)
        if (i_opened) close(chkptIO_ext)
        open(chkptIO_ext,action='write',position='rewind',status='replace',file=trove%chk_external_fname)
        !
        write(chkptIO_ext,"(4i9,10x,a)") trove%extF(1)%Npoints,trove%extF(1)%Orders,Ncoeff,extF%rank,"<- Npoints,Norder,Ncoeff,Rank"
        !
        do k1 = 1,extF%rank
          !
          do iterm = 1,trove%extF(k1)%Ncoeff
            !
            fl => trove%extF(k1)
            !
            do i = 0,trove%extF(1)%Npoints
              !
              !if (abs(fl%field(iterm,i))>small_) then
              !
              write(chkptIO_ext,"(i8,1x,i8,1x,i8,1x,e23.16)") k1,fl%ifromsparse(iterm),i,fl%field(iterm,i)
              !
              !endif
              !
            enddo
            !
          enddo
          !
        enddo
        !
        write(chkptIO_ext,"(i11,1x,i5,1x,i8,1x,e18.11,' <- End')") 987654321,0,0,0.0_ark
        !
        if (trove%sparse) then 
          write(chkptIO_ext,"(e15.8,5x,' <- sparse threshold used')") job%exp_coeff_thresh
        else
          write(chkptIO_ext,"(e15.8,5x,' <- no sparse was threshold used')") 0.0d0
        endif
        !
        write(chkptIO_ext,"(a)") 'End of external'
        !
        close(chkptIO_ext,status='keep')
        !
      end subroutine ExternalSave_ASCII
      !
      subroutine write_ascii(k1,k2,Ncoeff,Npoints,chkptIO_kin,ifromsparse,field)
        !
        integer(ik),intent(in)   :: k1,k2,Ncoeff,Npoints,chkptIO_kin
        real(ark),intent(in)     :: field(1:Ncoeff,0:Npoints)
        integer(ik),intent(in)   :: ifromsparse(1:Ncoeff)
        integer(ik) :: iterm,i
          !
          do iterm = 1,Ncoeff
            do i = 0,Npoints
              !
              if (abs(field(iterm,i))>job%exp_coeff_thresh) then 
                !
                write(chkptIO_kin,"(i5,1x,i5,1x,i8,1x,i8,1x,e26.18)") k1,k2,ifromsparse(iterm),i,real(field(iterm,i),rk)
                !
              endif
              !
             enddo
          enddo
          !
      end subroutine write_ascii
      !
      subroutine checkpointRestore_kinetic

        character(len=15) :: buf
        character(len=25) :: buf25
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,alloc,Tcoeff
        type(FLpolynomT),pointer    :: fl
        integer(ik)          :: Natoms,Nmodes,Npoints,k1,k2
        real(rk)             :: factor
        !
        integer(ik) :: n(2,8), m(2,8) , l(2,8), i, iterm, k(trove%Nmodes)
        !
        if (trove%separate_store.and.trim(trove%internal_coords)=='LOCAL') return
        !
        unitfname ='Check point of the Hamiltonian'
        call IOStart(trim(unitfname),chkptIO)
        open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=trove%chk_hamil_fname)
        !
        read(chkptIO) buf25
        if (buf25/='Start Hamiltonian objects') then
          write (out,"(' Checkpoint file ',a,' has bogus header: ',a)") trove%chk_hamil_fname, buf
          stop 'check_point_Hamiltonian - bogus file format (1)'
        end if
        !
        Natoms = trove%Natoms
        Nmodes = trove%Nmodes
        Npoints = trove%Npoints
        !
        if (.not.associated(trove%Amatrho).or..not.associated(trove%dAmatrho).or. &
            .not.associated(trove%Bmatrho).or..not.associated(trove%dBmatrho)) then 
           !
           write (out,"('basisRestore:  Amatrho-fields have to be allocated by now; maybe FLsetMolecule was no run yet')") 
           stop 'basisRestore, Amatrho-Bmatrho fields have to been alllocated before'
           !
        endif 
        !
        read(chkptIO) buf(1:7)
        if (buf(1:7)/='Amatrho') then
          write (out,"(' Checkpoint file ',a,' has bogus label Amatrho ',a)") trove%chk_fname, buf(1:7)
          stop 'check_point_Hamiltonian - bogus file format Amatrho'
        end if
        !
        read(chkptIO) trove%Amatrho
        !
        read(chkptIO) buf(1:8)
        if (buf(1:8)/='dAmatrho') then
          write (out,"(' Checkpoint file ',a,' has bogus label dAmatrho ',a)") trove%chk_fname, buf(1:8)
          stop 'check_point_Hamiltonian - bogus file format dAmatrho'
        end if
        !
        read(chkptIO) trove%dAmatrho
        !
        read(chkptIO) buf(1:7)
        if (buf(1:7)/='Bmatrho') then
          write (out,"(' Checkpoint file ',a,' has bogus label Bmatrho ',a)") trove%chk_fname, buf(1:7)
          stop 'check_point_Hamiltonian - bogus file format Bmatrho'
        end if
        !
        read(chkptIO) trove%Bmatrho
        !
        read(chkptIO) buf(1:8)
        if (buf(1:8)/='dBmatrho') then
          write (out,"(' Checkpoint file ',a,' has bogus label dBmatrho ',a)") trove%chk_fname, buf(1:8)
          stop 'check_point_Hamiltonian - bogus file format dBmatrho'
        end if
        !
        read(chkptIO) trove%dBmatrho
        !
        if (.not.associated(trove%g_vib).or..not.associated(trove%g_rot).or. &
            .not.associated(trove%g_cor)) then 
           !
           !write (out,"('basisRestore:  g-fields have to be allocated by now; maybe _checkpointRestore_kinetic_ was no run yet')") 
           !stop 'basisRestore, g-fields has to be alllocated before'
           !
           allocate (trove%g_vib(Nmodes,Nmodes),trove%g_rot(3,3),trove%g_cor(Nmodes,3),trove%pseudo,stat=alloc)
           if (alloc/=0) then
               write (out,"('chk_Restore_kin-Error ',i9,' trying to allocate g-fields')") alloc
               stop 'chk_Restore_kin, g-fields - out of memory'
           end if
           !
        endif 
        !
        if (FLl2_coeffs.and..not.associated(trove%L2_vib)) then 
           !
           allocate (trove%L2_vib(Nmodes,Nmodes),stat=alloc)
           if (alloc/=0) then
               write (out,"('chk_Restore_kin-Error ',i9,' trying to allocate L2_vib-field')") alloc
               stop 'chk_Restore_kin, L2_vib-field - out of memory'
           end if
           !
        endif
        !
        ! read from separate file in ASCII format if required 
        !
        if (trove%separate_store) return
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='g_vib') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_vib ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format g_vib'
        end if
        !
        read(chkptIO) Tcoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,Nmodes
              !
              fl => trove%g_vib(k1,k2)
              fl%Ncoeff = Tcoeff
              !
              call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_vib')
              read(chkptIO) fl%field     
              read(chkptIO) fl%iorder 
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='g_rot') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_rot ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format g_rot'
        end if
        !
        read(chkptIO) Tcoeff
        !
        do k1 = 1,3
           do k2 = 1,3
              !
              fl => trove%g_rot(k1,k2)
              fl%Ncoeff = Tcoeff
              !
              call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_rot')
              read(chkptIO) fl%field  
              read(chkptIO) fl%iorder
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='g_cor') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_cor ',a)") trove%chk_fname, buf(1:6)
          stop 'check_point_Hamiltonian - bogus file format g_cor'
        end if
        !
        read(chkptIO) Tcoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,3
              !
              fl => trove%g_cor(k1,k2)
              fl%Ncoeff = Tcoeff
              !
              call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_cor')
              read(chkptIO) fl%field  
              read(chkptIO) fl%iorder 
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:6)
        if (buf(1:6)/='pseudo') then
          write (out,"(' Checkpoint file ',a,' has bogus label poten ',a)") trove%chk_fname, buf(1:6)
          stop 'check_point_Hamiltonian - bogus file format poten'
        end if
        !
        read(chkptIO) Tcoeff
        !
        call polynom_initialization(trove%pseudo,trove%NKinOrder,Tcoeff,Npoints,'pseudo')
        !
        fl => trove%pseudo
        !
        read(chkptIO) fl%field
        read(chkptIO) fl%iorder
        !
        if (FLl2_coeffs) then 
          !
          read(chkptIO) buf(1:1)
          !
          backspace(chkptIO)
          !
          if (buf(1:1)/='L'.and..false.) then
            !
            if (trove%lincoord==0.or.trove%Nmodes/=7) then 
              !
              write (out,"(' Checkpoint file ',a,' has bogus label L2_vib ',a)") trove%chk_fname, buf(1:6)
              stop 'check_point_Hamiltonian - bogus file format L2_vib'
              !
            endif
            !
            ! THIS IS A TEMPORAL HACK FOR C2H2 CURVELINEAR COORDINATES 
            !
            write (out,"(' THIS IS A TEMPORAL HACK FOR C2H2 CURVELINEAR COORDINATES')")
            !
            do k1 = 1,Nmodes
               do k2 = 1,Nmodes
                  !
                  fl => trove%L2_vib(k1,k2)
                  fl%Ncoeff = Tcoeff
                  !
                  call polynom_initialization(fl,2,Tcoeff,Npoints,'L2_vib')
                  !
                  fl%field = 0
                  fl%iorder = 0
                  !
               enddo
            enddo
            !
            factor = 1.0_rk 
            !
            !real(planck,ark)*real(avogno,ark)*real(1.0d+16,kind=ark)/(4.0_ark*pi*pi*real(vellgt,ark))
            !
            n(1,1) = 4 ; n(2,1) = 5 ; m(1,1) = 5 ; m(2,1) = 5 ; l(1,1) = 2 ; l(2,1) = 0 
            n(1,2) = 4 ; n(2,2) = 5 ; m(1,2) = 4 ; m(2,2) = 4 ; l(1,2) = 0 ; l(2,2) = 2 
            n(1,3) = 4 ; n(2,3) = 5 ; m(1,3) = 4 ; m(2,3) = 5 ; l(1,3) = 1 ; l(2,3) = 1
                                      
            n(1,4) = 6 ; n(2,4) = 7 ; m(1,4) = 7 ; m(2,4) = 7 ; l(1,4) = 2 ; l(2,4) = 0 
            n(1,5) = 6 ; n(2,5) = 7 ; m(1,5) = 6 ; m(2,5) = 6 ; l(1,5) = 0 ; l(2,5) = 2
            n(1,6) = 6 ; n(2,6) = 7 ; m(1,6) = 6 ; m(2,6) = 7 ; l(1,6) = 1 ; l(2,6) = 1 
            !
            do iterm = 1,fl%Ncoeff
               !
               k(:) = FLIndexQ(:,iterm)
               !
               if ( sum(k(:))/=2 ) cycle
               !
               do i = 1,6 
                 ! 
                 if ( k(n(1,i))==l(1,i).and.k(n(2,i))==l(2,i) ) then 
                    trove%L2_vib(m(1,i),m(2,i))%field  = -factor
                    if ( l(1,i) == l(2,i) ) then 
                      trove%L2_vib(m(1,i),m(2,i))%field  = factor
                      trove%L2_vib(m(2,i),m(1,i))%field  = factor
                    endif
                 endif
                 !
               enddo
               !
            enddo
            !
          else
            !
            read(chkptIO) buf(1:6)
            !
            if (buf(1:6)/='L2_vib') then
               !
               write (out,"(' Checkpoint file ',a,' has bogus label L2_vib ',a)") trove%chk_fname, buf(1:6)
               stop 'check_point_Hamiltonian - bogus file format L2_vib'
               !
            endif
            !
            read(chkptIO) Tcoeff
            !
            do k1 = 1,Nmodes
               do k2 = 1,Nmodes
                  !
                  fl => trove%L2_vib(k1,k2)
                  fl%Ncoeff = Tcoeff
                  !
                  call polynom_initialization(fl,max(trove%NKinOrder,2),Tcoeff,Npoints,'L2_vib')
                  read(chkptIO) fl%field     
                  read(chkptIO) fl%iorder 
                  !
               enddo
            enddo
             !
          endif
          !
        endif
        !
        call MemoryReport
        !
      end subroutine checkpointRestore_kinetic


      subroutine checkpointRestore_kinetic_ascii

        character(len=14) :: buf
        character(len=25) :: buf25
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,alloc,Tcoeff
        type(FLpolynomT),pointer    :: fl
        integer(ik)          :: Natoms,Nmodes,Npoints,k1,k2,Tpoints,k1_,k2_,n,Torder,KinOrder
        real(rk)             :: factor
        real(ark)            :: field_
        real(rk)             :: exp_coeff_thresh

        !
        integer(ik) :: i, iterm, k(trove%Nmodes)
        !
        unitfname ='Check point of the kinetic'
        call IOStart(trim(unitfname),chkptIO)
        open(chkptIO,action='read',status='old',file=trove%chk_kinet_fname)
        !
        Natoms = trove%Natoms
        Nmodes = trove%Nmodes
        Npoints = trove%Npoints
        KinOrder = trove%NKinOrder
        !
        if (.not.associated(trove%g_vib).or..not.associated(trove%g_rot).or. &
            .not.associated(trove%g_cor)) then 
           !
           !write (out,"('basisRestore:  g-fields have to be allocated by now; maybe _checkpointRestore_kinetic_ascii_ was no run yet')") 
           !stop 'basisRestore, g-fields has to be alllocated before'
           !
           allocate (trove%g_vib(Nmodes,Nmodes),trove%g_rot(3,3),trove%g_cor(Nmodes,3),trove%pseudo,stat=alloc)
           if (alloc/=0) then
               write (out,"('chk_Restore_kin-Error ',i9,' trying to allocate g-fields')") alloc
               stop 'chk_Restore_kin, g-fields - out of memory'
           end if
           !
        endif 
        !
        if (FLl2_coeffs.and..not.associated(trove%L2_vib)) then 
           !
           allocate (trove%L2_vib(Nmodes,Nmodes),stat=alloc)
           if (alloc/=0) then
               write (out,"('chk_Restore_kin-Error ',i9,' trying to allocate L2_vib-field')") alloc
               stop 'chk_Restore_kin, L2_vib-field - out of memory'
           end if
           !
        endif
        !
        ! start reading 
        !
        read(chkptIO,*) Tpoints,Torder,Tcoeff
        !
        if (Tpoints/=Npoints) then
          write(out,"('Kinetic-ASCII-chk npoints is wrong:',2i8)") Tpoints,Npoints
          stop "Kinetic-ASCII-chk npoints is wrong"
        endif
        if (Torder/=KinOrder) then 
          write(out,"('Kinetic-ASCII-chk Norder is wrong:',2i8)") Torder,KinOrder
          stop "Kinetic-ASCII-chk Norder is wrong"
        endif
        !
        do k1 = 1,Nmodes
          do k2 = 1,Nmodes
            !
            fl => trove%g_vib(k1,k2)
            fl%Ncoeff = Tcoeff
            !
            call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_vib')
            !fl%ifromsparse = 0
            forall(n=1:Tcoeff) fl%ifromsparse(n) = n
            fl%sparse = .true.
            !
          enddo
        enddo
        !
        k1_= 0 ; k2_= 0 ; n = 0
        do_gvib : do 
           !
           read(chkptIO,*) k1,k2,iterm,i,field_
           !
           if (k1_/=k1.or.k2_/=k2) then
             k1_= k1 ; k2_= k2
             n = 0
           endif
           !
           if (k1==987654321) exit do_gvib
           !
           trove%g_vib(k1,k2)%field(iterm,i) = field_
           !
        enddo do_gvib
        !
        read(chkptIO,*) Tpoints,Torder,Tcoeff
        !
        if (Tpoints/=Npoints) then 
           print*,"grot-ASCII-chk npoints is wrong"
           stop "grot-ASCII-chk npoints is wrong"
        endif
        !
        if (Torder/=KinOrder) then
          print*,"grot-ASCII-chk Order is wrong"
          stop "grot-ASCII-chk Order is wrong"
        endif
        !
        do k1 = 1,3
          do k2 = 1,3
            !
            fl => trove%g_rot(k1,k2)
            fl%Ncoeff = Tcoeff
            !
            call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_rot')
            forall(n=1:Tcoeff) fl%ifromsparse(n) = n
            fl%sparse = .true.
            !
          enddo
        enddo
        !
        k1_= 0 ; k2_= 0 ; n = 0
        do_grot : do 
           !
           read(chkptIO,*) k1,k2,iterm,i,field_
           !
           if (k1_/=k1.or.k2_/=k2) then
             k1_= k1 ; k2_= k2
             n = 0
           endif
           !
           if (k1==987654321) exit do_grot
           !
           trove%g_rot(k1,k2)%field(iterm,i) = field_
           !
        enddo do_grot
        !
        read(chkptIO,*) Tpoints,Torder,Tcoeff
        !
        if (Tpoints/=Npoints) stop "gcor-ASCII-chk npoints is wrong"
        if (Torder/=KinOrder) stop "gcor-ASCII-chk Order is wrong"
        !
        do k1 = 1,Nmodes
          do k2 = 1,3
            !
            fl => trove%g_cor(k1,k2)
            fl%Ncoeff = Tcoeff
            !
            call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'g_cor')
            forall(n=1:Tcoeff) fl%ifromsparse(n) = n
            fl%sparse = .true.
            !
          enddo
        enddo
        !
        k1_= 0 ; k2_= 0 ; n = 0
        do_gcor : do 
           !
           read(chkptIO,*) k1,k2,iterm,i,field_
           !
           if (k1==987654321) exit do_gcor
           !
           if (k1_/=k1.or.k2_/=k2) then
             k1_= k1 ; k2_= k2
             n = 0
           endif
           !
           trove%g_cor(k1,k2)%field(iterm,i) = field_
           !
        enddo do_gcor
        !
        !
        read(chkptIO,*) Tpoints,Torder,Tcoeff
        !
        if (Tpoints/=Npoints) stop "pseudo-ASCII-chk npoints is wrong"
        if (Torder/=KinOrder) stop "pseudo-ASCII-chk Order is wrong"
        !
        fl => trove%pseudo
        fl%Ncoeff = Tcoeff
        !
        call polynom_initialization(fl,trove%NKinOrder,Tcoeff,Npoints,'pseudo')
        forall(n=1:Tcoeff) fl%ifromsparse(n) = n
        fl%sparse = .true.
        !
        n = 0
        do_pseu : do 
           !
           read(chkptIO,*) k1,k2,iterm,i,field_
           !
           if (k1==987654321) exit do_pseu
           !
           trove%pseudo%field(iterm,i) = field_
           !
        enddo do_pseu
        !
        if (FLl2_coeffs) then 
          !
          read(chkptIO,*) Tpoints,Torder,Tcoeff
          !
          if (Tpoints/=Npoints) stop "L2vib-ASCII-chk npoints is wrong"
          if (Torder/=2) then 
            write (out,"('L2vib-ASCII-chk Order is not 2',i8)") Torder
            stop "L2vib-ASCII-chk Order is wrong"
          endif
          !
          do k1 = 1,Nmodes
            do k2 = 1,Nmodes
              !
              fl => trove%L2_vib(k1,k2)
              fl%Ncoeff = Tcoeff
              !
              call polynom_initialization(fl,max(trove%NKinOrder,2),Tcoeff,Npoints,'L2_vib')
              forall(n=1:Tcoeff) fl%ifromsparse(1:n) = (/(n,n=1, Tcoeff)/)            
              fl%sparse = .true.
              !
            enddo
          enddo
          !
          k1_= 0 ; k2_= 0 ; n = 0
          do_L2vib : do 
             !
             read(chkptIO,*) k1,k2,iterm,i,field_
             !
             if (k1_/=k1.or.k2_/=k2) then
               k1_= k1 ; k2_= k2
               n = 0
             endif
             !
             if (k1==987654321) exit do_L2vib
             !
             trove%L2_vib(k1,k2)%field(iterm,i) = field_
             !
          enddo do_L2vib
          !
        endif
        !
        read(chkptIO,*) exp_coeff_thresh
        !
        if ( abs(exp_coeff_thresh-job%exp_coeff_thresh)>small_ ) then
           !
           write(out,"('WARNING: in kinetic.chk exp_coeff_thresh is inconsistent with used: ',2e18.10)") &
                     job%exp_coeff_thresh,exp_coeff_thresh
           !
        endif
        !
        read(chkptIO,"(a14)") buf
        !
        if (buf/='End of kinetic') then
          write (out,"(' Checkpoint file ',a,' has bogus label kinetic-ascii',a)") trove%chk_fname, buf
          stop 'check_point_Hamiltonian - bogus file format kinetic-ASCII'
        end if
        !
        call MemoryReport
        !
      end subroutine checkpointRestore_kinetic_ascii



      subroutine checkpointSkip_kinetic

        character(len=15) :: buf
        character(len=25) :: buf25
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,alloc,Tcoeff
        type(FLpolynomT),pointer    :: fl
        integer(ik)               :: Natoms,Nmodes,Npoints,k1,k2,isize
        real(rk),allocatable     :: array(:)
        integer(ik),allocatable  :: iarray(:)
        !
        unitfname ='Check point of the Hamiltonian'
        call IOStart(trim(unitfname),chkptIO)
        open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=trove%chk_hamil_fname)
        !
        read(chkptIO) buf25
        if (buf25/='Start Hamiltonian objects') then
          write (out,"(' Checkpoint file ',a,' has bogus header: ',a)") trove%chk_hamil_fname, buf
          stop 'check_point_Hamiltonian - bogus file format (1)'
        end if
        !
        Natoms = trove%Natoms
        Nmodes = trove%Nmodes
        Npoints = trove%Npoints
        !
        if (.not.associated(trove%Amatrho).or..not.associated(trove%dAmatrho).or. &
            .not.associated(trove%Bmatrho).or..not.associated(trove%dBmatrho)) then 
           !
           write (out,"('basisRestore:  Amatrho-fields have to be allocated by now; maybe FLsetMolecule was no run yet')") 
           stop 'basisRestore, Amatrho-Bmatrho fields have to been alllocated before'
           !
        endif 
        !
        read(chkptIO) buf(1:7)
        if (buf(1:7)/='Amatrho') then
          write (out,"(' Checkpoint file ',a,' has bogus label Amatrho ',a)") trove%chk_fname, buf(1:7)
          stop 'check_point_Hamiltonian - bogus file format Amatrho'
        end if
        !
        isize = size(trove%Amatrho)
        allocate(array(isize),stat=alloc)
        read(chkptIO) array
        deallocate(array)
        !
        read(chkptIO) buf(1:8)
        if (buf(1:8)/='dAmatrho') then
          write (out,"(' Checkpoint file ',a,' has bogus label dAmatrho ',a)") trove%chk_fname, buf(1:8)
          stop 'check_point_Hamiltonian - bogus file format dAmatrho'
        end if
        !
        isize = size(trove%dAmatrho)
        allocate(array(isize),stat=alloc)
        read(chkptIO) array
        deallocate(array)
        !
        read(chkptIO) buf(1:7)
        if (buf(1:7)/='Bmatrho') then
          write (out,"(' Checkpoint file ',a,' has bogus label Bmatrho ',a)") trove%chk_fname, buf(1:7)
          stop 'check_point_Hamiltonian - bogus file format Bmatrho'
        end if
        !
        isize = size(trove%Bmatrho)
        allocate(array(isize),stat=alloc)
        read(chkptIO) array
        deallocate(array)
        !
        read(chkptIO) buf(1:8)
        if (buf(1:8)/='dBmatrho') then
          write (out,"(' Checkpoint file ',a,' has bogus label dBmatrho ',a)") trove%chk_fname, buf(1:8)
          stop 'check_point_Hamiltonian - bogus file format dBmatrho'
        end if
        !
        isize = size(trove%dBmatrho)
        allocate(array(isize),stat=alloc)
        read(chkptIO) array
        deallocate(array)
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='g_vib') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_vib ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format g_vib'
        end if
        !
        isize = size(trove%g_vib(1,1)%field)
        allocate(array(isize),stat=alloc)
        !
        isize = size(trove%g_vib(1,1)%iorder)
        allocate(iarray(isize),stat=alloc)
        !
        read(chkptIO) Tcoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,Nmodes
              !
              read(chkptIO) array     
              read(chkptIO) iarray
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='g_rot') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_rot ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format g_rot'
        end if
        !
        read(chkptIO) Tcoeff
        !
        do k1 = 1,3
           do k2 = 1,3
              !
              read(chkptIO) array     
              read(chkptIO) iarray
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='g_cor') then
          write (out,"(' Checkpoint file ',a,' has bogus label g_cor ',a)") trove%chk_fname, buf(1:6)
          stop 'check_point_Hamiltonian - bogus file format g_cor'
        end if
        !
        read(chkptIO) Tcoeff
        !
        do k1 = 1,Nmodes
           do k2 = 1,3
              !
              read(chkptIO) array     
              read(chkptIO) iarray
              !
           enddo
        enddo
        !
        read(chkptIO) buf(1:6)
        if (buf(1:6)/='pseudo') then
          write (out,"(' Checkpoint file ',a,' has bogus label poten ',a)") trove%chk_fname, buf(1:6)
          stop 'check_point_Hamiltonian - bogus file format poten'
        end if
        !
        read(chkptIO) Tcoeff
        !
        read(chkptIO) array     
        read(chkptIO) iarray
        !
        if (FLl2_coeffs) then 
          !
          read(chkptIO) buf(1:1)
          !
          backspace(chkptIO)
          !
          if (buf(1:1)/='L') then
            !
            if (trove%lincoord==0.or.trove%Nmodes/=7) then 
              !
              write (out,"(' Checkpoint file ',a,' has bogus label L2_vib ',a)") trove%chk_fname, buf(1:6)
              stop 'check_point_Hamiltonian - bogus file format L2_vib'
              !
            endif
            !
          else
            !
            read(chkptIO) buf(1:6)
            !
            if (buf(1:6)/='L2_vib') then
               !
               write (out,"(' Checkpoint file ',a,' has bogus label L2_vib ',a)") trove%chk_fname, buf(1:6)
               stop 'check_point_Hamiltonian - bogus file format L2_vib'
               !
            endif
            !
            read(chkptIO) Tcoeff
            !
            do k1 = 1,Nmodes
               do k2 = 1,Nmodes
                  !
                  read(chkptIO) array     
                  read(chkptIO) iarray
                  !
               enddo
            enddo
             !
          endif
          !
        endif
        !
        deallocate(array)
        deallocate(iarray)
        !
        call MemoryReport
        !
      end subroutine checkpointSkip_kinetic



      !
      subroutine checkpointRestore_potential

        character(len=23) :: buf
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,alloc
        type(FLpolynomT),pointer    :: fl
        integer(ik)          :: Nmodes,Npoints,Ncoeff
        !
        unitfname ='Check point of the Hamiltonian'
        call IOStart(trim(unitfname),chkptIO)
        !
        read(chkptIO) buf(1:5)
        if (buf(1:5)/='poten') then
          write (out,"(' Checkpoint file ',a,' has bogus label poten ',a)") trove%chk_fname, buf(1:5)
          stop 'check_point_Hamiltonian - bogus file format poten'
        end if
        !
        allocate (trove%poten,stat=alloc)        
        !
        read(chkptIO) Ncoeff
        !
        call polynom_initialization(trove%poten,trove%NPotOrder,Ncoeff,trove%Npoints,'poten')
        !
        fl => trove%poten
        !
        read(chkptIO) fl%field
        read(chkptIO) fl%iorder
        !
        read(chkptIO) buf
        !
        if (buf/='End Hamiltonian objects') then
          write (out,"(' Checkpoint file ',a,' has bogus label poten ',a)") trove%chk_fname, buf
          stop 'check_point_Hamiltonian - bogus file format poten'
        end if
        !
        if (trim(trove%IO_ext_coeff)/='READ') close(chkptIO,status='keep')
        !
      end subroutine checkpointRestore_potential
      !

      !
      subroutine checkpointRestore_potential_ascii

        character(len=16) :: buf
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,alloc,k
        type(FLpolynomT),pointer    :: fl
        integer(ik)          :: Nmodes,Npoints,Ncoeff,iterm,i,Norder
        real(ark) :: field_
        real(rk)             :: exp_coeff_thresh
        !

        unitfname ='Check point of the potential'
        call IOStart(trim(unitfname),chkptIO)
        open(chkptIO,action='read',status='old',file=trove%chk_poten_fname)
        !
        if (.not.associated(trove%poten)) then 
           !
           allocate (trove%poten,stat=alloc) 
           if (alloc/=0) then
               write (out,"('chk_Restore_pot_ascii-Error ',i9,' trying to allocate poten-field')") alloc
               stop 'chk_Restore_pot_ascii, poten-field - out of memory'
           end if
           !
        endif
        !
        ! start reading 
        !
        read(chkptIO,*) Npoints,Norder,Ncoeff
        !
        if (Npoints/=trove%Npoints) then
          write(out,"('poten-ASCII-chk npoints is wrong:',2i8)") Npoints,Npoints
          stop "poten-ASCII-chk npoints is wrong"
        endif
        !
        if (Norder/=trove%NPotOrder) then 
          write(out,"('poten-ASCII-chk Norder is wrong:',2i8)") Norder,trove%NPotOrder
          stop "poten-ASCII-chk Norder is wrong"
        endif
        !
        fl => trove%poten
        fl%Ncoeff = Ncoeff
        !
        call polynom_initialization(trove%poten,trove%NPotOrder,Ncoeff,trove%Npoints,'poten')
        forall(k=1:Ncoeff) fl%ifromsparse(k) = k
        fl%sparse = .true.
        k = 0 
        !
        do_pot : do 
           !
           read(chkptIO,*) iterm,i,field_
           !
           if (iterm==987654321) exit do_pot
           !
           fl%field(iterm,i) = field_
           !
        enddo do_pot
        !
        read(chkptIO,*) exp_coeff_thresh
        !
        if ( abs(exp_coeff_thresh-job%exp_coeff_thresh)>small_ ) then
           !
           write(out,"('WARNING: in potential.chk exp_coeff_thresh is inconsistent with used: ',2e18.10)") &
                       job%exp_coeff_thresh,exp_coeff_thresh
           !
        endif
        !
        read(chkptIO,"(a16)") buf
        !
        if (buf/='End of potential') then
          write (out,"(' Checkpoint file ',a,' has bogus label poten-ascii',a)") trove%chk_fname, buf
          stop 'check_point_Hamiltonian - bogus file format poten-ASCII'
        end if
        !
      end subroutine checkpointRestore_potential_ascii





      !
      subroutine checkpointRestore_external

        character(len=19) :: buf
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,alloc,imu,nterms
        type(FLpolynomT),pointer    :: fl
        integer(ik)          :: Nmodes,Npoints,Ncoeff
        real(rk)             :: exp_coeff_thresh
        !
        unitfname ='Check point of the Hamiltonian'
        call IOStart(trim(unitfname),chkptIO)
        !
        read(chkptIO) buf(1:4)
        if (buf(1:4)/='extF') then
          write (out,"(' Checkpoint file ',a,' has bogus label extF ',a)") trove%chk_fname, buf(1:4)
          stop 'check_point_Hamiltonian - bogus file format extF'
        end if
        !
        if (.not.associated(trove%extF)) then 
           !
           allocate (trove%extF(extF%rank),stat=alloc)
           if (alloc/=0) then
               write (out,"(' Error ',i9,' trying to allocate extF-field')") alloc
               stop 'chk_Restore, extF-field - out of memory'
           end if
           !
        endif 
        !
        read(chkptIO) nterms
        !
        do imu = 1,extF%rank
           !
           fl => trove%extF(imu)
           fl%Ncoeff = nterms
           !
           call polynom_initialization(fl,extF%maxord(imu),nterms,trove%Npoints,'extF')
           read(chkptIO) fl%field     
           read(chkptIO) fl%iorder 
           !
        enddo
        !
        !read(chkptIO,*) exp_coeff_thresh
        !
        !if ( abs(exp_coeff_thresh-job%exp_coeff_thresh)>small_ ) then
        !  !
        !  if (.not.trove%sparse) then
        !    write(out,"('external.chk: A sparse option was used to save kinetic.chk, please switch SPARSE on by adding SPARSE to the inpiut')")
        !    stop 'external.chk error SPARSE should be switched on'
        !  endif
        !  !
        !  if (job%exp_coeff_thresh>exp_coeff_thresh.and.trim(trove%IO_basisset)=="READ") then
        !    write(out,"('external.chk: A larger exp_coeff_thresh ',e18.10,' is incompatible with BASIS ',e18.10,', change threshold or redo BASIS_SET SAVE')") & 
        !                job%exp_coeff_thresh>exp_coeff_thresh
        !    stop 'external.chk: A larger exp_coeff_thresh is incompatible with BASIS, change threshold or BASIS_SET SAVE'
        !  endif
        !  !
        !endif                
        !
        read(chkptIO) buf
        !
        if (buf/='End External object') then
          write (out,"(' Checkpoint file ',a,' has bogus label extF ',a)") trove%chk_fname, buf
          stop 'check_point_Hamiltonian - bogus file format extF'
        end if
        !
        close(chkptIO,status='keep')
        !
      end subroutine checkpointRestore_external
      !      

      subroutine checkpointRestore_external_ascii
      

        character(len=15) :: buf
        character(len=cl)  :: unitfname
        integer(ik)        :: chkptIO,alloc
        type(FLpolynomT),pointer    :: fl
        integer(ik)          :: Nmodes,Npoints,Ncoeff,iterm,i,imu,nterms,nrank,Norder
        real(ark) :: field_
        integer(ik),allocatable :: nterm(:)
        real(rk)             :: exp_coeff_thresh
        !
        unitfname ='Check point of the external'
        call IOStart(trim(unitfname),chkptIO)
        open(chkptIO,action='read',status='old',file=trove%chk_external_fname)
        !
        if (.not.associated(trove%extF)) then 
           !
           allocate (trove%extF(extF%rank),stat=alloc)
           if (alloc/=0) then
               write (out,"(' Error ',i9,' trying to allocate extF-field')") alloc
               stop 'chk_Restore, extF-field - out of memory'
           end if
           !
        endif 
        !
        ! start reading 
        !
        read(chkptIO,*) Npoints,Norder,nterms,nrank
        !
        if (nrank/=extF%rank) then 
           write(out,"(a,i8,1x,i8)") "external-ASCII-chk rank is wrong: " ,nrank,extF%rank
           stop "external-ASCII-chk rank is wrong"
        endif
        !
        if (Npoints/=trove%Npoints) stop "external-ASCII-chk npoints is wrong"
        !
        if (Npoints/=trove%Npoints) then
          write(out,"('external-ASCII-chk npoints is wrong:',2i8)") Npoints,Npoints
          stop "external-ASCII-chk npoints is wrong"
        endif
        !
        if (Norder/=trove%NExtOrder) then 
          write(out,"('external-ASCII-chk Norder is wrong:',2i8)") Norder,trove%NExtOrder
          stop "external-ASCII-chk Norder is wrong"
        endif
        !
        if (nrank/=extF%rank) then 
          write(out,"('external-ASCII-chk nrank is wrong:',2i8)") nrank,extF%rank
          stop "external-ASCII-chk nrank is wrong"
        endif
        !
        do imu = 1,extF%rank
           !
           fl => trove%extF(imu)
           fl%Ncoeff = nterms
           !
           call polynom_initialization(fl,extF%maxord(imu),nterms,trove%Npoints,'extF')
           forall(i=1:fl%Ncoeff) fl%ifromsparse(i) = i
           fl%sparse = .true.
           !
        enddo
        !
        allocate (nterm(extF%rank),stat=alloc)
        call ArrayStart('nterm(extF%rank)',alloc,size(nterm),kind(nterm))
        nterm = 0
        !
        do_ext : do 
           !
           read(chkptIO,*) imu,iterm,i,field_
           !
           if (imu==987654321) exit do_ext
           !
           trove%extF(imu)%field(iterm,i) = field_
           !
        enddo do_ext
        !
        read(chkptIO,*) exp_coeff_thresh
        !
        if ( abs(exp_coeff_thresh-job%exp_coeff_thresh)>small_ ) then
           !
           write(out,"('WARNING: in external.chk exp_coeff_thresh is inconsistent with used: ',2e18.10)") &
                       job%exp_coeff_thresh,exp_coeff_thresh
           !
        endif
        !
        read(chkptIO,"(a15)") buf
        !
        if (buf/='End of external') then
          write (out,"(' Checkpoint file ',a,' has bogus label external-ascii',a)") trove%chk_fname, buf
          stop 'check_point_Hamiltonian - bogus file format ext-ASCII'
        end if
        !
        deallocate (nterm)
        call ArrayStop('nterm(extF%rank)')
        !
      end subroutine checkpointRestore_external_ascii
      !
   end subroutine FLcheck_point_Hamiltonian
   !

   !
   ! A compact, sparse representation of a field 
   !
   subroutine FLCompact_a_field_sparse(fl,name)

     type(FLpolynomT),pointer  :: fl
     character(len=*),intent(in) :: name
     integer(ik)        :: Npoints,Ncoeff,iterm,i,icoeff,Nterms,alloc
     real(ark),allocatable    :: sfield(:,:)  ! Expansion parameters in the sparse representation
     integer(ik),allocatable  :: siorder(:)  ! iorder in sparse
     logical :: check = .true.
     !
     Ncoeff = fl%Ncoeff
     Npoints = fl%Npoints
     !
     ! Count large elements (using exp_coeff_thresh as threshold) and store in a sparse representation
     !
     iterm = 0
     !
     do icoeff = 1,Ncoeff
       if (any(abs(fl%field(icoeff,:))>job%exp_coeff_thresh)) then
          iterm = iterm + 1
       endif
     enddo
     !
     Nterms = iterm
     !
     ! Create a field in a sparse representaion
     !
     if (associated(fl%IndexQ)) then 
       deallocate(fl%IndexQ)
       call ArrayMinus(name//'IndexQ',isize=size(fl%IndexQ),ikind=kind(fl%IndexQ))
     endif
     !
     deallocate(fl%ifromsparse)
     call ArrayMinus(name//'ifromsparse',isize=size(fl%ifromsparse),ikind=kind(fl%ifromsparse))
     !
     if (Nterms==0) then
       !
       fl%Ncoeff = Nterms
       !
       call ArrayMinus(name,isize=size(fl%iorder),ikind=kind(fl%iorder))
       deallocate(fl%iorder)
       !
       call ArrayMinus(name,isize=size(fl%field),ikind=kind(fl%field))
       deallocate(fl%field,stat=alloc)
       !
       return
       !
     endif 
     !
     allocate(Sfield(Nterms,0:Npoints),fl%IndexQ(trove%Nmodes,Nterms),stat=alloc)
     call ArrayStart("Sfield",alloc,size(Sfield),kind(Sfield))
     call ArrayStart(name//'IndexQ',alloc,size(fl%IndexQ),kind(fl%IndexQ))
     !
     allocate(fl%ifromsparse(Nterms),stat=alloc)
     call ArrayStart(name//'ifromsparse',alloc,size(fl%ifromsparse),kind(fl%ifromsparse))
     !
     allocate(siorder(Nterms),stat=alloc)
     call ArrayStart("Sfield",alloc,size(siorder),kind(siorder))
     !
     iterm = 0
     Sfield = 0
     siorder = 0
     fl%ifromsparse = 0
     !
     do icoeff = 1,Ncoeff
       if (any(abs(fl%field(icoeff,:))>job%exp_coeff_thresh)) then
          !
          iterm = iterm + 1
          !
          Sfield(iterm,:) = fl%field(icoeff,:)
          fl%IndexQ(:,iterm) = FLIndexQ(:,icoeff)
          Siorder(iterm) = fl%iorder(icoeff)
          fl%ifromsparse(iterm) = icoeff
          !
       endif
     enddo
     !
     call ArrayMinus(name,isize=size(fl%iorder),ikind=kind(fl%iorder))
     deallocate(fl%iorder)
     !
     call ArrayMinus(name,isize=size(fl%field),ikind=kind(fl%field))
     deallocate(fl%field,stat=alloc)
     !
     allocate(fl%field(Nterms,0:Npoints),stat=alloc)
     call ArrayStart(name,alloc,size(fl%field),kind(fl%field))
     !
     allocate(fl%iorder(Nterms),stat=alloc)
     call ArrayStart(name,alloc,size(fl%iorder),kind(fl%iorder))
     !
     fl%field = Sfield
     fl%Ncoeff = Nterms
     fl%iorder = Siorder
     !
     deallocate(Sfield,siorder)
     !
     call ArrayStop("Sfield")
     !
   end subroutine FLCompact_a_field_sparse
   ! 
   !
   ! A compact, sparse representation of a field 
   !
   subroutine FLCompact_and_combine_fields_sparse(fl1,name1,fl2,name2)

     type(FLpolynomT),pointer  :: fl1,fl2
     character(len=*),intent(in) :: name1,name2
     integer(ik)        :: Npoints,Ncoeff1,Ncoeff2,iterm,i,icoeff,Nterms,alloc,Nterm1,Nterm2,Ncoeffmax
     real(ark),allocatable    :: sfield1(:,:),sfield2(:,:)  ! Expansion parameters in the sparse representation
     integer(ik),allocatable  :: siorder1(:),siorder2(:)      ! iorder in sparse
     integer(ik)   :: target_index(trove%Nmodes)
     logical :: check = .true.
     !
     Ncoeff1 = fl1%Ncoeff
     Ncoeff2 = fl2%Ncoeff
     !
     Npoints = fl1%Npoints
     !
     if (Npoints/=fl2%Npoints) then
       write(out,"('FLCompact_and_combine_fields_sparse: Illegal Npoints in two fields, should be the same',2i8)") & 
             fl1%Npoints,fl2%Npoints
       stop 'FLCompact_and_combine_fields_sparse: Illegal Npoints in two fields'
     endif
     !    
     target_index = 0
     target_index(1) = trove%NKinOrder 
     Ncoeffmax= FLQindex(trove%Nmodes_e,target_index)
     !
     ! Count large elements (using exp_coeff_thresh as threshold) and store in a sparse representation
     !
     iterm = 0
     !
     do icoeff = 1,Ncoeffmax
       if (any(abs(fl1%field(icoeff,:))>job%exp_coeff_thresh).or.any(abs(fl2%field(icoeff,:))>job%exp_coeff_thresh)) then
          iterm = iterm + 1
       endif
     enddo
     !
     nterms = iterm
     !
     if (Nterms==0) then
       !
       stop 'FLCompact_and_combine_fields_sparse is not implemented for Nterms=0' 
       !
       deallocate(fl1%IndexQ,fl2%IndexQ)
       call ArrayStop(name1//'IndexQ')
       call ArrayStop(name2//'IndexQ')
       !
       deallocate(fl1%ifromsparse,fl2%ifromsparse)
       call ArrayStop(name1//"ifromsparse")
       call ArrayStop(name2//"ifromsparse")
       !
       return
       !
     endif 
     !
     allocate(Sfield1(nterms,0:Npoints),Sfield2(nterms,0:Npoints),stat=alloc)
     call ArrayStart("Sfield",alloc,size(Sfield1),kind(Sfield1))
     call ArrayStart("Sfield",alloc,size(Sfield2),kind(Sfield2))
     !
     allocate(Siorder1(nterms),Siorder2(nterms),stat=alloc)
     call ArrayStart("Sfield",alloc,size(Siorder1),kind(Siorder1))
     call ArrayStart("Sfield",alloc,size(Siorder2),kind(Siorder2))
     !
     deallocate(fl1%IndexQ,fl2%IndexQ)
     call ArrayStop(name1//'IndexQ')
     call ArrayStop(name2//'IndexQ')
     !
     deallocate(fl1%ifromsparse,fl2%ifromsparse)
     call ArrayStop(name1//"ifromsparse")
     call ArrayStop(name2//"ifromsparse")
     !
     allocate(fl1%ifromsparse(nterms),fl1%IndexQ(trove%Nmodes,nterms),stat=alloc)
     allocate(fl2%ifromsparse(nterms),fl2%IndexQ(trove%Nmodes,nterms),stat=alloc)
     call ArrayStart(name1//"ifromsparse",alloc,size(fl1%ifromsparse),kind(fl1%ifromsparse))
     call ArrayStart(name1//"IndexQ",alloc,size(fl1%IndexQ),kind(fl1%IndexQ))
     call ArrayStart(name2//"ifromsparse",alloc,size(fl2%ifromsparse),kind(fl2%ifromsparse))
     call ArrayStart(name2//"IndexQ",alloc,size(fl2%IndexQ),kind(fl2%IndexQ))
     !
     iterm = 0
     !
     do icoeff = 1,Ncoeffmax
       if (any(abs(fl1%field(icoeff,:))>job%exp_coeff_thresh).or.any(abs(fl2%field(icoeff,:))>job%exp_coeff_thresh)) then
          !
          iterm = iterm + 1
          !
          Sfield1(iterm,:) = fl1%field(icoeff,:)
          Sfield2(iterm,:) = fl2%field(icoeff,:)
          siorder1(iterm) = fl1%iorder(icoeff)
          siorder2(iterm) = fl2%iorder(icoeff)
          !
          fl1%ifromsparse(iterm) = icoeff
          fl1%IndexQ(:,iterm) = FLIndexQ(:,icoeff)
          !
          fl2%ifromsparse(iterm) = icoeff
          fl2%IndexQ(:,iterm) = FLIndexQ(:,icoeff)
          !
       endif
     enddo
     !
     deallocate(fl1%iorder,fl2%iorder)
     call ArrayStop(name1)
     call ArrayStop(name2)
     !
     ! Create a field in a sparse representaion
     !
     allocate(fl1%iorder(nterms),fl2%iorder(nterms),stat=alloc)
     call ArrayStart(name1,alloc,size(fl1%iorder),kind(fl1%iorder))
     call ArrayStart(name2,alloc,size(fl2%iorder),kind(fl2%iorder))
     !
     deallocate(fl1%field,fl2%field,stat=alloc)
     call ArrayMinus(name1,isize=size(fl1%field),ikind=kind(fl1%field))
     call ArrayMinus(name2,isize=size(fl2%field),ikind=kind(fl2%field))
     !
     allocate(fl1%field(nterms,0:Npoints),fl2%field(nterms,0:Npoints),stat=alloc)
     call ArrayStart(name1,alloc,size(fl1%field),kind(fl1%field))
     call ArrayStart(name2,alloc,size(fl2%field),kind(fl2%field))
     !
     fl1%field = Sfield1
     fl1%Ncoeff = Nterms
     fl1%iorder = Siorder1
     !
     fl2%field = Sfield2
     fl2%Ncoeff = Nterms
     fl2%iorder = Siorder2
     !
     deallocate(Sfield1,Sfield2,siorder1,siorder2)
     !
     call ArrayStop("Sfield")
     !
   end subroutine FLCompact_and_combine_fields_sparse
   !
   !
   subroutine FLCompact_and_combine_three_fields_sparse(fl1,name1,fl2,name2,fl3,name3)

     type(FLpolynomT),pointer  :: fl1,fl2,fl3
     character(len=*),intent(in) :: name1,name2,name3
     integer(ik)        :: Npoints,Ncoeff1,Ncoeff2,Ncoeff3,iterm,i,icoeff,Nterms,alloc,Nterm1,Nterm2,Ncoeffmax
     real(ark),allocatable    :: sfield1(:,:),sfield2(:,:),sfield3(:,:)    ! Expansion parameters in the sparse representation
     integer(ik),allocatable  :: siorder1(:),siorder2(:),siorder3(:)       ! iorder in sparse
     integer(ik)   :: target_index(trove%Nmodes)
     logical :: check = .true.
     !
     Ncoeff1 = fl1%Ncoeff
     Ncoeff2 = fl2%Ncoeff
     Ncoeff3 = fl3%Ncoeff
     !
     Npoints = fl1%Npoints
     !
     if (Npoints/=fl2%Npoints.or.Npoints/=fl3%Npoints) then
       write(out,"('FLCompact_and_combine_three_fields_sparse: Illegal Npoints in 3 fields, should be the same',3i8)") & 
             fl1%Npoints,fl2%Npoints,fl3%Npoints
       stop 'FLCompact_and_combine_three_fields_sparse: Illegal Npoints in 3 fields'
     endif
     !    
     target_index = 0
     target_index(1) = trove%NKinOrder 
     Ncoeffmax= FLQindex(trove%Nmodes_e,target_index)
     !
     ! Count large elements (using exp_coeff_thresh as threshold) and store in a sparse representation
     !
     iterm = 0
     !
     do icoeff = 1,Ncoeffmax
       if (any(abs(fl1%field(icoeff,:))>job%exp_coeff_thresh).or.&
           any(abs(fl2%field(icoeff,:))>job%exp_coeff_thresh).or.&
           any(abs(fl3%field(icoeff,:))>job%exp_coeff_thresh)) then
          iterm = iterm + 1
       endif
     enddo
     !
     nterms = iterm
     !
     if (Nterms==0) then
       !
       stop 'FLCompact_and_combine_three_fields_sparse is not implemented for Nterms=0' 
       !
       deallocate(fl1%IndexQ,fl2%IndexQ,fl3%IndexQ)
       call ArrayStop(name1//'IndexQ')
       call ArrayStop(name2//'IndexQ')
       call ArrayStop(name3//'IndexQ')
       !
       deallocate(fl1%ifromsparse,fl2%ifromsparse,fl3%ifromsparse)
       call ArrayStop(name1//"ifromsparse")
       call ArrayStop(name2//"ifromsparse")
       call ArrayStop(name3//"ifromsparse")
       !
       return
       !
     endif 
     !
     allocate(Sfield1(nterms,0:Npoints),Sfield2(nterms,0:Npoints),Sfield3(nterms,0:Npoints),stat=alloc)
     call ArrayStart("Sfield",alloc,size(Sfield1),kind(Sfield1))
     call ArrayStart("Sfield",alloc,size(Sfield2),kind(Sfield2))
     call ArrayStart("Sfield",alloc,size(Sfield3),kind(Sfield3))
     !
     allocate(Siorder1(nterms),Siorder2(nterms),Siorder3(nterms),stat=alloc)
     call ArrayStart("Sfield",alloc,size(Siorder1),kind(Siorder1))
     call ArrayStart("Sfield",alloc,size(Siorder2),kind(Siorder2))
     call ArrayStart("Sfield",alloc,size(Siorder3),kind(Siorder3))
     !
     deallocate(fl1%IndexQ,fl2%IndexQ,fl3%IndexQ)
     call ArrayStop(name1//'IndexQ')
     call ArrayStop(name2//'IndexQ')
     call ArrayStop(name3//'IndexQ')
     !
     deallocate(fl1%ifromsparse,fl2%ifromsparse,fl3%ifromsparse)
     call ArrayStop(name1//"ifromsparse")
     call ArrayStop(name2//"ifromsparse")
     call ArrayStop(name3//"ifromsparse")
     !
     allocate(fl1%ifromsparse(nterms),fl1%IndexQ(trove%Nmodes,nterms),stat=alloc)
     allocate(fl2%ifromsparse(nterms),fl2%IndexQ(trove%Nmodes,nterms),stat=alloc)
     allocate(fl3%ifromsparse(nterms),fl3%IndexQ(trove%Nmodes,nterms),stat=alloc)
     call ArrayStart(name1//"ifromsparse",alloc,size(fl1%ifromsparse),kind(fl1%ifromsparse))
     call ArrayStart(name1//"IndexQ",alloc,size(fl1%IndexQ),kind(fl1%IndexQ))
     call ArrayStart(name2//"ifromsparse",alloc,size(fl2%ifromsparse),kind(fl2%ifromsparse))
     call ArrayStart(name2//"IndexQ",alloc,size(fl2%IndexQ),kind(fl2%IndexQ))
     call ArrayStart(name3//"ifromsparse",alloc,size(fl3%ifromsparse),kind(fl3%ifromsparse))
     call ArrayStart(name3//"IndexQ",alloc,size(fl3%IndexQ),kind(fl3%IndexQ))
     !
     iterm = 0
     !
     do icoeff = 1,Ncoeffmax
       if (any(abs(fl1%field(icoeff,:))>job%exp_coeff_thresh).or.&
           any(abs(fl2%field(icoeff,:))>job%exp_coeff_thresh).or.&
           any(abs(fl3%field(icoeff,:))>job%exp_coeff_thresh)) then
          !
          iterm = iterm + 1
          !
          Sfield1(iterm,:) = fl1%field(icoeff,:)
          Sfield2(iterm,:) = fl2%field(icoeff,:)
          Sfield3(iterm,:) = fl3%field(icoeff,:)
          !
          siorder1(iterm) = fl1%iorder(icoeff)
          siorder2(iterm) = fl2%iorder(icoeff)
          siorder3(iterm) = fl3%iorder(icoeff)
          !
          fl1%ifromsparse(iterm) = icoeff
          fl1%IndexQ(:,iterm) = FLIndexQ(:,icoeff)
          !
          fl2%ifromsparse(iterm) = icoeff
          fl2%IndexQ(:,iterm) = FLIndexQ(:,icoeff)
          !
          fl3%ifromsparse(iterm) = icoeff
          fl3%IndexQ(:,iterm) = FLIndexQ(:,icoeff)
          !
       endif
     enddo
     !
     deallocate(fl1%iorder,fl2%iorder,fl3%iorder)
     call ArrayStop(name1)
     call ArrayStop(name2)
     call ArrayStop(name3)
     !
     ! Create a field in a sparse representaion
     !
     allocate(fl1%iorder(nterms),fl2%iorder(nterms),fl3%iorder(nterms),stat=alloc)
     call ArrayStart(name1,alloc,size(fl1%iorder),kind(fl1%iorder))
     call ArrayStart(name2,alloc,size(fl2%iorder),kind(fl2%iorder))
     call ArrayStart(name3,alloc,size(fl3%iorder),kind(fl3%iorder))
     !
     deallocate(fl1%field,fl2%field,fl3%field,stat=alloc)
     call ArrayMinus(name1,isize=size(fl1%field),ikind=kind(fl1%field))
     call ArrayMinus(name2,isize=size(fl2%field),ikind=kind(fl2%field))
     call ArrayMinus(name3,isize=size(fl3%field),ikind=kind(fl3%field))
     !
     allocate(fl1%field(nterms,0:Npoints),fl2%field(nterms,0:Npoints),fl3%field(nterms,0:Npoints),stat=alloc)
     call ArrayStart(name1,alloc,size(fl1%field),kind(fl1%field))
     call ArrayStart(name2,alloc,size(fl2%field),kind(fl2%field))
     call ArrayStart(name3,alloc,size(fl3%field),kind(fl3%field))
     !
     fl1%field = Sfield1
     fl1%Ncoeff = Nterms
     fl1%iorder = Siorder1
     !
     fl2%field = Sfield2
     fl2%Ncoeff = Nterms
     fl2%iorder = Siorder2
     !
     fl3%field = Sfield3
     fl3%Ncoeff = Nterms
     fl3%iorder = Siorder3
     !
     deallocate(Sfield1,siorder1,Sfield2,siorder2,Sfield3,siorder3)
     !
     call ArrayStop("Sfield")
     !
   end subroutine FLCompact_and_combine_three_fields_sparse


   !
   subroutine FLfingerprint(action,chkptIO,PTorder,Npolyads,enercut)
    !
    character(len=*), intent(in) :: action ! 'read'  or 'write'
    integer(ik),intent(in) :: chkptIO,PTorder,Npolyads
    real(rk),intent(in   ) :: enercut(2)
     !
     if (job%verbose>=7) write(out,"(/'FLfingerprint/start')") 
     !
     select case (action)
       case default
         write (out,"(' FLfingerprint - action ',a,' is not valid')") trim(action)
         stop 'FLfingerprint - bogus command'
       case ('READ','read')
         call fingerprintRead
       case ('WRITE','write')
         call fingerprintWrite
     end select
     !
     if (job%verbose>=7) write(out,"('FLfingerprint/stop')") 
     !
   contains 
     !
     subroutine fingerprintWrite
      !
      integer(ik)  :: imode
      character(len=10) :: ifmt_modes,rfmt_atoms,rfmt_modes,rfmt_coords
      !
      write(chkptIO,"('Start Fingerprints')") 
      !
      write(ifmt_modes,"(i3,'i4')") trove%Nmodes+1
      ifmt_modes = adjustl(ifmt_modes)
      write(rfmt_modes,"(i3,'f18.8')") trove%Nmodes
      rfmt_modes = adjustl(rfmt_modes)
      write(rfmt_coords,"(i3,'f18.8')") trove%Ncoords
      rfmt_coords = adjustl(rfmt_coords)
      write(rfmt_atoms,"(i3,'f18.8')") trove%Natoms
      rfmt_atoms = adjustl(rfmt_atoms)
      !
      write(chkptIO,"(4i7,2f12.1,'  <= PTorder, Nmodes, Natoms,  Npolyads, enercut ')") &
                       PTorder,trove%Nmodes,trove%Natoms,Npolyads,enercut(:)
      write(chkptIO,'('//rfmt_modes//',''  <= dstep       '')')  trove%fdstep(:)
      write(chkptIO,'('//rfmt_atoms//',''  <= masses      '')')  trove%mass
      write(chkptIO,'('//rfmt_coords//',''  <= equilbrium  '')')  trove%local_eq
      !
      write(chkptIO,"(a10,' ',a10,' ',a10,' ',a10)") trove%internal_coords,trove%coords_transform,trove%symmetry,trove%Moltype
      !
      write(chkptIO,"('BASIS:   i  type     coord_kinet coord_poten model dim species class range dvrpoints',1x,&
                       &'res_coeffs npoints borders periodic period')") 
      !
      write(chkptIO,"(i8,'   <- Jrot, rotational angular momentum')") bset%dscr(0)%range(1)
      !
      do imode = 0,trove%Nmodes
        write(chkptIO,"(6x,i4,1x,3(a10,1x),i5,3x,a2,3x,i2,5x,i2,1x,2i4,2x,f6.1,2x,i9,1x,2f9.3,1x,i2,1x,i2,1x,a10,i9,i3,i3,i3)") &
                      imode, bset%dscr(imode)
      enddo
      !
      write(chkptIO,"('End Fingerprints')") 
      !
    end subroutine fingerprintWrite
    !  
    subroutine fingerprintRead
      !
      integer(ik)  :: imode,PTorder_t,Nmodes_t,Natoms_t,Npolyads_t,jrot_t,imode_
      real(ark)     :: enercut_t(1:2),f_t(1:trove%Nmodes), mass_t(1:trove%Natoms),g_t(1:trove%Ncoords)
      character(len=18) :: buf
      character(len=43) :: buf43
      character(len=10) :: char_t
      type(FLbasissetT) :: bs_
      !
      read(chkptIO,"(a18)") buf
      if (buf/='Start Fingerprints') then
        write (out,"(' fingerprintRead file has bogus header: ',a)") buf
        stop 'fingerprintRead - bogus file header'
      end if
      !
      read(chkptIO,*)  PTorder_t,Nmodes_t,Natoms_t,Npolyads_t,enercut_t(1:2)
      !
      if (PTorder_t/=PTorder.or.Nmodes_t/=trove%Nmodes.or. & 
          Natoms_t/=trove%Natoms.or.Npolyads_t/=Npolyads) then
        write (out,"(' fingerprintRead:  parameters mismatch: ')") 
        write (out,"('PTorder  : ',2i8)") PTorder_t,PTorder
        write (out,"('Nmodes   : ',2i8)") Nmodes_t,trove%Nmodes
        write (out,"('Natoms   : ',2i8)") Natoms_t,trove%Natoms
        write (out,"('Npolyads : ',2i8)") Npolyads_t,Npolyads
        stop 'fingerprintRead - parameters mismatch'
      end if
      !
      if (enercut_t(1)/=enercut(1).or.enercut_t(2)<enercut(2)) then
        write (out,"(' fingerprintRead:  enercut mismatch: ',4f18.4)") enercut_t(:),enercut(:)
        stop 'fingerprintRead - enercut mismatch'
      end if
      !
      read(chkptIO,*)  f_t(1:Nmodes_t)
      !
      if (any(abs(f_t(:)-trove%fdstep(:))>sqrt(small_))) then
        write (out,"(' fingerprintRead:  fdstep mismatch: ')")
        write (out,"(' fingerprintRead:  actual: ',40f18.8)") trove%fdstep(:)
        write (out,"(' fingerprintRead:  stored: ',40f18.8)") f_t(:)
        stop 'fingerprintRead - enercut mismatch'
      end if
      !
      read(chkptIO,*)  mass_t(1:Natoms_t)
      read(chkptIO,*)  g_t(1:trove%Ncoords)
      !
      if (any(abs(mass_t(:)-trove%mass(:))>1e-7)) then
        write (out,"(' fingerprintRead:  mass mismatch: ')")
        write (out,"(' fingerprintRead:  actual: ',40f18.4)") trove%mass(:)
        write (out,"(' fingerprintRead:  stored: ',40f18.4)") mass_t(:)
        stop 'fingerprintRead - mass mismatch'
      end if
      !
      if (any(abs(g_t(:)-trove%local_eq(:))>1e-7)) then
        write (out,"(' fingerprintRead:  local_eq mismatch: ')")
        write (out,"(' fingerprintRead:  actual: ',40f20.10)") trove%local_eq(:)
        write (out,"(' fingerprintRead:  stored: ',40f20.10)") g_t(:)
        stop 'fingerprintRead - local_eq mismatch'
      end if
      !
      read(chkptIO,"(a43)") buf43
      write(char_t,"(a10)") trove%internal_coords 
      !
      if (buf43(1:10)/=char_t) then
        write (out,"(' fingerprintRead:  internal_coords  mismatch: ',2a10)") buf43(1:10),trove%internal_coords
        stop 'fingerprintRead - internal_coords mismatch'
      end if
      !
      write(char_t,"(a10)") trove%coords_transform 
      if (buf43(12:21)/=char_t) then
        write (out,"(' fingerprintRead:  coords_transform  mismatch: ',2a10)") buf43(12:21),trove%coords_transform
        stop 'fingerprintRead - coords_transform mismatch'
      end if
      !
      write(char_t,"(a10)") trove%symmetry 
      if (buf43(23:32)/=char_t) then
        write (out,"(' fingerprintRead:  symmetry  mismatch: ',2a10)") buf43(23:32),trove%symmetry
        stop 'fingerprintRead - symmetry mismatch'
      end if
      !
      read(chkptIO,"(1x,a10)") char_t
      !
      read(chkptIO,*) jrot_t
      !
      !if (jrot_t/=bset%dscr(0)%range(1).and.job%verbose>=3) then
      !  write (out,"(' fingerprintRead:  Jrot mismatch ')")
      !  write (out,"(' fingerprintRead:  actual: ',i8)") bset%dscr(0)%range(1)
      !  write (out,"(' fingerprintRead:  stored: ',i8)") jrot_t
      !end if
      !
      do imode = 0,trove%Nmodes
        !
        !read(chkptIO,"(6x,i4,1x,3(a10,1x),i5,3x,a2,3x,i2,5x,i2,1x,3i4,2x,f6.1,2x,i9,1x,2f9.3,1x,i,i)") imode_,bs_
        !
        read(chkptIO,"(6x,i4,1x,3(a10,1x),i5,3x,a2,3x,i2,5x,i2,1x,2i4,2x,f6.1,2x,i9,1x,2f9.3,1x,i2,1x,i2,1x,a10,i9,i2,i2)") & 
        imode_,bs_%type,bs_%COORD_KINET,bs_%COORD_POTEN,bs_%MODEL,bs_%DIM,bs_%SPECIES,bs_%CLASS,bs_%RANGE,&
               bs_%RES_COEFFS,bs_%NPOINTS,bs_%BORDERS,bs_%PERIODIC,bs_%IPERIOD
        !
        if (bs_%range(2)/=job%bset(imode)%range(2)) then
          write (out,"('fingerprintRead:  parameters mismatch for  ',i9,'th mode:')") imode
          write (out,"('range2 (stored) /=  range (given)  : ',2i8)") bs_%range(2),job%bset(imode)%range(2)
          stop 'fingerprintRead - parameters mismatch:range'
        end if
        !
        if (bs_%IPERIOD/=job%bset(imode)%IPERIOD) then
          write (out,"('fingerprintRead:  parameters mismatch for  ',i9,'th mode:')") imode
          write (out,"('IPERIOD (stored) /=  IPERIOD (given)  : ',2i8)") bs_%IPERIOD,job%bset(imode)%IPERIOD
          stop 'fingerprintRead - parameters mismatch:IPERIOD'
        end if
        !
        !read(chkptIO,"(a10)") char_t
        !        
		!bs_%RES_COEFFS,bs_%NPOINTS,bs_%BORDERS,bs_%PERIODIC,bs_%IPERIOD,
		!bs_%DVR
		!bs_%DVRPOINTS
		!bs_%POSTPROCESS
		!bs_%LVIB
        !
      enddo
      !
      read(chkptIO,"(a16)") buf(1:16)
      if (buf(1:16)/='End Fingerprints') then
        write (out,"(' fingerprintRead file has bogus footer: ',a)") buf
        stop 'fingerprintRead - bogus file footer'
      end if
      !
    end subroutine fingerprintRead
    !
   end subroutine FLfingerprint


 !read eigenvalues and their assignment 
 !
 subroutine FLread_ZPE
    !
    implicit none
    !
    integer(ik)             :: nmodes, nroots, iroot, igamma,iounit, info, nroots_t, Npolyad_t,nsize,irec,ilevel,ideg
    real(rk)                :: energy
    character(cl)           :: filename, ioname, buf
    character(4)            :: jchar,gchar
    character(500)          :: buf500
    integer(ik)             :: gamma
    real(rk)                :: energy_
    !
    !type(PTeigenT)          :: eigen_t   ! temporal object used for sorting 'eigen'
    !
    if (job%ZPE>0) return
    ! 
    nmodes = FLNmodes
    !
    gamma = 1
    !
    write(jchar, '(i4)') 0
    write(gchar, '(i2)') gamma
    !
    filename = trim(job%eigenfile%filebase)//'_descr'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
    !
    write(ioname, '(a, i4,2x,i2)') 'eigenvalues for J,gamma = ', 0,gamma
    !
    call IOstart(trim(ioname), iounit)
    open(unit = iounit, action = 'read',status='old' , file = filename,err=14)
    !
    ! Check the fingerprint of the computed eigenvectors. 
    !
    call FLfingerprint('read',iounit,0,job%Npolyads_prim,(/job%enercutoff%primt,job%enercutoff%contr/))
    !
    buf500 = ''
    !
    read(iounit, '(a)') buf500
    read(iounit,"(i8,a4)") Npolyad_t,buf500(1:4)
    read(iounit, '(a)') buf500
    !
    ! Start reading the description of the eigensolution. 
    !
    read(iounit,*) nroots_t,nsize
    !
    read(iounit, '(a)') buf500
    if (buf500(1:3) == 'End') then 
      stop 'FLread_ZPE: no energies in J=0,gamma=1 file'
    endif 
    !
    read(buf500, *) irec, igamma, ilevel, ideg, energy
    !
    if (igamma/=gamma) then 
      write(out,"('FLread_ZPE error: igamma/=gamma: ',2i4)") igamma,gamma
      stop 'FLread_ZPE error: igamma/=gamma'
    endif
    !
    job%zpe = energy
    job%partfunc%zpe = energy
    !
    return
    !
14  if (job%verbose>=3) &
        write(out,"('Warning: the egine_descr0_1.chk files does nto exist to define the ZPE value; default will be used')")
    !
 end subroutine FLread_ZPE





!
! Assign PT-orders for different kinetic and potential expansion terms
!
  subroutine FLpt_orders_distribution

  integer(ik)                ::   iterm,imode,jmode
    !
    if (job%verbose>=4) write(out,"(/'FLpt_orders_distribution/start')") 
    !
    !
    if (job%verbose>=3) write(out,"('Distribution of polinomial coeff. over PT-orders')") 
    !
    ! Here we define the distribution of the field-coeffs with respect to different orders
    !
    if (job%verbose>=3) write(out,"('poten:'/'             #    PT  modes ')") 
    !
    call PT_orders_distribution(job%pot_pt_shift,trove%poten)
    !
    if (job%verbose>=3) then 
       !
       write(out,"('       #  pot     powers')") 
       !
       do iterm = 1,trove%poten%Ncoeff
          !
          write(out,"(6x,i8,i5,30i4)") iterm,trove%poten%iorder(iterm), & 
                              FLIndexQ(1:min(30,trove%Nmodes),iterm)
          !
       enddo
       !
    endif 
    !
    if (job%verbose>=3) write(out,"('g_vib:'/'  k1 k2      #    PT  modes ')") 
    !
    do imode =1,trove%Nmodes  
       do jmode =1,trove%Nmodes  
          !
          call PT_orders_distribution(0,trove%g_vib(imode,jmode),(/imode,jmode/))
          !
          if (job%verbose>=3) then 
             !
             do iterm = 1,trove%g_vib(1,1)%Ncoeff
               write(out,"(2i3,i8,i5,30i4)") imode,jmode,iterm,trove%g_vib (imode,jmode)%iorder(iterm), & 
                                             FLIndexQ(1:min(30,trove%Nmodes),iterm)
             enddo
             !
          endif 
          !
       enddo
    enddo
    !
    if (job%verbose>=3) write(out,"('g_rot:'/'  k1 k2      #    PT  modes ')") 
    !
    do imode =1,3
       do jmode =1,3
          call PT_orders_distribution(0,trove%g_rot(imode,jmode),(/-imode,-jmode/))
          !
          if (job%verbose>=3) then 
             !
             do iterm = 1,trove%g_rot(1,1)%Ncoeff
               write(out,"(2i3,i8,i5,30i4)") imode,jmode,iterm,trove%g_rot (imode,jmode)%iorder(iterm), & 
                                             FLIndexQ(1:min(30,trove%Nmodes),iterm)
             enddo
             !
          endif 
          !
       enddo
    enddo
    !
    if (job%verbose>=3) write(out,"('g_cor:'/'  k1 k2      #    PT  modes ')") 
    !
    do imode =1,trove%Nmodes  
       do jmode =1,3
          !
          call PT_orders_distribution(0,trove%g_cor(imode,jmode),(/imode,-jmode/))
          !
          if (job%verbose>=3) then 
             !
             do iterm = 1,trove%g_cor(1,1)%Ncoeff
               write(out,"(2i3,i8,i5,30i4)") imode,jmode,iterm,trove%g_cor (imode,jmode)%iorder(iterm), & 
                                             FLIndexQ(1:min(30,trove%Nmodes),iterm)
             enddo
             !
          endif 
          !
       enddo
    enddo
    !
    if (FLl2_coeffs) then
      !
      if (job%verbose>=3) write(out,"('L2_vib:'/'  k1 k2      #    PT  modes ')") 
      !
      do imode =1,trove%Nmodes  
         do jmode =1,trove%Nmodes  
            !
            call PT_orders_distribution(0,trove%L2_vib(imode,jmode),(/imode,jmode/))
            !
            if (job%verbose>=3) then 
               !
               do iterm = 1,trove%L2_vib(1,1)%Ncoeff
                 write(out,"(2i3,i8,i5,30i4)") imode,jmode,iterm,trove%L2_vib (imode,jmode)%iorder(iterm), & 
                                               FLIndexQ(1:min(30,trove%Nmodes),iterm)
               enddo
               !
            endif 
            !
         enddo
      enddo
      !
    endif
    !
    if (job%verbose>=4) write(out,"('FLpt_orders_distribution/end')") 
    !
  contains 

    !
    ! Here we define the distribution of the field-coeffs with respect to different orders
    !
    subroutine PT_orders_distribution(ordershift,fl,gk)
    !
    integer(ik),intent(in)  :: ordershift
    type(FLpolynomT),intent(inout)  :: fl
    integer(ik),optional       :: gk(2)
    !
    integer(ik) :: ki(trove%Nmodes),power_k,kdiag,iterm,iorder_t,imode
    character(len=cl)            :: char_t
      !
      do iterm = 1,fl%Ncoeff
         !
         ki(:) = FLIndexQ(:,iterm)
         !
         power_k = sum(ki)
         ! 
         !
         ! if kdiag<=1 - the term will be diagonal
         ! lets check it! 
         !
         kdiag = 0 
         !
         !
         ! if gk - present, we also take inot account the momenta_i 
         ! when checking for diagonality
         !
         if (present(gk)) then 
            !
            if (gk(1)/=gk(2)) then 
               kdiag = kdiag+2 
            endif 
            !
         endif 
         !
         do imode = 1,trove%Nmodes
            !
            if (ki(imode)/=0) then 
               kdiag = kdiag+1 
               !
               if (present(gk)) then 
                  !
                  if (any( imode/=gk(:) )) then 
                     kdiag = kdiag+1 
                  endif 
                  !
               endif 
               !
            endif 
            !
            char_t = trim(bset%dscr(imode)%type)
            !
            if (bset%dscr(imode)%type=='NUMEROV') then 
              !
              power_k = power_k-ki(imode) 
              !
              continue
              !
            endif 
            !
         enddo 
         !
         if (power_k<=ordershift.and.kdiag<=1) then 
            iorder_t = manifold
         elseif(power_k<=ordershift) then 
            iorder_t = 1 
         else
            iorder_t = power_k-ordershift
         endif 
         !
         fl%iorder(iterm) = iorder_t
         !
      enddo 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! fl%iorder = 0 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine PT_orders_distribution


  end subroutine FLpt_orders_distribution


!
! Add new 1D basis set type 
! 
  subroutine FLbset1DNew(ibs,BSsize)

    integer(ik),intent(in)      :: ibs         ! Index for the new 1D basis   
    integer(ik),intent(inout)   :: BSsize       ! Size of the 1D basis set 

    integer(ik)                 :: MatrixSize,imode,k,ipower,iterm,Nmodes,Tcoeff,ialloc,irho_eq,icoeff,jmode
    integer(ik)                 :: imu,alloc,alloc_p,nu_i,powers(trove%Nmodes),npoints,vl,vr,k1,k2,i,i_,isingular,jrot,krot,&
                                   kmax,nmax,krot1,krot2,krot11,krot21,k_l,k_r,i1,i2,j
    integer(ik)                 :: nl,nr,irho
    type(FLpolynomT),pointer    :: fl,gl
    type(Basis1DT), pointer     :: bs           ! 1D bset
    real(ark)                    :: f2(1:trove%Nmodes),g2(1:trove%Nmodes),f_t,g_t,rmk,amorse,f_m
    real(ark)                    :: f1d(0:trove%MaxOrder),p1d(0:trove%MaxOrder),g1d(0:trove%MaxOrder),chi(trove%Nmodes)
    real(ark)                   :: rho,rho_pot,rho_kin,rho_kin0,rho_pot0,rho_ext,rho_ext0
    real(ark)    :: ar_t(1:molec%ncoords)
    character(len=cl)     :: dir
    !
    real(ark),allocatable        :: f1drho(:),g1drho(:),drho(:,:),muzz(:),sinrho(:),cosrho(:),pseudo(:),mrho(:),xton(:,:)
    !
    integer(ik),parameter        ::  Nperiod_max = 10 
    real(ark)                    :: f_period(1:trove%Nmodes,Nperiod_max)
    !
    real(ark),allocatable       :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),phivphi_t(:),weight(:)
    real(ark),allocatable       :: phil_leg(:),phir_leg(:),dphil_leg(:),dphir_leg(:),phil_sin(:),phir_sin(:)
    real(ark),allocatable       :: dfunc(:,:),func(:,:)
    !
    real(ark)                   :: rho_b(2),step,rho_ref,mat_t,sqrt2,L
    real(ark)                   :: rho_range,rho_t
    integer(ik)                 :: io_slot       ! unit numeber to store the numerov eigenvectors and their derivatives
    integer(ik)                 :: iperiod=0,rec_len,iparity,numerpoints
    character(len=cl)    :: unitfname,char_
    !
    logical              :: reduced_model,periodic_model 
    real(ark)            :: poten_t,gvib_t(trove%Nmodes,trove%Nmodes),grot_t(3,3),gcor_t(trove%Nmodes,3),extF_t(extF%rank)
    real(ark)            :: rho_ref_,period
    !
    real(ark)   ::  rho_switch  = .0174532925199432957692369_ark       ! the value of abcisse rho of the switch between regions (1 deg)
    integer(ik) ::  iswitch                                 ! the grid point of switch
    real(ark)   :: g2_term
    !
    ! substitute for easier reference 
    !
    bs => bset%bs1D(ibs)

    if (job%verbose>=4) write(out,"(/'FLbset1DNew/start')") 
    !
    !  Initialize the basis set 
    !
    MatrixSize = 3*(trove%MaxOrder+1)*(BSsize+1)*(BSsize+1)
    bs%Size = BSsize
    bs%Order = trove%MaxOrder
    npoints = trove%Npoints
    rho_b  = trove%rho_border
    Nmodes = trove%Nmodes
    rho_range = rho_b(2)-rho_b(1)
    isingular = -1    !
    periodic_model = .false.
    period = 0
    !
    if (job%verbose>=6) then
      write (out,"(' Basis type ',a,' needs ',f12.5,' Mbytes of memory (plus a bit)')") &
             trim(bs%name), real(rk*MatrixSize,kind=rk)/(1024.0_rk**2)
    end if
    !
    allocate (bs%matelements(-1:3,0:trove%MaxOrder,0:BSsize,0:BSsize),bs%ener0(0:BSsize),stat=alloc)
    call ArrayStart('bs%matelements',alloc,1_ik,kind(bs%matelements),size(bs%matelements,kind=hik))
    call ArrayStart('bs%ener0',alloc,size(bs%ener0),kind(bs%ener0))
    !
    ! at the last mode we allocate the matrix elements arrays:
    !
    if (any(bs%mode(:)==Nmodes)) then 
        !
        ! These fields will be used for storing the Numerov-matrix elements
        !
        ialloc = 0
        !
        do k1 = 1,Nmodes
           do k2 = 1,Nmodes
              !
              fl => trove%g_vib(k1,k2)
              !
              Tcoeff = fl%Ncoeff
              !
              allocate (fl%me(Tcoeff,0:bs%Size,0:bs%Size),stat=alloc)
              call ArrayStart('trove%g_vib%me',alloc,size(fl%me),kind(fl%me))
              !
              ialloc = ialloc+abs(alloc)
              !
              fl%me  = 0
              !
           enddo
        enddo
        !
        do k1 = 1,3
           do k2 = 1,3
              !
              fl => trove%g_rot(k1,k2)
              !
              Tcoeff = fl%Ncoeff
              !
              allocate (fl%me(Tcoeff,0:bs%Size,0:bs%Size),stat=alloc)
              call ArrayStart('trove%g_rot%me',alloc,size(fl%me),kind(fl%me))
              !
              ialloc = ialloc+abs(alloc)
              !
              fl%me  = 0
              !
           enddo
        enddo
        !
        do k1 = 1,Nmodes
           do k2 = 1,3
              !
              fl => trove%g_cor(k1,k2)
              !
              Tcoeff = fl%Ncoeff
              !
              allocate (fl%me(Tcoeff,0:bs%Size,0:bs%Size),stat=alloc)
              call ArrayStart('trove%g_cor%me',alloc,size(fl%me),kind(fl%me))
              !
              ialloc = ialloc+abs(alloc)
              !
              fl%me  = 0
              !
           enddo
        enddo
        !
        Tcoeff = trove%poten%Ncoeff
        !
        allocate (trove%poten%me(Tcoeff,0:bs%Size,0:bs%Size),stat=alloc)
        call ArrayStart('trove%poten%me',alloc,size(fl%me),kind(fl%me))
        !
        trove%poten%me = 0
        !
        ialloc = ialloc+abs(alloc)
        !
        if (ialloc/=0) then
           write (out,"(' Error ',i9,' trying to allocate ME-arrays')") ialloc
           stop 'FLbset1DNew, me - out of memory'
        end if
        !
        if (FLextF_coeffs) then
          !
          do imu = 1,extF%rank
             !
             fl => trove%extF(imu)
             !
             Tcoeff = fl%Ncoeff
             !
             allocate (fl%me(Tcoeff,0:bs%Size,0:bs%Size),stat=alloc)
             call ArrayStart('trove%extF%me',alloc,size(fl%me),kind(fl%me))
             !
             ialloc = ialloc+abs(alloc)
             !
             fl%me  = 0
             !
          enddo
          !
        endif
        !
        if (FLl2_coeffs) then
           !
           do k1 = 1,Nmodes
              do k2 = 1,Nmodes
                 !
                 fl => trove%L2_vib(k1,k2)
                 !
                 Tcoeff = fl%Ncoeff
                 !
                 if (Tcoeff<1) cycle
                 !
                 allocate (fl%me(Tcoeff,0:bs%Size,0:bs%Size),stat=alloc)
                 call ArrayStart('trove%L2_vib%me',alloc,size(fl%me),kind(fl%me))
                 !
                 ialloc = ialloc+abs(alloc)
                 !
                 fl%me  = 0
                 !
              enddo
           enddo
           !
        endif
        !
    endif 
    !
    ! define the quadratic potential and kinetic constants f2 and g2
    !
    fl => trove%poten
    !
    if (size(fl%field,dim=1)/=fl%Ncoeff ) then
      write (out,"(' FLbset1DNew: fields poten not defined yet and cannot be used ')") 
      stop 'FLbset1DNew: fields poten not defined'
    end if
    !
    ! define the conversion factor to the normal coordinate
    !
    do imode = 1,bs%imodes
       !
       nu_i = bs%mode(imode)
       !
       if (trim(bset%dscr(nu_i)%dvr)=='HERMITE') then
         !
         irho_eq = 0 
         !
         if (manifold/=0) irho_eq = mod(nint( ( molec%chi_eq(trove%Nmodes)-trove%rho_border(1) )/&
                                                (trove%rhostep),kind=ik ),trove%npoints)
         !
         powers = 0 ; powers(nu_i) = 2
         k = FLQindex(trove%Nmodes_e,powers)
         !
         if (trove%sparse) then
           !
           call find_isparse_from_ifull(trove%poten%Ncoeff,trove%poten%ifromsparse,k,i)
           !
           f2(1) = 0 
           !
           if (i/=0) f2(1) = trove%poten%field(i,irho_eq)
           !
           g2(1) = trove%g_vib(nu_i,nu_i)%field(1,irho_eq)
           !
           if (abs(g2(1))<sqrt(small_).or.trove%g_vib(nu_i,nu_i)%ifromsparse(1)/=1) then 
             write(out,"('FLbset1DNew: g(2)=0 in the sparse-field or inconsistent sparse-recored/=1',i8)") &
                   trove%g_vib(nu_i,nu_i)%ifromsparse(1)
             stop 'FLbset1DNew: g(2)=0 in the sparse-field'
           endif
           !
         else
           !
           f2(1) = trove%poten%field(k,irho_eq)
           g2(1) = trove%g_vib(nu_i,nu_i)%field(1,irho_eq)
           !
         endif
         !
         trove%coord_f(nu_i) = sqrt( sqrt( g2(1)/( 2.0_ark*f2(1) ) ) )
         !
       endif
       !
    enddo 
    !
    ! Here we generate the 1D basis set matrix elements 
    ! the type of the basis set will define the type of the matrix elements generator 
    !
    !
    select case (trim(bs%type)) ! -----------------------------------------------------------------------
      case default
        write (out,"('FLbset1DNew: basis set ',a,' unknown')") trim(bs%name)
        stop 'FLbset1DNew: unknown basis set '
         !
      case('NORMAL','MORSE','HARMONIC')    ! -----------------------------------------------------------------------
        ! 
        ! harmonic  or morse bsets
        !
        ! for the harmonic-morse basis sets we will need the quadratic potential and zero-order kinetic parameters
        !
        select case(manifold)
        !
        case (0) ! rigid case, expansion around the rank=0 manifold
          !
          irho_eq = 0 
          !
        case (1) ! non-rigid case, expansion around a rank=1 manifold
          !
          irho_eq = mod(nint( ( molec%chi_eq(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
          !
        end select
        !
        if (irho_eq<0.or.irho_eq>trove%Npoints) then
           write(out,"('FLbset1DNew-numerov: equilibrium is ill defined:')")
           write(out,"('irho_eq is ',i6)") irho_eq
           stop 'FLbset1DNew:  equilibrium is ill defined'
        endif
        !
        f2 = 0 
        do imode = 1,bs%imodes
           nu_i = bs%mode(imode)
           powers = 0 ; powers(nu_i) = 2
           k = FLQindex(trove%Nmodes_e,powers)
           !
           if (trove%sparse) then
             !
             call find_isparse_from_ifull(trove%poten%Ncoeff,trove%poten%ifromsparse,k,i)
             !
             f2(imode) = 0 
             !
             if (i/=0) f2(imode) = trove%poten%field(i,irho_eq)
             !
           else
             !
             f2(imode) = fl%field(k,irho_eq)
             !
           endif
           !
        enddo 
        !
        ! Check if all quadratic parameters are equal for every considered mode
        !
        f_t= f2(1)
        if (any(abs(f2(1:bs%imodes)-f_t)>1e-7*abs(f_t))) then
           write(out,"('FLbset1DNew: not all quadratic parameters are equal')")
           if (bs%imodes<30) then
             write(out,"(30f18.8)") f2(1:bs%imodes)
           endif 
           if ( job%bset(nu_i)%check_sym ) then
             stop 'FLbset1DNew: not all quadratic parameters are equal'
           endif
        endif
        !
        ! for the kinetic part it is simpler - we just take the first member of the g_vib%coeffs
        !
        g2 = 0
        do imode =1,bs%imodes  
          nu_i = bs%mode(imode)
          fl => trove%g_vib(nu_i,nu_i)
          g2(imode) = fl%field(1,irho_eq)
          !
          if (trove%sparse) then
            !
            g2(imode) = trove%g_vib(nu_i,nu_i)%field(1,irho_eq)
            !
            if (abs(g2(1))<sqrt(small_).or.trove%g_vib(nu_i,nu_i)%ifromsparse(1)/=1) then 
              write(out,"('FLbset1DNew: g(2)=0 in the sparse-field or inconsistent sparse-recored/=1',i8)") &
                    trove%g_vib(nu_i,nu_i)%ifromsparse(1)
              stop 'FLbset1DNew: g(2)=0 in the sparse-field'
            endif
            !
          endif
          !
        enddo 
        !
        ! The same check for all g_0 kinetic parameters to be equal for every considered mode
        !
        f_t= g2(1) 
        if (any(abs(g2(1:bs%imodes)-f_t)>10000.0_rk*sqrt(small_)*abs(f_t))) then
           write(out,"('FLbset1DNew: not all zero-order kinetic parameters are equal')")
           if (bs%imodes<30) then
             write(out,"(30f18.8)") g2(1:bs%imodes)
             write(out,"('difference somewhat greater than ',f18.8)") 10000.0_rk*sqrt(small_)*abs(f_t)
           endif 
           if ( job%bset(nu_i)%check_sym ) then 
             stop 'FLbset1DNew: not all zero-order kinetic parameters are equal'
           endif
        endif
        !
        ! Define the conditional basis set parameters   
        !
        ! Check if defined f2 and g2 parameters are positive
        !
        if (g2(1)<=0.0_rk ) then
           write(out,"('FLbset1DNew: g2 is not positive ',f18.8)") g2(1)
           stop 'FLbset1DNew: g2 is not positive'
        endif

        if (f2(1)<=0.0_rk ) then
           write(out,"('FLbset1DNew: f2 is not positive ',f18.8,' use specparam')") f2(1)
           f2(1) = trove%specparam(bs%mode(1))**2*0.5_ark/g2(1)
           !stop 'FLbset1DNew: f2 or g2 are not positive'
        endif


        ! Here we define the difference between harmonic, normal, and morse bsets
        if (trim(bs%type)=='HARMONIC'.or.trim(bs%type)=='NORMAL') then 

          !
          ! use specparam if given
          !
          if (trove%specparam(bs%mode(1))>0d0 ) then
             !
             f2(1) = trove%specparam(bs%mode(1))**2*0.5_ark/g2(1)
             !
             if (job%verbose>=4) then
               write(out,"('the input special-parameter ',f18.8,' will be used as f2 for the Harmonic basis set')") & 
                     trove%specparam(bs%mode(1)) 
             endif
             !
          endif
          !
          bs%params    = 0
          bs%params(1) = sqrt( sqrt( g2(1)/( 2.0_ark*f2(1) ) ) )  ! conversion parameter to the normal coordinates
          bs%params(2) =       sqrt( 2.0_ark*g2(1)*f2(1)  )       ! omega (harmonic parameter)
          !
          f_t=bs%params(1)
          !
          call ME_harmonic(bs%Size,bs%order,f_t,bs%matelements,bs%ener0) 
          !
          ! energies have been found for a normalized Morse oscilator, here we compute the actual energies
          !
          bs%ener0 = 0.5_ark*g2(1)*bs%ener0
          !
          if (job%verbose>=2) then
             !
             write (out,"(/'Zero-order Harmonic energies:')")
             !
             do nu_i =0,bs%Size 
                write (out,"(i8,f16.4)") nu_i,bs%ener0(nu_i)
             enddo
             !
          endif 
         !
        elseif  (trim(bs%type)=='MORSE') then 
          !
          amorse = trove%specparam(bs%mode(1))
          rmk  =sqrt(2.0_ark*f2(1)/g2(1))/amorse
          !
          ! Morse basis set cannot go above De (morse parameter), which 
          ! is limited by the parameter rmk (the maximal number of bound states).
          ! If rmk is smaller than bs%Size we stop here and suggest to reduce the basis set range 
          !
          if (rmk<bs%Size) then
             write(out,"('FLbset1DNew: we cannot fit ',i8,' morse basis states into PES with rmk = ',f16.8)") bs%Size,rmk
             write(out,"('FLbset1DNew: consider increasing bset%range')")
             stop 'FLbset1DNew: rmk is too small - increase bset%range'
          endif
          !
          ! Check just in case if the morse basis set used for right coordinates 
          !
          if (any(trove%Coordinates(:,1:bs%imodes)/='MORSE')) then
             write(out,"('FLbset1DNew: morse bset cannot be used for coord-s ',30a)") trove%Coordinates(1,1:min(bs%imodes,30))
             stop 'FLbset1DNew: coordinates type not consistent with bset'
          endif
          !
          bs%params    = 0
          bs%params(1) =  rmk
          bs%params(2) =  amorse
          !
          call ME_morse(bs%Size,bs%order,  rmk, amorse, bs%matelements,bs%ener0) 
          !
          ! energies have been found for a normalized Morse oscilatro, here we compute the actula energies
          !
          bs%ener0 = 0.5_ark*g2(1)*bs%ener0
          !
          if (job%verbose>=2) then
             !
             write (out,"(/'Zero-order Morse energies:')")
             !
             do nu_i =0,bs%Size 
                write (out,"(i8,f16.4)") nu_i,bs%ener0(nu_i)
             enddo
             !
          endif 
          !
        else 
           !
           write(out,"('FLbset1DNew: only harmonic or morse bsets are assumed, not ',a)") trim(bs%type)
           stop 'FLbset1DNew: bad bset-type'
           !
        endif 
        !
     case('LAGUERRE') 
        ! 
        ! numerov bset
        if (trove%manifold_rank(bs%mode(1))/=0) then
           !
           ! Allocation of the potential and kinetic 1d matrixes
           !
           if (job%bset(nu_i)%iperiod>0.and.trim(bs%type)=='NUMEROV') then
             iperiod = job%bset(nu_i)%iperiod
             rho_b(2) = rho_b(2)/real(iperiod,ark)
             npoints = npoints/iperiod
           endif
           !
           if (job%bset(nu_i)%iperiod==-2.and.trim(bs%type)=='NUMEROV') then
             iperiod = abs(job%bset(nu_i)%iperiod)
             rho_b(2) = rho_b(1)+(rho_b(2)-rho_b(1))/real(iperiod,ark)
             npoints = npoints/iperiod
           endif
           !
           allocate (f1drho(0:Npoints),g1drho(0:Npoints),weight(0:Npoints),stat=alloc)
           if (alloc/=0) then
              write (out,"(' Error ',i9,' trying to allocate f1drho and g1drho')") alloc
              stop 'FLbset1DNew, f1drho and g1drho - out of memory'
           end if
           !
           if (trove%periodic.and.job%bset(Nmodes)%iperiod>Nperiod_max) then
              stop 'FLbset1DNew: Nperiod_max is too small'
           endif
           !    
           ! Double check:
           !
           if (bs%imodes/=1.or.bs%mode(bs%imodes)/=trove%Nmodes) then
              write (out,"(' FLbset1DNew-Numerov: it has to be 1d and the last bs,')") 
              write (out,"('                      you have bs%imodes, bs%mode: ',30i8)") bs%imodes,bs%mode(:)
              stop 'FLbset1DNew, f1drho and g1drho - out of memory'
           end if
           !
           ! for basis sets we will need 1d potential and kinetic enrgy part 
           ! in terms of the corresponding coordinate 
           !
           ! Potential and kinetic energy parts
           !
           if (trove%lincoord/=0) then
             !
             ! Well, if pseudo -> infinity, we better switch it off. Here it is done by removing it
             ! when we have a naught at rho=0. 
             isingular = 0 
             !
           endif 
           !
           ! standard case of a non-singular pseudo-function
           !
           f1drho(0:npoints) = trove%poten%field(1,0:npoints)+trove%pseudo%field(1,0:npoints)
           ! singular case is reconstrcuted assuming the stored pseudo is pseudo*rho**2
           if (isingular>=0) then 
             !        
             rho_switch = trove%specparam(nu_i)
             iswitch = mod(nint( ( rho_switch-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
             !
             do i = iswitch+1,npoints
                !
                rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
                f1drho(i) = trove%poten%field(1,i) +trove%pseudo%field(1,i)/rho**2
                !
             enddo
             !
           endif 
           !
           ! for the kinetic part - we just take the corresoinding diagonal member of the g_vib%field
           !
           nu_i = trove%Nmodes ; fl => trove%g_vib(nu_i,nu_i)
           !
           g1drho(0:npoints) = fl%field(1,0:npoints)
           !
           reduced_model = .false.
           !
           if (bset%dscr(nu_i)%model<trove%NPotOrder) then 
             !
             reduced_model = .true.
             !
           endif 
           !
           if (.not.trove%DVR.or.reduced_model) then
             !
             ! for Krot /=0 in the BASIS input we add the corresponding diaginal rotational angular momentum term 
             ! g_rot(z,z)%field(1,:)*Krot**2
             !
             krot = 0 
             !
             if ( job%bset(0)%range(2)/=0 ) then 
               !
               krot = job%bset(0)%range(2)
               !
               fl => trove%g_rot(3,3)
               !
               f1drho(0:npoints) = f1drho(0:npoints)+0.5_ark*fl%field(1,0:npoints)*real(krot,ark)**2
               !
             endif 
             !
             if ( job%bset(0)%range(1)/=0 ) then 
               !
               jrot = job%bset(0)%range(1) ! job%bset(0)%model
               !
               f1drho(0:npoints) = f1drho(0:npoints)+0.25_ark*(trove%g_rot(1,1)%field(1,0:npoints)+&
                                   trove%g_rot(2,2)%field(1,0:npoints))* &
               real(jrot*(jrot+1)-krot**2,ark)
               !
             endif
             !
           endif
           !
           bs%params    = 0
           !
           ! We have stored the numeber of points, rhomax, and rhomin as optional parameters of  "bset%dscr"
           ! now we need them:
           !
           ! This NUMEROV for the last (could be also non-rigid) coordinate
           !
           allocate (drho(0:Npoints,3),stat=alloc)
           if (alloc/=0) then
              write (out,"(' Error ',i9,' trying to allocate drho')") alloc
              stop 'FLbset1DNew, drho - out of memory'
           end if
           !
           do i = 0,npoints
              !
              rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
              !
              drho(i,1) = MLcoord_direct(rho,1,nu_i)
              drho(i,2) = MLcoord_direct(rho,2,nu_i)
              drho(i,3) = MLcoord_direct(rho,3,nu_i)
              !
           enddo
           !
           weight = 1.0_ark
           !
           numerpoints = trove%numerpoints
           !
           if (trove%numerpoints<0) numerpoints = npoints
           !
           call ME_laguerre(npoints+1,bs%Size,bs%order,rho_b,drho(0:npoints,1:3),nu_i,isingular,f1drho,g1drho,job%verbose,&
                            bs%matelements,bs%ener0)
           !
        else
           !
           stop 'RIGID LAGUERRE NOT IMPLEMENTED'
           !
        endif
        !
     case('NUMEROV','BOX','FOURIER','LEGENDRE','SINRHO','LAGUERRE-K') 
        ! 
        ! numerov bset
        if (trove%manifold_rank(bs%mode(1))/=0) then
           !
           !if (trove%sparse) then 
           !  !
           !  write(out,"('FLbset1DNew: NON-RIGID was not tested in combinatiion for SPARSE, try either RIGID or NO-SPARSE ')") 
           !  !stop 'FLbset1DNew: NON-RIGID is not working for SPARSE yet'
           !  !
           !endif
           !
           !
           ! Allocation of the potential and kinetic 1d matrixes
           !
           if (job%bset(nu_i)%iperiod>0.and.trim(bs%type)=='NUMEROV') then
             iperiod = job%bset(nu_i)%iperiod
             rho_b(2) = rho_b(2)/real(iperiod,ark)
             npoints = npoints/iperiod
           endif
           !
           if (job%bset(nu_i)%iperiod==-2.and.trim(bs%type)=='NUMEROV') then
             iperiod = abs(job%bset(nu_i)%iperiod)
             rho_b(2) = rho_b(1)+(rho_b(2)-rho_b(1))/real(iperiod,ark)
             npoints = npoints/iperiod
           endif
           !
           allocate (f1drho(0:Npoints),g1drho(0:Npoints),weight(0:Npoints),stat=alloc)
           if (alloc/=0) then
              write (out,"(' Error ',i9,' trying to allocate f1drho and g1drho')") alloc
              stop 'FLbset1DNew, f1drho and g1drho - out of memory'
           end if
           !    
           ! Double check:
           !
           if (bs%imodes/=1.or.bs%mode(bs%imodes)/=trove%Nmodes) then
              write (out,"(' FLbset1DNew-Numerov: it has to be 1d and the last bs,')") 
              write (out,"('                      you have bs%imodes, bs%mode: ',30i8)") bs%imodes,bs%mode(:)
              stop 'FLbset1DNew, f1drho and g1drho - out of memory'
           end if
           !
           ! for basis sets we will need 1d potential and kinetic enrgy part 
           ! in terms of the corresponding coordinate 
           !
           ! Potential and kinetic energy parts
           !
           if (trove%lincoord/=0) then
             !
             ! Well, if pseudo -> infinity, we better switch it off. Here it is done by removing it
             ! when we have a naught at rho=0. 
             isingular = 0 
             !
           endif 
           !
           !if (trim(molec%coords_transform)=='R-RHO'.and..true.) then 
           !  !
           !  do i = 0,Npoints
           !   !
           !   rho_t = ((rho_b(1)+trove%rhostep*real(i,ark)))
           !   !
           !   rho_t = sum( trove%mass(:)*( trove%b0(:,2,i)**2+ trove%b0(:,3,i)**2 ) )/(planck*avogno*real(1.0d+16,kind=rk)/(4.0_ark*pi*pi*vellgt))
           !   !
           !   trove%poten%field(:,i) = trove%poten%field(:,i)*rho_t
           !   !
           ! enddo
           !  !
           !endif
           !
           !
           ! standard case of a non-singular pseudo-function
           !
           f1drho(0:npoints) = trove%poten%field(1,0:npoints)+trove%pseudo%field(1,0:npoints)
           !
           if (isingular>=0.and.(trim(bs%type)=='SINRHO'.or.trim(bs%type)=='LAGUERRE-K')) then 
               f1drho(0:npoints) = trove%poten%field(1,0:npoints)
           endif
           ! singular case is reconstrcuted assuming the stored pseudo is pseudo*rho**2
           if (isingular>=0.and.(.not.trim(bs%type)=='LEGENDRE'.and..not.trim(bs%type)=='SINRHO'.and.&
                                 .not.trim(bs%type)=='LAGUERRE-K')) then 
             !        
             rho_switch = trove%specparam(nu_i)
             iswitch = mod(nint( ( rho_switch-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
             !
             do i = iswitch+1,npoints
                !
                rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
                f1drho(i) = trove%poten%field(1,i) +trove%pseudo%field(1,i)/rho**2
                !
             enddo
             !
           endif 
           !
           ! These are grid-based corfinates 
           !
           allocate (drho(0:Npoints,3),xton(0:Npoints,0:trove%MaxOrder),stat=alloc)
           if (alloc/=0) then
              write (out,"(' Error ',i9,' trying to allocate drho')") alloc
              stop 'FLbset1DNew, drho - out of memory'
           end if
           !
           do i = 0,npoints
              !
              rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
              !
              drho(i,1) = MLcoord_direct(rho,1,nu_i)
              drho(i,2) = MLcoord_direct(rho,2,nu_i)
              drho(i,3) = MLcoord_direct(rho,3,nu_i)
              !
              do ipower = 0, trove%NKinorder
                 xton(i,ipower) = MLcoord_direct(rho,1,nu_i,ipower)
              enddo
              !
           enddo
           !
           ! for the kinetic part - we just take the corresoinding diagonal member of the g_vib%field
           !
           nu_i = trove%Nmodes ; fl => trove%g_vib(nu_i,nu_i)
           !
           g1drho(0:npoints) = fl%field(1,0:npoints)
           !
           if (trove%sparse) then
             !
             call find_isparse_from_ifull(fl%Ncoeff,fl%ifromsparse,1,i)
             !
             if (i==0) then 
                !
                g1drho = 0 
                !
                do icoeff = 1, fl%Ncoeff 
                   f_t = 1.0_ark
                   do imode  = 1,Nmodes
                     rho =  trove%chi_eq(imode)
                     ipower = fl%IndexQ(imode,icoeff)
                     rho_kin0 = MLcoord_direct(rho,1,imode,ipower)
                     f_t = f_t*rho_kin0
                  enddo
                  g1drho = g1drho + f_t*fl%field(icoeff,0:npoints)
                  !
                enddo
                !
             endif
             !
           elseif (abs(g1drho(trove%ipotmin))<small_) then 
             !
             g1drho = 0 
             !
             do icoeff = 1, fl%Ncoeff 
                f_t = 1.0_ark
                do imode  = 1,Nmodes
                  rho =  trove%chi_eq(imode)
                  ipower = fl%IndexQ(imode,icoeff)
                  rho_kin0 = MLcoord_direct(rho,1,imode,ipower)
                  f_t = f_t*rho_kin0
               enddo
               g1drho = g1drho + f_t*fl%field(icoeff,0:npoints)
               !
             enddo
             !
           endif
           !
           reduced_model = .false.
           !
           if (bset%dscr(nu_i)%model<trove%NPotOrder) then 
             !
             reduced_model = .true.
             !
           endif 
           !
           if (trove%DVR.and..not.reduced_model) then
             !
             chi = 0
             !
             do i = 0,Npoints
               !
               chi(nu_i) = rho_b(1)+trove%rhostep*real(i,ark)-trove%chi0_ref(nu_i)
               !
               call FLcalc_poten_kinet_dvr(chi,i,poten_t,gvib_t,grot_t,gcor_t,extF_t,reduced_model)
               !
               if (job%verbose>=7) write(out,"(i8,f12.6,3g18.8)") i,chi(nu_i),poten_t,gvib_t(nu_i,nu_i)
               !
               g1drho(i) = gvib_t(nu_i,nu_i)
               f1drho(i) = poten_t
               !
             enddo
             !
           endif
           !
           if ( (.not.trim(bs%type)=='LEGENDRE'.and..not.trim(bs%type)=='SINRHO'.and..not.trim(bs%type)=='LAGUERRE-K').and.&
                 .not.trove%DVR.or.reduced_model) then
             !
             ! for Krot /=0 in the BASIS input we add the corresponding diaginal rotational angular momentum term 
             ! g_rot(z,z)%field(1,:)*Krot**2
             !
             krot = 0 
             !
             if ( job%bset(0)%range(2)/=0 ) then 
               !
               krot = job%bset(0)%range(2)
               !
               fl => trove%g_rot(3,3)
               !
               f1drho(0:npoints) = f1drho(0:npoints)+0.5_ark*fl%field(1,0:npoints)*real(krot,ark)**2
               !
             endif 
             !
             if ( job%bset(0)%range(1)/=0 ) then 
               !
               jrot = job%bset(0)%range(1) ! job%bset(0)%model
               !
               f1drho(0:npoints) = f1drho(0:npoints)+0.25_ark*(trove%g_rot(1,1)%field(1,0:npoints)+&
                                   trove%g_rot(2,2)%field(1,0:npoints))*real(jrot*(jrot+1)-krot**2,ark)
               !
             endif
             !
           endif
           !
           bs%params    = 0
           !
           ! We have stored the numeber of points, rhomax, and rhomin as optional parameters of  "bset%dscr"
           ! now we need them:
           !
           weight = 1.0_ark
           !
           numerpoints = trove%numerpoints
           !
           if (trove%numerpoints<0) numerpoints = npoints
           !
           if (trim(bs%type)=='NUMEROV') then
             !
             call ME_numerov(bs%Size,bs%order,rho_b,isingular,npoints,numerpoints,drho,xton,f1drho,g1drho,nu_i,&
                             job%bset(nu_i)%iperiod,job%verbose,bs%matelements,bs%ener0)
             !
           elseif (trim(bs%type)=='BOX') then 
             !
             call ME_box(bs%Size,bs%order,rho_b,isingular,npoints,drho,f1drho(0:npoints),g1drho(0:npoints),nu_i,&
                         job%bset(nu_i)%periodic,job%verbose,&
                         bs%matelements(-1:3,0:trove%MaxOrder,0:bs%Size,0:bs%Size),bs%ener0(0:bs%Size))
             !
           elseif (trim(bs%type)=='LEGENDRE') then 
             !
             kmax = job%bset(0)%range(2)
             nmax = bs%Size
             if ( kmax/=0 ) then 
               nmax = (bs%Size+1)/(kmax+1)-1
             endif
             !
             allocate (muzz(0:Npoints),sinrho(0:Npoints),cosrho(0:Npoints),stat=alloc)
             if (alloc/=0) then
                write (out,"(' Error ',i9,' trying to allocate muzz')") alloc
                stop 'FLbset1DNew, muzz - out of memory'
             end if
             !
             muzz = trove%g_rot(3,3)%field(1,0:npoints)
             !
             !call ME_Legendre(bs%Size,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,nu_i,job%verbose,bs%matelements,bs%ener0)
             !
             call ME_Associate_Legendre(bs%Size,kmax,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,muzz,nu_i,&
                                        job%verbose,bs%matelements,bs%ener0)
             !
             !call ME_sinrho_polynomial(bs%Size,kmax,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,muzz,nu_i,&
             !                          job%verbose,bs%matelements,bs%ener0)
             !
             do i = 0,npoints
                rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
                sinrho(i) = sin(rho)
                cosrho(i) = cos(rho)
             enddo
             !
             deallocate(muzz)
             !
           elseif (trim(bs%type)=='SINRHO') then 
             !
             kmax = job%bset(0)%range(2)
             nmax = bs%Size
             if ( kmax/=0 ) then 
               nmax = (bs%Size+1)/(kmax+1)-1
             endif
             !
             allocate (muzz(0:Npoints),sinrho(0:Npoints),cosrho(0:Npoints),pseudo(0:Npoints),stat=alloc)
             if (alloc/=0) then
                write (out,"(' Error ',i9,' trying to allocate muzz')") alloc
                stop 'FLbset1DNew, muzz - out of memory'
             end if
             !
             muzz = trove%g_rot(3,3)%field(1,0:npoints)
             pseudo = trove%pseudo%field(1,0:npoints)
             !
             !call ME_sinrho_polynomial_k(bs%Size,kmax,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,muzz,nu_i,&
             !                          job%verbose,bs%matelements,bs%ener0)
                                       !
             !call ME_sinrho_polynomial_k(bs%Size,kmax,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,muzz,pseudo,nu_i,&
             !                          job%verbose,bs%matelements,bs%ener0)
                                       !
             call ME_legendre_polynomial_k(bs%Size,kmax,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,muzz,pseudo,nu_i,&
                                       job%verbose,bs%matelements,bs%ener0)
                                       !
             !call ME_sinrho_polynomial_muzz(bs%Size,kmax,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,muzz,pseudo,nu_i,&
             !                          job%verbose,bs%matelements,bs%ener0)

             !
             do i = 0,npoints
                rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
                sinrho(i) = sin(rho)
                cosrho(i) = cos(rho)
             enddo
             !
             deallocate(muzz,pseudo)
             !
           elseif (trim(bs%type)=='LAGUERRE-K') then 
             !
             kmax = job%bset(0)%range(2)
             nmax = bs%Size
             if ( kmax/=0 ) then 
               nmax = (bs%Size+1)/(kmax+1)-1
             endif
             !
             allocate (muzz(0:Npoints),mrho(0:Npoints),pseudo(0:Npoints),stat=alloc)
             if (alloc/=0) then
                write (out,"(' Error ',i9,' trying to allocate muzz')") alloc
                stop 'FLbset1DNew, muzz - out of memory'
             end if
             !
             muzz = trove%g_rot(3,3)%field(1,0:npoints)
             pseudo = trove%pseudo%field(1,0:npoints)
             !
             fl => trove%g_rot(3,3)
             gl => trove%pseudo
             !
             if (trove%sparse) then
               call find_isparse_from_ifull(fl%Ncoeff,fl%ifromsparse,1,i1)
               call find_isparse_from_ifull(fl%Ncoeff,fl%ifromsparse,1,i2)
               !
               if (i1==0.or.i2==0) then 
                  !
                  muzz = 0 
                  pseudo = 0
                  !
                  do icoeff = 1, fl%Ncoeff 
                     f_t = 1.0_ark
                     do imode  = 1,Nmodes
                       rho =  trove%chi_eq(imode)
                       ipower = fl%IndexQ(imode,icoeff)
                       rho_kin0 = MLcoord_direct(rho,1,imode,ipower)
                       f_t = f_t*rho_kin0
                    enddo
                    muzz   = muzz   + f_t*fl%field(icoeff,0:npoints)
                    pseudo = pseudo + f_t*gl%field(icoeff,0:npoints)
                    !
                  enddo
                  !
               endif
               !
             elseif(abs(muzz(trove%ipotmin))<small_.or.abs(pseudo(trove%ipotmin))<small_) then
               !
               muzz = 0 
               pseudo = 0
               !
               do icoeff = 1, fl%Ncoeff 
                  f_t = 1.0_ark
                  do imode  = 1,Nmodes
                    rho =  trove%chi_eq(imode)
                    ipower = fl%IndexQ(imode,icoeff)
                    rho_kin0 = MLcoord_direct(rho,1,imode,ipower)
                    f_t = f_t*rho_kin0
                 enddo
                 muzz   = muzz   + f_t*fl%field(icoeff,0:npoints)
                 pseudo = pseudo + f_t*gl%field(icoeff,0:npoints)
                 !
               enddo
               !
             endif
             !
             g_t = g1drho(trove%ipotmin)
             !
             if (g_t<small_) then 
               write(out,"('At ME_laguerre_k: mu_rr(imin) cannot be zero ',g18.8)") g_t
               stop 'At ME_laguerre_k: illegal mu_rr(imin)'
             endif
             !
             if (trove%ipotmin>0.and.trove%ipotmin<npoints) then
               f_t = ( f1drho(trove%ipotmin+1)+f1drho(trove%ipotmin-1)-2.0_ark*f1drho(trove%ipotmin) )/trove%rhostep**2 
             elseif (trove%ipotmin == 0 ) then
               f_t = ( 2.0_ark*f1drho(2)-2.0_ark*f1drho(0) )/(trove%rhostep*2.0_ark)**2 
             elseif (trove%ipotmin == npoints ) then
               f_t = ( 2.0_ark*f1drho(npoints-1)-2.0_ark*f1drho(npoints) )/trove%rhostep**2 
             else 
               stop 'At ME_laguerre_k: illegal imin'
             endif
             !
             if (trove%specparam(bs%mode(1))>0d0 ) then
                !
                f_t = trove%specparam(bs%mode(1))**2*g_t
                !
                if (job%verbose>=5) then
                  write(out,"('the input special-parameter ',f18.8,' will be used to obtain m for laguarre basis')") & 
                        trove%specparam(bs%mode(1)) 
                endif
                !
             endif
             !
             f_m = sqrt(f_t/g_t)
             !
             call ME_laguerre_k(bs%Size,kmax,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,muzz,f_m,pseudo,nu_i,&
                                       job%verbose,bs%matelements,bs%ener0)
             !
             !call ME_laguerre_simple_k(bs%Size,kmax,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,muzz,f_m,pseudo,nu_i,&
             !                          job%verbose,bs%matelements,bs%ener0)
             !
             do i = 0,npoints
                rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
                mrho(i) = rho
             enddo
             !
             deallocate(muzz,pseudo)             
             !
           elseif (trim(bs%type)=='FOURIER') then 
             !
             call ME_Fourier(bs%Size,bs%order,rho_b,isingular,npoints,numerpoints,drho,f1drho,g1drho,nu_i,&
                             job%bset(nu_i)%iperiod,job%verbose,bs%matelements,bs%ener0)
             !
           endif
           !
           allocate(phil(0:trove%Npoints),phir(0:trove%Npoints),dphil(0:trove%Npoints),dphir(0:trove%Npoints), &
                 phivphi(0:trove%Npoints),stat=alloc)
           if (alloc/=0) then 
              write (out,"('FLbset1DNew: phi - out of memory')")
              stop 'FLbset1DNew: phi - out of memory'
           endif 
           !
           if (job%bset(nu_i)%iperiod/=0.and.trim(bs%type)=='NUMEROV') then
             !
             if (abs(job%bset(nu_i)%iperiod)/=2) then
               write (out,"(' FLbset1DNew warning: periodic copying procedure has been only tested for period of 2 ')") 
               !stop 'FLbset1DNew: illegal periodicity'
             endif 
             !
             write(unitfname,"('Numerov basis set # ',i6)") nu_i
             call IOStart(trim(unitfname),io_slot)
             !
             allocate (func(0:bs%Size+1,0:npoints),dfunc(0:bs%Size+1,0:npoints),stat=alloc)
             call ArrayStart('numerov-bs_funct' ,alloc,size(func),kind(func))
             call ArrayStart('numerov-bs_funct' ,alloc,size(dfunc),kind(dfunc))
             !
             do vl = 0,bs%Size+1
                !
                read (io_slot,rec=vl+1) (func(vl,k),k=0,npoints),(dfunc(vl,k),k=0,npoints)
                !
             enddo
             !
             close(io_slot)
             !
             inquire(iolength=rec_len) phil(:),dphil(:)
             !
             open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
             !
             if (abs(job%bset(nu_i)%iperiod)==2) then
               !
               sqrt2 = sqrt(0.5_ark)
               !
               do vl = 0,bs%Size+1
                  !
                  forall(i = 0:npoints) phil(i)  =  func(vl,i)*sqrt2
                  forall(i = 0:npoints) dphil(i) = dfunc(vl,i)*sqrt2
                  !
                  if (job%bset(nu_i)%iperiod>0) then
                    !
                    iparity = 2
                    !
                    if (mod((vl+1)/2,2)==1) iparity = 1
                    !
                    !if (iparity==1.and.abs(func(vl,npoints))>sqrt(small_)) then
                    !   !iparity = 1
                    !   write(out,'("basis_parity: parity = 1, but psi(npoints)/=0;  mode = ",i2," v = ",i3," psi(npoints) ")') nu_i,vl,func(vl,npoints)
                    !   stop "psi must be zero at its node for parity 1"
                    !endif
                    !
                    forall(i = 1:npoints) phil(npoints+i)  =  func(vl,npoints-i)*(-1.0_ark)**(iparity  )*sqrt2
                    forall(i = 1:npoints) dphil(npoints+i) = dfunc(vl,npoints-i)*(-1.0_ark)**(iparity+1)*sqrt2
                    !
                  else
                    !
                    iparity = 2
                    if (mod(vl,2)/=0) iparity = 1
                    !
                    forall(i = 1:npoints) phil(npoints+i)  =  func(vl,npoints-i)*(-1.0_ark)**(iparity  )*sqrt2
                    forall(i = 1:npoints) dphil(npoints+i) = dfunc(vl,npoints-i)*(-1.0_ark)**(iparity+1)*sqrt2
                    !
                  endif 
                  !
                  write (io_slot,rec=vl+1) (phil(i),i=0,trove%Npoints),(dphil(i),i=0,trove%Npoints)
                  !
               enddo
               !
             else
               !
               ! Bloch functions 
               !
               L = rho_range*0.5_ark
               !
               step = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
               !
               nl = 0 
               !
               if (.true.) then
                 !
                 sqrt2 = 1.0_ark/sqrt(real(job%bset(nu_i)%iperiod,ark))
                 !
                 loop_v : do vl = 0,bs%Size+1
                    !
                    !iparity = 2
                    !if (mod((vl+1)/2,2)==1) iparity = 1
                    !
                    !k = mod((vl+1)/2,job%bset(nu_i)%iperiod)
                    !
                    do i=0,trove%Npoints
                       !
                       irho = mod(i,Npoints)
                       !
                       rho = real(i,kind=ark)*step
                       phil(i)  = func(vl,irho)*sqrt2
                       dphil(i) = dfunc(vl,irho)*sqrt2
                       !
                    enddo
                    !
                    phivphi(:) = phil(:)*phil(:)
                    !
                    mat_t = simpsonintegral_ark(trove%Npoints,rho_range,phivphi)
                    !
                    phil(:) = 1.0_ark/sqrt(mat_t)*phil(:)
                    dphil(:) = 1.0_ark/sqrt(mat_t)*dphil(:)
                    !
                    write (io_slot,rec=vl+1) (phil(i),i=0,trove%Npoints),(dphil(i),i=0,trove%Npoints)
                    !
                 enddo loop_v
                 !
               else
                 !
                 loop_vl : do vl = 0,bs%Size+1
                    !
                    do k =0,job%bset(nu_i)%iperiod-1
                      !
                      nl = nl + 1
                      !
                      if (nl>bs%Size+1) cycle loop_vl 
                      !
                      !k = mod((vl+1)/2,job%bset(nu_i)%iperiod)
                      !
                      do i=0,trove%Npoints
                         !
                         irho = mod(i,Npoints)
                         !
                         rho = real(i,kind=ark)*step
                         !
                         if (k==0) then 
                           phil(i)  = func(vl,irho)
                           dphil(i) = dfunc(vl,irho)
                         elseif (mod(k,2)==0) then
                           phil(i)  = cos(real(k,ark)*pi*rho/L)*func(vl,irho)
                           dphil(i) = cos(real(k,ark)*pi*rho/L)*dfunc(vl,irho)-sin(real(k,ark)*pi*rho/L)*real(k,ark)*pi/&
                                      L*func(vl,irho)
                         else
                           phil(i)  = sin(real(k,ark)*pi*rho/L)*func(vl,irho)
                           dphil(i) = sin(real(k,ark)*pi*rho/L)*dfunc(vl,irho)+cos(real(k,ark)*pi*rho/L)*real(k,ark)*pi/&
                                      L*func(vl,irho)
                         endif
                         !
                      enddo
                      !
                      phivphi(:) = phil(:)*phil(:)
                      !
                      mat_t = simpsonintegral_ark(trove%Npoints,rho_range,phivphi)
                      !
                      phil(:) = 1.0_ark/sqrt(mat_t)*phil(:)
                      dphil(:) = 1.0_ark/sqrt(mat_t)*dphil(:)
                      !
                      write (io_slot,rec=nl+1) (phil(i),i=0,trove%Npoints),(dphil(i),i=0,trove%Npoints)
                      !
                    enddo
                    !
                 enddo loop_vl
                 !
               endif
               !
             endif 
             !
             npoints = trove%Npoints
             rho_b  = trove%rho_border
             !
             deallocate(func,dfunc)
             call ArrayStop('numerov-bs_funct')
             !
           endif 
           !
           if (job%verbose>=3) write(out,"(/'Primitive matrix elements...')")
           !
           call TimerStart('Primitive matrix elements')
           !
           deallocate(phil,phir,dphil,dphir)
           !
           if (trim(bs%type)=='LEGENDRE') then
             !
             !$omp parallel private(phivphi_t,phil,phir,dphil,dphir,alloc_p,phil_leg,phir_leg,&
             !$omp& dphil_leg,dphir_leg)
             allocate (phivphi_t(0:npoints),phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints),stat=alloc_p)
             if (alloc_p/=0) then
                write (out,"(' FLbset1DNew error ',i9,' trying to allocate 11 arrays phivphi_t,phi, etc')") alloc_p
                stop 'FLbset1DNew, phivphi_t  - out of memory'
             end if
             !
             allocate(phil_leg(0:trove%Npoints),phir_leg(0:trove%Npoints),&
                      dphil_leg(0:trove%Npoints),dphir_leg(0:trove%Npoints),stat=alloc_p)
             if (alloc_p/=0) then 
                write (out,"('FLbset1DNew: phil_leg - out of memory')")
                stop 'FLbset1DNew: phil_leg - out of memory'
             endif 
             !
             !$omp do private(vl,unitfname,io_slot,i,rho,k,vr,Tcoeff,iterm,k1,k2,imu,mat_t) schedule(dynamic)
             do vl = 0,bs%Size
                !
                write(unitfname,"('Numerov basis set # ',i6)") nu_i
                call IOStart(trim(unitfname),io_slot)
                !
                read (io_slot,rec=vl+1) (phil(k),k=0,npoints),(dphil(k),k=0,npoints)
                !
                krot1 = mod(vl,kmax+1)
                nl = (vl-krot1)/(kmax+1)
                !
                do i = 0,npoints
                   rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
                   !
                   !!phil_leg(i) = phil(i)
                   !!phil(i) = sqrt(sin(rho))*phil(i)*sin(rho)**krot1
                   !!dphil_leg(i) = cos(rho)*0.5_ark*phil_leg(i)+sin(rho)*dphil(i)
                   !
                   phil_leg(i) = phil(i)
                   phil(i) = sqrt(sin(rho))*phil(i)
                   dphil_leg(i) = cos(rho)*0.5_ark*phil_leg(i)-sin(rho)*dphil(i)-cos(rho)*phil_leg(i)*real(krot1,ark)
                   !
                   !
                   !dphil(i) = dphil(i)
                enddo
                !
                do vr = 0,bs%Size
                    !
                    read (io_slot,rec=vr+1) (phir(k),k=0,npoints),(dphir(k),k=0,npoints)
                    !
                    krot2 = mod(vr,kmax+1)
                    nr = (vr-krot2)/(kmax+1)
                    !
                    !if ( krot2/=krot1 ) cycle
                    !
                    !  write(out,"('inconsistent k in vr',2i8,' for vr,nr = ',2i5)") krot,(vr-nr-1)/(nmax+1),vr,nr
                    !  stop 'Inconsistent kl and kr'
                    !endif
                    !
                    do i = 0,npoints
                       rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
                       !
                       !!phir_leg(i) = phir(i)
                       !!phir(i) = sqrt(sin(rho))*phir(i)*sin(rho)**krot2
                       !!dphir_leg(i) = cos(rho)*0.5_ark*phir_leg(i)+sin(rho)*dphir(i)
                       !
                       phir_leg(i) = phir(i)
                       phir(i) = sqrt(sin(rho))*phir(i)
                       dphir_leg(i) = cos(rho)*0.5_ark*phir_leg(i)-sin(rho)*dphir(i)-cos(rho)*phir_leg(i)*real(krot2,ark)
                       !
                       !dphir(i) = sqrt(sin(rho))*dphir(i)
                    enddo
                    !
                    if (krot1==krot2) then
                      !
                      Tcoeff = trove%poten%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                        !
                        phivphi_t(:) = phil(:)*trove%poten%field(iterm,:)*phir(:)
                        !
                        trove%poten%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                        !
                      enddo
                      !
                      ! vibraional kinetic energy part
                      !
                      do k1 = 1,Nmodes
                         !
                         do k2 = 1,Nmodes
                           !
                           fl => trove%g_vib(k1,k2) 
                           !
                           Tcoeff = fl%Ncoeff
                           !
                           do iterm = 1,Tcoeff
                              !
                              if (k1==Nmodes.and.k2==Nmodes) then 
                                 !
                                 !phivphi_t(:) =-dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                                 phivphi_t(:) =-fl%field(iterm,:)*( dphil(:)*dphiR(:)*sinrho(:)-&
                                 (real(krot2,ark)*dphil(:)*phir_leg(:)+real(krot1,ark)*phil_leg(:)*dphir(:) )*cosrho(:) )
                                 !
                              elseif (k1==Nmodes) then 
                                 !
                                 !phivphi_t(:) =-dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                                 phivphi_t(:) = -dphil_leg(:)*fl%field(iterm,:)*phir_leg(:)
                                 !
                              elseif (k2==Nmodes) then 
                                 !
                                 !phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                                 phivphi_t(:) = phil_leg(:)*fl%field(iterm,:)*dphir_leg(:)
                                 !
                              else
                                 !
                                 phivphi_t(:) = phil(:)*fl%field(iterm,:)*phir(:)
                                 !
                              endif
                              !
                              trove%g_vib(k1,k2)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                              !
                           enddo
                           !
                        enddo
                        !
                      enddo
                      !
                      Tcoeff = trove%pseudo%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                         !
                         ! The gvib term can be combined with the pseudopotential part.
                         ! It will allow us to spare additional array and i.e. memory
                         ! and it is also necessary in case of the singular solution to do so.
                         !
                         !phivphi_t(:) =-2.0_ark*phil(:)*trove%pseudo%field(iterm,:)*phir(:)
                         !
                         phivphi_t(:) =-2.0_ark*phil_leg(:)*trove%pseudo%field(iterm,:)*phir_leg(:)
                         !
                         !!!phivphi_t(:) =-2.0_ark*phil_leg(:)*trove%pseudo%field(iterm,:)*phir_leg(:)*sinrho(:)**(krot1+krot2)
                         !
                         mat_t = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                         !
                         trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr)+mat_t
                         !
                      enddo
                      !
                      ! Vibraional Angular Momentum L2
                      !
                      if (FLl2_coeffs) then
                        !
                        stop 'FLl2_coeffs is not implemented for Legendre'
                        !
                        do k1 = 1,Nmodes
                          !
                          do k2 = 1,Nmodes
                             !
                             Tcoeff = trove%L2_vib(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                if (k1==Nmodes.and.k2==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                elseif (k1==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                elseif (k2==Nmodes) then 
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                else
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                endif
                                !
                                trove%L2_vib(k1,k2)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                          enddo
                          !
                        enddo
                        !
                      endif
                      !
                    endif
                    !
                    if (FLrotation) then
                       !
                       ! rotation kinetic energy part
                       !
                       !
                       ! Special basis set sqrt(sin(rho))L_n
                       !
                       do k1 = 1,3
                         do k2 = 1,3
                           !
                           fl => trove%g_rot(k1,k2)
                           !
                           if (k1==3.and.k2==3) then
                             !
                             do iterm = 1,fl%Ncoeff
                               !!phivphi_t(:) = phil_leg(:)*fl%field(iterm,:)*phir_leg(:)*sinrho(:)**(krot1+krot2-1)
                               !
                               phivphi_t(:) = -phil_leg(:)*fl%field(iterm,:)*phir_leg(:)
                               !
                               fl%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                               !
                             enddo
                             !
                           elseif (k1==3.or.k2==3) then
                             !
                             do iterm = 1,fl%Ncoeff
                               !!phivphi_t(:) =phil_leg(:)*fl%field(iterm,:)*phir_leg(:)*sinrho(:)**(krot1+krot2)
                               !
                               phivphi_t(:) =phil_leg(:)*fl%field(iterm,:)*phir_leg(:)
                               fl%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                             enddo
                             !
                           else
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                phivphi_t(:) = phil(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir(:)
                                !
                                trove%g_rot(k1,k2)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                           endif
                           !
                         enddo
                         !
                       enddo
                       !
                       !fl => trove%g_rot(3,3)
                       !!fl%me(:,vl,vr) = 0 
                       !do iterm = 1,fl%Ncoeff
                       !  phivphi_t(:) = phil_leg(:)*fl%field(iterm,:)*phir_leg(:)*sinrho(:)**(2*k)
                       !  fl%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                       !enddo
                       !
                       ! Coriolis part
                       !
                       do k1 = 1,Nmodes
                         !
                         do k2 = 1,3
                            !
                            fl => trove%g_cor(Nmodes,k2)
                            !
                            Tcoeff = fl%Ncoeff
                            !
                            do iterm = 1,Tcoeff
                               !
                               if (k1==Nmodes) then 
                                  !
                                  phivphi_t(:) = fl%field(iterm,:)*( phil_leg(:)*dphir_leg(:) - dphil_leg(:)*phir_leg(:) )
                                  !
                               else
                                  !
                                  phivphi_t(:) = phil(:)*fl%field(iterm,:)*phir(:)
                                  !
                               endif
                               !
                              fl%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                               !
                               !
                            enddo
                            !
                         enddo
                         !
                       enddo
                       !
                    endif
                    !
                    if (FLextF_coeffs) then
                      !
                      ! External field part
                      !
                      do imu = 1,extF%rank
                        !
                        Tcoeff = trove%extF(imu)%Ncoeff
                        !
                        do iterm = 1,Tcoeff
                           !
                           phivphi_t(:) = phil(:)*trove%extF(imu)%field(iterm,:)*phir(:)
                           !
                           trove%extF(imu)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                           !
                        enddo
                        !
                      enddo
                      !
                    endif 
                    !
                enddo 
                !
             enddo
             !$omp end do
             !
             deallocate (phivphi_t,phil,phir,dphil,dphir,phil_leg,phir_leg,phil_sin,phir_sin,dphil_leg,dphir_leg)
             !$omp end parallel 
             !

           elseif (trim(bs%type)=='SINRHO-MUZZ') then
             !
             if (job%verbose>=4) then 
                write(out,"('   Allocating 11 arrays of ',i8,', Ncore times ',g12.4,' ')") trove%Npoints+1,&
                                11.0_rk*real(trove%Npoints+1,rk)/1024.0_rk**3
             endif
             !
             !do i = 0,npoints
             !   rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
             !   !
             !   if (rho>0.5_ark*pi) then
             !     sinrho(i) = 1.0_ark
             !     cosrho(i) = 0
             !   else
             !     sinrho(i) = sin(rho)
             !     cosrho(i) = cos(rho)
             !   endif
             !   !
             !enddo
             !
             !$omp parallel private(phivphi_t,phil,phir,dphil,dphir,alloc_p,phil_leg,phir_leg,&
             !$omp& phil_sin,phir_sin,dphil_leg,dphir_leg)
             allocate (phivphi_t(0:npoints),phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints),&
                       phil_leg(0:npoints),phir_leg(0:npoints),phil_sin(0:npoints),phir_sin(0:npoints),&
                       dphil_leg(0:npoints),dphir_leg(0:npoints),stat=alloc_p)
             if (alloc_p/=0) then
                write (out,"(' FLbset1DNew error ',i9,' trying to allocate 11 arrays phivphi_t,phi, etc')") alloc_p
                stop 'FLbset1DNew, phivphi_t  - out of memory'
             end if
             !
             !$omp do private(vl,unitfname,io_slot,nl,krot1,vr,nr,krot2,Tcoeff,iterm,k1,k2,imu,mat_t) schedule(dynamic)
             do vl = 0,bs%Size
                !
                write(unitfname,"('Numerov basis set # ',i6)") nu_i
                call IOStart(trim(unitfname),io_slot)
                !
                read (io_slot,rec=vl+1) (phil_leg(k),k=0,npoints),(dphil_leg(k),k=0,npoints)
                !
                krot1 = mod(vl,kmax+1)
                nl = (vl-krot1)/(kmax+1)
                !krot11 = 0 ; if (krot1>0) krot11 = 1
                !
                phil_sin(:) = phil_leg(:)*sinrho(:)**krot1
                phil(:) = sqrt(sinrho(:))*phil_sin(:)
                !
                dphil(:) =cosrho(:)*0.5_ark*phil_leg(:)+sinrho(:)*dphil_leg(:)
                !
                dphil(:) = dphil(:)+cosrho(:)*phil_leg(:)*real(krot1,ark)
                !
                do vr = 0,bs%Size
                    !
                    krot2 = mod(vr,kmax+1)
                    nr = (vr-krot2)/(kmax+1)
                    if (abs(krot1-krot2)>2) cycle
                    read (io_slot,rec=vr+1) (phir_leg(k),k=0,npoints),(dphir_leg(k),k=0,npoints)
                    !
                    !krot21 = 0 ; if (krot2>0) krot21 = 1
                    !
                    !if ( krot2/=krot1 ) cycle
                    !
                    phir_sin(:) = phir_leg(:)*sinrho(:)**krot2
                    phir(:) = sqrt(sinrho(:))*phir_sin(:)
                    dphir(:) =cosrho(:)*0.5_ark*phir_leg(:)+sinrho(:)*dphir_leg(:)
                    !
                    dphir(:) = dphir(:)+cosrho(:)*phir_leg(:)*real(krot2,ark)
                    !
                    if (krot1==krot2) then
                      !
                      Tcoeff = trove%poten%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                        !
                        phivphi_t(:) = phil(:)*trove%poten%field(iterm,:)*phir(:)
                        !
                        trove%poten%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                        !
                      enddo
                      !
                      ! vibraional kinetic energy part
                      !
                      do k1 = 1,Nmodes
                         !
                         do k2 = 1,Nmodes
                           !
                           Tcoeff = trove%g_vib(k1,k2)%Ncoeff
                           !
                           do iterm = 1,Tcoeff
                              !
                              if (k1==Nmodes.and.k2==Nmodes) then 
                                 !
                                 phivphi_t(:) =-dphil_leg(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir_leg(:)*sinrho(:)
                                 !
                                 phivphi_t(:) = phivphi_t(:) - &
                                                trove%g_vib(k1,k2)%field(iterm,:)*cosrho(:)*real(krot1,ark)*&
                                                ( dphil_leg(:)*phir_sin(:)+phil_sin(:)*dphir_leg(:) )
                                 !
                              elseif (k1==Nmodes) then 
                                 !
                                 !phivphi_t(:) =-dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                                 phivphi_t(:) = -dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir_sin(:)
                                 !
                              elseif (k2==Nmodes) then 
                                 !
                                 !phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                                 phivphi_t(:) = phil_sin(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                              else
                                 !
                                 phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                              endif
                              !
                              trove%g_vib(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                              !
                           enddo
                           !
                        enddo
                        !
                      enddo
                      !
                      Tcoeff = trove%pseudo%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                         !
                         ! The gvib term can be combined with the pseudopotential part.
                         ! It will allow us to spare additional array and i.e. memory
                         ! and it is also necessary in case of the singular solution to do so.
                         !
                         !!phivphi_t(:) =-2.0_ark*phil(:)*trove%pseudo%field(iterm,:)*phir(:)
                         !
                         !!phivphi_t(:) =-2.0_ark*phil_leg(:)*trove%pseudo%field(iterm,:)*phir_leg(:)
                         !
                         phivphi_t(:) =-2.0_ark*phil_sin(:)*trove%pseudo%field(iterm,:)*phir_sin(:)
                         !
                         mat_t = integral_rect_ark(npoints,rho_range,phivphi_t)
                         !
                         trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr)+mat_t
                         !
                      enddo
                      !
                      ! Vibraional Angular Momentum L2
                      !
                      if (FLl2_coeffs) then
                        !
                        stop 'FLl2_coeffs is not implemented for Legendre'
                        !
                        do k1 = 1,Nmodes
                          !
                          do k2 = 1,Nmodes
                             !
                             Tcoeff = trove%L2_vib(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                if (k1==Nmodes.and.k2==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                elseif (k1==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                elseif (k2==Nmodes) then 
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                else
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                endif
                                !
                                trove%L2_vib(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                          enddo
                          !
                        enddo
                        !
                      endif
                      !
                    endif
                    !
                    if (FLrotation) then
                       !
                       ! rotation kinetic energy part
                       !
                       ! Special basis set sqrt(sin(rho))^(k+1/2) phi_n
                       !
                       do k1 = 1,3
                         do k2 = 1,3
                           !
                           if (k1==3.and.k2==3) then
                             !
                             if (krot1+krot2==0.or.krot1/=krot2) cycle
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                               !
                               phivphi_t(:) = phil_leg(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_leg(:)*sinrho(:)**(2*krot1-1)
                               !
                               !trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                               mat_t = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                               trove%poten%me(iterm,vl,vr) = trove%poten%me(iterm,vl,vr) + 0.5_ark*mat_t*real(krot1**2,ark)
                               !
                             enddo
                             !
                           elseif (k1==3.or.k2==3) then
                             !
                             if (krot1+krot2==0.or.abs(krot1-krot2)>1) cycle
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                               !
                               phivphi_t(:) =phil_sin(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_sin(:)
                               !
                               trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                             enddo
                             !
                           else
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                phivphi_t(:) = phil(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir(:)
                                !
                                trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                           endif
                           !
                         enddo
                         !
                       enddo
                       !
                       ! Coriolis part
                       !
                       do k1 = 1,Nmodes
                         !
                         do k2 = 1,3
                            !
                            Tcoeff = trove%g_cor(k1,k2)%Ncoeff
                            !
                            do iterm = 1,Tcoeff
                               !
                               if (k1==Nmodes) then 
                                  !
                                  phivphi_t(:) = trove%g_cor(k1,k2)%field(iterm,:)*( phil_sin(:)*dphir(:) - dphil(:)*phir_sin(:))
                                  !
                               else
                                  !
                                  phivphi_t(:) = phil(:)*trove%g_cor(k1,k2)%field(iterm,:)*phir(:)
                                  !
                               endif
                               !
                               trove%g_cor(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                               !
                            enddo
                            !
                         enddo
                         !
                       enddo
                       !
                    endif
                    !
                    if (FLextF_coeffs) then
                      !
                      ! External field part
                      !
                      do imu = 1,extF%rank
                        !
                        Tcoeff = trove%extF(imu)%Ncoeff
                        !
                        do iterm = 1,Tcoeff
                           !
                           phivphi_t(:) = phil(:)*trove%extF(imu)%field(iterm,:)*phir(:)
                           !
                           trove%extF(imu)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                           !
                        enddo
                        !
                      enddo
                      !
                    endif 
                    !
                enddo 
                !
             enddo
             !$omp end do
             !
             deallocate (phivphi_t,phil,phir,dphil,dphir,phil_leg,phir_leg,phil_sin,phir_sin,dphil_leg,dphir_leg)
             !$omp end parallel 
             !
             deallocate(sinrho,cosrho,stat=alloc_p)
             !
           elseif (trim(bs%type)=='SINRHO') then
             !
             if (job%verbose>=4) then 
                write(out,"('   Allocating 11 arrays of ',i8,', Ncore times ',g12.4,' ')") trove%Npoints+1,&
                                11.0_rk*real(trove%Npoints+1,rk)/1024.0_rk**3
             endif
             !
             !do i = 0,npoints
             !   rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
             !   !
             !   if (rho>0.5_ark*pi) then
             !     sinrho(i) = 1.0_ark
             !     cosrho(i) = 0
             !   else
             !     sinrho(i) = sin(rho)
             !     cosrho(i) = cos(rho)
             !   endif
             !   !
             !enddo
             !
             write(unitfname,"('Numerov basis set # ',i6)") nu_i
             call IOStart(trim(unitfname),io_slot)
             !
             !$omp parallel private(phivphi_t,phil,phir,dphil,dphir,alloc_p,phil_leg,phir_leg,&
             !$omp& phil_sin,phir_sin,dphil_leg,dphir_leg)
             allocate (phivphi_t(0:npoints),phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints),&
                       phil_leg(0:npoints),phir_leg(0:npoints),phil_sin(0:npoints),phir_sin(0:npoints),&
                       dphil_leg(0:npoints),dphir_leg(0:npoints),stat=alloc_p)
             if (alloc_p/=0) then
                write (out,"(' FLbset1DNew error ',i9,' trying to allocate 11 arrays phivphi_t,phi, etc')") alloc_p
                stop 'FLbset1DNew, phivphi_t  - out of memory'
             end if
             !
             !$omp do private(vl,nl,krot1,k_l,vr,nr,krot2,k_r,Tcoeff,iterm,k1,k2,imu,mat_t) schedule(dynamic)
             do vl = 0,bs%Size
                !
                read (io_slot,rec=vl+1) (phil_leg(k),k=0,npoints),(dphil_leg(k),k=0,npoints)
                !
                krot1 = mod(vl,kmax+1)
                nl = (vl-krot1)/(kmax+1)
                !                
                k_l = 0 ; if (krot1>0) k_l = 1
                !
                phil_sin(:) = phil_leg(:)*sinrho(:)**k_l
                phil(:) = sqrt(sinrho(:))*phil_sin(:)
                dphil(:) = cosrho(:)*0.5_ark*phil_leg(:)+sinrho(:)*dphil_leg(:)
                !dphil(:) = cosrho(:)*(k_l+0.5_ark)*phil_leg(:)+sinrho(:)*dphil_leg(:)
                !
                do vr = 0,bs%Size
                    !
                    read (io_slot,rec=vr+1) (phir_leg(k),k=0,npoints),(dphir_leg(k),k=0,npoints)
                    !
                    krot2 = mod(vr,kmax+1)
                    nr = (vr-krot2)/(kmax+1)
                    !
                    k_r = 0 ; if (krot2>0) k_r = 1
                    !
                    !if ( krot2/=krot1 ) cycle
                    !
                    phir_sin(:) = phir_leg(:)*sinrho(:)**k_r
                    phir(:) = sqrt(sinrho(:))*phir_sin(:)
                    dphir(:) = cosrho(:)*0.5_ark*phir_leg(:)+sinrho(:)*dphir_leg(:)
                    !dphir(:) = cosrho(:)*(k_r+0.5_ark)*phir_leg(:)+sinrho(:)*dphir_leg(:)
                    !
                    if (krot1==krot2) then
                      !
                      Tcoeff = trove%poten%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                        !
                        phivphi_t(:) = phil(:)*trove%poten%field(iterm,:)*phir(:)
                        !
                        trove%poten%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                        !
                      enddo
                      !
                      ! vibraional kinetic energy part
                      !
                      do k1 = 1,Nmodes
                         !
                         do k2 = 1,Nmodes
                           !
                           Tcoeff = trove%g_vib(k1,k2)%Ncoeff
                           !
                           do iterm = 1,Tcoeff
                              !
                              if (k1==Nmodes.and.k2==Nmodes) then 
                                 !
                                 phivphi_t(:) =-dphil_leg(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir_leg(:)*sinrho(:)
                                 !
                              elseif (k1==Nmodes) then 
                                 !
                                 !phivphi_t(:) =-dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                                 phivphi_t(:) = -dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir_sin(:)
                                 !
                              elseif (k2==Nmodes) then 
                                 !
                                 !phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                                 phivphi_t(:) = phil_sin(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                              else
                                 !
                                 phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                              endif
                              !
                              trove%g_vib(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                              !
                           enddo
                           !
                        enddo
                        !
                      enddo
                      !
                      Tcoeff = trove%pseudo%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                         !
                         ! The gvib term can be combined with the pseudopotential part.
                         ! It will allow us to spare additional array and i.e. memory
                         ! and it is also necessary in case of the singular solution to do so.
                         !
                         !!phivphi_t(:) =-2.0_ark*phil(:)*trove%pseudo%field(iterm,:)*phir(:)
                         !
                         !!phivphi_t(:) =-2.0_ark*phil_leg(:)*trove%pseudo%field(iterm,:)*phir_leg(:)
                         !
                         phivphi_t(:) =-2.0_ark*phil_sin(:)*trove%pseudo%field(iterm,:)*phir_sin(:)
                         !
                         mat_t = integral_rect_ark(npoints,rho_range,phivphi_t)
                         !
                         trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr)+mat_t
                         !
                      enddo
                      !
                      ! Vibraional Angular Momentum L2
                      !
                      if (FLl2_coeffs) then
                        !
                        stop 'FLl2_coeffs is not implemented for Legendre'
                        !
                        do k1 = 1,Nmodes
                          !
                          do k2 = 1,Nmodes
                             !
                             Tcoeff = trove%L2_vib(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                if (k1==Nmodes.and.k2==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                elseif (k1==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                elseif (k2==Nmodes) then 
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                else
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                endif
                                !
                                trove%L2_vib(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                          enddo
                          !
                        enddo
                        !
                      endif
                      !
                    endif
                    !
                    if (FLrotation) then
                       !
                       ! rotation kinetic energy part
                       !
                       !
                       ! Special basis set sqrt(sin(rho)) L_n
                       !
                       do k1 = 1,3
                         do k2 = 1,3
                           !
                           if (k1==3.and.k2==3) then
                             !
                             if (krot1+krot2==0.or.krot1/=krot2) cycle
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                               !
                               phivphi_t(:) = phil_leg(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_leg(:)*sinrho(:)
                               !
                               !trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                               trove%g_rot(k1,k2)%me(iterm,vl,vr) = 0
                               !
                               mat_t = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                               trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr)-&
                                                                            mat_t*real(krot1**2,ark)
                               !
                             enddo
                             !
                           elseif (k1==3.or.k2==3) then
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                               !
                               phivphi_t(:) =phil_sin(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_sin(:)
                               !
                               trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                             enddo
                             !
                           else
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                phivphi_t(:) = phil(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir(:)
                                !
                                trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                           endif
                           !
                         enddo
                         !
                       enddo
                       !
                       ! Coriolis part
                       !
                       do k1 = 1,Nmodes
                         !
                         do k2 = 1,3
                            !
                            Tcoeff = trove%g_cor(k1,k2)%Ncoeff
                            !
                            do iterm = 1,Tcoeff
                               !
                               if (k1==Nmodes) then 
                                  !
                                  !
                                  phivphi_t(:) = trove%g_cor(k1,k2)%field(iterm,:)*( phil_sin(:)*dphir(:) - dphil(:)*phir_sin(:))
                                  !
                               else
                                  !
                                  phivphi_t(:) = phil(:)*trove%g_cor(k1,k2)%field(iterm,:)*phir(:)
                                  !
                               endif
                               !
                               trove%g_cor(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                               !
                            enddo
                            !
                         enddo
                         !
                       enddo
                       !
                    endif
                    !
                    if (FLextF_coeffs) then
                      !
                      ! External field part
                      !
                      do imu = 1,extF%rank
                        !
                        Tcoeff = trove%extF(imu)%Ncoeff
                        !
                        do iterm = 1,Tcoeff
                           !
                           phivphi_t(:) = phil(:)*trove%extF(imu)%field(iterm,:)*phir(:)
                           !
                           trove%extF(imu)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                           !
                        enddo
                        !
                      enddo
                      !
                    endif 
                    !
                enddo 
                !
             enddo
             !$omp end do
             !
             deallocate (phivphi_t,phil,phir,dphil,dphir,phil_leg,phir_leg,phil_sin,phir_sin,dphil_leg,dphir_leg)
             !$omp end parallel 
             !
             deallocate(sinrho,cosrho,stat=alloc_p)
             !
           elseif (trim(bs%type)=='LAGUERRE-K') then
             !
             if (job%verbose>=4) then 
                write(out,"('   Allocating 11 arrays of ',i8,', Ncore times ',g12.4,' ')") trove%Npoints+1,&
                                11.0_rk*real(trove%Npoints+1,rk)/1024.0_rk**3
             endif
             !
             !do i = 0,npoints
             !   rho =  rho_b(1)+real(i,kind=ark)*trove%rhostep
             !   !
             !   if (rho>0.5_ark*pi) then
             !     sinrho(i) = 1.0_ark
             !     cosrho(i) = 0
             !   else
             !     sinrho(i) = sin(rho)
             !     cosrho(i) = cos(rho)
             !   endif
             !   !
             !enddo
             !
             write(unitfname,"('Numerov basis set # ',i6)") nu_i
             call IOStart(trim(unitfname),io_slot)
             !
             !$omp parallel private(phivphi_t,phil,phir,dphil,dphir,alloc_p,phil_leg,phir_leg,&
             !$omp& phil_sin,phir_sin,dphil_leg,dphir_leg)
             allocate (phivphi_t(0:npoints),phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints),&
                       phil_leg(0:npoints),phir_leg(0:npoints),phil_sin(0:npoints),phir_sin(0:npoints),&
                       dphil_leg(0:npoints),dphir_leg(0:npoints),stat=alloc_p)
             if (alloc_p/=0) then
                write (out,"(' FLbset1DNew error ',i9,' trying to allocate 11 arrays phivphi_t,phi, etc')") alloc_p
                stop 'FLbset1DNew, phivphi_t  - out of memory'
             end if
             !
             !$omp do private(vl,nl,krot1,k_l,vr,nr,krot2,k_r,Tcoeff,iterm,k1,k2,imu,mat_t) schedule(dynamic)
             do vl = 0,bs%Size
                !
                read (io_slot,rec=vl+1) (phil_leg(k),k=0,npoints),(dphil_leg(k),k=0,npoints)
                !
                krot1 = mod(vl,kmax+1)
                nl = (vl-krot1)/(kmax+1)
                !                
                k_l = 0 ; if (krot1>0) k_l = 1
                !
                phil_sin(:) = phil_leg(:)*mrho(:)**k_l
                phil(:) = sqrt(mrho(:))*phil_sin(:)
                !dphil(:) = 0.5_ark*phil_leg(:)+mrho(:)*dphil_leg(:)
                !
                dphil(:) = 0.5_ark*phil_sin(:)+mrho(:)*dphil_leg(:)
                !
                do vr = 0,bs%Size
                    !
                    read (io_slot,rec=vr+1) (phir_leg(k),k=0,npoints),(dphir_leg(k),k=0,npoints)
                    !
                    krot2 = mod(vr,kmax+1)
                    nr = (vr-krot2)/(kmax+1)
                    !
                    if (abs(krot1-krot2)>2) cycle 
                    !
                    k_r = 0 ; if (krot2>0) k_r = 1
                    !
                    phir_sin(:) = phir_leg(:)*mrho**k_r
                    phir(:) = sqrt(mrho(:))*phir_sin(:)
                    !dphir(:) = 0.5_ark*phir_leg(:)+mrho(:)*dphir_leg(:)
                    !
                    dphir(:) = 0.5_ark*phir_sin(:)+mrho(:)*dphir_leg(:)
                    !
                    if (krot1==krot2) then
                      !
                      Tcoeff = trove%poten%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                        !
                        phivphi_t(:) = phil(:)*trove%poten%field(iterm,:)*phir(:)
                        !
                        trove%poten%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                        !
                      enddo
                      !
                      ! vibraional kinetic energy part
                      !
                      do k1 = 1,Nmodes
                         !
                         do k2 = 1,Nmodes
                           !
                           Tcoeff = trove%g_vib(k1,k2)%Ncoeff
                           !
                           do iterm = 1,Tcoeff
                              !
                              if (k1==Nmodes.and.k2==Nmodes) then 
                                 !
                                 phivphi_t(:) =-dphil_leg(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir_leg(:)*mrho(:)
                                 !
                              elseif (k1==Nmodes) then 
                                 !
                                 !phivphi_t(:) =-dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                                 phivphi_t(:) = -dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir_sin(:)
                                 !
                              elseif (k2==Nmodes) then 
                                 !
                                 !phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                                 phivphi_t(:) = phil_sin(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                              else
                                 !
                                 phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                              endif
                              !
                              trove%g_vib(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                              !
                           enddo
                           !
                        enddo
                        !
                      enddo
                      !
                      Tcoeff = trove%pseudo%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                         !
                         ! The gvib term can be combined with the pseudopotential part.
                         ! It will allow us to spare additional array and i.e. memory
                         ! and it is also necessary in case of the singular solution to do so.
                         !
                         phivphi_t(:) =-2.0_ark*phil_sin(:)*trove%pseudo%field(iterm,:)*phir_sin(:)
                         !
                         mat_t = integral_rect_ark(npoints,rho_range,phivphi_t)
                         !
                         trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr)+mat_t
                         !
                      enddo
                      !
                      ! Vibraional Angular Momentum L2
                      !
                      if (FLl2_coeffs) then
                        !
                        stop 'FLl2_coeffs is not implemented for Legendre'
                        !
                        do k1 = 1,Nmodes
                          !
                          do k2 = 1,Nmodes
                             !
                             Tcoeff = trove%L2_vib(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                if (k1==Nmodes.and.k2==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                elseif (k1==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                elseif (k2==Nmodes) then 
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                else
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                endif
                                !
                                trove%L2_vib(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                          enddo
                          !
                        enddo
                        !
                      endif
                      !
                    endif
                    !
                    if (FLrotation) then
                       !
                       ! rotation kinetic energy part
                       !
                       !
                       ! Special basis set sqrt(sin(rho)) L_n
                       !
                       do k1 = 1,3
                         do k2 = 1,3
                           !
                           if (k1==3.and.k2==3) then
                             !
                             if (krot1+krot2==0.or.krot1/=krot2) cycle
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                               !
                               phivphi_t(:) = phil_leg(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_leg(:)*mrho(:)
                               !
                               if (trove%kmax<=trove%krot) then 
                                 !
                                 trove%g_rot(k1,k2)%me(iterm,vl,vr) = 0
                                 !
                                 mat_t = integral_rect_ark(npoints,rho_range,phivphi_t)
                                 !
                                 trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr)-&
                                                                              mat_t*real(krot1**2,ark)
                                 !
                               else
                                 !
                                 trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                 !
                               endif
                               !
                             enddo
                             !
                           elseif (k1==3.or.k2==3) then
                             !
                             if (abs(krot1-krot2)>1) cycle
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                               !
                               phivphi_t(:) =phil_sin(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_sin(:)
                               !
                               !!phivphi_t(:) =phil_leg(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_leg(:)
                               trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                             enddo
                             !
                           else
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                phivphi_t(:) = phil(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir(:)
                                !
                                trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                           endif
                           !
                         enddo
                         !
                       enddo
                       !
                       ! Coriolis part
                       !
                       do k1 = 1,Nmodes
                         !
                         do k2 = 1,3
                            !
                            if (k2==3.and.krot1/=krot2) cycle
                            if (k2/=3.and.krot1==krot2) cycle
                            if (abs(krot1-krot2)>1) cycle
                            !
                            Tcoeff = trove%g_cor(k1,k2)%Ncoeff
                            !
                            do iterm = 1,Tcoeff
                               !
                               if (k1==Nmodes) then 
                                  !
                                  phivphi_t(:) = trove%g_cor(k1,k2)%field(iterm,:)*( phil_sin(:)*dphir(:) - dphil(:)*phir_sin(:))
                                  !
                               else
                                  !
                                  phivphi_t(:) = phil(:)*trove%g_cor(k1,k2)%field(iterm,:)*phir(:)
                                  !
                               endif
                               !
                               trove%g_cor(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                               !
                            enddo
                            !
                         enddo
                         !
                       enddo
                       !
                    endif
                    !
                    if (FLextF_coeffs) then
                      !
                      ! External field part
                      !
                      if (abs(krot1-krot2)>2) cycle
                      !
                      do imu = 1,extF%rank
                        !
                        !if (imu==3.and.krot1/=krot2) cycle
                        !if (imu/=3.and.krot1==krot2) cycle
                        !
                        Tcoeff = trove%extF(imu)%Ncoeff
                        !
                        do iterm = 1,Tcoeff
                           !
                           phivphi_t(:) = phil(:)*trove%extF(imu)%field(iterm,:)*phir(:)
                           !
                           mat_t = integral_rect_ark(npoints,rho_range,phivphi_t)
                           !
                           if (abs(mat_t)>extF%matelem_threshold) then
                             trove%extF(imu)%me(iterm,vl,vr) = mat_t
                           endif
                           !
                        enddo
                        !
                      enddo
                      !
                    endif 
                    !
                enddo 
                !
             enddo
             !$omp end do
             !
             deallocate (phivphi_t,phil,phir,dphil,dphir,phil_leg,phir_leg,phil_sin,phir_sin,dphil_leg,dphir_leg)
             !$omp end parallel 
             !
             deallocate(mrho,stat=alloc_p)
             !
           elseif (trim(bs%type)=='SINRHO-1') then
             !
             if (job%verbose>=4) then 
                write(out,"('   Allocating 11 arrays of ',i8,', Ncore times ',g12.4,' ')") trove%Npoints+1,&
                                11.0_rk*real(trove%Npoints+1,rk)/1024.0_rk**3
             endif
             !
             !$omp parallel private(phivphi_t,phil,phir,dphil,dphir,alloc_p,phil_leg,phir_leg,&
             !$omp& phil_sin,phir_sin,dphil_leg,dphir_leg)
             allocate (phivphi_t(0:npoints),phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints),&
                       phil_leg(0:npoints),phir_leg(0:npoints),phil_sin(0:npoints),phir_sin(0:npoints),&
                       dphil_leg(0:npoints),dphir_leg(0:npoints),stat=alloc_p)
             if (alloc_p/=0) then
                write (out,"(' FLbset1DNew error ',i9,' trying to allocate 11 arrays phivphi_t,phi, etc')") alloc_p
                stop 'FLbset1DNew, phivphi_t  - out of memory'
             end if
             !
             !$omp do private(vl,unitfname,io_slot,nl,krot1,vr,nr,krot2,Tcoeff,iterm,k1,k2,imu,mat_t) schedule(dynamic)
             do vl = 0,bs%Size
                !
                write(unitfname,"('Numerov basis set # ',i6)") nu_i
                call IOStart(trim(unitfname),io_slot)
                !
                read (io_slot,rec=vl+1) (phil_leg(k),k=0,npoints),(dphil_leg(k),k=0,npoints)
                !
                krot1 = mod(vl,kmax+1)
                nl = (vl-krot1)/(kmax+1)
                krot11 = 0 ; if (krot1>0) krot11 = 1
                !
                phil_sin(:) = phil_leg(:)*sinrho(:)**krot1
                phil(:) = sqrt(sinrho(:))*phil_sin(:)
                dphil(:) = cosrho(:)*0.5_ark*phil_leg(:)+sinrho(:)*dphil_leg(:)
                !
                do vr = 0,bs%Size
                    !
                    read (io_slot,rec=vr+1) (phir_leg(k),k=0,npoints),(dphir_leg(k),k=0,npoints)
                    !
                    krot2 = mod(vr,kmax+1)
                    nr = (vr-krot2)/(kmax+1)
                    krot21 = 0 ; if (krot2>0) krot21 = 1
                    !
                    !if ( krot2/=krot1 ) cycle
                    !
                    phir_sin(:) = phir_leg(:)*sinrho**krot2
                    phir(:) = sqrt(sinrho(:))*phir_sin(:)
                    dphir(:) = cosrho(:)*0.5_ark*phir_leg(:)+sinrho(:)*dphir_leg(:)
                    !
                    if (krot1==krot2) then
                      !
                      Tcoeff = trove%poten%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                        !
                        phivphi_t(:) = phil(:)*trove%poten%field(iterm,:)*phir(:)
                        !
                        trove%poten%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                        !
                      enddo
                      !
                      ! vibraional kinetic energy part
                      !
                      do k1 = 1,Nmodes
                         !
                         do k2 = 1,Nmodes
                           !
                           Tcoeff = trove%g_vib(k1,k2)%Ncoeff
                           !
                           do iterm = 1,Tcoeff
                              !
                              if (k1==Nmodes.and.k2==Nmodes) then 
                                 !
                                 phivphi_t(:) =-dphil_leg(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir_leg(:)*sinrho(:)
                                 !
                              elseif (k1==Nmodes) then 
                                 !
                                 !phivphi_t(:) =-dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                                 phivphi_t(:) = -dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir_sin(:)
                                 !
                              elseif (k2==Nmodes) then 
                                 !
                                 !phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                                 phivphi_t(:) = phil_sin(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                              else
                                 !
                                 phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                              endif
                              !
                              trove%g_vib(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                              !
                           enddo
                           !
                        enddo
                        !
                      enddo
                      !
                      Tcoeff = trove%pseudo%Ncoeff
                      !
                      do iterm = 1,Tcoeff
                         !
                         ! The gvib term can be combined with the pseudopotential part.
                         ! It will allow us to spare additional array and i.e. memory
                         ! and it is also necessary in case of the singular solution to do so.
                         !
                         !!phivphi_t(:) =-2.0_ark*phil(:)*trove%pseudo%field(iterm,:)*phir(:)
                         !
                         !!phivphi_t(:) =-2.0_ark*phil_leg(:)*trove%pseudo%field(iterm,:)*phir_leg(:)
                         !
                         phivphi_t(:) =-2.0_ark*phil_sin(:)*trove%pseudo%field(iterm,:)*phir_sin(:)
                         !
                         mat_t = integral_rect_ark(npoints,rho_range,phivphi_t)
                         !
                         trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr)+mat_t
                         !
                      enddo
                      !
                      ! Vibraional Angular Momentum L2
                      !
                      if (FLl2_coeffs) then
                        !
                        stop 'FLl2_coeffs is not implemented for Legendre'
                        !
                        do k1 = 1,Nmodes
                          !
                          do k2 = 1,Nmodes
                             !
                             Tcoeff = trove%L2_vib(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                if (k1==Nmodes.and.k2==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                elseif (k1==Nmodes) then 
                                   !
                                   phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                elseif (k2==Nmodes) then 
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                   !
                                else
                                   !
                                   phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                   !
                                endif
                                !
                                trove%L2_vib(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                          enddo
                          !
                        enddo
                        !
                      endif
                      !
                    endif
                    !
                    if (FLrotation) then
                       !
                       ! rotation kinetic energy part
                       !
                       !
                       ! Special basis set sqrt(sin(rho)) L_n
                       !
                       do k1 = 1,3
                         do k2 = 1,3
                           !
                           if (k1==3.and.k2==3) then
                             !
                             if (krot1+krot2==0.or.krot1/=krot2) cycle
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                               !
                               phivphi_t(:) = phil_leg(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_leg(:)*sinrho(:)**(krot1+krot2-1)
                               !
                               !!phivphi_t(:) = -phil_leg(:)*fl%field(iterm,:)*phir_leg(:)
                               !
                               trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                             enddo
                             !
                           elseif (k1==3.or.k2==3) then
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                               !
                               phivphi_t(:) =phil_sin(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_sin(:)
                               !
                               !!phivphi_t(:) =phil_leg(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir_leg(:)
                               trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                             enddo
                             !
                           else
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                phivphi_t(:) = phil(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir(:)
                                !
                                trove%g_rot(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                           endif
                           !
                         enddo
                         !
                       enddo
                       !
                       ! Coriolis part
                       !
                       do k1 = 1,Nmodes
                         !
                         do k2 = 1,3
                            !
                            Tcoeff = trove%g_cor(k1,k2)%Ncoeff
                            !
                            do iterm = 1,Tcoeff
                               !
                               if (k1==Nmodes) then 
                                  !
                                  phivphi_t(:) = trove%g_cor(k1,k2)%field(iterm,:)*( phil_sin(:)*dphir(:) - dphil(:)*phir_sin(:))
                                  !
                               else
                                  !
                                  phivphi_t(:) = phil(:)*trove%g_cor(k1,k2)%field(iterm,:)*phir(:)
                                  !
                               endif
                               !
                               trove%g_cor(k1,k2)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                               !
                               !
                            enddo
                            !
                         enddo
                         !
                       enddo
                       !
                    endif
                    !
                    if (FLextF_coeffs) then
                      !
                      ! External field part
                      !
                      do imu = 1,extF%rank
                        !
                        Tcoeff = trove%extF(imu)%Ncoeff
                        !
                        do iterm = 1,Tcoeff
                           !
                           phivphi_t(:) = phil(:)*trove%extF(imu)%field(iterm,:)*phir(:)
                           !
                           trove%extF(imu)%me(iterm,vl,vr) = integral_rect_ark(npoints,rho_range,phivphi_t)
                           !
                        enddo
                        !
                      enddo
                      !
                    endif 
                    !
                enddo 
                !
             enddo
             !$omp end do
             !
             deallocate (phivphi_t,phil,phir,dphil,dphir,phil_leg,phir_leg,phil_sin,phir_sin,dphil_leg,dphir_leg)
             !$omp end parallel 
             !
             deallocate(sinrho,cosrho,stat=alloc_p)
             !
           else
             !
             !$omp parallel private(phivphi_t,phil,phir,dphil,dphir,alloc_p)
             allocate (phivphi_t(0:npoints),phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints),stat=alloc_p)
             if (alloc_p/=0) then
                write (out,"(' FLbset1DNew error ',i9,' trying to allocate 11 arrays phivphi_t,phi, etc')") alloc_p
                stop 'FLbset1DNew, phivphi_t  - out of memory'
             end if
             !
             !$omp do private(vl,unitfname,io_slot,i,rho,k,vr,Tcoeff,iterm,k1,k2,imu,mat_t) schedule(dynamic)
             do vl = 0,bs%Size
                !
                write(unitfname,"('Numerov basis set # ',i6)") nu_i
                call IOStart(trim(unitfname),io_slot)
                !
                read (io_slot,rec=vl+1) (phil(k),k=0,npoints),(dphil(k),k=0,npoints)
                !
                do vr = 0,bs%Size
                    !
                    !if (vl==vr) then
                    !    phir =  phil
                    !   dphir = dphil
                    !else
                    read (io_slot,rec=vr+1) (phir(k),k=0,npoints),(dphir(k),k=0,npoints)
                    !endif
                    !
                    !
                    Tcoeff = trove%poten%Ncoeff
                    !
                    do iterm = 1,Tcoeff
                      !
                      phivphi_t(:) = phil(:)*trove%poten%field(iterm,:)*phir(:)
                      !
                      trove%poten%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                      !
                    enddo
                    !
                    ! vibraional kinetic energy part
                    !
                    do k1 = 1,Nmodes
                       !
                       do k2 = 1,Nmodes
                         !
                         Tcoeff = trove%g_vib(k1,k2)%Ncoeff
                         !
                         do iterm = 1,Tcoeff
                            !
                            if (k1==Nmodes.and.k2==Nmodes) then 
                               !
                               phivphi_t(:) =-dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                               !
                            elseif (k1==Nmodes) then 
                               !
                               phivphi_t(:) =-dphil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                               !
                            elseif (k2==Nmodes) then 
                               !
                               phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*dphir(:)
                               !
                            else
                               !
                               phivphi_t(:) = phil(:)*trove%g_vib(k1,k2)%field(iterm,:)*phir(:)
                               !
                            endif
                            !
                            trove%g_vib(k1,k2)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                            !
                         enddo
                         !
                      enddo
                      !
                    enddo
                    !
                    !
                    Tcoeff = trove%pseudo%Ncoeff
                    !
                    do iterm = 1,Tcoeff
                       !
                       ! The gvib term can be combined with the pseudopotential part.
                       ! It will allow us to spare additional array and i.e. memory
                       ! and it is also necessary in case of the singular solution to do so.
                       !
                       phivphi_t(:) =-2.0_ark*phil(:)*trove%pseudo%field(iterm,:)*phir(:)
                       !
                       mat_t = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                       !
                       trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(iterm,vl,vr)+mat_t
                       !
                    enddo
                    !
                    ! Vibraional Angular Momentum L2
                    !
                    if (FLl2_coeffs) then
                      !
                      do k1 = 1,Nmodes
                        !
                        do k2 = 1,Nmodes
                           !
                           Tcoeff = trove%L2_vib(k1,k2)%Ncoeff
                           !
                           do iterm = 1,Tcoeff
                              !
                              if (k1==Nmodes.and.k2==Nmodes) then 
                                 !
                                 phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                              elseif (k1==Nmodes) then 
                                 !
                                 phivphi_t(:) =-dphil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                              elseif (k2==Nmodes) then 
                                 !
                                 phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*dphir(:)
                                 !
                              else
                                 !
                                 phivphi_t(:) = phil(:)*trove%L2_vib(k1,k2)%field(iterm,:)*phir(:)
                                 !
                              endif
                              !
                              trove%L2_vib(k1,k2)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                              !
                           enddo
                           !
                        enddo
                        !
                      enddo
                      !
                    endif
                    !
                    if (FLrotation) then
                       !
                       ! rotation kinetic energy part
                       !
                       do k1 = 1,3
                          !
                          do k2 = 1,3
                             !
                             Tcoeff = trove%g_rot(k1,k2)%Ncoeff
                             !
                             do iterm = 1,Tcoeff
                                !
                                phivphi_t(:) = phil(:)*trove%g_rot(k1,k2)%field(iterm,:)*phir(:)
                                !
                                trove%g_rot(k1,k2)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                                !
                             enddo
                             !
                          enddo
                          !
                       enddo
                       !
                       !
                       ! Coriolis part
                       !
                       do k1 = 1,Nmodes
                         !
                         do k2 = 1,3
                            !
                            Tcoeff = trove%g_cor(k1,k2)%Ncoeff
                            !
                            do iterm = 1,Tcoeff
                               !
                               if (k1==Nmodes) then 
                                  !
                                  phivphi_t(:) = trove%g_cor(k1,k2)%field(iterm,:)*( phil(:)*dphir(:) - dphil(:)*phir(:))
                                  !
                               else
                                  !
                                  phivphi_t(:) = phil(:)*trove%g_cor(k1,k2)%field(iterm,:)*phir(:)
                                  !
                               endif
                               !
                               trove%g_cor(k1,k2)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                               !
                               !
                            enddo
                            !
                         enddo
                         !
                       enddo
                       !
                    endif
                    !
                    if (FLextF_coeffs) then
                      !
                      ! External field part
                      !
                      do imu = 1,extF%rank
                        !
                        Tcoeff = trove%extF(imu)%Ncoeff
                        !
                        do iterm = 1,Tcoeff
                           !
                           phivphi_t(:) = phil(:)*trove%extF(imu)%field(iterm,:)*phir(:)
                           !
                           trove%extF(imu)%me(iterm,vl,vr) = simpsonintegral_ark(npoints,rho_range,phivphi_t)
                           !
                        enddo
                        !
                      enddo
                      !
                    endif 
                    !
                enddo 
                !
             enddo
             !$omp end do
             !
             deallocate (phivphi_t,phil,phir,dphil,dphir)
             !$omp end parallel 
             !
           endif
           !
           call TimerStop('Primitive matrix elements')
           !
           if (job%verbose>=3) call TimerReport
           !
           deallocate (f1drho,g1drho,phivphi,weight)
           !
        else
           !
           ! numerov bset
           !
           ! for basis sets we will need 1d potential and kinetic enrgy part 
           ! in terms of the corresponding coordinate 
           !
           select case(manifold)
             !
           case (0) ! rigid case, expansion around the rank=0 manifold
             !
             irho_eq = 0 
             !
           case (1) ! non-rigid case, expansion around a rank=1 manifold
             !
             irho_eq = mod(nint( ( molec%chi_eq(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
             !
             ! altenatively try to put equilibrimum at the minimum of the potential
             !
             irho_eq = minloc(trove%poten%field(1,:),dim=1)-1
             !
             ! for periodic case we need to choose the equilibrium at the 
             ! reference geometry; we apply this rule only if the periodicity is the same 
             ! as the size of the class 
             ! 
             if (trove%periodic) then 
               !.and.bs%imodes   ==job%bset(Nmodes)%iperiod) then
               period = (job%bset(Nmodes)%borders(2)-job%bset(Nmodes)%borders(1))/real(job%bset(Nmodes)%iperiod,ark)
               periodic_model = .true.
             endif
             !
           end select
           !
           if (irho_eq<0.or.irho_eq>trove%Npoints) then
              write(out,"('FLbset1DNew-numerov: equilibrium is ill defined:')")
              write(out,"('irho_eq is ',i6)") irho_eq
              stop 'FLbset1DNew:  equilibrium is ill defined'
           endif
           !
           nu_i = bs%mode(1)
           !
           npoints = bset%dscr(nu_i)%npoints
           rho_b(1:2)=bset%dscr(nu_i)%borders(1:2)
           step = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
           !rho_ref = molec%chi_eq(nu_i)
           !
           allocate (f1drho(0:Npoints),g1drho(0:Npoints),drho(0:Npoints,3),xton(0:Npoints,0:trove%MaxOrder),stat=alloc)
           if (alloc/=0) then
              write (out,"(' Error ',i9,' trying to allocate f1drho and g1drho')") alloc
              stop 'FLbset1DNew, f1drho and g1drho - out of memory'
           end if
           !
           xton = 0
           !
           ! Defining 1D potential and kinetic energy functions
           !
           reduced_model = .false.
           !
           if (bset%dscr(nu_i)%model<trove%NPotOrder) then 
             !
             reduced_model = .true.
             !
           endif 
           !
           if (.not.trove%DVR.or.reduced_model) then
             !
             f1d = 0 
             !
             do ipower = 0,min(trove%NPotOrder,max(bset%dscr(nu_i)%model,2))
                !
                f2 = 0 
                !
                do imode = 1,bs%imodes
                  !
                  nu_i = bs%mode(imode)
                  powers = 0 ; powers(nu_i) = ipower
                  k = FLQindex(trove%Nmodes_e,powers)
                  !
                  ! shift the minimum by the period of the last mode if present
                  if (periodic_model) then 
                    !
                    rho_ref_ = trove%rho_ref+period*real(imode-1,ark)
                    irho_eq = mod(nint( ( rho_ref_-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
                    !
                  endif
                  !
                  if (trove%sparse) then
                    !
                    call find_isparse_from_ifull(trove%poten%Ncoeff,trove%poten%ifromsparse,k,i)
                    !
                    f2(imode) = 0 
                    if (i/=0) f2(imode) = trove%poten%field(i,irho_eq)
                    !
                    if (ipower==2.and.abs(f2(imode))<sqrt(small_)) then
                      write(out,"('FLbset1DNew: f2=0 in the poten sparse-field for imode = ',i5,' ipower',i5)") imode,ipower
                      stop 'FLbset1DNew: f2=0 in the poten sparse-field'
                    endif
                    !
                    if (periodic_model) then 
                      !
                      do jmode = 1,bs%imodes
                        !
                        rho_ref_ = trove%rho_ref+period*real(jmode-1,ark)
                        irho_eq = mod(nint( ( rho_ref_-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
                        call find_isparse_from_ifull(trove%poten%Ncoeff,trove%poten%ifromsparse,k,i)
                        !
                        f_period(imode,jmode) = 0 
                        if (i/=0) f_period(imode,jmode)= trove%poten%field(i,irho_eq)
                        !
                      enddo
                      !
                    endif
                    !
                  else
                    !
                    f2(imode) = trove%poten%field(k,irho_eq)
                    !
                  endif
                  !
                enddo
                !
                ! Check if all potential parameters are equal for every considered mode
                !
                f_t= f2(1) 
                if (any(abs(f2(1:bs%imodes)-f_t)>1e-5*abs(f_t)).and.&
                    any(abs(f2(1:bs%imodes)-f_t)>1e6*sqrt(small_))) then
                   write(out,"('FLbset1DNew: not all numerov-pot parameters are equal')")
                   write(out,"('pot-ipower=',i6)") ipower
                   write(out,"(30f18.8)") (f2(imode),imode=1,min(bs%imodes,30)) 
                   if ( job%bset(nu_i)%check_sym ) then 
                     stop 'FLbset1DNew: not all numerov parameters are equal'
                   endif
                endif
                !
                f1d(ipower) = f2(1)
                !
             enddo 
             !
             p1d = 0 
             !
             do ipower = 0,min(trove%NKinOrder,max(bset%dscr(nu_i)%model-2,0))
                !
                f2 = 0 
                !
                do imode = 1,bs%imodes
                  !
                  nu_i = bs%mode(imode)
                  powers = 0 ; powers(nu_i) = ipower
                  k = FLQindex(trove%Nmodes_e,powers)
                  !
                  ! shift the minimum by the period of the last mode if present
                  if (periodic_model) then 
                    !
                    rho_ref_ = trove%rho_ref+period*real(imode-1,ark)
                    irho_eq = mod(nint( ( rho_ref_-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
                    !
                  endif
                  !
                  if (trove%sparse) then
                    !
                    call find_isparse_from_ifull(trove%pseudo%Ncoeff,trove%pseudo%ifromsparse,k,i)
                    !
                    f2(imode) = 0
                    if (i/=0) f2(imode) = trove%pseudo%field(i,irho_eq)
                    !
                    !if (ipower==0.and.abs(f2(imode))<sqrt(small_)) then
                    !  write(out,"('FLbset1DNew: f2=0 in the pseudo sparse-field for imode = ',i5,' ipower',i5)") imode,ipower
                    !  stop 'FLbset1DNew: f2=0 in the pseudo sparse-field'
                    !endif
                    !
                  else
                    !
                    f2(imode) = trove%pseudo%field(k,irho_eq)
                    !
                  endif
                  !
                enddo
                !
                ! Check if all potential parameters are equal for every considered mode
                !
                f_t= f2(1) 
                if (any(abs(f2(1:bs%imodes)-f_t)>1e-5*abs(f_t)).and.&
                    any(abs(f2(1:bs%imodes)-f_t)>1e6*sqrt(small_))) then
                   write(out,"('FLbset1DNew: not all numerov-pseudo parameters are equal')")
                   write(out,"('pseudo-ipower=',i6)") ipower
                   write(out,"(30f18.8)") (f2(imode),imode=1,min(bs%imodes,30)) 
                   if ( job%bset(nu_i)%check_sym ) then 
                     stop 'FLbset1DNew: not all numerov parameters are equal'
                   endif
                endif
                !
                p1d(ipower) = f2(1)
                !
             enddo 
             !
             ! for the kinetic part it is simpler - we just take the corresoinding diagonal member of the g_vib%coeffs
             !
             g1d = 0 
             !
             do ipower = 0,min(trove%NKinOrder,max(bset%dscr(nu_i)%model-2,0))
                !
                !
                do imode = 1,bs%imodes
                  !             
                  nu_i = bs%mode(imode)
                  fl => trove%g_vib(nu_i,nu_i) 
                  !
                  select case(job%bset(imode)%coord_kinet)
                    !
                    case('BOND-LENGTH', 'ANGLE', 'DIHEDRAL')
                      !
                      g2(imode) = 0
                      !
                      do j = 1, size(fl%ifromsparse)
                        !
                        g2_term = 0
                        !write(*,*) fl%ifromsparse(j)
                        !
                        powers = powers_from_index(Nmodes, fl%ifromsparse(j))
                        !
                        !write(*,*) "powers ", powers, " end"
                        !
                        if (powers(nu_i)/= ipower) cycle
                        g2_term  =  fl%field(j,irho_eq)
                        !
                        do i = 1, size(powers)
                          if(i == nu_i) cycle
                          !write(*,*) "eq ", trove%chi_eq(nu_i)   
                          !write(*,*) "powers(I) ", powers(i)
                          !write(*,*) "MLcoord", MLcoord_direct(trove%chi_eq(nu_i), 1, i, powers(i))
                          !
                          g2_term = g2_term*MLcoord_direct(trove%chi_eq(nu_i), 1, i, powers(i)) 
                          !
                        enddo
                        !
                        g2(imode) = g2(imode) + g2_term
                        !write(*,*) "nu_i ", nu_i, " g2(imode) ", g2(imode) 
                        !
                     enddo
                     !    
                  case default 
                     !
                     powers = 0 ; powers(nu_i) = ipower
                     k = FLQindex(trove%Nmodes_e,powers)
                     !
                     ! shift the minimum by the period of the last mode if present
                     if (periodic_model) then 
                       !
                       rho_ref_ = trove%rho_ref+period*real(imode-1,ark)
                       irho_eq = mod(nint( ( rho_ref_-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
                       !
                     endif
                     !
                     if (trove%sparse) then
                       !
                       call find_isparse_from_ifull(fl%Ncoeff,fl%ifromsparse,k,i)
                       !
                       g2(imode) = 0
                       !
                       if (i/=0) g2(imode) = fl%field(i,irho_eq)
                       !
                       if (ipower==0.and.abs(g2(imode))<sqrt(small_)) then
                         write(out,"('FLbset1DNew: g2=0 in the gvib sparse-field for imode = ',i5,' ipower',i5)") imode,ipower
                         stop 'FLbset1DNew: g2=0 in the gvib sparse-field'
                       endif
                       !
                     else
                       !
                       g2(imode) = fl%field(k,irho_eq)
                       !
                     endif
                     !
                  end select 
                  !
                enddo
                ! 
                !
                ! The same check for all g2 kinetic parameters to be equal for every considered mode
                !
                f_t= g2(1) 
                if (any(abs(g2(1:bs%imodes)-f_t)>100000.0_rk*sqrt(small_)*abs(f_t)).and.abs(f_t)>1000.0_rk*sqrt(small_)) then
                   write(out,"('FLbset1DNew-numerov: not all gvib-zero-order kinetic parameters are equal')")
                   write(out,"('gvib-ipower=',i6)") ipower
                   write(out,"(30f18.8)") (g2(imode),imode=1,min(bs%imodes,30)) 
                   if ( job%bset(nu_i)%check_sym ) then 
                     stop 'FLbset1DNew: not all zero-order kinetic parameters are equal'
                   endif
                endif
                !
                g1d(ipower) = g2(1)
                !
             enddo
             !
             !
             if (trim(molec%coords_transform)=='R-RHO'.and..false.) then 
                !
                rho_t = sum( trove%mass(:)*( trove%b0(:,2,irho_eq)**2+ trove%b0(:,3,irho_eq)**2 ) )/&
                                            (planck*avogno*real(1.0d+16,kind=rk)/(4.0_ark*pi*pi*vellgt))
                !
                f1d(:) = f1d(:)*rho_t
                !
             endif
             !
             ! Check if defined f2 and g2 parameters are not zero 
             !
             f_t = maxval(f1d,dim=1)
             !
             if (all(abs(f1d(:))<1e-5*abs(f_t)) ) then
                write(out,"('FLbset1DNew: all f1d  are zero ',f18.8)") f1d(:)
                stop 'FLbset1DNew: f1d are all zero'
             endif
             !
             f_t = maxval(g1d,dim=1)
             !
             if (all(abs(g1d(:))<1e-5*abs(f_t)) ) then
                write(out,"('FLbset1DNew: all g1d  are zero ',f18.8)") g1d(:)
                stop 'FLbset1DNew: g1d are all zero'
             endif
             !
             ! generating 1d tables
             ! for basis sets we will need 1d potential and kinetic enrgy part 
             ! in terms of the corresponding coordinate 
             !
             f1drho = 0 ; g1drho = 0 
             !
             nu_i = bs%mode(1)
             !
             rho =  trove%chi_eq(nu_i)
             !
             rho_kin0 = MLcoord_direct(rho,1,nu_i)
             rho_pot0 = MLcoord_direct(rho,2,nu_i)
             rho_ext0 = MLcoord_direct(rho,3,nu_i)
             !
             do i = 0,npoints
                !
                rho =  rho_b(1)+real(i,kind=ark)*step
                !
                rho_kin = MLcoord_direct(rho,1,nu_i)-rho_kin0
                rho_pot = MLcoord_direct(rho,2,nu_i)-rho_pot0
                rho_ext = MLcoord_direct(rho,3,nu_i)-rho_ext0
                !
                drho(i,1) = rho_kin
                drho(i,2) = rho_pot
                drho(i,3) = rho_ext
                !
                do ipower = 0,trove%NPotorder 
                   !
                   f1drho(i) = f1drho(i) + f1d(ipower)*rho_pot**ipower
                   !
                enddo 
                !
                do ipower = 0, trove%NKinorder 
                   xton(i,ipower) = MLcoord_direct(rho,1,nu_i,ipower)
                enddo
                !
                do ipower = 0,trove%NKinorder
                   !
                   !f1drho(i) = f1drho(i) + p1d(ipower)*rho_kin**ipower
                   !g1drho(i) = g1drho(i) + g1d(ipower)*rho_kin**ipower
                   !
                   f1drho(i) = f1drho(i) + p1d(ipower)*xton(i,ipower)
                   g1drho(i) = g1drho(i) + g1d(ipower)*xton(i,ipower)
                   !
                enddo 
                !
             enddo
             !
           else 
             !
             ! DVR 
             !
             reduced_model = .false.
             !
             chi = 0
             !
             nu_i = bs%mode(1)
             !
             rho =  trove%chi0_ref(nu_i)
             !
             rho_kin0 = MLcoord_direct(rho,1,nu_i)
             rho_pot0 = MLcoord_direct(rho,2,nu_i)
             rho_ext0 = MLcoord_direct(rho,3,nu_i)
             !
             do i = 0,Npoints
               !
               rho =  rho_b(1)+real(i,kind=ark)*step
               !
               !chi(nu_i) = rho_b(1)+real(i,kind=rk)*step-trove%chi0_ref(nu_i)
               !
               chi(nu_i) = rho_b(1)+real(i,kind=ark)*step-trove%chi_ref(nu_i,irho_eq)
               !
               call FLcalc_poten_kinet_dvr(chi,irho_eq,poten_t,gvib_t,grot_t,gcor_t,extF_t,reduced_model)
               !
               !if (any(bs%mode(:)==Nmodes).and.molec%meptype/='') then 
               !  !
               !  !ar_t(:) = ML_MEPfunc(molec%ncoords,rho)
               !  !
               !  dir = 'INVERSE'
               !  ar_t(:) = MLcoordinate_transform_func(chi,trove%Ncoords,dir)
               !  !
               !  poten_t = MLpotenfunc(ar_t)
               !  !
               !  !chi(1:3) = 1.173629909*cos(rho)-.2563444694*cos(rho)**2-.9172854395+.4764805833*sin(rho)-.1040729802*sin(rho)*cos(rho)
               !  !chi(4:5) = 0
               !  !chi(6) = -.6533902655*cos(rho)+.1427136269*cos(rho)**2+.5106766386+1.609277110*sin(rho)-.3514986145*sin(rho)*cos(rho)
               !  !
               !  !call FLcalc_poten_kinet_dvr(chi,0_ik,poten_t,gvib_t,grot_t,gcor_t,extF_t,reduced_model)
               !  !
               !endif 
               !
               g1drho(i) = gvib_t(nu_i,nu_i)
               f1drho(i) = poten_t
               !
               rho_kin = MLcoord_direct(rho,1,nu_i)-rho_kin0
               rho_pot = MLcoord_direct(rho,2,nu_i)-rho_pot0
               rho_ext = MLcoord_direct(rho,3,nu_i)-rho_ext0
               !
               drho(i,1) = rho_kin
               drho(i,2) = rho_pot
               drho(i,3) = rho_ext
               !
             enddo
             !
           endif
           !
           bs%params    = 0
           bs%params(1) = 0
           bs%params(2) = 0
           !
           ! We have stored the numeber of points, rhomax, and rhomin as optional parameters of  "bset%dscr"
           ! now we need them:
           !
           if (trim(bs%type)=='NUMEROV') then
             !
             call ME_numerov(bs%Size,bs%order,rho_b,isingular,npoints,npoints,drho,xton,f1drho,g1drho,nu_i,&
                             job%bset(nu_i)%iperiod,job%verbose,bs%matelements,bs%ener0)
             !
           elseif (trim(bs%type)=='BOX') then
             !
             call ME_box(bs%Size,bs%order,rho_b,isingular,npoints,drho,f1drho,g1drho,nu_i,job%bset(nu_i)%periodic,job%verbose,&
                         bs%matelements,bs%ener0)
             !
           elseif (trim(bs%type)=='FOURIER') then
             !
             call ME_fourier(bs%Size,bs%order,rho_b,isingular,npoints,numerpoints,drho,f1drho,g1drho,nu_i,&
                             job%bset(nu_i)%iperiod,job%verbose,bs%matelements,bs%ener0)
             !
           else
             !
             write(out,"('FLbset1DNew: illegal type for RIGID case ',a)") trim(bs%type)
             stop 'FLbset1DNew: illegal basis set type for RIGID'
            !
           endif
           !
           deallocate (f1drho,g1drho,drho,xton)
           !
        endif 
        !
     end select   
     !
     !
     ! For the rigid molecule (manifold=0) the potential and kinetic expansion are included 
     ! in the last-mode matrix elements, e.g.:
     ! <v|F_{ijk}( xi_{nmode}^i )|v'> =  force_{i,jk} <v| xi_{nmode}^i |v'>
     ! we thus integrate both the potential and kinetic functions over the last-mode coordinate
     !
     !
     if (job%verbose>=3) write(out,"(/'Primitive matrix elements...')")
     !
     if (any(bs%mode(:)==trove%Nmodes).and.manifold==0) then 
        !
        imode = trove%Nmodes
        !
        ! Potential energy part
        !
        fl => trove%poten
        !
        do icoeff = 1,fl%Ncoeff
           !
           ipower = trove%poten%IndexQ(imode,icoeff)
           !
           do vl = 0,bs%Size
              !
              do vr = 0,bs%Size
                  !
                  trove%poten%me(icoeff,vl,vr) = bs%matelements(0,ipower,vl,vr)*fl%field(icoeff,0)
                  !
               enddo
           enddo
        enddo
        !
        ! vibraional kinetic energy part 
        !
        do k1 = 1,Nmodes
           !
           do k2 = 1,Nmodes
              !
              fl => trove%g_vib(k1,k2)
              !
              do icoeff = 1,fl%Ncoeff
                 !
                 ipower = fl%IndexQ(imode,icoeff)
                 !
                 do vl = 0,bs%Size
                    !
                    do vr = 0,bs%Size
                       !
                       if (k1==imode.and.k2==imode) then 
                          !
                          ! this term will be combined with the pseudo-potential one. 
                          !
                          fl%me(icoeff,vl,vr) = bs%matelements( 2,ipower,vl,vr)*fl%field(icoeff,0)
                                      !-2.0_ark*bs%matelements(-1,ipower,vl,vr)*trove%pseudo%field(icoeff,0)
                          !
                       elseif (k1==imode) then 
                          !
                          fl%me(icoeff,vl,vr) = -bs%matelements(1,ipower,vr,vl)*fl%field(icoeff,0)
                          !
                       elseif (k2==imode) then 
                          !
                          fl%me(icoeff,vl,vr) =  bs%matelements(1,ipower,vl,vr)*fl%field(icoeff,0)
                          !
                       else 
                          !
                          fl%me(icoeff,vl,vr) = bs%matelements(-1,ipower,vl,vr)*fl%field(icoeff,0)
                          !
                       endif
                       !
                    enddo
                    !
                 enddo
                 !
              enddo
              !
           enddo
           !
        enddo
        !
        fl => trove%pseudo
        !
        do icoeff = 1,fl%Ncoeff
           !
           ipower = fl%IndexQ(imode,icoeff)
           !
           if (ipower/=trove%g_vib(imode,imode)%IndexQ(imode,icoeff)) then
             write(out,"('FLbset1DNew: ipower-s must be the same for pseudo and gvib(N,N): ',2i6)") &
                          ipower,trove%g_vib(imode,imode)%IndexQ(imode,icoeff)
             stop 'FLbset1DNew: ipower-s are no the same for pseudo and gvib(N,N)'
           endif
           !
           do vl = 0,bs%Size
              !
              do vr = 0,bs%Size
                    !
                    trove%g_vib(Nmodes,Nmodes)%me(icoeff,vl,vr) = trove%g_vib(Nmodes,Nmodes)%me(icoeff,vl,vr)&
                        -2.0_ark*bs%matelements(-1,ipower,vl,vr)*trove%pseudo%field(icoeff,0)
              !
              enddo
              !
           enddo
           !
        enddo
        !
        ! vibraional angular momentum L2
        !
        if (FLl2_coeffs) then
           !
           do k1 = 1,Nmodes
              !
              do k2 = 1,Nmodes
                 !
                 fl => trove%L2_vib(k1,k2)
                 !
                 do icoeff = 1,fl%Ncoeff
                    !
                    ipower = fl%IndexQ(imode,icoeff)
                    !
                    do vl = 0,bs%Size
                       !
                       do vr = 0,bs%Size
                         !
                         if (k1==imode.and.k2==imode) then 
                            !
                            fl%me(icoeff,vl,vr) = bs%matelements( 2,ipower,vl,vr)*fl%field(icoeff,0)
                            !
                         elseif (k1==imode) then 
                            !
                            fl%me(icoeff,vl,vr) = -bs%matelements(1,ipower,vr,vl)*fl%field(icoeff,0)
                            !
                         elseif (k2==imode) then 
                            !
                            fl%me(icoeff,vl,vr) =  bs%matelements(1,ipower,vl,vr)*fl%field(icoeff,0)
                            !
                         else 
                            !
                            fl%me(icoeff,vl,vr) = bs%matelements(-1,ipower,vl,vr)*fl%field(icoeff,0)
                            !
                         endif
                         !
                       enddo
                       !
                    enddo
                    !
                 enddo
                 !
              enddo
              !
           enddo
           !
        endif
        !
        if (FLrotation) then
           !
           ! rotation kinetic energy part
           !
           do k1 = 1,3
              !
              do k2 = 1,3
                 !
                 fl => trove%g_rot(k1,k2)
                 !
                 do icoeff = 1,fl%Ncoeff
                    !
                    ipower = fl%IndexQ(imode,icoeff)
                    !
                    do vl = 0,bs%Size
                       !
                       do vr = 0,bs%Size
                         !
                         fl%me(icoeff,vl,vr) =  bs%matelements(-1,ipower,vl,vr)*fl%field(icoeff,0)
                         !
                       enddo
                       !
                    enddo
                    !
                 enddo
                 !
              enddo
              !
           enddo
           !
           ! Coriolis part
           !
           do k1 = 1,Nmodes
              !
              do k2 = 1,3
                 !
                 fl => trove%g_cor(k1,k2)
                 !
                 do icoeff = 1,fl%Ncoeff
                    !
                    ipower = fl%IndexQ(imode,icoeff)
                    !
                    do vl = 0,bs%Size
                       !
                       do vr = 0,bs%Size
                         !
                         if (k1==Nmodes) then 
                            !
                            mat_t =  bs%matelements(1,ipower,vl,vr) -  bs%matelements(1,ipower,vr,vl)
                            !
                         else
                            !
                            mat_t = bs%matelements(-1,ipower,vl,vr)
                            !
                         endif
                         !
                         fl%me(icoeff,vl,vr) =  mat_t*fl%field(icoeff,0)
                         !
                       enddo
                       !
                    enddo
                    !
                 enddo
                 !
              enddo
              !
           enddo
           !
           !
        endif
        !
        if (FLextF_coeffs) then
          !
          ! external field part
          !
          do imu = 1,extF%rank
             !
             fl => trove%extF(imu)
             !
             do icoeff = 1,fl%Ncoeff
                !
                ipower = fl%IndexQ(imode,icoeff)
                !
                do vl = 0,bs%Size
                   !
                   do vr = 0,bs%Size
                     !
                     if (abs(bs%matelements(3,ipower,vl,vr)*fl%field(icoeff,0))>extF%matelem_threshold) then
                       fl%me(icoeff,vl,vr) =  bs%matelements(3,ipower,vl,vr)*fl%field(icoeff,0)
                     endif
                     !
                     !fl%me(icoeff,vl,vr) =  bs%matelements(3,ipower,vl,vr)*fl%field(icoeff,0)
                     !
                   enddo
                   !
                enddo
                !
             enddo
             !
          enddo
          !
        endif 
        !
     endif 
     !
     !  Be verbose!
     !
     call Printbset1DInfo (ibs)
     !
     if (job%verbose>=4) write(out,"('FLbset1DNew/end')") 
     !
     contains 
     !
     subroutine find_isparse_from_ifull(Nterms,ifromsparse,ifull,isparse)
       implicit none
       !
       integer(ik),intent(in)  :: Nterms
       integer(ik),intent(in)  :: ifromsparse(Nterms)
       integer(ik),intent(in)  :: ifull
       integer(ik),intent(out) :: isparse
       !
       integer(ik) :: i
       !
       i  = 1
       isparse = 0
       !
       do while(i<=Nterms)
         !
         if (ifromsparse(i)==ifull) then
           isparse = i
           exit 
         endif
         i = i + 1
       enddo
       !
       !if (isparse==0) then 
       !  write(6,"('FLbset1DNew: Could not localte isparse-term in the sparse-field for index  ',i6)") ifull
       !  stop 'FLbset1DNew: Could not localte a term in the sparse-field'
       !endif
       !
     end subroutine find_isparse_from_ifull
     !
  end subroutine FLbset1DNew



!
! Add new rotational basis set type
! 
  subroutine FLRotbset_new(Jmax)
  
    integer(ik),intent(in)      :: Jmax
    integer(ik)                 :: alloc,MatrixSize,k,J,Jk,iref,jx
    type(Basis1DT), pointer     :: bs  
    character(len=1)            :: char_
    real(ark)                    :: Bx,By,Bz,factor,Ix,Iy,Iz,Inertm(3)
    !
    if (job%verbose>=4) write(out,"(/'FLRotbset_new/start')") 
    !
    ! substitute for easier reference 
    !
    bs => bset%rot
    !
    !  Initialize the basis set 
    !
    bs%Size = (jmax+1)*(jmax+2)/2-1
    bs%Order = 0
    !
    MatrixSize = 7*(jmax+1)*(jmax+2)/2*5*2
    !
    ! Rotational constants computed at the equilibrium
    !
    Ix = sum(trove%mass(:)*( trove%a0(:,2)**2+ trove%a0(:,3)**2) )
    Iy = sum(trove%mass(:)*( trove%a0(:,1)**2+ trove%a0(:,3)**2) )
    Iz = sum(trove%mass(:)*( trove%a0(:,1)**2+ trove%a0(:,2)**2) )
    !
    factor = real(planck*avogno,ark)*real(1.0d+16,kind=rk)/(8.0_ark*pi*pi*real(vellgt,ark))
    !
    Bx = 0 ; By = 0 ; Bz = 0 
    !
    if (trove%lincoord/=1) Bx = factor/Ix
    if (trove%lincoord/=2) By = factor/Iy
    if (trove%lincoord/=3) Bz = factor/Iz
    !
    Inertm = (/Ix,Iy,Iz/)
    do jx = 1,3
      if (Inertm(jx)<sqrt(small_).and.trove%lincoord/=jx) then
        write (out,"('Too small Ix,Iy,Iz: ',3f14.6,'; if it is linear molecule change lincoord ',i0)") Ix,Iy,Iz,trove%lincoord
        stop 'Vanishing Ix,Iy,Iz'
      endif
    enddo
    !
    if (job%verbose>=1) then
      write (out,"(/'Rotational constants (Bx,By,Bz, cm-1): ',3f14.6)") Bx,By,Bz
    end if
    !
    !Bx = trove%g_rot(1,1)%field(1,trove%iPotmin)
    !By = trove%g_rot(2,2)%field(1,trove%iPotmin)
    !Bz = trove%g_rot(3,3)%field(1,trove%iPotmin)
    !
    if (job%verbose>=2) then
      write (out,"(' Basis type ',a,' needs ',f14.6,' Mbytes of memory (plus a bit)')") &
             trim(bs%name), real(rk*MatrixSize,kind=rk)/(1024.0_rk**2)
    end if
    !
    allocate (bs%matelements(7,(jmax+1)*(jmax+2)/2,-2:2,0:1),bs%ener0(0:bs%Size),stat=alloc)
    call ArrayStart('bs%matelements-rot',alloc,size(bs%matelements),kind(bs%matelements))
    call ArrayStart('bs%ener0-rot',alloc,size(bs%ener0),kind(bs%ener0))
    !
    ! Here we generate the rotational basis set matrix elements 
    !
    select case (trim(bs%type)) 
      case default
        write (out,"('FLRotbset_new: basis set ',a,' unknown')") trim(bs%type)
        stop 'FLRotbset_new: unknown basis set '
         !
      case('JKTAU')
          ! 
          bs%params    = 0
          !
          call ME_rotation(jmax,bs%matelements) 
          !
          if (job%verbose>=2) write (out,"(/'Zero-order Rotational energies')")
          if (job%verbose>=2) write (out,"( '       J       k          energy')")
          !
          do J =0,Jmax
             !
             do k =0,J
                !
                jk=k+(j*(j+1) )/2
                !
                bs%ener0(Jk) = 0.5_ark*real(J*(J+1_ik)-k**2,ark)*&
                              ( Bx+By ) + &
                              real(k**2,ark)*Bz
                !
                write (out,"(2i8,f16.4)") J,k,bs%ener0(Jk)
                !
             enddo
             !
          enddo
          !
    end select
    !
    !bs%Size = (jmax+1)*(jmax+2)/2
    bs%Order = jmax
    !
    !
    if (job%verbose>=4) write(out,"('FLRotbset_new/end')")   
    !
  end subroutine FLRotbset_new 
!
!
!
  subroutine Printbset1DInfo (ibs)
    integer(ik), intent(in)         :: ibs
    integer(ik)                     :: imode,i1,i2,iterm
    !
    type (Basis1DT), pointer :: bs
    !
    if (verbose<=0) return
    !
    !  Report grid parameters
    !
    bs => bset%bs1D(ibs)
    write (out,"(' 1D basis set #',i9,' (',a,')')") ibs,trim(bs%name)
    write (out,"(' Order    = ',i6)")    bs%Order
    write (out,"(' Size     = ',i6)")    bs%Size
    write (out,"(' Nmodes   = ',i6)")    bs%imodes

    if (bs%imodes<=30) then 
       write(out,"(' Modes: ',30i5)") (bs%mode(imode),imode=1,bs%imodes)
    endif
    
    write (out,"(' params  = ',3f18.8)") bs%params

    if (verbose>=6) then 
       !
       write (out,"(' mat.elements:')")
       !
       !
       do iterm = 0,trove%MaxOrder
           !
           do i1 = 0,bs%Size
              do i2 = 0,bs%Size
                 if (bset%dscr(ibs)%coord_kinet/=bset%dscr(ibs)%coord_kinet) &
                 write (out,"('mat-1(',3i4,')= ',f18.8)") iterm,i1,i2,bs%matelements(-1,iterm,i1,i2)
                !
                 write (out,"('mat0 (',3i4,')= ',f18.8)") iterm,i1,i2,bs%matelements(0,iterm,i1,i2)
                 write (out,"('mat1 (',3i4,')= ',f18.8)") iterm,i1,i2,bs%matelements(1,iterm,i1,i2)
                 write (out,"('mat2 (',3i4,')= ',f18.8)") iterm,i1,i2,bs%matelements(2,iterm,i1,i2)
                 write (out,"('mat3 (',3i4,')= ',f18.8)") iterm,i1,i2,bs%matelements(3,iterm,i1,i2)
              enddo
           enddo
       enddo
       !
    endif 
    !
  end subroutine Printbset1DInfo
!
!
! Primitive 1D harmonic energy 
!
  function FLHarmonicEnergy(imode,nu) result (f)
    integer(ik), intent(in)        :: imode
    integer(ik), intent(in)        :: nu

    integer(ik)                    :: nmodes,itype

    real(ark)                       :: f
    real(ark)                       :: omega
    

    nmodes = 0
    itype = 0 
    do while(itype<bset%n_bset1D_max.and.nmodes<imode) 
       itype = itype + 1 
       nmodes = nmodes+bset%bs1D(itype)%imodes
    enddo

    omega = bset%bs1d(itype)%params(2)

    f = (real(nu,kind=rk) + 0.5_ark)*omega

    if (job%verbose>=5) then 
        write(out,"('FLHarmonicEnergy: imode,nu = ',2i5,f18.8)") imode,nu,f
    endif


  end function FLHarmonicEnergy

   !
   ! Matrix elements calculations of the fields on the product of the primitive (1D) eigenfunctions
   ! <v1v2v3...| fields |w1w2w3...>, where fields = V, T_vib, T_rot, ot T_corr
   !
   function FLmatrixelements(norder,nu_i,nu_j,j) result(mat_elem)

      integer(ik),intent(in)        :: norder,nu_i(0:trove%Nmodes),nu_j(0:trove%Nmodes)
      integer(ik),intent(in),optional  :: j
      !character(len=cl),intent(in) :: fieldtype       ! Identifying name of the field
      real(rk)                      :: mat_elem

      integer(ik)                   :: imode,i,ibstype,iterm,k1,k2,pshift,k(trove%Nmodes)
      integer(ik)                   :: vl,vr,tau_i,tau_j,k_i,k_j,jk,dk
      real(rk)                      :: pot_t,gvib_t,mat(trove%Nmodes),grot_t,gcor_t
      type(FLpolynomT),pointer      :: fl
      type(Basis1DT),pointer        :: bs,bs_rot
      !real(rk)                      :: mat_legatee(trove%Nmodes,trove%Ncoeff)

         if (verbose>=6) write(out,"(/'FLmatrixelements/start: matrix elemnts for the hamiltonian ')") 
         !
         !call TimerStart('FLmatrixelements')
         
         ! Zero order for the potential function starts from the harmonic approximation (=2) 

         !if (norder<=1) then 
         !  if (norder+pshift-1>=0)             itmin = trove%RangeOrder(min(norder+pshift-1,trove%NPotOrder))+1
         !  if (norder+pshift<=trove%NPotOrder) itmax = trove%RangeOrder(min(norder+pshift  ,trove%NPotOrder))
         !endif 
         !
         pshift = job%pot_pt_shift
         !
         ! Extract large amplitude quantum numbers
         !
         vl = nu_i(trove%Nmodes) ; vr = nu_j(trove%Nmodes)
         !
         ! The J-free part has to be diagonal in terms of the rotational quanta
         !
         pot_t = 0
         gvib_t = 0 
         !
         if (nu_i(0)==nu_j(0)) then 
            ! 
            ! Potential part of the hamiltonian 
            !
            fl => trove%poten
            !
            !
            !mat_legatee = 0 
            !
            write(out,"('check FLIndexQ_legatee')")
            stop 'check FLIndexQ_legatee'
            !
            !do iterm = 1,fl%Ncoeff
            !   !
            !   k(:) = FLIndexQ(:,iterm)
            !   !
            !   ! For the zero order case we allow only diagonal terms
            !   !
            !   ! Check if the current iterm belongs to the present calculation case
            !   ! another words, if poten*xi^k belongs to the current perturb. order
            !   !
            !   !if (fl%iorder(iterm)==norder) then 
            !      !
            !      do ibstype = 1,bset%n_bset1D_max
            !         !
            !         imode = bset%bs1D(ibstype)%mode(1)
            !         !
            !         bs => bset%bs1D(ibstype)
            !         !
            !         do i = 1,bs%imodes
            !            !
            !            imode = bs%mode(i)
            !            !
            !            if (FLIndexQ_legatee(imode,iterm)/=0) then
            !              !
            !              if (FLIndexQ_legatee(imode,iterm)==-1) then 
            !                !
            !                mat_legatee(imode,iterm) = 1.0_rk
            !                !
            !              else
            !                !
            !                mat_legatee(imode,iterm) = mat_legatee(imode-1,FLIndexQ_legatee(imode,iterm))
            !                !
            !              endif
            !              !
            !              if (imode==trove%Nmodes) then
            !                !
            !                !mat(imode) = fl%me(iterm,vl,vr)
            !                mat_legatee(imode,iterm) = mat_legatee(imode,iterm) * fl%me(iterm,vl,vr)
            !                !
            !              else
            !                !
            !                !mat(imode) = bs%matelements(0,k(imode),nu_i(imode),nu_j(imode))
            !                mat_legatee(imode,iterm) = mat_legatee(imode,iterm) * bs%matelements(0,k(imode),nu_i(imode),nu_j(imode))
            !                !
            !              endif
            !              !
            !            endif   ! legatee/=0
            !            !
            !         enddo 
            !         ! 
            !      enddo
            !      !
!           !      pot_t = pot_t + product(mat(:))
            !      !
            !   !endif 
            !   !
            !enddo
            ! 
            !pot_t = 0
            !
            !do iterm = 1,fl%Ncoeff
            !   !
            !   if (fl%iorder(iterm)==norder) then 
            !      !
            !      pot_t = pot_t  + mat_legatee(trove%Nmodes,iterm)
            !      !
            !   endif 
            !   !
            !enddo


            ! 
            ! Pseudo-potential part of the hamiltonian 
            !
            !fl => trove%pseudo
            !
            !do iterm = 1,fl%Ncoeff
            !   !
            !   k(:) = FLIndexQ(:,iterm)
            !   !
            !   if (fl%iorder(iterm)==norder) then 
            !      !
            !      do ibstype = 1,bset%n_bset1D_max
            !         !
            !         imode = bset%bs1D(ibstype)%mode(1)
            !         !
            !         bs => bset%bs1D(ibstype)
            !         !
            !         do i = 1,bs%imodes
            !            !
            !            imode = bs%mode(i)
            !            !
            !            if (imode==trove%Nmodes) then
            !               !
            !               mat(imode) = fl%me(iterm,vl,vr)
            !               !
            !            else
            !               !
            !               mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
            !               !
            !            endif
            !            !
            !         enddo 
            !         ! 
            !      enddo
            !      pot_t = pot_t + product(mat(:))
            !      !
            !   endif 
            !   !
            !enddo 
            ! 
            ! Vibrational part of the kinetic operator g_vib
            !
            do k1 = 1,trove%Nmodes
               do k2 = 1,trove%Nmodes
                  !
                  fl => trove%g_vib(k1,k2)
                  !
                  do iterm = 1,fl%Ncoeff
                     !
                     k(:) = FLIndexQ(:,iterm)
                     !
                     ! Check if the current iterm belongs to the current perturb. order
                     !
                     if (fl%iorder(iterm)==norder) then 
                        !
                        do ibstype = 1,bset%n_bset1D_max
                           !
                           imode = bset%bs1D(ibstype)%mode(1)
                           !
                           bs => bset%bs1D(ibstype)
                           !
                           do i = 1,bs%imodes
                              !
                              imode = bs%mode(i)
                              !
                              if (imode==trove%Nmodes) then
                                 !
                                 mat(imode) = fl%me(iterm,vl,vr)
                                 !
                              else
                                 !
                                 if    (k1/=imode.and.k2/=imode) then 
                                   !
                                   mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                                   !
                                 elseif (k1==imode.and.k2/=imode) then
                                   !
                                   mat(imode) =-bs%matelements(1,k(imode),nu_j(imode),nu_i(imode))
                                   !
                                 elseif (k1/=imode.and.k2==imode) then
                                   !
                                   mat(imode) = bs%matelements(1,k(imode),nu_i(imode),nu_j(imode))
                                   !
                                 else !   if (k1==imode.and.k2==imode) then
                                   !
                                   mat(imode) = bs%matelements(2,k(imode),nu_i(imode),nu_j(imode))
                                   !
                                 endif
                                 !
                                 !
                              endif
                              !
                           enddo 
                           ! 
                        enddo
                        !
                        if (verbose>=7) write(out,"('k1,k2,iterm,iorder ',4i8)") k1,k2,iterm,fl%iorder(iterm)
                        if (verbose>=7) write(out,"('k1,k2,iterm,iorder ',31f18.8)") fl%field(iterm,1),mat(:)
                        !
                        gvib_t = gvib_t + product(mat(:))
                        !
                     endif
                     !
                  enddo
                  !
                  continue
                  ! 
               enddo
            enddo
            !
         endif 
         !
         ! Rotational and Coriolis parts of the kinetic energy operator are turned on when FLrotation is .true.
         !
         grot_t = 0 
         gcor_t = 0 
         !
         if (FLrotation.and.present(J)) then
          if (J/=0) then
            !
            ! Extract rotaitonal quantum numbers
            !
            if (nu_i(0)==0) then 
               !
               tau_i = mod(J,2)
               k_i   = 0
               !
            else
               !
               tau_i = mod(nu_i(0),2)
               k_i   = ( nu_i(0)+tau_i )/2
               !
            endif
            !
            if (nu_j(0)==0) then 
               !
               tau_j = mod(J,2)
               k_j   = 0
               !
            else
               !
               tau_j = mod(nu_j(0),2)
               k_j   = ( nu_j(0)+tau_j )/2
               !
            endif
            !
            Jk = 1+k_i+(j*(j+1) )/2
            dk = k_i - k_j
            !
            bs_rot => bset%rot
            ! 
            ! Rotational part of the kinetic operator g_rot
            !
            if (abs(dk)<=2) then 
               !
               do iterm = 1,trove%g_rot(1,1)%Ncoeff
                  !
                  k(:) = FLIndexQ(:,iterm)
                  !
                  ! Check if the current iterm belongs to the current perturb. order
                  !
                  if (trove%g_rot(1,1)%iorder(iterm)==norder) then 
                     !
                     do ibstype = 1,bset%n_bset1D_max
                        !
                        imode = bset%bs1D(ibstype)%mode(1)
                        !
                        bs => bset%bs1D(ibstype)
                        !
                        do i = 1,bs%imodes
                           !
                           imode = bs%mode(i)
                           !
                           if (imode/=trove%Nmodes) then
                              !
                              mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                              !
                           else
                              !
                              ! four terms that come together: 
                              ! 1/2(Jx^2 - Jy^2) and (JxJy+JyJx)
                              !     (JxJz+JzJx)  and (JzJy+JyJz)
                              !
                              ! first  - all even contributions 
                              !
                              if (tau_i==tau_j) then 
                                 !
                                 ! Jx2y2 
                                 !
                                 mat(imode) =  0.5_ark*&
                                             ( trove%g_rot(1,1)%me(iterm,vl,vr)-trove%g_rot(2,2)%me(iterm,vl,vr) )*&
                                               bs_rot%matelements(4,Jk,dk,tau_i)
                                 !
                                 ! Jxz
                                 !
                                 mat(imode) = mat(imode) +& 
                                 trove%g_rot(1,3)%me(iterm,vl,vr)*bs_rot%matelements(6,Jk,dk,tau_i)
                                 !
                                 ! and now odd ones 
                              else
                                 !
                                 ! Jxy
                                 !
                                 mat(imode) = trove%g_rot(1,2)%me(iterm,vl,vr)*bs_rot%matelements(5,Jk,dk,tau_i)
                                 !
                                 ! Jyz
                                 !
                                 mat(imode) = mat(imode) +& 
                                              trove%g_rot(2,3)%me(iterm,vl,vr)*bs_rot%matelements(7,Jk,dk,tau_i)
                                 !
                              endif 
                              !
                              if (tau_i==tau_j.and.k_i==k_j) then 
                                 !
                                 ! two terms A*Jx^2+BJy^2+C*Jz^2
                                 !
                                 mat(imode) = mat(imode) + 0.5_ark*real(J*(J+1_ik)-k_i**2,rk)*&
                                            ( trove%g_rot(1,1)%me(iterm,vl,vr)+trove%g_rot(2,2)%me(iterm,vl,vr) )
                                 mat(imode) = mat(imode) + real(k_i**2,rk)*trove%g_rot(3,3)%me(iterm,vl,vr) 
                                 !
                              endif 
                              !
                           endif
                           !
                        enddo 
                        ! 
                     enddo
                     !
                     grot_t = grot_t + product(mat(:))
                     !
                  endif
                  !
               enddo
               !
               continue
               !
            endif 
            ! 
            ! Coriolis part of the kinetic operator g_cor
            !
            if (abs(dk)<=1) then 
               !
               do k1 = 1,trove%Nmodes
                  !
                  do iterm = 1,trove%g_cor(k1,1)%Ncoeff
                     !
                     k(:) = FLIndexQ(:,iterm)
                     !
                     ! Check if the current iterm belongs to the current perturb. order
                     !
                     if (trove%g_cor(k1,1)%iorder(iterm)==norder) then 
                        !
                        do ibstype = 1,bset%n_bset1D_max
                           !
                           imode = bset%bs1D(ibstype)%mode(1)
                           !
                           bs => bset%bs1D(ibstype)
                           !
                           do i = 1,bs%imodes
                              !
                              imode = bs%mode(i)
                              !
                              if (imode/=trove%Nmodes) then
                                 !
                                 if (k1/=imode) then 
                                    !
                                    mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                                    !
                                 else
                                    !
                                    mat(imode) =bs%matelements(1,k(imode),nu_i(imode),nu_j(imode))-&
                                                bs%matelements(1,k(imode),nu_j(imode),nu_i(imode))
                                    !
                                 endif
                                 !
                              else
                                 !
                                 ! first  - all even contributions 
                                 !
                                 if (tau_i==tau_j) then 
                                    !
                                    mat(imode) =  trove%g_cor(k1,2)%me(iterm,vl,vr)*bs_rot%matelements(2,Jk,dk,tau_i)
                                    !
                                    ! and now odd ones 
                                 else
                                    !
                                    mat(imode) =  trove%g_cor(k1,1)%me(iterm,vl,vr)*bs_rot%matelements(1,Jk,dk,tau_i)+&
                                                  trove%g_cor(k1,3)%me(iterm,vl,vr)*bs_rot%matelements(3,Jk,dk,tau_i)
                                    !
                                 endif 
                                 !
                              endif
                              !
                           enddo 
                           ! 
                        enddo
                        !
                        gcor_t = gcor_t + product(mat(:))
                        !
                     endif
                     !
                  enddo
                  !
                  continue
                  !
               enddo 
               !
            endif 
            ! 
            !
          endif 
         endif 
         !
         mat_elem = pot_t-0.5_ark*(gvib_t-grot_t-gcor_t)
         !
         if (verbose>=7) write(out,"('FLmatrixelements:  pot_t,gvib_t,grot_t,gcor_t ',4f18.8)") pot_t,gvib_t,grot_t,gcor_t
         !
         !call TimerStart('FLmatrixelements')
         !
      if (verbose>=6) write(out,"('FLmatrixelements/end')") 
    
   end function FLmatrixelements



   !
   ! Matrix elements calculations of a term xi1^i1*xi2^i2*xi3^i3... 
   ! on the product of the primitive 1D eigenfunctions
   ! <v1v2v3...|  xi1^i1*xi2^i2*xi3^i3... |w1w2w3...>
   !
   function FLmatrixelements_single_term(nu_i,nu_j,term) result(mat_elem)

      integer(ik),intent(in)        :: nu_i(0:trove%Nmodes),nu_j(0:trove%Nmodes)
      integer(ik),intent(in)        :: term(0:trove%Nmodes)
      real(rk)                      :: mat_elem
      !
      integer(ik)                   :: imode,i,ibstype,iterm
      real(rk)                      :: mat(trove%Nmodes)
      type(FLpolynomT),pointer      :: fl
      type(Basis1DT),pointer        :: bs

         if (verbose>=6) write(out,"(/'FLmatrixelements_single_term/start: matrix elemnts for the hamiltonian ')") 
         !
         ! This is a J-free part and has to be diagonal in terms of the rotational quanta
         !
         if (nu_i(0)/=nu_j(0)) return 
         !
         if (nu_i(0)/=nu_j(0)) return 
         !
         fl => trove%poten
         ! 
         do ibstype = 1,bset%n_bset1D_max
            !
            imode = bset%bs1D(ibstype)%mode(1)
            !
            bs => bset%bs1D(ibstype)
            !
            do i = 1,bs%imodes
               !
               imode = bs%mode(i)
               !
               if (imode==trove%Nmodes) then
                  !
                  write(out,"('FLmatrixelements_single_term: ')")
                  write(out,"('This routine is not for the class that contains the last mode ',i8)") imode
                  stop 'The last mode cannot appear here'
                  !
               else
                  !
                  mat(imode) = bs%matelements(term(0),term(imode),nu_i(imode),nu_j(imode))
                  !
               endif
               !
            enddo 
            ! 
         enddo
         !
         mat_elem = product(mat(:))
         !
      if (verbose>=6) write(out,"('FLmatrixelements_single_term/end')") 
    
   end function FLmatrixelements_single_term



   subroutine FLread_fields_dimensions(poten_N,gvib_N,grot_N,gcor_N,potorder,kinorder,extForder,jmax,extF_N,L2vib_N)
     !
     integer(ik),intent(out)  :: poten_N,gvib_N,grot_N,gcor_N,L2vib_N,potorder,kinorder,extForder,jmax,extF_N(:)
     integer(ik)              :: i
       !
       poten_N = trove%poten%Ncoeff
       gvib_N = trove%g_vib (1,1)%Ncoeff
       grot_N = trove%g_rot (1,1)%Ncoeff
       gcor_N = trove%g_cor(1,1)%Ncoeff
       potorder = trove%Npotorder
       kinorder = trove%Nkinorder
       extForder = trove%NExtOrder
       L2vib_N = 0
       !
       if (FLL2_coeffs) L2vib_N = trove%L2_vib(1,1)%Ncoeff
       !
       if (FLextF_coeffs) then 
         !
         extF_N(1:extF%rank) = trove%extF(1:extF%rank)%Ncoeff
         !
       else
         !
         extF_N = 0 
         !
       endif 
       !
       jmax   = trove%jmax
       !
   end subroutine FLread_fields_dimensions

   !
   ! trasnfer the dimentsions of the fields to perturbation assuming sparse representation 
   !
   function FLread_fields_dimension_field(job_str,k1,k2) result (Ncoeff)
     !
     character(len=cl),intent(in) :: job_str
     integer(ik),intent(in)       :: k1,k2 
     integer(ik)                  :: Ncoeff
       !
      select case(trim(job_str)) 
          !
        case('poten')
          !
          Ncoeff = trove%poten%Ncoeff
          !
        case('gvib')
          !
          Ncoeff = trove%g_vib (k1,k2)%Ncoeff
          !
        case('grot')
          !
          Ncoeff = trove%g_rot(k1,k2)%Ncoeff
          !
        case('gcor')
          !
          Ncoeff = trove%g_cor(k1,k2)%Ncoeff
          !
        case('externalF')
          !
          Ncoeff = trove%extF(k1)%Ncoeff
          !
        case('L2_vib')
          !
          Ncoeff = trove%L2_vib(k1,k2)%Ncoeff
          !
        case default
          !
          write (out,"('FLread_fields_dimension_field: job_str ',a,' unknown')") job_str
          stop 'FLread_fields_dimension_field - bad job_str'
          !
      end select 
      !
   end function FLread_fields_dimension_field



   function FLread_extF_rank() result (rank)
     !
     integer(ik) :: rank
     !
     rank = extF%rank
     !
   end function FLread_extF_rank
   
   !
   ! This is to compute the values of all kinetic energy operator terms 
   ! as well as of the potential energy funciton at a given geometry. 
   !
   recursive subroutine FLcalc_poten_kinet_dvr(dchi,irho,poten,gvib,grot,gcor,extfield,reduced_model)
     !
     real(ark),intent(in) :: dchi(trove%Nmodes)
     integer(ik),intent(in) :: irho
     real(ark),intent(out):: poten
     real(ark),intent(out):: gvib(trove%Nmodes,trove%Nmodes)
     real(ark),intent(out):: grot(3,3)
     real(ark),intent(out):: gcor(trove%Nmodes,3)
     real(ark),intent(out):: extfield(:)
     logical,intent(in)  ::  reduced_model
     real(ark)           :: extfield_(extF%rank)
     !
     real(ark)            :: p_t,f_t,g_t,pseudo
     integer(ik)         :: k(trove%Nmodes),iterm,igeom,k1,k2,i
     real(ark)           :: xi(trove%Nmodes)
     type(FLpolynomT),pointer  :: fl
       !
       ! Potential energy contribution
       !
       if (.not.reduced_model) then
         ! 
         poten = FLpoten_linearized_dchi(dchi,irho)
         !
         call FLDVR_gmat_dvr(dchi,irho,gvib,grot,gcor,pseudo)
         !
         poten = poten + pseudo
         !
         if (FLextF_coeffs) then 
           !
           call dms4chi(irho,dchi,extfield_)
           !
           if (verbose>=6) write(out,"('rank = ',i7)") extF%rank
           !
           if (verbose>=6) write(out,"('extfield_ = ',50f18.8)") extfield_
           !
           extfield = extfield_
           !
         endif 
         !
       else
         !
         do i = 1,trove%Nmodes_e
            xi(i) = MLcoord_direct(dchi(i),2,i) 
         enddo
         !
         ! Contribution from the pseudo-potential part 
         !
         fl => trove%poten
         !
         p_t = 0
         !
         do iterm = 1,fl%Ncoeff
            !
            if (fl%iorder(iterm)/=0) cycle
            !
            k(:) = FLIndexQ(:,iterm)
            !
            f_t = product(xi(:)**k(:))
            !
            p_t = p_t + fl%field(iterm,irho)*f_t
            !
         enddo
         !
         poten = p_t
         !
         ! Switch to the kinetic-xi coordinate:
         !
         do i = 1,trove%Nmodes_e
            xi(i) = MLcoord_direct(dchi(i),1,i) 
         enddo
         !
         ! Contribution from the pseudo-potential part 
         !
         fl => trove%pseudo
         !
         g_t = 0
         !
         do iterm = 1,fl%Ncoeff
            !
            if (fl%iorder(iterm)/=0) cycle
            !
            k(:) = FLIndexQ(:,iterm)
            !
            f_t = product(xi(:)**k(:))
            !
            g_t = g_t + fl%field(iterm,irho)*f_t
            !
         enddo
         !
         ! Combine with the potential energy part 
         !
         poten = poten + g_t
         !
         ! Vibrational part of the kinetic operator g_vib
         !
         do k1 = 1,trove%Nmodes
            do k2 = 1,trove%Nmodes
               !
               fl => trove%g_vib(k1,k2)
               !
               g_t = 0
               !
               do iterm = 1,fl%Ncoeff
                  !
                  if (fl%iorder(iterm)/=0) cycle
                  !
                  k(:) = FLIndexQ(:,iterm)
                  !
                  f_t = product(xi(:)**k(:))
                  !
                  g_t = g_t + fl%field(iterm,irho)*f_t
                  !
               enddo
               !
               gvib(k1,k2) = g_t
               ! 
            enddo
         enddo
         ! 
         ! Rotational part and Coriolis parts of the kinetic operator g_rot
         !
         !if (treat_rotation) then
         !
         do k1 = 1,3
            do k2 = 1,3
               !
               fl => trove%g_rot(k1,k2)
               !
               g_t = 0
               !
               do iterm = 1,fl%Ncoeff
                  !
                  if (fl%iorder(iterm)/=0) cycle
                  !
                  k(:) = FLIndexQ(:,iterm)
                  !
                  f_t = product(xi(:)**k(:))
                  !
                  g_t = g_t + fl%field(iterm,irho)*f_t
                  !
               enddo
               !
               grot(k1,k2) = g_t
               ! 
            enddo
         enddo
         !
         !
         do k1 = 1,trove%Nmodes
            do k2 = 1,3
               !
               fl => trove%g_cor(k1,k2)
               !
               g_t = 0
               !
               do iterm = 1,fl%Ncoeff
                  !
                  if (fl%iorder(iterm)/=0) cycle
                  !
                  k(:) = FLIndexQ(:,iterm)
                  !
                  f_t = product(xi(:)**k(:))
                  !
                  g_t = g_t + fl%field(iterm,irho)*f_t
                  !
               enddo
               !
               gcor(k1,k2) = g_t
               ! 
            enddo
         enddo
         !
       endif 
       ! 
   end subroutine FLcalc_poten_kinet_dvr



   subroutine FLread_coeff_matelem(job_str,k1,k2,field)
     !
     character(len=cl),intent(in) :: job_str
     integer(ik),intent(in)       :: k1,k2 
     real(rk),intent(out)         :: field(:,:,:) !
      !
      !if (size(field,dim=1)/=size(bset%rot%matelements,dim=2)) 
      !
      select case(trim(job_str)) 
        !
        case('rot')
          !
          field(:,:,:) = bset%rot%matelements(k1,:,:,:)
          !
        case('poten')
          !
          field(:,:,:) = trove%poten%me(:,:,:)
          !
        case('gvib')
          !
          field(:,:,:) = trove%g_vib(k1,k2)%me(:,:,:)
          !
        case('grot')
          !
          field(:,:,:) = trove%g_rot(k1,k2)%me(:,:,:)
          !
        case('gcor')
          !
          field(:,:,:) = trove%g_cor(k1,k2)%me(:,:,:)
          !
        case('externalF')
          !
          field(:,:,:) = trove%extF(k1)%me(:,:,:)
          !
        case('vib')
          !
          field(:,:,:) = bset%bs1D(k1)%matelements(k2,:,:,:)
          !
        case('L2_vib')
          !
          field(:,:,:) = trove%L2_vib(k1,k2)%me(:,:,:)
          !
        case default
          !
          write (out,"('FLread_coeff_matelem: job_str ',a,' unknown')") job_str
          stop 'FLread_coeff_matelem - bad job_str'
          !
      end select 
      !
   end subroutine FLread_coeff_matelem
   !
   ! This is to transfer IndexQ to perturbation 
   !
   subroutine FLread_IndexQ_field(job_str,k1,k2,IndexQ)
     !
     character(len=cl),intent(in) :: job_str
     integer(ik),intent(in)       :: k1,k2 
     integer(ik),intent(out)         :: IndexQ(:,:)
     type(FLpolynomT),pointer        :: fl
      !
      !if (size(field,dim=1)/=size(bset%rot%matelements,dim=2)) 
      !
      select case(trim(job_str)) 
          !
        case('poten')
          !
          fl => trove%poten
          !
        case('gvib')
          !
          fl =>  trove%g_vib(k1,k2)
          !
        case('grot')
          !
          fl => trove%g_rot(k1,k2)
          !
        case('gcor')
          !
          fl => trove%g_cor(k1,k2)
          !
        case('externalF')
          !
          fl => trove%extF(k1)
          !
        case('L2_vib')
          !
          fl => trove%L2_vib(k1,k2)
          !
        case default
          !
          write (out,"('FLread_IndexQ_field: job_str ',a,' unknown')") job_str
          stop 'FLread_IndexQ_field - bad job_str'
          !
      end select 
      !
      IndexQ(:,:) = fl%IndexQ(:,:)
      !
   end subroutine FLread_IndexQ_field   
   !
   !
   subroutine FLread_iorder_send(job_str,k1,k2,field)
     !
     character(len=cl),intent(in) :: job_str
     integer(ik),intent(in)       :: k1,k2 
     integer(ik),intent(in)       :: field(:)
      !
      !
      !if (size(field,dim=1)/=size(bset%rot%matelements,dim=2)) 
      !
      select case(trim(job_str)) 
          !
        case('poten')
          !
          trove%poten%iorder(:) = field(:)
          !
        case('gvib')
          !
          trove%g_vib(k1,k2)%iorder(:) = field(:)
          !
        case('grot')
          !
          trove%g_rot(k1,k2)%iorder(:) = field(:)
          !
        case('gcor')
          !
          trove%g_cor(k1,k2)%iorder(:) = field(:)
          !
        case('L2vib')
          !
          trove%L2_vib(k1,k2)%iorder(:) = field(:)
          !
        case('externalF')
          !
          trove%extF(k1)%iorder(:) = field(:)
          !
        case default
          !
          write (out,"('FLread_iorder_send: job_str ',a,' unknown')") job_str
          stop 'FLread_iorder_send - bad job_str'
          !
      end select 
      !
   end subroutine FLread_iorder_send



   subroutine FLread_rot_matelem(k,rot_me)
     !
     integer(ik),intent(in) :: k             ! number of the element 
     real(rk),intent(out)   :: rot_me(:,:,:) !
      !
      if (k<1.or.k>7) then 
        !
        write(out,"('FLread_rot_matelem: illigal number for rot.mat.elem ',i8)") k
        stop 'FLread_rot_matelem: illigal number'
        !
      endif 
      !
      rot_me(:,:,:) = bset%rot%matelements(k,:,:,:)
      !
   end subroutine FLread_rot_matelem

  !
  ! Deallocate objects related to the primitive matrix elements in the 
  ! standard power representation of Hamiltonian expansion 
  !
   subroutine  FLfree_primitive_objects 
        !
        integer(ik) :: itype,alloc

        !
        do itype = 0,bset%n_bset1D_max
          !
          if (itype==0.and.FLrotation) then 
             !
             deallocate(bset%rot%matelements)
             deallocate(bset%rot%ener0)
             call ArrayStop('bs%matelements-rot')
             call ArrayStop('bs%ener0-rot')
             !
           elseif(itype/=0) then
             !
             deallocate (bset%bs1D(itype)%matelements,bset%bs1D(itype)%ener0,stat=alloc)
             !
           endif 
           !
        enddo
        !
        if (associated(trove%poten)) nullify(trove%poten) !deallocate(trove%poten)
        if (associated(trove%g_vib)) deallocate(trove%g_vib)
        if (associated(trove%g_rot)) deallocate(trove%g_rot)
        if (associated(trove%g_cor)) deallocate(trove%g_cor)
        if (associated(trove%pseudo)) then 
            !deallocate(trove%pseudo)
            nullify(trove%pseudo)
            call ArrayStop('pseudo')
        endif 
        if (associated(trove%extF)) then
           deallocate(trove%extF)
           call ArrayStop('extF')
        endif
        !
        call ArrayStop('g_vib')
        call ArrayStop('g_rot')
        call ArrayStop('g_cor')
        call ArrayStop('poten')
        call ArrayStop('bs%matelements')
        call ArrayStop('bs%ener0')
        !
   end subroutine  FLfree_primitive_objects 


   !
   ! Matrix elements calculations of the terms containing 
   ! values of the potential and kinetic energy expansions
   ! <v1v2v3...|  f |w1w2w3...>
   !
   subroutine FLmatrixelements_single_terms(nu_i,nu_j,poten,gvib,grot,gcor)

      integer(ik),intent(in)        :: nu_i(0:trove%Nmodes),nu_j(0:trove%Nmodes)
      real(rk),intent(out)           :: poten(trove%poten%Ncoeff)
      real(rk),intent(out)           :: gvib (trove%Nmodes,trove%Nmodes,trove%g_vib(1,1)%Ncoeff)
      real(rk),intent(out),optional  :: grot (3,3,trove%g_rot(1,1)%Ncoeff)
      real(rk),intent(out),optional  :: gcor(trove%Nmodes,3,trove%g_cor(1,1)%Ncoeff)
      !
      integer(ik)                   :: vl,vr,imode,i,ibstype,iterm,k1,k2,k(trove%Nmodes)
      real(rk)                      :: mat(trove%Nmodes)
      type(FLpolynomT),pointer      :: fl
      type(Basis1DT),pointer        :: bs

         if (verbose>=6) write(out,"(/'FLmatrixelements_single_terms/start: matrix elemnts for the hamiltonian ')") 
         !
         ! This is a J-free part and has to be diagonal in terms of the rotational quanta
         !
         if (nu_i(0)>0.or.nu_j(0)>0) then
            !
            write(out,"('FLmatrixelements_single_terms: ')")
            write(out,"('This routine is not for J/=0 ',2i8)") nu_i(0),nu_j(0)
            stop 'This routine is only for J=0'
            !
         endif 
         !
         ! Extract large amplitude quantum numbers
         !
         vl = nu_i(trove%Nmodes) ; vr = nu_j(trove%Nmodes)
         !
         ! Potential part of the hamiltonian 
         !
         fl => trove%poten
         !
         do iterm = 1,fl%Ncoeff
            !
            k(:) = FLIndexQ(:,iterm)
            !
            mat = 1
            !
            mode_loop_pot : do ibstype = 1,bset%n_bset1D_max
               !
               imode = bset%bs1D(ibstype)%mode(1)
               !
               if (nu_i(imode)<0.or.nu_i(imode)<0) cycle mode_loop_pot
               !
               bs => bset%bs1D(ibstype)
               !
               do i = 1,bs%imodes
                  !
                  imode = bs%mode(i)
                  !
                  if (imode==trove%Nmodes) then
                     !
                     mat(imode) = fl%me(iterm,vl,vr)
                     !
                  else
                     !
                     mat(imode) = bs%matelements(0,k(imode),nu_i(imode),nu_j(imode))
                     !
                  endif
                  !
               enddo 
               ! 
            enddo mode_loop_pot
            !
            poten(iterm) = product(mat(:))
            !
         enddo 
         ! 
         ! Vibrational part of the kinetic operator g_vib
         !
         do k1 = 1,trove%Nmodes
            do k2 = 1,trove%Nmodes
               !
               fl => trove%g_vib(k1,k2)
               !
               do iterm = 1,fl%Ncoeff
                  !
                  k(:) = FLIndexQ(:,iterm)
                  !
                  ! Check if the current iterm belongs to the current perturb. order
                  !
                  mode_loop_gvib : do ibstype = 1,bset%n_bset1D_max
                     !
                     imode = bset%bs1D(ibstype)%mode(1)
                     !
                     if (nu_i(imode)<0.or.nu_i(imode)<0) cycle mode_loop_gvib
                     !
                     bs => bset%bs1D(ibstype)
                     !
                     do i = 1,bs%imodes
                        !
                        imode = bs%mode(i)
                        !
                        if (imode==trove%Nmodes) then
                           !
                           mat(imode) = fl%me(iterm,vl,vr)
                           !
                        else
                           !
                           if    (k1/=imode.and.k2/=imode) then 
                             !
                             mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                             !
                           elseif (k1==imode.and.k2/=imode) then
                             !
                             mat(imode) =-bs%matelements(1,k(imode),nu_j(imode),nu_i(imode))
                             !
                           elseif (k1/=imode.and.k2==imode) then
                             !
                             mat(imode) = bs%matelements(1,k(imode),nu_i(imode),nu_j(imode))
                             !
                           else !   if (k1==imode.and.k2==imode) then
                             !
                             mat(imode) = bs%matelements(2,k(imode),nu_i(imode),nu_j(imode))
                             !
                           endif
                           !
                           !
                        endif
                        !
                     enddo 
                     ! 
                  enddo mode_loop_gvib
                  !
                  gvib(k1,k2,iterm) = product(mat(:))
                  !
               enddo
               ! 
            enddo
         enddo
         !
         ! Rotational and Coriolis parts of the kinetic energy operator are turned on when FLrotation is .true.
         !
         if (FLrotation.and.nu_i(0)>0.and.nu_j(0)>0) then
            !
            if (.not.present(grot).or..not.present(gcor)) then 
              !
              write(out,"('FLmatrixelements_single_terms: grot or gcor are not present for FLrotation')")
              stop 'FLmatrixelements_single_terms: grot or gcor is not present'
              !
            endif
            ! 
            ! Rotational part of the kinetic operator g_rot
            !
            do k1 = 1,3
               do k2 = 1,3
                  !
                  fl => trove%g_rot(k1,k2)
                  do iterm = 1,trove%g_rot(1,1)%Ncoeff
                    !
                    k(:) = FLIndexQ(:,iterm)
                    !
                    mode_loop_grot : do ibstype = 1,bset%n_bset1D_max
                       !
                       imode = bset%bs1D(ibstype)%mode(1)
                       !
                       if (nu_i(imode)<0.or.nu_i(imode)<0) cycle mode_loop_grot
                       !
                       bs => bset%bs1D(ibstype)
                       !
                       do i = 1,bs%imodes
                          !
                          imode = bs%mode(i)
                          !
                          if (imode/=trove%Nmodes) then
                             !
                             mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                             !
                          else
                             !
                             mat(imode) = fl%me(iterm,vl,vr)
                             !
                          endif
                          !
                       enddo 
                       ! 
                    enddo mode_loop_grot
                    !
                    grot(k1,k2,iterm) = product(mat(:))
                    !
                  enddo 
               enddo 
            enddo
            ! 
            ! Coriolis part of the kinetic operator g_cor
            !
            do k1 = 1,trove%Nmodes
               do k2 = 1,3
                  !
                  fl => trove%g_cor(k1,k2)
                  do iterm = 1,trove%g_rot(1,1)%Ncoeff
                    !
                    k(:) = FLIndexQ(:,iterm)
                    !
                    mode_loop_gcor : do ibstype = 1,bset%n_bset1D_max
                       !
                       imode = bset%bs1D(ibstype)%mode(1)
                       !
                       if (nu_i(imode)<0.or.nu_i(imode)<0) cycle mode_loop_gcor
                       !
                       bs => bset%bs1D(ibstype)
                       !
                       do i = 1,bs%imodes
                          !
                          imode = bs%mode(i)
                          !
                          if (imode/=trove%Nmodes) then
                             !
                             mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                             !
                          else
                             !
                             mat(imode) = fl%me(iterm,vl,vr)
                             !
                          endif
                          !
                       enddo 
                       ! 
                    enddo mode_loop_gcor
                    !
                    gcor(k1,k2,iterm) = product(mat(:))
                    !
                  enddo 
               enddo 
            enddo
            !
         endif 
         !
      if (verbose>=6) write(out,"('FLmatrixelements_single_terms/end')") 
    
   end subroutine FLmatrixelements_single_terms








   !
   ! Matrix elements calculations of the terms containing 
   ! values of the potential and kinetic energy expansions
   ! <v1v2v3...|  f |w1w2w3...>
   !
   subroutine FLmatrixelements_single_iterm(im1,im2,k,nu_i,nu_j,poten,gvib,grot,gcor)

      integer(ik),intent(in)        :: im1,im2
      integer(ik),intent(in)        :: k(trove%Nmodes)   ! term to be evaluated
      integer(ik),intent(in)        :: nu_i(0:trove%Nmodes),nu_j(0:trove%Nmodes)
      real(rk),intent(out),optional  :: poten
      real(rk),intent(out),optional  :: gvib(trove%Nmodes,trove%Nmodes)
      real(rk),intent(out),optional  :: grot(3,3)
      real(rk),intent(out),optional  :: gcor(trove%Nmodes,3)
      !
      integer(ik)                   :: vl,vr,imode,i,ibstype,k1,k2,iterm
      real(rk)                      :: mat(trove%Nmodes)
      type(FLpolynomT),pointer      :: fl
      type(Basis1DT),pointer        :: bs

         if (verbose>=6) write(out,"(/'FLmatrixelements_single_iterm/start: matrix elemnts for the hamiltonian ')") 
         !
         ! Extract large amplitude quantum numbers
         !
         vl = nu_i(trove%Nmodes) ; vr = nu_j(trove%Nmodes)
         !
         ! and the corresponding power 
         !
         iterm = k(trove%Nmodes)
         !
         poten  = 0
         gvib  = 0 
         !
         ! Potential part of the hamiltonian 
         !
         if (iterm<=trove%poten%Ncoeff) then 
            !
            fl => trove%poten
            !
            do imode = im1,im2
               !
               bs => bset%bs1D( bset%dscr(imode)%species )
               !
               if (imode==trove%Nmodes) then
                  !
                  mat(imode) = fl%me(iterm,vl,vr)
                  !
               else
                  !
                  mat(imode) = bs%matelements(0,k(imode),nu_i(imode),nu_j(imode))
                  !
               endif
               !
            enddo 
            !
            poten = product(mat(im1:im2))
            !
         endif 

         ! 
         ! Vibrational part of the kinetic operator g_vib
         !
         if (iterm<=trove%g_vib(1,1)%Ncoeff) then 
            !
            do k1 = 1,trove%Nmodes
               do k2 = 1,trove%Nmodes
                  !
                  fl => trove%g_vib(k1,k2)
                  !
                  do imode = im1,im2
                     !
                     bs => bset%bs1D( bset%dscr(imode)%species )
                     !
                     if (imode==trove%Nmodes) then
                        !
                        mat(imode) = fl%me(iterm,vl,vr)
                        !
                     else
                        !
                        if    (k1/=imode.and.k2/=imode) then 
                          !
                          mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                          !
                        elseif (k1==imode.and.k2/=imode) then
                          !
                          mat(imode) =-bs%matelements(1,k(imode),nu_j(imode),nu_i(imode))
                          !
                        elseif (k1/=imode.and.k2==imode) then
                          !
                          mat(imode) = bs%matelements(1,k(imode),nu_i(imode),nu_j(imode))
                          !
                        else !   if (k1==imode.and.k2==imode) then
                          !
                          mat(imode) = bs%matelements(2,k(imode),nu_i(imode),nu_j(imode))
                          !
                        endif
                        !
                        !
                     endif
                     !
                  enddo 
                  !
                  gvib(k1,k2) = product(mat(im1:im2))
                  ! 
               enddo
            enddo
            !
         endif 
         !
         ! Rotational and Coriolis parts of the kinetic energy operator are turned on when FLrotation is .true.
         !
         if (FLrotation.and.nu_i(0)>0.and.nu_j(0)>0) then
            !
            if (.not.present(grot).or..not.present(gcor)) then 
              !
              write(out,"('FLmatrixelements_single_iterm: grot or gcor are not present for FLrotation')")
              stop 'FLmatrixelements_single_iterm: grot or gcor is not present'
              !
            endif
            !
            grot  = 0 
            gcor = 0
            ! 
            ! Rotational part of the kinetic operator g_rot
            !
            if (iterm<=trove%g_rot(1,1)%Ncoeff) then 
               !
               do k1 = 1,3
                  do k2 = 1,3
                    !
                    fl => trove%g_rot(k1,k2)
                    !
                    mat = 1.0_ark
                    !
                    do imode = im1,im2
                       !
                       bs => bset%bs1D( bset%dscr(imode)%species )
                       !
                       if (imode/=trove%Nmodes) then
                          !
                          mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                          !
                       else
                          !
                          mat(imode) = fl%me(iterm,vl,vr)
                          !
                       endif
                       !
                    enddo 
                    !
                    grot(k1,k2) = product(mat(im1:im2))
                    !
                  enddo 
               enddo
               !
            endif 
            ! 
            ! Coriolis part of the kinetic operator g_cor
            !
            if (iterm<=trove%g_cor(1,1)%Ncoeff) then 
               !
               do k1 = 1,trove%Nmodes
                  do k2 = 1,3
                    !
                    fl => trove%g_cor(k1,k2)
                    !
                    do imode = im1,im2
                       !
                       bs => bset%bs1D( bset%dscr(imode)%species )
                       !
                       if (imode/=trove%Nmodes) then
                          !
                          if (k1/=imode) then 
                             !
                             mat(imode) = bs%matelements(-1,k(imode),nu_i(imode),nu_j(imode))
                             !
                          else
                             !
                             mat(imode) =bs%matelements(1,k(imode),nu_i(imode),nu_j(imode))-&
                                         bs%matelements(1,k(imode),nu_j(imode),nu_i(imode))
                             !
                          endif
                          !
                       else
                          !
                          mat(imode) = fl%me(iterm,vl,vr)
                          !
                       endif
                       !
                    enddo 
                    !
                    gcor(k1,k2) = product(mat(im1:im2))
                    !
                  enddo 
               enddo
               !
            endif 
            !
         endif 
         !
      if (verbose>=6) write(out,"('FLmatrixelements_single_iterm/end')") 
    
   end subroutine FLmatrixelements_single_iterm






   !
   ! Finite difference derivatives  
   ! d^n func / d q1^k1 d q2^k2 d q3^k3 .. | at q = x
   ! where (k1,k2,k3,..) = itarget(:)
   ! step(:) - finite differences spacing factor defined for every coordinate 
   ! 
   !
   function FLfinitediffs(itarget,func,ax,astep) result(f)

      integer(ik),intent(in) :: itarget(:)
      real(ark),external      :: func 
      real(ark),intent(in)    :: ax(:)
      real(ark),intent(in)    :: astep(:,:)
      integer(ik)            :: Nmodes
      real(ark)               :: f
      integer(ik) :: imode
      integer(ik) :: isearch(size(Itarget))


      Nmodes = size(Itarget)
      !
      if (verbose>=6) write(out,"(/'FLfinitediffs/start: finite derivatives for k=',30i4)") (itarget(imode),imode=1,min(30,Nmodes))
      !
      ! Itarget defines the total combination of the derivatives 
      ! while in isearch we store the current derivatives level. 
      ! In more details, we start with 
      ! Isearch =  Itarget
      ! then after every single derivative d / d q_i
      ! Isearch(i) = Isearch(i)-1
      !
      Isearch =  Itarget
      !
      ! here we go recursively untill we hit the bottom, i.e. untill isearch = 0 
      !
      f = fdiff(isearch,ax,func) 

      if (verbose>=6) write(out,"('FLfinitediffs/end')") 
      !
   contains 

     recursive function fdiff(isearch,x,func) result (f)

       integer(ik)            :: isearch(Nmodes)
       real(ark),intent(in)   :: x(Nmodes)
       real(ark),external     :: func 
       !
       real(ark)              :: f,h,f1,f2,f3,f4
       real(ark)              :: x1(size(Itarget)),x2(size(Itarget)),x3(size(Itarget)),x4(size(Itarget))
       integer(ik)            :: imode
       !
       !
       ! find first active coordinate 
       imode = 1
       do while(Isearch(imode)==0.and.imode<Nmodes)
         imode = imode+1   
       enddo
       !
       ! if all isearch(:) = 0 we are at the bottom level and go to the func
       !
       if (imode == Nmodes.and.isearch(imode)==0) then
          !
          f = func(x)
          !
       else 
          ! if not all derivatives have been evaluated, we continue to the next lower recursive level
          !
          ! here we go one level lower for the coordinate imode 
          ! 
          Isearch(imode) = Isearch(imode)-1
          !
          ! we use simple central formula for the finite difference derivatives  
          !
          select case(difftype)
          case default
             !
             write (out,"('FLfinitediffs: difftype ',i8,' unknown')") difftype
             stop 'FLfinitediffs - bad difftype'
             !
          case(2)
             !
             x2 = x ; x1 = x 
             x1(imode) = x(imode) - astep(1,imode)
             x2(imode) = x(imode) + astep(2,imode)
             !
             f1 = fdiff(isearch,x1,func)
             f2 = fdiff(isearch,x2,func)
             !
             h = astep(1,imode)+astep(2,imode)
             !
             ! debug !!!
             !
             !if (abs(f2-f1)>1e6) then 
             !   write(out,"('DEBUG SY: f1-f2 is too large',2e20.12)") f1,f2
             !   stop 'DEBUG SY: f1-f2 is too large'
             !endif 
             !
             f = ( f2-f1 )/(h)
             !
             !f_t = 0.5_ark*(f2-f1)
             !f=log(f_t+sqrt(f_t**2+1.0_ark))/astep(1,imode)
             !
          case(4)
             !
             x1 = x ; x2 = x ;  x3 = x ; x4 = x 
             !
             x1(imode) = x(imode) - astep(1,imode)
             x2(imode) = x(imode) + astep(2,imode)
             x3(imode) = x(imode) - astep(1,imode)*2.0_ark
             x4(imode) = x(imode) + astep(2,imode)*2.0_ark
             !
             h = astep(1,imode)+astep(2,imode)
             !
             f1 = fdiff(isearch,x1,func)
             f2 = fdiff(isearch,x2,func)
             f3 = fdiff(isearch,x3,func)
             f4 = fdiff(isearch,x4,func)
             !
             !write (out,"('FLfinitediffs-4: You should check the 4-point equation first!')") 
             !stop 'FLfinitediffs-4: Check the 4-point equation!'
             !
             f = (-f4+8.0_ark*f2+f3-8.0_ark*f1 )/(6.0_ark*h)

             !f = (-f4/12.0_ark+2.0_ark/3.0_ark*f2 & 
             !     +f3/12.0_ark-2.0_ark/3.0_ark*f1 )/h

             !
          end select
          !
          ! here we are back to the current level and the "imode"-coordinate is also back 
          !
          Isearch(imode) = Isearch(imode)+1
          !
       endif 
       !
     end function fdiff
     !
   end function FLfinitediffs


   !subroutine FLpoten4xi(rank,nmodes, ipoint, xi, pot) 


   !
   ! Finite difference derivatives  
   ! d^n func / d q1^k1 d q2^k2 d q3^k3 .. | at q = x
   ! where (k1,k2,k3,..) = itarget(:)
   ! step(:) - finite differences spacing factor defined for every coordinate 
   ! 
   !
   recursive function FLfinitediffs_2d(itarget,ipoint,get_func,x,step) result(f)

      interface
         subroutine get_func(rank,nmodes, ipoint, xi,f)
            use accuracy, only: ik, ark
            use moltype, only: extF
            integer(ik) :: rank,nmodes
            integer(ik) :: ipoint
            real(ark)   :: xi(nmodes)
            real(ark)   :: f(rank)
         end subroutine get_func
      end interface
      !
      integer(ik),intent(in) :: ipoint,itarget(:)
      real(rk),intent(in)    :: x(:)
      real(rk),intent(in)    :: step(:,:)
      integer(ik)            :: Nmodes
      real(ark)               :: f
      integer(ik) :: imode
      integer(ik) :: isearch(size(Itarget))
      real(ark)               :: ax(size(x))
      real(ark)               :: astep(2,size(x))


      Nmodes = size(Itarget)
      !
      !if (verbose>=6) write(out,"(/'FLfinitediffs_2d/start: finite derivatives for k=',30i4)") (itarget(imode),imode=1,min(30,Nmodes))
      !
      ! Itarget defines the total combination of the derivatives 
      ! while in isearch we store the current derivatives level. 
      ! In more details, we start with 
      ! Isearch =  Itarget
      ! then after every single derivative d / d q_i
      ! Isearch(i) = Isearch(i)-1
      !
      Isearch =  Itarget
      !
      ! Transform x and step to high precission
      ax = x
      astep = step
      !
      ! here we go recursively untill we hit the bottom, i.e. untill isearch = 0 
      !
      f = fdiff(Nmodes,isearch,ipoint,ax,astep,get_func) 

      !if (verbose>=6) write(out,"('FLfinitediffs_2d/end')") 
      !
   contains 
     !
     recursive function fdiff(Nmodes,isearch,ipoint,x,astep,get_func) result (f)
       interface
          subroutine get_func(rank,nmodes, ipoint, xi,f)
             use accuracy, only: ik, ark
             use moltype, only: extF
             integer(ik) :: rank,nmodes
             integer(ik) :: ipoint
             real(ark)   :: xi(nmodes)
             real(ark)   :: f(rank)
          end subroutine get_func
       end interface

       integer(ik)            :: Nmodes,isearch(Nmodes),ipoint
       real(ark),intent(in)   :: x(Nmodes)
       real(ark),intent(in)    :: astep(:,:)
       !
       real(ark)              :: f,h,f1,f2,fd(1)
       real(ark)              :: x1(Nmodes),x2(Nmodes)
       integer(ik)            :: imode
       !
       !
       ! find first active coordinate 
       imode = 1
       do while(Isearch(imode)==0.and.imode<Nmodes)
         imode = imode+1   
       enddo
       !
       ! if all isearch(:) = 0 we are at the bottom level and go to the func
       !
       if (imode == Nmodes.and.isearch(imode)==0) then
          !
          call get_func(1,nmodes,ipoint,x,fd)
          !
          f = fd(1)
          !
       else 
          ! if not all derivatives have been evaluated, we continue to the next lower recursive level
          !
          ! here we go one level lower for the coordinate imode 
          ! 
          Isearch(imode) = Isearch(imode)-1
          !
          ! we use simple central formula for the finite difference derivatives  
          !
          x2 = x ; x1 = x 
          x1(imode) = x(imode) - astep(1,imode)
          x2(imode) = x(imode) + astep(2,imode)
          !
          f1 = fdiff(Nmodes,isearch,ipoint,x1,astep,get_func)
          f2 = fdiff(Nmodes,isearch,ipoint,x2,astep,get_func)
          !
          h = astep(1,imode)+astep(2,imode)
          !
          f = ( f2-f1 )/(h)
          !
          ! here we are back to the current level and the "imode"-coordinate is also back 
          !
          Isearch(imode) = Isearch(imode)+1
          !
       endif 
       !
     end function fdiff

     !
   end function FLfinitediffs_2d


   !
   ! Collect all  points needed for the difference derivatives  
   !
   subroutine FLdiffsadresses(itarget,istep,ipoint,ipointaddress)

      integer(ik),intent(in)     :: itarget(:)
      integer(ik),intent(in)     :: istep(:,:)
      integer(ik),intent(out)    :: ipoint
      integer(ik),optional    :: ipointaddress(:,:)

      integer(ik)                :: iaddr(size(Itarget))
      integer(ik)            :: Nmodes
      integer(ik) :: imode,Nmax
      integer(ik) :: isearch(size(Itarget))


      Nmodes = size(Itarget)
      !
      Nmax = trove%NPotOrder
      !
      if (verbose>=6) write(out,"(/'FLdiffsadresses/start: finite derivatives for k=',30i4)")  & 
                                 (itarget(imode),imode=1,min(30,Nmodes))
      !
      ! Itarget defines the total combination of the derivatives 
      ! while in isearch we store the current derivatives level. 
      ! In more details, we start with 
      ! Isearch =  Itarget
      ! then after every single derivative d / d q_i
      ! Isearch(i) = Isearch(i)-1
      !
      Isearch =  0
      !
      ! here we go recursively untill we hit the bottom, i.e. untill isearch = 0 
      !
      iaddr = 0 
      !
      ipoint = 0 
      !
      call next_mode(isearch,1,iaddr,ipoint) 

      if (verbose>=6) write(out,"('FLdiffsadresses/end')") 
      !
   contains 

     recursive subroutine next_mode(isearch,imode,iaddr,ipoint)

       integer(ik),intent(inout)  :: isearch(Nmodes)
       integer(ik),intent(in) :: iaddr(Nmodes)
       integer(ik),intent(inout) :: ipoint
       !
       integer(ik)            :: imode,iaddr1(size(Itarget)),iaddr2(size(Itarget))
       integer                :: iaddr_t(size(Itarget)),iterm,jmode,N_t
       !
       if (imode>trove%Nmodes) return 
       !
       do iterm =-istep(1,imode)*itarget(imode),istep(2,imode)*itarget(imode)
        !
        isearch(imode) = iterm 
        !
        N_t = sum(abs(isearch))
        !
        if (N_t<=Nmax) then 
          !
          if (imode==trove%Nmodes) then 
            !
            ipoint = ipoint + 1 
            !
            if (present(ipointaddress)) then 
              !
              ipointaddress(ipoint,:) = isearch(:)
              !
            endif 
            !
          else
            !
            call next_mode(isearch,imode+1,iaddr,ipoint)
            !
          endif 
          !  
        endif 
        !
       enddo
       !
       isearch(imode) = 0  
       !
     end subroutine next_mode
     !
   end subroutine FLdiffsadresses





   !
   ! Collect all  points needed for the difference derivatives  
   !
   subroutine FLdiffsadresses_roman(itarget,istep,ipoint,par,ipointaddress)

      integer(ik),intent(in)     :: itarget(trove%Nmodes)
      integer(ik),intent(in)     :: istep(2,trove%Nmodes)
      integer(ik),intent(out)    :: ipoint
      integer(ik),intent(in)         :: par(0:trove%Nmodes+1) ! Parity of each mode in derivatives calculation
      integer(ik),optional       :: ipointaddress(:,:)

      integer(ik)                :: iaddr(trove%Nmodes)
      integer(ik) :: imode
      integer(ik) :: isearch(trove%Nmodes)


      !Nmodes = size(Itarget)
      !
      !Nmax = trove%NPotOrder
      !
      if (verbose>=6) write(out,"(/'FLdiffsadresses/start: finite derivatives for k=',30i4)")  & 
                                 (itarget(imode),imode=1,min(30,trove%Nmodes))
      !
      ! Itarget defines the total combination of the derivatives 
      ! while in isearch we store the current derivatives level. 
      ! In more details, we start with 
      ! Isearch =  Itarget
      ! then after every single derivative d / d q_i
      ! Isearch(i) = Isearch(i)-1
      !
      Isearch =  0
      !
      ! here we go recursively untill we hit the bottom, i.e. untill isearch = 0 
      !
      iaddr = 0 
      !
      ipoint = 0 
      !
      call next_mode(itarget,istep,par,isearch,1,iaddr,ipoint) 

      if (verbose>=6) write(out,"('FLdiffsadresses/end')") 
      !
     contains 
     !
     recursive subroutine next_mode(itarget,istep,par,isearch,imode,iaddr,ipoint)

       integer(ik),intent(in)     :: itarget(trove%Nmodes)
       integer(ik),intent(in)     :: istep(2,trove%Nmodes)
       integer(ik),intent(in)     :: par(0:trove%Nmodes+1)  ! Parity of each mode in derivatives calculation
       integer(ik),intent(inout)  :: isearch(trove%Nmodes)
       integer(ik),intent(in)     :: iaddr(trove%Nmodes)
       integer(ik),intent(inout)  :: ipoint
       !
       integer(ik)            :: imode,iaddr1(trove%Nmodes),iaddr2(trove%Nmodes)
       integer(ik)            :: iaddr_t(trove%Nmodes),iterm,jmode,N_t
       integer(ik)            :: iterm1,iterm2,step_iterm
       !
       if (imode>trove%Nmodes) return 
       !
       if (par(imode)==-1) then
        !
        iterm1 = -istep(1,imode)*itarget(imode)
        !
        iterm2 =  istep(2,imode)*itarget(imode)
        !
        step_iterm = 1
        !
       else
        !
        iterm1 = -itarget(imode) + mod(itarget(imode)+par(imode),2)
        !
        iterm2 =  itarget(imode) - mod(itarget(imode)+par(imode),2)
        !
        step_iterm = 2
        !
       endif
       !
       do iterm = iterm1,iterm2,step_iterm
        !
        isearch(imode) = iterm 
        !
        N_t = sum(abs(isearch))
        !
        if (N_t<=trove%MaxOrder) then 
          !
          if (imode==trove%Nmodes) then 
            !
            ipoint = ipoint + 1 
            !
            if (present(ipointaddress)) then 
              !
              ipointaddress(ipoint,:) = isearch(:)
              !
            endif 
            !
          else
            !
            call next_mode(itarget,istep,par,isearch,imode+1,iaddr,ipoint)
            !
          endif 
          !  
        endif 
        !
       enddo
       !
       isearch(imode) = 0  
       !
     end subroutine next_mode
     !
   end subroutine FLdiffsadresses_roman



  !
   ! Finite difference derivatives  
   ! d^n func / d q1^k1 d q2^k2 d q3^k3 .. | at q = x
   ! where (k1,k2,k3,..) = itarget(:)
   ! step(:) - finite differences spacing factor defined for every coordinate 
   ! 
   !
   function FLfinitediffs_imode(itarget,func,x,step,kmode) result(f)

      integer(ik),intent(in) :: itarget(:)
      real(ark),external      :: func 
      real(rk),intent(in)    :: x(:)
      real(rk),intent(in)    :: step(:,:)
      integer(ik),intent(in) :: kmode 
      integer(ik)            :: Nmodes
      real(ark)               :: f
      integer(ik) :: imode
      integer(ik) :: isearch(size(Itarget))
      real(ark)               :: ax(size(x))
      real(ark)               :: astep(2,size(x))


      Nmodes = size(Itarget)
      !
      if (verbose>=6) write(out,"(/'FLfinitediffs/start: finite derivatives for k=',30i4)") (itarget(imode),imode=1,min(30,Nmodes))
      !
      ! Itarget defines the total combination of the derivatives 
      ! while in isearch we store the current derivatives level. 
      ! In more details, we start with 
      ! Isearch =  Itarget
      ! then after every single derivative d / d q_i
      ! Isearch(i) = Isearch(i)-1
      !
      Isearch =  Itarget
      ! Transform x and step to high precission
      ax = x
      astep = step
      !
      ! here we go recursively untill we hit the bottom, i.e. untill isearch = 0 
      !
      f = fdiff(isearch,ax,func) 

      if (verbose>=6) write(out,"('FLfinitediffs/end')") 
      !
   contains 

     recursive function fdiff(isearch,x,func) result (f)

       integer(ik)            :: isearch(Nmodes)
       real(ark),intent(in)   :: x(Nmodes)
       real(ark),external     :: func 
       !
       real(ark)              :: f,h,f1,f2,f3,f4
       real(ark)              :: x1(size(Itarget)),x2(size(Itarget)),x3(size(Itarget)),x4(size(Itarget))
       integer(ik)            :: imode
       !
       !
       ! find first active coordinate 
       imode = 1
       do while(Isearch(imode)==0.and.imode<Nmodes)
         imode = imode+1   
       enddo
       !
       ! if all isearch(:) = 0 we are at the bottom level and go to the func
       !
       if (imode == Nmodes.and.isearch(imode)==0) then
          !
          f = func(x,kmode)
          !
       else 
          ! if not all derivatives have been evaluated, we continue to the next lower recursive level
          !
          ! here we go one level lower for the coordinate imode 
          ! 
          Isearch(imode) = Isearch(imode)-1
          !
          ! we use simple central formula for the finite difference derivatives  
          !
          select case(difftype)
          case default
             !
             write (out,"('FLfinitediffs: difftype ',i8,' unknown')") difftype
             stop 'FLfinitediffs - bad difftype'
             !
          case(2)
             !
             x2 = x ; x1 = x 
             x1(imode) = x(imode) - astep(1,imode)
             x2(imode) = x(imode) + astep(2,imode)
             !
             f1 = fdiff(isearch,x1,func)
             f2 = fdiff(isearch,x2,func)
             !
             h = astep(1,imode)+astep(2,imode)
             !
             f = ( f2-f1 )/(h)
             !
             !f_t = 0.5_ark*(f2-f1)
             !f=log(f_t+sqrt(f_t**2+1.0_ark))/astep(1,imode)
             !
          case(4)
             !
             x1 = x ; x2 = x ;  x3 = x ; x4 = x 
             !
             x1(imode) = x(imode) - astep(1,imode)
             x2(imode) = x(imode) + astep(2,imode)
             x3(imode) = x(imode) - astep(1,imode)*2.0_ark
             x4(imode) = x(imode) + astep(2,imode)*2.0_ark
             !
             h = astep(1,imode)+astep(2,imode)
             !
             f1 = fdiff(isearch,x1,func)
             f2 = fdiff(isearch,x2,func)
             f3 = fdiff(isearch,x3,func)
             f4 = fdiff(isearch,x4,func)
             !
             !write (out,"('FLfinitediffs-4: You should check the 4-point equation first!')") 
             !stop 'FLfinitediffs-4: Check the 4-point equation!'
             !
             f = (-f4+8.0_ark*f2+f3-8.0_ark*f1 )/(6.0_ark*h)

             !f = (-f4/12.0_ark+2.0_ark/3.0_ark*f2 & 
             !     +f3/12.0_ark-2.0_ark/3.0_ark*f1 )/h

             !
          end select
          !
          ! here we are back to the current level and the "imode"-coordinate is also back 
          !
          Isearch(imode) = Isearch(imode)+1
          !
       endif 
       !
     end function fdiff
     !
   end function FLfinitediffs_imode




   !
   ! Finite difference derivatives of a vector function func(:)
   ! d^n func / d q1^k1 d q2^k2 d q3^k3 .. | at q = x
   ! where (k1,k2,k3,..) = itarget(:)
   ! step(:) - finite differences spacing factor defined for every coordinate 
   ! 
   !
   subroutine FLfinitediffs_vect(job_str,Nmodes,Nsize_f,itarget,ax,astep,irho,f)
      !
      character(len=cl),intent(in)  :: job_str
      integer(ik),intent(in) :: Nmodes,Nsize_f
      integer(ik),intent(in) :: itarget(Nmodes)
      !
      real(ark),intent(in)    :: ax(Nmodes)
      real(ark),intent(in)    :: astep(2,Nmodes)
      integer,intent(in)     :: irho
      real(ark),dimension(Nsize_f) :: f
      integer(ik) :: imode
      integer(ik) :: isearch(Nmodes)
      character(len=cl) :: my_fmt !format for I/O specification
      !
      write(my_fmt,'(a,i0,a)') "(/a,",Nmodes,"i4)"
      !
      if (verbose>=6) write(out,my_fmt) 'Finite deriv-s for k=',(itarget(imode),imode=1,Nmodes)
      !
      ! Itarget defines the total combination of the derivatives 
      ! while in isearch we store the current derivatives level. 
      ! In more details, we start with 
      ! Isearch =  Itarget
      ! then after every single derivative d / d q_i
      ! Isearch(i) = Isearch(i)-1
      !
      Isearch =  Itarget
      ! Transform x and step to high precission
      !
      ! here we go recursively untill we hit the bottom, i.e. untill isearch = 0 
      !
      f =  fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,ax) 

      if (verbose>=6) write(out,"('FLfinitediffs_vect/end')") 
      !
      !contains 
      !
   end subroutine FLfinitediffs_vect


   !
   ! Finite difference derivatives of a vector function func(:)
   ! d^n func / d q1^k1 d q2^k2 d q3^k3 .. | at q = x
   ! where (k1,k2,k3,..) = itarget(:)
   ! step(:) - finite differences spacing factor defined for every coordinate 
   ! 
   !
   recursive function FLvect_finitediffs(job_str,Nsize_f,itarget,ax,astep,irho) result (f)

      character(len=cl),intent(in)  :: job_str
      integer(ik),intent(in) :: Nsize_f
      integer(ik),intent(in) :: itarget(:)
      !real(ark),external :: func

      real(ark),intent(in)    :: ax(:)
      real(ark),intent(in)    :: astep(:,:)
      integer,intent(in)     :: irho
      real(ark),dimension(Nsize_f) :: f
      integer(ik) :: imode
      integer(ik) :: isearch(size(Itarget))
      integer(ik)            :: Nmodes
      character(len=cl) :: my_fmt !format for I/O specification
      !
      Nmodes = size(Itarget)
      !
      write(my_fmt,'(a,i0,a)') "(/a,",Nmodes,"i4)"
      !
      if (verbose>=6) write(out,my_fmt) 'Finite deriv-s for k=',(itarget(imode),imode=1,Nmodes)
      !
      ! Itarget defines the total combination of the derivatives 
      ! while in isearch we store the current derivatives level. 
      ! In more details, we start with 
      ! Isearch =  Itarget
      ! then after every single derivative d / d q_i
      ! Isearch(i) = Isearch(i)-1
      !
      Isearch =  Itarget
      ! Transform x and step to high precission
      !
      ! here we go recursively untill we hit the bottom, i.e. untill isearch = 0 
      !
      f =  fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,ax) 

      if (verbose>=6) write(out,"('FLvect_finitediffs/end')") 
      !
      !contains 
      !
   end function FLvect_finitediffs


   recursive function fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,x) result (f)

     integer(ik),intent(in) :: irho,Nmodes,Nsize_f
     integer(ik)            :: isearch(Nmodes)
     real(ark),intent(in)   :: x(Nmodes)
     integer(ik),intent(in) :: itarget(Nmodes)
     character(len=cl),intent(in)  :: job_str
     real(ark),intent(in)  :: astep(2,Nmodes)
     !real(ark),external    ::  func
     real(ark),dimension(Nsize_f) :: f
     !
     real(ark)             :: f1(Nsize_f),f2(Nsize_f),f3(Nsize_f),f4(Nsize_f),h
     real(ark)             :: x1(Nmodes),x2(Nmodes),x3(Nmodes),x4(Nmodes)
     integer(ik)           :: imode
     logical     :: dir
     !
     ! find first active coordinate 
     imode = 1
     do while(Isearch(imode)==0.and.imode<Nmodes)
       imode = imode+1   
     enddo
     !
     ! if all isearch(:) = 0 we are at the bottom level and go to the calc_local2cartesian_vec
     !
     if (imode == Nmodes.and.isearch(imode)==0) then
        !
        select case(trim(job_str))
        case default
           !
           write (out,"('FLfinitediffs: job-type ',a,' unknown')") job_str
           stop 'FLfinitediffs - bad job-type'
           !
        case('local2cartesian')
           !
           f = calc_local2cartesian_vec(Nsize_f,1,x,irho)
           !
        case('cartesian2local')
           !
           f = calc_cartesian2local_vec(Nsize_f,x,irho)
           !
        case('xi2cartesian')
           !
           f = FLcalc_xi2cartesian_vec(Nsize_f,1,x,irho)
           !
        case('xi2dcartesian_drho')
           !
           f = calc_xi2dcartesian_drho_vec(Nsize_f,1,x,irho)
           !
        case('calc_xi2drho_dxna_vec')
           !
           f = calc_xi2drho_dxna_vec(Nsize_f,1,x,irho)
           !
        case('local2chi')
           !
           dir = .true.
           !
           f = MLcoordinate_transform_func(x,Nsize_f,dir)
           !
        !case('s_vib_s_rot_dvr')
        !   !
        !   f = s_vib_s_rot_dvr_vec(x,irho,Nsize_f)
           !
        end select
        !
     else 
        ! if not all derivatives have been evaluated, we continue to the next lower recursive level
        !
        ! here we go one level lower for the coordinate imode 
        !
        Isearch(imode) = Isearch(imode)-1
        !
        ! we use simple central formula for the finite difference derivatives  
        !
        select case(difftype)
        case default
           !
           write (out,"('FLfinitediffs: difftype ',i8,' unknown')") difftype
           stop 'FLfinitediffs - bad difftype'
           !
        case(2)
           !
           x1 = x ; x2 = x
           !
           x1(imode) = x(imode) - astep(1,imode)
           x2(imode) = x(imode) + astep(2,imode)
           !
           h = astep(1,imode)+astep(2,imode)
           !
           f1 = fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,x1) 
           f2 = fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,x2) 
           !
           f(:) = ( f2(:)-f1(:) )/(h)
           !
        case(4)
           !
           x1 = x ; x2 = x ;  x3 = x ; x4 = x 
           !
           x1(imode) = x(imode) - astep(1,imode)
           x2(imode) = x(imode) + astep(2,imode)
           x3(imode) = x(imode) - astep(1,imode)*2.0_ark
           x4(imode) = x(imode) + astep(2,imode)*2.0_ark
           !
           h = astep(1,imode)+astep(2,imode)
           !
           f1 = fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,x1)
           f2 = fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,x2)
           f3 = fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,x3)
           f4 = fdiff_vec(irho,Nmodes,Nsize_f,job_str,isearch,itarget,astep,x4)
           !
           !write (out,"('FLvect_finitediffs-4: You should check the 4-point equation first!')") 
           !stop 'FLvect_finitediffs-4: Check the 4-point equation!'
           !
           !f(:) = (-f4(:)/12.0_ark+2.0_ark/3.0_ark*f2(:)+f3(:)/12.0_ark-2.0_ark/3.0_ark*f1(:) )/h*2.0_ark
           !
           f(:) = (-f4(:)+8.0_ark*f2(:)+f3(:)-8.0_ark*f1(:) )/(6.0_ark*h)
           !
        end select
        !
        ! here we are back to the current level and the "imode"-coordinate is also back 
        !
        Isearch(imode) = Isearch(imode)+1
        !
     endif 
     !
   end function fdiff_vec

!
!
! Potential function in terms of the normal (linearized) coordinates 
!
   function poten_normal(xi) result (f) 

     real(ark),intent(in)    :: xi(trove%Nmodes)
     !
     real(ark)               :: f,r(trove%Ncoords),chi(trove%Nmodes)
     real(ark)               :: r_na(trove%Natoms,3)
     integer(ik)             :: n1,ix,i,irho
     !
     if (verbose>=6) write(out,"(/'poten_normal/start')") 
     !
     ! Build up geometrically defined coordinates from the cartesian coordinates, 
     ! which are defined in terms of the normal (linearized) coordinates 
     ! 
     !irho =     nint( ( molec%chi_eq(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik )
     !
     irho =  FLirho
     !
     do i = 1,trove%Nmodes_e
       chi(i) = MLcoord_invert(xi,2,i) 
     enddo
     !
     !chi = xi
     !
     do n1 = 1,trove%Natoms
        do ix = 1,3
          !
          r_na(n1,ix) = trove%b0(n1,ix,irho) + sum( ( trove%Amatrho(n1,ix,1:trove%Nmodes_e,irho) )* &
                                               (chi(1:trove%Nmodes_e)-trove%chi_ref(1:trove%Nmodes_e,irho)) )
          !
        enddo
     enddo
     !
     call FLfromcartesian2local(r_na,r)
     !
     ! The potential energy function value in terms of the GDC "r"
     f = MLpotenfunc(r,r_na)
     !
     if (verbose>=5) then
         if (irho==0) then
            write(out,"('r   = ',20f18.8)") (r(n1),n1=1,min(20,trove%Ncoords))
            write(out,"('xi  = ',20f18.8)") (xi(n1),n1=1,min(20,trove%Nmodes))
            write(out,"('chi = ',20f18.8)") (chi(n1),n1=1,min(20,trove%Nmodes))
            write(out,"('rnax = ',20f18.8)") (r_na(n1,1),n1=1,min(20,trove%Natoms))
            write(out,"('rnay = ',20f18.8)") (r_na(n1,2),n1=1,min(20,trove%Natoms))
            write(out,"('rnaz = ',20f18.8)") (r_na(n1,3),n1=1,min(20,trove%Natoms))
            write(out,"('f = ',f18.8)") f
         endif
     endif
     !
     !
     if (verbose>=6) write(out,"('poten_normal/end')") 
     !
   end function poten_normal


!
!
! Potential function in terms of the normal (linearized) coordinates 
!
   recursive function FLpoten_linearized_dchi(dchi,irho) result (f) 

     real(ark),intent(in)    :: dchi(trove%Nmodes)
     integer(ik),intent(in)  :: irho
     !
     real(ark)               :: f,r(trove%Ncoords)
     real(ark)               :: r_na(trove%Natoms,3)
     integer(ik)             :: n1,ix,i
     !
     if (verbose>=6) write(out,"(/'FLpoten_linearized_dchi/start')") 
     !
     ! Build up geometrically defined coordinates from the cartesian coordinates, 
     ! which are defined in terms of the normal (linearized) coordinates 
     !
     do n1 = 1,trove%Natoms
        do ix = 1,3
          !
          r_na(n1,ix) = trove%b0(n1,ix,irho) + sum( trove%Amatrho(n1,ix,1:trove%Nmodes_e,irho)*dchi(1:trove%Nmodes_e) )
          !
        enddo
     enddo
     !
     call FLfromcartesian2local(r_na,r)
     !
     ! The potential energy function value in terms of the GDC "r"
     f = MLpotenfunc(r,r_na)
     !
     if (verbose>=6) write(out,"('FLpoten_linearized_dchi/end')") 
     !
   end function FLpoten_linearized_dchi



!
!
! Potential function in terms of the normal (linearized) coordinates 
!
   recursive function FLpoten_linearized(xi,irho) result (f) 

     real(ark),intent(in)    :: xi(trove%Nmodes)
     integer(ik),intent(in)  :: irho
     !
     real(ark)               :: f,r(trove%Ncoords),chi(trove%Nmodes)
     real(ark)               :: r_na(trove%Natoms,3)
     integer(ik)             :: n1,ix,i
     !
     if (verbose>=6) write(out,"(/'FLpoten_linearized/start')") 
     !
     ! Build up geometrically defined coordinates from the cartesian coordinates, 
     ! which are defined in terms of the normal (linearized) coordinates 
     !
     !chi(1:trove%Nmodes_e) = trove%chi_ref(1:trove%Nmodes_e,irho)
     !
     do i = 1,trove%Nmodes_e
       chi(i) = MLcoord_invert(xi,2,i) 
     enddo
     !
     !chi = xi
     !
     do n1 = 1,trove%Natoms
        do ix = 1,3
          !
          r_na(n1,ix) = trove%b0(n1,ix,irho) + sum( ( trove%Amatrho(n1,ix,1:trove%Nmodes_e,irho) )* &
                                                    ( chi(1:trove%Nmodes_e)-trove%chi_ref(1:trove%Nmodes_e,irho)) )
          !
        enddo
     enddo
     !
     call FLfromcartesian2local(r_na,r)
     !
     ! The potential energy function value in terms of the GDC "r"
     f = MLpotenfunc(r,r_na)
     !
     if (verbose>=6) write(out,"('FLpoten_linearized/end')") 
     !
   end function FLpoten_linearized

!
!
! Potential function in terms of the xi coordinates = xi(geometrically defined)
!
   function poten_xi(xi) result (f) 

     real(ark),intent(in)    :: xi(trove%Nmodes)
     !
     real(ark)               :: f,r(trove%Ncoords),chi(trove%Nmodes),chi_eq(trove%Nmodes),r_(trove%Ncoords)
     real(ark)               :: r_na(trove%Natoms,3)
     !
     integer(ik)             :: i,pm,Ncoords
     logical                 :: dir
     character(len=cl) :: my_fmt !format for I/O specification

     if (verbose>=6) write(out,"(/'poten_xi/start')") 
     do i = 1,trove%Nmodes
        chi(i) = MLcoord_invert(xi,2,i)
     enddo
     !
     Ncoords = trove%Ncoords
     !
     chi_eq(:) = trove%chi_ref(:,FLirho)
     !
     ! reconstruct the TROVE-ccordinates from Zmat curvelinear coordinates 
     !
     dir = .false.
     !
     r = MLcoordinate_transform_func(chi,size(r),dir)
     !
     pm = 1
     !
     call MLfromlocal2cartesian(pm,r,r_na)
     !
     f = MLpotenfunc(r,r_na)
     !
     call FLfromcartesian2local(r_na,r_)
     !
     if ( any( abs( r(:)-r_(:) )>sqrt(small_) ) ) then 
       !
       write(my_fmt,'(a,i0,a)') "(4x,",Ncoords,"f18.6)"
       !
       write(out,'("poten_xi: Error in MLfromlocal2cartesian, r /= r_:")')
       write(out,my_fmt) r(:)
       write(out,my_fmt) r_(:)
       stop "poten_xi: Error in MLfromlocal2cartesian, r /= r_"
       !
     endif 
     !
     if (verbose>=6) write(out,"('poten_xi/end')") 
     !
   end function poten_xi


!
!
! Potential function in terms of the xi coordinates = xi(geometrically defined)
!
   function poten_chi(chi) result (f) 

     real(ark),intent(in)    :: chi(trove%Nmodes)
     real(ark)               :: r_na(trove%Natoms,3)
     !
     real(ark)               :: f,r(trove%Ncoords),chi_(trove%Nmodes),r_(trove%Ncoords)
     logical                 :: dir
     integer(ik)             :: pm = 1,Ncoords,i
     !
     if (verbose>=6) write(out,"(/'poten_chi/start')") 
     !
     Ncoords = trove%Ncoords
     !
     !chi_eq(:) = trove%chi_ref(:,FLirho)
     !
     dir = .false.
     r = MLcoordinate_transform_func(chi,size(r),dir)
     !
     ! The potential energy function value in terms of the GDC "r"
     !
     call MLfromlocal2cartesian(pm,r,r_na)
     !
     f = MLpotenfunc(r,r_na)
     !
     call FLfromcartesian2local(r_na,r_)
     dir = .true.
     chi_ = MLcoordinate_transform_func(r_,size(r_),dir)
     !
     do i = 1,trove%Nmodes
       !
       if ( abs( chi(i)-chi_(i) )>10.0*sqrt(small_).and.abs( chi(i)-chi_(i) )-2.0_ark*pi> 10.0*sqrt(small_)) then
         !
         write(out,'("poten_chi: chi /= chi_: ",i5,4x,2f18.6,1x,e10.3)') i,chi(i),chi_(i),chi(i)-chi_(i)
         !stop "poten_chi: Error in MLfromlocal2cartesian, r /= r_"
         !
        endif
        ! 
     enddo
     !
     if (verbose>=6) write(out,"('poten_chi/end')") 
     !
   end function poten_chi

!
!
! Potential function in terms of the local coordinates (geometrically defined)
!
   function poten_local(r) result (f) 

     real(ark),intent(in)    :: r(trove%Ncoords)
     !
     real(ark)               :: f
     real(ark)               :: r_na(trove%Natoms,3)
     integer(ik)             :: pm = 1
     !
     if (verbose>=6) write(out,"(/'poten_local/start')") 
     !
     ! The potential energy function value in terms of the GDC "r"
     !
     call MLfromlocal2cartesian(pm,r,r_na)
     !
     f = MLpotenfunc(r,r_na)
     !
     if (verbose>=6) write(out,"('poten_local/end')") 
     !
   end function poten_local

!
!
! Zero order energy calculated from the 1d zero order Hamiltonians
!
   function FLenergy_zero(nu) result (e_t) 

     integer(ik),intent(in)    :: nu(:)
     !
     real(rk)                  :: e_t
     integer(ik)               :: imode,ibstype,i,nu_t(0:trove%Nmodes),j,k,jk
     type(Basis1DT),pointer    :: bs
     !
     if (verbose>=6) write(out,"(/'FLenergy_zero/start')") 
       !
       nu_t(:) = nu(:)
       !
       e_t = 0 
       !
       if (bset_initialized) then 
           !
           if (FLrotation) then 
             !
             bs => bset%rot
             jk = nu_t(0) 
             !
             e_t = e_t + bs%ener0(jk)
             !
           endif 
           !
           do ibstype = 1,bset%n_bset1D_max
              !
              bs => bset%bs1D(ibstype)
              !
              do i = 1,bs%imodes
                 !
                 imode = bs%mode(i)
                 !
                 e_t = e_t + bs%ener0(nu_t(imode))
                 !
              enddo 
              ! 
           enddo
           !
       endif 
       !
     if (verbose>=6) write(out,"('FLenergy_zero/end')") 
     !
   end function FLenergy_zero


!
! Vector function to convert xi into cartesian coordinates 
!
   recursive function FLcalc_xi2cartesian_vec(Nsize_f,itype,xi,irho) result (r_na) 

     integer(ik),intent(in)  :: Nsize_f,itype
     real(ark),intent(in)    :: xi(trove%Nmodes)
     integer(ik),intent(in)  :: irho
     !
     real(ark),dimension(Nsize_f)  :: r_na
     !
     real(ark) :: chi(trove%Nmodes)
     integer(ik) :: i,numvar,x0,n0

     if (verbose>=8) write(out,"(/'FLcalc_xi2cartesian_vec/start')") 
     !
     do i = 1,trove%Nmodes
        !
        chi(i) = MLcoord_invert(xi,itype,i)
        !
     enddo
     !
     numvar = 0 
     do n0 = 1,trove%Natoms
        do x0 = 1,3
           numvar= numvar+1
           r_na(numvar) = trove%b0(n0,x0,irho) + & 
                          sum( trove%Amatrho(n0,x0,1:trove%Nmodes_e,irho)*&
                             (chi(1:trove%Nmodes_e)-trove%chi_ref(1:trove%Nmodes_e,irho)) )
           !
        enddo
     enddo
     !
     if (verbose>=8) write(out,"('FLcalc_xi2cartesian_vec/end')") 
     !
   end function FLcalc_xi2cartesian_vec

!
! Vector function to convert xi into d cartesian / d rho 
!
   function calc_xi2dcartesian_drho_vec(Nsize_f,itype,xi,irho) result (r_na) 

     integer(ik),intent(in)  :: Nsize_f,itype
     real(ark),intent(in)    :: xi(trove%Nmodes)
     integer(ik),intent(in)  :: irho
     !
     real(ark),dimension(Nsize_f)  :: r_na
     !
     real(ark) :: chi(trove%Nmodes)
     integer(ik) :: i,numvar,x0,n0

     if (verbose>=8) write(out,"(/'calc_xi2dcartesian_drho_vec/start')") 
     !
     do i = 1,trove%Nmodes
        !
        chi(i) = MLcoord_invert(xi,itype,i)
        !
     enddo
     !
     numvar = 0 
     do n0 = 1,trove%Natoms
        do x0 = 1,3
           numvar= numvar+1
           r_na(numvar) = trove%db0(n0,x0,irho,1) + &
                          sum( trove%dAmatrho(n0,x0,1:trove%Nmodes_e,irho,1)*&
                             (chi(1:trove%Nmodes_e)-trove%chi_ref(1:trove%Nmodes_e,irho)) )
        enddo
     enddo
     !
     if (verbose>=8) write(out,"('calc_xi2dcartesian_drho_vec/end')") 
     !
   end function calc_xi2dcartesian_drho_vec



!
! Vector function to convert local to into d local / d x_Na
!
   function calc_xi2drho_dxna_vec(Nsize_f,itype,xi,irho) result (drho) 

     integer(ik),intent(in)  :: Nsize_f,itype
     real(ark),intent(in)    :: xi(trove%Nmodes)
     integer(ik),intent(in)  :: irho
     !
     real(ark),dimension(Nsize_f)  :: drho
     !
     real(ark)   :: chi(trove%Nmodes),r_na0(Nsize_f),r_na(Nsize_f),dr,r(trove%Ncoords),chi_r(trove%Nmodes),chi_l(trove%Nmodes)
     integer(ik) :: i,numvar,x0,n0
     logical     :: dir

     !if (verbose>=8) write(out,"(/'calc_xi2drho_dxna_vec/start')") 
     !
     dr = trove%fdstep(1)
     !
     do i = 1,trove%Nmodes
        !
        chi(i) = MLcoord_invert(xi,itype,i)
        !
     enddo
     !
     numvar = 0 
     do n0 = 1,trove%Natoms
        do x0 = 1,3
           numvar= numvar+1
           r_na0(numvar) = trove%db0(n0,x0,irho,1) + &
                          sum( trove%dAmatrho(n0,x0,1:trove%Nmodes_e,irho,1)*&
                             (chi(1:trove%Nmodes_e)-trove%chi_ref(1:trove%Nmodes_e,irho)) )
        enddo
     enddo
     !
     r_na = r_na0
     !
     dir = .true.
     !
     numvar = 0
     do n0 = 1,trove%Natoms
        do x0 = 1,3
          numvar = numvar +1
          !
          r_na(numvar) = r_na0(numvar)+dr
          !
          call FLfromcartesian2local(r_na,r)
          !
          chi_r = MLcoordinate_transform_func(r,size(chi_r),dir)
          !
          r_na(numvar) = r_na0(numvar)-dr
          !
          call FLfromcartesian2local(r_na,r)
          !
          chi_l = MLcoordinate_transform_func(r,size(chi_l),dir)
          !
          r_na(numvar) = r_na0(numvar)
          !
          drho(numvar) = 0.5_ark*( chi_r(trove%Nmodes)-chi_l(trove%Nmodes) )/dr
          !
        enddo
     enddo
     !
     !if (verbose>=8) write(out,"('calc_xi2drho_dxna_vec/end')") 
     !
   end function calc_xi2drho_dxna_vec



!
!
! Vector function to convert local into cartesian coordinates 
!
   function calc_cartesian2local_vec(Nsize_f,r_na,irho) result (r) 

     integer(ik),intent(in)  :: Nsize_f

     real(ark),intent(in)    :: r_na(:)
     !
     integer(ik),intent(in)  :: irho
     !
     real(ark),dimension(Nsize_f)  :: r
     !
     real(ark)   :: r_na_t(trove%Natoms,3)
     integer(ik) :: ieq,ix,n1
     !
     ieq = 0
     do n1 = 1,trove%Natoms
        do ix = 1,3
          ieq = ieq +1 
          r_na_t(n1,ix) = r_na(ieq)
        enddo
     enddo
     !
     call FLfromcartesian2local(r_na_t,r)
     !
     if (verbose>=6) write(out,"('calc_local2cartesian_vec/end')") 
     !
   end function calc_cartesian2local_vec



!
!
! Vector function to convert local into cartesian coordinates 
!
   function calc_local2cartesian_vec(Nsize_f,itype,xi,irho) result (r_na) 

     integer(ik),intent(in)  :: Nsize_f,itype
     real(ark),intent(in)    :: xi(trove%Nmodes)
     integer(ik),intent(in)  :: irho
     !
     real(ark),dimension(Nsize_f)  :: r_na,r_na0
     !
     real(ark)                     :: r(trove%Ncoords)
     !
     real(ark) :: chi(trove%Nmodes),chi_eq(trove%Nmodes),r_na_t(trove%Natoms,3)
     logical   :: dir
     integer(ik) :: i,ieq,ix,n1

     if (verbose>=6) write(out,"(/'calc_local2cartesian_vec/start')") 
     !
     do i = 1,trove%Nmodes
        !
        chi(i) = MLcoord_invert(xi,itype,i)
        !
     enddo
     !
     dir = .false.
     !
     r = MLcoordinate_transform_func(chi,size(r),dir)
     !
     chi_eq(:) = trove%chi_ref(:,irho)
     !
     call fromlocal2cartesian(r,r_na_t)
     !
     ieq = 0
     do n1 = 1,trove%Natoms
        do ix = 1,3
          ieq = ieq +1 
          r_na0(ieq) = r_na_t(n1,ix)
        enddo
     enddo
     !cartesian_stored = r_na
     !
     !chi_eq = trove%chi_ref
     !
     r_na = r_na0
     !
     call from_local2cartesian_by_fit(chi,chi_eq,r_na,r_na0)
     
     !
     if (verbose>=6) write(out,"('calc_local2cartesian_vec/end')") 
     !
   end function calc_local2cartesian_vec



   subroutine fromlocal2cartesian(r,cartesian)

    real(ark),intent(in)    :: r(trove%Ncoords)
    real(ark),intent(out)   :: cartesian(trove%Natoms,3)

    real(ark)   :: x(trove%Natoms,3),y(trove%Natoms,3),chi(trove%Nmodes)
    real(ark)   :: n1(3),n2(3),n3(3),c1(3),c2(3),c3(3),sign_t
    !
    real(ark)   :: rot1(3,3),rot2(3,3), local(trove%Ncoords),CM_shift,dlocal(trove%Ncoords)
    real(ark)   :: a_t
    !
    real(ark)   :: rbond,alpha13,alpha23,alpha12,tau_2,tau,delta,theta0,phi0,chi0,theta2,r1

    integer(ik) ::  iangle,idihedral,iatom,j
    integer(ik) ::  p0,p1,p2,p3,p4,ix,jx,kx,n0,irho
    !
    logical     :: dir
     !
     if (verbose>=6) write(out,"(/'fromlocal2cartesian/start')") 
     !
     ! reference irho position 
     !
     !irho = mod(nint( ( molec%chi_eq(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
     !irho =     nint( ( molec%chi_eq(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik )
     irho = FLirho
     !if (irho<0.or.irho>trove%npoints) irho = mod(irho,trove%npoints)
     !
     dir = .true.
     !
     ! Reconstruct chi from r
     !
     chi = MLcoordinate_transform_func(r,size(chi),dir)
     !
     ! Construct the approximate cartesian positions from the first derivatives Amatrho
     ! as a reference 
     !
     do n0 = 1,trove%Natoms
        do ix = 1,3
           !
           y(n0,ix) = trove%b0(n0,ix,irho) + & 
                      sum( trove%Amatrho(n0,ix,1:trove%Nmodes_e,irho)*&
                         ( chi(1:trove%Nmodes_e)-trove%chi_ref(1:trove%Nmodes_e,irho) ) )
           !
        enddo
     enddo
     !
     ! cartesian coordinates
     !
     x = 0
     !
     ! #1 - at the reference coordinate of atom #1
     !
     x(1,:) = y(1,:)
     !
     ! #2 - on the vector between refernce atoms 1 and 2
     !
     r1 = sqrt( sum( ( y(2,:)-y(1,:) )**2 ) )
     !
     n1(:) = ( y(2,:)-y(1,:) )/r1
     !
     x(2,:) = x(1,:)+n1(:)*r(1)
     !
     iangle = 0
     idihedral = 0
     !
     do iatom = 3,trove%NAtoms
        !
        J = trove%zmatrix(iatom)%connect(4)
        !
        ! it is special when we work with the third atom
        ! it will be distinguished by "J"
        !
        if (iatom == 3) J = -1
        !
        rbond = r(iatom-1)
        !
        iangle = iangle + 1
        alpha13 = r(trove%Nbonds+iangle)
        !
        ! connections
        !
        p0 = iatom
        p1 = trove%zmatrix(iatom)%connect(1)
        p2 = trove%zmatrix(iatom)%connect(2)
        p3 = trove%zmatrix(iatom)%connect(3)
        !
        select case (J) 
          !
        case(-1:0)
           !
           p4 = p1
           !
        case(1)
           !
           p4 = p1
           !
        case(2)
           !
           p4 = p2
           p2 = p1
           p1 = p4
           !
        end select 
        !
        ! The rotationa back is defined by chi,theta
        ! We shift the coordinate system to have z' be coincided with the 1st reference vector n1 and to have 
        ! the 2nd ref. vector n2 be in the x'z'-plane.  
        !
        ! the 1st zmatrix ref. vector 
        !
        r1 = sqrt( sum( ( x(p2,:)-x(p1,:) )**2 ) )
        n1(:) = ( x(p2,:)-x(p1,:) )/r1
        !
        ! we keep track also of the reference vectors in similar way 
        !
        c1(:) = ( y(p2,:)-y(p1,:) )/sqrt( sum( ( y(p2,:)-y(p1,:) )**2 ) )
        !
        ! the 2d ref. vector is the second zmatrix-ref. vector, 
        ! unless it is the atom # 3, which has none. In the latter case 
        ! we use the reference y-atom #3, which we require to be in the same plane as the target atom #3
        !
        if (iatom == 3) then 
           n2(:)  = ( y(p0,:)-y(p4,:) )/sqrt( sum( ( y(p0,:)-y(p4,:) )**2 ) )
           c2(:)  = n2(:)
        else
           n2(:)  = ( x(p3,:)-x(p4,:) )/sqrt( sum( ( x(p3,:)-x(p4,:) )**2 ) )
           c2(:)  = ( y(p3,:)-y(p4,:) )/sqrt( sum( ( y(p3,:)-y(p4,:) )**2 ) )
        endif 
        !
        ! n3 vector will show the direction of the targeted atom # iatom
        !
        n3(:) = 0 
        !
        p1 = trove%zmatrix(iatom)%connect(1)
        !
        c3(:) = ( y(p0,:)-y(p1,:) )/sqrt( sum( ( y(p0,:)-y(p1,:) )**2 ) )
        !
        ! common solution
        !
        n3(3) = cos(alpha13) 
        alpha12 = acos(sum( n1(:)*n2(:) ))
        !
        ! The targeted vector depends on the J value:
        !
        select case (J) 
          !
        case(-1) 
           !
           ! this case for iatom=3 is similar to J = 1, i.e. will be defined using 
           ! the type 1 dihedral angle tau, which is zero
           !
           tau = 0
           !
           if (abs(alpha12)>small_) then
             !
             n3(2) = 0  ! tau/sin(alpha12)
             n3(1) = sqrt(cos(alpha13)**2*cos(alpha12)**2+1.0_ark-cos(alpha13)**2-cos(alpha12)**2)/sin(alpha12)
             !
           endif 
           !
        case(0)
           !
           iangle = iangle + 1 
           !
           alpha23 = r(trove%Nbonds+iangle)
           !
           tau_2 = 1.0_ark-cos(alpha13)**2-cos(alpha23)**2-cos(alpha12)**2 & 
                  +2.0_ark*cos(alpha13)*cos(alpha23)*cos(alpha12)
           !
           if ( tau_2<-sqrt(small_) ) then 
              !
              write (out,"('fromlocal2cartesian: tau**2<0: ',f18.8)") tau_2
              stop 'fromlocal2cartesian: tau**2<0'
              !
           elseif ( tau_2<0.0_ark) then 
              !
              tau = 0.0_ark
              !
           else 
              !
              tau = sqrt(tau_2)
              !
           endif 
           !
           if (abs(alpha12)>small_) then 
             !
             n3(2) = sqrt(tau_2)/sin(alpha12)
             n3(1) = ( cos(alpha23)-cos(alpha13)*cos(alpha12) )/sin(alpha12)
             !
           endif 
           !
        case(1)
           !
           idihedral = idihedral + 1 
           !
           tau = r(trove%Nbonds+trove%Nangles+idihedral)
           !
           if (abs(alpha12)>small_) then
             !
             n3(2) =-tau/sin(alpha12)
             n3(1) = sqrt(cos(alpha13)**2*cos(alpha12)**2+1.0_ark-tau**2-cos(alpha13)**2-cos(alpha12)**2)/sin(alpha12)
             n3(1) =-sign(n3(1),tau)
             !
           endif 
           !
        case (2)
           !
           idihedral = idihedral + 1 
           !
           delta = r(trove%Nbonds+trove%Nangles+idihedral)
           !
           if (abs(sin(alpha12))>small_) then 
             !
             !n3(1) = cos(delta)*sin(alpha13) 
             !n3(2) = sin(delta)**2+cos(delta)**2*cos(alpha13)**2-cos(alpha13)**2
             !n3(2) =-sqrt( sin(delta)**2+cos(delta)**2*cos(alpha13)**2-cos(alpha13)**2 )
             !
             n3(1) =-sin(alpha13)*cos(delta)
             n3(3) =-cos(alpha13) 
             n3(2) = sin(alpha13)*sin(delta)
             !
           endif 
           !
           !
        end select 
        !
        ! "Back rotation" #1 angles theta and chi 
        !
        theta0 = acos(n1(3))
        !
        chi0 = atan2(-n1(2),n1(1)) + pi 
        !
        !sign_t = sign(1.0_ark,n1(1)*sin(theta0)*cos(chi0))
        !if (sign_t>0) 
        !chi0 = chi0 + pi
        !
        !c_t=(/-sin(theta0)*cos(chi0),sin(theta0)*sin(chi0),cos(theta0)/)
        !
        rot1 = FL_euler_rotait(theta0,0.0_ark,chi0)
        !
        ! rotation of the vector n2 is needed to define the next back rotation,
        ! while rotation of n1 is only for test purposes 
        ! 
        n1(:) = matmul(transpose(rot1),n1(:))
        n2(:) = matmul(transpose(rot1),n2(:))
        !
        ! Second back rotation
        !
        phi0 = atan2(-n2(2),n2(1))
        !
        rot2 = FL_euler_rotait(0.0_ark,phi0,0.0_ark)
        !
        ! only for test porposes 
        !
        if (verbose>=4) n1(:) = matmul(transpose(rot2),n1(:))
        if (verbose>=4) n2(:) = matmul(transpose(rot2),n2(:))
        !
        ! here we found the back transformation
        !
        rot1 = matmul(rot1,rot2)
        !
        ! adjust vector n3 if we see some problem 
        !
        select case (J) 
          !
        case(-1) 
           !
        case(0)
           !
        case(1)
           !
           ! we constract the dihedral angle (akin direction) type 1 from n1,n2, and n3 and 
           ! compare its direction with a similar dihedral angle  constructed from the reference vectors 
           ! c1,c2, and c3. 
           ! if it turns out that the directions are different we take the signs from the transformed 
           ! vector c3 and assign them to the vector n3
           !
           tau = 0
           !
           do ix =1,3
             do jx =1,3
               do kx =1,3
                  tau = tau + epsil(ix,jx,kx)*n1(ix)*n2(jx)*n3(kx)
               enddo
             enddo
           enddo
           !
           sign_t = 0
           !
           do ix =1,3
             do jx =1,3
               do kx =1,3
                  sign_t = sign_t + epsil(ix,jx,kx)*c1(ix)*c2(jx)*c3(kx)
               enddo
             enddo
           enddo
           !
           ! in the transformation we keep the atom ordering, so to say, the left
           ! hand orientation has to be preserved, so the righ hand one
           !
           if (sign_t*tau<-small_) then
             !
             !c1(:) = matmul(transpose(rot1),c1(:))
             !c2(:) = matmul(transpose(rot1),c2(:))
             !
             c3(:) = matmul(transpose(rot1),c3(:))
             !
             n3(:) = sign(n3(:),c3(:))
             !
             !
           endif 
           !
           !
        case (2)
           !
           !
        end select 
        !
        !
        ! Applying the rotation 
        !
        n3(:) = matmul(rot1,n3(:))
        !
        ! only for test porposes 
        !
        if (verbose>=4) n1(:) = matmul(rot1,n1(:))
        if (verbose>=4) n2(:) = matmul(rot1,n2(:))
        !
        p1 = trove%zmatrix(iatom)%connect(1)
        !
        x(p0,:) = x(p1,:) + n3(:)*rbond
        !
        ! Check that the found vector x(p0,:) is not too far from the original vector y(p0,:)
        !
        n2(:) = ( y(p0,:)-y(p1,:) )/sqrt( sum( ( y(p0,:)-y(p1,:) )**2 ) )
        !
        theta2 = sum( n2(:)*n3(:) )
        !
        if ( abs(theta2)>1.0_ark+sqrt(small_) ) then 
           write (out,"('fromlocal2cartesian: cos(theta2)>1: ',f18.8)") theta2
           stop 'fromlocal2cartesian: cos(theta2)>1'
        elseif ( abs(theta2)>1.0_ark) then 
           !
           theta2 = acos(sign(1.0_ark,theta2))
           !
        else 
           !
           theta2 = acos( theta2 )
           !
        endif 
        !
        !theta2 = acos( sum( n2(:)*n3(:) ) )
        !
        if ( abs(theta2)>0.5_rk*pi ) then 
           write (out,"('fromlocal2cartesian: the vector x(',i4,') is too far from y, theta = ',f18.8)") p0,theta2
           stop 'fromlocal2cartesian: x too far from y'
        endif 
        !
        continue
        !
        !call FLfromcartesian2local(x,local)
        !
        !continue
        !
     enddo
     !
     ! Find center of mass and correct the vectors
     !
     do ix = 1,3 
       !
       CM_shift = sum( x(:,ix)*trove%mass(:) )/sum( trove%mass(:) )
       x(:,ix) = x(:,ix) - CM_shift
       !
     enddo 
     !
     ! Check if we found the silution
     !
     call FLfromcartesian2local(x,local)
     !
     !local = local - x
     !
     do n0 = 1,trove%Nmodes
         !
         if ( abs( abs( (local(n0)-r(n0) ) ) - 2.0_ark*pi )<sqrt(small_) ) then
            !
            local(n0) = r(n0)
            !
         endif
         !
     enddo
     !
     dlocal(:) = local(:)-r(:)
     !
     if ( any( abs(dlocal(:) )>100000.0*sqrt(small_) ) ) then 
          write(out,"(t3,'fromlocal2cartesian: found local coords are wrong: ',30d18.8)") r
          write(out,"(t3,'                                     compare with: ',30d18.8)") local 
          stop 'fromlocal2cartesian: found local coords are wrong'
     endif
     !
     ! Second Eckart equation
     ! 
     select case (trim(axis_system)) 
     !
     case default
     !
       write (out,"('fromlocal2cartesian: axis system ',a,' unknown')") trim(axis_system)
       stop 'fromlocal2cartesian - bad axis system'
       !
     case('Eckart') 
       !
       !call second_eckart_by_fit(x)
       !
       !y(:,:) =  trove%b0(:,:,irho)
       !
       ! estimate the rms of deviations of each atom from their reference positions
       !
       !rms = sqrt( sum( ( x(:,:)-y(:,:) )**2 )/12.0_rk )
       !
       !x(:,1:2) = x(:,2:1:-1)
       !
       !rms = sqrt( sum( ( x(:,:)-y(:,:) )**2 )/12.0_rk )
       !
       !continue
       !
       !case('Full_Eckart')
       !
       !cartesian = x
       !
       ! Second Eckart equation
       ! 
       !do ix = 1,3 
       !   do jx = 1,3 
       !      a(ix,jx) =  sum(trove%mass(:)*x(:,ix)*x(:,jx) )
       !   enddo
       !   !
       !enddo
       !
       !a = real(a,kind=rk)
       !
       !call lapack_syev(a,b(:,1))
       !
       ! Found coordinate transformation "c"
       !
       !c = real(a,kind=rk)
       !
       ! Transformation of a0 
       !
       !do ix = 1,4
       !   x(ix,:) = matmul(transpose(c),cartesian(ix,:))
       !enddo
       !
       ! estimate the rms of deviations of each atom from their reference positions
       !
       !rms = sqrt( sum( ( x(:,:)-y(:,:) )**2 )/12.0_rk )
       !
       !
     !case('PAS') 
       !
       !do kx = 1,3 
       !   !
       !   f_t = epsil(ix,jx,kx)**2
       !   !
       !   e_t = e_t + f_t*cartesian(n0,jx)*cartesian(n0,kx)*trove%mass(n0)
       !   !
       !enddo
       !
     end select 
     !
     ! for test purposes only 
     !
     if (verbose>=5) then
        !
        ! First Eckart equation
        ! 
        do ix = 1,3 
           !
           a_t = sum(trove%mass(:)*x(:,ix))
           !
           if (abs(a_t)>10.0_rk**(-rk)) then 
                write(out,"('fromlocal2cartesian: a0 is not a solution of Eckart 1 for  ix =',i4,d18.8)") ix,a_t
                   stop 'fromlocal2cartesian: a0 is not solution of Eckart1'
           endif
           !
        enddo
        !
        ! Second Eckart equation
        !
        do ix = 1,3 
           !
           select case (trim(axis_system)) 
             !
           case('Eckart') 
             !
             a_t = 0
             do jx = 1,3 
                do kx = 1,3 
                   a_t = a_t + epsil(ix,jx,kx)*sum(trove%mass(:)*y(:,jx)*x(:,kx) )
                enddo
             enddo
             !
             if (abs(a_t)>2000.0_rk*sqrt(small_)) then 
                 if (verbose>=5) write(out,"('fromlocal2cartes: x is not a solut of Eckart 2 for ix,jx =',2i4,d18.8)")&
                                       ix,jx,a_t
                !stop 'fromlocal2cartesian: x is not solution of Eckart2'
             endif
             !
           case('Full_Eckart')
             !
             do jx = 1,3 
                 if (ix/=jx) then  
                    !
                    a_t =  sum(trove%mass(:)*x(:,ix)*x(:,jx) )
                    !
                    if (abs(a_t)>100.0_rk*small_) then 
                        write(out,"('fromlocal2cartes: x is not a solution of Eckart 2 for ix,jx =',2i4,d18.8)") & 
                                     ix,jx,a_t
                        stop 'fromlocal2cartesian: x is not solution of Eckart2'
                    endif
                    !
                 endif
                 !
             enddo
             !
           end select 
           !
        enddo
        !
     endif
     !
     !
     !call FLfromcartesian2local(x,local)
     !
     cartesian = x
     !
     !call FLfromcartesian2local(cartesian,local)
     !
     !
     if (verbose>=6) write(out,"('fromlocal2cartesian/end')") 
     !
   end subroutine fromlocal2cartesian



   subroutine from_local2cartesian_by_fit(chi,chi_eq,xi,xi0)


    real(ark),intent(in)  :: chi(1:trove%Nmodes) 
    real(ark),intent(in)  :: chi_eq(1:trove%Nmodes) 
    real(ark),intent(out) :: xi(1:trove%Nmodes+6)
    ! 
    ! zero approximation
    !
    real(ark),optional,intent(in) :: xi0(1:trove%Nmodes+6)
    !
    real(ark)            :: cartesian(trove%Natoms,3)
    !
    real(ark)   :: a0(trove%Natoms,3)
    real(ark) :: eps(1:trove%Nmodes+6),parold(1:trove%Nmodes+6)
    real(ark) :: Xright(1:trove%Nmodes+6),Xleft(1:trove%Nmodes+6)
    integer(ik) :: ivar(1:trove%Nmodes+6)
    real(ark) :: rjacob(1:trove%Nmodes+6,1:trove%Nmodes+6)

    real(rk) :: am(1:trove%Nmodes+6,1:trove%Nmodes+6),bm(1:trove%Nmodes+6,1)

    real(ark) ::stadev_old,stability,stadev,ssq,stadev_best,tempx,deltax
    real(ark) :: Bmat(trove%Nmodes,trove%Natoms,3)
    !
    !real(ark) :: xi_eq(1:trove%Nmodes+6),chi_t(1:trove%Nmodes) 
    !
    integer(ik) :: iter,numpar,n0,itmax,i,parmax,irow,icolumn,ieq,x0,numvar,Natoms,x1,jmode,Nmodes
    integer(ik) :: irho,jacob_type
    !
    Natoms = trove%Natoms
    Nmodes = trove%Nmodes
    a0 = trove%a0
    !
    rjacob = 0 
    iter = 0
    stadev_old = 1.e10
    stability =  1.e10
    stadev    =  1.e10
    jacob_type = 1

    stadev_best = sqrt(small_)! 1e-14
    itmax = fititermax
    parmax = trove%Nmodes+6
    !factordeltax = 1e-5
    !
    deltax=sqrt(trove%fdstep(1))
    !
    ivar = 1
    !
    ! The address of the coordinate "rho" in the 1d rho-table 
    !
    !irho =     nint( ( trove%chi_ref(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik )
    !
    irho = FLirho
    !
    ! Initial values for xi = local
    !
    if (present(xi0)) then 
       !
       xi = xi0
       !
    else
       !
       numvar = 0 
       do n0 = 1,Natoms
          do x0 = 1,3
             numvar= numvar+1
             !xi(numvar) = trove%a0(n0,x0) + sum( trove%Amat(n0,x0,:)*(chi(:)-chi_eq(:)) )
             !
             xi(numvar) = trove%b0(n0,x0,irho) + &
                          sum( trove%Amatrho(n0,x0,1:trove%Nmodes_e,irho)*&
                             (chi(1:trove%Nmodes_e)-chi_eq(1:trove%Nmodes_e)) )
             !
             !xi_eq(numvar) = trove%b0(n0,x0,irho)
             !
          enddo
       enddo
       !
    endif
    !
    !call FLfromcartesian2local(trove%b0(:,:,irho),chi_t)
    !
    parold = xi
    !
    outer_loop: & 
    do while( iter<fititermax .and. stadev>stadev_best )   
       !
       numpar = sum(ivar)
       iter = iter + 1
       ssq=0
       !
       call calc_all_equations(xi,eps)
       !
       ! st. square deviation 
       !
       ssq=sum(eps(:)**2)
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       if (jacob_type==1) then 
          !
          ieq = 0
          do n0 = 1,trove%Natoms
            do x0 = 1,3
               !
               ieq = ieq +1 
               cartesian(n0,x0)=xi(ieq)
               !
            enddo
          enddo
          !
          call Bmat_generation(trove%b0(:,:,irho),Bmat)
          !
          numvar= 0
          do n0 = 1,Natoms
             !
             do x0 = 1,3
                !
                numvar= numvar+1
                !
                ! a) translational part (3 members: x,y,z)
                !
                rjacob(x0,numvar) = trove%mass(n0)
                !
                ieq = 3
                !
                ! b) rotational part (3 members: x,y,z)       
                !
                do x1 = 1,3
                   ieq = ieq +1 
                   rjacob(ieq,numvar) = sum( epsil(x0,x1,:)*trove%b0(n0,:,irho) )*trove%mass(n0)
                enddo
                !
                ! c) Vibrational part (Nmodes members)
                !
                do jmode = 1,Nmodes
                  !
                  ieq = ieq +1 
                  rjacob(ieq,numvar) = trove%Bmatrho(jmode,n0,x0,irho)
                  !
                enddo
                !
             enddo
          enddo
       else
          !
          do  i=1,parmax
              !
              tempx=xi(i)
              !
              xi(i)=tempx+deltax
              call calc_all_equations(xi,Xright)
              !
              xi(i)=tempx-deltax
              call calc_all_equations(xi,Xleft)
              rjacob(:,i)=(Xright(:)-Xleft(:))/(2.d0*deltax)
              xi(i)=tempx
              !
          enddo
          !
       endif
       !
       if (fititermax>=0) then
         !
         ! We constract a set of linear equations A x = B
         !
         ! form A matrix 
         !
         do irow=1,numpar       !==== row-...... ====!
           do icolumn=1,irow    !==== column-....====!
             am(irow,icolumn)=sum(rjacob(:,icolumn)*rjacob(:,irow))
             am(icolumn,irow)=am(irow,icolumn)
           enddo
         enddo
         !
         ! form B matrix 
         !
         do irow=1,numpar       !==== row-...... ====!
           bm(irow,1)=sum(eps(:)*rjacob(:,irow))
         enddo   
         !
         ! Solve the set of linear equations 
         !
         ! lapack_gelss - solves a linera equation by least squares method 
         ! 
         call lapack_gelss(am(:,:),bm(:,:))
         !
         xi(:)=xi(:)-bm(:,1)
         !
         stadev=sqrt(ssq/real(numpar,kind=rk))
         !
         if (stadev>1.e3.and.jacob_type==1) then 
            !
            jacob_type = 2
            xi = parold
            iter = 0 
            !
         endif
         !
      else   ! Itmax = 0, so there is no fit
             ! only straightforward calculations 
             !
         stadev_old=stadev
         stadev=sqrt(ssq/real(numpar,kind=rk))
         !
      endif 

      !
    enddo  outer_loop ! --- iter
   ! ##########################
   !
   if (iter==fititermax) then 
      write (out,"('from_local2cartesian_by_fit: no convergence after ',i8,' iterations')") iter
      write (out,"('stadev = ',f18.8)") stadev
      write (out,"('eps = ')")
      write (out,"(40f18.8)") eps(1:min(40,size(eps)))
      stop 'from_local2cartesian_by_fit: - bad axis system'
   endif 
   !
   ! Calculated Cartesian coordinates 
   ieq = 0
   do n0 = 1,trove%Natoms
     do x0 = 1,3
        !
        ieq = ieq +1 
        cartesian(n0,x0)=xi(ieq)
        !
     enddo
   enddo
   !
   ! We store the fitted cartesian coordinates xi here for the latter use
   !
   !cartesian_stored = xi
   !
   contains 

     !
     ! Caclulate the equations 
     !
     subroutine calc_all_equations(xi,equations)
      !
      real(ark),intent(in) :: xi(3*trove%Natoms)
      real(ark),intent(out) :: equations(3*trove%Natoms)
      !
      real(ark) :: cartesian(trove%Natoms,3)
      integer(ik) :: i,ix,n
      real(ark)   :: Eckart(6),r(trove%Ncoords),chi_t(trove%Nmodes)
      logical     :: dir

       !
       i = 0 
       !
       do n = 1,trove%Natoms
         !
         do ix =1,3
           !
           i = i +1 
           cartesian(n,ix) = xi(i)
           !
         enddo
         !
       enddo
       !
       call calc_eckart(cartesian,Eckart)
       !
       equations(1:6) = Eckart(1:6)
       !
       call FLfromcartesian2local(cartesian,r)
       !
       dir = .true.
       !
       chi_t = MLcoordinate_transform_func(r,size(chi_t),dir)
       !
       equations(7:trove%Nmodes+6) = chi_t(1:trove%Nmodes)-chi(1:trove%Nmodes)
       !
       if ( any( abs( mod(chi_t(:),pi) )<small_ ) ) then
          !
          do i = 1,trove%Nmodes
             !
             if ( abs( mod(chi_t(i),pi) )<small_ ) then
                !
                equations(6+i) = mod(equations(6+i)+2.0_ark*pi,2.0_ark*pi)
                !
             endif
             !
          enddo
          !
       endif
       !
     end subroutine calc_all_equations


   end subroutine from_local2cartesian_by_fit





   subroutine FLfromcartesian2local(cartesian,r)

    real(ark),intent(in)   :: cartesian(trove%Natoms,3)
    real(ark),intent(out)  :: r(trove%Ncoords)

    real(ark)   :: rcon(trove%Natoms,trove%Natoms),tau_sign
    real(ark)   :: tau,cosa1,cosa2,cosa3,sindelta,norm_2,cosa,&
                   a_t(3),a_t1(3),a_t2(3),a_t3(3),delta,B,vec1,vec2,dvec1(3),dvec2(3),r1,r2,r3,cosu,cosv,sinu,sinv,&
                   u(3),v(3),w(3),cosdelta,phi,rmat(3,3),fmod,u0(3),v0(3),sina1,sina2

    integer(ik) ::  ibond,iangle,kappa,zeta,k1,k2
    integer(ik) ::  n1,n2,n3,n4,n0,ix,iy,iz,J
     !
     ! geometrically defined coordinates 
     !
     ! 1. stretching coordinates (bonds)
     !
     rcon = 0
     do ibond = 1,trove%Nbonds
        !
        n1 = trove%bonds(ibond,1)
        n2 = trove%bonds(ibond,2)
        a_t(:) = cartesian(n1,:) - cartesian(n2,:)
        !r(ibond) = sqrt( sum( a_t(:)**2 ) )
        !
        rcon(n1,n2) = sqrt( sum( a_t(:)**2 )  )
        rcon(n2,n1) = rcon(n1,n2)
        r(ibond)    = rcon(n1,n2) ! -trove%req(ibond)
        !
     enddo
     !
     ! 2. bending coordinates (angles)
     !
     do iangle = 1,trove%Nangles
        !
        n1 = trove%angles(iangle,1)
        n0 = trove%angles(iangle,2)
        n2 = trove%angles(iangle,3)
        !
        a_t1(:) = cartesian(n1,:) - cartesian(n0,:)
        a_t2(:) = cartesian(n2,:) - cartesian(n0,:)
        !
        cosa = sum(a_t1(:)*a_t2(:) )/( rcon(n1,n0)*rcon(n2,n0) )
        !
        r(trove%Nbonds+iangle) = acos(cosa)
        !
     enddo
     !
     do iangle = 1,trove%Ndihedrals
       !
       J = trove%dihedtype(iangle)
       !
       select case (J)
       !
       case(1) ! type 1 
          !
          n1 = trove%dihedrals(iangle,1)
          n4 = trove%dihedrals(iangle,2)
          n3 = trove%dihedrals(iangle,3)
          n2 = trove%dihedrals(iangle,4)
          !
          a_t1(:) = cartesian(n1,:) - cartesian(n4,:)
          a_t2(:) = cartesian(n2,:) - cartesian(n4,:)
          a_t3(:) = cartesian(n3,:) - cartesian(n4,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          r3 =  sqrt(sum(a_t3(:)**2))
          !
          delta = 0 
          !
          do ix =1,3
            do iy =1,3
              do iz =1,3
                 delta = delta + epsil(ix,iy,iz)*a_t1(ix)*a_t2(iy)*a_t3(iz)
              enddo
            enddo
          enddo
          !
          tau = delta/( r1*r2*r3 )
          !
          cosa1 = sum(a_t2(:)*a_t3(:) )/( r2*r3 )
          cosa2 = sum(a_t1(:)*a_t2(:) )/( r1*r2 )
          cosa3 = sum(a_t1(:)*a_t3(:) )/( r1*r3 )
          !
          norm_2 = 3.0_ark-cosa3**2-cosa2**2-cosa1**2+2.0_ark*cosa3*cosa1-2.0_ark*cosa2+&
                   2.0_ark*cosa2*cosa3-2.0_ark*cosa1+2.0_ark*cosa2*cosa1-2.0_ark*cosa3
          !
          if (norm_2<small_.and.abs(tau)<small_) then 
            !
            sindelta = 1.0_ark
            !
          elseif(norm_2<small_) then 
            !
            write (out,"('FLfromcartesian2local: norm2 = ',f18.8', delta = infty!')") norm_2
            write (out,"('Consider change difftype ')")
            stop 'FLfromcartesian2local - bad norm2'
            !
          else
            !
            sindelta = tau/sqrt(norm_2)
            !
          endif
          !
          if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('FLfromcartesian2local: sindelta>1: ',f18.8)") sindelta
             write (out,"('Consider change difftype ')")
             stop 'FLfromcartesian2local - bad sindelta'
             !
          elseif ( sindelta>=1.0_ark) then 
             !
             r(trove%Nbonds+trove%Nangles+iangle) = 0.0_ark
             !
          else 
             r(trove%Nbonds+trove%Nangles+iangle) = asin(sindelta)
             !
          endif 
          !
          !r(trove%Nbonds+trove%Nangles+iangle) = delta
          !
       case(-2,2,-202,202,-302,302) ! type 2   B = (a*b)/(|a|*|b|), a = [y1 times y2]; b = [y2 times y3]
          !
          n1 = trove%dihedrals(iangle,1)
          n2 = trove%dihedrals(iangle,2)
          n3 = trove%dihedrals(iangle,3)
          n4 = trove%dihedrals(iangle,4)
          !
          u(:) = cartesian(n4,:) - cartesian(n3,:)
          w(:) = cartesian(n2,:) - cartesian(n3,:)
          v(:) = cartesian(n1,:) - cartesian(n2,:)
          !
          if (J<0) w = -w
          !
          !endif 
          !
          r1 =  sqrt(sum(u(:)**2))
          r2 =  sqrt(sum(w(:)**2))
          r3 =  sqrt(sum(v(:)**2))
          !
          u = u/r1
          w = w/r2
          v = v/r3
          !
          dvec1(:) = 0 
          dvec2(:) = 0 
          !
          do iy =1,3
            do iz =1,3
               dvec1(:) = dvec1(:) + epsil(:,iy,iz)*u(iy)*w(iz)
               dvec2(:) = dvec2(:) + epsil(:,iy,iz)*v(iy)*w(iz)
            enddo
          enddo
          !
          cosu = sum( u(:)*w(:) )
          cosv =-sum( v(:)*w(:) )
          !
          sinu = sqrt(1.0_ark-cosu**2)
          sinv = sqrt(1.0_ark-cosv**2)
          !
          B = sum( dvec1(:)*dvec2(:) )/( sinu*sinv )
          !
          tau_sign = 0
          !
          do iy =1,3
            do iz =1,3
               tau_sign = tau_sign - sum(epsil(:,iy,iz)*w(:)*dvec1(iy)*dvec2(iz))
            enddo
          enddo
          !
          dvec1(:) = MLvector_product(u(:),w(:))
          dvec2(:) = MLvector_product(v(:),w(:))
          !
          vec1 = sqrt( sum(dvec1(:)**2) )
          vec2 = sqrt( sum(dvec2(:)**2) )
          !
          dvec1 = dvec1/vec1
          dvec2 = dvec2/vec2
          !
          B = sum( dvec1(:)*dvec2(:) )
          !
          a_t = MLvector_product(dvec1(:),dvec2(:))
          !
          tau_sign = -sum( w(:)*a_t(:) )
          !
          cosdelta = B
          !
          a_t = MLvector_product(dvec1(:),dvec2(:))
          !
          sindelta = sqrt(sum(a_t(:)**2))
          !
          delta = atan2(sindelta,cosdelta)
          !
          !
          !if (abs(tau_sign)<(small_)) tau_sign = small_a
          !
          !if ( abs(B)>1.0_ark+sqrt(small_) ) then 
          !   !
          !   write (out,"('FLfromcartesian2local: costau>1: ',f18.8)") B
          !   stop 'FLfromcartesian2local - bad costau'
          !   !
          !elseif ( B>=1.0_ark) then 
          !   delta = 0.0_ark
          !elseif ( B<=-1.0_ark) then 
          !   delta = pi
          !else 
          !   delta = acos(B)
          !endif
          !
          fmod = 2.0_ark*pi
          !
          ! s[ecial case of (EM)-symmetry with tau=0..720 deg
          !if (trove%periodic.and.abs(job%bset(trove%Nmodes)%borders(2)-4.0_ark*pi)<sqrt(small_)) fmod = 4.0_ark*pi
          !
          if (tau_sign<-small_a) then 
             !
             delta = fmod-delta
             !
          endif
          !
          if ( delta<-small_.or.delta>fmod+small_ ) then 
            !
            delta  = mod(delta+fmod,fmod)
            !
          endif
          !
          r(trove%Nbonds+trove%Nangles+iangle) = delta
          !
       case(-402,402) ! type 2   B = (a*b)/(|a|*|b|), a = [y1 times y2]; b = [y2 times y3]
          !           ! special case of tau defined for the range 0..720
          !
          n1 = trove%dihedrals(iangle,1)
          n2 = trove%dihedrals(iangle,2)
          n3 = trove%dihedrals(iangle,3)
          n4 = trove%dihedrals(iangle,4)
          !
          u(:) = cartesian(n4,:) - cartesian(n3,:)
          w(:) = cartesian(n2,:) - cartesian(n3,:)
          v(:) = cartesian(n1,:) - cartesian(n2,:)
          !
          if (J<0) w = -w
          !
          !endif 
          !
          r1 =  sqrt(sum(u(:)**2))
          r2 =  sqrt(sum(w(:)**2))
          r3 =  sqrt(sum(v(:)**2))
          !
          u = u/r1
          w = w/r2
          v = v/r3
          !
          dvec1(:) = 0 
          dvec2(:) = 0 
          !
          do iy =1,3
            do iz =1,3
               dvec1(:) = dvec1(:) + epsil(:,iy,iz)*u(iy)*w(iz)
               dvec2(:) = dvec2(:) + epsil(:,iy,iz)*v(iy)*w(iz)
            enddo
          enddo
          !
          cosu = sum( u(:)*w(:) )
          cosv =-sum( v(:)*w(:) )
          !
          sinu = sqrt(1.0_ark-cosu**2)
          sinv = sqrt(1.0_ark-cosv**2)
          !
          B = sum( dvec1(:)*dvec2(:) )/( sinu*sinv )
          !
          tau_sign = 0
          !
          do iy =1,3
            do iz =1,3
               tau_sign = tau_sign - sum(epsil(:,iy,iz)*w(:)*dvec1(iy)*dvec2(iz))
            enddo
          enddo
          !
          dvec1(:) = MLvector_product(u(:),w(:))
          dvec2(:) = MLvector_product(v(:),w(:))
          !
          vec1 = sqrt( sum(dvec1(:)**2) )
          vec2 = sqrt( sum(dvec2(:)**2) )
          !
          dvec1 = dvec1/vec1
          dvec2 = dvec2/vec2
          !
          B = sum( dvec1(:)*dvec2(:) )
          !
          a_t = MLvector_product(dvec1(:),dvec2(:))
          !
          tau_sign = -sum( w(:)*a_t(:) )
          !
          cosdelta = B
          !
          a_t = MLvector_product(dvec1(:),dvec2(:))
          !
          sindelta = sqrt(sum(a_t(:)**2))
          !
          delta = atan2(sindelta,cosdelta)
          !
          u0(:) = trove%b0(n4,:,0) - trove%b0(n3,:,0) 
          v0(:) = trove%b0(n1,:,0) - trove%b0(n2,:,0) 
          !
          u0 = u0/sqrt(sum(u0(:)**2))
          v0 = v0/sqrt(sum(v0(:)**2))
          !
          cosa1 = sum(u0*u)
          cosa2 = sum(v0*v)
          !
          a_t = MLvector_product(u0(:),u(:))
          sina1 = sqrt(sum(a_t(:)**2))
          !          
          a_t = MLvector_product(v0(:),v(:))
          sina2 = sqrt(sum(a_t(:)**2))
          !
          fmod = 2.0_ark*pi
          if ( abs(dvec1(1))<sqrt(small_) ) then 
            if ( dvec1(2)>0.0_ark ) then 
                fmod = 4.0_ark*pi
            endif
          elseif (dvec1(1)<0.0_ark) then 
            fmod = 4.0_ark*pi
          endif
          !
          if (tau_sign<-sqrt(small_a)) then 
             !
             delta = fmod-delta
             !
          endif
          !
          ! special case of (EM)-symmetry with tau=0..720 deg
          if ( fmod > 2.0_ark*pi.and.tau_sign> small_a ) then 
            delta = delta + 2.0_ark*pi
          endif
          !
          if ( delta<-small_.or.delta>fmod+small_ ) then 
            !
            delta  = mod(delta+fmod,fmod)
            !
          endif
          !
          r(trove%Nbonds+trove%Nangles+iangle) = delta          
          !
       case(101) ! The special bond-angles for the linear molecule case 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = cartesian(n1,:) - cartesian(n0,:)
          a_t2(:) = cartesian(n2,:) - cartesian(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          !
          u =  a_t1(:)/r1
          v =  a_t2(:)/r2
          !
          w(:) = MLvector_product(u,v)
          !
          ! special angle is -arcsin( kappa . [uxv] ), 
          ! i.e. the scalar product kappa.w is a kappa component of w
          !
          sindelta = w(kappa)
          !
          if ( abs(sindelta)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('FLfromcartesian2local: sindelta>1: ',f18.8)") sindelta
             write (out,"('Consider change difftype ')")
             stop 'FLfromcartesian2local - bad sindelta'
             !
          elseif ( abs(sindelta)>=1.0_ark) then 
             !
             phi = 0.0_ark
             !
          else
             ! 
             phi = -asin(sindelta)
             !
          endif
          !
          r(trove%Nbonds+trove%Nangles+iangle) = phi
          !
          !if (kappa==2) then
          !   r(trove%Nbonds+trove%Nangles+iangle) = -r(trove%Nbonds+trove%Nangles+iangle)
          !endif 
          !
       case(102) ! The special bond-angles for the linear molecule case 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = cartesian(n1,:) - cartesian(n0,:)
          a_t2(:) = cartesian(n2,:) - cartesian(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          !
          u =  a_t1(:)/r1
          v =  a_t2(:)/r2
          !
          w(:) = MLvector_product(u,v)
          !
          ! special angle is -arcsin( kappa . [uxv] ), 
          ! i.e. the scalar product kappa.w is a kappa component of w
          !
          sindelta = w(kappa)
          !
          r(trove%Nbonds+trove%Nangles+iangle) = sindelta
          !
       case(103,105) ! The special bond-angles for the linear molecule case 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = cartesian(n1,:) - cartesian(n0,:)
          a_t2(:) = cartesian(n2,:) - cartesian(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          !
          u =  a_t1(:)/r1
          v =  a_t2(:)/r2
          !
          ! xyz
          !
          r(trove%Nbonds+trove%Nangles+iangle) = a_t1(kappa)
          !
          if (j==105) then
            !
            r(trove%Nbonds+trove%Nangles+iangle) = u(kappa)
            !
          endif
          !
          zeta = trove%zmatrix(n1)%connect(3)
          !
          do ibond = 1,trove%Nbonds
             !
             k1 = trove%bonds(ibond,1)
             k2 = trove%bonds(ibond,2)
             !
             if (k1/=n1.and.k2/=n1) cycle
             !
             r(ibond) = a_t1(zeta)
             !
             !r(ibond) = sum(a_t1(:)*a_t2(:))
             !
             !if (trove%b0(n1,zeta,0)<0) r(ibond) = -a_t1(zeta)
             !
          enddo
          !
       case(104) ! The special bond-angles for the linear molecule case 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = cartesian(n1,:) - cartesian(n0,:)
          a_t2(:) = cartesian(n2,:) - cartesian(n0,:)
          !
          ! xyz
          !
          r(trove%Nbonds+trove%Nangles+iangle) = a_t1(kappa)
          !
       case(106) ! The special bond-angles for the linear molecule case 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = cartesian(n1,:) - cartesian(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          !
          ! xyz
          !
          r(trove%Nbonds+trove%Nangles+iangle) = a_t1(kappa)/r1
          !
       case(107) ! The special bond-angles for the linear molecule case 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = cartesian(n1,:) - cartesian(n0,:)
          a_t2(:) = cartesian(n2,:) - cartesian(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          !
          !u =  (/0.0_ark,1.0_ark,0.0_ark/)
          !
          v =  a_t2(:)/r2
          !
          u(1) = 0 
          u(3) = v(2)/sqrt(v(2)**2+v(3)**2)
          u(2) =-v(3)/sqrt(v(2)**2+v(3)**2)
          !
          r2 =  sqrt(sum(u(:)**2))
          !
          if (v(3)<0) v = -v 
          !
          w(:) = MLvector_product(u,v)
          !
          rmat(1,:) = w(:)
          rmat(2,:) = u(:)
          rmat(3,:) = v(:)
          !
          a_t1(:) =matmul(rmat,a_t1(:))
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          !
          ! xyz
          !
          r(trove%Nbonds+trove%Nangles+iangle) = a_t1(kappa)
          !
       case(108) ! The special bond-angles for the linear molecule case 
          !
          ! 1 and 2 are the two outer atoms and n0 is the central atom
          !
          n1 = trove%dihedrals(iangle,1)
          n0 = trove%dihedrals(iangle,2)
          n2 = trove%dihedrals(iangle,3)
          !
          ! Cartesian component we build the two special angles for
          !
          kappa = trove%dihedrals(iangle,4)
          !
          a_t1(:) = cartesian(n1,:) - cartesian(n0,:)
          a_t2(:) = cartesian(n2,:) - cartesian(n0,:)
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          r2 =  sqrt(sum(a_t2(:)**2))
          !
          v =  a_t2(:)/r2
          !
          u(1) = 0 
          u(3) = v(2)/sqrt(v(2)**2+v(3)**2)
          u(2) =-v(3)/sqrt(v(2)**2+v(3)**2)
          !
          r2 =  sqrt(sum(u(:)**2))
          !
          if (v(3)<0) v = -v 
          !
          w(:) = MLvector_product(u,v)
          !
          rmat(1,:) = w(:)
          rmat(2,:) = u(:)
          rmat(3,:) = v(:)
          !
          a_t1(:) =matmul(rmat,a_t1(:))
          !
          r1 =  sqrt(sum(a_t1(:)**2))
          !
          ! xyz
          !
          r(trove%Nbonds+trove%Nangles+iangle) = a_t1(kappa)       
          !
          zeta = trove%zmatrix(n1)%connect(3)
          !
          do ibond = 1,trove%Nbonds
             !
             k1 = trove%bonds(ibond,1)
             k2 = trove%bonds(ibond,2)
             !
             if (k1/=n1.and.k2/=n1) cycle
             !
             r(ibond) = a_t1(zeta)
             !
          enddo          
          ! 
       end select 
       !
     enddo 
     !
   end subroutine FLfromcartesian2local


   !
   ! 3-dim rotation by three Euler angles
   !
   function FL_euler_rotait(theta,phi,chi) result (f)

      real(ark),intent(in)   :: theta,phi,chi
      real(ark),dimension(3,3) :: f

        f(1,:) = &
         (/cos(theta)*cos(phi)*cos(chi)-sin(phi)*sin(chi),  &
           cos(theta)*sin(phi)*cos(chi)+cos(phi)*sin(chi),  &
          -sin(theta)*cos(chi) /)
        f(2,:) = &
         (/-cos(theta)*cos(phi)*sin(chi)-sin(phi)*cos(chi),  &
           -cos(theta)*sin(phi)*sin(chi)+cos(phi)*cos(chi),  &
            sin(theta)*sin(chi)/)
        f(3,:) = &
         (/ sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta) /)

   end function FL_euler_rotait
   !


   !
   !
   subroutine calc_eckart(cartesian,Eckart)

    real(ark),intent(in)   :: cartesian(trove%Natoms,3)
    real(ark),intent(out)  :: Eckart(6)
    real(ark)              :: e_t,f_t
    !
    integer(ik) ::  ieq,ix,jx,kx,n0,irho
      !
      !irho = mod(nint( ( molec%chi_eq(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
      !irho =     nint( ( molec%chi_eq(trove%Nmodes)-trove%rho_border(1) )/(trove%rhostep),kind=ik )
      irho = FLirho
      !
      ! First Eckart equation
      ! 
      do ix = 1,3 
          !
          Eckart(ix) = sum(trove%mass(:)*cartesian(:,ix))
          !
      enddo
      !
      ! Second Eckart equation
      ! 
      ieq = 3
      !
      do ix = 1,3 
         !
         ieq = ieq +1 
         !
         e_t = 0 
         do n0 = 1,trove%Natoms
            !
            do jx = 1,3 
               !
               ! Eckart
               !
               select case (trim(axis_system)) 
               !
               case default
               !
                 write (out,"('Fcalc_eckart: axis system ',a,' unknown')") trim(axis_system)
                 stop 'Fcalc_eckart - bad axis system'
               !
               case('Eckart') 
                 !
                 e_t = e_t + sum(epsil(ix,:,jx)*trove%b0(n0,:,irho) )*cartesian(n0,jx)*trove%mass(n0)
                 !
               case('Full_Eckart') 
                 !
                 e_t = e_t + sum(epsil(ix,:,jx)*cartesian(n0,:) )*cartesian(n0,jx)*trove%mass(n0)
                 !
               case('PAS') 
                 !
                 do kx = 1,3 
                    !
                    f_t = epsil(ix,jx,kx)**2
                    !
                    e_t = e_t + f_t*cartesian(n0,jx)*cartesian(n0,kx)*trove%mass(n0)
                    !
                 enddo
                 !
               end select 
            enddo
            !
         enddo
         !
         Eckart(ieq) = e_t
         !
      enddo
      !
      !
   end subroutine calc_eckart


!
! Here we initialize a polynom
!
  subroutine polynom_initialization(polynom,orders,Ncoeff,Npoints,name)
     !
     integer(ik),intent(in) :: orders,Ncoeff,Npoints
     character(len=*),intent(in) :: name
     type(FLpolynomT),pointer :: polynom
     integer(ik) :: alloc,i
       !
       polynom%Orders = orders
       polynom%Ncoeff = Ncoeff
       polynom%Npoints = Npoints
       polynom%name   = name
       allocate (polynom%field(Ncoeff,0:Npoints),polynom%iorder(Ncoeff),&
                 polynom%IndexQ(trove%Nmodes,Ncoeff),stat=alloc)
       allocate(polynom%ifromsparse(Ncoeff),stat=alloc)
       call ArrayStart(name,alloc,size(polynom%field),kind(polynom%field))
       call ArrayStart(name,alloc,size(polynom%iorder),kind(polynom%iorder))
       call ArrayStart(name//"IndexQ",alloc,size(polynom%IndexQ),kind(polynom%IndexQ))
       call ArrayStart(name//'ifromsparse',alloc,size(polynom%ifromsparse),kind(polynom%ifromsparse))
       !
       !if (alloc/=0) then
       !   write (out,"(' Error ',i9,' trying to allocate fields of polynom ',a)") alloc,trim(name)
       !   stop 'polynom_initialization, polynom - out of memory'
       !end if
       polynom%field  = 0 
       polynom%iorder = 0
       polynom%IndexQ(:,1:Ncoeff) = FLIndexQ(:,1:Ncoeff)
       forall(i = 1:Ncoeff) polynom%ifromsparse(i) = i 
       !
       if (verbose>=5) write(out,"('field ',a,' initialized')") trim(name)
       ! 
  end subroutine polynom_initialization


  subroutine FL_exclude_specific_modes(imode1,imode2)
    !
    integer(ik),intent(in)  :: imode1,imode2
    !
    integer(ik)  :: excluded_power,iterm,k(trove%Nmodes),k1,k2
    !
    type(FLpolynomT),pointer      :: fl
    !
    fl => trove%poten
    !
    fl%iorder = 0 
    !
    do iterm = 1,fl%Ncoeff
       !
       k(:) = FLIndexQ(:,iterm)
       !
       ! For the zero order case we allow only diagonal terms
       !
       ! Check if the current iterm belongs to the present calculation case
       ! another words, if poten*xi^k belongs to the current perturb. order
       !
       excluded_power = sum(k(1:imode1-1)) + sum(k(imode2+1:trove%Nmodes_e))
       !
       if (excluded_power>0) fl%iorder(iterm) = 1
       !
     enddo
     !
     ! 
     ! Vibrational part of the kinetic operator g_vib
     !
     do k1 = 1,trove%Nmodes
        do k2 = 1,trove%Nmodes
           !
           fl => trove%g_vib(k1,k2)
           !
           fl%iorder = 0
           !
           do iterm = 1,fl%Ncoeff
              !
              k(:) = FLIndexQ(:,iterm)
              !
              excluded_power = sum(k(1:imode1-1)) + sum(k(imode2+1:trove%Nmodes_e))
              !
              if (excluded_power>0.or.imode1>k1.or.k1>imode2.or.&
                                      imode1>k2.or.k2>imode2)  &
                                                fl%iorder(iterm) = 1
              !
           enddo
           ! 
        enddo
     enddo
     !
     ! Rotational and Coriolis parts
     !
     if (FLrotation) then 
       !
       do k1 = 1,3
          do k2 = 1,3
             !
             fl => trove%g_rot(k1,k2)
             !
             fl%iorder = 0
             !
             do iterm = 1,fl%Ncoeff
                !
                k(:) = FLIndexQ(:,iterm)
                !
                excluded_power = sum(k(1:imode1-1)) + sum(k(imode2+1:trove%Nmodes_e))
                !
                if (excluded_power>0) fl%iorder(iterm) = 1
                !
             enddo
             ! 
          enddo
       enddo
       !
       do k1 = 1,trove%Nmodes
          do k2 = 1,3
             !
             fl => trove%g_cor(k1,k2)
             !
             fl%iorder = 0
             !
             do iterm = 1,fl%Ncoeff
                !
                k(:) = FLIndexQ(:,iterm)
                !
                excluded_power = sum(k(1:imode1-1)) + sum(k(imode2+1:trove%Nmodes_e))
                !
                if (excluded_power>0.or.imode1>k1.or.k1>imode2) fl%iorder(iterm) = 1
                !
             enddo
             ! 
          enddo
       enddo
       !
     endif 
     !
  end subroutine FL_exclude_specific_modes




 subroutine FLinit_External_field_andrey

    !expand external field (function) components in Taylor series on xi coordinates

    integer(ik) :: npoints, nmodes, imu, iterm, imode, jmode, ipoint, maxorder, alloc, nterms, &
                   imu_maxord, irank
    integer(ik) :: ipowers(trove%nmodes),Ncoeff,icoeff,i,lambda,jterm
    real(ark)   :: step(trove%nmodes)
    type(FLpolynomT),pointer     :: fl
    real(ark),allocatable        :: extF_t(:),rho_extF(:)
    real(ark)                    :: axi_eq(trove%Nmodes),rho
    real(ark)                    :: df_t(extF%rank),factor
    real(ark)                    :: aq_eq(trove%Nmodes),f_t
    integer(ik),allocatable      :: powers_(:,:)
    character(len=cl) :: my_fmt1,my_fmt2 !format for I/O specification
    character(len=cl) :: my_fmt !format for I/O specification
    !
    call TimerStart('External')
    !
    if (trim(trove%IO_ext_coeff)=='READ') then 
      !
      call FLcheck_point_Hamiltonian('EXTERNAL_READ') 
      !
      return
      !
    endif
    !
    npoints   = trove%Npoints
    nmodes    = trove%Nmodes
    !
    if (.not.FLextF_coeffs)  return
    !
    if (job%verbose>=1) then
      write(out, '(/a)') 'FLinit_External_field_andrey/start'
      write (out,"(' External fields need ',f9.3,' Mbytes of memory (plus a bit)')") &
             real(rk*product(extF%maxord(:)),kind=rk)/(1024.0_rk**2)
    end if
    !
    allocate (trove%extF(extF%rank),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate extF-field')") alloc
        stop 'FLinit_External_field_andrey, extF-field - out of memory'
    end if
    !
    if (job%verbose>=4)  write(out, '(/1x, a)') 'generate expansion terms'
    !
    nterms = maxval(extF%nterms(:))
    !
    write(my_fmt1,'(a,i0,a)') "(/2(1x, a)",Nmodes,"(1x, i3))"
    !
    write(my_fmt2,'(a,i0,a)') "(3x, i6, 8x,",Nmodes,"(1x, i3))"
    !
    do imu = 1,extF%rank
       !
       ! Number of the expansion coefficients
       !
       nterms  = trove%RangeOrder(extF%maxord(imu))
       !
       fl => trove%extF(imu)
       call polynom_initialization(fl,extF%maxord(imu),nterms,Npoints,'extF')
       !
       if (job%verbose>=7) then
          write(out, '(/1x, a, 25x, i3/1x, a, 4x, i3/1x, a, 1x, i5)') &
          'imu', imu, 'maximum expansion degree', trove%NExtOrder,'number of expansion terms', nterms
          write(out,my_fmt1) 'iterm', 'imode->', (imode, imode = 1, nmodes)
       end if
       !
       !print expansion terms
       !
       if (job%verbose>=5) then
          do iterm = 1, trove%extF(imu)%Ncoeff
             write(out,my_fmt2) iterm, (FLIndexQ(imode,iterm), imode = 1, nmodes)
          end do
       end if
       !
    enddo
    !
    !
    step(:) = extF%fdstep(:)
    !
    aq_eq(:)= trove%chi_ref(:,0)
    !
    do imode=1,trove%Nmodes_e
       !
       axi_eq(imode) = MLcoord_direct(aq_eq(imode),3,imode)
       !
    enddo
    !
    if (job%verbose>=4) then
       write(out, '(/1x, a)') 'generate expansion coefficients'
    end if
    !
    !component of external fied for which maximum expansion degree occured
    !
    imu_maxord = maxloc(trove%extF(1:extF%rank)%Ncoeff, dim=1)
    !
    !loop over nonrigid coordinate points
    !
    maxorder = trove%NExtOrder
    !
    do ipoint = 0, npoints
       !
       write(out,"('i = ',i8)") ipoint
       !
       if (manifold/=0) then 
          !
          aq_eq(:) = trove%chi_ref(:,ipoint)
          !
          do imode = 1,Nmodes
            !
            axi_eq(imode) = MLcoord_direct(aq_eq(imode),3,imode)
            !
          enddo
          !
       endif 
       !
       !loop over expansion terms
       !
       !omp parallel do private(iterm,ipowers,factor,imode,jmode,df_t) schedule(guided)
       do iterm = 1, trove%extF(imu_maxord)%Ncoeff
          !
          ipowers(1:nmodes) =  FLIndexQ(1:nmodes,iterm)
          !
          factor = 1.0_ark
          do imode = 1,Nmodes
             do jmode = 1,ipowers(imode)
                factor = factor*real(jmode,kind=ark)
             enddo
          enddo 
          !
          !here we compute derivatives for all components of external field...
          !
          df_t(1:extF%rank) = FL_fdf(nmodes,axi_eq(1:nmodes),ipoint,dms4xi,maxorder,ipowers,step,extF%rank)
          !
          !...but store only those which correspond to the expansion
          !of order defined for each component of external field
          !
          do irank = 1, extF%rank
             if (iterm <= trove%extF(irank)%Ncoeff) trove%extF(irank)%field(iterm,ipoint) = df_t(irank)/factor
          end do
          !
       end do
       !omp end parallel do
       !
    end do
    !
    !print expansion coefficients and expansion points
    !
    if (job%verbose>=4) then
       !
       Ncoeff = maxval(trove%extF(1:extF%rank)%Ncoeff, dim=1) 
       !
       write(my_fmt,'(a,i0,a)') "(/1x, a, 1x, a, 1x, a,",Ncoeff,"(10x, i3))"
       !
       write(out,my_fmt) 'ipoint', 'imu', 'iterm->', (iterm, iterm = 1, trove%extF(imu_maxord)%Ncoeff)
       !
       write(my_fmt,'(a,i0,a)') "(1x, i6, 2x, i2, 8x,",Ncoeff,"(1x, es12.4))"
       !
       do irank = 1, extF%rank
         do ipoint = 0, npoints
             !
             Ncoeff = trove%extF(imu_maxord)%Ncoeff
             !
             write(out,my_fmt) ipoint, irank, (trove%extF(irank)%field(iterm,ipoint), iterm = 1, trove%extF(irank)%Ncoeff)
          end do
       end do
    end if
    !
    if (job%verbose>=1) then
       write(out, '(/a)') 'FLinit_External_field_andrey/end'
    endif
    !
    call TimerStop('External')
    !
    call TimerReport
    !
 end subroutine FLinit_External_field_andrey




 subroutine FLinitilize_Potential_Andrey

 !expand external field (function) components in Taylor series on xi coordinates

    integer(ik) :: npoints, nmodes, imu, iterm, imode, jmode, ipoint, maxorder, alloc, nterms, &
                   imu_maxord, irank
    integer(ik) :: ipowers(trove%nmodes),Ncoeff,icoeff,i,lambda,jterm
    real(ark)   :: step(trove%nmodes)
    type(FLpolynomT),pointer     :: fl
    real(ark)                    :: axi_eq(trove%Nmodes),rho
    real(rk)                     :: df_t(1),factor
    real(ark)                    :: aq_eq(trove%Nmodes),f_t
    character(len=cl)            :: my_fmt !format for I/O specification
    !
    !
    if (job%verbose>=2) write(out,"(/'Expansion of the potential function...')")   
    !
    ! If the potentil function has been stored we can just read it from the hard disk and leave...
    !
    if (trim(trove%IO_hamiltonian)=='READ'.or.&
        trim(trove%IO_potential)=='READ') then 
        !
        call FLcheck_point_Hamiltonian('POTENTIAL_READ') 
        !
        return 
        !
    endif
    !
    call TimerStart('Potential')
    ! 
    ! Parameters for the internal use 
    !
    Nmodes  = trove%Nmodes
    Npoints = trove%Npoints
    !
    ! The potential energy function within the calculations must be given 
    ! in the normal coordinates. However at the input it is suposed to be 
    ! in the geometrically defined coordinates. 
    ! Therefore 
    ! To initialize potential energy function ("poten"-object) we perform 
    ! 1. Coordinate transformation from GDC to normal coordinates 
    ! 2. Calculate derivatives (normal force constants) using finite difference method
    !
    ! Allocation of the poten-field
    !
    allocate (trove%poten,stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i9,' trying to allocate poten-field')") alloc
        stop 'FLinitilize_Potential, poten-field - out of memory'
    end if
    !
    fl => trove%poten
    call polynom_initialization(fl,trove%NPotOrder,trove%RangeOrder(trove%NPotOrder),Npoints,'poten')
    !
    ! The potential function is expanded at its minimum, which is zero for the normal coordinates
    !
    aq_eq(:)= trove%chi_ref(:,0)
    !
    step(:) = trove%fdstep(:)
    !
    do imode=1,trove%Nmodes_e
       !
       axi_eq(imode) = MLcoord_direct(aq_eq(imode),2,imode)
       !
    enddo
    !
    if (job%verbose>=4) then
       write(out, '(/1x, a)') 'generate expansion coefficients'
    end if
    !
    maxorder = trove%NPotOrder
    !
    do ipoint = 0, npoints
       !
       write(out,"('i = ',i8)") ipoint
       !
       if (manifold/=0) then 
          aq_eq(Nmodes) = trove%rho_i(ipoint)
          axi_eq(Nmodes) = MLcoord_direct(aq_eq(Nmodes),2,Nmodes)
          !trove%chi_ref(Nmodes,ipoint) = aq_eq(Nmodes)
       endif 
       !
       !loop over expansion terms
       !
       !omp parallel do private(iterm,ipowers,factor,imode,df_t) schedule(guided)
       do iterm = 1, trove%poten%Ncoeff
          !
          ipowers(1:nmodes) =  FLIndexQ(1:nmodes,iterm)
          !
          factor = 1.0_ark
          do imode = 1,Nmodes
             do jmode = 1,ipowers(imode)
                factor = factor*real(jmode,kind=ark)
             enddo
          enddo 
          !
          !here we compute derivatives for all components of external field...
          !
          df_t(1:1) = FL_fdf(nmodes,axi_eq(1:nmodes),ipoint,FLpoten4xi,maxorder,ipowers,step,1)
          !
          !...but store only those which correspond to the expansion
          !of order defined for each component of external field
          !
          if (iterm <= trove%poten%Ncoeff) trove%poten%field(iterm,ipoint) = df_t(1)/factor
          !
       end do
       !omp end parallel do
       !
    end do
    !
    !print expansion coefficients and expansion points
    !
    if (job%verbose>=5) then
       !
       Ncoeff = trove%poten%Ncoeff
       !
       write(my_fmt,'(a,i0,a)') "(/1x, a, 1x, a, 1x, a,",Ncoeff,"(10x, i3))"
       !
       write(out,my_fmt) 'ipoint', 'imu', 'iterm->', (iterm, iterm = 1, trove%poten%Ncoeff)
       !
       write(my_fmt,'(a,i0,a)') "(1x, i6, 8x,",Ncoeff,"(1x, es12.4))"
       !
       do ipoint = 0, npoints
           !
           write(out,my_fmt) ipoint, (trove%poten%field(iterm,ipoint), iterm = 1, trove%poten%Ncoeff)
        end do
    end if
    !
    if (job%verbose>=2) write(out,"(/'...done!')")   
    !
 end subroutine FLinitilize_Potential_Andrey



 function FL_fdf(num_var, var_val, ivar, get_func, maxorder, der_ord, hh, rank)

 !compute partial derivative by the finite-difference technique

    interface
       subroutine get_func(rank,nmodes, ipoint, xi, func)
          use accuracy, only: ik, ark
          use moltype, only: extF
          integer(ik),intent(in) :: rank,nmodes,ipoint
          real(ark),intent(in)   :: xi(nmodes)
          real(ark),intent(out)  :: func(rank)
       end subroutine get_func
    end interface
    !
    integer(ik), intent(in) :: num_var, ivar, maxorder, der_ord(num_var), rank
    real(ark), intent(in)   :: var_val(num_var), hh(num_var)
    real(ark)               :: FL_fdf(rank)
    !
    real(ark)               :: h(num_var), x(2 ** maxorder + 1, num_var), arg(num_var), f(rank), frac
    integer(ik)             :: sgn(2 ** maxorder + 1, num_var), num_arg(num_var), no_arg(num_var)
    integer(ik)             :: sign, i, j, k, l, n
    real(ark)               :: func(rank)
    !
    !
    h = hh
    !
    do i = 1, num_var
       x(1, i) = var_val(i)
       x(2, i) = x(1, i)
       sgn(1, i) = 1
       sgn(2, i) = 1
       if (der_ord(i) == 0) then
          sgn(1, i) = -1
          h(i) = 0
       end if
       do j = 1, der_ord(i) - 1
          n = 2 ** (j + 1) + 1
          do k = 1, 2 ** j
             x(k, i) = x(k, i) + (-1) ** k * h(i)
             x(n - k, i) = x(k, i)
             sgn(k, i) = sgn(k, i) * (-1) ** k
             sgn(n - k, i) = sgn(k, i)
          end do
       end do
       num_arg(i) = 2 ** der_ord(i)
       do k = 1, num_arg(i)
          x(k, i) = x(k, i) + (-1) ** k * h(i)
          sgn(k, i) = sgn(k, i) * (-1) ** k
       end do
    end do
    !
    n = 1
    frac = 1.0_ark
    do i = 1, num_var
       n = n * num_arg(i)
       frac = frac * (2.0_ark * h(i)) ** der_ord(i)
    end do
    !
    no_arg = 1
    f = 0.0_ark
    do i = 1, n
       sign = 1
       do j = num_var, 2, -1
          if (no_arg(j) > (num_arg(j))) then
             no_arg(j) = 1
             no_arg(j - 1) = no_arg(j - 1) + 1
          end if
          arg(j) = x(no_arg(j), j)
          sign = sign * sgn(no_arg(j), j)
       end do
       arg(1) = x(no_arg(1), 1)
       sign = sign * sgn(no_arg(1), 1)
       !
       call get_func(rank,num_var, ivar, arg(1:num_var),func(1:rank))
       !
       f(1:rank) = f(1:rank) + func(1:rank) * sign
       !
       no_arg(num_var) = no_arg(num_var) + 1
    end do
    !
    FL_fdf(1:rank) = f(1:rank) / frac
    !
 end function FL_fdf



 subroutine dms4xi(rank,nmodes, ipoint, xi, mu_xyz)

 !return dipole moment cartesian component curr_imu for xi coordinates specified;

    integer(ik), intent(in) :: rank,nmodes,ipoint
    real(ark), intent(in)   :: xi(nmodes)
    real(ark),intent(out)   :: mu_xyz(1:rank)
    !
    real(ark)               :: r(molec%ncoords),xyz(molec%natoms, 3),chi(nmodes)
    integer(ik)             :: natoms,nlocals,iatom,icart,imode
    logical                 :: dir
    character(len=cl)       :: my_fmt !format for I/O specification
    !
    if (verbose>=6) write(out,"('dms4xi/start')")    
    !
    !if (rank/=3.and.rank/=molec%ncoords) then 
    !  write(out,"('dms4xi: is working only  for rank =3 or ncoords (not ',i5,') components ext-function')") extF%rank
    !  stop 'illegal rank for use in dms4xi'
    !endif
    !
    if (verbose>=6) then 
       write(out, '(/a, 1x, i5)') 'dms4xi run for ipoint', ipoint
    endif
    !
    !strategy:
    !(here pqr and xyz are space- and body-fixed coordinates of atoms)
    !
    !xi => xyz => r => (user-defined function) => pqr => (xyz = cosmat * pqr) => cosmat =\
    !             |                                                                       => mu_xyz = cosmat * mu_pqr
    !             r => (user-defined function) => mu_pqr ================================/
    !
    natoms  = molec%natoms
    nlocals = molec%ncoords
    !
    do imode = 1,trove%Nmodes
      chi(imode) = MLcoord_invert(xi,3,imode) 
    enddo
    !
    !chi = xi
    !
    do iatom = 1,trove%Natoms
       do icart = 1,3
         !
         xyz(iatom,icart) = trove%b0(iatom,icart,ipoint) + sum( ( trove%Amatrho(iatom,icart,1:trove%Nmodes_e,ipoint) )* &
                                                                ( chi(1:trove%Nmodes_e)-trove%chi_ref(1:trove%Nmodes_e,ipoint)) )
         !
       enddo
    enddo
    !
    if (trove%internal_coords=='LOCAL') then 
       !
       dir = .false.
       !
       r = MLcoordinate_transform_func(chi,size(r),dir)
       !
       xyz = 0
       !
    else
       !
       call FLfromcartesian2local(xyz,r)
       !
    endif
    !
    if (verbose>=6) then 
       write(my_fmt,'(a,i0,a,i0,a)') "(1x, a/",Nmodes,"(1x, es16.8)/1x, a/",nlocals,"(1x, es16.8))"
       write(out,my_fmt) 'xi', xi(1:nmodes), 'r', r(1:nlocals)
       write(my_fmt,'(a,i0,a)') "(1x, a,",Nmodes,"(/3(1x, es16.8)))"
       write(out,my_fmt) 'xyz', (xyz(iatom, 1:3), iatom = 1, natoms)
    endif
    !
    call MLextF_func(rank,molec%ncoords,molec%natoms,r,xyz,mu_xyz)
    !
    if (verbose>=6) write(out,"('dms4xi/end')")
    !
 end subroutine dms4xi



 recursive subroutine dms4chi(ipoint, dchi, mu_xyz)

 !return dipole moment cartesian component curr_imu for chi coordinates specified;

    integer(ik), intent(in) :: ipoint
    real(ark), intent(in)   :: dchi(:)
    real(ark),intent(out)   :: mu_xyz(:)
    !
    real(ark)               :: r(molec%ncoords),xyz(molec%natoms, 3)
    integer(ik)             :: natoms,nlocals,iatom,icart,rank,nmodes
    character(len=cl) :: my_fmt !format for I/O specification
    !
    if (verbose>=6) write(out,"('dms4chi/start')")    
    !
    !if (rank/=3.and.rank/=molec%ncoords) then 
    !  write(out,"('dms4chi: is working only  for rank =3 or ncoords (not ',i5,') components ext-function')") extF%rank
    !  stop 'illegal rank for use in dms4chi'
    !endif
    !
    if (verbose>=6) then 
       write(out, '(/a, 1x, i5)') 'dms4chi run for ipoint', ipoint
    endif
    !
    !strategy:
    !(here pqr and xyz are space- and body-fixed coordinates of atoms)
    !
    !chi => xyz => r => (user-defined function) => pqr => (xyz = cosmat * pqr) => cosmat =\
    !             |                                                                       => mu_xyz = cosmat * mu_pqr
    !             r => (user-defined function) => mu_pqr ================================/
    !
    natoms  = molec%natoms
    nlocals = molec%ncoords
    !
    !chi = xi
    !
    do iatom = 1,trove%Natoms
       do icart = 1,3
         !
         xyz(iatom,icart) = trove%b0(iatom,icart,ipoint) + &
                            sum(trove%Amatrho(iatom,icart,1:trove%Nmodes_e,ipoint)*dchi(1:trove%Nmodes_e))
         !
       enddo
    enddo
    !
    call FLfromcartesian2local(xyz,r)
    !
    nmodes = size(dchi)
    !
    if (verbose>=6) then 
       !
       write(my_fmt,'(a,i0,a,i0,a)') "(1x, a/",Nmodes,"(1x, es16.8)/1x, a/",nlocals,"(1x, es16.8))"
       write(out,my_fmt) 'xi', dchi(1:nmodes), 'r', r(1:nlocals)
       write(my_fmt,'(a,i0,a)') "(1x, a,",Nmodes,"(/3(1x, es16.8)))"
       write(out,my_fmt) 'xyz', (xyz(iatom, 1:3), iatom = 1, natoms)
       !
       !write(out, '(1x, a/<nmodes>(1x, es16.8)/1x, a/<nlocals>(1x, es16.8) )') 'xi', dchi(1:nmodes), 'r', r(1:nlocals)
       !write(out, '(1x, a, <natoms>(/3(1x, es16.8)))') 'xyz', (xyz(iatom, 1:3), iatom = 1, natoms)
    endif
    !
    rank = size(mu_xyz,dim=1)
    !
    call MLextF_func(rank,molec%ncoords,molec%natoms,r,xyz,mu_xyz)
    !
    if (verbose>=6) write(out,"('dms4chi/end')")
    !
 end subroutine dms4chi



 recursive subroutine FLpoten4xi(rank,nmodes, ipoint, xi, pot) 

 !return poten_linear for xi coordinates specified;

    integer(ik), intent(in) :: rank,nmodes, ipoint
    real(ark), intent(in)   :: xi(nmodes)
    real(ark),intent(out)   :: pot(1:rank)
    !
    real(ark)               :: chi,x(-2:2),f(-2:2),df_t
    integer(ik)             :: irho_eq,irho,i
    !
    !
    if (rank/=1) then 
      write(out,"('FLpoten4xi: is working only with a 1one (not ',i5,') component function, not')") extF%rank
      stop 'illegal rank for use in FLpoten4xi'
    endif
    !
    if (verbose>=6) then 
       write(out, '(/a, 1x, i5)') 'FLpoten4xi run for ipoint', ipoint
    endif
    !
    ! chi = MLcoord_invert(xi,2,nmodes)
    !
    pot(1) = FLpoten_linearized(xi,ipoint)
    !
    !irho_eq = 0
    !
    !if (ipoint>=2.and.ipoint<=trove%npoints-2) then 
    !  !
    !  irho_eq = mod(nint( ( chi-trove%rho_border(1) )/(trove%rhostep),kind=ik ),trove%npoints)
    !  !
    !  do i=-2,2 
    !   !
    !   x(i) = trove%rho_border(1)+(irho_eq+i)*trove%rhostep
    !   f(i) = FLpoten_linearized(xi,irho_eq+i)
    !   !
    !  enddo
    !  !
    !  call polintark(x(-2:2),f(-2:2),chi,pot(1),df_t)
    !  !
    !else
    !  !
    !  pot(1) = FLpoten_linearized(xi,irho_eq)
    !  !
    !endif
    !
 end subroutine FLpoten4xi


  subroutine FL_rotation_energy_surface
    !
    integer(ik)         :: jval,iphi,itheta,iter,ivar(2*trove%nmodes+2),Nmodes,parmax
    real(rk)            :: theta,phi,xi(2*trove%nmodes+2),rjacob(2*trove%nmodes+2),hess(2*trove%nmodes+2,2*trove%nmodes+2)
    real(rk)            :: ssq_old,stability,ssq,thetastep,phistep,taustep,ssq_best
    real(rk)            :: Energy,tempx,Fright,Fleft,Fplusminus,Fminusplus,am(2*trove%nmodes+2,2*trove%nmodes+2),&
                           bl(2*trove%nmodes+2,1)
    real(rk)            :: dx(2*trove%nmodes+2),t_xi(2*trove%nmodes+2),Tsing(2*trove%nmodes+2)
    integer(ik)         :: icol,i,irow,icolumn,numpar,ncol,ncol1,ncol2,i1,i2,rank,iwork, info
    real(rk)            :: stability_best,deltax(2*trove%nmodes+2)
    real(rk)            :: r_na(trove%natoms,3),largest,res_,tau,smallest,r(trove%Ncoords)
    integer(ik)         :: iphi_,itheta_,itau_,ncoords,iatom,naught_at(2),iphi_s,itheta_s,itau_s,alloc_p,Ntheta1

    !
    integer(ik)         :: Ntheta,Nphi,itmax,Ntau,itau
    !character(len=cl)     :: flag_res = 'OPTIM_COORDS'
    character(len=cl)     :: flag_res,comment
    !character(len=cl)     :: flag_res = 'THETA_TAU'
    !character(len=cl)     :: flag_res = 'GLOBAL_SEARCH'
    !
    real(rk) :: work((2*trove%Nmodes+2)*50)

    real(rk),allocatable  :: RES(:,:,:),xi_(:,:,:,:),r_(:,:,:),flag(:,:)
    character(len=cl)     :: my_fmt !format for I/O specification

      !
      jval = trove%jmax
      !
      Nmodes = trove%nmodes
      ncoords = trove%Ncoords
      !
      parmax = 2*trove%Nmodes+2
      !
      iwork= parmax*50
      !
      naught_at = 0 
      !
      if ( abs( trove%rho_border(2)-trove%rho_border(1)-2.0_ark*pi )<0.01.or.abs( trove%rho_border(2)-trove%rho_border(1)-&
           4.0_ark*pi )<0.01 ) then 
        naught_at = 1
      endif 
      !
      Ntheta = analysis%res%Ntheta
      Nphi = analysis%res%Nphi
      Ntau = analysis%res%Ntau
      flag_res = analysis%res%type
      !
      select case (trim(flag_res))
        !
      case ('THETA_TAU') 
         Nphi = 0
      case ('THETA_PHI') 
         Ntau = 0
      case ('PHI_TAU') 
         Ntheta = 0
      case ('GLOBAL_SEARCH') 
         Ntau = 0
      end select
      !
      thetastep = 0
      phistep = 0
      taustep = 0 
      !
      if (Ntheta/=0) thetastep = pi/Ntheta
      if (Nphi/=0) phistep   = 2.0_rk*pi/Nphi
      if (Ntau/=0) taustep   = 2.0_rk*pi/Ntau
      !
      Ntheta1 = 0
      !
      if (analysis%res%theta1>small_) then
         Ntheta1 = int(analysis%res%theta1/thetastep)
      endif
      !
      itmax = analysis%res%itermax
      !
      write(out,"(/'  Rotational energy surface'/)")
      !
      itheta_ = 0 ; iphi_ = 0 ;  itau_ = 0 ; largest = 0
      itheta_s= 0 ; iphi_s= 0 ;  itau_s= 0 ; smallest= 0
      !

      stability =  1.e10
      stability_best = 100.0_rk*sqrt(small_)
      ssq_best = 1e-16
      !
      allocate(RES(0:Ntheta,0:Nphi,0:Ntau),stat=info)
      call ArrayStart('RES',info,size(RES),kind(RES))
      allocate(xi_(parmax,0:Ntheta,0:Nphi,0:Ntau),stat=info)
      call ArrayStart('RES',info,size(xi_),kind(xi_))
      allocate(r_(trove%Ncoords,0:Nphi,0:Ntau),stat=info)
      call ArrayStart('RES',info,size(r_),kind(r_))
      allocate(flag(0:Nphi,0:Ntau),stat=info)
      call ArrayStart('RES',info,size(r_),kind(r_))
      !     
      !omp parallel private(work,alloc_p) shared(RES,xi_)
      !allocate(work(iwork),stat=alloc_p)
      !if (alloc_p/=0) then
      !   write (out,"(' Error ',i9,' trying to allocate array pot_points')") alloc_p
      !   stop 'FLinitilize_Potential, pot_points  - out of memory'
      !end if
      !
      do itheta=Ntheta1,Ntheta,1
        !
        theta=min(itheta*thetastep,pi)
        !
        !
        !$omp  parallel do private(iphi,phi,itau,tau,xi,ivar,numpar,&
        !$omp& rjacob,iter,ssq_old,ssq,dx,stability,ncol,i,r_na,r,Energy,deltax,t_xi,ncol1,i1,ncol2,i2,Fright,Fleft,Hess,&
        !$omp& Fplusminus,irow,icolumn,am,bl,Tsing,rank,info,comment) shared(RES,xi_,r_,flag) schedule(guided) 
        do iphi=0,Nphi,1
           !
           !if (job%verbose>=5) write(out,"('itheta,iphi = ',2i8)") itheta,iphi
           !
           phi=min(iphi*phistep,2.d0*pi)
           !
           do itau=0,Ntau,1
             !
             tau=min(itau*taustep,2.d0*pi)
             !
             xi(1:Nmodes) = trove%chi_ref(1:Nmodes,0)
             xi(Nmodes+1:2*Nmodes) = 0 
             xi(2*Nmodes+1:2*Nmodes+2) = (/theta,phi/)
             !
             if (itheta>0) xi(1:Nmodes) = xi_(1:Nmodes,itheta-1,iphi,itau) 
             !
             ivar = 1 
             !
             select case (flag_res)
               !
             case ('THETA_TAU')
               !
               xi(Nmodes)= tau
               !
               phi = analysis%res%phi/180.0_rk*pi
               !
               xi(2*Nmodes+1:2*Nmodes+2) = (/theta,phi/)
               ivar(Nmodes:2*Nmodes+2) = 0
               !
             case ('PHI_TAU')
               !
               xi(Nmodes)= tau
               !
               theta = analysis%res%theta/180.0_rk*pi
               !
               xi(2*Nmodes+1:2*Nmodes+2) = (/theta,phi/)
               ivar(Nmodes:2*Nmodes+2) = 0
               !
             case ('THETA_PHI')
               !
               tau = analysis%res%tau/180.0_rk*pi
               tau = max(trove%chi_ref(Nmodes,0),tau)
               tau = min(tau,trove%chi_ref(Nmodes,trove%npoints))
               !
               xi(Nmodes)= tau
               !
               xi(2*Nmodes+1:2*Nmodes+2) = (/theta,phi/)
               ivar(Nmodes:2*Nmodes+2) = 0
               !
             case ('THETA_PHI_TAU')
               !
               xi(Nmodes) = tau
               xi(2*Nmodes+1:2*Nmodes+2) = (/theta,phi/)
               ivar(Nmodes:2*Nmodes+2) = 0
               !
             case ('OPTIM_COORDS')
               !
               ivar(Nmodes+1:2*Nmodes+2) = 0
               !
               tau = analysis%res%tau/180.0_rk*pi
               !
               tau = max(trove%chi_ref(Nmodes,0),tau)
               tau = min(tau,trove%chi_ref(Nmodes,trove%npoints))
               !
               xi(Nmodes)= tau
               !
             case ('GLOBAL_SEARCH')
               !
               ivar(Nmodes+1:2*Nmodes) = 0
               !
               if (Ntheta==0) theta = analysis%res%theta/180.0_rk*pi
               !
               tau = analysis%res%tau/180.0_rk*pi
               !
               xi(Nmodes)= tau
               !
             end select
             !
             numpar = sum(ivar)
             !
             rjacob = 0 
             !hess = 0 
             iter = 0
             ssq_old = 1.e10
             ssq    =  1.e10
             dx = 0
             !
             stability = 1e8
             !
             ! start optimization here:
             !
             outer_loop: & 
             do while( iter<itmax .and. stability>stability_best ) 
                !
                !if (job%verbose>=5) write(out,"('iter = ',i8)") iter
                !
                iter = iter + 1
                ssq=0
                !   
                ! Update the pot. parameters to the new values 
                !
                ncol=0
                do i=1,parmax
                  if (ivar(i) /= 0) then
                     ncol=ncol+1
                     xi(i)=xi(i)+dx(ncol)*0.1
                  endif
                enddo
                !
                ! Caclulate the function 
                !
                !if (job%verbose>=5) write(out,"('classic_hamilt')") 
                !
                Energy=classic_hamilt(jval,xi,r,r_na)
                !
                !if (job%verbose>=5) write(out,"('...done!')") 
                !
                !if (job%verbose>=5) write(out,"('hessian:')") 
                !
                rjacob = 0
                !
                deltax=trove%fdstep(1) 
                if (manifold/=0) deltax(Nmodes) = trove%rhostep
                !
                if (itmax/=0 ) then
                  t_xi = xi
                  ncol1=0
                  do  i1=1,parmax
                    ! 
                    if (ivar(i1) /= 0) then
                      !
                      ncol1=ncol1+1
                      !if (job%verbose>=5) write(out,"('ncol1 = ',i8)") ncol1
                      !
                      ncol2=0
                      do  i2=1,i1
                        !
                        if (ivar(i2) /= 0) then
                          !
                          ncol2=ncol2+1
                          !
                          !if (job%verbose>=5) write(out,"('ncol2 = ',i8)") ncol2
                          !
                          t_xi = xi
                          !
                          t_xi(i1) = t_xi(i1)+deltax(i1)
                          t_xi(i2) = t_xi(i2)+deltax(i2)
                          Fright = classic_hamilt(jval,t_xi)
                          !
                          t_xi = xi
                          !
                          t_xi(i1) = t_xi(i1)-deltax(i1)
                          t_xi(i2) = t_xi(i2)-deltax(i2)
                          Fleft = classic_hamilt(jval,t_xi)
                          !
                          if (ncol1==ncol2) then 
                            !
                            Hess(ncol1,ncol2)=(Fright+Fleft-2.0d0*Energy)/(deltax(i1)*deltax(i2)*4.0_rk)
                            !
                          else
                            !
                            t_xi = xi
                            !
                            t_xi(i1) = t_xi(i1)+deltax(i1)
                            t_xi(i2) = t_xi(i2)-deltax(i2)
                            Fplusminus = classic_hamilt(jval,t_xi)
                            !
                            t_xi = xi
                            !
                            t_xi(i1) = t_xi(i1)-deltax(i1)
                            t_xi(i2) = t_xi(i2)+deltax(i2)
                            Fminusplus = classic_hamilt(jval,t_xi)
                            !
                            Hess(ncol1,ncol2)=(Fright+Fleft-Fplusminus-Fminusplus)/(deltax(i1)*deltax(i2)*4.0_rk)
                            !
                          endif
                          !
                          !if (ncol1==2*Nmodes+1.xor.ncol2==2*Nmodes+1) Hess(ncol1,ncol2) = -Hess(ncol1,ncol2)
                          !if (ncol1==2*Nmodes+2.xor.ncol2==2*Nmodes+2) Hess(ncol1,ncol2) = -Hess(ncol1,ncol2)
                          !
                          Hess(ncol2,ncol1)=Hess(ncol1,ncol2)
                          !
                          if (ncol1==ncol2) then 
                           !
                           rjacob(ncol1)=(Fright-Fleft)/(4*deltax(i1))
                           !
                           if (ncol1==2*Nmodes+1.or.ncol1==2*Nmodes+2) rjacob(ncol1) = -rjacob(ncol1)
                           !
                          endif 
                          !
                        endif
                      enddo ! --- ncol
                    endif
                  enddo ! --- ncol
                endif
                !
                !if (job%verbose>=5) write(out,"('...done!')") 
                !
                ssq = energy 
                !
                if (itmax.ne.0) then
                  ! form A matrix 
                  do irow=1,numpar         !==== row-...... ====!
                    do icolumn=1,irow      !==== column-....====!
                        am(irow,icolumn) = hess(irow,icolumn)
                        am(icolumn,irow) = am(irow,icolumn)
                    enddo
                  enddo
                  ! form B matrix 
                  do irow=1,numpar       !==== row-...... ====!
                    bl(irow,1)= rjacob(irow)
                  enddo   
                  !
                  ! Solve the set of linear equations 
                  !
                  !call linur(numpar,parmax,al,bl,dx,ierror)
                  !
                  !if (job%verbose>=5) write(out,"('dgelss')") 
                  !
                  call dgelss(numpar,numpar,1,am(1:numpar,1:numpar),numpar,bl(1:numpar,1),numpar,Tsing(1:numpar),&
                              -1.0d-12,rank,work,iwork,info)
                  !
                  !if (job%verbose>=5) write(out,"('...done!')") 
                  !
                  if (info/=0) then
                    write(6,"('dgelss:error',i7)") info
                    stop 'dgelss'
                  endif
                  !
                endif 
                !
                dx(1:numpar) = -bl(1:numpar,1)
                !
                !dx(2*Nmodes+1:2*Nmodes+2) = bl(:,1)
                !
                if (ssq>sqrt(small_)) then 
                  !
                  stability=abs( (ssq-ssq_old)/ssq )
                  ssq_old=ssq
                  !                    
                else   
                  ! 
                  ssq_old=ssq
                  !
                endif 
                !
             enddo  outer_loop ! --- iter
             !
             comment = ''
             flag(iphi,itau) = 1
             !
             if (iter==itmax) then
                !
                comment = 'Not found'
                flag(iphi,itau) = 0
                !
                !write(out,"('FL_rotation_energy_suface: could not find solution after ',i8,' iterations with stability = ',g14.6)") iter, stability
                !stop 'FL_rotation_energy_surface: could not find solution'
             endif 
             !
             RES(itheta,iphi,itau)=classic_hamilt(jval,xi,r,r_na)
             !
             xi_(:,itheta,iphi,itau) = xi
             r_(:,iphi,itau) = r
             !
             !
          enddo
          !
        enddo
        !$omp end parallel do
        !
        ! print out RES
        !
        write(my_fmt,'(a,i0,a)') "(2f12.4,4x,g16.8,2x,",Ncoords,"g13.6,2x,a10)"
        !
        do iphi=0,Nphi,1
           !
           !if (job%verbose>=5) write(out,"('itheta,iphi = ',2i8)") itheta,iphi
           !
           phi=min(iphi*phistep,2.d0*pi)
           !
           do itau=0,Ntau,1
             !
             tau=min(itau*taustep,2.d0*pi)
             r(:) = r_(:,iphi,itau)
             comment = ''
             if (flag(iphi,itau)==0) comment = 'Not found'
             !
             write(out,my_fmt) theta*180.0_rk/pi,phi*180.0_rk/pi,RES(itheta,iphi,itau),r(1:ncoords),comment
             !
           enddo
           !
        enddo
        !
      enddo
      !
      !deallocate(work)
      !omp end parallel
      !
      do itheta=Ntheta1,Ntheta,1
        !
        theta=min(itheta*thetastep,pi)
        !
        do iphi=0,Nphi,1
           !
           !if (job%verbose>=5) write(out,"('itheta,iphi = ',2i8)") itheta,iphi
           !
           phi=min(iphi*phistep,2.d0*pi)
           !
           do itau=0,Ntau,1
             !
             tau=min(itau*taustep,2.d0*pi)
             !
             if (RES(itheta,iphi,itau)>largest) then 
               !
               largest = RES(itheta,iphi,itau)
               itheta_ = itheta
               iphi_ = iphi
               itau_ = itau
               !
             endif 
             !
             if (RES(itheta,iphi,itau)<smallest) then 
               !
               smallest = RES(itheta,iphi,itau)
               itheta_s = itheta
               iphi_s = iphi
               itau_s = itau
               !
             endif 

             !
          enddo
          !
        enddo
        !
      enddo
      !
      RES_ = classic_hamilt(jval,xi_(:,itheta_,iphi_,itau_),r,r_na)
      theta=xi_(2*Nmodes+1,itheta_,iphi_,itau_) ! min(itheta_*thetastep,pi)
      phi=xi_(2*Nmodes+2,itheta_,iphi_,itau_) ! min(iphi_*phistep,2.d0*pi)
      write(out,"(/'  Highest point of RES:')")
      !
      write(my_fmt,'(a,i0,a)') "(/a,2f12.4,4x,f16.8,2x,",Ncoords,"f12.6)"
      write(out,my_fmt) 'st. point: ',theta*180.0_rk/pi,phi*180.0_rk/pi,RES_,r(1:ncoords)
      write(out,"(/i6/)") trove%natoms
      !
      do iatom = 1,trove%natoms
        !
        write(out,"('X',3x,3f14.8)") r_na(iatom,:)
        !
      enddo
      !
      RES_ = classic_hamilt(jval,xi_(:,itheta_s,iphi_s,itau_s),r,r_na)
      theta=xi_(2*Nmodes+1,itheta_s,iphi_s,itau_s) ! min(itheta_*thetastep,pi)
      phi=xi_(2*Nmodes+2,itheta_s,iphi_s,itau_s) ! min(iphi_*phistep,2.d0*pi)
      write(out,"(/'  Lowest poit of RES:')")
      !
      write(my_fmt,'(a,i0,a)') "(/a,2f12.4,4x,f16.8,2x,",Ncoords,"f12.6)"
      !
      write(out,my_fmt) 'min. point: ',theta*180.0_rk/pi,phi*180.0_rk/pi,RES_,r(1:ncoords)
      !
      do iatom = 1,trove%natoms
        !
        write(out,"('X',3x,3f14.8)") r_na(iatom,:)
        !
      enddo
      !
      deallocate(RES,xi_,r_,flag)
      call arraystop('RES')
      !
    contains 

    recursive function classic_hamilt(jval,xi,r,r_na) result (H)

     integer(ik),intent(in) :: jval
     real(rk),intent(in)    :: xi(:)
     real(rk),optional,intent(out) :: r_na(trove%natoms,3),r(trove%Ncoords)
     real(rk)  :: H,pot_t,gvib_t,grot_t,gcor_t,g(trove%Nmodes),f(trove%Nmodes)
     real(rk)  :: p(trove%nmodes),theta,phi,jrot(3)
     real(ark) :: r_na_t(trove%natoms,3),r_t(trove%Ncoords),chi(trove%nmodes),chi0(trove%nmodes)
     integer(ik)  :: Nmodes,n1,ix,k1,k2,irho,imode,iterm,k(trove%nmodes)
     type(FLpolynomT),pointer    :: fl
        !
        Nmodes = trove%Nmodes
        !
        chi(1:Nmodes) = xi(1:Nmodes)
        p(1:Nmodes) = xi(Nmodes+1:2*Nmodes)
        theta = xi(2*Nmodes+1)
        phi = xi(2*Nmodes+2)

        Jrot(1) = Jval*cos(phi)*sin(theta)
        Jrot(2) = Jval*sin(phi)*sin(theta)
        Jrot(3) = Jval*         cos(theta)
        !
        irho = 0 
        !
        if (manifold/=0) then 
          irho = nint( ( chi(Nmodes)-trove%rho_border(1) )/trove%rhostep,kind=ik )
          !
          if (all(naught_at(:)==1)) then 
            !
            chi(Nmodes) = mod(chi(Nmodes)+2.0_ark*pi,2.0_ark*pi)
            !
            irho = nint( ( chi(Nmodes)-trove%rho_border(1) )/trove%rhostep,kind=ik )
            !
            irho = mod(irho+trove%npoints,trove%npoints)
            !
          endif 
          !
          !do v = -Nr_t,Nr_t
          !  r_t(v)=trove%borders(1) + trove%rhostep*real(min(max(irho_t+v,0),npoints),ark)
          !enddo
          !
        endif 
        !
        ! Find the irho point and interpolate the rho value
        !
        !
        chi0(1:trove%Nmodes_e) = trove%chi_ref(1:trove%Nmodes_e,irho)
        !
        do n1 = 1,trove%Natoms
           do ix = 1,3
             !
             r_na_t(n1,ix) = trove%b0(n1,ix,irho) + sum( ( trove%Amatrho(n1,ix,1:trove%Nmodes_e,irho) )* &
                                                       ( chi(1:trove%Nmodes_e)-chi0(1:trove%Nmodes_e)) )
             !
           enddo
        enddo
        !
        call FLfromcartesian2local(r_na_t,r_t)
        !
        ! The potential energy function value in terms of the GDC "r"
        pot_t = MLpotenfunc(r_t,r_na_t)
        !
        do imode = 1,trove%Nmodes_e
           !
           g(imode) = MLcoord_direct(chi(imode),1,imode)-MLcoord_direct(chi0(imode),1,imode)
           !
        enddo 
        !
        gvib_t = 0
        !
        do k1 = 1,trove%Nmodes
           do k2 = 1,trove%Nmodes
              !
              fl => trove%g_vib(k1,k2)
              !
              do iterm = 1,fl%Ncoeff
                 !
                 k(:) = FLIndexQ(:,iterm)
                 !
                 f(:) = g(:)**k(:)
                 !
                 gvib_t = gvib_t + fl%field(iterm,irho)*p(k1)*product(f(:))*p(k2)
                 !
              enddo
              ! 
           enddo
        enddo

        !
        grot_t = 0
        !
        do k1 = 1,3
           do k2 = 1,3
              !
              fl => trove%g_rot(k1,k2)
              !
              do iterm = 1,fl%Ncoeff
                 !
                 k(:) = FLIndexQ(:,iterm)
                 !
                 f(:) = g(:)**k(:)
                 !
                 grot_t = grot_t + fl%field(iterm,irho)*Jrot(k1)*product(f(:))*Jrot(k2)
                 !
              enddo
              ! 
           enddo
        enddo

        gcor_t = 0
        !
        do k1 = 1,Nmodes
           do k2 = 1,3
              !
              fl => trove%g_cor(k1,k2)
              !
              do iterm = 1,fl%Ncoeff
                 !
                 k(:) = FLIndexQ(:,iterm)
                 !
                 f(:) = g(:)**k(:)
                 !
                 gcor_t = gcor_t + 2.0d0*fl%field(iterm,irho)*p(k1)*product(f(:))*Jrot(k2)
                 !
              enddo
              ! 
           enddo
        enddo
        !
        H = pot_t-0.5_rk*(gvib_t-grot_t-gcor_t)
        !
        if (present(r)) r = r_t
        if (present(r_na)) r_na = r_na_t
        !
    end function classic_hamilt
    !
  end subroutine FL_rotation_energy_surface


  function choose(n, k) result(res) 
    !
    implicit none
    integer(ik), intent (in) :: n
    integer(ik), intent (in) :: k
    integer(ik) :: res
    integer(ik) :: i, k_  
    
    res = 1
    k_ = k 
    if(k > n - k) then
      k_ = n - k
    endif
    
    do i = 0, k_ - 1
      res = res*(n - i)
      res = res/(i + 1)
    enddo
  end function choose
   
  function sum_choose(num_pos, high_order, pos_val) result(res)
    implicit none 
    integer(ik), intent(in) :: num_pos 
    integer(ik), intent(in) :: high_order
    integer(ik), intent(in) :: pos_val
    integer(ik) :: res  
    integer (ik):: i
       
    res = 0
    if(pos_val == 0) then
      res = res + choose(high_order + num_pos - 1, num_pos - 1)
    else 
      do i = 0, pos_val - 1 
        res = res + choose(high_order + num_pos - 1 - i, num_pos - 1) 
      end do
    endif 

  end function sum_choose
  !
  function powers_from_index(Nmodes, Nindex) result(powers) 
     implicit none 
     integer(ik), intent(in) :: Nmodes, Nindex
     integer(ik) :: powers(Nmodes)   
     integer(ik) :: polyad, max_index, diff, j, term_val, remaining_order 
     !
     polyad = -1
     max_index = 0
     powers = 0
     !
     do while(max_index < Nindex)
       polyad = polyad + 1
       max_index = choose(Nmodes + polyad, Nmodes)
     end do
     if(polyad == 0) then
       diff = Nindex - 1
     else 
       diff = Nindex - choose(Nmodes + polyad - 1, Nmodes)
     endif 
     remaining_order = polyad 
     do j = 1, Nmodes - 1 
       if(remaining_order == 0) then
         powers(j) = 0
         continue 
       else 
           term_val = 1
           do while(diff - sum_choose(Nmodes - j, remaining_order, term_val) > 0)
             term_val = term_val + 1
           end do
           term_val = term_val - 1 
           powers(j) = term_val
           if( term_val > 0) then 
             diff  = diff - sum_choose(Nmodes - j, remaining_order, term_val)
           endif
           remaining_order = remaining_order - term_val
       endif 
     end do 
     powers(Nmodes) = remaining_order
     !
  end function powers_from_index 


  end module fields


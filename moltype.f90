module moltype

  use accuracy
  use lapack

  implicit none
  private
  
  public MoleculeT,molec,MLinitialize_molec,MLequilibrium_chi,MLorienting_a0,&
         MLZmatrixT,MLfromlocal2cartesian,MLdiag_ulen,ML_check_steps,three_j,&
         intensity,MLIntensityT,MLthresholdsT,extF,MLext_locexp,MLvector_product,ML_sym_rotat,ML_euler_rotait,MLdiag_ulen_ark,&
         aacos,MLlinurark,MLlinur,faclog,aasin
  public MLtemplate_poten,MLtemplate_potential,MLtemplate_coord_transform,MLtemplate_b0,MLtemplate_extF,MLtemplate_kinetic
  public MLtemplate_symmetry_transformation,MLtemplate_rotsymmetry,ML_rjacobi_fit_ark,ML_splint,ML_splint_quint,ML_spline
  public MLorienting_a0_across_dadrho,manifold
         !
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level

  type  MLZmatrixT
     character(len=cl)    :: name         ! Identifying name of atom (no effect on anything)
     integer(ik)          :: connect(4)   ! z-matrix connections
  end type MLZmatrixT

  !
  character(len=cl),parameter :: axis_system = 'Eckart'

  abstract interface
    real(ark) function MLtemplate_poten (local,xyz)
      use accuracy
      real(ark),intent(in) ::  local(:)
      real(ark),intent(in) ::  xyz(:,:)
    end function MLtemplate_poten
    !
    real(ark) function MLtemplate_potential (ncoords,natoms,local,xyz,force)
      use accuracy
      integer(ik),intent(in) ::  ncoords,natoms
      real(ark),intent(in)   ::  local(ncoords)
      real(ark),intent(in)   ::  xyz(natoms,3)
      real(ark),intent(in)   ::  force(:)
    end function MLtemplate_potential
    !
    subroutine MLtemplate_extF(rank,ncoords,natoms,local,xyz,mu_xyz)
      use accuracy
      !
      integer(ik),intent(in) ::  rank,ncoords,natoms
      real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
      real(ark),intent(out)  ::  mu_xyz(rank)
      !
    end subroutine MLtemplate_extF
    !
    subroutine MLtemplate_kinetic(nmodes,nterms,rho,g_vib,g_rot,g_cor,pseudo)
      use accuracy
      !
      integer(ik),intent(in) ::  nmodes,nterms
      real(ark),intent(in)   ::  rho
      real(ark),intent(out)  ::  g_vib(nmodes,nmodes,nterms),g_rot(3,3,nterms),g_cor(nmodes,3,nterms),pseudo(Nterms)
      !
    end subroutine MLtemplate_kinetic
    !
    subroutine MLtemplate_symmetry_transformation(ioper,natoms,src,dst)
      use accuracy
      !
      integer(ik),intent(in)   :: ioper  ! group operation  
      integer(ik),intent(in)   :: natoms
      real(ark),intent(in)     :: src(natoms)
      real(ark),intent(out)    :: dst(natoms)
      !
    end subroutine MLtemplate_symmetry_transformation
    !
    subroutine MLtemplate_rotsymmetry(J,K,tau,gamma,ideg)
      use accuracy
      !
      integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
      integer(ik),intent(out) :: gamma,ideg
      !
    end subroutine MLtemplate_rotsymmetry
    !
  end interface 

  abstract interface
    subroutine MLtemplate_b0 (Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
      use accuracy
      !
      integer(ik), intent(in)   :: Npoints
      integer(ik), intent(in)   :: Natoms
      real(ark),   intent(out)  :: b0(Natoms,3,0:Npoints)
      real(ark),   intent(in) ,optional  :: rho_i(0:Npoints)
      real(ark),   intent(out),optional :: rho_ref
      real(ark),   intent(in) ,optional :: rho_borders(2)  ! rhomim, rhomax - borders
      !
    end subroutine MLtemplate_b0
    !
    function MLtemplate_coord_transform (src,ndst,direct) result (dst)
      use accuracy
      integer(ik),intent(in) :: ndst
      real(ark),intent(in)   :: src(:)
      logical,intent(in)     :: direct
      !
      real(ark),dimension(ndst) :: dst
      !
    end function MLtemplate_coord_transform
  end interface 




!procedure(interface1), pointer, nopass :: p
!Edit: In response to your comment, if you want to use the pass keyword, the interface would have to be changed as such:

!ABSTRACT INTERFACE 
!    integer function interface1(passed_object, a)
!        import :: type1
!        class(type1), intent(...) :: passed_object
!        real,         intent(in)  :: a
!    END function interface1
!END INTERFACE
  !
  type basic_function
      procedure(calc_func), pointer, nopass  :: func_pointer
      real(ark) :: coeff
      real(ark) :: inner_expon
      real(ark) :: outer_expon
  end type 
  !
  type ragged_array_lvl_1
      type(basic_function), allocatable :: func_set(:)
      integer :: num_terms
  end type ragged_array_lvl_1
  !
  type ragged_array_lvl_2
      type(ragged_array_lvl_1), allocatable :: mode_set(:)
  end type ragged_array_lvl_2
  !
  abstract interface
      subroutine calc_func(x, y)
        use accuracy
        implicit none
        real(ark) , intent(in) :: x 
        real(ark) , intent(inout) :: y
      end subroutine 
  end interface
  !
  type MoleculeT
     !
     character(len=cl) :: moltype   ! Identifying the Molecule type (e.g. XY3)
     real(ark),pointer :: AtomMasses(:)
     real(ark),pointer  :: req(:)
     real(ark),pointer  :: alphaeq(:)
     real(ark),pointer  :: taueq(:)
     real(ark),pointer  :: local_eq(:)
     real(ark),pointer  :: chi_eq(:)
     real(ark),pointer  :: specparam(:)
     real(ark)          :: rho_ref  
     real(ark)          :: rho_border(2)
     type(MLZmatrixT),pointer  :: zmatrix(:)       ! 
     !
     character(len=cl)         :: coords_transform ! type of the coordinate transformation 
     integer(ik),pointer       :: dihedtype(:)   ! dihedral type 
     character(len=cl),pointer :: coordinates(:,:)  ! Identifying the coordinate system, e.g. 'Cartesian', 'Bond-Angle'
     integer(ik)               :: Nmodes,Ndihedrals,Natoms,Nbonds,Nangles,Ncoords
     integer(ik)               :: parmax        ! number of pot. parameters 
     integer(ik)               :: N_meppars     ! number of MEP. parameters 
     character(len=cl)         :: potentype     ! type of potential function 
     character(len=cl)         :: kinetic_type  ! type of the kinetic function 
     character(len=cl)         :: meptype       ! type of MEP function
     real(ark),pointer         :: force(:)      ! force field
     real(ark),pointer         :: mep_params(:) ! MEP parameters
     integer(ik),pointer       :: mep_ind(:,:)  ! MEP powers-index
     character(len=16),pointer :: forcename(:)  ! characters to store the force constant labels
     integer(ik),pointer       :: pot_ind(:,:)  ! indexes with the powers for every expansion term
     integer(ik),pointer       :: ifit(:)       ! varying indexes to control fitting
     character(len=cl)         :: symmetry      ! molecular symmetry
     character(len=cl)         :: IO_primitive  ! control writing/reading of primitive basis sets 
     character(len=cl)         :: chk_primitive_fname  ! filename to store primitive functions on a grid
     !
     !procedure(MLtemplate_poten),pointer :: potenfunc => null ()
     !
     type(ragged_array_lvl_2), allocatable :: basic_function_list(:)
     logical  :: basic_function_set = .false.
  end type MoleculeT
  !
  type xyT
      real(ark) :: val
  end type xyT
  !
  type(MoleculeT),save :: molec
  !

  !external field (function)expansion in terms of local coordinates


  type MLext_locexp
     integer(ik)            :: rank    ! rank of the external field (number of components), 
                                       ! e.g. 3 for dipole, 9 for polarizability
     integer(ik), pointer   :: maxord(:)  ! maximal order of the xi-expansions for each component.
     integer(ik), pointer   :: nterms(:)
     integer(ik), pointer   :: term(:, :, :)
     integer(ik), pointer   :: ifit(:, :)
     real(ark), pointer     :: coef(:, :)
     character(cl), pointer :: name(:, :)
     real(ark), pointer     :: fdstep(:)
     character(cl), pointer :: intcoords(:)
     real(ark), pointer     :: geom_ref(:)
     integer(ik)            :: irho_ref
     character(cl)          :: ftype = 'GENERAL'  ! field type 
     real(rk)               :: matelem_threshold = -1e0   ! threshold to set the primitive matrix elements to zero, 
                                                          !required to reduce numerical noice in overtone intensities
  end type MLext_locexp

  !
  ! external field, e.g. dipole moment, palarizability etc.
  !
  type(MLext_locexp)         :: extF


  !type MLextf_xiexp
  !   integer(ik)          :: maxord(3)
  !   real(ark), pointer   :: fdstep(:)
  !   real(ark)            :: coef_thr
  ! 
  !   integer(ik)          :: nterms(3)
  !   integer(ik), pointer :: term(:, :, :)
  !   real(ark), pointer   :: coef(:, :, :)
  !   logical              :: man_terms = .false.
  !   logical              :: man_coefs = .false.
  !
  !   integer(ik), pointer :: npoints(:)
  !   real(ark), pointer   :: b(:, :)
  !end type MLextf_xiexp
  !type(MLextf_xiexp)        :: extF_xi


  !vars for transition intensities calculations


  type MLthresholdsT
     real(rk) :: intensity    = -1e0    ! threshold defining the output intensities
     real(rk) :: linestrength = -1e0    ! threshold defining the output linestrength
     real(rk) :: coeff        = -1e0    ! threshold defining the eigenfunction coefficients
                                        ! taken into account in the matrix elements evaluation.
    real(rk) :: leading_coeff = 0.1
  end type MLthresholdsT

  type MLIntensityT
     logical             :: do = .false.     ! process (.true.) or not (.false.) the intensity (or TM) calculations 
     character(cl)       :: action           ! type of the intensity calculations:
                                             ! absorption, emission, tm (transition moments),
                                             !  raman, and so on. 
     real(rk)            :: temperature      ! temperature in K
     real(rk)            :: part_func        ! partition function 
     real(rk)            :: ZPE              ! zero point energy
     type(MLthresholdsT) :: threshold        ! different thresholds
     real(rk),pointer    :: gns(:)           ! nuclear stat. weights
     integer(ik),pointer :: isym_pairs(:)    ! numbers defining symmetry pairs with allowed transitions, analogous to gns
     real(rk)            :: freq_window(1:2) ! frequency window (1/cm)
     real(rk)            :: erange_low(1:2)  ! energy range for the lower state
     real(rk)            :: erange_upp(1:2)  ! energy range for the upper state
     integer(ik)         :: J(1:2)           ! range of J-values, from..to; in order to avoid double counting of transitions
                                             ! in the calculations it is always assumed that 
                                             ! J1<=J_lower<=J2 and J1<=J_upper<J2;
                                             !
     integer(ik),pointer :: v_low(:,:)       ! lower state range of the quantun numbers employed 
     integer(ik),pointer :: v_upp(:,:)       ! upper state range of the quantun numbers employed 
                                             ! in intensity calculations; (imode,1:2), 
                                             ! where 1 stands for the beginning and 2 for the end. 
     !
     integer(ik)         :: swap_size    = 0 ! the number of vectors to keep in memory
     character(cl)       :: swap = "NONE"    ! whether save the compacted vectors or read
     character(cl)       :: swap_file  ="compress"   ! where swap the compacted eigenvectors to
     integer(ik)         :: int_increm = 1e9 ! used to print out the lower energies needed to select int_increm intensities
     integer(ik)         :: Ncache = 10000 ! used to cache intensities before prinout  to speed up 
     real(rk)         :: factor = 1.0d0   ! factor <1 to be applied the maxsize of the vector adn thus to be shrunk 
     real(rk)         :: wallclock    ! wallclock limit, needed to estmate how many transitions can be processed within one job
     logical          :: reduced      ! process intensity in a reduced symmetry adapted approach, only the (1,2) degenerate component
     logical          :: pruning = .false.    ! for the TM-based basis set pruning compute and store the max vib. intensity for each state
     logical          :: output_short = .false.    ! Long output is with all quantum numbers and energies; short is with indeces, energies and A-coef-s only
     logical          :: tdm_replace = .false.     ! Replace vibrational trandipole moments with experimental values
     !
     ! variables used in extfield module
     integer(ik)   :: tens_rank = 1
     integer(ik)   :: tens_ialpha = 3
     integer(ik)   :: tens_ibeta = 3
     character(cl) :: tens_oper = 'NA'
     integer(ik)   :: tens_oper_ielem = 1
     character(cl) :: tens_frame = 'NA'
     integer(ik)   :: tens_ibatch = 1
     real(rk)      :: tens_omega
     !
     integer(ik) :: nvib_quanta_upp = 0
     integer(ik) :: nvib_quanta_low = 0
     integer(ik) :: vib_quanta_upp(20,1000)
     integer(ik) :: vib_quanta_low(20,1000)
     !
 end type MLIntensityT

 type(MLIntensityT),save :: intensity
 !
 integer(ik) :: manifold

  contains



!
! The procedure defines  for a given molecule 
! the equilibrium geometry, potential function, masses, and something else.
! In fact it redirects the execution to the case related subroutine  
!
  subroutine MLinitialize_molec(Moltype,Coordinates,coords_transform,&
                                  Nbonds,Nangles,Ndihedrals,dihedtype_,&
                                  AtomMasses,local_eq, &
                                  force_,forcename_,ifit_,pot_ind_,specparam,potentype,kinetic_type,&
                                  IO_primitive,chk_numerov_fname,&
                                  symmetry_,rho_border,zmatrix_)


  character(len=cl),intent(in)  :: Moltype
  character(len=cl),intent(in)  :: Coordinates(:,:)  
  character(len=cl),intent(in)  :: coords_transform    ! coordinate transformation type
  integer(ik),   intent(in)  :: Nbonds,Nangles,Ndihedrals
  integer(ik),   intent(in)  :: dihedtype_(0:Ndihedrals)
  real(ark),   intent(in)  :: AtomMasses(:)
  real(ark),   intent(in)  :: local_eq(:)
  !
  real(ark),intent(in)         :: force_(:)
  integer(ik),intent(in)       :: ifit_(:)        ! indexes contrilloing fit
  integer(ik),intent(in)       :: pot_ind_(:,:)   ! indexes with the powers for every expansion term
  character(len=16),intent(in) :: forcename_(:)
  real(ark),   intent(in)      :: specparam(:)
  !
  character(len=cl),intent(in)  :: potentype,kinetic_type
  character(len=cl),intent(in)  :: symmetry_
  character(len=cl),intent(in)  :: IO_primitive,chk_numerov_fname
  real(ark)                     :: rho_border(2)     ! rhomim, rhomax - borders
  type(MLZmatrixT),intent(in)   :: zmatrix_(:)       ! 
  !
  integer(ik)              :: alloc,Ncoords
    !
    if (verbose>=4) write(out,"(/'MLinitialize_molec/start: molecular and potential parameters')") 
    !
    molec%moltype = Moltype
    molec%coords_transform = coords_transform
    !
    molec%Natoms = size(AtomMasses)
    molec%Ndihedrals = Ndihedrals
    molec%Nmodes = size(Coordinates,dim=2)
    molec%Nbonds = Nbonds
    molec%Nangles = Nangles
    !
    Ncoords = Nbonds+Nangles+Ndihedrals
    molec%Ncoords = Ncoords
    !
    molec%parmax = size(force_)
    !
    allocate (molec%force(molec%parmax),molec%forcename(molec%parmax),&
              molec%ifit(molec%parmax),molec%pot_ind(1:Ncoords,molec%parmax),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i8,' allocating matix force ')") alloc
        stop 'MLforce_and_geometry - cannot allocate force'
    end if
    !
    allocate (molec%coordinates(3,molec%Nmodes),molec%specparam(Ncoords),molec%AtomMasses(molec%Natoms),molec%req(Nbonds),&
              molec%alphaeq(Nangles),molec%taueq(Ndihedrals),molec%local_eq(Ncoords),&
              molec%dihedtype(0:Ndihedrals),stat=alloc)
    if (alloc/=0) then
       write (out,"(' Error ',i8,' allocating matixes for masses, reqs, or alphaeq-s ')") alloc
       stop 'MLforce_and_geometry - cannot allocate masses, reqs, or alphaeq-s'
    end if
    !
    allocate (molec%zmatrix(molec%Natoms),stat=alloc)
    if (alloc/=0) then
        write (out,"(' Error ',i8,' allocating matix zmatrix ')") alloc
        stop 'MLforce_and_geometry - cannot allocate zmatrix'
    end if
    !
    molec%zmatrix(1:molec%Natoms) = zmatrix_(1:molec%Natoms)
    !
    molec%force = force_
    molec%forcename = forcename_
    molec%pot_ind = pot_ind_
    molec%ifit = ifit_
    !
    ! store equilibrium parameters 
    !
    molec%req(1:Nbonds) = local_eq(1:Nbonds)
    molec%alphaeq(1:Nangles) = local_eq(Nbonds+1:Nbonds+Nangles)
    molec%taueq(1:Ndihedrals) = local_eq(Nbonds+Nangles+1:Ncoords)
    !
    molec%local_eq = local_eq
    molec%potentype = potentype
    molec%kinetic_type = kinetic_type
    molec%atomMasses = AtomMasses
    molec%specparam = specparam
    molec%dihedtype = dihedtype_
    molec%symmetry = symmetry_
    molec%IO_primitive = IO_primitive
    molec%chk_primitive_fname = chk_numerov_fname
    !
    !molec%potenfunc => MLpoten_xy2_morbid
    !
    ! Store the molecular structure parameters 
    !
    molec%coordinates(1,:) = Coordinates(1,:)
    molec%coordinates(2,:) = Coordinates(2,:)
    molec%coordinates(3,:) = Coordinates(3,:)
    !
    ! define the integration borders for each mode
    !
    molec%rho_border(:) = rho_border(:)
    !
    !call xy2_initialize(xy2)
    !
    if (verbose>=4) write(out,"('MLinitialize_molec/end')") 

  end subroutine MLinitialize_molec



!
! Here we define chi_eq - 
! equilibrium values of the internal coordinates chi
!
  subroutine MLequilibrium_chi(chi_eq)

     real(ark),intent(in)     ::  chi_eq(:)    ! Equilibrium values for the internal coordinates chi 
     integer(ik)  :: alloc


    if (verbose>=5) write(out,"(/'MLequilibrium_chi/start: chi_eq')") 

    allocate (molec%chi_eq(size(chi_eq)),stat=alloc)
    !
    molec%chi_eq = chi_eq
    !
    if (verbose>=5) write(out,"('MLequilibrium_chi/end')") 


  end subroutine MLequilibrium_chi


 subroutine MLorienting_a0(Natoms,AtomMasses,a0,transform,method)

     integer(ik),intent(in)  :: Natoms
     real(ark),   intent(in)  :: AtomMasses(Natoms)
     real(ark),   intent(inout) :: a0(Natoms,3)
     real(ark),   intent(out),optional :: transform(3,3)
     character(cl),optional :: method
     !
     real(ark)                :: a0_t(Natoms,3)
     real(ark)                :: CM_shift
     integer(ik)             :: ix,jx,imx
     !
     real(ark)                ::  b(3,1),a(3,3),c(3,3),a_t,c_t(3),Inert(3)
     real(rk)                 ::  a_rk(Natoms,3),b_rk(3,1)


      if (verbose>=5) write(out,"(/'MLorienting_a0/start')") 

      !
      ! Find center of mass
      !
      do ix = 1,3 
        CM_shift = sum(a0(:,ix)*AtomMasses(:))/sum(AtomMasses(:))
        a0(:,ix) = a0(:,ix) - CM_shift
      enddo 
      !
      Inert(1) = sum(AtomMasses(:)*( a0(:,2)**2+ a0(:,3)**2) )
      Inert(2) = sum(AtomMasses(:)*( a0(:,1)**2+ a0(:,3)**2) )
      Inert(3) = sum(AtomMasses(:)*( a0(:,1)**2+ a0(:,2)**2) )
      !
      ! First Eckart equation
      ! 
      do ix = 1,3 
         !
         a_t = sum(AtomMasses(:)*a0(:,ix))
         if (abs(a_t)>10.0_rk**(-rk)) then 
              write(out,"('MLorienting_a0: a0 is not a solution of Eckart 1 for  ix =',i4,d18.8)") ix,a_t
                 stop 'MLorienting_a0: a0 is not solution of Eckart1'
         endif
         !
      enddo


      !
      ! Second Eckart equation
      ! 
      a = 0 
      do ix = 1,3 
         !
         a(ix,ix) = Inert(ix)
         !
         do jx = 1,3
            if (ix/=jx) a(ix,jx) =-sum(AtomMasses(:)*a0(:,ix)*a0(:,jx) )
         enddo
         !
      enddo

      !
      !a = real(a,kind=rk)
      !
      if (present(method)) then 
        !
        select case (trim(method)) 
          !
        case('LAPACK')
          !
          a_rk = a
          !
          call lapack_syev(a_rk,b_rk(:,1))
          !
          c = a_rk
          !
        case('ULEN')
          !
          call MLdiag_ulen_ark(3,a,b(:,1),c)
          !
        case default
          write (out,"('MLorienting_a0: method  ',a,' unknown')") trim(method)
          stop 'MLorienting_a0: method unknown'
          !
        end select  
          !
      else 
          !
          call MLdiag_ulen_ark(3,a,b(:,1),c)
          !
      endif 
      !
      !
      !
      !c = a
      !
      !do ix=1,2
      !  do jx=ix,3
      !    !
      !    if (abs(a(ix,ix))<abs(a(jx,ix))) then 
      !        !
      !        a(:,ix) = c(:,jx)
      !        a(:,jx) = c(:,ix)
      !        !
      !    endif 
      !    !
      !  enddo
      !enddo
      !
      !c_t(:) = a(:,1)
      !a(:,1) = a(:,3)
      !a(:,3) = c_t(:)
      !!
      !do ix=1,3
      !  !
      !  if (a(ix,ix)<0.0_rk) a(:,ix) = -a(:,ix)
      !  !
      !enddo
      !
      ! Found coordinate transformation "c"
      !
      !c = real(a,kind=rk)
      !
      ! Transformation of a0 
      !
      do ix = 1,Natoms
         a0_t(ix,:) = matmul(transpose(c),a0(ix,:))
      enddo
      !
      Inert(1) = sum(AtomMasses(:)*( a0_t(:,2)**2+ a0_t(:,3)**2) )
      Inert(2) = sum(AtomMasses(:)*( a0_t(:,1)**2+ a0_t(:,3)**2) )
      Inert(3) = sum(AtomMasses(:)*( a0_t(:,1)**2+ a0_t(:,2)**2) )
      !
      imx = minloc(Inert(:),dim=1)
      !
      !c_t(:) = c(imx,:) ; c(imx,:) = c(3,:) ; c(3,:) = c_t(:)
      !c_t(:) = c(imn,:) ; c(imn,:) = c(1,:) ; c(1,:) = c_t(:)
      !
      !c_t(:) = c(:,imx) ; c(:,imx) = c(:,3) ; c(:,3) = c_t(:)
      !
      !do ix = 1,Natoms
      !   a0_t(ix,:) = matmul(transpose(c),a0(ix,:))
      !enddo
      !
      !Inert(1) = sum(AtomMasses(:)*( a0_t(:,2)**2+ a0_t(:,3)**2) )
      !Inert(2) = sum(AtomMasses(:)*( a0_t(:,1)**2+ a0_t(:,3)**2) )
      !Inert(3) = sum(AtomMasses(:)*( a0_t(:,1)**2+ a0_t(:,2)**2) )
      !
      !imn = minloc(Inert(:),dim=1)
      !c_t(:) = c(:,imn) ; c(:,imn) = c(:,1) ; c(:,1) = c_t(:)
      !
      do ix = 1,Natoms
         a0(ix,:) = matmul(transpose(c),a0(ix,:))
      enddo
      !
      Inert(1) = sum(AtomMasses(:)*( a0(:,2)**2+ a0(:,3)**2) )
      Inert(2) = sum(AtomMasses(:)*( a0(:,1)**2+ a0(:,3)**2) )
      Inert(3) = sum(AtomMasses(:)*( a0(:,1)**2+ a0(:,2)**2) )
      !
      if (verbose>=6) then 
        write(out,"('i0',3f12.8)") Inert(1:3)
      endif
      !
      ! Second Eckart equation
      ! 
      do ix = 1,3 
         do jx = 1,3 
            if (ix/=jx) then  
               a_t =  sum(AtomMasses(:)*a0(:,ix)*a0(:,jx) )

               if (abs(a_t)>1000.0_rk*small_) then 
                   write(out,"('MLorienting_a0: a0 is not a solution of Eckart 2 for ix,jx =',2i4,d18.8)") ix,jx,a_t
                   stop 'MLorienting_a0: a0 is not solution of Eckart2'
               endif
            endif
         enddo
         !
      enddo
      !
      ! we might need the transformation matrix out there
      !
      if (present(transform)) then 
         transform = c
      endif 
      !
      if (verbose>=5) write(out,"('MLorienting_a0/stop')") 


  end subroutine MLorienting_a0
  !


 subroutine MLorienting_a0_across_dadrho(Natoms,npoints,AtomMasses,arho_borders,a0,da0,periodic,lincoord,method,transform)

     integer(ik),intent(in)   :: Natoms,npoints
     real(ark),   intent(in)  :: AtomMasses(Natoms)
     real(ark),   intent(in)  :: arho_borders(2)  ! rhomim, rhomax - borders
     real(ark),   intent(inout) :: a0(Natoms,3,0:npoints),da0(Natoms,3,0:npoints)
     logical     ,intent(in)  :: periodic
     real(ark),   intent(out),optional :: transform(3,3)
     character(cl),optional :: method
     integer(ik),optional :: lincoord
     !
     real(ark) :: CM_shift,Eq(3),A(3,3),A1(3,3),A2(3,3),C1(3,3),C2(3,3),b1(3),b2(3),S(3,3),T(3,3)
     real(ark) :: b0(Natoms,3,0:npoints),db0(Natoms,3,0:npoints)
     real(rk)  :: a_rk(3,3),b_rk(3,1)
     integer(ik)              :: ix,jx,jy,imx,ipoint,i,ieq,jatom,ilincoord=0
      !
      if (verbose>=5) write(out,"(/'MLorienting_a0_across_dadrho/start')") 
      !
      if (present(lincoord)) ilincoord = lincoord
      !
      ! Find center of mass
      !
      do ipoint = 0,npoints
        !
        do ix = 1,3 
          CM_shift = sum(a0(:,ix,ipoint)*AtomMasses(:))/sum(AtomMasses(:))
          a0(:,ix,ipoint) = a0(:,ix,ipoint) - CM_shift
        enddo 
        !
      enddo
      !
      do ipoint = 0,npoints
         !
         ieq = 0
         do ix = 1,3
            ! 
            if (ix/=ilincoord) then 
               !
               ieq = ieq+1
               !
               Eq(ieq) = 0
               do jatom = 1,Natoms  
                  do jx = 1,3
                     do jy = 1,3
                       !
                       Eq(ieq) = Eq(ieq) + AtomMasses(jatom)*epsil(ix,jx,jy)*a0(jatom,jx,ipoint)*da0(jatom,jy,ipoint)
                       !
                     enddo
                  enddo
               enddo
               !
            endif
            !
         enddo
         !
         do jx = 1,3
            do jy = 1,3
              !
              A(jx,jy) = sum(AtomMasses(:)*a0(:,jx,ipoint)*da0(:,jy,ipoint))
              !
            enddo
         enddo
         !
         A1 = matmul(A,transpose(A))
         A2 = matmul(transpose(A),A)
         !
         call MLdiag_ulen_ark(3,A1,b1,C1)
         call MLdiag_ulen_ark(3,A2,b2,C2)
         !
         a_rk = A1
         call lapack_syev(a_rk,b_rk(:,1))
         C1 = a_rk
         !
         a_rk = A2
         call lapack_syev(a_rk,b_rk(:,1))
         C2 = a_rk
         !
         T = matmul(C2,transpose(C1))
         !
         A1 = matmul(S,transpose(S))
         A2 = matmul(transpose(S),S)
         !
         do jatom = 1,Natoms  
            !
            b0 (jatom,:,ipoint) = matmul(T, a0(jatom,:,ipoint))
            db0(jatom,:,ipoint) = matmul(T,da0(jatom,:,ipoint))
            !
         enddo
         !
         ieq = 0
         Eq = 0 
         do ix = 1,3
            ! 
            if (ix/=ilincoord) then 
               !
               ieq = ieq+1
               !
               Eq(ieq) = 0
               do jatom = 1,Natoms  
                  do jx = 1,3
                     do jy = 1,3
                       !
                       Eq(ieq) = Eq(ieq) + AtomMasses(jatom)*epsil(ix,jx,jy)*b0(jatom,jx,ipoint)*db0(jatom,jy,ipoint)
                       !
                     enddo
                  enddo
               enddo
               !
            endif
            !
         enddo
         !
         do jatom = 1,Natoms  
            !
            b0 (jatom,:,ipoint) = matmul(a0 (jatom,:,ipoint),T)
            db0(jatom,:,ipoint) = matmul(da0(jatom,:,ipoint),T)
            !
         enddo
         !
         ieq = 0
         Eq = 0 
         do ix = 1,3
            ! 
            if (ix/=ilincoord) then 
               !
               ieq = ieq+1
               !
               Eq(ieq) = 0
               do jatom = 1,Natoms  
                  do jx = 1,3
                     do jy = 1,3
                       !
                       Eq(ieq) = Eq(ieq) + AtomMasses(jatom)*epsil(ix,jx,jy)*b0(jatom,jx,ipoint)*db0(jatom,jy,ipoint)
                       !
                     enddo
                  enddo
               enddo
               !
            endif
            !
         enddo
         !
         S = matmul(T,A)
         !
         T = matmul(C1,transpose(C2))
         !
         S = matmul(A,T)
         !
         continue
         !
      enddo
      !
      if (present(transform)) then 
         !transform = c
         transform = 0
      endif 
      !
      if (verbose>=5) write(out,"('MLorienting_a0_across_dadrho/stop')") 


  end subroutine MLorienting_a0_across_dadrho



   !
   ! 6-dim rotations: (123),(312),(12)*,(23)*,(13)*
   !
   function ML_sym_rotat(ioper,a) result (f)

      integer(ik)   :: ioper
      real(ark),intent(in)   :: a(3,3)
      real(ark),dimension(3,3) :: f

      select case(ioper)
         !
      case (1)
         !
         f = a
         !
      case (2)
         !
         f = cshift(a,1,1)
         !
      case (3)
         !
         f = cshift(a,2,1)
         !
      case (4)
         !
         f(1,:) = a(1,:) ; f(2,:) = a(3,:) ; f(3,:) = a(2,:) ; f = transpose(f)
         !
      case (5)
         !
         f(1,:) = a(3,:) ; f(2,:) = a(2,:) ; f(3,:) = a(1,:) ; f = transpose(f)  
         !
      case (6)
         !
         f(1,:) = a(2,:) ; f(2,:) = a(1,:) ; f(3,:) = a(3,:) ; f = transpose(f)  
         !
      end select

   end function ML_sym_rotat


 subroutine ML_check_steps(xi,itype,imode,fstep)

   real(ark),intent(in) :: xi(1:molec%Nmodes)
   integer(ik),intent(in) :: itype,imode
   real(ark),intent(out) :: fstep(2)
   real(ark) :: rhoe

   if (verbose>=6) write(out,"(/'ML_check_steps/start')") 
   !
   fstep(1:2) = 1.0_ark
   !
   select case(trim(molec%coordinates(itype,imode)))
   case default
        write (out,"('ML_check_steps: coordinate type ',a,' unknown')") trim(molec%coordinates(itype,imode))
        stop 'ML_check_steps - bad coordinate-type'
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
        !
     case('COSX','SINX') 
        !
        if (xi(imode)-1.0_ark < small_) then 
            !
            fstep(2) = 0.0_ark
            !
        endif 
        if (xi(imode)+1.0_ark < small_) then 
            !
            fstep(1) = 0.0_ark
            !
        endif 
        !
     case('1-COSX') 
        !
        if (abs( xi(imode) )<small_) then 
            !
            fstep(1) = 0.0_ark
            !
        endif 
        if (abs( xi(imode) )>small_) then 
            !
            fstep(2) = 0.0_ark
            !
        endif 
        !
    end select 
    !
    if (fstep(1)+fstep(2)==0.0_ark) then 
       write (out,"('ML_check_steps: no numerical derivatives allowed around point ',f18.8)") xi(imode)
       write (out,"('imode -  ',i8)") imode
       stop 'ML_check_steps - bad point for derivatives'
    endif

   
   if (verbose>=6) write(out,"('ML_check_steps/end')") 
   
    
 end subroutine ML_check_steps


   subroutine MLfromlocal2cartesian(pm,r,cartesian)

    integer(ik),intent(in)  :: pm
    real(ark),intent(in)    :: r(molec%Ncoords)
    real(ark),intent(out)   :: cartesian(molec%Natoms,3)

    real(ark)   :: x(molec%Natoms,3)
    real(ark)   :: n1(3),n2(3),n3(3)
    !
    real(ark)   :: rot(3,3),CM_shift
    real(ark)   :: alpha,phi,v12(3),v23(3),f_t,beta,alpha3,cosa3,cosphi
    !
    real(ark)   :: rbond
    !
    integer(ik) ::  iangle,idihedral,iatom,j
    integer(ik) ::  p0,p1,p2,p3,p4,ix,kappa(3),idelta,ikappa,zeta
    character(len=cl)  :: txt = 'MLfromlocal2cartesian'
     !
     if (verbose>=6) write(out,"(/'MLfromlocal2cartesian/start')") 
     !
     if (abs(pm)/=1) then 
       !
       write(out,"(/'MLfromlocal2cartesian: illegal pm, not 1,-1: ',i8)") pm
       stop 'illegal pm'
       !
     endif 
     !
     ! cartesian coordinates
     !
     x = 0
     !
     ! #1 - at the reference coordinate of atom #1
     !
     x(1,:) = 0
     !
     ! #2 - on the vector between refernce atoms 1 and 2
     !
     n1(:) = (/0.0_ark,0.0_ark,1.0_ark/)
     !
     x(2,:) = x(1,:)+n1(:)*r(1)
     !
     iangle = 1
     idihedral = 0
     !
     if (molec%NAtoms>2) then
       !
       iatom = 3
       !
       J = molec%zmatrix(iatom)%connect(4)
       !
       select case (J)
         !
       case (0:101)
         !
         p1 = molec%zmatrix(iatom)%connect(1)
         p2 = molec%zmatrix(iatom)%connect(2)
         p3 = molec%zmatrix(iatom)%connect(3)
         !
         alpha = r(molec%Nbonds+1)
         !
         if (p1==1) then 
           n2(:) = (/sin(alpha),0.0_ark,cos(alpha)/)
         else 
           n2(:) = (/sin(alpha),0.0_ark,-cos(alpha)/)
         endif
         x(3,:) = x(p1,:)+n2(:)*r(2)
         !
       case (102)
         !
         idihedral = idihedral + 1
         !
         zeta = molec%zmatrix(iatom)%connect(3)
         !
         idelta = 0
         kappa(3) = zeta
         do ikappa = 1,3
           if (ikappa==zeta) cycle
           idelta = idelta + 1
           kappa(idelta) = ikappa
         enddo
         !
         n1 = 0
         !
         n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         idihedral = idihedral + 1
         n1(kappa(1)) =-r(molec%Nbonds+molec%Nangles+idihedral)
         !
         n1(kappa(3)) =-sqrt( 1.0_ark-( n1(kappa(1))**2+n1(kappa(2))**2 ) )
         !
         rbond = r(iatom-1)
         !
         p1 = molec%zmatrix(iatom)%connect(1)
         !
         x(iatom,:) = x(p1,:)+rbond*n1(:)
         !
       case (103)
         !
         idihedral = idihedral + 1
         !
         zeta = molec%zmatrix(iatom)%connect(3)
         !
         idelta = 0
         kappa(3) = zeta
         do ikappa = 1,3
           if (ikappa==zeta) cycle
           idelta = idelta + 1
           kappa(idelta) = ikappa
         enddo
         !
         n1 = 0
         !
         n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         idihedral = idihedral + 1
         n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         n1(kappa(3)) = r(iatom-1)
         !
         p1 = molec%zmatrix(iatom)%connect(1)
         !
         x(iatom,:) = x(p1,:)+n1(:)
         !
       case (104)
         !
         idihedral = idihedral + 1
         !
         zeta = molec%zmatrix(iatom)%connect(3)
         !
         idelta = 0
         kappa(3) = zeta
         do ikappa = 1,3
           if (ikappa==zeta) cycle
           idelta = idelta + 1
           kappa(idelta) = ikappa
         enddo
         !
         n1 = 0
         !
         n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         idihedral = idihedral + 1
         n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         rbond = r(iatom-1)
         !
         n1(kappa(3)) =-sqrt(rbond**2-(n1(kappa(1))**2+n1(kappa(2))**2))
         !
         p1 = molec%zmatrix(iatom)%connect(1)
         !
         x(iatom,:) = x(p1,:)+n1(:)
         !
       case (105)
         !
         idihedral = idihedral + 1
         !
         zeta = molec%zmatrix(iatom)%connect(3)
         !
         idelta = 0
         kappa(3) = zeta
         do ikappa = 1,3
           if (ikappa==zeta) cycle
           idelta = idelta + 1
           kappa(idelta) = ikappa
         enddo
         !
         n1 = 0
         !
         n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         idihedral = idihedral + 1
         n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         rbond = r(iatom-1)
         !
         n1(kappa(3)) =-sqrt(1.0_ark-(n1(kappa(1))**2+n1(kappa(2))**2))
         !
         p1 = molec%zmatrix(iatom)%connect(1)
         !
         x(iatom,:) = x(p1,:)+n1(:)*rbond
         !
       case (106)
         !
         idihedral = idihedral + 1
         !
         zeta = molec%zmatrix(iatom)%connect(3)
         !
         idelta = 0
         kappa(3) = zeta
         do ikappa = 1,3
           if (ikappa==zeta) cycle
           idelta = idelta + 1
           kappa(idelta) = ikappa
         enddo
         !
         n1 = 0
         !
         n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         idihedral = idihedral + 1
         n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral)
         !
         rbond = r(iatom-1)
         !
         n1(kappa(3)) =-sqrt(1.0_ark-(n1(kappa(1))**2+n1(kappa(2))**2))
         !
         p1 = molec%zmatrix(iatom)%connect(1)
         !
         x(iatom,:) = x(p1,:)+n1(:)*rbond
         !
       case (107)
         !
         idihedral = idihedral + 1
         !
         zeta = molec%zmatrix(iatom)%connect(3)
         !
         idelta = 0
         kappa(3) = zeta
         do ikappa = 1,3
           if (ikappa==zeta) cycle
           idelta = idelta + 1
           kappa(idelta) = ikappa
         enddo
         !
         n1 = 0
         !
         n1(kappa(1)) =-r(molec%Nbonds+molec%Nangles+idihedral)
         !
         idihedral = idihedral + 1
         n1(kappa(2)) =-r(molec%Nbonds+molec%Nangles+idihedral)
         !
         rbond = r(iatom-1)
         !
         n1(kappa(3)) =-sqrt(rbond**2-(n1(kappa(1))**2+n1(kappa(2))**2))
         !
         p1 = molec%zmatrix(iatom)%connect(1)
         !
         x(iatom,:) = x(p1,:)+n1(:)         
         !
       case (108)
         !
         idihedral = idihedral + 1
         !
         zeta = molec%zmatrix(iatom)%connect(3)
         !
         idelta = 0
         kappa(3) = zeta
         do ikappa = 1,3
           if (ikappa==zeta) cycle
           idelta = idelta + 1
           kappa(idelta) = ikappa
         enddo
         !
         n1 = 0
         !
         n1(kappa(1)) =-r(molec%Nbonds+molec%Nangles+idihedral)
         !
         idihedral = idihedral + 1
         n1(kappa(2)) =-r(molec%Nbonds+molec%Nangles+idihedral)
         !
         n1(kappa(3)) = r(iatom-1)
         !
         p1 = molec%zmatrix(iatom)%connect(1)
         !
         x(iatom,:) = x(p1,:)+n1(:)         
         !
       end select
       !
     endif
     !
     ! loop over the rest atoms
     !
     do iatom = 4,molec%NAtoms
        !
        J = molec%zmatrix(iatom)%connect(4)
        !
        ! it is special when we work with the third atom
        ! it will be distinguished by "J"
        !
        !if (iatom == 3) J = -1
        !
        rbond = r(iatom-1)
        !
        iangle = iangle + 1
        alpha = r(molec%Nbonds+iangle)
        !
        ! connections
        !
        p0 = iatom
        p1 = molec%zmatrix(iatom)%connect(1)
        p2 = molec%zmatrix(iatom)%connect(2)
        p3 = molec%zmatrix(iatom)%connect(3)
        !
        select case (J) 
          !
        case(-1,0)
           !
           iangle = iangle + 1
           !
           beta = r(molec%Nbonds+iangle)
           !
           v12 = x(p2,:)-x(p1,:)
           v23 = x(p3,:)-x(p1,:)
           n2 = v12/sqrt(sum(v12(:)**2)) 
           n3 = MLvector_product(v23,v12)
           !
           f_t = sqrt(sum(n3(:)**2))
           !
           n3 = n3/f_t
           !
           n1 = MLvector_product(n2,n3)
           !
           cosa3 = sum(n2(:)*v23(:))/sqrt(sum(v23(:)**2))
           !
           alpha3 = acos(cosa3)
           !
           alpha3 = aacos(cosa3,txt)
           !
           cosphi =  ( cos(beta)-cos(alpha)*cos(alpha3) )/( sin(alpha)*sin(alpha3) )
           !
           !phi = acos(cosphi)
           !
           phi = aacos(cosphi,txt)
           !
           x(iatom,:) = x(p1,:)+rbond*( cos(alpha)*n2(:) &
                                      + sin(alpha)*cos(phi)*n1(:) &
                                      + sin(alpha)*sin(phi)*n3(:) )
           !
           continue
           !
        case(1)
           !
           idihedral = idihedral + 1
           !
           phi = r(molec%Nbonds+molec%Nangles+idihedral)
           !
           v12 = x(p2,:)-x(p1,:)
           v23 = x(p3,:)-x(p1,:)
           n2 = v12/sqrt(sum(v12(:)**2)) 
           n3 = MLvector_product(v12,v23)
           !
           f_t = sqrt(sum(n3(:)**2))
           !
           n3 = n3/f_t
           !
           n1 = MLvector_product(n2,n3)
           !
           x(iatom,:) = x(p1,:)+rbond*( cos(alpha)*n2(:) &
                                      + sin(alpha)*cos(phi)*n1(:) &
                                      - sin(alpha)*sin(phi)*n3(:) )
           !
        case(-2,2,-202,202,-302,302,-402,402,-502,502,-602,602)
           !
           idihedral = idihedral + 1
           !
           phi = r(molec%Nbonds+molec%Nangles+idihedral)
           !
           v12 = x(p2,:)-x(p1,:)
           v23 = x(p3,:)-x(p2,:)
           !
           !if (J<0) v12 = -v12
           !
           n2 = v12/sqrt(sum(v12(:)**2)) 
           n3 = MLvector_product(v23,v12)
           !
           f_t = sqrt(sum(n3(:)**2))
           !
           n3 = n3/f_t
           !
           n1 = MLvector_product(n2,n3)
           !
           if (J<0) n3 = -n3
           !
           x(iatom,:) = x(p1,:)+rbond*( cos(alpha)*n2(:) &
                                      + sin(alpha)*cos(phi)*n1(:) &
                                      - sin(alpha)*sin(phi)*n3(:) )
        case(4:100)
           !
           !write(out,"(/'MLfromlocal2cartesian: case(4:100)')")
           !stop 'not implemented case(4:100)'
           !
           iangle = iangle + 1
           !
           !p3 = molec%zmatrix(iatom)%connect(4)
           !
           beta = r(molec%Nbonds+iangle)
           !
           v12 = x(p2,:)-x(p1,:)
           v23 = x(p3,:)-x(p1,:)
           n2 = v12/sqrt(sum(v12(:)**2)) 
           n3 = MLvector_product(v12,v23)
           !
           f_t = sqrt(sum(n3(:)**2))
           !
           n3 = n3/f_t
           !
           n1 = MLvector_product(n3,n2)
           !
           cosa3 = sum(n2(:)*v23(:))/sqrt(sum(v23(:)**2))
           !
           alpha3 = acos(cosa3)
           !
           cosphi =  ( cos(beta)-cos(alpha)*cos(alpha3) )/( sin(alpha)*sin(alpha3) )
           !
           phi = acos(cosphi)
           !
           x(iatom,:) = x(p1,:)+rbond*( cos(alpha)*n2(:) &
                                      + sin(alpha)*cos(phi)*n1(:) &
                                      + sin(alpha)*sin(phi)*n3(:) )
        case(101)
           !
           idihedral = idihedral + 1
           !
           zeta = molec%zmatrix(iatom)%connect(3)
           !
           idelta = 0
           kappa(3) = zeta
           do ikappa = 1,3
             if (ikappa==zeta) cycle
             idelta = idelta + 1
             kappa(idelta) = ikappa
           enddo
           !
           n1 = 0
           !
           n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
           n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral)
           !
           n1(kappa(3)) =sqrt(1.0_ark-(n1(kappa(1))**2+n1(kappa(2))**2))
           !
           x(iatom,:) = x(p1,:)+rbond*n1(:)
           !
           continue
           !  
        case(102)
           !
           idihedral = idihedral + 1
           !
           zeta = molec%zmatrix(iatom)%connect(3)
           !
           idelta = 0
           kappa(3) = zeta
           do ikappa = 1,3
             if (ikappa==zeta) cycle
             idelta = idelta + 1
             kappa(idelta) = ikappa
           enddo
           !
           n1 = 0
           !
           n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
           n1(kappa(2)) =-r(molec%Nbonds+molec%Nangles+idihedral)
           !
           n1(kappa(3)) =sqrt( 1.0_ark-( n1(kappa(1))**2+n1(kappa(2))**2 ) )
           !
           x(iatom,:) = x(p1,:)+rbond*n1(:)
           !
           continue
           !
        case(103)
           !
           idihedral = idihedral + 1
           !
           zeta = molec%zmatrix(iatom)%connect(3)
           !
           idelta = 0
           kappa(3) = zeta
           do ikappa = 1,3
             if (ikappa==zeta) cycle
             idelta = idelta + 1
             kappa(idelta) = ikappa
           enddo
           !
           n1 = 0
           !
           n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
           n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
           !
           n1(kappa(3)) = r(iatom-1)
           !
           x(iatom,:) = x(p1,:)+n1(:)
           !
           continue
           !
        case(105)
           !
           idihedral = idihedral + 1
           !
           zeta = molec%zmatrix(iatom)%connect(3)
           !
           idelta = 0
           kappa(3) = zeta
           do ikappa = 1,3
             if (ikappa==zeta) cycle
             idelta = idelta + 1
             kappa(idelta) = ikappa
           enddo
           !
           n1 = 0
           !
           n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
           n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
           !
           n1(kappa(3)) = r(iatom-1)
           !
           x(iatom,:) = x(p1,:)+n1(:)
           !
        case(104)
           !
           idihedral = idihedral + 1
           !
           zeta = molec%zmatrix(iatom)%connect(3)
           !
           idelta = 0
           kappa(3) = zeta
           do ikappa = 1,3
             if (ikappa==zeta) cycle
             idelta = idelta + 1
             kappa(idelta) = ikappa
           enddo
           !
           n1 = 0
           !
           n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
           n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
           !
           n1(kappa(3)) = sqrt(rbond**2-(n1(kappa(1))**2+n1(kappa(2))**2))
           !
           x(iatom,:) = x(p1,:)+n1(:)
           !
        case(106)
           !
           idihedral = idihedral + 1
           !
           zeta = molec%zmatrix(iatom)%connect(3)
           !
           idelta = 0
           kappa(3) = zeta
           do ikappa = 1,3
             if (ikappa==zeta) cycle
             idelta = idelta + 1
             kappa(idelta) = ikappa
           enddo
           !
           n1 = 0
           !
           n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
           n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
           !
           n1(kappa(3)) = sqrt(1.0_ark-(n1(kappa(1))**2+n1(kappa(2))**2))
           !
           x(iatom,:) = x(p1,:)+rbond*n1(:)
           !
           continue           
           !
        case(107)
           !
           idihedral = idihedral + 1
           !
           zeta = molec%zmatrix(iatom)%connect(3)
           !
           idelta = 0
           kappa(3) = zeta
           do ikappa = 1,3
             if (ikappa==zeta) cycle
             idelta = idelta + 1
             kappa(idelta) = ikappa
           enddo
           !
           n1 = 0
           !
           n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
           n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
           !
           n1(kappa(3)) = sqrt(rbond**2-(n1(kappa(1))**2+n1(kappa(2))**2))
           !
           x(iatom,:) = x(p1,:)+n1(:)
           !
        case(108)
           !
           idihedral = idihedral + 1
           !
           zeta = molec%zmatrix(iatom)%connect(3)
           !
           idelta = 0
           kappa(3) = zeta
           do ikappa = 1,3
             if (ikappa==zeta) cycle
             idelta = idelta + 1
             kappa(idelta) = ikappa
           enddo
           !
           n1 = 0
           !
           n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
           n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
           !
           n1(kappa(3)) = r(iatom-1)
           !
           x(iatom,:) = x(p1,:)+n1(:)
           !
        end select 
        !
     enddo
     !
     ! Find center of mass and correct the vectors
     !
     do ix = 1,3 
       !
       CM_shift = sum( x(:,ix)*molec%AtomMasses(:) )/sum( molec%AtomMasses(:) )
       x(:,ix) = x(:,ix) - CM_shift
       !
     enddo 
     !
     cartesian = x
     !
     if (verbose>=6) write(out,"('MLfromlocal2cartesian/end')") 
     !
   end subroutine MLfromlocal2cartesian
   !
   !  
   function MLvector_product(v1,v2) result (v)
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
   end function MLvector_product



  subroutine MLlinur(dimen,npar,coeff,constant,solution,error)

  integer(ik),intent(in)  :: dimen,npar
  integer(ik),intent(out) :: error 
  real(rk),intent(in)  :: coeff(npar,npar),constant(npar)
  real(rk),intent(out) :: solution(npar)
  real(rk)          :: a0(npar,npar)
  real(rk)          :: c
  integer(ik)                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    error=0
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
       error=i
       return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      if (a0(i,i).eq.0) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1
         c=0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine MLlinur


  subroutine MLlinurark(dimen,coeff,constant,solution,error)

  integer(ik),intent(in)  :: dimen
  real(ark),intent(in)  :: coeff(dimen,dimen),constant(dimen)
  real(ark),intent(out) :: solution(dimen)
  integer(ik),intent(out) :: error 
  !
  real(ark)          :: a0(dimen,dimen)
  real(ark)          :: c
  integer(ik)        :: i1,i2,i,k8,k,l,k9

  !----- begin ----!

    if (verbose>=6) write(out,"('MLlinurark...')")
    !
    error=0
    !
    if (dimen==1) then 
     !
     if (abs(coeff(1,1))<small_) then
        !write(out,"(i,'-th parameter is wrong, small a0(i,i) =  ',g18.8)") i,coeff(i,i)
        error=i
        return
     endif
     !
     solution(1) = constant(1)/coeff(1,1)
     !
     return
     !
    endif
    !
    a0 = coeff
    solution = constant
    !
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
       !write(out,"(i,'-th parameter is wrong, c >= a0(i,i) =  ',g18.8)") i,a0(i,i)
       error=i
       return
      endif
      !
      a0(i,i)=sqrt(a0(i,i)-c)
      if (abs(a0(i,i))<small_) then
         !write(out,"(i,'-th parameter is wrong, small a0(i,i) =  ',g18.8)") i,a0(i,i)
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1
         c=0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    !
    if (verbose>=6) write(out,"('...MLlinurark done!')")   
    !
  end subroutine MLlinurark




   subroutine MLdiag_ulen_II(n,a,d,ve,niter_,err_)
      !
      integer(ik),intent(in),optional  ::  niter_
      real(rk),intent(in),optional     ::  err_
      !
      integer(ik)  ::  p,q,p1,q1,irot,i,n2,kp,kp1,ip1,j,iq1,kq1,ii,n
      real(rk)     ::  a(n,n),d(n),ve(n,n)
      real(rk)     ::  sm,tresh,g,h,s,c,t,tau,theta,ff1,ff2,ff3,ff4
      real(rk),allocatable  ::  b(:),z(:)
      integer(ik)  ::  niter = 100
      real(rk)     ::  err
      !
      !err = job%tolerance
      !
      err = sqrt(small_)
      !
      if (present(niter_)) niter = niter_
      if (present(err_)) err = err_
      !
      allocate(b(n+10),z(n+10))
      !
      ve = 0 
      !
      do p=1,n
        ve(p,p)=0.0_rk
      enddo 
      !
      z = 0
      !
      do p=1,n
        d(p)=a(p,p)
        b(p)=d(p)
      enddo
      !
      irot=0
      do i=1,niter
        !
        sm=0.0_rk
        do p=1,n-1
           !
           sm = sm + sum(dabs(a(p,p+1:n))) 
           !
           !do q=p+1,n
           !  sm=sm+dabs(a(p,q))
           !enddo
           !
        enddo
        !
        if(sm<=err) exit

        tresh=0
        if(i<4) tresh=0.2_rk*sm/(real(n*n,rk))
        !
        do p1=1,n-1
           do q1=p1+1,n
              !
              g=100.0_rk*dabs(a(p1,q1))
              !
              ff1=dabs(d(p1)+g)
              ff2=dabs(d(p1))
              ff3=dabs(d(q1)+g)
              ff4=dabs(d(q1))
              !
              if(i<4.or.(ff1/=ff2).or.ff3/=ff4) then 
                 !
                 ff1=dabs(a(p1,q1))
                 if(ff1<=tresh) exit
                 !
                 h=d(q1)-d(p1)
                 !
                 ff1=dabs(h)+g
                 !
                 ff2=dabs(h)
                 !
                 if(ff1==ff2) then 
                   !
                   t=a(p1,q1)/h
                   !
                 else
                   !
                   theta=0.5_rk*h/a(p1,q1)
                   t=1.0_rk/(dabs(theta)+sqrt(1.0_rk+theta*theta))
                   !
                   if(theta<0.0_rk)  t=-t
                   !
                 endif
                 !
                 c=1.0_rk/sqrt(1.0_rk+t*t)
                 !
                 s=t*c
                 !
                 tau=s/(1.0_rk+c)
                 !
                 h=t*a(p1,q1)
                 !
                 z(p1)=z(p1)-h
                 z(q1)=z(q1)+h
                 d(p1)=d(p1)-h
                 d(q1)=d(q1)+h
                 !
                 a(p1,q1)=0.0_rk
                 !
                 do j=1,p1-1
                   g=a(j,p1)
                   h=a(j,q1)
                   a(j,p1)=g-s*(h+g*tau)
                   a(j,q1)=h+s*(g-h*tau)
                 enddo
                 !
                 do j=p1+1,q1-1
                   g=a(p1,j)
                   h=a(j,q1)
                   a(p1,j)=g-s*(h+g*tau)
                   a(j,q1)=h+s*(g-h*tau)
                 enddo
                 !
                 do j=q1+1,n
                   g=a(p1,j)
                   h=a(q1,j)
                   a(p1,j)=g-s*(h+g*tau)
                   a(q1,j)=h+s*(g-h*tau)
                 enddo
                 !
                 do j=1,n
                   g=ve(j,p1)
                   h=ve(j,q1)
                   ve(j,p1)=g-s*(h+g*tau)
                   ve(j,q1)=h+s*(g-h*tau)
                 enddo
                 irot=irot+1
              endif 
              a(p1,q1)=0.0_rk
              !
           enddo
        enddo
        !
        d = b + z
        b = d
        z = 0
        !
      enddo
      !
      deallocate (b,z)
      !


  deallocate(b,z)

  end subroutine MLdiag_ulen_II




   subroutine MLdiag_ulen(n,a,d,ve)
      !
      integer(ik)  ::  p,q,p1,q1,irot,i,n2,kp,kp1,ip1,j,iq1,kq1,ii,n
      real(rk)     ::  a(n,n),d(n),ve(n,n)
      real(rk)     ::  sm,tresh,g,h,s,c,t,tau,theta,ff1,ff2,ff3,ff4
      real(rk)     ::  err
      real(rk),allocatable  ::  b(:),z(:)
      !
      err = small_
      !
      allocate(b(n+10),z(n+10))


 101  format(5e14.5)
      do p=1,n
        do q=1,n
          ve(p,q)=0.0_rk
          if(p.eq.q) ve(p,q)=1.0_rk
        enddo
      enddo
      do p=1,n
        z(p)=0.0_rk
        d(p)=a(p,p)
        b(p)=d(p)
      enddo
      irot=0
      do 50 i=1,50
      sm=0.0_rk
      n2=n-1
      do p=1,n2
        kp=p+1
        do q=kp,n
          sm=sm+dabs(a(p,q))
        enddo
      enddo
      if(sm.le.err) goto 50
      tresh=0.0_rk
      if(i-4) 3,4,4
  3   tresh=0.2_rk*sm/(n*n)
  4     do 33 p1=1,n2
        kp1=p1+1
      do 33 q1=kp1,n
      g=100*dabs(a(p1,q1))
      ff1=dabs(d(p1)+g)
      ff2=dabs(d(p1))
      ff3=dabs(d(q1)+g)
      ff4=dabs(d(q1))
      if(i.gt.4.and.(ff1.eq.ff2).and.ff3.eq.ff4) goto 7
      ff1=dabs(a(p1,q1))
        if(ff1.le.tresh) goto 33
      h=d(q1)-d(p1)
      ff1=dabs(h)+g
      ff2=dabs(h)
      if(ff1.ne.ff2) goto 13
      t=a(p1,q1)/h
      goto 6
  13    theta=0.5_rk*h/a(p1,q1)
        t=1._rk/(dabs(theta)+sqrt(1._rk+theta*theta))
      if(theta) 5,6,6
  5     t=-t
  6     c=1._rk/sqrt(1._rk+t*t)
        s=t*c
      tau=s/(1._rk+c)
      h=t*a(p1,q1)
      z(p1)=z(p1)-h
      z(q1)=z(q1)+h
      d(p1)=d(p1)-h
      d(q1)=d(q1)+h
      a(p1,q1)=0.0_rk
      ip1=p1-1
        do j=1,ip1
          g=a(j,p1)
          h=a(j,q1)
          a(j,p1)=g-s*(h+g*tau)
          a(j,q1)=h+s*(g-h*tau)
        enddo
        iq1=q1-1
        do 21 j=kp1,iq1
        g=a(p1,j)
        h=a(j,q1)
        a(p1,j)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  21    continue
        kq1=q1+1
        do 26 j=kq1,n
        g=a(p1,j)
        h=a(q1,j)
        a(p1,j)=g-s*(h+g*tau)
        a(q1,j)=h+s*(g-h*tau)
  26    continue
          do 29 j=1,n
        g=ve(j,p1)
        h=ve(j,q1)
        ve(j,p1)=g-s*(h+g*tau)
        ve(j,q1)=h+s*(g-h*tau)
  29      continue
        irot=irot+1
  7     a(p1,q1)=0.0_rk
  33    continue
        do 44 ii=1,n
      d(ii)=b(ii)+z(ii)
      b(ii)=d(ii)
      z(ii)=0.0d0
  44  continue
  50  continue



  deallocate(b,z)

  end subroutine MLdiag_ulen





   subroutine MLdiag_ulen_ark(n,a,d,ve)
      !
      integer(ik)  ::  p,q,p1,q1,irot,i,n2,kp,kp1,ip1,j,iq1,kq1,ii,n
      real(ark)     ::  a(n,n),d(n),ve(n,n)
      real(ark)     ::  sm,tresh,g,h,s,c,t,tau,theta,ff1,ff2,ff3,ff4
      real(ark)     ::  err
      real(ark),allocatable  ::  b(:),z(:)
      !
      err = small_
      !
      allocate(b(n+10),z(n+10))


 101  format(5e14.5)
      do 10 p=1,n
      do 10 q=1,n
      ve(p,q)=0.0_ark
      if(p.eq.q) ve(p,q)=1.0_ark
  10  continue
      do 99 p=1,n
      z(p)=0.0_ark
      d(p)=a(p,p)
      b(p)=d(p)
 99   continue
      irot=0
      do 50 i=1,50
      sm=0.0_ark
      n2=n-1
      do 30 p=1,n2
      kp=p+1
      do 30 q=kp,n
      sm=sm+abs(a(p,q))
  30  continue
      if(sm.le.err) goto 50
      tresh=0.0_ark
      if(i-4) 3,4,4
  3   tresh=0.2_ark*sm/(n*n)
  4     do 33 p1=1,n2
        kp1=p1+1
      do 33 q1=kp1,n
      g=100*abs(a(p1,q1))
      ff1=abs(d(p1)+g)
      ff2=abs(d(p1))
      ff3=abs(d(q1)+g)
      ff4=abs(d(q1))
      if(i.gt.4.and.(ff1.eq.ff2).and.ff3.eq.ff4) goto 7
      ff1=abs(a(p1,q1))
        if(ff1.le.tresh) goto 33
      h=d(q1)-d(p1)
      ff1=abs(h)+g
      ff2=abs(h)
      if(ff1.ne.ff2) goto 13
      t=a(p1,q1)/h
      goto 6
  13    theta=0.5_ark*h/a(p1,q1)
        t=1._ark/(abs(theta)+sqrt(1._ark+theta*theta))
      if(theta) 5,6,6
  5     t=-t
  6     c=1._ark/sqrt(1._ark+t*t)
        s=t*c
      tau=s/(1._ark+c)
      h=t*a(p1,q1)
      z(p1)=z(p1)-h
      z(q1)=z(q1)+h
      d(p1)=d(p1)-h
      d(q1)=d(q1)+h
      a(p1,q1)=0.0_ark
      ip1=p1-1
        do 20 j=1,ip1
        g=a(j,p1)
        h=a(j,q1)
        a(j,p1)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  20      continue
        iq1=q1-1
        do 21 j=kp1,iq1
        g=a(p1,j)
        h=a(j,q1)
        a(p1,j)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  21    continue
        kq1=q1+1
        do 26 j=kq1,n
        g=a(p1,j)
        h=a(q1,j)
        a(p1,j)=g-s*(h+g*tau)
        a(q1,j)=h+s*(g-h*tau)
  26    continue
          do 29 j=1,n
        g=ve(j,p1)
        h=ve(j,q1)
        ve(j,p1)=g-s*(h+g*tau)
        ve(j,q1)=h+s*(g-h*tau)
  29      continue
        irot=irot+1
  7     a(p1,q1)=0.0_ark
  33    continue
        do 44 ii=1,n
      d(ii)=b(ii)+z(ii)
      b(ii)=d(ii)
      z(ii)=0.0d0
  44  continue
  50  continue



  deallocate(b,z)

  end subroutine MLdiag_ulen_ark



  function three_j0(j1,j2,j3,k1,k2,k3)

      real(rk) :: three_j0
      !
      integer(ik) :: j1,j2,j3,k1,k2,k3,newmin,newmax,new,iphase
      real(rk)   :: a,b,c,al,be,ga,delta,clebsh,minus
      real(rk)   :: term,term1,term2,term3,summ,dnew,term4,term5,term6

!
      a = j1
      b = j2
      c = j3
      al= k1
      be= k2
      ga= k3

      three_j0=0
!
!     (j1+j2).ge.j and j.ge.abs(a-b)    -m=m1+m2    j1,j2,j.ge.0
!     abs(m1).le.j1    abs(m2).le.j2   abs(m).le.j
!
      if(c.gt.a+b) return
      if(c.lt.abs(a-b)) return
      if(a.lt.0.or.b.lt.0.or.c.lt.0) return
      if(a.lt.abs(al).or.b.lt.abs(be).or.c.lt.abs(ga)) return
      if(-1.0_rk*ga.ne.al+be) return
!
!
!     compute delta(abc)
!
      delta=sqrt(fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.0_rk))
!
!
      term1=fakt(a+al)*fakt(a-al)
      term2=fakt(b-be)*fakt(b+be)
      term3=fakt(c+ga)*fakt(c-ga)
      term=sqrt((2.0_rk*c+1.0_rk)*term1*term2*term3)
!
!
!     now compute summation term
!
!     sum to get summation in eq(2.34) of brink and satchler.  sum until
!     a term inside factorial goes negative.  new is index for summation
!     .  now find what the range of new is.
!
!
      newmin=idnint(max((a+be-c),(b-c-al),0.0_rk))
      newmax=idnint(min((a-al),(b+be),(a+b-c)))
!
!
      summ=0
!
!
      do new=newmin,newmax
        dnew=real(new,rk)
        term4=fakt(a-al-dnew)*fakt(c-b+al+dnew)
        term5=fakt(b+be-dnew)*fakt(c-a-be+dnew)
        term6=fakt(dnew)*fakt(a+b-c-dnew)
        summ=summ+(-1.0_rk)**new/(term4*term5*term6)
      enddo
!
!     so clebsch-gordon <j1j2m1m2ljm> is clebsh
!
      clebsh=delta*term*summ/sqrt(10.0_rk)
!
!     convert clebsch-gordon to three_j0
!
      iphase=idnint(a-b-ga)
      minus = -1.0_rk
      if (mod(iphase,2).eq.0) minus = 1.0_rk
      three_j0=minus*clebsh/sqrt(2.0_rk*c+1.0_rk)

!     threej=(-1.d0)**(iphase)*clebsh/sqrt(2.0_rk*c+1.d0)
!
!
   end function three_j0

   function three_j(j1,j2,j3,k1,k2,k3)

      real(rk) :: three_j
      integer(ik),intent(in) :: j1,j2,j3,k1,k2,k3
      real(rk) :: a,b,c,al,be,ga
      !
      integer(ik):: newmin,newmax,new,iphase
      real(rk)   :: delta,clebsh,minus
      real(rk)   :: term,term1,term2,term3,summ,dnew,term4,term5,term6,delta_log,term16,termlog

      a = j1
      b = j2
      c = j3
      al= k1
      be= k2
      ga= k3

      three_j=0
!
!     (j1+j2).ge.j and j.ge.abs(a-b)    -m=m1+m2    j1,j2,j.ge.0
!     abs(m1).le.j1    abs(m2).le.j2   abs(m).le.j
!
      if(c.gt.a+b) return
      if(c.lt.abs(a-b)) return
      if(a.lt.0.or.b.lt.0.or.c.lt.0) return
      if(a.lt.abs(al).or.b.lt.abs(be).or.c.lt.abs(ga)) return
      if(-1.0_rk*ga.ne.al+be) return
!
!
!     compute delta(abc)
!
!     delta=sqrt(fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.0_rk))
      delta_log = faclogf(a+b-c)+faclogf(a+c-b)+faclogf(b+c-a)-faclogf(a+b+c+1.0_rk)
      !
      delta=sqrt(exp(delta_log)) 
!
!
      !term1=fakt(a+al)*fakt(a-al)
      !term2=fakt(b-be)*fakt(b+be)
      !term3=fakt(c+ga)*fakt(c-ga)
      !
      !term=sqrt( (2.0_rk*c+1.0_rk)*term1*term2*term3 )
      !
      !
      term1=faclogf(a+al)+faclogf(a-al)
      term2=faclogf(b-be)+faclogf(b+be)
      term3=faclogf(c+ga)+faclogf(c-ga)
      !
      termlog = ( term1+term2+term3+delta_log )*0.5_rk
 
      term=sqrt( (2.0_rk*c+1.0_rk) )
!
!
!     now compute summation term
!
!     sum to get summation in eq(2.34) of brink and satchler.  sum until
!     a term inside factorial goes negative.  new is index for summation
!     .  now find what the range of new is.
!
!
      newmin=idnint(max((a+be-c),(b-c-al),0.0_rk))
      newmax=idnint(min((a-al),(b+be),(a+b-c)))
!
!
      summ=0
!
!
      do new=newmin,newmax
        !
        dnew=real(new,rk)
        !
        term4=faclogf(a-al-dnew)+faclogf(c-b+al+dnew)
        term5=faclogf(b+be-dnew)+faclogf(c-a-be+dnew)
        term6=faclogf(dnew)+faclogf(a+b-c-dnew)
        !
        term16=termlog-(term4+term5+term6)
        !
        summ=summ+(-1.0_rk)**new*exp(term16)
        !
      enddo
!
!     so clebsch-gordon <j1j2m1m2ljm> is clebsh
!
      clebsh=term*summ ! /sqrt(10.0_rk)
!
!     convert clebsch-gordon to three_j
!
      iphase=idnint(a-b-ga)
      minus = -1.0_rk
      if (mod(iphase,2).eq.0) minus = 1.0_rk
      three_j=minus*clebsh/term

!     threej=(-1.d0)**(iphase)*clebsh/sqrt(2.0_rk*c+1.d0)
!
!
   end function three_j



   function fakt(a) result (f)

      real(rk),intent(in) :: a
      real(rk)            :: ax,f
      integer(ik)         :: i,ic
!

!
      ax=a
      f=1.0_rk
      if(abs(ax)<1.d-24) return
      f=.1d0
      if(ax.lt.0.d0) then 
         write (*,"(' fkt.err  negative argument for functi on fakt. argument = ',e12.5)") ax
         stop 'fkt.err  negative argument'
      endif 
      !
      ic=idnint(ax)
      ax=ax/10.0_rk
      f=ax
      do  i=1,ic-1
        f=f*(ax-real(i,rk)*0.1_rk)
      enddo

    end function fakt

   !
   ! 3-dim rotation by three Euler angles
   !
   function ML_euler_rotait(theta,phi,chi) result (f)

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

   end function ML_euler_rotait

   
   function aacos(x,txt) result (f)
      !
      real(ark),intent(in)     :: x
      character(cl),intent(in),optional :: txt
      real(ark)                :: f
      !
      if ( abs(x)>1.0_ark+100.0_ark*sqrt(small_) ) then 
         !
         write (out,"('|cos(x)|>1: ',f18.8,a)") x,txt
         stop 'aacos - bad cosalpha'
         !
      elseif ( x>=1.0_ark) then
         !
         f = 0.0_ark
         !
      elseif ( x<=-1.0_ark) then 
         f = pi
      else 
         f = acos(x)
      endif
      !
   end function aacos


   function aasin(x,txt) result (f)
      !
      real(ark),intent(in)     :: x
      character(len=cl),intent(in),optional :: txt
      real(ark)                :: f
      !
      if ( abs(x)>1.0_ark+100.0_ark*sqrt(small_) ) then 
         !
         write (out,"('|cos(x)|>1: ',f18.8,a)") x,txt
         stop 'aasin - bad sinalpha'
         !
      elseif ( x>=1.0_ark) then
         !
         f = pi*0.5_ark
         !
      elseif ( x<=-1.0_ark) then 
         f = pi*1.5_ark
      else 
         f = asin(x)
      endif
      !
   end function aasin

  !
  !
  ! calculate factorial by log function 
  ! 

  function faclog(k)   result (v)
    integer(ik),intent(in) ::  k
    real(ark)              :: v 
    integer(ik) j

    v=0
    if(k>=2) then 
      do j=2,k
         v=v+log(real(j,ark))
      enddo 
    endif 
    
  end function faclog



  function faclogf(a)   result (v)
    real(rk),intent(in) ::  a
    real(ark)              :: v 
    integer(ik) j,k

    v=0
    k=nint(a)
    if(k>=2) then 
      do j=2,k
         v=v+log(real(j,ark))
      enddo 
    endif 
    
  end function faclogf


  subroutine ML_rjacobi_fit_ark(N,M,a,b,x,tolerance)
    !
    integer(ik),intent(in) :: N,M
    real(ark),intent(in)  :: a(N,M),b(N)
    real(ark),intent(inout) :: x(M)
    real(rk),intent(in) :: tolerance
    !
    real(ark) :: x0(m),eps(n)
    !
    real(ark),allocatable :: rjacob(:,:),b_(:),am(:,:),bm(:),b_r(:),b_l(:),cm(:)
    !
    real(ark) :: stadev_old,stability,stadev,ssq,stadev_best,h
    !
    double precision,allocatable :: ad(:,:),bd(:,:)
    !
    integer(ik),parameter :: itmax=100
    integer(ik) :: iter,k,irow,icolumn,ierror,alloc
    !
    if (verbose>=5) then
      write(out,"('ML_rjacobi_fit_ark start ')") 
    endif
    !
    rjacob = 0 
    iter = 0
    stadev_old = 2.e10
    stability =  1.e10
    stadev    =  1.e10
    !
    x0 = x
    !
    stadev_best = tolerance ! (small_)*1e-2
    !
    allocate(rjacob(n,m),stat=alloc)
    call ArrayStart('ML_rjacobi_fit_ark',alloc,1_ik,ark,size(rjacob,kind=hik))
    allocate(am(m,m),stat=alloc)
    call ArrayStart('ML_rjacobi_fit_ark',alloc,1_ik,ark,size(am,kind=hik))
    allocate(bm(m),stat=alloc)
    call ArrayStart('ML_rjacobi_fit_ark',alloc,1_ik,ark,size(bm,kind=hik))
    allocate(cm(m),stat=alloc)
    call ArrayStart('ML_rjacobi_fit_ark',alloc,1_ik,ark,size(cm,kind=hik))
    !
    allocate(b_(n),b_r(n),b_l(n),stat=alloc)
    call ArrayStart('ML_rjacobi_fit_ark',alloc,1_ik,ark,3*size(b_,kind=hik))
    !
    allocate(ad(m,m),stat=alloc)
    call ArrayStart('ML_rjacobi_fit_ark',alloc,1_ik,rk,size(ad,kind=hik))
    allocate(bd(m,1),stat=alloc)
    call ArrayStart('ML_rjacobi_fit_ark',alloc,1_ik,rk,size(bd,kind=hik))
    !
    iter = 0 
    !
    outer_loop: & 
    do while( iter<itmax .and. stadev>stadev_best .and. abs(stadev-stadev_old)>small_ )
      !
      stadev_old=stadev
      !
      iter = iter + 1
      ssq=0
      !
      ! Calculate the target 
      !
      b_ = matmul(a,x)
      !
      eps(:) = b(:)-b_(:)
      !
      do k = 1,m
        !
        h = 1.e-5*abs(x(k)) ; if (h<sqrt(small_)) h = sqrt(small_)
        !
        x(k) = x(k) + h 
        !
        b_r = matmul(a,x)
        !
        x(k) = x(k) - h - h
        !
        b_l = matmul(a,x)
        !
        rjacob(:,k)  = ( b_r(:)-b_l(:))/h*0.5_ark
        !
        x(k) = x(k) + h
        !
      enddo 
      !
      ssq=sqrt(sum(eps(:)**2))
      !
      ! We constract a set of linear equations A x = B
      !
      ! form A matrix 
      !
      do irow=1,m       !==== row-...... ====!
        do icolumn=1,irow    !==== column-....====!
          am(irow,icolumn)=sum(rjacob(:,icolumn)*rjacob(:,irow))
          am(icolumn,irow)=am(irow,icolumn)
        enddo
      enddo
      !
      ! form B matrix 
      !
      do irow=1,m       !==== row-...... ====!
        bm(irow)=sum(eps(:)*rjacob(:,irow))
      enddo   
      !
      ! Solve the set of linear equations 
      !
      call MLlinurark(m,am,bm,cm,ierror)
      !
      if (ierror>0) then
        !
        if (verbose>=5) then
          write(out,"('ML_rjacobi_fit_ark error: ',3i7,g12.4)") n,m,ierror,tolerance
        endif
        !
        ad = am ; bd(:,1) = bm(:) 
        !
        call lapack_gelss(ad(:,:),bd(:,:))
        !
        cm(:) = bd(:,1)
        !
      endif
      !
      x(:)=x(:)+cm(:)
      !
      stadev=ssq/sqrt(real(n,ark))
      !
    enddo  outer_loop ! --- iter
    !
    b_ = matmul(a,x)
    !
    eps(:) = b(:)-b_(:)
    !
    ssq=sqrt(sum(eps(:)**2))
    !
    stadev=ssq/sqrt(real(n,ark))
    !
    if (iter==itmax.and.verbose>=5.and.stadev>sqrt(small_)) then
       write(6,"('ML_rjacobi_fit_ark: could not find solution after ',i8,' iterations')") iter
       stop 'ML_rjacobi_fit_ark: could not find solution'
    endif 
    !
    deallocate(rjacob,am,bm,cm,b_,b_r,b_l,ad,bd)
    !
    call ArrayStop('ML_rjacobi_fit_ark')
    !
    if (verbose>=5) then
      write(out,"('ML_rjacobi_fit_ark end, rms = ',g18.4)") stadev  
    endif
    !
  end subroutine ML_rjacobi_fit_ark

  subroutine ML_spline(x,y,n,yp1,ypn,y2)
  ! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
  ! x1 < x2 < ... < xN, and given values yp1 and ypn for the 1st derivative of the
  ! interpolating function at points 1 and n, respectively, this routine returns an array y2(1:n)
  ! of length n which contains the second derivatives of the interpolating function at
  ! the tabulated points xi. If yp1 and/or ypn are equal to 1e30 or larger, the routine is
  ! signaled to set the corresponding boundary condition for a natural spline, with zero second
  ! derivative on that boundary.
  ! From Numerical Recipes in Fortran 77
  !
  integer(ik)         :: n,i,k
  real(kind=ark)      :: yp1,ypn,x(n),y(n),y2(n)
  real(kind=ark)      :: p,qn,sig,un,u(n)

  ! if (yp1.gt.safe_max) then   ! Lorenzo Lodi, 13 February 2014; modified so that natural splines are always used
    y2(1)=0._ark
    u(1)=0._ark
  ! else
  !   y2(1)=-0.5_ark
  !   u(1)=(3._ark/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  ! endif
  
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2._ark
    y2(i)=(sig-1._ark)/p
    !
    if (abs(x(i+1)-x(i))<sqrt(epsilon(1.0_ark))) then
      write(out, "(A,2i5,' x(i),x(i+1) = ',2g16.8,' y(i),y(i+1) = ',2g16.8)") &
                'Error in spline: identical grid points: i,i+1 = ', i,i+1, x(i),x(i+1),y(i),y(i+1)
      stop 'Error in spline: identical grid points'
    endif
    !
    u(i)=(6._ark*( ( y(i+1)-y(i) )/( x(i+1)-x(i) )-( y(i)-y(i-1) )/( x(i)-x(i-1)) )/( x(i+1)-x(i-1) )-sig*u(i-1) )/p
  enddo
  ! if (ypn.gt.safe_max) then   ! Lorenzo Lodi, 13 February 2014; modified so that natural splines are always used
    qn=0._ark
    un=0._ark
  ! else
  !   qn=0.5_ark
  !   un=(3._ark/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  ! endif
  !
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1._ark)
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
  return
  end  subroutine ML_spline

  subroutine ML_splint(xa,ya,y2a,n,x,y)
  !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
  ! xai's in order), and given the array y2a(1:n), which is the output from spline above,
  ! and given a value of x, this routine returns a cubic-spline interpolated value y.
  ! From Numerical Recipes in Fortran 77
  
  integer :: n,klo,khi,k
  real(kind=ark) :: x,y,xa(n),y2a(n),ya(n)
  real(kind=ark) :: a,b,h
  klo=1
  khi=n
  
  1 if (khi-klo.gt.1) then
    k=(khi+klo)/2
    if(xa(k).gt.x)then
      khi=k
    else
      klo=k
    endif
  goto 1
  endif
  h=xa(khi)-xa(klo)
  if (h.eq.0._ark) stop 'bad xa input in splint'
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6._ark
  return
  end subroutine ML_splint
  
  subroutine ML_splint_quint(xa,ya,n,x,y, b,c,d,e,f)
  !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
  ! xai's in order), and given the arrays B,C,D,E,F (1:n), which is the output from spline above,
  ! and given a value of x, this routine returns a quintic-spline interpolated value y.
  ! Adapted from Numerical Recipes in Fortran 77
  
  integer, intent(in) :: n
  real(kind=ark) , intent(in) :: xa(n),ya(n)
  real(kind=ark) , intent(in)  :: x
  real(kind=ark) , intent(out) :: y
  real(kind=ark) ::  b(n),c(n),d(n),e(n), f(n)
  real(kind=ark) ::  p, q
  integer :: klo,khi,k, i
  real(kind=ark) :: h
  klo=1
  khi=n
  
  1 if (khi-klo > 1) then
    k=(khi+klo)/2
    if(xa(k) > x)then
      khi=k
    else
      klo=k
    endif
    goto 1
  endif
  
  h=xa(khi)-xa(klo)
  if (h == 0._ark) stop 'bad xa input in splint'
  
  i = klo
  p = x -xa(klo)
  q = xa(khi) - x
  y = ((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+ya(I)
  
  return
  end subroutine ML_splint_quint
  !
end module moltype
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

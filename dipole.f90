module dipole
! set dipole_debug > 2 with small expansions only
#define dipole_debug 0
                              ! 0 - none
                              ! 1 - some checkings only (no printing)
                              ! 2 - minimal printing
                              ! 3 - extendent printing

 use accuracy,     only : hik, ik, rk, ark, cl, wl, out, vellgt, planck, avogno, boltz, pi, small_, rad
 use fields,       only : job,analysis,bset
 use timer,        only : IOstart,IOStop,Arraystart,Arraystop,Timerstart,Timerstop,MemoryReport,TimerReport,&
                          TimerProbe,memory_limit,memory_now
 use molecules,    only : MLcoord_direct,MLrotsymmetry_generate,ddlmn_conj,dlmn,Phi_rot,calc_phirot
 use moltype,      only : manifold,molec, extF, intensity, three_j
 use symmetry,     only : sym


 use tran,         only : TReigenvec_unit, bset_contr, read_contrind, & 
                          read_eigenval,index_correlation,eigen,Neigenlevels,istate2ilevel,TReigenvec_unit_external

 private
 public dm_tranint,dm_analysis_density,dm_analysis_contribution

 real(rk),allocatable,save  :: dipole_me(:,:,:)


 type DmatT
      real(rk),pointer   :: mat(:,:) 
 end type DmatT

 type Dmat3T
      real(rk),pointer   :: mat(:,:,:) 
 end type Dmat3T

 type DimatT
      integer(ik),pointer   :: imat(:,:,:)
 end type DimatT

 type DkmatT
      integer(ik),pointer   :: kmat(:,:)
 end type DkmatT
 ! 
 type dipoleT                      
    real(rk),pointer    :: rot(:,:,:,:,:)
 end type dipoleT
 !
 type(dipoleT),allocatable     :: wigner(:,:) ! Rotational component of the dipole moment matrix elements 

contains





 subroutine dm_tranint

 ! compute  partition function, electric dipole transition moments, linestrength, and intensities
    !
    implicit none
    !
    integer(ik)          :: info

    integer(ik)              :: Jmin, Jmax, nJ, jind, j, ierror
    integer(ik), allocatable :: Jval(:)

    real(rk)             :: exp_en, part, beta, energy, q_part(20,0:100)

    integer(ik)          :: ilevel, irrep




    ! initialize array of J values
    ! J=0 is always present no matter if we use it in intensity calcs or not. 


    Jmin = minval(intensity%J(1:2))
    Jmax = maxval(intensity%J(1:2))
    nJ   = Jmax - Jmin + 1
    !
    if (Jmin>0) nJ = nJ + 1
    !
    allocate(Jval(nJ), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: Jval - out of memory'
    !
    Jval = 0 
    jind = 1
    do j = max(Jmin,1), Jmax
       jind = jind + 1
       Jval(jind) = j
    end do
    !
    select case (trim(intensity%action))
    !
    case('ABSORPTION', 'EMISSION', 'TM')
       !
       ! read contraction indexes (icase,ilambda) for all J specified
       !
       call read_contrind(nJ, Jval(1:nJ))
       !
       ! find correlation between indexes for J = 0 and J > 0
       ! We need this because the contracted matrix elements of the 
       ! dipole moment are computed on the J=0 conatracted basis functions. 
       ! When J/=0 the numbering (bookkeeping) of these basis set has changed and
       ! we need to find the correlation between the bookkeepings.  
       !
       call index_correlation(nJ, Jval(1:nJ))
       !
       ! read eigenvalues and their labeling, i.e. description;
       ! initialize file-units for reading eigenvectors
       !
       call read_eigenval(nJ, Jval(1:nJ),ierror)
       !
       if (ierror/=0) then 
           write(out,"('dm_tranint: read_eigenval error, some eigen-files do not exists')")
           stop 'dm_tranint: read_eigenval error, some eigenfiles are missing'
       endif 
       !
       !restore vibrational contracted matrix elements
       !for all  dipole moment vector components
       !
       call restore_vib_matrix_elements
       !
       if (intensity%tdm_replace) then
         !
         call replace_vib_trans_dipole_moments
         !
       endif
       !
       ! Run intensity simulations 
       !
       if (job%IOvector_symm) then
          call dm_intensity_symmvec(Jval)
       else
          call dm_intensity(Jval)
       endif
       !
       continue 
       !
       write(out, '(/a)') 'done'
       !

    !compute partition function
    !
    case('PARTFUNC')

       write(out, '(/a)') 'compute partition function'


       !read contraction indexes for all J specified;
       !read eigenvalues


       call read_contrind(nJ, Jval(1:nJ))
       call read_eigenval(nJ, Jval(1:nJ),ierror)
       !
       if (ierror/=0) then 
           write(out,"('dm_tranint: read_eigenval error, some eigen-files do not exists')")
           stop 'dm_tranint: read_eigenval error, some eigenfiles are missing'
       endif 


       write(out, '(/3x, a, 1x, a, (/1x, i3, 11x, i6))') 'J', 'number of levels',                                  &
                                                             Jval(nJ), Neigenlevels
       !
       beta = planck * vellgt / (boltz * intensity%temperature)
       !
       !loop ove J values
       !
       q_part  = 0 
       !
       do ilevel = 1, Neigenlevels

             j = eigen(ilevel)%jval
             energy = eigen(ilevel)%energy
             irrep  = eigen(ilevel)%igamma
             exp_en = exp(-(energy - intensity%ZPE) * beta)
             !q_part = q_part + intensity%gns(irrep) * real(2*j+1,rk) * exp_en
             q_part(irrep,j) = q_part(irrep,j) + intensity%gns(irrep) * real(2*j+1,rk) * exp_en

       end do

       part = sum(q_part(:,:))

       do j=0,Jmax
         do irrep = 1,sym%Nrepresen
            write(out, '(2i4,es16.8)') irrep,j,q_part(irrep,j)
         enddo
       enddo 
       !
       write(out, '(es16.8)') part
       !
       !irrep = sym%Nrepresen
       !write(out, '(/1x, a, 3x, es16.8/<irrep>(/1x, a, i1, a, 1x, es16.8)/1x, a, 4x, es16.8//a)') &
       !'zpe', intensity%ZPE, ('qpart(', irrep, ')', q_part(irrep), irrep = 1, sym%Nrepresen), 'qpart', sum(q_part), 'done'

       !write(out, '(/1x, a, 3x, es16.8/1x, a, 1x, es16.8//a)') 'zpe', intensity%ZPE, 'qpart', q_part, 'done'


    end select

    !
    call MemoryReport
    !
    call TimerReport



 end subroutine dm_tranint




 subroutine dm_analysis_density
 !
 ! Analysis of the eigenfunctions and of the prob. densities

    !
    implicit none
    !
    integer(ik)          :: info

    integer(ik)              :: Jmin, Jmax, nJ, jind, j, ierror
    integer(ik), allocatable :: Jval(:)

    integer(ik)          :: i,ilist,nlist,imode,kmode,jmode

    ! initialize array of J values
    ! J=0 is always present no matter if we use it in intensity calcs or not. 
    !
    i = 1
    Jmin = 1000 ; Jmax = 0
    !
    do while (i<100.and.analysis%j_list(i)/=-1)
      !
      j = analysis%j_list(i)
      !
      Jmin = min(Jmin,j)
      Jmax = max(Jmax,j)
      i = i + 1
      !
    enddo
    !
    nJ   = Jmax - Jmin + 1
    !
    if (Jmin>0) nJ = nJ + 1
    !
    allocate(Jval(nJ), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: Jval - out of memory'
    !
    Jval = 0 
    jind = 1
    do j = max(Jmin,1), Jmax
       jind = jind + 1
       Jval(jind) = j
    end do
    !
    !
    call read_contrind(nJ, Jval(1:nJ))
    !
    call index_correlation(nJ, Jval(1:nJ))
    !
    ! read eigenvalues and their labeling, i.e. description;
    ! initialize file-units for reading eigenvectors
    !
    call read_eigenval(nJ, Jval(1:nJ),ierror)
    !
    if (ierror/=0) then 
        write(out,"('dm_tranint: read_eigenval error, some eigen-files do not exists')")
        stop 'dm_tranint: read_eigenval error, some eigenfiles are missing'
    endif 
    !
    ilist = 1 
    !
    do while (ilist<size(analysis%dens_list).and.analysis%dens_list(ilist)/=-1)
      !
      ilist = ilist + 1 
      !
    enddo
    !
    nlist = ilist-4 ! minloc(analysis%dens_list(:),dim=1)-4
    !
    if (nlist>size(analysis%dens_list)-3) then 
      !
      write(out,"('dm_analysis_density: too many entries in dens_list: ',i8,' max - ',i8)") nlist,size(analysis%dens_list)-3
      stop 'dm_analysis_density: too many entries in dens_list'
      !
    endif 
    !
    if (analysis%reducible_eigen_contribution.or.analysis%population) then 
      !
      call DM_contribution_symmvec(Jval)
      !
    else 
      !
      imode = analysis%dens_list(nlist+1) ; jmode = analysis%dens_list(nlist+2) ; kmode = analysis%dens_list(nlist+3)
      !
      call DM_density_symmvec(Jval)
      !
      call DM_contribution_symmvec(Jval)
      !
    endif
    !
    write(out, '(/a)') 'done'
    !
    call MemoryReport
    !
    call TimerReport
    !
 end subroutine dm_analysis_density



 subroutine dm_analysis_contribution
 !
 ! Analysis of the reducible functions' contribution to the eigenfunction
    !
    implicit none
    !
    integer(ik)          :: info
    !
    integer(ik)              :: Jmin, Jmax, nJ, jind, j, ierror
    integer(ik), allocatable :: Jval(:)

    integer(ik)          :: i,ilist,nlist,imode,kmode,jmode

    ! initialize array of J values
    ! J=0 is always present no matter if we use it in intensity calcs or not. 
    !
    i = 1
    Jmin = 1000 ; Jmax = 0
    !
    do while (i<100.and.analysis%j_list(i)/=-1)
      !
      j = analysis%j_list(i)
      !
      Jmin = min(Jmin,j)
      Jmax = max(Jmax,j)
      i = i + 1
      !
    enddo
    !
    nJ   = Jmax - Jmin + 1
    !
    if (Jmin>0) nJ = nJ + 1
    !
    allocate(Jval(nJ), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: Jval - out of memory'
    !
    Jval = 0 
    jind = 1
    do j = max(Jmin,1), Jmax
       jind = jind + 1
       Jval(jind) = j
    end do
    !
    !
    call read_contrind(nJ, Jval(1:nJ))
    !
    call index_correlation(nJ, Jval(1:nJ))
    !
    ! read eigenvalues and their labeling, i.e. description;
    ! initialize file-units for reading eigenvectors
    !
    call read_eigenval(nJ, Jval(1:nJ),ierror)
    !
    if (ierror/=0) then 
        write(out,"('dm_tranint: read_eigenval error, some eigen-files do not exists')")
        stop 'dm_tranint: read_eigenval error, some eigenfiles are missing'
    endif
    !
    call DM_contribution_symmvec(Jval)
    !
    write(out, '(/a)') 'done'
    !
    call MemoryReport
    !
    call TimerReport
    !
 end subroutine dm_analysis_contribution


 subroutine DM_contribution_symmvec(Jval)
    implicit none
    !
    integer(ik),intent(in)  :: Jval(:)
    !
    integer(ik)    :: nJ,dimenmax
    integer(ik)    :: ilevelI, ilevelF, ndegI, ndegF, idegI, idegF, irec, idimen, nsizeF,nsizeI
    integer(ik)    :: irep
    integer(ik)    :: nmodes,info,indI,indF,jI,jF,Nrepresen
    integer(ik)    :: igammaI,igammaF,quantaI(0:molec%nmodes),quantaF(0:molec%nmodes),normalI(0:molec%nmodes)
    integer(ik)    :: dimenI,dimenF,unitI,unitF,jmax,kI,kF,external_vec_unit
    real(rk)       :: energyI, energyF, f_t,ener_,population
    logical        :: passed
    !
    integer(ik), allocatable :: icoeffI(:),istored(:),isave(:),isaved(:),ilevel_analyse(:)
    integer(ik), allocatable :: nlevelsG(:,:),ram_size(:,:),ilevelsG(:,:),iram(:,:)
    real(rk),    allocatable :: vecI(:), vecF(:),vec(:),vec_(:),rot_density(:,:,:,:), highest_contribution(:), ener_top_contri(:)
    complex(rk), allocatable :: dens_coord(:,:,:,:)
    real(rk),allocatable     :: density(:,:,:),Jrot_mat(:,:,:),vecE(:),vecE_(:)
    integer(ik), pointer :: Jeigenvec_unit(:,:)
    type(DmatT),pointer  :: vec_ram(:,:)
    type(DkmatT),pointer :: ijterm(:)
    !
    integer(ik)  :: jind,nlevels,igamma,nclasses,nsize_need,nlevelI,dimenmax_,nsizemax,itauI,itauF,termI,termF
    !
    integer(ik)  :: ram_size_
    !
    integer(ik)  :: iram_,Nterms,iterm
    !
    integer(hik) :: matsize,ram_size_max,dimenmax_ram,dimenmax_ram_
    integer(ik)  :: ilist,nlist,imode,kmode,jmode,dens_size(3),Nstored,i,k,i1,i2,i3,ipp,ip,npp,itheta,iphi
    real(rk)     :: rot(3),dens_,dtheta,dphi
    !
    integer(ik)  :: irootI,irow,ib,nelem,ielem,isrooti
    real(rk)     :: dtemp0 
    !
    double precision,allocatable  :: rw(:)
    double complex,allocatable  :: wd(:),rotmat(:,:),vl(:,:),vr(:,:),w(:),angular_m(:,:),crot_density(:,:,:,:)
    character(len=1)        :: jobvl,jobvr
    integer(ik)             :: cnu_i(0:molec%nmodes),ktau,tauI
    character(len=cl)       :: my_fmt !format for I/O specification
    character(len=wl)       :: my_fmt1 !format for I/O specification
    integer(ik), allocatable :: basis_set_contraction(:,:)
    integer(ik)                 :: iounit
    character(cl)               :: ioname
    !
    call TimerStart('Contribution analysis')
    !
    nmodes = molec%nmodes
    nclasses = size(eigen(1)%cgamma)-1
    Nrepresen = sym%Nrepresen
    !
    nJ = size(Jval)
    !
    ! Prepare the list of units with the stored eigenvectors
    !
    allocate(Jeigenvec_unit(nJ,sym%Nrepresen), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: Jeigenvec_unit - out of memory'
    !
    allocate(nlevelsG(nJ,sym%Nrepresen),ilevelsG(nJ,sym%Nrepresen),ram_size(nJ,sym%Nrepresen),iram(nJ,sym%Nrepresen),stat = info)
    call ArrayStart('nlevelsG',info,size(nlevelsG),kind(nlevelsG))
    call ArrayStart('ilevelsG',info,size(ilevelsG),kind(ilevelsG))
    call ArrayStart('ram_size',info,size(ram_size),kind(ram_size))
    call ArrayStart('iram',info,size(iram),kind(iram))
    !
    do jind=1,nJ
      do igamma = 1,sym%Nrepresen
        Jeigenvec_unit(jind,igamma) = TReigenvec_unit(jind,Jval,igamma)
      enddo
    enddo
    !
    ! maximal size of basis functions 
    !
    dimenmax = 0
    nsizemax = 0
    !
    !loop over J quantities
    !
    do jind = 1, nJ
        !
        ! Estimate the maximal size of the basis functions. 
        !
        dimenmax = max(dimenmax,bset_contr(jind)%Maxcontracts)
        nsizemax = max(nsizemax,maxval(bset_contr(jind)%nsize(:)))
        !
    enddo 
    !
    jmax = Jval(nJ)
    !
    nlevelI = 0
    !
    !number of initial states
    !
    nlevels = Neigenlevels 
    !
    allocate (ijterm(nJ),stat=info)
    !
    do jind=1,nJ
      !
      idimen = bset_contr(jind)%Maxcontracts
      !
      allocate (ijterm(jind)%kmat(bset_contr(jind)%Maxsymcoeffs,sym%Nrepresen),stat=info)
      call ArrayStart('ijterm',info,size(ijterm(jind)%kmat),kind(ijterm(jind)%kmat))
      !
      do igammaI = 1,sym%Nrepresen
         !
         Nterms = 0 
         !
         do iterm = 1,bset_contr(jind)%Maxsymcoeffs
           !
           ijterm(jind)%kmat(iterm,igammaI) = Nterms
           !
           Nterms = Nterms + bset_contr(jind)%irr(igammaI)%N(iterm) 
           !
         enddo
         !
      enddo
    enddo
    !
    ! loop over initial states
    !
    ener_ = 0
    !
    matsize = 0
    !
    nsizemax = 1
    nlevelsG = 0
    !
    ilist = 1 
    do while (ilist<nlevels)
      ilist = ilist + 1 
    enddo
    nlist = ilist-4
    f_t = 0
    !
    do ilevelF = 1, nlevels
       !
       jF = eigen(ilevelF)%jval
       indF = eigen(ilevelF)%jind
       ndegF   = eigen(ilevelF)%ndeg
       !
       !energy and and quanta of the final state
       !
       energyF = eigen(ilevelF)%energy
       igammaF = eigen(ilevelF)%igamma        
       quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes)
       !
       passed = .false.
       loop_ilist : do ilist=1,nJ
         !
         if (jF==analysis%j_list(ilist)) then 
           passed = .true.
           exit loop_ilist
         endif 
         !
       enddo loop_ilist
       !
       if (.not.passed) cycle
       !
       nlevelsG(indF,igammaF) = nlevelsG(indF,igammaF) + 1
       !
       nsize_need = max(bset_contr(nJ)%nsize(igammaF),1)
       !
       f_t = f_t + real(nsize_need,rk)*( real(rk,rk) )/1024_rk**3
       !
       matsize = matsize + nsize_need
       nsizemax = max(nsizemax,nsize_need)
       !
    enddo
    !
    write(my_fmt,'(a,i0,a)') "(a,i4,a,",Nrepresen,"i8)"
    !
    do jind = 1, nJ
      jI = Jval(jind) 
      write(out,my_fmt) 'Number of states for J = ',jI,' each symm = ',nlevelsG(jind,:)
      !
    enddo
    !
    dimenmax_ram = max(nsizemax,1)
    !
    write(out,"('Max size of vectors: sym = ',i0,' prim =  ',i0)") dimenmax_ram,matsize
    !
    allocate(vec_ram(nJ,sym%Nrepresen),stat = info)
    call ArrayStart('vec_ram',0,1,rk)
    !
    ! estimate the memory requirements 
    !
    if (job%verbose>=4) write(out,"(a,f12.5,a,i0,' per each vector = ',f12.5'Gb')") &
        'All vectors will require approx ',f_t,' Gb; number of vectors = ',nlevels,f_t/real(nlevels,rk)
    !
    ! memory per vector
    !
    f_t = f_t/real(nlevels,rk)
    !
    ram_size_max = max(0,min( nlevels,int((memory_limit-memory_now)/f_t,hik)))
    !
    if (job%verbose>=4) write(out,"(a,f12.5,a,f12.5,'Gb; number of vectors to be kept in RAM =  ',i0/)") &
                               'Memory needed = ',f_t*real(nlevels,rk),'Gb;  memory  available = ',&
                               memory_limit-memory_now,ram_size_max
    !
    do jind = 1, nJ
       !
       !if (job%verbose>=4) write(out,"(' J    irep  N-RAM-irep   Nreps/Ntotal   ')") 
       !
       do irep = 1,Nrepresen
          !
          f_t = real( nlevelsG(jind,irep),rk )/real(nlevels, rk )
          ram_size_ = int(f_t*ram_size_max,ik)
          !
          ram_size_ = min(ram_size_,nlevelsG(jind,irep))
          !
          if (abs(ram_size_-nlevelsG(jind,irep))<=1) ram_size_=nlevelsG(jind,irep)
          !
          ram_size(jind,irep)  = ram_size_
          !
          if (ram_size_>0) then
            !
            if (job%verbose>=5) write(out,"(i3,4x,i2,4x,i0,4x,f14.5)") jval(jind),irep,ram_size_,f_t
            !
            !dimenmax_ram_ = max(min(int(bset_contr(jind)%nsize(irep)*intensity%factor),bset_contr(jind)%nsize(irep)),1)
            !
            dimenmax_ram_ = max(bset_contr(jind)%nsize(irep),1)
            !
            matsize = int(ram_size_,hik)*int(dimenmax_ram_,hik)
            !
            if (job%verbose>=4) write(out,"('Allocate (J=',i3,',irep=',i2,') matrix   ',i8,' x',i8,' = ',i0,', ',&
                                f14.4,'Gb')") jval(jind),irep,ram_size_,dimenmax_ram_,matsize,&
                                real(matsize*8,rk)/1024.0_rk**3
            !
            allocate(vec_ram(jind,irep)%mat(dimenmax_ram_,ram_size_),stat=info)
            !
            call ArrayStart('vec_ram',info,1,rk,matsize)
            !
          endif
          !
       enddo
    enddo
    !
    allocate(istored(nlevels),stat=info)
    !
    call ArrayStart('istored',info,size(istored),kind(istored))
    !
    istored = 0
    !
    allocate(isaved(nlevels),stat=info)
    !
    call ArrayStart('isaved',info,size(isaved),kind(isaved))
    !
    isaved = 0
    !
    allocate(isave(nJ),stat=info)
    !
    if (job%verbose>=5) call MemoryReport
    !
    if (job%verbose>=6) write(out,"(/'Pre-screening, compacting and storing...')")
    !
    ! prescreen all eigenfunctions, compact and store on the disk
    !
    call TimerStart('Distributing memory')
    !
    iram = 0
    isave = 0
    ilevelsG = 0
    !
    Nstored  = 0
    !
    write(*,*) "test"
    allocate(ilevel_analyse(nlevels))
    do ilevelF = 1, nlevels
      !
      if (job%verbose>=6.and.mod(ilevelF,nlevels/min(500,ilevelF))==0 ) write(out,"('ilevel = ',i0)") ilevelF
      !
      jF = eigen(ilevelF)%jval 
      !
      !energy and and quanta of the final state
      !
      energyF = eigen(ilevelF)%energy
      igammaF = eigen(ilevelF)%igamma        
      quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
      !
      passed = .false.
      loop_ilist1 : do ilist=1,nJ
        !
        if (jF==analysis%j_list(ilist)) then 
          passed = .true.
          exit loop_ilist1
        endif 
        !
      enddo loop_ilist1
      !
      if (.not.passed) cycle
      !
      Nstored  = Nstored + 1
      !
      ilevel_analyse(Nstored) = ilevelF
      !
      indF = eigen(ilevelF)%jind
      !
      ilevelsG(indF,igammaF) = ilevelsG(indF,igammaF) + 1
      !
      ndegF   = eigen(ilevelF)%ndeg
      dimenF = bset_contr(indF)%Maxcontracts
      nsizeF = bset_contr(indF)%nsize(igammaF)
      !
      unitF = Jeigenvec_unit(indF,igammaF)
      !
      !
      if (nsizeF<=dimenmax_ram.and.iram(indF,igammaF)<ram_size(indF,igammaF)) then
        !
        iram(indF,igammaF) = iram(indF,igammaF) + 1
        !
        istored(ilevelF) = iram(indF,igammaF)
        !
        irec = eigen(ilevelF)%irec(1)
        !
        read(unitF,rec=irec) vec_ram(indF,igammaF)%mat(1:nsizeF,iram(indF,igammaF))
        !
      else
        !
        isave(indF) = isave(indF) + 1
        !
        isaved(ilevelF) = eigen(ilevelF)%irec(1)
        !
      endif
      !
    enddo
    !
    !deallocate(vecI,vecF)
    !call ArrayStop('intensity-vec')
    !
    call TimerStop('Distributing memory')
    !
    write(out,"(/'...done!')")
    !
    dimenmax_ = 1
    do jind = 1, nJ
       do irep = 1,Nrepresen
         dimenmax_ = max(bset_contr(jind)%nsize(irep),dimenmax_)
       enddo
    enddo
    !
    allocate(vecI(dimenmax_),vecF(dimenmax_),vec(dimenmax),vec_(dimenmax),vecE(dimenmax_),vecE_(dimenmax), stat = info)
    call ArrayStart('density-vectors',info,size(vecI),kind(vecI))
    call ArrayStart('density-vectors',info,size(vecF),kind(vecF))
    call ArrayStart('density-vectors',info,size(vec),kind(vec))
    call ArrayStart('density-vectors',info,size(vec_),kind(vec_))
    call ArrayStart('density-vectors',info,size(vecE),kind(vecE))
    call ArrayStart('density-vectors',info,size(vecE_),kind(vecE_))
    !
    allocate(icoeffI(dimenmax), stat = info)
    call ArrayStart('density-vectors',info,size(icoeffI),kind(icoeffI))
    !
    !write(out,"('Number of states for each symm = ',<Nrepresen>i8)") nlevelsG(:)
    !
    write(my_fmt,'(a,i0,a)') "(a,i4,a,",Nrepresen,"i8)"
    !
    do jind = 1, nJ
      jI = Jval(jind) 
      write(out,my_fmt) 'Number of states for ',jI,' each symm = ',nlevelsG(jind,:)
    enddo
    !
    if (all(nlevelsG(:,:)<1)) then 
      !
      write(out,"('DM_contribution_symmvec error: number of states to analyse is zero!')") 
      return 
      !
    endif
    !
    !  The matrix where some of the eigenvectors will be stored
    !
    !
    if (job%verbose>=5) call MemoryReport
    !
    if (analysis%population) then
       !
       ! prepare the external vector
       !
       passed = .false.
       !
       loop_Jref: do indI = 1, nJ
         !
         jI = Jval(indI) 
         !
         if (jI==analysis%ref%J) then 
           passed = .true.
           exit loop_Jref
         endif 
         !
       enddo loop_Jref
       !
       if (.not.passed) then
         write(out,"('DM_contribution_symmvec: The reference state is not in the analyss list',i4)") analysis%ref%J
         stop 'DM_contribution_symmvec: The reference state is not in the analyss list'
       endif
       !
       igammaI = analysis%ref%iGamma
       ilevelI = analysis%ref%iLevel
       !
       external_vec_unit = TReigenvec_unit_external(indI,Jval,igammaI)
       !
       irec = eigen(ilevelI)%irec(1)
       idegI = 1
       nsizeI = bset_contr(indI)%nsize(igammaI)
       dimenI = bset_contr(indI)%Maxcontracts
       read(external_vec_unit,rec=irec) vecE(1:nsizeI)
       call convert_symvector_to_contrvector(1,dimenI,igammaI,idegI,ijterm(indI)%kmat(:,:),vecE,vecE_)
       !
       do jind = 1, nJ 
         indI = eigen(ilevel_analyse(1))%jind
         dimenI = bset_contr(eigen(ilevel_analyse(1))%jind)%Maxcontracts
         !
         levels_loop: do i =1,Nstored
           !
           ilevelI = ilevel_analyse(i)
           !
           ! start measuring time per line
           !
           indI = eigen(ilevelI)%jind
           ! 
           !dimension of the bases for the initial states
           !
           dimenI = bset_contr(indI)%Maxcontracts
           !
           !energy, quanta, and gedeneracy order of the initial state
           !
           jI = eigen(ilevelI)%jval
           !
           if(jI/=Jval(jind)) then
             cycle
           endif
           !
           energyI = eigen(ilevelI)%energy
           igammaI  = eigen(ilevelI)%igamma
           quantaI(0:nmodes) = eigen(ilevelI)%quanta(0:nmodes)
           normalI(0:nmodes) = eigen(ilevelI)%normal(0:nmodes)
           ndegI   = eigen(ilevelI)%ndeg
           nsizeI = bset_contr(indI)%nsize(igammaI)
           !
           ! where the eigenvector is stored 
           !
           unitI = Jeigenvec_unit(indI,igammaI) 
           !
           !read eigenvector of initial state
           !
           if (istored(ilevelI)==0) then
               !
               call TimerStart('Reading eigenfuncs')
               !
               irec = isaved(ilevelI)
               !
               read(unitI,rec=irec) vecI(1:nsizeI)
               !
               call TimerStop('Reading eigenfuncs')
               !
           else
               !
               iram_ = istored(ilevelI)
               !
               vecI(1:nsizeI) = vec_ram(indI,igammaI)%mat(1:nsizeI,iram_)
               !
           endif
           !
           ndegF  = ndegI
           !
           do idegI = 1, ndegI
              !
              call convert_symvector_to_contrvector(indI,dimenI,igammaI,idegI,ijterm(indI)%kmat(:,:),vecI,vec)
              !
              population = sum(vec(1:dimenI)*vecE_(1:dimenI))
              !
              write(my_fmt1,'(a,i0,a)') "(i7,2x,a4,i4,f13.6,1x,a1,",nmodes,"(i3),a1,1x,e16.8)"
              !
              write(out,my_fmt1) Jval(jind),sym%label(igammaI),ilevelI,energyI,"(",quantaI(1:nmodes),")",population
              !
           enddo
           !
         end do levels_loop
         !
       enddo
       !
    endif
    !
    !
    if (job%verbose>=5) call MemoryReport
    !
    if (analysis%reducible_eigen_contribution) then
       !
       ioname = "highest eigen contributions" 
       call IOstart(trim(ioname),iounit)
       open(unit = iounit, action = 'write',status='replace',file = "contribution.chk")
       do jind = 1, nJ 
       indI = eigen(ilevel_analyse(1))%jind
       dimenI = bset_contr(eigen(ilevel_analyse(1))%jind)%Maxcontracts
       allocate(highest_contribution(bset_contr(indI)%icontr2icase(dimenI,1)))
       allocate(ener_top_contri(bset_contr(indI)%icontr2icase(dimenI,1)))
       highest_contribution = 0.0_ark
       allocate(basis_set_contraction(bset_contr(indI)%icontr2icase(dimenI,1),0:Nclasses + 1))
       ipp = 0 
       !
       Ilevels_loop: do i =1,Nstored
         !
         ilevelI = ilevel_analyse(i)
         !
         ! start measuring time per line
         !
         indI = eigen(ilevelI)%jind
         ! 
         !dimension of the bases for the initial states
         !
         dimenI = bset_contr(indI)%Maxcontracts
         !
         !energy, quanta, and gedeneracy order of the initial state
         !
         jI = eigen(ilevelI)%jval
         !
         if(jI/=Jval(jind)) then
           cycle
         endif
         energyI = eigen(ilevelI)%energy
         igammaI  = eigen(ilevelI)%igamma
         quantaI(0:nmodes) = eigen(ilevelI)%quanta(0:nmodes)
         normalI(0:nmodes) = eigen(ilevelI)%normal(0:nmodes)
         ndegI   = eigen(ilevelI)%ndeg
         nsizeI = bset_contr(indI)%nsize(igammaI)
         !
         ! where the eigenvector is stored 
         !
         unitI = Jeigenvec_unit(indI,igammaI) 
         !
         !read eigenvector of initial state
         !
         !
         if (istored(ilevelI)==0) then
             !
             call TimerStart('Reading eigenfuncs')
             !
             irec = isaved(ilevelI)
             !
             read(unitI,rec=irec) vecI(1:nsizeI)
             !
             call TimerStop('Reading eigenfuncs')
             !
         else
             !
             iram_ = istored(ilevelI)
             !
             vecI(1:nsizeI) = vec_ram(indI,igammaI)%mat(1:nsizeI,iram_)
             !
         endif
         !
         ndegF  = ndegI
         !
         do idegI = 1, ndegI
            !
            call convert_symvector_to_contrvector(indI,dimenI,igammaI,idegI,ijterm(indI)%kmat(:,:),vecI,vec)
            !
            write(my_fmt1,'(a,i0,a)') "(2x,i4,i7,e16.8,3x,a1,3i3,1x,i4,i2,a1,1x,a1,",Nclasses,"(i3),a1)"
            !
            do irootI = 1, dimenI
                 !
                 irow = bset_contr(indI)%icontr2icase(irootI,1)
                 ib   = bset_contr(indI)%icontr2icase(irootI,2)
                 !
                 cnu_i(0:Nclasses) = bset_contr(indI)%contractive_space(0:Nclasses, irow)
                 !
                 ktau = bset_contr(indI)%ktau(irootI)
                 tauI  = mod(ktau,2_ik)
                 kI = bset_contr(indI)%k(irootI)
                 !
                 !ndeg = bset_contr(jind)%index_deg(irow)%size1
                 !
                 if(abs(vec(irootI))>highest_contribution(irow)) then 
                    highest_contribution(irow) = abs(vec(irootI))
                    ener_top_contri(irow) = energyI
                 endif
                 basis_set_contraction(irow,0:Nclasses) = cnu_i(0:Nclasses)
                 basis_set_contraction(irow,Nclasses+1) = jI
                 if (.false.) then  
                   !
                   write(out,my_fmt1) & 
                              igammaI,irootI,&
                              vec(irootI),"(", &
                              jI,kI,tauI,irow,ib,")", &
                              "(",cnu_i(1:Nclasses),")"
                   !
                 endif
                 !
            end do
            !
         enddo
         !
       end do Ilevels_loop
       !
       write(iounit,"(a)") "Start-contract"
       !
       write(my_fmt1,'(a,i0,a)') "(2x,i7,e16.8,3x,f13.6,3x,i3,1x,",Nclasses,"(1x,i5))"
       !
       do i = 1, size(highest_contribution)
         write(iounit,my_fmt1) i, highest_contribution(i), ener_top_contri(i) -intensity%ZPE,  &
                             basis_set_contraction(i, Nclasses+1), &
                              basis_set_contraction(i,1:Nclasses) 
       enddo
       !
       write(iounit,"(a)") "End-contract"
       !
       deallocate(highest_contribution)
       deallocate(basis_set_contraction)
       enddo
       close(iounit)
       !
    endif
    !
    deallocate(vecI, vecF, vec, vec_,icoeffI,vecE,vecE_)
    call ArrayStop('density-vectors')
    !
    do irep = 1,Nrepresen
     do jind = 1,nJ
      !
      nullify(vec_ram(jind,irep)%mat)
      !
      if (associated(vec_ram(jind,irep)%mat)) deallocate(vec_ram(jind,irep)%mat)
      !
     enddo
    enddo
    !
    call ArrayStop('vec_ram')
    !
    do jind = 1,nJ
      !
      if (associated(ijterm(jind)%kmat)) deallocate(ijterm(jind)%kmat)
      !
    enddo
    call ArrayStop('ijterm')
    !
    deallocate(nlevelsG,ilevelsG,ram_size,iram)
    call ArrayStop('nlevelsG')
    call ArrayStop('ilevelsG')
    call ArrayStop('ram_size')
    call ArrayStop('iram')
    !
    deallocate(istored)
    call ArrayStop('istored')
    !
    deallocate(isaved)
    call ArrayStop('isaved')
    !
    call TimerStop('Contribution analysis')
    !
  end subroutine DM_contribution_symmvec



 !
 ! Electric dipole moment intensities and transition moments calculations for vectors in symmetrised representation
 !
 subroutine dm_intensity_symmvec(Jval)
    implicit none
    !
    integer(ik),intent(in)  :: Jval(:)
    !
    integer(ik)    :: nJ,dimenmax
    integer(ik)    :: ilevelI, ilevelF, ndegI, ndegF, idegI, idegF, irec, idimen, nsizeF,nsizeI,ilevelI_
    integer(ik)    :: cdimenI,irep
    integer(ik)    :: nmodes,info,indI,indF,itransit,Ntransit,jI,jF,Nrepresen,Ntransit_
    integer(ik)    :: igammaI,igammaF,quantaI(0:molec%nmodes),quantaF(0:molec%nmodes),normalF(0:molec%nmodes),&
                      normalI(0:molec%nmodes)
    integer(ik)    :: dimenI,dimenF,unitI,unitF,irootF,irootI,imu,jmax,kI,kF,icontrF,irow,ib,int_increm,ktauI
    real(rk)       :: energyI, energyF, nu_if,linestr,linestr_deg(sym%Maxdegen,sym%Maxdegen),f_t,ener_
    real(rk)       :: time_per_line,real_time,total_time_predict,cpu_time
    real(rk)       :: tm_deg(3,sym%Maxdegen,sym%Maxdegen)
    logical        :: passed,passed_
    logical        :: vector_diagonal = .false.
    integer(ik)    :: ID_I
    !
    integer(ik), allocatable :: icoeffI(:), istored(:),isave(:),isaved(:),iID(:)
    integer(ik), allocatable :: nlevelsG(:,:),ram_size(:,:),ilevelsG(:,:),iram(:,:)
    real(rk),    allocatable :: vecI(:), vecF(:),vec(:),vecPack(:),vec_(:)
    real(rk),    allocatable :: vecPack_kblock(:,:)
    integer(ik), allocatable :: icoeff_kblock(:,:),cdimen_kblock(:),itau_kblock(:,:)
    real(rk),allocatable     :: half_linestr(:,:,:,:),threej(:,:,:,:),max_intens_state(:)
    integer(ik), pointer :: Jeigenvec_unit(:,:)
    type(DmatT),pointer  :: vec_ram(:,:)
    type(DkmatT),pointer :: ijterm(:)
    !
    integer(ik)  :: jind,nlevels,igamma,nformat,nclasses,nsize_need,nlevelI,dimenmax_,nsizemax
    !
    integer(ik)  :: ram_size_,alloc_p
    !
    integer(ik)  :: iram_,Nterms,iterm,ielem,isrootF,isrootI,nelem
    !
    integer(ik)  :: igamma_pair(sym%Nrepresen)
    integer(hik) :: matsize,ram_size_max,dimenmax_ram,dimenmax_ram_
    !
    !
    character(len=1) :: branch
    character(len=wl) :: my_fmt,my_fmt_tm,my_fmt1
    !
    character(cl)           :: filename, ioname
    character(4)            :: jchar,gchar
    integer(ik)             :: iounit
    !
    real(rk)     :: ddot,boltz_fc,beta,intens_cm_mol,A_coef_s_1,A_einst,absorption_int,dtemp0,dmu(3),&
                    intens_cm_molecule
    !
    character(len=cl) :: my_fmt0 !format for I/O specification
    !
    call TimerStart('Intensity calculations')
    !
    if (sym%maxdegen>2) then 
      !
      write(out,"('dm_intensity_symmvec: this procedure has not been tested for the symmetries with degeneracies higher than 2')")
      !stop 'dm_intensity_symmvec was not tested for present symmetry'
      !
    endif
    !
    ! In case the vector obtained without diagonalization, as diagonal solution of the Hamiltonian matrix
    ! the reading of the vectors can be omitted 
    !
    if (trim(job%diagonalizer)=='NO-DIAGONALIZATION') vector_diagonal = .true.
    !
    nmodes = molec%nmodes
    nclasses = size(eigen(1)%cgamma)-1
    Nrepresen = sym%Nrepresen
    !
    beta = planck * vellgt / (boltz * intensity%temperature)
    intens_cm_mol  = 8.0d-36 * pi**3* avogno / (3.0_rk * planck * vellgt)
    intens_cm_molecule  = 8.0d-36 * pi**3/ (3.0_rk * planck * vellgt)
    A_coef_s_1     =64.0d-36 * pi**4  / (3.0_rk * planck)
    !
    nJ = size(Jval)
    !
    ! Prepare the list of units with the stored eigenvectors
    !
    allocate(Jeigenvec_unit(nJ,sym%Nrepresen), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: Jeigenvec_unit - out of memory'
    !
    allocate(nlevelsG(nJ,sym%Nrepresen),ilevelsG(nJ,sym%Nrepresen),ram_size(nJ,sym%Nrepresen),&
             iram(nJ,sym%Nrepresen),stat = info)
    call ArrayStart('nlevelsG',info,size(nlevelsG),kind(nlevelsG))
    call ArrayStart('ilevelsG',info,size(ilevelsG),kind(ilevelsG))
    call ArrayStart('ram_size',info,size(ram_size),kind(ram_size))
    call ArrayStart('iram',info,size(iram),kind(iram))
    !
    do jind=1,nJ
      do igamma = 1,sym%Nrepresen
        Jeigenvec_unit(jind,igamma) = TReigenvec_unit(jind,Jval,igamma)
      enddo
    enddo
    !
    ! maximal size of basis functions 
    !
    dimenmax = 0
    nsizemax = 0
    !
    !loop over J quantities
    !
    do jind = 1, nJ
        !
        ! Estimate the maximal size of the basis functions. 
        !
        dimenmax = max(dimenmax,bset_contr(jind)%Maxcontracts)
        nsizemax = max(nsizemax,maxval(bset_contr(jind)%nsize(:)))
        !
    enddo 
    !
    jmax = Jval(nJ)
    !
    ! Precompute the 3-j symbol values 
    !
    allocate(threej(0:jmax,0:jmax,-1:1,-1:1), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: threej - out of memory'
    !
    threej = 0 
    !
    do jI = 0,jmax
      do jF = max(jI-1,0),min(jI+1,jmax)
        do kI = 0,jI
          do kF = max(kI-1,0),min(kI+1,jF)
            !
            !threej(jI, kI, jF - jI, kF - kI) = cg(jI, kI, jF - jI, kF - kI) *                                                 &
            !               (-1.0_rk) ** (jI - 1 + kF) / sqrt( real(2*jF + 1,rk) )
            !
            !dp = three_j(jI, 1, jF, kI, kF - kI, -kF)
            !
            !if (abs(dp-threej(jI, kI, jF - jI, kF - kI))>sqrt(small_)) then 
            !  write(out,"('a problem with threej for jI,jF,kI,kF = ',4i4,3es16.8)") jI,jF,kI,kF,dp,threej(jI, kI, jF - jI, kF - kI),dp-threej(jI, kI, jF - jI, kF - kI)
            !  stop 'a problem with threej'
            !endif
            !
            threej(jI, kI, jF - jI, kF - kI) = three_j(jI, 1, jF, kI, kF - kI, -kF)
            !
          enddo
        enddo
      enddo
    enddo
    !
    ! in the case of the proper rotatonal symmetry the matrix elements of the 
    ! wigner matrix have to be computed using the correspoding symmetrization 
    ! transformation derived here. they will be used together with the 
    ! vibrational matrix elements of the dipole moment when evaluating 
    ! line strengthes using eigenvectos.  
    !
    if (job%rotsym_do) call contraced_rotational_dipole(nJ, jval, jmax, threej)
    !
    !allocate arrays for eigenvectors and coefficients
    !
    ! First we count transitions to be calculated. It will help us to keep track of the 
    ! calculation progress.
    !
    Ntransit = 0 
    !
    nlevelI = 0
    !
    !number of initial states
    !
    nlevels = Neigenlevels
    !
    allocate (ijterm(nJ),stat=info)
    !
    do jind=1,nJ
      !
      idimen = bset_contr(jind)%Maxcontracts
      !
      allocate (ijterm(jind)%kmat(bset_contr(jind)%Maxsymcoeffs,sym%Nrepresen),stat=info)
      call ArrayStart('ijterm',info,size(ijterm(jind)%kmat),kind(ijterm(jind)%kmat))
      !
      do igammaI = 1,sym%Nrepresen
         !
         Nterms = 0 
         !
         do iterm = 1,bset_contr(jind)%Maxsymcoeffs
           !
           ijterm(jind)%kmat(iterm,igammaI) = Nterms
           !
           Nterms = Nterms + bset_contr(jind)%irr(igammaI)%N(iterm) 
           !
         enddo
         !
      enddo
    enddo
    !
    ! For a given symmetry igamma with some gns(igamma) we find its counterpart jgamma/=igamma
    ! having the same gns(jgamma). We assume that there is only one such pair 
    ! in case of absorption or emission calcs. 
    !
    call find_igamma_pair(igamma_pair)
    !
    call TimerStart('Intens_Filter-1')
    !
    if (intensity%int_increm<1e6) write(out,"(/'Intensity increments:')")
    !
    ! loop over initial states
    !
    ener_ = 0
    !
    !omp parallel do private(ilevelI,jI,energyI,igammaI,quantaI,ilevelF,jF,energyF,igammaF,quantaF,passed) schedule(guided) reduction(+:Ntransit,nlevelI)
    do ilevelI = 1, nlevels
      !
      ! rotational quantum number 
      !
      jI = eigen(ilevelI)%jval
      !
      !energy energy and and quanta of the initial state
      !
      energyI = eigen(ilevelI)%energy
      !
      igammaI  = eigen(ilevelI)%igamma
      quantaI(0:nmodes) = eigen(ilevelI)%quanta(0:nmodes) 
      !
      call energy_filter_lower(jI,energyI,quantaI,eigen(ilevelI)%normal(0),passed)
      !
      if (.not.passed) cycle
      !
      nlevelI = nlevelI + 1
      !
      !loop over final states
      !
      do ilevelF = 1, nlevels
         !
         jF = eigen(ilevelF)%jval 
         !
         !energy and and quanta of the final state
         !
         energyF = eigen(ilevelF)%energy
         igammaF = eigen(ilevelF)%igamma        
         quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
         !
         call intens_filter(jI,jF,energyI,energyF,igammaI,igammaF,quantaI,quantaF,igamma_pair,passed)
         !
         if (passed) then 
           !
           Ntransit = Ntransit + 1
           !
           if (mod(Ntransit,intensity%int_increm)==0.and.energyI>ener_) then
              ener_ = energyI
              write(out,"(i0,f16.4)") Ntransit,energyI-intensity%ZPE
           endif
           !
         endif 
         !
      enddo
      !
    enddo
    !omp end parallel do
    !
    call TimerStop('Intens_Filter-1')
    !
    !loop over final states -> count states for each symmetry
    !
    dimenmax_ = 1
    do jind = 1, nJ
       do irep = 1,Nrepresen
         dimenmax_ = max(bset_contr(jind)%nsize(irep),dimenmax_)
       enddo
    enddo
    !
    allocate(vecI(dimenmax_),vecF(dimenmax_), stat = info)
    call ArrayStart('intensity-vec',info,size(vecI),kind(vecI))
    call ArrayStart('intensity-vec',info,size(vecF),kind(vecF))
    !
    nlevelsG = 0
    !
    f_t = 0
    matsize = 0
    !
    intensity%factor = 1
    nsizemax = 1
    !
    do ilevelF = 1, nlevels
       !
       jF = eigen(ilevelF)%jval
       indF = eigen(ilevelF)%jind
       ndegF   = eigen(ilevelF)%ndeg
       !
       !energy and and quanta of the final state
       !
       energyF = eigen(ilevelF)%energy
       igammaF = eigen(ilevelF)%igamma        
       quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes)
       normalF(0:nmodes) = eigen(ilevelF)%normal(0:nmodes)
       !
       call energy_filter_upper(jF,energyF,quantaF,normalF(0),passed)
       !
       call energy_filter_lower(jF,energyF,quantaF,normalF(0),passed_)
       !
       if (.not.passed.and..not.passed_) cycle
       !
       nlevelsG(indF,igammaF) = nlevelsG(indF,igammaF) + 1
       !
       nsize_need = max(bset_contr(nJ)%nsize(igammaF),1)
       !
       f_t = f_t + real(nsize_need,rk)*( real(rk,rk) )/1024_rk**3
       !
       !nsize_ram = nsize_need + nsize_need
       !
       matsize = matsize + nsize_need
       nsizemax = max(nsizemax,nsize_need)
       !
    enddo
    !
    write(my_fmt0,'(a,i0,a)') "(a,i4,a,",Nrepresen,"i8)"
    !
    do jind = 1, nJ
      jI = Jval(jind) 
      write(out,my_fmt0) 'Number of states for J = ',jI,' each symm = ',nlevelsG(jind,:)
      !
    enddo
    !
    dimenmax_ram = max(nsizemax,1)
    !
    write(out,"('Max size of vectors: sym = ',i0,' prim =  ',i0)") dimenmax_ram,matsize
    !
    allocate(vec_ram(nJ,sym%Nrepresen),stat = info)
    call ArrayStart('vec_ram',0,1,rk)
    !
    ! estimate the memory requirements 
    !
    if (job%verbose>=4) write(out,"(a,f12.5,a,i0,' per each vector = ',f12.5'Gb')") &
                             'All vectors will require approx ',f_t,' Gb; number of vectors = ',&
                             nlevels,f_t/real(nlevels,rk)
    !
    ! memory per vector
    !
    f_t = f_t/real(nlevels,rk)
    !
    ram_size_max = max(0,min( nlevels,int((memory_limit-memory_now*1.2_rk)/f_t,hik)))
    !
    if (job%verbose>=4) write(out,"('Memory needed = ',f12.5,'Gb;  memory  available = ',f12.5,&
                        'Gb; number of vectors to be kept in RAM =  ',i0/)") f_t*real(nlevels,rk),&
                        memory_limit-memory_now,ram_size_max
    !
    do jind = 1, nJ
       !
       !if (job%verbose>=4) write(out,"(' J    irep  N-RAM-irep   Nreps/Ntotal   ')") 
       !
       do irep = 1,Nrepresen
          !
          f_t = real( nlevelsG(jind,irep),rk )/real(nlevels, rk )
          ram_size_ = int(f_t*ram_size_max,ik)
          !
          ram_size_ = min(ram_size_,nlevelsG(jind,irep))
          !
          if (abs(ram_size_-nlevelsG(jind,irep))<=1) ram_size_=nlevelsG(jind,irep)
          !
          ram_size(jind,irep)  = ram_size_
          !
          if (ram_size_>0) then
            !
            if (job%verbose>=5) write(out,"(i3,4x,i2,4x,i0,4x,f14.5)") jval(jind),irep,ram_size_,f_t
            !
            !dimenmax_ram_ = max(min(int(bset_contr(jind)%nsize(irep)*intensity%factor),bset_contr(jind)%nsize(irep)),1)
            !
            dimenmax_ram_ = max(bset_contr(jind)%nsize(irep),1)
            !
            matsize = int(ram_size_,hik)*int(dimenmax_ram_,hik)
            !
            if (job%verbose>=4) write(out,"('Allocate (J=',i3,',irep=',i2,') matrix   ',i8,' x',i8,' = ',i0,', ',f14.4,'Gb')") &
                                jval(jind),irep,ram_size_,dimenmax_ram_,matsize,real(matsize*8,rk)/1024.0_rk**3
            !
            allocate(vec_ram(jind,irep)%mat(dimenmax_ram_,ram_size_),stat=info)
            !
            call ArrayStart('vec_ram',info,1,rk,matsize)
            !
          endif
          !
       enddo
    enddo
    !
    allocate(istored(nlevels),stat=info)
    !
    call ArrayStart('istored',info,size(istored),kind(istored))
    !
    istored = 0
    !
    allocate(isaved(nlevels),stat=info)
    !
    call ArrayStart('isaved',info,size(isaved),kind(isaved))
    !
    isaved = 0
    !
    allocate(isave(nJ),stat=info)
    !
    if (job%verbose>=5) call MemoryReport
    !
    if (job%verbose>=4) write(out,"(/'Pre-screening, compacting and storing...')")
    !
    ! prescreen all eigenfunctions, compact and store on the disk
    !
    call TimerStart('Distributing memory')
    !
    iram = 0
    isave = 0
    ilevelsG = 0
    !
	!cdimenmax = 0
	!cfactor = 0
	!
	!dimenmax_ram = max(min(int(matsize),nsizemax),1)
    !
    !
    if (Ntransit==0) then 
         write(out,"('dm_intensity_symmvec: the transition filters are too tight: no entry')") 
         write(out,"(' window = ',2f12.4,'; Elow = ',2f12.4,'; Eupp = ',2f12.4,' cm-1')") &
              intensity%freq_window(1:2),intensity%erange_low(1:2),intensity%erange_upp(1:2)
         stop 'dm_intensity_symmvec: the filters are too tight' 
    endif 
    !
    !omp do private(ilevelF,jF,energyF,igammaF,quantaF,indF,ndegF,dimenF,unitF,unitO,unitC,idegF,irec,cdimen,idimen) schedule(guided)
    do ilevelF = 1, nlevels
      !
      if (job%verbose>=5.and.mod(ilevelF,nlevels/min(500,ilevelF))==0 ) write(out,"('ilevel = ',i0)") ilevelF
      !
      jF = eigen(ilevelF)%jval 
      !
      !energy and and quanta of the final state
      !
      energyF = eigen(ilevelF)%energy
      igammaF = eigen(ilevelF)%igamma        
      quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
      normalF(0:nmodes) = eigen(ilevelF)%normal(0:nmodes)
      !
      call energy_filter_upper(jF,energyF,quantaF,normalF(0),passed)
      !
      call energy_filter_lower(jF,energyF,quantaF,normalF(0),passed_)
      !
      if (.not.passed.and..not.passed_) cycle
      !
      indF = eigen(ilevelF)%jind
      !
      ilevelsG(indF,igammaF) = ilevelsG(indF,igammaF) + 1
      !
      ndegF   = eigen(ilevelF)%ndeg
      dimenF = bset_contr(indF)%Maxcontracts
      nsizeF = bset_contr(indF)%nsize(igammaF)
      !
      unitF = Jeigenvec_unit(indF,igammaF)
      !
      if (unitF==-1) stop 'This file is not supposed to be accessed'
      !
      indI = eigen(ilevelF)%icoeff
      !
      if (nsizeF<=dimenmax_ram.and.iram(indF,igammaF)<ram_size(indF,igammaF)) then
        !
        iram(indF,igammaF) = iram(indF,igammaF) + 1
        !
        istored(ilevelF) = iram(indF,igammaF)
        !
        irec = eigen(ilevelF)%irec(1)
        !
        if (.not.vector_diagonal) then
           !
           read(unitF,rec=irec) vec_ram(indF,igammaF)%mat(1:nsizeF,iram(indF,igammaF))
           !
        else
           !
           vec_ram(indF,igammaF)%mat(:,iram(indF,igammaF)) = 0 
           vec_ram(indF,igammaF)%mat(eigen(ilevelF)%icoeff,iram(indF,igammaF)) = 1.0_rk
           !
        endif
        !
      else
        !
        isave(indF) = isave(indF) + 1
        !
        isaved(ilevelF) = eigen(ilevelF)%irec(1)
        !
      endif
      !
    enddo
    !
    deallocate(vecI,vecF)
    call ArrayStop('intensity-vec')
    !
    call TimerStop('Distributing memory')
    !
    write(out,"(/'...done!')")
    !
    dimenmax_ = 1
    do jind = 1, nJ
       do irep = 1,Nrepresen
         dimenmax_ = max(bset_contr(jind)%nsize(irep),dimenmax_)
       enddo
    enddo
    !
    allocate(vecI(dimenmax_),vec(dimenmax),vecPack(dimenmax), stat = info)
    !
    call ArrayStart('intensity-vectors',info,size(vecI),kind(vecI))
    !call ArrayStart('intensity-vectors',info,size(vecF),kind(vecF))
    call ArrayStart('intensity-vectors',info,size(vec),kind(vec))
    call ArrayStart('intensity-vectors',info,size(vecPack),kind(vecpack))
    !
    allocate(icoeffI(dimenmax), stat = info)
    call ArrayStart('intensity-vectors',info,size(icoeffI),kind(icoeffI))

    !
    if (.not.job%rotsym_do) then 
      allocate(icoeff_kblock(dimenmax,0:jmax),cdimen_kblock(0:jmax),itau_kblock(dimenmax,0:jmax),stat = info)
      call ArrayStart('intensity-vectors-kblock',info,size(icoeff_kblock),kind(icoeff_kblock))
      call ArrayStart('intensity-vectors-kblock',info,size(cdimen_kblock),kind(cdimen_kblock))
      call ArrayStart('intensity-vectors-kblock',info,size(itau_kblock),kind(itau_kblock))
      allocate(vecPack_kblock(dimenmax,0:jmax), stat = info)
      call ArrayStart('intensity-vectors-kblock',info,size(vecPack_kblock),kind(vecPack_kblock))
    endif

    !
    ! loop over final states -> count states for each symmetry
    !
    nlevelsG = 0
    !
    do ilevelF = 1, nlevels
       !
       jF = eigen(ilevelF)%jval 
       !
       !energy and and quanta of the final state
       !
       energyF = eigen(ilevelF)%energy
       igammaF = eigen(ilevelF)%igamma        
       quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes)
       normalF(0:nmodes) = eigen(ilevelF)%normal(0:nmodes)
       !
       indF = eigen(ilevelF)%jind
       !
       call energy_filter_upper(jF,energyF,quantaF,normalF(0),passed)
       !
       if (.not.passed) cycle
       !
       nlevelsG(indF,igammaF) = nlevelsG(indF,igammaF) + 1
       !
    enddo
    !
    !write(out,"('Number of states for each symm = ',<Nrepresen>i8)") nlevelsG(:)
    !
    write(my_fmt0,'(a,i0,a)') "(a,i4,a,",Nrepresen,"i8)"
    !
    do jind = 1, nJ
      jI = Jval(jind) 
      write(out,my_fmt0) 'Number of states for ',jI,' each symm = ',nlevelsG(jind,:)
      !
    enddo
    !
#if (dipole_debug >= 0)
       write(out,"(' Total number of lower states = ',i8)") nlevelI
       write(out,"(' Total number of transitions  = ',i10)") Ntransit
#endif
    !
    !
    ! In order to speed up the line strength evaluation, 
    ! S_if = | <i|mu|f> |^2  = | \sum_{nm} C_in C_fm <n|mu|m> |^2
    ! in three steps:
    ! 1. Evaluating the expansion of the i-state 
    !    s_{im} =  \sum_{n} C_in  <n|mu|m> 
    ! 2. performing the expansion of the final state f:
    !    s_{if} =  \sum_{m} C_fm  s_{im}. 
    ! 3. Building S_{if}
    !    S_{if} = s_{if}^2
    !
    !  The temporaly object s_{im} will be referted to as 
    !  a half-linestrength "half_linestr"
    !
    allocate(half_linestr(dimenmax,nJ,sym%Maxdegen,3),stat=info)
    !
    call ArrayStart('half_linestr',info,size(half_linestr),kind(half_linestr))
    !
    ! in case of the TM-pruning the maximal intensity of each state will be estimated and stored
    !
    if (intensity%pruning) then
      !
      allocate(max_intens_state(nlevels),stat=info)
      call ArrayStart('max_intens_state',info,size(max_intens_state),kind(max_intens_state))
      !
      max_intens_state = 0
      !
      do indI = 1, nJ
         !
         do igammaI = 1,sym%Nrepresen
           !
           write(jchar, '(i4)') jval(indI)
           write(gchar, '(i3)') igammaI
           !
           filename = trim(job%eigenfile%filebase)//'_intens'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
           !
           ioname = 'max state intensities'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))
           !
           call IOstart(trim(ioname),iounit)
           !
           open(unit = iounit, action = 'write',status='replace',file = filename)
           !
         enddo
         !
      enddo
      !
    else 
      ! 
      allocate(max_intens_state(1),stat=info)
      !
    endif
    !
    ! We can construct a unique ID for exomol format for all states using the formula 
    ! ID = ilevel(J) + Nvib*(2(J-1)+1), 
    ! where Nvib = bset_contr(1)%Maxcontracts is the size of the vibrational basis and the number of states at J=0.
    ! Thus we only need to introduce Nstates_before(J) = Nvib*(2(J-1)+1)
    !
    !  - prepare the states file:
    !
    ! loop over initial states for current J-values for the exomol-format only
    !
    if (job%exomol_format) then
      !
      !write(out,"(/a/)") 'States file in the Exomol format'
      !
      allocate(iID(nlevels),stat=info)
      call ArrayStart('iID',info,size(iID),kind(iID))
      !
      do ilevelI = 1,nlevels
         !
         indI = eigen(ilevelI)%jind
         igammaI  = eigen(ilevelI)%igamma
         !
         !if (indI/=jind.or.igamma/=igammaI) cycle
         !
         !dimenI = bset_contr(indI)%Maxcontracts
         !
         !energy, quanta, and gedeneracy order of the initial state
         !
         !jI = eigen(ilevelI)%jval
         !energyI = eigen(ilevelI)%energy-intensity%ZPE
         !quantaI(0:nmodes) = eigen(ilevelI)%quanta(0:nmodes)
         !normalI(0:nmodes) = eigen(ilevelI)%normal(0:nmodes)
         !
         ID_I = eigen(ilevelI)%ilevel + bset_contr(indI)%nsize_base(igammaI)
         !
         iID(ilevelI) = ID_I
         !
         !write(out,"(i12,1x,f12.6,1x,i6,1x,i7,2x,a3,2x,<nmodes>i3,1x,<nclasses>(1x,a3),1x,2i4,1x,a3,2x,f5.2,' ::',1x,i9,1x,<nmodes>i3)") & 
         !ID_I,energyI,int(intensity%gns(igammaI),4)*(2*jI+1),jI,sym%label(igammaI),normalI(1:nmodes),eigen(ilevelI)%cgamma(1:nclasses),&
         !eigen(ilevelI)%krot,eigen(ilevelI)%taurot,eigen(ilevelI)%cgamma(0),&
         !eigen(ilevelI)%largest_coeff,eigen(ilevelI)%ilevel,quantaI(1:nmodes)
         !
      enddo
      !
    endif
    !
    !
    !  The matrix where some of the eigenvectors will be stored
    !
    if (job%verbose>=5) call MemoryReport
    !
    write(out,"(/a,a,a,a)") 'Linestrength S(f<-i) [Debye**2],',' Transition moments [Debye],'& 
                          &,'Einstein coefficient A(if) [1/s],','and Intensities [cm/mol]'
    !
    ! prepare the format to write intensities
    !
    !
    nformat = sym%Maxdegen**2
    !
    write(my_fmt,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') '(i4,1x,a4,3x,"<-",i4,1x,a4,3x,a1,2x,f11.4,1x,"<-",1x,f11.4,1x,&
                   f11.4,2x,"(",1x,a3,x,i3,1x,")",1x,"(",1x,',&
                   nclasses,'(x,a3),1x,',nmodes,'(1x, i3),1x,")",1x,"<- ","(",1x,a3,x,i3,1x,")",1x,"(",1x,',&
                   nclasses,'(x,a3),1x,',nmodes,'(1x, i3),1x,")",1x,3(1x, es16.8),2x,1x,i6,1x,"<-",&
                   &1x,i6,1x,i8,1x,i8,1x,"(",1x,',&
                   nmodes,'(1x, i3),1x,")",1x,"<- ",1x,"(",1x,',nmodes,'(1x, i3),1x,")",1x,',nformat,'(1x, es16.8))'
    !
    !
    ! Prepare the table header
    !
    select case (trim(intensity%action))
      !
      case('ABSORPTION','EMISSION')
       !
       if (.not.job%exomol_format) then
         !
         !write(my_fmt_tm,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
         !                   '(i4,1x,a3,3x,"<-",i4,1x,a3,3x,a1,2x,f13.6,1x,"<-",1x, f13.6,1x,f12.6,2x,"(",a3,";",i3,")",1x,"(",',&
         !                   nclasses,'a3,";",',nmodes,'(1x,i3),")",2x,"<- ","(",a3,";",i3,")",1x,"(",',&
         !                   nclasses,'a3,";",',nmodes,'(1x,i3),")",2(1x,es15.8),i8,2x,"(",',&
         !                   nmodes,'(1x, i3),")",2x,"<- ",1x,"(",',nmodes,'(1x, i3),")",3(',nformat,'(1x,f16.8,1x,3i1)))' 
         !
         write(my_fmt1,'(a,i0,a,i0,a,i0,a,i0,a)') "(/t4,a1,t6,a8,t17,a1,t19,a5,t25,a3,t35,a2,t42,a2,t50,a2,t62,a5,t85,",nclasses,&
                                   "(4x),1x,",nmodes,"(4x),3x,a2,14x,",nclasses,"(4x),1x,",nmodes,&
                                   "(4x),8x,a7,10x,a5,12x,a7,12x,a2,8x,a2,8x,a1)"

         !
         !write(my_fmt1,'(a)') "(/t4,a1,t6,a8,t17,a1,t19,a5)"
         !write(out,"(/t4a1,t6a8,t17a1,t19a5,t25a3,t35a1,t42a2,t50a1,t62a5,t85,<nclasses>(4x),1x,<nmodes>(4x),3x,a2,14x,<nclasses>(4x),1x,<nmodes>(4x),8x,a7,10x,a5,12x,a7,12x,a1,8x,a1,8x,a1)") 'J','Gamma <-','J','Gamma','Typ','Ef','<-','Ei','nu_if','<-','S(f<-i)','A(if)','I(f<-i)','Ni','Nf','N'
         !
         write(out,my_fmt1) 'J','Gamma <-','J','Gamma','Typ','Ef','<-','Ei' ,'nu_if','<-','S(f<-i)','A(if)','I(f<-i)','Ni','Nf','N'
         !
      else
         !
         write(out,"(t11,a4,t18,a2,t25,a4,t34,a5,t55,a5)") 'ID_f','<-','ID_i','A(if)','nu_if'
         !
      endif
      !
      case('TM')
       !
       write(my_fmt_tm,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
                               '(i4,1x,a3,3x,"<-",i4,1x,a3,3x,a1,2x,f13.6,1x,"<-",1x, f13.6,1x,f12.6,2x,"(",a3,";",i3,")",1x,"(",',&
                               nclasses,'a3,";",',nmodes,'(1x,i3),")",2x,"<- ","(",a3,";",i3,")",1x,"(",',&
                               nclasses,'a3,";",',nmodes,'(1x,i3),")",2(1x,es15.8),i8,2x,"(",',&
                               nmodes,'(1x, i3),")",2x,"<- ",1x,"(",',nmodes,'(1x, i3),")",3(',nformat,'(1x,f16.8)))' 
                               !'(1x,f16.8,1x,3i1)))' 
                               !
       write(my_fmt1,'(a,i0,a,i0,a,i0,a,i0,a)') &
                                  "(/t4,a1,t6,a8,t17,a1,t19,a5,t2,5a3,t35,a2,t42,a2,t52,a2,t65,a5,t84,",nclasses,"(3x),1x",nmodes,&
                                  &"(3x),9x,a2,8x,",nclasses,"(3x),1x,",nmodes,"(3x),17x,a8,8x,a1,12x,a1,18x,a1,18x,a1)"
       !
       write(out,my_fmt1) 'J','Gamma <-','J','Gamma','Typ','Ef','<-','Ei','nu_if','<-','TM(f<-i)','N','x','y','z'
       !
    end select 
    !
    ! ---------------------------------
    ! the actual intensity calculations
    ! --------------------------------- 
    !
    itransit = 0
    !
    call TimerStart('Intensity loop')
    !
    ! loop over initial states
    !
    Ilevels_loop: do ilevelI = 1,nlevels
      !
      ! start measuring time per line
      !
      indI = eigen(ilevelI)%jind
      !
      !dimension of the bases for the initial states
      !
      dimenI = bset_contr(indI)%Maxcontracts
      !
      !energy, quanta, and gedeneracy order of the initial state
      !
      jI = eigen(ilevelI)%jval
      energyI = eigen(ilevelI)%energy
      igammaI  = eigen(ilevelI)%igamma
      quantaI(0:nmodes) = eigen(ilevelI)%quanta(0:nmodes)
      normalI(0:nmodes) = eigen(ilevelI)%normal(0:nmodes)
      ndegI   = eigen(ilevelI)%ndeg
      nsizeI = bset_contr(indI)%nsize(igammaI)
      !
      call energy_filter_lower(jI,energyI,quantaI,normalI(0),passed)
      !
      if (.not.passed) cycle
      !
      ! where the eigenvector is stored 
      !
      unitI = Jeigenvec_unit(indI,igammaI) 
      !
      if (unitI==-1) stop 'This file is not supposed to be accessed'
      !
      !read eigenvector of initial state
      !
      call TimerStart('Reading eigenfuncs-I')
      !
      if (istored(ilevelI)==0) then
          !
          irec = isaved(ilevelI)
          !
          if (.not.vector_diagonal) then
             !
             read(unitI,rec=irec) vecI(1:nsizeI)
             !
          else
             !
             vecI(:) = 0 
             vecI(eigen(ilevelI)%icoeff) = 1.0_rk
             !
          endif
          !
      else
          !
          iram_ = istored(ilevelI)
          !
          vecI(1:nsizeI) = vec_ram(indI,igammaI)%mat(1:nsizeI,iram_)
          !
      endif
      !
      call TimerStop('Reading eigenfuncs-I')
      !
      ! Compute the half-linestrength
      !
      half_linestr = 0
      !
      call TimerStart('Half_linestrength')
      !
      do indF = 1, nJ
        !
        jF = Jval(indF)
        !
        dimenF = bset_contr(indF)%Maxcontracts
        !
        ! Check if it is really necessary to start the calculations for the given levelI -> jF, 
        ! i.e. one can skip the rest if no transitions will start from the given ilevelI and 
        ! finish anywehere at J= jF. 
        !
        call TimerStart('Intens_Filter-2')
        !
        passed = .false.
        !
        !omp parallel do private(ilevelF,energyF,igammaF,quantaF,passed_) schedule(guided) reduction(+:passed)
        do ilevelF = 1, nlevels
          !
          if (eigen(ilevelF)%jval/=jF) cycle 
          !
          energyF = eigen(ilevelF)%energy
          igammaF = eigen(ilevelF)%igamma        
          quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
          !
          call intens_filter(jI,jF,energyI,energyF,igammaI,igammaF,quantaI,quantaF,igamma_pair,passed)
          !
          if (passed) exit
          !
          !passed = passed_
          !
        enddo
        !omp end parallel do 
        !
        call TimerStop('Intens_Filter-2')
        !
        if (.not.passed) cycle
        !
        select case (trim(intensity%action))
          !
        case('ABSORPTION','EMISSION')
          !
          igammaF = igamma_pair(igammaI)
          !
          ! The degeneracy of the basis vector have to be the same for the intial and final states:
          !
          ndegF  = ndegI
          idegF = 1
          !
          if ((jI/=intensity%J(1).or.jF/=intensity%J(1)).and. abs(jI-jF)<=1.and.jI+jF>=1) then 
            !
            do idegI = 1, ndegI
              !
              call degeneracy_filter(igammaI,igammaF,idegI,idegF,passed)
              if (.not.passed) cycle
              !
              call TimerStart('Unfolding-to-primitives')
              !
              vec = 0
              !
              !$omp parallel do private(irootI,irow,ib,iterm,dtemp0,nelem,ielem,isrootI) shared(vec) schedule(static)
              do irootI = 1,dimenI
                 !
                 irow = bset_contr(indI)%icontr2icase(irootI,1)
                 ib   = bset_contr(indI)%icontr2icase(irootI,2)
                 !
                 iterm = ijterm(indI)%kmat(irow,igammaI)
                 !
                 nelem = bset_contr(indI)%irr(igammaI)%N(irow)
                 !
                 if (.not.vector_diagonal.or..false.) then
                   !
                   dtemp0 = 0
                   !
                   do ielem = 1,nelem
                      !
                      isrootI = iterm+ielem 
                      !
                      dtemp0 = dtemp0 + vecI(isrootI)*bset_contr(indI)%irr(igammaI)%repres(isrootI,idegI,ib)
                      !
                   enddo
                   !
                 else
                    !
                    isrootI = eigen(ilevelI)%icoeff
                    !
                    if (isrootI-iterm<1.or.isrootI-iterm>nelem) cycle
                    !
                    dtemp0 = vecI(isrootI)*bset_contr(indI)%irr(igammaI)%repres(isrootI,idegI,ib)
                    !
                 endif 
                 !
                 vec(irootI) = dtemp0
                 !
              enddo
              !$omp end parallel do
              !
              call TimerStop('Unfolding-to-primitives')
              !
              call TimerStart('Pre-screening')
              !
              !vecPack = pack(vec,abs(vec)>intensity%threshold%coeff)
              !cdimenI = count(abs(vec)>intensity%threshold%coeff)
              !
              cdimenI    = 0 
              icoeffI(:) = 0
              cdimen_kblock  = 0
              icoeff_kblock  = 0  
              itau_kblock    = 0
              vecPack_kblock = 0 
              do idimen = 1, dimenI
                 if (abs(vec(idimen)) > intensity%threshold%coeff) then
                    cdimenI = cdimenI + 1
                    icoeffI(cdimenI) = idimen
                    vecPack(cdimenI) = vec(idimen)
                    !
                    kI = bset_contr(indI)%k(idimen)
                    cdimen_kblock(kI) = cdimen_kblock(kI) + 1
                    icoeff_kblock(cdimen_kblock(kI),kI) = idimen
                    vecPack_kblock(cdimen_kblock(kI),kI) = vec(idimen)
                    !
                    ktauI = bset_contr(indI)%ktau(idimen)
                    itau_kblock(cdimen_kblock(kI),kI) = mod(ktauI,2_ik)
                    !
                 end if
              end do
              !
              call TimerStop('Pre-screening')
              !
              if (job%rotsym_do) then 
                !
                call do_1st_half_linestrength_rotsym_symmvec(jI,jF,indI,indF,cdimenI,icoeffI,vecPack,&
                                              half_linestr(:,indF,idegI,1))
                                              !
              else
                !
                call do_1st_half_linestrength_III_symmvec(jI,jF,indI,indF,jmax,dimenmax,&
                cdimen_kblock,icoeff_kblock,vecPack_kblock,itau_kblock,threej,half_linestr(:,indF,idegI,1))
                !call do_1st_half_linestrength_II_symmvec(jI,jF,indI,indF,cdimenI,icoeffI,vecPack,&
                !                              jmax,threej,half_linestr(:,indF,idegI,1))
              endif 
              !
            enddo
            !
          endif 
          !
        case('TM')
            !
            do idegI = 1, ndegI
              !
              !$omp parallel do private(irootI,irow,ib,iterm,dtemp0,nelem,ielem,isrootI) shared(vec) schedule(static)
              do irootI = 1,dimenI
                 !
                 irow = bset_contr(indI)%icontr2icase(irootI,1)
                 ib   = bset_contr(indI)%icontr2icase(irootI,2)
                 !
                 iterm = ijterm(indI)%kmat(irow,igammaI)
                 !
                 dtemp0 = 0
                 !
                 nelem = bset_contr(indI)%irr(igammaI)%N(irow)
                 !
                 do ielem = 1,nelem
                    !
                    isrootI = iterm+ielem 
                    !
                    dtemp0 = dtemp0 + vecI(isrootI)*bset_contr(indI)%irr(igammaI)%repres(isrootI,idegI,ib)
                    !
                 enddo
                 !
                 vec(irootI) = dtemp0
                 !
              enddo
              !$omp end parallel do
              !
              call do_1st_half_tm_symmvec(indI,indF,vec,&
                                     half_linestr(:,indF,idegI,:))
            enddo
            !
        end select
        !
      enddo
      !
      call TimerStop('Half_linestrength')
      !
      !loop over final states
      !
      !$omp parallel private(vecF,vec_,alloc_p) reduction(max:max_intens_state)
      allocate(vecF(dimenmax_),vec_(dimenmax),stat = alloc_p)
      if (alloc_p/=0) then
          write (out,"(' dipole: ',i9,' trying to allocate array -vecF')") alloc_p
          stop 'dipole-vecF - out of memory'
      end if
      !
      !$omp do private(ilevelF,indF,dimenF,jF,energyF,igammaF,quantaF,normalF,ndegF,nsizeF,unitF,passed,branch,&
      !$omp& nu_if,irec,iram_,linestr_deg,idegF,irootF,irow,ib,iterm,nelem,dtemp0,ielem,isrootF,idegI,linestr,A_einst,boltz_fc,&
      !$omp& absorption_int,tm_deg,dmu,icontrF) reduction(+:itransit) schedule(static) 
      Flevels_loop: do ilevelF = 1,nlevels
         !
         indF = eigen(ilevelF)%jind
         !
         !dimension of the bases for the final state
         !
         dimenF = bset_contr(indF)%Maxcontracts
         !
         !energy, quanta, and gedeneracy order of the final state
         !
         jF = eigen(ilevelF)%jval
         energyF = eigen(ilevelF)%energy
         igammaF  = eigen(ilevelF)%igamma
         quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes)
         normalF(0:nmodes) = eigen(ilevelF)%normal(0:nmodes)
         ndegF   = eigen(ilevelF)%ndeg
         nsizeF = bset_contr(indF)%nsize(igammaF)
         !
         ! where the eigenvector is stored 
         !
         unitF = Jeigenvec_unit(indF,igammaF)
         !
         if (unitF==-1) stop 'This file is not supposed to be accessed'
         !
         call energy_filter_upper(jF,energyF,quantaF,normalF(0),passed)
         !
         if (.not.passed) cycle Flevels_loop
         !
         call intens_filter(jI,jF,energyI,energyF,igammaI,igammaF,quantaI,quantaF,igamma_pair,passed)
         !
         if (.not.passed) cycle Flevels_loop
         !
         ! Which PQR branch this transition belong to ->
         ! 
         branch = PQR_branch(jI,jF)
         !
         nu_if = energyF - energyI 
         if (trim(intensity%action)=='EMISSION') nu_if = -nu_if 
         !
         ! Count the processed transitions 
         !
         itransit = itransit + 1
         !
         !read eigenvector of the final state
         !
         if (istored(ilevelF)==0) then
           !
           if (.not.vector_diagonal) then
              !
              !call TimerStart('Reading eigenfuncs')
              !
              irec = isaved(ilevelF)
              !
              !read(unitF,rec=irec) vec_ram(indF,igammaF)%mat(1:nsizeF,iram(indF,igammaF))
              !
              !$omp critical
              read(unitF, rec = irec) vecF(1:nsizeF)
              !$omp end critical 
              !
              !call TimerStop('Reading eigenfuncs')
              !
           else
              !
              vecF = 0 
              vecF(eigen(ilevelF)%icoeff) = 1.0_rk
              !
           endif
           !
         else
           !
           iram_ = istored(ilevelF)
           !
           vecF(1:nsizeF) = vec_ram(indF,igammaF)%mat(1:nsizeF,iram_)
           !
         endif
         !
         select case (trim(intensity%action))
           !
           case default
             !
             stop 'only ABSORPTION and TM are properly coded so far'
             !
           case('ABSORPTION')
             !
             linestr_deg = 0 
             !
             !loops over degenerate states
             !
             !call TimerStart('Linestrength')
             !
             linestr_deg = 0
             !
             Idegs_loop: do idegF = 1, ndegF
                !
                if (intensity%reduced.and.idegF/=1) cycle
                !
                if (.not.vector_diagonal) then
                   !
                   vec_ = 0
                   !
                   !omp parallel do private(irootF,irow,ib,iterm,dtemp0,nelem,ielem,isrootF) shared(vec) schedule(static)
                   do irootF = 1,dimenF
                      !
                      irow = bset_contr(indF)%icontr2icase(irootF,1)
                      ib   = bset_contr(indF)%icontr2icase(irootF,2)
                      !
                      !icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
                      !
                      iterm = ijterm(indF)%kmat(irow,igammaF)
                      !
                      nelem = bset_contr(indF)%irr(igammaF)%N(irow)
                      !
                      dtemp0 = 0
                      !
                      do ielem = 1,nelem
                         !
                         isrootF = iterm+ielem 
                         !
                         dtemp0 = dtemp0 + vecF(isrootF)*bset_contr(indF)%irr(igammaF)%repres(isrootF,idegF,ib)
                         !
                      enddo
                      !
                      !else
                      !   !
                      !   isrootF = eigen(ilevelF)%icoeff
                      !   !
                      !   if (isrootF-iterm<1.or.isrootF-iterm>nelem) cycle
                      !   !
                      !   dtemp0 = vecF(isrootF)*bset_contr(indF)%irr(igammaF)%repres(isrootF,idegF,ib)
                      !   !
                      !endif 
                      !
                      vec_(irootF) = dtemp0
                      !
                   enddo
                   !omp end parallel do
                   !
                   Fdegs_loop: do idegI = 1,ndegI
                     !
                     !if (intensity%reduced.and.ndegF/=1.and.ndegI/=1.and.idegI/=1) cycle
                     !
                     call degeneracy_filter(igammaI,igammaF,idegI,idegF,passed)
                     !
                     if (.not.passed) cycle
                     !
                     ! complete the line strength 
                     !
                     linestr_deg(idegI,idegF) = ddot(dimenF,half_linestr(1:dimenF,indF,idegI,1),1,vec_(1:dimenF),1)
                     !
                     !linestr_deg(idegI,idegF) = dtemp
                     !
                   end do Fdegs_loop
                   !
                else 
                   !
                   !omp parallel do private(irootF,idegi,hf_ls,irow,ib,iterm,dtemp0,nelem,ielem,isrootF) shared(vec) schedule(static)
                   do irootF = 1,dimenF
                      !
                      do idegI = 1,ndegI
                        !
                        call degeneracy_filter(igammaI,igammaF,idegI,idegF,passed)
                        !
                        if (.not.passed) cycle
                        !
                        linestr = half_linestr(irootF,indF,idegI,1)
                        !
                        if (abs(linestr)<intensity%threshold%linestrength) cycle
                        !
                        irow = bset_contr(indF)%icontr2icase(irootF,1)
                        ib   = bset_contr(indF)%icontr2icase(irootF,2)
                        !
                        !icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
                        !
                        iterm = ijterm(indF)%kmat(irow,igammaF)
                        !
                        nelem = bset_contr(indF)%irr(igammaF)%N(irow)
                        !
                        dtemp0 = 0
                        !
                        isrootF = eigen(ilevelF)%icoeff
                        !
                        if (isrootF-iterm<1.or.isrootF-iterm>nelem) cycle
                        !
                        dtemp0 = vecF(isrootF)*bset_contr(indF)%irr(igammaF)%repres(isrootF,idegF,ib)
                        !                        
                        linestr_deg(idegI,idegF) = linestr_deg(idegI,idegF) + dtemp0*linestr
                        !
                     enddo
                     !
                   enddo
                   !omp end parallel do  
                   !
                endif
                !
             end do Idegs_loop
             !
             !call TimerStop('Linestrength')
             !
             !sum up all degenerate components
             !
             linestr = sum( linestr_deg(1:ndegI,1:ndegF)**2 )/real(ndegI,rk)
             !
             ! swtich to a reduced C3v/D3h/Td case by commenting this line and undcommenting the next one
             !
             if (intensity%reduced.and.ndegF/=1.and.ndegI/=1) linestr  = linestr*real(ndegI,rk)
             !
             !linestr = sum( linestr_deg(1:ndegI,1:ndegF)**2 )
             !
             ! calculate the intensity 
             !
             A_einst = A_coef_s_1*real(2*jI+1,rk)*linestr*abs(nu_if)**3
             !
             linestr = linestr * intensity%gns(igammaI) * real( (2*jI + 1)*(2 * jF + 1),rk )
             !
             boltz_fc = abs(nu_if) * exp(-(energyI-intensity%ZPE) * beta) * (1.0_rk - exp(-abs(nu_if) * beta))&
                        / intensity%part_func
             !
             ! intensity in cm/mol
             !
             absorption_int = linestr * intens_cm_molecule * boltz_fc
             !
             !
             if (absorption_int>=intensity%threshold%intensity.and.linestr>=intensity%threshold%linestrength) then 
               !
               iram_ = istored(ilevelF)
               !
               !
               !$omp critical
               if (job%exomol_format) then
                 !
                 write(out, "( i12,1x,i12,1x,1x,es16.8,1x,f16.6,' ||')")&
                              iID(ilevelF),iID(ilevelI),A_einst,nu_if     
                              !
               elseif (intensity%output_short) then 
                 !
                 write(out, "( i4,1x,i2,9x,i4,1x,i2,3x,&
                              &2x, f11.4,3x,f11.4,1x,f11.4,2x,&
                              &2(1x,es16.8),3x,i6,1x,2x,i6,1x,'||')")&
                              !
                              jF,igammaF,jI,igammaI, & 
                              energyF-intensity%ZPE,energyI-intensity%ZPE,nu_if,                 &
                              A_einst,absorption_int,&
                              eigen(ilevelF)%ilevel,eigen(ilevelI)%ilevel
               else

                 !
                !write(out, "( (i4, 1x, a4, 3x),'<-', (i4, 1x, a4, 3x),a1,&
                !             &(2x, f11.4,1x),'<-',(1x, f11.4,1x),f11.4,2x,&
                !             &'(',1x,a3,x,i3,1x,')',1x,'(',1x,<nclasses>(x,a3),1x,<nmodes>(1x, i3),1x,')',1x,'<- ',   &
                !             &'(',1x,a3,x,i3,1x,')',1x,'(',1x,<nclasses>(x,a3),1x,<nmodes>(1x, i3),1x,')',1x,   &
                !             & 3(1x, es16.8),2x,(1x,i6,1x),'<-',(1x,i6,1x),i8,1x,i8,&
                !             1x,'(',1x,<nmodes>(1x, i3),1x,')',1x,'<- ',1x,'(',1x,<nmodes>(1x, i3),1x,')',1x,& 
                !             <nformat>(1x, es16.8))")  &
                              !
                 write(out,my_fmt) &
                              jF,sym%label(igammaF),jI,sym%label(igammaI),branch, &
                              energyF-intensity%ZPE,energyI-intensity%ZPE,nu_if,                 &
                              eigen(ilevelF)%cgamma(0),eigen(ilevelF)%krot,&
                              eigen(ilevelF)%cgamma(1:nclasses),eigen(ilevelF)%quanta(1:nmodes), &
                              eigen(ilevelI)%cgamma(0),eigen(ilevelI)%krot,&
                              eigen(ilevelI)%cgamma(1:nclasses),eigen(ilevelI)%quanta(1:nmodes), &
                              linestr,A_einst,absorption_int,&
                              eigen(ilevelF)%ilevel,eigen(ilevelI)%ilevel,&
                              itransit,istored(ilevelF),normalF(1:nmodes),normalI(1:nmodes),&
                              linestr_deg(1:ndegI,1:ndegF)
                              !            


               endif
               !$omp end critical
               !             
             endif
             !
           case('TM')
             !
             if (jI/=0.or.jF/=0) then
               write(out,"(' TM: can be used only for J=0, not for ',2i4)") jI,jF
               stop 'TM: can be used only for J=0' 
             endif
             !
             tm_deg = 0
             !
             !loops over degenerate states
             !
             do idegI = 1, ndegI
               do idegF = 1, ndegF
                !
                ! complete the tm
                !
                dmu = 0
                !
                !omp parallel do private(icontrF,irow,ib,iterm,dtemp0,nelem,ielem,isrootF) reduction(+:dmu) schedule(guided)
                do icontrF = 1,bset_contr(1)%Maxcontracts
                   !
                   irow = bset_contr(1)%icontr2icase(icontrF,1)
                   ib   = bset_contr(1)%icontr2icase(icontrF,2)
                   !
                   !icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
                   !
                   iterm = ijterm(1)%kmat(irow,igammaF)
                   !
                   dtemp0 = 0
                   !
                   nelem = bset_contr(1)%irr(igammaF)%N(irow)
                   !
                   do ielem = 1,nelem
                      !
                      isrootF = iterm+ielem 
                      !
                      dtemp0 = dtemp0 + vecF(isrootF)*bset_contr(indF)%irr(igammaF)%repres(isrootF,idegF,ib)
                      !
                   enddo
                   !
                   !dtemp = dtemp + half_linestr(idegI,1,icontrF)*dtemp0
                   !
                   dmu(:) = dmu(:) + half_linestr(icontrF,indF,idegI,:)*dtemp0
                   !
                enddo
                !omp end parallel do
                !
                tm_deg(:,idegI,idegF) = dmu(:) 
                !
               end do
             end do
             !
             !sum up all degenerate components
             !
             linestr = sqrt( sum( tm_deg(1:3,1:ndegI,1:ndegF)**2 ) )
             !
             ! Einstein A (vib) coeff
             !
             A_einst = A_coef_s_1*linestr**2*abs(nu_if)**3
             !
             boltz_fc = abs(nu_if) * exp(-(energyI-intensity%ZPE) * beta) * (1.0_rk - exp(-abs(nu_if) * beta))&
                        / intensity%part_func
             !
             ! intensity in cm/mol
             !
             absorption_int = linestr**2 * intens_cm_molecule * boltz_fc
             !
             if (linestr>=intensity%threshold%intensity) then 
               !
               if (job%exomol_format) then
                 !
                 write(out, "( i12,1x,i12,1x,1x,es16.8,1x,f16.6,1x,es16.8,' ||')")&
                              iID(ilevelF),iID(ilevelI),A_einst,nu_if,absorption_int
                              !
               else
                 !
                 !nformat = ndegI*ndegF
                 !
                 !$omp critical
                 !write(out, "( (i4, 1x, a3, 3x),'<-', (i4, 1x, a3, 3x),a1,&
                 !             &(2x, f13.6,1x),'<-',(1x, f13.6,1x),f12.6,2x,&
                 !             &'(',a3,';',i3,')',1x,'(',<nclasses>a3,';',<nmodes>(1x, i3),')',2x,'<- ',   &
                 !             &'(',a3,';',i3,')',1x,'(',<nclasses>a3,';',<nmodes>(1x, i3),')',   &
                 !             &2(1x,es15.8),i8,&
                 !             &2x,'(',<nmodes>(1x, i3),')',2x,'<- ',   &
                 !             &1x,'(',<nmodes>(1x, i3),')',   &
                 !             &3(<nformat>( 1x,f16.8,1x,3i1) ) )") &
                 write(out,my_fmt_tm) &
                              jF,sym%label(igammaF),jI,sym%label(igammaI),branch, &
                              energyF-intensity%ZPE,energyI-intensity%ZPE,nu_if ,&
                              eigen(ilevelF)%cgamma(0),eigen(ilevelF)%krot,&
                              eigen(ilevelF)%cgamma(1:nclasses),eigen(ilevelF)%quanta(1:nmodes), &
                              eigen(ilevelI)%cgamma(0),eigen(ilevelI)%krot,&
                              eigen(ilevelI)%cgamma(1:nclasses),eigen(ilevelI)%quanta(1:nmodes), &
                              linestr,absorption_int,itransit,&
                              normalF(1:nmodes),normalI(1:nmodes),&
                              !( ( ( ( tm_deg(imu,idegI,idegF),idegF,idegI,imu ), idegF=1,ndegF ),idegI=1,ndegI),imu=1,3 )
                              ( ( ( ( tm_deg(imu,idegI,idegF) ), idegF=1,ndegF ),idegI=1,ndegI),imu=1,3 )
                 !$omp end critical
                 !
               endif 
               !'i4,1x,a3,3x,"<-",i4,1x,a3,3x,a1,2x,f13.6,1x,"<-",1x, f13.6,1x,f12.6,2x,"(",a3,";",i3,")",1x,"(",2a3,";",3(1x,i3),")",2x,"<- ","(",a3,";",i3,")",1x,"(",2a3,";",3(1x,i3),")",2(1x,es15.8),i8,2x,"(",3(1x, i3),")",2x,"<- ",1x,"(",3(1x, i3),")",3(1(1x,f16.8,1x,3i1))'
               !
            endif 
            !
            ! estiimate the maximal intensity for each state needed for the TM-pruning of the basis set
            !
            if (intensity%pruning) then 
              !
              max_intens_state(ilevelF) = max(max_intens_state(ilevelF),absorption_int)
              max_intens_state(ilevelI) = max(max_intens_state(ilevelI),absorption_int)
              !
            endif
            !
         end select
         !
      end do Flevels_loop
      !$omp enddo
      !
      deallocate(vecF,vec_)    
      !$omp end parallel
      !
      ! stop counting time per block and estimate seconds per line 
      !
      call TimerProbe('Intensity loop',real_time,cpu_time)
      !
      time_per_line = real_time/real(min(itransit,1),rk)
      total_time_predict = time_per_line*real(Ntransit,rk)/3600.0
      !
      !if (job%verbose>=5) then
      !   write(out,"(/t2,i,' lines ',f12.2,' s,  ',g12.2,'s per line; total estimate for ',i,' lines:',f12.2,'h')") itransit,real_time,time_per_line,Ntransit,total_time_predict
      !endif
      !
      if (job%verbose>=6) then
          write(out,"('--- ',t4,f18.6,2x,i0,' l ',f12.2,' s, (',g12.2,' l/s ); Ttot= ',f12.2,'hrs.')") &
                    energyI-intensity%ZPE,itransit,real_time,1.0_rk/time_per_line,total_time_predict
      endif
      !
      if (mod(ilevelI,min(100000,nlevelI))==0.and.(int(total_time_predict/intensity%wallclock)/=0).and.job%verbose>=5) then
         !
         write(out,"(/'Recommended energy distribution for ',f12.2,' h limit:')") intensity%wallclock
         !
         int_increm = max(Ntransit/int(total_time_predict/intensity%wallclock),1)
         !
         ! print energy ranges
         !
         ener_ = intensity%ZPE
         Ntransit_ = 0
         !
         do ilevelI_ = 1, nlevels
           !
           jI = eigen(ilevelI_)%jval
           energyI = eigen(ilevelI_)%energy
           igammaI  = eigen(ilevelI_)%igamma
           quantaI(0:nmodes) = eigen(ilevelI_)%quanta(0:nmodes)
           normalI(0:nmodes) = eigen(ilevelI_)%normal(0:nmodes)
           !
           call energy_filter_lower(jI,energyI,quantaI,normalI(0),passed)
           !
           if (.not.passed) cycle
           !
           do ilevelF = 1, nlevels
              !
              jF = eigen(ilevelF)%jval 
              energyF = eigen(ilevelF)%energy
              igammaF = eigen(ilevelF)%igamma        
              quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
              !
              call intens_filter(jI,jF,energyI,energyF,igammaI,igammaF,quantaI,quantaF,igamma_pair,passed)
              !
              if (passed) then 
                !
                Ntransit_ = Ntransit_ + 1
                !
                if (mod(Ntransit_,int_increm)==0.and.energyI>ener_) then
                   !
                   write(out,"('|  ',2f16.0,i0)") ener_-intensity%ZPE,energyI-intensity%ZPE,Ntransit_
                   ener_ = energyI
                endif
                !
              endif 
              !
           enddo
           !
         enddo
         !
         write(out,"('|  ',2f16.0,i0)") ener_-intensity%ZPE,intensity%erange_low(2),Ntransit_
         !
      endif
      !
      if (intensity%pruning) then 
        !
        do indF = 1, nJ
           !
           do igammaF = 1,sym%Nrepresen
             !
             write(jchar, '(i4)') jval(indF)
             write(gchar, '(i3)') igammaF
             !
             filename = trim(job%eigenfile%filebase)//'_intens'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
             !
             ioname = 'max state intensities'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))
             !
             call IOstart(trim(ioname),iounit)
             !
             rewind(iounit)
             !
             do irootF = 1, bset_contr(indF)%nsize(igammaF)
               !
               ilevelF = istate2ilevel(indF,igammaF,irootF)
               !
               energyF = 0 ; absorption_int = 0
               if (ilevelF/=0) then 
                 energyF =  eigen(ilevelF)%energy
                 absorption_int = max_intens_state(ilevelF)
               endif
               !
               write(iounit,"(i5,i3,1x,i9,2x,f20.12,2x,e15.8)") jval(indF),igammaF,irootF,energyF,absorption_int
               !
             enddo
             !
           enddo
           !
        enddo
        !
      endif 
      !
      if (job%verbose>=5) call TimerReport
      !
      call flush(out) 
      !
    end do Ilevels_loop
    !
    call TimerStop('Intensity loop')
    !
    ! write out the list of states with maximal intensities
    !
    if (intensity%pruning) then 
      !
      ioname = 'max state intensities'
      !
      call IOstart(trim(ioname),iounit)
      !
      do indI = 1, nJ
         !
         do igammaI = 1,sym%Nrepresen
           !
           write(jchar, '(i4)') jval(indI)
           write(gchar, '(i3)') igammaI
           !
           filename = trim(job%eigenfile%filebase)//'_intens'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
           !
           ioname = 'max state intensities'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))
           !
           call IOstart(trim(ioname),iounit)
           !
           close(unit = iounit,status='keep')
           !
           call IOstop(trim(ioname))
           !
         enddo
         !
      enddo
      !
    endif 
    !
    deallocate(vecI, vecPack, vec, icoeffI)
    call ArrayStop('intensity-vectors')

    if (.not.job%rotsym_do) then
      deallocate(vecPack_kblock,icoeff_kblock,cdimen_kblock,itau_kblock)
      call ArrayStop('intensity-vectors-kblock')
    endif
    !
    deallocate(half_linestr)
    call ArrayStop('half_linestr')
    !
    do irep = 1,Nrepresen
     do jind = 1,nJ
      !
      nullify(vec_ram(jind,irep)%mat)
      !
      if (associated(vec_ram(jind,irep)%mat)) deallocate(vec_ram(jind,irep)%mat)
      !
     enddo
    enddo
    !
    call ArrayStop('vec_ram')
    !
    do jind = 1,nJ
      !
      if (associated(ijterm(jind)%kmat)) deallocate(ijterm(jind)%kmat)
      !
    enddo
    call ArrayStop('ijterm')
    !
    deallocate(nlevelsG,ilevelsG,ram_size,iram)
    call ArrayStop('nlevelsG')
    call ArrayStop('ilevelsG')
    call ArrayStop('ram_size')
    call ArrayStop('iram')
    !
    deallocate(istored)
    call ArrayStop('istored')
    !
    if (allocated(iID) ) then 
      deallocate(iID)
      call ArrayStop('iID')
    endif
    !
    deallocate(isaved)
    call ArrayStop('isaved')
    !
    if (intensity%pruning) then
      !
      deallocate(max_intens_state)
      call ArrayStop('max_intens_state')
      !
    endif
    !
    call TimerStop('Intensity calculations')
    !
  end subroutine dm_intensity_symmvec




 !
 !
 ! Here we restore the vibrational (J=0) contracted matrix elements of the dipole moment
 !
 subroutine restore_vib_matrix_elements 
   !
   implicit none 
   !
   integer(ik)        :: chkptIO,alloc
   character(len=cl)  :: job_is

   character(len=20)  :: buf20
   character(len=4)   :: buf4
   character(len=4)   :: jchar
   
   integer(ik)        :: ncontr_t,imu,imu_t,i,j
   integer(hik)       :: matsize,rootsize,rootsize2
   !
   real(rk),allocatable  :: dipole_(:,:)
   !
   character(len=cl) :: filename
   !
   if (job%verbose>=4) then 
     if (.not.job%IOextF_divide) then
       !
       write(out, '(/a, 1x, a)') 'read vibrational contracted matrix elements from chk', trim(job%extFmat_file)
       !
     else
       !
       write(out, '(/a, 1x, a)') 'read vibrational contracted matrix elements from files', trim(job%extmat_suffix)
       !
     endif
     !
   endif 
   !
   job_is ='external field contracted matrix elements for J=0'
   call IOStart(trim(job_is),chkptIO)
   !
   if (.not.job%IOextF_divide) then
      !
      if (job%verbose>=4) write(out, '(a, 1x, a)') 'chk = file', trim(job%extFmat_file)
      !
      open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=job%extFmat_file)
      !
      read(chkptIO) buf20
      if (buf20/='Start external field') then
        write (out,"(' restore_vib_matrix_elements ',a,' has bogus header: ',a)") job%extFmat_file,buf20
        stop 'restore_vib_matrix_elements - bogus file format'
      end if
      !
      read(chkptIO) ncontr_t
      !
      if (bset_contr(1)%Maxcontracts/=ncontr_t) then
        write (out,"(' Dipole moment checkpoint file ',a)") job%extFmat_file
        write (out,"(' Actual and stored basis sizes at J=0 do not agree  ',2i8)") bset_contr(1)%Maxcontracts,ncontr_t
        stop 'restore_vib_matrix_elements - in file - illegal ncontracts '
      end if
      !
   endif
   !
   ncontr_t = bset_contr(1)%Maxcontracts
   !
   if (job%verbose>=5) write(out,"(/'restore_vib_matrix_elements...: Number of elements: ',i8)") ncontr_t
   !
   rootsize  = int(ncontr_t*(ncontr_t+1)/2,hik)
   rootsize2 = int(ncontr_t*ncontr_t,hik)
   !
   matsize = rootsize2*3
   !
   if (job%verbose>=5) write(out,"(/'allocate 4 dipole_me matrices with ',i22,' elements each...')") rootsize2
   allocate(dipole_me(ncontr_t,ncontr_t,3),stat=alloc)
   call ArrayStart('dipole_me',alloc,1,kind(dipole_me),matsize)
   !
   allocate(dipole_(ncontr_t,ncontr_t),stat=alloc)
   call ArrayStart('dipole_',alloc,1,kind(dipole_),rootsize2)
   !
   do imu = 1,3
     !
     if (job%verbose>=4) write(out,'("imu = ",i3)') imu 
     !
     if (job%IOextF_divide) then
       !
       write(jchar, '(i4)') imu
       !
       filename = trim(job%extmat_suffix)//trim(adjustl(jchar))//'.chk'
       !
       open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=filename)
       !
       read(chkptIO) buf4
       if (buf4/='extF') then
         write (out,"(' restore_vib_matrix_elements ',a,' has bogus header: ',a)") filename,buf4
         stop 'restore_vib_matrix_elements - bogus file format'
       end if
       !
     else
       !
       read(chkptIO) imu_t
       if (imu_t/=imu) then
         write (out,"(' restore_vib_matrix_elements ',a,' has bogus imu - restore_vib_matrix_elements: ',i8,'/=',i8)") imu_t,imu
         stop 'restore_vib_matrix_elements - bogus imu restore_vib_matrix_elements'
       end if
       !
     endif
     !
     if (job%verbose>=5) write(out,"('read dipole_ ...')") !,advance='NO'
     !
     read(chkptIO) dipole_
     !
     if (job%verbose>=5) write(out,"('copy to dipole_me ...')")
     !
     dipole_me(1:ncontr_t,1:ncontr_t,imu) = dipole_(1:ncontr_t,1:ncontr_t)
     !
     !if (extF%matelem_threshold>0) then 
     !  !
     !  do i = 1,ncontr_t
     !    do j = 1,ncontr_t
     !       if ( abs( dipole_me(j,i,imu) )<extF%matelem_threshold) then 
     !          dipole_me(j,i,imu) = 0
     !       endif
     !    enddo
     !  enddo
     !  !
     !endif
     !
     if (job%verbose>=5) write(out,"(' done')")
     !
     if (job%IOextF_divide) then
       !
       read(chkptIO) buf4
       if (buf4/='extF') then
         write (out,"(' restore_vib_matrix_elements ',a,' has bogus footer: ',a)") job%kinetmat_file,buf4
         stop 'restore_vib_matrix_elements - bogus file format'
       end if
       !
       close(chkptIO,status='keep')
       !
     endif
     !
   enddo
   !
   deallocate(dipole_)
   call ArrayStop('dipole_')
   !
   if (.not.job%IOextF_divide) then
     !
     read(chkptIO) buf20(1:18)
     if (buf20(1:18)/='End external field') then
       write (out,"(' restore_vib_matrix_elements ',a,' has bogus footer: ',a)") job%kinetmat_file,buf20(1:17)
       stop 'restore_vib_matrix_elements - bogus file format'
     end if
     !
     close(chkptIO,status='keep')
     !
   endif
   !
   if (job%verbose>=4)   write(out, '(a/)') '...done'
   !
 end subroutine restore_vib_matrix_elements
 !
 !
 ! Here we reaplce the vibrational (J=0) dipole moment elements
 !
 subroutine replace_vib_trans_dipole_moments
   !
   implicit none 
   !
   integer(ik)        :: chkptIO
   integer(ik)        :: ncontr_t,imu,nclasses,i1,i2
   logical            :: eof
   character(len=cl)  :: job_is
   real(rk)           :: f_tdm
   !
   nclasses = size(eigen(1)%cgamma)-1
   !
   if (nclasses/=1) then
     !
     write(out,"('replace_vib_trans_dipole_moments error: only works for J=0, nclasse=1 not ',i0)") nclasses
     stop 'replace_vib_trans_dipole_moments is only for model J=0'
     !
   endif
   !
   job_is ='replace J=0 dipole moment'
   call IOStart(trim(job_is),chkptIO)
   !
   if (job%verbose>=4) write(out, '(a, 1x, a)') 'chk = file', trim(job%tdm_file)//'.chk'
   !
   open(chkptIO,action='read',status='old',file=trim(job%tdm_file)//'.chk')
   !
   read(chkptIO,*) ncontr_t
   !
   if (bset_contr(1)%Maxcontracts/=ncontr_t) then
     write (out,"(' Vib TDM file ',a)") job%tdm_file
     write (out,"(' Actual and stored basis sizes at J=0 do not agree  ',2i8)") bset_contr(1)%Maxcontracts,ncontr_t
     stop 'replace_vib_trans_dipole_moments - in file - illegal ncontracts '
   end if
   !
   ncontr_t = bset_contr(1)%Maxcontracts
   !
   if (job%verbose>=5) write(out,"(/'replace_vib_trans_dipole_moments...: Number of elements: ',i8)") ncontr_t
   !
   if (job%verbose>=4) write(out,"(a)") "Replace dipole moments:"
   !
   eof = .false.
   !
   loop_tdm: do
     !
     if (eof) exit loop_tdm
     read(chkptIO,*,end=111) i1,i2,f_tdm
     !
     if (job%verbose>=4) write(out,"(2i8,2x,i3,2x,3e15.7,2x,e15.7,1x,'//')") i1,i2,imu,dipole_me(i1,i2,:),f_tdm
     !
     dipole_me(i1,i2,:) = dipole_me(i1,i2,:)*f_tdm
     dipole_me(i2,i1,:) = dipole_me(i2,i1,:)*f_tdm
     !
     cycle
     111 continue
       eof = .true.
     exit
     !
   enddo loop_tdm
   !
   if (job%verbose>=4)   write(out, '(a/)') '...done'
   !
 end subroutine replace_vib_trans_dipole_moments
 !
 !
 function cg(j0, k0, dj, dk)

 ! return values for some CG coeeficients
 ! j0      and k0      are j" and k" (initial state)
 ! j0 + dj and k0 + dk are j' and k' (final state)
 ! compute <j0 k0, 1 dk|j0+dj k0+dk>

    integer(ik), intent(in) :: j0, k0, dj, dk
    real(rk)                :: cg
    real(rk)                :: j, m

    j = real(j0, kind = rk)
    m = real(k0 + dk, kind = rk)

    cg = 0.0_rk

    if (dj ==  0 .and. dk ==  0) cg = m / sqrt(j * (j + 1.0_rk))

    if (dj ==  0 .and. dk ==  1) cg = -sqrt((j + m) * (j - m + 1.0_rk) / (2.0_rk * j * (j + 1.0_rk)))

    if (dj ==  0 .and. dk == -1) cg = sqrt((j - m) * (j + m + 1.0_rk) / (2.0_rk * j * (j + 1.0_rk)))

    if (dj ==  1 .and. dk ==  0) cg = sqrt((j - m + 1.0_rk) * (j + m + 1.0_rk) / ((2.0_rk * j + 1.0_rk) * (j + 1.0_rk)))

    if (dj ==  1 .and. dk ==  1) cg = sqrt((j + m) * (j + m + 1.0_rk) / ((2.0_rk * j + 1.0_rk) * (2.0_rk * j + 2.0_rk)))

    if (dj ==  1 .and. dk == -1) cg = sqrt((j - m) * (j - m + 1.0_rk) / ((2.0_rk * j + 1.0_rk) * (2.0_rk * j + 2.0_rk)))

    if (dj == -1 .and. dk ==  0) cg =-sqrt((j - m) * (j + m) / (j * (2.0_rk * j + 1.0_rk)))

    if (dj == -1 .and. dk ==  1) cg = sqrt((j - m) * (j - m + 1.0_rk) / (2.0_rk * j * (2.0_rk * j + 1.0_rk)))

    if (dj == -1 .and. dk == -1) cg = sqrt((j + m + 1.0_rk) * (j + m) / (2.0_rk * j * (2.0_rk * j + 1.0_rk)))

    return

 end function cg


 !
 ! Analysis density calculations for vectors in symmetrised representation
 !
 subroutine DM_density_symmvec(Jval)
    implicit none
    !
    integer(ik),intent(in)  :: Jval(:)
    !
    integer(ik)    :: nJ,dimenmax
    integer(ik)    :: ilevelI, ilevelF, ndegI, ndegF, idegI, idegF, irec, idimen, nsizeF,nsizeI
    integer(ik)    :: irep
    integer(ik)    :: nmodes,info,indI,indF,jI,jF,Nrepresen
    integer(ik)    :: igammaI,igammaF,quantaI(0:molec%nmodes),quantaF(0:molec%nmodes),normalI(0:molec%nmodes)
    integer(ik)    :: dimenI,dimenF,unitI,unitF,jmax,kI,kF
    real(rk)       :: energyI, energyF, f_t,ener_
    logical        :: passed
    !
    integer(ik), allocatable :: icoeffI(:),istored(:),isave(:),isaved(:)
    integer(ik), allocatable :: nlevelsG(:,:),ram_size(:,:),ilevelsG(:,:),iram(:,:)
    real(rk),    allocatable :: vecI(:), vecF(:),vec(:),vec_(:),rot_density(:,:,:,:)
    complex(rk), allocatable :: dens_coord(:,:,:,:)
    real(rk),allocatable     :: density(:,:,:),Jrot_mat(:,:,:)
    integer(ik), pointer :: Jeigenvec_unit(:,:)
    type(DmatT),pointer  :: vec_ram(:,:)
    type(DkmatT),pointer :: ijterm(:)
    !
    integer(ik)  :: jind,nlevels,igamma,nclasses,nsize_need,nlevelI,dimenmax_,nsizemax,itauI,itauF,termI,termF
    !
    integer(ik)  :: ram_size_
    !
    integer(ik)  :: iram_,Nterms,iterm
    !
    integer(hik) :: matsize,ram_size_max,dimenmax_ram,dimenmax_ram_
    integer(ik)  :: ilist,nlist,imode,kmode,jmode,dens_size(3),ilevel_analyse(100),Nstored,i,k,i1,i2,i3,ipp,ip,npp,itheta,iphi
    real(rk)     :: rot(3),dens_,dtheta,dphi
    !
    integer(ik)  :: irootI,irow,ib,nelem,ielem,isrooti
    real(rk)     :: dtemp0 
    !
    double precision,allocatable  :: rw(:)
    double complex,allocatable  :: wd(:),rotmat(:,:),vl(:,:),vr(:,:),w(:),angular_m(:,:),crot_density(:,:,:,:)
    character(len=1)        :: jobvl,jobvr
    integer(ik)             :: cnu_i(0:molec%nmodes),ktau,tauI
    character(len=cl)       :: my_fmt !format for I/O specification
    character(len=wl)       :: my_fmt1 !format for I/O specification
    !
    call TimerStart('Density analysis')
    !
    nmodes = molec%nmodes
    nclasses = size(eigen(1)%cgamma)-1
    Nrepresen = sym%Nrepresen
    !
    nJ = size(Jval)
    !
    ! Prepare the list of units with the stored eigenvectors
    !
    allocate(Jeigenvec_unit(nJ,sym%Nrepresen), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: Jeigenvec_unit - out of memory'
    !
    allocate(nlevelsG(nJ,sym%Nrepresen),ilevelsG(nJ,sym%Nrepresen),ram_size(nJ,sym%Nrepresen),iram(nJ,sym%Nrepresen),stat = info)
    call ArrayStart('nlevelsG',info,size(nlevelsG),kind(nlevelsG))
    call ArrayStart('ilevelsG',info,size(ilevelsG),kind(ilevelsG))
    call ArrayStart('ram_size',info,size(ram_size),kind(ram_size))
    call ArrayStart('iram',info,size(iram),kind(iram))
    !
    do jind=1,nJ
      do igamma = 1,sym%Nrepresen
        Jeigenvec_unit(jind,igamma) = TReigenvec_unit(jind,Jval,igamma)
      enddo
    enddo
    !
    ! maximal size of basis functions 
    !
    dimenmax = 0
    nsizemax = 0
    !
    !loop over J quantities
    !
    do jind = 1, nJ
        !
        ! Estimate the maximal size of the basis functions. 
        !
        dimenmax = max(dimenmax,bset_contr(jind)%Maxcontracts)
        nsizemax = max(nsizemax,maxval(bset_contr(jind)%nsize(:)))
        !
    enddo 
    !
    jmax = Jval(nJ)
    !
    nlevelI = 0
    !
    !number of initial states
    !
    nlevels = Neigenlevels
    !
    allocate (ijterm(nJ),stat=info)
    !
    do jind=1,nJ
      !
      idimen = bset_contr(jind)%Maxcontracts
      !
      allocate (ijterm(jind)%kmat(bset_contr(jind)%Maxsymcoeffs,sym%Nrepresen),stat=info)
      call ArrayStart('ijterm',info,size(ijterm(jind)%kmat),kind(ijterm(jind)%kmat))
      !
      do igammaI = 1,sym%Nrepresen
         !
         Nterms = 0 
         !
         do iterm = 1,bset_contr(jind)%Maxsymcoeffs
           !
           ijterm(jind)%kmat(iterm,igammaI) = Nterms
           !
           Nterms = Nterms + bset_contr(jind)%irr(igammaI)%N(iterm) 
           !
         enddo
         !
      enddo
    enddo
    !
    ! loop over initial states
    !
    ener_ = 0
    !
    matsize = 0
    !
    nsizemax = 1
    nlevelsG = 0
    !
    ilist = 1 
    do while (ilist<size(analysis%dens_list).and.analysis%dens_list(ilist)/=-1)
      ilist = ilist + 1 
    enddo
    nlist = ilist-4
    f_t = 0
    !
    do ilevelF = 1, nlevels
       !
       jF = eigen(ilevelF)%jval
       indF = eigen(ilevelF)%jind
       ndegF   = eigen(ilevelF)%ndeg
       !
       !energy and and quanta of the final state
       !
       energyF = eigen(ilevelF)%energy
       igammaF = eigen(ilevelF)%igamma        
       quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes)
       !
       passed = .false.
       loop_ilist : do ilist=1,nlist
         !
         if (jF==analysis%j_list(ilist).and.igammaF==analysis%sym_list(ilist).and.&
             eigen(ilevelF)%irec(1)==analysis%dens_list(ilist)) then 
           passed = .true.
           exit loop_ilist
         endif 
         !
       enddo loop_ilist
       !
       if (.not.passed) cycle
       !
       nlevelsG(indF,igammaF) = nlevelsG(indF,igammaF) + 1
       !
       nsize_need = max(bset_contr(nJ)%nsize(igammaF),1)
       !
       f_t = f_t + real(nsize_need,rk)*( real(rk,rk) )/1024_rk**3
       !
       matsize = matsize + nsize_need
       nsizemax = max(nsizemax,nsize_need)
       !
    enddo
    !
    write(my_fmt,'(a,i0,a)') "(a,i4,a,",Nrepresen,"i8)"
    !
    do jind = 1, nJ
      jI = Jval(jind) 
      write(out,my_fmt) 'Number of states for J = ',jI,' each symm = ',nlevelsG(jind,:)
      !
    enddo
    !
    dimenmax_ram = max(nsizemax,1)
    !
    write(out,"('Max size of vectors: sym = ',i0,' prim =  ',i0)") dimenmax_ram,matsize
    !
    allocate(vec_ram(nJ,sym%Nrepresen),stat = info)
    call ArrayStart('vec_ram',0,1,rk)
    !
    ! estimate the memory requirements 
    !
    if (job%verbose>=4) write(out,"(a,f12.5,a,i0,' per each vector = ',f12.5'Gb')") &
        'All vectors will require approx ',f_t,' Gb; number of vectors = ',nlevels,f_t/real(nlevels,rk)
    !
    ! memory per vector
    !
    f_t = f_t/real(nlevels,rk)
    !
    ram_size_max = max(0,min( nlevels,int((memory_limit-memory_now)/f_t,hik)))
    !
    if (job%verbose>=4) write(out,"(a,f12.5,a,f12.5,'Gb; number of vectors to be kept in RAM =  ',i0/)") &
                               'Memory needed = ',f_t*real(nlevels,rk),'Gb;  memory  available = ',&
                               memory_limit-memory_now,ram_size_max
    !
    do jind = 1, nJ
       !
       !if (job%verbose>=4) write(out,"(' J    irep  N-RAM-irep   Nreps/Ntotal   ')") 
       !
       do irep = 1,Nrepresen
          !
          f_t = real( nlevelsG(jind,irep),rk )/real(nlevels, rk )
          ram_size_ = int(f_t*ram_size_max,ik)
          !
          ram_size_ = min(ram_size_,nlevelsG(jind,irep))
          !
          if (abs(ram_size_-nlevelsG(jind,irep))<=1) ram_size_=nlevelsG(jind,irep)
          !
          ram_size(jind,irep)  = ram_size_
          !
          if (ram_size_>0) then
            !
            if (job%verbose>=5) write(out,"(i3,4x,i2,4x,i0,4x,f14.5)") jval(jind),irep,ram_size_,f_t
            !
            !dimenmax_ram_ = max(min(int(bset_contr(jind)%nsize(irep)*intensity%factor),bset_contr(jind)%nsize(irep)),1)
            !
            dimenmax_ram_ = max(bset_contr(jind)%nsize(irep),1)
            !
            matsize = int(ram_size_,hik)*int(dimenmax_ram_,hik)
            !
            if (job%verbose>=4) write(out,"('Allocate (J=',i3,',irep=',i2,') matrix   ',i8,' x',i8,' = ',i0,', ',&
                                f14.4,'Gb')") jval(jind),irep,ram_size_,dimenmax_ram_,matsize,&
                                real(matsize*8,rk)/1024.0_rk**3
            !
            allocate(vec_ram(jind,irep)%mat(dimenmax_ram_,ram_size_),stat=info)
            !
            call ArrayStart('vec_ram',info,1,rk,matsize)
            !
          endif
          !
       enddo
    enddo
    !
    allocate(istored(nlevels),stat=info)
    !
    call ArrayStart('istored',info,size(istored),kind(istored))
    !
    istored = 0
    !
    allocate(isaved(nlevels),stat=info)
    !
    call ArrayStart('isaved',info,size(isaved),kind(isaved))
    !
    isaved = 0
    !
    allocate(isave(nJ),stat=info)
    !
    if (job%verbose>=5) call MemoryReport
    !
    if (job%verbose>=6) write(out,"(/'Pre-screening, compacting and storing...')")
    !
    ! prescreen all eigenfunctions, compact and store on the disk
    !
    call TimerStart('Distributing memory')
    !
    iram = 0
    isave = 0
    ilevelsG = 0
    !
    Nstored  = 0
    !
    do ilevelF = 1, nlevels
      !
      if (job%verbose>=6.and.mod(ilevelF,nlevels/min(500,ilevelF))==0 ) write(out,"('ilevel = ',i0)") ilevelF
      !
      jF = eigen(ilevelF)%jval 
      !
      !energy and and quanta of the final state
      !
      energyF = eigen(ilevelF)%energy
      igammaF = eigen(ilevelF)%igamma        
      quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
      !
      passed = .false.
      loop_ilist1 : do ilist=1,nlist
        !
        if (jF==analysis%j_list(ilist).and.igammaF==analysis%sym_list(ilist).and.&
            eigen(ilevelF)%irec(1)==analysis%dens_list(ilist)) then 
          passed = .true.
          exit loop_ilist1
        endif 
        !
      enddo loop_ilist1
      !
      if (.not.passed) cycle
      !
      Nstored  = Nstored + 1
      !
      ilevel_analyse(Nstored) = ilevelF
      !
      indF = eigen(ilevelF)%jind
      !
      ilevelsG(indF,igammaF) = ilevelsG(indF,igammaF) + 1
      !
      ndegF   = eigen(ilevelF)%ndeg
      dimenF = bset_contr(indF)%Maxcontracts
      nsizeF = bset_contr(indF)%nsize(igammaF)
      !
      unitF = Jeigenvec_unit(indF,igammaF)
      !
      if (nsizeF<=dimenmax_ram.and.iram(indF,igammaF)<ram_size(indF,igammaF)) then
        !
        iram(indF,igammaF) = iram(indF,igammaF) + 1
        !
        istored(ilevelF) = iram(indF,igammaF)
        !
        irec = eigen(ilevelF)%irec(1)
        !
        read(unitF,rec=irec) vec_ram(indF,igammaF)%mat(1:nsizeF,iram(indF,igammaF))
        !
      else
        !
        isave(indF) = isave(indF) + 1
        !
        isaved(ilevelF) = eigen(ilevelF)%irec(1)
        !
      endif
      !
    enddo
    !
    !deallocate(vecI,vecF)
    !call ArrayStop('intensity-vec')
    !
    call TimerStop('Distributing memory')
    !
    write(out,"(/'...done!')")
    !
    dimenmax_ = 1
    do jind = 1, nJ
       do irep = 1,Nrepresen
         dimenmax_ = max(bset_contr(jind)%nsize(irep),dimenmax_)
       enddo
    enddo
    !
    allocate(vecI(dimenmax_),vecF(dimenmax_),vec(dimenmax),vec_(dimenmax), stat = info)
    !
    call ArrayStart('density-vectors',info,size(vecI),kind(vecI))
    call ArrayStart('density-vectors',info,size(vecF),kind(vecF))
    call ArrayStart('density-vectors',info,size(vec),kind(vec))
    !
    allocate(icoeffI(dimenmax), stat = info)
    call ArrayStart('density-vectors',info,size(icoeffI),kind(icoeffI))
    !
    !write(out,"('Number of states for each symm = ',<Nrepresen>i8)") nlevelsG(:)
    !
    write(my_fmt,'(a,i0,a)') "(a,i4,a,",Nrepresen,"i8)"
    !
    do jind = 1, nJ
      jI = Jval(jind) 
      write(out,my_fmt) 'Number of states for ',jI,' each symm = ',nlevelsG(jind,:)
    enddo
    !
    if (all(nlevelsG(:,:)<1)) then 
      !
      write(out,"('DM_density_symmvec error: number of states to analyse is zero!')") 
      return 
      !
    endif
    !
    !  The matrix where some of the eigenvectors will be stored
    !
    if (job%verbose>=5) call MemoryReport
    !
    if (analysis%reduced_density) then
      !
      imode = min(analysis%dens_list(nlist+1),nclasses) ; jmode = min(analysis%dens_list(nlist+2),nclasses)
      kmode = min(analysis%dens_list(nlist+3),nclasses)
      !
      dens_size(:) = 1
      do jind=1,nJ
        !
        dens_size(1) = max(maxval(bset_contr(jind)%contractive_space(imode,:),dim=1),dens_size(1))
        dens_size(2) = max(maxval(bset_contr(jind)%contractive_space(jmode,:),dim=1),dens_size(2))
        dens_size(3) = max(maxval(bset_contr(jind)%contractive_space(kmode,:),dim=1),dens_size(3))
        !
      enddo
      !
      if (any(dens_size(:)<1)) then 
        !
        write(out,"('DM_density_symmvec error: dens_size=0 :',3i8)") dens_size
        stop 'DM_density_symmvec error: dens_size=0'
        !
      endif
      !
      if (job%verbose>=3) write(out,"('Allocate density matrix of dims = ',3i9)") dens_size(:)
      !
      allocate (density(dens_size(1),dens_size(2),dens_size(3)),stat=info)
      call ArrayStart('density',info,size(density),kind(density))
      !
      density = 0
      !
    endif
    !
    if (analysis%rotation_matrix) then 
      !
      ! count all states incl. degeneracy
      npp = 0 
      do i =1,Nstored
        ilevelI = ilevel_analyse(i)
        ndegI   = eigen(ilevelI)%ndeg
        npp = npp + ndegI 
      enddo
      allocate(Jrot_mat(3,npp,npp),stat=info)
      call ArrayStart('Jrot_mat',info,size(Jrot_mat),kind(Jrot_mat))
      !
      if (analysis%rotation_density) then
        allocate(rot_density(2*jmax+1,2*jmax+1,npp,npp),stat=info)
        call ArrayStart('rot_density',info,size(rot_density),kind(rot_density))
        allocate(crot_density(2*jmax+1,2*jmax+1,npp,npp),stat=info)
        call ArrayStart('rot_density',info,size(crot_density),kind(crot_density))
      endif 
      !
    endif
    !
    ! loop over initial states
    !
    ipp = 0 
    !
    Ilevels_loop: do i =1,Nstored
      !
      ilevelI = ilevel_analyse(i)
      !
      ! start measuring time per line
      !
      indI = eigen(ilevelI)%jind
      !
      !dimension of the bases for the initial states
      !
      dimenI = bset_contr(indI)%Maxcontracts
      !
      !energy, quanta, and gedeneracy order of the initial state
      !
      jI = eigen(ilevelI)%jval
      energyI = eigen(ilevelI)%energy
      igammaI  = eigen(ilevelI)%igamma
      quantaI(0:nmodes) = eigen(ilevelI)%quanta(0:nmodes)
      normalI(0:nmodes) = eigen(ilevelI)%normal(0:nmodes)
      ndegI   = eigen(ilevelI)%ndeg
      nsizeI = bset_contr(indI)%nsize(igammaI)
      !
      ! where the eigenvector is stored 
      !
      unitI = Jeigenvec_unit(indI,igammaI) 
      !
      !read eigenvector of initial state
      !
      !
      if (istored(ilevelI)==0) then
          !
          call TimerStart('Reading eigenfuncs')
          !
          irec = isaved(ilevelI)
          !
          read(unitI,rec=irec) vecI(1:nsizeI)
          !
          call TimerStop('Reading eigenfuncs')
          !
      else
          !
          iram_ = istored(ilevelI)
          !
          vecI(1:nsizeI) = vec_ram(indI,igammaI)%mat(1:nsizeI,iram_)
          !
      endif
      !
      ndegF  = ndegI
      !
      do idegI = 1, ndegI
         !
         call convert_symvector_to_contrvector(indI,dimenI,igammaI,idegI,ijterm(indI)%kmat(:,:),vecI,vec)
         !
         if (analysis%reduced_density) then
           !
           call do_reduced_density(indI,vec,imode,jmode,kmode,density)
           !
           !print*, vec
           !
           write(my_fmt1,'(a,i0,a,i0,a,i0,a)') "(i4,1x, a3,2x,f13.6,1x,a1,a3,a1,i3,a1,1x,a1,",nclasses,"a3,a1,",&
                                               nmodes,"(1x, i3),a1,2x,a1,",nmodes,"(1x, i3),a1)"
           !
           !write(out,"(i4,1x, a3,2x,f13.6,1x,a1,a3,a1,i3,a1,1x,a1,<nclasses>a3,a1,<nmodes>(1x, i3),a1,2x,a1,<nmodes>(1x, i3),a1)") & 
           !
           write(out,my_fmt1) & 
                            jI,sym%label(igammaI),&
                            energyI-intensity%ZPE,'(',&
                            eigen(ilevelI)%cgamma(0),';',eigen(ilevelI)%krot,')',&
                            '(',eigen(ilevelI)%cgamma(1:nclasses),';',eigen(ilevelI)%quanta(1:nmodes),')', &
                            '(',normalI(1:nmodes),')'
           !
           do i1 = 1,dens_size(1)
             do i2 = 1,dens_size(2)
               do i3 = 1,dens_size(3)
                 !
                 if (density(i1,i2,i3)>small_) then
                   !
                   write(out,"(4i8,x,f20.14)") & 
                            i,i1,i2,i3,density(i1,i2,i3)
                 endif
                 !
               enddo
             enddo
           enddo
           !
         endif
         !
         if (analysis%print_vector) then
            !
            write(my_fmt1,'(a,i0,a)') "(2x,i4,i7,e16.8,3x,a1,3i3,1x,i4,i2,a1,1x,a1,",Nclasses,"(i3),a1)"
            !
            do irootI = 1, dimenI
                 !
                 irow = bset_contr(indI)%icontr2icase(irootI,1)
                 ib   = bset_contr(indI)%icontr2icase(irootI,2)
                 !
                 cnu_i(0:Nclasses) = bset_contr(indI)%contractive_space(0:Nclasses, irow)
                 !
                 ktau = bset_contr(indI)%ktau(irootI)
                 tauI  = mod(ktau,2_ik)
                 kI = bset_contr(indI)%k(irootI)
                 !
                 !ndeg = bset_contr(jind)%index_deg(irow)%size1
                 !
                 if (abs(vec(irootI))>analysis%threshold) then  
                   !
                   write(out,my_fmt1) & 
                              igammaI,irootI,&
                              vec(irootI),"(", &
                              jI,kI,tauI,irow,ib,")", &
                              "(",cnu_i(1:Nclasses),")"
                   !
                 endif
                 !
            end do
            !
         endif 
         !
         if (analysis%rotation_matrix) then 
            !
            ipp = ipp + 1
            !
            ! store this vector to be used for plotting the rotational density 
            !if (analysis%rotation_density) clustervec(:,ipp)  = vec
            !
            ip = 0 
            !
            do k =1,Nstored
              !
              ilevelF = ilevel_analyse(k)
              jF = eigen(ilevelF)%jval
              indF = eigen(ilevelF)%jind
              !
              if (jI/=jF) then
                stop 'rotation_matrix: Something wrong with Ji/jF'
              endif
              !
              energyF = eigen(ilevelF)%energy
              igammaF = eigen(ilevelF)%igamma        
              quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
              normalI(0:nmodes) = eigen(ilevelI)%normal(0:nmodes)
              ndegF   = eigen(ilevelF)%ndeg
              nsizeF = bset_contr(indF)%nsize(igammaF)
              dimenF = bset_contr(indF)%Maxcontracts
              !
              ! where the eigenvector is stored 
              !
              unitF = Jeigenvec_unit(indF,igammaF) 
              !
              !read eigenvector of initial state
              !
              if (istored(ilevelF)==0) then
                  !
                  irec = isaved(ilevelF)
                  !
                  read(unitF,rec=irec) vecF(1:nsizeF)
                  !
              else
                  !
                  iram_ = istored(ilevelF)
                  !
                  vecF(1:nsizeF) = vec_ram(indF,igammaF)%mat(1:nsizeF,iram_)
                  !
              endif
              !
              do idegF = 1, ndegF
                  ip = ip + 1
                  !
                  call convert_symvector_to_contrvector(indF,dimenF,igammaF,idegF,ijterm(indF)%kmat(:,:),vecF,vec_)
                  !
                  call do_Angular_momentum_average(jI,indI,vec,vec_,jrot_mat(:,ip,ipp))
                  !
                  if (analysis%rotation_density) call do_rotational_density(indI,vec,vec_,rot_density(:,:,ip,ipp))
                  !
              enddo
              !
            enddo
            !
         endif 
         !
      enddo
      !
    end do Ilevels_loop
    !
    if (analysis%reduced_density) then
      deallocate (density)
      call ArrayStop('density')
    endif
    !
    if (analysis%rotation_matrix) then
       !
       if (job%verbose>=2) then
         ! 
         write(out,"(/'   Expectation values of the rotational (angular momentum) matrix: ')")
         !
         do ipp = 1,npp
           do ip = 1,npp
              write(out,"(3x,3f20.14,3i8)") jrot_mat(:,ipp,ip),ipp,ip
           enddo
         enddo
         !
       endif
       !
       if (analysis%rotation_density) then
         !
         ! rotation of the ang. momentum by the user-defined angles to the axis with the maximal projection
         !
         rot(1) = sin( analysis%res%theta )*cos( analysis%res%phi )
         rot(2) = sin( analysis%res%theta )*sin( analysis%res%phi )
         rot(3) = cos( analysis%res%theta )
         !
         allocate(angular_m(npp,npp),stat=info)
         !
         ! angular momentum in the rotated frame
         !
         do ipp = 1,npp
            do ip = 1,npp
               angular_m(ipp,ip) = cmplx( 0.0_rk,-sum(rot(:)*jrot_mat(:,ipp,ip)),rk )
            enddo
         enddo
         !
         jobvr='V'
         jobvl='V'
         !
         allocate(rw(npp),stat=info)
         deallocate(rw,stat=info)
         !
         allocate(wd(2*npp+5),vl(npp,npp),vr(npp,npp),w(npp),rotmat(npp,npp),stat=info)
         allocate(rw(2*npp+5),stat=info)
         call ArrayStart('zgeev-mat',info,size(wd),kind(wd))
         call ArrayStart('zgeev-mat',info,size(vr),kind(vr))
         call ArrayStart('zgeev-mat',info,size(vl),kind(vl))
         call ArrayStart('zgeev-mat',info,size(rotmat),kind(rotmat))
         call ArrayStart('zgeev-mat',info,size(angular_m),kind(angular_m))
         call ArrayStart('zgeev-mat',info,size(rw),kind(rw))
         !
         rotmat = angular_m
         !
         ! diagonalize the angular momenum matrix 
         !
         call ZGEEV(JOBVL,JOBVR,npp,rotmat,npp,W,VL,npp,VR,npp,wd,size(wd),rw,info)
         !
         if (info /= 0) stop 'ZGEEV angular_m error'
         !
         deallocate(rw,stat=info)
         !
         if (job%verbose>=2) then 
           !
           write(out,"(/'   Expectation values of the angular momentum in the rotated frame: ')")
           !
           do ipp = 1,npp
              write(out,"(3x,i8,1x,2f16.6)") ipp,w(ipp)
           enddo
           !
         endif
         !
         ! check the solution and compare with the eigenvalues
         !
         !rotmat = matmul(conjg(transpose(Vr)),matmul(angular_m,vr))
         !rotmat = matmul(conjg(transpose(Vr)),vr)
         !
         ! now we have a transformation to maximize the projecion of the angular momentum 
         ! which we will apply to the cluster eigenfunctions
         !
         if (job%verbose>=3) write(out,"(/'   Rotational density in the rotated frame...')")
         !
         do kI = 0,jmax
           do itauI = 0,min(kI,1)
             !
             termI = 2*kI+itauI ; if (kI==0) termI = termI + 1
             !
             do kF = 0,jmax
               do itauF = 0,1
                 !
                 termF = 2*kF+itauF ; if (kF==0) termF = termF + 1
                 !
                 ! rotational density in the rotated (primitve cluster) frame
                 !
                 crot_density(termI,termF,:,:) = matmul(conjg(transpose(Vr)),matmul(cmplx(rot_density(termI,termF,:,:),0.0_rk),vr))
                 !
               enddo
             enddo
           enddo
         enddo
         !
         if (job%verbose>=3) write(out,"(/'   ...done!')")
         !
         ! transform to the primitive |J,k,m=0> (k=-J..J) and coordinate (theta,phi) representations 
         !
         allocate(dens_coord(0:analysis%res%ntheta,0:analysis%res%nphi,2*jmax+1,2*jmax+1),stat=info)
         call ArrayStart('dens_coord',info,size(dens_coord),kind(dens_coord))
         !
         if (job%verbose>=3) write(out,"(/'   Generate coordinate representation... ')")
         !
         call Transform_rotdens_to_coord_repres(jmax,dens_coord)
         !
         if (job%verbose>=3) write(out,"(/'   ...done!')")
         !
         write(out,"(' Reduced rotational density of the cluster states:')") 
         !
         do ip=1,npp
           !
           write(out,"(' state =  ',i6)") ip
           !
           dtheta = pi/(analysis%res%ntheta)
           dphi   = 2.0_rk*pi/(analysis%res%nphi)
           !
           do itheta = 0,analysis%res%ntheta
             do iphi = 0,analysis%res%nphi
                 !
                 ! Rotational probabily density value at the instanteneous (theta,phi) coordinate in the primitive cluster state (PCS) representaion
                 !
                 dens_ = 0 
                 do i1 = 1,2*jmax+1
                   do i2 = 1,2*jmax+1
                     !
                     !dens_ =sum(crot_density(:,:,ip,ip)*dens_coord(itheta,iphi,:,:))
                     !
                     dens_ = dens_ + real(crot_density(i1,i2,ip,ip)*dens_coord(itheta,iphi,i1,i2),rk)
                     !
                   enddo
                 enddo
                 !
                 !write(out,"(2f9.2,3x,e18.7)") itheta*dtheta*rad,iphi*dphi*rad,dens_
                 !
                 write(out,"(2f9.2,3x,e18.7)") itheta*dtheta*rad,iphi*dphi*rad,abs(dens_)
                 !
             enddo
           enddo
         enddo
         !
         !clustervec_ = matmul(clustervec,vr) 
         !
         deallocate(wd,vl,vr,w,stat=info)
         deallocate(rotmat,angular_m,stat=info)
         call ArrayStop('zgeev-mat')
         !
         deallocate(rot_density,crot_density)
         call ArrayStop('rot_density')
         !
         deallocate(dens_coord)
         call ArrayStop('dens_coord')
         !
       endif
       !
       deallocate(Jrot_mat)
       call ArrayStop('Jrot_mat')
    endif
    !
    deallocate(vecI, vecF, vec, vec_,icoeffI)
    call ArrayStop('density-vectors')
    !
    do irep = 1,Nrepresen
     do jind = 1,nJ
      !
      nullify(vec_ram(jind,irep)%mat)
      !
      if (associated(vec_ram(jind,irep)%mat)) deallocate(vec_ram(jind,irep)%mat)
      !
     enddo
    enddo
    !
    call ArrayStop('vec_ram')
    !
    do jind = 1,nJ
      !
      if (associated(ijterm(jind)%kmat)) deallocate(ijterm(jind)%kmat)
      !
    enddo
    call ArrayStop('ijterm')
    !
    deallocate(nlevelsG,ilevelsG,ram_size,iram)
    call ArrayStop('nlevelsG')
    call ArrayStop('ilevelsG')
    call ArrayStop('ram_size')
    call ArrayStop('iram')
    !
    deallocate(istored)
    call ArrayStop('istored')
    !
    deallocate(isaved)
    call ArrayStop('isaved')
    !
    call TimerStop('Density analysis')
    !
  end subroutine DM_density_symmvec


 !
 ! Electric dipole moment intensities and transition moments calculations 
 !
 subroutine dm_intensity(Jval)

    integer(ik),intent(in)  :: Jval(:)

    integer(ik)    :: nJ,dimenmax
    integer(ik)    :: ilevelI, ilevelF, ndegI, ndegF, idegI, idegF, irec, idimen, cdimen, cdimenmax
    integer(ik)    :: cdimenI(sym%Maxdegen),cdimenF(sym%Maxdegen),nlevelsG(sym%Nrepresen),ilevelsG(sym%Nrepresen),irep
    integer(ik)    :: nmodes,info,indI,indF,itransit,Ntransit,jI,jF,Nrepresen
    integer(ik)    :: igammaI,igammaF,quantaI(0:molec%nmodes),quantaF(0:molec%nmodes),normalF(0:molec%nmodes),&
                      normalI(0:molec%nmodes)
    integer(ik)    :: dimenI,dimenF,unitI,unitF,imu,jmax,kI,kF,unitO,unitC
    real(rk)       :: energyI, energyF, nu_if,linestr,linestr_deg(sym%Maxdegen,sym%Maxdegen),f_t,cfactor,ener_
    real(rk)       :: tm_deg(3,sym%Maxdegen,sym%Maxdegen)
    logical        :: passed,passed_

    integer(ik), allocatable :: icoeffI(:,:), icoeffF(:,:),Jfuncs_unit(:),Jindex_unit(:),cdimens(:,:),istored(:),&
                                isave(:),isaved(:,:),idimenmax(:)
    real(rk),    allocatable :: vecI(:,:), vecF(:,:)
    real(rk),allocatable     :: half_linestr(:,:,:,:),threej(:,:,:,:)
    integer(ik), pointer :: Jeigenvec_unit(:)
    type(Dmat3T)  :: vec_swap(sym%Nrepresen)
    type(DkmatT) :: cdimen_swap(sym%Nrepresen)
    type(DimatT) :: icoeff_swap(sym%Nrepresen)
    !
    integer(ik)  :: jind,nlevels,nformat,nclasses,dimenmax_swap,nlevelI,dimenmax_
    !
    integer(ik)  :: swap_size_,swap_size(sym%Nrepresen)
    !
    integer(ik)  :: iswap_,iswap(sym%Nrepresen),M,MP1
    !
    integer(ik)  :: igamma_pair(sym%Nrepresen),rec_len,irec_len
    integer(hik) :: matsize,swap_size_max
    !
    !
    character(len=1) :: branch
    character(len=cl):: unitfname,filename
    character(4)     :: jchar
    !
    real(rk)         :: boltz_fc, beta, intens_cm_molecule, A_coef_s_1, A_einst, absorption_int, dtemp
    character(len=cl):: my_fmt !format for I/O specification
    character(len=wl) :: my_fmt_a,my_fmt_tm
    !
    call TimerStart('Intensity calculations')
    !
    if (sym%maxdegen>2) then 
      !
      write(out,"('dm_intensity: this procedure has not been tested for the symmetries with degeneracies higher than 2')")
      !stop 'dm_intensity was not tested for present symmetry'
      !
    endif
    !
    nmodes = molec%nmodes
    nclasses = size(eigen(1)%cgamma)-1
    Nrepresen = sym%Nrepresen
    !
    beta = planck * vellgt / (boltz * intensity%temperature)
    intens_cm_mol  = 8.0d-36 * pi**3* avogno / (3.0_rk * planck * vellgt)
    intens_cm_molecule  = 8.0d-36 * pi**3 / (3.0_rk * planck * vellgt)
    A_coef_s_1     =64.0d-36 * pi**4  / (3.0_rk * planck)
    !
    nJ = size(Jval)
    !
    ! Prepare the list of units with the stored eigenvectors
    !
    allocate(Jeigenvec_unit(nJ), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: Jeigenvec_unit - out of memory'
    !
    do jind=1,nJ 
      Jeigenvec_unit(jind) = TReigenvec_unit(jind,Jval)
    enddo
    !
    ! maximal size of basis functions 
    !
    dimenmax = 0 
    !
    !loop over J quantities
    !
    do jind = 1, nJ
        !
        ! Estimate the maximal size of the basis functions. 
        !
        dimenmax = max(dimenmax,bset_contr(jind)%Maxcontracts)
        !
    enddo 
    !
    jmax = Jval(nJ)
    !
    ! Precompute the 3-j symbol values 
    !
    allocate(threej(0:jmax,0:jmax,-1:1,-1:1), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: threej - out of memory'
    !
    threej = 0 
    !
    do jI = 0,jmax
      do jF = max(jI-1,0),min(jI+1,jmax)
        do kI = 0,jI
          do kF = max(kI-1,0),min(kI+1,jF)
            !
            !threej(jI, kI, jF - jI, kF - kI) = cg(jI, kI, jF - jI, kF - kI) *                                                 &
            !               (-1.0_rk) ** (jI - 1 + kF) / sqrt( real(2*jF + 1,rk) )
            !
            !dp = three_j(jI, 1, jF, kI, kF - kI, -kF)
            !
            !if (abs(dp-threej(jI, kI, jF - jI, kF - kI))>sqrt(small_)) then 
            !  write(out,"('a problem with threej for jI,jF,kI,kF = ',4i4,3es16.8)") jI,jF,kI,kF,dp,threej(jI, kI, jF - jI, kF - kI),dp-threej(jI, kI, jF - jI, kF - kI)
            !  stop 'a problem with threej'
            !endif
            !
            threej(jI, kI, jF - jI, kF - kI) = three_j(jI, 1, jF, kI, kF - kI, -kF)
            !
          enddo
        enddo
      enddo
    enddo
    !
    ! in the case of the proper rotatonal symmetry the matrix elements of the 
    ! wigner matrix have to be computed using the correspoding symmetrization 
    ! transformation derived here. they will be used together with the 
    ! vibrational matrix elements of the dipole moment when evaluating 
    ! line strengthes using eigenvectos.  
    !
    if (job%rotsym_do) call contraced_rotational_dipole(nJ, jval, jmax, threej)
    !
    !allocate arrays for eigenvectors and coefficients
    !
    ! First we count transitions to be calculated. It will help us to keep track of the 
    ! calculation progress.
    !
    Ntransit = 0 
    !
    nlevelI = 0
    !
    !number of initial states
    !
    nlevels = Neigenlevels
    !
    ! For a given symmetry igamma with some gns(igamma) we find its counterpart jgamma/=igamma
    ! having the same gns(jgamma). We assume that there is only one such pair 
    ! in case of absorption or emission calcs. 
    !
    call find_igamma_pair(igamma_pair)
    !
    call TimerStart('Intens_Filter-1')
    !
    if (intensity%int_increm<1e6) write(out,"(/'Intensity increments:')")
    !
    ! loop over initial states
    !
    ener_ = 0
    !
    !omp parallel do private(ilevelI,jI,energyI,igammaI,quantaI,ilevelF,jF,energyF,igammaF,quantaF,passed) schedule(guided) reduction(+:Ntransit,nlevelI)
    do ilevelI = 1, nlevels
      !
      ! rotational quantum number 
      !
      jI = eigen(ilevelI)%jval
      !
      !energy energy and and quanta of the initial state
      !
      energyI = eigen(ilevelI)%energy
      !
      igammaI  = eigen(ilevelI)%igamma
      quantaI(0:nmodes) = eigen(ilevelI)%quanta(0:nmodes) 
      normalI(0:nmodes) = eigen(ilevelI)%normal(0:nmodes)
      !
      call energy_filter_lower(jI,energyI,quantaI,normalI(0),passed)
      !
      if (.not.passed) cycle
      !
      nlevelI = nlevelI + 1
      !
      !loop over final states
      !
      do ilevelF = 1, nlevels
         !
         jF = eigen(ilevelF)%jval 
         !
         !energy and and quanta of the final state
         !
         energyF = eigen(ilevelF)%energy
         igammaF = eigen(ilevelF)%igamma        
         quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
         !
         call intens_filter(jI,jF,energyI,energyF,igammaI,igammaF,quantaI,quantaF,igamma_pair,passed)
         !
         if (passed) then 
           !
           Ntransit = Ntransit + 1
           !
           if (mod(Ntransit,intensity%int_increm)==0.and.energyI>ener_) then
              ener_ = energyI
              write(out,"(i0,f16.4)") Ntransit,energyI-intensity%ZPE
           endif
           !
         endif 
         !
      enddo
      !
    enddo
    !omp end parallel do
    !
    call TimerStop('Intens_Filter-1')
    !
    !loop over final states -> count states for each symmetry
    !
    nlevelsG = 0
    !
    do ilevelF = 1, nlevels
       !
       jF = eigen(ilevelF)%jval
       indF = eigen(ilevelF)%jind
       ndegF   = eigen(ilevelF)%ndeg
       !
       !energy and and quanta of the final state
       !
       energyF = eigen(ilevelF)%energy
       igammaF = eigen(ilevelF)%igamma        
       quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
       normalF(0:nmodes) = eigen(ilevelF)%normal(0:nmodes)
       !
       call energy_filter_upper(jF,energyF,quantaF,normalF(0),passed)
       !
       call energy_filter_lower(jF,energyF,quantaF,normalF(0),passed_)
       !
       if (.not.passed.and..not.passed_) cycle
       !
       nlevelsG(igammaF) = nlevelsG(igammaF) + 1
       !
    enddo
    !
    write(my_fmt,'(a,i0,a)') "(a,",Nrepresen,"i8)"
    !
    write(out,my_fmt) 'Number of states for each symm = ',nlevelsG(:)
    !
    dimenmax_swap = max(min(int(dimenmax*intensity%factor),dimenmax),1)
    !
    write(out,"('Size of the reduced vector = ',i0,' total =  ',i0)") dimenmax_swap,dimenmax
    !
    call ArrayStart('vec_swap',0,1,rk)
    call ArrayStart('icoeff_swap',0,1,ik)
    !
    ! estimate the memory requirements 
    !
    f_t = 0
    do irep = 1,Nrepresen
       !
       f_t = f_t + real(nlevelsG(irep),rk)*real(sym%degen(irep),rk)*real(dimenmax_swap,rk)*(real(rk,rk)+&
                   real(ik,rk)+1.0_rk)/1024_rk**3
       !
       !f_t = f_t + real(sym%degen(irep),rk)*(real(dimenmax_swap,rk)+1.0_rk)*real(ik,rk)/1024_rk**3
       !
    enddo
    !
    matsize = int(sum( nlevelsG(:) ),hik)
    !
    if (job%verbose>=4) write(out,"('All vectors will require approx ',f12.5,' Gb; number of vectors = ',i0,&
                                  ' per each vector = ',f12.5'Gb')") f_t,matsize,f_t/real(matsize,rk)
    !
    ! memory per vector
    !
    f_t = f_t/real(matsize,rk)
    !
    !
    swap_size_max = max(0,min( matsize,int((memory_limit-memory_now)/f_t,hik)))
    !
    if (job%verbose>=4) write(out,"('Memory needed = ',f12.5,'Gb;  memory  available = ',f12.5,&
                                  'Gb number of vectors will be shared in RAM =  ',i0)") f_t*real(matsize,rk),&
                                  memory_limit-memory_now,swap_size_max
    !
    do irep = 1,Nrepresen
       !
       vec_swap(irep)%mat => null()
       icoeff_swap(irep)%imat => null()
       cdimen_swap(irep)%kmat => null()
       !
       f_t = real( nlevelsG(irep),rk )/real( sum( nlevelsG(:) ), rk )
       swap_size_ = int(f_t*swap_size_max,ik)
       !
       if (job%verbose>=4) write(out,"('irep = ',i2,';  total swap ',i0,' swap_irep =  ',i0,' f_t = ',f14.5)") &
                                       irep,swap_size_max,swap_size_,f_t
       !
       swap_size_ = min(swap_size_,nlevelsG(irep))
       !
       if (abs(swap_size_-nlevelsG(irep))<=1) swap_size_=nlevelsG(irep)
       !
       swap_size(irep)  = swap_size_
       !
       if (swap_size_>0) then 
         !
         matsize = int(swap_size_,hik)*int(sym%degen(irep),hik)*int(dimenmax_swap,hik)
         !
         if (job%verbose>=4) write(out,"('Allocate ',i2,'-th  swap matrix   ',i8,' x',i3,' x',i8,' = ',i0,', ',f14.4,'Gb')") &
                                   irep,swap_size_,sym%degen(irep),dimenmax_swap,matsize,real(matsize*8,rk)/1024.0_rk**3
         !
         allocate(vec_swap(irep)%mat(swap_size_,sym%degen(irep),dimenmax_swap),stat=info)
         !
         call ArrayStart('vec_swap',info,1,rk,matsize)
         !
         matsize = int(swap_size_,hik)*int(sym%degen(irep),hik)*int(dimenmax_swap,hik)+int(swap_size_,hik)*int(sym%degen(irep),hik)
         !
         if (job%verbose>=4) write(out,"('Allocate ',i2,'-th iswap matrix   ',i8,' x',i3,' x',i8,' = ',i0,', ',f14.4,'Gb')")&
                                        irep,swap_size_,sym%degen(irep),dimenmax_swap,matsize,real(matsize*4,rk)/1024.0_rk**3
         !
         allocate(icoeff_swap(irep)%imat(swap_size_,sym%degen(irep),dimenmax_swap),stat=info)
         allocate(cdimen_swap(irep)%kmat(swap_size_,sym%degen(irep)),stat=info)
         !
         call ArrayStart('icoeff_swap',info,1,ik,matsize)
         !
       endif
       !
       !f_t = real( nlevelsG(irep),rk )/real( sum( nlevelsG(:) ), rk )
       !swap_size_ = int(f_t*intensity%swap_size(1),ik)+1
       !
       !swap_size_ = max(min(swap_size_,nlevelsG(irep)),swap_size(irep))
       !
       !swap_size_i(irep)  = swap_size_
       !
    enddo
    !
    allocate(istored(nlevels),stat=info)
    !
    call ArrayStart('istored',info,size(istored),kind(istored))
    !
    istored = 0
    !
    allocate(isaved(nlevels,sym%Maxdegen),stat=info)
    !
    call ArrayStart('isaved',info,size(isaved),kind(isaved))
    !
    isaved = 0
    !
    allocate(isave(nJ),stat=info)
    allocate(Jindex_unit(nJ), stat = info)
    allocate(Jfuncs_unit(nJ), stat = info)
    allocate(cdimens(nlevels,sym%Maxdegen),stat=info)
    allocate(idimenmax(nJ),stat=info)
    !
    call ArrayStart('cdimens',info,size(cdimens),ik)
    !
    if (job%verbose>=5) call MemoryReport
    !
    !inquire(iolength=rec_len) beta
    !rec_len = rec_len*dimenmax_swap
    !
    !inquire(iolength=irec_len) cdimen
    !irec_len = irec_len*(dimenmax_swap+1)
    !
    do jind = 1, nJ
        !
        write(jchar, '(i4)') jval(jind)
        !
        write(unitfname, '(a, i4)') 'compacted eigenfuncs for ', jval(jind)
        call iostart(trim(unitfname), unitO)
        !
        Jfuncs_unit(jind) = unitO
        !
        write(unitfname, '(a, i4)') 'compacted coeffs for ', jval(jind)
        call iostart(trim(unitfname), unitC)
        !
        Jindex_unit(jind) = unitC
        !
        dimenmax_ = max(min(int(bset_contr(jind)%Maxcontracts*intensity%factor),bset_contr(jind)%Maxcontracts),1)
        !
        idimenmax(jind) = dimenmax_
        !
        !bset_contr(jind)%Maxcontracts
        !
        inquire(iolength=rec_len) beta
        rec_len = rec_len*dimenmax_
        !
        inquire(iolength=irec_len) cdimen
        irec_len = irec_len*(dimenmax_+1)
        !
        ! use sctratch 
        !
        if (trim(intensity%swap)=="NONE") then 
          !
          open(unitO,status='scratch',access='direct',recl=rec_len ) 
          open(unitC,status='scratch',access='direct',recl=irec_len) 
          !
        elseif(trim(intensity%swap)=="SAVE") then
          !
          if (jval(jind)==0.and.intensity%J(1)/=0) cycle
          !
          ! use the predefined file to store the compacted vectors
          !
          filename = trim(intensity%swap_file)//'_vect'//trim(adjustl(jchar))//'.tmp'
          open(unitO,access='direct',action = 'readwrite',status='replace',file=filename,recl=rec_len) 
          !
          filename = trim(intensity%swap_file)//'_coef'//trim(adjustl(jchar))//'.tmp'
          open(unitC,access='direct',action = 'readwrite',status='replace',file=filename,recl=irec_len) 
          !
        elseif(trim(intensity%swap)=="READ") then
          !
          ! use the predefined file to store the compacted vectors
          !
          filename = trim(intensity%swap_file)//'_vect'//trim(adjustl(jchar))//'.tmp'
          open(unitO,access='direct',action = 'read',status='old',file=filename,recl=rec_len) 
          !
          filename = trim(intensity%swap_file)//'_coef'//trim(adjustl(jchar))//'.tmp'
          open(unitC,access='direct',action = 'read',status='old',file=filename,recl=irec_len) 
          !
        endif
        !
    enddo 
    !
    if (job%verbose>=4) write(out,"(/'Pre-screening, compacting and storing...')")
    !
    ! prescreen all eigenfunctions, compact and store on the disk
    !
    call TimerStart('Prescreening eigenfuncs')
    !
    iswap = 0
    isave = 0
    ilevelsG = 0
	cdimenmax = 0
	cfactor = 0
    !
    allocate(vecI(sym%Maxdegen,dimenmax_swap),vecF(sym%Maxdegen,dimenmax), stat = info)
    if (info/=0)  stop 'vecI,vecF,icoeffF - out of memory'
    !
    allocate(icoeffF(sym%Maxdegen,dimenmax), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error:  icoeffF - out of memory'
    !
    if (trim(intensity%swap)=="SAVE".or.trim(intensity%swap)=="NONE") then 
       !omp parallel private(vecI,vecF,icoeffF,info) shared(cdimens) 
       !
       !omp do private(ilevelF,jF,energyF,igammaF,quantaF,indF,ndegF,dimenF,unitF,unitO,unitC,idegF,irec,cdimen,idimen) schedule(guided)
       do ilevelF = 1, nlevels
         !
         if (job%verbose>=5.and.mod(ilevelF,nlevels/min(500,ilevelF))==0 ) write(out,"('ilevel = ',i0)") ilevelF
         !
         jF = eigen(ilevelF)%jval
         !
         if (trim(intensity%swap)=="SAVE".and.jF==0.and.intensity%J(1)/=0) cycle
         !
         !energy and and quanta of the final state
         !
         energyF = eigen(ilevelF)%energy
         igammaF = eigen(ilevelF)%igamma        
         quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes)
         !
         indF = eigen(ilevelF)%jind
         !
         ndegF   = eigen(ilevelF)%ndeg
         dimenF = bset_contr(indF)%Maxcontracts
         !
         unitF = Jeigenvec_unit(indF) 
         unitO = Jfuncs_unit(indF)
         unitC = Jindex_unit(indF)
         !
         ! prescreen all eigenfunctions, compact and store on the disk
         !
         do idegF = 1, ndegF
           !
           !call TimerStart('Reading eigenfuncs-0')
           !
           irec = eigen(ilevelF)%irec(idegF)
           read(unitF, rec = irec) vecF(idegF,1:dimenF)
           !
           !call TimerStop('Reading eigenfuncs-0')
           !
           ! delete smallest basis contributions to the eigenvector of initial state
           !
           !call TimerStart('Compacting eigenfuncs')
           !
           cdimen = 0 
           icoeffF(idegF,:) = 0
           !
           do idimen = 1, dimenF
              if (abs(vecF(idegF,idimen)) > intensity%threshold%coeff) then
                 cdimen = cdimen + 1
                 icoeffF(idegF,cdimen) = idimen
                 vecI(idegF,cdimen) = vecF(idegF,idimen)
              end if
           end do
		   !
		   cdimenmax = max(cdimen,cdimenmax)
		   cfactor = max(real(cdimen)/real(idimenmax(indF)),cfactor)
           !
           cdimens(ilevelF,idegF) = cdimen
           !
           if (cdimen>idimenmax(indF)) then 
              !
              write(out,"(a,i0,' vs ',i0,a,f12.4)") 'intens: dimenmax_swap is too small (',&
                         idimenmax(indF),cdimen,') reduce the coeff value or increase swap_factor to minimum',&
                         real(cdimen,rk)/real(bset_contr(jind)%Maxcontracts,rk)
			  write(out,"('dimenmax_ = ',i0,' intensity%factor = ',i0)")  dimenmax_,intensity%factor
              stop 'intens: increase swap_factor!'
              !
           endif
           !
           write(unitO,rec=irec) vecI(idegF,1:cdimen)
           write(unitC,rec=irec) cdimen,icoeffF(idegF,1:cdimen)
           !
         enddo
         !
      enddo
      !
      if(trim(intensity%swap)=="SAVE") then 
        !
        close(unitO,status='keep')
        close(unitC,status='keep')
        !
        return 
      endif
	  !
	  write(out,"(/a,i0,' out of ',i0,' suggested compression  = ',f12.5/)") &
	               'Maximal number of non-zero values after vector compression  = ',cdimenmax,dimenmax,cfactor
      !
    endif
    !
    if (Ntransit==0) then 
         write(out,"('dm_intensity: the transition filters are too tight: no entry')") 
         stop 'dm_intensity: the filters are too tight' 
    endif 
    !
    !omp do private(ilevelF,jF,energyF,igammaF,quantaF,indF,ndegF,dimenF,unitF,unitO,unitC,idegF,irec,cdimen,idimen) schedule(guided)
    do ilevelF = 1, nlevels
      !
      if (job%verbose>=5.and.mod(ilevelF,nlevels/min(500,ilevelF))==0 ) write(out,"('ilevel = ',i0)") ilevelF
      !
      jF = eigen(ilevelF)%jval 
      !
      !energy and and quanta of the final state
      !
      energyF = eigen(ilevelF)%energy
      igammaF = eigen(ilevelF)%igamma        
      quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
      normalF(0:nmodes) = eigen(ilevelF)%normal(0:nmodes) 
      !
      call energy_filter_upper(jF,energyF,quantaF,normalF(0),passed)
      !
      call energy_filter_lower(jF,energyF,quantaF,normalF(0),passed_)
      !
      if (.not.passed.and..not.passed_) cycle
      !
      ilevelsG(igammaF) = ilevelsG(igammaF) + 1
      !
      indF = eigen(ilevelF)%jind
      !
      ndegF   = eigen(ilevelF)%ndeg
      dimenF = bset_contr(indF)%Maxcontracts
      !
      unitF = Jeigenvec_unit(indF) 
      unitO = Jfuncs_unit(indF)
      unitC = Jindex_unit(indF)
      !
      ! prescreen all eigenfunctions, compact and store on the disk
      !
      do idegF = 1, ndegF
        !
        !call TimerStart('Reading eigenfuncs-0')
        !
        irec = eigen(ilevelF)%irec(idegF)
        !read(unitF, rec = irec) vecF(idegF,1:dimenF)
        !
        read(unitC,rec=irec) cdimen
        cdimens(ilevelF,idegF) = cdimen
        !
        !
        !call TimerStop('Compacting eigenfuncs')
        !
        !call TimerStart('Store compacted eigenfuncs')
        !
        !call TimerStop('Store compacted eigenfuncs')
        !
      enddo
      !
      if (all(cdimens(ilevelF,1:ndegF)<=dimenmax_swap).and.iswap(igammaF)<swap_size(igammaF)) then
        !
        iswap(igammaF) = iswap(igammaF) + 1
        !
        istored(ilevelF) = iswap(igammaF)
        !
        !omp parallel do private(idegF,cdimen) schedule(guided)
        do idegF = 1, ndegF
          !
          cdimen = cdimens(ilevelF,idegF)
          !
          cdimen_swap(igammaF)%kmat(iswap(igammaF),idegF) = cdimen
          !
          irec = eigen(ilevelF)%irec(idegF)
          !
          read(unitO,rec=irec) vec_swap(igammaF)%mat(iswap(igammaF),idegF,1:cdimen) ! vecI(idegF,1:cdimen)
          read(unitC,rec=irec) cdimen,icoeff_swap(igammaF)%imat(iswap(igammaF),idegF,1:cdimen) ! icoeffF(idegF,1:cdimen)
          !
          !if (cdimen>dimenmax_swap) then 
          !   !
          !   write(out,"('intens: dimenmax_swap is too smal (',i,' vs ',i,') increase swap_factor minimun to ',f12.4)") dimenmax_swap,cdimen,real(cdimen,rk)/real(dimenmax,rk)
          !   stop 'intens: increase swap_factor!'
          !   !
          !endif
          !
          !icoeff_swap(igammaF)%imat(iswap(igammaF),idegF,1:cdimen) = icoeffF(idegF,1:cdimen)
          !
          !call dcopy(cdimen,vecI(idegF,1:cdimen),1,vec_swap(igammaF)%mat(iswap(igammaF),idegF,1:cdimen),1)
          !
          !vec_swap(igammaF)%mat(iswap(igammaF),idegF,1:cdimen) = vecF(idegF,icoeffF(idegF,1:cdimen))
          !
        enddo
        !omp end parallel do
        !
      else
        !
        do idegF = 1, ndegF
          !
          isave(indF) = isave(indF) + 1
          !
          isaved(ilevelF,idegF) = eigen(ilevelF)%irec(idegF)
          !
          !irec = isave(indF)
          !
          !write(unitO,rec=irec) vecI(idegF,1:cdimen)
          !write(unitC,rec=irec) icoeffF(idegF,1:cdimen)
          !
        enddo
        !
      endif
      !
    enddo
    !omp end do
    !
    deallocate(vecI,vecF)
    deallocate(icoeffF)
    !omp end parallel 
    !
    call TimerStop('Prescreening eigenfuncs')
    !
    write(out,"(/'...done!')")
    !
    allocate(vecI(sym%Maxdegen,dimenmax),vecF(sym%Maxdegen,dimenmax), stat = info)
    !
    call ArrayStart('intensity-vectors',info,size(vecI),kind(vecI))
    call ArrayStart('intensity-vectors',info,size(vecF),kind(vecF))
    !
    allocate(icoeffI(sym%Maxdegen,dimenmax),icoeffF(sym%Maxdegen,dimenmax), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error:  icoeffI,icoeffF - out of memory'
    !
    ! loop over final states -> count states for each symmetry
    !
    nlevelsG = 0
    !
    do ilevelF = 1, nlevels
       !
       jF = eigen(ilevelF)%jval 
       !
       !energy and and quanta of the final state
       !
       energyF = eigen(ilevelF)%energy
       igammaF = eigen(ilevelF)%igamma        
       quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
       normalF(0:nmodes) = eigen(ilevelF)%normal(0:nmodes) 
       !
       call energy_filter_upper(jF,energyF,quantaF,normalF(0),passed)
       !
       if (.not.passed) cycle
       !
       nlevelsG(igammaF) = nlevelsG(igammaF) + 1
       !
    enddo
    !
    write(my_fmt,'(a,i0,a)') "(a,",Nrepresen,"i8)"
    !
    write(out,my_fmt) 'Number of states for each symm = ',nlevelsG(:)
    !
#if (dipole_debug >= 0)
       write(out,"(' Total number of lower states = ',i8)") nlevelI
       write(out,"(' Total number of transitions  = ',i8)") Ntransit
#endif
    !
    !
    ! In order to speed up the line strength evaluation, 
    ! S_if = | <i|mu|f> |^2  = | \sum_{nm} C_in C_fm <n|mu|m> |^2
    ! in three steps:
    ! 1. Evaluating the expansion of the i-state 
    !    s_{im} =  \sum_{n} C_in  <n|mu|m> 
    ! 2. performing the expansion of the final state f:
    !    s_{if} =  \sum_{m} C_fm  s_{im}. 
    ! 3. Building S_{if}
    !    S_{if} = s_{if}^2
    !
    !  The temporaly object s_{im} will be referted to as 
    !  a half-linestrength "half_linestr"
    !
    allocate(half_linestr(nJ,sym%Maxdegen,3,dimenmax),stat=info)
    !
    call ArrayStart('half_linestr',info,size(half_linestr),kind(half_linestr))
    !
    !  The matrix where some of the eigenvectors will be stored
    !
    if (job%verbose>=5) call MemoryReport
    !
    write(out,"(/a,a,a,a)") 'Linestrength S(f<-i) [Debye**2],',' Transition moments [Debye],'& 
                          &,'Einstein coefficient A(if) [1/s],','and Intensities [cm/mol]'
    !
    ! Prepare the table header
    !
    select case (trim(intensity%action))
      !
      case('ABSORPTION','EMISSION')

       write(my_fmt_a,'(a,i0,a,i0,a,i0,a,i0,a)') "(/t4a1,t6a8,t17a1,t19a5,t25a3,t35a1,t42a2,t50a2,t62a5,t85,",&
                      nclasses,"(4x),1x,",nmodes,"(4x),3x,a2,14x,",nclasses,"(4x),1x,",nmodes,&
                      "(4x),8x,a7,10x,a7,12x,a7,12x,a2,8x,a2,8x,a1)"
       !
       write(out,my_fmt_a) 'J','Gamma <-','J','Gamma','Typ','Ei','<-','Ef','nu_if','<-','S(f<-i)','A(if)','I(f<-i)','Ni','Nf','N'
       !
      case('TM')
       !
       write(my_fmt_a,'(a,i0,a,i0,a,i0,a,i0,a)') &
                    "(t4a1,t6a8,t17a1,t19a5,t25a3,t35a2,t42a2,t52a2,t65a5,t84,",nclasses,"(3x),1x,",nmodes,&
                    "(3x),9x,a2,8x,",nclasses,"(3x),1x,",nmodes,"(3x),17x,a8,8x,a1,12x,a1,18x,a1,18x,a1)"
       !
       write(out,my_fmt_a) 'J','Gamma <-','J','Gamma','Typ','Ei','<-','Ef','nu_if','<-','TM(f<-i)','N','x','y'

       !
    end select 
    !
    ! ---------------------------------
    ! the actual intensity calculations
    ! --------------------------------- 
    !
    itransit = 0
    !
    ! loop over initial states
    !
    Ilevels_loop: do ilevelI = 1,nlevels
      !
      indI = eigen(ilevelI)%jind
      !
      !dimension of the bases for the initial states
      !
      dimenI = bset_contr(indI)%Maxcontracts
      !
      !energy, quanta, and gedeneracy order of the initial state
      !
      jI = eigen(ilevelI)%jval
      energyI = eigen(ilevelI)%energy
      igammaI  = eigen(ilevelI)%igamma
      quantaI(0:nmodes) = eigen(ilevelI)%quanta(0:nmodes)
      normalI(0:nmodes) = eigen(ilevelI)%normal(0:nmodes)
      ndegI   = eigen(ilevelI)%ndeg
      !
      call energy_filter_lower(jI,energyI,quantaI,normalI(0),passed)
      !
      if (.not.passed) cycle
      !
      ! where the eigenvector is stored 
      !
      unitI = Jeigenvec_unit(indI) 
      unitO = Jfuncs_unit(indI)
      unitC = Jindex_unit(indI)
      !
      !read eigenvector of initial state
      !
      !
      if (istored(ilevelI)==0) then
        !
        do idegI = 1, ndegI
          !
          call TimerStart('Reading eigenfuncs')
          !
          irec = isaved(ilevelI,idegI)
          !
          cdimen = cdimens(ilevelI,idegI)
          !
          read(unitC,rec=irec) cdimen,icoeffI(idegI,1:cdimen)
          read(unitO,rec=irec) vecI(idegI,1:cdimen)
          !
          !read(unitC, rec = irec) icoeffI(idegI,1:cdimen)
          !write(unitO,rec=irec) vecI(idegF,1:cdimen)
          !
          cdimenI(idegI) = cdimen
          !
          call TimerStop('Reading eigenfuncs')
          !
        enddo
        !
      else
        !
        do idegI = 1, ndegI
          !
          iswap_ = istored(ilevelI)
          !
          cdimen = cdimen_swap(igammaI)%kmat(iswap_,idegI) 
          !
          cdimenI(idegI) = cdimen
          !
          icoeffI(idegI,1:cdimen) = icoeff_swap(igammaI)%imat(iswap_,idegI,1:cdimen)
          !
          call dcopy(cdimen,vec_swap(igammaI)%mat(iswap_,idegI,1:cdimen),1,vecI(idegI,1:cdimen),1)
          !
        enddo
        !
      endif
      !
      ! Compute the half-linestrength
      !
      half_linestr = 0
      !
      do indF = 1, nJ
        !
        jF = Jval(indF)
        !
        dimenF = bset_contr(indF)%Maxcontracts
        !
        ! Check if it is really necessary to start the calculations for the given levelI -> jF, 
        ! i.e. one can skip the rest if no transitions will start from the given ilevelI and 
        ! finish anywehere at J= jF. 
        !
        call TimerStart('Intens_Filter-2')
        !
        passed = .false.
        !
        !omp parallel do private(ilevelF,energyF,igammaF,quantaF,passed_) schedule(guided) reduction(+:passed)
        do ilevelF = 1, nlevels
          !
          if (eigen(ilevelF)%jval/=jF) cycle 
          !
          energyF = eigen(ilevelF)%energy
          igammaF = eigen(ilevelF)%igamma        
          quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes) 
          !
          call intens_filter(jI,jF,energyI,energyF,igammaI,igammaF,quantaI,quantaF,igamma_pair,passed)
          !
          if (passed) exit
          !
          !passed = passed_
          !
        enddo
        !omp end parallel do 
        !
        call TimerStop('Intens_Filter-2')
        !
        if (.not.passed) cycle
        !
        select case (trim(intensity%action))
          !
        case('ABSORPTION','EMISSION')
          !
          igammaF = igamma_pair(igammaI)
          !
          ! The degeneracy of the basis vector have to be the same for the intial and final states:
          !
          ndegF  = ndegI
          !
          if ((jI/=intensity%J(1).or.jF/=intensity%J(1)).and. abs(jI-jF)<=1.and.jI+jF>=1) then 
            !
            idegI = 1
            !
            ! swtich to a reduced C3v/D3h case by commenting the next loop
            !
            do idegI = 1, ndegI
              !
              !ideg_range = 1
              !
              ! In case of real calculations we might reduce the problem to only one degeneracy
              !
              !if (ndegI>1) ideg_range = mod(idegI,2)+1
              !
              !#if (dipole_debug >=3)
              !  !
              !  ideg_range(1) = 1; ideg_range(2) = ndegF 
              !  !
              !#endif
              !
              if (job%rotsym_do) then 
                !
                call do_1st_half_linestrength_rotsym(jI,jF,indI,indF,cdimenI(idegI),dimenF,&
                                              icoeffI(idegI,1:dimenI),vecI(idegI,1:cdimenI(idegI)),&
                                              half_linestr(indF,idegI,1,:))
                                              !
              else
                !
                call do_1st_half_linestrength_II(jI,jF,indI,indF,cdimenI(idegI),dimenF,&
                                              icoeffI(idegI,1:dimenI),vecI(idegI,1:cdimenI(idegI)),&
                                              jmax,threej,half_linestr(indF,idegI,1,:))
              endif 
              !
            enddo
            !
          endif 
          !
        case('TM')
            !
            !do igammaF = 1,Nrepresen
            !
            !ndegF  = sym%degen(igammaF)
            !
            do idegI = 1, ndegI
              !
              !do idegF = 1, sym%maxdegen
                !
                call do_1st_half_tm_II(indI,indF,cdimenI(idegI),dimenF,&
                                    icoeffI(idegI,1:dimenI),vecI(idegI,1:cdimenI(idegI)),&
                                    half_linestr(indF,idegI,:,:))
                !
              !enddo
            enddo
            !
            !enddo
            !
        end select
        !
      enddo 
      !
      !loop over final states
      !
      Flevels_loop: do ilevelF = 1,nlevels
         !
         indF = eigen(ilevelF)%jind
         !
         !dimension of the bases for the final state
         !
         dimenF = bset_contr(indF)%Maxcontracts
         !
         !energy, quanta, and gedeneracy order of the final state
         !
         jF = eigen(ilevelF)%jval
         energyF = eigen(ilevelF)%energy
         igammaF  = eigen(ilevelF)%igamma
         quantaF(0:nmodes) = eigen(ilevelF)%quanta(0:nmodes)
         normalF(0:nmodes) = eigen(ilevelF)%normal(0:nmodes)
         ndegF   = eigen(ilevelF)%ndeg
         !
         ! where the eigenvector is stored 
         !
         unitF = Jeigenvec_unit(indF)
         unitO = Jfuncs_unit(indF)
         unitC = Jindex_unit(indF)
         !
         call energy_filter_upper(jF,energyF,quantaF,normalF(0),passed)
         !
         if (.not.passed) cycle Flevels_loop
         !
         call TimerStart('Intens_Filter-3')
         !
         call intens_filter(jI,jF,energyI,energyF,igammaI,igammaF,quantaI,quantaF,igamma_pair,passed)
         !
         call TimerStop('Intens_Filter-3')
         !
         if (.not.passed) cycle Flevels_loop
         !
         ! Which PQR branch this transition belong to ->
         ! 
         branch = PQR_branch(jI,jF)
         !
         nu_if = energyF - energyI 
         if (trim(intensity%action)=='EMISSION') nu_if = -nu_if 
         !
         ! Count the processed transitions 
         !
         itransit = itransit + 1
         !
         !read eigenvector of the final state
         !
         if (istored(ilevelF)==0) then
           !
           call TimerStart('Reading eigenfuncs')
           !
           !omp parallel do private(idegF,cdimen,irec)  shared(icoeffF,vecF,cdimenF) schedule(guided)
           !
           idegF = 1
           if (ndegF>1) idegF = 2
           !
           ! swtich to a reduced C3v/D3h case by commenting the next loop
           !
           do idegF = 1, ndegF
             !
             !irec = eigen(ilevelF)%irec(idegF)
             !
             irec = isaved(ilevelF,idegF)
             !
             cdimen = cdimens(ilevelF,idegF)
             !
             read(unitC, rec = irec) cdimen,icoeffF(idegF,1:cdimen)
             read(unitO, rec = irec) vecF(idegF,1:cdimen)
             !
             cdimenF(idegF) = cdimen
             !
           enddo
           !omp end parallel do
           !
           call TimerStop('Reading eigenfuncs')
           !
         !else
         !  !
         !  idegF = 1
         !  if (ndegI>1) idegF = 2
         !  !
         !  iswap_ = istored(ilevelF)
         !  !
         !  cdimenF(idegF) = cdimen_swap(igammaF)%kmat(iswap_,idegF)
         !  !
         !  cdimen = cdimenF(idegF)
         !  !
         !  M = mod(cdimen,5)
         !  !
         !  if (M/=0) then 
         !    !$omp parallel do private(idimen) shared(icoeffF,vecF) schedule(guided)
         !    do idimen = 1,M
         !      !
         !      icoeffF(idegF,idimen) = icoeff_swap(igammaF)%imat(iswap_,idegF,idimen)
         !      vecF(idegF,idimen)    = vec_swap(igammaF)%mat(iswap_,idegF,idimen)
         !      !
         !    enddo
         !    !$omp end parallel do
         !  endif
         !  !
         !  MP1 = M + 1
         !  !
         !  !$omp parallel do private(idimen) shared(icoeffF,vecF) schedule(guided)
         !  do idimen = MP1,cdimen,5
         !    !
         !    icoeffF(idegF,idimen  ) = icoeff_swap(igammaF)%imat(iswap_,idegF,idimen  )
         !    icoeffF(idegF,idimen+1) = icoeff_swap(igammaF)%imat(iswap_,idegF,idimen+1)
         !    icoeffF(idegF,idimen+2) = icoeff_swap(igammaF)%imat(iswap_,idegF,idimen+2)
         !    icoeffF(idegF,idimen+3) = icoeff_swap(igammaF)%imat(iswap_,idegF,idimen+3)
         !    icoeffF(idegF,idimen+4) = icoeff_swap(igammaF)%imat(iswap_,idegF,idimen+4)
         !    !
         !    vecF(idegF,idimen  )    = vec_swap(igammaF)%mat(iswap_,idegF,idimen  )
         !    vecF(idegF,idimen+1)    = vec_swap(igammaF)%mat(iswap_,idegF,idimen+1)
         !    vecF(idegF,idimen+2)    = vec_swap(igammaF)%mat(iswap_,idegF,idimen+2)
         !    vecF(idegF,idimen+3)    = vec_swap(igammaF)%mat(iswap_,idegF,idimen+3)
         !    vecF(idegF,idimen+4)    = vec_swap(igammaF)%mat(iswap_,idegF,idimen+4)
         !    !
         !  enddo
         !  !$omp end parallel do
         !  !
         endif
         !
         select case (trim(intensity%action))
           !
           case default
             !
             stop 'only ABSORPTION and TM are properly coded up to now'
             !
           case('ABSORPTION')
             !
             linestr_deg = 0 
             !
             !loops over degenerate states
             !
             call TimerStart('Linestrength')
             !
             idegI = 1 ; idegF = 1
             if (ndegF>1) idegF = 2 ! mod(idegI,2)+1
             !
             ! swtich to a reduced C3v/D3h case by commenting the next loop
             !
             Idegs_loop: do idegI = 1, ndegI
                !!
                !ideg_range = 1
                !!
                !if (ndegI>1) ideg_range = mod(idegI,2)+1
                !!
                !#if (dipole_debug >=3)
                !  !
                !  ideg_range(1) = 1; ideg_range(2) = ndegF 
                !  !
                !#endif
                !
                ! swtich to a reduced C3v/D3h case by commenting the next loop
                !
                Fdegs_loop: do idegF = 1,ndegF
                 !
                 ! complete the line strength 
                 !
                 if (istored(ilevelF)==0) then
                   !
                   cdimen = cdimenF(idegF)
                   !
                   !linestr_deg(idegI,idegF) = ddot(cdimen,half_linestr( indF,idegI,1,icoeffF(idegF,1:cdimen ) ),1,vecF(idegF,1:cdimen),1)
                   !
                   M = mod(cdimen,5)
                   !
                   dtemp = 0
                   !
                   if (M/=0) then 
                     !$omp parallel do private(idimen) reduction(+:dtemp) schedule(guided)
                     do idimen = 1,M
                       dtemp = dtemp + half_linestr(indF,idegI,1,icoeffF(idegF,idimen))*vecF(idegF,idimen)
                     enddo
                     !$omp end parallel do
                   endif
                   !
                   MP1 = M + 1
                   !
                   !$omp parallel do private(idimen) reduction(+:dtemp) schedule(guided)
                   do idimen = MP1,cdimen,5
                     !
                     dtemp = dtemp + & 
                                     half_linestr(indF,idegI,1,icoeffF(idegF,idimen  ))*vecF(idegF,idimen  ) + &
                                     half_linestr(indF,idegI,1,icoeffF(idegF,idimen+1))*vecF(idegF,idimen+1) + &
                                     half_linestr(indF,idegI,1,icoeffF(idegF,idimen+2))*vecF(idegF,idimen+2) + &
                                     half_linestr(indF,idegI,1,icoeffF(idegF,idimen+3))*vecF(idegF,idimen+3) + &
                                     half_linestr(indF,idegI,1,icoeffF(idegF,idimen+4))*vecF(idegF,idimen+4)
                     !
                   enddo
                   !$omp end parallel do
                   !
                   linestr_deg(idegI,idegF) = dtemp

                   !
                 else
                   !
                   iswap_ = istored(ilevelF)
                   !
                   cdimen = cdimen_swap(igammaF)%kmat(iswap_,idegF)
                   !
                   M = mod(cdimen,5)
                   !
                   dtemp = 0
                   !
                   if (M/=0) then 
                     !$omp parallel do private(idimen) reduction(+:dtemp) schedule(guided)
                     do idimen = 1,M
                       dtemp = dtemp + half_linestr( indF,idegI,1,icoeff_swap(igammaF)%imat(iswap_,idegF,idimen))*&
                                       vec_swap(igammaF)%mat(iswap_,idegF,idimen)
                     enddo
                     !$omp end parallel do
                   endif
                   !
                   MP1 = M + 1
                   !
                   !$omp parallel do private(idimen) reduction(+:dtemp) schedule(guided)
                   do idimen = MP1,cdimen,5
                     !
                     dtemp = dtemp + & 
                     half_linestr( indF,idegI,1,icoeff_swap(igammaF)%imat(iswap_,idegF,idimen  ))*&
                                   vec_swap(igammaF)%mat(iswap_,idegF,idimen  ) + &
                     half_linestr( indF,idegI,1,icoeff_swap(igammaF)%imat(iswap_,idegF,idimen+1))*&
                                   vec_swap(igammaF)%mat(iswap_,idegF,idimen+1) + &
                     half_linestr( indF,idegI,1,icoeff_swap(igammaF)%imat(iswap_,idegF,idimen+2))*&
                                   vec_swap(igammaF)%mat(iswap_,idegF,idimen+2) + &
                     half_linestr( indF,idegI,1,icoeff_swap(igammaF)%imat(iswap_,idegF,idimen+3))*&
                                   vec_swap(igammaF)%mat(iswap_,idegF,idimen+3) + &
                     half_linestr( indF,idegI,1,icoeff_swap(igammaF)%imat(iswap_,idegF,idimen+4))*&
                                   vec_swap(igammaF)%mat(iswap_,idegF,idimen+4)
                     !
                   enddo
                   !$omp end parallel do
                   !
                   linestr_deg(idegI,idegF) = dtemp
                   !
                   !linestr_deg(idegI,idegF) = ddot(cdimen,half_linestr( indF,idegI,1,icoeff_swap(igammaF)%imat(iswap_,idegF,1:cdimen)),1,vecF(idegF,1:cdimen),1)
                   !
                 endif
                 !
               end do Fdegs_loop
             !
             end do Idegs_loop
             !
             call TimerStop('Linestrength')
             !
             !sum up all degenerate components
             !
             ! swtich to a reduced C3v/D3h case by commenting this line and undcommenting the next one
             !
             linestr = sum( linestr_deg(1:ndegI,1:ndegF)**2 )/real(ndegI,rk)
             !
             !linestr = sum( linestr_deg(1:ndegI,1:ndegF)**2 )
             !
             ! calculate the intensity 
             !
             A_einst = A_coef_s_1*real(2*jI+1,rk)*linestr*abs(nu_if)**3
             !
             linestr = linestr * intensity%gns(igammaI) * real( (2*jI + 1)*(2 * jF + 1),rk )
             !
             boltz_fc = abs(nu_if) * exp(-(energyI-intensity%ZPE) * beta) * (1.0_rk - exp(-abs(nu_if) * beta))&
                        / intensity%part_func
             !
             ! intensity in cm/mol
             !
             absorption_int = linestr * intens_cm_molecule * boltz_fc
             !
             !
             if (absorption_int>=intensity%threshold%intensity.and.linestr>=intensity%threshold%linestrength) then 
               !
               nformat = ndegI*ndegF
               !
               iswap_ = istored(ilevelF)
               !
               cdimen = maxval(cdimens(ilevelF,1:ndegF),dim=1)

               write(my_fmt_a,'(a,i0,a,i0,a,i0,a,i0,a,a,i0,a,i0,a,i0,a)') &
                     "((i4,1x,a4,3x),a2,(i4,1x,a4,3x),a1,(2x,f11.4,1x),a2,(1x,f11.4,1x),f11.4,2x,a1,1x,a3,x,i3,1x,a1,1x,a1,1x,",&
                     nclasses,"(x,a3),1x,",nmodes,"(1x, i3),1x,a1,1x,a3,a1,1x,a3,x,i3,1x,a1,1x,a1,1x,",&
                     nclasses,"(x,a3),1x,",nmodes,"(1x, i3),1x,a1,1x,3(1x, es16.8),2x,(1x,i6,1x),",&
                     "a2,(1x,i6,1x),i8,1x,i12,1x,i8,x,a1,1x,",&
                     nmodes,"(1x, i3),1x,a1,1x,a3,1x,a1,1x,",nmodes,"(1x, i3),1x,a1,1x,",nformat,"(1x, es16.8))"

               !
               write(out,my_fmt_a)  &
                            jF,sym%label(igammaF),'<-',jI,sym%label(igammaI),branch, &
                            energyF-intensity%ZPE,'<-',energyI-intensity%ZPE,nu_if,                 &
                            '(',eigen(ilevelF)%cgamma(0),eigen(ilevelF)%krot,')',&
                            '(',eigen(ilevelF)%cgamma(1:nclasses),eigen(ilevelF)%quanta(1:nmodes),')','<- ', &
                            '(',eigen(ilevelI)%cgamma(0),eigen(ilevelI)%krot,')',&
                            eigen(ilevelI)%cgamma(1:nclasses),eigen(ilevelI)%quanta(1:nmodes), &
                            linestr,A_einst,absorption_int,&
                            eigen(ilevelF)%ilevel,'<-',eigen(ilevelI)%ilevel,&
                            itransit,cdimen,istored(ilevelF),'(',normalF(1:nmodes),')','<- ','(',normalI(1:nmodes),')',&
                            linestr_deg(1:ndegI,1:ndegF)
             endif
             !
           case('TM')
             !
             !loops over degenerate states
             !
             if (istored(ilevelF)==0) then
               !
               do idegI = 1, ndegI
                 do idegF = 1, ndegF
                  !
                  cdimen = cdimenF(idegF)
                  !
                  ! complete the tm
                  !
                  forall(imu=1:3) tm_deg(imu,idegI,idegF) = dot_product( half_linestr(indF,idegI,imu,icoeffF(idegF,1:cdimen)),&
                                                                         vecF(idegF,1:cdimen) )
                  !
                 end do
               end do
               !
             else
               !
               iswap_ = istored(ilevelF)
               !
               do idegI = 1, ndegI
                 do idegF = 1, ndegF
                  !
                  cdimen = cdimen_swap(igammaF)%kmat(iswap_,idegF)
                  !
                  ! complete the tm
                  !
                  forall(imu=1:3) tm_deg(imu,idegI,idegF) = dot_product( &
                                  half_linestr(indF,idegI,imu,icoeff_swap(igammaF)%imat(iswap_,idegF,1:cdimen)),&
                                  vec_swap(igammaF)%mat(iswap_,idegF,1:cdimen) )
                  !
                 end do
               end do
               !
             endif
             !
             !do idegI = 1, ndegI
             !  do idegF = 1, ndegF
             !   !
             !   ! complete the tm
             !   !
             !   forall(imu=1:3) tm_deg(imu,idegI,idegF) = dot_product( half_linestr(indF,idegI,imu,1:dimenF),vecF( idegF,1:dimenF ) )
             !   !
             !  end do
             !end do
             !
             !sum up all degenerate components
             !
             linestr = sqrt( sum( tm_deg(1:3,1:ndegI,1:ndegF)**2 ) )
             !
             boltz_fc = abs(nu_if) * exp(-(energyI-intensity%ZPE) * beta) * (1.0_rk - exp(-abs(nu_if) * beta))&
                        / intensity%part_func
             !
             ! intensity in cm/mol
             !
             absorption_int = linestr**2 * intens_cm_molecule * boltz_fc
             !
             if (linestr>=intensity%threshold%intensity) then 
               !
               nformat = ndegI*ndegF
               !
               write(my_fmt_a,'(a,i0,a,i0,a,i0,a,i0,a,a,i0,a,i0,a,i0,a)') &
                       "(i4,1x,a3,3x,a2,i4,1x,a3,3x,a1,(2x,f13.6,1x),a2,(1x, f13.6,1x),f12.6,2x,a1,a3,a1,i3,a1,1x,a1,",&
                       nclasses,"a3,a1,",nmodes,"(1x, i3),a1,2x,a3,a1,a3,a1,i3,,1x,a1,",nclasses,"a3,a1,",nmodes,&
                       "(1x,i3),a1,2(1x,es15.8),i8,2x,a1,",nmodes,"(1x, i3),a1,2x,a3,1x,a1,",nmodes,&
                       "(1x, i3),a1,3(",nformat,"(1x,f16.8)))"
                       !"(1x, i3),a1,3(",nformat,"(1x,f16.8,1x,3i1)))"

               !write(out, "((i4, 1x, a3, 3x),a2, (i4, 1x, a3, 3x),a1,&
               !             &(2x, f13.6,1x),a2,(1x, f13.6,1x),f12.6,2x,&
               !             &a1,a3,a1,i3,a1,1x,a1,<nclasses>a3,a1,<nmodes>(1x, i3),a1,2x,a3,   &
               !             &a1,a3,a1,i3,,1x,a1,<nclasses>a3,a1,<nmodes>(1x, i3),a1,   &
               !             &2(1x,es15.8),i8,&
               !             &2x,a1,<nmodes>(1x, i3),a1,2x,a3,   &
               !             &1x,a1,<nmodes>(1x, i3),a1,   &
               !             &3(<nformat>( 1x,f16.8,1x,3i1) ) )") &


               !
               write(out,my_fmt_a) &
                            !
                            jF,sym%label(igammaF),'<-',jI,sym%label(igammaI),branch, &
                            energyF-intensity%ZPE,'<-',energyI-intensity%ZPE,nu_if,                              &
                            '(',eigen(ilevelF)%cgamma(0),';',eigen(ilevelF)%krot,')',&
                            '(',eigen(ilevelF)%cgamma(1:nclasses),';',eigen(ilevelF)%quanta(1:nmodes),')','<- ', &
                            '(',eigen(ilevelI)%cgamma(0),';',eigen(ilevelI)%krot,')',&
                            '(',eigen(ilevelI)%cgamma(1:nclasses),';',eigen(ilevelI)%quanta(1:nmodes),')', &
                            linestr,absorption_int,itransit,&
                            '(',normalF(1:nmodes),')','<- ','(',normalI(1:nmodes),')',&
                            !( ( ( (tm_deg(imu,idegI,idegF),idegF,idegI,imu),idegF=1,ndegF ),idegI=1,ndegI ),imu=1,3 )
                            ( ( ( (tm_deg(imu,idegI,idegF)),idegF=1,ndegF ),idegI=1,ndegI ),imu=1,3 )
             endif 
             !
         end select
         !
      end do Flevels_loop
      !
      if (job%verbose>=5) call TimerReport
      !
    end do Ilevels_loop

    deallocate(vecI, vecF)
    call ArrayStop('intensity-vectors')
    !
    deallocate(half_linestr)
    call ArrayStop('half_linestr')
    !
    do irep = 1,Nrepresen
      !
      if (associated(vec_swap(irep)%mat)) deallocate(vec_swap(irep)%mat)
      if (associated(icoeff_swap(irep)%imat)) deallocate(icoeff_swap(irep)%imat)
      if (associated(cdimen_swap(irep)%kmat)) deallocate(cdimen_swap(irep)%kmat)
      !
    enddo
    !
    call ArrayStop('vec_swap')
    call ArrayStop('icoeff_swap')
    !
    deallocate(istored)
    call ArrayStop('istored')
    !
    deallocate(isaved)
    call ArrayStop('isaved')
    !
    deallocate(idimenmax)
    !
    deallocate(icoeffI, icoeffF)
    !
    call TimerStop('Intensity calculations')
    !
  end subroutine dm_intensity



  !
  function PQR_branch(jI,jF) result (X)
    !
    implicit none 
    integer(ik),intent(in)  :: jI,jF
    character(len=1)        :: X
    !
    select case(trim(intensity%action))
        !
      case default 
        !
        X = 'D'
        !
      case ('ABSORPTION')
        !
        if (jI>jF) then
          X = 'P'
        elseif(jI/=jF) then
          X = 'R'
        else
          X = 'Q'
        endif
        !
      case ('EMISSION')
        !
        if (jI>jF) then
          X = 'R'
        elseif(jI/=jF) then
          X = 'P'
        else
          X = 'Q'
        endif
        !
    end select
    !
  end function PQR_branch


     subroutine find_igamma_pair(igamma_pair)
      !
      implicit none 
      integer(ik),intent(out) :: igamma_pair(:)
      integer(ik)     :: igammaI,igammaF,ngamma
      !
      if (trim(intensity%action)=='TM') then 
        !
        igamma_pair = 1
        return 
        ! 
      endif
      !
      do igammaI = 1,sym%Nrepresen
        !
        ! count number of hits
        !
        ngamma = 0
        igamma_pair(igammaI) = igammaI
        !
        do igammaF = 1,sym%Nrepresen
          !
          if (igammaI/=igammaF.and.intensity%isym_pairs(igammaI)==intensity%isym_pairs(igammaF)) then 
            !
            igamma_pair(igammaI) = igammaF
            !
            ngamma = ngamma + 1 
            !
            if (ngamma>1) then 
              !
              write(out,"('find_igamma_pair: Assumption that selection rules come in pairs is not fulfilled!')")
              stop 'find_igamma_pair: Assumption that all selection rules work in pairs is not fulfilled!'
              !
            endif   
            !
          endif
          !
        enddo
        !
        if ( intensity%gns(igammaI)/=intensity%gns(igamma_pair(igammaI)) ) then 
          !
          write(out,"('find_igamma_pair: selection rules do not agree with Gns')")
          stop 'find_igamma_pair: selection rules do not agree with Gns!'
          !
        endif   
        !
      enddo 
      !
     end subroutine find_igamma_pair


     subroutine energy_filter_lower(J,energy,quanta,normal_0,passed)
       !
       implicit none 
       integer(ik),intent(in) :: J
       real(rk),intent(in)    :: energy
       integer(ik),intent(in) :: quanta(0:molec%nmodes),normal_0
       logical,intent(out)    :: passed
         !
         ! passed = .true.
         !
         ! if (.not.intensity%do) return
         !
         passed = .false.
         !
         if (                                                             &
             ! nuclear stat.weight: 
             !
             J>=intensity%J(1).and.                                       &
             J<=intensity%J(2).and.                                       &
             !
             energy-intensity%ZPE>=intensity%erange_low(1).and.           &
             energy-intensity%ZPE<=intensity%erange_low(2).and.           &
             !
             all(quanta(1:molec%nmodes) >= intensity%v_low(1:molec%nmodes, 1)).and.   &
             all(quanta(1:molec%nmodes) <= intensity%v_low(1:molec%nmodes, 2)) ) then 
             !
             passed = .true.
             !
         endif
         !
         if (job%triatom_sing_resolve) then
            if (J==0.and.normal_0/=0) then
              passed = .false.
            endif
         endif 
         !
     end subroutine energy_filter_lower



     subroutine energy_filter_upper(J,energy,quanta,normal_0,passed)
       !
       implicit none 
       integer(ik),intent(in) :: J
       real(rk),intent(in)    :: energy
       integer(ik),intent(in) :: quanta(0:molec%nmodes),normal_0
       logical,intent(out)    :: passed
         !
         ! passed = .true.
         !
         ! if (.not.intensity%do) return
         !
         passed = .false.
         !
         if (                                                             &
             ! nuclear stat.weight: 
             !
             J>=intensity%J(1).and.                                       &
             J<=intensity%J(2).and.                                       &
             !
             energy-intensity%ZPE>=intensity%erange_upp(1).and.           &
             energy-intensity%ZPE<=intensity%erange_upp(2).and.           &
             !
             all(quanta(1:molec%nmodes) >= intensity%v_upp(1:molec%nmodes, 1)).and.   &
             all(quanta(1:molec%nmodes) <= intensity%v_upp(1:molec%nmodes, 2)) ) then 
             !
             passed = .true.
             !
         endif 
         !
         if (job%triatom_sing_resolve) then
            if (J==0.and.normal_0/=0) then
              passed = .false.
            endif
         endif 
         !
     end subroutine energy_filter_upper



     subroutine degeneracy_filter(gammaI,gammaF,idegI,idegF,passed)
       !
       implicit none 
       integer(ik),intent(in) :: gammaI,gammaF,idegI,idegF
       logical,intent(out)    :: passed
         !
         ! passed = .true.
         !
         ! if (.not.intensity%do) return
         !
         passed = .true.
         !
         if ( sym%degen(gammaI)==1.or.sym%degen(gammaF)==1 ) return
         if (.not.intensity%reduced) return
         !
         select case (trim(sym%group))
           !
         case ("C3V(M)","C3V","D3H(M)","D3H")
           !
           passed = .false.
           if (idegI==2.and.idegF==1) passed = .true.
           !
         case("TD(M)","TD")
           !
           passed = .false.
           !
           if ( (gammaI==4.and.gammaF==5.and.idegI==2.and.idegF==1).or.&
                (gammaI==5.and.gammaF==4.and.idegI==3.and.idegF==1) ) passed = .true.
           !
           if ( (gammaI==3.and.gammaF==3.and.idegI==2.and.idegF==1) ) passed = .true.
           !
         case default
           !
           passed = .true.
           !
         end select 
         !
     end subroutine degeneracy_filter





     subroutine intens_filter(jI,jF,energyI,energyF,igammaI,igammaF,quantaI,quantaF,igamma_pair,passed)
        implicit none 
        !
        integer(ik),intent(in) :: jI,jF,igammaI,igammaF,quantaI(0:molec%nmodes),quantaF(0:molec%nmodes)
        real(rk),intent(in)    :: energyI,energyF
        integer(ik),intent(in) :: igamma_pair(sym%Nrepresen)
        real(rk)               :: nu_if
        logical,intent(out)    :: passed
        integer(ik)            :: nmodes

          passed = .false.
          !
          nmodes = molec%nmodes
          !
          nu_if = energyF - energyI 
          if (trim(intensity%action)=='EMISSION') nu_if = -nu_if 
          !
          !
          if (                                                             &
              ! nuclear stat.weight: 
              !
              intensity%gns(igammaI)>small_.and.                           &
              !
              ! absorption/emission go only in one direction
              !
              nu_if>-small_.and.                                           &
              !
              ! spectroscopic window
              !
              nu_if>=intensity%freq_window(1).and.                         &
              nu_if<=intensity%freq_window(2).and.                         &
              !
              jI>=intensity%J(1).and.                                      &
              jI<=intensity%J(2).and.                                      &
              !
              jF>=intensity%J(1).and.                                      &
              jF<=intensity%J(2).and.                                      &
              !
              energyI-intensity%ZPE>=intensity%erange_low(1).and.          &
              energyI-intensity%ZPE<=intensity%erange_low(2).and.          &
              !
              energyF-intensity%ZPE>=intensity%erange_upp(1).and.          &
              energyF-intensity%ZPE<=intensity%erange_upp(2).and.          &
              !
              all(quantaI(1:nmodes) >= intensity%v_low(1:nmodes, 1)).and.  &
              all(quantaI(1:nmodes) <= intensity%v_low(1:nmodes, 2)).and.  &
              !
              all(quantaF(1:nmodes) >= intensity%v_upp(1:nmodes, 1)).and.  &
              all(quantaF(1:nmodes) <= intensity%v_upp(1:nmodes, 2))) then 
              !
              passed = .true.
              !
          endif 
          !
          !
          if (trim(intensity%action)=='ABSORPTION'.or.trim(intensity%action)=='EMISSION') then 
             !
             ! In order to avoid double counting of transitions
             ! we exclude jI=jF==intensity%J(1), i.e. Q branch for the highest J is never considered:
             !
             passed = passed.and.                                              &
             !
             !(jF/=intensity%J(2).or.jF/=jI).and.                               &
             !
             (jF/=intensity%J(1).or.jI/=intensity%J(1)).and.                    &
             !
             ! selection rules: 
             !
             intensity%isym_pairs(igammaI)==intensity%isym_pairs(igammaF).and.  &
             !
             igamma_pair(igammaI)==igammaF.and.                                 &
             !
             ! selection rules from the 3j-symbols
             !
             abs(jI-jF)<=1.and.jI+jF>=1
             !
          endif
          !
     end subroutine intens_filter



      !
      ! In order to speed up the line strength evaluation, 
      ! S_if = | <i|mu|f> |^2  = | \sum_{nm} C_in C_fm <n|mu|m> |^2
      ! in three steps:
      ! 1. Evaluating the expansion of the i-state 
      !    s_{im} =  \sum_{n} C_in  <n|mu|m> 
      ! 2. performing the expansion of the final state f:
      !    s_{if} =  \sum_{m} C_fm  s_{im}. 
      ! 3. Building S_{if}
      !    S_{if} = s_{if}^2
      ! 
      !  This routine performs the first half of <i|mu|f>. 
      !
      subroutine do_1st_half_linestrength(jI,jF,indI,indF,cdimenI,dimenF,icoeff,vector,jmax,threej,half_ls)
        !
        implicit none 
        integer(ik),intent(in)  :: jI,jF,indI,indF,cdimenI,dimenF,jmax
        integer(ik),intent(in)  :: icoeff(:)
        real(rk),intent(in)     :: vector(:),threej(0:jmax,0:jmax,-1:1,-1:1)
        real(rk),intent(out)    :: half_ls(:)
        integer(ik)             :: irootF, cirootI, icontrF, icontrI, & 
                                   kF, kI, tauF, tauI,sigmaF, sigmaI, ktau, irootI
        real(rk)                :: ls, f3j, sq2

          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_1st_half_linestr')
          !
          half_ls    = 0
          !
          sq2 = 1.0_ark/sqrt(2.0_rk)
          !
          !loop over final state basis components
          !
          !$omp parallel do private(irootF,icontrF,ktau,kF,tauF,cirootI,irootI,icontrI,tauI,sigmaI,sigmaF,kI,f3j,ls) shared(half_ls) schedule(guided)
          loop_F : do irootF = 1, dimenF
               !
               icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
               ktau = bset_contr(indF)%ktau(irootF)
               tauF  = mod(ktau,2_ik)
               kF = bset_contr(indF)%k(irootF)
               !
               sigmaF = mod(kF, 3)*tauF
               !
               !loop over initial state basis components
               !
               loop_I : do cirootI = 1, cdimenI
                  !
                  irootI = icoeff(cirootI)
                  kI = bset_contr(indI)%k(irootI)
                  !
                  if (abs(kF - kI)>1) cycle loop_I
                  !
                  icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
                  ktau = bset_contr(indI)%ktau(irootI)
                  tauI  = mod(ktau,2_ik)
                  !
                  sigmaI = mod(kI, 3)*tauI
                  !
                  f3j  =  threej(jI, kI, jF - jI, kF - kI)                 
                  ! 
                  ! 3j-symbol selection rule
                  !
                  if (abs(f3j)<intensity%threshold%coeff) cycle loop_I
                  !
                  !index of the corresponding vibrational contracted matrix element (cind)
                  !
                  !irow = max(icontrF, icontrI)
                  !icol = min(icontrF, icontrI)
                  !cind = irow * (irow - 1) / 2 + icol
                  !
                  !compute line strength
                  !
                  ls = 0 
                  !
                  if (kF == kI) then
                      !
                      ls  =  real(tauF-tauI,rk) * dipole_me(icontrI,icontrF,3)*f3j*vector(irootI)
                      !
                  elseif(tauF/=tauI) then 
                      !
                      ls =  real((kF-kI)*(tauF-tauI),rk)*dipole_me(icontrI,icontrF,1)*f3j*vector(irootI) 
                      !
                      if (kI*kF /= 0) ls = ls*sq2
                      !
                  elseif(tauF==tauI) then 
                      !
                      ls =  -dipole_me(icontrI,icontrF,2)*f3j*vector(irootI)
                      !
                      if (kI*kF /= 0) ls = ls*sq2
                      !
                  endif
                  !
                  !if (kI*kF /= 0.and.kF/=kI) ls = ls*sq2
                  !
                  ! The factor I**(tauF-tauI) is equivalent to I*(tauF-tauI)
                  !
                  half_ls(irootF) = half_ls(irootF) + (-1.0_rk)**(sigmaI+kI)*ls
                  !
               end do  loop_I
               !
               half_ls(irootF) = half_ls(irootF)*(-1.0_rk)**(sigmaF)
               !
            end do   loop_F
            !$omp end parallel do
            !
            call TimerStop('do_1st_half_linestr')
            !
      end subroutine do_1st_half_linestrength



      subroutine do_1st_half_linestrength_II_symmvec(jI,jF,indI,indF,cdimenI,icoeffI,vector,jmax,threej,half_ls)

        implicit none 
        integer(ik),intent(in)  :: jI,jF,indI,indF,jmax,cdimenI,icoeffI(:)
        real(rk),intent(in)     :: vector(:),threej(0:jmax,0:jmax,-1:1,-1:1)
        real(rk),intent(out)    :: half_ls(:)
        integer(ik)             :: irootF, icontrF, icontrI, & 
                                   kF, kI, tauF, tauI,sigmaF, sigmaI, ktau,&
                                   irootI,dimenI, dimenF,cirootI, startK,endK
        real(rk)                :: ls, f3j, sq2

          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_1st_half_linestr')
          !
          half_ls    = 0
          !
          sq2 = 1.0_ark/sqrt(2.0_rk)
          !
          dimenI = bset_contr(indI)%Maxcontracts
          dimenF = bset_contr(indF)%Maxcontracts
          !
          !loop over final state basis components
          !
          !$omp parallel do private(irootF,icontrF,ktau,tauF,kF,sigmaF,cirootI,irootI,kI,icontrI,tauI,sigmaI,f3j,ls) shared(half_ls) schedule(static)
          loop_F : do irootF = 1, dimenF
               !
               icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
               !irlevelF = bset_contr(indF)%ktau(irootF)
               !irdegF   = bset_contr(indF)%k(irootF)
               !
               ktau = bset_contr(indF)%ktau(irootF)
               tauF  = mod(ktau,2_ik)
               kF = bset_contr(indF)%k(irootF)
               startK = max(kF-1,0)
               endK   = max(kF+1,0)
               !
               sigmaF = mod(kF, 3)*tauF
               !
               !loop over initial state basis components
               !
               !loop_I : do irootI = 1, dimenI
               !
               loop_I : do cirootI = 1, cdimenI
                  !
                  irootI = icoeffI(cirootI)
                  !
                  kI = bset_contr(indI)%k(irootI)
                  !
                  if (abs(kF - kI)>1) cycle loop_I
                  !if ( kI< kF - 1) cycle loop_I
                  !if ( kI> kF + 1) then 
                  !  half_ls(irootF) = half_ls(irootF)*(-1.0_rk)**(sigmaF)
                  !  cycle loop_F
                  !endif
                  !
                  icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
                  !
                  !irlevelI = bset_contr(indI)%ktau(irootI)
                  !irdegI   = bset_contr(indI)%k(irootI)
                  !
                  ktau = bset_contr(indI)%ktau(irootI)
                  tauI  = mod(ktau,2_ik)
                  !
                  sigmaI = mod(kI, 3)*tauI
                  !
                  f3j  =  threej(jI, kI, jF - jI, kF - kI)                 
                  ! 
                  ! 3j-symbol selection rule
                  !
                  if (abs(f3j)<intensity%threshold%coeff) cycle loop_I
                  !
                  !compute line strength
                  !
                  ls = 0 
                  !
                  if (kF == kI) then
                      !
                      ls  =  real(tauF-tauI,rk) * dipole_me(icontrI,icontrF, 3)
                      !
                  elseif(tauF/=tauI) then 
                      !
                      ls =  real((kF-kI)*(tauF-tauI),rk)*dipole_me(icontrI,icontrF,1) 
                      !
                      if (kI*kF /= 0) ls = ls*sq2
                      !
                  elseif(tauF==tauI) then 
                      !
                      ls =  -dipole_me(icontrI,icontrF,2)
                      !
                      if (kI*kF /= 0) ls = ls*sq2
                      !
                  endif
                  !
                  !if (kI*kF /= 0.and.kF/=kI) ls = ls*sq2
                  !
                  ! The factor I**(tauF-tauI) is equivalent to I*(tauF-tauI)
                  !
                  half_ls(irootF) = half_ls(irootF) + (-1.0_rk)**(sigmaI+kI)*ls*f3j*vector(cirootI)
                  !
               end do  loop_I
               !
               half_ls(irootF) = half_ls(irootF)*(-1.0_rk)**(sigmaF)
               !
            end do   loop_F
            !$omp end parallel do
            !
            call TimerStop('do_1st_half_linestr')
            !
      end subroutine do_1st_half_linestrength_II_symmvec


      subroutine do_1st_half_linestrength_III_symmvec(jI,jF,indI,indF,jmax,dimenmax,cdimen_blk,icoeff_blk,&
                                                       vector_blk,itau_blk,threej,half_ls)

        implicit none 
        integer(ik),intent(in)  :: jI,jF,indI,indF,jmax,dimenmax,cdimen_blk(0:jmax),icoeff_blk(dimenmax,0:jmax),&
                                   itau_blk(dimenmax,0:jmax)
        real(rk),intent(in)     :: vector_blk(dimenmax,0:jmax),threej(0:jmax,0:jmax,-1:1,-1:1)
        real(rk),intent(out)    :: half_ls(:)
        integer(ik)             :: irootF, icontrF, icontrI, & 
                                   kF, kI, tauF, tauI,sigmaF, sigmaI, ktau,&
                                   irootI,dimenI, dimenF,cirootI, startK,endK
        real(rk)                :: ls, f3j, sq2

          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_1st_half_linestr')
          !
          half_ls    = 0
          !
          sq2 = 1.0_ark/sqrt(2.0_rk)
          !
          dimenI = bset_contr(indI)%Maxcontracts
          dimenF = bset_contr(indF)%Maxcontracts
          !
          !loop over final state basis components
          !
          !$omp parallel do private(irootF,icontrF,ktau,tauF,kF,startK,endK,sigmaF,kI,f3j,cirootI,irootI,icontrI,tauI,sigmaI,ls) &
          !$omp& shared(half_ls) schedule(static)
          loop_F : do irootF = 1, dimenF
               !
               icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
               !irlevelF = bset_contr(indF)%ktau(irootF)
               !irdegF   = bset_contr(indF)%k(irootF)
               !
               ktau = bset_contr(indF)%ktau(irootF)
               tauF  = mod(ktau,2_ik)
               kF = bset_contr(indF)%k(irootF)
               startK = max(kF-1,0)
               endK   = min(kF+1,jI)
               !
               sigmaF = mod(kF, 3)*tauF
               !
               !loop over initial state basis components
               !
               !loop_I : do irootI = 1, dimenI
               loop_k : do kI = startK,endK
                 !
                 f3j  =  threej(jI, kI, jF - jI, kF - kI)                 
                 ! 
                 ! 3j-symbol selection rule
                 !
                 if (abs(f3j)<intensity%threshold%coeff) cycle loop_K
                 !
                 loop_I : do cirootI = 1, cdimen_blk(kI)
                    !
                    irootI = icoeff_blk(cirootI,kI)
                    !
                    !kI = bset_contr(indI)%k(irootI)
                    !if (abs(kF - kI)>1) cycle loop_I
                    !if ( kI< kF - 1) cycle loop_I
                    !
                    if ( irootI==0 ) then 
                      write(out,"('do_1st_half_linestrength_III_symmvec error: illegal index irootI =0 ')")
                      stop 'do_1st_half_linestrength_III_symmvec error: irootI = 0 :( '
                    endif
                    !
                    icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
                    !
                    !irlevelI = bset_contr(indI)%ktau(irootI)
                    !irdegI   = bset_contr(indI)%k(irootI)
                    !
                    !ktau = bset_contr(indI)%ktau(irootI)
                    !tauI  = mod(ktau,2_ik)
                    !
                    tauI  = itau_blk(cirootI,kI)
                    !
                    sigmaI = mod(kI, 3)*tauI
                    !
                    !compute line strength
                    !
                    ls = 0 
                    !
                    if (kF == kI) then
                        !
                        ls  =  real(tauF-tauI,rk) * dipole_me(icontrI,icontrF, 3)
                        !
                    elseif(tauF/=tauI) then 
                        !
                        ls =  real((kF-kI)*(tauF-tauI),rk)*dipole_me(icontrI,icontrF,1) 
                        !
                        if (kI*kF /= 0) ls = ls*sq2
                        !
                    elseif(tauF==tauI) then 
                        !
                        ls =  -dipole_me(icontrI,icontrF,2)
                        !
                        if (kI*kF /= 0) ls = ls*sq2
                        !
                    endif
                    !
                    ! for magnetic moments
                    !ls  =   sum(dipole_me(icontrI,icontrF, :))
                    !
                    !if (kI*kF /= 0.and.kF/=kI) ls = ls*sq2
                    !
                    ! The factor I**(tauF-tauI) is equivalent to I*(tauF-tauI)
                    !
                    half_ls(irootF) = half_ls(irootF) + (-1.0_rk)**(sigmaI+kI)*ls*f3j*vector_blk(cirootI,kI)
                    !
                 end do  loop_I
               end do  loop_K
               !
               half_ls(irootF) = half_ls(irootF)*(-1.0_rk)**(sigmaF)
               !
            end do   loop_F
            !$omp end parallel do
            !
            call TimerStop('do_1st_half_linestr')
            !
      end subroutine do_1st_half_linestrength_III_symmvec


      !
      ! symmilar procedure of evaluating the (left) half-transformation of the line strength 
      ! in case the proper treatment of the rotational symmetrization is performed. 
      ! this is the only way to treat Td(M) spectra
      !
      subroutine do_1st_half_linestrength_rotsym_symmvec(jI,jF,indI,indF,cdimenI,icoeffI,vector,half_ls)
        !
        implicit none 
        integer(ik),intent(in)  :: jI,jF,indI,indF,cdimenI,icoeffI(:)
        real(rk),intent(in)     :: vector(:)
        real(rk),intent(out)    :: half_ls(:)
        integer(ik)             :: irootF, cirootI, icontrF, icontrI, & 
                                   irlevelI, irlevelF, irdegI, irdegF, irootI,dJ, dimenI, dimenF
        real(rk)                :: f_w(3),dip

          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_1st_half_linestr')
          !
          half_ls    = 0
          !
          dJ = jF-jI
          !
          dimenI = bset_contr(indI)%Maxcontracts
          dimenF = bset_contr(indF)%Maxcontracts
          !
          !loop over final state basis components
          !
          !$omp parallel do private(irootF,icontrF,irlevelF,irdegF,cirootI,irootI,icontrI,irlevelI,irdegI,f_w,dip) &
          !$omp& shared(half_ls) schedule(static)
          loop_F : do irootF = 1, dimenF
               !
               icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
               irlevelF = bset_contr(indF)%ktau(irootF)
               irdegF   = bset_contr(indF)%k(irootF)
               !
               !loop over initial state basis components
               !
               loop_I : do cirootI = 1, cdimenI
                  !
                  irootI = icoeffI(cirootI)
                  icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
                  irlevelI = bset_contr(indI)%ktau(irootI)
                  irdegI   = bset_contr(indI)%k(irootI)
                  !
                  f_w(:) = wigner(indI,dJ)%rot(:,irlevelI,irlevelF,irdegI,irdegF)
                  !
                  dip = sum(dipole_me(icontrI,icontrF,:)*f_w(:))
                  !
                  half_ls(irootF) = half_ls(irootF) + vector(cirootI)*dip
                  !
               end do  loop_I
               !
            end do   loop_F
            !$omp end parallel do
            !
            call TimerStop('do_1st_half_linestr')
            !
      end subroutine do_1st_half_linestrength_rotsym_symmvec




      subroutine do_1st_half_tm(jI,jF,indI,igammaI,idegI,ijterm,vector,half_tm)
        !
        implicit none 
        integer(ik),intent(in)  :: jI,jF,indI,igammaI,idegI
        integer(ik),intent(in)  :: ijterm(:,:)
        real(rk),intent(in)     :: vector(:)
        real(rk),intent(out)    :: half_tm(:,:)
        integer(ik)             :: icontrF, icontrI, dimen, irow, ib, iterm, ielem, isrootI
        real(rk)                :: vec
          !
          if (jI/=0.or.jF/=0) then
            write(out,"(' do_1st_half_tm: can be used only for J=0, not for ',2i4)") jI,jF
            stop 'do_1st_half_tm: can be used only for J=0' 
          endif
          !
          half_tm    = 0
          !
          dimen = bset_contr(1)%Maxcontracts
          !
          !loop over final state basis components
          !
          !$omp parallel do private(icontrF,icontrI,irow,ib,iterm,ielem,isrootI,vec) shared(half_tm) schedule(guided)
          loop_F : do icontrF = 1, dimen
               !
               !loop over initial state basis components
               !
               loop_I : do icontrI = 1, dimen
                  !
                  !icontrI = bset_contr(1)%iroot_correlat_j0(irootI)
                  !irlevelI = bset_contr(1)%ktau(irootI)
                  !irdegI   = bset_contr(1)%k(irootI)
                  !
                  irow = bset_contr(1)%icontr2icase(icontrF,1)
                  ib   = bset_contr(1)%icontr2icase(icontrF,2)
                  !
                  iterm = ijterm(irow,igammaI) 
                  !
                  do ielem = 1,bset_contr(indI)%irr(igammaI)%N(irow)
                     !
                     isrootI = iterm+ielem 
                     !
                     vec = vector(isrootI)*bset_contr(indI)%irr(igammaI)%repres(isrootI,idegI,ib)
                     !
                     half_tm(:,icontrF) = half_tm(:,icontrF) + dipole_me(icontrI,icontrF, :) *  vec
                     !
                  enddo
                  !
               end do  loop_I
               !
            end do   loop_F
            !$omp end parallel do
            !
      end subroutine do_1st_half_tm


      subroutine do_1st_half_tm_symmvec(indI,indF,vector,half_tm)
        
        implicit none 
        integer(ik),intent(in)  :: indI,indF
        real(rk),intent(in)     :: vector(:)
        real(rk),intent(out)    :: half_tm(:,:)
        integer(ik)             :: irootF, icontrF, icontrI, dimenI, dimenF, irootI, icase,ilambda
          !
          !if (jI/=0.or.jF/=0) then
          !  write(out,"(' do_1st_half_tm_II: can be used only for J=0, not for ',2i4)") jI,jF
          !  stop 'do_1st_half_tm_II: can be used only for J=0' 
          !endif
          !
          half_tm    = 0
          !
          !loop over final state basis components
          dimenF = bset_contr(indF)%Maxcontracts
          dimenI = bset_contr(indI)%Maxcontracts
          !
          !loop over final state basis components
          !
          !$omp parallel do private(irootF,icase,ilambda,icontrF,irootI,icontrI) shared(half_tm) schedule(static)
          loop_F : do irootF = 1, dimenF

               icase   = bset_contr(indF)%icontr2icase(irootF, 1)
               ilambda = bset_contr(indF)%icontr2icase(irootF, 2)

               icontrF = bset_contr(indF)%icontr_correlat_j0(icase, ilambda)
               !
               !loop over initial state basis components
               !
               loop_I : do irootI = 1, dimenI
                  !
                  icase   = bset_contr(indI)%icontr2icase(irootI, 1)
                  ilambda = bset_contr(indI)%icontr2icase(irootI, 2)
                  !
                  icontrI = bset_contr(indI)%icontr_correlat_j0(icase, ilambda)
                  !
                  !icontrI = bset_contr(1)%iroot_correlat_j0(irootI)
                  !irlevelI = bset_contr(1)%ktau(irootI)
                  !irdegI   = bset_contr(1)%k(irootI)
                  !
                  half_tm(irootF,:) = half_tm(irootF,:) + dipole_me(icontrI,icontrF, :) *  vector(irootI)
                  !
               end do  loop_I
               !
            end do   loop_F
            !$omp end parallel do


      end subroutine do_1st_half_tm_symmvec



      subroutine do_1st_half_linestrength_II(jI,jF,indI,indF,cdimenI,dimenF,icoeff,vector,jmax,threej,half_ls)

        integer(ik),intent(in)  :: jI,jF,indI,indF,cdimenI,dimenF,jmax
        integer(ik),intent(in)  :: icoeff(:)
        real(rk),intent(in)     :: vector(:),threej(0:jmax,0:jmax,-1:1,-1:1)
        real(rk),intent(out)    :: half_ls(:)
        integer(ik)             :: irootF, cirootI, icontrF, icontrI, & 
                                   kF, kI, tauF, tauI,sigmaF, sigmaI, ktau, irootI
        real(rk)                :: ls, f3j, sq2

          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_1st_half_linestr')
          !
          half_ls    = 0
          !
          sq2 = 1.0_ark/sqrt(2.0_rk)
          !
          !loop over final state basis components
          !
          !$omp parallel do private(irootF,icontrF,ktau,kF,tauF,cirootI,irootI,icontrI,tauI,sigmaI,sigmaF,kI,f3j,ls) shared(half_ls) schedule(guided)
          loop_F : do irootF = 1, dimenF
               !
               icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
               ktau = bset_contr(indF)%ktau(irootF)
               tauF  = mod(ktau,2_ik)
               kF = bset_contr(indF)%k(irootF)
               !
               sigmaF = mod(kF, 3)*tauF
               !
               !loop over initial state basis components
               !
               loop_I : do cirootI = 1, cdimenI
                  !
                  irootI = icoeff(cirootI)
                  kI = bset_contr(indI)%k(irootI)
                  !
                  if (abs(kF - kI)>1) cycle loop_I
                  !
                  icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
                  ktau = bset_contr(indI)%ktau(irootI)
                  tauI  = mod(ktau,2_ik)
                  !
                  sigmaI = mod(kI, 3)*tauI
                  !
                  f3j  =  threej(jI, kI, jF - jI, kF - kI)                 
                  ! 
                  ! 3j-symbol selection rule
                  !
                  if (abs(f3j)<intensity%threshold%coeff) cycle loop_I
                  !
                  !index of the corresponding vibrational contracted matrix element (cind)
                  !
                  !irow = max(icontrF, icontrI)
                  !icol = min(icontrF, icontrI)
                  !cind = irow * (irow - 1) / 2 + icol
                  !
                  !compute line strength
                  !
                  ls = 0 
                  !
                  if (kF == kI) then
                      !
                      ls  =  real(tauF-tauI,rk) * dipole_me(icontrI,icontrF, 3)*f3j*vector(cirootI)
                      !
                  elseif(tauF/=tauI) then 
                      !
                      ls =  real((kF-kI)*(tauF-tauI),rk)*dipole_me(icontrI,icontrF,1)*f3j*vector(cirootI) 
                      !
                      if (kI*kF /= 0) ls = ls*sq2
                      !
                  elseif(tauF==tauI) then 
                      !
                      ls =  -dipole_me(icontrI,icontrF,2)*f3j*vector(cirootI)
                      !
                      if (kI*kF /= 0) ls = ls*sq2
                      !
                  endif
                  !
                  !if (kI*kF /= 0.and.kF/=kI) ls = ls*sq2
                  !
                  ! The factor I**(tauF-tauI) is equivalent to I*(tauF-tauI)
                  !
                  half_ls(irootF) = half_ls(irootF) + (-1.0_rk)**(sigmaI+kI)*ls
                  !
               end do  loop_I
               !
               half_ls(irootF) = half_ls(irootF)*(-1.0_rk)**(sigmaF)
               !
            end do   loop_F
            !$omp end parallel do
            !
            call TimerStop('do_1st_half_linestr')
            !
      end subroutine do_1st_half_linestrength_II

      !
      ! symmilar procedure of evaluating the (left) half-transformation of the line strength 
      ! in case the proper treatment of the rotational symmetrization is performed. 
      ! this is the only way to treat Td(M) spectra
      !
      subroutine do_1st_half_linestrength_rotsym(jI,jF,indI,indF,cdimenI,dimenF,icoeff,vector,half_ls)

        integer(ik),intent(in)  :: jI,jF,indI,indF,cdimenI,dimenF
        integer(ik),intent(in)  :: icoeff(:)
        real(rk),intent(in)     :: vector(:)
        real(rk),intent(out)    :: half_ls(:)
        integer(ik)             :: irootF, cirootI, icontrF, icontrI, & 
                                   irlevelI, irlevelF, irdegI, irdegF, irootI,dJ
        real(rk)                :: ls, f_w(3)

          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_1st_half_linestr')
          !
          half_ls    = 0
          !
          dJ = jF-jI
          !
          !loop over final state basis components
          !
          !$omp  parallel do private(irootF,icontrF,irlevelF,irdegF,cirootI,irootI,icontrI,irlevelI,irdegI,f_w,ls) shared(half_ls)&
          !$omp& schedule(guided)
          loop_F : do irootF = 1, dimenF
               !
               icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
               irlevelF = bset_contr(indF)%ktau(irootF)
               irdegF   = bset_contr(indF)%k(irootF)
               !
               !loop over initial state basis components
               !
               loop_I : do cirootI = 1, cdimenI
                  !
                  irootI = icoeff(cirootI)
                  !kI = bset_contr(indI)%k(irootI)
                  !
                  !if (abs(kF - kI)>1) cycle loop_I
                  !
                  icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
                  irlevelI = bset_contr(indI)%ktau(irootI)
                  irdegI   = bset_contr(indI)%k(irootI)
                  ! 
                  !irow = max(icontrF, icontrI)
                  !icol = min(icontrF, icontrI)
                  !cind = irow * (irow - 1) / 2 + icol
                  !
                  !compute line strength
                  !
                  f_w(:) = wigner(indI,dJ)%rot(:,irlevelI,irlevelF,irdegI,irdegF)
                  !
                  ls = sum(dipole_me(icontrI,icontrF,:)*f_w(:))*vector(cirootI)
                  !
                  half_ls(irootF) = half_ls(irootF) + ls
                  !
               end do  loop_I
               !
            end do   loop_F
            !$omp end parallel do
            !
            call TimerStop('do_1st_half_linestr')
            !
      end subroutine do_1st_half_linestrength_rotsym

      subroutine do_1st_half_tm_II(indI,indF,cdimenI,dimenF,icoeff,vector,half_tm)

        integer(ik),intent(in)  :: indI,indF,cdimenI,dimenF
        integer(ik),intent(in)  :: icoeff(:)
        real(rk),intent(in)     :: vector(:)
        real(rk),intent(out)    :: half_tm(:,:)
        integer(ik)             :: irootF, cirootI, icase, ilambda, icontrF, icontrI
          !
          half_tm    = 0
          !
          !loop over final state basis components
          !
          do irootF = 1, dimenF

               icase   = bset_contr(indF)%icontr2icase(irootF, 1)
               ilambda = bset_contr(indF)%icontr2icase(irootF, 2)

               icontrF = bset_contr(indF)%icontr_correlat_j0(icase, ilambda)

               !loop over initial state basis components

               do cirootI = 1, cdimenI

                  icase   = bset_contr(indI)%icontr2icase(icoeff(cirootI), 1)
                  ilambda = bset_contr(indI)%icontr2icase(icoeff(cirootI), 2)

                  icontrI = bset_contr(indI)%icontr_correlat_j0(icase, ilambda)


                  !index of the corresponding vibrational contracted matrix element (cind)


                  !irow = max(icontrF, icontrI)
                  !icol = min(icontrF, icontrI)
                  !cind = irow * (irow - 1) / 2 + icol

                  !compute TM


                  half_tm(:,irootF) = half_tm(:,irootF) + dipole_me(icontrI,icontrF, :) *  vector(cirootI)

               end do

            end do


      end subroutine do_1st_half_tm_II








  ! This routine generates the rotational basis functions set in the symmetry adapted contracted representation. 
  ! Each basis function gets a symmetry label here. 
  !
  subroutine contraced_rotational_dipole(njval, jval, jmax, threej)
    !
    implicit none 
    integer(ik), intent(in) :: njval, jval(njval),jmax
    real(rk),intent(in)     :: threej(0:jmax,0:jmax,-1:1,-1:1)
    !
    integer(ik)         :: alloc
    integer(ik)         :: dimenI,nrootsI,idegI,ndeg,dimenf,nrootsF,idegF,irootI,irootF,dk
    integer(ik)         :: jI,jF,icountI,icountF,kI,kF,tauI,tauF,NcountI,NtotalI(sym%Nrepresen),NcountF,NtotalF(sym%Nrepresen)
    integer(ik)         :: nlevelsI,nlevelsF,isym,i,k0,tau0,sigmaI,sigmaF,n,indI,indF,dJ
    integer(ik),allocatable   :: count_indexI(:,:),rot_kI(:),rot_tauI(:),count_indexF(:,:),rot_kF(:),rot_tauF(:),NdegI(:),NdegF(:)
    real(rk),allocatable      :: eigenvectsI(:,:),eigenvectsF(:,:)
    real(rk),allocatable   :: mat(:,:,:),mat_(:,:)

    real(rk)           :: f3j,f_IF,f_I


    ! In order to find the correlation of the vibrational quantum numbers at J=0 and J/=0 
    ! the J=0 quanta destibution has to be defined before. We assume here that this is the case, 
    ! bset_contr(1) is always reserved for J=0
    !
    allocate(wigner(njval,-1:1),stat=alloc)
    !
    Ndeg = min(sym%Maxdegen,2)
    !
    do indI = 1, njval
      !
      jI = jval(indI)
      !
      dimenI = 2*jI+1 
      nrootsI = dimenI
      !
      allocate(count_indexI(nrootsI,nrootsI),eigenvectsI(nrootsI,nrootsI),rot_kI(dimenI),rot_tauI(dimenI),NdegI(dimenI),stat=alloc)
      call ArrayStart('rotational_dipoleI',alloc,size(count_indexI),kind(count_indexI))
      call ArrayStart('rotational_dipoleI',alloc,size(eigenvectsI),kind(eigenvectsI))
      !
      rot_kI(1)   = 0
      rot_tauI(1) = mod(jI,2)
      !
      irootI = 1
      do k0 = 1,jI
        do tau0 = 0,1
          !
          irootI = irootI + 1
          !
          rot_kI(irootI)  = k0
          rot_tauI(irootI)= tau0
          !
        enddo 
      enddo
      ! 
      call MLrotsymmetry_generate(jI,job%verbose,count_indexI,eigenvectsI,NcountI,NtotalI)
      icountI = 0
      do isym = 1,sym%Nrepresen
        !
        do i = 1,NtotalI(isym)
          !
          icountI = icountI + 1
          !
          NdegI(icountI) = sym%degen(isym)
          !
        enddo
      enddo
      nlevelsI  = icountI
      !
      do indF = 1, njval
         !
         jF = jval(indF)
         !
         dJ = jF-jI
         !
         if (abs(dJ)>1) cycle
         !
         dimenF = 2*jF+1 
         nrootsF = dimenF
         !
         allocate(count_indexF(nrootsF,nrootsF),eigenvectsF(nrootsF,nrootsF),rot_kF(dimenF),rot_tauF(dimenF),&
                  NdegF(dimenF),stat=alloc)
         call ArrayStart('rotational_dipoleF',alloc,size(count_indexF),kind(count_indexF))
         call ArrayStart('rotational_dipoleF',alloc,size(eigenvectsF),kind(eigenvectsF))
         ! 
         call MLrotsymmetry_generate(jF,job%verbose,count_indexF,eigenvectsF,NcountF,NtotalF)
         !
         rot_kF(1)   = 0
         rot_tauF(1) = mod(jF,2)
         !
         irootF = 1
         do k0 = 1,jF
           do tau0 = 0,1
             !
             irootF = irootF + 1
             !
             rot_kF(irootF)  = k0
             rot_tauF(irootF)= tau0
             !
           enddo 
         enddo
         !
         icountF = 0
         do isym = 1,sym%Nrepresen
           !
           do i = 1,NtotalF(isym)
             !
             icountF = icountF + 1
             NdegF(icountF) = sym%degen(isym)
             !
             !do ideg = 1,Ndeg
             !  !
             !  iroot = count_index(icount,ideg)
             !  !
             !enddo
             !
           enddo
         enddo
         nlevelsF  = icountF
         !
         allocate(mat(3,dimenI,dimenF),mat_(dimenI,dimenF),stat=alloc)
         call ArrayStart('rotational_dipoleF',alloc,size(mat),kind(mat))
         call ArrayStart('rotational_dipoleF',alloc,size(mat_),kind(mat_))
         !
         allocate(wigner(indI,dJ)%rot(3,nlevelsI,nlevelsF,sym%maxdegen,sym%maxdegen),stat=alloc)
         call ArrayStart('wigner%rot',alloc,size(wigner(indI,dJ)%rot),kind(wigner(indI,dJ)%rot))
         !
         wigner(indI,dJ)%rot = 0
         !
         mat = 0
         !
         do irootI  = 1,dimenI
            !
            kI   = rot_kI(irootI)
            tauI = rot_tauI(irootI)
            sigmaI = mod(kI, 3)*tauI
            !
            f_I = (-1.0_rk)**(sigmaI+kI)
            !
            do irootF  = 1,dimenF
               !
               kF   = rot_kF(irootF)
               tauF = rot_tauF(irootF)
               sigmaF = mod(kf, 3)*tauF
               !
               f_IF = f_I*(-1.0_rk)**(sigmaF)
               !
               dk = kI - kF
               !
               if (abs(dk)>1) cycle
               !
               f3j  =  threej(jI, kI, jF - jI, kF - kI)
               !
               if (abs(f3j)<small_) cycle
               !
               if (kF == kI) then
                   !
                   mat(3,irootI,irootF) = real(tauF-tauI,rk)*f3j*f_IF
                   !
                   !mat(3,irootI,irootF) = f3j*f_IF
                   !
               elseif(tauF/=tauI) then 
                   !
                   mat(1,irootI,irootF) = real((kF-kI)*(tauF-tauI),rk)*f3j*f_IF
                   !
                   !mat(1,irootI,irootF) = f3j*f_IF
                   !
                   if (kI*kF /= 0) mat(1,irootI,irootF) = mat(1,irootI,irootF)/sqrt(2.0_rk)
                   !
               elseif(tauF==tauI) then 
                   !
                   mat(2,irootI,irootF) = -f3j*f_IF
                   !
                   !mat(2,irootI,irootF) = f3j
                   !
                   if (kI*kF /= 0) mat(2,irootI,irootF) = mat(2,irootI,irootF)/sqrt(2.0_rk)
                   !
               endif
               !
            enddo
            !
         enddo
         !
         do n = 1,3
           !
           mat_(:,:) = matmul( matmul( transpose(eigenvectsI),mat(n,:,:) ),eigenvectsF )
           !
           do icountI = 1,nlevelsI
             do idegI =1,NdegI(icountI)
               irootI = count_indexI(icountI,idegI)
               !
               do icountF = 1,nlevelsF
                 do idegF =1,NdegF(icountF)
                   irootF = count_indexF(icountF,idegF)
                   wigner(indI,dJ)%rot(n,icountI,icountF,idegI,idegF) = mat_(irootI,irootF)
                 enddo
               enddo
               !
             enddo
           enddo
           !
         enddo
         !
         deallocate(mat,mat_,count_indexF,eigenvectsF,rot_kF,rot_tauF,NdegF)
         call ArrayStop('rotational_dipoleF')
         !
      enddo
      !
      deallocate(count_indexI,eigenvectsI,rot_kI,rot_tauI,NdegI)
      call ArrayStop('rotational_dipoleI')
      !
    enddo
    !
  end subroutine contraced_rotational_dipole
  !
  !
  !
  subroutine do_reduced_density(indI,vector,imode,jmode,kmode,density)
  
    implicit none 
    integer(ik),intent(in)  :: indI,imode,jmode,kmode
    real(rk),intent(in)     :: vector(:)
    real(rk),intent(out)    :: density(:,:,:)
    integer(ik)             :: icontrI, & 
                               kI, tauI,ktau,ib,&
                               irootI,dimenI, irow,cnu_i(0:bset_contr(indI)%nclasses)
      !
      density    = 0
      !
      dimenI = bset_contr(indI)%Maxcontracts
      !
      !omp parallel do private(irootI,irow,ib,cnu_i,icontrI,ktau,tauI,kI) shared(density) schedule(static)
      do irootI = 1, dimenI
           !
           irow = bset_contr(indI)%icontr2icase(irootI,1)
           ib   = bset_contr(indI)%icontr2icase(irootI,2)
           !
           cnu_i(:) = bset_contr(indI)%contractive_space(:, irow)
           !
           icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
           !
           ktau = bset_contr(indI)%ktau(irootI)
           tauI  = mod(ktau,2_ik)
           kI = bset_contr(indI)%k(irootI)
           !
           density(cnu_i(imode),cnu_i(jmode),cnu_i(kmode)) = density(cnu_i(imode),cnu_i(jmode),cnu_i(kmode)) &
                                                        + vector(irootI)**2
        end do
        !omp end parallel do
        !
  end subroutine do_reduced_density



      subroutine do_Angular_momentum_average(jI,indI,vectorI,vectorF,jrot_mat)
      
        implicit none 
        integer(ik),intent(in)  :: jI,indI
        real(rk),intent(in)     :: vectorI(:),vectorF(:)
        real(rk),intent(out)    :: jrot_mat(:)
        integer(ik)             :: irootF, icontrF, icontrI, & 
                                   kF, kI, tauF, tauI, ktau,&
                                   irootI,dimenI, dimenF, irow, ib, indF,jrow,jb,Jk,dk,i1
          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_Angular_momentum_average')
          !
          jrot_mat    = 0
          indF = indI
          !
          dimenI = bset_contr(indI)%Maxcontracts
          dimenF = bset_contr(indF)%Maxcontracts
          !
          !loop over final state basis components
          !
          !omp parallel do private(irootF,icontrF,ktau,tauF,kF,irootI,kI,icontrI,tauI,Jk,dk,i1) shared(jrot_mat) schedule(static)
          loop_F : do irootF = 1, dimenF
               !
               icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
               !
               !irow = bset_contr(indI)%icontr2icase(irootF,1)
               !ib   = bset_contr(indI)%icontr2icase(irootF,2)
               !
               ktau = bset_contr(indF)%ktau(irootF)
               tauF  = mod(ktau,2_ik)
               kF = bset_contr(indF)%k(irootF)
               !
               !cnu_i(:) = bset_contr(indI)%contractive_space(:, irow)
               !deg_i(:) = bset_contr(indI)%index_deg(icase)%icoeffs(:,ib)
               !
               !loop over initial state basis components
               !
               loop_I : do irootI = 1, dimenI
                  !
                  kI = bset_contr(indI)%k(irootI)
                  !
                  if (abs(kF - kI)>1) cycle loop_I
                  !
                  icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
                  !
                  !jrow = bset_contr(indI)%icontr2icase(irootF,1)
                  !jb   = bset_contr(indI)%icontr2icase(irootF,2)
                  !
                  !cnu_f(:) = bset_contr(indI)%contractive_space(:, jrow)
                  !
                  !irlevelI = bset_contr(indI)%ktau(irootI)
                  !irdegI   = bset_contr(indI)%k(irootI)
                  !
                  ktau = bset_contr(indI)%ktau(irootI)
                  tauI  = mod(ktau,2_ik)
                  !
                  Jk = 1+ki+(ji*(ji+1) )/2
                  dk = ki - kf
                  !
                  if (abs(dk)<=1.and.icontrI==icontrF ) then 
                    !.and.all(deg_i(1:PT%Nclasses)==deg_j(1:PT%Nclasses))) then 
                    !
                    do i1=1,3
                      !
                      if ( (taui==tauf.and.i1==2).or.(taui/=tauf.and.i1/=2) ) then
                        !
                        !f_t = contr(0)%rot(i1)%coeff3d(Jk,dk,tau_i)
                        !
                        jrot_mat(i1) = jrot_mat(i1) &
                                                 + vectorI(irootI)*vectorF(irootF)*bset%rot%matelements(i1,Jk,dk,taui)
                                                 !
                                                 !ideg    = bset_contr(jind)%index_deg(icase)%icoeffs(0,ilambda)
                                                 !
                                                 !contr(0)%rot(i1)%coeff3d(Jk,dk,tau_i)
                        !
                      endif
                      !
                    enddo
                    !
                  endif 

                 !
               end do  loop_I
               !
            end do   loop_F
            !omp end parallel do
            !
            call TimerStop('do_Angular_momentum_average')
            !
      end subroutine do_Angular_momentum_average



      subroutine do_rotational_density(indI,vectorI,vectorF,density)
      
        implicit none 
        integer(ik),intent(in)  :: indI
        real(rk),intent(in)     :: vectorI(:),vectorF(:)
        real(rk),intent(out)    :: density(:,:)
        integer(ik)             :: irootF, icontrF, icontrI, & 
                                   kF, kI, tauF, tauI, ktau,&
                                   irootI,dimenI, dimenF, indF,termI,termF

          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_rotational_density')
          !
          density    = 0
          indF = indI
          !
          dimenI = bset_contr(indI)%Maxcontracts
          dimenF = bset_contr(indF)%Maxcontracts
          !
          !loop over final state basis components
          !
          !omp parallel do private(irootF,icontrF,ktau,tauF,kF,termF,irootI,kI,icontrI,tauI,termI) shared(density) schedule(static)
          loop_F : do irootF = 1, dimenF
               !
               icontrF = bset_contr(indF)%iroot_correlat_j0(irootF)
               !
               !irow = bset_contr(indI)%icontr2icase(irootF,1)
               !ib   = bset_contr(indI)%icontr2icase(irootF,2)
               !
               ktau = bset_contr(indF)%ktau(irootF)
               tauF  = mod(ktau,2_ik)
               kF = bset_contr(indF)%k(irootF)
               termF = 2*kF+tauF ; if (kF==0) termF = termF + 1
               !
               !cnu_i(:) = bset_contr(indI)%contractive_space(:, irow)
               !deg_i(:) = bset_contr(indI)%index_deg(icase)%icoeffs(:,ib)
               !
               !loop over initial state basis components
               !
               if (abs(vectorF(irootF))<small_) cycle
               !
               loop_I : do irootI = 1, dimenI
                  !
                  kI = bset_contr(indI)%k(irootI)
                  !
                  if (abs(kF - kI)>1) cycle loop_I
                  !
                  icontrI = bset_contr(indI)%iroot_correlat_j0(irootI)
                  !
                  !jrow = bset_contr(indI)%icontr2icase(irootF,1)
                  !jb   = bset_contr(indI)%icontr2icase(irootF,2)
                  !
                  !cnu_f(:) = bset_contr(indI)%contractive_space(:, jrow)
                  !
                  !irlevelI = bset_contr(indI)%ktau(irootI)
                  !irdegI   = bset_contr(indI)%k(irootI)
                  !
                  ktau = bset_contr(indI)%ktau(irootI)
                  tauI  = mod(ktau,2_ik)
                  termI = 2*kI+tauI ; if (kI==0) termI = termI + 1
                  !
                  if (icontrI==icontrF ) then 
                     !
                     density(termI,termF) = density(termI,termF) + vectorI(irootI)*vectorF(irootF)
                     !
                  endif 
                 !
               end do  loop_I
               !
            end do   loop_F
            !omp end parallel do
            !
            call TimerStop('do_rotational_density')
            !
      end subroutine do_rotational_density




      subroutine Transform_rotdens_to_coord_repres(jI,density)
      
        implicit none 
        integer(ik),intent(in)  :: jI
        complex(rk),intent(out) :: density(0:,0:,:,:)
        integer(ik)             :: k1,k2,i1,i2,iphi,itheta,n1,n2,sigma,nphi,ntheta,tau1,tau2,info,alloc_p,jsize,mI
        double complex          :: R(2*jI+1,2*jI+1)
        double complex,parameter:: alpha = (1.0d0,0.0d0),beta=(0.0d0,0.0d0)
        real(ark)                :: phi,theta,chi,dtheta,dphi,normfactor
        complex(rk),allocatable :: Dfunc(:,:,:)
        !complex(rk)             :: f_t,g_t
        double complex,allocatable :: T(:,:),RT(:,:),F(:,:)

          !
          !dms_tmp = dipole_me
          !
          call TimerStart('transform_rotdens_to_coord_repres')
          !
          jsize = 2*jI+1
          !
          ! We assume m to maximize the projection
          !
          mI = jI
          !
          ! R is a transformation to the Wang functions
          !
          R = 0
          !
          if (mod(jI,2)==0) then 
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
          !
          do k1 = 1,jI
             !
             sigma = mod(k1,3)
             !
             R(i1  ,i1  ) = cmplx(1.0_ark/sqrt(2.0_ark)           ,0.0_ark                                )
             R(i1  ,i1+1) = cmplx(0.0_ark                         ,-(-1.0_ark)**sigma/sqrt(2.0_ark)       )
             R(i1+1,i1  ) = cmplx((-1.0_ark)**(jI+k1)/sqrt(2.0_ark),0.0_ark                               )
             R(i1+1,i1+1) = cmplx(0.0_ark                         ,(-1.0_ark)**(jI+k1+sigma)/sqrt(2.0_ark))
             !
             i1 = i1 + 2
             !
          enddo
          !
          chi = 0
          !
          density    = 0
          !
          ntheta = analysis%res%ntheta
          nphi = analysis%res%nphi
          !
          dtheta = real(pi,rk)/(Ntheta)
          dphi   = 2.0_rk*real(pi,rk)/(Nphi)
          !
          allocate(Dfunc(-jI:jI,0:ntheta,0:nphi),stat = info)
          call ArrayStart('transform_rotdens_to_coord_repres:Dfunc',info,size(Dfunc),kind(Dfunc))
          !
          normfactor=sqrt(real(2*jI+1,rk)/(8.0_rk*pi**2))
          !
          !$omp parallel do private(itheta,theta,iphi,phi,i1,k1,tau1,n1) shared(Dfunc) schedule(static)
          do itheta = 0,Ntheta
            !
            if (job%verbose>=6) write(out,"('itheta = ',i8)") itheta
            !
            theta = itheta*dtheta
            !
            do iphi = 0,Nphi
               !
               if (job%verbose>=6) write(out,"(2i8)") itheta,iphi
               !
               phi = iphi*dphi
               !
               i1 = 0
               do k1 = 0,jI
                  do tau1 = 0,min(1,k1)
                    i1 = i1 + 1
                    !
                    n1=k1*(-1)**tau1
                    !
                    !g_t = ( Dfunc(n1,itheta,iphi) - Phi_rot(jI,n1,theta,phi) )
                    !
                    !f_t = abs( Dfunc(n1,itheta,iphi) - Phi_rot(jI,n1,theta,phi) )
                    !
                    !if ( any( (/abs(real(Dfunc(n1,itheta,iphi))),abs(imag(Dfunc(n1,itheta,iphi)))/)>sqrt(small_) ) ) then
                    !  continue
                    !endif 
                    !
                    !if ( any( (/real(f_t),imag(f_t)/)>sqrt(small_) ) ) then
                    !  continue
                    !endif 
                    !
                    !Dfunc(n1,itheta,iphi) =  Phi_rot(jI,n1,theta,phi)
                    !
                    Dfunc(n1,itheta,iphi) = calc_phirot(jI,mI,n1,theta,phi)
                    !
                    !Dfunc(n1,itheta,iphi) =  ddlmn_conj(jI,n1,0,phi,theta,chi,info)*normfactor
                    !
                    !Dfunc(n1,itheta,iphi) =  ddlmn_conj(jI,0,n1,chi,theta,phi,info)*normfactor
                    !
                    !if (info/=0) then 
                    !  !
                    !  Dfunc(n1,itheta,iphi) =  Phi_rot(jI,n1,theta,phi)
                    !  !
                    !endif
                    !
                    if (job%verbose>=6) write(out,"(i8,2x,2f14.6,2x,2e16.8)") n1,theta,phi,Dfunc(n1,itheta,iphi)
                    !
                  enddo
               enddo
               !
            enddo
            !
          enddo
          !$omp end parallel do
          !
          !$omp parallel private(T,F,RT,alloc_p) shared(density) 
          allocate (T(jsize,jsize),RT(jsize,jsize),F(jsize,jsize),stat=alloc_p)
          if (alloc_p/=0)  then 
          write(out,"('transform_rotdens_to_coord_repres: T - out of memory')") 
             stop 'transform_rotdens_to_coord_repres - T out of memory'
          endif 
          !
          !$omp do private(itheta,theta,iphi,phi,i1,k1,tau1,n1,i2,k2,tau2,n2) schedule(static)
          do itheta = 0,Ntheta
            !
            if (job%verbose>=6) write(out,"('itheta = ',i8,'!')") itheta
            !
            theta = itheta*dtheta
            !
            do iphi = 0,Nphi
               !
               if (job%verbose>=6) write(out,"(2i8)") itheta,iphi
               !
               phi = iphi*dphi
               !
               i1 = 0
               do k1 = 0,jI
                  do tau1 = 0,min(1,k1)
                    i1 = i1 + 1
                    !
                    n1=k1*(-1)**tau1
                    !
                    i2 = 0
                    do k2 = 0,jI
                       do tau2 = 0,min(1,k2)
                         i2 = i2 + 1
                         !
                         n2=k2*(-1)**tau2
                         !
                         ! prob. density of the two (i2,i1) primitive rotational funcions with m = 0 anc chi=0
                         !
                         T(i2,i1) =  conjg(Dfunc(n2,itheta,iphi))*Dfunc(n1,itheta,iphi)
                         !
                       enddo
                    enddo
                    !
                  enddo
               enddo
               !
               ! Rotational prob. density matrix in the symmetrized (Wang) representaion
               !
               if (job%verbose>=6) write(out,"('zgemm...')") 
               !
               !F(:,:) = matmul( transpose( conjg(R) ),matmul( T,R ) )
               !
               !F(:,:) = matmul( conjg( R ) , matmul( T,transpose(R) ) )
               !F(:,:) = matmul( transpose( R ) , matmul( T,conjg(R) ) )
               !
               !F(:,:) = matmul( ( R ) , matmul( T,conjg(transpose(R)) ) )
               !
               call zgemm('C','N',jsize,jsize,jsize,alpha,R ,jsize,T,jsize,beta,RT,jsize)
               call zgemm('N','N',jsize,jsize,jsize,alpha,RT,jsize,R,jsize,beta,F ,jsize)
               !
               if (job%verbose>=6) write(out,"('...done')") 
               !
               density(itheta,iphi,:,:) = F(:,:)*real(sin(theta),rk)
               !
            enddo
            !
          enddo
          !$omp end do
          !
          deallocate(T,RT,F)
          !$omp end parallel
          !
          deallocate(Dfunc)
          call ArrayStop('transform_rotdens_to_coord_repres:Dfunc')
          !
          call TimerStop('transform_rotdens_to_coord_repres')
          !
      end subroutine transform_rotdens_to_coord_repres


      subroutine convert_symvector_to_contrvector(indI,dimenI,igammaI,idegI,ijterm,symvec,contrvec)
       !
       integer(ik),intent(in) :: indI,dimenI,igammaI,idegI
       real(rk), intent(in)   :: symvec(:)
       integer(ik),intent(in) :: ijterm(:,:)
       real(rk), intent(out)  :: contrvec(:)
       integer(ik)   :: irootI,irow,ib,nelem,ielem,isrootI,iterm
       real(rk)      :: dtemp0
         !
         !$omp parallel do private(irootI,irow,ib,iterm,dtemp0,nelem,ielem,isrootI) shared(contrvec) schedule(static)
         do irootI = 1,dimenI
            !
            irow = bset_contr(indI)%icontr2icase(irootI,1)
            ib   = bset_contr(indI)%icontr2icase(irootI,2)
            !
            iterm = ijterm(irow,igammaI)
            !
            dtemp0 = 0
            !
            nelem = bset_contr(indI)%irr(igammaI)%N(irow)
            !
            do ielem = 1,nelem
               !
               isrootI = iterm+ielem 
               !
               dtemp0 = dtemp0 + symvec(isrootI)*bset_contr(indI)%irr(igammaI)%repres(isrootI,idegI,ib)
               !
            enddo
            !
            contrvec(irootI) = dtemp0
            !
         enddo
         !$omp end parallel do


     end subroutine convert_symvector_to_contrvector


end module dipole

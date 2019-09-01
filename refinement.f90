module refinement

 use accuracy,     only : ik, hik, rk, ark, cl, wl, out, small_
 use fields,       only : job,fitting,j0fit,FLNmodes,FLindexQ,FLQindex,FL_fdf,FLpoten4xi,&
                          FLfinitediffs_2d,FLpoten_linearized,analysis,action
 use timer,        only : IOstart,Arraystart,Arraystop,Arrayminus,Timerstart,Timerstop,MemoryReport,TimerReport
 use molecules,    only : MLcoord_direct,MLinvmat,MLinvmatark,MLcoordinate_transform_func,MLpotentialfunc
 use moltype,      only : manifold,molec, extF,ML_check_steps,MLdiag_ulen,MLlinur,MLfromlocal2cartesian
 use symmetry,     only : sym
 use lapack,       only : lapack_syevr,lapack_syev,lapack_syevd


 use tran,         only : TReigenvec_unit, bset_contr, read_contrind, & 
                          read_eigenval,index_correlation,eigen,Neigenlevels

 private


 real(rk) :: stadev_best=1e-04,stab_best=1e-12  ! best standard error and srability 
 character(len=cl) :: deriv_type = 'hellman'    ! hellman or direct - how we calculate derivatives. 
                                                ! direct means the finite diferentiation of energies wrt parameters.
 !character(len=6)   :: fit_type='dgelss'       ! to switch between fitting methods. 
 integer(ik) :: maxiter_as = 3                  ! maximal number of iterations to find a match for assignement


 public refinement_by_fit,external_expectation_values

 real(rk),allocatable       :: poten_me(:,:,:)
 integer(ik), allocatable   :: Jindex(:)
 !
 integer(ik), pointer :: Jeigenvec_unit(:,:)
 integer(ik),parameter :: fit_debug = 1

 type coeffT                      
    real(rk),pointer    :: coeff(:)
    integer(ik),pointer :: icoeff(:)
 end type coeffT
 
  type  FfitT
     integer(ik)  :: Nentries
     integer(ik)  :: IOunit
     integer(ik),pointer  :: ilevel(:)
     integer(ik),pointer  :: ideg(:)
  end type FfitT

  type(FfitT),allocatable,save :: fit(:,:)

 type FkmatT
      integer(ik),pointer   :: kmat(:,:)
 end type FkmatT



  type(FkmatT),pointer :: ijterm(:)
  

contains



 subroutine refinement_by_fit
   !
   implicit none
   !
   ! to perform the refinement of PES through the fitting to the experimental and the ab initio data. 

    integer(ik)          :: info

    integer(ik)              :: Jmin, Jmax, nJ, jind, j, igamma, idimen, Nterms, iterm, ierror, Nentries
    integer(ik), allocatable :: Jval(:)

       !
       if (job%verbose>=2) write(out,"(/'The refinement by the Simultaneous Fitting')")   
       !
       !
       ! Count J-s
       !
       nj = 0 
       !
       do j = 1, size(fitting%J_list)
          !
          if (fitting%J_list(j)/=-1) nj = nj + 1
          !
       end do
       !
       Jmax = maxval(fitting%J_list(:))
       Jmin = minval(fitting%J_list(:),mask=fitting%J_list.gt.-1)
       !
       if (Jmin>0) nJ = nJ + 1
       !
       allocate(Jval(nJ),Jindex(0:jmax),stat = info)
       if (info /= 0) stop 'refinement_by_fit allocation error: Jval - out of memory'
       !
       Jval = 0 
       jind = 1
       !
       do j = 1, size(fitting%J_list)
          !
          if (fitting%J_list(j)>0) then
            !
            jind = jind + 1
            Jval(jind) = fitting%J_list(j)
            Jindex(Jval(jind)) = jind
            !
          endif 
          !
       end do
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
       allocate(Jeigenvec_unit(nJ,sym%Nrepresen), stat = info)
       if (info /= 0) stop 'dm_tranint allocation error: Jeigenvec_unit - out of memory'
       !
       do jind=1,nJ
         do igamma = 1,sym%Nrepresen
           !
           Nentries = 0
           do j = 1, Neigenlevels
             if (eigen(j)%jval==Jval(jind).and.eigen(j)%igamma==igamma.and.job%isym_do(igamma)) Nentries = Nentries + 1
           enddo
           !
           if (Nentries==0) cycle
           !
           Jeigenvec_unit(jind,igamma) = TReigenvec_unit(jind,Jval,igamma)
         enddo
       enddo
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
         do igamma = 1,sym%Nrepresen
            !
            Nterms = 0 
            !
            do iterm = 1,bset_contr(jind)%Maxsymcoeffs
              !
              ijterm(jind)%kmat(iterm,igamma) = Nterms
              !
              Nterms = Nterms + bset_contr(jind)%irr(igamma)%N(iterm) 
              !
            enddo
            !
         enddo
       enddo
       !
       if (job%verbose>=2) write(out,"('...done!')")
       !
       if (action%band_fitting) then 
         !
         ! the fitting of the band centers to the experimental ro-vib. term values
         !
         ! precalculate all derivatives  of the Hamiltonian matrix wrt J=0 term values 
         !
         call prepare_diff_evib(nJ,Jval(1:nJ))
         !
         ! do the fitting
         !
         if (trim(job%IOfitpot_action)/='SPLIT'.and.trim(job%IOfitpot_action)/='JOIN') then 
           !
           call bandcentres_fitting(Jval)
           !
         endif
         !
       else
         !
         ! simultaneous fitting of the potential parameters to the experimental ro-vib. term values
         !
         ! precalculate all derivatives  of the Hamiltonian matrix wrt potential parameters
         !
         call prepare_pot_matrix(nJ,Jval(1:nJ))
         !
         ! do the fitting
         !
         if (trim(job%IOfitpot_action)=='READ') then 
           !
           call sf_fitting(Jval)
           !
         endif
         !
       endif
       !
       if (allocated(poten_me)) deallocate(poten_me)
       if (allocated(Jindex)) deallocate(Jindex)
       !
       if (job%verbose>=2) call MemoryReport
       !
       if (job%verbose>=2) call TimerReport
       !
 end subroutine refinement_by_fit
 !


 subroutine external_expectation_values

    ! to perform the refinement of PES through the fitting to the experimental and the ab initio data. 
    !
    integer(ik)              :: info
    !
    integer(ik)              :: Jmin, Jmax, nJ, jind, j, igamma, idimen, Nterms, iterm
    integer(ik), allocatable :: Jval(:)
       !
       if (job%verbose>=2) write(out,"(/'The refinement by the Simultaneous Fitting')")   
       !
       !
       ! Count J-s
       !
       nj = 0 
       !
       do j = 1, size(analysis%J_list)
          !
          if (analysis%J_list(j)/=-1.and.all(analysis%J_list(j)/=analysis%J_list(1:j-1))) nj = nj + 1
          !
       end do
       !
       Jmax = maxval(analysis%J_list(:))
       Jmin = minval(analysis%J_list(:),mask=analysis%J_list.gt.-1)
       !
       if (Jmin>0) nJ = nJ + 1
       !
       allocate(Jval(nJ),Jindex(0:jmax),stat = info)
       if (info /= 0) stop 'refinement_by_fit allocation error: Jval - out of memory'
       !
       Jval = 0 
       jind = 1
       !
       do j = 1, size(analysis%J_list)
          !
          if (analysis%J_list(j)>0.and.all(analysis%J_list(j)/=analysis%J_list(1:j-1))) then
            !
            jind = jind + 1
            Jval(jind) = analysis%J_list(j)
            Jindex(Jval(jind)) = jind
            !
          endif 
          !
       end do
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
       call read_eigenval(nJ, Jval(1:nJ))
       !
       allocate(Jeigenvec_unit(nJ,sym%Nrepresen), stat = info)
       if (info /= 0) stop 'dm_tranint allocation error: Jeigenvec_unit - out of memory'
       !
       do jind=1,nJ
         do igamma = 1,sym%Nrepresen
           Jeigenvec_unit(jind,igamma) = TReigenvec_unit(jind,Jval,igamma)
         enddo
       enddo
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
         do igamma = 1,sym%Nrepresen
            !
            Nterms = 0 
            !
            do iterm = 1,bset_contr(jind)%Maxsymcoeffs
              !
              ijterm(jind)%kmat(iterm,igamma) = Nterms
              !
              Nterms = Nterms + bset_contr(jind)%irr(igamma)%N(iterm) 
              !
            enddo
            !
         enddo
       enddo
       if (job%verbose>=2) write(out,"('...done!')")   
       !
       call calc_exp_values(nJ,Jval(1:nJ))
       !
       if (allocated(poten_me)) deallocate(poten_me)
       if (allocated(Jindex)) deallocate(Jindex)
       do jind=1,nJ
         deallocate(ijterm(jind)%kmat)
       enddo
       call ArrayStop('ijterm')
       !
       if (job%verbose>=2) call MemoryReport
       !
       if (job%verbose>=2) call TimerReport
       !
 end subroutine external_expectation_values


 !
 !
 !
 subroutine sf_fitting(Jval)

      integer(ik),intent(in)  :: Jval(:)

      integer(ik)    :: npts,en_npts,pot_npts,nused      ! number of data points: current, all obs. energies, all pot. points, and actuall used. 
      integer(ik)       :: parmax                  ! total number of potential parameters 
      integer(ik)       :: jrot                    ! current value of the rotational quantum number 
      integer(ik)       :: nmodes,ncoords,lwork,fititer,nclasses
      integer           :: alloc,info, ierror      ! error state variables 
      real(rk)          :: wtsum

      real(rk),allocatable :: wt_bit(:),wtall(:) ! weight factors - only for the energies and total.
      real(rk),allocatable :: deriv0(:),rjacob(:,:),eps(:),energy_(:),enercalc(:)
      real(rk),allocatable :: local(:,:),pot_values(:),wspace(:),Tsing(:),chi_(:,:)
      real(rk),allocatable :: al(:,:),bl(:),dx(:),ai(:,:),sterr(:),sigma(:)
      real(rk),allocatable :: am(:,:),bm(:)
      real(ark),allocatable :: al_ark(:,:),ai_ark(:,:)
      real(rk),allocatable :: mat(:,:),dmat(:,:),vector(:)
      character(len=cl),allocatable :: nampar(:)    ! parameter names 
      real(ark),allocatable :: pot_terms(:)
      !
      integer(ik),allocatable :: ivar(:),ilargest(:),ifitparam(:)
      !
      logical      :: still_run,do_deriv,deriv_recalc
      real(rk)     :: stadev_old,stadev_before,stability,stadev,tempx,deltax,v,potright,potleft,sum_sterr,conf_int
      real(rk)     :: ssq,rms,ssq1,ssq2,rms1,rms2,fit_factor,Smallest,Largest
      real(rk)     :: ezero_,a_wats = 1.0_rk, lambda = 0.01_rk,nu = 10.0_rk
      integer(ik)  :: iener1,iener2,i,numpar,itmax,j,jlistmax,rank,ncol,nroots
      integer(ik)  :: iJ,isym,jsym,iener,jener,irow,icolumn,l,ndigits,imu_t
      integer(ik)  :: potunit,enunit,abinitunit,i0,iter_th
      integer(ik)  :: nlevels,ilevel,Nentries,ientry,jentry,k,imode,iunit,irange(2),quanta(0:FLNmodes)
      real(rk)     :: vrange(2),f_t,v_l,v_r
      real(ark)    :: ar_t(1:molec%ncoords)
      real(rk)     :: xi(1:FLNmodes)
      real(ark)    :: chi(1:FLNmodes),dummy_xyz(molec%natoms,3)
      character(cl) :: filename, ioname, buf,dir,pot_suffix,symchar,jchar
      integer(ik)   :: ipowers(1:FLNmodes),nterms,igamma,iroot
      integer(hik)  :: matsize,matsize2
      real(rk),allocatable  ::  deriv(:,:)   ! to store the calculated derivatives 
      real(rk),allocatable  :: pot_matrix(:,:)
      real(rk),allocatable :: potparam(:)
      logical       :: do_symmetry
      character(len=1)   :: rng
      character(len=1),allocatable  :: mark(:)
      character(len=cl) :: my_fmt,my_fmt_pot1,my_fmt_pot2 !format for I/O specification
      character(len=cl) :: my_fmt_en2,my_fmt_par1,my_fmt_par2 !format for I/O specification
      character(len=wl) :: my_fmt_en1 !wider format for I/O specification
       !
       if (job%verbose>=2) write(out,"(/'The least-squares fitting ...')")
       !
       dummy_xyz = 0
       !
       call TimerStart('Simultaneous Fitting')
       !
       en_npts = fitting%Nenergies
       !
       filename = fitting%geom_file
       write(ioname, '(a, i4)') 'ab initio points for the fit '
       call IOstart(trim(ioname), potunit)
       open(unit = potunit, action = 'read',status='old',file = trim(filename))
       !
       ! file (temp) where all computed energies will be printed out 
       !
       filename =  trim(fitting%output_file)//'.en'
       write(ioname, '(a, i4)') 'calc. energies from the fit '
       call IOstart(trim(ioname), enunit)
       open(unit = enunit, action = 'write',status='replace' , file = filename)
       !
       ! file (temp) where all computed potential energy corrections are printed out 
       !
       filename =  trim(fitting%output_file)//'.pot'
       write(ioname, '(a, i4)') 'pot. energies at the ab initio geometries'
       call IOstart(trim(ioname), abinitunit)
       open(unit = abinitunit, action = 'write',status='replace' , file = filename)
       !
       pot_npts = 0
       buf(1:3)=""
       do while( buf(1:3)/="---" )
         !
         pot_npts=pot_npts+1
         read (potunit,"(a80)") buf
         !
       enddo 
       !
       pot_npts=pot_npts-1
       !
       write(out,"(/'Number of obs. energies: ',i9)") en_npts
       write(out,"(/'Number of potential energy data points: ',i9)") pot_npts
       !
       npts = en_npts + pot_npts
       !
       allocate (local(molec%ncoords,pot_npts),stat=info)
       call ArrayStart('coords-mat',info,size(local),kind(local))
       !
       parmax = extF%rank
       !
       if (parmax>molec%parmax) then 
         !
         write(out,"('fitting: potential and fitting parameters do not agree (parmax>molec%parmax):',2i0)") parmax,molec%parmax
         stop 'potential and fitting parameters do not agree'
         !
       endif
       !
       nmodes  = FLNmodes
       !
       ncoords = molec%ncoords
       !
       nclasses = size(eigen(1)%cgamma)-1
       !
       itmax = fitting%itermax
       !
       fit_factor = fitting%factor
       !
       jlistmax = size(Jval)

       ipowers = 0 
       ipowers(1) = extF%maxord(1)
       nterms = FLQindex(FLNmodes,ipowers)
       !
       allocate (enercalc(en_npts),stat=info)
       call ArrayStart('obs-en-mat',info,size(enercalc),kind(enercalc))
       allocate (mark(en_npts),stat=info)
       allocate (wtall(npts),stat=info)
       call ArrayStart('weight-mat',info,size(wtall),kind(wtall))
       allocate (wt_bit(npts),stat=info)
       call ArrayStart('weight-mat',info,size(wt_bit),kind(wt_bit))
       !
       allocate (potparam(parmax),nampar(parmax),ivar(parmax),ifitparam(parmax),&
                 al(parmax,parmax),ai(parmax,parmax),bl(parmax),dx(parmax),sterr(parmax),&
                 Tsing(parmax),pot_terms(extF%rank),al_ark(parmax,parmax),ai_ark(parmax,parmax),stat=info)
       call ArrayStart('potparam-mat',info,size(potparam),kind(potparam))
       call ArrayStart('potparam-mat',info,size(ivar),kind(ivar))
       call ArrayStart('potparam-dx',info,size(dx),kind(dx))
       call ArrayStart('potparam-sterr',info,size(sterr),kind(sterr))
       call ArrayStart('potparam-Tsing',info,size(Tsing),kind(Tsing))
       call ArrayStart('a-b-mat',info,size(al),kind(al))
       call ArrayStart('a-b-mat',info,size(bl),kind(bl))
       call ArrayStart('pot_terms',info,size(pot_terms),kind(pot_terms))
       !
       allocate (am(parmax,parmax),bm(parmax),stat=info)
       call ArrayStart('potparam-mat',info,size(am),kind(am))
       call ArrayStart('potparam-mat',info,size(bm),kind(bm))
       !
       allocate (pot_values(pot_npts),stat=info)
       call ArrayStart('pot_values-mat',info,size(pot_values),kind(pot_values))
       !
       ! Allocate objects, that will be used for the fitting procedure:
       !
       allocate (deriv0(parmax),rjacob(npts,parmax),eps(npts),stat=info)
       call ArrayStart('deriv0',info,size(deriv0),kind(deriv0))
       call ArrayStart('rjacob',info,size(rjacob),kind(rjacob))
       call ArrayStart('eps',info,size(eps),kind(eps))
       !
       if (fitting%robust>small_) then
         !
         allocate (sigma(npts),stat=info)
         call ArrayStart('eps',info,size(sigma),kind(sigma))
         !
         a_wats = fitting%watson
         !
       endif 
       !
       ! The last object to allocate - the lapack related work array
       !
       lwork = 50*parmax
       !
       allocate (wspace(lwork),stat=info)
       call ArrayStart('wspace',info,size(wspace),kind(wspace))
       !
       wtall = 0 
       !
       rewind(potunit)
       if (job%verbose>=2)  write (out,"(/'Read geometries and pot. energies...')") 
       !
       do i=1,pot_npts
         !
         if (fit_debug > 6) then
           write (out,"('i = ',i0)") i
         endif
         !
         read (potunit,*) ar_t(1:molec%ncoords),pot_values(i),wtall(en_npts+i)
         local(:,i) = ar_t(:)
         !
       enddo 
       !
       if (job%verbose>=2)  write (out,"('...done!')") 
       !
       close(potunit,status="keep")
       !
       ! Objects with energies and derivatives at the given fitting iteration 
       !
       potparam(:) = extF%coef(1,:)
       ivar(:)     = extF%ifit(1,:)
       nampar(:)   = extF%name(1,:)
       !
       forall(i=1:en_npts) wtall(i) = fitting%obs(i)%weight
       wtsum = sum(fitting%obs(1:en_npts)%weight) ; wtsum = max(wtsum,sqrt(small_))
       wtall(1:en_npts) = wtall(1:en_npts)/max(wtsum,sqrt(small_))
       !
       wtsum = sum(wtall(en_npts+1:npts))
       wtall(en_npts+1:npts) = wtall(en_npts+1:npts)/wtsum
       !
       ! 2. Factorizing the obs. weights by the factor "fit_factor":
       !
       wtall(1:en_npts) = wtall(1:en_npts)*fitting%factor
       !
       ! 3. And normilizing the all weight factors. 
       !
       wtsum = sum(wtall(1:npts))
       wtall(1:npts) = wtall(1:npts)/wtsum
       ! 
       ! Count how many data points actually will be fitted. 
       !
       nused=0
       wt_bit = 0 
       !
       do i=1,npts
         if (wtall(i)>small_) then 
           !
           nused=nused+1
           wt_bit(i) = 1.0_rk
           !
         endif
       enddo
       !
       ! sigma = exp. data precision for robust fit 
       !
       if (fitting%robust>small_) then
         !
         sigma = 1.0_rk
         !
         do i=1,npts
           if (wtall(i)>small_) sigma(i) = sigma(i)/sqrt(wtall(i))
           !if (wtall(i)>small_) sigma(i) = sigma(i)/sqrt(fitting%obs(i)%weight)
         enddo
         !
         wtsum = 1.0_rk ! sqrt(sum(sigma(1:en_npts)**2))
         !
         sigma(:) = sigma(:)*fitting%robust/wtsum
         !
       endif 
       !
       ! assigning mark with null
       !
       mark(:) = ' '
       !
       write(out,"('Number of data points used in the fit: ',i9)") nused
       !
       ! The fitting loop is about to start. 
       !
       ! fititer will count the iterations. 
       !
       fititer = 0
       !
       !
       ! define printing formats
       write(my_fmt,'(a,i0,a)') "(3i5,2x,a3,1x,3f13.4,2x,e9.2,2x,a1,i3,a1,1x,a1,",nmodes,"(i3),a1,a)"
       write(my_fmt_en1,'(a,i0,a,i0,a,i0,a)') "(3i5,2x,a3,1x,3f13.4,2x,e9.2,2x,a2,a3,a1,i3,a2,1x,a2,",&
                                              nclasses,"a3,a1,",nmodes,"(i3),a2,1x,a1,",nmodes+1,"(i3),a1,a)"
                                              !
       write(my_fmt_en2,'(a,i0,a,i0,a)') "(3i5,2x,a3,1x,3f13.4,2x,e9.2,2x,a2,a3,a1,i3,a2,1x,a2,",&
                                               nclasses,"a3,a1,",nmodes,"(i3),a2)"

       write(my_fmt_pot1,'(a,i0,a)') "(1h1,5x,a,a,a//4x,",ncoords,"(7x),a,7x,a,3x,a,3x,a/)"

       write(my_fmt_pot2,'(a1,i0,a)') "(",ncoords,"(2x,f18.9),3(x,g18.10),x,e12.4)"
       write(my_fmt_par1,'(a1,i0,a)') "(a8,4x,",Ncoords,"i3,1x,i2,e22.14)"
       !
       nlevels = Neigenlevels
       !
       ! The outer fitting loop - allows to restart the fitting in case 
       ! we decide to remove some of the varying parameters. 
       ! This option is working together with fit_type ='linur'
       !
       still_run = .true.
       outer_loop: do while (still_run)  
          !
          ! Initial values for the standard error and  stability.
          !
          stadev_old = 1e10
          stadev_before = 1e10
          stability = 1e10
          stadev    = 1e10
          !
          ! Count the actually varying number of parameters:
          !
          numpar  = 0
          ifitparam = 1
          do i=1,parmax
            if (ivar(i) > 0) then 
              numpar=numpar+1
              ifitparam(numpar) = i
            endif
          enddo 
          !
          if (numpar==0.and.itmax/=0) then 
            !
            write (out,"('No varying paramters, check input!')") 
            stop 'sf_fitting: No varying paramters'
            !
          endif 
          !
          rjacob = 0 
          !
          ! The loop starts here. 
          !
          do while( fititer.le.itmax .and. stadev.ge.fitting%target_rms.and. stability.ge.stab_best)
            !
            fititer = fititer + 1
            !
            ! If fit_factor is set zero - there would be no fit to the experimental energies. 
            !
            if (fit_factor<small_) rjacob(1:en_npts,:) = 0
            !
            write(out,"(/'Iteration = ',i8)") fititer
            write(enunit,"(/'Iteration = ',i8)") fititer
            !
            ! Reconstruct the potential expansion from the local to linearized coords.
            !
            pot_terms = potparam
            !
            eps = 0
            !
            ! Print out the calc. and obs.-calc., i.e. result of the fit. 
            ! Only fitted energies are printed. 
            !
            write(out,"(/1X,100('-'),/a,/1X,100('-'))") &
                       '|## |  N |  J | sym|     Obs.    |    Calc.   | Obs.-Calc. |   Weight |    K     vib. quanta'
            !
            write(enunit,"(/1X,100('-'),/a,/1X,100('-'))") &
                         '|## |  N |  J | Sym|     Obs.    |    Calc.   | Obs.-Calc. |   Weight |    K    quanta   (Calc./Obs.) '
            !
            do j = 1,jlistmax
               !
               jrot = jval(j)
               !
               deriv_recalc = .true.
               !
               loop_isym : do isym = 1,sym%Nrepresen
                 !
                 if (job%IOfitpot_divide) then
                   write(jchar, '(i4)') jrot
                   write(symchar, '(i4)') isym
                   pot_suffix = trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'_'
                 else
                   iunit = fit(isym,j)%IOunit
                 endif
                 !
                 Nentries = fit(isym,j)%Nentries
                 !
                 if (Nentries<1) cycle
                 !
                 ! Check if we don't need to process this pair of J,isym 
                 !
                 do_symmetry = .false.
                 do i = 1,en_npts
                   if (fitting%obs(i)%Jrot==Jrot.and.fitting%obs(i)%symmetry==isym) then 
                     do_symmetry = .true.
                     exit
                   endif
                 enddo
                 if (.not.do_symmetry) cycle
                 !
                 if (job%verbose>=6) then 
                    write (out,"(/'Matrix size = ',i8)") Nentries
                 endif
                 !
                 matsize = int(Nentries*(Nentries+1)/2,hik)
                 matsize2= int(Nentries*Nentries,hik)
                 !
                 allocate(pot_matrix(Nentries,Nentries),stat=alloc)
                 call ArrayStart('pot_matrix',alloc,1,kind(pot_matrix),matsize2)
                 !
                 if (.not.job%IOfitpot_divide) then
                   !
                   rewind(iunit) 
                   read(iunit) buf(1:12)
                   !
                 endif
                 !
                 call TimerStart('Energy calculations')
                 !
                 allocate (mat(Nentries,Nentries),stat=info)
                 call ArrayStart('fitting-mat',info,1,kind(mat),matsize2)
                 !
                 mat = 0
                 !
                 do ientry = 1, Nentries
                   !
                   ilevel = fit(isym,j)%ilevel(ientry)
                   mat(ientry,ientry) = eigen(ilevel)%energy
                   !
                 enddo 
                 !
                 do i=1,extF%rank
                   !
                   !read (iunit,rec=i) pot_matrix(1:matsize)
                   !
                   if (job%IOfitpot_divide) then
                     !
                     ! we can sckip the parameter if it is zero and not used in refinement for the sliced-representation
                     ! 
                     if (ivar(i)<=0.and.abs(pot_terms(i))<small_) cycle
                     !
                     call divided_slice_read(i,'potF',pot_suffix,Nentries,pot_matrix)
                     !
                   else
                     !
                     read(iunit) imu_t
                     !
                     if ( i/=imu_t ) then
                       write (out,"(' pot_matrix at J = ',i4,', isym = ',i4,' i <> j: ',2i8)") jrot,isym,i,imu_t
                       stop 'pot_matrix - i<>j'
                     end if
                     !
                     read (iunit) pot_matrix
                     !
                   endif
                   !
                   !$omp parallel do private(ientry,jentry) shared(mat) schedule(guided)
                   do ientry = 1, Nentries
                     !
                     do jentry = 1, ientry
                       !
                       !iroot = ientry*(ientry-1)/2+jentry
                       !
                       mat(ientry,jentry) = mat(ientry,jentry) + pot_terms(i)*pot_matrix(ientry,jentry)
                       !
                       mat(jentry,ientry) = mat(ientry,jentry)
                       !
                       if (fit_debug > 3) then
                         write (out,"('mat (',i0,',',i0,')= ',es14.7)") ientry,jentry,mat(ientry,jentry)
                       endif
                       !
                     enddo
                     !
                   enddo
                   !$omp end parallel do
                   !
                 enddo
                 !
                 ! find the smallest absolute diagonal value of mat
                 !
                 Smallest  = huge(1.0_rk)
                 k  = 1
                 !
                 do ientry=1,Nentries
                    if (abs(mat(ientry,ientry))<=Smallest) then 
                        Smallest = abs( mat(ientry,ientry) )
                        k = ientry
                    endif
                 enddo
                 !
                 if (fit_debug > 2) then
                    !
                    write (out,"(/'Smallest diag. value of mat = ',es14.7,'at k = ',i8,' upper_range = ',es14.7/)") &
                                  mat(k,k),k,job%upper_ener+mat(k,k)
                    !
                 endif
                 !
                 if (allocated(energy_)) deallocate(energy_)
                 !
                 allocate (energy_(Nentries),stat=info)
                 call ArrayStart('fit-ener-deriv',info,size(energy_),kind(energy_))
                 !
                 vrange(1) = -0.0_rk ; vrange(2) = job%upper_ener+mat(k,k)
                 irange(1) = 1 ; irange(2) = min(job%nroots(isym),Nentries)
                 !
                 if (job%nroots(isym)/=1e6) then
                    rng = 'I'
                 elseif (job%upper_ener/=1e9) then 
                    rng = 'V'
                 else
                    rng = 'A'
                 endif
                 !
                 !energy_(1) = 0
                 !nroots = 1
                 !
                 !call lapack_syevr(mat(1:Nentries,1:Nentries),energy_(1:Nentries),rng,iroots=nroots,vrange=vrange,irange=irange)
                 !
                 select case (trim(job%diagonalizer)) 

                 case default
                   !
                   write (out,"('sf_fitting: type of the diagonalizer  ',a,' unknown')") trim(job%diagonalizer)
                   stop 'bandcentres_fitting - wrong diagonalizer '
                   !
                 case('SYEV') 
                   !
                   call lapack_syev(mat(1:Nentries,1:Nentries),energy_(1:Nentries))
                   nroots = Nentries
                   !
                 case('SYEVD') 
                   !
                   call lapack_syevd(mat(1:Nentries,1:Nentries),energy_(1:Nentries))
                   nroots = Nentries
                   !
                 case('SYEVR') 
                   !
                   vrange(1) = -0.0_rk ; vrange(2) = job%upper_ener+mat(k,k)
                   irange(1) = 1 ; irange(2) = min(job%nroots(isym),Nentries)
                   !
                   if (job%nroots(isym)/=1e6) then
                      rng = 'I'
                   elseif (job%upper_ener/=1e9) then 
                      rng = 'V'
                   else
                      rng = 'A'
                   endif 
                   !
                   call lapack_syevr(mat(1:Nentries,1:Nentries),energy_(1:Nentries),rng,iroots=nroots,vrange=vrange,irange=irange)
                   !
                 end select
                 !
                 ! Zero point energy:
                 !
                 if (jrot==0.and.isym==1) ezero_ = energy_(1)
                 !
                 ! find the largest coefficient for each root for the assignment
                 !
                 allocate(ilargest(nroots),stat=info)
                 call ArrayStart('ilargest',info,size(ilargest),kind(ilargest))
                 !
                 do i = 1,nroots
                    !
                    Largest  = -1.0
                    k  = 1
                    !
                    do ientry=1,Nentries
                       if (abs(mat(ientry,i))>=Largest) then 
                           Largest = abs( mat(ientry,i) )
                           k = ientry
                       endif
                    enddo
                    !
                    ilargest(i) = k
                    !
                 enddo
                 !
                 call TimerStop('Energy calculations')
                 !
                 ! Check we if we need to computed the derivatives 
                 !
                 do_deriv = .false.
                 !
                 if (itmax.ge.1.and.fit_factor>1e-12.and.trim(deriv_type)=='hellman'.and.mod(fititer+2,3)==0) do_deriv = .true.
                 !
                 if (do_deriv) then 
                   !
                   if (.not.job%IOfitpot_divide) then
                     rewind(iunit) 
                     read(iunit) buf(1:12)
                   endif 
                   !
                   call TimerStart('Energy derivatives')
                   !
                   allocate (deriv(nroots,numpar),stat=info)
                   call ArrayStart('fit-ener-deriv',info,size(deriv),kind(deriv))
                   !
                   deriv = 0
                   !
                   allocate (dmat(Nentries,Nentries),stat=info)
                   call ArrayStart('dmat',info,size(dmat),kind(dmat))
                   !
                   ncol = 0
                   !
                   do i=1,extF%rank
                      !
                      if (.not.job%IOfitpot_divide) then
                        !
                        read(iunit) imu_t
                        !
                        read (iunit) pot_matrix
                        !
                      elseif(ivar(i)>0) then !.and.abs(pot_terms(i))>small_) then
                        !
                        call divided_slice_read(i,'potF',pot_suffix,Nentries,pot_matrix)
                        !
                      endif
                      !
                      if (ivar(i) <=0) cycle
                      !
                      ncol = ncol + 1 
                      !
                      !do ncol=1,numpar 
                      !i = ifitparam(ncol)
                      !read (iunit,rec=i) pot_matrix
                      !
                      !$omp parallel do private(ientry,jentry) shared(dmat) schedule(guided)
                      do ientry = 1, Nentries
                        !
                        do jentry = 1,ientry
                          !
                          !iroot = ientry*(ientry-1)/2+jentry
                          !
                          dmat(ientry,jentry) = pot_matrix(ientry,jentry)
                          dmat(jentry,ientry) = pot_matrix(ientry,jentry)
                          !
                        enddo
                        !
                      enddo
                      !$omp end parallel do
                      !
                      !$omp parallel private(vector,info) shared(deriv,deriv0)
                      allocate (vector(Nentries),stat=info)
                      if (info/=0)  stop 'fitting-vector - out of memory'
                      !
                      !$omp do private(iener,iJ,igamma,ientry) schedule(guided)
                      do iener = 1,en_npts
                        !
                        iJ = fitting%obs(iener)%Jrot ; igamma = fitting%obs(iener)%symmetry
                        !
                        if (iJ==Jrot.and.igamma==isym) then 
                          !
                          ientry  = fitting%obs(iener)%N ;  if (ientry>nroots.or.isym>sym%Nrepresen) cycle 
                          !
                          call dgemv('N',Nentries,Nentries,1.0d0,dmat,Nentries,mat(:,ientry),1,0.0d0,vector,1)
                          !
                          deriv(ientry,ncol)= dot_product(mat(:,ientry),vector(:))
                          !
                        endif
                        !
                      enddo
                      !$omp end do
                      !
                      deallocate(vector)
                      !$omp end parallel
                      !
                      if (jrot==0.and.isym==1) then 
                        !
                        if (fitting%obs(1)%N==1.and.1==fitting%obs(1)%Jrot.and.fitting%obs(1)%symmetry==0) then
                           !
                           deriv0(ncol) = deriv(1,ncol)
                           !
                        else
                           !
                           allocate (vector(Nentries),stat=info)
                           if (info/=0)  stop 'fitting-vector - out of memory'
                           !
                           call dgemv('N',Nentries,Nentries,1.0d0,dmat,Nentries,mat(:,1),1,0.0d0,vector,1)
                           deriv0(ncol)= dot_product(mat(:,1),vector(:))
                           !
                           deallocate(vector)
                           !
                        endif
                      endif
                      !
                   enddo 
                   !
                   deriv_recalc = .false.
                   !
                   call TimerStop('Energy derivatives')
                   !
                   deallocate(dmat)
                   call ArrayStop('dmat')
                   !
                 endif 
                 !
                 deallocate(pot_matrix)
                 call ArrayStop('pot_matrix')
                 !
                 deallocate(mat)
                 call ArrayStop('fitting-mat')
                 !
                 ! correct the addresses of the energy
                 !
                 if (abs(fitting%threshold_lock)>0) then 
                   !
                   if (fitting%threshold_lock<0)  maxiter_as = 1
                   !
                   i0 = 1
                   !
                   iJ = fitting%obs(1)%Jrot ; igamma = fitting%obs(1)%symmetry
                   !
                   do iener = 1,en_npts
                     !
                     mark(iener) = '*'
                     !
                     !if (iJ/=fitting%obs(iener)%jrot.or.igamma /= fitting%obs(iener)%symmetry) i0 = 1
                     !
                     iJ = fitting%obs(iener)%Jrot ; igamma = fitting%obs(iener)%symmetry
                     !
                     if (iJ==Jrot.and.igamma==isym) then
                       !
                       loop_thresh : do iter_th = 1,maxiter_as
                         !
                         loop_i : do i=i0,nroots
                            !
                            do jener = max(iener-10,1),iener-1
                              !
                              if (fitting%obs(jener)%jrot/=Jrot.or.fitting%obs(jener)%symmetry/=isym) cycle
                              !
                              if (fitting%obs(jener)%N == i) cycle loop_i
                              !
                            enddo
                            !
                            quanta(1:nmodes) = eigen(fit(isym,j)%ilevel(ilargest(i)))%quanta(1:nmodes)
                            quanta(0) = eigen(fit(isym,j)%ilevel(ilargest(i)))%krot
                            !
                            if (abs( fitting%obs(iener)%energy-(energy_(i)-ezero_) )<=real(iter_th,8)*abs(fitting%threshold_lock) &
                                .and.all(quanta(0:nmodes)==fitting%obs(iener)%quanta(0:nmodes))) then 
                                !
                                fitting%obs(iener)%N = i
                                mark(iener) = ' '
                                exit loop_thresh
                                !
                            endif 
                            !
                         enddo loop_i
                         !
                       enddo loop_thresh
                       !
                       if (mark(iener) == ' ')  i0 = max(fitting%obs(iener)%N-10,1)
                       !
                     endif 
                     !
                   enddo
                   !
                 endif
                 !
                 do iener = 1,en_npts
                   !
                   iJ = fitting%obs(iener)%Jrot ; igamma = fitting%obs(iener)%symmetry
                   !
                   if (iJ==Jrot.and.igamma==isym) then 
                     !
                     i  = fitting%obs(iener)%N ;  if (i>nroots.or.isym>sym%Nrepresen) cycle 
                     !
                     enercalc(iener) = energy_(i)-ezero_
                     !
                     eps(iener) = fitting%obs(iener)%energy-enercalc(iener)
                     !
                     if (do_deriv) rjacob(iener,1:numpar) = deriv(i,1:numpar)-deriv0(1:numpar)
                     !
                     ilevel = fit(isym,j)%ilevel(ilargest(i))
                     !
                     write (out,my_fmt) &
                            iener,i,iJ,sym%label(isym),enercalc(iener)+eps(iener),enercalc(iener),eps(iener),&
                            wtall(iener),&
                            '(',eigen(ilevel)%krot,')','(',eigen(ilevel)%quanta(1:nmodes),')',mark(iener)
                            !fitting%obs(i)%quanta(1:nmodes)
                     !
                   endif 
                   !
                 enddo ! --- i
                 !
                 ! Printing all calculated term values. If the obs. counterpats exist, 
                 ! the obs.-calc. are printed as well. 
                 ! This list can be used to identify the obs-calc pairs, i.e. the number of 
                 ! a term value as it appear in the list of calculated energies. 
                 !
                 ! write(out,"(/'Iteration = ',i4)") fititer
                 !
                 do i = 1,nroots
                    !
                    if (i==1.or.energy_(i)-ezero_>sqrt(small_)) then 
                       !
                       ilevel = fit(isym,j)%ilevel(ilargest(i))                       
                       !
                       jener = 1
                       !
                       do while (jener/=en_npts+1) 
                         !
                         if ( i==fitting%obs(jener)%N.and.Jrot==fitting%obs(jener)%Jrot.and.isym==fitting%obs(jener)%symmetry ) exit
                         !
                         jener = jener+1
                         !
                       enddo
                       !
                       if (jener<=en_npts) then 
                         !
                         write(enunit,my_fmt_en1) &
                            i,fitting%obs(jener)%N,Jrot,&
                            sym%label(isym),enercalc(jener)+eps(jener),&
                            enercalc(jener),eps(jener),wtall(jener),&
                            '( ',eigen(ilevel)%cgamma(0),';',eigen(ilevel)%krot,' )',&
                            '( ',eigen(ilevel)%cgamma(1:nclasses),';',eigen(ilevel)%quanta(1:nmodes),' )',&
                            '(',fitting%obs(jener)%quanta(0:nmodes),')',mark(jener)
                         !
                       else
                         !
                         !write(enunit,"(4i5,' ',3f13.4,2x,e8.2,5x,<nmodes>(i3))") i,0,Jrot,isym,0.0,energy_(i)-ezero_,0.0,0.0
                         !
                         write(enunit,my_fmt_en2) &
                            i,0,Jrot,sym%label(isym),0.0,energy_(i)-ezero_,0.0,0.0,&
                            '( ',eigen(ilevel)%cgamma(0),';',eigen(ilevel)%krot,' )',&
                            '( ',eigen(ilevel)%cgamma(1:nclasses),';',eigen(ilevel)%quanta(1:nmodes),' )'
                         !
                       endif 
                       !
                    endif 
                    !
                 enddo
                 !
                 ! switch off energies with large obs-calc assuming unwanted swapping 
                 if (fitting%threshold_obs_calc>small_) then
                     !
                     do iener = 1,en_npts
                       !
                       iJ = fitting%obs(iener)%Jrot ; igamma = fitting%obs(iener)%symmetry
                       !
                       if (iJ==Jrot.and.igamma==isym) then
                         !
                         if (abs(eps(iener))>fitting%threshold_obs_calc.and.wtall(iener)>small_) then
                           wtall(iener) = 0
                           if (fitting%robust>small_) sigma(iener) = 0
                           eps(iener) = 0
                           nused = nused - 1
                         endif
                         !
                       endif
                     enddo
                     nused = max(nused,1)
                     wtsum = sum(wtall(1:npts))
                     wtall(1:npts) = wtall(1:npts)/wtsum
                     !
                 endif
                 !
                 deallocate(energy_)
                 if (allocated(deriv)) deallocate(deriv )
                 call ArrayStop('fit-ener-deriv')
                 !
                 deallocate(ilargest)
                 call ArrayStop('ilargest')
                 !
               enddo loop_isym
               !
            enddo 
            !
            ! Alternative way of calculating the derivatives  - with the finite 
            ! differences. It is essentially slower and we use it only for 
            ! the testing of the xpect3 derivativies.
            !
            if (trim(deriv_type)/='hellman'.and.itmax.ge.1.and.fit_factor>1e-12) then
              !
              stop 'FINITE-DIFF->NOT ACTIVATED ' ! call finite_diff_of_energy(rjacob,potparam)
              !
            endif 
            !
            !
            ! Here the potential energy section starts. 
            !
            call TimerStart('Potential energy points')
            !
            do_deriv = .false.
            !
            if (itmax.ge.1.and.fit_factor<1e12.and.trim(deriv_type)=='hellman'.and.mod(fititer+2,3)==0) do_deriv = .true.
            !
            dummy_xyz = 0
            !
            do ientry=1,pot_npts
              !
              ar_t = local(:,ientry)
              !
              pot_terms(1:parmax) = molec%force(1:parmax) + potparam(1:parmax)
              !
              call MLfromlocal2cartesian(1,ar_t,dummy_xyz)
              !
              v =  MLpotentialfunc(ncoords,molec%natoms,ar_t,dummy_xyz,pot_terms)
              !
              !
              ! eps - epsilon = ab initio energies - calculated pot. energies,
              ! where we comntinue counting the fitting data points starting with en_npts - 
              ! obs. data. 
              !
              eps(en_npts+ientry) = pot_values(ientry)-v
              !
              !pot_terms = potparam
              !
              ! Calculate derivatives with respect to parameters 
              ! using the finite diff. method. 
              !
              if (do_deriv) then
                !
                do i = 1,parmax
                  !
                  if (ivar(i)<0) then 
                    pot_terms(i) = molec%force(i)
                  else
                    pot_terms(i) = 0
                  endif
                  !
                enddo
                !
                !omp do private(ncol,i) schedule(guided)
                do ncol=1,numpar 
                   !
                   i = ifitparam(ncol)
                   !
                   pot_terms(i) = 1.0_ark
                   rjacob(en_npts+ientry,ncol) =  MLpotentialfunc(ncoords,molec%natoms,ar_t,dummy_xyz,pot_terms)
                   pot_terms(i) = 0
                   !
                   !tempx=pot_terms(i)
                   !deltax=1e-3*abs(tempx)
                   !if (deltax<small_) deltax=1e-3
                   !
                   !pot_terms(i)=tempx+deltax
                   !v_r =  MLpotentialfunc(ar_t,pot_terms)
                   !
                   !pot_terms(i)=tempx-deltax
                   !v_l =  MLpotentialfunc(ar_t,pot_terms)
                   !
                   !pot_terms(i)=tempx
                   !rjacob(en_npts+ientry,ncol)=(v_r-v_l)/(2.0_rk*deltax)
                   !
                enddo ! --- ncol
                !omp end do
                !
              endif     
              !
            enddo  ! ---  ientry
            !
            call TimerStop('Potential energy points')
            !
            ! print out the potential energy corrections at the 'ab initio' geometries
            !
            rewind(abinitunit)
            !
            write(abinitunit,my_fmt_pot1) '************',' ab initio points','************','zero','calc.','zero-calc','weight'
            !
            do i=1,pot_npts
               !
               v = pot_values(i)-eps(i+en_npts)
               !
               write (abinitunit,my_fmt_pot2),& 
                      local(:,i), &
                      pot_values(i),v, &
                      eps(i+en_npts),wtall(i+en_npts) 
              !
            enddo

            !
            ! ssq  - weighted rms**2, rms  - root mean square deviation. 
            !
            ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
            rms=sqrt(sum(eps(1:npts)*eps(1:npts))/npts)
            !
            !
            ! Prepare the linear system a x = b as in the Newton fitting approach.  
            !
            if (itmax>=1) then
               !----- form the a and b matrix ------c
               ! form A matrix 
               do irow=1,numpar       
                 do icolumn=1,irow    
                   al(irow,icolumn)=sum(rjacob(1:npts,icolumn)*rjacob(1:npts,irow)*wtall(1:npts))
                   al(icolumn,irow)=al(irow,icolumn)
                   if (fit_debug > 2) then
                     write (out,"('al (',i0,',',i0,')= ',es14.7)") irow,icolumn,al(irow,icolumn)
                   endif
                 enddo
               enddo
               !
               ! form B matrix 
               do irow=1,numpar      
                 bl(irow)=sum(eps(1:npts)*rjacob(1:npts,irow)*wtall(1:npts))
                 if (fit_debug > 2) then
                   write (out,"('bl (',i0,')= ',es14.7)") irow,bl(irow)
                 endif
               enddo
               !
               !
               ! Using Marquardt's fitting method
               !
               ! solve the linear equatins for two values of lambda and lambda/10
               !
               ! Defining scalled (with covariance) A and b
               ! 
               ! form A matrix 
               do irow=1,numpar       
                 do icolumn=1,irow    
                   Am(irow,icolumn)=al(irow,icolumn)/sqrt( al(irow,irow)*al(icolumn,icolumn) )
                   Am(icolumn,irow)=Am(irow,icolumn)
                 enddo
                 bm(irow) = bl(irow)/sqrt(al(irow,irow))
               enddo
               !
               ! define shifted A as A =  A+lambda I
               ! lambda is Marquard's scaling factor
               !
               do irow=1,numpar       
                   Am(irow,irow)=Am(irow,irow)*(1.0_rk+lambda)
               enddo
               !
               ! Two types of the linear solver are availible: 
               ! 1. linur (integrated into the program, from Ulenikov Oleg)
               ! 2. dgelss - Lapack routine (recommended).
               !
               select case (trim(fitting%fit_type)) 
               !
               case default
                 !
                 write (out,"('fit_type ',a,' unknown')") trim(fitting%fit_type)
                 stop 'fit_type unknown'
                 !
               case('LINUR') 
                 !
                 call MLlinur(numpar,numpar,am(1:numpar,1:numpar),bm(1:numpar),dx(1:numpar),ierror)
                 !
                 ! In case of dependent parameters  "linur" reports an error = ierror, 
                 ! which is a number of the dependent parameter. We can remove this paramter 
                 ! from the fit and set its value to zero. And start the iteration again. 
                 !
                 if (ierror.ne.0) then 
                   do ncol=1,numpar 
                      !
                      i = ifitparam(ncol)
                      if  ( ncol.eq.ierror ) then 
                          ivar(i) = 0
                          potparam(i) = 0
                          !
                          write(out,"(i0,'-th is out - ',a8)") i,nampar(i)
                          !
                      endif 
                      !
                   enddo 
                   !
                   numpar  = 0
                   ifitparam = 1
                   do i=1,parmax
                     if (ivar(i) > 0) then 
                       numpar=numpar+1
                       ifitparam(numpar) = i
                     endif
                   enddo 
                   !
                   cycle outer_loop    
                  endif 
                  !
               case ('DGELSS')
                 !
                 ai = am 
                 call dgelss(numpar,numpar,1,ai(1:numpar,1:numpar),numpar,bm(1:numpar),numpar,Tsing,-1.D-12,rank,wspace,lwork,info)
                 !
                 if (info/=0) then
                   write(out,"('dgelss:error',i0)") info
                   stop 'dgelss'
                 endif
                 !
                 dx = bm
                 !
               end select 
               !
               ! convert back from Marquardt's representation
               !
               do ncol=1,numpar
                  dx(ncol) =  dx(ncol)/sqrt(al(ncol,ncol))
               enddo
               !
               !----- update the parameter values ------!
               !
               do ncol=1,numpar 
                  i = ifitparam(ncol)
                  potparam(i)=potparam(i)+dx(ncol)*fitting%fit_scale
               enddo
               !
               ! Robust fit: adjust the fitting weights
               !
               if (fitting%robust>small_) then
                 !
                 call robust_fit(numpar,fitting%robust,sigma(1:npts),eps(1:npts),wtall(1:npts))
                 !
                 ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
				 !
				 !fitting%robust = -1.0_rk
                 !
               endif 
               !
               ! Estimate standard deviation error. 
               !
               if ( nused.ne.numpar ) then 
                 stadev=sqrt(ssq/float(nused-numpar))
               else 
                 stadev=sqrt(ssq/nused)
               endif
               !               
               if (stadev<stadev_old) then
                 lambda = lambda/nu
               else 
                 lambda = min(lambda*nu,10000.0)
               endif
               !
               if (job%verbose>=5) write(out,"(/'Marquardts parameter lambda = ',g15.7)") lambda
               !
               ! Estimate the standard errors for each parameter using 
               ! the inverse matrix of a. 
               !
               al_ark = am
               call MLinvmatark(al_ark,ai_ark,numpar,info)
               ai = ai_ark
               !
               sum_sterr=0.0_rk
               do ncol=1,numpar 
                  i = ifitparam(ncol)
                  if (nused.eq.numpar) then  
                     sterr(ncol)=0
                  else
                     sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev/sqrt(al(ncol,ncol))
                     sum_sterr=sum_sterr+abs(sterr(ncol)/potparam(i))
                  endif
               enddo    
               !
               sum_sterr=sum_sterr/numpar 
               !
               ! This is how we define stability of the fit:
               ! as a relative change of stadev comparing with the step before. 
               !
               stability=abs( (stadev-stadev_old)/stadev )
               stadev_old=stadev
               stadev_before  = stadev_old
               !
            else
               !
               stadev=sqrt(ssq/nused)
               !
            endif
            !
            ! Print the updated parameters. 
            !
            write(out,"(/'Correction to potential parameters:')")
            !
            ! 'list'
            !
            if (all(extF%term(:,1,1)==-1)) then
               !
               do i=1,parmax
                  write (out,"(a8,4x,i2,e22.14)") nampar(i),ivar(i),potparam(i)
               enddo
               !
               write(out,"(/'Potential parameters:')")
               !
               do i=1,parmax
                 write (out,"(a8,4x,i2,e22.14)") nampar(i),ivar(i),molec%force(i)+potparam(i)
               enddo
               !
            else
               !
               ! 'powers'
               !
               do i=1,parmax
                  write (out,my_fmt_par1) nampar(i),(extF%term(l,1,i),l=1,Ncoords),ivar(i),potparam(i)
               enddo
               !
               write(out,"(/'Potential parameters:')")
               !
               do i=1,parmax
                 write (out,my_fmt_par1) nampar(i),(extF%term(l,1,i),l=1,Ncoords),ivar(i),molec%force(i)+potparam(i)
               enddo
               !
            endif 
            !
            !
            ! Output some statistics and results 
            !
            !  only if we are fitting:  
            !
            if (itmax.ne.0) then
              ! 
              write (out,"(//'Potential parameters rounded in accord. with their standard errors'/)")
              l = 0 
              do i=1,parmax
                if (ivar(i) > 0) then
                   !
                   l=l+1
                   ndigits = 0
                   conf_int = sterr(l)
                   do while (conf_int.le.10.0.and.ndigits<12)
                     ndigits = ndigits +1 
                     conf_int = conf_int*10
                   enddo
                   !
                   if (conf_int>1e10) conf_int = 0
                   !
                   write(my_fmt_par2,'(a,i0,a)') "(a8,i4,2x,f22.",ndigits,",a1,i14,a1)"
                   !
                   write (out,my_fmt_par2) nampar(i),ivar(i),potparam(i),'(',nint(conf_int),')'
                   !
                   !write (out,"(a8,i4,2x,f22.<ndigits>,'(',i14,')')") nampar(i),ivar(i),potparam(i),nint(conf_int)
                   !
                else 
                   !
                   ndigits =2
                   if (potparam(i).ne.0.0) ndigits = 8
                   !
                   write(my_fmt_par2,'(a,i0,a)') "(a8,i4,2x,f22.",ndigits,")"
                   !
                   write (out,my_fmt_par2) nampar(i),ivar(i),potparam(i)
                   !
                endif
              enddo  ! --- i
              !
            endif 
            !
            still_run = .false.
            !
            ! Print out the ssq for the rovib. energies and pot. data points separetely:
            !
            ssq1 = 0 ; ssq2 = 0 
            !
            wtsum = sum(wt_bit(1:en_npts))
            !
            if (wtsum/=0) ssq1 = sqrt( sum(eps(1:en_npts)**2*dble(wt_bit(1:en_npts)))/wtsum )
            !
            wtsum = sum(wt_bit(1+en_npts:npts))
            !
            if (wtsum/=0) ssq2 = sqrt( sum(eps(1+en_npts:npts)**2*dble(wt_bit(1+en_npts:npts)))/wtsum )
            !
            rms1=sqrt(sum(eps(1:en_npts)**2)/en_npts)
            rms2=sqrt(sum(eps(1+en_npts:npts)**2)/max(pot_npts,1))
            !
            write (out,6552) fititer,nused,numpar,stadev,ssq1,ssq2,stability
            !
            ! Print the potential energy points into a separate unit. 
            !
            !
            if (job%verbose>=6) call TimerReport
            !
          enddo  ! --- fititer
          !
       enddo outer_loop
       !
       if (allocated(potparam)) deallocate(potparam)
       !
       deallocate (nampar,ivar,ifitparam,al,ai,bl,dx,sterr,Tsing,pot_terms,al_ark,ai_ark)
       deallocate (am,bm)
       call ArrayStop('potparam-mat')
       call ArrayStop('potparam-dx')
       call ArrayStop('potparam-sterr')
       call ArrayStop('potparam-Tsing')
       call ArrayStop('a-b-mat')
       call ArrayStop('pot_terms')
       !
       deallocate (pot_values)
       call ArrayStop('pot_values-mat')
       !
       ! Allocate objects, that will be used for the fitting procedure:
       !
       deallocate (deriv0,rjacob,eps,wspace)
       call ArrayStop('deriv0')
       call ArrayStop('rjacob')
       call ArrayStop('eps')
       call ArrayStop('wspace')
       !
       if (allocated(sigma)) deallocate(sigma)
       !
       close(enunit,status="keep")
       close(abinitunit,status="keep")
       !
       call TimerStop('Simultaneous Fitting')
       !


6552   format(/3X,88('-')/'   |  Iter  | Points | Params |    Deviat     |',&
       '     ssq_ener  |    ssq_pot  | Convergence |'/&
       3X,88('-')/,&
       '   | ',i6,' | ',i6,' | ',i6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' | ',E10.3,'  |',/3X,88('-')/)

6553   format(/3X,88('-')/'   |  Iter  | Points | Params |    Deviat     |',&
       '     rms_ener  |    rms_pot  | Convergence |'/&
       3X,88('-')/,&
       '   | ',i6,' | ',i6,' | ',i6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' | ',E10.3,'  |',/3X,88('-')/)


  contains 




    subroutine finite_diff_of_energy(rjacob,potparam)
      !
      real(rk) :: rjacob(:,:),potparam(:)
      real(rk),allocatable :: enerleft(:),enerright(:),energy_(:)
      real(rk)  :: ezero_,deltax
      integer(ik)  :: i,ncol,k,gamma
      character(len=1)   :: rng
      !
      ! calculate derivatives with respect to parameters
      !
      allocate (enerright(en_npts),enerleft(en_npts),stat=alloc)
      call ArrayStart('enerright-enerleft',info,size(enerright),kind(enerright))
      call ArrayStart('enerright-enerleft',info,size(enerleft),kind(enerleft))
      !
      rjacob = 0 
      !
      ncol=0
      !
      do  i=1,parmax
        if (ivar(i) > 0) then
          !
          write(out,"('jacob:',i4,'-th param')") i
          !
          ncol=ncol+1
          tempx=potparam(i)
          deltax=1e-4*abs(tempx)
          if (deltax .le. 1e-15) deltax=1e-6
          !
          potparam(i)=tempx+deltax
          !
          do isym = 1,sym%Nrepresen
            !
            ! Determine the size of the Hamiltonian
            !
            Nentries = 0
            !
            do ilevel = 1, nlevels
              !
              if (eigen(ilevel)%jval==jrot.and.eigen(ilevel)%igamma==isym) Nentries = Nentries + 1
              !
            enddo
            !
            ! Construct the Hamiltonian matrix vector by vector
            !
            if (Nentries<1) cycle
            !
            allocate (mat(Nentries,Nentries),stat=info)
            call ArrayStart('tmp-mat',info,size(mat),kind(mat))
            !
            ientry = 0
            !
            do ilevel = 1, nlevels
              !
              if (eigen(ilevel)%jval/=jrot.or.eigen(ilevel)%igamma/=isym) cycle
              !
              ientry = ientry + 1
              !
              !call Hamiltonian_vector(0,ilevel,jrot,isym,nlevels,poten_me,mat(ientry,ientry:Nentries))
              !
              mat(ientry,1:ientry-1) = mat(1:ientry-1,ientry)
              !
            enddo
            !
            ! find the smallest absolute diagonal value of mat
            !
            Smallest  = huge(1.0_rk)
            k  = 1
            !
            do ientry=1,Nentries
               if (abs(mat(ientry,ientry))<=Smallest) then 
                   Smallest = abs( mat(ientry,ientry) )
                   k = ientry
               endif
            enddo
            !
            if (allocated(energy_) ) deallocate(energy_)
            allocate (energy_(Nentries),stat=alloc)
            !
            vrange(1) = -0.0_rk ; vrange(2) = job%upper_ener+mat(k,k)
            !
            rng  = 'V'
            call lapack_syevr(mat(1:Nentries,1:Nentries),energy_(1:Nentries),rng,iroots=nroots,vrange=vrange,irange=irange)
            !
            deallocate(mat)
            !
            ! Zero point energy:
            !
            if (jrot==0.and.isym==1) ezero_ = energy_(1)
            !
            enerright = 0 
            !
            do k = 1,en_npts
              !
              iJ = fitting%obs(k)%Jrot ; gamma = fitting%obs(k)%symmetry
              !
              if (iJ==Jrot.and.gamma==isym) then 
                !
                iener  = fitting%obs(k)%N ;  if (iener>nroots.or.gamma>sym%Nrepresen) cycle ! iener  = 1
                !
                enerright(k) = energy_(iener)-ezero_
                !
              endif 
              !
            enddo ! --- k
            ! 
          enddo 
          !
          potparam(i)=tempx-deltax

          do isym = 1,sym%Nrepresen
            !
            ! Determine the size of the Hamiltonian
            !
            Nentries = 0
            !
            do ilevel = 1, nlevels
              !
              if (eigen(ilevel)%jval==jrot.and.eigen(ilevel)%igamma==isym) Nentries = Nentries + 1
              !
            enddo
            !
            ! Construct the Hamiltonian matrix vector by vector
            !
            if (Nentries<1) cycle
            !
            allocate (mat(Nentries,Nentries),stat=alloc)
            !
            ientry = 0
            !
            do ilevel = 1, nlevels
              !
              if (eigen(ilevel)%jval/=jrot.or.eigen(ilevel)%igamma/=isym) cycle
              !
              ientry = ientry + 1
              !
              !call Hamiltonian_vector(0,ilevel,jrot,isym,nlevels,poten_me,mat(ientry,ientry:Nentries))
              !
              mat(ientry,1:ientry-1) = mat(1:ientry-1,ientry)
              !
            enddo
            !
            ! find the smallest absolute diagonal value of mat
            !
            Smallest  = huge(1.0_rk)
            k  = 1
            !
            do ientry=1,Nentries
               if (abs(mat(ientry,ientry))<=Smallest) then 
                   Smallest = abs( mat(ientry,ientry) )
                   k = ientry
               endif
            enddo
            !
            if (allocated(energy_) ) deallocate(energy_)
            allocate (energy_(Nentries),stat=alloc)
            !
            vrange(1) = -0.0_rk ; vrange(2) = job%upper_ener+mat(k,k)
            !
            rng  = 'V'
            call lapack_syevr(mat(1:Nentries,1:Nentries),energy_(1:Nentries),rng,iroots=nroots,vrange=vrange,irange=irange)
            !
            deallocate(mat)
            call ArrayStop('tmp-mat')
            !
            ! Zero point energy:
            !
            if (jrot==0.and.isym==1) ezero_ = energy_(1)
            !
            enerleft = 0
            !
            do k = 1,en_npts
              !
              iJ = fitting%obs(k)%Jrot ; gamma = fitting%obs(k)%symmetry
              !
              if (iJ==Jrot.and.isym==gamma) then 
                !
                iener  = fitting%obs(k)%N ;  if (iener>nroots.or.isym>sym%Nrepresen) cycle ! iener  = 1
                !
                enerleft(k) = energy_(iener)-ezero_
                !
              endif 
              !
            enddo ! --- k
            ! 
          enddo 
          !
          rjacob(1:en_npts,ncol)=(enerright(1:en_npts)-enerleft(1:en_npts))/(2.0d0*deltax)
          !
          potparam(i)=tempx
          !
        endif
      enddo ! --- ncol
      !
      deallocate (enerright,enerleft,energy_)
      call ArrayStop('enerright-enerleft')
      call ArrayStop('enerright-enerleft')
      !
    end subroutine finite_diff_of_energy
    !
 end subroutine sf_fitting


 subroutine bandcentres_fitting(Jval)

      integer(ik),intent(in)  :: Jval(:)

      integer(ik)    :: npts,en_npts,nused      ! number of data points: current, all obs. energies, all pot. points, and actuall used. 
      integer(ik)       :: parmax                  ! total number of potential parameters 
      integer(ik)       :: jrot                    ! current value of the rotational quantum number 
      integer(ik)       :: nmodes,ncoords,lwork,fititer,nclasses
      integer           :: alloc,info, ierror      ! error state variables 
      integer(ik)       :: j0ener_total            ! total number of the J=0 term values 
      real(rk)          :: wtsum

      real(rk),allocatable :: wt_bit(:),wtall(:) ! weight factors - only for the energies and total.
      real(rk),allocatable :: deriv0(:),rjacob(:,:),eps(:),energy_(:),enercalc(:),ener_j0(:,:)
      real(rk),allocatable :: local(:,:),wspace(:),Tsing(:),chi_(:,:)
      real(rk),allocatable :: al(:,:),bl(:),dx(:),ai(:,:),sterr(:),sigma(:)
      real(rk),allocatable :: mat(:,:),dmat(:,:),vector(:)
      character(len=cl),allocatable :: nampar(:)    ! parameter names 
      !
      integer(ik),allocatable :: ivar(:),ilargest(:),ifitparam(:),ilargest_j0(:,:)
      !
      logical      :: still_run,do_deriv,deriv_recalc
      real(rk) :: stadev_old,stability,stadev,tempx,deltax,v,potright,potleft,sum_sterr,conf_int
      real(rk) :: ssq,rms,ssq1,ssq2,rms1,rms2,Smallest,Largest
      real(rk) :: ezero_,a_wats = 1.0_rk
      integer(ik)  :: iener1,iener2,i,numpar,itmax,j,jlistmax,rank,ncol,nroots
      integer(ik)  :: iJ,isym,jsym,iener,jener,irow,icolumn,l,ndigits,jlevel
      integer(ik)  :: potunit,enunit,quanta(1:FLNmodes)
      integer(ik)  :: nlevels,ilevel,Nentries,ientry,jentry,k,imode,iunit,irange(2),icontr,icase,icase1,j0ener_total_plus,imu_t
      integer(ik)  :: jcase,jlambda,jcontr,j0ener_max,nroots_j0(sym%Nrepresen),j0Neners,ilambda,j0ener_fit_plus,i0,iter_th
      real(rk)     :: vrange(2),f_t,v_l,v_r
      real(ark)    :: ar_t(1:molec%ncoords)
      real(rk)     :: xi(1:FLNmodes)
      real(ark)    :: chi(1:FLNmodes)
      character(cl) :: filename, ioname, buf,dir
      integer(ik)   :: ipowers(1:FLNmodes),nterms,igamma,iterm
      integer(hik)  :: matsize,matsize2
      real(rk),allocatable  ::  deriv(:,:)   ! to store the calculated derivatives 
      real(rk),allocatable  :: pot_matrix(:,:)
      real(rk),allocatable :: xparam(:),cnu_j(:)
      logical       :: do_symmetry
      character(len=1)   :: rng
      character(len=1),allocatable  :: mark(:)
      character(len=cl)             :: my_fmt0 !format for I/O specification
      character(len=cl) :: my_fmt,my_fmt_en3,my_fmt_en4 !format for I/O specification
      character(len=cl) :: my_fmt_en1,my_fmt_en2,my_fmt_par1,my_fmt_par2 !format for I/O specification
       !
       if (job%verbose>=2) write(out,"(/'The least-squares fitting ...')")
       !
       call TimerStart('Band centres fitting')
       !
       en_npts = fitting%Nenergies
       !
       ! file (temp) where all computed energies will be printed out 
       !
       filename =  trim(fitting%output_file)//'.en'
       write(ioname, '(a, i4)') 'calc. energies from the fit '
       call IOstart(trim(ioname), enunit)
       open(unit = enunit, action = 'write',status='replace' , file = filename)
       !
       write(out,"(/'Number of obs. energies: ',i9)") en_npts
       !
       npts = en_npts
       !
       nmodes  = FLNmodes
       !
       ncoords = molec%ncoords
       !
       nclasses = size(eigen(1)%cgamma)-1
       !
       itmax = fitting%itermax
       !
       jlistmax = size(Jval)
       !
       ! number fitting j0energies
       !
       parmax = j0fit%Nenergies
       !
       ! number of fitting j0energies incl. degeneracies
       !
       j0ener_fit_plus = 0 
       do i=1,j0fit%Nenergies 
          do ilambda = 1, bset_contr(1)%index_deg(i)%size1
            j0ener_fit_plus = j0ener_fit_plus + 1
          enddo
       enddo
       !
       ! total number of j=0 energies incl. degeneracies
       !
       j0ener_total_plus = bset_contr(1)%Maxcontracts
       !
       ! the maximal number of energies for the given symmetry
       !
       !j0ener_total = 0 
       j0ener_max = 0
       do isym= 1,sym%Nrepresen
         !j0ener_total = j0ener_total + fit(isym,1)%Nentries
         j0ener_max = max(j0ener_max,fit(isym,1)%Nentries)
       enddo
       !
       allocate (ener_j0(sym%Nrepresen,j0ener_max),ilargest_j0(sym%Nrepresen,j0ener_max),stat=info)
       call ArrayStart('ener_j0',info,size(ener_j0),kind(ener_j0))
       !
       ener_j0 = 0 
       ilargest_j0 = 1
       !
       allocate (enercalc(en_npts),stat=info)
       call ArrayStart('obs-en-mat',info,size(enercalc),kind(enercalc))
       allocate (mark(en_npts),stat=info)
       allocate (wtall(npts),stat=info)
       call ArrayStart('weight-mat',info,size(wtall),kind(wtall))
       allocate (wt_bit(npts),stat=info)
       call ArrayStart('weight-mat',info,size(wt_bit),kind(wt_bit))
       !
       allocate (xparam(parmax),nampar(parmax),ivar(parmax),ifitparam(parmax),&
                 al(parmax,parmax),ai(parmax,parmax),bl(parmax),dx(parmax),sterr(parmax),&
                 Tsing(parmax),cnu_j(nclasses),stat=info)
       call ArrayStart('potparam-mat',info,size(xparam),kind(xparam))
       call ArrayStart('potparam-mat',info,size(ivar),kind(ivar))
       call ArrayStart('potparam-dx',info,size(dx),kind(dx))
       call ArrayStart('potparam-sterr',info,size(sterr),kind(sterr))
       call ArrayStart('potparam-Tsing',info,size(Tsing),kind(Tsing))
       call ArrayStart('a-b-mat',info,size(al),kind(al))
       call ArrayStart('a-b-mat',info,size(bl),kind(bl))
       !
       ! Allocate objects, that will be used for the fitting procedure:
       !
       allocate (deriv0(parmax),rjacob(npts,parmax),eps(npts),stat=info)
       call ArrayStart('deriv0',info,size(deriv0),kind(deriv0))
       call ArrayStart('rjacob',info,size(rjacob),kind(rjacob))
       call ArrayStart('eps',info,size(eps),kind(eps))
       !
       if (fitting%robust>small_) then
         !
         allocate (sigma(npts),stat=info)
         call ArrayStart('eps',info,size(sigma),kind(sigma))
         !
       endif 
       !
       ! The last object to allocate - the lapack related work array
       !
       lwork = 50*parmax
       !
       allocate (wspace(lwork),stat=info)
       call ArrayStart('wspace',info,size(wspace),kind(wspace))
       !
       wtall = 0 
       !
       ! Objects with energies and derivatives at the given fitting iteration 
       !
       write(my_fmt0,'(a,i0,a)') "(,",nmodes,"i3)"
       write(my_fmt,'(a,i0,a)') "(3i5,2x,a3,1x,3f13.4,2x,e9.2,2x,a1,i3,a1,1x,a1,",nmodes,"(1x, i3),a1,a)"
       write(my_fmt_en1,'(a,i0,a,i0,a,i0,a)') "(3i5,2x,a3,1x,3f13.4,2x,e9.2,2x,a3,a1,i3,a2,1x,a2,",&
                                              nclasses,"a3,a1,",nmodes,"(1x, i3),a2,1x,a1,",nmodes,"(1x, i3),a1,a)"
                                              !
       write(my_fmt_en2,'(a,i0,a,i0,a)') "(3i5,2x,a3,1x,3f13.4,2x,e9.2,2x,a2,a3,a1,i3,a2,1x,a2,",&
                                               nclasses,"a3,a1,",nmodes,"(1x, i3),a2)"

       write(my_fmt_en3,'(a,i0,a,i0,a,i0,a)') "(3i5,2x,f18.10,2x,",nmodes+1,"(1x, i3),f8.2,1x,a1,",nmodes+1,&
                         "(1x, i3),a1,2x,a1,",nmodes+1,"(1x, i3),a1,f12.4,2x,f12.4)"

       write(my_fmt_en4,'(a,i0,a)') "(3i5,2x,f18.10,2x,",nmodes+1,"(1x, i3),f8.2)"

       !write(my_fmt_pot2,'(a1,i0,a)') "(",ncoords,"(2x,f18.9),3(x,g18.10),x,e12.4)"
       !write(my_fmt_par1,'(a1,i0,a)') "(a8,4x,",Ncoords,"i3,1x,i2,e22.14)"

       !
       do i = 1,j0fit%Nenergies
         !
         xparam(i) = j0fit%obs(i)%energy
         ivar(i)   = j0fit%obs(i)%weight
         !
         write(nampar(i),my_fmt0) j0fit%obs(i)%quanta(1:nmodes)
         !
       enddo
       !
       forall(i=1:en_npts) wtall(i) = fitting%obs(i)%weight
       wtsum = sum(fitting%obs(1:en_npts)%weight) ; wtsum = max(wtsum,sqrt(small_))
       wtall(1:en_npts) = wtall(1:en_npts)/max(wtsum,sqrt(small_))
       ! 
       ! Count how many data points actually will be fitted. 
       !
       nused=0
       wt_bit = 0 
       !
       do i=1,npts
         if (wtall(i)>small_) then 
           !
           nused=nused+1
           wt_bit(i) = 1.0_rk
           !
         endif
       enddo
       !
       ! sigma = exp. data precision for robust fit 
       !
       if (fitting%robust>small_) then
         !
         sigma = 1.0_rk
         !
         do i=1,npts
           if (wtall(i)>small_) sigma(i) = sigma(i)/sqrt(wtall(i))
         enddo
         !
         wtsum = 1.0_rk ! sqrt(sum(sigma(1:en_npts)**2))
         !
         sigma(:) = sigma(:)*fitting%robust/wtsum
         !
       endif 
       !
       ! assigning mark with null
       !
       mark(:) = ' '
       !
       write(out,"('Number of data points used in the fit: ',i9)") nused
       !
       ! The fitting loop is about to start. 
       !
       ! fititer will count the iterations. 
       !
       fititer = 0
       !
       nlevels = Neigenlevels
       !
       ! The outer fitting loop - allows to restart the fitting in case 
       ! we decide to remove some of the varying parameters. 
       ! This option is working together with fit_type ='linur'
       !
       still_run = .true.
       outer_loop: do while (still_run)  
          !
          ! Initial values for the standard error and  stability.
          !
          stadev_old = 1e10
          stability = 1e10
          stadev    = 1e10
          !
          ! Count the actually varying number of parameters:
          !
          numpar  = 0
          ifitparam = 1
          do i=1,parmax
            if (ivar(i) > 0) then 
              numpar=numpar+1
              ifitparam(numpar) = i
            endif
          enddo 
          !
          if (numpar==0.and.itmax/=0) then 
            !
            write (out,"('No varying paramters, check input!')") 
            stop 'bandcentres_fitting: No varying paramters'
            !
          endif 
          !
          rjacob = 0 
          !xparam = 0
          !
          ! The loop starts here. 
          !
          do while( fititer.le.itmax .and. stadev.ge.fitting%target_rms.and. stability.ge.stab_best)
            !
            fititer = fititer + 1
            !
            write(out,"(/'Iteration = ',i8)") fititer
            write(enunit,"(/'Iteration = ',i8)") fititer
            !
            eps = 0
            !
            ! Print out the calc. and obs.-calc., i.e. result of the fit. 
            ! Only fitted energies are printed. 
            !
            write(out,"(/1X,100('-'),/a,/1X,100('-'))") &
                       '|## |  N |  J | sym|     Obs.    |    Calc.   | Obs.-Calc. |   Weight |    K     vib. quanta'
            !
            write(enunit,"(/1X,100('-'),/a,/1X,100('-'))") &
                         '|## |  N |  J | Sym|     Obs.    |    Calc.   | Obs.-Calc. |   Weight |    K    quanta   (Calc./Obs.) '
            !
            do j = 1,jlistmax
               !
               jrot = jval(j)
               !
               deriv_recalc = .true.
               !
               loop_isym : do isym = 1,sym%Nrepresen
                 !
                 iunit = fit(isym,j)%IOunit
                 !
                 Nentries = fit(isym,j)%Nentries
                 !
                 if (Nentries<1) cycle
                 !
                 ! Check if we don't need to process this pair of J,isym 
                 !
                 do_symmetry = .false.
                 do i = 1,en_npts
                   if (fitting%obs(i)%Jrot==Jrot.and.fitting%obs(i)%symmetry==isym) then 
                     do_symmetry = .true.
                     exit
                   endif
                 enddo
                 if (.not.do_symmetry) cycle
                 !
                 if (job%verbose>=6) then 
                    write (out,"(/'Matrix size = ',i8)") Nentries
                 endif
                 !
                 matsize = int(Nentries*(Nentries+1)/2,hik)
                 matsize2= int(Nentries*Nentries,hik)
                 !
                 allocate(pot_matrix(Nentries,Nentries),stat=alloc)
                 call ArrayStart('pot_matrix',alloc,1,kind(pot_matrix),matsize2)
                 !
                 call TimerStart('Energy calculations')
                 !
                 allocate (mat(Nentries,Nentries),stat=info)
                 call ArrayStart('fitting-mat',info,1,kind(mat),matsize2)
                 !
                 mat = 0
                 !
                 do ientry = 1, Nentries
                   !
                   ilevel = fit(isym,j)%ilevel(ientry)
                   mat(ientry,ientry) = eigen(ilevel)%energy
                   !
                 enddo 
                 !
                 rewind(iunit) 
                 read(iunit) buf(1:13)
                 !
                 do icontr=1,j0ener_fit_plus
                   !
                   icase   = bset_contr(1)%icontr2icase(icontr,1)
                   ilambda = bset_contr(1)%icontr2icase(icontr,2)
                   !
                   if (ilambda==1) icase1 = icase
                   !
                   if (icase1>j0fit%Nenergies) cycle
                   !
                   read(iunit) imu_t
                   !
                   if (icontr/=imu_t) then
                     write (out,"(' evib_matrix at J = ',i4,', isym = ',i4,' i <> j: ',2i8)") jrot,isym,icontr,imu_t
                     stop 'evib_matrix - i<>j'
                   end if
                   !
                   read (iunit) pot_matrix
                   !
                   !read (iunit,rec=icontr) pot_matrix(1:matsize)
                   !
                   !$omp parallel do private(ientry,jentry) shared(mat) schedule(guided)
                   do ientry = 1, Nentries
                     !
                     do jentry = 1, ientry
                       !
                       mat(ientry,jentry) = mat(ientry,jentry) + xparam(icase1)*pot_matrix(ientry,jentry)
                       !
                       mat(jentry,ientry) = mat(ientry,jentry)
                       !
                       if (fit_debug > 3) then
                          write (out,"('mat (',i0,',',i0,')= ',es14.7)") ientry,jentry,mat(ientry,jentry)
                       endif
                       !
                     enddo
                     !
                   enddo
                   !$omp end parallel do
                   !
                 enddo
                 !
                 !
                 ! find the smallest absolute diagonal value of mat
                 !
                 Smallest  = huge(1.0_rk)
                 k  = 1
                 !
                 do ientry=1,Nentries
                    if (abs(mat(ientry,ientry))<=Smallest) then 
                        Smallest = abs( mat(ientry,ientry) )
                        k = ientry
                    endif
                 enddo
                 !
                 if (fit_debug > 2) then
                    !
                    write (out,"(/'Smallest diag. value of mat = ',es14.7,'at k = ',i8,' upper_range = ',es14.7/)") &
                                   mat(k,k),k,job%upper_ener+mat(k,k)
                    !
                 endif
                 !
                 if (allocated(energy_)) deallocate(energy_)
                 !
                 allocate (energy_(Nentries),stat=info)
                 call ArrayStart('fit-ener-deriv',info,size(energy_),kind(energy_))
                 !
                 select case (trim(job%diagonalizer)) 

                 case default
                   !
                   write (out,"('bandcentres_fitting: type of the diagonalizer  ',a,' unknown')") trim(job%diagonalizer)
                   stop 'bandcentres_fitting - wrong diagonalizer '
                   !
                 case('SYEV') 
                   !
                   call lapack_syev(mat(1:Nentries,1:Nentries),energy_(1:Nentries))
                   nroots = Nentries
                   !
                 case('SYEVD') 
                   !
                   call lapack_syevd(mat(1:Nentries,1:Nentries),energy_(1:Nentries))
                   nroots = Nentries
                   !
                 case('SYEVR') 
                   !
                   vrange(1) = -0.0_rk ; vrange(2) = job%upper_ener+mat(k,k)
                   irange(1) = 1 ; irange(2) = min(job%nroots(isym),Nentries)
                   !
                   if (job%nroots(isym)/=1e6) then
                      rng = 'I'
                   elseif (job%upper_ener/=1e9) then 
                      rng = 'V'
                   else
                      rng = 'A'
                   endif 
                   !
                   call lapack_syevr(mat(1:Nentries,1:Nentries),energy_(1:Nentries),rng,iroots=nroots,vrange=vrange,irange=irange)
                   !
                 end select
                 !
                 ! Zero point energy:
                 !
                 if (jrot==0.and.isym==1) ezero_ = energy_(1)
                 !
                 ! find the largest coefficient for each root for the assignment
                 !
                 allocate(ilargest(nroots),stat=info)
                 call ArrayStart('ilargest',info,size(ilargest),kind(ilargest))
                 !
                 do i = 1,nroots
                    !
                    Largest  = -1.0
                    k  = 1
                    !
                    do ientry=1,Nentries
                       if (abs(mat(ientry,i))>=Largest) then 
                           Largest = abs( mat(ientry,i) )
                           k = ientry
                       endif
                    enddo
                    !
                    ilargest(i) = k
                    !
                    if (jrot==0) then 
                      !
                      ener_j0(isym,i) = energy_(i) - ezero_
                      ilargest_j0(isym,i) = k
                      nroots_j0(isym) = nroots
                      !
                    endif
                    !
                 enddo
                 !
                 call TimerStop('Energy calculations')
                 !
                 ! Check we if we need to computed the derivatives 
                 !
                 do_deriv = .false.
                 !
                 if (itmax.ge.1.and.trim(deriv_type)=='hellman'.and.mod(fititer+2,3)==0) do_deriv = .true.
                 !
                 if (do_deriv) then 
                   !
                   call TimerStart('Energy derivatives')
                   !
                   allocate (deriv(nroots,numpar),stat=info)
                   call ArrayStart('fit-ener-deriv',info,size(deriv),kind(deriv))
                   !
                   deriv = 0
                   !
                   allocate (dmat(Nentries,Nentries),stat=info)
                   call ArrayStart('dmat',info,size(dmat),kind(dmat))
                   !
                   ncol = 0
                   rewind(iunit) 
                   read(iunit) buf(1:13)
                   !
                   do icontr=1,j0ener_fit_plus
                      !
                      icase   = bset_contr(1)%icontr2icase(icontr,1)
                      ilambda = bset_contr(1)%icontr2icase(icontr,2)
                      !
                      !if (ilambda==1) icase1 = icase
                      !
                      if (icase>j0fit%Nenergies) cycle
                      !
                      read(iunit) imu_t
                      !
                      read(iunit) pot_matrix
                      !
                      if (ivar(icase) <= 0) cycle
                      !
                      if (ilambda==1) ncol = ncol + 1
                      !
                      !do ncol=1,numpar 
                      !
                      !icase = ifitparam(ncol)
                      !
                      !do ilambda = 1,bset_contr(1)%index_deg(icase)%size1
                      !
                      !icontr = bset_contr(1)%icontr_correlat_j0(icase,ilambda)
                      !
                      !if (icontr>j0ener_fit_plus) cycle
                      !
                      !read(iunit,rec=icontr) pot_matrix(1:matsize)
                      !
                      !$omp parallel do private(ientry,jentry) shared(dmat) schedule(guided)
                      do ientry = 1, Nentries
                        !
                        do jentry = 1,ientry
                          !
                          dmat(ientry,jentry) = dmat(ientry,jentry)+pot_matrix(ientry,jentry)
                          dmat(jentry,ientry) = dmat(ientry,jentry)
                          !
                        enddo
                        !
                      enddo
                      !$omp end parallel do
                      !
                      !enddo
                      !
                      !$omp parallel private(vector,info) shared(deriv,deriv0)
                      allocate (vector(Nentries),stat=info)
                      if (info/=0)  stop 'fitting-vector - out of memory'
                      !
                      !$omp do private(iener,iJ,igamma,ientry) schedule(guided)
                      do iener = 1,en_npts
                        !
                        iJ = fitting%obs(iener)%Jrot ; igamma = fitting%obs(iener)%symmetry
                        !
                        if (iJ==Jrot.and.igamma==isym) then 
                          !
                          ientry  = fitting%obs(iener)%N ;  if (ientry>nroots.or.isym>sym%Nrepresen) cycle 
                          !
                          call dgemv('N',Nentries,Nentries,1.0d0,dmat,Nentries,mat(:,ientry),1,0.0d0,vector,1)
                          !
                          deriv(ientry,ncol)= dot_product(mat(:,ientry),vector(:))
                          !
                        endif
                        !
                      enddo
                      !$omp end do
                      !
                      deallocate(vector)
                      !$omp end parallel
                      !
                      if (jrot==0.and.isym==1) then 
                        !
                        if (fitting%obs(1)%N==1.and.1==fitting%obs(1)%Jrot.and.fitting%obs(1)%symmetry==0) then
                           !
                           deriv0(ncol) = deriv(1,ncol)
                           !
                        else
                           !
                           allocate (vector(Nentries),stat=info)
                           if (info/=0)  stop 'fitting-vector - out of memory'
                           !
                           call dgemv('N',Nentries,Nentries,1.0d0,dmat,Nentries,mat(:,1),1,0.0d0,vector,1)
                           deriv0(ncol)= dot_product(mat(:,1),vector(:))
                           !
                           deallocate(vector)
                           !
                        endif
                      endif
                      !
                   enddo 
                   !
                   deriv_recalc = .false.
                   !
                   call TimerStop('Energy derivatives')
                   !
                   deallocate(dmat)
                   call ArrayStop('dmat')
                   !
                 endif 
                 !
                 deallocate(pot_matrix)
                 call ArrayStop('pot_matrix')
                 !
                 deallocate(mat)
                 call ArrayStop('fitting-mat')
                 !
                 ! correct the addresses of the energy
                 !
                 if (fitting%threshold_lock>0) then 
                   !
                   i0 = 1
                   !
                   do iener = 1,en_npts
                     !
                     iJ = fitting%obs(iener)%Jrot ; igamma = fitting%obs(iener)%symmetry
                     mark(iener) = '*'
                     !
                     if (iJ==Jrot.and.igamma==isym) then
                       !
                       loop_thresh : do iter_th = 1,maxiter_as
                         !
                         loop_i : do i=i0,nroots
                            !
                            do jener = max(iener-10,1),iener-1
                              !
                              if (fitting%obs(jener)%N == i) cycle loop_i
                              !
                            enddo
                            !
                            quanta(1:nmodes) = eigen(fit(isym,j)%ilevel(ilargest(i)))%quanta(1:nmodes)
                            !
                            if (abs( fitting%obs(iener)%energy-(energy_(i)-ezero_) )<=real(iter_th,8)*fitting%threshold_lock.and.&
                                all(quanta(1:nmodes)==fitting%obs(iener)%quanta(1:nmodes))) then 
                                !
                                fitting%obs(iener)%N = i
                                mark(iener) = ' '
                                exit loop_thresh
                                !
                            endif 
                            !
                         enddo loop_i
                         !
                       enddo loop_thresh
                       !
                       if (mark(iener) == ' ')  i0 = max(fitting%obs(iener)%N-10,1)
                       !
                     endif 
                     !
                   enddo
                   !
                 endif
                 !
                 write(my_fmt,'(a,i0,a)') "(3i5,2x,a3,1x,3f13.4,2x,e9.2,2x,a1,i3,a1,1x,a1,",nmodes,"(1x, i3),a1,a)"
                 !
                 do iener = 1,en_npts
                   !
                   iJ = fitting%obs(iener)%Jrot ; igamma = fitting%obs(iener)%symmetry
                   !
                   if (iJ==Jrot.and.igamma==isym) then 
                     !
                     i  = fitting%obs(iener)%N ;  if (i>nroots.or.isym>sym%Nrepresen) cycle 
                     !
                     enercalc(iener) = energy_(i)-ezero_
                     !
                     eps(iener) = fitting%obs(iener)%energy-enercalc(iener)
                     !
                     if (do_deriv) rjacob(iener,1:numpar) = deriv(i,1:numpar)-deriv0(1:numpar)
                     !
                     ilevel = fit(isym,j)%ilevel(ilargest(i))
                     !
                     write (out,my_fmt) &
                            iener,i,iJ,sym%label(isym),enercalc(iener)+eps(iener),enercalc(iener),-eps(iener),&
                            wtall(iener),&
                            eigen(ilevel)%krot,eigen(ilevel)%quanta(1:nmodes),mark(iener)
                            !fitting%obs(i)%quanta(1:nmodes)
                     !
                   endif 
                   !
                 enddo ! --- i
                 !
                 ! Printing all calculated term values. If the obs. counterpats exist, 
                 ! the obs.-calc. are printed as well. 
                 ! This list can be used to identify the obs-calc pairs, i.e. the number of 
                 ! a term value as it appear in the list of calculated energies. 
                 !
                 ! write(out,"(/'Iteration = ',i4)") fititer
                 !
                 do i = 1,nroots
                    !
                    if (i==1.or.energy_(i)-ezero_>sqrt(small_)) then 
                       !
                       ilevel = fit(isym,j)%ilevel(ilargest(i))                       
                       !
                       jener = 1
                       !
                       do while (jener/=en_npts+1.and..not.(i==fitting%obs(jener)%N.and.Jrot==fitting%obs(jener)%Jrot.and.&
                                 isym==fitting%obs(jener)%symmetry) )
                         !
                         jener = jener+1
                         !
                       enddo
                       !
                       if (jener<en_npts) then 
                         !
                         write(enunit,my_fmt_en1) &
                            i,fitting%obs(jener)%N,Jrot,&
                            sym%label(isym),enercalc(jener)+eps(jener),&
                            enercalc(jener),eps(jener),wtall(jener),&
                            '( ',eigen(ilevel)%cgamma(0),';',eigen(ilevel)%krot,' )',&
                            '( ',eigen(ilevel)%cgamma(1:nclasses),';',eigen(ilevel)%quanta(1:nmodes),' )',&
                            '(',fitting%obs(jener)%quanta(1:nmodes),')',mark(jener)
                         !
                       else
                         !
                         !write(enunit,"(4i5,' ',3f13.4,2x,e8.2,5x,<nmodes>(i3))") i,0,Jrot,isym,0.0,energy_(i)-ezero_,0.0,0.0
                         !
                         write(enunit,my_fmt_en2) &
                            i,0,Jrot,sym%label(isym),0.0,energy_(i)-ezero_,0.0,0.0,&
                            '( ',eigen(ilevel)%cgamma(0),';',eigen(ilevel)%krot,' )',&
                            '( ',eigen(ilevel)%cgamma(1:nclasses),';',eigen(ilevel)%quanta(1:nmodes),' )'
                         !
                       endif 
                       !
                    endif 
                    !
                 enddo
                 !
                 deallocate(energy_)
                 if (allocated(deriv)) deallocate(deriv )
                 call ArrayStop('fit-ener-deriv')
                 !
                 deallocate(ilargest)
                 call ArrayStop('ilargest')
                 !
               enddo loop_isym
               !
            enddo 
            !
            ! ssq  - weighted rms**2, rms  - root mean square deviation. 
            !
            ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
            rms=sqrt(sum(eps(1:npts)*eps(1:npts))/npts)
            !
            ! Prepare the linear system a x = b as in the Newton fitting approach.  
            !
            if (itmax>=1) then
               !----- form the a and b matrix ------c
               ! form A matrix 
               do irow=1,numpar       
                 do icolumn=1,irow    
                   al(irow,icolumn)=sum(rjacob(1:npts,icolumn)*rjacob(1:npts,irow)*wtall(1:npts))
                   al(icolumn,irow)=al(irow,icolumn)
                   if (fit_debug > 2) then
                     write (out,"('al (',i0,',',i0,')= ',es14.7)") irow,icolumn,al(irow,icolumn)
                   endif
                 enddo
               enddo
               !
               ! form B matrix 
               do irow=1,numpar      
                 bl(irow)=sum(eps(1:npts)*rjacob(1:npts,irow)*wtall(1:npts))
                 if (fit_debug > 2) then
                   write (out,"('bl (',i0,')= ',es14.7)") irow,bl(irow)
                 endif
               enddo  
               !
               ! Two types of the linear solver are availible: 
               ! 1. linur (integrated into the program, from Ulenikov Oleg)
               ! 2. dgelss - Lapack routine (recommended).
               !
               select case (trim(fitting%fit_type)) 
               !
               case default
                 !
                 write (out,"('fit_type ',a,' unknown')") trim(fitting%fit_type)
                 stop 'fit_type unknown'
                 !
               case('LINUR') 
                 !
                 call MLlinur(numpar,numpar,al(1:numpar,1:numpar),bl(1:numpar),dx(1:numpar),ierror)
                 !
                 ! In case of dependent parameters  "linur" reports an error = ierror, 
                 ! which is a number of the dependent parameter. We can remove this paramter 
                 ! from the fit and set its value to zero. And start the iteration again. 
                 !
                 if (ierror.ne.0) then 
                   do ncol=1,numpar 
                      !
                      i = ifitparam(ncol)
                      if  ( ncol.eq.ierror ) then 
                          ivar(i) = 0
                          xparam(i) = 0
                          !
                          write(out,"(i0,'-th is out - ',a8)") i,nampar(i)
                          !
                      endif 
                      !
                   enddo 
                   !
                   numpar  = 0
                   ifitparam = 1
                   do i=1,parmax
                     if (ivar(i) > 0) then 
                       numpar=numpar+1
                       ifitparam(numpar) = i
                     endif
                   enddo 
                   !
                   cycle outer_loop    
                  endif 
                  !
               case ('DGELSS')
                 !
                 ai = al 
                 call dgelss(numpar,numpar,1,ai(1:numpar,1:numpar),numpar,bl(1:numpar),numpar,Tsing,1.D-12,rank,wspace,lwork,info)
                 !
                 if (info/=0) then
                   write(out,"('dgelss:error',i0)") info
                   stop 'dgelss'
                 endif
                 !
                 dx = bl
                 !
               end select 
               !
               !----- update the parameter values ------!
               !
               do ncol=1,numpar 
                  i = ifitparam(ncol)
                  xparam(i)=xparam(i)+dx(ncol)
               enddo
               !
               ! Robust fit: adjust the fitting weights
               !
               if (fitting%robust>small_) then
                 !
                 call robust_fit(numpar,a_wats,sigma(1:npts),eps(1:npts),wtall(1:npts))
                 !
                 ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
				 !
				 !fitting%robust = -1.0_rk
                 !
               endif 
               !
               ! Estimate standard deviation error. 
               !
               if ( nused.ne.numpar ) then 
                 stadev=sqrt(ssq/float(nused-numpar))
               else 
                 stadev=sqrt(ssq/nused)
               endif
               !
               ! Estimate the standard errors for each parameter using 
               ! the inverse matrix of a. 
               !
               call MLinvmat(al(1:numpar,1:numpar),ai(1:numpar,1:numpar),numpar,info)
               !
               sum_sterr=0.d0
               do ncol=1,numpar 
                  i = ifitparam(ncol)
                  if (nused.eq.numpar) then  
                     sterr(ncol)=0
                  else
                     sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
                     sum_sterr=sum_sterr+abs(sterr(ncol)/xparam(i))
                  endif
               enddo    
               !
               sum_sterr=sum_sterr/numpar 
               !
               ! This is how we define stability of the fit:
               ! as a relative change of stadev comparing with the step before. 
               !
               stability=abs( (stadev-stadev_old)/stadev )
               stadev_old=stadev
               !
            else
               !
               stadev=sqrt(ssq/nused)
               !
            endif
            !
            ! Print the updated parameters. 
            !
            write(out,"(/'Correction to band centers:')")
            !
            do i=1,j0fit%Nenergies
              !
              ilevel = j0fit%obs(i)%N ; if (ilevel>nroots.or.j0fit%obs(i)%symmetry>sym%Nrepresen) ilevel = 1e6
              !
              !if (ilevel>nroots.or.j0fit%obs(i)%symmetry>sym%Nrepresen) then
              !  write(out,"('bandcentres_fitting: wrong j0ener-parameters: i,ilevel,isym = ',3i)") i,ilevel,j0fit%obs(i)%symmetry
              !  stop 'bandcentres_fitting: wrong j0ener-parameters'
              !endif
              !
              isym = j0fit%obs(i)%symmetry
              !
              !iener  = fitting%obs(ilevel)%N ;  if (iener>nroots.or.isym>sym%Nrepresen) cycle 
              !
              iener = 1
              !
              do while (iener/=en_npts+1.and..not.(ilevel==fitting%obs(iener)%N.and.fitting%obs(iener)%Jrot==0.and.&
                        fitting%obs(iener)%symmetry==isym) )
                !
                iener = iener+1
                !
              enddo
              !
              !if (iener>en_npts) then
              !  write(out,"('bandcentres_fitting: wrong j0ener- or obs-parameters: i,ilevel,iener = ',3i)") i,ilevel,iener
              !  stop 'bandcentres_fitting: wrong j0ener- or obs-parameters'
              !endif
              !
              if (iener<=en_npts) then 
                !
                jlevel = fit(isym,1)%ilevel(ilargest_j0(isym,ilevel))
                !
                write(out,my_fmt_en3) &
                   Jrot,j0fit%obs(i)%symmetry,j0fit%obs(i)%N,&
                   xparam(i),j0fit%obs(i)%quanta(0:nmodes),real(ivar(i),rk), &
                   '(',eigen(jlevel)%quanta(1:nmodes),')',&
                   '(',fitting%obs(iener)%quanta(0:nmodes),')',&
                   ener_j0(isym,ilevel),fitting%obs(iener)%energy
                   !
              else if (ilevel<=nroots_j0(isym)) then
                !
                write(out,my_fmt_en3) &
                   Jrot,j0fit%obs(i)%symmetry,j0fit%obs(i)%N,&
                   xparam(i),j0fit%obs(i)%quanta(0:nmodes),real(ivar(i),rk), &
                   '(',eigen(jlevel)%quanta(0:nmodes),')',&
                   ener_j0(isym,ilevel)
                   !
              else 
                !
                write(out,my_fmt_en4) &
                   Jrot,j0fit%obs(i)%symmetry,j0fit%obs(i)%N,&
                   xparam(i),j0fit%obs(i)%quanta(0:nmodes),real(ivar(i),rk)
               
              endif
              !
            enddo
            !
            write(out,"(/'Band centers:')")
            !
            do i=1,parmax
              !
              write (out,"(a12,4x,i2,f18.8)") nampar(i),ivar(i),eigen(i)%energy+xparam(i)
              !
            enddo
            !
            ! Output some statistics and results 
            !
            !  only if we are fitting:  
            !
            if (itmax.ne.0) then
              !
              !still_run = .false.
              ! 
              write (out,"(//'Correction to band centers rounded in accord. with their standart errors'/)")
              l = 0 
              do i=1,parmax
                if (ivar(i) .ne. 0) then
                   l=l+1
                   ndigits = 0
                   conf_int = sterr(l)
                   do while (conf_int.le.10.0.and.ndigits<10)
                     ndigits = ndigits +1 
                     conf_int = conf_int*10
                   enddo
                   !
                   if (conf_int>1e8) conf_int = 0 
                   !
                   !write (out,"(a8,i4,2x,f22.<ndigits>,'(',i14,')')") nampar(i),ivar(i),xparam(i),nint(conf_int)
                   !
                   write(my_fmt_par2,'(a1,i0,a)') "(a8,i4,2x,f22.,",ndigits,"a1,i14,a1)"
                   !
                   write (out,my_fmt_par2) nampar(i),ivar(i),xparam(i),'(',nint(conf_int),')'
                   !
                else 
                   ndigits =2
                   if (xparam(i).ne.0.0) ndigits = 8
                   !
                   !write (out,"(a8,i4,2x,f22.<ndigits>)") nampar(i),ivar(i),xparam(i)
                   write(my_fmt_par2,'(a1,i0,a)') "(a8,i4,2x,f22.,",ndigits,")"
                   !
                   write (out,my_fmt_par2) nampar(i),ivar(i),xparam(i)
                   !
                endif
              enddo  ! --- i
              !
            endif 
            !
            still_run = .false.
            !
            ! Print out the ssq for the rovib. energies and pot. data points separetely:
            !
            ssq1 = 0 ; ssq2 = 0 
            !
            wtsum = sum(wt_bit(1:en_npts))
            !
            if (wtsum/=0) ssq1 = sqrt( sum(eps(1:en_npts)**2*dble(wt_bit(1:en_npts)))/wtsum )
            !
            wtsum = sum(wt_bit(1+en_npts:npts))
            !
            if (wtsum/=0) ssq2 = sqrt( sum(eps(1+en_npts:npts)**2*dble(wt_bit(1+en_npts:npts)))/wtsum )
            !
            rms1=sqrt(sum(eps(1:en_npts)**2)/en_npts)
            rms2=rms1
            !
            write (out,6552) fititer,nused,numpar,stadev,ssq1,ssq2,stability
            !
            ! Print the potential energy points into a separate unit. 
            !
            if (job%verbose>=6) call TimerReport
            !
          enddo  ! --- fititer
          !
       enddo outer_loop
       !
       if (allocated(xparam)) deallocate(xparam)
       !
       deallocate (nampar,ivar,ifitparam,al,ai,bl,dx,sterr,Tsing)
       call ArrayStop('potparam-mat')
       call ArrayStop('potparam-dx')
       call ArrayStop('potparam-sterr')
       call ArrayStop('potparam-Tsing')
       call ArrayStop('a-b-mat')
       !
       ! Allocate objects, that will be used for the fitting procedure:
       !
       deallocate (deriv0,rjacob,eps,wspace,ener_j0)
       call ArrayStop('deriv0')
       call ArrayStop('rjacob')
       call ArrayStop('eps')
       call ArrayStop('wspace')
       call ArrayStop('ener_j0')
       !
       if (allocated(sigma)) deallocate(sigma)
       !
       close(enunit,status="keep")
       !
       call TimerStop('Band centres fitting')
       !


6552   format(/3X,80('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    ssq_ener |   ssq_pot | Stability |'/&
       3X,80('-')/,&
       '   | ',i6,' | ',i6,' | ',i6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',E10.3,' |',/3X,80('-')/)

6553   format(/3X,80('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    rms_ener |   rms_pot | Stability |'/&
       3X,80('-')/,&
       '   | ',i6,' | ',i6,' | ',i6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',E10.3,' |',/3X,80('-')/)

    !
 end subroutine bandcentres_fitting



    subroutine Robust_fit(numpar,a_wats,sigma,eps,wt)

      integer(ik),intent(in) :: numpar
      real(rk),intent(inout) :: a_wats
      real(rk),intent(in)    :: sigma(:),eps(:)
      real(rk),intent(inout) :: wt(:)
      !
      integer(ik)            :: npts,i,nrow,nused
      real(rk)               :: da1,da2,wtsum,da
      !
      if (job%verbose>=4) write(out,"(/'Robust fitting ...')")
      !
      npts = size(sigma)
      !
      nused = 0 
      do i=1,npts
        if (wt(i)>small_) nused=nused+1
      enddo
      !
      !Watson alpha-parameter
      ! 
      !do i = 1,-1
      !  !
      !  da1 = 0
      !  da2 = 0
      !  !
      !  do nrow=1,npts
      !    if (wt(nrow)>small_) then 
      !      da1 = da1+eps(nrow)**2/( sigma(nrow)**2+a_wats*eps(nrow)**2 )
      !      da2 = da2+eps(nrow)**4/( sigma(nrow)**2+a_wats*eps(nrow)**2 )**2
      !    endif 
      !  enddo 
      !  !
      !  da =( da1 -  real(nused-numpar,rk) )/da2
      !  a_wats = a_wats + da
      !  !
      !  if (a_wats<sqrt(small_)) a_wats = 1e-3+real(i,rk)*1e-2
      !  !
      !enddo
      !
      !a_wats = 0.001_rk
      !
      if (job%verbose>=4) write(out,"('Watson parameter =',f18.8)") a_wats
      ! 
      !  adjusting the weights  
      ! 
      do nrow=1,npts
         if (wt(nrow)>small_) wt(nrow) = 1.0_rk/( sigma(nrow)**2 + a_wats*eps(nrow)**2 )
      enddo 
      ! 
      ! "re-normalizing" the weight factors
      !
      wtsum = sum(wt(1:npts))
      wt(1:npts) = wt(1:npts)/wtsum
      !
      if (job%verbose>=4) write(out,"('... done!')")
      !
    end subroutine Robust_fit
 
 !
 !
 ! Here we restore the vibrational (J=0) contracted matrix elements of the dipole moment
 !
 subroutine restore_vib_matrix_elements 
   !
   integer(ik)        :: chkptIO,alloc
   character(len=cl)  :: job_is

   character(len=20)  :: buf20
   integer(ik)        :: ncontr_t,rootsize_t,imu,imu_t
   integer(hik)       :: rootsize,rootsize2


   if (fit_debug > 1) then
      write(out, '(/a, 1x, a)') 'read vibrational contracted matrix elements from file', trim(job%extFmat_file)
   endif
   !
   job_is ='external field contracted matrix elements for J=0'
   call IOStart(trim(job_is),chkptIO)
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
   if (trim(job%kinetmat_format)=='OLD') read(chkptIO) rootsize_t
   !
   if (bset_contr(1)%Maxcontracts/=ncontr_t) then
     write (out,"(' Dipole moment checkpoint file ',a)") job%extFmat_file
     write (out,"(' Actual and stored basis sizes at J=0 do not agree  ',2i0)") bset_contr(1)%Maxcontracts,ncontr_t
     stop 'restore_vib_matrix_elements - in file - illegal ncontracts '
   end if
   !
   rootsize = int(ncontr_t*(ncontr_t+1)/2,hik)
   rootsize2= int(ncontr_t*ncontr_t,hik)
   !
   if (fit_debug > 2) then
      write(out,"(/'restore_vib_matrix_elements...: Number of elements: ',i8)") ncontr_t
   endif
   !
   allocate(poten_me(ncontr_t,ncontr_t,extF%rank),stat=alloc)
   call ArrayStart('poten_me',alloc,1,kind(poten_me),rootsize2)

   !if (imu_t/=1) then
   !  write (out,"(' restore_vib_matrix_elements ',a,' has bogus imu - restore_vib_matrix_elements: ',i8,'/=',i8)") imu_t,1
   !  stop 'restore_vib_matrix_elements - bogus imu restore_vib_matrix_elements'
   !end if
   !
   do imu = 1,extF%rank
     !
     read(chkptIO) imu_t
     !
     read(chkptIO) poten_me(:,:,imu)
     !
   enddo
   !
   read(chkptIO) buf20(1:18)
   if (buf20(1:18)/='End external field') then
     write (out,"(' restore_vib_matrix_elements ',a,' has bogus footer: ',a)") job%kinetmat_file,buf20(1:17)
     stop 'restore_vib_matrix_elements - bogus file format'
   end if
   !
   close(chkptIO,status='keep')

   if (fit_debug > 2) then
      write(out, '(/a)') 'done'
   endif

   !
 end subroutine restore_vib_matrix_elements
 !
 ! 
 !
 subroutine Hamiltonian_vector(ientry,indJ,nentries,tmat,cdimen,poten,mat)
   !
   integer(ik),intent(in)  :: ientry,indJ,nentries
   real(rk),intent(in)     :: poten(:,:)
   type(coeffT),intent(in) :: tmat(nentries)
   integer(ik),intent(in)  :: cdimen(nentries)
   real(rk),intent(out)    :: mat(nentries)
   !
   integer(ik)  :: dimen,info
   integer(ik)  :: jentry
   integer(ik)             :: i,irootF,cirootI,icontrF,icontrI,ktauF,irootI
   real(rk),allocatable    :: half_matelem(:)
    !
    !
    !ilevelI = fit(gamma,indJ)%ilevel(ientry)
    !
    ! the maximal size of the basis functions. 
    !
    dimen = bset_contr(indJ)%Maxcontracts
    !
    allocate(half_matelem(dimen),stat=info)
    call ArrayStart('half_matelem',info,size(half_matelem),kind(half_matelem))
    !
    ! Compute the half_matelem
    !
    half_matelem    = 0
    !
    !loop over final state basis components
    !
    !$omp parallel do private(irootF,icontrF,ktauF,cirootI,irootI,icontrI) schedule(guided) ! shared(half_matelem) 
    do irootF = 1, dimen
       !
       icontrF = bset_contr(indJ)%iroot_correlat_j0(irootF)
       !
       ktauF = bset_contr(indJ)%ktau(irootF)
       !
       !loop over initial state basis components
       !
       loop_I : do cirootI = 1, cdimen(ientry)
          !
          irootI = tmat(ientry)%icoeff(cirootI)
          !
          if (ktauF/=bset_contr(indJ)%ktau(irootI)) cycle
          !
          icontrI = bset_contr(indJ)%iroot_correlat_j0(irootI)
          !
          !index of the corresponding vibrational contracted matrix element (cind)
          !
          !irow = max(icontrF, icontrI)
          !icol = min(icontrF, icontrI)
          !cind = irow * (irow - 1) / 2 + icol
          !
          !compute me
          !
          !#if (fit_debug >= 3)
          !  write (out,"('irootF,icontrF,cirootI,icontrI,irow,icol,cindtmat(ientry)%icoef = ',8i,' poten,tmat =  ',2(2x,es18.9))") irootF,icontrF,cirootI,icontrI,irow,icol,cind,tmat(ientry)%icoeff(cirootI),&
          !                            poten(cind),tmat(ientry)%coeff(cirootI)
          !#endif
          !
          half_matelem(irootF) = half_matelem(irootF) + poten(icontrI,icontrF)*tmat(ientry)%coeff(cirootI)
          !
       end do loop_I
       !
    end do
    !$omp end parallel do
    !
    if (fit_debug > 2) then
      write (out,"('ientry = ',i0,'; dimen = ',i0)") ientry,dimen
    endif
    !
    !loop over final states
    !
    mat = 0
    !
    !$omp parallel do private(jentry,cirootI) schedule(guided) ! shared(mat) 
    Flevels_loop: do jentry = 1,ientry
       !
       !ilevelF = fit(gamma,indJ)%ilevel(jentry)
       !
       !mat(jentry) = mat(jentry) + sum( half_matelem(tmat(jentry)%icoeff(:))*tmat(jentry)%coeff(:) )
	   !
	   !cirootI = cdimen(jentry)
	   !
	   !mat(jentry) = dot_product( half_matelem( tmat(jentry)%icoeff( 1:cirootI ) ),tmat(jentry)%coeff( 1:cirootI ) )
       !
       do cirootI = 1, cdimen(jentry)
        !
        mat(jentry) = mat(jentry) + half_matelem(tmat(jentry)%icoeff(cirootI))*tmat(jentry)%coeff(cirootI)
        !
        !#if (fit_debug >= 3)
        !  write (out,"('jentry = ',2i,'; mat,half_matelem,tmat-coef,icoef = ',3es18.8,i0)") & 
        !      jentry,cirootI,mat(jentry),&
        !      half_matelem(tmat(jentry)%icoeff(cirootI)),&
        !      tmat(jentry)%coeff(cirootI),tmat(jentry)%icoeff(cirootI)
        !#endif
        !
       enddo
       !
       !#if (fit_debug >= 3)
       !  write (out,"('jentry = ',i0,'; mat(jentry) = ',es18.8)") jentry,mat(jentry)
       !#endif
       !
    end do Flevels_loop
    !$omp end parallel do
    !
    deallocate(half_matelem)
    call ArrayStop('half_matelem')
    !
    !
 end subroutine Hamiltonian_vector



 subroutine calc_exp_values(nJ,Jval)

     implicit none
     integer(ik),intent(in)  :: nJ,Jval(:)

     integer(ik)  :: jind,Jrot,isym,Nentries,ientry,ilevel,alloc,rec_len,i,parmax,iunit,i1,i2,matsize,irow,ib,ktau,icontr,iterm
     real(rk)     :: f_t,dtemp0
     character(len=cl) :: unitfname,job_file

     integer(ik)        :: chkptIO,nelem,ielem,isrooti,k,itau
     character(len=cl)  :: job_is,filename,iostatus,symchar,jchar,parchar

     character(len=20)  :: buf20
     integer(hik)       :: matsize2
     integer(ik)        :: rootsize_t,imu,imu_t,dimen,irec,cdimen_,idimen,nsize,ncontr_t
     integer(hik)       :: rootsize2,rootsize
     integer(ik)        :: ilist,nlist,ideg,iroot,num
     real(rk),allocatable       :: pot_matrix(:,:)
     real(rk),allocatable       :: poten_(:,:),eigenval(:),eigenvec(:,:)
     logical            :: do_calc = .true.
     real(rk),allocatable     :: vec(:)
     type(coeffT),allocatable :: tmat(:)
     integer(ik),allocatable  :: cdimen(:)
     real(rk),allocatable     :: psi(:,:,:),mat_t(:,:)
     !double precision :: alpha,beta0,beta
     double precision,parameter :: alpha = 1.0d0,beta0=0.0d0,beta=1.0d0
     !
     if (job%verbose>=2) write(out,"(/'Prepare the matrix elements for the potential energy function')")   
     !
     !stop 'calc_exp_values is not working yet after implementing vec_sym'
     !
     call TimerStart('Prepare the potential matrix')
     !
     ! the matrix elements in the eigen-representaion  can be saved for the later use:
     ! or can be restored from the previous save
     !
     allocate(fit(sym%Nrepresen,nJ))
     !
     parmax = extF%rank
     !
     if (job%verbose>=4) then 
        write(out, '(/a, 1x, a)') 'Reading vibrational contracted matrix elements from file', trim(job%extFmat_file)
        write(out, '(a)') 'and converting them to the representaion of the eigenfunctions ... '
     end if
     !
     job_is ='external field contracted matrix elements for J=0'
     call IOStart(trim(job_is),chkptIO)
     !
     job_file = job%extFmat_file
     !
     open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=job%extFmat_file)
     !
     read(chkptIO) buf20
     if (buf20/='Start external field') then
       write (out,"(' calc_exp_values ',a,' has bogus header: ',a)") job%extFmat_file,buf20
       stop 'calc_exp_values - bogus file format'
     end if
     !
     read(chkptIO) ncontr_t
     if (trim(job%kinetmat_format)=='OLD') read(chkptIO) rootsize_t
     !
     if (bset_contr(1)%Maxcontracts/=ncontr_t) then
       write (out,"(' Dipole moment checkpoint file ',a)") job%extFmat_file
       write (out,"(' Actual and stored basis sizes at J=0 do not agree  ',i0,1x,i0)") bset_contr(1)%Maxcontracts,ncontr_t
       stop 'calc_exp_values - in file - illegal ncontracts '
     end if
     !
     rootsize = int(ncontr_t*(ncontr_t+1)/2,hik)
     rootsize2 = int(ncontr_t*ncontr_t,hik)
     !
     if (job%verbose>=4) then 
        write(out,"(/'calc_exp_values...: Number of elements: ',i8)") ncontr_t
     end if
     !
     allocate(poten_(ncontr_t,ncontr_t),stat=alloc)
     call ArrayStart('poten_',alloc,1,kind(poten_),rootsize2)
     !
     rewind(chkptIO)
     !
     ilist = 1 
     !
     do while (ilist<size(analysis%dens_list).and.analysis%dens_list(ilist)/=-1)
       !
       ilist = ilist + 1 
       !
     enddo
     nlist = max(ilist-4,0)
     !
     do jind = 1,nJ
       !
       jrot = jval(jind) 
       !
       ! Determine the size of the Hamiltonian
       !
       do isym=1,sym%Nrepresen
         !
         Nentries = 0
         !
         do ilevel = 1, Neigenlevels
             !
             irec = eigen(ilevel)%irec(1)
             !
             if (jrot/=eigen(ilevel)%jval.or.isym/=eigen(ilevel)%igamma) cycle
             !
             !jrot = eigen(ilevel)%igamma  
             !isym = eigen(ilevel)%jval 
             !
             if (eigen(ilevel)%jval==jrot) then 
                do ilist = 1,nlist 
                   if (analysis%J_list(ilist)==Jrot.and.isym==analysis%sym_list(ilist).and.analysis%dens_list(ilist)==irec)  then 
                      Nentries = Nentries + 1
                   endif
                enddo
             endif
             !
           !enddo 
           !
         enddo
         !
         fit(isym,jind)%Nentries=Nentries 
         !
         if (Nentries<1) cycle
         !
         allocate(fit(isym,jind)%ilevel(Nentries),fit(isym,jind)%ideg(Nentries),stat=alloc)
         call ArrayStart('fit%ilevel',alloc,2*size(fit(isym,jind)%ilevel),kind(fit(isym,jind)%ilevel))
         !
         Nentries = 0
         !
         do ilevel = 1, Neigenlevels
           !
           !do ideg = 1,eigen(ilevel)%ndeg
             !
             irec = eigen(ilevel)%irec(1)
             !
             if (jrot/=eigen(ilevel)%jval.or.isym/=eigen(ilevel)%igamma) cycle
             !
             if (eigen(ilevel)%jval==jrot) then 
                do ilist = 1,nlist 
                   !
                   if (analysis%J_list(ilist)==Jrot.and.isym==analysis%sym_list(ilist).and.analysis%dens_list(ilist)==irec)  then 
                     Nentries = Nentries + 1
                     fit(isym,jind)%ilevel(Nentries) = ilevel
                   endif
                enddo
             endif
             !
           !enddo 
           !
         enddo
         !
       enddo
       !
     end do
     !
     if (job%verbose>=4) write (out,"(/'Transformation to eigensolution presentaion...'/)")
     !
     do jind = 1,nJ
        !
        jrot = Jval(jind)
        !
        do isym=1,sym%Nrepresen
          !
          if (job%verbose>=2) write (out,"('jrot = ',i0,'; sym = ',i0)") Jrot,isym
          !write (out,"(/'iparam#      ilist     matrix   ')")
          !
          Nentries = fit(isym,jind)%Nentries
          !
          if (Nentries<1) cycle
          !
          matsize = Nentries*(Nentries+1)/2
          !
          allocate(pot_matrix(Nentries,Nentries),stat=alloc)
          call ArrayStart('pot_matrix',alloc,size(pot_matrix),kind(pot_matrix))
          !
          allocate(tmat(Nentries),stat=alloc)
          !
          dimen = bset_contr(jind)%Maxcontracts
          nsize = bset_contr(jind)%nsize(isym)
          !
          allocate(cdimen(Nentries),stat=alloc)
          call ArrayStart('cdimen',alloc,size(cdimen),kind(cdimen))
          !
          iunit = Jeigenvec_unit(jind,isym)
          !
          ! Prepare the transformational matrix
          !
          cdimen = 0
          !
          if (job%verbose>=5) call TimerStart('Prepare tmat for J0-convertion')
          !
          !omp parallel private(vec,alloc)  shared(cdimen,tmat)
          !
          matsize2= int(ncontr_t*Nentries,hik)
          allocate(psi(ncontr_t,Nentries,0:2*Jrot+1),mat_t(Nentries,ncontr_t),stat=alloc)
          call ArrayStart('psi',alloc,1,kind(psi),matsize2)
          call ArrayStart('mat_t',alloc,1,kind(mat_t),matsize2)
          !
          psi = 0
          !
          allocate(vec(nsize),stat = alloc)
          if (alloc /= 0) stop 'fitting-vec allocation error: vec - out of memory'
          !
          do ientry = 1,Nentries
             !
             ilevel = fit(isym,jind)%ilevel(ientry)
             !
             irec = eigen(ilevel)%irec(1)
             read(iunit, rec = irec) vec(1:nsize)
             !
             !omp parallel do private(idimen,ktau,icontr) shared(psi) schedule(guided)
             !do idimen = 1, dimen
             !   !
             !   ktau = bset_contr(jind)%ktau(idimen)
             !   icontr = bset_contr(jind)%iroot_correlat_j0(idimen)
             !   !
             !   psi(icontr,ientry,ktau) = vec(idimen)
             !   !
             !end do
             !omp end parallel do
             !
             !$omp parallel do private(idimen,irow,ib,ktau,icontr,iterm,dtemp0,nelem,ielem,isrootI) shared(psi) schedule(guided)
             do idimen = 1, dimen
               !
               irow = bset_contr(jind)%icontr2icase(idimen,1)
               ib   = bset_contr(jind)%icontr2icase(idimen,2)
               !
               ktau = bset_contr(jind)%ktau(idimen)
               icontr = bset_contr(jind)%iroot_correlat_j0(idimen)
               !
               iterm = ijterm(jind)%kmat(irow,isym)
               !
               dtemp0 = 0
               !
               nelem = bset_contr(jind)%irr(isym)%N(irow)
               !
               do ielem = 1,nelem
                  !
                  isrootI = iterm+ielem 
                  !
                  dtemp0 = dtemp0 + vec(isrootI)*bset_contr(jind)%irr(isym)%repres(isrootI,1,ib)
                  !
               enddo
               !
               psi(icontr,ientry,ktau) = dtemp0
               !
             enddo
             !$omp end parallel do
             !
          end do
          !
          deallocate(vec)
          !
          if (job%verbose>=5) call TimerStop('Prepare tmat for J0-convertion')
          !
          if (job%verbose>=5)  call MemoryReport
          !
          read(chkptIO) buf20
          read(chkptIO) ncontr_t
          !
          do i=1,parmax
            !
            if (job%verbose>=6) write (out,"('iparam = ',i0)") i
            !
            !write (out,"(' ')")
            !
            read(chkptIO) imu_t
            read(chkptIO) poten_
            !
            pot_matrix = 0
            !
            do k = 0,Jrot
              !
              do itau = 0,1
                 !
                 !if (k==0.and.mod(jrot,2)/=itau) cycle
                 !
                 ktau = 2*k+itau
                 !
                 mat_t = 0
                 !
                 call dgemm('T','N',Nentries,ncontr_t,ncontr_t,alpha,psi(:,:,ktau),ncontr_t,& 
                             poten_,ncontr_t,beta0,mat_t,Nentries)
                 call dgemm('N','N',Nentries,Nentries,ncontr_t,alpha,mat_t,Nentries,& 
                             psi(:,:,ktau),ncontr_t,beta,pot_matrix,Nentries)
                 !
                 !
              enddo
            enddo
            !
            allocate (eigenval(Nentries),eigenvec(Nentries,Nentries),stat=alloc)
            call ArrayStart('eigenval',alloc,size(eigenval),kind(eigenval))
            !
            call MLdiag_ulen(Nentries,pot_matrix,eigenval,eigenvec)
            !
            !write (out,"(' ')")
            !
            do ientry = 1, Nentries
              !
              ilevel = fit(isym,jind)%ilevel(ientry)
              !
              write (out,"(2i7,f18.6,(2x,es18.9))") i,ientry,eigen(ilevel)%energy,eigenval(ientry)
              !
            enddo
            !
            deallocate(eigenval,eigenvec)
            !
            call ArrayStop('eigenval')
            !
         enddo
         !
         if (jind==1.and.isym==1) then 
           !
           read(chkptIO) buf20(1:18)
           if (buf20(1:18)/='End external field') then
             write (out,"(' calc_exp_values ',a,' has bogus footer: ',a)") job%kinetmat_file,buf20(1:17)
             stop 'calc_exp_values - bogus file format'
           end if
           !
         endif
         !
         rewind(chkptIO)
         ! 
         deallocate(pot_matrix,stat=alloc)
         call ArrayStop('pot_matrix')
         !
         deallocate(mat_t)
         call ArrayStop('mat_t')
         !
         deallocate(psi)
         call ArrayStop('psi')
         !
         if (job%verbose>=5) call TimerReport
         !
       enddo
       !
     end do
     !
     deallocate(poten_)
     !
     close(chkptIO,status='keep')
     !
     if (job%verbose>=2) write (out,"('...done!')")
     !
     call TimerStop('Prepare the potential matrix')
     !
     if (job%verbose>=4) call TimerReport
     !
 end subroutine calc_exp_values



 subroutine prepare_diff_evib(nJ,Jval)

     integer(ik),intent(in)  :: nJ,Jval(:)

     integer(ik)  :: jind,Jrot,isym,Nentries,ientry,ilevel,alloc,rec_len,i,j0ener_fit_plus,iunit,junit,i1,i2
     real(rk)     :: f_t
     character(len=cl) :: unitfname,job_file

     character(len=cl)  :: job_is,filename,iostatus,symchar,jchar,parchar

     character(len=20)  :: buf20
     integer(ik)        :: ncontr_t,rootsize_t,imu,imu_t,dimen,nsize,irec,cdimen_,idimen,jentry,&
                           jlevel,ij,ktau_i,ktau_j,cjroot,jroot,jcontr,icase_j0,ilambda
                           !
     integer(hik)       :: hmatsize
     integer(hik)       :: rootsize,rootsize2,matsize,matsize2
     integer(ik)        :: ciroot,iroot,ktauF,icontr,irow,ib
     real(rk),allocatable       :: deriv_matrix(:,:,:)
     real(rk),allocatable       :: poten_(:,:)
     logical            :: do_calc = .true.
     real(rk),allocatable  :: vec(:)
     type(coeffT),allocatable :: tmat(:)
     integer(ik),allocatable  :: cdimen(:)
     !
     integer(ik)         :: icontr_max,ktaumax,kroot
     integer(ik),allocatable  :: kdimen(:,:,:),kdimen_(:,:)
     type(coeffT),allocatable :: kmat(:,:,:)
     logical       :: do_symmetry     
     !
     if (job%verbose>=2) write(out,"(/'Prepare the matrix elements for the potential energy function')")   
     !
     call TimerStart('Prepare derivatives wrt evib')
     !
     ! the matrix elements in the eigen-representaion  can be saved for the later use:
     ! or can be restored from the previous save
     !
     select case(trim(job%IOfitpot_action))
       !
     case('READ')
       !
       iostatus = 'old'
       do_calc = .false.
       !
     case('JOIN')
       !
       iostatus = 'replace'
       do_calc = .false.
       !
     case('SAVE','SPLIT')
       !
       iostatus = 'replace'
       !
     case default 
       !
       iostatus = 'scratch'
       do_calc = .false.
       !
     end select
     !
     ! DIVIDE can be used only with SLOW:
     !
     if (job%IOfitpot_action=='SPLIT'.and.trim(fitting%method)/='SLOW') then
       !
       write(out,"('SPLIT can be used only with SLOW, not with ',a)") fitting%method
       stop 'SPLIT can be used only with SLOW'
       !
     endif
     !
     allocate(fit(sym%Nrepresen,nJ))
     !
     !j0Neners = bset_contr(1)%Maxcontracts
     !
     j0ener_fit_plus = 0 
     !
     do i=1,j0fit%Nenergies 
        !
        do ilambda = 1, bset_contr(1)%index_deg(i)%size1
          !
          j0ener_fit_plus = j0ener_fit_plus + 1
          !
        enddo
        !
     enddo
     !
     fitting%iparam(2) = min(fitting%iparam(2),j0ener_fit_plus)
     !
     if (job%IOfitpot_action/='SPLIT') fitting%iparam = (/1,j0ener_fit_plus/)
     !
     if (job%verbose>=5) write (out,"('Total number of J=0 energies (incl. degeneracies) = ',i0/)") j0ener_fit_plus
     !
     if (job%verbose>=4) then 
        write(out, '(/a)') 'Reading the eigenvalues'
     end if
     !
     do jind = 1, nJ
       !
       Jrot = Jval(jind)
       !
       do isym= 1,sym%Nrepresen
         !
         ! Determine the size of the Hamiltonian
         !
         Nentries = 0
         !
         do ilevel = 1, Neigenlevels
           !
           if (eigen(ilevel)%jval==jrot.and.eigen(ilevel)%igamma==isym.and.job%isym_do(isym)) Nentries = Nentries + 1
           !
         enddo
         !
         if (job%verbose>=5) write (out,"('Nentries = ',i0/)") Nentries
         !
         fit(isym,jind)%Nentries=Nentries 
         !
         allocate(fit(isym,jind)%ilevel(Nentries),stat=alloc)
         call ArrayStart('fit%ilevel',alloc,size(fit(isym,jind)%ilevel),kind(fit(isym,jind)%ilevel))
         !
         Nentries = 0
         !
         do ilevel = 1, Neigenlevels
           !
           if (eigen(ilevel)%jval==jrot.and.eigen(ilevel)%igamma==isym.and.job%isym_do(isym)) then
             Nentries = Nentries + 1
             fit(isym,jind)%ilevel(Nentries) = ilevel
           endif 
           !
         enddo
         !
         if (jind==1.and.fitting%J_list(1)/=0)    cycle
         !
         if (Nentries<1)                          cycle
         !
         if (trim(job%IOfitpot_action)=='SPLIT') cycle
         !
         write(unitfname,"('deriv_matrix for j = ',i6,' sym = ',i4)") jrot,isym
         !
         call IOStart(trim(unitfname),fit(isym,jind)%IOunit)
         !
         inquire(iolength=rec_len) f_t
         !
         matsize = Nentries*(Nentries+1)/2
         !
         rec_len = rec_len*matsize
         !
         write(jchar, '(i4)') jrot
         write(symchar, '(i4)') isym
         !
         filename = trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'.chk'
         !
         !open(fit(isym,jind)%IOunit,access='direct',recl=rec_len,status=iostatus,file=filename)
         !
         open(fit(isym,jind)%IOunit,form='unformatted',position='rewind',status=iostatus,file=filename)
         !
       enddo
       !
     end do
     !
     if (do_calc) then 
       !
       if (job%verbose>=4) write (out,"(/'Transformation to the eigensolution presentaion...')")
       !
       do jind = 1, nJ
          !
          if (jind==1.and.fitting%J_list(1)/=0) cycle
          !
          Jrot = Jval(jind)
          !
          do isym= 1,sym%Nrepresen
            !
            if (job%verbose>=4) write (out,"('jrot = ',i0,'; sym = ',i0)") Jrot,isym
            !
            Nentries = fit(isym,jind)%Nentries
            !
            if (Nentries<1) cycle
            !
            matsize = int(Nentries*(Nentries+1)/2,hik)
            matsize2= int(Nentries*(Nentries+1)/2,hik)
            !
            if (job%verbose>=5) write (out,"('matsize,j0ener_fit_plus = ',2i0)") matsize,j0ener_fit_plus
            !
            idimen = j0ener_fit_plus
            !
            if (trim(fitting%method)=='SLOW') idimen = 1
            !
            hmatsize = int(matsize*idimen,hik)
            !
            allocate(deriv_matrix(Nentries,Nentries,idimen),stat=alloc)
            call ArrayStart('deriv_matrix',alloc,1,kind(deriv_matrix),hmatsize)
            !
            deriv_matrix = 0
            !
            allocate(tmat(Nentries),stat=alloc)
            !
            dimen = bset_contr(jind)%Maxcontracts
            nsize = bset_contr(jind)%nsize(isym)
            !
            if (trim(fitting%method)=='SLOW'.or.trim(fitting%method)=='FAST') then 
              !
              if (job%verbose>=5) write(out,"('not very fast approach...')")
              !
              !
              allocate(cdimen(Nentries),stat=alloc)
              call ArrayStart('cdimen',alloc,size(cdimen),kind(cdimen))
              !
              iunit = Jeigenvec_unit(jind,isym)
              !
              ! Prepare the transformational matrix
              !
              cdimen = 0
              !
              if (job%verbose>=5) call TimerStart('Prepare tmat for J0-convertion')
              !
              !omp parallel private(vec,alloc)  ! shared(cdimen,tmat)
              allocate(vec(nsize),stat = alloc)
              if (alloc /= 0) stop 'fitting-vec allocation error: vec - out of memory'
              !
              !omp do private(ientry,ilevel,irec,cdimen_,idimen) schedule(guided)
              do ientry = 1,Nentries
                 !
                 ilevel = fit(isym,jind)%ilevel(ientry)
                 !
                 irec = eigen(ilevel)%irec(1)
                 read(iunit, rec = irec) vec(1:nsize)
                 !
                 ! Determine number of nonzero terms 
                 !
                 cdimen_ = 0 
                 do idimen = 1, dimen
                    if (abs(vec(idimen)) > fitting%threshold_coeff) then
                       cdimen_ = cdimen_ + 1
                    end if
                 end do
                 !
                 if (job%verbose>=6) then 
                   write (out,"(' ientry  =  ',i0,' cdimen_ = ',i0)") ientry,cdimen_
                 end if
                 !
                 cdimen(ientry) = cdimen_
                 !
                 allocate(tmat(ientry)%coeff(cdimen_),stat=alloc)
                 if (alloc /= 0) stop 'fitting-vec allocation error: tmat%coeff - out of memory'
                 !
                 allocate(tmat(ientry)%icoeff(cdimen_),stat=alloc)
                 if (alloc /= 0) stop 'fitting-vec allocation error: tmat%icoeff - out of memory'
                 !
                 cdimen_ = 0
                 do idimen = 1, dimen
                    if (abs(vec(idimen)) > fitting%threshold_coeff) then
                       !
                       cdimen_ = cdimen_ + 1
                       tmat(ientry)%icoeff(cdimen_) = idimen
                       tmat(ientry)%coeff(cdimen_)  = vec(idimen)
                       !
                    end if
                 end do
                 !
              end do
              !omp end do
              !
              deallocate(vec)
              !omp end parallel
              !
              do ientry = 1,Nentries
                 call ArrayStart('tmat',0,size(tmat(ientry)%coeff),kind(tmat(ientry)%coeff))
                 call ArrayStart('tmat',0,size(tmat(ientry)%icoeff),kind(tmat(ientry)%icoeff))
              end do

              ! 
            elseif (trim(fitting%method)=='VERYFAST') then 
              !
              if (job%verbose>=5) write(out,"('very fast (and very expensive) approach...')")     
              !
              icontr_max = maxval(bset_contr(jind)%iroot_correlat_j0(:),dim=1)
              ktaumax = 2*jrot+1
              !
              if (job%verbose>=5) write(out,"('ktaumax = ',i0,'; icontr_max = ',i0,' dimen = ',i0)") ktaumax,icontr_max,dimen
              !
              allocate(kmat(Nentries,icontr_max,0:ktaumax),stat=alloc)
              !
              allocate(kdimen_(icontr_max,0:ktaumax),stat=alloc)
              call ArrayStart('kdimen_',alloc,size(kdimen_),kind(kdimen_))
              !
              allocate(kdimen(Nentries,icontr_max,0:ktaumax),stat=alloc)
              call ArrayStart('kdimen',alloc,size(kdimen),kind(kdimen))
              !
              allocate(cdimen(Nentries),stat=alloc)
              call ArrayStart('cdimen',alloc,size(cdimen),kind(cdimen))
              !
              iunit = Jeigenvec_unit(jind,isym)
              !
              ! Prepare the transformational matrix
              !
              kdimen = 0
              !
              if (job%verbose>=5) call TimerStart('Prepare tmat for J0-convertion')
              !
              !omp parallel private(vec,alloc)  ! shared(cdimen,tmat)
              allocate(vec(dimen),stat = alloc)
              if (alloc /= 0) stop 'fitting-vec allocation error: vec - out of memory'
              !
              !omp do private(ientry,ilevel,irec,cdimen_,idimen) schedule(guided)
              do ientry = 1,Nentries
                 !
                 if (job%verbose>=5.and.mod(ientry,Nentries/50)==0) write(out,"('ientry = ',i0)") ientry
                 !
                 ilevel = fit(isym,jind)%ilevel(ientry)
                 !
                 irec = eigen(ilevel)%irec(1)
                 read(iunit, rec = irec) vec(1:dimen)
                 !
                 ! Determine number of nonzero terms 
                 !
                 kdimen_ = 0 
                 cdimen_ = 0 
                 !
                 do idimen = 1, dimen
                    if (abs(vec(idimen)) > fitting%threshold_coeff) then
                       !
                       cdimen_ = cdimen_ + 1
                       !
                       icontr = bset_contr(jind)%iroot_correlat_j0(idimen)
                       !
                       ktau_i = bset_contr(jind)%ktau(idimen)
                       !
                       kdimen_(icontr,ktau_i) = kdimen_(icontr,ktau_i) + 1
                       !
                    end if
                 end do
                 !
                 allocate(tmat(ientry)%coeff(cdimen_),stat=alloc)
                 if (alloc /= 0) stop 'fitting-vec allocation error: tmat%coeff - out of memory'
                 !
                 allocate(tmat(ientry)%icoeff(cdimen_),stat=alloc)
                 if (alloc /= 0) stop 'fitting-vec allocation error: tmat%icoeff - out of memory'
                 !
                 call ArrayStart('tmat',0,size(tmat(ientry)%coeff),kind(tmat(ientry)%coeff))
                 call ArrayStart('tmat',0,size(tmat(ientry)%icoeff),kind(tmat(ientry)%icoeff))
                 !
                 cdimen(ientry) = cdimen_
                 kdimen(ientry,:,:) = kdimen_(:,:)
                 !
                 do ktau_i = 0,ktaumax
                   !
                   do icontr = 1,icontr_max
                     !
                     cdimen_ = max(kdimen(ientry,icontr,ktau_i),1)
                     !
                     if (job%verbose>=6) write (out,"(' ientry  =  ',i0,'ktau_i = ',i4,' icontr = ',i0,' cdimen_ = ',i0)") &
                                         ientry,ktau_i,icontr,cdimen_
                     !
                     !allocate(kmat(ientry,icontr,ktau_i)%coeff(cdimen_),stat=alloc)
                     !if (alloc /= 0) stop 'fitting-vec allocation error: kmat%coeff - out of memory'
                     !
                     allocate(kmat(ientry,icontr,ktau_i)%icoeff(cdimen_),stat=alloc)
                     if (alloc /= 0) stop 'fitting-vec allocation error: kmat%icoeff - out of memory'
                     !
                     !call ArrayStart('kmat',0,size(kmat(ientry,icontr,ktau_i)%coeff),kind(kmat(ientry,icontr,ktau_i)%coeff))
                     call ArrayStart('kmat',0,size(kmat(ientry,icontr,ktau_i)%icoeff),kind(kmat(ientry,icontr,ktau_i)%icoeff))
                     !
                   enddo 
                   !
                 enddo 
                 !
                 kdimen_ = 0
                 cdimen_ = 0 
                 !
                 do idimen = 1, dimen
                    if (abs(vec(idimen)) > fitting%threshold_coeff) then
                       !
                       icontr = bset_contr(jind)%iroot_correlat_j0(idimen)
                       !
                       ktau_i = bset_contr(jind)%ktau(idimen)
                       !
                       kdimen_(icontr,ktau_i) = kdimen_(icontr,ktau_i) + 1
                       !
                       cdimen_ = cdimen_ + 1
                       !
                       tmat(ientry)%icoeff(cdimen_)  = idimen
                       tmat(ientry)%coeff(cdimen_)   = vec(idimen)
                       !
                       kmat(ientry,icontr,ktau_i)%icoeff(kdimen_(icontr,ktau_i)) = cdimen_
                       !
                    end if
                 end do
                 !
              end do
              !omp end do
              !
              deallocate(vec)
              !omp end parallel
              !
              deallocate(kdimen_)
              !
            else
              !
              write(out,"('unknown j0fit-method: ',a)") fitting%method
              stop 'unknown j0fit-method'
              !
            endif
            !
            if (job%verbose>=5) call TimerStop('Prepare tmat for J0-convertion')
            if (job%verbose>=5)  call MemoryReport
            !
            iunit = fit(isym,jind)%IOunit
            !
            if (job%IOfitpot_action/='SPLIT') write (iunit) 'start fitener'
            !
            if (trim(fitting%method)=='SLOW') then 
              !
              if (job%verbose>=5) write(out,"('slow and less expensive approach...')")
              !
              do i=fitting%iparam(1),fitting%iparam(2)
                !
                if (job%verbose>=5) write (out,"('iparam = ',i0)") i
                !
                do ientry = 1, Nentries
                  !
                  if (job%verbose>=5.and.mod(ientry,Nentries/50)==0) write(out,"('ientry = ',i0)") ientry
                  !
                  ilevel = fit(isym,jind)%ilevel(ientry)
                  !
                  do jentry = 1,ientry
                     !
                     jlevel = fit(isym,jind)%ilevel(jentry)
                     !
                     !ij = ientry*(ientry-1)/2+jentry
                     !
                     !$omp  parallel do private(ciroot,iroot,icontr,ktau_i,cjroot,jroot,jcontr,ktau_j) &
                     !$omp& shared(deriv_matrix) schedule(guided) 
                     do ciroot = 1, cdimen(ientry)
                       !
                       iroot = tmat(ientry)%icoeff(ciroot)
                       !
                       icontr = bset_contr(jind)%iroot_correlat_j0(iroot)
                       !
                       if (icontr/=i) cycle 
                       !
                       ktau_i = bset_contr(jind)%ktau(iroot)
                       !
                       do cjroot = 1, cdimen(jentry)
                         !
                         jroot = tmat(jentry)%icoeff(cjroot)
                         !
                         jcontr = bset_contr(jind)%iroot_correlat_j0(jroot)
                         !
                         ktau_j = bset_contr(jind)%ktau(jroot)
                         !
                         if (icontr>j0ener_fit_plus.or.icontr/=jcontr.or.ktau_i/=ktau_j) cycle
                         !
                         deriv_matrix(ientry,jentry,1) = deriv_matrix(ientry,jentry,1) + &
                                              tmat(ientry)%coeff(ciroot)*tmat(jentry)%coeff(cjroot)
                         !
                       enddo
                       !
                     enddo
                     !$omp end parallel do
                     !
                  enddo
                  !
                enddo
                !
                if (job%IOfitpot_action=='SPLIT') then 
                  !
                  write(unitfname,"('single pot_matrix')")
                  !
                  call IOStart(trim(unitfname),junit)
                  !
                  write(jchar, '(i4)') jrot
                  write(symchar, '(i4)') isym
                  write(parchar, '(i4)') i
                  !
                  filename = &
                     trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'_'//trim(adjustl(parchar))//'.tmp'
                  !
                  open(junit,form='unformatted',action='write',position='rewind',status='replace',file=filename)
                  !
                  write (junit) deriv_matrix(:,:,1)
                  !
                  close(junit)
                  !
                else
                  !
                  write (iunit) i
                  write (iunit) deriv_matrix(:,:,1)
                  !
                endif
                !
              enddo
              !
              if (job%IOfitpot_action/='SPLIT') then 
                 write (iunit) 'end fitener'
                 close(iunit)
              endif 
              !
              deallocate(tmat)
              call ArrayStop('tmat')
              !
              deallocate(cdimen)
              call ArrayStop('cdimen')
              ! 
            elseif (trim(fitting%method)=='FAST') then 
              !
              if (job%verbose>=4) write(out,"('fast approach...')")     
              !
              do ientry = 1, Nentries
                !
                if (job%verbose>=5.and.mod(ientry,Nentries/50)==0) write(out,"('ientry = ',i0)") ientry
                !
                ilevel = fit(isym,jind)%ilevel(ientry)
                !
                do jentry = 1,ientry
                   !
                   jlevel = fit(isym,jind)%ilevel(jentry)
                   !
                   !ij = ientry*(ientry-1)/2+jentry
                   !
                   !$omp parallel do private(ciroot,iroot,icontr,icase_j0,ktau_i,cjroot,jroot,jcontr,ktau_j) &
                   !$omp& shared(deriv_matrix) schedule(guided) 
                   l_iroot: do ciroot = 1, cdimen(ientry)
                     !
                     iroot = tmat(ientry)%icoeff(ciroot)
                     !
                     icontr = bset_contr(jind)%iroot_correlat_j0(iroot)
                     !
                     icase_j0  = bset_contr(1)%icontr2icase(icontr, 1)
                     !
                     if (icontr>j0ener_fit_plus.or.icontr>bset_contr(1)%Maxcontracts.or.icase_j0>j0fit%Nenergies) cycle
                     !
                     ktau_i = bset_contr(jind)%ktau(iroot)
                     !
                     l_jroot : do cjroot = 1, cdimen(jentry)
                       !
                       jroot = tmat(jentry)%icoeff(cjroot)
                       !
                       jcontr = bset_contr(jind)%iroot_correlat_j0(jroot)
                       !
                       ktau_j = bset_contr(jind)%ktau(jroot)
                       !
                       if (icontr/=jcontr.or.ktau_i/=ktau_j) cycle l_jroot
                       !
                       deriv_matrix(ientry,jentry,icontr) = deriv_matrix(ientry,jentry,icontr) + &
                                                 tmat(ientry)%coeff(ciroot)*tmat(jentry)%coeff(cjroot)
                       !
                       cycle l_iroot
                       !
                     enddo l_jroot
                     !
                   enddo l_iroot
                   !$omp end parallel do
                   !
                enddo
                !
              enddo
		      !
              iunit = fit(isym,jind)%IOunit
              !write (iunit) 'start fitener'
              !
              do i=1,j0ener_fit_plus
                !
                write (iunit) i
                write (iunit) deriv_matrix(:,:,i)
                !
                !write (iunit,rec=i) deriv_matrix(1:matsize,i)
                !
              enddo
              !
              write (iunit) 'end fitener'
              close(iunit)
              !
              deallocate(tmat)
              call ArrayStop('tmat')
              !
              deallocate(cdimen)
              call ArrayStop('cdimen')
              !
            elseif (trim(fitting%method)=='VERYFAST') then 
              !
              if (job%verbose>=5) write(out,"('fast and very expensive approach...')")     
              !
              do ientry = 1, Nentries
                !
                if (job%verbose>=5.and.mod(ientry,Nentries/50)==0) write(out,"('ientry = ',i0)") ientry
                !
                ilevel = fit(isym,jind)%ilevel(ientry)
                !
                do jentry = 1,ientry
                   !
                   jlevel = fit(isym,jind)%ilevel(jentry)
                   !
                   !ij = ientry*(ientry-1)/2+jentry
                   !
                   !$omp parallel do private(ciroot,iroot,icontr,icase_j0,ktau_i,kroot,cjroot) shared(deriv_matrix) schedule(guided) 
                   do ciroot = 1, cdimen(ientry)
                     !
                     iroot = tmat(ientry)%icoeff(ciroot)
                     !
                     icontr = bset_contr(jind)%iroot_correlat_j0(iroot)
                     !
                     icase_j0  = bset_contr(1)%icontr2icase(icontr, 1)
                     !
                     if (icontr>j0ener_fit_plus.or.icontr>bset_contr(1)%Maxcontracts.or.icase_j0>j0fit%Nenergies) cycle
                     !
                     ktau_i = bset_contr(jind)%ktau(iroot)
                     !
                     do kroot = 1, kdimen(jentry,icontr,ktau_i)
                       !
                       cjroot = kmat(jentry,icontr,ktau_i)%icoeff(kroot)
                       !
                       deriv_matrix(ientry,jentry,icontr) = deriv_matrix(ientry,jentry,icontr) + &
                                                            tmat(ientry)%coeff(ciroot)*tmat(jentry)%coeff(cjroot)
                       !
                     enddo
                     !
                   enddo
                   !$omp end parallel do
                   !
                enddo
                !
              enddo
		      !
              iunit = fit(isym,jind)%IOunit
              !write (iunit) 'start fitener'
              !
              do i=1,j0ener_fit_plus
                !
                write (iunit) i
                write (iunit) deriv_matrix(:,:,i)
                !
                !write (iunit,rec=i) deriv_matrix(1:matsize,i)
                !
              enddo
              !
              write (iunit) 'end fitener'
              close(iunit)
              !
              deallocate(tmat)
              call ArrayStop('tmat')
              !
              deallocate(kmat)
              call ArrayStop('kmat')
              !
              deallocate(kdimen)
              call ArrayStop('kdimen')
              !
              deallocate(cdimen)
              call ArrayStop('cdimen')
              !
            else
              !
              write(out,"('unknown j0fit-method: ',a)") fitting%method
              stop 'unknown j0fit-method'
              !
            endif
            !
            deallocate(deriv_matrix)
            call ArrayStop('deriv_matrix')
            !
            if (job%verbose>=5) call TimerReport
            !
            if (job%IOfitpot_action/='SPLIT') then 
              !
              write(jchar, '(i4)') jrot
              write(symchar, '(i4)') isym
              filename = trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'.chk'
              open(fit(isym,jind)%IOunit,form='unformatted',position='rewind',status='old',action='read',file=filename)
            endif 
            !
          enddo
          !
       end do
       !
     else if (job%IOfitpot_action=='JOIN') then 
       !
       if (job%verbose>=4) write (out,"(/'Join computed mat. elements ...')")
       !
       do jind = 1, nJ
          !
          if (jind==1.and.fitting%J_list(1)/=0) cycle
          !
          Jrot = Jval(jind)
          !
          do isym= 1,sym%Nrepresen
            !
            if (job%verbose>=4) write (out,"('jrot = ',i0,'; sym = ',i0)") Jrot,isym
            !
            Nentries = fit(isym,jind)%Nentries
            !
            if (Nentries<1) cycle
            !
            matsize = int(Nentries*(Nentries+1)/2,hik)
            matsize2= int(Nentries*Nentries,hik)
            !
            if (job%verbose>=5) write (out,"('matsize,j0ener_fit_plus = ',2i0)") matsize,j0ener_fit_plus
            !
            if (trim(fitting%method)=='SLOW') idimen = 1
            !
            hmatsize = matsize2
            !
            allocate(deriv_matrix(Nentries,Nentries,1),stat=alloc)
            call ArrayStart('deriv_matrix',alloc,1,kind(deriv_matrix),matsize2)
            !
            iunit = fit(isym,jind)%IOunit
            write (iunit) 'start fitener'
            !
            do i=1,j0ener_fit_plus
              !
              if (job%verbose>=5) write (out,"('iparam = ',i0)") i
              !
              write(unitfname,"('single deriv_matrix')")
              !
              call IOStart(trim(unitfname),junit)
              !
              write(jchar, '(i4)') jrot
              write(symchar, '(i4)') isym
              write(parchar, '(i4)') i
              !
              filename = &
                  trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'_'//trim(adjustl(parchar))//'.chk'
              !
              open(junit,form='unformatted',action='read',position='rewind',status='old',file=filename)
              !
              read (junit) deriv_matrix(:,:,1)
              !
              close(junit)
              !
              write(iunit) i
              write (iunit) deriv_matrix(:,:,1)
              !
           enddo
           !
           write (iunit) 'end fitener'
           close(iunit)
           ! 
           deallocate(deriv_matrix,stat=alloc)
           call ArrayStop('deriv_matrix')
           !
           if (job%verbose>=5) call TimerReport
           !
         enddo
         !
       end do
       !
    else if (job%IOfitpot_action=='READ') then
       ! 
       do jind = 1, nJ
          !
          if (jind==1.and.fitting%J_list(1)/=0) cycle
          !
          Jrot = Jval(jind)
          !
          do isym= 1,sym%Nrepresen
            !
            Nentries = fit(isym,jind)%Nentries
            !
            if (Nentries<1) cycle
            iunit = fit(isym,jind)%IOunit
            !
            do_symmetry = .false.
            do i = 1,fitting%Nenergies
              if (fitting%obs(i)%Jrot==Jrot.and.fitting%obs(i)%symmetry==isym) then 
                do_symmetry = .true.
                exit
              endif
            enddo
            !
            if (.not.do_symmetry) cycle
            !
            read(iunit) buf20(1:13)
            !
            if (buf20(1:13)/='start fitener') then
              write (out,"(' evib_matrix at J = ',i4,', isym = ',i4,' has bogus footer: ',a)") jrot,isym,buf20(1:13)
              stop 'evib - bogus file format'
            end if
            !
         enddo
         !
       end do
       !
     endif
     !
     if (job%verbose>=2) write (out,"('...done!')")
     !
     call TimerStop('Prepare derivatives wrt evib')
     !
     if (job%verbose>=4) call TimerReport
     !
 end subroutine prepare_diff_evib
 !
 subroutine prepare_pot_matrix(nJ,Jval)
 !
     implicit none
     integer(ik),intent(in)  :: nJ,Jval(:)
     !
     integer(ik)       :: jind,Jrot,isym,Nentries,ientry,ilevel,alloc,rec_len,i,parmax,iunit,junit,i1,i2
     real(rk)          :: f_t,dtemp0
     character(len=cl) :: unitfname,job_file
     !
     integer(ik)        :: chkptIO
     character(len=cl)  :: job_is,filename,iostatus,symchar,jchar,parchar,pot_suffix
     !
     character(len=20)  :: buf20
     integer(ik)        :: ncontr_t,iterm,nelem,ielem,isrootI,irow,ib
     integer(ik)        :: rootsize_t,imu,imu_t,dimen,nsize,irec,cdimen_,idimen,j,iterm1,iterm2,ktau,icontr,k,itau
     integer(hik)       :: rootsize,rootsize2,matsize,matsize2
     real(rk),allocatable  :: pot_matrix(:,:)
     real(rk),allocatable  :: poten_(:,:)
     logical               :: do_calc = .true.
     real(rk),allocatable  :: vec(:)
     type(coeffT),allocatable :: tmat(:)
     integer(ik),allocatable  :: cdimen(:)
     real(rk),allocatable     :: psi(:,:,:),mat_t(:,:)
     double precision,parameter :: alpha = 1.0d0,beta0=0.0d0,beta=1.0d0
     !
     if (job%verbose>=2) write(out,"(/'Prepare the matrix elements for the potential energy function')")   
     !
     call TimerStart('Prepare the potential matrix')
     !
     ! the matrix elements in the eigen-representaion  can be saved for the later use:
     ! or can be restored from the previous save
     !
     select case(trim(job%IOfitpot_action))
       !
     case('READ')
       !
       iostatus = 'old'
       do_calc = .false.
       !
     case('REWRITE')
       !
       iostatus = 'old'
       do_calc = .false.
       !
     case('JOIN')
       !
       iostatus = 'replace'
       do_calc = .false.
       !
     case('SAVE','SPLIT')
       !
       iostatus = 'replace'
       !
     case default 
       !
       iostatus = 'scratch'
       !
     end select
     !
     allocate(fit(sym%Nrepresen,nJ))
     !
     parmax = extF%rank
     !
     fitting%iparam(2) = min(fitting%iparam(2),parmax)
     !
     if (.not.job%IOfitpot_divide) fitting%iparam = (/1,parmax/)
     !
     if (do_calc) then 
       !
       if (job%verbose>=4) then 
          write(out, '(/a, 1x, a)') 'Reading vibrational contracted matrix elements from file', trim(job%extFmat_file)
          write(out, '(a)') 'and converting them to the representaion of the eigenfunctions ... '
       end if
       !
       if (job%IOextF_divide) then
         !
         ncontr_t = bset_contr(1)%Maxcontracts
         !
       else
         !
         job_is ='external field contracted matrix elements for J=0'
         call IOStart(trim(job_is),chkptIO)
         !
         job_file = job%extFmat_file
         !
         open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=job_file)
         !
         read(chkptIO) buf20
         if (buf20/='Start external field') then
           write (out,"(' prepare_pot_matrix ',a,' has bogus header: ',a)") job_file,buf20
           stop 'prepare_pot_matrix - bogus file format'
         end if
         !
         read(chkptIO) ncontr_t
         !
         if (bset_contr(1)%Maxcontracts/=ncontr_t) then
           write (out,"(' Dipole moment checkpoint file ',a)") job_file
           write (out,"(' Actual and stored basis sizes at J=0 do not agree  ',2i0)") bset_contr(1)%Maxcontracts,ncontr_t
           stop 'prepare_pot_matrix - in file - illegal ncontracts '
         end if
         !
         rewind(chkptIO)
         !
       endif
       !
       if (job%IOfitpot_divide) then
         !
         iterm1 = max(fitting%iparam(1),1)
         iterm2 = min(fitting%iparam(2),parmax)
         !
         if (job%verbose>=4) write(out,"('  The fitpot-chks will be divided into ',i6,' chk-slices')")extF%rank
         if (job%verbose>=4) write(out,"('  This run is for the checkpoint slices from ',i4,' to ',i4)") iterm1,iterm2
         !
       endif
       !
       rootsize  = int(ncontr_t*(ncontr_t+1)/2,hik)
       rootsize2 = int(ncontr_t*ncontr_t,hik)
       !
       if (job%verbose>=4) then 
          write(out,"(/'prepare_pot_matrix...: Number of elements: ',i8)") ncontr_t
       end if
       !
       allocate(poten_(ncontr_t,ncontr_t),stat=alloc)
       call ArrayStart('poten_',alloc,1,kind(poten_),rootsize2)
       !
       !if (job%IOfitpot_action=='SPLIT') close(chkptIO)
       !
     endif
     !
     do jind = 1, nJ
       !
       Jrot = Jval(jind)
       !
       do isym= 1,sym%Nrepresen
         !
         ! Determine the size of the Hamiltonian
         !
         Nentries = 0
         !
         do ilevel = 1, Neigenlevels
           !
           if (eigen(ilevel)%jval==jrot.and.eigen(ilevel)%igamma==isym.and.job%isym_do(isym)) Nentries = Nentries + 1
           !
         enddo
         !
         fit(isym,jind)%Nentries=Nentries 
         !
         allocate(fit(isym,jind)%ilevel(Nentries),stat=alloc)
         call ArrayStart('fit%ilevel',alloc,size(fit(isym,jind)%ilevel),kind(fit(isym,jind)%ilevel))
         !
         Nentries = 0
         !
         do ilevel = 1, Neigenlevels
           !
           if (eigen(ilevel)%jval==jrot.and.eigen(ilevel)%igamma==isym.and.job%isym_do(isym)) then
             Nentries = Nentries + 1
             fit(isym,jind)%ilevel(Nentries) = ilevel
           endif 
           !
         enddo
         !
         fit(isym,jind)%IOunit = 0
         !
         if (jind==1.and.fitting%J_list(1)/=0)    cycle
         !
         if (Nentries<1)                          cycle
         !
         if (job%IOfitpot_divide) cycle
         !
         !matsize = Nentries*(Nentries+1)/2
         !
         write(unitfname,"('pot_matrix for j = ',i6,' sym = ',i4)") jrot,isym
         !
         call IOStart(trim(unitfname),fit(isym,jind)%IOunit)
         !
         !inquire(iolength=rec_len) f_t
         !
         !rec_len = rec_len*matsize
         !
         !fitting%JGamma_do(jrot,isym) = .true.
         !
         write(jchar, '(i4)') jrot
         write(symchar, '(i4)') isym
         !
         pot_suffix = trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))
         !
         filename = trim(pot_suffix)//'.chk'
         !
         open(fit(isym,jind)%IOunit,form='unformatted',position='rewind',status=iostatus,file=filename,err=22)
         !
         cycle
         !
         22 fit(isym,jind)%Nentries = 0 ! fitting%JGamma_do(jrot,isym) = .false.
         !
         !open(fit(isym,jind)%IOunit,access='direct',recl=rec_len,status=iostatus,file=filename)
         !
       enddo
       !
     end do
     !
     if (do_calc) then 
       !
       if (job%verbose>=4) write (out,"(/'Transformation to eigensolution presentaion...')")
          !
          if (job%verbose>=4) then 
             select case (trim(fitting%method))
             case default
               write (out,"(/'Transformation to eigensolution presentaion (fast)...')")
             case ('SLOW')
               write (out,"(/'Transformation to eigensolution presentaion (slow)...')")
             end select
          end if
          !
       do jind = 1, nJ
          !
          if (jind==1.and.fitting%J_list(1)/=0) cycle
          !
          Jrot = Jval(jind)
          !
          do isym= 1,sym%Nrepresen
            !
            if (job%verbose>=4) write (out,"('jrot = ',i0,'; sym = ',i0)") Jrot,isym
            !
            Nentries = fit(isym,jind)%Nentries
            !
            if (Nentries<1) cycle
            !
            matsize = int(Nentries*(Nentries+1)/2,hik)
            matsize2= int(Nentries*Nentries,hik)
            !
            allocate(pot_matrix(Nentries,Nentries),stat=alloc)
            call ArrayStart('pot_matrix',alloc,1,kind(pot_matrix),matsize2)
            !
            allocate(tmat(Nentries),stat=alloc)
            !
            dimen = bset_contr(jind)%Maxcontracts
            nsize = bset_contr(jind)%nsize(isym)
            !
            !iunit = TReigenvec_unit(jind,Jval)
            !
            iunit = Jeigenvec_unit(jind,isym)
            !
            if (trim(fitting%method)=='FAST') then
              !
              matsize2= int(ncontr_t*Nentries,hik)
              allocate(psi(ncontr_t,Nentries,0:2*Jrot+1),mat_t(Nentries,ncontr_t),stat=alloc)
              call ArrayStart('psi',alloc,1,kind(psi),matsize2)
              call ArrayStart('mat_t',alloc,1,kind(mat_t),matsize2)
              !
              psi = 0
              !
              allocate(vec(nsize),stat = alloc)
              if (alloc /= 0) stop 'fitting-vec allocation error: vec - out of memory'
              !
              do ientry = 1,Nentries
                 !
                 ilevel = fit(isym,jind)%ilevel(ientry)
                 !
                 irec = eigen(ilevel)%irec(1)
                 read(iunit, rec = irec) vec(1:nsize)
                 !
                 !$omp parallel do private(idimen,irow,ib,ktau,icontr,iterm,dtemp0,nelem,ielem,isrootI) shared(psi) schedule(guided)
                 do idimen = 1, dimen
                   !
                   irow = bset_contr(jind)%icontr2icase(idimen,1)
                   ib   = bset_contr(jind)%icontr2icase(idimen,2)
                   !
                   ktau = bset_contr(jind)%ktau(idimen)
                   icontr = bset_contr(jind)%iroot_correlat_j0(idimen)
                   !
                   iterm = ijterm(jind)%kmat(irow,isym)
                   !
                   dtemp0 = 0
                   !
                   nelem = bset_contr(jind)%irr(isym)%N(irow)
                   !
                   do ielem = 1,nelem
                      !
                      isrootI = iterm+ielem 
                      !
                      dtemp0 = dtemp0 + vec(isrootI)*bset_contr(jind)%irr(isym)%repres(isrootI,1,ib)
                      !
                   enddo
                   !
                   psi(icontr,ientry,ktau) = dtemp0
                   !
                 enddo
                 !$omp end parallel do
                 !
              end do
              !
              deallocate(vec)
              !
            else
              !
              write(out,"('SLOW option in prepare_pot_matrix is not working after adding sym_vec')")
              !
              stop 'SLOW option in prepare_pot_matrix is not working'
              !
              ! Prepare the transformational matrix
              !
              allocate(cdimen(Nentries),stat=alloc)
              call ArrayStart('cdimen',alloc,size(cdimen),kind(cdimen))
              !
              cdimen = 0
              !
              if (job%verbose>=5) call TimerStart('Prepare tmat for J0-convertion')
              !
              !omp parallel private(vec,alloc)  ! shared(cdimen,tmat)
              allocate(vec(nsize),stat = alloc)
              if (alloc /= 0) stop 'fitting-vec allocation error: vec - out of memory'
              !
              !omp do private(ientry,ilevel,irec,cdimen_,idimen) schedule(guided)
              do ientry = 1,Nentries
                 !
                 ilevel = fit(isym,jind)%ilevel(ientry)
                 !
                 irec = eigen(ilevel)%irec(1)
                 read(iunit, rec = irec) vec(1:nsize)
                 !
                 ! Determine number of nonzero terms 
                 !
                 cdimen_ = 0 
                 do idimen = 1, dimen
                    if (abs(vec(idimen)) > fitting%threshold_coeff) then
                       cdimen_ = cdimen_ + 1
                    end if
                 end do
                 !
                 if (job%verbose>=6) then 
                   write (out,"(' ientry  =  ',i0,' cdimen_ = ',i0)") ientry,cdimen_
                 end if
                 !
                 cdimen(ientry) = cdimen_
                 !
                 allocate(tmat(ientry)%coeff(cdimen_),stat=alloc)
                 if (alloc /= 0) stop 'fitting-vec allocation error: tmat%coeff - out of memory'
                 !
                 allocate(tmat(ientry)%icoeff(cdimen_),stat=alloc)
                 if (alloc /= 0) stop 'fitting-vec allocation error: tmat%icoeff - out of memory'
                 !
                 cdimen_ = 0
                 do idimen = 1, dimen
                    if (abs(vec(idimen)) > fitting%threshold_coeff) then
                       !
                       cdimen_ = cdimen_ + 1
                       tmat(ientry)%icoeff(cdimen_) = idimen
                       tmat(ientry)%coeff(cdimen_) = vec(idimen)
                       !
                    end if
                 end do
                 !
                 !tmat(ientry)%coeff(:) = vec(tmat(ientry)%icoeff(:))
                 !
              end do
              !omp end do
              !
              deallocate(vec)
              !omp end parallel
              !
              do ientry = 1,Nentries
                 call ArrayStart('tmat',0,size(tmat(ientry)%coeff),kind(tmat(ientry)%coeff))
                 call ArrayStart('tmat',0,size(tmat(ientry)%icoeff),kind(tmat(ientry)%icoeff))
              end do
              !
              if (job%verbose>=5) call TimerStop('Prepare tmat for J0-convertion')
              !
            endif
            !
            if (job%verbose>=6)  call MemoryReport
            !
            if (.not.job%IOextF_divide) then
              !
              read(chkptIO) buf20
              read(chkptIO) ncontr_t
              !
            endif
            !
            if (.not.job%IOfitpot_divide) then 
              iunit = fit(isym,jind)%IOunit
              write(iunit) 'start potfit'
            endif
            !
            do i=1,fitting%iparam(2)
              !
              if (.not.job%IOextF_divide) then 
                !
                read(chkptIO) imu_t
                !
                !if (imu_t>fitting%iparam(2)) exit
                !
                read(chkptIO) poten_
                !
              endif
              !
              !if (i<fitting%iparam(1).or.extF%ifit(i)==0) cycle
              !
              if (i<fitting%iparam(1)) cycle
              !
              if (job%IOextF_divide) then
                !
                if (extF%ifit(1,i)<1) cycle
                !
                call divided_slice_read(i,'extF',job%extmat_suffix,ncontr_t,poten_)
                !
              endif
              !
              if (job%verbose>=5) write (out,"('  iparam = 'i4)") i
              !
              if (trim(fitting%method)=='FAST') then
                !
                pot_matrix = 0
                !
                do k = 0,Jrot
                  !
                  do itau = 0,1
                     !
                     !if (k==0.and.mod(jrot,2)/=itau) cycle
                     !
                     ktau = 2*k+itau
                     !
                     call dgemm('T','N',Nentries,ncontr_t,ncontr_t,alpha,psi(:,:,ktau),ncontr_t,& 
                                 poten_,ncontr_t,beta0,mat_t,Nentries)
                     call dgemm('N','N',Nentries,Nentries,ncontr_t,alpha,mat_t,Nentries,& 
                                 psi(:,:,ktau),ncontr_t,beta,pot_matrix,Nentries)
                     !
                     !
                  enddo
                enddo
                !
              else 
                !
                do ientry = 1, Nentries
                  !
                  ilevel = fit(isym,jind)%ilevel(ientry)
                  !
                  call Hamiltonian_vector(ientry,jind,Nentries,tmat,cdimen,poten_,pot_matrix(:,ientry))
                  !
                enddo
                !
              endif 
              !
              if (job%IOfitpot_divide) then
                !
                write(jchar, '(i4)') jrot
                write(symchar, '(i4)') isym
                !
                pot_suffix = trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'_'
                !
                call divided_slice_write(i,'potF',pot_suffix,Nentries,pot_matrix)
                !
                !
              else
                !
                write(iunit) i
                write (iunit) pot_matrix
                !
              endif
              !
           enddo
           !
           if (.not.job%IOfitpot_divide) then
              ! 
              write(iunit) 'end potfit'
              close(iunit)
              !
              write(jchar, '(i4)') jrot
              write(symchar, '(i4)') isym
              !
              pot_suffix = trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))
              !
              filename = trim(pot_suffix)//'.chk'
              !
              open(fit(isym,jind)%IOunit,form='unformatted',position='rewind',status='old',action='read',file=filename)
              !
           endif
           !
           if (.not.job%IOextF_divide) then
              !
              if (jind==1.and.isym==1) then 
                !
                read(chkptIO) buf20(1:18)
                if (buf20(1:18)/='End external field') then
                  write (out,"(' prepare_pot_matrix ',a,' has bogus footer: ',a)") job%kinetmat_file,buf20(1:17)
                  stop 'prepare_pot_matrix - bogus file format'
                end if
                !
              endif
              !
              rewind(chkptIO)
              !
           endif
           ! 
           deallocate(pot_matrix,stat=alloc)
           call ArrayStop('pot_matrix')
           !
           if (trim(fitting%method)=='FAST') then
             !
             deallocate(mat_t)
             call ArrayStop('mat_t')
             !
             deallocate(psi)
             call ArrayStop('psi')
             !
           else
             !
             deallocate(tmat)
             call ArrayStop('tmat')
             !
             deallocate(cdimen)
             call ArrayStop('cdimen')
             !
           endif
           !
           if (job%verbose>=5) call TimerReport
           !
         enddo
         !
       end do
       !
       deallocate(poten_)
       !
       if (.not.job%IOextF_divide) then
         close(chkptIO,status='keep')
       endif
       !
    else if (job%IOfitpot_action=='READ') then
       ! 
       do jind = 1, nJ
          !
          if (jind==1.and.fitting%J_list(1)/=0) cycle
          !
          Jrot = Jval(jind)
          !
          do isym= 1,sym%Nrepresen
            !
            if (fit(isym,jind)%Nentries<1) cycle
            !
            if (.not.job%IOfitpot_divide) then
              !
              iunit = fit(isym,jind)%IOunit
              !
              read(iunit) buf20(1:12)
              !
              if (buf20(1:12)/='start potfit') then
                write (out,"(' pot_matrix at J = ',i4,', isym = ',i4,' has bogus footer: ',a)") jrot,isym,buf20(1:12)
                stop 'pot_matrix - bogus file format'
              end if
              !
            endif
            !
         enddo
         !
       end do
       !
    else if (job%IOfitpot_action=='REWRITE') then 
       !
       if (job%verbose>=4) write (out,"(/'rewrite computed mat. elements ...')")
       !
       do jind = 1, nJ
          !
          if (jind==1.and.fitting%J_list(1)/=0) cycle
          !
          Jrot = Jval(jind)
          !
          do isym= 1,sym%Nrepresen
            !
            if (job%verbose>=4) write (out,"('jrot = ',i0,'; sym = ',i0)") Jrot,isym
            !
            Nentries = fit(isym,jind)%Nentries
            !
            if (Nentries<1) cycle
            !
            matsize = int(Nentries*(Nentries+1)/2,hik)
            matsize2= int(Nentries*Nentries,hik)
            !
            dimen = bset_contr(jind)%Maxcontracts
            !
            allocate(pot_matrix(Nentries,Nentries),poten_(Nentries,Nentries),stat=alloc)
            call ArrayStart('pot_matrix',alloc,1,kind(pot_matrix),matsize2)
            call ArrayStart('pot_matrix',alloc,1,kind(poten_),matsize2)
            !
            j = 0
            i = 15
            !
            do 
              !
              i = i + 1
              if (job%verbose>=5) write (out,"('f-'30i4)") molec%pot_ind(:,i)
              !
              j = j + 1
              !
              iunit = fit(isym,jind)%IOunit
              !
              read (iunit,rec=i) pot_matrix
              !
              if (molec%pot_ind(2,i)/=molec%pot_ind(3,i).or.molec%pot_ind(4,i)/=molec%pot_ind(5,i)) then
                !
                i = i + 1
                if (job%verbose>=5) write (out,"('g-'30i4)") molec%pot_ind(:,i)
                !
                read (iunit,rec=i) poten_
                !
                pot_matrix = pot_matrix+poten_
                !
              endif 
              !
              write(unitfname,"('single pot_matrix')")
              !
              call IOStart(trim(unitfname),junit)
              !
              write(jchar, '(i4)') jrot
              write(symchar, '(i4)') isym
              write(parchar, '(i4)') j
              !
              filename = &
                  trim(job%fitpot_file)//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'_'//trim(adjustl(parchar))//'.tmp'
              !
              open(junit,form='unformatted',action='write',position='rewind',status='replace',file=filename)
              !
              write (junit) pot_matrix
              !
              close(junit)
              !
              if (i>=parmax) exit
              !
           enddo
           ! 
           deallocate(pot_matrix,poten_,stat=alloc)
           call ArrayStop('pot_matrix')
           !
           if (job%verbose>=5) call TimerReport
           !
         enddo
         !
       end do
       !
    endif
    !
    if (job%verbose>=2) write (out,"('...done!')")
    !
    call TimerStop('Prepare the potential matrix')
    !
    if (job%verbose>=4) call TimerReport
    !
    !
 end subroutine prepare_pot_matrix


 subroutine divided_slice_write(islice,name,suffix,N,field)
   !
   integer(ik),intent(in)      :: islice
   character(len=*),intent(in) :: name,suffix
   integer(ik),intent(in)      :: N
   real(rk),intent(in)         :: field(N,N)
   character(len=4)            :: jchar
   integer(ik)                 :: chkptIO
   character(len=cl)           :: buf,filename,job_is
   !
   write(job_is,"('single swap_matrix')")
   !
   call IOStart(trim(job_is),chkptIO)
   !
   write(jchar, '(i4)') islice
   !
   filename = trim(suffix)//trim(adjustl(jchar))//'.chk'
   !
   open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=filename)
   !
   write(chkptIO) trim(name)
   !
   write(chkptIO) field
   !
   write(chkptIO) trim(name)
   !
   close(chkptIO)
   !
 end subroutine divided_slice_write
 !       
 !
 subroutine divided_slice_read(islice,name,suffix,N,field)
   !
   integer(ik),intent(in)      :: islice
   character(len=*),intent(in) :: name,suffix
   integer(ik),intent(in)      :: N
   real(rk),intent(out)        :: field(N,N)
   character(len=4)            :: jchar
   integer(ik)                 :: chkptIO
   character(len=cl)           :: buf,filename,job_is
   integer(ik)                 :: ilen
   logical                     :: ifopened
   !
   write(job_is,"('single swap_matrix')")
   !
   call IOStart(trim(job_is),chkptIO)
   !
   write(jchar, '(i4)') islice
   !
   filename = trim(suffix)//trim(adjustl(jchar))//'.chk'
   !
   open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=filename)
   !
   ilen = LEN_TRIM(name)
   !
   read(chkptIO) buf(1:ilen)
   if ( trim(buf(1:ilen))/=trim(name) ) then
     write (out,"(' kinetic divided_slice_read in slice ',a20,': header is missing or wrong',a)") filename,buf(1:ilen)
     stop 'divided_slice_read - in slice -  header missing or wrong'
   end if
   !
   read(chkptIO) field
   !
   read(chkptIO) buf(1:ilen)
   if ( trim(buf(1:ilen))/=trim(name) ) then
     write (out,"(' kinetic divided_slice_read in slice ',a20,': header is missing or wrong',a)") filename,buf(1:ilen)
     stop 'divided_slice_read - in slice -  header missing or wrong'
   end if
   !
   close(chkptIO)
   !
 end subroutine divided_slice_read

 !
end module refinement


!
! defines different transformations of the eigenvectors given 
! in terms of the contracted basis state representaion. 
!
module tran
! set tran_debug > 2 with small vibrational bases and small expansions only
!#define tran_debug  1 

 use accuracy,     only : ik, rk, hik, ark, cl, wl, out, small_
 use timer,        only : IOstart,IOstop,arraystart,arraystop,arrayminus,Timerstart,Timerstop,TimerReport,MemoryReport
 use me_numer,     only : simpsonintegral_ark
 use molecules,    only : MLcoord_direct
 use fields,       only : manifold, FLfingerprint, job,FLNmodes,FLextF_coeffs,FLread_extF_rank,FLextF_matelem,fitting
 use moltype,      only : intensity,extF
 use symmetry,     only : sym

 use perturbation, only : PTintcoeffsT,PTrotquantaT,PTNclasses,PTstore_icontr_cnu,PTeigenT,PTdefine_contr_from_eigenvect,PTrepresT

 private
 public read_contrind,read_eigenval, TReigenvec_unit, bset_contrT, & 
        bset_contr,eigen, index_correlation,Neigenlevels, &
        TRconvert_matel_j0_eigen,TRconvert_repres_J0_to_contr,istate2ilevel

 type bset_contrT
    integer(ik)                 :: jval            ! rotational quantum number that correspond to the contr. basis set
    integer(ik)                 :: Maxsymcoeffs
    integer(ik)                 :: max_deg_size
    integer(ik)                 :: Maxcontracts
    integer(ik)                 :: Nclasses
    integer(ik),        pointer :: icontr2icase(:, :)
    integer(ik),        pointer :: icase2icontr(:, :)
    integer(ik),        pointer :: ktau(:)
    integer(ik),        pointer :: k(:)
    integer(ik),        pointer :: contractive_space(:, :)
    type(PTintcoeffsT), pointer :: index_deg(:)
    type(PTrotquantaT), pointer :: rot_index(:, :)
    integer(ik),        pointer :: icontr_correlat_j0(:, :)
    integer(ik),        pointer :: iroot_correlat_j0(:)
    integer(ik),        pointer :: nsize(:)  ! number of levels in a J,Gamma-block, can be zero for not used symmetry
    integer(ik),        pointer :: nsize_base(:) ! number of levels before current J,gamma-vbloxk, needed to count exomol-ID states
    type(PTrepresT),    pointer :: irr(:)
 end type bset_contrT
 ! 
 type coeffT                      
    real(rk),pointer    :: coeff(:)
    integer(ik),pointer :: icoeff(:)
 end type coeffT
 !
 integer(ik),parameter :: tran_debug = 1
 !
 type(bset_contrT), allocatable,save :: bset_contr(:)  ! information on the contracted basis set
                                                       ! note: bset_contr(1) is always reserved for J=0
 integer(ik),  allocatable           :: i2d_to_1d(:,:) ! a 2d upper diagonal matrix index is transformed into a 1d matrix

 type(PTeigenT), allocatable   :: eigen(:)    ! Eigensolution: description of  eigenvectors 
                                              ! used in the intensity calculations. 
 integer(ik)                 :: Neigenlevels  ! number of 'eigen' levels to be processed 
 integer(ik)                 :: Neigenroots   ! number of 'eigen' roots to be processed 
 integer(ik)                 :: extF_rank     ! 
 integer(ik),pointer         :: TRicontr_cnu(:,:)     ! store cnu  vs icontr index
 integer(ik),pointer         :: Ncontr02icase0(:,:)   ! from icontr0 to icase0, limits of the contracted active space: (i,1) = icontr1, (i,2) = icontr2
 integer(ik),pointer         :: istate2ilevel(:,:,:)  !  correlate the state number in the eigen...chk with ilevel 
 !
 !character(cl), parameter :: eigen_vectr = 'eigen_vectors', eigen_quant = 'eigen_quanta', eigen_descr = 'eigen_descr'

contains

 !
 ! read contraction indexes for J given
 !
 subroutine read_contrind(njval, jval)
    !
    implicit none 
    !
    integer(ik), intent(in) :: njval, jval(njval)

    integer(ik)             :: nclasses, jind, ncontr, ncases, nlambdas, icontr, icase, ilambda, iroot
    integer(ik)             :: iounit, info,Ntotal(sym%Nrepresen),mat_size,igamma,icoeff,ideg

    character(4)            :: jchar
    character(cl)           :: filename, ioname,  buf
    character(len=cl)       :: my_fmt !format for I/O specification
    !
    if (job%verbose>=2) write(out,"(/'Reading the information on the contr. basis...')")
    !
    if (allocated(bset_contr)) deallocate(bset_contr)
    allocate(bset_contr(njval), stat = info)
    if (info /= 0) stop 'read_contr_ind allocation error: bset_contr - out of memory'
    !
    do jind = 1, njval
       !
       !
       bset_contr(jind)%jval = jval(jind)
       !
       !init filename
       !
       write(jchar, '(i4)') jval(jind)
       filename = trim(job%eigenfile%filebase)//'_quanta'//trim(adjustl(jchar))//'.chk'


       if (tran_debug > 2) then
          write(out, '(/a, 1x, i2, 2(1x, a))') 'read contraction indexes for J =', jval(jind), &
                                               'from file', trim(filename)
       end if


       write(ioname, '(a, i4)') 'contraction indexes for J=', jval(jind)
       call iostart(trim(ioname), iounit)
       open(unit = iounit, action = 'read',status='old', file = filename)

       !read vibrational quanta


       read(iounit, '(a)') buf
       if (buf(1:25) /= 'Start Primitive basis set') stop 'read_contrind error: wrong header'


       read(iounit, '(4i8)') ncases, nlambdas, ncontr, nclasses

       bset_contr(jind)%Maxsymcoeffs = ncases
       bset_contr(jind)%max_deg_size = nlambdas
       bset_contr(jind)%Maxcontracts = ncontr
       bset_contr(jind)%nclasses     = nclasses


       if (tran_debug > 2) then
          write(out, '(3(/1x, a, 1x, i5))') 'ncases ', ncases, 'ndegmax', nlambdas, 'ncontr ', ncontr
          write(out, '(/1x, a, 1x, a, 4x, a/40x, 1x, a, 1x, a, 6x, a)') 'icase', 'ndeg', 'ilevel(0:nclasses)', 'ideg', &
                                                                        'iroot', 'ideg(0:nclasses)'
       end if
       !
       allocate(bset_contr(jind)%index_deg(ncases),bset_contr(jind)%contractive_space(0:nclasses, ncases),stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%contractive_space),kind(bset_contr(jind)%contractive_space))
       !
       allocate(bset_contr(jind)%nsize(sym%Nrepresen),stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%nsize),kind(bset_contr(jind)%nsize))
       !
       allocate(bset_contr(jind)%nsize_base(sym%Nrepresen),stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%nsize_base),kind(bset_contr(jind)%nsize_base))
       !
       allocate(bset_contr(jind)%icontr2icase(ncontr, 2),bset_contr(jind)%icase2icontr(ncases, nlambdas),stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%icontr2icase),kind(bset_contr(jind)%icontr2icase))
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%icase2icontr),kind(bset_contr(jind)%icase2icontr))
       !
       icontr = 0
       iroot = 0 
       !
       do icase = 1, ncases
          !
          read(iounit,*) nlambdas, bset_contr(jind)%contractive_space(0:nclasses, icase)
          !
          bset_contr(jind)%index_deg(icase)%size1 = nlambdas
          !
          if (tran_debug > 2) then
             write(my_fmt,'(a,i0,a)') "(x, i5, 2x, i3, 1x,",nclasses,"(1x, i3), 1x, i3)"
             write(out,my_fmt)   &
             icase, nlambdas, bset_contr(jind)%contractive_space(0:nclasses, icase)
          end if
          !
          allocate(bset_contr(jind)%index_deg(icase)%icoeffs(0:nclasses, nlambdas))
          call ArrayStart('bset_contr',info,size(bset_contr(jind)%index_deg(icase)%icoeffs),&
                          kind(bset_contr(jind)%index_deg(icase)%icoeffs))
          !
          do ilambda = 1, nlambdas
             !
             write(my_fmt,'(a,i0,a)') "(i8,i6,",nclasses,"i6)"
             read(iounit,my_fmt) iroot, bset_contr(jind)%index_deg(icase)%icoeffs(0:nclasses, ilambda)
             !
             if (tran_debug > 2) then
                write(my_fmt,'(a,i0,a)') "(42x, i3, 1x, i5, 3x,",nclasses,"(1x, i3), 1x, i3)"
                write(out,my_fmt) ilambda, iroot, bset_contr(jind)%index_deg(icase)%icoeffs(0:nclasses, ilambda)
             end if
             !
             icontr = icontr + 1
             !
             bset_contr(jind)%icontr2icase(icontr, 1)      = icase
             bset_contr(jind)%icontr2icase(icontr, 2)      = ilambda
             bset_contr(jind)%icase2icontr(icase, ilambda) = icontr

             if (iroot /= icontr) stop 'read_contrind error: wrong indexing'

          end do

       end do

       if (icontr /= ncontr) stop 'read_contrind error: wrong indexing'


       !read rotational quanta


       read(iounit, '(2i8)') ncases, nlambdas


       if (tran_debug > 2) then
          write(out, '(2(/1x, a, 1x, i2))') 'nrot      ', ncases, 'nrotdegmax', nlambdas
          write(out, '(/1x, a, 1x, a, 8x, a, 3x, a, 1x, a)') 'irot', 'irotdeg', 'J', 'K', 'Tau'
       end if


       allocate(bset_contr(jind)%rot_index(ncases, nlambdas), stat = info)
       if (info /= 0) stop 'read_contrind allocation error: bset_contr%rot_index - out of memory'

       do
          read(iounit, '(a)') buf
          if (buf(:3) == 'End') exit

          read(buf, *) icase, ilambda, bset_contr(jind)%rot_index(icase, ilambda)%j,                                       &
                                       bset_contr(jind)%rot_index(icase, ilambda)%k,                                       &
                                       bset_contr(jind)%rot_index(icase, ilambda)%tau

          if (tran_debug > 2) then
             write(out, '(2x, i3, 5x, i3, 5x, 3(1x, i3))') icase, ilambda, bset_contr(jind)%rot_index(icase, ilambda)%j,   &
                                                                           bset_contr(jind)%rot_index(icase, ilambda)%k,   &
                                                                           bset_contr(jind)%rot_index(icase, ilambda)%tau
          end if

       end do

       !read(iounit, '(24a)') buf(1:24)
       if (buf(1:24) /= 'End Primitive basis set') stop 'read_contrind error: wrong footer'
       !
       if (job%IOvector_symm) then
         !
         read(iounit, '(32a)') buf(1:32)
         if (buf(1:32) /= 'Start irreducible transformation') stop 'read_contrind error: wrong sym-footer'
         !
         allocate(bset_contr(jind)%irr(sym%Nrepresen),stat = info)
         call ArrayStart('bset_contr',info,1,1)
         !
         read(iounit,*) mat_size
         !
         read(iounit,*) Ntotal(1:sym%Nrepresen)
         !
         do igamma = 1,sym%Nrepresen
           !
           allocate(bset_contr(jind)%irr(igamma)%N(bset_contr(jind)%Maxsymcoeffs),&
                    bset_contr(jind)%irr(igamma)%repres(Ntotal(igamma),sym%degen(igamma),mat_size),stat = info)
           call ArrayStart('bset_contr',info,size(bset_contr(jind)%irr(igamma)%repres),&
                    kind(bset_contr(jind)%irr(igamma)%repres))
           !
           do icoeff = 1,bset_contr(jind)%Maxsymcoeffs
              !
              read(iounit,*) bset_contr(jind)%irr(igamma)%N(icoeff)
              !
           enddo
           !
           do icoeff = 1,Ntotal(igamma)
              !
              do ideg = 1,sym%degen(igamma)
                 !
                 read(iounit,*) bset_contr(jind)%irr(igamma)%repres(icoeff,ideg,1:mat_size)
                 !
              enddo
              !
           enddo
           !
         enddo
         !
         read(iounit, '(30a)') buf(1:30)
         !
         if (buf(1:30) /= 'End irreducible transformation') stop 'read_contrind error: wrong irrep-footer'
         !
       endif 
       !
       close(unit = iounit)
       !
    end do
    !
    if (job%verbose>=2) write(out,"('...done!')")

 end subroutine read_contrind

 !
 !find correlation between contracted basis indexes for J = 0 and J > 0
 !
 subroutine index_correlation(njval, jval)
    implicit none 
    !

    integer(ik), intent(in) :: njval, jval(njval)

    integer(ik)             :: jind, icase, jcase, ilambda, jlambda,icontr,jcontr, nclasses,info, iroot,jroot
    integer(ik)             :: ilevel,ideg,k,tau,dimen,irow,icol
    !
    integer(ik),allocatable  :: cnu_i(:),cnu_j(:)
    !
    logical                 :: found
    character(len=cl)       :: my_fmt   !format for I/O specification
    !
    if (job%verbose>=2) write(out,"(/'Establish the correlation between the indexes of J=0 and J>0 contr. basis funct.')")
    !

    if (.not. allocated(bset_contr)) stop 'index_correlation error: allocated(bset_contr) = .false.'

    !
    ! In order to find the correlation of the vibrational quantum numbers at J=0 and J/=0 
    ! the J=0 quanta destibution has to be defined before. We assume here that this is the case, 
    ! bset_contr(1) is always reserved for J=0
    !
    if (bset_contr(1)%jval/=0) then 
      write(out,"('index_correlation: bset_contr(1) is not for J=0')")
      stop 'index_correlation: bset_contr(1) has to be reserved for J=0'
    endif 
    !
    do jind = 1, njval
       !
       nclasses = bset_contr(jind)%nclasses
       if (bset_contr(jind)%nclasses/=bset_contr(1)%nclasses) then 
         write(out,"('index_correlation: Nclasses are different for diff. J:',2i0)") &
                     bset_contr(1)%nclasses,bset_contr(jind)%nclasses
         stop 'index_correlation: Nclasses cannot be different for diff. J'
       endif
       !
       allocate(cnu_i(1:nclasses),cnu_j(1:nclasses),stat = info)
       if (info /= 0) stop 'index_correlation: cnu_i allocation error'
       !
       if (tran_debug > 2) then
          write(out, '(/a, 1x, i2)') 'find correlation between contraction indexes for J = 0 and J =', jval(jind)
       end if
       !
       if (jind > size(bset_contr)) stop 'index_correlation error: jind > size(bset_contr)'
       !
       allocate(bset_contr(jind)%icontr_correlat_j0(bset_contr(jind)%Maxsymcoeffs,bset_contr(jind)%max_deg_size), stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%icontr_correlat_j0),kind(bset_contr(jind)%icontr_correlat_j0))
       !
       if (tran_debug > 2) then
          write(out, '(/7x, a, i3, 17x, a/1x, a, 1x, a, 1x, a, 4x, a, 1x, a, 1x, a)')                                  &
          'J =', jval(jind), 'J =  0', 'icase', 'ilambda', 'iroot', 'icase', 'ilambda', 'iroot'
       end if
       !
       l_icase : do icase = 1, bset_contr(jind)%Maxsymcoeffs
          !
          cnu_i(1:nclasses) = bset_contr(jind)%contractive_space(1:nclasses, icase)
          !
          l_ilambda : do ilambda = 1, bset_contr(jind)%index_deg(icase)%size1
             !
             found = .false. 
             !
             l_jcase : do jcase = 1, bset_contr(1)%Maxsymcoeffs
               !
               cnu_j(1:nclasses) = bset_contr(1)%contractive_space(1:nclasses, jcase)
               !
               l_jlambda : do jlambda = 1, bset_contr(1)%index_deg(jcase)%size1
                 !
                 if (all(cnu_i(:) == cnu_j(:))  .and.   &
                     all(bset_contr(   1)%index_deg(jcase)%icoeffs(1:nclasses,jlambda) == &
                         bset_contr(jind)%index_deg(icase)%icoeffs(1:nclasses,ilambda))) then 
                         !
                         found = .true.
                         exit l_jcase
                         !
                 endif 
                 !
               end do l_jlambda
             end do l_jcase
             !
             if (.not.found) then 
               write(my_fmt,'(a,i0,a)') "(a,2i6,",nclasses,"i4,i4)"
               write(out,"('index_correlation: not found for J = ',i8,' -> problems with checkpoints?')") jval(jind)
               write(out,my_fmt) 'J,icase,cnu,ideg:',jval(jind),icase,cnu_i(:),ilambda
               stop 'index_correlation: not found'
             endif 
             !
             jcontr = bset_contr(1)%icase2icontr(jcase,jlambda)
             !
             bset_contr(jind)%icontr_correlat_j0(icase, ilambda) = jcontr
             !
             if (tran_debug > 2) then
                write(out, '(i6, 6x, i2, i6, 3x, i6, 6x, i2, i6)')   &
                !                      case,                      lambda, contr
                icase,ilambda,bset_contr(jind)%icase2icontr(icase,ilambda), &
                bset_contr(1)%icontr2icase(jcontr,1), bset_contr(1)%icontr2icase(jcontr,2), jcontr
                !
             end if
             !
          end do l_ilambda
          !
       end do l_icase
       !
       ! The same but for a 1d index iroot in place of (icase,lambda)
       !
       allocate(bset_contr(jind)%iroot_correlat_j0(bset_contr(jind)%Maxcontracts), stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%iroot_correlat_j0),kind(bset_contr(jind)%iroot_correlat_j0))
       allocate(bset_contr(jind)%ktau(bset_contr(jind)%Maxcontracts), stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%ktau),kind(bset_contr(jind)%ktau))
       allocate(bset_contr(jind)%k(bset_contr(jind)%Maxcontracts), stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%k),kind(bset_contr(jind)%k))
       !
       do iroot = 1, bset_contr(jind)%Maxcontracts
          !
          icase   = bset_contr(jind)%icontr2icase(iroot, 1)
          ilambda = bset_contr(jind)%icontr2icase(iroot, 2)
          !
          jcontr = bset_contr(jind)%icontr_correlat_j0(icase, ilambda)
          !
          bset_contr(jind)%iroot_correlat_j0(iroot) = jcontr
          !
          ilevel  = bset_contr(jind)%contractive_space(0,icase)
          ideg    = bset_contr(jind)%index_deg(icase)%icoeffs(0,ilambda)
          !
          k      = bset_contr(jind)%rot_index(ilevel,ideg)%k
          tau    = bset_contr(jind)%rot_index(ilevel,ideg)%tau
          !
          bset_contr(jind)%ktau(iroot) = 2*k+tau
          bset_contr(jind)%k(iroot)    = k
          !
          if (job%rotsym_do) then 
            bset_contr(jind)%ktau(iroot) = ilevel
            bset_contr(jind)%k(iroot)    = ideg
          endif 
          !
          if (tran_debug >= 3) then
            write (out,"('iroot,icase,ilambda,jcontr = ',4i7)") iroot,icase,ilambda,jcontr 
          end if
          !
       enddo
       !
       if (tran_debug > 2) then
          write(out, '(/a)') 'done'
       end if
       !
       deallocate(cnu_i,cnu_j)
       !
    end do
    !
    ! Introduce the address-index matrix for correspondence between 1d and 2d (i<j) index arrays
    !
    dimen  = bset_contr(1)%Maxcontracts
    !
    allocate(i2d_to_1d(dimen,dimen), stat = info)
    call ArrayStart('i2d_to_1d',info,1,kind(i2d_to_1d),size(i2d_to_1d,kind=hik))
    !
    do icontr = 1,dimen
      !
      do jcontr = 1,dimen
        !
        irow = max(icontr, jcontr)
        icol = min(icontr, jcontr)
        !
        i2d_to_1d(icontr,jcontr) = ( irow * (irow - 1) )/ 2 + icol
        !
      enddo
    enddo
    !
    if (job%verbose>=2) write(out,"('...done!')")
    !
 end subroutine index_correlation


 ! 
 !read eigenvalues and their assignment 
 !
 subroutine read_eigenval(njval, jval, error)
    !
    implicit none
    !
    integer(ik), intent(in) :: njval, jval(njval)
    integer(ik),optional,intent(out) :: error

    integer(ik)             :: jind, nmodes, nroots, ndeg, nlevels,  iroot, irec, igamma, ilevel, jlevel, &
                               ideg, ilarge_coef,k0,tau0,nclasses,nsize,nsize_base,id_,j_,   &
                               iounit, jounit, info, quanta(0:FLNmodes), iline, nroots_t, nu(0:FLNmodes),&
                               normal(0:FLNmodes),Npolyad_t
    integer(ik),allocatable :: ktau_rot(:,:),isym(:)
    !
    real(rk)                :: energy,energy_t,largest_coeff
    !
    character(cl)           :: filename, ioname, buf
    character(4)            :: jchar,gchar
    character(500)          :: buf500
    !
    logical                 :: passed
    logical                 :: normalmode_input = .false.,largest_coeff_input = .false.
    integer(ik)             :: jind_t,maxdeg,gamma,jval_,irec_, igamma_, ilevel_
    real(rk)                :: energy_, state_intensity
    !
    type(PTeigenT)          :: eigen_t   ! temporal object used for sorting 'eigen'
    character(len=wl)       :: my_fmt !format for I/O specification
    ! 
    if (job%verbose>=2) write(out,"(/'Read and sort eigenvalues in increasing order...')")
    !
    call TimerStart('Read eigenvalues')

    if (.not. allocated(bset_contr)) stop 'read_eigenval error: associated(bset_contr) = .false.'
    !
    nmodes = FLNmodes
    !nclasses = PTNclasses
    !
    nclasses = bset_contr(1)%nclasses
    !
    ! In practice we do not need all stored rootes, but only those that pass the 
    ! filters given in the input. Here we perform a pre-selection of the states, 
    ! which are actually required in the calculations. Towards this we read the unit twice:
    ! a) to count and b) to read.  
    ! After reading all states will  be combined together and sorted wrt energy increasing 
    ! to produce just one list of states. 
    !
    ! Astimate the maximal number of energy records
    !
    nroots_t = 0
    !
    do jind = 1, njval
      nroots_t = max(bset_contr(jind)%Maxcontracts,nroots_t)
    enddo 
    !
    ! create a temp. array needed for filtering out levels
    !
    allocate(isym(0:nclasses),ktau_rot(0:2*maxval( jval(:),dim=1 ),2)) 
    !
    allocate(istate2ilevel(njval,sym%Nrepresen,nroots_t),stat=info) 
    !
    call ArrayStart('istate2ilevel',info,size(istate2ilevel),kind(istate2ilevel))
    !
    istate2ilevel = 0
    !
    ! initialize the error code which will indicate if an eigen-file exists
    if (present(error)) error = 0 
    !
    nroots  = 0
    nlevels = 0 
    maxdeg = 1
    !
    do jind = 1, njval
       !
       nsize_base = 0
       do gamma = 1,sym%Nrepresen
          !
          if (.not.job%select_gamma(gamma).and.(jval(jind)/=0.or.gamma/=1)) cycle
          !
          write(jchar, '(i4)') jval(jind)
          write(gchar, '(i3)') gamma
          !
          filename = trim(job%eigenfile%filebase)//'_descr'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
          !
          if (tran_debug > 2) then
             write(out, '(/a, 1x, i2, 2(1x, a))') 'read eigenvalues for J =', jval(jind),&
                         ', gamma = ',gamma,' from file', trim(filename)
          end if
          !
          if (jind > size(bset_contr)) stop 'read_eigenval error: jind > size(bset_contr)'
          !
          write(ioname, '(a, i4,2x,i4)') 'eigenvalues for J,gamma = ', jval(jind),gamma
          !
          call IOstart(trim(ioname), iounit)
          open(unit = iounit, action = 'read',status='old' , file = filename, err=15)
          !
          ! for the TM-based basis set pruning open the chk-file with the vibrational intensities 
          !
          if (job%TMpruning)  then 
             !
             filename = trim(job%eigenfile%filebase)//'_intens'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
             !
             write(ioname, '(a, i4,2x,i4)') 'J=0 intensities for J,gamma = ', jval(jind),gamma
             !
             call IOstart(trim(ioname), jounit)
             open(unit = jounit, action = 'read',status='old' , file = filename)
             !
          endif
          !
          ! Check the fingerprint of the computed eigenvectors. 
          !
          call FLfingerprint('read',iounit,0,job%Npolyads_prim,(/job%enercutoff%primt,job%enercutoff%contr/))
          !
          buf500 = ''
          !
          read(iounit, '(a)') buf500
          read(iounit,"(i8,a4)") Npolyad_t,buf500(1:4)
          if (trim(buf500(1:3))==' <=')  normalmode_input = .true.
          if (trim(buf500(1:4))==' <==') largest_coeff_input = .true.
          read(iounit, '(a)') buf500
          !
          !do while (buf500(1:5) /= 'Class')
          !   read(iounit, '(a)') buf500
          !end do
          !
          ! Start reading the description of the eigensolution. 
          !
          read(iounit,*) nroots_t,nsize
          !
          bset_contr(jind)%nsize_base(gamma) = nsize_base + bset_contr(1)%Maxcontracts*jval(jind)**2
          !
          nsize_base = nsize_base + nsize
          !
          if (.not.job%select_gamma(gamma).and.(jval(jind)/=0.or.gamma/=1)) then
            close(iounit)
            cycle
          endif
          !
          if ( .not.job%IOvector_symm.and.nroots_t /= bset_contr(jind)%Maxcontracts) then
             write(out,"('read_eigenval error: wrong number of contracted solutions')")
             stop 'read_eigenval error: wrong number of contracted solutions'
          endif
          !
          do
             !
             read(iounit, '(a)') buf500
             if (buf500(1:3) == 'End') exit
             !
             read(buf500, *) irec, igamma, ilevel, ideg, energy, nu(0:nmodes)
             !
             if (igamma/=gamma) then 
               write(out,"('read_eigenval error: igamma/=gamma: ',2i4)") igamma,gamma
               stop 'read_eigenval error: igamma/=gamma'
             endif
             !
             if (job%ZPE<0.and.igamma==1.and.Jval(jind)==0) job%zpe = energy
             if (intensity%ZPE<0.and.igamma==1.and.Jval(jind)==0) intensity%zpe = energy
             !
             call filter(energy,igamma,passed)
             !
             if (job%TMpruning)  then 
                !
                read(jounit, *) jval_, igamma_, ilevel_, energy_, state_intensity
                !
                if (jval_/=jval(jind).or.igamma_/=igamma.or.ilevel_/=ilevel) then 
                  !
                  write(out,"('eigen_intens.chk error: ineteger records do not agree with eigen_descr:' )")
                  write(out,"('jind =  ',2i5,', igamma = ',2i4,'  ilevel =  ',2i8 )") jval_,jval(jind),igamma_,igamma,ilevel_,ilevel
                  stop 'eigen_intens.chk error: records do not agree with eigen_descr'
                  !
                endif
                !
                if (energy_>small_.and.abs(energy-energy_)>sqrt(small_)) then 
                  !
                  write(out,"('eigen_intens.chk error: energy  does not agree with eigen_descr:' )")
                  write(out,"('energy =  ',2f20.12)") energy_,energy
                  stop 'eigen_intens.chk error: energy  does not agree with eigen_descr'
                  !
                endif
                !
                if (state_intensity<job%TMcutoff.and.energy-job%ZPE>job%TMenermin) passed = .false.
                !
             endif
             !
             if (passed) then 
                if (ideg==1) then 
                    !
                    nlevels = nlevels + 1
                    istate2ilevel(jind,igamma,ilevel) = nlevels
                    !
                endif
                !
                if (job%IOvector_symm) then
                  !
                  nroots = nroots + sym%degen(igamma)
                  maxdeg = max(maxdeg,sym%degen(igamma))
                  !
                else
                  !
                  nroots = nroots + 1
                  maxdeg = max(maxdeg,ideg)
                  !
                endif
                !
             endif 
             !
          end do
          !
          close(iounit)
          if (job%TMpruning) close(jounit)
          !
          cycle
          !
          15 continue 
          !
          call IOStop(trim(ioname))
          !
          write(out,"('read_eigenval warninig: eigenfilefile ',a,'is missing')") filename
          if (present(error)) error = 1 
          !
       enddo
    enddo
    !
    if (nroots==0) then 
       write(out,"('read_eigenval: the filters are too tight: no entry')") 
       stop 'read_eigenval: the filters are too tight' 
    endif 
    !
    ! total number of levels to be processed, a global variable 
    !
    Neigenlevels = nlevels
    Neigenroots = nroots
    !
    allocate(eigen(Neigenlevels),stat = info)
    if (info /= 0) stop 'read_eigenval allocation error: eigen - out of memory'
    !
    if (job%verbose>=4) then 
       write(out,"('   Number of selected eigenvectors: ',i8)") Neigenlevels
    end if
    !
    do ilevel = 1,Neigenlevels
      !
      allocate(eigen(ilevel)%irec(maxdeg),eigen(ilevel)%iroot(maxdeg),eigen(ilevel)%quanta(0:nmodes),&
               eigen(ilevel)%normal(0:nmodes),eigen(ilevel)%cgamma(0:nclasses), stat = info)
      if (info /= 0) stop 'read_eigenval allocation error: eigen%irec, eigen%quanta - out of memory'
      eigen(ilevel)%ndeg   = 0
      eigen(ilevel)%iroot = 0
      eigen(ilevel)%quanta = 0
      !
    enddo
    !
    if (job%exomol_format) then
      write(out,"(/a/)") 'States file in the Exomol format'
    endif
    !
    ! Now we can actually read and store the energies and their description. 
    !
    nroots  = 0
    nlevels = 0 
    !
    do jind = 1, njval
       !
       ! reconstruct k_rot and tau_rot from a 1d distribution
       !
       ktau_rot(0,1) = 0
       ktau_rot(0,2) = mod(Jval(jind),2)
       !
       iroot = 0
       !
       do k0 = 1,Jval(jind)
         !
         do tau0 = 0,1
           !
           iroot = iroot + 1
           ktau_rot(iroot,1) = k0
           ktau_rot(iroot,2) = tau0
           !
         enddo 
         !
       enddo
       !
       do gamma = 1,sym%Nrepresen
          !
          bset_contr(jind)%nsize(gamma) = 0
          !
          if (.not.job%select_gamma(gamma).and.(jval(jind)/=0.or.gamma/=1)) cycle
          !
          write(jchar, '(i4)') jval(jind)
          write(gchar, '(i3)') gamma
          !
          filename = trim(job%eigenfile%filebase)//'_descr'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
          !
          if (jind > size(bset_contr)) stop 'read_eigenval error: jind > size(bset_contr)'
          !
          write(ioname, '(a, i4,2x,i4)') 'eigenvalues for J,gamma = ', jval(jind),gamma
          !
          call IOstart(trim(ioname), iounit)
          open(unit = iounit, action = 'read',status='old' , file = filename, err=16)
          !
          buf500 = ''
          !
          do while (buf500(1:5) /= 'Class')
             read(iounit, '(a)') buf500
          end do
          !
          ! Number of roots and the size 
          !
          read(iounit,*) nroots_t,nsize
          !
          bset_contr(jind)%nsize(gamma) = nsize
          !
          do
             !
             read(iounit, '(a)') buf500
             if (buf500(1:3) == 'End') exit
             !
             largest_coeff = 1.0_rk
             !
             if (normalmode_input.and.largest_coeff_input) then
               !
               read(buf500, *) irec, igamma, ilevel, ideg, energy, quanta(0:nmodes), ilarge_coef,isym(0:nclasses),&
                               normal(0:nmodes),largest_coeff
               !
             elseif (normalmode_input) then 
               !
               read(buf500, *) irec,igamma,ilevel,ideg,energy,quanta(0:nmodes),ilarge_coef,isym(0:nclasses),normal(0:nmodes)
               !
             elseif (largest_coeff_input) then 
               !
               ! obsolete
               !
               read(buf500, *) irec, igamma, ilevel, ideg, energy, quanta(0:nmodes), ilarge_coef,isym(0:nclasses),largest_coeff
               !
               normal = quanta
               !
             else
               !
               ! obsolete
               !
               read(buf500, *) irec, igamma, ilevel, ideg, energy, quanta(0:nmodes), ilarge_coef,isym(0:nclasses)
               !
               normal = quanta
               !
             endif 
             !
             !if (job%ZPE<0.and.igamma==1.and.Jval(jind)==0) job%zpe = energy
             !
             passed = .true.
             if (job%triatom_sing_resolve) then
                if (Jval(jind)==0.and.normal(0)/=0) then
                  passed = .false.
                endif
             endif 
             !
             if (job%exomol_format.and.(jind/=1.or.intensity%J(1)==0).and.passed) then
               !
               !write(out,"(/a/)") 'States file in the Exomol format'
               !
               ID_ = ilevel + bset_contr(jind)%nsize_base(gamma)
               !
               J_ = Jval(jind)
               !
               !write(out,"(i12,1x,f12.6,1x,i6,1x,i7,2x,a3,2x,<nmodes>i3,1x,<nclasses>(1x,a3),1x,2i4,1x,a3,2x,f5.2,a3,1x,i9,1x,<nmodes>i3)") & 
               !
               write(my_fmt,'(a,i0,a,i0,a,i0,a)') &
                     "(i12,1x,f12.6,1x,i6,1x,i7,2x,a3,2x,",nmodes,"i3,1x",nclasses,"(1x,a3),2i4,1x,a3,2x,f5.2,a3,1x,i9,1x",nmodes,"i3)"
               !
               write(out,my_fmt) & 
               ID_,energy-intensity%ZPE,int(intensity%gns(gamma),4)*(2*J_+1),J_,sym%label(gamma),&
               normal(1:nmodes),sym%label(isym(1:nclasses)),&
               ktau_rot(quanta(0),1),ktau_rot(quanta(0),2),sym%label(isym(0)),&
               largest_coeff,' ::',ilevel,quanta(1:nmodes)
               !
             endif
             !
             passed = .true.
             !
             if (istate2ilevel(jind,igamma,ilevel)==0) passed = .false.
             !
             !call filter(energy,igamma,passed)
             !
             if (passed) then
                !
                ! 'nlevels' runs over the levels (i.e. one 'ilevel' for all denerate components)
                !
                nlevels = istate2ilevel(jind,igamma,ilevel)
                !
                eigen(nlevels)%irec(ideg)  = irec
                !
                if (job%IOvector_symm) then
                  !
                  ! 'nroots' runs over all degerate states 
                  !
                  nroots = nroots + sym%degen(igamma)
                  !
                  eigen(nlevels)%ndeg  = sym%degen(igamma)
                  !
                else
                  !
                  ! 'nroots' runs over all degerate states 
                  !
                  nroots = nroots + 1
                  !
                  eigen(nlevels)%ndeg        = eigen(nlevels)%ndeg + 1
                  !
                endif
                !
                if (ideg==1) then  
                  !
                  eigen(nlevels)%jind       = jind
                  eigen(nlevels)%jval       = Jval(jind)
                  eigen(nlevels)%ilevel     = ilevel
                  eigen(nlevels)%krot       = ktau_rot(quanta(0),1)
                  eigen(nlevels)%taurot     = ktau_rot(quanta(0),2)
                  eigen(nlevels)%energy     = energy
                  eigen(nlevels)%igamma     = igamma
                  eigen(nlevels)%quanta(:)  = quanta(:)
                  eigen(nlevels)%normal(:)  = normal(:)
                  eigen(nlevels)%cgamma(:)  = sym%label(isym(:))
                  eigen(nlevels)%icoeff     = ilarge_coef
                  eigen(nlevels)%largest_coeff = largest_coeff
                  !
                endif 
                !
             endif 
             !
          end do
          !
          close(unit = iounit)
          !
          !
          !print energies
          if (tran_debug > 2) then
             !
             write(my_fmt,'(a,i0,a)') "(/1x, a, 11x, a, 1x, a, 1x, a, 8x, a,",nmodes,"(2x), 1x, a)"
             write(out, '(/1x, a, 2x, i8/1x, a, 1x, i8)') 'number of roots', nroots, 'number of levels', nlevels
             write(out,my_fmt) 'ilevel', 'energy', 'ndeg','igamma', 'nu(0:nmodes)', 'irec'
             !
             write(my_fmt,'(a,i0,a,i0,a)') "(2x,i8,1x,f14.8,2x,i3,4x,i3, 5x,",nmodes,"(1x, i3),1x,i7,5x",ndeg,"(1x, i7)))"
             !
             do ilevel = 1, nlevels
                !
                ndeg = eigen(ilevel)%ndeg
                write(out,my_fmt)  &
                ilevel, eigen(ilevel)%energy, eigen(ilevel)%ndeg,        &
                eigen(ilevel)%igamma, eigen(ilevel)%quanta(1:nmodes),    &
                eigen(ilevel)%irec(1:ndeg)
                !
             end do
             write(out, '(a)') '...done!'
          end if
          !
          cycle
          16 continue 
          if (present(error)) error = 1 
          !
       enddo
       !
    end do
    !
    ! Sort according the energy increasing 
    !
    allocate(eigen_t%irec(maxdeg),eigen_t%iroot(maxdeg),eigen_t%quanta(0:nmodes),eigen_t%cgamma(0:nmodes), stat = info)
    if (info /= 0) stop 'read_eigenval allocation error: eigen_t%irec,quanta - out of memory'
    !
    do ilevel =1,nlevels
      !
      energy = eigen(ilevel)%energy
      !
      do jlevel =ilevel+1,nlevels
        !
        if (energy>eigen(jlevel)%energy) then 
          !
          eigen_t       = eigen(jlevel)
          eigen(jlevel) = eigen(ilevel)
          eigen(ilevel) = eigen_t
          energy        = eigen(ilevel)%energy
          !
          istate2ilevel(eigen(ilevel)%jind,eigen(ilevel)%igamma,eigen(ilevel)%ilevel) = ilevel
          istate2ilevel(eigen(jlevel)%jind,eigen(jlevel)%igamma,eigen(jlevel)%ilevel) = jlevel
          !
        endif 
        !
      enddo
      !
    enddo
    !
    iroot = 0 
    !
    do ilevel =1,nlevels
      do ideg = 1,eigen(ilevel)%ndeg
        !
        iroot = iroot +1 
        eigen(ilevel)%iroot(ideg) = iroot 
        !
      enddo
    enddo
    !
    deallocate(ktau_rot,isym) 
    !
    if (job%verbose>=2) write(out,"('...done!')")
    !
    call TimerStop('Read eigenvalues')
    !
    contains 
      !
      subroutine filter(energy,igamma,passed)
        !
        implicit none
        !
        real(rk),intent(in)    :: energy
        integer(ik),intent(in) :: igamma
        logical,intent(out)    :: passed
          !
          ! passed = .true.
          !
          ! if (.not.intensity%do) return
          !
          passed = .false.
          !
          if (job%isym_do(igamma).and.energy-job%ZPE>=job%erange(1)  &
                                 .and.energy-job%ZPE<=job%erange(2)) then 
              !
              passed = .true.
              !
          endif 

      end subroutine filter
      !
 end subroutine read_eigenval





 function TReigenvec_unit(jind,jval,igamma)

 !return file-units for reading eigenvectors

    integer(ik), intent(in)          :: jind,jval(:)
    integer(ik),optional, intent(in) :: igamma
    !
    integer(ik)             :: TReigenvec_unit
    !
    integer(ik)             :: kind,ilevel,jlevel, ncontr, iounit, info, reclen, irec
    !
    real(rk), pointer       :: vec1(:),vec2(:)
    real(rk)                :: f_t
    !
    character(4)            :: jchar,gchar = 'xxxx'
    character(cl)           :: filename, ioname
    logical                 :: exists,hasopened
    !
    if (present(igamma)) then
       !
       if (.not.job%select_gamma(igamma).and.(jval(jind)/=0.or.igamma/=1)) then 
         TReigenvec_unit = -1
         return
       endif
       !
       write(ioname, '(a, i4,2x,i4)') 'eigenvalues for J,gamma = ', jval(jind),igamma
       !
    else
       write(ioname, '(a, i4)') 'eigenvectors for J=', jval(jind)
       !
       if (job%IOvector_symm) then
          write(out,"(' TReigenvec_unit: igamma must be present for EIGEN_SYMM')")
          stop 'TReigenvec_unit: igamma must be present for EIGEN_SYMM'
       endif
       !
    endif
    !
    call iostart(trim(ioname), iounit)
    !
    TReigenvec_unit = iounit
    !
    write(jchar, '(i4)') jval(jind)
    if (present(igamma)) write(gchar, '(i3)') igamma
    !
    !filename = trim(job%eigenfile%filebase)//'_vectors'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
    !if (job%IOvector_symm) 
    !
    filename = trim(job%eigenfile%filebase)//'_vectors'//trim(adjustl(jchar))//'_'//trim(adjustl(gchar))//'.chk'
    !
    inquire (file = filename, exist = exists)
    if (.not.exists) then 
      write (out,"('TReigenvec_unit: Cannot find file for j= ',i3,4x,a)") jval(jind),filename
      stop 'TReigenvec_unit: file with eigenvectors does not exist' 
    endif 
    !
    inquire (unit = iounit, opened = hasopened)
    !
    if (hasopened) return
    !
    !write(jchar, '(i4)') jval(jind)
    !filename = trim(job%eigenfile%filebase)//'_vectors'//trim(adjustl(jchar))//'.chk'
    !
    !if (job%IOvector_symm) then 
    !  mat_size = Ntotal(igamma)
    !  ib = Ntotal(igamma)
    !else
    !  mat_size = PT%Maxprimitive
    !  ib = sum(PT%Max_sym_levels(:)*sym%degen(:))
    !endif 
    !
    if (job%IOvector_symm) then 
      ncontr = max(bset_contr(jind)%nsize(igamma),1)
    else
      ncontr = bset_contr(jind)%Maxcontracts
    endif 
    !
    inquire(iolength = reclen) f_t
    reclen = ncontr*reclen
    !
    open(unit = iounit, access = 'direct', recl = reclen, action='read',status='old' , file = filename,err=22)
    !
    return
    !
    22 continue
    !
    write (out,"('TReigenvec_unit: Error opening  eigenvectors-file for j= ',i4,a)") jval(jind),filename
    stop 'TReigenvec_unit: Error opening  eigenvectors-file' 
    !
    ! test reading
    !
    if (job%verbose >= 7) then 
      !
      if (job%verbose>=2) write(out,"('Checking orthogonality of eigenfunctions for J = ',i4,' ...')") jval(jind)
      !
      ncontr = bset_contr(jind)%Maxcontracts
      !
      !omp parallel private(vec,info)
      allocate(vec1(ncontr),vec2(ncontr), stat = info)
      if (info /= 0) stop 'TReigenvec_unit allocation error: vec - out of memory'
      !
      !omp do private(ilevel,kind,irec,f_t) schedule(guided)
      do ilevel = 1,Neigenlevels
         !
         kind = eigen(ilevel)%jind
         !
         if (jind/=kind) cycle
         !
         irec = eigen(ilevel)%irec(1)
         !
         read(iounit, rec = irec) vec1(1:ncontr)
         !
         !check normalization
         !
         f_t = dot_product(vec1(:),vec1(:))
         if(abs(f_t - 1.0_rk) > 0.1_rk ** (rk - 1)) then
            write(out, '(1x, a, 1x, 1(1h(), i5, 1(1h,), i3, 1(1h)), 1x, a, 1x, es16.8)')      &
            'TReigenvec_unit error: initial state vector', ilevel, 1, 'is not normalized:', f_t
            stop 'TReigenvec_unit error: initial state vector is not normalized'
         end if
         !
         do jlevel = ilevel+1,Neigenlevels
            !
            kind = eigen(jlevel)%jind
            !
            if (jind/=kind) cycle
            !
            irec = eigen(jlevel)%irec(1)
            !
            read(iounit, rec = irec) vec2(1:ncontr)
            !
            !check normalization
            !
            f_t = dot_product(vec1(:),vec2(:))
            if(abs(f_t) > 0.1_rk ** (rk - 1)) then
               write(out, '(1x, a, 1x, 1(1h(), i5, i5, i3, 1(1h)), 1x, a, 1x, es16.8)')      &
               'TReigenvec_unit error: vectoras', ilevel,jlevel, 'are not orthogonal: ', f_t
               stop 'TReigenvec_unit error: vectors are not orthogonal'
            end if
            !
         end do
         !
      end do
      !omp end do
      !
      deallocate(vec1,vec2)
      !omp end parallel
      !
      if (job%verbose>=2) write(out,"('...done!')")
      !
    end if
    !
 end function TReigenvec_unit


 subroutine TRconvert_repres_J0_to_contr(Jrot)
    !
    implicit none

    integer(ik),intent(in) :: Jrot
    integer(ik)        :: alloc,nroots,i
    integer(ik)        :: iroot,ilevel,gamma,Jval(1)

       !
       if (job%verbose>=2) write(out,"(/'Convert the J=0 eigenvec. to contracted representaion ')")
       !
       call MemoryReport
       !
       Jval(1) = jrot
       !
       if(jrot/=0) then
          write(out,"('TRconvert_repres_J0_to_contr: illegal jrot (not 0): ',i0)") jrot 
          stop 'TRconvert_repres_J0_to_contr: illegal jrot'
       end if
       !
       ! make all modes to be one class  
       if (job%convert_model_j0) then 
         PTNclasses = 1
         do i=1,FLNmodes
            job%bset(i)%class = 1
         enddo
       endif
       !
       if(PTNclasses/=1) then
          write(out,"('TRconvert_repres_J0_to_contr: illegal number of classes (not 1): ',i0)") PTNclasses 
       !   stop 'TRconvert_repres_J0_to_contr: illegal PTNclasses'
       end if
       !
       ! read contraction indexes (icase,ilambda) for all J specified
       !
       call read_contrind(1,Jval)
       !
       ! find correlation between indexes for J = 0 and J > 0
       ! We need this because the contracted matrix elements of the 
       ! dipole moment are computed on the J=0 conatracted basis functions. 
       ! When J/=0 the numbering (bookkeeping) of these basis set has changed and
       ! we need to find the correlation between the bookkeepings.  
       !
       call index_correlation(1,Jval)
       !
       ! read eigenvalues and their labeling, i.e. description;
       ! initialize file-units for reading eigenvectors
       !
       call read_eigenval(1,Jval)
       !
       job%contrfile%dscr       = 'j0'//trim(job%contrfile%dscr)
       job%contrfile%primitives = 'j0'//trim(job%contrfile%primitives)
       job%contrfile%vectors    = 'j0'//trim(job%contrfile%vectors)
       job%contrfile%dvr        = 'j0'//trim(job%contrfile%dvr)
       !
       call PTdefine_contr_from_eigenvect(Neigenroots,Neigenlevels,eigen(:))
       !
       ! 
 end subroutine TRconvert_repres_J0_to_contr

 !
 ! Compute matrix elements of the vibrational part of the kinetic energy 
 ! operator on the J=0 eigenvectors
 !

 subroutine TRconvert_matel_j0_eigen(jrot)
    implicit none
    integer(ik),intent(in) :: Jrot
    integer(ik)        :: info,imu
    logical            :: treat_vibration =.true.  ! switch off/on the vibration
    integer(ik)        :: k1,k2,ilevel,ideg,iroot,i1,i2,chkptIO,extF_rank,dimen,irec,iunit,jroot
    integer(ik)        :: idimen,imu_t,irow,ib,Nsize,mat_size,igamma,ielem,ierror
    character(len=cl)  :: job_is
    character(len=cl)  :: task
    real(rk),allocatable  :: vec(:)
    !real(rk),allocatable  :: fcoeff(:,:)
    !integer(ik),allocatable  :: cdimen(:),icoeff(:,:)
    real(rk),allocatable :: extF_me(:,:)
    character(len=20)  :: buf20
    integer(ik)        :: ncontr_t,rootsize_t,junit,iterm1=0,iterm2=1e6,islice,Nterms,iterm,icoeff
    integer(hik)       :: matsize2,matsize,rootsize,rootsize2
    real(rk),allocatable :: gmat(:,:),psi(:,:)
    real(rk),allocatable :: mat_s(:,:),mat_t(:,:)
    integer(ik),allocatable :: ijterm(:,:)
    double precision,parameter :: alpha = 1.0d0,beta=0.0d0
    character(len=cl)  :: jchar,filename
      !
      if (job%verbose>=2) write(out,"(/'Compute J=0 vib. matrix elements of the kinetic energy operator...')")
      !
      if (job%verbose>=2) call TimerStart('Convert J0-mat.elems to contr. repres.')
      !
      if(jrot/=0) then
          write(out,"('TRconvert_matel_j0_eigen: illegal jrot (not 0): ',i0)") jrot 
          stop 'TRconvert_matel_j0_eigen: illegal jrot'
      end if
      !
      if(PTNclasses/=1) then
         write(out,"('TRconvert_matel_j0_eigen: illegal number of classes (not 1): ',i0)") PTNclasses 
         stop 'TRconvert_matel_j0_eigen: illegal PTNclasses'
      end if
      !
      if (trim(job%IOkinet_action)/='CONVERT'.and.trim(job%IOextF_action)/='CONVERT'.AND..not.job%convert_model_j0) then
          write(out,"(a,a)") &
                    'TRconvert_matel_j0_eigen: Illegal MATELEM or EXTMATELEM',&
                    ' at least one must be set to CONVERT or EIGENfunc SAVE CONVERT'
          stop 'TRconvert_matel_j0_eigen: illegal MATELEM or EXTMATELEM <> CONVERT'
      end if
      !
      matsize  = int(Neigenroots*(Neigenroots+1)/2,hik)
      matsize2 = int(Neigenroots*Neigenroots,hik)
      !
      if (job%IOvector_symm) then 
        mat_size = maxval(bset_contr(1)%nsize(:))
      else
        mat_size = bset_contr(1)%Maxcontracts
      endif 
      !
      dimen = bset_contr(1)%Maxcontracts
      !
      allocate(mat_s(Neigenroots,Neigenroots),stat=info)
      call ArrayStart('mat_s',info,1,kind(mat_s),matsize2)
      !
      matsize = int(dimen*Neigenroots,hik)
      !
      if (job%verbose>=3) write(out,"(/' Allocate two matrices of ',i8,'x',i8,' = ',i0,' elements.')") & 
                          Neigenroots,Neigenroots,matsize
      !
      allocate(psi(dimen,Neigenroots),mat_t(Neigenroots,dimen),stat=info)
      call ArrayStart('psi',info,1,kind(psi),matsize)
      call ArrayStart('mat_t',info,1,kind(mat_t),matsize)
      !
      psi = 0
      !
      if (job%convert_model_j0.and.trim(job%IOkinet_action)=='SAVE') then 
          job%IOj0matel_action = 'SAVE'
      end if 
      !
      if (job%convert_model_j0.and.trim(job%IOextF_action)=='SAVE') then 
          FLextF_matelem = .true.
      end if 
      !
      if (job%IOvector_symm) then 
        !
        allocate(vec(mat_size),stat = info)
        call ArrayStart('vec',info,size(vec),kind(vec))
        !
        allocate (ijterm(bset_contr(1)%Maxsymcoeffs,sym%Nrepresen),stat=info)
        call ArrayStart('ijterm',info,size(ijterm),kind(ijterm))
        !
        do igamma = 1,sym%Nrepresen
           !
           Nterms = 0 
           !
           do iterm = 1,bset_contr(1)%Maxsymcoeffs
             !
             ijterm(iterm,igamma) = Nterms
             !
             Nterms = Nterms + bset_contr(1)%irr(igamma)%N(iterm) 
             !
           enddo
           !
        enddo
        !
        !allocate(cdimen(Neigenroots),stat=info)
        !call ArrayStart('cdimen',info,size(cdimen),kind(cdimen))
        !
        !iunit = TReigenvec_unit(1,(/0/))
        !
        ! Prepare the transformational matrix
        !
        !cdimen = 0
        !
        !if (job%verbose>=5) call TimerStart('Prepare fcoeff for J0-convertion')
        !
        !if (job%verbose>=3) write(out,"(/' Compress the eigenvectors using the thresh = ',g8.1,' ...')") job%coeff_thresh
        !
        ! Estimate the maxiaml size of the basis after compression
        !
        !cdimenmax = 0
        !
        iroot = 0
        !
        do ilevel = 1,Neigenlevels
           !
           igamma = eigen(ilevel)%igamma
           iunit = TReigenvec_unit(1,(/0/),igamma)
           !
           irec = eigen(ilevel)%irec(1)
           !
           Nsize = bset_contr(1)%nsize(igamma)
           !
           read(iunit, rec = irec) vec(1:Nsize)
           !
           do ideg = 1, eigen(ilevel)%ndeg
             !
             iroot = iroot + 1
             !
             !$omp parallel do private(icoeff,irow,ib,iterm,ielem) shared(vec) schedule(dynamic)
             do icoeff = 1,dimen
                !
                psi(icoeff,iroot) = 0 
                !
                irow = bset_contr(1)%icontr2icase(icoeff,1)
                ib   = bset_contr(1)%icontr2icase(icoeff,2)
                !
                iterm = ijterm(irow,igamma) 
                !
                do ielem = 1,bset_contr(1)%irr(igamma)%N(irow)
                   !
                   psi(icoeff,iroot) = psi(icoeff,iroot) + vec(iterm+ielem)*bset_contr(1)%irr(igamma)%repres(iterm+ielem,ideg,ib)
                   !
                enddo
                !
             enddo 
             !$omp end parallel do
             !
           end do
           !
        end do
        !
        deallocate(vec)
        call ArrayStop('vec')
        !
        deallocate(ijterm)
        call ArrayStop('ijterm')
        !
      else
        !
        !iunit = TReigenvec_unit(1,(/0/))
        !
        iroot = 0
        !
        do ilevel = 1,Neigenlevels
           !
           igamma = eigen(ilevel)%igamma
           iunit = TReigenvec_unit(1,(/0/),igamma)
           !
           do ideg = 1, eigen(ilevel)%ndeg
             !
             iroot = iroot + 1
             !
             iroot = eigen(ilevel)%iroot(ideg)
             !
             irec = eigen(ilevel)%irec(ideg)
             !
             read(iunit, rec = irec) psi(1:dimen,iroot)
             !
           enddo
           !
        enddo
        !
      endif
      !
      !if (job%verbose>=5) write(out,"(/'Maximal number of non-zero values after vector compression  = ',i,' out of ',i/)") cdimenmax,dimen
      !
      if (job%verbose>=3) write(out,"(' ...done!')") 
      !
      !if (job%verbose>=5) call TimerStop('Prepare fcoeff for J0-convertion')
      !
      if (job%verbose>=4)  call MemoryReport
      !
      ! Generate the main matelem file with the indexes stored at the top. 
      ! In case of "DIVIDE"-processing only for islice = 1 
      !
      if (trim(job%IOj0matel_action)=='SAVE') then
        !
        if (.not.job%IOmatelem_split.or.job%iswap(1)<=1) then 
          !
          job_is ='Eigen-vib. matrix elements of the rot. kinetic part'
          call IOStart(trim(job_is),chkptIO)
          !
          open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=job%kineteigen_file)
          write(chkptIO) 'Start Kinetic part'
          !
          treat_vibration = .false.
          !
          call PTstore_icontr_cnu(Neigenroots,chkptIO,job%IOj0matel_action)
          !
          if (job%vib_rot_contr) then
            write(chkptIO) 'vib-rot'
          endif
          !
        endif
        !
        rootsize = int(bset_contr(1)%Maxcontracts*(bset_contr(1)%Maxcontracts+1)/2,hik)
        rootsize2= int(bset_contr(1)%Maxcontracts*bset_contr(1)%Maxcontracts,hik)
        rootsize2= int(bset_contr(1)%Maxcontracts,hik)
        rootsize2 = rootsize2*rootsize2
        !
        allocate(gmat(dimen,dimen),stat=info)
        call ArrayStart('gmat-fields',info,1,kind(gmat),rootsize2)
        !
        ! Preparing slicing 
        !
        !
        if (job%IOmatelem_split) then
          !
          iterm1 = max(job%iswap(1),1)
          iterm2 = min(job%iswap(2),(FLNmodes+3)*3)
          !
          if (job%verbose>=4) write(out,"('  The j0_matelem.chk will be divided into 3 x 3 + ',i4,'x 3  = ',i5,' chk-slices')") &
                                    FLNmodes,9+3*FLNmodes
          if (job%verbose>=4) write(out,"('  islice = 1-9 (Grot), and 10-',i4,'(Gcor). ')") 9+3*FLNmodes
          if (job%verbose>=4) write(out,"('  This run is for the checkpoint slices from ',i4,' to ',i4)") iterm1,iterm2
          !
        endif
        !

        ! The eigen-vibrational (J=0) matrix elements of the rotational and coriolis 
        ! kinetic parts are being computed here. 
        !
        if (job%verbose>=3) write(out,"(/' Transform grot to J0-repres...')")
        !
        if (.not.job%IOmatelem_split) then
          !
          task = 'rot'
          !
          write(chkptIO) 'g_rot'
          !
          call restore_rot_kinetic_matrix_elements(jrot,treat_vibration,task,iunit)
          !
        else
          !
          task = 'top'
          !
          call restore_rot_kinetic_matrix_elements(jrot,treat_vibration,task,iunit)
          !
        endif
        !
        if (job%verbose>=5) call TimerStart('J0-convertion for g_rot')
        !
        ! Run the loop over all term of the expansion of the Hamiltonian 
        !
        islice = 0
        !
        do k1 = 1,3
          !
          do k2 = 1,3
            !
            islice = islice + 1
            !
            if (job%IOmatelem_split.and.(islice<iterm1.or.iterm2<islice)) cycle
            !
            if (job%IOmatelem_split.and..not.job%vib_rot_contr) then 
              !
              call divided_slice_read(islice,'g_rot',job%matelem_suffix,dimen,gmat,ierror)
              !
            elseif (job%IOmatelem_split.and.job%vib_rot_contr) then 
              !
              call divided_slice_read_vibrot(islice,job%matelem_suffix,dimen,gmat)
              !
            else
              !
              read(iunit) gmat
              !
            endif
            !
            !
            call dgemm('T','N',Neigenroots,dimen,dimen,alpha,psi,dimen,& 
                        gmat,dimen,beta,mat_t,Neigenroots)
            call dgemm('N','N',Neigenroots,Neigenroots,dimen,alpha,mat_t,Neigenroots,& 
                        psi,dimen,beta,mat_s,Neigenroots)
            !
            if (job%IOmatelem_split.and..not.job%vib_rot_contr) then 
              !
              call divided_slice_write(islice,'g_rot',job%j0matelem_suffix,Neigenroots,mat_s)
              !
            elseif (job%IOmatelem_split.and.job%vib_rot_contr) then 
              !
              call divided_slice_write_vibrot(islice,job%j0matelem_suffix,Neigenroots,mat_s)
              !
            else
              !
              write (chkptIO) mat_s
              !
            endif
            !
          enddo
          !
        enddo
        !
        if (job%verbose>=5) call TimerStop('J0-convertion for g_rot')
        !
        if (job%verbose>=3) write(out,"(' ...done!')")
        if (job%verbose>=3) write(out,"(/' Transform gcor to J0-repres...')")
        !
        if (.not.job%IOmatelem_split) then
          !
          task = 'cor'
          !
          call restore_rot_kinetic_matrix_elements(jrot,treat_vibration,task,iunit)
          !
          write(chkptIO) 'g_cor'
          !
        endif
        !
        if (job%verbose>=5) call TimerStart('J0-convertion for g_cor')
        !
        ! Run the loop over all term of the expansion of the Hamiltonian 
        !
        islice = 9
        !
        !do k1 = 1,FLNmodes
          !
          !if (job%contrci_me_fast.and.k1>1) cycle
          !
          do k2 = 1,3
            !
            islice = islice + 1
            !
            if (job%IOmatelem_split.and.(islice<iterm1.or.iterm2<islice)) cycle
            !
            if (job%IOmatelem_split.and..not.job%vib_rot_contr) then 
              !
              call divided_slice_read(islice,'g_cor',job%matelem_suffix,dimen,gmat,ierror)
              !
            elseif (job%IOmatelem_split.and.job%vib_rot_contr) then 
              !
              call divided_slice_read_vibrot(islice,job%matelem_suffix,dimen,gmat)
              !
            else
              !
              read(iunit) gmat
              !
            endif
            !
            call dgemm('T','N',Neigenroots,dimen,dimen,alpha,psi,dimen,& 
                        gmat,dimen,beta,mat_t,Neigenroots)
            call dgemm('N','N',Neigenroots,Neigenroots,dimen,alpha,mat_t,Neigenroots,& 
                        psi,dimen,beta,mat_s,Neigenroots)
            !
            !
            if (job%IOmatelem_split.and..not.job%vib_rot_contr) then 
              !
              call divided_slice_write(islice,'g_cor',job%j0matelem_suffix,Neigenroots,mat_s)
              !
            elseif (job%IOmatelem_split.and.job%vib_rot_contr) then 
              !
              call divided_slice_write_vibrot(islice,job%j0matelem_suffix,Neigenroots,mat_s)
              !
            else
              !
              write (chkptIO) mat_s
              !
            endif
            !
          enddo
          ! 
        !enddo
        !
        if (job%verbose>=5) call TimerStop('J0-convertion for g_cor')
        !
        if (.not.job%IOmatelem_split.or.job%iswap(1)==1) write(chkptIO) 'End Kinetic part'
        !
        if (.not.job%vib_rot_contr) then 
          close(chkptIO,status='keep')
        endif
        !
        task = 'end'
        !
        call restore_rot_kinetic_matrix_elements(jrot,treat_vibration,task,iunit)
        !
        if (allocated(gmat)) deallocate(gmat)
        call ArrayStop('gmat-fields')
        !
        if (job%verbose>=3) write(out,"(' ...done!')")
        !
      endif 
      !
      ! External field part 
      !
      if (FLextF_matelem) then
        !
        if (job%verbose>=3) write(out,"(/' Transform extF to J0-representation...')")
        !
        if (job%verbose>=2) write(out,"(/' External function...')")
        !
        extF_rank = FLread_extF_rank()
        !
        if (job%IOextF_divide) then
          !
          iterm1 = max(fitting%iparam(1),1)
          iterm2 = min(fitting%iparam(2),extF_rank)
          !
          if (job%verbose>=4) write(out,"('  The j0_external.chk will be divided into ',i6,' chk-slices')") extF_rank
          if (job%verbose>=4) write(out,"('  This run is for the checkpoint slices from ',i4,' to ',i4)") iterm1,iterm2
          !
          ncontr_t = bset_contr(1)%Maxcontracts
          !
        else
          !
          job_is ='external field contracted matrix elements for J=0'
          call IOStart(trim(job_is),iunit)
          !
          filename = job%extFmat_file
          !
          open(iunit,form='unformatted',action='read',position='rewind',status='old',file=filename)
          !
          read(iunit) buf20
          if (buf20/='Start external field') then
            write (out,"(' restore_vib_matrix_elements ',a,' has bogus header: ',a)") filename,buf20
            stop 'restore_vib_matrix_elements - bogus file format'
          end if
          !
          read(iunit) ncontr_t
          !
          if (bset_contr(1)%Maxcontracts/=ncontr_t) then
            write (out,"(' Dipole moment checkpoint file ',a)") filename
            write (out,"(' Actual and stored basis sizes at J=0 do not agree  ',2i8)") bset_contr(1)%Maxcontracts,ncontr_t
            stop 'restore_Extvib_matrix_elements - in file - illegal ncontracts '
          end if
          !
        endif
        !
        if (.not.job%IOextF_divide.or.job%IOextF_stitch) then 
          !
          ! Prepare the checkpoint file
          !
          job_is ='external field contracted matrix elements for J=0'
          call IOStart(trim(job_is),chkptIO)
          !
          open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=job%exteigen_file)
          write(chkptIO) 'Start external field'
          !
          ! store the matrix elements 
          !
          write(chkptIO) Neigenroots
          !
        endif
        !
        !if (job%IOextF_divide) close(iunit)
        !
        rootsize = int(ncontr_t*(ncontr_t+1)/2,hik)
        rootsize2= int(ncontr_t*ncontr_t,hik)
        !
        if (job%verbose>=4) write(out,"(/'restore Extvib...: Number of elements: ',i0)") ncontr_t
        !
        allocate(extF_me(ncontr_t,ncontr_t),stat=info)
        call ArrayStart('extF_me',info,1,kind(extF_me),rootsize2)
        !
        do imu = fitting%iparam(1),fitting%iparam(2)
          !
          if (job%verbose>=4) write(out,"('  imu = ',i8)",advance='NO') imu
          !
          if (job%IOextF_divide) then
            !
            call divided_slice_read(imu,'extF',job%extmat_suffix,dimen,extF_me,ierror)
            !
            if (ierror==1) cycle
            !
            if (job%verbose>=4) write(out,"('.')",advance='YES')
            !
          else
            !
            read(iunit) imu_t
            !
            read(iunit) extF_me
            !
          endif
          !
          call dgemm('T','N',Neigenroots,dimen,dimen,alpha,psi,dimen,& 
                      extF_me,dimen,beta,mat_t,Neigenroots)
          call dgemm('N','N',Neigenroots,Neigenroots,dimen,alpha,mat_t,Neigenroots,& 
                      psi,dimen,beta,mat_s,Neigenroots)
          !
          !mat_s = 0 
          !
          !do ilevel = 1, Neigenlevels
          !  !
          !  do ideg = 1,eigen(ilevel)%ndeg
          !    !
          !    iroot = eigen(ilevel)%iroot(ideg)
          !    !
          !    !i1 = iroot*(iroot-1)/2+1
          !    !i2 = i1-1+iroot
          !    !
          !    call eigen_vib_matelem_vector(0,ilevel,iroot,Neigenlevels,Neigenroots,cdimenmax,icoeff,fcoeff,cdimen,extF_me,mat_s(iroot,1:iroot))
          !    !
          !  enddo
          !  !
          !enddo
          !
          !!$omp parallel do private(iroot,jroot) schedule(dynamic) shared(mat_s)
          !do iroot=1,Neigenroots
          !  do jroot=1,iroot-1
          !    mat_s(jroot,iroot) = mat_s(iroot,jroot)
          !  enddo
          !enddo
          !!$omp end parallel do
          !
          !
          if (.not.job%IOextF_divide.or.job%IOextF_stitch) then
            !
            write(chkptIO) imu
            write(chkptIO) mat_s
            !
          else
            !
            call divided_slice_write(imu,'extF',job%j0extmat_suffix,Neigenroots,mat_s)
            !
          endif
          !
        enddo
        !
        if (allocated(extF_me)) deallocate(extF_me)
        call ArrayStop('extF_me')
        !
        if (.not.job%IOextF_divide) then
          !
          read(iunit) buf20(1:18)
          if (buf20(1:18)/='End external field') then
            write (out,"(' restore_Extvib_matrix_elements ',a,' has bogus footer: ',a)") job%kinetmat_file,buf20(1:17)
            stop 'restore_Extvib_matrix_elements - bogus file format'
          end if
          !
          close(iunit,status='keep')
          !
          !job_is ='external field contracted matrix elements for J=0'
          !call IOStart(trim(job_is),iunit)
          !
        endif
        !
        if (.not.job%IOextF_divide.or.job%IOextF_stitch) then
          !
          write(chkptIO) 'End external field'
          close(chkptIO,status='keep')
          !
        endif
        !
        if (job%verbose>=3) write(out,"(/' ...done!')")
        !
      endif 
      !
      deallocate(mat_s)
      call ArrayStop('mat_s')
      !
      deallocate(mat_t)
      call ArrayStop('mat_t')
      !
      deallocate(psi)
      call ArrayStop('psi')
      !
      !deallocate(icoeff,fcoeff)
      !call ArrayStop('fcoeff')
      !
      call MemoryReport
      !
      if (job%verbose>=2) call TimerStop('Convert J0-mat.elems to contr. repres.')
      !
      call TimerReport
      !
      if (job%verbose>=2) write(out,"('...done!')")   
      !
    contains
      !
      subroutine divided_slice_write(islice,name,suffix,N,field)
        !
        implicit none
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
      subroutine divided_slice_read(islice,name,suffix,N,field,ierror)
        !
        implicit none
        !
        integer(ik),intent(in)      :: islice
        character(len=*),intent(in) :: name,suffix
        integer(ik),intent(in)      :: N
        real(rk),intent(out)        :: field(N,N)
        integer(ik),intent(out)     :: ierror
        character(len=4)            :: jchar
        integer(ik)                 :: chkptIO
        character(len=cl)           :: buf,filename,job_is
        integer(ik)                 :: ilen
        logical                     :: ifopened
        !
        ierror = 0
        !
        write(job_is,"('single swap_matrix')")
        !
        call IOStart(trim(job_is),chkptIO)
        !
        write(jchar, '(i4)') islice
        !
        filename = trim(suffix)//trim(adjustl(jchar))//'.chk'
        !
        open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=filename,err=15)
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
        return
        !
        ! This error code will allow simply skipping the corresponding record/file without crashing the program 
        !
   15   ierror = 1
        !
        ! we allow to skip opening the file only for the external matrix elements
        !
        if (trim(name)/="extF") then 
          write (out,"(' kinetic divided_slice_read in slice ',a20,': file does not exist')") filename
          stop 'divided_slice_read - in slice -  file does not exist'
        endif
        !
        if (job%verbose>=4) write (out,"(' (skipped).')",advance='YES') 
        !
      end subroutine divided_slice_read
      !
      subroutine divided_slice_read_vibrot(islice,suffix,N,field)
        !
        implicit none
        !
        integer(ik),intent(in)      :: islice
        character(len=*),intent(in) :: suffix
        integer(ik),intent(in)      :: N
        real(rk),intent(out)        :: field(N,N)
        character(len=4)            :: jchar
        integer(ik)                 :: chkptIO
        character(len=cl)           :: buf,filename,job_is
        integer(ik)                 :: rec_len,icontr,icontr1,icontr2
        logical                     :: ifopened
        real(rk)                    :: f_t
        !
        if (.not.job%IOmatelem_split) return
        !
        write(jchar, '(i4)') islice
        !
        filename = trim(suffix)//trim(adjustl(jchar))//'.chk'
        !
        write(job_is,"('single swap_matrix #',i8)") islice
        !
        call IOStart(trim(job_is),chkptIO)
        !
        inquire(iolength=rec_len) f_t
        rec_len = rec_len*N*bset_contr(1)%max_deg_size
        !
        ! use direct access format for vib-rot
        !
        open(chkptIO,access='direct',action='read',status='old',recl=rec_len,file=filename)
        !
        do icontr=1,N
          !
          icontr1 = Ncontr02icase0(icontr,1)
          icontr2 = Ncontr02icase0(icontr,2)
          !
          read(chkptIO,rec=icontr) field(1:N,icontr1:icontr2)
          !
        enddo
        !
        close(chkptIO)
        call IOStop(trim(job_is))
        !
      end subroutine divided_slice_read_vibrot
      !
      !
      subroutine divided_slice_write_vibrot(islice,suffix,N,field)
        !
        implicit none
        !
        integer(ik),intent(in)      :: islice
        character(len=*),intent(in) :: suffix
        integer(ik),intent(in)      :: N
        real(rk),intent(in)         :: field(N,N)
        character(len=4)            :: jchar
        integer(ik)                 :: chkptIO
        character(len=cl)           :: buf,filename,job_is
        integer(ik)                 :: rec_len,icontr,icontr1,icontr2
        logical                     :: ifopened
        real(rk)                    :: f_t
        !
        if (.not.job%IOmatelem_split) return
        !
        write(jchar, '(i4)') islice
        !
        filename = trim(suffix)//trim(adjustl(jchar))//'.chk'
        !
        write(job_is,"('single swap_j0_matrix #',i8)") islice
        !
        call IOStart(trim(job_is),chkptIO)
        !
        inquire(iolength=rec_len) f_t
        rec_len = rec_len*N*bset_contr(1)%max_deg_size
        !
        ! use direct access format for vib-rot
        !
        open(chkptIO,access='direct',action='write',status='replace',recl=rec_len,file=filename)
        !
        do icontr=1,N
          !
          icontr1 = Ncontr02icase0(icontr,1)
          icontr2 = Ncontr02icase0(icontr,2)
          !
          write(chkptIO,rec=icontr) field(1:N,icontr1:icontr2)
          !
        enddo
        !
        close(chkptIO,status='keep')
        call IOStop(trim(job_is))
        !
      end subroutine divided_slice_write_vibrot

      !
 end subroutine TRconvert_matel_j0_eigen
 !

 ! Here we restore the vibrational (J=0) matrix elements of the rotational kinetic part G_rot and G_cor
 !
 subroutine restore_rot_kinetic_matrix_elements(jrot,treat_vibration,kinetic_part,chkptIO)
   !
   implicit none
   !
   integer(ik),intent(in)       :: jrot
   logical,intent(in)           :: treat_vibration
   character(len=3),intent(in)  :: kinetic_part
   integer(ik),intent(inout)    :: chkptIO
   !
   integer(ik)        :: iclasses,alloc,ilevel,ib,max_deg_size,nclasses,islice
   character(len=cl)  :: job_is
   character(len=18)  :: buf18
   integer(ik)        :: ncontr_t,rootsize
   integer(ik),allocatable :: imat_t(:,:)
   real(rk),allocatable :: mat_t(:)
   !
   nclasses = bset_contr(1)%nclasses
   !
   select case (kinetic_part)
   !
   case('rot','top')
     !
     job_is ='Vib. matrix elements of the rot. kinetic'
     call IOStart(trim(job_is),chkptIO)
     !
     open(chkptIO,form='unformatted',action='read',position='rewind',status='old',file=job%kinetmat_file)
     !
     read(chkptIO) buf18
     if (buf18/='Start Kinetic part') then
       write (out,"(' Vib. kinetic checkpoint file ',a,' has bogus header: ',a)") job%kinetmat_file,buf18
       stop 'PTcontracted_matelem_class - bogus file format'
     end if
     !
     read(chkptIO) ncontr_t
     !
     if (bset_contr(1)%Maxcontracts/=ncontr_t) then
       write (out,"(' Vib. kinetic checkpoint file ',a)") job%kinetmat_file
       write (out,"(' Actual and stored basis sizes at J=0 do not agree  ',2i0)") bset_contr(1)%Maxcontracts,ncontr_t
       stop 'PTcontracted_matelem_class - in file - illegal nroots '
     end if
     !
     rootsize = bset_contr(1)%Maxcontracts*(bset_contr(1)%Maxcontracts+1)/2
     !
     if (job%verbose>=6) write(out,"(/'Restore_rot_kin...: Number of elements: ',i8)") bset_contr(1)%Maxcontracts
     !
     ! Read the indexes of the J=0 contracted basis set. 
     !
     allocate (imat_t(0:nclasses,ncontr_t),stat=alloc)
     call ArrayStart('contractive_space',alloc,size(imat_t),kind(imat_t))
     !
     read(chkptIO) buf18(1:10)
     if (buf18(1:10)/='icontr_cnu') then
       write (out,"(' Vib. kinetic checkpoint file ',a,': icontr_cnu is missing ',a)") job%kinetmat_file,buf18(1:10)
       stop 'restore_rot_kinetic_matrix_elements - in file -  icontr_cnu missing'
     end if
     !
     read(chkptIO) imat_t(0:nclasses,1:ncontr_t)
     !
     if (job%vib_rot_contr) then
        !
        allocate (TRicontr_cnu(0:nclasses,ncontr_t),stat=alloc)
        call ArrayStart('TRicontr_cnu',alloc,size(TRicontr_cnu),kind(TRicontr_cnu))
        TRicontr_cnu = imat_t
        !
     endif
     !
     read(chkptIO) buf18(1:11)
     if (buf18(1:11)/='icontr_ideg') then
       write (out,"(' Vib. kinetic checkpoint file ',a,': icontr_ideg is missing ',a)") job%kinetmat_file,buf18(1:11)
       stop 'restore_rot_kinetic_matrix_elements - in file -  icontr_ideg missing'
     end if
     !
     read(chkptIO) imat_t(0:nclasses,1:ncontr_t)
     !
     deallocate(imat_t)
     !
     call arraystop('contractive_space')
     !
     ! Read the rotational part only if its really needed. 
     !
     if (trim(kinetic_part)=='rot') then
       !
       read(chkptIO) buf18(1:4)
       if (buf18(1:4)/='g_ro') then
         write (out,"(' Vib. kinetic checkpoint file ',a,': g_rot is missing ',a)") trim(job%kinetmat_file),buf18(1:5)
         !
         if (buf18(1:4)=='hvib'.or.buf18(1:3)=='End') &
               write (out,"(a,a)") &
               ' Most likely the split chk-points are supposed to be used.',&
               'Re-do MATELEM SAVE or use SPLIT in MATELEM !'
         stop 'restore_rot_kinetic_matrix_elements - in file -  g_rot missing'
       end if
       !
     endif
     !
     if (job%vib_rot_contr) then
       !
       !inquire(iolength=rec_len) f_t
       !rec_len = rec_len*ncontr*PT%max_deg_size
       !!
       !do islice = 1,9+3*PT%Nmodes
       !    call divided_slice_open_vib_rot(islice,rec_len,job%matelem_suffix)
       !enddo
       !
       call find_groundstate_icontr(bset_contr(1)%nclasses)
       !
     endif
     !
   case('cor')
     !
     read(chkptIO) buf18(1:5)
     if (buf18(1:5)/='g_cor') then
       write (out,"(' Vib. kinetic checkpoint file ',a,': g_cor is missing ',a)") job%kinetmat_file,buf18(1:5)
       stop 'restore_rot_kinetic_matrix_elements - in file -  g_cor missing'
     end if
     !
   case('end')
     !
     if (job%vib_rot_contr) then
        !
        do islice = 1,9+3*FLNmodes
            call divided_slice_close_vib_rot(islice)
        enddo
        !
     endif
     !
     read(chkptIO) buf18(1:4)
     if (buf18(1:4)/='hvib'.and.buf18(1:3)/='End'.and.buf18(1:4)/='vib-') then
       write (out,"(' Vib. kinetic checkpoint file ',a,': hvib, vib-rot, End is missing ',a)") job%kinetmat_file,buf18(1:4)
       stop 'restore_rot_kinetic_matrix_elements - in file -  hvib or End missing'
     end if
     !
     close(chkptIO,status='keep')
     !
   end select
   !
   contains
    !
    subroutine divided_slice_open_vib_rot(islice,rec_len,suffix)
      !
      implicit none
      integer(ik),intent(in)      :: islice,rec_len
      character(len=*),intent(in) :: suffix
      integer(ik)                 :: chkptIO
      character(len=4)            :: jchar
      character(len=cl)           :: buf,filename,job_is
      integer(ik)                 :: ilen
      logical                     :: ifopened
      !
      if (.not.job%IOmatelem_split) return
      !
      write(jchar, '(i4)') islice
      !
      filename = trim(suffix)//trim(adjustl(jchar))//'.chk'
      !
      write(job_is,"('single swap_matrix #',i8)") islice
      !
      call IOStart(trim(job_is),chkptIO)
      !
      ! use direct access format for vib-rot
      !
      open(chkptIO,access='direct',action='read',status='old',recl=rec_len,file=filename)
      !
    end subroutine divided_slice_open_vib_rot
    !
    subroutine divided_slice_close_vib_rot(islice)
      !
      implicit none
      integer(ik),intent(in)      :: islice
      integer(ik)                 :: chkptIO
      character(len=4)            :: jchar
      character(len=cl)           :: filename,job_is
      integer(ik)                 :: ilen
      !
      if (.not.job%IOmatelem_split) return
      !
      write(jchar, '(i4)') islice
      !
      write(job_is,"('single swap_matrix #',i8)") islice
      !
      call IOStart(trim(job_is),chkptIO)
      !
      ! use direct access format for vib-rot
      !
      close(chkptIO,status='keep')
      !
    end subroutine divided_slice_close_vib_rot    
    !
    ! find correspondence between contracted quantum numbers: current and for J=0 
    !
    subroutine find_groundstate_icontr(Nclasses)
       !
       integer(ik), intent(in) :: Nclasses
       integer(ik) :: maxcontr
       !
       integer(ik)  :: icontr,iterm,cnu_(1:bset_contr(1)%nclasses),icontr0,alloc,maxcontr0
       !
       if (job%verbose>=4) write(out,"('   Find correlation between J=0 and J/=0 contr. basis functions...')")
       !
       maxcontr = size(TRicontr_cnu,dim=2)
       !
       ! count icontr0-s
       icontr0 = 0
       cnu_ = 0
       !
       do icontr = 1,maxcontr
         !
         if (any(TRicontr_cnu(1:Nclasses,icontr)/=cnu_(1:Nclasses))) then 
           !
           icontr0 = icontr0 + 1
           cnu_ = TRicontr_cnu(1:Nclasses,icontr)
           !
         endif
         !
       enddo
       !
       maxcontr0 = icontr0
       !
       allocate(Ncontr02icase0(maxcontr0,2),stat=alloc)
       call ArrayStart('Ncontr02icase0',alloc,size(Ncontr02icase0),kind(Ncontr02icase0))
       !
       ! count ideg0
       icontr0 = 0
       cnu_ = 0
       !
       Ncontr02icase0(1,1) = 1
       Ncontr02icase0(maxcontr0,2) = maxcontr
       !
       do icontr = 1,maxcontr
         !
         if (any(TRicontr_cnu(1:Nclasses,icontr)/=cnu_(1:Nclasses))) then 
           !
           icontr0 = icontr0 + 1
           cnu_ = TRicontr_cnu(1:Nclasses,icontr)
           Ncontr02icase0(icontr0,1) = icontr
           if (icontr0>1) Ncontr02icase0(icontr0-1,2) = icontr-1
           !
         endif
         !
       enddo
       !
    end subroutine find_groundstate_icontr
    !
 end subroutine restore_rot_kinetic_matrix_elements
 !
 !
 !
 subroutine eigen_vib_matelem_vector(iparity,ilevelI,irootI,nlevels,nroots,cdimenmax,icoeff,fcoeff,cdimen,field,mat)
   !
   implicit none
   !
   integer(ik),intent(in)  :: iparity,ilevelI,irootI,nlevels,nroots,cdimenmax
   integer(ik),intent(in)  :: icoeff(cdimenmax,nroots)
   real(rk),intent(in)     :: fcoeff(cdimenmax,nroots)
   integer(ik),intent(in)  :: cdimen(nroots)
   real(rk),intent(in)     :: field(:,:)
   real(rk),intent(out)    :: mat(:)
   !
   integer(ik)             :: ilevelF,idegF,info,dimen,cdimen_,M,MP1,idimen
   real(rk),allocatable    :: half_matelem(:)
   integer(ik),allocatable :: contr(:)
   integer(ik)             :: i,indJ,irootF, cirootI, icase, ilambda, icontrF, icontrI
   real(rk)                :: f_t,dtemp
    !
    !if (job%verbose>=5) call TimerStart('Eigen-vib-me-vec')
    !
    ! Compute the half_matelem
    !
    dimen = bset_contr(1)%Maxcontracts
    !
    allocate(half_matelem(dimen),stat=info)
    call ArrayStart('half_matelem',info,size(half_matelem),kind(half_matelem))
    !
    indJ = 1
    !
    !loop over final state basis components
    !
    if (job%verbose>=5) call TimerStart('J0-convertion: half_matelem')
    !
    if (job%verbose>=7) write(out,"('  J0-convertion: half_matelem ...')")
    !
    !half_matelem = 0
    !
    cdimen_ = cdimen(irootI)
    !
    ! precompute the indexing
    !
    allocate(contr(cdimen_),stat=info)
    call ArrayStart('half_matelem',info,size(contr),kind(contr))
    !
    !$omp parallel do private(cirootI) schedule(dynamic) shared(contr)
    do cirootI = 1,cdimen_
      contr(cirootI) = bset_contr(indJ)%iroot_correlat_j0(icoeff(cirootI,irootI))
    enddo
    !$omp end parallel do
    !
    if (iparity==1) then 
      !
      !$omp parallel do private(irootF,icontrF,dtemp,cirootI,icontrI,f_t) schedule(dynamic) ! shared(half_matelem)
      do irootF = 1, dimen
         !
         icontrF = bset_contr(indJ)%iroot_correlat_j0(irootF)
         !
         dtemp = 0
         !
         !loop over initial state basis components
         !
         do cirootI = 1,cdimen_
            !
            icontrI = contr(cirootI)
            !
            !icontrI = bset_contr(indJ)%iroot_correlat_j0(icoeff(cirootI,irootI))
            !
            !index of the corresponding vibrational contracted matrix element (cind)
            !
            !irow = max(icontrF, icontrI)
            !icol = min(icontrF, icontrI)
            !cind = ( irow * (irow - 1) )/ 2 + icol
            !
            !cind = i2d_to_1d(icontrI,icontrF)
            !
            !f_t = sign(1,icontrI-icontrF) 
            !
            f_t = fcoeff(cirootI,irootI)*field(icontrI,icontrF)
            !
            if (icontrI<icontrF) f_t = -f_t
            !
            dtemp = dtemp + f_t ! *fcoeff(cirootI,irootI)*field(cind)
            !
         end do
         !
         half_matelem(irootF) = dtemp
         !
      end do
      !$omp end parallel do
      !
    else
      !
      !$omp parallel do private(irootF,icontrF,dtemp,cirootI,icontrI) schedule(dynamic) ! shared(half_matelem) 
      do irootF = 1, dimen
         !
         icontrF = bset_contr(indJ)%iroot_correlat_j0(irootF)
         !
         dtemp = 0
         !
         !loop over initial state basis components
         do cirootI = 1, cdimen_
            !
            icontrI = contr(cirootI)
            !
            !icontrI = bset_contr(indJ)%iroot_correlat_j0(icoeff(cirootI,irootI))
            !
            !index of the corresponding vibrational contracted matrix element (cind)
            !
            !irow = max(icontrF, icontrI)
            !icol = min(icontrF, icontrI)
            !cind = ( irow * (irow - 1) )/ 2 + icol
            !
            !cind = i2d_to_1d(icontrI,icontrF)
            !
            dtemp = dtemp + field(icontrI,icontrF)*fcoeff(cirootI,irootI)
            !
         end do
         !
         half_matelem(irootF) = dtemp
         !
      end do
      !$omp end parallel do
      !
    endif 
    !
    if (job%verbose>=7) write(out,"('  ...done!')")
    if (job%verbose>=5) call TimerStop('J0-convertion: half_matelem')
    !
    if (job%verbose>=5) call TimerStart('J0-convertion: step 2')
    !
    if (job%verbose>=7) write(out,"('  J0-convertion: step2 ...')")
    !
    !loop over final states
    !
    !ilevel = eigen(irootI)%ilevel
    !
    !$omp parallel do private(ilevelF,idegF,irootF,cdimen_,m,dtemp,idimen,mp1) shared(mat) schedule(dynamic)
    Flevels_loop: do ilevelF = 1,ilevelI
       !
       do idegF = 1,eigen(ilevelF)%ndeg
         !
         irootF = eigen(ilevelF)%iroot(idegF)
         !
         if (irootF>irootI) cycle
         !
         ! complete the me 
         !
         cdimen_ = cdimen(irootF)
         !
         M = mod(cdimen_,5)
         !
         dtemp = 0
         !
         if (M/=0) then 
           do idimen = 1,M
             !
             dtemp = dtemp + half_matelem( icoeff( idimen,irootF ) )*fcoeff(idimen,irootF)
             !
           enddo
         endif
         !
         MP1 = M + 1
         !
         do idimen = MP1,cdimen_,5
           !
           dtemp = dtemp + & 
                           half_matelem( icoeff( idimen  ,irootF ) )*fcoeff(idimen  ,irootF) + &
                           half_matelem( icoeff( idimen+1,irootF ) )*fcoeff(idimen+1,irootF) + &
                           half_matelem( icoeff( idimen+2,irootF ) )*fcoeff(idimen+2,irootF) + &
                           half_matelem( icoeff( idimen+3,irootF ) )*fcoeff(idimen+3,irootF) + &
                           half_matelem( icoeff( idimen+4,irootF ) )*fcoeff(idimen+4,irootF)
           !
         enddo
         !
         !
         mat(irootF) = mat(irootF) + dtemp
         !
         ! mat(irootF) = mat(irootF) + sum( half_matelem(tmat(irootF)%icoeff(1:cdimen(irootF)) )*tmat(irootF)%coeff(1:cdimen(irootF)))
         !
       end do
       !
    end do Flevels_loop
    !$omp end parallel do
    !
    if (job%verbose>=7) write(out,"('  ...done!')")
    !
    if (job%verbose>=5) call TimerStop('J0-convertion: step 2')
    !
    deallocate(half_matelem,contr)
    call ArrayStop('half_matelem')
    !
    !if (job%verbose>=5) call TimerStop('Eigen-vib-me-vec')
    !
 end subroutine eigen_vib_matelem_vector




function to_upper(str) result (string)

  character(*), intent(in) :: str
  character(len(str)) :: string

  integer :: ic, i

  character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

  string = str
  do i = 1, len_trim(str)
    ic = index(low, str(i:i))
    if (ic>0) string(i:i) = cap(ic:ic)
  end do

end function to_upper


function to_lower(str) result (string)

  character(*), intent(in) :: str
  character(len(str)) :: string

  integer :: ic, i

  character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

  string = str
  do i = 1, len_trim(str)
    ic = index(cap, str(i:i))
    if (ic>0) string(i:i) = low(ic:ic)
  end do

end function to_lower


subroutine tokenize(str, n, w)

  character(*), intent(in) :: str
  integer(ik), intent(out) :: n
  character(cl), intent(out) :: w(:)

  integer(ik) :: pos2, pos1, maxn, ln
  character(1000) :: s

  maxn = size(w)
  n = 0
  pos1 = 1

  DO
    s = trim(adjustl(str))
    ln = len(trim(adjustl(str)))
    if (pos1>ln) exit
    pos2 = INDEX(s(pos1:ln), ' ')
    IF (pos2 == 0) THEN
       n = n + 1
       if (n>maxn) then
         write(out, '(/a,1x,i4)') 'timer/tokenize error: number of words exceeds maximum =', maxn
         stop 'STOP'
       endif
       w(n) = s(pos1:ln)
       EXIT
    END IF
    if (s(pos1:pos1+pos2-1)=='') then
      pos1 = pos2+pos1
      cycle
    endif
    n = n + 1
    if (n>maxn) then
      write(out, '(/a,1x,i4)') 'timer/tokenize error: number of words exceeds maximum =', maxn
      stop 'STOP'
    endif
    w(n) = s(pos1:pos1+pos2-1)
    pos1 = pos2+pos1
 END DO

end subroutine tokenize


end module tran

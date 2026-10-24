 module FITPES_C2H6
  !
  public ik,rk,do_fitting,pi

  private
 

 integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
 integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
 integer, parameter :: ark         = selected_real_kind(25,32)  ! "Normal" reals and complex (complexi? :-)
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10
!, nparam_eq = 5
!  TAKING THIS OUT FOR NOW TO TRY AND MAKE MORE GENERAL (BPM)
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  real(rk)  :: pi


  contains
  
 subroutine do_fitting
 !
 implicit none 
 !
 ! Parameters:  
 ! 
 ! number of parameters and maximum number of ab initio points   
 !
 !
 integer,parameter          ::  enermax =  300000, nmax0=8,taumax0=8, nparam = 1000
 !
 ! where the output goes to 
 ! 
 integer,parameter          :: f_inp=5, f_out=6, f_res=9, f_new =10
 ! , nparam_eq = 5
 ! 
 ! Dinite difference differentiation with 2 or 3 points
 !
 integer,parameter          ::  findif = 3 
 !
 ! the iteration will be stoped when the best standard deviation reached 
 !
 !double precision,parameter :: stadev_best = 0.001d0   
 !
 integer,parameter          ::  verbose = 5
 ! 
 ! parameter to control the correlation when we through out the parameters with large st. error (0.0..1.0)
 !
 double precision,parameter :: fact_corr=0.6, too_big = 1e8
 ! 
 ! Watson's alpha parameter to control the Robust fit
 double precision,parameter :: a_wats = 0.1
 ! parameter for the finite differncies differention 
 !
 double precision,parameter :: factordeltax=0.001
 !
 integer,parameter          ::  ireplace_by = 0
 !
 ! Universal constants 
 ! 
 double precision,parameter :: bohr = 0.52917715e+00,hartre=2.194746354d+05
 ! Avareging region for the robustfit = 3 
 integer,parameter          :: ishift = 10
  character(len=40),parameter  :: fit_type= 'linur'  !'dgelss'
 ! Variables:  
 double precision  :: maxdeviat  ! st. deviation at the previos iteration 
 !
 integer           :: i,npts,nused,itmax,numpar,ncol,iter,nrow,l,itertypemax,irow,icolumn ! just integers  
 ! 
 character(len=80) :: title(4),longlabel ! input title of the job from the first four input lines 
 character(len=1)  :: mark           ! used at the parameters output to mark the dependent or independent ones 
 character(len=10) :: label          ! Temporary label 
 double precision  :: wtsum,ssq      ! 
 !
 ! some variables 
 !
 double precision  :: stadev,epsil,stadev_old,v,dd
 double precision  :: tempx,deltax,potright,potleft,stability
 double precision  :: r14,r24,r34,alpha1,alpha2,alpha3
 double precision  :: corr,ve
 integer           :: ierror,n0,i1,i2,i3,i4,i5,i6,itertype,last_i,robustfit
 double precision  :: stab_best,acc_rms,conf_int
 integer           :: maxivar,minivar,tivar
 double precision  :: last_stadev,factor_back,factor_out
 double precision  :: rms,s,da1,da2,da,ZPE,stadev_best,TEMPE,gamma,vup,emax
 integer           :: last_ivar,ndigits,npoints_delta,rand_itmin,alloc,alloc_p
 integer           :: max_pnts,ipar,i0,l0,parmax,fivar_t,nmax,taumax
 logical           :: yes,still_run,ifopen,find_best_params
 integer           :: flag_borh, flag_hartre,ivartmp,it0,imin,imax,ilarge_sterr,ilow_sterr,Nout,parmax_t,i_t
 double precision  :: conv_length,conv_energy,paramtmp,relat_sterr,relat_sterr_min,fparam_t,beta312,alpha34,V0
 character(len=68)  :: fmt1,fmt2
 character(len=2)   :: fmt0
 character(len=130)  ::  char_tmp 
 real(ark)           :: w_t
 integer           ::  iparam,NDEG,NBONDS,NANGS,nparam_eq
! MADE NPARAM_EQ A 
 DOUBLE PRECISION, ALLOCATABLE :: IT(:), LOCAL(:),local_p(:)
 ! NDEG ADDED BY BPM TO MAKE CODE MORE GENERAL, NOW READ IN NUMBER OF DEGREES
 ! AND ALLLOCATE MATRICES
 ! ALSO ADDED NBONDS
 !
 ! some matrices, to be allocated  
 !
 double precision, allocatable :: rjacob(:,:),last_param(:),param(:)
 double precision, allocatable :: param_tmp(:)
 double precision, allocatable :: crr(:),coord(:,:)
 double precision, allocatable :: energy(:)
 double precision, allocatable :: eps(:),wt(:),parold(:),sigma(:),wt_tmp(:),parinitial(:)
 double precision, allocatable :: al(:,:),bl(:),dx(:)
 double precision, allocatable :: ai(:,:),sterr(:),dV(:,:),f(:)
 double precision, allocatable :: sens(:)
 character(len=10), allocatable:: parnam(:),parnam_tmp(:)
 integer, allocatable :: ivar(:),icorr(:),ivar0(:),ivar_tmp(:),ind_tmp(:),ipower_tmp(:,:),ipower(:,:),icol(:),IPOWER_(:)
 !
 double precision, allocatable :: Tsing(:,:) 
 integer                       :: rank0,info
 integer                       :: lwork,kcol
 double precision,allocatable  :: wspace(:)
 double precision              :: alpha13,alpha12,alpha14,alpha23,alpha24,cosbeta,aacos,beta412,cosa34
 character(len=40)        ::  txt
 double precision t12,t13,t23,t45,t46,t56,t24,t15,t16,t34,t35,t26,t14,tbar,f_t,t43,t25,t36
 double precision A1,A2,A3,A4,redundancy



 !
 ! Common block is to reduce the number of assignments within the potential routine
 !
 ! Here we go!

 !
 !  initial constants 
 !
 last_stadev = 1000000.
 maxivar = 8
 maxdeviat = 1e8
 !
 pi = 4.0d0 * datan2(1.0d0,1.0d0)
 !
 ! Input data from file: 
 !
 call skiplines(f_inp,1)
 !
 !  input the job title 
 !
 do i=1,4
   read  (f_inp,"(a80)") title(i)
 enddo
 !
 !  output the job title 
 !
 write (f_out,"(3(132('*')/' '),16('*'),100x,16('*'))") 
 write (f_out,"(4(' ',16('*'),10x,a80,10x,16('*')/))") (title(i), i=1,4)
 write (f_out,"(' ',16('*'),100x,16('*')/' ',3(132('*')/' '))") 
 !
 call skiplines(f_inp,4)
 !
 READ (f_inp,*) NDEG
 READ (f_inp,*) NBONDS
 READ (f_inp,*) NANGS
! ADDED IN BY BPM TO MAKE MORE GENERAL
 ALLOCATE(IPOWER_(NDEG)) 
 ALLOCATE(IT(NDEG))
 ALLOCATE(LOCAL(NDEG))
 !
 read(f_inp,*) gamma,emax,vup
!
 read (f_inp,*) itmax
 read (f_inp,*) stab_best,stadev_best
 read (f_inp,*) conv_length
 read (f_inp,*) conv_energy
 read (f_inp,*) factor_out
 read (f_inp,*) factor_back
 read (f_inp,*) robustfit
 !
 WRITE(f_out,"(i10,' <= number of degrees of freedom')")         NDEG
 WRITE(f_out,"(i10,' <= number of bonds')")         NBONDS
 WRITE(f_out,"(i10,' <= number of bond angles')")         NANGS
 WRITE(f_out,"(f10.6,' <= Schwenkes gamma,emax,Vup')")    gamma,emax,vup
 write (f_out,"(i10,' <= number of iterations in the fit')")         itmax
 write (f_out,"(f10.6,' <= The fit is over if the stability is reached ')") stab_best
 write (f_out,"(f10.6,' <= The fit is over if the st.dev. best is reached ')") stadev_best
 write (f_out,"(f10.6,' <= sterr parameter out ')")                  factor_out
 write (f_out,"(f10.6,' <= stdev parameter out ')")                  factor_back
 write (f_out,"(i10,  ' <= watson robust on (1) or of (0) ')")       robustfit
 !
 if (robustfit<0) then 
   find_best_params = .false.
 else
   find_best_params = .true.
 endif 
 !
 ! Converiosn factors:
 ! bond length 
 !
 !if (flag_borh ==1) then 
 !  conv_lenght = bohr 
 !  write (f_out,"(' The bond lengths are in bohr')") 
 !else 
 !  conv_lenght = 1.0
 !  write (f_out,"(' The bond lengths are in Angstrom')") 
 !endif   
 !
 ! energies 
 !if (flag_hartre==1) then 
 !  conv_energy = hartre
 !  write (f_out,"(' The ab initio energies are in hartre')") 
 !else 
 !  conv_energy = 1.0 
 !  write (f_out,"(' The ab initio energies are in 1/cm')") 
 !endif   


 !
 ! Array's allocation:
 !

 parmax_t = nparam
 lwork =  parmax_t*8
 !
 allocate (param_tmp(parmax_t),parnam_tmp(parmax_t),ivar_tmp(parmax_t),ipower_tmp(1:NDEG,parmax_t),stat=alloc)
 if (alloc/=0) then 
  write(6,"('param_tmp - out of memory')")  
    stop 'param_tmp - out of memory'
 endif 

 ipower_tmp = 0


 !
 ! input  parameters of structur
 !
 parnam_tmp = 'xxxxxx'
 param_tmp = 0 
 ivar_tmp = 0 
 !
 call skiplines(f_inp,3)
 READ (f_inp,*) NPARAM_EQ
 WRITE(f_out,"(i10,' <= number of non-linear parameters')")  NPARAM_Eq
 do i = 1,nparam_eq
   read (f_inp,"(a80)") longlabel
   read (longlabel,*) label,ivartmp,paramtmp
     parnam_tmp(i) = label 
     ivar_tmp(i)   = ivartmp
     param_tmp(i)  = paramtmp
 enddo 

 if (verbose>=5) then 
   write (f_out,"(' parameters of structure ')")
   do i = 1,nparam_eq
      write (f_out,"(a10,'=',f16.3)") parnam_tmp(i),param_tmp(i)
   enddo 
 endif
 !
 Yes = .True. 
 i = nparam_eq
 nmax = 0
 taumax=0
 read_loop: do while ( Yes )
    read (f_inp,"(a130)") char_tmp
    if ( char_tmp(1:3)/='end' ) then
       !
       i = i + 1
       !
       read (char_tmp(1:130),*) parnam_tmp(i),ipower_(1:NDEG),fivar_t,fparam_t
!      THIS PART COULD DO WITH A FORMAT STATEMENT BPM OR JUST DON'T PRINT
       if (verbose>=3) then 
          write (f_out,"(a10,19i3,1x,g16.9)") trim(parnam_tmp(i)),ipower_(1:NDEG),fivar_t,fparam_t
       endif
       !
       param_tmp(i) = fparam_t
       ivar_tmp(i) = fivar_t
       ipower_tmp(1:NDEG,i) = ipower_(1:NDEG)
       !
    else
       Yes = .False.
    endif
 enddo read_loop
 !
 parmax = i 
 !
 allocate (param(parmax),parnam(parmax),ivar(parmax),icorr(parmax),ivar0(parmax),ipower(1:NDEG,parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('param - out of memory')")  
    stop 'param - out of memory'
 endif 
 !
 param(1:parmax) = param_tmp(1:parmax)
 ivar(1:parmax) = ivar_tmp(1:parmax)
 ipower(1:NDEG,1:parmax) = ipower_tmp(1:NDEG,1:parmax)
 parnam(1:parmax) = parnam_tmp(1:parmax)
 !
 deallocate (param_tmp,parnam_tmp,ipower_tmp,ivar_tmp)
 !
 allocate (parold(parmax),parinitial(parmax),last_param(parmax),sterr(parmax),&
           Tsing(parmax,parmax),Wspace(lwork),crr(parmax),al(parmax,parmax),bl(parmax),&
           dx(parmax),ai(parmax,parmax),coord(NDEG,enermax),energy(enermax),&
           wt_tmp(enermax),sens(parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('sens - out of memory')")  
  stop
 endif 

 parinitial = param 

! BPM MAKE FIRST LINEAR PARAMATER (ALL 0 COEFFCIENTS ) V0 FOR SUBTRACTING FROM
! ENERGY LATER (CURRENTLY NOT USING)

 V0 = PARAM(NPARAM_EQ+1) 

!---input energy values------!
 npts = 0
 yes = .true. 
 Nout = 1
 do while( yes )
    npts=npts+1
    if (npts>enermax) then
      write(f_out,"('Too many ab initio points, increase enermax:',I6,' vs ',I6)") npts,enermax
      endif
    read (f_inp,*) coord(1:NDEG,npts),energy(npts),wt_tmp(npts)
    if ( wt_tmp(npts).lt.0.0 ) wt_tmp(npts) = abs(wt_tmp(npts))*factor_back
    if ( coord(1,npts).lt.0.0 ) yes = .false.
 enddo 
 !
 ! Number of ab initio points  
 npts=npts-1
 !
 allocate (rjacob(npts,parmax),eps(npts),wt(npts),sigma(npts),dV(npts,parmax),stat=alloc)
 if (alloc/=0) then 
  write(6,"('rjacob - out of memory')")  
  stop
 endif 
 !
 wt(1:npts) = wt_tmp(1:npts)
 !
 !
 ! BPM modify coordinates from MOLRPO TO COORDINATES WANTED FOR FIT HERE
 !
 do i = 1, npts
    call transform_coordiates_tau_11(ndeg,coord(:,i))
 end do
 !
 !
 allocate(icol(parmax))
 ncol=0
 do  i=1,parmax
     if (ivar(i) .ne. 0) then
       ncol = ncol + 1
       icol(ncol) = i
     endif
 enddo
 !
 numpar  = 0
 do i=1,parmax
   if (ivar(i) .gt. 0) numpar=numpar+1
 enddo
 !
 !
 ! We introduce a "zero point energy" shift to avoid loosing accuracy 
 ! because of too too large numbers
 ZPE = 0.0 !energy(1)
 ! Convert the energies to the internal dimension, e.g. 1/cm 
 energy(:) = ( energy(:)-ZPE )*conv_energy
 !
 ! For the Robust fit we need to define accuracy of the abinitio data "sigma"
 !
 !
 !      ADD IN WEIGHTING FUNCTION HERE. NEED TO CHANGE MAX ENERGIES ETC FOR
 !      DIFFERENT WEIGHTING. COULD ADD AS INPUT PARAMETER
 !
 do i = 1, npts
   tempe = energy(i) - param(nparam_eq+1)
   w_t = (tanh(-real(gamma,ark)*( real(tempe,ark) - real(Vup,ark))) + 1.000000200000020000002_ark)/2.000000200000020000002_ark
   wt(i) = w_t
   !
   if (tempe>emax) wt(i) = 0 
   !wt(i) = (tanh(-gamma*( max(tempe,emax) - Vup)) +1.000000200000020000002)/2.000002000002000002
 !  write(6,*) energy(i), tempe, wt(i)
 end do
 ! 
 ! "normalising" the weight factors
 nused=0
 wtsum=0.0d0
 do i=1,npts
   if (wt(i) > 0.0) nused=nused+1
   wtsum=wtsum+wt(i)
 enddo 
 wtsum=wtsum/nused
 wt(:)=wt(:)/wtsum
 !
 do nrow=1,npts
   if (wt(nrow)>0.0) then 
     sigma(nrow) = 1.0d0/wt(nrow)
   else 
     sigma(nrow) = 1000.0 
   endif 
 enddo
 !
 ! a_wats = 1.0d0
 !
 ! Sometimes pot. parameters are dependent. They are nor to be in the fit 
 ! If we meet a dependent parameter, we get rid of it and start the fit again 
 ! This is the outer loop for
 !
 !
 still_run = .true.
 outer_loop: do while (still_run)  
   !
   !  make a table header 
   !
   if (itmax.ne.0) then 
     write(f_out,"(/3x,66('-'))") 
     write(f_out,"('   |   iter  |  points |    deviat     |     rms       |  stability  |')")
     write(f_out,"(3x,66('-'))") 
   endif 
   !
   ! Parameters to control current and previous st.deviations and stability
   !
   stadev_old = 1.e10
   stability =  1.e10
   stadev    =  1.e10
   ! 
   ! Parameters to control excluding dependent parameters   
   ! We need to know the maximum and minimum values of the pot. parameters weights "ivar" 
   !
   ! numpar is to count the number fitted varying pot. parameters  
   !
   numpar  = 0
   maxivar = 0
   minivar = 1000000
   do i=1,parmax
     if (ivar(i) .gt. 0) numpar=numpar+1
     if (ivar(i) .gt. maxivar) maxivar = ivar(i)
     if (ivar(i) .lt. minivar .and. ivar(i).ne.0) minivar = ivar(i)
   enddo 
   !
   if (nused<numpar) then 
       write(f_out,"('worning! the number of parameter <  the number of data points',2I6)") nused,numpar
       !stop 'worning! the number of parameter <  the number of data points'
    endif
   if (nused==numpar) then 
      write(f_out,"('worning! the same number of parameters and data points')")
   endif 
   ncol=0
   do  i=1,parmax
       if (ivar(i) .ne. 0) then
         ncol = ncol + 1
         icol(ncol) = i
       endif
   enddo
   !
   !##################################################c
   !############ start of the fitting ################c
   !##################################################c
   rjacob = 0 
   iter = 0
   sens = 0.d0 
   do while( iter<=itmax .and. stadev>stadev_best .and. (  stability.ge.stab_best) )   
     iter = iter + 1
     ssq=0.0d+00
     rms=0.0d+00
     !
     parold = param
     !
     dV = 0 
     !
     if (verbose>=6) print*,"Generate PES derivatives"
     !
     !$omp parallel private(f,local_p) shared(dV)
     allocate(f(parmax),local_p(ndeg),stat=alloc_p)
     if (alloc_p/=0) then 
       write (f_out,"('f - out of memory')")
       stop 'f - out of memory'
     endif
     !
     !$omp do private(nrow) schedule(dynamic)
     do nrow=1,npts
       !
       !$omp critical
       if (verbose>=6) write(6,"('n = ',i9)") nrow
       !$omp end critical
       !
       ! MODIFIED BY BPM TO MAKE MORE GENERAL WITH NBONDS/NANGS. NEED TO
       ! TEST/CHECK
       !
       ! Running bond lengths
       !
       local_p(1:NBONDS)=coord(1:NBONDS,nrow)*conv_length
       !
       ! Running interbond angles
       !
       local_p(NBONDS+1:NDEG)=coord(NBONDS+1:NBONDS+NANGS,nrow)*pi/180.0d0
       !
       !do iparam = nparam_eq+1,parmax
         !
         call diff_V(ipower,local_p,param(1:nparam_eq),NDEG,NBONDS,NANGS,nparam_eq,parmax,f)
         dV(nrow,1:parmax) = f(1:parmax)
         !
       !enddo
       !
     enddo
     !$omp end do
     !
     !deallocate (f)
     deallocate(f,local_p,stat=alloc_p)
     if (alloc_p/=0) then 
       write (f_out,"('f - out of memory')")
       stop 'f - out of memory'
     endif
     !$omp end parallel 
     !
     if (verbose>=6) print*,"...done!"
     !
     if (verbose>=6) print*,"Check redundancy ..."
     do ipar = nparam_eq+1,parmax
       !
       redundancy = sum(dv(:,ipar)**2)
       !
       ! exclude from the fit if derivative at all geometetries is zero
       !
       if (redundancy<1e-8) then 
         !
         !
         if (verbose>=6) then
           !
           !$omp critical
           if (ivar(ipar)/=0.and.abs(param(ipar))>1e-15) write(f_out,"(i8,'-th is out - ',a10)") ipar,parnam(ipar)
           !$omp end critical
           !
         endif
         !
         !param(ipar) = 0 
         ivar(ipar) = 0 
         !
       endif
       !
     enddo
     !
     if (verbose>=6) print*,"...done!"
     !
     if (verbose>=6) print*,"...Generate the fitting Jacobi matrix"
     !
     do nrow=1,npts
       !
       !print*,"nrow = ",nrow
       !
       ! Running bond lengths
       !
       local(1:NBONDS)=coord(1:NBONDS,nrow)*conv_length
       !
       ! Running interbond angles
       !
       local(NBONDS+1:NBONDS+NANGS)=coord(NBONDS+1:NANGS+NBONDS,nrow)*pi/180.0d0
       ! 
       ! Value of the potential energy function at the current geometry 
       !
       v = potential(parmax,dV(nrow,1:parmax),local,param,ipower,0,NDEG,NBONDS,NANGS,NPARAM_EQ)
       !
       ! calculate derivatives with respect to parameters (RJACOB)
       !
       if (itmax/=0 .and. wt(nrow)>0.0) then
         !$omp parallel do private(kcol,i) shared(rjacob) schedule(dynamic)
         do kcol=1,ncol                
         !
         i = icol(kcol)
         !

      !   ncol=0
      !   do  i=1,parmax
      !     if (ivar(i) .ne. 0) then
      !       !
      !       ncol=ncol+1
             rjacob(nrow,kcol)=potential(parmax,dV(nrow,1:parmax),local,param,ipower,i,NDEG,NBONDS,NANGS,nparam_eq)
             !
       !    endif
         enddo ! --- ncol
        !$omp end parallel do
       endif
       !
       eps(nrow) = energy(nrow)-v
       !
       ! Weighted st. square deviation 
       !
       ssq=ssq+eps(nrow)*eps(nrow)*wt(nrow)
       !
       ! mean square 
       !
       rms=rms+eps(nrow)*eps(nrow)
     enddo  ! ---  nrow 
     !
     if (verbose>=6) print*,"Construct and solve a linear system..."
     !
     ! We constract a set of linear equations A x = B
     !
     if (itmax.ne.0) then
        ! form A matrix 
        do irow=1,numpar       !==== row-...... ====!
          do icolumn=1,irow    !==== column-....====!
            al(irow,icolumn)=sum(rjacob(:,icolumn)*rjacob(:,irow)*wt(:))
            al(icolumn,irow)=al(irow,icolumn)
          enddo
        enddo
        ! form B matrix 
        do irow=1,numpar       !==== row-...... ====!
          bl(irow)=sum(eps(:)*rjacob(:,irow)*wt(:))
        enddo   
        !
        ! Solve the set of linear equations 
        !
       select case (trim(fit_type))
       !
       case default
         write (6,"('fit_type ',a,' unknown')") trim(fit_type)
         stop 'fit_type unknown'
       case('linur')
         !
         call linur(numpar,parmax,al,bl,dx,ierror)
         !
          if (ierror.ne.0) then
            ncol=0
            do i=1,parmax
              if (ivar(i) .ne. 0) then
                ncol=ncol+1
                if ( ncol.eq.ierror ) then
                    ivar(i) = 0
                    param(i) = 0 
                    if (ireplace_by/=0) param(i) = parinitial(i)
                    write(f_out,"(i8,'-th is out - ',a7)") i,parnam(i)
                endif
              endif
            enddo
            cycle outer_loop
          endif
       case ('dgelss')
         !
         write(6,"('DGELSS: lapack was not presented ')")
         stop 'dgelss not compiled'
         !
         ai = al
         !
         !
         !call DGELSS(numpar,numpar,1,AL,numpar,BL,numpar,Tsing,1.D-12,RANK0,wspace,lwork,ierror)
         !
         if (info/=0) then
           write(6,"('DGELSS:error',i8)") ierror
           stop 'dgelss'
         endif
         !
         dx = bl
         !
       end select
       !
       !
       if ( nused.ne.numpar ) then 
          stadev=dsqrt(ssq/float(nused-numpar))
       else 
          stadev=dsqrt(ssq/nused)
       endif
       !
        ncol=0
        do i=1,parmax
          if (ivar(i) /= 0) then
             ncol=ncol+1
             param(i)=param(i)+dx(ncol)
          endif
        enddo
        !
        ! Calcualte the inverse to A matrix, needed for the standard error of the fitted pot. parameters
        !
        call invmat(al,ai,numpar,parmax)
        !
        !
        ncol = 0 
        do i=1,parmax
          if (ivar(i) /= 0) then
             ncol=ncol+1
            !  
            ! If nused = numpar, the problem is degenerated, Ax = B is a system of ordinary equations
            ! st. error doesn't make sense  
            !
            if (nused==numpar) then  
              sterr(ncol)=0
            else
              !
              ! Standard definition of the st. error 
              !
              sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
              !
            endif
          endif
        enddo    
        !
        ! Here we define stability to see if the fit is converged 
        !      
        stability=abs( (stadev-stadev_old)/stadev )
        stadev_old=stadev
        !
      if (RobustFit>0) then 
        !
        !Watson alpha-parameter
        ! 
        !do it0 = 1,20 
        ! da1 = 0
        ! da2 = 0
        ! do nrow=1,npts
        !   if (wt(nrow) .gt. 0.0d0) then 
        !     da1 = da1+eps(nrow)**2/( sigma(nrow)**2+a_wats*eps(nrow)**2 )
        !     da2 = da2+eps(nrow)**4/( sigma(nrow)**2+a_wats*eps(nrow)**2 )**2
        !   endif 
        ! enddo 
        ! !
        ! da =( da1 -  float(nused-numpar) )/da2
        ! a_wats = a_wats + da
        ! !
        ! if (a_wats<0.0000001) a_wats = 0.001+it0*0.01
        !enddo
        !
        ! 
        !  adjusting the weights  
        ! 
        ! We try a_wats = 0.25        
        ! a_wats = 0.25        
        do nrow=1,npts
           if (wt(nrow) .gt. 0.0d0) then 
             wt(nrow) = 1.d0/( sigma(nrow)**2 + a_wats*eps(nrow)**2 )
           endif 
        enddo 
        ! 
        ! "re-normalising" the weight factors
        !
        nused=0
        wtsum=0.0d0
        do i=1,npts
          if (wt(i) > 0.0) nused=nused+1
          wtsum=wtsum+wt(i)
        enddo 
        wtsum=wtsum/nused
        wt(1:npts)=wt(1:npts)/wtsum

        ssq = sum(eps(1:npts)**2*wt(1:npts))
        !
        ! We can recalculated the stand. deviation with new weigh factors
        ! otherwise stadev is defined as a weighted RMS 
        !
        if (nused.eq.numpar) then 
          stadev=dsqrt(ssq/float(nused))
        else 
          stadev=dsqrt(ssq/float(nused-numpar))
        endif 
      endif ! --- robust fit
            
     else   ! Itermax = 0, so there is no fit
            ! only straightforward calculations 
        !
        stadev_old=stadev
        stadev=dsqrt(ssq/float(nused))
     endif 
     !
     if (verbose>=6) print*,"...done!"
     !
     ! Here we define the root mean square 
     !
     rms = dsqrt(rms/npts)
     !
     ! Do some output 
     !
     if (itmax/=0 .and. stadev < maxdeviat ) then 
         write (f_out,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |')") &
                iter,nused,stadev,rms,stability
     endif 
     !
   enddo  ! --- iter
   !################ end of iterations  #######################
   !
   ! update the energies and residuals based on the last change 
   !
   if (verbose>=6) print*,"Update rms..."
   !
   do nrow=1,npts
     !
     ! Running bond lengths
     !
     local(1:NBONDS)=coord(1:NBONDS,nrow)*conv_length
     !
     ! Running interbond angles
     !
     local(NBONDS+1:NBONDS+NANGS)=coord(NBONDS+1:NANGS+NBONDS,nrow)*pi/180.0d0
     ! 
     ! Value of the potential energy function at the current geometry 
     !
     v = potential(parmax,dV(nrow,1:parmax),local,param,ipower,0,NDEG,NBONDS,NANGS,NPARAM_EQ)
     !
     eps(nrow) = energy(nrow)-v
     !
     ! Weighted st. square deviation 
     !
     ssq=ssq+eps(nrow)*eps(nrow)*wt(nrow)
     !
     ! mean square 
     !
     rms=rms+eps(nrow)*eps(nrow)
     !
   enddo  ! ---  nrow 
   !
   rms = dsqrt(rms/npts)
   !
   if ( nused.ne.numpar ) then 
      stadev=dsqrt(ssq/float(nused-numpar))
   else 
      stadev=dsqrt(ssq/nused)
   endif
   !
   if (RobustFit>0) write(f_out,"(/'Watson Alpha parameter',f8.4,' fitted with rms2 = ',f8.4)") &
                                     a_wats !,da1/dfloat((nused-numpar))
   !
   ! if RobustFit==2,3 after first iteration it converts to the standard fit RobustFit=0, but with robust weights
   !
   if ( RobustFit==2 ) then
     RobustFit = 0 
   else if ( RobustFit==3 ) then ! Averaging after the first iteration 
     do i=1,npts
       imin = max(   1,i-npts/ishift) 
       imax = min(npts,i+npts/ishift) 
       ssq = sum(wt(imin:imax))/(imax-imin)
       wt_tmp(i) = ssq
     end do
     wtsum=sum(wt_tmp(1:npts))
     wtsum=wtsum/nused
     wt(1:npts)=wt_tmp(1:npts)/wtsum
     RobustFit = 0 
     cycle outer_loop    
   end if 
   !
   ! We need to know the equilibrium value of the pot. function to substract it from the 
   ! output energies 
   !
   ve = param(nparam_eq+1) ! -ZPE*conv_energy

   ! Output some staqtistics and results 
   !
   !  only if we are fitting:  
   !
   if (itmax.ne.0.and.find_best_params) then
      !
      ! 
      ! It's time to do the output of pot. parameters 
      ! We print them out rounded with their st.errors
      ! 
      write (f_out,"(//'Potential parameters rounded in accord. with their standart errors'/)")
      l = 0 
      do i=1,parmax
        if (ivar(i) .ne. 0) then
           l=l+1
           ndigits = 0
           conf_int = sterr(l)
           do while (conf_int.le.10.0)
             ndigits = ndigits +1 
             conf_int = conf_int*10
           enddo
           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           fmt1 = '(a2,I4,2x,f14.'//fmt0//'    ''    ('',i8,'')'')'
           write (f_out,FMT=fmt1) trim(parnam(i)),ivar(i),param(i),nint(conf_int)
           !
        else 
           ndigits =2
           if (param(i).ne.0.0) ndigits = 8

           write(fmt0,"(i2)") ndigits
           fmt0 = adjustl(fmt0)
           !
           fmt1 = '(a2,18i3,I4,2x,f14.'//fmt0//')'
           if (i<=nparam_eq) then 
             fmt1 = '(a10,21x,I4,2x,f14.'//fmt0//')'
           endif
           !
           fmt1 = '(a2,18i3,I4,2x,f14.'//fmt0//')'
           write (f_out,FMT=fmt1) trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)


           !do i=1,nparam_eq
           !   write (f_out,"(a10,27x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
           !enddo 
           !do i=nparam_eq+1,parmax
           !   write (f_out,"(a10,9i3,i4,2x,f18.8)") parnam(i),ipower(1:9,i),ivar(i),param(i)
           !enddo 

           !
        endif
      enddo  ! --- i


     write (f_out,"('------------------------------')")
    
     write (f_out,"(a80)") (title(i), i=1,4)

     do i=1,nparam_eq
        write (f_out,"(a10,21x,i4,2x,f18.8)") trim(parnam(i)),ivar(i),param(i)
     enddo 
     do i=nparam_eq+1,parmax
        write (f_out,"(a2,18i3,i4,2x,f21.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
     enddo 

     !
      7007  format(a10,I4,2x,f14.4,'    (',i8,')')
      7008  format(a10,I4,2x,f14.4)
      !
   endif
   !
   if (itmax.ne.0.and.find_best_params) then
     !
     if (verbose>=6) print*,"Optimize the fitting set by exluding correlated parameters ..."
     !
     ! If any of the pot. parameters (last_i) has been removed from the fit at the previous iteration step
     ! at this moment we check whether the st. deviation is worse than it was before we got rid of the pot. parameter 
     ! in this case we take the last_i-pot. parameter back (for one parameter)
     ! this very tricky process is controled by factor_back 
     ! 
     if ( (stadev-last_stadev)/stadev > Nout*factor_back .and.last_i/=0) then 
       ivar = ivar0     
       ivar(last_i)  =  min(ivar(last_i)+1,9)
       param = last_param
       last_stadev = 10000000.
       write(f_out,"(I6,'-th is back - ',A10)") last_i,parnam(last_i)
       cycle outer_loop    
     endif 

      !
      ! Another statistics of the fit: 
      ! correlations between pot.parameters 
      !
      write (f_out,"(/'Correlation coefficients larger than 0.9: ')")
      do i1=1,numpar
        do i2=i1+1,numpar
          corr=ai(i1,i2)/sqrt(abs(ai(i1,i1)*ai(i2,i2)))
           if (abs(corr) .gt. 0.9) write (f_out,"(10x,'corr(',i2,',',i2,') = ',f10.7)") i1,i2,abs(corr)
        enddo
      enddo
      ! 
      ! sensitivity
      !
      ncol = 0 
      do i=1,parmax
        if (ivar(i) /=  0) then
          ncol=ncol+1
          if (nused /= numpar) then  
            ssq = sum(rjacob(:,ncol)**2*wt(:))
            sens(ncol)=stadev*0.1d0/dfloat(numpar)/dsqrt( ssq/dfloat(nused) )
          endif
        endif
      enddo    
      ! 
      ! Robust fit secion 
      !
      !
      ! Pot. parameters output with their st. errors 
      !
      write (f_out,"(/15x,80('-')/15x,':     old parm   :     new parm   :  delta parm  :   std. err.  : sensitility'/15x,80('-'))")
      l = 0
      do i=1,parmax
        if (ivar(i) /= 0) then
          l=l+1
          !
          epsil=param(i)-parold(i)
          crr(l)=(sterr(l)*al(l,l)/stadev)**2
          !
          mark  = ' '
          if (abs(sens(l)) >= abs(epsil) ) mark = '!'
          !
          if (abs(param(i)) > abs(sterr(l))) then
             write (f_out,"(i3,2x,a10,1x,2(':',f14.2,2x),':',3(f14.3,':'),a1)") &
                    l,parnam(i),parold(i),param(i),epsil,sterr(l),sens(l),mark
          else
              write (f_out,"(i3,2('#'),a10,1x,2(':',f14.2,2x),':',3(f14.3,':'),a1,'#')") & 
                      l,parnam(i),parold(i),param(i),epsil,sterr(l),sens(l),mark
          endif
        endif
      enddo



    inquire(f_res,opened=ifopen)
    if ( ifopen ) then
      rewind(f_res)
    else
      open  (f_res,file='00.res',status='replace' )
    endif
    !
    if (itmax.ne.0) then 
     write (f_res,"(a80)") (title(i), i=1,4)
      do i=1,nparam_eq
        write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
      enddo 
      do i=nparam_eq+1,parmax
       if (abs(param(i))>1e8) then
         write (f_res,"(a2,18i3,i4,2x,g19.9)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       else
         write (f_res,"(a2,18i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       endif
      enddo 
    endif  ! ---- itmax/=0
    !
    write (f_res,"('------------------------------')")
    !
    write (f_res,"(a80)") (title(i), i=1,4)
    do i=1,nparam_eq
       write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
    enddo 
    do i=nparam_eq+1,parmax
      if (abs(param(i))>1e-12) then 
       if (abs(param(i))>1e8) then
         write (f_res,"(a2,18i3,i4,2x,g19.9)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       else
         write (f_res,"(a2,18i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       endif
      endif
    enddo 
    !
    do nrow=1,npts
      !
       v = energy(nrow)-eps(nrow)
       write (f_res,"(7f15.6,11f11.4,3f18.7,2x,e9.2)")            & 
             coord(1:NDEG,nrow),                                  &
             energy(nrow)-ve,v-ve,eps(nrow),                   &
             wt(nrow)

    enddo

      !
      !  Here we check if there is a pot. parameter with too big st. error  
      !  We consider a parameter to be unneccesary and remove from the fit 
      !  keeping in mind that we might have to take it back in case the st.deviation will get worse
      !
      !
      ! Search for the largest relative sterr
      !
      tivar = minivar
      ivar0 = ivar
      icorr = 0
      last_param = param
      Nout = 0
      do while ( tivar.lt.maxivar .and.factor_out>0)
        if (verbose>=5) write (f_out,"(' tivar, factor_out =',i5,f18.8)") tivar,factor_out
        relat_sterr_min = 1.0e20
        ilow_sterr = -1

        relat_sterr = 0.0
        ilarge_sterr = -1

        l = numpar+1
        do i=parmax,1,-1
          if (ivar(i) /= 0) then
            l=l-1
            !
            ! Get read of parameters with large rel.st.error and their dependencies 
            !
            if (verbose>=5) write (f_out,"(' param,number, ivar, st.error ',f18.7,i5,i5,f18.7)") &
                 param(i),i,ivar(i),abs(sterr(l)/param(i))
            !
            if ( i>nparam_eq+1.and.abs(param(i))>too_big.and.ivar(i)/=maxivar ) then
              last_param = param
              last_stadev = stadev
              last_i = i
              ivar(i) = 0
              !
              ! The parameter gets zero value  
              !
              param(i) = 0.d0
              write(f_out,"(' parameter ',i6,' is  out - ',A10)") i,parnam(i) 
              cycle outer_loop    
            endif 
            !
            if (param(i)/=0.0d0.and.ivar(i) == tivar) then
              if (icorr(i)==0.and.factor_out<abs(sterr(l)/param(i))) then
                if (verbose>=5) write (f_out,"(' param, number, st.err., ivar',f18.7,i5,f18.7,i5)") &
                 param(i),i,abs(sterr(l)/param(i)),ivar(i)
                ! relat_sterr=abs(sterr(l)/param(i))
                ! ilarge_sterr = i
                icorr(i) = 1
                l0=0
                relat_sterr = abs(sterr(l)/param(i))
                ilarge_sterr = i
                do i0=1,parmax,1
                  if (ivar(i0) /= 0.and.ivar(i0) /= maxivar) then
                    l0=l0+1
                    corr=abs(ai(l,l0))/sqrt(abs(ai(l,l)*ai(l0,l0)))
                    !   if (verbose>=4) then 
                    !     write (f_out,"('corr(',i,',',i,') = ',f14.7)") l,l0,corr
                    !   endif 
                    if (i0/=i.and.param(i0)/=0.0d0.and.ivar(i0) == tivar) then
                      if (icorr(i0)==0.and.corr>fact_corr ) then
                        if (verbose>=4) write (f_out,"('corr, param, (st. err.)  number',2f18.7,f18.7,i4)") & 
                                      corr,param(i0),sterr(l0),i0
                        icorr(i0) = 1
                        !
                        ! Here one finds the (lowest) largest st. error and which par. it belongs to
                        if ( relat_sterr<abs(sterr(l0)/param(i0))) then
                            relat_sterr=abs(sterr(l0)/param(i0))
                            ilarge_sterr = i0
                        endif
                      endif
                    endif
                  endif
                enddo
                ivar(ilarge_sterr)  = 0
                param(ilarge_sterr) = 0.0
                if (ireplace_by/=0) param(ilarge_sterr) = parinitial(ilarge_sterr)
                if (verbose>2) then 
                  write(6,"('relat_sterr,ilarge_sterr, tivar:',f18.8,2i5)") relat_sterr,ilarge_sterr, tivar
                endif 
                Nout = Nout+1
                !
                ! Here one finds the (lowest) largest st. error and which par. it belongs to
                if ( relat_sterr_min>relat_sterr) then
                    relat_sterr_min=relat_sterr
                   ilow_sterr = ilarge_sterr
                endif
              endif
            endif
          endif
        enddo
         if (verbose>2) then 
           write(f_out,"('relat_sterr_min,ilow_sterr, tivar,rms:',d18.8,2i5,e11.3)") relat_sterr_min,ilow_sterr, tivar,rms
        endif 
        if (verbose>=4) then 
          do i=1,nparam_eq
             write (f_out,"(a10,22x,i4,2x,f19.8)") parnam(i),ivar(i),param(i)
          enddo 
          do i=nparam_eq+1,parmax
             write (f_out,"(a2,18i3,i4,2x,f19.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
          enddo 
        endif 
        !
        if (relat_sterr_min > factor_out .and. ilow_sterr/=-1) then
          last_ivar = ilow_sterr
          ! last_param = param
          last_stadev = stadev
          last_i = ilow_sterr
          ! ivar(ilarge_sterr) = 0
          !
          ! The parameter gets zero value  
          !
          ! param(ilarge_sterr) = 0.d0
          write(f_out,"(I6,' parameters are out - ',A10)") Nout
          cycle outer_loop    
        endif 
        tivar = tivar + 1
      enddo

      still_run = .false.


!     tivar = minivar
!     !
!     do while ( tivar.lt.maxivar )
!       l = numpar+1
!       !
!       do i=parmax,1,-1
!         if (ivar(i) .ne. 0) then
!           l=l-1
!           epsil=param(i)-parold(i)
!           crr(l)=(sterr(l)*al(l,l)/stadev)**2
!           if (factor_out*abs(param(i)) < abs(sterr(l)) .and. ivar(i) == tivar) then
!              last_ivar = ivar(i)
!              last_param = param
!              last_stadev = stadev
!              last_i = i
!              ivar(i) = 0
!              !
!              ! The parameter gets zero value  
!              !
!              param(i) = 0.d0
!              write(f_out,"(I6,'-th is out - ',A10)") i,parnam(i)
!              cycle outer_loop    
!           endif
!         endif
!       enddo
!       tivar = tivar + 1
!     enddo

    endif 
    still_run = .false.

!--  printing  out the resulted information 




 ! Here we can stop the fit 
 enddo outer_loop   


!--  printing  out the resulted information 

   inquire(f_res,opened=ifopen)
   if ( ifopen ) then
     rewind(f_res)
   else
     open  (f_res,file='00.res',status='replace' )
   endif

  if (itmax.ne.0) then 
    write (f_res,"(a80)") (title(i), i=1,4)
    do i=1,nparam_eq
       write (f_res,"(a10,21x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
    enddo 
    do i=nparam_eq+1,parmax
      if (abs(param(i))>1e-12) then 
       if (abs(param(i))>1e8) then
         write (f_res,"(a2,18i3,i4,2x,g19.9)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       else
         write (f_res,"(a2,18i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       endif
      endif
    enddo 
  endif  ! ---- itmax/=0
  
  write (f_res,"('------------------------------')")
  !
  write (f_res,"(a80)") (title(i), i=1,4)
  do i=1,nparam_eq
     write (f_res,"(a10,27x,i4,2x,f18.8)") parnam(i),ivar(i),param(i)
  enddo 
  do i=nparam_eq+1,parmax
       if (abs(param(i))>1e8) then
         write (f_res,"(a2,18i3,i4,2x,g19.9)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       else
         write (f_res,"(a2,18i3,i4,2x,f18.8)") trim(parnam(i)),ipower(1:NDEG,i),ivar(i),param(i)
       endif
  enddo 
  
  

do nrow=1,npts
   !
    v = energy(nrow)-eps(nrow)
    write (f_res,"(7f15.6,1x,11f11.4,3f18.7,2x,e9.2)")            &
          coord(1:NDEG,nrow),                                  &
          energy(nrow)-ve,v-ve,eps(nrow),                   &
          wt(nrow)
   !
 enddo

write (f_res,"('------------------------------')")
!
!--  printing  out the fitted PES in the form of the original ab initio
!--  for refitting with the lineariyed coordinates later 

 !open  (f_new,file='new.en',status='replace' )


 !ve = param(10)
 !write (f_new,"('ve=',e20.8)") ve


 !do nrow=1,npts
 !  !
 !  v = (energy(nrow)-eps(nrow))+ZPE
 !  !energy(:) = ( energy(:)-ZPE )*conv_energy
 !  if (v<20000.) &
 !  write (f_res,"(4f9.6,5f10.4,3g18.5,2x,e8.2)")    &
 !       coord(1:9,nrow),                           &
 !       v/conv_energy,-wt(nrow)
 !enddo


 write (f_res,"(/3x,59('-')/'   |  param  |  points |   deviat    |    rms      | stability |')")
 write (f_res,"(3x,59('-'))")
 write (f_res,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e11.3,' |')") &
               numpar,nused,stadev,rms,stability

 write (f_res,"(/92('#'))")



 write (f_out,"(/3x,59('-')/'   |  param  |  points |   deviat    |    rms      | stability |')")
 write (f_out,"(3x,59('-'))")
 write (f_out,"('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e11.3,' |')") &
               numpar,nused,stadev,rms,stability

 write (f_out,"(/92('#'))")

 end subroutine do_fitting
 !
 !

 !
 recursive function potential(parmax,dF,local,param,ipower,ipar,NDEG,NBONDS,NANGS,NPARAM_EQ) result (V)

 integer,intent(in) ::           parmax,NDEG,NBONDS,NANGS
 double precision,intent(in) ::  dF(parmax),local(NDEG),param(parmax)
 integer,intent(in) ::           ipar,ipower(NDEG,parmax)
 
 double precision         ::  param_(parmax)
 !
 double precision         ::  deltax,potright,potleft,V
 double precision         ::  dF_(parmax)
 integer                  ::  ierror,i,k
 character(len=40)        ::  txt
 !
 ! parameter for the finite differncies differention 
 !
 double precision,parameter :: factordeltax=0.001
  !
  V = 0
  !
  if (ipar==0) then
    !
    do i = nparam_eq+1,parmax
      !
      !k = ind(i)
      !
      V = V + param(i)*dF(i)
      !
    enddo
    !
  elseif (ipar>nparam_eq) then
    !
    !k = ind(ipar)
    !
    V = dF(ipar)
    !
  else
    !
    deltax=factordeltax*abs(param(ipar))
    if (deltax .le. 1e-15) deltax=1e-6
    !
    param_ = param
    !
    param_(ipar) = param(ipar)+deltax
    !
    !call diff_V(nmax,local,param_(1:2),dF_)
    !
    call diff_V(ipower,local,param_(1:nparam_eq),NDEG,NBONDS,NANGS,NPARAM_EQ,PARMAX,dF_)
    !
    !do iparam = nparam_eq+1,parmax
    !  !
    !  dF_(iparam) = diff_V(ipower(1:NDEG,iparam),local,param_(1:parmax),NDEG,NBONDS,NANGS,NPARAM_EQ,PARMAX)
    !  !
    !enddo
    !
    potright=potential(parmax,dF_,local,param_,ipower,0,NDEG,NBONDS,NANGS,NPARAM_EQ)
    !
    param_(ipar) = param(ipar)-deltax
    !
    !call diff_V(nmax,local,param_(1:2),dF_)
    call diff_V(ipower,local,param_(1:nparam_eq),NDEG,NBONDS,NANGS,NPARAM_EQ,PARMAX,dF_)
    !
    !do iparam = nparam_eq+1,parmax                                                      
    !  !
    !  dF_(iparam) = diff_V(ipower(1:NDEG,iparam),local,param_(1:parmax),NDEG,NBONDS,NANGS,NPARAM_EQ,parmax)
    !  !
    !enddo
    !
    potleft=potential(parmax,dF_,local,param_,ipower,0,NDEG,NBONDS,NANGS,NPARAM_EQ)
    !
    V=(potright-potleft)/(2.d0*deltax)
    !
  endif
  !
 end function potential
 !
 function aacos(x,ierror,txt) result (f)
    !
    double precision,intent(in)   :: x
    character(len=40),intent(in),optional :: txt
    integer,intent(out)    :: ierror
    double precision              :: f,small_,pi
    !
    small_ = epsilon(1.0d0)
    pi = 4.0d0 * datan2(1.0d0,1.0d0)
    !
    ierror = 0
    !
    if ( abs(x)>1.0d0+10.0*sqrt(small_) ) then 
       !
       write (6,"('|cos(x)|>1: ',f18.8,2x,a)") x,txt
       !
       ierror = 1
       !
       !stop 'aacos - bad cosalpha'
       !
    elseif ( x>=1.0d0) then 
       f = 0.0d0
    elseif ( x<=-1.0d0) then 
       f = pi
    else 
       f = acos(x)
    endif
    !
 end function aacos
 !

  recursive subroutine diff_V(ipower,local,param,NDEG,NBONDS,NANGS,nparam_eq,parmax,f)
  ! NEED TO MODIFY FOR C2H6/MAKE GENERAL
  ! SPECIFY FOR C2H6 FOR NOW, MAKE GENERAL IN FUTURE OR LEAVE AS A BIT OF CODE
  ! WHICH NEEDS MODIFIED FOR EACH MOLECULE
  !
  implicit none
  !
  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
  integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  !
  integer, intent(in)        :: parmax
  integer,intent(in)  :: ipower(1:ndeg,parmax)
  double precision,intent(in)  :: local(ndeg),param(nparam_eq)
  double precision,intent(out)  :: f(parmax)
  double precision             :: s1,s2,s3,s4,dst(ndeg)
  double precision         ::  v,r1,r2,r3,r4,alpha12,alpha13,alpha14,alpha23,alpha24,alpha34,deg
  double precision         ::  r0,r12,r13,r14,r23,r24,r34,t12,t13,t14,t23,t24,t34
  double precision         ::  r1e,r2e,a,y1,y2,y3,y4,y5,y6,y7,y8,y9,xi(ndeg),chi(ndeg),betae,pi,b,tau12,tau13,tau23
  DOUBLE PRECISION         ::  R5,R6,R7,A1,A2,A3,A4,TAU
  integer                  ::  ierror,ioper,NDEG,NBONDS,NANGS,nparam_eq,nsym,iparam,k(1:ndeg)
  character(len=40)        ::  txt


  NSYM = 36
  !
  txt = 'diff_V'
  !
  pi = 4.0d0 * datan2(1.0d0,1.0d0)
  !
  deg=pi/180.0d0

  r1e     = param(1)
  r2e     = param(2)
  betae   = param(3)*deg
  !
  a       = param(4)
  b       = param(5)
  !
  r1      = local(1)
  r2      = local(2)
  r3      = local(3)
  r4      = local(4)
  r5      = local(5)
  r6      = local(6)
  r7      = local(7)

  f = 0
  !
  do ioper =1,nsym
    !
    !call ML_symmetry_transformation_XY3_II(ioper,xi,chi,NDEG)
    !
    call ML_symmetry_transformation_C2H6(ioper, ndeg, local, dst)
    !
    ! NOW MODIFIED TO INCLUDE TAU IN THE SYMMETRY COORDINATES  BPM
    !
    xi(1)=1.0d+00-exp(-a*(dst(1)-r1e))
    xi(2)=1.0d+00-exp(-b*(dst(2)-r2e))
    xi(3)=1.0d+00-exp(-b*(dst(3)-r2e))
    xi(4)=1.0d+00-exp(-b*(dst(4)-r2e))
    xi(5)=1.0d+00-exp(-b*(dst(5)-r2e))
    xi(6)=1.0d+00-exp(-b*(dst(6)-r2e))
    xi(7)=1.0d+00-exp(-b*(dst(7)-r2e))
    !
    xi(8)  = dst(8)  - betae
    xi(9)  = dst(9)  - betae
    xi(10) = dst(10) - betae
    xi(11) = dst(11) - betae
    xi(12) = dst(12) - betae
    xi(13) = dst(13) - betae
    !
    xi(14) = dst(14)
    xi(15) = dst(15)
    xi(16) = dst(16)
    xi(17) = dst(17) 
    !
    !  XI(18) = LOCAL(18) - pi
    xi(18) = 1.d0 + cos(3.d0*dst(18) )
    !
    do iparam = nparam_eq+1,parmax
      !
      !dV(nrow,iparam) = diff_V(ipower(1:NDEG,iparam),local,param(1:parmax),NDEG,NBONDS,NANGS,nparam_eq,parmax)
      !
      k = ipower(:,iparam)
      !
      f(iparam) = f(iparam) + product(xi(1:18)**k(1:18))
      !
    enddo
    !
  enddo
  !
  f(:) = f(:)/real(nsym,8)

end subroutine diff_V

  !
  !
  function find_alpha34(alpha12,alpha13,alpha23,alpha14,alpha24,ierror) result (s10_f)
    !
    implicit none
    !
    real(8),intent(in)  ::alpha12,alpha13,alpha14,alpha23,alpha24
    integer,intent(out) :: ierror
    real(8) s10_f
    real(8) xtmp(4),ytmp(4),ztmp(4),temp

     ierror = 0
     !
     s10_f = 0
     !
     xtmp(1)=1.0d0
     ytmp(1)=0
     ztmp(1)=0
     xtmp(2)=cos(alpha12)
     ytmp(2)=sin(alpha12)
     ztmp(2)=0


     xtmp(3)=cos(alpha13)
     ytmp(3)=(cos(alpha23)-xtmp(2)*xtmp(3))/ytmp(2)
     !
     temp = 1.0d0-xtmp(3)**2-ytmp(3)**2
     !
     if (temp<0) then 
       !
       ierror = 1 ; return 
       !
     endif
     !
     ztmp(3)=sqrt(temp)
     !
     xtmp(4)=cos(alpha14)
     ytmp(4)=(cos(alpha24)-xtmp(2)*xtmp(4))/ytmp(2)
     !
     temp = 1.0d0-xtmp(4)**2-ytmp(4)**2
     !
     if (temp<0) then 
       !
       ierror = 1 ; return 
       !
     endif
     !
     ztmp(4)=-sqrt(temp)
     !
     temp  = xtmp(3)*xtmp(4)+ytmp(3)*ytmp(4)+ztmp(3)*ztmp(4)
     !
     if (abs(temp)>1.d00) then 
       !
       ierror = 1 ; return 
       !
     endif

     s10_f=acos(temp)

 end function find_alpha34


 function calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24) result (alpha34)
   !
   implicit none
   !
   real(8),intent(in)  ::alpha12,alpha13,alpha14,alpha23,alpha24
   real(8) :: alpha34,ca12,ca13,ca23,ca14,ca24
   real(8) :: cosalpha,small_

       small_ = epsilon(1.0d0)
       ca12 = cos(alpha12)
       ca13 = cos(alpha13)
       ca14 = cos(alpha14)
       ca23 = cos(alpha23)
       ca24 = cos(alpha24)

       cosalpha = 0.5d0*( &
                  2.0d0*ca12*ca23*ca14-     &
                  2.0d0*ca13*ca14+     &
                  2.0d0*ca12*ca24*ca13-     &
                  2.0d0*ca23*ca24+     &
                  2.0d0*sqrt(     &
                  1.0d0-ca24**2+     &
                  2.0d0*ca12*ca24*ca14+     &
                  2.0d0*ca12*ca23*ca13-     &
                  2.0d0*ca12**3*ca23*ca13-     &
                  2.0d0*ca12**2-     &
                  ca14**2-     &
                  ca13**2+     &
                  ca14**2*ca23**2+     &
                  ca13**2*ca24**2+     &
                  4.0d0*ca12**2*ca23*ca14*ca24*ca13+     &
                  ca12**4-     &
                  2.0d0*ca12*ca23*ca14**2*ca13-     &
                  2.0d0*ca12*ca23**2*ca14*ca24-     &
                  2.0d0*ca13**2*ca14*ca12*ca24-     &
                  2.0d0*ca12*ca24**2*ca13*ca23+     &
                  ca13**2*ca14**2+     &
                  ca23**2*ca24**2+     &
                  ca12**2*ca13**2+     &
                  ca12**2*ca14**2+     &
                  ca12**2*ca23**2+     &
                  ca12**2*ca24**2-     &
                  2.0d0*ca12**3*ca24*ca14-     &
                  ca23**2))/( -sin(alpha12)**2 )

       if ( abs(cosalpha)>1.0d0+sqrt(small_) ) then 
          !
          write (6,"('calc_alpha34: cosalpha>1: ',f18.8)") cosalpha
          stop 'calc_alpha34 - bad cosalpha'
          !
       elseif ( cosalpha>=1.0d0) then 
          alpha34 = 0.0d0
       else 
          alpha34 = acos(cosalpha)
       endif
       !

 end function calc_alpha34


!
!   Skip n lines  in the input file 
!
  subroutine skiplines( inpunit,n )
  integer,intent(in) :: n,inpunit
  character(len=80)  :: label
  integer            :: i0

    do i0=1,n
       read  (inpunit,"(a80)") label 
    enddo

  end subroutine skiplines


 subroutine linur(dimen,npar,coeff,constant,solution,error)
   !
   implicit none

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  double precision,intent(in)  :: coeff(npar,npar),constant(npar)
  double precision,intent(out) :: solution(npar)
  double precision          :: a0(npar,npar)
  double precision          :: c
  integer                   :: i1,i2,i,k8,k,l,k9

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
         c=0.0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0.0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0.0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine linur

!------------------------------------------!
  subroutine invmat(al,ai,dimen,npar)
  implicit none
  integer,intent(in)           :: npar,dimen
  double precision,intent(in)  :: al(npar,npar)
  double precision,intent(out) :: ai(npar,npar)
  double precision             :: h(npar),p,q
  integer                      :: i1,i2,k,i,j,k8,k9
      

    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
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
 end subroutine invmat
!------------------------------------------!



  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  subroutine ML_symmetry_transformation_XY3(ioper,src,dst)
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    double precision,intent(in)      :: src(1:9)
    double precision,intent(out)     :: dst(1:9)
    !
    double precision         :: repres(12,2,2),a,b,e,o
    !
    if (verbose>=5) write(f_out,"('ML_symmetry_transformation_XY3/start')") 
    !
    dst(1) = src(1)
    !
    select case(ioper)
    !
    case (1) ! identity 

      dst = src

    case (3) ! (123)

      dst(2) = src(3)
      dst(3) = src(4)
      dst(4) = src(2)
      !
      dst(5) = src(6)
      dst(6) = src(7)
      dst(7) = src(5)

    case (2) ! (321)

      dst(2) = src(4)
      dst(3) = src(2)
      dst(4) = src(3)
      !
      dst(5) = src(7)
      dst(6) = src(5)
      dst(7) = src(6)
 
    case (6) ! (12)

      dst(2) = src(3)
      dst(3) = src(2)
      dst(4) = src(4)
      !
      dst(5) = src(6)
      dst(6) = src(5)
      dst(7) = src(7)

    case (5) ! (13)

      dst(2) = src(4)
      dst(3) = src(3)
      dst(4) = src(2)
      !
      dst(5) = src(7)
      dst(6) = src(6)
      dst(7) = src(5)

    case (4) ! (23)

      dst(2) = src(2)
      dst(3) = src(4)
      dst(4) = src(3)
      !
      dst(5) = src(5)
      dst(6) = src(7)
      dst(7) = src(6)

    case default

      write (f_out,"('ML_symmetry_transformation_XY3: operation ',i8,' unknown')") ioper
      stop 'ML_symmetry_transformation_XY3 - bad operation. type'
 
    end select 
    !
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
    !
    repres ( 1,:,:)= reshape((/ e, o,  & 
                                o, e/),(/2,2/))
    !
    repres ( 3,:,:)= reshape((/-a,-b,  &
                                b,-a/),(/2,2/))
    !
    repres ( 2,:,:)= reshape((/-a, b,  &
                               -b,-a/),(/2,2/))
    !
    repres ( 4,:,:)= reshape((/ e, o,  &
                                o,-e/),(/2,2/))
    !
    repres ( 6,:,:)= reshape((/-a, b,  &
                                b, a/),(/2,2/))
    !
    repres ( 5,:,:)= reshape((/-a,-b,  &
                               -b, a/),(/2,2/))
    !
    dst(8) = repres(ioper,1,1)*src(8)+repres(ioper,1,2)*src(9)
    dst(9) = repres(ioper,2,1)*src(8)+repres(ioper,2,2)*src(9)
           
  end subroutine ML_symmetry_transformation_XY3        
  

  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  subroutine ML_symmetry_transformation_XY3_II(ioper,src,dst,NDEG)
!       THIS COULD ALSO BE MADE AS GENERAL AS POSSIBLE BUT WILL ALWAYS NEED TO
!       DEFINE MATRICES FOR EACH MOLECULE BPM
    implicit none 
    !
    integer,intent(in)    :: ioper,NDEG  ! group operation  
    double precision,intent(in)      :: src(1:NDEG)
    double precision,intent(out)     :: dst(1:NDEG)
    !
    double precision         :: repres(12,NDEG,NDEG),a,b,e,o
    INTEGER NSYM
    !
    if (verbose>=5) write(f_out,"('ML_symmetry_transformation_XY3/start')") 
   !
   a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
   NSYM = 12
   !
   repres = 0 
   !
   ! E
   repres(1,1,1) = 1.0d0
   repres(1,2,2) = 1.0d0
   repres(1,3,3) = 1.0d0
   repres(1,4,4) = 1.0d0
   repres(1,5,5) = 1.0d0
   repres(1,6,6) = 1.0d0
   repres(1,7,7) = 1.0d0
   repres(1,8,8) = 1.0d0
   repres(1,9,9) = 1.0d0
   REPRES(1,10,10) = 1.0d0
   REPRES(1,11,11) = 1.0d0
   REPRES(1,12,12) = 1.0d0
   REPRES(1,13,13) = 1.0d0
   REPRES(1,14,14) = 1.0d0
   REPRES(1,15,15) = 1.0d0
   REPRES(1,16,16) = 1.0d0
   REPRES(1,17,17) = 1.0d0
   REPRES(1,18,18)  = 1.D0
   !
   !C3+/(123)(465)
   repres(2,1,1) = 1.0d0
   repres(2,4,2) = 1.0d0
   repres(2,2,3) = 1.0d0
   repres(2,3,4) = 1.0d0
   repres(2,6,5) = 1.0d0
   repres(2,7,6) = 1.0d0
   repres(2,5,7) = 1.0d0
   repres(2,10,8) = 1.0d0
   repres(2,8,9) = 1.0d0
   repres(2,9,10) =  1.0d0
   repres(2,12,11) = 1.0d0
   repres(2,13,12) = 1.0d0
   repres(2,11,13) = 1.0d0
   repres(2,14,14) = -a
   repres(2,15,14) = b
   repres(2,14,15) = -b
   repres(2,15,15) = -a
   repres(2,16,16) = -a
   repres(2,17,16) = b
   repres(2,16,17) = -b
   repres(2,17,17) = -a
   REPRES(2,18,18)  = 1.D0
   !
   !C3-/(123)(465)
   repres(3,1,1) = 1.0d0
   repres(3,3,2) = 1.0d0
   repres(3,4,3) = 1.0d0
   repres(3,2,4) = 1.0d0
   repres(3,7,5) = 1.0d0
   repres(3,5,6) = 1.0d0
   repres(3,6,7) = 1.0d0
   repres(3,9,8) = 1.0d0
   repres(3,10,9) = 1.0d0
   repres(3,8,10) =  1.0d0
   repres(3,13,11) = 1.0d0
   repres(3,11,12) = 1.0d0
   repres(3,12,13) = 1.0d0
   repres(3,14,14) = -a
   repres(3,15,14) = -b
   repres(3,14,15) = b
   repres(3,15,15) = -a
   repres(3,16,16) = -a
   repres(3,17,16) = -b
   repres(3,16,17) = b
   repres(3,17,17) = -a
   REPRES(3,18,18)  = 1.D0
   !
   !C2/(16)(24)(35)(78)
   repres(4,1,1) = 1.0d0
   repres(4,7,2) = 1.0d0
   repres(4,5,3) = 1.0d0
   repres(4,6,4) = 1.0d0
   repres(4,3,5) = 1.0d0
   repres(4,4,6) = 1.0d0
   repres(4,2,7) = 1.0d0
   repres(4,13,8) = 1.0d0
   repres(4,11,9) = 1.0d0
   repres(4,12,10) =  1.0d0
   repres(4,9,11) = 1.0d0
   repres(4,10,12) = 1.0d0
   repres(4,8,13) = 1.0d0
   repres(4,16,14) = 1.d0
   repres(4,17,15) = -1.d0
   repres(4,14,16) = 1.d0
   repres(4,15,17) = -1.d0
   REPRES(4,18,18)  = 1.D0
   !
   !C2'/(16)(24)(35)(78)
   repres(5,1,1) = 1.0d0
   repres(5,6,2) = 1.0d0
   repres(5,7,3) = 1.0d0
   repres(5,5,4) = 1.0d0
   repres(5,4,5) = 1.0d0
   repres(5,2,6) = 1.0d0
   repres(5,3,7) = 1.0d0
   repres(5,12,8) = 1.0d0
   repres(5,13,9) = 1.0d0
   repres(5,11,10) =  1.0d0
   repres(5,10,11) = 1.0d0
   repres(5,8,12) = 1.0d0
   repres(5,9,13) = 1.0d0
   repres(5,16,14) = -a
   repres(5,17,14) = -b
   repres(5,16,15) = b
   repres(5,17,15) = a
   repres(5,14,16) = -a
   repres(5,15,16) = -b
   repres(5,14,17) = -b
   repres(5,15,17) =  a
   REPRES(5,18,18)  = 1.D0
   !
   !C2''/(16)(24)(35)(78)
   repres(6,1,1) = 1.0d0
   repres(6,5,2) = 1.0d0
   repres(6,6,3) = 1.0d0
   repres(6,7,4) = 1.0d0
   repres(6,2,5) = 1.0d0
   repres(6,3,6) = 1.0d0
   repres(6,4,7) = 1.0d0
   repres(6,11,8) = 1.0d0
   repres(6,12,9) = 1.0d0
   repres(6,13,10) =  1.0d0
   repres(6,8,11) = 1.0d0
   repres(6,9,12) = 1.0d0
   repres(6,10,13) = 1.0d0
   repres(6,16,14) = -a
   repres(6,17,14) = b
   repres(6,16,15) = b
   repres(6,17,15) = a
   repres(6,14,16) = -a
   repres(6,14,16) = -a
   repres(6,15,16) = b
   repres(6,14,17) = b
   repres(6,15,17) = a
   REPRES(6,18,18)  = 1.D0
   !
   !i/(14)(26)(35)(78)*
   repres(7,1,1) = 1.0d0
   repres(7,5,2) = 1.0d0
   repres(7,7,3) = 1.0d0
   repres(7,6,4) = 1.0d0
   repres(7,2,5) = 1.0d0
   repres(7,4,6) = 1.0d0
   repres(7,3,7) = 1.0d0
   repres(7,11,8) = 1.0d0
   repres(7,13,9) = 1.0d0
   repres(7,12,10) =  1.0d0
   repres(7,8,11) = 1.0d0
   repres(7,10,12) = 1.0d0
   repres(7,9,13) = 1.0d0
   repres(7,16,14) = 1.d0
   repres(7,17,15) = 1.d0
   repres(7,14,16) = 1.d0
   repres(7,15,17) = 1.d0
   REPRES(7,18,18)  = 1.D0
   !
   !S6/(163425)(78)*
   repres(8,1,1) = 1.0d0
   repres(8,6,2) = 1.0d0
   repres(8,5,3) = 1.0d0
   repres(8,7,4) = 1.0d0
   repres(8,4,5) = 1.0d0
   repres(8,3,6) = 1.0d0
   repres(8,2,7) = 1.0d0
   repres(8,12,8) = 1.0d0
   repres(8,11,9) = 1.0d0
   repres(8,13,10) =  1.0d0
   repres(8,10,11) = 1.0d0
   repres(8,9,12) = 1.0d0
   repres(8,8,13) = 1.0d0
   repres(8,16,14) = -a
   repres(8,17,14) = b
   repres(8,16,15) = -b
   repres(8,17,15) = -a
   repres(8,14,16) = -a
   repres(8,15,16) = b
   repres(8,14,17) = -b
   repres(8,15,17) = -a
   REPRES(8,18,18)  = 1.D0
   !
   !S6'/(14)(26)(35)(78)*
   repres(9,1,1) = 1.0d0
   repres(9,7,2) = 1.0d0
   repres(9,6,3) = 1.0d0
   repres(9,5,4) = 1.0d0
   repres(9,3,5) = 1.0d0
   repres(9,2,6) = 1.0d0
   repres(9,4,7) = 1.0d0
   repres(9,13,8) = 1.0d0
   repres(9,12,9) = 1.0d0
   repres(9,11,10) =  1.0d0
   repres(9,9,11) = 1.0d0
   repres(9,8,12) = 1.0d0
   repres(9,10,13) = 1.0d0
   repres(9,16,14) = -a
   repres(9,17,14) = -b
   repres(9,16,15) = b
   repres(9,17,15) = -a
   repres(9,14,16) = -a
   repres(9,15,16) = -b
   repres(9,14,17) = b
   repres(9,15,17) = -a    
   REPRES(9,18,18)  = 1.D0
   !
   !sigmad/(12)(46)*
   repres(10,1,1) = 1.0d0
   repres(10,4,2) = 1.0d0
   repres(10,3,3) = 1.0d0
   repres(10,2,4) = 1.0d0
   repres(10,6,5) = 1.0d0
   repres(10,5,6) = 1.0d0
   repres(10,7,7) = 1.0d0
   repres(10,10,8) = 1.0d0
   repres(10,9,9) = 1.0d0
   repres(10,8,10) =  1.0d0
   repres(10,12,11) = 1.0d0
   repres(10,11,12) = 1.0d0
   repres(10,13,13) = 1.0d0
   repres(10,14,14) = -a
   repres(10,15,14) = -b
   repres(10,14,15) = -b
   repres(10,15,15) = a
   repres(10,16,16) = -a
   repres(10,17,16) = -b
   repres(10,16,17) = -b
   repres(10,17,17) = a
   REPRES(10,18,18)  = 1.D0
   !
   !sigmad'/(12)(46)*
   repres(11,1,1) = 1.0d0
   repres(11,2,2) = 1.0d0
   repres(11,4,3) = 1.0d0
   repres(11,3,4) = 1.0d0
   repres(11,5,5) = 1.0d0
   repres(11,7,6) = 1.0d0
   repres(11,6,7) = 1.0d0
   repres(11,8,8) = 1.0d0
   repres(11,10,9) = 1.0d0
   repres(11,9,10) =  1.0d0
   repres(11,11,11) = 1.0d0
   repres(11,13,12) = 1.0d0
   repres(11,12,13) = 1.0d0
   repres(11,14,14) =  -a
   repres(11,15,14) =  b
   repres(11,14,15) =  b
   repres(11,15,15) =  a
   repres(11,16,16) = -a
   repres(11,17,16) =  b
   repres(11,16,17) =  b
   repres(11,17,17) =  a
   REPRES(11,18,18)  = 1.D0
   !
   !sigmad''/(12)(46)*
   repres(12,1,1) = 1.0d0
   repres(12,3,2) = 1.0d0
   repres(12,2,3) = 1.0d0
   repres(12,4,4) = 1.0d0
   repres(12,7,5) = 1.0d0
   repres(12,6,6) = 1.0d0
   repres(12,5,7) = 1.0d0
   repres(12,9,8) = 1.0d0
   repres(12,8,9) = 1.0d0
   repres(12,10,10) =  1.0d0
   repres(12,13,11) = 1.0d0
   repres(12,12,12) = 1.0d0
    repres(12,11,13) = 1.0d0
    repres(12,14,14) = 1.d0
    repres(12,15,15) =  -1.d0
    repres(12,16,16) = 1.d0
    repres(12,17,17) =  -1.d0
    REPRES(12,18,18)  = 1.D0
    !
    if (ioper<0.or.ioper>NSYM) then
      write (6,"('symmetry_transformation_local: operation ',i8,' unknown')") ioper
      stop 'symmetry_transformation_local - bad operation. type'
    endif
    !
    ! BPM : REVERSED MULTIPLICATION ORDER TO AGREE WITH MY NOTES
    dst = matmul(src,repres(ioper,:,:)) 
           
  end subroutine ML_symmetry_transformation_XY3_II      
  
  
 ! subroutine ML_symmetry_transformation_XY3_II(ioper,src,dst,NDEG)
!!       THIS COULD ALSO BE MADE AS GENERAL AS POSSIBLE BUT WILL ALWAYS NEED TO
!!       DEFINE MATRICES FOR EACH MOLECULE BPM
 !
!!       THIS IS THE MODIFIED MATRICES FOR SERGEY'S COORDINATES 
 !   implicit none 
 !   !
 !   integer,intent(in)    :: ioper,NDEG  ! group operation  
 !   double precision,intent(in)      :: src(1:NDEG)
 !   double precision,intent(out)     :: dst(1:NDEG)
 !   !
 !   double precision         :: repres(12,NDEG,NDEG),a,b,e,o
 !   INTEGER NSYM
 !   !
 !   if (verbose>=5) write(f_out,"('ML_symmetry_transformation_XY3/start')") 
 !   !
 !   a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
 !   NSYM = 12
 !   !
 !   repres = 0 
 !   !
 !   ! E
 !   repres(1,1,1) = 1.0d0
 !   repres(1,2,2) = 1.0d0
 !   repres(1,3,3) = 1.0d0
 !   repres(1,4,4) = 1.0d0
 !   repres(1,5,5) = 1.0d0
 !   repres(1,6,6) = 1.0d0
 !   repres(1,7,7) = 1.0d0
 !   repres(1,8,8) = 1.0d0
 !   repres(1,9,9) = 1.0d0
 !   REPRES(1,10,10) = 1.0d0
 !   REPRES(1,11,11) = 1.0d0
 !   REPRES(1,12,12) = 1.0d0
 !   REPRES(1,13,13) = 1.0d0
 !   REPRES(1,14,14) = 1.0d0
 !   REPRES(1,15,15) = 1.0d0
 !   REPRES(1,16,16) = 1.0d0
 !   REPRES(1,17,17) = 1.0d0
 !   REPRES(1,18,18)  = 1.D0
 !   !
 !   !C3+/(123)(465)
 !   repres(2,1,1) = 1.0d0
 !   repres(2,4,2) = 1.0d0
 !   repres(2,2,3) = 1.0d0
 !   repres(2,3,4) = 1.0d0
 !   repres(2,7,5) = 1.0d0
 !   repres(2,5,6) = 1.0d0
 !   repres(2,6,7) = 1.0d0
 !   repres(2,10,8) = 1.0d0
 !   repres(2,8,9) = 1.0d0
 !   repres(2,9,10) =  1.0d0
 !   repres(2,13,11) = 1.0d0
 !   repres(2,11,12) = 1.0d0
 !   repres(2,12,13) = 1.0d0
 !   repres(2,14,14) = -a
 !   repres(2,15,14) = -b
 !   repres(2,14,15) =  b
 !   repres(2,15,15) = -a
 !   repres(2,16,16) = -a
 !   repres(2,17,16) = -b
 !   repres(2,16,17) =  b
 !   repres(2,17,17) = -a
 !   REPRES(2,18,18)  = 1.D0
 !   !
 !   !C3-/(123)(465)
 !   repres(3,1,1) = 1.0d0
 !   repres(3,3,2) = 1.0d0
 !   repres(3,4,3) = 1.0d0
 !   repres(3,2,4) = 1.0d0
 !   repres(3,6,5) = 1.0d0
 !   repres(3,7,6) = 1.0d0
 !   repres(3,5,7) = 1.0d0
 !   repres(3,9,8) = 1.0d0
 !   repres(3,10,9) = 1.0d0
 !   repres(3,8,10) =  1.0d0
 !   repres(3,12,11) = 1.0d0
 !   repres(3,13,12) = 1.0d0
 !   repres(3,11,13) = 1.0d0
 !   repres(3,14,14) = -a
 !   repres(3,15,14) =  b
 !   repres(3,14,15) = -b
 !   repres(3,15,15) = -a
 !   repres(3,16,16) = -a
 !   repres(3,17,16) =  b
 !   repres(3,16,17) = -b
 !   repres(3,17,17) = -a
 !   REPRES(3,18,18)  = 1.D0
 !   !
 !   !C2/(16)(24)(35)(78)
 !   repres(4,1,1) = 1.0d0
 !   repres(4,6,2) = 1.0d0
 !   repres(4,5,3) = 1.0d0
 !   repres(4,7,4) = 1.0d0
 !   repres(4,3,5) = 1.0d0
 !   repres(4,2,6) = 1.0d0
 !   repres(4,4,7) = 1.0d0
 !   repres(4,12,8) = 1.0d0
 !   repres(4,11,9) = 1.0d0
 !   repres(4,13,10) =  1.0d0
 !   repres(4,9,11) = 1.0d0
 !   repres(4,8,12) = 1.0d0
 !   repres(4,10,13) = 1.0d0
 !   repres(4,16,14) = -a
 !   repres(4,17,14) =  b
 !   repres(4,16,15) =  b
 !   repres(4,17,15) =  a
 !   repres(4,14,16) = -a
 !   repres(4,15,16) =  b
 !   repres(4,14,17) =  b
 !   repres(4,15,17) =  a
 !   REPRES(4,18,18)  = 1.D0
 !   !
 !   !C2'/(16)(24)(35)(78)
 !   repres(5,1,1) = 1.0d0
 !   repres(5,5,2) = 1.0d0
 !   repres(5,7,3) = 1.0d0
 !   repres(5,6,4) = 1.0d0
 !   repres(5,2,5) = 1.0d0
 !   repres(5,4,6) = 1.0d0
 !   repres(5,3,7) = 1.0d0
 !   repres(5,11,8) = 1.0d0
 !   repres(5,13,9) = 1.0d0
 !   repres(5,12,10) =  1.0d0
 !   repres(5,8,11) = 1.0d0
 !   repres(5,10,12) = 1.0d0
 !   repres(5,9,13) = 1.0d0
 !   repres(5,16,14) = 1.d0
 !   repres(5,17,15) = -1.d0
 !   repres(5,14,16) = 1.d0
 !   repres(5,15,17) = -1.d0
 !   REPRES(5,18,18)  = 1.D0
 !   !
 !   !C2''/(16)(24)(35)(78)
 !   repres(6,1,1) = 1.0d0
 !   repres(6,7,2) = 1.0d0
 !   repres(6,6,3) = 1.0d0
 !   repres(6,5,4) = 1.0d0
 !   repres(6,4,5) = 1.0d0
 !   repres(6,3,6) = 1.0d0
 !   repres(6,2,7) = 1.0d0
 !   repres(6,13,8) = 1.0d0
 !   repres(6,12,9) = 1.0d0
 !   repres(6,11,10) =  1.0d0
 !   repres(6,10,11) = 1.0d0
 !   repres(6,9,12) = 1.0d0
 !   repres(6,8,13) = 1.0d0
 !   repres(6,16,14) = -a
 !   repres(6,17,14) = -b
 !   repres(6,16,15) = -b
 !   repres(6,17,15) =  a
 !   repres(6,14,16) = -a
 !   repres(6,15,16) = -b
 !   repres(6,14,17) = -b
 !   repres(6,15,17) =  a
 !   REPRES(6,18,18)  = 1.D0
 !   !
 !   !i/(14)(26)(35)(78)*
 !   repres(7,1,1) = 1.0d0
 !   repres(7,5,2) = 1.0d0
 !   repres(7,6,3) = 1.0d0
 !   repres(7,7,4) = 1.0d0
 !   repres(7,2,5) = 1.0d0
 !   repres(7,3,6) = 1.0d0
 !   repres(7,4,7) = 1.0d0
 !   repres(7,11,8) = 1.0d0
 !   repres(7,12,9) = 1.0d0
 !   repres(7,13,10) =  1.0d0
 !   repres(7,8,11) = 1.0d0
 !   repres(7,9,12) = 1.0d0
 !   repres(7,10,13) = 1.0d0
 !   repres(7,16,14) = 1.d0
 !   repres(7,17,15) = 1.d0
 !   repres(7,14,16) = 1.d0
 !   repres(7,15,17) = 1.d0
 !   REPRES(7,18,18)  = 1.D0
 !   !
 !   !S6/(163425)(78)*
 !   repres(8,1,1) = 1.0d0
 !   repres(8,6,2) = 1.0d0
 !   repres(8,7,3) = 1.0d0
 !   repres(8,5,4) = 1.0d0
 !   repres(8,3,5) = 1.0d0
 !   repres(8,4,6) = 1.0d0
 !   repres(8,2,7) = 1.0d0
 !   repres(8,12,8) = 1.0d0
 !   repres(8,13,9) = 1.0d0
 !   repres(8,11,10) =  1.0d0
 !   repres(8,9,11) = 1.0d0
 !   repres(8,10,12) = 1.0d0
 !   repres(8,8,13) = 1.0d0
 !   repres(8,16,14) = -a
 !   repres(8,17,14) =  b
 !   repres(8,16,15) = -b
 !   repres(8,17,15) = -a
 !   repres(8,14,16) = -a
 !   repres(8,15,16) = b
 !   repres(8,14,17) = -b
 !   repres(8,15,17) = -a
 !   REPRES(8,18,18)  = 1.D0
 !   !
 !   !S6'/(14)(26)(35)(78)*
 !   repres(9,1,1) = 1.0d0
 !   repres(9,7,2) = 1.0d0
 !   repres(9,5,3) = 1.0d0
 !   repres(9,6,4) = 1.0d0
 !   repres(9,4,5) = 1.0d0
 !   repres(9,2,6) = 1.0d0
 !   repres(9,3,7) = 1.0d0
 !   repres(9,13,8) = 1.0d0
 !   repres(9,11,9) = 1.0d0
 !   repres(9,12,10) =  1.0d0
 !   repres(9,10,11) = 1.0d0
 !   repres(9,8,12) = 1.0d0
 !   repres(9,9,13) = 1.0d0
 !   repres(9,16,14) = -a
 !   repres(9,17,14) = -b
 !   repres(9,16,15) = b
 !   repres(9,17,15) = -a
 !   repres(9,14,16) = -a
 !   repres(9,15,16) = -b
 !   repres(9,14,17) = b
 !   repres(9,15,17) = -a    
 !   REPRES(9,18,18)  = 1.D0
 !   !
 !   !sigmad/(12)(46)*
 !   repres(10,1,1) = 1.0d0
 !   repres(10,4,2) = 1.0d0
 !   repres(10,3,3) = 1.0d0
 !   repres(10,2,4) = 1.0d0
 !   repres(10,7,5) = 1.0d0
 !   repres(10,6,6) = 1.0d0
 !   repres(10,5,7) = 1.0d0
 !   repres(10,10,8) = 1.0d0
 !   repres(10,9,9) = 1.0d0
 !   repres(10,8,10) =  1.0d0
 !   repres(10,13,11) = 1.0d0
 !   repres(10,12,12) = 1.0d0
 !   repres(10,11,13) = 1.0d0
 !   repres(10,14,14) = -a
 !   repres(10,15,14) = -b
 !   repres(10,14,15) = -b
 !   repres(10,15,15) = a
 !   repres(10,16,16) = -a
 !   repres(10,17,16) = -b
 !   repres(10,16,17) = -b
 !   repres(10,17,17) = a
 !   REPRES(10,18,18)  = 1.D0
 !   !
 !   !sigmad'/(12)(46)*
 !   repres(11,1,1) = 1.0d0
 !   repres(11,2,2) = 1.0d0
 !   repres(11,4,3) = 1.0d0
 !   repres(11,3,4) = 1.0d0
 !   repres(11,5,5) = 1.0d0
 !   repres(11,7,6) = 1.0d0
 !   repres(11,6,7) = 1.0d0
 !   repres(11,8,8) = 1.0d0
 !   repres(11,10,9) = 1.0d0
 !   repres(11,9,10) =  1.0d0
 !   repres(11,11,11) = 1.0d0
 !   repres(11,13,12) = 1.0d0
 !   repres(11,12,13) = 1.0d0
 !   repres(11,14,14) =  1.d0
 !   repres(11,15,15) =  -1.d0
 !   repres(11,16,16) = 1.d0
 !   repres(11,17,17) =  -1.d0
 !   REPRES(11,18,18)  = 1.D0
 !   !
 !   !sigmad''/(12)(46)*
 !   repres(12,1,1) = 1.0d0
 !   repres(12,3,2) = 1.0d0
 !   repres(12,2,3) = 1.0d0
 !   repres(12,4,4) = 1.0d0
 !   repres(12,6,5) = 1.0d0
 !   repres(12,5,6) = 1.0d0
 !   repres(12,7,7) = 1.0d0
 !   repres(12,9,8) = 1.0d0
 !   repres(12,8,9) = 1.0d0
 !   repres(12,10,10) =  1.0d0
 !   repres(12,12,11) = 1.0d0
 !   repres(12,11,12) = 1.0d0
 !   repres(12,13,13) = 1.0d0
 !   repres(12,14,14) = -a
 !   repres(12,15,14) =  b
 !   repres(12,14,15) =  b
 !   repres(12,15,15) =  a
 !   repres(12,16,16) = -a
 !   repres(12,17,16) =  b
 !   repres(12,16,17) =  b
 !   repres(12,17,17) =  a
 !   REPRES(12,18,18)  = 1.D0
 !   !
 !   if (ioper<0.or.ioper>NSYM) then
 !     write (6,"('symmetry_transformation_local: operation ',i8,' unknown')") ioper
 !     stop 'symmetry_transformation_local - bad operation. type'
 !   endif
 !   !
 !   ! BPM : REVERSED MULTIPLICATION ORDER TO AGREE WITH MY NOTES
 !   dst = matmul(src,repres(ioper,:,:)) 
 !          
 ! end subroutine ML_symmetry_transformation_XY3_II      
 ! 
 !
 
 
   recursive subroutine ML_symmetry_transformation_C2H6(ioper, nmodes, src, dst)
    !
    ! Symmetry transformation rules of coordinates
    !
    integer, intent(in)  ::  ioper
    integer, intent(in)  ::  nmodes
    real(8), intent(in)    ::  src(1:nmodes)
    real(8), intent(out)   ::  dst(1:nmodes)
    real(8) :: a,b,e,o,g(1:4,1:4),c123(2,2),c132(2,2),a123(3,3),a132(3,3),sxy(2,2),i(3,3),i2(3,3)
 !
    real(8),dimension(size(src)) :: tmp
    !
    integer(4)  :: tn(72,2), temp(144)
    integer(4) :: nsrc
    character(80) :: coords_transform
    !
    coords_transform = 'TAU-11'
    !
    temp(1:36)   = (/0, 0, 2, 0, 6, 4, 0, 7, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 0,21,19, 2, 3, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6/)
    temp(73:108) = (/0, 0, 2, 0, 2, 2, 0, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8, 0, 7, 7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
    !
    !temp(1:36)   = (/0, 0, 2, 0, 4, 5, 0, 7, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 0,21,19, 2, 3, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6/)
    !temp(73:108) = (/0, 0, 2, 0, 2, 2, 0, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8, 0, 7, 7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
    ! 
    temp(37:72)   = (/ 0, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37/)
    temp(109:144) = (/ 0,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36/)
    !
    tn = reshape( temp, (/ 72, 2/))
    !
    a = 0.5d0
    b = 0.5d0*sqrt(3.0d0)
    !
    e = 1.0d0
    o = 0.0d0
    !
    c123 = transpose(reshape( (/ -a, -b, &
                                  b, -a/), (/ 2, 2/)))
    c132  = matmul(c123,c123)
    !
    a123 = transpose(reshape( (/ o, o, e,&
                                 e, o, o,&
                                 o, e, o/), (/ 3, 3/)))
    a132  = matmul(a123,a123)
    !
    i = transpose(reshape( (/ e, o, o,&
                              o, e, o,&
                              o, o, e/), (/ 3, 3/)))
    !
    i2= transpose(reshape( (/ e, o, o,&
                              o, o, e,&
                              o, e, o/), (/ 3, 3/)))
    !
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))
    !
    select case(trim(coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'MLpoten_c2h6_88 error', trim(coords_transform), 'is unknown'
      stop 'ML_symmetry_transformation_C2H6 error: bad coordinate type'
      !
    case('TAU-11')    
      !
      select case(ioper)
        !
      case (1) ! E
        !
        dst(1:18) = src(1:18)
        !
      case (2) !C(+)/(123)(456)
        !
        dst(1) = src(1)
        !
        dst(2:4) = matmul(a123,src(2:4))
        dst(5:7) = matmul(a123,src(5:7))
        !
        dst(8:10) = matmul(a123,src(8:10))
        dst(11:13) = matmul(a123,src(11:13))
        !
        !! a07
        !dst(14:15) = matmul(c132,src(14:15))
        !dst(16:17) = matmul(c132,src(16:17))
        !!
        dst(14:15) = matmul(c123,src(14:15))
        dst(16:17) = matmul(c123,src(16:17))
        !
        ! a02
        !dst(14:15) = matmul(c132,src(14:15))
        !dst(16:17) = matmul(c132,src(16:17))
        !!
        !dst(18) = src(18)  + 4.0d0/3.0d0*pi
        !
        ! a05
        dst(18) = src(18)  - 4.0d0/3.0d0*pi
        !
        ! a04
        !dst(18) = src(18)
        !
        do while(dst(18) < 0.0d0) 
              dst(18) = dst(18) + 4.0d0*pi
        enddo
        !
        do while(dst(18) > 4.0d0*pi) 
              dst(18) = dst(18) - 4.0d0*pi
        enddo
        !    
      case (4) !sxy(+)/(14)(26)(35)(ab)* 
        !
        dst(1) = src(1)
        !
        dst(2:4) = matmul(i2,src(5:7))
        dst(5:7) = matmul(i2,src(2:4))
        !
        dst(8:10) = matmul(i2,src(11:13))
        dst(11:13) = matmul(i2,src(8:10))
        !
        ! a02
        !dst(14:15) = src(16:17)
        !dst(16:17) = src(14:15)
        !!
        dst(14:15) = matmul(sxy,src(16:17))
        dst(16:17) = matmul(sxy,src(14:15))
        !!
        dst(18) =  4.0d0*pi - src(18)
        !
        ! a12
        !dst(18) = -2.0d0*pi - src(18)
        !
        do while(dst(18) < 0.0d0) 
              dst(18) = dst(18) + 4.0d0*pi
        enddo
        do while(dst(18) > 4.0d0*pi) 
              dst(18) = dst(18) - 4.0d0*pi
        enddo
       !
      case (7) ! C(-)/(132)(456)
        !
        dst(1) = src(1)
        !
        dst(2:4) = matmul(a132,src(2:4))
        dst(5:7) = matmul(a123,src(5:7))
        !
        dst(8:10) = matmul(a132,src(8:10))
        dst(11:13) = matmul(a123,src(11:13))
        !!
        dst(14:15) = matmul(c132,src(14:15))
        dst(16:17) = matmul(c123,src(16:17))
        !
        !a07,a08
        !dst(14:15) = matmul(c123,src(14:15))
        !dst(16:17) = matmul(c132,src(16:17))
        !!
        dst(18) = src(18)
        ! a04
        !dst(18) = src(18)  + 4.0d0/3.0d0*pi
        !
      case (19) !sxy(-)/(14)(25)(36)(ab)
        !
        dst(1) = src(1)
        !
        dst(2:4) = matmul(i,src(5:7))
        dst(5:7) = matmul(i,src(2:4))
        !
        dst(8:10) = matmul(i,src(11:13))
        dst(11:13) = matmul(i,src(8:10))
        !!
        dst(14:15) = src(16:17)
        dst(16:17) = src(14:15)
        !
        !a02
        !dst(14:15) = matmul(sxy,src(16:17))
        !dst(16:17) = matmul(sxy,src(14:15))
        !!
        dst(18) = src(18)
        !
      case(37) !E'
         dst(1:17) = src(1:17)
         dst(18) = src(18) + 2.0d0*pi
         do while(dst(18) < 0.0d0) 
              dst(18) = dst(18) + 4.0d0*pi
         enddo
         do while(dst(18) > 4.0d0*pi) 
              dst(18) = dst(18) - 4.0d0*pi
         enddo
         !
      end select
      !
      if (all(tn(ioper,:)/=0)) then
          call ML_symmetry_transformation_C2H6(tn(ioper,1),nmodes,src,tmp)
          call ML_symmetry_transformation_C2H6(tn(ioper,2),nmodes,tmp,dst)
      endif 
      !
    end select 
    !
  end subroutine ML_symmetry_transformation_C2H6



 
  subroutine transform_coordiates_tau_11(nmodes,coords)
    !
    implicit none
    integer, intent(in) :: nmodes
    real(8), intent(inout)    ::  coords(1:nmodes)
    real(8) :: t14,t12,t13,t46,t45,t23,t56,t16,t34,t15,t43,t25,t36,t24,t35,tbar,A1,A2,A3,A4
     !
     t14 = coords(14)
     t24 = coords(15)
     t25 = coords(16) ! negative
     t35 = coords(17)
     t36 = coords(18) ! negative 
     !
     ! don't swap H 5 and 6!
     !t45 =  dabs(coord(18,i))
     !t46 =  coord(17,i)
     !
     t12 = mod(t14-t24+360.0d0,360.0d0)
     t23 = mod(t25-t35+360.0d0,360.0d0)
     t13 = mod(360.0d0-t12-t23+360.0d0,360.0d0)
     !
     t56 = mod(t36-t35+360.0d0,360.0d0)
     t45 = mod(t25-t24+360.0d0,360.0d0)
     t46 = mod(360.0d0-t56-t45+360.0d0,360.0d0)
     !
     tbar = (t14 + t25 + t36)/3.0d0
     !
     ! the coordinate chouce corresponds to tau-11 in trove 
     A1 = 1.0d0/sqrt(6.0d0)*(2.0d0*T23 - T13 - T12)
     A2 = 1.0d0/sqrt(2.0d0)*(T13 - T12)
     A3 = 1.0d0/sqrt(6.0d0)*(2.d0*T56 - T46 - T45)
     A4 = 1.0d0/sqrt(2.0d0)*(T46 - T45) 
     !
     !
     coords(14) = A1
     coords(15) = A2
     coords(16) = A3
     coords(17) = A4
     coords(18) = tbar
     !
   
   end subroutine transform_coordiates_tau_11
 
 
 
  
  end module FITPES_C2H6
  
  
  
  
  
  
  program driver
  use FITPES_C2H6
  
  call do_fitting
  
  end program driver

module me_numer
! 
! Numerical matrix elements based on the Numerov procedure
!
  use accuracy
  use molecules
  use moltype
  use lapack
  use timer
  implicit none
  private  

  public ik, rk, out
  public ME_numerov,simpsonintegral,simpsonintegral_ark,integral_rect_ark

  !
  real(ark), parameter :: enermax=  700000.0_rk   !  Upper limit for the energy search (1/cm)
  real(ark), parameter :: enerdev=  0.1_rk        !  In the i-th solution has been already defined at 
                                                  !  previouse steps, we solve the eigen-problem anyway
                                                  !  but within +/- "enerdev" (cm-1) interval from 
                                                  !  the found solution
  real(ark), parameter :: enerstep = 100.0_ark    !  Energy interval needed for the step-wise search of the 
                                                  !  solution, when the iterative procedure does not help 
                                                  !
  real(ark), parameter :: enersmall= 0.000001_ark !  Small energy value to check whether the search interval not collapsed 
  !
  real(ark), parameter :: thrsh_int = 10.0_ark**(-ark) ! Threshold in Numerov-integration 
  real(ark), parameter :: thrsh_dif = 0.1e-2_ark       ! Threshold in numerical differentiation
  real(ark), parameter :: thrsh_upper  = 1.0e20_ark    ! Upper limit for fncts in out-numerov integr.
  !
  integer(ik), parameter :: itermax= 500000      ! Maximal possible number of iterations (tries) for Numerov 
  integer(ik), parameter :: maxslots = 500       ! Maximal possible number of slots for found energies 
  integer(ik), parameter :: epoints = 40         ! N of points used for extrapolation in lsq_fit
  !
  real(ark) :: rhostep                           ! mesh size
  !
  real(ark),parameter :: wave_init=1.0_ark       ! initial arbitrary value at rho = rho_ref
  !
  !
  integer(ik), parameter :: iterlimit=5000       ! Iteration limit in numerov integration 
  !
  integer(ik) :: npoints                         ! number of grid points 
  real(ark)   :: rho_b(2)                        ! rhomin..rhomax
  integer(ik) :: isingular                       ! point with singularity, none if <0
  integer(ik) :: iperiod                         ! the periodicity (can be negative for the reflecation)
  !
  integer(ik) :: verbose = 6                     ! Verbosity level
  !
  logical     :: periodic                        ! periodic boundary conditions: true or false
  !
  integer(ik) :: vmax,maxorder,imin
  !
  integer(ik) :: io_slot                         ! unit numeber to store the numerov eigenvectors and their derivatives
  !
  integer(ik) :: imode                           ! the current mode
  !
  integer(ik) :: Nr = 4                          ! 2*Nr+1 is the number of interpolation points
  !
  integer(ik) :: iparity                         ! 
  !
  character(len=cl),parameter :: deriv_method = 'ML_diffs' !  We use this method to estimate d pvi_v / d rho  
                                                         ! where phi_v is a numerical eigenfunction from numerov
                                                         ! deriv_method can be either 'd04aaf' of '5 points'
                                                         ! '5 points' 
  real(ark) ::  rho_switch  = .0174532925199432957692369_ark       ! the value of abcisse rho of the switch between regions (1 deg)

  integer(ik) :: iswitch                                 ! the grid point of switch

  !
  contains

  !
  ! Matrix elements with Numerov-eigenfunctions 
  !
  subroutine ME_numerov(vmax_,maxorder_,rho_b_,isingular_,npoints_,numerpoints_,drho_,xton_,poten_,mu_rr_,icoord,iperiod_,&
                        verbose_,g_numerov,energy)
   !
   integer(ik),intent(in)   :: vmax_,maxorder_,npoints_,isingular_,numerpoints_,iperiod_
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder_,0:vmax_,0:vmax_)
   real(ark),intent(out)    :: energy(0:vmax_)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten_(0:npoints_),mu_rr_(0:npoints_),drho_(0:npoints_,3),xton_(0:npoints_,0:maxorder_)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose_   ! Verbosity level
   !
   real(ark)            :: rho,cross_prod,factor
   real(ark)            :: h_t,sigma,sigma_t,rms,psipsi_t,characvalue,rhostep_,step_scale,fval,df_t
   !
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,k,i_,i1,i2
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:)
   real(ark),allocatable :: f(:),poten(:),mu_rr(:),d2fdr2(:),dfdr(:),rho_(:),xton(:,:)
   character(len=cl)     :: unitfname 
   real(ark),allocatable :: enerslot(:),enerslot_(:)
   real(ark),allocatable :: psi(:,:),dpsi(:,:)
    !
    if (verbose>=1) write (out,"(/'Numerov matrix elements calculations')")
     !
     ! global variables
     !
     vmax      = vmax_
     maxorder  = maxorder_
     npoints   = numerpoints_
     rho_b = rho_b_
     isingular = isingular_
     iperiod = iperiod_
     verbose = verbose_
     imode=icoord
     rho_switch = molec%specparam(icoord)
     !
     periodic = .false.
     if (iperiod>0) periodic = .true.
     !
     allocate(phil(0:npoints_),phir(0:npoints_),dphil(0:npoints_),dphir(0:npoints_), &
              phivphi(0:npoints_),rho_kinet(0:npoints_),rho_poten(0:npoints_),rho_extF(0:npoints_),enerslot(0:maxslots), &
              f(0:npoints),dfdr(0:npoints),d2fdr2(0:npoints),poten(0:npoints),mu_rr(0:npoints),&
              xton(0:npoints,0:maxorder_),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     enerslot = 0 
     characvalue = 50000.0
     !
     ! numerov step size 
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
     !
     ! real step size 
     rhostep_ = (rho_b(2)-rho_b(1))/real(npoints_,kind=ark)
     !
     ! define the rho-type coordinate 
     !
     rho_kinet(:) = drho_(:,1)
     rho_poten(:) = drho_(:,2)
     rho_extF(:)  = drho_(:,3)
     !
     step_scale = rhostep/rhostep_
     !
     ! interpolation 
     !
     if (npoints==npoints_) then 
       !
       poten = poten_
       mu_rr = mu_rr_
       !
       do lambda  = 0,maxorder_
         !
         xton(:,lambda) = xton_(:,lambda)
         !
       enddo
       !
     else
       !
       allocate(rho_(0:npoints_),stat=alloc)
       if (alloc/=0) stop 'rho_ - out of memory'
       !
       forall(i_ = 0:npoints_) rho_(i_)  =  rho_b(1)+real(i_,kind=ark)*rhostep_
       !
       do i = 0,npoints 
          !
          rho =  rho_b(1)+real(i,kind=ark)*rhostep
          !
          i_ = int( real(i,ark)*step_scale )
          !
          i1 = max(0,i_-Nr) ; i2 = min(npoints_,i_+Nr)
          !
          call polintark(rho_(i1:i2),poten_(i1:i2),rho,fval,df_t)
          poten(i) = fval
          !
          call polintark(rho_(i1:i2),mu_rr_(i1:i2),rho,fval,df_t)
          mu_rr(i) = fval
          !
          do lambda = 0,maxorder
            !
            call polintark(rho_(i1:i2),xton_(i1:i2,lambda),rho,fval,df_t)
            xton(i,lambda) = fval
            !
          enddo
          !
       enddo
       !
       deallocate(rho_)
       !
     endif
     !
     ! Do some reporting
     !
     if (verbose>=3) then 
         write (out,"('vmax        = ',i8)") vmax
         write (out,"('maxorder    = ',i8)") maxorder
         write (out,"('icoord      = ',i4)") icoord
         write (out,"('rho_b (x)   = ',2f12.4)") rho_b(1:2) !*180.0_rk/pi
         write (out,"('rhostep (x) = ',2f12.4)") rhostep  !*180.0_rk/pi
     endif 
     !
     call diff_2d_4points_ark(npoints,rho_b,mu_rr,periodic,0_ik,dfdr,d2fdr2)
     !
     f(:) = (d2fdr2(:)/mu_rr(:)-dfdr(:)**2/mu_rr(:)**2*0.5_ark)*0.5_ark
     !
     ! Print out
     !
     if (verbose>=3) then 
        write(out,"('grid values (i,rho,rho_kinet,rho_poten,poten, mu_rr, f): ')") 
        do i_=0,npoints_,2
          i = int( real(i_,ark)/step_scale )
          rho = rho_b(1)+real(i_,kind=ark)*rhostep_
          write(out,"(i8,3f14.6,3g14.6)") i_,rho,rho_kinet(i_),rho_poten(i_),poten(i),mu_rr(i),f(i)
        enddo 
     endif 
     !
     inquire(iolength=rec_len) phil(:),dphil(:)
     !
     write(unitfname,"('Numerov basis set # ',i6)") icoord
     call IOStart(trim(unitfname),io_slot)
     !
     open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
     !
     ! solve 1D Schroedinger equation with Numerov algorithm
     !
     iparity = 0
     !
     iswitch = max(1,int( rho_switch/rhostep_))
     !
     if (isingular>-1.and.(iswitch<=0.or.iswitch>=npoints_)) then
       write(out,"('me_numerov: illegal iswitch or npoints',2i8)") iswitch,npoints_
       stop 'iswitch inconsistent with npoints' 
     endif
     !
     if (iperiod/=0) vmax = vmax/2
     !
     if ( trim(molec%IO_primitive)/='READ') then
       call numerov(npoints,npoints_,step_scale,poten,mu_rr,f,enerslot)
     endif
     !
     if (iperiod/=0) then
       allocate(enerslot_(0:maxslots),stat=alloc)
       if (alloc/=0) then 
         write (out,"('phi - out of memory')")
         stop 'phi - out of memory'
       endif 
       enerslot_ =0
       ! 
       iparity = 1
       !
       if ( trim(molec%IO_primitive)/='READ') then
          call numerov(npoints,npoints_,step_scale,poten,mu_rr,f,enerslot_)
       endif
       !
       do vl = vmax,0,-1
         !
         enerslot(2*vl  ) = enerslot(vl)
         if (2*vl+1<=vmax_) enerslot(2*vl+1) = enerslot_(vl)
         !
       enddo
       !
       vmax = vmax_
       !
       deallocate(enerslot_)
       !
     endif
     !
     ! renormalising-re-orthogonalasing the basis set 
     !
     allocate(psi(0:npoints_,0:vmax),dpsi(0:npoints_,0:vmax),stat=alloc)
     if (alloc/=0) then 
       write (out,"('psi - out of memory')")
       stop 'psi - out of memory'
     endif 
     !
     do vl = 0,vmax
       !
       read (io_slot,rec=vl+1) (psi(i,vl),i=0,npoints_),(dpsi(i,vl),i=0,npoints_)
       !
     enddo
     !
     !omp parallel do private(vl,cross_prod,factor,vr) shared(psi,dpsi) schedule(dynamic)
     do vl =  0,vmax
       !
       phivphi(:) = psi(:,vl)*psi(:,vl)
       cross_prod = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
       !
       factor = 1.0_ark/sqrt(cross_prod)
       !
       psi(:,vl) = psi(:,vl)*factor
       dpsi(:,vl) = dpsi(:,vl)*factor
       !
       do vr = 0,vl-1
         !
         phivphi(:) = psi(:,vl)*psi(:,vr)
         cross_prod = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
         !
         psi(:,vl) = psi(:,vl)-cross_prod*psi(:,vr)
         dpsi(:,vl) = dpsi(:,vl)-cross_prod*dpsi(:,vr)
         !
         phivphi(:) = psi(:,vl)*psi(:,vl)
         cross_prod = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
         !
         factor = 1.0_ark/sqrt(cross_prod)
         !
         psi(:,vl) = psi(:,vl)*factor
         dpsi(:,vl) = dpsi(:,vl)*factor
         ! 
       enddo
       !
       write (io_slot,rec=vl+1) (psi(i,vl),i=0,npoints_),(dpsi(i,vl),i=0,npoints_)
       !
     enddo
     !omp end parallel do
     !
     deallocate(psi,dpsi)
     !
     ! Matrix elements 
     !
     sigma = 0.0_ark 
     rms   = 0.0_ark 
     characvalue = max(maxval(enerslot(0:vmax)),characvalue)
     !
     do vl = 0,vmax
        !
        read (io_slot,rec=vl+1) (phil(i),i=0,npoints_),(dphil(i),i=0,npoints_)
        !
        do vr = vl,vmax
            !
            if (iperiod/=0.and.mod(abs(vl-vr),2)==1) cycle
            !
            if (vl==vr) then
                phir =  phil
               dphir = dphil
            else
               read (io_slot,rec=vr+1) (phir(i),i=0,npoints_),(dphir(i),i=0,npoints_)
            endif
            !
            ! Here we prepare integrals of the potential 
            ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
            ! obtained above by the Numerov
            !
            phivphi(:) = phil(:)*poten_(:)*phir(:)
            !
            h_t = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
            !
            ! momenta-quadratic part 
            !
            phivphi(:) =-dphil(:)*mu_rr_(:)*dphir(:)
            !
            psipsi_t = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
            !
            ! Add the diagonal kinetic part to the tested mat. elem-s
            !
            h_t = h_t - 0.5_ark*psipsi_t
            !
            psipsi_t = 0 
            !
            do lambda = 0,maxorder
               !
               ! momenta-free part in potential part
               !
               if (lambda==0) then 
                  phivphi(:) = phil(:)*phir(:)
               else
                  phivphi(:) = phil(:)*rho_poten(:)**lambda*phir(:)
               endif
               !
               g_numerov(0,lambda,vl,vr) = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
               !
               ! external field expansion
               !
               if (lambda==0) then 
                  phivphi(:) = phil(:)*phir(:)
               else
                  phivphi(:) = phil(:)*rho_extF(:)**lambda*phir(:)
               endif
               !
               g_numerov(3,lambda,vl,vr) = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
               if (vl/=vr) g_numerov(3,lambda,vr,vl) = g_numerov(3,lambda,vl,vr)
               !
               ! momenta-free in kinetic part 
               !
               if (lambda==0) then 
                  phivphi(:) = phil(:)*phir(:)
               else
                  phivphi(:) = phil(:)*rho_kinet(:)**lambda*phir(:)
               endif
               !
               phivphi(:) = phil(:)*xton(:,lambda)*phir(:)
               !
               g_numerov(-1,lambda,vl,vr) = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
               !
               ! We also control the orthogonality of the basis set 
               !
               if (lambda==0) psipsi_t = g_numerov(0,lambda,vl,vr)
               !
               if (vl/=vr) g_numerov(-1:0,lambda,vr,vl) = g_numerov(-1:0,lambda,vl,vr)
               !
               ! momenta-quadratic part 
               !
               if (lambda==0) then 
                  phivphi(:) =-dphil(:)*dphir(:)
               else
                  phivphi(:) =-dphil(:)*rho_kinet(:)**lambda*dphir(:)
               endif
               !
               phivphi(:) =-dphil(:)*xton(:,lambda)*dphir(:)
               !
               g_numerov(2,lambda,vl,vr) = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
               !
               if (vl/=vr) g_numerov(2,lambda,vr,vl) = g_numerov(2,lambda,vl,vr)
               !
               ! momenta-linear part:
               ! < vl | d/dx g(x) | vr > = - < vr | g(x) d/dx | vl >
               !
               !
               if (lambda==0) then 
                  phivphi(:) = phil(:)*dphir(:)
               else
                  phivphi(:) = phil(:)*rho_kinet(:)**lambda*dphir(:)
               endif
               !
               phivphi(:) = phil(:)*xton(:,lambda)*dphir(:)
               !
               g_numerov(1,lambda,vl,vr) = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
               !
               if (vl/=vr) then
                  !
                  if (lambda==0) then 
                     phivphi(:) = dphil(:)*phir(:)
                  else
                     phivphi(:) = dphil(:)*rho_kinet(:)**lambda*phir(:)
                  endif
                  !
                  phivphi(:) = dphil(:)*xton(:,lambda)*phir(:)
                  !
                  g_numerov(1,lambda,vr,vl) = integral_rect_ark(npoints_,rho_b(2)-rho_b(1),phivphi)
                  !
               endif 
               !
               if(trim(extF%ftype)=='XY2_G-COR-ELEC'.or.trim(extF%ftype)=='XY2_G-TENS-RANK3') then
                  !
                  g_numerov(3,lambda,vl,vr) = g_numerov(1,lambda,vl,vr)
                  g_numerov(3,lambda,vr,vl) = g_numerov(1,lambda,vr,vl)
                  !
               endif               
               !
               !
               if (verbose>=7) then 
                   write(out,"('g_numerov(0,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(0,lambda,vl,vr)
                   write(out,"('g_numerov(1,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(1,lambda,vl,vr)
                   write(out,"('g_numerov(2,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(2,lambda,vl,vr)
                   write(out,"('g_numerov(3,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(3,lambda,vl,vr)
                   if (vl/=vr) then 
                     write(out,"('g_numerov(0,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(0,lambda,vr,vl)
                     write(out,"('g_numerov(1,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(1,lambda,vr,vl)
                     write(out,"('g_numerov(2,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(2,lambda,vr,vl)
                     write(out,"('g_numerov(3,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(3,lambda,vr,vl)
                   endif 
               endif 
               !
            enddo 
            !
            ! Count the error, as a maximal deviation sigma =  | <i|H|j>-E delta_ij |
            !
            sigma_t =  abs(h_t)
            if (vl==vr) then 
              enerslot(vl) = h_t
              sigma_t =  abs(h_t-enerslot(vl))
            endif
            !
            sigma = max(sigma,sigma_t)
            rms = rms + sigma_t**2
            !
            ! Now we test the h_t = <vl|h|vr> matrix elements and check if Numerov cracked
            ! the Schroedinger all right
            if (vl/=vr.and.abs(h_t)>sqrt(small_)*abs(characvalue)*1e4) then 
               write(out,"('ME_numerov: Integration differs from Numerovs solution for <',i4,'|H|',i4,'> = ',f20.10)") vl,vr,h_t
               write(out,"('            Try increasing the integration range.')")
               if ( trim(molec%IO_primitive)/='READ') then
                 stop 'ME_numerov: bad Numerov solution'
               endif
            endif 
            !
            if (vl==vr.and.abs(h_t-enerslot(vl))>sqrt(small_)*abs(characvalue)*1e4) then 
               write(out,"('ME_numerov: wrong <',i4,'|H|',i4,'> (',f16.6,') =/= energy (',f16.6,')')") vl,vr,h_t,enerslot(vl)
               write(out,"('            Try increasing the integration range.')")
               stop 'ME_numerov: bad Numerov solution'
            endif 
            !
            ! Reporting the quality of the matrix elemenst 
            !
            if (verbose>=3) then 
              if (vl/=vr) then 
               write(out,"('<',i4,'|H|',i4,'> = ',e16.2,'<-',8x,'0.0',5x,'; <',i4,'|',i4,'> = ',e16.2,'<-',8x,'0.0')") & 
                                vl,vr,h_t,vl,vr,psipsi_t
              else
               write(out,"('<',i4,'|H|',i4,'> = ',f16.6,'<-',f16.6,'; <',i4,'|',i4,'> = ',f16.6)")& 
                              vl,vr,h_t,enerslot(vl),vl,vr,psipsi_t
              endif 
            endif 
            !
        enddo
     enddo
     !
     energy(0:vmax) = enerslot(0:vmax)-enerslot(0)
     !
     rms = sqrt(rms/real((vmax+1)*(vmax+2)/2,kind=ark))
     !
     if (verbose>=1) then 
        write(out,"('Maximal deviation sigma =  | <i|H|j>-E delta_ij | is ',f18.8)") sigma
        write(out,"('rms deviation is ',f18.8)") sqrt(rms)
     endif 
     !
     deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,rho_poten,rho_extF,enerslot,f,poten,mu_rr,d2fdr2,dfdr,xton)
     !
     !
  end subroutine ME_numerov

!
!********************************************************************!
!                                                                    !
! numerov calculates bending functions and stores them               !
!                                                                    !
!                                                                    !
!********************************************************************!
  subroutine numerov(npoints,npoints_,step_scale,poten,mu_rr,f,enerslot) 
   !
   integer(ik),intent(in) :: npoints,npoints_
   real(ark),intent(in) :: step_scale
   real(ark),intent(in) :: poten(0:npoints),mu_rr(0:npoints),f(0:npoints)
   !
   real(ark),intent(out) :: enerslot(0:maxslots)

   integer(ik) :: v,i_
   integer(ik) :: ierr,numnod,ntries,vrec,k0,ipoint,mtries
   !
   real(ark) :: eold,tsum
   !
   real(ark) :: pcout
   !
   real(ark),allocatable :: pot_eff(:),i0(:),phi_f(:),phi_t(:),psi_t(:),phi_f_(:),phi_d_(:)
   !
   real(ark) :: eguess,enerlow,enerupp,enerdelta,enershift,enermid
   !
   !real(rk) :: simpsonintegral
   !
   logical :: notfound
   !
   real(ark) :: potmin,oldphi,newphi,v_t(-2:2),der(4),erest(4),v_tt(-4:4)
   !
   real(ark) :: r_t(1:4),r_tt,f_t(1:4),f_tt,df_t1,df_t2
   !
   integer(ik) :: alloc,i,j,jp2,iter,nm2,ic,imin_l,imin_r
   integer(ik) :: icslots(0:maxslots),iright,ileft,ireflect


   allocate(pot_eff(0:npoints),i0(0:npoints),phi_f(0:npoints),phi_t(0:npoints),psi_t(-30:npoints+30),phi_f_(0:npoints),&
            phi_d_(0:npoints),stat=alloc)
   if (alloc/=0) then 
      write (6,"('numerov: allocation is faild - out of memory')")
      stop 'numerov: allocation is faild - out of memory'
   endif 


   !fcoef = planck*avogno*1.0d+16/( 4.0d+00*pi*pi*vellgt )
   !
   ! Check for the potential "1/0" probelm in mu_rr
   !
   do i=0,npoints
     if ( abs(mu_rr(i))<small_*100.0_rk ) then 
        write(out,"('numerov: mu_rr is zero at i = ',i8)") i
        stop 'numerov: mu_rr has a zero element'
     endif
   enddo 
   !
   ! before commencing the numerov-cooley numerical integration we
   ! store the following function in the f1 array:
   !
   !                       0                 
   !     w  = f  (p ) + 2 i  v  (p )   
   !      i    1   i          0   i
   !
   !      in the renner-teller version of morbid, this function is
   !      extended with a term originating in the non-zero
   !      electronic angular momentum.
   !
   !                  -1       0                 -1 
   ! v  is given in cm  , and i   is given in  cm
   !  0                         
   !
   ! calculate effective potential w_i
   !
   pot_eff(:) = f(:) + 2.0_ark*( poten(:) )/mu_rr(:)
   !
   ! Effective inertia moment
   !
   i0(:) = 2.0_ark/mu_rr(:)
   !
   ! find minimum : position 
   !
   potmin = safe_max
   !
   do i=0,npoints
      !
      if (pot_eff(i)<potmin) then 
         imin = i
         potmin = pot_eff(i)
      endif
      !
   enddo
   ! and value 
   !
   potmin = pot_eff(imin)
   !
   potmin = max(0.0_ark,potmin)
   !
   !if (imin==0) imin = npoints/2
   !
   !imin = minloc(pot_eff(:),dim=1)
   !
   if (imin<0.or.imin>npoints) then 
       write(out,"('numerov: pot_eff has no minimum',i8)") 
       stop 'numerov: pot_eff has no minimum'
   endif 
   !
   imin_l = minloc(pot_eff(0:imin-1),dim=1)
   imin_r = minloc(pot_eff(imin+1:npoints),dim=1)+imin
   !
   !
   ileft  = npoints
   iright = 0
   !
   ! Print out
   !
   if (verbose>=5) then 
      write(out,"('grid values (poten, mu_rr, poten_eff): ')") 
      do i=0,npoints,2
        write(out,"(i8,3f18.8)") i,poten(i),mu_rr(i),pot_eff(i)
      enddo 
      write(out,"('potmin(eff) =  ',f18.8,' at i = ',i8)") potmin,imin
   endif 
   !
   ! assign initial values 
   !
   enerslot(0:maxslots)=poten(imin)-sqrt(safe_max)
   icslots(0:maxslots) = npoints
   !
   ! Parameters for the Simson rule integration 
   !
   ! facodd=2.0_rk*rhostep/3.0_rk
   ! faceve=2.0_rk*facodd
   !
   ! Here we start solving Schroedinger equation for the v-th eigenvalue
   !
   enerlow=potmin/2.0_ark*mu_rr(imin)
   enerupp=enermax
   !
   do v=0,vmax+1
     !
     if (v/=0) then 
       enerlow = enerslot(v-1)
     endif
     !
     ! Determine the highest undefined energy slot 
     i = v
     do while (enerslot(i)<poten(imin).and.i<maxslots)
       i = i + 1
     enddo 
     !
     if (enerslot(i)>=poten(imin)) then
        enerupp=enerslot(i)
     else 
        enerupp=enermax
     endif 
     !
     if (enerslot(v)>poten(imin)) then
        !
        enerlow=max(enerslot(v)-enerdev,poten(imin))
        enerupp=min(enerslot(v)+enerdev,enermax)
        ic  = icslots(v)
        !
     endif
     !
     enerdelta=  enerstep
     enershift = enerstep
     ntries=0
     mtries = 0

     !
     ! start itererative procedure for searching the current solution 
     !
     iter = 0
     !
     notfound = .true.
     !
     search_loop : do while(notfound.and.iter<itermax) 
       !
       iter = iter + 1
       !
       enermid=(enerlow+enerupp)/2.0
       !
       !enerdelta = enerdelta*0.5_ark
       !
       eguess=enermid
       !
       eold=eguess
       !
       ! Numerov procedure 
       !
       ic = icslots(v)
       !
       if (ic/=npoints.and.enerslot(v)>=poten(imin).and.iter==1) eguess = enerslot(v)
       !
       if (verbose>=5) write(out,"('eguess = ',e14.7)") eguess
       !
       ! special case of a singular point i = 0
       !
       if (isingular==0) then 
          !
          forall(i = 0:npoints) phi_f(i)  = sqrt(rhostep*real(i,ark))!*sqrt(mu_rr(i))
          !
       endif
       !
       call numcoo ( v, pot_eff, i0, eguess, enerlow, ic, phi_f, pcout, ierr)
       ! 
       if (ierr>1) then 
           !
           write(out,"('numerov: no solution found in numcoo, ierr = ',i8)") ierr
           !
           do i = 0,maxslots
             if (enerslot(i)>poten(imin)) write (out,"('        v,ener = ',i8,f20.10)") i,enerslot(i)
           enddo 
           !
           stop 'numerov: no solution found in numcoo'
           !
       endif 
       ! 
       !
       if (ierr == 0) then  

        !if (trim(molec%coords_transform)=='R-EXPRHO'.and.imode==molec%Nmodes) then 
        !  numnod=0
        !  nm2=npoints-2
        !  oldphi=phi_f(0)*exp((rho_b(1)+rhostep*real(0,rk)))
        !  newphi=phi_f(2)*exp((rho_b(1)+rhostep*real(2,rk)))
        !  do i=2,nm2
        !     
        !     if (oldphi*phi_f(i)*exp((rho_b(1)+rhostep*real(i,rk)))<0.0_ark .or.  & 
        !          ( abs(phi_f(i)*exp((rho_b(1)+rhostep*real(i,rk))))==0.0_ark .and. newphi*oldphi<0.0_ark) ) then 
        !          if (exp( (rho_b(1)+rhostep*real(i,rk)) )>1e-4) numnod=numnod+1
        !     endif 
        !     ! 
        !     newphi=phi_f(i+2)*exp((rho_b(1)+rhostep*real(i+2,rk)))
        !     oldphi=phi_f(i)*exp((rho_b(1)+rhostep*real(i,rk)))
        !  enddo
        !else

          numnod=0
          nm2=npoints-2
          oldphi=phi_f(0)
          newphi=phi_f(2)
          do i=2,nm2
             if ( oldphi*phi_f(i)<0.0 .or.  & 
                  ( phi_f(i)==0.0_ark  .and. newphi*oldphi<0.0_ark) ) numnod=numnod+1
             ! 
             newphi=phi_f(i+2)
             oldphi=phi_f(i)
          enddo
          !
       endif 
       !
       if (ierr==0.and.periodic) then 
           !
           phi_t(:) = phi_f(:)/sqrt(mu_rr(:))
           !
           tsum = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phi_t(:)**2)
           !
           phi_t(:)=phi_t(:)/sqrt(tsum)
           !
           ireflect = 0
           !
           if (iparity/=0) ireflect = -1
            !
           call diff_2d_4points_ark(Npoints,rho_b,phi_t,periodic,ireflect,psi_t(0:Npoints))
           !
           do i=2,5
            r_t(i-1) = rho_b(1)+real(i,kind=ark)*rhostep
           enddo
           r_tt = rho_b(1)
           f_t(1:4) = psi_t(2:5)
           call polintark(r_t(1:4),f_t,r_tt,df_t1,f_tt)
           !
           do i=1,4
            r_t(i) = rho_b(1)+real(npoints-6+i,kind=ark)*rhostep
           enddo
           r_tt = rho_b(2)
           f_t(1:4) = psi_t(npoints-5:npoints-2)
           call polintark(r_t(1:4),f_t,r_tt,df_t2,f_tt)
           !
           if (abs(phi_t(0)-(-1.0_ark)**iparity*phi_t(npoints  ))>sqrt(thrsh_int).or.&
               abs(df_t1-(-1.0_ark)**iparity*df_t2)>1000.*sqrt(thrsh_int)) then 
              !
              if (verbose>=5) write(out,"(/'f(0),f(1),.,f(N-1),f(N),f(0)-tau*f(N),energy: ',5g18.8,', energy= ',f12.4,',n=',i9)")&
                              phi_t(0),phi_t(npoints),phi_t(1),phi_t(npoints-1),& 
                                             phi_t(0)-(-1.0_rk)**iparity*phi_t(npoints ),eguess,numnod
              if (verbose>=5) write(out,"(/'df(0)-df(N): ',3g18.8)") df_t1,df_t2,df_t1-(-1.0_ark)**iparity*df_t2
              !
              if ( numnod+1<maxslots.and.numnod>v) then 
                  enerslot(numnod)=eguess
                  icslots(numnod) = ic 
              endif 
              !
              numnod = maxslots+1
              !
           endif 
       endif
       !
       !
       if (ierr==0.and.iperiod<0) then 
           !
           phi_t(:) = phi_f(:)/sqrt(mu_rr(:))
           !
           tsum = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phi_t(:)**2)
           !
           phi_t(:)=phi_t(:)/sqrt(tsum)
           !
           ireflect = 1
           if (iparity/=0) ireflect = -1
           !
           call diff_2d_4points_ark(Npoints,rho_b,phi_t,periodic,ireflect,psi_t(0:Npoints))
           !
           if ( (iparity==0.and.abs(psi_t(npoints  ))>sqrt(thrsh_int)).or.&
                (iparity==1.and.abs(phi_t(npoints  ))>sqrt(thrsh_int))) then 
              !
              if (verbose>=5) write(out,"(/'phi(N) and psi(N) = ',2g18.8,', energy = ',f12.4,', n = ',9i7)") &
                              phi_t(npoints),psi_t(npoints ),eguess,numnod
              !
              if ( numnod+1<maxslots.and.numnod>v) then 
                  enerslot(numnod)=eguess
                  icslots(numnod) = ic 
              endif 
              !
              numnod = maxslots+1
              !
           endif 
           !
       endif
       ! 

       !if (iperiod<0) then 
       !   ! 
       !   if (abs(phi_f(npoints))<small_.and.mod(v,2)/=0) then 
       !     !
       !     numnod = numnod*2+1
       !    !
       !   else
       !     !
       !     numnod = numnod*2
       !     !
       !   endif
       !   !
       !endif 
       !
       if ( ierr==0.and.numnod+1<maxslots.and.numnod>=v ) then 
          enerslot(numnod)=eguess
          icslots(numnod) = ic 
          if (verbose>=6) then 
             write (out,"('v,numnod,ener = ',2i8,f20.10)") v,numnod,enerslot(numnod)
             !
             do i=0,npoints 
                !
                write(out,"(i8,2f18.8)") i,phi_f(i),exp((rho_b(1)+rhostep*real(i,rk)))
                !
             enddo
             !
          endif 
       endif 
       !
       if (verbose>=5) write(out,"('v = ',i5,'; numnod = ',i8,', efound = ',f16.7,', ierr = ',i9)") v,numnod,eguess,ierr
       !
       if ( ierr == 0.and.v==numnod ) notfound = .false.
       !
       ! if it cannot find a solution by Numerov procedure and ierr = 1, we can still try to
       ! devide the searching interval and repeate
       ! the Numerov procedure with a closer initial guess
       ! We do the same thing if the found solution is not what we need (/=v)
       !
       if (abs(enerupp-enerlow) < enersmall.and.(ierr == 1 .or. numnod  /= v)) then 
          !
          ntries=0
          !
          !enerdelta = enerdelta*0.5_ark
          !
          ! Are we still at the begining?
          !
          if (v.eq.0) then
             !
             enerlow=poten(imin)
             !
          !elseif (iperiod<0) then
          !   !
          !   enerlow=poten(imin)
          !   !
          else
            !
             enerlow=enerslot(v-1)-0.001
             !
             !enerlow=poten(imin)
          endif
          !
          ! Determine the highest undefined energy slot 
          i = v+1
          do while (enerslot(i)<poten(imin).and.i<maxslots)
            i = i + 1
          enddo 
          !
          if (enerslot(i)>=poten(imin)) then
             enerupp=enerslot(i)
          else 
             enerupp=enermax
          endif 
          !
          if (enermid>=enerupp) then 
             enerdelta = enerdelta*0.5_ark
             enershift = 0.0_ark
          endif 
          !
          enerlow=enerlow+enershift
          enershift=enershift+enerdelta
          enerupp=enerlow
          !
          cycle search_loop

       endif

       !
       if (ierr <=1 .and. numnod  /= v) then 
          !
          ! We need to do something if the seacrh interval is collapsed
          ! this must have happend at ierr=1, and no solution found 
          !
          if (mtries==0) then
             !
             mtries = 1
             ! 
             ! Are we still at the begining?
             !
             if (v==0) then
                enerlow=poten(imin)
             else
                enerlow=enerslot(v-1) 
                !enerlow=poten(imin)
             endif
             !
             ! Determine the highest undefined energy slot 
             i = v+1
             do while (enerslot(i)<poten(imin).and.i<maxslots)
               i = i + 1
             enddo 
             !
             if (enerslot(i)>=poten(imin)) then
                enerupp=enerslot(i)
             else 
                enerupp=enermax
             endif 
             !
          endif 
          !
          if (numnod<v.or.ierr==1) then
             enerlow=enermid
          else
             enerupp=enermid
          endif
          !
       endif
       !
     enddo search_loop
     !
     if( notfound ) then
       !
        write(out,"('numerov: no solution found after ',i8,' iterations')") itermax
        write (out,"('        v = ',i8,', eguess= ',e20.10)") v,eguess
        !
        do i = 0,maxslots
          if (enerslot(i)>poten(imin)) write (out,"('        v,ener = ',i8,f20.10)") i,enerslot(i)
        enddo 
        !
        stop 'numerov: no solution found'
        !
     endif 
     !
     if (ierr==0.and.periodic) then 
        if (verbose>=5) write(out,"(/'phi_t(0)-phi_t(npoints  )    : ',3g18.8)") phi_t(0),phi_t(npoints  ),phi_t(0)-phi_t(npoints  )
        if (verbose>=5) write(out,"(/'dphi_t(0)-dphi_t(npoints)    : ',5g18.8)") phi_t(1),phi_t(npoints-1),df_t1,df_t2,df_t1-df_t2
     endif 
     !
     ! multiply the wavefunction with sqrt(irr) (eq. (6.4) of jensen)
     ! 
     phi_f(:) = phi_f(:)/sqrt(mu_rr(:))
     !
     phi_t(:) = phi_f(:)*phi_f(:)
     !
     !if (trim(molec%coords_transform)=='R-RHO'.and.imode==molec%Nmodes) then 
     !  !
     !  do i = 0,Npoints
     !    !
     !   phi_t(i) = phi_t(i)*sqrt(rho_b(1)+rhostep*real(i,rk) )**1
     !   !
     !  enddo
     !  !
     !endif 
     !
     ! normalize the wavefunction using the Simpson rule integration 
     !
     !tsum = sum(phi_f(:)*phi_f(:):0:2)*facodd + sum(phi_f(:)*phi_f(:):0:2)*faceve
     !
     !do i=0,npoints-1,2
     !   tsum=tsum+facodd*phi_f(i)*phi_f(i)+faceve*phi_f(i+1)*phi_f(i+1)
     !enddo             ! --- i 
     !
     !   numerical intagration with simpson's rule #2
     !
     tsum = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phi_t)
     !
     phi_f(:)=phi_f(:)/sqrt(tsum)
     !
     !if (iperiod<0) phi_f(:)=phi_f(:)/sqrt(2.0_ark)
     !
     !
     if (verbose>=5) then 
        write (out,"('v,ener = ',i8,f20.10)") v,enerslot(numnod)
     endif 
     !
     if (verbose>=6) then 
        write(out,"(f18.8)") phi_f
     endif 
     !
     ! direvatives of the wavefunctions method I
     !
     if(deriv_method == '5 points') then 
        !
        do i=0,npoints      ! --- np
           !
           do k0 = -2, 2
              !
              if (periodic) then 
                 ipoint = min(max(i+k0,0),npoints)
                 v_t(k0) = phi_f(ipoint)
              else
                 if (i+k0<0) then 
                    v_t(k0) = phi_f(-(i+k0))
                 elseif(i+k0>npoints) then 
                    v_t(k0) = phi_f(npoints)
                 else
                    v_t(k0) = phi_f(i+k0)
                 endif
              endif
              !
           enddo
           !
           ! 5-points expression
           !
           phi_t(i) = (-v_t( 2)/12.0_ark+2.0_ark/3.0_ark*v_t(1) & 
                       +v_t(-2)/12.0_ark-2.0_ark/3.0_ark*v_t(-1) )/rhostep
           !
           !phi_t(i) = (v_t(1)-v_t(-1) )/(rhostep*2.0_ark)
           !
        enddo 
        !
        ! direvatives of the wavefunctions method II
        !
     elseif(deriv_method == 'ML_diffs') then 
        !
        if (iperiod<0) then 
          !
          ireflect = 1
          if (iparity/=0) ireflect = -1
          !
          call diff_2d_4points_ark(Npoints,rho_b,phi_f,periodic,ireflect,phi_t)
          !
        else 
          !
          ireflect = 0
          if (iparity/=0) ireflect = -1
          !
          call diff_2d_4points_ark(Npoints,rho_b,phi_f,periodic,ireflect,phi_t)
          !
        endif
        !
        ! direvatives of the wavefunctions method III
        !
     else
         write(out,"('numerov: bad deriv_method ',a)") trim(deriv_method)
               stop 'numerov: bad deriv_method'

     endif
     !
     !
     if (verbose>=5) then 
        !
        do i=0,npoints 
           !
           write(out,"(i8,2f18.8,' || ',1x,i8)") i,phi_f(i),phi_t(i),v
           !
        enddo
        !
     endif 
     !
     ! Find the outmost points where the function is not zero
     !
     i = 0 
     !
     do while (abs(phi_f(i))<100.0*sqrt(small_).and.i<npoints)
      i = i + 1
     enddo
     !
     ileft = min(ileft,i)
     !
     i = npoints 
     !
     do while (abs(phi_f(i))<100.0*sqrt(small_).and.i>0)
      i = i - 1
     enddo
     !
     iright = max(iright,i)
     !
     ! dump the eigenfunction
     !
     if (npoints_==npoints) then 
          phi_f_ = phi_f
          phi_d_ = phi_t
     else
       do i_ = 0,npoints_
          !
          i = nint( real(i_,ark)/step_scale )
          phi_f_(i_) = phi_f(i)
          phi_d_(i_) = phi_t(i)
          !
       enddo
     endif
     !
     vrec=v+1
     !
     if (iperiod/=0) then 
       !
       vrec = 2*v+1 
       !
       if (iparity/=0) vrec = 2*v+2
       !
     endif
     !
     write (io_slot,rec=vrec) (phi_f_(i_),i_=0,npoints_),(phi_d_(i_),i_=0,npoints_)
     !
   enddo
   !
   write(out,"(/' Outmost points of the nonvanishing wave-functions are: [',i6,'...',i6,'], i.e. [',f12.6,'...',f12.6,'] ')") &
                 ileft,iright,rho_b(1)+rhostep*real(ileft,rk),rho_b(1)+rhostep*real(iright,rk)
   !
   ! report the found Numerov energies 
   !
   write (out,"(/' Numerov-energies are:')") 
   !
   enerlow = minval(enerslot(0:vmax),dim=1) !mask=(enerslot>=poten(imin)
   !
   do v=0,vmax   
     write (out,"(i8,f18.8)") v,enerslot(v)-enerlow
   enddo
   !
   write (out,"('Zero-point-energy is :',f18.6)") enerlow
   !
   deallocate(pot_eff,i0,phi_f,phi_t,phi_f_,phi_d_)
   !
  end subroutine numerov



!     subroutine solves eqn 6.5 , for the given value of kqua ,
!     by the numerov - cooley technique.
!
!                  (                       0                  )
!      2           (                    2 i                   )
!     d            (                2      pp                 )
!     ---phi (p) = ( f (p) + f (p) k  + ----- ( v   (p) - e ) ) phi (p)
!       2   b      (  1       2         _2       eff          )    b
!     dp           (                    h                     )
!
!     for this the following defintions (equations 6.14 to 6.17) are
!     needed :-
!
!     ph  = phi  (p )
!       i      b   i
!
!             0         _2
!     i  = 2 i   (p ) / h
!      i      pp   i
!
!                             2      0                   _2
!     u  = f  (p ) + f  (p ) k  + 2 i   (p ) v    (p ) / h
!      i    1   i     2   i          pp   i   eff   i
!
!                 2
!     y  = ( 1 - h  ( u  - i  e ) ) ph
!      i         --    i    i         i
!                12
!
!     starting from an initial guess for the energy , given in wavfun
!     an iterative scheme is followed. this is performed by using the
!     recursion relation (eqn 6.18) :-
!
!                           2
!     y    + y    - 2 y  = h  ( u  - i  e ) ph
!      i+1    i-1      i         i    i       i
!
!     h is here the rhostep length used to generate the f1 , f2 , i0's
!
!     to calculate values of y(i) and phi_f(i) both from rho = 0
!     and from rho = rhomax. these two parts are called the outward
!     and inward integrations respectively.
!     the outward integration is carried out until a maximum is reached
!     in the wavefunction. this crossing point (ic , yc , pc) is then
!     the stopping point used for the inward integration.
!     both sets of wavefunctions are scaled so that pc(out) = pc(in) = 1
!     an error for the energy can be calculated from eqn 6.24 :-
!
!                                                         npoints
!               out           in       2                    \--      2
!     d(e) = (-y    + 2 y  - y    ) / h + u  - i  e ) y  /   \  i  ph
!               c-1      c    c+1          c    c      c     /   i   i
!                                                           /--
!                                                           i=1
!
!     and this is added to the eguess at each iteration , when d(e)
!     is less than thrsh3 the iteration has converged.
!

  subroutine numcoo( v, pot_eff, i0, eguess, enerlow, iref, phi_f, pcout, ierr )

     integer(ik),intent(in) :: v
     integer(ik),intent(inout) :: iref
     real(ark),intent(in)  ::  i0(0:npoints),pot_eff(0:npoints)
     real(ark),intent(in)  ::  enerlow
     real(ark),intent(inout)  :: eguess
     real(ark),intent(out) :: phi_f(0:npoints)
     real(ark),intent(out) :: pcout
     integer(ik),intent(out) :: ierr

     real(ark) :: hh,dx,sumout,sumin,tsum,ycm1,pcin,ycp1,yc,phi_t,x(0:epoints)
     integer(ik) :: niter,ic,i,istart
     !
     !y(i)=(1.0_rk-hh*(pot_eff(i)-i0(i)*eguess)/12.0_rk)*phi_f(i)
     !
     hh=rhostep*rhostep
     !
     !
     niter=0
     dx = safe_max
     !
     pcout = 1.0_ark
     pcin  = 1.0_ark
     !
     if (mod(v,2)/=0)  pcout = -pcout
     !
     if ( .not.periodic.and.iperiod<0.and.iparity==0) pcout = 1.0_ark
     !
     if (periodic.and.iperiod>0) pcout = pcout*(-1.0_ark)**iparity
     !
     ierr = 0
     istart = 1
     !
     do while ( ierr==0 .and. abs(dx)>thrsh_int  )   ! --- niter 
        !
        if (.not.periodic.and.iperiod<0) then 
            !
            !if (mod(v,2)==0.and.pcin*pcout>small_) then 
            if (iparity==0) then 
              !
              phi_t = sqrt(small_)
              !
              phi_f(npoints  )  = phi_t
              phi_f(npoints-1)  = phi_f(npoints)*(12.0_ark+5.0_ark*hh*(pot_eff(npoints  )-i0(npoints  )*eguess))/ &
                                                 (12.0_ark-        hh*(pot_eff(npoints-1)-i0(npoints-1)*eguess))
            else
              !            
              phi_f(npoints-1)  = (small_) ! small_
              phi_f(npoints  )  = 0
              !
            endif 
            !
            phi_f(0)  = 0
            phi_f(1)  = small_
            !
         elseif (isingular==0) then 
            !
            istart = iswitch+2
            phi_f(npoints-1)  = safe_min ! small_
            phi_f(npoints  )  = 0.0_ark
            !
         elseif (abs(rho_b(1))<=0.01_ark.and..not.periodic) then 
            !
            phi_f(1)  = sqrt(small_) ! small_
            phi_f(0)  = small_
            !
          !endif 
          !
          !phi_f(0)  = wave_init
          !phi_f(1)  = phi_f(0)*(12.0_rk+5.0_rk*hh*(pot_eff(0)-i0(0)*eguess))/ &
          !                     (12.0_rk-       hh*(pot_eff(1)-i0(1)*eguess))
          !
          !p1 := -P0*(12+5*h^2*UIE0)/(-12+h^2*UIE1)
          !
          !
          phi_f(npoints-1)  = sqrt(small_) ! small_
          phi_f(npoints  )  = small_
          !  !
          !if (imode==3) then ! .and.mod(v,2)==0.and.v/=100) then 
          !   !
          !   phi_t = wave_init
          !   !
          !   phi_f(0)  = phi_t
          !   phi_f(1)  = phi_f(0)*(12.0_ark+5.0_ark*hh*(pot_eff(0)-i0(0)*eguess))/ &
          !                        (12.0_ark-        hh*(pot_eff(1)-i0(1)*eguess))
          !   !
          !endif 
          !
        !elseif (naught_at(1)==0.and.naught_at(2)==1) then 
        !  !
        !  phi_f(npoints  )  = wave_init
        !  phi_f(npoints-1)  = phi_f(npoints)*(12.0_ark+5.0_ark*hh*(pot_eff(npoints)-i0(npoints)*eguess))/ &
        !                                     (12.0_ark-        hh*(pot_eff(npoints-1)-i0(npoints-1)*eguess))
        !  !
        !  phi_f(1)  = safe_min ! small_
        !  phi_f(0)  = 0.0_ark

            !
        elseif (periodic) then 
          !
          if (pcin*pcout>small_.or.(mod(v,2)/=0.and.iperiod>0.and.iparity==1)) then 
            !
            phi_t = wave_init
            !
            !
            phi_f(npoints  )  = phi_t
            phi_f(npoints-1)  = phi_f(npoints)*(12.0_ark+5.0_ark*hh*(pot_eff(npoints)-i0(npoints)*eguess))/ &
                                               (12.0_ark-        hh*(pot_eff(npoints-1)-i0(npoints-1)*eguess))
            !
            phi_f(0)  = phi_t

            if ((mod(v,2)/=0.and.iperiod>0.and.iparity==1)) phi_f(0) = -phi_f(0)

            phi_f(1)  = phi_f(0)*(12.0_ark+5.0_ark*hh*(pot_eff(0)-i0(0)*eguess))/ &
                                 (12.0_ark-        hh*(pot_eff(1)-i0(1)*eguess))



           !
          else
            !
            phi_f(npoints-1)  = small_ !sqrt(small_) ! small_
            phi_f(npoints  )  = 0
            !
            phi_f(1)  = -small_
            phi_f(0)  = 0
            !
          endif 
          !
          !endif 
          !
          !phi_f(npoints  )  = small_ !sqrt(small_) ! small_
          !phi_f(npoints-1)  = sqrt(small_)
          !
          !phi_f(1)  = sqrt(small_) ! small_
          !phi_f(0)  = small_
          !
          !endif 
          !
          !
        else
          !
          phi_f(0)  = 0.0_ark
          phi_f(1)  = small_
          phi_f(npoints-1)  = safe_min ! small_
          phi_f(npoints  )  = 0.0_ark
          !
          !if (pcin*pcout>small_) then 
          !  !
          !if (imode==3.and.mod(v,2)==0.and.v/=0) then 
          !   !
          !   phi_t = wave_init
          !   !
          !   phi_f(0)  = phi_t
          !   phi_f(1)  = phi_f(0)*(12.0_ark+5.0_ark*hh*(pot_eff(0)-i0(0)*eguess))/ &
          !                        (12.0_ark-        hh*(pot_eff(1)-i0(1)*eguess))
          !   !
          !endif 
          !endif 
          !
        endif 
        !
        ! Outer integration 
        !
        !
        call intout ( v, pot_eff, i0, eguess, phi_f, sumout, istart, iref,  ic, pcout,  ierr)
        !
        if (ierr/=0) return
        !
        ! Inner integration 
        !
        call intin ( pot_eff, i0, eguess, phi_f, sumin, ic, pcin)
        !
        tsum=sumin/(pcin*pcin)+sumout/(pcout*pcout)
        !
        if (tsum > safe_max ) then  ! thrsh_upper
           ierr=2 ! 3 
           write (out,"(' numerov:  no solution, tsum > thrsh_upper:',2g18.8 )") tsum,thrsh_upper
           return
        endif 
        !
        yc=y(ic)/pcout
        ycm1=y(ic-1)/pcout
        ycp1=y(ic+1)/pcin
        !
        dx=((-ycm1+yc+yc-ycp1)/hh+(pot_eff(ic)-i0(ic)*eguess))*yc/tsum
        !
        eguess=eguess+dx
        !
        niter=niter+1
        !
        if (niter==iterlimit) then 
           ierr=2
           write (out,"(' numerov:  iteration limit reached in numerical integration (niter=',i8)") iterlimit
           write (out,"('           energy difference =',d13.4)") dx
           write (out,"('           convergence threshold =',d13.4,/)") thrsh_int
           return
        endif
        !if (eguess<enerlow) then 
        !   ierr=1
        !endif 
        !
     enddo
     if (eguess<enerlow) then 
        ierr=1
     endif 
     !
     iref = ic
     !
     phi_f(0:ic)=phi_f(0:ic)/pcout
     phi_f(ic+1:npoints)=phi_f(ic+1:npoints)/pcin
     !
  contains 
     !
     function y(i) result (v)
     integer(ik),intent(in) :: i
     real(ark) :: v
        !
        v=(1.0_ark-hh*(pot_eff(i)-i0(i)*eguess)/12.0_ark)*phi_f(i)
        !
     end function y
     !
  end subroutine numcoo

!
!     subroutine performs the outward integration of numcoo ,
!     until the first maximum in the wave function is found.
!     the sum i0(i)*phi_f(i)**2 is saved for the outward integration
!     and will later br divided by pc(out)**2.

  subroutine intout ( v, pot_eff, i0, eguess, phi_f, sumout,istart, iref, ic, pcout, ierr)
      !
     integer(ik),intent(in ) :: v,istart,iref
     real(ark),intent(in ) :: eguess
     integer(ik),intent(inout) :: ierr
     integer(ik),intent(out) :: ic


     real(ark),intent(in)  ::  i0(0:npoints),pot_eff(0:npoints)
     real(ark),intent(inout) :: phi_f(0:npoints)
     real(ark),intent(out) :: sumout,pcout
     !
     integer(ik) :: i,iend,imin_ref
     logical :: notfound
     !
     real(ark) :: hh,const,tsum,redfac,yi,yim1,yip1,phi_t
     !
     hh=rhostep*rhostep
     !
     !istart=1
     !
     yim1=phi_f(istart-1)*(1.0_ark-hh*(pot_eff(istart-1)-i0(istart-1)*eguess)/12.0_ark)
     !
     yi=phi_f(istart)*(1.0_ark-hh*(pot_eff(istart)-i0(istart)*eguess)/12.0_ark)
     !
     const=hh*(pot_eff(istart)-i0(istart)*eguess)
     !
     yip1=const*phi_f(istart)+yi+yi-yim1
     !
     tsum=0.0_ark
     !
     notfound = .true.
     !
     i = istart-1
     iend = npoints-2 
     if (v==0) iend = imin
     !
     ! if minimum is at i= 0 choose iend at the moddle 
     !
     !
     imin_ref = imin
     !
     if (isingular==0) then 
       !
       !if (imin<=5) 
       !imin_ref = npoints/2
       !if (iend<=5) 
       !iend = npoints/2
       !
     endif
     !
     do while(notfound.and.i<=iend)
        !
        i = i + 1
        !
        const=hh*(pot_eff(i)-i0(i)*eguess)
        !
        phi_t = yi/(1.0_ark-const/12.0_ark)
        !
        phi_f(i)=phi_t
        !
        yip1=const*phi_t+yi+yi-yim1
        !
        if (i>=imin_ref) then 
          !
          if ( sign( 1.0_ark,yi-yim1 )/=sign(1.0_ark,yip1-yi).and.i.ne.1) then 
          !if ( sign( 1.0_rk,yim1 )/=sign(1.0_rk,yip1 ).and.i.ne.1) then 
              notfound = .false.
              cycle 
          endif 
          !
          if ( i==iref ) then 
          !if ( sign( 1.0_rk,yim1 )/=sign(1.0_rk,yip1 ).and.i.ne.1) then 
              notfound = .false.
              cycle 
          endif 
          !
        endif 
        !
        yim1=yi
        !
        yi=yip1
        !
        tsum=tsum+i0(i)*phi_f(i)*phi_f(i)
        !
        ! renormalizing phi_f
        !
        if (tsum > thrsh_upper) then 
            redfac=sqrt(thrsh_upper)
            !
            phi_f(0:i)=phi_f(0:i)/redfac
            yi=yi/redfac
            yim1=yim1/redfac
            tsum=tsum/thrsh_upper
        endif 
        !
     enddo
     !
     sumout=tsum
     !
     if (notfound) then 
        if (v/=0) ierr=1
        ic=imin_ref
        if (verbose>=7) then 
          write (out,"('intout: no turning point found in outward integration')")
        endif 
     else
        sumout=tsum
        ic=i-1
        pcout=phi_f(ic)
     endif


     pcout=phi_f(ic)


  end subroutine intout


!     subroutine performs the inward integration from rho=rhomax
!     stopping at the maximum in the wavefunction determined in
!     intout. the recursion relations given in numcoo are used
!     in conjunction with the two starting values :-
!
!     phi(rhomax) = 0
!
!     phi(rhomax-rhostep) = small_
!
!
  subroutine intin ( pot_eff, i0, eguess, phi_f , sumin , ic , pcin )
      
     real(ark),intent(in ) :: eguess
     integer(ik),intent(in) :: ic

     real(ark),intent(in)  ::  i0(0:npoints),pot_eff(0:npoints)
     real(ark),intent(inout) :: phi_f(0:npoints)
     real(ark),intent(out) :: sumin,pcin
     !
     real(ark) :: hh,const,tsum,yi,yim1,yip1,redfac,phi_t
     integer(ik) :: i,ist,nm1,iend,kend,km1
     !
     hh=rhostep**2
     !
     ist=ic+1
     iend=npoints-1
     nm1=iend
     !
     yip1=phi_f(nm1+1)*(1.0_ark-hh*(pot_eff(nm1+1)-i0(nm1+1)*eguess)/12.0_ark)
     yi  =phi_f(nm1  )*(1.0_ark-hh*(pot_eff(nm1  )-i0(nm1  )*eguess)/12.0_ark)
     !
     if (periodic.and.abs(phi_f(0))>sqrt(small_).and..false.) then 
       !
       kend = 2
       km1 = kend
       !
       yip1=phi_f(km1+1)*(1.0_ark-hh*(pot_eff(km1+1)-i0(km1+1)*eguess)/12.0_ark)
       yi  =phi_f(km1  )*(1.0_ark-hh*(pot_eff(km1  )-i0(km1  )*eguess)/12.0_ark)
       !
      ! do i=kend,0,-1
      !   !
      !   const=hh*(pot_eff(i)-i0(i)*eguess)
      !   !
      !   phi_t=yi/(1.0_ark-const/12.0_ark)
      !   !
      !   phi_f(i)=phi_t
      !   !
      !   yim1=const*phi_t+yi+yi-yip1
      !   !
      !   yip1=yi
      !   yi=yim1
      !   !
      ! enddo
       !
       !iend=npoints-2
       !nm1=iend
       !
       !phi_f(npoints-2) = 2.0_ark*phi_f(npoints)-phi_f(2)
       !
       !yip1=phi_f(nm1+1)*(1.0_ark-hh*(pot_eff(nm1+1)-i0(nm1+1)*eguess)/12.0_ark)
       !yi  =phi_f(nm1  )*(1.0_ark-hh*(pot_eff(nm1  )-i0(nm1  )*eguess)/12.0_ark)
       !
     endif
     !
     !yip1=0.0_rk
     !
     !yi=phi_f(nm1)*(1.0_rk-hh*(pot_eff(nm1)-i0(nm1)*eguess)/12.0_rk)
     !
     tsum=0.0_ark
     !
     do i=iend,ist,-1
        !i=ist+iend-ii
        const=hh*(pot_eff(i)-i0(i)*eguess)
        !
        phi_t=yi/(1.0_ark-const/12.0_ark)
        !
        phi_f(i)=phi_t
        !
        yim1=const*phi_t+yi+yi-yip1
        !
        yip1=yi
        yi=yim1
        !
        tsum=tsum+i0(i)*phi_f(i)*phi_f(i)
        !
        ! renormalizing phi_f
        !
        if (tsum > thrsh_upper) then 
            redfac=sqrt(thrsh_upper)
            !
            phi_f(i:iend)=phi_f(i:iend)/redfac
            yi=yi/redfac
            yip1=yip1/redfac
            tsum=tsum/thrsh_upper
        endif 
        !
     enddo
     !
     sumin=tsum
     !
     pcin=yi/(1.0_ark-hh*(pot_eff(ic)-i0(ic)*eguess)/12.0_ark)
     !
  end subroutine intin

!
! integration with Simpson rules 
!                                      
  function simpsonintegral_ark(npoints,xmax,f) result (si) 
    integer(ik),intent(in) :: npoints
    !
    real(ark),intent(in) :: xmax,f(0:npoints)
    !
    real(ark) :: si
    !
    integer(ik) :: i
    !
    real(ark) ::  feven,fodd,f0,fmax,h
      !
      h = xmax/real(Npoints,kind=ark)  !   integration step   
      feven=0         
      fodd =0
      f0   =f(0)
      fmax =f(Npoints)

     !
     !  sum of odd and even contributions 
     !
     do i = 1,npoints-2,2
        fodd   = fodd  + f(i  )
        feven  = feven + f(i+1)
     enddo
     !
     fodd   = fodd  + f(npoints-1)
     !
     si =  h/3.0_ark*( 4.0_ark*fodd + 2.0_ark*feven + f0 + fmax)

  end function  simpsonintegral_ark


!
! integration with Simpson rules 
!                                      
  function simpsonintegral(npoints,xmax,f) result (si) 
    integer(ik),intent(in) :: npoints
    !
    real(rk),intent(in) :: xmax,f(0:npoints)
    !
    real(rk) :: si
    !
    integer(ik) :: i
    !
    real(ark) ::  feven,fodd,f0,fmax,h
      !
      h = xmax/real(Npoints,kind=ark)  !   integration step   
      feven=0         
      fodd =0
      f0   =f(0)
      fmax =f(Npoints)

     !
     !  sum of odd and even contributions 
     !
     do i = 1,npoints-2,2
        fodd   = fodd  + f(i  )
        feven  = feven + f(i+1)
     enddo
     !
     fodd   = fodd  + f(npoints-1)
     !
     si =  h/3.0_ark*( 4.0_ark*fodd + 2.0_ark*feven + f0 + fmax)

  end function  simpsonintegral


! integration with rectangular rules 
!                                      
  function integral_rect_ark(npoints,xmax,f) result (si) 
    integer(ik),intent(in) :: npoints
    !
    real(ark),intent(in) :: xmax,f(0:npoints)
    !
    real(ark) :: si,st
    !
    integer(ik) :: i
    !
    real(ark) ::  feven,fodd,f0,fmax,h
     !
     h = xmax/real(Npoints,kind=ark)  !   integration step   
     !
     !  sum of odd and even contributions 
     !
     si = sum(f)*h
     !
     !st  = simpsonintegral_ark(npoints,xmax,f)
     !
     !if (abs(si-st)>0.1) then 
     !   continue
     !endif
     !
  end function  integral_rect_ark


     !
     ! extrapolation at borders
     !
     subroutine interapolate_at_center(N,V,dv)
     !
     integer,intent(in) :: N
     real(rk),intent(inout) ::dv,V(-N:N)
     !
     integer            :: i1,i2,ipoints,fact
     real(rk)           :: x1,x2,d1t,a(N,N),b(N,1)
        !
        fact = 2
        if (abs(v(-2)/v(2)+1.0_rk)<sqrt(small_)) then 
           dv = 0 
           return 
        endif 
        !
        do i1 = 1,N
           !
           x1 = i1*rhostep
           !x2 = i1*rhostep
           !
           !
           b(i1,1) = V( i1)
           !
           do i2 = 2,N
             !
             a(i1 ,i2) = x1**((i2-1)*fact)
             !
           enddo
        enddo
        !
        !b(1,1) = V(0)
        !a(1,1) = 1.0_rk
        a(:,1) = 1.0_rk
        !a(1,2:N) = 0
        !
        !  lapack_gelss 
        ! 
        call lapack_gelss(a(:,:),b(:,:))
        !
        dv = b(1,1)

        !
   end subroutine interapolate_at_center



  subroutine numerov_lsqfit(npoints,x,f,phi)

    !
    integer(ik),intent(in) :: npoints
    real(rk),intent(in) :: f(0:npoints),x(0:npoints)
    real(rk),intent(out) :: phi(0:npoints)

    integer(ik),parameter :: pparam  = 6
    integer(ik),parameter :: qparam  = 9
    integer(ik),parameter :: itmax   = 200
    !
    real(rk)  :: stadev_best,ssq,stadev,f_t,alpha
    integer   :: iter,i,k,irow,icolumn
    integer   :: ierror

    real(rk) :: a(pparam,pparam),b(pparam,1),c(1:pparam),rho
    real(rk) :: acoef(0:qparam),eps(1:npoints)
    !
    real(rk) :: dx(pparam),rjacob(npoints,pparam)
      !
      c = 0 ; c(1) = -0.5_rk
      !
      stadev_best = 100.0_rk
      stadev = safe_max
      !
      iter = 0 
      !
      outer_loop: do while( iter<itmax .and. stadev>stadev_best )   
        !
        iter = iter + 1
        ssq=0
        !
        do i = 1,npoints
           !
           rho=x(i)
           !
           eps(i) = 0 ! c(1)/rho+c(2)
           !
           do k = -2,pparam-3
             eps(i) = eps(i)+c(3+k)*rho**k  !(2*k-4)
           enddo 
           !
           rjacob(i,1) = 1.0_rk/rho
           rjacob(i,2) = 1.0_rk
           !
           do k = -2,pparam-3
             rjacob(i,3+k) =rho**k   !(2*k-4)
           enddo 
           !
        enddo
        !
        eps(1:npoints) = eps(1:npoints) - f(1:npoints)
        !
        ssq = sum(eps(2:npoints)**2)
        !
        !----- form the a and b matrix ------
        !
        do irow=1,pparam   
          do icolumn=1,pparam
            a(irow,icolumn) = sum( rjacob(2:npoints,icolumn)*rjacob(2:npoints,irow) )
          enddo
            b(irow,1)=sum(eps(2:npoints)*rjacob(2:npoints,irow))
        enddo

        call MLlinur(pparam,pparam,a(1:pparam,1:pparam),b(1:pparam,1),dx,ierror)
        !
        if (ierror/=0) then 
          !
          write(out,"('numerov_lsqfit-MLlinur error ',i8)") ierror
          stop 'No fitted solution in numerov_lsqfit-MLlinur error'
          ! 
        endif 
        !
        c(:)=c(:)-dx(:)
        !
        !call lapack_gelss(a(:,:),b(:,:))
        !
        !c(:)=c(:)-b(:,1)
        !
        stadev=sqrt(ssq/real(npoints-pparam,kind=rk))
        !
      enddo outer_loop
      !
      if (iter==itmax) then 
         write (out,"('numerov_lsqfit: no convergence after ',i8,' iterations')") iter
         write (out,"('stadev = ',f18.8)") stadev
         write (out,"('eps = ')")
         write (out,"(40f18.8)") eps(1:min(40,size(eps)))
         stop 'numerov_lsqfit: - no convergence'
      endif 
      !
      f_t = max(1.0_rk+4.0_rk*c(1),1.0_rk-4.0_rk*c(1))
      !
      if (f_t<small_) then 
        !
        write(out,"('me_numer: sqrt(-1)')")
        stop 'me_numer: sqrt(-1)'
        !
      endif
      !
      alpha = 0.5_rk*(1.0_rk+sqrt(f_t))
      !
      acoef(0) = 1.0_rk
      !
      do i = 0,qparam
        !
        acoef(i) = 0 
        f_t = (alpha+real(i+1,rk))*(alpha+real(i+2,rk))-c(1)
        !
        do k = 1,pparam 
          !
          acoef(i) = acoef(i) + c(k)*acoef(i-k)
          !
        enddo
        !
        acoef(i) = acoef(i)/f_t
        !
      enddo
      !
      phi(0) = 0 
      !
      do i = 1,npoints
        !
        rho = x(i)
        !
        phi(i)= rho**alpha*acoef(1)
        !
        do k = 1,qparam-1
          phi(i) = phi(i)+acoef(k)*rho**(2*k)
        enddo 
        !
      enddo
      !

  end subroutine numerov_lsqfit


 end module me_numer



  subroutine ME_Fourier(vmax,maxorder,rho_b_,isingular,npoints_,numerpoints_,drho_,xton_,poten_,mu_rr_,icoord,&
                        iperiod,verbose,g_numerov,energy,KinOrder,PotOrder,ExtOrder)
   !
   integer(ik),intent(in) :: vmax,maxorder,npoints_,numerpoints_,isingular
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:vmax,0:vmax)
   real(ark),intent(out)    :: energy(0:vmax)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten_(0:npoints_),mu_rr_(0:npoints_),drho_(0:npoints_,3),xton_(0:npoints_,0:maxorder)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose   ! Verbosity level
   integer(ik),intent(in) :: iperiod
   integer(ik),intent(in) :: KinOrder, PotOrder, ExtOrder 
   !
   integer(ik),parameter  :: Factor_FF=30 ! factor to increase the Fourier basis set size 
   !
   real(ark)            :: rho,L,rhostep,potmin,rhostep_
   real(ark)            :: psipsi_t,characvalue,rho_b(2),cross_prod,factor,fval,df_t,step_scale
   !
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,n,imin,io_slot,kl,kr,p,fmax,npoints,i_,i1,i2,alloc_p, KineticOrder
   !
   real(ark),allocatable :: poten(:),mu_rr(:),rho_(:),xton(:,:)
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:),phi(:),&
                            dphi(:),vect(:,:)
   real(ark),allocatable :: phil_(:),phir_(:),dphil_(:),dphir_(:),phivphi_(:)
   real(ark),allocatable :: psi(:,:),dpsi(:,:),psi_(:,:),dpsi_(:,:)
   real(rk),allocatable ::  h(:,:),ener(:)
   !
   
   character(len=cl)    :: unitfname 
    !
    if (verbose>=2) write (out,"(/20('*'),' Fourier real functions primitive matrix elements calculations')")
     !
     call TimerStart('ME_Fourier')
     !
     ! global variables 
     !
     ! large grid 
     npoints   = numerpoints_
     !
     allocate(rho_kinet(0:npoints_),rho_poten(0:npoints_),rho_extF(0:npoints_),stat=alloc)
     call ArrayStart('rho-grids',alloc,size(rho_kinet),kind(rho_kinet))
     call ArrayStart('rho-grids',alloc,size(rho_poten),kind(rho_poten))
     call ArrayStart('rho-grids',alloc,size(rho_extF),kind(rho_extF))
     !
     fmax = vmax*factor_ff*iperiod
     !
     allocate(poten(0:npoints),mu_rr(0:npoints),stat=alloc)
     call ArrayStart('mu-poten',alloc,size(poten),kind(poten))
     call ArrayStart('mu-poten',alloc,size(mu_rr),kind(mu_rr))
     !
     allocate(xton(0:npoints,0:maxorder),stat=alloc)
     call ArrayStart('xton',alloc,size(xton),kind(xton))
     !
     allocate(psi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Fourier',alloc,size(psi),kind(psi))
     allocate(dpsi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Fourier',alloc,size(dpsi),kind(dpsi))
     !
     allocate(psi_(vmax+1,0:npoints_),stat=alloc)
     call ArrayStart('psi-phi-Fourier',alloc,size(psi_),kind(psi_))
     allocate(dpsi_(vmax+1,0:npoints_),stat=alloc)
     call ArrayStart('psi-phi-Fourier',alloc,size(dpsi_),kind(dpsi_))
     !
     allocate(h(fmax+1,fmax+1),ener(fmax+1),stat=alloc)
     call ArrayStart('h-Fourier',alloc,size(h),kind(h))
     call ArrayStart('h-Fourier',alloc,size(ener),kind(ener))
     !
     rho_b = rho_b_
     !
     ! step size 
     rhostep  = (rho_b(2)-rho_b(1))/real(npoints ,kind=ark)
     rhostep_ = (rho_b(2)-rho_b(1))/real(npoints_,kind=ark)
     !
     ! interpolation: mapping to the larger grid
     !
     step_scale = rhostep/rhostep_
     !
     if (npoints==npoints_) then 
       !
       poten = poten_
       mu_rr = mu_rr_
       !
       do lambda  = 0,maxorder
         xton(:,lambda) = xton_(:,lambda)
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
     ! periodic factor to produce cos( k nx ) and sin(k nx ) functions
     p = 1
     if (iperiod>0) then 
       p = iperiod
     endif
     !
     ! Do some reporting
     !
     if (verbose>=3) then 
         write (out,"('vmax = ',i8)") vmax
         write (out,"('fmax (Fourier size) = ',i8)") fmax
         write (out,"('maxorder = ',i8)") maxorder
         write (out,"('icoord = ',i4)") icoord
         write (out,"('rho_b (x) = ',2f12.4)") rho_b(1:2) !*180.0_ark/pi
         write (out,"('rhostep (x) = ',2f12.4)") rhostep  !*180.0_ark/pi
         write (out,"('periodicity = ',i8)") p
     endif 
     !
     
       
     potmin = huge(1.0_ark)
     !
     do i=0,npoints
        !
        if (poten(i)<potmin) then 
           imin = i
           potmin = poten(i)
        endif
        !
     enddo
     !
     if (imin<0.or.imin>npoints) then 
         write(out,"('ME_Fourier: pot_eff has no minimum',i8)") 
         stop 'ME_Fourier: pot_eff has no minimum'
     endif 
     !
     ! define the rho-type coordinate 
     !
     rho_kinet(:) = drho_(:,1)
     rho_poten(:) = drho_(:,2)
     rho_extF(:)  = drho_(:,3)
     !
     inquire(iolength=rec_len) rho_kinet(:),rho_poten(:)
     !
     ! Print out
     !
     if (verbose>=3) then 
        write(out,"('grid values (i,rho,rho_kinet,rho_poten,poten, mu_rr): ')") 
        do i_=0,npoints_,2
          i = int( real(i_,ark)/step_scale )
          rho = rho_b(1)+real(i_,kind=ark)*rhostep_
          write(out,"(i8,3f14.6,2g14.6)") i_,rho,rho_kinet(i_),rho_poten(i_),poten(i),mu_rr(i)
        enddo 
     endif     
     write(unitfname,"('Numerov basis set # ',i6)") icoord
     call IOStart(trim(unitfname),io_slot)
     !
     open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
     !
     !
     ! Matrix elements
     !
     ! box size: 
     !
     L = (rho_b(2)-rho_b(1))*0.5_ark
     !
     ! Build primitive Fourier functions on a grid
     !
     if (verbose>=4) write(out,"('   Build primitive Fourier functions and 1D Hamiltonian on a grid ...')")
     !
     !$omp parallel private(phil,phir,dphil,dphir,alloc_p,phivphi) shared(h) 
     allocate(phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints),phivphi(0:npoints),stat=alloc_p)
     if (alloc_p/=0) then 
       write (out,"('phil,phir - out of memory')")
       stop 'phil,phir - out of memory'
     endif 
     !
     !$omp do private(vl,kl,i,rho,vr,kr,psipsi_t) schedule(dynamic)
     do vl = 0,fmax
        kl = (vl+1)/2*p
        !
        do i=0,npoints
           !
           rho = real(i,kind=ark)*rhostep
           !
           if (vl==0) then 
             phil(i)  = sqrt(0.5_ark/L)
             dphil(i) = 0
           elseif (mod(vl,2)==0) then
             phil(i)  = sqrt(1.0_ark/L)*cos(real(kl,ark)*pi*rho/L)
             dphil(i) =-sqrt(1.0_ark/L)*sin(real(kl,ark)*pi*rho/L)*real(kl,ark)*pi/L
           else
             phil(i)  = sqrt(1.0_ark/L)*sin(real(kl,ark)*pi*rho/L)
             dphil(i) = sqrt(1.0_ark/L)*cos(real(kl,ark)*pi*rho/L)*real(kl,ark)*pi/L
           endif
           !
        enddo
        !
        do vr = vl,fmax
            kr = (vr+1)/2*p
            !
            do i=0,npoints
               !
               rho = real(i,kind=ark)*rhostep
               !
               if (vr==0) then 
                 phir(i)  = sqrt(0.5_ark/L)
                 dphir(i) = 0
               elseif (mod(vr,2)==0) then
                 phir(i)  = sqrt(1.0_ark/L)*cos(real(kr,ark)*pi*rho/L)
                 dphir(i) =-sqrt(1.0_ark/L)*sin(real(kr,ark)*pi*rho/L)*real(kr,ark)*pi/L
               else
                 phir(i)  = sqrt(1.0_ark/L)*sin(real(kr,ark)*pi*rho/L)
                 dphir(i) = sqrt(1.0_ark/L)*cos(real(kr,ark)*pi*rho/L)*real(kr,ark)*pi/L
               endif
               !
            enddo
            !
            ! check orthogonality (normalisation)
            phivphi(:) = phil(:)*phir(:)
            psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
            !
            ! Here we prepare integrals of the potential 
            ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
            ! obtained above by the Numerov
            !
            phivphi(:) = phil(:)*poten(:)*phir(:)
            !
            h(vl+1,vr+1) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
            !
            ! momenta-quadratic part 
            !
            phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)
            !
            psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
            !
            ! Add the diagonal kinetic part to the tested mat. elem-s
            !
            h(vl+1,vr+1) = h(vl+1,vr+1) - 0.5_ark*psipsi_t
            !
            h(vr+1,vl+1) = h(vl+1,vr+1)
            !
        enddo
     enddo
     !$omp end do
     !
     deallocate (phil,phir,dphil,dphir,phivphi)
     !$omp end parallel 
     !
     if (verbose>=4) write(out,"('   Diagonalize the 1D primitive Hamiltonian ...')")
     call lapack_syev(h,ener)
     !
     if (verbose>=4) write(out,"('ZPE = ',f16.6)") ener(1)
     !
     energy(0:vmax) = ener(1:vmax+1)-ener(1)
     !
     ! Schmidt orthogonalization to make eigenvectors orthogonal in ark
     !
     allocate(vect(fmax+1,fmax+1),stat=alloc)
     call ArrayStart('h-vect',alloc,size(vect),kind(vect))
     !
     vect = h
     !
     !omp parallel do private(vl,cross_prod,factor,vr) shared(vect) schedule(dynamic)
     do vl =  1,vmax+1
       !
       cross_prod = sum(vect(:,vl)*vect(:,vl))
       !
       factor = 1.0_ark/sqrt(cross_prod)
       !
       vect(:,vl) = vect(:,vl)*factor
       !
       do vr = 1,vl-1
         !
         cross_prod = sum(vect(:,vl)*vect(:,vr))
         !
         vect(:,vl) = vect(:,vl)-cross_prod*vect(:,vr)
         !
         cross_prod = sum(vect(:,vl)*vect(:,vl))
         !
         factor = 1.0_ark/sqrt(cross_prod)
         !
         vect(:,vl) = vect(:,vl)*factor
         ! 
       enddo
       !
     enddo
     !omp end parallel do
     !
     !
     if (verbose>=4) then 
       write (out,"(/' Fourier-optimized energies are:')") 
       !
       do vl=0,vmax   
         write (out,"(i8,f18.8)") vl,energy(vl)
       enddo
     endif
     !
     if (verbose>=4) write(out,"('   Transform primitive to contracted on the grid ...')")
     !
     !$omp parallel private(phi,dphi,alloc_p) shared(Psi,dPsi)
     allocate(phi(fmax+1),dphi(fmax+1),stat=alloc_p)
     if (alloc_p/=0) then 
       write (out,"('phi,dphi - out of memory')")
       stop 'phi,dphi - out of memory'
     endif 
     !
     !$omp do private(i,rho,vl,kl) schedule(dynamic)
     do i=0,npoints
        !
        rho = real(i,kind=ark)*rhostep
        do vl = 0,fmax
           kl = (vl+1)/2*p
           !
           if (vl==0) then 
             phi(vl+1)  = sqrt(0.5_ark/L)
             dphi(vl+1) = 0
           elseif (mod(vl,2)==0) then
             phi(vl+1)  = sqrt(1.0_ark/L)*cos(real(kl,ark)*pi*rho/L)
             dphi(vl+1) =-sqrt(1.0_ark/L)*sin(real(kl,ark)*pi*rho/L)*real(kl,ark)*pi/L
           else
             phi(vl+1)  = sqrt(1.0_ark/L)*sin(real(kl,ark)*pi*rho/L)
             dphi(vl+1) = sqrt(1.0_ark/L)*cos(real(kl,ark)*pi*rho/L)*real(kl,ark)*pi/L
           endif
           !
        enddo
        !
        Psi (1:vmax+1,i)  = matmul(transpose(vect(1:fmax+1,1:vmax+1)),phi(1:fmax+1))
        DPsi(1:vmax+1,i)  = matmul(transpose(vect(1:fmax+1,1:vmax+1)),dphi(1:fmax+1))
        !
     enddo
     !$omp end do
     !
     deallocate (phi,dphi)
     !$omp end parallel 
     !
     !
     ! dump the eigenfunction
     if (npoints_==npoints) then 
          Psi_  = Psi
          DPsi_ = DPsi
     else
       !
       !$omp parallel do private(i_,i) shared(Psi_,dPsi_) schedule(dynamic)
       do i_ = 0,npoints_
          !
          i = nint( real(i_,ark)/step_scale )
          !
          Psi_(1:vmax+1,i_) = Psi(1:vmax+1,i)
          DPsi_(1:vmax+1,i_) = DPsi(1:vmax+1,i)
          !
       enddo
       !$omp end parallel do
       !
     endif
     !
     deallocate(vect)
     call ArrayStop('h-vect')
     !
     ! store basis functios 
     do vl = 0,vmax
        !
        write (io_slot,rec=vl+1) (Psi_(vl+1,i),i=0,npoints_),(dPsi_(vl+1,i),i=0,npoints_)
        !
     enddo
     !
     ! Build contracted Fourier functions on a grid
     !
     if (verbose>=4) write(out,"('   Generate matrix elements of elements of 1D Hamiltonian, x^n, p x^n, p^2 x^n  ...')")
     !
     !$omp parallel private(phil_,phir_,dphil_,dphir_,alloc_p,phivphi_) shared(h,g_numerov)
     allocate(phil_(0:npoints_),phir_(0:npoints_),dphil_(0:npoints_),dphir_(0:npoints_),phivphi_(0:npoints_),stat=alloc_p)
     if (alloc_p/=0) then 
       write (out,"('phi_ - out of memory')")
       stop 'phi_ - out of memory'
     endif
     !
     ! 
     !$omp do private(vl,i,rho,vr,psipsi_t,lambda) schedule(dynamic)
     do vl = 0,vmax
        !
        phil_(:)  =  Psi_(vl+1,:)
        dphil_(:) = dPsi_(vl+1,:)
        !
        !write (io_slot,rec=vl+1) (phil_(i),i=0,npoints_),(dphil_(i),i=0,npoints_)
        !
        do vr = vl,vmax
            !
            phir_  = Psi_(vr+1,:)
            dphir_ = dPsi_(vr+1,:)
            !
            ! check orthogonality (normalisation)
            phivphi_(:) = phil_(:)*phir_(:)
            psipsi_t = simpsonintegral_ark(npoints_,rho_b_(2)-rho_b_(1),phivphi_)
            !
            ! Here we prepare integrals of the potential 
            ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
            ! obtained above by the Numerov
            !
            phivphi_(:) = phil_(:)*poten_(:)*phir_(:)
            !
            h(vl+1,vr+1) = simpsonintegral_ark(npoints_,rho_b_(2)-rho_b_(1),phivphi_)
            !
            ! momenta-quadratic part 
            !
            phivphi_(:) =-dphil_(:)*mu_rr(:)*dphir_(:)
            !
            psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
            !
            ! Add the diagonal kinetic part to the tested mat. elem-s
            !
            h(vl+1,vr+1) = h(vl+1,vr+1) - 0.5_ark*psipsi_t
            !
            h(vr+1,vl+1) = h(vl+1,vr+1)
            !
            psipsi_t = 0 
            !
            write(*,*) "test line 1164 maxorder ", maxorder
            do lambda = 0,maxorder
              !
              ! momenta-free part in potential part
              !
              if (lambda==0) then 
                 phivphi_(:) = phil_(:)*phir_(:)
              else
                 phivphi_(:) = phil_(:)*rho_poten(:)**lambda*phir_(:)
              endif
              !
              if(lambda > PotOrder) then 
                 g_numerov(0,lambda,vl,vr) = 0.0_ark
              else 
                 g_numerov(0,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
              endif
              !
              ! external field expansion
              !
              if (lambda==0) then 
                 phivphi_(:) = phil_(:)*phir_(:)
              else
                 phivphi_(:) = phil_(:)*rho_extF(:)**lambda*phir_(:)
              endif
              !
              if(lambda > ExtOrder) then
                g_numerov(3,lambda,vl,vr) = 0.0_ark
              else
                g_numerov(3,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
              endif
              !  
              if (vl/=vr) g_numerov(3,lambda,vr,vl) = g_numerov(3,lambda,vl,vr)
              !
              ! momenta-free in kinetic part 
              !
              phivphi_(:) = phil_(:)*xton(:,lambda)*phir_(:)
              !
              if(lambda > KinOrder) then
                g_numerov(-1,lambda,vl,vr) = 0.0_ark  
              else
                g_numerov(-1,lambda,vl,vr) = simpsonintegral_ark(npoints_,rho_b(2)-rho_b(1),phivphi_)
              endif
              !
              ! We also control the orthogonality of the basis set 
              !
              if (lambda==0) psipsi_t = g_numerov(0,lambda,vl,vr)
              !
              if (vl/=vr) g_numerov(-1:0,lambda,vr,vl) = g_numerov(-1:0,lambda,vl,vr)
              !
              ! momenta-quadratic part 
              !
              phivphi_(:) =-dphil_(:)*xton(:,lambda)*dphir_(:)
              !
              if(lambda > KinOrder) then
                g_numerov(2,lambda,vl,vr) = 0.0_ark
              else
                g_numerov(2,lambda,vl,vr) = simpsonintegral_ark(npoints_,rho_b(2)-rho_b(1),phivphi_)
              endif
              !
              if (vl/=vr) g_numerov(2,lambda,vr,vl) = g_numerov(2,lambda,vl,vr)
              !
              ! momenta-linear part:
              ! < vl | d/dx g(x) | vr > = - < vr | g(x) d/dx | vl >
              !
              phivphi_(:) = phil_(:)*xton(:,lambda)*dphir_(:)
              !
              if(lambda > KinOrder) then 
                g_numerov(1,lambda,vl,vr) = 0.0_ark
              else
                g_numerov(1,lambda,vl,vr) = simpsonintegral_ark(npoints_,rho_b(2)-rho_b(1),phivphi_)
              endif 
              !
              if (vl/=vr) then
                !
                phivphi_(:) = dphil_(:)*xton(:,lambda)*phir_(:)
                !
                if(lambda > KinOrder) then
                  g_numerov(1,lambda,vr,vl) = 0.0_ark
                else 
                  g_numerov(1,lambda,vr,vl) = simpsonintegral_ark(npoints_,rho_b(2)-rho_b(1),phivphi_)
                endif
                !
              endif 
              !
              if (verbose>=5) then 
                  write(out,"('g_numerov(-1,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(-1,lambda,vl,vr)
                  write(out,"('g_numerov(0,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(0,lambda,vl,vr)
                  write(out,"('g_numerov(1,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(1,lambda,vl,vr)
                  write(out,"('g_numerov(2,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(2,lambda,vl,vr)
                  write(out,"('g_numerov(3,',i4,i4,i4,') = ',f18.8)") lambda,vl,vr,g_numerov(3,lambda,vl,vr)
                  if (vl/=vr) then 
                    write(out,"('g_numerov(-1,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(-1,lambda,vr,vl)
                    write(out,"('g_numerov(0,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(0,lambda,vr,vl)
                    write(out,"('g_numerov(1,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(1,lambda,vr,vl)
                    write(out,"('g_numerov(2,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(2,lambda,vr,vl)
                    write(out,"('g_numerov(3,',i4,i4,i4,') = ',f18.8)") lambda,vr,vl,g_numerov(3,lambda,vr,vl)
                  endif 
              endif 
              !
            enddo 
            !
        enddo
        !
        if (verbose>=5) then 
           !
           !write (out,"('v = ',i8,f18.8)") vl,h(vl+1,vl+1)-h(1,1)
           !$omp critical
           do i=0,npoints 
              write(out,"(i8,2f18.8,' || ',1x,i8)") i,phil_(i),dphil_(i),vl
           enddo
           !$omp end critical
           !
        endif 
        !
     enddo
     !$omp end do
     !
     deallocate (phil_,phir_,dphil_,dphir_,phivphi_)
     !$omp end parallel 
     !
     deallocate(rho_kinet,rho_poten,rho_extF)
     !
     call ArrayStop('rho-grids')
     !
     deallocate(poten,mu_rr)
     call ArrayStop('mu-poten')
     !
     deallocate(xton)
     call ArrayStop('xton')
     !
     deallocate(h,ener)
     call ArrayStop('h-Fourier')
     deallocate(psi,dpsi,psi_,dpsi_)
     call ArrayStop('psi-phi-Fourier')
     !
     if (verbose>=2) write (out,"(/40('*')/)")
     !
     call TimerStop('ME_Fourier')
     !
  end subroutine ME_Fourier

module me_str
  use accuracy
  implicit none
  private

  public ik, rk, out
  public ME_morse,ME_harmonic


  integer(ik), parameter :: verbose     = 4                          ! Verbosity level

!
  contains
!


  !********************************************************************!
  !                                                                    
  !     Stretching matrix elements                                    
  !     based on the Efremov recursion procedure                       
  !                                                                    
  !********************************************************************!
  !********************************************************************!
  !                                                                    !
  subroutine ME_morse(vstrmax, maxorder, rmk, amorse, g_str, energy ) 
  !                                                                    !
  !********************************************************************!

    integer(ik), intent(in)       ::  vstrmax
    integer(ik), intent(in)       ::  maxorder
    real(ark), intent(in)         ::  rmk    ! Morse parameter RMK
    real(ark), intent(in)         ::  amorse ! Morse parameter a   
    real(ark), intent(out)        ::  g_str(-1:3,0:maxorder,0:vstrMax,0:vstrmax)
    real(ark), intent(out)        ::  energy(0:vstrmax)

    real(ark)                     ::  over(0:vstrmax,0:vstrmax)
    real(ark) , allocatable       ::  jei_mat0(:,:,:),jei_mat1(:,:,:)
    real(ark) , allocatable       ::  y(:,:)
    integer(ik)                   ::  vn,vm,j,lambda,nu,alloc
    real(ark)                     ::  delta,sumterm,alpha

    !
    ! the routine starts here 
    !

    if (verbose>=1) write (out,"(20('*'),' Morse matrix elements calculations')")

    !
    ! calculate one-dimensional overlap integrals
    !
    call ovrlap(vstrmax, rmk, amorse, over)


    allocate( jei_mat0(0:maxorder+2,0:vstrmax,0:vstrmax), &
              jei_mat1(0:maxorder+2,0:vstrmax,0:vstrmax), & 
              y(0:maxorder+2,0:maxorder+2),stat=alloc)
    if (alloc/=0) then 
      write (6,"('ME_morse: jei_mat,y - out of memory')")
      stop 'ME_morse: jei_mat,y - out of memory'
    endif  


    ! Calculate  matrix elements < morse (a nbond  v_i) | y^k | morse (c nbond  v_j) > 
    ! and                        < morse (b nbond  v_i) | y^k | morse (c nbond  v_j) > 
    ! y:=(k)->simplify((-1)^k*n!/kappa^k/(n-k)!/k!/2^k);

     call calc_jmatrix_xton      (maxorder,vstrmax,rmk,jei_mat0)
     call calc_jmatrix_xton_dxtom(maxorder,vstrmax,rmk,jei_mat1)

     call y2x_expantion(maxorder,rmk,y)

     do vn = 0,vstrmax 
       alpha = 2.0_ark*rmk-2.0_ark*vn-1.0_ark 
       do vm = vn,vstrmax
         j  = vm-vn
         delta  = 0
         if ( vn==vm ) delta = 1.0_ark 
         !
         do lambda = 0, maxorder
           sumterm = 0  
           do nu = 0,lambda-1 
             sumterm = sumterm + y(lambda,nu+1)*jei_mat0(nu,vn,j)/amorse 
           enddo 
           !
           g_str(0,lambda,vm,vn) = delta+sumterm*over(vm,vn)
           if (vn/=vm) g_str(0,lambda,vn,vm) = g_str(0,lambda,vm,vn)

         enddo ! lambda 
         !
         if ( vm==vn ) then 
           g_str(1,0,vn,vn) = 0
           !
           do lambda = 1, maxorder
             g_str(1,lambda,vn,vn) = lambda*amorse*0.5_ark*( g_str(0,lambda,vn,vn) - g_str(0,lambda-1,vn,vn) )
           enddo 
           !
         else 
           g_str(1,0,vm,vn) = 0.5_ark*over(vm,vn)
           g_str(1,0,vn,vm) = -g_str(1,0,vm,vn)

           do lambda = 1, maxorder
             sumterm = 0.5_ark
             do nu = 1,lambda 
               sumterm = sumterm + y(lambda,nu)*(                            &
                     jei_mat0(nu,vn,j)*0.5_ark                                &
                   - jei_mat0(nu-1,vn,j)*( real(vn,rk)+0.5_ark*alpha)         &
                   + jei_mat1(nu-1,vn,j)*( real(vn,rk)+       alpha) )
             enddo 
             !
             g_str(1,lambda,vm,vn) = sumterm*over(vm,vn)
             g_str(1,lambda,vn,vm) =-g_str(1,lambda,vm,vn)                   &
                                   -lambda*amorse*( g_str(0,lambda-1,vm,vn)     &
                                   -g_str(0,lambda,vm,vn) )

           enddo ! lambda 

         endif 
         !
         do lambda = 0,maxorder
            sumterm = 0
            do nu = 0,lambda+1
               sumterm = sumterm + ( y(lambda+2,nu+1)             &
                         - (1.0_ark-(0.5_ark*alpha/rmk)**2)*y(lambda,nu+1) )*  &
                        jei_mat0(nu,vn,j)
            enddo
            !
            g_str(2,lambda,vm,vn) = delta*0.25_ark*(amorse*alpha)**2 + rmk*rmk*amorse*sumterm*over(vm,vn)
            !
            if ( lambda>=1 )         &
                 g_str(2,lambda,vm,vn) = g_str(2,lambda,vm,vn) &
               + real(lambda,rk)*amorse*( g_str(1,lambda-1,vm,vn)-g_str(1,lambda,vm,vn) )
            !
            if (vn.ne.vm) g_str(2,lambda,vn,vm) = g_str(2,lambda,vm,vn)
            ! 
            !if ( lambda.le.maxorder-2 ) then 
            !   sumterm = rmk*rmk*amorse*amorse*g_str(0,lambda+2,vn,vm)-                 &
            !            (4.0_ark*rmk**2-alpha**2)*amorse**2*0.25_ark*g_str(0,lambda,vn,vm)
            !endif 
         enddo ! lambda 
        enddo      !----- vm   -----
      enddo      !----- vn   -----
      !
      ! p=-1 is an optional dimemnsion needed to distinguish different variable used in the kinetic and potential 
      ! parts. It is absolete in ME_morse, therefore we copy p=0 to p=-1
      !
      g_str(-1,:,:,:) = g_str(0,:,:,:)
      g_str( 3,:,:,:) = g_str(0,:,:,:)
      !
      do nu = 0,vstrmax
         !
         energy(nu) = -g_str(2,0,nu,nu)+amorse**2*rmk**2*g_str(0,2,nu,nu)
         !
      enddo 
      !
      energy(:) = energy(:) - energy(0)
      !
      if (allocated(jei_mat0)) deallocate( jei_mat0 )
      if (allocated(jei_mat1)) deallocate( jei_mat1 )
      if (allocated(y       )) deallocate( y )
     !
  end subroutine  ME_morse

  !
  ! j-matrix
  ! based on the efremov's recursion procedure 
  !
  subroutine calc_jmatrix_xton(maxorder,vstrmax,rmk, jei_mat)

    integer(ik), intent(in)      :: vstrmax
    integer(ik), intent(in)      :: maxorder
    real(ark), intent(in)         :: rmk    ! Morse parameter RMK
    real(ark), intent(out)        :: jei_mat(0:maxorder+2,0:vstrmax,0:vstrmax)

    real(ark)                     :: sumterm,expon,alpha,jlamterm,gterm
    integer(ik)                  :: lambda,j,vn,vm,m0


    do vn =  0, vstrmax 
       do vm = vn, vstrmax
          !
          j = vm-vn
          jei_mat(0,vn,j) = 1.0_ark
          do lambda = 1,maxorder+2
             alpha = 2.0_ark*rmk-2.0_ark*real(vn,kind=rk)-1.0_ark
             sumterm = 0.0_ark
             if ( lambda<j ) then 
                jlamterm=exp(faclog(j+lambda)-faclog(j-lambda-1_ik)+faclog(vn))
                do m0 = 0, min(lambda,vn)
                   gterm = gfunctlog(alpha+real(vn-j+1_ik,kind=ark),lambda-m0)
                   !   
                   expon= gterm + faclog(j-lambda+m0-1_ik)    &
                                - faclog(lambda-m0)-faclog(m0)-faclog(j+m0)-faclog(vn-m0)
                   sumterm = sumterm + real((-1)**m0,rk)*exp(expon)
                enddo 
             else    
                jlamterm=exp(faclog(lambda+j)+faclog(lambda-j)+faclog(vn))
                do m0 = 0, min(lambda-j,vn)
                   !
                   gterm = gfunctlog(alpha+real(vn-j+1_ik,kind=ark),lambda-m0)
                   expon= gterm-faclog(m0)-faclog(lambda-j-m0)-faclog(lambda-m0)  &
                               -faclog(j+m0)-faclog(vn-m0)
                   sumterm = sumterm + exp(expon)
                   !
                enddo ! --- m0
             endif 
             !
             jei_mat(lambda,vn,j) = sumterm*jlamterm
          enddo 
       enddo    
    enddo    

  end subroutine calc_jmatrix_xton



  !
  !  j-matrix for derivatives
  !  based on the efremov's recursion procedure 
  !
  subroutine calc_jmatrix_xton_dxtom(maxorder,vstrmax, rmk, jei_mat)

    integer(ik), intent(in)      ::  vstrmax
    integer(ik), intent(in)      ::  maxorder
    real(ark), intent(in)         :: rmk    ! Morse parameter RMK
    real(ark), intent(out)        :: jei_mat(0:maxorder+2,0:vstrmax,0:vstrmax)

    real(ark)                     ::  sumterm,expon,alpha,jlamterm,gterm
    integer(ik)                  ::  lambda,j,vn,vm,m0

    do vn =  0, vstrmax 
      do vm = vn, vstrmax
         j = vm-vn
         jei_mat(0,vn,j) = 0
         do lambda = 1,maxorder+2
            alpha = 2.0_ark*rmk-2.0_ark*vn-1.0_ark
            sumterm = 0
            if ( lambda<j ) then 
               !
               jlamterm =exp(faclog(j+lambda)-faclog(j-lambda-1_ik)+faclog(vn))
               !
               do m0 = 0, min(lambda-1,vn-1)
                  gterm = gfunctlog(alpha+real(vn-j+1_ik,ark),lambda-m0-1_ik)
                  expon= gterm+faclog(j-lambda+m0-1_ik)  &
                              -faclog(lambda-m0-1_ik)-faclog(m0)-faclog(j+m0+1_ik)-faclog(vn-m0-1_ik)
                  !
                  sumterm = sumterm  - real((-1)**m0,rk)*exp(expon)
                  !
               enddo 
               !
            else 
               jlamterm=-exp(faclog(lambda+j)+faclog(lambda-j)+faclog(vn))
               do m0 = 0, min(lambda-j,vn-1)
                  gterm = gfunctlog(alpha+real(vn-j+1_ik,ark),lambda-m0-1_ik)
                  expon= gterm - faclog(m0) - faclog(lambda-j-m0) - faclog(lambda-m0-1_ik) &
                               - faclog(j+m0+1_ik) - faclog(vn-m0-1_ik)
                  sumterm = sumterm + exp(expon)
               enddo
            endif 
           jei_mat(lambda,vn,j) = sumterm*jlamterm
         enddo 
      enddo  
    enddo    

  end subroutine calc_jmatrix_xton_dxtom
  !
  function gfunctlog(x,n)  result (v)
    real(ark)  ,intent(in) ::  x
    integer(ik),intent(in) ::  n
    real(ark)              :: v 
    integer(ik)            :: l
    real(ark)              :: gterm

    gterm  = 0 
    if ( n>=0 ) then 
       do l = 1,n 
         gterm  = gterm+log(x+real(l,ark)-1.0_ark)
       enddo 
    else
       do l = 1,abs(n) 
          gterm  = gterm-log(x+real(n+l,ark)-1.0_ark)
       enddo
     endif 

     v = gterm

  end function gfunctlog

  !
  subroutine y2x_expantion(maxorder,rmk,y)
    integer(ik), intent(in) :: maxorder
    real(ark), intent(in)    :: rmk    ! Morse parameter RMK
    real(ark),intent(out)    :: y(0:maxorder+2,0:maxorder+2)
    real(ark)                :: expon
    integer(ik)             :: lambda,nu

    y = 0
    !
    !      y:=(k)->simplify((-1)^k*n!/kappa^k/(n-k)!/k!/2^k);
    !
     do lambda = 0, maxorder+2 
        y(lambda,0) = 1.0_ark
        do nu = 1,lambda
           expon  = faclog(lambda)-faclog(nu)-faclog(lambda-nu)
           y(lambda,nu) = 1.0_ark/(-2.0_ark*rmk)**nu*exp(expon)
        enddo
     enddo 
  end subroutine y2x_expantion
  !



  !
  ! calculate one-dimensional overlap integrals
  !
  subroutine ovrlap(vstrmax, rmk, amorse, over)

  integer(ik),intent(in) :: vstrmax  
  real(ark),intent(in) ::  rmk,amorse
  real(ark),intent(out) ::  over(0:vstrmax,0:vstrmax)

  real(ark) ::  a,rifm,rifn,rpl,rxl,rp
  integer(ik) :: n,mm,j,isi,i

    a=amorse
   
    do n = 0,vstrmax
       do mm = n+1,vstrmax
          j = mm-n
          isi=(-1)**j
          rifm=faclog(mm)
          rifn=faclog(n)
          rpl= rifm - rifn  &
             + log(2.0_ark*rmk-2.0_ark*real(mm,kind=rk)-1.0_ark) &
             + log(2.0_ark*rmk-2.0_ark*real(n ,kind=rk)-1.0_ark)

          rxl=0.0_ark
          do i=1,j
             rxl=rxl+log( 2.0_ark*rmk-real(mm-i+1,kind=rk) )
          enddo 

          rpl=log(a)+0.5_ark*(rpl - rxl)
          rp =real(isi,kind=rk)*exp(rpl)
          over(n,mm) = rp
          over(mm,n) = rp
       enddo 
       over(n,n) = a*(2.0_ark*rmk-2.0_ark*real(n,kind=rk)-1.0_ark)
    enddo 
  end subroutine ovrlap


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
  !
  !********************************************************************!
  !                                                                    
  !     Stretching matrix elements                                    
  !     for harmonic eigenfunctions 
  !                                                                    
  !********************************************************************!
  !********************************************************************!
  !                                                                    !
  subroutine ME_harmonic(vstrmax, maxorder, coeff_norm, g_str, energy) 
  !                                                                    !
  !********************************************************************!

    integer(ik), intent(in)          ::  vstrmax
    integer(ik), intent(in)          ::  maxorder
    real(ark),   intent(in)          ::  coeff_norm

    real(ark), intent(out)           ::  g_str(-1:3,0:maxorder,0:vstrmax,0:vstrmax)
    real(ark), intent(out)           ::  energy(0:vstrmax)

    real(ark),allocatable            ::  g_t(:,:,:,:)
 
    integer(ik)                      ::  v1,v2,v0,n,p,alloc
    real(ark)                        ::  matrel,u,f,fp

    !
    ! the routine starts here 
    !

    if (verbose>=1) write (out,"(20('*'),' Harmonic matrix elements calculations')")

    ! Temporaly matrix g_t similar to g_str, but with extended (+1) size
    !
    allocate( g_t(0:2,0:maxorder,0:vstrmax+maxorder+2,0:vstrmax+maxorder+2),stat=alloc)
    if (alloc/=0) then 
      write (6,"('ME_harmonic: g_t - out of memory')")
      stop 'ME_harmonic: g_t - out of memory'
    endif  


    g_str = 0
    g_t   = 0
    !
    ! First we define basic matrix elements <v1|q|v2+/-1> and <v1|p|v2+/-1> 
    !
    g_t(0,0,0,0) = 1.0_ark
    g_t(0,1,0,1) = sqrt( 0.5_ark)
    g_t(1,0,0,1) = sqrt( 0.5_ark)
    !
    do v0 = 1,vstrmax+maxorder
       !
       u = sqrt( 0.5_ark*real(v0,kind=ark) )

       ! <v0|1|v0> 
       !
       g_t(0,0,v0,v0  ) = 1.0_ark
       ! <v0|q|v0-1> 
       !
       g_t(0,1,v0  ,v0-1) = u
       g_t(0,1,v0-1,v0  ) = u


       ! <v0|p|v0-1> 
       !
       g_t(1,0,v0  ,v0-1) =-u
       g_t(1,0,v0-1,v0  ) = u

    enddo  ! --- v0

    ! coordinate part
    !
    do n = 2,maxorder
       do v1 = 0,vstrmax+maxorder
          do v2 = max(v1-n,0),min(v1+n,vstrmax+maxorder+2)
             !
             ! applying the summation rule
             !
             matrel = 0.0_ark
             do v0 = max(v1-1,0),min(v1+1,vstrmax+maxorder+2)
                matrel= matrel + g_t(0,1,v1,v0)*g_t(0,n-1,v0,v2)
                continue
             enddo 
             g_t(0,n,v1,v2) = matrel
          enddo 
       enddo 
    enddo 

    ! momenta part
    !
    do n = 1,maxorder
       do v1 = 0,vstrmax+maxorder+1
          do v2 = max(v1-1-n,0),min(v1+1+n,vstrmax+maxorder+2)
             !
             ! applying the summation rule
             !
             matrel = 0.0_ark
             do v0 = max(v2-1,0),min(v2+1,vstrmax+maxorder+2)
                matrel=   matrel + g_t(0,n,v1,v0)*g_t(1,0,v0,v2)
             enddo 
             g_t(1,n,v1,v2) = matrel
          enddo 
       enddo 
    enddo 
   
    do n = 0,maxorder
       do v1 = 0,vstrmax+maxorder+2
          do v2 = max(v1-2-n,0),min(v1+2+n,vstrmax+maxorder+2)
             !
             ! applying the summation rule
             !
             matrel = 0.0_ark
             do v0 = max(v1-1,0),min(v1+1,vstrmax+maxorder+2)
                matrel=   matrel + g_t(1,0,v1,v0)*g_t(1,n,v0,v2)
             enddo 
             g_t(2,n,v1,v2) = matrel
          enddo 
       enddo 
    enddo 
    !
    ! if the coordinates are not normal then we need to apply a normalization coefficient coeff_norm
    !
    if (verbose>=4) write (out,"('p,n,v1,v2,g_str: ')")
    !
    g_str(0,0,0,0) = 1.0_ark
    !
    do p = 0,2
       !
       do n = 0,maxorder
          !
          f = 1.0_ark
          !
          if (n/=p) f = coeff_norm**(n-p)
          !
          g_str(p,n,0:vstrmax,0:vstrmax) = g_t(p,n,0:vstrmax,0:vstrmax)*f
          !
          if (verbose>=4) then 
              do v1 = 0,vstrmax
                 do v2 = 0,vstrmax
                    !
                    ! g_str(p,n,v1,v2) = g_t(p,n,v1,v2)*f
                    !
                    write (out,"(4i8,f20.8)") p,n,v1,v2,g_str(p,n,v1,v2)
                    !
                 enddo
              enddo
          endif 
          !
       enddo 
    enddo 
    !
    ! p=-1 stands for the kinetic part, but without momenta
    !
    g_str(-1,:,:,:) = g_str(0,:,:,:)
    g_str( 3,:,:,:) = g_str(0,:,:,:)
    !
    do n = 0,vstrmax
       !
       energy(n) = -g_str(2,0,n,n)+1.0_ark/coeff_norm**4*g_str(0,2,n,n)
       !
    enddo 
    !
    energy(:) = energy(:) - energy(0)
    !
    !
    ! Clean up !
    !
    deallocate(g_t)  

    if (verbose>=1) write (out,"(20('*'),' Harmonic matrix elements calculations/end')")



  end subroutine  ME_harmonic

end module me_str

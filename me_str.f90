module me_str
  use accuracy
  use timer
  use me_numer
  implicit none
  private

  public ik, rk, out
  public ME_morse,ME_harmonic,ME_Laguerre


  integer(ik), parameter :: verbose     = 3                          ! Verbosity level

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
    !
    if (verbose>=1) write (out,"(20('*'),' Harmonic matrix elements calculations/end')")
  !
  end subroutine  ME_harmonic



  subroutine lf_function ( m, n, alpha, x, cx )
  !*****************************************************************************80
  !
  !LF_FUNCTION evaluates the Laguerre function Lf(n,alpha,x).
  !
  !  Recursion:
  !
  !    Lf(0,ALPHA,X) = 1
  !    Lf(1,ALPHA,X) = 1+ALPHA-X
  !
  !    Lf(N,ALPHA,X) = (2*N-1+ALPHA-X)/N * Lf(N-1,ALPHA,X) 
  !                      - (N-1+ALPHA)/N * Lf(N-2,ALPHA,X)
  !
  !  Restrictions:
  !
  !    -1 < ALPHA
  !
  !  Special values:
  !
  !    Lf(N,0,X) = L(N,X).
  !    Lf(N,ALPHA,X) = LM(N,ALPHA,X) for ALPHA integral.
  !
  !  Norm:
  !
  !    Integral ( 0 <= X < +oo ) exp ( - X ) * Lf(N,ALPHA,X)^2 dX
  !    = Gamma ( N + ALPHA + 1 ) / N!
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    10 March 2012
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Milton Abramowitz, Irene Stegun,
  !    Handbook of Mathematical Functions,
  !    National Bureau of Standards, 1964,
  !    ISBN: 0-486-61272-4,
  !    LC: QA47.A34.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) M, the number of evaluation points.
  !
  !    Input, integer ( kind = 4 ) N, the highest order function to compute.
  !
  !    Input, real ( kind = 8 ) ALPHA, the parameter.  -1 < ALPHA is required.
  !
  !    Input, real ( kind = 8 ) X(M), the evaluation points.
  !
  !    Output, real ( kind = 8 ) CX(1:M,0:N), the functions of 
  !    degrees 0 through N evaluated at the points X.
  !
    implicit none
    !
    integer ( kind = ik ) :: m
    integer ( kind = ik ) :: n
    !
    real ( kind = ark )  :: alpha
    real ( kind = ark )  :: cx(1:m,0:n)
    integer ( kind = 4 ) :: i
    real ( kind = ark )  :: x(1:m)
    !
    if ( alpha <= -1.0+small_ ) then
      write (out, '(a)' ) ' '
      write (out, '(a)' ) 'LF_FUNCTION - Fatal error!'
      write (out, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
      write (out, '(a)' ) '  but ALPHA must be greater than -1.'
      stop
    end if
    !
    if ( n < 0 ) then
      return
    end if
    !
    cx(1:m,0) = 1.0_ark
    !
    if ( n == 0 ) then
      return
    end if
    !
    cx(1:m,1) = 1.0_ark + alpha - x(1:m)
    !
    do i = 2, n
      !
      cx(1:m,i) = ( &
        ( real ( 2 * i - 1, kind = ark ) + alpha - x(1:m) ) * cx(1:m,i-1)   &
      + ( real (   - i + 1, kind = ark ) - alpha          ) * cx(1:m,i-2) ) &
        / real (     i,     kind = ark )
      !
    end do
   !
  return
  !
  end subroutine lf_function
  !

  subroutine lm_polynomial( mm, n, m, x, cx )
  !
  !*****************************************************************************80
  !
  !! LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
  !
  !  First terms:
  !
  !    M = 0
  !
  !    Lm(0,0,X) =   1
  !    Lm(1,0,X) =  -X   +  1
  !    Lm(2,0,X) =   X^2 -  4 X   +  2
  !    Lm(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
  !    Lm(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
  !    Lm(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
  !    Lm(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
  !
  !    M = 1
  !
  !    Lm(0,1,X) =    0
  !    Lm(1,1,X) =   -1,
  !    Lm(2,1,X) =    2 X - 4,
  !    Lm(3,1,X) =   -3 X^2 + 18 X - 18,
  !    Lm(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
  !
  !    M = 2
  !
  !    Lm(0,2,X) =    0
  !    Lm(1,2,X) =    0,
  !    Lm(2,2,X) =    2,
  !    Lm(3,2,X) =   -6 X + 18,
  !    Lm(4,2,X) =   12 X^2 - 96 X + 144
  !
  !    M = 3
  !
  !    Lm(0,3,X) =    0
  !    Lm(1,3,X) =    0,
  !    Lm(2,3,X) =    0,
  !    Lm(3,3,X) =   -6,
  !    Lm(4,3,X) =   24 X - 96
  !
  !    M = 4
  !
  !    Lm(0,4,X) =    0
  !    Lm(1,4,X) =    0
  !    Lm(2,4,X) =    0
  !    Lm(3,4,X) =    0
  !    Lm(4,4,X) =   24
  !
  !  Recursion:
  !
  !    Lm(0,M,X)   = 1 
  !    Lm(1,M,X)   = (M+1-X)
  !
  !    if 2 <= N:
  !
  !      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X) 
  !                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
  !
  !  Special values:
  !
  !    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal 
  !    to the Laguerre polynomials L(N,X).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    08 February 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Milton Abramowitz, Irene Stegun,
  !    Handbook of Mathematical Functions,
  !    National Bureau of Standards, 1964,
  !    ISBN: 0-486-61272-4,
  !    LC: QA47.A34.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) MM, the number of evaluation points.
  !
  !    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
  !    Note that polynomials 0 through N will be computed.
  !
  !    Input, integer ( kind = 4 ) M, the parameter.  M must be nonnegative.
  !
  !    Input, real ( kind = 8 ) X(MM), the evaluation points.
  !
  !    Output, real ( kind = 8 ) CX(MM,0:N), the associated Laguerre polynomials 
  !    of degrees 0 through N evaluated at the evaluation points.
  !
    implicit none
    !
    integer(ik),intent(in) :: mm
    integer(ik),intent(in) :: n
    integer(ik),intent(in) :: m
    real(ark),intent(in)   :: x(mm)
    real(ark),intent(out)  :: cx(mm,0:n)
    !
    integer (ik) i
    !
    if ( m < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LM_POLYNOMIAL - Fatal error!'
      write ( *, '(a,i8)' ) '  Input value of M = ', m
      write ( *, '(a)' ) '  but M must be nonnegative.'
      stop
    end if
    !
    if ( n < 0 ) then
      return
    end if
    !
    cx(:,0) = 1.0_ark
    !
    if ( n == 0 ) then
      return
    end if
    !
    cx(:,1) = real ( m + 1,ark) - x(:)
    !
    do i = 2, n
      cx(:,i) = &
        ( ( real(   m + 2 * i - 1,ark ) - x(:) ) * cx(:,i-1)   &
          + real( - m     - i + 1,ark )          * cx(:,i-2) ) &
          / real(           i,    ark )
    end do
    !
  end subroutine lm_polynomial
  !


  subroutine me_laguerre(nn,n,maxorder,xrange,drho,icoord,isingular,poten,mu_rr,verbose,g_numerov,energy)
    !
    implicit none
    !
    integer(ik),intent(in) :: nn
    integer(ik),intent(in) :: n
    integer(ik) :: m
    integer(ik),intent(in)   :: icoord,maxorder
    real(ark),intent(in) :: xrange(2)
    !
    integer(ik),intent(in) :: isingular
    real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:n,0:n)
    real(ark),intent(out)    :: energy(0:n)
    !
    real(ark),intent(in) :: poten(nn),mu_rr(nn),drho(nn,3)
    integer(ik),intent(in) :: verbose   ! Verbosity level
    !
    real(ark)  :: x(nn)
    real(ark)  :: f(nn,0:n)
    real(ark)  :: df(nn,0:n)
    real(ark)  :: cx(nn,0:n)
    real(ark)  :: f_m,x_(1:nn),factor,fnorm,alpha,f_by_sqrtx(1:nn),f_by_sqrtx_1(1:nn),dx,potmin,L,rho,psipsi_t
    integer(ik) :: i,k,imin
    !
    integer(ik) :: alloc,rec_len,io_slot
    character(len=cl)    :: unitfname 
    !
    integer(ik) :: vl,vr,lambda
    !
    real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:)

    allocate(phil(nn),phir(nn),dphil(nn),dphir(nn), &
             phivphi(nn),rho_kinet(nn),rho_poten(nn),rho_extF(nn),stat=alloc)
    if (alloc/=0) then 
      write (out,"('phi - out of memory')")
      stop 'phi - out of memory'
    endif 
    !
    ! step size 
    dx = (xrange(2)-xrange(1))/real(nn-1,kind=ark)
    !
    ! Do some reporting
    !
    if (verbose>=3) then 
        write (out,"('vmax = ',i8)") n
        write (out,"('maxorder = ',i8)") maxorder
        write (out,"('icoord = ',i4)") icoord
        write (out,"('xrange (x) = ',2f12.4)") xrange(1:2) !*180.0_ark/pi
        write (out,"('dx (x) = ',2f12.4)") dx 
    endif 
    !
    potmin = huge(1.0_ark)
    !
    do i=1,nn
       !
       if (poten(i)<potmin) then 
          imin = i
          potmin = poten(i)
       endif
       !
    enddo
    !
    if (imin<1.or.imin>nn) then 
        write(out,"('ML_box: pot_eff has no minimum',i8)") 
        stop 'ML_box: pot_eff has no minimum'
    endif 
    !
    ! define the rho-type coordinate 
    !
    rho_kinet(:) = drho(:,1)
    rho_poten(:) = drho(:,2)
    rho_extF(:)  = drho(:,3)

    !
    inquire(iolength=rec_len) f(:,0),df(:,0)
    !
    write(unitfname,"('Laguerre basis set # ',i6)") icoord
    call IOStart(trim(unitfname),io_slot)
    !
    open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
    !
    if ( n < 0 ) then
      return
    end if
    !
    f_m = 1.0_ark
    !
    m = 0
    !
    x_(:) = sqrt(x(:))/sqrt(f_m)
    !
    alpha = real(2*m+1,ark)/2.0_ark
    !
    call lm_polynomial(nn,n, m, x_,cx)
    !
    do i = 0, n
      !
      factor = exp(faclog(i+m))
      !
      fnorm = sqrt(2.0_ark*sqrt(f_m)/factor)
      !
      f(:,i) = sqrt(x_(:)**(alpha))*exp(-0.5_ark*x_(:))*cx(:,i)/fnorm
      !
      f_by_sqrtx(:) = sqrt(f_m**alpha)*sqrt(x(:)**(m))*exp(-0.5_ark*x_(:))*cx(:,i)/fnorm
      !
      f_by_sqrtx_1 = 0
      if (i>0) then
        f_by_sqrtx_1(:) = sqrt(f_m**alpha)*sqrt(x(:)**m)*exp(-0.5_ark*x_(:))*cx(:,i-1)*fnorm
      endif
      !
      ! diff(f,rho)*sqrt(rho)
      df(:,i) = ( (alpha+real(2*n,ark))-x_(:) )*f_by_sqrtx(:) -f_by_sqrtx_1(:)*real(2*(i+m),ark)
      !
      write (io_slot,rec=i+1) (f(k,i),k=1,nn),(df(k,i),k=1,nn)
      !
      phivphi(:) = f(:,i)**2
      factor = simpsonintegral_ark(nn-1,xrange(2)-xrange(1),phivphi)
      !
      continue
      !
    end do
    !
    do vl = 0,n
      energy(vl) = (vl+0.5_ark)
    enddo
    g_numerov = 0
    !
    !characvalue = maxval(enerslot(0:vmax))
    !
    do vl = 0,n
       !
       phil(:)  = f(:,vl)
       dphil(:) = df(:,vl)
       !
       do vr = vl,n
           !
           phir(:)  = f(:,vr)
           dphir(:) = df(:,vr)
           !
           ! Here we prepare integrals of the potential 
           ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
           ! obtained above by the Numerov
           !
           !phivphi(:) = phil(:)*poten(:)*phir(:)
           !
           psipsi_t = 0 
           phivphi(:) = phil(:)*phir(:)
           !
           psipsi_t = simpsonintegral_ark(nn-1,xrange(2)-xrange(1),phivphi)    
           !
           continue  
           !
       enddo
    enddo
    !
    deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,rho_poten,rho_extF)


    !
  end subroutine me_laguerre
  !



  !
  ! modified version of r8_gamma
  ! 
  function ark_gamma(x) result(f)
  !  
  !*****************************************************************************80
  !
  !! R8_GAMMA evaluates Gamma(X) for a real argument.
  !
  !  Discussion:
  !
  !    This routine calculates the gamma function for a real argument X.
  !
  !    Computation is based on an algorithm outlined in reference 1.
  !    The program uses rational functions that approximate the gamma
  !    function to at least 20 significant decimal digits.  Coefficients
  !    for the approximation over the interval (1,2) are unpublished.
  !    Those for the approximation for 12 <= X are from reference 2.
  !
  !  Modified:
  !
  !    11 February 2008
  !
  !  Author:
  !
  !    Original FORTRAN77 version by William Cody, Laura Stoltz.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    William Cody,
  !    An Overview of Software Development for Special Functions,
  !    in Numerical Analysis Dundee, 1975,
  !    edited by GA Watson,
  !    Lecture Notes in Mathematics 506,
  !    Springer, 1976.
  !
  !    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
  !    Charles Mesztenyi, John Rice, Henry Thatcher,
  !    Christoph Witzgall,
  !    Computer Approximations,
  !    Wiley, 1968,
  !    LC: QA297.C64.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) X, the argument of the function.
  !
  !    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
  !
    implicit none
  !
  !  Coefficients for minimax approximation over (12, INF).
  !
    real (ark),dimension ( 7 ) :: c = (/ &
     -1.910444077728e-03_ark, &
      8.4171387781295e-04_ark, &
     -5.952379913043012e-04_ark, &
      7.93650793500350248e-04_ark, &
     -2.777777777777681622553e-03_ark, &
      8.333333333333333331554247e-02_ark, &
      5.7083835261e-03_ark /)
      !
    real (ark), parameter :: eps = 2.22e-16_ark
    real (ark) fact
    integer (ik) :: i
    integer (ik) :: n
    real (ark), dimension ( 8 ) :: p = (/ &
      -1.71618513886549492533811e+00_ark, &
       2.47656508055759199108314e+01_ark, &
      -3.79804256470945635097577e+02_ark, &
       6.29331155312818442661052e+02_ark, &
       8.66966202790413211295064e+02_ark, &
      -3.14512729688483675254357e+04_ark, &
      -3.61444134186911729807069e+04_ark, &
       6.64561438202405440627855e+04_ark /)
    !
    logical parity
    real (ark), parameter :: pi = 3.1415926535897932384626434_ark
    real (ark), dimension ( 8 ) :: q = (/ &
      -3.08402300119738975254353e+01_ark, &
       3.15350626979604161529144e+02_ark, &
      -1.01515636749021914166146e+03_ark, &
      -3.10777167157231109440444e+03_ark, &
       2.25381184209801510330112e+04_ark, &
       4.75584627752788110767815e+03_ark, &
      -1.34659959864969306392456e+05_ark, &
      -1.15132259675553483497211e+05_ark /)
      !
    real (ark) :: f
    real (ark) ::  res
    real (ark),parameter :: sqrtpi = 0.9189385332046727417803297_ark
    real (ark) :: sum
    real (ark),intent(in) :: x
    real (ark), parameter :: xbig = 171.624D+00
    real (ark) :: xden
    real (ark), parameter :: xinf = 1.0D+30
    real (ark), parameter :: xminin = 2.23D-308
    real (ark) :: xnum
    real (ark) :: y
    real (ark) :: y1
    real (ark) :: ysq
    real (ark) :: z
    !  
    parity = .false.
    fact = 1.0_ark
    n = 0
    y = x
    !
    !  Argument is negative.
    !
    if ( y <= small_ ) then
      !
      y = - x
      y1 = aint(y)
      res = y - y1
      !
      if ( res /= 0.0_ark ) then
        !
        if ( y1 /= aint ( y1 * 0.5_ark ) * 2.0_ark ) then
          parity = .true.
        end if
        !
        fact = - pi / sin ( pi * res )
        y = y + 1.0_ark
        !
      else
        !
        res = xinf
        f = res
        return
        !
      end if
        !
    end if
    !
    !  Argument is positive.
    !
    if ( y < eps ) then
      !
      !  Argument < EPS.
      !
      if ( xminin <= y ) then
        res = 1.0_ark / y
      else
        res = xinf
        f = res
        return
      end if
      !
    else if ( y < 12.0_ark ) then
      !
      y1 = y
      !
    !  0.0 < argument < 1.0.
    if ( y < 1.0_ark ) then
      !
      z = y
      y = y + 1.0_ark
      !
    !  1.0 < argument < 12.0.
    !  Reduce argument if necessary.
    !
    else
      !
      n = int ( y ) - 1
      y = y - real ( n,ark)
      z = y - 1.0_ark
      !
    end if
    !
    !  Evaluate approximation for 1.0 < argument < 2.0.
    !
    xnum = 0
    xden = 1.0_ark
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do
    !
    res = xnum / xden + 1.0_ark
    !
    !  Adjust result for case  0.0 < argument < 1.0.
    !
    if ( y1 < y ) then
      !
      res = res / y1
    !
    !  Adjust result for case 2.0 < argument < 12.0.
    !
    else if ( y < y1 ) then
      !
      do i = 1, n
        res = res * y
        y = y + 1.0_ark
      end do
      !
    end if
    !
  else
    !
    !  Evaluate for 12.0 <= argument.
    !
    if ( y <= xbig ) then
      !
      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5_ark ) * log(y)
      res = exp ( sum )
      !
    else
      !
      res = xinf
      f = res
      return
      !
    end if
    !
  end if
  !
  !  Final adjustments and return.
  !
  if ( parity ) then
    res = - res
  end if
  !
  if ( fact /= 1.0_ark ) then
    res = fact / res
  end if
  !
  f = res
  !
  return
  !
  end function ark_gamma


end module me_str

module me_bnd
  use accuracy
  use me_numer
  use molecules
  use timer
  use lapack
  !
  implicit none

  public ik, rk, out
  public degener_harm_q,ME_box,ME_Fourier,ME_Legendre,ME_Associate_Legendre,ME_sinrho_polynomial

  integer(ik), parameter :: verbose     = 1                       ! Verbosity level
  integer(ik) :: Nr = 4                          ! 2*Nr+1 is the number of interpolation points

  private

!
  contains
!

!********************************************************************!
!                                                                    
!     Bending matrix elements                                    
!     2D Isotropic double degenerated harmonic oscilator
!                                                                    
!********************************************************************!
!                                                                    
!     Bending  matrix elements for q4^k4*q5^k5  term                 
!     on the double degenerated harmonic oscilator                   
!                                                                    
  subroutine degener_harm_q(maxorder,vbndmaxtmp,me) 


    integer(ik),intent(in) :: maxorder,vbndmaxtmp
    real(rk),intent(out)   :: me(0:maxorder,0:maxorder,-maxorder:maxorder,-maxorder:maxorder,0:vbndmaxtmp,0:vbndmaxtmp,4)

    integer(ik) :: n1,n2,v,l,dv,dl,n,k,kl,k1,k2
    real(rk) :: u,p

    ! 
    !              here we  go! 
    !

    if (verbose>1) write (out,"(20('*'),' calculation of the degenerated oscilator matrix elements.')")
    !
    me = 0
    do v = 0,vbndmaxtmp 
       do kl= 0, int(v/2)
          l = v-2*kl
          u = real(v,rk)
          p = real(l,rk)
          !
          me(0,0,0,0,v,l,1) = 1.0_rk
          if ( l/=0 ) me(0,0,0,0,v,l,4) = 1.0_rk
 
          if ( l==0 ) then 
             me(1,0,-1, 1,v,0,1) = -0.5_rk*sqrt( u        )
             me(1,0, 1, 1,v,0,1) =  0.5_rk*sqrt( u+2.0_rk )

             me(0,1,-1, 1,v,0,3) = -0.5_rk*sqrt( u        )
             me(0,1, 1, 1,v,0,3) =  0.5_rk*sqrt( u+2.0_rk )

          elseif ( l.eq.1 ) then 

             me(1,0,-1,-1,v,1,1) =  0.5_rk*sqrt( u+1.0_rk          )
             me(1,0,-1, 1,v,1,1) = -0.5_rk*sqrt( 0.5_rk*(u-1.0_rk  ) )
             me(1,0, 1,-1,v,1,1) = -0.5_rk*sqrt( u+1._rk           )
             me(1,0, 1, 1,v,1,1) =  0.5_rk*sqrt( 0.5_rk*(u+3.0_rk) )

             me(0,1,-1,-1,v,1,2) =  0.5_rk*sqrt( u+1.0_rk )
             me(0,1,-1, 1,v,1,2) =  0.5_rk*sqrt( 0.5_rk*(u-1.0_rk  ) )
             me(0,1, 1,-1,v,1,2) = -0.5_rk*sqrt( u+1.0_rk )
             me(0,1, 1, 1,v,1,2) = -0.5_rk*sqrt( 0.5_rk*(u+3.0_rk) )

             me(1,0,-1, 1,v,1,4) =  me(1,0,-1, 1,v,1,1)
             me(1,0, 1, 1,v,1,4) =  me(1,0, 1, 1,v,1,1)

             me(0,1,-1, 1,v,1,3) = -me(0,1,-1, 1,v,1,2)
             me(0,1, 1, 1,v,1,3) = -me(0,1, 1, 1,v,1,2)

          else

             me(1,0,-1,-1,v,l,1) =  0.5_rk*sqrt( 0.5_rk*(u+p     ) )
             me(1,0,-1, 1,v,l,1) = -0.5_rk*sqrt( 0.5_rk*(u-p     ) )
             me(1,0, 1,-1,v,l,1) = -0.5_rk*sqrt( 0.5_rk*(u-p+2.0_rk) )
             me(1,0, 1, 1,v,l,1) =  0.5_rk*sqrt( 0.5_rk*(u+p+2.0_rk) )

             me(0,1,-1,-1,v,l,2) =  0.5_rk*sqrt( 0.5_rk*(u+p     ) )
             me(0,1,-1, 1,v,l,2) =  0.5_rk*sqrt( 0.5_rk*(u-p     ) )
             me(0,1, 1,-1,v,l,2) = -0.5_rk*sqrt( 0.5_rk*(u-p+2.0_rk) )
             me(0,1, 1, 1,v,l,2) = -0.5_rk*sqrt( 0.5_rk*(u+p+2.0_rk) )

             me(1,0,-1,-1,v,l,4) =  me(1,0,-1,-1,v,l,1)
             me(1,0,-1, 1,v,l,4) =  me(1,0,-1, 1,v,l,1)
             me(1,0, 1,-1,v,l,4) =  me(1,0, 1,-1,v,l,1)
             me(1,0, 1, 1,v,l,4) =  me(1,0, 1, 1,v,l,1)

             me(0,1,-1,-1,v,l,3) = -me(0,1,-1,-1,v,l,2)
             me(0,1,-1, 1,v,l,3) = -me(0,1,-1, 1,v,l,2)
             me(0,1, 1,-1,v,l,3) = -me(0,1, 1,-1,v,l,2)
             me(0,1, 1, 1,v,l,3) = -me(0,1, 1, 1,v,l,2)

          endif

       enddo  ! --- kl
    enddo  ! --- v


    call me_tobe_zero_for_nonphys_region(maxorder,vbndmaxtmp,me)

    do n = 2,maxorder
       do n1= 0,n
          n2 = n - n1
          do dv = -n,n
             do dl = -n,n
                do v = 0,vbndmaxtmp-1    !-n 
                   do l = 0, v
                      do k1 = 1, 2
                         do k2 = 1, 2
                            k = k1*k2+dim(k1,k2)
                           if( (v+dv).ge.0  .and. (l+dl).ge.0 ) then
                             me(n1,n2,dv,dl,v,l,k) = matrrec_q(n1,n2,dv,dl,v,l,k1,k2)
                           endif 

                         enddo 
                      enddo 
                   enddo 
                enddo 
             enddo 
          enddo 
       enddo 
    enddo 


  !
  contains  
  !
    function matrrec_q(n1,n2,dv,dl,v,l,k1,k2) result (matrel)

      integer(ik),intent(in) :: n1,n2,dv,dl,v,l,k1,k2
      integer(ik)            :: k0,l0
      real(rk)               :: matrel

      matrel = 0
      !
      if (n1 .eq. 0 ) then 
         do k0 = -1,1,2
            do l0 = -1,1,2
               if ( v+k0.ge.0 .and. l+l0.ge.0 .and. iabs(dv-k0).le.maxorder .and. iabs(dl-l0).le.maxorder ) then
                  if ( k2.eq.1 ) then 
                     matrel=   matrel + me(0,n2-1,dv-k0,dl-l0,v+k0,l+l0,2*k1)*me(0,1,k0,l0,v,l,3)
                  else 
                     matrel=   matrel + me(0,n2-1,dv-k0,dl-l0,v+k0,l+l0,k1+dim(k1,1))*me(0,1,k0,l0,v,l,2)
                  endif 
               endif
            enddo 
         enddo
      else 
         do k0 = -1,1,2
            do l0 = -1,1,2
               if ( v+k0.ge.0 .and. l+l0.ge.0 .and. iabs(dv-k0).le.maxorder .and. iabs(dl-l0).le.maxorder ) then
                  if ( k2.eq.1 ) then 
                     matrel=   matrel + me(n1-1,n2,dv-k0,dl-l0,v+k0,l+l0,k1+dim(k1,1))*me(1,0,k0,l0,v,l,1 )
                  else 
                     matrel=   matrel + me(n1-1,n2,dv-k0,dl-l0,v+k0,l+l0,2*k1)*me(1,0,k0,l0,v,l,4 )
                  endif
               endif
            enddo 
         enddo 
      endif 

    end function matrrec_q

  end subroutine degener_harm_q
  !
  !
  !
  !   Bending  matrix elements for p4^k4*p5^k5  term
  !   on the double degenerated harmonic oscilator  
  !
  subroutine degener_harm_p( vbndmaxtmp,me ) 

    integer(ik),intent(in) :: vbndmaxtmp
    real(rk),intent(out) :: me(0:1,0:1,-1:1,-1:1,0:vbndmaxtmp,0:vbndmaxtmp,4)
    integer(ik)          :: v,l,kl
    real(rk)             :: u,p

    me = 0
    do v = 0,vbndmaxtmp 
       do kl= 0, int(v/2)
          l = v-2*kl 
          u = real(v,rk)
          p = real(l,rk)
          !
          me(0,0,0,0,v,l,1) = 1.0_rk
          if ( l.ne.0 ) me(0,0,0,0,v,l,4) = 1.0_rk

          if ( l.eq.0 ) then 

             me(1,0,-1, 1,v,0,1) = -0.5_rk*sqrt( u )
             me(1,0, 1, 1,v,0,1) = -0.5_rk*sqrt( u+2.0_rk )

             me(0,1,-1, 1,v,0,3) = -0.5_rk*sqrt( u )
             me(0,1, 1, 1,v,0,3) = -0.5_rk*sqrt( u+2.0_rk )

          elseif ( l.eq.1 ) then 

             me(1,0,-1,-1,v,1,1) =  0.5_rk*sqrt( u+1.0_rk )
             me(1,0,-1, 1,v,1,1) = -0.5_rk*sqrt( 0.5_rk*(u-1.0_rk     ) )
             me(1,0, 1,-1,v,1,1) =  0.5_rk*sqrt( u+1.0_rk )
             me(1,0, 1, 1,v,1,1) = -0.5_rk*sqrt( 0.5_rk*(u+3.0_rk) )

             me(0,1,-1,-1,v,1,2) =  0.5_rk*sqrt( u+1.0_rk )
             me(0,1,-1, 1,v,1,2) =  0.5_rk*sqrt( 0.5_rk*(u-1.0_rk  ) )
             me(0,1, 1,-1,v,1,2) =  0.5_rk*sqrt( u+1.0_rk )
             me(0,1, 1, 1,v,1,2) =  0.5_rk*sqrt( 0.5_rk*(u+3.0_rk) )

             me(0,1,-1, 1,v,1,3) = -me(0,1,-1, 1,v,1,2)
             me(0,1, 1, 1,v,1,3) = -me(0,1, 1, 1,v,1,2)

             me(1,0,-1, 1,v,1,4) =  me(1,0,-1, 1,v,1,1)
             me(1,0, 1, 1,v,1,4) =  me(1,0, 1, 1,v,1,1)


          else

             me(1,0,-1,-1,v,l,1) =  0.5_rk*sqrt( 0.5_rk*(u+p     ) )
             me(1,0,-1, 1,v,l,1) = -0.5_rk*sqrt( 0.5_rk*(u-p     ) )
             me(1,0, 1,-1,v,l,1) =  0.5_rk*sqrt( 0.5_rk*(u-p+2.0_rk) )
             me(1,0, 1, 1,v,l,1) = -0.5_rk*sqrt( 0.5_rk*(u+p+2.0_rk) )

             me(0,1,-1,-1,v,l,2) =  0.5_rk*sqrt( 0.5_rk*(u+p     ) )
             me(0,1,-1, 1,v,l,2) =  0.5_rk*sqrt( 0.5_rk*(u-p     ) )
             me(0,1, 1,-1,v,l,2) =  0.5_rk*sqrt( 0.5_rk*(u-p+2.0_rk) )
             me(0,1, 1, 1,v,l,2) =  0.5_rk*sqrt( 0.5_rk*(u+p+2.0_rk) )

             me(1,0,-1,-1,v,l,4) =  me(1,0,-1,-1,v,l,1)
             me(1,0,-1, 1,v,l,4) =  me(1,0,-1, 1,v,l,1)
             me(1,0, 1,-1,v,l,4) =  me(1,0, 1,-1,v,l,1)
             me(1,0, 1, 1,v,l,4) =  me(1,0, 1, 1,v,l,1)

             me(0,1,-1,-1,v,l,3) = -me(0,1,-1,-1,v,l,2)
             me(0,1,-1, 1,v,l,3) = -me(0,1,-1, 1,v,l,2)
             me(0,1, 1,-1,v,l,3) = -me(0,1, 1,-1,v,l,2)
             me(0,1, 1, 1,v,l,3) = -me(0,1, 1, 1,v,l,2)

          endif 

      enddo  ! --- k
    enddo    ! --- v


    call me_tobe_zero_for_nonphys_region(1_ik,vbndmaxtmp,me)


  end subroutine degener_harm_p

!
!     the bending  matrix elements of pa* q4^k4*q5^k5* pb term  
!     on the double degenerated harmonic oscilator 
!
  subroutine degener_harm_g( maxorder,vbndmax,lbndmax,vbndmaxtmp,me_q, me_p,coeff_norm,g_ben ) 

    integer(ik),intent(in) :: maxorder,vbndmaxtmp
    integer(ik),intent(in) :: vbndmax,lbndmax
    real(rk),intent(in)    :: me_q(0:maxorder,0:maxorder,-maxorder:maxorder,-maxorder:maxorder,0:vbndmaxtmp,0:vbndmaxtmp,4)
    real(rk),intent(in)    :: me_p(0:1,0:1,-1:1,-1:1,0:vbndmaxtmp,0:vbndmaxtmp,4)
    real(rk),intent(in)    :: coeff_norm
    real(rk),intent(out)    :: g_ben(9,int((maxorder+2)*(maxorder+1)/2),0:vbndmax,0:lbndmax,0:vbndmax,0:lbndmax,4)

    integer(ik)            :: n1,n2,dv,dl,n,k,k1,k2,pl,pr
    real(rk)               :: g_tmp,minus,normfact
    integer(ik)            :: p,pmax
    integer(ik)            :: a(2,2,2),deltap(2,2)
    integer(ik)            :: vl,vr,ll,lr,kl,kr,n0


    deltap = 0
    deltap(1,1) = 1
    deltap(2,2) = 1

    a(1,1,1) = 2
    a(2,1,1) = 0

    a(1,1,2) = 1
    a(2,1,2) = 1

    a(1,2,1) = 1
    a(2,2,1) = 1

    a(1,2,2) = 0
    a(2,2,2) = 2


    ! concor=planck*avogno*1.0d+16/(4.0d+00*pi*pi*vellgt)

    ! cosr = cos(rhoe)
    ! g044 = 3.0_rk/2.0_rk*(3.0_rk*m1-6.0_rk*m1*cosr**2+3.0_rk*m1*cosr**4+2.0_rk*m4+2.0_rk*m4*cosr**2)/m4/m1/re14**2/(1.0_rk+3.0_rk*cosr**2)*concor
    ! coeff_norm = sqrt(sqrt(g044*0.5_rk/f44))

    g_ben = 0

    do n = 0,maxorder
       do n1= 0,n
          n2 = n - n1
          n0= int( n*(n+1)/2 )+n1+1
          do pr = 0, 2
             do pl = 0, 2
                if ( pl+pr.eq.0 ) then 
                   pmax = 0 
                elseif ( pl*pr.eq.0 ) then   
                   pmax = 1 
                else 
                   pmax = 2 
                endif 
                p = pr*3+pl+1
                normfact = 1.0_rk
                if ( n1+n2-pmax.ne.0 ) normfact = coeff_norm**(n1+n2-pmax)
                do vl = 0,vbndmax
                   do vr = 0,vbndmax
                      dv = vl-vr
                      do kl = 0, int(vl/2)
                         do kr = 0, int(vr/2)
                            ll = vl-2*kl
                            lr = vr-2*kr
                            dl = ll-lr
                            do k1 = 1, 2
                               do k2 = 1, 2
                                  if ( iabs(dv).le.n+pmax .and. iabs(dl).le.n+pmax ) then
                                     k = k1*k2+dim(k1,k2)
                                     !
                                     select case  (pmax)
                                     case(0)
                                       g_tmp=me_q(n1,n2,dv,dl,vr,lr,k)
                                     case(1)
                                       g_tmp= matrelem_g0(pl,pr,n1,n2,dv,dl,vr,lr,k1,k2)
                                     case(2)
                                       if ( pl.eq.pr ) then 
                                           g_tmp= matrelem_g(pl,pl,n1,n2,dv,dl,vr,lr,k1,k2)
                                       else 
                                           g_tmp= matrelem_g(pl,pr,n1,n2,dv,dl,vr,lr,k1,k2) &
                                                 +matrelem_g(pr,pl,n1,n2,dv,dl,vr,lr,k1,k2)
                                       endif

                                     end select
                                     !
                                     minus = 1.0_rk
                                     if ( k1.eq.2 .and. mod(ll+3,3).eq.2 ) minus = -minus
                                     if ( k2.eq.2 .and. mod(lr+3,3).eq.2 ) minus = -minus
                                     !
                                     g_ben(p,n0,vl,kl,vr,kr,k) = g_tmp*minus*normfact
                                  endif 
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo




  !
  contains  
  !
    function matrelem_g0(p1,p2,n1,n2,dv,dl,v,l,k1,k2) result (matrel)

      integer(ik),intent(in) :: p1,p2,n1,n2,dv,dl,v,l,k1,k2
      integer(ik)            :: v0,l0,v1,l1,ind1,ind2,ind0,il1,il2,ir1,ir2
      real(rk)               ::  matrel
      
      matrel = 0
      if ( p2.eq.0 ) then 
         il1 = deltap(1,p1)
         il2 = deltap(2,p1)
         do v0 = -1,1,2
            do l0 = -1,1,2
               if ( v+dv+v0.ge.0 .and. l+dl+l0.ge.0 .and.l+dl+l0.le.v+dv+v0 .and. &
                    iabs(dv+v0).le.maxorder .and. iabs(dl+l0).le.maxorder ) then
                  ind1 = mod(k1+p1,2)+1
                  ind0 = ind1*k2   + dim(ind1,k2)
                  ind1 =   k1*ind1 + dim(  k1,ind1)
                  matrel=   matrel + me_p(il1,il2,-v0  ,-l0  ,v+dv+v0,l+dl+l0,ind1)*  &
                                     me_q(n1 ,n2 ,dv+v0,dl+l0,v      ,l      ,ind0)

               endif 
            enddo 
         enddo 
      else

         ir1= deltap(1,p2)
         ir2= deltap(2,p2)

         do v0 = -1,1,2
            do l0 = -1,1,2
               if ( v+v0.ge.0.and.l+l0.ge.0.and.l+l0.le.v+v0.and.   &
                   iabs(dv-v0).le.maxorder .and. iabs(dl-l0).le.maxorder ) then
                 ind1 = mod(k2+p2,2)+1
                 ind0 = ind1*k2   + dim(ind1,k2)
                 ind1 =   k1*ind1 + dim(  k1,ind1)
                 matrel=   matrel +  me_q(n1,n2,dv-v0,dl-l0,v+v0,l+l0,ind1)*me_p(ir1,ir2,v0,l0,v,l,ind0)
               endif
            enddo 
         enddo 
      endif

    end function matrelem_g0
    !
    !
    function matrelem_g(p1,p2,n1,n2,dv,dl,v,l,k1,k2)  result (matrel)

      integer(ik),intent(in) :: p1,p2,n1,n2,dv,dl,v,l,k1,k2
      integer(ik)            :: v0,l0,v1,l1,ind1,ind2,ind0,il1,il2,ir1,ir2
      real(rk)               ::  matrel

      matrel = 0

      il1 = deltap(1,p1)
      il2 = deltap(2,p1)
      ir1 = deltap(1,p2)
      ir2 = deltap(2,p2)
      do v0 = -1,1,2
         do l0 = -1,1,2
            if ( v+dv+v0.ge.0 .and. l+dl+l0.ge.0 .and. l+dl+l0.le.v+dv+v0 ) then
               do v1 = -1,1,2
                  do l1 = -1,1,2
                     if ( v+v1.ge.0 .and. l+l1.ge.0  .and. l+l1.le.v+v1.and. &
                         iabs(dv+v0-v1).le.maxorder .and. iabs(dl+l0-l1).le.maxorder ) then
                         !
                         ind1 = mod(k1+p1,2)+1
                         ind2 = mod(k2+p2,2)+1
                         !
                         ind0 = ind1*ind2 + dim(ind1,ind2)
                         ind1 =   k1*ind1 + dim(  k1,ind1)
                         ind2 = ind2*k2   + dim(ind2,k2  )
                         !
                         matrel=   matrel + me_p(il1,il2,-v0,-l0,v+dv+v0,l+dl+l0,ind1)*   &
                                   me_q(n1 ,n2 ,dv+v0-v1,dl+l0-l1,v+v1   ,l+l1   ,ind0)*  &
                                   me_p(ir1,ir2,v1      ,l1      ,v      ,l      ,ind2)


                     endif
                  enddo 
               enddo 
            endif
         enddo 
      enddo 
    !
    end function matrelem_g
  !
  end subroutine degener_harm_g
  !
  !
  subroutine me_tobe_zero_for_nonphys_region(nmax,vbndmaxtmp,me) 

    integer(ik),intent(in) :: nmax 
    integer(ik),intent(in) :: vbndmaxtmp 
    real(rk),intent(inout) :: me(0:nmax,0:nmax,-nmax:nmax,-nmax:nmax,0:vbndmaxtmp,0:vbndmaxtmp,4)
    integer(ik) :: n1,n2,dv,dl,v,k,l

    do n1 = 0,nmax
       n2 = nmax - n1
       do dv = -1,1
          do dl = -1,1
             do v = 0,vbndmaxtmp
                do l = 0, v
                   do k = 1, 4
                      if ( (v+dv).lt.0 .or. (l+dl).lt.0 .or. l+dl.gt.v+dv)  me(n1,n2,dv,dl,v,l,k) = 0.0_rk
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
 
  end subroutine me_tobe_zero_for_nonphys_region 
  !



  !
  ! Matrix elements with box-eigenfunctions 
  !
  subroutine ME_box(vmax,maxorder,rho_b_,isingular,npoints,drho,poten,mu_rr,icoord,periodic,verbose,g_numerov,energy)
   !
   integer(ik),intent(in) :: vmax,maxorder,npoints,isingular
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:vmax,0:vmax)
   real(ark),intent(out)    :: energy(0:vmax)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten(0:npoints),mu_rr(0:npoints),drho(0:npoints,3)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose   ! Verbosity level
   logical,intent(in)     :: periodic

   real(ark)            :: rho,L,rhostep,potmin,rhostep_
   real(ark)            :: psipsi_t,characvalue,rho_b(2)
   !
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,n,imin,io_slot,kl,kr,p
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:)
   !real(ark),allocatable :: f(:),poten(:),mu_rr(:)

   character(len=cl)    :: unitfname 
    !
    if (verbose>=1) write (out,"(/20('*'),' Particle in a box matrix elements calculations')")
     !
     ! global variables 
     !
     allocate(phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints), &
              phivphi(0:npoints),rho_kinet(0:npoints),rho_poten(0:npoints),rho_extF(0:npoints),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     rho_b = rho_b_
     !
     ! step size 
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
     !
     ! periodic factor to produce cos( k nx ) and sin(k nx ) functions
     p = 1
     !if (periodic>0) then 
     !  p = periodic
     !endif
     !
     ! Do some reporting
     !
     if (verbose>=3) then 
         write (out,"('vmax = ',i8)") vmax
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
     inquire(iolength=rec_len) phil(:),dphil(:)
     !
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
     L = rho_b(2)-rho_b(1)
     !
     do vl = 0,vmax
       kl = vl*p
       energy(vl) = 0.5_ark*(kl+1)**2*mu_rr(imin)*pi**2/L**2
     enddo
     !
     !characvalue = maxval(enerslot(0:vmax))
     !
     do vl = 0,vmax
        kl = vl*p
        !
        do i=0,npoints
           !
           rho = real(i,kind=ark)*rhostep
           !
           phil(i)  = sqrt(2.0_ark/L)*sin(real(kl+1,ark)*pi*rho/L)
           dphil(i) = sqrt(2.0_ark/L)*cos(real(kl+1,ark)*pi*rho/L)*real(kl+1,ark)*pi/L
           !
        enddo
        !
        write (io_slot,rec=vl+1) (phil(i),i=0,npoints),(dphil(i),i=0,npoints)
        !
        do vr = vl,vmax
            kr = vr*p
            !
            if (vl==vr) then
                phir =  phil
               dphir = dphil
            else
              do i=0,npoints
                 !
                 rho = real(i,kind=ark)*rhostep
                 !
                 phir(i)  = sqrt(2.0_ark/L)*sin(real(kr+1,ark)*pi*rho/L)
                 dphir(i) = sqrt(2.0_ark/L)*cos(real(kr+1,ark)*pi*rho/L)*real(kr+1,ark)*pi/L
                 !
              enddo
            endif
            !
            ! Here we prepare integrals of the potential 
            ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
            ! obtained above by the Numerov
            !
            !phivphi(:) = phil(:)*poten(:)*phir(:)
            !
            !h_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
            !
            ! momenta-quadratic part 
            !
            !phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)
            !
            !psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
            !
            ! Add the diagonal kinetic part to the tested mat. elem-s
            !
            !h_t = h_t - 0.5_ark*psipsi_t
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
               g_numerov(0,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
               !
               ! external field expansion
               !
               if (lambda==0) then 
                  phivphi(:) = phil(:)*phir(:)
               else
                  phivphi(:) = phil(:)*rho_extF(:)**lambda*phir(:)
               endif
               !
               g_numerov(3,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
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
               g_numerov(-1,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
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
               g_numerov(2,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
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
               g_numerov(1,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
               !
               if (vl/=vr) then
                  !
                  if (lambda==0) then 
                     phivphi(:) = dphil(:)*phir(:)
                  else
                     phivphi(:) = dphil(:)*rho_kinet(:)**lambda*phir(:)
                  endif
                  !
                  g_numerov(1,lambda,vr,vl) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
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
        enddo
     enddo
     !
     deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,rho_poten,rho_extF)
     !
     !
  end subroutine ME_box


  !
  ! Matrix elements with Fourier-eigenfunctions 
  !
  subroutine ME_Fourier(vmax,maxorder,rho_b_,isingular,npoints_,numerpoints_,drho_,poten_,mu_rr_,icoord,iperiod,verbose,g_numerov,energy)
   !
   integer(ik),intent(in) :: vmax,maxorder,npoints_,numerpoints_,isingular
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:vmax,0:vmax)
   real(ark),intent(out)    :: energy(0:vmax)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten_(0:npoints_),mu_rr_(0:npoints_),drho_(0:npoints_,3)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose   ! Verbosity level
   integer(ik),intent(in) :: iperiod
   !
   integer(ik),parameter  :: Factor_FF=30 ! factor to increase the Fourier basis set size 
   !
   real(ark)            :: rho,L,rhostep,potmin,rhostep_
   real(ark)            :: psipsi_t,characvalue,rho_b(2),cross_prod,factor,fval,df_t,step_scale
   !
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,n,imin,io_slot,kl,kr,p,fmax,npoints,i_,i1,i2
   !
   real(ark),allocatable :: poten(:),mu_rr(:),rho_(:)
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:),phi(:),dphi(:),vect(:,:)
   real(ark),allocatable :: phil_(:),phir_(:),dphil_(:),dphir_(:),phivphi_(:)
   real(ark),allocatable :: psi(:,:),dpsi(:,:),psi_(:,:),dpsi_(:,:)
   real(rk),allocatable ::  h(:,:),ener(:)
   !
   character(len=cl)    :: unitfname 
    !
    if (verbose>=1) write (out,"(/20('*'),' Fourier real functions primitive matrix elements calculations')")
     !
     ! global variables 
     !
     ! large grid 
     npoints   = numerpoints_
     !
     allocate(phil_(0:npoints_),phir_(0:npoints_),dphil_(0:npoints_),dphir_(0:npoints_), &
              phivphi_(0:npoints_),rho_kinet(0:npoints_),rho_poten(0:npoints_),rho_extF(0:npoints_),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif
     !
     allocate(phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints), &
              phivphi(0:npoints),poten(0:npoints),mu_rr(0:npoints),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     fmax = vmax*factor_ff*iperiod
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
     allocate(h(fmax+1,fmax+1),ener(fmax+1),phi(fmax+1),dphi(fmax+1),stat=alloc)
     call ArrayStart('h-Fourier',alloc,size(h),kind(h))
     call ArrayStart('h-Fourier',alloc,size(ener),kind(ener))
     call ArrayStart('h-Fourier-phi',alloc,size(phi),kind(phi))
     call ArrayStart('h-Fourier-phi',alloc,size(dphi),kind(dphi))
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
         write(out,"('ML_box: pot_eff has no minimum',i8)") 
         stop 'ML_box: pot_eff has no minimum'
     endif 
     !
     ! define the rho-type coordinate 
     !
     rho_kinet(:) = drho_(:,1)
     rho_poten(:) = drho_(:,2)
     rho_extF(:)  = drho_(:,3)
     !
     inquire(iolength=rec_len) phil_(:),dphil_(:)
     !
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
     !characvalue = maxval(enerslot(0:vmax))
     !
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
     !
     call lapack_syev(h,ener)
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
     do vl =  1,vmax+1
       !
       cross_prod = sum(vect(:,vl)*vect(:,vl))
       !
       factor = 1.0_ark/sqrt(cross_prod)
       !
       vect(:,vl) = vect(:,vl)*factor
       !
       do vr = vl+1,vmax+1
         !
         cross_prod = sum(vect(:,vl)*vect(:,vr))
         !
         vect(:,vr) = vect(:,vr)-cross_prod*vect(:,vl)
         ! 
       enddo
       !
     enddo
     !
     write (out,"(/' Fourier-optimized energies are:')") 
     !
     do vl=0,vmax   
       write (out,"(i8,f18.8)") vl,energy(vl)
     enddo
     !
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
     !
     ! dump the eigenfunction
     !
     if (npoints_==npoints) then 
          Psi_  = Psi
          DPsi_ = DPsi
     else
       do i_ = 0,npoints_
          !
          i = nint( real(i_,ark)/step_scale )
          !
          Psi_(1:vmax+1,i_) = Psi(1:vmax+1,i)
          DPsi_(1:vmax+1,i_) = DPsi(1:vmax+1,i)
          !
       enddo
     endif
     !
     deallocate(vect)
     call ArrayStop('h-vect')
     !
     do vl = 0,vmax
        kl = (vl+1)/2*p
        !
        phil_(:)  =  Psi_(vl+1,:)
        dphil_(:) = dPsi_(vl+1,:)
        !
        write (io_slot,rec=vl+1) (phil_(i),i=0,npoints_),(dphil_(i),i=0,npoints_)
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
               g_numerov(0,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
               !
               ! external field expansion
               !
               if (lambda==0) then 
                  phivphi_(:) = phil_(:)*phir_(:)
               else
                  phivphi_(:) = phil_(:)*rho_extF(:)**lambda*phir_(:)
               endif
               !
               g_numerov(3,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
               if (vl/=vr) g_numerov(3,lambda,vr,vl) = g_numerov(3,lambda,vl,vr)
               !
               ! momenta-free in kinetic part 
               !
               if (lambda==0) then 
                  phivphi_(:) = phil_(:)*phir_(:)
               else
                  phivphi_(:) = phil_(:)*rho_kinet(:)**lambda*phir_(:)
               endif
               !
               g_numerov(-1,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
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
                  phivphi_(:) =-dphil_(:)*dphir_(:)
               else
                  phivphi_(:) =-dphil_(:)*rho_kinet(:)**lambda*dphir_(:)
               endif
               !
               g_numerov(2,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
               !
               if (vl/=vr) g_numerov(2,lambda,vr,vl) = g_numerov(2,lambda,vl,vr)
               !
               ! momenta-linear part:
               ! < vl | d/dx g(x) | vr > = - < vr | g(x) d/dx | vl >
               !
               !
               if (lambda==0) then 
                  phivphi_(:) = phil_(:)*dphir_(:)
               else
                  phivphi_(:) = phil_(:)*rho_kinet(:)**lambda*dphir_(:)
               endif
               !
               g_numerov(1,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
               !
               if (vl/=vr) then
                  !
                  if (lambda==0) then 
                     phivphi_(:) = dphil_(:)*phir_(:)
                  else
                     phivphi_(:) = dphil_(:)*rho_kinet(:)**lambda*phir_(:)
                  endif
                  !
                  g_numerov(1,lambda,vr,vl) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi_)
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
        enddo
        !
        if (verbose>=6) then 
           !
           write (out,"('v = ',i8,f18.8)") vl,h(vl+1,vl+1)-h(1,1)
           !
           do i=0,npoints 
              write(out,"(i8,2f18.8,' ||')") i,phil_(i),dphil_(i)
           enddo
        endif 
        !
     enddo
     !
     deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,rho_poten,rho_extF)
     deallocate(phil_,phir_,dphil_,dphir_,phivphi_)
     deallocate(poten,mu_rr)
     deallocate(h,ener)
     deallocate(psi,dpsi,psi_,dpsi_,phi,dphi)
     call ArrayStop('psi-phi-Fourier')
     call ArrayStop('h-Fourier')
     call ArrayStop('h-Fourier-phi')
     !
  end subroutine ME_Fourier

  !
  ! Matrix elements with Legendre-eigenfunctions 
  !
  subroutine ME_Legendre(vmax,maxorder,rho_b_,isingular,npoints,drho,poten,mu_rr,icoord,verbose,g_numerov,energy)
   !
   integer(ik),intent(in) :: vmax,maxorder,npoints,isingular
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:vmax,0:vmax)
   real(ark),intent(out)    :: energy(0:vmax)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten(0:npoints),mu_rr(0:npoints),drho(0:npoints,3)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose   ! Verbosity level
   !
   real(ark)            :: rho,rhostep,potmin,C_l,C_r,phi_rho(vmax+1),dphi_rho(vmax+1),ddphi_rho(vmax+1)
   real(ark)            :: psipsi_t,characvalue,rho_b(2),h_t,sigma_t,sigma,rms
   !
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,n,imin,io_slot
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:),phi(:),dphi(:),x(:),sinrho(:)
   real(ark),allocatable :: L(:,:),dL(:,:)
   real(rk),allocatable  :: h(:,:),ener(:),psi(:,:),dpsi(:,:),ddpsi(:,:)
   !
   character(len=cl)    :: unitfname 
     !
     if (verbose>=3) write (out,"(/20('*'),' Legendre real functions primitive matrix elements calculations')")
     !
     ! global variables 
     !
     allocate(phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints), &
              phivphi(0:npoints),rho_kinet(0:npoints),rho_poten(0:npoints),rho_extF(0:npoints),&
              x(0:npoints),sinrho(0:npoints),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     allocate(psi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Legendre',alloc,size(psi),kind(psi))
     allocate(dpsi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Legendre',alloc,size(dpsi),kind(dpsi))
     allocate(ddpsi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Legendre',alloc,size(ddpsi),kind(ddpsi))
     !
     allocate(h(vmax+1,vmax+1),ener(vmax+1),phi(vmax+1),dphi(vmax+1),stat=alloc)
     call ArrayStart('h-Legendre',alloc,size(h),kind(h))
     call ArrayStart('h-Legendre',alloc,size(ener),kind(ener))
     call ArrayStart('h-Legendre-phi',alloc,size(phi),kind(phi))
     call ArrayStart('h-Legendre-phi',alloc,size(dphi),kind(dphi))
     !
     allocate(L(0:npoints,0:vmax),dL(0:npoints,0:vmax),stat=alloc)
     call ArrayStart('Legendre',alloc,size(L),kind(L))
     call ArrayStart('Legendre',alloc,size(dL),kind(dL))
     !
     rho_b = rho_b_
     !
     ! step size 
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
     !
     ! Do some reporting
     !
     if (verbose>=3) then 
         write (out,"('vmax = ',i8)") vmax
         write (out,"('maxorder = ',i8)") maxorder
         write (out,"('icoord = ',i4)") icoord
         write (out,"('rho_b (x) = ',2f12.4)") rho_b(1:2) !*180.0_ark/pi
         write (out,"('rhostep (x) = ',2f12.4)") rhostep  !*180.0_ark/pi
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
         write(out,"('ML_box: pot_eff has no minimum',i8)") 
         stop 'ML_box: pot_eff has no minimum'
     endif 
     !
     ! define the x = cos(phi) coordinate 
     !
     do i=0,npoints
        !
        rho = rho_b(1)+real(i,kind=ark)*rhostep
        x(i) = cos(rho)
        sinrho(i) = sin(rho)
        !
     enddo
     !
     ! Evaluate the Legendre polynomials
     !
     call p_polynomial_value(npoints+1,vmax,x(0:),L(0:,0:))
     call p_polynomial_prime(npoints+1,vmax,x(0:),dL(0:,0:))
     !
     ! define the rho-type coordinate 
     !
     rho_kinet(:) = drho(:,1)
     rho_poten(:) = drho(:,2)
     rho_extF(:)  = drho(:,3)
     !
     inquire(iolength=rec_len) phil(:),dphil(:)
     !
     write(unitfname,"('Numerov basis set # ',i6)") icoord
     call IOStart(trim(unitfname),io_slot)
     !
     open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
     !
     !
     ! Matrix elements
     !
     !characvalue = maxval(enerslot(0:vmax))
     !
     do vl = 0,vmax
        !
        ! normalisation constant
        C_l = sqrt( real(2*vl+1,ark)/2.0_ark )
        !
        do i=0,npoints
           !
           rho = rho_b(1)+real(i,kind=ark)*rhostep
           !
           phil(i)  = L(i,vl)*C_l*sqrt(sin(rho))
           !dphil(i) = ( cos(rho)*L(i,vl)-sin(rho)**2*dL(i,vl) )*C_l
           !
           dphil(i) = -sin(rho)*dL(i,vl)*sqrt(sin(rho))*C_l
           !
        enddo
        !
        do vr = vl,vmax
            !
            ! norm constant
            C_r = sqrt( real(2*vr+1,ark)/2.0_ark )
            !
            do i=0,npoints
               !
               rho = rho_b(1)+real(i,kind=ark)*rhostep
               !
               phir(i)  = L(i,vr)*C_r*sqrt(sin(rho))
               !dphir(i) = ( cos(rho)*L(i,vr)-sin(rho)**2*dL(i,vr) )*C_r
               dphir(i) = -sin(rho)*dL(i,vr)*sqrt(sin(rho))*C_r
               !
            enddo
            !
            ! check orthagonality and noralisation
            !
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
            !phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)
            !
            phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)
            !
            psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
            !
            ! correction due to the derivatives at rho=0 and rho = Pi
            !
            !psipsi_t = psipsi_t - ( phil(npoints)*mu_rr(npoints)*phir(npoints) - phil(0)*mu_rr(0)*phir(0) )
            !
            ! Add the diagonal kinetic part to the tested mat. elem-s
            !
            h(vl+1,vr+1) = h(vl+1,vr+1) - 0.5_ark*psipsi_t
            !
            h(vr+1,vl+1) = h(vl+1,vr+1)
            !
        enddo
     enddo
     !
     call lapack_syev(h,ener)
     !
     energy(0:vmax) = ener(1:vmax+1)-ener(1)
     !
     write (out,"(/' Legendre-optimized energies are:')") 
     !
     do vl=0,vmax   
       write (out,"(i8,f18.8)") vl,energy(vl)
     enddo
     !
     do i=0,npoints
        !
        rho = rho_b(1)+real(i,kind=ark)*rhostep
        !
        do vl = 0,vmax
           !
           C_l = sqrt( real(2*vl+1,ark)/2.0_ark )
           !
           phi_rho(vl+1)  = L(i,vl)*C_l !*sqrt(sin(rho))
           !dphi_rho(vl+1) = ( cos(rho)*L(i,vl)*0.5_ark-sin(rho)**2*dL(i,vl) )*C_l
           dphi_rho(vl+1) = -sin(rho)*dL(i,vl)*C_l
           !
           ddphi_rho(vl+1) = -sin(rho)*dL(i,vl)*sqrt(sin(rho))*C_l
           !
        enddo
        !
        Psi (1:vmax+1,i)  = matmul(transpose(h(1:vmax+1,1:vmax+1)),  phi_rho(1:vmax+1))
        DPsi(1:vmax+1,i)  = matmul(transpose(h(1:vmax+1,1:vmax+1)), dphi_rho(1:vmax+1))
        DDPsi(1:vmax+1,i) = matmul(transpose(h(1:vmax+1,1:vmax+1)),ddphi_rho(1:vmax+1))
        !
     enddo
     !
     rms = 0
     sigma = 0.0_ark 
     rms   = 0.0_ark 
     characvalue = maxval(ener(:))
     !
     do vl = 0,vmax
        !
        phil(:)  =  Psi(vl+1,:)
        dphil(:) = dPsi(vl+1,:)
        !
        write (io_slot,rec=vl+1) (phil(i),i=0,npoints),(dphil(i),i=0,npoints)
        !
        do vr = vl,vmax
            !
            phir = Psi(vr+1,:)
            dphir = dPsi(vr+1,:)
            !
            ! Here we prepare integrals of the potential 
            ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
            ! obtained above by the Numerov
            !
            phivphi(:) = phil(:)*poten(:)*phir(:)*sinrho(:)
            !
            h_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
            !
            ! momenta-quadratic part 
            !
            phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)*sinrho(:)
            !
            psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
            !
            !psipsi_t = psipsi_t - (phil(npoints)*mu_rr(npoints)*phir(npoints)- phil(0)*mu_rr(0)*phir(0))
            !
            ! Add the diagonal kinetic part to the tested mat. elem-s
            !
            h_t = h_t - 0.5_ark*psipsi_t
            !
            ! check the solution
            !
            sigma_t =  abs(h_t)
            if (vl==vr) sigma_t =  abs(h_t-ener(vl+1))
            !
            sigma = max(sigma,sigma_t)
            rms = rms + sigma_t**2
            !
            ! Now we test the h_t = <vl|h|vr> matrix elements and check if Numerov cracked
            ! the Schroedinger all right
            if (vl/=vr.and.abs(h_t)>sqrt(small_)*abs(characvalue)*1e4) then 
               write(out,"('ME_numerov: wrong Numerovs solution for <',i4,'|H|',i4,'> = ',f20.10)") vl,vr,h_t
               stop 'ME_numerov: bad Numerov solution'
            endif 
            !
            if (vl==vr.and.abs(h_t-ener(vl+1))>sqrt(small_)*abs(characvalue)*1e4) then 
               write(out,"('ME_numerov: wrong <',i4,'|H|',i4,'> (',f16.6,') =/= energy (',f16.6,')')") vl,vr,h_t,ener(vl)
               stop 'ME_numerov: bad Numerov solution'
            endif 
            !
            ! Reporting the quality of the matrix elemenst 
            !
            if (verbose>=3) then 
              if (vl/=vr) then 
               write(out,"('<',i4,'|H|',i4,'> = ',e16.2,'<-',8x,'0.0',5x,'; <',i4,'|',i4,'> = ',e16.2,'<-',8x,'0.0')") & 
                                vl,vr,h_t,vl,vr,sigma_t
              else
               write(out,"('<',i4,'|H|',i4,'> = ',f16.6,'<-',f16.6,'; <',i4,'|',i4,'> = ',f16.6)")& 
                              vl,vr,h_t,ener(vl+1),vl,vr,sigma_t
              endif 
            endif 
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
               g_numerov(0,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
               !
               ! external field expansion
               !
               if (lambda==0) then 
                  phivphi(:) = phil(:)*phir(:)
               else
                  phivphi(:) = phil(:)*rho_extF(:)**lambda*phir(:)
               endif
               !
               g_numerov(3,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
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
               g_numerov(-1,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
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
               g_numerov(2,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
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
               g_numerov(1,lambda,vl,vr) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
               !
               if (vl/=vr) then
                  !
                  if (lambda==0) then 
                     phivphi(:) = dphil(:)*phir(:)
                  else
                     phivphi(:) = dphil(:)*rho_kinet(:)**lambda*phir(:)
                  endif
                  !
                  g_numerov(1,lambda,vr,vl) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
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
        enddo
     enddo
     !
     deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,rho_poten,rho_extF,x,sinrho,L,dL)
     deallocate(h,ener)
     deallocate(psi,dpsi,ddpsi,phi,dphi)
     call ArrayStop('psi-phi-Legendre')
     call ArrayStop('h-Legendre')
     call ArrayStop('h-Legendre-phi')
     call ArrayStop('Legendre')
     !
  end subroutine ME_Legendre



  !
  ! Matrix elements with Associate Legendre-eigenfunctions 
  !
  subroutine ME_Associate_Legendre(vmax,kmax,maxorder,rho_b_,isingular,npoints,drho,poten,mu_rr,mu_zz,icoord,verbose,g_numerov,energy)
   !
   implicit none
   integer(ik),intent(in) :: vmax,kmax,maxorder,npoints,isingular
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:vmax,0:vmax)
   real(ark),intent(out)    :: energy(0:vmax)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten(0:npoints),mu_rr(0:npoints),drho(0:npoints,3),mu_zz(0:npoints)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose   ! Verbosity level
   !
   real(ark)            :: rho,rhostep,potmin,C_l,C_r,ddphi_rho(vmax+1),zpe
   real(ark)            :: psipsi_t,characvalue,rho_b(2),h_t,sigma_t,sigma,rms,C1,C2,C3,C4,cross_prod,factor,mu_zz_t,mu_rr_t
   !
   integer(ik) :: vl,vr,nl,nr,il,ir,nmax,lambda,alloc,i,k,rec_len,n,imin,io_slot,lmax,nmax1
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:),phi(:)
   real(ark),allocatable :: L(:,:),dL(:,:),dphi(:),x(:),sinrho(:),cosrho(:),vect(:,:),psi(:,:),dpsi(:,:),phi_rho(:),dphi_rho(:)
   real(rk),allocatable  :: h(:,:),ener(:)
   !
   character(len=cl)    :: unitfname 
     !
     if (verbose>=3) write (out,"(/20('*'),' Legendre real functions primitive matrix elements calculations')")
     !
     ! global variables 
     !
     ! vibrational size is basis_size/(kmax+1)-1
     !
     nmax = (vmax+1)/(kmax+1)-1
     lmax = kmax + nmax
     !
     allocate(phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints), &
              phivphi(0:npoints),rho_kinet(0:npoints),rho_poten(0:npoints),rho_extF(0:npoints),&
              x(0:npoints),sinrho(0:npoints),cosrho(0:npoints),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     rho_b = rho_b_
     !
     ! step size 
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
     !
     ! Do some reporting
     !
     if (verbose>=3) then 
         write (out,"('vmax = ',i8)") vmax
         write (out,"('kmax = ',i8)") kmax
         write (out,"('lmax = ',i8)") lmax
         write (out,"('maxorder = ',i8)") maxorder
         write (out,"('icoord = ',i4)") icoord
         write (out,"('rho_b (x) = ',2f12.4)") rho_b(1:2)*180.0_ark/pi
         write (out,"('rhostep (x) = ',2f12.4)") rhostep  !*180.0_ark/pi
     endif 
     !
     if (kmax>lmax) then
       write(out,"('ME_Associate_Legendre error: illegal kmax>max ',2i8)") kmax,lmax
       stop 'ME_Associate_Legendre error: illegal kmax>lmax'
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
         write(out,"('ML_box: pot_eff has no minimum',i8)") 
         stop 'ML_box: pot_eff has no minimum'
     endif 
     !
     ! define the x = cos(phi) coordinate 
     !
     do i=0,npoints
        !
        rho = rho_b(1)+real(i,kind=ark)*rhostep
        x(i) = cos(rho)
        sinrho(i) = sin(rho)
        cosrho(i) = cos(rho)
        !
     enddo
     !
     !
     ! define the rho-type coordinate 
     !
     rho_kinet(:) = drho(:,1)
     rho_poten(:) = drho(:,2)
     rho_extF(:)  = drho(:,3)
     !
     inquire(iolength=rec_len) phil(:),dphil(:)
     !
     write(unitfname,"('Numerov basis set # ',i6)") icoord
     call IOStart(trim(unitfname),io_slot)
     !
     open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
     !
     nmax1 = nmax+1
     !
     allocate(h(nmax1,nmax1),ener(nmax1),phi(nmax1),dphi(nmax1),vect(nmax1,nmax1),stat=alloc)
     call ArrayStart('h-Legendre',alloc,size(h),kind(h))
     call ArrayStart('h-Legendre',alloc,size(ener),kind(ener))
     call ArrayStart('h-Legendre',alloc,size(vect),kind(vect))
     call ArrayStart('Legendre-phi',alloc,size(phi),kind(phi))
     call ArrayStart('Legendre-phi',alloc,size(dphi),kind(dphi))
     !
     allocate(psi(nmax1,0:npoints),stat=alloc)
     call ArrayStart('psi-Legendre',alloc,size(psi),kind(psi))
     allocate(dpsi(nmax1,0:npoints),stat=alloc)
     call ArrayStart('psi-Legendre',alloc,size(dpsi),kind(dpsi))     
     allocate(phi_rho(nmax1),stat=alloc)
     call ArrayStart('psi-Legendre',alloc,size(phi_rho),kind(phi_rho))     
     allocate(dphi_rho(nmax1),stat=alloc)
     call ArrayStart('psi-Legendre',alloc,size(dphi_rho),kind(dphi_rho))     
     !
     ! start a large loop over k
     !
     loop_k : do k = 0,kmax
       !
       if (verbose>=4) write(out,"(' K = ',i8)") k
       !
       allocate(L(0:npoints,0:k+nmax),dL(0:npoints,0:k+nmax),stat=alloc)
       call ArrayStart('Legendre',alloc,size(L),kind(L))
       call ArrayStart('Legendre',alloc,size(dL),kind(dL))
       !
       ! Evaluate the Legendre polynomials
       !
       ! normalised Associated Legendre polynomials
       call pmn_polynomial_value ( npoints+1, k+nmax, k, x(0:), L(0:,0:) )
       !
       dL = 0
       if (k>0) then
         call pmn_polynomial_value ( npoints+1, k+nmax, k-1, x(0:), dL(0:,0:) )
       else
         !
         ! un-normalised Legendre polyomials 
         !call p_polynomial_prime(npoints+1,nmax,x(0:),dL(0:,0:))
         !C_l = sqrt( real(2*vl+1,ark)/2.0_ark )
         !dL(:,vl) = dL(:,vl)*C_l
         !do vl = k,k+nmax
         !   dL(:,vl) = dL(:,vl)*C_l
         !enddo 
         if (k+nmax>0) then
           call pmn_polynomial_value ( npoints+1, k+nmax, 1_ik, x(0:), dL(0:,0:) )
           dL(0:,0:) = -dL(0:,0:)
         endif
         !
       endif
       !
       !do vl = k,vmax
       !   do i=0,npoints
       !     dL(:,:) = -int( (vl+k)*(vl-k+1),ik )*dL(i,vl)
       !   enddo
       !enddo
       !
       do nl = 0,nmax
          !
          vl = k+nl 
          !
          do i=0,npoints
             !
             rho = rho_b(1)+real(i,kind=ark)*rhostep
             phil(i)  = L(i,vl)*sqrt(sin(rho))
             dphil(i) = dL(i,vl)*sqrt(real((vl+k)*(vl-k+1),ark))
             !
          enddo
          !
          do nr = nl,nmax
              !
              vr = k+nr
              !
              do i=0,npoints
                 !
                 rho = rho_b(1)+real(i,kind=ark)*rhostep
                 phir(i)  = L(i,vr)*sqrt(sin(rho))
                 dphir(i) = dL(i,vr)*sqrt(real((vr+k)*(vr-k+1),ark))
                 !
              enddo
              !
              ! check orthagonality and noralisation
              !
              phivphi(:) = phil(:)*phir(:)
              psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              ! Here we prepare integrals of the potential 
              ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
              ! obtained above by the Numerov
              !
              phivphi(:) = phil(:)*poten(:)*phir(:)
              !
              h(nl+1,nr+1) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              ! momenta-quadratic part 
              !
              !phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)
              !
              !phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)
              !
              !C1 = real((vl+k)*(vl-k+1)*(vr+k)*(vr-k+1),ark)
              !C2 = real((vl+k)*(vl-k+1)*k,ark)
              !C3 = real((vr+k)*(vr-k+1)*k,ark)
              !C4 = real(k*k,ark)             
              !
              phivphi(:) =-mu_rr(:)*( dphil(:)*dphir(:)*sinrho(:)- &
                                      cosrho(:)*real(k,ark)*( dphil(:)*L(:,vr)+L(:,vl)*dphir(:) ) )
              !
              mu_rr_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              phivphi = real(k*k,ark) *mu_zz(:)*L(:,vl)*L(:,vr)
              !
              mu_zz_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              ! correction due to the derivatives at rho=0 and rho = Pi
              !
              !psipsi_t = psipsi_t - ( phil(npoints)*mu_rr(npoints)*phir(npoints) - phil(0)*mu_rr(0)*phir(0) )
              !
              ! Add the diagonal kinetic part to the tested mat. elem-s
              !
              h(nl+1,nr+1) = h(nl+1,nr+1) - 0.5_ark*mu_rr_t+0.5_ark*mu_zz_t
              !
              h(nr+1,nl+1) = h(nl+1,nr+1)
              !
          enddo
       enddo
       !
       call lapack_syev(h,ener)
       !
       write (out,"(/' Legendre-optimized energies are:')") 
       !
       zpe = ener(1)
       !
       do nl=0,nmax
         i = k*(nmax+1)+nl
         energy(i) = ener(nl+1)-zpe
         write (out,"(2i8,f18.8)") k,nl,energy(i)
       enddo

       !
       ! Schmidt orthogonalization to make eigenvectors orthogonal in ark
       !
       vect = h
       !
       do nl =  1,nmax1
         !
         cross_prod = sum(vect(:,nl)*vect(:,nl))
         !
         factor = 1.0_ark/sqrt(cross_prod)
         !
         vect(:,nl) = vect(:,nl)*factor
         !
         do nr = nl+1,nmax1
           !
           cross_prod = sum(vect(:,nl)*vect(:,nr))
           !
           vect(:,nr) = vect(:,nr)-cross_prod*vect(:,nl)
           ! 
         enddo
         !
         cross_prod = sum(vect(:,nl)*vect(:,nl))
         !
         factor = 1.0_ark/sqrt(cross_prod)
         vect(:,nl) = vect(:,nl)*factor
         !
       enddo 
       !
       do i=0,npoints
          !
          rho = rho_b(1)+real(i,kind=ark)*rhostep
          !
          do nl = 0,nmax
             !
             vl = k+nl
             !
             phil(nl)  =   L(i,vl)
             dphil(nl) =   dL(i,vl)*sqrt( real((vl+k)*(vl-k+1),ark) )
             !
             phi_rho(nl+1)  = L(i,vl)
             dphi_rho(nl+1) = dL(i,vl)*sqrt( real((vl+k)*(vl-k+1),ark) )
             !
             !ddphi_rho(vl+1) = dL(i,vl)*sqrt(sin(rho))
             !
          enddo
          !
          Psi (1:nmax1,i)  = matmul(transpose(vect), phi_rho)
          DPsi(1:nmax1,i)  = matmul(transpose(vect),dphi_rho)
          !
       enddo
       !
       sigma = 0
       rms   = 0
       characvalue = maxval(ener(:))
       !
       do nl = 0,nmax
          !
          vl = k+nl
          il = k*(nmax+1)+nl
          !
          phil(:)  =  Psi(nl+1,:)
          dphil(:) = dPsi(nl+1,:)
          !
          !phil(:)  = L(:,vl)!*sqrt(sinrho(:))
          !dphlr(:) = -dL(:,vl)*sqrt(real((vl+k)*(vl-k+1),ark))
          !
          write (io_slot,rec=il+1) (phil(i),i=0,npoints),(dphil(i),i=0,npoints)
          !
          do nr = nl,nmax
              !
              vr = k+nr
              ir = k*(nmax+1)+nr
              !
              phir = Psi(nr+1,:)
              dphir = dPsi(nr+1,:)
              !
              ! check orthagonality and noralisation
              !
              phivphi(:) = phil(:)*phir(:)*sinrho(:)
              psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              !phir(:)  = L(:,vr)!*sqrt(sinrho(:))
              !dphir(:) = -dL(:,vr)*sqrt(real((vr+k)*(vr-k+1),ark))
              !
              ! Here we prepare integrals of the potential 
              ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
              ! obtained above by the Numerov
              !
              phivphi(:) = phil(:)*poten(:)*phir(:)*sinrho(:)
              !
              h_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              ! momenta-quadratic part 
              !
              !phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)*sinrho(:)
              !
              !C1 = int((vl+k)*(vl-k+1)*(vr+k)*(vr-k+1),ik)
              !C2 = int((vl+k)*(vl-k+1)*k,ik)
              !C3 = int((vr+k)*(vr-k+1)*k,ik)
              !C4 = real(k*k,ark)
              !
              phivphi(:) =-mu_rr(:)*( dphil(:)*dphir(:)*sinrho(:) - & 
                          cosrho(:)*real(k,ark)*( dphil(:)*phir(:)+phil(:)*dphir(:) ) )
              !
              mu_rr_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              phivphi = real(k*k,ark)*mu_zz(:)*phil(:)*phir(:)
              !
              mu_zz_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              !psipsi_t = psipsi_t - (phil(npoints)*mu_rr(npoints)*phir(npoints)- phil(0)*mu_rr(0)*phir(0))
              !
              ! Add the diagonal kinetic part to the tested mat. elem-s
              !
              h_t = h_t - 0.5_ark*mu_rr_t+0.5_ark*mu_zz_t
              !
              ! check the solution
              !
              sigma_t =  abs(h_t)
              if (vl==vr) sigma_t =  abs(h_t-ener(nl+1))
              !
              sigma = max(sigma,sigma_t)
              rms = rms + sigma_t**2
              !
              ! Now we test the h_t = <vl|h|vr> matrix elements and check if Numerov cracked
              ! the Schroedinger all right
              if (nl/=nr.and.abs(h_t)>sqrt(small_)*abs(characvalue)*1e4) then 
                 write(out,"('ME_numerov: wrong Numerovs solution for <',i4,'|H|',i4,'> = ',f20.10)") nl,nr,h_t
                 stop 'ME_numerov: bad Numerov solution'
              endif 
              !
              if (nl==nr.and.abs(h_t-ener(nl+1))>sqrt(small_)*abs(characvalue)*1e4) then 
                 write(out,"('ME_numerov: wrong <',i4,'|H|',i4,'> (',f16.6,') =/= energy (',f16.6,')')") nl,nr,h_t,ener(vl)
                 stop 'ME_numerov: bad Numerov solution'
              endif 
              !
              ! Reporting the quality of the matrix elemenst 
              !
              if (verbose>=5) then 
                if (vl/=vr) then 
                 write(out,"('<',i4,'|H|',i4,'> = ',e16.2,'<-',8x,'0.0',5x,'; <',i4,'|',i4,'> = ',e16.2,'<-',8x,'0.0')") & 
                                  nl,nr,h_t,nl,nr,sigma_t
                else
                 write(out,"('<',i4,'|H|',i4,'> = ',f16.6,'<-',f16.6,'; <',i4,'|',i4,'> = ',f16.6)")& 
                                nl,nr,h_t,ener(nl+1),nl,nr,sigma_t
                endif 
              endif 
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
                 g_numerov(0,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 !
                 ! external field expansion
                 !
                 if (lambda==0) then 
                    phivphi(:) = phil(:)*phir(:)
                 else
                    phivphi(:) = phil(:)*rho_extF(:)**lambda*phir(:)
                 endif
                 !
                 g_numerov(3,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 if (il/=ir) g_numerov(3,lambda,ir,il) = g_numerov(3,lambda,il,ir)
                 !
                 ! momenta-free in kinetic part 
                 !
                 if (lambda==0) then 
                    phivphi(:) = phil(:)*phir(:)
                 else
                    phivphi(:) = phil(:)*rho_kinet(:)**lambda*phir(:)
                 endif
                 !
                 g_numerov(-1,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 !
                 ! We also control the orthogonality of the basis set 
                 !
                 if (lambda==0) psipsi_t = g_numerov(0,lambda,vl,vr)
                 !
                 if (il/=ir) g_numerov(-1:0,lambda,ir,il) = g_numerov(-1:0,lambda,il,ir)
                 !
                 ! momenta-quadratic part 
                 !
                 if (lambda==0) then 
                    phivphi(:) =-dphil(:)*dphir(:)
                 else
                    phivphi(:) =-dphil(:)*rho_kinet(:)**lambda*dphir(:)
                 endif
                 !
                 g_numerov(2,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 !
                 if (nl/=nr) g_numerov(2,lambda,ir,il) = g_numerov(2,lambda,il,ir)
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
                 g_numerov(1,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 !
                 if (vl/=vr) then
                    !
                    if (lambda==0) then 
                       phivphi(:) = dphil(:)*phir(:)
                    else
                       phivphi(:) = dphil(:)*rho_kinet(:)**lambda*phir(:)
                    endif
                    !
                    g_numerov(1,lambda,ir,il) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                    !
                 endif 
                 !
                 if (verbose>=7) then 
                     write(out,"('g_numerov(0,',i4,i4,i4,') = ',f18.8)") lambda,nl,nr,g_numerov(0,lambda,nl,nr)
                     write(out,"('g_numerov(1,',i4,i4,i4,') = ',f18.8)") lambda,nl,nr,g_numerov(1,lambda,nl,nr)
                     write(out,"('g_numerov(2,',i4,i4,i4,') = ',f18.8)") lambda,nl,nr,g_numerov(2,lambda,nl,nr)
                     write(out,"('g_numerov(3,',i4,i4,i4,') = ',f18.8)") lambda,nl,nr,g_numerov(3,lambda,nl,nr)
                     if (nl/=nr) then 
                       write(out,"('g_numerov(0,',i4,i4,i4,') = ',f18.8)") lambda,nr,nl,g_numerov(0,lambda,nr,nl)
                       write(out,"('g_numerov(1,',i4,i4,i4,') = ',f18.8)") lambda,nr,nl,g_numerov(1,lambda,nr,nl)
                       write(out,"('g_numerov(2,',i4,i4,i4,') = ',f18.8)") lambda,nr,nl,g_numerov(2,lambda,nr,nl)
                       write(out,"('g_numerov(3,',i4,i4,i4,') = ',f18.8)") lambda,nr,nl,g_numerov(3,lambda,nr,nl)
                     endif 
                 endif 
                 !
              enddo 
              !
          enddo
       enddo
       !
       deallocate(L,dL)
       call ArrayStop('Legendre')
       !
     enddo loop_k
     !
     ! cleanup
     !
     deallocate(h,ener,phi,dphi,vect)
     call ArrayStop('h-Legendre')
     call ArrayStop('Legendre-phi')
     deallocate(psi,dpsi,phi_rho,dphi_rho)
     call ArrayStop('psi-Legendre')
     !
     deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,rho_poten,rho_extF,x,sinrho,cosrho)
     !
  end subroutine ME_Associate_Legendre


  !
  ! Matrix elements with Associate Legendre-eigenfunctions 
  !
  subroutine ME_sinrho_polynomial(vmax,kmax,maxorder,rho_b_,isingular,npoints,drho,poten,mu_rr,mu_zz,icoord,verbose,g_numerov,energy)
   !
   implicit none
   integer(ik),intent(in) :: vmax,kmax,maxorder,npoints,isingular
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:vmax,0:vmax)
   real(ark),intent(out)    :: energy(0:vmax)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten(0:npoints),mu_rr(0:npoints),drho(0:npoints,3),mu_zz(0:npoints)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose   ! Verbosity level
   !
   real(ark)            :: rho_,rhostep,potmin,C_l,C_r,zpe
   real(ark)            :: psipsi_t,characvalue,rho_b(2),h_t,sigma_t,sigma,rms,C1,C2,C3,C4,cross_prod,factor,mu_zz_t,mu_rr_t
   !
   integer(ik) :: vl,vr,nl,nr,il,ir,nmax,lambda,alloc,i,k,rec_len,n,imin,io_slot,lmax,nmax1
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:),phi(:)
   real(ark),allocatable :: L(:,:),dL(:,:),dphi(:),x(:),sinrho(:),cosrho(:),vect(:,:),rho(:),psi(:,:),dpsi(:,:),phi_rho(:),dphi_rho(:)
   real(rk),allocatable  :: h(:,:),ener(:)
   !
   character(len=cl)    :: unitfname 
     !
     if (verbose>=3) write (out,"(/20('*'),' Pseudo Associate Legendre functions primitive matrix elements calculations')")
     !
     ! global variables 
     !
     ! vibrational size is basis_size/(kmax+1)-1
     !
     nmax = (vmax+1)/(kmax+1)-1
     lmax = kmax + nmax
     !
     allocate(phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints), &
              phivphi(0:npoints),rho_kinet(0:npoints),rho_poten(0:npoints),rho_extF(0:npoints),&
              x(0:npoints),sinrho(0:npoints),cosrho(0:npoints),rho(0:npoints),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     rho_b = rho_b_
     !
     ! step size 
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
     !
     ! Do some reporting
     !
     if (verbose>=3) then 
         write (out,"('vmax = ',i8)") vmax
         write (out,"('kmax = ',i8)") kmax
         write (out,"('lmax = ',i8)") lmax
         write (out,"('maxorder = ',i8)") maxorder
         write (out,"('icoord = ',i4)") icoord
         write (out,"('rho_b (x) = ',2f12.4)") rho_b(1:2)*180.0_ark/pi
         write (out,"('rhostep (x) = ',2f12.4)") rhostep  !*180.0_ark/pi
     endif 
     !
     if (kmax>lmax) then
       write(out,"('ME_sinrho_polynomial error: illegal kmax>max ',2i8)") kmax,lmax
       stop 'ME_sinrho_polynomial error: illegal kmax>lmax'
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
         write(out,"('ML_box: pot_eff has no minimum',i8)") 
         stop 'ML_box: pot_eff has no minimum'
     endif 
     !
     ! define the x = cos(phi) coordinate 
     !
     do i=0,npoints
        !
        rho_ = rho_b(1)+real(i,kind=ark)*rhostep
        x(i) = cos(rho_)
        sinrho(i) = sin(rho_)
        cosrho(i) = cos(rho_)
        rho(i) = rho_
        !
     enddo
     !
     !
     ! define the rho-type coordinate 
     !
     rho_kinet(:) = drho(:,1)
     rho_poten(:) = drho(:,2)
     rho_extF(:)  = drho(:,3)
     !
     inquire(iolength=rec_len) phil(:),dphil(:)
     !
     write(unitfname,"('Numerov basis set # ',i6)") icoord
     call IOStart(trim(unitfname),io_slot)
     !
     open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
     !
     nmax1 = nmax+1
     !
     allocate(h(nmax1,nmax1),ener(nmax1),phi(nmax1),dphi(nmax1),vect(nmax1,nmax1),stat=alloc)
     call ArrayStart('h-Legendre',alloc,size(h),kind(h))
     call ArrayStart('h-Legendre',alloc,size(ener),kind(ener))
     call ArrayStart('h-Legendre',alloc,size(vect),kind(vect))
     call ArrayStart('Legendre-phi',alloc,size(phi),kind(phi))
     call ArrayStart('Legendre-phi',alloc,size(dphi),kind(dphi))
     !
     allocate(psi(nmax1,0:npoints),stat=alloc)
     call ArrayStart('psi-Legendre',alloc,size(psi),kind(psi))
     allocate(dpsi(nmax1,0:npoints),stat=alloc)
     call ArrayStart('psi-Legendre',alloc,size(dpsi),kind(dpsi))     
     allocate(phi_rho(nmax1),stat=alloc)
     call ArrayStart('psi-Legendre',alloc,size(phi_rho),kind(phi_rho))     
     allocate(dphi_rho(nmax1),stat=alloc)
     call ArrayStart('psi-Legendre',alloc,size(dphi_rho),kind(dphi_rho))     
     !
     ! start a large loop over k
     !
     loop_k : do k = 0,kmax
       !
       if (verbose>=4) write(out,"(' K = ',i8)") k
       !
       allocate(L(0:npoints,0:nmax),dL(0:npoints,0:nmax),stat=alloc)
       call ArrayStart('Legendre',alloc,size(L),kind(L))
       call ArrayStart('Legendre',alloc,size(dL),kind(dL))
       !
       ! Generate polynomial sqrt(sin(rho))*sin(rho)^k*L^k_n by orthogonalising L^k_n = cos(rho)^n
       !
       ! for the expansion coefficients of the polynomial wrt x = cos(rho) and we start with a diagonal form
       !
       L = 0
       dL = 0
       !
       do vl =  0,nmax
         !
         L(:,vl) = x(:)**vl
         !
         if (vl>0) dL(:,vl) = -real(vl,ark)*x(:)**(vl-1)*sinrho(:)
         !
         Psi(vl+1,:) = L(:,vl)*sqrt(sinrho(:))*sinrho(:)**k
         !
         !phivphi(:) = psi(vl+1,:)*psi(vl+1,:)
         !
         !factor = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
         !Psi(vl+1,:) = 1.0_ark/sqrt(factor)*Psi(vl+1,:)
         !L(:,vl) = L(:,vl)/sqrt(factor)
         !dL(:,vl) = dL(:,vl)/sqrt(factor)
         !
       enddo
       !
       ! building the overlap matrix and diagonalizing it
       !
       !do vl = 0,nmax
       !   !
       !   do vr = vl,nmax
       !       !
       !       ! check orthagonality and noralisation
       !       !
       !       phivphi(:) = psi(vl+1,:)*psi(vr+1,:)
       !       !
       !       h(vl+1,vr+1) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
       !       h(vr+1,vl+1) = h(vl+1,vr+1)
       !       !
       !   enddo
       !enddo
       !
       ! orthogonalisation using the weight sqrt(sin(rho))*sin(rho)^k
       !
       do vl =  0,nmax
         !
         phivphi(:) = psi(vl+1,:)*psi(vl+1,:)
         cross_prod = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
         !
         !cross_prod = sum(psi(vl+1,:)*psi(vl+1,:))*rhostep
         !
         factor = 1.0_ark/sqrt(cross_prod)
         !
         psi(vl+1,:) = psi(vl+1,:)*factor
         L(:,vl)  =  L(:,vl)*factor
         dL(:,vl) = dL(:,vl)*factor
         !
         do vr = vl+1,nmax
           !
           phivphi(:) = psi(vl+1,:)*psi(vr+1,:)
           cross_prod = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
           !
           !cross_prod = sum(psi(vl+1,:)*psi(vr+1,:))*rhostep
           !
           psi(vr+1,:) = psi(vr+1,:)-cross_prod*psi(vl+1,:)
           !
           phivphi(:) = psi(vr+1,:)*psi(vr+1,:)
           factor = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
           !
           !factor = sum(psi(vr+1,:)*psi(vr+1,:))*rhostep
           !
           factor = 1.0_ark/sqrt(factor)
           psi(vr+1,:) = psi(vr+1,:)*factor
           L(:,vr)   = ( L(:,vr)-cross_prod* L(:,vl))*factor
           dL(:,vr)  = (dL(:,vr)-cross_prod*dL(:,vl))*factor
           ! 
         enddo
         !
       enddo 
       !
       !call lapack_syev(h,ener)
       !
       !vect = h
       !
       ! orthogonalisation using the weight sqrt(sin(rho))*sin(rho)^k
       !
       !do vl =  1,nmax1
       !  !
       !  cross_prod = sum(vect(:,vl)*vect(:,vl))
       !  !
       !  factor = 1.0_ark/sqrt(cross_prod)
       !  !
       !  vect(:,vl) = vect(:,vl)*factor
       !  !
       !  do vr = vl+1,nmax1
       !    !
       !    cross_prod = sum(vect(:,vl)*vect(:,vr))
       !    !
       !    vect(:,vr) = vect(:,vr)-cross_prod*vect(:,vl)
       !    ! 
       !  enddo
       !  !
       !  cross_prod = sum(vect(:,vl)*vect(:,vl))
       !  !
       !  factor = 1.0_ark/sqrt(cross_prod)
       !  vect(:,vl) = vect(:,vl)*factor
       !  !
       !enddo 
       !
       !
       !do i=0,npoints
       !   !
       !   do vl = 0,nmax
       !      !
       !      phi_rho(vl+1)  = L(i,vl)
       !      dphi_rho(vl+1) = dL(i,vl)
       !      !
       !      !ddphi_rho(vl+1) = dL(i,vl)*sqrt(sin(rho))
       !      !
       !   enddo
       !   !
       !   !Psi (1:nmax1,i)  = matmul(transpose(vect), phi_rho)
       !   !DPsi(1:nmax1,i)  = matmul(transpose(vect),dphi_rho)
       !   !
       !   psi (1:nmax1,i)  = matmul(transpose(vect), phi_rho)
       !   dpsi(1:nmax1,i)  = matmul(transpose(vect),dphi_rho)
       !   !
       !enddo
       !!
       !do vl = 0,nmax
       !   !
       !   phivphi(:) = psi(vl+1,:)*psi(vl+1,:)*sinrho(:)**(2*k+1)
       !   !
       !   factor = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
       !   L(:,vl) = psi(vl+1,:)/sqrt(factor)
       !   dL(:,vl)= dpsi(vl+1,:)/sqrt(factor)
       !   !
       !enddo
       !
       do vl = 0,nmax
          !
          phil(:)  = L(:,vl)*sqrt(sinrho(:))*sinrho(:)**k
          dphil(:) = dL(:,vl)*sinrho(:)**k
          if (k>0) dphil(:) = dphil(:) + real(k,ark)*sinrho(:)**(k-1)*L(:,vl)*cosrho(:)
          !
          do vr = vl,nmax
              !
              phir(:)  = L(:,vr)*sqrt(sinrho(:))*sinrho(:)**k
              dphir(:) = dL(:,vr)*sinrho(:)**k
              if (k>0) dphir(:) = dphir(:) + real(k,ark)*sinrho(:)**(k-1)*L(:,vr)*cosrho(:)
              !
              ! check orthagonality and noralisation
              !
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
              phivphi(:) =-dphil(:)*mu_rr(:)*dphir(:)*sinrho(:)
              !
              !phivphi(:) =-mu_rr(:)*( dphil(:)*dphir(:)*sinrho(:)- &
              !                        cosrho(:)*real(k,ark)*( dphil(:)*L(:,vr)+L(:,vl)*dphir(:) ) )
              !
              mu_rr_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              mu_zz_t = 0
              !
              if (k>0) then 
                !
                phivphi = real(k*k,ark)*mu_zz(:)*L(:,vl)*L(:,vr)*sinrho(:)**(2*k-1)
                !
                mu_zz_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                !
              endif
              !
              ! correction due to the derivatives at rho=0 and rho = Pi
              !
              !psipsi_t = psipsi_t - ( phil(npoints)*mu_rr(npoints)*phir(npoints) - phil(0)*mu_rr(0)*phir(0) )
              !
              ! Add the diagonal kinetic part to the tested mat. elem-s
              !
              h(vl+1,vr+1) = h(vl+1,vr+1) - 0.5_ark*mu_rr_t+0.5_ark*mu_zz_t
              !
              h(vr+1,vl+1) = h(vl+1,vr+1)
              !
          enddo
       enddo
       !
       call lapack_syev(h,ener)
       !
       write (out,"(/' Optimized energies are:')") 
       !
       zpe = ener(1)
       !
       do vl=0,nmax
         i = k*(nmax+1)+vl
         energy(i) = ener(vl+1)-zpe
         write (out,"(2i8,f18.8)") k,vl,energy(i)
       enddo
       !
       ! Schmidt orthogonalization to make eigenvectors orthogonal in ark
       !
       vect = h
       !
       do vl =  1,nmax1
         !
         cross_prod = sum(vect(:,vl)*vect(:,vl))
         !
         factor = 1.0_ark/sqrt(cross_prod)
         !
         vect(:,vl) = vect(:,vl)*factor
         !
         do vr = vl+1,nmax1
           !
           cross_prod = sum(vect(:,vl)*vect(:,vr))
           !
           vect(:,vr) = vect(:,vr)-cross_prod*vect(:,vl)
           ! 
         enddo
         !
         cross_prod = sum(vect(:,vl)*vect(:,vl))
         !
         factor = 1.0_ark/sqrt(cross_prod)
         vect(:,vl) = vect(:,vl)*factor
         !
       enddo 
       !
       do i=0,npoints
          !
          do vl = 0,nmax
             !
             phi_rho(vl+1)  = L(i,vl)
             !dphi_rho(vl+1) = dL(i,vl)
             dphi_rho(vl+1) = dL(i,vl)*sinrho(i)**k
             if (k>0) dphi_rho(vl+1)= dphi_rho(vl+1) + real(k,ark)*sinrho(i)**(k-1)*L(i,vl)*cosrho(i)
             !
          enddo
          !
          Psi (1:nmax1,i)  = matmul(transpose(vect), phi_rho)
          DPsi(1:nmax1,i)  = matmul(transpose(vect),dphi_rho)
          !
       enddo
       !
       sigma = 0
       rms   = 0
       characvalue = maxval(ener(:))
       !
       do vl = 0,nmax
          !
          il = k*(nmax+1)+vl
          !
          phil(:)  =  Psi(vl+1,:)
          dphil(:) = dPsi(vl+1,:)
          !
          write (io_slot,rec=il+1) (phil(i),i=0,npoints),(dphil(i),i=0,npoints)
          !
          do vr = vl,nmax
              !
              ir = k*(nmax+1)+vr
              !
              phir = Psi(vr+1,:)
              dphir = dPsi(vr+1,:)
              !
              ! check orthagonality and noralisation
              !
              phivphi(:) = phil(:)*phir(:)*sinrho(:)**(2*k+1)
              psipsi_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              ! Here we prepare integrals of the potential 
              ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
              ! obtained above by the Numerov
              !
              phivphi(:) = phil(:)*poten(:)*phir(:)*sinrho(:)**(2*k+1)
              !
              h_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              ! momenta-quadratic part 
              !
              phivphi(:) =-mu_rr(:)*dphil(:)*dphir(:)*sinrho(:)
              !
              mu_rr_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
              !
              mu_zz_t = 0
              !
              if (k>0) then 
                !
                phivphi = real(k*k,ark)*mu_zz(:)*phil(:)*phir(:)*sinrho(:)**(2*k-1)
                !
                mu_zz_t = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                !
              endif
              !
              ! Add the diagonal kinetic part to the tested mat. elem-s
              !
              h_t = h_t - 0.5_ark*mu_rr_t+0.5_ark*mu_zz_t
              !
              ! check the solution
              !
              sigma_t =  abs(h_t)
              if (vl==vr) sigma_t =  abs(h_t-ener(vl+1))
              !
              sigma = max(sigma,sigma_t)
              rms = rms + sigma_t**2
              !
              ! Now we test the h_t = <vl|h|vr> matrix elements and check if Numerov cracked
              ! the Schroedinger all right
              if (vl/=vr.and.abs(h_t)>sqrt(small_)*abs(characvalue)*1e4) then 
                 write(out,"('ME_numerov: wrong Numerovs solution for <',i4,'|H|',i4,'> = ',f20.10)") vl,vr,h_t
                 stop 'ME_numerov: bad Numerov solution'
              endif 
              !
              if (vl==vr.and.abs(h_t-ener(vl+1))>sqrt(small_)*abs(characvalue)*1e4) then 
                 write(out,"('ME_numerov: wrong <',i4,'|H|',i4,'> (',f16.6,') =/= energy (',f16.6,')')") vl,vr,h_t,ener(vl+1)
                 stop 'ME_numerov: bad Numerov solution'
              endif 
              !
              ! Reporting the quality of the matrix elemenst 
              !
              if (verbose>=5) then 
                if (vl/=vr) then 
                 write(out,"('<',i4,'|H|',i4,'> = ',e16.2,'<-',8x,'0.0',5x,'; <',i4,'|',i4,'> = ',e16.2,'<-',8x,'0.0')") & 
                                  vl,vr,h_t,vl,vr,sigma_t
                else
                 write(out,"('<',i4,'|H|',i4,'> = ',f16.6,'<-',f16.6,'; <',i4,'|',i4,'> = ',f16.6)")& 
                                vl,vr,h_t,ener(vl+1),vl,vr,sigma_t
                endif 
              endif 
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
                 g_numerov(0,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 !
                 ! external field expansion
                 !
                 if (lambda==0) then 
                    phivphi(:) = phil(:)*phir(:)
                 else
                    phivphi(:) = phil(:)*rho_extF(:)**lambda*phir(:)
                 endif
                 !
                 g_numerov(3,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 if (il/=ir) g_numerov(3,lambda,ir,il) = g_numerov(3,lambda,il,ir)
                 !
                 ! momenta-free in kinetic part 
                 !
                 if (lambda==0) then 
                    phivphi(:) = phil(:)*phir(:)
                 else
                    phivphi(:) = phil(:)*rho_kinet(:)**lambda*phir(:)
                 endif
                 !
                 g_numerov(-1,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 !
                 ! We also control the orthogonality of the basis set 
                 !
                 if (lambda==0) psipsi_t = g_numerov(0,lambda,vl,vr)
                 !
                 if (il/=ir) g_numerov(-1:0,lambda,ir,il) = g_numerov(-1:0,lambda,il,ir)
                 !
                 ! momenta-quadratic part 
                 !
                 if (lambda==0) then 
                    phivphi(:) =-dphil(:)*dphir(:)
                 else
                    phivphi(:) =-dphil(:)*rho_kinet(:)**lambda*dphir(:)
                 endif
                 !
                 g_numerov(2,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 !
                 if (vl/=vr) g_numerov(2,lambda,ir,il) = g_numerov(2,lambda,il,ir)
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
                 g_numerov(1,lambda,il,ir) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                 !
                 if (vl/=vr) then
                    !
                    if (lambda==0) then 
                       phivphi(:) = dphil(:)*phir(:)
                    else
                       phivphi(:) = dphil(:)*rho_kinet(:)**lambda*phir(:)
                    endif
                    !
                    g_numerov(1,lambda,ir,il) = simpsonintegral_ark(npoints,rho_b(2)-rho_b(1),phivphi)
                    !
                 endif 
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
          enddo
       enddo
       !
       deallocate(L,dL)
       call ArrayStop('Legendre')
       !
     enddo loop_k
     !
     ! cleanup
     !
     deallocate(h,ener,phi,dphi,vect)
     call ArrayStop('h-Legendre')
     call ArrayStop('Legendre-phi')
     deallocate(psi,dpsi,phi_rho,dphi_rho)
     call ArrayStop('psi-Legendre')
     !
     deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,rho_poten,rho_extF,x,sinrho,cosrho,rho)
     !
     if (verbose>=3) write (out,"(/20('*'),' ... done!')")
     !
  end subroutine ME_sinrho_polynomial




subroutine p_polynomial_value ( m, n, x, v )

!*****************************************************************************80
!
!! P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(n,1) = 1.
!    P(n,-1) = (-1)^N.
!    | P(n,x) | <= 1 in [-1,1].
!
!    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!    The Legendre polynomials are orthogonal under the inner product defined
!    as integration from -1 to 1:
!
!      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
!        = 0 if I =/= J
!        = 2 / ( 2*I+1 ) if I = J.
!
!    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
!    where
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
!
!    The formula is:
!
!      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Differential equation:
!
!    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,x) =      1
!    P( 1,x) =      1 X
!    P( 2,x) = (    3 X^2 -       1)/2
!    P( 3,x) = (    5 X^3 -     3 X)/2
!    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
!
!  Recursion:
!
!    P(0,x) = 1
!    P(1,x) = x
!    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
!
!    P'(0,x) = 0
!    P'(1,x) = 1
!    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
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
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = rk ) X(M), the evaluation points.
!
!    Output, real ( kind = rk ) V(M,0:N), the values of the Legendre polynomials 
!    of order 0 through N at the points X.
!
  implicit none
  !
  integer ( kind = ik ) m
  integer ( kind = ik ) n
  !
  integer ( kind = ik ) i
  real ( kind = ark ) v(m,0:n)
  real ( kind = ark ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0_ark

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = rk ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = rk ) *          v(1:m,i-2) ) &
               / real (     i,     kind = rk )
 
  end do
 
  return
  !
end subroutine p_polynomial_value



subroutine p_polynomial_prime ( m, n, x, vp )

!*****************************************************************************80
!
!! P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2012
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
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) VP(M,0:N), the values of the derivatives of the
!    Legendre polynomials of order 0 through N.
!
  implicit none

  integer ( kind = ik ) m
  integer ( kind = ik ) n

  integer ( kind = ik ) i
  real ( kind = ark ) v(m,0:n)
  real ( kind = ark ) vp(m,0:n)
  real ( kind = ark ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0_ark
  vp(1:m,0) = 0

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
  vp(1:m,1) = 1.0_ark
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = ark ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = ark ) *          v(1:m,i-2) ) &
               / real (     i,     kind = ark )
 
    vp(1:m,i) = ( real ( 2 * i - 1, kind = ark ) * ( v(1:m,i-1) &
                                                   + x(1:m) * vp(1:m,i-1) ) &
                - real (     i - 1, kind = ark ) *   vp(1:m,i-2)               ) &
                / real (     i,     kind = ark )
 
  end do
 
  return
end subroutine p_polynomial_prime


subroutine pmn_polynomial_value ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! PMN_POLYNOMIAL_VALUE: normalized Legendre polynomial Pmn(n,m,x).
!
!  Discussion:
!
!    The unnormalized associated Legendre functions P_N^M(X) have
!    the property that
!
!      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX 
!      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
!
!    By dividing the function by the square root of this term,
!    the normalized associated Legendre functions have norm 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2005
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
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
!
  implicit none

  integer ( kind = ik ) mm
  integer ( kind = ik ) n

  real ( kind = ark ) cx(mm,0:n)
  real ( kind = ark ) factor
  integer ( kind = ik ) j
  integer ( kind = ik ) m
  !real ( kind = ark ) ark_factorial
  real ( kind = ark ) x(mm)

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMN_POLYNOMIAL_VALUE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop 1
  end if
 
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMN_POLYNOMIAL_VALUE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    stop 1
  end if

  cx(1:mm,0:n) = 0.0_ark

  if ( m <= n ) then
    cx(1:mm,m) = 1.0_ark
    factor = 1.0_ark
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * factor * sqrt ( 1.0_ark - x(1:mm)**2 )
      factor = factor + 2.0_ark
    end do
  end if

  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = ark ) * cx(1:mm,m)
  end if

  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = ark ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = ark ) *           cx(1:mm,j-2) ) &
                 / real (     j - m,     kind = ark )
  end do
!
!  Normalization.
!
  do j = m, n
    factor = sqrt ( ( real ( 2 * j + 1, kind = ark ) * ark_factorial ( j - m ) ) &
      / ( 2.0_ark * ark_factorial ( j + m ) ) )
    cx(1:mm,j) = cx(1:mm,j) * factor
  end do

  return
end subroutine pmn_polynomial_value



function ark_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N. => ark
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = ark ) ark_factorial
  integer ( kind = ik ) i
  integer ( kind = ik ) n

  ark_factorial = 1.0_ark

  do i = 1, n
    ark_factorial = ark_factorial * real ( i, kind = ark )
  end do

  return
end function ark_factorial



  ! 
end module me_bnd

module me_bnd
  use accuracy
  use me_numer
  use timer
  use lapack
  !
  implicit none

  public ik, rk, out
  public degener_harm_q,ME_box,ME_Fourier

  integer(ik), parameter :: verbose     = 1                       ! Verbosity level

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

   real(ark)            :: rho,L,rhostep,potmin
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
  subroutine ME_Fourier(vmax,maxorder,rho_b_,isingular,npoints,drho,poten,mu_rr,icoord,iperiod,verbose,g_numerov,energy)
   !
   integer(ik),intent(in) :: vmax,maxorder,npoints,isingular
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:vmax,0:vmax)
   real(ark),intent(out)    :: energy(0:vmax)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten(0:npoints),mu_rr(0:npoints),drho(0:npoints,3)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose   ! Verbosity level
   integer(ik),intent(in) :: iperiod
   !
   integer(ik),parameter  :: Factor_FF=10 ! factor to increase the Fourier basis set size 
   !
   real(ark)            :: rho,L,rhostep,potmin
   real(ark)            :: psipsi_t,characvalue,rho_b(2)
   !
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,n,imin,io_slot,kl,kr,p,fmax
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:),phi(:),dphi(:)
   real(rk),allocatable :: h(:,:),ener(:),psi(:,:),dpsi(:,:)
   !
   character(len=cl)    :: unitfname 
    !
    if (verbose>=1) write (out,"(/20('*'),' Fourier real functions primitive matrix elements calculations')")
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
     fmax = vmax*factor_ff
     !
     allocate(psi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Fourier',alloc,size(psi),kind(psi))
     allocate(dpsi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Fourier',alloc,size(dpsi),kind(dpsi))
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
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=ark)
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
        Psi (1:vmax+1,i)  = matmul(transpose(h(1:fmax+1,1:vmax+1)),phi(1:fmax+1))
        DPsi(1:vmax+1,i)  = matmul(transpose(h(1:fmax+1,1:vmax+1)),dphi(1:fmax+1))
        !
     enddo
     !
     do vl = 0,vmax
        kl = (vl+1)/2*p
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
     deallocate(h,ener)
     deallocate(psi,dpsi,phi,dphi)
     call ArrayStop('psi-phi-Fourier')
     call ArrayStop('h-Fourier')
     call ArrayStop('h-Fourier-phi')
     !
  end subroutine ME_Fourier

  !
  ! Matrix elements with Fourier-eigenfunctions 
  !
  subroutine ME_Legendre_sqrt_sinrho(vmax,maxorder,rho_b_,npoints,drho,poten,mu_rr,icoord,verbose,g_numerov,energy)
   !
   integer(ik),intent(in) :: vmax,maxorder,npoints
   real(ark),intent(out)    :: g_numerov(-1:3,0:maxorder,0:vmax,0:vmax)
   real(ark),intent(out)    :: energy(0:vmax)
   !
   real(ark),intent(in) :: rho_b_(2)
   real(ark),intent(in) :: poten(0:npoints),mu_rr(0:npoints),drho(0:npoints,3)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   integer(ik),intent(in) :: verbose   ! Verbosity level
   !
   real(ark)            :: rho,rhostep,potmin
   real(ark)            :: psipsi_t,characvalue,rho_b(2)
   !
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,n,imin,io_slot,kl,kr,p,fmax
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:),phi(:),dphi(:),x(:)
   real(ark),allocatable :: L(:,:),dL(:,:)
   real(rk),allocatable  :: h(:,:),ener(:),psi(:,:),dpsi(:,:)
   !
   character(len=cl)    :: unitfname 
     !
     if (verbose>=3) write (out,"(/20('*'),' Legendre real functions primitive matrix elements calculations')")
     !
     ! global variables 
     !
     allocate(phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints), &
              phivphi(0:npoints),rho_kinet(0:npoints),rho_poten(0:npoints),rho_extF(0:npoints),&
              x(0:npoints),L(vmax,0:npoints),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     allocate(psi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Legendre',alloc,size(psi),kind(psi))
     allocate(dpsi(vmax+1,0:npoints),stat=alloc)
     call ArrayStart('psi-phi-Legendre',alloc,size(dpsi),kind(dpsi))
     !
     allocate(h(fmax+1,fmax+1),ener(fmax+1),phi(fmax+1),dphi(fmax+1),stat=alloc)
     call ArrayStart('h-Legendre',alloc,size(h),kind(h))
     call ArrayStart('h-Legendre',alloc,size(ener),kind(ener))
     call ArrayStart('h-Legendre-phi',alloc,size(phi),kind(phi))
     call ArrayStart('h-Legendre-phi',alloc,size(dphi),kind(dphi))
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
     ! define the x = cos(phi) coordinate 
     !
     do i=0,npoints
        !
        rho = real(i,kind=ark)*rhostep
        x(i) = cos(rho)
        !
     enddo
     !
     ! Evaluate the Legendre polynomials
     !
     call p_polynomial_value(npoints,vmax,x,L)
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
           phil(i)  = L(vl,i)
           dphil(i) = dL(vl,i)
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
               phir(i)  = L(vr,i)
               dphir(i) = dL(vr,i)
               !
            enddo
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
     write (out,"(/' Legendre-optimized energies are:')") 
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
           phi(i)  = L(vl,i)
           dphi(i) = dL(vl,i)
           !
        enddo
        !
        Psi (1:vmax+1,i)  = matmul(transpose(h(1:fmax+1,1:vmax+1)),phi(1:fmax+1))
        DPsi(1:vmax+1,i)  = matmul(transpose(h(1:fmax+1,1:vmax+1)),dphi(1:fmax+1))
        !
     enddo
     !
     do vl = 0,vmax
        kl = (vl+1)/2*p
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
     deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,rho_poten,rho_extF,x,L,dL)
     deallocate(h,ener)
     deallocate(psi,dpsi,phi,dphi)
     call ArrayStop('psi-phi-Legendre')
     call ArrayStop('h-Legendre')
     call ArrayStop('h-Legendre-phi')
     !
  end subroutine ME_Legendre_sqrt_sinrho




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

  v(1:m,0) = 1.0D+00

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





  ! 
end module me_bnd

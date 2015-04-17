module me_bnd
  use accuracy
  use me_numer
  use timer
  use lapack
  !
  implicit none

  public ik, rk, out
  public degener_harm_q,ME_box

  integer(ik), parameter :: verbose     = 1                       ! Verbosity level

  private

   character(len=cl):: solution_method  = '5POINTDIFFERENCES'
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
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,n,imin,io_slot
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
       energy(vl) = 0.5_ark*(vl+1)**2*mu_rr(imin)*pi**2/L**2
     enddo
     !
     !characvalue = maxval(enerslot(0:vmax))
     !
     do vl = 0,vmax
        !
        do i=0,npoints
           !
           rho = real(i,kind=ark)*rhostep
           !
           phil(i)  = sqrt(2.0_ark/L)*sin(real(vl+1,ark)*pi*rho/L)
           dphil(i) = sqrt(2.0_ark/L)*cos(real(vl+1,ark)*pi*rho/L)*real(vl+1,ark)*pi/L
           !
        enddo
        !
        write (io_slot,rec=vl+1) (phil(i),i=0,npoints),(dphil(i),i=0,npoints)
        !
        do vr = vl,vmax
            !
            if (vl==vr) then
                phir =  phil
               dphir = dphil
            else
              do i=0,npoints
                 !
                 rho = real(i,kind=ark)*rhostep
                 !
                 phir(i)  = sqrt(2.0_ark/L)*sin(real(vr+1,ark)*pi*rho/L)
                 dphir(i) = sqrt(2.0_ark/L)*cos(real(vr+1,ark)*pi*rho/L)*real(vr+1,ark)*pi/L
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
  ! Matrix elements with DVR-grid 
  !
  subroutine ME_DVR_grid(vmax,maxorder,rho_b_,isingular,npoints,drho,poten,mu_rr,icoord,periodic,verbose,g_numerov,energy)
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
   integer(ik) :: vl,vr,lambda,alloc,i,rec_len,n,imin,io_slot,ngrid
   !
   real(ark),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:),rho_poten(:),rho_extF(:)
   !real(ark),allocatable :: f(:),poten(:),mu_rr(:)
   real(rk),allocatable    :: vec(:),vibmat(:,:),vibener(:)
   real(rk)                :: epot,d2dr
   integer(ik) :: igrid,jgrid
   character(len=1)        :: rng,jobz
   real(rk)                :: vrange(2),zpe
   integer(ik)             :: irange(2),nroots

   character(len=cl)    :: unitfname 
    !
    if (verbose>=1) write (out,"(/20('*'),' DVR-grid calculations')")
     !
     ! global variables 
     !
     allocate(phil(0:npoints),phir(0:npoints),dphil(0:npoints),dphir(0:npoints), &
              phivphi(0:npoints),rho_kinet(0:npoints),rho_poten(0:npoints),rho_extF(0:npoints),stat=alloc)
     !
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     ngrid = npoints+1
     !
     allocate(vibmat(ngrid,ngrid),vibener(ngrid),vec(ngrid),stat=alloc)
     call ArrayStart('vibmat',alloc,size(vibmat),kind(vibmat))
     call ArrayStart('vibener',alloc,size(vibener),kind(vibmat))
     call ArrayStart('vec',alloc,size(vec),kind(vec))
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
       energy(vl) = 0.5_ark*(vl+1)**2*mu_rr(imin)*pi**2/L**2
     enddo
     !
     !characvalue = maxval(enerslot(0:vmax))
     !
     ! build the kinetic matrix elements

     !$omp parallel do private(igrid,f_rot,epot,f_l2,iL2,erot) shared(vibmat) schedule(guided)
     do igrid =0,npoints
       !
       if (verbose>=6) write(out,'("igrid = ",i0)') igrid
       !
       ! the diagonal term with the potential function
       !
       epot = rho_poten(igrid)
       !
       vibmat(igrid+1,igrid+1) = epot
       !
       d2dr = 30.0_ark
       !
       method_choice: select case(solution_method)
         case ("5POINTDIFFERENCES") ! default 
           vibmat(igrid,igrid) = vibmat(igrid,igrid) + d2dr
           !
           ! The nondiagonal matrix elemenets are:
           ! The vibrational kinetic energy operator will connect only the
           ! neighbouring grid points igrid+/1 and igrid+/2.
           !
           ! Comment by Lorenzo Lodi
           ! The following method corresponds to approximating the second derivative of the wave function
           ! psi''  by the 5-point finite difference formula:
           !
           ! f''(0) = [-f(-2h) +16*f(-h) - 30*f(0) +16*f(h) - f(2h) ] / (12 h^2)  + O( h^4 )
           !
          if (igrid>1) then
            vibmat(igrid+1,igrid+1-1) = -16.0_ark
            vibmat(igrid+1-1,igrid+1) = vibmat(igrid+1,igrid+1-1)
          endif
           !
          if (igrid>2) then
            vibmat(igrid+1,igrid+1-2) = 1.0_ark
            vibmat(igrid+1-2,igrid+1) = vibmat(igrid,igrid-2)
          endif

          case("SINC")   ! Experimental part using Colbert Miller sinc DVR (works only for uniform grids at the moment)
            vibmat(igrid+1,igrid+1) = vibmat(igrid+1,igrid+1) +(12._ark)* pi**2 / 3._ark

             do jgrid =igrid+1, ngrid
               vibmat(igrid+1,jgrid+1) = +(12._ark)*2._rk* real( (-1)**(igrid+jgrid), ark) / real(igrid - jgrid, ark)**2
               vibmat(jgrid+1,igrid+1) = vibmat(igrid+1,jgrid+1)
             enddo

          case default
           write(out, '(A)') 'Error: unrecognized solution method' // trim(solution_method)
           write(out, '(A)') 'Possible options are: '
           write(out, '(A)') '                      5POINTDIFFERENCES'
           write(out, '(A)') '                      SINC'
          end select method_choice

     enddo
     !$omp end parallel do
     !
     !
     ! diagonalize the vibrational hamiltonian using the DSYEVR routine from LAPACK
     !
     jobz = 'V'
     vrange(1) = -0.0_rk ; vrange(2) = (vmax)
     irange(1) = 1 ; irange(2) = min(vmax+1,Ngrid)
     nroots = Ngrid
     rng = 'A'
     if (vmax<1e8) then
        rng = 'V'
     elseif (vmax/=1e8) then
        rng = 'I'
     endif
     !
     call lapack_syevr(vibmat,vibener,rng=rng,jobz=jobz,iroots=nroots,vrange=vrange,irange=irange)
     !
     if (nroots<1) then
       nroots = 1
       vibener = 0
       vibmat = 0
       vibmat(1,1) = 1.0_rk
     endif
     !
     zpe = vibener(1)


     !
     ! print out the vibrational fields in the J=0 representaion
     if (verbose>=4) then
        write(out,'(/"Vibrational (contracted) energies: ")')
        write(out,'("    N        Energy/cm    State v"/)')
        do i = 1,nroots
          write(out,'(i5,f18.6," [ ",i4," ] ")') i,vibener(i)-zpe,i-1
        enddo
     endif

     !
     ! check the orthogonality of the basis
     !
     if (verbose>=6) then
       !
       write(out,'(/"Check the contracted basis for ortho-normality")')
       !
       !$omp parallel do private(ilevel,jlevel,psipsi_t) schedule(guided)
       do ilevel = 1,nroots
         do jlevel = 1,ilevel
           !
           psipsi_t  = sum(vibmat(:,ilevel)*vibmat(:,jlevel))
           !
           if (ilevel/=jlevel.and.abs(psipsi_t)>sqrt(small_)) then
              write(out,"('orthogonality is brocken : <',i4,'|',i4,'> (',f16.6,')')") ilevel,jlevel,psipsi_t
              stop 'Brocken orthogonality'
           endif
           !
           if (ilevel==jlevel.and.abs(psipsi_t-1.0_rk)>sqrt(small_)) then
              write(out,"('normalization is brocken:  <',i4,'|',i4,'> (',f16.6,')')") ilevel,jlevel,psipsi_t
              stop 'Brocken normalization'
           endif
           !
           ! Reporting the quality of the matrix elemenst
           !
           if (verbose>=6) then
             if (ilevel/=jlevel) then
               write(out,"('<',i4,'|',i4,'> = ',e16.2,'<-',8x,'0.0')") ilevel,jlevel,psipsi_t
             else
               write(out,"('<',i4,'|',i4,'> = ',f16.2,'<-',8x,'1.0')") ilevel,jlevel,psipsi_t
             endif
           endif
           !
         enddo
       enddo
       !$omp end parallel do
       !
     endif
     !
     ! dealocate some objects
     !
     !$omp parallel do private(ilevel,jlevel) schedule(guided)
     do ilevel = 1,nroots
       do jlevel = 1,ilevel
         !
         ! in the grid representation of the vibrational basis set
         ! the matrix elements are evaluated simply by a sumation of over the grid points
         !
         field%matelem(ilevel,jlevel)  = sum(vibmat(:,ilevel)*(field%gridvalue(:))*vibmat(:,jlevel))
         field%matelem(jlevel,ilevel) = field%matelem(ilevel,jlevel)
         !
       enddo
     enddo
     !$omp end parallel do


     !
     deallocate(vibmat,vibener)
     call ArrayStop('vibmat')
     call ArrayStop('vibener')
     deallocate(vec)
     call ArrayStop('vec')




     !
     do vl = 0,vmax
        !
        do i=0,npoints
           !
           rho = real(i,kind=ark)*rhostep
           !
           phil(i)  = sqrt(2.0_ark/L)*sin(real(vl+1,ark)*pi*rho/L)
           dphil(i) = sqrt(2.0_ark/L)*cos(real(vl+1,ark)*pi*rho/L)*real(vl+1,ark)*pi/L
           !
        enddo
        !
        write (io_slot,rec=vl+1) (phil(i),i=0,npoints),(dphil(i),i=0,npoints)
        !
        do vr = vl,vmax
            !
            if (vl==vr) then
                phir =  phil
               dphir = dphil
            else
              do i=0,npoints
                 !
                 rho = real(i,kind=ark)*rhostep
                 !
                 phir(i)  = sqrt(2.0_ark/L)*sin(real(vr+1,ark)*pi*rho/L)
                 dphir(i) = sqrt(2.0_ark/L)*cos(real(vr+1,ark)*pi*rho/L)*real(vr+1,ark)*pi/L
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
  end subroutine ME_DVR_grid


  ! 
end module me_bnd

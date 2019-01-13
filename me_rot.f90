module me_rot
  use accuracy
  use timer
  implicit none
  private

  public ME_rotation


  integer(ik),parameter :: verbose     = 5                          ! Verbosity level
  !
  integer(ik) :: Jmax


!
  contains
!

!                                                                    
! Rotational matrix elements  for symmetrized rigid rotor
!                                                                    
!
  subroutine ME_rotation(Jmax_,me) 
  !
  integer(ik),intent(in) :: Jmax_
  real(ark),intent(out)   ::  me(7,(jmax_+1)*(jmax_+2)/2,-2:2,0:1)
  !
  integer(ik) :: j,k1,k2,tau1,tau2,pm,jk,dk,alloc
  real(ark)    :: fjk,gjk,minus1,mejx,mejz,mejxy,mejyz,mejx2y2,mejy,mejxz,fjk_t
  !
  real(ark) ,allocatable :: jzmat(:,:,:)
  real(ark) ,allocatable :: jx2y2mat(:,:,:),jxymat(:,:,:),jxmat(:,:,:)
  real(ark) ,allocatable :: jymat(:,:,:),jxzmat(:,:,:),jyzmat(:,:,:)
    !
    !
    if (verbose>=5) write (out,"(20('*'),'  Rotational matrix elements calculations')")
    !
    Jmax = Jmax_
    !
    allocate(jzmat(0:jmax,-jmax:jmax,-jmax:jmax   ),stat=alloc)
    call ArrayStart('ME_rotation-jzmat',alloc,size(jzmat),kind(jzmat))
    allocate(jx2y2mat(0:jmax,-jmax:jmax,-jmax:jmax),stat=alloc)
    call ArrayStart('ME_rotation-jx2y2mat',alloc,size(jx2y2mat),kind(jx2y2mat))
    allocate(jxymat(0:jmax,-jmax:jmax,-jmax:jmax  ),stat=alloc)
    call ArrayStart('ME_rotation-jxymat',alloc,size(jxymat),kind(jxymat))
    allocate(jxmat(0:jmax,-jmax:jmax,-jmax:jmax   ),stat=alloc)
    call ArrayStart('ME_rotation-jxmat',alloc,size(jxmat),kind(jxmat))
    allocate(jymat(0:jmax,-jmax:jmax,-jmax:jmax   ),stat=alloc)
    call ArrayStart('ME_rotation-jymat',alloc,size(jymat),kind(jymat))
    allocate(jxzmat(0:jmax,-jmax:jmax,-jmax:jmax  ),stat=alloc)
    call ArrayStart('ME_rotation-jxzmat',alloc,size(jxzmat),kind(jxzmat))
    allocate(jyzmat(0:jmax,-jmax:jmax,-jmax:jmax  ),stat=alloc)
    call ArrayStart('ME_rotation-jyzmat',alloc,size(jyzmat),kind(jyzmat))
    !
    jxmat = 0 
    jymat = 0 
    jzmat = 0 
    jx2y2mat = 0 
    jxymat = 0 
    jxzmat = 0 
    jyzmat = 0 
    !
    do j =0,jmax 
       !
       do k1 = -j,j
          !
          do k2 = -j,j
             !
             pm = sign(1,k2-k1)
             !
             if (k1==k2) then 
                jzmat(j,k1,k2) = real(k1,ark)
             endif 
             !
             if (abs(k1-k2)==2) then
                !
                fjk_t = int( j*(j+1)-(k1+pm*2)*(k1+pm*1) ,hik)*int( j*(j+1)-k1*(k1+pm*1 ) ,hik)
                !
                fjk = 0.5_ark*sqrt( fjk_t )
                !
                jx2y2mat(j,k1,k2) = fjk
                jxymat(j,k1,k2) = real(pm,ark)*fjk
                !
                if (verbose>=7) write(out,"('j,k1,k2,pm,fjk = ',4i7,3e18.8)") j,k1,k2,pm,fjk,jx2y2mat(j,k1,k2),jxymat(j,k1,k2)
                !
             endif 
             !
             if (iabs(k1-k2).eq.1) then 
                !
                gjk = 0.5_ark*sqrt( real( j*(j+1) - k1*(k1+pm) , kind=ark ) )
                jxmat(j,k1,k2) = gjk
                jymat(j,k1,k2) = real(pm,ark)*gjk
                jxzmat(j,k1,k2) = gjk*real( 2*k1+pm ,ark)
                jyzmat(j,k1,k2) = gjk*real( ( 2*k1+pm )*pm ,ark) 
                !
             endif
             !
          enddo
       enddo
    enddo
    !
    me = 0 
    !
    do j =0,jmax 
       !
       do k1 = 0,j
          !
          jk=1+k1+(j*(j+1) )/2
          !
          do tau1 = 0,1
             !
             if (k1.ne.0.or.mod(j+tau1,2).eq.0) then
                !
                minus1 = -1.0_ark
                !
                if (mod(tau1,2).eq.0) minus1 = 1.0_ark
                !
                do k2=max(0,k1-2),min(j,k1+2)
                   do tau2 = 0,1
                      !
                      if (k2/=0.or.mod(j+tau2,2)==0) then
                         !
                         dk = k1-k2
                         !
                         mejx2y2 =         mult_rot_me(j,k1,k2,tau1,tau2,jx2y2mat)
                         mejy    =        -mult_rot_me(j,k1,k2,tau1,tau2,jymat   )
                         mejxz   =         mult_rot_me(j,k1,k2,tau1,tau2,jxzmat  )
                         mejx    =  minus1*mult_rot_me(j,k1,k2,tau1,tau2,jxmat   )
                         mejz    =  minus1*mult_rot_me(j,k1,k2,tau1,tau2,jzmat   )
                         mejxy   =  minus1*mult_rot_me(j,k1,k2,tau1,tau2,jxymat  )
                         mejyz   =  minus1*mult_rot_me(j,k1,k2,tau1,tau2,jyzmat  )
                         !
                         ! Here we check for the consistency 
                         !
                         if (tau1.ne.tau2) then 
                            !
                            if (mejx2y2/=0) then
                               write(out,"('ME_rotation: mejx2y2 <>0',5i7,e18.8)") j,k1,k2,tau1,tau2,mejx2y2
                               stop 'ME_rotation: mejx2y2 <>0' 
                            endif
                            ! 
                            if (mejy/=0) then 
                               write(out,"('ME_rotation: mejy <>0',5i7,e18.8)") j,k1,k2,tau1,tau2,mejy
                               stop 'ME_rotation: mejy <>0'
                            endif
                            ! 
                            if (mejxz/=0) then 
                               write(out,"('ME_rotation: mejxz <>0',5i7,e18.8)") j,k1,k2,tau1,tau2,mejxz
                               stop 'ME_rotation: mejxz <>0'
                            endif
                            !
                            ! here we store for different taus, but using only the first tau1 as the reference
                            !
                            me(1,jk,dk,tau1) =  mejx
                            me(3,jk,dk,tau1) =  mejz
                            me(5,jk,dk,tau1) =  mejxy
                            me(7,jk,dk,tau1) =  mejyz
                            ! 
                         else
                            ! 
                            if (mejx.ne.0) then 
                               write(out,"('ME_rotation: mejx <>0',5i7,e18.8)") j,k1,k2,tau1,tau2,mejx
                               stop 'ME_rotation: mejx <>0'
                            endif 
                            !
                            if (mejz.ne.0) then 
                               write(out,"('ME_rotation: mejz <>0',5i7,e18.8)") j,k1,k2,tau1,tau2,mejz
                               stop 'ME_rotation: mejz <>0'
                            endif 
                            !
                            if (mejxy.ne.0) then 
                               write(out,"('ME_rotation: mejxy <>0',5i7,e18.8)") j,k1,k2,tau1,tau2,mejxy
                               stop 'ME_rotation: mejxy <>0'
                            endif 
                            !
                            if (mejyz.ne.0) then 
                               write(out,"('ME_rotation: mejyz <>0',5i7,e18.8)") j,k1,k2,tau1,tau2,mejyz
                               stop 'ME_rotation: mejyz <>0'
                            endif
                            !
                            ! here we store for the same taus, but using only the first tau1 as the reference
                            !
                            me(2,jk,dk,tau1) =  mejy
                            me(4,jk,dk,tau1) =  mejx2y2
                            me(6,jk,dk,tau1) =  mejxz
                            !
                         endif
                         !
                      endif 
                   enddo
                enddo
             endif 
          enddo
       enddo
    enddo
    !
    deallocate(jzmat,jx2y2mat,jxymat,jxmat,jymat,jxzmat,jyzmat)
    !
    call ArrayStop('ME_rotation-jzmat')
    call ArrayStop('ME_rotation-jx2y2mat')
    call ArrayStop('ME_rotation-jxymat')
    call ArrayStop('ME_rotation-jxmat')
    call ArrayStop('ME_rotation-jymat')
    call ArrayStop('ME_rotation-jxzmat')
    call ArrayStop('ME_rotation-jyzmat')
    !
    if (verbose>=5) write (out,"(20('*'),'  Rotational matrix elements calculations/end')")
    !
  end subroutine  ME_rotation


!
! Miscelan. routime 
!
 function mult_rot_me(j,k1,k2,tau1,tau2,me) result (f)
   !

   integer(ik),intent(in) :: j,k1,k2,tau1,tau2
   real(ark),intent(in)    :: me(0:jmax,-jmax:jmax,-jmax:jmax)
   !
   real(ark) ::  sigma1,sigma2,minus0,f
   real(ark) ::  minus1,minus2,minus12,fk1,fk2
      !
      fk1 = 1.0_ark/sqrt(2.0_ark)
      if (k1.eq.0) fk1 = 0.5_ark
      !
      fk2 = 1.0_ark/sqrt(2.0_ark)
      if (k2.eq.0) fk2 = 0.5_ark
      !
      sigma1 = 1.0_ark
      if (tau1.eq.1.and.mod(k1,3).eq.1) sigma1 = -1.0_ark
      !
      sigma2 = 1.0_ark
      if (tau2.eq.1.and.mod(k2,3).eq.1) sigma2 = -1.0_ark
      !
      minus0 = sigma1*sigma2
      !
      minus1 = -1.0_ark
      if (mod(j+tau1+k1,2).eq.0) minus1 = 1.0_ark
      !
      minus2 = -1.0_ark
      if (mod(j+tau2+k2,2).eq.0) minus2 = 1.0_ark
      !
      minus12 = -1.0_ark
      if (mod(tau1+tau2+k1+k2,2).eq.0) minus12 = 1.0_ark
      !
      if (verbose>=6) write(out,"('j,k1,k2,minus1,minus1,minus12 = ',3i6,5e18.8)") j,k1,k2,minus1,minus1,minus12
      if (verbose>=6) write(out,"('me(j,pm k1, pm k2) = ',4e18.8)") me(j,k1, k2),me(j,-k1, k2),me(j,k1, -k2),me(j,-k1, -k2)
      !
      f = minus0*fk1*fk2*(        me(j,k1, k2) + minus1 *me(j,-k1,k2)   &
                         + minus2*me(j,k1,-k2) + minus12*me(j,-k1,-k2) )
      !                    
 end function mult_rot_me
 !

end module me_rot

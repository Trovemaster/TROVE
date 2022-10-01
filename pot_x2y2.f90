!
!  This unit defines all specific routines for a fouratomic molecule of ABCD type
!
module pot_x2y2
  use accuracy
  use moltype
  use lapack

  implicit none

  public MLdms_X2Y2_MB
  !
  private

  integer(ik), parameter :: verbose     = 4                        ! Verbosity level


  contains



 !===============================================================================
 !                   electric dipole moment section
 !===============================================================================

! HCCH Dipole using R,r1,r2,alpha1,alpha2,tau coordinates
!define cartesian components of the dipole moment in space-fixed system
 !
 ! mu_x has transformation properties of cos(tau) and r2+r3, a1+a2
 ! mu_y has transformation properties of cos(tau) and r2-r3  a1-a2
 ! mu_z has transformation properties of sin(tau) and r2-r3  a1-a2
 !
 recursive subroutine MLdms_X2Y2_MB(rank,ncoords,natoms,r,xyz,f)
    !
    implicit none
    !
    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  r(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: imu, iterm, ind(1:molec%ncoords),dimen
    real(ark)             :: mu_t,f_t,r21,r31,r1, r2, r3, alpha1, alpha2, tau, e1(3),e2(3),e3(3),tmat(3,3),&
                             tmatout(3,3),n1(3),n2(3),n3(3),xyz_(natoms,3)
    real(ark)             :: re1(1:3),re2(1:3),alphae(1:3),e0(3),costau, &
                             beta1(1:3),beta2(1:3),y(molec%ncoords,1:3), mu(3),xi(molec%ncoords),tau_,sintau,r0,tau1,tau2,tau_sign
    !
    integer(ik),parameter :: lspace = 150
    integer(ik) :: ierror,rank0,info
    real(rk)    :: dip_rk(3, 1), tmat_rk(3, 3), tsing(3), wspace(lspace),tol = -1.0d-12
    real(ark)    :: cosalpha1,cosalpha2,x1(3),x2(3),x3(3),n30,v2(3),v3(3),v20,v30,n0(3),e20,e30
    character(len=cl)  :: txt
    !
    !
    txt = 'Error: MLdms_X2Y2_MB'
    !
    x1(:) = xyz(1,:)-xyz(2,:)
    x2(:) = xyz(3,:)-xyz(1,:)
    x3(:) = xyz(4,:)-xyz(2,:)
    !
    r1 = sqrt(sum(x1(:)**2))
    r2 = sqrt(sum(x2(:)**2))
    r3 = sqrt(sum(x3(:)**2))
    !
    cosalpha1 = sum(-x1(:)*x2(:))/(r1*r2)
    cosalpha2 = sum( x1(:)*x3(:))/(r1*r3)
    !
    alpha1 = aacos(cosalpha1,txt)
    alpha2 = aacos(cosalpha2,txt)
    !
    n3(:) = -x1(:)/r1
    !
    n30 = sum(n3*n3)
    !
    n3(:) = n3/sqrt(n30)
    !
    v2(:) = MLvector_product(x2(:),n3(:))
    !
    v3(:) = MLvector_product(n3(:),x3(:))
    !
    v20 = sum(v2*v2)
    v30 = sum(v3*v3)
    !
    !setting tau
    if (v20>sqrt(small_).and.v30>sqrt(small_)) then !both v2 and v3 defined
       !
       v2 = v2/sqrt(v20)
       v3 = v3/sqrt(v30)
       !
       costau =-sum(v2*v3)
       !
       tau = aacos(costau,txt)
       !
    elseif (v20<sqrt(small_).and.v30>sqrt(small_)) then !v2 undefined
       !
       v3 = v3/sqrt(v30)
       tau=0
       !
    elseif (v20>sqrt(small_).and.v30<sqrt(small_)) then !v3 undefined
       v2 = v2/sqrt(v20)
       tau=0
       !
    else !both undefined
       !
       tau = 0
       !
    endif
    !End of setting tau
    !
    e2(:) = -MLvector_product(v2(:),n3(:))
    !
    e3(:) = MLvector_product(v3(:),n3(:))
    !
    e20 = sum(e2*e2)
    e30 = sum(e3*e3)
    !
    n1 = 0
    n2 = 0
    !
    !setting n1 and n2
    if (abs(tau-pi)<sqrt(small_).or.abs(tau+pi)<sqrt(small_)) then !tau=pi planar
      !
      n2(:)=-e2(:)
      !
    elseif(e30>small_.and.e20>small_) then
      !
      n1(:) =  (e2(:) + e3(:))
      n2(:) =  (e2(:) - e3(:))
      !
    elseif (e20<small_.and.e30>small_) then !v3 defined
      !
      n1(:) = e3(:)
      !
    elseif (e30<small_.and.e20>small_) then !v2 defined
      !
      n1(:) = e2(:)
      !
    elseif (e30<small_.and.e20<small_) then  !both v2 and v3 undefined
      !
      !define arbitrary direction for n1
      !
      n1 = (/ 1.0_ark,0.0_ark,0.0_ark /)
      !
    else
      !
      stop 'MLdms_X2Y2_MB: you should not be here'
      !
    endif
    !
    !if (e2(2)>0.or.e3(2)>0) then
    !  n1 = -n1
    !endif
    !
    !if (v2(2)>0.or.v3(2)>0) then
    !  n2 = -n2
    !endif
    !
    if (sum(n1(:)*n1(:))>small_) then !n1 defined
      !
      n1(:) = n1(:)/SQRT(sum(n1(:)*n1(:)))
      n2(:) = MLvector_product(n3(:),n1(:))
      !
    elseif (sum(n2(:)**2)>small_) then !n2 defined
      !
      n2(:) = n2(:)/SQRT(sum(n2(:)*n2(:)))
      n1(:) = MLvector_product(n2(:),n3(:))
      !
    else !both n1 and n2 are still undefined
      !
      write(6,"(' both n1 and n2 are  0 : n1 = ',3g12.5,' n2 = ',3g12.5)") n1,n2
      stop ' n1 = n2 = 0 '
      !
    endif
    !
    tmat(1,:) = n1(:)
    tmat(2,:) = n2(:)
    tmat(3,:) = n3(:)
    !
    tau_ = tau
    !
    re1(1:3)     = extF%coef(1,1:3)
    re2(1:3)     = extF%coef(2,1:3)
    alphae(1:3) = extF%coef(3,1:3)/rad
    !
    beta1(1:3)   = extF%coef(4,1:3)
    beta2(1:3)   = extF%coef(5,1:3)
    !
    y(1,:) = (r1 - re1(:)) * exp(-beta1(:) * (r1 - re1(:)) ** 2)
    y(2,:) = (r2 - re2(:)) * exp(-beta2(:) * (r2 - re2(:)) ** 2)
    y(3,:) = (r3 - re2(:)) * exp(-beta2(:) * (r3 - re2(:)) ** 2)
    !
    y(4,:) = (alphae(:) - alpha1)
    y(5,:) = (alphae(:) - alpha2)
    !
    !y(4,:) = sin(alphae(:)) - sin(alpha1)
    !y(5,:) = sin(alphae(:)) - sin(alpha2)
    !
    y(6,:) = cos(tau_)
    !
    mu=0
    mu_t=0
    !
    do imu = 1, 3
       !
       do iterm =  7, extF%nterms(imu)
          !
          ind(1:6) = extF%term(1:6, iterm, imu)
          xi(1:6) = y(1:6,imu) ** ind(1:6)
          !
          mu_t = product(xi(1:6))
          !
          if (ind(2)/=ind(3).or.ind(4)/=ind(5)) then
            !
            ind(2) = extF%term(3, iterm, imu)
            ind(3) = extF%term(2, iterm, imu)
            ind(4) = extF%term(5, iterm, imu)
            ind(5) = extF%term(4, iterm, imu)
            !
            xi(2:5) = y(2:5,imu) ** ind(2:5)
            !
            f_t = 1.0_ark
            if (imu/=1)  f_t = -1.0_ark
            !
            mu_t = mu_t + f_t*product(xi(1:6))
            !
            !
          endif
          !
          !if (imu==1) mu(imu) = mu(1)*cos(tau_*0.5_ark)
          !if (imu==2) mu(imu) = mu(2)*sin(tau_*0.5_ark)
          !
          mu(imu) = mu(imu) + extF%coef(iterm, imu)*mu_t
          !
          !if (imu==1) mu(imu) = mu(1)*cos(tau_*0.5_ark)
          !if (imu==2) mu(imu) = mu(2)*sin(tau_*0.5_ark)
          !
       end do
       !
    end do
    !
    mu(1) = mu(1)*cos(tau_*0.5_ark)
    mu(2) = mu(2)*sin(tau_*0.5_ark)
    !
    !write(out,"('tmat start')")
    !
    !if (verbose>=6) write(out,"('tmat')")
    !if (verbose>=6) write(out,*) tmat

    tmat = transpose(tmat)

    !if (verbose>=6) write(out,"('tmat trans')")
    !if (verbose>=6) write(out,*) tmat
    !if (verbose>=6) write(out,"('trans tmat')") tmat
    !tmat = inv(tmat)
    !if (verbose>=6) write(out,"('inv tmat')") tmat
    !
    dimen=3
    !
    !call MLinvmatark(tmat(1:dimen,1:dimen),tmatout(1:dimen,1:dimen),dimen,info)
    !
    f(1:3) = matmul(tmat,mu)
    !
    continue
    !
    if (abs(r(4))>1e-4   ) then
      continue
    endif
    !
    if (abs(r(5))>1e-4   ) then
      continue
    endif
    !
    if (abs(r(6))>1e-4  ) then
      continue
    endif
    !
    if (abs(r(4))>1e-4 .and. abs(r(6))>1e-4  ) then
      continue
    endif
    !
    if (abs(r(5))>1e-4 .and. abs(r(7))>1e-4  ) then
      continue
    endif
    !
    !if (verbose>=6) write(out,"('f is tmat times mu')") f
    !if (verbose>=6) write(out,"('mu')") mu
    !
    !f(1:2) = f(2:1:-1)
    !
  end subroutine MLdms_X2Y2_MB
  !
  !
end module pot_x2y2

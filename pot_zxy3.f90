!
!  This unit defines all specific routines for a fouratomic molecule of ZXY2 type
!
module pot_zxy3
  use accuracy
  use moltype
  !
  implicit none
  !
  public MLpoten_zxy3_sym,MLdms2xyz_zxy3_sym,MLpoten_zxy3_Nikitin
  !
  private
  !
  integer(ik), parameter :: verbose     = 3                          ! Verbosity level
  !
  contains




 function MLpoten_zxy3_sym(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          :: i,k_ind(9),ioper
   real(ark)            :: xi(9),chi(9,6),r1,r2,r3,r4,r1e,r2e,beta2,a,b,betae,tau12,tau23,tau13
   real(ark)            :: alpha12,alpha13,alpha14,alpha23,alpha24,alpha34,phi1,phi2,phi3
   real(ark)            :: cosbeta,beta312,beta412,beta413,cosa34,cosa12,cosa13,cosa23,cosa14,cosa24
   real(ark)            :: dx(4,3),txi(10),s1,s2,s3,s4,s5,s6,s7,s8,s9
   real(ark)            :: term
   character(len=cl)    :: txt
   !
   !
   r1e   = force(1)
   r2e   = force(2)
   betae = force(3)*pi/180.0_ark
   a     = force(4)
   b     = force(5)
   !
   r1      = local(1)
   r2      = local(2)
   r3      = local(3)
   r4      = local(4)
   !
   xi(1)=1.0_ark-exp(-a*(r1-r1e))
   xi(2)=1.0_ark-exp(-b*(r2-r2e))
   xi(3)=1.0_ark-exp(-b*(r3-r2e))
   xi(4)=1.0_ark-exp(-b*(r4-r2e))
   !
   select case(trim(molec%coords_transform))
   case default
      write (out,"('MLcoordinate_transform_func: coord. type ',a,' unknown')") trim(molec%coords_transform)
      stop 'MLcoordinate_transform_func - bad coord. type'
      !
   case('R-BETA-TAU')
      !
      !stop 'not tested R-BETA-TAU'
      !
      !
      if (size(local)==10) then 
        !
        dx(1,:)=xyz(2,:)-xyz(1,:)
        dx(2,:)=xyz(3,:)-xyz(1,:)
        dx(3,:)=xyz(4,:)-xyz(1,:)
        dx(4,:)=xyz(5,:)-xyz(1,:)
        !
        r1=sqrt(sum(dx(1,:)**2))
        r2=sqrt(sum(dx(2,:)**2))
        r3=sqrt(sum(dx(3,:)**2))
        r4=sqrt(sum(dx(4,:)**2))
        !
        xi(1)=1.0_ark-exp(-a*(r1-r1e))
        xi(2)=1.0_ark-exp(-b*(r2-r2e))
        xi(3)=1.0_ark-exp(-b*(r3-r2e))
        xi(4)=1.0_ark-exp(-b*(r4-r2e))
        !
        !  now get the angles
        !
        cosa12=(sum(dx(1,:)*dx(2,:)))/(r1*r2)
        cosa13=(sum(dx(1,:)*dx(3,:)))/(r1*r3)
        cosa14=(sum(dx(1,:)*dx(4,:)))/(r1*r4)
        cosa23=(sum(dx(2,:)*dx(3,:)))/(r2*r3)
        cosa24=(sum(dx(2,:)*dx(4,:)))/(r2*r4)
        cosa34=(sum(dx(3,:)*dx(4,:)))/(r3*r4)
        !
        alpha12=acos(cosa12)
        alpha13=acos(cosa13)
        alpha14=acos(cosa14)
        alpha23=acos(cosa23)
        alpha24=acos(cosa24)
        alpha34=acos(cosa34)
        !
        alpha34 = local(10)
        !
        cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
        beta312 = aacos(cosbeta,txt)
        !
        cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
        beta412 = aacos(cosbeta,txt)
        !
        xi(5)=alpha12-betae
        xi(6)=alpha13-betae
        xi(7)=alpha14-betae
        !
        phi1 = 2.0_ark*pi-(beta312+beta412)
        phi3 = beta312
        phi2 = beta412
        !
        xi(8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
        xi(9) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
        !
      else
        !
        !cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
        !beta312 = aacos(cosbeta,txt)
        !
        !cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
        !beta412 = aacos(cosbeta,txt)
        !
        !cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
        !alpha34 = aacos(cosa34,txt)
        !
        xi(1)=1.0_ark-exp(-a*(r1-r1e))
        xi(2)=1.0_ark-exp(-b*(r2-r2e))
        xi(3)=1.0_ark-exp(-b*(r3-r2e))
        xi(4)=1.0_ark-exp(-b*(r4-r2e))
        !
        xi(5) = local(5)-betae
        xi(6) = local(6)-betae
        xi(7) = local(7)-betae
        phi3 = 2.0_ark*pi-local(8)
        phi2 = local(9)
        phi1 = 2.0_ark*pi-phi3-phi2
        !
        !xi(8)=(2.0_ark*tau23-tau12-tau13)/sqrt(6.0_ark)
        !xi(9)=(tau12-tau13)/sqrt(2.0_ark)
        !
      endif 
      !
   case('R-BETA-SYM')
      !
      alpha12 = local(5)
      alpha13 = local(6)
      alpha23 = local(7)
      alpha14 = local(8)
      alpha24 = local(9)
      !
      if (size(local)==10) then 
        !
        dx(1,:)=xyz(2,:)-xyz(1,:)
        dx(2,:)=xyz(3,:)-xyz(1,:)
        dx(3,:)=xyz(4,:)-xyz(1,:)
        dx(4,:)=xyz(5,:)-xyz(1,:)
        !
        r1=sqrt(sum(dx(1,:)**2))
        r2=sqrt(sum(dx(2,:)**2))
        r3=sqrt(sum(dx(3,:)**2))
        r4=sqrt(sum(dx(4,:)**2))
        !
        xi(1)=1.0_ark-exp(-a*(r1-r1e))
        xi(2)=1.0_ark-exp(-b*(r2-r2e))
        xi(3)=1.0_ark-exp(-b*(r3-r2e))
        xi(4)=1.0_ark-exp(-b*(r4-r2e))
        !
        !  now get the angles
        !
        cosa12=(sum(dx(1,:)*dx(2,:)))/(r1*r2)
        cosa13=(sum(dx(1,:)*dx(3,:)))/(r1*r3)
        cosa14=(sum(dx(1,:)*dx(4,:)))/(r1*r4)
        cosa23=(sum(dx(2,:)*dx(3,:)))/(r2*r3)
        cosa24=(sum(dx(2,:)*dx(4,:)))/(r2*r4)
        cosa34=(sum(dx(3,:)*dx(4,:)))/(r3*r4)
        !
        alpha12=acos(cosa12)
        alpha13=acos(cosa13)
        alpha14=acos(cosa14)
        alpha23=acos(cosa23)
        alpha24=acos(cosa24)
        alpha34=acos(cosa34)
        !
        !alpha34 = local(10)
        !
        cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
        beta312 = aacos(cosbeta,txt)
        !
        cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
        beta412 = aacos(cosbeta,txt)
        !
        cosbeta = (cos(alpha34)-cos(alpha13)*cos(alpha14) )/(sin(alpha13)*sin(alpha14))
        beta413 = aacos(cosbeta,txt)
        !
        if (abs(beta312+beta412+beta413-2.0_rk*pi)>sqrt(small_)) then 
          write(out,"('MLpoten_zxy3_sym: beta312+beta412+beta413/=2.0_rk*pi for t1,t2,t3',4f16.8)") beta312,beta412,beta413,beta312+beta412+beta413
          stop
        endif
        !
      else
        !
        cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
        beta312 = aacos(cosbeta,txt)
        !
        cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
        beta412 = aacos(cosbeta,txt)
        !
        cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
        alpha34 = aacos(cosa34,txt)
        !
      endif 
      !
      xi(5)=alpha12-betae
      xi(6)=alpha13-betae
      xi(7)=alpha14-betae
      !
      phi1 = beta413 !2.0_ark*pi-(beta312+beta412)
      phi3 = beta312
      phi2 = beta412
      !
      xi(8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
      xi(9) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
      !
   end select 
   !
   txi(1:7) = xi(1:7)
   txi(8) = phi1
   txi(9) = phi3
   txi(10) = phi2
   !
   f = 0
   !
   do ioper =1,6
      !
      call symmetry_transformation_local(ioper,xi,chi(:,ioper))
      !
      !call symmetry_transformation_local10(ioper,txi,chi(:,ioper))
      !
   enddo
   !
   do i = 6,molec%parmax
     !
     k_ind(1:9) = molec%pot_ind(1:9,i)
     !
     term = 0
     !
     do ioper =1,6
        !
        !call symmetry_transformation_local(ioper,xi,chi)
        !
        term = term + product(chi(1:9,ioper)**k_ind(1:9))
        !
     enddo
     !
     f = f + term*force(i)
     !
   enddo
   !
  end function MLpoten_zxy3_sym
  !


  !Dipole moment for a CH3Cl-like molecule with the symmetry-adapted expansion generated on the fly 
  recursive subroutine MLdms2xyz_zxy3_sym(rank,ncoords,natoms,local,xyz0,f)
    !
    implicit none
    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz0(natoms,3)
    real(ark), intent(out) ::  f(rank)
    !
    integer(ik)            :: i,k_ind(9),ioper,imu,iterm,lwork,info,nsv
    real(ark)              :: xi(9),chi(9,6),r1,r2,r3,r4,r1e,r2e,beta2,a,b,betae,tau12,tau23,tau13
    real(ark)              :: alpha12,alpha13,alpha14,alpha23,alpha24,alpha34,phi1,phi2,phi3
    real(ark)              :: cosbeta,beta312,beta412,beta413,cosa34,cosa12,cosa13,cosa23,cosa14,cosa24
    real(ark)              :: dx(4,3),mu(3),nu(3)
    real(ark)              :: term(3),dip(3),tmat(4,3),amat(3,3),xyz(natoms,3)
    double precision       :: dmat(3,3),dipd(3,1),work(64*3),sv(3),svtol
    character(len=cl)      :: txt
    !
    r1e  = extF%coef(1,1)
    r2e  = extF%coef(2,1)
    betae = extF%coef(3,1)*pi/180.0_ark
    !
    r1      = local(1)
    r2      = local(2)
    r3      = local(3)
    r4      = local(4)
    !
    xi(1)=r1-r1e
    xi(2)=r2-r2e
    xi(3)=r3-r2e
    xi(4)=r4-r2e
    !
    xyz = xyz0
    !
!xyz(	1	,	1	)=	0.0000000000
!xyz(	1	,	2	)=	0.0000000000
!xyz(	1	,	3	)=	0.0000000000
!xyz(	2	,	1	)=	0.0000000000
!xyz(	2	,	2	)=	0.0000000000
!xyz(	2	,	3	)=	1.7760000000
!xyz(	3	,	1	)=	0.0000000000
!xyz(	3	,	2	)=	1.0411169188
!xyz(	3	,	3	)=	-0.3469390595
!xyz(	4	,	1	)=	0.9866556263
!xyz(	4	,	2	)=	-0.2880130220
!xyz(	4	,	3	)=	-0.3425123857
!xyz(	5	,	1	)=	-0.8901295519
!xyz(	5	,	2	)=	-0.5139165364
!xyz(	5	,	3	)=	-0.3425123857
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('MLcoordinate_transform_func: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'MLcoordinate_transform_func - bad coord. type'
       !
    case('R-BETA-TAU')
       !
       !stop 'not tested R-BETA-TAU'
       !
       if (size(local)==10) then 
         !
         dx(1,:)=xyz(2,:)-xyz(1,:)
         dx(2,:)=xyz(3,:)-xyz(1,:)
         dx(3,:)=xyz(4,:)-xyz(1,:)
         dx(4,:)=xyz(5,:)-xyz(1,:)
         !
         r1=sqrt(sum(dx(1,:)**2))
         r2=sqrt(sum(dx(2,:)**2))
         r3=sqrt(sum(dx(3,:)**2))
         r4=sqrt(sum(dx(4,:)**2))
         !
         xi(1)=r1-r1e
         xi(2)=r2-r2e
         xi(3)=r3-r2e
         xi(4)=r4-r2e
         !
         !  now get the angles
         !
         cosa12=(sum(dx(1,:)*dx(2,:)))/(r1*r2)
         cosa13=(sum(dx(1,:)*dx(3,:)))/(r1*r3)
         cosa14=(sum(dx(1,:)*dx(4,:)))/(r1*r4)
         cosa23=(sum(dx(2,:)*dx(3,:)))/(r2*r3)
         cosa24=(sum(dx(2,:)*dx(4,:)))/(r2*r4)
         cosa34=(sum(dx(3,:)*dx(4,:)))/(r3*r4)
         !
         alpha12=acos(cosa12)
         alpha13=acos(cosa13)
         alpha14=acos(cosa14)
         alpha23=acos(cosa23)
         alpha24=acos(cosa24)
         alpha34=acos(cosa34)
         !
         alpha34 = local(10)
         !
         cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
         beta312 = aacos(cosbeta,txt)
         !
         cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
         beta412 = aacos(cosbeta,txt)
         !
         xi(5)=alpha12-betae
         xi(6)=alpha13-betae
         xi(7)=alpha14-betae
         !
         phi1 = 2.0_ark*pi-(beta312+beta412)
         phi3 = beta312
         phi2 = beta412
         !
         xi(8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
         xi(9) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
         !
       else
         !
         xi(1)=r1-r1e
         xi(2)=r2-r2e
         xi(3)=r3-r2e
         xi(4)=r4-r2e
         !
         xi(5) = local(5)-betae
         xi(6) = local(6)-betae
         xi(7) = local(7)-betae
         phi3 = 2.0_ark*pi-local(8)
         phi2 = local(9)
         phi1 = 2.0_ark*pi-phi3-phi2
         !
         !xi(8)=(2.0_ark*tau23-tau12-tau13)/sqrt(6.0_ark)
         !xi(9)=(tau12-tau13)/sqrt(2.0_ark)
         !
       endif 
       !
    case('R-BETA-SYM')
       !
       alpha12 = local(5)
       alpha13 = local(6)
       alpha23 = local(7)
       alpha14 = local(8)
       alpha24 = local(9)
       !
       if (size(local)==10) then 
         !
         dx(1,:)=xyz(2,:)-xyz(1,:)
         dx(2,:)=xyz(3,:)-xyz(1,:)
         dx(3,:)=xyz(4,:)-xyz(1,:)
         dx(4,:)=xyz(5,:)-xyz(1,:)
         !
         r1=sqrt(sum(dx(1,:)**2))
         r2=sqrt(sum(dx(2,:)**2))
         r3=sqrt(sum(dx(3,:)**2))
         r4=sqrt(sum(dx(4,:)**2))
         !
         xi(1)=r1-r1e
         xi(2)=r2-r2e
         xi(3)=r3-r2e
         xi(4)=r4-r2e
         !
         !  now get the angles
         !
         cosa12=(sum(dx(1,:)*dx(2,:)))/(r1*r2)
         cosa13=(sum(dx(1,:)*dx(3,:)))/(r1*r3)
         cosa14=(sum(dx(1,:)*dx(4,:)))/(r1*r4)
         cosa23=(sum(dx(2,:)*dx(3,:)))/(r2*r3)
         cosa24=(sum(dx(2,:)*dx(4,:)))/(r2*r4)
         cosa34=(sum(dx(3,:)*dx(4,:)))/(r3*r4)
         !
         alpha12=acos(cosa12)
         alpha13=acos(cosa13)
         alpha14=acos(cosa14)
         alpha23=acos(cosa23)
         alpha24=acos(cosa24)
         alpha34=acos(cosa34)
         !
         !alpha34 = local(10)
         !
         cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
         beta312 = aacos(cosbeta,txt)
         !
         cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
         beta412 = aacos(cosbeta,txt)
         !
         cosbeta = (cos(alpha34)-cos(alpha13)*cos(alpha14) )/(sin(alpha13)*sin(alpha14))
         beta413 = aacos(cosbeta,txt)
         !
         if (abs(beta312+beta412+beta413-2.0_rk*pi)>sqrt(small_)) then 
           write(out,"('MLpoten_zxy3_sym: beta312+beta412+beta413/=2.0_rk*pi for t1,t2,t3',4f16.8)") beta312,beta412,beta413,beta312+beta412+beta413
           stop
         endif
         !
       else
         !
         cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
         beta312 = aacos(cosbeta,txt)
         !
         cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
         beta412 = aacos(cosbeta,txt)
         !
         cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
         alpha34 = aacos(cosa34,txt)
         !
       endif 
       !
       xi(5)=alpha12-betae
       xi(6)=alpha13-betae
       xi(7)=alpha14-betae
       !
       phi1 = beta413 !2.0_ark*pi-(beta312+beta412)
       phi3 = beta312
       phi2 = beta412
       !
       xi(8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
       xi(9) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
       !
    end select 
    !
    do ioper =1,6
       !
       !call symmetry_transformation_local2(ioper,xi,chi(:,ioper))
       !
       call symmetry_transformation_local(ioper,xi,chi(:,ioper))
       !
       !call symmetry_transformation_local10(ioper,txi,chi(:,ioper))
       !
    enddo
    !
    dip(1:3)=0
    !
    do imu=1,3
     !
     mu = 0
     !
     do iterm=4,extF%nterms(imu) !(iterm=1:3 are equilibrium constants)
       !
       k_ind(1:9) = extF%term(1:9,iterm,imu)
       !
       term = 0
       !
       do ioper =1,6
          !
          mu(imu) = product(chi(1:9,ioper)**k_ind(1:9))
          !
          call symmetry_transformation_dipole(ioper,mu,nu)
          !
          term(:) = term(:) + nu(:)
          !
       enddo
       !
       dip(:)=dip(:)+extF%coef(iterm,imu)*term(:)
       !
     enddo
    enddo
    ! 
    ! transformation of DMS to the TROVE Cartesian coordinate system 
    !
    !
    tmat(1, :) = dx(1,:) / r1
    tmat(2, :) = dx(2,:) / r2
    tmat(3, :) = dx(3,:) / r3
    tmat(4, :) = dx(4,:) / r4
    !
    amat(1,:) = 1._ark/sqrt(6._ark)*(2.0_ark*tmat(2,:)-tmat(3,:)-tmat(4,:))
    amat(2,:) = 1._ark/sqrt(2._ark)*(                  tmat(3,:)-tmat(4,:))
    amat(3,:) = tmat(1, :)
    !
    !forall(i=1:3) amat(i,:) = amat(i,:)/sqrt(sum(amat(i,:)**2))
    !
    mu = dip
    !
    call MLlinurark(3,amat,dip(:),f,info)
    !
    info = 1
    !
    if (info/=0) then 
      !
      dipd(1:3,1)= mu(1:3)
      dmat = amat
      !
      lwork=size(work)
      svtol=-1.0d-12
      call dgelss(3,3,1,dmat,3,dipd,3,sv,svtol,nsv,work,lwork,info)
      if (info/=0) then
       write(out,'(/a,1x,i3)') 'MLdms2xyz_zxy3_sym error: dgelss failed, info=',info
       stop 'MLdms2xyz_zxy3_sym error: dgelss failed'
      endif
      !
      f(1:3)=real(dipd(1:3,1),kind=ark)
      !
    endif
    !
  end subroutine MLdms2xyz_zxy3_sym

  function MLpoten_zxy3_Nikitin(ncoords,natoms,local,xyz,force) result(f) 
    !
    integer(ik),intent(in) ::  ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords)
    real(ark),intent(in)   ::  xyz(natoms,3)
    real(ark),intent(in)   ::  force(:)
    real(ark)              ::  f
    !
    real(ark) :: r1,r2,r3,r4,cosa12,alpha12,cosa13,alpha13,cosa23,alpha23,cosa14,alpha14,cosa24,alpha24,cosa34,alpha34,beta312,beta412,cosbeta
    !
    real(ark) :: dx(4,3)
    !
    integer(ik) :: k_ind(14),iterm,n
    !
    !real(ark) ::   betac, betaa,rc2,r02,a12,a13,a14,a23,a24,a34,t12,t13,t14,t23,t24,t34,vdump,cosae
    !
    real(ark) :: v
    real(ark) :: r1e,r2e,betae,a,b,y(14),xi(14)
    !
    character(len=cl)         :: txt
    !
    real(ark) :: ae= 1.9106332362490185563277142050315_ark
      !
      r1e   = force(1)
      r2e   = force(2)
      betae = force(3)*pi/180.0_ark
      a     = force(4)
      b     = force(5)
      !
      r1      = local(1)
      r2      = local(2)
      r3      = local(3)
      r4      = local(4)
      !
      select case(trim(molec%coords_transform))
      case default
         write (out,"('MLcoordinate_transform_func: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLcoordinate_transform_func - bad coord. type'
         !
      case('R-BETA-TAU')
         !
         stop 'R-BETA-TAU is not working yet'
         !
      case('R-BETA-SYM')
         !
         alpha12 = local(5)
         alpha13 = local(6)
         alpha23 = local(7)
         alpha14 = local(8)
         alpha24 = local(9)
         !
         if (size(local)==10) then 
           !
           dx(1,:)=xyz(2,:)-xyz(1,:)
           dx(2,:)=xyz(3,:)-xyz(1,:)
           dx(3,:)=xyz(4,:)-xyz(1,:)
           dx(4,:)=xyz(5,:)-xyz(1,:)
           !
           r1=sqrt(sum(dx(1,:)**2))
           r2=sqrt(sum(dx(2,:)**2))
           r3=sqrt(sum(dx(3,:)**2))
           r4=sqrt(sum(dx(4,:)**2))
           !
           !  now get the angles
           !
           cosa12=(sum(dx(1,:)*dx(2,:)))/(r1*r2)
           cosa13=(sum(dx(1,:)*dx(3,:)))/(r1*r3)
           cosa14=(sum(dx(1,:)*dx(4,:)))/(r1*r4)
           cosa23=(sum(dx(2,:)*dx(3,:)))/(r2*r3)
           cosa24=(sum(dx(2,:)*dx(4,:)))/(r2*r4)
           cosa34=(sum(dx(3,:)*dx(4,:)))/(r3*r4)
           !
           alpha12=acos(cosa12)
           alpha13=acos(cosa13)
           alpha14=acos(cosa14)
           alpha23=acos(cosa23)
           alpha24=acos(cosa24)
           alpha34=acos(cosa34)
           !
           alpha34 = local(10)
           !
           cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
           beta312 = aacos(cosbeta,txt)
           !
           cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
           beta412 = aacos(cosbeta,txt)
           !
         else
           !
           cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
           beta312 = aacos(cosbeta,txt)
           !
           cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
           beta412 = aacos(cosbeta,txt)
           !
           cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
           alpha34 = aacos(cosa34,txt)
           !
         endif 
         !
      end select 
      !
      y(1)=1.0_ark-exp(-a*(r1-r1e))
      y(2)=1.0_ark-exp(-b*(r2-r2e))
      y(3)=1.0_ark-exp(-b*(r3-r2e))
      y(4)=1.0_ark-exp(-b*(r4-r2e))
      !
      y(5)=cos(alpha12)-cos(betae)
      y(6)=cos(alpha13)-cos(betae)
      y(7)=cos(alpha14)-cos(betae)
      !
      y(8)=cos(beta312)-cos(2.0_ark*pi/3.0_ark)
      y(9)=cos(beta412)-cos(2.0_ark*pi/3.0_ark)
      !
      y(10)=sin(alpha12)-sin(betae)
      y(11)=sin(alpha13)-sin(betae)
      y(12)=sin(alpha14)-sin(betae)
      !
      y(13)=sin(beta312)-sin(2.0_ark*pi/3.0_ark)
      y(14)=sin(beta412)-sin(2.0_ark*pi/3.0_ark)
      !
      f = 0
      !
      n = size(force)
      !
      do iterm=6,n
        !
        k_ind(1:4) = molec%pot_ind(1:4,iterm)
        k_ind(5:9) = molec%pot_ind(5:9,iterm)/10
        k_ind(10:14) = mod(molec%pot_ind(5:9,iterm),10)
        !
        xi(1:14) = y(1:14)**k_ind(1:14)
        !
        f = f + force(iterm)*product(xi(1:14))
        !
      enddo
      !
      f = f*219474.6313705_rk
      !
  end function MLpoten_zxy3_Nikitin
  !

  subroutine symmetry_transformation_local2(ioper,src,dst)
    !
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    real(ark),intent(in)  :: src(1:9)
    real(ark),intent(out) :: dst(1:9)
    !
    real(ark)             :: repres(6,9,9),a,b,e,o
    !
    if (verbose>=5) write(out,"('symmetry_transformation_local/start')") 
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    dst(:) = src(:)
    !
    repres = 0 
    !
    repres(1,1,1) = 1.0_ark
    repres(1,2,2) = 1.0_ark
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    repres(1,6,6) = 1.0_ark
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !
    repres(2,1,1) = 1.0_ark
    repres(2,2,4) = 1.0_ark
    repres(2,3,2) = 1.0_ark
    repres(2,4,3) = 1.0_ark
    repres(2,5,7) = 1.0_ark
    repres(2,6,5) = 1.0_ark
    repres(2,7,6) = 1.0_ark
    repres(2,8,8) = -a
    repres(2,8,9) = -b
    repres(2,9,8) =  b
    repres(2,9,9) = -a
    !
    repres(3,1,1) = 1.0_ark
    repres(3,2,3) = 1.0_ark
    repres(3,3,4) = 1.0_ark
    repres(3,4,2) = 1.0_ark
    repres(3,5,6) = 1.0_ark
    repres(3,6,7) = 1.0_ark
    repres(3,7,5) = 1.0_ark
    repres(3,8,8) = -a
    repres(3,8,9) =  b
    repres(3,9,8) = -b
    repres(3,9,9) = -a
    !
    repres(4,1,1) = 1.0_ark
    repres(4,2,2) = 1.0_ark
    repres(4,3,4) = 1.0_ark
    repres(4,4,3) = 1.0_ark
    repres(4,5,5) = 1.0_ark
    repres(4,6,7) = 1.0_ark
    repres(4,7,6) = 1.0_ark
    repres(4,8,8) = 1.0_ark
    repres(4,8,9) = 0.0_ark
    repres(4,9,8) = 0.0_ark
    repres(4,9,9) = -1.0_ark
    !
    repres(5,1,1) = 1.0_ark
    repres(5,2,4) = 1.0_ark
    repres(5,3,3) = 1.0_ark
    repres(5,4,2) = 1.0_ark
    repres(5,5,7) = 1.0_ark
    repres(5,6,6) = 1.0_ark
    repres(5,7,5) = 1.0_ark
    repres(5,8,8) = -a
    repres(5,8,9) = -b
    repres(5,9,8) = -b
    repres(5,9,9) =  a
    !
    repres(6,1,1) = 1.0_ark
    repres(6,2,3) = 1.0_ark
    repres(6,3,2) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,6) = 1.0_ark
    repres(6,6,5) = 1.0_ark
    repres(6,7,7) = 1.0_ark
    repres(6,8,8) = -a
    repres(6,8,9) = b
    repres(6,9,8) = b
    repres(6,9,9) = a
    !
    if (ioper<0.or.ioper>6) then
      write (out,"('symmetry_transformation_local: operation ',i8,' unknown')") ioper
      stop 'symmetry_transformation_local - bad operation. type'
    endif
    !
    dst = matmul(repres(ioper,:,:),src) 
    !
    !dst(8) = repres(ioper,1,1)*src(8)+repres(ioper,1,2)*src(9)
    !dst(9) = repres(ioper,2,1)*src(8)+repres(ioper,2,2)*src(9)
           
  end subroutine symmetry_transformation_local2        
  

  subroutine symmetry_transformation_local10(ioper,src,dst)
    !
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    real(ark),intent(in)  :: src(1:10)
    real(ark),intent(out) :: dst(1:9)
    !
    real(ark)             :: repres(6,10,10),a,b,e,o,coord(10),phi1,phi2,phi3
    !
    if (verbose>=5) write(out,"('symmetry_transformation_local10/start')") 
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    repres = 0 
    !
    repres(1,1,1) = 1.0_ark
    repres(1,2,2) = 1.0_ark
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    repres(1,6,6) = 1.0_ark
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    repres(1,10,10) = 1.0_ark
    !
    repres(2,1,1) = 1.0_ark
    repres(2,2,4) = 1.0_ark
    repres(2,3,2) = 1.0_ark
    repres(2,4,3) = 1.0_ark
    repres(2,5,7) = 1.0_ark
    repres(2,6,5) = 1.0_ark
    repres(2,7,6) = 1.0_ark
    repres(2,8,10) = 1.0_ark
    repres(2,9,8) = 1.0_ark
    repres(2,10,9) = 1.0_ark
    !
    repres(3,1,1) = 1.0_ark
    repres(3,2,3) = 1.0_ark
    repres(3,3,4) = 1.0_ark
    repres(3,4,2) = 1.0_ark
    repres(3,5,6) = 1.0_ark
    repres(3,6,7) = 1.0_ark
    repres(3,7,5) = 1.0_ark
    repres(3,8,9) = 1.0_ark
    repres(3,9,10) = 1.0_ark
    repres(3,10,8) = 1.0_ark
    !
    repres(4,1,1) = 1.0_ark
    repres(4,2,2) = 1.0_ark
    repres(4,3,4) = 1.0_ark
    repres(4,4,3) = 1.0_ark
    repres(4,5,5) = 1.0_ark
    repres(4,6,7) = 1.0_ark
    repres(4,7,6) = 1.0_ark
    repres(4,8,8) = 1.0_ark
    repres(4,9,10) = 1.0_ark
    repres(4,10,9) = 1.0_ark
    !
    repres(5,1,1) = 1.0_ark
    repres(5,2,4) = 1.0_ark
    repres(5,3,3) = 1.0_ark
    repres(5,4,2) = 1.0_ark
    repres(5,5,7) = 1.0_ark
    repres(5,6,6) = 1.0_ark
    repres(5,7,5) = 1.0_ark
    repres(5,8,10) = 1.0_ark
    repres(5,9,9) = 1.0_ark
    repres(5,10,8) = 1.0_ark
    !
    repres(6,1,1) = 1.0_ark
    repres(6,2,3) = 1.0_ark
    repres(6,3,2) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,6) = 1.0_ark
    repres(6,6,5) = 1.0_ark
    repres(6,7,7) = 1.0_ark
    repres(6,8,9) = 1.0_ark
    repres(6,9,8) = 1.0_ark
    repres(6,10,10) = 1.0_ark
    !
    if (ioper<0.or.ioper>6) then
      write (out,"('symmetry_transformation_local10: operation ',i8,' unknown')") ioper
      stop 'symmetry_transformation_local10 - bad operation. type'
    endif
    !
    coord = matmul(repres(ioper,:,:),src) 
    !
    dst(1:7) = coord(1:7)
    !
    phi1 = coord(8)
    phi3 = coord(9)
    phi2 = coord(10)
    !
    dst(8) = 1.0_ark/sqrt(6.0_ark)*( 2.0_ark*phi1-phi2-phi3 )
    dst(9) = 1.0_ark/sqrt(2.0_ark)*(              phi2-phi3 )
    !
    !dst(8) = repres(ioper,1,1)*src(8)+repres(ioper,1,2)*src(9)
    !dst(9) = repres(ioper,2,1)*src(8)+repres(ioper,2,2)*src(9)
           
  end subroutine symmetry_transformation_local10 
  

  subroutine symmetry_transformation_local(ioper,src,dst)
    !
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    real(ark),intent(in)  :: src(1:9)
    real(ark),intent(out) :: dst(1:9)
    !
    real(ark)             :: repres(12,2,2),a,b,e,o
    !
    if (verbose>=5) write(out,"('symmetry_transformation_local/start')") 
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

      write (out,"('symmetry_transformation_local: operation ',i8,' unknown')") ioper
      stop 'symmetry_transformation_local - bad operation. type'
 
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
           
  end subroutine symmetry_transformation_local        


  
  subroutine symmetry_transformation_dipole(ioper,src,dst)
    implicit none 
    !
    integer,intent(in)    :: ioper  ! group operation  
    real(ark),intent(in)      :: src(1:3)
    real(ark),intent(out)     :: dst(1:3)
    !
    real(ark)         :: repres(6,3,3),a,b,e,o
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    repres = 0
    !
    repres(:,3,3) = 1.0_ark 
    !
    repres(1,1,1) = 1.0_ark
    repres(1,2,2) = 1.0_ark
    !
    repres(2,1,1) = -a
    repres(2,1,2) = -b
    repres(2,2,1) =  b
    repres(2,2,2) = -a
    !
    repres(3,1,1) = -a
    repres(3,1,2) =  b
    repres(3,2,1) = -b
    repres(3,2,2) = -a
    !
    repres(4,1,1) = 1.0_ark
    repres(4,1,2) = 0.0_ark
    repres(4,2,1) = 0.0_ark
    repres(4,2,2) = -1.0_ark
    !
    repres(5,1,1) = -a
    repres(5,1,2) = -b
    repres(5,2,1) = -b
    repres(5,2,2) =  a
    !
    repres(6,1,1) = -a
    repres(6,1,2) = b
    repres(6,2,1) = b
    repres(6,2,2) = a
    !
    if (ioper<0.or.ioper>6) then
      write (6,"('symmetry_transformation_dipole: operation ',i8,' unknown')") ioper
      stop 'symmetry_transformation_dipole - bad operation. type'
    endif
    !
    dst = matmul(transpose(repres(ioper,:,:)),src) 
    !
  end subroutine symmetry_transformation_dipole  


  !
end module pot_zxy3

!
!  This unit is for a user defined potential: N2O PES 
!
module pot_user
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xyz

  implicit none

  public MLdipole,MLpoten,ML_MEP

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  !
  integer(ik) :: npol(5),order(5,5)
  real(ark)   :: gam1(5,5),x0(5,5),r0(5,5),beta(5,5)
  integer(ik) :: np(5)
  !
  contains
  !
  function ML_MEP(dim,rho)  result(f)

   integer(ik),intent(in) ::  dim
   real(ark),intent(in)   ::  rho
   real(ark)              ::  f(dim)
   !
   if (dim/=3) stop 'Illegal size of the function - must be 3'
   !
   f(:) = molec%local_eq(:)
   f(molec%Ncoords) = rho

  end function ML_MEP



 !
 ! Defining potential energy function (built for SO2)

 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f,f1,f2
   integer(ik)            ::  nparam
   !
   f1 = MLpoten_n2opotlongrange(ncoords,natoms,local,xyz,force)
   !
   nparam = int(force(1))
   !
   f2 = MLpoten_xyz_N2O_Zobov(ncoords,natoms,local,xyz,force(nparam+1:))
   !
   f = f1+f2
   !
   !
 end function MLpoten
 !


 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
    !
    real(ark) :: xyz_(natoms,3), r1, r2,  n1(3), n2(3), n3(3), tmat(3,3), CN_CM(3),&
                 x(natoms,3),xyz0(natoms,3),cos_theta,alpha,mu(3),u1(3),u2(3),u3(3),bigr,smallr

    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      xyz0 = MLloc2pqr_xyz(local)
      !
      x(1,:) = xyz0(2,:) - xyz0(1,:)
      x(2,:) = xyz0(3,:) - xyz0(1,:)
      !
    else
      !
      x(1,:) = xyz(2,:) - xyz(1,:)
      x(2,:) = xyz(3,:) - xyz(1,:)
      !
    endif
    !
    r1 = sqrt(sum(x(1,:)**2))
    r2 = sqrt(sum(x(2,:)**2))
    !
    stop 'MLdipole N2O is undefined'
    !
  end subroutine MLdipole




  function MLpoten_xyz_N2O_Zobov(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f
   !
   real(ark)            :: r1,r2,alpha,xcos,v,v1,v2,v3,v4,v5,v6,v7,v8

   real(ark)            :: aa1,aa2,re1,re2,alphae,xst,xs1,xs2,v0
   integer(ik)          :: N
   real(ark)            :: rho,xp,vp1
   real(ark)            :: ReNN_ang,ReNO_ang
   !
   if (verbose>=6) write(out,"('MLpoten_xyz_tyuterev/start')")
   !
   r1  = local(1) ;  r2    = local(2) ;  alpha = local(3)
   !
   if (molec%Nangles==0) then
     stop 'MLpoten_xyz_N2O_Zobov: ilegal number of bond anlges 0'
   endif 
   !
   ReNN_ang    = force( 1)
   ReNO_ang    = force( 2)
   alphae      = force( 3)*pi/180.0_ark
   aa1         = force( 4)
   aa2         = force( 5)
   !
   ! calculate potential energy function values
   !    
   !
   !
   !
   rho=pi-alpha
   !   
   xst= 1.0_ark - cos(rho)
   xs1= 1.0_ark - exp(-aa1*(r1-ReNN_ang))
   xs2= 1.0_ark - exp(-aa2*(r2-ReNO_ang))
   !
   N = 5
   !
  v0= force(N+1)  *xs1**0*xs2**0*xst**0
 vp1= force(N+2)  *xs1**1*xs2**0*xst**0& 
     +force(N+3)  *xs1**0*xs2**1*xst**0& 
     +force(N+4)  *xs1**0*xs2**0*xst**1& 
     +force(N+5)  *xs1**2*xs2**0*xst**0& 
     +force(N+6)  *xs1**1*xs2**1*xst**0& 
     +force(N+7)  *xs1**1*xs2**0*xst**1& 
     +force(N+8)  *xs1**0*xs2**2*xst**0& 
     +force(N+9)  *xs1**0*xs2**1*xst**1& 
     +force(N+10) *xs1**0*xs2**0*xst**2&    ! end of 2
     +force(N+11) *xs1**3*xs2**0*xst**0& 
     +force(N+12) *xs1**2*xs2**1*xst**0& 
     +force(N+13) *xs1**2*xs2**0*xst**1& 
     +force(N+14) *xs1**1*xs2**2*xst**0& 
     +force(N+15) *xs1**1*xs2**1*xst**1& 
     +force(N+16) *xs1**1*xs2**0*xst**2& 
     +force(N+17) *xs1**0*xs2**3*xst**0& 
     +force(N+18) *xs1**0*xs2**2*xst**1& 
     +force(N+19) *xs1**0*xs2**1*xst**2& 
     +force(N+20) *xs1**0*xs2**0*xst**3     !  end of 3

       f=v0+vp1

  end function MLpoten_xyz_N2O_Zobov


  function MLpoten_n2opotlongrange(ncoords,natoms,local,xyz,force) result(f)

     integer(ik),intent(in) ::  ncoords,natoms
     real(ark),intent(in)   ::  local(ncoords)
     real(ark),intent(in)   ::  xyz(natoms,3)
     real(ark),intent(in)   ::  force( :)
     real(ark)              ::  f
     integer(ik) :: iopttmp(3)

     real(ark) :: coe,emin
     real(ark) :: icoe

     real(ark) :: De1,De_1,De2,De_2,De3,De_3
     real(ark) :: edp1,edp2,edp3,edp4,edp5,edp6
     real(ark) :: rnnref,rnoref,alpha1,alpha1b,alpha2,alpha2b
     !
     integer(ik) :: i,k,Ncoe
     real(ark)   :: r12ref,Ae1,Ae2,edamp2,edamp4,edamp5,edamp6,alphae
     !
     real(ark)   :: rmin,rminbohr,alpha,rref,a2b
     !
     real(ark)   :: str1,str2,enetmp1,atp1,enetmp2,edamp11,edamp12,edamp1,V
     real(ark)   :: v0,rrco1,rrco2,r1,r2,a3,xx1,xx2,rco1,rco2,ang,dstr1,dstr2,sumstr2,sumstr4,angref,angx,bdamp2,bdamp4
     !
     real(ark)   :: v1(3),v2(3),x1,x2,cosalpha1,rnn,rno,ang2,sumx0,rno2,enetmp2B
     real(ark)   :: etmp2,anno,ax,enetmp3,etmp1
        !
        !
        rnn  = local(1) ;  rno    = local(2) ;  alpha = local(3)
        !
        alphae  = molec%alphaeq(1)
        !
        if (molec%Nangles==0) then
          stop 'MLpoten_xyz_N2O_Zobov: ilegal number of bond anlges 0'
        endif 
	    !
        Ncoe    = int(force(1))
        !
        rnnref  = force( 2)
        rnoref  = force( 3)
        alphae  = force( 4)
        alpha1  = force( 5)
        alpha2  = force( 6)
        De1     = force( 7) 
        De_1    = force( 8) 
        De2     = force( 9) 
        De_2    = force(10) 
        De3     = force(11) 
        De_3    = force(12) 
        edamp2  = force(13)
        edamp4  = force(14)
        edamp5  = force(15)
        edamp6  = force(16)
        edp1    = force(17)
        edp2    = force(18)
        edp3    = force(19)
        edp4    = force(20)
        edp5    = force(21)
        edp6    = force(22)
        !
        Emin    = force(23)
        !
        r1=rnn-rnnref
        r2=rno-rnoref
        !
        a3=sin(alpha)
        !
	    a3=a3*a3
        !
        v0=0
        do i=22,Ncoe
          v0=v0+force(i)*r1**molec%pot_ind(1,i)*r2**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
        end do
        !
        sumstr2=(rnn-rnnref)**2+(rno-rnoref)**2
        sumstr4=(rnn-rnnref)**4+(rno-rnoref)**4
        ang2=pi-alpha
        angref=150.0_ark*pi/180.0_ark-pi
        !angx=min(-abs(alpha*180.0_ark/pi)-angref,0.0_ark)
        !bdamp2=angx**2;bdamp4=angx**4        !
        bdamp2=ang2**2; bdamp4=ang2**4
        etmp2=exp(edp1*sumstr2+edp2*sumstr4+edp3*bdamp2+edp4*bdamp4)

	    sumx0=V0
	    !
        V0=V0*etmp2
        !
! calc etmp1
        sumstr2=(rnn-rnnref)**2+(rno-rnoref)**2
        sumstr4=(rnn-rnnref)**4+(rno-rnoref)**4
        enetmp1=De1*(1.0_ark-exp(-alpha1*(rnn-rnnref)))**2 + De_1*(1.0_ark-exp(-alpha1*(rnn-rnnref)))**4
        enetmp2=De2*(1.0_ark-exp(-alpha2*(rno-rnoref)))**2 + De_2*(1.0_ark-exp(-alpha2*(rno-rnoref)))**4
        !
        rno2=rnn*rnn+rno*rno-2.0_ark*rnn*rno*cos(alpha)
        rno2=sqrt(rno2)
        enetmp2B=De2*(1-exp(-alpha2*(rno2-rnoref)))**2 + De_2*(1.0_ark-exp(-alpha2*(rno2-rnoref)))**4
        anno=alpha; ax=(pi-anno)/2.0_ark
        enetmp3=De3*sin(ax)**2 + De_3*sin(ax)**4  !bending simulation
        edamp1=exp(edp5*sumstr2+0.0_ark*edp6*sumstr4)
        enetmp3=edamp1*enetmp3
        !
        etmp1=enetmp1+enetmp2+enetmp3 !+ enetmp2B
        !
        V0=V0+etmp1
	    !
        f=V0-emin*219474.63067d0
        !
      end function MLpoten_n2opotlongrange





end module pot_user

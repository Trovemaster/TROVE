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
   integer(ik)            ::  nparam,ntot
   !
   !
   nparam = int(force(1))
   ntot = size(force)
   !
   f = MLpoten_co2_ames1_morse_powers_isotopologue(ncoords,natoms,local,xyz,force)
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





  function MLpoten_co2_ames1_morse_powers_isotopologue(ncoords,natoms,local,xyz,force) result(f)
   !
   implicit none
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          :: i,k,Ncoe
   real(ark)            :: r12ref,alpha2,De1,De2,Ae1,Ae2,edamp2,edamp4,edamp5,edamp6
   !
   real(ark)            :: emin,rmin,rminbohr,alpha,rref,a2b
   !
   real(ark)            :: str1,str2,enetmp1,atp1,enetmp2,edamp11,edamp12,edamp1,V
   real(ark)            :: v0,rrco1,rrco2,r1,r2,a3,xx1,xx2,rco1,rco2,ang,dstr1,dstr2,sumstr2,sumstr4,angref,angx,bdamp2,bdamp4
   !
   real(ark)            ::   v1(3),v2(3),x1,x2,cosalpha1
   integer(ik)          :: Ncoeff2,n1,n2,n3
   real(ark)            :: r12,r32,y1,y2,y3,aa1,re12,alphae,carbon12Mass,carbonIsotopeMass,carbonIsotopologueFactor,f1,f2
   
   character(len=cl)    :: txt = 'MLpoten_c3_R_theta'
        !
        Ncoe = force(molec%parmax)
        !
        r12ref  = force(1)
        alpha2  = force(2)
        De1     = force(3)
        De2     = force(4)
        Ae1     = force(5)
        Ae2     = force(6)
        edamp2  = force(7)
        edamp4  = force(8)
        edamp5  = force(9)
        edamp6  = force(10)
        Emin    = force(11)
        Rmin    = force(12)
        rminbohr= force(13)
        alpha   = force(14)
        rref    = force(15)
        !
        rco1 = local(1)
        rco2 = local(2)
        !
        rrco1=rco1/0.529177249_ark
        rrco2=rco2/0.529177249_ark
        !
        ang = pi-local(3)
        if  (molec%Ndihedrals>1) then
          !
          ang = asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
          !
        endif 
        !
        if  (molec%Ndihedrals>1) then
          !
          ang = asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
          !
          v1(:) = xyz(2,:)-xyz(1,:)
          v2(:) = xyz(3,:)-xyz(1,:)
          !
          x1 = sqrt(sum(v1(:)**2))
          x2 = sqrt(sum(v2(:)**2))
          !
          cosalpha1 = sum( v1(:)*v2(:))/(x1*x2)
          !
          !alpha = aacos(cosalpha1,txt)
          ang  = aacos(cosalpha1,txt)-pi
          !
        endif
        !
        r1=1._ark-exp(-alpha*(rco1-r12ref))
        r2=1._ark-exp(-alpha*(rco2-r12ref))
        a3=cos(ang)
        !
        v0=0
        do i=16,Ncoe
           !
           v0=v0+force(i)*r1**molec%pot_ind(1,i)*r2**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
           if(molec%pot_ind(1,i).ne.molec%pot_ind(2,i))then
             v0=v0+force(i)*r2**molec%pot_ind(1,i)*r1**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
           end if
           !
        end do
        !
        xx1=rco1; xx2=rco2
        dstr1=(xx1-r12ref)
        dstr2=(xx2-r12ref)
        sumstr2=dstr1**2+dstr2**2
        sumstr4=dstr1**4+dstr2**4
        !
        angref=150.0_ark*pi/180.0_ark-pi
        !
        angx=min(-abs(ang)-angref,0.0_ark)
        bdamp2=angx**2;bdamp4=angx**4
        bdamp2=ang**2; bdamp4=ang**4
        !
        !angx=min(-abs(ang)+angref,0.d0)
        !bdamp2=angx**2;bdamp4=angx**4
        !bdamp2=ang**2; bdamp4=ang**4
        !
        V0=V0*exp(edamp2*sumstr2+edamp4*sumstr4+edamp5*bdamp2+edamp6*bdamp4)
        !
        a2b=0.529177249_ark
        str1=1._ark-exp(-alpha2*(xx1-r12ref))
        str2=1._ark-exp(-alpha2*(xx2-r12ref))
        enetmp1=De1*(str1**2+str2**2)+De2*(str1**4+str2**4)
        enetmp1=enetmp1/219474.63067_ark
        !
        atp1=acos(-1._ark)-abs(ang)
        enetmp2=Ae1*(1._ark+cos(atp1))**2+Ae2*(1._ark+cos(atp1))**4
        enetmp2=enetmp2/219474.63067_ark
        !
        sumstr2=(xx1-r12ref)**2+(xx2-r12ref)**2
        sumstr4=(xx1-r12ref)**4+(xx2-r12ref)**4
        edamp11=-0.5_ark; edamp12=-0.5_ark
        edamp1=exp(0.2_ark*edamp11*sumstr2+0*edamp12*sumstr4)
        enetmp2=edamp1*enetmp2
        !
        V0=V0+enetmp1+enetmp2
        !
        V=V0-emin
        !
        f1 = V*219474.63067_ark
        !
        ! part 2: isotopologues
        !
        r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
        !
        re12    = force(Ncoe+1)
        alphae  = force(Ncoe+2)*pi/180.0_ark
        aa1     = force(Ncoe+3)
        carbon12Mass = force(Ncoe+4)
        carbonIsotopeMass = force(Ncoe+5)
        !
        ! calculate potential energy function values
        !
        y1=1.0_ark-exp(-aa1*(r12-re12))
        y2=1.0_ark-exp(-aa1*(r32-re12))
        !
        y3=cos(alpha)-cos(alphae)
        !
        carbonIsotopologueFactor = (carbonIsotopeMass - carbon12Mass)/carbon12Mass
        !
        Ncoeff2 = molec%parmax-1
        !
        v0=0
        do i=Ncoe+6,Ncoeff2
           !
           n1 = molec%pot_ind(1,i)
           n2 = molec%pot_ind(2,i)
           n3 = molec%pot_ind(3,i)
           !
           ! write(out, '(A, I0)') "The value of i is:", i
           ! write(out, '(A, f18.8)') "The value of force(i) is:", force(i)
           v0=v0+force(i)*y1**n1*y2**n2*y3**n3
           if( n1/=n2 )then
             v0=v0+force(i)*y2**n1*y1**n2*y3**n3
           end if
           !
        end do
        !
        f2 = v0*carbonIsotopologueFactor
        ! write(out, '(A, f18.8)') "The isotopologue correction to the potential is:", v0 
        f=f1+f2
        !
  end function MLpoten_co2_ames1_morse_powers_isotopologue


end module pot_user

!
!  This unit defines all specific routines for a fiveatomic molecule of XY4 type
!
module pot_xy4
  use accuracy
  use moltype
  use lapack
  use symmetry,     only : sym

  implicit none

  public MLpoten_xy4_ZZZ,MLpoten_xy4_Bowman2000,ML_XY4_calc_alpha34,MLpoten_xy4_Nikitin
  !
  private
  !
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  !
  contains



    function ML_XY4_calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24) result (alpha34)
      !
      real(ark),intent(in)  ::alpha12,alpha13,alpha14,alpha23,alpha24
      real(ark) alpha34,ca12,ca13,ca23,ca14,ca24
      real(ark) cosalpha

          ca12 = cos(alpha12)
          ca13 = cos(alpha13)
          ca14 = cos(alpha14)
          ca23 = cos(alpha23)
          ca24 = cos(alpha24)

          cosalpha = 0.5_ark*( &
                     2.0_ark*ca12*ca23*ca14-     &
                     2.0_ark*ca13*ca14+     &
                     2.0_ark*ca12*ca24*ca13-     &
                     2.0_ark*ca23*ca24+     &
                     2.0_ark*sqrt(     &
                     1.0_ark-ca24**2+     &
                     2.0_ark*ca12*ca24*ca14+     &
                     2.0_ark*ca12*ca23*ca13-     &
                     2.0_ark*ca12**3*ca23*ca13-     &
                     2.0_ark*ca12**2-     &
                     ca14**2-     &
                     ca13**2+     &
                     ca14**2*ca23**2+     &
                     ca13**2*ca24**2+     &
                     4.0_ark*ca12**2*ca23*ca14*ca24*ca13+     &
                     ca12**4-     &
                     2.0_ark*ca12*ca23*ca14**2*ca13-     &
                     2.0_ark*ca12*ca23**2*ca14*ca24-     &
                     2.0_ark*ca13**2*ca14*ca12*ca24-     &
                     2.0_ark*ca12*ca24**2*ca13*ca23+     &
                     ca13**2*ca14**2+     &
                     ca23**2*ca24**2+     &
                     ca12**2*ca13**2+     &
                     ca12**2*ca14**2+     &
                     ca12**2*ca23**2+     &
                     ca12**2*ca24**2-     &
                     2.0_ark*ca12**3*ca24*ca14-     &
                     ca23**2))/( -sin(alpha12)**2 )

          if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
             !
             write (out,"('ML_XY4_calc_alpha34: cosalpha>1: ',f18.8)") cosalpha
             stop 'ML_XY4_calc_alpha34 - bad cosalpha'
             !
          elseif ( cosalpha>=1.0_ark) then 
             alpha34 = 0.0_ark
          else 
             alpha34 = acos(cosalpha)
          endif
          !
          !alpha34 =molec%local_eq(5) + ( 5.0_ark*molec%local_eq(5) - (alpha12+alpha13+alpha23+alpha14+alpha24) )
          !

    end function ML_XY4_calc_alpha34



!
! CH4 PES from Bowman: 
! S. Carter and J. M. Bowman, J. Phys. Chem. A 104, 2355 (2000).
! 

  function MLpoten_xy4_Bowman2000(ncoords,natoms,local,xyz,force) result(f) 
    !
    integer(ik),intent(in) ::  ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords)
    real(ark),intent(in)   ::  xyz(natoms,3)
    real(ark),intent(in)   ::  force(:)
    real(ark)              ::  f
    !
    real(ark) :: s(9),f2(9,9),f3(9,9,9),f4(9,9,9,9),fact(5)
    real(ark) :: dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4
    real(ark) :: r1,r2,r3,r4,cosa12,a12,cosa13,a13,cosa23,a23,cosa14,a14,cosa24,a24,cosa34,a34
    real(ark) :: dr1,dr2,dr3,dr4,da12,da13,da23,da14,da24,da34
    real(ark) :: s1,s2a,s2b,s3x,s3y,s3z,s4x,s4y,s4z
    real(ark) :: ff2,ff3,ff4,srt2,srt12
    integer(ik) :: ind(9),nn,i,j,k,l,m,n,n2,n3,n4
    real(ark),parameter ::  re = 1.0890_ark, ae = 1.91063324_ark
     !
     srt2=1.0_ark/sqrt(2.0_ark)
     srt12=1.0_ark/sqrt(12.0_ark)
     !
     fact(1)=1.0_ark
     fact(2)=1.0_ark
     fact(3)=2.0_ark
     fact(4)=6.0_ark
     fact(5)=24.0_ark
     !
     f2 = 0
     f3 = 0 
     f4 = 0 
     !
     n2 = 12
     !
     do i=1,n2
       !
       ind = 1
       !
       m = molec%pot_ind(1,i)
       n = molec%pot_ind(2,i)
       !
       f2(m,n) =  force(i)
       !
       ind(m)=ind(m)+1
       ind(n)=ind(n)+1
       f2(m,n)=f2(m,n)/( fact(ind(1))*fact(ind(2))*fact(ind(3))*fact(ind(4))*&
                         fact(ind(5))*fact(ind(6))*fact(ind(7))*fact(ind(8))*fact(ind(9)) )
       !
     enddo
     !
     n3 = n2+37
     !
     do i=n2+1,n3
       !
       ind = 1
       !
       j = molec%pot_ind(1,i)
       m = molec%pot_ind(2,i)
       n = molec%pot_ind(3,i)
       !
       f3(j,m,n) =  force(i)
       !
       ind(j)=ind(j)+1
       ind(m)=ind(m)+1
       ind(n)=ind(n)+1
       f3(j,m,n)=f3(j,m,n)/(fact(ind(1))*fact(ind(2))*fact(ind(3))*fact(ind(4))*&
                            fact(ind(5))*fact(ind(6))*fact(ind(7))*fact(ind(8))*fact(ind(9)))
     enddo
     !
     n4 = n3 + 116
     !
     do i=n3+1,n4
       !
       ind = 1
       !
       k = molec%pot_ind(1,i)
       j = molec%pot_ind(2,i)
       m = molec%pot_ind(3,i)
       n = molec%pot_ind(4,i)
       !
       f4(k,j,m,n) =  force(i)
       !
       ind(k)=ind(k)+1
       ind(j)=ind(j)+1
       ind(m)=ind(m)+1
       ind(n)=ind(n)+1
       !
       f4(k,j,m,n)=f4(k,j,m,n)/(fact(ind(1))*fact(ind(2))*fact(ind(3))*fact(ind(4))*&
                   fact(ind(5))*fact(ind(6))*fact(ind(7))*fact(ind(8))*fact(ind(9)))
     enddo
     !
     !  set up for ch4
     !  c is atom 1 and then come the four h atoms
     !  get the four ch distances
     !
     dx1=xyz(2,1)-xyz(1,1)
     dy1=xyz(2,2)-xyz(1,2)
     dz1=xyz(2,3)-xyz(1,3)
     dx2=xyz(3,1)-xyz(1,1)
     dy2=xyz(3,2)-xyz(1,2)
     dz2=xyz(3,3)-xyz(1,3)
     dx3=xyz(4,1)-xyz(1,1)
     dy3=xyz(4,2)-xyz(1,2)
     dz3=xyz(4,3)-xyz(1,3)
     dx4=xyz(5,1)-xyz(1,1)
     dy4=xyz(5,2)-xyz(1,2)
     dz4=xyz(5,3)-xyz(1,3)
     !
     r1=sqrt(dx1**2+dy1**2+dz1**2)
     r2=sqrt(dx2**2+dy2**2+dz2**2)
     r3=sqrt(dx3**2+dy3**2+dz3**2)
     r4=sqrt(dx4**2+dy4**2+dz4**2)
     !
     !  now get the angles
     !
     cosa12=(dx1*dx2+dy1*dy2+dz1*dz2)/(r1*r2)
     a12=acos(cosa12)
     cosa13=(dx1*dx3+dy1*dy3+dz1*dz3)/(r1*r3)
     a13=acos(cosa13)
     cosa14=(dx1*dx4+dy1*dy4+dz1*dz4)/(r1*r4)
     a14=acos(cosa14)
     cosa23=(dx2*dx3+dy2*dy3+dz2*dz3)/(r2*r3)
     a23=acos(cosa23)
     cosa24=(dx2*dx4+dy2*dy4+dz2*dz4)/(r2*r4)
     a24=acos(cosa24)
     cosa34=(dx3*dx4+dy3*dy4+dz3*dz4)/(r3*r4)
     a34=acos(cosa34)
     !
     ! convert to angst
     !
     !r1=r1*bohr
     !r2=r2*bohr
     !r3=r3*bohr
     !r4=r4*bohr
     !
     !  get deltas
     !
     dr1=r1-re
     dr2=r2-re
     dr3=r3-re
     dr4=r4-re
     da13=a13-ae
     da12=a12-ae
     da14=a14-ae
     da23=a23-ae
     da24=a24-ae
     da34=a34-ae
     !
     !  now form the symmetry displacement coordinates
     s1=(dr1+dr2+dr3+dr4)/2.0_ark
     s(1)=s1
     s2a=(2.0_ark*da12-da13-da14-da23-da24+2.0_ark*da34)*srt12
     s(2)=s2a
     s2b=(da13-da14-da23+da24)/2.0_ark
     s(3)=s2b
     s3x=(dr1-dr2+dr3-dr4)/2.0_ark
     s(4)=s3x
     s3y=(dr1-dr2-dr3+dr4)/2.0_ark
     s(5)=s3y
     s3z=(dr1+dr2-dr3-dr4)/2.0_ark
     s(6)=s3z
     s4x=(da24-da13)*srt2
     s(7)=s4x
     s4y=(da23-da14)*srt2
     s(8)=s4y
     s4z=(da34-da12)*srt2
     s(9)=s4z
     !
     !  quadratic ff first
     !
     ff2=0
     do i=1,9
       do j=1,i
        !
        !  take account of index symmetry
        !
        ff2=ff2+f2(i,j)*s(i)*s(j)
       enddo
     enddo
     !
     ! now do cubic
     ff3=0
     do i=1,9
       do j=1,i
         do k=1,j
           ff3=ff3+f3(i,j,k)*s(i)*s(j)*s(k)
         enddo
       enddo
     enddo
     !  now do quartic
     ff4=0
     do i=1,9
       do j=1,i
         do k=1,j
           do l=1,k
             ff4=ff4+f4(i,j,k,l)*s(i)*s(j)*s(k)*s(l)
           enddo
         enddo
       enddo
     enddo
     !
     !  convert to hartrees
     !
     f=(ff2+ff3+ff4)*1.0e-11/planck/vellgt
     !
  end function MLpoten_xy4_Bowman2000
  !


  function MLpoten_xy4_Nikitin(ncoords,natoms,local,xyz,force) result(f) 
    !
    integer(ik),intent(in) ::  ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords)
    real(ark),intent(in)   ::  xyz(natoms,3)
    real(ark),intent(in)   ::  force(:)
    real(ark)              ::  f
    !
    real(ark) :: r1,r2,r3,r4,cosa12,alpha12,cosa13,alpha13,cosa23,alpha23,cosa14,alpha14,cosa24,alpha24,cosa34    !
    real(ark) :: dx(4,3),g(40612),alpha34,beta312,beta412,cosbeta
    !
    integer(ik) :: k_ind(14),iterm,n
    !
    !real(ark) ::   betac, betaa,rc2,r02,a12,a13,a14,a23,a24,a34,t12,t13,t14,t23,t24,t34,vdump,cosae
    !
    real(ark) :: v
    real(ark) :: re,alphae,a,y(14),xi(14)
    !
    character(len=cl)         :: txt
    !
    real(ark) :: ae= 1.9106332362490185563277142050315_ark
      !
      if (size(local)==9.and.molec%NDihedrals==0) then
        !
        r1  = local(1) 
        r2  = local(2) 
        r3  = local(3) 
        r4  = local(4) 
        !
        alpha12 = local(5)
        alpha13 = local(6)
        alpha23 = local(7)
        alpha14 = local(8)
        alpha24 = local(9)
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
      else
        !
        !  set up for ch4
        !  c is atom 1 and then come the four h atoms
        !  get the four ch distances
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
        cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
        beta312 = aacos(cosbeta,txt)
        !
        cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
        beta412 = aacos(cosbeta,txt)
        !
        !cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
        !alpha34 = aacos(cosa34,txt)
        !
      endif
      !
      re = molec%req(1)
      !
      alphae  = ae
      !
      a = molec%specparam(1)
      !
      y(1)=1.0_ark-exp(-a*(r1-re))
      y(2)=1.0_ark-exp(-a*(r2-re))
      y(3)=1.0_ark-exp(-a*(r3-re))
      y(4)=1.0_ark-exp(-a*(r4-re))
      !
      y(5)=cosa12-cos(ae)
      y(6)=cosa13-cos(ae)
      y(7)=cosa14-cos(ae)
      !
      y(8)=cos(beta312)
      y(9)=cos(beta412)
      !
      y(10)=sin(alpha12)-sin(ae)
      y(11)=sin(alpha13)-sin(ae)
      y(12)=sin(alpha14)-sin(ae)
      !
      y(13)=sin(beta312)
      y(14)=sin(beta412)
      !
      f = 0
      !
      n = size(force)
      !
      do iterm=1,n
        !
        k_ind(1:4) = molec%pot_ind(1:4,iterm)
        k_ind(5:9) = molec%pot_ind(5:9,iterm)/10
        k_ind(10:14) = mod(molec%pot_ind(5:9,iterm),10)
        !
        xi(1:14) = y(1:14)**k_ind(1:14)
        !
        f = f + g(iterm)*product(xi(1:14))
        !
      enddo
      !
      f = f*219474.6313705_rk
      !
  end function MLpoten_xy4_Nikitin
  !

  !
  function MLpoten_xy4_ZZZ(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)    ::  v0,v2,v3,v4,chi(10)

   real(ark)    :: s1,s2a,s2b,s3x,s3y,s3z,s4x,s4y,s4z
   real(ark)    :: alpha12,alpha13,alpha14,alpha23,alpha24,alpha34,beta312,beta412,cosa23,cosa24,cosa34
   
   real(ark)    ::F0,F11,F22,F33,F43,F44,F111,F12a2a,F13x3x,F13x4x,    &
                F14x4x,F2a2a2a,F2a3z3z,F2a3z4z,F2a4z4z,F3x3y3z,        &
                F3x3y4z,F3x4y4z,F4x4y4z,F1111,F112a2a,F113x4x,         &
                F114x4x,F12a2a2a,F12a3z3z,F12a3z4z,F12a4z4z,           &
                F13x3y3z,F13x3y4z,F13x4y4z,F14x4y4z,F2a2a2a2a,         &
                F2a2a3z3z,F2b2b3z3z,F2a2a3z4z,F2b2b3z4z,F2a2a4z4z,     &
                F2b2b4z4z,F2a3x3y4z,F2a3z4x4y,F3x3x3x3x,F3x3x3y3y,     &
                F3x3x3x4x,F3x3x3y4y,F3x3x4x4x,F3x3y4x4y,F3x3x4y4y,     &
                F3x4x4x4x,F3x4x4y4y,F4x4x4x4x,F4x4x4y4y,F113x3x
                !
   character(len=cl)         :: txt

    if (verbose>=6) write(out,"(/'MLpoten_xy4_zzz/start')") 
    !
    txt = 'MLpoten_xy4_ZZZ'
    !
    F0          =  force(  1)
    F11         =  force(  2)
    F22         =  force(  3)
    F33         =  force(  4)
    F43         =  force(  5)
    F44         =  force(  6)
    
    F111        =  force(  7)
    F12a2a      =  force(  8)
    F13x3x      =  force(  9)
    F13x4x      =  force( 10)
    F14x4x      =  force( 11)
    F2a2a2a     =  force( 12)
    F2a3z3z     =  force( 13)
    F2a3z4z     =  force( 14)
    F2a4z4z     =  force( 15)
    F3x3y3z     =  force( 16)
    F3x3y4z     =  force( 17)
    F3x4y4z     =  force( 18)
    F4x4y4z     =  force( 19)

    F1111       =  force( 20)
    F112a2a     =  force( 21)
    F113x3x     =  force( 22)
    F113x4x     =  force( 23)
    F114x4x     =  force( 24)
    F12a2a2a    =  force( 25)
    F12a3z3z    =  force( 26)
    F12a3z4z    =  force( 27)
    F12a4z4z    =  force( 28)
    F13x3y3z    =  force( 29)
    F13x3y4z    =  force( 30)
    F13x4y4z    =  force( 31)
    F14x4y4z    =  force( 32)
    F2a2a2a2a   =  force( 33)
    F2a2a3z3z   =  force( 34)
    F2b2b3z3z   =  force( 35)
    F2a2a3z4z   =  force( 36)
    F2b2b3z4z   =  force( 37)
    F2a2a4z4z   =  force( 38)
    F2b2b4z4z   =  force( 39)
    F2a3x3y4z   =  force( 40)
    F2a3z4x4y   =  force( 41)
    F3x3x3x3x   =  force( 42)
    F3x3x3y3y   =  force( 43)
    F3x3x3x4x   =  force( 44)
    F3x3x3y4y   =  force( 45)
    F3x3x4x4x   =  force( 46)
    F3x3y4x4y   =  force( 47)
    F3x3x4y4y   =  force( 48)
    F3x4x4x4x   =  force( 49)
    F3x4x4y4y   =  force( 50)
    F4x4x4x4x   =  force( 51)
    F4x4x4y4y   =  force( 52)

    !
    select case(trim(molec%coords_transform))
       !
       case default
       write (out,"('MLpoten_xy4_ZZZ: coord_transf ',a,' unknown')") trim(molec%coords_transform)
       stop 'MLpoten_xy4_ZZZ - bad coord. type'
       !
    case('R-ALPHA')
       !
       chi(1:10) = local(1:10) - molec%local_eq(1:10)
       !
       s1  =(chi(1)+chi(2)+chi(3)+chi(4))*0.5_ark
       s3x =(chi(1)-chi(2)+chi(3)-chi(4))*0.5_ark
       s3y =(chi(1)-chi(2)-chi(3)+chi(4))*0.5_ark
       s3z =(chi(1)+chi(2)-chi(3)-chi(4))*0.5_ark
       !
       alpha12 = local(5)
       alpha13 = local(6)
       alpha23 = local(7)
       alpha14 = local(8)
       alpha24 = local(9)
       alpha34 = local(10)
       !
       s2a=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
       s2b=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
       s4x=(alpha24-alpha13)/sqrt(2.0_ark)
       s4y=(alpha23-alpha14)/sqrt(2.0_ark)
       s4z=(alpha34-alpha12)/sqrt(2.0_ark)
       !
    case('R-SYM','R-SYM-F-E')
       !
       chi(1:9) = local(1:9) - molec%local_eq(1:9)
       !
       !alpha=calc_alpha34(local(5),local(6),local(7),local(8),local(9))
       !alpha=find_alpha34(local(5),local(6),local(7),local(8),local(9))
       !
       !alpha=alpha-molec%local_eq(5)
       !
       s1  =(chi(1)+chi(2)+chi(3)+chi(4))*0.5_ark
       !s2a =(2.0_ark*chi(5)-chi(6)-chi(8)-chi(7)-chi(9)+2.0_ark*alpha)/sqrt(12.0_ark)
       !s2b =(chi(6)-chi(8)-chi(7)+chi(9))*0.5_ark
       s3x =(chi(1)-chi(2)+chi(3)-chi(4))*0.5_ark
       s3y =(chi(1)-chi(2)-chi(3)+chi(4))*0.5_ark
       s3z =(chi(1)+chi(2)-chi(3)-chi(4))*0.5_ark
       !
       !s4x =(chi(9)-chi(6))/sqrt(2.0_ark)
       !s4y =(chi(7)-chi(8))/sqrt(2.0_ark)
       !s4z =(alpha-chi(5))/sqrt(2.0_ark)
       !
       alpha12 = local(5)
       alpha13 = local(6)
       alpha23 = local(7)
       alpha14 = local(8)
       alpha24 = local(9)
       !
       !alpha34 = calc_alpha34(alpha12,alpha13,alpha23,alpha14,alpha24)
       !
       if (size(local)==10) then 
         !
         alpha34 = local(10)
         !
        else
         !
         alpha34 = ML_XY4_calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
         !
       endif 
       !
       !alpha34 = find_alpha34(local(5),local(6),local(7),local(8),local(9))
       !
       s2a=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
       s2b=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
       s4x=(alpha24-alpha13)/sqrt(2.0_ark)
       s4y=(alpha23-alpha14)/sqrt(2.0_ark)
       s4z=(alpha34-alpha12)/sqrt(2.0_ark)
       !
    case('R-2-SYM')
       !
       chi(1:9) = local(1:9) - molec%local_eq(1:9)
       !
       alpha12 = local(5)
       alpha13 = local(6)
       alpha14 = local(7)
       beta312 = local(8)
       beta412 = local(9)
       !
       cosa23 = cos(alpha12)*cos(alpha13)+cos(beta312)*sin(alpha12)*sin(alpha13)
       cosa24 = cos(alpha12)*cos(alpha14)+cos(beta412)*sin(alpha12)*sin(alpha14)
       cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
       !
       alpha23 = aacos(cosa23,txt)
       alpha24 = aacos(cosa24,txt)
       alpha34 = aacos(cosa34,txt)
       !
       !alpha34=find_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
       !alpha34=calc_alpha34(alpha12,alpha13,alpha14,alpha23,alpha24)
       !
       s1  =(chi(1)+chi(2)+chi(3)+chi(4))*0.5_ark
       !
       s3x =(chi(1)-chi(2)+chi(3)-chi(4))*0.5_ark
       s3y =(chi(1)-chi(2)-chi(3)+chi(4))*0.5_ark
       s3z =(chi(1)+chi(2)-chi(3)-chi(4))*0.5_ark
       !
       s2a=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
       s2b=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
       s4x=(alpha24-alpha13)/sqrt(2.0_ark)
       s4y=(alpha23-alpha14)/sqrt(2.0_ark)
       s4z=(alpha34-alpha12)/sqrt(2.0_ark)
       !
    case('SYM-S')
       !
       CHI(1)=(local(1)+local(2)+local(3)+local(4))*0.5_ark-2.0_ark*molec%local_eq(1)
       CHI(2)=(2.0_ark*local(5)-local(6)-local(7)-local(8)-local(9)+2.0_ark*local(10))/sqrt(12.0_ark)
       CHI(3)=(local(6)-local(7)-local(8)+local(9))*0.5_ark
       CHI(4)=(local(1)-local(2)+local(3)-local(4))*0.5_ark
       CHI(5)=(local(1)-local(2)-local(3)+local(4))*0.5_ark
       CHI(6)=(local(1)+local(2)-local(3)-local(4))*0.5_ark
       CHI(7)=(local(9)-local(6))/sqrt(2.0_ark)
       CHI(8)=(local(8)-local(7))/sqrt(2.0_ark)
       CHI(9)=(local(10)-local(5))/sqrt(2.0_ark)
       !
       s1 =  CHI(1) ;  s2a    = CHI(2) ;  s2b    = CHI(3)
       s3x = CHI(4) ;  s3y    = CHI(5) ;  s3z    = CHI(6)
       s4x = CHI(7) ;  s4y    = CHI(8) ;  s4z    = CHI(9)
       !
    end select
    !
    v0 = F0

    v2 =(F11*s1**2 &
          +F22*(s2a**2+s2b**2) & 
          +F33*(s3x**2+s3y**2+s3z**2)  &
          +2.0_ark*F43*(s3x*s4x+s3y*s4y+s3z*s4z)        &
          +F44*(s4x**2+s4y**2+s4z**2)    &
          )/2.0_ark



    v3 = (F111*s1**3                                                        &
           +3.0_ark*F12a2a*s1*(s2a**2+s2b**2)                                &
           +3.0_ark*F13x3x*s1*(s3x**2+s3y**2+s3z**2)                         &
           +6.0_ark*F13x4x*s1*(s3x*s4x+s3y*s4y+s3z*s4z)                      &
           +3.0_ark*F14x4x*s1*(s4x**2+s4y**2+s4z**2)                         &
           +       F2a2a2a*S2a*(s2a**2-3.0_ark*s2b**2)                       &
           +3.0_ark*F2a3z3z*(s2a*(s3z**2-0.5_ark*s3x**2-0.5_ark*s3y**2)        &
                           +s2b*sqrt(3.0_ark)*0.5_ark*(s3x**2-s3y**2))        &
           +6.0_ark*F2a3z4z*(s2a*(s3z*s4z-0.5_ark*s3x*s4x-0.5_ark*s3y*s4y)     &
                           +s2b*sqrt(3.0_ark)*0.5_ark*(s3x*s4x-s3y*s4y))      &
           +3.0_ark*F2a4z4z*(s2a*(s4z**2-0.5_ark*s4x**2-0.5_ark*s4y**2)        &
                    +0.5_ark*sqrt(3.0_ark)*s2b*(s4x**2-s4y**2)   )            &
           +6.0_ark*F3x3y3z*s3x*s3y*s3z                                      &
           +6.0_ark*F3x3y4z*(s3x*s3y*s4z+s3x*s3z*s4y+s3y*s3z*s4x)            &
           +6.0_ark*F3x4y4z*(s3x*s4y*s4z+s3y*s4x*s4z+s3z*s4x*s4y)            &
           +6.0_ark*F4x4y4z*s4x*s4y*s4z   &
           )/6.0_ark

    v4 = (F1111*s1**4 & 
           +6.0_ark*F112a2a*s1**2*(s2a**2+s2b**2)                &
           +6.0_ark*F113x3x*s1**2*(s3x**2+s3y**2+s3z**2)                     &
           +12.0_ark*F113x4x*s1**2*(s3x*s4x+s3y*s4y+s3z*s4z)                 &
           +6.0_ark*F114x4x*s1**2*(s4x**2+s4y**2+s4z**2)                     &
           +4.0_ark*F12a2a2a*s1*s2a*(s2a**2-3.0_ark*s2b**2)                   &
           +12.0_ark*F12a3z3z*s1*(s2a*(s3z**2-0.5_ark*s3x**2-0.5_ark*s3y**2)   &
                           +s2b*sqrt(3.0_ark)*0.5_ark*(s3x**2-s3y**2))        &
           +24.0_ark*F12a3z4z*s1*(s2a*(s3z*s4z-0.5_ark*s3x*s4x-0.5_ark*s3y*s4y)&
                           +s2b*sqrt(3.0_ark)*0.5_ark*(s3x*s4x-s3y*s4y))      &
           +12.0_ark*F12a4z4z*s1*(s2a*(s4z**2-0.5_ark*s4x**2-0.5_ark*s4y**2)   &
                           +s2b*sqrt(3.0_ark)*0.5_ark*(s4x**2-s4y**2))        &
           +24.0_ark*F13x3y3z*s1*s3x*s3y*s3z                                 &
           +24.0_ark*F13x3y4z*s1*(s3x*s3y*s4z+s3x*s3z*s4y+s3y*s3z*s4x)       &
           +24.0_ark*F13x4y4z*s1*(s3x*s4y*s4z+s3y*s4x*s4z+s3z*s4x*s4y)       &
           +24.0_ark*F14x4y4z*s1*s4x*s4y*s4z                                 &
           +F2a2a2a2a*(s2a**2+s2b**2)**2                                    &
           +6.0_ark*F2a2a3z3z*s2a**2*s3z**2                                  &
           +6.0_ark*F2b2b3z3z*s2b**2*s3z**2                                  &
           +6.0_ark*0.25_ark*(s3x**2+s3y**2)*((F2a2a3z3z+3.0_ark*F2b2b3z3z)*s2a**2+(3.0_ark*F2a2a3z3z+F2b2b3z3z)*s2b**2)  &
           +12.0_ark*0.25_ark*sqrt(3.0_ark)*(s3y**2-s3x**2)*s2a*s2b*(F2a2a3z3z-F2b2b3z3z)                            &
           +12.0_ark*F2a2a3z4z*s2a**2*s3z*s4z                                &
           +12.0_ark*F2b2b3z4z*s2b**2*s3z*s4z                                &
           +12.0_ark*0.25_ark*(s3x*s4x+s3y*s4y)*((F2a2a3z4z+3.0_ark*F2b2b3z4z)*s2a**2+(3.0_ark*F2a2a3z4z+F2b2b3z4z)*s2b**2)  &
           +24.0_ark*0.25_ark*sqrt(3.0_ark)*(s3y*s4y-s3x*s4x)*s2a*s2b          &
                          *(F2a2a3z4z-F2b2b3z4z)                            &
           +6.0_ark*F2a2a4z4z*s2a**2*s4z**2                                  &
           +6.0_ark*F2b2b4z4z*s2b**2*s4z**2                                  &
           +6.0_ark*0.25_ark*(s4x**2+s4y**2)*((F2a2a4z4z+3.0_ark*F2b2b4z4z)*s2a**2+(3.0_ark*F2a2a4z4z+F2b2b4z4z)*s2b**2)  &
           +12.0_ark*0.25_ark*sqrt(3.0_ark)*(s4y**2-s4x**2)*s2a*s2b            &
                          *(F2a2a4z4z-F2b2b4z4z)                            &
           +24.0_ark*F2a3x3y4z*(s2a*(s3x*s3y*s4z-0.5_ark*s3x*s3z*s4y          &
                            -0.5_ark*s3y*s3z*s4x)                            &
                      +s2b*sqrt(3.0_ark)*0.5_ark*s3z*(s3y*s4x-s3x*s4y))     &
           +24.0_ark*F2a3z4x4y*(s2a*(s3z*s4x*s4y-0.5_ark*s3x*s4y*s4z          &
                   -0.5_ark*s3y*s4x*s4z)                                     &
                      +sqrt(3.0_ark)*s2b*0.5_ark*s4z*(s3x*s4y-s3y*s4x))       &
           +F3x3x3x3x*(s3x**4+s3y**4+s3z**4)                                &
           +6.0_ark*F3x3x3y3y*(s3x**2*s3y**2+s3x**2*s3z**2+s3y**2*s3z**2)    &
           +4.0_ark*F3x3x3x4x*(s3x**3*s4x+s3y**3*s4y+s3z**3*s4z)             &
           +12.0_ark*F3x3x3y4y*(s3x**2*s3y*s4y+s3x**2*s3z*s4z+s3y**2*s3x*s4x &
                      +s3y**2*s3z*s4z+s3z**2*s3x*s4x+s3z**2*s3y*s4y)        &
           +6.0_ark*F3x3x4x4x*(s3x**2*s4x**2+s3y**2*s4y**2+s3z**2*s4z**2)    &
           +24.0_ark*F3x3y4x4y*(s3x*s3y*s4x*s4y+s3x*s3z*s4x*s4z              &
                                 +s3y*s3z*s4y*s4z)                          &
           +6.0_ark*F3x3x4y4y*(s3x**2*s4y**2+s3x**2*s4z**2+s3y**2*s4x**2     &
                      +s3y**2*s4z**2+s3z**2*s4x**2+s3z**2*s4y**2)           &
           +4.0_ark*F3x4x4x4x*(s3x*s4x**3+s3y*s4y**3+s3z*s4z**3)             &
           +12.0_ark*F3x4x4y4y*(s3x*s4x*s4y**2+s3x*s4x*s4z**2+s3y*s4y*s4x**2 &
                      +s3y*s4y*s4z**2+s3z*s4z*s4x**2+s3z*s4z*s4y**2)        &
           +F4x4x4x4x*(s4x**4+s4y**4+s4z**4)                                &
           +6.0_ark*F4x4x4y4y*(s4x**2*s4y**2+s4x**2*s4z**2+s4y**2*s4z**2)    &
            )/24.0_ark
            !
          f =  (v0+v2+v3+v4)*1.0e-11/planck/vellgt

      if (verbose>=6) write(out,"('MLpoten_xy4_zzz/end')") 
      !
 end function MLpoten_xy4_ZZZ



end module pot_xy4

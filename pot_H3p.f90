!
!  This unit is for a user defined potential 
!
module pot_user
  use accuracy
  use moltype

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


 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
       !
       f = 0
       !
  end subroutine MLdipole

 !
 ! Defining potential energy function (built for SO2)

 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   f = MLpoten_H3p_singlet(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !


      function MLpoten_H3p_singlet(ncoords,natoms,local,xyz,force) result(f) 
      !
      !subroutine MLpoten_H3p_singlet(r1,r2,r3,pot,idimp)
      implicit none
      !
      integer(ik),intent(in) ::  ncoords,natoms
      real(ark),intent(in)   ::  local(ncoords)
      real(ark),intent(in)   ::  xyz(natoms,3)
      real(ark),intent(in)   ::  force(:)
      real(ark),parameter :: tocm = 219474.63067_ark
      real(ark)   :: x(297),xmass(3),x1,x2,x3
      real(ark)   ::  f


      !integer(ik) :: nparam,npol,order
      integer(ik) :: inst,nst,idimp,idimp2
      real(ark) ::  r1,r2,r3,pot(3),dx(3,3),ac,Vas,vr
      real(ark) ::  h(3,3),tb(5),dim(3,3),xdiag,xodiag
      !real(ark) ::  gam1,x0,r0,beta,x
      real(ark) ::  minad,sumad,rmax
      !    
      !     if idimp=0  the three-body terms are used
      !     else        pure DIM potential

      !  isolate points with two large distances > rmax
      idimp = 0
      !
      dx(1,:)=xyz(2,:)-xyz(1,:)
      dx(2,:)=xyz(3,:)-xyz(1,:)
      dx(3,:)=xyz(3,:)-xyz(2,:)
      !
      r1=sqrt(sum(dx(1,:)**2))/bohr
      r2=sqrt(sum(dx(2,:)**2))/bohr
      r3=sqrt(sum(dx(3,:)**2))/bohr
      !
      !r1=local(1)/bohr
      !r2=local(2)/bohr
      !r3=local(3)/bohr
      !
      rmax=9.0_ark
      minad=min(r1,r2,r3)
      sumad=r1+r2+r3
      if (sumad-minad.gt.2.0_ark*rmax.or.minad.gt.3.5_ark) then
      !  pure DIM
         idimp2=1
      else
      !  if idimp=0 potential with three-body terms 
      !  if idimp=1 pure DIM
         idimp2=idimp
      end if

      !
      !  parameters of polorder and switch are set in block data co3bdnl
      !  parameters of coef3 are set in block data co3bdl
      !  attention: the parameter "nn" should be set as in co3bdl
      !
      !  "co3bdl" may be created running the code in mblockd.f
      !
      !  determine number of parameters
      !
      call co3bdnl()
      call co3bdl(x)
      !
      np(1)=nparam(order(1,1))+nparam(order(2,1))
      np(2)=nparam(order(1,2))
      np(3)=nparam(order(1,3)) 
      np(4)=nparam(order(1,4)) 
      np(5)=nparam(order(1,5)) 

      tb(1)=0
      tb(2)=0
      tb(3)=0
      tb(4)=0
      tb(5)=0

      if (idimp2.eq.0) then
      inst=1

      do nst=1,5
         call thrbody(x,r1,r2,r3,tb(nst),nst,inst)
         inst=inst+np(nst)
      enddo
      end if

      xdiag=tb(1)
      xodiag=tb(2)

      call dimpot(r1,r2,r3,xdiag,xodiag,h,dim,pot)

        pot(1)=pot(1) + tb(3)
        pot(2)=pot(2) + tb(4)
        pot(3)=pot(3) + tb(5)
        !
        f = (pot(1) + 1.343835625028_ark)*tocm
        !
        !f = ((r1-1.910377821)**2+(r2-1.910377821)**2+(r3-1.910377821)**2)*40000.0
        !
        ac =  0
        !
        call potvAC(ac,r1,r2,r3)
        if(ac.lt.-200.0_ark) ac=-200.0_ark
        if(ac.gt.0.0_ark) ac=0
        !
        x1 = 1.0_ark
        x2 = 2.0_ark
        x3 = 3.0_ark
        xmass = 1.00782505_ark
        !
        ac = -ac*( x1/xmass(1)+x1/xmass(2)+x1/xmass(3) )/x3/tocm*1836.15_ark/1822.89_ark
        !
        ! Asymetric part of AC
        Vas = 0
        !call potvACasym(Vas,r1,r2,r3)
        !Vas = Vas*(x1/xmass(1)-x1/xmass(2))/x3/cmtoau*1836.15/1822.89
        !
        ! Relativistic correction
        vr = 0
        call potvRCb(vr,r1,r2,r3)
        if(abs(vr).gt.10.0_ark) vr=0
        vr = vr/tocm
        !
        ! BO+AC+Rel
        f = f + ac + Vas + vr
        !
      end function MLpoten_H3p_singlet


      subroutine thrbody(x,d1,d2,d3,v,nst,inst)
      implicit none
      real(ark) ::  x(:),d1,d2,d3,v
      integer(ik) :: nst,inst
      real(ark) ::  R(3),A(3,3),Q(3),q1,q2,q3
      real(ark) ::  g1p,g2p,g3p,rho
      integer(ik) :: i,j,k,l,num,nr
      integer(ik) :: idimp
      !common/polorder/npol(5),order(5,5)
      !common/switch/gam1(5,5),x0(5,5),r0(5,5),beta(5,5)
      !common/param/np(5)
      R(1)=d1
      R(2)=d2
      R(3)=d3
      
      rho=sqrt((R(1)**2+R(2)**2+R(3)**2)/sqrt(3.0_ark))
      A(1,1)=sqrt(1._ark/3._ark)
      A(1,2)=A(1,1)
      A(1,3)=A(1,1)
      A(2,1)=0._ark
      A(2,2)=sqrt(1._ark/2._ark)
      A(2,3)=-A(2,2)
      A(3,1)=sqrt(2._ark/3._ark)
      A(3,2)=-sqrt(1._ark/6._ark)
      A(3,3)=A(3,2)
      
      num=inst-1
      
      V=0
      
      do nr=1,npol(nst)
        
        do i=1,3
          Q(i)=0
          do j=1,3
            Q(i)=Q(i)+A(i,j)*(1._ark-exp(-beta(nr,nst)*(R(j)/R0(nr,nst)-1._ark)))&
                 /beta(nr,nst)
          enddo
        enddo
                   
        q1=Q(1)
        q2=Q(2)**2+Q(3)**2
        q3=Q(3)**3-3._ark*Q(3)*Q(2)**2
          
        do l=0,order(nr,nst)
          do i=0,l
            g1p=pow(q1,i)
            do j=0,(l-i),2
              g2p=pow(q2,j/2)
              k=l-i-j
              if(mod(k,3).EQ.0) then
                g3p=pow(q3,k/3)
                num=num+1
                  V=V+x(num)*g1p*g2p*g3p             &
                    *damp(gam1(nr,nst),d1,x0(nr,nst))&
                    *damp(gam1(nr,nst),d2,x0(nr,nst))&
                    *damp(gam1(nr,nst),d3,x0(nr,nst))&
                    *damp(gam1(nr,nst),q1,x0(nr,nst))&     ! 16/1/2011
                    *damp(gam1(nr,nst),q2,x0(nr,nst))&     ! 16/1/2011
                    *damp(gam1(nr,nst),q3,x0(nr,nst))     ! 16/1/2011
!     &               *damp(gam1(nr,nst),rho,x0(nr,nst))
                endif

            enddo
          enddo
        enddo
      enddo  ! end of npol loop
!
!  avoid instabilities and small distances
!
            V=V*(1.0_ark-damp(10.0_ark,d1,0.8_ark))&
               *(1.0_ark-damp(10.0_ark,d2,0.8_ark))&
               *(1.0_ark-damp(10.0_ark,d3,0.8_ark))

!
!   the three-body term may misbehave, exclude it if the conditions
!   below are satified. The potential will then be pure DIM, which is
!   good at large distances.
!
!      if(d1+d2+d3.gt.13.0d0.and.min(d1,d2,d3).gt.2.5d0) v=0.0d0   !
!      not used

      !
      end subroutine thrbody
      !
      function pow(x,i)
      implicit none
      !
      real(ark) ::  pow
      real(ark),intent(in) ::  x
      integer(ik),intent(in) :: i
      ! calculates x**i, with 0**0=1
      if (i.eq.0) then
        pow=1.0_ark
      else
        pow=x**i
      end if
      end function pow
          
      function damp(gam,q,qq)
      real(ark),intent(in) :: q,gam,qq
      real(ark) :: damp
      !
      damp=1._ark/(1._ark+exp(gam*(q-qq)))
      !
      end function damp
        
      subroutine dimpot(r1,r2,r3,xdiag,xodiag,h,dim,pot)
      !
      implicit none
      real(ark) ::  h(3,3),r1,r2,r3,pot(3),z(3,3)
      real(ark) ::  xdiag,xodiag,dim(3,3)
      integer(ik) :: nrot,i,j

!  a. alijah 17/02/2010
!  xdiag:   three-body term to be added to the diagonals
!  xodiag:  three-body term to be added to the off-diagonals
!  simplification to get rid of tb(1-7)
!    e --> tb(1)
!    x --> tb(2)
!    a --> tb(3)
!    b --> tb(4)
!    c --> tb(5)
!    d --> tb(6)
!   xx --> tb(7)

      H(1,1)=xh2pd(r2)+ah2pd(r2)+xh2pd(r3)+ah2pd(r3)
      H(1,1)=0.5_ark*H(1,1)+pothhx(r1)+xdiag-1.0_ark

      H(2,2)=xh2pd(r1)+ah2pd(r1)+xh2pd(r3)+ah2pd(r3)
      H(2,2)=0.5_ark*H(2,2)+pothhx(r2)+xdiag-1.0_ark

      H(3,3)=xh2pd(r1)+ah2pd(r1)+xh2pd(r2)+ah2pd(r2)
      H(3,3)=0.5_ark*H(3,3)+pothhx(r3)+xdiag-1.0_ark

      H(1,2)=0.5_ark*(xh2pd(r3)-ah2pd(r3)-xodiag**2)
      H(2,1)=H(1,2)

      H(1,3)=0.5_ark*(xh2pd(r2)-ah2pd(r2)-xodiag**2)
      H(3,1)=H(1,3)

      H(2,3)=0.5_ark*(xh2pd(r1)-ah2pd(r1)-xodiag**2)
      H(3,2)=H(2,3)

      do i=1,3
        do j=1,3
          dim(i,j)=h(i,j)
        enddo
      enddo

      call jacobi2(h,3,3,pot,z,nrot)
      call piksrt(3,pot)
      return
      end subroutine dimpot
        
      !  calculates number of parameters up to order "order"
      function nparam (order)
      implicit none
      integer(ik) :: nparam
      integer(ik) :: i,j,k,l,order
      nparam=0
      do l=0,order
        do i=0,l
          do j=0,(l-i),2
            k=l-i-j
            if(mod(k,3).eq.0) then
              nparam=nparam+1
            endif
            end do
          end do
        end do
        return
      end function nparam


      subroutine co3bdnl
!  sets non-linear parameters of three-body term
      implicit none
      !common/polorder/npol(5),order(5,5)
      !common/switch/gam1(5,5),x0(5,5),r0(5,5),beta(5,5)
               
!  orders of polynomials:
!    if order=-1, no polynomial
!    if order= 0, only constant term
!    if order= 1, terms of order 0 and 1
!    if order= 2, terms of order 0, 1 and 2
!    etc.
!
! number of polynomials of the diagonal
      npol(1) =2
! first polynomial of the diagonal
      gam1(1,1)=0.3_ark
      r0(1,1)  =1.65_ark
      beta(1,1)=1.3_ark
      x0(1,1)  =12.0_ark
      order(1,1)=-1
! second polynomial of the diagonal
      gam1(2,1)=0.3_ark
      r0(2,1)  =2.5_ark
      beta(2,1)=1.0_ark
      x0(2,1)  =14.0_ark
      order(2,1)=15

! polynomial of the off-diagonal
      npol(2)  =1
      gam1(1,2)=0.3_ark
      r0(1,2)  =2.5_ark
      beta(1,2)=1.0_ark
      x0(1,2)  =14.0_ark
      order(1,2)=13
      
! polynomial outside DIM matrix for ground state
      npol(3)  =1
      gam1(1,3)=0.3_ark
      r0(1,3)  =2.5_ark
      beta(1,3)=1.0_ark
      x0(1,3)  =10.0_ark
      order(1,3)=-1
      
! polynomial outside DIM matrix for first excited state
      npol(4)  =1
      gam1(1,4)=0.3_ark
      r0(1,4)  =2.5_ark
      beta(1,4)=1.0_ark
      x0(1,4)  =12.0_ark
      order(1,4)=-1
      
! polynomial outside DIM matrix for second excited state
      npol(5)  =1
      gam1(1,5)=0.3_ark
      r0(1,5)  =2.5_ark
      beta(1,5)=1.0_ark
      x0(1,5)  =12.0_ark
      order(1,5)=-1
      
      end subroutine co3bdnl

      
      subroutine co3bdl(x)
      real(ark) ::  x(297)
      x =  (/ 0.2038135752082E-01 ,& 
             -0.4502653609961E-02 , &
             -0.2049733139575E-01 , &
              0.1388978585601E-01 , &
              0.2624774724245E-01 , &
             -0.5493988282979E-03 , &
              0.2584048919380E-01 , &
              0.3755598887801E-01 , &
              0.1371347345412E-01 , &
             -0.1439545862377E-01 , &
             -0.2458681911230E-01 , &
             -0.8049705065787E-02 , &
              0.3551620244980E-01 , &
              0.5522616580129E-01 , &
             -0.1107459589839E+00 , &
             -0.2481598639861E-02 , &
             -0.7615165784955E-02 , &
             -0.1556981354952E-02 , &
             -0.2044748142362E-01 , &
              0.3450882650213E-04 , &
             -0.2876006998122E-01 , &
              0.3585752798244E-02 , &
              0.5190435051918E-01 , &
             -0.1708707539365E-02 , &
              0.2274964749813E-01 , &
              0.1964472758118E-03 , &
              0.2397969597951E-02 , &
              0.2068273315672E-03 , &
             -0.5980014428496E-01 , &
             -0.8253587409854E-02 , &
             -0.4225731641054E-01 , &
             -0.1474103238434E-01 , &
             -0.9330268949270E-02 , &
             -0.6015709415078E-01 , &
              0.1426827609539E+00 , &
              0.1501764813838E-05 , &
             -0.7821730524302E-01 , &
             -0.9919969737530E-01 , &
              0.8972355723381E-01 , &
             -0.1051881723106E-01 , &
             -0.1506378054619E+00 , &
             -0.1030388195068E-01 , &
             -0.1700490526855E-01 , &
             -0.2058380693197E+00 , &
             -0.8895791321993E-01 , &
             -0.1030693203211E+00 , &
              0.1231466904283E+00 , &
              0.2777793817222E-01 , &
              0.1696047693258E-03 , &
              0.2040153890848E+00 , &
              0.3578749597073E+00 , &
             -0.5582627654076E-01 , &
             -0.4295857623219E-01 , &
              0.1307931524934E-03 , &
             -0.3089929744601E-01 , &
             -0.3914486244321E-01 , &
             -0.2266425266862E-01 , &
             -0.3675921559334E+00 , &
              0.1150405555964E+00 , &
              0.8040013313293E+00 , &
              0.8343981951475E-01 , &
              0.9832771420479E+00 , &
              0.1498659372330E+01 , &
              0.1561497449875E+01 , &
              0.4981927871704E+00 , &
              0.3854436874390E+00 , &
              0.5157361552119E-01 , &
             -0.5597009789199E-02 , &
             -0.8208463899791E-02 , &
              0.2799553871155E+00 , &
              0.1264895796776E+00 , &
             -0.8795589674264E-03 , &
              0.4569451212883E+00 , &
              0.3874383270741E+00 , &
              0.2551791965961E+00 , &
              0.1009184122086E+01 , &
              0.2745246887207E+00 , &
              0.1018496394157E+01 , &
              0.1291464686394E+01 , &
              0.1277048707008E+01 , &
              0.5290365219116E-01 , &
              0.9185030460358E+00 , &
              0.3099000081420E-01 , &
              0.1510376925580E-02 , &
              0.9264567052014E-04 , &
              0.4045693203807E-01 , &
              0.3528299555182E-01 , &
              0.1624061316252E+00 , &
              0.7760862112045E+00 , &
              0.1307412516326E-01 , &
              0.1503127068281E+00 , &
             -0.1806402653456E+00 , &
              0.1454018354416E+00 , &
             -0.8989291787148E+00 , &
             -0.2468922138214E+01 , &
             -0.2391477674246E+00 , &
             -0.2874656677246E+01 , &
             -0.3005650281906E+01 , &
             -0.2573872089386E+01 , &
             -0.3544934093952E+00 , &
              0.1113087534904E+00 , &
             -0.1227685366757E-02 , &
              0.1945724524558E-01 , &
              0.2700309082866E-01 , &
             -0.1965133659542E-01 , &
             -0.2602054774761E+00 , &
             -0.1566321551800E+00 , &
              0.1117114797235E+00 , &
             -0.3906232118607E+00 , &
             -0.4097225964069E+00 , &
             -0.4356949925423E+00 , &
              0.4296246469021E+00 , &
             -0.2132443428040E+01 , &
             -0.1110988020897E+01 , &
             -0.1931788206100E+01 , &
             -0.3352178335190E+01 , &
             -0.2033104300499E+00 , &
             -0.3341197252274E+01 , &
             -0.2821970701218E+01 , &
             -0.2959844112396E+01 , &
             -0.3598988354206E+00 , &
             -0.6201191544533E+00 , &
              0.8133341680150E-06 , &
              0.1234306255355E-02 , &
              0.1268017664552E-01 , &
             -0.5187008995563E-02 , &
             -0.1626813709736E+00 , &
             -0.3367761671543E+00 , &
              0.1109622512013E-01 , &
             -0.8243380188942E+00 , &
             -0.4655699804425E-01 , &
             -0.3448127508163E+00 , &
             -0.6386389732361E+00 , &
             -0.7093213200569E+00 , &
             -0.3560272976756E-02 , &
              0.2475702613592E+00 , &
              0.3179576396942E+00 , &
              0.7398064136505E+00 , &
              0.1610691308975E+01 , &
              0.4166177272797E+01 , &
              0.9598708748817E+00 , &
              0.3224317550659E+01 , &
              0.2147500514984E+01 , &
              0.7551590204239E+00 , &
             -0.2166554182768E+00 , &
             -0.2505114376545E+00 , &
              0.9207403287292E-02 , &
              0.1126132156060E-04 , &
             -0.2862761961296E-02 , &
             -0.1534178009024E-03 , &
              0.9977171197534E-02 , &
             -0.1740714460611E+00 , &
              0.1246334798634E-01 , &
             -0.2196481227875E+00 , &
             -0.1066186577082E+00 , &
             -0.6484547257423E-01 , &
             -0.2513209283352E+00 , &
              0.3015027567744E-01 , &
             -0.3334693908691E+00 , &
              0.5404766201973E+00 , &
              0.8669694662094E+00 , &
              0.5779470801353E+00 , &
              0.1180241554976E+00 , &
              0.2760182857513E+01 , &
              0.1963732719421E+01 , &
              0.2369036197662E+01 , &
              0.4611462593079E+01 , &
              0.8045859336853E+00 , &
              0.3192297935486E+01 , &
              0.2113111972809E+01 , &
              0.1158948063850E+01 , &
             -0.2899852395058E-01 , &
              0.3210568428040E-01 , &
              0.4643308464438E-02 , &
              0.5794972926378E-01 , &
              0.2163177281618E+00 , &
              0.3456968755700E-07 , &
             -0.7665733825490E-08 , &
              0.3650718182325E-01 , &
             -0.6945794820786E-01 , &
             -0.5244998261333E-01 , &
             -0.1063322369009E-01 , &
              0.1378841549158E+00 , &
             -0.3282044827938E+00 , &
             -0.5219735205173E-01 , &
             -0.1209799572825E+00 , &
             -0.4702143371105E-02 , &
             -0.6626180559397E-01 , &
              0.2407051852060E-05 , &
              0.1468833833933E+00 , &
             -0.9982571646105E-04 , &
              0.2347189188004E-01 , &
              0.4199467948638E-03 , &
              0.1984704583883E+00 , &
             -0.1109232529998E+00 , &
             -0.1109396368265E+00 , &
             -0.2261913381517E-01 , &
              0.2691740883165E-04 , &
              0.2081657201052E+00 , &
              0.1839612275362E+00 , &
              0.5080552101135E+00 , &
             -0.1613639108837E-01 , &
              0.3844989538193E+00 , &
              0.1302360445261E+00 , &
             -0.3850560188293E+00 , &
             -0.1574946641922E+00 , &
             -0.1519335508347E+00 , &
             -0.2564406394958E+00 , &
              0.4450597465038E+00 , &
             -0.2165480405092E+00 , &
             -0.3114369213581E+00 , &
              0.2186658978462E+00 , &
              0.6476427316666E+00 , &
              0.3730938732624E+00 , &
             -0.1092024073005E+00 , &
             -0.9590282291174E-01 , &
             -0.1192676718347E-03 , &
             -0.6356548666954E+00 , &
              0.1947915554047E+00 , &
              0.6025955080986E+00 , &
              0.2409428507090E+00 , &
              0.8307718038559E+00 , &
              0.7525328397751E+00 , &
              0.1109323501587E+01 , &
              0.6152945756912E+00 , &
              0.1256046380149E-04 , &
              0.1692902892828E+00 , &
              0.6652852892876E-01 , &
             -0.2153564430773E-01 , &
             -0.3678647056222E-01 , &
              0.8088537305593E-01 , &
              0.5307267885655E-02 , &
              0.6883329153061E+00 , &
              0.2174600362778E+01 , &
              0.6245110630989E+00 , &
              0.2579672574997E+01 , &
              0.4459482669830E+01 , &
              0.2777974843979E+01 , &
              0.3691339790821E+00 , &
              0.5397902727127E+00 , &
              0.1189048052765E-02 , &
             -0.2816260792315E-01 , &
             -0.3072485029697E+00 , &
              0.7174132466316E+00 , &
             -0.5216906666756E+00 , &
             -0.4354291260242E+00 , &
             -0.3040901422501E+00 , &
             -0.3592750523239E-02 , &
             -0.9823914766312E+00 , &
             -0.6712844371796E+00 , &
             -0.1427285722457E-02 , &
              0.6134739518166E+00 , &
              0.6469817757607E+00 , &
              0.2023749053478E+00 , &
              0.2512786984444E+00 , &
              0.8785905838013E+00 , &
             -0.1337261050940E+00 , &
             -0.6628892384470E-02 , &
             -0.3364108800888E+00 , &
             -0.5055920407176E-01 , &
              0.1961876749992E+00 , &
              0.1384985866025E-02 , &
             -0.7081563025713E-01 , &
             -0.7799867987633E+00 , &
              0.7183402776718E-01 , &
             -0.1918944001198E+01 , &
             -0.1496339321136E+01 , &
             -0.2778274774551E+01 , &
             -0.5837337970734E+01 , &
             -0.1115664243698E+01 , &
             -0.3906371831894E+01 , &
             -0.4546409130096E+01 , &
             -0.3241823196411E+01 , &
             -0.6220843270421E-01 , &
              0.2666761279106E+00 , &
             -0.6738465279341E-01 , &
             -0.6620140373707E-01 , &
             -0.5357830226421E-01 , &
              0.2084175050259E+00 , &
             -0.3408109247684E+00 , &
              0.2126413285732E+00 , &
              0.4915786087513E+00 , &
             -0.5661582350731E+00 , &
              0.2188280783594E-01 , &
              0.5355895336834E-04 , &
              0.3853269517422E+00 , &
             -0.1088703632355E+01 , &
             -0.7743591070175E+00 , &
             -0.2267065197229E+00 , &
             -0.2393134385347E+00 , &
             -0.3650400936604E+00 , &
              0.1053227926604E-02 , &
             -0.1342126488686E+01 , &
             -0.1051973462105E+01 , &
             -0.1317051351070E+00 , &
              0.1048775576055E-01 , &
             -0.6959872320294E-02 /)
      end subroutine co3bdl 



!===============================================================
! EHFACE2U POTENTIAL CURVE FOR H2+ ( X 2^sigma^(+)_(g) )
! ## 2006-06-22 ##
! USING AB INITIO POINTS FROM:
! D. M. BISHOP AND R. W. WETMORE, MOL. PHYS. 26(1),145 (1972)
! RANGE of R: 0.6 to 10 a0
!
!  Initial fit  
!  RMS(m= 95 )= 67.6787429191265630  cm-1
!  Final fit  
!  RMS(m= 95 )= 0.675527097585787786E-02  cm-1
!===============================================================     
      function XH2PD(R)
      IMPLICIT NONE
      real(ark),intent(in) :: r
      real(ark) :: dd(20),XH2PD,re
      integer(ik) :: nlow,nupp,ii,iexp,i
      real(ark) :: D1(20),D2(20),C(20),GAMMA,agp,AI(11),g0,g1,g2,r0,rm,x,pol1,VHF,gam,rhh
      real(ark) :: DISP,DAMPI
      !
      !COMMON/POTEN1/D1(20),D2(20),C(20),GAMMA,AGP,AI(11),R0,RM
      !COMMON/LIM1/NLOW,NUPP
      !COMMON/TH1/G0,G1,G2
      !COMMON/EXPOE1/RE,IEXP
      !COMMON/ASYEXC1/CATILD,ATILD(2),ALPHT,GTILD
      !
      NLOW=4
      NUPP=11
      DO II=NLOW,NUPP
        D1(II)=AN(II)
        D2(II)=BN(II)
      ENDDO
      IEXP=1
      RE=2.0_ark
      C(4)=2.250_ark
      C(5)=0.0_ark
      C(6)=7.5000_ark
      C(7)=53.25_ark
      C(8)=65.625_ark
      C(9)=886.5_ark
      C(10)=1063.125_ark
      C(11)=21217.5_ark
      GAMMA=2.5_ark
      AGP=-0.180395506614312_ark
      AI(1)=-0.897676265677028185_ark
      AI(2)=-0.771599228853606545_ark
      AI(3)=-0.245766669963638718_ark
      AI(4)=-0.788889284685244524E-01
      AI(5)=-0.252032558464952844E-01
      AI(6)=0.681894227468654839E-02
      AI(7)=0.655940943163255976E-03
      AI(8)=-0.531288172311992135E-03
      AI(9)=0.890418306898401330E-04
      AI(10)=-0.666314834544138477E-05
      AI(11)=0.194182431833699709E-06
      G0=0.960151039243191562_ark
      G1=-0.353946173857859037_ark
      G2=-0.496213155382123572_ark
      R0=3.4641_ark
      RM=0.11000000D+02
      X=R-RE

      pol1=ai(11)
      do i=10,1,-1
         pol1=pol1*x + ai(i)
      end do
      pol1=pol1*x + 1.0_ark

      GAM=G0*(1.0_ark+G1*TANH(G2*X))
      VHF=-AGP/(R**IEXP)*POL1*EXP(-GAM*X)
!
      RHH=0.5_ark*(RM+GAMMA*R0)
      X=R/RHH    
!
      DISP=0
      DO I=NLOW,NUPP
        DAMPI=(1.0_ark-EXP(-D1(I)*X-D2(I)*X**2))**I
        DD(I)=DAMPI
        DISP=DISP-C(I)*DAMPI*R**(-I)
      enddo
      XH2PD=VHF+DISP
      
      end  function XH2PD


!===============================================================
! POTENTIAL CURVE FOR H2+ ( A 2^sigma^(+)_(u) )
! ## 2006-06-22 ##
! USING AB INITIO POINTS FROM:
! J. M. PEEK, JCP 43(9), 3004 (1965)
! RANGE of R: 3.5 to 15 a0 
!
!  Final fit  
!  RMS(m= 24 )= 0.125196743261995869  cm-1
!===============================================================
      real*8 function ah2pd(r)
      implicit none
      real(ark) :: r,coef(0:7),v
      integer i
      !
      coef( 0 )=  1.11773285795729826_ark
      coef( 1 )= -1.27592697554394174_ark
      coef( 2 )=  0.235612064424508216_ark
      coef( 3 )= -0.500203729467869895E-01
      coef( 4 )=  0.568627052480373801E-02
      coef( 5 )= -0.382978465312642114E-03
      coef( 6 )=  0.149149267032670154E-04
      coef( 7 )= -0.267518847221239873E-06
!
      v=coef(7)
      do i=6,0,-1
        v=v*r+coef(i)
      enddo
      ah2pd=exp(v)+xh2pd(r)
      return
      end function ah2pd


      FUNCTION POTHHX(R)
!===============================================================
!  NEWEST FIT FOR H2 GROUND STATE POTENTIAL CURVE ##2006-06-22##
!  USING AB INITIO POINTS FROM:
!  L. WOLNIEWICZ JCP 99(3),1851 (1993) ---> FIRST POINTS
!  L. WOLNIEWICZ JCP 103(5),1792 (1995) ---> CORRECTIONS    
!  RANGE of R: 0.6 to 12 a0
!
!  Initial fit  
!  RMS(m= 52 )= 35.2987773353410503  cm-1
!  Final fit
!  RMS(m= 52 )= 0.967169906620670844E-01  cm-1
!===============================================================
      IMPLICIT NONE
      real(ark) :: POTHHX,r
      real(ark) ::  DD(20)
      integer(ik) :: ii,NLOW,nupp,iexp,i
      real(ark) :: D1(20),D2(20),CATILD,GTILD,re,C(20),GAMMA,AGP,AI(9),R0,RM,DAMPI,DEXC
      real(ark) :: gam,pol1,VHF,ATILD(2),ASEXC,ALPHT,g0,g1,g2,x,rhh,DISP
      !COMMON/POTEN/D1(20),D2(20),C(20),GAMMA,AGP,AI(9),R0,RM
      !COMMON/LIM/NLOW,NUPP
      !COMMON/TH/G0,G1,G2
      !COMMON/EXPOE/RE,IEXP
      !COMMON/ASYEXC/CATILD,ATILD(2),ALPHT,GTILD
      !
      NLOW=6
      NUPP=16
      DO II=NLOW,NUPP
        D1(II)=AN(II)
        D2(II)=BN(II)
      ENDDO
      CATILD=-0.8205_ark
      ATILD(1)=0
      ATILD(2)=0
      ALPHT=2.5_ark
      GTILD=2.0_ark
      IEXP=1
      RE=0.14010000E+01
      C(6)=0.64990000E+01
      C(7)=0.0D0
      C(8)=0.12440000E+03
      C(9)=0
      C(10)=0.32858000E+04
      C(11)=-0.34750000E+04
      C(12)=0.12150000E+06
      C(13)=-0.29140000E+06
      C(14)=0.60610000E+07
      C(15)=-0.23050000E+08
      C(16)=0.39380000E+09
      GAMMA=2.5_ark
      AGP=0.229794389784158_ark
      AI(1)=1.74651398886700093_ark
      AI(2)=0.631036031819560028_ark
      AI(3)=0.747363488024733624_ark
      AI(4)=0.956724297662875783E-01
      AI(5)=0.131320504483065703_ark
      AI(6)=-0.812200084994067194E-07
      AI(7)=0.119803887928360935E-01
      AI(8)=-0.212584227748381302E-02
      AI(9)=0.509125901134908042E-03
      G0=1.02072511539524680_ark
      G1=1.82599688484061118_ark
      G2=0.269916332495592104_ark
      R0=0.69282032E+01
      RM=0.11000000E+02
      X=R-RE

       pol1=ai(9)
       do i=8,1,-1              
         pol1=pol1*x + ai(i)
       end do
       pol1=pol1*x + 1.0d0 

      GAM=G0*(1.0D0+G1*TANH(G2*X))
      VHF=-AGP/(R**IEXP)*POL1*EXP(-GAM*X)
    
      ASEXC=10.0_ark
      DO I=1,2
        ASEXC=ASEXC+ATILD(I)*R**I
      ENDDO
      ASEXC=ASEXC*CATILD*R**ALPHT*EXP(-GTILD*R)
    
      RHH=0.5_ark*(RM+GAMMA*R0)
      X=R/RHH
      DEXC=(1.0_ark-EXP(-D1(NLOW)*X-D2(NLOW)*X**2))**NLOW
    
      ASEXC=ASEXC*DEXC
      VHF=VHF+ASEXC
    
      DISP=0
      DO 1 I=NLOW,NUPP
        DAMPI=(1.0_ark-EXP(-D1(I)*X-D2(I)*X**2))**I
        DD(I)=DAMPI
        DISP=DISP-C(I)*DAMPI*R**(-I)
 1    CONTINUE
      POTHHX=VHF+DISP
      RETURN
      END FUNCTION POTHHX


      FUNCTION AN(N)
      IMPLICIT NONE
      real(ark) :: an
      integer(ik) :: N
      real(ark) :: ALPH0,ALPH1
      
      ALPH0=16.36606_ark
      ALPH1=0.70172_ark
      AN=ALPH0/(real(N,ark))**ALPH1
      RETURN
      END FUNCTION AN

      FUNCTION BN(N)
      IMPLICIT NONE
      real(ark) :: bn
      integer(ik) :: N
      real(ark) :: bet0,bet1
      

      BET0=17.19338_ark
      BET1=0.09574_ark
      BN=BET0*EXP(-BET1*real(N,ark))
      RETURN  
      END FUNCTION BN



      SUBROUTINE jacobi2(a,n,np,d,v,nrot)
      implicit none
      integer(ik) :: n,np,nrot,NMAX
      real(ark) ::  a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      integer(ik) :: i,ip,iq,j
      real(ark) ::  c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do ip=1,n
        do iq=1,n
          v(ip,iq)=0.d0
        enddo
        v(ip,ip)=1.d0
      enddo
      do ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
      enddo
      nrot=0
      do i=1,50
        sm=0.d0
        do ip=1,n-1
          do iq=ip+1,n
            sm=sm+abs(a(ip,iq))
          enddo
        enddo
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do ip=1,n-1
          do iq=ip+1,n
            g=100.d0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.d0
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/sqrt(1+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
              enddo
              do j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
              enddo
              nrot=nrot+1
            endif
          enddo
        enddo
        do ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
        enddo
      enddo
      stop 'too many iterations in jacobi'    ! AA
      !      pause 'too many iterations in jacobi'   ! original
      !      return                                  ! original
      END SUBROUTINE jacobi2
      !  (C) Copr. 1986-92 Numerical Recipes Software <%8=$j,)]#1_%.


      SUBROUTINE piksrt(n,arr)
      implicit none
      integer(ik) n
      real(ark) :: arr(n)
      integer(ik) ::  i,j
      real(ark) ::  a
      loop2: do j=2,n
        a=arr(j)
        do i=j-1,1,-1
          if(arr(i).le.a) then
            arr(i+1)=a
            cycle loop2
          endif
          arr(i+1)=arr(i)
        enddo
        i=0
        arr(i+1)=a
      enddo loop2
      return
      END SUBROUTINE piksrt
      !  (C) Copr. 1986-92 Numerical Recipes Software



      subroutine potvRCb(Vrel,P1,P2,P3)

!  Relativistic correction, fit of Barchorz 2009 data
!     UNITS: HARTREE & BOHR
!     symmetric part

      IMPLICIT real(ark) (A-H,O-Z)
      integer(ik),parameter:: ncf=31
!      dimension ft(ncf)
!      common /potential/ cv(ncf), der(ncf)
      dimension cv(ncf), ft(ncf)
      integer(ik) :: nv,nfmax,npot,norder,n,k,m,i
      DATA C0/0.0D0/,RE/1.65000/,BET/1.300D0/
      DATA cv/-2.58247_ark, 0.88025_ark, -0.58043_ark, -0.84413_ark, 0.31196_ark,&
     1.17078_ark, 0.31370_ark, -0.14455_ark, -0.85299_ark, -0.54332_ark, -0.22522_ark,&
     0.04732_ark, 0.34433_ark, 0.37552_ark, 0.27914_ark, 0.07864_ark, -0.01041_ark,&
     -0.06251_ark, -0.11575_ark, -0.12251_ark, -0.06413_ark, -0.00581_ark, -0.00078_ark,&
     0.00126_ark, 0.00271_ark, 0.01333_ark, 0.01845_ark, 0.01204_ark, 0.00203_ark,&
     0.00025, 0.00160_ark/

      data nv/ ncf/,nfmax/7/
      DATA ZERO/0.0_ark/,ONE/1.0_ark/,TWO/2.0_ark/,THREE/3.0_ark/

      SQ3=SQRT(THREE)
      SQ2=SQRT(TWO)
!      FACTOR=BET/RE
      DR1= (P1-RE)
      DR2= (P2-RE)
      DR3= (P3-RE)

! Displacement coordinates
       Y1=DR1
       Y2=DR2
       Y3=DR3
      SA=(Y1+Y2+Y3)/SQ3
      SX1=(DR2+DR2-DR1-DR3)/(SQ2*SQ3)
      SY1=(DR1-DR3)/SQ2
      QUAD1=SX1**2+SY1**2
      SE1=sqrt(QUAD1)
      if(abs(se1).lt.1.0e-10) then
       phi=acos(0.0_ark)
      else
       phi=acos(sx1/se1)
      endif

      npot=1
      ft(1)=1.0_ark
      do 100 norder=1,nfmax
      do 100 n=norder,0,-1
      do 100 k=0,norder-n,3
      if(mod(k,3).ne.0) goto 100
      m=norder-k-n
      if (mod(m,2) .ne. 0) goto 100
      npot=npot+1
      if(npot.gt.ncf) goto 100
      ft(npot)=sa**n * se1**(m+k) * cos(dble(k)*phi)
  100 continue
      V=ZERO
      DO 40 I=1,NV
   40 V=V+CV(I)*FT(I)
      Vrel = V
      RETURN
      END subroutine potvRCb

!  Potential from Ludwik's calculations of BO and AC
!  Calculates AC
!  input bond lengths in Bohr
      subroutine potvAC(V,p1,p2,p3)     

      implicit real(ark)(A-H,O-Z)
      real(ark) Ves, p1,p2,p3
      integer(ik),parameter :: ncf=98
      integer(ik) :: nfmax,npot,norder,k,m,n,i
      dimension ft(ncf), cv(ncf)
      data RE/1.65_ark/,BET/1.300_ark/
      data nfmax/12/
      !save icall, cv

      !icall = icall + 1
      !
      !if(icall.eq.1) then
      !  open(unit=31,status='old',file='f.31.symAC')
      !do i=1,ncf
      !  read(31,2) cv(i)
      !2       format(f20.5)
      !end do
      !close(31)
      !end if
      !
      cv( 1) =      -115.13592_ark     
      cv( 2) =        34.95581_ark     
      cv( 3) =        -8.94854_ark     
      cv( 4) =       -19.22486_ark     
      cv( 5) =        -0.47363_ark     
      cv( 6) =        20.42779_ark     
      cv( 7) =        -1.08335_ark     
      cv( 8) =        -1.18039_ark     
      cv( 9) =        -0.97418_ark     
      cv(10) =         1.61562_ark     
      cv(11) =         0.84843_ark     
      cv(12) =        -0.76828_ark     
      cv(13) =        -1.84770_ark     
      cv(14) =        -1.50647_ark     
      cv(15) =        -0.27232_ark     
      cv(16) =         0.26785_ark     
      cv(17) =        -1.12614_ark     
      cv(18) =        -0.28722_ark     
      cv(19) =         2.75205_ark     
      cv(20) =        -2.17235_ark     
      cv(21) =         1.90988_ark     
      cv(22) =        -0.38403_ark     
      cv(23) =         0.17021_ark     
      cv(24) =         0.00890_ark     
      cv(25) =        -1.87033_ark     
      cv(26) =         1.56787_ark     
      cv(27) =        -1.46615_ark     
      cv(28) =        -2.39546_ark     
      cv(29) =        -0.45831_ark     
      cv(30) =        -0.38466_ark     
      cv(31) =         0.00960_ark     
      cv(32) =        -0.00415_ark     
      cv(33) =         5.12902_ark     
      cv(34) =        -7.51230_ark     
      cv(35) =         3.68779_ark     
      cv(36) =        -1.19697_ark     
      cv(37) =         2.50128_ark     
      cv(38) =         0.59069_ark     
      cv(39) =        -1.04071_ark     
      cv(40) =         0.17463_ark     
      cv(41) =         0.02009_ark     
      cv(42) =        -1.97854_ark     
      cv(43) =        -1.28897_ark     
      cv(44) =        -1.74433_ark     
      cv(45) =         0.32559_ark     
      cv(46) =         0.94534_ark     
      cv(47) =         0.32374_ark     
      cv(48) =         0.04094_ark     
      cv(49) =         1.81070_ark     
      cv(50) =         0.03173_ark     
      cv(51) =         0.07976_ark     
      cv(52) =        -0.03801_ark     
      cv(53) =         0.00163_ark     
      cv(54) =         1.62614_ark     
      cv(55) =        -7.38655_ark     
      cv(56) =        12.42367_ark     
      cv(57) =        -2.21440_ark     
      cv(58) =         7.99644_ark     
      cv(59) =        -4.24962_ark     
      cv(60) =        -1.12576_ark     
      cv(61) =        -1.08647_ark     
      cv(62) =        -1.43964_ark     
      cv(63) =        -0.61741_ark     
      cv(64) =         0.31037_ark     
      cv(65) =        -0.00046_ark     
      cv(66) =        -0.03831_ark     
      cv(67) =        -0.01450_ark     
      cv(68) =         1.61879_ark     
      cv(69) =         3.07243_ark     
      cv(70) =        -2.27573_ark     
      cv(71) =         4.67040_ark     
      cv(72) =       -13.97697_ark     
      cv(73) =         1.39648_ark     
      cv(74) =         0.49974_ark     
      cv(75) =         0.80347_ark     
      cv(76) =         1.97912_ark     
      cv(77) =         0.91585_ark     
      cv(78) =        -0.53555_ark     
      cv(79) =        -0.01358_ark     
      cv(80) =         0.12127_ark     
      cv(81) =         0.05223_ark     
      cv(82) =        -0.00432_ark     
      cv(83) =         0.00048_ark     
      cv(84) =        -2.34873_ark     
      cv(85) =         2.20312_ark     
      cv(86) =        -3.73879_ark     
      cv(87) =        -3.41559_ark     
      cv(88) =         5.92814_ark     
      cv(89) =         0.63597_ark     
      cv(90) =         0.23146_ark     
      cv(91) =        -0.13923_ark     
      cv(92) =        -0.66705_ark     
      cv(93) =        -0.38679_ark     
      cv(94) =         0.19273_ark     
      cv(95) =         0.00978_ark     
      cv(96) =        -0.09287_ark     
      cv(97) =        -0.04204_ark     
      cv(98) =         0.00989_ark     
      !
      SQ3=SQRT(3.0_ark)
      SQ2=SQRT(2.0_ark)
      FACTOR=BET/RE
      DR1= (P1-RE)
      DR2= (P2-RE)
      DR3= (P3-RE)
       Y1=(1.0_ark-EXP(-FACTOR*DR1))/BET
       Y2=(1.0_ark-EXP(-FACTOR*DR2))/BET
       Y3=(1.0_ark-EXP(-FACTOR*DR3))/BET
      SA=(Y1+Y2+Y3)/SQ3
      SX=(Y3+Y3-Y1-Y2)/(SQ2*SQ3)
      SY=(Y2-Y1)/SQ2
      SX1=(DR3+DR3-DR1-DR2)/(SQ2*SQ3)
      SY1=(DR2-DR1)/SQ2
      QUAD1=SX1**2+SY1**2
      SE1=SQRT(QUAD1)
      if(abs(se1).lt.1.0e-10) then
        phi=acos(0.0_ark)
      else
        phi=acos(sx1/se1)
      end if
     !
      npot=1
      ft(1)=1.0_ark
      do 100 norder=1,nfmax
      do 100 n=norder,0,-1
      do 100 k=0,norder-n,3
      m=norder-k-n
      if (mod(m,2) .ne. 0) goto 100
      npot=npot+1
      if(npot.gt.ncf) go to 100
      ft(npot)=sa**n * se1**(m+k) * cos(dble(k)*phi)
  100 continue
      Vad=0.0_ark

      do i=1,ncf
       Vad=Vad + cv(i)*ft(i)
      end do
      V = Vad 
   
      return
      end subroutine potvAC


end module pot_user

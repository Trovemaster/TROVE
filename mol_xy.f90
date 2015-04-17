!
!  This unit defines all specific routines for a diatomic molecule of XY type
!
module mol_xy
  use accuracy
  use moltype

  implicit none

  public ML_b0_XY,MLpoten_xy_gener

  private
 
  integer(ik), parameter :: verbose     = 2                          ! Verbosity level

  !type(MoleculeT),save :: molec

  contains


  !
  ! Here we define structural parameters a0 for rigid XY twoatomic molecule.
  !
  subroutine ML_b0_XY(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)

     integer(ik),intent(in)  :: Npoints,Natoms

     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in) ,optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional :: rho_ref
     real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders

     real(ark)                :: r0,m1,m2,a0(molec%Natoms,3),rho
     integer(ik)             :: Nbonds,i


      if (verbose>=4) write(out,"('ML_b0_XY/start')") 

      Nbonds = size(molec%req)

      if (Nbonds/=1) then
        write(out,"('ML_b0_XY: Nbonds must be 1 in this routine, not  ',i8)") Nbonds
        stop 'ML_b0_XY: wrong Nbonds '
      endif 

      if (molec%Natoms/=2) then
        write(out,"('ML_b0_XY: Natoms must be 2 in this routine, not  ',i8)") molec%Natoms
        stop 'ML_b0_XY: wrong Natoms '
      endif 

      r0 = molec%req(1)

      m1 = molec%AtomMasses(1) ; m2 = molec%AtomMasses(2)

      a0(1,1) = 0.0_ark
      a0(1,2) = 0.0_ark
      a0(1,3) = -r0*m2/(m1+m2)
      !
      a0(2,1) = 0.0_ark
      a0(2,2) = 0.0_ark
      a0(2,3) =  r0*m1/(m1+m2)
      !
      b0(:,:,0) = a0(:,:)
      !
      !
      if (Npoints/=0) then
         ! 
         if (.not.present(rho_borders).or..not.present(rho_ref)) then  
            write(out,"('ML_b0_XY: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_XY: rho_borders or rho_ref not specified '
         endif
         !
         rho_ref = 0.0_ark
         r0 = molec%req(1)

         do i = 0,npoints
            !
            rho = r0+rho_i(i)
            !
            b0(1,1,i) = 0.0_ark
            b0(1,2,i) = 0.0_ark
            b0(1,3,i) = -rho*m2/(m1+m2)
            !
            b0(2,1,i) = 0.0_ark
            b0(2,2,i) = 0.0_ark
            b0(2,3,i) =  rho*m1/(m1+m2)
            !
            if (verbose>=4) then 
              write(out,"('b0',i4,12f12.8)") i,b0(:,:,i)
            endif
            !
         enddo
         !
      endif
      !
      if (verbose>=4) write(out,"('ML_b0_XY/end')") 

  end subroutine ML_b0_XY

  !
  ! Defining potential energy function: generalaized 
  !
  ! This type is for XY- twoatomic molecule, Morse, Dunham, and SPF types of expansion 
  !
  function MLpoten_xy_gener(ncoords,natoms,local,xyz,force) result(v) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  v

   real(ark)    :: req,a,r,y,vh2
   !
   integer :: N,k

      if (verbose>=6) write(out,"('MLpoten_xy_gener/start')") 
      !
      N = size(force(:))
      !
      req   =  molec%req(1)
      !
      r    = local(1)
      !
      !v = f(2)*y**2+f(3)*y**3+f(4)*y**4+f(5)*y**5+f(6)*y**6+f(7)*y**7+f(8)*y**8+f(9)*y**9+f(10)*y**10

    select case(trim(molec%potentype))
    case default
         write (out,"('MLpotenfunc: potential type ',a,' unknown')") trim(molec%potentype)
         stop 'MLpoten_xy_gener - bad potential'

    case('POTEN_XY_DUNHAM') 
         !
         y = (r-req)/req
         !v = f(2)*y**2*(1.0_rk+f(3)*y**1+f(4)*y**2+f(5)*y**3+f(6)*y**4+f(7)*y**5+f(8)*y**6+f(9)*y**7+f(10)*y**8)
         !
         v = 1.0_ark
         !
         do k=2,N
            v = v+force(k)*y**(k-1)
         enddo
         !
         v = force(1)*y**2*v
         !
    case('POTEN_XY_SPF') 
         !
         y = (r-req)/r
         !v = f(2)*y**2*(1.0_rk+f(3)*y**1+f(4)*y**2+f(5)*y**3+f(6)*y**4+f(7)*y**5+f(8)*y**6+f(9)*y**7+f(10)*y**8)
         !
         v = 1.0_ark
         !
         do k=2,N
            v = v+force(k)*y**(k-1)
         enddo
         !
         v = force(1)*y**2*v
         !
    case('POTEN_XY_SPF_H2') 
         !
         y = (r-req)/r
         !v = f(2)*y**2*(1.0_rk+f(3)*y**1+f(4)*y**2+f(5)*y**3+f(6)*y**4+f(7)*y**5+f(8)*y**6+f(9)*y**7+f(10)*y**8)
         !
         v = 1.0_ark
         !
         do k=2,N-4
            v = v+force(k)*y**(k-1)
         enddo
         !
         v = force(1)*y**2*v
         !
         vh2 = 0
         !
         do k=1,4
            vh2 = vh2+force(N-4+k)*y**k
         enddo
         !
         v = v + vh2
         !
    case('POTEN_XY_MORSE') 
         !
         a  =  molec%specparam( 1)
         !
         y = 1.0_rk-exp(-a*(r-req)) 
         !
         v = 0
         !
         do k=1,N
            v = v+force(k)*y**(k+1)
         enddo
         !
         !v = force(2)*y**2+force(3)*y**3+force(4)*y**4+force(5)*y**5+force(6)*y**6+force(7)*y**7+f(8)*y**8+f(9)*y**9+f(10)*y**10
         !
    end select 
    !
    if (verbose>=6) write(out,"('MLpoten_xy_gener/end')") 
 
 end function MLpoten_xy_gener


end module mol_xy

!
!  This unit defines all specific routines for a fouratomic molecule of ABCD type
!
module mol_x2y2
  use accuracy
  use moltype
  use lapack
  use pot_abcd
  use symmetry,only : sym

  implicit none

  public ML_b0_X2Y2,ML_coordinate_transform_X2Y2,ML_rotsymmetry_X2Y2,ML_symmetry_transformation_X2Y2
  !
  private
  !
  integer(ik), parameter :: verbose = 3  ! Verbosity level
  !
  contains
  !
  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter 
  ! as conjugate momenta coordinates
  !
  function ML_coordinate_transform_X2Y2(src,ndst,direct) result (dst)
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct

    !
    real(ark),dimension(ndst) :: dst
    real(ark)                 :: dsrc(size(src))
    real(ark)                 :: rho,r_eq(6)
    !
    integer(ik) :: nsrc
    !
    if (verbose>=6) write(out,"('ML_coordinate_transform_X2Y2/start')") 
    !
    if (direct) then 
       !
       dsrc(:) = src(:) - molec%local_eq(:)
       !
    else 
       !
       dsrc(:) = src(:)
       !
    endif
    !
    nsrc = size(src)
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('MLcoordinate_transform_func: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'MLcoordinate_transform_func - bad coord. type'
       !
    case('NORMAL')
       !
       if (direct) then 
           dst = src
       else
           dst = src
       endif
       !
    case('R-ALPHA-TAU')
       !
       if (direct) then
          ! 
          dst(1:5) = dsrc(1:5)
          dst(6) = src(6)
          !
      else ! not direct
          !
          dst(1:5) = dsrc(1:5)+molec%local_eq(1:5)
          dst(6) = src(6)
          !
      endif
      !
    end select
    !
    !
    if (verbose>=6) write(out,"('ML_coordinate_transform_X2Y2/end')") 
    !
    !
  end function ML_coordinate_transform_X2Y2


  ! Here we define structural parameters a0 for ABCD molecule,
  !
  subroutine ML_b0_X2Y2(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
     !
     !
     integer(ik),intent(in) :: Npoints,Natoms
     real(ark),   intent(out) :: b0(Natoms,3,0:Npoints)
     real(ark),   intent(in) ,optional  :: rho_i(0:Npoints)
     real(ark),   intent(out),optional :: rho_ref
     real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
     !
     real(ark)    ::  rho,transform(3,3),c(3,3),a0(molec%Natoms,3),CM_shift,re1,re2,re3,ae1,ae2,r_eq(6),phi,ayz,ayy,azz
     integer(ik)  ::  n,i,ix,jx,i0,in,i1,ipar,istep,j0,Nangles,Nbonds
     !
     real(ark)    :: Inert0(3),Inert(3),Inert1(3),a(3,3),b(3,1),x(5),a0_tt,a0_t(molec%Natoms,3)
     real(ark)                        :: rho_ark
      !
      if (verbose>=4) write(out,"('ML_b0_X2Y2/start')") 
      !
      if (size(molec%req)/=3) then
        write(out,"('ML_b0_X2Y2 Nbonds must be 3 in this routine, not  ',i9)") size(molec%req)
        stop 'ML_b0_X2Y2: wrong Nbonds '
      endif 

      if (molec%Natoms/=4) then
        write(out,"('ML_b0_X2Y2: Natoms must be 4 in this routine, not  ',i9)") molec%Natoms
        stop 'ML_b0_X2Y2: wrong Natoms '
      endif 
      !
      Nbonds = molec%Nbonds
      Nangles = molec%Nangles
      !
      re1   = molec%req(1)
      re2   = abs(molec%req(2))
      re3   = abs(molec%req(3))
      !
      rho = 0 
      !
      if (Nangles==2.and.molec%Ndihedrals==1) then
        !
        ae1   = molec%alphaeq(1)
        ae2   = molec%alphaeq(2)
        !
        rho = molec%taueq(1)
        !
      else
        write(out,"('Error ML_b0_X2Y2: wrong number of angles')")
        stop 'Error ML_b0_X2Y2: wrong number of angles'
      endif 
      !
      if (Npoints/=0.and.(.not.present(rho_borders).or..not.present(rho_ref))) then  
            write(out,"('ML_b0_X2Y2: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
            stop 'ML_b0_X2Y2: rho_borders or rho_ref not specified '
      endif
      !
      !
      select case(trim(molec%coords_transform))
         !
      case default
         !
         write (out,"('MLcoordinate_transform_func: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLcoordinate_transform_func - bad coord. type'
         !
      case('R-ALPHA-TAU','NORMAL')
        !
        a0(1,1) = 0.0_ark
        a0(1,2) = 0.0_ark
        a0(1,3) = 0.0_ark
        !
        a0(2,1) = 0.0_ark
        a0(2,2) = 0.0_ark
        a0(2,3) = re1
        !
        a0(3,1) = re2*sin(ae2)*cos(rho*0.5_ark)
        a0(3,2) = re2*sin(ae1)*sin(rho*0.5_ark)
        a0(3,3) = re2*cos(ae1)
        !
        a0(4,1) = re3*sin(ae2)*cos(rho*0.5_ark)
        a0(4,2) =-re3*sin(ae2)*sin(rho*0.5_ark)
        a0(4,3) = re1-re3*cos(ae2)        
        !
        do n = 1,3 
          CM_shift = sum(a0(:,n)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
          a0(:,n) = a0(:,n) - CM_shift
        enddo 
        !
      end select
      !
      b0(:,:,0) = a0(:,:)
      !
      ! We can also simulate the reference structure at each point of rho, if npoints presented
      !
      if (Npoints/=0) then
         !
         rho_ref = molec%taueq(1)
         !
         do i = 0,npoints
            !
            rho = rho_i(i)
            !
            re1   = molec%req(1)
            re2   = molec%req(2)
            re3   = molec%req(3)
            ae1   = molec%alphaeq(1)
            ae2   = molec%alphaeq(2)
            !
            i0 = npoints/2
            !
            b0(3,2,i) = re2*sin(ae1)
            b0(3,1,i) = 0.0_ark
            b0(3,3,i) = re2*cos(ae1)
            !
            b0(1,2,i) = 0.0_ark
            b0(1,1,i) = 0.0_ark
            b0(1,3,i) = 0.0_ark
            !
            b0(2,2,i) = 0.0_ark
            b0(2,1,i) = 0.0_ark
            b0(2,3,i) = re1
            !
            b0(4,2,i) = re3*sin(ae2)*cos(rho)
            b0(4,1,i) = re3*sin(ae2)*sin(rho)
            b0(4,3,i) = re1-re3*cos(ae2)
            !
            ! Find center of mass
            !
            do n = 1,3 
              CM_shift = sum(b0(:,n,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
              b0(:,n,i) = b0(:,n,i) - CM_shift
            enddo 
            !
            !do n = 1,molec%Natoms
            !   b0(n,:,i) = matmul(transpose(transform),b0(n,:,i))
            !enddo
            !
         enddo
         !
         ! build the non-rigid reference structure point-by-point 
         ! and make sure it does not flip acex
         !
         do ipar = 0,1
           !
           istep = (-1)**(ipar+2)
           !
           i = npoints/2
           !
           a0(:,:) = b0(:,:,i)
           !
           call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,transform)
           !
           Inert0(1) = sum(molec%AtomMasses(:)*( a0(:,2)**2+ a0(:,3)**2) )
           Inert0(2) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,3)**2) )
           Inert0(3) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,2)**2) )
           !
           do ix = 1,molec%Natoms
              a0(ix,:) = matmul(transpose(transform),b0(ix,:,i))
           enddo
           !
           b0(:,:,i) = a0(:,:)
           !
           do j0 = 1,npoints/2
              !
              !write(out,'("j0 = ",i8)') j0
              !
              i = i + istep
              !
              do ix = 1,molec%Natoms
                 a0(ix,:) = matmul(transpose(transform),b0(ix,:,i))
              enddo
              !
              Inert(1) = sum(molec%AtomMasses(:)*( b0(:,2,i)**2+ b0(:,3,i)**2) )
              Inert(2) = sum(molec%AtomMasses(:)*( b0(:,1,i)**2+ b0(:,3,i)**2) )
              Inert(3) = sum(molec%AtomMasses(:)*( b0(:,1,i)**2+ b0(:,2,i)**2) )
              !
              ! Second Eckart equation
              ! 
              do ix = 1,3 
                 do jx = 1,3 
                    a(ix,jx) =  sum(molec%AtomMasses(:)*a0(:,ix)*a0(:,jx) )
                 enddo
                 !
              enddo
              !
              call MLorienting_a0(molec%Natoms,molec%AtomMasses,a0,a)
              !
              ! Found coordinate transformation "c"
              !
              transform = matmul(transform,real(a,kind=rk))
              !
              ! Transformation of a0 
              !
              do ix = 1,molec%Natoms
                 a0(ix,:) = matmul(transpose(transform),b0(ix,:,i))
              enddo
              !
              do ix = 1,3 
                 do jx = 1,3 
                    a(ix,jx) =  sum(molec%AtomMasses(:)*a0(:,ix)*a0(:,jx) )
                 enddo
                 !
              enddo
              !
              Inert(1) = sum(molec%AtomMasses(:)*( a0(:,2)**2+ a0(:,3)**2) )
              Inert(2) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,3)**2) )
              Inert(3) = sum(molec%AtomMasses(:)*( a0(:,1)**2+ a0(:,2)**2) )
              !
              Inert1 = 1
              !
              if (abs(i-i0)>2) then 
                !
                do in = 1,molec%Natoms
                   do ix = 1,3
                     !
                     N = min(4,abs(i-i0))
                     !
                     x(1:N+1) = rho_i(i-istep*N:i:istep)
                     !
                     call extrapolate(N,x,b0(in,ix,i-istep*N:i-istep:istep),a0_tt)
                     a0_t(in,ix) = a0_tt
                   enddo 
                enddo
                !
                Inert1(1) = sum(molec%AtomMasses(:)*( a0_t(:,2)**2+ a0_t(:,3)**2) )
                Inert1(2) = sum(molec%AtomMasses(:)*( a0_t(:,1)**2+ a0_t(:,3)**2) )
                Inert1(3) = sum(molec%AtomMasses(:)*( a0_t(:,1)**2+ a0_t(:,2)**2) )
                !
              endif 
              !
              if (all(abs(Inert(:)-Inert0(:))<5.0).or.all(abs(Inert(:)-Inert1(:))<5.0)) then 
                Inert0 = Inert
                b0(:,:,i) = a0(:,:)
                cycle 
              endif 
              !
              if (verbose>=4) then 
                write(out,"('b0',i4,12f12.8)") i,b0(:,:,i)
              endif
              !
           enddo
           !    
         enddo
         !
      endif
      !
      if (verbose>=4) then
        !
        do i = 0,npoints
           !
           write(out,"(i6)") molec%natoms
           !
           write(out,"(/a,3x,3f14.8)") trim(molec%zmatrix(1)%name),b0(1,:,i) 
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(2)%name),b0(2,:,i)
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(3)%name),b0(3,:,i)
           write(out,"( a,3x,3f14.8)") trim(molec%zmatrix(4)%name),b0(4,:,i)
           !
           !
        enddo 
        !
        !write(out,"('Coordinates:')")
        !
        !do i = 0,npoints
        !  write(out,"('b0',i4,12f12.8)") i,b0(:,:,i)
        !enddo
        !
      endif
      !
      if (verbose>=4) write(out,"('ML_b0_X2Y2/end')") 
      !
    contains 
     !
     ! extrapolation at borders
     !
     subroutine extrapolate(N,x,src,dst)
     !
     integer(ik),intent(in)   :: N
     real(ark),intent(inout)  :: src(1:N),x(1:N+1)
     real(ark),intent(out)    :: dst

     integer(ik)        :: i1,i2
     real(rk)           :: a(N,N),b(N,1)
        !
        !
        do i1 = 1,N
           !
           a(i1,1) = 1.0_ark
           !
           b(i1,1) = src(i1)
           !
           do i2 = 2,N
             !
             a(i1,i2) = x(i1)**(i2-1)
             !
           enddo
        enddo
        !
        !  lapack_gelss 
        ! 
        call lapack_gelss(a(:,:),b(:,:))
        !
        dst = b(1,1)
        !
        do i1 = 2,N
          !
          dst = dst + real(b(i1,1),ark)*x(N+1)**(i1-1)
          !
        enddo
        !
   end subroutine extrapolate
   !
  end subroutine ML_b0_X2Y2


  ! Here we define the coordinate transformation from the local coordinates (r,alpha) to 
  ! the internal coordinates chi used in the jacobian transformation and thereafter in the 
  ! as conjugate momenta coordinates
  !
  subroutine ML_symmetry_transformation_X2Y2(ioper,natoms,src,dst)
    !
    integer(ik),intent(in)    :: ioper  ! group operation  
    integer(ik),intent(in)    :: natoms
    real(ark),intent(in)      :: src(natoms)
    real(ark),intent(out)     :: dst(natoms)
    logical                   :: extended = .false.
    !
    integer(ik) :: nsrc,Nrot,N_Cn,ioper_,irot,NC2
    real(ark)   :: q1x,q2x,q1y,q2y,phi_n,phi,repres(sym%Noper,2,2)
    !
    if (verbose>=7) write(out,"('ML_symmetry_transformation_X2Y2/start')") 
    !
    nsrc = size(src)
    !
    if(trim(molec%symmetry)=='C'.or.trim(molec%symmetry)=='C(M)') then 
       dst = src
       return 
    endif
    !
    if (molec%rho_border(2)>2.0_ark*pi) extended = .true.
    !
    select case(trim(molec%coords_transform))
       !
       case default
       write (out,"('ML_symmetry_transformation_X2Y2: coord_transf ',a,' unknown')") trim(molec%coords_transform)
       stop 'ML_symmetry_transformation_X2Y2 - bad coord. type'
       !
    case('LINEAR','X-XE')
       !
       select case(trim(molec%symmetry))
       case default
          write (out,"('ML_symmetry_transformation_X2Y2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_X2Y2 - bad symm. type'
          !
       case('CS','CS(M)')
         !
         select case(ioper)
         !
         case (1) ! identity 
           !
           dst = src
           !
         case (2) ! (E*)
           !
           dst(1:5) = src(1:5)
           dst(6) =-src(6)
           !
         case default

           write (out,"('ML_symmetry_transformation_X2Y2: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_X2Y2 - bad operation. type'
 
         end select 
         !
       case('C2H','C2H(M)')
             !
         select case(ioper)

         case (1) ! E 

           dst = src

         case (2) ! (12)(34)

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)

         case (3) ! (E*)

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) =-src(6)
           !
         case (4) ! (12)(34)*

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) =-src(6)

         case default

           write (out,"('ML_symmetry_transformation_X2Y2: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_X2Y2 - bad operation. type'
 
         end select 
         !
       end select
        !
    case('R-ALPHA-TAU')
       !
       select case(trim(molec%symmetry))
       !
       case default
          write (out,"('ML_symmetry_transformation_X2Y2: symmetry ',a,' unknown')") trim(molec%symmetry)
          stop 'ML_symmetry_transformation_X2Y2 - bad symm. type'
          !
       case('C','C(M)')
         !
         dst = src
         !
       case('CS','CS(M)')
         !
         select case(ioper)
         !
         case (1) ! identity 

           dst = src

         case (2) ! (E*)

           dst = src
           dst(6) = 2.0_ark*pi-src(6)

         case default

           write (out,"('ML_symmetry_transformation_X2Y2: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_X2Y2 - bad operation. type'
 
         end select 
         !
       case('D2H(M)','G4(M)')
           !
         select case(ioper)
           !
         case (1) ! E 
           !
           dst = src
           !
         case (2) ! (C2z) -> Rz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi+src(6) 
           if (extended) then
              if(dst(6)>4.0_ark*pi) dst(6) = dst(6) - 4.0_ark*pi
           else
              dst(6) = src(6)
           endif
           !
         case (3) ! (C2x)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = src(6)
           !
         case (4) ! (C2y)
           !
           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi+src(6) ! ; if (dst(6)>0) dst(6) = dst(6) - 4.0_ark*pi
           if (extended) then 
              if(dst(6)>4.0_ark*pi) dst(6) = dst(6) - 4.0_ark*pi
           else
              dst(6) = src(6)
           endif
           !
         case (5) ! E* -> sigma_xz
           !
           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 4.0_ark*pi-src(6)
           if (extended) then 
              dst(6) = dst(6)
           else
              dst(6) = 2.0_ark*pi-src(6)
           endif

         case (6) ! sigma_yz

           dst(1) = src(1)
           dst(2) = src(2)
           dst(3) = src(3)
           dst(4) = src(4)
           dst(5) = src(5)
           dst(6) = 2.0_ark*pi-src(6) !; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           if (extended) then 
              if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           else
              dst(6) = dst(6)
           endif

         case (7) ! sigma_xy

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 4.0_ark*pi-src(6)
           if (extended) then 
              dst(6) = dst(6)
           else
              dst(6) = 2.0_ark*pi-src(6)
           endif

         case (8) ! i

           dst(1) = src(1)
           dst(2) = src(3)
           dst(3) = src(2)
           dst(4) = src(5)
           dst(5) = src(4)
           dst(6) = 2.0_ark*pi-src(6) !; if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           if (extended) then 
              if (dst(6)<0) dst(6) = dst(6) + 4.0_ark*pi
           else
              dst(6) = dst(6)
           endif

         case default

           write (out,"('ML_symmetry_transformation_X2Y2: operation ',i8,' unknown')") ioper
           stop 'ML_symmetry_transformation_X2Y2 - bad operation. type'
 
         end select 
         !
       end select
       !
    end select
    !
    if (verbose>=7) write(out,"('ML_symmetry_transformation_X2Y2/end')") 
    !
    !
  end subroutine ML_symmetry_transformation_X2Y2


  ! Here we find the rotational symmetry from J,K
  !
  subroutine ML_rotsymmetry_X2Y2(J,K,tau,gamma,ideg)
    !
    integer(ik),intent(in)  :: J,K,tau   ! rot. quanta
    integer(ik),intent(out) :: gamma,ideg
    !
    integer(ik) :: N,N_Cn,K_,L
    !
    if (verbose>=7) write(out,"('ML_rotsymmetry_X2Y2/start')") 
    !
    select case(trim(molec%symmetry))
    case default
       write (out,"('ML_rotsymmetry_X2Y2: symmetry ',a,' unknown')") trim(molec%symmetry)
       stop 'ML_rotsymmetry_X2Y2 - bad symm. type'
       !
    case('C','C(M)')
       !
       gamma = 1
       ideg = 1 
       !
    case('CS','CS(M)')
       !
       gamma = 0 
       ideg = 1 
       !
       if (mod(K+tau+2,2)==0) gamma = 1 !; return
       if (mod(K+tau+2,2)/=0) gamma = 2 !; return
       !
    case('C2H(M)','C2H')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma =1 !1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma =3 !3 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma =3 !4 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma =1 !2 !; return
       !
    case('D2H(M)','G4(M)')
       !
       gamma = 0 
       ideg = 1
       if (mod(K+2,2)==0.and.tau==0) gamma = 1 !1 !; return
       if (mod(K+2,2)==0.and.tau==1) gamma = 3 !3 !; return
       if (mod(K+2,2)/=0.and.tau==0) gamma = 5 !7 !; return
       if (mod(K+2,2)/=0.and.tau==1) gamma = 7 !5 !; return
       !
    case('DNH','DNH(M)')
       !
       gamma = 0 
       ideg = 1
       !
       N = sym%N
       N_Cn = sym%N/2
       k_ = mod(K+N_Cn,N_Cn)
       l = k_ ; if (k_>N_Cn) l = sym%N-k_
       !
       if (mod(sym%N,2)==1) then
          !
          if (mod(K+N_Cn,N_Cn)==0) then
             !
             if     (tau==0.and.mod(k+2,2)==0) then 
                gamma = 1 
             elseif (tau==1.and.mod(k+2,2)==0) then 
                gamma = 2
             elseif (tau==0.and.mod(k+2,2)/=0) then 
                gamma = 4
             elseif (tau==1.and.mod(k+2,2)/=0) then 
                gamma = 3
             else
                stop 'ML_rotsymmetry_X2Y2-Dnh: illegal k,tau (K mod N  = 0)'
             endif
             !
          elseif (tau<=1.and.k<=j) then
             !
             ideg = 1 ! tau +1
             if (mod(k+tau,2)/=0) ideg = 2
             !
             if     (mod(k+2,2)==0) then 
                 gamma = 4+2*l-1
             else
                gamma = 4+2*l
             endif
             !
          else
               stop 'ML_rotsymmetry_X2Y2-Dnh: illegal k,tau (K mod N  /= 0)'
          endif
          !
       else ! even Dnh
          !
          if (mod(K+N_Cn,N_Cn)==0) then
             !
             if     (tau==0.and.mod(k+2,2)==0) then 
                gamma = 1 
             elseif (tau==1.and.mod(k+2,2)==0) then 
                gamma = 2
             elseif (tau==0.and.mod(k+2,2)/=0.and.mod(N_Cn,2)/=0) then 
                gamma = 4
             elseif (tau==1.and.mod(k+2,2)/=0.and.mod(N_Cn,2)/=0) then 
                gamma = 3
             else
                stop 'ML_rotsymmetry_X2Y2-Dnh: illegal k,tau (K mod N  = 0)'
             endif
             !
          elseif (tau<=1.and.k<=j) then
             !
             !ideg = tau +1
             !
             ideg = 1
             !
             if (mod(k+tau,2)/=0) ideg = 2
             !
             gamma = 8+2*l-1
             !
             !if     (mod(k+2,2)==0) then 
             !    gamma = 8+2*l
             !else
             !    gamma = 8+2*l-1
             !endif
             !
          else
               stop 'ML_rotsymmetry_X2Y2-Dnh: illegal k,tau (K mod N  /= 0)'
          endif
          !
       endif
       !
    end select
    !
    if (verbose>=7) write(out,"('ML_rotsymmetry_X2Y2/end')") 
    !
    !
  end subroutine ML_rotsymmetry_X2Y2


end module mol_x2y2

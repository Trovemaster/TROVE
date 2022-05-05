!
module symmetry
  use accuracy
  use timer

  implicit none
  public SymmetryInitialize,sym


  type  ScIIT
     integer(ik)          :: Noper     ! Number of operations in the CII operator
     integer(ik)          :: Nzeta     ! Number of symmetric elements taking into account the degeneracies 
     integer(ik),pointer    :: ioper(:)  ! the operation number in the MS group
     integer(ik),pointer  :: coeff(:)  ! coefficients of the CII operator
     integer(ik),pointer  :: izeta(:)  ! symmetry indentification as a eigenvalues of the CII operator
  end type ScIIT


  type  SymmetryT
     character(len=cl)    :: group = 'C' ! The symmetry group 
     integer(ik)          :: Nrepresen = 1  ! Number of irreduc. represent.
     integer(ik)          :: Noper     = 1  ! Number of operations
     integer(ik)          :: Nclasses  = 1  ! Number of classes
     integer(ik),pointer  :: Nelements(:)   ! Number of elements in a class
     real(rk),pointer     :: characters(:,:)! Character table
     type(SrepresT),pointer :: irr(:,:)     ! irreducible representaion 
     integer(ik),pointer  :: degen(:)       ! degeneracy
     character(len=4),pointer  :: label(:)  ! The symmetry label 
     integer(ik)          :: Maxdegen  = 1  ! Maximal degeneracy order
     integer(ik),pointer  :: igenerator(:)  ! address of the class generator in the sym%Ngroup list
     type(ScIIT)          :: CII            ! the elements of the CII operator 
     real(ark),pointer    :: euler(:,:)     ! rotational angles equivalent to the group operations
     integer(ik)          :: class_size_max = 8 ! current maximal class size 
     integer(ik)          :: N  = 1         ! The group order, currently desgined for Dnh where N is odd 
     integer(ik),pointer  :: lquant(:)      ! Store the value of the (vibrational) angular momentum 
     integer(ik), allocatable :: product_table(:,:)        ! Stores information on obtaining all group elements from the generators
     logical :: product_table_set = .false. ! Whether or not the product table has been set  
     !
  end type SymmetryT

  type  SrepresT
     real(ark),pointer  :: repres(:,:)      ! matrix representation of the group 
  end type SrepresT


  type(SymmetryT) , save  :: sym
  integer(ik),parameter   :: max_irreps=100
  integer(ik),parameter   :: verbose_ = 3

contains 


  subroutine SymmetryInitialize(sym_group)
  character(len=cl),intent(inout) :: sym_group
  integer(ik):: alloc,iclass,gamma,ioper,ielem,irepr,Nrot,irep,k,irot,N_Cn,ioper_,icn,NC2,joper,jclass
  real(ark)  :: a,b,e,o,p2,p3,p4,p23,p43,phi,phi_n,factor,f_t,mat_t(2,2),repres_(2,2), m_one,mat_tt(4,4)
  character(len=4) :: Kchar
  !
  integer(ik),allocatable :: iclass_of(:)
  real(ark),allocatable :: characters_(:,:)
  
  real(ark), dimension(2 , 2) :: i, c, c2, sxy, sxy_, s2, s2_, s3, s3_
  real(ark), dimension(4, 4) :: gi, gcp, gc2p, gsp, gs2p, gs3p, gcm, gc2m, gsm
  real(ark), dimension(4, 4) :: gs2m, gs3m, g1,g2,g4,g7,g19 
  real(ark), dimension(6, 4, 4) :: G_rep_1, G_rep_2
  real(ark), dimension(6, 2, 2) :: E_rep_1, E_rep_2
  integer(ik), dimension(6) :: A2_char  
  integer(ik), dimension(6 , 6) :: pos_array
  integer(ik) :: j,r,s,n_c2v
  real :: log_n  
  integer(ik), dimension(144) :: tempG36  
  integer, dimension(2,2) :: z2_mat
  character(len = 5) :: k_num 
  character(len=3) :: sym_sub_label
  character(len=10) :: sym_cur_label
  !   
  !
  sym%group=sym_group
  !
  select case(trim(sym_group))

  case("C(M)","C")

    sym%Nrepresen=1
    sym%Noper=1
    sym%Nclasses=1
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters(1,1)=1.0_ark
    sym%degen=(/1/)
    sym%Nelements=(/1/)
    sym%label=(/'A'/)

    do gamma = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem = 1,sym%Nelements(iclass)
          !
          ioper = ioper + 1
          !
          allocate(sym%irr(gamma,ioper)%repres(sym%degen(gamma),sym%degen(gamma)),stat=alloc)
          !
          if (sym%degen(gamma)==1) then
             sym%irr(gamma,ioper)%repres(1,1)= sym%characters(gamma,iclass)
          endif 
          !
        enddo 
      enddo 
    enddo 
    !
  case("CS(M)","CS")
    !
    sym%Nrepresen=2
    sym%Noper=2
    sym%Nclasses=2
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                 (/ 1, 1,  & 
                                    1,-1  /),(/2,2/))
    sym%degen=(/1,1/)
    sym%Nelements=(/1,1/)
    sym%label=(/'A''','A"'/)
    !
    call irr_allocation
    !
  case("C2V(M)","C2V")
    !
    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                  (/ 1, 1, 1, 1, &   ! A1
                                     1, 1,-1,-1, &   ! A2
                                     1,-1,-1, 1, &   ! B1
                                     1,-1, 1,-1 /),(/4,4/)) ! B2
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'A1','A2','B1','B2'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o /)
    sym%euler( 3,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/p2,pi,p3/)
    !
    call irr_allocation
    !
  case("C2VN","C2V(N)","C2VN(M)")
    !
    sym_group = "C2VN"
    !
    log_n = log(4.0_ark*( real(sym%N,ark)+1.0_ark ))/log(2.0_ark)
    n_c2v = ceiling(log_n)

    sym%Nrepresen = 2**n_c2v
    sym%Noper = 2**n_c2v
    sym%Nclasses = 2**n_c2v
    sym%CII%Noper = 0

    call simple_arrays_allocation

    z2_mat = reshape( (/1, 1, &
                      1,-1/),(/2,2/))


    sym%characters(1:2,1:2) = z2_mat
    do j = 2, n_c2v
      sym%characters(1:2**(j-1), 2**(j-1)+1:2**j) = sym%characters(1:2**(j-1),1:2**(j-1))*z2_mat(2,1)
      sym%characters(2**(j-1)+1:2**j, 1:2**(j-1)) = sym%characters(1:2**(j-1),1:2**(j-1))*z2_mat(1,2)
      sym%characters(2**(j-1)+1:2**j, 2**(j-1)+1:2**j) = sym%characters(1:2**(j-1),1:2**(j-1))*z2_mat(2,2)
    enddo

    !do j=1, 2**n_c2v
    !  do k=1, 2**n_c2v
    !    !write(*,*) sym%characters(k,j)
    !  enddo
    !enddo
    
    do j=1, 2**n_c2v
      sym%degen(j) = 1
      sym%Nelements(j) = 1
      if(mod(j,4) == 1) then
        sym_sub_label = 'e'
      elseif(mod(j,4) == 2) then
        sym_sub_label = 'f' 
      elseif(mod(j,4) == 3) then
        sym_sub_label = 'g' 
      elseif(mod(j,4) == 0) then
        sym_sub_label = 'h'
      endif
      write(k_num,'(i5)') (j-1-mod(j-1,4))/4
      sym_cur_label = trim(sym_sub_label)//trim(adjustl(k_num))
      sym%label(j) =  sym_cur_label
      !write(*,*) sym%label(j)
    enddo
    !
    sym%label(1:4)=(/'A1','B2','A2','B1'/)
    !
    ! generators and the product table  
    !
    allocate(sym%product_table(sym%Noper,2),stat=alloc)
    call ArrayStart('sym%product_table',alloc,size(sym%product_table),kind(sym%product_table))
    !
    sym%product_table(1:4,1:2) = reshape((/0,0,0,0, &
                                           0,0,0,0/), (/4,2/))
    !
    do j = 1,2**(n_c2v-2)-1
      sym%product_table(4*j+1: 4*j+4,1:2) = reshape((/1,2,3,4, &
                                                      1,1,1,1/), (/4,2/))
    enddo
    !
    sym%product_table_set = .true.

    call irr_allocation
    !
  case("C2H(M)","C2H")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                  (/ 1, 1, 1, 1, &   ! Ag
                                     1, 1,-1,-1, &   ! Au
                                     1,-1,-1, 1, &   ! Bg
                                     1,-1, 1,-1 /),(/4,4/)) ! Bu
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'Ag','Au','Bg','Bu'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 2,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    !
    ! 1324 and 1234 are working
    !
    call irr_allocation
    !
  case("CS(EM)")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2a  E* sab i
                                  (/ 1,  1,  1,  1, &   ! Ag 
                                     1,  1, -1, -1, &   ! Au 
                                     1, -1,  1, -1, &   ! Bg
                                     1, -1, -1,  1/),(/4 ,4/)) ! Bu
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'Ag','Au','Bg','Bu'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 4,:) = (/o,pi,o/)
    !
    call irr_allocation
    !
  case("G4(M)","G4")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &       ! E (12)(34)  E*  (12)(34)*
                                  (/ 1, 1, 1, 1, &   ! A+
                                     1, 1,-1,-1, &   ! A-
                                     1,-1,-1, 1, &   ! B+
                                     1,-1, 1,-1 /),(/4,4/)) ! B-
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'A+','A-','B-','B+'/)
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p2,pi,p3/)
    sym%euler( 3,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    !
    call irr_allocation
    !
  case("G4(EM)")

    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E   a   b  ab   E' E'a E'b E'ab 
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ags
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! Aus
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! Bgs
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! Bus
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! Agd 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Aud
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! Bgd
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! Bud
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ags','Aus','Bgs','Bus','Agd','Aud','Bgd','Bud'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 2,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    sym%euler( 5,:) = (/pi,o,o/)
    sym%euler( 7,:) = (/o,pi,o/)
    sym%euler( 6,:) = (/p2,pi,p3/)
    sym%euler( 8,:) = 0


    !sym%euler( 1,:) = 0
    !sym%euler( 3,:) = (/p2,pi,p3/)
    !sym%euler( 4,:) = (/o,pi,o/)
    !sym%euler( 2,:) = (/pi,o,o/)
    !sym%euler( 7,:) = (/pi,o,o/)
    !sym%euler( 5,:) = (/o,pi,o/)
    !sym%euler( 6,:) = (/p2,pi,p3/)
    !sym%euler( 8,:) = 0


    !
    call irr_allocation

    case("G36(M)", "G36")
      
      write(*,*) "test 1"
      sym%Nrepresen = 9
      sym%Noper = 36
      sym%Nclasses = 9
      sym%CII%Noper = 0
                               
      call simple_arrays_allocation
    
      !RE = 2C2+E-   SE = 3C3+E-  ER = 2E+C2-  RR = 4C2+C2- SR = 6C3+C2- ER = 3E+C3- RS = 6C2+C3- SS = 9C3+C3-
      sym%characters = reshape( & 
        ! EE  RE  SE  ER  RR  SR  ER  RS  SS  
        (/ 1,  1,  1,  1,  1,  1,  1,  1,  1, & ! A1s
           1,  1,  1,  1,  1,  1, -1, -1, -1, & ! A2s
           1,  1, -1,  1,  1, -1,  1,  1, -1, & ! A3s
           1,  1, -1,  1,  1, -1, -1, -1,  1, & ! A4s
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0, & ! E1s
           2,  2, -2, -1, -1,  1,  0,  0,  0, & ! E2s
           2, -1,  0,  2, -1,  0,  2, -1,  0, & ! E3s
           2, -1,  0,  2, -1,  0, -2,  1,  0, & ! E4s
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0 /),(/9,9/)) ! Gs  i
      !
      sym%characters = transpose(sym%characters)
      sym%degen = (/1, 1, 1, 1, 2, 2, 2, 2, 4/)
      sym%Nelements = (/1, 2, 3, 2, 4, 6, 3, 6, 9/)
      sym%label=(/'A1', 'A2', 'A3', 'A4', 'E1', 'E2', 'E3', 'E4', 'G ' /)
      a = 0.5_ark
      b = 0.5_ark*sqrt(3.0_ark)
      e = 1.0_ark
      o = 0.0_ark
      m_one = -1.0_ark       
 
      !E rep of C_3v 

      i = transpose(reshape( (/ e, o, &
                                o, e /), (/ 2, 2/)))
  
      c = transpose(reshape( (/ -a, -b, &
                                 b, -a/), (/ 2, 2/)))
  
    c2 = matmul(c,c)
  
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))
  
    s3 = matmul(c, sxy)

    s2 = matmul(c,s3)
  
    
 
    E_rep_1(1,:,:) = i
    E_rep_1(2,:,:) = c
    E_rep_1(3,:,:) = c2
    E_rep_1(4,:,:) = sxy
    E_rep_1(5,:,:) = s2
    E_rep_1(6,:,:) = s3
    
    A2_char(:3) = 1
    A2_char(4:6) = -1
  
    call irr_allocation 
  
    pos_array = transpose( reshape((/  1,  7,  8, 19, 20, 21, &
                                       2,  9, 11, 22, 24, 26, &
                                       3, 10, 12, 23, 25, 27, &
                                       4, 13, 16, 28, 31, 34, &
                                       5, 14, 17, 29, 32, 35, &
                                       6, 15, 18, 30, 33, 36/), (/6,6/))) 
    ! E1 and E2
 
    do j = 1, 6
      do k = 1, 6 
        sym%irr(5, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(6, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
      end do
    end do

    !E3s and E4s
    
    do j = 1, 6
      do k = 1, 6
        sym%irr(7,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(8,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
      end do    
    end do
    
    ! Gs and Gd

    ! Implements outer product,
    ! consider two matrices
    !   a b    and e f
    !   c d        g h
    ! then the outer product is 
    !   ae af be bf 
    !   ag ah bg bh
    !   ce cf de df
    !   cg ch dg dh

    do j = 1, 6
      do k = 1, 6
        do r = 1, 2
          do s = 1, 2
            sym%irr(9,pos_array(k,j))%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)
          end do
        end do
      end do 
    end do
    write(*,*) "test 1"


    !
    ! sy16: Per's matrices, did not improve 
    ! sy24: as part of sap 7 and 8 
    !
    case("G36(EMI)") 

      sym%Nrepresen = 18
      sym%Noper = 72
      sym%Nclasses = 18
      sym%CII%Noper = 0
                               
      call simple_arrays_allocation
      !
      tempG36(1:36)   =(/0,0,2,0,6,4,0,7,2,3,2,3,4,5,6,4,5,6,0,21,19, 2, 3, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6/)
      tempG36(73:108) =(/0,0,2,0,2,2,0,7,7,7,8,8,7,7,7,8,8,8,0, 7, 7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
      tempG36(37:72)  =(/0,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,&
                         37,37,37,37,37,37/)
      tempG36(109:144)=(/0, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,&
                         31,32,33,34,35,36/)
      !
      allocate(sym%product_table(72,2),stat=alloc)
      call ArrayStart('sym%product_table',alloc,size(sym%product_table),kind(sym%product_table))
      !
      sym%product_table = reshape( tempG36, (/ 72, 2/))
      sym%product_table_set = .true.
      !
      !RE = 2C2+E-   SE = 3C3+E-  ER = 2E+C2-  RR = 4C2+C2- SR = 6C3+C2- ER = 3E+C3- RS = 6C2+C3- SS = 9C3+C3-
      sym%characters = reshape( & 
        ! EE  RE  SE  ER  RR  SR  ER  RS  SS EE' RE' SE' ER' RR' SR' ER' RS' SS' 
        (/ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! A1s
           1,  1,  1,  1,  1,  1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, & ! A2s
           1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1, & ! A3s
           1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1, -1,  1,  1, -1, -1, -1,  1, & ! A4s
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0,  2,  2,  2, -1, -1, -1,  0,  0,  0, & ! E1s
           2,  2, -2, -1, -1,  1,  0,  0,  0,  2,  2, -2, -1, -1,  1,  0,  0,  0, & ! E2s
           2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0, & ! E3s
           2, -1,  0,  2, -1,  0, -2,  1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, & ! E4s
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0,  4, -2,  0, -2,  1,  0,  0,  0,  0, & ! Gs
           !
           1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, & ! A1d
           1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, & ! A2d 
           1,  1, -1,  1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1, & ! A3d
           1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1,  1,  1, -1, & ! A4d
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0, -2, -2, -2,  1,  1,  1,  0,  0,  0, & ! E1d
           2,  2, -2, -1, -1,  1,  0,  0,  0, -2, -2,  2,  1,  1, -1,  0,  0,  0, & ! E2d
           2, -1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0, & ! E3d
           2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0,  2, -1,  0, & ! E4d
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0, -4,  2,  0,  2, -1,  0,  0,  0,  0 /),(/18,18/)) ! Gd
      sym%characters = transpose(sym%characters)
      sym%degen = (/1, 1, 1, 1, 2, 2, 2, 2, 4, 1, 1, 1, 1, 2, 2, 2, 2, 4/)
      sym%Nelements = (/1, 2, 3, 2, 4, 6, 3, 6, 9, 1, 2, 3, 2, 4, 6, 3, 6, 9 /)
      sym%label=(/'A1s', 'A2s', 'A3s', 'A4s', 'E1s', 'E2s', 'E3s', 'E4s', 'Gs ', &
                  'A1d', 'A2d', 'A3d', 'A4d', 'E1d', 'E2d', 'E3d', 'E4d', 'Gd ' /)

      p2 = 0.5_ark*pi
      p3 = 1.5_ark*pi
      o  = 0.0_ark
      !
      sym%euler( 1,:) = 0
      sym%euler( 2,:) = (/o,o,o/)
      sym%euler( 4,:) = (/o,o,pi/)
      sym%euler( 7,:) = (/o,o,2.0_ark/3.0_ark*pi/)
      sym%euler( 19,:) = (/p2,pi,p3/)
      !
      a = 0.5_ark
      b = 0.5_ark*sqrt(3.0_ark)
      e = 1.0_ark
      o = 0.0_ark
      m_one = -1.0_ark       
 
      i = transpose(reshape( (/ e, o, &
                                o, e /), (/ 2, 2/)))
  
      c = transpose(reshape( (/ -a, -b, &
                                 b, -a/), (/ 2, 2/)))
  
    c2 = matmul(c,c)
  
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))

    !sxy = transpose(reshape( (/ o, e, &
    !                            e, o /), (/ 2, 2/)))
    !
    s3 = matmul(c, sxy)
    !
    s2 = matmul(c,s3)
    ! 
    E_rep_1(1,:,:) = i
    E_rep_1(2,:,:) = c
    E_rep_1(3,:,:) = c2
    E_rep_1(4,:,:) = sxy
    E_rep_1(5,:,:) = s2
    E_rep_1(6,:,:) = s3
    !
    A2_char(:3) = 1.0_ark
    A2_char(4:6) = -1.0_ark
    !
    call irr_allocation 
    !
    pos_array = transpose( reshape((/  1,  7,  8, 19, 20, 21, &
                                       2,  9, 11, 22, 24, 26, &
                                       3, 10, 12, 23, 25, 27, &
                                       4, 13, 16, 28, 31, 34, &
                                       5, 14, 17, 29, 32, 35, &
                                       6, 15, 18, 30, 33, 36/), (/6,6/))) 
    ! E1s, E2s, E1d, and E2d 
    !
    do j = 1, 6
      do k = 1, 6 
        sym%irr(5, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(5, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)
        
        sym%irr(6, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
        sym%irr(6, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*A2_char(k)
       
        sym%irr(5+9, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(5+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*m_one
        
        sym%irr(6+9, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
        sym%irr(6+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*A2_char(k)*m_one
      end do
    end do
    !
    !E3s and E4s
    !
    do j = 1, 6
      do k = 1, 6
        sym%irr(7,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(7,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)

        sym%irr(8,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
        sym%irr(8,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)

        sym%irr(7+9,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(7+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*m_one

        sym%irr(8+9,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
        sym%irr(8+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)*m_one

      end do    
    end do
  
    g1 = transpose(reshape( (/ e, o, o ,o, &
                               o, e, o, o, &
                               o, o, e, o, &
                               o, o, o, e/), (/4,4/)))
    !
    ! q01 try this 
    ! u07 
    ! b08
    g2 = transpose(reshape( (/ -a,-b, o, o, &
                                b,-a, o, o, &
                                o, o,-a,-b, &
                                o, o, b,-a  /), (/4,4/))) 

    !c01
    !
    !g2 = transpose(reshape( (/ -a,-b, o, o, &
    !                            b,-a, o, o, &
    !                            o, o,-a, b, &
    !                            o, o,-b,-a  /), (/4,4/))) 

    !!
    !g2 = transpose(reshape( (/ -a, b, o, o, &
    !                           -b,-a, o, o, &
    !                            o, o,-a, b, &
    !                            o, o,-b,-a  /), (/4,4/))) 

    !t02
    !g4 = transpose(reshape( (/o, o,-e, o, &
    !                           o, o, o, e, &
    !                          -e, o, o, o, &
    !                           o, e, o, o  /), (/4,4/)))
    !!
    g4 = transpose(reshape( (/ e, o, o, o, &
                               o,-e, o, o, &
                               o, o,-e, o, &
                               o, o, o, e  /), (/4,4/)))
    !
    ! t03
    !g4 = transpose(reshape( (/-e, o, o, o, &
    !                           o, e, o, o, &
    !                           o, o, e, o, &
    !                           o, o, o,-e  /), (/4,4/)))

    ! a10
    g7 = transpose(reshape( (/ -a, o, o,-b, &
                                o,-a, b, o, &
                                o,-b,-a, o, &
                                b, o, o,-a  /), (/4,4/))) 
    ! c02
    !g7 = transpose(reshape( (/ -a, o, o,-b, &
    !                            o,-a,-b, o, &
    !                            o, b,-a, o, &
    !                            b, o, o,-a  /), (/4,4/))) 
    !!
    !g7 = transpose(reshape( (/ -a, o, o, b, &
    !                            o,-a,-b, o, &
    !                            o, b,-a, o, &
    !                           -b, o, o,-a  /), (/4,4/))) 

    ! sy5: try this instd of !! wrong character 
    !g7 = transpose(reshape( (/-a, b, o, o, &
    !                           -b,-a, o, o, &
    !                            o, o,-a,-b, &
    !                            o, o, b,-a  /), (/4,4/))) 
    ! t06 
    !g19= transpose(reshape( (/-e, o, o, o, &
    !                           o, e, o, o, &
    !                           o, o,-e, o, &
    !                           o, o, o, e  /), (/4,4/)))

    !!
    g19= transpose(reshape( (/ e, o, o, o, &
                               o, e, o, o, &
                               o, o,-e, o, &
                               o, o, o,-e  /), (/4,4/)))

    ! sy7: try this for !! completely wrong 
    !g19= transpose(reshape( (/o, o, e, o, &
    !                           o, o, o, e, &
    !                           e, o, o, o, &
    !                           o, e, o, o  /), (/4,4/)))


    i = transpose(reshape( (/ e, o, &
                              o, e /), (/ 2, 2/)))
    !!
    c = transpose(reshape( (/ -a, -b, &
                               b, -a/), (/ 2, 2/)))
    !
    !t08
    !u09
    ! b07
    ! b16: still same problem 
    !c = transpose(reshape( (/ -a,  b, &
    !                          -b, -a/), (/ 2, 2/)))

  
    c2 = matmul(c,c)
    !!
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))
    !
    ! t07 
    ! u10
    ! b10
    !sxy = transpose(reshape( (/ o, e, &
    !                            e, o /), (/ 2, 2/)))
    !
    s3 = matmul(c, sxy)
    !
    s2 = matmul(c,s3)
    !
    sym%irr( 5, 1)%repres = i
    sym%irr( 6, 1)%repres = i
    sym%irr( 7, 1)%repres = i
    sym%irr( 8, 1)%repres = i
    !
    sym%irr( 5, 2)%repres = i
    sym%irr( 6, 2)%repres = i
    sym%irr( 7, 2)%repres = c
    sym%irr( 8, 2)%repres = c
    !!
    sym%irr( 5, 4)%repres = i
    sym%irr( 6, 4)%repres = i*m_one
    sym%irr( 7, 4)%repres = sxy
    sym%irr( 8, 4)%repres = sxy
    !
    !u11
    !sym%irr( 5, 4)%repres = i
    !sym%irr( 6, 4)%repres = i
    !sym%irr( 7, 4)%repres = sxy
    !sym%irr( 8, 4)%repres = sxy
    !
    ! u12
    !
    !sym%irr( 5, 4)%repres = i
    !sym%irr( 6, 4)%repres = i*m_one
    !sym%irr( 7, 4)%repres = sxy
    !sym%irr( 8, 4)%repres = sxy*m_one
    !
    !!
    sym%irr( 5, 7)%repres = c2
    sym%irr( 6, 7)%repres = c2
    sym%irr( 7, 7)%repres = i
    sym%irr( 8, 7)%repres = i
    !
    ! a09 
    !
    !sym%irr( 5, 7)%repres = c
    !sym%irr( 6, 7)%repres = c
    !sym%irr( 7, 7)%repres = i
    !sym%irr( 8, 7)%repres = i
    !
    ! p06 try changing sign for oper 6  did not work 
    !
    !sym%irr( 5, 7)%repres = c2
    !sym%irr( 6, 7)%repres = c2*m_one
    !sym%irr( 7, 7)%repres = i
    !sym%irr( 8, 7)%repres = i
    !!
    sym%irr( 5,19)%repres = sxy
    sym%irr( 6,19)%repres = sxy
    sym%irr( 7,19)%repres = i
    sym%irr( 8,19)%repres = i*m_one
    !
    ! b15
    !sym%irr( 5,19)%repres = sxy
    !sym%irr( 6,19)%repres = sxy
    !sym%irr( 7,19)%repres = i
    !sym%irr( 8,19)%repres = i
    !
    ! p07 change E2 sign according with character of A2, did not help
    ! p08 the same but with sxy from alternative choice , did not work
    !sym%irr( 5,19)%repres = sxy
    !sym%irr( 6,19)%repres = sxy*m_one
    !sym%irr( 7,19)%repres = i
    !sym%irr( 8,19)%repres = i*m_one
    ! p10 another test, did not work 
    !sym%irr( 5,19)%repres = sxy*m_one
    !sym%irr( 6,19)%repres = sxy
    !sym%irr( 7,19)%repres = i
    !sym%irr( 8,19)%repres = i*m_one
    !
    sym%irr( 5,37)%repres = i
    sym%irr( 6,37)%repres = i
    sym%irr( 7,37)%repres = i
    sym%irr( 8,37)%repres = i
    !
    do irep = 5,8
      sym%irr(irep+9, 1)%repres = sym%irr(irep, 1)%repres
      sym%irr(irep+9, 2)%repres = sym%irr(irep, 2)%repres
      sym%irr(irep+9, 4)%repres = sym%irr(irep, 4)%repres
      sym%irr(irep+9, 7)%repres = sym%irr(irep, 7)%repres
      sym%irr(irep+9,19)%repres = sym%irr(irep,19)%repres
    enddo
    !
    sym%irr(5+9,37)%repres = i*m_one
    sym%irr(6+9,37)%repres = i*m_one
    sym%irr(7+9,37)%repres = i*m_one
    sym%irr(8+9,37)%repres = i*m_one
    !
    sym%irr( 9, 1)%repres = g1
    sym%irr( 9, 2)%repres = g2
    sym%irr( 9, 4)%repres = g4
    sym%irr( 9, 7)%repres = g7
    sym%irr( 9,19)%repres = g19
    !
    sym%irr(18, 1)%repres = g1
    sym%irr(18, 2)%repres = g2
    sym%irr(18, 4)%repres = g4
    sym%irr(18, 7)%repres = g7
    sym%irr(18,19)%repres = g19
    !
    sym%irr( 9, 1+36)%repres = g1
    sym%irr(18, 1+36)%repres = g1*m_one
    !
    do irep = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem =1,sym%Nelements(iclass)
          ioper = ioper + 1
          !
          call do_g36_transform(irep,ioper,sym%degen(irep),sym%irr(irep,ioper)%repres)
          !
          f_t = 0
          do k = 1,sym%degen(irep)
              f_t = f_t + (sym%irr(irep,ioper)%repres(k,k))
          enddo
          !
          !sym%characters(irep,iclass) = f_t
          !
        enddo
      enddo
    enddo
    !

    !
    case("G36(EM1)") 

      sym%Nrepresen = 18
      sym%Noper = 72
      sym%Nclasses = 18
      sym%CII%Noper = 0
                               
      call simple_arrays_allocation
      !
      tempG36(1:36)   = (/0,0,2,0,6,4,0,7,2,3,2,3,4,5,6,4,5,6,0,21,19, 2, 3, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6/)
      tempG36(73:108) = (/0,0,2,0,2,2,0,7,7,7,8,8,7,7,7,8,8,8,0, 7, 7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
      tempG36(37:72)  = (/0,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,&
                          37,37,37,37,37,37/)
      tempG36(109:144)= (/0, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,&
                         31,32,33,34,35,36/)
      !
      allocate(sym%product_table(72,2),stat=alloc)
      call ArrayStart('sym%product_table',alloc,size(sym%product_table),kind(sym%product_table))
      !
      sym%product_table = reshape( tempG36, (/ 72, 2/))
      sym%product_table_set = .true.
      !
      !RE = 2C2+E-   SE = 3C3+E-  ER = 2E+C2-  RR = 4C2+C2- SR = 6C3+C2- ER = 3E+C3- RS = 6C2+C3- SS = 9C3+C3-
      sym%characters = reshape( & 
        ! EE  RE  SE  ER  RR  SR  ER  RS  SS EE' RE' SE' ER' RR' SR' ER' RS' SS' 
        (/ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! A1s
           1,  1,  1,  1,  1,  1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, & ! A2s
           1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1, & ! A3s
           1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1, -1,  1,  1, -1, -1, -1,  1, & ! A4s
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0,  2,  2,  2, -1, -1, -1,  0,  0,  0, & ! E1s
           2,  2, -2, -1, -1,  1,  0,  0,  0,  2,  2, -2, -1, -1,  1,  0,  0,  0, & ! E2s
           2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0, & ! E3s
           2, -1,  0,  2, -1,  0, -2,  1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, & ! E4s
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0,  4, -2,  0, -2,  1,  0,  0,  0,  0, & ! Gs
           !
           1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, & ! A1d
           1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, & ! A2d 
           1,  1, -1,  1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1, & ! A3d
           1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1,  1,  1, -1, & ! A4d
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0, -2, -2, -2,  1,  1,  1,  0,  0,  0, & ! E1d
           2,  2, -2, -1, -1,  1,  0,  0,  0, -2, -2,  2,  1,  1, -1,  0,  0,  0, & ! E2d
           2, -1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0, & ! E3d
           2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0,  2, -1,  0, & ! E4d
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0, -4,  2,  0,  2, -1,  0,  0,  0,  0 /),(/18,18/)) ! Gd 
      sym%characters = transpose(sym%characters)
      sym%degen = (/1, 1, 1, 1, 2, 2, 2, 2, 4, 1, 1, 1, 1, 2, 2, 2, 2, 4/)
      sym%Nelements = (/1, 2, 3, 2, 4, 6, 3, 6, 9, 1, 2, 3, 2, 4, 6, 3, 6, 9 /)
      sym%label=(/'A1s', 'A2s', 'A3s', 'A4s', 'E1s', 'E2s', 'E3s', 'E4s', 'Gs ', &
                  'A1d', 'A2d', 'A3d', 'A4d', 'E1d', 'E2d', 'E3d', 'E4d', 'Gd ' /)

      p2 = 0.5_ark*pi
      p3 = 1.5_ark*pi
      o  = 0.0_ark
      !
      sym%euler( 1,:) = 0
      sym%euler( 2,:) = (/o,o,o/)
      sym%euler( 4,:) = (/o,o,pi/)
      sym%euler( 7,:) = (/o,o,2.0_ark/3.0_ark*pi/)
      sym%euler( 19,:) = (/p2,pi,p3/)
      !
      a = 0.5_ark
      b = 0.5_ark*sqrt(3.0_ark)
      e = 1.0_ark
      o = 0.0_ark
      m_one = -1.0_ark       
 
      i = transpose(reshape( (/ e, o, &
                                o, e /), (/ 2, 2/)))
  
      c = transpose(reshape( (/ -a, -b, &
                                 b, -a/), (/ 2, 2/)))
  
    c2 = matmul(c,c)
  
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))

    !sxy = transpose(reshape( (/ o, e, &
    !                            e, o /), (/ 2, 2/)))
    !
    s3 = matmul(c, sxy)
    !
    s2 = matmul(c,s3)
    ! 
    E_rep_1(1,:,:) = i
    E_rep_1(2,:,:) = c
    E_rep_1(3,:,:) = c2
    E_rep_1(4,:,:) = sxy
    E_rep_1(5,:,:) = s2
    E_rep_1(6,:,:) = s3
    !
    A2_char(:3) = 1.0_ark
    A2_char(4:6) = -1.0_ark
    !
    call irr_allocation 
    !
    pos_array = transpose( reshape((/  1,  7,  8, 19, 20, 21, &
                                       2,  9, 11, 22, 24, 26, &
                                       3, 10, 12, 23, 25, 27, &
                                       4, 13, 16, 28, 31, 34, &
                                       5, 14, 17, 29, 32, 35, &
                                       6, 15, 18, 30, 33, 36/), (/6,6/))) 
    ! E1s, E2s, E1d, and E2d 
    !
    do j = 1, 6
      do k = 1, 6 
        sym%irr(5, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(5, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)
        
        sym%irr(6, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
        sym%irr(6, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*A2_char(k)
       
        sym%irr(5+9, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(5+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*m_one
        
        sym%irr(6+9, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
        sym%irr(6+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*A2_char(k)*m_one
      end do
    end do
    !
    !E3s and E4s
    !
    do j = 1, 6
      do k = 1, 6
        sym%irr(7,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(7,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)

        sym%irr(8,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
        sym%irr(8,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)

        sym%irr(7+9,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(7+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*m_one

        sym%irr(8+9,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
        sym%irr(8+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)*m_one

      end do    
    end do
  
    g1 = transpose(reshape( (/ e, o, o ,o, &
                               o, e, o, o, &
                               o, o, e, o, &
                               o, o, o, e/), (/4,4/)))
    !
    ! q01 try this 
    ! u07 
    ! b08
    !g2 = transpose(reshape( (/ -a,-b, o, o, &
    !                            b,-a, o, o, &
    !                            o, o,-a,-b, &
    !                            o, o,-b,-a  /), (/4,4/))) 

    ! y01
    g7 = transpose(reshape( (/ -a, o, o,-b, &
                                o,-a, b, o, &
                                o,-b,-a, o, &
                                b, o, o,-a  /), (/4,4/))) 




    !c01
    !
    !g2 = transpose(reshape( (/ -a,-b, o, o, &
    !                            b,-a, o, o, &
    !                            o, o,-a, b, &
    !                            o, o,-b,-a  /), (/4,4/))) 

    !!
    !g2 = transpose(reshape( (/ -a, b, o, o, &
    !                           -b,-a, o, o, &
    !                            o, o,-a, b, &
    !                            o, o,-b,-a  /), (/4,4/))) 

    !t02
    !g4 = transpose(reshape( (/o, o,-e, o, &
    !                           o, o, o, e, &
    !                          -e, o, o, o, &
    !                           o, e, o, o  /), (/4,4/)))
    !!
    g4 = transpose(reshape( (/ e, o, o, o, &
                               o,-e, o, o, &
                               o, o,-e, o, &
                               o, o, o, e  /), (/4,4/)))
    !
    ! t03
    !g4 = transpose(reshape( (/-e, o, o, o, &
    !                           o, e, o, o, &
    !                           o, o, e, o, &
    !                           o, o, o,-e  /), (/4,4/)))

    ! a10
    !g7 = transpose(reshape( (/ -a, o, o,-b, &
    !                            o,-a, b, o, &
    !                            o,-b,-a, o, &
    !                            b, o, o,-a  /), (/4,4/))) 

    !y01
    g7 = transpose(reshape( (/ -a,-b, o, o, &
                                b,-a, o, o, &
                                o, o,-a,-b, &
                                o, o,-b,-a  /), (/4,4/))) 


    ! c02
    !g7 = transpose(reshape( (/ -a, o, o,-b, &
    !                            o,-a,-b, o, &
    !                            o, b,-a, o, &
    !                            b, o, o,-a  /), (/4,4/))) 
    !!
    !g7 = transpose(reshape( (/ -a, o, o, b, &
    !                            o,-a,-b, o, &
    !                            o, b,-a, o, &
    !                           -b, o, o,-a  /), (/4,4/))) 

    ! sy5: try this instd of !! wrong character 
    !g7 = transpose(reshape( (/-a, b, o, o, &
    !                           -b,-a, o, o, &
    !                            o, o,-a,-b, &
    !                            o, o, b,-a  /), (/4,4/))) 
    ! t06 
    !g19= transpose(reshape( (/-e, o, o, o, &
    !                           o, e, o, o, &
    !                           o, o,-e, o, &
    !                           o, o, o, e  /), (/4,4/)))

    !!
    g19= transpose(reshape( (/ e, o, o, o, &
                               o, e, o, o, &
                               o, o,-e, o, &
                               o, o, o,-e  /), (/4,4/)))

    ! sy7: try this for !! completely wrong 
    !g19= transpose(reshape( (/o, o, e, o, &
    !                           o, o, o, e, &
    !                           e, o, o, o, &
    !                           o, e, o, o  /), (/4,4/)))


    i = transpose(reshape( (/ e, o, &
                              o, e /), (/ 2, 2/)))
    !!
    c = transpose(reshape( (/ -a, -b, &
                               b, -a/), (/ 2, 2/)))
    !
    !t08
    !u09
    ! b07
    ! b16: still same problem 
    !c = transpose(reshape( (/ -a,  b, &
    !                          -b, -a/), (/ 2, 2/)))

  
    c2 = matmul(c,c)
    !!
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))
    !
    ! t07 
    ! u10
    ! b10
    !sxy = transpose(reshape( (/ o, e, &
    !                            e, o /), (/ 2, 2/)))
    !
    s3 = matmul(c, sxy)
    !
    s2 = matmul(c,s3)
    !
    sym%irr( 5, 1)%repres = i
    sym%irr( 6, 1)%repres = i
    sym%irr( 7, 1)%repres = i
    sym%irr( 8, 1)%repres = i
    !
    sym%irr( 5, 2)%repres = i
    sym%irr( 6, 2)%repres = i
    sym%irr( 7, 2)%repres = c
    sym%irr( 8, 2)%repres = c
    !!
    sym%irr( 5, 4)%repres = i
    sym%irr( 6, 4)%repres = i*m_one
    sym%irr( 7, 4)%repres = sxy
    sym%irr( 8, 4)%repres = sxy
    !
    !u11
    !sym%irr( 5, 4)%repres = i
    !sym%irr( 6, 4)%repres = i
    !sym%irr( 7, 4)%repres = sxy
    !sym%irr( 8, 4)%repres = sxy
    !
    ! u12
    !
    !sym%irr( 5, 4)%repres = i
    !sym%irr( 6, 4)%repres = i*m_one
    !sym%irr( 7, 4)%repres = sxy
    !sym%irr( 8, 4)%repres = sxy*m_one
    !
    !!
    sym%irr( 5, 7)%repres = c2
    sym%irr( 6, 7)%repres = c2
    sym%irr( 7, 7)%repres = i
    sym%irr( 8, 7)%repres = i
    !
    ! a09 
    !
    !sym%irr( 5, 7)%repres = c
    !sym%irr( 6, 7)%repres = c
    !sym%irr( 7, 7)%repres = i
    !sym%irr( 8, 7)%repres = i
    !
    ! p06 try changing sign for oper 6  did not work 
    !
    !sym%irr( 5, 7)%repres = c2
    !sym%irr( 6, 7)%repres = c2*m_one
    !sym%irr( 7, 7)%repres = i
    !sym%irr( 8, 7)%repres = i
    !!
    sym%irr( 5,19)%repres = sxy
    sym%irr( 6,19)%repres = sxy
    sym%irr( 7,19)%repres = i
    sym%irr( 8,19)%repres = i*m_one
    !
    ! b15
    !sym%irr( 5,19)%repres = sxy
    !sym%irr( 6,19)%repres = sxy
    !sym%irr( 7,19)%repres = i
    !sym%irr( 8,19)%repres = i
    !
    ! p07 change E2 sign according with character of A2, did not help
    ! p08 the same but with sxy from alternative choice , did not work
    !sym%irr( 5,19)%repres = sxy
    !sym%irr( 6,19)%repres = sxy*m_one
    !sym%irr( 7,19)%repres = i
    !sym%irr( 8,19)%repres = i*m_one
    ! p10 another test, did not work 
    !sym%irr( 5,19)%repres = sxy*m_one
    !sym%irr( 6,19)%repres = sxy
    !sym%irr( 7,19)%repres = i
    !sym%irr( 8,19)%repres = i*m_one
    !
    sym%irr( 5,37)%repres = i
    sym%irr( 6,37)%repres = i
    sym%irr( 7,37)%repres = i
    sym%irr( 8,37)%repres = i
    !
    do irep = 5,8
      sym%irr(irep+9, 1)%repres = sym%irr(irep, 1)%repres
      sym%irr(irep+9, 2)%repres = sym%irr(irep, 2)%repres
      sym%irr(irep+9, 4)%repres = sym%irr(irep, 4)%repres
      sym%irr(irep+9, 7)%repres = sym%irr(irep, 7)%repres
      sym%irr(irep+9,19)%repres = sym%irr(irep,19)%repres
    enddo
    !
    sym%irr(5+9,37)%repres = i*m_one
    sym%irr(6+9,37)%repres = i*m_one
    sym%irr(7+9,37)%repres = i*m_one
    sym%irr(8+9,37)%repres = i*m_one
    !
    sym%irr( 9, 1)%repres = g1
    sym%irr( 9, 2)%repres = g2
    sym%irr( 9, 4)%repres = g4
    sym%irr( 9, 7)%repres = g7
    sym%irr( 9,19)%repres = g19
    !
    sym%irr(18, 1)%repres = g1
    sym%irr(18, 2)%repres = g2
    sym%irr(18, 4)%repres = g4
    sym%irr(18, 7)%repres = g7
    sym%irr(18,19)%repres = g19
    !
    sym%irr( 9, 1+36)%repres = g1
    sym%irr(18, 1+36)%repres = g1*m_one
    !
    do irep = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem =1,sym%Nelements(iclass)
          ioper = ioper + 1
          !
          call do_g36_transform(irep,ioper,sym%degen(irep),sym%irr(irep,ioper)%repres)
          !
          f_t = 0
          do k = 1,sym%degen(irep)
              f_t = f_t + (sym%irr(irep,ioper)%repres(k,k))
          enddo
          !
          !sym%characters(irep,iclass) = f_t
          !
        enddo
      enddo
    enddo
    !




    case("G36(EM2)") ! sy 

      sym%Nrepresen = 18
      sym%Noper = 72
      sym%Nclasses = 18
      sym%CII%Noper = 0
                               
      call simple_arrays_allocation
    
      !RE = 2C2+E-   SE = 3C3+E-  ER = 2E+C2-  RR = 4C2+C2- SR = 6C3+C2- ER = 3E+C3- RS = 6C2+C3- SS = 9C3+C3-
      sym%characters = reshape( & 
        ! EE  RE  SE  ER  RR  SR  ER  RS  SS EE' RE' SE' ER' RR' SR' ER' RS' SS' 
        (/ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! A1s
           1,  1,  1,  1,  1,  1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, & ! A2s
           1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1, & ! A3s
           1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1, -1,  1,  1, -1, -1, -1,  1, & ! A4s
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0,  2,  2,  2, -1, -1, -1,  0,  0,  0, & ! E1s
           2,  2, -2, -1, -1,  1,  0,  0,  0,  2,  2, -2, -1, -1,  1,  0,  0,  0, & ! E2s
           2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0, & ! E3s
           2, -1,  0,  2, -1,  0, -2,  1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, & ! E4s
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0,  4, -2,  0, -2,  1,  0,  0,  0,  0, & ! Gs
           !
           1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, & ! A1d
           1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, & ! A2d 
           1,  1, -1,  1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1, & ! A3d
           1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1,  1,  1, -1, & ! A4d
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0, -2, -2, -2,  1,  1,  1,  0,  0,  0, & ! E1d
           2,  2, -2, -1, -1,  1,  0,  0,  0, -2, -2,  2,  1,  1, -1,  0,  0,  0, & ! E2d
           2, -1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0, & ! E3d
           2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0,  2, -1,  0, & ! E4d
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0, -4,  2,  0,  2, -1,  0,  0,  0,  0 /),(/18,18/)) ! Gd
      sym%characters = transpose(sym%characters)
      sym%degen = (/1, 1, 1, 1, 2, 2, 2, 2, 4, 1, 1, 1, 1, 2, 2, 2, 2, 4/)
      sym%Nelements = (/1, 2, 3, 2, 4, 6, 3, 6, 9, 1, 2, 3, 2, 4, 6, 3, 6, 9 /)
      sym%label=(/'A1s', 'A2s', 'A3s', 'A4s', 'E1s', 'E2s', 'E3s', 'E4s', 'Gs ', &
                  'A1d', 'A2d', 'A3d', 'A4d', 'E1d', 'E2d', 'E3d', 'E4d', 'Gd ' /)
      a = 0.5_ark
      b = 0.5_ark*sqrt(3.0_ark)
      e = 1.0_ark
      o = 0.0_ark
      m_one = -1.0_ark       
 
      i = transpose(reshape( (/ e, o, &
                                o, e /), (/ 2, 2/)))
  
      c = transpose(reshape( (/ -a, -b, &
                                 b, -a/), (/ 2, 2/)))
  
    c2 = matmul(c,c)
  
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))

    !sxy = transpose(reshape( (/ o, e, &
    !                            e, o /), (/ 2, 2/)))
    !
    s3 = matmul(c, sxy)
    !
    s2 = matmul(c,s3)
    ! 
    E_rep_1(1,:,:) = i
    E_rep_1(2,:,:) = c
    E_rep_1(3,:,:) = c2
    E_rep_1(4,:,:) = sxy
    E_rep_1(5,:,:) = s2
    E_rep_1(6,:,:) = s3
    
    A2_char(:3) = 1.0_ark
    A2_char(4:6) = -1.0_ark
  
    call irr_allocation 
  
    pos_array = transpose( reshape((/  1,  7,  8, 19, 20, 21, &
                                       2,  9, 11, 22, 24, 26, &
                                       3, 10, 12, 23, 25, 27, &
                                       4, 13, 16, 28, 31, 34, &
                                       5, 14, 17, 29, 32, 35, &
                                       6, 15, 18, 30, 33, 36/), (/6,6/))) 
    ! E1s, E2s, E1d, and E2d 
 
    do j = 1, 6
      do k = 1, 6 
        sym%irr(5, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(5, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)
        
        sym%irr(6, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
        sym%irr(6, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*A2_char(k)
       
        sym%irr(5+9, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(5+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*m_one
        
        sym%irr(6+9, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
        sym%irr(6+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*A2_char(k)*m_one
      end do
    end do

    !E3s and E4s
    
    do j = 1, 6
      do k = 1, 6
        sym%irr(7,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(7,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)

        sym%irr(8,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
        sym%irr(8,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)

        sym%irr(7+9,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(7+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*m_one

        sym%irr(8+9,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
        sym%irr(8+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)*m_one

      end do    
    end do
  
    g1 = transpose(reshape( (/ e, o, o ,o, &
                               o, e, o, o, &
                               o, o, e, o, &
                               o, o, o, e/), (/4,4/)))
    !!
    !g2 = transpose(reshape( (/ -a,-b, o, o, &
    !                            b,-a, o, o, &
    !                            o, o,-a,-b, &
    !                            o, o, b,-a  /), (/4,4/))) 

    !sy01: repalce g2 from !!! to this, did not work
    !
    ! q03 try this 
    g2 = transpose(reshape( (/ -a, b, o, o, &
                               -b,-a, o, o, &
                                o, o,-a, b, &
                                o, o,-b,-a  /), (/4,4/))) 

    !sy3: try this in place of !!  did not work even with one class                          !
    !g4 = transpose(reshape( (/o, o,-e, o, &
    !                           o, o, o, e, &
    !                          -e, o, o, o, &
    !                           o, e, o, o  /), (/4,4/)))

    ! sy2: use this for g4 in place of !! did not help
    ! q01: change again for g4 !!
    ! q03 
    g4 = transpose(reshape( (/ e, o, o, o, &
                               o,-e, o, o, &
                               o, o,-e, o, &
                               o, o, o, e  /), (/4,4/)))
    !!
    !g4 = transpose(reshape( (/-e, o, o, o, &
    !                           o, e, o, o, &
    !                           o, o, e, o, &
    !                           o, o, o,-e  /), (/4,4/)))


    ! sy4: try this instead of !! did not work 
    ! sy17: try this instead of !!  again to make it consistent with sy15, did not help
    ! sy18: use this together with e7-alternative choice, did not work either
    ! q01: try his and also in g19 insetad of !!, there is a sign difference in G2 and G3
    ! q02: try this 
    !g7 = transpose(reshape( (/ -a, o, o,-b, &
    !                            o,-a, b, o, &
    !                            o,-b,-a, o, &
    !                            b, o, o,-a  /), (/4,4/))) 
    !!
    g7 = transpose(reshape( (/ -a, o, o, b, &
                                o,-a,-b, o, &
                                o, b,-a, o, &
                               -b, o, o,-a  /), (/4,4/))) 

    ! sy5: try this instd of !! wrong character 
    !g7 = transpose(reshape( (/-a, b, o, o, &
    !                           -b,-a, o, o, &
    !                            o, o,-a,-b, &
    !                            o, o, b,-a  /), (/4,4/))) 
    ! sy6: try this for !! does not work even for one class
    !g19= transpose(reshape( (/-e, o, o, o, &
    !                           o, e, o, o, &
    !                           o, o,-e, o, &
    !                           o, o, o, e  /), (/4,4/)))

    !!
    g19= transpose(reshape( (/ e, o, o, o, &
                               o, e, o, o, &
                               o, o,-e, o, &
                               o, o, o,-e  /), (/4,4/)))

    ! sy7: try this for !! completely wrong 
    !g19= transpose(reshape( (/o, o, e, o, &
    !                           o, o, o, e, &
    !                           e, o, o, o, &
    !                           o, e, o, o  /), (/4,4/)))


      i = transpose(reshape( (/ e, o, &
                                o, e /), (/ 2, 2/)))
      !!
      c = transpose(reshape( (/ -a, -b, &
                                 b, -a/), (/ 2, 2/)))
      !
      !sy8: try this for !! did not work
      !sy20: try this for !! agian, did not help
      !c = transpose(reshape( (/ -a,  b, &
      !                          -b, -a/), (/ 2, 2/)))

  
    c2 = matmul(c,c)
    !!
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))
    !
    ! sy9 try this for !! did not work
    ! sy21 try this for !! again, made it even worse 
    ! p04 try changing for tau-5 and inverting A-B for (7), did nit work
    ! p05, everythung else is usual, did nit help 
    ! p09 
    !sxy = transpose(reshape( (/ o, e, &
    !                            e, o /), (/ 2, 2/)))
    !
    s3 = matmul(c, sxy)
    !
    s2 = matmul(c,s3)
    !
    sym%irr( 5, 1)%repres = i
    sym%irr( 6, 1)%repres = i
    sym%irr( 7, 1)%repres = i
    sym%irr( 8, 1)%repres = i
    !
    sym%irr( 5, 2)%repres = i
    sym%irr( 6, 2)%repres = i
    sym%irr( 7, 2)%repres = c
    sym%irr( 8, 2)%repres = c
    !!
    sym%irr( 5, 4)%repres = i
    sym%irr( 6, 4)%repres = i*m_one
    sym%irr( 7, 4)%repres = sxy
    sym%irr( 8, 4)%repres = sxy
    !
    ! p09 change as for A1,..A4
    !
    !sym%irr( 5, 4)%repres = i
    !sym%irr( 6, 4)%repres = i*m_one
    !sym%irr( 7, 4)%repres = sxy
    !sym%irr( 8, 4)%repres = sxy*m_one
    !
    !!
    sym%irr( 5, 7)%repres = c2
    sym%irr( 6, 7)%repres = c2
    sym%irr( 7, 7)%repres = i
    sym%irr( 8, 7)%repres = i
    !
    ! sy18: use this together with g7-alternative choice, did no help
    ! sy19: only this option without g7-swap, did not help
    !
    !sym%irr( 5, 7)%repres = c
    !sym%irr( 6, 7)%repres = c
    !sym%irr( 7, 7)%repres = i
    !sym%irr( 8, 7)%repres = i
    !
    ! p06 try changing sign for oper 6  did not work 
    !
    !sym%irr( 5, 7)%repres = c2
    !sym%irr( 6, 7)%repres = c2*m_one
    !sym%irr( 7, 7)%repres = i
    !sym%irr( 8, 7)%repres = i
    !!
    sym%irr( 5,19)%repres = sxy
    sym%irr( 6,19)%repres = sxy
    sym%irr( 7,19)%repres = i
    sym%irr( 8,19)%repres = i*m_one
    !
    ! p07 change E2 sign according with character of A2, did not help
    ! p08 the same but with sxy from alternative choice , did not work
    !sym%irr( 5,19)%repres = sxy
    !sym%irr( 6,19)%repres = sxy*m_one
    !sym%irr( 7,19)%repres = i
    !sym%irr( 8,19)%repres = i*m_one
    ! p10 another test, did not work 
    !sym%irr( 5,19)%repres = sxy*m_one
    !sym%irr( 6,19)%repres = sxy
    !sym%irr( 7,19)%repres = i
    !sym%irr( 8,19)%repres = i*m_one
    !
    sym%irr( 5,37)%repres = i
    sym%irr( 6,37)%repres = i
    sym%irr( 7,37)%repres = i
    sym%irr( 8,37)%repres = i
    !
    do irep = 5,8
      sym%irr(irep+9, 1)%repres = sym%irr(irep, 1)%repres
      sym%irr(irep+9, 2)%repres = sym%irr(irep, 2)%repres
      sym%irr(irep+9, 4)%repres = sym%irr(irep, 4)%repres
      sym%irr(irep+9, 7)%repres = sym%irr(irep, 7)%repres
      sym%irr(irep+9,19)%repres = sym%irr(irep,19)%repres
    enddo
    !
    sym%irr(5+9,37)%repres = i*m_one
    sym%irr(6+9,37)%repres = i*m_one
    sym%irr(7+9,37)%repres = i*m_one
    sym%irr(8+9,37)%repres = i*m_one
    !
    sym%irr( 9, 1)%repres = g1
    sym%irr( 9, 2)%repres = g2
    sym%irr( 9, 4)%repres = g4
    sym%irr( 9, 7)%repres = g7
    sym%irr( 9,19)%repres = g19
    !
    sym%irr(18, 1)%repres = g1
    sym%irr(18, 2)%repres = g2
    sym%irr(18, 4)%repres = g4
    sym%irr(18, 7)%repres = g7
    sym%irr(18,19)%repres = g19
    !
    sym%irr( 9, 1+36)%repres = g1
    sym%irr(18, 1+36)%repres = g1*m_one
    !
    do irep = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem =1,sym%Nelements(iclass)
          ioper = ioper + 1
          !
          call do_g36_transform(irep,ioper,sym%degen(irep),sym%irr(irep,ioper)%repres)
          !
          f_t = 0
          do k = 1,sym%degen(irep)
              f_t = f_t + (sym%irr(irep,ioper)%repres(k,k))
          enddo
          !
          !sym%characters(irep,iclass) = f_t
          !
        enddo
      enddo
    enddo





    case("G36(EM0)") !! Tom's

      sym%Nrepresen = 18
      sym%Noper = 72
      sym%Nclasses = 18
      sym%CII%Noper = 0
                               
      call simple_arrays_allocation
    
      !RE = 2C2+E-   SE = 3C3+E-  ER = 2E+C2-  RR = 4C2+C2- SR = 6C3+C2- ER = 3E+C3- RS = 6C2+C3- SS = 9C3+C3-
      sym%characters = reshape( & 
        ! EE  RE  SE  ER  RR  SR  ER  RS  SS EE' RE' SE' ER' RR' SR' ER' RS' SS' 
        (/ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! A1s
           1,  1,  1,  1,  1,  1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, & ! A2s
           1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1, & ! A3s
           1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1, -1,  1,  1, -1, -1, -1,  1, & ! A4s
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0,  2,  2,  2, -1, -1, -1,  0,  0,  0, & ! E1s
           2,  2, -2, -1, -1,  1,  0,  0,  0,  2,  2, -2, -1, -1,  1,  0,  0,  0, & ! E2s
           2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0, & ! E3s
           2, -1,  0,  2, -1,  0, -2,  1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, & ! E4s
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0,  4, -2,  0, -2,  1,  0,  0,  0,  0, & ! Gs
           !
           1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, & ! A1d
           1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, & ! A2d 
           1,  1, -1,  1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1, & ! A3d
           1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1,  1,  1, -1, & ! A4d
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0, -2, -2, -2,  1,  1,  1,  0,  0,  0, & ! E1d
           2,  2, -2, -1, -1,  1,  0,  0,  0, -2, -2,  2,  1,  1, -1,  0,  0,  0, & ! E2d
           2, -1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0, & ! E3d
           2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0,  2, -1,  0, & ! E4d
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0, -4,  2,  0,  2, -1,  0,  0,  0,  0 /),(/18,18/)) ! Gd
      sym%characters = transpose(sym%characters)
      sym%degen = (/1, 1, 1, 1, 2, 2, 2, 2, 4, 1, 1, 1, 1, 2, 2, 2, 2, 4/)
      sym%Nelements = (/1, 2, 3, 2, 4, 6, 3, 6, 9, 1, 2, 3, 2, 4, 6, 3, 6, 9 /)
      sym%label=(/'A1s', 'A2s', 'A3s', 'A4s', 'E1s', 'E2s', 'E3s', 'E4s', 'Gs ', &
                  'A1d', 'A2d', 'A3d', 'A4d', 'E1d', 'E2d', 'E3d', 'E4d', 'Gd ' /)
      a = 0.5_ark
      b = 0.5_ark*sqrt(3.0_ark)
      e = 1.0_ark
      o = 0.0_ark
      m_one = -1.0_ark       
 
      !E rep of C_3v 

      i = transpose(reshape( (/ e, o, &
                                o, e /), (/ 2, 2/)))
  
      c = transpose(reshape( (/ -a, -b, &
                                 b, -a/), (/ 2, 2/)))
  
    c2 = matmul(c,c)
    !!
    sxy = transpose(reshape( (/ -e,  o, &
                                 o,  e /), (/ 2, 2/)))
    !
    ! sy22 try this for !! made it worse
    !sxy = transpose(reshape( (/ o, e, &
    !                            e, o /), (/ 2, 2/)))
    !
    s3 = matmul(c, sxy)
    !
    s2 = matmul(c,s3)
    ! 
    E_rep_1(1,:,:) = i
    E_rep_1(2,:,:) = c
    E_rep_1(3,:,:) = c2
    E_rep_1(4,:,:) = sxy
    E_rep_1(5,:,:) = s2
    E_rep_1(6,:,:) = s3
    
    A2_char(:3) = 1.0_ark
    A2_char(4:6) = -1.0_ark
  
    call irr_allocation 
  
    pos_array = transpose( reshape((/  1,  7,  8, 19, 20, 21, &
                                       2,  9, 11, 22, 24, 26, &
                                       3, 10, 12, 23, 25, 27, &
                                       4, 13, 16, 28, 31, 34, &
                                       5, 14, 17, 29, 32, 35, &
                                       6, 15, 18, 30, 33, 36/), (/6,6/))) 
    ! E1s, E2s, E1d, and E2d 
 
    do j = 1, 6
      do k = 1, 6 
        sym%irr(5, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(5, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)
        
        sym%irr(6, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
        sym%irr(6, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*A2_char(k)
       
        sym%irr(5+9, pos_array(k,j))%repres = E_rep_1(j,:,:)
        sym%irr(5+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*m_one
        
        sym%irr(6+9, pos_array(k,j))%repres = E_rep_1(j,:,:)*A2_char(k)
        sym%irr(6+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*A2_char(k)*m_one
      end do
    end do

    !E3s and E4s
    
    do j = 1, 6
      do k = 1, 6
        sym%irr(7,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(7,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)

        sym%irr(8,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
        sym%irr(8,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)

        sym%irr(7+9,pos_array(k,j))%repres = E_rep_1(k,:,:)
        sym%irr(7+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*m_one

        sym%irr(8+9,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
        sym%irr(8+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)*m_one

      end do    
    end do
    
    ! Gs and Gd

    ! Implements outer product,
    ! consider two matrices
    !   a b    and e f
    !   c d        g h
    ! then the outer product is 
    !   ae af be bf 
    !   ag ah bg bh
    !   ce cf de df
    !   cg ch dg dh

    do j = 1, 6
      do k = 1, 6
        do r = 1, 2
          do s = 1, 2
            sym%irr(9,pos_array(k,j))%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)
            sym%irr(9,pos_array(k,j)+36)%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)

            sym%irr(9+9,pos_array(k,j))%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)
            sym%irr(9+9,pos_array(k,j)+36)%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)*m_one
          end do
        end do
      end do 
    end do

!  do j = 5, 8
!    do p= 1,36
!       write(*,*) "j = ", j, ", p= ", p  
!       do r = 1,2
!            write(*, "(2(f12.6,1x))") sym%irr(j,p)%repres(r,1),sym%irr(j,p)%repres(r,2)
!        end do
!       write(*,*) char(10)
!    end do
!  end do

!  do p=1,36
!    write(*,*) "j=", j, " p= ", p
!    do r=1,4
!        write(*, "(4(f12.6,1x))") sym%irr(9,p)%repres(r,1), sym%irr(9,p)%repres(r,2), sym%irr(9,p)%repres(r,3),&
!                                   sym%irr(9,p)%repres(r,4)
!    end do 
!    write(*,*) char(10)
!  end do 


       ! characters as traces of the corresponding representations 
       
       !do irep = 1,sym%Nrepresen
       !  ioper = 0
       !  do iclass = 1,sym%Nclasses
       !    do ielem =1,sym%Nelements(iclass)
       !      ioper = ioper + 1
       !      f_t = 0
       !      do k = 1,sym%degen(irep)
       !          f_t = f_t + (sym%irr(irep,ioper)%repres(k,k))
       !      enddo
       !      sym%characters(irep,iclass) = f_t
       !    enddo
       !  enddo
       !enddo
      !
 case("G36(EM)") !! Tom's
      !
      sym%Nrepresen = 18
      sym%Noper = 72
      sym%Nclasses = 18
      sym%CII%Noper = 0
      !                               
      call simple_arrays_allocation
      !
      tempG36(1:36)   = (/0,0,2,0,6,4,0,7,2,3,2,3,4,5,6,4,5,6,0,21,19, 2, 3, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6/)
      tempG36(73:108) = (/0,0,2,0,2,2,0,7,7,7,8,8,7,7,7,8,8,8,0, 7, 7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
      tempG36(37:72)  = (/0,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,&
                          37,37,37,37,37,37/)
      tempG36(109:144)= (/0, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,&
                          31,32,33,34,35,36/)
      !
      allocate(sym%product_table(72,2),stat=alloc)
      call ArrayStart('sym%product_table',alloc,size(sym%product_table),kind(sym%product_table))
      !
      sym%product_table = reshape( tempG36, (/ 72, 2/))
      sym%product_table_set = .true.
      !    
      !RE = 2C2+E-   SE = 3C3+E-  ER = 2E+C2-  RR = 4C2+C2- SR = 6C3+C2- ER = 3E+C3- RS = 6C2+C3- SS = 9C3+C3-
      sym%characters = reshape( & 
        ! EE  RE  SE  ER  RR  SR  ER  RS  SS EE' RE' SE' ER' RR' SR' ER' RS' SS' 
        (/ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, & ! A1s
           1,  1,  1,  1,  1,  1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, & ! A2s
           1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1,  1,  1, -1, & ! A3s
           1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1, -1,  1,  1, -1, -1, -1,  1, & ! A4s
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0,  2,  2,  2, -1, -1, -1,  0,  0,  0, & ! E1s
           2,  2, -2, -1, -1,  1,  0,  0,  0,  2,  2, -2, -1, -1,  1,  0,  0,  0, & ! E2s
           2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0,  2, -1,  0, & ! E3s
           2, -1,  0,  2, -1,  0, -2,  1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, & ! E4s
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0,  4, -2,  0, -2,  1,  0,  0,  0,  0, & ! Gs
           !
           1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, & ! A1d
           1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, & ! A2d 
           1,  1, -1,  1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1, & ! A3d
           1,  1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1, -1, -1,  1,  1,  1, -1, & ! A4d
           !
           2,  2,  2, -1, -1, -1,  0,  0,  0, -2, -2, -2,  1,  1,  1,  0,  0,  0, & ! E1d
           2,  2, -2, -1, -1,  1,  0,  0,  0, -2, -2,  2,  1,  1, -1,  0,  0,  0, & ! E2d
           2, -1,  0,  2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0, & ! E3d
           2, -1,  0,  2, -1,  0, -2,  1,  0, -2,  1,  0, -2,  1,  0,  2, -1,  0, & ! E4d
           !
           4, -2,  0, -2,  1,  0,  0,  0,  0, -4,  2,  0,  2, -1,  0,  0,  0,  0 /),(/18,18/)) ! Gd 
           !
      sym%characters = transpose(sym%characters)
      sym%degen = (/1, 1, 1, 1, 2, 2, 2, 2, 4, 1, 1, 1, 1, 2, 2, 2, 2, 4/)
      sym%Nelements = (/1, 2, 3, 2, 4, 6, 3, 6, 9, 1, 2, 3, 2, 4, 6, 3, 6, 9 /)
      sym%label=(/'A1s', 'A2s', 'A3s', 'A4s', 'E1s', 'E2s', 'E3s', 'E4s', 'Gs ', &
                  'A1d', 'A2d', 'A3d', 'A4d', 'E1d', 'E2d', 'E3d', 'E4d', 'Gd ' /)
      !
      a = 0.5_ark
      b = 0.5_ark*sqrt(3.0_ark)
      e = 1.0_ark
      o = 0.0_ark
      m_one = -1.0_ark       
      ! 
      !E rep of C_3v 
      i = transpose(reshape( (/ e, o, &
                                o, e /), (/ 2, 2/)))
  
      c = transpose(reshape( (/ -a, -b, &
                                 b, -a/), (/ 2, 2/)))
      !
      c2 = matmul(c,c)
      !
      sxy = transpose(reshape( (/ e,  o, &
                                  o, -e /), (/ 2, 2/)))
      !
      s3 = matmul(c, sxy)
      !
      s2 = matmul(c,s3)
      ! 
      E_rep_1(1,:,:) = i
      E_rep_1(2,:,:) = c
      E_rep_1(3,:,:) = c2
      E_rep_1(4,:,:) = sxy
      E_rep_1(5,:,:) = s2
      E_rep_1(6,:,:) = s3
      
      E_rep_2(1,:,:) = i
      E_rep_2(2,:,:) = c
      E_rep_2(3,:,:) = c2
      E_rep_2(4,:,:) = -sxy
      E_rep_2(5,:,:) = -s2
     E_rep_2(6,:,:) = -s3
    
      A2_char(:3) = 1.0_ark
     A2_char(4:6) = -1.0_ark
    
     call irr_allocation 
    
      pos_array = transpose( reshape((/  1,  7,  8, 19, 20, 21, &
                                         2,  9, 11, 22, 24, 26, &
                                         3, 10, 12, 23, 25, 27, &
                                         4, 13, 16, 28, 31, 34, &
                                         5, 14, 17, 29, 32, 35, &
                                         6, 15, 18, 30, 33, 36/), (/6,6/))) 
     ! E1s, E2s, E1d, and E2d 
    
      do j = 1, 6
        do k = 1, 6 
          sym%irr(5, pos_array(k,j))%repres = E_rep_1(j,:,:)
          sym%irr(5, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)
          
          sym%irr(6, pos_array(k,j))%repres = E_rep_2(j,:,:)*A2_char(k)
          sym%irr(6, pos_array(k,j)+36)%repres = E_rep_2(j,:,:)*A2_char(k)
         
          sym%irr(5+9, pos_array(k,j))%repres = E_rep_1(j,:,:)
          sym%irr(5+9, pos_array(k,j)+36)%repres = E_rep_1(j,:,:)*m_one
          
          sym%irr(6+9, pos_array(k,j))%repres = E_rep_2(j,:,:)*A2_char(k)
          sym%irr(6+9, pos_array(k,j)+36)%repres = E_rep_2(j,:,:)*A2_char(k)*m_one
        end do
     end do
    
      !E3s and E4s
      
      do j = 1, 6
        do k = 1, 6
          sym%irr(7,pos_array(k,j))%repres = E_rep_1(k,:,:)
         sym%irr(7,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)
    
          sym%irr(8,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
         sym%irr(8,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)
    
          sym%irr(7+9,pos_array(k,j))%repres = E_rep_1(k,:,:)
         sym%irr(7+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*m_one
    
          sym%irr(8+9,pos_array(k,j))%repres = E_rep_1(k,:,:)*A2_char(j)
         sym%irr(8+9,pos_array(k,j)+36)%repres = E_rep_1(k,:,:)*A2_char(j)*m_one
    
        end do    
      end do
      
     ! Gs and Gd
    
      ! Implements outer product,
      ! consider two matrices
      !   a b    and e f
      !   c d        g h
      ! then the outer product is 
      !   ae af be bf 
      !   ag ah bg bh
      !   ce cf de df
     !   cg ch dg dh
    
      do j = 1, 6
        do k = 1, 6
          do r = 1, 2
            do s = 1, 2
              sym%irr(9,pos_array(k,j))%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)
             sym%irr(9,pos_array(k,j)+36)%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)
    
              sym%irr(9+9,pos_array(k,j))%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)
              sym%irr(9+9,pos_array(k,j)+36)%repres(2*r-1:2*r,2*s-1:2*s) = E_rep_1(k, r, s)*E_rep_1(j, :, :)*m_one
            end do
          end do
        end do 
     end do




  case("D2H(M)")
    !
    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2a  C2b C2c E* sab sac sbc i
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ag 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Au 
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! B1g
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! B1u
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! B2g 
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! B2u
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! B3g
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! B3u
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ag ','Au ','B1g','B1u','B2g','B2u','B3g','B3u'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/o,pi,o/) ! Ry
    sym%euler( 4,:) = (/p2,pi,p3/) ! Rx
    !
    sym%euler( 5,:) = (/p2,pi,p3/)
    sym%euler( 6,:) = (/o,pi,o/)
    sym%euler( 7,:) = (/pi,o,o/)
    sym%euler( 8,:) = 0
    !
    call irr_allocation
    !
  case("D2H")

    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2z  C2y C2x i  sxy sxz syz  
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ag 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Au 
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! B1g
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! B1u
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! B2g 
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! B2u
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! B3g
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! B3u
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ag ','Au ','B1g','B1u','B2g','B2u','B3g','B3u'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/o,pi,o/) ! Ry
    sym%euler( 4,:) = (/p2,pi,p3/) ! Rx
    !
    sym%euler( 7,:) = (/p2,pi,p3/)
    sym%euler( 8,:) = (/o,pi,o/)
    sym%euler( 5,:) = (/pi,o,o/)
    sym%euler( 6,:) = 0
    !
    call irr_allocation    
    !
  case("C3V(M)","C3V")
    !
    sym%Nrepresen=3
    sym%Noper=6
    sym%Nclasses=3
    sym%CII%Noper = 0
    !
    call simple_arrays_allocation


    sym%characters= reshape( &      !A1 A2 E
                                  (/ 1, 1, 2, &  
                                     1, 1,-1, &  
                                     1,-1, 0 /),(/3,3/)) 
    sym%degen=(/1,1,2/)
    sym%Nelements=(/1,2,3/)
    sym%label=(/'A1','A2','E '/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 2.0_ark/3.0_ark*pi
    p4 = 4.0_ark/3.0_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p4,o ,o /)
    sym%euler( 3,:) = (/p3,o ,o /)
    sym%euler( 4,:) = (/o ,pi,o /)
    sym%euler( 5,:) = (/p3,pi,o /)
    sym%euler( 6,:) = (/p4,pi,o /)
    !
    sym%lquant(3) = 1
    !
    call irr_allocation
       !
       sym%irr(3,1)%repres = reshape((/1.0_ark,               0.0_ark, &
                                       0.0_ark,               1.0_ark/),(/2,2/))
       !
       sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                        0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
       !
       sym%irr(3,3)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
                                       -0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
       !
       sym%irr(3,4)%repres = reshape((/ 1.0_ark,               0.0_ark, &
                                        0.0_ark,              -1.0_ark/),(/2,2/))
       !
       sym%irr(3,5)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
                                        0.5_ark*sqrt(3.0_ark), 0.5_ark            /),(/2,2/))
       !
       sym%irr(3,6)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                       -0.5_ark*sqrt(3.0_ark), 0.5_ark             /),(/2,2/))
       !
  case("D3H(M)","D3H")

    sym%Nrepresen=6
    sym%Noper=12
    sym%Nclasses=6
    sym%CII%Noper = 0

    call simple_arrays_allocation
    !
    sym%characters= reshape( &      !A1 A2 E  A1 A2 E
                                  (/ 1, 1, 2, 1, 1, 2, &  
                                     1, 1,-1, 1, 1,-1, &  
                                     1,-1, 0, 1,-1, 0, &          
                                     1, 1, 2,-1,-1,-2, &  
                                     1, 1,-1,-1,-1, 1, &  
                                     1,-1, 0,-1, 1, 0 /),(/6,6/)) 
    sym%degen=(/1,1,2,1,1,2/)
    sym%Nelements=(/1,2,3,1,2,3/)
    sym%label=(/'A1''','A2''','E'' ','A1"','A2"','E" '/)
    !
    o   = 0.0_ark
    p2  = 0.5_ark*pi
    p3  = pi/3.0_ark
    p23 = 2.0_ark/3.0_ark*pi
    p43 = 4.0_ark/3.0_ark*pi
    !
    sym%euler( 1,:) = 0               ! (E)
    sym%euler( 2,:) = (/p23,o ,o /)   ! (132)
    sym%euler( 3,:) = (/p43,o ,o /)   ! (123)
    sym%euler( 4,:) = (/ pi,pi,o /)   ! (23)
    sym%euler( 5,:) = (/ p3,pi,o /)   ! (12)
    sym%euler( 6,:) = (/-p3,pi,o /)   ! (13)
    sym%euler( 7,:) = (/ pi, o,o /)   ! (E)*
    sym%euler( 8,:) = (/-p3, o,o /)   ! (132)*
    sym%euler( 9,:) = (/ p3, o,o /)   ! (123)*
    sym%euler(10,:) = (/  o,pi,o /)   ! (23)*
    sym%euler(11,:) = (/p43,pi,o /)   ! (12)*
    sym%euler(12,:) = (/p23,pi,o /)   ! (13)*
    !
    call irr_allocation
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    sym%irr(3,1)%repres = reshape((/ e, o,  &
                                     o, e/),(/2,2/))
    !
    sym%irr(3,2)%repres = reshape((/-a, b,  &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(3,3)%repres = reshape((/-a,-b,  &
                                     b,-a/),(/2,2/))
    !
    sym%irr(3,4)%repres = reshape((/ e, o,  &
                                     o,-e/),(/2,2/))
    !
    sym%irr(3,5)%repres = reshape((/-a,-b,  &
                                    -b, a /),(/2,2/))
    !
    sym%irr(3,6)%repres = reshape((/-a, b,  &
                                     b, a  /),(/2,2/))
    !
    !
    sym%irr(3,7)%repres = reshape((/ e, o, &
                                     o, e/),(/2,2/))
    !
    sym%irr(3,8)%repres = reshape((/-a, b, &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(3,9)%repres = reshape((/-a,-b, &
                                     b,-a/),(/2,2/))
    !
    sym%irr(3,10)%repres= reshape((/ e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(3,11)%repres= reshape((/-a,-b, &
                                    -b, a /),(/2,2/))
    !
    sym%irr(3,12)%repres= reshape((/-a, b, &
                                     b, a  /),(/2,2/))
    !
    sym%irr(6,1)%repres = reshape((/e, o,  &
                                    o, e/),(/2,2/))
    !
    sym%irr(6,2)%repres = reshape((/-a, b, &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(6,3)%repres = reshape((/-a,-b, &
                                     b,-a/),(/2,2/))
    !
    sym%irr(6,4)%repres = reshape((/-e, o, &
                                     o, e/),(/2,2/))
    !
    sym%irr(6,5)%repres = reshape((/ a, b, &
                                     b,-a /),(/2,2/))
    !
    sym%irr(6,6)%repres = reshape((/ a,-b, &
                                    -b,-a  /),(/2,2/))
    !
    sym%irr(6,7)%repres = reshape((/-e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(6,8)%repres = reshape((/ a,-b, &
                                     b, a/),(/2,2/))
    !
    sym%irr(6,9)%repres = reshape((/ a, b, &
                                    -b, a/),(/2,2/))
    !
    sym%irr(6,10)%repres= reshape((/ e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(6,11)%repres= reshape((/-a,-b, &
                                    -b, a /),(/2,2/))
    !
    sym%irr(6,12)%repres= reshape((/-a, b, &
                                     b, a  /),(/2,2/))
    !
    sym%lquant(1:2) = 0
    sym%lquant(4:5) = 0
    sym%lquant(3) = 1
    sym%lquant(6) = 1
    !
  case("D3D(M)","D3D")
    !
    sym%Nrepresen=6
    sym%Noper=12
    sym%Nclasses=6
    sym%CII%Noper = 0
    !
    call simple_arrays_allocation
    !
    sym%characters= reshape( &      !E  2C3  3C2'  i  2S6  3sd
                                  (/ 1,  1,   1,   1,  1,   1,  &  !A1g
                                     1,  1,  -1,   1,  1,  -1,  &  !A2g
                                     2, -1,   0,   2, -1,   0,  &  !Eg         
                                     1,  1,   1,  -1, -1,  -1,  &  !A1u
                                     1,  1,  -1,  -1, -1,   1,  &  !A2u
                                     2, -1,   0,  -2,  1,   0 /),(/6,6/))   !Eu
    !
    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,2,1,1,2/)
    sym%Nelements=(/1,2,3,1,2,3/)
    sym%label=(/'A1g''','A2g''','Eg'' ','A1u"','A2u"','Eu" '/)
    !
    o   = 0.0_ark
    p2  = 0.5_ark*pi
    p3  = pi/3.0_ark
    p23 = 2.0_ark/3.0_ark*pi
    p43 = 4.0_ark/3.0_ark*pi
    !
    sym%euler( 1,:) = 0               ! (E)
    sym%euler( 2,:) = (/p23,o ,o /)   ! (132)
    sym%euler( 3,:) = (/p43,o ,o /)   ! (123)
    sym%euler( 4,:) = (/ pi,pi,o /)   ! (23)
    sym%euler( 5,:) = (/ p3,pi,o /)   ! (12)
    sym%euler( 6,:) = (/-p3,pi,o /)   ! (13)
    sym%euler( 7,:) = (/ pi, o,o /)   ! (E)*
    sym%euler( 8,:) = (/-p3, o,o /)   ! (132)*
    sym%euler( 9,:) = (/ p3, o,o /)   ! (123)*
    sym%euler(10,:) = (/  o,pi,o /)   ! (23)*
    sym%euler(11,:) = (/p43,pi,o /)   ! (12)*
    sym%euler(12,:) = (/p23,pi,o /)   ! (13)*
    !
    call irr_allocation
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    sym%irr(3,1)%repres = reshape((/e, o,  &
                                    o, e/),(/2,2/))
    !
    sym%irr(3,2)%repres = reshape((/-a,  b, &
                                    -b, -a/),(/2,2/))
    !
    sym%irr(3,3)%repres = reshape((/-a, -b, &
                                     b, -a/),(/2,2/))
    !
    sym%irr(3,4)%repres = reshape((/ -a,  b, &
                                      b,  a /),(/2,2/))
    !
    sym%irr(3,5)%repres = reshape((/e, o, &
                                    o,-e/),(/2,2/))
    !
    sym%irr(3,6)%repres = reshape((/ -a,-b, &
                                     -b, a  /),(/2,2/))
    !
    sym%irr(3,7)%repres = reshape((/ e, o, &
                                     o, e/),(/2,2/))
    !
    sym%irr(3,8)%repres = reshape((/ -a, -b, &
                                      b, -a/),(/2,2/))
    !
    sym%irr(3,9)%repres = reshape((/ -a,  b, &
                                     -b, -a/),(/2,2/))
    !
    sym%irr(3,10)%repres= reshape((/ -a, -b, &
                                     -b,  a /),(/2,2/))
    !
    sym%irr(3,11)%repres= reshape((/  e, o, &
                                      o,-e/),(/2,2/))
    !
    sym%irr(3,12)%repres= reshape((/ -a,  b, &
                                      b,  a  /),(/2,2/))
    !
    sym%irr(6,1)%repres = reshape((/e, o,  &
                                    o, e/),(/2,2/))
    !
    sym%irr(6,2)%repres = reshape((/-a,  b, &
                                    -b, -a/),(/2,2/))
    !
    sym%irr(6,3)%repres = reshape((/-a, -b, &
                                     b, -a/),(/2,2/))
    !
    sym%irr(6,5)%repres = reshape((/  e, o, &
                                      o,-e/),(/2,2/))
    !
    sym%irr(6,4)%repres = reshape((/ -a,  b, &
                                      b,  a /),(/2,2/))
    !
    sym%irr(6,6)%repres = reshape((/ -a, -b, &
                                     -b,  a  /),(/2,2/))
    !
    sym%irr(6,7)%repres = reshape((/-e, o, &
                                     o, -e/),(/2,2/))
    !
    sym%irr(6,8)%repres = reshape((/ a,  b, &
                                    -b,  a/),(/2,2/))
    !
    sym%irr(6,9)%repres = reshape((/ a, -b, &
                                     b,  a/),(/2,2/))
    !
    sym%irr(6,11)%repres= reshape((/ -e,  o, &
                                      o,  e/),(/2,2/))
    !
    sym%irr(6,10)%repres= reshape((/  a,  b, &
                                      b, -a /),(/2,2/))
    !
    sym%irr(6,12)%repres= reshape((/  a,  -b, &
                                     -b,  -a  /),(/2,2/))
    !
    sym%lquant(1:2) = 0
    sym%lquant(4:5) = 0
    sym%lquant(3) = 1
    sym%lquant(6) = 1

    allocate (characters_(sym%Nrepresen,sym%Nclasses),stat=alloc)
    !
    if (alloc/=0) stop 'characters_ - out of memory'
    !
    characters_ = sym%characters
    !
    ! characters as traces of the corresponding representations 
    !
    do irep = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem =1,sym%Nelements(iclass)
          ioper = ioper + 1
          f_t = 0
          do k = 1,sym%degen(irep)
              f_t = f_t + (sym%irr(irep,ioper)%repres(k,k))
          enddo
          sym%characters(irep,iclass) = f_t
          !
          if (abs(sym%characters(irep,iclass)-characters_(irep,iclass))>small_) then 
            !
            write(out,"(' characters do not agree: irep,iclass = ',4i8)") irep,iclass,sym%characters(irep,iclass), &
                  characters_(irep,iclass)
            stop 'simple_arrays_allocation - out of memory'
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    deallocate(characters_)
    !
  case("TD(M)","TD")

    sym%Nrepresen=5
    sym%Noper=24
    sym%Nclasses=5
    sym%CII%Noper = 6

    call simple_arrays_allocation

    sym%CII%coeff = (/1.0_ark,1.0_ark,4.0_ark,10.0_ark,1.0_ark,1.0_ark/)
    sym%CII%ioper = (/19,20,21,22,23,24/)



    sym%characters= reshape( &      !A1 A2 E  F1 F2
                                      (/ 1, 1, 2, 3, 3, &  
                                         1, 1,-1, 0, 0, &  
                                         1, 1, 2,-1,-1, &
                                         1,-1, 0, 1,-1, &
                                         1,-1, 0,-1, 1  /),(/5,5/)) 
    sym%degen=(/1,1,2,3,3/)
    sym%Nelements=(/1,8,3,6,6/)
    sym%label=(/'A1','A2','E ','F1','F2'/)

    sym%CII%Nzeta = sum(sym%degen(:))
    allocate(sym%CII%izeta(sym%CII%Nzeta),stat=alloc)
    !
    sym%CII%izeta = (/18,-18,12,-12,-14,-8,4,14,-4,8/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p3,p2,o /)
    sym%euler( 3,:) = (/pi,p2,p3/)
    sym%euler( 4,:) = (/ o,p2,p3/)
    sym%euler( 5,:) = (/p3,p2,pi/)
    sym%euler( 6,:) = (/p2,p2,o /)
    sym%euler( 7,:) = (/pi,p2,p2/)
    sym%euler( 8,:) = (/ o,p2,p2/)
    sym%euler( 9,:) = (/p2,p2,pi/)
    sym%euler(10,:) = (/p2,pi,p3/)
    sym%euler(11,:) = (/ o,pi,o /)
    sym%euler(12,:) = (/pi, o,o /)
    sym%euler(13,:) = (/ o,p3,o /)
    sym%euler(14,:) = (/ o,p2,o /)
    sym%euler(15,:) = (/p2, o,o /)
    sym%euler(16,:) = (/p3, o,o /)
    sym%euler(17,:) = (/p2,p3,p3/)
    sym%euler(18,:) = (/p2,p2,p3/)
    sym%euler(19,:) = (/p2,p2,p2/)
    sym%euler(20,:) = (/ o,p2,pi/)
    sym%euler(21,:) = (/ o,pi,p2/)
    sym%euler(22,:) = (/pi,pi,p2/)
    sym%euler(23,:) = (/pi,p2,o /)
    sym%euler(24,:) = (/p3,p2,p3/)
    !
    !sym%euler(:,1) = sym%euler(:,1)+pi*0.25_ark
    !sym%euler(:,2) = sym%euler(:,2)
    !sym%euler(:,3) = sym%euler(:,3)-pi*0.25_ark
    !
    call irr_allocation
    !
    sym%irr(3,1)%repres = reshape((/ 1.0_ark,               0.0_ark, &
                                     0.0_ark,               1.0_ark/),(/2,2/))
    !
    !sym%irr(3,2)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
    !                                -0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
    ! working 
    sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                     0.5_ark*sqrt(3.0_ark),-0.5_ark             /),(/2,2/))
    !
    !sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
    !                                 0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
    !
    ! 
    !sym%irr(3,21)%repres = reshape((/-0.5_ark,              0.5_ark*sqrt(3.0_ark),&
    !                                  0.5_ark*sqrt(3.0_ark),0.5_ark             /),(/2,2/))
    !working    
    sym%irr(3,21)%repres = reshape((/-0.5_ark,              0.5_ark*sqrt(3.0_ark),&
                                      0.5_ark*sqrt(3.0_ark),0.5_ark             /),(/2,2/))

    !
    sym%irr(4,1)%repres = reshape((/ 1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark, 1.0_ark, 0.0_ark, &
                                     0.0_ark, 0.0_ark, 1.0_ark/),(/3,3/))
    !
    !sym%irr(4,2)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    !!sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !!                                 1.0_ark, 0.0_ark, 0.0_ark, &
    !!                                 0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
    ! working ->
    !
    sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
                                     1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))

    !sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark, 1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark, &
    !                                 0.0_ark, 1.0_ark, 0.0_ark/),(/3,3/))


    !sym%irr(4,2)%repres = transpose(reshape((/ 0.0_ark, 0.5_ark*sqrt(2.0_ark), 0.5_ark*sqrt(2.0_ark), &
    !                                 0.5_ark*sqrt(2.0_ark), 0.5_ark,-0.5_ark, &
    !                                -0.5_ark*sqrt(2.0_ark), 0.5_ark,-0.5_ark/),(/3,3/)))
    !
    !
    !sym%irr(4,21)%repres = reshape((/-1.0_ark, 0.0_ark, 0.0_ark, &
    !                                  0.0_ark, 0.0_ark,-1.0_ark, &
    !                                  0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))

    ! working->
    sym%irr(4,21)%repres = reshape((/  0.0_ark, 1.0_ark, 0.0_ark, &
                                       1.0_ark, 0.0_ark, 0.0_ark, &
                                       0.0_ark, 0.0_ark,-1.0_ark/),(/3,3/))


    !
    ! roman
    !

     !  sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
     !                                   1.0_ark, 0.0_ark, 0.0_ark, &
     !                                   0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
     !  !
     !  sym%irr(4,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
     !                                    0.0_ark,-1.0_ark, 0.0_ark, &
     !                                   -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))




    !
    !sym%irr(4,21)%repres = reshape((/-1.0_ark, 0.0_ark, 0.0_ark, &
    !                                  0.0_ark, 1.0_ark, 0.0_ark, &
    !                                  0.0_ark, 0.0_ark,-1.0_ark/),(/3,3/))
    
    !sym%irr(4,21)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 -1.0_ark, 0.0_ark, 0.0_ark, &
    !!                                  0.0_ark, 0.0_ark, -1.0_ark/),(/3,3/))
    !


    !
    !!!sym%irr(4,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !!!                                  0.0_ark,-1.0_ark, 0.0_ark, &
    !!!                                 -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,1)%repres = reshape((/ 1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark, 1.0_ark, 0.0_ark, &
                                     0.0_ark, 0.0_ark, 1.0_ark/),(/3,3/))
    ! working
    !sym%irr(5,2)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,2)%repres = reshape((/ 0.0_ark, 0.0_ark, 1.0_ark, &
                                    -1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
    ! working
    !sym%irr(5,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                  0.0_ark, 1.0_ark, 0.0_ark, &
    !                                 -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
                                      0.0_ark, 1.0_ark, 0.0_ark, &
                                     -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    !sym%irr(5,2)%repres = reshape((/ 0.0_ark, 0.0_ark,  1.0_ark, &
    !                                -1.0_ark, 0.0_ark, 0.0_ark, &
    !                                 0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))


    !
    ! (2)  = (123)
    ! (21) = (14)*
    do irepr=3,5
      sym%irr(irepr,3)%repres = matmul(sym%irr(irepr,2)%repres,sym%irr(irepr,2)%repres)   ! (132)
      sym%irr(irepr,13)%repres= matmul(sym%irr(irepr,21)%repres,sym%irr(irepr,2)%repres)  ! (1423)*
      sym%irr(irepr,8)%repres = matmul(sym%irr(irepr,13)%repres,sym%irr(irepr,21)%repres) ! (234)
      sym%irr(irepr,5)%repres = matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,3)%repres)   ! (134)
      sym%irr(irepr,4)%repres = matmul(sym%irr(irepr,5)%repres,sym%irr(irepr,5)%repres)   ! (143)
      sym%irr(irepr,6)%repres = matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,3)%repres)   ! (142)
      sym%irr(irepr,7)%repres = matmul(sym%irr(irepr,6)%repres,sym%irr(irepr,6)%repres)   ! (124)
      sym%irr(irepr,9)%repres = matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,8)%repres)   ! (243)
      sym%irr(irepr,10)%repres= matmul(sym%irr(irepr,7)%repres,sym%irr(irepr,2)%repres)   ! (13)(24)
      sym%irr(irepr,11)%repres= matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,2)%repres)   ! (12)(34)
      sym%irr(irepr,12)%repres= matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,2)%repres)   ! (14)(23)
      sym%irr(irepr,14)%repres= matmul(sym%irr(irepr,3)%repres,sym%irr(irepr,21)%repres)  ! (1324)*
      sym%irr(irepr,15)%repres= matmul(sym%irr(irepr,11)%repres,sym%irr(irepr,21)%repres) ! (1243)*
      sym%irr(irepr,16)%repres= matmul(sym%irr(irepr,10)%repres,sym%irr(irepr,21)%repres) ! (1342)*
      sym%irr(irepr,17)%repres= matmul(sym%irr(irepr,9)%repres,sym%irr(irepr,21)%repres)  ! (1432)*
      sym%irr(irepr,18)%repres= matmul(sym%irr(irepr,2)%repres,sym%irr(irepr,21)%repres)  ! (1234)*
      sym%irr(irepr,19)%repres= matmul(sym%irr(irepr,5)%repres,sym%irr(irepr,21)%repres)  ! (13)* 
      sym%irr(irepr,20)%repres= matmul(sym%irr(irepr,7)%repres,sym%irr(irepr,21)%repres)  ! (12)*
      sym%irr(irepr,22)%repres= matmul(sym%irr(irepr,12)%repres,sym%irr(irepr,21)%repres) ! (23)*
      sym%irr(irepr,23)%repres= matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,21)%repres)  ! (34)*
      sym%irr(irepr,24)%repres= matmul(sym%irr(irepr,6)%repres,sym%irr(irepr,21)%repres)  ! (24)*
    enddo
    sym%lquant(1:2) = 0
    sym%lquant(3) = 1
    sym%lquant(4:5) = 2
    !
  case("DNH(M)","DNH") ! D_infinity_H(M)
    !
    if (mod(sym%N,2)==1) then
       !  write(out,"('symmetry: currently Dnh is only for an even ggroup order N ',i8)") sym%N
       !  stop 'symmetry: illegal order of Dnh group '
       !endif
       !
       ! Number of rotations to test for < infinity 
       !
       Nrot = sym%N ! must be >=1
       !
       ! Number of Cn classes 
       N_Cn = sym%N/2
       !
       sym%Noper=2+4*N_Cn+2*Nrot
       sym%Nclasses=4+N_Cn*2
       sym%Nrepresen= 4+N_Cn*2
       sym%CII%Noper = 0
       !
       phi = 2.0_ark*pi/real(Nrot,ark)
       !
       call simple_arrays_allocation
       !
       ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
       !
       allocate(iclass_of(sym%Noper),stat=alloc)
       if (alloc/=0) stop 'symmetry: iclass_ alloc error'
       iclass_of = 0
       !
       sym%label(1:4)=(/'A1''','A2''','A1"','A2"'/)
       !
       sym%characters(:,:) = 0
       !
       ! A1g,A1u,A2g,A2u:
       ! E
       sym%characters(1:4,1) = 1.0_ark
       ! Cinf
       sym%characters(1:4,1+1:1+N_Cn) = 1.0_ark
       ! C'inf
       sym%characters(1,1+N_Cn+1) = 1._ark
       sym%characters(2,1+N_Cn+1) =-1._ark
       sym%characters(3,1+N_Cn+1) = 1._ark
       sym%characters(4,1+N_Cn+1) =-1._ark
       ! sigmah
       sym%characters(1,1+N_Cn+2) = 1._ark
       sym%characters(2,1+N_Cn+2) = 1._ark
       sym%characters(3,1+N_Cn+2) =-1._ark
       sym%characters(4,1+N_Cn+2) =-1._ark
       ! Sinf
       sym%characters(1,1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 1._ark
       sym%characters(2,1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 1._ark
       sym%characters(3,1+N_Cn+2+1:1+N_Cn+2+N_Cn) =-1._ark
       sym%characters(4,1+N_Cn+2+1:1+N_Cn+2+N_Cn) =-1._ark
       ! sigmav
       sym%characters(1,1+N_Cn+2+N_Cn+1) = 1._ark
       sym%characters(2,1+N_Cn+2+N_Cn+1) =-1._ark
       sym%characters(3,1+N_Cn+2+N_Cn+1) =-1._ark
       sym%characters(4,1+N_Cn+2+N_Cn+1) = 1._ark
       !
       !sym%label(1:4)=(/'A1g','A2g','A1u','A2u'/)
       !
       !sym%label=(/'A1''','A2''','E'' ','A1"','A2"','E" '/)
       !
       ! E1' E1" E2' E2" E3' E3" ....
       !
       sym%lquant(1:4) = 0 
       !
       irep = 4
       do k = 1,(sym%Nrepresen-4)/2
         !
         irep = irep + 1
         !
         sym%lquant(irep  ) = k
         sym%lquant(irep+1) = k
         !
         write(Kchar, '(i4)') K
         !
         sym%label(irep  ) = 'E'//trim(adjustl(Kchar))//''''
         sym%label(irep+1) = 'E'//trim(adjustl(Kchar))//'"'
         !
         ! E 
         !
         sym%characters(irep  ,1) = 2.0_ark
         sym%characters(irep+1,1) = 2.0_ark
         !
         ! Cn
         !
         do irot = 1,N_Cn
           !
           sym%characters(irep  ,1+irot)          = 2.0_ark*cos(phi*irot*k)
           sym%characters(irep+1,1+irot)          = 2.0_ark*cos(phi*irot*k)
           !
         enddo
         !
         ! C2'
         !
         sym%characters(irep  ,1+N_Cn+2) = 0
         sym%characters(irep+1,1+N_Cn+2) = 0
         !
         ! sigmah
         !
         sym%characters(irep  ,1+N_Cn+2) = 2.0_ark
         sym%characters(irep+1,1+N_Cn+2) =-2.0_ark
         !
         do irot = 1,N_Cn
           !
           sym%characters(irep  ,1+N_Cn+2+irot)   = 2.0_ark*cos(phi*irot*k)
           sym%characters(irep+1,1+N_Cn+2+irot)   =-2.0_ark*cos(phi*irot*k)
           !
         enddo
         !
         sym%characters(irep+1,1+2*N_Cn+1) = 0 
         !
         irep = irep + 1
         !
       enddo
       !
       sym%degen(:)   = 2
       sym%degen(1:4) = 1
       !
       sym%Nelements(1) = 1
       sym%Nelements(1+ 1:1+ N_Cn) = 2
       sym%Nelements(1+N_Cn+1) = Nrot
       sym%Nelements(1+N_Cn+2) = 1
       sym%Nelements(1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 2
       sym%Nelements(1+N_Cn+2+N_Cn+1) = Nrot
       !
       o  = 0.0_ark
       p2 = 0.5_ark*pi
       p3 = 1.5_ark*pi
       !
       sym%euler(:,:) = 0
       !
       ioper = 1
       do irot = 1,N_Cn
         !
         sym%euler(1+ioper  ,:) = (/o, phi*irot,o/) ! Rz
         sym%euler(1+ioper+1,:) = (/o,-phi*irot,o/) ! Rz
         sym%euler(1+2*N_Cn+Nrot+1+ioper,:)   = (/o,pi+phi*irot,o/) ! Rz
         sym%euler(1+2*N_Cn+Nrot+1+ioper+1,:) = (/o,pi-phi*irot,o/) ! Rz
         !
         ioper = ioper + 2
       enddo
       !
       call irr_allocation
       !
       ! Generate irr-representations
       !
       do ioper = 1,sym%Noper
         !
         factor = 1.0_ark
         !
         if (ioper==1) then ! E 
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) = 1.0_ark
           sym%irr(3,ioper)%repres(1,1) = 1.0_ark
           sym%irr(4,ioper)%repres(1,1) = 1.0_ark
           !
         elseif (ioper<=1+2*N_Cn) then !  Cinf
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) = 1.0_ark
           sym%irr(3,ioper)%repres(1,1) = 1.0_ark
           sym%irr(4,ioper)%repres(1,1) = 1.0_ark
           !
         elseif (ioper<=1+2*N_Cn+Nrot) then !  C2'
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) =-1.0_ark
           sym%irr(3,ioper)%repres(1,1) = 1.0_ark
           sym%irr(4,ioper)%repres(1,1) =-1.0_ark
           !
         elseif (ioper==1+2*N_Cn+Nrot+1) then ! sigmah
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) = 1.0_ark
           sym%irr(3,ioper)%repres(1,1) =-1.0_ark
           sym%irr(4,ioper)%repres(1,1) =-1.0_ark
           !
         elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn) then !  Sinf
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) = 1.0_ark
           sym%irr(3,ioper)%repres(1,1) =-1.0_ark
           sym%irr(4,ioper)%repres(1,1) =-1.0_ark
           !
         elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then ! sigmav
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) =-1.0_ark
           sym%irr(3,ioper)%repres(1,1) =-1.0_ark
           sym%irr(4,ioper)%repres(1,1) = 1.0_ark
           !
         else
           !
           stop  'symmetry: illegal ioper'
           !
         endif
         !
       enddo
       !
       irep = 4
       do k = 1,(sym%Nrepresen-4)/2
         !
         irep = irep + 1
         !
         ioper = 1
         do ioper = 1,sym%Noper
           !
           factor = 1.0_ark
           !
           if (ioper==1) then ! E 
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) = 1.0_ark
             !
           elseif (ioper<=1+2*N_Cn) then !  Cinf
             !
             ioper_ =ioper-1 
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           elseif (ioper<=1+2*N_Cn+Nrot) then !  C2'
             !
             irot = ioper-(1+2*N_Cn)
             !
             phi_n = phi*irot*k*2.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper==1+2*N_Cn+Nrot+1) then ! sigmah
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) =-1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) =-1.0_ark
             !
           elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn) then !  Sinf
             !
             ioper_ = ioper-(1+2*N_Cn+Nrot+1)
             !
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if (mod(ioper_,2)==0)  phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) =-cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then ! sigmav
             !
             irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
             !
             phi_n  = phi*irot*k*2.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) =-cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           else
             !
             stop  'symmetry: illegal ioper'
             !
           endif
           !
         enddo
         !
         irep = irep + 1
         !
       enddo
       ! characters as traces of the corresponding representations 
       !
       do irep = 1,sym%Nrepresen
         ioper = 0
         do iclass = 1,sym%Nclasses
           do ielem =1,sym%Nelements(iclass)
             ioper = ioper + 1
             f_t = 0
             do k = 1,sym%degen(irep)
                 f_t = f_t + (sym%irr(irep,ioper)%repres(k,k))
             enddo
             sym%characters(irep,iclass) = f_t
           enddo
         enddo
       enddo
       !
       deallocate(iclass_of)
       !
       ! do the even part in the same order of operations as the odd  part
       !
    elseif (mod(sym%N,2)==0) then
       !
       Nrot = sym%N  ! Number of equivalent rotations
       NC2 = sym%N/2 ! Number of orthog. C2' axes
       !
       ! Number of Cn classes without C2
       N_Cn = sym%N/2-1
       !
       ! 1xE, 2xN_CnxCn, C2, C2', C2" ...
       !
       sym%Noper=2*(1+2*N_Cn+1+2*NC2)
       !
       ! we could alos assume Nclasses = Noper in order to be more flexible in the order of the operarations below, 
       ! which are more natuarally groupped by the similar actions <- not now
       !
       sym%Nclasses= 8+N_Cn*2 !<-classes
       sym%Nrepresen= 8+N_Cn*2
       sym%CII%Noper = 0
       !
       f_t = (-1.0_ark)**(Nrot/2)
       !
       phi = 2.0_ark*pi/real(Nrot,ark)
       !
       call simple_arrays_allocation
       !
       allocate(iclass_of(sym%Noper),stat=alloc)
       if (alloc/=0) stop 'symmetry: iclass_ alloc error'
       iclass_of = 0
       !
       ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
       !
       !sym%label(1:4)=(/'A1g','A2g','A1u','A2u'/)
       !
       sym%label(1:8)=(/'A1g','A2g','B1g','B2g','A1u','A2u','B1u','B2u'/)
       !
       ! Characters
       !
       sym%characters(:,:) = 0.0_ark
       !
       ! A1g,A2g,B1g,B2g,A1u,A2u,B1u,B2u:
       ! E
       sym%characters(1:8,1) = 1.0_ark
       !
       ! Cn - rotations
       !
       do icn=1,N_Cn
         !
         sym%characters(1,1+icn) = 1.0_ark
         sym%characters(2,1+icn) = 1.0_ark
         sym%characters(3,1+icn) =(-1.0_ark)**icn
         sym%characters(4,1+icn) =(-1.0_ark)**icn
         sym%characters(5,1+icn) = 1.0_ark
         sym%characters(6,1+icn) = 1.0_ark
         sym%characters(7,1+icn) =(-1.0_ark)**icn
         sym%characters(8,1+icn) =(-1.0_ark)**icn
         !
       enddo
       !
       ! C2
       sym%characters(1,2+N_Cn) = 1.0_ark
       sym%characters(2,2+N_Cn) = 1.0_ark
       sym%characters(3,2+N_Cn) = f_t
       sym%characters(4,2+N_Cn) = f_t
       sym%characters(5,2+N_Cn) = 1.0_ark
       sym%characters(6,2+N_Cn) = 1.0_ark
       sym%characters(7,2+N_Cn) = f_t
       sym%characters(8,2+N_Cn) = f_t
       ! nxC2'
       sym%characters(1,3+N_Cn) = 1.0_ark
       sym%characters(2,3+N_Cn) =-1.0_ark
       sym%characters(3,3+N_Cn) = 1.0_ark
       sym%characters(4,3+N_Cn) =-1.0_ark
       sym%characters(5,3+N_Cn) = 1.0_ark
       sym%characters(6,3+N_Cn) =-1.0_ark
       sym%characters(7,3+N_Cn) = 1.0_ark
       sym%characters(8,3+N_Cn) =-1.0_ark
       ! nxC2"
       sym%characters(1,4+N_Cn) = 1.0_ark
       sym%characters(2,4+N_Cn) =-1.0_ark
       sym%characters(3,4+N_Cn) =-1.0_ark
       sym%characters(4,4+N_Cn) = 1.0_ark
       sym%characters(5,4+N_Cn) = 1.0_ark
       sym%characters(6,4+N_Cn) =-1.0_ark
       sym%characters(7,4+N_Cn) =-1.0_ark
       sym%characters(8,4+N_Cn) = 1.0_ark
       ! i
       sym%characters(1,5+N_Cn) = 1.0_ark
       sym%characters(2,5+N_Cn) = 1.0_ark
       sym%characters(3,5+N_Cn) = 1.0_ark
       sym%characters(4,5+N_Cn) = 1.0_ark
       sym%characters(5,5+N_Cn) =-1.0_ark
       sym%characters(6,5+N_Cn) =-1.0_ark
       sym%characters(7,5+N_Cn) =-1.0_ark
       sym%characters(8,5+N_Cn) =-1.0_ark
       !
       ! Sn
       do icn=1,N_Cn
         sym%characters(1,5+N_Cn+icn) = 1.0_ark
         sym%characters(2,5+N_Cn+icn) = 1.0_ark
         sym%characters(3,5+N_Cn+icn) =(-1.0_ark)**icn*f_t
         sym%characters(4,5+N_Cn+icn) =(-1.0_ark)**icn*f_t
         sym%characters(5,5+N_Cn+icn) =-1.0_ark
         sym%characters(6,5+N_Cn+icn) =-1.0_ark
         sym%characters(7,5+N_Cn+icn) =-(-1.0_ark)**icn*f_t
         sym%characters(8,5+N_Cn+icn) =-(-1.0_ark)**icn*f_t
       enddo
       ! sigmah
       sym%characters(1,6+2*N_Cn) = 1.0_ark
       sym%characters(2,6+2*N_Cn) = 1.0_ark
       sym%characters(3,6+2*N_Cn) = f_t
       sym%characters(4,6+2*N_Cn) = f_t
       sym%characters(5,6+2*N_Cn) =-1.0_ark
       sym%characters(6,6+2*N_Cn) =-1.0_ark
       sym%characters(7,6+2*N_Cn) =-f_t
       sym%characters(8,6+2*N_Cn) =-f_t
       ! sigmav
       sym%characters(1,7+2*N_Cn) = 1.0_ark
       sym%characters(2,7+2*N_Cn) =-1.0_ark
       sym%characters(3,7+2*N_Cn) = f_t
       sym%characters(4,7+2*N_Cn) =-f_t
       sym%characters(5,7+2*N_Cn) =-1.0_ark
       sym%characters(6,7+2*N_Cn) = 1.0_ark
       sym%characters(7,7+2*N_Cn) =-f_t
       sym%characters(8,7+2*N_Cn) = f_t
       ! sigmad
       sym%characters(1,8+2*N_Cn) = 1.0_ark
       sym%characters(2,8+2*N_Cn) =-1.0_ark
       sym%characters(3,8+2*N_Cn) =-f_t
       sym%characters(4,8+2*N_Cn) = f_t
       sym%characters(5,8+2*N_Cn) =-1.0_ark
       sym%characters(6,8+2*N_Cn) = 1.0_ark
       sym%characters(7,8+2*N_Cn) = f_t
       sym%characters(8,8+2*N_Cn) =-f_t
       !
       !
       irep = 8
       do k = 1,(sym%Nrepresen-8)/2
         !
         irep = irep + 1
         !
         write(Kchar, '(i4)') K
         !
         sym%label(irep  ) = 'E'//trim(adjustl(Kchar))//'g'
         sym%label(irep+1) = 'E'//trim(adjustl(Kchar))//'u'
         !
         if (k>=10) then
           sym%label(irep  ) = trim(adjustl(Kchar))//'g'
           sym%label(irep+1) = trim(adjustl(Kchar))//'u'
         endif
         !
         sym%characters(irep  ,1) = 2.0_ark
         sym%characters(irep+1,1) = 2.0_ark
         !
         sym%characters(irep  ,5+N_Cn) = 2.0_ark
         sym%characters(irep+1,5+N_Cn) =-2.0_ark
         !
         sym%lquant(irep  ) = k
         sym%lquant(irep+1) = k
         !
         do icn = 1,N_Cn
           !
           sym%characters(irep  ,1+icn)          = 2.0_ark*cos(phi*icn*k)
           sym%characters(irep+1,1+icn)          = 2.0_ark*cos(phi*icn*k)
           !
           sym%characters(irep  ,5+N_Cn+icn)   =-2.0_ark*cos(phi*icn*k)
           sym%characters(irep+1,5+N_Cn+icn)   = 2.0_ark*cos(phi*icn*k)
           !
         enddo
         !
         sym%characters(irep  ,2+N_Cn) = 2.0_ark*(-1)**irep
         sym%characters(irep+1,2+N_Cn) = 2.0_ark*(-1)**irep
         !
         sym%characters(irep  ,6+2*N_Cn) = 2.0_ark*(-1)**irep
         sym%characters(irep+1,6+2*N_Cn) =-2.0_ark*(-1)**irep
         !
         irep = irep + 1
         !
       enddo
       !
       !sym%label=(/'A1''','A2''','E'' ','A1"','A2"','E" '/)
       !
       ! E1' E1" E2' E2" E3' E3" ....
       !
       sym%lquant(1:4) = 0 
       !
       !
       sym%degen(:)   = 2
       sym%degen(1:8) = 1
       !
       sym%Nelements(1) = 1
       sym%Nelements(1+ 1:1+ N_Cn) = 2
       sym%Nelements(2+N_Cn) = 1
       sym%Nelements(3+N_Cn) = NC2
       sym%Nelements(4+N_Cn) = NC2
       sym%Nelements(5+N_Cn) = 1
       sym%Nelements(5+N_Cn+1:5+2*N_Cn) = 2
       sym%Nelements(6+2*N_Cn) = 1
       sym%Nelements(7+2*N_Cn) = NC2
       sym%Nelements(8+2*N_Cn) = NC2
       !
       ! Define elements in  classes
       !
       sym%Nelements = 0
       !
       ! E
       sym%Nelements(1) = 1
       !
       ! Cn
       do ioper = 1,N_Cn
         iclass = 1 + ioper
         sym%Nelements(iclass) = 2
       enddo 
       !
       ! C2
       sym%Nelements(1+N_Cn+1) = 1
       !
       ! C2'
       sym%Nelements(1+N_Cn+2) = NC2
       !
       ! C2"
       sym%Nelements(1+N_Cn+3) = NC2
       !
       ! i
       sym%Nelements(1+N_Cn+4) = 1
       !
       ! Sn
       do ioper = 1,N_Cn
         iclass = 1+N_Cn+4+ioper
         sym%Nelements(iclass) = 2
       enddo 
       !
       iclass = 5+N_Cn+N_Cn
       !
       ! sigmah
       sym%Nelements(iclass+1) = 1
       !
       ! sigmav
       sym%Nelements(iclass+2) = NC2
       !
       ! sigmad
       sym%Nelements(iclass+3) = NC2
       !
       !
       o  = 0.0_ark
       p2 = 0.5_ark*pi
       p3 = 1.5_ark*pi
       !
       sym%euler(:,:) = 0
       !
       ioper = 1
       do irot = 1,N_Cn
         !
         sym%euler(1+ioper  ,:) = (/o, phi*irot,o/) ! Rz
         sym%euler(1+ioper+1,:) = (/o,-phi*irot,o/) ! Rz
         sym%euler(1+2*N_Cn+Nrot+1+ioper,:)   = (/o,pi+phi*irot,o/) ! Rz
         sym%euler(1+2*N_Cn+Nrot+1+ioper+1,:) = (/o,pi-phi*irot,o/) ! Rz
         !
         ioper = ioper + 2
       enddo
       !
       call irr_allocation
       !
       do ioper = 1,sym%Noper
          !
          if (ioper==1) then ! E 
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) = 1.0_ark
            sym%irr(3,ioper)%repres(1,1) = 1.0_ark
            sym%irr(4,ioper)%repres(1,1) = 1.0_ark
            sym%irr(5,ioper)%repres(1,1) = 1.0_ark
            sym%irr(6,ioper)%repres(1,1) = 1.0_ark
            sym%irr(7,ioper)%repres(1,1) = 1.0_ark
            sym%irr(8,ioper)%repres(1,1) = 1.0_ark
            !
          elseif (ioper<=1+2*N_Cn) then ! Cn x 2 x(n/2-1)
            !
            ioper_ =(ioper)/2 
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) = 1.0_ark
            sym%irr(3,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            sym%irr(4,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            sym%irr(5,ioper)%repres(1,1) = 1.0_ark
            sym%irr(6,ioper)%repres(1,1) = 1.0_ark
            sym%irr(7,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            sym%irr(8,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            !
          elseif (ioper<=1+2*N_Cn+1) then !  C2 only once
            !
            ioper_ =(ioper)/2 
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) = 1.0_ark
            sym%irr(3,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            sym%irr(4,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            sym%irr(5,ioper)%repres(1,1) = 1.0_ark
            sym%irr(6,ioper)%repres(1,1) = 1.0_ark
            sym%irr(7,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            sym%irr(8,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            !
          elseif (ioper<=2+2*N_Cn+NC2) then !  C2'
            !
            ioper_ = ioper-(1+2*N_Cn+1)
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) =-1.0_ark
            sym%irr(3,ioper)%repres(1,1) = 1.0_ark
            sym%irr(4,ioper)%repres(1,1) =-1.0_ark
            sym%irr(5,ioper)%repres(1,1) = 1.0_ark
            sym%irr(6,ioper)%repres(1,1) =-1.0_ark
            sym%irr(7,ioper)%repres(1,1) = 1.0_ark
            sym%irr(8,ioper)%repres(1,1) =-1.0_ark
            !
          elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
            !
            ioper_ = ioper-(1+2*N_Cn+NC2)
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) =-1.0_ark
            sym%irr(3,ioper)%repres(1,1) =-1.0_ark
            sym%irr(4,ioper)%repres(1,1) = 1.0_ark
            sym%irr(5,ioper)%repres(1,1) = 1.0_ark
            sym%irr(6,ioper)%repres(1,1) =-1.0_ark
            sym%irr(7,ioper)%repres(1,1) =-1.0_ark
            sym%irr(8,ioper)%repres(1,1) = 1.0_ark
            !
          elseif (ioper==3+2*N_Cn+2*NC2) then ! i
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) = 1.0_ark
            sym%irr(3,ioper)%repres(1,1) = 1.0_ark
            sym%irr(4,ioper)%repres(1,1) = 1.0_ark
            sym%irr(5,ioper)%repres(1,1) =-1.0_ark
            sym%irr(6,ioper)%repres(1,1) =-1.0_ark
            sym%irr(7,ioper)%repres(1,1) =-1.0_ark
            sym%irr(8,ioper)%repres(1,1) =-1.0_ark
            !
          elseif (ioper<=3+4*N_Cn+2*NC2) then !  2xnxSn
            !
            ioper_ =ioper-(3+2*N_Cn+2*NC2)
            !
            ioper_ =(ioper_+1)/2+NC2
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) = 1.0_ark
            sym%irr(3,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            sym%irr(4,ioper)%repres(1,1) = (-1.0_ark)**(ioper_)
            sym%irr(5,ioper)%repres(1,1) =-1.0_ark
            sym%irr(6,ioper)%repres(1,1) =-1.0_ark
            sym%irr(7,ioper)%repres(1,1) = (-1.0_ark)**(ioper_+1)
            sym%irr(8,ioper)%repres(1,1) = (-1.0_ark)**(ioper_+1)
            !
          elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigma_h
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) = 1.0_ark
            sym%irr(3,ioper)%repres(1,1) =(-1.0_ark)**(NC2)
            sym%irr(4,ioper)%repres(1,1) =(-1.0_ark)**(NC2)
            sym%irr(5,ioper)%repres(1,1) =-1.0_ark
            sym%irr(6,ioper)%repres(1,1) =-1.0_ark
            sym%irr(7,ioper)%repres(1,1) =(-1.0_ark)**(NC2+1)
            sym%irr(8,ioper)%repres(1,1) =(-1.0_ark)**(NC2+1)
            !
          elseif (ioper<=4+4*N_Cn+3*NC2) then !  sigmav: odd rotations
            !
            ioper_ =ioper-(4+4*N_Cn+2*NC2)+NC2
            !
            ioper_ = NC2
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) =-1.0_ark
            sym%irr(3,ioper)%repres(1,1) = (-1.0_ark)**(ioper_  )
            sym%irr(4,ioper)%repres(1,1) = (-1.0_ark)**(ioper_+1)
            sym%irr(5,ioper)%repres(1,1) =-1.0_ark
            sym%irr(6,ioper)%repres(1,1) = 1.0_ark
            sym%irr(7,ioper)%repres(1,1) = (-1.0_ark)**(ioper_+1)
            sym%irr(8,ioper)%repres(1,1) = (-1.0_ark)**(ioper_  )
            !
          elseif (ioper<=4+4*N_Cn+4*NC2) then !  sigmad
            !
            ioper_ =ioper-(4+4*N_Cn+3*NC2)+NC2
            !
            ioper_ = NC2+1
            !
            sym%irr(1,ioper)%repres(1,1) = 1.0_ark
            sym%irr(2,ioper)%repres(1,1) =-1.0_ark
            sym%irr(3,ioper)%repres(1,1) = (-1.0_ark)**(ioper_  )
            sym%irr(4,ioper)%repres(1,1) = (-1.0_ark)**(ioper_+1)
            sym%irr(5,ioper)%repres(1,1) =-1.0_ark
            sym%irr(6,ioper)%repres(1,1) = 1.0_ark
            sym%irr(7,ioper)%repres(1,1) = (-1.0_ark)**(ioper_+1)
            sym%irr(8,ioper)%repres(1,1) = (-1.0_ark)**(ioper_  )
            !
          else
            !
            stop  'symmetry: illegal ioper'
            !
          endif
          !
       enddo
       !
       ! Generate irr-representations
       !
       irep = 8
       do k = 1,(sym%Nrepresen-8)/2
         !
         irep = irep + 1
         !
         sym%lquant(irep  ) = k
         sym%lquant(irep+1) = k
         !
         do ioper = 1,sym%Noper
           !
           factor = 1.0_ark
           !
           if (ioper==1) then ! E !! 
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) = 1.0_ark
             !
           elseif (ioper<=1+2*N_Cn) then ! Cn x 2 x(n/2-1)
             !
             ioper_ =ioper-1 
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           elseif (ioper<=2+2*N_Cn) then !  C2 only once
             !
             ioper_ =ioper-1
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
             !
           elseif (ioper<=2+2*N_Cn+NC2) then !  C2'
             !
             irot =ioper-(2+2*N_Cn)-1
             !
             phi_n = phi*irot*2.0_ark*k
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper<=2+2*N_Cn+2*NC2) then !  C2"
             !
             irot =ioper-(2+2*N_Cn+NC2)-1
             !
             phi_n = phi*(2*irot+1)*k
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
             !
           elseif (ioper==3+2*N_Cn+2*NC2) then ! i
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
             sym%irr(irep+1,ioper)%repres(1,1) =-1.0_ark
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) =-1.0_ark
             !
           elseif (ioper<=3+4*N_Cn+2*NC2) then !  2xnxSn
             !
             ioper_ =ioper-(3+2*N_Cn+2*NC2)
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)*(-1.0_ark)**k
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)*(-1.0_ark)**k
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)*(-1.0_ark)**k
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)*(-1.0_ark)**k
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)*(-1.0_ark)**(k+1)
             sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)*(-1.0_ark)**(k+1)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)*(-1.0_ark)**(k+1)
             sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)*(-1.0_ark)**(k+1)
             !
           elseif (ioper<=4+4*N_Cn+2*NC2) then !  sigma_h
             !
             ! (Snxn) for phin=pi
             !
             sym%irr(irep,ioper)%repres(1,1) = (-1.0_ark)**k
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = (-1.0_ark)**k
             !
             sym%irr(irep+1,ioper)%repres(1,1) = (-1.0_ark)**(k+1)
             sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep+1,ioper)%repres(2,2) = (-1.0_ark)**(k+1)
             !
           elseif (ioper<=4+4*N_Cn+3*NC2) then !  sigmav -> reflection through a line making angle phi is equivalent to
             ! cos(2phi) sin(2phi)
             ! sin(2phi) -cos(2phi)
             !
             irot = ioper-(4+4*N_Cn+2*NC2)-1
             !
             phi_n = phi*irot*k*2.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)*(-1.0_ark)**k
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)*(-1.0_ark)**k
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)*(-1.0_ark)**k
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)*(-1.0_ark)**k
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)*(-1.0_ark)**(k+1)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)*(-1.0_ark)**(k+1)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)*(-1.0_ark)**(k+1)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)*(-1.0_ark)**(k+1)
             !
           elseif (ioper<=4+4*N_Cn+4*NC2) then ! sigmad
             !
             irot = ioper-(4+4*N_Cn+3*NC2)-1
             !
             !phi_n = (-phi*0.5+phi*irot)*2.0_ark
             !
             phi_n = phi*(2*irot+1)*k
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)*(-1.0_ark)**k
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)*(-1.0_ark)**k
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)*(-1.0_ark)**k
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)*(-1.0_ark)**k
             !
             sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)*(-1.0_ark)**(k+1)
             sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)*(-1.0_ark)**(k+1)
             !
             sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)*(-1.0_ark)**(k+1)
             sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)*(-1.0_ark)**(k+1)
             !
           else
             !
             stop  'symmetry: illegal ioper'
             !
           endif
           !
         enddo
         !
         irep = irep + 1
         !
       enddo
       !
       ! characters as traces of the corresponding representations 
       !
       do irep = 1,sym%Nrepresen
         ioper = 0
         do iclass = 1,sym%Nclasses
           do ielem =1,sym%Nelements(iclass)
             ioper = ioper + 1
             f_t = 0
             do k = 1,sym%degen(irep)
                 f_t = f_t + (sym%irr(irep,ioper)%repres(k,k))
             enddo
             sym%characters(irep,iclass) = f_t
           enddo
         enddo
       enddo
       !
       !do ioper =1,sym%Noper
       !  !
       !  iclass = iclass_of(ioper)
       !  !
       !  do joper =ioper+1,sym%Noper
       !    !
       !    jclass = iclass_of(joper)
       !    !
       !    if (iclass>jclass) then
       !      !
       !      iclass_of(joper) = iclass
       !      iclass_of(ioper) = jclass
       !      !
       !      do irep = 1,8
       !        !
       !        f_t= sym%irr(irep,joper)%repres(1,1)
       !        sym%irr(irep,joper)%repres(1,1) = sym%irr(irep,ioper)%repres(1,1)
       !        sym%irr(irep,ioper)%repres(1,1) = f_t
       !        !
       !      enddo
       !      !
       !      do irep = 9,sym%Nrepresen
       !        !
       !        repres_= sym%irr(irep,joper)%repres
       !        sym%irr(irep,joper)%repres = sym%irr(irep,ioper)%repres
       !        sym%irr(irep,ioper)%repres = repres_
       !        !
       !      enddo
       !      !
       !    endif
       !    !
       !  enddo
       !enddo
       !
       deallocate(iclass_of)
       !
       !do irep = 9,sym%Nrepresen
       !  !
       !  mat_t  = matmul(sym%irr(irep,2)%repres,sym%irr(irep,2)%repres) ! C4xC4 = C2
       !  mat_t  = matmul(sym%irr(irep,2)%repres,mat_t)                  ! C4xC2 = C4-2
       !  mat_t  = matmul(sym%irr(irep,2)%repres,sym%irr(irep,5)%repres) ! C2" = C4(2)xC2'(1)
       !  mat_t  = matmul(sym%irr(irep,3)%repres,sym%irr(irep,5)%repres) !
       !  mat_t  = matmul(sym%irr(irep,4)%repres,sym%irr(irep,12)%repres)-sym%irr(irep, 9)%repres ! i
       !  mat_t  = matmul(sym%irr(irep,2)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,10)%repres ! Sn(2)
       !  mat_t  = matmul(sym%irr(irep,3)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,11)%repres ! Sn(2)!
       !  mat_t  = matmul(sym%irr(irep,3)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,10)%repres ! Sn(2)
       !  mat_t  = matmul(sym%irr(irep,2)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,11)%repres ! Sn(2)
       !  mat_t  = matmul(sym%irr(irep,5)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,13)%repres ! sigmav1
       !!  mat_t  = matmul(sym%irr(irep,6)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,14)%repres ! sigmav2 
       !  mat_t  = matmul(sym%irr(irep,7)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,15)%repres ! sigmad1
       !  mat_t  = matmul(sym%irr(irep,8)%repres,sym%irr(irep,12)%repres)-sym%irr(irep,16)%repres ! sigmad2
       !  !
       !  continue
       !  !
       !enddo
       !
    endif
    !
  case("CNV(M)","CNV") ! C_infinity_V(M)
    !
    if (mod(sym%N,2)==1) then
       !
       ! Number of rotations 
       !
       Nrot = sym%N ! must be >=1
       !
       ! Number of Cn classes 
       N_Cn = sym%N/2
       !
       sym%Noper=1+2*N_Cn+Nrot
       sym%Nclasses=2+N_Cn
       sym%Nrepresen= 2+N_Cn
       sym%CII%Noper = 0
       !
       phi = 2.0_ark*pi/real(Nrot,ark)
       !
       call simple_arrays_allocation
       !
       ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
       !
       allocate(iclass_of(sym%Noper),stat=alloc)
       if (alloc/=0) stop 'symmetry: iclass_ alloc error'
       iclass_of = 0
       !
       sym%label(1:2)=(/'A1','A2'/)
       !
       sym%characters(:,:) = 0
       !
       ! A1g,A1u,A2g,A2u:
       ! E
       sym%characters(1:2,1) = 1.0_ark
       ! Cinf
       sym%characters(1:2,1+N_Cn) = 1.0_ark
       ! sigmav
       sym%characters(1,1+N_Cn+1) = 1._ark
       sym%characters(2,1+N_Cn+1) =-1._ark
       !
       ! E1' E1" E2' E2" E3' E3" ....
       !
       sym%lquant(1:2) = 0 
       !
       irep = 2
       do k = 1,sym%Nrepresen-2
         !
         irep = k + 2
         !
         sym%lquant(irep  ) = k
         !
         write(Kchar, '(i4)') K
         !
         sym%label(irep  ) = 'E'//trim(adjustl(Kchar))
         !
         ! E 
         !
         sym%characters(irep  ,1) = 2.0_ark
         !
         ! Cn
         !
         do irot = 1,N_Cn
           !
           sym%characters(irep  ,1+irot)          = 2.0_ark*cos(phi*irot*k)
           !
         enddo
         !
         sym%characters(irep,1+N_Cn+1) = 0 
         !
       enddo
       !
       sym%degen(:)   = 2
       sym%degen(1:2) = 1
       !
       sym%Nelements(1) = 1
       sym%Nelements(1+1:1+ N_Cn) = 2
       sym%Nelements(1+N_Cn+1) = Nrot
       !
       o  = 0.0_ark
       p2 = 0.5_ark*pi
       p3 = 1.5_ark*pi
       !
       sym%euler(:,:) = 0
       !
       ioper = 1
       do irot = 1,N_Cn
         !
         sym%euler(1+ioper  ,:) = (/o, phi*irot,o/) ! Rz
         sym%euler(1+ioper+1,:) = (/o,-phi*irot,o/) ! Rz
         !
         ioper = ioper + 2
       enddo
       !
       call irr_allocation
       !
       ! Generate irr-representations
       !
       do ioper = 1,sym%Noper
         !
         factor = 1.0_ark
         !
         if (ioper==1) then ! E 
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) = 1.0_ark
           !
         elseif (ioper<=1+2*N_Cn) then !  Cinf
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) = 1.0_ark
           !
         elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then ! sigmav
           !
           sym%irr(1,ioper)%repres(1,1) = 1.0_ark
           sym%irr(2,ioper)%repres(1,1) =-1.0_ark
           !
         else
           !
           stop  'symmetry: illegal ioper'
           !
         endif
         !
       enddo
       !
       irep = 2
       do k = 1,sym%Nrepresen-2
         !
         irep = k+2
         !
         ioper = 1
         do ioper = 1,sym%Noper
           !
           factor = 1.0_ark
           !
           if (ioper==1) then ! E 
             !
             sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
             sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
             !
             sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
             sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
             !
           elseif (ioper<=1+2*N_Cn) then !  Cinf
             !
             ioper_ =ioper-1 
             irot = (ioper_+1)/2
             !
             phi_n = phi*irot*k
             !
             ! Second oper in a class is with negative phi
             if ( mod(ioper_,2)==0 ) phi_n = -phi_n
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
             !
           elseif (ioper<=1+2*N_Cn+Nrot) then ! sigmav
             !
             irot = ioper-(1+2*N_Cn)
             !
             phi_n  = phi*irot*k*2.0_ark
             !
             sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
             sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
             !
             sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
             sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
             !
           else
             !
             stop  'symmetry: illegal ioper'
             !
           endif
           !
         enddo
         !
       enddo
       ! characters as traces of the corresponding representations 
       !
       do irep = 1,sym%Nrepresen
         ioper = 0
         do iclass = 1,sym%Nclasses
           do ielem =1,sym%Nelements(iclass)
             ioper = ioper + 1
             f_t = 0
             do k = 1,sym%degen(irep)
                 f_t = f_t + (sym%irr(irep,ioper)%repres(k,k))
             enddo
             sym%characters(irep,iclass) = f_t
           enddo
         enddo
       enddo
       !
       deallocate(iclass_of)
       !
       ! do the even part in the same order of operations as the odd  part
       !
    elseif (mod(sym%N,2)==0) then
       !
       write(out,"('CNV for N-even has not been implemented yet, not working')")
       stop 'CNV for N-even has not been implemented yet, not working'
       !
    endif
    !
  case("DINFTYH(M)") ! D_infinity_H(M)

    ! Number of rotations to test for < infinity 
    !
    Nrot = 37 ! must be >=1
    !
    sym%Noper=6+2*Nrot
    sym%Nclasses=6
    sym%Nrepresen= max_irreps  ! must be even and >=4
    sym%CII%Noper = 0
    !
    phi = 2.0_ark*pi/real(Nrot,ark)
    !
    call simple_arrays_allocation
    !
    ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
    !
    sym%characters(:,:) = 0
    !
    ! A1g,A1u,A2g,A2u:
    ! E
    sym%characters(1:4,1) = 1
    ! Cinf
    sym%characters(1:4,2) = 1
    ! sigmav
    sym%characters(1,3) = 1
    sym%characters(2,3) = 1
    sym%characters(3,3) =-1
    sym%characters(4,3) =-1
    ! i
    sym%characters(1,4) = 1
    sym%characters(2,4) =-1
    sym%characters(3,4) = 1
    sym%characters(4,4) =-1
    ! Sinf
    sym%characters(1,5) = 1
    sym%characters(2,5) =-1
    sym%characters(3,5) = 1
    sym%characters(4,5) =-1
    ! C'inf
    sym%characters(1,6) = 1
    sym%characters(2,6) =-1
    sym%characters(3,6) =-1
    sym%characters(4,6) = 1
    !
    sym%label(1:4)=(/'A1g','A2g','A1u','A2u'/)
    !
    ! E_1g E_1u E_2g E_2u E_3g E_3u .... 
    !
    irep = 4
    do k = 1,(sym%Nrepresen-4)/2
      !
      irep = irep + 1
      !
      write(Kchar, '(i4)') K
      !
      sym%label(irep  ) = 'E'//trim(adjustl(Kchar))//'g'
      sym%label(irep+1) = 'E'//trim(adjustl(Kchar))//'u'
      !
      sym%characters(irep  ,1) = 2.0_ark
      sym%characters(irep+1,1) = 2.0_ark
      !
      sym%characters(irep  ,4) = 2.0_ark
      sym%characters(irep+1,4) =-2.0_ark
      !
      !do irot = 1,Nrot
      !  sym%characters(irep  ,2)          = 2.0_ark*cos(phi*irot*k)
      !  sym%characters(irep+1,2)          = 2.0_ark*cos(phi*irot*k)
      !  sym%characters(irep  ,5) = 2.0_ark*cos(phi*irot*k)*(-1)**k
      !  sym%characters(irep+1,5) = 2.0_ark*cos(phi*irot*k)*(-1)**(k+1)
      !enddo
      !
      sym%characters(irep  ,2)          = 2.0_ark*cos(phi*k)
      sym%characters(irep+1,2)          = 2.0_ark*cos(phi*k)
      sym%characters(irep  ,5) = 2.0_ark*cos(phi*k)*(-1)**k
      sym%characters(irep+1,5) = 2.0_ark*cos(phi*k)*(-1)**(k+1)
      !
      irep = irep + 1
      !
    enddo
    !
    sym%degen(:)   = 2
    sym%degen(1:4) = 1
    !
    sym%Nelements = (/1,2,Nrot,1,2,Nrot/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler(:,:) = 0
    !
    !do irot = 1,Nrot
      sym%euler(2,:) = (/o, phi,o/) ! Rz
      sym%euler(3,:) = (/o,-phi,o/) ! Rz
      sym%euler(3+Nrot+2,:) = (/o,pi+phi,o/) ! Rz
      sym%euler(3+Nrot+2,:) = (/o,pi-phi,o/) ! Rz
    !enddo
    !
    call irr_allocation
    !
  case default

    write(out,"('symmetry: undefined symmetry group ',a)") trim(sym_group)
    stop 'symmetry: undefined symmetry group '

  end select
  !
  !if (max_irreps<sym%Nrepresen) then 
  !  !
  !  write(out,"('symmetry: max_irreps is too small: ',i5,' increase to > ',i5)") max_irreps,sym%Nrepresen
  !  stop 'symmetry: size of _select_gamma_ is too small'
  !  !
  !endif 
  !
  sym%maxdegen = maxval(sym%degen(:),dim=1)

  !
  ! store the address of the group generator from ioper = 1..Noper list 
  !
  ioper = 1
  !
  do iclass = 1,sym%Nclasses
    !
    sym%igenerator(iclass) = ioper
    ioper = ioper + sym%Nelements(iclass)
    !
  enddo
  
  !
  ! check 
  call check_characters_and_representation
  !

  contains

   subroutine simple_arrays_allocation

    integer(ik) :: alloc,nCII

    nCII = max(1,sym%CII%Noper)
    !
    allocate (sym%characters(sym%Nrepresen,sym%Nclasses),sym%irr(sym%Nrepresen,sym%Noper),&  
              sym%degen(sym%Nrepresen),sym%Nelements(sym%Nclasses),sym%label(sym%Nrepresen),&
              sym%igenerator(sym%Nclasses),&
              sym%CII%ioper(nCII),sym%CII%coeff(nCII),sym%euler(sym%Noper,3),sym%lquant(sym%Nrepresen),stat=alloc)

    if (alloc/=0) stop 'simple_arrays_allocation - out of memory'
    !
    sym%CII%coeff = 0
    sym%CII%ioper = 1
    sym%euler = 0
    sym%lquant = 0 
    !
   end subroutine simple_arrays_allocation



   subroutine irr_allocation

    integer(ik) :: gamma,ioper,iclass,ielem,alloc

    do gamma = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem = 1,sym%Nelements(iclass)
          !
          ioper = ioper + 1
          !
          allocate(sym%irr(gamma,ioper)%repres(sym%degen(gamma),sym%degen(gamma)),stat=alloc)
          !
          if (sym%degen(gamma)==1) then
             sym%irr(gamma,ioper)%repres(1,1)= sym%characters(gamma,iclass)
          endif 
          !
        enddo 
      enddo 
    enddo 


   if (alloc/=0) then
       write (out,"(' symmetryInitialize ',i9,' error trying to allocate symmetry')") alloc
      stop 'symmetryInitialize, symmetries - out of memory'
   end if


   end subroutine irr_allocation


 end subroutine symmetryInitialize



   subroutine check_characters_and_representation

    integer(ik) :: igamma,jgamma,ioper,ideg,iclass,ielem,k1,k2,m1,m2
    real(ark)   :: temp,f_t
    character(len=10) :: my_fmt !format for I/O specification

    do igamma = 1,sym%Nrepresen
      do jgamma = igamma,sym%Nrepresen
        !
        temp = sum(real(sym%Nelements(:),ark)*sym%characters(igamma,:)*sym%characters(jgamma,:))
        !
        if (igamma/=jgamma.and.abs(temp)>sqrt(small_)) then 
          write (out,"(' check_char_and_repres: not orhogonal for igamma,jgamma = ',2i4,' -> ',g18.6)") igamma,jgamma,temp
          !stop 'check_characters_and_representation: not orhogonal'
        endif
        !
        if (igamma==jgamma.and.abs(temp-sym%Noper)>sqrt(small_)) then 
          write (out,"(' check_charac_and_repres: dot product ',f16.2,' for isym = ',i4,' /= size of the group ',f16.0)") & 
                temp,igamma,sym%Noper
          !stop 'check_characters_and_representation: not orhogonal'
        endif
        !
      enddo 
    enddo 

    !
    ! Check characters
    !
    if (verbose_>=5) write(out,"('Irrep matrices:')")
    !
    ioper = 0
    do iclass = 1,sym%Nclasses
      !
      do ielem = 1,sym%Nelements(iclass)
        !
        ioper = ioper + 1
        !
        do igamma = 1,sym%Nrepresen
          !
          temp = 0
          !
          do ideg = 1,sym%degen(igamma)
            !
            temp = temp + sym%irr(igamma,ioper)%repres(ideg,ideg)
            !
          enddo
          !
          if (verbose_>=5) then
            write(out,"('igamma,iclass,ioper = ',3i6)") igamma,iclass,ioper
            do ideg = 1,sym%degen(igamma)
              write(my_fmt,'(a,i0,a)') "(",sym%degen(igamma),"f18.8)"
              write(out,my_fmt) sym%irr(igamma,ioper)%repres(ideg,:)
            enddo
          endif
          !
          if (abs(temp-sym%characters(igamma,iclass))>sqrt(small_)) then 
            write (out,"(' symmetry: character and representation do not agree ',2f18.8,', igamma,iclass,ioper = ',3i5)") & 
                  sym%characters(igamma,iclass),temp,igamma,iclass,ioper
            stop 'symmetry: character and representation do not agree'
          endif
          !
        enddo
        !
      enddo
      !
    enddo

    do igamma = 1,sym%Nrepresen
      do jgamma = igamma,sym%Nrepresen
        do k1 = 1,sym%degen(igamma)
          do k2 = 1,sym%degen(igamma)
            do m1 = 1,sym%degen(jgamma)
              do m2 = 1,sym%degen(jgamma)
                ioper = 0
                f_t = 0
                do iclass = 1,sym%Nclasses
                  do ielem =1,sym%Nelements(iclass)
                   ioper = ioper + 1
                   f_t = f_t + sym%irr(igamma,ioper)%repres(k1,k2)*sqrt(real(sym%degen(igamma))/real(sym%Noper,ark))*&
                               sym%irr(jgamma,ioper)%repres(m1,m2)*sqrt(real(sym%degen(jgamma))/real(sym%Noper,ark))
                  enddo          
                enddo
                !
                if ((k1/=m1.or.k2/=m2.or.igamma/=jgamma).and.abs(f_t)>small_) then
                  write(out,"('Non orthogonal irreps for igamma,jgamma,k1,k2,m1,m2 = ',6i7,f16.7)") igamma,jgamma,k1,k2,m1,m2,f_t
                  stop 'Non orthogonal irreps'
                endif
                !
                if ((k1==m1.and.k2==m2.and.igamma==jgamma).and.abs(f_t-1.0_ark)>small_) then
                  write(out,"('Non normalized irreps for igamma,jgamma,k1,k2,m1,m2 = ',6i7,f16.7)") igamma,jgamma,k1,k2,m1,m2,f_t
                  stop 'Non normailized irreps'
                endif
                !          
              enddo          
            enddo          
          enddo          
        enddo          
      enddo          
    enddo          


   end subroutine check_characters_and_representation
   !
   recursive subroutine do_g36_transform(irep,ioper,dim,dst)

    integer(ik), intent(in)  ::  ioper,irep,dim
    real(ark), intent(out)   ::  dst(dim,dim)
    !
    integer(ik)  :: tn(72,2), temp(144)
    integer(ik) :: nsrc
    real(ark)   ::  mat1(dim,dim),mat2(dim,dim)
     !
     temp(1:36)   = (/0,0,2,0,6,4,0,7,2,3,2,3,4,5,6,4,5,6,0,21,19,2,3,2,3,2,3,4,5,6,4,5,6,4,5,6/)
     temp(73:108) = (/0,0,2,0,2,2,0,7,7,7,8,8,7,7,7,8,8,8,0,7,7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
     !
     !temp(1:36)   = (/0, 0, 2, 0, 4, 5, 0, 7, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 0,21,19, 2, 3, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6/)
     !temp(73:108) = (/0, 0, 2, 0, 2, 2, 0, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8, 0, 7, 7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
     ! 
     temp(37:72)  = (/0,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37/)
     temp(109:144)= (/0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36/)
     !
     tn = reshape( temp, (/ 72, 2/))
     !
        
     select case (ioper)
       !
     case (1,2,4,7,19,37) 
       !
       dst = sym%irr(irep,ioper)%repres
       !
     case default 
       !
       !
       call do_g36_transform(irep,tn(ioper,1),dim,mat1)
       call do_g36_transform(irep,tn(ioper,2),dim,mat2)
       !
       dst = matmul(mat2,mat1)
       !
     end select
     !
   end subroutine do_g36_transform
   !

end module symmetry




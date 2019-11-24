#define fwig  0

module rotme_cart_tens

use accuracy
use timer
use fwigxjpf
use moltype
implicit none


integer(ik), parameter :: maxnelem = 30


type kmat_type
  real(rk), allocatable :: me(:,:,:,:)
end type kmat_type


type rotme_cart_tens_type
  integer(ik)                   ::  nelem = 0
  integer(ik)                   ::  nirrep = 0
  character(cl)                 ::  name = 'N/A'
  character(cl), allocatable    ::  selem(:)
  character(cl), allocatable    ::  sirrep(:)
  complex(rk), allocatable      ::  tmat_s(:,:)
  complex(rk), allocatable      ::  tmat_x(:,:)
  logical                       ::  sym_a = .true.
  integer(ik)                   ::  jmin = 0
  integer(ik)                   ::  jmax = 0
  integer(ik), allocatable      ::  nktau(:)
  integer(ik), allocatable      ::  ktau(:,:,:)
  integer(ik), allocatable      ::  ktau_ind(:,:)
  integer(ik)                   ::  wang_ksign(2) = (/-1,1/)
  type(kmat_type), allocatable  ::  mmat(:,:)
  type(kmat_type), allocatable  ::  kmat(:,:)
  integer(ik)                   ::  dj = 1000
  integer(ik)                   ::  dm = 1000
  integer(ik)                   ::  dk = 1000
  integer(ik)                   ::  kmat_cmplx = 0
  integer(ik), allocatable      ::  mmat_cmplx(:)
  real(rk)                      ::  zero_tol = 1.0d-14
  procedure(prim_me_func), pointer, nopass :: func => null()
  contains
  procedure, public  ::  init => init_rotme_cart_tens_type
  procedure, public  ::  wang_coefs => wang_jktau_coefs
end type rotme_cart_tens_type


abstract interface
  subroutine prim_me_func(q1,q2,name,nelem,nirrep,mf,lf,sirrep,selem)
    use accuracy
    integer(ik), intent(in) :: q1(:), q2(:)
    character(cl), intent(out) :: name
    integer(ik), intent(inout) :: nirrep, nelem
    complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
    character(cl), intent(out), optional :: sirrep(:), selem(:)
  end subroutine prim_me_func
end interface


contains


!#include 'rotme_vzz.f90'
subroutine rotme_vzz_trace0(q1, q2, name, nelem, nirrep, mf, lf, sirrep, selem)

  integer(ik), intent(in) :: q1(:), q2(:)
  character(cl), intent(out) :: name
  integer(ik), intent(inout) :: nelem, nirrep
  complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
  character(cl), intent(out), optional :: sirrep(:), selem(:)

  integer(ik) :: rank(2), isigma, sigma, irrep, j1, j2, k1, k2, m1, m2, nelem_s
  complex(rk) :: tmat_s(6,6), tmat_x(6,6)

  j1 = q1(1)
  k1 = q1(2)
  m1 = q1(3)
  j2 = q2(1)
  k2 = q2(2)
  m2 = q2(3)

  name = 'Vzz'
  nelem = 6
  nirrep = 1
  nelem_s = 5
  rank(1:nirrep) = (/2/)

  if (present(sirrep)) then
    sirrep(1:nirrep) = (/'T(2)'/)
  endif

  if (present(selem)) then
    selem(1:nelem) = (/'xx','xy','xz','yy','yz','zz'/)
  endif


  ! sigma = -2
  tmat_s(1,1:nelem) = (/ cmplx( 0.5_rk/sqrt(6.0_rk), 0.0_rk ), &!    xx
                         cmplx( 0.0_rk, -1.0_rk/sqrt(6.0_rk) ), &!   xy
                         cmplx( 0.0_rk, 0.0_rk), &!                  xz
                         cmplx( -0.5_rk/sqrt(6.0_rk), 0.0_rk ), &!   yy
                         cmplx( 0.0_rk, 0.0_rk), &!                  yz
                         cmplx( 0.0_rk, 0.0_rk) /) !                 zz
  ! sigma = -1
  tmat_s(2,1:nelem) = (/ cmplx(0.0_rk, 0.0_rk), &!                   xx
                         cmplx(0.0_rk, 0.0_rk), &!                   xy
                         cmplx(1.0_rk/sqrt(6.0_rk), 0.0_rk), &!      xz
                         cmplx(0.0_rk, 0.0_rk), &!                   yy
                         cmplx(0.0_rk, -1.0_rk/sqrt(6.0_rk)), &!     yz
                         cmplx(0.0_rk, 0.0_rk) /) !                  zz

  ! sigma = 0
  tmat_s(3,1:nelem) = (/ cmplx(0.0_rk, 0.0_rk), &!                   xx
                         cmplx(0.0_rk, 0.0_rk), &!                   xy
                         cmplx(0.0_rk, 0.0_rk), &!                   xz
                         cmplx(0.0_rk, 0.0_rk), &!                   yy
                         cmplx(0.0_rk, 0.0_rk), &!                   yz
                         cmplx(0.5_rk, 0.0_rk) /) !                  zz

  ! sigma = +1
  tmat_s(4,1:nelem) = (/ cmplx(0.0_rk, 0.0_rk), &!                   xx
                         cmplx(0.0_rk, 0.0_rk), &!                   xy
                         cmplx(-1.0_rk/sqrt(6.0_rk), 0.0_rk), &!     xz
                         cmplx(0.0_rk, 0.0_rk), &!                   yy
                         cmplx(0.0_rk, -1.0_rk/sqrt(6.0_rk)), &!     yz
                         cmplx(0.0_rk, 0.0_rk) /) !                  zz

  ! sigma = +2
  tmat_s(5,1:nelem) = (/ cmplx( 0.5_rk/sqrt(6.0_rk), 0.0_rk ), &!    xx
                         cmplx( 0.0_rk, 1.0_rk/sqrt(6.0_rk) ), &!    xy
                         cmplx( 0.0_rk, 0.0_rk), &!                  xz
                         cmplx( -0.5_rk/sqrt(6.0_rk), 0.0_rk ), &!   yy
                         cmplx(0.0_rk, 0.0_rk), &!                   yz
                         cmplx(0.0_rk, 0.0_rk) /) !                  zz

  if (present(mf)) then
    mf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        mf(1:nelem,irrep) = mf(1:nelem,irrep) + threej_symbol(j1,rank(irrep),j2,k1,sigma,-k2) * (-1)**(k2) * tmat_s(isigma,1:nelem)
      enddo
    enddo
  endif

  if (present(lf)) then
    call pseudoinverse(nelem_s, nelem, tmat_s(1:nelem_s,1:nelem), tmat_x(1:nelem,1:nelem_s))
    lf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        lf(1:nelem,irrep) = lf(1:nelem,irrep) + tmat_x(1:nelem,isigma) * threej_symbol(j1,rank(irrep),j2,m1,sigma,-m2)
      enddo
    enddo
    lf = lf * (-1)**m2 * sqrt(real((2*j1+1)*(2*j2+1),rk))
  endif

end subroutine rotme_vzz_trace0


!#include 'rotme_alpha.f90'

subroutine rotme_alpha(q1, q2, name, nelem, nirrep, mf, lf, sirrep, selem)

  integer(ik), intent(in) :: q1(:), q2(:)
  character(cl), intent(out) :: name
  integer(ik), intent(inout) :: nelem, nirrep
  complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
  character(cl), intent(out), optional :: sirrep(:), selem(:)

  integer(ik) :: rank(2), isigma, sigma, irrep, j1, j2, k1, k2, m1, m2, nelem_s
  complex(rk) :: tmat_s(6,6), tmat_x(6,6)

  j1 = q1(1)
  k1 = q1(2)
  m1 = q1(3)
  j2 = q2(1)
  k2 = q2(2)
  m2 = q2(3)

  name = 'ALPHA'
  nelem = 6
  nirrep = 2
  nelem_s = 6
  rank(1:nirrep) = (/0,2/)

  if (present(sirrep)) then
    sirrep(1:nirrep) = (/'T(0)','T(2)'/)
  endif

  if (present(selem)) then
    selem(1:nelem) = (/'xx','xy','xz','yy','yz','zz'/)
  endif

  ! Polarizability
  !  0 {0: -sqrt(3)*xx/3 - sqrt(3)*yy/3 - sqrt(3)*zz/3}
  !  1 {0: 0, 1: 0, -1: 0}
  !  2 {0: -sqrt(6)*xx/6 - sqrt(6)*yy/6 + sqrt(6)*zz/3, 1: -xz - I*yz, 2: xx/2 + I*xy - yy/2, -1: xz - I*yz, -2: xx/2 - I*xy - yy/2}
  ! number of elements: 6
  ! order of Cartesian components: [xx, xy, xz, yy, yz, zz]
  ! order of spherical components: [[0, 0], [2, -2], [2, -1], [2, 0], [2, 1], [2, 2]]
  ! Cartesian to spherical transformation matrix:
  tmat_s( 1 , 1 )= -1.0_rk/3.0_rk*sqrt(3.0_rk)
  tmat_s( 1 , 2 )= 0
  tmat_s( 1 , 3 )= 0
  tmat_s( 1 , 4 )= -1.0_rk/3.0_rk*sqrt(3.0_rk)
  tmat_s( 1 , 5 )= 0
  tmat_s( 1 , 6 )= -1.0_rk/3.0_rk*sqrt(3.0_rk)
  tmat_s( 2 , 1 )= 1.0_rk/2.0_rk
  tmat_s( 2 , 2 )= cmplx(0,-1)
  tmat_s( 2 , 3 )= 0
  tmat_s( 2 , 4 )= -1.0_rk/2.0_rk
  tmat_s( 2 , 5 )= 0
  tmat_s( 2 , 6 )= 0
  tmat_s( 3 , 1 )= 0
  tmat_s( 3 , 2 )= 0
  tmat_s( 3 , 3 )= 1
  tmat_s( 3 , 4 )= 0
  tmat_s( 3 , 5 )= cmplx(0,-1)
  tmat_s( 3 , 6 )= 0
  tmat_s( 4 , 1 )= -1.0_rk/6.0_rk*sqrt(6.0_rk)
  tmat_s( 4 , 2 )= 0
  tmat_s( 4 , 3 )= 0
  tmat_s( 4 , 4 )= -1.0_rk/6.0_rk*sqrt(6.0_rk)
  tmat_s( 4 , 5 )= 0
  tmat_s( 4 , 6 )= (1.0_rk/3.0_rk)*sqrt(6.0_rk)
  tmat_s( 5 , 1 )= 0
  tmat_s( 5 , 2 )= 0
  tmat_s( 5 , 3 )= -1
  tmat_s( 5 , 4 )= 0
  tmat_s( 5 , 5 )= cmplx(0,-1)
  tmat_s( 5 , 6 )= 0
  tmat_s( 6 , 1 )= 1.0_rk/2.0_rk
  tmat_s( 6 , 2 )= cmplx(0,1)
  tmat_s( 6 , 3 )= 0
  tmat_s( 6 , 4 )= -1.0_rk/2.0_rk
  tmat_s( 6 , 5 )= 0
  tmat_s( 6 , 6 )= 0


  if (present(mf)) then
    mf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        mf(1:nelem,irrep) = mf(1:nelem,irrep) + threej_symbol(j1,rank(irrep),j2,k1,sigma,-k2) * (-1)**(k2) * tmat_s(isigma,1:nelem)
      enddo
    enddo
  endif

  if (present(lf)) then
    call pseudoinverse(nelem_s, nelem, tmat_s(1:nelem_s,1:nelem), tmat_x(1:nelem,1:nelem_s))
    lf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        lf(1:nelem,irrep) = lf(1:nelem,irrep) + tmat_x(1:nelem,isigma) * threej_symbol(j1,rank(irrep),j2,m1,sigma,-m2)
      enddo
    enddo
    lf = lf * (-1)**m2 * sqrt(real((2*j1+1)*(2*j2+1),rk))
  endif

end subroutine rotme_alpha



!#include 'rotme_beta.f90'

subroutine rotme_beta(q1, q2, name, nelem, nirrep, mf, lf, sirrep, selem)

  integer(ik), intent(in) :: q1(:), q2(:)
  character(cl), intent(out) :: name
  integer(ik), intent(inout) :: nelem, nirrep
  complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
  character(cl), intent(out), optional :: sirrep(:), selem(:)

  integer(ik) :: rank(2), isigma, sigma, irrep, j1, j2, k1, k2, m1, m2, nelem_s
  complex(rk) :: tmat_s(10,10), tmat_x(10,10)

  j1 = q1(1)
  k1 = q1(2)
  m1 = q1(3)
  j2 = q2(1)
  k2 = q2(2)
  m2 = q2(3)

  name = 'BETA'
  nelem = 10
  nirrep = 2
  nelem_s = 10
  rank(1:nirrep) = (/1, 3/)

  if (present(sirrep)) then
    sirrep(1:nirrep) = (/'T(1)','T(3)'/)
  endif

  if (present(selem)) then
    selem(1:nelem) = (/'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', 'yyy', 'yyz', 'yzz', 'zzz'/)
  endif

  ! Hyperpolarizability
  !  0 {0: 0}
  !  1 {0: -2*sqrt(15)*xxz/15 - 2*sqrt(15)*yyz/15 - 2*sqrt(15)*zzz/15, 1: sqrt(30)*xxx/15 + sqrt(30)*I*xxy/15 + sqrt(30)*xyy/15 + sqrt(30)*xzz/15 + sqrt(30)*I*yyy/15 + sqrt(30)*I*yzz/15, -1: -sqrt(30)*xxx/15 + sqrt(30)*I*xxy/15 - sqrt(30)*xyy/15 - sqrt(30)*xzz/15 + sqrt(30)*I*yyy/15 + sqrt(30)*I*yzz/15}
  !  2 {0: 0, 1: 0, 2: 0, -1: 0, -2: 0}
  !  3 {0: -3*sqrt(10)*xxz/10 - 3*sqrt(10)*yyz/10 + sqrt(10)*zzz/5, 1: sqrt(30)*xxx/20 + sqrt(30)*I*xxy/20 + sqrt(30)*xyy/20 - sqrt(30)*xzz/5 + sqrt(30)*I*yyy/20 - sqrt(30)*I*yzz/5, 2: sqrt(3)*xxz/2 + sqrt(3)*I*xyz - sqrt(3)*yyz/2, 3: -sqrt(2)*xxx/4 - 3*sqrt(2)*I*xxy/4 + 3*sqrt(2)*xyy/4 + sqrt(2)*I*yyy/4, -2: sqrt(3)*xxz/2 - sqrt(3)*I*xyz - sqrt(3)*yyz/2, -3: sqrt(2)*xxx/4 - 3*sqrt(2)*I*xxy/4 - 3*sqrt(2)*xyy/4 + sqrt(2)*I*yyy/4, -1: -sqrt(30)*xxx/20 + sqrt(30)*I*xxy/20 - sqrt(30)*xyy/20 + sqrt(30)*xzz/5 + sqrt(30)*I*yyy/20 - sqrt(30)*I*yzz/5}
  ! number of elements: 10
  ! order of Cartesian components: [xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz]
  ! order of spherical components: [[1, -1], [1, 0], [1, 1], [3, -3], [3, -2], [3, -1], [3, 0], [3, 1], [3, 2], [3, 3]]
  ! Cartesian to spherical transformation matrix:
  tmat_s( 1 , 1 )= -1.0_rk/15.0_rk*sqrt(30.0_rk)
  tmat_s( 1 , 2 )= cmplx(0,(1.0_rk/15.0_rk)*sqrt(30.0_rk))
  tmat_s( 1 , 3 )= 0
  tmat_s( 1 , 4 )= -1.0_rk/15.0_rk*sqrt(30.0_rk)
  tmat_s( 1 , 5 )= 0
  tmat_s( 1 , 6 )= -1.0_rk/15.0_rk*sqrt(30.0_rk)
  tmat_s( 1 , 7 )= cmplx(0,(1.0_rk/15.0_rk)*sqrt(30.0_rk))
  tmat_s( 1 , 8 )= 0
  tmat_s( 1 , 9 )= cmplx(0,(1.0_rk/15.0_rk)*sqrt(30.0_rk))
  tmat_s( 1 , 10 )= 0
  tmat_s( 2 , 1 )= 0
  tmat_s( 2 , 2 )= 0
  tmat_s( 2 , 3 )= -2.0_rk/15.0_rk*sqrt(15.0_rk)
  tmat_s( 2 , 4 )= 0
  tmat_s( 2 , 5 )= 0
  tmat_s( 2 , 6 )= 0
  tmat_s( 2 , 7 )= 0
  tmat_s( 2 , 8 )= -2.0_rk/15.0_rk*sqrt(15.0_rk)
  tmat_s( 2 , 9 )= 0
  tmat_s( 2 , 10 )= -2.0_rk/15.0_rk*sqrt(15.0_rk)
  tmat_s( 3 , 1 )= (1.0_rk/15.0_rk)*sqrt(30.0_rk)
  tmat_s( 3 , 2 )= cmplx(0,(1.0_rk/15.0_rk)*sqrt(30.0_rk))
  tmat_s( 3 , 3 )= 0
  tmat_s( 3 , 4 )= (1.0_rk/15.0_rk)*sqrt(30.0_rk)
  tmat_s( 3 , 5 )= 0
  tmat_s( 3 , 6 )= (1.0_rk/15.0_rk)*sqrt(30.0_rk)
  tmat_s( 3 , 7 )= cmplx(0,(1.0_rk/15.0_rk)*sqrt(30.0_rk))
  tmat_s( 3 , 8 )= 0
  tmat_s( 3 , 9 )= cmplx(0,(1.0_rk/15.0_rk)*sqrt(30.0_rk))
  tmat_s( 3 , 10 )= 0
  tmat_s( 4 , 1 )= (1.0_rk/4.0_rk)*sqrt(2.0_rk)
  tmat_s( 4 , 2 )= cmplx(0,-3.0_rk/4.0_rk*sqrt(2.0_rk))
  tmat_s( 4 , 3 )= 0
  tmat_s( 4 , 4 )= -3.0_rk/4.0_rk*sqrt(2.0_rk)
  tmat_s( 4 , 5 )= 0
  tmat_s( 4 , 6 )= 0
  tmat_s( 4 , 7 )= cmplx(0,(1.0_rk/4.0_rk)*sqrt(2.0_rk))
  tmat_s( 4 , 8 )= 0
  tmat_s( 4 , 9 )= 0
  tmat_s( 4 , 10 )= 0
  tmat_s( 5 , 1 )= 0
  tmat_s( 5 , 2 )= 0
  tmat_s( 5 , 3 )= (1.0_rk/2.0_rk)*sqrt(3.0_rk)
  tmat_s( 5 , 4 )= 0
  tmat_s( 5 , 5 )= cmplx(0,-sqrt(3.0_rk))
  tmat_s( 5 , 6 )= 0
  tmat_s( 5 , 7 )= 0
  tmat_s( 5 , 8 )= -1.0_rk/2.0_rk*sqrt(3.0_rk)
  tmat_s( 5 , 9 )= 0
  tmat_s( 5 , 10 )= 0
  tmat_s( 6 , 1 )= -1.0_rk/20.0_rk*sqrt(30.0_rk)
  tmat_s( 6 , 2 )= cmplx(0,(1.0_rk/20.0_rk)*sqrt(30.0_rk))
  tmat_s( 6 , 3 )= 0
  tmat_s( 6 , 4 )= -1.0_rk/20.0_rk*sqrt(30.0_rk)
  tmat_s( 6 , 5 )= 0
  tmat_s( 6 , 6 )= (1.0_rk/5.0_rk)*sqrt(30.0_rk)
  tmat_s( 6 , 7 )= cmplx(0,(1.0_rk/20.0_rk)*sqrt(30.0_rk))
  tmat_s( 6 , 8 )= 0
  tmat_s( 6 , 9 )= cmplx(0,-1.0_rk/5.0_rk*sqrt(30.0_rk))
  tmat_s( 6 , 10 )= 0
  tmat_s( 7 , 1 )= 0
  tmat_s( 7 , 2 )= 0
  tmat_s( 7 , 3 )= -3.0_rk/10.0_rk*sqrt(10.0_rk)
  tmat_s( 7 , 4 )= 0
  tmat_s( 7 , 5 )= 0
  tmat_s( 7 , 6 )= 0
  tmat_s( 7 , 7 )= 0
  tmat_s( 7 , 8 )= -3.0_rk/10.0_rk*sqrt(10.0_rk)
  tmat_s( 7 , 9 )= 0
  tmat_s( 7 , 10 )= (1.0_rk/5.0_rk)*sqrt(10.0_rk)
  tmat_s( 8 , 1 )= (1.0_rk/20.0_rk)*sqrt(30.0_rk)
  tmat_s( 8 , 2 )= cmplx(0,(1.0_rk/20.0_rk)*sqrt(30.0_rk))
  tmat_s( 8 , 3 )= 0
  tmat_s( 8 , 4 )= (1.0_rk/20.0_rk)*sqrt(30.0_rk)
  tmat_s( 8 , 5 )= 0
  tmat_s( 8 , 6 )= -1.0_rk/5.0_rk*sqrt(30.0_rk)
  tmat_s( 8 , 7 )= cmplx(0,(1.0_rk/20.0_rk)*sqrt(30.0_rk))
  tmat_s( 8 , 8 )= 0
  tmat_s( 8 , 9 )= cmplx(0,-1.0_rk/5.0_rk*sqrt(30.0_rk))
  tmat_s( 8 , 10 )= 0
  tmat_s( 9 , 1 )= 0
  tmat_s( 9 , 2 )= 0
  tmat_s( 9 , 3 )= (1.0_rk/2.0_rk)*sqrt(3.0_rk)
  tmat_s( 9 , 4 )= 0
  tmat_s( 9 , 5 )= cmplx(0,sqrt(3.0_rk))
  tmat_s( 9 , 6 )= 0
  tmat_s( 9 , 7 )= 0
  tmat_s( 9 , 8 )= -1.0_rk/2.0_rk*sqrt(3.0_rk)
  tmat_s( 9 , 9 )= 0
  tmat_s( 9 , 10 )= 0
  tmat_s( 10 , 1 )= -1.0_rk/4.0_rk*sqrt(2.0_rk)
  tmat_s( 10 , 2 )= cmplx(0,-3.0_rk/4.0_rk*sqrt(2.0_rk))
  tmat_s( 10 , 3 )= 0
  tmat_s( 10 , 4 )= (3.0_rk/4.0_rk)*sqrt(2.0_rk)
  tmat_s( 10 , 5 )= 0
  tmat_s( 10 , 6 )= 0
  tmat_s( 10 , 7 )= cmplx(0,(1.0_rk/4.0_rk)*sqrt(2.0_rk))
  tmat_s( 10 , 8 )= 0
  tmat_s( 10 , 9 )= 0
  tmat_s( 10 , 10 )= 0

  if (present(mf)) then
    mf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        mf(1:nelem,irrep) = mf(1:nelem,irrep) + threej_symbol(j1,rank(irrep),j2,k1,sigma,-k2) * (-1)**(k2) * tmat_s(isigma,1:nelem)
      enddo
    enddo
  endif

  if (present(lf)) then
    call pseudoinverse(nelem_s, nelem, tmat_s(1:nelem_s,1:nelem), tmat_x(1:nelem,1:nelem_s))
    lf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        lf(1:nelem,irrep) = lf(1:nelem,irrep) + tmat_x(1:nelem,isigma) * threej_symbol(j1,rank(irrep),j2,m1,sigma,-m2)
      enddo
    enddo
    lf = lf * (-1)**m2 * sqrt(real((2*j1+1)*(2*j2+1),rk))
  endif

end subroutine rotme_beta


!#include 'rotme_mu.f90'

subroutine rotme_mu(q1, q2, name, nelem, nirrep, mf, lf, sirrep, selem)

  integer(ik), intent(in) :: q1(:), q2(:)
  character(cl), intent(out) :: name
  integer(ik), intent(inout) :: nelem, nirrep
  complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
  character(cl), intent(out), optional :: sirrep(:), selem(:)

  integer(ik) :: rank(1), isigma, sigma, irrep, j1, j2, k1, k2, m1, m2, nelem_s
  complex(rk) :: tmat_s(3,3), tmat_x(3,3)

  j1 = q1(1)
  k1 = q1(2)
  m1 = q1(3)
  j2 = q2(1)
  k2 = q2(2)
  m2 = q2(3)

  name = 'MU'
  nelem = 3
  nirrep = 1
  nelem_s = 3
  rank(1:nirrep) = (/1/)

  if (present(sirrep)) then
    sirrep(1:nirrep) = (/'T(1)'/)
  endif

  if (present(selem)) then
    selem(1:nelem) = (/'x','y','z'/)
  endif

  ! Dipole moment
  !  0 {0: 0}
  !  1 {0: z, 1: sqrt(2)*(-x - I*y)/2, -1: sqrt(2)*(x - I*y)/2}
  ! number of elements: 3
  ! order of Cartesian components: [x, y, z]
  ! order of spherical components: [[1, -1], [1, 0], [1, 1]]
  ! Cartesian to spherical transformation matrix:
  tmat_s( 1 , 1 )= (1.0_rk/2.0_rk)*sqrt(2.0_rk)
  tmat_s( 1 , 2 )= cmplx(0,-1.0_rk/2.0_rk*sqrt(2.0_rk))
  tmat_s( 1 , 3 )= 0
  tmat_s( 2 , 1 )= 0
  tmat_s( 2 , 2 )= 0
  tmat_s( 2 , 3 )= 1
  tmat_s( 3 , 1 )= -1.0_rk/2.0_rk*sqrt(2.0_rk)
  tmat_s( 3 , 2 )= cmplx(0,-1.0_rk/2.0_rk*sqrt(2.0_rk))
  tmat_s( 3 , 3 )= 0


  if (present(mf)) then
    mf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        mf(1:nelem,irrep) = mf(1:nelem,irrep) + threej_symbol(j1,rank(irrep),j2,k1,sigma,-k2) * (-1)**(k2) * tmat_s(isigma,1:nelem)
      enddo
    enddo
  endif

  if (present(lf)) then
    call pseudoinverse(nelem_s, nelem, tmat_s(1:nelem_s,1:nelem), tmat_x(1:nelem,1:nelem_s))
    lf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        lf(1:nelem,irrep) = lf(1:nelem,irrep) + tmat_x(1:nelem,isigma) * threej_symbol(j1,rank(irrep),j2,m1,sigma,-m2)
      enddo
    enddo
    lf = lf * (-1)**m2 * sqrt(real((2*j1+1)*(2*j2+1),rk))
  endif

end subroutine rotme_mu


!#include 'rotme_costheta.f90'

! <j',k',m'|cos(theta)|j,k,m> = <j',k',m'|d_{0,0}^{1}|j,k,m>

subroutine rotme_costheta(q1, q2, name, nelem, nirrep, mf, lf, sirrep, selem)

  integer(ik), intent(in) :: q1(:), q2(:)
  character(cl), intent(out) :: name
  integer(ik), intent(inout) :: nelem, nirrep
  complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
  character(cl), intent(out), optional :: sirrep(:), selem(:)

  integer(ik) :: j1, j2, k1, k2, m1, m2

  j1 = q1(1)
  k1 = q1(2)
  m1 = q1(3)
  j2 = q2(1)
  k2 = q2(2)
  m2 = q2(3)

  name = 'COSTHETA'
  nelem = 1
  nirrep = 1

  if (present(sirrep)) then
    sirrep(1:nirrep) = (/'T(0)'/)
  endif

  if (present(selem)) then
    selem(1:nelem) = (/''/)
  endif

  if (present(mf)) then
    if (k1==k2) then
      mf = threej_symbol(j1,1,j2,k1,0,-k2) * (-1)**(k2)
    else
      mf = 0
    endif
  endif

  if (present(lf)) then
    if (m1==m2) then
      lf = threej_symbol(j1,1,j2,m1,0,-m2) * (-1)**m2 * sqrt(real((2*j1+1)*(2*j2+1),rk))
    else
      lf = 0
    endif
  endif

end subroutine rotme_costheta



!#include 'rotme_3cos2theta_min1.f90'

! <j',k',m'|3cos^2(theta)-1|j,k,m> = 2 <j',k',m'|d_{0,0}^{2}|j,k,m>

subroutine rotme_3cos2theta_min1(q1, q2, name, nelem, nirrep, mf, lf, sirrep, selem)

  integer(ik), intent(in) :: q1(:), q2(:)
  character(cl), intent(out) :: name
  integer(ik), intent(inout) :: nelem, nirrep
  complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
  character(cl), intent(out), optional :: sirrep(:), selem(:)

  integer(ik) :: j1, j2, k1, k2, m1, m2

  j1 = q1(1)
  k1 = q1(2)
  m1 = q1(3)
  j2 = q2(1)
  k2 = q2(2)
  m2 = q2(3)

  name = '3COSTHETA2_MIN1'
  nelem = 1
  nirrep = 1

  if (present(sirrep)) then
    sirrep(1:nirrep) = (/'T(0)'/)
  endif

  if (present(selem)) then
    selem(1:nelem) = (/''/)
  endif

  if (present(mf)) then
    if (k1==k2) then
      mf = threej_symbol(j1,2,j2,k1,0,-k2) * (-1)**(k2)
    else
      mf = 0
    endif
  endif

  if (present(lf)) then
    if (m1==m2) then
      lf(1,1) = threej_symbol(j1,2,j2,m1,0,-m2) * (-1)**m2 * sqrt(real((2*j1+1)*(2*j2+1),rk))
    else
      lf = 0
    endif
  endif

end subroutine rotme_3cos2theta_min1


!#include 'rotme_mf.f90'

subroutine rotme_mf(q1, q2, name, nelem, nirrep, mf, lf, sirrep, selem)

  integer(ik), intent(in) :: q1(:), q2(:)
  character(cl), intent(out) :: name
  integer(ik), intent(inout) :: nelem, nirrep
  complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
  character(cl), intent(out), optional :: sirrep(:), selem(:)

  integer(ik) :: j1, j2, k1, k2, m1, m2

  j1 = q1(1)
  k1 = q1(2)
  m1 = q1(3)
  j2 = q2(1)
  k2 = q2(2)
  m2 = q2(3)

  name = ''
  nirrep = 1
  nelem = 1

  if (present(sirrep)) then
    sirrep(1:nirrep) = 'T(0)'
  endif

  if (present(selem)) then
    selem(1:nelem) = 'MF'
  endif

  if (present(mf)) then
    mf(:,:) = 0.0
    if (j1==j2 .and. k1==k2) then
      mf = 1.0
    endif
  endif

  if (present(lf)) then
    lf(:,:) = 0.0
    if (j1==j2 .and. m1==m2) then
      lf = 1.0
    endif
  endif

end subroutine rotme_mf


!#include 'rotme_j.f90'

subroutine rotme_j(q1, q2, name, nelem, nirrep, mf, lf, sirrep, selem)

  integer(ik), intent(in) :: q1(:), q2(:)
  character(cl), intent(out) :: name
  integer(ik), intent(inout) :: nelem, nirrep
  complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
  character(cl), intent(out), optional :: sirrep(:), selem(:)

  integer(ik) :: rank(1), isigma, sigma, irrep, j1, j2, k1, k2, m1, m2, nelem_s, ncoef(3), kval(2,3), fac, ialpha, icoef
  real(ark) :: jm_coef, jp_coef
  complex(rk) :: tmat_s(3,3), tmat_x(3,3), coef(2,3), cx1, cx2, cy1, cy2, cz1, res(3)

  j1 = q1(1)
  k1 = q1(2)
  m1 = q1(3)
  j2 = q2(1)
  k2 = q2(2)
  m2 = q2(3)

  name = 'J'
  nelem = 3
  nirrep = 1
  nelem_s = 3
  rank(1:nirrep) = (/1/)

  if (present(sirrep)) then
    sirrep(1:nirrep) = (/'T(1)'/)
  endif

  if (present(selem)) then
    selem(1:nelem) = (/'x','y','z'/)
  endif


  fac = j1*(j1+1)-k1*(k1-1)
  if (fac<0) fac = 0
  jp_coef = sqrt(real(fac,ark))

  fac = j1*(j1+1)-k1*(k1+1)
  if (fac<0) fac = 0
  jm_coef = sqrt(real(fac,ark))

  cx1 =  jp_coef*0.5_ark
  cx2 =  jm_coef*0.5_ark
  cy1 =  jm_coef*cmplx(0.0_ark,0.5_ark)
  cy2 = -jp_coef*cmplx(0.0_ark,0.5_ark)
  cz1 = k1

  ncoef(1:3) = (/2,2,1/)
  coef(1:2,1) = (/cx1,cx2/);  kval(1:2,1) = (/1,-1/)
  coef(1:2,2) = (/cy1,cy2/);  kval(1:2,2) = (/-1,1/)
  coef(1:1,3) = (/cz1/);      kval(1:1,3) = (/0/)


  !  0 {0: 0}
  !  1 {0: z, 1: sqrt(2)*(-x - I*y)/2, -1: sqrt(2)*(x - I*y)/2}
  ! number of elements: 3
  ! order of Cartesian components: [x, y, z]
  ! order of spherical components: [[1, -1], [1, 0], [1, 1]]
  ! Cartesian to spherical transformation matrix:
  tmat_s( 1 , 1 )= (1.0_rk/2.0_rk)*sqrt(2.0_rk)
  tmat_s( 1 , 2 )= cmplx(0,-1.0_rk/2.0_rk*sqrt(2.0_rk))
  tmat_s( 1 , 3 )= 0
  tmat_s( 2 , 1 )= 0
  tmat_s( 2 , 2 )= 0
  tmat_s( 2 , 3 )= 1
  tmat_s( 3 , 1 )= -1.0_rk/2.0_rk*sqrt(2.0_rk)
  tmat_s( 3 , 2 )= cmplx(0,-1.0_rk/2.0_rk*sqrt(2.0_rk))
  tmat_s( 3 , 3 )= 0


  if (present(mf)) then
    mf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1

        res = 0
        do ialpha=1, nelem
          do icoef=1, ncoef(ialpha)
            res(ialpha) = res(ialpha) + coef(icoef,ialpha) * threej_symbol(j1,rank(irrep),j2,k1-kval(icoef,ialpha),sigma,-k2)
          enddo
        enddo

        mf(1:nelem,irrep) = mf(1:nelem,irrep) + res(1:nelem) * (-1)**(k2) * tmat_s(isigma,1:nelem)
      enddo
    enddo
  endif

  if (present(lf)) then
    call pseudoinverse(nelem_s, nelem, tmat_s(1:nelem_s,1:nelem), tmat_x(1:nelem,1:nelem_s))
    lf(:,:) = 0.0
    isigma = 0
    do irrep=1, nirrep
      do sigma=-rank(irrep), rank(irrep)
        isigma = isigma + 1
        lf(1:nelem,irrep) = lf(1:nelem,irrep) + tmat_x(1:nelem,isigma) * threej_symbol(j1,rank(irrep),j2,m1,sigma,-m2)
      enddo
    enddo
    lf = lf * (-1)**m2 * sqrt(real((2*j1+1)*(2*j2+1),rk))
  endif

end subroutine rotme_j


!###################################################################################################################################



subroutine init_rotme_cart_tens_type(tens, jmin, jmax, dj, verbose)

  class(rotme_cart_tens_type) :: tens
  integer(ik), intent(in) :: jmin, jmax, dj
  logical, optional, intent(in) :: verbose

  integer(ik) :: info, j, ktau, tau_range(2), tau, nimg_j, j1, j2, nktau1, nktau2, nimg_k, ktau1, ktau2, ncoefs1, ncoefs2, icoef, jcoef, nimg_elem, &!
                 irrep, ielem, dk_max, dj_max, k, k1, k2, tau1, tau2, nelem, nirrep, dm_max, m1, m2, m1min, m1max, m2min, m2max, mmin, mmax, &!
                 nimg_elem_lf(maxnelem), nimg_k_lf(maxnelem), nimg_j_lf(maxnelem)
  real(rk) :: zero_tol
  complex(rk) :: coefs1(2), coefs2(2), me(maxnelem,maxnelem), dens, res(maxnelem,maxnelem)
  logical :: verbose_, img_j((jmax+1)**2), img_k((2*jmax+1)**2), img_elem(maxnelem**2), img_elem_lf(maxnelem**2,maxnelem), &!
             img_k_lf((2*jmax+1)**2,maxnelem), img_j_lf((jmax+1)**2,maxnelem)
  character(cl) :: name

  write(out, '(///a)') 'Compute rotational matrix elements of Cartesian tensor operator (rotme_cart_tens/init_rotme_cart_tens_type)'

  verbose_ = .false.
  if (present(verbose)) verbose_ = verbose

  zero_tol = tens%zero_tol
  tens%jmin = jmin
  tens%jmax = jmax

  write(out, '(/1x,a,2(1x,i4)/1x,a,1x,i4)') 'min and max J quanta:', jmin, jmax, "user-defined |j-j'| selection rules:", dj

  call tens%func( (/0,0,0/), (/0,0,0/), name,nelem,nirrep)
  tens%name = name
  tens%nelem = nelem
  tens%nirrep = nirrep

  write(out, '(/1x,a,1x,a/1x,a,1x,i4/1x,a,1x,i4)') 'tensor:', tens%name, 'no. Cartesian components:', tens%nelem, 'no. irreps:', tens%nirrep

  if (allocated(tens%sirrep)) deallocate(tens%sirrep)
  if (allocated(tens%selem)) deallocate(tens%selem)
  allocate(tens%sirrep(tens%nirrep),tens%selem(tens%nelem), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%sirrep(tens%nirrep),tens%selem(tens%nelem)', &!
    'tens%nelem, tens%nirrep =', tens%nelem, tens%nirrep
    stop 'STOP'
  endif

  call tens%func( (/0,0,0/), (/0,0,0/), name,nelem,nirrep, sirrep=tens%sirrep(1:tens%nirrep), selem=tens%selem(1:tens%nelem))

  write(out, '(1x,a,100(3x,a))') 'Cartesian components:', (trim(tens%selem(ielem)), ielem=1, tens%nelem)
  write(out, '(1x,a,100(3x,a))') 'irreps:', (trim(tens%sirrep(irrep)), irrep=1, tens%nirrep)


  ! initialize rotational quanta K and tau

  if (allocated(tens%nktau)) deallocate(tens%nktau)
  if (allocated(tens%ktau)) deallocate(tens%ktau)
  if (allocated(tens%ktau_ind)) deallocate(tens%ktau_ind)
  allocate( tens%nktau(jmin:jmax), tens%ktau(2,2*jmax+1,jmin:jmax), tens%ktau_ind(0:jmax,0:1), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%nktau(jmin:jmax), tens%ktau(2,2*jmax+1,jmin:jmax), , tens%ktau_ind(0:jmax,0:1)', &!
    'jmin, jmax =', jmin, jmax
    stop 'STOP'
  endif
  tens%nktau = 0
  tens%ktau_ind = 0

  do j=jmin, jmax
    ktau = 0
    do k=0, j
      tau_range = (/0,1/)
      if (k==0) tau_range = mod(j,2)
      do tau=tau_range(1), tau_range(2)
        ktau = ktau + 1
        tens%ktau(1:2,ktau,j) = (/k,tau/)
        tens%ktau_ind(k,tau) = ktau
      enddo
    enddo
    tens%nktau(j) = ktau
  enddo
  tens%wang_ksign = (/1,-1/)


  ! init module "fwigxjpf" for computing 3j symbols

  write(out, '(/1x,a,1x,i3)') 'initialize "fwigxjpf" for computing 3-j symbols with max angular momentum =', jmax
#if (fwig>0)
  call fwig_table_init(int(2*(jmax+1),kind=4), int(3,kind=4))
  call fwig_temp_init(int(2*(jmax+1),kind=4))
#endif

  dj_max = 0
  dk_max = 0
  dm_max = 0


  ! compute rotaitonal matrix elements, molecule-fixed K-tensor

  if (allocated(tens%kmat)) deallocate(tens%kmat)
  allocate(tens%kmat(jmin:jmax,jmin:jmax), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%kmat(jmin:jmax,jmin:jmax)', &!
    'jmin, jmax =', jmin, jmax
    stop 'STOP'
  endif

  if (verbose_) then
    write(out, '(/1x,a/2x,a,2x,a,4x,a,1x,a,4x,a,1x,a,9x,a,8x,a)') &!
    'molecule-fixed frame matrix elements', 'j1', 'j2', 'k1', 'tau1', 'k2', 'tau2', 'elem', 'irrep'
  endif

  nimg_j = 0
  img_j(:) = .false.

  do j1=jmin, jmax
    do j2=jmin, jmax

      if (abs(j1-j2)>dj) cycle

      nktau1 = tens%nktau(j1)
      nktau2 = tens%nktau(j2)

      if (allocated(tens%kmat(j1,j2)%me)) deallocate(tens%kmat(j1,j2)%me)
      allocate(tens%kmat(j1,j2)%me(tens%nelem,tens%nirrep,nktau1,nktau2), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%kmat(j1,j2)%me(tens%nelem,tens%nirrep,nktau1,nktau2)', &!
        'j1, j2, tens%nelem, tens%nirrep, nktau1, nktau2 =', j1, j2, tens%nelem, tens%nirrep, nktau1, nktau2
        stop 'STOP'
      endif

      tens%kmat(j1,j2)%me = 0

      nimg_k = 0
      img_k(:) = .false.

      do ktau1=1, nktau1

        k1 = tens%ktau(1,ktau1,j1)
        tau1 = tens%ktau(2,ktau1,j1)
        call tens%wang_coefs(j1, k1, tau1, ncoefs1, coefs1)

        do ktau2=1, nktau2

          k2 = tens%ktau(1,ktau2,j2)
          tau2 = tens%ktau(2,ktau2,j2)
          call tens%wang_coefs(j2, k2, tau2, ncoefs2, coefs2)

          me(:,:) = 0
          do icoef=1, ncoefs1
            do jcoef=1, ncoefs2
              dens = coefs1(icoef) * conjg(coefs2(jcoef))
              call tens%func( (/j1,k1*tens%wang_ksign(icoef),0/), (/j2,k2*tens%wang_ksign(jcoef),0/), name,nelem,nirrep, mf=res )
              me(1:nelem,1:nirrep) = me(1:nelem,1:nirrep) + res(1:nelem,1:nirrep) * dens
            enddo
          enddo

          nimg_elem = 0
          img_elem(:) = .false.

          do irrep=1, tens%nirrep
            do ielem=1, tens%nelem

              if ( abs(real(me(ielem,irrep),rk))>zero_tol .and. abs(aimag(me(ielem,irrep)))<=zero_tol ) then
                nimg_elem = nimg_elem + 1
                img_elem(nimg_elem) = .false.
                tens%kmat(j1,j2)%me(ielem,irrep,ktau1,ktau2) = real(me(ielem,irrep),rk)
              elseif ( abs(real(me(ielem,irrep),rk))<=zero_tol .and. abs(aimag(me(ielem,irrep)))>zero_tol ) then
                nimg_elem = nimg_elem + 1
                img_elem(nimg_elem) = .true.
                tens%kmat(j1,j2)%me(ielem,irrep,ktau1,ktau2) = aimag(me(ielem,irrep))
              elseif ( abs(real(me(ielem,irrep),rk))<=zero_tol .and. abs(aimag(me(ielem,irrep)))<=zero_tol ) then
                cycle
              else
                write(out, '(/a)') 'rotme_cart_tens/init_rotme_cart_tens_type error: element of K-tensor is a complex-valued number'
                write(out, '(1x,a,2(1x,i3),1x,a,2(1x,es16.8),1x,a)') 'ielem/irrep/value:', ielem,irrep,'(', me(ielem,irrep), 'i)'
                stop 'STOP'
              endif

              if (verbose_) then
                if (img_elem(nimg_elem)) then
                  write(out, '(1x,i3,1x,i3,3x,i3,2x,i3,3x,i3,2x,i3,3x,a10,3x,a10,3x,a,1x,es16.8,1x,es16.8,1x,a)') &!
                  j1, j2, k1, tau1, k2, tau2, trim(tens%selem(ielem)), trim(tens%sirrep(irrep)), &!
                  '(', 0.0, tens%kmat(j1,j2)%me(ielem,irrep,ktau1,ktau2), ' i )'
                else
                  write(out, '(1x,i3,1x,i3,3x,i3,2x,i3,3x,i3,2x,i3,3x,a10,3x,a10,3x,a,1x,es16.8,1x,es16.8,1x,a)') &!
                  j1, j2, k1, tau1, k2, tau2, trim(tens%selem(ielem)), trim(tens%sirrep(irrep)), &!
                  '(', tens%kmat(j1,j2)%me(ielem,irrep,ktau1,ktau2), 0.0, ' i )'
                endif
              endif

            enddo
          enddo

          if (nimg_elem>0) then
            nimg_k = nimg_k + 1
            if (all(img_elem(1:nimg_elem))) then
              img_k(nimg_k) = .true.
            elseif (all(.not.img_elem(1:nimg_elem))) then
              img_k(nimg_k) = .false.
            else
              write(out, '(/a/a,6(1x,i3))') 'rotme_cart_tens/init_rotme_cart_tens_type error: K-tensor has mixed real- and imaginary-valued elements', &!
              'j1,k1,tau1/j2,k2,tau2 =', j1,k1,tau1, j2,k2,tau2
              stop 'STOP'
            endif
          endif

          if (any(abs(tens%kmat(j1,j2)%me(:,:,ktau1,ktau2))>zero_tol)) then
            dk_max = max(dk_max,abs(ktau1-ktau2))
            dj_max = max(dj_max,abs(j1-j2))
          endif

        enddo ! ktau2
      enddo ! ktau1

      if (nimg_k>0) then
        nimg_j = nimg_j + 1
        if (all(img_k(1:nimg_k))) then
          img_j(nimg_j) = .true.
        elseif (all(.not.img_k(1:nimg_k))) then
          img_j(nimg_j) = .false.
        else
          write(out, '(/a/a,2(1x,i3))') &!
          'rotme_cart_tens/init_rotme_cart_tens_type error: K-tensor has mixed real- and imaginary-valued elements for different pairs of k,tau-quanta', &!
          ' j1/j2 =', j1, j2
          stop 'STOP'
        endif
      endif

    enddo ! j2
  enddo ! j1

  if (all(img_j(1:nimg_j))) then
    tens%kmat_cmplx = -1 ! imaginary numbers
  elseif (all(.not.img_j(1:nimg_j))) then
    tens%kmat_cmplx = 0 ! real numbers
  else
    write(out, '(/a)') &!
    'rotme_cart_tens/init_rotme_cart_tens_type error: K-tensor has mixed real- and imaginary-valued elements for different pairs of J-quanta'
    stop 'STOP'
  endif


  ! compute rotaitonal matrix elements, laboratory-fixed M-tensor

  mmin = -jmax
  mmax = jmax

  if (allocated(tens%mmat)) deallocate(tens%mmat)
  if (allocated(tens%mmat_cmplx)) deallocate(tens%mmat_cmplx)
  allocate(tens%mmat(jmin:jmax,jmin:jmax), tens%mmat_cmplx(tens%nelem), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%mmat(jmin:jmax,jmin:jmax), tens%mmat_cmplx(tens%nelem)', &!
    'jmin, jmax, tens%nelem =', jmin, jmax, tens%nelem
    stop 'STOP'
  endif

  if (verbose_) then
    write(out, '(/1x,a/2x,a,2x,a,4x,a,2x,a,8x,a,9x,a)') &!
    'laboratory-fixed frame matrix elements', 'j1', 'j2', 'm1', 'm2', 'irrep', 'elem'
  endif

  nimg_j_lf(:) = 0
  img_j_lf(:,:) = .false.

  do j1=jmin, jmax
    do j2=jmin, jmax

      if (abs(j1-j2)>dj) cycle

      m1min = max(mmin,-j1)
      m1max = min(mmax,j1)
      m2min = max(mmin,-j2)
      m2max = min(mmax,j2)

      if (allocated(tens%mmat(j1,j2)%me)) deallocate(tens%mmat(j1,j2)%me)
      allocate(tens%mmat(j1,j2)%me(tens%nirrep,tens%nelem,m1min:m1max,m2min:m2max), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%mmat(j1,j2)%me(tens%nirrep,tens%nelem,m1min:m1max,m2min:m2max)', &!
        'j1, j2, tens%nirrep, tens%nelem, m1min, m1max, m2min, m2max =', j1, j2, tens%nirrep, tens%nelem, m1min, m1max, m2min, m2max
        stop 'STOP'
      endif

      tens%mmat(j1,j2)%me = 0

      nimg_k_lf(:) = 0
      img_k_lf(:,:) = .false.

      do m1=m1min, m1max
        do m2=m2min, m2max

          call tens%func( (/j1,0,m1/), (/j2,0,m2/), name,nelem,nirrep, lf=res)
          me(1:nelem,1:nirrep) = res(1:nelem,1:nirrep)


          nimg_elem_lf(:) = 0
          img_elem_lf(:,:) = .false.

          do ielem=1, tens%nelem
            do irrep=1, tens%nirrep

              if ( abs(real(me(ielem,irrep),rk))>zero_tol .and. abs(aimag(me(ielem,irrep)))<=zero_tol ) then
                nimg_elem_lf(ielem) = nimg_elem_lf(ielem) + 1
                img_elem_lf(nimg_elem_lf(ielem),ielem) = .false.
                tens%mmat(j1,j2)%me(irrep,ielem,m1,m2) = real(me(ielem,irrep),rk)
              elseif ( abs(real(me(ielem,irrep),rk))<=zero_tol .and. abs(aimag(me(ielem,irrep)))>zero_tol ) then
                nimg_elem_lf(ielem) = nimg_elem_lf(ielem) + 1
                img_elem_lf(nimg_elem_lf(ielem),ielem) = .true.
                tens%mmat(j1,j2)%me(irrep,ielem,m1,m2) = aimag(me(ielem,irrep))
              elseif ( abs(real(me(ielem,irrep),rk))<=zero_tol .and. abs(aimag(me(ielem,irrep)))<=zero_tol ) then
                cycle
              else
                write(out, '(/a)') 'rotme_cart_tens/init_rotme_cart_tens_type error: element of M-tensor is a complex-valued number'
                write(out, '(1x,a,2(1x,i3),1x,a,2(1x,es16.8),1x,a)') 'ielem/irrep/value:', ielem,irrep,'(', me(ielem,irrep), 'i)'
                stop 'STOP'
              endif

              if (verbose_) then
                if (img_elem_lf(nimg_elem_lf(ielem),ielem)) then
                  write(out, '(1x,i3,1x,i3,3x,i3,1x,i3,3x,a10,3x,a10,3x,a,1x,es16.8,1x,es16.8,1x,a)') &!
                  j1, j2, m1, m2, trim(tens%sirrep(irrep)), trim(tens%selem(ielem)), &!
                  '(', 0.0, tens%mmat(j1,j2)%me(irrep,ielem,m1,m2), ' i )'
                else
                  write(out, '(1x,i3,1x,i3,3x,i3,1x,i3,3x,a10,3x,a10,3x,a,1x,es16.8,1x,es16.8,1x,a)') &!
                  j1, j2, m1, m2, trim(tens%sirrep(irrep)), trim(tens%selem(ielem)), &!
                  '(', tens%mmat(j1,j2)%me(irrep,ielem,m1,m2), 0.0, ' i )'
                endif
              endif

            enddo ! irrep

            if (nimg_elem_lf(ielem)>0) then
              nimg_k_lf(ielem) = nimg_k_lf(ielem) + 1
              if (all(img_elem_lf(1:nimg_elem_lf(ielem),ielem))) then
                img_k_lf(nimg_k_lf(ielem),ielem) = .true.
              elseif (all(.not.img_elem_lf(1:nimg_elem_lf(ielem),ielem))) then
                img_k_lf(nimg_k_lf(ielem),ielem) = .false.
              else
                write(out, '(/a/a,4(1x,i3))') 'rotme_cart_tens/init_rotme_cart_tens_type error: M-tensor has mixed real- and imaginary-valued elements', &!
                'j1,m1/j2,m2 =', j1,m1, j2,m2
                stop 'STOP'
              endif
            endif

          enddo ! ielem

          if (any(abs(tens%mmat(j1,j2)%me(:,:,m1,m2))>zero_tol)) then
            dm_max = max(dm_max,abs(m1-m2))
            dj_max = max(dj_max,abs(j1-j2))
          endif

        enddo ! m2
      enddo ! m1

      do ielem=1, nelem
        if (nimg_k_lf(ielem)>0) then
          nimg_j_lf(ielem) = nimg_j_lf(ielem) + 1
          if (all(img_k_lf(1:nimg_k_lf(ielem),ielem))) then
            img_j_lf(nimg_j_lf(ielem),ielem) = .true.
          elseif (all(.not.img_k_lf(1:nimg_k_lf(ielem),ielem))) then
            img_j_lf(nimg_j_lf(ielem),ielem) = .false.
          else
            write(out, '(/a/a,2(1x,i3))') &!
            'rotme_cart_tens/init_rotme_cart_tens_type error: M-tensor has mixed real- and imaginary-valued elements for different pairs of m-quanta', &!
            'j1/j2 =', j1, j2
            stop 'STOP'
          endif
        endif
      enddo

    enddo ! j2
  enddo ! j1

  do ielem=1, nelem
    if (all(img_j_lf(1:nimg_j_lf(ielem),ielem))) then
      tens%mmat_cmplx(ielem) = -1
    elseif (all(.not.img_j_lf(1:nimg_j_lf(ielem),ielem))) then
      tens%mmat_cmplx(ielem) = 0
    else
      write(out, '(/a)') &!
      'rotme_cart_tens/init_rotme_cart_tens_type error: M-tensor has mixed real- and imaginary-valued elements for different pairs of J-quanta'
      stop 'STOP'
    endif
  enddo


  tens%dj = dj_max
  tens%dm = dm_max
  tens%dk = dk_max

  write(out, '(/1x,a,1x,i4)') "|ktau-ktau'| selection rules:", tens%dk
  write(out, '(1x,a,1x,i4)') "|m-m'| selection rules:", tens%dm
  write(out, '(1x,a,1x,i4)') "|j-j'| selection rules:", tens%dj

  write(out,'(a)') 'done (rotme_cart_tens/init_rotme_cart_tens_type)'

end subroutine init_rotme_cart_tens_type



!###################################################################################################################################



subroutine wang_jktau_coefs(tens, j, k, tau, ncoefs, coefs)

  class(rotme_cart_tens_type), intent(in) :: tens
  integer(ik), intent(in) :: j, k, tau
  integer(ik), intent(out) :: ncoefs
  complex(rk), intent(out), optional :: coefs(2)

  integer(ik) :: sigma
  real(rk) :: fac1, fac2

  sigma = mod(k,3)*tau
  fac1 = real((-1)**sigma,rk)/sqrt(2.0_rk)
  fac2 = fac1 * (-1)**(j+k)

  if (tau==0) then

    if (k==0) then
      ncoefs = 1
      if (present(coefs)) coefs = (/ cmplx(1.0_rk,0.0_rk,kind=rk), cmplx(0.0_rk,0.0_rk,kind=rk) /)
    elseif (k>0) then
      ncoefs = 2
      if (present(coefs)) coefs = (/ cmplx(fac1,0.0_rk,kind=rk), cmplx(fac2,0.0_rk,kind=rk) /)
    else
      write(out, '(/a,1x,i3,1x,a)') 'rotme_cart_tens/wang_jktau_coefs error: rotational quantum number "k" =', k, 'is less than zero (expected >=0)'
      stop 'STOP'
    endif

  elseif (tau==1) then

    if (k==0) then
      ncoefs = 1
      if (present(coefs)) coefs = (/ cmplx(0.0_rk,1.0_rk,kind=rk), cmplx(0.0_rk,0.0_rk,kind=rk) /)
    elseif (k>0) then
      ncoefs = 2
      if (present(coefs)) coefs = (/ cmplx(0.0_rk,fac1,kind=rk), cmplx(0.0_rk,-fac2,kind=rk) /)
    else
      write(out, '(/a,1x,i3,1x,a)') 'rotme_cart_tens/wang_jktau_coefs error: rotational quantum number "k" =', k, 'is less than zero (expected >=0)'
      stop 'STOP'
    endif

  else

    write(out, '(/a,1x,i3,1x,a)') 'rotme_cart_tens/wang_jktau_coefs error: rotational quantum number "tau" =', tau, '(expected 0 or 1)'
    stop 'STOP'

  endif

end subroutine wang_jktau_coefs



!###################################################################################################################################



subroutine pseudoinverse(nrows, ncols, mat, invmat, tol_)

  integer(ik), intent(in) :: nrows, ncols
  complex(rk), intent(in) :: mat(nrows,ncols)
  complex(rk), intent(out) :: invmat(ncols,nrows)
  real(rk), intent(in), optional :: tol_

  integer(ik) lwork, info, i
  complex(rk) :: work1(1)
  complex(rk), allocatable :: work(:), matu(:,:), matvt(:,:), mat_d(:,:), matt(:,:)
  double precision :: tol, rwork(5*max(nrows,ncols))
  double precision, allocatable :: sv(:)

  if (present(tol_)) then
    tol = tol_
  else
    tol = 1.0d-14 !epsilon(1.0_rk)
  endif

  allocate(matt(nrows,ncols), sv(min(ncols,nrows)), matu(nrows,nrows), matvt(ncols,ncols), mat_d(ncols,nrows), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'rotbas/pseudoinverse error: failed to allocate matt(nrows,ncols), sv(min(ncols,nrows)), matu(nrows,nrows), matvt(ncols,ncols), mat_d(ncols,nrows)', &!
    'ncols, nrows =', ncols, nrows
    stop 'STOP'
  endif
  matt = 0
  sv = 0
  matu = 0
  matvt = 0
  mat_d = 0

  matt = mat

  ! estimate size of the workspace array (by setting lwork<0 to ZGESVD)
  lwork = -1
  call zgesvd('A', 'A', nrows, ncols, matt, nrows, sv, matu, nrows, matvt, ncols, work1, lwork, rwork, info)
  lwork = int(work1(1), kind=ik)

  ! allocate workspace array
  allocate(work(lwork), stat=info)
  if (info/=0) then
    write(out, '(/a,1x,i6)') 'rotbas/pseudoinverse error: failed to allocate workspace array required for SVD, size =', lwork
    stop 'STOP'
  endif

  ! perform actual SVD
  call zgesvd('A', 'A', nrows, ncols, matt, nrows, sv, matu, nrows, matvt, ncols, work, lwork, rwork, info)

  if (info/=0) then
    write(out, '(/a,1x,i6)') 'rotbas/pseudoinverse error: SVD failed, info =', info
    stop 'STOP'
  endif

  ! inverse diagonal matrix
  mat_d = 0.0
  do i=1, min(nrows,ncols)
    if (sv(i)>=tol) then
      mat_d(i,i) = 1.0_rk/sv(i)
    !else
    !  write(out, '(/a)') 'rotbas/pseudoinverse error: matrix is singular'
    !  write(out, '(a)') 'singular elements:'
    !  do j=1, min(nrows,ncols)
    !    write(out, '(1x,i3,1x,f)') j, sv(j)
    !  enddo
    !  stop 'STOP'
    endif
  enddo

  ! compute pseudoinverse
  mat_d = matmul(mat_d, conjg(transpose(matu)))
  invmat(1:ncols,1:nrows) = matmul(conjg(transpose(matvt)), mat_d)

  ! clear workspace
  deallocate(work,matt,mat_d,sv,matu,matvt)

end subroutine pseudoinverse



!###################################################################################################################################



function threej_symbol(j1,j2,j3, m1,m2,m3) result(f)

  integer(ik), intent(in) :: j1,j2,j3,m1,m2,m3
  real(rk) :: f

  integer(4) :: two, two_j1, two_j2, two_j3, two_m1, two_m2, two_m3

  two = 2_ik

  two_j1 = two*j1
  two_j2 = two*j2
  two_j3 = two*j3

  two_m1 = two*m1
  two_m2 = two*m2
  two_m3 = two*m3

#if (fwig>0)
  f = fwig3jj(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3)
#else
  f = three_j(j1,j2,j3, m1,m2,m3)
#endif

end function threej_symbol


end module rotme_cart_tens

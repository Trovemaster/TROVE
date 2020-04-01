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

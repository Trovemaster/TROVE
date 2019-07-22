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

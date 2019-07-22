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

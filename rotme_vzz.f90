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

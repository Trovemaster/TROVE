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

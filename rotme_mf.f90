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

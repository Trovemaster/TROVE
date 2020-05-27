! Contains external functions for computing spin-spin tensor of XY2-type molecule.
!
! 1. prop_xy2_spinspin_dipoleYY - computes dipolar spin-spin coupling tensor
!    between magnetic dipoles of Y1 and Y2 atoms. Tnesor depends only on
!    molecular geometry, thus no input expansion coefficients are required.
!    The magnetic dipole moments of Y atoms together with spins are provided in
!    the first expansion coefficient (of the first tensor element, i.e.
!    extF%coef(1,1)) as a prefactor mu0*(m1*muN/s1)*(m2*muN/s2)*1e30/(4*np.pi)/kHz_to_joule

module prop_xy2_spinspin
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xy2
  use timer

  implicit none

  public prop_xy2_spinspin_dipoleYY

contains


! Dipolar spin-spin coupling tensor between magnetic dipoles of Y1 and Y2 atoms in XY2 molecule

subroutine prop_xy2_spinspin_dipoleYY(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: i, j
  real(ark) :: x(natoms,3), tmat(3,3), dtens(3,3), r, n(3), prefac

  if (rank/=9) then
    write(out, '(/a,1x,i3,1x,a)') &
      'prop_xy2_spinspin_dipoleYY: rank of the dipole moment vector =', rank, ', expected 9'
    stop
  endif

  ! xyz are undefined for the local case
  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_spinspin_dipoleYY: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_spinspin_dipoleYY - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  prefac = extF%coef(1,1) ! mu0*(m1*muN/s1)*(m2*muN/s2)*1e30/(4*np.pi)/kHz_to_joule (final units are kHz*Angstrom**3)

  r = sqrt(sum( (x(2,:)-x(3,:))**2 )) ! distance between Y1 and Y2
  n = (x(2,:)-x(3,:))/r ! unit vector from Y1 to Y2

  do i=1, 3
    do j=1, 3
      if (i==j) then
        tmat(i,j) = 1.0_ark
      else
        tmat(i,j) = 0.0_ark
      endif
      tmat(i,j) = tmat(i,j) - 3.0_ark*n(i)*n(j)
    enddo
  enddo

  dtens = prefac/r**3*tmat

  f = (/dtens(1,1), dtens(1,2), dtens(1,3), &
        dtens(2,1), dtens(2,2), dtens(2,3), &
        dtens(3,1), dtens(3,2), dtens(3,3)/)

end subroutine prop_xy2_spinspin_dipoleYY


end module prop_xy2_spinspin

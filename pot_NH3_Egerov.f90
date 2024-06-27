!
!  This unit is for a user defined potential 
!
module pot_user
    use accuracy
    use moltype
  
    implicit none
  
    public MLdipole,MLpoten,ML_MEP,MLpoten_name
  
    private
   
    integer(ik), parameter :: verbose     = 4                          ! Verbosity level
    !
    contains
    !
    !
    function ML_MEP(dim,rho)  result(f)
  
     integer(ik),intent(in) ::  dim
     real(ark),intent(in)   ::  rho
     real(ark)              ::  f(dim)
     !
     if (dim/=3) stop 'Illegal size of the function - must be 3'
     !
     f(:) = molec%local_eq(:)
     f(molec%Ncoords) = rho
  
    end function ML_MEP
  
   ! Check the potential name 
   subroutine MLpoten_name(name)
     !
     character(len=cl),intent(in) ::  name
     character(len=cl),parameter ::  poten_name = 'GENERAL'
     ! 
     if (poten_name/=trim(name)) then
       write(out,"('a,a,a,a')") 'Wrong Potential ',trim(name),'; should be ',trim(poten_name)
     endif
     !
     write(out,"('a')") '  Using USER-tpye PES ',trim(poten_name)
     !
   end subroutine MLpoten_name
   !
  
   recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
     !
     integer(ik),intent(in) ::  rank,ncoords,natoms
     real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
     real(ark),intent(out)  ::  f(rank)
         !
         f = 0
         !
    end subroutine MLdipole
  
   !
   ! Defining potential energy function (built for SO2)
  
   function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
     !
     integer(ik),intent(in) ::  ncoords,natoms
     integer(ik)            ::  i, j
     real(ark),intent(in)   ::  local(ncoords)
     real(ark),intent(in)   ::  xyz(natoms,3)
     real(ark),intent(in)   ::  force(:)
     real(ark)              ::  S(ncoords), expansionIndices(940, 6), coefficients(940)
     real(ark)              ::  y1,y2,y3,y4,y5, r1, r2, r3, alpha1, alpha2, alpha3, phi1, phi2, phi3, addToPotential
     real(ark)              ::  f, b, rho_eq, r_eq, alpha_eq
     !
     ! f = MLpoten_xy3_morbid_Roman_10(ncoords,natoms,local,xyz,force)
     !
     b = 1.9_ark
     r_eq = 1.01099_ark
     alpha_eq = 106.7473616524_ark*2.0_ark*pi/360.0_ark
     rho_eq = 67.92572*pi/180.0_ark
     r1 = sqrt(dot_product(xyz(2, :) - xyz(1, :), xyz(2, :) - xyz(1, :)))
     r2 = sqrt(dot_product(xyz(3, :) - xyz(1, :), xyz(3, :) - xyz(1, :)))
     r3 = sqrt(dot_product(xyz(4, :) - xyz(1, :), xyz(4, :) - xyz(1, :)))
     y1 = 1.0_ark - exp(-b*(r1 - r_eq))
     y2 = 1.0_ark - exp(-b*(r2 - r_eq))
     y3 = 1.0_ark - exp(-b*(r3 - r_eq))
     S(1) = (y1 + y2 + y3)/sqrt(3.0_ark)
     S(3) = (2.0_ark*y1 - y2 - y3)/sqrt(6.0_ark)
     S(4) = (y2 - y3)/(sqrt(2.0_ark))
   !   phi1 = dot_product(xyz(1, :) - xyz(2, :), xyz(3, :) - xyz(1, :))/sqrt(abs(dot_product(xyz(1, :) - xyz(2, :), xyz(3, :) - xyz(1, :)))) - cos(alpha_eq)
   !   phi2 = dot_product(xyz(1, :) - xyz(3, :), xyz(4, :) - xyz(1, :))/sqrt(abs(dot_product(xyz(1, :) - xyz(3, :), xyz(4, :) - xyz(1, :)))) - cos(alpha_eq)
   !   phi3 = dot_product(xyz(1, :) - xyz(4, :), xyz(2, :) - xyz(1, :))/sqrt(abs(dot_product(xyz(1, :) - xyz(4, :), xyz(2, :) - xyz(1, :)))) - cos(alpha_eq)
     phi1 = dot_product(xyz(1, :) - xyz(3, :), xyz(4, :) - xyz(1, :))/(r2*r3) - cos(alpha_eq)
     phi2 = dot_product(xyz(1, :) - xyz(2, :), xyz(4, :) - xyz(1, :))/(r1*r3) - cos(alpha_eq)
     phi3 = dot_product(xyz(1, :) - xyz(3, :), xyz(4, :) - xyz(1, :))/(r2*r3) - cos(alpha_eq)
     S(5) = (2.0_ark*phi1 - phi2 - phi3)/sqrt(6.0_ark)
     S(6) = (phi2 - phi3)/sqrt(2.0_ark)
     alpha1 = acos(phi1 + cos(alpha_eq))
     alpha2 = acos(phi2 + cos(alpha_eq))
     alpha3 = acos(phi3 + cos(alpha_eq))
     S(2) = sin(rho_eq) - 2.0_ark*sin((alpha1 + alpha2 + alpha3)/6.0_ark)/sqrt(3.0_ark)

   !   write(*, '(A, F20.10)') "point ", sqrt(dot_product(xyz(1, :) - xyz(2, :), xyz(3, :) - xyz(1, :)))
     
     expansionIndices(1, :) = (/ 0, 0, 0, 0, 0, 2 /)
     coefficients(1) = 8.35221017607176990000e-002_ark
     expansionIndices(2, :) = (/ 0, 0, 0, 0, 2, 0 /)
     coefficients(2) = 8.35221017607176990000e-002_ark
     expansionIndices(3, :) = (/ 0, 0, 0, 1, 0, 1 /)
     coefficients(3) = 2.41234392625650910000e-002_ark
     expansionIndices(4, :) = (/ 0, 0, 1, 0, 1, 0 /)
     coefficients(4) = 2.41234392625650910000e-002_ark
     expansionIndices(5, :) = (/ 0, 0, 0, 2, 0, 0 /)
     coefficients(5) = 2.26714821838658600000e-001_ark
     expansionIndices(6, :) = (/ 0, 0, 2, 0, 0, 0 /)
     coefficients(6) = 2.26714821838658600000e-001_ark
     expansionIndices(7, :) = (/ 0, 2, 0, 0, 0, 0 /)
     coefficients(7) = 1.47975045171455700000e+000_ark
     expansionIndices(8, :) = (/ 1, 1, 0, 0, 0, 0 /)
     coefficients(8) = -3.03009271334584720000e-001_ark
     expansionIndices(9, :) = (/ 2, 0, 0, 0, 0, 0 /)
     coefficients(9) = 2.23104722925976220000e-001_ark
     expansionIndices(10, :) = (/ 0, 0, 0, 0, 1, 2 /)
     coefficients(10) = 2.38114496268895460000e-002_ark
     expansionIndices(11, :) = (/ 0, 0, 0, 0, 3, 0 /)
     coefficients(11) = -7.93714987562855160000e-003_ark
     expansionIndices(12, :) = (/ 0, 0, 0, 1, 1, 1 /)
     coefficients(12) = -5.97264673511189220000e-003_ark
     expansionIndices(13, :) = (/ 0, 0, 1, 0, 0, 2 /)
     coefficients(13) = -2.98632336755636810000e-003_ark
     expansionIndices(14, :) = (/ 0, 0, 1, 0, 2, 0 /)
     coefficients(14) = 2.98632336755636810000e-003_ark
     expansionIndices(15, :) = (/ 0, 0, 0, 2, 1, 0 /)
     coefficients(15) = -1.72121867972136020000e-002_ark
     expansionIndices(16, :) = (/ 0, 0, 1, 1, 0, 1 /)
     coefficients(16) = -3.44243735944223390000e-002_ark
     expansionIndices(17, :) = (/ 0, 0, 2, 0, 1, 0 /)
     coefficients(17) = 1.72121867972136020000e-002_ark
     expansionIndices(18, :) = (/ 0, 0, 1, 2, 0, 0 /)
     coefficients(18) = 3.26518978459615570000e-002_ark
     expansionIndices(19, :) = (/ 0, 0, 3, 0, 0, 0 /)
     coefficients(19) = -1.08839659486520760000e-002_ark
     expansionIndices(20, :) = (/ 0, 1, 0, 0, 0, 2 /)
     coefficients(20) = 1.90778666226051910000e-001_ark
     expansionIndices(21, :) = (/ 0, 1, 0, 0, 2, 0 /)
     coefficients(21) = 1.90778666226051910000e-001_ark
     expansionIndices(22, :) = (/ 0, 1, 0, 1, 0, 1 /)
     coefficients(22) = 9.22752008111396290000e-002_ark
     expansionIndices(23, :) = (/ 0, 1, 1, 0, 1, 0 /)
     coefficients(23) = 9.22752008111396290000e-002_ark
     expansionIndices(24, :) = (/ 0, 1, 0, 2, 0, 0 /)
     coefficients(24) = -1.07064550937149200000e-001_ark
     expansionIndices(25, :) = (/ 0, 1, 2, 0, 0, 0 /)
     coefficients(25) = -1.07064550937149200000e-001_ark
     expansionIndices(26, :) = (/ 0, 3, 0, 0, 0, 0 /)
     coefficients(26) = -1.86136692190152430000e+000_ark
     expansionIndices(27, :) = (/ 1, 0, 0, 0, 0, 2 /)
     coefficients(27) = -2.27143334754858520000e-002_ark
     expansionIndices(28, :) = (/ 1, 0, 0, 0, 2, 0 /)
     coefficients(28) = -2.27143334754858520000e-002_ark
     expansionIndices(29, :) = (/ 1, 0, 0, 1, 0, 1 /)
     coefficients(29) = -7.44935893754538870000e-004_ark
     expansionIndices(30, :) = (/ 1, 0, 1, 0, 1, 0 /)
     coefficients(30) = -7.44935893754538870000e-004_ark
     expansionIndices(31, :) = (/ 1, 0, 0, 2, 0, 0 /)
     coefficients(31) = -4.74405762626982500000e-002_ark
     expansionIndices(32, :) = (/ 1, 0, 2, 0, 0, 0 /)
     coefficients(32) = -4.74405762626982500000e-002_ark
     expansionIndices(33, :) = (/ 1, 2, 0, 0, 0, 0 /)
     coefficients(33) = 3.69804301781042450000e-001_ark
     expansionIndices(34, :) = (/ 2, 1, 0, 0, 0, 0 /)
     coefficients(34) = -6.45101124148464990000e-002_ark
     expansionIndices(35, :) = (/ 3, 0, 0, 0, 0, 0 /)
     coefficients(35) = -1.82072337155755080000e-002_ark
     expansionIndices(36, :) = (/ 0, 0, 0, 0, 0, 4 /)
     coefficients(36) = 1.88670719821665170000e-002_ark
     expansionIndices(37, :) = (/ 0, 0, 0, 0, 2, 2 /)
     coefficients(37) = 3.77341439643330350000e-002_ark
     expansionIndices(38, :) = (/ 0, 0, 0, 0, 4, 0 /)
     coefficients(38) = 1.88670719821665170000e-002_ark
     expansionIndices(39, :) = (/ 0, 0, 0, 1, 0, 3 /)
     coefficients(39) = 9.28173533522482640000e-003_ark
     expansionIndices(40, :) = (/ 0, 0, 0, 1, 2, 1 /)
     coefficients(40) = 9.28173533522482640000e-003_ark
     expansionIndices(41, :) = (/ 0, 0, 1, 0, 1, 2 /)
     coefficients(41) = 9.28173533522482640000e-003_ark
     expansionIndices(42, :) = (/ 0, 0, 1, 0, 3, 0 /)
     coefficients(42) = 9.28173533522482640000e-003_ark
     expansionIndices(43, :) = (/ 0, 0, 0, 2, 0, 2 /)
     coefficients(43) = -1.48794477288687680000e-002_ark
     expansionIndices(44, :) = (/ 0, 0, 0, 2, 2, 0 /)
     coefficients(44) = -1.48794477288687680000e-002_ark
     expansionIndices(45, :) = (/ 0, 0, 2, 0, 0, 2 /)
     coefficients(45) = -1.48794477288687680000e-002_ark
     expansionIndices(46, :) = (/ 0, 0, 2, 0, 2, 0 /)
     coefficients(46) = -1.48794477288687680000e-002_ark
     expansionIndices(47, :) = (/ 0, 0, 0, 2, 0, 2 /)
     coefficients(47) = 3.92202035243651660000e-003_ark
     expansionIndices(48, :) = (/ 0, 0, 0, 2, 2, 0 /)
     coefficients(48) = -3.92202035243651660000e-003_ark
     expansionIndices(49, :) = (/ 0, 0, 1, 1, 1, 1 /)
     coefficients(49) = 1.56880814097482870000e-002_ark
     expansionIndices(50, :) = (/ 0, 0, 2, 0, 0, 2 /)
     coefficients(50) = -3.92202035243651660000e-003_ark
     expansionIndices(51, :) = (/ 0, 0, 2, 0, 2, 0 /)
     coefficients(51) = 3.92202035243651660000e-003_ark
     expansionIndices(52, :) = (/ 0, 0, 0, 3, 0, 1 /)
     coefficients(52) = 1.13655408164072300000e-003_ark
     expansionIndices(53, :) = (/ 0, 0, 1, 2, 1, 0 /)
     coefficients(53) = 1.13655408164072300000e-003_ark
     expansionIndices(54, :) = (/ 0, 0, 2, 1, 0, 1 /)
     coefficients(54) = 1.13655408164072300000e-003_ark
     expansionIndices(55, :) = (/ 0, 0, 3, 0, 1, 0 /)
     coefficients(55) = 1.13655408164072300000e-003_ark
     expansionIndices(56, :) = (/ 0, 0, 0, 4, 0, 0 /)
     coefficients(56) = 9.67539818117143030000e-003_ark
     expansionIndices(57, :) = (/ 0, 0, 2, 2, 0, 0 /)
     coefficients(57) = 1.93507963623428610000e-002_ark
     expansionIndices(58, :) = (/ 0, 0, 4, 0, 0, 0 /)
     coefficients(58) = 9.67539818117143030000e-003_ark
     expansionIndices(59, :) = (/ 0, 1, 0, 0, 1, 2 /)
     coefficients(59) = 2.28179346733276730000e-002_ark
     expansionIndices(60, :) = (/ 0, 1, 0, 0, 3, 0 /)
     coefficients(60) = -7.60597822444131470000e-003_ark
     expansionIndices(61, :) = (/ 0, 1, 0, 1, 1, 1 /)
     coefficients(61) = -2.81702582722839110000e-002_ark
     expansionIndices(62, :) = (/ 0, 1, 1, 0, 0, 2 /)
     coefficients(62) = -1.40851291361439470000e-002_ark
     expansionIndices(63, :) = (/ 0, 1, 1, 0, 2, 0 /)
     coefficients(63) = 1.40851291361439470000e-002_ark
     expansionIndices(64, :) = (/ 0, 1, 0, 2, 1, 0 /)
     coefficients(64) = -4.35996392504400740000e-002_ark
     expansionIndices(65, :) = (/ 0, 1, 1, 1, 0, 1 /)
     coefficients(65) = -8.71992785008678250000e-002_ark
     expansionIndices(66, :) = (/ 0, 1, 2, 0, 1, 0 /)
     coefficients(66) = 4.35996392504400740000e-002_ark
     expansionIndices(67, :) = (/ 0, 1, 1, 2, 0, 0 /)
     coefficients(67) = 4.73032161489105040000e-002_ark
     expansionIndices(68, :) = (/ 0, 1, 3, 0, 0, 0 /)
     coefficients(68) = -1.57677387163009250000e-002_ark
     expansionIndices(69, :) = (/ 0, 2, 0, 0, 0, 2 /)
     coefficients(69) = -8.52697837059654410000e-002_ark
     expansionIndices(70, :) = (/ 0, 2, 0, 0, 2, 0 /)
     coefficients(70) = -8.52697837059654410000e-002_ark
     expansionIndices(71, :) = (/ 0, 2, 0, 1, 0, 1 /)
     coefficients(71) = 3.47690021700984790000e-001_ark
     expansionIndices(72, :) = (/ 0, 2, 1, 0, 1, 0 /)
     coefficients(72) = 3.47690021700984790000e-001_ark
     expansionIndices(73, :) = (/ 0, 2, 0, 2, 0, 0 /)
     coefficients(73) = 1.44456055465740810000e-001_ark
     expansionIndices(74, :) = (/ 0, 2, 2, 0, 0, 0 /)
     coefficients(74) = 1.44456055465740810000e-001_ark
     expansionIndices(75, :) = (/ 0, 4, 0, 0, 0, 0 /)
     coefficients(75) = 4.49998542662941500000e+000_ark
     expansionIndices(76, :) = (/ 1, 0, 0, 1, 1, 1 /)
     coefficients(76) = -6.06389954985904560000e-003_ark
     expansionIndices(77, :) = (/ 1, 0, 1, 0, 0, 2 /)
     coefficients(77) = -3.03194977492995130000e-003_ark
     expansionIndices(78, :) = (/ 1, 0, 1, 0, 2, 0 /)
     coefficients(78) = 3.03194977492995130000e-003_ark
     expansionIndices(79, :) = (/ 1, 0, 0, 2, 1, 0 /)
     coefficients(79) = -9.67464731475196730000e-003_ark
     expansionIndices(80, :) = (/ 1, 0, 1, 1, 0, 1 /)
     coefficients(80) = -1.93492946295012010000e-002_ark
     expansionIndices(81, :) = (/ 1, 0, 2, 0, 1, 0 /)
     coefficients(81) = 9.67464731475196730000e-003_ark
     expansionIndices(82, :) = (/ 1, 0, 1, 2, 0, 0 /)
     coefficients(82) = -4.94473842191388010000e-002_ark
     expansionIndices(83, :) = (/ 1, 0, 3, 0, 0, 0 /)
     coefficients(83) = 1.64824614063769080000e-002_ark
     expansionIndices(84, :) = (/ 1, 1, 0, 0, 0, 2 /)
     coefficients(84) = -4.71920319746126850000e-002_ark
     expansionIndices(85, :) = (/ 1, 1, 0, 0, 2, 0 /)
     coefficients(85) = -4.71920319746126850000e-002_ark
     expansionIndices(86, :) = (/ 1, 1, 0, 1, 0, 1 /)
     coefficients(86) = 4.79039953908221360000e-003_ark
     expansionIndices(87, :) = (/ 1, 1, 1, 0, 1, 0 /)
     coefficients(87) = 4.79039953908221360000e-003_ark
     expansionIndices(88, :) = (/ 1, 1, 0, 2, 0, 0 /)
     coefficients(88) = -7.27055389042235570000e-002_ark
     expansionIndices(89, :) = (/ 1, 1, 2, 0, 0, 0 /)
     coefficients(89) = -7.27055389042235570000e-002_ark
     expansionIndices(90, :) = (/ 1, 3, 0, 0, 0, 0 /)
     coefficients(90) = -2.67620394533707630000e+000_ark
     expansionIndices(91, :) = (/ 2, 0, 0, 0, 0, 2 /)
     coefficients(91) = -5.86646280290414970000e-003_ark
     expansionIndices(92, :) = (/ 2, 0, 0, 0, 2, 0 /)
     coefficients(92) = -5.86646280290414970000e-003_ark
     expansionIndices(93, :) = (/ 2, 0, 0, 1, 0, 1 /)
     coefficients(93) = -9.38124255797328530000e-003_ark
     expansionIndices(94, :) = (/ 2, 0, 1, 0, 1, 0 /)
     coefficients(94) = -9.38124255797328530000e-003_ark
     expansionIndices(95, :) = (/ 2, 0, 0, 2, 0, 0 /)
     coefficients(95) = 3.05587036649660820000e-002_ark
     expansionIndices(96, :) = (/ 2, 0, 2, 0, 0, 0 /)
     coefficients(96) = 3.05587036649660820000e-002_ark
     expansionIndices(97, :) = (/ 2, 2, 0, 0, 0, 0 /)
     coefficients(97) = 4.77934301196786200000e-001_ark
     expansionIndices(98, :) = (/ 3, 1, 0, 0, 0, 0 /)
     coefficients(98) = -7.73238497870463920000e-003_ark
     expansionIndices(99, :) = (/ 4, 0, 0, 0, 0, 0 /)
     coefficients(99) = 2.68948822731941840000e-003_ark
     expansionIndices(100, :) = (/ 0, 0, 0, 0, 1, 4 /)
     coefficients(100) = 1.20326902903677680000e-002_ark
     expansionIndices(101, :) = (/ 0, 0, 0, 0, 3, 2 /)
     coefficients(101) = 8.02179352691708560000e-003_ark
     expansionIndices(102, :) = (/ 0, 0, 0, 0, 5, 0 /)
     coefficients(102) = -4.01089676345697030000e-003_ark
     expansionIndices(103, :) = (/ 0, 0, 0, 1, 1, 3 /)
     coefficients(103) = 7.43267664786763630000e-003_ark
     expansionIndices(104, :) = (/ 0, 0, 0, 1, 3, 1 /)
     coefficients(104) = 7.43267664786763630000e-003_ark
     expansionIndices(105, :) = (/ 0, 0, 1, 0, 0, 4 /)
     coefficients(105) = 3.71633832393563870000e-003_ark
     expansionIndices(106, :) = (/ 0, 0, 1, 0, 4, 0 /)
     coefficients(106) = -3.71633832393563870000e-003_ark
     expansionIndices(107, :) = (/ 0, 0, 0, 1, 1, 3 /)
     coefficients(107) = 1.96326688064112420000e-003_ark
     expansionIndices(108, :) = (/ 0, 0, 0, 1, 3, 1 /)
     coefficients(108) = -1.96326688064112420000e-003_ark
     expansionIndices(109, :) = (/ 0, 0, 1, 0, 0, 4 /)
     coefficients(109) = -4.90816720160329090000e-004_ark
     expansionIndices(110, :) = (/ 0, 0, 1, 0, 2, 2 /)
     coefficients(110) = 2.94490032096216730000e-003_ark
     expansionIndices(111, :) = (/ 0, 0, 1, 0, 4, 0 /)
     coefficients(111) = -4.90816720160329090000e-004_ark
     expansionIndices(112, :) = (/ 0, 0, 0, 2, 1, 2 /)
     coefficients(112) = 2.74108134794312250000e-003_ark
     expansionIndices(113, :) = (/ 0, 0, 0, 2, 3, 0 /)
     coefficients(113) = -9.13693782648006010000e-004_ark
     expansionIndices(114, :) = (/ 0, 0, 2, 0, 1, 2 /)
     coefficients(114) = 2.74108134794312250000e-003_ark
     expansionIndices(115, :) = (/ 0, 0, 2, 0, 3, 0 /)
     coefficients(115) = -9.13693782648006010000e-004_ark
     expansionIndices(116, :) = (/ 0, 0, 0, 2, 1, 2 /)
     coefficients(116) = -5.44442401884904020000e-003_ark
     expansionIndices(117, :) = (/ 0, 0, 0, 2, 3, 0 /)
     coefficients(117) = -5.44442401884904020000e-003_ark
     expansionIndices(118, :) = (/ 0, 0, 1, 1, 0, 3 /)
     coefficients(118) = -1.08888480376980800000e-002_ark
     expansionIndices(119, :) = (/ 0, 0, 1, 1, 2, 1 /)
     coefficients(119) = -1.08888480376980800000e-002_ark
     expansionIndices(120, :) = (/ 0, 0, 2, 0, 1, 2 /)
     coefficients(120) = 5.44442401884904020000e-003_ark
     expansionIndices(121, :) = (/ 0, 0, 2, 0, 3, 0 /)
     coefficients(121) = 5.44442401884904020000e-003_ark
     expansionIndices(122, :) = (/ 0, 0, 0, 3, 1, 1 /)
     coefficients(122) = -9.15903033558053070000e-003_ark
     expansionIndices(123, :) = (/ 0, 0, 1, 2, 0, 2 /)
     coefficients(123) = -4.57951516779026530000e-003_ark
     expansionIndices(124, :) = (/ 0, 0, 1, 2, 2, 0 /)
     coefficients(124) = 4.57951516779026530000e-003_ark
     expansionIndices(125, :) = (/ 0, 0, 2, 1, 1, 1 /)
     coefficients(125) = -9.15903033558053070000e-003_ark
     expansionIndices(126, :) = (/ 0, 0, 3, 0, 0, 2 /)
     coefficients(126) = -4.57951516779026530000e-003_ark
     expansionIndices(127, :) = (/ 0, 0, 3, 0, 2, 0 /)
     coefficients(127) = 4.57951516779026530000e-003_ark
     expansionIndices(128, :) = (/ 0, 0, 1, 2, 0, 2 /)
     coefficients(128) = 1.49127291814066790000e-002_ark
     expansionIndices(129, :) = (/ 0, 0, 1, 2, 2, 0 /)
     coefficients(129) = 1.49127291814066790000e-002_ark
     expansionIndices(130, :) = (/ 0, 0, 3, 0, 0, 2 /)
     coefficients(130) = -4.97090972713718280000e-003_ark
     expansionIndices(131, :) = (/ 0, 0, 3, 0, 2, 0 /)
     coefficients(131) = -4.97090972713718280000e-003_ark
     expansionIndices(132, :) = (/ 0, 0, 0, 4, 1, 0 /)
     coefficients(132) = 1.65846315337923680000e-003_ark
     expansionIndices(133, :) = (/ 0, 0, 1, 3, 0, 1 /)
     coefficients(133) = 3.31692630675684860000e-003_ark
     expansionIndices(134, :) = (/ 0, 0, 3, 1, 0, 1 /)
     coefficients(134) = 3.31692630675684860000e-003_ark
     expansionIndices(135, :) = (/ 0, 0, 4, 0, 1, 0 /)
     coefficients(135) = -1.65846315337923680000e-003_ark
     expansionIndices(136, :) = (/ 0, 0, 0, 4, 1, 0 /)
     coefficients(136) = 5.48233995274915080000e-004_ark
     expansionIndices(137, :) = (/ 0, 0, 1, 3, 0, 1 /)
     coefficients(137) = -2.19293598109944560000e-003_ark
     expansionIndices(138, :) = (/ 0, 0, 2, 2, 1, 0 /)
     coefficients(138) = -3.28940397164970580000e-003_ark
     expansionIndices(139, :) = (/ 0, 0, 3, 1, 0, 1 /)
     coefficients(139) = 2.19293598109944560000e-003_ark
     expansionIndices(140, :) = (/ 0, 0, 4, 0, 1, 0 /)
     coefficients(140) = 5.48233995274915080000e-004_ark
     expansionIndices(141, :) = (/ 0, 0, 1, 4, 0, 0 /)
     coefficients(141) = -4.62101220783683090000e-003_ark
     expansionIndices(142, :) = (/ 0, 0, 3, 2, 0, 0 /)
     coefficients(142) = -3.08067480522656670000e-003_ark
     expansionIndices(143, :) = (/ 0, 0, 5, 0, 0, 0 /)
     coefficients(143) = 1.54033740261267970000e-003_ark
     expansionIndices(144, :) = (/ 0, 1, 0, 0, 0, 4 /)
     coefficients(144) = 3.75204396236362290000e-002_ark
     expansionIndices(145, :) = (/ 0, 1, 0, 0, 2, 2 /)
     coefficients(145) = 7.50408792472724570000e-002_ark
     expansionIndices(146, :) = (/ 0, 1, 0, 0, 4, 0 /)
     coefficients(146) = 3.75204396236362290000e-002_ark
     expansionIndices(147, :) = (/ 0, 1, 0, 1, 0, 3 /)
     coefficients(147) = 1.51640749002263930000e-002_ark
     expansionIndices(148, :) = (/ 0, 1, 0, 1, 2, 1 /)
     coefficients(148) = 1.51640749002263930000e-002_ark
     expansionIndices(149, :) = (/ 0, 1, 1, 0, 1, 2 /)
     coefficients(149) = 1.51640749002263930000e-002_ark
     expansionIndices(150, :) = (/ 0, 1, 1, 0, 3, 0 /)
     coefficients(150) = 1.51640749002263930000e-002_ark
     expansionIndices(151, :) = (/ 0, 1, 0, 2, 0, 2 /)
     coefficients(151) = -3.47714981852756540000e-002_ark
     expansionIndices(152, :) = (/ 0, 1, 0, 2, 2, 0 /)
     coefficients(152) = -3.47714981852756540000e-002_ark
     expansionIndices(153, :) = (/ 0, 1, 2, 0, 0, 2 /)
     coefficients(153) = -3.47714981852756540000e-002_ark
     expansionIndices(154, :) = (/ 0, 1, 2, 0, 2, 0 /)
     coefficients(154) = -3.47714981852756540000e-002_ark
     expansionIndices(155, :) = (/ 0, 1, 0, 2, 0, 2 /)
     coefficients(155) = 1.48906250300433820000e-002_ark
     expansionIndices(156, :) = (/ 0, 1, 0, 2, 2, 0 /)
     coefficients(156) = -1.48906250300433820000e-002_ark
     expansionIndices(157, :) = (/ 0, 1, 1, 1, 1, 1 /)
     coefficients(157) = 5.95625001201819520000e-002_ark
     expansionIndices(158, :) = (/ 0, 1, 2, 0, 0, 2 /)
     coefficients(158) = -1.48906250300433820000e-002_ark
     expansionIndices(159, :) = (/ 0, 1, 2, 0, 2, 0 /)
     coefficients(159) = 1.48906250300433820000e-002_ark
     expansionIndices(160, :) = (/ 0, 1, 0, 4, 0, 0 /)
     coefficients(160) = -3.67918674138864070000e-003_ark
     expansionIndices(161, :) = (/ 0, 1, 2, 2, 0, 0 /)
     coefficients(161) = -7.35837348277728130000e-003_ark
     expansionIndices(162, :) = (/ 0, 1, 4, 0, 0, 0 /)
     coefficients(162) = -3.67918674138864070000e-003_ark
     expansionIndices(163, :) = (/ 0, 2, 0, 0, 1, 2 /)
     coefficients(163) = -6.26217553583841010000e-001_ark
     expansionIndices(164, :) = (/ 0, 2, 0, 0, 3, 0 /)
     coefficients(164) = 2.08739184527947010000e-001_ark
     expansionIndices(165, :) = (/ 0, 2, 0, 1, 1, 1 /)
     coefficients(165) = -7.95543942460977260000e-003_ark
     expansionIndices(166, :) = (/ 0, 2, 1, 0, 0, 2 /)
     coefficients(166) = -3.97771971230488630000e-003_ark
     expansionIndices(167, :) = (/ 0, 2, 1, 0, 2, 0 /)
     coefficients(167) = 3.97771971230488630000e-003_ark
     expansionIndices(168, :) = (/ 0, 2, 0, 2, 1, 0 /)
     coefficients(168) = -1.74470147481012190000e-001_ark
     expansionIndices(169, :) = (/ 0, 2, 1, 1, 0, 1 /)
     coefficients(169) = -3.48940294962024390000e-001_ark
     expansionIndices(170, :) = (/ 0, 2, 2, 0, 1, 0 /)
     coefficients(170) = 1.74470147481012190000e-001_ark
     expansionIndices(171, :) = (/ 0, 2, 1, 2, 0, 0 /)
     coefficients(171) = -1.65132562509422210000e-001_ark
     expansionIndices(172, :) = (/ 0, 2, 3, 0, 0, 0 /)
     coefficients(172) = 5.50441875031407390000e-002_ark
     expansionIndices(173, :) = (/ 0, 3, 0, 0, 0, 2 /)
     coefficients(173) = 2.11043588551879550000e+000_ark
     expansionIndices(174, :) = (/ 0, 3, 0, 0, 2, 0 /)
     coefficients(174) = 2.11043588551879550000e+000_ark
     expansionIndices(175, :) = (/ 0, 3, 0, 1, 0, 1 /)
     coefficients(175) = -1.06732536055868920000e+000_ark
     expansionIndices(176, :) = (/ 0, 3, 1, 0, 1, 0 /)
     coefficients(176) = -1.06732536055868920000e+000_ark
     expansionIndices(177, :) = (/ 0, 3, 0, 2, 0, 0 /)
     coefficients(177) = -1.40927195695630060000e+000_ark
     expansionIndices(178, :) = (/ 0, 3, 2, 0, 0, 0 /)
     coefficients(178) = -1.40927195695630060000e+000_ark
     expansionIndices(179, :) = (/ 0, 5, 0, 0, 0, 0 /)
     coefficients(179) = -6.76831809371480730000e+000_ark
     expansionIndices(180, :) = (/ 1, 0, 0, 1, 0, 3 /)
     coefficients(180) = 4.29066736487255030000e-003_ark
     expansionIndices(181, :) = (/ 1, 0, 0, 1, 2, 1 /)
     coefficients(181) = 4.29066736487255030000e-003_ark
     expansionIndices(182, :) = (/ 1, 0, 1, 0, 1, 2 /)
     coefficients(182) = 4.29066736487255030000e-003_ark
     expansionIndices(183, :) = (/ 1, 0, 1, 0, 3, 0 /)
     coefficients(183) = 4.29066736487255030000e-003_ark
     expansionIndices(184, :) = (/ 1, 0, 0, 2, 0, 2 /)
     coefficients(184) = -1.60883041647021660000e-002_ark
     expansionIndices(185, :) = (/ 1, 0, 0, 2, 2, 0 /)
     coefficients(185) = -1.60883041647021660000e-002_ark
     expansionIndices(186, :) = (/ 1, 0, 2, 0, 0, 2 /)
     coefficients(186) = -1.60883041647021660000e-002_ark
     expansionIndices(187, :) = (/ 1, 0, 2, 0, 2, 0 /)
     coefficients(187) = -1.60883041647021660000e-002_ark
     expansionIndices(188, :) = (/ 1, 0, 0, 2, 0, 2 /)
     coefficients(188) = 4.75881172240139990000e-003_ark
     expansionIndices(189, :) = (/ 1, 0, 0, 2, 2, 0 /)
     coefficients(189) = -4.75881172240139990000e-003_ark
     expansionIndices(190, :) = (/ 1, 0, 1, 1, 1, 1 /)
     coefficients(190) = 1.90352468896082920000e-002_ark
     expansionIndices(191, :) = (/ 1, 0, 2, 0, 0, 2 /)
     coefficients(191) = -4.75881172240139990000e-003_ark
     expansionIndices(192, :) = (/ 1, 0, 2, 0, 2, 0 /)
     coefficients(192) = 4.75881172240139990000e-003_ark
     expansionIndices(193, :) = (/ 1, 0, 0, 3, 0, 1 /)
     coefficients(193) = -4.28661992390491630000e-003_ark
     expansionIndices(194, :) = (/ 1, 0, 1, 2, 1, 0 /)
     coefficients(194) = -4.28661992390491630000e-003_ark
     expansionIndices(195, :) = (/ 1, 0, 2, 1, 0, 1 /)
     coefficients(195) = -4.28661992390491630000e-003_ark
     expansionIndices(196, :) = (/ 1, 0, 3, 0, 1, 0 /)
     coefficients(196) = -4.28661992390491630000e-003_ark
     expansionIndices(197, :) = (/ 1, 0, 0, 4, 0, 0 /)
     coefficients(197) = 5.72252091552913170000e-003_ark
     expansionIndices(198, :) = (/ 1, 0, 2, 2, 0, 0 /)
     coefficients(198) = 1.14450418310582630000e-002_ark
     expansionIndices(199, :) = (/ 1, 0, 4, 0, 0, 0 /)
     coefficients(199) = 5.72252091552913170000e-003_ark
     expansionIndices(200, :) = (/ 1, 1, 0, 0, 1, 2 /)
     coefficients(200) = 4.56915255044300640000e-002_ark
     expansionIndices(201, :) = (/ 1, 1, 0, 0, 3, 0 /)
     coefficients(201) = -1.52305085014742010000e-002_ark
     expansionIndices(202, :) = (/ 1, 1, 0, 1, 1, 1 /)
     coefficients(202) = -2.89326564454336390000e-002_ark
     expansionIndices(203, :) = (/ 1, 1, 1, 0, 0, 2 /)
     coefficients(203) = -1.44663282227188630000e-002_ark
     expansionIndices(204, :) = (/ 1, 1, 1, 0, 2, 0 /)
     coefficients(204) = 1.44663282227188630000e-002_ark
     expansionIndices(205, :) = (/ 1, 1, 0, 2, 1, 0 /)
     coefficients(205) = -3.50687144335446800000e-002_ark
     expansionIndices(206, :) = (/ 1, 1, 1, 1, 0, 1 /)
     coefficients(206) = -7.01374288670794520000e-002_ark
     expansionIndices(207, :) = (/ 1, 1, 2, 0, 1, 0 /)
     coefficients(207) = 3.50687144335446800000e-002_ark
     expansionIndices(208, :) = (/ 1, 1, 1, 2, 0, 0 /)
     coefficients(208) = 4.52845539622235420000e-002_ark
     expansionIndices(209, :) = (/ 1, 1, 3, 0, 0, 0 /)
     coefficients(209) = -1.50948513207387160000e-002_ark
     expansionIndices(210, :) = (/ 1, 2, 0, 0, 0, 2 /)
     coefficients(210) = -9.43595784761765010000e-002_ark
     expansionIndices(211, :) = (/ 1, 2, 0, 0, 2, 0 /)
     coefficients(211) = -9.43595784761765010000e-002_ark
     expansionIndices(212, :) = (/ 1, 2, 0, 1, 0, 1 /)
     coefficients(212) = 5.17944544099733890000e-001_ark
     expansionIndices(213, :) = (/ 1, 2, 1, 0, 1, 0 /)
     coefficients(213) = 5.17944544099733890000e-001_ark
     expansionIndices(214, :) = (/ 1, 2, 0, 2, 0, 0 /)
     coefficients(214) = 3.25352154185234110000e-001_ark
     expansionIndices(215, :) = (/ 1, 2, 2, 0, 0, 0 /)
     coefficients(215) = 3.25352154185234110000e-001_ark
     expansionIndices(216, :) = (/ 1, 4, 0, 0, 0, 0 /)
     coefficients(216) = 8.94354508604022410000e+000_ark
     expansionIndices(217, :) = (/ 2, 0, 0, 0, 1, 2 /)
     coefficients(217) = 1.03217268954211070000e-003_ark
     expansionIndices(218, :) = (/ 2, 0, 0, 0, 3, 0 /)
     coefficients(218) = -3.44057563180703610000e-004_ark
     expansionIndices(219, :) = (/ 2, 0, 0, 1, 1, 1 /)
     coefficients(219) = -9.97451557803276210000e-003_ark
     expansionIndices(220, :) = (/ 2, 0, 1, 0, 0, 2 /)
     coefficients(220) = -4.98725778901638100000e-003_ark
     expansionIndices(221, :) = (/ 2, 0, 1, 0, 2, 0 /)
     coefficients(221) = 4.98725778901638100000e-003_ark
     expansionIndices(222, :) = (/ 2, 0, 1, 2, 0, 0 /)
     coefficients(222) = -8.61422098240238050000e-003_ark
     expansionIndices(223, :) = (/ 2, 0, 3, 0, 0, 0 /)
     coefficients(223) = 2.87140699413412700000e-003_ark
     expansionIndices(224, :) = (/ 2, 1, 0, 0, 0, 2 /)
     coefficients(224) = -5.54050756036391580000e-002_ark
     expansionIndices(225, :) = (/ 2, 1, 0, 0, 2, 0 /)
     coefficients(225) = -5.54050756036391580000e-002_ark
     expansionIndices(226, :) = (/ 2, 1, 0, 1, 0, 1 /)
     coefficients(226) = -1.67210431511704700000e-002_ark
     expansionIndices(227, :) = (/ 2, 1, 1, 0, 1, 0 /)
     coefficients(227) = -1.67210431511704700000e-002_ark
     expansionIndices(228, :) = (/ 2, 1, 0, 2, 0, 0 /)
     coefficients(228) = -1.97204314402776890000e-002_ark
     expansionIndices(229, :) = (/ 2, 1, 2, 0, 0, 0 /)
     coefficients(229) = -1.97204314402776890000e-002_ark
     expansionIndices(230, :) = (/ 2, 3, 0, 0, 0, 0 /)
     coefficients(230) = -2.90811756581810780000e+000_ark
     expansionIndices(231, :) = (/ 3, 0, 0, 1, 0, 1 /)
     coefficients(231) = -9.51954960134395250000e-003_ark
     expansionIndices(232, :) = (/ 3, 0, 1, 0, 1, 0 /)
     coefficients(232) = -9.51954960134395250000e-003_ark
     expansionIndices(233, :) = (/ 3, 0, 0, 2, 0, 0 /)
     coefficients(233) = 2.80476830353233950000e-003_ark
     expansionIndices(234, :) = (/ 3, 0, 2, 0, 0, 0 /)
     coefficients(234) = 2.80476830353233950000e-003_ark
     expansionIndices(235, :) = (/ 3, 2, 0, 0, 0, 0 /)
     coefficients(235) = 2.92380005840328360000e-001_ark
     expansionIndices(236, :) = (/ 0, 0, 0, 0, 0, 6 /)
     coefficients(236) = 4.87517249418921090000e-003_ark
     expansionIndices(237, :) = (/ 0, 0, 0, 0, 2, 4 /)
     coefficients(237) = 1.46255174825699710000e-002_ark
     expansionIndices(238, :) = (/ 0, 0, 0, 0, 4, 2 /)
     coefficients(238) = 1.46255174825699710000e-002_ark
     expansionIndices(239, :) = (/ 0, 0, 0, 0, 6, 0 /)
     coefficients(239) = 4.87517249418921090000e-003_ark
     expansionIndices(240, :) = (/ 0, 0, 0, 0, 0, 6 /)
     coefficients(240) = -2.97439814636483420000e-004_ark
     expansionIndices(241, :) = (/ 0, 0, 0, 0, 2, 4 /)
     coefficients(241) = 4.46159721954725130000e-003_ark
     expansionIndices(242, :) = (/ 0, 0, 0, 0, 4, 2 /)
     coefficients(242) = -4.46159721954725130000e-003_ark
     expansionIndices(243, :) = (/ 0, 0, 0, 0, 6, 0 /)
     coefficients(243) = 2.97439814636483420000e-004_ark
     expansionIndices(244, :) = (/ 0, 0, 0, 1, 0, 5 /)
     coefficients(244) = 4.70981426626130610000e-003_ark
     expansionIndices(245, :) = (/ 0, 0, 0, 1, 2, 3 /)
     coefficients(245) = 9.41962853252630540000e-003_ark
     expansionIndices(246, :) = (/ 0, 0, 0, 1, 4, 1 /)
     coefficients(246) = 4.70981426626130610000e-003_ark
     expansionIndices(247, :) = (/ 0, 0, 1, 0, 1, 4 /)
     coefficients(247) = 4.70981426626130610000e-003_ark
     expansionIndices(248, :) = (/ 0, 0, 1, 0, 3, 2 /)
     coefficients(248) = 9.41962853252630540000e-003_ark
     expansionIndices(249, :) = (/ 0, 0, 1, 0, 5, 0 /)
     coefficients(249) = 4.70981426626130610000e-003_ark
     expansionIndices(250, :) = (/ 0, 0, 0, 1, 0, 5 /)
     coefficients(250) = -1.94345861175421570000e-004_ark
     expansionIndices(251, :) = (/ 0, 0, 0, 1, 2, 3 /)
     coefficients(251) = 1.94345861175421550000e-003_ark
     expansionIndices(252, :) = (/ 0, 0, 0, 1, 4, 1 /)
     coefficients(252) = -9.71729305876866860000e-004_ark
     expansionIndices(253, :) = (/ 0, 0, 1, 0, 1, 4 /)
     coefficients(253) = 9.71729305876866860000e-004_ark
     expansionIndices(254, :) = (/ 0, 0, 1, 0, 3, 2 /)
     coefficients(254) = -1.94345861175421550000e-003_ark
     expansionIndices(255, :) = (/ 0, 0, 1, 0, 5, 0 /)
     coefficients(255) = 1.94345861175421570000e-004_ark
     expansionIndices(256, :) = (/ 0, 0, 0, 2, 0, 4 /)
     coefficients(256) = -8.72960181818371110000e-004_ark
     expansionIndices(257, :) = (/ 0, 0, 0, 2, 4, 0 /)
     coefficients(257) = 8.72960181818371110000e-004_ark
     expansionIndices(258, :) = (/ 0, 0, 1, 1, 1, 3 /)
     coefficients(258) = -3.49184072727314270000e-003_ark
     expansionIndices(259, :) = (/ 0, 0, 1, 1, 3, 1 /)
     coefficients(259) = -3.49184072727314270000e-003_ark
     expansionIndices(260, :) = (/ 0, 0, 2, 0, 0, 4 /)
     coefficients(260) = 8.72960181818371110000e-004_ark
     expansionIndices(261, :) = (/ 0, 0, 2, 0, 4, 0 /)
     coefficients(261) = -8.72960181818371110000e-004_ark
     expansionIndices(262, :) = (/ 0, 0, 0, 4, 0, 2 /)
     coefficients(262) = -3.76737521525158510000e-003_ark
     expansionIndices(263, :) = (/ 0, 0, 0, 4, 2, 0 /)
     coefficients(263) = -3.76737521525158510000e-003_ark
     expansionIndices(264, :) = (/ 0, 0, 2, 2, 0, 2 /)
     coefficients(264) = -7.53475043050317020000e-003_ark
     expansionIndices(265, :) = (/ 0, 0, 2, 2, 2, 0 /)
     coefficients(265) = -7.53475043050317020000e-003_ark
     expansionIndices(266, :) = (/ 0, 0, 4, 0, 0, 2 /)
     coefficients(266) = -3.76737521525158510000e-003_ark
     expansionIndices(267, :) = (/ 0, 0, 4, 0, 2, 0 /)
     coefficients(267) = -3.76737521525158510000e-003_ark
     expansionIndices(268, :) = (/ 0, 0, 0, 4, 0, 2 /)
     coefficients(268) = 1.39873652931037270000e-003_ark
     expansionIndices(269, :) = (/ 0, 0, 0, 4, 2, 0 /)
     coefficients(269) = -1.39873652931037270000e-003_ark
     expansionIndices(270, :) = (/ 0, 0, 1, 3, 1, 1 /)
     coefficients(270) = 5.59494611724094350000e-003_ark
     expansionIndices(271, :) = (/ 0, 0, 3, 1, 1, 1 /)
     coefficients(271) = 5.59494611724094350000e-003_ark
     expansionIndices(272, :) = (/ 0, 0, 4, 0, 0, 2 /)
     coefficients(272) = -1.39873652931037270000e-003_ark
     expansionIndices(273, :) = (/ 0, 0, 4, 0, 2, 0 /)
     coefficients(273) = 1.39873652931037270000e-003_ark
     expansionIndices(274, :) = (/ 0, 0, 0, 5, 0, 1 /)
     coefficients(274) = -1.22031941228223770000e-003_ark
     expansionIndices(275, :) = (/ 0, 0, 1, 4, 1, 0 /)
     coefficients(275) = -1.22031941228223770000e-003_ark
     expansionIndices(276, :) = (/ 0, 0, 2, 3, 0, 1 /)
     coefficients(276) = -2.44063882456543220000e-003_ark
     expansionIndices(277, :) = (/ 0, 0, 3, 2, 1, 0 /)
     coefficients(277) = -2.44063882456543220000e-003_ark
     expansionIndices(278, :) = (/ 0, 0, 4, 1, 0, 1 /)
     coefficients(278) = -1.22031941228223770000e-003_ark
     expansionIndices(279, :) = (/ 0, 0, 5, 0, 1, 0 /)
     coefficients(279) = -1.22031941228223770000e-003_ark
     expansionIndices(280, :) = (/ 0, 0, 0, 6, 0, 0 /)
     coefficients(280) = 7.31589322065984500000e-004_ark
     expansionIndices(281, :) = (/ 0, 0, 2, 4, 0, 0 /)
     coefficients(281) = 2.19476796619830430000e-003_ark
     expansionIndices(282, :) = (/ 0, 0, 4, 2, 0, 0 /)
     coefficients(282) = 2.19476796619830430000e-003_ark
     expansionIndices(283, :) = (/ 0, 0, 6, 0, 0, 0 /)
     coefficients(283) = 7.31589322065984500000e-004_ark
     expansionIndices(284, :) = (/ 0, 0, 0, 4, 0, 2 /)
     coefficients(284) = -1.41489343435579930000e-003_ark
     expansionIndices(285, :) = (/ 0, 0, 0, 4, 2, 0 /)
     coefficients(285) = 1.41489343435579930000e-003_ark
     expansionIndices(286, :) = (/ 0, 0, 1, 3, 1, 1 /)
     coefficients(286) = 1.13191474748452860000e-002_ark
     expansionIndices(287, :) = (/ 0, 0, 2, 2, 0, 2 /)
     coefficients(287) = 8.48936060613257780000e-003_ark
     expansionIndices(288, :) = (/ 0, 0, 2, 2, 2, 0 /)
     coefficients(288) = -8.48936060613257780000e-003_ark
     expansionIndices(289, :) = (/ 0, 0, 3, 1, 1, 1 /)
     coefficients(289) = -1.13191474748452860000e-002_ark
     expansionIndices(290, :) = (/ 0, 0, 4, 0, 0, 2 /)
     coefficients(290) = -1.41489343435579930000e-003_ark
     expansionIndices(291, :) = (/ 0, 0, 4, 0, 2, 0 /)
     coefficients(291) = 1.41489343435579930000e-003_ark
     expansionIndices(292, :) = (/ 0, 1, 0, 1, 1, 3 /)
     coefficients(292) = -7.77628177255017170000e-003_ark
     expansionIndices(293, :) = (/ 0, 1, 0, 1, 3, 1 /)
     coefficients(293) = -7.77628177255017170000e-003_ark
     expansionIndices(294, :) = (/ 0, 1, 1, 0, 0, 4 /)
     coefficients(294) = -3.88814088627699100000e-003_ark
     expansionIndices(295, :) = (/ 0, 1, 1, 0, 4, 0 /)
     coefficients(295) = 3.88814088627699100000e-003_ark
     expansionIndices(296, :) = (/ 0, 1, 0, 1, 1, 3 /)
     coefficients(296) = -1.85197325318839920000e-002_ark
     expansionIndices(297, :) = (/ 0, 1, 0, 1, 3, 1 /)
     coefficients(297) = 1.85197325318839920000e-002_ark
     expansionIndices(298, :) = (/ 0, 1, 1, 0, 0, 4 /)
     coefficients(298) = 4.62993313297145170000e-003_ark
     expansionIndices(299, :) = (/ 0, 1, 1, 0, 2, 2 /)
     coefficients(299) = -2.77795987978305260000e-002_ark
     expansionIndices(300, :) = (/ 0, 1, 1, 0, 4, 0 /)
     coefficients(300) = 4.62993313297145170000e-003_ark
     expansionIndices(301, :) = (/ 0, 1, 0, 3, 1, 1 /)
     coefficients(301) = -2.68828546428699320000e-002_ark
     expansionIndices(302, :) = (/ 0, 1, 1, 2, 0, 2 /)
     coefficients(302) = -1.34414273214349660000e-002_ark
     expansionIndices(303, :) = (/ 0, 1, 1, 2, 2, 0 /)
     coefficients(303) = 1.34414273214349660000e-002_ark
     expansionIndices(304, :) = (/ 0, 1, 2, 1, 1, 1 /)
     coefficients(304) = -2.68828546428699320000e-002_ark
     expansionIndices(305, :) = (/ 0, 1, 3, 0, 0, 2 /)
     coefficients(305) = -1.34414273214349660000e-002_ark
     expansionIndices(306, :) = (/ 0, 1, 3, 0, 2, 0 /)
     coefficients(306) = 1.34414273214349660000e-002_ark
     expansionIndices(307, :) = (/ 0, 1, 1, 2, 0, 2 /)
     coefficients(307) = 2.77935425655557000000e-002_ark
     expansionIndices(308, :) = (/ 0, 1, 1, 2, 2, 0 /)
     coefficients(308) = 2.77935425655557000000e-002_ark
     expansionIndices(309, :) = (/ 0, 1, 3, 0, 0, 2 /)
     coefficients(309) = -9.26451418852159220000e-003_ark
     expansionIndices(310, :) = (/ 0, 1, 3, 0, 2, 0 /)
     coefficients(310) = -9.26451418852159220000e-003_ark
     expansionIndices(311, :) = (/ 0, 2, 0, 0, 0, 4 /)
     coefficients(311) = -2.81878818479547220000e-002_ark
     expansionIndices(312, :) = (/ 0, 2, 0, 0, 2, 2 /)
     coefficients(312) = -5.63757636958966900000e-002_ark
     expansionIndices(313, :) = (/ 0, 2, 0, 0, 4, 0 /)
     coefficients(313) = -2.81878818479547220000e-002_ark
     expansionIndices(314, :) = (/ 0, 2, 0, 1, 0, 3 /)
     coefficients(314) = 1.06666365635360930000e-001_ark
     expansionIndices(315, :) = (/ 0, 2, 0, 1, 2, 1 /)
     coefficients(315) = 1.06666365635360930000e-001_ark
     expansionIndices(316, :) = (/ 0, 2, 1, 0, 1, 2 /)
     coefficients(316) = 1.06666365635360930000e-001_ark
     expansionIndices(317, :) = (/ 0, 2, 1, 0, 3, 0 /)
     coefficients(317) = 1.06666365635360930000e-001_ark
     expansionIndices(318, :) = (/ 0, 2, 0, 2, 0, 2 /)
     coefficients(318) = 6.51537750069275390000e-002_ark
     expansionIndices(319, :) = (/ 0, 2, 0, 2, 2, 0 /)
     coefficients(319) = 6.51537750069275390000e-002_ark
     expansionIndices(320, :) = (/ 0, 2, 2, 0, 0, 2 /)
     coefficients(320) = 6.51537750069275390000e-002_ark
     expansionIndices(321, :) = (/ 0, 2, 2, 0, 2, 0 /)
     coefficients(321) = 6.51537750069275390000e-002_ark
     expansionIndices(322, :) = (/ 0, 2, 0, 2, 0, 2 /)
     coefficients(322) = 2.53180306077172910000e-002_ark
     expansionIndices(323, :) = (/ 0, 2, 0, 2, 2, 0 /)
     coefficients(323) = -2.53180306077172910000e-002_ark
     expansionIndices(324, :) = (/ 0, 2, 1, 1, 1, 1 /)
     coefficients(324) = 1.01272122430869170000e-001_ark
     expansionIndices(325, :) = (/ 0, 2, 2, 0, 0, 2 /)
     coefficients(325) = -2.53180306077172910000e-002_ark
     expansionIndices(326, :) = (/ 0, 2, 2, 0, 2, 0 /)
     coefficients(326) = 2.53180306077172910000e-002_ark
     expansionIndices(327, :) = (/ 0, 2, 0, 3, 0, 1 /)
     coefficients(327) = 1.94575941401654130000e-001_ark
     expansionIndices(328, :) = (/ 0, 2, 1, 2, 1, 0 /)
     coefficients(328) = 1.94575941401654130000e-001_ark
     expansionIndices(329, :) = (/ 0, 2, 2, 1, 0, 1 /)
     coefficients(329) = 1.94575941401654130000e-001_ark
     expansionIndices(330, :) = (/ 0, 2, 3, 0, 1, 0 /)
     coefficients(330) = 1.94575941401654130000e-001_ark
     expansionIndices(331, :) = (/ 0, 2, 0, 4, 0, 0 /)
     coefficients(331) = 1.07145418624116880000e-001_ark
     expansionIndices(332, :) = (/ 0, 2, 2, 2, 0, 0 /)
     coefficients(332) = 2.14290837248185260000e-001_ark
     expansionIndices(333, :) = (/ 0, 2, 4, 0, 0, 0 /)
     coefficients(333) = 1.07145418624116880000e-001_ark
     expansionIndices(334, :) = (/ 0, 3, 0, 0, 1, 2 /)
     coefficients(334) = 4.49752810261089040000e+000_ark
     expansionIndices(335, :) = (/ 0, 3, 0, 0, 3, 0 /)
     coefficients(335) = -1.49917603420339040000e+000_ark
     expansionIndices(336, :) = (/ 0, 3, 0, 1, 1, 1 /)
     coefficients(336) = -8.59671871601152770000e-001_ark
     expansionIndices(337, :) = (/ 0, 3, 1, 0, 0, 2 /)
     coefficients(337) = -4.29835935800576390000e-001_ark
     expansionIndices(338, :) = (/ 0, 3, 1, 0, 2, 0 /)
     coefficients(338) = 4.29835935800576390000e-001_ark
     expansionIndices(339, :) = (/ 0, 3, 0, 2, 1, 0 /)
     coefficients(339) = 7.56836420607795660000e-001_ark
     expansionIndices(340, :) = (/ 0, 3, 1, 1, 0, 1 /)
     coefficients(340) = 1.51367284121559130000e+000_ark
     expansionIndices(341, :) = (/ 0, 3, 2, 0, 1, 0 /)
     coefficients(341) = -7.56836420607795660000e-001_ark
     expansionIndices(342, :) = (/ 0, 3, 1, 2, 0, 0 /)
     coefficients(342) = 1.79063770706788650000e+000_ark
     expansionIndices(343, :) = (/ 0, 3, 3, 0, 0, 0 /)
     coefficients(343) = -5.96879235689200090000e-001_ark
     expansionIndices(344, :) = (/ 0, 4, 0, 0, 0, 2 /)
     coefficients(344) = -1.19371239857027600000e+001_ark
     expansionIndices(345, :) = (/ 0, 4, 0, 0, 2, 0 /)
     coefficients(345) = -1.19371239857027600000e+001_ark
     expansionIndices(346, :) = (/ 0, 4, 0, 1, 0, 1 /)
     coefficients(346) = 4.96944974747279300000e+000_ark
     expansionIndices(347, :) = (/ 0, 4, 1, 0, 1, 0 /)
     coefficients(347) = 4.96944974747279300000e+000_ark
     expansionIndices(348, :) = (/ 0, 4, 0, 2, 0, 0 /)
     coefficients(348) = 5.70557570936992780000e+000_ark
     expansionIndices(349, :) = (/ 0, 4, 2, 0, 0, 0 /)
     coefficients(349) = 5.70557570936992780000e+000_ark
     expansionIndices(350, :) = (/ 0, 6, 0, 0, 0, 0 /)
     coefficients(350) = 1.53719741416340090000e+001_ark
     expansionIndices(351, :) = (/ 1, 0, 0, 0, 1, 4 /)
     coefficients(351) = 1.42137325116302600000e-002_ark
     expansionIndices(352, :) = (/ 1, 0, 0, 0, 3, 2 /)
     coefficients(352) = 9.47582167442636430000e-003_ark
     expansionIndices(353, :) = (/ 1, 0, 0, 0, 5, 0 /)
     coefficients(353) = -4.73791083721132510000e-003_ark
     expansionIndices(354, :) = (/ 1, 0, 0, 1, 1, 3 /)
     coefficients(354) = 1.25754970812496660000e-002_ark
     expansionIndices(355, :) = (/ 1, 0, 0, 1, 3, 1 /)
     coefficients(355) = 1.25754970812496660000e-002_ark
     expansionIndices(356, :) = (/ 1, 0, 1, 0, 0, 4 /)
     coefficients(356) = 6.28774854062791410000e-003_ark
     expansionIndices(357, :) = (/ 1, 0, 1, 0, 4, 0 /)
     coefficients(357) = -6.28774854062791410000e-003_ark
     expansionIndices(358, :) = (/ 1, 0, 0, 1, 1, 3 /)
     coefficients(358) = 2.76508349420317730000e-003_ark
     expansionIndices(359, :) = (/ 1, 0, 0, 1, 3, 1 /)
     coefficients(359) = -2.76508349420317730000e-003_ark
     expansionIndices(360, :) = (/ 1, 0, 1, 0, 0, 4 /)
     coefficients(360) = -6.91270873550861980000e-004_ark
     expansionIndices(361, :) = (/ 1, 0, 1, 0, 2, 2 /)
     coefficients(361) = 4.14762524130544290000e-003_ark
     expansionIndices(362, :) = (/ 1, 0, 1, 0, 4, 0 /)
     coefficients(362) = -6.91270873550861980000e-004_ark
     expansionIndices(363, :) = (/ 1, 0, 0, 2, 1, 2 /)
     coefficients(363) = -8.62985729100212810000e-003_ark
     expansionIndices(364, :) = (/ 1, 0, 0, 2, 3, 0 /)
     coefficients(364) = -8.62985729100212810000e-003_ark
     expansionIndices(365, :) = (/ 1, 0, 1, 1, 0, 3 /)
     coefficients(365) = -1.72597145820042560000e-002_ark
     expansionIndices(366, :) = (/ 1, 0, 1, 1, 2, 1 /)
     coefficients(366) = -1.72597145820042560000e-002_ark
     expansionIndices(367, :) = (/ 1, 0, 2, 0, 1, 2 /)
     coefficients(367) = 8.62985729100212810000e-003_ark
     expansionIndices(368, :) = (/ 1, 0, 2, 0, 3, 0 /)
     coefficients(368) = 8.62985729100212810000e-003_ark
     expansionIndices(369, :) = (/ 1, 0, 0, 3, 1, 1 /)
     coefficients(369) = -1.08215995783738730000e-002_ark
     expansionIndices(370, :) = (/ 1, 0, 1, 2, 0, 2 /)
     coefficients(370) = -5.41079978918693630000e-003_ark
     expansionIndices(371, :) = (/ 1, 0, 1, 2, 2, 0 /)
     coefficients(371) = 5.41079978918693630000e-003_ark
     expansionIndices(372, :) = (/ 1, 0, 2, 1, 1, 1 /)
     coefficients(372) = -1.08215995783738730000e-002_ark
     expansionIndices(373, :) = (/ 1, 0, 3, 0, 0, 2 /)
     coefficients(373) = -5.41079978918693630000e-003_ark
     expansionIndices(374, :) = (/ 1, 0, 3, 0, 2, 0 /)
     coefficients(374) = 5.41079978918693630000e-003_ark
     expansionIndices(375, :) = (/ 1, 0, 1, 2, 0, 2 /)
     coefficients(375) = 2.82942323124511680000e-002_ark
     expansionIndices(376, :) = (/ 1, 0, 1, 2, 2, 0 /)
     coefficients(376) = 2.82942323124511680000e-002_ark
     expansionIndices(377, :) = (/ 1, 0, 3, 0, 0, 2 /)
     coefficients(377) = -9.43141077082013630000e-003_ark
     expansionIndices(378, :) = (/ 1, 0, 3, 0, 2, 0 /)
     coefficients(378) = -9.43141077082013630000e-003_ark
     expansionIndices(379, :) = (/ 1, 0, 0, 4, 1, 0 /)
     coefficients(379) = 3.34607804593474360000e-003_ark
     expansionIndices(380, :) = (/ 1, 0, 1, 3, 0, 1 /)
     coefficients(380) = 6.69215609186620850000e-003_ark
     expansionIndices(381, :) = (/ 1, 0, 3, 1, 0, 1 /)
     coefficients(381) = 6.69215609186620850000e-003_ark
     expansionIndices(382, :) = (/ 1, 0, 4, 0, 1, 0 /)
     coefficients(382) = -3.34607804593474360000e-003_ark
     expansionIndices(383, :) = (/ 1, 0, 1, 4, 0, 0 /)
     coefficients(383) = -5.46776733046294530000e-003_ark
     expansionIndices(384, :) = (/ 1, 0, 3, 2, 0, 0 /)
     coefficients(384) = -3.64517822031101180000e-003_ark
     expansionIndices(385, :) = (/ 1, 0, 5, 0, 0, 0 /)
     coefficients(385) = 1.82258911015479140000e-003_ark
     expansionIndices(386, :) = (/ 1, 1, 0, 0, 0, 4 /)
     coefficients(386) = -1.55958538137003190000e-002_ark
     expansionIndices(387, :) = (/ 1, 1, 0, 0, 2, 2 /)
     coefficients(387) = -3.11917076274006370000e-002_ark
     expansionIndices(388, :) = (/ 1, 1, 0, 0, 4, 0 /)
     coefficients(388) = -1.55958538137003190000e-002_ark
     expansionIndices(389, :) = (/ 1, 1, 0, 2, 0, 2 /)
     coefficients(389) = -5.50758555598481750000e-002_ark
     expansionIndices(390, :) = (/ 1, 1, 0, 2, 2, 0 /)
     coefficients(390) = -5.50758555598481750000e-002_ark
     expansionIndices(391, :) = (/ 1, 1, 2, 0, 0, 2 /)
     coefficients(391) = -5.50758555598481750000e-002_ark
     expansionIndices(392, :) = (/ 1, 1, 2, 0, 2, 0 /)
     coefficients(392) = -5.50758555598481750000e-002_ark
     expansionIndices(393, :) = (/ 1, 1, 0, 2, 0, 2 /)
     coefficients(393) = 1.53849524437796900000e-002_ark
     expansionIndices(394, :) = (/ 1, 1, 0, 2, 2, 0 /)
     coefficients(394) = -1.53849524437796900000e-002_ark
     expansionIndices(395, :) = (/ 1, 1, 1, 1, 1, 1 /)
     coefficients(395) = 6.15398097751274700000e-002_ark
     expansionIndices(396, :) = (/ 1, 1, 2, 0, 0, 2 /)
     coefficients(396) = -1.53849524437796900000e-002_ark
     expansionIndices(397, :) = (/ 1, 1, 2, 0, 2, 0 /)
     coefficients(397) = 1.53849524437796900000e-002_ark
     expansionIndices(398, :) = (/ 1, 1, 0, 4, 0, 0 /)
     coefficients(398) = 1.10263035699957080000e-002_ark
     expansionIndices(399, :) = (/ 1, 1, 2, 2, 0, 0 /)
     coefficients(399) = 2.20526071399914160000e-002_ark
     expansionIndices(400, :) = (/ 1, 1, 4, 0, 0, 0 /)
     coefficients(400) = 1.10263035699957080000e-002_ark
     expansionIndices(401, :) = (/ 1, 2, 0, 0, 1, 2 /)
     coefficients(401) = 9.21212522452078820000e-002_ark
     expansionIndices(402, :) = (/ 1, 2, 0, 0, 3, 0 /)
     coefficients(402) = -3.07070840817359610000e-002_ark
     expansionIndices(403, :) = (/ 1, 2, 0, 1, 1, 1 /)
     coefficients(403) = 3.37384794614262800000e-001_ark
     expansionIndices(404, :) = (/ 1, 2, 1, 0, 0, 2 /)
     coefficients(404) = 1.68692397307131400000e-001_ark
     expansionIndices(405, :) = (/ 1, 2, 1, 0, 2, 0 /)
     coefficients(405) = -1.68692397307131400000e-001_ark
     expansionIndices(406, :) = (/ 1, 2, 0, 2, 1, 0 /)
     coefficients(406) = -5.11357115487325430000e-001_ark
     expansionIndices(407, :) = (/ 1, 2, 1, 1, 0, 1 /)
     coefficients(407) = -1.02271423097465090000e+000_ark
     expansionIndices(408, :) = (/ 1, 2, 2, 0, 1, 0 /)
     coefficients(408) = 5.11357115487325430000e-001_ark
     expansionIndices(409, :) = (/ 1, 2, 1, 2, 0, 0 /)
     coefficients(409) = -7.35239007903031760000e-001_ark
     expansionIndices(410, :) = (/ 1, 2, 3, 0, 0, 0 /)
     coefficients(410) = 2.45079669301010580000e-001_ark
     expansionIndices(411, :) = (/ 1, 3, 0, 0, 0, 2 /)
     coefficients(411) = 9.25959359296515520000e-001_ark
     expansionIndices(412, :) = (/ 1, 3, 0, 0, 2, 0 /)
     coefficients(412) = 9.25959359296515520000e-001_ark
     expansionIndices(413, :) = (/ 1, 3, 0, 1, 0, 1 /)
     coefficients(413) = -2.66734278176963800000e+000_ark
     expansionIndices(414, :) = (/ 1, 3, 1, 0, 1, 0 /)
     coefficients(414) = -2.66734278176963800000e+000_ark
     expansionIndices(415, :) = (/ 1, 3, 0, 2, 0, 0 /)
     coefficients(415) = -3.15487875028764670000e+000_ark
     expansionIndices(416, :) = (/ 1, 3, 2, 0, 0, 0 /)
     coefficients(416) = -3.15487875028764670000e+000_ark
     expansionIndices(417, :) = (/ 1, 5, 0, 0, 0, 0 /)
     coefficients(417) = -2.94162970489296180000e+001_ark
     expansionIndices(418, :) = (/ 2, 0, 0, 0, 0, 4 /)
     coefficients(418) = -1.01342951641511960000e-003_ark
     expansionIndices(419, :) = (/ 2, 0, 0, 0, 2, 2 /)
     coefficients(419) = -2.02685903282978030000e-003_ark
     expansionIndices(420, :) = (/ 2, 0, 0, 0, 4, 0 /)
     coefficients(420) = -1.01342951641511960000e-003_ark
     expansionIndices(421, :) = (/ 2, 0, 0, 2, 0, 2 /)
     coefficients(421) = -1.64609807677887280000e-002_ark
     expansionIndices(422, :) = (/ 2, 0, 0, 2, 2, 0 /)
     coefficients(422) = -1.64609807677887280000e-002_ark
     expansionIndices(423, :) = (/ 2, 0, 2, 0, 0, 2 /)
     coefficients(423) = -1.64609807677887280000e-002_ark
     expansionIndices(424, :) = (/ 2, 0, 2, 0, 2, 0 /)
     coefficients(424) = -1.64609807677887280000e-002_ark
     expansionIndices(425, :) = (/ 2, 0, 0, 2, 0, 2 /)
     coefficients(425) = 3.77511276019006780000e-003_ark
     expansionIndices(426, :) = (/ 2, 0, 0, 2, 2, 0 /)
     coefficients(426) = -3.77511276019006780000e-003_ark
     expansionIndices(427, :) = (/ 2, 0, 1, 1, 1, 1 /)
     coefficients(427) = 1.51004510407602710000e-002_ark
     expansionIndices(428, :) = (/ 2, 0, 2, 0, 0, 2 /)
     coefficients(428) = -3.77511276019006780000e-003_ark
     expansionIndices(429, :) = (/ 2, 0, 2, 0, 2, 0 /)
     coefficients(429) = 3.77511276019006780000e-003_ark
     expansionIndices(430, :) = (/ 2, 0, 0, 3, 0, 1 /)
     coefficients(430) = -6.61410213487317790000e-003_ark
     expansionIndices(431, :) = (/ 2, 0, 1, 2, 1, 0 /)
     coefficients(431) = -6.61410213487317790000e-003_ark
     expansionIndices(432, :) = (/ 2, 0, 2, 1, 0, 1 /)
     coefficients(432) = -6.61410213487317790000e-003_ark
     expansionIndices(433, :) = (/ 2, 0, 3, 0, 1, 0 /)
     coefficients(433) = -6.61410213487317790000e-003_ark
     expansionIndices(434, :) = (/ 2, 0, 0, 4, 0, 0 /)
     coefficients(434) = 5.62841546187971710000e-003_ark
     expansionIndices(435, :) = (/ 2, 0, 2, 2, 0, 0 /)
     coefficients(435) = 1.12568309237568880000e-002_ark
     expansionIndices(436, :) = (/ 2, 0, 4, 0, 0, 0 /)
     coefficients(436) = 5.62841546187971710000e-003_ark
     expansionIndices(437, :) = (/ 2, 1, 0, 0, 1, 2 /)
     coefficients(437) = -4.16729820549751330000e-002_ark
     expansionIndices(438, :) = (/ 2, 1, 0, 0, 3, 0 /)
     coefficients(438) = 1.38909940183250430000e-002_ark
     expansionIndices(439, :) = (/ 2, 1, 0, 1, 1, 1 /)
     coefficients(439) = -2.66066314114886990000e-002_ark
     expansionIndices(440, :) = (/ 2, 1, 1, 0, 0, 2 /)
     coefficients(440) = -1.33033157057443500000e-002_ark
     expansionIndices(441, :) = (/ 2, 1, 1, 0, 2, 0 /)
     coefficients(441) = 1.33033157057443500000e-002_ark
     expansionIndices(442, :) = (/ 2, 2, 0, 1, 0, 1 /)
     coefficients(442) = 5.77914943306844900000e-001_ark
     expansionIndices(443, :) = (/ 2, 2, 1, 0, 1, 0 /)
     coefficients(443) = 5.77914943306844900000e-001_ark
     expansionIndices(444, :) = (/ 2, 2, 0, 2, 0, 0 /)
     coefficients(444) = 4.73169528191837830000e-001_ark
     expansionIndices(445, :) = (/ 2, 2, 2, 0, 0, 0 /)
     coefficients(445) = 4.73169528191837830000e-001_ark
     expansionIndices(446, :) = (/ 2, 4, 0, 0, 0, 0 /)
     coefficients(446) = 1.44405094692787510000e+001_ark
     expansionIndices(447, :) = (/ 3, 1, 0, 0, 0, 2 /)
     coefficients(447) = 1.37665732835014440000e-002_ark
     expansionIndices(448, :) = (/ 3, 1, 0, 0, 2, 0 /)
     coefficients(448) = 1.37665732835014440000e-002_ark
     expansionIndices(449, :) = (/ 3, 1, 0, 2, 0, 0 /)
     coefficients(449) = 4.06946609777004560000e-002_ark
     expansionIndices(450, :) = (/ 3, 1, 2, 0, 0, 0 /)
     coefficients(450) = 4.06946609777004560000e-002_ark
     expansionIndices(451, :) = (/ 3, 3, 0, 0, 0, 0 /)
     coefficients(451) = -2.77137099472502780000e+000_ark
     expansionIndices(452, :) = (/ 4, 0, 0, 0, 0, 2 /)
     coefficients(452) = -5.80331540857252550000e-003_ark
     expansionIndices(453, :) = (/ 4, 0, 0, 0, 2, 0 /)
     coefficients(453) = -5.80331540857252550000e-003_ark
     expansionIndices(454, :) = (/ 4, 0, 0, 1, 0, 1 /)
     coefficients(454) = -4.56929793416770750000e-003_ark
     expansionIndices(455, :) = (/ 4, 0, 1, 0, 1, 0 /)
     coefficients(455) = -4.56929793416770750000e-003_ark
     expansionIndices(456, :) = (/ 5, 1, 0, 0, 0, 0 /)
     coefficients(456) = 8.93608045246492430000e-002_ark
     expansionIndices(457, :) = (/ 6, 0, 0, 0, 0, 0 /)
     coefficients(457) = -3.16182867599765260000e-003_ark
     expansionIndices(458, :) = (/ 0, 0, 0, 0, 1, 6 /)
     coefficients(458) = 1.36452019599825750000e-002_ark
     expansionIndices(459, :) = (/ 0, 0, 0, 0, 3, 4 /)
     coefficients(459) = 2.27420032666319880000e-002_ark
     expansionIndices(460, :) = (/ 0, 0, 0, 0, 5, 2 /)
     coefficients(460) = 4.54840065332696130000e-003_ark
     expansionIndices(461, :) = (/ 0, 0, 0, 0, 7, 0 /)
     coefficients(461) = -4.54840065332696130000e-003_ark
     expansionIndices(462, :) = (/ 0, 0, 0, 1, 1, 5 /)
     coefficients(462) = 9.42417406579784850000e-003_ark
     expansionIndices(463, :) = (/ 0, 0, 0, 1, 3, 3 /)
     coefficients(463) = 1.88483481316035350000e-002_ark
     expansionIndices(464, :) = (/ 0, 0, 0, 1, 5, 1 /)
     coefficients(464) = 9.42417406579784850000e-003_ark
     expansionIndices(465, :) = (/ 0, 0, 1, 0, 0, 6 /)
     coefficients(465) = 4.71208703289761890000e-003_ark
     expansionIndices(466, :) = (/ 0, 0, 1, 0, 2, 4 /)
     coefficients(466) = 4.71208703289761890000e-003_ark
     expansionIndices(467, :) = (/ 0, 0, 1, 0, 4, 2 /)
     coefficients(467) = -4.71208703289761890000e-003_ark
     expansionIndices(468, :) = (/ 0, 0, 1, 0, 6, 0 /)
     coefficients(468) = -4.71208703289761890000e-003_ark
     expansionIndices(469, :) = (/ 0, 0, 0, 1, 1, 5 /)
     coefficients(469) = 2.30331044517275520000e-003_ark
     expansionIndices(470, :) = (/ 0, 0, 0, 1, 5, 1 /)
     coefficients(470) = -2.30331044517275520000e-003_ark
     expansionIndices(471, :) = (/ 0, 0, 1, 0, 0, 6 /)
     coefficients(471) = -5.75827611293188810000e-004_ark
     expansionIndices(472, :) = (/ 0, 0, 1, 0, 2, 4 /)
     coefficients(472) = 2.87913805646644850000e-003_ark
     expansionIndices(473, :) = (/ 0, 0, 1, 0, 4, 2 /)
     coefficients(473) = 2.87913805646644850000e-003_ark
     expansionIndices(474, :) = (/ 0, 0, 1, 0, 6, 0 /)
     coefficients(474) = -5.75827611293188810000e-004_ark
     expansionIndices(475, :) = (/ 0, 0, 0, 2, 1, 4 /)
     coefficients(475) = 1.57838201696235900000e-002_ark
     expansionIndices(476, :) = (/ 0, 0, 0, 2, 3, 2 /)
     coefficients(476) = 1.05225467797559360000e-002_ark
     expansionIndices(477, :) = (/ 0, 0, 0, 2, 5, 0 /)
     coefficients(477) = -5.26127338987590530000e-003_ark
     expansionIndices(478, :) = (/ 0, 0, 2, 0, 1, 4 /)
     coefficients(478) = 1.57838201696235900000e-002_ark
     expansionIndices(479, :) = (/ 0, 0, 2, 0, 3, 2 /)
     coefficients(479) = 1.05225467797559360000e-002_ark
     expansionIndices(480, :) = (/ 0, 0, 2, 0, 5, 0 /)
     coefficients(480) = -5.26127338987590530000e-003_ark
     expansionIndices(481, :) = (/ 0, 0, 0, 2, 1, 4 /)
     coefficients(481) = -2.00306239586244390000e-003_ark
     expansionIndices(482, :) = (/ 0, 0, 0, 2, 3, 2 /)
     coefficients(482) = -4.00612479172488770000e-003_ark
     expansionIndices(483, :) = (/ 0, 0, 0, 2, 5, 0 /)
     coefficients(483) = -2.00306239586244390000e-003_ark
     expansionIndices(484, :) = (/ 0, 0, 1, 1, 0, 5 /)
     coefficients(484) = -4.00612479172488770000e-003_ark
     expansionIndices(485, :) = (/ 0, 0, 1, 1, 2, 3 /)
     coefficients(485) = -8.01224958345291700000e-003_ark
     expansionIndices(486, :) = (/ 0, 0, 1, 1, 4, 1 /)
     coefficients(486) = -4.00612479172488770000e-003_ark
     expansionIndices(487, :) = (/ 0, 0, 2, 0, 1, 4 /)
     coefficients(487) = 2.00306239586244390000e-003_ark
     expansionIndices(488, :) = (/ 0, 0, 2, 0, 3, 2 /)
     coefficients(488) = 4.00612479172488770000e-003_ark
     expansionIndices(489, :) = (/ 0, 0, 2, 0, 5, 0 /)
     coefficients(489) = 2.00306239586244390000e-003_ark
     expansionIndices(490, :) = (/ 0, 1, 0, 0, 0, 6 /)
     coefficients(490) = 1.28166697127569150000e-002_ark
     expansionIndices(491, :) = (/ 0, 1, 0, 0, 2, 4 /)
     coefficients(491) = 3.84500091382768970000e-002_ark
     expansionIndices(492, :) = (/ 0, 1, 0, 0, 4, 2 /)
     coefficients(492) = 3.84500091382768970000e-002_ark
     expansionIndices(493, :) = (/ 0, 1, 0, 0, 6, 0 /)
     coefficients(493) = 1.28166697127569150000e-002_ark
     expansionIndices(494, :) = (/ 0, 1, 0, 0, 0, 6 /)
     coefficients(494) = -2.33212529256537380000e-003_ark
     expansionIndices(495, :) = (/ 0, 1, 0, 0, 2, 4 /)
     coefficients(495) = 3.49818793884806020000e-002_ark
     expansionIndices(496, :) = (/ 0, 1, 0, 0, 4, 2 /)
     coefficients(496) = -3.49818793884806020000e-002_ark
     expansionIndices(497, :) = (/ 0, 1, 0, 0, 6, 0 /)
     coefficients(497) = 2.33212529256537380000e-003_ark
     expansionIndices(498, :) = (/ 0, 1, 0, 1, 0, 5 /)
     coefficients(498) = 2.71008526668017580000e-002_ark
     expansionIndices(499, :) = (/ 0, 1, 0, 1, 2, 3 /)
     coefficients(499) = 5.42017053336247700000e-002_ark
     expansionIndices(500, :) = (/ 0, 1, 0, 1, 4, 1 /)
     coefficients(500) = 2.71008526668017580000e-002_ark
     expansionIndices(501, :) = (/ 0, 1, 1, 0, 1, 4 /)
     coefficients(501) = 2.71008526668017580000e-002_ark
     expansionIndices(502, :) = (/ 0, 1, 1, 0, 3, 2 /)
     coefficients(502) = 5.42017053336247700000e-002_ark
     expansionIndices(503, :) = (/ 0, 1, 1, 0, 5, 0 /)
     coefficients(503) = 2.71008526668017580000e-002_ark
     expansionIndices(504, :) = (/ 0, 1, 0, 2, 0, 4 /)
     coefficients(504) = 3.25357997530843350000e-002_ark
     expansionIndices(505, :) = (/ 0, 1, 0, 2, 2, 2 /)
     coefficients(505) = 6.50715995061686700000e-002_ark
     expansionIndices(506, :) = (/ 0, 1, 0, 2, 4, 0 /)
     coefficients(506) = 3.25357997530843350000e-002_ark
     expansionIndices(507, :) = (/ 0, 1, 2, 0, 0, 4 /)
     coefficients(507) = 3.25357997530843350000e-002_ark
     expansionIndices(508, :) = (/ 0, 1, 2, 0, 2, 2 /)
     coefficients(508) = 6.50715995061686700000e-002_ark
     expansionIndices(509, :) = (/ 0, 1, 2, 0, 4, 0 /)
     coefficients(509) = 3.25357997530843350000e-002_ark
     expansionIndices(510, :) = (/ 0, 1, 0, 2, 0, 4 /)
     coefficients(510) = 5.76328492483375860000e-003_ark
     expansionIndices(511, :) = (/ 0, 1, 0, 2, 4, 0 /)
     coefficients(511) = -5.76328492483375860000e-003_ark
     expansionIndices(512, :) = (/ 0, 1, 1, 1, 1, 3 /)
     coefficients(512) = 2.30531396993327760000e-002_ark
     expansionIndices(513, :) = (/ 0, 1, 1, 1, 3, 1 /)
     coefficients(513) = 2.30531396993327760000e-002_ark
     expansionIndices(514, :) = (/ 0, 1, 2, 0, 0, 4 /)
     coefficients(514) = -5.76328492483375860000e-003_ark
     expansionIndices(515, :) = (/ 0, 1, 2, 0, 4, 0 /)
     coefficients(515) = 5.76328492483375860000e-003_ark
     expansionIndices(516, :) = (/ 0, 2, 0, 0, 1, 4 /)
     coefficients(516) = -6.82445132937132160000e-001_ark
     expansionIndices(517, :) = (/ 0, 2, 0, 0, 3, 2 /)
     coefficients(517) = -4.54963421958487440000e-001_ark
     expansionIndices(518, :) = (/ 0, 2, 0, 0, 5, 0 /)
     coefficients(518) = 2.27481710979149130000e-001_ark
     expansionIndices(519, :) = (/ 0, 2, 0, 1, 1, 3 /)
     coefficients(519) = -1.78240237873709790000e-001_ark
     expansionIndices(520, :) = (/ 0, 2, 0, 1, 3, 1 /)
     coefficients(520) = -1.78240237873709790000e-001_ark
     expansionIndices(521, :) = (/ 0, 2, 1, 0, 0, 4 /)
     coefficients(521) = -8.91201189368548970000e-002_ark
     expansionIndices(522, :) = (/ 0, 2, 1, 0, 4, 0 /)
     coefficients(522) = 8.91201189368548970000e-002_ark
     expansionIndices(523, :) = (/ 0, 2, 0, 1, 1, 3 /)
     coefficients(523) = -1.92130155157054590000e-001_ark
     expansionIndices(524, :) = (/ 0, 2, 0, 1, 3, 1 /)
     coefficients(524) = 1.92130155157054590000e-001_ark
     expansionIndices(525, :) = (/ 0, 2, 1, 0, 0, 4 /)
     coefficients(525) = 4.80325387892636470000e-002_ark
     expansionIndices(526, :) = (/ 0, 2, 1, 0, 2, 2 /)
     coefficients(526) = -2.88195232735581900000e-001_ark
     expansionIndices(527, :) = (/ 0, 2, 1, 0, 4, 0 /)
     coefficients(527) = 4.80325387892636470000e-002_ark
     expansionIndices(528, :) = (/ 0, 2, 0, 3, 1, 1 /)
     coefficients(528) = 1.65617219594152430000e-001_ark
     expansionIndices(529, :) = (/ 0, 2, 1, 2, 0, 2 /)
     coefficients(529) = 8.28086097970949510000e-002_ark
     expansionIndices(530, :) = (/ 0, 2, 1, 2, 2, 0 /)
     coefficients(530) = -8.28086097970949510000e-002_ark
     expansionIndices(531, :) = (/ 0, 2, 2, 1, 1, 1 /)
     coefficients(531) = 1.65617219594152430000e-001_ark
     expansionIndices(532, :) = (/ 0, 2, 3, 0, 0, 2 /)
     coefficients(532) = 8.28086097970949510000e-002_ark
     expansionIndices(533, :) = (/ 0, 2, 3, 0, 2, 0 /)
     coefficients(533) = -8.28086097970949510000e-002_ark
     expansionIndices(534, :) = (/ 0, 2, 1, 2, 0, 2 /)
     coefficients(534) = -1.90484221071887700000e-001_ark
     expansionIndices(535, :) = (/ 0, 2, 1, 2, 2, 0 /)
     coefficients(535) = -1.90484221071887700000e-001_ark
     expansionIndices(536, :) = (/ 0, 2, 3, 0, 0, 2 /)
     coefficients(536) = 6.34947403572871110000e-002_ark
     expansionIndices(537, :) = (/ 0, 2, 3, 0, 2, 0 /)
     coefficients(537) = 6.34947403572871110000e-002_ark
     expansionIndices(538, :) = (/ 0, 2, 0, 4, 1, 0 /)
     coefficients(538) = -2.61357158396862630000e-001_ark
     expansionIndices(539, :) = (/ 0, 2, 1, 3, 0, 1 /)
     coefficients(539) = -5.22714316793725260000e-001_ark
     expansionIndices(540, :) = (/ 0, 2, 3, 1, 0, 1 /)
     coefficients(540) = -5.22714316793725260000e-001_ark
     expansionIndices(541, :) = (/ 0, 2, 4, 0, 1, 0 /)
     coefficients(541) = 2.61357158396862630000e-001_ark
     expansionIndices(542, :) = (/ 0, 2, 1, 4, 0, 0 /)
     coefficients(542) = -3.68394934673988370000e-001_ark
     expansionIndices(543, :) = (/ 0, 2, 3, 2, 0, 0 /)
     coefficients(543) = -2.45596623116207810000e-001_ark
     expansionIndices(544, :) = (/ 0, 2, 5, 0, 0, 0 /)
     coefficients(544) = 1.22798311558052850000e-001_ark
     expansionIndices(545, :) = (/ 0, 3, 0, 0, 0, 4 /)
     coefficients(545) = 3.30778062908504690000e+000_ark
     expansionIndices(546, :) = (/ 0, 3, 0, 0, 2, 2 /)
     coefficients(546) = 6.61556125817268730000e+000_ark
     expansionIndices(547, :) = (/ 0, 3, 0, 0, 4, 0 /)
     coefficients(547) = 3.30778062908504690000e+000_ark
     expansionIndices(548, :) = (/ 0, 3, 0, 1, 0, 3 /)
     coefficients(548) = 1.78915795993339130000e+000_ark
     expansionIndices(549, :) = (/ 0, 3, 0, 1, 2, 1 /)
     coefficients(549) = 1.78915795993339130000e+000_ark
     expansionIndices(550, :) = (/ 0, 3, 1, 0, 1, 2 /)
     coefficients(550) = 1.78915795993339130000e+000_ark
     expansionIndices(551, :) = (/ 0, 3, 1, 0, 3, 0 /)
     coefficients(551) = 1.78915795993339130000e+000_ark
     expansionIndices(552, :) = (/ 0, 3, 0, 2, 0, 2 /)
     coefficients(552) = 1.47405937086355430000e+000_ark
     expansionIndices(553, :) = (/ 0, 3, 0, 2, 2, 0 /)
     coefficients(553) = 1.47405937086355430000e+000_ark
     expansionIndices(554, :) = (/ 0, 3, 2, 0, 0, 2 /)
     coefficients(554) = 1.47405937086355430000e+000_ark
     expansionIndices(555, :) = (/ 0, 3, 2, 0, 2, 0 /)
     coefficients(555) = 1.47405937086355430000e+000_ark
     expansionIndices(556, :) = (/ 0, 3, 0, 2, 0, 2 /)
     coefficients(556) = 2.30433583977541350000e-001_ark
     expansionIndices(557, :) = (/ 0, 3, 0, 2, 2, 0 /)
     coefficients(557) = -2.30433583977541350000e-001_ark
     expansionIndices(558, :) = (/ 0, 3, 1, 1, 1, 1 /)
     coefficients(558) = 9.21734335910165400000e-001_ark
     expansionIndices(559, :) = (/ 0, 3, 2, 0, 0, 2 /)
     coefficients(559) = -2.30433583977541350000e-001_ark
     expansionIndices(560, :) = (/ 0, 3, 2, 0, 2, 0 /)
     coefficients(560) = 2.30433583977541350000e-001_ark
     expansionIndices(561, :) = (/ 0, 3, 0, 3, 0, 1 /)
     coefficients(561) = -1.86815665388909040000e+000_ark
     expansionIndices(562, :) = (/ 0, 3, 1, 2, 1, 0 /)
     coefficients(562) = -1.86815665388909040000e+000_ark
     expansionIndices(563, :) = (/ 0, 3, 2, 1, 0, 1 /)
     coefficients(563) = -1.86815665388909040000e+000_ark
     expansionIndices(564, :) = (/ 0, 3, 3, 0, 1, 0 /)
     coefficients(564) = -1.86815665388909040000e+000_ark
     expansionIndices(565, :) = (/ 0, 3, 0, 4, 0, 0 /)
     coefficients(565) = -1.56552127007293350000e+000_ark
     expansionIndices(566, :) = (/ 0, 3, 2, 2, 0, 0 /)
     coefficients(566) = -3.13104254014709450000e+000_ark
     expansionIndices(567, :) = (/ 0, 3, 4, 0, 0, 0 /)
     coefficients(567) = -1.56552127007293350000e+000_ark
     expansionIndices(568, :) = (/ 0, 4, 0, 0, 1, 2 /)
     coefficients(568) = -7.10904902740490030000e+001_ark
     expansionIndices(569, :) = (/ 0, 4, 0, 0, 3, 0 /)
     coefficients(569) = 2.36968300913458770000e+001_ark
     expansionIndices(570, :) = (/ 0, 4, 0, 1, 1, 1 /)
     coefficients(570) = -5.92137068509013440000e+000_ark
     expansionIndices(571, :) = (/ 0, 4, 1, 0, 0, 2 /)
     coefficients(571) = -2.96068534254506720000e+000_ark
     expansionIndices(572, :) = (/ 0, 4, 1, 0, 2, 0 /)
     coefficients(572) = 2.96068534254506720000e+000_ark
     expansionIndices(573, :) = (/ 0, 4, 0, 2, 1, 0 /)
     coefficients(573) = -5.25076052139603360000e+000_ark
     expansionIndices(574, :) = (/ 0, 4, 1, 1, 0, 1 /)
     coefficients(574) = -1.05015210427920670000e+001_ark
     expansionIndices(575, :) = (/ 0, 4, 2, 0, 1, 0 /)
     coefficients(575) = 5.25076052139603360000e+000_ark
     expansionIndices(576, :) = (/ 0, 4, 1, 2, 0, 0 /)
     coefficients(576) = -1.25510243978137210000e+001_ark
     expansionIndices(577, :) = (/ 0, 4, 3, 0, 0, 0 /)
     coefficients(577) = 4.18367479927057140000e+000_ark
     expansionIndices(578, :) = (/ 0, 5, 0, 0, 0, 2 /)
     coefficients(578) = 1.30398540330004860000e+002_ark
     expansionIndices(579, :) = (/ 0, 5, 0, 0, 2, 0 /)
     coefficients(579) = 1.30398540330004860000e+002_ark
     expansionIndices(580, :) = (/ 0, 5, 0, 1, 0, 1 /)
     coefficients(580) = -1.62540568667730020000e+001_ark
     expansionIndices(581, :) = (/ 0, 5, 1, 0, 1, 0 /)
     coefficients(581) = -1.62540568667730020000e+001_ark
     expansionIndices(582, :) = (/ 0, 5, 0, 2, 0, 0 /)
     coefficients(582) = -2.66948633047460720000e+001_ark
     expansionIndices(583, :) = (/ 0, 5, 2, 0, 0, 0 /)
     coefficients(583) = -2.66948633047460720000e+001_ark
     expansionIndices(584, :) = (/ 0, 7, 0, 0, 0, 0 /)
     coefficients(584) = -2.50636283694460820000e+001_ark
     expansionIndices(585, :) = (/ 1, 0, 0, 0, 0, 6 /)
     coefficients(585) = -6.00039232477717980000e-004_ark
     expansionIndices(586, :) = (/ 1, 0, 0, 0, 2, 4 /)
     coefficients(586) = -1.80011769743344190000e-003_ark
     expansionIndices(587, :) = (/ 1, 0, 0, 0, 4, 2 /)
     coefficients(587) = -1.80011769743344190000e-003_ark
     expansionIndices(588, :) = (/ 1, 0, 0, 0, 6, 0 /)
     coefficients(588) = -6.00039232477717980000e-004_ark
     expansionIndices(589, :) = (/ 1, 0, 0, 1, 0, 5 /)
     coefficients(589) = 6.07065202707867020000e-003_ark
     expansionIndices(590, :) = (/ 1, 0, 0, 1, 2, 3 /)
     coefficients(590) = 1.21413040541621000000e-002_ark
     expansionIndices(591, :) = (/ 1, 0, 0, 1, 4, 1 /)
     coefficients(591) = 6.07065202707867020000e-003_ark
     expansionIndices(592, :) = (/ 1, 0, 1, 0, 1, 4 /)
     coefficients(592) = 6.07065202707867020000e-003_ark
     expansionIndices(593, :) = (/ 1, 0, 1, 0, 3, 2 /)
     coefficients(593) = 1.21413040541621000000e-002_ark
     expansionIndices(594, :) = (/ 1, 0, 1, 0, 5, 0 /)
     coefficients(594) = 6.07065202707867020000e-003_ark
     expansionIndices(595, :) = (/ 1, 0, 0, 2, 0, 4 /)
     coefficients(595) = -3.82558297548990250000e-003_ark
     expansionIndices(596, :) = (/ 1, 0, 0, 2, 2, 2 /)
     coefficients(596) = -7.65116595097980500000e-003_ark
     expansionIndices(597, :) = (/ 1, 0, 0, 2, 4, 0 /)
     coefficients(597) = -3.82558297548990250000e-003_ark
     expansionIndices(598, :) = (/ 1, 0, 2, 0, 0, 4 /)
     coefficients(598) = -3.82558297548990250000e-003_ark
     expansionIndices(599, :) = (/ 1, 0, 2, 0, 2, 2 /)
     coefficients(599) = -7.65116595097980500000e-003_ark
     expansionIndices(600, :) = (/ 1, 0, 2, 0, 4, 0 /)
     coefficients(600) = -3.82558297548990250000e-003_ark
     expansionIndices(601, :) = (/ 1, 0, 0, 2, 0, 4 /)
     coefficients(601) = -2.87693193944079090000e-003_ark
     expansionIndices(602, :) = (/ 1, 0, 0, 2, 4, 0 /)
     coefficients(602) = 2.87693193944079090000e-003_ark
     expansionIndices(603, :) = (/ 1, 0, 1, 1, 1, 3 /)
     coefficients(603) = -1.15077277577620360000e-002_ark
     expansionIndices(604, :) = (/ 1, 0, 1, 1, 3, 1 /)
     coefficients(604) = -1.15077277577620360000e-002_ark
     expansionIndices(605, :) = (/ 1, 0, 2, 0, 0, 4 /)
     coefficients(605) = 2.87693193944079090000e-003_ark
     expansionIndices(606, :) = (/ 1, 0, 2, 0, 4, 0 /)
     coefficients(606) = -2.87693193944079090000e-003_ark
     expansionIndices(607, :) = (/ 1, 0, 0, 4, 0, 2 /)
     coefficients(607) = -1.26359215447782740000e-002_ark
     expansionIndices(608, :) = (/ 1, 0, 0, 4, 2, 0 /)
     coefficients(608) = -1.26359215447782740000e-002_ark
     expansionIndices(609, :) = (/ 1, 0, 2, 2, 0, 2 /)
     coefficients(609) = -2.52718430895565480000e-002_ark
     expansionIndices(610, :) = (/ 1, 0, 2, 2, 2, 0 /)
     coefficients(610) = -2.52718430895565480000e-002_ark
     expansionIndices(611, :) = (/ 1, 0, 4, 0, 0, 2 /)
     coefficients(611) = -1.26359215447782740000e-002_ark
     expansionIndices(612, :) = (/ 1, 0, 4, 0, 2, 0 /)
     coefficients(612) = -1.26359215447782740000e-002_ark
     expansionIndices(613, :) = (/ 1, 1, 0, 0, 1, 4 /)
     coefficients(613) = 2.20809589806972660000e-001_ark
     expansionIndices(614, :) = (/ 1, 1, 0, 0, 3, 2 /)
     coefficients(614) = 1.47206393204744600000e-001_ark
     expansionIndices(615, :) = (/ 1, 1, 0, 0, 5, 0 /)
     coefficients(615) = -7.36031966023434490000e-002_ark
     expansionIndices(616, :) = (/ 1, 2, 0, 0, 0, 4 /)
     coefficients(616) = -7.95897552165634510000e-001_ark
     expansionIndices(617, :) = (/ 1, 2, 0, 0, 2, 2 /)
     coefficients(617) = -1.59179510433090890000e+000_ark
     expansionIndices(618, :) = (/ 1, 2, 0, 0, 4, 0 /)
     coefficients(618) = -7.95897552165634510000e-001_ark
     expansionIndices(619, :) = (/ 1, 2, 0, 2, 0, 2 /)
     coefficients(619) = 2.99615881604856970000e-001_ark
     expansionIndices(620, :) = (/ 1, 2, 0, 2, 2, 0 /)
     coefficients(620) = 2.99615881604856970000e-001_ark
     expansionIndices(621, :) = (/ 1, 2, 2, 0, 0, 2 /)
     coefficients(621) = 2.99615881604856970000e-001_ark
     expansionIndices(622, :) = (/ 1, 2, 2, 0, 2, 0 /)
     coefficients(622) = 2.99615881604856970000e-001_ark
     expansionIndices(623, :) = (/ 1, 2, 0, 3, 0, 1 /)
     coefficients(623) = 1.21737867028125520000e+000_ark
     expansionIndices(624, :) = (/ 1, 2, 1, 2, 1, 0 /)
     coefficients(624) = 1.21737867028125520000e+000_ark
     expansionIndices(625, :) = (/ 1, 2, 2, 1, 0, 1 /)
     coefficients(625) = 1.21737867028125520000e+000_ark
     expansionIndices(626, :) = (/ 1, 2, 3, 0, 1, 0 /)
     coefficients(626) = 1.21737867028125520000e+000_ark
     expansionIndices(627, :) = (/ 1, 2, 0, 4, 0, 0 /)
     coefficients(627) = 7.54245263811392670000e-001_ark
     expansionIndices(628, :) = (/ 1, 2, 2, 2, 0, 0 /)
     coefficients(628) = 1.50849052762244410000e+000_ark
     expansionIndices(629, :) = (/ 1, 2, 4, 0, 0, 0 /)
     coefficients(629) = 7.54245263811392670000e-001_ark
     expansionIndices(630, :) = (/ 1, 3, 0, 0, 1, 2 /)
     coefficients(630) = 1.21775291528821090000e+001_ark
     expansionIndices(631, :) = (/ 1, 3, 0, 0, 3, 0 /)
     coefficients(631) = -4.05917638429338720000e+000_ark
     expansionIndices(632, :) = (/ 1, 3, 0, 1, 1, 1 /)
     coefficients(632) = -4.14685090448236700000e-001_ark
     expansionIndices(633, :) = (/ 1, 3, 1, 0, 0, 2 /)
     coefficients(633) = -2.07342545224118350000e-001_ark
     expansionIndices(634, :) = (/ 1, 3, 1, 0, 2, 0 /)
     coefficients(634) = 2.07342545224118350000e-001_ark
     expansionIndices(635, :) = (/ 1, 3, 0, 2, 1, 0 /)
     coefficients(635) = 4.18910740005011560000e+000_ark
     expansionIndices(636, :) = (/ 1, 3, 1, 1, 0, 1 /)
     coefficients(636) = 8.37821480010023120000e+000_ark
     expansionIndices(637, :) = (/ 1, 3, 2, 0, 1, 0 /)
     coefficients(637) = -4.18910740005011560000e+000_ark
     expansionIndices(638, :) = (/ 1, 3, 1, 2, 0, 0 /)
     coefficients(638) = 9.56925482371995880000e+000_ark
     expansionIndices(639, :) = (/ 1, 3, 3, 0, 0, 0 /)
     coefficients(639) = -3.18975160790614300000e+000_ark
     expansionIndices(640, :) = (/ 1, 4, 0, 0, 0, 2 /)
     coefficients(640) = -3.04401376563397470000e+001_ark
     expansionIndices(641, :) = (/ 1, 4, 0, 0, 2, 0 /)
     coefficients(641) = -3.04401376563397470000e+001_ark
     expansionIndices(642, :) = (/ 1, 4, 0, 1, 0, 1 /)
     coefficients(642) = 2.33800624908670120000e+001_ark
     expansionIndices(643, :) = (/ 1, 4, 1, 0, 1, 0 /)
     coefficients(643) = 2.33800624908670120000e+001_ark
     expansionIndices(644, :) = (/ 1, 4, 0, 2, 0, 0 /)
     coefficients(644) = 3.06779558454284360000e+001_ark
     expansionIndices(645, :) = (/ 1, 4, 2, 0, 0, 0 /)
     coefficients(645) = 3.06779558454284360000e+001_ark
     expansionIndices(646, :) = (/ 1, 6, 0, 0, 0, 0 /)
     coefficients(646) = 8.34020464049663520000e+000_ark
     expansionIndices(647, :) = (/ 2, 0, 0, 0, 1, 4 /)
     coefficients(647) = 4.19923894990918960000e-003_ark
     expansionIndices(648, :) = (/ 2, 0, 0, 0, 3, 2 /)
     coefficients(648) = 2.79949263327524970000e-003_ark
     expansionIndices(649, :) = (/ 2, 0, 0, 0, 5, 0 /)
     coefficients(649) = -1.39974631663704310000e-003_ark
     expansionIndices(650, :) = (/ 2, 0, 0, 2, 1, 2 /)
     coefficients(650) = -1.09285161185470740000e-002_ark
     expansionIndices(651, :) = (/ 2, 0, 0, 2, 3, 0 /)
     coefficients(651) = -1.09285161185470740000e-002_ark
     expansionIndices(652, :) = (/ 2, 0, 1, 1, 0, 3 /)
     coefficients(652) = -2.18570322370892010000e-002_ark
     expansionIndices(653, :) = (/ 2, 0, 1, 1, 2, 1 /)
     coefficients(653) = -2.18570322370892010000e-002_ark
     expansionIndices(654, :) = (/ 2, 0, 2, 0, 1, 2 /)
     coefficients(654) = 1.09285161185470740000e-002_ark
     expansionIndices(655, :) = (/ 2, 0, 2, 0, 3, 0 /)
     coefficients(655) = 1.09285161185470740000e-002_ark
     expansionIndices(656, :) = (/ 2, 1, 0, 0, 0, 4 /)
     coefficients(656) = 1.96829766369592540000e-002_ark
     expansionIndices(657, :) = (/ 2, 1, 0, 0, 2, 2 /)
     coefficients(657) = 3.93659532739095990000e-002_ark
     expansionIndices(658, :) = (/ 2, 1, 0, 0, 4, 0 /)
     coefficients(658) = 1.96829766369592540000e-002_ark
     expansionIndices(659, :) = (/ 2, 2, 0, 0, 1, 2 /)
     coefficients(659) = -5.76380990884725390000e-001_ark
     expansionIndices(660, :) = (/ 2, 2, 0, 0, 3, 0 /)
     coefficients(660) = 1.92126996961449640000e-001_ark
     expansionIndices(661, :) = (/ 2, 2, 0, 1, 1, 1 /)
     coefficients(661) = 5.39698983843832190000e-001_ark
     expansionIndices(662, :) = (/ 2, 2, 1, 0, 0, 2 /)
     coefficients(662) = 2.69849491921839820000e-001_ark
     expansionIndices(663, :) = (/ 2, 2, 1, 0, 2, 0 /)
     coefficients(663) = -2.69849491921839820000e-001_ark
     expansionIndices(664, :) = (/ 2, 2, 0, 2, 1, 0 /)
     coefficients(664) = -9.22587492189204660000e-001_ark
     expansionIndices(665, :) = (/ 2, 2, 1, 1, 0, 1 /)
     coefficients(665) = -1.84517498437893110000e+000_ark
     expansionIndices(666, :) = (/ 2, 2, 2, 0, 1, 0 /)
     coefficients(666) = 9.22587492189204660000e-001_ark
     expansionIndices(667, :) = (/ 2, 2, 1, 2, 0, 0 /)
     coefficients(667) = -2.09988838181697450000e+000_ark
     expansionIndices(668, :) = (/ 2, 2, 3, 0, 0, 0 /)
     coefficients(668) = 6.99962793938534330000e-001_ark
     expansionIndices(669, :) = (/ 2, 3, 0, 0, 0, 2 /)
     coefficients(669) = 1.91635921216184820000e+000_ark
     expansionIndices(670, :) = (/ 2, 3, 0, 0, 2, 0 /)
     coefficients(670) = 1.91635921216184820000e+000_ark
     expansionIndices(671, :) = (/ 2, 3, 0, 1, 0, 1 /)
     coefficients(671) = -8.34244675670895930000e+000_ark
     expansionIndices(672, :) = (/ 2, 3, 1, 0, 1, 0 /)
     coefficients(672) = -8.34244675670895930000e+000_ark
     expansionIndices(673, :) = (/ 2, 3, 0, 2, 0, 0 /)
     coefficients(673) = -1.07290870077245690000e+001_ark
     expansionIndices(674, :) = (/ 2, 3, 2, 0, 0, 0 /)
     coefficients(674) = -1.07290870077245690000e+001_ark
     expansionIndices(675, :) = (/ 2, 5, 0, 0, 0, 0 /)
     coefficients(675) = -8.26514866825268940000e+001_ark
     expansionIndices(676, :) = (/ 3, 0, 0, 0, 0, 4 /)
     coefficients(676) = -5.96489871112884300000e-003_ark
     expansionIndices(677, :) = (/ 3, 0, 0, 0, 2, 2 /)
     coefficients(677) = -1.19297974222623610000e-002_ark
     expansionIndices(678, :) = (/ 3, 0, 0, 0, 4, 0 /)
     coefficients(678) = -5.96489871112884300000e-003_ark
     expansionIndices(679, :) = (/ 3, 0, 0, 2, 0, 2 /)
     coefficients(679) = -2.96337067234089600000e-002_ark
     expansionIndices(680, :) = (/ 3, 0, 0, 2, 2, 0 /)
     coefficients(680) = -2.96337067234089600000e-002_ark
     expansionIndices(681, :) = (/ 3, 0, 2, 0, 0, 2 /)
     coefficients(681) = -2.96337067234089600000e-002_ark
     expansionIndices(682, :) = (/ 3, 0, 2, 0, 2, 0 /)
     coefficients(682) = -2.96337067234089600000e-002_ark
     expansionIndices(683, :) = (/ 3, 2, 0, 1, 0, 1 /)
     coefficients(683) = 1.19391543874391280000e+000_ark
     expansionIndices(684, :) = (/ 3, 2, 1, 0, 1, 0 /)
     coefficients(684) = 1.19391543874391280000e+000_ark
     expansionIndices(685, :) = (/ 3, 2, 0, 2, 0, 0 /)
     coefficients(685) = 1.23851414213596690000e+000_ark
     expansionIndices(686, :) = (/ 3, 2, 2, 0, 0, 0 /)
     coefficients(686) = 1.23851414213596690000e+000_ark
     expansionIndices(687, :) = (/ 3, 4, 0, 0, 0, 0 /)
     coefficients(687) = 2.97495251237368910000e+001_ark
     expansionIndices(688, :) = (/ 4, 1, 0, 0, 0, 2 /)
     coefficients(688) = 3.98605322132695310000e-002_ark
     expansionIndices(689, :) = (/ 4, 1, 0, 0, 2, 0 /)
     coefficients(689) = 3.98605322132695310000e-002_ark
     expansionIndices(690, :) = (/ 4, 3, 0, 0, 0, 0 /)
     coefficients(690) = -4.92083189346669860000e+000_ark
     expansionIndices(691, :) = (/ 5, 0, 0, 0, 0, 2 /)
     coefficients(691) = -1.23998546496072290000e-002_ark
     expansionIndices(692, :) = (/ 5, 0, 0, 0, 2, 0 /)
     coefficients(692) = -1.23998546496072290000e-002_ark
     expansionIndices(693, :) = (/ 5, 2, 0, 0, 0, 0 /)
     coefficients(693) = 1.40465057790493010000e-001_ark
     expansionIndices(694, :) = (/ 6, 1, 0, 0, 0, 0 /)
     coefficients(694) = 9.06219972550804660000e-002_ark
     expansionIndices(695, :) = (/ 7, 0, 0, 0, 0, 0 /)
     coefficients(695) = -2.61062652617386840000e-003_ark
     expansionIndices(696, :) = (/ 0, 0, 0, 0, 0, 8 /)
     coefficients(696) = 5.43852981543006470000e-003_ark
     expansionIndices(697, :) = (/ 0, 0, 0, 0, 2, 6 /)
     coefficients(697) = 2.17541192617244360000e-002_ark
     expansionIndices(698, :) = (/ 0, 0, 0, 0, 4, 4 /)
     coefficients(698) = 3.26311788925762110000e-002_ark
     expansionIndices(699, :) = (/ 0, 0, 0, 0, 6, 2 /)
     coefficients(699) = 2.17541192617244360000e-002_ark
     expansionIndices(700, :) = (/ 0, 0, 0, 0, 8, 0 /)
     coefficients(700) = 5.43852981543006470000e-003_ark
     expansionIndices(701, :) = (/ 0, 0, 0, 0, 0, 8 /)
     coefficients(701) = -7.02219236942119810000e-004_ark
     expansionIndices(702, :) = (/ 0, 0, 0, 0, 2, 6 /)
     coefficients(702) = 9.83106931718967690000e-003_ark
     expansionIndices(703, :) = (/ 0, 0, 0, 0, 6, 2 /)
     coefficients(703) = -9.83106931718967690000e-003_ark
     expansionIndices(704, :) = (/ 0, 0, 0, 0, 8, 0 /)
     coefficients(704) = 7.02219236942119810000e-004_ark
     expansionIndices(705, :) = (/ 0, 0, 0, 1, 0, 7 /)
     coefficients(705) = 4.79200607723380670000e-003_ark
     expansionIndices(706, :) = (/ 0, 0, 0, 1, 2, 5 /)
     coefficients(706) = 1.43760182316975160000e-002_ark
     expansionIndices(707, :) = (/ 0, 0, 0, 1, 4, 3 /)
     coefficients(707) = 1.43760182316975160000e-002_ark
     expansionIndices(708, :) = (/ 0, 0, 0, 1, 6, 1 /)
     coefficients(708) = 4.79200607723380670000e-003_ark
     expansionIndices(709, :) = (/ 0, 0, 1, 0, 1, 6 /)
     coefficients(709) = 4.79200607723380670000e-003_ark
     expansionIndices(710, :) = (/ 0, 0, 1, 0, 3, 4 /)
     coefficients(710) = 1.43760182316975160000e-002_ark
     expansionIndices(711, :) = (/ 0, 0, 1, 0, 5, 2 /)
     coefficients(711) = 1.43760182316975160000e-002_ark
     expansionIndices(712, :) = (/ 0, 0, 1, 0, 7, 0 /)
     coefficients(712) = 4.79200607723380670000e-003_ark
     expansionIndices(713, :) = (/ 0, 0, 0, 2, 0, 6 /)
     coefficients(713) = 1.76970393050034030000e-003_ark
     expansionIndices(714, :) = (/ 0, 0, 0, 2, 2, 4 /)
     coefficients(714) = 5.30911179149932190000e-003_ark
     expansionIndices(715, :) = (/ 0, 0, 0, 2, 4, 2 /)
     coefficients(715) = 5.30911179149932190000e-003_ark
     expansionIndices(716, :) = (/ 0, 0, 0, 2, 6, 0 /)
     coefficients(716) = 1.76970393050034030000e-003_ark
     expansionIndices(717, :) = (/ 0, 0, 2, 0, 0, 6 /)
     coefficients(717) = 1.76970393050034030000e-003_ark
     expansionIndices(718, :) = (/ 0, 0, 2, 0, 2, 4 /)
     coefficients(718) = 5.30911179149932190000e-003_ark
     expansionIndices(719, :) = (/ 0, 0, 2, 0, 4, 2 /)
     coefficients(719) = 5.30911179149932190000e-003_ark
     expansionIndices(720, :) = (/ 0, 0, 2, 0, 6, 0 /)
     coefficients(720) = 1.76970393050034030000e-003_ark
     expansionIndices(721, :) = (/ 0, 0, 0, 8, 0, 0 /)
     coefficients(721) = -3.52956993646548750000e-004_ark
     expansionIndices(722, :) = (/ 0, 0, 2, 6, 0, 0 /)
     coefficients(722) = -1.41182797458646630000e-003_ark
     expansionIndices(723, :) = (/ 0, 0, 4, 4, 0, 0 /)
     coefficients(723) = -2.11774196187902160000e-003_ark
     expansionIndices(724, :) = (/ 0, 0, 6, 2, 0, 0 /)
     coefficients(724) = -1.41182797458646630000e-003_ark
     expansionIndices(725, :) = (/ 0, 0, 8, 0, 0, 0 /)
     coefficients(725) = -3.52956993646548750000e-004_ark
     expansionIndices(726, :) = (/ 0, 1, 0, 0, 1, 6 /)
     coefficients(726) = -9.96749730819395250000e-002_ark
     expansionIndices(727, :) = (/ 0, 1, 0, 0, 3, 4 /)
     coefficients(727) = -1.66124955136524700000e-001_ark
     expansionIndices(728, :) = (/ 0, 1, 0, 0, 5, 2 /)
     coefficients(728) = -3.32249910273090530000e-002_ark
     expansionIndices(729, :) = (/ 0, 1, 0, 0, 7, 0 /)
     coefficients(729) = 3.32249910273090530000e-002_ark
     expansionIndices(730, :) = (/ 0, 1, 0, 1, 1, 5 /)
     coefficients(730) = -3.44247818330508940000e-002_ark
     expansionIndices(731, :) = (/ 0, 1, 0, 1, 5, 1 /)
     coefficients(731) = 3.44247818330508940000e-002_ark
     expansionIndices(732, :) = (/ 0, 1, 1, 0, 0, 6 /)
     coefficients(732) = 8.60619545826272340000e-003_ark
     expansionIndices(733, :) = (/ 0, 1, 1, 0, 2, 4 /)
     coefficients(733) = -4.30309772913211530000e-002_ark
     expansionIndices(734, :) = (/ 0, 1, 1, 0, 4, 2 /)
     coefficients(734) = -4.30309772913211530000e-002_ark
     expansionIndices(735, :) = (/ 0, 1, 1, 0, 6, 0 /)
     coefficients(735) = 8.60619545826272340000e-003_ark
     expansionIndices(736, :) = (/ 0, 1, 0, 2, 1, 4 /)
     coefficients(736) = 2.00860887553197170000e-001_ark
     expansionIndices(737, :) = (/ 0, 1, 0, 2, 3, 2 /)
     coefficients(737) = 1.33907258368885590000e-001_ark
     expansionIndices(738, :) = (/ 0, 1, 0, 2, 5, 0 /)
     coefficients(738) = -6.69536291844165520000e-002_ark
     expansionIndices(739, :) = (/ 0, 1, 2, 0, 1, 4 /)
     coefficients(739) = 2.00860887553197170000e-001_ark
     expansionIndices(740, :) = (/ 0, 1, 2, 0, 3, 2 /)
     coefficients(740) = 1.33907258368885590000e-001_ark
     expansionIndices(741, :) = (/ 0, 1, 2, 0, 5, 0 /)
     coefficients(741) = -6.69536291844165520000e-002_ark
     expansionIndices(742, :) = (/ 0, 2, 0, 0, 0, 6 /)
     coefficients(742) = 2.45063315990961290000e-001_ark
     expansionIndices(743, :) = (/ 0, 2, 0, 0, 2, 4 /)
     coefficients(743) = 7.35189947972883840000e-001_ark
     expansionIndices(744, :) = (/ 0, 2, 0, 0, 4, 2 /)
     coefficients(744) = 7.35189947972883840000e-001_ark
     expansionIndices(745, :) = (/ 0, 2, 0, 0, 6, 0 /)
     coefficients(745) = 2.45063315990961290000e-001_ark
     expansionIndices(746, :) = (/ 0, 2, 0, 0, 0, 6 /)
     coefficients(746) = -1.77603585804302460000e-002_ark
     expansionIndices(747, :) = (/ 0, 2, 0, 0, 2, 4 /)
     coefficients(747) = 2.66405378706491780000e-001_ark
     expansionIndices(748, :) = (/ 0, 2, 0, 0, 4, 2 /)
     coefficients(748) = -2.66405378706491780000e-001_ark
     expansionIndices(749, :) = (/ 0, 2, 0, 0, 6, 0 /)
     coefficients(749) = 1.77603585804302460000e-002_ark
     expansionIndices(750, :) = (/ 0, 2, 0, 1, 0, 5 /)
     coefficients(750) = 2.35881636079534830000e-001_ark
     expansionIndices(751, :) = (/ 0, 2, 0, 1, 2, 3 /)
     coefficients(751) = 4.71763272159069660000e-001_ark
     expansionIndices(752, :) = (/ 0, 2, 0, 1, 4, 1 /)
     coefficients(752) = 2.35881636079534830000e-001_ark
     expansionIndices(753, :) = (/ 0, 2, 1, 0, 1, 4 /)
     coefficients(753) = 2.35881636079534830000e-001_ark
     expansionIndices(754, :) = (/ 0, 2, 1, 0, 3, 2 /)
     coefficients(754) = 4.71763272159069660000e-001_ark
     expansionIndices(755, :) = (/ 0, 2, 1, 0, 5, 0 /)
     coefficients(755) = 2.35881636079534830000e-001_ark
     expansionIndices(756, :) = (/ 0, 2, 0, 2, 0, 4 /)
     coefficients(756) = -1.75765790682822280000e-001_ark
     expansionIndices(757, :) = (/ 0, 2, 0, 2, 2, 2 /)
     coefficients(757) = -3.51531581365644550000e-001_ark
     expansionIndices(758, :) = (/ 0, 2, 0, 2, 4, 0 /)
     coefficients(758) = -1.75765790682822280000e-001_ark
     expansionIndices(759, :) = (/ 0, 2, 2, 0, 0, 4 /)
     coefficients(759) = -1.75765790682822280000e-001_ark
     expansionIndices(760, :) = (/ 0, 2, 2, 0, 2, 2 /)
     coefficients(760) = -3.51531581365644550000e-001_ark
     expansionIndices(761, :) = (/ 0, 2, 2, 0, 4, 0 /)
     coefficients(761) = -1.75765790682822280000e-001_ark
     expansionIndices(762, :) = (/ 0, 2, 0, 5, 0, 1 /)
     coefficients(762) = 2.61769348782680710000e-001_ark
     expansionIndices(763, :) = (/ 0, 2, 1, 4, 1, 0 /)
     coefficients(763) = 2.61769348782680710000e-001_ark
     expansionIndices(764, :) = (/ 0, 2, 2, 3, 0, 1 /)
     coefficients(764) = 5.23538697565361420000e-001_ark
     expansionIndices(765, :) = (/ 0, 2, 3, 2, 1, 0 /)
     coefficients(765) = 5.23538697565361420000e-001_ark
     expansionIndices(766, :) = (/ 0, 2, 4, 1, 0, 1 /)
     coefficients(766) = 2.61769348782680710000e-001_ark
     expansionIndices(767, :) = (/ 0, 2, 5, 0, 1, 0 /)
     coefficients(767) = 2.61769348782680710000e-001_ark
     expansionIndices(768, :) = (/ 0, 3, 0, 0, 1, 4 /)
     coefficients(768) = -3.63425308596538600000e-001_ark
     expansionIndices(769, :) = (/ 0, 3, 0, 0, 3, 2 /)
     coefficients(769) = -2.42283539064514140000e-001_ark
     expansionIndices(770, :) = (/ 0, 3, 0, 0, 5, 0 /)
     coefficients(770) = 1.21141769532198920000e-001_ark
     expansionIndices(771, :) = (/ 0, 3, 0, 1, 1, 3 /)
     coefficients(771) = -8.10927018212042160000e-001_ark
     expansionIndices(772, :) = (/ 0, 3, 0, 1, 3, 1 /)
     coefficients(772) = -8.10927018212042160000e-001_ark
     expansionIndices(773, :) = (/ 0, 3, 1, 0, 0, 4 /)
     coefficients(773) = -4.05463509105972450000e-001_ark
     expansionIndices(774, :) = (/ 0, 3, 1, 0, 4, 0 /)
     coefficients(774) = 4.05463509105972450000e-001_ark
     expansionIndices(775, :) = (/ 0, 3, 0, 2, 1, 2 /)
     coefficients(775) = 5.29302350772943790000e+000_ark
     expansionIndices(776, :) = (/ 0, 3, 0, 2, 3, 0 /)
     coefficients(776) = -1.76434116924286390000e+000_ark
     expansionIndices(777, :) = (/ 0, 3, 2, 0, 1, 2 /)
     coefficients(777) = 5.29302350772943790000e+000_ark
     expansionIndices(778, :) = (/ 0, 3, 2, 0, 3, 0 /)
     coefficients(778) = -1.76434116924286390000e+000_ark
     expansionIndices(779, :) = (/ 0, 3, 0, 4, 1, 0 /)
     coefficients(779) = 1.71918054392190190000e+000_ark
     expansionIndices(780, :) = (/ 0, 3, 1, 3, 0, 1 /)
     coefficients(780) = 3.43836108784421630000e+000_ark
     expansionIndices(781, :) = (/ 0, 3, 3, 1, 0, 1 /)
     coefficients(781) = 3.43836108784421630000e+000_ark
     expansionIndices(782, :) = (/ 0, 3, 4, 0, 1, 0 /)
     coefficients(782) = -1.71918054392190190000e+000_ark
     expansionIndices(783, :) = (/ 0, 3, 1, 4, 0, 0 /)
     coefficients(783) = 2.61927872002551390000e+000_ark
     expansionIndices(784, :) = (/ 0, 3, 3, 2, 0, 0 /)
     coefficients(784) = 1.74618581335146010000e+000_ark
     expansionIndices(785, :) = (/ 0, 3, 5, 0, 0, 0 /)
     coefficients(785) = -8.73092906675311030000e-001_ark
     expansionIndices(786, :) = (/ 0, 4, 0, 0, 0, 4 /)
     coefficients(786) = -9.65855223957085580000e+000_ark
     expansionIndices(787, :) = (/ 0, 4, 0, 0, 2, 2 /)
     coefficients(787) = -1.93171044791417120000e+001_ark
     expansionIndices(788, :) = (/ 0, 4, 0, 0, 4, 0 /)
     coefficients(788) = -9.65855223957085580000e+000_ark
     expansionIndices(789, :) = (/ 0, 4, 0, 1, 0, 3 /)
     coefficients(789) = -6.00026473221909920000e+000_ark
     expansionIndices(790, :) = (/ 0, 4, 0, 1, 2, 1 /)
     coefficients(790) = -6.00026473221909920000e+000_ark
     expansionIndices(791, :) = (/ 0, 4, 1, 0, 1, 2 /)
     coefficients(791) = -6.00026473221909920000e+000_ark
     expansionIndices(792, :) = (/ 0, 4, 1, 0, 3, 0 /)
     coefficients(792) = -6.00026473221909920000e+000_ark
     expansionIndices(793, :) = (/ 0, 4, 0, 2, 0, 2 /)
     coefficients(793) = -8.82776507100306240000e+000_ark
     expansionIndices(794, :) = (/ 0, 4, 0, 2, 2, 0 /)
     coefficients(794) = -8.82776507100306240000e+000_ark
     expansionIndices(795, :) = (/ 0, 4, 2, 0, 0, 2 /)
     coefficients(795) = -8.82776507100306240000e+000_ark
     expansionIndices(796, :) = (/ 0, 4, 2, 0, 2, 0 /)
     coefficients(796) = -8.82776507100306240000e+000_ark
     expansionIndices(797, :) = (/ 0, 4, 0, 3, 0, 1 /)
     coefficients(797) = 8.58118877419019200000e+000_ark
     expansionIndices(798, :) = (/ 0, 4, 1, 2, 1, 0 /)
     coefficients(798) = 8.58118877419019200000e+000_ark
     expansionIndices(799, :) = (/ 0, 4, 2, 1, 0, 1 /)
     coefficients(799) = 8.58118877419019200000e+000_ark
     expansionIndices(800, :) = (/ 0, 4, 3, 0, 1, 0 /)
     coefficients(800) = 8.58118877419019200000e+000_ark
     expansionIndices(801, :) = (/ 0, 4, 0, 4, 0, 0 /)
     coefficients(801) = 8.64863397339527220000e+000_ark
     expansionIndices(802, :) = (/ 0, 4, 2, 2, 0, 0 /)
     coefficients(802) = 1.72972679467905440000e+001_ark
     expansionIndices(803, :) = (/ 0, 4, 4, 0, 0, 0 /)
     coefficients(803) = 8.64863397339527220000e+000_ark
     expansionIndices(804, :) = (/ 0, 5, 0, 0, 1, 2 /)
     coefficients(804) = 2.92594958533547360000e+002_ark
     expansionIndices(805, :) = (/ 0, 5, 0, 0, 3, 0 /)
     coefficients(805) = -9.75316528445157760000e+001_ark
     expansionIndices(806, :) = (/ 0, 5, 0, 1, 1, 1 /)
     coefficients(806) = 9.69069979207997710000e+000_ark
     expansionIndices(807, :) = (/ 0, 5, 1, 0, 0, 2 /)
     coefficients(807) = 4.84534989604073910000e+000_ark
     expansionIndices(808, :) = (/ 0, 5, 1, 0, 2, 0 /)
     coefficients(808) = -4.84534989604073910000e+000_ark
     expansionIndices(809, :) = (/ 0, 5, 0, 2, 1, 0 /)
     coefficients(809) = 1.39429242505006620000e+001_ark
     expansionIndices(810, :) = (/ 0, 5, 1, 1, 0, 1 /)
     coefficients(810) = 2.78858485009970000000e+001_ark
     expansionIndices(811, :) = (/ 0, 5, 2, 0, 1, 0 /)
     coefficients(811) = -1.39429242505006620000e+001_ark
     expansionIndices(812, :) = (/ 0, 5, 1, 2, 0, 0 /)
     coefficients(812) = 3.95255215299596120000e+001_ark
     expansionIndices(813, :) = (/ 0, 5, 3, 0, 0, 0 /)
     coefficients(813) = -1.31751738433198700000e+001_ark
     expansionIndices(814, :) = (/ 0, 6, 0, 0, 0, 2 /)
     coefficients(814) = -5.13283570883237530000e+002_ark
     expansionIndices(815, :) = (/ 0, 6, 0, 0, 2, 0 /)
     coefficients(815) = -5.13283570883237530000e+002_ark
     expansionIndices(816, :) = (/ 0, 6, 0, 1, 0, 1 /)
     coefficients(816) = 5.83343210358891890000e+001_ark
     expansionIndices(817, :) = (/ 0, 6, 1, 0, 1, 0 /)
     coefficients(817) = 5.83343210358891890000e+001_ark
     expansionIndices(818, :) = (/ 0, 6, 0, 2, 0, 0 /)
     coefficients(818) = 2.47049246919140340000e+001_ark
     expansionIndices(819, :) = (/ 0, 6, 2, 0, 0, 0 /)
     coefficients(819) = 2.47049246919140340000e+001_ark
     expansionIndices(820, :) = (/ 0, 8, 0, 0, 0, 0 /)
     coefficients(820) = 2.38169246969607970000e+001_ark
     expansionIndices(821, :) = (/ 1, 0, 0, 0, 1, 6 /)
     coefficients(821) = -5.13014887438405490000e-003_ark
     expansionIndices(822, :) = (/ 1, 0, 0, 0, 3, 4 /)
     coefficients(822) = -8.55024812397130530000e-003_ark
     expansionIndices(823, :) = (/ 1, 0, 0, 0, 5, 2 /)
     coefficients(823) = -1.71004962479447280000e-003_ark
     expansionIndices(824, :) = (/ 1, 0, 0, 0, 7, 0 /)
     coefficients(824) = 1.71004962479447280000e-003_ark
     expansionIndices(825, :) = (/ 1, 0, 0, 3, 1, 3 /)
     coefficients(825) = -3.25260441910132230000e-002_ark
     expansionIndices(826, :) = (/ 1, 0, 0, 3, 3, 1 /)
     coefficients(826) = -3.25260441910132230000e-002_ark
     expansionIndices(827, :) = (/ 1, 0, 1, 2, 0, 4 /)
     coefficients(827) = -1.62630220955043560000e-002_ark
     expansionIndices(828, :) = (/ 1, 0, 1, 2, 4, 0 /)
     coefficients(828) = 1.62630220955043560000e-002_ark
     expansionIndices(829, :) = (/ 1, 0, 2, 1, 1, 3 /)
     coefficients(829) = -3.25260441910132230000e-002_ark
     expansionIndices(830, :) = (/ 1, 0, 2, 1, 3, 1 /)
     coefficients(830) = -3.25260441910132230000e-002_ark
     expansionIndices(831, :) = (/ 1, 0, 3, 0, 0, 4 /)
     coefficients(831) = -1.62630220955043560000e-002_ark
     expansionIndices(832, :) = (/ 1, 0, 3, 0, 4, 0 /)
     coefficients(832) = 1.62630220955043560000e-002_ark
     expansionIndices(833, :) = (/ 1, 1, 0, 0, 0, 6 /)
     coefficients(833) = 4.12618650198663530000e-002_ark
     expansionIndices(834, :) = (/ 1, 1, 0, 0, 2, 4 /)
     coefficients(834) = 1.23785595059618860000e-001_ark
     expansionIndices(835, :) = (/ 1, 1, 0, 0, 4, 2 /)
     coefficients(835) = 1.23785595059618860000e-001_ark
     expansionIndices(836, :) = (/ 1, 1, 0, 0, 6, 0 /)
     coefficients(836) = 4.12618650198663530000e-002_ark
     expansionIndices(837, :) = (/ 1, 2, 0, 0, 1, 4 /)
     coefficients(837) = -1.21058362144396670000e+000_ark
     expansionIndices(838, :) = (/ 1, 2, 0, 0, 3, 2 /)
     coefficients(838) = -8.07055747630019420000e-001_ark
     expansionIndices(839, :) = (/ 1, 2, 0, 0, 5, 0 /)
     coefficients(839) = 4.03527873814841950000e-001_ark
     expansionIndices(840, :) = (/ 1, 2, 0, 2, 1, 2 /)
     coefficients(840) = 4.56105250181033040000e-001_ark
     expansionIndices(841, :) = (/ 1, 2, 0, 2, 3, 0 /)
     coefficients(841) = 4.56105250181033040000e-001_ark
     expansionIndices(842, :) = (/ 1, 2, 1, 1, 0, 3 /)
     coefficients(842) = 9.12210500361859690000e-001_ark
     expansionIndices(843, :) = (/ 1, 2, 1, 1, 2, 1 /)
     coefficients(843) = 9.12210500361859690000e-001_ark
     expansionIndices(844, :) = (/ 1, 2, 2, 0, 1, 2 /)
     coefficients(844) = -4.56105250181033040000e-001_ark
     expansionIndices(845, :) = (/ 1, 2, 2, 0, 3, 0 /)
     coefficients(845) = -4.56105250181033040000e-001_ark
     expansionIndices(846, :) = (/ 1, 3, 0, 0, 0, 4 /)
     coefficients(846) = 1.94571823955353840000e+000_ark
     expansionIndices(847, :) = (/ 1, 3, 0, 0, 2, 2 /)
     coefficients(847) = 3.89143647910860220000e+000_ark
     expansionIndices(848, :) = (/ 1, 3, 0, 0, 4, 0 /)
     coefficients(848) = 1.94571823955353840000e+000_ark
     expansionIndices(849, :) = (/ 1, 3, 0, 1, 0, 3 /)
     coefficients(849) = 1.74345538050975430000e+000_ark
     expansionIndices(850, :) = (/ 1, 3, 0, 1, 2, 1 /)
     coefficients(850) = 1.74345538050975430000e+000_ark
     expansionIndices(851, :) = (/ 1, 3, 1, 0, 1, 2 /)
     coefficients(851) = 1.74345538050975430000e+000_ark
     expansionIndices(852, :) = (/ 1, 3, 1, 0, 3, 0 /)
     coefficients(852) = 1.74345538050975430000e+000_ark
     expansionIndices(853, :) = (/ 1, 3, 0, 3, 0, 1 /)
     coefficients(853) = -8.32661687320075930000e+000_ark
     expansionIndices(854, :) = (/ 1, 3, 1, 2, 1, 0 /)
     coefficients(854) = -8.32661687320075930000e+000_ark
     expansionIndices(855, :) = (/ 1, 3, 2, 1, 0, 1 /)
     coefficients(855) = -8.32661687320075930000e+000_ark
     expansionIndices(856, :) = (/ 1, 3, 3, 0, 1, 0 /)
     coefficients(856) = -8.32661687320075930000e+000_ark
     expansionIndices(857, :) = (/ 1, 3, 0, 4, 0, 0 /)
     coefficients(857) = -5.68097284773724990000e+000_ark
     expansionIndices(858, :) = (/ 1, 3, 2, 2, 0, 0 /)
     coefficients(858) = -1.13619456954789550000e+001_ark
     expansionIndices(859, :) = (/ 1, 3, 4, 0, 0, 0 /)
     coefficients(859) = -5.68097284773724990000e+000_ark
     expansionIndices(860, :) = (/ 1, 4, 0, 0, 1, 2 /)
     coefficients(860) = -5.24605586835362130000e+001_ark
     expansionIndices(861, :) = (/ 1, 4, 0, 0, 3, 0 /)
     coefficients(861) = 1.74868528945092730000e+001_ark
     expansionIndices(862, :) = (/ 1, 4, 0, 2, 1, 0 /)
     coefficients(862) = -1.55759993795930940000e+001_ark
     expansionIndices(863, :) = (/ 1, 4, 1, 1, 0, 1 /)
     coefficients(863) = -3.11519987591861880000e+001_ark
     expansionIndices(864, :) = (/ 1, 4, 2, 0, 1, 0 /)
     coefficients(864) = 1.55759993795930940000e+001_ark
     expansionIndices(865, :) = (/ 1, 4, 1, 2, 0, 0 /)
     coefficients(865) = -4.53109947361196530000e+001_ark
     expansionIndices(866, :) = (/ 1, 4, 3, 0, 0, 0 /)
     coefficients(866) = 1.51036649120374700000e+001_ark
     expansionIndices(867, :) = (/ 1, 5, 0, 0, 0, 2 /)
     coefficients(867) = 1.27490821203075300000e+002_ark
     expansionIndices(868, :) = (/ 1, 5, 0, 0, 2, 0 /)
     coefficients(868) = 1.27490821203075300000e+002_ark
     expansionIndices(869, :) = (/ 1, 5, 0, 1, 0, 1 /)
     coefficients(869) = -1.01983956628515860000e+002_ark
     expansionIndices(870, :) = (/ 1, 5, 1, 0, 1, 0 /)
     coefficients(870) = -1.01983956628515860000e+002_ark
     expansionIndices(871, :) = (/ 1, 5, 0, 2, 0, 0 /)
     coefficients(871) = -1.30764525885980080000e+002_ark
     expansionIndices(872, :) = (/ 1, 5, 2, 0, 0, 0 /)
     coefficients(872) = -1.30764525885980080000e+002_ark
     expansionIndices(873, :) = (/ 1, 7, 0, 0, 0, 0 /)
     coefficients(873) = 2.18815695005428270000e+002_ark
     expansionIndices(874, :) = (/ 2, 0, 0, 0, 0, 6 /)
     coefficients(874) = 2.20335855856287200000e-003_ark
     expansionIndices(875, :) = (/ 2, 0, 0, 0, 2, 4 /)
     coefficients(875) = 6.61007567568861610000e-003_ark
     expansionIndices(876, :) = (/ 2, 0, 0, 0, 4, 2 /)
     coefficients(876) = 6.61007567568861610000e-003_ark
     expansionIndices(877, :) = (/ 2, 0, 0, 0, 6, 0 /)
     coefficients(877) = 2.20335855856287200000e-003_ark
     expansionIndices(878, :) = (/ 2, 1, 0, 0, 1, 4 /)
     coefficients(878) = 8.07203702574379730000e-002_ark
     expansionIndices(879, :) = (/ 2, 1, 0, 0, 3, 2 /)
     coefficients(879) = 5.38135801716725440000e-002_ark
     expansionIndices(880, :) = (/ 2, 1, 0, 0, 5, 0 /)
     coefficients(880) = -2.69067900858250860000e-002_ark
     expansionIndices(881, :) = (/ 2, 2, 0, 0, 0, 4 /)
     coefficients(881) = -6.72320501628424490000e-002_ark
     expansionIndices(882, :) = (/ 2, 2, 0, 0, 2, 2 /)
     coefficients(882) = -1.34464100325684900000e-001_ark
     expansionIndices(883, :) = (/ 2, 2, 0, 0, 4, 0 /)
     coefficients(883) = -6.72320501628424490000e-002_ark
     expansionIndices(884, :) = (/ 2, 2, 0, 2, 0, 2 /)
     coefficients(884) = 5.82129781512261910000e-001_ark
     expansionIndices(885, :) = (/ 2, 2, 0, 2, 2, 0 /)
     coefficients(885) = 5.82129781512261910000e-001_ark
     expansionIndices(886, :) = (/ 2, 2, 2, 0, 0, 2 /)
     coefficients(886) = 5.82129781512261910000e-001_ark
     expansionIndices(887, :) = (/ 2, 2, 2, 0, 2, 0 /)
     coefficients(887) = 5.82129781512261910000e-001_ark
     expansionIndices(888, :) = (/ 2, 2, 0, 3, 0, 1 /)
     coefficients(888) = 8.43937480509686200000e-001_ark
     expansionIndices(889, :) = (/ 2, 2, 1, 2, 1, 0 /)
     coefficients(889) = 8.43937480509686200000e-001_ark
     expansionIndices(890, :) = (/ 2, 2, 2, 1, 0, 1 /)
     coefficients(890) = 8.43937480509686200000e-001_ark
     expansionIndices(891, :) = (/ 2, 2, 3, 0, 1, 0 /)
     coefficients(891) = 8.43937480509686200000e-001_ark
     expansionIndices(892, :) = (/ 2, 3, 0, 0, 1, 2 /)
     coefficients(892) = 4.92445289109022790000e+000_ark
     expansionIndices(893, :) = (/ 2, 3, 0, 0, 3, 0 /)
     coefficients(893) = -1.64148429703016880000e+000_ark
     expansionIndices(894, :) = (/ 2, 3, 0, 1, 1, 1 /)
     coefficients(894) = -3.46099522304571570000e+000_ark
     expansionIndices(895, :) = (/ 2, 3, 1, 0, 0, 2 /)
     coefficients(895) = -1.73049761152370560000e+000_ark
     expansionIndices(896, :) = (/ 2, 3, 1, 0, 2, 0 /)
     coefficients(896) = 1.73049761152370560000e+000_ark
     expansionIndices(897, :) = (/ 2, 3, 0, 2, 1, 0 /)
     coefficients(897) = 3.74250448385992130000e+000_ark
     expansionIndices(898, :) = (/ 2, 3, 1, 1, 0, 1 /)
     coefficients(898) = 7.48500896771617620000e+000_ark
     expansionIndices(899, :) = (/ 2, 3, 2, 0, 1, 0 /)
     coefficients(899) = -3.74250448385992130000e+000_ark
     expansionIndices(900, :) = (/ 2, 3, 1, 2, 0, 0 /)
     coefficients(900) = 1.40860390148465800000e+001_ark
     expansionIndices(901, :) = (/ 2, 3, 3, 0, 0, 0 /)
     coefficients(901) = -4.69534633828245870000e+000_ark
     expansionIndices(902, :) = (/ 2, 4, 0, 0, 0, 2 /)
     coefficients(902) = -1.01913105645860020000e+001_ark
     expansionIndices(903, :) = (/ 2, 4, 0, 0, 2, 0 /)
     coefficients(903) = -1.01913105645860020000e+001_ark
     expansionIndices(904, :) = (/ 2, 4, 0, 1, 0, 1 /)
     coefficients(904) = 4.14927922003984480000e+001_ark
     expansionIndices(905, :) = (/ 2, 4, 1, 0, 1, 0 /)
     coefficients(905) = 4.14927922003984480000e+001_ark
     expansionIndices(906, :) = (/ 2, 4, 0, 2, 0, 0 /)
     coefficients(906) = 6.18419285914721680000e+001_ark
     expansionIndices(907, :) = (/ 2, 4, 2, 0, 0, 0 /)
     coefficients(907) = 6.18419285914721680000e+001_ark
     expansionIndices(908, :) = (/ 2, 6, 0, 0, 0, 0 /)
     coefficients(908) = 2.43967059737864760000e+002_ark
     expansionIndices(909, :) = (/ 3, 0, 0, 1, 1, 3 /)
     coefficients(909) = -4.36934042314379430000e-002_ark
     expansionIndices(910, :) = (/ 3, 0, 0, 1, 3, 1 /)
     coefficients(910) = -4.36934042314379430000e-002_ark
     expansionIndices(911, :) = (/ 3, 0, 1, 0, 0, 4 /)
     coefficients(911) = -2.18467021157163520000e-002_ark
     expansionIndices(912, :) = (/ 3, 0, 1, 0, 4, 0 /)
     coefficients(912) = 2.18467021157163520000e-002_ark
     expansionIndices(913, :) = (/ 3, 1, 0, 0, 0, 4 /)
     coefficients(913) = -8.47535869737384570000e-002_ark
     expansionIndices(914, :) = (/ 3, 1, 0, 0, 2, 2 /)
     coefficients(914) = -1.69507173947543360000e-001_ark
     expansionIndices(915, :) = (/ 3, 1, 0, 0, 4, 0 /)
     coefficients(915) = -8.47535869737384570000e-002_ark
     expansionIndices(916, :) = (/ 3, 1, 0, 2, 0, 2 /)
     coefficients(916) = -2.99249643944080100000e-001_ark
     expansionIndices(917, :) = (/ 3, 1, 0, 2, 2, 0 /)
     coefficients(917) = -2.99249643944080100000e-001_ark
     expansionIndices(918, :) = (/ 3, 1, 2, 0, 0, 2 /)
     coefficients(918) = -2.99249643944080100000e-001_ark
     expansionIndices(919, :) = (/ 3, 1, 2, 0, 2, 0 /)
     coefficients(919) = -2.99249643944080100000e-001_ark
     expansionIndices(920, :) = (/ 3, 3, 0, 0, 0, 2 /)
     coefficients(920) = -3.25101230880359360000e-001_ark
     expansionIndices(921, :) = (/ 3, 3, 0, 0, 2, 0 /)
     coefficients(921) = -3.25101230880359360000e-001_ark
     expansionIndices(922, :) = (/ 3, 3, 0, 1, 0, 1 /)
     coefficients(922) = -7.38342976908927810000e+000_ark
     expansionIndices(923, :) = (/ 3, 3, 1, 0, 1, 0 /)
     coefficients(923) = -7.38342976908927810000e+000_ark
     expansionIndices(924, :) = (/ 3, 3, 0, 2, 0, 0 /)
     coefficients(924) = -1.01886604238886530000e+001_ark
     expansionIndices(925, :) = (/ 3, 3, 2, 0, 0, 0 /)
     coefficients(925) = -1.01886604238886530000e+001_ark
     expansionIndices(926, :) = (/ 3, 5, 0, 0, 0, 0 /)
     coefficients(926) = -1.36862445572916610000e+002_ark
     expansionIndices(927, :) = (/ 4, 0, 0, 0, 0, 4 /)
     coefficients(927) = 8.21389940686711310000e-003_ark
     expansionIndices(928, :) = (/ 4, 0, 0, 0, 2, 2 /)
     coefficients(928) = 1.64277988137342260000e-002_ark
     expansionIndices(929, :) = (/ 4, 0, 0, 0, 4, 0 /)
     coefficients(929) = 8.21389940686711310000e-003_ark
     expansionIndices(930, :) = (/ 4, 2, 0, 0, 0, 2 /)
     coefficients(930) = 5.76616601093827840000e-001_ark
     expansionIndices(931, :) = (/ 4, 2, 0, 0, 2, 0 /)
     coefficients(931) = 5.76616601093827840000e-001_ark
     expansionIndices(932, :) = (/ 4, 4, 0, 0, 0, 0 /)
     coefficients(932) = 3.61553020077958960000e+001_ark
     expansionIndices(933, :) = (/ 5, 1, 0, 0, 0, 2 /)
     coefficients(933) = -2.25824652086142720000e-001_ark
     expansionIndices(934, :) = (/ 5, 1, 0, 0, 2, 0 /)
     coefficients(934) = -2.25824652086142720000e-001_ark
     expansionIndices(935, :) = (/ 5, 3, 0, 0, 0, 0 /)
     coefficients(935) = -3.16511700980540760000e+000_ark
     expansionIndices(936, :) = (/ 6, 0, 0, 0, 0, 2 /)
     coefficients(936) = 1.61016073233605990000e-002_ark
     expansionIndices(937, :) = (/ 6, 0, 0, 0, 2, 0 /)
     coefficients(937) = 1.61016073233605990000e-002_ark
     expansionIndices(938, :) = (/ 6, 2, 0, 0, 0, 0 /)
     coefficients(938) = 3.70868601613380310000e-001_ark
     expansionIndices(939, :) = (/ 7, 1, 0, 0, 0, 0 /)
     coefficients(939) = -3.72451834151799720000e-001_ark
     expansionIndices(940, :) = (/ 8, 0, 0, 0, 0, 0 /)
     coefficients(940) = 3.71248982774002480000e-003_ark


     f = 0.0_ark
     do i = 1, 940
        addToPotential = 1.0_ark
        do j = 1, 6
           if (expansionIndices(i, j) /= 0) then
              addToPotential = addToPotential*S(j)**expansionIndices(i, j)
           endif
        enddo
        addToPotential = addToPotential*coefficients(i)
        f = f + addToPotential
     enddo
   ! write(*, '(A, F20.10, F20.10, F20.10, F20.10)') "sum of coefficients = ", force(1),  force(2), force(3), force(4)
   ! write(*, '(A, F20.10, F20.10, F20.10, F20.10, F20.10, F20.10, F20.10)') "point ", S(1), S(2), S(3), S(4), S(5), S(6), f
   ! write(*, '(A, F20.10)') "size = ", size(force)
   f = f*219474.631370213_ark
   ! write(*, '(A, F20.10, F20.10, F20.10, F20.10, F20.10, F20.10, A, F20.10, F20.10, F20.10, F20.10, F20.10, F20.10, F20.10)') "local ", r1, r2, r3, alpha1, alpha2, alpha3, " symmeterised ", S(1), S(2), S(3), S(4), S(5), S(6), f
   end function MLpoten
   !
   !!!!start change
    function ML_MEP_OH3P(x)  result(dst)
  
     real(ark),intent(in)   ::  x
     real(ark)              ::  dst
     real(ark)              ::  rf(0:6),y
     integer(ik)          ::  k(0:6) = (/0,1,2,3,4,5,6/)
       !
       rf(0:6)  = molec%mep_params(1:7) 
       !
       y=sin(x)-1.0_ark
       !
       dst =  sum(rf(:)*y**k(:))
       ! 
    end function ML_MEP_OH3P
   !!!!end change
   !
    ! Defining potential energy function 
    !
    ! This type is for XY3-molecules, MORBID type of expansion #10
    ! see, e.g., JCP 117, 11265 (2002) for the defienition 
    !
    function MLpoten_xy3_morbid_Roman_10(ncoords,natoms,local,xyz,force) result(f) 
     !
     integer(ik),intent(in) ::  ncoords,natoms
     real(ark),intent(in)   ::  local(ncoords)
     real(ark),intent(in)   ::  xyz(natoms,3)
     real(ark),intent(in)   ::  force(:)
     real(ark)              ::  force2(304)
     real(ark)              ::  f,f1,f2
     !
     integer(ik)          :: N
  
     real(ark)            ::  r14,r24,r34,alpha1,alpha2,alpha3,alpha0
     real(ark)            ::  y1,y2,y3,y4,y5
     real(ark)            ::  v0,rhoe,v1,v2,v3,v4,v5,v6,rho,re14,aa1,tau
     real(ark)            ::  sinrho,alpha,coro,tau_2,cosalpha,sinphi,phi1,phi2,phi3,delta
        !
        if (verbose>=6) write(out,"('MLpoten_xy3_morbid_10/start')") 
        !
        re14    = molec%req(1)
        !
        alpha0  = molec%alphaeq(1)
        rhoe    = pi-asin(2.0_ark*sin(alpha0*0.5_ark)/sqrt(3.0_ark))
        !
        aa1     = molec%specparam(1)
        !
        r14    = local(1) ;  r24    = local(2) ;  r34    = local(3)
        !
        if (size(local)==6.and.molec%Ndihedrals==1) then 
           !
           select case(trim(molec%coords_transform))
           case default
              !
              write (out,"('MLpoten_xy3_morbid_10: coord. type ',a,' unknown')") trim(molec%coords_transform)
              stop 'MLpoten_xy3_morbid_10 - bad coord. type'
              !
  
           case('R-SYMPHI-TAU')
              !
              !
              alpha3 = local(4)
              alpha2 = local(5)
              !
              tau    = local(6)
              !
              cosalpha = cos(alpha2)*cos(alpha3) + sin(alpha2)*sin(alpha3)*cos(tau)
              !
              if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLpoten_xy3_morbid_45760: cosalpha>1: ',f18.8)") cosalpha
                 stop 'MLpoten_xy3_morbid_45760 - bad cosalpha'
                 !
              endif 
              !
              alpha1 = acos(cosalpha)
              !
           case('R-S-DELTA','R-2D-DELTA','R-S-RHO','SYM-DELTA','R-SYMPHI-DELTA','R-PHI-DELTA','R-S-DELTA-MEP','R-PHI-DELTA-MEP','R-SYMPHI-DELTA-MEP','R-A2A3-DELTA')
              !
              alpha3 = local(4)
              alpha2 = local(5)
              !
              delta = local(6)
              !
              sinphi = sin(alpha2*0.5_ark)/cos(delta)
              phi2 = asin(sinphi)*2.0_ark
              phi2 = mod(phi2+2.0_ark*pi,2.0_ark*pi)
              !
              sinphi = sin(alpha3*0.5_ark)/cos(delta)
              phi3 = asin(sinphi)*2.0_ark
              phi3 = mod(phi3+2.0_ark*pi,2.0_ark*pi)
              phi1 = 2.0_ark*pi-phi2-phi3
              !     
              cosalpha = cos( delta )**2+sin(delta)**2*cos(phi1)
              !
              if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLpoten_xy3_morbid_45760: cosalpha>1: ',f18.8)") cosalpha
                 stop 'MLpoten_xy3_morbid_45760 - bad cosalpha'
                 !
              elseif ( cosalpha>=1.0_ark) then 
                 alpha1 = 0.0_ark
              else 
                 alpha1 = acos(cosalpha)
              endif
              !
           case('R-A2-A3-TAU','R-THETA-TAU')
              !
              !
              alpha2 = local(4)
              alpha3 = local(5)
              !
              tau    = local(6)
              !
              cosalpha = cos(alpha2)*cos(alpha3) + sin(alpha2)*sin(alpha3)*cos(tau)
              !
              if ( abs(cosalpha)>1.0_ark+sqrt(small_) ) then 
                 !
                 write (out,"('MLpoten_xy3_morbid_45760: cosalpha>1: ',f18.8)") cosalpha
                 stop 'MLpoten_xy3_morbid_45760 - bad cosalpha'
                 !
              endif 
              !
              alpha1 = acos(cosalpha)
  
              !
           end select
  
  
        elseif (size(local)==7.and.molec%Ndihedrals==1) then
           !
           alpha3 = local(4)
           alpha2 = local(5)
           alpha1 = local(6)
           tau    = local(7)
           ! 
        elseif (size(local)==6.and.molec%Ndihedrals==0) then 
           !
           alpha3 = local(4)
           alpha2 = local(5)
           alpha1 = local(6)
           !
           tau = sqrt(1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 & 
                     +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3) )
           !
        else
          write (out,"('MLpoten_xy3_morbid_delta: bad locals; local = ',30f18.8)") local(:)
          stop 'MLpoten_xy3_morbid_delta - bad local coordinates'
        endif
        !
        alpha=(alpha1+alpha2+alpha3)/3.0_ark
        !
        if ( 2.0_ark*sin(alpha*0.5_ark)/sqrt(3.0_ark).ge.1.0_ark ) then 
          sinrho=1.0_ark ; rho = 0.5_ark*pi
        else 
          sinrho = 2.0_ark*sin(alpha*0.5_ark)/sqrt(3.0_ark)
          !
          rho=pi-asin(sinrho)
          !
        endif
        !
        coro=(sin(rhoe)-sinrho)
        !
        if (trim(molec%potentype)=='POTEN_XY3_MORBID_10_MEP') re14 = ML_mep_oh3p(rho)
        !
        y1=1.0_ark-exp(-aa1*(r14-re14))
        y2=1.0_ark-exp(-aa1*(r24-re14))
        y3=1.0_ark-exp(-aa1*(r34-re14))
        !
        y4=(2.0_ark*alpha1-alpha2-alpha3)/sqrt(6.0_ark)
        y5=(alpha2-alpha3)/sqrt(2.0_ark)
        !
        f1 = poten_xy3_morbid_10(y1,y2,y3,y4,y5,coro,force)
  
        force2(:)=0.0_ark
  
  force2(1)=-4501.81203917_ark
  force2(2)=-19881.26329389_ark
  force2(3)=-23974.57727407_ark
  force2(4)=-9334.08769568_ark
  force2(6)=-322.99732624_ark
  force2(12)=5732.64717612_ark
  force2(13)=13470.10273810_ark
  force2(14)=10831.73624903_ark
  force2(15)=2696.17194994_ark
  force2(19)=-390.35384268_ark
  force2(26)=-199529.16354149_ark
  force2(27)=-744261.41099013_ark
  force2(28)=-1006934.10723694_ark
  force2(29)=-576639.05043630_ark
  force2(30)=-114345.27125559_ark
  force2(31)=-28897.81099237_ark
  force2(32)=-123077.58667199_ark
  force2(33)=-196547.81666239_ark
  force2(34)=-139124.46991904_ark
  force2(35)=-36788.63210179_ark
  force2(36)=1959.66818301_ark
  force2(37)=1306.96225957_ark
  force2(38)=-3189.04443438_ark
  force2(39)=-2357.13560627_ark
  force2(40)=-1180.58134135_ark
  force2(41)=-1405.50631046_ark
  force2(42)=624.03401645_ark
  force2(43)=853.31266981_ark
  force2(44)=42446.96444091_ark
  force2(45)=139604.39540308_ark
  force2(46)=152760.79415013_ark
  force2(47)=55622.89083767_ark
  force2(48)=1928.88517275_ark
  force2(49)=3306.00669222_ark
  force2(50)=1337.21817462_ark
  force2(52)=14899.53219932_ark
  force2(53)=49070.51150982_ark
  force2(54)=53871.52094367_ark
  force2(55)=19720.93341366_ark
  force2(56)=8002.31666938_ark
  force2(57)=26001.07795769_ark
  force2(58)=28166.98161352_ark
  force2(59)=10165.90720340_ark
  force2(60)=-7373.48068441_ark
  force2(61)=-23883.57111100_ark
  force2(62)=-25628.40548148_ark
  force2(63)=-9138.11334295_ark
  force2(64)=-7446.34665931_ark
  force2(65)=-24161.39395965_ark
  force2(66)=-26037.86030725_ark
  force2(67)=-9324.29441909_ark
  force2(68)=-2170.24940951_ark
  force2(69)=-4331.72359507_ark
  force2(70)=-2131.17374574_ark
  force2(71)=939.09024088_ark
  force2(72)=1915.74690059_ark
  force2(73)=985.54062584_ark
  force2(74)=-4.02427755_ark
  
       coro=-sinrho
        y1=r14!-re14
        y2=r24!-re14
        y3=r34!-re14
  !
        f2 = poten_xy3_morbid_10(y1,y2,y3,y4,y5,coro,force2)
  !      !
        f=f1+f2*force(305)
  !      !
        if (verbose>=6) write(out,"('MLpoten_xy3_morbid_10/end')") 
   
   end function MLpoten_xy3_morbid_Roman_10
  
  
  
   function poten_xy3_morbid_10(y1,y2,y3,y4,y5,coro,force) result (f)
     !
     real(ark),intent(in) ::  y1,y2,y3,y4,y5,coro
     real(ark),intent(in) ::  force(:)
     real(ark)            ::  f
     !
     integer(ik)          :: N
  
     real(ark)            ::  v0,v1,v2,v3,v4,v5,v6
  
  
     real(ark)            ::  &
        fea1  ,                                      &
        fea11  ,fea12  ,fea14  ,fea44  ,                     &
        fea111 ,fea112 ,fea114 ,fea123 ,                     &
        fea124 ,fea144 ,fea155 ,fea455 ,                     &
        fea1111,fea1112,fea1114,fea1122,                     &
        fea1123,fea1124,fea1125,fea1144,                     &
        fea1155,fea1244,fea1255,fea1444,                     &
        fea1455,fea4444
     real(ark)            ::  &
        fea44444 ,fea33455 ,fea33445 ,fea33345 ,fea33344 ,&
        fea33334 ,fea33333 ,fea25555 ,fea24455 ,fea24445 ,fea23333 ,&
        fea13455 ,fea13445 ,fea13345 ,fea12355 ,fea11334 ,fea11333 ,&
        fea11255 ,fea11245 ,fea11234 ,fea11233 ,fea11135 ,fea11134 ,&
        fea11123 ,fea555555,fea444444,fea335555,fea334455,fea334445,&
        fea333555,fea333333,fea244555,fea244455,fea233445,fea233444,&
        fea233345,fea233344,fea233335,fea223355,fea222335,fea222334,&
        fea222333,fea222255,fea222245,fea222233,fea222224,fea145555,&
        fea134444,fea133444,fea133345,fea133334,fea133333,fea124555,&
        fea124455,fea123455,fea123345,fea113555,fea113345,fea112355,&
        fea112335,fea112233,fea111444,fea111234,fea111233,fea111123
        !
     real(ark)            ::  s1,s2,s3,s4,s5,tau
        !
     real(ark)            ::  &      
        ve  ,            &
        f0a1,f1a,f2a,f3a,f4a,f5a,f6a,f7a,f8a, &
        f1a1,f2a1,f3a1,f4a1,f5a1,f6a1,  &
        f0a11,f1a11,f2a11,f3a11,f4a11, &
        f0a12,f1a12,f2a12,f3a12,f4a12, &
        f0a14,f1a14,f2a14,f3a14,f4a14, &
        f0a44,f1a44,f2a44,f3a44,f4a44, &
        f0a111,f1a111,f2a111,f3a111  , &
        f0a112,f1a112,f2a112,f3a112  , &
        f0a114,f1a114,f2a114,f3a114  , &
        f0a123,f1a123,f2a123,f3a123  , &
        f0a124,f1a124,f2a124,f3a124  , &
        f0a144,f1a144,f2a144,f3a144  , &
        f0a155,f1a155,f2a155,f3a155  , &
        f0a455,f1a455,f2a455,f3a455  , &
        f0a1111,f1a1111,f2a1111      , &
        f0a1112,f1a1112,f2a1112      , &
        f0a1114,f1a1114,f2a1114      , &
        f0a1122,f1a1122,f2a1122      , &
        f0a1123,f1a1123,f2a1123      , &
        f0a1124,f1a1124,f2a1124      , &
        f0a1125,f1a1125,f2a1125      , &
        f0a1144,f1a1144,f2a1144      , &
        f0a1155,f1a1155,f2a1155      , &
        f0a1244,f1a1244,f2a1244      , &
        f0a1255,f1a1255,f2a1255      , &
        f0a1444,f1a1444,f2a1444      , &
        f0a1455,f1a1455,f2a1455      , &
        f0a4444,f1a4444,f2a4444      , &
        f0a44444 ,f1a44444 
  
     real(ark)            ::  &      
        f2a44444 ,f0a33455 ,f1a33455 ,f2a33455 ,f0a33445 ,f1a33445 ,&
        f2a33445 ,f0a33345 ,f1a33345 ,f2a33345 ,f0a33344 ,f1a33344 ,&
        f2a33344 ,f0a33334 ,f1a33334 ,f2a33334 ,f0a33333 ,f1a33333 ,&
        f2a33333 ,f0a25555 ,f1a25555 ,f2a25555 ,f0a24455 ,f1a24455 ,&
        f2a24455 ,f0a24445 ,f1a24445 ,f2a24445 ,f0a23333 ,f1a23333 ,&
        f2a23333 ,f0a13455 ,f1a13455 ,f2a13455 ,f0a13445 ,f1a13445 ,&
        f2a13445 ,f0a13345 ,f1a13345 ,f2a13345 ,f0a12355 ,f1a12355 ,&
        f2a12355 ,f0a11334 ,f1a11334 ,f2a11334 ,f0a11333 ,f1a11333 ,&
        f2a11333 ,f0a11255 ,f1a11255 ,f2a11255 ,f0a11245 ,f1a11245 ,&
        f2a11245 ,f0a11234 ,f1a11234 ,f2a11234 ,f0a11233 ,f1a11233 ,&
        f2a11233 ,f0a11135 ,f1a11135 ,f2a11135 ,f0a11134 ,f1a11134 ,&
        f2a11134 ,f0a11123 ,f1a11123 ,f2a11123 ,f0a555555,f1a555555
  
  
     real(ark)            ::  &      
        f2a555555,f0a444444,f1a444444,f2a444444,f0a335555,f1a335555,&
        f2a335555,f0a334455,f1a334455,f2a334455,f0a334445,f1a334445,&
        f2a334445,f0a333555,f1a333555,f2a333555,f0a333333,f1a333333,&
        f2a333333,f0a244555,f1a244555,f2a244555,f0a244455,f1a244455,&
        f2a244455,f0a233445,f1a233445,f2a233445,f0a233444,f1a233444,&
        f2a233444,f0a233345,f1a233345,f2a233345,f0a233344,f1a233344,&
        f2a233344,f0a233335,f1a233335,f2a233335,f0a223355,f1a223355,&
        f2a223355,f0a222335,f1a222335,f2a222335,f0a222334,f1a222334,&
        f2a222334,f0a222333,f1a222333,f2a222333,f0a222255,f1a222255,&
        f2a222255,f0a222245,f1a222245,f2a222245,f0a222233,f1a222233,&
        f2a222233,f0a222224,f1a222224,f2a222224,f0a145555,f1a145555,&
        f2a145555,f0a134444,f1a134444,f2a134444,f0a133444,f1a133444,&
        f2a133444,f0a133345,f1a133345,f2a133345,f0a133334,f1a133334,&
        f2a133334,f0a133333,f1a133333,f2a133333,f0a124555,f1a124555,&
        f2a124555,f0a124455,f1a124455,f2a124455,f0a123455,f1a123455,&
        f2a123455,f0a123345,f1a123345,f2a123345,f0a113555,f1a113555,&
        f2a113555,f0a113345,f1a113345,f2a113345,f0a112355,f1a112355,&
        f2a112355,f0a112335,f1a112335,f2a112335,f0a112233,f1a112233,&
        f2a112233,f0a111444,f1a111444,f2a111444,f0a111234,f1a111234,&
        f2a111234,f0a111233,f1a111233,f2a111233,f0a111123,f1a111123,&
        f2a111123
  
  
        N = size(force)
        !
        ve         = force(  1)
        f1a        = force(  2)
        f2a        = force(  3)
        f3a        = force(  4)
        f4a        = force(  5)
        f5a        = force(  6)
        f6a        = force(  7)
        f7a        = force(  8)
        f0a1       = force(  9)
        f1a1       = force( 10)
        f2a1       = force( 11)
        f3a1       = force( 12)
        f4a1       = force( 13)
        f5a1       = force( 14)
        f6a1       = force( 15)
        f0a11      = force( 16)
        f1a11      = force( 17)
        f2a11      = force( 18)
        f3a11      = force( 19)
        f4a11      = force( 20)
        f0a12      = force( 21)
        f1a12      = force( 22)
        f2a12      = force( 23)
        f3a12      = force( 24)
        f4a12      = force( 25)
        f0a14      = force( 26)
        f1a14      = force( 27)
        f2a14      = force( 28)
        f3a14      = force( 29)
        f4a14      = force( 30)
        f0a44      = force( 31)
        f1a44      = force( 32)
        f2a44      = force( 33)
        f3a44      = force( 34)
        f4a44      = force( 35)
  
        v0=ve+f1a*coro+f2a*coro**2+f3a*coro**3+f4a*coro**4+f5a*coro**5 &
              +f6a*coro**6+f7a*coro**7  !  +f8a*coro**8
  
        fea1= f0a1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4+f5a1*coro**5+f6a1*coro**6
        !
        fea11=   f0a11+f1a11*coro+f2a11*coro**2+f3a11*coro**3+f4a11*coro**4
        fea12=   f0a12+f1a12*coro+f2a12*coro**2+f3a12*coro**3+f4a12*coro**4
        fea14=   f0a14+f1a14*coro+f2a14*coro**2+f3a14*coro**3+f4a14*coro**4
        fea44=   f0a44+f1a44*coro+f2a44*coro**2+f3a44*coro**3+f4a44*coro**4
        !
        !
        v1 = (y3+y2+y1)*fea1
        !
        v2 = (y2*y3+y1*y3+y1*y2)*fea12                                                                 &
         +(y2**2+y3**2+y1**2)*fea11                                                                    &
         +(-sqrt(3.0_ark)*y3*y5/2.0_ark-y3*y4/2.0_ark+y1*y4+sqrt(3.0_ark)*y2*y5/2.0_ark-y2*y4/2.0_ark)*fea14 &
         +(y5**2+y4**2)*fea44
        !
        v3 = 0 ; v4 = 0 ; v5 = 0 ; v6 = 0 
        !
        if (N>35) then  
          f0a111     = force( 36)
          f1a111     = force( 37)
          f2a111     = force( 38)
          f3a111     = force( 39)
          f0a112     = force( 40)
          f1a112     = force( 41)
          f2a112     = force( 42)
          f3a112     = force( 43)
          f0a114     = force( 44)
          f1a114     = force( 45)
          f2a114     = force( 46)
          f3a114     = force( 47)
          f0a123     = force( 48)
          f1a123     = force( 49)
          f2a123     = force( 50)
          f3a123     = force( 51)
          f0a124     = force( 52)
          f1a124     = force( 53)
          f2a124     = force( 54)
          f3a124     = force( 55)
          f0a144     = force( 56)
          f1a144     = force( 57)
          f2a144     = force( 58)
          f3a144     = force( 59)
          f0a155     = force( 60)
          f1a155     = force( 61)
          f2a155     = force( 62)
          f3a155     = force( 63)
          f0a455     = force( 64)
          f1a455     = force( 65)
          f2a455     = force( 66)
          f3a455     = force( 67)
          !
          fea111= f0a111+f1a111*coro+f2a111*coro**2+f3a111*coro**3
          fea112= f0a112+f1a112*coro+f2a112*coro**2+f3a112*coro**3
          fea114= f0a114+f1a114*coro+f2a114*coro**2+f3a114*coro**3
          fea123= f0a123+f1a123*coro+f2a123*coro**2+f3a123*coro**3
          fea124= f0a124+f1a124*coro+f2a124*coro**2+f3a124*coro**3
          fea144= f0a144+f1a144*coro+f2a144*coro**2+f3a144*coro**3
          fea155= f0a155+f1a155*coro+f2a155*coro**2+f3a155*coro**3
          fea455= f0a455+f1a455*coro+f2a455*coro**2+f3a455*coro**3
          !
          v3 = (y1*y3*y4+y1*y2*y4-2.0_ark*y2*y3*y4+sqrt(3.0_ark)*y1*y2*y5-sqrt(3.0_ark)*y1*y3*y5)*fea124   &
           +(3.0_ark/4.0_ark*y3*y4**2-sqrt(3.0_ark)*y3*y4*y5/2.0_ark+y1*y5**2+y2*y5**2/4.0_ark               & 
           +3.0_ark/4.0_ark*y2*y4**2+sqrt(3.0_ark)*y2*y4*y5/2.0_ark+y3*y5**2/4.0_ark)*fea155                 &
           +(y2*y3**2+y1*y3**2+y1**2*y3+y1*y2**2+y2**2*y3+y1**2*y2)*fea112+                             &
           (-y4**3/3.0_ark+y4*y5**2)*fea455+fea123*y1*y2*y3                                              &
           +(y1*y4**2+3.0_ark/4.0_ark*y3*y5**2+3.0_ark/4.0_ark*y2*y5**2+y2*y4**2/4.0_ark                     &
           -sqrt(3.0_ark)*y2*y4*y5/2.0_ark+sqrt(3.0_ark)*y3*y4*y5/2.0_ark+y3*y4**2/4.0_ark)*fea144           &
           +(y3**3+y2**3+y1**3)*fea111+(-y2**2*y4/2.0_ark-y3**2*y4/2.0_ark+sqrt(3.0_ark)*y2**2*y5/2.0_ark   & 
           +y1**2*y4-sqrt(3.0_ark)*y3**2*y5/2.0_ark)*fea114
           !
        endif
        !
        if (N>67) then  
          !
          f0a1111    = force( 68)
          f1a1111    = force( 69)
          f2a1111    = force( 70)
          f0a1112    = force( 71)
          f1a1112    = force( 72)
          f2a1112    = force( 73)
          f0a1114    = force( 74)
          f1a1114    = force( 75)
          f2a1114    = force( 76)
          f0a1122    = force( 77)
          f1a1122    = force( 78)
          f2a1122    = force( 79)
          f0a1123    = force( 80)
          f1a1123    = force( 81)
          f2a1123    = force( 82)
          f0a1124    = force( 83)
          f1a1124    = force( 84)
          f2a1124    = force( 85)
          f0a1125    = force( 86)
          f1a1125    = force( 87)
          f2a1125    = force( 88)
          f0a1144    = force( 89)
          f1a1144    = force( 90)
          f2a1144    = force( 91)
          f0a1155    = force( 92)
          f1a1155    = force( 93)
          f2a1155    = force( 94)
          f0a1244    = force( 95)
          f1a1244    = force( 96)
          f2a1244    = force( 97)
          f0a1255    = force( 98)
          f1a1255    = force( 99)
          f2a1255    = force(100)
          f0a1444    = force(101)
          f1a1444    = force(102)
          f2a1444    = force(103)
          f0a1455    = force(104)
          f1a1455    = force(105)
          f2a1455    = force(106)
          f0a4444    = force(107)
          f1a4444    = force(108)
          f2a4444    = force(109)
          !
          fea1111= f0a1111+f1a1111*coro+f2a1111*coro**2
          fea1112= f0a1112+f1a1112*coro+f2a1112*coro**2
          fea1114= f0a1114+f1a1114*coro+f2a1114*coro**2
          fea1122= f0a1122+f1a1122*coro+f2a1122*coro**2
          fea1123= f0a1123+f1a1123*coro+f2a1123*coro**2
          fea1124= f0a1124+f1a1124*coro+f2a1124*coro**2
          fea1125= f0a1125+f1a1125*coro+f2a1125*coro**2
          fea1144= f0a1144+f1a1144*coro+f2a1144*coro**2
          fea1155= f0a1155+f1a1155*coro+f2a1155*coro**2
          fea1244= f0a1244+f1a1244*coro+f2a1244*coro**2
          fea1255= f0a1255+f1a1255*coro+f2a1255*coro**2
          fea1444= f0a1444+f1a1444*coro+f2a1444*coro**2
          fea1455= f0a1455+f1a1455*coro+f2a1455*coro**2
          fea4444= f0a4444+f1a4444*coro+f2a4444*coro**2
          !
          s2 = (y4**4+y5**4+2.0_ark*y4**2*y5**2)*fea4444+(3.0_ark/8.0_ark*sqrt(3.0_ark)*&
           y2*y5**3-3.0_ark/8.0_ark*sqrt(3.0_ark)*y3*y4**2*y5-3.0_ark/8.0_ark*sqrt(3.0_ark)*y3*&
           y5**3-9.0_ark/8.0_ark*y2*y4*y5**2-y3*y4**3/8.0_ark-y2*y4**3/8.0_ark-9.0_ark/8.0_ark*&
           y3*y4*y5**2+y1*y4**3+3.0_ark/8.0_ark*sqrt(3.0_ark)*y2*y4**2*y5)*fea1444 &
           +(3.0_ark/4.0_ark*y2**2*y4**2+3.0_ark/4.0_ark*y3**2*y4**2+y1**2*y5**2+y3**2*y5**2/4.0_ark &
           -sqrt(3.0_ark)*y3**2*y4*y5/2.0_ark+sqrt(3.0_ark)*y2**2*y4*y5/2.0_ark+y2**2&
           *y5**2/4.0_ark)*fea1155 
           s1 = s2+(y3**2*y4**2/4.0_ark+3.0_ark/4.0_ark*y3**2*y5**2+y1**2*y4**2+y2**2*&
           y4**2/4.0_ark+sqrt(3.0_ark)*y3**2*y4*y5/2.0_ark-sqrt(3.0_ark)*y2**2*y4*y5/2.0_ark&
           +3.0_ark/4.0_ark*y2**2*y5**2)*fea1144+(y1**3*y4+sqrt(3.0_ark)*y2**3*y5/2.0_ark&
           -sqrt(3.0_ark)*y3**3*y5/2.0_ark-y2**3*y4/2.0_ark-y3**3*y4/2.0_ark)*fea1114&
           +(y2**4+y1**4+y3**4)*fea1111+(sqrt(3.0_ark)*y1*y3*y4*y5+3.0_ark/2.0_ark*y2*y3*y5**2&
           -y2*y3*y4**2/2.0_ark+y1*y2*y4**2-sqrt(3.0_ark)*y1*y2*y4*y5+y1*y3*y4**2)*fea1244 
           !
          s2 = s1+(y1*y3*y5**2+y1*y2*y5**2-sqrt(3.0_ark)*y1*y3*y4*y5-y2*y3*y5**&
           2/2.0_ark+3.0_ark/2.0_ark*y2*y3*y4**2+sqrt(3.0_ark)*y1*y2*y4*y5)*fea1255+&
           (-y1*y3**2*y4/2.0_ark+y1**2*y3*y4-sqrt(3.0_ark)*y1*y3**2*y5/2.0_ark-sqrt(3.0_ark)*y2&
           *y3**2*y5/2.0_ark+y1**2*y2*y4+sqrt(3.0_ark)*y2**2*y3*y5/2.0_ark-y2**2*y3*y4&
           /2.0_ark+sqrt(3.0_ark)*y1*y2**2*y5/2.0_ark-y2*y3**2*y4/2.0_ark-y1*y2**2*y4/2.0_ark&
           )*fea1124+(y1**2*y2*y5+sqrt(3.0_ark)*y1*y3**2*y4/2.0_ark+sqrt(3.0_ark)*y1*&
           y2**2*y4/2.0_ark-sqrt(3.0_ark)*y2*y3**2*y4/2.0_ark-sqrt(3.0_ark)*y2**2*y3*y4/2.0_ark&
           -y2**2*y3*y5/2.0_ark+y2*y3**2*y5/2.0_ark-y1*y3**2*y5/2.0_ark+y1*y2**2*y5&
           /2.0_ark-y1**2*y3*y5)*fea1125 
           !
          v4 = s2+(y2*y3**3+y1**3*y3+y1**3*y2+y1*y2**3+y1*y3**3+y2**3*y3)*fea1112+&
           (y2**2*y3**2+y1**2*y3**2+y1**2*y2**2)*fea1122+(y1*y2**2*y3&
           +y1**2*y2*y3+y1*y2*y3**2)*fea1123+(5.0_ark/8.0_ark*y2*y4*y5**2+sqrt(3.0_ark)*&
           y2*y5**3/8.0_ark-sqrt(3.0_ark)*y3*y4**2*y5/8.0_ark+sqrt(3.0_ark)*y2*y4**2*y5/8.0_ark&
           -3.0_ark/8.0_ark*y2*y4**3+y1*y4*y5**2-sqrt(3.0_ark)*y3*y5**3/8.0_ark&
           +5.0_ark/8.0_ark*y3*y4*y5**2-3.0_ark/8.0_ark*y3*y4**3)*fea1455
          !
        endif
        !
        if (N>109) then  
          !
          f0a44444   = force(110)
          f1a44444   = force(111)
          f2a44444   = force(112)
          f0a33455   = force(113)
          f1a33455   = force(114)
          f2a33455   = force(115)
          f0a33445   = force(116)
          f1a33445   = force(117)
          f2a33445   = force(118)
          f0a33345   = force(119)
          f1a33345   = force(120)
          f2a33345   = force(121)
          f0a33344   = force(122)
          f1a33344   = force(123)
          f2a33344   = force(124)
          f0a33334   = force(125)
          f1a33334   = force(126)
          f2a33334   = force(127)
          f0a33333   = force(128)
          f1a33333   = force(129)
          f2a33333   = force(130)
          f0a25555   = force(131)
          f1a25555   = force(132)
          f2a25555   = force(133)
          f0a24455   = force(134)
          f1a24455   = force(135)
          f2a24455   = force(136)
          f0a24445   = force(137)
          f1a24445   = force(138)
          f2a24445   = force(139)
          f0a23333   = force(140)
          f1a23333   = force(141)
          f2a23333   = force(142)
          f0a13455   = force(143)
          f1a13455   = force(144)
          f2a13455   = force(145)
          f0a13445   = force(146)
          f1a13445   = force(147)
          f2a13445   = force(148)
          f0a13345   = force(149)
          f1a13345   = force(150)
          f2a13345   = force(151)
          f0a12355   = force(152)
          f1a12355   = force(153)
          f2a12355   = force(154)
          f0a11334   = force(155)
          f1a11334   = force(156)
          f2a11334   = force(157)
          f0a11333   = force(158)
          f1a11333   = force(159)
          f2a11333   = force(160)
          f0a11255   = force(161)
          f1a11255   = force(162)
          f2a11255   = force(163)
          f0a11245   = force(164)
          f1a11245   = force(165)
          f2a11245   = force(166)
          f0a11234   = force(167)
          f1a11234   = force(168)
          f2a11234   = force(169)
          f0a11233   = force(170)
          f1a11233   = force(171)
          f2a11233   = force(172)
          f0a11135   = force(173)
          f1a11135   = force(174)
          f2a11135   = force(175)
          f0a11134   = force(176)
          f1a11134   = force(177)
          f2a11134   = force(178)
          f0a11123   = force(179)
          f1a11123   = force(180)
          f2a11123   = force(181)
  
          fea44444 = f0a44444  + f1a44444 *coro+ f2a44444 *coro**2
          fea33455 = f0a33455  + f1a33455 *coro+ f2a33455 *coro**2
          fea33445 = f0a33445  + f1a33445 *coro+ f2a33445 *coro**2
          fea33345 = f0a33345  + f1a33345 *coro+ f2a33345 *coro**2
          fea33344 = f0a33344  + f1a33344 *coro+ f2a33344 *coro**2
          fea33334 = f0a33334  + f1a33334 *coro+ f2a33334 *coro**2
          fea33333 = f0a33333  + f1a33333 *coro+ f2a33333 *coro**2
          fea25555 = f0a25555  + f1a25555 *coro+ f2a25555 *coro**2
          fea24455 = f0a24455  + f1a24455 *coro+ f2a24455 *coro**2
          fea24445 = f0a24445  + f1a24445 *coro+ f2a24445 *coro**2
          fea23333 = f0a23333  + f1a23333 *coro+ f2a23333 *coro**2
          fea13455 = f0a13455  + f1a13455 *coro+ f2a13455 *coro**2
          fea13445 = f0a13445  + f1a13445 *coro+ f2a13445 *coro**2
          fea13345 = f0a13345  + f1a13345 *coro+ f2a13345 *coro**2
          fea12355 = f0a12355  + f1a12355 *coro+ f2a12355 *coro**2
          fea11334 = f0a11334  + f1a11334 *coro+ f2a11334 *coro**2
          fea11333 = f0a11333  + f1a11333 *coro+ f2a11333 *coro**2
          fea11255 = f0a11255  + f1a11255 *coro+ f2a11255 *coro**2
          fea11245 = f0a11245  + f1a11245 *coro+ f2a11245 *coro**2
          fea11234 = f0a11234  + f1a11234 *coro+ f2a11234 *coro**2
          fea11233 = f0a11233  + f1a11233 *coro+ f2a11233 *coro**2
          fea11135 = f0a11135  + f1a11135 *coro+ f2a11135 *coro**2
          fea11134 = f0a11134  + f1a11134 *coro+ f2a11134 *coro**2
          fea11123 = f0a11123  + f1a11123 *coro+ f2a11123 *coro**2
          !    
          s3 = (y4**5-2.0_ark*y4**3*y5**2-3.0_ark*y4*y5**4)*fea44444+(-4.0_ark*y3*y4*&
          y5**3*sqrt(3.0_ark)+9.0_ark*y1*y4**2*y5**2-3.0_ark/2.0_ark*y1*y4**4+4.0_ark*y2*y4&
          *y5**3*sqrt(3.0_ark)+3.0_ark*y2*y4**4+5.0_ark/2.0_ark*y1*y5**4+3.0_ark*y3*y4**4+&
          y2*y5**4+y3*y5**4)*fea25555+(-y2*y4**4+y3*y4**2*y5**2-2.0_ark*y2*y4*y5&
          **3*sqrt(3.0_ark)-y3*y4**4-7.0_ark/2.0_ark*y1*y4**2*y5**2-3.0_ark/4.0_ark*y1*y5**4&
          +2.0_ark*y3*y4*y5**3*sqrt(3.0_ark)+y2*y4**2*y5**2+5.0_ark/4.0_ark*y1*y4**4)*fea24455 
          !
          s2 = s3+(y2*y4**3*y5-3.0_ark*y3*y4*y5**3+2.0_ark/3.0_ark*y3*y4**4*sqrt(3.0_ark&
          )+3.0_ark/4.0_ark*y1*y5**4*sqrt(3.0_ark)+3.0_ark*y2*y4*y5**3-&
          7.0_ark/12.0_ark*y1*y4**4*sqrt(3.0_ark)+3.0_ark/2.0_ark*y1*y4**2*y5**2*sqrt(3.0_ark)-y3*y4**3*y5&
          +2.0_ark/3.0_ark*y2*y4**4*sqrt(3.0_ark))*fea24445+(-y2**2*y5**3+y3**2*y4**2*y5+ &
          y3**2*y5**3+4.0_ark/9.0_ark*y2**2*y4**3*sqrt(3.0_ark)-5.0_ark/9.0_ark*y1**2*y4**3*&
          sqrt(3.0_ark)+4.0_ark/9.0_ark*y3**2*y4**3*sqrt(3.0_ark)-y2**2*y4**2*y5&
          -y1**2*y4*y5**2*sqrt(3.0_ark))*fea33445+(y3**2*y4*y5**2-y1**2*y4**3/3.0_ark&
          -y3**2*y4**3/3.0_ark+y1**2*y4*y5**2+y2**2*y4*y5**2-y2**2*y4**3/3.0_ark)*fea33455 
          !
          s1 = s2+(-y2**3*y4*y5+y3**3*y4*y5+y2**3*y5**2*sqrt(3.0_ark)/3.0_ark+y1**&
          3*y4**2*sqrt(3.0_ark)/2.0_ark+y3**3*y5**2*sqrt(3.0_ark)/3.0_ark- & 
          y1**3*y5**2*sqrt(3.0_ark)/6.0_ark)*fea33345+(y3**3*y4**2+y3**3*y5**2+y2**3*y4**2+y2**3&
          *y5**2+y1**3*y5**2+y1**3*y4**2)*fea33344+(y3**4*y4+sqrt(3.0_ark)*y3**&
          4*y5+y2**4*y4-2.0_ark*y1**4*y4-sqrt(3.0_ark)*y2**4*y5)*fea33334+(y2**5+ &
          y3**5+y1**5)*fea33333+(-4.0_ark/9.0_ark*y1*y2*y4**3*sqrt(3.0_ark)-y1*y2*y5**3+ &
          y1*y3*y4**2*y5+y2*y3*y4*y5**2*sqrt(3.0_ark)-y1*y2*y4**2*y5+5.0_ark/9.0_ark&
          *y2*y3*y4**3*sqrt(3.0_ark)-4.0_ark/9.0_ark*y1*y3*y4**3*sqrt(3.0_ark)+y1*y3*y5&
          **3)*fea13445+(y2*y3*y4*y5**2+y1*y2*y4*y5**2-y2*y3*y4**3/3.0_ark- & 
          y1*y2*y4**3/3.0_ark-y1*y3*y4**3/3.0_ark+y1*y3*y4*y5**2)*fea13455 
          
          s3 = s1+(y1**2*y3*y5**2+y2**2*y3*y4**2+y2**2*y3*y5**2+y1*y2**2*y5**2+&
          y1**2*y2*y5**2+y1*y2**2*y4**2+y2*y3**2*y4**2+y1*y3**2*y4**2+&
          y1**2*y3*y4**2+y1**2*y2*y4**2+y1*y3**2*y5**2+y2*y3**2*y5**2)*fea11255&
          +(2.0_ark/3.0_ark*y1**2*y3*y4**2*sqrt(3.0_ark)+y1*y3**2*y5**2*sqrt(3.0_ark)/2.0_ark+&
          y1*y2**2*y5**2*sqrt(3.0_ark)/2.0_ark+y2**2*y3*y5**2*sqrt(3.0_ark)/2.0_ark-&
          y1*y2**2*y4*y5+y2*y3**2*y4*y5+y1*y3**2*y4*y5-y2**2*y3*y4*y5+y2*y3**&
          2*y4**2*sqrt(3.0_ark)/6.0_ark+y1*y3**2*y4**2*sqrt(3.0_ark)/6.0_ark+y1*y2**2*y4&
          **2*sqrt(3.0_ark)/6.0_ark+2.0_ark/3.0_ark*y1**2*y2*y4**2*sqrt(3.0_ark)+&
          y2*y3**2*y5**2*sqrt(3.0_ark)/2.0_ark+y2**2*y3*y4**2*sqrt(3.0_ark)/6.0_ark)*fea13345 
          s4 = s3+(y1**2*y2*y4*y5+y1**2*y3*y4**2*sqrt(3.0_ark)/3.0_ark+y1**2*y2*y4&
          **2*sqrt(3.0_ark)/3.0_ark-y1*y2**2*y4**2*sqrt(3.0_ark)/6.0_ark+y2*y3**2*y4*y5-&
          y2**2*y3*y4*y5-y1**2*y3*y4*y5+y2*y3**2*y4**2*sqrt(3.0_ark)/3.0_ark+y1*y2&
          **2*y5**2*sqrt(3.0_ark)/2.0_ark-y1*y3**2*y4**2*sqrt(3.0_ark)/6.0_ark+y2**2*y3*&
          y4**2*sqrt(3.0_ark)/3.0_ark+y1*y3**2*y5**2*sqrt(3.0_ark)/2.0_ark)*fea11245 
          s2 = s4+(-y1**3*y2*y5+y1**3*y3*y5+y2**3*y3*y5/2.0_ark-y1*y2**3*y4*sqrt(3.0_ark)/2.0_ark-&
          y1*y2**3*y5/2.0_ark-y2*y3**3*y5/2.0_ark+y1*y3**3*y5/2.0_ark+y2&
          **3*y3*y4*sqrt(3.0_ark)/2.0_ark+y2*y3**3*y4*sqrt(3.0_ark)/2.0_ark-y1*y3**3*y4*&
          sqrt(3.0_ark)/2.0_ark)*fea11135+(y1**3*y3*y4-y2**3*y3*y4/2.0_ark+y1**3*y2*y4-&
          y2*y3**3*y4/2.0_ark-y1*y3**3*y4/2.0_ark+y1*y2**3*y5*sqrt(3.0_ark)/2.0_ark+y2&
          **3*y3*y5*sqrt(3.0_ark)/2.0_ark-y2*y3**3*y5*sqrt(3.0_ark)/2.0_ark-y1*y2**3*y4/&
          2.0_ark-y1*y3**3*y5*sqrt(3.0_ark)/2.0_ark)*fea11134 
          
          v5 = s2+(y1*y2**4+y1**4*y3+y1**4*y2+y2**4*y3+y2*y3**4+y1*y3**4)*fea23333+&
          (-2.0_ark*y2**2*y3**2*y4+y1**2*y2**2*y4-sqrt(3.0_ark)*y1**2*y3**2&
          *y5+sqrt(3.0_ark)*y1**2*y2**2*y5+y1**2*y3**2*y4)*fea11334+(y1**2*y3**&
          3+y1**3*y3**2+y2**2*y3**3+y1**2*y2**3+y1**3*y2**2+y2**3*y3**2)*fea11333+&
          (y1*y2*y3*y4**2+y1*y2*y3*y5**2)*fea12355+(-y1*y2*y3**2*y4/2.0_ark-&
          y1*y2**2*y3*y4/2.0_ark-sqrt(3.0_ark)*y1*y2*y3**2*y5/2.0_ark+y1**2*y2*y3*y4+&
          sqrt(3.0_ark)*y1*y2**2*y3*y5/2.0_ark)*fea11234+(y1*y2**3*y3+y1*y2*y3**3+&
          y1**3*y2*y3)*fea11123+(y1**2*y2**2*y3+y1*y2**2*y3**2+y1**2*y2*y3**2)*fea11233
          !
        endif
        !
        if (N>181) then  
          !
          f0a555555  = force(182)
          f1a555555  = force(183)
          f2a555555  = force(184)
          f0a444444  = force(185)
          f1a444444  = force(186)
          f2a444444  = force(187)
          f0a335555  = force(188)
          f1a335555  = force(189)
          f2a335555  = force(190)
          f0a334455  = force(191)
          f1a334455  = force(192)
          f2a334455  = force(193)
          f0a334445  = force(194)
          f1a334445  = force(195)
          f2a334445  = force(196)
          f0a333555  = force(197)
          f1a333555  = force(198)
          f2a333555  = force(199)
          f0a333333  = force(200)
          f1a333333  = force(201)
          f2a333333  = force(202)
          f0a244555  = force(203)
          f1a244555  = force(204)
          f2a244555  = force(205)
          f0a244455  = force(206)
          f1a244455  = force(207)
          f2a244455  = force(208)
          f0a233445  = force(209)
          f1a233445  = force(210)
          f2a233445  = force(211)
          f0a233444  = force(212)
          f1a233444  = force(213)
          f2a233444  = force(214)
          f0a233345  = force(215)
          f1a233345  = force(216)
          f2a233345  = force(217)
          f0a233344  = force(218)
          f1a233344  = force(219)
          f2a233344  = force(220)
          f0a233335  = force(221)
          f1a233335  = force(222)
          f2a233335  = force(223)
          f0a223355  = force(224)
          f1a223355  = force(225)
          f2a223355  = force(226)
          f0a222335  = force(227)
          f1a222335  = force(228)
          f2a222335  = force(229)
          f0a222334  = force(230)
          f1a222334  = force(231)
          f2a222334  = force(232)
          f0a222333  = force(233)
          f1a222333  = force(234)
          f2a222333  = force(235)
          f0a222255  = force(236)
          f1a222255  = force(237)
          f2a222255  = force(238)
          f0a222245  = force(239)
          f1a222245  = force(240)
          f2a222245  = force(241)
          f0a222233  = force(242)
          f1a222233  = force(243)
          f2a222233  = force(244)
          f0a222224  = force(245)
          f1a222224  = force(246)
          f2a222224  = force(247)
          f0a145555  = force(248)
          f1a145555  = force(249)
          f2a145555  = force(250)
          f0a134444  = force(251)
          f1a134444  = force(252)
          f2a134444  = force(253)
          f0a133444  = force(254)
          f1a133444  = force(255)
          f2a133444  = force(256)
          f0a133345  = force(257)
          f1a133345  = force(258)
          f2a133345  = force(259)
          f0a133334  = force(260)
          f1a133334  = force(261)
          f2a133334  = force(262)
          f0a133333  = force(263)
          f1a133333  = force(264)
          f2a133333  = force(265)
          f0a124555  = force(266)
          f1a124555  = force(267)
          f2a124555  = force(268)
          f0a124455  = force(269)
          f1a124455  = force(270)
          f2a124455  = force(271)
          f0a123455  = force(272)
          f1a123455  = force(273)
          f2a123455  = force(274)
          f0a123345  = force(275)
          f1a123345  = force(276)
          f2a123345  = force(277)
          f0a113555  = force(278)
          f1a113555  = force(279)
          f2a113555  = force(280)
          f0a113345  = force(281)
          f1a113345  = force(282)
          f2a113345  = force(283)
          f0a112355  = force(284)
          f1a112355  = force(285)
          f2a112355  = force(286)
          f0a112335  = force(287)
          f1a112335  = force(288)
          f2a112335  = force(289)
          f0a112233  = force(290)
          f1a112233  = force(291)
          f2a112233  = force(292)
          f0a111444  = force(293)
          f1a111444  = force(294)
          f2a111444  = force(295)
          f0a111234  = force(296)
          f1a111234  = force(297)
          f2a111234  = force(298)
          f0a111233  = force(299)
          f1a111233  = force(300)
          f2a111233  = force(301)
          f0a111123  = force(302)
          f1a111123  = force(303)
          f2a111123  = force(304)
  
          fea555555= f0a555555 + f1a555555*coro+ f2a555555*coro**2
          fea444444= f0a444444 + f1a444444*coro+ f2a444444*coro**2
          fea335555= f0a335555 + f1a335555*coro+ f2a335555*coro**2
          fea334455= f0a334455 + f1a334455*coro+ f2a334455*coro**2
          fea334445= f0a334445 + f1a334445*coro+ f2a334445*coro**2
          fea333555= f0a333555 + f1a333555*coro+ f2a333555*coro**2
          fea333333= f0a333333 + f1a333333*coro+ f2a333333*coro**2
          fea244555= f0a244555 + f1a244555*coro+ f2a244555*coro**2
          fea244455= f0a244455 + f1a244455*coro+ f2a244455*coro**2
          fea233445= f0a233445 + f1a233445*coro+ f2a233445*coro**2
          fea233444= f0a233444 + f1a233444*coro+ f2a233444*coro**2
          fea233345= f0a233345 + f1a233345*coro+ f2a233345*coro**2
          fea233344= f0a233344 + f1a233344*coro+ f2a233344*coro**2
          fea233335= f0a233335 + f1a233335*coro+ f2a233335*coro**2
          fea223355= f0a223355 + f1a223355*coro+ f2a223355*coro**2
          fea222335= f0a222335 + f1a222335*coro+ f2a222335*coro**2
          fea222334= f0a222334 + f1a222334*coro+ f2a222334*coro**2
          fea222333= f0a222333 + f1a222333*coro+ f2a222333*coro**2
          fea222255= f0a222255 + f1a222255*coro+ f2a222255*coro**2
          fea222245= f0a222245 + f1a222245*coro+ f2a222245*coro**2
          fea222233= f0a222233 + f1a222233*coro+ f2a222233*coro**2
          fea222224= f0a222224 + f1a222224*coro+ f2a222224*coro**2
          fea145555= f0a145555 + f1a145555*coro+ f2a145555*coro**2
          fea134444= f0a134444 + f1a134444*coro+ f2a134444*coro**2
          fea133444= f0a133444 + f1a133444*coro+ f2a133444*coro**2
          fea133345= f0a133345 + f1a133345*coro+ f2a133345*coro**2
          fea133334= f0a133334 + f1a133334*coro+ f2a133334*coro**2
          fea133333= f0a133333 + f1a133333*coro+ f2a133333*coro**2
          fea124555= f0a124555 + f1a124555*coro+ f2a124555*coro**2
          fea124455= f0a124455 + f1a124455*coro+ f2a124455*coro**2
          fea123455= f0a123455 + f1a123455*coro+ f2a123455*coro**2
          fea123345= f0a123345 + f1a123345*coro+ f2a123345*coro**2
          fea113555= f0a113555 + f1a113555*coro+ f2a113555*coro**2
          fea113345= f0a113345 + f1a113345*coro+ f2a113345*coro**2
          fea112355= f0a112355 + f1a112355*coro+ f2a112355*coro**2
          fea112335= f0a112335 + f1a112335*coro+ f2a112335*coro**2
          fea112233= f0a112233 + f1a112233*coro+ f2a112233*coro**2
          fea111444= f0a111444 + f1a111444*coro+ f2a111444*coro**2
          fea111234= f0a111234 + f1a111234*coro+ f2a111234*coro**2
          fea111233= f0a111233 + f1a111233*coro+ f2a111233*coro**2
          fea111123= f0a111123 + f1a111123*coro+ f2a111123*coro**2
          !
          s3 = (y2**3*y4**3*sqrt(3.0_ark)-y2**3*y4**2*y5+y3**3*y4**2*y5-&
          5.0_ark/3.0_ark*y2**3*y4*y5**2*sqrt(3.0_ark)+y3**3*y4**3*sqrt(3.0_ark)-5.0_ark/3.0_ark*y3**&
          3*y4*y5**2*sqrt(3.0_ark)-y2**3*y5**3+y3**3*y5**3-8.0_ark/3.0_ark*y1**3*y4*y5**2*sqrt(3.0_ark))*fea333555+&
          (y1**4*y5**2*sqrt(3.0_ark)/2.0_ark+y2**4*y4*y5+y2**4*y4**2*sqrt(3.0_ark)/3.0_ark+&
          y3**4*y4**2*sqrt(3.0_ark)/3.0_ark-y3**4*y4&
          *y5-y1**4*y4**2*sqrt(3.0_ark)/6.0_ark)*fea222245+(y1*y3**5+y1*y2**5+y2**&
          5*y3+y1**5*y3+y1**5*y2+y2*y3**5)*fea133333+(y1**4*y3*y4-2.0_ark*y2**4&
          *y3*y4+y1**4*y2*y4+y1*y2**4*y5*sqrt(3.0_ark)+y1*y3**4*y4-2.0_ark*y2*y3**&
          4*y4+y1**4*y2*y5*sqrt(3.0_ark)-y1*y3**4*y5*sqrt(3.0_ark)-y1**4*y3*y5*sqrt(3.0_ark)+&
          y1*y2**4*y4)*fea133334+(-y1*y2*y3*y4**3/3.0_ark+y1*y2*y3*y4*y5**2)*fea123455 
          
          s4 = s3+(2.0_ark/3.0_ark*sqrt(3.0_ark)*y1*y2**2*y3**2*y4-y1**2*y2**2*y3*y5-&
          sqrt(3.0_ark)*y1**2*y2**2*y3*y4/3.0_ark+y1**2*y2*y3**2*y5-&
          sqrt(3.0_ark)*y1**2*y2*y3**2*y4/3.0_ark)*fea112335+(y1*y2**2*y3*y5**2+y1*y2*y3**2*y5**2+&
          y1*y2*y3**2*y4**2+y1*y2**2*y3*y4**2+y1**2*y2*y3*y4**2+y1**2*y2*y3*y5**2)*fea112355 
          
          s2 = s4+(y2**3*y3**2*y5-y1**3*y2**2*y5/2.0_ark-y1**2*y3**3*y5/2.0_ark-y2&
          **2*y3**3*y5+y1**3*y2**2*y4*sqrt(3.0_ark)/2.0_ark-y1**2*y2**3*y4*sqrt(3.0_ark)/2.0_ark+&
          y1**3*y3**2*y5/2.0_ark+y1**2*y2**3*y5/2.0_ark+y1**3*y3**2*y4*sqrt(3.0_ark)/2.0_ark-&
          y1**2*y3**3*y4*sqrt(3.0_ark)/2.0_ark)*fea222335+(-y1**2*y2&
          **2*y5**2*sqrt(3.0_ark)/2.0_ark-y1**2*y3**2*y5**2*sqrt(3.0_ark)/2.0_ark-y1**2*&
          y2**2*y4**2*sqrt(3.0_ark)/6.0_ark-y1**2*y2**2*y4*y5-2.0_ark/3.0_ark*y2**2*y3**&
          2*y4**2*sqrt(3.0_ark)+y1**2*y3**2*y4*y5-y1**2*y3**2*y4**2*sqrt(3.0_ark)/&
          6.0_ark)*fea113345+(y2**2*y3**2*y5**2+y2**2*y3**2*y4**2+y1**2*y2**2*y5**2+&
          y1**2*y3**2*y4**2+y1**2*y3**2*y5**2+y1**2*y2**2*y4**2)*fea223355 
          
          s3 = s2+(y1*y2*y3**2*y4**2*sqrt(3.0_ark)/6.0_ark+y1*y2*y3**2*y4*y5+y1*y2&
          *y3**2*y5**2*sqrt(3.0_ark)/2.0_ark+2.0_ark/3.0_ark*y1**2*y2*y3*y4**2*sqrt(3.0_ark&
          )-y1*y2**2*y3*y4*y5+y1*y2**2*y3*y4**2*sqrt(3.0_ark)/6.0_ark+y1*y2**2*y3*&
          y5**2*sqrt(3.0_ark)/2.0_ark)*fea123345+(-y1**3*y2**2*y5*sqrt(3.0_ark)/2.0_ark-&
          y1**3*y2**2*y4/2.0_ark-y1**3*y3**2*y4/2.0_ark-y1**2*y2**3*y4/2.0_ark+y1**3*&
          y3**2*y5*sqrt(3.0_ark)/2.0_ark-y1**2*y3**3*y4/2.0_ark+y2**3*y3**2*y4-y1**2*&
          y2**3*y5*sqrt(3.0_ark)/2.0_ark+y2**2*y3**3*y4+y1**2*y3**3*y5*sqrt(3.0_ark)/&
          2.0_ark)*fea222334+(3.0_ark*y3**2*y4**4+5.0_ark/2.0_ark*y1**2*y5**4+y2**2*y5**&
          4+3.0_ark*y2**2*y4**4-4.0_ark*y3**2*y4*y5**3*sqrt(3.0_ark)+y3**2*y5**4+9.0_ark&
          *y1**2*y4**2*y5**2-3.0_ark/2.0_ark*y1**2*y4**4+4.0_ark*y2**2*y4*y5**3*sqrt(&
          3.0_ark))*fea335555+(y1**3*y2**3+y1**3*y3**3+y2**3*y3**3)*fea222333 
          
          s4 = s3+(y3*y4**5/5.0_ark-y2*y4**4*y5*sqrt(3.0_ark)/2.0_ark-2.0_ark/5.0_ark*y1*y4&
          **5-2.0_ark*y1*y4**3*y5**2-3.0_ark/10.0_ark*y2*y5**5*sqrt(3.0_ark)+y3*y4**3*y5&
          **2+y3*y4**4*y5*sqrt(3.0_ark)/2.0_ark+y2*y4**3*y5**2+3.0_ark/10.0_ark*y3*y5**5&
          *sqrt(3.0_ark)+y2*y4**5/5.0_ark)*fea244455+(y2**5*y4-2.0_ark*y1**5*y4-sqrt(&
          3.0_ark)*y2**5*y5+y3**5*y4+sqrt(3.0_ark)*y3**5*y5)*fea222224 
          
          s5 = s4+(-y3*y5**5*sqrt(3.0_ark)/5.0_ark+y2*y5**5*sqrt(3.0_ark)/5.0_ark+y1*y4*&
          y5**4-7.0_ark/15.0_ark*y2*y4**5+y2*y4**4*y5*sqrt(3.0_ark)/3.0_ark-y3*y4**4*y5*&
          sqrt(3.0_ark)/3.0_ark+y3*y4*y5**4+y2*y4*y5**4+2.0_ark*y1*y4**3*y5**2-7.0_ark/15.0_ark*y3*y4**5-&
          y1*y4**5/15.0_ark)*fea145555 
          
          s1 = s5+(-sqrt(3.0_ark)*y1*y2*y3**3*y5/2.0_ark+y1**3*y2*y3*y4+sqrt(3.0_ark)&
          *y1*y2**3*y3*y5/2.0_ark-y1*y2**3*y3*y4/2.0_ark-y1*y2*y3**3*y4/2.0_ark)*fea111234+&
          (y3*y4**4*y5/3.0_ark+y3*y4**5*sqrt(3.0_ark)/18.0_ark-y2*y4**4*y5/3.0_ark&
          -y2*y4*y5**4*sqrt(3.0_ark)/2.0_ark-y3*y4**2*y5**3+2.0_ark/9.0_ark*y1*y4**5*sqrt(3.0_ark)+&
          y2*y4**5*sqrt(3.0_ark)/18.0_ark+y2*y4**2*y5**3-2.0_ark/3.0_ark*y1*y4**&
          3*y5**2*sqrt(3.0_ark)-y3*y4*y5**4*sqrt(3.0_ark)/2.0_ark)*fea244555+(y1*y2*y4**2*y5**2-&
          3.0_ark/4.0_ark*y2*y3*y4**4-y1*y2*y5**4-y1*y3*y5**4+5.0_ark/4.0_ark*y2*y3*y5**4+&
          y1*y3*y4**2*y5**2-7.0_ark/2.0_ark*y2*y3*y4**2*y5**2-2.0_ark*y1&
          *y2*y4**3*y5*sqrt(3.0_ark)+2.0_ark*y1*y3*y4**3*y5*sqrt(3.0_ark))*fea124455 
          
          s3 = s1+(y2**6+y1**6+y3**6)*fea333333+(y1*y2**4*y3+y1**4*y2*y3+y1*&
          y2*y3**4)*fea111123+fea112233*y1**2*y2**2*y3**2+(y1**4*y4**2+y2**4&
          *y4**2+y2**4*y5**2+y3**4*y4**2+y1**4*y5**2+y3**4*y5**2)*fea222255 
          s4 = s3+(3.0_ark*y1*y3*y5**4+y1*y3*y4**4+9.0_ark*y2*y3*y4**2*y5**2-3.0_ark/&
          2.0_ark*y2*y3*y5**4-4.0_ark*y1*y3*y4**3*y5*sqrt(3.0_ark)+y1*y2*y4**4+&
          4.0_ark*y1*y2*y4**3*y5*sqrt(3.0_ark)+3.0_ark*y1*y2*y5**4+5.0_ark/2.0_ark*y2*y3*y4**4)*fea134444+&
          (-y1*y3**2*y5**3*sqrt(3.0_ark)/3.0_ark-7.0_ark/3.0_ark*y1**2*y3*y4*y5&
          **2+5.0_ark/3.0_ark*y1*y2**2*y4**2*y5*sqrt(3.0_ark)-13.0_ark/3.0_ark*y2**2*y3*y4*&
          y5**2-4.0_ark/3.0_ark*y2*y3**2*y5**3*sqrt(3.0_ark)-7.0_ark/3.0_ark*y1**2*y2*y4*y5&
          **2-16.0_ark/3.0_ark*y1*y3**2*y4*y5**2+4.0_ark/3.0_ark*y1**2*y3*y4**2*y5*sqrt(&
          3.0_ark)+4.0_ark/3.0_ark*y2**2*y3*y5**3*sqrt(3.0_ark)+3.0_ark*y1**2*y2*y4**3+&
          y2*y3**2*y4**3+y1*y2**2*y5**3*sqrt(3.0_ark)/3.0_ark+y2**2*y3*y4**3-13.0_ark/3.0_ark&
          *y2*y3**2*y4*y5**2-5.0_ark/3.0_ark*y1*y3**2*y4**2*y5*sqrt(3.0_ark)-&
          4.0_ark/3.0_ark*y1**2*y2*y4**2*y5*sqrt(3.0_ark)+3.0_ark*y1**2*y3*y4**3-16.0_ark/3.0_ark*y1*&
          y2**2*y4*y5**2)*fea233444 
          
          s5 = s4+(2.0_ark*y1*y3**2*y5**3+4.0_ark*y2*y3**2*y5**3+4.0_ark*y2**2*y3*y4*&
          y5**2*sqrt(3.0_ark)-2.0_ark*y1*y2**2*y5**3+y1**2*y3*y4*y5**2*sqrt(3.0_ark)+&
          6.0_ark*y1*y3**2*y4**2*y5-6.0_ark*y1*y2**2*y4**2*y5-3.0_ark*y1**2*y3*y4**2*&
          y5+y1**2*y2*y4*y5**2*sqrt(3.0_ark)+4.0_ark*y1*y3**2*y4*y5**2*sqrt(3.0_ark)-&
          3.0_ark*y1**2*y2*y4**3*sqrt(3.0_ark)-4.0_ark*y2**2*y3*y5**3+3.0_ark*y1**2*y2*y4**2*y5-&
          y1**2*y2*y5**3+y1**2*y3*y5**3-3.0_ark*y1**2*y3*y4**3*sqrt(3.0_ark&
          )+4.0_ark*y2*y3**2*y4*y5**2*sqrt(3.0_ark)+4.0_ark*y1*y2**2*y4*y5**2*sqrt(3.0_ark))*fea113555 
          
          s2 = s5+(-2.0_ark/3.0_ark*y3**2*y4**4*sqrt(3.0_ark)-3.0_ark/2.0_ark*y1**2*y4**2*y5**2*sqrt(3.0_ark)-&
          3.0_ark/4.0_ark*y1**2*y5**4*sqrt(3.0_ark)-y2**2*y4**3*y5+&
          7.0_ark/12.0_ark*y1**2*y4**4*sqrt(3.0_ark)+y3**2*y4**3*y5+3.0_ark*y3**2*y4*y5**3&
          -2.0_ark/3.0_ark*y2**2*y4**4*sqrt(3.0_ark)-3.0_ark*y2**2*y4*y5**3)*fea334445+(&
          -3.0_ark*y1*y3*y4**3*y5+2.0_ark/3.0_ark*y1*y2*y5**4*sqrt(3.0_ark)-y1*y3*y4*y5**3+&
          2.0_ark/3.0_ark*y1*y3*y5**4*sqrt(3.0_ark)+3.0_ark*y1*y2*y4**3*y5-7.0_ark/12.0_ark&
          *y2*y3*y5**4*sqrt(3.0_ark)+3.0_ark/2.0_ark*y2*y3*y4**2*y5**2*sqrt(3.0_ark)+y1*&
          y2*y4*y5**3+3.0_ark/4.0_ark*y2*y3*y4**4*sqrt(3.0_ark))*fea124555+(2.0_ark*y3**&
          2*y4*y5**3*sqrt(3.0_ark)-7.0_ark/2.0_ark*y1**2*y4**2*y5**2+y2**2*y4**2*y5**&
          2-y2**2*y4**4-y3**2*y4**4-2.0_ark*y2**2*y4*y5**3*sqrt(3.0_ark)-3.0_ark/4.0_ark&
          *y1**2*y5**4+5.0_ark/4.0_ark*y1**2*y4**4+y3**2*y4**2*y5**2)*fea334455 
          s3 = s2+(-6.0_ark*y4**2*y5**4+9.0_ark*y4**4*y5**2+y5**6)*fea555555+(y2*y3**3*y4**2+&
          y2*y3**3*y5**2+y1*y3**3*y4**2+y1*y2**3*y4**2+y1**3*y2*y4**2+&
          y1*y2**3*y5**2+y1**3*y3*y5**2+y1**3*y3*y4**2+y1**3*y2*y5**2+y2**3*y3*y4**2+&
          y1*y3**3*y5**2+y2**3*y3*y5**2)*fea233344+(y1*y2**3*y5**2*sqrt(3.0_ark)/6.0_ark-&
          y2**3*y3*y5**2*sqrt(3.0_ark)/3.0_ark-y2*y3**3*y5**2&
          *sqrt(3.0_ark)/3.0_ark+y1**3*y2*y4*y5-y1**3*y2*y5**2*sqrt(3.0_ark)/3.0_ark-&
          y1**3*y3*y4*y5-y1**3*y3*y5**2*sqrt(3.0_ark)/3.0_ark-y1*y3**3*y4**2*sqrt(3.0_ark&
          )/2.0_ark+y1*y3**3*y5**2*sqrt(3.0_ark)/6.0_ark-y2**3*y3*y4*y5+y2*y3**3*y4*&
          y5-y1*y2**3*y4**2*sqrt(3.0_ark)/2.0_ark)*fea233345+(-3.0_ark*y2**3*y4*y5**2&
          +y3**3*y4**3-3.0_ark*y3**3*y4*y5**2-3.0_ark*y1**3*y4*y5**2+y2**3*y4**3+&
          y1**3*y4**3)*fea111444+(y1*y2**3*y3**2+y1**3*y2**2*y3+y1**2*y2**3*y3+&
          y1*y2**2*y3**3+y1**2*y2*y3**3+y1**3*y2*y3**2)*fea111233 
          
          s4 = s3+(9.0_ark*y4**2*y5**4-6.0_ark*y4**4*y5**2+y4**6)*fea444444+(-5.0_ark&
          /3.0_ark*y1*y2**2*y4**2*y5*sqrt(3.0_ark)+y1*y2**2*y4**3-4.0_ark/3.0_ark*y1**2*&
          y3*y4**2*y5*sqrt(3.0_ark)-2.0_ark*y1**2*y2*y4**3-y1*y2**2*y5**3*sqrt(3.0_ark&
          )/3.0_ark+4.0_ark/3.0_ark*y2**2*y3*y4*y5**2-4.0_ark/3.0_ark*y2**2*y3*y5**3*sqrt(&
          3.0_ark)-2.0_ark*y1**2*y3*y4**3+7.0_ark/3.0_ark*y1*y2**2*y4*y5**2-2.0_ark/3.0_ark*y1&
          **2*y3*y4*y5**2+y1*y3**2*y4**3+4.0_ark/3.0_ark*y2*y3**2*y5**3*sqrt(3.0_ark)&
          +y1*y3**2*y5**3*sqrt(3.0_ark)/3.0_ark+4.0_ark/3.0_ark*y1**2*y2*y4**2*y5*sqrt(3.0_ark)&
          +4.0_ark/3.0_ark*y2*y3**2*y4*y5**2+5.0_ark/3.0_ark*y1*y3**2*y4**2*y5*sqrt(&
          3.0_ark)-2.0_ark/3.0_ark*y1**2*y2*y4*y5**2+7.0_ark/3.0_ark*y1*y3**2*y4*y5**2)*fea133444 
          
          s5 = s4+(-y1**3*y2*y4*y5+2.0_ark/3.0_ark*y2**3*y3*y5**2*sqrt(3.0_ark)+y1*y3&
          **3*y4**2*sqrt(3.0_ark)/2.0_ark+y1**3*y3*y4**2*sqrt(3.0_ark)/2.0_ark+y1**3*y3*&
          y5**2*sqrt(3.0_ark)/6.0_ark+y1**3*y2*y5**2*sqrt(3.0_ark)/6.0_ark+y1**3*y3*y4*y5+&
          y1*y2**3*y5**2*sqrt(3.0_ark)/6.0_ark+y1**3*y2*y4**2*sqrt(3.0_ark)/2.0_ark+&
          2.0_ark/3.0_ark*y2*y3**3*y5**2*sqrt(3.0_ark)-y1*y2**3*y4*y5+y1*y2**3*y4**2*sqrt(3.0_ark)/2.0_ark+&
          y1*y3**3*y5**2*sqrt(3.0_ark)/6.0_ark+y1*y3**3*y4*y5)*fea133345 
  
          v6 = s5+(-y2**2*y3*y4**2*y5+y1**2*y3*y4*y5**2*sqrt(3.0_ark)/3.0_ark+y2*y3**2*y4**2*y5+&
          y2*y3**2*y5**3-y1*y2**2*y5**3+4.0_ark/3.0_ark*y2**2*y3*y4*&
          y5**2*sqrt(3.0_ark)+4.0_ark/3.0_ark*y2*y3**2*y4*y5**2*sqrt(3.0_ark)-y1*y2**2*y4**2*y5+&
          4.0_ark/3.0_ark*y1*y3**2*y4*y5**2*sqrt(3.0_ark)-y2**2*y3*y5**3+y1*y3**2*y5**3+&
          y1**2*y2*y4*y5**2*sqrt(3.0_ark)/3.0_ark-y1**2*y2*y4**3*sqrt(3.0_ark)&
          +y1*y3**2*y4**2*y5-y1**2*y3*y4**3*sqrt(3.0_ark)+4.0_ark/3.0_ark*y1*y2**&
          2*y4*y5**2*sqrt(3.0_ark))*fea233445+(y2*y3**4*y4*sqrt(3.0_ark)-y1**4*y2*&
          y5+y2**4*y3*y4*sqrt(3.0_ark)-y1**4*y3*y4*sqrt(3.0_ark)+y2*y3**4*y5-2.0_ark*&
          y1*y2**4*y5+2.0_ark*y1*y3**4*y5-y1**4*y2*y4*sqrt(3.0_ark)+y1**4*y3*y5-y2&
          **4*y3*y5)*fea233335+(y2**2*y3**4+y1**4*y3**2+y1**2*y2**4+y2**4*y3&
          **2+y1**2*y3**4+y1**4*y2**2)*fea222233
          !
        endif
        !
        f =  v0+v1+v2+v3+v4+v5+v6
  
  
    end function poten_xy3_morbid_10
  
  end module pot_user
  
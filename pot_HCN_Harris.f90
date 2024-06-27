!
!  This unit is for a user defined potential 
!
module pot_user
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xyz

  implicit none

  public MLdipole,MLpoten,ML_MEP,MLpoten_name

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  !
  integer(ik) :: npol(5),order(5,5)
  real(ark)   :: gam1(5,5),x0(5,5),r0(5,5),beta(5,5)
  integer(ik) :: np(5)
  !
  contains
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
     write(out,"(a,a,a,a)") 'Wrong Potential ',trim(name),'; should be ',trim(poten_name)
   endif
   !
   write(out,"(a,a)") '  Using USER-type PES ',trim(poten_name)
   !
 end subroutine MLpoten_name
 !


 !
 ! Defining potential energy function (built for SO2)

 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f
   !
   f = MLpoten_HCN_Harris(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !


 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
    !
    real(ark) :: xyz_(natoms,3), r1, r2,  n1(3), n2(3), n3(3), tmat(3,3), CN_CM(3),&
                 x(natoms,3),xyz0(natoms,3),cos_theta,alpha,mu(3),u1(3),u2(3),u3(3),bigr,smallr

    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      xyz0 = MLloc2pqr_xyz(local)
      !
      x(1,:) = xyz0(2,:) - xyz0(1,:)
      x(2,:) = xyz0(3,:) - xyz0(1,:)
      !
    else
      !
      x(1,:) = xyz(2,:) - xyz(1,:)
      x(2,:) = xyz(3,:) - xyz(1,:)
      !
    endif
    !
    r1 = sqrt(sum(x(1,:)**2))
    r2 = sqrt(sum(x(2,:)**2))
    !
    n1 = x(1,:) / r1
    n2 = x(2,:) / r2
    !
    alpha = acos(sum(n1*n2))
    !
    cos_theta = sum(n1*n2)
    !
    bigr  = r1/bohr
    smallr = r2/bohr
    !
    mu(1)=x_comp_dpl(cos_theta, bigr, smallr)
    mu(3)=z_comp_dpl(cos_theta, bigr, smallr)
    !
    mu(2) = 0
    !
    CN_CM = (x(2,:)*molec%AtomMasses(3))/(sum(molec%AtomMasses(:)))
    !
    u3 = x(1,:)-CN_CM
    !
    u3 = u3 / sqrt(sum(u3(:)**2))
    !
    u2 = (/0._ark,1.0_ark,0._ark/)
    !
    !u2 = n2
    !u2 = u2 / sqrt(sum(u2(:)**2))
    !
    u1 = MLvector_product(u2,u3)
    !
    tmat(1, :) = u1
    tmat(2, :) = u2
    tmat(3, :) = u3
    !
    f(1:3) = matmul(mu,tmat)
    !
    f = f/0.393430_ark ! to debye
    !
  end subroutine MLdipole





! VQZANO+ PES second version. 25th Feb 2004.
!
!    This is the VQZANO1 potential with 3 point mappings,
!      HNC morphing, relativistic (relcor2) and 
!      adiabatic (DBOCcor) corrections OF:-
!
! T. van Mourik, G. J. Harris, O. L. Polyansky, J. Tennyson,
! A. G. Csaszar and P. J. Knowles, J. Chem. Phys. 115, 3706 (2001).
!
! HISTORY OF CODE
!
!  14th Sept 2000  Origonal version by GJH 
!  16th Mar  2001  Revised version  by GJH
!  25th Feb  2004  Corrected for coefficient truncation errors and
!                  resulting stationary point mapping errors by GJH.
!
! ===============================================================
! ***************************************************************
! ===============================================================

! *************************************************

!subroutine MLpoten_HCN_Harris(pot, smallr, bigR, xgamma)

function MLpoten_HCN_Harris(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f

!
! Jacobi coordinates are used.
! Units of lengths in a0 and energy in Hartree.
!
! INPUT
!
! smallr is the C to N bond length.

! bigr is the H to CN center of mass distance.

! xgamma is the cosine of the angle between the bigr and smallr,
!      an angle of 0  is linear HCN
!      and an angle of Pi is CNH
!
! OUTPUT
!
! pot the potential energy in Hartree
!
! HCN minimum is at:-
! R=3.1855
! r=2.1785

      real(ark) ::  coeff(252), bigReC(4), smallreC(3),BsmallrC(1), BbigRC(3), R_co(5), m_co(9)
      real(ark) :: CN_CM(3),xyz0(natoms,3),x(2,3),r1,r2,smallr,bigr,gamma,xgamma,n1(3), n2(3),pot
      real(ark) :: EH,sr,br,cosgamma,rmin,tempr,shift1,rlebigr,rlesmallr,bigmorse,smallmorse,betabigr,betasmallr
      real(ark) :: apolynom,temp,gamorg,escale,aprxrbig,rcor,DBOCcor
      integer(ik) :: icoefref, i,j,k


      !
      ! xyz are undefined for the local case
      if (all(abs(xyz)<small_)) then 
        !
        xyz0 = MLloc2pqr_xyz(local)
        !
        x(1,:) = xyz0(2,:) - xyz0(1,:)
        x(2,:) = xyz0(3,:) - xyz0(1,:)
        !
      else
        !
        x(1,:) = xyz(2,:) - xyz(1,:)
        x(2,:) = xyz(3,:) - xyz(1,:)
        !
      endif
      !
      CN_CM = (x(1,:)*molec%AtomMasses(2))/(sum(molec%AtomMasses(:)))
      !
      r1 = sqrt(sum(x(1,:)**2))
      smallr = r1/bohr
      !
      r2 = sqrt(sum((x(2,:)-CN_CM(:))**2))
      bigr  = r2/bohr
      !
      n1 = x(1,:) / r1
      n2 = (x(2,:)-CN_CM(:)) / r2
      !
      gamma = acos(sum(n1*n2))
      !
      xgamma = cos(gamma)
      !
      coeff(1)=0.14940421045032473d0
      coeff(2)=0.217900654326183458d0
      coeff(3)=-0.327392028888569882d0
      coeff(4)=0.835149682488571383d0
      coeff(5)=-0.661914417213925382d0
      coeff(6)=0.459971320763026929d0
      coeff(7)=-0.273866303732020775d0
      coeff(8)=-4.44011996836376689d0
      coeff(9)=11.092776718502913d0
      coeff(10)=-13.9524491007917121d0
      coeff(11)=11.012961975027332d0
      coeff(12)=-7.26313371323174033d0
      coeff(13)=2.22448398431496846d0
      coeff(14)=-0.567318455569261296d0
      coeff(15)=-11.9640827597094234d0
      coeff(16)=32.8653804821333301d0
      coeff(17)=-34.0404066146125822d0
      coeff(18)=26.4703151786085983d0
      coeff(19)=-10.2448811590431805d0
      coeff(20)=4.32881481015119248d0
      coeff(21)=0.630240096316927718d0
      coeff(22)=-12.8165844409942811d0
      coeff(23)=35.0125449481178856d0
      coeff(24)=-27.870513172896129d0
      coeff(25)=24.1683798317469438d0
      coeff(26)=-3.07770877307937299d0
      coeff(27)=2.18326384079993883d0
      coeff(28)=2.46114066751500652d0
      coeff(29)=-5.09487215547354142d0
      coeff(30)=17.597963501201852d0
      coeff(31)=-7.99277500160802139d0
      coeff(32)=10.4884249553040064d0
      coeff(33)=2.66301409632960294d0
      coeff(34)=0.996764643446198058d0
      coeff(35)=2.09946959256891876d0
      coeff(36)=-0.703740212541100076d0
      coeff(37)=3.27905614092978895d0
      coeff(38)=-0.306076218649636794d0
      coeff(39)=1.87568055546004295d0
      coeff(40)=0.933946587021956015d0
      coeff(41)=0.278425171895907322d0
      coeff(42)=0.65792221319690957d0
      coeff(43)=4.79670129410093857d0
      coeff(44)=-13.1800600884948861d0
      coeff(45)=13.4840914735171191d0
      coeff(46)=-11.9580853105867104d0
      coeff(47)=6.42577079074903925d0
      coeff(48)=-2.49123707146380389d0
      coeff(49)=0.607472235122679231d0
      coeff(50)=-7.89849133799544658d0
      coeff(51)=25.1589820523135178d0
      coeff(52)=-22.2888265767286886d0
      coeff(53)=23.1929301464181875d0
      coeff(54)=-7.30009638305873811d0
      coeff(55)=4.1865822852996565d0
      coeff(56)=1.46836639793191452d0
      coeff(57)=-50.5330802775276632d0
      coeff(58)=114.406187033876383d0
      coeff(59)=-133.041030039099562d0
      coeff(60)=80.6240317074125061d0
      coeff(61)=-41.4342872689681454d0
      coeff(62)=5.32954162529848163d0
      coeff(63)=-1.05787093721719226d0
      coeff(64)=-60.2851896270153013d0
      coeff(65)=108.032425862793044d0
      coeff(66)=-144.568192050097044d0
      coeff(67)=52.2262947988976505d0
      coeff(68)=-38.5488354429517188d0
      coeff(69)=-6.43052804920155153d0
      coeff(70)=-0.247080415183886974d0
      coeff(71)=-29.720097128681807d0
      coeff(72)=24.865213190357119d0
      coeff(73)=-68.4487383668746743d0
      coeff(74)=-0.0251946841118654547d0
      coeff(75)=-24.2008817437959262d0
      coeff(76)=-9.85201778163083597d0
      coeff(77)=1.55839015596968458d0
      coeff(78)=-5.0047676699041352d0
      coeff(79)=-3.43858035544414142d0
      coeff(80)=-12.2039861533350778d0
      coeff(81)=-2.25814109580532698d0
      coeff(82)=-8.13283230977259383d0
      coeff(83)=-3.99117790275047103d0
      coeff(84)=1.25279897438148715d0
      coeff(85)=20.8036047668437169d0
      coeff(86)=-47.5279621339239009d0
      coeff(87)=60.2300636622413646d0
      coeff(88)=-43.9793450716278635d0
      coeff(89)=27.1038177803201726d0
      coeff(90)=-9.40221067114798413d0
      coeff(91)=3.13434570006633654d0
      coeff(92)=-11.1163350222532166d0
      coeff(93)=11.7792174719772819d0
      coeff(94)=-36.6636632117557018d0
      coeff(95)=10.9186893345684749d0
      coeff(96)=-20.3231652982557202d0
      coeff(97)=-0.113407046159886923d0
      coeff(98)=-3.91167981684250297d0
      coeff(99)=-108.38148996906446d0
      coeff(100)=249.010701131692587d0
      coeff(101)=-284.970366750339567d0
      coeff(102)=185.408326698166527d0
      coeff(103)=-109.874292113305107d0
      coeff(104)=33.7297262519754285d0
      coeff(105)=-19.1135372402709026d0
      coeff(106)=-98.1310653318854648d0
      coeff(107)=257.18609219016259d0
      coeff(108)=-250.524805629300039d0
      coeff(109)=189.957003395679738d0
      coeff(110)=-102.598875479190147d0
      coeff(111)=24.4409064516133971d0
      coeff(112)=-22.4419254439761535d0
      coeff(113)=-4.23469120278782542d0
      coeff(114)=85.5239709588279186d0
      coeff(115)=-49.9962123601848539d0
      coeff(116)=114.159219655855833d0
      coeff(117)=-45.3333409447954793d0
      coeff(118)=-1.45266945791890963d0
      coeff(119)=-10.2409789131196286d0
      coeff(120)=14.0046926224722434d0
      coeff(121)=0.228389578366991434d0
      coeff(122)=0.592139306963308471d0
      coeff(123)=39.0750420262782737d0
      coeff(124)=-9.49856667268513727d0
      coeff(125)=-6.4021519862004659d0
      coeff(126)=-0.451235768817357659d0
      coeff(127)=15.5677024289479361d0
      coeff(128)=-42.8526104991173929d0
      coeff(129)=37.5848122640195081d0
      coeff(130)=-28.3655343569062971d0
      coeff(131)=5.19581160700064373d0
      coeff(132)=-2.38014608559048141d0
      coeff(133)=-2.52584705312179486d0
      coeff(134)=-61.7257431194142384d0
      coeff(135)=195.435978585348013d0
      coeff(136)=-190.552420901440047d0
      coeff(137)=195.786024883592538d0
      coeff(138)=-92.9815589243333036d0
      coeff(139)=51.0648967354336233d0
      coeff(140)=-12.3539356971417142d0
      coeff(141)=-226.05258154251356d0
      coeff(142)=654.334735212995248d0
      coeff(143)=-562.970034173071587d0
      coeff(144)=534.327791258369423d0
      coeff(145)=-209.946943298766379d0
      coeff(146)=132.124388708316738d0
      coeff(147)=-19.7837592205861324d0
      coeff(148)=-167.856639772833117d0
      coeff(149)=567.436843943118731d0
      coeff(150)=-426.493408394709266d0
      coeff(151)=524.714925632482882d0
      coeff(152)=-125.400907528039235d0
      coeff(153)=136.816986685874855d0
      coeff(154)=-19.4639342702480096d0
      coeff(155)=6.92965254059776948d0
      coeff(156)=120.353749310575446d0
      coeff(157)=-158.21889462617053d0
      coeff(158)=287.134527967565371d0
      coeff(159)=-15.721708590422863d0
      coeff(160)=59.6370591713209051d0
      coeff(161)=-9.65069174642444399d0
      coeff(162)=37.5086851831762387d0
      coeff(163)=-17.2820167523175888d0
      coeff(164)=-43.4876566190562057d0
      coeff(165)=81.670913965896217d0
      coeff(166)=7.98616268395208555d0
      coeff(167)=3.5915178003836471d0
      coeff(168)=-1.31084971183005057d0
      coeff(169)=-8.49338556661102412d0
      coeff(170)=28.4499305965348798d0
      coeff(171)=-38.0055315627361919d0
      coeff(172)=44.7421989870746499d0
      coeff(173)=-28.5021543517553095d0
      coeff(174)=16.4652733581730212d0
      coeff(175)=-3.9702076094434707d0
      coeff(176)=-154.26878299291816d0
      coeff(177)=337.283260656237265d0
      coeff(178)=-446.478135789661496d0
      coeff(179)=295.150552411495791d0
      coeff(180)=-181.669803091897394d0
      coeff(181)=53.1698460567440208d0
      coeff(182)=-6.68633327469084964d0
      coeff(183)=-392.559366501732743d0
      coeff(184)=697.148899455349974d0
      coeff(185)=-950.292339165103396d0
      coeff(186)=453.846744674231419d0
      coeff(187)=-327.950331825650804d0
      coeff(188)=66.7420024502840966d0
      coeff(189)=-5.60423992400915793d0
      coeff(190)=-341.855063951868785d0
      coeff(191)=435.685071018598849d0
      coeff(192)=-896.920052715331572d0
      coeff(193)=209.915919841349364d0
      coeff(194)=-282.977076772405832d0
      coeff(195)=55.3632684162478713d0
      coeff(196)=-5.49938828126457749d0
      coeff(197)=-67.0029167154898878d0
      coeff(198)=15.4353939881417699d0
      coeff(199)=-487.506148168839836d0
      coeff(200)=26.6401499628514139d0
      coeff(201)=-112.609855324789777d0
      coeff(202)=23.8420323212550498d0
      coeff(203)=-3.8907626415968107d0
      coeff(204)=24.4021330289898155d0
      coeff(205)=-50.0559916368805536d0
      coeff(206)=-143.002841315594935d0
      coeff(207)=8.15036150283440724d0
      coeff(208)=-6.25375080821558273d0
      coeff(209)=1.30917207719287929d0
      coeff(210)=-1.36035423758027847d0
      coeff(211)=-11.4223709641714344d0
      coeff(212)=30.2573863334748862d0
      coeff(213)=-38.6980457616891423d0
      coeff(214)=32.5476707701009811d0
      coeff(215)=-19.6140444249989411d0
      coeff(216)=7.70543274621848194d0
      coeff(217)=-0.259930351577392869d0
      coeff(218)=-73.8642627982317255d0
      coeff(219)=213.266441880770153d0
      coeff(220)=-199.017958243210714d0
      coeff(221)=167.152490923235094d0
      coeff(222)=-60.4156462971610874d0
      coeff(223)=14.4818513843296322d0
      coeff(224)=1.54144204733481677d0
      coeff(225)=-140.528331123614929d0
      coeff(226)=465.913406208203455d0
      coeff(227)=-285.722301436150354d0
      coeff(228)=260.989603824715535d0
      coeff(229)=-68.6547591363068364d0
      coeff(230)=14.6775390967586618d0
      coeff(231)=5.32361174537634982d0
      coeff(232)=-67.1037226424881053d0
      coeff(233)=467.555923121107807d0
      coeff(234)=-130.191942382774234d0
      coeff(235)=188.706186094361021d0
      coeff(236)=-42.6384789693953906d0
      coeff(237)=15.784078931513617d0
      coeff(238)=6.04338941414585238d0
      coeff(239)=49.314791519219286d0
      coeff(240)=234.004630307177036d0
      coeff(241)=-28.3706394384382075d0
      coeff(242)=71.6637535185337326d0
      coeff(243)=-9.17026849397892469d0
      coeff(244)=9.93611164655507751d0
      coeff(245)=2.6669742312744035d0
      coeff(246)=40.506711919631071d0
      coeff(247)=51.9245838131337154d0
      coeff(248)=-16.4333079845281006d0
      coeff(249)=11.1712491707103627d0
      coeff(250)=3.95529588005568118d0
      coeff(251)=1.78716569984616024d0
      coeff(252)=0.0708944545076405778d0


      bigReC(1)=3.3088409915851508813d0
      bigReC(2)=-2.4563416322017199711d0
      bigReC(3)=-1.0989512458800640982d0
      bigReC(4)=0.8804973060435498100d0

      BbigRC(1)=0.83389452953072795705d0
      BbigRC(2)=-0.063550404733614548891d0
      BbigRC(3)=0.28431597260435503838d0

      smallreC(1)=2.2535051699746539988d0
      smallreC(2)=-0.89102535433330598558d0
      smallreC(3)=2.228022014935958417d0

      BsmallrC(1)=1.6606332738405520377d0
      pot = 0.0
      icoefref = 1

! LEI PATH COEFFICENTS

      R_co(1) = 2.06021279535746d0
      R_co(2) = 0.206149293036705d0
      R_co(3) = 1.38850960996447d0
      R_co(4) = -0.0678082821258994d0
      R_co(5) = -0.412935160444408d0

! HNC R MOPRHING COEFFIECNTS

      m_co(1) = -0.000303673553151895d0
      m_co(2) = 0.00435465925587477d0
      m_co(3) = 0.00402462914869603d0
      m_co(4) = 0.0139187573949336d0
      m_co(5) = 0.0416697727891190d0
      m_co(6) = 0.0264801392783952d0
      m_co(7) = -0.0216034082004957d0
      m_co(8) = -0.0493168970474053d0
      m_co(9) = -0.0297172762227854d0

      EH = 219474.624_ark ! Hartree to cm-1 conversion factor

      !gamma = acos(xgamma)

      sr=smallr
      br=bigr

! 3 POINT MORPHING, coordinate mappings

      sr = smallr + 1.65e-3_ark + (1.80738759962837e-3_ark * ((sin(gamma)+&
          (0.131325005d0*sin(gamma*2.0_ark)))**4))

      br = bigR + 2.1656e-3_ark - (1.4252441453928e-2_ark * ((sin(gamma)+&
         (0.131325005_ark*sin(gamma*2.0_ark)))**4))

      gamma = gamma+ (1.05737198462503e-2_ark*& 
         ((sin(gamma)+(0.131325005_ark*sin(gamma*2.0_ark)))**16))

      cosgamma = cos(gamma)

      ! BEGIN HNC R MORPHING

      Rmin = 0
      do j=0, 4
        Rmin = Rmin + (R_co(j+1)*(cosgamma**j))
      enddo

      tempR = bR-Rmin
      shift1=0

      do i=0, 2
         do j=0, 2

            shift1 = shift1 + (m_co(j+(i*3)+1) * (tempR**i) * (cosgamma**j)) 
         enddo
      enddo

      bR=bR+(shift1*exp(-((0.6_ark*(acos(cosgamma)-3.1415927_ark))**6)))

! NEXT THE MAIN FIT CALCULATION

      rLEbigR =  bigReC(1)+(bigReC(2)*cosgamma)&
      +(bigReC(3)*cosgamma*cosgamma)           &
      +(bigReC(4)*cosgamma*cosgamma*sr)

      rLEsmallr = smallreC(1)+(smallreC(2)*cosgamma) &
      +(smallreC(3)*cosgamma*cosgamma)

      betabigR= BbigRC(1) + (BbigRC(2) * cosgamma)  &
      + (BbigRC(3) * cosgamma * cosgamma)

      betasmallr = BsmallrC(1)

      do i=0, 5
               
         call rmorsecoord(bR, betabigR, rLEbigR, i, bigmorse)

         do j=0, 5

      call rmorsecoord(sr, betasmallr, rLEsmallr, j, smallmorse)

            do k=0, 6

               call plgndr(k, 0, cosgamma, apolynom)

               temp = coeff(icoefref)*bigmorse*smallmorse*apolynom

               pot = pot + temp

               icoefref=icoefref+1
            enddo
         enddo
      enddo

! CRITICAL POINT ENERGY SCALING NEXT

      gamorg=acos(xgamma)

      escale=(8.18864726379888d0 *  &
            (atan(15.0d0*(gamorg-1.2d0))+1.570796327d0)) + &
      (11.9947385784680d0 *                                &
      ((sin(gamorg)+(0.131325005d0*sin(gamorg*2.0d0)))**16.0))

      pot = pot + (escale/EH)

!     lower limits to potential (hole patches)


! general low R patch
      if (bigR .le. 1.0_rk) then
         pot=100.0_ark
      endif

!     Patch for irritating holes at higher low R
      if (cosgamma .le. 0.969d0 .and. cosgamma .gt. 0.071d0) then
!       calculate approxomate LEI
      aprxRbig = 2.113058664103_ark &
      + (cosgamma * 0.1544661511693d0) &
      + (cosgamma * cosgamma * 0.9706585547594d0)
      !
      if (bigR .LE. (aprxRbig-0.8d0)) then
         pot=100.0
      endif
      endif

! small r patch
      if (smallr .LE. 1.7) then 
       pot = 100.0
      endif

!     patch for large smallr
      if (smallr .GE. 3.5) then 
       pot = 100.0
      endif

!     patch for large R
! general upper limit
      if(bigR .GE. 7.0) then
       pot = 100.0
      endif

! angular conditions upper limit patchs
      if (bigR .GE. 6.5 .AND. cosgamma .le. 0.35) then
       pot = 100.0
      endif
      if (bigR .GE. 6.0 .AND. cosgamma .le. 0.3) then
       pot = 100.0
      endif
      if (bigR .GE. 5.5 .AND. cosgamma .le. 0.23173) then
       pot = 100.0
      endif

      if (bigR .GE. 5.0 .AND. cosgamma .le. 0.12945) then
       pot = 100.0
      endif

      if (bigR .GE. 4.8 .AND. cosgamma .le. -0.12944) then
       pot = 100.0
      endif

      pot = pot-3.061537d-5 ! zero potential at HCN minimum

!	call relativistic  correction

      call rel_corr(rcor, smallr, bigR, cosgamma)
      pot = pot + rcor

! 	call DBOC correction

      call DBOC_corr(DBOCcor, smallr, bigR, cosgamma)
      pot = pot + DBOCcor

!      write(6,*) "return ", pot
      !
      f  = pot*EH
      !
      return
      end function MLpoten_HCN_Harris


! ************************************************************************

      subroutine plgndr(l, m, x, pmm)

! subroutine to calculate the asociated legenre polynomials, translated drectly from a C 
! version to fortran version. The C version comes from Numerical recipies in C,  Press et al. 
      integer(ik),intent(in) :: l,m
      real(ark),intent(in) :: x
      real(ark),intent(out) :: pmm
      real(ark) :: fact, pll, pmmp1, somx2,temp,fm,fll
      integer(ik) :: i, ll

      if (m .LT. 0 .OR. m .GT. l .OR. (x*x) .GT. 1.0) then
         write(out,*)  "!ERROR! bad arguments in routine plgndr"
         stop "!ERROR! bad arguments in routine plgndr"
      endif

      pmm = 1.0_ark

      if (m .GT. 0) then
   
         somx2=sqrt((1.0_ark - x) * (1.0_ark + x))
         fact = 1.0_rk
         do i=1, m
            pmm = pmm*(-fact * somx2)
            fact = fact+ 2.0_ark
         enddo
      endif

      if (l .NE. m) then
         pmmp1 = x * (2.0_ark * m + 1.0_ark) * pmm

         if (l .EQ. (m + 1)) then
            pmm=pmmp1
         else
            do ll = m+2, l
           fm = m
           fll = ll
           
          pll=((x*((2.0_ark*fll)-1.0_ark)*pmmp1)-((fll+fm-1.0_ark)*pmm))
          temp = (fll-fm)
          pll = pll / temp
            pmm=pmmp1
            pmmp1=pll
            enddo
            pmm=pll
         endif
      endif

!      write(6,*) l, x, pmm

      return
      end subroutine plgndr

! ********************************************************************

      subroutine rmorsecoord(R, beta, requalib, n, rmor)

      real(ark),intent(in) :: R,beta,requalib
      real(ark),intent(out) :: rmor
      real(ark) :: temp,ratio
      integer(ik) :: n

      ratio = -beta * ((R - requalib) /requalib)
      temp = (1 - exp(ratio))

      rmor = temp**n

!      write (6,*) temp, rmor, n

      return
      end subroutine rmorsecoord

! *************************************************

      subroutine DBOC_corr(pot, smallr, bigR, cosgamma)

      real(ark),intent(out) :: pot
      real(ark),intent(in) :: smallr, bigR, cosgamma

! Jacobi coordinates are used.
! Units of lengths in a0 and energy in Hartree.
!
! adiabatic correction surface.
!
! INPUT
! smallr is the C to N bond length.

! bigr is the H to CN center of mass distance.

! cosgamma is the cosine of the angle between the bigr and smallr,
!      an angle of 0 coresponds to a linear molecule of configuration
!      HCN and an angle of Pi is CNH
!
! OUTPUT
!
! pot the adiabatic correction in Hartree
!
      real(ark) ::  coeff(60),EH,bigmorse,smallmorse,apolynom,temp
      integer(ik) :: icoef,i,j,k
      !
      EH = 219474.624_ark

      coeff(1)=61300.7135899d0
      coeff(2)=895.886162423d0
      coeff(3)=-70565.8990549d0
      coeff(4)=4773.9817209d0
      coeff(5)=8856.65590136d0
      coeff(6)=-80797.3926551d0
      coeff(7)=-1682.1616785d0
      coeff(8)=92565.9091745d0
      coeff(9)=-5602.42892746d0
      coeff(10)=-9198.09369058d0
      coeff(11)=35939.9687435d0
      coeff(12)=929.150781931d0
      coeff(13)=-40169.8574301d0
      coeff(14)=2130.09338513d0
      coeff(15)=2704.67838662d0
      coeff(16)=-5313.95221068d0
      coeff(17)=-157.566257215d0
      coeff(18)=5780.36035164d0
      coeff(19)=-260.350322598d0
      coeff(20)=-167.44015843d0
      coeff(21)=-55370.316994d0
      coeff(22)=4754.39164204d0
      coeff(23)=83030.3122802d0
      coeff(24)=-8281.26766355d0
      coeff(25)=-26969.0665816d0
      coeff(26)=74094.2408883d0
      coeff(27)=-5925.76456913d0
      coeff(28)=-110370.495912d0
      coeff(29)=10527.0804951d0
      coeff(30)=34834.6549419d0
      coeff(31)=-33001.6947725d0
      coeff(32)=2471.70897939d0
      coeff(33)=48649.3467107d0
      coeff(34)=-4437.80010694d0
      coeff(35)=-14790.4934647d0
      coeff(36)=4887.37761221d0
      coeff(37)=-345.374981145d0
      coeff(38)=-7117.25664686d0
      coeff(39)=620.415816354d0
      coeff(40)=2062.69994319d0
      coeff(41)=12554.0937029d0
      coeff(42)=-1936.70511436d0
      coeff(43)=-21736.5239009d0
      coeff(44)=2472.52194725d0
      coeff(45)=9092.30843229d0
      coeff(46)=-16808.5214793d0
      coeff(47)=2498.58479994d0
      coeff(48)=29077.903879d0
      coeff(49)=-3205.61349595d0
      coeff(50)=-12071.3921973d0
      coeff(51)=7489.67926539d0
      coeff(52)=-1077.32941536d0
      coeff(53)=-12914.6546977d0
      coeff(54)=1383.34822742d0
      coeff(55)=5304.19453882d0
      coeff(56)=-1109.80874249d0
      coeff(57)=155.135287705d0
      coeff(58)=1905.01555022d0
      coeff(59)=-198.680122409d0
      coeff(60)=-771.435653363d0


      pot = -840.6019631_ark ! averge value of correction
      icoef = 1

      do i=0, 2
               
         bigmorse=bigR**i

         do j=0, 3

         smallmorse=smallr**j

            do k=0, 4

               apolynom=cosgamma**k

               temp = coeff(icoef)*bigmorse*smallmorse*apolynom

               pot = pot + temp

               icoef=icoef+1
            enddo
         enddo
      enddo

      pot = pot / EH ! CONVERT FORM CM-1 TO HARTREE

      return
      end subroutine DBOC_corr


! ************************************************************************
 
      subroutine rel_corr(pot, smallr, bigR, cosgamma)
      
      real(ark),intent(in)  :: smallr, bigR, cosgamma
      real(ark),intent(out) :: pot
      
!
! The relativistic correction surface.
!
! Jacobi coordinates are used.
! Units of lengths in a0 and energy in Hartree.
!
! INPUT
! smallr is the C to N bond length.
!
! bigr is the H to CN center of mass distance.
!
! cosgamma is the cosine of the angle between the bigr and smallr,
!      an angle of 0 corresponds to a linear molecule of configuration
!      HCN and an angle of Pi is CNH
!
! OUTPUT
!
! pot the relativistic correction in Hartree
!

      real(ark) ::  coeff(80),bigmorse,smallmorse,apolynom,temp
      integer(ik) :: icoef,i,j,k

      coeff(1)=-0.0805585205815292d0
      coeff(2)=0.0950673303720766d0
      coeff(3)=-0.253364906126069d0
      coeff(4)=-0.084129451710227d0
      coeff(5)=-0.0745419772100667d0
      coeff(6)=0.0477227356064123d0
      coeff(7)=-0.130026022694324d0
      coeff(8)=0.335629886457965d0
      coeff(9)=0.107905706720457d0
      coeff(10)=0.0993733461559782d0
      coeff(11)=-0.0237612599129882d0
      coeff(12)=0.0595297964464055d0
      coeff(13)=-0.149509867823676d0
      coeff(14)=-0.0460734392042994d0
      coeff(15)=-0.0439643487346239d0
      coeff(16)=0.00388436986158473d0
      coeff(17)=-0.00906898001283882d0
      coeff(18)=0.0221602015192701d0
      coeff(19)=0.0065691932071385d0
      coeff(20)=0.00645621168784551d0
      coeff(21)=0.0266090487964399d0
      coeff(22)=-0.115552391578968d0
      coeff(23)=0.252922948126617d0
      coeff(24)=0.0870075200551087d0
      coeff(25)=0.0506809331899964d0
      coeff(26)=-0.041492313771614d0
      coeff(27)=0.156584563160005d0
      coeff(28)=-0.332335667002259d0
      coeff(29)=-0.1116697071046650d0
      coeff(30)=-0.06770235689405d0
      coeff(31)=0.0220056679462455d0
      coeff(32)=-0.0710958741686161d0
      coeff(33)=0.146735654274853d0
      coeff(34)=0.047768960584727d0
      coeff(35)=0.0298863818957252d0
      coeff(36)=-0.00375993578145236d0
      coeff(37)=0.0107485671294735d0
      coeff(38)=-0.021583118481192d0
      coeff(39)=-0.00682587066702194d0
      coeff(40)=-0.00437935957298068d0
      coeff(41)=-0.0131085901581377d0
      coeff(42)=0.0460835147852451d0
      coeff(43)=-0.0826635630451006d0
      coeff(44)=-0.0312833158935659d0
      coeff(45)=-0.00904915849992159d0
      coeff(46)=0.0193658221974871d0
      coeff(47)=-0.0618956123071639d0
      coeff(48)=0.107774673538991d0
      coeff(49)=0.040184308696244d0
      coeff(50)=0.0123015451000135d0
      coeff(51)=-0.00973451656341789d0
      coeff(52)=0.0278651342014908d0
      coeff(53)=-0.0471541574047327d0
      coeff(54)=-0.0172250650594388d0
      coeff(55)=-0.00546984649469844d0
      coeff(56)=0.001602231058765d0
      coeff(57)=-0.00418032677019256d0
      coeff(58)=0.00687897584843289d0
      coeff(59)=0.00246635113824095d0
      coeff(60)=0.000805450828032454d0
      coeff(61)=0.00229338173792963d0
      coeff(62)=-0.00610190383854349d0
      coeff(63)=0.00863846660156985d0
      coeff(64)=0.00393643841031911d0
      coeff(65)=0.000202814342030121d0
      coeff(66)=-0.00323246774517316d0
      coeff(67)=0.00812760212687042d0
      coeff(68)=-0.011178025856714d0
      coeff(69)=-0.00505844398525825d0
      coeff(70)=-0.000330146921869386d0
      coeff(71)=0.00154631457145571d0
      coeff(72)=-0.0036279978431523d0
      coeff(73)=0.00484658611284714d0
      coeff(74)=0.00217108938922934d0
      coeff(75)=0.000160901000351014d0
      coeff(76)=-0.000244921320434022d0
      coeff(77)=0.000539999643014766d0
      coeff(78)=-0.000700947801549237d0
      coeff(79)=-0.000311210335395826d0
      coeff(80)=-2.53086031415751d-5

      pot = 0.0472411971073_ark
      icoef = 1

      do i=0, 3
               
         bigmorse=bigR**i

         do j=0, 3

         smallmorse=smallr**j

            do k=0, 4

               apolynom=cosgamma**k

               temp = coeff(icoef)*bigmorse*smallmorse*apolynom

               pot = pot + temp

               icoef=icoef+1
            enddo
         enddo
      enddo

      return
      end subroutine rel_corr

! -------
      !subroutine pot(v0,vl,r1,r2)
      !return
      !
      !end subroutine rel_corr



  
! This program calculates the dipole of:-
! T. van Mourik, G. J. Harris, O. L. Polyansky, J. Tennyson, 
! A. G. Csaszar and P. J. Knowles, J. Chem. Phys. 115, 3706 (2001).
! Written by by GJH March 16th 2001.
!
! 
! Fits are 3*3*5 and 3*3*7 fits to x and z respectively

!      A test main program
 
!      implicit double precision (a-h,o-z)
!      cosg=cos(1.3962634015952)
!      bigr= 2.16667
!      smallr=2.13333 
!      icomp=0

!      should then give dipd=0.499620862d0 for icomp=1
!      for icomp=0 dipd=-0.00166465088d0


!      smallr=2.2 
!      bigr=3.05172d0
!      cosg=-1.0

!      do 10 i=1, 201

!      icomp=0
!      call dipd(dipole0, smallr, bigr, cosg,  icomp )
!      icomp=1
!      call dipd(dipole1, smallr, bigr, cosg,  icomp )

!      write(6,20) cosg, bigr, smallr, dipole0*2.5417662, dipole1*2.5417662
!      cosg=cosg+0.009999d0

! 10    continue
! 20    format(5d12.4)

!      end


! Jacobi coordinates are used.
! Units of lengths in a0 and dipole in au.
!
! smallr is the C to N bond length. 

! bigr is the H to CN center of mass distance.

! cosg is the cosine of the angle between the bigr and smallr, 
!      an angle of 0 corresponds to a linear molecule of configuration
!      HCN and an angle of Pi is CNH

! icomp selects the x or z component of the dipole, 
!     if icomp=0 then returns z component of dipole
!     if icomp=1 then the x component of the dipole is returned. Dipole is the output 

! Z dipole component is along C to H bond.
! Xdipole component is perpendicular to z component, H bends in xz plane.



! *************************************************************************
      function z_comp_dpl(cosg, bigr, smallr) result (f)

      implicit none
      real(ark),intent(in) :: cosg, bigr, smallr
      real(ark) :: coeffz(63)
      real(ark) :: B(3), C(3) , f , temp,exprbig,exprsmall,apolynom
      integer(ik) :: icoefref,i,j,k


! code to calculate the z component of the dipole moment.
! Standard Deviation = 0.008 D
           

      coeffz(1)=-1039.93729474448763d0
      coeffz(2)=353.933029028493396d0
      coeffz(3)=-469.273752106125326d0
      coeffz(4)=5053.36372859932299d0
      coeffz(5)=-174.791252669703362d0
      coeffz(6)=-1790.40179998532812d0
      coeffz(7)=-1737.30128936150733d0
      coeffz(8)=1036.43863494317905d0
      coeffz(9)=-351.840317729231717d0
      coeffz(10)=466.654729016545031d0
      coeffz(11)=-5035.72026003557988d0
      coeffz(12)=174.888204603534292d0
      coeffz(13)=1783.68882790358092d0
      coeffz(14)=1730.16793139302162d0
      coeffz(15)=4.13996041131474449d0
      coeffz(16)=-1.77817321839038074d0
      coeffz(17)=2.46984600203225505d0
      coeffz(18)=-20.8620301949723831d0
      coeffz(19)=0.248557957243695296d0
      coeffz(20)=7.69363017734038415d0
      coeffz(21)=7.82369848101642729d0
      coeffz(22)=-170.934815933980384d0
      coeffz(23)=1038.00304354550415d0
      coeffz(24)=350.15335220496472d0
      coeffz(25)=-1945.5046107913207d0
      coeffz(26)=-192.883669686284673d0
      coeffz(27)=646.020865093223846d0
      coeffz(28)=464.968280321535981d0
      coeffz(29)=170.495969427259945d0
      coeffz(30)=-1035.59570071766628d0
      coeffz(31)=-348.571134507640448d0
      coeffz(32)=1938.59278642478573d0
      coeffz(33)=191.925410021956726d0
      coeffz(34)=-643.632109747355521d0
      coeffz(35)=-463.047969228713786d0
      coeffz(36)=0.60555784179134213d0
      coeffz(37)=-3.61463889815628707d0
      coeffz(38)=-1.65622524995037699d0
      coeffz(39)=8.06544948486997637d0
      coeffz(40)=0.985792239355321115d0
      coeffz(41)=-2.74472229985483116d0
      coeffz(42)=-2.10370116762671985d0
      coeffz(43)=11.6158010664658819d0
      coeffz(44)=-1.34518644018171973d0
      coeffz(45)=-30.6940205962599487d0
      coeffz(46)=12.6697949751099632d0
      coeffz(47)=20.7241897553344583d0
      coeffz(48)=-4.33083732019109644d0
      coeffz(49)=-9.04791261638706548d0
      coeffz(50)=-11.5726193461914656d0
      coeffz(51)=1.34352911201902514d0
      coeffz(52)=30.5762810875374454d0
      coeffz(53)=-12.6246369560126679d0
      coeffz(54)=-20.6427638916496952d0
      coeffz(55)=4.31545198368715627d0
      coeffz(56)=9.01190335420644576d0
      coeffz(57)=-0.0494259888418515414d0
      coeffz(58)=0.00383241424113694307d0
      coeffz(59)=0.132747818418275613d0
      coeffz(60)=-0.0526642790873265692d0
      coeffz(61)=-0.0909424933648184257d0
      coeffz(62)=0.0179744240449498035d0
      coeffz(63)=0.0400004758413293306d0

      B(1) = -0.33194847405555984d0
      B(2) = 0.26846609653383304d0
      B(3) = 1.4784268008920333d0

      C(1) = 0.41155209135581288d0
      C(2) = 0.41038551093500483d0
      C(3) = 0.58913645613028966d0

      temp=0
      icoefref=1
      
      do i=0, 2
         exprbig=exp(bigR *  B(i+1))
         do j=0, 2
            exprsmall=exp(smallr * C(j+1))
            do k=0, 6
               call plgndr(k, 0, cosg, apolynom)
                  temp=temp+(coeffz(icoefref)*exprbig*exprsmall*apolynom)
               icoefref=icoefref+1
            enddo
         enddo
      enddo
      f = temp
      return
      
      end function z_comp_dpl


! ************************************************************************

      function x_comp_dpl(cosg, bigr, smallr) result (f)


      implicit none
      !
      real(ark),intent(in) :: cosg, bigr, smallr
      real(ark) :: coeffx(45)
      real(ark) :: B(3), C(3) , f, temp,exprbig,exprsmall,apolynom
      integer(ik) :: icoefref,i,j,k


      coeffx(1)=3307.26744745775115d0
      coeffx(2)=31514.9042500608293d0
      coeffx(3)=-20148.594904970726d0
      coeffx(4)=546.081750420754701d0
      coeffx(5)=7723.26526967422282d0
      coeffx(6)=-584.587960975932364d0
      coeffx(7)=-32963.3776253884116d0
      coeffx(8)=21191.7974142975256d0
      coeffx(9)=-590.053047331642632d0
      coeffx(10)=-7990.27903779592029d0
      coeffx(11)=-0.031796628743699731d0
      coeffx(12)=-0.00940828117939172576d0
      coeffx(13)=0.00955414571501609397d0
      coeffx(14)=-0.00072965350671948823d0
      coeffx(15)=-0.00340116044617136415d0
      coeffx(16)=-793.754144943871578d0
      coeffx(17)=-3486.98023111115815d0
      coeffx(18)=2421.86358491171674d0
      coeffx(19)=-146.240552443951966d0
      coeffx(20)=-908.472444995532274d0
      coeffx(21)=521.168738604204326d0
      coeffx(22)=3646.28887952151323d0
      coeffx(23)=-2546.78588888870825d0
      coeffx(24)=154.568740785581287d0
      coeffx(25)=939.567725418272133d0
      coeffx(26)=0.00359417943193839101d0
      coeffx(27)=0.00108249483506098995d0
      coeffx(28)=-0.00108048813729769784d0
      coeffx(29)=0.000123617108271128898d0
      coeffx(30)=0.000402299503589065202d0
      coeffx(31)=-2556.9099117606792d0
      coeffx(32)=-28050.9552366028232d0
      coeffx(33)=17760.1943565608219d0
      coeffx(34)=-413.178836891090924d0
      coeffx(35)=-6829.82090935623486d0
      coeffx(36)=106.994087651491977d0
      coeffx(37)=29341.1779083554478d0
      coeffx(38)=-18680.3377835709173d0
      coeffx(39)=449.374983259753014d0
      coeffx(40)=7066.18946387678389d0
      coeffx(41)=0.0282356308831904958d0
      coeffx(42)=0.0083400647971135194d0
      coeffx(43)=-0.00848293720114831184d0
      coeffx(44)=0.000614212002770376261d0
      coeffx(45)=0.00300582095992624047d0




! Standard Deviation = 0.002 D

      temp=0
      icoefref=1

      B(1) = -0.47194808781916352d0
      B(2) = -0.54677773064856535d0
      B(3) = -0.46378831330917923d0

      C(1) = 0.0060793434403400061d0
      C(2) = -0.011738771795116187d0
      C(3) = 3.9847070586924296d0

      do i=0, 2
         exprbig=exp(bigR *  B(i+1))
         do j=0, 2
            exprsmall=exp(smallr * C(j+1))
            do k=1, 5
               call plgndr(k, 1, cosg, apolynom)
               temp=temp+(coeffx(icoefref)*exprbig*exprsmall*apolynom)
               icoefref=icoefref+1
            enddo
         enddo
      enddo
      f = temp

      return
      end function x_comp_dpl




end module pot_user

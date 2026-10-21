!***********************************************************************
!
!  CH3Br_PEF.f90 — Standalone evaluation of the CH3Br potential energy
!                  function (PEF) in C3v(M) symmetry-adapted coordinates
!                  (modeled after HOPO_PEF.f90 & TROVE poten_zxy3_sym).
!
!  Default embedded PEF:
!    188-parameter (sty1542 topology) PEF fitted to 1,169 CCSD(T)/aug-cc-pVQZ
!    1D + 2D geometries of CH3Br (98 <=2D parameters varied, RMS = 1.3835 cm^-1).
!
!  Can also load any TROVE POTEN file (.trove) from the command line:
!    ./CH3Br_PEF.x [optional_poten_file.trove] [optional_grid_file.dat]
!
!***********************************************************************

module ch3br_pef_mod
  implicit none
  integer, parameter :: rk = selected_real_kind(15, 307)
  integer, parameter :: nparmax = 1000
  real(rk), parameter :: pi = 4.0_rk * atan2(1.0_rk, 1.0_rk)
  real(rk), parameter :: deg2rad = pi / 180.0_rk

  integer  :: npar_active = 183
  real(rk) :: Re_p    = 1.93400000_rk
  real(rk) :: re_ch_p = 1.08500000_rk
  real(rk) :: betae_p = 107.70000000_rk
  real(rk) :: a0_p    = 1.65000000_rk
  real(rk) :: b0_p    = 1.75000000_rk
  integer  :: ipow_p(9, nparmax) = 0
  real(rk) :: coeff_p(nparmax)   = 0.0_rk
  logical  :: initialized = .false.

contains

  subroutine init_default_ch3br_pef()
    implicit none
    if (initialized) return
    npar_active = 183
    Re_p    = 1.93400000_rk
    re_ch_p = 1.08500000_rk
    betae_p = 107.70000000_rk
    a0_p    = 1.65000000_rk
    b0_p    = 1.75000000_rk
    ipow_p(1:9,   1) = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  1) =         1.11084528_rk
    ipow_p(1:9,   2) = (/ 1, 0, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  2) =      -161.92297767_rk
    ipow_p(1:9,   3) = (/ 2, 0, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  3) =      4741.13249799_rk
    ipow_p(1:9,   4) = (/ 3, 0, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  4) =       300.27362571_rk
    ipow_p(1:9,   5) = (/ 4, 0, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  5) =       352.50326927_rk
    ipow_p(1:9,   6) = (/ 5, 0, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  6) =       176.00788504_rk
    ipow_p(1:9,   7) = (/ 6, 0, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  7) =        47.29355107_rk
    ipow_p(1:9,   8) = (/ 0, 1, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  8) =        46.37809248_rk
    ipow_p(1:9,   9) = (/ 0, 2, 0, 0, 0, 0, 0, 0, 0 /); coeff_p(  9) =     22527.91241481_rk
    ipow_p(1:9,  10) = (/ 0, 3, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 10) =     -2305.91575390_rk
    ipow_p(1:9,  11) = (/ 0, 4, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 11) =      1714.22157199_rk
    ipow_p(1:9,  12) = (/ 0, 5, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 12) =       357.04480325_rk
    ipow_p(1:9,  13) = (/ 0, 6, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 13) =        95.52345739_rk
    ipow_p(1:9,  14) = (/ 0, 0, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 14) =      -137.63761249_rk
    ipow_p(1:9,  15) = (/ 0, 0, 0, 0, 2, 0, 0, 0, 0 /); coeff_p( 15) =     11124.14577864_rk
    ipow_p(1:9,  16) = (/ 0, 0, 0, 0, 3, 0, 0, 0, 0 /); coeff_p( 16) =     -1125.39394732_rk
    ipow_p(1:9,  17) = (/ 0, 0, 0, 0, 4, 0, 0, 0, 0 /); coeff_p( 17) =       716.44309182_rk
    ipow_p(1:9,  18) = (/ 0, 0, 0, 0, 6, 0, 0, 0, 0 /); coeff_p( 18) =      -181.78093943_rk
    ipow_p(1:9,  19) = (/ 0, 0, 0, 0, 0, 0, 0, 2, 0 /); coeff_p( 19) =      3232.96888381_rk
    ipow_p(1:9,  20) = (/ 0, 0, 0, 0, 0, 0, 0, 4, 0 /); coeff_p( 20) =       -24.00940647_rk
    ipow_p(1:9,  21) = (/ 0, 0, 0, 0, 0, 0, 0, 6, 0 /); coeff_p( 21) =        10.77963391_rk
    ipow_p(1:9,  22) = (/ 1, 1, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 22) =       288.13690869_rk
    ipow_p(1:9,  23) = (/ 1, 2, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 23) =       849.27485494_rk
    ipow_p(1:9,  24) = (/ 2, 1, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 24) =      -426.68780758_rk
    ipow_p(1:9,  25) = (/ 1, 3, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 25) =       838.64236433_rk
    ipow_p(1:9,  26) = (/ 2, 2, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 26) =      -432.10804463_rk
    ipow_p(1:9,  27) = (/ 1, 4, 0, 0, 0, 0, 0, 0, 0 /); coeff_p( 27) =       621.43757287_rk
    ipow_p(1:9,  28) = (/ 1, 0, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 28) =      5551.98886124_rk
    ipow_p(1:9,  29) = (/ 1, 0, 0, 0, 2, 0, 0, 0, 0 /); coeff_p( 29) =     -6071.54771916_rk
    ipow_p(1:9,  30) = (/ 2, 0, 0, 0, 2, 0, 0, 0, 0 /); coeff_p( 30) =      -944.68817784_rk
    ipow_p(1:9,  31) = (/ 1, 0, 0, 0, 3, 0, 0, 0, 0 /); coeff_p( 31) =      3026.50337874_rk
    ipow_p(1:9,  32) = (/ 3, 0, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 32) =      -210.79866807_rk
    ipow_p(1:9,  33) = (/ 1, 0, 0, 0, 4, 0, 0, 0, 0 /); coeff_p( 33) =      -668.10640076_rk
    ipow_p(1:9,  34) = (/ 4, 0, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 34) =      -553.25979488_rk
    ipow_p(1:9,  35) = (/ 2, 0, 0, 0, 3, 0, 0, 0, 0 /); coeff_p( 35) =       305.02685647_rk
    ipow_p(1:9,  36) = (/ 1, 0, 0, 0, 0, 0, 0, 2, 0 /); coeff_p( 36) =       172.16562373_rk
    ipow_p(1:9,  37) = (/ 0, 1, 1, 0, 0, 0, 0, 0, 0 /); coeff_p( 37) =       198.23617489_rk
    ipow_p(1:9,  38) = (/ 0, 2, 1, 0, 0, 0, 0, 0, 0 /); coeff_p( 38) =       346.67631801_rk
    ipow_p(1:9,  39) = (/ 0, 1, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 39) =     -1406.29770554_rk
    ipow_p(1:9,  40) = (/ 0, 1, 0, 0, 2, 0, 0, 0, 0 /); coeff_p( 40) =     -1386.32298067_rk
    ipow_p(1:9,  41) = (/ 0, 2, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 41) =      -639.84703462_rk
    ipow_p(1:9,  42) = (/ 0, 2, 0, 0, 2, 0, 0, 0, 0 /); coeff_p( 42) =     -1923.01868539_rk
    ipow_p(1:9,  43) = (/ 0, 1, 0, 0, 3, 0, 0, 0, 0 /); coeff_p( 43) =      -754.15080158_rk
    ipow_p(1:9,  44) = (/ 0, 3, 0, 0, 2, 0, 0, 0, 0 /); coeff_p( 44) =     -1594.16234999_rk
    ipow_p(1:9,  45) = (/ 0, 2, 0, 0, 3, 0, 0, 0, 0 /); coeff_p( 45) =      -165.23860651_rk
    ipow_p(1:9,  46) = (/ 0, 1, 0, 0, 5, 0, 0, 0, 0 /); coeff_p( 46) =      1158.79042442_rk
    ipow_p(1:9,  47) = (/ 0, 3, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 47) =       -13.16131752_rk
    ipow_p(1:9,  48) = (/ 0, 1, 0, 0, 0, 1, 0, 0, 0 /); coeff_p( 48) =     -1413.76695661_rk
    ipow_p(1:9,  49) = (/ 0, 1, 0, 0, 0, 2, 0, 0, 0 /); coeff_p( 49) =      -994.02760493_rk
    ipow_p(1:9,  50) = (/ 0, 2, 0, 0, 0, 1, 0, 0, 0 /); coeff_p( 50) =     -1277.42802343_rk
    ipow_p(1:9,  51) = (/ 0, 2, 0, 0, 0, 2, 0, 0, 0 /); coeff_p( 51) =     -1105.66212064_rk
    ipow_p(1:9,  52) = (/ 0, 1, 0, 0, 0, 3, 0, 0, 0 /); coeff_p( 52) =      -322.26447567_rk
    ipow_p(1:9,  53) = (/ 0, 3, 0, 0, 0, 1, 0, 0, 0 /); coeff_p( 53) =      -511.43795256_rk
    ipow_p(1:9,  54) = (/ 0, 3, 0, 0, 0, 2, 0, 0, 0 /); coeff_p( 54) =      -481.65631704_rk
    ipow_p(1:9,  55) = (/ 0, 2, 0, 0, 0, 3, 0, 0, 0 /); coeff_p( 55) =      -306.82498623_rk
    ipow_p(1:9,  56) = (/ 0, 1, 0, 0, 0, 0, 0, 1, 0 /); coeff_p( 56) =     -1451.76332200_rk
    ipow_p(1:9,  57) = (/ 0, 1, 0, 0, 0, 0, 0, 2, 0 /); coeff_p( 57) =      -154.63871999_rk
    ipow_p(1:9,  58) = (/ 0, 3, 0, 0, 0, 0, 0, 1, 0 /); coeff_p( 58) =      -522.01017961_rk
    ipow_p(1:9,  59) = (/ 0, 1, 0, 0, 0, 0, 0, 3, 0 /); coeff_p( 59) =      -211.24089667_rk
    ipow_p(1:9,  60) = (/ 0, 2, 0, 0, 0, 0, 0, 1, 0 /); coeff_p( 60) =     -1120.27507234_rk
    ipow_p(1:9,  61) = (/ 0, 2, 0, 0, 0, 0, 0, 3, 0 /); coeff_p( 61) =       -37.24841623_rk
    ipow_p(1:9,  62) = (/ 0, 2, 0, 0, 0, 0, 0, 2, 0 /); coeff_p( 62) =      -118.77530961_rk
    ipow_p(1:9,  63) = (/ 0, 0, 1, 0, 0, 0, 0, 2, 0 /); coeff_p( 63) =     -1606.17149569_rk
    ipow_p(1:9,  64) = (/ 0, 0, 2, 0, 0, 0, 0, 2, 0 /); coeff_p( 64) =     -1486.86289132_rk
    ipow_p(1:9,  65) = (/ 0, 0, 3, 0, 0, 0, 0, 2, 0 /); coeff_p( 65) =     -1583.91136564_rk
    ipow_p(1:9,  66) = (/ 0, 0, 4, 0, 0, 0, 0, 2, 0 /); coeff_p( 66) =     -1251.55451446_rk
    ipow_p(1:9,  67) = (/ 0, 0, 0, 0, 1, 1, 0, 0, 0 /); coeff_p( 67) =      3091.44519631_rk
    ipow_p(1:9,  68) = (/ 0, 0, 0, 0, 2, 1, 0, 0, 0 /); coeff_p( 68) =      6399.95982083_rk
    ipow_p(1:9,  69) = (/ 0, 0, 0, 0, 2, 2, 0, 0, 0 /); coeff_p( 69) =      1114.09159258_rk
    ipow_p(1:9,  70) = (/ 0, 0, 0, 0, 4, 1, 0, 0, 0 /); coeff_p( 70) =       776.58015707_rk
    ipow_p(1:9,  71) = (/ 0, 0, 0, 0, 3, 3, 0, 0, 0 /); coeff_p( 71) =       952.39106563_rk
    ipow_p(1:9,  72) = (/ 0, 0, 0, 0, 1, 0, 0, 1, 0 /); coeff_p( 72) =      4261.20078636_rk
    ipow_p(1:9,  73) = (/ 0, 0, 0, 0, 1, 0, 0, 2, 0 /); coeff_p( 73) =       912.61430500_rk
    ipow_p(1:9,  74) = (/ 0, 0, 0, 0, 1, 0, 0, 3, 0 /); coeff_p( 74) =       123.65356728_rk
    ipow_p(1:9,  75) = (/ 0, 0, 0, 0, 1, 0, 0, 4, 0 /); coeff_p( 75) =       -89.86130314_rk
    ipow_p(1:9,  76) = (/ 0, 0, 0, 0, 2, 0, 0, 1, 0 /); coeff_p( 76) =      1089.39159315_rk
    ipow_p(1:9,  77) = (/ 0, 0, 0, 0, 2, 0, 0, 2, 0 /); coeff_p( 77) =       407.01529087_rk
    ipow_p(1:9,  78) = (/ 0, 0, 0, 0, 2, 0, 0, 3, 0 /); coeff_p( 78) =      -192.65137127_rk
    ipow_p(1:9,  79) = (/ 0, 0, 0, 0, 3, 0, 0, 1, 0 /); coeff_p( 79) =      -808.36860470_rk
    ipow_p(1:9,  80) = (/ 0, 0, 0, 0, 4, 0, 0, 1, 0 /); coeff_p( 80) =       264.16104016_rk
    ipow_p(1:9,  81) = (/ 0, 0, 0, 0, 5, 0, 0, 1, 0 /); coeff_p( 81) =      -379.19808039_rk
    ipow_p(1:9,  82) = (/ 0, 0, 0, 0, 0, 1, 0, 2, 0 /); coeff_p( 82) =     -1594.58685038_rk
    ipow_p(1:9,  83) = (/ 0, 0, 0, 0, 0, 2, 0, 2, 0 /); coeff_p( 83) =     -2768.20757006_rk
    ipow_p(1:9,  84) = (/ 0, 0, 0, 0, 0, 1, 0, 4, 0 /); coeff_p( 84) =       -82.36054388_rk
    ipow_p(1:9,  85) = (/ 0, 0, 0, 0, 0, 0, 0, 3, 0 /); coeff_p( 85) =      -970.99891407_rk
    ipow_p(1:9,  86) = (/ 0, 0, 0, 0, 0, 0, 0, 3, 2 /); coeff_p( 86) =       212.23166473_rk
    ipow_p(1:9,  87) = (/ 0, 0, 0, 0, 0, 0, 0, 4, 2 /); coeff_p( 87) =       -69.90408513_rk
    ipow_p(1:9,  88) = (/ 1, 1, 1, 0, 0, 0, 0, 0, 0 /); coeff_p( 88) =      1527.08611005_rk
    ipow_p(1:9,  89) = (/ 1, 1, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 89) =     -3196.62497485_rk
    ipow_p(1:9,  90) = (/ 1, 2, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 90) =     -1595.95028235_rk
    ipow_p(1:9,  91) = (/ 2, 1, 0, 0, 1, 0, 0, 0, 0 /); coeff_p( 91) =       414.89321722_rk
    ipow_p(1:9,  92) = (/ 1, 1, 0, 0, 2, 0, 0, 0, 0 /); coeff_p( 92) =       654.09088851_rk
    ipow_p(1:9,  93) = (/ 1, 2, 0, 0, 2, 0, 0, 0, 0 /); coeff_p( 93) =      1229.44297086_rk
    ipow_p(1:9,  94) = (/ 1, 1, 0, 0, 0, 1, 0, 0, 0 /); coeff_p( 94) =       533.29788702_rk
    ipow_p(1:9,  95) = (/ 1, 1, 0, 0, 0, 0, 0, 1, 0 /); coeff_p( 95) =      -153.61245100_rk
    ipow_p(1:9,  96) = (/ 1, 0, 0, 0, 1, 1, 0, 0, 0 /); coeff_p( 96) =      1563.55481206_rk
    ipow_p(1:9,  97) = (/ 2, 0, 0, 0, 1, 1, 0, 0, 0 /); coeff_p( 97) =       797.13404799_rk
    ipow_p(1:9,  98) = (/ 3, 0, 0, 0, 1, 1, 0, 0, 0 /); coeff_p( 98) =      1271.54141889_rk
    ipow_p(1:9,  99) = (/ 1, 0, 0, 0, 1, 0, 0, 1, 0 /); coeff_p( 99) =      -987.29606070_rk
    ipow_p(1:9, 100) = (/ 1, 0, 0, 0, 1, 0, 0, 2, 0 /); coeff_p(100) =      -367.41207578_rk
    ipow_p(1:9, 101) = (/ 1, 0, 0, 0, 2, 0, 0, 1, 0 /); coeff_p(101) =      -457.74667809_rk
    ipow_p(1:9, 102) = (/ 2, 0, 0, 0, 1, 0, 0, 1, 0 /); coeff_p(102) =      -911.64906641_rk
    ipow_p(1:9, 103) = (/ 0, 1, 1, 0, 1, 0, 0, 0, 0 /); coeff_p(103) =      3910.06191003_rk
    ipow_p(1:9, 104) = (/ 0, 1, 1, 0, 2, 0, 0, 0, 0 /); coeff_p(104) =      1535.82006527_rk
    ipow_p(1:9, 105) = (/ 0, 2, 1, 0, 1, 0, 0, 0, 0 /); coeff_p(105) =       949.56386642_rk
    ipow_p(1:9, 106) = (/ 0, 1, 2, 0, 1, 0, 0, 0, 0 /); coeff_p(106) =       967.21547119_rk
    ipow_p(1:9, 107) = (/ 0, 1, 1, 0, 0, 0, 1, 0, 0 /); coeff_p(107) =      2398.38497384_rk
    ipow_p(1:9, 108) = (/ 0, 1, 1, 0, 0, 0, 0, 1, 0 /); coeff_p(108) =      2979.67286505_rk
    ipow_p(1:9, 109) = (/ 0, 2, 1, 0, 0, 0, 0, 1, 0 /); coeff_p(109) =       575.58138455_rk
    ipow_p(1:9, 110) = (/ 0, 2, 1, 0, 0, 0, 0, 0, 1 /); coeff_p(110) =     -1550.22920739_rk
    ipow_p(1:9, 111) = (/ 0, 1, 0, 0, 1, 1, 0, 0, 0 /); coeff_p(111) =     -4764.90175837_rk
    ipow_p(1:9, 112) = (/ 0, 1, 0, 0, 2, 1, 0, 0, 0 /); coeff_p(112) =     -1608.41191150_rk
    ipow_p(1:9, 113) = (/ 0, 2, 0, 0, 1, 1, 0, 0, 0 /); coeff_p(113) =     -2030.37312355_rk
    ipow_p(1:9, 114) = (/ 0, 2, 0, 0, 2, 1, 0, 0, 0 /); coeff_p(114) =     -3080.02440981_rk
    ipow_p(1:9, 115) = (/ 0, 1, 0, 0, 1, 2, 0, 0, 0 /); coeff_p(115) =     -1728.61549640_rk
    ipow_p(1:9, 116) = (/ 0, 2, 0, 0, 1, 2, 0, 0, 0 /); coeff_p(116) =     -2095.27646137_rk
    ipow_p(1:9, 117) = (/ 0, 1, 0, 0, 1, 0, 0, 1, 0 /); coeff_p(117) =     -1991.04526461_rk
    ipow_p(1:9, 118) = (/ 0, 2, 0, 0, 1, 0, 0, 1, 0 /); coeff_p(118) =     -1116.80590747_rk
    ipow_p(1:9, 119) = (/ 0, 1, 0, 0, 1, 0, 0, 2, 0 /); coeff_p(119) =      -179.61175350_rk
    ipow_p(1:9, 120) = (/ 0, 3, 0, 0, 1, 0, 0, 1, 0 /); coeff_p(120) =     -1139.90739618_rk
    ipow_p(1:9, 121) = (/ 0, 1, 0, 0, 1, 0, 0, 0, 2 /); coeff_p(121) =       790.36542762_rk
    ipow_p(1:9, 122) = (/ 0, 1, 0, 0, 0, 1, 1, 0, 0 /); coeff_p(122) =       517.79200214_rk
    ipow_p(1:9, 123) = (/ 0, 1, 0, 0, 0, 2, 1, 0, 0 /); coeff_p(123) =      1248.98842999_rk
    ipow_p(1:9, 124) = (/ 0, 2, 0, 0, 0, 1, 1, 0, 0 /); coeff_p(124) =       831.60543112_rk
    ipow_p(1:9, 125) = (/ 0, 1, 0, 0, 0, 1, 0, 1, 0 /); coeff_p(125) =     -2481.27232302_rk
    ipow_p(1:9, 126) = (/ 0, 2, 0, 0, 0, 1, 0, 1, 0 /); coeff_p(126) =     -1381.83693611_rk
    ipow_p(1:9, 127) = (/ 0, 1, 0, 0, 0, 2, 0, 1, 0 /); coeff_p(127) =       999.48729838_rk
    ipow_p(1:9, 128) = (/ 0, 1, 0, 0, 0, 3, 0, 1, 0 /); coeff_p(128) =       922.25334674_rk
    ipow_p(1:9, 129) = (/ 0, 1, 0, 0, 0, 1, 0, 2, 0 /); coeff_p(129) =       266.29894429_rk
    ipow_p(1:9, 130) = (/ 0, 1, 0, 0, 0, 1, 0, 0, 1 /); coeff_p(130) =      1266.16517067_rk
    ipow_p(1:9, 131) = (/ 0, 1, 0, 0, 0, 4, 0, 0, 1 /); coeff_p(131) =     -2061.97869514_rk
    ipow_p(1:9, 132) = (/ 0, 1, 0, 0, 0, 0, 0, 1, 2 /); coeff_p(132) =      -572.05805206_rk
    ipow_p(1:9, 133) = (/ 0, 2, 0, 0, 0, 0, 0, 1, 2 /); coeff_p(133) =      -942.51690152_rk
    ipow_p(1:9, 134) = (/ 0, 0, 1, 1, 0, 0, 0, 2, 0 /); coeff_p(134) =       924.26799457_rk
    ipow_p(1:9, 135) = (/ 0, 0, 2, 1, 0, 0, 0, 2, 0 /); coeff_p(135) =       832.36811174_rk
    ipow_p(1:9, 136) = (/ 0, 0, 2, 0, 1, 0, 0, 1, 0 /); coeff_p(136) =     -2085.65937281_rk
    ipow_p(1:9, 137) = (/ 0, 0, 3, 0, 1, 0, 0, 1, 0 /); coeff_p(137) =     -1687.52312499_rk
    ipow_p(1:9, 138) = (/ 0, 0, 1, 0, 2, 0, 0, 2, 0 /); coeff_p(138) =      1228.06542303_rk
    ipow_p(1:9, 139) = (/ 0, 0, 0, 0, 1, 1, 1, 0, 0 /); coeff_p(139) =     -1208.66373493_rk
    ipow_p(1:9, 140) = (/ 0, 0, 0, 0, 2, 1, 1, 0, 0 /); coeff_p(140) =     -1806.34257401_rk
    ipow_p(1:9, 141) = (/ 0, 0, 0, 0, 3, 1, 1, 0, 0 /); coeff_p(141) =       841.99174974_rk
    ipow_p(1:9, 142) = (/ 0, 0, 0, 0, 2, 2, 0, 1, 0 /); coeff_p(142) =      1305.13606472_rk
    ipow_p(1:9, 143) = (/ 0, 0, 0, 0, 3, 1, 0, 1, 0 /); coeff_p(143) =      1044.41277095_rk
    ipow_p(1:9, 144) = (/ 0, 0, 0, 0, 1, 1, 0, 2, 0 /); coeff_p(144) =      1205.22519003_rk
    ipow_p(1:9, 145) = (/ 0, 0, 0, 0, 2, 2, 0, 2, 0 /); coeff_p(145) =      1267.52704928_rk
    ipow_p(1:9, 146) = (/ 0, 0, 0, 0, 1, 1, 0, 3, 0 /); coeff_p(146) =      -862.05936278_rk
    ipow_p(1:9, 147) = (/ 0, 0, 0, 0, 2, 1, 0, 0, 1 /); coeff_p(147) =      3235.68203073_rk
    ipow_p(1:9, 148) = (/ 0, 0, 0, 0, 1, 0, 0, 1, 2 /); coeff_p(148) =      -444.91090662_rk
    ipow_p(1:9, 149) = (/ 0, 0, 0, 0, 2, 0, 0, 1, 2 /); coeff_p(149) =     -1583.30023805_rk
    ipow_p(1:9, 150) = (/ 0, 0, 0, 0, 1, 0, 0, 2, 2 /); coeff_p(150) =       441.49467172_rk
    ipow_p(1:9, 151) = (/ 0, 0, 0, 0, 0, 3, 1, 1, 0 /); coeff_p(151) =      1825.31388183_rk
    ipow_p(1:9, 152) = (/ 0, 0, 0, 0, 0, 4, 1, 1, 0 /); coeff_p(152) =      1705.46992240_rk
    ipow_p(1:9, 153) = (/ 1, 1, 0, 0, 1, 1, 0, 0, 0 /); coeff_p(153) =     -1917.14851483_rk
    ipow_p(1:9, 154) = (/ 1, 0, 0, 0, 1, 1, 1, 0, 0 /); coeff_p(154) =      2129.92504604_rk
    ipow_p(1:9, 155) = (/ 1, 0, 0, 0, 1, 1, 0, 1, 0 /); coeff_p(155) =      3261.46357616_rk
    ipow_p(1:9, 156) = (/ 0, 1, 1, 0, 1, 1, 0, 0, 0 /); coeff_p(156) =      6759.28487368_rk
    ipow_p(1:9, 157) = (/ 0, 1, 1, 0, 0, 0, 1, 1, 0 /); coeff_p(157) =     -1875.32875477_rk
    ipow_p(1:9, 158) = (/ 0, 1, 0, 0, 1, 1, 1, 0, 0 /); coeff_p(158) =     -7136.58097308_rk
    ipow_p(1:9, 159) = (/ 0, 1, 0, 0, 1, 1, 0, 1, 0 /); coeff_p(159) =       922.02959037_rk
    ipow_p(1:9, 160) = (/ 0, 1, 0, 0, 0, 2, 0, 1, 1 /); coeff_p(160) =     -1986.26360638_rk
    ipow_p(1:9, 161) = (/ 0, 0, 0, 0, 2, 1, 0, 1, 1 /); coeff_p(161) =     -1863.99930491_rk
    ipow_p(1:9, 162) = (/ 0, 0, 0, 0, 5, 0, 0, 0, 0 /); coeff_p(162) =      -133.06088668_rk
    ipow_p(1:9, 163) = (/ 2, 0, 0, 0, 1, 0, 0, 0, 0 /); coeff_p(163) =       -20.44445112_rk
    ipow_p(1:9, 164) = (/ 5, 0, 0, 0, 1, 0, 0, 0, 0 /); coeff_p(164) =      -495.86722852_rk
    ipow_p(1:9, 165) = (/ 1, 0, 0, 0, 0, 0, 0, 4, 0 /); coeff_p(165) =       -19.68116765_rk
    ipow_p(1:9, 166) = (/ 0, 3, 1, 0, 0, 0, 0, 0, 0 /); coeff_p(166) =      -359.06380321_rk
    ipow_p(1:9, 167) = (/ 0, 4, 0, 0, 2, 0, 0, 0, 0 /); coeff_p(167) =      -839.66703532_rk
    ipow_p(1:9, 168) = (/ 0, 4, 0, 0, 0, 0, 0, 1, 0 /); coeff_p(168) =      -133.46410289_rk
    ipow_p(1:9, 169) = (/ 0, 0, 1, 0, 0, 0, 0, 4, 0 /); coeff_p(169) =       202.17417315_rk
    ipow_p(1:9, 170) = (/ 0, 0, 0, 0, 4, 2, 0, 0, 0 /); coeff_p(170) =       736.00167859_rk
    ipow_p(1:9, 171) = (/ 0, 0, 0, 0, 0, 2, 0, 4, 0 /); coeff_p(171) =      -246.22874846_rk
    ipow_p(1:9, 172) = (/ 0, 0, 0, 0, 0, 4, 0, 2, 0 /); coeff_p(172) =      -401.66866392_rk
    ipow_p(1:9, 173) = (/ 2, 0, 0, 0, 3, 1, 0, 0, 0 /); coeff_p(173) =      2192.62347837_rk
    ipow_p(1:9, 174) = (/ 2, 0, 0, 0, 0, 2, 0, 2, 0 /); coeff_p(174) =      -881.42398072_rk
    ipow_p(1:9, 175) = (/ 1, 0, 0, 0, 0, 0, 0, 1, 2 /); coeff_p(175) =        66.58244331_rk
    ipow_p(1:9, 176) = (/ 0, 1, 0, 0, 0, 0, 0, 2, 2 /); coeff_p(176) =      -360.14055280_rk
    ipow_p(1:9, 177) = (/ 0, 0, 0, 0, 1, 1, 0, 1, 0 /); coeff_p(177) =      3070.87572060_rk
    ipow_p(1:9, 178) = (/ 0, 0, 0, 0, 0, 3, 1, 2, 0 /); coeff_p(178) =     -1594.29925199_rk
    ipow_p(1:9, 179) = (/ 0, 0, 0, 0, 1, 1, 0, 2, 1 /); coeff_p(179) =     -2352.82144491_rk
    ipow_p(1:9, 180) = (/ 0, 0, 0, 0, 1, 2, 0, 2, 1 /); coeff_p(180) =      1388.31262684_rk
    ipow_p(1:9, 181) = (/ 0, 0, 0, 0, 1, 1, 0, 3, 1 /); coeff_p(181) =      -900.04612840_rk
    ipow_p(1:9, 182) = (/ 0, 0, 0, 0, 1, 1, 0, 1, 1 /); coeff_p(182) =     -3415.91973797_rk
    ipow_p(1:9, 183) = (/ 0, 0, 0, 0, 2, 2, 0, 1, 1 /); coeff_p(183) =     -2422.90322633_rk
    initialized = .true.
  end subroutine init_default_ch3br_pef

  subroutine load_trove_poten(filename)
    character(len=*), intent(in) :: filename
    character(len=256) :: line, label
    integer :: iu, ios, p(11), k
    real(rk) :: val
    open(newunit=iu, file=trim(filename), status="old", action="read")
    k = 0
    do
      read(iu, "(a)", iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:5) == "POTEN" .or. line(1:6) == "NPARAM" .or. &
          line(1:8) == "POT_TYPE" .or. line(1:5) == "COEFF" .or. &
          line(1:3) == "end" .or. line(1:3) == "END") cycle
      read(line, *) label
      if (trim(label) == "Re") then
        read(line, *) label, p(1:11), Re_p
      else if (trim(label) == "re") then
        read(line, *) label, p(1:11), re_ch_p
      else if (trim(label) == "betae") then
        read(line, *) label, p(1:11), betae_p
      else if (trim(label) == "a0") then
        read(line, *) label, p(1:11), a0_p
      else if (trim(label) == "b0") then
        read(line, *) label, p(1:11), b0_p
      else if (trim(label) == "f") then
        k = k + 1
        read(line, *) label, ipow_p(1:9, k), p(1:2), val
        coeff_p(k) = val
      end if
    end do
    close(iu)
    npar_active = k
    initialized = .true.
  end subroutine load_trove_poten

  subroutine poten_CH3Br(V, r0, r1, r2, r3, beta1_deg, beta2_deg, beta3_deg, tau12_deg, tau13_deg)
    real(rk), intent(out) :: V
    real(rk), intent(in)  :: r0, r1, r2, r3, beta1_deg, beta2_deg, beta3_deg, tau12_deg, tau13_deg
    real(rk) :: xi(9), chi(9, 6), betae_rad, tau12, tau13, tau23, term
    integer  :: ioper, m

    if (.not. initialized) call init_default_ch3br_pef()

    betae_rad = betae_p * deg2rad
    xi(1) = 1.0_rk - exp(-a0_p * (r0 - Re_p))
    xi(2) = 1.0_rk - exp(-b0_p * (r1 - re_ch_p))
    xi(3) = 1.0_rk - exp(-b0_p * (r2 - re_ch_p))
    xi(4) = 1.0_rk - exp(-b0_p * (r3 - re_ch_p))
    xi(5) = beta1_deg * deg2rad - betae_rad
    xi(6) = beta2_deg * deg2rad - betae_rad
    xi(7) = beta3_deg * deg2rad - betae_rad
    tau12 = tau12_deg * deg2rad
    tau13 = tau13_deg * deg2rad
    tau23 = 2.0_rk * pi - tau12 - tau13
    xi(8) = (2.0_rk * tau23 - tau13 - tau12) / sqrt(6.0_rk)
    xi(9) = (tau12 - tau13) / sqrt(2.0_rk)

    do ioper = 1, 6
      call sym_transform_xy3(ioper, xi, chi(:, ioper))
    end do

    V = 0.0_rk
    do m = 1, npar_active
      term = 0.0_rk
      do ioper = 1, 6
        term = term + product(chi(1:9, ioper)**ipow_p(1:9, m))
      end do
      V = V + coeff_p(m) * term
    end do
  end subroutine poten_CH3Br

  subroutine sym_transform_xy3(ioper, src, dst)
    integer,  intent(in)  :: ioper
    real(rk), intent(in)  :: src(9)
    real(rk), intent(out) :: dst(9)
    real(rk) :: repres(6, 9, 9), a, b
    a = 0.5_rk
    b = 0.5_rk * sqrt(3.0_rk)
    repres = 0.0_rk
    repres(1,1,1)=1; repres(1,2,2)=1; repres(1,3,3)=1; repres(1,4,4)=1
    repres(1,5,5)=1; repres(1,6,6)=1; repres(1,7,7)=1; repres(1,8,8)=1; repres(1,9,9)=1

    repres(2,1,1)=1; repres(2,2,4)=1; repres(2,3,2)=1; repres(2,4,3)=1
    repres(2,5,7)=1; repres(2,6,5)=1; repres(2,7,6)=1
    repres(2,8,8)=-a; repres(2,8,9)=b; repres(2,9,8)=-b; repres(2,9,9)=-a

    repres(3,1,1)=1; repres(3,2,3)=1; repres(3,3,4)=1; repres(3,4,2)=1
    repres(3,5,6)=1; repres(3,6,7)=1; repres(3,7,5)=1
    repres(3,8,8)=-a; repres(3,8,9)=-b; repres(3,9,8)=b; repres(3,9,9)=-a

    repres(4,1,1)=1; repres(4,2,2)=1; repres(4,3,4)=1; repres(4,4,3)=1
    repres(4,5,5)=1; repres(4,6,7)=1; repres(4,7,6)=1
    repres(4,8,8)=1; repres(4,9,9)=-1

    repres(5,1,1)=1; repres(5,2,4)=1; repres(5,3,3)=1; repres(5,4,2)=1
    repres(5,5,7)=1; repres(5,6,6)=1; repres(5,7,5)=1
    repres(5,8,8)=-a; repres(5,8,9)=b; repres(5,9,8)=b; repres(5,9,9)=a

    repres(6,1,1)=1; repres(6,2,3)=1; repres(6,3,2)=1; repres(6,4,4)=1
    repres(6,5,6)=1; repres(6,6,5)=1; repres(6,7,7)=1
    repres(6,8,8)=-a; repres(6,8,9)=-b; repres(6,9,8)=-b; repres(6,9,9)=a

    dst = matmul(repres(ioper, :, :), src)
  end subroutine sym_transform_xy3

end module ch3br_pef_mod

program test_ch3br_pef
  use ch3br_pef_mod
  implicit none
  character(len=256) :: arg1, arg2, line
  integer :: nargs, iu, ios, npts
  real(rk) :: c(9), e_ha, e_ab, v_calc, v_ref, diff, ssq

  nargs = command_argument_count()
  if (nargs >= 1) then
    call get_command_argument(1, arg1)
    call load_trove_poten(arg1)
    write(*, "(a,a,a,i0,a)") "Loaded TROVE POTEN from ", trim(arg1), " (", npar_active, " linear terms)"
  else
    call init_default_ch3br_pef()
    write(*, "(a,i0,a)") "Using embedded sty1542 CH3Br PEF (", npar_active, " linear terms)"
  end if

  call poten_CH3Br(v_ref, Re_p, re_ch_p, re_ch_p, re_ch_p, betae_p, betae_p, betae_p, 120.0_rk, 120.0_rk)
  write(*, "(a,3f12.6,a,f14.6,a)") "Reference geometry (Re, re, betae) = ", Re_p, re_ch_p, betae_p, " -> V(ref) = ", v_ref, " cm^-1"

  if (nargs >= 2) then
    call get_command_argument(2, arg2)
  else
    arg2 = "../grids/ch3br_pes_1d_2d_aug-cc-pVQZ.dat"
  end if

  open(newunit=iu, file=trim(arg2), status="old", action="read", iostat=ios)
  if (ios == 0) then
    npts = 0
    ssq  = 0.0_rk
    do
      read(iu, "(a)", iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
      read(line, *) c(1:9), e_ha, e_ab
      call poten_CH3Br(v_calc, c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9))
      diff = e_ab - v_calc
      ssq  = ssq + diff * diff
      npts = npts + 1
    end do
    close(iu)
    write(*, "(a,i0,a,a,a,f12.6,a)") "Evaluated ", npts, " geometries from ", trim(arg2), ": unweighted RMS = ", sqrt(ssq / real(npts, rk)), " cm^-1"
  end if
end program test_ch3br_pef

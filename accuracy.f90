module accuracy
  implicit none
  private
  public sik, ik, hik, rk, ark, out, inp, safe_max,safe_min,max_exp, pi, twopi, cl, wl
  public accuracyInitialize
  public planck,avogno,vellgt,boltz,bohr
  public epsil,small_,sqrt2,sqrt3,rad,fititermax,small_a
  !
  integer, parameter :: sik         = selected_int_kind(4)       ! Small integers
  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
                                                                 ! C "int" type, or the DX interface won't
                                                                 ! work correctly.
  integer, parameter :: hik         = selected_int_kind(16)      ! "Pointer" integers - sufficient to store
                                                                 ! memory address
  integer, parameter :: drk         = selected_real_kind(12,25)  ! "Double" reals and complex (complexi? :-)
  integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  integer, parameter :: ark         = selected_real_kind(25,32)  ! "Accurate" reals and complex (complexi? :-)
  integer, parameter :: inp         = 5                          ! Input I/O channel
  integer            :: out         = 6                          ! Output I/O channel
  integer, parameter :: nfilelegendre = 101                      ! Dump-outout channel for eigenfunction 

  ! universal constants
  real(drk), parameter :: planck     =  6.62606896e-27             ! Planck constant
  real(drk), parameter :: avogno     =  6.0221415E+23             ! Avogadro constant
  real(drk), parameter :: vellgt     =  2.99792458E+10            ! Speed of light constant
  real(drk), parameter :: boltz      =  1.380658E-16              ! Boltzmann constant
  real(drk), parameter :: bohr       =  0.529177249               ! a.u.
  
  real(rk)           :: safe_max                                 ! Largest number we want to work with
  real(rk)           :: safe_min                                 ! Smalles number we want to work with
  real(rk)           :: max_exp                                  ! Largest number OK for exponentiating
  real(ark)          :: small_a                                  ! a positive model number that is almost 
  real(rk)           :: small_                                   ! a positive model number that is almost 
                                                                 ! negligible compared to unity in the current model 
  real(ark)          :: pi,twopi                                 ! Pi, 2*Pi
  real(ark)          :: epsil(3,3,3)                             ! epsil - antisymmetric tensor
  real(rk),parameter :: sqrt2 = 1.414213562373095048801689_rk    ! \sqrt{2}
  real(rk),parameter :: sqrt3 = 1.732050807568877293527446_rk    ! \sqrt{3}
  real(rk),parameter :: rad   = 57.295779513082320875_rk    ! radian = 180/pi
  integer, parameter :: cl          = 80                         ! Max character string length
  integer, parameter :: wl          = 500                       ! Very large max character string length 
  integer, parameter :: fititermax  = 200                        ! Max number of iteratons in different fittings 

  contains

    subroutine accuracyInitialize

      safe_max = huge(1.0_rk)**(0.25_rk)
      safe_min = tiny(1.0_rk)**(0.25_rk)
      max_exp  = log(safe_max)
      pi       = 4.0_ark * atan2(1.0_ark,1.0_ark)
      twopi    = 2.0_ark * pi
      small_ = epsilon(1.0_rk )
      small_a= epsilon(1.0_ark)

      epsil = 0 
      epsil(1,2,3) = 1.0_ark
      epsil(1,3,2) =-1.0_ark
      epsil(2,1,3) =-1.0_ark
      epsil(2,3,1) = 1.0_ark
      epsil(3,1,2) = 1.0_ark
      epsil(3,2,1) =-1.0_ark


    end subroutine accuracyInitialize

end module accuracy


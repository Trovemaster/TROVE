module timer
!
!  Timing routines. Because of the possibility for parallel
!  execution, potentially on two scales at the same time
!  (OpenMP and MPI), all timings are done using real time.
!
!  The externally visible routines are:
!
!    TimerStart  - Begins a timing region
!    TimerStop   - Ends a timing region
!    TimerReport - Produces timing report
!
!  Internally, timers are implemented using a hash table and
!  should have constant cost regardless of the number of
!  timers. Hash implementation follows Algorthim L (linear
!  insertion with linear search) of chapter 6.4 of Knuth's
!  vol 3. Hash function is taken from Aho, Sethi, and Ullman,
!  pp. 434-438. Because timers are persistent, there is no
!  deletion.
!
  use accuracy
  implicit none
  private
  public TimerStart, TimerStop, TimerReport, TimerProbe, IOStart, IOStop , ArrayStart, ArrayStop, ArrayMinus, MemoryReport
  public memory_limit,maxmemory,memory_now
  public to_upper, to_lower, tokenize
  !
  integer, parameter :: trk        = selected_real_kind(12)
  integer, parameter :: table_size = 10000 ! Max number of entries to track
  integer, parameter :: tarray_size = 10000 ! Max number of entries to track
  integer, parameter :: name_len   =   40 ! Max length of timer name
  !
  real(rk)           :: maxmemory   =   0 ! Maximal memory allocated 
  real(rk)           :: memory_limit =   100000.0 ! Gb
  !
  type tim
    logical                 :: used       ! Slot used?
    logical                 :: active     ! Currently active?
    character(len=name_len) :: name       ! Timer name
    real(trk)               :: calls      ! Number of times the timer was invoked
                                          !
                                          ! All times below are in seconds
                                          !
    real(trk)               :: real_time  ! Total real time on this timer
    real(trk)               :: cpu_time   ! ditto for CPU time
    real(trk)               :: real_kids  ! Real time spent in nested timers
    real(trk)               :: cpu_kids   ! ditto for CPU time
    real(trk)               :: real_start ! For active timers, time of activation
    real(trk)               :: cpu_start  ! ditto for CPI time
    integer(ik)             :: stack_p    ! For active timers, position in the stack
  end type tim
  !
  !  Timers and sundry data
  !
  type(tim), target :: t_table (table_size) ! All of our timers
  integer(ik)       :: t_nested(table_size) ! Stack of currently active timers
  integer(ik)       :: t_appear(table_size) ! Appearance order for the timers
  integer(ik)       :: t_count  = 0         ! Number of defined timers
  integer(ik)       :: t_active = 0         ! Number of currently active timers
  real(trk)         :: prog_start           ! Timebase for the real time


  ! -----  io-units 
  type tio_unit
    logical                 :: used       ! Slot used?
    logical                 :: active     ! Currently active?
    character(len=name_len) :: name       ! Timer name
    integer(ik)             :: stack_p    ! For active units, position in the stack
    integer(ik)             :: slot     ! Number of times the unit was invoked
  end type tio_unit


  ! -----  io-units 
  type tarray_unit
    logical                 :: used       ! Slot used?
    logical                 :: active     ! Currently active?
    character(len=cl)       :: name       ! Timer name
    integer(ik)             :: stack_p    ! For active arrays, position in the stack
    integer(ik)             :: slot       ! Number of times the unit was invoked
    real(trk)             :: size       ! total size
  end type tarray_unit


  type(tarray_unit), target :: array_table (tarray_size) ! All of our array uints
  integer(ik)       :: array_appear(tarray_size) ! Appearance order for the i/o-unit
  integer(ik)       :: array_count =  0         ! Number of defined io units (we will start counting from unit=10)
  integer(ik)       :: array_active = 0         ! Number of currently active units
  real(trk)         :: memory_now   = 0         ! Number of currently active units


  type(tio_unit), target :: io_table (table_size) ! All of our i/o uints
  integer(ik)       :: io_appear(table_size) ! Appearance order for the i/o-unit
  integer(ik)       :: io_count =  0         ! Number of defined io units (we will start counting from unit=10)
  integer(ik)       :: io_active = 0         ! Number of currently active units


  !
  contains
  !
  !  Public interfaces
  !
    !
    !  Start timing region.
    !
    subroutine TimerStart(region)
      character(len=*), intent(in) :: region  ! Timer name
      !
      integer(ik)          :: pos  ! Timer position
      type(tim), pointer   :: t    ! Current timer (for convenience)
      !
      !  One-time initialization
      !
      if (t_count==0) then
        prog_start      = get_real_time()
        t_table(:)%used = .false.
      end if
      !
      pos =  insert_item(region)
      t   => t_table(pos)
      !
      if (t%active) then
        write (out,"('TimerStart: timer ',a,' is already active')") trim(region)
        stop 'TimerStart - nested timer'
      end if
      !
      !  Push the new timer to the timer stack
      !
      t_active = t_active + 1
      t_nested(t_active) = pos
      !
      t%active     = .true.
      t%stack_p    = t_active
      t%calls      = t%calls + 1
      t%real_start = get_real_time()
      t%cpu_start  = get_cpu_time ()
    end subroutine TimerStart
    !
    !  End timing region.
    !
    subroutine TimerStop(region,real_time_out)
      character(len=*), intent(in) :: region  ! Timer name
      !
      integer(ik)        :: pos        ! Timer position
      type(tim), pointer :: t          ! Current timer (for convenience)
      type(tim), pointer :: pt         ! Parent timer
      real(trk)          :: real_time  ! Real time for this invocation
      real(trk)          :: cpu_time   ! ditto CPU
      real(rk),optional,intent(out)  :: real_time_out  ! Real time for this invocation
      !
      pos =  insert_item(region)
      t   => t_table(pos)
      !
      if (.not.t%active) then
        write (out,"('TimerStop: timer ',a,' is not running')") trim(region)
        stop 'TimerStop - inactive timer'
      end if
      !
      !  Get timings for this invocation
      !
      real_time = get_real_time() - t%real_start
      cpu_time  = get_cpu_time()  - t%cpu_start
      !
      !  Update time counts for the leaving timer
      !
      t%real_time = t%real_time + real_time
      t%cpu_time  = t%cpu_time + cpu_time
      !
      !  Mark as inactive
      !
      t%active    = .false.
      !
      !  Pop the timer from stack, and update counts for the parent
      !
      t_active = t_active - 1
      if (t_active>0) then
        pt => t_table(t_nested(t_active))
        pt%real_kids = pt%real_kids + real_time
        pt%cpu_kids  = pt%cpu_kids  + cpu_time
      end if
      !
      if (present(real_time_out)) real_time_out = real_time
      !
    end subroutine TimerStop


    !
    !  Check time for a region.
    !
    subroutine TimerProbe(region,real_time_out,cpu_time_out)
      character(len=*), intent(in) :: region  ! Timer name
      !
      integer(ik)        :: pos        ! Timer position
      type(tim), pointer :: t          ! Current timer (for convenience)
      type(tim), pointer :: pt         ! Parent timer
      real(rk),intent(out) :: real_time_out  ! Real time for this invocation
      real(rk),intent(out) :: cpu_time_out   ! ditto CPU
      !
      pos =  insert_item(region)
      t   => t_table(pos)
      !
      if (.not.t%active) then
        write (out,"('TimerStop: timer ',a,' is not running')") trim(region)
        !stop 'TimerStop - inactive timer'
      end if
      !
      !  Get timings for this invocation
      !
      real_time_out = get_real_time() - t%real_start
      cpu_time_out  = get_cpu_time()  - t%cpu_start
      !
    end subroutine TimerProbe
    !

    !
    !  Produce timing report
    !
    subroutine TimerReport
      real(trk)          :: real_now
      real(trk)          :: cpu_now
      real(trk)          :: real_time, cpu_time
      real(trk)          :: real_kids, cpu_kids
      real(trk)          :: real_threshold
      real(trk)          :: cpu_threshold
      character(len=19)  :: timestamp
      integer(ik)        :: ord
      integer(ik)        :: pos, kid_pos
      type(tim), pointer :: t, k
      character(len=1)   :: active
      integer(ik)        :: omitted
      !
      real_now  = get_real_time()
      cpu_now   = get_cpu_time()
      timestamp = get_timestamp()
      !
      real_threshold = 0.01_trk * (real_now - prog_start)
      cpu_threshold  = 0.01_trk * cpu_now
      !
      write (out,"(/t15,'Timing data at ',a/)") timestamp
      write (out,"(t2,'     ',t38,'     ',t45,'Total time (seconds)',t67,'Self time (seconds)')")
      write (out,"(t2,'Timer',t38,'Calls',t45,'--------------------',t67,'-------------------')")
      write (out,"(t2,'-----',t38,'-----',t50,'Real',t61,'CPU',t72,'Real',t83,'CPU')")
      write (out,"()")
      !
      omitted = 0
      scan: do ord=1,t_count
        pos = t_appear(ord)
        t => t_table(pos)
        if (.not.t%used) then
          write (out,"('Timer ',i4,' in slot ',i5,' is defined but unused?!')") ord, pos
          stop 'TimerReport - logic error'
        end if
        !
        ! Calculate active-timer corrections
        !
        real_time = 0 ; real_kids = 0 ;
        cpu_time  = 0 ; cpu_kids  = 0 ;
        active     = ' '
        if (t%active) then
          real_time = real_now - t%real_start
          cpu_time  = cpu_now  - t%cpu_start
          if (t_active/=t%stack_p) then
            !
            ! If we are not at the top of the stack, adjust
            ! cumulative children time.
            !
            kid_pos   = t_nested(t%stack_p+1)
            k         => t_table(kid_pos)
            real_kids = real_now - k%real_start
            cpu_kids  = cpu_now  - k%cpu_start
          end if
          active     = '*'
        end if
        !
        real_time = real_time + t%real_time
        cpu_time  = cpu_time  + t%cpu_time
        real_kids = real_kids + t%real_kids
        cpu_kids  = cpu_kids  + t%cpu_kids
        !
        !  Output needed?
        !
        if (real_time<real_threshold .and. cpu_time<cpu_threshold) then
          omitted = omitted + 1
          cycle scan
        end if
        !
        !  Output
        !
        write (out,"(t2,a30,t33,a1,t35,f8.0,t45,2(f9.1,1x,f9.1,3x))") &
               t%name, active, t%calls, real_time, cpu_time, &
               real_time - real_kids, cpu_time - cpu_kids
      end do scan
      if (omitted>0) then
        write (out,"(/' (',i3,' timers contributing less than 1% are not shown)')") &
               omitted
      end if
      write (out,"()")
    end subroutine TimerReport
  !
  !  Support routines
  !
    !
    !  Get real time in seconds. SYSTEM_CLOCK is allowed to roll over,
    !  so this code gets a little messy to handle this.
    !
    function get_real_time() result(t)
      real(trk) :: t
      !
      integer         :: count, count_rate, count_max
      real(trk), save :: overflow   =  0
      integer, save   :: last_count = -1
      !
      call system_clock(count,count_rate,count_max)
      !
      ! Try to detect a rollover
      !
      if (count<last_count) then
        overflow = overflow + count_max
      end if
      last_count = count
      !
      ! Convert to seconds
      !
      t = (overflow+count)/count_rate
    end function get_real_time
    !
    !  Get CPU time, whatever this means
    !
    function get_cpu_time() result(t)
      real(trk) :: t
      !
      call cpu_time(t)
    end function get_cpu_time
    !
    !  Return a (reasonably) pretty time stamp
    !
    function get_timestamp() result(s)
      character(len=19) :: s
      !
      character(len=8)  :: date
      character(len=10) :: time
      !
      call date_and_time(date=date,time=time)
      s = 'YYYY/MM/DD HH:MM:SS'
      s( 1: 4) = date(1:4)
      s( 6: 7) = date(5:6)
      s( 9:10) = date(7:8)
      s(12:13) = time(1:2)
      s(15:16) = time(3:4)
      s(18:19) = time(5:6)
    end function get_timestamp
    !
    !  Linear insertion with linear search (Algorithm L)
    !
    function insert_item(name) result(pos)
      character(len=*), intent(in) :: name
      integer(ik)                  :: pos
      !
      pos = string_hash(name)
      search: do
        if (.not.t_table(pos)%used) then
          !
          ! This is a new key, insert it
          !
          t_count = t_count + 1
          if (t_count>=table_size/5) then
            write (out,"('Too many timers. Increase table_size in "// &
                       "timer.f90 to at least ',i5)") t_count*5
            stop 'timer%insert_item'
          end if
          t_appear(t_count)      = pos
          t_table(pos)%used      = .true.
          t_table(pos)%active    = .false.
          t_table(pos)%name      = name
          t_table(pos)%calls     = 0
          t_table(pos)%real_time = 0
          t_table(pos)%cpu_time  = 0
          t_table(pos)%real_kids = 0
          t_table(pos)%cpu_kids  = 0
          exit search
        end if
        if (t_table(pos)%name==name) then
          !
          ! This is an existing key, simply return the location
          !
          exit search
        end if
        pos = 1 + modulo(pos-2,table_size)
      end do search
      !
    end function insert_item
    !
    integer function string_hash(str) result(h)
      character(len=*), intent(in) :: str
!
      integer :: i, chr, g
      integer :: mask
      data mask/Z"1FFFFFF"/
!
!    This hash assumes at least 29-bit integers. It is supposedly
!    documented in Aho, Sethi, and Ullman, pp. 434-438
!
      h = 0
      do i=1,len_trim(str)
        chr = ichar(str(i:i))
        h   = ishft(h,  4) + chr
        g   = ishft(h,-24)
        h   = iand(ieor(h,g),mask)
      end do
      h = 1 + modulo(h,table_size)
    end function string_hash
   !

    !
    !  Start new i/o unit .
    !
    subroutine IOStart(name,slot)
      character(len=*), intent(in) :: name  ! Unit name
      integer(ik), intent(out)     :: slot  ! Unit slot
      !
      integer(ik)          :: pos  ! Unit position
      type(tio_unit), pointer   :: t    ! Current io_unit (for convenience)
      logical                   :: ifopen
      !
      !  One-time initialization
      !
      if (io_count==0) then
        io_count = 10
        io_table(:)%used = .false.
      end if
      !
      pos =  insert_iounit(name)
      t   => io_table(pos)
      !
      if (t%active) then
        slot  = t%slot
        return
      end if
      !
      inquire (t%slot,OPENED=ifopen)
      !
      if (ifopen) then
        write (out,"('IOStart: unit ',a,' is already open')") trim(name)
        stop 'IOStart - unit was open before'
      end if
      !
      !  Push the new timer to the timer stack
      !
      io_active = io_active + 1
      slot  = t%slot
      !
      t%active     = .true.

    end subroutine IOStart


    !
    !  End an i/o-unit
    !
    subroutine IOStop(name)
      character(len=*), intent(in) :: name  ! Unit name
      !
      integer(ik)        :: pos        ! unit position
      type(tio_unit), pointer :: t          ! Current unit (for convenience)
      !
      pos =  insert_iounit(name)
      t   => io_table(pos)
      !
      if (.not.t%active) then
        write (out,"('IOStop: timer ',a,' is not running')") trim(name)
        stop 'IOStop - inactive IO-slot'
      end if
      !
      !  Get timings for this invocation
      !
      !  Mark as inactive
      !
      t%active    = .false.
      !
      !  Pop the timer from stack, and update counts for the parent
      !
      io_active = io_active - 1
      !if (t_active>0) then
      !  pt => io_table(t_nested(t_active))
      !end if
      !
    end subroutine IOStop


    !
    !  Linear insertion with linear search (Algorithm L)
    !
    function insert_iounit(name) result(pos)
      character(len=*), intent(in) :: name
      integer(ik)                  :: pos
      !
      pos = string_hash(name)
      search: do
        if (.not.io_table(pos)%used) then
          !
          ! This is a new key, insert it
          !
          io_count = io_count + 1
          if (io_count>=table_size/5) then
            write (out,"('Too many io_units. Increase table_size in "// &
                       "timer.f90 to at least ',i5)") io_count*5
            stop 'timer%insert_item'
          end if
          io_appear(io_count)      = pos
          io_table(pos)%used      = .true.
          io_table(pos)%active    = .false.
          io_table(pos)%name      = name
          io_table(pos)%slot      = io_count
          exit search
        end if
        if (io_table(pos)%name==name) then
          !
          ! This is an existing key, simply return the location
          !
          exit search
        end if
        pos = 1 + modulo(pos-2,table_size)
      end do search
      !
    end function insert_iounit



    !
    !  Start new array.
    !
    subroutine ArrayStart(name,alloc,isize,ikind,hsize)
      character(len=*), intent(in) ::  name  ! Unit name
      integer(ik), intent(in)      :: alloc  ! allocation error
      integer(ik), intent(in)      :: isize  ! Unit size
      integer(ik), intent(in)      :: ikind  ! Unit kind
      integer(hik), intent(in),optional   :: hsize  ! Unit size
      !
      integer(ik)          :: pos  ! Unit position
      type(tarray_unit), pointer   :: t    ! Current io_unit (for convenience)
      logical                   :: ifopen
      real(rk)                  :: size_
      integer(hik)              :: hsize_  
      !
      !  One-time initialization
      !
      hsize_ = int(isize,hik)
      if (present(hsize)) hsize_ = hsize
      !
      size_ = (ikind*real(hsize_))/real(1024**3) ! size in GByte

      if (alloc/=0) then
          write(out,"(/' Error ',i8,' trying to allocate array ',a)") alloc,name
          write(out,"( ' Array dimension = ',i0,' array size =  ',f20.2,' Gb')") hsize_,size_
          call MemoryReport
          stop 'ArrayStart - allocation error'
      end if
      !
      if (array_count==0) then
        array_table(:)%used = .false.
      end if
      !
      pos =  insert_arrayunit(name)
      t   => array_table(pos)
      !
      if (size_<0.0_trk) then
          write (out,"(' Size of array ',a,' is negative ',f20.8)") name,size_
          size_ = 1
          call MemoryReport
          !
          !stop 'ArrayStart - negative size'
      end if
      !
      if (.not.t%active) then
        !write (out,"('ArrayStart: ArrayStart ',a,' is already active')") trim(name)
        !stop 'ArrayStart - nested ArrayStart'
        !
        !  Push the new array to the array stack
        !
        array_active = array_active + 1
        t%active     = .true.
        !
      end if
      !
      memory_now   = memory_now + size_
      !
      maxmemory = max(memory_now,maxmemory)
      !
      if (memory_now>memory_limit) then 
        !
        write(out,"('Run out of memory: check your input memory')")
        call MemoryReport
        stop 'Run out of memory'
        !
      endif 
      !
      if (trim(array_table(pos)%name)==name) then
          !
          t%size = t%size + size_ 
          !
      end if
      !
      !slot  = t%slot
      !
    end subroutine ArrayStart


    !
    !  End an array-unit
    !
    subroutine ArrayStop(name)
      character(len=*), intent(in) :: name  ! Unit name
      !
      integer(ik)        :: pos        ! unit position
      type(tarray_unit), pointer :: t          ! Current unit (for convenience)
      real(rk)           :: mem
      !
      pos =  insert_arrayunit(name)
      t   => array_table(pos)
      !
      if (.not.t%active) then
        write (out,"('ArrayStop: array ',a,' is not running')") trim(name)
        stop 'ArrayStop - inactive array counter'
      end if
      !
      memory_now   = memory_now - t%size
      t%size = 0 
      mem = memory_now
      !
      !if (memory_now<=0) then 
        !
        !  Mark as inactive
        !
        t%active    = .false.
        ! 
        !  Pop the timer from stack, and update counts for the parent
        !
        array_active = array_active - 1
        !
      !endif 
      !
    end subroutine ArrayStop


    !
    !  allocation
    !
    subroutine ArrayMinus(name,isize,ikind,hsize)
      character(len=*), intent(in) :: name  ! Unit name
      integer(ik), intent(in)      :: isize  ! Unit size
      integer(ik), intent(in)      :: ikind  ! Unit kind
      integer(hik), intent(in),optional   :: hsize  ! Unit size
      !
      real(hik)                    :: hsize_
      !
      integer(ik)        :: pos        ! unit position
      type(tarray_unit), pointer :: t          ! Current unit (for convenience)
      real(rk)           :: mem,size_
      !
      pos =  insert_arrayunit(name)
      t   => array_table(pos)
      !
      if (.not.t%active) then
        write (out,"('ArrayStop: array ',a,' is not running')") trim(name)
        stop 'ArrayStop - inactive array counter'
      end if
      !
      hsize_ = int(isize,hik)
      if (present(hsize)) hsize_ = hsize
      !
      size_ = (ikind*real(hsize_))/real(1024**3) ! size in GByte
      !
      memory_now   = memory_now - size_
      t%size = max(t%size-size_,0.0_rk)
      mem = memory_now
      !
      !if (t%size<=0) then 
      !  !
      !  !  Mark as inactive
      !  !
      !  t%active    = .false.
      !  ! 
      !  !  Pop the counter from stack, and update counts for the parent
      !  !
      !  array_active = array_active - 1
      !  !
      !endif 
      !
    end subroutine ArrayMinus


    !
    !  Linear insertion with linear search (Algorithm L)
    !
    function insert_arrayunit(name) result(pos)
      character(len=*), intent(in)  :: name
      integer(ik)                  :: pos,ord
      type(tarray_unit), pointer :: t
      !
      pos = string_hash(name)
      search: do
        if (.not.array_table(pos)%used) then
          !
          ! This is a new key, insert it
          !
          array_count = array_count + 1
          if (array_count>=tarray_size/4) then
            write (out,"('Too many array_units. Increase tarray_size in "// &
                       "timer.f90 to at least ',i5)") array_count*4+1

            scan: do ord=1,array_count
              pos = array_appear(ord)
              t => array_table(pos)
              !
              ! Calculate active-array corrections
              !
              !  Output
              !
              write (out,"(t2,a,t50,t55,e11.4)") &
                     t%name, t%size
            end do scan


            stop 'timer%insert_item'
          end if
          array_appear(array_count)      = pos
          array_table(pos)%used      = .true.
          array_table(pos)%active    = .false.
          array_table(pos)%name      = name
          array_table(pos)%slot      = array_count
          array_table(pos)%size      = 0
          exit search
        endif 
        if (trim(array_table(pos)%name)==name) then
          !
          ! This is an existing key, simply return the location
          !
          !array_table(pos)%size = array_table(pos)%size + size_ 
          !
          exit search
          !
        end if
        pos = 1 + modulo(pos-2,tarray_size)
      end do search
      !
    end function insert_arrayunit


    !
    !  Produce timing report
    !
    subroutine MemoryReport
      real(trk)          :: mem
      real(trk)          :: mem_threshold
      integer(ik)        :: ord
      integer(ik)        :: pos, kid_pos
      type(tarray_unit), pointer :: t, k
      character(len=1)   :: active
      integer(ik)        :: omitted
      !
      ! memory_now  = array_size
      !
      maxmemory = max(memory_now,maxmemory)
      !
      mem_threshold = 0.01_trk * memory_now
      !
      write (out,"(//t2,'Memory Report:')")
      write (out,"(t2,'Active Arrays',t55,'size (Gb)')")
      !write (out,"()")
      !
      omitted = 0
      scan: do ord=1,array_count
        pos = array_appear(ord)
        t => array_table(pos)
        if (.not.t%used) then
          write (out,"('Array ',i4,' in slot ',i5,' is defined but unused?!')") ord, pos
          stop 'MemoryReport - logic error'
        end if
        !
        ! Calculate active-array corrections
        !
        mem= 0 
        !active     = ' '
        !if (t%active) then
        !  mem = memory_now
        !  active     = '*'
        !end if
        !
        mem = mem   + t%size
        !
        !  Output needed?
        !
        if (mem<mem_threshold.or.mem<1e-7.or..not.t%active) then
          omitted = omitted + 1
          cycle scan
        end if
        !
        !  Output
        !
        write (out,"(t2,a,t50,t55,e11.4)") &
               t%name, t%size
      end do scan
      write (out,"(t2,'Total memory   = ',t47,f18.8,' Gb')") memory_now
      write (out,"(t2,'Maximal memory = ',t47,f18.8,' Gb (',f16.1,')')") maxmemory,memory_limit


      if (omitted>0) then
        write (out,"(/' (',i9,' arrays contributing less than 1% are not shown)')") &
               omitted
      end if
      write (out,"()")
    end subroutine MemoryReport


function to_upper(str) result (string)

  character(*), intent(in) :: str
  character(len(str)) :: string

  integer :: ic, i

  character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

  string = str
  do i = 1, len_trim(str)
    ic = index(low, str(i:i))
    if (ic>0) string(i:i) = cap(ic:ic)
  end do

end function to_upper


function to_lower(str) result (string)

  character(*), intent(in) :: str
  character(len(str)) :: string

  integer :: ic, i

  character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

  string = str
  do i = 1, len_trim(str)
    ic = index(cap, str(i:i))
    if (ic>0) string(i:i) = low(ic:ic)
  end do

end function to_lower


subroutine tokenize(str, n, w)

  character(*), intent(in) :: str
  integer(ik), intent(out) :: n
  character(cl), intent(out) :: w(:)

  integer(ik) :: pos2, pos1, maxn, ln
  character(1000) :: s

  maxn = size(w)
  n = 0
  pos1 = 1

  DO
    s = trim(adjustl(str))
    ln = len(trim(adjustl(str)))
    if (pos1>ln) exit
    pos2 = INDEX(s(pos1:ln), ' ')
    IF (pos2 == 0) THEN
       n = n + 1
       if (n>maxn) then
         write(out, '(/a,1x,i4)') 'timer/tokenize error: number of words exceeds maximum =', maxn
         stop 'STOP'
       endif
       w(n) = s(pos1:ln)
       EXIT
    END IF
    if (s(pos1:pos1+pos2-1)=='') then
      pos1 = pos2+pos1
      cycle
    endif
    n = n + 1
    if (n>maxn) then
      write(out, '(/a,1x,i4)') 'timer/tokenize error: number of words exceeds maximum =', maxn
      stop 'STOP'
    endif
    w(n) = s(pos1:pos1+pos2-1)
    pos1 = pos2+pos1
 END DO

end subroutine tokenize


end module timer

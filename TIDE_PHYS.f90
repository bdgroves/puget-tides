!=====================================================================
! TIDE_PHYS.f90   SHARED BY PUGET-TIDES AND NORMALS
!   tide_util  CSV SPLITTING, DATES, PACIFIC TIME, NUMBER FORMATTING
!   tide_phys  THE TIDE FROM ITS HARMONIC CONSTITUENTS, HIGHS AND
!              LOWS, THE SUN, PERCENTILE CLASSES
!=====================================================================
module tide_util
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: miss = -9999.0_dp
  integer, parameter :: maxf = 24
  character(len=3), parameter :: mon(12) = [character(len=3) :: &
    'JAN','FEB','MAR','APR','MAY','JUN', &
    'JUL','AUG','SEP','OCT','NOV','DEC']
contains

  logical function ok(x)
    real(dp), intent(in) :: x
    ok = x > miss + 1.0_dp
  end function ok

  ! split a comma-separated line; empty fields stay empty
  subroutine split(line, f, n)
    character(len=*), intent(in) :: line
    character(len=64), intent(out) :: f(maxf)
    integer, intent(out) :: n
    integer :: i, s, l
    f = ' '
    n = 1
    s = 1
    l = len_trim(line)
    do i = 1, l
      if (line(i:i) == ',') then
        if (i > s) f(n) = adjustl(line(s:i-1))
        if (n < maxf) n = n + 1
        s = i + 1
      end if
    end do
    if (l >= s) f(n) = adjustl(line(s:l))
  end subroutine split

  real(dp) function rval(s)
    character(len=*), intent(in) :: s
    integer :: ios
    rval = miss
    if (len_trim(s) == 0) return
    read(s, *, iostat=ios) rval
    if (ios /= 0) rval = miss
  end function rval

  ! julian day number of a calendar date
  integer function jdn(y, m, d)
    integer, intent(in) :: y, m, d
    integer :: a, yy, mm
    a = (14 - m) / 12
    yy = y + 4800 - a
    mm = m + 12*a - 3
    jdn = d + (153*mm + 2)/5 + 365*yy + yy/4 - yy/100 &
        + yy/400 - 32045
  end function jdn

  subroutine ymd(j, y, m, d)
    integer, intent(in) :: j
    integer, intent(out) :: y, m, d
    integer :: a, b, c, e, f
    a = j + 32044
    b = (4*a + 3) / 146097
    c = a - 146097*b/4
    e = (4*c + 3) / 1461
    f = c - 1461*e/4
    m = (5*f + 2) / 153
    d = f - (153*m + 2)/5 + 1
    y = 100*b + e - 4800 + m/10
    m = m + 3 - 12*(m/10)
  end subroutine ymd

  integer function jdn_s(s)          ! 'YYYY-MM-DD'
    character(len=*), intent(in) :: s
    integer :: y, m, d, ios
    jdn_s = 0
    if (len_trim(s) < 10) return
    read(s(1:4), *, iostat=ios) y
    if (ios /= 0) return
    read(s(6:7), *, iostat=ios) m
    if (ios /= 0) return
    read(s(9:10), *, iostat=ios) d
    if (ios /= 0) return
    jdn_s = jdn(y, m, d)
  end function jdn_s

  character(len=10) function dstr(j)
    integer, intent(in) :: j
    integer :: y, m, d
    call ymd(j, y, m, d)
    write(dstr, '(I4.4,A,I2.2,A,I2.2)') y, '-', m, '-', d
  end function dstr

  character(len=6) function mstr(j)  ! 'APR 02'
    integer, intent(in) :: j
    integer :: y, m, d
    call ymd(j, y, m, d)
    write(mstr, '(A3,1X,I2.2)') mon(m), d
  end function mstr

  ! fixed-width number for the report; '--' when missing
  function rf(x, w, d) result(s)
    real(dp), intent(in) :: x
    integer, intent(in) :: w, d
    character(len=w) :: s
    character(len=16) :: fm
    if (.not. ok(x)) then
      s = repeat(' ', w-2) // '--'
      return
    end if
    if (d == 0) then
      write(fm, '(A,I0,A)') '(I', w, ')'
      write(s, fm) nint(x)
    else
      write(fm, '(A,I0,A,I0,A)') '(F', w, '.', d, ')'
      write(s, fm) x
    end if
  end function rf

  ! trimmed number for csv; empty when missing
  function cs(x, d) result(s)
    real(dp), intent(in) :: x
    integer, intent(in) :: d
    character(len=:), allocatable :: s
    character(len=24) :: b
    character(len=16) :: fm
    if (.not. ok(x)) then
      s = ''
      return
    end if
    if (d == 0) then
      write(b, '(I24)') nint(x)
    else
      write(fm, '(A,I0,A)') '(F24.', d, ')'
      if (abs(x) < 0.5_dp * 10.0_dp**(-d)) then
        write(b, fm) 0.0_dp                ! never '-0.00'
      else
        write(b, fm) x
      end if
    end if
    s = trim(adjustl(b))
    if (s(1:1) == '.') s = '0' // s
    if (len(s) > 1) then
      if (s(1:2) == '-.') s = '-0' // s(2:)
    end if
  end function cs

  character(len=20) function itoa(i)
    integer, intent(in) :: i
    write(itoa, '(I0)') i
  end function itoa

  ! day of year 1-366 on a leap-year calendar, so Feb 29 has its
  ! own slot and Mar 1 is always day 61
  integer function doy366(m, d)
    integer, intent(in) :: m, d
    integer, parameter :: c(12) = [0, 31, 60, 91, 121, 152, &
                                   182, 213, 244, 274, 305, 335]
    doy366 = c(m) + d
  end function doy366

  ! heapsort, ascending, in place
  subroutine hsort(n, a)
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(n)
    integer :: i, k
    real(dp) :: t
    do i = n / 2, 1, -1
      call sift(i, n)
    end do
    do k = n, 2, -1
      t = a(1)
      a(1) = a(k)
      a(k) = t
      call sift(1, k - 1)
    end do
  contains
    subroutine sift(lo, hi)
      integer, intent(in) :: lo, hi
      integer :: r, c
      real(dp) :: v
      r = lo
      v = a(r)
      do
        c = 2 * r
        if (c > hi) exit
        if (c < hi) then
          if (a(c+1) > a(c)) c = c + 1
        end if
        if (a(c) <= v) exit
        a(r) = a(c)
        r = c
      end do
      a(r) = v
    end subroutine sift
  end subroutine hsort

  ! Weibull plotting position: the value at P*(N+1) in a sorted
  ! sample, interpolated; clamped to the ends
  real(dp) function pctl(n, a, p)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), p
    real(dp) :: r
    integer :: k
    pctl = miss
    if (n < 1) return
    r = p * (n + 1)
    k = int(r)
    if (k < 1) then
      pctl = a(1)
    else if (k >= n) then
      pctl = a(n)
    else
      pctl = a(k) + (r - k) * (a(k+1) - a(k))
    end if
  end function pctl


  ! a time that couldn't be read. Not miss: hours before 2000 are
  ! negative, and -9999 is a real hour in November 1998
  logical function tok(t)
    real(dp), intent(in) :: t
    tok = t > -1.0e29_dp
  end function tok

  ! hours from J2000 (2000-01-01 12:00 UTC) for 'YYYY-MM-DD HH:MM'
  real(dp) function hrs(s)
    character(len=*), intent(in) :: s
    integer :: hh, mi, ios
    hrs = -1.0e30_dp
    if (jdn_s(s) == 0 .or. len_trim(s) < 16) return
    read(s(12:13), *, iostat=ios) hh
    if (ios /= 0) return
    read(s(15:16), *, iostat=ios) mi
    if (ios /= 0) return
    hrs = (jdn_s(s) - 2451545) * 24.0_dp - 12.0_dp + hh + mi / 60.0_dp
  end function hrs

  ! julian day number of the UTC day an hour falls in
  integer function jday(t)
    real(dp), intent(in) :: t
    jday = 2451545 + floor((t + 12.0_dp) / 24.0_dp + 1.0e-9_dp)
  end function jday

  ! 'YYYY-MM-DDTHH:MMZ', rounded to the minute
  character(len=17) function tstr(t)
    real(dp), intent(in) :: t
    real(dp) :: tm
    integer :: m
    tm = nint(t * 60.0_dp) / 60.0_dp
    m = nint((tm + 12.0_dp - 24.0_dp * floor((tm + 12.0_dp) &
        / 24.0_dp + 1.0e-9_dp)) * 60.0_dp)
    write(tstr, '(A10,A,I2.2,A,I2.2,A)') dstr(jday(tm)), 'T', &
      m / 60, ':', mod(m, 60), 'Z'
  end function tstr

  ! day of the week, 0 = Sunday
  integer function dow(j)
    integer, intent(in) :: j
    dow = mod(j + 1, 7)
  end function dow

  ! Pacific time's offset from UTC, hours: daylight time from 2 a.m.
  ! on the second Sunday of March to 2 a.m. on the first Sunday of
  ! November
  real(dp) function pac_off(t)
    real(dp), intent(in) :: t
    integer :: y, m, d, j1, j2
    real(dp) :: a, b
    call ymd(jday(t), y, m, d)
    j1 = jdn(y, 3, 8)
    j1 = j1 + mod(7 - dow(j1), 7)
    j2 = jdn(y, 11, 1)
    j2 = j2 + mod(7 - dow(j2), 7)
    a = (j1 - 2451545) * 24.0_dp - 12.0_dp + 10.0_dp
    b = (j2 - 2451545) * 24.0_dp - 12.0_dp + 9.0_dp
    pac_off = -8.0_dp
    if (t >= a .and. t < b) pac_off = -7.0_dp
  end function pac_off

  ! Pacific local time, 'YYYY-MM-DD HH:MM', and its zone
  character(len=20) function lstr(t)
    real(dp), intent(in) :: t
    character(len=17) :: u
    u = tstr(t + pac_off(t))
    lstr = u(1:10) // ' ' // u(12:16)
    if (pac_off(t) > -7.5_dp) then
      lstr = trim(lstr) // ' PDT'
    else
      lstr = trim(lstr) // ' PST'
    end if
  end function lstr

end module tide_util

!---------------------------------------------------------------------
! tide_phys: THE TIDE FROM ITS HARMONIC CONSTITUENTS
!   h(t) = Z0 + SUM f H cos(V + u - kappa), after Schureman (1958),
!   Manual of Harmonic Analysis and Prediction of Tides, the method
!   NOAA's own tide-predicting machines and programs used.
!   t is in hours from 2000-01-01 12:00 UTC (J2000).
!---------------------------------------------------------------------
module tide_phys
  use tide_util
  implicit none
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: d2r = pi / 180.0_dp, r2d = 180.0_dp / pi
  integer, parameter :: nc = 37
  ! NOAA's order, so harcon.json maps straight across
  character(len=4), parameter :: cname(nc) = [character(len=4) :: &
    'M2','S2','N2','K1','M4','O1','M6','MK3','S4','MN4', &
    'NU2','S6','MU2','2N2','OO1','LAM2','S1','M1','J1','MM', &
    'SSA','SA','MSF','MF','RHO','Q1','T2','R2','2Q1','P1', &
    '2SM2','M3','L2','2MK3','K2','M8','MS4']
  ! the constituents of one station
  type station_tide
    real(dp) :: z0 = 0.0_dp            ! MSL above the chart datum
    real(dp) :: amp(nc) = 0.0_dp       ! H, feet
    real(dp) :: kap(nc) = 0.0_dp       ! kappa, Greenwich phase, deg
  end type station_tide
  ! node factors and corrections, kept for one year at a time
  integer, private :: fu_year = -1
  real(dp), private :: fcache(nc), ucache(nc)
contains

  ! hours from J2000 for a calendar date and UTC hour
  real(dp) function hours_j2000(y, m, d, hr)
    integer, intent(in) :: y, m, d
    real(dp), intent(in) :: hr
    hours_j2000 = (jdn(y, m, d) - 2451545) * 24.0_dp - 12.0_dp + hr
  end function hours_j2000

  ! the calendar year an hour falls in
  integer function year_of(t)
    real(dp), intent(in) :: t
    integer :: y, m, d
    call ymd(2451545 + floor((t + 12.0_dp) / 24.0_dp), y, m, d)
    year_of = y
  end function year_of

  ! mean longitudes, degrees: moon s, sun h, lunar perigee p,
  ! moon's node N, solar perigee p1 (Meeus 1998, linear terms)
  subroutine astro(t, s, h, p, en, p1)
    real(dp), intent(in) :: t
    real(dp), intent(out) :: s, h, p, en, p1
    real(dp) :: c
    c = t / 876600.0_dp                       ! julian centuries
    s  = modulo(218.3164477_dp + 481267.88123421_dp * c, 360.0_dp)
    h  = modulo(280.46646_dp + 36000.76983_dp * c, 360.0_dp)
    p  = modulo(83.3532465_dp + 4069.0137287_dp * c, 360.0_dp)
    en = modulo(125.04452_dp - 1934.136261_dp * c, 360.0_dp)
    p1 = modulo(282.93768_dp + 1.71946_dp * c, 360.0_dp)
  end subroutine astro

  ! speeds of T s h p N p1, degrees per hour
  subroutine rates(r)
    real(dp), intent(out) :: r(6)
    r = [15.0_dp, 481267.88123421_dp, 36000.76983_dp, &
         4069.0137287_dp, -1934.136261_dp, 1.71946_dp]
    r(2:6) = r(2:6) / 876600.0_dp
  end subroutine rates

  ! equilibrium argument V, node correction u and node factor f
  ! of every constituent at hour t; speeds in deg/hour
  subroutine args(t, v, u, f, spd)
    real(dp), intent(in) :: t
    real(dp), intent(out) :: v(nc), u(nc), f(nc), spd(nc)
    real(dp) :: s, h, p, en, p1, tt, r(6), w
    real(dp) :: ci, ai, xi, nu, nup, nupp2, pp, ar, ra, aq, qa
    real(dp) :: fm2, fo1, fk1, fk2, fj1, foo1, fmm, fmf
    real(dp) :: um2, uo1, uk1, uk2
    real(dp), parameter :: om = 23.452_dp * d2r   ! obliquity
    real(dp), parameter :: ei = 5.145_dp * d2r    ! moon's orbit
    real(dp) :: a1, a2, n
    call astro(t, s, h, p, en, p1)
    call rates(r)
    tt = modulo(15.0_dp * t, 360.0_dp)  ! hour angle of mean sun
    n = en * d2r
    ! inclination of the moon's orbit to the equator, I
    ci = cos(om) * cos(ei) - sin(om) * sin(ei) * cos(n)
    ai = acos(ci)
    ! nu and xi from the half-angle formulas
    a1 = atan2(cos(0.5_dp*(om-ei)) * sin(0.5_dp*n), &
               cos(0.5_dp*(om+ei)) * cos(0.5_dp*n))
    a2 = atan2(sin(0.5_dp*(om-ei)) * sin(0.5_dp*n), &
               sin(0.5_dp*(om+ei)) * cos(0.5_dp*n))
    nu = a1 - a2
    xi = n - a1 - a2
    xi = atan2(sin(xi), cos(xi))
    nup = atan2(sin(2*ai) * sin(nu), sin(2*ai) * cos(nu) + 0.3347_dp)
    nupp2 = atan2(sin(ai)**2 * sin(2*nu), &
                  sin(ai)**2 * cos(2*nu) + 0.0727_dp)
    pp = p * d2r - xi                     ! P = p - xi
    ! L2: R and 1/Ra
    ar = atan2(sin(2*pp), 1.0_dp/(6.0_dp*tan(0.5_dp*ai)**2) &
               - cos(2*pp))
    ra = sqrt(1.0_dp - 12.0_dp*tan(0.5_dp*ai)**2*cos(2*pp) &
              + 36.0_dp*tan(0.5_dp*ai)**4)
    ! M1: Q and 1/Qa
    aq = atan2((5*cos(ai) - 1) * sin(pp), (7*cos(ai) + 1) * cos(pp))
    qa = sqrt(0.25_dp + 1.5_dp*cos(ai)/cos(0.5_dp*ai)**2*cos(2*pp) &
              + 2.25_dp*cos(ai)**2/cos(0.5_dp*ai)**4)
    ! node factors
    fm2 = cos(0.5_dp*ai)**4 / 0.9154_dp
    fo1 = sin(ai) * cos(0.5_dp*ai)**2 / 0.3800_dp
    fk1 = sqrt(0.8965_dp*sin(2*ai)**2 + &
               0.6001_dp*sin(2*ai)*cos(nu) + 0.1006_dp)
    fk2 = sqrt(19.0444_dp*sin(ai)**4 + &
               2.7702_dp*sin(ai)**2*cos(2*nu) + 0.0981_dp)
    fj1 = sin(2*ai) / 0.7214_dp
    foo1 = sin(ai) * sin(0.5_dp*ai)**2 / 0.0164_dp
    fmm = (2.0_dp/3.0_dp - sin(ai)**2) / 0.5021_dp
    fmf = sin(ai)**2 / 0.1578_dp
    ! node corrections, degrees
    um2 = (2*xi - 2*nu) * r2d
    uo1 = (2*xi - nu) * r2d
    uk1 = -nup * r2d
    uk2 = -nupp2 * r2d
    ! the base constituents: Doodson T s h p N p1, extra degrees
    call base(1,  2,-2, 2, 0, 0, 0,   0.0_dp, um2, fm2)        ! M2
    call base(2,  2, 0, 0, 0, 0, 0,   0.0_dp, 0.0_dp, 1.0_dp)  ! S2
    call base(3,  2,-3, 2, 1, 0, 0,   0.0_dp, um2, fm2)        ! N2
    call base(4,  1, 0, 1, 0, 0, 0, -90.0_dp, uk1, fk1)        ! K1
    call base(6,  1,-2, 1, 0, 0, 0,  90.0_dp, uo1, fo1)        ! O1
    call base(11, 2,-3, 4,-1, 0, 0,   0.0_dp, um2, fm2)        ! NU2
    call base(13, 2,-4, 4, 0, 0, 0,   0.0_dp, um2, fm2)        ! MU2
    call base(14, 2,-4, 2, 2, 0, 0,   0.0_dp, um2, fm2)        ! 2N2
    call base(15, 1, 2, 1, 0, 0, 0, -90.0_dp, &
              (-2*xi - nu)*r2d, foo1)                          ! OO1
    call base(16, 2,-1, 0, 1, 0, 0, 180.0_dp, um2, fm2)        ! LAM2
    call base(17, 1, 0, 0, 0, 0, 0,   0.0_dp, 0.0_dp, 1.0_dp)  ! S1
    call base(18, 1,-1, 1, 1, 0, 0, -90.0_dp, &
              (xi - nu)*r2d + aq*r2d, fo1*qa)                  ! M1
    call base(19, 1, 1, 1,-1, 0, 0, -90.0_dp, -nu*r2d, fj1)    ! J1
    call base(20, 0, 1, 0,-1, 0, 0,   0.0_dp, 0.0_dp, fmm)     ! MM
    call base(21, 0, 0, 2, 0, 0, 0,   0.0_dp, 0.0_dp, 1.0_dp)  ! SSA
    call base(22, 0, 0, 1, 0, 0, 0,   0.0_dp, 0.0_dp, 1.0_dp)  ! SA
    call base(24, 0, 2, 0, 0, 0, 0,   0.0_dp, -2*xi*r2d, fmf)  ! MF
    call base(25, 1,-3, 3,-1, 0, 0,  90.0_dp, uo1, fo1)        ! RHO
    call base(26, 1,-3, 1, 1, 0, 0,  90.0_dp, uo1, fo1)        ! Q1
    call base(27, 2, 0,-1, 0, 0, 1,   0.0_dp, 0.0_dp, 1.0_dp)  ! T2
    call base(28, 2, 0, 1, 0, 0,-1, 180.0_dp, 0.0_dp, 1.0_dp)  ! R2
    call base(29, 1,-4, 1, 2, 0, 0,  90.0_dp, uo1, fo1)        ! 2Q1
    call base(30, 1, 0,-1, 0, 0, 0,  90.0_dp, 0.0_dp, 1.0_dp)  ! P1
    call base(32, 3,-3, 3, 0, 0, 0,   0.0_dp, 1.5_dp*um2, &
              fm2**1.5_dp)                                     ! M3
    call base(33, 2,-1, 2,-1, 0, 0, 180.0_dp, um2 - ar*r2d, &
              fm2*ra)                                          ! L2
    call base(35, 2, 0, 2, 0, 0, 0,   0.0_dp, uk2, fk2)        ! K2
    ! compound tides: sums of the base ones
    call comb(5,  1, 2.0_dp, 1, 0.0_dp)         ! M4   = 2 M2
    call comb(7,  1, 3.0_dp, 1, 0.0_dp)         ! M6   = 3 M2
    call comb(8,  1, 1.0_dp, 4, 1.0_dp)         ! MK3  = M2 + K1
    call comb(9,  2, 2.0_dp, 2, 0.0_dp)         ! S4   = 2 S2
    call comb(10, 1, 1.0_dp, 3, 1.0_dp)         ! MN4  = M2 + N2
    call comb(12, 2, 3.0_dp, 2, 0.0_dp)         ! S6   = 3 S2
    call comb(23, 2, 1.0_dp, 1, -1.0_dp)        ! MSF  = S2 - M2
    call comb(31, 2, 2.0_dp, 1, -1.0_dp)        ! 2SM2 = 2 S2 - M2
    call comb(34, 1, 2.0_dp, 4, -1.0_dp)        ! 2MK3 = 2 M2 - K1
    call comb(36, 1, 4.0_dp, 1, 0.0_dp)         ! M8   = 4 M2
    call comb(37, 1, 1.0_dp, 2, 1.0_dp)         ! MS4  = M2 + S2
    v = modulo(v, 360.0_dp)
    u = modulo(u + 180.0_dp, 360.0_dp) - 180.0_dp
  contains
    subroutine base(k, a, b, c, d, e, g, x, uu, ff)
      integer, intent(in) :: k, a, b, c, d, e, g
      real(dp), intent(in) :: x, uu, ff
      v(k) = a*tt + b*s + c*h + d*p + e*en + g*p1 + x
      spd(k) = a*r(1) + b*r(2) + c*r(3) + d*r(4) + e*r(5) + g*r(6)
      u(k) = uu
      f(k) = ff
    end subroutine base
    subroutine comb(k, i, a, j, b)
      integer, intent(in) :: k, i, j
      real(dp), intent(in) :: a, b
      v(k) = a*v(i) + b*v(j)
      spd(k) = a*spd(i) + b*spd(j)
      u(k) = a*u(i) + b*u(j)
      w = abs(a)
      f(k) = f(i)**w * f(j)**abs(b)
    end subroutine comb
  end subroutine args

  ! node factors and corrections for the middle of t's year, the
  ! convention NOAA's predictions use
  subroutine nodal(t, f, u)
    real(dp), intent(in) :: t
    real(dp), intent(out) :: f(nc), u(nc)
    real(dp) :: v(nc), spd(nc)
    real(dp) :: s, h, p, en, p1
    integer :: y
    y = year_of(t)
    if (y /= fu_year) then
      call args(hours_j2000(y, 7, 2, 0.0_dp), v, ucache, fcache, spd)
      ! NOAA starts M1 each year from T - s + h - 90, without the
      ! perigee p, but runs it at a speed that includes p's motion;
      ! taking p at January 1 out reproduces their tables
      call astro(hours_j2000(y, 1, 1, 0.0_dp), s, h, p, en, p1)
      ucache(18) = ucache(18) - p
      fu_year = y
    end if
    f = fcache
    u = ucache
  end subroutine nodal

  ! the predicted tide, and its rate of change (ft per hour)
  subroutine tide_at(st, t, h, dh)
    type(station_tide), intent(in) :: st
    real(dp), intent(in) :: t
    real(dp), intent(out) :: h
    real(dp), intent(out), optional :: dh
    real(dp) :: v(nc), u(nc), f(nc), spd(nc), fu(nc), uu(nc), a
    integer :: k
    call args(t, v, u, f, spd)
    call nodal(t, fu, uu)
    ! sum the constituents: h = Z0 + SUM f H cos(V + u - kappa)
    h = st%z0
    if (present(dh)) dh = 0.0_dp
    do k = 1, nc
      if (st%amp(k) == 0.0_dp) cycle
      a = (v(k) + uu(k) - st%kap(k)) * d2r
      h = h + fu(k) * st%amp(k) * cos(a)
      if (present(dh)) dh = dh - fu(k) * st%amp(k) * &
                                 spd(k) * d2r * sin(a)
    end do
  end subroutine tide_at

  ! the turning point between t0 and t1, where the rate changes sign:
  ! bisection on the analytic derivative, to about a second
  real(dp) function turn(st, t0, t1)
    type(station_tide), intent(in) :: st
    real(dp), intent(in) :: t0, t1
    real(dp) :: a, b, m, h, da, dm
    integer :: i
    a = t0
    b = t1
    call tide_at(st, a, h, da)
    do i = 1, 24
      m = 0.5_dp * (a + b)
      call tide_at(st, m, h, dm)
      if (sign(1.0_dp, dm) == sign(1.0_dp, da)) then
        a = m
        da = dm
      else
        b = m
      end if
    end do
    turn = 0.5_dp * (a + b)
  end function turn

  ! every high and low between t0 and t1: scan the rate of change
  ! every 6 minutes and refine each turn. ty is 'H' or 'L'
  subroutine events(st, t0, t1, mx, n, te, he, ty)
    type(station_tide), intent(in) :: st
    real(dp), intent(in) :: t0, t1
    integer, intent(in) :: mx
    integer, intent(out) :: n
    real(dp), intent(out) :: te(mx), he(mx)
    character(len=1), intent(out) :: ty(mx)
    real(dp), parameter :: dt = 0.1_dp
    real(dp) :: t, h, d0, d1, tt
    n = 0
    t = t0
    call tide_at(st, t, h, d0)
    do while (t < t1 .and. n < mx)
      call tide_at(st, t + dt, h, d1)
      if (d0 > 0.0_dp .neqv. d1 > 0.0_dp) then
        tt = turn(st, t, t + dt)
        n = n + 1
        te(n) = tt
        call tide_at(st, tt, he(n))
        ty(n) = 'L'
        if (d0 > 0.0_dp) ty(n) = 'H'
      end if
      d0 = d1
      t = t + dt
    end do
    call prune(n, te, he, ty)
  end subroutine events

  ! drop a high and low next to each other that are within 0.1 ft:
  ! the tide only stands for an hour there, and NOAA's tables leave
  ! such pairs out
  subroutine prune(n, te, he, ty)
    integer, intent(inout) :: n
    real(dp), intent(inout) :: te(:), he(:)
    character(len=1), intent(inout) :: ty(:)
    integer :: a, best
    real(dp) :: dmin
    do
      best = 0
      dmin = 0.1_dp
      do a = 1, n - 1
        if (abs(he(a+1) - he(a)) < dmin) then
          dmin = abs(he(a+1) - he(a))
          best = a
        end if
      end do
      if (best == 0) exit
      te(best:n-2) = te(best+2:n)
      he(best:n-2) = he(best+2:n)
      ty(best:n-2) = ty(best+2:n)
      n = n - 2
    end do
  end subroutine prune

  ! the sun's elevation in degrees at hour t (NOAA's solar position
  ! equations, good to about a minute of time for sunrise and set)
  real(dp) function sun_elev(t, lat, lon)
    real(dp), intent(in) :: t, lat, lon
    real(dp) :: c, l0, m, eqc, sl, ob, ra, dec, gm, ha
    c = t / 876600.0_dp
    l0 = modulo(280.46646_dp + 36000.76983_dp * c, 360.0_dp)
    m = (357.52911_dp + 35999.05029_dp * c) * d2r
    eqc = (1.914602_dp - 0.004817_dp * c) * sin(m) &
        + 0.019993_dp * sin(2 * m) + 0.000289_dp * sin(3 * m)
    sl = (l0 + eqc) * d2r
    ob = (23.439291_dp - 0.0130042_dp * c) * d2r
    ra = atan2(cos(ob) * sin(sl), cos(sl))
    dec = asin(sin(ob) * sin(sl))
    ! Greenwich mean sidereal time, then the local hour angle
    gm = modulo(280.46061837_dp + 360.98564736629_dp * t / 24.0_dp, &
                360.0_dp)
    ha = (gm + lon) * d2r - ra
    sun_elev = asin(sin(lat * d2r) * sin(dec) + &
               cos(lat * d2r) * cos(dec) * cos(ha)) * r2d
  end function sun_elev

  ! a value against its percentiles, the way the USGS classes
  ! streamflow
  character(len=10) function versus(x, p10, p25, p75, p90)
    real(dp), intent(in) :: x, p10, p25, p75, p90
    versus = ' '
    if (.not. ok(x) .or. .not. ok(p10)) return
    if (x < p10) then
      versus = 'MUCH BELOW'
    else if (x < p25) then
      versus = 'BELOW'
    else if (x <= p75) then
      versus = 'NORMAL'
    else if (x <= p90) then
      versus = 'ABOVE'
    else
      versus = 'MUCH ABOVE'
    end if
  end function versus

  ! sea-level trend from monthly means: least squares on time with
  ! one mean per calendar month, so the seasons don't bias it. The
  ! standard error is widened for the months' autocorrelation (AR1),
  ! the way NOAA reports its trends. Returns ft per year.
  subroutine trend(n, yr, mo, x, b, se, rho)
    integer, intent(in) :: n, mo(n)
    real(dp), intent(in) :: yr(n), x(n)
    real(dp), intent(out) :: b, se, rho
    real(dp) :: mx(12), my(12), sxy, sxx, r(n), s2, a, c, neff
    integer :: k, cnt(12), i
    b = miss
    se = miss
    rho = miss
    if (n < 120) return
    mx = 0
    my = 0
    cnt = 0
    do i = 1, n
      mx(mo(i)) = mx(mo(i)) + yr(i)
      my(mo(i)) = my(mo(i)) + x(i)
      cnt(mo(i)) = cnt(mo(i)) + 1
    end do
    do k = 1, 12
      if (cnt(k) > 0) then
        mx(k) = mx(k) / cnt(k)
        my(k) = my(k) / cnt(k)
      end if
    end do
    sxy = 0
    sxx = 0
    do i = 1, n
      sxy = sxy + (yr(i) - mx(mo(i))) * (x(i) - my(mo(i)))
      sxx = sxx + (yr(i) - mx(mo(i)))**2
    end do
    b = sxy / sxx
    s2 = 0
    do i = 1, n
      r(i) = x(i) - my(mo(i)) - b * (yr(i) - mx(mo(i)))
      s2 = s2 + r(i)**2
    end do
    ! lag-one autocorrelation of the residuals, months in a row only
    a = 0
    c = 0
    do i = 2, n
      if (abs(yr(i) - yr(i-1) - 1.0_dp/12.0_dp) < 0.01_dp) then
        a = a + r(i) * r(i-1)
        c = c + r(i-1)**2
      end if
    end do
    rho = 0
    if (c > 0) rho = max(0.0_dp, min(0.95_dp, a / c))
    neff = n * (1.0_dp - rho) / (1.0_dp + rho)
    se = sqrt(s2 / (n - 13) / sxx) * sqrt(n / neff)
  end subroutine trend

end module tide_phys

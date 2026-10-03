!=====================================================================
! PUGET-TIDES.f90   TIDES OF PUGET SOUND
!   EIGHT NOAA TIDE STATIONS FROM NEAH BAY TO TACOMA, TWICE A DAY.
!   PREDICTS THE TIDE FROM NOAA'S HARMONIC CONSTITUENTS, CHECKS ITS
!   HIGHS AND LOWS AGAINST NOAA'S OWN TABLES, MEASURES THE SURGE
!   (WATER LEVEL MINUS TIDE), FINDS THE KING TIDES AND THE DAYLIGHT
!   MINUS TIDES OF THE COMING YEAR, AND FITS THE SEA-LEVEL TREND.
!
!   IN:   fetch_status.csv stations.csv harcon.csv datums.csv
!         water_level.csv pressure.csv monthly_mean.csv
!         noaa_hilo.csv noaa_trend.csv normals.csv (if built)
!   OUT:  now.csv series.csv hilo.csv kingtides.csv minustides.csv
!         travel.csv trend.csv annual.csv summary.csv
!         puget-tides-report.txt
!   RC:   0 NORMAL  4 WARNING (A STATION LATE OR MISSING, NO NORMALS,
!         OR FORTRAN AND NOAA DISAGREE)  8 NO CONSTITUENTS  12 NO FILE
!
!   BUILD gfortran -O2 -o puget-tides TIDE_PHYS.f90 PUGET-TIDES.f90
!=====================================================================
program puget_tides
  use tide_util
  use tide_phys
  implicit none
  integer, parameter :: maxs = 12, maxo = 9000, maxm = 2000
  integer, parameter :: maxe = 64, maxy = 1600, maxk = 10
  character(len=*), parameter :: ver = 'PUGET-TIDES V1.0'
  real(dp), parameter :: m2spd = 28.984104_dp
  real(dp), parameter :: ft2mm = 304.8_dp, hpa_ft = 0.0328084_dp

  ! stations
  integer :: ns = 0
  character(len=8) :: sid(maxs)
  character(len=16) :: sname(maxs)
  real(dp) :: slat(maxs), slon(maxs)
  type(station_tide) :: st(maxs)
  integer :: nhc(maxs) = 0
  ! datums, feet above MLLW
  real(dp) :: dmhhw(maxs) = miss, dmhw(maxs) = miss, dmsl(maxs) = miss
  real(dp) :: dmlw(maxs) = miss, dgt(maxs) = miss, dhat(maxs) = miss
  real(dp) :: dlat(maxs) = miss, dmax(maxs) = miss, dmin(maxs) = miss
  character(len=16) :: dmaxw(maxs) = ' ', dminw(maxs) = ' '
  ! water levels, predictions and surge
  integer :: no(maxs) = 0
  real(dp) :: ot(maxo, maxs), oh(maxo, maxs), op(maxo, maxs)
  ! air pressure
  integer :: np(maxs) = 0
  real(dp) :: pt(maxo, maxs), pv(maxo, maxs)
  ! monthly mean sea level
  integer :: nm(maxs) = 0, my(maxm, maxs), mmo(maxm, maxs)
  real(dp) :: mv(maxm, maxs)
  ! NOAA's predicted highs and lows
  integer :: nn(maxs) = 0
  real(dp) :: nt(maxe, maxs), nh(maxe, maxs)
  character(len=1) :: nty(maxe, maxs)
  ! normals by day of the year: surge and daily highest water
  real(dp) :: nrs(5, 366, maxs), nmx(5, 366, maxs)
  logical :: have_nrm(maxs) = .false.
  ! NOAA's published trend
  real(dp) :: ntr(maxs) = miss, nte(maxs) = miss
  ! results
  integer :: ne(maxs) = 0, nyr(maxs) = 0
  real(dp) :: et(maxe, maxs), eh(maxe, maxs)
  character(len=1) :: ety(maxe, maxs)
  real(dp) :: yt(maxy, maxs), yh(maxy, maxs)
  character(len=1) :: yty(maxy, maxs)
  real(dp) :: tl(maxs) = miss, lv(maxs) = miss, lp(maxs) = miss
  real(dp) :: res24(maxs) = miss, ib24(maxs) = miss, hpa(maxs) = miss
  real(dp) :: rmax(maxs) = miss, rmin(maxs) = miss, rmean(maxs) = miss
  real(dp) :: trmax(maxs) = miss, trmin(maxs) = miss
  real(dp) :: ymax(maxs) = miss
  integer :: ydoy(maxs) = 0, tdoy(maxs) = 0
  character(len=10) :: res_vs(maxs) = ' ', ymax_vs(maxs) = ' '
  integer :: vn(maxs) = 0, vtot(maxs) = 0, vext(maxs) = 0
  real(dp) :: vdt(maxs) = miss, vdh(maxs) = miss
  real(dp) :: tb(maxs) = miss, tse(maxs) = miss, trho(maxs) = miss
  integer :: ty0(maxs) = 0, ty1(maxs) = 0, tnm(maxs) = 0
  integer :: kidx(maxk, maxs), nki(maxs) = 0
  integer :: nminus(maxs) = 0, iminus1(maxs) = 0, iminlow(maxs) = 0

  character(len=512) :: line
  character(len=64) :: f(maxf)
  character(len=20) :: runutc = ' ', pdate = ' '
  character(len=200) :: msgs(200)
  integer :: nmsg = 0, rc = 0, ios, nf, i, k, u, t0i
  real(dp) :: tnow, x
  integer :: current = 0

  call read_status()
  write(*,'(A)') 'TIDE000I ' // ver // ' STARTED ' // trim(runutc)
  call read_stations()
  call read_harcon()
  call read_datums()
  call read_levels()
  call read_pressure()
  call read_monthly()
  call read_noaa_hilo()
  call read_noaa_trend()
  call read_normals()
  if (sum(nhc(1:ns)) == 0) then
    write(*,'(A)') 'TIDE008E NO HARMONIC CONSTITUENTS. NOTHING WRITTEN.'
    call exit(8)
  end if

  do i = 1, ns
    if (nhc(i) < 30) then
      call note('TIDE011W ' // trim(sname(i)) // &
                ': NO CONSTITUENTS, SKIPPED', 4)
      cycle
    end if
    call water_now(i)
    call week(i)
    call verify(i)
    call year_ahead(i)
    call sea_level(i)
  end do
  current = count(tl(1:ns) > tnow - 3.0_dp)

  call write_now()
  call write_series()
  call write_hilo()
  call write_year()
  call write_travel()
  call write_trend()
  call write_summary()
  call write_report()
  write(*,'(A,I4)') 'TIDE001I STATIONS ...................', ns
  write(*,'(A,I4)') 'TIDE002I WITH WATER LEVELS ...........', &
    count(no(1:ns) > 0)
  write(*,'(A,I4)') 'TIDE003I CURRENT (UNDER 3 HOURS OLD) .', current
  write(*,'(A,I4)') 'TIDE004I WITH NORMALS ................', &
    count(have_nrm(1:ns))
  write(*,'(A)') 'TIDE006I REPORT WRITTEN: puget-tides-report.txt'
  write(*,'(A)') 'TIDE007I RESULTS WRITTEN: now.csv series.csv ' // &
    'hilo.csv kingtides.csv minustides.csv travel.csv ' // &
    'trend.csv annual.csv summary.csv'
  do k = 1, nmsg
    write(*,'(A)') trim(msgs(k))
  end do
  write(*,'(A,I2)') 'TIDE999I ' // ver // ' ENDED  RC=', rc
  call exit(rc)

contains

  subroutine note(m, level)
    character(len=*), intent(in) :: m
    integer, intent(in) :: level
    if (nmsg < size(msgs)) then
      nmsg = nmsg + 1
      msgs(nmsg) = m
    end if
    rc = max(rc, level)
  end subroutine note

  integer function stn(s)
    character(len=*), intent(in) :: s
    integer :: a
    stn = 0
    do a = 1, ns
      if (trim(sid(a)) == trim(s)) stn = a
    end do
  end function stn

  subroutine open_in(name, unit, need)
    character(len=*), intent(in) :: name
    integer, intent(out) :: unit
    logical, intent(in) :: need
    open(newunit=unit, file=name, status='old', iostat=ios)
    if (ios /= 0) then
      unit = 0             ! newunit numbers are negative
      if (need) then
        write(*,'(A)') 'TIDE012E CANNOT OPEN ' // name
        call exit(12)
      end if
      return
    end if
    read(unit, '(A)', iostat=ios) line          ! header
  end subroutine open_in

  !-- input ---------------------------------------------------------
  subroutine read_status()
    integer :: v(8)
    call open_in('fetch_status.csv', u, .false.)
    if (u /= 0) then
      do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        call split(line, f, nf)
        if (f(1) == 'run_utc') runutc = f(2)
        if (f(1) == 'pacific_date') pdate = f(2)
      end do
      close(u)
    end if
    if (len_trim(runutc) < 16) then
      call date_and_time(values=v)
      write(runutc, '(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') v(1), '-', &
        v(2), '-', v(3), ' ', v(5), ':', v(6)
      tnow = hrs(runutc) - v(4) / 60.0_dp
      runutc = tstr(tnow)
      runutc = runutc(1:10) // ' ' // runutc(12:16)
      call note('TIDE010W NO fetch_status.csv: USING THE CLOCK', 0)
    end if
    tnow = hrs(runutc)
    if (len_trim(pdate) < 10) pdate = lstr(tnow)
  end subroutine read_status

  subroutine read_stations()
    call open_in('stations.csv', u, .true.)
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      if (nf < 5 .or. ns >= maxs) cycle
      ns = ns + 1
      sid(ns) = f(1)
      sname(ns) = f(2)
      slat(ns) = rval(f(4))
      slon(ns) = rval(f(5))
    end do
    close(u)
  end subroutine read_stations

  subroutine read_harcon()
    call open_in('harcon.csv', u, .true.)
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      i = stn(f(1))
      if (i == 0) cycle
      do k = 1, nc
        if (trim(f(2)) == trim(cname(k))) then
          st(i)%amp(k) = rval(f(3))
          st(i)%kap(k) = rval(f(4))
          nhc(i) = nhc(i) + 1
        end if
      end do
    end do
    close(u)
  end subroutine read_harcon

  ! datums.csv is relative to the station datum; keep everything as
  ! feet above mean lower low water
  subroutine read_datums()
    real(dp) :: raw(12, maxs)
    character(len=16) :: w(2, maxs)
    integer, parameter :: imllw = 1
    raw = miss
    w = ' '
    call open_in('datums.csv', u, .true.)
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      i = stn(f(1))
      if (i == 0) cycle
      x = rval(f(3))
      select case (trim(f(2)))
      case ('MLLW')
        raw(imllw, i) = x
      case ('MHHW')
        raw(2, i) = x
      case ('MHW')
        raw(3, i) = x
      case ('MSL')
        raw(4, i) = x
      case ('MLW')
        raw(5, i) = x
      case ('GT')
        raw(6, i) = x
      case ('HAT')
        raw(7, i) = x
      case ('LAT')
        raw(8, i) = x
      case ('max')
        raw(9, i) = x
        w(1, i) = f(4)
      case ('min')
        raw(10, i) = x
        w(2, i) = f(4)
      end select
    end do
    close(u)
    do i = 1, ns
      if (.not. ok(raw(imllw, i))) then
        call note('TIDE013W ' // trim(sname(i)) // &
                  ': NO DATUMS, TIDE LEFT ABOUT MEAN SEA LEVEL', 4)
        cycle
      end if
      dmhhw(i) = rel(raw(2, i), raw(imllw, i))
      dmhw(i) = rel(raw(3, i), raw(imllw, i))
      dmsl(i) = rel(raw(4, i), raw(imllw, i))
      dmlw(i) = rel(raw(5, i), raw(imllw, i))
      dgt(i) = raw(6, i)
      dhat(i) = rel(raw(7, i), raw(imllw, i))
      dlat(i) = rel(raw(8, i), raw(imllw, i))
      dmax(i) = rel(raw(9, i), raw(imllw, i))
      dmin(i) = rel(raw(10, i), raw(imllw, i))
      dmaxw(i) = w(1, i)
      dminw(i) = w(2, i)
      if (ok(dmsl(i))) st(i)%z0 = dmsl(i)
    end do
  end subroutine read_datums

  real(dp) function rel(y, base)
    real(dp), intent(in) :: y, base
    rel = miss
    if (ok(y) .and. ok(base)) rel = y - base
  end function rel

  subroutine read_levels()
    call open_in('water_level.csv', u, .false.)
    if (u == 0) then
      call note('TIDE014W NO water_level.csv', 4)
      return
    end if
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      i = stn(f(1))
      if (i == 0) cycle
      x = rval(f(3))
      if (.not. ok(x) .or. no(i) >= maxo) cycle
      if (.not. tok(hrs(f(2)))) cycle
      no(i) = no(i) + 1
      ot(no(i), i) = hrs(f(2))
      oh(no(i), i) = x
    end do
    close(u)
  end subroutine read_levels

  subroutine read_pressure()
    call open_in('pressure.csv', u, .false.)
    if (u == 0) return
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      i = stn(f(1))
      if (i == 0) cycle
      x = rval(f(3))
      if (.not. ok(x) .or. np(i) >= maxo) cycle
      if (x < 900.0_dp .or. x > 1090.0_dp) cycle
      if (.not. tok(hrs(f(2)))) cycle
      np(i) = np(i) + 1
      pt(np(i), i) = hrs(f(2))
      pv(np(i), i) = x
    end do
    close(u)
  end subroutine read_pressure

  subroutine read_monthly()
    integer :: yy, mo
    call open_in('monthly_mean.csv', u, .false.)
    if (u == 0) then
      call note('TIDE015W NO monthly_mean.csv: NO SEA-LEVEL TREND', 4)
      return
    end if
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      i = stn(f(1))
      if (i == 0) cycle
      x = rval(f(4))
      read(f(2), *, iostat=ios) yy
      if (ios /= 0) cycle
      read(f(3), *, iostat=ios) mo
      if (ios /= 0 .or. .not. ok(x) .or. nm(i) >= maxm) cycle
      nm(i) = nm(i) + 1
      my(nm(i), i) = yy
      mmo(nm(i), i) = mo
      mv(nm(i), i) = x
    end do
    close(u)
  end subroutine read_monthly

  subroutine read_noaa_hilo()
    call open_in('noaa_hilo.csv', u, .false.)
    if (u == 0) return
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      i = stn(f(1))
      if (i == 0 .or. nn(i) >= maxe) cycle
      if (.not. tok(hrs(f(2)))) cycle
      nn(i) = nn(i) + 1
      nt(nn(i), i) = hrs(f(2))
      nh(nn(i), i) = rval(f(3))
      nty(nn(i), i) = f(4)(1:1)
    end do
    close(u)
  end subroutine read_noaa_hilo

  subroutine read_noaa_trend()
    call open_in('noaa_trend.csv', u, .false.)
    if (u == 0) return
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      i = stn(f(1))
      if (i == 0) cycle
      ntr(i) = rval(f(2))
      nte(i) = rval(f(3))
    end do
    close(u)
  end subroutine read_noaa_trend

  subroutine read_normals()
    integer :: d
    nrs = miss
    nmx = miss
    call open_in('normals.csv', u, .false.)
    if (u /= 0) then
      do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        call split(line, f, nf)
        i = stn(f(1))
        if (i == 0 .or. nf < 15) cycle
        read(f(2), *, iostat=ios) d
        if (ios /= 0 .or. d < 1 .or. d > 366) cycle
        do k = 1, 5
          nrs(k, d, i) = rval(f(4 + k))
          nmx(k, d, i) = rval(f(10 + k))
        end do
        if (ok(nrs(3, d, i))) have_nrm(i) = .true.
      end do
      close(u)
    end if
    ! no normals.csv at all is a warning; a gauge missing from it has
    ! no 1991-2020 record (Bremerton), which is a fact, not a fault
    do i = 1, ns
      if (have_nrm(i)) cycle
      if (u == 0) then
        call note('TIDE018W ' // trim(sname(i)) // &
          ': NO NORMALS ON FILE', 4)
      else
        call note('TIDE018I ' // trim(sname(i)) // &
          ': NO 1991-2020 RECORD, NO NORMALS', 0)
      end if
    end do
  end subroutine read_normals

  !-- the water now -------------------------------------------------
  integer function doy_of(t)
    real(dp), intent(in) :: t
    integer :: y, m, d
    call ymd(jday(t), y, m, d)
    doy_of = doy366(m, d)
  end function doy_of

  subroutine water_now(i)
    integer, intent(in) :: i
    integer :: a, n24, nym
    real(dp) :: s24, mid, off, ps, dum
    if (no(i) == 0) then
      call note('TIDE015W ' // trim(sname(i)) // &
                ': NO WATER LEVELS', 4)
      return
    end if
    do a = 1, no(i)
      call tide_at(st(i), ot(a, i), op(a, i))
    end do
    tl(i) = ot(no(i), i)
    lv(i) = oh(no(i), i)
    lp(i) = op(no(i), i)
    if (tl(i) < tnow - 3.0_dp) call note('TIDE016W ' // &
      trim(sname(i)) // ': LATEST WATER LEVEL ' // &
      trim(lstr(tl(i))), 4)
    ! surge: the last 24 hours, and the whole 30 days
    s24 = 0
    n24 = 0
    rmean(i) = 0
    rmax(i) = -huge(1.0_dp)
    rmin(i) = huge(1.0_dp)
    do a = 1, no(i)
      x = oh(a, i) - op(a, i)
      rmean(i) = rmean(i) + x
      if (x > rmax(i)) then
        rmax(i) = x
        trmax(i) = ot(a, i)
      end if
      if (x < rmin(i)) then
        rmin(i) = x
        trmin(i) = ot(a, i)
      end if
      if (ot(a, i) > tl(i) - 24.0_dp) then
        s24 = s24 + x
        n24 = n24 + 1
      end if
    end do
    rmean(i) = rmean(i) / no(i)
    if (n24 >= 200) res24(i) = s24 / n24
    ! air pressure: the inverse barometer, a foot of water for every
    ! 30 hPa, about a centimetre per hectopascal
    ps = 0
    a = 0
    do k = 1, np(i)
      if (pt(k, i) > tl(i) - 24.0_dp .and. pt(k, i) <= tl(i)) then
        ps = ps + pv(k, i)
        a = a + 1
      end if
    end do
    if (a >= 100) then
      hpa(i) = ps / a
      ib24(i) = (1013.25_dp - hpa(i)) * hpa_ft
    end if
    ! yesterday's highest water, Pacific time
    off = pac_off(tnow)
    mid = (jday(tnow + off) - 2451545) * 24.0_dp - 12.0_dp - off
    nym = 0
    dum = -huge(1.0_dp)
    do a = 1, no(i)
      if (ot(a, i) >= mid - 24.0_dp .and. ot(a, i) < mid) then
        nym = nym + 1
        dum = max(dum, oh(a, i))
      end if
    end do
    if (nym >= 200) ymax(i) = dum
    ydoy(i) = doy_of(mid - 12.0_dp + off)
    tdoy(i) = doy_of(tnow + off)
    if (have_nrm(i)) then
      res_vs(i) = vsn(res24(i), nrs(:, tdoy(i), i))
      ymax_vs(i) = vsn(ymax(i), nmx(:, ydoy(i), i))
    end if
  end subroutine water_now

  character(len=10) function vsn(v, q)
    real(dp), intent(in) :: v, q(5)
    real(dp) :: r
    ! compare at the precision the files store
    r = miss
    if (ok(v)) r = nint(v * 100.0_dp) / 100.0_dp
    vsn = versus(r, q(1), q(2), q(4), q(5))
  end function vsn

  ! highs and lows from a day ago to eight days ahead
  subroutine week(i)
    integer, intent(in) :: i
    call events(st(i), tnow - 24.0_dp, tnow + 192.0_dp, maxe, ne(i), &
                et(:, i), eh(:, i), ety(:, i))
  end subroutine week

  ! FORTRAN's highs and lows against NOAA's tables for the same week
  subroutine verify(i)
    integer, intent(in) :: i
    integer :: a, b, best
    real(dp) :: d
    vtot(i) = 0
    if (nn(i) == 0) return
    vdt(i) = 0
    vdh(i) = 0
    do a = 1, nn(i)
      ! only the week both lists cover
      if (nt(a, i) < et(1, i) + 0.5_dp .or. &
          nt(a, i) > et(ne(i), i) - 0.5_dp) cycle
      vtot(i) = vtot(i) + 1
      best = 0
      do b = 1, ne(i)
        if (ety(b, i) /= nty(a, i)) cycle
        if (abs(et(b, i) - nt(a, i)) > 1.0_dp) cycle
        if (best == 0) then
          best = b
        else if (abs(et(b, i) - nt(a, i)) < &
                 abs(et(best, i) - nt(a, i))) then
          best = b
        end if
      end do
      if (best == 0) cycle
      vn(i) = vn(i) + 1
      d = abs(et(best, i) - nt(a, i)) * 60.0_dp
      vdt(i) = max(vdt(i), d)
      vdh(i) = max(vdh(i), abs(eh(best, i) - nh(a, i)))
    end do
    ! and FORTRAN's own events in that window that NOAA doesn't list
    do b = 1, ne(i)
      if (et(b, i) < nt(1, i) + 0.5_dp .or. &
          et(b, i) > nt(nn(i), i) - 0.5_dp) cycle
      best = 0
      do a = 1, nn(i)
        if (nty(a, i) == ety(b, i) .and. &
            abs(nt(a, i) - et(b, i)) < 1.0_dp) best = a
      end do
      if (best == 0) vext(i) = vext(i) + 1
    end do
    if (vext(i) > 0) call note('TIDE021W ' // trim(sname(i)) // &
      ': ' // trim(itoa(vext(i))) // ' HIGHS OR LOWS NOAA DOES ' // &
      'NOT LIST', 4)
    if (vn(i) < vtot(i) .or. vdt(i) > 3.0_dp .or. vdh(i) > 0.05_dp) &
      call note('TIDE020W ' // trim(sname(i)) // ': FORTRAN AND ' // &
        'NOAA DIFFER: ' // trim(itoa(vn(i))) // ' OF ' // &
        trim(itoa(vtot(i))) // ' MATCHED, MAX' // rf(vdt(i), 5, 1) &
        // ' MIN,' // rf(vdh(i), 6, 3) // ' FT', 4)
  end subroutine verify

  ! a year of highs and lows: king tides and daylight minus tides
  subroutine year_ahead(i)
    integer, intent(in) :: i
    integer :: a, b, c
    logical :: used(maxy)
    call events(st(i), tnow, tnow + 365.25_dp * 24.0_dp, maxy, &
                nyr(i), yt(:, i), yh(:, i), yty(:, i))
    ! the highest highs, no two within three days of each other, so
    ! each one is a different spell of king tides
    used = .false.
    nki(i) = 0
    do c = 1, maxk
      b = 0
      do a = 1, nyr(i)
        if (yty(a, i) /= 'H' .or. used(a)) cycle
        if (b == 0) then
          b = a
        else if (yh(a, i) > yh(b, i)) then
          b = a
        end if
      end do
      if (b == 0) exit
      nki(i) = nki(i) + 1
      kidx(nki(i), i) = b
      do a = 1, nyr(i)
        if (abs(yt(a, i) - yt(b, i)) < 72.0_dp) used(a) = .true.
      end do
    end do
    ! lows below zero while the sun is up
    nminus(i) = 0
    iminus1(i) = 0
    iminlow(i) = 0
    do a = 1, nyr(i)
      if (.not. minus_day(i, a)) cycle
      nminus(i) = nminus(i) + 1
      if (iminus1(i) == 0) iminus1(i) = a
      if (iminlow(i) == 0) then
        iminlow(i) = a
      else if (yh(a, i) < yh(iminlow(i), i)) then
        iminlow(i) = a
      end if
    end do
  end subroutine year_ahead

  logical function minus_day(i, a)
    integer, intent(in) :: i, a
    minus_day = .false.
    if (yty(a, i) /= 'L' .or. yh(a, i) >= 0.0_dp) return
    minus_day = sun_elev(yt(a, i), slat(i), slon(i)) > -0.833_dp
  end function minus_day

  ! the sea-level trend over the whole record
  subroutine sea_level(i)
    integer, intent(in) :: i
    real(dp) :: yr(maxm)
    integer :: a
    if (nm(i) < 240) then
      if (nm(i) > 0) call note('TIDE019I ' // trim(sname(i)) // &
        ': ' // trim(itoa(nm(i))) // ' MONTHLY MEANS, TOO FEW ' // &
        'FOR A TREND', 0)
      return
    end if
    do a = 1, nm(i)
      yr(a) = my(a, i) + (mmo(a, i) - 0.5_dp) / 12.0_dp
    end do
    call trend(nm(i), yr(1:nm(i)), mmo(1:nm(i), i), mv(1:nm(i), i), &
               tb(i), tse(i), trho(i))
    ty0(i) = my(1, i)
    ty1(i) = my(nm(i), i)
    tnm(i) = nm(i)
  end subroutine sea_level

  !-- output --------------------------------------------------------
  subroutine csv_open(name, head)
    character(len=*), intent(in) :: name, head
    open(newunit=u, file=name // '.tmp', status='replace')
    write(u, '(A)') head
  end subroutine csv_open

  subroutine csv_close(name)
    character(len=*), intent(in) :: name
    close(u)
    call rename(name // '.tmp', name)
  end subroutine csv_close

  ! the next event of a type after tnow, as an index into et
  integer function next_ev(i, c)
    integer, intent(in) :: i
    character(len=1), intent(in) :: c
    integer :: a
    next_ev = 0
    do a = 1, ne(i)
      if (et(a, i) > tnow .and. ety(a, i) == c) then
        next_ev = a
        return
      end if
    end do
  end function next_ev

  subroutine write_now()
    integer :: nh_, nl_
    call csv_open('now.csv', 'station,name,latest_utc,level_ft,' // &
      'tide_ft,surge_ft,surge24_ft,surge24_p10,surge24_p50,' // &
      'surge24_p90,surge24_vs_normal,pressure_hpa,barometer_ft,' // &
      'max_surge30_ft,max_surge30_utc,min_surge30_ft,' // &
      'min_surge30_utc,mean_surge30_ft,yesterday_high_ft,' // &
      'yesterday_high_p50,yesterday_high_p90,' // &
      'yesterday_high_vs_normal,' &
      // 'next_high_utc,next_high_ft,next_low_utc,next_low_ft,mhhw_ft')
    do i = 1, ns
      nh_ = next_ev(i, 'H')
      nl_ = next_ev(i, 'L')
      line = trim(sid(i)) // ',' // trim(sname(i)) // ','
      if (ok(tl(i))) line = trim(line) // tstr(tl(i))
      line = trim(line) // ',' // cs(lv(i), 2) // ',' // cs(lp(i), 2) &
        // ',' // cs(dif(lv(i), lp(i)), 2) // ',' // cs(res24(i), 2) &
        // ',' // cs(nrs(1, max(1, tdoy(i)), i), 2) // ',' // &
        cs(nrs(3, max(1, tdoy(i)), i), 2) // ',' // &
        cs(nrs(5, max(1, tdoy(i)), i), 2) // ',' // trim(res_vs(i)) &
        // ',' // cs(hpa(i), 1) // ',' // cs(ib24(i), 2) // ',' // &
        cs(rmax(i), 2) // ',' // tstr_ok(trmax(i)) // ',' // &
        cs(rmin(i), 2) // ',' // tstr_ok(trmin(i)) // ',' // &
        cs(rmean(i), 2) // ',' // cs(ymax(i), 2) // ',' // &
        cs(nmx(3, max(1, ydoy(i)), i), 2) // ',' // &
        cs(nmx(5, max(1, ydoy(i)), i), 2) // ',' // trim(ymax_vs(i)) &
        // ',' // ev_t(i, nh_) // ',' // ev_h(i, nh_) // ',' // &
        ev_t(i, nl_) // ',' // ev_h(i, nl_) // ',' // cs(dmhhw(i), 2)
      write(u, '(A)') trim(line)
    end do
    call csv_close('now.csv')
  end subroutine write_now

  real(dp) function dif(a, b)
    real(dp), intent(in) :: a, b
    dif = miss
    if (ok(a) .and. ok(b)) dif = a - b
  end function dif

  function tstr_ok(t) result(s)
    real(dp), intent(in) :: t
    character(len=:), allocatable :: s
    s = ''
    if (ok(t)) s = tstr(t)
  end function tstr_ok

  function ev_t(i, a) result(s)
    integer, intent(in) :: i, a
    character(len=:), allocatable :: s
    s = ''
    if (a > 0) s = tstr(et(a, i))
  end function ev_t

  function ev_h(i, a) result(s)
    integer, intent(in) :: i, a
    character(len=:), allocatable :: s
    s = ''
    if (a > 0) s = cs(eh(a, i), 2)
  end function ev_h

  ! hourly water level and tide: a week back, three days ahead
  subroutine write_series()
    integer :: hh, a
    real(dp) :: t, h, ob
    call csv_open('series.csv', &
      'station,time_utc,level_ft,tide_ft,surge_ft')
    t0i = floor(tnow) - 7 * 24
    do i = 1, ns
      if (nhc(i) < 30) cycle
      a = 1
      do hh = t0i, t0i + 10 * 24
        t = real(hh, dp)
        call tide_at(st(i), t, h)
        ob = miss
        do while (a < no(i))
          if (ot(a, i) >= t - 0.01_dp) exit
          a = a + 1
        end do
        if (a <= no(i)) then
          if (abs(ot(a, i) - t) < 0.01_dp) ob = oh(a, i)
        end if
        write(u, '(A)') trim(sid(i)) // ',' // tstr(t) // ',' // &
          cs(ob, 2) // ',' // cs(h, 2) // ',' // cs(dif(ob, h), 2)
      end do
    end do
    call csv_close('series.csv')
  end subroutine write_series

  subroutine write_hilo()
    integer :: a, b, best
    character(len=20) :: ls
    call csv_open('hilo.csv', 'station,time_utc,local,type,ft,' // &
      'noaa_ft,noaa_minutes')
    do i = 1, ns
      do a = 1, ne(i)
        if (et(a, i) < tnow - 12.0_dp) cycle
        best = 0
        do b = 1, nn(i)
          if (nty(b, i) == ety(a, i) .and. &
              abs(nt(b, i) - et(a, i)) < 1.0_dp) best = b
        end do
        ls = lstr(et(a, i))
        line = trim(sid(i)) // ',' // tstr(et(a, i)) // ',' // &
          trim(ls) // ',' // ety(a, i) // ',' // cs(eh(a, i), 2)
        if (best > 0) then
          line = trim(line) // ',' // cs(nh(best, i), 2) // ',' // &
            cs((et(a, i) - nt(best, i)) * 60.0_dp, 1)
        else
          line = trim(line) // ',,'
        end if
        write(u, '(A)') trim(line)
      end do
    end do
    call csv_close('hilo.csv')
  end subroutine write_hilo

  subroutine write_year()
    integer :: a, b
    call csv_open('kingtides.csv', 'station,rank,time_utc,local,ft,' &
      // 'above_mhhw_ft')
    do i = 1, ns
      do a = 1, nki(i)
        b = kidx(a, i)
        write(u, '(A)') trim(sid(i)) // ',' // trim(itoa(a)) // ',' &
          // tstr(yt(b, i)) // ',' // trim(lstr(yt(b, i))) // ',' // &
          cs(yh(b, i), 2) // ',' // cs(dif(yh(b, i), dmhhw(i)), 2)
      end do
    end do
    call csv_close('kingtides.csv')
    call csv_open('minustides.csv', 'station,time_utc,local,ft')
    do i = 1, ns
      do a = 1, nyr(i)
        if (minus_day(i, a)) write(u, '(A)') trim(sid(i)) // ',' // &
          tstr(yt(a, i)) // ',' // trim(lstr(yt(a, i))) // ',' // &
          cs(yh(a, i), 2)
      end do
    end do
    call csv_close('minustides.csv')
  end subroutine write_year

  ! how the tide travels in from the ocean: the main lunar tide's
  ! phase and size at each station against Neah Bay's
  subroutine write_travel()
    call csv_open('travel.csv', 'station,name,lat,lon,m2_ft,' // &
      'm2_phase,lag_hours,m2_ratio,k1_ft,great_range_ft,mhhw_ft,' // &
      'record_high_ft,record_high_when,record_low_ft,record_low_when')
    do i = 1, ns
      write(u, '(A)') trim(sid(i)) // ',' // trim(sname(i)) // ',' // &
        cs(slat(i), 4) // ',' // cs(slon(i), 4) // ',' // &
        cs(st(i)%amp(1), 3) // ',' // cs(st(i)%kap(1), 1) // ',' // &
        cs(lag(i), 2) // ',' // cs(st(i)%amp(1) / st(1)%amp(1), 2) &
        // ',' // cs(st(i)%amp(4), 3) // ',' // cs(dgt(i), 2) // ',' &
        // cs(dmhhw(i), 2) // ',' // cs(dmax(i), 2) // ',' // &
        trim(dmaxw(i)) // ',' // cs(dmin(i), 2) // ',' // trim(dminw(i))
    end do
    call csv_close('travel.csv')
  end subroutine write_travel

  real(dp) function lag(i)
    integer, intent(in) :: i
    lag = modulo(st(i)%kap(1) - st(1)%kap(1), 360.0_dp) / m2spd
  end function lag

  subroutine write_trend()
    integer :: a, y, c
    real(dp) :: s
    call csv_open('trend.csv', 'station,name,first_year,last_year,' &
      // 'months,mm_per_yr,ci95_mm_per_yr,in_per_century,lag1,' // &
      'noaa_mm_per_yr,noaa_err_mm_per_yr,line_first_ft,line_last_ft')
    do i = 1, ns
      write(u, '(A)') trim(sid(i)) // ',' // trim(sname(i)) // ',' // &
        trim(itoa(ty0(i))) // ',' // trim(itoa(ty1(i))) // ',' // &
        trim(itoa(tnm(i))) // ',' // cs(sc(tb(i), ft2mm), 2) // ',' &
        // cs(sc(tse(i), 1.96_dp * ft2mm), 2) // ',' // &
        cs(sc(tb(i), 1200.0_dp), 1) // ',' // cs(trho(i), 2) // ',' &
        // cs(ntr(i), 2) // ',' // cs(nte(i), 2) // ',' // &
        cs(fit(i, ty0(i)), 3) // ',' // cs(fit(i, ty1(i)), 3)
    end do
    call csv_close('trend.csv')
    ! annual means, for the chart: years with at least 10 months
    call csv_open('annual.csv', 'station,year,msl_ft,months')
    do i = 1, ns
      a = 1
      do while (a <= nm(i))
        y = my(a, i)
        s = 0
        c = 0
        do while (a <= nm(i))
          if (my(a, i) /= y) exit
          s = s + mv(a, i)
          c = c + 1
          a = a + 1
        end do
        if (c >= 10) write(u, '(A)') trim(sid(i)) // ',' // &
          trim(itoa(y)) // ',' // cs(s / c, 3) // ',' // trim(itoa(c))
      end do
    end do
    call csv_close('annual.csv')
  end subroutine write_trend

  ! the trend line at the middle of a year: through the mean of all
  ! the months, at the trend's slope
  real(dp) function fit(i, y)
    integer, intent(in) :: i, y
    real(dp) :: mt, mx_
    integer :: a
    fit = miss
    if (.not. ok(tb(i)) .or. nm(i) == 0) return
    mt = 0
    mx_ = 0
    do a = 1, nm(i)
      mt = mt + my(a, i) + (mmo(a, i) - 0.5_dp) / 12.0_dp
      mx_ = mx_ + mv(a, i)
    end do
    fit = mx_ / nm(i) + tb(i) * (y + 0.5_dp - mt / nm(i))
  end function fit

  real(dp) function sc(v, k)
    real(dp), intent(in) :: v, k
    sc = miss
    if (ok(v)) sc = v * k
  end function sc

  subroutine write_summary()
    integer :: s, a, kmin
    real(dp) :: dtm, dhm
    s = stn('9447130')
    dtm = 0
    dhm = 0
    do i = 1, ns
      if (ok(vdt(i))) dtm = max(dtm, vdt(i))
      if (ok(vdh(i))) dhm = max(dhm, vdh(i))
    end do
    call csv_open('summary.csv', 'key,value')
    call kv('version', ver)
    call kv('run_utc', trim(runutc))
    call kv('pacific_date', pdate(1:10))
    call kv('return_code', trim(itoa(rc)))
    call kv('stations', trim(itoa(ns)))
    call kv('stations_current', trim(itoa(current)))
    call kv('stations_with_normals', trim(itoa(count(have_nrm(1:ns)))))
    call kv('verify_events', trim(itoa(sum(vn(1:ns)))) // ' of ' // &
      trim(itoa(sum(vtot(1:ns)))))
    call kv('verify_extra_events', trim(itoa(sum(vext(1:ns)))))
    call kv('verify_max_minutes', cs(dtm, 1))
    call kv('verify_max_ft', cs(dhm, 3))
    if (s > 0) then
      call kv('seattle_latest_utc', tstr_ok(tl(s)))
      call kv('seattle_level_ft', cs(lv(s), 2))
      call kv('seattle_tide_ft', cs(lp(s), 2))
      call kv('seattle_surge_ft', cs(dif(lv(s), lp(s)), 2))
      call kv('seattle_surge24_ft', cs(res24(s), 2))
      call kv('seattle_surge24_vs_normal', trim(res_vs(s)))
      a = next_ev(s, 'H')
      call kv('seattle_next_high_utc', ev_t(s, a))
      call kv('seattle_next_high_ft', ev_h(s, a))
      a = next_ev(s, 'L')
      call kv('seattle_next_low_utc', ev_t(s, a))
      call kv('seattle_next_low_ft', ev_h(s, a))
      if (nki(s) > 0) then
        a = kidx(1, s)
        call kv('seattle_king_tide_utc', tstr(yt(a, s)))
        call kv('seattle_king_tide_ft', cs(yh(a, s), 2))
      end if
      if (iminus1(s) > 0) then
        call kv('seattle_next_daylight_minus_utc', &
                tstr(yt(iminus1(s), s)))
        call kv('seattle_next_daylight_minus_ft', &
                cs(yh(iminus1(s), s), 2))
      end if
      kmin = iminlow(s)
      if (kmin > 0) then
        call kv('seattle_lowest_daylight_utc', tstr(yt(kmin, s)))
        call kv('seattle_lowest_daylight_ft', cs(yh(kmin, s), 2))
      end if
      call kv('seattle_daylight_minus_tides', trim(itoa(nminus(s))))
      call kv('seattle_trend_mm_per_yr', cs(sc(tb(s), ft2mm), 2))
      call kv('seattle_trend_ci95', cs(sc(tse(s), 1.96_dp*ft2mm), 2))
      call kv('seattle_trend_since', trim(itoa(ty0(s))))
    end if
    call kv('messages', trim(itoa(nmsg)))
    call csv_close('summary.csv')
  end subroutine write_summary

  subroutine kv(k_, v_)
    character(len=*), intent(in) :: k_, v_
    write(u, '(A)') k_ // ',' // trim(v_)
  end subroutine kv

  !-- the printed report, 132 columns ---------------------------------
  subroutine write_report()
    character(len=132) :: rule, dash
    character(len=20) :: ls
    integer :: a, b, nh_, nl_, s
    rule = repeat('=', 132)
    dash = repeat('-', 132)
    open(newunit=u, file='puget-tides-report.tmp', status='replace')
    write(u, '(A)') rule
    write(u, '(A,T45,A,T106,A)') ver, &
      'TIDES OF PUGET SOUND   NEAH BAY TO TACOMA', &
      'RUN ' // trim(runutc) // ' UTC'
    write(u, '(A,T45,A)') 'NOAA CO-OPS WATER LEVELS', &
      'FEET ABOVE MEAN LOWER LOW WATER (MLLW)   TIMES PACIFIC'
    write(u, '(A)') rule
    write(u, '(A)')

    write(u, '(A)') 'SECTION I    THE WATER NOW    ' // &
      'MEASURED LEVEL, ' // &
      'PREDICTED TIDE, AND THE DIFFERENCE: THE SURGE'
    write(u, '(A)')
    write(u, '(A)') 'STATION          LATEST              LEVEL' // &
      '    TIDE   SURGE  24-HR SURGE  VS NORMAL     PRESSURE  ' // &
      'BAROMETER   NEXT HIGH                NEXT LOW'
    write(u, '(A)') '                                        FT' // &
      '      FT      FT           FT                     HPA  ' // &
      '       FT'
    write(u, '(A)') dash
    do i = 1, ns
      nh_ = next_ev(i, 'H')
      nl_ = next_ev(i, 'L')
      ls = '--'
      if (ok(tl(i))) ls = lstr(tl(i))
      write(u, '(A16,1X,A20,3A8,A13,2X,A10,A13,A11,3X,A25,A25)') &
        sname(i), ls, rf(lv(i), 8, 2), rf(lp(i), 8, 2), &
        rf(dif(lv(i), lp(i)), 8, 2), rf(res24(i), 13, 2), res_vs(i), &
        rf(hpa(i), 13, 1), rf(ib24(i), 11, 2), evs(i, nh_), evs(i, nl_)
    end do
    write(u, '(A)') dash
    write(u, '(A)') '  SURGE IS WHAT THE WEATHER ADDS: WIND, AIR ' // &
      'PRESSURE AND RIVER FLOW. BAROMETER IS THE PART ' // &
      'LOW PRESSURE ' // &
      'EXPLAINS, ABOUT A FOOT PER 30 HPA BELOW 1013.'
    write(u, '(A)') '  VS NORMAL: THE 24-HOUR SURGE AGAINST EVERY ' // &
      'DAY WITHIN 15 DAYS OF THE DATE, 1991-2020. TIDE PREDICTED ' // &
      'BY THIS PROGRAM FROM NOAA''S CONSTITUENTS.'
    write(u, '(A)')

    write(u, '(A)') 'SECTION II   HIGHS AND LOWS, NEXT TWO DAYS'
    write(u, '(A)')
    do i = 1, ns
      line = '  ' // sname(i)
      b = 0
      do a = 1, ne(i)
        if (et(a, i) < tnow .or. et(a, i) > tnow + 48.0_dp) cycle
        b = b + 1
        if (b > 6) exit
        ls = lstr(et(a, i))
        line(19 + (b - 1) * 19:) = ety(a, i) // ' ' // ls(6:7) // &
          '/' // ls(9:10) // ' ' // ls(12:16) // rf(eh(a, i), 5, 1)
      end do
      write(u, '(A)') trim(line)
    end do
    write(u, '(A)')

    write(u, '(A)') 'SECTION III  FORTRAN AGAINST NOAA    THIS ' // &
      'PROGRAM''S HIGHS AND LOWS AGAINST NOAA''S PUBLISHED ' // &
      'PREDICTIONS, NEXT 8 DAYS'
    write(u, '(A)')
    write(u, '(A)') 'STATION           MATCHED   LARGEST TIME ' // &
      'DIFFERENCE   LARGEST HEIGHT DIFFERENCE   NOT IN NOAA'
    write(u, '(A)') '                                       MIN' // &
      '                          FT'
    write(u, '(A)') dash(1:104)
    do i = 1, ns
      write(u, '(A16,I6,A,I3,A21,A28,I14)') sname(i), vn(i), ' OF', &
        vtot(i), rf(vdt(i), 21, 1), rf(vdh(i), 28, 3), vext(i)
    end do
    write(u, '(A)') '  NOAA ROUNDS TO THE MINUTE AND THE ' // &
      'THOUSANDTH OF A FOOT.'
    write(u, '(A)')

    write(u, '(A)') 'SECTION IV   THE TIDE COMES IN    THE MAIN ' // &
      'LUNAR TIDE (M2) FROM THE OCEAN TO THE END OF THE SOUND'
    write(u, '(A)')
    write(u, '(A)') 'STATION             M2 FT   PHASE  HOURS ' // &
      'AFTER NEAH BAY   SIZE VS NEAH BAY   GREAT RANGE FT   ' // &
      'RECORD HIGH FT  WHEN'
    write(u, '(A)') dash
    do i = 1, ns
      write(u, '(A16,A9,A8,A21,A19,A17,A18,2X,A16)') sname(i), &
        rf(st(i)%amp(1), 9, 2), rf(st(i)%kap(1), 8, 1), &
        rf(lag(i), 21, 2), rf(st(i)%amp(1) / st(1)%amp(1), 19, 2), &
        rf(dgt(i), 17, 2), rf(dmax(i), 18, 2), dmaxw(i)
    end do
    write(u, '(A)') '  GREAT RANGE: MEAN HIGHER HIGH TO MEAN ' // &
      'LOWER LOW WATER, 1983-2001 NATIONAL TIDAL DATUM EPOCH.'
    write(u, '(A)')

    write(u, '(A)') 'SECTION V    KING TIDES    THE HIGHEST ' // &
      'PREDICTED TIDES OF THE NEXT 12 MONTHS, ONE PER SPELL'
    write(u, '(A)')
    do i = 1, ns
      line = '  ' // sname(i)
      do a = 1, min(5, nki(i))
        b = kidx(a, i)
        ls = lstr(yt(b, i))
        line(19 + (a - 1) * 21:) = mstr(jdn_s(ls(1:10))) // ' ' // &
          ls(12:16) // rf(yh(b, i), 6, 1)
      end do
      write(u, '(A)') trim(line)
    end do
    write(u, '(A)')

    write(u, '(A)') 'SECTION VI   MINUS TIDES IN DAYLIGHT    LOWS ' // &
      'BELOW ZERO WHILE THE SUN IS UP, NEXT 12 MONTHS'
    write(u, '(A)')
    write(u, '(A)') 'STATION          COUNT   NEXT ONE' // &
      '                       LOWEST'
    write(u, '(A)') dash(1:90)
    do i = 1, ns
      line = sname(i) // rf(real(nminus(i), dp), 6, 0)
      if (iminus1(i) > 0) then
        ls = lstr(yt(iminus1(i), i))
        line = line(1:25) // ls // rf(yh(iminus1(i), i), 7, 1)
      end if
      if (iminlow(i) > 0) then
        ls = lstr(yt(iminlow(i), i))
        line = line(1:56) // ls // rf(yh(iminlow(i), i), 7, 1)
      end if
      write(u, '(A)') trim(line)
    end do
    write(u, '(A)') '  SUNRISE AND SUNSET FROM THE SUN''S ' // &
      'POSITION, COMPUTED HERE. IN FALL AND WINTER THE ' // &
      'LOWEST TIDES COME AT NIGHT.'
    write(u, '(A)')

    write(u, '(A)') 'SECTION VII  SEA LEVEL    TREND IN MONTHLY ' // &
      'MEAN SEA LEVEL OVER EACH STATION''S WHOLE RECORD'
    write(u, '(A)')
    write(u, '(A)') 'STATION          FROM    TO   MONTHS   ' // &
      'MM/YR   +/- 95%   INCHES/CENTURY   NOAA MM/YR   +/-'
    write(u, '(A)') dash(1:100)
    do i = 1, ns
      write(u, '(A16,2I6,I9,A8,A10,A17,A13,A6)') sname(i), ty0(i), &
        ty1(i), tnm(i), rf(sc(tb(i), ft2mm), 8, 2), &
        rf(sc(tse(i), 1.96_dp * ft2mm), 10, 2), &
        rf(sc(tb(i), 1200.0_dp), 17, 1), rf(ntr(i), 13, 2), &
        rf(nte(i), 6, 2)
    end do
    write(u, '(A)') '  LEAST SQUARES WITH A MEAN FOR EACH ' // &
      'CALENDAR MONTH. THE UNCERTAINTY IS WIDENED FOR ONE MONTH ' // &
      'FOLLOWING ANOTHER (AR1); NOAA''S MODEL GIVES A NARROWER ONE.'
    write(u, '(A)') '  THE LAND MOVES TOO: THESE ARE RELATIVE SEA ' // &
      'LEVEL, THE WATER AGAINST THE PIER.'
    write(u, '(A)')

    write(u, '(A)') rule
    write(u, '(A)') 'MESSAGES'
    do k = 1, nmsg
      write(u, '(A)') '  ' // trim(msgs(k))
    end do
    if (nmsg == 0) write(u, '(A)') '  NONE'
    s = count(no(1:ns) > 0)
    write(u, '(A)') 'SOURCES      NOAA CO-OPS WATER LEVELS, ' // &
      'HARMONIC CONSTITUENTS, DATUMS AND MONTHLY MEANS'
    write(u, '(A)') 'END OF REPORT   ' // trim(itoa(ns)) // &
      ' STATIONS   ' // trim(itoa(s)) // ' WITH WATER LEVELS' // &
      '   RETURN CODE ' // trim(itoa(rc))
    write(u, '(A)') rule
    close(u)
    call rename('puget-tides-report.tmp', 'puget-tides-report.txt')
  end subroutine write_report

  function evs(i, a) result(s)
    integer, intent(in) :: i, a
    character(len=25) :: s
    character(len=20) :: ls
    s = '--'
    if (a == 0) return
    ls = lstr(et(a, i))
    s = ls(6:16) // rf(eh(a, i), 6, 1) // ' FT'
  end function evs

end program puget_tides

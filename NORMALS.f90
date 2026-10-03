!=====================================================================
! NORMALS.f90   WHAT'S NORMAL FOR THE DATE, 1991-2020
!   READS EVERY HOURLY WATER LEVEL FROM 1991 TO 2020 AT EACH STATION,
!   SUBTRACTS THE TIDE THIS PROGRAM PREDICTS FOR THAT HOUR, AND KEEPS
!   EACH DAY'S MEAN SURGE AND HIGHEST WATER. THEN, FOR EVERY DAY OF
!   THE YEAR, THE 10TH TO 90TH PERCENTILES: SURGE FROM EVERY DAY
!   WITHIN 15 DAYS OF THE DATE, HIGHEST WATER FROM WITHIN 7 DAYS.
!   DAYS RUN MIDNIGHT TO MIDNIGHT PACIFIC STANDARD TIME.
!
!   IN:   stations.csv harcon.csv datums.csv history/files.txt
!         history/<station>_<year>.csv
!   OUT:  normals.csv
!   RC:   0 NORMAL  4 A STATION HAS FEWER THAN 20 YEARS  8 NOTHING
!
!   BUILD gfortran -O2 -o normals TIDE_PHYS.f90 NORMALS.f90
!=====================================================================
program normals
  use tide_util
  use tide_phys
  implicit none
  integer, parameter :: maxs = 12, maxd = 40 * 366, maxw = 1000
  integer, parameter :: wres = 15, wmax = 7, minn = 30
  type(station_tide) :: st(maxs)
  character(len=8) :: sid(maxs)
  character(len=16) :: sname(maxs)
  integer :: ns = 0, nhc(maxs) = 0
  real(dp) :: mllw(maxs) = miss, msl(maxs) = miss
  ! one value per station-day: day of year, mean surge, highest water
  integer :: nd(maxs) = 0, ddoy(maxd, maxs)
  real(dp) :: dres(maxd, maxs), dmax(maxd, maxs)
  logical :: yr(1991:2020, maxs) = .false.
  integer :: hours_read = 0, hours_used = 0, files = 0
  character(len=512) :: line
  character(len=64) :: f(maxf)
  character(len=40) :: fname
  integer :: u, uf, ios, nf, i, k, rc = 0
  real(dp) :: x

  write(*,'(A)') 'NRML000I NORMALS V1.0 STARTED'
  call read_inputs()
  open(newunit=uf, file='history/files.txt', status='old', iostat=ios)
  if (ios /= 0) then
    write(*,'(A)') 'NRML012E CANNOT OPEN history/files.txt'
    call exit(12)
  end if
  do
    read(uf, '(A)', iostat=ios) fname
    if (ios /= 0) exit
    if (len_trim(fname) < 12) cycle
    call one_file(trim(fname))
  end do
  close(uf)
  write(*,'(A,I9)') 'NRML001I FILES READ ...............', files
  write(*,'(A,I9)') 'NRML002I HOURS READ ...............', hours_read
  write(*,'(A,I9)') 'NRML003I HOURS USED ...............', hours_used
  do i = 1, ns
    write(*,'(A,A16,I7,A,I3,A)') 'NRML005I ', sname(i), nd(i), &
      ' DAYS IN ', count(yr(:, i)), ' YEARS'
    if (count(yr(:, i)) < 20) then
      write(*,'(A)') 'NRML010W ' // trim(sname(i)) // &
        ': FEWER THAN 20 YEARS, NO NORMALS'
      rc = 4
    end if
  end do
  if (sum(nd(1:ns)) == 0) then
    write(*,'(A)') 'NRML008E NO DATA. NOTHING WRITTEN.'
    call exit(8)
  end if
  call write_normals()
  write(*,'(A)') 'NRML006I NORMALS WRITTEN: normals.csv'
  write(*,'(A,I2)') 'NRML999I NORMALS ENDED  RC=', rc
  call exit(rc)

contains

  integer function stn(s)
    character(len=*), intent(in) :: s
    integer :: a
    stn = 0
    do a = 1, ns
      if (trim(sid(a)) == trim(s)) stn = a
    end do
  end function stn

  subroutine read_inputs()
    open(newunit=u, file='stations.csv', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'NRML012E CANNOT OPEN stations.csv'
      call exit(12)
    end if
    read(u, '(A)') line
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      if (nf < 5 .or. ns >= maxs) cycle
      ns = ns + 1
      sid(ns) = f(1)
      sname(ns) = f(2)
    end do
    close(u)
    open(newunit=u, file='harcon.csv', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'NRML012E CANNOT OPEN harcon.csv'
      call exit(12)
    end if
    read(u, '(A)') line
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
    open(newunit=u, file='datums.csv', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'NRML012E CANNOT OPEN datums.csv'
      call exit(12)
    end if
    read(u, '(A)') line
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      i = stn(f(1))
      if (i == 0) cycle
      if (f(2) == 'MLLW') mllw(i) = rval(f(3))
      if (f(2) == 'MSL') msl(i) = rval(f(3))
    end do
    close(u)
    do i = 1, ns
      if (ok(mllw(i)) .and. ok(msl(i))) st(i)%z0 = msl(i) - mllw(i)
    end do
  end subroutine read_inputs

  ! one station-year: <station>_<year>.csv of time,ft hourly
  subroutine one_file(name)
    character(len=*), intent(in) :: name
    integer :: s, day, cur, cnt
    real(dp) :: t, h, p, sres, hmax
    s = stn(name(1:index(name, '_') - 1))
    if (s == 0) return
    if (nhc(s) < 30) return
    open(newunit=u, file='history/' // name, status='old', iostat=ios)
    if (ios /= 0) return
    files = files + 1
    read(u, '(A)') line
    cur = -1
    cnt = 0
    sres = 0
    hmax = -huge(1.0_dp)
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call split(line, f, nf)
      t = hrs(f(1))
      x = rval(f(2))
      hours_read = hours_read + 1
      if (.not. tok(t) .or. .not. ok(x)) cycle
      day = jday(t - 8.0_dp)               ! Pacific standard time
      if (day /= cur) then
        call close_day(s, cur, cnt, sres, hmax)
        cur = day
        cnt = 0
        sres = 0
        hmax = -huge(1.0_dp)
      end if
      call tide_at(st(s), t, p)
      h = x
      sres = sres + (h - p)
      hmax = max(hmax, h)
      cnt = cnt + 1
      hours_used = hours_used + 1
    end do
    call close_day(s, cur, cnt, sres, hmax)
    close(u)
  end subroutine one_file

  ! a day counts if it has at least 20 of its 24 hours
  subroutine close_day(s, day, cnt, sres, hmax)
    integer, intent(in) :: s, day, cnt
    real(dp), intent(in) :: sres, hmax
    integer :: yy, mm, dd
    if (day < 0 .or. cnt < 20 .or. nd(s) >= maxd) return
    call ymd(day, yy, mm, dd)
    if (yy < 1991 .or. yy > 2020) return
    nd(s) = nd(s) + 1
    ddoy(nd(s), s) = doy366(mm, dd)
    dres(nd(s), s) = sres / cnt
    dmax(nd(s), s) = hmax
    yr(yy, s) = .true.
  end subroutine close_day

  subroutine write_normals()
    real(dp) :: a(maxw), b(maxw)
    integer :: d, na, nb, j, dd
    character(len=:), allocatable :: row
    open(newunit=u, file='normals.tmp', status='replace')
    write(u, '(A)') 'station,doy,years,n_surge,surge_p10,' // &
      'surge_p25,' // &
      'surge_p50,surge_p75,surge_p90,n_high,high_p10,high_p25,' // &
      'high_p50,high_p75,high_p90'
    do i = 1, ns
      if (count(yr(:, i)) < 20) cycle
      do d = 1, 366
        na = 0
        nb = 0
        do j = 1, nd(i)
          dd = abs(ddoy(j, i) - d)
          dd = min(dd, 366 - dd)            ! the year wraps
          if (dd <= wres .and. na < maxw) then
            na = na + 1
            a(na) = dres(j, i)
          end if
          if (dd <= wmax .and. nb < maxw) then
            nb = nb + 1
            b(nb) = dmax(j, i)
          end if
        end do
        row = trim(sid(i)) // ',' // trim(itoa(d)) // ',' // &
          trim(itoa(count(yr(:, i)))) // ',' // trim(itoa(na))
        row = row // pcts(na, a)
        row = row // ',' // trim(itoa(nb)) // pcts(nb, b)
        write(u, '(A)') row
      end do
    end do
    close(u)
    call rename('normals.tmp', 'normals.csv')
  end subroutine write_normals

  function pcts(n, v) result(s)
    integer, intent(in) :: n
    real(dp), intent(inout) :: v(:)
    character(len=:), allocatable :: s
    real(dp), parameter :: q(5) = [0.10_dp, 0.25_dp, 0.50_dp, &
                                   0.75_dp, 0.90_dp]
    integer :: j
    s = ''
    if (n < minn) then
      s = ',,,,,'
      return
    end if
    call hsort(n, v(1:n))
    do j = 1, 5
      s = s // ',' // cs(pctl(n, v(1:n), q(j)), 2)
    end do
  end function pcts

end program normals

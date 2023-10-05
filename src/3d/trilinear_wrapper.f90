program trilinear_wrapper
  double precision, parameter :: one      = 1.0d0
  double precision, parameter :: four     = 4.0d0
  double precision, parameter :: f14   = one / four

  ! time
  double precision :: tstart, tstop, tcalc
  ! globals for trilinear
  double precision :: lower(3), dxi(3), extent(3), volg(-5:2,-5:2,-5:2)
  double precision, parameter :: nx=40.0d0, ny=40.0d0, nz=40.0d0
  ! trilinear out
  double precision :: weights(0:1,0:1,0:1), points(3,4)
  double precision :: pvol, temp
  integer :: is, js, ks
  ! others
  integer :: i, p

  ! extent(3) => [0, 40[
  call random_number(extent)
  extent = extent * 40.0d0

  dxi = one / ( extent / dble((/nx, ny, nz/)) )
  ! dxi = 1 / ( extent / dble((/nx, ny, nz/))
  call random_number(lower)
  volg = 0.0d0

  tstart = gettime()

  do i = 1, 100000000
    call random_number(pvol)
    call random_number(points)
    temp = f14 * pvol
    do p = 1, 4
      call trilinear(points(:, p), is, js, ks, weights)
      ! volg(ks:ks+1, js:js+1, is:is+1) = volg(ks:ks+1, js:js+1, is:is+1) + temp
      volg(ks:ks+1, js:js+1, is:is+1) = volg(ks:ks+1, js:js+1, is:is+1) + weights * temp
    end do
  end do

  tstop = gettime()
  tcalc = tstop - tstart
  
  print *, tcalc
  ! print *, volg

contains

  pure subroutine trilinear_vec(pos, ii, jj, kk, ww)
    double precision, intent(in)  :: pos(3)
    integer,          intent(out) :: ii, jj, kk
    double precision, intent(out) :: ww(0:1,0:1,0:1)
    double precision              :: xyz(3)
    double precision              :: wt(0:1, 0:1)
    double precision              :: px, py, pz, pxc, pyc, pzc


    ! (i, j, k)
    xyz = (pos - lower) * dxi
    ii = floor(xyz(1))
    jj = floor(xyz(2))
    kk = floor(xyz(3))

    px = xyz(1) - dble(ii)
    pxc = one - px

    py = xyz(2) - dble(jj)
    pyc = one - py

    pz = xyz(3) - dble(kk)
    pzc = one - pz

    w00 = pyc * pxc
    w10 = pyc * px
    w01 = py * pxc
    w11 = py * px

    ! Note order of indices is k,j,i
    ww(0,0,0) = pzc * w00
    ww(0,0,1) = pzc * w10
    ww(0,1,0) = pzc * w01
    ww(0,1,1) = pzc * w11
    ww(1,0,0) = pz * w00
    ww(1,0,1) = pz * w10
    ww(1,1,0) = pz * w01
    ww(1,1,1) = pz * w11

    ! can we do:
    ! wt(0, 0) = pzc * pyc
    ! wt(0, 1) = pzc * py
    ! wt(1, 0) = pz * pyc
    ! wt(1, 1) = pz * py
    ! 
    ! ww(:,:,0) = pxc * wt
    ! ww(:,:,1) = px * wt
    ! end

  end subroutine trilinear_vec

  pure subroutine trilinear(pos, ii, jj, kk, ww)
    double precision, intent(in)  :: pos(3)
    integer,          intent(out) :: ii, jj, kk
    double precision, intent(out) :: ww(0:1,0:1,0:1)
    double precision              :: xyz(3)
    double precision              :: wt(0:1, 0:1)
    double precision              :: px, py, pz, pxc, pyc, pzc


    ! (i, j, k)
    xyz = (pos - lower) * dxi
    ii = floor(xyz(1))
    jj = floor(xyz(2))
    kk = floor(xyz(3))

    px = xyz(1) - dble(ii)
    pxc = one - px

    py = xyz(2) - dble(jj)
    pyc = one - py

    pz = xyz(3) - dble(kk)
    pzc = one - pz

    w00 = pyc * pxc
    w10 = pyc * px
    w01 = py * pxc
    w11 = py * px

    ! Note order of indices is k,j,i
    ww(0,0,0) = pzc * w00
    ww(0,0,1) = pzc * w10
    ww(0,1,0) = pzc * w01
    ww(0,1,1) = pzc * w11
    ww(1,0,0) = pz * w00
    ww(1,0,1) = pz * w10
    ww(1,1,0) = pz * w01
    ww(1,1,1) = pz * w11

    ! can we do:
    ! wt(0, 0) = pzc * pyc
    ! wt(0, 1) = pzc * py
    ! wt(1, 0) = pz * pyc
    ! wt(1, 1) = pz * py
    ! 
    ! ww(:,:,0) = pxc * wt
    ! ww(:,:,1) = px * wt
    ! end

  end subroutine trilinear


  double precision function gettime()

    logical, save :: firstcall = .true.

    integer, parameter :: int32kind = selected_int_kind( 9)
    integer, parameter :: int64kind = selected_int_kind(18)

    integer, parameter :: intkind = int64kind

    integer(kind = intkind) :: count,rate

    double precision, save :: ticktime

    if (firstcall) then

       firstcall = .false.

       call system_clock(count, rate)

       ticktime = 1.0d0/dble(rate)
       gettime  = dble(count)*ticktime

       !          write(*,*) 'Clock resolution is ', ticktime*1.0e6, ', usecs'

    else

       call system_clock(count)

       gettime = dble(count)*ticktime

    end if

  end function gettime

end program trilinear_wrapper

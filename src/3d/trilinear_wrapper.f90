program trilinear_wrapper
  use constants, only : one

  double precision :: lower(3)
  double precision, parameter, dimension=3 :: extent = (\10.0, 12.0, 14.0\)
  double precision, parameter :: dxi=


  one

  !double precision, intent(in)  :: pos(3)
  !integer,          intent(out) :: ii, jj, kk
  !double precision, intent(out) :: ww(0:1,0:1,0:1)
  !double precision              :: xyz(3)
  !double precision              :: wt(0:1, 0:1)
  !double precision              :: px, py, pz, pxc, pyc, pzc
  call trilinear(pos, ii, jj, kk, ww)

end program trilinear_wrapper

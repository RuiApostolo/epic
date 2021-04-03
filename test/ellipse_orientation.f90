program ellipse_orientation
    use hdf5
    use constants, only : pi
    use parcel_container
    use writer
    implicit none

    integer(hid_t) :: group
    double precision :: origin(2) = (/0.0, 0.0/)
    double precision :: extent(2) =  (/0.2, 0.2/)
    integer :: iter
    integer :: grid(2) = (/2, 2/)
    double precision :: angle
    double precision, parameter :: lam = 3.0

    n_parcels = 2
    call parcel_alloc(2)

    parcels%position = 0.0
    parcels%volume = 0.25 * product(extent / (grid - 1))

    do iter = 0, 359

        angle = dble(iter) * pi / 180.0d0

        print *, angle

        ! b11
        parcels%B(1, 1) = lam * cos(angle) ** 2 + 1.0 / lam * sin(angle) ** 2

        ! b12
        parcels%B(1, 2) = 0.5 * (lam - 1.0 / lam) * sin(2.0 * angle)


        call open_h5_file('ellipse_orientation.hdf5')
        call write_h5_scalar_step_attrib(iter, "t", dble(iter))
        call write_h5_scalar_step_attrib(iter, "dt", dble(iter))


        print *, iter

        call write_h5_parcels(iter)
        call close_h5_file

    enddo


    call open_h5_file('ellipse_orientation.hdf5')

    group = open_h5_group("/")
    call write_h5_integer_scalar_attrib(group, "nsteps", 360)
    call h5gclose_f(group, h5err)

    group = open_h5_group("mesh")
    call write_h5_double_vector_attrib(group, "extent", extent)
    call write_h5_double_vector_attrib(group, "origin", origin)
    call write_h5_integer_vector_attrib(group, "grid", grid)
    call h5gclose_f(group, h5err)

    group = open_h5_group("parcel")
    call write_h5_integer_scalar_attrib(group, "n_per_cell", 1)
    call write_h5_logical_attrib(group, "is_random", .false.)
    call write_h5_integer_scalar_attrib(group, "seed", 42)
    call write_h5_logical_attrib(group, "is_elliptic", .true.)
    call write_h5_double_scalar_attrib(group, "lambda", 3.0d0)
    call h5gclose_f(group, h5err)
    call close_h5_file

    call parcel_dealloc()

end program ellipse_orientation

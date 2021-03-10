module parcels
    implicit none

    integer :: n_parcels

    double precision, allocatable, dimension(:) :: &
        x, y,        & ! positions
        dxdt, dydt,  & ! velocitues
        stretch

    contains

        subroutine split(threshold)
            double precision, intent(in) :: threshold


        end subroutine split

        subroutine alloc_parcel_mem(num)
            integer, intent(in) :: num

            allocate(x(num))
            allocate(y(num))
            allocate(dxdt(num))
            allocate(dydt(num))
            allocate(stretch(num))
        end subroutine alloc_parcel_mem

        subroutine dealloc_parcel_mem
            deallocate(x)
            deallocate(y)
            deallocate(dxdt)
            deallocate(dydt)
            deallocate(stretch)
        end subroutine dealloc_parcel_mem

        subroutine create(num)
            integer, intent(in) :: num

        end subroutine create

end module parcels

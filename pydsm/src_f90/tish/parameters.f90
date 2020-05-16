module parameters
    integer, parameter :: dp = kind(1.d0)
    real(dp), parameter :: pi=4.d0*atan(1.d0)

    integer, parameter :: maxnzone = 20
    integer, parameter :: maxnr = 1000
    integer, parameter :: maxlmax = 80000
    integer, parameter :: maxnlay = 88300
    integer, parameter :: ilog = 0
    real(dp), parameter :: lmaxdivf = 2.d4
    real(dp), parameter :: shallowdepth = 100.d0

    real(dp), parameter :: flattening = 1.d0 / 298.25d0
    ! specfem flattening is f = 1/299.8

end module parameters
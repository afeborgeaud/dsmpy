module parameters
    integer, parameter :: dp = kind(1.d0) !selected_real_kind=12
    real(dp), parameter :: pi=3.1415926535897932d0

    integer, parameter :: maxnzone = 20
    integer, parameter :: maxnr = 1000
    integer, parameter :: maxlmax = 80000
    integer, parameter :: maxnlay = 80880
    integer, parameter :: maxnslay = 48840
    integer, parameter :: maxnllay = 32040
    integer, parameter :: ilog = 0
    real(dp), parameter :: lmaxdivf = 2.d4
    real(dp), parameter :: shallowdepth = 100.d0

end module parameters
